/*
PRICE was written by Graham Ruby.
This file is part of PRICE.

PRICE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PRICE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PRICE.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ASSEMBLYJOBSUBSET_CPP
#define ASSEMBLYJOBSUBSET_CPP


#include "AssemblyJobSubset.h"
#include "ScoredSeqFlip.h"
#include "AssemblyException.h"
#include "AlignmentUngapped.h"
#include <queue>
#include <string>
#include <iostream>
#include <typeinfo>
#include <algorithm>
#include <omp.h>
#include "MatchSeqTest.h"
using namespace::std;

AssemblyJobSubset::AssemblyJobSubset(){}
AssemblyJobSubset::AssemblyJobSubset(set<ScoredSeq*>* seqs,
				     float minFractId,
				     AssemblyJob::AssemblyStrandedness strandedness) :
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _strandedness(strandedness),
  _minFractId(minFractId){
  _inputSeqs.insert( seqs->begin(), seqs->end() );
}


AssemblyJobSubset::~AssemblyJobSubset(){
}


void AssemblyJobSubset::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJSubset::makeShallow can't be called after the job was run.");
  }
  set<ScoredSeq*> shallowCopies;
  // buffer all input, then copy
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    (*it)->bottomBuffer();
  }
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    shallowCopies.insert( (*it)->shallowCopy() );
    (*it)->deepUnbuffer();
  }
  if ( shallowDelete){
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      delete (*it);
    }
  }
  // replace the input contents with the shallow copies
  _inputSeqs.clear();
  _inputSeqs.insert( shallowCopies.begin(), shallowCopies.end() );
}


bool AssemblyJobSubset::wasRun(){ return _wasRun; }

void AssemblyJobSubset::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}

void AssemblyJobSubset::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobSubset::remainingInputSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobSubset::remainingInputSeqs; seqs were already removed.");
  }
  outsideSet->insert( _remainSeqs.begin(), _remainSeqs.end() );
}

void AssemblyJobSubset::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}

void AssemblyJobSubset::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobSubset::newlyAssembledSeqs; job hasn't run yet.");
  }
  if ( _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobSubset::newlyAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}

void AssemblyJobSubset::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobSubset::allAssembledSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved or _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobSubset::allAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
  outsideSet->insert( _remainSeqs.begin(), _remainSeqs.end() );
}



void AssemblyJobSubset::OK(){
  if (! _wasRun ){
    if ( _novelSeqs.size() != 0 or _allAssembledSeqs.size() != 0 ){ //or _remainingSeqs.size() != 0 ){
      throw AssemblyException::LogicError("AssemblyJobSubset::OK; there should be no product seqs if job hasn't run.");
    }
    if ( _novelRemoved or _inputRemoved ){
      throw AssemblyException::LogicError("AssemblyJobSubset::OK; seqs can't have been removed if job hasn't run.");
    }
  }
}




bool AssemblyJobSubset::SortSsnByLength::operator() (ScoredSeqNormalized* ssnA, ScoredSeqNormalized* ssnB){
  return ssnA->size() < ssnB->size();
}



void AssemblyJobSubset::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){
    runJobMapping(threadedness);
    runJobCombine();
    _wasRun = true;
  }
}


void AssemblyJobSubset::runJobMapping(AssemblyJob::AssemblyThreadedness threadedness){

  // set up short-seq-indexed data structures
  _numInput = _inputSeqs.size();
  // fill these in during the step 1 loop
  _shortSeqArray = new ScoredSeqNormalized*[ _numInput+1 ];
  _shortIndexedAlSets = new vector<Alignment*>*[ _numInput+1 ];

  // STEP 1: create a BWT collection from all the sufficient-length sequences
  // and set up collections that will allow me to keep track of all the match
  // combinations even as I start consolidating sequences
  set<ScoredSeq*> inputForBwt; // cast differently for SSCBwt input
  long shortIndex = 0;
  for (set<ScoredSeq*>::iterator inputIt = _inputSeqs.begin(); inputIt != _inputSeqs.end(); ++inputIt){
    ScoredSeqNormalized* newSsn = new ScoredSeqNormalized( (*inputIt) );
    // this is the only truly necessary log-time addition
    inputForBwt.insert( newSsn );
    _shortSeqArray[shortIndex] = newSsn;
    _shortIndexedAlSets[shortIndex] = new vector<Alignment*>;
    shortIndex++;
  }
  // even though there is not a pre-specified min overlap, the alignments will all be run using
  // an input-specified min overlap, so penalize edge gaps is false
  ScoredSeqCollectionBwt* longSeqBwt = new ScoredSeqCollectionBwt( &inputForBwt, _minFractId, false );
  longSeqBwt->disableEdgeScaling();

  // STEP 2: query all the sufficient-length sequences to that collection
  // make a threading array of the seq test
  MatchSeqTest** matchTestArray = new MatchSeqTest*[ _numInput+1 ];
  long* minOvlArray = new long[ _numInput+1 ];
  ScoredSeq** recastSeqArray = new ScoredSeq*[ _numInput+1 ];
  for (long n = 0; n < _numInput; ++n){
    matchTestArray[n] = new MatchSeqTestMinLength( _shortSeqArray[n]->size() + 1);
    minOvlArray[n] = _shortSeqArray[n]->size();
    recastSeqArray[n] = _shortSeqArray[n];
  }

  // forward threadedness property to the search
  ScoredSeqCollectionBwt::AlignmentThreadedness alThread;
  if (threadedness == AssemblyJob::NOTTHREADED){ alThread = ScoredSeqCollectionBwt::notThreaded; }
  else { alThread = ScoredSeqCollectionBwt::threaded; }
  if (_strandedness == AssemblyJob::DOUBLESTRANDED){
    longSeqBwt->getBestMatches(_numInput, _shortIndexedAlSets, recastSeqArray, matchTestArray, minOvlArray,
			       alThread, ScoredSeqCollectionBwt::hardMinOvl);
  } else {
    longSeqBwt->getBestMatches(_numInput, _shortIndexedAlSets, recastSeqArray, '+', matchTestArray, minOvlArray,
			       alThread, ScoredSeqCollectionBwt::hardMinOvl);
  }

  for (long n = 0; n < _numInput; ++n){ delete matchTestArray[n]; }
  delete [] matchTestArray;
  delete [] minOvlArray;
  delete [] recastSeqArray;

  // i will use read counts in the SSNormalized to point to its index in an array of alignments
  // the index has to start at 1, the 0 index cannot be used because that means no counts
  // the counts in the long sequences will encode the index in the array
  vector<Alignment*>** tempLongAlignArray = new vector<Alignment*>*[ _numInput + 2 ];
  long longCountP1 = 1;

  // now make sure everything is organized properly in a NON-THREADED loop
  // DO NOT delete non-empty array elements because those sets are now values in _shortToAlignments
  // use the counts to figure out if a query has already been used as a long sequence: this cuts the
  // log-time searches in half when the match is non-redundant
  for (long n = 0; n < _numInput; ++n){
    // delete all elements at once, at the end
    if (_shortIndexedAlSets[n]->size() != 0){

      // now I must create or add to the long table entry
      for (vector<Alignment*>::iterator alIt = _shortIndexedAlSets[n]->begin(); alIt != _shortIndexedAlSets[n]->end(); ++alIt){
	Alignment* al = *alIt;
	ScoredSeqNormalized* longSeq = dynamic_cast<ScoredSeqNormalized*>( al->seqB() );

	// if there are no counts, then a new element must be added, otherwise the alignment just needs to be
	// added to the existing array
	if (longSeq->getCount() == 0){
	  longSeq->addCounts(longCountP1);
	  tempLongAlignArray[longCountP1] = new vector<Alignment*>;
	  tempLongAlignArray[longCountP1]->push_back( al );
	  ++longCountP1;
	} else {
	  tempLongAlignArray[ longSeq->getCount() ]->push_back( al );
	}
      }
    }
  }
  delete longSeqBwt;

  // before changing the counts, figure out the shortSeqArray index for each of the long seqs - also
  // put any sequence with a count into the _usedSeqs bin
  // the NESTED sequences here are correctly sorted for fast set addition
  long* tempLongToShort = new long[ longCountP1 ];
  for (long n = 0; n < _numInput; ++n){
    if ( _shortSeqArray[n]->getCount() > 0 ){
      tempLongToShort[ _shortSeqArray[n]->getCount() ] = n;
      _shortSeqArray[n]->resetCount();
      _shortSeqArray[n]->addCount();
    }
  }
  // this vector must initially be sorted in the same order as the temp arrays
  vector<ScoredSeqNormalized*> initialLongs;
  for (long n = 1; n < longCountP1; ++n){ initialLongs.push_back( _shortSeqArray[ tempLongToShort[n] ] ); }

  // add normalization counts to the normalized seqs based on their shortSeq alignment count
  for (long n = 0; n < _numInput; ++n){ _shortSeqArray[n]->addCounts( _shortIndexedAlSets[n]->size() ); }

  // add the nested seqs found in alignments to _usedSeqs - it is appropriately sorted
  set<ScoredSeq*>::iterator uIt = _usedSeqs.begin();
  for (long n = 0; n < _numInput; ++n){
    if ( _shortSeqArray[n]->getCount() > 0 ){
      uIt = _usedSeqs.insert( uIt, _shortSeqArray[n]->getNested() );
    }
  }


  // SECOND
  // now, filter through _inputSeqs for the retained seqs
  set<ScoredSeq*>::iterator inputIt = _inputSeqs.begin();
  set<ScoredSeq*>::iterator remainIt = _remainSeqs.begin();
  set<ScoredSeq*>::iterator usedIt = _usedSeqs.begin();
  set<ScoredSeq*>::iterator usedEnd = _usedSeqs.end();
  while (usedIt != usedEnd){
    if (*inputIt == *usedIt){ ++usedIt; }
    else{ remainIt = _remainSeqs.insert(remainIt, *inputIt); }
    ++inputIt;
  }
  _remainSeqs.insert(inputIt, _inputSeqs.end());

  // remove the count that was added above: only shortSeq counts are used for normalization
  // replace the count with the position in _initialLongs
  long* tempLongToCount = new long[ longCountP1 ];
  long lN = 0;
  for (vector<ScoredSeqNormalized*>::iterator it = initialLongs.begin(); it != initialLongs.end(); ++it){
    tempLongToCount[lN] = (*it)->getCount() - 1;
    (*it)->resetCount();
    (*it)->addCounts(lN);
    ++lN;
  }

  // now sort the longs, shortest to longest
  SortSsnByLength ssnSorter; // sort by seq length
  sort( initialLongs.begin(), initialLongs.end(), ssnSorter );

  // create the properly-ordered arrays and reset the normalization counts to their actual values
  _longAlignArray = new vector<Alignment*>*[ lN + 1 ];
  _longToShortIndexes = new long[ lN + 1 ];
  _initialLongs = new ScoredSeqNormalized*[ lN + 1 ];
  lN = 0;
  for (vector<ScoredSeqNormalized*>::iterator it = initialLongs.begin(); it != initialLongs.end(); ++it){
    long oldN = (*it)->getCount();
    _longAlignArray[lN] = tempLongAlignArray[oldN+1];
    _longToShortIndexes[lN] = tempLongToShort[oldN+1];
    (*it)->resetCount();
    (*it)->addCounts( tempLongToCount[oldN] );
    _initialLongs[lN] = *it;
    ++lN;
  }
  _longCount = lN;

  // clean up
  delete [] tempLongToShort;
  delete [] tempLongAlignArray;
  delete [] tempLongToCount;
}



// this cannot be safely parallelized, as far as i can see
void AssemblyJobSubset::runJobCombine(){

  // things that are created as intermediates
  vector<ScoredSeqNormalized*> transientCarriers;

  for (long olN = 0; olN < _longCount; ++olN){
    ScoredSeqNormalized* oldLong = _initialLongs[olN];
    vector<Alignment*>* longSet = _longAlignArray[olN];

    // combine the sequences to yield a new normalized seq
    ScoredSeqNormalized* newSeq = combineHelper(longSet, oldLong);

    if (oldLong == newSeq){ throw AssemblyException::LogicError("AJSubset: alignment set cannot be empty."); }

    // if the sequence does not exist in the short set, that means the end of the line has been
    // reached and the product can be added to the output set - i can know if it exists based on its
    // normalization count (long entries do not add to the count, so i only need to look for the short seq
    // if the count > 0)
    if ( oldLong->getCount() == 0 ){
      //_novelSeqs.insert(newSeq);
      _novelSeqs.insert(newSeq->shallowCopy());
      newSeq->deepDelete();
    } else {
      vector<Alignment*>* shortSet = _shortIndexedAlSets[ _longToShortIndexes[olN] ];
      transientCarriers.push_back(newSeq);
      for (vector<Alignment*>::iterator shortAlIt = shortSet->begin(); shortAlIt != shortSet->end(); ++shortAlIt){
	Alignment* shortAl = (*shortAlIt);
	ScoredSeqNormalized* localLong = dynamic_cast<ScoredSeqNormalized*>( shortAl->seqB() );
	shortAl->seqReplace(newSeq, localLong);
      }
    }
    delete longSet;
  }

  // short-to-long was not cleared during the collapse
  for (long n = 0; n < _numInput; ++n){
    delete _shortIndexedAlSets[n];
    delete _shortSeqArray[n];
  }
  delete [] _shortIndexedAlSets;
  delete [] _shortSeqArray;
  delete [] _initialLongs;
  delete [] _longToShortIndexes;
  delete [] _longAlignArray;

  // deep-delete the transient sequences
  for (vector<ScoredSeqNormalized*>::iterator trashIt = transientCarriers.begin(); trashIt != transientCarriers.end(); ++trashIt){
    (*trashIt)->deepDelete();
  }
}



// alSet members all have the same seqB
// seqB is longer than seqA
// this function deletes the alignments in alSet
ScoredSeqNormalized* AssemblyJobSubset::combineHelper(vector<Alignment*>* alSet, ScoredSeqNormalized* oldSeq){
  long numAl = alSet->size();
  if (numAl == 0){ return oldSeq; }
  else if (numAl == 1){
    Alignment* al = *(alSet->begin());
    ScoredSeq* combined = al->combine();
    ScoredSeqNormalized* carrier;
    switch (al->orientation()){
    case '+': carrier = new ScoredSeqNormalized( combined ); break;
    case '-': carrier = new ScoredSeqNormalized( combined->flipCopy() ); delete combined; break;
    default: throw AssemblyException::LogicError("AJSubset bad sense this won't happen");
    }
    //ScoredSeqNormalized* carrier = new ScoredSeqNormalized( al->combine() );
    carrier->addCounts( oldSeq->getCount() );
    delete al;
    return carrier;
  } else {
    long numAlM1 = numAl - 1;
    Alignment** alArray = new Alignment*[ numAl ];
    vector<Alignment*>::iterator alIt = alSet->begin();
    ScoredSeq* longSeq = (*alIt)->seqB();
    long sizeB = longSeq->size();
    AlignmentUngapped* seedAl;
    for (long alN = 0; alN <= numAlM1; ++alN){
      AlignmentUngapped* newAl;
      Alignment* oldAl = *alIt;
      switch (oldAl->orientation()){
      case '+': newAl = new AlignmentUngapped(longSeq, oldAl->seqA(), '+',
					      0 - oldAl->getConstantOffset()); break;
      case '-': newAl = new AlignmentUngapped(longSeq, oldAl->seqA(), '-',
					      sizeB - oldAl->seqA()->size() + oldAl->getConstantOffset()); break;
      default: throw AssemblyException::LogicError("AJSubset bad sense this won't happen");
      }
      if (alN == numAlM1){ seedAl = newAl; }
      else { alArray[alN] = newAl; }
      delete oldAl;
      ++alIt;
    }

    // collapse the alignments
    ScoredSeqNormalized* carrier;
    if (numAl == 1){ carrier = new ScoredSeqNormalized( seedAl->combine() ); }
    else { carrier = new ScoredSeqNormalized( seedAl->multiCombine(numAlM1, alArray) ); }
    carrier->addCounts( oldSeq->getCount() );

    // clean up
    delete seedAl;
    for (long n = 0; n < numAlM1; ++n){ delete alArray[n]; }
    delete [] alArray;

    return carrier;
  }
}

void AssemblyJobSubset::printAlignNw(Alignment * al){
  if ( al->isNull() ){
    throw AssemblyException::ArgError("in printAlignNw, alignment may not be null.");
  }
  ScoredSeq * seqA = al->seqA();
  ScoredSeq * seqB = al->seqB();
  char sense = al->orientation();
  long posA = 0;
  long posB = 0;
  long maxA = seqA->size();
  long maxB = seqB->size();

  string stringA = "";
  string stringB = "";
  string matchLines = "";

  while ( posA < maxA or posB < maxB ){

    if ( posA < maxA and al->isGapped(posA,seqA) ){
      stringA += seqA->nucAtPosPlus(posA);
      stringB += '-';
      matchLines += ' ';
      posA++;
    } else if ( posB < maxB and al->isGapped(posB,seqB) ){
      stringA += '-';
      stringB += seqB->nucAtPosition(posB,sense);
      matchLines += ' ';
      posB++;
    } else {
      char nucA = seqA->nucAtPosPlus(posA);
      char nucB = seqB->nucAtPosition(posB,sense);
      stringA += nucA;
      stringB += nucB;
      if ( nucA == nucB ){ matchLines += '|'; }
      else { matchLines += ' '; }
      posA++;
      posB++;
    }

  }

  cout << stringA << endl;
  cout << matchLines << endl;
  cout << stringB << endl;
}




#endif
