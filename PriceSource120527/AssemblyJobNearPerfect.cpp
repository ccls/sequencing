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

#ifndef ASSEMBLYJOBNEARPERFECT_CPP
#define ASSEMBLYJOBNEARPERFECT_CPP


#include "AssemblyJobNearPerfect.h"

#include "AssemblyException.h"
#include "AlignmentFull.h"
#include <queue>
#include <string>
#include <iostream>
#include <math.h>
using namespace::std;

AssemblyJobNearPerfect::AssemblyJobNearPerfect(){}
AssemblyJobNearPerfect::AssemblyJobNearPerfect(set<ScoredSeq*>* seqs,
					       float minFractId,
					       AssemblyJob::AssemblyStrandedness strandedness) :
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _strandedness(strandedness),
  _maxFractMis(1.0 - minFractId){
  _inputSeqs.insert(seqs->begin(),seqs->end());
  // true if mismatches <= maxMisNum
  _mismatchCountAsm = new AlignmentScoreMatrix(0,1,0,0); // mismatch counts as 1
  _matchCountAsm = new AlignmentScoreMatrix(1,0,0,0); // mismatch counts as 1
}


AssemblyJobNearPerfect::~AssemblyJobNearPerfect(){
  delete _mismatchCountAsm;
  delete _matchCountAsm;
}


void AssemblyJobNearPerfect::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJNP::makeShallow can't be called after the job was run.");
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


bool AssemblyJobNearPerfect::wasRun(){ return _wasRun; }

void AssemblyJobNearPerfect::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}

void AssemblyJobNearPerfect::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobNP::remainingInputSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobNP::remainingInputSeqs; seqs were already removed.");
  }
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
}

void AssemblyJobNearPerfect::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}

void AssemblyJobNearPerfect::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobNP::newlyAssembledSeqs; job hasn't run yet.");
  }
  if ( _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobNP::newlyAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}

void AssemblyJobNearPerfect::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobNP::allAssembledSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved or _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobNP::allAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}



void AssemblyJobNearPerfect::OK(){
}







void AssemblyJobNearPerfect::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){

    // STEP 1: organize input data by sequence length
    // these vectors are ordered appropriately for set addition - ordered by the original input sequence!
    map<long, vector<ScoredSeq*>* > lenToSeqSet;
    for (set<ScoredSeq*>::iterator inputIt = _inputSeqs.begin(); inputIt != _inputSeqs.end(); ++inputIt){
      ScoredSeq* inSeq = (*inputIt);
      map<long, vector<ScoredSeq*>* >::iterator lenToSeqsIt = lenToSeqSet.find( inSeq->size() );
      if ( lenToSeqsIt == lenToSeqSet.end() ){
	vector<ScoredSeq*>* newEntry = new vector<ScoredSeq*>;
	newEntry->push_back( inSeq );
	lenToSeqSet.insert( pair<long, vector<ScoredSeq*>* > (inSeq->size(), newEntry) );
      } else {
	lenToSeqsIt->second->push_back( inSeq );
      }
    }

    // STEP 2: run each set through a subroutine that consolidates its contents based on near-matches

    // organize into arrays to facilitate threading
    long indexToSeqLen[ lenToSeqSet.size() + 1 ];
    vector<ScoredSeq*>** indexToSeqSet = new vector<ScoredSeq*>*[ lenToSeqSet.size() + 1 ];
    // this will reduce the complexity in the end of looking for input sequences
    set<ScoredSeq*>* indexToSeqInput = new set<ScoredSeq*>[ lenToSeqSet.size() + 1 ];
    long numLengthSets = 0;
    for (map<long, vector<ScoredSeq*>* >::iterator lenToSeqsIt = lenToSeqSet.begin(); lenToSeqsIt != lenToSeqSet.end(); ++lenToSeqsIt){
      indexToSeqLen[numLengthSets] = lenToSeqsIt->first;
      indexToSeqSet[numLengthSets] = lenToSeqsIt->second;
      set<ScoredSeq*>::iterator isiIt = indexToSeqInput[numLengthSets].begin();
      for (vector<ScoredSeq*>::iterator vIt = lenToSeqsIt->second->begin(); vIt != lenToSeqsIt->second->end(); ++vIt){
	isiIt = indexToSeqInput[numLengthSets].insert(isiIt, *vIt);
      }
      //indexToSeqInput[numLengthSets].insert(indexToSeqSet[numLengthSets]->begin(), indexToSeqSet[numLengthSets]->end());
      numLengthSets++;
    }

    // run the consolidation threaded if appropriate
    if (threadedness == AssemblyJob::THREADED){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numLengthSets; ++n){
	subAssembly( indexToSeqLen[n], indexToSeqSet[n] );
      }
    } else {
      for (long n = 0; n < numLengthSets; ++n){
	subAssembly( indexToSeqLen[n], indexToSeqSet[n] );
      }
    }

    // STEP 3: clean up
    // re-consolidate and sort the sequences
    if (_usedSeqs.size() != 0){ throw AssemblyException::LogicError("AJPerfect::runJob, _novelSeqs not empty near end of job."); }
    if (_usedSeqs.size() != 0){ throw AssemblyException::LogicError("AJPerfect::runJob, _retainedSeqs not empty near end of job."); }
    for (long n = 0; n < numLengthSets; ++n){
      // temporary data structures to make addition time constant
      set<ScoredSeq*> goToNovel;
      set<ScoredSeq*>::iterator nsIt = goToNovel.begin();
      set<ScoredSeq*> goToRetained;;
      set<ScoredSeq*>::iterator rsIt = goToRetained.begin();
      for (vector<ScoredSeq*>::iterator it = indexToSeqSet[n]->begin(); it != indexToSeqSet[n]->end(); ++it){
	if (indexToSeqInput[n].count(*it) == 0){ nsIt = goToNovel.insert(nsIt, *it); }
	else { rsIt = goToRetained.insert(rsIt, *it); }
      }
      _novelSeqs.insert(goToNovel.begin(), goToNovel.end());
      _retainedSeqs.insert(goToRetained.begin(), goToRetained.end());
      delete indexToSeqSet[n];
    }
    delete [] indexToSeqSet;
    delete [] indexToSeqInput;

    // figure out which of the input seqs was used
    // here, the seqs ARE sorted for efficient set addition
    if (_usedSeqs.size() != 0){ throw AssemblyException::LogicError("AJPerfect::runJob, _usedSeqs not empty near end of job."); }
    set<ScoredSeq*>::iterator usIt = _usedSeqs.begin();
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      if (_retainedSeqs.count(*it) == 0){ usIt = _usedSeqs.insert(usIt, *it); }
    }

    _wasRun = true;
  }
  //OK();
}


AssemblyJobNearPerfect::WindowStartAndLenCarrier::WindowStartAndLenCarrier(){}
AssemblyJobNearPerfect::WindowStartAndLenCarrier::WindowStartAndLenCarrier(long seqLen, long maxMisNum){
  // calculate the number of windows
  // scale the number of windows with seqLen (also considers best-case scenario with maxMisNum)
  long maxNumWindows = maxMisNum + 1;
  long baseCase = 50;
  long minNumWindows;
  if (seqLen < baseCase){ minNumWindows = 2; }
  else{ minNumWindows = long( log10( double(seqLen) / baseCase ) / log10( 2 ) ) + 2; }
  if (maxNumWindows < minNumWindows){ _numWindows = maxNumWindows; }
  else { _numWindows = minNumWindows; }

  // find the window dimensions
  _winStart = new long[_numWindows];
  _winLen = new long[_numWindows];
  // delineate the windows as evenly as possible
  _winStart[0] = 0;
  for (int n = 1; n < _numWindows; ++n){
    _winStart[n] = n * seqLen / _numWindows;
    _winLen[n-1] = _winStart[n] - _winStart[n-1];
  }
  _winLen[_numWindows-1] = seqLen - _winStart[_numWindows-1];
}
AssemblyJobNearPerfect::WindowStartAndLenCarrier::~WindowStartAndLenCarrier(){
  delete [] _winStart;
  delete [] _winLen;
}

// THE NEW...
void AssemblyJobNearPerfect::subAssembly2(long seqLen, vector<ScoredSeq*>* monoLenSeqSet){
  /*
  long maxMisNum = long(float(seqLen) * _maxFractMis);
  WindowStartAndLenCarrier* wslc = new WindowStartAndLenCarrier(seqLen, maxMisNum);
  vector< vector<SeqCarrier*>* > wToStrIndToSeqs[ wslc->_numWindows + 1 ];
  map<string, long> wToStrToIndex[ wslc->_numWindows + 1 ];

  for (vector<ScoredSeq*>::iterator seqIt = monoLenSeqSet->begin(); seqIt != monoLenSeqSet->end(); ++seqIt){

    // get the sub-seq strings
    char** subSeqArray = new char*[ wslc->_numWindows + 1 ];
    map<string, long>::iterator stiIt[ wslc->_numWindows + 1 ];

    // find the index of each window's substring in the vector
    for (int n = 0; n < wslc->_numWindows; ++n){
      subSeqArray[n] = inputSeq->getSubseq(wslc->_winStart[n], wslc->_winLen[n], sense);
      map<string, long>::iterator stiIt = wToStrToIndex[n].find( subSeqArray[n] );
      if (stiIt == wToStrToIndex[n].end()){ subSeqIndexes[n] = wToStrToIndex[n].size(); }
      else { subSeqIndexes[n] = stiIt->second; }
    }

    // get the potential matches - as they are found, remove them from
  }
  */
}



// THE OLD...
void AssemblyJobNearPerfect::subAssembly(long seqLen, vector<ScoredSeq*>* monoLenSeqSet){
  // it only needs to modify the set if there are multiple entries, otherwise the seq will not be changed
  if (monoLenSeqSet->size() > 1){
    ScoredSeqNormalized** withWrappers = new ScoredSeqNormalized*[ monoLenSeqSet->size() + 1 ];
    long ssnIndex = 0;
    for (vector<ScoredSeq*>::iterator it = monoLenSeqSet->begin(); it != monoLenSeqSet->end(); ++it){
      withWrappers[ssnIndex] = new ScoredSeqNormalized( *it );
      ssnIndex++;
    }
    withWrappers[ssnIndex] = NULL;

    subAssembly( seqLen, withWrappers );

    monoLenSeqSet->clear();
    ssnIndex = 0;
    while ( withWrappers[ssnIndex] != NULL ){
      monoLenSeqSet->push_back( withWrappers[ssnIndex]->getNested() );
      delete withWrappers[ssnIndex];
      ssnIndex++;
    }
    delete [] withWrappers;
  }
}


// MODIFIES: monoLenSeqInput
void AssemblyJobNearPerfect::subAssembly(long seqLen, ScoredSeqNormalized** monoLenSeqInput){
  long maxMisNum = long(float(seqLen) * _maxFractMis);


  // local so that it can be emptied
  // these normalized wrappers are NOT sorted for fast set addition
  set<ScoredSeqNormalized*> monoLenSeqSet;
  long ssnCount = 0;
  while ( monoLenSeqInput[ssnCount] != NULL ){
    monoLenSeqSet.insert( monoLenSeqInput[ssnCount] );
    ssnCount++;
  }
  
  WindowStartAndLenCarrier* wslc = new WindowStartAndLenCarrier(seqLen, maxMisNum);
  map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs = new map<string, set<ScoredSeqNormalized*>* >[ wslc->_numWindows ];

  for (long ssnIndex = 0; ssnIndex < ssnCount; ++ssnIndex){
    ScoredSeqNormalized* query = monoLenSeqInput[ssnIndex];
    monoLenSeqSet.erase( query );

    // get all the legit matches for each strand of the query
    map<ScoredSeqNormalized*,AlignmentFull*> legitMatches; // key is seqB: seqA is the query seq; seqB is the found match
    // NOTE: this may NOT contain both (+) and (-) orientation matches to a given seqB

    if (_strandedness == AssemblyJob::DOUBLESTRANDED){
      getLegitMatchesHelper( query, '+', wslc, maxMisNum, windowToSubstrToSeqs, &legitMatches );
      getLegitMatchesHelper( query, '-', wslc, maxMisNum, windowToSubstrToSeqs, &legitMatches );
    } else {
      getLegitMatchesHelper( query, '+', wslc, maxMisNum, windowToSubstrToSeqs, &legitMatches );
    }

    // define the set of sequences to be removed and add the new sequences to the _novel set
    // NOTE: the query is not yet in the search structure
    set<ScoredSeqNormalized*> seqsToRemove;
    // either the query or the newly-consolidated sequence will be added
    ScoredSeqNormalized* consolidatedSeq;
    if ( legitMatches.size() == 1 ){
      map<ScoredSeqNormalized*,AlignmentFull*>::iterator legitIt = legitMatches.begin();
      seqsToRemove.insert( legitIt->first );
      consolidatedSeq = new ScoredSeqNormalized( legitIt->second->combine() );
      delete legitIt->second;
      consolidatedSeq->buffer();
      ScoredSeq* innerQuery = query->getNested();
      delete query;
      if (_inputSeqs.count(innerQuery) == 0){ delete innerQuery; }
    } else if ( legitMatches.size() > 1 ){
      // set up extra alignments
      long extraAlCount = legitMatches.size() - 1;
      AlignmentFull** extraAlArray = new AlignmentFull*[ extraAlCount ];
      // skip the first alignment; that will be the seed alignment for combining the others

      map<ScoredSeqNormalized*,AlignmentFull*>::iterator legitIt = legitMatches.begin();
      AlignmentFull* seedAlignment = legitIt->second;
      seqsToRemove.insert( legitIt->first );
      long extraAlIndex = 0;
      ++legitIt;
      while (legitIt != legitMatches.end()){
	seqsToRemove.insert( legitIt->first );
	extraAlArray[extraAlIndex] = legitIt->second;
	++extraAlIndex;
	++legitIt;
      }
      //use the first alignment for the combine
      consolidatedSeq = new ScoredSeqNormalized( seedAlignment->multiCombine(extraAlCount, extraAlArray) );
      delete seedAlignment;
      consolidatedSeq->buffer();
      ScoredSeq* innerQuery = query->getNested();
      delete query;
      if (_inputSeqs.count(innerQuery) == 0){ delete innerQuery; }
      // clean up extra alignments
      for (long n = 0; n < extraAlCount; ++n){ delete extraAlArray[n]; }
      delete [] extraAlArray;
    } else {
      consolidatedSeq = query;
    }

    // insert the new seq (there will be only one by definition)
    monoLenSeqSet.insert( consolidatedSeq );
    insertSeqHelper( consolidatedSeq, wslc, windowToSubstrToSeqs );

    // delete the old seqs if that is allowed
    for (set<ScoredSeqNormalized*>::iterator oldSeqIt = seqsToRemove.begin(); oldSeqIt != seqsToRemove.end(); ++oldSeqIt){
      monoLenSeqSet.erase( *oldSeqIt );
      removeSeqHelper( (*oldSeqIt), wslc, windowToSubstrToSeqs );
    }
  }

  for (long n = 0; n < wslc->_numWindows; ++n){
    for (map<string, set<ScoredSeqNormalized*>* >::iterator it = windowToSubstrToSeqs[n].begin(); it != windowToSubstrToSeqs[n].end(); ++it){
      delete it->second;
    }
  }
  delete wslc;
  delete [] windowToSubstrToSeqs;

  // modify the input set to just have the now legit contig seqs
  {
    long ssnCount = 0;
    for (set<ScoredSeqNormalized*>::iterator it = monoLenSeqSet.begin(); it != monoLenSeqSet.end(); ++it){
      monoLenSeqInput[ssnCount] = *it;
      ssnCount++;
    }
    monoLenSeqInput[ssnCount] = NULL;
  }
}


// add the combined seq to the front/back special map collections
void AssemblyJobNearPerfect::insertSeqHelper(ScoredSeqNormalized* newContig,
					     WindowStartAndLenCarrier* wslc,
					     map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs){
  for (int n = 0; n < wslc->_numWindows; ++n){
    char* subSeq = newContig->getSubseq(wslc->_winStart[n], wslc->_winLen[n], '+');
    map<string,set<ScoredSeqNormalized*>* >::iterator matchIt = windowToSubstrToSeqs[n].find(subSeq);
    if ( matchIt == windowToSubstrToSeqs[n].end() ){
      set<ScoredSeqNormalized*>* newSeqSet = new set<ScoredSeqNormalized*>;
      newSeqSet->insert( newContig );
      windowToSubstrToSeqs[n].insert( pair<string, set<ScoredSeqNormalized*>* > (subSeq, newSeqSet) );
    } else {
      matchIt->second->insert( newContig );
    }
    delete [] subSeq;
  }
}


// gets rid of the seq from the provided collections, moves it around the class-level
// collections as appropriate, deletes if that is ok.
void AssemblyJobNearPerfect::removeSeqHelper(ScoredSeqNormalized* oldSeq,
					     WindowStartAndLenCarrier* wslc,
					     map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs){
  // remove the sequence from the front/back maps
  // make sure that the seq has been put in the table; if it has, remove it
  //char* oldSeqString = oldSeq->getSeq('+');
  for (int n = 0; n < wslc->_numWindows; ++n){
    char* subSeq = oldSeq->getSubseq(wslc->_winStart[n], wslc->_winLen[n], '+');
    map<string,set<ScoredSeqNormalized*>* >::iterator matchIt = windowToSubstrToSeqs[n].find(subSeq);
    if ( matchIt != windowToSubstrToSeqs[n].end() ){
      matchIt->second->erase(oldSeq);
      if ( matchIt->second->size() == 0 ){
	delete matchIt->second;
	windowToSubstrToSeqs[n].erase( matchIt );
      }
    }
    delete [] subSeq;
  }

  // delete the seq and its carrier
  ScoredSeq* oldInnerSeq = oldSeq->getNested();
  delete oldSeq;
  if (_inputSeqs.count(oldInnerSeq) == 0){ delete oldInnerSeq; }
}


// gets the matches satisfying the requirements as efficiently as possible
// using the front/back data structures and ungapped alignment technique
// does not delete anything, does not modify class data structures
// MODIFIES: legitMatches
void AssemblyJobNearPerfect::getLegitMatchesHelper(ScoredSeq* inputSeq, char sense,
						   WindowStartAndLenCarrier* wslc, long maxMisNum,
						   map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs,
						   map<ScoredSeqNormalized*,AlignmentFull*>* legitMatches){
  // establish the query sequences
  long seqLen = inputSeq->size();
  //string* subSeqArray = new string[ wslc->_numWindows ];
  char** subSeqArray = new char*[ wslc->_numWindows + 1 ];
  for (int n = 0; n < wslc->_numWindows; ++n){
    subSeqArray[n] = inputSeq->getSubseq(wslc->_winStart[n], wslc->_winLen[n], sense);
  }

  set<ScoredSeqNormalized*> allCandidateMatches;
  long allMatchSum = long( pow(double(2), double(wslc->_numWindows) ) ) - 1;
  for (long n = 0; n < wslc->_numWindows; ++n){
    map<string,set<ScoredSeqNormalized*>* >::iterator matchIt = windowToSubstrToSeqs[n].find( subSeqArray[n] );
    if ( matchIt != windowToSubstrToSeqs[n].end() ){
      for (set<ScoredSeqNormalized*>::iterator it = matchIt->second->begin(); it != matchIt->second->end(); ++it){
	(*it)->addCounts( long( pow(double(2), double(n)) ) );
	allCandidateMatches.insert( (*it) );
      }
    }
  }

  // check for legitimacy of candidate matches
  for (set<ScoredSeqNormalized*>::iterator candIt = allCandidateMatches.begin(); candIt != allCandidateMatches.end(); ++candIt){

    ScoredSeqNormalized* candSeq = (*candIt);
    bool legitMatchFound = false;

    // figure out if there is a legit match
    long candSeqCount = candSeq->getCount();
    if ( candSeqCount == allMatchSum ){ legitMatchFound = true; }
    else {
      // see if the number of windows with matches was sufficient for a quality match
      long numGoodWindows = 0;
      bool* isWindowGood = new bool[wslc->_numWindows];
      for (int n = 0; n < wslc->_numWindows; ++n){
	long binaryCount = ( candSeqCount % long( pow(double(2), double(n+1)) ) ) / long( pow(double(2), double(n) ) );
	numGoodWindows += binaryCount;
	isWindowGood[n] = bool(binaryCount);
      }
      if (wslc->_numWindows - numGoodWindows <= maxMisNum){
	// do the real comparison so that it will stop doing comparisons if the errors have over-accumulated
	legitMatchFound = true;
	long numRemainingMisses = maxMisNum;
	long n = 0;
	while (legitMatchFound and n < wslc->_numWindows){
	  // don't re-check perfect windows
	  if (! isWindowGood[n] ){
	    char* compareSeq = candSeq->getSubseq(wslc->_winStart[n], wslc->_winLen[n], '+');
	    long numMis = compareSeqHelper(subSeqArray[n], compareSeq, wslc->_winLen[n], sense, numRemainingMisses);
	    delete [] compareSeq;
	    if (numMis > numRemainingMisses){ legitMatchFound = false; }
	    else { numRemainingMisses -= numMis; }
	  }
	  n++;
	}
      }
      delete [] isWindowGood;
    }

    // so that the opposite-strand matches can be found clearly
    candSeq->resetCount();

    // if there is a legit match, add the alignment to the collection of legit matches
    if (legitMatchFound){
      AlignmentFull* legitMatch = new AlignmentFull( inputSeq, candSeq, sense );
      // check if there is already a legit match to the other strand
      map<ScoredSeqNormalized*,AlignmentFull*>::iterator foundIt = legitMatches->find( candSeq );
      if ( foundIt == legitMatches->end() ){
	legitMatches->insert( pair<ScoredSeqNormalized*,AlignmentFull*> (candSeq,legitMatch) );
      } else if ( foundIt->second->score( _matchCountAsm, false ) < legitMatch->score( _matchCountAsm, false ) ){
	delete foundIt->second; // delete the old alignment, replace with the better alignment
	foundIt->second = legitMatch;
      }
    }
  }

  for (int n = 0; n < wslc->_numWindows; ++n){ delete [] subSeqArray[n]; }
  delete [] subSeqArray;
}


inline long AssemblyJobNearPerfect::compareSeqHelper(char* seqA, char* seqB, long seqSize, char sense, long maxMisNum){
  long numMissed = 0;
  long pos = 0;
  while ( numMissed <= maxMisNum and pos < seqSize ){
    numMissed += _mismatchCountAsm->match( seqA[pos], seqB[pos] );
    ++pos;
  }
  return numMissed;
}


#endif
