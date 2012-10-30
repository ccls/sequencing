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

#ifndef ASSEMBLYJOBPERFECT_CPP
#define ASSEMBLYJOBPERFECT_CPP


#include "AssemblyJobPerfect.h"

#include "AssemblyJob.h"
#include "AssemblyException.h"
#include <queue>
#include <string>
#include <iostream>
using namespace::std;

AssemblyJobPerfect::AssemblyJobPerfect(){}
AssemblyJobPerfect::AssemblyJobPerfect(set<ScoredSeq*>* seqs,
				       AssemblyJob::AssemblyStrandedness strandedness) :
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _strandedness(strandedness){
  _inputSeqs.insert(seqs->begin(),seqs->end());
}


AssemblyJobPerfect::~AssemblyJobPerfect(){
}



bool AssemblyJobPerfect::wasRun(){ return _wasRun; }

void AssemblyJobPerfect::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}

void AssemblyJobPerfect::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobP::remainingInputSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobP::remainingInputSeqs; seqs were already removed.");
  }
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
}

void AssemblyJobPerfect::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}

void AssemblyJobPerfect::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobP::newlyAssembledSeqs; job hasn't run yet.");
  }
  if ( _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobP::newlyAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}

void AssemblyJobPerfect::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobP::allAssembledSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved or _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobP::allAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}



void AssemblyJobPerfect::OK(){
  if (! _wasRun ){
    if ( _novelSeqs.size() != 0 or _allAssembledSeqs.size() != 0 or _remainingSeqs.size() != 0 ){
      throw AssemblyException::LogicError("AssemblyJobPerfect::OK; there should be no product seqs if job hasn't run.");
    }
    if ( _novelRemoved or _inputRemoved ){
      throw AssemblyException::LogicError("AssemblyJobPerfect::OK; seqs can't have been removed if job hasn't run.");
    }
  }
}







void AssemblyJobPerfect::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJNull::makeShallow can't be called after the job was run.");
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

void AssemblyJobPerfect::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){

    // get all of the lengths; only run the collapse for sets of similar-length sequences
    // STEP 1: organize input data by sequence length
    // note: the orders of seqs in the vectors are set-sorted - i am using vectors
    // here so i don't have to keep track of iterators
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

    // organize into arrays to facilitate threading
    // NOTE: the vectors are still sorted appropriately for efficient set addition
    vector<ScoredSeq*>** indexToSeqSet = new vector<ScoredSeq*>*[ lenToSeqSet.size() + 1 ];
    // this will reduce the complexity in the end of looking for input sequences
    set<ScoredSeq*>* indexToRetained = new set<ScoredSeq*>[ lenToSeqSet.size() + 1 ];
    set<ScoredSeq*>* indexToUsed = new set<ScoredSeq*>[ lenToSeqSet.size() + 1 ];
    set<ScoredSeq*>* indexToNovel = new set<ScoredSeq*>[ lenToSeqSet.size() + 1 ];
    long numLengthSets = 0;
    for (map<long, vector<ScoredSeq*>* >::iterator lenToSeqsIt = lenToSeqSet.begin(); lenToSeqsIt != lenToSeqSet.end(); ++lenToSeqsIt){
      indexToSeqSet[numLengthSets] = lenToSeqsIt->second;
      numLengthSets++;
    }

    // do the collapsing
    if (threadedness == AssemblyJob::THREADED){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numLengthSets; ++n){
	subAssembly(indexToSeqSet[n], &indexToRetained[n], &indexToUsed[n], &indexToNovel[n]);
      }
    } else {
      for (long n = 0; n < numLengthSets; ++n){
	subAssembly(indexToSeqSet[n], &indexToRetained[n], &indexToUsed[n], &indexToNovel[n]);
      }
    }

    // STEP 3: clean up
    // re-consolidate and sort the sequences
    for (long n = 0; n < numLengthSets; ++n){
      _retainedSeqs.insert(indexToRetained[n].begin(), indexToRetained[n].end());
      _usedSeqs.insert(indexToUsed[n].begin(), indexToUsed[n].end());
      _novelSeqs.insert(indexToNovel[n].begin(), indexToNovel[n].end());
      delete indexToSeqSet[n];
    }
    delete [] indexToSeqSet;
    delete [] indexToRetained;
    delete [] indexToUsed;
    delete [] indexToNovel;
    _wasRun = true;
  }
  //OK();
}




AssemblyJobPerfect::SameSizeSeqCombiner::SameSizeSeqCombiner(ScoredSeq* original) :
  _original(original),
  _isOriginal(true),
  _seqGotten(false)
{}
AssemblyJobPerfect::SameSizeSeqCombiner::~SameSizeSeqCombiner(){
  if (! _seqGotten){ throw AssemblyException::LogicError("AJPerfect::SameSizeSeqCombiner shouldn't be deleted if seq wasn't gotten."); }
}
ScoredSeq* AssemblyJobPerfect::SameSizeSeqCombiner::getSeq(){
  _seqGotten = true;
  if ( _isOriginal ){ return _original; }
  else { return ScoredSeq::repExposedSeq( _original->getSeq('+'), _scores, _links, _original->size() ); }
}
bool AssemblyJobPerfect::SameSizeSeqCombiner::isOriginal(){ return _isOriginal; }
ScoredSeq* AssemblyJobPerfect::SameSizeSeqCombiner::getOriginal(){ return _original; }

void AssemblyJobPerfect::SameSizeSeqCombiner::addSeq(ScoredSeq* seq, char sense){
  // make sure that the carrier arrays exist
  if ( _isOriginal ){
    _isOriginal = false;
    _scores = _original->getScores('+');
    _links = _original->getLinks('+');
  }
  // add to the scores
  float* tempScores = seq->getScores(sense);
  long seqSize = _original->size();
  for (long n = 0; n < seqSize; ++n){ _scores[n] += tempScores[n]; }
  delete [] tempScores;
  // add to the links
  float* tempLinks = seq->getLinks(sense);
  long seqSizeM1 = seqSize - 1;
  for (long n = 0; n < seqSizeM1; ++n){ _links[n] += tempLinks[n]; }
  delete [] tempLinks;
}



void AssemblyJobPerfect::subAssembly(vector<ScoredSeq*>* monoLenSeqSet, set<ScoredSeq*>* retainedSeqs, 
				     set<ScoredSeq*>* usedSeqs, set<ScoredSeq*>* novelSeqs){

  long setSize = monoLenSeqSet->size();
  if (setSize == 1){ retainedSeqs->insert( *(monoLenSeqSet->begin()) ); }
  else if (setSize > 1){

    // retains the set-appropriate order for retained seqs
    vector<SameSizeSeqCombiner*> orderedCombiners;
    // the condensed set, indexed by seqeuence
    map<string,SameSizeSeqCombiner*> stringToSeq;

    // accelerates addition to the old seq sets
    set<ScoredSeq*>::iterator rIt = retainedSeqs->begin();
    set<ScoredSeq*>::iterator uIt = usedSeqs->begin();

    for (vector<ScoredSeq*>::iterator inputIt = monoLenSeqSet->begin(); inputIt != monoLenSeqSet->end(); ++inputIt){

      // try for a sense match
      char* keySeq = (*inputIt)->getSeq('+');
      map<string,SameSizeSeqCombiner*>::iterator matchIt = stringToSeq.find(keySeq);

      if ( matchIt!=stringToSeq.end() ){
	matchIt->second->addSeq( *inputIt, '+' );
	uIt = usedSeqs->insert(uIt, *inputIt);
      } else if (_strandedness == AssemblyJob::DOUBLESTRANDED){
	char* keySeqRc = (*inputIt)->getSeq('-');
	matchIt = stringToSeq.find(keySeqRc);
	if ( matchIt != stringToSeq.end() ){
	  matchIt->second->addSeq( *inputIt, '-' );
	  uIt = usedSeqs->insert(uIt, *inputIt);
	}
	delete [] keySeqRc;
      }

      // if no match was found, add the seq
      if (matchIt == stringToSeq.end()){
	SameSizeSeqCombiner* combiner = new SameSizeSeqCombiner(*inputIt);
	stringToSeq.insert( pair<string,SameSizeSeqCombiner*> (keySeq, combiner) );
	orderedCombiners.push_back( combiner );
      }
      delete [] keySeq;
    }

    // now add the leftovers and novel sequences to the appropriate bins
    // retained are added quickly; novel require log time
    set<ScoredSeq*> tempUsed;
    set<ScoredSeq*>::iterator tuIt = tempUsed.begin();
    for (vector<SameSizeSeqCombiner*>::iterator combIt = orderedCombiners.begin(); combIt != orderedCombiners.end(); ++combIt){
      if ( (*combIt)->isOriginal() ){ rIt = retainedSeqs->insert(rIt, (*combIt)->getSeq()); }
      else {
	novelSeqs->insert( (*combIt)->getSeq() );
	tuIt = tempUsed.insert(tuIt, (*combIt)->getOriginal());
      }
      delete *combIt;
    }
    usedSeqs->insert(tempUsed.begin(), tempUsed.end());
  }
}



#endif
