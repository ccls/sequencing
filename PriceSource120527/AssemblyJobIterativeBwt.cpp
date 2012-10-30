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


#ifndef ASSEMBLYJOBITERATIVEBWT_CPP
#define ASSEMBLYJOBITERATIVEBWT_CPP


#include "AssemblyJobIterativeBwt.h"
#include "AssemblyException.h"
#include "MatchAlignTest.h"

#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace::std;

AssemblyJobIterativeBwt::AssemblyJobIterativeBwt(){}
// THESE CONSTRUCTORS SHOULD BE CHANGED

// HERE, NON-THREADED IS DEFAULT
AssemblyJobIterativeBwt::AssemblyJobIterativeBwt(set<ScoredSeq*>* seqs,
						 AssemblyJob::AssemblyStrandedness strandedness,
						 AssemblyJob::GapStatus gapStatus,
						 DynamicProgrammingAligner * als,
						 AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _maxOverlap(-1),
  _gapped(gapStatus==AssemblyJob::GAPPED),
  _usingMaxOverlap(false),
  _strandedness(strandedness){
  _als = new DynamicProgrammingAligner(als);
  constructorHelper(seqs);
    //OK();
}
AssemblyJobIterativeBwt::AssemblyJobIterativeBwt(set<ScoredSeq*>* seqs,
						 AssemblyJob::AssemblyStrandedness strandedness,
						 AssemblyJob::GapStatus gapStatus,
						 DynamicProgrammingAligner * als,
						 long maxOverlap,
						 AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _maxOverlap(maxOverlap),
  _gapped(gapStatus==AssemblyJob::GAPPED),
  _usingMaxOverlap(true),
  _strandedness(strandedness){
  _als = new DynamicProgrammingAligner(als);
  constructorHelper(seqs);
  //OK();
}





// separates out the seqs that are too short
void AssemblyJobIterativeBwt::constructorHelper(set<ScoredSeq*>* inputSeqs){
  _minOverlap = _als->getMinOverlap();
  set<ScoredSeq*>::iterator tssIt = _tooShortSeqs.begin();
  set<ScoredSeq*>::iterator isIt = _inputSeqs.begin();
  for (set<ScoredSeq*>::iterator it = inputSeqs->begin(); it != inputSeqs->end(); ++it) {
    if ( (*it)->size() < _minOverlap ){ tssIt = _tooShortSeqs.insert(tssIt, *it); }
    else { isIt = _inputSeqs.insert(isIt, *it); }
  }
}

AssemblyJobIterativeBwt::~AssemblyJobIterativeBwt(){
  //OK();
  delete _als;
  //delete _seqCollectionBwt;
}



void AssemblyJobIterativeBwt::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJItBwt::makeShallow can't be called after the job was run.");
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


bool AssemblyJobIterativeBwt::wasRun(){
  //OK();
  return _wasRun;
}

void AssemblyJobIterativeBwt::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  //OK();
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
  //outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
  outsideSet->insert( _tooShortSeqs.begin(), _tooShortSeqs.end() );
}

void AssemblyJobIterativeBwt::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  //OK();
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::remainingInputSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::remainingInputSeqs; seqs were already removed.");
  }
  //outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
  outsideSet->insert( _tooShortSeqs.begin(), _tooShortSeqs.end() );
}

void AssemblyJobIterativeBwt::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}

void AssemblyJobIterativeBwt::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  //OK();
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::newlyAssembledSeqs; job hasn't run yet.");
  }
  if ( _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::newlyAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}

void AssemblyJobIterativeBwt::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  //OK();
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::allAssembledSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved or _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobIterativeBwt::allAssembledSeqs; seqs were already removed.");
  }
  //outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
  outsideSet->insert( _retainedSeqs.begin(), _retainedSeqs.end() );
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
  outsideSet->insert( _tooShortSeqs.begin(), _tooShortSeqs.end() );
}



void AssemblyJobIterativeBwt::OK(){
  throw AssemblyException::ImplementationError("not updated yet.");
  if (int(_wasRun)!=0 and int(_wasRun)!=1){
    cout << int(_wasRun) << endl;
    throw AssemblyException::LogicError("_wasRun is a nonsense value.");
  }
  if (int(_inputRemoved)!=0 and int(_inputRemoved)!=1){
    cout << int(_inputRemoved) << endl;
    throw AssemblyException::LogicError("_inputRemoved is a nonsense value.");
  }
  if (int(_novelRemoved)!=0 and int(_novelRemoved)!=1){
    cout << int(_novelRemoved) << endl;
    throw AssemblyException::LogicError("_novelRemoved is a nonsense value.");
  }

  if (! _wasRun ){
    if ( _novelSeqs.size() != 0 ){
      throw AssemblyException::LogicError("AssemblyJobIterativeBwt::OK; there should be no product seqs if job hasn't run.");
    }
    if ( _novelRemoved or _inputRemoved ){
      throw AssemblyException::LogicError("AssemblyJobIterativeBwt::OK; seqs can't have been removed if job hasn't run.");
    }
  }
}



void AssemblyJobIterativeBwt::runJobGetMatchesHelper(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
						     ScoredSeqCollectionBwt* seqCollectionBwt,
						     MatchSeqTest** seqTestArray, AssemblyJob::AssemblyThreadedness threadedness){
  ScoredSeqCollectionBwt::AlignmentThreadedness useThreads;
  if (threadedness == AssemblyJob::THREADED){ useThreads = ScoredSeqCollectionBwt::threaded; }
  else { useThreads = ScoredSeqCollectionBwt::notThreaded; }

  long* minOvlArray = new long[ numQueries+1 ];
  ScoredSeqCollectionBwt::MinOvlStringency ovlStringency;

  if (_assemblyMode == AssemblyJob::FULLASSEMBLY){
    ovlStringency = ScoredSeqCollectionBwt::hardMinOvl;
    long minOvl = seqCollectionBwt->getMinOverlap();
    for (long n = 0; n < numQueries; ++n){ minOvlArray[n] = minOvl; }
  } else {
    ovlStringency = ScoredSeqCollectionBwt::softMinOvl;
    for (long n = 0; n < numQueries; ++n){ minOvlArray[n] = seqArray[n]->size(); }
  }

  if (_strandedness == AssemblyJob::DOUBLESTRANDED){
    seqCollectionBwt->getBestMatches(numQueries,matchesArray,seqArray,seqTestArray,minOvlArray,useThreads,ovlStringency);
  } else {
    seqCollectionBwt->getBestMatches(numQueries,matchesArray,seqArray,'+',seqTestArray,minOvlArray,useThreads,ovlStringency);
  }

  delete [] minOvlArray;
}



void AssemblyJobIterativeBwt::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){

    // make the queue and buffer the sequences
    set<ScoredSeq*> tempSet;
    ScoredSeqCollectionBwt * seqCollectionBwt;
    tempSet.insert( _inputSeqs.begin(), _inputSeqs.end() );
    for (set<ScoredSeq*>::iterator it = tempSet.begin(); it != tempSet.end(); ++it){ (*it)->buffer(); }

    // make the collection - it doesn't already exist
    if (_gapped){
      if (! _usingMaxOverlap){ seqCollectionBwt = new ScoredSeqCollectionBwt(&_inputSeqs,_als); }
      else { seqCollectionBwt = new ScoredSeqCollectionBwt(&_inputSeqs,_als,_maxOverlap); }
    } else {
      if (! _usingMaxOverlap){ seqCollectionBwt = new ScoredSeqCollectionBwt(&_inputSeqs,_als->getFractId(),_minOverlap); }
      else { seqCollectionBwt = new ScoredSeqCollectionBwt(&_inputSeqs,_als->getFractId(),_minOverlap,_maxOverlap); }
    }

    AlignmentScoreMatrix* alScoreMatrix = _als->getScoreMatrix();

    // I only need one "thread" if threading is turned off, which will
    // keep the memory load down
    int numThreads;
    if (threadedness == AssemblyJob::NOTTHREADED){ numThreads = 1; }
    else { numThreads = omp_get_max_threads() * 10; }

    // iterate through the queries
    while ( tempSet.size() > 0 ){

      // get a number of query seqs to match the number of threads
      int numUsedThreads;
      if (numThreads > tempSet.size()){ numUsedThreads = tempSet.size(); }
      else { numUsedThreads = numThreads; }
      ScoredSeq** querySeqs = new ScoredSeq*[ numUsedThreads ];
      for (int n = 0; n < numUsedThreads; ++n){
	set<ScoredSeq*>::iterator tempIt = tempSet.begin();
	querySeqs[n] = (*tempIt);
	tempSet.erase(tempIt);
      }

      // don't use bestMatches to prevent a good but invalid match
      // from blocking a valid one.
      vector<Alignment*>** testMatches = new vector<Alignment*>*[ numUsedThreads ];
      for (int n = 0; n < numUsedThreads; ++n){ testMatches[n] = new vector<Alignment*>; }

      // THIS IS A THREAD-SAFE BLOCK (GIVEN THAT getMatches IS THREAD-SAFE)
      // it appears twice so that a nested thread will not be created if this AJ is already threaded
      // Sequences already used will be excluded from the results.
      MatchSeqTest* matchSeqTest = new MatchSeqTestNotInSet(&_usedSeqs);
      MatchSeqTest** seqTestArray = new MatchSeqTest*[ numUsedThreads+1 ];
      for (int n = 0; n < numUsedThreads; ++n){ seqTestArray[n] = matchSeqTest; }
      runJobGetMatchesHelper(numUsedThreads, testMatches, querySeqs, seqCollectionBwt, seqTestArray, threadedness);
      delete [] seqTestArray;
      delete matchSeqTest;

      // THIS IS A NON-THREADED BLOCK
      // it is set up to prepare for a future threaded block
      vector<Alignment*>** stillValidMatches = new vector<Alignment*>*[ numUsedThreads ];
      for (int n = 0; n < numUsedThreads; ++n){ stillValidMatches[n] = new vector<Alignment*>; }
      for (int n = 0; n < numUsedThreads; ++n){
        // the query seq must not be null and must not have been consumed already
	if ( _usedSeqs.find(querySeqs[n]) != _usedSeqs.end() ){
  	  // delete the alignments, they must be discarded since one of their sequences was consumed
	  for (vector<Alignment*>::iterator hitIt = testMatches[n]->begin(); hitIt != testMatches[n]->end(); ++hitIt){ delete *hitIt; }
	} else {
	  // make sure targets weren't already removed, and pick only one of the remaining valid alignments
	  bool isFirst = true;
	  for (vector<Alignment*>::iterator hitIt = testMatches[n]->begin(); hitIt != testMatches[n]->end(); ++hitIt){
	    // check that the input seq was not already consumed
	    if ( isFirst ){
	      ScoredSeq* localSeqB = (*hitIt)->seqB();
	      if ( _usedSeqs.find(localSeqB ) == _usedSeqs.end() ){
		stillValidMatches[n]->push_back( (*hitIt) );
		// this guarantees that the target will also be consumed
		tempSet.erase( localSeqB );
		_usedSeqs.insert( localSeqB );
		isFirst = false;
	      } else { delete *hitIt; }
	    } else { delete *hitIt; }
	  }
	  // now move the sequence according to its status
	  if (stillValidMatches[n]->size() > 0){
	    _usedSeqs.insert( querySeqs[n] );
	  } else if (testMatches[n]->size() > 0){
	    // in this case, the best matches were cancelled because those sequences were used up in parallel
	    tempSet.insert( querySeqs[n] );
	  }
	}
      }

      // this is once again not thread safe
      vector<Alignment*> finalAlignments;
      vector<ScoredSeq*> normalizedToDelete;
      for (int n = 0; n < numUsedThreads; ++n){
        // create a normalized query and use it to replace seqA in the alignments
        // while moving the query/targets from _input to _used

	if (stillValidMatches[n]->size() == 1){
	  finalAlignments.push_back( *(stillValidMatches[n]->begin()) );
	} else if (stillValidMatches[n]->size() > 1){
	  // THIS IS STILL SET UP FOR THERE TO BE MULTIPLE ALIGNMENTS TO THE SAME QUERY, EVEN THOUGH THAT WILL NO LONGER BE THE CASE
	  ScoredSeqNormalized* normalizedQuery = new ScoredSeqNormalized(querySeqs[n]);
	  normalizedToDelete.push_back(normalizedQuery);
	  normalizedQuery->addCounts( stillValidMatches[n]->size() );
	  for (vector<Alignment*>::iterator hitIt = stillValidMatches[n]->begin(); hitIt != stillValidMatches[n]->end(); ++hitIt){
	    // re-check in case best alignments match the same seq
	    (*hitIt)->seqReplace( normalizedQuery, (*hitIt)->seqB() );
	    finalAlignments.push_back( (*hitIt) );
	  }
	}
      }

      // organize a threading-ready set of data structures
      Alignment** toBeCombined = new Alignment*[ finalAlignments.size() + 1 ];
      long numToCombine = 0;
      for (vector<Alignment*>::iterator hitIt = finalAlignments.begin(); hitIt != finalAlignments.end(); ++hitIt){
	toBeCombined[numToCombine] = *hitIt;
	numToCombine++;
      }
      ScoredSeq** wasCombined = new ScoredSeq*[ numToCombine + 1 ];

      // run the combine methods
      if (threadedness == AssemblyJob::NOTTHREADED){
	for (long n = 0; n < numToCombine; ++n){ wasCombined[n] = toBeCombined[n]->combine(); }
      } else {
        #pragma omp parallel for schedule(dynamic)
	for (long n = 0; n < numToCombine; ++n){ wasCombined[n] = toBeCombined[n]->combine(); }
      }

      // place the results and clean up
      for (long n = 0; n < numToCombine; ++n){
	_novelSeqs.insert( wasCombined[n] );
	delete toBeCombined[n];
      }
      delete [] toBeCombined;
      delete [] wasCombined;

      for (vector<ScoredSeq*>::iterator delIt = normalizedToDelete.begin(); delIt != normalizedToDelete.end(); ++delIt){ delete *delIt; }

      delete [] querySeqs;
      for (int n = 0; n < numUsedThreads; ++n){
	delete testMatches[n];
	delete stillValidMatches[n];
      }
      delete [] testMatches;
      delete [] stillValidMatches;
    }
    delete alScoreMatrix;
    delete seqCollectionBwt;

    // finally, sort out the used versus retained seqs
    set<ScoredSeq*>::iterator rIt = _retainedSeqs.begin();
    for (set<ScoredSeq*>::iterator inIt = _inputSeqs.begin(); inIt != _inputSeqs.end(); ++inIt){
      if (_usedSeqs.count(*inIt) == 0){ rIt = _retainedSeqs.insert(rIt, *inIt); }
    }

    _wasRun = true;
  }
}

#endif
