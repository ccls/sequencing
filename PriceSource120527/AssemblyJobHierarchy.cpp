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


#ifndef ASSEMBLYJOBHIERARCHY_CPP
#define ASSEMBLYJOBHIERARCHY_CPP



#include "AssemblyJobHierarchy.h"

#include "AssemblyJobNull.h"
#include "AssemblyJobPerfect.h"
#include "AssemblyJobNearPerfect.h"
#include "AssemblyJobSubset.h"
#include "AssemblyJobIterativeBwt.h"
#include "AssemblyException.h"
#include "AssemblerListenerNull.h"
#include "AssemblyJobSsDeBruijn.h"
#include <queue>
#include <limits>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <algorithm>
using namespace::std;

AssemblyJobHierarchy::AssemblyJobHierarchy(){}

// NOT THREADED
AssemblyJobHierarchy::AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
					   ParamsMinOverlap * pmo,
					   ParamsMinFractId * pmf,
					   AlignmentScoreMatrix * sm,
					   ParamsDeBruijn* pmdb,
					   AssemblyJob::AssemblyStrandedness strandedness,
					   AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _readCountFactor(1),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1),
  _maxRepeatFract(0.9),
  _strandedness(strandedness){
  _sm = new AlignmentScoreMatrix(sm);
  _pmo = new ParamsMinOverlap(pmo);
  _pmf = new ParamsMinFractId(pmf);
  _inputSeqs.insert( seqs->begin(), seqs->end() );
  _nullListener = new AssemblerListenerNull();
  _listener = _nullListener;
  _paramsDeBruijn = new ParamsDeBruijn(pmdb);
}
AssemblyJobHierarchy::AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
					   ParamsMinOverlap * pmo,
					   ParamsMinFractId * pmf,
					   AlignmentScoreMatrix * sm,
					   ParamsDeBruijn* pmdb,
					   AssemblyJob::AssemblyStrandedness strandedness,
					   AssemblerListener* listener,
					   AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _listener(listener),
  _readCountFactor(1),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1),
  _maxRepeatFract(0.9),
  _strandedness(strandedness){
  _sm = new AlignmentScoreMatrix(sm);
  _pmo = new ParamsMinOverlap(pmo);
  _pmf = new ParamsMinFractId(pmf);
  _inputSeqs.insert( seqs->begin(), seqs->end() );
  _nullListener = new AssemblerListenerNull();
  _paramsDeBruijn = new ParamsDeBruijn(pmdb);
}

AssemblyJobHierarchy::AssemblyJobHierarchy(float readCountFactor,
					   set<ScoredSeq*>* seqs,
					   ParamsMinOverlap * pmo,
					   ParamsMinFractId * pmf,
					   AlignmentScoreMatrix * sm,
					   ParamsDeBruijn* pmdb,
					   AssemblyJob::AssemblyStrandedness strandedness,
					   AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _readCountFactor(readCountFactor),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1),
  _maxRepeatFract(0.9),
  _strandedness(strandedness){
  _sm = new AlignmentScoreMatrix(sm);
  _pmo = new ParamsMinOverlap(pmo);
  _pmf = new ParamsMinFractId(pmf);
  _inputSeqs.insert( seqs->begin(), seqs->end() );
  _nullListener = new AssemblerListenerNull();
  _listener = _nullListener;
  _paramsDeBruijn = new ParamsDeBruijn(pmdb);
}
AssemblyJobHierarchy::AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
					   ParamsMinOverlap * pmo,
					   ParamsMinFractId * pmf,
					   AlignmentScoreMatrix * sm,
					   ParamsDeBruijn* pmdb,
					   AssemblyJob::AssemblyStrandedness strandedness,
					   long maxOverlap,
					   AssemblerListener* listener,
					   AssemblyJob::AssemblyMode assemblyMode) :
  _assemblyMode(assemblyMode),
  _listener(listener),
  _readCountFactor(1),
  _wasRun(false),
  _inputRemoved(false),
  _novelRemoved(false),
  _strandedness(strandedness),
  _maxRepeatFract(0.9),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap){
  _sm = new AlignmentScoreMatrix(sm);
  _pmo = new ParamsMinOverlap(pmo);
  _pmf = new ParamsMinFractId(pmf);
  _inputSeqs.insert( seqs->begin(), seqs->end() );
  _nullListener = new AssemblerListenerNull();
  _paramsDeBruijn = new ParamsDeBruijn(pmdb);
}




AssemblyJobHierarchy::~AssemblyJobHierarchy(){
  delete _sm;
  delete _pmo;
  delete _pmf;
  delete _nullListener;
  delete _paramsDeBruijn;
}



void AssemblyJobHierarchy::makeShallow(bool shallowDelete){
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


long AssemblyJobHierarchy::getN50(set<ScoredSeq*>* seqs){
  vector<long> contigLengths;
  long totalSize = 0;
  for (set<ScoredSeq*>::iterator seqIt = seqs->begin(); seqIt != seqs->end(); ++seqIt){
    contigLengths.push_back( (*seqIt)->size() );
    totalSize += (*seqIt)->size();
  }
  sort( contigLengths.begin(), contigLengths.end() );
  long halfSize = totalSize / 2;

  long n50 = 0;
  vector<long>::iterator sizeIt = contigLengths.begin();
  while ( halfSize > 0 and sizeIt != contigLengths.end() ){
    n50 = (*sizeIt); // current n50; will only be updated if half the total size is not reached
    halfSize -= (*sizeIt);
    ++sizeIt;
  }
  return n50;
}


void AssemblyJobHierarchy::executeAndCleanUpJob(AssemblyJob* aj, AssemblyJob::AssemblyThreadedness threadedness){
  aj->runJob(threadedness);
  _seqCarrier.clear();
  aj->allAssembledSeqs( &_seqCarrier );
  aj->discardedInputSeqs( &_localUsedSeqs );
  delete aj;
}



void AssemblyJobHierarchy::runJobHelpPerfect(AssemblyJob::AssemblyThreadedness threadedness){
  // do the perfect assembly
  AssemblyJob * asj1 = new AssemblyJobPerfect(&_seqCarrier,_strandedness);
  executeAndCleanUpJob(asj1,threadedness);

  // do the ungapped near-perfect assembly
  AssemblyJob * asj2 = new AssemblyJobNearPerfect(&_seqCarrier, _pmf->globalFractId(), _strandedness);
  executeAndCleanUpJob(asj2,threadedness);
}


void AssemblyJobHierarchy::runJobHelpDeBruijn(AssemblyJob::AssemblyThreadedness threadedness){
  // i only have a method for single-stranded de bruijn
  if (_strandedness == AssemblyJob::SINGLESTRANDED){
    set<ScoredSeq*> process;
    set<ScoredSeq*>::iterator prIt = process.begin();
    set<ScoredSeq*> dontProcess;
    set<ScoredSeq*>::iterator dprIt = dontProcess.begin();
    // these are already sorted for a set, so they can be inserted using an iterator
    for (set<ScoredSeq*>::iterator it = _seqCarrier.begin(); it != _seqCarrier.end(); ++it){
      if ( _paramsDeBruijn->inputToDeBruijn(*it) ){ prIt = process.insert(prIt, *it); }
      else { dprIt = dontProcess.insert(dprIt, *it); }
    }
    // i don't need to do anything if there are not processable seqs.
    if (process.size() >= _paramsDeBruijn->minNumSeqs()){

      // all of the input seqs are discarded
      _localUsedSeqs.insert(process.begin(), process.end());
      AssemblyJob* asjDb = new AssemblyJobSsDeBruijn(&process, _paramsDeBruijn);
      asjDb->runJob(threadedness);
      process.clear();
      asjDb->allAssembledSeqs( &process );
      delete asjDb;

      // collapse the output seqs and add back the unprocessed seqs
      AssemblyJob* asjSub = new AssemblyJobSubset(&process,
						  _pmf->calculateMinFractId( &process, _readCountFactor ),
						  _strandedness);
      executeAndCleanUpJob(asjSub,threadedness);

      _seqCarrier.insert(dontProcess.begin(), dontProcess.end());
    }
  }
}


void AssemblyJobHierarchy::runJobHelpIterativeBwt(AssemblyJob::GapStatus gapStatus, AssemblyJob::AssemblyThreadedness threadedness){
  // calculate the initial params
  long correctedSeqNum = long(float(_seqCarrier.size()) * _readCountFactor);
  long ogCorrectedSeqNum = correctedSeqNum;
  //float minFractId = _pmf->calculateMinFractId( correctedSeqNum );

  float minFractId = _pmf->calculateMinFractId(&_seqCarrier, _readCountFactor);
  long minOverlap = _pmo->calculateMinOverlap( correctedSeqNum );

  bool doRound = true;
  while ( doRound ){

    // create the job
    AssemblyJob* asj;
    DynamicProgrammingAligner * ast = new DynamicProgrammingAligner(minFractId,minOverlap,_sm);
    if (! _usingMaxOverlap){ asj = new AssemblyJobIterativeBwt(&_seqCarrier, _strandedness, gapStatus, ast, _assemblyMode); }
    else { asj = new AssemblyJobIterativeBwt(&_seqCarrier, _strandedness, gapStatus, ast, _maxOverlap, _assemblyMode); }
    delete ast;

    // run the job
    executeAndCleanUpJob(asj,threadedness);

    // clean up newly redundant seqs if anything has happened
    long localCorrectedSeqNum = long(float(_seqCarrier.size()) * _readCountFactor);

    // determine if another round is called for and, if so, update the param values
    if ( localCorrectedSeqNum < correctedSeqNum ){
      //minFractId = _pmf->calculateMinFractId( localCorrectedSeqNum );
      minFractId = _pmf->calculateMinFractId( &_seqCarrier, _readCountFactor );
      minOverlap = _pmo->calculateMinOverlap( localCorrectedSeqNum );
    } else { doRound = false; }
    correctedSeqNum = localCorrectedSeqNum;
  }

  // do the subset job only at the end and only if progress was made above
  if ( correctedSeqNum < ogCorrectedSeqNum ){
    AssemblyJob* asjSub1 = new AssemblyJobSubset(&_seqCarrier,
						 _pmf->calculateMinFractId( &_seqCarrier, _readCountFactor ),
						 _strandedness);
    executeAndCleanUpJob(asjSub1,threadedness);
  }
}


void AssemblyJobHierarchy::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if ( _inputSeqs.size()==0 ){ _wasRun = true; }
  else if (! _wasRun){
    _seqCarrier.insert(_inputSeqs.begin(), _inputSeqs.end());

    // first, do the perfect jobs
    runJobHelpPerfect(threadedness);
    runJobHelpDeBruijn(threadedness);

    // second, give the sequences chances to combine with their subseqs (if there is even a point)
    if (_seqCarrier.size() > 1){
      long correctedSeqNum = long(float(_seqCarrier.size()) * _readCountFactor);
      //float minFractId = _pmf->calculateMinFractId( correctedSeqNum );
      float minFractId = _pmf->calculateMinFractId( &_seqCarrier, _readCountFactor );
      long minOverlap = _pmo->calculateMinOverlap( correctedSeqNum );
      AssemblyJob* asjSub1 = new AssemblyJobSubset(&_seqCarrier, minFractId, _strandedness);
      executeAndCleanUpJob(asjSub1,threadedness);
    }

    // do the iterative bwt jobs (if there is even a point)
    //if (_seqCarrier.size() > 1){ runJobHelpIterativeBwt( AssemblyJob::ungapped, threadedness ); }
    if (_seqCarrier.size() > 1){ runJobHelpIterativeBwt( AssemblyJob::GAPPED, threadedness ); }

    // clean up
    set<ScoredSeq*> tempUsedSeqs;
    set<ScoredSeq*>::iterator tusIt = tempUsedSeqs.begin();
    for (set<ScoredSeq*>::iterator it = _localUsedSeqs.begin(); it != _localUsedSeqs.end(); ++it){
      if ( _inputSeqs.find( (*it) ) == _inputSeqs.end() ){ (*it)->deepDelete(); }
      else { tusIt = tempUsedSeqs.insert(tusIt, *it); }
    }
    _usedSeqs.insert(tempUsedSeqs.begin(), tempUsedSeqs.end());
    _localUsedSeqs.clear();

    // sort the output and delete those seqs that were not in the official input nor the assembled output
    _allAssembledSeqs.insert( _seqCarrier.begin(), _seqCarrier.end() );

    set<ScoredSeq*> tempNovelSeqs;
    set<ScoredSeq*>::iterator tnsIt = tempNovelSeqs.begin();
    set<ScoredSeq*> tempRemainSeqs;
    set<ScoredSeq*>::iterator trsIt = tempRemainSeqs.begin();

    for (set<ScoredSeq*>::iterator it = _allAssembledSeqs.begin(); it != _allAssembledSeqs.end(); ++it){
      if ( _inputSeqs.find( (*it) ) == _inputSeqs.end() ){
	tnsIt = tempNovelSeqs.insert(tnsIt, *it);
	//_novelSeqs.insert( (*it) );
      } else {
	trsIt = tempRemainSeqs.insert(trsIt, *it);
	//_remainingSeqs.insert( (*it) );
      }
    }
    _novelSeqs.insert(tempNovelSeqs.begin(), tempNovelSeqs.end());
    _remainingSeqs.insert(tempRemainSeqs.begin(), tempRemainSeqs.end());

    _wasRun = true;
  }

  //OK();
}





bool AssemblyJobHierarchy::wasRun(){ return _wasRun; }

void AssemblyJobHierarchy::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}

void AssemblyJobHierarchy::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::remainingInputSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::remainingInputSeqs; seqs were already removed.");
  }
  outsideSet->insert( _remainingSeqs.begin(), _remainingSeqs.end() );
}

void AssemblyJobHierarchy::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}

void AssemblyJobHierarchy::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::newlyAssembledSeqs; job hasn't run yet.");
  }
  if ( _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::newlyAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}

void AssemblyJobHierarchy::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  if (! _wasRun){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::allAssembledSeqs; job hasn't run yet.");
  }
  if ( _inputRemoved or _novelRemoved ){
    throw AssemblyException::LogicError("in AssemblyJobHierarchy::allAssembledSeqs; seqs were already removed.");
  }
  outsideSet->insert( _allAssembledSeqs.begin(), _allAssembledSeqs.end() );
}



void AssemblyJobHierarchy::OK(){
  if (! _wasRun ){
    if ( _novelSeqs.size() != 0 or _allAssembledSeqs.size() != 0 or _remainingSeqs.size() != 0 ){
      throw AssemblyException::LogicError("AssemblyJobHierarchy::OK; there should be no product seqs if job hasn't run.");
    }
    if ( _novelRemoved or _inputRemoved ){
      throw AssemblyException::LogicError("AssemblyJobHierarchy::OK; seqs can't have been removed if job hasn't run.");
    }
  }
}



#endif
