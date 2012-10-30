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


#ifndef ASSEMBLYJOBGRAPH_CPP
#define ASSEMBLYJOBGRAPH_CPP

#include "AssemblyJobGraph.h"
#include "AssemblyJobNullCopy.h"

#include "AssemblyException.h"
#include <queue>
#include <iostream>
using namespace::std;


AssemblyJob* AssemblyJobGraph::AssemblyJobGraphFactory(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode,
						       AssemblyJobFactory* jobFactory){
  AssemblyJob* aj = new AssemblyJobGraph(beginNode, endNode, jobFactory);
  set<ScoredSeq*> input;
  aj->allInputSeqs( &input );
  if (input.size() < 2){
    delete aj;
    aj = new AssemblyJobNullCopy( &input );
  }
  return aj;
}


AssemblyJobGraph::AssemblyJobGraph(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode,
				   AssemblyJobFactory* jobFactory){
  constructorHelper(beginNode,endNode,jobFactory);
}


void AssemblyJobGraph::constructorHelper(AssemblyGraphNode* beginNode, AssemblyGraphNode* endNode, AssemblyJobFactory* ajf){
  _inputRemoved = false;
  _wasRun = false;
  _jobFactory = new AssemblyJobFactory(ajf);
  _isShallow = false;

  // check this requirement - if it is not met, bad things will happen
  if (_jobFactory->strandedness() != AssemblyJob::SINGLESTRANDED){
    throw AssemblyException::ArgError("AJGraph, factory cannot be ds");
  }

  // figure out how many nodes there are
  _numNodes = 1;
  AssemblyGraphNode* markNode = beginNode;
  while (markNode != endNode){
    markNode = markNode->getPlusLink();
    _numNodes++;
  }

  // make an input set for each edge (including tips) - include just the nested seqs, delete flips
  _edgeAssembledPlus = new set<ScoredSeq*>[_numNodes + 1];
  _edgeAssembledMinus = new set<ScoredSeq*>[_numNodes + 1];
  _minOverlaps = new long[_numNodes];
  _needsToBeAssembled = new bool[_numNodes + 1];
  for (long n = 0; n <= _numNodes; ++n){ _needsToBeAssembled[n] = true; }

  // FIRST, organize the input seqs (delete the flips)
  _rawInputSetsPlus = new set<ScoredSeq*>[_numNodes + 1];
  _rawInputSetsMinus = new set<ScoredSeq*>[_numNodes + 1];

  // and set these up
  _normInputSetsPlus = new set<ScoredSeq*>[_numNodes + 1];
  _normInputSetsMinus = new set<ScoredSeq*>[_numNodes + 1];

  // take care of the leading edge
  set<ScoredSeqFlip*> localCollector;
  beginNode->getMinusEdgeSeqs(&localCollector);
  for (set<ScoredSeqFlip*>::iterator it = localCollector.begin(); it != localCollector.end(); ++it){
    switch( (*it)->getSense() ){
    case '+': _rawInputSetsPlus[0].insert( (*it)->getNested() ); break;
    case '-': _rawInputSetsMinus[0].insert( (*it)->getNested() ); break;
    default: throw AssemblyException::ArgError("AJG, bad sense char");
    }
    delete *it;
  }
  _rawInputTotalPlus.insert( _rawInputSetsPlus[0].begin(), _rawInputSetsPlus[0].end() );
  _rawInputTotalMinus.insert( _rawInputSetsMinus[0].begin(), _rawInputSetsMinus[0].end() );

  // take care of the linking edges
  AssemblyGraphNode* firstNode = beginNode;
  _minOverlaps[0] = firstNode->getCenterContig()->size();
  long setIndex = 1;
  while (firstNode != endNode){
    AssemblyGraphNode* secondNode = firstNode->getPlusLink();
    _minOverlaps[setIndex] = secondNode->getCenterContig()->size();
    set<ScoredSeqFlip*> localCollector;
    firstNode->getPlusEdgeSeqs(&localCollector);
    secondNode->getMinusEdgeSeqs(&localCollector);
    for (set<ScoredSeqFlip*>::iterator it = localCollector.begin(); it != localCollector.end(); ++it){
      switch( (*it)->getSense() ){
      case '+': _rawInputSetsPlus[setIndex].insert( (*it)->getNested() ); break;
      case '-': _rawInputSetsMinus[setIndex].insert( (*it)->getNested() ); break;
      default: throw AssemblyException::ArgError("AJG, bad sense char");
      }
      delete *it;
    }
    _rawInputTotalPlus.insert( _rawInputSetsPlus[setIndex].begin(), _rawInputSetsPlus[setIndex].end() );
    _rawInputTotalMinus.insert( _rawInputSetsMinus[setIndex].begin(), _rawInputSetsMinus[setIndex].end() );
    setIndex++;
    firstNode = secondNode;
  }

  // take care of the trailing edge
  localCollector.clear();
  endNode->getPlusEdgeSeqs(&localCollector);
  for (set<ScoredSeqFlip*>::iterator it = localCollector.begin(); it != localCollector.end(); ++it){
    switch( (*it)->getSense() ){
    case '+': _rawInputSetsPlus[_numNodes].insert( (*it)->getNested() ); break;
    case '-': _rawInputSetsMinus[_numNodes].insert( (*it)->getNested() ); break;
    default: throw AssemblyException::ArgError("AJG, bad sense char");
    }
    delete *it;
  }
  _rawInputTotalPlus.insert( _rawInputSetsPlus[_numNodes].begin(), _rawInputSetsPlus[_numNodes].end() );
  _rawInputTotalMinus.insert( _rawInputSetsMinus[_numNodes].begin(), _rawInputSetsMinus[_numNodes].end() );
}


AssemblyJobGraph::~AssemblyJobGraph(){
  delete _jobFactory;
  delete [] _edgeAssembledPlus;
  delete [] _edgeAssembledMinus;
  delete [] _minOverlaps;
  delete [] _needsToBeAssembled;
  delete [] _rawInputSetsPlus;
  delete [] _rawInputSetsMinus;
  delete [] _normInputSetsPlus;
  delete [] _normInputSetsMinus;
}




void AssemblyJobGraph::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJGraph::makeShallow can't be called after the job was run.");
  }

  // THE OLD WAY STARTS HERE
  set<ScoredSeq*> totalInput;
  totalInput.insert(_rawInputTotalPlus.begin(), _rawInputTotalPlus.end());
  totalInput.insert(_rawInputTotalMinus.begin(), _rawInputTotalMinus.end());
  for (set<ScoredSeq*>::iterator it = totalInput.begin(); it != totalInput.end(); ++it){ (*it)->bottomBuffer(); }

  long numNodesP1 = _numNodes + 1;
  set<ScoredSeq*> addToTotalPlus;
  set<ScoredSeq*>* addToPlusSets = new set<ScoredSeq*>[numNodesP1];
  set<ScoredSeq*> addToTotalMinus;
  set<ScoredSeq*>* addToMinusSets = new set<ScoredSeq*>[numNodesP1];

  for (set<ScoredSeq*>::iterator it = totalInput.begin(); it != totalInput.end(); ++it){
    ScoredSeq* copy = (*it)->shallowCopy();
    (*it)->deepUnbuffer();

    // replace if found in the plus sets
    set<ScoredSeq*>::iterator foundPlusTotalIt = _rawInputTotalPlus.find( *it );
    if (foundPlusTotalIt != _rawInputTotalPlus.end()){
      _rawInputTotalPlus.erase( foundPlusTotalIt );
      addToTotalPlus.insert( copy );
      // now check all of the subsets
      for (long n = 0; n < numNodesP1; ++n){
	set<ScoredSeq*>::iterator foundItN = _rawInputSetsPlus[n].find( *it );
	if (foundItN != _rawInputSetsPlus[n].end()){
	  _rawInputSetsPlus[n].erase( foundItN );
	  addToPlusSets[n].insert( copy );
	}
      }
    }

    // replace if found in the minus sets
    set<ScoredSeq*>::iterator foundMinusTotalIt = _rawInputTotalMinus.find( *it );
    if (foundMinusTotalIt != _rawInputTotalMinus.end()){
      _rawInputTotalMinus.erase( foundMinusTotalIt );
      addToTotalMinus.insert( copy );
      // now check all of the subsets
      for (long n = 0; n < numNodesP1; ++n){
	set<ScoredSeq*>::iterator foundItN = _rawInputSetsMinus[n].find( *it );
	if (foundItN != _rawInputSetsMinus[n].end()){
	  _rawInputSetsMinus[n].erase( foundItN );
	  addToMinusSets[n].insert( copy );
	}
      }
    }
  }

  // replace the contents of the now-empty total input
  if (_rawInputTotalPlus.size() != 0){ throw AssemblyException::LogicError("AJGraph::makeShallow, plus set should be empty."); }
  if (_rawInputTotalMinus.size() != 0){ throw AssemblyException::LogicError("AJGraph::makeShallow, minus set should be empty."); }
  _rawInputTotalPlus.insert(addToTotalPlus.begin(), addToTotalPlus.end());
  _rawInputTotalMinus.insert(addToTotalMinus.begin(), addToTotalMinus.end());

  // replace the contents of the now-empty individual subsets
  for (long n = 0; n < numNodesP1; ++n){
    if (_rawInputSetsPlus[n].size() != 0){ throw AssemblyException::LogicError("AJGraph::makeShallow, plus subset should be empty."); }
    if (_rawInputSetsMinus[n].size() != 0){ throw AssemblyException::LogicError("AJGraph::makeShallow, minus subset should be empty."); }
    _rawInputSetsPlus[n].insert(addToPlusSets[n].begin(), addToPlusSets[n].end());
    _rawInputSetsMinus[n].insert(addToMinusSets[n].begin(), addToMinusSets[n].end());
  }
  delete [] addToPlusSets;
  delete [] addToMinusSets;

  _isShallow = true;
}


void AssemblyJobGraph::createNormSetsHelper(long seqN, char sense, ScoredSeq** seqArray, long* numInstances, bool* inArray, long* indexArray){
  bool onThisStrand;
  set<ScoredSeq*>* rawInputSets;
  switch (sense){
  case '+':
    onThisStrand = _rawInputTotalPlus.count( seqArray[seqN] ) > 0;
    rawInputSets = _rawInputSetsPlus;
    break;
  case '-':
    onThisStrand = _rawInputTotalMinus.count( seqArray[seqN] ) > 0;
    rawInputSets = _rawInputSetsMinus;
    break;
  default:
    throw AssemblyException::ArgError("bad sense char in AJG::createNormSetsHelper");
  }
  if (onThisStrand){
    inArray[seqN] = true;
    for (long n = 0; n < _numNodes+1; ++n){
      if (rawInputSets[n].count( seqArray[seqN] ) > 0){
	numInstances[seqN]++;
	indexArray[seqN] = n;
      }
    }
  } else { inArray[seqN] = false; }
}

void AssemblyJobGraph::createNormSets(AssemblyJob::AssemblyThreadedness threadedness){

  // create a non-redundant array so that the process below can be threaded
  set<ScoredSeq*> totalInput;
  totalInput.insert(_rawInputTotalPlus.begin(), _rawInputTotalPlus.end());
  totalInput.insert(_rawInputTotalMinus.begin(), _rawInputTotalMinus.end());
  ScoredSeq** seqArray = new ScoredSeq*[totalInput.size() + 1];
  long inputCount = 0;
  for (set<ScoredSeq*>::iterator it = totalInput.begin(); it != totalInput.end(); ++it){
    seqArray[inputCount] = *it;
    inputCount++;
  }
  long* numInstances = new long[inputCount + 1];
  bool* inPlus = new bool[inputCount + 1];
  long* indexPlus = new long[inputCount + 1];
  bool* inMinus = new bool[inputCount + 1];
  long* indexMinus = new long[inputCount + 1];

  // this can be threaded
  if (threadedness == AssemblyJob::THREADED){
    #pragma omp parallel for schedule(static)
    for (long seqN = 0; seqN < inputCount; ++seqN){
      numInstances[seqN] = 0;
      createNormSetsHelper(seqN, '+', seqArray, numInstances, inPlus, indexPlus);
      createNormSetsHelper(seqN, '-', seqArray, numInstances, inMinus, indexMinus);
    }
  } else {
    for (long seqN = 0; seqN < inputCount; ++seqN){
      numInstances[seqN] = 0;
      createNormSetsHelper(seqN, '+', seqArray, numInstances, inPlus, indexPlus);
      createNormSetsHelper(seqN, '-', seqArray, numInstances, inMinus, indexMinus);
    }
  }

  // this should not be threaded because it modifies sets
  for (long seqN = 0; seqN < inputCount; ++seqN){
    if (numInstances[seqN]==1){
      if (inPlus[seqN]){ _normInputSetsPlus[indexPlus[seqN]].insert( seqArray[seqN] ); }
      else { _normInputSetsMinus[indexMinus[seqN]].insert( seqArray[seqN] ); }
    } else if (numInstances[seqN]==2 and inPlus[seqN] and inMinus[seqN]){
      ScoredSeqNormalized* normShell = new ScoredSeqNormalized( seqArray[seqN] );
      normShell->addCounts( 2 );
      _normToLivingCount.insert(pair<ScoredSeq*,long>(normShell,2));
      _justNormPlus.insert( normShell );
      _justNormMinus.insert( normShell );
      _normInputSetsPlus[indexPlus[seqN]].insert( normShell );
      _normInputSetsMinus[indexMinus[seqN]].insert( normShell );
    } else if (numInstances[seqN] < 1){
      throw AssemblyException::LogicError("AJGraph::createNormSets seq cannot be used zero times");
    } else {
      long localCount = 0;
      ScoredSeqNormalized* normShell = new ScoredSeqNormalized( seqArray[seqN] );
      normShell->addCounts( numInstances[seqN] );
      _normToLivingCount.insert(pair<ScoredSeq*,long>(normShell,numInstances[seqN]));
      if (inPlus[seqN]){
	_justNormPlus.insert( normShell );
	for (long n = 0; n < _numNodes+1; ++n){
	  if (_rawInputSetsPlus[n].count( seqArray[seqN] ) > 0){ _normInputSetsPlus[n].insert( normShell ); localCount++; }
	}
      }
      if (inMinus[seqN]){
	_justNormMinus.insert( normShell );
	for (long n = 0; n < _numNodes+1; ++n){
	  if (_rawInputSetsMinus[n].count( seqArray[seqN] ) > 0){ _normInputSetsMinus[n].insert( normShell ); localCount++; }
	}
      }
      if (localCount !=numInstances[seqN]){
        throw AssemblyException::LogicError("AJGraph::createNormSets seq is included an inappropriate num of times.");
      }
    }
  }
  delete [] seqArray;
  delete [] numInstances;
  delete [] inPlus;
  delete [] indexPlus;
  delete [] inMinus;
  delete [] indexMinus;
}

void AssemblyJobGraph::insertFullDiscard(ScoredSeq* seq, bool delSeq){
  map<ScoredSeq*,long>::iterator foundIt = _normToLivingCount.find(seq);
  if (foundIt == _normToLivingCount.end()){
    if (delSeq){ seq->deepDelete(); }
    else { _fullDiscardSet.insert( seq ); }
  } else { foundIt->second--; }
}



void AssemblyJobGraph::runCenterAssemblyJob(long tipIndex, long edgeIndex, AssemblyJob::AssemblyThreadedness threadedness){

  // collect all of the input sequences and flip them appropriately
  set<ScoredSeq*> jobSet;

  // these set up for threading of the SSFlip generation and set addition for the four sets;
  // the number of sets is not variable, so these are all declared in scope (no deletion necessary)
  set<ScoredSeq*>* sourceJobSets[] = {&_edgeAssembledPlus[tipIndex], &_edgeAssembledMinus[tipIndex],
				      &_edgeAssembledPlus[edgeIndex], &_edgeAssembledMinus[edgeIndex]};
  char flipDirections[] = {'+','-','+','-'};
  set<ScoredSeq*> destinationJobSets[4];

  if (threadedness == THREADED){
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < 4; ++n){
      for (set<ScoredSeq*>::iterator it = sourceJobSets[n]->begin(); it != sourceJobSets[n]->end(); ++it){
	destinationJobSets[n].insert( ScoredSeqFlip::getFlip(*it,flipDirections[n]) );
      }
      sourceJobSets[n]->clear();
    }
  } else {
    for (int n = 0; n < 4; ++n){
      for (set<ScoredSeq*>::iterator it = sourceJobSets[n]->begin(); it != sourceJobSets[n]->end(); ++it){
	destinationJobSets[n].insert( ScoredSeqFlip::getFlip(*it,flipDirections[n]) );
      }
      sourceJobSets[n]->clear();
    }
  }

  // combining the four sets will now have a linear order-of-growth and is not threaded
  for (int n = 0; n < 4; ++n){ jobSet.insert(destinationJobSets[n].begin(), destinationJobSets[n].end()); }

  // switch out the min overlap requirement and combine the two edge results, to generate
  // a new tip edge set of seqs.
  long minOverlap;
  if (tipIndex > edgeIndex){ minOverlap = _minOverlaps[edgeIndex]; }
  else { minOverlap = _minOverlaps[tipIndex]; }

  // raise the bar on the min overlap requirement and nullify any max overlap
  ParamsMinOverlap* oldPmo = _jobFactory->getParamsMinOverlap();
  ParamsMinOverlap* newPmo = new ParamsMinOverlap(oldPmo,  minOverlap);
  AlignmentScoreMatrix* asMatrix = _jobFactory->getScoreMatrix();
  ParamsMinFractId* pmf = _jobFactory->getParamsMinFractId();
  ParamsDeBruijn* pdb = _jobFactory->getParamsDeBruijn();
  AssemblyJobFactory* metaFactory = new AssemblyJobFactory(newPmo, pmf, asMatrix, pdb, AssemblyJob::SINGLESTRANDED, _jobFactory->getListener());
  delete oldPmo;
  delete newPmo;
  delete asMatrix;
  delete pmf;
  delete pdb;

  // make and run the job and collect output into the same set that was passed in as an arg
  AssemblyJob* aj = metaFactory->assemblyJob( &jobSet );
  delete metaFactory;
  aj->runJob(threadedness);
  set<ScoredSeq*> leftoverSeqs;
  aj->discardedInputSeqs( &leftoverSeqs );
  // collect these separately to avoid lookup below
  set<ScoredSeq*> novelSeqs;
  aj->newlyAssembledSeqs( &novelSeqs );
  set<ScoredSeq*> preservedSeqs;
  aj->remainingInputSeqs( &preservedSeqs );
  delete aj;

  for (set<ScoredSeq*>::iterator it = leftoverSeqs.begin(); it != leftoverSeqs.end(); ++it){
    ScoredSeq* innerSeq = dynamic_cast<ScoredSeqFlip*>(*it)->getNested();
    // i need to check the entire input set because a sequence may have propagated over from
    // another edge job if it had a similar sequence to a sequence from this edge job
    if (_normToLivingCount.find(innerSeq) == _normToLivingCount.end()){ innerSeq->deepDelete(); }
    else { insertFullDiscard( innerSeq, false ); }
  }

  // place the assembled sequences appropriately
  for (set<ScoredSeq*>::iterator it = preservedSeqs.begin(); it != preservedSeqs.end(); ++it){
    ScoredSeqFlip* recast = dynamic_cast<ScoredSeqFlip*>( *it );
    switch ( recast->getSense() ){
    case '+': _edgeAssembledPlus[edgeIndex].insert( recast->getNested() ); break;
    case '-': _edgeAssembledMinus[edgeIndex].insert( recast->getNested() ); break;
    default: throw AssemblyException::LogicError("AJG::assembleOneEdge, bad sense char");
    }
  }
  // this operation should come second since it is linear, but would increase the log time
  // for the addition of individual sequence above
  _edgeAssembledPlus[edgeIndex].insert(novelSeqs.begin(), novelSeqs.end());

  // delete the wrappers
  for (set<ScoredSeq*>::iterator it = jobSet.begin(); it != jobSet.end(); ++it){ delete *it; }
}


void AssemblyJobGraph::assembleOneEdge(long edgeIndex, AssemblyJob::AssemblyThreadedness threadedness){
  // create strand-appropriate wrappers
  set<ScoredSeq*> jobSet;

  // these set up for threading of the SSFlip generation and set addition for the four sets;
  // the number of sets is not variable, so these are all declared in scope (no deletion necessary)
  set<ScoredSeq*>* sourceJobSets[] = {&_normInputSetsPlus[edgeIndex], &_normInputSetsMinus[edgeIndex]};
  char flipDirections[] = {'+','-'};
  set<ScoredSeq*> destinationJobSets[2];

  if (threadedness == THREADED){
    #pragma omp parallel for schedule(static)
    for (int n = 0; n < 2; ++n){
      for (set<ScoredSeq*>::iterator it = sourceJobSets[n]->begin(); it != sourceJobSets[n]->end(); ++it){
	destinationJobSets[n].insert( ScoredSeqFlip::getFlip(*it,flipDirections[n]) );
      }
    }
  } else {
    for (int n = 0; n < 2; ++n){
      for (set<ScoredSeq*>::iterator it = sourceJobSets[n]->begin(); it != sourceJobSets[n]->end(); ++it){
	destinationJobSets[n].insert( ScoredSeqFlip::getFlip(*it,flipDirections[n]) );
      }
    }
  }

  // combining the four sets will now have a linear order-of-growth and is not threaded
  for (int n = 0; n < 2; ++n){ jobSet.insert(destinationJobSets[n].begin(), destinationJobSets[n].end()); }

  // run the job
  AssemblyJob* aj = _jobFactory->assemblyJob( &jobSet );
  aj->runJob(threadedness);
  set<ScoredSeq*> leftoverSeqs;
  aj->discardedInputSeqs( &leftoverSeqs );
  // collect these separately to avoid lookup below
  set<ScoredSeq*> novelSeqs;
  aj->newlyAssembledSeqs( &novelSeqs );
  set<ScoredSeq*> preservedSeqs;
  aj->remainingInputSeqs( &preservedSeqs );
  delete aj;

  for (set<ScoredSeq*>::iterator it = leftoverSeqs.begin(); it != leftoverSeqs.end(); ++it){
    insertFullDiscard( dynamic_cast<ScoredSeqNested*>(*it)->getNested(), false );
  }

  // place the assembled sequences appropriately
  for (set<ScoredSeq*>::iterator it = preservedSeqs.begin(); it != preservedSeqs.end(); ++it){
    ScoredSeqFlip* recast = dynamic_cast<ScoredSeqFlip*>( *it );
    switch ( recast->getSense() ){
    case '+': _edgeAssembledPlus[edgeIndex].insert( recast->getNested() ); break;
    case '-': _edgeAssembledMinus[edgeIndex].insert( recast->getNested() ); break;
    default: throw AssemblyException::LogicError("AJG::assembleOneEdge, bad sense char");
    }
  }
  // this operation should come second since it is linear, but would increase the log time
  // for the addition of individual sequence above
  _edgeAssembledPlus[edgeIndex].insert(novelSeqs.begin(), novelSeqs.end());

  // place the assembled sequences appropriately

  // delete the wrappers
  for (set<ScoredSeq*>::iterator it = jobSet.begin(); it != jobSet.end(); ++it){ delete *it; }
  _needsToBeAssembled[edgeIndex] = false;
}



void AssemblyJobGraph::collapseOneNode(long tipIndex, long edgeIndex, set<ScoredSeq*>* localProducts,
				       AssemblyJob::AssemblyThreadedness threadedness){

  // make sure that the two input sets have been dealt with
  if ( _needsToBeAssembled[tipIndex] ){ assembleOneEdge(tipIndex, threadedness); }
  if ( _needsToBeAssembled[edgeIndex] ){ assembleOneEdge(edgeIndex, threadedness); }

  // I want to retain just those seqs that derive in some way from the edgeIndex
  // set for the next tip set.  So I will make a collection for mapping before I
  // empty that set.  These are copies so I can delete the originals during assembly.
  set<ScoredSeq*> anchorSet;
  // making these copies could be time-consuming if there are a lot of assembly output sequences and/or
  // they are long - so it should be threaded when possible.  creating the copies can be threaded:

  long numAnchorPlus = _edgeAssembledPlus[edgeIndex].size();
  ScoredSeq** sourceAnchorPlus = new ScoredSeq*[ numAnchorPlus + 1 ];
  long apN = 0;
  for(set<ScoredSeq*>::iterator it = _edgeAssembledPlus[edgeIndex].begin();
      it != _edgeAssembledPlus[edgeIndex].end(); ++it){ sourceAnchorPlus[apN] = *it; ++apN; }
  ScoredSeq** copyAnchorPlus = new ScoredSeq*[ numAnchorPlus + 1 ];
  if (threadedness == THREADED){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numAnchorPlus; ++n){ copyAnchorPlus[n] = sourceAnchorPlus[n]->shallowCopy(); }
  } else {
    for (long n = 0; n < numAnchorPlus; ++n){ copyAnchorPlus[n] = sourceAnchorPlus[n]->shallowCopy(); }
  }
  delete [] sourceAnchorPlus;

  long numAnchorMinus = _edgeAssembledMinus[edgeIndex].size();
  ScoredSeq** sourceAnchorMinus = new ScoredSeq*[ numAnchorMinus + 1 ];
  long amN = 0;
  for(set<ScoredSeq*>::iterator it = _edgeAssembledMinus[edgeIndex].begin();
      it != _edgeAssembledMinus[edgeIndex].end(); ++it){ sourceAnchorMinus[amN] = *it; ++amN; }
  ScoredSeq** copyAnchorMinus = new ScoredSeq*[ numAnchorMinus + 1 ];
  if (threadedness == THREADED){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numAnchorMinus; ++n){ copyAnchorMinus[n] = sourceAnchorMinus[n]->shallowCopy(); }
  } else {
    for (long n = 0; n < numAnchorMinus; ++n){ copyAnchorMinus[n] = sourceAnchorMinus[n]->shallowCopy(); }
  }
  delete [] sourceAnchorMinus;

  // the two sets can also be sorted in a threaded manner
  ScoredSeq** copyAnchors[] = {copyAnchorPlus, copyAnchorMinus};
  long numAnchors[] = {numAnchorPlus, numAnchorMinus};
  set<ScoredSeq*> preAnchorSets[2];
  if (threadedness == THREADED){
    #pragma omp parallel for schedule(static)
    for (int aN = 0; aN < 2; ++aN){
      for (long n = 0; n < numAnchors[aN]; ++n){ preAnchorSets[aN].insert(copyAnchors[aN][n]); }
      delete [] copyAnchors[aN];
    }
  } else {
    for (int aN = 0; aN < 2; ++aN){
      for (long n = 0; n < numAnchors[aN]; ++n){ preAnchorSets[aN].insert(copyAnchors[aN][n]); }
      delete [] copyAnchors[aN];
    }
  }

  // now the two sets can be combined in linear time - also, create an array version for later
  for (int aN = 0; aN < 2; ++aN){ anchorSet.insert(preAnchorSets[aN].begin(), preAnchorSets[aN].end()); }
  long totalAnchors = anchorSet.size();
  ScoredSeq** fullAnchorArray = new ScoredSeq*[ totalAnchors + 1 ];
  long* fullAnchorLengths = new long[ totalAnchors + 1 ];
  long faaN = 0;
  for (set<ScoredSeq*>::iterator it = anchorSet.begin(); it != anchorSet.end(); ++it){
    fullAnchorArray[faaN] = *it; fullAnchorLengths[faaN] = (*it)->size(); ++faaN;
  }

  // collapse the anchor set - this won't take long and will save mapping time
  /* THIS IS PROBABLY NOT WORTH IT BECAUSE OF THE OVERHEAD OF BUILDING THE BWT
  AssemblyJob* anchorJob = _jobFactory->redundancyJob(&anchorSet);
  anchorJob->runJob(threadedness);
  anchorSet.clear();
  anchorJob->discardedInputSeqs(&anchorSet);
  for (set<ScoredSeq*>::iterator it = anchorSet.begin(); it != anchorSet.end(); ++it){ (*it)->deepDelete(); }
  anchorSet.clear();
  anchorJob->allAssembledSeqs(&anchorSet);
  delete anchorJob;
  */
  // run the job, then create one set of the products for the purpose of calculating parameters
  runCenterAssemblyJob(tipIndex, edgeIndex, threadedness);
  set<ScoredSeq*> justForParams;
  justForParams.insert( _edgeAssembledPlus[edgeIndex].begin(), _edgeAssembledPlus[edgeIndex].end() );
  justForParams.insert( _edgeAssembledMinus[edgeIndex].begin(), _edgeAssembledMinus[edgeIndex].end() );

  // use the properties of the assembled set to determine the alignment parameters;
  // min overlap will be specified for each search based on the query length
  ParamsMinFractId* pmf = _jobFactory->getParamsMinFractId();
  ParamsMinOverlap* pmo = _jobFactory->getParamsMinOverlap();
  // gapped alignment to maximize retention
  AlignmentScoreMatrix* asMatrix = _jobFactory->getScoreMatrix();
  DynamicProgrammingAligner* dpa = new DynamicProgrammingAligner(pmf->calculateMinFractId(&justForParams),
								 pmo->calculateMinOverlap(justForParams.size()),
								 asMatrix);
  delete asMatrix;
  ScoredSeqCollectionBwt* anchorSeqCollection = new ScoredSeqCollectionBwt(&anchorSet, dpa);
  // don't delete dpa yet, it will be used to construct the contig collections below
  delete pmf;
  delete pmo;

  char senseChoices[] = {'+','-'};
  set<ScoredSeq*>* edgeAssembledChoices[] = { _edgeAssembledPlus, _edgeAssembledMinus };
  for (int senseN = 0; senseN < 2; ++senseN){
    // establish some loop variables
    char sense = senseChoices[senseN];
    set<ScoredSeq*>* edgeAssembled = edgeAssembledChoices[senseN];

    // set up an array for threaded searching (same structure whether threaded or not), then clear the set
    long numQueries = edgeAssembled[edgeIndex].size();
    // NOTE: these are sorted appropriately for a set
    ScoredSeq** cbiArray = new ScoredSeq*[ numQueries + 1 ];
    long* cbiLengths = new long[ numQueries + 1 ];
    long cbiIndex = 0;
    // also create a set of the normalized wrappers for finding alignments
    for (set<ScoredSeq*>::iterator seqIt = edgeAssembled[edgeIndex].begin(); seqIt != edgeAssembled[edgeIndex].end(); ++seqIt){
      cbiArray[cbiIndex] = *seqIt;
      cbiLengths[cbiIndex] = (*seqIt)->size();
      cbiIndex++;
    }
    edgeAssembled[edgeIndex].clear();


    ScoredSeqCollectionBwt::AlignmentThreadedness mapThreads;
    if (threadedness == THREADED){ mapThreads = ScoredSeqCollectionBwt::threaded; }
    else { mapThreads = ScoredSeqCollectionBwt::notThreaded; }
    //bool* cbiKeepBool = anchorSeqCollection->hasFullMatch(numQueries, cbiArray, sense, mapThreads);

    // search using the contigs as queries - the array tracks the results
    bool* cbiKeepBool = anchorSeqCollection->hasMatch(numQueries, cbiArray, sense, cbiLengths, mapThreads,
						      ScoredSeqCollectionBwt::softMinOvl);

    // search using the contigs as database - the norm counts track the results
    // create a collection of wrappers for the new contigs. also, save time by only
    // building a database for the seqs that didn't already find a match
    set<ScoredSeq*> cbiNormSet;
    // this array is only populated at elements where the value of cbiKeepBool is false
    ScoredSeqNormalized** cbiNormArray = new ScoredSeqNormalized*[ numQueries + 1 ];
    for (long n = 0; n < numQueries; ++n){
      if (! cbiKeepBool[n]){
	ScoredSeqNormalized* cbiNorm = new ScoredSeqNormalized( cbiArray[n] );
	cbiNormArray[n] = cbiNorm;
	cbiNormSet.insert( cbiNorm );
      }
    }
    ScoredSeqCollectionBwt* newContigCollection = new ScoredSeqCollectionBwt(&cbiNormSet, dpa);
    vector<Alignment*>** matchesArray = new vector<Alignment*>*[ totalAnchors + 1 ];
    for (long n = 0; n < totalAnchors; ++n){ matchesArray[n] = new vector<Alignment*>; }
    newContigCollection->getMatches(totalAnchors, matchesArray, fullAnchorArray, sense, fullAnchorLengths, mapThreads, 
				    ScoredSeqCollectionBwt::softMinOvl);
    for (long n = 0; n < totalAnchors; ++n){
      for (vector<Alignment*>::iterator alIt = matchesArray[n]->begin(); alIt != matchesArray[n]->end(); ++alIt){
	// the found contig is seqB
	ScoredSeqNormalized* ssn = dynamic_cast<ScoredSeqNormalized*>( (*alIt)->seqB() );
	ssn->addCount();
        delete *alIt;
      }
      delete matchesArray[n];
    }
    delete [] matchesArray;
    delete newContigCollection;

    // place the contents according to the results (modifies the collections in shared memory; not threaded)
    // REMEMBER: these array elements are sorted appropriately for a set.
    set<ScoredSeq*> addToKeep;
    set<ScoredSeq*>::iterator atkIt = addToKeep.begin();
    set<ScoredSeq*> addToDiscard;
    set<ScoredSeq*>::iterator atdIt = addToDiscard.begin();
    for (long n = 0; n < numQueries; ++n){
      if ( cbiKeepBool[n] ){ atkIt = addToKeep.insert(atkIt, cbiArray[n]); }
      else {
	if (cbiNormArray[n]->getCount() > 0){ atkIt = addToKeep.insert(atkIt, cbiArray[n]); }
        else { atdIt = addToDiscard.insert(atdIt, cbiArray[n]); }
	delete cbiNormArray[n];
      }
    }
    edgeAssembled[edgeIndex].insert(addToKeep.begin(), addToKeep.end());
    localProducts->insert(addToDiscard.begin(), addToDiscard.end());
    delete [] cbiArray;
    delete [] cbiNormArray;
    delete [] cbiLengths;
    delete [] cbiKeepBool;
  }

  delete [] fullAnchorArray;
  delete [] fullAnchorLengths;
  for (set<ScoredSeq*>::iterator it = anchorSet.begin(); it != anchorSet.end(); ++it){ (*it)->deepDelete(); }
  delete anchorSeqCollection;
  delete dpa;
}



void AssemblyJobGraph::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){

    // set up the normalization sets
    createNormSets(threadedness);
    // keeps track of the products from edge jobs that don't get propagated to the
    // next edge over
    set<ScoredSeq*> allLocalProducts;

    // approach the center job from the two tips moving inwards
    long beginEdgeIndex = 0;
    long endEdgeIndex = _numNodes;
    while (beginEdgeIndex != endEdgeIndex){
      collapseOneNode(beginEdgeIndex, beginEdgeIndex + 1, &allLocalProducts, threadedness);
      beginEdgeIndex++;
      // re-check to make sure that the forward step didn't already collapse things
      if (beginEdgeIndex != endEdgeIndex){
	collapseOneNode(endEdgeIndex, endEdgeIndex - 1, &allLocalProducts, threadedness);
	endEdgeIndex--;
      }
    }

    // the assembly has finished, and the products are distributed in two places.  the ones
    // that passed the filter by matching sequences from the adjacent input bin are still in
    // _edgeSets[endEdgeIndex], aka _edgeSets[beginEdgeIndex].  the ones that were generated
    // but cast aside because they didn't match the targets are in _fullRedundantProductSet.
    // they should be combined and made non-redundant before being released.
    allLocalProducts.insert( _edgeAssembledPlus[beginEdgeIndex].begin(), _edgeAssembledPlus[beginEdgeIndex].end() );
    allLocalProducts.insert( _edgeAssembledMinus[beginEdgeIndex].begin(), _edgeAssembledMinus[beginEdgeIndex].end() );
    //_edgeSets[beginEdgeIndex].clear();

    // this set will be given nested objects, so they will no longer be sorted. however,
    // time can still be saved by having the set in which the sequences are sought and removed
    // be temporarily different from the one to which they are added.
    set<ScoredSeq*> newLocalProducts;
    for (map<ScoredSeq*,long>::iterator it = _normToLivingCount.begin(); it != _normToLivingCount.end(); ++it){
      ScoredSeqNormalized* recast = dynamic_cast<ScoredSeqNormalized*>( it->first );
      if (it->second == recast->getCount()){
	// leftover
	allLocalProducts.erase( recast );
	newLocalProducts.insert( recast->getNested() );
      } else {
	// inner seq was discarded
	_fullDiscardSet.insert( recast->getNested() );
	if (it->second != 0) {
	  allLocalProducts.erase( recast );
	  recast->resetCount();
	  recast->addCounts( it->second );
	  newLocalProducts.insert( recast->shallowCopy() );
	}
      }
      delete recast;
    }
    allLocalProducts.insert(newLocalProducts.begin(), newLocalProducts.end());

    // make the finished set non-redundant
    AssemblyJob* ajNr = _jobFactory->redundancyJob( &allLocalProducts );
    ajNr->runJob(threadedness);

    // sort the discarded sequences, deleting the removed ones
    set<ScoredSeq*> discardedSeqs;
    ajNr->discardedInputSeqs( &discardedSeqs );
    // the contents are already sorted, so can be added to this temporary set in constant time
    // and then merged to the _fullDiscardSet in linear time
    set<ScoredSeq*> addToFullDiscard;
    set<ScoredSeq*>::iterator atfdIt = addToFullDiscard.begin();
    for (set<ScoredSeq*>::iterator it = discardedSeqs.begin(); it != discardedSeqs.end(); ++it){
      if (_rawInputTotalPlus.count(*it)==0 and _rawInputTotalMinus.count(*it)==0){ (*it)->deepDelete(); }
      else { atfdIt = addToFullDiscard.insert(atfdIt, *it); }
    }
    _fullDiscardSet.insert(addToFullDiscard.begin(), addToFullDiscard.end());

    // get the final output set and define the subset that is novel - note that the "novel" sequences
    // are guaranteed to be novel, but the "remaining" sequences could be novel if they were novel
    // products from a prior assembly step
    set<ScoredSeq*> couldBeNovelSet;
    ajNr->remainingInputSeqs( &couldBeNovelSet );
    // again, the contents are already sorted, so can be added to these temporary sets in constant time
    // and then merged to the _fullDiscardSet in linear time
    set<ScoredSeq*> addToFullNovel;
    set<ScoredSeq*>::iterator atfnIt = addToFullNovel.begin();
    set<ScoredSeq*> addToFullRemain;
    set<ScoredSeq*>::iterator atfrIt = addToFullRemain.begin();
    for (set<ScoredSeq*>::iterator it = couldBeNovelSet.begin(); it != couldBeNovelSet.end(); ++it){
      if (_rawInputTotalPlus.count(*it)==0 and _rawInputTotalMinus.count(*it)==0){ atfnIt = addToFullNovel.insert(atfnIt, *it); }
      else { atfrIt = addToFullRemain.insert(atfrIt, *it); }
    }
    _fullRemainingSet.insert(addToFullRemain.begin(), addToFullRemain.end());
    // the novel set gets those plus the truly novel sequences
    ajNr->newlyAssembledSeqs( &_fullNovelSet );
    _fullNovelSet.insert(addToFullNovel.begin(), addToFullNovel.end());
    delete ajNr;
    _wasRun = true;
  }
}

bool AssemblyJobGraph::wasRun(){ return _wasRun; }

void AssemblyJobGraph::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _rawInputTotalPlus.begin(), _rawInputTotalPlus.end() );
  outsideSet->insert( _rawInputTotalMinus.begin(), _rawInputTotalMinus.end() );
}
void AssemblyJobGraph::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _fullRemainingSet.begin(), _fullRemainingSet.end() );
}
void AssemblyJobGraph::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _fullDiscardSet.begin(), _fullDiscardSet.end() );
}
void AssemblyJobGraph::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _fullNovelSet.begin(), _fullNovelSet.end() );
}
void AssemblyJobGraph::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _fullRemainingSet.begin(), _fullRemainingSet.end() );
  outsideSet->insert( _fullNovelSet.begin(), _fullNovelSet.end() );
}



void AssemblyJobGraph::OK(){}



#endif
