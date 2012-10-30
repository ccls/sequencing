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

#ifndef ASSEMBLYJOBSSDEBRUIJN_CPP
#define ASSEMBLYJOBSSDEBRUIJN_CPP

#include "AssemblyJobSsDeBruijn.h"
// just for tests!
#include "ScoredSeqShallow.h"
#include "ScoredSeqMonoScore.h"

#include "AssemblyException.h"
#include <queue>
#include <iostream>
#include <omp.h>
using namespace::std;

AssemblyJobSsDeBruijn::AssemblyJobSsDeBruijn(){}

AssemblyJobSsDeBruijn::AssemblyJobSsDeBruijn(set<ScoredSeq*>* seqs, ParamsDeBruijn* paramsDeBruijn) :
  _kmerSize(paramsDeBruijn->kmerSize()),
  _wasRun(false){
  _paramsDeBruijn = new ParamsDeBruijn(paramsDeBruijn);
  _inputSeqs.insert( seqs->begin(), seqs->end() );
  for (set<ScoredSeq*>::iterator it = seqs->begin(); it != seqs->end(); ++it){
    if ( paramsDeBruijn->inputToDeBruijn(*it) ){ _usedSeqs.insert(*it); }
    else { _remainingSeqs.insert(*it); }
  }
}
AssemblyJobSsDeBruijn::~AssemblyJobSsDeBruijn(){
  delete _paramsDeBruijn;
}



void AssemblyJobSsDeBruijn::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJNull::makeShallow can't be called after the job was run.");
  }
  set<ScoredSeq*> shallowCopies;
  // sequences will have to be re-tested for paramsDeBruijn compliance
  _usedSeqs.clear();
  _remainingSeqs.clear();

  // buffer all input, then copy
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    (*it)->bottomBuffer();
  }
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    shallowCopies.insert( (*it)->shallowCopy() );
    (*it)->deepUnbuffer();
  }
  if (shallowDelete){
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      delete (*it);
    }
  }
  // replace the input contents with the shallow copies
  _inputSeqs.clear();
  _inputSeqs.insert( shallowCopies.begin(), shallowCopies.end() );
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    if ( _paramsDeBruijn->inputToDeBruijn(*it) ){ _usedSeqs.insert(*it); }
    else { _remainingSeqs.insert(*it); }
  }
}
void AssemblyJobSsDeBruijn::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){

    // FIRST build the tree structure (this may be thread-able)
    pair<AssemblyJobSsDeBruijn::NucNode*,long> topNodeAndCount;
    if (threadedness == AssemblyJob::THREADED){
      // break up the query set
      int numThreads = omp_get_max_threads();
      set<ScoredSeq*>* threadSets = new set<ScoredSeq*>[ numThreads + 1 ];
      long countIndex = 0;
      for (set<ScoredSeq*>::iterator it = _usedSeqs.begin(); it != _usedSeqs.end(); ++it){
	threadSets[countIndex % numThreads].insert( *it );
	countIndex++;
      }
      pair<AssemblyJobSsDeBruijn::NucNode*,long>* tnacArray = new pair<AssemblyJobSsDeBruijn::NucNode*,long>[ numThreads + 1 ];
      #pragma omp parallel for schedule(dynamic)
      for (int n = 0; n < numThreads; ++n){
	tnacArray[n] = makeGraphTree(&threadSets[n], _kmerSize);
      }

      // thread the combining of trees as much as possible (log collapse)
      int numOldTrees = numThreads;
      int numNewTrees = (numThreads + 1) / 2;
      while (numOldTrees > 1){
	// first, set up the new tree array with the even nodes
	pair<AssemblyJobSsDeBruijn::NucNode*,long>* tempTnacArray = new pair<AssemblyJobSsDeBruijn::NucNode*,long>[ numNewTrees + 1 ];
	for (int n = 0; n < numOldTrees; ++n){
	  if (n % 2 == 0){ tempTnacArray[n/2] = tnacArray[n]; }
	}
	// now, add the contents of the odd nodes to the even nodes (threaded)
        #pragma omp parallel for schedule(dynamic)
	for (int n = 0; n < numNewTrees; ++n){
	  int secondIndex = n*2 + 1;
	  if (secondIndex < numOldTrees){
	    tempTnacArray[n].second += combineTwoTrees(tempTnacArray[n].first, tnacArray[secondIndex].first, _kmerSize);
	    (tnacArray[secondIndex].first)->deepDelete();
	  }
	}
	// adjust the count and re-set the pointer to the new array
	numOldTrees = numNewTrees;
	numNewTrees = (numOldTrees + 1) / 2;
	delete [] tnacArray;
	tnacArray = tempTnacArray;
      }

      topNodeAndCount = tnacArray[0];
      delete [] threadSets;
      delete [] tnacArray;
    } else {
      topNodeAndCount = makeGraphTree(&_usedSeqs, _kmerSize);
    }

    // SECOND isolate the de Bruijn graph that exists at the tips of the tree branches
    DeBruijnGraph* graph = new DeBruijnGraph(topNodeAndCount.first, topNodeAndCount.second, _kmerSize);

    // THIRD pick out a best edge on each side of each kmer node
    findBestLinks(graph);

  // FIFTH construct paths to form contigs
  // and concurrently
  // SIXTH post-filter the output seqs according to the params
    makeSeqsFromGraph(graph, &_novelSeqs);
    delete graph;
    topNodeAndCount.first->deepDelete();
    _wasRun = true;
  }
}




bool AssemblyJobSsDeBruijn::wasRun(){ return _wasRun; }

void AssemblyJobSsDeBruijn::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}
void AssemblyJobSsDeBruijn::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _remainingSeqs.begin(), _remainingSeqs.end() );
}
void AssemblyJobSsDeBruijn::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _usedSeqs.begin(), _usedSeqs.end() );
}
void AssemblyJobSsDeBruijn::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
}
void AssemblyJobSsDeBruijn::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _novelSeqs.begin(), _novelSeqs.end() );
  outsideSet->insert( _remainingSeqs.begin(), _remainingSeqs.end() );
}


// this top method is just an organized launcher
long AssemblyJobSsDeBruijn::combineTwoTrees(SeqNode* topNodeA, SeqNode* topNodeB, int kmerSize){
  return combineTreesMidHelp(topNodeA, kmerSize, topNodeA, topNodeB, kmerSize);
}
long AssemblyJobSsDeBruijn::combineTreesMidHelp(SeqNode* topNodeA, int kmerSize,
						SeqNode* branchA, SeqNode* branchB, int remainingKmer){
  char nucs[] = {'A','T','C','G'};
  long newNodeSum = 0;
  if (remainingKmer == 1){
    for (int n = 0; n < 4; ++n){
      DeBruijnNode* newNodeB = dynamic_cast<DeBruijnNode*>(branchB->getBranch(nucs[n]));
      if (newNodeB != 0){
	DeBruijnNode* newNodeA = dynamic_cast<DeBruijnNode*>(branchA->getBranch(nucs[n]));
	if (newNodeA == 0){
	  newNodeA = new DeBruijnNode(branchA,nucs[n]);
	  branchA->addBranch( newNodeA );
	  newNodeSum += 1;
	}
	newNodeA->addKmerScore( newNodeB->getKmerScore() );
	newNodeSum += combineTreesBottomHelp(topNodeA, kmerSize, newNodeA, newNodeB);
      }
    }
  } else {
    for (int n = 0; n < 4; ++n){
      SeqNode* newNodeB = branchB->getBranch(nucs[n]);
      if (newNodeB != 0){
	SeqNode* newNodeA = branchA->getBranch(nucs[n]);
	if (newNodeA == 0){
	  newNodeA = new NucNode(branchA,nucs[n]);
	  branchA->addBranch( newNodeA );
	}
	newNodeSum += combineTreesMidHelp(topNodeA, kmerSize, newNodeA, newNodeB, remainingKmer-1);
      }
    }
  }
  return newNodeSum;
}
long AssemblyJobSsDeBruijn::combineTreesBottomHelp(SeqNode* topNodeA, int kmerSize,
						   DeBruijnNode* tipA, DeBruijnNode* tipB){
  long newNodeSum = 0;
  char nucs[] = {'A','T','C','G'};
  // perform this operation for 5p links...
  for (int n = 0; n < 4; ++n){
    DeBruijnNode* linkB = tipB->get5pLink(nucs[n]);
    if (linkB == 0){} // nothing needs to happen, linkA stays the same
    else {
      DeBruijnNode* linkA = tipA->get5pLink(nucs[n]);
      // make sure that the node exists, then add the counts
      if (linkA == 0){
	pair<DeBruijnNode*,long> newLink = combineTreesAddNodeHelp(topNodeA, linkB, kmerSize);
	linkA = newLink.first;
	newNodeSum += newLink.second;
      }
      tipA->add5pLink(linkA, tipB->get5pLinkScore(nucs[n]));
    }
  }
  // ...and now the same for 3p links
  for (int n = 0; n < 4; ++n){
    DeBruijnNode* linkB = tipB->get3pLink(nucs[n]);
    if (linkB == 0){} // nothing needs to happen, linkA stays the same
    else {
      DeBruijnNode* linkA = tipA->get3pLink(nucs[n]);
      // make sure that the node exists, then add the counts
      if (linkA == 0){
	pair<DeBruijnNode*,long> newLink = combineTreesAddNodeHelp(topNodeA, linkB, kmerSize);
	linkA = newLink.first;
	newNodeSum += newLink.second;
      }
      tipA->add3pLink(linkA, tipB->get3pLinkScore(nucs[n]));
    }
  }
  return newNodeSum;
}
pair<AssemblyJobSsDeBruijn::DeBruijnNode*,long> AssemblyJobSsDeBruijn::combineTreesAddNodeHelp(SeqNode* topNodeA, DeBruijnNode* bottomNodeB, int kmerSize){
  char* kmer = new char[kmerSize+1];
  bottomNodeB->getKmer(kmer,0);
  SeqNode* currentNode = topNodeA;
  for (int n = 0; n < kmerSize-1; ++n){
    SeqNode* newNode = currentNode->getBranch(kmer[n]);
    if (newNode == 0){
      newNode = new NucNode(currentNode,kmer[n]);
      currentNode->addBranch( newNode );
    }
    currentNode = newNode;
  }
  DeBruijnNode* finishedNode = dynamic_cast<DeBruijnNode*>( currentNode->getBranch(kmer[kmerSize-1]) );
  long newNodeCount = 0;
  if (finishedNode == 0){
    finishedNode = new DeBruijnNode(currentNode, kmer[kmerSize-1]);
    currentNode->addBranch(finishedNode);
    newNodeCount = 1;
  }
  delete [] kmer;
  return pair<AssemblyJobSsDeBruijn::DeBruijnNode*,long>(finishedNode,newNodeCount);
}



void AssemblyJobSsDeBruijn::testTreeConstruction(set<ScoredSeq*>* input, int kmerSize, set<ScoredSeq*>* foundKmers){
  pair<NucNode*,long> topNodeAndNum = makeGraphTree(input,kmerSize);
  SeqNode* topNode = topNodeAndNum.first;
  // in case the input set is not empty
  long numKmers = foundKmers->size();
  testTreeConstructionHelper( foundKmers, topNode, kmerSize );
  // make sure that the correct number of base nodes were reported back
  numKmers = foundKmers->size() - numKmers;
  if (numKmers != topNodeAndNum.second){ throw AssemblyException::LogicError("AJSDB::testTC wrong number of nodes reported."); }
  topNode->deepDelete();
}
void AssemblyJobSsDeBruijn::testTreeConstructionHelper(set<ScoredSeq*>* foundKmers, SeqNode* node, int stepsToBottom){
  if (stepsToBottom == 0){
    DeBruijnNode* dbNode = dynamic_cast<DeBruijnNode*>(node);
    foundKmers->insert( new ScoredSeqMonoScore(dbNode->getKmer(), dbNode->getKmerScore()) );
    //foundKmers->insert( ScoredSeq::getScoredSeq(dbNode->getKmer(), dbNode->getKmerScore()) );
  } else {
    char nucList[] = {'A','T','C','G'};
    for (int n = 0; n < 4; ++n){
      SeqNode* nextNode = node->getBranch( nucList[n] );
      if (nextNode != 0){ testTreeConstructionHelper( foundKmers, nextNode, stepsToBottom-1 ); }
    }
  }
}


void AssemblyJobSsDeBruijn::findBestLinks(DeBruijnGraph* graph){
  char nucs[] = {'A','T','C','G'};
  for (DeBruijnGraph::Iterator nodeIt = graph->begin(); nodeIt != graph->end(); ++nodeIt){
    DeBruijnNode* currentNode = *nodeIt;
    DeBruijnNode* best5pLink = 0;
    float best5pLinkScore = 0;
    DeBruijnNode* best3pLink = 0;
    float best3pLinkScore = 0;
    for (int n = 0; n < 4; ++n){
      DeBruijnNode* cand5pLink = currentNode->get5pLink(nucs[n]);
      if (cand5pLink != 0){
	// take the worse of the two linkage scores
	float next5pLinkScore;
	if (currentNode->get5pLinkScore(nucs[n]) > cand5pLink->get3pLinkScore(currentNode->get3pNuc()) ){
	  next5pLinkScore = cand5pLink->get3pLinkScore(currentNode->get3pNuc());
	} else { next5pLinkScore = currentNode->get5pLinkScore(nucs[n]); }
	if (next5pLinkScore > best5pLinkScore){
	  best5pLinkScore = next5pLinkScore;
	  best5pLink = cand5pLink;
	}
      }
      DeBruijnNode* cand3pLink = currentNode->get3pLink(nucs[n]);
      if (cand3pLink != 0){
	// take the worse of the two linkage scores
	float next3pLinkScore;
	if (currentNode->get3pLinkScore(nucs[n]) > cand3pLink->get5pLinkScore(currentNode->get5pNuc()) ){
	  next3pLinkScore = cand3pLink->get5pLinkScore(currentNode->get5pNuc());
	} else { next3pLinkScore = currentNode->get3pLinkScore(nucs[n]); }
	if (next3pLinkScore > best3pLinkScore){
	  best3pLinkScore = next3pLinkScore;
	  best3pLink = cand3pLink;
	}
      }
    }
    currentNode->setBest5pLink(best5pLink);
    currentNode->setBest3pLink(best3pLink);
  }

  // FOURTH resolve best-edge conflicts
  for (DeBruijnGraph::Iterator nodeIt = graph->begin(); nodeIt != graph->end(); ++nodeIt){
    DeBruijnNode* currentNode = *nodeIt;
    DeBruijnNode* best5pLink = currentNode->best5pLink();
    if (best5pLink != 0 and best5pLink->best3pLink() != currentNode){
      DeBruijnNode* bestRecip = best5pLink->best3pLink();
      if (bestRecip == 0){ best5pLink->setBest3pLink(currentNode); }
      else { currentNode->setBest5pLink( 0 ); }
    }
    DeBruijnNode* best3pLink = currentNode->best3pLink();
    if (best3pLink != 0 and best3pLink->best5pLink() != currentNode){
      DeBruijnNode* bestRecip = best3pLink->best5pLink();
      if (bestRecip == 0){ best3pLink->setBest5pLink(currentNode); }
      else { currentNode->setBest3pLink( 0 ); }
    }
  }
}



void AssemblyJobSsDeBruijn::makeSeqsFromGraph(DeBruijnGraph* graph, set<ScoredSeq*>* productSeqs){
  graph->sortNodes();
  for (DeBruijnGraph::Iterator nodeIt = graph->begin(); nodeIt != graph->end(); ++nodeIt){
    // some nodes will already have been traversed
    if (! (*nodeIt)->wasVisited() ){
      (*nodeIt)->visit();
      DeBruijnNode* legitNode5p = *nodeIt;
      DeBruijnNode* legitNode3p = *nodeIt;
      DeBruijnNode* tempNode5p;
      DeBruijnNode* tempNode3p;
      bool active5p = true;
      bool active3p = true;
      long numNodes = 1;
      while (active5p or active3p){
	if (active5p){
	  tempNode5p = legitNode5p->best5pLink();
	  if (tempNode5p == 0 or tempNode5p->wasVisited()){ active5p = false; }
	  else {
	    legitNode5p = tempNode5p;
	    legitNode5p->visit();
	    numNodes++;
	  }
	}
	if (active3p){
	  tempNode3p = legitNode3p->best3pLink();
	  if (tempNode3p == 0 or tempNode3p->wasVisited()){ active3p = false; }
	  else {
            legitNode3p = tempNode3p;
	    legitNode3p->visit();
	    numNodes++;
	  }
	}
      }
      ScoredSeq* novelSeq = makeSeqFromPath(legitNode5p,legitNode3p,numNodes);
      if ( _paramsDeBruijn->retainDeBruijnOutput(novelSeq) ){ productSeqs->insert( novelSeq ); }
      else { delete novelSeq; }
    }
  }
}


ScoredSeq* AssemblyJobSsDeBruijn::makeSeqFromPath(DeBruijnNode* node5p, DeBruijnNode* node3p, long numNodes){
  long seqLen = _kmerSize + numNodes - 1;
  char* cSeq = new char[seqLen + 1];
  cSeq[seqLen] = '\0';
  float* scores = new float[seqLen];
  float* links = new float[seqLen-1];

  // do the initial kmer and any overlapping positions
  node5p->getKmer(cSeq,0);
  float firstScore = node5p->getKmerScore();
  for (int n = 0; n < _kmerSize; ++n){
    scores[n] = firstScore;
    if (n < _kmerSize - 1){ links[n] = firstScore; }
  }

  // fill in the rest
  DeBruijnNode* currentNode = node5p;
  DeBruijnNode* priorNode = 0;
  long pos5p = 0;
  long pos3p = _kmerSize - 1;
  while (currentNode != node3p){
    // increment
    pos5p++;
    pos3p++;
    priorNode = currentNode;
    currentNode = currentNode->best3pLink();
    // fill in values
    cSeq[pos3p] = currentNode->get3pNuc();
    scores[pos3p] = currentNode->getKmerScore();
    links[pos3p-1] = priorNode->get3pLinkScore( currentNode->get3pNuc() );
    if (links[pos5p-1] < currentNode->get5pLinkScore( cSeq[pos5p-1] )){
      links[pos5p-1] = currentNode->get5pLinkScore( cSeq[pos5p-1] );
    }
    for (long pos = pos5p; pos < pos3p; ++pos){
      if (scores[pos] < currentNode->getKmerScore()){ scores[pos] = currentNode->getKmerScore(); }
    }
  }

  ScoredSeq* novelSeq = new ScoredSeqShallow(true, cSeq, scores, links, seqLen);
  //ScoredSeq* novelSeq = ScoredSeq::getScoredSeq(cSeq, scores, links, seqLen);
  /*
  delete [] cSeq;
  delete [] scores;
  delete [] links;
  */
  return novelSeq;
}



pair<AssemblyJobSsDeBruijn::NucNode*,long> AssemblyJobSsDeBruijn::makeGraphTree(set<ScoredSeq*>* seqs, int kmerSize){
  NucNode* topNode = new NucNode(0, '\0');
  long numBaseNodes = 0;
  if (kmerSize > 0){
    for (set<ScoredSeq*>::iterator seqIt = seqs->begin(); seqIt != seqs->end(); ++seqIt){
      ScoredSeq* seq = *seqIt;

      // these are the nodes for which the current nuc is at the Nth position in the kmer;
      // the elements earlier in the array represent kmers that come later in the sequence
      SeqNode* currentNodes[kmerSize-1];
      float worstScores[kmerSize-1];
      // for making the de Bruijn graph at the tips
      DeBruijnNode* finishedNode = 0;
      DeBruijnNode* priorNode = 0;

      // used to determine how far along the kmer construction has made it
      long maxCurrentNodeIndex = 0;
      bool createFinishedNode = false;

      char* nucs = seq->getSeq('+');
      float* scores = seq->getScores('+');
      float* links = seq->getLinks('+');
      long seqSize = seq->size();
      long seqSizeM1 = seqSize - 1;

      for (long pos = 0; pos < seqSize; ++pos){
	char nuc = nucs[pos];
	float nucScore = scores[pos];
	// will only be used up to the second-to-last position
	float linkScore;
	if (pos < seqSizeM1){ linkScore = links[pos]; }
	else { linkScore = 0; }

	// re-set the derivation of kmers if the nucleotide is an N
	if (nuc == 'N'){
	  maxCurrentNodeIndex = 0;
	  createFinishedNode = false;
	  finishedNode = 0;
	  priorNode = 0;
	} else {
	  priorNode = finishedNode;
	  if (createFinishedNode){
	    SeqNode* oldNode = currentNodes[maxCurrentNodeIndex];
	    float worstScore = worstScores[maxCurrentNodeIndex];
	    finishedNode = dynamic_cast<DeBruijnNode*>( oldNode->getBranch(nuc) );
	    if (finishedNode == 0){
	      finishedNode = new DeBruijnNode(oldNode, nuc);
	      oldNode->addBranch(finishedNode);
	      numBaseNodes++;
	    }
	    // add the kmer score to the node; the worst score possible!
	    if ( nucScore < worstScore ){ worstScore = nucScore; }
	    finishedNode->addKmerScore(worstScore);

	    // add the 3p-directed link score
	    if (priorNode != 0){
	      float link3pScore = links[pos-1];
	      float link5pScore = links[pos-kmerSize];
	      float worseLinkScore = link3pScore;
	      if (link5pScore < worseLinkScore){ worseLinkScore = link5pScore; }
	      priorNode->add3pLink(finishedNode,worseLinkScore);
	      finishedNode->add5pLink(priorNode,worseLinkScore);
	    }
	  }
	  for (int n = maxCurrentNodeIndex; n >= 0; --n){
	    SeqNode* oldNode;
	    float worstScore;
	    if (n==0){
	      oldNode = topNode;
	      worstScore = nucScore;
	    }
	    else {
	      oldNode = currentNodes[n-1];
	      worstScore = worstScores[n-1];
	      if ( nucScore < worstScore ){ worstScore = nucScore; }
	    }
	    if ( linkScore < worstScore ){ worstScore = linkScore; }
	    SeqNode* newNode = oldNode->getBranch(nuc);
	    if (newNode == 0){
	      newNode = new NucNode(oldNode,nuc);
	      oldNode->addBranch( newNode );
	    }
	    currentNodes[n] = newNode;
	    worstScores[n] = worstScore;
	  }
	  if (! createFinishedNode){
	    if (maxCurrentNodeIndex < kmerSize-2){ maxCurrentNodeIndex++; }
	    else { createFinishedNode = true; }
	  }
	}
      }
      delete [] nucs;
      delete [] scores;
      delete [] links;
    }
  }
  return pair<NucNode*,long>(topNode,numBaseNodes);
}





void AssemblyJobSsDeBruijn::OK(){}



// implement the virtual destructor
AssemblyJobSsDeBruijn::SeqNode::~SeqNode(){}

// DOES NOT REFERENCE DATA FIELDS; JUST USES PUBLIC METHODS
string AssemblyJobSsDeBruijn::SeqNode::getKmer(){
  // first, figure out how big the kmer is
  int kmerSize = 0;
  SeqNode* currentNode = this;
  while (currentNode != 0 and currentNode->get3pNuc() != '\0'){
    kmerSize++;
    currentNode = currentNode->getParent();
  }
  // second, fill in the kmer
  char seq[kmerSize+1];
  seq[kmerSize] = '\0';
  currentNode = this;
  for (int n = kmerSize-1; n >= 0; --n){
    seq[n] = currentNode->get3pNuc();
    currentNode = currentNode->getParent();
  }
  return string(seq);
}
void AssemblyJobSsDeBruijn::SeqNode::getKmer(char* seqToFillIn, long startCoord){
  // first, figure out how big the kmer is
  int kmerSize = 0;
  SeqNode* currentNode = this;
  while (currentNode != 0 and currentNode->get3pNuc() != '\0'){
    kmerSize++;
    currentNode = currentNode->getParent();
  }
  // second, fill in the kmer
  currentNode = this;
  for (int n = startCoord+kmerSize-1; n >= startCoord; --n){
    seqToFillIn[n] = currentNode->get3pNuc();
    currentNode = currentNode->getParent();
  }
}
char AssemblyJobSsDeBruijn::SeqNode::get5pNuc(){
  // first, figure out how big the kmer is
  char currentNuc = '\0';
  SeqNode* currentNode = this;
  while (currentNode != 0 and currentNode->get3pNuc() != '\0'){
    currentNuc = currentNode->get3pNuc();
    currentNode = currentNode->getParent();
  }
  return currentNuc;
}



AssemblyJobSsDeBruijn::NucNode::NucNode(SeqNode* parent, char nuc) :
  _parent(parent),
  _nuc(nuc)
{
  for (int n = 0; n < 4; ++n){ _branches[n] = 0; }
}
AssemblyJobSsDeBruijn::NucNode::~NucNode(){}
void AssemblyJobSsDeBruijn::NucNode::deepDelete(){
  for (int n = 0; n < 4; ++n){ if(_branches[n] != 0){ _branches[n]->deepDelete(); } }
  delete this;
}

AssemblyJobSsDeBruijn::SeqNode* AssemblyJobSsDeBruijn::NucNode::getParent(){ return _parent; }
char AssemblyJobSsDeBruijn::NucNode::get3pNuc(){ return _nuc; }

bool AssemblyJobSsDeBruijn::NucNode::isBranchable(){ return true; }
int AssemblyJobSsDeBruijn::NucNode::getNucIndex(char nuc){
  int nucIndex;
  switch ( nuc ){
  case 'A': nucIndex = 0; break;
  case 'C': nucIndex = 1; break;
  case 'G': nucIndex = 2; break;
  case 'T': nucIndex = 3; break;
  case '\0':
    throw AssemblyException::ArgError("NucNode cannot have a null nuc char");
    break;
  default:
    throw AssemblyException::ArgError("invalid nuc input for NucNode");
  }
  return nucIndex;
}
void AssemblyJobSsDeBruijn::NucNode::addBranch(SeqNode* newNode){
  int nucIndex = getNucIndex( newNode->get3pNuc() );
  if (_branches[nucIndex] != 0){
    throw AssemblyException::CallingError("can't add a branch that already exists in NucNode");
  } else {
    _branches[nucIndex] = newNode;
  }
}
AssemblyJobSsDeBruijn::SeqNode* AssemblyJobSsDeBruijn::NucNode::getBranch(char nuc){
  return _branches[ getNucIndex( nuc ) ]; }





AssemblyJobSsDeBruijn::DeBruijnNode::DeBruijnNode(SeqNode* parent, char nuc) :
  _parent(parent),
  _nuc(nuc),
  _wasVisited(false),
  _kmerScore(0),
  _best5pLink(0),
  _best3pLink(0)
{
  for (int n = 0; n < 4; ++n){
    _5pLinks[n] = 0;
    _3pLinks[n] = 0;
    _5pLinkScores[n] = 0;
    _3pLinkScores[n] = 0;
  }
}
AssemblyJobSsDeBruijn::DeBruijnNode::~DeBruijnNode(){}
void AssemblyJobSsDeBruijn::DeBruijnNode::deepDelete(){
  delete this;
}

AssemblyJobSsDeBruijn::SeqNode* AssemblyJobSsDeBruijn::DeBruijnNode::getParent(){ return _parent; }
char AssemblyJobSsDeBruijn::DeBruijnNode::get3pNuc(){ return _nuc; }
bool AssemblyJobSsDeBruijn::DeBruijnNode::isBranchable(){ return false; }
void AssemblyJobSsDeBruijn::DeBruijnNode::addBranch(SeqNode* newNode){
  throw AssemblyException::ArgError("DeBruijnNode cannot have a branch added, it is not branchable");
}
AssemblyJobSsDeBruijn::SeqNode* AssemblyJobSsDeBruijn::DeBruijnNode::getBranch(char nuc){
  throw AssemblyException::ArgError("DeBruijnNode cannot return a branch, it is not branchable");
}
int AssemblyJobSsDeBruijn::DeBruijnNode::getNucIndex(char nuc){
  int nucIndex;
  switch ( nuc ){
  case 'A': nucIndex = 0; break;
  case 'C': nucIndex = 1; break;
  case 'G': nucIndex = 2; break;
  case 'T': nucIndex = 3; break;
  case '\0':
    throw AssemblyException::ArgError("NucNode cannot have a null nuc char");
    break;
  default:
    throw AssemblyException::ArgError("invalid nuc input for NucNode");
  }
  return nucIndex;
}
void AssemblyJobSsDeBruijn::DeBruijnNode::add5pLink(DeBruijnNode* node, float linkScore){
  int nucIndex = getNucIndex( node->get5pNuc() );
  if (_5pLinks[nucIndex] == 0){ _5pLinks[nucIndex] = node; }
  _5pLinkScores[nucIndex] += linkScore;
}
void AssemblyJobSsDeBruijn::DeBruijnNode::add3pLink(DeBruijnNode* node, float linkScore){
  int nucIndex = getNucIndex( node->get3pNuc() );
  if (_3pLinks[nucIndex] == 0){ _3pLinks[nucIndex] = node; }
  _3pLinkScores[nucIndex] += linkScore;
}
AssemblyJobSsDeBruijn::DeBruijnNode* AssemblyJobSsDeBruijn::DeBruijnNode::get5pLink(char nuc){
  return _5pLinks[ getNucIndex(nuc) ]; }
AssemblyJobSsDeBruijn::DeBruijnNode* AssemblyJobSsDeBruijn::DeBruijnNode::get3pLink(char nuc){
  return _3pLinks[ getNucIndex(nuc) ]; }
float AssemblyJobSsDeBruijn::DeBruijnNode::get5pLinkScore(char nuc){
  return _5pLinkScores[ getNucIndex(nuc) ]; }
float AssemblyJobSsDeBruijn::DeBruijnNode::get3pLinkScore(char nuc){
  return _3pLinkScores[ getNucIndex(nuc) ]; }
    // properties of the kmer overall
float AssemblyJobSsDeBruijn::DeBruijnNode::getKmerScore(){ return _kmerScore; }
void AssemblyJobSsDeBruijn::DeBruijnNode::addKmerScore(float score){ _kmerScore += score; }

void AssemblyJobSsDeBruijn::DeBruijnNode::visit(){ _wasVisited = true; }
bool AssemblyJobSsDeBruijn::DeBruijnNode::wasVisited(){ return _wasVisited; }


void AssemblyJobSsDeBruijn::DeBruijnNode::setBest5pLink(DeBruijnNode* link){ _best5pLink = link; }
void AssemblyJobSsDeBruijn::DeBruijnNode::setBest3pLink(DeBruijnNode* link){ _best3pLink = link; }
AssemblyJobSsDeBruijn::DeBruijnNode* AssemblyJobSsDeBruijn::DeBruijnNode::best5pLink(){ return _best5pLink; }
AssemblyJobSsDeBruijn::DeBruijnNode* AssemblyJobSsDeBruijn::DeBruijnNode::best3pLink(){ return _best3pLink; }




AssemblyJobSsDeBruijn::DeBruijnGraph::DeBruijnGraph(){}
AssemblyJobSsDeBruijn::DeBruijnGraph::DeBruijnGraph(SeqNode* topNode, long nodeCount, int kmerSize) :
  _topNode(topNode),
  _nodeCount(nodeCount),
  _kmerSize(kmerSize)
{
  if (nodeCount == 0){ _nodeArray = new DeBruijnNode*[1]; }
  else {
    _nodeArray = new DeBruijnNode*[nodeCount];
    long currentIndex = 0;
    currentIndex = nodeArrayHelper(topNode, kmerSize, currentIndex);
    if (currentIndex != nodeCount){
      throw AssemblyException::LogicError("AJSDB::DBG wrong num nodes found.");
    }
  }
}

long AssemblyJobSsDeBruijn::DeBruijnGraph::nodeArrayHelper(SeqNode* node, int stepsToBottom, long currentIndex){
  if (stepsToBottom == 0){
    _nodeArray[currentIndex] = dynamic_cast<DeBruijnNode*>(node);
    return currentIndex + 1;
  } else {
    char nucList[] = {'A','T','C','G'};
    for (int n = 0; n < 4; ++n){
      SeqNode* nextNode = node->getBranch( nucList[n] );
      if (nextNode != 0){ currentIndex = nodeArrayHelper(nextNode, stepsToBottom-1, currentIndex); }
    }
    return currentIndex;
  }
}


AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator AssemblyJobSsDeBruijn::DeBruijnGraph::begin(){ return Iterator(this,0); }
AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator AssemblyJobSsDeBruijn::DeBruijnGraph::end(){ return Iterator(this,_nodeCount); }

AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::Iterator(){}
AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::Iterator(DeBruijnGraph* source, long nodeNum) : _nodeNum(nodeNum), _source(source){}
AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::Iterator(const Iterator& it) : _source(it._source), _nodeNum(it._nodeNum){}
AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::operator++(){ _nodeNum++; return *this; } 
bool AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::operator==(const Iterator& it){ return _nodeNum==it._nodeNum and _source==it._source; }
bool AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::operator!=(const Iterator& it){ return _nodeNum!=it._nodeNum or _source!=it._source; }
AssemblyJobSsDeBruijn::DeBruijnNode* AssemblyJobSsDeBruijn::DeBruijnGraph::Iterator::operator*() {
  return _source->_nodeArray[ _nodeNum ]; }

AssemblyJobSsDeBruijn::DeBruijnGraph::~DeBruijnGraph(){
  delete [] _nodeArray;
}

void AssemblyJobSsDeBruijn::DeBruijnGraph::sortNodes(){
  if (_nodeCount > 1){
    // this is going to make sortNodes O( n * log(m) ) where n in the num of nodes but m is the num
    // of unique scores, which is smaller.
    map<float,long> scoreToCount;
    // figure out how many instances of each score there are (many scores will be 1 or 0)
    for (long n = 0; n < _nodeCount; ++n){
      float score = _nodeArray[n]->getKmerScore();
      map<float,long>::iterator scoreBin = scoreToCount.find( score );
      if (scoreBin == scoreToCount.end()){ scoreToCount.insert( pair<float,long>(score,1) ); }
      else { scoreBin->second++; }
    }

    // figure out what index each score will start at and set up an array of counts for the score bin
    long numScores = scoreToCount.size();
    long rankToStart[numScores];
    long rankToCount[numScores];
    map<float,long> scoreToRank;
    long rank = numScores;
    long start = _nodeCount;
    for (map<float,long>::iterator it = scoreToCount.begin(); it != scoreToCount.end(); ++it){
      rank--;
      start -= it->second;
      rankToStart[rank] = start;
      rankToCount[rank] = 0;
      scoreToRank.insert( pair<float,long>( it->first, rank ) );
    }

    // now place the nodes appropriately in a new array
    DeBruijnNode** sortedArray = new DeBruijnNode*[_nodeCount];
    for (long n = 0; n < _nodeCount; ++n){
      long rank = scoreToRank.find( _nodeArray[n]->getKmerScore() )->second;
      sortedArray[rankToStart[rank] + rankToCount[rank]] = _nodeArray[n];
      rankToCount[rank]++;
    }

    // now replace the old array with the new array
    delete [] _nodeArray;
    _nodeArray = sortedArray;
  }
}


#endif
