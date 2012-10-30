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

#ifndef ASSEMBLYGRAPHNODE_CPP
#define ASSEMBLYGRAPHNODE_CPP

#include "AssemblyGraphNode.h"
#include "AssemblyException.h"
#include <iostream>
using namespace::std;



//AssemblyGraphNode::AssemblyGraphNode(){}
// REQUIRES: SSFlip senses must be worked out so that the two seq sets are concordant in terms of their sense.
// (would jive if combined into a single-stranded assembly job)
// REQUIRES: one of the SSFlip objects in each edge collection wraps the centerContig
AssemblyGraphNode::AssemblyGraphNode(ScoredSeq* centerContig, set<ScoredSeqFlip*>* plusEdgeSeqs, set<ScoredSeqFlip*>* minusEdgeSeqs) :
  _centerContig(centerContig),
  _centerIsPlus(true),
  _wasVisited(false)
{
  for (set<ScoredSeqFlip*>::iterator it = plusEdgeSeqs->begin(); it != plusEdgeSeqs->end(); ++it){
    _plusEdgeSeqs.insert( ScoredSeqFlip::getFlip( (*it)->getNested(), (*it)->getSense() ) );
  }
  for (set<ScoredSeqFlip*>::iterator it = minusEdgeSeqs->begin(); it != minusEdgeSeqs->end(); ++it){
    _minusEdgeSeqs.insert( ScoredSeqFlip::getFlip( (*it)->getNested(), (*it)->getSense() ) );
  }

  // initial linkage values are null
  _plusLink = NULL;
  _minusLink = NULL;
  _plusLinkScore = 0.0;
  _minusLinkScore = 0.0;
  _plusLinkAgrees = true;
  _minusLinkAgrees = true;

  // run a check; this can be commented out for speed
  int plusCount = 0;
  int minusCount = 0;
  for (set<ScoredSeqFlip*>::iterator it = _plusEdgeSeqs.begin(); it != _plusEdgeSeqs.end(); ++it){
    if ( (*it)->getNested()==centerContig ){
      if ( (*it)->getSense()=='+' ){ plusCount++; }
      else { minusCount++; }
    }
  }
  for (set<ScoredSeqFlip*>::iterator it = _minusEdgeSeqs.begin(); it != _minusEdgeSeqs.end(); ++it){
    if ( (*it)->getNested()==centerContig ){
      if ( (*it)->getSense()=='+' ){ plusCount++; }
      else { minusCount++; }
    }
  }
  if ( plusCount+minusCount != 2 or (plusCount != 2 and minusCount != 2) ){
    throw AssemblyException::ArgError("AGN constructor, problem with center contig status in edge sets");
  }
}


AssemblyGraphNode::~AssemblyGraphNode(){
  for (set<ScoredSeqFlip*>::iterator it = _plusEdgeSeqs.begin(); it != _plusEdgeSeqs.end(); ++it){ delete *it; }
  for (set<ScoredSeqFlip*>::iterator it = _minusEdgeSeqs.begin(); it != _minusEdgeSeqs.end(); ++it){ delete *it; }
}


void AssemblyGraphNode::visit(){ _wasVisited = true; }
bool AssemblyGraphNode::wasVisited(){ return _wasVisited; }
void AssemblyGraphNode::unvisit(){ _wasVisited = false; }


void AssemblyGraphNode::setPlusLink(AssemblyGraphNode* plusLink, float linkScore, bool sameDir){
  _plusLink = plusLink;
  _plusLinkScore = linkScore;
  _plusLinkAgrees = sameDir;
}
void AssemblyGraphNode::setMinusLink(AssemblyGraphNode* minusLink, float linkScore, bool sameDir){
  _minusLink = minusLink;
  _minusLinkScore = linkScore;
  _minusLinkAgrees = sameDir;
}


AssemblyGraphNode* AssemblyGraphNode::getPlusLink(){ return _plusLink; }
AssemblyGraphNode* AssemblyGraphNode::getMinusLink(){ return _minusLink; }
float AssemblyGraphNode::getPlusLinkScore(){ return _plusLinkScore; }
float AssemblyGraphNode::getMinusLinkScore(){ return _minusLinkScore; }
bool AssemblyGraphNode::plusLinkAgrees(){ return _plusLinkAgrees; }
bool AssemblyGraphNode::minusLinkAgrees(){ return _minusLinkAgrees; }
void AssemblyGraphNode::setPlusAgreement(bool sameDir){ _plusLinkAgrees = sameDir; }
void AssemblyGraphNode::setMinusAgreement(bool sameDir){ _minusLinkAgrees = sameDir; }


ScoredSeq* AssemblyGraphNode::getCenterContig(){ return _centerContig; }
bool AssemblyGraphNode::centerContigIsPlus(){ return _centerIsPlus; }


void AssemblyGraphNode::flipOrientation(){
  // flip the orientations of all the sequences; avoid redundant deletes
  set<ScoredSeqFlip*> nrShallowDeletes;
  flipOrientationHelper(&_plusEdgeSeqs, &nrShallowDeletes);
  flipOrientationHelper(&_minusEdgeSeqs,&nrShallowDeletes);
  for (set<ScoredSeqFlip*>::iterator it = nrShallowDeletes.begin(); it != nrShallowDeletes.end(); ++it){ delete (*it); }
  // flip the link orientations
  AssemblyGraphNode* linkHolder = _plusLink;
  _plusLink = _minusLink;
  _minusLink = linkHolder;
  float linkScoreHolder = _plusLinkScore;
  _plusLinkScore = _minusLinkScore;
  _minusLinkScore = linkScoreHolder;
  // and the values of the linked contigs have to be inverted
  if (_plusLink != NULL){
    if (_plusLinkAgrees){ _plusLink->setPlusAgreement( false ); }
    else { _plusLink->setMinusAgreement( true ); }
  }
  if (_minusLink != NULL){
    if (_minusLinkAgrees){ _minusLink->setMinusAgreement( false ); }
    else { _minusLink->setPlusAgreement( true ); }
  }
  // now, these have to be both switched and logically inverted
  bool agreeHolder = _plusLinkAgrees;
  _plusLinkAgrees = (! _minusLinkAgrees);
  _minusLinkAgrees = (! agreeHolder);
  // record the orientation of the center contig
  _centerIsPlus = (! _centerIsPlus);
}
void AssemblyGraphNode::flipOrientationHelper(set<ScoredSeqFlip*>* setToFlip, set<ScoredSeqFlip*>* garbageBin){
  set<ScoredSeqFlip*> newSet;
  for (set<ScoredSeqFlip*>::iterator it = setToFlip->begin(); it != setToFlip->end(); ++it){
    ScoredSeq* nestedSeq = (*it)->getNested();
    switch( (*it)->getSense() ) {
    case '+': newSet.insert( ScoredSeqFlip::getFlip(nestedSeq, '-') ); break;
    case '-': newSet.insert( ScoredSeqFlip::getFlip(nestedSeq, '+') ); break;
    default: throw AssemblyException::LogicError("bad sense found in AGN::flipOrientation");
    }
  }
  garbageBin->insert(setToFlip->begin(), setToFlip->end());
  setToFlip->clear();
  setToFlip->insert(newSet.begin(), newSet.end());
}

void AssemblyGraphNode::getPlusEdgeSeqs(set<ScoredSeqFlip*>* edgeSeqs){
  for (set<ScoredSeqFlip*>::iterator it = _plusEdgeSeqs.begin(); it != _plusEdgeSeqs.end(); ++it){
    edgeSeqs->insert( ScoredSeqFlip::getFlip( (*it)->getNested(), (*it)->getSense() ) );
  }
}
void AssemblyGraphNode::getMinusEdgeSeqs(set<ScoredSeqFlip*>* edgeSeqs){
  for (set<ScoredSeqFlip*>::iterator it = _minusEdgeSeqs.begin(); it != _minusEdgeSeqs.end(); ++it){
    edgeSeqs->insert( ScoredSeqFlip::getFlip( (*it)->getNested(), (*it)->getSense() ) );
  }
}


#endif
