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


#ifndef SCOREDSEQWITHCARRIERS_CPP
#define SCOREDSEQWITHCARRIERS_CPP

#include "ScoredSeqWithCarriers.h"

#include <typeinfo>
using namespace::std;

long ScoredSeqWithCarriers::_readCount = 0;

ScoredSeqWithCarriers::ScoredSeqWithCarriers(){} //default constructor

ScoredSeqWithCarriers::ScoredSeqWithCarriers(ScoredSeq* innerSeq, int numOfSets) :
  _thisIsAlive( true ),
  _numOfSets(numOfSets),
  _innerSeq(innerSeq){
  _pe = 0;
  _carriedSeqSets = new set<ScoredSeq*>*[ numOfSets ];
  for (int n = 0; n < _numOfSets; ++n){ _carriedSeqSets[n] = new set<ScoredSeq*>; }
  _okToDeletePe = true;
  //_readCount++;
}


ScoredSeqWithCarriers::~ScoredSeqWithCarriers(){
  //_readCount--;
  makeDead();
  if ( _pe != 0 and _pe->isAlive() ) {
    // create a "copy" with the same references using function calls
    ScoredSeqWithCarriers * copyOfThis = new ScoredSeqWithCarriers( _innerSeq, _numOfSets );
    // delete the carriedSeqSets in the copy to make room for the re-assignment of this's carriedSeqSets
    for (long n = 0; n < _numOfSets; ++n){ delete copyOfThis->_carriedSeqSets[n]; }
    delete [] copyOfThis->_carriedSeqSets;
    copyOfThis->_carriedSeqSets = _carriedSeqSets; // copy the reference directly
    copyOfThis->addPairedEnd( _pe );
    copyOfThis->makeDead();
    _pe->replaceDeadPe( copyOfThis ); // pass that copy to the paired end
  } else {
    if ( _okToDeletePe and _pe != 0 ){ _pe->shallowDelete(); } // <- this needs to be able to actually delete without recursing
    unbuffer();
    for (long n = 0; n < _numOfSets; ++n){ delete _carriedSeqSets[n]; }
    delete [] _carriedSeqSets;
  }
}

void ScoredSeqWithCarriers::shallowDelete(){
  _okToDeletePe = false;
  delete this;
}


void ScoredSeqWithCarriers::deepDelete(){
  makeDead();
  if ( _pe != 0 and _pe->isAlive() ){
    // (deep) delete the inner seq; assume this operation is legit
    _innerSeq->deepDelete();
    // get pe's inner seq's paired end without bringing it back to life
    _innerSeq = dynamic_cast<ScoredSeqPaired*>( _pe->getNested() )->getTempPairedEnd();
  } else if ( _pe != 0 ){
    // shallow delete the paired end (no copy will be made)
    _pe->shallowDelete();
    _pe = 0;
    // this will take care of the inner seq's paired end as well
    _innerSeq->deepDelete();
    _innerSeq = 0;
  } else {
    _innerSeq->deepDelete();
    _innerSeq = 0;
  }
  delete this;
}


void ScoredSeqWithCarriers::addPairedEnd(ScoredSeqWithCarriers * ss){ _pe = ss; }
void ScoredSeqWithCarriers::makeDead(){ _thisIsAlive = false; }
void ScoredSeqWithCarriers::makeAlive(){
  _thisIsAlive = true;
  // this is ok because in order for this to be paired, this's inner
  // seq must also be paired.
  dynamic_cast<ScoredSeqPaired*>(_innerSeq)->makeAlive();
}
bool ScoredSeqWithCarriers::isAlive(){ return _thisIsAlive; }
void ScoredSeqWithCarriers::replaceDeadPe(ScoredSeqWithCarriers * deadReplacement){
  if (_pe == deadReplacement) {
    throw AssemblyException::ArgError("PE original and copy have same memory address.");
  }
  deadReplacement->makeDead();
  _pe = deadReplacement;
}


float ScoredSeqWithCarriers::scoreAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->scoreAtPosPlus(position);
  case '-': return _innerSeq->scoreAtPosMinus(position);
  default: throw AssemblyException::ArgError("SSWithCarriers:scoreAtPosition bad sense char");
  }
}
float ScoredSeqWithCarriers::scoreAtPosPlus(long position) { return _innerSeq->scoreAtPosPlus(position); }
float ScoredSeqWithCarriers::scoreAtPosMinus(long position) { return _innerSeq->scoreAtPosMinus(position); }

float ScoredSeqWithCarriers::linkAfterPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->linkAfterPosPlus(position);
  case '-': return _innerSeq->linkAfterPosMinus(position);
  default: throw AssemblyException::ArgError("SSWithCarriers:linkAfterPosition bad sense char");
  }
}
float ScoredSeqWithCarriers::linkAfterPosPlus(long position){ return _innerSeq->linkAfterPosPlus(position); }
float ScoredSeqWithCarriers::linkAfterPosMinus(long position){ return _innerSeq->linkAfterPosMinus(position); }


char ScoredSeqWithCarriers::nucAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->nucAtPosPlus(position);
  case '-': return _innerSeq->nucAtPosMinus(position);
  default: throw AssemblyException::ArgError("SSWithCarriers:nucAtPosition bad sense char");
  }
}
char ScoredSeqWithCarriers::nucAtPosPlus(long position) { return _innerSeq->nucAtPosPlus(position); }
char ScoredSeqWithCarriers::nucAtPosMinus(long position) { return _innerSeq->nucAtPosMinus(position); }

char* ScoredSeqWithCarriers::getSeq(char sense) { return _innerSeq->getSeq(sense); }
float* ScoredSeqWithCarriers::getScores(char sense){ return _innerSeq->getScores(sense); }
float* ScoredSeqWithCarriers::getLinks(char sense){ return _innerSeq->getLinks(sense); }

char* ScoredSeqWithCarriers::getSubseq(long position, long length, char sense){
  return _innerSeq->getSubseq(position, length, sense);
}
float* ScoredSeqWithCarriers::getSubScores(long position, long length, char sense){
  return _innerSeq->getSubScores(position, length, sense);
}
float* ScoredSeqWithCarriers::getSubLinks(long position, long length, char sense){
  return _innerSeq->getSubLinks(position, length, sense);
}



long ScoredSeqWithCarriers::size() {
  return _innerSeq->size();
}

void ScoredSeqWithCarriers::buffer(){
  //if (_innerSeq != 0){ _innerSeq->buffer(); }
}
void ScoredSeqWithCarriers::unbuffer(){
  //if (_innerSeq != 0){ _innerSeq->unbuffer(); }
}
void ScoredSeqWithCarriers::bottomBuffer(){
  _innerSeq->bottomBuffer();
}

bool ScoredSeqWithCarriers::hasPairedEnd() {
  return _pe != 0;
}

ScoredSeq * ScoredSeqWithCarriers::getPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  _pe->makeAlive();
  return _pe;
}
ScoredSeq * ScoredSeqWithCarriers::getTempPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}

ScoredSeq * ScoredSeqWithCarriers::shallowCopy(){
  return ScoredSeq::copyShallowSeq(_innerSeq, '+');
}
ScoredSeq * ScoredSeqWithCarriers::flipCopy(){
  return ScoredSeq::copyShallowSeq(_innerSeq, '-');
}

bool ScoredSeqWithCarriers::isNested(){ return true; }
ScoredSeq* ScoredSeqWithCarriers::getNested(){ return _innerSeq; }


void ScoredSeqWithCarriers::deepUnbuffer(){
  unbuffer();
  _innerSeq->deepUnbuffer();
}


void ScoredSeqWithCarriers::addCarriedSeq(ScoredSeq* seq, int setNum){
  _carriedSeqSets[setNum]->insert(seq);
}
void ScoredSeqWithCarriers::addCarriedSeqs(set<ScoredSeq*>* seqs, int setNum){
  _carriedSeqSets[setNum]->insert(seqs->begin(), seqs->end());
}
void ScoredSeqWithCarriers::removeCarriedSeq(ScoredSeq* seq, int setNum){
  _carriedSeqSets[setNum]->erase(seq);
}
void ScoredSeqWithCarriers::getCarriedSeqs(set<ScoredSeq*>* seqs, int setNum){
  seqs->insert(_carriedSeqSets[setNum]->begin(), _carriedSeqSets[setNum]->end());
}
void ScoredSeqWithCarriers::clearCarriedSeqs(int setNum){
  _carriedSeqSets[setNum]->clear();
}
long ScoredSeqWithCarriers::numCarriedSeqs(int setNum){
  return _carriedSeqSets[setNum]->size();
}
int ScoredSeqWithCarriers::numCarriers(){
  return _numOfSets;
}
bool ScoredSeqWithCarriers::hasCarriedSeq(ScoredSeq* seq, int setNum){
  return _carriedSeqSets[setNum]->count(seq) > 0;
}


#endif
