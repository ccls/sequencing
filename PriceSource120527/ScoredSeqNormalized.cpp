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


#ifndef SCOREDSEQNORMALIZED_CPP
#define SCOREDSEQNORMALIZED_CPP

#include "ScoredSeqNormalized.h"

long ScoredSeqNormalized::_readCount = 0;

ScoredSeqNormalized::ScoredSeqNormalized(){}; //default constructor

ScoredSeqNormalized::ScoredSeqNormalized(ScoredSeq* innerSeq) :
  _normCount(0),
  _innerSeq(innerSeq){
  _bufferSeq = 0;
  _thisIsAlive = true;
  _pe = 0;
  _okToDeletePe = true;
  //_readCount++;
  //OK();
}

ScoredSeqNormalized::~ScoredSeqNormalized(){
  //_readCount--;
  makeDead();
  if ( _pe != 0 and _pe->isAlive() ) {
    unbuffer();
    // create a "copy" with the same references using function calls
    ScoredSeqNormalized * copyOfThis = new ScoredSeqNormalized( _innerSeq );
    copyOfThis->addCounts( _normCount );
    copyOfThis->addPairedEnd( _pe );
    copyOfThis->makeDead();
    _pe->replaceDeadPe( copyOfThis ); // pass that copy to the paired end
  } else {
    if ( _okToDeletePe and _pe != 0 ){ _pe->shallowDelete(); } // <- this needs to be able to actually delete without recursing
    unbuffer();
  }
}

void ScoredSeqNormalized::shallowDelete(){
  _okToDeletePe = false;
  delete this;
}


void ScoredSeqNormalized::deepDelete(){
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

void ScoredSeqNormalized::OK(){
  if (_innerSeq == 0){
    throw AssemblyException::LogicError("inner seq should not be null in SSNormalized");
  }
  _innerSeq->getSeq('+'); // call a method
}


void ScoredSeqNormalized::addPairedEnd(ScoredSeqNormalized * ss){ _pe = ss; }
void ScoredSeqNormalized::makeDead(){ _thisIsAlive = false; }
void ScoredSeqNormalized::makeAlive(){
  // this is ok because in order for this to be paired, this's inner
  // seq must also be paired.
  _thisIsAlive = true;
  dynamic_cast<ScoredSeqPaired*>(_innerSeq)->makeAlive();
  //OK();
}
bool ScoredSeqNormalized::isAlive(){ return _thisIsAlive; }
void ScoredSeqNormalized::replaceDeadPe(ScoredSeqNormalized * deadReplacement){
  if (_pe == deadReplacement) {
    throw AssemblyException::ArgError("PE original and copy have same memory address.");
  }
  deadReplacement->makeDead();
  _pe = deadReplacement;
  //OK();
}



float ScoredSeqNormalized::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+':
    if (_normCount == 0){ return _innerSeq->scoreAtPosPlus(position); }
    else { return _innerSeq->scoreAtPosPlus(position) / _normCount; }
  case '-':
    if (_normCount == 0){ return _innerSeq->scoreAtPosMinus(position); }
    else { return _innerSeq->scoreAtPosMinus(position) / _normCount; }
  default: throw AssemblyException::ArgError("SSNormalized:scoreAtPosition bad sense char");
  }
}
float ScoredSeqNormalized::scoreAtPosPlus(long position) {
  //OK();
  if (_normCount == 0){ return _innerSeq->scoreAtPosPlus(position); }
  else { return _innerSeq->scoreAtPosPlus(position) / _normCount; }
}
float ScoredSeqNormalized::scoreAtPosMinus(long position) {
  //OK();
  if (_normCount == 0){ return _innerSeq->scoreAtPosMinus(position); }
  else { return _innerSeq->scoreAtPosMinus(position) / _normCount; }
}



float ScoredSeqNormalized::linkAfterPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+':
    if (_normCount == 0){ return _innerSeq->linkAfterPosPlus(position); }
    else { return _innerSeq->linkAfterPosPlus(position) / _normCount; }
  case '-':
    if (_normCount == 0){ return _innerSeq->linkAfterPosMinus(position); }
    else { return _innerSeq->linkAfterPosMinus(position) / _normCount; }
  default: throw AssemblyException::ArgError("SSNormalized:scoreAtPosition bad sense char");
  }
}
float ScoredSeqNormalized::linkAfterPosPlus(long position) {
  //OK();
  if (_normCount == 0){ return _innerSeq->linkAfterPosPlus(position); }
  else { return _innerSeq->linkAfterPosPlus(position) / _normCount; }
}
float ScoredSeqNormalized::linkAfterPosMinus(long position) {
  //OK();
  if (_normCount == 0){ return _innerSeq->linkAfterPosMinus(position); }
  else { return _innerSeq->linkAfterPosMinus(position) / _normCount; }
}


char ScoredSeqNormalized::nucAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->nucAtPosPlus(position);
  case '-': return _innerSeq->nucAtPosMinus(position);
  default: throw AssemblyException::ArgError("SSNormalized:nucAtPosition bad sense char");
  }
}
char ScoredSeqNormalized::nucAtPosPlus(long position) { return _innerSeq->nucAtPosPlus(position); }
char ScoredSeqNormalized::nucAtPosMinus(long position) { return _innerSeq->nucAtPosMinus(position); }

char* ScoredSeqNormalized::getSeq(char sense) {
  //OK();
  return _innerSeq->getSeq(sense);
}
float* ScoredSeqNormalized::getScores(char sense) {
  float* scores = _innerSeq->getScores(sense);
  if (_normCount > 1){
    long size = _innerSeq->size();
    for (long n = 0; n < size; ++n){ scores[n] /= _normCount; }
  }
  return scores;
}
float* ScoredSeqNormalized::getLinks(char sense) {
  float* links = _innerSeq->getLinks(sense);
  if (_normCount > 1){
    long sizeM1 = _innerSeq->size() - 1;
    for (long n = 0; n < sizeM1; ++n){ links[n] /= _normCount; }
  }
  return links;
}



char* ScoredSeqNormalized::getSubseq(long position, long length, char sense){
  //OK();
  return _innerSeq->getSubseq(position, length, sense);
}
float* ScoredSeqNormalized::getSubScores(long position, long length, char sense){
  float* scores = _innerSeq->getSubScores(position, length, sense);
  if (_normCount > 1){
    for (long n = 0; n < length; ++n){ scores[n] /= _normCount; }
  }
  return scores;
}
float* ScoredSeqNormalized::getSubLinks(long position, long length, char sense){
  float* links = _innerSeq->getSubLinks(position, length, sense);
  if (_normCount > 1){
    for (long n = 0; n < length; ++n){ links[n] /= _normCount; }
  }
  return links;
}


long ScoredSeqNormalized::size() {
  //OK();
  return _innerSeq->size();
}

void ScoredSeqNormalized::buffer(){}
void ScoredSeqNormalized::unbuffer(){}
void ScoredSeqNormalized::bottomBuffer(){
  //OK();
  _innerSeq->bottomBuffer();
}


bool ScoredSeqNormalized::hasPairedEnd() {
  //OK();
  return _pe != 0;
}

ScoredSeq * ScoredSeqNormalized::getPairedEnd() {
  //OK();
  _pe->makeAlive();
  return _pe;
}
ScoredSeq * ScoredSeqNormalized::getTempPairedEnd() {
  //OK();
  return _pe;
}

ScoredSeq * ScoredSeqNormalized::shallowCopy(){
  //OK();
  return ScoredSeq::copyShallowSeq(this, '+');
}
ScoredSeq * ScoredSeqNormalized::flipCopy(){
  //OK();
  return ScoredSeq::copyShallowSeq(this, '-');
}

// these necessitate a re-buffer
void ScoredSeqNormalized::addCount(){
  //OK();
  _normCount++;
  if (_bufferSeq != 0){
    unbuffer();
    buffer();
  }
}
void ScoredSeqNormalized::addCounts(long counts){
  //OK();
  _normCount += counts;
}
long ScoredSeqNormalized::getCount(){ return _normCount; }
void ScoredSeqNormalized::resetCount(){
  _normCount = 0;
}

bool ScoredSeqNormalized::isNested(){ return true; }
ScoredSeq* ScoredSeqNormalized::getNested(){ return _innerSeq; }

void ScoredSeqNormalized::deepUnbuffer(){
  unbuffer();
  _innerSeq->deepUnbuffer();
}


#endif
