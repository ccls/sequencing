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


#ifndef SCOREDSEQSUBSEQ_CPP
#define SCOREDSEQSUBSEQ_CPP

#include "ScoredSeqSubseq.h"

long ScoredSeqSubseq::_readCount = 0;

ScoredSeqSubseq::ScoredSeqSubseq(){} //default constructor

ScoredSeqSubseq::ScoredSeqSubseq(ScoredSeq* innerSeq, long offset, long length) :
  _offset(offset),
  _size(length),
  _innerSeq(innerSeq){
  if (_offset < 0 or _size < 1 or _offset + _size > _innerSeq->size()){
    throw AssemblyException::ArgError("ScoredSeqSubseq is beyond the scope of its inner sequence");
  }
  _minusOffset = _innerSeq->size() - _size - _offset;
  //_readCount++;
}

ScoredSeqSubseq::~ScoredSeqSubseq(){
  //_readCount--;
  unbuffer();
}



float ScoredSeqSubseq::scoreAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->scoreAtPosPlus(position + _offset);
  case '-': return _innerSeq->scoreAtPosMinus(position + _minusOffset);
  default: throw AssemblyException::ArgError("SSSubseq:scoreAtPosition, sense must be (+) or (-)");
  }
}
float ScoredSeqSubseq::scoreAtPosPlus(long position){ return _innerSeq->scoreAtPosPlus(position + _offset); }
float ScoredSeqSubseq::scoreAtPosMinus(long position){ return _innerSeq->scoreAtPosMinus(position + _minusOffset); }


float ScoredSeqSubseq::linkAfterPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->linkAfterPosPlus(position + _offset);
  case '-': return _innerSeq->linkAfterPosMinus(position + _minusOffset);
  default: throw AssemblyException::ArgError("SSSubseq:linkAfterPosition, sense must be (+) or (-)");
  }
}
float ScoredSeqSubseq::linkAfterPosPlus(long position){ return _innerSeq->linkAfterPosPlus(position + _offset); }
float ScoredSeqSubseq::linkAfterPosMinus(long position){ return _innerSeq->linkAfterPosMinus(position + _minusOffset); }


char ScoredSeqSubseq::nucAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->nucAtPosPlus(position + _offset);
  case '-': return _innerSeq->nucAtPosMinus(position + _minusOffset);
  default: throw AssemblyException::ArgError("SSSubseq:nucAtPosition, sense must be (+) or (-)");
  }
}
char ScoredSeqSubseq::nucAtPosPlus(long position) { return _innerSeq->nucAtPosPlus(position + _offset); }
char ScoredSeqSubseq::nucAtPosMinus(long position) { return _innerSeq->nucAtPosMinus(position + _minusOffset); }

char* ScoredSeqSubseq::getSeq(char sense) {
  switch (sense){
  case '+': return _innerSeq->getSubseq(_offset, _size, sense);
  case '-': return _innerSeq->getSubseq(_minusOffset, _size, sense);
  default: throw AssemblyException::ArgError("SSSubseq::getSeq, sense must be (+) or (-)");
  }
}
float* ScoredSeqSubseq::getScores(char sense) {
  //float* scores = new float[_size + 1];
  switch (sense){
  case '+': return _innerSeq->getSubScores(_offset, _size, sense);
    //for (long n = 0; n < _size; ++n){ scores[n] = _innerSeq->scoreAtPosPlus(n + _offset); }
    //return scores;
  case '-': return _innerSeq->getSubScores(_minusOffset, _size, sense);
    //for (long n = 0; n < _size; ++n){ scores[n] = _innerSeq->scoreAtPosMinus(n + _minusOffset); }
    //return scores;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
float* ScoredSeqSubseq::getLinks(char sense) {
  switch (sense){
  case '+': return _innerSeq->getSubLinks(_offset, _size - 1, sense);
  case '-': return _innerSeq->getSubLinks(_minusOffset, _size - 1, sense);
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}


char* ScoredSeqSubseq::getSubseq(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubseq(position + _offset, length, sense);
  case '-': return _innerSeq->getSubseq(position + _minusOffset, length, sense);
  default: throw AssemblyException::ArgError("SSSubseq::getSubseq, sense must be (+) or (-)");
  }
}
float* ScoredSeqSubseq::getSubScores(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubScores(position + _offset, length, sense);
  case '-': return _innerSeq->getSubScores(position + _minusOffset, length, sense);
  default: throw AssemblyException::ArgError("SSSubseq::getSubScores, sense must be (+) or (-)");
  }
}
float* ScoredSeqSubseq::getSubLinks(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubLinks(position + _offset, length, sense);
  case '-': return _innerSeq->getSubLinks(position + _minusOffset, length, sense);
  default: throw AssemblyException::ArgError("SSSubseq::getSubLinks, sense must be (+) or (-)");
  }
}



long ScoredSeqSubseq::size() {
  return _size;
}

void ScoredSeqSubseq::buffer(){}
void ScoredSeqSubseq::unbuffer(){}
void ScoredSeqSubseq::bottomBuffer(){
  _innerSeq->bottomBuffer();
}


bool ScoredSeqSubseq::hasPairedEnd() {
  return false;
}

ScoredSeq * ScoredSeqSubseq::getPairedEnd() {
  throw AssemblyException::CallingError("subseq has no paired end");
}
ScoredSeq * ScoredSeqSubseq::getTempPairedEnd() {
  throw AssemblyException::CallingError("subseq has no paired end");
}

ScoredSeq * ScoredSeqSubseq::shallowCopy(){
  return ScoredSeq::copyShallowSeq(this, '+');
}
ScoredSeq * ScoredSeqSubseq::flipCopy(){
  return ScoredSeq::copyShallowSeq(this, '-');
}


bool ScoredSeqSubseq::isNested(){ return true; }
ScoredSeq* ScoredSeqSubseq::getNested(){ return _innerSeq; }
void ScoredSeqSubseq::deepDelete(){
  _innerSeq->deepDelete();
  delete this;
}
void ScoredSeqSubseq::deepUnbuffer(){
  unbuffer();
  _innerSeq->deepUnbuffer();
}


#endif
