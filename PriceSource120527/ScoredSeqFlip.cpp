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


#ifndef SCOREDSEQFLIP_CPP
#define SCOREDSEQFLIP_CPP

#include "ScoredSeqFlip.h"


ScoredSeqFlip::~ScoredSeqFlip(){};
ScoredSeqFlip* ScoredSeqFlip::getFlip(ScoredSeq* innerSeq, char sense){
  switch (sense){
  case '+': return new ScoredSeqFlipPlus(innerSeq);
  case '-': return new ScoredSeqFlipMinus(innerSeq);
  default: throw AssemblyException::ArgError("SSFlip Factory: bad sense");
  }
}






ScoredSeqFlipPlus::ScoredSeqFlipPlus(ScoredSeq* innerSeq) : _innerSeq(innerSeq){}
ScoredSeqFlipPlus::~ScoredSeqFlipPlus(){}
float ScoredSeqFlipPlus::scoreAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->scoreAtPosPlus(position); 
  case '-': return _innerSeq->scoreAtPosMinus(position); 
  default: throw AssemblyException::ArgError("SSFlipPlus:scoreAtPosition bad sense char");
  }
}
float ScoredSeqFlipPlus::scoreAtPosPlus(long position){ return _innerSeq->scoreAtPosPlus(position); }
float ScoredSeqFlipPlus::scoreAtPosMinus(long position){ return _innerSeq->scoreAtPosMinus(position); }

float ScoredSeqFlipPlus::linkAfterPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->linkAfterPosPlus(position); 
  case '-': return _innerSeq->linkAfterPosMinus(position); 
  default: throw AssemblyException::ArgError("SSFlipPlus:linkAfterPosition bad sense char");
  }
}
float ScoredSeqFlipPlus::linkAfterPosPlus(long position){ return _innerSeq->linkAfterPosPlus(position); }
float ScoredSeqFlipPlus::linkAfterPosMinus(long position){ return _innerSeq->linkAfterPosMinus(position); }

char ScoredSeqFlipPlus::nucAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->nucAtPosPlus(position); 
  case '-': return _innerSeq->nucAtPosMinus(position); 
  default: throw AssemblyException::ArgError("SSFlipPlus:nucAtPosition bad sense char");
  }
}
char ScoredSeqFlipPlus::nucAtPosPlus(long position) { return _innerSeq->nucAtPosPlus(position); }
char ScoredSeqFlipPlus::nucAtPosMinus(long position) { return _innerSeq->nucAtPosMinus(position); }

char* ScoredSeqFlipPlus::getSeq(char sense) { return _innerSeq->getSeq(sense); }
float* ScoredSeqFlipPlus::getScores(char sense) { return _innerSeq->getScores(sense); }
float* ScoredSeqFlipPlus::getLinks(char sense) { return _innerSeq->getLinks(sense); }

char* ScoredSeqFlipPlus::getSubseq(long position, long length, char sense){ return _innerSeq->getSubseq(position,length,sense); }
float* ScoredSeqFlipPlus::getSubScores(long position, long length, char sense){ return _innerSeq->getSubScores(position,length,sense); }
float* ScoredSeqFlipPlus::getSubLinks(long position, long length, char sense){ return _innerSeq->getSubLinks(position,length,sense); }

char ScoredSeqFlipPlus::getSense() { return '+'; }
long ScoredSeqFlipPlus::size() { return _innerSeq->size(); }
void ScoredSeqFlipPlus::buffer(){}
void ScoredSeqFlipPlus::unbuffer(){}
void ScoredSeqFlipPlus::bottomBuffer(){ _innerSeq->bottomBuffer(); }
bool ScoredSeqFlipPlus::hasPairedEnd() { return false; }
ScoredSeq * ScoredSeqFlipPlus::getPairedEnd() {
  throw AssemblyException::CallingError("Paired end issues still need to be worked out for ScoredSeqFlipPlus");
}
ScoredSeq * ScoredSeqFlipPlus::getTempPairedEnd() {
  throw AssemblyException::CallingError("Paired end issues still need to be worked out for ScoredSeqFlipPlus");
}
ScoredSeq * ScoredSeqFlipPlus::shallowCopy(){ return ScoredSeq::copyShallowSeq(_innerSeq, '+'); }
ScoredSeq * ScoredSeqFlipPlus::flipCopy(){ return ScoredSeq::copyShallowSeq(_innerSeq, '-'); }
bool ScoredSeqFlipPlus::isNested(){ return true; }
ScoredSeq* ScoredSeqFlipPlus::getNested(){ return _innerSeq; }
void ScoredSeqFlipPlus::deepDelete(){
  deepUnbuffer();
  _innerSeq->deepDelete();
  delete this;
}
void ScoredSeqFlipPlus::deepUnbuffer(){
  unbuffer();
  _innerSeq->deepUnbuffer();
}






ScoredSeqFlipMinus::ScoredSeqFlipMinus(ScoredSeq* innerSeq) : _innerSeq(innerSeq){}
ScoredSeqFlipMinus::~ScoredSeqFlipMinus(){}
float ScoredSeqFlipMinus::scoreAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->scoreAtPosMinus(position);
  case '-': return _innerSeq->scoreAtPosPlus(position);
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float ScoredSeqFlipMinus::scoreAtPosPlus(long position){ return _innerSeq->scoreAtPosMinus(position); }
float ScoredSeqFlipMinus::scoreAtPosMinus(long position){ return _innerSeq->scoreAtPosPlus(position); }

float ScoredSeqFlipMinus::linkAfterPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->linkAfterPosMinus(position);
  case '-': return _innerSeq->linkAfterPosPlus(position);
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float ScoredSeqFlipMinus::linkAfterPosPlus(long position){ return _innerSeq->linkAfterPosMinus(position); }
float ScoredSeqFlipMinus::linkAfterPosMinus(long position){ return _innerSeq->linkAfterPosPlus(position); }

char ScoredSeqFlipMinus::nucAtPosition(long position, char sense) {
  switch (sense){
  case '+': return _innerSeq->nucAtPosMinus(position);
  case '-': return _innerSeq->nucAtPosPlus(position);
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
char ScoredSeqFlipMinus::nucAtPosPlus(long position) { return _innerSeq->nucAtPosMinus(position); }
char ScoredSeqFlipMinus::nucAtPosMinus(long position) { return _innerSeq->nucAtPosPlus(position); }

char* ScoredSeqFlipMinus::getSeq(char sense) {   //OK();
  switch (sense){
  case '+': return _innerSeq->getSeq('-');
  case '-': return _innerSeq->getSeq('+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float* ScoredSeqFlipMinus::getScores(char sense) {   //OK();
  switch (sense){
  case '+': return _innerSeq->getScores('-');
  case '-': return _innerSeq->getScores('+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float* ScoredSeqFlipMinus::getLinks(char sense) {   //OK();
  switch (sense){
  case '+': return _innerSeq->getLinks('-');
  case '-': return _innerSeq->getLinks('+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}


char* ScoredSeqFlipMinus::getSubseq(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubseq(position,length,'-');
  case '-': return _innerSeq->getSubseq(position,length,'+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float* ScoredSeqFlipMinus::getSubScores(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubScores(position,length,'-');
  case '-': return _innerSeq->getSubScores(position,length,'+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}
float* ScoredSeqFlipMinus::getSubLinks(long position, long length, char sense){
  switch (sense){
  case '+': return _innerSeq->getSubLinks(position,length,'-');
  case '-': return _innerSeq->getSubLinks(position,length,'+');
  default: throw AssemblyException::ArgError("SSFM: bad sense");
  }
}

char ScoredSeqFlipMinus::getSense() { return '-'; }
long ScoredSeqFlipMinus::size() { return _innerSeq->size(); }
void ScoredSeqFlipMinus::buffer(){}
void ScoredSeqFlipMinus::unbuffer(){}
void ScoredSeqFlipMinus::bottomBuffer(){ _innerSeq->bottomBuffer(); }
bool ScoredSeqFlipMinus::hasPairedEnd(){ return false; }
ScoredSeq * ScoredSeqFlipMinus::getPairedEnd() {
  throw AssemblyException::CallingError("Paired end issues still need to be worked out for ScoredSeqFlipMinus");
}
ScoredSeq * ScoredSeqFlipMinus::getTempPairedEnd() {
  throw AssemblyException::CallingError("Paired end issues still need to be worked out for ScoredSeqFlipMinus");
}
ScoredSeq * ScoredSeqFlipMinus::shallowCopy(){ return ScoredSeq::copyShallowSeq(_innerSeq, '-'); }
ScoredSeq * ScoredSeqFlipMinus::flipCopy(){ return ScoredSeq::copyShallowSeq(_innerSeq, '+'); }
bool ScoredSeqFlipMinus::isNested(){ return true; }
ScoredSeq* ScoredSeqFlipMinus::getNested(){ return _innerSeq; }
void ScoredSeqFlipMinus::deepDelete(){
  deepUnbuffer();
  _innerSeq->deepDelete();
  delete this;
}
void ScoredSeqFlipMinus::deepUnbuffer(){
  unbuffer();
  _innerSeq->deepUnbuffer();
}




#endif
