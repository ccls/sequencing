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


#ifndef SCOREDSEQMONOSCORE_CPP
#define SCOREDSEQMONOSCORE_CPP

#include "ScoredSeqMonoScore.h"

#include "ScoredSeq.h"
#include "AssemblyException.h"
#include <cstring>
#include <omp.h>
using namespace::std;


void ScoredSeqMonoScore::constructorSeqHelper(string seq){
  _size = seq.size();
  _seq = new char [_size+1];
  strcpy (_seq, seq.c_str());
  _rcSeq = reverseComplement(_seq,_size);
}

// this one doesn't use the constructor helper
ScoredSeqMonoScore::ScoredSeqMonoScore(ScoredSeq* seq, float score) :
  _score(score){
  _size = seq->size();
  _seq = seq->getSeq('+');
  _rcSeq = seq->getSeq('-');
  _pe = 0;
  _thisIsAlive = true;
  //OK();
}

ScoredSeqMonoScore::ScoredSeqMonoScore(string seq, float score) :
  _score(score) {
  constructorSeqHelper(seq);
  _pe = 0;
  _thisIsAlive = true;
  //OK();
}

// CHAR ARRAYS
ScoredSeqMonoScore::ScoredSeqMonoScore(char* seq, float score, long size, bool expRep) :
  _size(size),
  _score(score) {
  if (expRep){ _seq = seq; }
  else { constructorSeqCopyHelper(seq); }
  _rcSeq = reverseComplement(_seq,_size);
  _pe = 0;
  _thisIsAlive = true;
  //OK();
}

// helper method - _size is already specified
void ScoredSeqMonoScore::constructorSeqCopyHelper(char* seq){
  _seq = new char[_size+1];
  _seq[_size] = '\0';
  for (long n = 0; n < _size; ++n){ _seq[n] = seq[n]; }
}

ScoredSeqMonoScore::ScoredSeqMonoScore(ScoredSeqMonoScore* seq, char sense) :
  _size(seq->_size),
  _score(seq->_score),
  _pe(NULL),
  _thisIsAlive(true)
{
  switch(sense){
  case '+':
    _seq = seq->getSeq('+');
    _rcSeq = seq->getSeq('-');
    break;
  case '-':
    _seq = seq->getSeq('-');
    _rcSeq = seq->getSeq('+');
    break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}

ScoredSeqMonoScore::~ScoredSeqMonoScore(){
  bool deleteSelf = true; // keeps track if the reference fields need to be destroyed
  if (_thisIsAlive){
    makeDead();
    if ( _pe != 0 ){
      if ( _pe->isAlive() ) {
	// create a "copy" with the same references using function calls
	ScoredSeqMonoScore * copyOfThis = new ScoredSeqMonoScore( _seq, _score, _size, true );
	copyOfThis->addPairedEnd( _pe );
	copyOfThis->makeDead();
	_pe->replaceDeadPe( copyOfThis ); // pass that copy to the paired end
	deleteSelf = false;
      } else {
	delete _pe;
	_pe = 0;
	delete [] _seq;
      }
    } else {
      delete [] _seq;
    }
  } else {
    delete [] _seq;
  }
  delete [] _rcSeq;
}

void ScoredSeqMonoScore::OK(){
}

float ScoredSeqMonoScore::scoreAtPosition(long position, char sense) {
  //OK();
  if ( nucAtPosition(position,sense) == 'N' ){ return 0; }
  else{ return _score; }
}
float ScoredSeqMonoScore::scoreAtPosPlus(long position) {
  //OK();
  if ( _seq[position] == 'N' ){ return 0; }
  else{ return _score; }
}
float ScoredSeqMonoScore::scoreAtPosMinus(long position) {
  //OK();
  if ( _rcSeq[position] == 'N' ){ return 0; }
  else{ return _score; }
}


float ScoredSeqMonoScore::linkAfterPosition(long position, char sense){ return _score; }
float ScoredSeqMonoScore::linkAfterPosPlus(long position){ return _score; }
float ScoredSeqMonoScore::linkAfterPosMinus(long position){ return _score; }


char ScoredSeqMonoScore::nucAtPosition(long position, char sense) {
  //OK();
  //if (position < 0 or position >= _size){ throw AssemblyException::LengthError("position is outside of the scope of the sequence"); }
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char ScoredSeqMonoScore::nucAtPosPlus(long position) { return _seq[ position ]; }
char ScoredSeqMonoScore::nucAtPosMinus(long position) { return _rcSeq[ position ]; }

char* ScoredSeqMonoScore::getSeq(char sense) {
  //OK();
  char* seq = new char[_size + 1];
  switch (sense){
  case '+':
    strcpy(seq, _seq);
    return seq;
  case '-':
    strcpy(seq, _rcSeq);
    return seq;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
float* ScoredSeqMonoScore::getScores(char sense) {
  //OK();
  float* scores = new float[_size + 1];
  switch (sense){
  case '+':
    for (long n = 0; n < _size; ++n){
      if (_seq[n] == 'N'){ scores[n] = 0; }
      else { scores[n] = _score; }
    }
    break;
  case '-':
    for (long n = 0; n < _size; ++n){
      if (_rcSeq[n] == 'N'){ scores[n] = 0; }
      else { scores[n] = _score; }
    }
    break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
  return scores;
}
float* ScoredSeqMonoScore::getLinks(char sense) {
  //OK();
  float* links = new float[_size + 1];
  long sizeM1 = _size - 1;
  for (long n = 0; n < sizeM1; ++n){ links[n] = _score; }
  return links;
}




char* ScoredSeqMonoScore::getSubseq(long position, long length, char sense){
  //OK();
  char* subseq = new char[ length + 1 ];
  subseq[length] = '\0';
  switch (sense){
  case '+':
    for (long n = 0; n < length; ++n){ subseq[n] = _seq[position + n]; }
    return subseq;
  case '-':
    for (long n = 0; n < length; ++n){ subseq[n] = _rcSeq[position + n]; }
    return subseq;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
float* ScoredSeqMonoScore::getSubScores(long position, long length, char sense){
  //OK();
  float* scores = new float[length + 1];
  switch (sense){
  case '+':
    for (long n = 0; n < length; ++n){
      if (_seq[position + n] == 'N'){ scores[n] = 0; }
      else { scores[n] = _score; }
    }
    break;
  case '-':
    for (long n = 0; n < length; ++n){
      if (_rcSeq[position + n] == 'N'){ scores[n] = 0; }
      else { scores[n] = _score; }
    }
    break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
  return scores;
}
float* ScoredSeqMonoScore::getSubLinks(long position, long length, char sense){
  //OK();
  float* links = new float[length + 1];
  for (long n = 0; n < length; ++n){ links[n] = _score; }
  return links;
}




long ScoredSeqMonoScore::size() {
  //OK();
  return _size;
}

void ScoredSeqMonoScore::buffer(){}
void ScoredSeqMonoScore::unbuffer(){}
void ScoredSeqMonoScore::bottomBuffer(){}
void ScoredSeqMonoScore::deepDelete(){ delete this; }
void ScoredSeqMonoScore::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it


bool ScoredSeqMonoScore::hasPairedEnd() { return _pe != 0; }
ScoredSeq * ScoredSeqMonoScore::getPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  _pe->makeAlive();
  return _pe;
}
ScoredSeq * ScoredSeqMonoScore::getTempPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}

void ScoredSeqMonoScore::makeDead(){ _thisIsAlive = false; }
void ScoredSeqMonoScore::makeAlive(){ _thisIsAlive = true; }
bool ScoredSeqMonoScore::isAlive(){ return _thisIsAlive; }

void ScoredSeqMonoScore::replaceDeadPe(ScoredSeqMonoScore * deadReplacement){
  if (_pe == deadReplacement) {
    throw AssemblyException::ArgError("PE original and copy have same memory address.");
  }
  deadReplacement->makeDead();
  _pe = deadReplacement;
}


void ScoredSeqMonoScore::addPairedEnd(ScoredSeqMonoScore * ss){
  _pe = ss;
}

ScoredSeq * ScoredSeqMonoScore::shallowCopy(){
  return new ScoredSeqMonoScore(this,'+');
}
ScoredSeq * ScoredSeqMonoScore::flipCopy(){
  return new ScoredSeqMonoScore(this,'-');
}


bool ScoredSeqMonoScore::isNested(){ return false; }



#endif
