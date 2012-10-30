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


#ifndef READ_CPP
#define READ_CPP


#include "Read.h"

#include "ScoredSeq.h"
#include "ReadFileIndex.h"
#include "ReadFile.h"
#include "AssemblyException.h"

#include <omp.h>

# include <iostream>
using namespace::std;

long Read::_readCount = 0;

Read::Read(){}; //default constructor
Read::Read(ReadFile *readFile, ReadFileIndex * rfi) :
  _requestedBuffer( false ),
  _readFile( readFile ),
  _rfi( rfi ),
  _pe( 0 ),
  _buffer( 0 ),
  _thisIsAlive( true )
{
  _okToDeletePe = true;
  //_readCount++;
}

Read::Read(const Read& r){
  throw AssemblyException::CopyConstructorError("Read");
}
Read& Read::operator=(const Read& that){
  throw AssemblyException::CopyAssignmentOperatorError("Read");
}

Read::~Read(){
  //_readCount--;
  makeDead();
  if ( _pe != 0 and _pe->isAlive() ) {
    unbufferPrivate();
    // create a "copy" with the same references using function calls
    Read * copyOfThis = new Read( _readFile, _rfi );
    copyOfThis->addPairedEnd( _pe );
    // no need for this to be added critically since nobody else has access
    // to this new copy yet.
    if ( _buffer != 0 ){ copyOfThis->addScoredSeqImp( _buffer ); }
    copyOfThis->makeDead();
    _pe->replaceDeadPe( copyOfThis ); // pass that copy to the paired end
    // this does not need to be critical because nobody should have
    // legit access to this object as it is being deleted.
    unqueuePrivate();
  } else {
    if ( _okToDeletePe and _pe != 0 ){ _pe->shallowDelete(); } // <- this needs to be able to actually delete without recursing
    // tell _readFile that this member is dead, then find out if it is ok to delete _readFile

    // this was not thread-safe
    //ReadFileCommunicator::readIsBeingDeleted( _readFile, this );
    // delete the non-singleton references
    // this does not need to be critical because nobody should have
    // legit access to this object as it is being deleted.
    unbufferPrivate();
    delete _rfi;
  }
}
void Read::shallowDelete(){
  _okToDeletePe = false;
  delete this;
}

// BEGIN NEW VERSION

float Read::scoreAtPosition(long position, char sense){
  float score;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    switch (sense){
    case '+': score = _buffer->scoreAtPosPlus(position); break;
    case '-': score = _buffer->scoreAtPosMinus(position); break;
    default: throw AssemblyException::ArgError("Read:scoreAtPosition bad sense char");
    }
  }
  return score;
}
float Read::scoreAtPosPlus(long position){
  float score;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    score = _buffer->scoreAtPosPlus(position);
  }
  return score;
}
float Read::scoreAtPosMinus(long position){
  float score;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    score = _buffer->scoreAtPosMinus(position);
  }
  return score;
}


float Read::linkAfterPosition(long position, char sense){
  float link;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    switch (sense){
    case '+': link =  _buffer->linkAfterPosPlus(position); break;
    case '-': link =  _buffer->linkAfterPosMinus(position); break;
    default: throw AssemblyException::ArgError("Read:linkAfterPosition bad sense char");
    }
  }
  return link;
}
float Read::linkAfterPosPlus(long position){
  float link;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    link =  _buffer->linkAfterPosPlus(position);
  }
  return link;
}
float Read::linkAfterPosMinus(long position){
  float link;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    link =  _buffer->linkAfterPosMinus(position);
  }
  return link;
}


char Read::nucAtPosition(long position, char sense){
  char nuc;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    switch (sense){
    case '+': nuc = _buffer->nucAtPosPlus(position); break;
    case '-': nuc = _buffer->nucAtPosMinus(position); break;
    default: throw AssemblyException::ArgError("Read:nucAtPosition bad sense char");
    }
  }
  return nuc;
}
char Read::nucAtPosPlus(long position) {
  char nuc;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    nuc = _buffer->nucAtPosPlus(position);
  }
  return nuc;
}
char Read::nucAtPosMinus(long position) {
  char nuc;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    nuc = _buffer->nucAtPosMinus(position);
  }
  return nuc;
}

char* Read::getSeq(char sense){
  char* seq;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    seq = _buffer->getSeq(sense);
  }
  return seq;
}
float* Read::getScores(char sense){
  float* scores;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    scores = _buffer->getScores(sense);
  }
  return scores;
}
float* Read::getLinks(char sense){
  float* links;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    links = _buffer->getLinks(sense);
  }
  return links;
}



char* Read::getSubseq(long position, long length, char sense){
  char* seq;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    seq = _buffer->getSubseq(position,length,sense);
  }
  return seq;
}
float* Read::getSubScores(long position, long length, char sense){
  float* scores;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    scores = _buffer->getSubScores(position,length,sense);
  }
  return scores;
}
float* Read::getSubLinks(long position, long length, char sense){
  float* links;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    links = _buffer->getSubLinks(position,length,sense);
  }
  return links;
}



long Read::size(){
  long size;
  #pragma omp critical (Read)
  {
    fillBufferPrivate();
    size = _buffer->size();
  }
  return size;
}


void Read::buffer(){
  #pragma omp critical (Read)
  {
    bufferPrivate();
  }
}
void Read::unbuffer(){
  #pragma omp critical (Read)
  {
    unbufferPrivate();
  }
}
void Read::bufferPrivate(){
  if ( _buffer == 0 and !_requestedBuffer ){
    ReadFileCommunicator::bufferQueue(_readFile,this);
    _requestedBuffer = true;
  }
}
void Read::unbufferPrivate(){
  delete _buffer;
  _buffer = 0;
  unqueuePrivate();
}

void Read::bottomBuffer(){ buffer(); }
void Read::deepDelete(){ delete this; }
void Read::deepUnbuffer(){ unbuffer(); }

void Read::fillBufferPrivate(){
  if ( _buffer == 0 ){
    bufferPrivate();
    ReadFileCommunicator::bufferFill(_readFile);
    _requestedBuffer = false;
  }
}
void Read::unqueuePrivate(){
  if (_requestedBuffer){
    ReadFileCommunicator::bufferUnqueue(_readFile,this);
    _requestedBuffer = false;
  }
}

ScoredSeq * Read::shallowCopy(){
  return ScoredSeq::copyShallowSeq(this, '+');
}

ScoredSeq * Read::flipCopy(){
  return ScoredSeq::copyShallowSeq(this, '-');
}



void Read::addScoredSeqImp(ScoredSeq* ssi){
  if (_buffer == 0){ _buffer = ssi; }
  else { delete ssi; }
  _requestedBuffer = false;
}

ReadFileIndex * Read::getRfiRef(){ return _rfi; }
ReadFile * Read::getReadFile(){ return _readFile; }


bool Read::hasPairedEnd() { return _pe != 0; }

ScoredSeq * Read::getPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  _pe->makeAlive();
  return _pe;
}
ScoredSeq * Read::getTempPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::LogicError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}


void Read::makeDead(){ _thisIsAlive = false; }
void Read::makeAlive(){ _thisIsAlive = true; }
bool Read::isAlive(){ return _thisIsAlive; }


void Read::replaceDeadPe(Read * deadReplacement){
  if (_pe == deadReplacement) {
    throw AssemblyException::ArgError("PE original and copy have same memory address.");
  }
  deadReplacement->makeDead();
  _pe = deadReplacement;
}


void Read::addPairedEnd(Read * ss){
  _pe = ss;
}

bool Read::isNested(){ return false; }



#endif
