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


#ifndef READFILEDUAL_CPP
#define READFILEDUAL_CPP


#include "ReadFileDual.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include "ScoredSeqWithCarriers.h"
#include <vector>
using namespace std;

ReadFileDual::ReadFileDual(){}

ReadFileDual::ReadFileDual(ReadFile* fileA, ReadFile* fileB, long ampliconSize) :
  _ampliconSize( ampliconSize ) {

  // check the files for agreement before starting to use them
  _readFileA = fileA->copy();
  _readFileB = fileB->copy();

  _getBothEnds = true;
  _isOpen = false;
}


// private constructor
ReadFileDual::ReadFileDual(ReadFileDual* parent, ReadFile* rfA, ReadFile* rfB) :
  _readFileA(rfA),
  _readFileB(rfB)
{
  _ampliconSize = parent->_ampliconSize;
  _getBothEnds = parent->_getBothEnds;
  _isOpen = false;
}


ReadFileDual::~ReadFileDual(){
  if (_isOpen){ close(); }
  delete _readFileA;
  delete _readFileB;
}



void ReadFileDual::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
  // determine the number of reads to include in each file
  long currentTally = 0;
  long totalHalfReads = numReads() / 2;
  for (int fileNum = 0; fileNum < numFiles; ++fileNum){
    long priorTally = currentTally;
    currentTally = (fileNum + 1) * totalHalfReads / numFiles;
    long halfToAdd = currentTally - priorTally;
    if (halfToAdd > 0){
      // create the file
      putFilesHere->insert( subFile(priorTally * 2, halfToAdd * 2) );
    }
  }
}

ReadFile* ReadFileDual::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip % 2 != 0 or numReadsKeep % 2 != 0){
    throw AssemblyException::LogicError("RFFDual subfile must subdivide even numbers of reads.");
  } else if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  long halfSkip = numReadsSkip / 2;
  long halfKeep = numReadsKeep / 2;
  ReadFile* subA = _readFileA->subFile(halfSkip, halfKeep);
  ReadFile* subB = _readFileB->subFile(halfSkip, halfKeep);
  return new ReadFileDual(this, subA, subB);
}

ReadFileDual* ReadFileDual::copy(){
  return new ReadFileDual(this, _readFileA->copy(), _readFileB->copy());
}


string ReadFileDual::getName(){ return _readFileA->getName() + "::" + _readFileB->getName(); }

long ReadFileDual::numReads(){
  if (_readFileA->numReads() != _readFileB->numReads()){
    cout << _readFileA->numReads() << endl;
    cout << _readFileB->numReads() << endl;
    throw AssemblyException::LogicError("RFFDual's two component read files must have the same num entries.");
  }
  return _readFileA->numReads() * 2;
}



long ReadFileDual::ampliconSize(){
  return _ampliconSize;
}

// open the files
void ReadFileDual::open(){
  if (_isOpen){
    throw AssemblyException::CallingError("don't open RFFDual, it is already open.");
  }
  _readFileA->open();
  _readFileB->open();
  _openWithCarriers = false;
  _openNormalized = false;
  _getBothEnds = true;
  _isOpen = true;
}
void ReadFileDual::open(bool getBothEnds){
  open();
  _getBothEnds = getBothEnds;
}

void ReadFileDual::openWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFileDual::openWithCarriers(int numSets, bool getBothEnds){
  openWithCarriers(numSets);
  _getBothEnds = getBothEnds;
}

void ReadFileDual::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFileDual::openNormalized(bool getBothEnds){
  openNormalized();
  _getBothEnds = getBothEnds;
}

void ReadFileDual::openNormWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFileDual::openNormWithCarriers(int numSets, bool getBothEnds){
  openNormWithCarriers(numSets);
  _getBothEnds = getBothEnds;
}


bool ReadFileDual::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFDual while it is closed.");
  }
  // this check ensures that only one paired-end is stored at a time
  if ( _getReadsQueue.size() > 0 ){ return true; }
  else {
    getReadPairHelper();
    return _getReadsQueue.size() > 0;
  }
}

ScoredSeq* ReadFileDual::getRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't try to get reads from RFFDual while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq* readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    return readToReturn;
  } else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}


void ReadFileDual::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFDual while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq * readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    delete readToReturn;
  } else {
    throw AssemblyException::LogicError("cannot skip read if there is no read to get.");
  }
}

void ReadFileDual::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFDual while it is closed.");
  }
  if (numToSkip < 1){}
  if (numToSkip == 1){ skipRead(); }
  else {
    while ( _getReadsQueue.size() > 0 and numToSkip > 0 ){
      ScoredSeq * readToReturn = _getReadsQueue.front();
      _getReadsQueue.pop();
      delete readToReturn;
      numToSkip--;
    }
    long numToSkipPerFile = numToSkip / 2;
    _readFileA->skipReads( numToSkipPerFile );
    _readFileB->skipReads( numToSkipPerFile );
    if ( numToSkip % 2 == 1 ){ skipRead(); }
  }
}


void ReadFileDual::close(){
  _readFileA->close();
  _readFileB->close();
  while (_getReadsQueue.size() > 0){
    ScoredSeq * readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    delete readToReturn;
  }
  _getBothEnds = true;
  _isOpen = false;
}


void ReadFileDual::getReadPairHelper(){

  if ( _readFileA->hasRead() and _readFileB->hasRead() ){

    // type-independent carriers
    ScoredSeq* seqA;
    ScoredSeq* seqB;

    // create the reads themselves and pair them
    Read * newReadA = dynamic_cast<Read*>( _readFileA->getRead() );
    Read * newReadB = dynamic_cast<Read*>( _readFileB->getRead() );
    seqA = newReadA;
    seqB = newReadB;
    ReadFileCommunicator::pairEnds(newReadA, newReadB);

    if (_openNormalized){
      ScoredSeqNormalized* ssnA = new ScoredSeqNormalized( seqA );
      ScoredSeqNormalized* ssnB = new ScoredSeqNormalized( seqB );
      seqA = ssnA;
      seqB = ssnB;
      ReadFileCommunicator::pairEnds(ssnA, ssnB);
    }
    if (_openWithCarriers){
      ScoredSeqWithCarriers* sswcA = new ScoredSeqWithCarriers( seqA, _numCarrierSets );
      ScoredSeqWithCarriers* sswcB = new ScoredSeqWithCarriers( seqB, _numCarrierSets );
      seqA = sswcA;
      seqB = sswcB;
      ReadFileCommunicator::pairEnds(sswcA, sswcB);
    }

    // pair the reads and add them to the queue
    _getReadsQueue.push( seqA );
    // if the pair will not be returned, make it dead
    if (_getBothEnds){ _getReadsQueue.push( seqB ); }
    else { seqB->deepDelete(); }

  //} else {
    //_getReadsFile.close();
  }
}



void ReadFileDual::bufferQueue(Read * r){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
}
void ReadFileDual::bufferUnqueue(Read * r){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
}


void ReadFileDual::bufferFill(){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
  ReadFileCommunicator::bufferFill(_readFileA);
  ReadFileCommunicator::bufferFill(_readFileB);
}
void ReadFileDual::bufferFill(BufferThreadedness threadedness){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
  ReadFileCommunicator::bufferFill(_readFileA, threadedness);
  ReadFileCommunicator::bufferFill(_readFileB, threadedness);
}

void ReadFileDual::updateBuffer(){}
void ReadFileDual::emptyBuffer(){}

#endif
