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


#ifndef READFILEPAIRED_CPP
#define READFILEPAIRED_CPP


#include "ReadFilePaired.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include "ScoredSeqWithCarriers.h"
#include <vector>
using namespace std;

ReadFilePaired::ReadFilePaired(){}

ReadFilePaired::ReadFilePaired(ReadFile* file, long ampliconSize) :
  _ampliconSize( ampliconSize ) {

  // check the files for agreement before starting to use them
  _readFile = file->copy();

  _getBothEnds = true;
  _isOpen = false;
}


// private constructor
ReadFilePaired::ReadFilePaired(ReadFilePaired* parent, ReadFile* rf) :
  _readFile(rf)
{
  _ampliconSize = parent->_ampliconSize;
  _getBothEnds = parent->_getBothEnds;
  _isOpen = false;
}


ReadFilePaired::~ReadFilePaired(){
  if (_isOpen){ close(); }
  delete _readFile;
}



void ReadFilePaired::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
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

ReadFile* ReadFilePaired::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip % 2 !=0 or numReadsKeep % 2 != 0){
    throw AssemblyException::LogicError("RFFPaired subfile must subdivide even numbers of reads.");
  } else if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  
  ReadFile* sub = _readFile->subFile(numReadsSkip, numReadsKeep);
  return new ReadFilePaired(this, sub);
}

ReadFilePaired* ReadFilePaired::copy(){
  return new ReadFilePaired(this, _readFile->copy());
}


string ReadFilePaired::getName(){ return _readFile->getName(); }

long ReadFilePaired::numReads(){
  long numReads = _readFile->numReads();
  if (numReads % 2 != 0){ throw AssemblyException::FileError("RFPaired has an uneven number of sequences."); }
  return numReads;
}



long ReadFilePaired::ampliconSize(){
  return _ampliconSize;
}

// open the files
void ReadFilePaired::open(){
  if (_isOpen){
    throw AssemblyException::CallingError("don't open RFFPaired, it is already open.");
  }
  _readFile->open();
  _openWithCarriers = false;
  _openNormalized = false;
  _getBothEnds = true;
  _isOpen = true;
}
void ReadFilePaired::open(bool getBothEnds){
  open();
  _getBothEnds = getBothEnds;
}

void ReadFilePaired::openWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFilePaired::openWithCarriers(int numSets, bool getBothEnds){
  openWithCarriers(numSets);
  _getBothEnds = getBothEnds;
}

void ReadFilePaired::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFilePaired::openNormalized(bool getBothEnds){
  openNormalized();
  _getBothEnds = getBothEnds;
}

void ReadFilePaired::openNormWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFilePaired::openNormWithCarriers(int numSets, bool getBothEnds){
  openNormWithCarriers(numSets);
  _getBothEnds = getBothEnds;
}


bool ReadFilePaired::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFPaired while it is closed.");
  }
  // this check ensures that only one paired-end is stored at a time
  if ( _getReadsQueue.size() > 0 ){ return true; }
  else {
    getReadPairHelper();
    return _getReadsQueue.size() > 0;
  }
}

ScoredSeq* ReadFilePaired::getRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't try to get reads from RFFPaired while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq* readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    return readToReturn;
  } else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}


void ReadFilePaired::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFPaired while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq * readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    delete readToReturn;
  } else {
    throw AssemblyException::LogicError("cannot skip read if there is no read to get.");
  }
}

void ReadFilePaired::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFPaired while it is closed.");
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
    // divide then multiply to ensure that number skipped is even
    long numToSkipPerFile = numToSkip / 2;
    _readFile->skipReads( numToSkipPerFile * 2 );
    if ( numToSkip % 2 == 1 ){ skipRead(); }
  }
}


void ReadFilePaired::close(){
  _readFile->close();
  while (_getReadsQueue.size() > 0){
    ScoredSeq * readToReturn = _getReadsQueue.front();
    _getReadsQueue.pop();
    delete readToReturn;
  }
  _getBothEnds = true;
  _isOpen = false;
}


void ReadFilePaired::getReadPairHelper(){

  if ( _readFile->hasRead() ){

    // type-independent carriers
    ScoredSeq* seqA;
    ScoredSeq* seqB;

    // create the reads themselves and pair them
    Read * newReadA = dynamic_cast<Read*>( _readFile->getRead() );
    if (! _readFile->hasRead() ){
      throw AssemblyException::FileError("RFPaired has an uneven number of sequences (found late)");
    }
    Read * newReadB = dynamic_cast<Read*>( _readFile->getRead() );
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



void ReadFilePaired::bufferQueue(Read * r){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
}
void ReadFilePaired::bufferUnqueue(Read * r){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
}


void ReadFilePaired::bufferFill(){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
  ReadFileCommunicator::bufferFill(_readFile);
}
void ReadFilePaired::bufferFill(BufferThreadedness threadedness){
  // this method is functionless since the inner read files are the ones
  // that are referenced for buffering.
  ReadFileCommunicator::bufferFill(_readFile, threadedness);
}

void ReadFilePaired::updateBuffer(){}
void ReadFilePaired::emptyBuffer(){}

#endif
