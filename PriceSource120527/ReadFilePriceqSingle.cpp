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


#ifndef READFILEPRICEQSINGLE_CPP
#define READFILEPRICEQSINGLE_CPP


#include "ReadFilePriceqSingle.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <algorithm>
#include "ScoredSeqWithCarriers.h"
#include "OutputFilePriceq.h"
#include <vector>
#include <limits.h>
#include <cstring>
using namespace::std;

long ReadFilePriceqSingle::_maxBlockSize = LONG_MAX;


ReadFilePriceqSingle::ReadFilePriceqSingle(string filename) :
  _readStart( 0 ),
  _readLength( 0 ),
  _invert(false) {
  constructorHelper(filename);
}
ReadFilePriceqSingle::ReadFilePriceqSingle(bool invert, string filename) :
  _readStart( 0 ),
  _readLength( 0 ),
  _invert(invert) {
  constructorHelper(filename);
}

ReadFilePriceqSingle::ReadFilePriceqSingle(string filename, int readStart, int readLength) :
  _readStart( readStart ),
  _readLength( readLength ),
  _invert(false) {
  constructorHelper(filename);
}
ReadFilePriceqSingle::ReadFilePriceqSingle(bool invert, string filename, int readStart, int readLength) :
  _readStart( readStart ),
  _readLength( readLength ),
  _invert(invert) {
  constructorHelper(filename);
}


void ReadFilePriceqSingle::constructorHelper(string filename){
  _okToSkip.insert('\r');
  _okToSkip.insert('\n');
  _okToSkip.insert('\0');

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::ArgError("RFPQ priceq file does not exist."); }

  if ( _readLength < 0 ){
    throw AssemblyException::ArgError("read length cannot be less than zero.");
  }

  _numReads = 0; // this initial value will be replaced if the true num immediately
  _numReadsDetermined = false;

  _initialBlock = 0;
  _initialPos = 0;

  _posBlocks = new long[1];
  _posBlocks[0] = 0; // the first block starts with the zero coord
  _numPosBlocks = 1;

  _isOpen = false;
}


// PRIVATE CONTSTRUCTOR
ReadFilePriceqSingle::ReadFilePriceqSingle(ReadFilePriceqSingle* precursor, int initialBlock, long initialPos, long numReads, bool numReadsDetermined) :
  _initialBlock( initialBlock ),
  _initialPos( initialPos ),
  _numReads( numReads ),
  _numReadsDetermined( numReadsDetermined )
 {
   _readStart = precursor->_readStart;
   _readLength = precursor->_readLength;
   _invert = precursor->_invert;

  _hasReadWaiting = false;

  // run checks on filename
  string filename = precursor->getName();
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  _numPosBlocks = precursor->_numPosBlocks;
  _posBlocks = new long[_numPosBlocks];
  for (int n = 0; n < _numPosBlocks; ++n){ _posBlocks[n] = precursor->_posBlocks[n]; }
  _isOpen = false;
}




ReadFilePriceqSingle::~ReadFilePriceqSingle(){
  if (_isOpen){ close(); }
  delete[] _posBlocks;
  delete[] _filename;
}





void ReadFilePriceqSingle::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
  // determine the number of reads to include in each file
  long readsPerFile[numFiles];
  long currentTally = 0;
  long totalReads = numReads();
  for (int fileNum = 0; fileNum < numFiles; ++fileNum){
    long priorTally = currentTally;
    currentTally = (fileNum + 1) * totalReads / numFiles;
    readsPerFile[fileNum] = currentTally - priorTally;
  }
  // now create the files
  open();
  for (int fileNum = 0; fileNum < numFiles; ++fileNum){
    if (readsPerFile[fileNum] > 0){
      ReadFile* subFile = new ReadFilePriceqSingle(this, _currentBlock, _currentPos, readsPerFile[fileNum], true);
      putFilesHere->insert( subFile );
      skipReads( readsPerFile[fileNum] );
    }
  }
  close();
}


ReadFile* ReadFilePriceqSingle::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  open();
  skipReads( numReadsSkip );
  ReadFile* subFile = new ReadFilePriceqSingle(this, _currentBlock, _currentPos, numReadsKeep, true);
  close();
  return subFile;
}

ReadFilePriceqSingle* ReadFilePriceqSingle::copy(){
  ReadFilePriceqSingle* subFile = new ReadFilePriceqSingle(this, _initialBlock, _initialPos, _numReads, _numReadsDetermined);
  return subFile;
}



string ReadFilePriceqSingle::getName(){ return string(_filename); }


long ReadFilePriceqSingle::numReads(){
  // if the count is zero, then there is no guarantee that a count has taken place;
  // if zero is the true count, then re-counting will take little time to finish.
  #pragma omp critical (RFPQS)
  {
    if (! _numReadsDetermined){
      open();
      while ( hasRead() ){
	skipRead();
	_numReads++;
      }
      close();
      _numReadsDetermined = true;
    }
  }
  return _numReads;
}


long ReadFilePriceqSingle::ampliconSize(){
  throw AssemblyException::CallingError("A single-read fastq file does not have an amplicon size by definition.");
}

// open the file
void ReadFilePriceqSingle::open(){
  if (_isOpen){
    throw AssemblyException::CallingError("don't open RFFSingle, it is already open.");
  }
  if (_getReadsFile.is_open()){
    throw AssemblyException::CallingError("don't open RFFSingle, it is already open.");
  }

  _getReadsFile.open(_filename,ifstream::in);
  _openWithCarriers = false;
  _openNormalized = false;
  _hasReadWaiting = false;
  _currentBlock = 0;
  _currentPos = 0;
  _numReadsRead = 0;
  // move to the specified initial position
  while (_currentBlock < _initialBlock){
    _currentBlock++;
    _getReadsFile.seekg( _posBlocks[_currentBlock], ios::cur );
  }
  _getReadsFile.seekg( _initialPos, ios::cur );
  _currentPos = _initialPos;
  _isOpen = true;
}
// PE reads are not returned by this file type, so the bool is meaningless
// and the method is just forwarded.
void ReadFilePriceqSingle::open(bool getBothEnds){ open(); }

void ReadFilePriceqSingle::openWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFilePriceqSingle::openWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }

void ReadFilePriceqSingle::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFilePriceqSingle::openNormalized(bool getBothEnds){ openNormalized(); }

void ReadFilePriceqSingle::openNormWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFilePriceqSingle::openNormWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }


bool ReadFilePriceqSingle::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFSingle while it is closed.");
  }
  if (! _hasReadWaiting){ readFromFileHelper(); }
  return _hasReadWaiting;
}

ScoredSeq* ReadFilePriceqSingle::getRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't try to get reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq* newRead = getReadHelper();
    if (newRead == NULL){
      throw AssemblyException::LogicError("cannot get read if there is no read to get.");
    }
    return newRead;
  } else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFilePriceqSingle::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){ _hasReadWaiting = false; }
  else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFilePriceqSingle::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  for (long n = 0; n < numToSkip; ++n){
    skipRead();
  }
}



void ReadFilePriceqSingle::close(){
  // these haven't been accessed, so this is the only way to delete them
  if (_getReadsFile.is_open() ){
    _getReadsFile.close();
    _getReadsFile.clear();
  }
  while ( hasRead() ){ delete getRead(); }
  _hasReadWaiting = false;
  _isOpen = false;
}



void ReadFilePriceqSingle::readFromFileHelper(){
  if ( _getReadsFile.is_open() ){
    long localPos = 0; // all will be defined by this and then adjusted for the block pos
    while ( (! _getReadsFile.eof() ) and _getReadsFile.peek() != '@' ){
      string uselessLine;
      getline(_getReadsFile,uselessLine);
      localPos += uselessLine.size() + 1;
    }


    if ( (! _getReadsFile.eof() ) and _getReadsFile.peek() == '@' and ( (! _numReadsDetermined) or _numReadsRead < _numReads)  ){
      _numReadsRead++;

      long stepSize;

      string uselessLine; // just something that getline needs to write to; gets discarded.
      getline(_getReadsFile,uselessLine); // the first name line
      stepSize = uselessLine.size() + 1;
      _seqStart = localPos + stepSize;
      localPos = _seqStart;

      getline(_getReadsFile,uselessLine); // the sequence line
      _seqLen = uselessLine.size();
      localPos += _seqLen + 1;

      getline(_getReadsFile,uselessLine); // the second name line
      _scoreStart = localPos + uselessLine.size() + 1;
      localPos = _scoreStart;

      getline(_getReadsFile,uselessLine); // the quality score line
      _scoreLen = uselessLine.size();
      localPos += _scoreLen + 1;

      getline(_getReadsFile,uselessLine); // the third name line
      _linkStart = localPos + uselessLine.size() + 1;
      localPos = _linkStart;

      getline(_getReadsFile,uselessLine); // the link score line
      _linkLen = uselessLine.size();
      localPos += _linkLen + 1;

      // check if the block has run out of space; if so
      // create or move on to a new block
      long blockSpaceLeft = _maxBlockSize - _currentPos;
      if (localPos >= blockSpaceLeft){
	_currentBlock++;
	if (_currentBlock == _numPosBlocks){
	  long* newPB = new long[ _currentBlock + 1 ];
	  for (int n = 0; n < _numPosBlocks; ++n){ newPB[n] = _posBlocks[n]; }
	  newPB[_currentBlock] = _currentPos;
	  delete _posBlocks;
	  _posBlocks = newPB;
	  _numPosBlocks++;
	}
	_currentPos = 0;
      }
      // adjust the positions to global block positions
      _seqStart += _currentPos;
      _scoreStart += _currentPos;
      _linkStart += _currentPos;
      _currentPos += localPos;
      _hasReadWaiting = true;

    } else {
      _getReadsFile.close();
      _hasReadWaiting = false;
    }
  }
}



ScoredSeq* ReadFilePriceqSingle::getReadHelper(){
  if (! _hasReadWaiting){ readFromFileHelper(); }
  if ( _hasReadWaiting ){

    // create the reads and add them to the carrier
    if ( _seqLen != _scoreLen ){
      throw AssemblyException::LogicError("The sequence and score list are different lengths.");
    } else if ( _seqLen < _readStart + _readLength ){
      throw AssemblyException::LogicError("The read's specified dimensions extend off the edge of the available sequence.");
    }

    // create the read itself
    ReadFileIndex * newIndex;
    if (_readLength == 0){
      newIndex = new ReadFileIndex3start1len(_seqStart + _readStart,
					     _scoreStart + _readStart,
					     _linkStart + _readStart,
					     _seqLen - _readStart,
					     _currentBlock);
    } else {
      newIndex = new ReadFileIndex3start1len(_seqStart + _readStart,
					     _scoreStart + _readStart,
					     _linkStart + _readStart,
					     _readLength,
					     _currentBlock);
    }
    //ScoredSeq* seq = new Read(this, newIndex);
    ScoredSeq* seq = ReadFileCommunicator::makeRead(this, newIndex);

    // make the type adjustment if appropriate
    if (_openNormalized){ seq = new ScoredSeqNormalized( seq ); }
    if (_openWithCarriers){ seq = new ScoredSeqWithCarriers( seq, _numCarrierSets ); }
    _hasReadWaiting = false;

    return seq;

  } else { return 0; }
}


void ReadFilePriceqSingle::bufferQueue(Read * r){
  #pragma omp critical (RFPQSbufferFill)
  { _bufferSet.insert(r); }
}
void ReadFilePriceqSingle::bufferUnqueue(Read * r){
  #pragma omp critical (RFPQSbufferFill)
  { _bufferSet.erase(r); }
}


void ReadFilePriceqSingle::bufferFill(){
  bufferFill(bufferNotThreaded);
}

void ReadFilePriceqSingle::bufferFill(BufferThreadedness threadedness){
  char** rawSeqArray;
  char** rawScoresArray;
  char** rawLinksArray;
  long* rfiReadSize;
  Read** readArray;
  ScoredSeq** bufferArray;

  if (threadedness == bufferThreaded){
    rawSeqArray = new char*[ _bufferSet.size() + 1 ];
    rawScoresArray = new char*[ _bufferSet.size() + 1 ];
    rawLinksArray = new char*[ _bufferSet.size() + 1 ];
    rfiReadSize = new long[ _bufferSet.size() + 1 ];
    readArray = new Read*[ _bufferSet.size() + 1 ];
    bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];

    long seqCount = bufferFillSetupHelper(rawSeqArray, rawScoresArray, rawLinksArray, rfiReadSize, readArray);
    #pragma omp parallel for schedule(static)
    for (long seqN = 0; seqN < seqCount; ++seqN){
      bufferFillInterpretHelper(seqN, rawSeqArray, rawScoresArray, rawLinksArray, rfiReadSize, readArray, bufferArray);
    }
    for (long seqN = 0; seqN < seqCount; ++seqN){
      ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
    }
  } else {
    #pragma omp critical (RFPQSbuffer)
    {
      rawSeqArray = new char*[ _bufferSet.size() + 1 ];
      rawScoresArray = new char*[ _bufferSet.size() + 1 ];
      rawLinksArray = new char*[ _bufferSet.size() + 1 ];
      rfiReadSize = new long[ _bufferSet.size() + 1 ];
      readArray = new Read*[ _bufferSet.size() + 1 ];
      bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];

      long seqCount = bufferFillSetupHelper(rawSeqArray, rawScoresArray, rawLinksArray, rfiReadSize, readArray);
      for (long seqN = 0; seqN < seqCount; ++seqN){
	bufferFillInterpretHelper(seqN, rawSeqArray, rawScoresArray, rawLinksArray, rfiReadSize, readArray, bufferArray);
      }
      for (long seqN = 0; seqN < seqCount; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
    }
  }
  delete [] rawSeqArray;
  delete [] rawScoresArray;
  delete [] rawLinksArray;
  delete [] rfiReadSize;
  delete [] readArray;
  delete [] bufferArray;
}


bool ReadFilePriceqSingle::SortReadByRfiStart::operator() (Read* readA, Read* readB){
  return ReadFileCommunicator::getRfiRef(readA)->readStart() > ReadFileCommunicator::getRfiRef(readB)->readStart();
}

long ReadFilePriceqSingle::bufferFillSetupHelper(char** rawSeqArray, char** rawScoresArray, char** rawLinksArray,
						 long* rfiReadSize, Read** readArray){
  long seqCount = 0;
  int* rfiBlockArray = new int[ _bufferSet.size() + 1 ];
  long* countByBlock = new long[ _currentBlock+1 ];
  for (int n = 0; n < _currentBlock+1; ++n){ countByBlock[n] = 0; }
  for (set<Read*>::iterator buffIt = _bufferSet.begin(); buffIt != _bufferSet.end(); ++buffIt){
    Read* read = (*buffIt);
    readArray[seqCount] = read;
    int block = ReadFileCommunicator::getRfiRef(read)->blockNum();
    rfiBlockArray[seqCount] = block;
    countByBlock[ block ]++;
    seqCount++;
  }

  // sort the reads into bins by file block number
  long** blockToIndexes = new long*[ _currentBlock+1 ];
  long* indexByBlock = new long[ _currentBlock+1 ];
  for (int n = 0; n < _currentBlock+1; ++n){
    blockToIndexes[n] = new long[ countByBlock[n]+1 ];
    indexByBlock[n] = 0;
  }
  for (long sc = 0; sc < seqCount; ++sc){
    int block = rfiBlockArray[sc];
    blockToIndexes[block][ indexByBlock[block] ] = sc;
    indexByBlock[block]++;
  }

  // get the raw sequence input
  BlockPosFileCarrier* bpc = new BlockPosFileCarrier(this);
  long readArrayIndex = 0;
  for (int block = 0; block < _currentBlock+1; ++block){

    // sort the reads within each block by file start position
    vector<Read*> readVector;
    for (long index = 0; index < countByBlock[block]; ++index){ readVector.push_back( readArray[ blockToIndexes[block][index] ] ); }
    SortReadByRfiStart sorter; // sort by seq length, longest to shortest
    sort( readVector.begin(), readVector.end(), sorter );

    for (vector<Read*>::iterator readIt = readVector.begin(); readIt != readVector.end(); ++readIt){
      Read* read = *readIt;
      readArray[readArrayIndex] = read;
      ReadFileIndex* rfi = ReadFileCommunicator::getRfiRef(read);
      rfiReadSize[readArrayIndex] = rfi->readSize();
      // get the raw sequence
      bpc->seekBlockPos(rfi->blockNum(), rfi->readStart());
      rawSeqArray[readArrayIndex] = new char[ rfi->readSize() + 1];
      bpc->getString(rawSeqArray[readArrayIndex], streamsize( rfi->readSize() + 1 ) );
      // get the raw scores
      bpc->seekBlockPos(rfi->blockNum(), rfi->scoreStart());
      rawScoresArray[readArrayIndex] = new char[ rfi->scoreSize() + 1];
      bpc->getString(rawScoresArray[readArrayIndex], streamsize( rfi->scoreSize() + 1 ) );
      // get the raw links
      bpc->seekBlockPos(rfi->blockNum(), rfi->linkStart());
      rawLinksArray[readArrayIndex] = new char[ rfi->linkSize() + 1];
      bpc->getString(rawLinksArray[readArrayIndex], streamsize( rfi->linkSize() + 1 ) );
      readArrayIndex++;
    }
  }
  _bufferSet.clear();
  delete bpc;
  delete [] rfiBlockArray;
  delete [] countByBlock;
  for (int n = 0; n < _currentBlock+1; ++n){ delete [] blockToIndexes[n]; }
  delete [] blockToIndexes;
  delete [] indexByBlock;
  return seqCount;
}


void ReadFilePriceqSingle::bufferFillInterpretHelper(long seqN, char** rawSeqArray,
						     char** rawScoresArray, char** rawLinksArray,
						     long* rfiReadSize, Read** readArray, ScoredSeq** bufferArray){
  StringInterpreter* si = new StringInterpreter(rawSeqArray[seqN], rfiReadSize[seqN], rfiReadSize[seqN], &_okToSkip, _invert);
  delete [] rawSeqArray[seqN];
  // use this count for read size; the other strings may include bogus chars;
  // this is especially true if the file has carriage returns.
  float* scores = cstringToScores(rawScoresArray[seqN], si->_seqLen);
  delete [] rawScoresArray[seqN];
  float* links = cstringToScores(rawLinksArray[seqN], si->_seqLen - 1);
  delete [] rawLinksArray[seqN];
  // make the ScoredSeqShallow
  bufferArray[seqN] = ScoredSeq::repExposedSeq( si->_seqString, scores, links, si->_seqLen );
  delete si;
}


float* ReadFilePriceqSingle::cstringToScores(char* scoreCstring, long readSize){
  float* scores = new float[ readSize + 1 ];
  if (_invert){
    for (long n = 0; n < readSize; n++){
      scores[readSize - n - 1] = OutputFilePriceq::charToScore( scoreCstring[n] );
    }
  } else {
    for (long n = 0; n < readSize; n++){
      scores[n] = OutputFilePriceq::charToScore( scoreCstring[n] );
    }
  }
  return scores;
}


void ReadFilePriceqSingle::updateBuffer(){}
void ReadFilePriceqSingle::emptyBuffer(){}

// set up some initial conditions
ReadFilePriceqSingle::BlockPosFileCarrier::BlockPosFileCarrier(){}
ReadFilePriceqSingle::BlockPosFileCarrier::BlockPosFileCarrier(ReadFilePriceqSingle* host) :
  _host(host){
  if (_file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is already open."); }
  _file.open(host->_filename,ifstream::in);
  if (! _file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is still closed."); }
  _block = 0;
  _pos = 0;
}
ReadFilePriceqSingle::BlockPosFileCarrier::~BlockPosFileCarrier(){
  _file.close();
  _file.clear();
}


void ReadFilePriceqSingle::BlockPosFileCarrier::getString(char* cstringToFill, streamsize length){
  _file.get(cstringToFill, length );
  _pos += length - 1;
}

void ReadFilePriceqSingle::BlockPosFileCarrier::seekBlockPos(int newBlock, long newPos){
  // increment up or down until the correct block is reached
  _file.clear();
  while (_block > newBlock){
    _file.seekg( 0 - _pos, ios::cur );
    _pos = _host->_posBlocks[_block];
    _block--;
  }
  while (_block < newBlock){
    _block++;
    _file.seekg( _host->_posBlocks[_block] - _pos, ios::cur );
    _pos = 0;
  }
  // go to the correct spot in the correct block
  _file.seekg( newPos - _pos, ios::cur );
  _pos = newPos;
}




// THE RFI CLASS
ReadFilePriceqSingle::ReadFileIndex3start1len::ReadFileIndex3start1len(){}
ReadFilePriceqSingle::ReadFileIndex3start1len::ReadFileIndex3start1len(long readStart, long scoreStart, long linkStart, long readSize, int blockNum) :
  _readStart(readStart),
  _scoreStart(scoreStart),
  _linkStart(linkStart),
  _readSize(readSize),
  _blockNum(blockNum) {}
ReadFilePriceqSingle::ReadFileIndex3start1len::~ReadFileIndex3start1len(){}
long ReadFilePriceqSingle::ReadFileIndex3start1len::readStart(){ return _readStart; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::readSize(){ return _readSize; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::scoreStart(){ return _scoreStart; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::scoreSize(){ return _readSize; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::linkStart(){ return _linkStart; }
// there are one fewer links than nucleotides
long ReadFilePriceqSingle::ReadFileIndex3start1len::linkSize(){ return _readSize - 1; }
int ReadFilePriceqSingle::ReadFileIndex3start1len::blockNum(){ return _blockNum; }



#endif
