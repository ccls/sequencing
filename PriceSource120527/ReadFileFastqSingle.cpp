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


#ifndef READFILEFASTQSINGLE_CPP
#define READFILEFASTQSINGLE_CPP


#include "ReadFileFastqSingle.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <algorithm>
#include "ScoredSeqWithCarriers.h"
#include <vector>
#include <limits.h>
#include <cstring>
#include <cmath>
using namespace::std;

long ReadFileFastqSingle::_maxBlockSize = LONG_MAX;

ReadFileFastqSingle::ReadFileFastqSingle(){}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, Encoding encoding) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( 1.0 ),
  _invert(false) {

  setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, Encoding encoding, float countFactor) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(false) {
  setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}

ReadFileFastqSingle::ReadFileFastqSingle(bool invert, string filename, Encoding encoding, float countFactor) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(invert) {
   setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, Encoding encoding, int readStart, int readLength) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( 1.0 ),
  _invert(false) {
   setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, Encoding encoding, int readStart, int readLength, float countFactor) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( countFactor),
  _invert(false) {
   setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}


ReadFileFastqSingle::ReadFileFastqSingle(bool invert, string filename, Encoding encoding, int readStart, int readLength, float countFactor) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( countFactor),
  _invert(invert) {
   setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}


ReadFileFastqSingle::ReadFileFastqSingle(bool invert, long readLength, string filename, Encoding encoding, float countFactor) :
  _splitMode(true),
  _maxReadSize( readLength ),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(invert) {
   setEncodingHelper(encoding, _countFactor);
   constructorHelper(filename);
}


void ReadFileFastqSingle::constructorConstHelper(){
  _okToSkip.insert('\r');
  _okToSkip.insert('\n');
  _okToSkip.insert('\0');
}
void ReadFileFastqSingle::constructorHelper(string filename){
  constructorConstHelper();

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::FileError("RFFQ fastq file does not exist."); }

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



void ReadFileFastqSingle::setEncodingHelper(Encoding encoding, float countFactor){
  if (encoding == illuminaSystem){
    _encodingConverter = new ScoreConverterIllumina(countFactor);
  } else if (encoding == fastqSystem){
    _encodingConverter = new ScoreConverterFastq(countFactor);
  } else {
    throw AssemblyException::ArgError("RFFQS::setEncoding, encoding enum type is unknown.");
  }
}


// PRIVATE CONTSTRUCTOR
ReadFileFastqSingle::ReadFileFastqSingle(ReadFileFastqSingle* precursor, int initialBlock, long initialPos, long numReads, bool numReadsDetermined) :
  _initialBlock( initialBlock ),
  _initialPos( initialPos ),
  _numReads( numReads ),
  _numReadsDetermined( numReadsDetermined )
 {
   constructorConstHelper();
   setEncodingHelper(precursor->_encodingConverter->getEncoding(), precursor->_countFactor);
   _readStart = precursor->_readStart;
   _readLength = precursor->_readLength;
   _countFactor = precursor->_countFactor;
   _invert = precursor->_invert;
   _splitMode = precursor->_splitMode;
   if (_splitMode){ _maxReadSize = precursor->_maxReadSize; }

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




ReadFileFastqSingle::~ReadFileFastqSingle(){
  if (_isOpen){ close(); }
  delete [] _posBlocks;
  delete [] _filename;
  delete _encodingConverter;
}





void ReadFileFastqSingle::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
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
      ReadFile* subFile = new ReadFileFastqSingle(this, _currentBlock, _currentPos, readsPerFile[fileNum], true);
      putFilesHere->insert( subFile );
      skipReads( readsPerFile[fileNum] );
    }
  }
  close();
}


ReadFile* ReadFileFastqSingle::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  open();
  skipReads( numReadsSkip );
  ReadFile* subFile = new ReadFileFastqSingle(this, _currentBlock, _currentPos, numReadsKeep, true);
  close();
  return subFile;
}

ReadFileFastqSingle* ReadFileFastqSingle::copy(){
  ReadFileFastqSingle* subFile = new ReadFileFastqSingle(this, _initialBlock, _initialPos, _numReads, _numReadsDetermined);
  return subFile;
}



string ReadFileFastqSingle::getName(){ return string(_filename); }


long ReadFileFastqSingle::numReads(){
  // if the count is zero, then there is no guarantee that a count has taken place;
  // if zero is the true count, then re-counting will take little time to finish.
  #pragma omp critical (RFFQS)
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


long ReadFileFastqSingle::ampliconSize(){
  throw AssemblyException::CallingError("A single-read fastq file does not have an amplicon size by definition.");
}

// open the file
void ReadFileFastqSingle::open(){
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
void ReadFileFastqSingle::open(bool getBothEnds){ open(); }

void ReadFileFastqSingle::openWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFileFastqSingle::openWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }

void ReadFileFastqSingle::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFileFastqSingle::openNormalized(bool getBothEnds){ openNormalized(); }

void ReadFileFastqSingle::openNormWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFileFastqSingle::openNormWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }


bool ReadFileFastqSingle::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFSingle while it is closed.");
  }
  if (! _hasReadWaiting){ readFromFileHelper(); }
  return _hasReadWaiting;
}

ScoredSeq* ReadFileFastqSingle::getRead(){
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

void ReadFileFastqSingle::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){ _hasReadWaiting = false; }
  else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFileFastqSingle::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  for (long n = 0; n < numToSkip; ++n){
    skipRead();
  }
}



void ReadFileFastqSingle::close(){
  // these haven't been accessed, so this is the only way to delete them
  if (_getReadsFile.is_open() ){
    _getReadsFile.close();
    _getReadsFile.clear();
  }
  while ( hasRead() ){ delete getRead(); }
  _hasReadWaiting = false;
  _isOpen = false;
}



void ReadFileFastqSingle::readFromFileHelper(){
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
      _currentPos += localPos;
      _hasReadWaiting = true;

    } else {
      _getReadsFile.close();
      _hasReadWaiting = false;
    }
  }
}



ScoredSeq* ReadFileFastqSingle::getReadHelper(){
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
      newIndex = new ReadFileIndex2start1len(_seqStart + _readStart, _scoreStart + _readStart, _seqLen - _readStart, _currentBlock);
    } else {
      newIndex = new ReadFileIndex2start1len(_seqStart + _readStart, _scoreStart + _readStart, _readLength, _currentBlock);
    }
    ScoredSeq* seq = ReadFileCommunicator::makeRead(this, newIndex);
    //ScoredSeq* seq = new Read(this, newIndex);

    // make the type adjustment if appropriate
    if (_openNormalized){ seq = new ScoredSeqNormalized( seq ); }
    if (_openWithCarriers){ seq = new ScoredSeqWithCarriers( seq, _numCarrierSets ); }
    _hasReadWaiting = false;

    return seq;

  } else { return 0; }
}


void ReadFileFastqSingle::bufferQueue(Read * r){
  #pragma omp critical (RFFQSbuffer)
  { _bufferSet.insert(r); }
}
void ReadFileFastqSingle::bufferUnqueue(Read * r){
  #pragma omp critical (RFFQSbuffer)
  { _bufferSet.erase(r); }
}


void ReadFileFastqSingle::bufferFill(){
  bufferFill(bufferNotThreaded);
}

/* OLD VERSION
void ReadFileFastqSingle::bufferFill(BufferThreadedness threadedness){
  if (threadedness == bufferNotThreaded){
    //#pragma omp critical (RFFQSbuffer)
    { bufferFillHelper(threadedness); }
  } else { bufferFillHelper(threadedness); }
}
*/

void ReadFileFastqSingle::bufferFill(BufferThreadedness threadedness){
  char** rawSeqArray;
  char** rawScoresArray;
  long* rfiReadSize;
  Read** readArray;
  ScoredSeq** bufferArray;
  if (threadedness == bufferThreaded){
    rawSeqArray = new char*[ _bufferSet.size() + 1 ];
    rawScoresArray = new char*[ _bufferSet.size() + 1 ];
    rfiReadSize = new long[ _bufferSet.size() + 1 ];
    readArray = new Read*[ _bufferSet.size() + 1 ];
    bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];

    long seqCount = bufferFillSetupHelper(rawSeqArray, rawScoresArray, rfiReadSize, readArray);
    if ( _splitMode ){
      #pragma omp parallel for schedule(static)
      for (long seqN = 0; seqN < seqCount; ++seqN){
	bufferFillInterpretSplitHelper(seqN, rawSeqArray, rawScoresArray, rfiReadSize, readArray, bufferArray);
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (long seqN = 0; seqN < seqCount; ++seqN){
	bufferFillInterpretHelper(seqN, rawSeqArray, rawScoresArray, rfiReadSize, readArray, bufferArray);
      }
    }
    for (long seqN = 0; seqN < seqCount; ++seqN){
      ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
    }
  } else {
    long seqCount;
    #pragma omp critical (RFFQSbuffer)
    {
      rawSeqArray = new char*[ _bufferSet.size() + 1 ];
      rawScoresArray = new char*[ _bufferSet.size() + 1 ];
      rfiReadSize = new long[ _bufferSet.size() + 1 ];
      readArray = new Read*[ _bufferSet.size() + 1 ];
      bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];
      seqCount = bufferFillSetupHelper(rawSeqArray, rawScoresArray, rfiReadSize, readArray);

      if ( _splitMode ){
	for (long seqN = 0; seqN < seqCount; ++seqN){
	  bufferFillInterpretSplitHelper(seqN, rawSeqArray, rawScoresArray, rfiReadSize, readArray, bufferArray);
	}
      } else {
	for (long seqN = 0; seqN < seqCount; ++seqN){
	  bufferFillInterpretHelper(seqN, rawSeqArray, rawScoresArray, rfiReadSize, readArray, bufferArray);
	}
      }

      for (long seqN = 0; seqN < seqCount; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
    }
  }
  delete [] rawSeqArray;
  delete [] rawScoresArray;
  delete [] rfiReadSize;
  delete [] readArray;
  delete [] bufferArray;
}


bool ReadFileFastqSingle::SortReadByRfiStart::operator() (Read* readA, Read* readB){
  return ReadFileCommunicator::getRfiRef(readA)->readStart() > ReadFileCommunicator::getRfiRef(readB)->readStart();
}

long ReadFileFastqSingle::bufferFillSetupHelper(char** rawSeqArray,
						char** rawScoresArray,
						long* rfiReadSize,
						Read** readArray){
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

void ReadFileFastqSingle::bufferFillInterpretHelper(long seqN,
						    char** rawSeqArray,
						    char** rawScoresArray,
						    long* rfiReadSize,
						    Read** readArray,
						    ScoredSeq** bufferArray){
  StringInterpreter* si = new StringInterpreter(rawSeqArray[seqN], rfiReadSize[seqN], rfiReadSize[seqN], &_okToSkip, _invert);
  delete [] rawSeqArray[seqN];
  // use this count for read size; the other strings may include bogus chars;
  // this is especially true if the file has carriage returns.
  float* scores = _encodingConverter->cstringToScores(rawScoresArray[seqN], si->_seqLen, _invert, si->_seqString);
  delete [] rawScoresArray[seqN];
  // create a links array too
  float* links = new float[ si->_seqLen ];
  for (long n = 0; n < si->_seqLen; ++n){ links[n] = _countFactor; }
  // make the ScoredSeqShallow
  bufferArray[seqN] = ScoredSeq::repExposedSeq( si->_seqString, scores, links, si->_seqLen );
  delete si;
}
void ReadFileFastqSingle::bufferFillInterpretSplitHelper(long seqN,
							 char** rawSeqArray,
							 char** rawScoresArray,
							 long* rfiReadSize,
							 Read** readArray,
							 ScoredSeq** bufferArray){
  StringInterpreter* si = new StringInterpreter(rawSeqArray[seqN], rfiReadSize[seqN], rfiReadSize[seqN], &_okToSkip, _invert);
  delete [] rawSeqArray[seqN];

  // now adjsut the parameters to get only part of the sequence and adjusted scores
  long seqLen = si->_seqLen;
  if (seqLen > _maxReadSize){
    // null-char-terminate at the end of the legit seq
    si->_seqString[ _maxReadSize ] = '\0';
    seqLen = _maxReadSize;
  }

  // figure out where the scores would change to half their value and adjust as necessary
  long scoreChange = si->_seqLen - seqLen;
  if (scoreChange < 0){ scoreChange = 0; }
  // i will iterate to scoreChange below, so make sure it isn't bigger than seqLen
  if (scoreChange > seqLen){ scoreChange = seqLen; }

  // get the scores; the ones after scoreChange will still need to be halved
  float* scores = _encodingConverter->cstringToScores(rawScoresArray[seqN], seqLen, _invert, si->_seqString);
  delete [] rawScoresArray[seqN];
  // links has an extra element just because that is faster below than excluding
  // the last element and treating it differently
  float* links = new float[ seqLen ];

  // first the full values...
  float halfCount = _countFactor / 2;
  for (long n = 0; n < scoreChange; ++n){ links[n] = _countFactor; }
  // ...now the half values
  for (long n = scoreChange; n < seqLen; ++n){
    scores[n] = scores[n] / 2;
    links[n] = halfCount;
  }

  // make the ScoredSeqShallow
  bufferArray[seqN] = ScoredSeq::repExposedSeq( si->_seqString, scores, links, seqLen );
  delete si;
}



ReadFileFastqSingle::ScoreConverter::~ScoreConverter(){}

ReadFileFastqSingle::ScoreConverterFastq::ScoreConverterFastq(float countFactor){
  // calclulate the scores ahead of time
  for (int n = 0; n < 94; ++n){
    _scores[n] = countFactor * (float(1) - pow( float(10), float(0 - n) / float(10) ) );
  }
}
ReadFileFastqSingle::ScoreConverterFastq::~ScoreConverterFastq(){}
ReadFileFastqSingle::Encoding ReadFileFastqSingle::ScoreConverterFastq::getEncoding(){ return fastqSystem; }

float* ReadFileFastqSingle::ScoreConverterFastq::cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString){
  float* scores = new float[ readSize + 1 ];
  if (invert){
    long readSizeM1 = readSize - 1;
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[readSizeM1 - n], seqString, n); }
  } else {
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[n], seqString, n); }
  }
  return scores;
}
inline float ReadFileFastqSingle::ScoreConverterFastq::charToScore(char nuc, char* seqString, long seqPos){
  switch(nuc){
  case '!':
    seqString[seqPos] = 'N';
    //return _scores[0];
    return 0;
  case '"':
    seqString[seqPos] = 'N';
    //return _scores[1];
    return 0;
  case '#':
    seqString[seqPos] = 'N';
    //return _scores[2];
    return 0;
  case '$': return _scores[3];
  case '%': return _scores[4];
  case '&': return _scores[5];
  case '\'': return _scores[6];
  case '(': return _scores[7];
  case ')': return _scores[8];
  case '*': return _scores[9];
  case '+': return _scores[10];
  case ',': return _scores[11];
  case '-': return _scores[12];
  case '.': return _scores[13];
  case '/': return _scores[14];
  case '0': return _scores[15];
  case '1': return _scores[16];
  case '2': return _scores[17];
  case '3': return _scores[18];
  case '4': return _scores[19];
  case '5': return _scores[20];
  case '6': return _scores[21];
  case '7': return _scores[22];
  case '8': return _scores[23];
  case '9': return _scores[24];
  case ':': return _scores[25];
  case ';': return _scores[26];
  case '<': return _scores[27];
  case '=': return _scores[28];
  case '>': return _scores[29];
  case '?': return _scores[30];
  case '@': return _scores[31];
  case 'A': return _scores[32];
  case 'B': return _scores[33];
  case 'C': return _scores[34];
  case 'D': return _scores[35];
  case 'E': return _scores[36];
  case 'F': return _scores[37];
  case 'G': return _scores[38];
  case 'H': return _scores[39];
  case 'I': return _scores[40];
  case 'J': return _scores[41];
  case 'K': return _scores[42];
  case 'L': return _scores[43];
  case 'M': return _scores[44];
  case 'N': return _scores[45];
  case 'O': return _scores[46];
  case 'P': return _scores[47];
  case 'Q': return _scores[48];
  case 'R': return _scores[49];
  case 'S': return _scores[50];
  case 'T': return _scores[51];
  case 'U': return _scores[52];
  case 'V': return _scores[53];
  case 'W': return _scores[54];
  case 'X': return _scores[55];
  case 'Y': return _scores[56];
  case 'Z': return _scores[57];
  case '[': return _scores[58];
  case '\\': return _scores[59];
  case ']': return _scores[60];
  case '^': return _scores[61];
  case '_': return _scores[62];
  case '`': return _scores[63];
  case 'a': return _scores[64];
  case 'b': return _scores[65];
  case 'c': return _scores[66];
  case 'd': return _scores[67];
  case 'e': return _scores[68];
  case 'f': return _scores[69];
  case 'g': return _scores[70];
  case 'h': return _scores[71];
  case 'i': return _scores[72];
  case 'j': return _scores[73];
  case 'k': return _scores[74];
  case 'l': return _scores[75];
  case 'm': return _scores[76];
  case 'n': return _scores[77];
  case 'o': return _scores[78];
  case 'p': return _scores[79];
  case 'q': return _scores[80];
  case 'r': return _scores[81];
  case 's': return _scores[82];
  case 't': return _scores[83];
  case 'u': return _scores[84];
  case 'v': return _scores[85];
  case 'w': return _scores[86];
  case 'x': return _scores[87];
  case 'y': return _scores[88];
  case 'z': return _scores[89];
  case '{': return _scores[90];
  case '|': return _scores[91];
  case '}': return _scores[92];
  case '~': return _scores[93];
  default:
    cerr << "Bad char: " << nuc << endl;
    throw AssemblyException::FileError("Fastq-encoded file contained an illegal score character.");
  }
}




ReadFileFastqSingle::ScoreConverterIllumina::ScoreConverterIllumina(float countFactor){
  // calclulate the scores ahead of time
  for (int n = 0; n < 63; ++n){
    _scores[n] = countFactor * (float(1) - pow( float(10), float(0 - n) / float(10) ) );
  }
}
ReadFileFastqSingle::ScoreConverterIllumina::~ScoreConverterIllumina(){}
ReadFileFastqSingle::Encoding ReadFileFastqSingle::ScoreConverterIllumina::getEncoding(){ return illuminaSystem; }

float* ReadFileFastqSingle::ScoreConverterIllumina::cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString){
  float* scores = new float[ readSize + 1 ];
  if (invert){
    long readSizeM1 = readSize - 1;
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[readSizeM1 - n], seqString, n); }
  } else {
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[n], seqString, n); }
  }
  return scores;
}
inline float ReadFileFastqSingle::ScoreConverterIllumina::charToScore(char nuc, char* seqString, long seqPos){
  switch(nuc){
  case '@':
    seqString[seqPos] = 'N';
    //return _scores[0];
    return 0;
  case 'A':
    seqString[seqPos] = 'N';
    //return _scores[1];
    return 0;
  case 'B':
    seqString[seqPos] = 'N';
    //return _scores[2];
    return 0;
  case 'C': return _scores[3];
  case 'D': return _scores[4];
  case 'E': return _scores[5];
  case 'F': return _scores[6];
  case 'G': return _scores[7];
  case 'H': return _scores[8];
  case 'I': return _scores[9];
  case 'J': return _scores[10];
  case 'K': return _scores[11];
  case 'L': return _scores[12];
  case 'M': return _scores[13];
  case 'N': return _scores[14];
  case 'O': return _scores[15];
  case 'P': return _scores[16];
  case 'Q': return _scores[17];
  case 'R': return _scores[18];
  case 'S': return _scores[19];
  case 'T': return _scores[20];
  case 'U': return _scores[21];
  case 'V': return _scores[22];
  case 'W': return _scores[23];
  case 'X': return _scores[24];
  case 'Y': return _scores[25];
  case 'Z': return _scores[26];
  case '[': return _scores[27];
  case '\\': return _scores[28];
  case ']': return _scores[29];
  case '^': return _scores[30];
  case '_': return _scores[31];
  case '`': return _scores[32];
  case 'a': return _scores[33];
  case 'b': return _scores[34];
  case 'c': return _scores[35];
  case 'd': return _scores[36];
  case 'e': return _scores[37];
  case 'f': return _scores[38];
  case 'g': return _scores[39];
  case 'h': return _scores[40];
  case 'i': return _scores[41];
  case 'j': return _scores[42];
  case 'k': return _scores[43];
  case 'l': return _scores[44];
  case 'm': return _scores[45];
  case 'n': return _scores[46];
  case 'o': return _scores[47];
  case 'p': return _scores[48];
  case 'q': return _scores[49];
  case 'r': return _scores[50];
  case 's': return _scores[51];
  case 't': return _scores[52];
  case 'u': return _scores[53];
  case 'v': return _scores[54];
  case 'w': return _scores[55];
  case 'x': return _scores[56];
  case 'y': return _scores[57];
  case 'z': return _scores[58];
  case '{': return _scores[59];
  case '|': return _scores[60];
  case '}': return _scores[61];
  case '~': return _scores[62];
  default:
    cerr << "Bad char: " << nuc << endl;
    throw AssemblyException::FileError("Illumina-encoded (Phred+64) file contained an illegal score character.");
  }
}




void ReadFileFastqSingle::updateBuffer(){}
void ReadFileFastqSingle::emptyBuffer(){}

// set up some initial conditions
ReadFileFastqSingle::BlockPosFileCarrier::BlockPosFileCarrier(){}
ReadFileFastqSingle::BlockPosFileCarrier::BlockPosFileCarrier(ReadFileFastqSingle* host) :
  _host(host){
  if (_file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is already open."); }
  _file.open(host->_filename,ifstream::in);
  if (! _file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is still closed."); }
  _block = 0;
  _pos = 0;
}
ReadFileFastqSingle::BlockPosFileCarrier::~BlockPosFileCarrier(){
  _file.close();
  _file.clear();
}


void ReadFileFastqSingle::BlockPosFileCarrier::getString(char* cstringToFill, streamsize length){
  _file.get(cstringToFill, length );
  _pos += length - 1;
}

void ReadFileFastqSingle::BlockPosFileCarrier::seekBlockPos(int newBlock, long newPos){
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
ReadFileFastqSingle::ReadFileIndex2start1len::ReadFileIndex2start1len(){}
ReadFileFastqSingle::ReadFileIndex2start1len::ReadFileIndex2start1len(long readStart, long scoreStart, long readSize, int blockNum) :
  _readStart(readStart),
  _readSize(readSize),
  _scoreStart(scoreStart),
  _blockNum(blockNum) {
}
ReadFileFastqSingle::ReadFileIndex2start1len::~ReadFileIndex2start1len(){}
long ReadFileFastqSingle::ReadFileIndex2start1len::readStart(){ return _readStart; }
long ReadFileFastqSingle::ReadFileIndex2start1len::readSize(){ return _readSize; }
long ReadFileFastqSingle::ReadFileIndex2start1len::scoreStart(){ return _scoreStart; }
long ReadFileFastqSingle::ReadFileIndex2start1len::scoreSize(){ return _readSize; }
long ReadFileFastqSingle::ReadFileIndex2start1len::linkStart(){
  throw AssemblyException::CallingError("linkStart is irrelevant to the RFI imp for RFFastqSingle.");
}
long ReadFileFastqSingle::ReadFileIndex2start1len::linkSize(){
  throw AssemblyException::CallingError("linkSize is irrelevant to the RFI imp for RFFastqSingle.");
}
int ReadFileFastqSingle::ReadFileIndex2start1len::blockNum(){ return _blockNum; }





#endif
