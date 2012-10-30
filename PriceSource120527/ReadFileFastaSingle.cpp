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

#ifndef READFILEFASTASINGLE_CPP
#define READFILEFASTASINGLE_CPP


#include "ReadFileFastaSingle.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <algorithm>
#include "ScoredSeqWithCarriers.h"
#include <vector>
#include <limits.h>
#include <cstring>
#include <omp.h>
using namespace::std;

long ReadFileFastaSingle::_maxBlockSize = LONG_MAX;

ReadFileFastaSingle::ReadFileFastaSingle(string filename) :
  _nucScore( 1.0 ),
  _invert(false),
  _splitMode(false),
  _maxReadSize(0)
 {
   constructorHelp();

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::FileError("RFFA fasta file does not exist."); }

  _numReads = 0; // this initial value will be replaced if the true num immediately
  _numReadsDetermined = false;

  _initialBlock = 0;
  _initialPos = 0;

  _posBlocks = new long[1];
  _posBlocks[0] = 0; // the first block starts with the zero coord
  _numPosBlocks = 1;

  _isOpen = false;
}


ReadFileFastaSingle::ReadFileFastaSingle(string filename, float nucScore, bool invert) :
  _nucScore( nucScore ),
  _invert(invert),
  _splitMode(false),
  _maxReadSize(0) {
   constructorHelp();

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::FileError("RFFA fasta file does not exist."); }

  _numReads = 0; // this initial value will be replaced if the true num immediately
  _numReadsDetermined = false;

  _initialBlock = 0;
  _initialPos = 0;

  _posBlocks = new long[1];
  _posBlocks[0] = 0; // the first block starts with the zero coord
  _numPosBlocks = 1;

  _isOpen = false;
}

ReadFileFastaSingle::ReadFileFastaSingle(string filename, float nucScore, long readSize, bool invert) :
  _nucScore( nucScore ),
  _invert(invert),
  _splitMode(true),
  _maxReadSize(readSize) {
   constructorHelp();

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::FileError("RFFA fasta file does not exist."); }

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
ReadFileFastaSingle::ReadFileFastaSingle(ReadFileFastaSingle* precursor, int initialBlock, long initialPos, long numReads, bool numReadsDetermined) :
  _initialBlock( initialBlock ),
  _initialPos( initialPos),
  _numReads( numReads ),
  _numReadsDetermined( numReadsDetermined )
 {
   constructorHelp();
  _nucScore = precursor->_nucScore;
  _invert = precursor->_invert;
  _splitMode = precursor->_splitMode;
  _maxReadSize = precursor->_maxReadSize;

  _hasReadWaiting = false;

  _minCountedScore = precursor->_minCountedScore;

  // run checks on filename
  string filename = precursor->getName();
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  _numPosBlocks = precursor->_numPosBlocks;
  _posBlocks = new long[_numPosBlocks];
  for (int n = 0; n < _numPosBlocks; ++n){ _posBlocks[n] = precursor->_posBlocks[n]; }
  _isOpen = false;
}


void ReadFileFastaSingle::constructorHelp(){
  _okToSkip.insert('\r');
  _okToSkip.insert('\n');
  _okToSkip.insert('\0');
  _okToSkip.insert(' ');
  _okToSkip.insert('-');
}


ReadFileFastaSingle::~ReadFileFastaSingle(){
  if (_isOpen){ close(); }
  delete[] _posBlocks;
  delete[] _filename;
}





void ReadFileFastaSingle::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
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
      ReadFile* subFile = new ReadFileFastaSingle(this, _currentBlock, _currentPos, readsPerFile[fileNum], true);
      putFilesHere->insert( subFile );
      skipReads( readsPerFile[fileNum] );
    }
  }
  close();
}


ReadFile* ReadFileFastaSingle::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  open();
  skipReads( numReadsSkip );
  ReadFile* subFile = new ReadFileFastaSingle(this, _currentBlock, _currentPos, numReadsKeep, true);
  close();
  return subFile;
}

ReadFileFastaSingle* ReadFileFastaSingle::copy(){
  ReadFileFastaSingle* subFile = new ReadFileFastaSingle(this, _initialBlock, _initialPos, _numReads, _numReadsDetermined);
  return subFile;
}



string ReadFileFastaSingle::getName(){ return string(_filename); }


long ReadFileFastaSingle::numReads(){
  // if the count is zero, then there is no guarantee that a count has taken place;
  // if zero is the true count, then re-counting will take little time to finish.
  #pragma omp critical (RFFAS)
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


long ReadFileFastaSingle::ampliconSize(){
  throw AssemblyException::CallingError("A single-read fastq file does not have an amplicon size by definition.");
}

// open the file
void ReadFileFastaSingle::open(){
  if (_isOpen){
    throw AssemblyException::CallingError("don't open RFFaSingle, it is already open.");
  }
  if (_getReadsFile.is_open()){
    throw AssemblyException::CallingError("don't open RFFaSingle, it is already open.");
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

  //find the first entry
  string uselessLine;
  while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '>' ){
    getline(_getReadsFile,uselessLine);
    _currentPos += uselessLine.size();
  }
}
// PE reads are not returned by this file type, so the bool is meaningless
// and the method is just forwarded.
void ReadFileFastaSingle::open(bool getBothEnds){ open(); }

void ReadFileFastaSingle::openWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFileFastaSingle::openWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }

void ReadFileFastaSingle::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFileFastaSingle::openNormalized(bool getBothEnds){ openNormalized(); }

void ReadFileFastaSingle::openNormWithCarriers(int numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFileFastaSingle::openNormWithCarriers(int numSets, bool getBothEnds){ openWithCarriers(numSets); }


bool ReadFileFastaSingle::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFSingle while it is closed.");
  }
  if (! _hasReadWaiting){ readFromFileHelper(); }
  return _hasReadWaiting;
}

ScoredSeq* ReadFileFastaSingle::getRead(){
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

void ReadFileFastaSingle::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){ _hasReadWaiting = false; }
  else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFileFastaSingle::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  for (long n = 0; n < numToSkip; ++n){
    skipRead();
  }
}



void ReadFileFastaSingle::close(){
  // these haven't been accessed, so this is the only way to delete them
  if (_getReadsFile.is_open() ){
    _getReadsFile.close();
    _getReadsFile.clear();
  }
  while ( hasRead() ){ delete getRead(); }
  _hasReadWaiting = false;
  _isOpen = false;
}



void ReadFileFastaSingle::readFromFileHelper(){
  if ( _getReadsFile.is_open() and (! _getReadsFile.eof() ) and 
       ( (! _numReadsDetermined) or _numReadsRead < _numReads) ){
    _numReadsRead++;

    long localPos = 0; // all will be defined by this and then adjusted for the block pos
    long stepSize;

    string uselessLine; // just something that getline needs to write to; gets discarded.
    getline(_getReadsFile,uselessLine); // the name line
    stepSize = uselessLine.size() + 1;
    _seqStart = localPos + stepSize;
    localPos = _seqStart;

    // find the next entry's carrot or the end of the file
    bool endOfSeq = false;
    _seqLen = 0;
    while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '>' ){
      string uselessLine; // must be defined in this scope to initialize as null
      getline(_getReadsFile,uselessLine); // the sequence lines
      _seqLen += uselessLine.size() + 1;
    }
    localPos += _seqLen;

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
    _currentPos += localPos;
    _hasReadWaiting = true;

  } else { 
    _getReadsFile.close();
    _hasReadWaiting = false;
  }
}



ScoredSeq* ReadFileFastaSingle::getReadHelper(){
  if (! _hasReadWaiting){ readFromFileHelper(); }
  if ( _hasReadWaiting ){

    // create the read itself
    ReadFileIndex * newIndex;
    newIndex = new ReadFileIndexJustSeq(_seqStart, _seqLen, _currentBlock);
    ScoredSeq* seq = ReadFileCommunicator::makeRead(this, newIndex);
    //ScoredSeq* seq = new Read(this, newIndex);

    // make the type adjustment if appropriate
    if (_openNormalized){ seq = new ScoredSeqNormalized( seq ); }
    if (_openWithCarriers){ seq = new ScoredSeqWithCarriers( seq, _numCarrierSets ); }
    _hasReadWaiting = false;

    return seq;

  } else { return NULL; }
}


void ReadFileFastaSingle::bufferQueue(Read * r){
  #pragma omp critical (RFFASbuffer)
  { _bufferSet.insert(r); }
}
void ReadFileFastaSingle::bufferUnqueue(Read * r){
  #pragma omp critical (RFFASbuffer)
  { _bufferSet.erase(r); }
}



void ReadFileFastaSingle::bufferFill(){
  bufferFill(bufferNotThreaded);
}

///*
void ReadFileFastaSingle::bufferFill(BufferThreadedness threadedness){
  char** rawSeqArray;
  long* rfiReadSize;
  Read** readArray;
  ScoredSeq** bufferArray;
  if (threadedness == bufferThreaded){
    rawSeqArray = new char*[ _bufferSet.size() + 1 ];
    rfiReadSize = new long[ _bufferSet.size() + 1 ];
    readArray = new Read*[ _bufferSet.size() + 1 ];
    bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];

    long seqCount = bufferFillSetupHelper(rawSeqArray, rfiReadSize, readArray);
    if ( _splitMode ){
      #pragma omp parallel for schedule(static)
      for (long seqN = 0; seqN < seqCount; ++seqN){
	bufferFillInterpretSplitHelper(seqN, rawSeqArray, rfiReadSize, readArray, bufferArray);
      }
    } else {
      #pragma omp parallel for schedule(static)
      for (long seqN = 0; seqN < seqCount; ++seqN){
	bufferFillInterpretHelper(seqN, rawSeqArray, rfiReadSize, readArray, bufferArray);
      }
    }
    for (long seqN = 0; seqN < seqCount; ++seqN){
      ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
    }
  } else {
    long seqCount;
    #pragma omp critical (RFFASbuffer)
    {
      rawSeqArray = new char*[ _bufferSet.size() + 1 ];
      rfiReadSize = new long[ _bufferSet.size() + 1 ];
      readArray = new Read*[ _bufferSet.size() + 1 ];
      bufferArray = new ScoredSeq*[ _bufferSet.size() + 1 ];
      seqCount = bufferFillSetupHelper(rawSeqArray, rfiReadSize, readArray);

      if ( _splitMode ){
	for (long seqN = 0; seqN < seqCount; ++seqN){
	  bufferFillInterpretSplitHelper(seqN, rawSeqArray, rfiReadSize, readArray, bufferArray);
	}
      } else {
	for (long seqN = 0; seqN < seqCount; ++seqN){
	  bufferFillInterpretHelper(seqN, rawSeqArray, rfiReadSize, readArray, bufferArray);
	}
      }

      for (long seqN = 0; seqN < seqCount; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
    }
      //}
  }
  delete [] rawSeqArray;
  delete [] rfiReadSize;
  delete [] readArray;
  delete [] bufferArray;
}

bool ReadFileFastaSingle::SortReadByRfiStart::operator() (Read* readA, Read* readB){
  return ReadFileCommunicator::getRfiRef(readA)->readStart() > ReadFileCommunicator::getRfiRef(readB)->readStart();
}

long ReadFileFastaSingle::bufferFillSetupHelper(char** rawSeqArray, long* rfiReadSize, Read** readArray){
  // set up for sorting the reads by file position
  long seqCount = 0;
  int* rfiBlockArray = new int[ _bufferSet.size() + 1 ];
  // get a count for the number of reads in each block
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

  // get the raw sequence input; the readArray will be re-sorted in the process
  BlockPosFileCarrier* bpc = new BlockPosFileCarrier(this);
  long readArrayIndex = 0;
  for (int block = 0; block < _currentBlock+1; ++block){

    // sort the reads within each block by file start position
    vector<Read*> readVector;
    for (long index = 0; index < countByBlock[block]; ++index){ readVector.push_back( readArray[ blockToIndexes[block][index] ] ); }
    SortReadByRfiStart sorter; // sort by ascending start position in the file
    sort( readVector.begin(), readVector.end(), sorter );

    for (vector<Read*>::iterator readIt = readVector.begin(); readIt != readVector.end(); ++readIt){
      Read* read = *readIt;
      readArray[readArrayIndex] = read;
      ReadFileIndex* rfi = ReadFileCommunicator::getRfiRef(read);
      rfiReadSize[readArrayIndex] = rfi->readSize();
      bpc->seekBlockPos(rfi->blockNum(), rfi->readStart());
      rawSeqArray[readArrayIndex] = new char[ rfi->readSize() + 1];
      bpc->getString(rawSeqArray[readArrayIndex], streamsize( rfi->readSize() + 1 ) );
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


void ReadFileFastaSingle::bufferFillInterpretHelper(long seqN, char** rawSeqArray, long* rfiReadSize,
						    Read** readArray, ScoredSeq** bufferArray){
  StringInterpreter* si = new StringInterpreter(rawSeqArray[seqN], rfiReadSize[seqN], rfiReadSize[seqN], &_okToSkip, _invert);
  delete [] rawSeqArray[seqN];
  bufferArray[seqN] = ScoredSeq::repExposedSeq(si->_seqString, _nucScore, si->_seqLen);
  delete si;
}
void ReadFileFastaSingle::bufferFillInterpretSplitHelper(long seqN, char** rawSeqArray, long* rfiReadSize,
							 Read** readArray, ScoredSeq** bufferArray){
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
  if (scoreChange >= seqLen){
    bufferArray[seqN] = ScoredSeq::repExposedSeq(si->_seqString, _nucScore, seqLen);
  } else if (scoreChange <= 0) {
    bufferArray[seqN] = ScoredSeq::repExposedSeq(si->_seqString, _nucScore / 2, seqLen);
  } else {
    float halfNucScore = _nucScore / 2;
    float* scores = new float[ seqLen ];
    // it is easier to just fill in a useless value at the end of links
    float* links = new float[ seqLen ];
    // first the full values...
    for (long n = 0; n < scoreChange; ++n){
      scores[n] = _nucScore;
      links[n] = _nucScore;
    }
    // ...now the half values
    for (long n = scoreChange; n < seqLen; ++n){
      scores[n] = halfNucScore;
      links[n] = halfNucScore;
    }
    bufferArray[seqN] = ScoredSeq::repExposedSeq(si->_seqString, scores, links, seqLen);
  }
  delete si;
}



void ReadFileFastaSingle::updateBuffer(){}
void ReadFileFastaSingle::emptyBuffer(){}

// set up some initial conditions
ReadFileFastaSingle::BlockPosFileCarrier::BlockPosFileCarrier(){}
ReadFileFastaSingle::BlockPosFileCarrier::BlockPosFileCarrier(ReadFileFastaSingle* host) :
  _host(host){
  if (_file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is already open."); }
  _file.open(host->_filename,ifstream::in);
  if (! _file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is still closed."); }
  _block = 0;
  _pos = 0;
}
ReadFileFastaSingle::BlockPosFileCarrier::~BlockPosFileCarrier(){
  _file.close();
  _file.clear();
}


void ReadFileFastaSingle::BlockPosFileCarrier::getString(char* cstringToFill, streamsize length){
  _file.get(cstringToFill, length, '\0');
  _pos += length - 1;
}

void ReadFileFastaSingle::BlockPosFileCarrier::seekBlockPos(int newBlock, long newPos){
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


ReadFileFastaSingle::ReadFileIndexJustSeq::ReadFileIndexJustSeq(){}
ReadFileFastaSingle::ReadFileIndexJustSeq::ReadFileIndexJustSeq(long readStart, long readSize, int blockNum) :
  _readStart(readStart),
  _readSize(readSize),
  _blockNum(blockNum) {
}
ReadFileFastaSingle::ReadFileIndexJustSeq::~ReadFileIndexJustSeq(){}
long ReadFileFastaSingle::ReadFileIndexJustSeq::readStart(){ return _readStart; }
long ReadFileFastaSingle::ReadFileIndexJustSeq::readSize(){ return _readSize; }
long ReadFileFastaSingle::ReadFileIndexJustSeq::scoreStart(){
  throw AssemblyException::CallingError("scoreStart is irrelevant to the RFI imp for RFFastaSingle.");
}
long ReadFileFastaSingle::ReadFileIndexJustSeq::scoreSize(){
  throw AssemblyException::CallingError("scoreSize is irrelevant to the RFI imp for RFFastaSingle.");
}
long ReadFileFastaSingle::ReadFileIndexJustSeq::linkStart(){
  throw AssemblyException::CallingError("linkStart is irrelevant to the RFI imp for RFFastaSingle.");
}
long ReadFileFastaSingle::ReadFileIndexJustSeq::linkSize(){
  throw AssemblyException::CallingError("linkSize is irrelevant to the RFI imp for RFFastaSingle.");
}
int ReadFileFastaSingle::ReadFileIndexJustSeq::blockNum(){ return _blockNum; }




#endif
