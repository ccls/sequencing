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



#ifndef MAPREADSETUPHELPER_CPP
#define MAPREADSETUPHELPER_CPP

#include "MapReadSetupHelper.h"
#include "Read.h"
#include "ScoredSeqSubseq.h"
#include "ScoredSeqNested.h"
#include "ScoredSeqMonoScore.h"
#include <typeinfo>
#include <sstream>
#include <omp.h>
#include <limits.h>
using namespace::std;



MapReadSetupHelper::~MapReadSetupHelper(){}


MrsMainHelper::MrsMainHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf, set<ScoredSeqWithCarriers*>* doNotExtendContigs) :
  _numReadFiles(numReadFiles),
  _rpf(rpf) {// = rpf->copy();
  // figure out the total file sizes for stdout update and set up for threaded file reading
  _numActiveFiles = 0;
  _totalFileSize = 0;
  _doNotExtendContigs.insert(doNotExtendContigs->begin(), doNotExtendContigs->end());

  _readFileArray = new ParameterizedReadFile*[ _numReadFiles + 1 ];
  for (int n = 0; n < _numReadFiles; ++n){ _readFileArray[n] = readFileArray[n]; }
  _fileIsActive = new bool[ _numReadFiles ]; // keeps track of how many files are still being read
  _fileToPmap = new ParamsMapping*[ _numReadFiles ];

  #pragma omp parallel for schedule(dynamic)
  for (int n = 0; n < _numReadFiles; ++n){
    long numReads = _readFileArray[n]->numReads();
    //#pragma omp critical (ECmapReads)
    {
      #pragma omp atomic
      _totalFileSize += numReads;
      // pre-open the files that have reads
      if (numReads == 0){
        #pragma omp critical (ECmapReads)
	{
	  _fileIsActive[n] = false;
	  _fileToPmap[n] = 0;
	}
      } else {
        #pragma omp critical (ECmapReads)
	{
	  _fileToPmap[n] = _readFileArray[n]->getMapParams();
	  _fileIsActive[n] = true;
	}
        #pragma omp atomic
	_numActiveFiles++;
	_readFileArray[n]->openNormWithCarriers(2,false);
      }
    }
  }
}
MrsMainHelper::~MrsMainHelper(){
  for (int n = 0; n < _numReadFiles; ++n){ 
    _readFileArray[n]->close();
    delete _fileToPmap[n];
  }
  delete [] _readFileArray;
  delete [] _fileIsActive;
  delete [] _fileToPmap;
}

int MrsMainHelper::numFilesTotal(){ return _numReadFiles; }
int MrsMainHelper::numFilesActive(){ return _numActiveFiles; }
bool MrsMainHelper::isFileActive(int fileNum){ return _fileIsActive[fileNum]; }
void MrsMainHelper::triggerFiles(){
  for (int n = 0; n < _numReadFiles; ++n){
    ReadFileCommunicator::bufferFill(_readFileArray[n], ReadFile::bufferThreaded);
  }
}


long MrsMainHelper::totalReadCount(){ return _totalFileSize; }
ExtendJobMapper* MrsMainHelper::getNewEjm(set<ScoredSeqWithCarriers*>* parallelMapSet, int fileNum){
  return new ExtendJobMapperNoLinks( parallelMapSet, _fileToPmap[fileNum]->fractId(), _readFileArray[fileNum]->ampliconSize() );
}


void MrsMainHelper::getPairedReadArrays(ScoredSeqWithCarriers*** readArrays, ScoredSeqWithCarriers*** pairArrays,
					long* numPairArray, int* fileNumArray, int numFiles){

  // thread the reading from files by file
  #pragma omp parallel for schedule(dynamic)
  for (int fileIndex = 0; fileIndex < numFiles; ++fileIndex){
    int fileNum = fileNumArray[fileIndex];
    ParameterizedReadFile* readFile = _readFileArray[fileNum];
    long arrayTargetSize = _fileToPmap[fileNum]->stepSize() / 2;
    if (arrayTargetSize == 0){ arrayTargetSize = 1; }
    ScoredSeqWithCarriers** readArray = new ScoredSeqWithCarriers*[ arrayTargetSize ];
    ScoredSeqWithCarriers** pairArray = new ScoredSeqWithCarriers*[ arrayTargetSize ];
    readArrays[fileIndex] = readArray;
    pairArrays[fileIndex] = pairArray;

    long pairCount = 0;
    while ( pairCount < arrayTargetSize and readFile->hasRead() ){
      ScoredSeqWithCarriers* newRead = dynamic_cast<ScoredSeqWithCarriers*>( readFile->getRead() );
      // only paired-end reads will be mapped, if the file has a mixture
      if ( newRead->hasPairedEnd() ){
	ScoredSeqWithCarriers* newPair = dynamic_cast<ScoredSeqWithCarriers*>( newRead->getPairedEnd() );
	readArray[pairCount] = newRead;
	pairArray[pairCount] = newPair;
	newRead->bottomBuffer();
	newPair->bottomBuffer();
	++pairCount;
      } else { newRead->deepDelete(); }
    }
    numPairArray[fileIndex] = pairCount;
  }
}


long MrsMainHelper::decideMappability(bool* shouldReadsMap, bool* shouldPairsMap, int fileNum,
				      ScoredSeqWithCarriers** reads, ScoredSeqWithCarriers** pairs, long rpCount){
  long numUnmappable = 0;
  for (long n = 0; n < rpCount; ++n){
    // assumes that testKillsPair is true for this test
    if (_rpf->isReadOk(reads[n]) and _rpf->isReadOk(pairs[n])){
      shouldReadsMap[n] = true;
      shouldPairsMap[n] = true;
    } else {
      shouldReadsMap[n] = false;
      shouldPairsMap[n] = false;
      numUnmappable += 2;
    }
    /*
    shouldReadsMap[n] = _rpf->isPairOk(reads[n], pairs[n]);
    shouldPairsMap[n] = shouldReadsMap[n];
    */
  }
  return numUnmappable;
}


// TEMPORARILY, THIS IS A FORWARDING METHOD - EVENTUALLY, I WILL DELETE IT ENTRIELY
ScoredSeqWithCarriers*** MrsMainHelper::getNewReadArrays(int* fileNumArray, int numFiles){
  ScoredSeqWithCarriers*** readArrays = new ScoredSeqWithCarriers**[ numFiles + 1 ];
  ScoredSeqWithCarriers*** pairArrays = new ScoredSeqWithCarriers**[ numFiles + 1 ];
  long* pairNumArray = new long[ numFiles + 1 ];

  getPairedReadArrays(readArrays, pairArrays, pairNumArray, fileNumArray, numFiles);

  ScoredSeqWithCarriers*** allSeqsToMapArrays = new ScoredSeqWithCarriers**[ numFiles + 1 ];
  for (int n = 0; n < numFiles; ++n){
    long readCounter = 0;
    ScoredSeqWithCarriers** seqsToMapArray = new ScoredSeqWithCarriers*[ pairNumArray[n] * 2 + 1 ];
    for (long readN = 0; readN < pairNumArray[n]; ++readN){
      seqsToMapArray[ readCounter ] = readArrays[n][readN];
      seqsToMapArray[ readCounter+1 ] = pairArrays[n][readN];
      readCounter += 2;
    }
    // end each array with a null
    seqsToMapArray[readCounter] = NULL;
    allSeqsToMapArrays[n] = seqsToMapArray;
    delete [] readArrays[n];
    delete [] pairArrays[n];
  }
  delete [] readArrays;
  delete [] pairArrays;
  delete [] pairNumArray;

  return allSeqsToMapArrays;
}

void MrsMainHelper::updateActiveFiles(){
  _numActiveFiles = 0;
  for (int n = 0; n < _numReadFiles; ++n){
    if (_fileIsActive[n]){
      if ( _readFileArray[n]->hasRead() ){ _numActiveFiles++; }
      else { _fileIsActive[n] = false; }
    }
  }
}
bool MrsMainHelper::includeContigCarrier(ScoredSeqWithCarriers* carrier){ return _doNotExtendContigs.count(carrier) == 0; }






MrsAddendumHelper::MrsAddendumHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf,
						  vector<ScoredSeqWithCarriers*>** hitReadsByFile) :
  _rpf(rpf),
  _numReadFiles(numReadFiles) {
  _numActiveFiles = 0;
  _totalNumQueries = 0;

  _readFileArray = new ParameterizedReadFile*[ _numReadFiles + 1 ];
  for (int n = 0; n < _numReadFiles; ++n){ _readFileArray[n] = readFileArray[n]; }

  // make sets of the paired-ends of the mapped reads, THOSE will be mapped
  _readQueries = new ScoredSeqWithCarriers**[ _numReadFiles ];
  _pairQueries = new ScoredSeqWithCarriers**[ _numReadFiles ];
  _numPairs = new long[ _numReadFiles ];
  _readPairCounter = new long[ _numReadFiles ];
  _fileIsActive = new bool[ _numReadFiles ];
  _fileToPmap = new ParamsMapping*[ _numReadFiles ];
  for (int n = 0; n < _numReadFiles; ++n){
    _numPairs[n] = hitReadsByFile[n]->size();
    _readQueries[n] = new ScoredSeqWithCarriers*[_numPairs[n]+1];
    _pairQueries[n] = new ScoredSeqWithCarriers*[_numPairs[n]+1];
    long pN = 0;
    for (vector<ScoredSeqWithCarriers*>::iterator it = hitReadsByFile[n]->begin(); it != hitReadsByFile[n]->end(); ++it){
      // I also want to re-map the read to the entire dataset so that its repetitive nature (ability
      // to map to many different contigs) can be fairly evaluated
      _readQueries[n][pN] = *it;
      _pairQueries[n][pN] = dynamic_cast<ScoredSeqWithCarriers*>( (*it)->getPairedEnd() );
      ++pN;
    }
    _numPairs[n] = pN;
    _readPairCounter[n] = 0;
    _totalNumQueries += pN * 2;
    if (pN == 0){ 
      _fileIsActive[n] = false;
    } else {
      _fileIsActive[n] = true;
      _numActiveFiles++;
    }
    _fileToPmap[n] = readFileArray[n]->getMapParams();
  }
}


MrsAddendumHelper::MrsAddendumHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf,
						  vector<ScoredSeqWithCarriers*>** hitReadsByFile,
						  vector<ScoredSeqWithCarriers*>** hitPairsByFile) :
  _rpf(rpf),
  _numReadFiles(numReadFiles) {
  _numActiveFiles = 0;
  _totalNumQueries = 0;

  int numFilesP1 = _numReadFiles + 1;
  _readFileArray = new ParameterizedReadFile*[ numFilesP1 ];
  for (int n = 0; n < _numReadFiles; ++n){ _readFileArray[n] = readFileArray[n]; }

  // make sets of the paired-ends of the mapped reads, THOSE will be mapped
  _readQueries = new ScoredSeqWithCarriers**[ numFilesP1 ];
  _pairQueries = new ScoredSeqWithCarriers**[ numFilesP1 ];
  _numPairs = new long[ numFilesP1 ];
  _readPairCounter = new long[ numFilesP1 ];
  _fileIsActive = new bool[ numFilesP1 ];
  _fileToPmap = new ParamsMapping*[ numFilesP1 ];
  for (int n = 0; n < _numReadFiles; ++n){
    if (hitReadsByFile[n]->size() != hitPairsByFile[n]->size()){
      throw AssemblyException::ArgError("EC::MRSAH constructor, the same number of reads and pairs must be input");
    }
    _numPairs[n] = hitReadsByFile[n]->size();
    _readQueries[n] = new ScoredSeqWithCarriers*[_numPairs[n]+1];
    _pairQueries[n] = new ScoredSeqWithCarriers*[_numPairs[n]+1];
    vector<ScoredSeqWithCarriers*>::iterator readIt = hitReadsByFile[n]->begin();
    vector<ScoredSeqWithCarriers*>::iterator pairIt = hitPairsByFile[n]->begin();
    for (long pN = 0; pN < _numPairs[n]; ++pN){
      _readQueries[n][pN] = *readIt;
      _pairQueries[n][pN] = *pairIt;
      ++readIt;
      ++pairIt;
    }
    _readPairCounter[n] = 0;
    _totalNumQueries += _numPairs[n] * 2;
    if (_numPairs[n] == 0){ 
      _fileIsActive[n] = false;
    } else {
      _fileIsActive[n] = true;
      _numActiveFiles++;
    }
    _fileToPmap[n] = readFileArray[n]->getMapParams();
  }
}



MrsAddendumHelper::~MrsAddendumHelper(){
  for (int n = 0; n < _numReadFiles; ++n){ 
    delete [] _readQueries[n];
    delete [] _pairQueries[n];
    delete _fileToPmap[n];
  }
  delete [] _readQueries;
  delete [] _pairQueries;
  delete [] _numPairs;
  delete [] _readPairCounter;
  delete [] _fileIsActive;
  delete [] _fileToPmap;
  delete [] _readFileArray;
}
int MrsAddendumHelper::numFilesTotal(){ return _numReadFiles; }
int MrsAddendumHelper::numFilesActive(){ return _numActiveFiles; }
long MrsAddendumHelper::totalReadCount(){ return _totalNumQueries; }
bool MrsAddendumHelper::isFileActive(int fileNum){ return _fileIsActive[fileNum]; }
ExtendJobMapper* MrsAddendumHelper::getNewEjm(set<ScoredSeqWithCarriers*>* parallelMapSet, int fileNum){
  return new ExtendJobMapperWithLinks( _readFileArray[fileNum]->ampliconSize(), parallelMapSet, _fileToPmap[fileNum]->fractId() );
}




void MrsAddendumHelper::getPairedReadArrays(ScoredSeqWithCarriers*** readArrays, ScoredSeqWithCarriers*** pairArrays,
					    long* numPairArray, int* fileNumArray, int numFiles){

  long numFilesP1 = numFiles + 1;

  // thread the reading from files by file
  #pragma omp parallel for schedule(dynamic)
  for (int fileIndex = 0; fileIndex < numFiles; ++fileIndex){
    int fileNum = fileNumArray[fileIndex];
    ParameterizedReadFile* readFile = _readFileArray[fileNum];
    long arrayTargetSize = _fileToPmap[fileNum]->stepSizeAddendum() / 2;
    if (arrayTargetSize == 0){ arrayTargetSize = 1; }
    ScoredSeqWithCarriers** readArray = new ScoredSeqWithCarriers*[ arrayTargetSize ];
    ScoredSeqWithCarriers** pairArray = new ScoredSeqWithCarriers*[ arrayTargetSize ];
    readArrays[fileIndex] = readArray;
    pairArrays[fileIndex] = pairArray;

    long pairCount = 0;
    while ( pairCount < arrayTargetSize and _readPairCounter[ fileNum ] < _numPairs[ fileNum ] ){
      readArray[pairCount] = _readQueries[ fileNum ][ _readPairCounter[ fileNum ] ];
      pairArray[pairCount] = _pairQueries[ fileNum ][ _readPairCounter[ fileNum ] ];
      readArray[pairCount]->bottomBuffer();
      pairArray[pairCount]->bottomBuffer();
      ++pairCount;
      ++_readPairCounter[ fileNum ];
    }
    numPairArray[fileIndex] = pairCount;
  }
}


long MrsAddendumHelper::decideMappability(bool* shouldReadsMap, bool* shouldPairsMap, int fileNum,
					  ScoredSeqWithCarriers** reads, ScoredSeqWithCarriers** pairs, long rpCount){
  long numUnmappable = 0;
  for (long n = 0; n < rpCount; ++n){
    if (_rpf->isReadOk(reads[n])){
      shouldReadsMap[n] = true;
    } else {
      shouldReadsMap[n] = false;
      ++numUnmappable;
    }
    if (_rpf->isReadOk(pairs[n])){
      shouldPairsMap[n] = true;
    } else {
      shouldPairsMap[n] = false;
      ++numUnmappable;
    }
    //shouldReadsMap[n] = _rpf->isReadOk(reads[n]);
    //shouldPairsMap[n] = _rpf->isReadOk(pairs[n]);
  }
  return numUnmappable;
}


ScoredSeqWithCarriers*** MrsAddendumHelper::getNewReadArrays(int* fileNumArray, int numFiles){
  int numFilesP1 = numFiles + 1;
  ScoredSeqWithCarriers*** readArrays = new ScoredSeqWithCarriers**[ numFilesP1 ];
  ScoredSeqWithCarriers*** pairArrays = new ScoredSeqWithCarriers**[ numFilesP1 ];
  long* pairNumArray = new long[ numFilesP1 ];

  getPairedReadArrays(readArrays, pairArrays, pairNumArray, fileNumArray, numFiles);

  ScoredSeqWithCarriers*** allSeqsToMapArrays = new ScoredSeqWithCarriers**[ numFilesP1 ];
  for (int n = 0; n < numFiles; ++n){
    bool* shouldReadsMap = new bool[ pairNumArray[n]+1 ];
    bool* shouldPairsMap = new bool[ pairNumArray[n]+1 ];
    decideMappability(shouldReadsMap, shouldPairsMap, n, readArrays[n], pairArrays[n], pairNumArray[n]);
    long readCounter = 0;
    ScoredSeqWithCarriers** seqsToMapArray = new ScoredSeqWithCarriers*[ pairNumArray[n] * 2 + 1 ];
    for (long readN = 0; readN < pairNumArray[n]; ++readN){
      if ( shouldReadsMap[readN] ){
	seqsToMapArray[ readCounter ] = readArrays[n][readN]; ++readCounter;
      } else { readArrays[n][readN]->deepDelete(); }
      if ( shouldPairsMap[readN] ){
	seqsToMapArray[ readCounter ] = pairArrays[n][readN]; ++readCounter;
      } else { pairArrays[n][readN]->deepDelete(); }
    }
    // end each array with a null
    seqsToMapArray[readCounter] = NULL;
    allSeqsToMapArrays[n] = seqsToMapArray;
    delete [] readArrays[n];
    delete [] pairArrays[n];
    delete [] shouldReadsMap;
    delete [] shouldPairsMap;
  }
  delete [] readArrays;
  delete [] pairArrays;
  delete [] pairNumArray;

  return allSeqsToMapArrays;
}


void MrsAddendumHelper::triggerFiles(){
  for (int n = 0; n < _numReadFiles; ++n){
    ReadFileCommunicator::bufferFill(_readFileArray[n], ReadFile::bufferThreaded);
  }
}


void MrsAddendumHelper::updateActiveFiles(){
  _numActiveFiles = 0;
  for (int fileN = 0; fileN < _numReadFiles; ++fileN){
    if (_fileIsActive[fileN]){
      if (_readPairCounter[fileN] < _numPairs[fileN]){ ++_numActiveFiles; }
      else { _fileIsActive[fileN] = false; }
    }
  }
}
bool MrsAddendumHelper::includeContigCarrier(ScoredSeqWithCarriers* carrier){ return true; }

#endif



