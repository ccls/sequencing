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

#ifndef ASSEMBLER_CPP
#define ASSEMBLER_CPP

#include "Assembler.h"


#include "Read.h"
#include "ScoredSeqSubseq.h"
#include "ScoredSeqNested.h"
#include "ScoredSeqMonoScore.h"
#include "OutputFileNull.h"
#include "OutputFileFasta.h"
#include "OutputFilePriceq.h"
#include "OutputFileMulti.h"
#include "AssemblerListenerNull.h"
#include "AssemblerListenerStd.h"
#include "AssemblerListenerVerbose.h"
#include "AssemblerListenerCarrier.h"
#include "RepeatDetector.h"
#include <sstream>
#include <omp.h>
#include <limits.h>

ParamsMinOverlap* Assembler::getDefaultParamsMinOverlap(){ return new ParamsMinOverlap(35,20,5); }
ParamsMinOverlap* Assembler::getNullParamsMinOverlap(){ return new ParamsMinOverlap(0,LONG_MAX,0); }


Assembler::Assembler(int totalCycles){
  _cycleInfo = new ParamsCycleManagement(totalCycles);
  _maxThreadsPerFile = 1;

  _miniMinOverlapParams = getDefaultParamsMinOverlap();
  _miniMinFractIdParams = ParamsMinFractId::getDefault();
  _metaMinFractIdParams = ParamsMinFractId::getDefault();
  _alScoreMatrix = AlignmentScoreMatrix::getDefault();

  _paramsDeBruijn = new ParamsDeBruijn();
  _stdOutListener = new AssemblerListenerStd();
  _fileListener = NULL;
  EcoFilter* defaultEcoFilterMinCount = EcoFilterMinCount::defaultFilter();
  _cycleInfo->addEcoFilter(defaultEcoFilterMinCount);
  delete defaultEcoFilterMinCount;
  _waitingTargetMode = NULL;
  _dynamicRepeatFilter = NULL;
}

void Assembler::setMiniMinOverlap(long minOvl){
  ParamsMinOverlap* pmoMiniHolder = new ParamsMinOverlap(minOvl,
							 _miniMinOverlapParams->minOvlThreshold(),
							 _miniMinOverlapParams->minOvlSlope());
  delete _miniMinOverlapParams;
  _miniMinOverlapParams = pmoMiniHolder;
}
void Assembler::setMiniMinOvlContigThreshold(long contigThreshold){
  ParamsMinOverlap* pmoMiniHolder = new ParamsMinOverlap(_miniMinOverlapParams->globalMinOverlap(),
							 contigThreshold,
							 _miniMinOverlapParams->minOvlSlope());
  delete _miniMinOverlapParams;
  _miniMinOverlapParams = pmoMiniHolder;
}
void Assembler::setMiniMinOvlScalingSlope(float slope){
  ParamsMinOverlap* pmoMiniHolder = new ParamsMinOverlap(_miniMinOverlapParams->globalMinOverlap(),
							 _miniMinOverlapParams->minOvlThreshold(),
							 slope);
  delete _miniMinOverlapParams;
  _miniMinOverlapParams = pmoMiniHolder;
}
void Assembler::setMiniMinFractId(float minFractId){
  ParamsMinFractId* pmfMiniHolder = new ParamsMinFractId(minFractId,
							 _miniMinFractIdParams->fractIdThreshold(),
							 _miniMinFractIdParams->denominatorLogStep());
  delete _miniMinFractIdParams;
  _miniMinFractIdParams = pmfMiniHolder;
}
void Assembler::setMiniFractIdContigThreshold(long contigThreshold){
  ParamsMinFractId* pmfMiniHolder = new ParamsMinFractId(_miniMinFractIdParams->globalFractId(),
							 contigThreshold,
							 _miniMinFractIdParams->denominatorLogStep());
  delete _miniMinFractIdParams;
  _miniMinFractIdParams = pmfMiniHolder;
}
void Assembler::setMiniFractIdScalingSlope(float slope){
  ParamsMinFractId* pmfMiniHolder = new ParamsMinFractId(_miniMinFractIdParams->globalFractId(),
							 _miniMinFractIdParams->fractIdThreshold(),
							 slope);
  delete _miniMinFractIdParams;
  _miniMinFractIdParams = pmfMiniHolder;
}
void Assembler::setMetaMinFractId(float minFractId){
  ParamsMinFractId* pmfMetaHolder = new ParamsMinFractId(minFractId,
							 _metaMinFractIdParams->fractIdThreshold(),
							 _metaMinFractIdParams->denominatorLogStep());
  delete _metaMinFractIdParams;
  _metaMinFractIdParams = pmfMetaHolder;
}
void Assembler::setMetaFractIdContigThreshold(long contigThreshold){
  ParamsMinFractId* pmfMetaHolder = new ParamsMinFractId(_metaMinFractIdParams->globalFractId(),
							 contigThreshold,
							 _metaMinFractIdParams->denominatorLogStep());
  delete _metaMinFractIdParams;
  _metaMinFractIdParams = pmfMetaHolder;
}
void Assembler::setMetaFractIdScalingSlope(float slope){
  ParamsMinFractId* pmfMetaHolder = new ParamsMinFractId(_metaMinFractIdParams->globalFractId(),
							 _metaMinFractIdParams->fractIdThreshold(),
							 slope);
  delete _metaMinFractIdParams;
  _metaMinFractIdParams = pmfMetaHolder;
}


void Assembler::setNucMatchScore(long score){
  _alScoreMatrix->_match = score;
}
void Assembler::setNucMismatchPenalty(long penalty){
  _alScoreMatrix->_mismatch = penalty;
}
void Assembler::setOpenGapPenalty(long penalty){
  _alScoreMatrix->_newGap = penalty;
}
void Assembler::setExtendGapPenalty(long penalty){
  _alScoreMatrix->_extendGap = penalty;
}


Assembler::~Assembler(){
  for (set<ScoredSeq*>::iterator it = _currentCycleContigs.begin(); it != _currentCycleContigs.end(); ++it){ (*it)->deepDelete(); }

  delete _miniMinOverlapParams;
  delete _miniMinFractIdParams;
  delete _metaMinFractIdParams;
  delete _alScoreMatrix;

  delete _cycleInfo;
  delete _paramsDeBruijn;
  for (set<ParameterizedReadFile*>::iterator it = _inputPairedReadFiles.begin(); it != _inputPairedReadFiles.end(); ++it){ delete (*it); }
  for (set<OutputFile*>::iterator it = _outfileSet.begin(); it != _outfileSet.end(); ++it){ delete *it; }
  delete _stdOutListener;
  delete _fileListener;
  for (vector<ReadPairFilter*>::iterator it = _waitingRpf.begin(); it != _waitingRpf.end(); ++it){ delete *it; }
  for (vector<ReadPairFilter*>::iterator it = _waitingIcf.begin(); it != _waitingIcf.end(); ++it){ delete *it; }
  for (set<ParameterizedInitialFile*>::iterator it = _inputInitialContigFiles.begin(); it != _inputInitialContigFiles.end(); ++it){
    delete *it;
  }
  delete _waitingTargetMode;
  if (_dynamicRepeatFilter != NULL){ delete _dynamicRepeatFilter; }
}



void Assembler::runSetup(){
  // consolidate the parameters
  _miniAlignParams = new ParamsAlignment(_miniMinOverlapParams, _miniMinFractIdParams, _alScoreMatrix);
  ParamsMinOverlap* dummyPmo = getNullParamsMinOverlap();
  _metaAlignParams = new ParamsAlignment(dummyPmo, _metaMinFractIdParams, _alScoreMatrix);
  delete dummyPmo;

  // initial contig files
  _initialContigFiles.insert( _inputInitialContigFiles.begin(), _inputInitialContigFiles.end() );
  _initialContigFiles.insert( _inputInitialContigFilesNoTarget.begin(), _inputInitialContigFilesNoTarget.end() );
  // and set up the target mode, including synchronization with the alignment score matrix
  if (_waitingTargetMode != NULL){
    _waitingTargetMode->addFiles( &_inputInitialContigFiles );
    _waitingTargetMode->setAlignmentScoreMatrix(_alScoreMatrix);
    _cycleInfo->addEcoFilter(_waitingTargetMode);
  }

  // read files
  ParameterizedReadFile* rfArray[ _inputPairedReadFiles.size() + 1 ];
  int prfCount = 0;
  for (set<ParameterizedReadFile*>::iterator it = _inputPairedReadFiles.begin();
       it != _inputPairedReadFiles.end(); ++it){
    rfArray[prfCount] = *it;
    prfCount++;
  }
  long rfSizes[ prfCount + 1 ];

  #pragma omp parallel for schedule(dynamic)
  for (int n = 0; n < prfCount; ++n){ rfSizes[n] = rfArray[n]->numReads(); }

  int numThreads = omp_get_max_threads();
  long readTotal = 0;
  for (int n = 0; n < prfCount; ++n){ readTotal += rfSizes[n]; }
  long threadDenom = readTotal / numThreads;
  if (threadDenom == 0){ threadDenom = 1; }

  #pragma omp parallel for schedule(dynamic)
  for (int n = 0; n < prfCount; ++n){
    int rfThreads = rfSizes[n] / threadDenom;
    if (rfThreads == 0){ rfThreads = 1; }
    else if (rfThreads > _maxThreadsPerFile){ rfThreads = _maxThreadsPerFile; }
    set<ParameterizedReadFile*> splitFiles;
    rfArray[n]->splitFile(&splitFiles, rfThreads);
    #pragma omp critical (Assembler_addRfsToSet)
    { _pairedReadFiles.insert( splitFiles.begin(), splitFiles.end() ); }
  }

  // output file(s)
  if (_outfileSet.size() == 0){ _outfile = new OutputFileNull(); }
  else {
    OutputFile** outfileArray = new OutputFile*[ _outfileSet.size() ];
    int arrayIndex = 0;
    for (set<OutputFile*>::iterator it = _outfileSet.begin(); it != _outfileSet.end(); ++it){
      outfileArray[arrayIndex] = *it;
      arrayIndex++;
    }
    _outfile = new OutputFileMulti(outfileArray, arrayIndex);
    delete [] outfileArray;
  }

  // the listener
  if (_fileListener == NULL){ _listener = _stdOutListener; }
  else {
    AssemblerListener* listenerArray[] = {_stdOutListener, _fileListener};
    _listener = new AssemblerListenerCarrier(listenerArray, 2);
  }
}


void Assembler::runCleanup(){
  delete _miniAlignParams;
  delete _metaAlignParams;

  delete _outfile;
  if (_fileListener != NULL){ delete _listener; }
  for (set<ParameterizedReadFile*>::iterator it = _pairedReadFiles.begin(); it != _pairedReadFiles.end(); ++it){
    if (_inputPairedReadFiles.count(*it) == 0){ delete *it; }
  }
  _pairedReadFiles.clear();
  _initialContigFiles.clear();
}


void Assembler::addInitialFile(string filename, int numSteps, int cyclesPerStep, float countFactor){
  ReadFile* file = ReadFile::makeBasicReadFile(filename, countFactor);
  ParameterizedInitialFile* pif = new ParameterizedInitialFile( file, numSteps, cyclesPerStep );
  delete file;
  _inputInitialContigFiles.insert(pif);
}
void Assembler::addInitialFile(long numSeqsToUse, string filename, int numSteps, int cyclesPerStep, float countFactor){
  ReadFile* file = ReadFile::makeBasicReadFile(filename, countFactor);
  ParameterizedInitialFile* pif = new ParameterizedInitialFile( numSeqsToUse, file, numSteps, cyclesPerStep );
  delete file;
  _inputInitialContigFiles.insert(pif);
}
void Assembler::addInitialFileNoTarget(string filename, int numSteps, int cyclesPerStep, float countFactor){
  ReadFile* file = ReadFile::makeBasicReadFile(filename, countFactor);
  ParameterizedInitialFile* pif = new ParameterizedInitialFile( file, numSteps, cyclesPerStep );
  delete file;
  _inputInitialContigFilesNoTarget.insert(pif);
}
void Assembler::addInitialFileNoTarget(long numSeqsToUse, string filename, int numSteps, int cyclesPerStep, float countFactor){
  ReadFile* file = ReadFile::makeBasicReadFile(filename, countFactor);
  ParameterizedInitialFile* pif = new ParameterizedInitialFile( numSeqsToUse, file, numSteps, cyclesPerStep );
  delete file;
  _inputInitialContigFilesNoTarget.insert(pif);
}



void Assembler::addReadFile(string filename, long ampSize){
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFile(string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int numCycles, string filename, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int numCycles, string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}


void Assembler::addReadFile(string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFile(string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow paired-end read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makePairedReadFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}


// MATE PAIRS BEGIN

void Assembler::addMatePairFile(string filename, long ampSize){
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFile(string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int numCycles, string filename, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int numCycles, string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filename, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}


void Assembler::addMatePairFile(string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFile(string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize, float fractId){
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("Assembler does not allow mate-pair read files to be the same file.");
  }
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeMatePairFile(filenameA, filenameB, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}

// MATE PAIRS END


// FALSE PAIRS BEGIN

void Assembler::addFalsePairFile(string filename, long readLength, long ampSize){
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addFalsePairFile(string filename, long readLength, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams ) );
  delete rf;
  delete mapParams;
}
void Assembler::addFalsePairFileLimitUse(int startCycle, int numCycles, string filename, long readLength, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles ) );
  delete rf;
  delete mapParams;
}
void Assembler::addFalsePairFileLimitUse(int startCycle, int numCycles, string filename, long readLength, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, numCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addFalsePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long readLength, long ampSize){
  ParamsMapping* mapParams = ParamsMapping::getDefault();
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}
void Assembler::addFalsePairFileLimitUse(int startCycle, int onCycles, int skipCycles,
					 string filename, long readLength, long ampSize, float fractId){
  if (fractId < .25 or fractId > 1){ throw AssemblyException::ArgError("% ID required for mapping is out of range 25-100"); }
  ParamsMapping* defMapParams = ParamsMapping::getDefault();
  ParamsMapping* mapParams = new ParamsMapping(defMapParams->stepSize(), fractId );
  delete defMapParams;
  ReadFile* rf = ReadFile::makeFalsePairFile(filename, readLength, ampSize);
  _inputPairedReadFiles.insert( new ParameterizedReadFile( rf, mapParams, startCycle, onCycles, skipCycles) );
  delete rf;
  delete mapParams;
}

// FALSE PAIRS END


void Assembler::addOutputFile(string filename){
  _outfileSet.insert( OutputFile::makeOutputFile(filename) );
}


// FILTER READS

void Assembler::addBadSequenceFilter(string badSeqFile, float fractId){
  ReadFile* rf = ReadFile::makeBasicReadFile(badSeqFile);
  set<ScoredSeq*> badSeqs;
  rf->open();
  while ( rf->hasRead() ){ badSeqs.insert( rf->getRead() ); }
  rf->close();
  _waitingRpf.push_back( new ReadPairFilterAvoidSeqs(&badSeqs, fractId) );
  for (set<ScoredSeq*>::iterator badIt = badSeqs.begin(); badIt != badSeqs.end(); ++badIt){ delete *badIt; }
  delete rf;
}
void Assembler::addHomopolymerFilter(long maxLength){
  _waitingRpf.push_back( new ReadPairFilterHomoPolymer(maxLength) );
}
void Assembler::addDinucRepeatFilter(long maxLength){
  _waitingRpf.push_back( new ReadPairFilterDinucleotide(maxLength) );
}
void Assembler::addReadQualityFilter(float minFractGood, float minProbCorrect){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect) );
}
void Assembler::addReadQualityFilter(float minFractGood, float minProbCorrect, int numSkipCycles, int numRunCycles){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect, numSkipCycles, numRunCycles) );
}
void Assembler::addReadCalledBasesFilter(float minFractGood){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterUncalledBases(maxFractBad) );
}
void Assembler::addReadCalledBasesFilter(float minFractGood, int numSkipCycles, int numRunCycles){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterUncalledBases(maxFractBad, numSkipCycles, numRunCycles) );
}



// FILTER INITIAL CONTIGS

void Assembler::icBadSequenceFilter(string badSeqFile, float fractId){
  icBadSequenceFilterHelper(badSeqFile, fractId, false, 0);
}
void Assembler::icBadSequenceFilter(string badSeqFile, float fractId, long maxSeqLength){
  icBadSequenceFilterHelper(badSeqFile, fractId, true, maxSeqLength);
}
void Assembler::icBadSequenceFilterHelper(string badSeqFile, float fractId, bool useMaxLength, long maxSeqLength){
  ReadFile* rf = ReadFile::makeBasicReadFile(badSeqFile);
  set<ScoredSeq*> badSeqs;
  rf->open();
  while ( rf->hasRead() ){ badSeqs.insert( rf->getRead() ); }
  rf->close();
  if (useMaxLength){ _waitingIcf.push_back( new ReadPairFilterAvoidSeqs(&badSeqs, fractId, maxSeqLength) ); }
  else { _waitingIcf.push_back( new ReadPairFilterAvoidSeqs(&badSeqs, fractId) ); }
  for (set<ScoredSeq*>::iterator badIt = badSeqs.begin(); badIt != badSeqs.end(); ++badIt){ delete *badIt; }
  delete rf;
}
void Assembler::icHomopolymerFilter(long maxLength){
  _waitingIcf.push_back( new ReadPairFilterHomoPolymer(maxLength) );
}
void Assembler::icHomopolymerFilter(long maxLength, long maxSeqLength){
  _waitingIcf.push_back( new ReadPairFilterHomoPolymer(maxLength, maxSeqLength) );
}
void Assembler::icDinucRepeatFilter(long maxLength){
  _waitingIcf.push_back( new ReadPairFilterDinucleotide(maxLength) );
}
void Assembler::icDinucRepeatFilter(long maxLength, long maxSeqLength){
  _waitingIcf.push_back( new ReadPairFilterDinucleotide(maxLength, maxSeqLength) );
}
void Assembler::icReadQualityFilter(float minFractGood, float minProbCorrect){
  float maxFractBad = 1.0 - minFractGood;
  _waitingIcf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect) );
}
void Assembler::icReadQualityFilter(float minFractGood, float minProbCorrect, long maxSeqLength){
  float maxFractBad = 1.0 - minFractGood;
  _waitingIcf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect, maxSeqLength) );
}
void Assembler::icReadCalledBasesFilter(float minFractGood){
  float maxFractBad = 1.0 - minFractGood;
  _waitingIcf.push_back( new ReadPairFilterUncalledBases(maxFractBad) );
}
void Assembler::icReadCalledBasesFilter(float minFractGood, long maxSeqLength){
  float maxFractBad = 1.0 - minFractGood;
  _waitingIcf.push_back( new ReadPairFilterUncalledBases(maxFractBad, maxSeqLength) );
}


// FILTER NEW CONTIGS

void Assembler::addLengthFilter(long minLength, int cyclesToSkip){
  _cycleInfo->addSizeFilter(minLength,cyclesToSkip);
}
void Assembler::addTargetMode(float minFractId, int cyclesToSkip, bool fullFile){
  delete _waitingTargetMode;
  _waitingTargetMode = new EcoFilterInitialContigMatch(minFractId,cyclesToSkip,fullFile);
}
void Assembler::addTargetMode(float minFractId, int cyclesToSkip,
			      int numFilterCycles, int numSkipCycles, bool fullFile){
  delete _waitingTargetMode;
  _waitingTargetMode = new EcoFilterInitialContigMatch(minFractId,cyclesToSkip,numFilterCycles,numSkipCycles,fullFile);
}



void Assembler::verboseLogFile(string filename){
  delete _fileListener;
  _fileListener = new AssemblerListenerVerbose( filename );
}
void Assembler::makeLogNull(){
  delete _stdOutListener;
  _stdOutListener = new AssemblerListenerNull();
}
void Assembler::makeLogVerbose(){
  delete _stdOutListener;
  _stdOutListener = new AssemblerListenerVerbose();
}


void Assembler::setMaxThreadsPerFile(int maxThreadsPerFile){
  _maxThreadsPerFile = maxThreadsPerFile;
}

void Assembler::setDeBruijnKmerSize(int kmerSize){
  _paramsDeBruijn->setKmerSize(kmerSize);
}
void Assembler::setDeBruijnMaxSeqLength(long maxSeqLength){
  _paramsDeBruijn->setMaxSeqLength(maxSeqLength);
}
void Assembler::setDeBruijnMinSeqNum(long minSeqNum){
  _paramsDeBruijn->setMinNumSeqs(minSeqNum);
}



// METHODS THAT FORWARD TO PARAMSCYCLEMANAGEMENT

long Assembler::getLinkMax(){
  return _cycleInfo->maxContigLinkages();
}
void Assembler::setLinkMax(long maxNumLinks){
  _cycleInfo->setMaxContigLinkages(maxNumLinks);
}

void Assembler::getResetCycles(set<int>* resetCycles){
  int totalCycles = _cycleInfo->totalCycles();
  for (int n = 1; n <= totalCycles; ++n){
    if (_cycleInfo->shouldContigsReset(n)){ resetCycles->insert(n); }
  }
}
void Assembler::addResetCycle(int cycleNum){
  _cycleInfo->addContigReset(cycleNum);
}

int Assembler::getOutputInterval(){
  return _cycleInfo->getOutputInterval();
}
void Assembler::setOutputInterval(int interval){
  _cycleInfo->setOutputInterval(interval);
}

void Assembler::findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
				    long minRepeatSize, float matchFractId){
  _cycleInfo->addRepeatDetection(cycleNum, numStdDevs, minFoldIncrease, ParamsCycleManagement::BEFORECYCLE,
				 minRepeatSize, matchFractId);
}
void Assembler::findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
				    long minRepeatSize, float matchFractId, string outfileName){
  _cycleInfo->addRepeatDetection(cycleNum, numStdDevs, minFoldIncrease, ParamsCycleManagement::BEFORECYCLE,
				 minRepeatSize, matchFractId, outfileName);
}
void Assembler::findRepsAfterCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
				   long minRepeatSize, float matchFractId){
  _cycleInfo->addRepeatDetection(cycleNum, numStdDevs, minFoldIncrease, ParamsCycleManagement::AFTERCYCLE,
				 minRepeatSize, matchFractId);
}
void Assembler::findRepsAfterCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
				   long minRepeatSize, float matchFractId, string outfileName){
  _cycleInfo->addRepeatDetection(cycleNum, numStdDevs, minFoldIncrease, ParamsCycleManagement::AFTERCYCLE,
				 minRepeatSize, matchFractId, outfileName);
}


void Assembler::addTrim(int cycleNum, float coverageLevel){
  _cycleInfo->addTrim(cycleNum, coverageLevel);
}
void Assembler::addTrim(int cycleNum, float coverageLevel, long minLength){
  _cycleInfo->addTrim(cycleNum, coverageLevel, minLength);
}
void Assembler::addBasalTrim(int numPriorCycles, float coverageLevel){
  _cycleInfo->addBasalTrim(numPriorCycles, coverageLevel);
}
void Assembler::addBasalTrim(int numPriorCycles, float coverageLevel, long minLength){
  _cycleInfo->addBasalTrim(numPriorCycles, coverageLevel, minLength);
}
void Assembler::addInitialTrim(float coverageLevel){
  _cycleInfo->addInitialTrim(coverageLevel);
}
void Assembler::addInitialTrim(float coverageLevel, long minLength){
  _cycleInfo->addInitialTrim(coverageLevel, minLength);
}


void Assembler::runAssembly(){
  runSetup();
  _listener->beginAssembly("BEGIN ASSEMBLY");
  for (int n = 0; n < _cycleInfo->totalCycles(); ++n){
    _listener->beginCycle(n+1);
    generateInitialContigs(n);
    extendOneCycle(n);

    // make the number into an addendum string and write the file
    if ( _cycleInfo->outputThisCycle(n+1) ){
      stringstream addendum;
      addendum << "cycle" << n+1;
      writeOutfile( addendum.str() );
    }

    _listener->endCycle(&_currentCycleContigs);
  }
  // delete the contigs from the prior cycle in the end
  for (set<ScoredSeq*>::iterator it = _priorCycleContigs.begin(); it != _priorCycleContigs.end(); ++it){
    (*it)->deepDelete();
  }
  _listener->endAssembly("END ASSEMBLY");
  runCleanup();
}

void Assembler::writeOutfile(string addendum){
  if (_outfile->getName(addendum) != ""){
    _listener->writeOutfile( _outfile->getName(addendum) );
  }
  _outfile->open(addendum);
  ScoredSeq** sortedContigs = OutputFile::makeSortedContigArray( &_currentCycleContigs );
  _outfile->writeContigs(sortedContigs, &_currentCycleLeftovers);
  delete [] sortedContigs;
  _outfile->close();
}



void Assembler::getContigs(set<ScoredSeq*>* outputContigs){
  for (set<ScoredSeq*>::iterator it = _currentCycleContigs.begin(); it != _currentCycleContigs.end(); ++it){
    outputContigs->insert( (*it)->shallowCopy() );
  }
}



void Assembler::findLeftoverContigs(set<ScoredSeq*>* current, set<ScoredSeq*>* prior, set<ScoredSeq*>* leftover){
  // collect the sequences of the contigs being input into this assembly round
  set<string> seqsFromPrior;
  for (set<ScoredSeq*>::iterator it = prior->begin(); it != prior->end(); ++it){
    char* tempSeq = (*it)->getSeq('+');
    seqsFromPrior.insert( string(tempSeq) );
    delete [] tempSeq;
  }
  // now identify the contigs from the current cycle that didn't make progress from the last
  for (set<ScoredSeq*>::iterator it = current->begin(); it != current->end(); ++it){
    char* tempPlus = (*it)->getSeq('+');
    char* tempMinus =  (*it)->getSeq('-');
    if (seqsFromPrior.count( tempPlus ) > 0 or 
	seqsFromPrior.count( tempMinus ) > 0 ){
      leftover->insert( (*it) );
    }
    delete [] tempPlus;
    delete [] tempMinus;
  }
}


void Assembler::generateInitialContigs(int cycleNum){
  set<ScoredSeq*> initialSeqs;

  // put the current contigs in with the new initial seqs
  initialSeqs.insert( _currentCycleContigs.begin(), _currentCycleContigs.end() );
  _currentCycleContigs.clear();

  // the new set of input sequences will need to be filtered
  long initialCount = initialSeqs.size();
  set<ScoredSeq*> unfilteredInputSeqs;
  for (set<ParameterizedInitialFile*>::iterator rfIt = _initialContigFiles.begin(); rfIt != _initialContigFiles.end(); ++rfIt){
    (*rfIt)->getContigs(&unfilteredInputSeqs, cycleNum);
  }
  _listener->collectedInitialContigs(&unfilteredInputSeqs);

  // this is the set of filtered sequences; iterator for fast addition since
  // the sequences are already sorted for sets
  set<ScoredSeq*> filteredInputSeqs;

  // only access the filtering code and print filter output if there are actually filters
  if (_waitingIcf.size() == 0){
    filteredInputSeqs.insert(unfilteredInputSeqs.begin(), unfilteredInputSeqs.end());
  } else {
    set<ScoredSeq*>::iterator fisIt = filteredInputSeqs.begin();
    for (set<ScoredSeq*>::iterator ufisIt = unfilteredInputSeqs.begin(); ufisIt != unfilteredInputSeqs.end(); ++ ufisIt){
      bool seqIsOk = true;
      vector<ReadPairFilter*>::iterator wfIt = _waitingIcf.begin();
      while (seqIsOk and wfIt != _waitingIcf.end()){
	if ( (*wfIt)->applyToLength( (*ufisIt)->size() ) ){
	  seqIsOk = (*wfIt)->isReadOk( *ufisIt );
	}
	++wfIt;
      }
      if (seqIsOk){ fisIt = filteredInputSeqs.insert(fisIt, *ufisIt); }
      else { (*ufisIt)->deepDelete(); }
    }
    _listener->filteredInitialContigs(&filteredInputSeqs);
  }

  // now add the filtered sequences to the initial sequence set
  initialSeqs.insert(filteredInputSeqs.begin(), filteredInputSeqs.end());

  // the input contigs
  _cycleInfo->executeInitialTrim(&initialSeqs);

  if ( initialCount == initialSeqs.size() ){
    _currentCycleContigs.insert( initialSeqs.begin(), initialSeqs.end() );
  } else {
    _listener->beginAddInitialContigs(&initialSeqs);

    // condense (assemble) the initial contigs (assembly is double-stranded)
    AssemblyJobFactory * contigFactory = new AssemblyJobFactory(_metaAlignParams,_paramsDeBruijn,AssemblyJob::DOUBLESTRANDED,_listener);
    AssemblyJob * ajr2 = contigFactory->redundancyJob( &initialSeqs );
    delete contigFactory;
    ajr2->runJob(AssemblyJob::THREADED);
    // add the output of the initial assembly to the collection
    // currentCycleContigs is the initial input; priorCycleContigs is initially empty
    ajr2->allAssembledSeqs( &_currentCycleContigs );

    set<ScoredSeq*> usedSeqs;
    ajr2->discardedInputSeqs( &usedSeqs );
    for (set<ScoredSeq*>::iterator it = usedSeqs.begin(); it != usedSeqs.end(); ++it){ (*it)->deepDelete(); }
    delete ajr2;

    _listener->endAddInitialContigs(&_currentCycleContigs);
  }
}


void Assembler::extendOneCycle(int cycleNum){

  // since the params consider cycles as indexed from 1
  int cycleNumP1 = cycleNum + 1;

  // identify the contigs that have not changed since the previous cycle (unless this is a "reset" cycle)
  set<ScoredSeq*> leftoverContigs;
  if ( ! _cycleInfo->shouldContigsReset(cycleNumP1) ){
    findLeftoverContigs(&_currentCycleContigs, &_priorCycleContigs, &leftoverContigs);
  }
  for (set<ScoredSeq*>::iterator it = _priorCycleContigs.begin(); it != _priorCycleContigs.end(); ++it){ (*it)->deepDelete(); }
  _priorCycleContigs.clear();

  ostringstream outMessage;
  outMessage << "num contigs leftover: " << leftoverContigs.size();
  _listener->message(outMessage.str());

  // isolate just the read files that are appropriate for this cycle
  set<ParameterizedReadFile*> thisCycleReadFiles;
  for (set<ParameterizedReadFile*>::iterator it = _pairedReadFiles.begin(); it != _pairedReadFiles.end(); ++it){
    if ( (*it)->useThisCycle(cycleNum) ){ thisCycleReadFiles.insert( (*it) ); }
  }

  // seek repeats (maybe)
  detectRepeatsDynamically(cycleNumP1, ParamsCycleManagement::BEFORECYCLE);

  // filter the ReadPairFilters for those being used this cycle, and sort by
  // whether or not the test will knock out both reads of a pair (in which case
  // it will be applied prior to the first mapping step) or just one read of a
  // pair (in which case it will wait until the second mapping step to be applied);
  // remember to consider the _dynamicRepeatFilter
  ReadPairFilter** rpfSingleUsed = new ReadPairFilter*[ _waitingRpf.size() + 2 ];
  ReadPairFilter** rpfDoubleUsed = new ReadPairFilter*[ _waitingRpf.size() + 2 ];
  int numRpfSingleUsed = 0;
  int numRpfDoubleUsed = 0;
  // the repeat filter will be a single filter
  if (_dynamicRepeatFilter != NULL){
    rpfSingleUsed[numRpfSingleUsed] = _dynamicRepeatFilter;
    ++numRpfSingleUsed;
  }
  for (vector<ReadPairFilter*>::iterator rpfIt = _waitingRpf.begin(); rpfIt != _waitingRpf.end(); ++rpfIt){
    if ( (*rpfIt)->useThisCycle(cycleNum) ){
      if ( (*rpfIt)->testKillsPair() ){
	rpfDoubleUsed[numRpfDoubleUsed] = *rpfIt;
	++numRpfDoubleUsed;
      } else {
	rpfSingleUsed[numRpfSingleUsed] = *rpfIt;
	++numRpfSingleUsed;
      }
    }
  }

  // create single and double ReadPairFilters that include those tests being used this cycle
  ReadPairFilter* rpfSingle;
  if (numRpfSingleUsed == 0){ rpfSingle = new ReadPairFilterNull(); }
  else if (numRpfSingleUsed == 1){ rpfSingle = rpfSingleUsed[0]; }
  else { rpfSingle = new ReadPairFilterMulti(rpfSingleUsed, numRpfSingleUsed); }
  delete [] rpfSingleUsed;

  ReadPairFilter* rpfDouble;
  if (numRpfDoubleUsed == 0){ rpfDouble = new ReadPairFilterNull(); }
  else if (numRpfDoubleUsed == 1){ rpfDouble = rpfDoubleUsed[0]; }
  else { rpfDouble = new ReadPairFilterMulti(rpfDoubleUsed, numRpfDoubleUsed); }
  delete [] rpfDoubleUsed;


  ExtendCycle* ec = new ExtendCycle(cycleNum, &thisCycleReadFiles, rpfSingle, rpfDouble, &_currentCycleContigs, &leftoverContigs, _listener);

  // this was either a null or a multi-carrier - unless there was just one
  if (numRpfSingleUsed != 1){ delete rpfSingle; }
  if (numRpfDoubleUsed != 1){ delete rpfDouble; }

  _priorCycleContigs.insert(_currentCycleContigs.begin(), _currentCycleContigs.end());
  _currentCycleContigs.clear();

  // run the cycle and collect, filter, and trim the output
  EcoFilter* sizeFilter = _cycleInfo->getSizeFilter();
  EcoFilter* fullFilter = _cycleInfo->getFullEcoFilter();
  ec->runCycle(_miniAlignParams, _metaAlignParams, _paramsDeBruijn, _cycleInfo->maxContigLinkages(), sizeFilter, fullFilter);
  delete sizeFilter;
  delete fullFilter;
  ec->getContigs(&_currentCycleContigs);
  delete ec;

  _cycleInfo->executeTrim( &_currentCycleContigs, cycleNumP1 );

  // identify those contigs that did not change during this cycle
  _currentCycleLeftovers.clear();
  findLeftoverContigs(&_currentCycleContigs, &_priorCycleContigs, &_currentCycleLeftovers);

  // seek repeats (maybe again)
  detectRepeatsDynamically(cycleNumP1, ParamsCycleManagement::AFTERCYCLE);
}


void Assembler::detectRepeatsDynamically(int cycleNumP1, ParamsCycleManagement::WhenInCycle whenInCycle){
  if ( _cycleInfo->seekDynamicRepeats(cycleNumP1, whenInCycle) ){

    _listener->beginRepeatDetection();

    RepeatDetector* repDetector = new RepeatDetector(_cycleInfo->repeatMinSize(cycleNumP1, whenInCycle), _metaAlignParams);
    ScoredSeq** repeats = repDetector->findRepeats(&_currentCycleContigs,
						   _cycleInfo->repeatStandDevs(cycleNumP1, whenInCycle),
						   _cycleInfo->repeatMinFoldIncrease(cycleNumP1, whenInCycle),
						   RepeatDetector::THREADED);

    delete repDetector;
    set<ScoredSeq*> repSet;
    long repNum = 0;
    long totalRepLength = 0;
    while (repeats[repNum] != NULL){
      repSet.insert( repeats[repNum] );
      totalRepLength += repeats[repNum]->size();
      ++repNum;
    }
    _listener->endRepeatDetection(repNum,totalRepLength);

    // out with the old, in with the new
    if (_dynamicRepeatFilter != NULL){ delete _dynamicRepeatFilter; }
    _dynamicRepeatFilter = new ReadPairFilterAvoidSeqs(&repSet, _cycleInfo->repeatMatchFractId(cycleNumP1, whenInCycle));

    // write a file?
    if ( _cycleInfo->makeRepeatOutfile(cycleNumP1, whenInCycle) ){
      OutputFile* outfile = OutputFile::makeOutputFile( _cycleInfo->repeatOutfileName(cycleNumP1, whenInCycle) );
      outfile->open();
      outfile->writeContigs(repeats);
      outfile->close();
      delete outfile;
    }

    // delete the original copies of the repeat sequences
    for (set<ScoredSeq*>::iterator badIt = repSet.begin(); badIt != repSet.end(); ++badIt){ delete *badIt; }
    delete [] repeats;
  }
}


long Assembler::getTotalSize(){
  long count = 0;
  for (set<ScoredSeq*>::iterator itS = _currentCycleContigs.begin(); itS != _currentCycleContigs.end(); ++itS){
    count += (*itS)->size();
  }
  return count;
}

long Assembler::getMaxSize(){
  long maxSize = 0;
  for (set<ScoredSeq*>::iterator itS = _currentCycleContigs.begin(); itS != _currentCycleContigs.end(); ++itS){
    if ( (*itS)->size() > maxSize ){ maxSize = (*itS)->size(); }
  }
  return maxSize;
}


long Assembler::getN50(){
  vector<long> contigLengths;
  for (set<ScoredSeq*>::iterator seqIt = _currentCycleContigs.begin(); seqIt != _currentCycleContigs.end(); ++seqIt){
    contigLengths.push_back( (*seqIt)->size() );
  }
  sort( contigLengths.begin(), contigLengths.end() );
  long halfSize = getTotalSize() / 2;

  long n50 = 0;
  vector<long>::iterator sizeIt = contigLengths.begin();
  while ( halfSize > 0 and sizeIt != contigLengths.end() ){
    n50 = (*sizeIt); // current n50; will only be updated if half the total size is not reached
    halfSize -= (*sizeIt);
    ++sizeIt;
  }
  return n50;
}
bool reverseN50Helper(long i,long j) { return (i>j); }


#endif



/*
THE AUTHOR: generated by www.degraeve.com/img2txt.php
::::;::::;::::i:,,,tiiiitt;;,;,,,;,,,,,,
::::;::::;:,i:i;i;;jfGLfji:,,;,,,;,,,,,,
::::;::::,:;,jGGfLGLLGLffLtit;,,,;,,,,,,
::::,::::;;jLGGGEEEDGGGLjfLLt;::,;,,,,,,
::::,::::ijLDGDGLKWKEDGGLffLt;,::;:,,,,,
::::,:::ifDGDDGGDEKWKEDLLLGGft:::;::,,,:
::::,::;tiLDDEDDDDDEKKKKEDLLLLi,:;,::,,,
:..:,::jitEDKKDGLLLLGGDEEDDfGLf,:;::,,,,
:::::::L,LDEWDLLffffLLGDDDDfLGLi:,:::,,:
..::::;jiDDKKLLLfffffffLGDEGLLLf:,:.::::
...:,.;:tfEWDGLLffffjjjjjLEDLDLf::...:.:
:...:..;iEEKGGLLfffjfjtttjLDDDLf::...::.
....:t:.LEEKLGLLfffjjjttiitDEEGLL:....:.
:...:LGLEDDDGLLffLLffjttiiijEKGj;:......
:...:jKKKKWELDDDKEDLfjjjjjiiLEEj::......
..:.tiEEKKWDWWKKWKKDffDKDGj;DEDG;:....:.
..:;iGKEKKKLWKEKKWKELEWEEGG#EEDEKtj;,;;,
ii;iiLEEDKDGKDEEEKEDjEEEGfGt#EDDEDi,,,,,
i;;iiiKKKKDLGGGDEEELtWLEEGjtGGGDDL,,,,,,
iiiiiiGfGEDLLLLGDDELiGfLLftfLGEDLfi;;iii
iiiiiiitKDDLLGEEEEGfiiGLLji;LGEGijiiiiii
iiiitiitLfDLGGDDDLGLiitELfi;jDGiijiiiiii
iiiiiiiittGLGDEEEEEDfLfDDLtiGLiiijiiiitt
tttttiiiitEGGGEDDEDDGffLGDtiGLi;itiiiiii
fffffffLGfEDGGDDWEEGGLLKGLjLjfjiitiiiiii
ttttttttitLEGGLGDEGLjDftftiDjiiiitiiiiii
#########WL#GGGGGDDGGjitftDKjiiiitiiiitt
WWWWWWWWKKfWLGGGGGGGLjtttftiiiiiijtttttt
WWWWWDWWWfEKLGGDLGLLfjtittiiitiiijtttttt
WWWWWKKKfGEDGGGDDGLLfttjitiiitiiitiiiitt
WWWWWKLfLfLDGGGDEDGGLjfjjiittttttftttttj
WWWLLGGLLLGEGGGGDDDDDLjtKG;;tttttftttjjj
jjfjjLLfLDKKGGGGGDGGjtttWD;jt;tttftjjjjj
ttjLLLffLLDKDGDGGGGfjtttKf;Lfi,;tjjjjjjj
tjDLDfGLtjGEDGDGDDDfjjjGWL;;tL;,;:tjjjjf
fGLDLLtjjtfDDDDDDDGLjfLEED,iGfLt,,;::;tj
jfLDfjGGGLtLEDGDDGGLLGELfGt;ijGD;ti;;,,j
LjfDfLGKWWWjfDGGDDDDGDjGjft;;fDLij,ti,:D
fDLLLDEKKWKGLjGLGGDGGtLGjfD:itfjfjt,j;tD
LLDGGEDWDEDGDjjtijLitfDjfjKKitDDjff,fiiL
DjLDfDGGGGEGGGjjttjtGfffGDKWKfGDiGL;tiji
jDLjDfGLLLGDGEftiiijtfLLfGEEWttLfLLf;jit
fLEGLjGjffGGGGEftitjfjGjtjjLKiLGttGftjit
GjfDjffjjfLLLLjDfjfffGLti;jjtiGL;fELfijj
LfjGfffjjfLLfffjDLLfGGjji;ttfLLLtEDLLjGj
fLGGffjjjjfGjjjjfGjLGLttiijtttjGjjfDiLLi
LffLjfjjjffGjjjfLGfGDGjtttttjtjjtfLLfff;
LfffjjjjfffGjjffLtfEDDjtfjjtjjffjffjtfjt
GLffjjfjfffGjjfftiGDDGLfffjtjjfLjfGGjjLj
*/
