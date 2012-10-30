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

#ifndef PARAMETERIZEDREADFILE_CPP
#define PARAMETERIZEDREADFILE_CPP

#include "ParameterizedReadFile.h"
#include "AssemblyException.h"
#include "ReadFileCommunicator.h"
#include <iostream>
using namespace::std;

// constructors

ParameterizedReadFile::ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm) :
  _initialCycle(0),
  _isFinite(false),
  _totalCycles(0),
  _isCyclic(false),
  _cyclesOn(0),
  _cyclesSkipped(0)
{
  _readFile = readFile->copy();
  _pm = new ParamsMapping(pm);
}
ParameterizedReadFile::ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle) :
  _initialCycle(initialCycle),
  _isFinite(false),
  _totalCycles(0),
  _isCyclic(false),
  _cyclesOn(0),
  _cyclesSkipped(0)
{
  _readFile = readFile->copy();
  _pm = new ParamsMapping(pm);
}
ParameterizedReadFile::ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle, int numCycles) :
  _initialCycle(initialCycle),
  _isFinite(true),
  _totalCycles(numCycles),
  _isCyclic(false),
  _cyclesOn(0),
  _cyclesSkipped(0)
{
  _readFile = readFile->copy();
  _pm = new ParamsMapping(pm);
}
ParameterizedReadFile::ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle, int cyclesOn, int cyclesSkipped) :
  _initialCycle(initialCycle),
  _isFinite(false),
  _totalCycles(0),
  _isCyclic(true),
  _cyclesOn(cyclesOn),
  _cyclesSkipped(cyclesSkipped)
{
  _readFile = readFile->copy();
  _pm = new ParamsMapping(pm);
}
ParameterizedReadFile::ParameterizedReadFile(ParameterizedReadFile* readFile, ReadFile* replacementFile){
  _readFile = replacementFile->copy();
  _initialCycle = readFile->_initialCycle;
  _isFinite = readFile->_isFinite;
  _totalCycles = readFile->_totalCycles;
  _isCyclic = readFile->_isCyclic;
  _cyclesOn = readFile->_cyclesOn;
  _cyclesSkipped = readFile->_cyclesSkipped;
  _pm = new ParamsMapping(readFile->_pm);
}

ParameterizedReadFile::~ParameterizedReadFile(){
  delete _readFile;
  delete _pm;
}


void ParameterizedReadFile::splitFile(set<ReadFile*>* putFilesHere, int numFiles){
  set<ReadFile*> splitFiles;
  _readFile->splitFile(&splitFiles, numFiles);
  for (set<ReadFile*>::iterator it = splitFiles.begin(); it != splitFiles.end(); ++it){
    putFilesHere->insert( new ParameterizedReadFile( this, *it ) );
    delete *it;
  }
}
void ParameterizedReadFile::splitFile(set<ParameterizedReadFile*>* putFilesHere, int numFiles){
  set<ReadFile*> splitFiles;
  _readFile->splitFile(&splitFiles, numFiles);
  for (set<ReadFile*>::iterator it = splitFiles.begin(); it != splitFiles.end(); ++it){
    putFilesHere->insert( new ParameterizedReadFile( this, *it ) );
    delete *it;
  }
}
ReadFile* ParameterizedReadFile::subFile(long numReadsSkip, long numReadsKeep){
  ReadFile* innerSubFile = _readFile->subFile(numReadsSkip, numReadsKeep);
  ReadFile* subFile = new ParameterizedReadFile( this, innerSubFile );
  delete innerSubFile;
  return subFile;
}
ParameterizedReadFile* ParameterizedReadFile::copy(){
  return new ParameterizedReadFile( this, _readFile );
}


bool ParameterizedReadFile::useThisCycle(int cycleNum){
  if (cycleNum < _initialCycle){ return false; }
  else if (_isCyclic){ return ((cycleNum - _initialCycle) % (_cyclesOn + _cyclesSkipped) < _cyclesOn); }
  else if (_isFinite) { return (cycleNum - _initialCycle < _totalCycles); }
  else { return true; }
}

ParamsMapping* ParameterizedReadFile::getMapParams(){
  return new ParamsMapping(_pm);
}


// inherited interface from ReadFile; these are forwarded to the carried file
string ParameterizedReadFile::getName(){ return _readFile->getName(); }
void ParameterizedReadFile::open(){ _readFile->open(); }
void ParameterizedReadFile::open(bool getBothEnds){ _readFile->open(getBothEnds); }
void ParameterizedReadFile::openWithCarriers(int numSets){ _readFile->openWithCarriers(numSets); }
void ParameterizedReadFile::openWithCarriers(int numSets, bool getBothEnds){ _readFile->openWithCarriers(numSets,getBothEnds); }
void ParameterizedReadFile::openNormalized(){ _readFile->openNormalized(); }
void ParameterizedReadFile::openNormalized(bool getBothEnds){ _readFile->openNormalized(getBothEnds); }
void ParameterizedReadFile::openNormWithCarriers(int numSets){ _readFile->openNormWithCarriers(numSets); }
void ParameterizedReadFile::openNormWithCarriers(int numSets, bool getBothEnds){ _readFile->openNormWithCarriers(numSets,getBothEnds); }
bool ParameterizedReadFile::hasRead(){ return _readFile->hasRead(); }
ScoredSeq* ParameterizedReadFile::getRead(){ return _readFile->getRead(); }
void ParameterizedReadFile::skipRead(){ _readFile->skipRead(); }
void ParameterizedReadFile::skipReads(long numToSkip){ _readFile->skipReads(numToSkip); }
void ParameterizedReadFile::close(){ _readFile->close(); }
long ParameterizedReadFile::numReads(){ return _readFile->numReads(); }
long ParameterizedReadFile::ampliconSize(){ return _readFile->ampliconSize(); }

void ParameterizedReadFile::bufferQueue(Read * r){ throw AssemblyException::CallingError("PRF::bufferQueue"); }
void ParameterizedReadFile::bufferUnqueue(Read * r){ throw AssemblyException::CallingError("PRF::bufferUnqueue"); }
void ParameterizedReadFile::updateBuffer(){ throw AssemblyException::CallingError("PRF::updateBuffer"); }
void ParameterizedReadFile::emptyBuffer(){ throw AssemblyException::CallingError("PRF::emptyBuffer"); }

void ParameterizedReadFile::bufferFill(){ ReadFileCommunicator::bufferFill(_readFile); }
void ParameterizedReadFile::bufferFill(BufferThreadedness threadedness){ ReadFileCommunicator::bufferFill(_readFile, threadedness); }


#endif
