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

#ifndef PARAMETERIZEDINITIALFILE_CPP
#define PARAMETERIZEDINITIALFILE_CPP

#include "ParameterizedInitialFile.h"
#include "AssemblyException.h"
#include <iostream>
using namespace::std;

// constructors
ParameterizedInitialFile::ParameterizedInitialFile(){}

ParameterizedInitialFile::ParameterizedInitialFile(ReadFile* readFile) :
  _initialCycle(0)
{
  _readFile = readFile->copy();
  _totalInput = _readFile->numReads();
  _totalOutput = _totalInput;
  makeInitDeploymentArray(1, 1);
}
ParameterizedInitialFile::ParameterizedInitialFile(ReadFile* readFile, int numSteps, int cyclesPerStep) :
  _initialCycle(0)
{
  _readFile = readFile->copy();
  _totalInput = _readFile->numReads();
  _totalOutput = _totalInput;
  makeInitDeploymentArray(numSteps, cyclesPerStep);
}
ParameterizedInitialFile::ParameterizedInitialFile(ReadFile* readFile, int numSteps, int cyclesPerStep, int initialCycle) :
  _initialCycle(initialCycle)
{
  _readFile = readFile->copy();
  _totalInput = _readFile->numReads();
  _totalOutput = _totalInput;
  makeInitDeploymentArray(numSteps, cyclesPerStep);
}

ParameterizedInitialFile::ParameterizedInitialFile(long totalOutput, ReadFile* readFile, int numSteps, int cyclesPerStep) :
  _totalOutput(totalOutput),
  _initialCycle(0)
{
  _readFile = readFile->copy();
  _totalInput = _readFile->numReads();
  makeInitDeploymentArray(numSteps, cyclesPerStep);
}
ParameterizedInitialFile::ParameterizedInitialFile(long totalOutput, ReadFile* readFile, int numSteps, int cyclesPerStep, int initialCycle) :
  _totalOutput(totalOutput),
  _initialCycle(initialCycle)
{
  _readFile = readFile->copy();
  _totalInput = _readFile->numReads();
  makeInitDeploymentArray(numSteps, cyclesPerStep);
}

ParameterizedInitialFile::~ParameterizedInitialFile(){
  delete [] _cycleToInputNum;
  delete [] _cycleToOutputNum;
  delete [] _cycleToPastInputCount;
  delete _readFile;
}


void ParameterizedInitialFile::makeInitDeploymentArray(int numSteps, int cyclesPerStep){
  _totalCycles = numSteps * cyclesPerStep + _initialCycle;

  _cycleToInputNum = new long[ _totalCycles ];
  _cycleToOutputNum = new long[ _totalCycles ];
  _cycleToPastInputCount = new long[ _totalCycles ];

  long alreadyAdded = 0;
  long alreadyRead = 0;
  int stepsTaken = 0;
  for (int n = 0; n < _initialCycle; ++n){
    _cycleToPastInputCount[n] = 0;
    _cycleToInputNum[n] = 0;
    _cycleToOutputNum[n] = 0;
  }
  for (int n  = _initialCycle; n < _totalCycles; ++n){
    _cycleToPastInputCount[n] = alreadyRead;

    if ( (n - _initialCycle) % cyclesPerStep == 0 and stepsTaken < numSteps){
      stepsTaken++;

      long priorAlreadyRead = alreadyRead;
      alreadyRead = stepsTaken * _totalInput / numSteps;
      _cycleToInputNum[n] = alreadyRead - priorAlreadyRead;

      long priorAlreadyAdded = alreadyAdded;
      alreadyAdded = stepsTaken * _totalOutput / numSteps;
      _cycleToOutputNum[n] = alreadyAdded - priorAlreadyAdded;

    } else {
      _cycleToInputNum[n] = 0;
      _cycleToOutputNum[n] = 0;
    }
  }

  // make sure that the tally adds up correctly
  if (alreadyAdded != _totalOutput){
    cout << alreadyAdded << endl;
    cout << _totalOutput << endl;
    throw AssemblyException::LogicError("the number of initial contigs scheduled to be added does not match the total.");
  }
  if (alreadyRead != _totalInput){
    cout << alreadyRead << endl;
    cout << _totalInput << endl;
    throw AssemblyException::LogicError("the number of initial contigs scheduled to be read does not match the total.");
  }
}



// REQUIRES len(selectionArray) == _cycleToInputNum[cycle]
// MODIFIES selectionArray
void ParameterizedInitialFile::fillCycleSelectionArray(bool* selectionArray, int cycle){

  // used to decide when a read should be collected;
  // will be one less than the number of "True" entries
  int intOutputCount = -1;

  // when this reaches a new int value, enter "True"
  float outCountPerRead = float(_cycleToOutputNum[cycle]) / _cycleToInputNum[cycle];

  for (long n = 0; n < _cycleToInputNum[cycle]; ++n){
    float floatOutputCount = outCountPerRead * n;
    if ( int(floatOutputCount) > intOutputCount ){
      selectionArray[n] = true;
      intOutputCount++;
    } else {
      selectionArray[n] = false;
    }
  }
}






int ParameterizedInitialFile::totalCycles(){ return _totalCycles; }

void ParameterizedInitialFile::getContigs(set<ScoredSeq*>* initialContigs, int cycleNum){

  if ( cycleNum < _totalCycles and _cycleToOutputNum[cycleNum] > 0 ){

    // skip the already-read reads
    _readFile->open();
    _readFile->skipReads( _cycleToPastInputCount[cycleNum] );

    // get the new contigs
    bool* selectionArray = new bool[ _cycleToInputNum[cycleNum] ];
    fillCycleSelectionArray(selectionArray, cycleNum);
    for (long n = 0; n < _cycleToInputNum[cycleNum]; ++n){
      if ( selectionArray[n] ){
	ScoredSeq* newRead = _readFile->getRead();
	initialContigs->insert( newRead->shallowCopy() );
	newRead->deepDelete();
      } else {
	_readFile->skipRead();
      }
    }
    delete [] selectionArray;

    _readFile->close();
  }
}

#endif
