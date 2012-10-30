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


#ifndef ASSEMBLERLISTENERCARRIER_CPP
#define ASSEMBLERLISTENERCARRIER_CPP
#include <string>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <cstring>
#include "AssemblerListenerCarrier.h"
#include "AssemblerListenerNull.h"
using namespace::std;

AssemblerListenerCarrier::AssemblerListenerCarrier(){}
AssemblerListenerCarrier::AssemblerListenerCarrier(AssemblerListener** listeners, int numListeners) :
  _numListeners(numListeners)
{
  _listeners = new AssemblerListener*[numListeners];
  for (int n = 0; n < _numListeners; ++n){ _listeners[n] = listeners[n]; }
}
AssemblerListenerCarrier::~AssemblerListenerCarrier(){
  //for (int n = 0; n < _numListeners; ++n){ delete _listeners[n]; }
  delete [] _listeners;
}


void AssemblerListenerCarrier::showParameters(){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->showParameters(); }
}
void AssemblerListenerCarrier::assemblyJobBeginStats(AssemblyJob* aj){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->assemblyJobBeginStats(aj); }
}
void AssemblerListenerCarrier::assemblyJobEndStats(AssemblyJob* aj){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->assemblyJobEndStats(aj); }
}

void AssemblerListenerCarrier::beginAssembly(string message){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginAssembly(message); }
}
void AssemblerListenerCarrier::endAssembly(string message){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endAssembly(message); }
}

void AssemblerListenerCarrier::beginCycle(int cycleNum){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginCycle(cycleNum); }
}
void AssemblerListenerCarrier::writeOutfile(string outfileName){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->writeOutfile(outfileName); }
}
void AssemblerListenerCarrier::endCycle(set<ScoredSeq*>* contigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endCycle(contigs); }
}

void AssemblerListenerCarrier::collectedInitialContigs(set<ScoredSeq*>* contigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->collectedInitialContigs(contigs); }
}
void AssemblerListenerCarrier::filteredInitialContigs(set<ScoredSeq*>* contigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->filteredInitialContigs(contigs); }
}
void AssemblerListenerCarrier::beginAddInitialContigs(set<ScoredSeq*>* contigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginAddInitialContigs(contigs); }
}
void AssemblerListenerCarrier::endAddInitialContigs(set<ScoredSeq*>* contigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endAddInitialContigs(contigs); }
}

void AssemblerListenerCarrier::beginMapping(long totalNumReads){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginMapping(totalNumReads); }
}
void AssemblerListenerCarrier::updateMapping(long addReadsToTally){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->updateMapping(addReadsToTally); }
}
void AssemblerListenerCarrier::endMapping(long readsFiltered){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endMapping(readsFiltered); }
}

void AssemblerListenerCarrier::beginEdgeAssembly(long totalNumJobs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginEdgeAssembly(totalNumJobs); }
}
void AssemblerListenerCarrier::updateEdgeAssembly(long addJobsToTally){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->updateEdgeAssembly(addJobsToTally); }
}
void AssemblerListenerCarrier::endEdgeAssembly(){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endEdgeAssembly(); }
}

void AssemblerListenerCarrier::beginMetaAssembly(long totalNumContigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginMetaAssembly(totalNumContigs); }
}
void AssemblerListenerCarrier::endMetaAssembly(long finalNumContigs){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endMetaAssembly(finalNumContigs); }
}

void AssemblerListenerCarrier::beginRepeatDetection(){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->beginRepeatDetection(); }
}
void AssemblerListenerCarrier::endRepeatDetection(long numRepeats, long totalLength){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->endRepeatDetection(numRepeats, totalLength); }
}

void AssemblerListenerCarrier::message(string message){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->message(message); }
}
void AssemblerListenerCarrier::verboseMessage(string message){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->verboseMessage(message); }
}
void AssemblerListenerCarrier::errorMessage(string message){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->errorMessage(message); }
}

void AssemblerListenerCarrier::timeStamp(){
  for (int n = 0; n < _numListeners; ++n){ _listeners[n]->timeStamp(); }
}

#endif
