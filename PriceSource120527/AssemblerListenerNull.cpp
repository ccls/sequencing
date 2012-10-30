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


#ifndef ASSEMBLERLISTENERNULL_CPP
#define ASSEMBLERLISTENERNULL_CPP

#include "AssemblerListenerNull.h"
#include <string>
using namespace::std;

AssemblerListenerNull::AssemblerListenerNull(){}
AssemblerListenerNull::~AssemblerListenerNull(){}

void AssemblerListenerNull::showParameters(){}
void AssemblerListenerNull::assemblyJobBeginStats(AssemblyJob* aj){}
void AssemblerListenerNull::assemblyJobEndStats(AssemblyJob* aj){}

void AssemblerListenerNull::beginAssembly(string message){}
void AssemblerListenerNull::endAssembly(string message){}

void AssemblerListenerNull::beginCycle(int cycleNum){}
void AssemblerListenerNull::writeOutfile(string outfileName){}
void AssemblerListenerNull::endCycle(set<ScoredSeq*>* contigs){}

void AssemblerListenerNull::collectedInitialContigs(set<ScoredSeq*>* contigs){}
void AssemblerListenerNull::filteredInitialContigs(set<ScoredSeq*>* contigs){}
void AssemblerListenerNull::beginAddInitialContigs(set<ScoredSeq*>* contigs){}
void AssemblerListenerNull::endAddInitialContigs(set<ScoredSeq*>* contigs){}

void AssemblerListenerNull::beginMapping(long totalNumReads){}
void AssemblerListenerNull::updateMapping(long addReadsToTally){}
void AssemblerListenerNull::endMapping(long readsFiltered){}

void AssemblerListenerNull::beginEdgeAssembly(long totalNumJobs){}
void AssemblerListenerNull::updateEdgeAssembly(long addJobsToTally){}
void AssemblerListenerNull::endEdgeAssembly(){}

void AssemblerListenerNull::beginMetaAssembly(long totalNumContigs){}
void AssemblerListenerNull::endMetaAssembly(long finalNumContigs){}

void AssemblerListenerNull::beginRepeatDetection(){}
void AssemblerListenerNull::endRepeatDetection(long numRepeats, long totalLength){}

void AssemblerListenerNull::message(string message){}
void AssemblerListenerNull::verboseMessage(string message){}
void AssemblerListenerNull::errorMessage(string message){}

void AssemblerListenerNull::timeStamp(){}

#endif
