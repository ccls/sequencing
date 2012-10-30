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


#ifndef ASSEMBLERLISTENERSTD_CPP
#define ASSEMBLERLISTENERSTD_CPP
#include <string>
#include <iostream>
#include <omp.h>
#include "AssemblerListenerStd.h"
using namespace::std;

AssemblerListenerStd::AssemblerListenerStd(){
  _edgeUnderway = false;
  _mapUnderway = false;
}
AssemblerListenerStd::~AssemblerListenerStd(){}

void AssemblerListenerStd::showParameters(){}
void AssemblerListenerStd::assemblyJobBeginStats(AssemblyJob* aj){}
void AssemblerListenerStd::assemblyJobEndStats(AssemblyJob* aj){}

void AssemblerListenerStd::beginAssembly(string message){}
void AssemblerListenerStd::endAssembly(string message){}

void AssemblerListenerStd::beginCycle(int cycleNum){
  #pragma omp critical (ALStd)
  { cout << endl << "****** CYCLE " << cycleNum << " ****** " << endl; }
}
void AssemblerListenerStd::writeOutfile(string outfileName){
  #pragma omp critical (ALStd)
  { cout << "Writing outfile:  " << outfileName << endl; }
}
void AssemblerListenerStd::endCycle(set<ScoredSeq*>* contigs){
  #pragma omp critical (ALStd)
  {
    cout << "End Cycle: num contigs/total size/N50/max size = " << contigs->size()
         << " / " << getTotalSize(contigs) 
         << " / " << getN50(contigs)
         << " / " << getMaxSize(contigs) << endl;
  }
}

void AssemblerListenerStd::collectedInitialContigs(set<ScoredSeq*>* contigs){
  #pragma omp critical (ALStd)
  {
    cout << "Collected new initial contigs: num contigs/total size/N50/max size = " << contigs->size()
         << " / " << getTotalSize(contigs) 
         << " / " << getN50(contigs)
         << " / " << getMaxSize(contigs) << endl;
  }
}

void AssemblerListenerStd::filteredInitialContigs(set<ScoredSeq*>* contigs){
  #pragma omp critical (ALStd)
  {
    cout << "Filtered initial contigs remaining: num contigs/total size/N50/max size = " << contigs->size()
         << " / " << getTotalSize(contigs) 
         << " / " << getN50(contigs)
         << " / " << getMaxSize(contigs) << endl;
  }
}

void AssemblerListenerStd::beginAddInitialContigs(set<ScoredSeq*>* contigs){
  #pragma omp critical (ALStd)
  {
    cout << "Adding initial contigs: num contigs/total size/N50/max size = " << contigs->size()
         << " / " << getTotalSize(contigs) 
         << " / " << getN50(contigs)
         << " / " << getMaxSize(contigs) << endl;
  }
}
void AssemblerListenerStd::endAddInitialContigs(set<ScoredSeq*>* contigs){
  #pragma omp critical (ALStd)
  {
    cout << "Added initial contigs: num contigs/total size/N50/max size = " << contigs->size()
         << " / " << getTotalSize(contigs) 
         << " / " << getN50(contigs)
         << " / " << getMaxSize(contigs) << endl;
  }
}

void AssemblerListenerStd::beginMapping(long totalNumReads){
  #pragma omp critical (ALStd)
  {
    _mapUnderway = true;
    _mapTotal = totalNumReads;
    _mapTally = 0;
    _mapPercent = 0;
    cout << "Percent of reads (out of " << _mapTotal << ") mapped: " << _mapPercent;
    cout.flush();
  }
}
void AssemblerListenerStd::updateMapping(long addReadsToTally){
  #pragma omp critical (ALStd)
  { 
    _mapTally += addReadsToTally;
    if (_mapTotal > 0){
      int newPercent = _mapTally * long(100) / _mapTotal;
      if (newPercent != _mapPercent){
	if (_mapPercent < 10){ cout << "\b"; }
	else if (_mapPercent < 100){ cout << "\b\b"; }
	else { cout << "\b\b\b"; }
	_mapPercent = newPercent;
	cout << _mapPercent;
	cout.flush();
      }
    }
  }
}
void AssemblerListenerStd::endMapping(long readsFiltered){
  #pragma omp critical (ALStd)
  { 
    cout << endl;
    _mapUnderway = false;
    cout << "Num. reads filtered: " << readsFiltered << endl;
  }
}

void AssemblerListenerStd::beginEdgeAssembly(long totalNumJobs){
  #pragma omp critical (ALStd)
  { 
    _edgeUnderway = true;
    _edgeTotal = totalNumJobs;
    _edgeTally = 0;
    _edgePercent = 0;
    cout << "Percent of jobs (out of " << _edgeTotal << ") done: " << _edgePercent;
    cout.flush();
  }
}
void AssemblerListenerStd::updateEdgeAssembly(long addJobsToTally){
  #pragma omp critical (ALStd)
  {
    _edgeTally += addJobsToTally;
    if (_edgeTotal > 0){
      int newPercent = _edgeTally * long(100) / _edgeTotal;
      if (newPercent != _edgePercent){
	if (_edgePercent < 10){ cout << "\b"; }
	else if (_edgePercent < 100){ cout << "\b\b"; }
	else { cout << "\b\b\b"; }
	_edgePercent = newPercent;
	cout << _edgePercent;
	cout.flush();
      }
    }
  }
}
void AssemblerListenerStd::endEdgeAssembly(){
  #pragma omp critical (ALStd)
  { 
    cout << endl;
    _edgeUnderway = false;
  }
}

void AssemblerListenerStd::beginMetaAssembly(long totalNumContigs){
  #pragma omp critical (ALStd)
  { 
    cout << "meta-assembly of " << totalNumContigs << " contigs... ";
    cout.flush();
  }
}
void AssemblerListenerStd::endMetaAssembly(long finalNumContigs){
  #pragma omp critical (ALStd)
  { 
    cout << "collapsed into " << finalNumContigs << " contigs." << endl;
  }
}


void AssemblerListenerStd::beginRepeatDetection(){
  #pragma omp critical (ALStd)
  { 
    cout << "searching for repeats... ";
    cout.flush();
  }
}
void AssemblerListenerStd::endRepeatDetection(long numRepeats, long totalLength){
  #pragma omp critical (ALStd)
  { 
    cout << "detected " << numRepeats << " repeats (" << totalLength << " nt)." << endl;
  }
}


void AssemblerListenerStd::message(string message){
  #pragma omp critical (ALStd)
  {
    cout << "Msg: " << message << endl;
  }
}
void AssemblerListenerStd::verboseMessage(string message){}
void AssemblerListenerStd::errorMessage(string message){
  #pragma omp critical (ALStd)
  {
    cout << "Error: " << message << endl;
  }
}

void AssemblerListenerStd::timeStamp(){}

#endif
