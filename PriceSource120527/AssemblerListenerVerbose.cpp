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


#ifndef ASSEMBLERLISTENERVERBOSE_CPP
#define ASSEMBLERLISTENERVERBOSE_CPP
#include <string>
#include <sstream>
#include <iostream>
#include <omp.h>
#include <time.h>
#include "AssemblerListenerVerbose.h"
using namespace::std;

AssemblerListenerVerbose::AssemblerListenerVerbose(){
  _edgeUnderway = false;
  _mapUnderway = false;
  _useCout = true;
  _cycleNum = 0;
}
AssemblerListenerVerbose::AssemblerListenerVerbose(string filename){
  _edgeUnderway = false;
  _mapUnderway = false;
  // test creation of the file
  ensureWritability(filename.c_str());
  // now open the file for real
  _logfile.open(filename.c_str());
  _useCout = false;
  _cycleNum = 0;
}

AssemblerListenerVerbose::~AssemblerListenerVerbose(){
  if (! _useCout){ _logfile.close(); }
}

void AssemblerListenerVerbose::printThis(const char* message){
  if (_useCout){
    cout << message;
    cout.flush();
  } else {
    _logfile << message;
    _logfile.flush();
  }
}

void AssemblerListenerVerbose::showParameters(){}
void AssemblerListenerVerbose::assemblyJobBeginStats(AssemblyJob* aj){
  set<ScoredSeq*> inputSeqs;
  aj->allInputSeqs( &inputSeqs );
  if (inputSeqs.size() > 1){
    ostringstream outMessage;
    int threadNum = omp_get_thread_num();
    long totalLen = 0;
    long maxLen = 0;
    long nCount = 0;
    for (set<ScoredSeq*>::iterator it = inputSeqs.begin(); it != inputSeqs.end(); ++it){
      long seqLen = (*it)->size();
      if (seqLen > maxLen){ maxLen = seqLen; }
      totalLen += seqLen;
      char* seq = (*it)->getSeq('+');
      for (long n = 0; n < seqLen; ++n){
	if (seq[n] == 'N'){ nCount++; }
      }
      delete [] seq;
    }
    outMessage << "Thread " << threadNum << ": begin num. seqs:\t" << inputSeqs.size() << endl;
    outMessage << "Thread " << threadNum << ": begin total length:\t" << totalLen << endl;
    outMessage << "Thread " << threadNum << ": begin max length:\t" << maxLen << endl;
    outMessage << "Thread " << threadNum << ": begin num. Ns:\t" << nCount << endl;
    #pragma omp critical (ALStd)
    {
      printThis( outMessage.str().c_str() );
    }
  }
}
void AssemblerListenerVerbose::assemblyJobEndStats(AssemblyJob* aj){
  set<ScoredSeq*> inputSeqs;
  aj->allInputSeqs( &inputSeqs );
  if (inputSeqs.size() > 1){
    ostringstream outMessage;
    int threadNum = omp_get_thread_num();
    set<ScoredSeq*> outputSeqs;
    aj->allAssembledSeqs( &outputSeqs );
    long totalLen = 0;
    long maxLen = 0;
    long nCount = 0;
    for (set<ScoredSeq*>::iterator it = outputSeqs.begin(); it != outputSeqs.end(); ++it){
      long seqLen = (*it)->size();
      if (seqLen > maxLen){ maxLen = seqLen; }
      totalLen += seqLen;
      char* seq = (*it)->getSeq('+');
      for (long n = 0; n < seqLen; ++n){
	if (seq[n] == 'N'){ nCount++; }
      }
      delete [] seq;
    }
    outMessage << "Thread " << threadNum << ": end num. seqs:\t" << outputSeqs.size() << endl;
    outMessage << "Thread " << threadNum << ": end total length:\t" << totalLen << endl;
    outMessage << "Thread " << threadNum << ": end max length:\t" << maxLen << endl;
    outMessage << "Thread " << threadNum << ": end num. Ns:\t" << nCount << endl;
    #pragma omp critical (ALStd)
    {
      printThis( outMessage.str().c_str() );
    }
  }
}


void AssemblerListenerVerbose::beginAssembly(string message){}
void AssemblerListenerVerbose::endAssembly(string message){}

void AssemblerListenerVerbose::beginCycle(int cycleNum){
  _cycleNum = cycleNum;
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << endl << "****** CYCLE " << _cycleNum << " ****** " << endl;
    outMessage << "$START_CYCLE " << asctime(timeinfo);
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::writeOutfile(string outfileName){
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "Writing outfile:  " << outfileName << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::endCycle(set<ScoredSeq*>* contigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "cycle num/num contigs/total size/N50/max size:" << endl
	       << "$END_CYCLE_STATS\t" << _cycleNum 
	       << "\t" << contigs->size()
               << "\t" << getTotalSize(contigs) 
               << "\t" << getN50(contigs)
               << "\t" << getMaxSize(contigs) << endl;
    outMessage << "$END_CYCLE " << asctime(timeinfo); 
    printThis( outMessage.str().c_str() );
  }
}


void AssemblerListenerVerbose::collectedInitialContigs(set<ScoredSeq*>* contigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "$COLLECTED_INITIAL_CONTIGS " << asctime(timeinfo); 
    outMessage << "cycle num/num contigs/total size/N50/max size:" << endl
	       << "$COLLECTED_INITIAL_CONTIGS_STATS\t" << _cycleNum 
	       << "\t" << contigs->size()
               << "\t" << getTotalSize(contigs) 
               << "\t" << getN50(contigs)
               << "\t" << getMaxSize(contigs) << endl;
    printThis( outMessage.str().c_str() );
  }
}

void AssemblerListenerVerbose::filteredInitialContigs(set<ScoredSeq*>* contigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "$FILTERED_INITIAL_CONTIGS " << asctime(timeinfo); 
    outMessage << "cycle num/num contigs/total size/N50/max size:" << endl
	       << "$FILTERED_INITIAL_CONTIGS_STATS\t" << _cycleNum 
	       << "\t" << contigs->size()
               << "\t" << getTotalSize(contigs) 
               << "\t" << getN50(contigs)
               << "\t" << getMaxSize(contigs) << endl;
    printThis( outMessage.str().c_str() );
  }
}

void AssemblerListenerVerbose::beginAddInitialContigs(set<ScoredSeq*>* contigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "$BEGIN_ADD_INITIAL_CONTIGS " << asctime(timeinfo); 
    outMessage << "cycle num/num contigs/total size/N50/max size:" << endl
	       << "$BEGIN_ADD_INITIAL_CONTIGS_STATS\t" << _cycleNum 
	       << "\t" << contigs->size()
               << "\t" << getTotalSize(contigs) 
               << "\t" << getN50(contigs)
               << "\t" << getMaxSize(contigs) << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::endAddInitialContigs(set<ScoredSeq*>* contigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "cycle num/num contigs/total size/N50/max size:" << endl
	       << "$END_ADD_INITIAL_CONTIGS_STATS\t" << _cycleNum 
	       << "\t" << contigs->size()
               << "\t" << getTotalSize(contigs) 
               << "\t" << getN50(contigs)
               << "\t" << getMaxSize(contigs) << endl;
    outMessage << "$END_ADD_INITIAL_CONTIGS " << asctime(timeinfo); 
    printThis( outMessage.str().c_str() );
  }
}


void AssemblerListenerVerbose::beginMapping(long totalNumReads){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    _mapUnderway = true;
    _mapTotal = totalNumReads;
    _mapTally = 0;
    _mapPercent = 0;
    ostringstream outMessage;
    outMessage << "$START_MAP " << asctime(timeinfo);
    outMessage << "Percent of reads (out of " << _mapTotal << ") mapped: " << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::updateMapping(long addReadsToTally){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    _mapTally += addReadsToTally;
    if (_mapTotal > 0){
      int newPercent = _mapTally * long(100) / _mapTotal;
      if (newPercent != _mapPercent){
	_mapPercent = newPercent;
	ostringstream outMessage;
	outMessage << _mapPercent << "% $PROGRESS_MAP " << asctime(timeinfo);
	printThis( outMessage.str().c_str() );
      }
    }
  }
}
void AssemblerListenerVerbose::endMapping(long readsFiltered){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "$END_MAP " << asctime(timeinfo);
    outMessage << "Num. reads filtered: " << readsFiltered << endl;
    printThis( outMessage.str().c_str() );
    _mapUnderway = false;
  }
}

void AssemblerListenerVerbose::beginEdgeAssembly(long totalNumJobs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    _edgeUnderway = true;
    _edgeTotal = totalNumJobs;
    _edgeTally = 0;
    _edgePercent = 0;
    ostringstream outMessage;
    outMessage << "$START_MINI " << asctime(timeinfo);
    outMessage << "Percent of jobs (out of " << _edgeTotal << ") done: " << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::updateEdgeAssembly(long addJobsToTally){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  {
    _edgeTally += addJobsToTally;
    if (_edgeTotal > 0){
      int newPercent = _edgeTally * long(100) / _edgeTotal;
      if (newPercent != _edgePercent){
	_edgePercent = newPercent;
	ostringstream outMessage;
	outMessage << _edgePercent << "% $PROGRESS_MINI " << asctime(timeinfo);
	printThis( outMessage.str().c_str() );
      }
    }
  }
}
void AssemblerListenerVerbose::endEdgeAssembly(){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "$END_MINI " << asctime(timeinfo);
    printThis( outMessage.str().c_str() );
    _edgeUnderway = false;
  }
}

void AssemblerListenerVerbose::beginMetaAssembly(long totalNumContigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "$START_META " << asctime(timeinfo);
    outMessage << "meta-assembly of " << totalNumContigs << " contigs... ";
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::endMetaAssembly(long finalNumContigs){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "collapsed into " << finalNumContigs << " contigs." << endl;
    outMessage << "$END_META " << asctime(timeinfo);
    printThis( outMessage.str().c_str() );
  }
}


void AssemblerListenerVerbose::beginRepeatDetection(){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "$START_REPDETECT " << asctime(timeinfo);
    outMessage << "searching for repeats... ";
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::endRepeatDetection(long numRepeats, long totalLength){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  #pragma omp critical (ALStd)
  { 
    ostringstream outMessage;
    outMessage << "detected " << numRepeats << " repeats (" << totalLength << " nt)." << endl;
    outMessage << "$END_REPDETECT " << asctime(timeinfo);
    printThis( outMessage.str().c_str() );
  }
}


void AssemblerListenerVerbose::message(string message){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  // get the cstr time and "remove" the newline
  char* timeString = asctime(timeinfo);
  int pos = 0;
  while (timeString[pos] != '\n' and timeString[pos] != '\0'){ ++pos; }
  timeString[pos] = '\0';
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "Msg (" << timeString << "): " << message << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::verboseMessage(string message){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  // get the cstr time and "remove" the newline
  char* timeString = asctime(timeinfo);
  int pos = 0;
  while (timeString[pos] != '\n' and timeString[pos] != '\0'){ ++pos; }
  timeString[pos] = '\0';
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "Vmsg (" << timeString << "): " << message << endl;
    printThis( outMessage.str().c_str() );
  }
}
void AssemblerListenerVerbose::errorMessage(string message){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  // get the cstr time and "remove" the newline
  char* timeString = asctime(timeinfo);
  int pos = 0;
  while (timeString[pos] != '\n' and timeString[pos] != '\0'){ ++pos; }
  timeString[pos] = '\0';
  #pragma omp critical (ALStd)
  {
    ostringstream outMessage;
    outMessage << "Error (" << timeString << "): " << message << endl;
    printThis( outMessage.str().c_str() );
  }
}

void AssemblerListenerVerbose::timeStamp(){}

#endif
