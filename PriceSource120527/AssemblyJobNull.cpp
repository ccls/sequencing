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


#ifndef ASSEMBLYJOBNULL_CPP
#define ASSEMBLYJOBNULL_CPP

#include "AssemblyJobNull.h"

#include "AssemblyException.h"
#include <queue>
#include <iostream>
using namespace::std;

AssemblyJobNull::AssemblyJobNull(){}

AssemblyJobNull::AssemblyJobNull(set<ScoredSeq*>* seqs) :
  _inputRemoved(false),
  _wasRun(false){
  _inputSeqs.insert( seqs->begin(), seqs->end() );
}
AssemblyJobNull::~AssemblyJobNull(){}



void AssemblyJobNull::makeShallow(bool shallowDelete){
  if (_wasRun){
    throw AssemblyException::CallingError("AJNull::makeShallow can't be called after the job was run.");
  }
  set<ScoredSeq*> shallowCopies;
  // buffer all input, then copy
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    (*it)->bottomBuffer();
  }
  for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
    shallowCopies.insert( (*it)->shallowCopy() );
    (*it)->deepUnbuffer();
  }
  if (shallowDelete){
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      delete (*it);
    }
  }
  // replace the input contents with the shallow copies
  _inputSeqs.clear();
  _inputSeqs.insert( shallowCopies.begin(), shallowCopies.end() );
}
void AssemblyJobNull::runJob(AssemblyJob::AssemblyThreadedness threadedness){ 
  _wasRun = true;
}
bool AssemblyJobNull::wasRun(){ return _wasRun; }

void AssemblyJobNull::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}
void AssemblyJobNull::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}
void AssemblyJobNull::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  // nothing happens because no seqs were discarded
}
void AssemblyJobNull::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  // nothing happens because no seqs were created
}
void AssemblyJobNull::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}



void AssemblyJobNull::OK(){}



#endif
