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


#ifndef ASSEMBLYJOBNULLCOPY_CPP
#define ASSEMBLYJOBNULLCOPY_CPP

#include "AssemblyJobNullCopy.h"

#include "AssemblyException.h"
#include <queue>
#include <iostream>
using namespace::std;

AssemblyJobNullCopy::AssemblyJobNullCopy(){}

AssemblyJobNullCopy::AssemblyJobNullCopy(set<ScoredSeq*>* seqs) :
  _inputRemoved(false),
  _isShallow(false),
  _wasRun(false){
  _inputSeqs.insert( seqs->begin(), seqs->end() );
}
AssemblyJobNullCopy::~AssemblyJobNullCopy(){}



void AssemblyJobNullCopy::makeShallow(bool shallowDelete){
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

void AssemblyJobNullCopy::runJob(AssemblyJob::AssemblyThreadedness threadedness){
  if (! _wasRun){
    // buffer all input, then copy
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      (*it)->bottomBuffer();
    }
    for (set<ScoredSeq*>::iterator it = _inputSeqs.begin(); it != _inputSeqs.end(); ++it){
      _copySeqs.insert( (*it)->shallowCopy() );
      (*it)->deepUnbuffer();
    }
    _wasRun = true;
  }
}

bool AssemblyJobNullCopy::wasRun(){ return _wasRun; }

void AssemblyJobNullCopy::allInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}
void AssemblyJobNullCopy::remainingInputSeqs( set<ScoredSeq*>* outsideSet ){
}
void AssemblyJobNullCopy::discardedInputSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _inputSeqs.begin(), _inputSeqs.end() );
}
void AssemblyJobNullCopy::newlyAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _copySeqs.begin(), _copySeqs.end() );
}
void AssemblyJobNullCopy::allAssembledSeqs( set<ScoredSeq*>* outsideSet ){
  outsideSet->insert( _copySeqs.begin(), _copySeqs.end() );
}



void AssemblyJobNullCopy::OK(){}



#endif
