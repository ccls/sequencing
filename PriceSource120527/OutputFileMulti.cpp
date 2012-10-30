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


#ifndef OUTPUTFILEMULTI_CPP
#define OUTPUTFILEMULTI_CPP

#include <string>
#include <iostream>
#include "OutputFileMulti.h"
#include "AssemblyException.h"
using namespace::std;

OutputFileMulti::OutputFileMulti(OutputFile** fileArray, int numFiles){
  _isOpen = false;
  _numFiles = numFiles;
  _fileArray = new OutputFile*[ numFiles + 1 ];
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n] = fileArray[n]; }
}
OutputFileMulti::~OutputFileMulti(){
  delete [] _fileArray;
}

void OutputFileMulti::open(){
  _isOpen = true;
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n]->open(); }
}
void OutputFileMulti::open(string addendum){
  _isOpen = true;
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n]->open(addendum); }
}
bool OutputFileMulti::isOpen(){ return _isOpen; }
void OutputFileMulti::writeContigs(ScoredSeq** contigs){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to null file when it isn't open.");
  }
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n]->writeContigs(contigs); }
}
void OutputFileMulti::writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchanged){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to null file when it isn't open.");
  }
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n]->writeContigs(contigs, unchanged); }
}
void OutputFileMulti::close(){
  _isOpen = false;
  for (int n = 0; n < _numFiles; ++n){ _fileArray[n]->close(); }
}
string OutputFileMulti::getName(){ 
  string currentString = "";
  for (int n = 0; n < _numFiles; ++n){
    if (n > 0){ currentString += ", "; }
    currentString += _fileArray[n]->getName();
  }
  return currentString;
}
string OutputFileMulti::getName(string addendum){ 
  string currentString = "";
  for (int n = 0; n < _numFiles; ++n){
    if (n > 0){ currentString += ", "; }
    currentString += _fileArray[n]->getName(addendum);
  }
  return currentString;
}



#endif
