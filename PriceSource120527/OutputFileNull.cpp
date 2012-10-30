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


#ifndef OUTPUTFILENULL_CPP
#define OUTPUTFILENULL_CPP

#include <string>
#include <iostream>
#include "OutputFileNull.h"
#include "AssemblyException.h"
using namespace::std;

OutputFileNull::OutputFileNull(){ _isOpen = false; }
OutputFileNull::~OutputFileNull(){}

void OutputFileNull::open(){ _isOpen = true; }
void OutputFileNull::open(string addendum){ _isOpen = true; }
bool OutputFileNull::isOpen(){ return _isOpen; }
void OutputFileNull::writeContigs(ScoredSeq** contigs){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to null file when it isn't open.");
  }
}
void OutputFileNull::writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchanged){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to null file when it isn't open.");
  }
}
void OutputFileNull::close(){ _isOpen = false; }
string OutputFileNull::getName(){ return ""; }
string OutputFileNull::getName(string addendum){ return ""; }



#endif
