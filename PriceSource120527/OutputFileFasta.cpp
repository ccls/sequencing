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


#ifndef OUTPUTFILEFASTA_CPP
#define OUTPUTFILEFASTA_CPP

#include <string>
#include <iostream>
#include <algorithm>
#include <cstring>
#include "OutputFileFasta.h"
#include "AssemblyException.h"
using namespace::std;

OutputFileFasta::OutputFileFasta(){}
OutputFileFasta::OutputFileFasta(string filename){
  _charPerLine = 60; // standard
  _isOpen = false;

  // parse the name
  size_t dotPos = filename.rfind('.');
  if (dotPos == string::npos){
    throw AssemblyException::ArgError("fasta filename should be .fa or .fasta; '.' is missing");
  }
  string filenameEnd = filename.substr(dotPos);
  if (filenameEnd != ".fa" and filenameEnd != ".fasta"){
    cout << filenameEnd << endl;
    throw AssemblyException::ArgError("fasta filename should be .fa or .fasta; addendum is wrong");
  }
  _filenameEnd = new char[ filenameEnd.size()+1 ];
  strcpy (_filenameEnd, filenameEnd.c_str());
  _filenameEnd[ filenameEnd.size() ] = '\0';

  string filenameBegin = filename.substr(0, dotPos);
  _filenameBegin = new char[ filenameBegin.size()+1 ];
  strcpy (_filenameBegin, filenameBegin.c_str());
  _filenameBegin[ filenameBegin.size() ] = '\0';

  // ensure writability of the base name
  ensureWritability( getName().c_str() );
}
OutputFileFasta::~OutputFileFasta(){
  delete [] _filenameBegin;
  delete [] _filenameEnd;
}

void OutputFileFasta::open(){
  if (! _isOpen){
    _outputFileObj.open( getName().c_str() );
    _isOpen = true;
    _contigCounter = 0;
  }
}
void OutputFileFasta::open(string addendum){
  if (! _isOpen){
    _outputFileObj.open( getName(addendum).c_str() );
    _isOpen = true;
    _contigCounter = 0;
  }
}
bool OutputFileFasta::isOpen(){ return _isOpen; }

void OutputFileFasta::writeContigs(ScoredSeq** contigs){
  set<ScoredSeq*> noTerminal;
  writeContigs(contigs, &noTerminal);
}

void OutputFileFasta::writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchangedContigs){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to fasta file when it isn't open.");
  }

  long index = 0;
  while (contigs[index] != NULL){
    _contigCounter++;
    ScoredSeq* seq = contigs[index];
    index++;
    long seqSize = seq->size();
    _outputFileObj << ">contig_" << _contigCounter << " (" << seqSize << "nt)";
    if (unchangedContigs->count(seq) > 0){ _outputFileObj << " unchanged"; }
    _outputFileObj << "\n";

    long startPos = 0;
    while (startPos < seqSize){
      int lineLen;
      if (startPos + _charPerLine > seqSize){ lineLen = seqSize - startPos; }
      else { lineLen = _charPerLine; }
      char* tempSeq = seq->getSubseq(startPos, lineLen, '+');
      _outputFileObj << tempSeq << "\n";
      delete [] tempSeq;
      startPos += _charPerLine;
    }
  }
}

void OutputFileFasta::close(){
  if (_isOpen){
    _outputFileObj.close();
    _isOpen = false;
  }
}
string OutputFileFasta::getName(){
  return string(_filenameBegin) + string(_filenameEnd);
}
string OutputFileFasta::getName(string addendum){
  return string(_filenameBegin) + "." + addendum + string(_filenameEnd);
}


#endif
