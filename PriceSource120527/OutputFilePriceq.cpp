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


#ifndef OUTPUTFILEPRICEQ_CPP
#define OUTPUTFILEPRICEQ_CPP

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "OutputFilePriceq.h"
#include "AssemblyException.h"
using namespace::std;

OutputFilePriceq::OutputFilePriceq(){}
OutputFilePriceq::OutputFilePriceq(string filename){
  _charPerLine = 60; // standard
  _isOpen = false;

  // parse the name
  size_t dotPos = filename.rfind('.');
  if (dotPos == string::npos){
    throw AssemblyException::ArgError("PriceQ filename should be .priceq; '.' is missing");
  }
  _filenameEnd = filename.substr(dotPos);
  if (_filenameEnd != ".priceq"){
    cout << _filenameEnd << endl;
    throw AssemblyException::ArgError("PriceQ filename should be .priceq.txt; append is wrong");
  }
  _filenameBegin = filename.substr(0, dotPos);

  // ensure writability of the base name
  ensureWritability( getName().c_str() );
}
OutputFilePriceq::~OutputFilePriceq(){}

void OutputFilePriceq::open(){
  if (! _isOpen){
    _outputFileObj.open( getName().c_str() );
    _isOpen = true;
    _contigCounter = 0;
  }
}
void OutputFilePriceq::open(string addendum){
  if (! _isOpen){
    _outputFileObj.open( getName(addendum).c_str() );
    _isOpen = true;
    _contigCounter = 0;
  }
}
bool OutputFilePriceq::isOpen(){ return _isOpen; }

void OutputFilePriceq::writeContigs(ScoredSeq** contigs){
  set<ScoredSeq*> noTerminal;
  writeContigs(contigs, &noTerminal);
}

void OutputFilePriceq::writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchangedContigs){
  if(! _isOpen){
    throw AssemblyException::CallingError("cannot write contigs to fasta file when it isn't open.");
  }

  long index = 0;
  while (contigs[index] != NULL){
    _contigCounter++;
    ScoredSeq* seq = contigs[index];
    index++;
    long seqSize = seq->size();
    char* tempSeq;

    // the seq
    _outputFileObj << "@contig_" << _contigCounter << " (" << seqSize << "nt)";
    if (unchangedContigs->count(seq) > 0){ _outputFileObj << " unchanged"; }
    _outputFileObj << "\n";
    tempSeq = seq->getSeq('+');
    _outputFileObj << tempSeq << "\n";
    delete [] tempSeq;

    // the nuc scores
    _outputFileObj << "+\n";
    for (long n = 0; n < seqSize; ++n){ _outputFileObj << scoreToChar( seq->scoreAtPosPlus(n) ); }
    _outputFileObj << "\n";

    // the link scores
    _outputFileObj << "~\n";
    for (long n = 0; n < seqSize - 1; ++n){ _outputFileObj << scoreToChar( seq->linkAfterPosPlus(n) ); }
    _outputFileObj << "\n";
  }
}

char OutputFilePriceq::scoreToChar(float score){
  // 126 is the max value; 126 - 35 = 91
  int ascInt;
  if (score >= 1.0){
    float logged = log(score) / log(float(1.5));
    if (logged > float(91)){ ascInt = 126; }
    else { ascInt = int(logged) + 35; }
  }
  else if (score < 0.3){ ascInt = 33; }
  else if (score < 0.9){ ascInt = 34; }
  else { ascInt = 35; }
  return char(ascInt);
}
float OutputFilePriceq::charToScore(char scoreChar){
  int ascInt = int(scoreChar);
  switch( ascInt ){
  case 33: return float(0.0);
  case 34: return float(0.3);
  default: return pow(float(1.5), float(ascInt - 35));
  }
}

void OutputFilePriceq::close(){
  if (_isOpen){
    _outputFileObj.close();
    _isOpen = false;
  }
}
string OutputFilePriceq::getName(){
  return _filenameBegin + _filenameEnd;
}
string OutputFilePriceq::getName(string addendum){
  return _filenameBegin + "." + addendum + _filenameEnd;
}

#endif
