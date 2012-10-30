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



#ifndef ALIGNMENTNULL_CPP
#define ALIGNMENTNULL_CPP

#include "AlignmentNull.h"

#include "AssemblyException.h"
#include <iostream>
#include <typeinfo>
#include <set>
using namespace std;


AlignmentNull::~AlignmentNull(){}
AlignmentNull::AlignmentNull(ScoredSeq * seqA, ScoredSeq* seqB, char orientation) :
  _seqA(seqA), _seqB(seqB), _orientation(orientation)
{}

AlignmentNull* AlignmentNull::copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentNull::copySeqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentNull::copySeqReplace seqsB aren't equal"); }
  return new AlignmentNull(seqA, seqB, _orientation);
}
void AlignmentNull::seqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentNull::seqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentNull::seqReplace seqsB aren't equal"); }
  _seqA = seqA;
  _seqB = seqB;
}
AlignmentNull* AlignmentNull::copyRcSeqA(ScoredSeq* newSeqA){
  if (newSeqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentNull::copyRcSeqA seqsA aren't equal"); }
  switch (_orientation){
  case '+': return new AlignmentNull(newSeqA,_seqB,'-');
  case '-': return new AlignmentNull(newSeqA,_seqB,'+');
  default: throw AssemblyException::StateError("AlNull::copyRcSeqA bad sense");
  }
}

ScoredSeq * AlignmentNull::seqA(){ return _seqA; }
ScoredSeq* AlignmentNull::seqB(){ return _seqB; }
char AlignmentNull::orientation(){ return _orientation; }
bool AlignmentNull::isNull(){ return true; }
bool AlignmentNull::isLinked(long pos, ScoredSeq * seq){
  throw AssemblyException::CallingError("AlignmentNull is always null (isLinked)");
} 
bool AlignmentNull::isGapped(long pos, ScoredSeq * seq){ 
  throw AssemblyException::CallingError("AlignmentNull is always null (isGapped)");
} 
long AlignmentNull::getLinkage(long pos, ScoredSeq * seq){
  throw AssemblyException::CallingError("AlignmentNull is always null (getLinkage)");
} 
long AlignmentNull::gapPairedAfter(long pos, ScoredSeq * seq){
  throw AssemblyException::CallingError("AlignmentNull is always null (gapPairedAfter)");
}
long AlignmentNull::score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps){
  throw AssemblyException::CallingError("AlignmentNull is always null (score)");
} 
long AlignmentNull::scoreOverhangA(AlignmentScoreMatrix* scoreMatrix){
  throw AssemblyException::CallingError("AlignmentNull is always null (score)");
} 


bool AlignmentNull::isLocked(){ return true; }


ScoredSeq* AlignmentNull::combine(){
  throw AssemblyException::CallingError("AlignmentNull is always null (combine1)");
} 
ScoredSeq* AlignmentNull::combine(int denominator){
  throw AssemblyException::CallingError("AlignmentNull is always null (combine2)");
} 
bool AlignmentNull::hasConstantOffset(){
  throw AssemblyException::CallingError("AlignmentNull is always null (hasConstantOffset)");
}
long AlignmentNull::getConstantOffset(){
  throw AssemblyException::CallingError("AlignmentNull is always null (getConstantOffset)");
}


void AlignmentNull::OK(){}



#endif
