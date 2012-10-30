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

#ifndef PARAMSDEBRUIJN_CPP
#define PARAMSDEBRUIJN_CPP

#include "ParamsDeBruijn.h"
#include "AssemblyException.h"
#include "ScoredSeqSubseq.h"
#include <iostream>
using namespace::std;



float ParamsDeBruijn::defaultMinCount(){ return 1.5; }
int ParamsDeBruijn::defaultKmerSize(){ return 20; }
long ParamsDeBruijn::defaultMaxSeqLength(){ return 100; }
long ParamsDeBruijn::defaultMinNumSeqs(){ return 3; }
// constructors
ParamsDeBruijn::ParamsDeBruijn(){
  _minCount = defaultMinCount();
  _kmerSize = defaultKmerSize();
  _maxSeqLength = defaultMaxSeqLength();
  _minNumSeqs = defaultMinNumSeqs();
}
ParamsDeBruijn::ParamsDeBruijn(int kmerSize) :
  _kmerSize(kmerSize)
{
  _minCount = defaultMinCount();
  _maxSeqLength = defaultMaxSeqLength();
  _minNumSeqs = defaultMinNumSeqs();
}
ParamsDeBruijn::ParamsDeBruijn(int kmerSize, float minCount) :
  _kmerSize(kmerSize),
  _minCount(minCount)
{
  _maxSeqLength = defaultMaxSeqLength();
  _minNumSeqs = defaultMinNumSeqs();
}
ParamsDeBruijn::ParamsDeBruijn(int kmerSize, float minCount, long maxSeqLength) :
  _kmerSize(kmerSize),
  _minCount(minCount),
  _maxSeqLength(maxSeqLength)
{
  _minNumSeqs = defaultMinNumSeqs();
}
ParamsDeBruijn::ParamsDeBruijn(int kmerSize, float minCount, long maxSeqLength, long minNumSeqs) :
  _kmerSize(kmerSize),
  _minCount(minCount),
  _maxSeqLength(maxSeqLength),
  _minNumSeqs(minNumSeqs)
{
}


ParamsDeBruijn::~ParamsDeBruijn(){
}


ParamsDeBruijn::ParamsDeBruijn(ParamsDeBruijn* pdb) :
  _kmerSize(pdb->_kmerSize),
  _minCount(pdb->_minCount),
  _maxSeqLength(pdb->_maxSeqLength),
  _minNumSeqs(pdb->_minNumSeqs)
{
}
ParamsDeBruijn* ParamsDeBruijn::copy(){ return new ParamsDeBruijn(this); }


// check to see if a sequence should go into DeBruijn assembly
bool ParamsDeBruijn::inputToDeBruijn(ScoredSeq* seq){
  if (seq->size() > _maxSeqLength or seq->size() < _kmerSize){ return false; }
  else { return true; }
}


// check to see if an output seq should be kept
bool ParamsDeBruijn::retainDeBruijnOutput(ScoredSeq* seq){
  int pos = 0;
  while (pos < seq->size() and seq->scoreAtPosPlus(pos) < _minCount and
	 (pos == seq->size() - 1 or seq->linkAfterPosPlus(pos) < _minCount) ){ pos++; }
  return pos != seq->size();
}


// get the values
int ParamsDeBruijn::kmerSize(){ return _kmerSize; }
float ParamsDeBruijn::minCount(){ return _minCount; }
long ParamsDeBruijn::maxSeqLength(){ return _maxSeqLength; }
long ParamsDeBruijn::minNumSeqs(){ return _minNumSeqs; }

// set the values
void ParamsDeBruijn::setMinCount(float minCount){ _minCount = minCount; }
void ParamsDeBruijn::setKmerSize(int kmerSize){ _kmerSize = kmerSize; }
void ParamsDeBruijn::setMaxSeqLength(long maxSeqLength){ _maxSeqLength = maxSeqLength; }
void ParamsDeBruijn::setMinNumSeqs(long minNumSeqs){ _minNumSeqs = minNumSeqs; }


#endif
