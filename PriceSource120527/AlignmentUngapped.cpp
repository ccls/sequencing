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



#ifndef ALIGNMENTUNGAPPED_CPP
#define ALIGNMENTUNGAPPED_CPP

#include "AlignmentUngapped.h"

#include "AssemblyException.h"
#include <iostream>
#include <typeinfo>
#include <set>
using namespace std;




AlignmentUngapped::~AlignmentUngapped(){}
AlignmentUngapped::AlignmentUngapped(ScoredSeq * seqA, ScoredSeq* seqB, char orientation, long offset) :
  _seqA(seqA), _seqB(seqB), _orientation(orientation), _offset(offset)
{ 
  if (offset > 0){ _startA = offset; _startB = 0; }
  else { _startA = 0; _startB = 0 - offset; }
  if ( _seqA->size() - _startA > _seqB->size() - _startB ){
    _endA = _startA + _seqB->size() - _startB;
    _endB = _seqB->size();
  } else {
    _endA = _seqA->size();
    _endB = _startB + _seqA->size() - _startA;
  }
  /*
  _endA = offset + _seqB->size();
  if (_endA <= _seqA->size()){ _endB = _seqB->size(); }
  else { _endA = _seqA->size(); _endB = _seqA->size() - offset; }
  */
}

AlignmentUngapped* AlignmentUngapped::copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentUngapped::copySeqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentUngapped::copySeqReplace seqsB aren't equal"); }
  return new AlignmentUngapped(seqA, seqB, _orientation, _offset);
}
void AlignmentUngapped::seqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentUngapped::seqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentUngapped::seqReplace seqsB aren't equal"); }
  _seqA = seqA;
  _seqB = seqB;
}
AlignmentUngapped* AlignmentUngapped::copyRcSeqA(ScoredSeq* newSeqA){
  if (newSeqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentNull::seqReplace seqsA aren't equal"); }
  long newOffset = _seqA->size() - _offset - _seqB->size();
  switch (_orientation){
  case '+': return new AlignmentUngapped(newSeqA,_seqB,'-',newOffset);
  case '-': return new AlignmentUngapped(newSeqA,_seqB,'+',newOffset);
  default: throw AssemblyException::StateError("AlUngapped::copyRcSeqA bad sense");
  }
}

ScoredSeq * AlignmentUngapped::seqA(){ return _seqA; }
ScoredSeq* AlignmentUngapped::seqB(){ return _seqB; }
bool AlignmentUngapped::hasConstantOffset(){ return true; }
long AlignmentUngapped::getConstantOffset(){ return _offset; }
char AlignmentUngapped::orientation(){ return _orientation; }
bool AlignmentUngapped::isNull(){ return false; }
bool AlignmentUngapped::isLinked(long pos, ScoredSeq * seq){
  if (seq == _seqA){ return pos >= _startA and pos < _endA; }
  else if (seq == _seqB){ return pos >= _startB and pos < _endB; }
  else { throw AssemblyException::ArgError("AlignmentUngapped::isLinked bad seq in"); }
}
bool AlignmentUngapped::isGapped(long pos, ScoredSeq * seq){ return (! isLinked(pos, seq) ); }
long AlignmentUngapped::getLinkage(long pos, ScoredSeq * seq){
  if (seq == _seqA){ return pos - _offset; }
  else if (seq == _seqB){ return pos + _offset; }
  else { throw AssemblyException::ArgError("AlignmentUngapped::getLinkage bad seq in"); }
}
long AlignmentUngapped::gapPairedAfter(long pos, ScoredSeq * seq){
  long startPos;
  long endAfter;
  if (seq == _seqA){ startPos = _startA; endAfter = _endB - 1; }
  else if (seq == _seqB){ startPos = _startB; endAfter = _endA - 1; }
  else { throw AssemblyException::ArgError("AlignmentUngapped::isLinked bad seq in"); }
  if (pos < startPos){ return -1; }
  else { return endAfter; }
}
long AlignmentUngapped::score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps){
  long misCount = 0;
  long length = _endA - _startA;
  for (long posA = _startA; posA < _endA; ++posA){
    if (! scoreMatrix->isMatch(_seqA->nucAtPosPlus(posA),
			       _seqB->nucAtPosition(posA - _offset,_orientation)) ){ ++misCount; }
  }
  if (penalizeTerminalGaps){
    long overhang = _startA + _startB + _seqA->size() - _endA + _seqB->size() - _endB;
    length += overhang;
    misCount += overhang;
  }
  return (misCount * scoreMatrix->_mismatch) + ( (length - misCount) * scoreMatrix->_match );
}
long AlignmentUngapped::scoreOverhangA(AlignmentScoreMatrix* scoreMatrix){
  long misCount = 0;
  long length = _seqA->size();
  for (long posA = _startA; posA < _endA; ++posA){
    if (! scoreMatrix->isMatch(_seqA->nucAtPosPlus(posA),
			       _seqB->nucAtPosition(posA - _offset,_orientation)) ){ ++misCount; }
  }
  misCount += _startA + length - _endA;
  return (misCount * scoreMatrix->_mismatch) + ( (length - misCount) * scoreMatrix->_match );
} 
bool AlignmentUngapped::isLocked(){ return true; }
ScoredSeq* AlignmentUngapped::combine(){
  ScoreCalculator* calc = new ScoreCalcAdd();
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
} 
ScoredSeq* AlignmentUngapped::combine(int denominator){
  ScoreCalculator* calc = new ScoreCalcNormalizeA(denominator);
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
}
ScoredSeq* AlignmentUngapped::combine(ScoreCalculator* calc){
  long sizeA = _seqA->size();
  long sizeB = _seqB->size();
  // alignment length    5p overhang           3p overhang
  // (_endA - _startA) + (_startA + _startB) + (sizeA - _endA + sizeB - _endB)
  // collapsed:
  long alignSize = _startB + sizeA + sizeB - _endB;
  char* newSeq = new char[ alignSize + 1 ];
  float* newScores = new float[ alignSize + 1 ];
  float* newLinks = new float[ alignSize + 1 ];
  newSeq[ alignSize ] = '\0';

  if (alignSize > sizeA + sizeB or alignSize < sizeA or alignSize < sizeB){
    throw AssemblyException::LogicError("AlUngapped::combine has an invalid alignment length given sequence lengths");
  }

  char* nucsA = _seqA->getSeq('+');
  float* scoresA = _seqA->getScores('+');
  float* linksA = _seqA->getLinks('+');

  char* nucsB = _seqB->getSeq(_orientation);
  float* scoresB = _seqB->getScores(_orientation);
  float* linksB = _seqB->getLinks(_orientation);

  // Fill in the leading overhang
  long posAlign;
  long sizeAm1 = sizeA - 1;
  long sizeBm1 = sizeB - 1;
  if (_startB == 0){
    if (_startA != 0){
      for (long alPos = 0; alPos < _startA; ++alPos){
	newSeq[alPos] = nucsA[alPos];
	newScores[alPos] = calc->adjustA( scoresA[alPos] );
      }
      if (sizeA != 1){
	long numLinks;
	if (_startA == sizeA){ numLinks = _startA - 1; }
	else { numLinks = _startA; }
	for (long alPos = 0; alPos < numLinks; ++alPos){ newLinks[alPos] = calc->adjustA( linksA[alPos] ); }
      }
      if (_startA == sizeA){ newLinks[sizeAm1] = 0; }
    }
    posAlign = _startA;
  } else {
    if (_startB != 0){
      for (long alPos = 0; alPos < _startB; ++alPos){
	newSeq[alPos] = nucsB[alPos];
	newScores[alPos] = calc->adjustA( scoresB[alPos] );
      }
      if (sizeB != 1){
	long numLinks;
	if (_startB == sizeB){ numLinks = _startB - 1; }
	else { numLinks = _startB; }
	for (long alPos = 0; alPos < numLinks; ++alPos){ newLinks[alPos] = calc->adjustA( linksB[alPos] ); }
      }
      if (_startB == sizeB){ newLinks[sizeBm1] = 0; }
    }
    posAlign = _startB;
  }

  // Fill in the aligned portion
  // update posB and posAlign at the end of the loop
  long posB = _startB;
  for (long posA = _startA; posA < _endA; ++posA){
    // decide about the nuc
    char nucA = nucsA[posA];
    char nucB = nucsB[posB];
    float scoreA = scoresA[posA];
    float scoreB = scoresB[posB];
    if ( nucA == nucB ){
      newSeq[posAlign] = nucA;
      newScores[posAlign] = calc->aPlusB(scoreA,scoreB);
    } else if ( calc->adjustA(scoreA) >= calc->adjustB(scoreB) ){ // A is better or the same
      if ( nucA == 'N' ){
	// this is recovery from rounding error to a negative value
	newSeq[posAlign] = nucB;
	newScores[posAlign] = 0;
      } else {
	newSeq[posAlign] = nucA;
	newScores[posAlign] = calc->aMinusB(scoreA,scoreB);
      }
    } else {
      newSeq[posAlign] = nucB;
      newScores[posAlign] = calc->bMinusA(scoreB,scoreA);
    }
    // the links just get added no matter what
    if (posA + 1 < _endA){
      newLinks[posAlign] = calc->aPlusB(linksA[posA], linksB[posB]);
    }
    ++posB;
    ++posAlign;
  }

  // Fill in the trailing overhang
  if (_endB == sizeB){
    for (long pos = _endA; pos < sizeA; ++pos){
      newSeq[posAlign] = nucsA[pos];
      newScores[posAlign] = calc->adjustA( scoresA[pos] );
      if (pos > 0){ newLinks[posAlign-1] = calc->adjustA( linksA[pos-1] ); }
      posAlign++;
    }
  } else {
    for (long pos = _endB; pos < sizeB; ++pos){
      newSeq[posAlign] = nucsB[pos];
      newScores[posAlign] = calc->adjustB( scoresB[pos] );
      if (pos > 0){ newLinks[posAlign-1] = calc->adjustB( linksB[pos-1] ); }
      posAlign++;
    }
  }

  delete [] nucsA;
  delete [] scoresA;
  delete [] linksA;
  delete [] nucsB;
  delete [] scoresB;
  delete [] linksB;

  ScoredSeq * returnVal = ScoredSeq::repExposedSeq(newSeq, newScores, newLinks, alignSize);
  return returnVal;
}


ScoredSeq* AlignmentUngapped::multiCombine(long numAls, Alignment** alArray){

  // the full scope of the multi-alignment must be known from the beginning
  // initial values are from this
  long minOffset;
  if (_offset < 0){ minOffset = _offset; }
  else { minOffset = 0; }
  long maxCoordA;
  if (_seqA->size() >= _seqB->size() + _offset){ maxCoordA = _seqA->size(); }
  else { maxCoordA = _seqB->size() + _offset; }
  // but the scope may have to be expanded based on the other alignments
  for (long alN = 0; alN < numAls; ++alN){
    // checks for the method requirements
    if (! alArray[alN]->hasConstantOffset()){
      throw AssemblyException::ArgError("AlUngapped::multiCombine was given an alignment w/out a constant offset.");
    }
    if (alArray[alN]->seqA() != _seqA){
      throw AssemblyException::ArgError("AlFull:multiCombine requires that all alignments have the same seqA");
    }
    long otherOffset = alArray[alN]->getConstantOffset();
    // adjust the alignment scopes if appropriate
    if (otherOffset < minOffset){ minOffset = otherOffset; }
    if (alArray[alN]->seqB()->size() + otherOffset > maxCoordA){ maxCoordA = alArray[alN]->seqB()->size() + otherOffset; }
  }

  // set up the final sequence (the rest of the function will modify it
  long totalLen = maxCoordA - minOffset;
  char* finalSeq = new char[totalLen+1];
  finalSeq[totalLen] = '\0';
  float* finalScores = new float[totalLen+1];
  float* finalLinks = new float[totalLen+1];

  // the initial values are the ones from seqA; positions that seqA does not include get null Ns
  long sizeAm1 = _seqA->size() - 1;
  long startA = 0 - minOffset;
  char* tempSeqA = _seqA->getSeq('+');
  float* tempScoresA = _seqA->getScores('+');
  float* tempLinksA = _seqA->getLinks('+');
  long fN = 0;
  // initial null
  while (fN < startA){
    finalSeq[fN] = 'N';
    finalScores[fN] = 0;
    finalLinks[fN] = 0;
    ++fN;
  }
  // real values including link scores
  for (long aN = 0; aN < sizeAm1; ++aN){
    finalSeq[fN] = tempSeqA[aN];
    finalScores[fN] = tempScoresA[aN];
    finalLinks[fN] = tempLinksA[aN];
    ++fN;
  }
  // no link at this position
  if (sizeAm1 >= 0){
    finalSeq[fN] = tempSeqA[sizeAm1];
    finalScores[fN] = tempScoresA[sizeAm1];
    finalLinks[fN] = 0;
    ++fN;
  }
  // trailing null (i will include a zero link score even though it is superfluous
  while (fN < totalLen){
    finalSeq[fN] = 'N';
    finalScores[fN] = 0;
    finalLinks[fN] = 0;
    ++fN;
  }
  delete [] tempSeqA;
  delete [] tempScoresA;
  delete [] tempLinksA;

  // go through all of the input als; "this" will be the final al collapsed
  for (long alN = 0; alN <= numAls; ++alN){

    Alignment* localAl;
    if (alN == numAls){ localAl = this; }
    else { localAl = alArray[alN]; }

    ScoredSeq* localSeqB = localAl->seqB();
    char localOrientation = localAl->orientation();

    char* nucsB = localSeqB->getSeq(localOrientation);
    float* scoresB = localSeqB->getScores(localOrientation);
    float* linksB = localSeqB->getLinks(localOrientation);
    long sizeB = localSeqB->size();
    long sizeBm1 = sizeB - 1;

    // the starting position for seqB aligned to the full consensus
    long fN = localAl->getConstantOffset() - minOffset;
    for (long bN = 0; bN < sizeB; ++bN){
      // decide about the nuc
      char nucA = finalSeq[fN];
      char nucB = nucsB[bN];
      float scoreA = finalScores[fN];
      float scoreB = scoresB[bN];
      if ( nucA == nucB ){ finalScores[fN] += scoreB; }
      else if ( scoreA >= scoreB ){ // A is better or the same
	if ( nucA == 'N' ){
	  // this is recovery from rounding error to a negative value
	  finalSeq[fN] = nucB;
	  finalScores[fN] = 0;
	} else { finalScores[fN] -= scoreB; }
      } else {
	finalSeq[fN] = nucB;
	finalScores[fN] = scoreB - scoreA;
      }
      // the links just get added no matter what
      if (bN < sizeBm1){ finalLinks[fN] += linksB[bN]; }
      ++fN;
    }
    delete [] nucsB;
    delete [] scoresB;
    delete [] linksB;
  }

  ScoredSeq * returnVal = ScoredSeq::repExposedSeq(finalSeq, finalScores, finalLinks, totalLen);
  return returnVal;
}

void AlignmentUngapped::OK(){}







#endif
