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



#ifndef ALIGNMENTFULL_CPP
#define ALIGNMENTFULL_CPP

#include "AlignmentFull.h"

#include "AssemblyException.h"
#include <iostream>
#include <typeinfo>
using namespace std;




AlignmentFull::~AlignmentFull(){}
AlignmentFull::AlignmentFull(ScoredSeq * seqA, ScoredSeq* seqB, char orientation) :
  _seqA(seqA), _seqB(seqB), _orientation(orientation)
{ 
  if (_seqA->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentNull::constructor seqs aren't equal"); }
}

AlignmentFull* AlignmentFull::copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentFull::copySeqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentFull::copySeqReplace seqsB aren't equal"); }
  return new AlignmentFull(seqA, seqB, _orientation);
}
void AlignmentFull::seqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentFull::seqReplace seqsA aren't equal"); }
  if (seqB->size() != _seqB->size()){ throw AssemblyException::ArgError("AlignmentFull::seqReplace seqsB aren't equal"); }
  _seqA = seqA;
  _seqB = seqB;
}
AlignmentFull* AlignmentFull::copyRcSeqA(ScoredSeq* newSeqA){
  if (newSeqA->size() != _seqA->size()){ throw AssemblyException::ArgError("AlignmentFull::copyRcSeqA seqsA aren't equal"); }
  switch (_orientation){
  case '+': return new AlignmentFull(newSeqA,_seqB,'-');
  case '-': return new AlignmentFull(newSeqA,_seqB,'+');
  default: throw AssemblyException::StateError("AlFull::copyRcSeqA bad sense");
  }
}

ScoredSeq * AlignmentFull::seqA(){ return _seqA; }
ScoredSeq* AlignmentFull::seqB(){ return _seqB; }
char AlignmentFull::orientation(){ return _orientation; }
bool AlignmentFull::isNull(){ return false; }
bool AlignmentFull::isLinked(long pos, ScoredSeq * seq){ return true; }
bool AlignmentFull::isGapped(long pos, ScoredSeq * seq){ return false; }
long AlignmentFull::getLinkage(long pos, ScoredSeq * seq){ return pos; }
long AlignmentFull::gapPairedAfter(long pos, ScoredSeq * seq){
  throw AssemblyException::CallingError("AlignmentFull has no gaps (gapPairedAfter)");
}
bool AlignmentFull::hasConstantOffset(){ return true; }
long AlignmentFull::getConstantOffset(){ return 0; }

long AlignmentFull::score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps){
  return scoreOverhangA(scoreMatrix);
}
long AlignmentFull::scoreOverhangA(AlignmentScoreMatrix* scoreMatrix){
  long misCount = 0;
  long length = _seqA->size();
  for (long pos = 0; pos < length; ++pos){
    if (! scoreMatrix->isMatch(_seqA->nucAtPosPlus(pos), _seqB->nucAtPosition(pos,_orientation)) ){ ++misCount; }
  }
  return (misCount * scoreMatrix->_mismatch) + ( (length - misCount) * scoreMatrix->_match );
} 
bool AlignmentFull::isLocked(){ return true; }
ScoredSeq* AlignmentFull::combine(){
  ScoreCalculator* calc = new ScoreCalcAdd();
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
} 
ScoredSeq* AlignmentFull::combine(int denominator){
  ScoreCalculator* calc = new ScoreCalcNormalizeA(denominator);
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
}
ScoredSeq* AlignmentFull::combine(ScoreCalculator* calc){
  long alignSize = _seqA->size();
  char* newSeq = new char[ alignSize + 1 ];
  float* newScores = new float[ alignSize + 1 ];
  float* newLinks = new float[ alignSize + 1 ];
  newSeq[ alignSize ] = '\0';

  char* nucsA = _seqA->getSeq('+');
  float* scoresA = _seqA->getScores('+');
  float* linksA = _seqA->getLinks('+');

  char* nucsB = _seqB->getSeq(_orientation);
  float* scoresB = _seqB->getScores(_orientation);
  float* linksB = _seqB->getLinks(_orientation);

  for (long pos = 0; pos < alignSize; ++pos){
    // decide about the nuc
    char nucA = nucsA[pos];
    char nucB = nucsB[pos];
    float scoreA = scoresA[pos];
    float scoreB = scoresB[pos];
    if ( nucA == nucB ){
      newSeq[pos] = nucA;
      newScores[pos] = calc->aPlusB(scoreA,scoreB);
    } else if ( calc->adjustA(scoreA) >= calc->adjustB(scoreB) ){ // A is better or the same
      if ( nucA == 'N' ){
	newSeq[pos] = nucB;
	newScores[pos] = 0;
      } else {
	newSeq[pos] = nucA;
	newScores[pos] = calc->aMinusB(scoreA,scoreB);
      }
    } else {
      newSeq[pos] = nucB;
      newScores[pos] = calc->bMinusA(scoreB,scoreA);
    }
    // the links just get added no matter what
    /*
    if (pos + 1 < alignSize){
      newLinks[pos] = calc->aPlusB(linksA[pos], linksB[pos]);
    }
    */
  }
  long alSizeM1 = alignSize - 1;
  for (long pos = 0; pos < alSizeM1; ++pos){
    newLinks[pos] = calc->aPlusB(linksA[pos], linksB[pos]);
  }
  delete [] nucsA;
  delete [] scoresA;
  delete [] linksA;
  delete [] nucsB;
  delete [] scoresB;
  delete [] linksB;

  //ScoredSeq * returnVal = ScoredSeq::getScoredSeq(newSeq, newScores, newLinks, alignSize);
  ScoredSeq * returnVal = ScoredSeq::repExposedSeq(newSeq, newScores, newLinks, alignSize);
  /*
  delete [] newSeq;
  delete [] newScores;
  delete [] newLinks;
  */
  return returnVal;
}


ScoredSeq* AlignmentFull::multiCombine(long numAls, AlignmentFull** alArray){
  long alignSize = _seqA->size();

  // initially, these have the same values as seqA
  char* newSeq = _seqA->getSeq('+');
  float* newScores = _seqA->getScores('+');
  float* newLinks = _seqA->getLinks('+');

  // go through all of the input als; "this" will be the final al collapsed
  for (long alN = 0; alN <= numAls; ++alN){
    char* nucsB;
    float* scoresB;
    float* linksB;
    if (alN == numAls){
      nucsB = _seqB->getSeq(_orientation);
      scoresB = _seqB->getScores(_orientation);
      linksB = _seqB->getLinks(_orientation);
    } else {
      if (alArray[alN]->_seqA != _seqA){
	throw AssemblyException::ArgError("AlFull:multiCombine requires that all alignments have the same seqA");
      }
      nucsB = alArray[alN]->_seqB->getSeq(_orientation);
      scoresB = alArray[alN]->_seqB->getScores(_orientation);
      linksB = alArray[alN]->_seqB->getLinks(_orientation);
    }
    for (long pos = 0; pos < alignSize; ++pos){
      // decide about the nuc
      char nucA = newSeq[pos];
      char nucB = nucsB[pos];
      float scoreA = newScores[pos];
      float scoreB = scoresB[pos];
      // don't need to change A identity
      if ( nucA == nucB ){
	newScores[pos] += scoreB;
      } else if ( scoreA >= scoreB ){ // A is better or the same
	if ( nucA == 'N' ){
	  newSeq[pos] = nucB;
	  newScores[pos] = 0;
	} else {
	  newScores[pos] -= scoreB;
	}
      } else {
	newSeq[pos] = nucB;
	newScores[pos] = scoreB - scoreA;
      }
      // the links just get added no matter what
      if (pos + 1 < alignSize){
      }
    }
    long alSizeM1 = alignSize - 1;
    for (long pos = 0; pos < alSizeM1; ++pos){  newLinks[pos] += linksB[pos]; }

    delete [] nucsB;
    delete [] scoresB;
    delete [] linksB;
  }

  //ScoredSeq * returnVal = ScoredSeq::getScoredSeq(newSeq, newScores, newLinks, alignSize);
  ScoredSeq * returnVal = ScoredSeq::repExposedSeq(newSeq, newScores, newLinks, alignSize);
  /*
  delete [] newSeq;
  delete [] newScores;
  delete [] newLinks;
  */
  return returnVal;
}


void AlignmentFull::OK(){}







#endif
