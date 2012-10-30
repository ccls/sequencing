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

/*
REP INVARIANT:
The whole time:
- _aMinLink==len(seqA) if there are no links
- _aMaxLink==-1 if there are no links
- if there are any links, _aMinLink <= _aMaxLink
- _isNull if there are no links or gaps
For a non-null alignment (when alignment is locked):
- every index of _aToB and _aGaps must be filled in from _aMin to _aMax, inclusively
- _aMin <= _aMax, NO MATTER WHAT -> UNLESS ALIGNMENT IS NULL, in which case 
  _aMax is -5 (below any possible value) and _aMin is len(seqA) + 5 (above any possible value)

- any linked position must have a value of -5 in _aGaps, and any gapped position in
  the range _aMin <= x <= _aMax must have the value -5 in _aToB (the link array)
- all linked positions must be within the range bounded by _aMin and _aMax
- only terminal gaps may be outside the range bounded by _aMin and _aMax, and
  those are not recorded in the array
- _aMin must be the position immediately following a terminal gap, and _aMax must
  be the position immediately preceeding it
- all of the above also applies to b coords
*/

#ifndef ALIGNMENTGAPPED_CPP
#define ALIGNMENTGAPPED_CPP

#include "AlignmentGapped.h"

#include "AssemblyException.h"
#include <iostream>
#include <typeinfo>
#include <set>
using namespace std;





ScoredSeq* AlignmentGapped::combine(){
  ScoreCalculator* calc = new ScoreCalcAdd();
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
}

ScoredSeq* AlignmentGapped::combine(int denominator){
  ScoreCalculator* calc = new ScoreCalcNormalizeA(denominator);
  ScoredSeq* outputSeq =  combine(calc);
  delete calc;
  return outputSeq;
}
bool AlignmentGapped::hasConstantOffset(){
  if (_isNull){ throw AssemblyException::CallingError("A null AlignmentGapped has no constant offset (hasConstantOffset)"); }
  return false;
}
long AlignmentGapped::getConstantOffset(){
  throw AssemblyException::CallingError("AlignmentGapped has no constant offset (getConstantOffset)");
}



AlignmentGapped::AlignmentGapped(ScoredSeq * seqA, ScoredSeq* seqB, char orientation) :
  _locked(false),
  _seqA(seqA),
  _seqB(seqB),
  _orientation(orientation) {

  _sizeA = _seqA->size();
  _sizeB = _seqB->size();
  _sizeAm1 = _sizeA - 1;
  _sizeBm1 = _sizeB - 1;
  _sizeAp5 = _sizeA + 5;
  _sizeBp5 = _sizeB + 5;

  _aToB = new long[_sizeA + 1];
  _bToA = new long[_sizeB + 1];

  // values are the positions in the other sequence that preceed the gaps
  _aGaps = new long[_sizeA + 1];
  _bGaps = new long[_sizeB + 1];


  if (_orientation != '+' and _orientation != '-'){ throw AssemblyException::ArgError("Alignment: invalid orientation char"); }
  if (_seqA == _seqB){ throw AssemblyException::ArgError("Alignment: two input ScoredSeqs cannot be the same"); }
  _isNull = true; // a holder value

  _aMin = _sizeAp5;
  _aMax = -5;
  _bMin = _sizeBp5;
  _bMax = -5;
  // b min/max links are derived fields based on these
  _aMinLink = _sizeA;
  _aMaxLink = -1;
}


AlignmentGapped* AlignmentGapped::copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  return new AlignmentGapped(this, seqA, seqB);
}
void AlignmentGapped::seqReplace(ScoredSeq* seqA, ScoredSeq* seqB){
  if (seqA->size() != _sizeA){ throw AssemblyException::ArgError("AlignmentGapped::seqReplace seqsA aren't equal"); }
  if (seqB->size() != _sizeB){ throw AssemblyException::ArgError("AlignmentGapped::seqReplace seqsB aren't equal"); }
  _seqA = seqA;
  _seqB = seqB;
}


inline long AlignmentGapped::reversePosCoord(long pos, long seqLenM1){
  if (pos == -5){ return pos; }
  else { return seqLenM1 - pos; }
}
inline long AlignmentGapped::reversePosCoordGap(long pos, long seqLenM1){
  if (pos == -5){ return pos; }
  else { return seqLenM1 - pos - 1; }
}

AlignmentGapped* AlignmentGapped::copyRcSeqA(ScoredSeq* newSeqA){
  //throw AssemblyException::ImplementationError("AlGapped::copyRcSeqA is not implemented for the current rep inv");

  // make the orientation opposite what it was
  if (newSeqA->size() != _sizeA){ throw AssemblyException::ArgError("AlignmentGapped::copyRcSeqA seqsA aren't equal"); }
  AlignmentGapped* newAl;
  switch (_orientation){
  case '+': newAl = new AlignmentGapped(newSeqA,_seqB,'-'); break;
  case '-': newAl = new AlignmentGapped(newSeqA,_seqB,'+'); break;
  default: throw AssemblyException::StateError("AlGapped::copyRcSeqA bad sense");
  }

  newAl->_isNull = _isNull;
  newAl->_aMin = reversePosCoord(_aMax,_sizeAm1);
  newAl->_aMax = reversePosCoord(_aMin,_sizeAm1);
  newAl->_bMin = reversePosCoord(_bMax,_sizeBm1);
  newAl->_bMax = reversePosCoord(_bMin,_sizeBm1);
  newAl->_aMinLink = reversePosCoord(_aMaxLink,_sizeAm1);
  newAl->_aMaxLink = reversePosCoord(_aMinLink,_sizeAm1);

  // for the two updates below, i do not need to check if the sequence
  // is null because if it is, _abMin > _abMax so nothing happens

  // reverse all of the linkage/gap positions
  for (long pos = _aMin; pos <= _aMax; ++pos){
    newAl->_aToB[ reversePosCoord(pos,_sizeAm1) ] = reversePosCoord(_aToB[pos], _sizeBm1);
    newAl->_aGaps[ reversePosCoord(pos,_sizeAm1) ] = reversePosCoordGap(_aGaps[pos], _sizeBm1);
  }
  // for seqB, just cover the gaps, the links have been handled
  for (long pos = _bMin; pos <= _bMax; ++pos){
    newAl->_bToA[ reversePosCoord(pos,_sizeBm1) ] = reversePosCoord(_bToA[pos], _sizeAm1);
    newAl->_bGaps[ reversePosCoord(pos,_sizeBm1) ] = reversePosCoordGap(_bGaps[pos], _sizeAm1);
  }

  newAl->_locked = _locked;
  //newAl->OK();
  return newAl;
}

AlignmentGapped::AlignmentGapped(AlignmentGapped* oldAlignment, ScoredSeq* newSeqA, ScoredSeq* newSeqB) :
  _locked(false),
  _seqA(newSeqA),
  _seqB(newSeqB) {

  _sizeA = _seqA->size();
  _sizeB = _seqB->size();
  _sizeAm1 = _sizeA - 1;
  _sizeBm1 = _sizeB - 1;
  _sizeAp5 = _sizeA + 5;
  _sizeBp5 = _sizeB + 5;

  // copy all of the values from the old alignment
  if (! oldAlignment->_locked ){
    throw AssemblyException::ArgError("Alignment: you shouldn't be copying an unlocked alignment");
  }
  _orientation = oldAlignment->_orientation;
  // _isNull will be changed by lock() at the end of the constructor if appropriate
  _isNull = oldAlignment->_isNull;

  _aMin = oldAlignment->_aMin;
  _aMax = oldAlignment->_aMax;
  _bMin = oldAlignment->_bMin;
  _bMax = oldAlignment->_bMax;
  _aMinLink = oldAlignment->_aMinLink;
  _aMaxLink = oldAlignment->_aMaxLink;

  // check the length requirement
  if ( _sizeA != oldAlignment->_sizeA ) {
    throw AssemblyException::ArgError("seqA sizes don't match in Alignment copying"); }
  if ( _sizeB != oldAlignment->_sizeB ) {
    cerr << endl;
    cerr << _sizeB << endl;
    cerr << oldAlignment->_sizeB << endl;
    throw AssemblyException::ArgError("seqB sizes don't match in Alignment copying"); }
  // copy the arrays
  _aToB = new long[_sizeA + 1];
  _aGaps = new long[_sizeA + 1];
  _bToA = new long[_sizeB + 1];
  _bGaps = new long[_sizeB + 1];
  // for the two updates below, i do not need to check if the sequence
  // is null because if it is, _abMin > _abMax so nothing happens
  for (long n = _aMin; n <= _aMax; ++n){
    _aToB[n] = oldAlignment->_aToB[n];
    _aGaps[n] = oldAlignment->_aGaps[n];
  }
  for (long n = _bMin; n <= _bMax; ++n){
    _bToA[n] = oldAlignment->_bToA[n];
    _bGaps[n] = oldAlignment->_bGaps[n];
  }
  if (oldAlignment->_locked){ lock(); }
  //OK();
}


AlignmentGapped::~AlignmentGapped(){
  //OK();
  delete [] _aToB;
  delete [] _bToA;
  delete [] _aGaps;
  delete [] _bGaps;
}


inline void AlignmentGapped::helpAddLinkSeqA(long posSeq, long posPair){
  _aToB[posSeq] = posPair;
  _bToA[posPair] = posSeq;
  _aGaps[posSeq] = -5;
  _bGaps[posPair] = -5;
  updateMinMaxA(posSeq);
  updateMinMaxB(posPair);
}
inline void AlignmentGapped::helpAddLinkSeqB(long posSeq, long posPair){
  _aToB[posPair] = posSeq;
  _bToA[posSeq] = posPair;
  _aGaps[posPair] = -5;
  _bGaps[posSeq] = -5;
  updateMinMaxA(posPair);
  updateMinMaxB(posSeq);
}
inline void AlignmentGapped::updateMinMaxA(long pos){
  if (pos < _aMin){ _aMin = pos; }
  if (pos > _aMax){ _aMax = pos; }
}
inline void AlignmentGapped::updateMinMaxB(long pos){
  if (pos < _bMin){ _bMin = pos; }
  if (pos > _bMax){ _bMax = pos; }
}

void AlignmentGapped::addLinkage(ScoredSeq * seq, long posSeq, long posPair){
  if ( _locked ){ throw AssemblyException::ArgError("mutator method cannot be called because object is locked."); }
  _isNull = false;
  if ( seq == _seqA ){
    helpAddLinkSeqA(posSeq,posPair);
    if (_aMinLink > posSeq){ _aMinLink = posSeq; }
    if (_aMaxLink < posSeq){ _aMaxLink = posSeq; }
  } else if ( seq == _seqB ){
    helpAddLinkSeqB(posSeq,posPair);
    if (_aMinLink > posPair){ _aMinLink = posPair; }
    if (_aMaxLink < posPair){ _aMaxLink = posPair; }
  }
  else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}

void AlignmentGapped::addLinkageBlock(ScoredSeq * seq, long posSeq, long posPair, long blockLength){
  _isNull = false;
  if ( seq == _seqA ){
    for (long n = 0; n < blockLength; n++){ helpAddLinkSeqA(posSeq + n, posPair + n); }
    if (_aMinLink > posSeq){ _aMinLink = posSeq; }
    if (_aMaxLink < posSeq + blockLength - 1){ _aMaxLink = posSeq + blockLength - 1; }
  } else if ( seq == _seqB ){
    for (long n = 0; n < blockLength; n++){ helpAddLinkSeqB(posSeq + n, posPair + n); }
    if (_aMinLink > posPair){ _aMinLink = posPair; }
    if (_aMaxLink < posPair + blockLength - 1){ _aMaxLink = posPair + blockLength - 1; }
  } else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}


inline void AlignmentGapped::addGapSeqA(long posIns, long posGap){
  addGapNoUpdateSeqA(posIns,posGap);
  updateMinMaxA(posIns);
}
inline void AlignmentGapped::addGapSeqB(long posIns, long posGap){
  addGapNoUpdateSeqB(posIns,posGap);
  updateMinMaxB(posIns);
}



inline void AlignmentGapped::addGapNoUpdateSeqA(long posIns, long posGap){
  _isNull = false;
  if (posGap >= _sizeB){ throw AssemblyException::ArgError("AlGapped:addGapSeqA, highest gap pos after can be is seq length - 1."); }
  _aGaps[posIns] = posGap;
  _aToB[posIns] = -5;
}
inline void AlignmentGapped::addGapNoUpdateSeqB(long posIns, long posGap){
  _isNull = false;
  if (posGap >= _sizeA){ throw AssemblyException::ArgError("AlGapped:addGapSeqB, highest gap pos after can be is seq length - 1."); }
  _bGaps[posIns] = posGap;
  _bToA[posIns] = -5;
}

void AlignmentGapped::addGap(ScoredSeq* inSeq, long posIns, long posGap){
  if ( _locked ){ throw AssemblyException::ArgError("mutator method cannot be called because object is locked."); }
  if (posGap < -1){ throw AssemblyException::ArgError("lowest gap pos after can be is -1."); }
  _isNull = false;
  if ( inSeq == _seqA ){
    addGapSeqA(posIns,posGap);
  } else if ( inSeq == _seqB ){
    addGapSeqB(posIns,posGap);
  } else {
    throw AssemblyException::ArgError("AlGapped:addGap, provided seq is not part of this alignment.");
  }
}

void AlignmentGapped::helpAddGapBlockSeqA(long posIns, long posGap, long blockLength){
  long posInsEnd = posIns + blockLength;
  if (posGap == -1){
    // first special case: only the last nuc of the block needs to be marked
    // since it is the end of the leading gap
    addGapSeqA(posIns + blockLength - 1, posGap);
  } else if (posInsEnd == _sizeA and posGap == _sizeBm1){
    // second special case: only the first position of the block needs to
    // marked since the rest of the alignment is a trailing gap
    addGapSeqA(posIns, posGap);
  } else {
    for (long n = posIns; n < posInsEnd; ++n){ addGapNoUpdateSeqA(n, posGap); }
    if (posIns < _aMin){ _aMin = posIns; }
    if (posInsEnd - 1 > _aMax){ _aMax = posInsEnd - 1; }
  }
}
void AlignmentGapped::helpAddGapBlockSeqB(long posIns, long posGap, long blockLength){
  long posInsEnd = posIns + blockLength;
  if (posGap == -1){
    // first special case: only the last nuc of the block needs to be marked
    // since it is the end of the leading gap
    addGapSeqB(posIns + blockLength - 1, posGap);
  } else if (posInsEnd == _sizeB and posGap == _sizeAm1){
    // second special case: only the first position of the block needs to
    // marked since the rest of the alignment is a trailing gap
    addGapSeqB(posIns, posGap);
  } else {
    for (long n = posIns; n < posInsEnd; ++n){ addGapNoUpdateSeqB(n, posGap); }
    if (posIns < _bMin){ _bMin = posIns; }
    if (posInsEnd - 1 > _bMax){ _bMax = posInsEnd - 1; }
  }
}

void AlignmentGapped::addGapBlock(ScoredSeq* inSeq, long posIns, long posGap, long blockLength){
  _isNull = false;
  //OK();
  if ( inSeq == _seqA ){
    helpAddGapBlockSeqA(posIns,posGap,blockLength);
    // implied leading
    if (posIns == 0 and posGap != -1){ helpAddGapBlockSeqB(0,-1,posGap+1); }
    // implied trailing
    if (posIns + blockLength == _sizeA and posGap != _sizeBm1){
      helpAddGapBlockSeqB(posGap+1, _sizeAm1, _sizeBm1 - posGap);
    }
  } else if ( inSeq == _seqB ){
    helpAddGapBlockSeqB(posIns,posGap,blockLength);
    // now worry about implied blocks
    // implied leading
    if (posIns == 0 and posGap != -1){ helpAddGapBlockSeqA(0, -1, posGap + 1); }
    if (posIns + blockLength == _sizeB and posGap != _sizeAm1){
      long bl = _sizeAm1 - posGap;
      helpAddGapBlockSeqA(posGap+1, _sizeBm1, _sizeAm1 - posGap);
    }
  } else {
    throw AssemblyException::ArgError("provided seq is not part of this alignment.");
  }
}

ScoredSeq * AlignmentGapped::seqA(){ return _seqA; }
ScoredSeq* AlignmentGapped::seqB(){ return _seqB; }
char AlignmentGapped::orientation(){ return _orientation; }

bool AlignmentGapped::isNull(){ 
  if (! _locked ){ throw AssemblyException::CallingError("cannot call 'isNull' if alignment is not locked."); }
  return _isNull;
}

bool AlignmentGapped::isLinked(long pos, ScoredSeq * seq){
  if (seq == _seqA) {
    if ( pos > _sizeAm1 ){ throw AssemblyException::ArgError("AlignmentGapped::isLinked, pos was outside of seqA"); }
    return (pos >= _aMin and pos <= _aMax and _aToB[pos] != -5);
  } else if (seq == _seqB) {
    if ( pos > _sizeBm1 ){ throw AssemblyException::ArgError("AlignmentGapped::isLinked, pos was outside of seqB"); }
    return (pos >= _bMin and pos <= _bMax and _bToA[pos] != -5);
  }
  else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}
inline bool AlignmentGapped::isLinkedSeqA(long pos){ return (pos >= _aMin and pos <= _aMax and _aToB[pos] != -5); }
inline bool AlignmentGapped::isLinkedSeqB(long pos){ return (pos >= _bMin and pos <= _bMax and _bToA[pos] != -5); }

bool AlignmentGapped::isGapped(long pos, ScoredSeq * seq){
  if (seq == _seqA) {
    if ( pos > _sizeAm1 ){ throw AssemblyException::ArgError("AlignmentGapped::isGapped, pos was outside of seqA"); }
    return (pos < _aMin or pos > _aMax or _aGaps[pos] != -5);
  }
  else if (seq == _seqB) {
    if ( pos > _sizeBm1 ){ throw AssemblyException::ArgError("AlignmentGapped::isGapped, pos was outside of seqB"); }
    return (pos < _bMin or pos > _bMax or _bGaps[pos] != -5);
  }
  else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}
inline bool AlignmentGapped::isGappedSeqA(long pos){ return (pos < _aMin or pos > _aMax or _aGaps[pos] != -5); }
inline bool AlignmentGapped::isGappedSeqB(long pos){ return (pos < _bMin or pos > _bMax or _bGaps[pos] != -5); }


long AlignmentGapped::getLinkage(long pos, ScoredSeq * seq){
  if (! isLinked(pos,seq) ){
    throw AssemblyException::ArgError("provided pos is not linked.");
  }
  if (seq == _seqA) { 
    if (! isLinkedSeqA(pos) ){ throw AssemblyException::ArgError("AlignmentGapped::getLink - provided seqA pos is not linked."); }
    if (pos < _aMin){ return -1; }
    else if (pos > _aMax){ return _seqB->size()-1; }
    else { return _aToB[pos]; }
  } else if (seq == _seqB) { 
    if (! isLinkedSeqB(pos) ){ throw AssemblyException::ArgError("AlignmentGapped::getLink - provided seqA pos is not linked."); }
    if (pos < _bMin){ return -1; }
    else if (pos > _bMax){ return _seqA->size()-1; }
    else { return _bToA[pos]; }
  } else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}
inline long AlignmentGapped::getLinkageSeqA(long pos){ return _aToB[pos]; }
inline long AlignmentGapped::getLinkageSeqB(long pos){ return _bToA[pos]; }


long AlignmentGapped::gapPairedAfter(long pos, ScoredSeq * seq){
  if (seq == _seqA) {
    if (! isGappedSeqA(pos)){ throw AssemblyException::ArgError("AlignmentGapped, cant get gap pair for an ungapped position seqA."); }
    if (pos < _aMin){ return -1; }
    else if (pos > _aMax){ return _sizeBm1; }
    else { return _aGaps[pos]; }
  } else if (seq == _seqB) {
    if (! isGappedSeqB(pos)){ throw AssemblyException::ArgError("AlignmentGapped, cant get gap pair for an ungapped position seqB."); }
    if (pos < _bMin){ return -1; }
    else if (pos > _bMax){ return _sizeAm1; }
    else { return _bGaps[pos]; }
  } else { throw AssemblyException::ArgError("provided seq is not part of this alignment."); }
}
inline long AlignmentGapped::gapPairedAfterSeqA(long pos){ 
  if (pos < _aMin){ return -1; }
  else if (pos > _aMax){ return _sizeBm1; }
  else { return _aGaps[pos]; }
}
inline long AlignmentGapped::gapPairedAfterSeqB(long pos){ 
  if (pos < _bMin){ return -1; }
  else if (pos > _bMax){ return _sizeAm1; }
  else { return _bGaps[pos]; }
}


long AlignmentGapped::scoreOverhangA(AlignmentScoreMatrix* scoreMatrix){
  // first, get the score for the aligned segment
  long mainScore = score(scoreMatrix, false);
  // now measure the seqA overhangs
  mainScore += _aMin * scoreMatrix->_extendGap;
  mainScore += (_sizeAm1 - _aMax) * scoreMatrix->_extendGap;
  return mainScore;
}

long AlignmentGapped::score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps){
  //throw AssemblyException::ImplementationError("ALGapped::score not implemented for current rep inv");

  //OK();
  if (! _locked){
    cerr << int(_locked) << endl;
    throw AssemblyException::CallingError("cannot call 'calculateScore' if alignment is not locked.");
  }
  if ( _isNull ){
    cerr << int(_locked) << endl;
    throw AssemblyException::CallingError("cannot call 'calculateScore' if alignment is null.");
  }

  // the tally variable
  long scoreSum = 0;

  // establish where the beginning and end of each alignment are
  long posA = _aMin;
  long posB = _bMin;
  long sizeA = _aMax + 1;
  long sizeB = _bMax + 1;

  // make the decisions for the aligned segment
  while ( posA < sizeA or posB < sizeB ){
    if ( posA==sizeA or isGappedSeqA(posA) or posB==sizeB or isGappedSeqB(posB) ){
      long aStartPos = posA;    // unlike the combine method, here i start with
      long bStartPos = posB;    // the gapped position
      // bring both positions to the point after the alignment block
      while ( posA < sizeA and isGappedSeqA(posA) ){ ++posA; }
      while ( posB < sizeB and isGappedSeqB(posB) ){ ++posB; }

      // check that block is legit
      if ( posA!=sizeA and isLinkedSeqA(posA) and getLinkageSeqA(posA)==posB ){} // this is legit
      else if ( posA==sizeA and posB==sizeB ){} // this is also legit
      else {
	cerr << endl << "ERR: " << posA << " " << posB << " " << _aMin << " " << _bMin << endl;
	cerr << _seqA->getSeq('+') << endl;
	cerr << _seqB->getSeq(_orientation) << endl;
	throw AssemblyException::LogicError("AlignGapped::score, conditions after a gap block are not legit.");
      }
      scoreSum += scoreGapBlock(aStartPos, posA, bStartPos, posB, scoreMatrix, penalizeTerminalGaps);
    } else { // MATCH
      scoreSum += scoreMatrix->match( _seqA->nucAtPosPlus(posA),
				      _seqB->nucAtPosition(posB,_orientation) );
      posA++;
      posB++;
    }
  }

  if (penalizeTerminalGaps){
    if (_aMin > 0){ scoreSum += _aMin * scoreMatrix->_extendGap; }
    if (_aMax + 1 < _seqA->size()){
      scoreSum += (_seqA->size() - _aMax - 1) * scoreMatrix->_extendGap;
    }
    if (_bMin > 0){ scoreSum += _bMin * scoreMatrix->_extendGap; }
    if (_bMax + 1 < _seqB->size()){
      scoreSum += (_seqB->size() - _bMax - 1) * scoreMatrix->_extendGap;
    }
  }
  return scoreSum;
}


// this is a helper for score
long AlignmentGapped::scoreGapBlock(long aStartPos, long aEndPos,
				    long bStartPos, long bEndPos,
				    AlignmentScoreMatrix* scoreMatrix,
				    bool penalizeTerminalGaps){

  long scoreTally = 0;

  // collect the sizes of the gaps that will be penalized by defining gap segments within the block
  GapToSize* aGapToSize = getSegmentsFromBlockSeqA(aStartPos,aEndPos,_sizeBm1,penalizeTerminalGaps);
  for (long n = 0; n < aGapToSize->_gCount; ++n){
    scoreTally += scoreMatrix->_newGap + (aGapToSize->_gLengths[n] - 1) * scoreMatrix->_extendGap;
  }
  delete aGapToSize;

  GapToSize* bGapToSize = getSegmentsFromBlockSeqB(bStartPos,bEndPos,_sizeAm1,penalizeTerminalGaps);
  for (long n = 0; n < bGapToSize->_gCount; ++n){
    scoreTally += scoreMatrix->_newGap + (bGapToSize->_gLengths[n] - 1) * scoreMatrix->_extendGap;
  }
  delete bGapToSize;
  return scoreTally;
}

AlignmentGapped::GapToSize::GapToSize(long maxLen) : _gCount(0), _gIndex(-1){
  _gStarts = new long[maxLen+1];
  _gLengths = new long[maxLen+1];
}
AlignmentGapped::GapToSize::~GapToSize(){
  delete [] _gStarts;
  delete [] _gLengths;
}
void AlignmentGapped::GapToSize::addPosition(long pos){
  _gCount++;
  _gIndex++;
  _gStarts[_gIndex] = pos;
  _gLengths[_gIndex] = 1;
}
void AlignmentGapped::GapToSize::addCountToCurrent(){
  _gLengths[_gIndex]++;
}

// modifies "gapToSize"; this is a helper for scoreGapBlock
AlignmentGapped::GapToSize* AlignmentGapped::getSegmentsFromBlockSeqA(long startPos, long endPos,
								      long maxOpposingPosition, bool penalizeTerminalGaps){
  // the value associated with pos will be changed during the run of this method;
  // startPos will not
  long pos = startPos;
  // counts will be added to this gapToSize bin - is only initialized
  // when a gap is found
  long currentGapPos; 

  // used to make sure that discontinuous gaps are separated;
  // initially, this is set so that it is below the theorietical
  // minimum position, therefore initial gaps will be treated
  // appropriately (as new gaps)
  long priorGapAfter = -10;

  GapToSize* gapToSize = new GapToSize(endPos - startPos + 1);
  while (pos < endPos){
    long currentGapAfter = gapPairedAfterSeqA(pos);
    if ( penalizeTerminalGaps or (currentGapAfter!=-1 and currentGapAfter!=maxOpposingPosition) ){
      // the dummy priorGapAfter value above makes it unneccesary to check if this is the startPos
      if (currentGapAfter!=priorGapAfter ){
	priorGapAfter = currentGapAfter;
	currentGapPos = pos;
	gapToSize->addPosition(currentGapPos);
      } else {
	gapToSize->addCountToCurrent();
      }
    }
    ++pos;
  }
  return gapToSize;
}
// modifies "gapToSize"; this is a helper for scoreGapBlock
AlignmentGapped::GapToSize* AlignmentGapped::getSegmentsFromBlockSeqB(long startPos, long endPos,
								      long maxOpposingPosition, bool penalizeTerminalGaps){
  long pos = startPos;
  long currentGapPos; // counts will be added to this gapToSize bin
  long priorGapAfter; // used to make sure that discontinuous gaps are separated

  GapToSize* gapToSize = new GapToSize(endPos - startPos + 1);
  while (pos < endPos){
    long currentGapAfter = gapPairedAfterSeqB(pos);
    if ( penalizeTerminalGaps or (currentGapAfter!=-1 and currentGapAfter!=maxOpposingPosition) ){
      if ( pos==startPos or currentGapAfter!=priorGapAfter ){
	priorGapAfter = currentGapAfter;
	currentGapPos = pos;
	gapToSize->addPosition(currentGapPos);
      } else {
	gapToSize->addCountToCurrent();
      }
    }
    ++pos;
  }
  return gapToSize;
}




void AlignmentGapped::lock(){
  //OK();
  _isNull = (_aMax==-5 and _bMax==-5);
  if (! _isNull){

    // fill in the 3p terminal gaps (or don't, but make sure they are considered)
    if (_aMax == _sizeAm1 and _bMax == _sizeBm1){
      // fine, the whole thing is filled in
    } else if (_aMax == _sizeAm1) {
      // make sure that the position after the gap in B points to the end of A
      if (_bMax != -5){ addGapSeqB(_bMax + 1, _sizeAm1); }
    } else if (_bMax == _sizeBm1) {
      // make sure that the position after the gap in A points to the end of B
      if (_aMax != -5){ addGapSeqA(_aMax + 1, _sizeBm1); }
    } else {
      throw AssemblyException::LogicError("AlGapped::lock, alignment doesn't go to the 3p end of either sequence.");
    }

    // fill in the 5p terminal gaps (or don't, but make sure they are considered)
    if (_aMin == 0 and _bMin == 0){
      // fine, the whole thing is filled in
    } else if (_aMin == 0) {
      // make sure that the position after the gap in B points to the end of A
      if (_bMax != _sizeBp5){ addGapSeqB(_bMin - 1, -1); }
    } else if (_bMin == 0) {
      // make sure that the position after the gap in A points to the end of B
      if (_aMax != _sizeAp5){ addGapSeqA(_aMin - 1, -1); }
    } else {
      throw AssemblyException::LogicError("AlGapped::lock, alignment doesn't go to the 5p end of either sequence.");
    }
  }

  // lock down and check rep
  _locked = true;
  //OK();
}


bool AlignmentGapped::isLocked(){ return _locked; }

void AlignmentGapped::OK(){
  if ( int(_locked)!=0 and int(_locked)!=1 ){
    cerr << int(_locked) << endl;
    throw AssemblyException::LogicError("in Alignment, _locked is a nonsense value.");
  }
  if (_orientation!='+' and _orientation!='-'){
    cerr << _orientation << endl;
    throw AssemblyException::LogicError("in Alignment, _orientation is a nonsense value.");
  }
  if ( _locked ){
    if (! _isNull ){
      if (_aMin != 0 and _bMin != 0){
	throw AssemblyException::LogicError("AlGapped::OK, either _aMin or _bMin should be 0.");
      }
      if (_aMax != _sizeAm1 and _bMax != _sizeBm1){
	throw AssemblyException::LogicError("AlGapped::OK, either _aMax or _bMax should be the last pos of the sequence.");
      }

      long aLinkCount = 0;
      long aGapCount = 0;
      for (long n = _aMin; n <= _aMax; ++n){
	if (_aToB[n] != -5){ ++aLinkCount; }
	if (_aGaps[n] != -5){ ++aGapCount; }
      }
      long bLinkCount = 0;
      long bGapCount = 0;
      for (long n = _bMin; n <= _bMax; ++n){
	if (_bToA[n] != -5){ ++bLinkCount; }
	if (_bGaps[n] != -5){ ++bGapCount; }
      }
      if ( aLinkCount != bLinkCount ){
	cerr << endl;
	cerr << "Link counts: " << aLinkCount << " " << bLinkCount << endl;
	cerr << endl;
	cerr << "A min/max: " << _aMin << " " << _aMax << endl;
	for (long n = _aMin; n <= _aMax; ++n){ cerr << _aToB[n] << '\t' << _aGaps[n] << endl; }
	cerr << endl;
	cerr << "B min/max: " << _bMin << " " << _bMax << endl;
	for (long n = _bMin; n <= _bMax; ++n){ cerr << _bToA[n] << '\t' << _bGaps[n] << endl; }
	cerr << aLinkCount << " " << bLinkCount << endl;
	throw AssemblyException::LogicError("linkage maps are of unequal size.");
      }
      if ( aLinkCount + aGapCount != _aMax - _aMin + 1 ){
	cerr << "A " << aLinkCount << " " << aGapCount << " " << _aMin << " " << _aMax << " " << _seqA->size() << endl;
	cerr << "B " << bLinkCount << " " << bGapCount << " " << _bMin << " " << _bMax << " " << _seqB->size() << endl;
	cerr << _bGaps[_bMin] << " " << _bGaps[_bMax] << endl;
	throw AssemblyException::LengthError("not all seqA positions are accounted for.");
      }
      if ( bLinkCount + bGapCount != _bMax - _bMin + 1 ){
	throw AssemblyException::LengthError("not all seqB positions are accounted for.");
      }
    }
    for (long aKey = _aMin; aKey < _aMax+1; ++aKey){
      long bKey = _aToB[ aKey ];
      if ( bKey != -5 ){
	if ( _bToA[ bKey ] != aKey ){
	  throw AssemblyException::LogicError("linkage is not reciprocal");
	}
	if ( _aGaps[aKey] != -5 ){
	  throw AssemblyException::LogicError("linked position in A is also gapped");
	}
	if ( _bGaps[bKey] != -5 ){
	  throw AssemblyException::LogicError("linked position in B is also gapped");
	}
      }
    // check that everything flows in the right direction
    }
  }
}

AlignmentGapped::SeqInfoCarrier::SeqInfoCarrier(ScoredSeq* seq, char sense) :
  _seq(seq),
  _sense(sense)
{
  _nucs = seq->getSeq(sense);
  _scores = seq->getScores(sense);
  _links = seq->getLinks(sense);
}
AlignmentGapped::SeqInfoCarrier::~SeqInfoCarrier(){
  delete [] _nucs;
  delete [] _scores;
  delete [] _links;
}


ScoredSeq* AlignmentGapped::combine(ScoreCalculator* calc){
  if (! _locked){
    throw AssemblyException::ArgError("in combine, unlocked alignment should not be combined");
  }

  // these variables are defined just for ease of access
  //char senseA = '+';
  char senseB = _orientation;
  SeqInfoCarrier* infoA = new SeqInfoCarrier(_seqA,'+');
  SeqInfoCarrier* infoB = new SeqInfoCarrier(_seqB,_orientation);

  // the values for the constructor of the new ScoredSeq
  CombineDataCarrier * carrier = new CombineDataCarrier(_seqA->size() + _seqB->size());

  // establish where the beginning and end of each alignment are
  long posA = 0;
  long posB = 0;
  bool lastPosWasMatch = false;

  // make the decisions for the aligned segment
  while ( posA < _sizeA or posB < _sizeB ){

    if ( posA==_sizeA or isGappedSeqA(posA) or posB==_sizeB or isGappedSeqB(posB) ){
      long aStartPos = posA-1;
      long bStartPos = posB-1;
      // bring both positions to the point after the alignment block
      while ( posA < _sizeA and isGappedSeqA(posA) ){ ++posA; }
      while ( posB < _sizeB and isGappedSeqB(posB) ){ ++posB; }

      // check that block is legit
      if ( posA!=_sizeA and isLinkedSeqA(posA) and getLinkageSeqA(posA)==posB ){} // this is legit
      else if ( posA==_sizeA and posB==_sizeB ){} // this is also legit
      else {
	throw AssemblyException::LogicError("in combine, conditions after a gap block are not legit.");
      }

      // if a linkage was already added to the end, it must be removed so it can be replaced
      if (lastPosWasMatch){ carrier->decrementLink(); }
      combineGapBlock(infoA, infoB, aStartPos, posA, bStartPos, posB, carrier, calc);
      lastPosWasMatch = false;

    } else { // MATCH
      char nucA = infoA->_nucs[posA];
      char nucB = infoB->_nucs[posB];
      float scoreA = infoA->_scores[posA];
      float scoreB = infoB->_scores[posB];

      float linkA = 0.0;
      float linkB = 0.0;
      // neither, one or both of these may be the case; a zero value will be added if neither is
      // the case, which is fine since it will not be used for ScoredSeq construction
      if ( posA + 1 < _sizeA ){ linkA = infoA->_links[posA]; }
      if ( posB + 1 < _sizeB ){ linkB = infoB->_links[posB]; }
      float localLinkScore = calc->aPlusB(linkA,linkB);

      if ( nucA == nucB ){
	carrier->addNuc(nucA, calc->aPlusB(scoreA,scoreB), localLinkScore);
      } else if ( calc->adjustA(scoreA) >= calc->adjustB(scoreB) ){ // A is better or the same
	if ( nucA == 'N' ){
	  carrier->addNuc(nucB, 0.0, localLinkScore);
	} else {
	  carrier->addNuc(nucA, calc->aMinusB(scoreA,scoreB), localLinkScore);
	}
      } else { // here, B is better
	carrier->addNuc(nucB, calc->bMinusA(scoreB,scoreA), localLinkScore);
      }

      posA++;
      posB++;
      lastPosWasMatch = true;
    }
  }

  delete infoA;
  delete infoB;

  ScoredSeq * returnVal = carrier->getSeq();
  delete carrier;
  return returnVal;
}



// modifies "links"; this is a helper for combineGapBlock
void AlignmentGapped::getSegmentsFromBlock(map<long,float>* links, long startPos, long endPos, char sense,
				     ScoredSeq* seq, ScoreCalcXy* calc){
  long pos = startPos;
  long seqSize = seq->size();
  while (pos < endPos){
    if (pos==-1 or pos+1==seqSize){ links->insert( pair<long,float> ( pos, calc->adjustX(0) ) ); }
    else if ( pos==startPos or pos+1==endPos or 
	      gapPairedAfter(pos,seq) != gapPairedAfter(pos+1,seq) ){
      links->insert( pair<long,float> ( pos, calc->adjustX( seq->linkAfterPosition(pos,sense) ) ) );
    }
    ++pos;
  }
}


// modifies "carrier"; this is a helper for combineGapBlock
void AlignmentGapped::fillInGaps(CombineDataCarrier* carrier,
				 SeqInfoCarrier* seqInfo, long startPos, long endPos,
				 float penaltyConstant, map<long,float>* links, ScoreCalcXy* calc){
  long seqSize = seqInfo->_seq->size();
  for (long n = startPos; n < endPos; ++n){
    // fill in non-link values
    if ( n != startPos ){
      carrier->addNuc(seqInfo->_nucs[n], calc->adjustX( seqInfo->_scores[n] ));
    }
    // the before and after links are imaginary and don't get added to the carrier
    if (n != -1 and n+1 != seqSize){
      map<long,float>::iterator compScoreIt = links->find(n);
      if (compScoreIt==links->end()){ // this was a non-competitive position
	carrier->addLink( calc->adjustX( seqInfo->_links[n] ) );
      } else { // this was a competitive position, so the value is modified by the constant
	carrier->addLink( compScoreIt->second * penaltyConstant );
      }
    }
  }
}

void AlignmentGapped::combineGapBlock(SeqInfoCarrier* infoA, SeqInfoCarrier* infoB,
				      long aStartPos, long aEndPos,
				      long bStartPos, long bEndPos,
				      CombineDataCarrier* carrier, ScoreCalculator* calc){

  char senseA = '+';
  char senseB = _orientation;


  // collect the relevant (competing) link scores by defining gap segments within the block
  map<long,float> aLinks;
  map<long,float> bLinks;
  ScoreCalcXy * tempCalcA = ScoreCalcXy::getScoreCalcXy(calc,true);
  ScoreCalcXy * tempCalcB = ScoreCalcXy::getScoreCalcXy(calc,false);
  getSegmentsFromBlock(&aLinks,aStartPos,aEndPos,senseA,_seqA,tempCalcA);
  getSegmentsFromBlock(&bLinks,bStartPos,bEndPos,senseB,_seqB,tempCalcB);


  // compute the means and totals
  float aTotal = 0;
  for (map<long,float>::iterator it = aLinks.begin(); it != aLinks.end(); ++it){
    aTotal += it->second;
  }
  float aMean = aTotal / aLinks.size();
  if (aMean < 0.0){ aMean = 0.0; }
  float bTotal = 0;
  for (map<long,float>::iterator it = bLinks.begin(); it != bLinks.end(); ++it){
    bTotal += it->second;
  }
  float bMean = bTotal / bLinks.size();
  if (bMean < 0.0){ bMean = 0.0; }


  // pick the winner and add to the contig
  float penaltyConstant;
  if (aMean > bMean){ // seqA won
    if (aTotal==0){ penaltyConstant = 1; } // avoid divide-by-zero error
    else { penaltyConstant = 1.0 - bMean * aLinks.size() / aTotal; }
    fillInGaps(carrier, infoA, aStartPos, aEndPos, penaltyConstant, &aLinks, tempCalcA);
  } else if (aMean < bMean) { // seqB won
    if (bTotal==0){ penaltyConstant = 1; } // avoid divide-by-zero error
    else { penaltyConstant = 1.0 - aMean * bLinks.size() / bTotal; }
    fillInGaps(carrier, infoB, bStartPos, bEndPos, penaltyConstant, &bLinks, tempCalcB);
  } else { // there was a tie; this could be due to an edge-of-sequence issue
    // we are at the edge of A, so B is chosen
    if (bTotal==0){ penaltyConstant = 1; } // avoid divide-by-zero error
    else { penaltyConstant = 1.0 - aMean * bLinks.size() / bTotal; }
    if ( (aStartPos == -1 and aEndPos != _sizeA and 
	  isGappedSeqB(0) and 
	  gapPairedAfterSeqB(0) == -1 ) or
	 (aEndPos == _sizeA and aStartPos != -1 and
	  isGappedSeqB(_sizeBm1) and
	  gapPairedAfterSeqB(bEndPos-1) == _sizeAm1) ){
      fillInGaps(carrier, infoB, bStartPos, bEndPos, penaltyConstant, &bLinks, tempCalcB);
    }
    // we are at the edge of B, so A is chosen - or, it is just a tie so A is chosen
    else { fillInGaps(carrier, infoA, aStartPos, aEndPos, penaltyConstant, &aLinks, tempCalcA); }
  }
  delete tempCalcA;
  delete tempCalcB;
}




AlignmentGapped::CombineDataCarrier::CombineDataCarrier(){}
AlignmentGapped::CombineDataCarrier::CombineDataCarrier(long maxSize) : _maxSize(maxSize) {
  _seq = new char[_maxSize];
  _scores = new float[_maxSize];
  _links = new float[_maxSize];
  _seq[0] = '\0';
  _size = 0;
  _linkOffset = 0;
}
AlignmentGapped::CombineDataCarrier::~CombineDataCarrier(){
  delete [] _seq;
  delete [] _scores;
  delete [] _links;
}
void AlignmentGapped::CombineDataCarrier::addNuc(char nuc, float score){
  _seq[_size] = nuc;
  _scores[_size] = score;
  _linkOffset--;
  _size++;
  _seq[_size] = '\0';
}
void AlignmentGapped::CombineDataCarrier::addNuc(char nuc, float score, float link){
  _seq[_size] = nuc;
  _scores[_size] = score;
  _links[_size + _linkOffset] = link;
  _size++;
  _seq[_size] = '\0';
}
void AlignmentGapped::CombineDataCarrier::addLink(float link){
  _links[_size + _linkOffset] = link;
  _linkOffset++;
}
void AlignmentGapped::CombineDataCarrier::decrementLink(){
  _linkOffset--;
}
ScoredSeq* AlignmentGapped::CombineDataCarrier::getSeq(){
  return new ScoredSeqShallow(_seq, _scores, _links, _size);
}

inline bool AlignmentGapped::isPositionFinal(ScoredSeq * seqX, long posX, ScoredSeq * seqY, long posY){
  return isPositionFinal(seqX,posX) or isPositionFinal(seqY,posY);
}
inline bool AlignmentGapped::isPositionFinal(ScoredSeq * seq, long pos){
  return pos + 1 == seq->size();
}



#endif
