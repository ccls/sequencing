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

#ifndef DYNAMICPROGRAMMINGALIGNER_CPP
#define DYNAMICPROGRAMMINGALIGNER_CPP


#include "DynamicProgrammingAligner.h"
#include "AlignmentScoreMatrix.h"
#include "ScoredSeq.h"
#include "AlignmentNull.h"
#include "AssemblyException.h"
#include <map>
#include <iostream>
#include <cmath>
using namespace::std;

void DynamicProgrammingAligner::OK(){
  if (_fractId < 0 or _fractId > 1){
    cerr << _fractId << endl;
    throw AssemblyException::LogicError("_fractId is a nonsense value.");
  }
}

DynamicProgrammingAligner::DynamicProgrammingAligner(){}
DynamicProgrammingAligner::DynamicProgrammingAligner(float fractId,
						     long minOverlap,
						     AlignmentScoreMatrix * scoreMatrix) :
  _fractId(fractId),
  _minOverlap(minOverlap),
  _maxOverlap(0),
  _usingMaxOverlap(false),
  _minScore(0) { 
  _scoreMatrix = new AlignmentScoreMatrix(scoreMatrix);
  //OK(); 
}
DynamicProgrammingAligner::DynamicProgrammingAligner(float fractId,
						     long minOverlap,
						     long maxOverlap,
						     AlignmentScoreMatrix * scoreMatrix) :
  _fractId(fractId),
  _minOverlap(minOverlap),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap),
  _minScore(0) { 
  _scoreMatrix = new AlignmentScoreMatrix(scoreMatrix);
  //OK(); 
}
DynamicProgrammingAligner::DynamicProgrammingAligner(float fractId,
						     long minOverlap,
						     AlignmentScoreMatrix * scoreMatrix,
						     long minScore) :
  _fractId(fractId),
  _minOverlap(minOverlap),
  _usingMaxOverlap(false),
  _maxOverlap(0),
  _minScore(minScore) {
  _scoreMatrix = new AlignmentScoreMatrix(scoreMatrix);
  //OK();
}
DynamicProgrammingAligner::DynamicProgrammingAligner(float fractId,
						     long minOverlap,
						     long maxOverlap,
						     AlignmentScoreMatrix * scoreMatrix,
						     long minScore) :
  _fractId(fractId),
  _minOverlap(minOverlap),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap),
  _minScore(minScore) {
  _scoreMatrix = new AlignmentScoreMatrix(scoreMatrix);
  //OK();
}
DynamicProgrammingAligner::DynamicProgrammingAligner(DynamicProgrammingAligner* dpa){
  _fractId = dpa->_fractId;
  _minOverlap = dpa->_minOverlap;
  _usingMaxOverlap = dpa->_usingMaxOverlap;
  _maxOverlap = dpa->_maxOverlap;
  _minScore = dpa->_minScore;
  _scoreMatrix = new AlignmentScoreMatrix(dpa->_scoreMatrix);
}
DynamicProgrammingAligner::~DynamicProgrammingAligner(){
  delete _scoreMatrix;
}


float DynamicProgrammingAligner::getFractId(){ return _fractId; }
long DynamicProgrammingAligner::getMinOverlap(){ return _minOverlap; }
long DynamicProgrammingAligner::getMinScore(){ return _minScore; }
AlignmentScoreMatrix * DynamicProgrammingAligner::getScoreMatrix(){ return new AlignmentScoreMatrix(_scoreMatrix); }

DynamicProgrammingAligner * DynamicProgrammingAligner::minScoreCopy(long minScore){
  //OK();
  DynamicProgrammingAligner* newDpa;
  if (! _usingMaxOverlap){
    newDpa = new DynamicProgrammingAligner( _fractId, _minOverlap, _scoreMatrix, _minScore );
  } else {
    newDpa = new DynamicProgrammingAligner( _fractId, _minOverlap, _maxOverlap, _scoreMatrix, _minScore );
  }
  return newDpa;
}


// REQUIRES: seqA is as long or longer than seqB
DynamicProgrammingAligner::EdgeBoundaryCarrierCollection* DynamicProgrammingAligner::getOverlapLimits(ScoredSeq* seqA, ScoredSeq* seqB,
												      long minOverlap, long maxOverlap){
  if (seqB->size() > seqA->size()){
    throw AssemblyException::ImplementationError("DPA::getOverlapLimits seq sizes inappropriate.");
  }

  EdgeBoundaryCarrierCollection* ebcc;

  // this condition ensures that the full semi-global alignment will
  // be performed if it is allowed.
  if ( maxOverlap >= seqB->size() ){
    ebcc = new EdgeBoundaryCarrierCollection(1);
    EdgeBoundaryCarrier* ebc = new EdgeBoundaryCarrier(0 - seqB->size() + _minOverlap, seqA->size() - _minOverlap);
    ebcc->_ebcArray[0] = ebc;
  } else {
    ebcc = new EdgeBoundaryCarrierCollection(2);
    // the boundary conditions will differ for each of the two alignments
    // for each of the two alignments
    EdgeBoundaryCarrier* ebc;

    // this alignment has the 3p end of A overhanging
    long initialPosA = seqA->size() - maxOverlap - 1;
    ebc = new EdgeBoundaryCarrier(initialPosA + 1, seqA->size() - _minOverlap);
    ebcc->_ebcArray[0] = ebc;

    // this alignment has the 3p end of B overhanging
    long initialPosB = seqB->size() - maxOverlap - 1;
    ebc = new EdgeBoundaryCarrier(0 - seqB->size() + _minOverlap, 0 - initialPosB - 1);
    ebcc->_ebcArray[1] = ebc;
  }
  return ebcc;
}


Alignment * DynamicProgrammingAligner::align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation){
  //OK();
  // do the alignment job with a value that is clearly above the upper bound
  // for maxOverlap.
  return align(extSeqA,extSeqB,orientation,extSeqA->size()+extSeqB->size());
}


Alignment * DynamicProgrammingAligner::align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long maxOverlap){
  //OK();

  // previousLine prevents long look-up times for seqA positions and
  // scales memory usage to the length of seqB, so the longer of the two
  // input seqs should be selected as seqA.
  ScoredSeq * seqA = 0;
  ScoredSeq * seqB = 0;
  char senseA;
  char senseB;
  if ( extSeqA->size() >= extSeqB->size() ){
    seqA = extSeqA;
    seqB = extSeqB;
    senseA = '+';
    senseB = orientation;
  } else {
    seqA = extSeqB;
    seqB = extSeqA;
    senseA = orientation;
    senseB = '+';
  }

  // adjust max overlap to conform to the global value
  if (_usingMaxOverlap and maxOverlap > _maxOverlap){ maxOverlap = _maxOverlap; }

  // this will be replaced if a satisfactory square is found
  Square* bestSquare = 0;

  // get the best square from all of the possible windows
  EdgeBoundaryCarrierCollection* overlapLimits = getOverlapLimits(seqA, seqB, _minOverlap, maxOverlap);

  for (long limitN = 0; limitN < overlapLimits->_size; ++limitN){
    EdgeBoundaryCarrier* limitEbc = overlapLimits->_ebcArray[limitN];
    Square* bestCand = alignHelp(seqA,seqB,senseA,senseB,limitEbc);
    bestSquare = getNewBestRef( bestSquare, bestCand );
    if ( bestCand != bestSquare ){ delete bestCand; }
    delete limitEbc;
  }
  delete overlapLimits;

  // the alignment, which will (maybe) be mutated by alignmentFrom Square.
  Alignment* alignment;
  if (bestSquare == 0){
    alignment = new AlignmentNull(extSeqA,extSeqB,orientation);
  } else {
    AlignmentGapped* alignmentValid = new AlignmentGapped(extSeqA,extSeqB,orientation);
    alignment = alignmentValid;
    DynamicProgrammingAligner::alignmentFromSquare(bestSquare, alignmentValid, seqA, seqB);
  }
  return alignment;
}


// offset is the position in extSeqA where extSeqB starts (negative if extSeqB has a 5p overhang)
Alignment * DynamicProgrammingAligner::align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long minOffset, long maxOffset){
  //OK();

  // adjust the min/max offsets based on the fractId
  long offsets[] = {minOffset, maxOffset};
  for (int n = 0; n < 2; ++n){
    long startOverlap;
    if (offsets[n] <= 0){ startOverlap = 0; }
    else { startOverlap = offsets[n]; }
    long endOverlap;
    if (extSeqA->size() <= extSeqB->size() + offsets[n]){ endOverlap = extSeqA->size(); }
    else { endOverlap = extSeqB->size() + offsets[n]; }
    long overlap = endOverlap - startOverlap;
    long ovlVariance = long(float(overlap) * (1.0 - _fractId));
    if (n==0){ minOffset -= ovlVariance; }
    else { maxOffset += ovlVariance; }
  }


  // previousLine prevents long look-up times for seqA positions and
  // scales memory usage to the length of seqB, so the longer of the two
  // input seqs should be selected as seqA.
  ScoredSeq * seqA = 0;
  ScoredSeq * seqB = 0;
  char senseA;
  char senseB;
  if ( extSeqA->size() >= extSeqB->size() ){
    seqA = extSeqA;
    seqB = extSeqB;
    senseA = '+';
    senseB = orientation;
  } else {
    seqA = extSeqB;
    seqB = extSeqA;
    senseA = orientation;
    senseB = '+';
    long oldMinOffset = minOffset;
    minOffset = 0 - maxOffset;
    maxOffset = 0 - oldMinOffset;
  }

  // apply the min/max overlap limits
  long maxOverlap;
  if( _usingMaxOverlap ){ maxOverlap = _maxOverlap; }
  else { maxOverlap = seqA->size() + seqB->size(); }
  EdgeBoundaryCarrierCollection* overlapLimits = getOverlapLimits(seqA, seqB, _minOverlap, maxOverlap);

  // generate an initial ebc that will then constrain the min/max overlaps
  vector<EdgeBoundaryCarrier*> ebcSet;
  EdgeBoundaryCarrier* initialEbc = new EdgeBoundaryCarrier(minOffset, maxOffset);
  for (long limitN = 0; limitN < overlapLimits->_size; ++limitN){
    EdgeBoundaryCarrier* limitEbc = overlapLimits->_ebcArray[limitN];
    if ( initialEbc->overlaps(limitEbc) ){ ebcSet.push_back( initialEbc->getOverlap(limitEbc) ); }
    delete limitEbc;
  }
  delete initialEbc;
  delete overlapLimits;

  // this will be replaced if a satisfactory square is found
  Square* bestSquare = 0;

  for (vector<EdgeBoundaryCarrier*>::iterator ebcIt = ebcSet.begin(); ebcIt != ebcSet.end(); ++ebcIt){
    Square* bestCand = alignHelp(seqA,seqB,senseA,senseB,*ebcIt);
    bestSquare = getNewBestRef( bestSquare, bestCand );
    if ( bestCand != bestSquare ){ delete bestCand; }
    delete *ebcIt;
  }

  // the alignment, which will (maybe) be mutated by alignmentFrom Square.
  Alignment* alignment;
  if (bestSquare == 0){
    alignment = new AlignmentNull(extSeqA,extSeqB,orientation);
  } else {
    AlignmentGapped* alignmentValid = new AlignmentGapped(extSeqA,extSeqB,orientation);
    alignment = alignmentValid;
    DynamicProgrammingAligner::alignmentFromSquare(bestSquare, alignmentValid, seqA, seqB);
  }
  return alignment;
}




// offset is the position in extSeqA where extSeqB starts (negative if extSeqB has a 5p overhang)
Alignment * DynamicProgrammingAligner::align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, set<long>* inputOffsets){
  long* offsetArray = new long[ inputOffsets->size() + 1 ];
  long offsetIndex = 0;
  for (set<long>::iterator it = inputOffsets->begin(); it != inputOffsets->end(); ++it){
    offsetArray[offsetIndex] = *it;
    offsetIndex++;
  }
  Alignment* al = align(extSeqA, extSeqB, orientation, offsetArray, long(inputOffsets->size()));
  delete [] offsetArray;
  return al;
}



// offset is the position in extSeqA where extSeqB starts (negative if extSeqB has a 5p overhang)
Alignment * DynamicProgrammingAligner::align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long* offsetArray, long numOffsets){
  //OK();

  // this will replace the input offsets (they may need to be adjusted if the input sequences are switched
  long* offsets = new long[ numOffsets + 1 ];

  // previousLine prevents long look-up times for seqA positions and
  // scales memory usage to the length of seqB, so the longer of the two
  // input seqs should be selected as seqA.
  ScoredSeq * seqA = 0;
  ScoredSeq * seqB = 0;
  char senseA;
  char senseB;
  if ( extSeqA->size() >= extSeqB->size() ){
    seqA = extSeqA;
    seqB = extSeqB;
    senseA = '+';
    senseB = orientation;
    for (long n = 0; n < numOffsets; ++n){ offsets[n] = offsetArray[n]; }
  } else {
    seqA = extSeqB;
    seqB = extSeqA;
    senseA = orientation;
    senseB = '+';
    // switch around the offsets
    for (long n = 0; n < numOffsets; ++n){ offsets[n] = 0 - offsetArray[n]; }
  }

  // apply the min/max overlap limits
  long maxOverlap;
  if( _usingMaxOverlap ){ maxOverlap = _maxOverlap; }
  else { maxOverlap = seqA->size() + seqB->size(); }
  EdgeBoundaryCarrierCollection* overlapLimits = getOverlapLimits(seqA, seqB, _minOverlap, maxOverlap);

  // remove the offsets that do not lie in the allowed overlap limits                                                                                  
  long* validOffsets = new long[numOffsets + 1];
  long numValidOffsets = 0;
  for (long inputIndex = 0; inputIndex < numOffsets; ++inputIndex){
    bool ovlIsOk = false;
    for (long limitN = 0; limitN < overlapLimits->_size; ++limitN){
      if (offsets[inputIndex] >= overlapLimits->_ebcArray[limitN]->_minOffset and
          offsets[inputIndex] <= overlapLimits->_ebcArray[limitN]->_maxOffset){
        ovlIsOk = true;
      } 
    }
    if ( ovlIsOk ){
      validOffsets[numValidOffsets] = offsets[inputIndex];
      numValidOffsets++;
    }   
  }
  // make the initial offset boundary carriers
  EdgeBoundaryCarrierCollection* inputEbcc = new EdgeBoundaryCarrierCollection(numValidOffsets);
  for (long inputIndex = 0; inputIndex < numValidOffsets; ++inputIndex){
    long offset = validOffsets[inputIndex];
    // adjust the min/max offsets based on the fractId
    long startOverlap;
    if (offset <= 0){ startOverlap = 0; }
    else { startOverlap = offset; }
    long endOverlap;
    if (extSeqA->size() <= seqB->size() + offset){ endOverlap = seqA->size(); }
    else { endOverlap = seqB->size() + offset; }
    long overlap = endOverlap - startOverlap;
    long ovlVariance = adjustMaxMismatch(overlap - long(_fractId * float(overlap)));
    // generate an initial ebc that will then constrain the min/max overlaps
    EdgeBoundaryCarrier* newEbc = new EdgeBoundaryCarrier(offset - ovlVariance, offset + ovlVariance, offset);
    inputEbcc->_ebcArray[inputIndex] = newEbc;
  }

  // make the set non-redundant
  inputEbcc = makeEbccNr(inputEbcc);

  // remove the portions that don't jive with the min/max overlap limits
  EdgeBoundaryCarrierCollection* limitedEbcc = new EdgeBoundaryCarrierCollection(inputEbcc->_size * overlapLimits->_size);
  long limitedIndex = 0;
  for (long ebcN = 0; ebcN < inputEbcc->_size; ++ebcN){
    EdgeBoundaryCarrier* inputEbc = inputEbcc->_ebcArray[ebcN];
    for (long limitN = 0; limitN < overlapLimits->_size; ++limitN){
      EdgeBoundaryCarrier* limitEbc = overlapLimits->_ebcArray[limitN];
      if ( inputEbc->overlaps(limitEbc) ){
	limitedEbcc->_ebcArray[limitedIndex] = inputEbc->getOverlap(limitEbc);
	limitedIndex++;
      }
    }
    delete inputEbc;
  }
  for (long limitN = 0; limitN < overlapLimits->_size; ++limitN){ delete overlapLimits->_ebcArray[limitN]; }
  delete overlapLimits;

  limitedEbcc->_size = limitedIndex;
  delete inputEbcc;

  // make the constrained portions non-redundant and sorted
  limitedEbcc = makeEbccNr(limitedEbcc);
  sortEbccByOffsets(limitedEbcc, validOffsets, numValidOffsets);

  // a null alignment is the base case; score zero is the base case
  long bestScoreSoFar = 0;
  Alignment* bestAlignment = NULL;
  DynamicProgrammingAligner* localExecutor = this->minScoreCopy(bestScoreSoFar);

  for (long ebcN = 0; ebcN < limitedEbcc->_size; ++ebcN){
    // this will be replaced if a satisfactory square is found
    Square* bestSquare = NULL;
    bestSquare = localExecutor->alignHelp(seqA,seqB,senseA,senseB,limitedEbcc->_ebcArray[ebcN]);
    delete limitedEbcc->_ebcArray[ebcN];

    // the alignment, which will (maybe) be mutated by alignmentFrom Square.
    if (bestSquare != NULL){
      AlignmentGapped * alignment = new AlignmentGapped(extSeqA,extSeqB,orientation);
      DynamicProgrammingAligner::alignmentFromSquare(bestSquare, alignment, seqA, seqB);
      if ( alignment->isNull() ){ delete alignment; }
      else {
	long newScore = alignment->score(_scoreMatrix, false);
	if ( bestAlignment == NULL or newScore > bestScoreSoFar){
	  bestScoreSoFar = newScore;
	  DynamicProgrammingAligner* newLocalExecutor = this->minScoreCopy(bestScoreSoFar);
	  delete localExecutor;
	  localExecutor = newLocalExecutor;
	  delete bestAlignment;
	  bestAlignment = alignment;
	} else { delete alignment; }
      }
    }
  }
  delete localExecutor;
  delete limitedEbcc;
  delete [] offsets;
  delete [] validOffsets;
  if (bestAlignment==NULL){
    bestAlignment = new AlignmentNull(extSeqA,extSeqB,orientation);
  }
  return bestAlignment;
}

DynamicProgrammingAligner::EdgeBoundaryCarrierCollection* DynamicProgrammingAligner::makeEbccNr(EdgeBoundaryCarrierCollection* oldEbcc){

  // initial values
  long nrEbcCount = 0;
  EdgeBoundaryCarrier** nrEbcArray = new EdgeBoundaryCarrier*[1]; // 1 instead of 0 so it is not null

  // iterate through elements of the old collection
  for (long oldEbccIndex = 0; oldEbccIndex < oldEbcc->_size; ++oldEbccIndex){

    // generate an initial ebc that will then constrain the min/max overlaps
    EdgeBoundaryCarrier* newEbc = oldEbcc->_ebcArray[oldEbccIndex];

    // keeps track of which ones were redundant
    long newNrEbcCount = 1; // starts with the count of newEbc
    bool* notUsed = new bool[nrEbcCount+1]; // +1 so that it won't be null
    for (long n = 0; n < nrEbcCount; ++n){
      if ( newEbc->isContinuous(nrEbcArray[n]) ){
	EdgeBoundaryCarrier* replacement = newEbc->combine( nrEbcArray[n] );
	delete nrEbcArray[n];
	delete newEbc;
	newEbc = replacement;
	notUsed[n] = false;
      } else {
	notUsed[n] = true;
	newNrEbcCount++;
      }
    }
    EdgeBoundaryCarrier** newNrEbcArray = new EdgeBoundaryCarrier*[newNrEbcCount];
    newNrEbcArray[0] = newEbc;
    long currentIndex = 1;
    for (long n = 0; n < nrEbcCount; ++n){
      if ( notUsed[n] ){
	newNrEbcArray[currentIndex] = nrEbcArray[n];
	currentIndex++;
      }
    }
    delete [] nrEbcArray;
    nrEbcArray = newNrEbcArray;
    nrEbcCount = newNrEbcCount;
    delete [] notUsed;
  }
  delete oldEbcc;
  EdgeBoundaryCarrierCollection* newEbcc = new EdgeBoundaryCarrierCollection();
  newEbcc->_size = nrEbcCount;
  newEbcc->_ebcArray = nrEbcArray;
  return newEbcc;
}


// MODIFIES ebcc
void DynamicProgrammingAligner::sortEbccByOffsets(EdgeBoundaryCarrierCollection* ebcc, long* sortedOffsets, long numOffsets){

  // the new sorted array
  EdgeBoundaryCarrier** sortedEbcArray = new EdgeBoundaryCarrier*[ ebcc->_size + 1 ];
  long nSorted = 0;

  bool* ebcAvailable = new bool[ ebcc->_size + 1 ];
  for (long n = 0; n < ebcc->_size; ++n){ ebcAvailable[n] = true; }

  long nOffset = 0;
  while (nOffset < numOffsets and nSorted < ebcc->_size){
    long offset = sortedOffsets[nOffset];
    nOffset++;

    for (long ebcIndex = 0; ebcIndex < ebcc->_size; ++ebcIndex){
      if (ebcAvailable[ebcIndex] and
          offset >= ebcc->_ebcArray[ebcIndex]->_minOffset and
	  offset <= ebcc->_ebcArray[ebcIndex]->_maxOffset){
        sortedEbcArray[nSorted] = ebcc->_ebcArray[ebcIndex];
	ebcAvailable[ebcIndex] = false;
	nSorted++;
      }
    }
  }
  if (nSorted < ebcc->_size){ throw AssemblyException::LogicError("DPA::sort not all ebc's have been transferred."); }

  // replace unsroted order with sorted order
  for (long n = 0; n < ebcc->_size; ++n){ ebcc->_ebcArray[n] = sortedEbcArray[n]; }
  delete [] sortedEbcArray;
  delete [] ebcAvailable;
}


DynamicProgrammingAligner::Square * DynamicProgrammingAligner::alignHelp(ScoredSeq * seqA, ScoredSeq * seqB, char senseA, char senseB, EdgeBoundaryCarrier* ebc){

  // these are for method optimizations, since these values are frequently obtained
  long seqAsizeM1 = seqA->size() - 1;
  long seqBsizeM1 = seqB->size() - 1;

  // figure out the maximum number of mismatches that will be allowed; this is a meta-statistic
  // that can be varied from this position to affect the alignment efficiency
  long maxOverlapLength = getMaxOverlapLength(seqA->size(), seqB->size(), ebc);
  // subtract the num matches so that mismatch value will end up being rounded up
  //long maxNumMismatch = adjustMaxMismatch(maxOverlapLength - long(float(maxOverlapLength) * _fractId));
  long maxNumMismatch = maxOverlapLength - long(float(maxOverlapLength) * _fractId);

  // calculate initial positions based on offsets
  long initialPosA = -1;
  long initialPosB = -1;
  if (ebc->_minOffset > 0){ initialPosA = ebc->_minOffset - 1; }
  else if (ebc->_maxOffset < 0){ initialPosB = 0 - ebc->_maxOffset - 1; }

  // figure out the max value of posB and the max length of the row
  long finalPosB = seqBsizeM1;
  if ( ebc->_maxOffset + seqBsizeM1 > seqAsizeM1 ){ finalPosB = ebc->_maxOffset + seqBsizeM1 - seqAsizeM1; }
  long addedWidth = ebc->_maxOffset - ebc->_minOffset + 1 - finalPosB + initialPosB;
  if (addedWidth < 0){ addedWidth = 0; }

  // temporary holders that will be used in the loops below
  // remember that these "lines" are upward-slanting diagonals across the grid
  DpaLine* currentLine = new DpaLine(initialPosB, finalPosB, addedWidth);
  DpaLine* oneNextLine = new DpaLine(initialPosB, finalPosB, addedWidth);
  DpaLine* twoNextLine = new DpaLine(initialPosB, finalPosB, addedWidth);

  Square * bestSquare = NULL; // the best square in the grid so far

  // establish the initial square and the current line
  Square * initialSquare = new Square(initialPosA,initialPosB,0);
  currentLine->addSquare(initialSquare, initialPosB);
  initialSquare->_isInGrid = true;

  // do the rest of the lines of the matrix
  long seqSizeBoundary = seqAsizeM1 + seqBsizeM1 + 2;
  for (long nA = initialPosA;
       // first condition defines the extent of the grid; second allows early termination
       // if there are no possible solutions.
       (nA < seqSizeBoundary) and (currentLine->_numActive + oneNextLine->_numActive + twoNextLine->_numActive > 0);
       ++nA){ // the -1 line was already made

    // iterate through the priorLine to generate entries on the currentLine
    for (DpaLine::Iterator it = currentLine->begin(); it!=currentLine->end(); ++it){

      Square* sourceSquare = *it;
      long posA = sourceSquare->_posA;
      long posB = sourceSquare->_posB;
      // so that they don't have to be repeatedly added below
      long posAp1 = posA+1;
      long posBp1 = posB+1;

      // gap in B (same posB as sourceSquare)
      if ( posA < seqAsizeM1 ) {
	Square* candNewSquare = new Square(posAp1,posB,posAp1,seqA,senseA,sourceSquare,_scoreMatrix);
	bestSquare = evaluateNewPosition(candNewSquare, oneNextLine, seqAsizeM1, seqBsizeM1, bestSquare, ebc, maxNumMismatch);
      }
      // gap in A (same posA as sourceSquare)
      if ( posB < seqBsizeM1 ) {
	Square* candNewSquare = new Square(posA,posBp1,posBp1,seqB,senseB,sourceSquare,_scoreMatrix);
	bestSquare = evaluateNewPosition(candNewSquare, oneNextLine, seqAsizeM1, seqBsizeM1, bestSquare, ebc, maxNumMismatch);
      }
      // match (posAp1 and posBp1 versus sourceSquare)
      if ( posA < seqAsizeM1 and posB < seqBsizeM1 ) {
	Square* candNewSquare = new Square(posAp1,posBp1,seqA,seqB,senseA,senseB,sourceSquare,_scoreMatrix);
	bestSquare = evaluateNewPosition(candNewSquare, twoNextLine, seqAsizeM1, seqBsizeM1, bestSquare, ebc, maxNumMismatch);
      }
    }

    // delete currentLine contents
    for (DpaLine::Iterator it = currentLine->begin(); it!=currentLine->end(); ++it){
      Square* candForDestruction = *it;
      candForDestruction->_isInGrid = false;
      if ( candForDestruction->okToDelete() ){ delete candForDestruction; }
    }

    // rotate the lines
    DpaLine* currentHolder = currentLine;
    currentLine = oneNextLine;
    oneNextLine = twoNextLine;
    long currentMin;
    long currentMax;
    if (currentLine->_currentSet and oneNextLine->_currentSet){
      if (currentLine->_currentMin < oneNextLine->_currentMin){ currentMin = currentLine->_currentMin; }
      else { currentMin = oneNextLine->_currentMin; }
      if (currentLine->_currentMax > oneNextLine->_currentMax){ currentMax = currentLine->_currentMax; }
      else { currentMax = oneNextLine->_currentMax; }
    } else if (currentLine->_currentSet){
      currentMin = currentLine->_currentMin;
      currentMax = currentLine->_currentMax;
    } else if (oneNextLine->_currentSet){
      currentMin = oneNextLine->_currentMin;
      currentMax = oneNextLine->_currentMax;
    } else {
      currentMin = 0;
      currentMax = 0;
    }
    twoNextLine = currentHolder;
    currentHolder->recycle(currentMin, currentMax);
  }

  // delete the lines
  delete currentLine;
  delete oneNextLine;
  delete twoNextLine;

  if ( bestSquare != 0 ){
    bestSquare->_isInGrid = false;
    bestSquare->_isBest = false;
  }
  return bestSquare;
}


// use the current best square to generate linkages within the alignment.  lock the alignment.
void DynamicProgrammingAligner::alignmentFromSquare(Square* bestSquare, AlignmentGapped* alignment, ScoredSeq* seqA, ScoredSeq* seqB){

  if ( bestSquare != 0 ){ // if there is no bestSquare, the alignment will remain empty (null)
    Square * currentSquare = bestSquare;
    while ( currentSquare != 0 and currentSquare->_origin != 0 ){  //and currentSquare->posA() > -1 and currentSquare->posB() > -1 ) {
      if ( currentSquare->isMatch() ){
	alignment->addLinkage( seqA, currentSquare->_posA, currentSquare->_posB );
      } else {
	// this is confusing; i should change the alignment construction interface
	if ( currentSquare->_extraNucSeq == seqA ){
	  alignment->addGap( seqA, currentSquare->_posA, currentSquare->_posB );
	} else if ( currentSquare->_extraNucSeq == seqB ){
	  alignment->addGap( seqB, currentSquare->_posB, currentSquare->_posA );
	} else {
	  throw AssemblyException::LogicError("in DynamicProgrammingAligner::alignmentFromSquare; bad seq with ins was returned.");
	}
      }
      if ( currentSquare->_origin != 0 ){
	currentSquare = currentSquare->_origin;
      } else { currentSquare = 0; }
    }
    bestSquare->_isBest = false;
    delete bestSquare;
  }

  alignment->lock();
}



// returns the current bestSquare based on its own evaluation
// MUST DELETE either candNewSquare OR the old square it was supposed to replace
// MUST DELETE old bestSquare IF old bestSquare != new bestSquare AND old bestSquare is okToDelete
DynamicProgrammingAligner::Square* DynamicProgrammingAligner::evaluateNewPosition(Square* candNewSquare, DpaLine* nextLine,
										  long seqAsizeM1, long seqBsizeM1,
										  Square* bestSquare, EdgeBoundaryCarrier* ebc,
										  long maxNumMismatch){
  Square* priorNewSquare = nextLine->getSquare( candNewSquare->_posB );

  // this was optimized to ask the smallest number of questions non-repetitively, at the expense
  // of the code being more verbose

  if (priorNewSquare == NULL){
    // here, there was not a square to begin with, so the new square is added if it is OK
    if ( isSquareLegit(candNewSquare,seqAsizeM1,seqBsizeM1,ebc,maxNumMismatch) ){
      // add the new square to the grid
      nextLine->addSquare(candNewSquare, candNewSquare->_posB);
      candNewSquare->_isInGrid = true;
      // check to see if the new square is the best square - it might be!
      if ( candNewSquare->_posA == seqAsizeM1 or candNewSquare->_posB == seqBsizeM1 ){
	bestSquare = getNewBestRef(bestSquare, candNewSquare);
      }
    } else { delete candNewSquare; }

    // otherwise, the new square must be compared to the old square (and be legit)
  } else  if (priorNewSquare->_score < candNewSquare->_score and isSquareLegit(candNewSquare,seqAsizeM1,seqBsizeM1,ebc,maxNumMismatch) ){
    // it is better, so I add the new square to the grid
    nextLine->addSquare(candNewSquare, candNewSquare->_posB);
    candNewSquare->_isInGrid = true;
    // check to see if the new square is the best square - it might be! (if it is replacing
    // a square that was the best square, then it is by definition because it has a higher score)
    if (priorNewSquare->_isBest){
      candNewSquare->_isBest = true;
      bestSquare = candNewSquare;
    } else if ( candNewSquare->_posA == seqAsizeM1 or candNewSquare->_posB == seqBsizeM1 ){
      bestSquare = getNewBestRef(bestSquare, candNewSquare);
    }
    delete priorNewSquare;

    // in this case, the new square will not be used, so it is deleted
  } else { delete candNewSquare; }

  return bestSquare;
}



DynamicProgrammingAligner::Square* DynamicProgrammingAligner::getNewBestRef(Square* bestSquare, Square* bestCand){
  if ( bestCand != NULL ){
    // more repetitive code, but fewer operations
    if ( bestSquare == NULL ){
      bestSquare = bestCand;
      bestSquare->_isBest = true;
    } else if (bestCand->_score > bestSquare->_score) {
      // the self-equality check that used to be here is superfluous because the
      // this is the only place where it could cause a problem and if the squares were
      // the same then their scores would too so this deletion would never be reached
      bestSquare->_isBest = false;
      if ( bestSquare->okToDelete() ){ delete bestSquare; }
      bestSquare = bestCand;
      bestSquare->_isBest = true;
    }
  }
  return bestSquare;
}



// figures out if, in the best case scenario, the alignment could succeed
// the percent ID is measured as the percentage of positions in the alignment
// that are matches (so the denominator includes the total number of gaps)
bool DynamicProgrammingAligner::isSquareLegit(Square* newSquare, long seqAsizeM1, long seqBsizeM1, EdgeBoundaryCarrier* ebc,
					      long maxNumMismatch){

  // first, check for the total num misses (also easy and fast and likely to fail)
  if (maxNumMismatch < newSquare->_sumGapAndMis ){ return false; }

  // second, check the offset boundaries, since this is easy and fast
  long offset = newSquare->_posA - newSquare->_posB;
  if (offset < ebc->_minOffset or offset > ebc->_maxOffset){ return false; }

  // i will add 1 to whichever one makes it
  long ovlSoFar; // given by the position in the sequence
  if ( newSquare->_posA < newSquare->_posB ){ ovlSoFar = newSquare->_posA + 1; }
  else { ovlSoFar = newSquare->_posB + 1; }

  // the remaining overlap could be defined by the short end of either sequence
  long ovlRemain;
  long ovlRemainA = seqAsizeM1 - newSquare->_posA;
  long ovlRemainB = seqBsizeM1 - newSquare->_posB;
  if ( ovlRemainA < ovlRemainB ){ ovlRemain = ovlRemainA; }
  else { ovlRemain = ovlRemainB; }

  // see if the projected overlap is sufficient
  if ( ovlSoFar + ovlRemain - newSquare->_sumSkip < _minOverlap ){ return false; }

  // this is split up so that calculations for stepsToOffset do not need to be performed
  // once an OK offset has been visited (it is still more efficient to calculate using
  // the value stepsToOffset==0 the first time an OK offset is visited)
  if (newSquare->_visitedTargetOffset){

    // see if the projected best-case fractId is sufficient
    long denominator = ovlRemain + newSquare->_denominator;
    if ( denominator <= 0 ){ return false; }
    else {
      if (float(newSquare->_sumMatch + ovlRemain) < _fractId * float(denominator)){ return false; }
    }

    // see if the projected best-case score is sufficient
    if ( newSquare->_score + (ovlRemain * _scoreMatrix->_match) < _minScore){ return false; }

  } else {
    // introduce a penalty for the number of gaps that would be
    // required for the alignment to visit a seed offset value
    long stepsToOffset = ebc->_distToOkOffset[ offset - ebc->_minOffset ];
    if (stepsToOffset == 0){ newSquare->_visitedTargetOffset = true; }

    // see if the projected best-case fractId is sufficient
    long denominator = ovlRemain + newSquare->_denominator + stepsToOffset;
    if ( denominator <= 0 ){ return false; }
    else {
      if (float(newSquare->_sumMatch + ovlRemain) < _fractId * float(denominator)){ return false; }
    }

    // see if the projected best-case score is sufficient
    if ( newSquare->_score + (ovlRemain * _scoreMatrix->_match) + (stepsToOffset * _scoreMatrix->_extendGap) < _minScore){ return false; }
  }

  // all tests have passed, so...
  return true;
}



// helper to simplify code dealing with overlap lengths
// this is more verbose than it needs to be to make it more efficient
long DynamicProgrammingAligner::getOverlapLength(long seqLenA, long seqLenB, long offset){
  long aEnd = offset + seqLenB;
  if (offset > 0){
    if (aEnd > seqLenA){ return seqLenA - offset; }
    else { return aEnd - offset; }
  } else {
    if (aEnd > seqLenA){ return seqLenA; }
    else { return aEnd; }
  }
}


// helper for defining the longest possible overlap for the window; useful
// for defining the max number of allowed mismatches for the corridor defined
// by the provided EBC
long DynamicProgrammingAligner::getMaxOverlapLength(long seqLenA, long seqLenB, EdgeBoundaryCarrier* ebc){
  if (ebc->_minOffset == ebc->_maxOffset){ return getOverlapLength( seqLenA, seqLenB, ebc->_minOffset ); }
  else {
    // start in the middle; this is the best position to start looking for
    // a maximum overlap value
    long currentOffset = (ebc->_minOffset + ebc->_maxOffset) / 2;
    long nextOffset = currentOffset + 1;
    long currentOvl = getOverlapLength(seqLenA, seqLenB, currentOffset);
    long nextOvl = getOverlapLength(seqLenA, seqLenB, nextOffset);
    // advance the offset as far a possible until it either reaches the edge of the corridor
    // or reaches a maximum.  because the overlaps are diagonal segments through the rectangle
    // of the alignment grid, the value has reached a maximum if the adjacent value is the same.
    if (currentOvl < nextOvl){
      while (nextOffset < ebc->_maxOffset and currentOvl < nextOvl){
	currentOffset = nextOffset;
	currentOvl = nextOvl;
	nextOffset++;
	nextOvl = getOverlapLength(seqLenA, seqLenB, nextOffset);
      }
    } else if (nextOvl < currentOvl){
      // switch the current/next identities
      long holderOffset = currentOffset;
      long holderOvl = currentOvl;
      currentOffset = nextOffset;
      currentOvl = nextOvl;
      nextOffset = holderOffset;
      nextOvl = holderOvl;
      // the same as above, but decrementing instead of incrementing
      while (nextOffset > ebc->_minOffset and currentOvl < nextOvl){
	currentOffset = nextOffset;
	currentOvl = nextOvl;
	nextOffset--;
	nextOvl = getOverlapLength(seqLenA, seqLenB, nextOffset);
      }
    }
    // the remaining possibility was that the two overlaps were equal, in which
    // case the maximum was already found and can be returned.
    return nextOvl;
  }
}


// this is used at multiple points (corridor width, maxMismatch)
long DynamicProgrammingAligner::adjustMaxMismatch(long maxMismatch){
  float rootedMismatch = sqrt(float(maxMismatch));
  float roundedDown = float(long(rootedMismatch));
  if (rootedMismatch > roundedDown){ return long(rootedMismatch) + long(1); }
  else { return long(rootedMismatch); }
}


//DynamicProgrammingAligner::DpaLine::DpaLine(){}
DynamicProgrammingAligner::DpaLine::DpaLine(long minIndex, long maxIndex, long extraCells) : 
  _numActive(0),
  _currentMin(0),
  _currentMax(0),
  _currentSet(false),
  _minIndex(minIndex),
  _maxIndex(maxIndex){
  // rowSize must be at least 1 for the arrays to be correctly made
  _rowSize = _maxIndex - _minIndex + 3 + extraCells;
  if (_rowSize < 1){ _rowSize = 1; }
  // zero index is grid position -1
  _squareArray = new Square*[ _rowSize ];
  // empty grid positions are null
  for (long n = 0; n < _rowSize; ++n){ _squareArray[n] = NULL; }
  // stores the indexes of active cells (cells cannot be de-activated
  _activeCells = new long[ _rowSize ];
}
DynamicProgrammingAligner::DpaLine::~DpaLine(){
  delete [] _squareArray;
  delete [] _activeCells;
}
void DynamicProgrammingAligner::DpaLine::recycle(long minIndex, long maxIndex) {
  long newRowSize = maxIndex - minIndex + 3;
  if (newRowSize > _rowSize){
    delete [] _squareArray;
    delete [] _activeCells;
    _squareArray = new Square*[ newRowSize ];
    _activeCells = new long[ newRowSize ];
    _rowSize = newRowSize;
    for (long n = 0; n < _rowSize; ++n){ _squareArray[n] = NULL; }
  } else {
    for (long n = 0; n < _numActive; ++n){ _squareArray[ _activeCells[n] ] = NULL; }
  }
  _numActive = 0;
  _currentMin = 0;
  _currentMax = 0;
  _currentSet = false;
  _minIndex = minIndex;
  _maxIndex = maxIndex; 
}


DynamicProgrammingAligner::DpaLine::Iterator DynamicProgrammingAligner::DpaLine::begin(){ return Iterator(this,0); }
DynamicProgrammingAligner::DpaLine::Iterator DynamicProgrammingAligner::DpaLine::end(){ return Iterator(this,_numActive); }

bool DynamicProgrammingAligner::DpaLine::hasSquare(long posB){ return _squareArray[posB - _minIndex] != NULL; }
void DynamicProgrammingAligner::DpaLine::addSquare(Square* square, long posB){
  long newIndex = posB - _minIndex;
  if (_squareArray[newIndex] == NULL){
    _activeCells[_numActive] = newIndex;
    _numActive++;
    if (_currentSet){
      // else because it can't be both
      if (posB < _currentMin){ _currentMin = posB; }
      else if (posB > _currentMax){ _currentMax = posB; }
    } else {
      _currentMin = posB;
      _currentMax = posB;
      _currentSet = true;
    }
  }
  _squareArray[newIndex] = square;
}
DynamicProgrammingAligner::Square* DynamicProgrammingAligner::DpaLine::getSquare(long posB){ return _squareArray[posB - _minIndex]; }
long DynamicProgrammingAligner::DpaLine::size(){ return _numActive; }

//DynamicProgrammingAligner::DpaLine::Iterator::Iterator(){}
DynamicProgrammingAligner::DpaLine::Iterator::Iterator(DpaLine* source, long itemNum) : _itemNum(itemNum), _source(source){}
DynamicProgrammingAligner::DpaLine::Iterator::Iterator(const Iterator& it) : _source(it._source), _itemNum(it._itemNum){}
DynamicProgrammingAligner::DpaLine::Iterator DynamicProgrammingAligner::DpaLine::Iterator::operator++(){ _itemNum++; return *this; } 
bool DynamicProgrammingAligner::DpaLine::Iterator::operator==(const Iterator& it){ return _itemNum==it._itemNum and _source==it._source; }
bool DynamicProgrammingAligner::DpaLine::Iterator::operator!=(const Iterator& it){ return _itemNum!=it._itemNum or _source!=it._source; }
DynamicProgrammingAligner::Square* DynamicProgrammingAligner::DpaLine::Iterator::operator*() { return _source->_squareArray[ _source->_activeCells[ _itemNum ] ]; }




DynamicProgrammingAligner::Square::~Square(){
  if (_deleteDeep){
    // both could be -1, and the array has to be at least 1 long
    Square** sqToDelete = new Square*[_posA + _posB + 3];
    long delIndex = 0;
    Square* currentSq = _origin;
    while (currentSq != NULL){
      currentSq->_hasChildSum--;
      if ( currentSq->okToDelete() ){
	currentSq->_deleteDeep = false;
	sqToDelete[delIndex] = currentSq;
	currentSq = currentSq->_origin;
	delIndex++;
      } else { currentSq = NULL; }
    }
    for (long n = 0; n < delIndex; ++n){ delete sqToDelete[n]; }
    delete [] sqToDelete;
  }
}


DynamicProgrammingAligner::Square::Square(long posA, long posB, long score) : 
  _deleteDeep(true),
  _isBest(false), 
  _isInGrid(false),
  _posA(posA),
  _posB(posB),
  _score(score), 
  _origin(0),
  _isMatch(false),
  _sumMatch(0),
  _sumGapAndMis(0),
  _sumSkip(0),
  _denominator(0),
  _extraNucSeq(0),
  _hasChildSum(0),
  _visitedTargetOffset(false)
 {}

// this constructor is for gapped squares only.  extraNucSeq is the sequence
// from which the nucleotide is being included (extraNucPos should equal either
// posA or posB).
DynamicProgrammingAligner::Square::Square(long posA, long posB, long extraNucPos, ScoredSeq* extraNucSeq, char extraNucSense, Square * origin, AlignmentScoreMatrix * scoreMatrix) : 
  _deleteDeep(true),
  _isBest(false), 
  _isInGrid(false),
  _posA(posA),
  _posB(posB),
  _isMatch(false),
  _origin(origin),
  _extraNucSeq(extraNucSeq),
  _hasChildSum(0),
  // these values are derived from the previous square and will be modified below
  _sumMatch(origin->_sumMatch),
  _sumGapAndMis(origin->_sumGapAndMis + 1),
  _sumSkip(origin->_sumSkip),
  _denominator(origin->_denominator + 1),
  _visitedTargetOffset(origin->_visitedTargetOffset)
 {
   origin->_hasChildSum++;

   if ( origin->_extraNucSeq == extraNucSeq ){ // determines if this is a gap extension
     _score = origin->_score + scoreMatrix->_extendGap;
   } else { _score = origin->_score + scoreMatrix->_newGap; }
   if ( extraNucSeq->scoreAtPosition(extraNucPos,extraNucSense)==0 ){
     _sumSkip++;
     _denominator--;
   }

   // provides the "to-edge-ness" on the front end
   if ( posA == -1 or posB == -1 ){
     _score = 0;
     _sumMatch = 0;
     _sumGapAndMis = 0;
     _sumSkip = 0;
     _denominator = 0;
   }
}
// this constructor is for matched (or mismatched) positions ONLY
DynamicProgrammingAligner::Square::Square(long posA, long posB, ScoredSeq* seqA, ScoredSeq* seqB, char senseA, char senseB, Square * origin, AlignmentScoreMatrix * scoreMatrix) : 
  _deleteDeep(true),
  _isBest(false), 
  _isInGrid(false),
  _posA(posA),
  _posB(posB),
  _isMatch(true),
  _origin(origin),
  _extraNucSeq(0),
  _hasChildSum(0),
  // these values are derived from the previous square and will be modified below
  _sumMatch(origin->_sumMatch),
  _sumGapAndMis(origin->_sumGapAndMis),
  _sumSkip(origin->_sumSkip),
  _denominator(origin->_denominator + 1),
  _visitedTargetOffset(origin->_visitedTargetOffset)
{
  origin->_hasChildSum++;
  if ( scoreMatrix->isMatch(seqA->nucAtPosition(posA,senseA), seqB->nucAtPosition(posB,senseB)) ){
    _score = origin->_score + scoreMatrix->_match;
    _sumMatch++;
  } else {
    _score = origin->_score + scoreMatrix->_mismatch;
    _sumGapAndMis++;
    if ( seqA->scoreAtPosition(posA,senseA)==0 or seqB->scoreAtPosition(posB,senseB)==0 ){
      _sumSkip++;
      _denominator--;
    }
  }
  // I do not check for the -1 pos match and re-set scores to zero because, since
  // this is a match, neither seqA nor seqB can be at pos==-1.
}
bool DynamicProgrammingAligner::Square::isMatch(){ return _isMatch; }
bool DynamicProgrammingAligner::Square::isGap(){ return (! _isMatch); }
bool DynamicProgrammingAligner::Square::okToDelete(){
  return ( _hasChildSum == 0 and (! _isInGrid) and (! _isBest) );
}



// the simpler case - all offsets are OK
DynamicProgrammingAligner::EdgeBoundaryCarrier::EdgeBoundaryCarrier(long minOffset, long maxOffset) :
  _minOffset(minOffset),
  _maxOffset(maxOffset)
{
  long arraySize = _maxOffset - _minOffset + 1;
  if (arraySize < 1){ arraySize = 1; }
  _distToOkOffset = new long[ arraySize ];
  for (long n = 0; n < arraySize; ++n){ _distToOkOffset[n] = 0; }
}
// the more complex case, with offset targets
DynamicProgrammingAligner::EdgeBoundaryCarrier::EdgeBoundaryCarrier(long minOffset, long maxOffset, long targetOffset) :
  _minOffset(minOffset),
  _maxOffset(maxOffset)
{
  long arraySize = _maxOffset - _minOffset + 1;
  if (arraySize < 1){ arraySize = 1; }
  _distToOkOffset = new long[ arraySize ];
  long offsetTarget = targetOffset - _minOffset;
  long borderMark = offsetTarget;
  if (borderMark > arraySize){ borderMark = arraySize; }
  for (long n = 0; n < borderMark; ++n){ _distToOkOffset[n] = offsetTarget - n; }
  for (long n = borderMark; n < arraySize; ++n){ _distToOkOffset[n] = n - offsetTarget; }
}
DynamicProgrammingAligner::EdgeBoundaryCarrier::~EdgeBoundaryCarrier(){
  delete [] _distToOkOffset;
}



bool DynamicProgrammingAligner::EdgeBoundaryCarrier::isContinuous(EdgeBoundaryCarrier* otherEbc){
  if (_minOffset > otherEbc->_maxOffset + 1){ return false; }
  if (_maxOffset + 1 < otherEbc->_minOffset){ return false; }
  return true;
}

DynamicProgrammingAligner::EdgeBoundaryCarrier* DynamicProgrammingAligner::EdgeBoundaryCarrier::combine(EdgeBoundaryCarrier* otherEbc){
  long newMinOffset;
  long newMaxOffset;
  if (_minOffset <= otherEbc->_minOffset){ newMinOffset = _minOffset; }
  else { newMinOffset = otherEbc->_minOffset; }
  if (_maxOffset >= otherEbc->_maxOffset){ newMaxOffset = _maxOffset; }
  else { newMaxOffset = otherEbc->_maxOffset; }
  EdgeBoundaryCarrier* newEbc = new EdgeBoundaryCarrier(newMinOffset,newMaxOffset);
  // now fill in the array
  for (long os = newMinOffset; os <= newMaxOffset; ++os){
    long index = os - newMinOffset;
    long thisIndex = os - _minOffset;
    long otherIndex = os - otherEbc->_minOffset;
    if (os < _minOffset or os > _maxOffset){
      newEbc->_distToOkOffset[index] = otherEbc->_distToOkOffset[ otherIndex ];
    } else if (os < otherEbc->_minOffset or os > otherEbc->_maxOffset){
      newEbc->_distToOkOffset[index] = _distToOkOffset[ thisIndex ];
    } else {
      if (_distToOkOffset[thisIndex] < otherEbc->_distToOkOffset[ otherIndex ]){
        newEbc->_distToOkOffset[index] = _distToOkOffset[ thisIndex ];
      } else {
        newEbc->_distToOkOffset[index] = otherEbc->_distToOkOffset[ otherIndex ];
      }
    }
  }
  return newEbc;
}

bool DynamicProgrammingAligner::EdgeBoundaryCarrier::overlaps(EdgeBoundaryCarrier* otherEbc){
  if (_minOffset > otherEbc->_maxOffset){ return false; }
  if (_maxOffset < otherEbc->_minOffset){ return false; }
  return true;
}
// REQUIRES: EBC's overlap
DynamicProgrammingAligner::EdgeBoundaryCarrier* DynamicProgrammingAligner::EdgeBoundaryCarrier::getOverlap(EdgeBoundaryCarrier* otherEbc){
  if (! overlaps(otherEbc) ){ throw AssemblyException::CallingError("can't call DPA::EBC::getOverlap for non-overlapping EBC's"); }
  long newMinOffset;
  long newMaxOffset;
  if (_minOffset >= otherEbc->_minOffset){ newMinOffset = _minOffset; }
  else { newMinOffset = otherEbc->_minOffset; }
  if (_maxOffset <= otherEbc->_maxOffset){ newMaxOffset = _maxOffset; }
  else { newMaxOffset = otherEbc->_maxOffset; }
  EdgeBoundaryCarrier* newEbc = new EdgeBoundaryCarrier(newMinOffset,newMaxOffset);
  // now fill in the array
  for (long os = newMinOffset; os <= newMaxOffset; ++os){
    long index = os - newMinOffset;
    long thisIndex = os - _minOffset;
    long otherIndex = os - otherEbc->_minOffset;
    if (_distToOkOffset[thisIndex] < otherEbc->_distToOkOffset[ otherIndex ]){
      newEbc->_distToOkOffset[index] = _distToOkOffset[ thisIndex ];
    } else {
      newEbc->_distToOkOffset[index] = otherEbc->_distToOkOffset[ otherIndex ];
    }
  }
  return newEbc;
}



DynamicProgrammingAligner::EdgeBoundaryCarrierCollection::EdgeBoundaryCarrierCollection(){}
DynamicProgrammingAligner::EdgeBoundaryCarrierCollection::EdgeBoundaryCarrierCollection(long size) : _size(size){
  if (_size > 0){ _ebcArray = new EdgeBoundaryCarrier*[_size]; }
  else { _ebcArray = new EdgeBoundaryCarrier*[1]; }
}
DynamicProgrammingAligner::EdgeBoundaryCarrierCollection::~EdgeBoundaryCarrierCollection(){ delete [] _ebcArray; }


#endif




/*
THE ORIGINAL PYTHON PROTOTYPE, WRITTEN BY ME BACK IN THE BARTEL LAB DAYS

class _Square:
   def __init__(self,source,score,charString):
      self._score = score
      self._source = source
      self._charString = charString
   def score(self): return self._score
   def source(self): return self._source
   def charString(self): return self._charString
   def __hash__(self): return int((score+source.__hash__()+charString.__hash__())/3) 
   def __eq__(self,other):
      return other!=None and self.score()==other.score() and self.source()==other.source() and\
             self.charString()==other.charString()
   def __ne__(self,other): return not(self.__eq__(other))


def localAlign(seqA,seqB,scoreMatrix):
  for k in scoreMatrix.keys():
    if not(_masterScores.has_key(k)):
    raise ValueError(k+' is not an acceptable score type')
    for k in filter(lambda k: _masterScores[k], _masterScores.keys()):
      if not(scoreMatrix.has_key(k)):
      raise ValueError('score matrix is missin '+k)
   
      lastRow = [_Square(None,0,'')]
    bestScoringSquare = lastRow[-1]

      for n in range(len(seqB)):
	best = _Square(lastRow[-1],
		       max([lastRow[-1].score()+scoreMatrix['egap'],0]),
		       '-'+seqB[n])
	if best.score() > bestScoringSquare.score(): bestScoringSquare = best
	lastRow.append(best)

	extendGapDict = {True:scoreMatrix['egap'], False:scoreMatrix['ngap']}
	extendGapScore = lambda extend,endOfSeq: extendGapDict[extend or endOfSeq]
  matchScore = {True:scoreMatrix['match'], False:scoreMatrix['mis']}
  for nA in range(len(seqA)):
    newRow = [_Square(lastRow[0],
		      max([lastRow[0].score() + scoreMatrix['egap'],0]),
		      seqA[nA]+'-')]
    if newRow[0].score() > bestScoringSquare.score(): bestScoringSquare = newRow[0]
    for nB in range(len(seqB)):
      aGap = _Square(newRow[-1],
		     newRow[-1].score() + extendGapScore(newRow[-1].charString()[0]=='-',
							 nA+1==len(seqA)),
		     '-'+seqB[nB])
      bGap = _Square(lastRow[nB+1],
		     lastRow[nB+1].score() + extendGapScore(lastRow[nB+1].charString()[1]=='-',
							    nB+1==len(seqB)),
		     seqA[nA]+'-')
      pair = _Square(lastRow[nB],
		     lastRow[nB].score() + matchScore[seqA[nA]==seqB[nB]],
		     seqA[nA]+seqB[nB])
          ## set all scores to zero if they are negative
    squareList = map(lambda i: _Square(i.source(),max([i.score(),0]),i.charString()), [aGap,bGap,pair])
    maxScore = max(map(lambda i: i.score(), squareList))
    for sq in squareList:
    if sq.score()==maxScore: best = sq
      if best.score() > bestScoringSquare.score(): bestScoringSquare = best
      newRow.append(best)
       lastRow = newRow

    lastSquare = bestScoringSquare
      alignmentScore = lastSquare.score()
    aBackwards = []
    bBackwards = []
      while lastSquare.source()!=None and lastSquare.score()>0:
      aBackwards.append(lastSquare.charString()[0])
      bBackwards.append(lastSquare.charString()[1])
      lastSquare = lastSquare.source()
      aBackwards.reverse()
      bBackwards.reverse()
      return Alignment(alignmentScore,''.join(aBackwards),''.join(bBackwards))
*/


