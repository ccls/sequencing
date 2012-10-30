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

#ifndef SCOREDSEQCOLLECTIONBWT_CPP
#define SCOREDSEQCOLLECTIONBWT_CPP

#include "ScoredSeqCollectionBwt.h"

#include "AssemblyException.h"
#include "AlignmentNull.h"
#include "AlignmentUngapped.h"
#include "ScoredSeqFlip.h"
#include <iostream>
#include <omp.h>
#include <algorithm>
using namespace::std;


long ScoredSeqCollectionBwt::_defaultMinScore = 0;

ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, bool penalizeEdgeGaps) : 
  _asi(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(0),
  _penalizeEdgeGaps(penalizeEdgeGaps),
  _usingMaxOverlap(false),
  _maxOverlap(-1) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap) : 
  _asi(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(minOverlap),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap, long maxOverlap) : 
  _asi(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(minOverlap),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DynamicProgrammingAligner * asi) :
  _gapMode(gapped),
  _fractId(asi->getFractId()),
  _minOverlap(asi->getMinOverlap()),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1) {
  _asi = new DynamicProgrammingAligner(asi);
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DynamicProgrammingAligner * asi, long maxOverlap) :
  _gapMode(gapped),
  _fractId(asi->getFractId()),
  _minOverlap(asi->getMinOverlap()),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap) {
  _asi = new DynamicProgrammingAligner(asi);
  constructorHelper(inputSeqs);
}


void ScoredSeqCollectionBwt::constructorHelper(set<ScoredSeq*>* inputSeqs){
  // fill _allSeqs now so that the log of the collection doesn't need to be faced again
  _allSeqs.insert( inputSeqs->begin(), inputSeqs->end() );

  // by default, edge scaling will happen.  it can be disabled to save
  // time during read mapping.
  _edgeScalingEnabled = true;

  // minOverlap cannot be zero for this class; sub-string words must be defined
  // and statistically scaled for this overlap minimum, so a hard bottom line
  // needs to be defined here
  long minOverlapLimit = 4;
  if (_minOverlap < minOverlapLimit){ _minOverlap = minOverlapLimit; }

  // fractId must be at least that expected for random sequence matches
  if (_fractId < 0.25){ _fractId = 0.25; }

  // figure out the windows and make the bins
  //long seqCount = 0;
  long seqLenSum = 0;
  // determine the length of the shortest sequence and whether or not it is
  // smaller than _minOverlap
  long minSeqLength = 0;
  long seqCount = inputSeqs->size();

  set<ScoredSeq*>::iterator seqIt = inputSeqs->begin();
  if (seqIt != inputSeqs->end()){
    minSeqLength = (*seqIt)->size();
    seqLenSum = minSeqLength + 1;
    ++seqIt;
    while (seqIt != inputSeqs->end()){
      long seqSize = (*seqIt)->size();
      if (seqSize < minSeqLength){ minSeqLength = seqSize; }
      seqLenSum += seqSize + 1;
      ++seqIt;
    }
  }
  if (minSeqLength < _minOverlap){ _minFragmentSize = _minOverlap; }
  else { _minFragmentSize = long(float(minSeqLength) * _fractId); }

  // some data that will be useful to have pre-computed
  _maxFractMis = 1.0 - _fractId;
  long matchScore = 1;
  long misScore = -1;
  _scoreAsm = new AlignmentScoreMatrix(matchScore,misScore,misScore,misScore);
  _matchCountAsm = new AlignmentScoreMatrix(1,0,0,0);
  _misCountAsm = new AlignmentScoreMatrix(0,1,1,1); // gaps count as mismatches (terminal gaps will exist)

  // figure out the windows and make the bins
  if (seqCount==0){
    _binSize = 0; // this value is only possible in this case
    _maxBin = 0;
  } else {
    _binSize = (seqLenSum / seqCount) + 1; // +1 to round up
    _maxBin = seqLenSum / _binSize;
  }
  _binToSeqs = new vector<ScoredSeqLocus*>[_maxBin+1];

  // convert the full text to a numeric array, filling in the bins
  if (seqLenSum==0){ seqLenSum = 1; } // at least an end char will be added
  _textSize = seqLenSum;

  // a local array of text length
  long* textAsNums = new long[seqLenSum];

  long currentIndex = 0;
  _numContigLoci = seqCount;
  if (_numContigLoci > 0){
    _orderedContigLoci = new ScoredSeqLocus*[ _numContigLoci ];
  }

  long currentOffset = 0;
  for (set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it){
    // add the sequence, preceeded by a separator if necessary
    if ( it != _allSeqs.begin() ){
      char* separator = new char[2];
      separator[0] = '&';
      separator[1] = '\0';
      burrowsWheelerUtilities::makeNumString( separator, long(1), currentOffset, textAsNums );
      delete [] separator;
      currentOffset++;
    }
    ScoredSeq* seq = (*it);
    char* seqString = seq->getSeq('+');
    burrowsWheelerUtilities::makeNumString( seqString, seq->size(), currentOffset, textAsNums );
    delete [] seqString;

    // fill the overlapping bins with this ScoredSeq in Locus form
    long firstBin = currentOffset / _binSize;
    long lastBinP1 = ( currentOffset + seq->size() ) / _binSize + 1;
    ScoredSeqLocus* seqAsLocus = new ScoredSeqLocus(seq, currentOffset, currentOffset + seq->size() - 1, currentIndex);
    for (long bin = firstBin; bin < lastBinP1; ++bin){ _binToSeqs[bin].push_back( seqAsLocus ); }

    _orderedContigLoci[currentIndex] = seqAsLocus;
    currentIndex++;

    currentOffset += seq->size();
  }

  // end the sequence with an end char (even if there were no sequence entries)
  char* endChar = new char[2];
  endChar[0] = '#';
  endChar[1] = '\0';
  burrowsWheelerUtilities::makeNumString( endChar, long(1), currentOffset, textAsNums );
  delete [] endChar;

  // make the burrows-wheeler transform
  _sortedSuffixes = new long[seqLenSum];
  long alphaSize = burrowsWheelerUtilities::getAlphabetSize();
  burrowsWheelerUtilities::sortSuffixes( textAsNums, _sortedSuffixes, seqLenSum, alphaSize );
  _bwTransform = new long[seqLenSum];
  burrowsWheelerUtilities::bwTransform(textAsNums, _sortedSuffixes, seqLenSum, _bwTransform);

  delete [] textAsNums;

  // make the accessory data structures for alignment
  _occCounts = new long[seqLenSum];
  _tableC = burrowsWheelerUtilities::fillBwtHelperArrays(_bwTransform, seqLenSum, alphaSize, _occCounts);
  //OK();
}



void ScoredSeqCollectionBwt::OK(){
  if ( _fractId > 1 ){
    throw AssemblyException::ArgError("fraction ID is a fraction; cannot be greater than 1");
  }
  if (_allSeqs.size()==0){
    if (_textSize != 1){ throw AssemblyException::LogicError("BWT text len is 1 if there are no seqs"); }
  } else {
    long seqLenSum = 0;
    for (set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it){
      seqLenSum += (*it)->size() + 1; // +1 because of the separator/end char
    }
    if (_textSize != seqLenSum){ throw AssemblyException::LogicError("BWT text len is wrong"); }
  }
}



ScoredSeqCollectionBwt::~ScoredSeqCollectionBwt(){
  delete _asi;
  delete _matchCountAsm;
  delete _misCountAsm;
  delete _scoreAsm;
  delete [] _sortedSuffixes;
  delete [] _bwTransform;
  delete [] _occCounts;
  delete [] _tableC;
  for (long n = 0; n < _numContigLoci; ++n){ delete _orderedContigLoci[n]; }
  if (_numContigLoci > 0){ delete [] _orderedContigLoci; }

  delete [] _binToSeqs;
}




ScoredSeqCollectionBwt * ScoredSeqCollectionBwt::copy(){
  // uses a private constructor to efficiently copy stuff
  throw AssemblyException::ImplementationError("not implemented");
}

float ScoredSeqCollectionBwt::getFractId(){ return _fractId; }
long ScoredSeqCollectionBwt::getMinOverlap(){ return _minOverlap; }



// DEFAULT: edge scaling is enabled
void ScoredSeqCollectionBwt::enableEdgeScaling(){ _edgeScalingEnabled = true; }
void ScoredSeqCollectionBwt::disableEdgeScaling(){ _edgeScalingEnabled = false; }



bool ScoredSeqCollectionBwt::contains(ScoredSeq* s){
  //OK();
  return ( _allSeqs.find(s) != _allSeqs.end() );
}



void ScoredSeqCollectionBwt::getSeqs( set<ScoredSeq*>* allSeqs ){
  //OK();
  allSeqs->insert( _allSeqs.begin(), _allSeqs.end() );
}



// UNTHREADED VERSIONS
// min overlap is specified - these four methods cascade down to the bottom one
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, long minOverlap,
					MinOvlStringency ovlStringency){
  getMatches(matches, seq, '.', minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, long minOverlap,
					MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  getMatches(matches, seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq,
					MatchSeqTest* matchTest, long minOverlap,
					MinOvlStringency ovlStringency){
  getMatches(matches, seq, '.', matchTest, minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
					MatchSeqTest* matchTest, long minOverlap,
					MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense,
		       matchTest, minOverlap,
		       semiGlobal, allMatches, ovlStringency, notThreaded);
  matches->insert( matches->end(), tempMatches.begin(), tempMatches.end() );
}



// THREADED VERSIONS
// these four methods cascade down to the one on the bottom
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getMatches(numQueries, matchesArray, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** matchTestArray = new MatchSeqTest*[  numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ matchTestArray[n] = seqTest; }
  getMatches(numQueries, matchesArray, seqArray, sense, matchTestArray, minOvlArray, threadedness, ovlStringency);
  delete seqTest;
  delete [] matchTestArray;
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
					MatchSeqTest** matchTestArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getMatches(numQueries, matchesArray, seqArray, '.', matchTestArray, minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
					MatchSeqTest** matchTestArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, matchTestArray, minOvlArray,
		     semiGlobal, allMatches, ovlStringency, threadedness);

  for (long n = 0; n < numQueries; ++n){
    matchesArray[n]->insert( matchesArray[n]->end(), tempMatches[n]->begin(), tempMatches[n]->end() );
    delete tempMatches[n];
  }
  delete [] tempMatches;
}







void ScoredSeqCollectionBwt::setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
						     ScoredSeq** seqFlipArray){
  // get the numeric representation of the query
  seqFlipArray[index] = ScoredSeqFlip::getFlip(seqArray[index], sense);
}
void ScoredSeqCollectionBwt::setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
						     AlignmentMaker** alMakerArray, ScoredSeq** seqFlipArray){
  if (_penalizeEdgeGaps){
    alMakerArray[index] = new AlMakerPenalizeEdgeA(_scoreAsm, _defaultMinScore, _maxFractMis); // value in alMaker
  } else {
    alMakerArray[index] = new AlMakerNoEdge(_scoreAsm, _defaultMinScore, _maxFractMis); // value in alMaker
  }
  setupPrivateMultiHelper(index, seqArray, sense,
			  seqFlipArray);
}
long ScoredSeqCollectionBwt::setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
						     long* minOvlArray, MatchSeqTest** matchTestArray,
						     ScoredSeq** seqFlipArray,
						     OffsetTracker*** conToOffByQuery, long* otCountByQuery, AlignmentType alType){
  // for each matched contig, keeps track of the num hits with any given offset
  conToOffByQuery[index] = getOffsets(seqArray[index],sense,minOvlArray[index],matchTestArray[index],alType);
  // count the number of offset trackers for each query
  long hitIndex = 0;
  while (conToOffByQuery[index][hitIndex] != NULL){ ++hitIndex; }
  otCountByQuery[index] = hitIndex;
  setupPrivateMultiHelper(index, seqArray, sense,
			  seqFlipArray);
  return hitIndex;
}
long ScoredSeqCollectionBwt::setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
						     long* minOvlArray, MatchSeqTest** matchTestArray,
						     AlignmentMaker** alMakerArray, ScoredSeq** seqFlipArray,
						     OffsetTracker*** conToOffByQuery, long* otCountByQuery, AlignmentType alType){
  // for each matched contig, keeps track of the num hits with any given offset
  conToOffByQuery[index] = getOffsets(seqArray[index],sense,minOvlArray[index],matchTestArray[index],alType);
  // count the number of offset trackers for each query
  long hitIndex = 0;
  while (conToOffByQuery[index][hitIndex] != NULL){ ++hitIndex; }
  otCountByQuery[index] = hitIndex;
  setupPrivateMultiHelper(index, seqArray, sense,
			  alMakerArray, seqFlipArray);
  return hitIndex;
}




long ScoredSeqCollectionBwt::makeAlignmentHelper(OffsetTracker* tracker, ScoredSeq* seq, vector<Alignment*>* matches,
						 long minOverlap, AlignmentMaker* alMaker, ScoredSeq* seqFlip, char sense,
						 long bestScoreSoFar, bool onlyBestAlignment, GapMode gapMode,
						 MatchesSought matchesSought){
  bool justFirstMatch = (matchesSought == firstMatch);
  ScoredSeq* contig = tracker->targetContig();
  long numOffsets = tracker->numOffsets();
  OffsetAndQueryCarrier** allOffsets = tracker->getOffsets(gapMode == ungapped);

  OffsetAndQueryCarrier** adequateOffsets = new OffsetAndQueryCarrier*[ numOffsets + 1 ];
  long numAdequateOffsets = 0;
  for (long allIndex = 0; allIndex < numOffsets; ++allIndex){
    // determine how long the alignment will be??  i will need to do that to make sure
    // that only legit alignments are kept.
    long alLength = getOverlapFromOffset(allOffsets[allIndex]->_offset, seq->size(), contig->size());
    if (alLength < 0){ throw AssemblyException::LogicError("SSCBwt::getMatches, alLength should be > 0"); }
    if ( (_penalizeEdgeGaps or alLength >= minOverlap) and ( (! _usingMaxOverlap) or alLength <= _maxOverlap) ){
      adequateOffsets[numAdequateOffsets] = allOffsets[allIndex];
      ++numAdequateOffsets;
    }
  }

  // this will be different for gapped versus ungapped alignments
  if (gapMode == ungapped){ // UNGAPPED ALIGNMENT
    vector<Alignment*> tempMatches;
    long inputBestScore = bestScoreSoFar;

    // this will reduce the number of seq replacements and (more time-consuming) RcFlips
    // that need to be performed on alignments (only in the ungapped case; deal with gapped differently?)
    for (long offN = 0; offN < numAdequateOffsets; ++offN){
      Alignment* alignCand = alMaker->makeAlignment(seqFlip, contig, '+', adequateOffsets[offN]);
      if ( alignCand->isNull() ){ delete alignCand; }
      else {
	if (justFirstMatch){ offN = numAdequateOffsets; }
	bool returnMatch = true;
	if (onlyBestAlignment){
	  long newScore = alMaker->getScore(alignCand);
	  if (newScore < bestScoreSoFar){ returnMatch = false; }
	  // if they are equal, then the match is just added
	  else if (newScore > bestScoreSoFar){
	    bestScoreSoFar = newScore;
	    for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
	    tempMatches.clear();
	    alMaker->setMinScore(newScore);
	  }
	}
	if (returnMatch){ tempMatches.push_back(alignCand); }
	else { delete alignCand; }
      }
    }
    // now compare to and deal with the input set's contents
    if (onlyBestAlignment and bestScoreSoFar > inputBestScore){
      for (vector<Alignment*>::iterator it = matches->begin(); it != matches->end(); ++it){ delete *it; }
      matches->clear();
    }
    // and modify the output alignments appropriately according to the seq's orientation
    matches->insert(matches->end(), tempMatches.begin(), tempMatches.end());

  } else { // GAPPED ALIGNMENT
    if (numAdequateOffsets > 0){
      DynamicProgrammingAligner * localAsi = _asi->minScoreCopy(bestScoreSoFar);
      long* justOffsetNums = new long[ numAdequateOffsets ];
      for (long offN = 0; offN < numAdequateOffsets; ++offN){ justOffsetNums[offN] = adequateOffsets[offN]->_offset; }
      Alignment* alignCand = localAsi->align(seqFlip,contig,'+',justOffsetNums,numAdequateOffsets);
      delete [] justOffsetNums;
      if ( alignCand->isNull() ){ delete alignCand; }
      else {
	bool returnMatch = true;
	if (onlyBestAlignment){
	  long newScore = alignCand->score(_scoreAsm,_penalizeEdgeGaps);
	  if (newScore < bestScoreSoFar){ returnMatch = false; }
	  // if they are equal, then the match is just added
	  else if (newScore > bestScoreSoFar){
	    bestScoreSoFar = newScore;
	    for (vector<Alignment*>::iterator it = matches->begin(); it != matches->end(); ++it){ delete *it; }
	    matches->clear();
	  }
	}
	if (returnMatch){ matches->push_back(alignCand); }
	else { delete alignCand; }
      }
      delete localAsi;
    }
  }

  for (long allIndex = 0; allIndex < numOffsets; ++allIndex){
    for (long wqN = 0; wqN < allOffsets[allIndex]->_numQueries; ++wqN){ delete allOffsets[allIndex]->_queries[wqN]; }
    delete allOffsets[allIndex];
  }
  delete [] allOffsets;
  delete [] adequateOffsets;

  return bestScoreSoFar;
}

/*
void ScoredSeqCollectionBwt::getBestFullMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* seqTest){
  MatchSeqTest* useSeqTest;
  if (seqTest == NULL){ useSeqTest = MatchSeqTest::getNullTest(); }
  else { useSeqTest = seqTest; }
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, useSeqTest, _minFragmentSize, fullAlignment, bestMatches, hardMinOvl, notThreaded);
  if (seqTest == NULL){ delete useSeqTest; }
  matches->insert( matches->end(), tempMatches.begin(), tempMatches.end() );
}

// THREADED
void ScoredSeqCollectionBwt::getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
						AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  getBestFullMatches(numQueries, matchesArray, seqArray, '.', threadedness, seqTestArray);
}
void ScoredSeqCollectionBwt::getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
						AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  MatchSeqTest** useSeqTestArray;
  MatchSeqTest* dummyTest = NULL;
  if (seqTestArray == NULL){
    dummyTest = MatchSeqTest::getNullTest();
    useSeqTestArray = new MatchSeqTest*[ numQueries+1 ];
    for (long n = 0; n < numQueries; ++n){ useSeqTestArray[n] = dummyTest; }
  } else { useSeqTestArray = seqTestArray; }
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  long* minOvlArray = new long[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ minOvlArray[n] = _minFragmentSize; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, useSeqTestArray, minOvlArray, fullAlignment, firstMatch, hardMinOvl, threadedness);

  delete [] minOvlArray;
  for (long n = 0; n < numQueries; ++n){
    matchesArray[n]->insert( matchesArray[n]->end(), tempMatches[n]->begin(), tempMatches[n]->end() );
    delete tempMatches[n];
  }
  delete [] tempMatches;
  if (seqTestArray == NULL){
    delete dummyTest;
    delete useSeqTestArray;
  }
}
*/



// min overlap is specified - these four methods cascade down
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq,
					    long minOverlap, MinOvlStringency ovlStringency){
  getBestMatches(matches, seq, '.', minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
					    long minOverlap, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  getBestMatches(matches, seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* seqTest,
					    long minOverlap, MinOvlStringency ovlStringency){
  getBestMatches(matches, seq, '.', seqTest, minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* seqTest,
					    long minOverlap, MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, seqTest, minOverlap, semiGlobal, bestMatches, ovlStringency, notThreaded);
  matches->insert( matches->end(), tempMatches.begin(), tempMatches.end() );
}

// THREADED VERSIONS
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getBestMatches(numQueries, matchesArray, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** seqTestArray = new MatchSeqTest*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ seqTestArray[n] = seqTest; }
  getBestMatches(numQueries, matchesArray, seqArray, sense, seqTestArray, minOvlArray, threadedness, ovlStringency);
  delete seqTest;
  delete [] seqTestArray;
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** seqTestArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getBestMatches(numQueries, matchesArray, seqArray, '.', seqTestArray, minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, seqTestArray, minOvlArray, semiGlobal, bestMatches, ovlStringency, threadedness);

  for (long n = 0; n < numQueries; ++n){
    matchesArray[n]->insert( matchesArray[n]->end(), tempMatches[n]->begin(), tempMatches[n]->end() );
    delete tempMatches[n];
  }
  delete [] tempMatches;
}



// THE BEGINNING OF HASMATCH
bool ScoredSeqCollectionBwt::hasFullMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest){
  MatchSeqTest* useSeqTest;
  if (seqTest == NULL){ useSeqTest = MatchSeqTest::getNullTest(); }
  else { useSeqTest = seqTest; }
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, useSeqTest, _minFragmentSize, fullAlignment, firstMatch, hardMinOvl, notThreaded);
  bool hasMatch = (tempMatches.size() > 0);
  for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
  if (seqTest == NULL){ delete useSeqTest; }
  return hasMatch;
}


bool* ScoredSeqCollectionBwt::hasFullMatch(long numQueries, ScoredSeq** seqArray,
					  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  return hasFullMatch(numQueries, seqArray, '.', threadedness, seqTestArray);
}

bool* ScoredSeqCollectionBwt::hasFullMatch(long numQueries, ScoredSeq** seqArray, char sense,
					  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  MatchSeqTest** useSeqTestArray;
  MatchSeqTest* dummyTest = NULL;
  if (seqTestArray == NULL){
    dummyTest = MatchSeqTest::getNullTest();
    useSeqTestArray = new MatchSeqTest*[ numQueries+1 ];
    for (long n = 0; n < numQueries; ++n){ useSeqTestArray[n] = dummyTest; }
  } else { useSeqTestArray = seqTestArray; }
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  long* minOvlArray = new long[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ minOvlArray[n] = _minFragmentSize; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, useSeqTestArray, minOvlArray, fullAlignment, firstMatch, hardMinOvl, threadedness);

  delete [] minOvlArray;

  bool* answers = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    answers[n] = (tempMatches[n]->size() > 0);
    for (vector<Alignment*>::iterator it = tempMatches[n]->begin(); it != tempMatches[n]->end(); ++it){ delete *it; }
  }
  for (long n = 0; n < numQueries; ++n){ delete tempMatches[n]; }
  delete [] tempMatches;
  if (seqTestArray == NULL){
    delete dummyTest;
    delete [] useSeqTestArray;
  }
  return answers;
}



// min overlap is specified - these four methods cascade down
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, long minOverlap, MinOvlStringency ovlStringency){
  return hasMatch(seq, '.', minOverlap, ovlStringency);
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, char sense, long minOverlap, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  bool returnVal = hasMatch(seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
  return returnVal;
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency){
  return hasMatch(seq, '.', seqTest, minOverlap, ovlStringency);
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, seqTest, minOverlap, semiGlobal, firstMatch, ovlStringency, notThreaded);
  bool hasMatch = (tempMatches.size() > 0);
  for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
  return hasMatch;
}


// THREADED
// min overlap is specified - these four methods cascade down
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  return hasMatch(numQueries, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, char sense, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** seqTestArray = new MatchSeqTest*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ seqTestArray[n] = seqTest; }
  bool* returnVal = hasMatch(numQueries, seqArray, sense, seqTestArray, minOvlArray, threadedness, ovlStringency);
  delete [] seqTestArray;
  delete seqTest;
  return returnVal;
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, MatchSeqTest** seqTestArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  return hasMatch(numQueries, seqArray, '.', seqTestArray, minOvlArray, threadedness, ovlStringency);
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, seqTestArray, minOvlArray, semiGlobal, firstMatch, ovlStringency, threadedness);

  bool* answers = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    answers[n] = (tempMatches[n]->size() > 0);
    for (vector<Alignment*>::iterator it = tempMatches[n]->begin(); it != tempMatches[n]->end(); ++it){ delete *it; }
  }
  for (long n = 0; n < numQueries; ++n){ delete tempMatches[n]; }
  delete [] tempMatches;
  return answers;
}




// THE END OF HASMATCH


// ASSUMES: matches is empty; otherwise, it will be modified
void ScoredSeqCollectionBwt::runSearchesConverter(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
						  MatchSeqTest* seqTest, long minOverlap,
						  AlignmentType alType, MatchesSought matchesSought,
						  MinOvlStringency ovlStringency, AlignmentThreadedness threadedness){

  vector<Alignment*>** matchesArray = new vector<Alignment*>*[ 1 ];
  matchesArray[ 0 ] = matches;

  ScoredSeq** seqArray = new ScoredSeq*[ 1 ];
  seqArray[0] = seq;
  MatchSeqTest** seqTestArray = new MatchSeqTest*[ 1 ];
  seqTestArray[0] = seqTest;

  long* minOvlArray = new long[ 1 ];
  minOvlArray[0] = minOverlap;

  runSearchesPrivate(1, matchesArray, seqArray, sense,
		     seqTestArray, minOvlArray,
		     alType, matchesSought, ovlStringency, threadedness);

  delete [] matchesArray;
  delete [] seqArray;
  delete [] seqTestArray;
  delete [] minOvlArray;
}



// ASSUMES: matches is empty; otherwise, it will be modified
void ScoredSeqCollectionBwt::runSearchesPrivate(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
						MatchSeqTest** seqTestArray, long* minOvlArray,
						AlignmentType alType, MatchesSought matchesSought,
						MinOvlStringency ovlStringency, AlignmentThreadedness threadedness){
  int numSenses;
  char senses[2];
  if (sense == '.'){
    numSenses = 2;
    senses[0] = '+';
    senses[1] = '-';
  } else if (sense == '+' or sense == '-'){
    numSenses = 1;
    senses[0] = sense;
  } else {
    cerr << '"' << sense << '"' << endl;
    throw AssemblyException::ArgError("SSCBwt:getBestMatchesPrivate bad sense char when '.' is allowed.");
  }

  long* correctedMinOvls = new long[ numQueries+1 ];
  if (ovlStringency == softMinOvl){
    for (long n = 0; n < numQueries; ++n){ correctedMinOvls[n] = long(float(minOvlArray[n]) * _fractId); }
  } else {
    for (long n = 0; n < numQueries; ++n){ correctedMinOvls[n] = minOvlArray[n]; }
  }

  OffsetTracker*** otArraysBySense[2];
  for (int sn = 0; sn < numSenses; ++sn){
    otArraysBySense[sn] = new OffsetTracker**[ numQueries+1 ];
    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long n2 = 0; n2 < numQueries; ++n2){
	otArraysBySense[sn][n2] = getOffsets(seqArray[n2],senses[sn],minOvlArray[n2],seqTestArray[n2],alType);
      }
    } else {
      for (long n2 = 0; n2 < numQueries; ++n2){
	otArraysBySense[sn][n2] = getOffsets(seqArray[n2],senses[sn],minOvlArray[n2],seqTestArray[n2],alType);
      }
    }
  }

  long* bestScoreArray = new long[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ bestScoreArray[n] = _defaultMinScore; }

  if ((matchesSought != allMatches) or (_gapMode == ungapped)){
    for (int n = 0; n < numSenses; ++n){
      runAlignmentsHelper(numQueries, otArraysBySense[n], matchesArray, seqArray, senses[n], ungapped, seqTestArray,
			  correctedMinOvls, bestScoreArray, matchesSought, threadedness);
    }
  }
  if (_gapMode == gapped){
    for (long n = 0; n < numQueries; ++n){
      for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
      matchesArray[n]->clear();
    }
    for (int n = 0; n < numSenses; ++n){
      runAlignmentsHelper(numQueries, otArraysBySense[n], matchesArray, seqArray, senses[n], gapped, seqTestArray,
			  correctedMinOvls, bestScoreArray, matchesSought, threadedness);
    }
  }

  // delete the offset tracker arrays and their contents
  for (int sn = 0; sn < numSenses; ++sn){
    for (int qn = 0; qn < numQueries; ++qn){
      long hitIndex = 0;
      while (otArraysBySense[sn][qn][hitIndex] != NULL){
	delete otArraysBySense[sn][qn][hitIndex];
	hitIndex++;
      }
      delete [] otArraysBySense[sn][qn];
    }
    delete [] otArraysBySense[sn];
  }
  //OK();

  delete [] correctedMinOvls;
  delete [] bestScoreArray;
}



// this code looks redundant because it is optimized for speed
inline long ScoredSeqCollectionBwt::getOverlapFromOffset(long offset, long seqSizeA, long seqSizeB){
  if (offset < 0){
    if (offset + seqSizeB >= seqSizeA){
      return seqSizeA;
    } else {
      return offset + seqSizeB;
    }
  } else {
    if (offset + seqSizeB >= seqSizeA){
      return seqSizeA - offset;
    } else {
      return seqSizeB;
    }
  }
}


// REQUIRES: all members of "matches" (if any) have the same score
void ScoredSeqCollectionBwt::runAlignmentsHelper(long numQueries, OffsetTracker*** contigToOffsetsArray,
						 vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, GapMode gapMode,
						 MatchSeqTest** seqTestArray, long* minOvlArray, long* bestScoreArray,
						 MatchesSought matchesSought, AlignmentThreadedness threadedness){
  bufferSeqs();

  bool justFirstMatch = (matchesSought == firstMatch);
  bool getAllMatches = (matchesSought == allMatches);

  long* initBestScoreArray = new long[ numQueries+1 ];
  long* otCountByQuery = new long[ numQueries+1 ];
  long totalOtCount = 0;

  bool* matchesStillSought = new bool[ numQueries+1 ];
  ScoredSeq** seqFlipArray = new ScoredSeq*[ numQueries+1 ];

  // this is so that all of the alignments produced here can be modified in the end
  vector<Alignment*>* tempMatchesArray = new vector<Alignment*>[ numQueries+1 ];
  // this is threaded by query, the best I can do at this point, but sets up for threading by hit
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numQueries; ++n){
      initBestScoreArray[n] = bestScoreArray[n];
      matchesStillSought[n] = ( (! justFirstMatch) or matchesArray[n]->size() == 0 );
      otCountByQuery[n] = 0;
      if (matchesStillSought[n]){
	while (contigToOffsetsArray[n][otCountByQuery[n]] != NULL){ ++otCountByQuery[n]; }
      }
      setupPrivateMultiHelper(n, seqArray, sense, seqFlipArray);
      #pragma omp atomic
      totalOtCount += otCountByQuery[n];
    }
  } else {
    for (long n = 0; n < numQueries; ++n){
      initBestScoreArray[n] = bestScoreArray[n];
      matchesStillSought[n] = ( (! justFirstMatch) or matchesArray[n]->size() == 0 );
      otCountByQuery[n] = 0;
      while (contigToOffsetsArray[n][otCountByQuery[n]] != NULL){ ++otCountByQuery[n]; }
      setupPrivateMultiHelper(n, seqArray, sense, seqFlipArray);
      totalOtCount += otCountByQuery[n];
    }
  }


  // this sets things up for threading across multiple hits per query by defining
  // the query and the target for each array element across an array that encompasses
  // all of the hits for all queries
  long totalOtCountP1 = totalOtCount + 1;
  long* useQueryIndex = new long[ totalOtCountP1 ];
  OffsetTracker** contigToOffsets = new OffsetTracker*[ totalOtCountP1 ];
  ScoredSeq** useSeqArray = new ScoredSeq*[ totalOtCountP1 ];
  vector<Alignment*>** useMatchesArray = new vector<Alignment*>*[ totalOtCountP1 ];
  long* useMinOvlArray = new long[ totalOtCountP1 ];
  ScoredSeq** useSeqFlipArray = new ScoredSeq*[ totalOtCountP1 ];
  long totalIndex = 0;


  // figure out the max number of hits so that hits from the same query can be maximally spread out
  long maxOtCount = 0;
  for (long n = 0; n < numQueries; ++n){
    if (otCountByQuery[n] > maxOtCount){ maxOtCount = otCountByQuery[n]; }
  }
  // i am unsure about which is the better way to order these jobs (n outer or n2 outer) for threading
  for (long n2 = 0; n2 < maxOtCount; ++n2){
    for (long n = 0; n < numQueries; ++n){
      if (n2 < otCountByQuery[n]){
	useQueryIndex[totalIndex] = n;
	contigToOffsets[totalIndex] = contigToOffsetsArray[n][n2];
	useSeqArray[totalIndex] = seqArray[n];
	useMatchesArray[totalIndex] = &tempMatchesArray[n];
	useMinOvlArray[totalIndex] = minOvlArray[n];
	useSeqFlipArray[totalIndex] = seqFlipArray[n];
	++totalIndex;
      }
    }
  }
  delete [] otCountByQuery;

  // make alignments
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long hitIndex = 0; hitIndex < totalIndex; ++hitIndex){
      // set up the alMaker to reflect the most up-to-date min score
      long oldBestScore;
      bool matchSought;
      #pragma omp critical (SSCBwtScoreStatus)
      {
	oldBestScore = bestScoreArray[useQueryIndex[hitIndex]];
	matchSought = matchesStillSought[useQueryIndex[hitIndex]];
      }
      if (matchSought){
	AlignmentMaker * alMaker = makeAlMaker(gapMode, oldBestScore);
	vector<Alignment*> tempMatches;
	long newBestScore = makeAlignmentHelper(contigToOffsets[hitIndex], useSeqArray[hitIndex], &tempMatches,
						useMinOvlArray[hitIndex], alMaker, useSeqFlipArray[hitIndex],
						sense, oldBestScore, true, gapMode, matchesSought);
        #pragma omp critical (SSCBwtScoreStatus)
	{
	dealWithMatches(&tempMatches, useMatchesArray[hitIndex],
			getAllMatches, justFirstMatch,
			useQueryIndex[hitIndex], bestScoreArray, newBestScore, matchesStillSought);
	}
	delete alMaker;
      }
    }
  } else {
    for (long hitIndex = 0; hitIndex < totalIndex; ++hitIndex){
      // set up the alMaker to reflect the most up-to-date min score
      long oldBestScore = bestScoreArray[useQueryIndex[hitIndex]];
      bool matchSought = matchesStillSought[useQueryIndex[hitIndex]];
      if (matchSought){
	AlignmentMaker * alMaker = makeAlMaker(gapMode, oldBestScore);
	vector<Alignment*> tempMatches;
	long newBestScore = makeAlignmentHelper(contigToOffsets[hitIndex], useSeqArray[hitIndex], &tempMatches,
						useMinOvlArray[hitIndex], alMaker, useSeqFlipArray[hitIndex],
						sense, oldBestScore, true, gapMode, matchesSought);
	dealWithMatches(&tempMatches, useMatchesArray[hitIndex],
			getAllMatches, justFirstMatch,
			useQueryIndex[hitIndex], bestScoreArray, newBestScore, matchesStillSought);
	delete alMaker;
      }
    }
  }

  // replace the sequences in the alignment
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numQueries; ++n){
      if (bestScoreArray[n] > initBestScoreArray[n]){
	for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
	matchesArray[n]->clear();
      }
    }
  } else {
    for (long n = 0; n < numQueries; ++n){
      if (bestScoreArray[n] > initBestScoreArray[n]){
	for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
	matchesArray[n]->clear();
      }
    }
  }
  delete [] initBestScoreArray;

  if (sense == '+'){
    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  (*it)->seqReplace(seqArray[n], (*it)->seqB());
	}
	matchesArray[n]->insert(matchesArray[n]->end(), tempMatchesArray[n].begin(), tempMatchesArray[n].end());
      }
    } else {
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  (*it)->seqReplace(seqArray[n], (*it)->seqB());
	}
	matchesArray[n]->insert(matchesArray[n]->end(), tempMatchesArray[n].begin(), tempMatchesArray[n].end());
      }
    }
  } else if (sense == '-'){
    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  matchesArray[n]->push_back( (*it)->copyRcSeqA(seqArray[n]) );
	  delete *it;
	}
      }
    } else {
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  matchesArray[n]->push_back( (*it)->copyRcSeqA(seqArray[n]) );
	  delete *it;
	}
      }
    }
  } else { throw AssemblyException::ArgError("SSCBWT:alignment method under construction, bad sense char"); }

  delete [] tempMatchesArray;
  delete [] matchesStillSought;

  delete [] useQueryIndex;
  delete [] contigToOffsets;
  delete [] useSeqArray;
  delete [] useMatchesArray;
  delete [] useMinOvlArray;
  delete [] useSeqFlipArray;

  for (long n = 0; n < numQueries; ++n){ delete seqFlipArray[n]; }
  delete [] seqFlipArray;
}



ScoredSeqCollectionBwt::AlignmentMaker* ScoredSeqCollectionBwt::makeAlMaker(GapMode gapMode, long bestScoreSoFar){
  if (gapMode == gapped){ return NULL; }
  else if (_penalizeEdgeGaps){ return new AlMakerPenalizeEdgeA(_scoreAsm, bestScoreSoFar, _maxFractMis); }
  else { return new AlMakerNoEdge(_scoreAsm, bestScoreSoFar, _maxFractMis); }
}

void ScoredSeqCollectionBwt::dealWithMatches(vector<Alignment*>* tempMatches, vector<Alignment*>* realMatches,
					     bool getAllMatches, bool justFirstMatch,
					     long index, long* bestScoreArray, long newBestScore, bool* matchSoughtArray){
  if (getAllMatches){
    realMatches->insert(realMatches->end(), tempMatches->begin(), tempMatches->end());
  } else {
    if (newBestScore < bestScoreArray[index]){
      for (vector<Alignment*>::iterator it = tempMatches->begin(); it != tempMatches->end(); ++it){ delete *it; }
    } else {
      if (newBestScore > bestScoreArray[index]){
	for (vector<Alignment*>::iterator it = realMatches->begin(); it != realMatches->end(); ++it){ delete *it; }
	realMatches->clear();
	bestScoreArray[index] = newBestScore;
      }
      realMatches->insert(realMatches->end(), tempMatches->begin(), tempMatches->end());
      if (justFirstMatch and tempMatches->size() > 0){ matchSoughtArray[index] = false; }
    }
  }
}

bool ScoredSeqCollectionBwt::SortWqByLength::operator() (WordQuery* wqA, WordQuery* wqB){
  return wqA->size() < wqB->size();
}


ScoredSeqCollectionBwt::OffsetTracker** ScoredSeqCollectionBwt::getOffsets(ScoredSeq* seq, char sense, long minOverlap,
									   MatchSeqTest* matchTest, AlignmentType alType){
  // this is what I am collecting now so that I can sort the return values
  map<ScoredSeq*,long> hitToIndex;
  long currentIndex = 0;
  OffsetTracker** hitArray = new OffsetTracker*[ _allSeqs.size() ];

  long querySize = seq->size();
  char* seqString = seq->getSeq(sense);  
  // get the size and position of the sub-strings for which
  // matches will be sought
  vector<WordQuery*> sortedWordQueries;
  if (alType == semiGlobal){
    getSubtextWordQueries(seq, seqString, &sortedWordQueries, sense, minOverlap);
  } else if (alType == fullAlignment){
    getSubtextWordQueriesFullAlignment(seq, seqString, &sortedWordQueries, sense, minOverlap);
  } else {
    throw AssemblyException::ArgError("SSCBwt:getOffsets got a bad alType");
  }
  // i need to create this array, but i don't need to fill it if there are no good word queries
  long* queryNumArray = new long[querySize];
  if (sortedWordQueries.size() > 0){ burrowsWheelerUtilities::makeNumString( seqString, querySize, 0, queryNumArray ); }
  delete [] seqString;
  if (sortedWordQueries.size() > 0){

    // create the length-sorted version of this collection, shortest to longest
    SortWqByLength wqSorter; // sort by size, shortest to longest
    sort( sortedWordQueries.begin(), sortedWordQueries.end(), wqSorter );

    // here is the block data structure that will be used to keep track of redundant queries
    WordQueryBinOrganizer* wqbo = new WordQueryBinOrganizer(&sortedWordQueries);

    long totalQueryCount = wqbo->numQueries();
    for (long queryIndex = 0; queryIndex < totalQueryCount; ++queryIndex){

      // make sure that the query is still valid before doing the search
      if (wqbo->queryStillValid( queryIndex )){
	WordQuery* wordQuery = wqbo->getQuery(queryIndex);

	// these will get used a lot
	long stSize = wordQuery->size();
	long startNum = wordQuery->start();

	// remove the query from the collection
	//wqbo->removeQuery( queryIndex );

	// do the search with the sub-query word
	// deal with possible redundant matchLoci if multiple words match
	long* matchArray = burrowsWheelerUtilities::bwtCountAndFind(queryNumArray, startNum, stSize,
								    _bwTransform, _sortedSuffixes, _textSize,
								    _tableC, _occCounts);

	// keep track of the starts AND how many non-overlapping windows have matched
	long numMatches = matchArray[0];
	if (numMatches == 0){ wqbo->removeSupersetQueries(queryIndex); }
	for (long n = 1; n <= numMatches; ++n){
	  long hitStart = matchArray[n];
	  ScoredSeqLocus* contigLocus = _orderedContigLoci[ getOverlapContigIndex(hitStart, hitStart + stSize - 1) ];
	  ScoredSeq* contig = contigLocus->_seq;
	  if (seq != contig and matchTest->seqIsOk(contig) ){ // don't align to self!!!
	    // determine the offset for the alignment of the query to this seq;
	    // number is the start position of target in seq (i.e. a seq coord)
	    long contigStart = hitStart - contigLocus->_start;

	    // make sure that the hit is on an appropriate region of the target contig
	    if ( wordQuery->acceptableOverlap(contigStart, contig) ){
	      long offset = startNum - contigStart;
	      long contigIndex;
	      map<ScoredSeq*,long>::iterator hitIt = hitToIndex.find( contig );
	      if (hitIt == hitToIndex.end()){
		hitToIndex.insert( pair<ScoredSeq*,long>(contig,currentIndex) );
		contigIndex = currentIndex;
		hitArray[contigIndex] = new OffsetTracker(contig);
		currentIndex++;
	      } else {
		contigIndex = hitIt->second;
	      }
	      hitArray[contigIndex]->addOffset( offset, wordQuery, numMatches );
	    }
	  }
	}
	delete [] matchArray;
      }
    }
    delete wqbo;
    // delete them afterwards, separately
    for (vector<WordQuery*>::iterator it = sortedWordQueries.begin(); it != sortedWordQueries.end(); ++it){ delete *it; }
  }

  OffsetTracker** sortedHitWithOffsets = new OffsetTracker*[ currentIndex + 1 ];
  // the last index is a pair of nulls so that I don't have to return a length
  sortedHitWithOffsets[currentIndex] = NULL;

  // OffsetTrackers will be sorted, but that only needs to be gone through if there ARE OffsetTrackers
  if (currentIndex > 0){

    // now sort by the seed match counts.  first, figure out how many of each number of counts there are
    map<long,long> hitCountsToCount;
    for (long n = 0; n < currentIndex; ++n){
      long maxOffsetCount = hitArray[n]->maxOffsetCount();
      map<long,long>::iterator it = hitCountsToCount.find( maxOffsetCount );
      if (it == hitCountsToCount.end()){ hitCountsToCount.insert( pair<long,long>(maxOffsetCount,1) ); }
      else { it->second++; }
    }

    // figure out the max hit count - it is the last element in the map (since the map keys are ordered)
    // this is safe because I already ensured that the map is not empty   
    // note that the last hit count could have multiple entries, so that count must be added (and its min
    // is 1, so it also guarantees that this number is >= 1)
    map<long,long>::reverse_iterator lastEntryIt = hitCountsToCount.rbegin();
    long maxHitCount = lastEntryIt->first + lastEntryIt->second;

    // now, figure out the index of the sorted list where each seed count #'s sublist will begin
    // decrement since the indexes (num counts) of the map will be retrieved in increasing order
    long* hitCountsToStart = new long[maxHitCount];
    long currentStart = currentIndex;
    for (map<long,long>::iterator it = hitCountsToCount.begin(); it != hitCountsToCount.end(); ++it){
      currentStart -= it->second;
      hitCountsToStart[it->first] = currentStart;
    }

    // now, create the new sorted list and fill it
    for (long n = 0; n < currentIndex; ++n){
      sortedHitWithOffsets[hitCountsToStart[hitArray[n]->maxOffsetCount()]] = hitArray[n];
      hitCountsToStart[hitArray[n]->maxOffsetCount()]++;
    }
    delete [] hitCountsToStart;
  }
			     
  delete [] queryNumArray;
  delete [] hitArray;
  return sortedHitWithOffsets;
}


void ScoredSeqCollectionBwt::getSubtextWordQueries(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet, char sense, long minOverlap){
  // local so that filterN's can be done in here but not repeated
  vector<WordQuery*> localQueries;

  // this is the whole-sequence case
  long blockStartPos = 0;
  long blockEndPos = querySeq->size();
  getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, false);

  // these are the edge cases
  if (_edgeScalingEnabled){
    long fragmentLen = minOverlap;
    long qSize = querySeq->size();
    while (fragmentLen < qSize){
      // the front end
      blockStartPos = 0;
      blockEndPos = fragmentLen;
      getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);
      // the back end
      blockStartPos = qSize - fragmentLen;
      blockEndPos = qSize;
      getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);
      // scale up the window
      fragmentLen *= 2;
    }
  }

  subtextWordQueriesFilterNs(querySeq, seqSeq, &localQueries);
  wordQuerySet->insert(wordQuerySet->end(), localQueries.begin(), localQueries.end());
}


void ScoredSeqCollectionBwt::getSubtextWordQueriesFullAlignment(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet,
								char sense, long minTargetLength){
  // local so that filterN's can be done in here but not repeated
  vector<WordQuery*> localQueries;

  // this is the whole-sequence case; it will allow hits to the entire query to be aligned
  long blockStartPos = 0;
  long blockEndPos = querySeq->size();
  getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);

  long fragmentLen = minTargetLength / 2;
  while (fragmentLen < querySeq->size()){
    getSubtextWordQueriesHelper(querySeq, &localQueries, sense, fragmentLen);
    // scale up the window
    fragmentLen *= 2;
  }
  subtextWordQueriesFilterNs(querySeq, seqSeq, &localQueries);
  wordQuerySet->insert(wordQuerySet->end(), localQueries.begin(), localQueries.end());
}


void ScoredSeqCollectionBwt::subtextWordQueriesFilterNs(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet){
  if (wordQuerySet->size() > 0){
    WordQueryBinOrganizer* wqbo = new WordQueryBinOrganizer(wordQuerySet);

    // now do the filtering
    long numBlocks = wqbo->numBlocks();
    for (long block = 0; block < numBlocks; ++block){
      if (wqbo->blockHasQueries(block)){
	// scan for N
	long pos = wqbo->blockStart(block);
	long blockEnd = wqbo->blockEnd(block);
	while (pos < blockEnd and seqSeq[pos] != 'N'){ ++pos; }
	if (pos != blockEnd){
	  // an N was found, so all overlapping windows can be gotten rid of!
	  long* queriesInBin = wqbo->getBlockQueryIndexes(block);
	  long qibCount = wqbo->getBlockQueryIndexCount(block);
	  for (long qN = 0; qN < qibCount; ++qN){
	    wqbo->removeQuery( queriesInBin[qN], block );
	    delete wqbo->getQuery( queriesInBin[qN] );
	  }
	  delete [] queriesInBin;
	}
      }
    }
    wordQuerySet->clear();
    wqbo->getQueries(wordQuerySet);
    delete wqbo;
  }
}



ScoredSeqCollectionBwt::WordQueryBinOrganizer::WordQueryBinOrganizer(vector<WordQuery*>* wordQueries){
  _numQueries = wordQueries->size();
  long numQp1 = _numQueries + 1; // just an optimization; it is used more later in this constructor
  _queries = new WordQuery*[ numQp1 ];
  _queryIsValid = new bool[ numQp1 ];
  vector<WordQuery*>::iterator inputIt = wordQueries->begin();
  for (long n = 0; n < _numQueries; ++n){
    _queryIsValid[n] = true;
    _queries[n] = *inputIt;
    ++inputIt;
  }

  // get all the non-redundant edges sorted
  set<long> edgeSetStarts;
  set<long> edgeSetEnds;
  long maxEdge = 0;
  for (long n = 0; n < _numQueries; ++n){
    edgeSetStarts.insert( _queries[n]->start() );
    long end = _queries[n]->end();
    if (end > maxEdge){ maxEdge = end; }
    edgeSetEnds.insert( end );
  }

  // only the elements that are actually starts will be occupied
  _startToBlock = new long[ maxEdge+1 ];

  // this can save time because there should be a lot of redundancy between
  // starts, and also between ends, but not combining starts and ends
  set<long> edgeSet;
  edgeSet.insert(edgeSetStarts.begin(), edgeSetStarts.end());
  edgeSet.insert(edgeSetEnds.begin(), edgeSetEnds.end());

  // use the edges to construct blocks; they will be ordered in the set
  long edgeSetSizeP1 = edgeSet.size() + 1;
  _blockToStart = new long[ edgeSetSizeP1 ];
  _blockToEnd = new long[ edgeSetSizeP1 ];
  _blockToBinCount = new long[ edgeSetSizeP1 ];
  _blockToQisInBin = new long*[ edgeSetSizeP1 ];
  _blockToQiBinIndex = new long*[ edgeSetSizeP1 ];
  _numBlocks = edgeSet.size() - 1;
  if (_numBlocks < 0){ _numBlocks = 0; }

  long block = 0;
  for (set<long>::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it){
    _blockToBinCount[block] = 0;
    _blockToStart[block] = *it;
    _blockToQisInBin[block] = new long[ numQp1 ];
    _blockToQiBinIndex[block] = new long[ numQp1 ];
    if (block > 0){ _blockToEnd[block-1] = *it; }
    _startToBlock[*it] = block;
    ++block;
  }

  // insert each word query into all of its blocks
  for (long n = 0; n < _numQueries; ++n){
    long nb = _startToBlock[ _queries[n]->start() ];
    long endEdge = _startToBlock[ _queries[n]->end() ];
    while (nb < endEdge){
      _blockToQisInBin[nb][ _blockToBinCount[nb] ] = n;
      _blockToQiBinIndex[nb][n] = _blockToBinCount[nb];
      _blockToBinCount[nb]++;
      ++nb;
    }
  }
}
ScoredSeqCollectionBwt::WordQueryBinOrganizer::~WordQueryBinOrganizer(){
  delete [] _blockToStart;
  delete [] _blockToEnd;
  delete [] _startToBlock;

  delete [] _blockToBinCount;
  for (long n = 0; n <= _numBlocks; ++n){
    delete [] _blockToQisInBin[n];
    delete [] _blockToQiBinIndex[n];
  }
  delete [] _blockToQisInBin;
  delete [] _blockToQiBinIndex;

  delete [] _queries;
  delete [] _queryIsValid;
}

long ScoredSeqCollectionBwt::WordQueryBinOrganizer::numQueries(){ return _numQueries; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::numBlocks(){ return _numBlocks; }
bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockHasQueries(long block){
  return _blockToBinCount[block] > 0;
}
bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockHasQuery(long block, long qIndex){
  return _queryIsValid[qIndex] and _blockToStart[block] >= _queries[qIndex]->start() and _blockToEnd[block] <= _queries[qIndex]->end();
}
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockStart(long block){ return _blockToStart[block]; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockEnd(long block){ return _blockToEnd[block]; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::getBlockQueryIndexCount(long block){
  return _blockToBinCount[block];
}
long* ScoredSeqCollectionBwt::WordQueryBinOrganizer::getBlockQueryIndexes(long block){
  long* indexes = new long[ _blockToBinCount[block] + 1 ];
  for (long n = 0; n < _blockToBinCount[block]; ++n){ indexes[n] = _blockToQisInBin[block][n]; }
  return indexes;
}


bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::queryStillValid(long qIndex){
  return _queryIsValid[qIndex];
}
ScoredSeqCollectionBwt::WordQuery* ScoredSeqCollectionBwt::WordQueryBinOrganizer::getQuery(long qIndex){
  return _queries[qIndex];
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeQuery(long qIndex){
  long nb = _startToBlock[ _queries[qIndex]->start() ];
  long endEdge = _startToBlock[ _queries[qIndex]->end() ];
  while (nb < endEdge){
    _blockToBinCount[nb]--;
    if (_blockToBinCount[nb] > 0){
      // this line moves the last element in the array of qIndexes to the position of the one
      // being removed, so that all active qIndexes are now found at the front of the array
      long movedQi = _blockToQisInBin[nb][ _blockToBinCount[nb] ];
      long oldIndex = _blockToBinCount[nb];
      long newIndex = _blockToQiBinIndex[nb][qIndex];
      _blockToQisInBin[nb][newIndex] = movedQi;
      // since the qIndex above was moved, the array that helps find it must also be updated
      _blockToQiBinIndex[nb][movedQi] = newIndex;
    }
    ++nb;
  }
  _queryIsValid[ qIndex ] = false;
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeQuery(long qIndex, long currentBlock){
  _queryIsValid[ qIndex ] = false;
  long nb = _startToBlock[ _queries[qIndex]->start() ];
  if (nb <= currentBlock){ nb = currentBlock + 1; }
  long endEdge = _startToBlock[ _queries[qIndex]->end() ];
  while (nb < endEdge){
    _blockToBinCount[nb]--;
    if (_blockToBinCount[nb] > 0){
      long movedQi = _blockToQisInBin[nb][ _blockToBinCount[nb] ];
      long oldIndex = _blockToBinCount[nb];
      long newIndex = _blockToQiBinIndex[nb][qIndex];
      _blockToQisInBin[nb][newIndex] = movedQi;
      // since the qIndex above was moved, the array that helps find it must also be updated
      _blockToQiBinIndex[nb][movedQi] = newIndex;
    }
    ++nb;
  }
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeSupersetQueries(long qIndex){
  long firstBin = _startToBlock[ _queries[qIndex]->start() ];
  long lastBin = _startToBlock[ _queries[qIndex]->end() ] - 1;
  long* inFirstBin = getBlockQueryIndexes(firstBin);
  long qiCount = getBlockQueryIndexCount(firstBin);
  if (firstBin==lastBin){
    for (long n = 0; n < qiCount; ++n){ removeQuery( inFirstBin[n] ); }
  } else {
    for (long n = 0; n < qiCount; ++n){
      if ( blockHasQuery(lastBin, inFirstBin[n]) ){ removeQuery( inFirstBin[n] ); }
    }
  }
  delete [] inFirstBin;
}

void ScoredSeqCollectionBwt::WordQueryBinOrganizer::getQueries(vector<WordQuery*>* wqSet){
  for (long n = 0; n < _numQueries; ++n){
    if (_queryIsValid[n]){ wqSet->push_back( _queries[n] ); }
  }
}



  // num allowed errors:        0  1  2  3  4  5  6  7  8  9  10
int ScoredSeqCollectionBwt::_specialCaseLenDenom[6] =  {1, 2, 3, 3, 4, 5};//, 6, 7, 8, 9, 10};
int ScoredSeqCollectionBwt::_specialCaseStepDenom[6] = {1, 1, 1, 2, 2, 1};//, 1, 1, 1, 1, 1};
int ScoredSeqCollectionBwt::_specialCaseWindowNum[6] = {1, 2, 3, 5, 7, 5};//, 6, 7, 8, 9, 10};

inline long ScoredSeqCollectionBwt::getNumWindowsHelper(long numAllowedMis){
  if (numAllowedMis < 6){ return _specialCaseWindowNum[numAllowedMis]; }
  else {
    // see notebook7 10/26/10 for derivation
    switch(numAllowedMis){
    case 6: return 6;
    case 7: return 6;
    case 8: return 7;
    case 9: return 7;
    case 10: return 8;
    case 11: return 8;
    case 12: return 8;
    case 13: return 9;
    case 14: return 9;
    case 15: return 10;
    case 16: return 10;
    case 17: return 11;
    case 18: return 11;
    case 19: return 11;
    case 20: return 12;
    case 21: return 12;
    case 22: return 13;
    case 23: return 13;
    case 24: return 13;
    case 25: return 14;
    case 26: return 14;
    case 27: return 14;
    case 28: return 15;
    case 29: return 15;
    case 30: return 15;
    case 31: return 16;
    case 32: return 16;
    case 33: return 16;
    case 34: return 17;
    case 35: return 17;
    case 36: return 17;
    case 37: return 18;
    case 38: return 18;
    case 39: return 18;
    case 40: return 19;
    case 41: return 19;
    case 42: return 19;
    case 43: return 20;
    case 44: return 20;
    case 45: return 20;
    case 46: return 21;
    case 47: return 21;
    case 48: return 21;
    case 49: return 21;
    case 50: return 22;
    case 51: return 22;
    case 52: return 22;
    case 53: return 23;
    case 54: return 23;
    case 55: return 23;
    case 56: return 24;
    case 57: return 24;
    case 58: return 24;
    case 59: return 24;
    case 60: return 25;
    case 61: return 25;
    case 62: return 25;
    case 63: return 26;
    case 64: return 26;
    case 65: return 26;
    case 66: return 27;
    case 67: return 27;
    case 68: return 27;
    case 69: return 27;
    case 70: return 28;
    case 71: return 28;
    case 72: return 28;
    case 73: return 29;
    case 74: return 29;
    case 75: return 29;
    case 76: return 29;
    case 77: return 30;
    case 78: return 30;
    case 79: return 30;
    case 80: return 31;
    case 81: return 31;
    case 82: return 31;
    case 83: return 31;
    case 84: return 32;
    case 85: return 32;
    case 86: return 32;
    case 87: return 32;
    case 88: return 33;
    case 89: return 33;
    case 90: return 33;
    case 91: return 34;
    case 92: return 34;
    case 93: return 34;
    case 94: return 34;
    case 95: return 35;
    case 96: return 35;
    case 97: return 35;
    case 98: return 35;
    case 99: return 36;
    case 100: return 36;
    case 101: return 36;
    case 102: return 36;
    case 103: return 37;
    case 104: return 37;
    case 105: return 37;
    case 106: return 38;
    case 107: return 38;
    case 108: return 38;
    case 109: return 38;
    case 110: return 39;
    case 111: return 39;
    case 112: return 39;
    case 113: return 39;
    case 114: return 40;
    case 115: return 40;
    case 116: return 40;
    case 117: return 40;
    case 118: return 41;
    case 119: return 41;
    case 120: return 41;
    case 121: return 41;
    case 122: return 42;
    case 123: return 42;
    case 124: return 42;
    case 125: return 42;
    case 126: return 43;
    case 127: return 43;
    case 128: return 43;
    case 129: return 43;
    case 130: return 44;
    case 131: return 44;
    case 132: return 44;
    case 133: return 45;
    case 134: return 45;
    case 135: return 45;
    case 136: return 45;
    case 137: return 46;
    case 138: return 46;
    case 139: return 46;
    case 140: return 46;
    default: return numAllowedMis / 3;
    }
  }
}


void ScoredSeqCollectionBwt::getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
							 long blockStartPos, long blockEndPos, bool useLimits){
  long seqLen = blockEndPos - blockStartPos;
  bool fiveEdge;
  if (blockStartPos == 0){ fiveEdge = true; }
  else if (blockEndPos == querySeq->size()){ fiveEdge = false; }
  else { throw AssemblyException::ArgError("SSCBwt::getSWQHelper block is at neither the front nor back of the query."); }

  long numAllowedMis = long(float(seqLen) * _maxFractMis);
  // a safety precaution for a rediculous situation
  if (numAllowedMis==seqLen){ numAllowedMis = seqLen - 1; }

  // figure these out for the special or general cases
  long startUpperLimit;
  long startFactor;
  long subLen;
  long fragLenT2;

  // special cases
  if (numAllowedMis < 6){
    long startDenom = _specialCaseLenDenom[numAllowedMis] * _specialCaseStepDenom[numAllowedMis];
    subLen = seqLen / _specialCaseLenDenom[numAllowedMis];
    long numWindows = getNumWindowsHelper(numAllowedMis);
    startFactor = seqLen / startDenom;
    startUpperLimit = blockStartPos + (numWindows * startFactor);
  } else { // general cases
    long numWindows = getNumWindowsHelper(numAllowedMis);
    subLen = seqLen / numWindows;
    if (subLen < 2){ subLen = 2; }
    startFactor =  seqLen / numWindows;
    // i am putting the booleans outside the loop to save time
    // the input args below with the equation calculates the start position
    startUpperLimit = blockStartPos + (numWindows * startFactor);
  }

  // now create the word queries using the parameters calculated above
  if (! useLimits){
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryNoEdgeLimit(querySeq, start, subLen, sense));
    }
  } else if (fiveEdge){
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryFiveEdgeLimit(querySeq, start, subLen, sense, seqLen));
    }
  } else {
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryThreeEdgeLimit(querySeq, start, subLen, sense, seqLen));
    }
  }
}



void ScoredSeqCollectionBwt::getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
							 long fragmentLength){

  long numAllowedMis = long(float(fragmentLength) * _maxFractMis);
  // a safety precaution for a rediculous situation
  if (numAllowedMis==fragmentLength){ numAllowedMis = fragmentLength - 1; }
  long querySize = querySeq->size();

  // figure these out for the special or general cases
  long startUpperLimit;
  long startFactor;
  long subLen;
  long fragLenT2;

  // special cases
  if (numAllowedMis < 6){
    long startDenom = _specialCaseLenDenom[numAllowedMis] * _specialCaseStepDenom[numAllowedMis];
    subLen = fragmentLength / _specialCaseLenDenom[numAllowedMis];
    long numWindows = getNumWindowsHelper(numAllowedMis) * querySize / fragmentLength;
    startFactor = fragmentLength / startDenom;
    fragLenT2 = fragmentLength * 2;
    // increment up using the startFactor
    startUpperLimit = numWindows * startFactor;
    // the highest word cannot read off the edge of the query sequence
    if (startUpperLimit + subLen > querySize + 1){ startUpperLimit = querySize + 1 - subLen; }
  } else { // general cases
    long numWindows = getNumWindowsHelper(numAllowedMis);
    // adjust for the query length
    subLen = fragmentLength / numWindows;
    if (subLen < 2){ subLen = 2; }
    numWindows = numWindows * querySize / fragmentLength;
    // other constants that only need be calculated once
    fragLenT2 = fragmentLength * 2;
    startFactor = fragmentLength / numWindows;
    // increment up using the startFactor
    startUpperLimit = numWindows * startFactor;
    // the highest word cannot read off the edge of the query sequence
    if (startUpperLimit + subLen > querySize + 1){ startUpperLimit = querySize + 1 - subLen; }
  }

  // now create the windows using the variables calculated above
  for (long start = 0; start < startUpperLimit; start += startFactor){
    wordQuerySet->push_back(new WordQueryFullTarget(querySeq, start, subLen, sense, _fractId, fragLenT2));
  }

}





void ScoredSeqCollectionBwt::bufferSeqs(){
  //OK();
  for ( set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it ){ (*it)->buffer(); }
  //OK();
}

void ScoredSeqCollectionBwt::unbufferSeqs(){
  //OK();
  for ( set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it ){ (*it)->unbuffer(); }
  //OK();
}



long ScoredSeqCollectionBwt::size(){
  //OK();
  return _allSeqs.size();
}


// overlaps is the return value, locus is the query
// ASSUMES that the locus only overlaps one contig
long ScoredSeqCollectionBwt::getOverlapContigIndex(long start, long end){

  // determine the bin to look in
  long startBin = start / _binSize;
  if (startBin < 0){ startBin = 0; }

  // return the legit match's index
  vector<ScoredSeqLocus*>::iterator firstLocusIt = _binToSeqs[startBin].begin();
  while (firstLocusIt != _binToSeqs[startBin].end() and (! (*firstLocusIt)->overlaps(start,end))){ ++firstLocusIt; }
  if ( firstLocusIt == _binToSeqs[startBin].end() ){ throw AssemblyException::LogicError("SSCBwt::gOCI can't find match."); }
  return (*firstLocusIt)->_index;
}


ScoredSeqCollectionBwt::ScoredSeqLocus::ScoredSeqLocus(){}
ScoredSeqCollectionBwt::ScoredSeqLocus::ScoredSeqLocus(ScoredSeq* seq, long start, long end, long index) :
  _start(start),
  _end(end),
  _seq(seq),
  _index(index) {
  }
bool ScoredSeqCollectionBwt::ScoredSeqLocus::overlaps(long otherStart, long otherEnd){
  return (! (_start > otherEnd or otherStart > _end) );
}
bool ScoredSeqCollectionBwt::ScoredSeqLocus::contains(long otherStart, long otherEnd){
  return (! (_start > otherStart or otherEnd > _end) );
}



ScoredSeqCollectionBwt::AlignmentMaker::~AlignmentMaker(){}

bool ScoredSeqCollectionBwt::AlignmentMaker::makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart,
								     long maxNumMis, long exBlockCount,
								     long* examineStarts, long* examineEnds){


  // gaps are irrelevant
  AlignmentScoreMatrix* misCountAsm = new AlignmentScoreMatrix(0,1,1,1);
  long misCount = 0;

  while ( (maxNumMis >= misCount) and (exBlockCount > 0) ){
    --exBlockCount;
    long posA = examineEnds[exBlockCount];
    long posB = bStart + (posA - aStart);
    while ( (maxNumMis >= misCount) and (posA > examineStarts[exBlockCount]) ){
      --posA;
      --posB;
      misCount += misCountAsm->match(seqA->nucAtPosPlus(posA), seqB->nucAtPosPlus(posB));
    }
  }

  delete misCountAsm;
  return maxNumMis >= misCount;
}
bool ScoredSeqCollectionBwt::AlignmentMaker::makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart,
								      long maxNumMis, long exBlockCount,
								      long* examineStarts, long* examineEnds){


  // gaps are irrelevant
  AlignmentScoreMatrix* misCountAsm = new AlignmentScoreMatrix(0,1,1,1);
  long misCount = 0;

  while ( (maxNumMis >= misCount) and (exBlockCount > 0) ){
    --exBlockCount;
    long posA = examineEnds[exBlockCount];
    long posB = bStart + (posA - aStart);
    while ( (maxNumMis >= misCount) and (posA > examineStarts[exBlockCount]) ){
      --posA;
      --posB;
      misCount += misCountAsm->match(seqA->nucAtPosPlus(posA), seqB->nucAtPosMinus(posB));
    }
  }

  delete misCountAsm;
  return maxNumMis >= misCount;
}




ScoredSeqCollectionBwt::AlMakerNoEdge::AlMakerNoEdge(){}
ScoredSeqCollectionBwt::AlMakerNoEdge::AlMakerNoEdge(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis) :
  _scoreMatrix(scoreMatrix),
  _minScore(minScore),
  _maxFractMis(maxFractMis)
{
  if (_scoreMatrix->_match <= _scoreMatrix->_mismatch){
    throw AssemblyException::ArgError("SSCBwt::AlMakerNoEdge needs a lower mismatch score than match score (i.e. mismatch should be negative)");
  }
}
ScoredSeqCollectionBwt::AlMakerNoEdge::~AlMakerNoEdge(){}
Alignment* ScoredSeqCollectionBwt::AlMakerNoEdge::makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense,
								OffsetAndQueryCarrier* offsetCarrier){
  long offset = offsetCarrier->_offset;

  // get the starts and alignment length
  long aStart;
  long bStart;
  if (offset < 0){
    aStart = 0;
    bStart = 0 - offset;
  } else {
    aStart = offset;
    bStart = 0;
  }
  long aEnd;
  if (offset + seqB->size() >= seqA->size()){ aEnd = seqA->size(); }
  else { aEnd = offset + seqB->size(); }
  long alLength = aEnd - aStart;

  // i can maybe take out these tests
  if (aStart != 0 and bStart != 0){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither start is zero");
  }
  if ( aStart + alLength != seqA->size() and bStart + alLength != seqB->size() ){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither end is the end of the sequence");
  }

  // the min score must be beaten, and the fractId must also be satisfied
  long maxNumMisByScore = (_minScore - _scoreMatrix->_match * alLength) / (_scoreMatrix->_mismatch - _scoreMatrix->_match );
  long maxNumMisByFid = long(float(alLength) * _maxFractMis);
  long maxNumMis;
  if (maxNumMisByScore < maxNumMisByFid){ maxNumMis = maxNumMisByScore; }
  else { maxNumMis = maxNumMisByFid; }

  // lowest to highest
  long exBlockCount = 1 + offsetCarrier->_numQueries;
  long* examineStarts = new long[exBlockCount];
  long* examineEnds = new long[exBlockCount];
  examineStarts[0] = aStart;
  examineEnds[exBlockCount-1] = aStart + alLength;
  for (long n = 0; n < offsetCarrier->_numQueries; ++n){
    examineEnds[n] = offsetCarrier->_queries[n]->start();
    examineStarts[n+1] = offsetCarrier->_queries[n]->end();
  }

  bool alIsGood;
  switch (sense){
  case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
  case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
  default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
  }

  delete [] examineStarts;
  delete [] examineEnds;

  if (alIsGood){ return new AlignmentUngapped(seqA, seqB, sense, offset); }
  else { return new AlignmentNull(seqA, seqB, sense); }
}

long ScoredSeqCollectionBwt::AlMakerNoEdge::getScore(Alignment* al){
  return al->score(_scoreMatrix, false);
}
long ScoredSeqCollectionBwt::AlMakerNoEdge::getMinScore(){ return _minScore; }
void ScoredSeqCollectionBwt::AlMakerNoEdge::setMinScore(long minScore){ _minScore = minScore; }
AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerNoEdge::getScoreMatrix(){ return _scoreMatrix; }



ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::AlMakerPenalizeEdgeA(){}
ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::AlMakerPenalizeEdgeA(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis) :
  _scoreMatrix(scoreMatrix),
  _minScore(minScore),
  _maxFractMis(maxFractMis)
{
  if (_scoreMatrix->_match <= _scoreMatrix->_mismatch){
    throw AssemblyException::ArgError("SSCBwt::AlMakerNoEdge needs a lower mismatch score than match score (i.e. mismatch should be negative)");
  }
}
ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::~AlMakerPenalizeEdgeA(){}

Alignment* ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense,
								       OffsetAndQueryCarrier* offsetCarrier){
  long offset = offsetCarrier->_offset;

  // get the starts and alignment length
  long aStart;
  long bStart;
  if (offset < 0){
    aStart = 0;
    bStart = 0 - offset;
  } else {
    aStart = offset;
    bStart = 0;
  }
  long aEnd;
  if (offset + seqB->size() >= seqA->size()){ aEnd = seqA->size(); }
  else { aEnd = offset + seqB->size(); }
  long alLength = aEnd - aStart;

  // i can maybe take out these tests
  if (aStart != 0 and bStart != 0){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither start is zero");
  }
  if ( aStart + alLength != seqA->size() and bStart + alLength != seqB->size() ){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither end is the end of the sequence");
  }

  // the min score must be beaten, and the fractId must also be satisfied
  long maxNumMisByScore = (_minScore - _scoreMatrix->_match * seqA->size()) / (_scoreMatrix->_mismatch - _scoreMatrix->_match );
  long maxNumMisByFid = long(float(seqA->size()) * _maxFractMis);
  long maxNumMis;
  if (maxNumMisByScore < maxNumMisByFid){ maxNumMis = maxNumMisByScore; }
  else { maxNumMis = maxNumMisByFid; }
  // adjust for the edges
  maxNumMis -= seqA->size() - alLength;

  // lowest to highest
  long exBlockCount = 1 + offsetCarrier->_numQueries;
  long* examineStarts = new long[exBlockCount];
  long* examineEnds = new long[exBlockCount];
  examineStarts[0] = aStart;
  examineEnds[exBlockCount-1] = aStart + alLength;
  for (long n = 0; n < offsetCarrier->_numQueries; ++n){
    examineEnds[n] = offsetCarrier->_queries[n]->start();
    examineStarts[n+1] = offsetCarrier->_queries[n]->end();
  }

  bool alIsGood;
  switch (sense){
  case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
  case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
  default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
  }

  delete [] examineStarts;
  delete [] examineEnds;

  if (alIsGood){ return new AlignmentUngapped(seqA, seqB, sense, offset); }
  else { return new AlignmentNull(seqA, seqB, sense); }
}
long ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getScore(Alignment* al){
  return al->scoreOverhangA(_scoreMatrix);
}
long ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getMinScore(){ return _minScore; }
void ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::setMinScore(long minScore){ _minScore = minScore; }
AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getScoreMatrix(){ return _scoreMatrix; }




ScoredSeqCollectionBwt::OffsetTracker::OffsetTracker(ScoredSeq* target) : _targetContig(target){
  _maxOffsetCount = 0;
  // this value is meaningless unless the above values are >0
  _maxCountOffset = 0;
}
ScoredSeqCollectionBwt::OffsetTracker::~OffsetTracker(){
  for (map<long,OffsetFields*>::iterator it = _offsetToFields.begin(); it != _offsetToFields.end(); ++it){
    delete it->second;
  }
}

ScoredSeqCollectionBwt::OffsetTracker::OffsetFields::OffsetFields(float count, WordQuery* maxWord) : 
  _count(count),
  _maxWord(new WordQuerySimple(maxWord)){}
ScoredSeqCollectionBwt::OffsetTracker::OffsetFields::~OffsetFields(){
  delete _maxWord;
  for (vector<WordQuerySimple*>::iterator wqIt = _otherWords.begin(); wqIt != _otherWords.end(); ++wqIt){ delete *wqIt; }
}
void ScoredSeqCollectionBwt::OffsetTracker::OffsetFields::addOffset(float localAddition, WordQuery* seedingQuery){
  _count += localAddition;
  if ( seedingQuery->notContinuous(_maxWord) ){
    // figure out which word is more important (longer)
    if (seedingQuery->size() > _maxWord->size()){
      _otherWords.push_back(_maxWord);
      _maxWord = new WordQuerySimple(seedingQuery);
    } else { _otherWords.push_back(new WordQuerySimple(seedingQuery)); }
  } else {
    _maxWord->merge(seedingQuery);
  }
}

void ScoredSeqCollectionBwt::OffsetTracker::addOffset(long offset, WordQuery* seedingQuery, long normFactor){
  float localCount;
  if (normFactor < 1){ throw AssemblyException::ArgError("SSCBwt::OffsetTracker::addOffset, normFactor cannot be less than 1"); }
  float localAddition = float(1) / float(normFactor);
  // now find or create the offset entry
  map<long,OffsetFields*>::iterator offsetEntry = _offsetToFields.find(offset);
  if (offsetEntry == _offsetToFields.end()){
    localCount = localAddition;
    _offsetToFields.insert( pair<long,OffsetFields*>(offset, new OffsetFields(localAddition, seedingQuery)) );
  } else {
    offsetEntry->second->addOffset(localAddition, seedingQuery);
    localCount = offsetEntry->second->_count;
  }
  if (localCount > _maxOffsetCount){
    _maxOffsetCount = localCount;
    _maxCountOffset = offset;
  }
}
long ScoredSeqCollectionBwt::OffsetTracker::maxOffsetCount(){ return _maxOffsetCount; }
long ScoredSeqCollectionBwt::OffsetTracker::maxCountOffset(){ return _maxCountOffset; }
ScoredSeqCollectionBwt::OffsetAndQueryCarrier** ScoredSeqCollectionBwt::OffsetTracker::getOffsets(bool includeQueries){
  // figure out how many offsets there are with each count
  map<float,long> countToNum;
  for (map<long,OffsetFields*>::iterator it = _offsetToFields.begin(); it != _offsetToFields.end(); ++it){
    map<float,long>::iterator foundCtm = countToNum.find( it->second->_count );
    if (foundCtm == countToNum.end()){ countToNum.insert(pair<float,long>(it->second->_count,1)); }
    else { foundCtm->second += 1; }
  }
  // figure out the index in the sorted array at which each offset count will begin
  map<float,long> countToIndex;
  map<float,long>::iterator ctiIt = countToIndex.begin();
  long index = _offsetToFields.size();
  for (map<float,long>::iterator it = countToNum.begin(); it != countToNum.end(); ++it){
    index -= it->second;
    // i can skip to the end because the keys were already sorted in the other map
    ctiIt = countToIndex.insert(ctiIt, pair<float,long>(it->first,index));
  }

  long numOutput = _offsetToFields.size();
  OffsetAndQueryCarrier** output = new OffsetAndQueryCarrier*[ numOutput + 1];
  long outIndex = 0;
  map<long,OffsetFields*>::iterator fieldIt = _offsetToFields.begin();
  while (outIndex < numOutput){
    // this is a test - i can comment it out later
    //if (countIt->first != queryIt->first){ throw AssemblyException::LogicError("SSCBwt:OT:getOffsets maps are messed up"); }

    OffsetAndQueryCarrier* newOut;
    if (includeQueries){
      // determine if there are any other queries to deal with
      if (fieldIt->second->_otherWords.size() == 0){
	newOut = new OffsetAndQueryCarrier(fieldIt->first, 1);
	newOut->_queries[0] = new WordQuerySimple(fieldIt->second->_maxWord);
      } else {
	// sort the queries and make them non-redundant while creating the output
	vector<WordQuery*> sortQueries;
	sortQueries.insert(sortQueries.end(), fieldIt->second->_otherWords.begin(), fieldIt->second->_otherWords.end());
	sortQueries.push_back( fieldIt->second->_maxWord );
	SortWqByStart wqSorter;
	sort( sortQueries.begin(), sortQueries.end(), wqSorter );

	WordQuerySimple** wqArray = new WordQuerySimple*[ sortQueries.size() ];
	vector<WordQuery*>::iterator sortIt = sortQueries.begin();
	long wqIndex = 0;
	wqArray[wqIndex] = new WordQuerySimple( *sortIt );
	++sortIt;
	while (sortIt != sortQueries.end()){
	  if (wqArray[wqIndex]->notContinuous( *sortIt )){
	    ++wqIndex;
	    wqArray[wqIndex] = new WordQuerySimple( *sortIt );
	  } else {
	    wqArray[wqIndex]->merge( *sortIt );
	  }
	  ++sortIt;
	}
	// advance wqIndex one beyond the end of the array's last active element
	++wqIndex;
	newOut = new OffsetAndQueryCarrier(fieldIt->first, wqIndex);
	for (long n = 0; n < wqIndex; ++n){ newOut->_queries[n] = wqArray[n]; }
	delete [] wqArray;
      }
    } else {
      newOut = new OffsetAndQueryCarrier(fieldIt->first, 0);
    }
    map<float,long>::iterator foundCtm = countToIndex.find( fieldIt->second->_count );
    output[ foundCtm->second ] = newOut;
    foundCtm->second += 1;
    ++outIndex;
    ++fieldIt;
  }
  return output;
}
bool ScoredSeqCollectionBwt::OffsetTracker::SortWqByStart::operator() (WordQuery* wqA, WordQuery* wqB){
  return wqA->start() < wqB->start();
}


ScoredSeqCollectionBwt::OffsetAndQueryCarrier::OffsetAndQueryCarrier(long offset, long numQueries) :
  _offset(offset), _numQueries(numQueries)
{
  _queries = new WordQuery*[ numQueries+1 ];
}
ScoredSeqCollectionBwt::OffsetAndQueryCarrier::~OffsetAndQueryCarrier(){
  delete [] _queries;
}



long ScoredSeqCollectionBwt::OffsetTracker::numOffsets(){ return long(_offsetToFields.size()); }

ScoredSeq* ScoredSeqCollectionBwt::OffsetTracker::targetContig(){ return _targetContig; }


ScoredSeqCollectionBwt::WordQuery::~WordQuery(){}


ScoredSeqCollectionBwt::WordQuerySimple::WordQuerySimple(long start, long size) :
  _start(start), _size(size), _end(start+size)
{}
ScoredSeqCollectionBwt::WordQuerySimple::WordQuerySimple(WordQuery* wq){
  _start = wq->start();
  _size = wq->size();
  _end = wq->end();
}
void ScoredSeqCollectionBwt::WordQuerySimple::merge(WordQuery* wq){
  if (wq->start() < _start){ _start = wq->start(); }
  if (wq->end() > _end){ _end = wq->end(); }
  _size = _end - _start;

}

ScoredSeqCollectionBwt::WordQuerySimple::~WordQuerySimple(){}
long ScoredSeqCollectionBwt::WordQuerySimple::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQuerySimple::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQuerySimple::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQuerySimple::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQuerySimple::acceptableOverlap(long targetStart, ScoredSeq* target){ return true; }


ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::WordQueryNoEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense) :
  _sourceSeq(sourceSeq), _start(start), _size(size), _sense(sense), _end(start+size)
{}
ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::~WordQueryNoEdgeLimit(){}
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){ return true; }


ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::WordQueryFiveEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense, long edgeLimit) :
  _sourceSeq(sourceSeq), _start(start), _size(size), _sense(sense), 
  _edgeLimit(edgeLimit), _end(start+size) {}
ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::~WordQueryFiveEdgeLimit(){}
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){
  long distance = target->size() - targetStart;
  return distance <= _edgeLimit * 2;
}
ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::WordQueryThreeEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense, long edgeLimit) :
  _sourceSeq(sourceSeq), _start(start), _size(size), _sense(sense), 
  _edgeLimit(edgeLimit), _end(start+size) {}
ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::~WordQueryThreeEdgeLimit(){}
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){
  long distance = targetStart + _size;
  return distance <= _edgeLimit * 2;
}


ScoredSeqCollectionBwt::WordQueryFullTarget::WordQueryFullTarget(ScoredSeq* sourceSeq, long start, long size, char sense, float fractId, long maxTargetSize) :
  _sourceSeq(sourceSeq), _start(start), _end(start+size), _size(size), _sense(sense), _fractId(fractId), _maxTargetSize(maxTargetSize)
{}
ScoredSeqCollectionBwt::WordQueryFullTarget::~WordQueryFullTarget(){}
long ScoredSeqCollectionBwt::WordQueryFullTarget::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryFullTarget::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryFullTarget::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryFullTarget::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryFullTarget::acceptableOverlap(long targetStart, ScoredSeq* target){
  //the target must be smaller than the source sequence, excpting for the fractId compensation
  if (float(_sourceSeq->size()) < float(target->size()) * _fractId){ return false; }
  else { 
    long targetTrim = long((float(1) - _fractId) * float(target->size()) );
    return (targetStart - targetTrim <= _start) and 
      (_start + target->size() - targetStart <= _sourceSeq->size() + targetTrim);
  }
}


#endif


