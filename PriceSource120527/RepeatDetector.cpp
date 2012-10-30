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


#ifndef REPEATDETECTOR_CPP
#define REPEATDETECTOR_CPP


#include "RepeatDetector.h"

#include "ScoredSeq.h"
#include "AssemblyException.h"
#include "AssemblyJob.h"
#include "AssemblyJobFactory.h"
#include "ParamsDeBruijn.h"
#include "AssemblerListenerNull.h"

#include <omp.h>
#include <math.h>
#include <algorithm>

# include <iostream>
using namespace::std;

RepeatDetector::RepeatDetector(long minRepeatSize, ParamsAlignment* alignParams) : _minRepeatSize(minRepeatSize) {
  _windowSize = _minRepeatSize / 2;
  if (_windowSize < 10){ _windowSize = 10; }
  _alignParams = new ParamsAlignment(alignParams);
}
RepeatDetector::~RepeatDetector(){
  delete _alignParams;
}


ScoredSeq** RepeatDetector::findRepeats(set<ScoredSeq*>* contigs, float numStdDevs, float minFoldIncrease,
					RepeatDetectorThreadedness threadedness){
  ScoredSeq** resultsArray;
  if (contigs->size() == 0){
    resultsArray = new ScoredSeq*[ 1 ];
    resultsArray[0] = NULL;
  } else {

    // generate the stats for the individual contigs
    ContigWithStats** cwsArray = new ContigWithStats*[ contigs->size() + 1 ];
    long cCount = 0;
    long wCount = 0;
    // don't thread this loop because the individual constructor calls are threaded
    for (set<ScoredSeq*>::iterator cIt = contigs->begin(); cIt != contigs->end(); ++cIt){
      cwsArray[cCount] = new ContigWithStats(*cIt, _windowSize, threadedness);
      wCount += cwsArray[cCount]->_numWindows;
      ++cCount;
    }

    // now generate the meta-stats
    ContigFragmentWindow** windowArray = new ContigFragmentWindow*[ wCount + 1 ];
    long windowIndex = 0;
    // probably not worth threading
    for (long cN = 0; cN < cCount; ++cN){
      long windowStart = windowIndex;
      windowIndex += cwsArray[cN]->_numWindows;
      long windowEnd = windowIndex;
      for (long wN = windowStart; wN < windowEnd; ++wN){ windowArray[wN] = cwsArray[cN]->getWindow(wN - windowStart); }
    }
    sortWindowArray(windowArray, wCount);

    float totalNormCount = 0;
    for (long n = 0; n < wCount; ++n){ totalNormCount += windowArray[n]->normCount(); }

    // get the median and percentile-derived "std dev"
    float scoreMedian = getWeightedPercentileScore(0.5, windowArray, wCount, totalNormCount);
    float cutoffLimit = scoreMedian * minFoldIncrease;
    // the values are sorted highest-to-lowest, so this will be one std dev ABOVE the median
    // i am using a method here that is valid for the definitions of mean and std dev in the
    // case of a perfect normal distribution.  i am NEVER expecting such a distribution, so
    // these percentile-based metrics will differ substantially from the calculated mean/stdDev.
    // but this measure will better capture what i am really looking for regardless of the value
    // distribution.
    //float scoreStdDev = getWeightedPercentileScore(0.159, windowArray, wCount, totalNormCount) - scoreMedian;
    float scoreStdDev = getWeightedPercentileScore(0.841, windowArray, wCount, totalNormCount) - scoreMedian;
    float scoreCutoff = scoreMedian + (scoreStdDev * numStdDevs);
    if (scoreCutoff < cutoffLimit){ scoreCutoff = cutoffLimit; }

    // sort the windows initially by their status as above or below the coverage cutoff
    long hcWindowCount = 0;
    long bgWindowCount = 0;
    for (long wN = 0; wN < wCount; ++wN){
      // i am using the length-normalized score so that short segments with high scores can still
      // seed the detection of a repeat (if it is too small a segment, the lenght filter will catch it)
      if (windowArray[wN]->normScore() > scoreCutoff){ ++hcWindowCount; }
      //if (windowArray[wN]->score() > scoreCutoff){ ++hcWindowCount; }
      else { ++bgWindowCount; }
    }

    // only sort if there actually are two groups (presumably, bgWindows will never be empty)
    EmissionCarrier* scoreEm;
    EmissionCarrier* linkEm;
    if (hcWindowCount == 0){
      // these are essentially null carriers
      scoreEm = new EmissionCarrier(scoreCutoff, 0.5, 0.5);
      linkEm = new EmissionCarrier(scoreCutoff, 0.5, 0.5);
    } else {
      // collect some stats about the probabilities of scores in the two contexts
      // high-coverage windows; initial values are laplace pseudocounts
      // index 0 => background (lc), index 1 => foreground (hc)
      long hcNucCount[2] = {1, 1};
      long totalNucCount[2] = {2, 2};
      long hcLinkCount[2] = {1, 1};
      long totalLinkCount[2] = {2, 2};

      // THREAD BY WINDOW
      if (threadedness == THREADED){
        #pragma omp parallel for schedule(dynamic)
	for (long wN = 0; wN < wCount; ++wN){
	  long* nucThenLinkThenBg = new long[4];
	  getHcCounts(nucThenLinkThenBg, windowArray[wN], scoreCutoff);
	  // increment the appropriate values protected as atomic
	  int cI = 0;
	  if (windowArray[wN]->score() > scoreCutoff){ cI = 1; }
          #pragma omp atomic
	  hcNucCount[cI] += nucThenLinkThenBg[0];
          #pragma omp atomic
	  hcLinkCount[cI] += nucThenLinkThenBg[1];
          #pragma omp atomic
	  totalNucCount[cI] += nucThenLinkThenBg[2];
          #pragma omp atomic
	  totalLinkCount[cI] += nucThenLinkThenBg[3];
	  delete [] nucThenLinkThenBg;
	}
      } else {
	long* nucThenLinkThenBg = new long[4];
	for (long wN = 0; wN < wCount; ++wN){
	  getHcCounts(nucThenLinkThenBg, windowArray[wN], scoreCutoff);
	  // increment the appropriate values protected as atomic
	  int cI = 0;
	  if (windowArray[wN]->score() > scoreCutoff){ cI = 1; }
	  hcNucCount[cI] += nucThenLinkThenBg[0];
	  hcLinkCount[cI] += nucThenLinkThenBg[1];
	  totalNucCount[cI] += nucThenLinkThenBg[2];
	  totalLinkCount[cI] += nucThenLinkThenBg[3];
	}
	delete [] nucThenLinkThenBg;
      }
      scoreEm = new EmissionCarrier(scoreCutoff,
				    float(hcNucCount[0]) / float(totalNucCount[0]),
				    float(hcNucCount[1]) / float(totalNucCount[1]) );
      linkEm = new EmissionCarrier(scoreCutoff,
				   float(hcLinkCount[0]) / float(totalLinkCount[0]),
				   float(hcLinkCount[1]) / float(totalLinkCount[1]) );
    }

    // get the windows that exceed the above score cutoff and collapse them
    // into a set of contiguous segments
    ScoredSeq*** hcSeqsByContig = new ScoredSeq**[ cCount+1 ];
    if (threadedness == THREADED){
      #pragma omp parallel for schedule(dynamic)
      for (long cN = 0; cN < cCount; ++cN){ hcSeqsByContig[cN] = cwsArray[cN]->getHighCoverageFragments(scoreCutoff, scoreEm, linkEm); }
    } else {
      for (long cN = 0; cN < cCount; ++cN){ hcSeqsByContig[cN] = cwsArray[cN]->getHighCoverageFragments(scoreCutoff, scoreEm, linkEm); }
    }
    delete scoreEm;
    delete linkEm;

    // consolidate into one set (not thread-able)
    set<ScoredSeq*> hcSeqs;
    for (long cN = 0; cN < cCount; ++cN){
      long sN = 0;
      while (hcSeqsByContig[cN][sN] != NULL){
	if( hcSeqsByContig[cN][sN]->size() < _minRepeatSize ){ delete hcSeqsByContig[cN][sN]; }
      else { hcSeqs.insert( hcSeqsByContig[cN][sN] ); }
	++sN;
      }
      delete [] hcSeqsByContig[cN];
    }
    delete [] hcSeqsByContig;

    // collapse redundant repeat elements with a redundancy assembly
    AssemblerListener* nullListener = new AssemblerListenerNull();
    // the last zero makes the job null by specifying the de bruijn will only be applied to sequences
    // with a max length of zero
    ParamsDeBruijn* nullPdb = new ParamsDeBruijn(0, 0, 0);
    AssemblyJobFactory* contigFactory = new AssemblyJobFactory(_alignParams, nullPdb, AssemblyJob::DOUBLESTRANDED, nullListener);
    AssemblyJob* aj = contigFactory->redundancyJob( &hcSeqs );
    if (threadedness==THREADED){ aj->runJob(AssemblyJob::THREADED); }
    else { aj->runJob(AssemblyJob::NOTTHREADED); }
    hcSeqs.clear();
    aj->allAssembledSeqs( &hcSeqs );
    // get the input data that was combined into the output contigs and delete it
    set<ScoredSeq*> discardedInput;
    aj->discardedInputSeqs( &discardedInput );
    for (set<ScoredSeq*>::iterator disIt = discardedInput.begin(); disIt != discardedInput.end(); ++disIt){ (*disIt)->deepDelete(); }
    // clean up the redundancy collapse
    delete aj;
    delete contigFactory;
    delete nullListener;
    delete nullPdb;

    // clean up the overall process
    for (long n = 0; n < cCount; ++n){ delete cwsArray[n]; }
    for (long n = 0; n < wCount; ++n){ delete windowArray[n]; }
    delete [] cwsArray;
    delete [] windowArray;

    // return results array
    resultsArray = new ScoredSeq*[ hcSeqs.size() + 1 ];
    long rIndex = 0;
    for (set<ScoredSeq*>::iterator it = hcSeqs.begin(); it != hcSeqs.end(); ++it){
      resultsArray[rIndex] = *it;
      ++rIndex;
    }
    resultsArray[rIndex] = NULL;
  }

  return resultsArray;
}

// a helper method for the threading in the above method
void RepeatDetector::getHcCounts(long* nucThenLinkThenBg, ContigFragmentWindow* window, float scoreCutoff){
  ScoredSeq* sourceSeq = window->sourceSeq();
  long startPos = window->startPos();
  long nucEndPos = startPos + window->size();
  long linkEndPos = nucEndPos;
  if (linkEndPos == sourceSeq->size()){ --linkEndPos; }

  // do both operations here; nucEndPos is always at least as big as linkEndPos
  long hcNucToAdd = 0;
  long hcLinkToAdd = 0;
  for (long pN = startPos; pN < linkEndPos; ++pN){
    if (sourceSeq->scoreAtPosPlus(pN) > scoreCutoff){ ++hcNucToAdd; }
    if (sourceSeq->linkAfterPosPlus(pN) > scoreCutoff){ ++hcLinkToAdd; }
  }
  // get the last pos if appropriate
  if (linkEndPos < nucEndPos and sourceSeq->scoreAtPosPlus(linkEndPos) > scoreCutoff){ ++hcNucToAdd; }
  nucThenLinkThenBg[0] = hcNucToAdd;
  nucThenLinkThenBg[1] = hcLinkToAdd;
  nucThenLinkThenBg[2] = nucEndPos - startPos;
  nucThenLinkThenBg[3] = linkEndPos - startPos;
}



float RepeatDetector::getWeightedPercentileScore(float percentile, ContigFragmentWindow** sortedWindows, long numWindows, float normNumWindows){
  if (numWindows == 0){ return 0; }
  else {

    float countLimit = (float(1) - percentile) * normNumWindows;
    long wIndex = 0;
    float countTally = sortedWindows[0]->normCount();
    while (countTally < countLimit and wIndex + 1 < numWindows){
      ++wIndex;
      countTally += sortedWindows[wIndex]->normCount();
    }
    // account for intermediate values (percentile is at the edge between two bins)
    if (countTally == countLimit and wIndex < numWindows+1){
      return (sortedWindows[wIndex]->normScore() + sortedWindows[wIndex+1]->normScore()) / 2;
    } else {
      return sortedWindows[wIndex]->normScore();
    }
  }
}


void RepeatDetector::sortWindowArray(ContigFragmentWindow** windowArray, long wCount){
  vector<ContigFragmentWindow*> sortedWindows;
  for (long n = 0; n < wCount; ++n){ sortedWindows.push_back( windowArray[n] ); }
  SortWindowsByCount sorter; // sort by seq length, longest to shortest
  sort( sortedWindows.begin(), sortedWindows.end(), sorter );
  long index = 0;
  for (vector<ContigFragmentWindow*>::iterator it = sortedWindows.begin(); it != sortedWindows.end(); ++it){
    windowArray[index] = *it;
    index++;
  }
}

bool RepeatDetector::SortWindowsByCount::operator() (ContigFragmentWindow* wA, ContigFragmentWindow* wB){
  return wA->normScore() > wB->normScore();
}



RepeatDetector::ContigWithStats::ContigWithStats(ScoredSeq* seq, long windowSize, RepeatDetectorThreadedness threadedness) :
  _contig(seq) {

  // i am going to use non-overlapping windows
  long cSize = seq->size();
  _numWindows = cSize / windowSize;
  // include incomplete windows, but their contribution to mean/median
  // will be normalized to their size
  if ( cSize % windowSize != 0 ){ ++_numWindows; }

  float* scores = seq->getScores('+');
  float* links = seq->getLinks('+');

  _scoreMeans = new float[ _numWindows + 1 ];
  _windowNorm = new float[ _numWindows + 1 ];
  _windowStarts = new long[ _numWindows + 1 ];
  _windowSizes = new long[ _numWindows + 1 ];

  // i.e. times two (scores plus linkages, including the trailing link for non-terminal windows)
  float maxLengthDenom = float(windowSize + windowSize);

  if (threadedness == THREADED){
    #pragma omp parallel for schedule(dynamic)
    for (long winNum = 0; winNum < _numWindows; ++winNum){
      windowStatHelper(winNum,windowSize,cSize,scores,links,maxLengthDenom);
    }
  } else {
    for (long winNum = 0; winNum < _numWindows; ++winNum){
      windowStatHelper(winNum,windowSize,cSize,scores,links,maxLengthDenom);
    }
  }

  // this is looped separately to avoid threading issues (these would have required atomic
  // statements within the helper function which would have separated the threading declaration
  // from its safety features
  _scoreSum = 0;
  _windowCount = 0;
  for (long winNum = 0; winNum < _numWindows; ++winNum){
    _scoreSum += _scoreMeans[winNum];
    _windowCount += _windowNorm[winNum];
  }

  delete [] scores;
  delete [] links;
}
RepeatDetector::ContigWithStats::~ContigWithStats(){
  delete [] _scoreMeans;
  delete [] _windowNorm;
  delete [] _windowStarts;
  delete [] _windowSizes;
}



RepeatDetector::ContigFragmentWindow* RepeatDetector::ContigWithStats::getWindow(long windowNum){
  if (windowNum >= _numWindows){ throw AssemblyException::ArgError("RD::CWS::getWindow index was too high"); }
  return new ContigFragmentWindow(this, windowNum);
}
inline void RepeatDetector::ContigWithStats::windowStatHelper(long winNum, long windowSize, long cSize, float* scores, float *links, float maxLengthDenom){

  // determine the edges of the window and determine if there is a terminal link
  long winStart = winNum * windowSize;
  long winEnd = winStart + windowSize;
  bool includeFinalLink = true;
  if (cSize <= winEnd){
    winEnd = cSize;
    includeFinalLink = false;
  }
  _windowStarts[winNum] = winStart;
  _windowSizes[winNum] = winEnd - winStart;

  // tally the appropriate scores
  long winEndM1 = winEnd - 1;
  float scoreTally = 0;
  long lengthDenom = (winEnd - winStart) * 2 + 1;
  for (long wPos = winStart; wPos < winEndM1; ++wPos){ scoreTally += scores[wPos] + links[wPos]; }
  scoreTally += scores[winEndM1];
  if (includeFinalLink){
    scoreTally += links[winEndM1];
    ++lengthDenom;
  }

  // calculate the actual values and update the meta-values
  _scoreMeans[winNum] = scoreTally / float(lengthDenom);
  _windowNorm[winNum] = float(lengthDenom) / maxLengthDenom;
}



ScoredSeq** RepeatDetector::ContigWithStats::getHighCoverageFragments(float windowScoreCutoff,
								      EmissionCarrier* scoreEmissions,
								      EmissionCarrier* linkEmissions){
  bool lastWindowKept = false;
  long fragCount = 0;
  // these will mark the windows with high counts to be included
  long* hcFragIndexStarts = new long[ _numWindows+1 ];
  long* hcFragIndexLengths = new long[ _numWindows+1 ];

  // define continuous high-coverage blocks
  for (long winNum = 0; winNum < _numWindows; ++winNum){
    if (_scoreMeans[winNum] > windowScoreCutoff){
      // the previous index's count must be modified, so -1
      if (lastWindowKept){ ++hcFragIndexLengths[fragCount-1]; }
      else {
	lastWindowKept = true;
	hcFragIndexStarts[fragCount] = winNum;
	hcFragIndexLengths[fragCount] = 1;
	++fragCount;
      }
    } else { lastWindowKept = false; }
  }

  // now define the regions across which the HMM should be parsed and use the HMM
  // to define the precise high-coverage boundaries
  vector<ContigFragment*> allLegitFragments;

  // do not thread this, it is externally threaded by contig
  for (long fragNum = 0; fragNum < fragCount; ++fragNum){
    CoverageState beginState = LOWCOVERAGE;
    long beginIndex = hcFragIndexStarts[fragNum];
    if (beginIndex == 0){ beginState = AMBIGUOUSCOVERAGE; }
    else { --beginIndex; }

    CoverageState endState = LOWCOVERAGE;
    long endIndex = hcFragIndexStarts[fragNum] + hcFragIndexLengths[fragNum];
    if (endIndex == _numWindows){ endState = AMBIGUOUSCOVERAGE; }
    else { ++endIndex; }
    long endIndexM1 = endIndex - 1;

    if (hcFragIndexLengths[fragNum] <= 2){
      // parse the whole block all at once
      long startCoord = _windowStarts[beginIndex];
      long endCoord = _windowStarts[endIndexM1] + _windowSizes[endIndexM1];
      float transitionProb;
      if (startCoord==endCoord){ transitionProb = 0; }
      else { transitionProb = float(2) / float(endCoord-startCoord); }
      parseCoverageHmm(&allLegitFragments,startCoord,endCoord,beginState,endState,transitionProb,
		       scoreEmissions,linkEmissions);
    } else {
      // parse the beginning of the block
      long numWinBegin = 2;
      if (beginState == AMBIGUOUSCOVERAGE){ numWinBegin = 1; }
      long startCoordBegin = _windowStarts[beginIndex];
      long endCoordBegin = _windowStarts[beginIndex+numWinBegin] + _windowSizes[beginIndex+numWinBegin];
      float transitionProbBegin;
      if (startCoordBegin==endCoordBegin){ transitionProbBegin = 0; }
      else { transitionProbBegin = float(1) / float(endCoordBegin-startCoordBegin); }
      parseCoverageHmm(&allLegitFragments, startCoordBegin, endCoordBegin, beginState, HIGHCOVERAGE, transitionProbBegin,
		       scoreEmissions,linkEmissions);

      // parse the end of the block
      long numWinEnd = 2;
      if (endState == AMBIGUOUSCOVERAGE){ numWinEnd = 1; }
      long startCoordEnd = _windowStarts[endIndex-numWinEnd];
      long endCoordEnd = _windowStarts[endIndexM1] + _windowSizes[endIndexM1];
      float transitionProbEnd;
      if (startCoordEnd==endCoordEnd){ transitionProbEnd = 0; }
      else { transitionProbEnd = float(1) / float(endCoordEnd-startCoordEnd); }
      parseCoverageHmm(&allLegitFragments, startCoordEnd, endCoordEnd, HIGHCOVERAGE, endState, transitionProbEnd,
		       scoreEmissions,linkEmissions);

      // if there is an unparsed intervening segment, add in the complete region as one high-coverage block
      if (startCoordEnd > endCoordBegin){
	allLegitFragments.push_back( new ContigFragmentSimple(_contig, endCoordBegin, startCoordEnd - endCoordBegin) );
      }
    }
  }
  delete [] hcFragIndexStarts;
  delete [] hcFragIndexLengths;


  // collapse redundant HC blocks
  // 1) sort by start position
  SortFragmentsByStart sortByStart;
  sort( allLegitFragments.begin(), allLegitFragments.end(), sortByStart );
  // 2) copy to a new vector, removing overlap redundancy
  vector<ContigFragment*> nrHcFragments;
  vector<ContigFragment*>::iterator fragIt = allLegitFragments.begin();
  if (fragIt != allLegitFragments.end()){
    nrHcFragments.push_back( *fragIt );
    ++fragIt;
  }
  while (fragIt != allLegitFragments.end()){
    vector<ContigFragment*>::iterator currentIt = nrHcFragments.end();
    currentIt--;
    long endCurrent = (*currentIt)->startPos() + (*currentIt)->size();
    if ( endCurrent < (*fragIt)->startPos() ){ nrHcFragments.push_back( *fragIt ); }
    else {
      long start = (*currentIt)->startPos();
      long endNext = (*fragIt)->startPos() + (*fragIt)->size();
      if (endNext > endCurrent){
	delete *currentIt;
	*currentIt = new ContigFragmentSimple(_contig, start, endNext - start);
      }
      delete *fragIt;
    }
    ++fragIt;
  }

  // convert the fragment objects into sequences for output
  ScoredSeq** returnSeqs = new ScoredSeq*[ nrHcFragments.size()+1 ];
  long seqIndex = 0;
  for (vector<ContigFragment*>::iterator fIt = nrHcFragments.begin(); fIt != nrHcFragments.end(); ++fIt){
    ContigFragment* frag = *fIt;
    long startPos = frag->startPos();
    long size = frag->size();
    delete frag;
    char* subseq = _contig->getSubseq(startPos, size, '+');
    float* scores = _contig->getSubScores(startPos, size, '+');
    float* links = _contig->getSubLinks(startPos, size-1, '+');
    returnSeqs[seqIndex] = ScoredSeq::repExposedSeq(subseq, scores, links, size);
    ++seqIndex;
  }

  returnSeqs[seqIndex] = NULL;
  return returnSeqs;
}



bool RepeatDetector::SortFragmentsByStart::operator() (ContigFragment* fragA, ContigFragment* fragB){
  return fragA->startPos() < fragB->startPos();
}


RepeatDetector::EmissionCarrier::EmissionCarrier(float scoreCutoff, float hcProbLcRegion, float hcProbHcRegion) :
  _scoreCutoff(scoreCutoff)
{
  _highStateEmission[1] = log(hcProbHcRegion);
  _lowStateEmission[1] = log(hcProbLcRegion);
  _highStateEmission[0] = log(float(1) - hcProbHcRegion);
  _lowStateEmission[0] = log(float(1) - hcProbLcRegion);
}
RepeatDetector::EmissionCarrier::~EmissionCarrier(){}

void RepeatDetector::ContigWithStats::parseCoverageHmm(vector<ContigFragment*>* outputFragments, long beginIndex, long endIndex,
						       CoverageState beginState, CoverageState endState,
						       float changeProb, EmissionCarrier* scoreEmissions, EmissionCarrier* linkEmissions){

  // if this is not the case, the returned vector will exist but be empty
  if (endIndex > beginIndex){
    // zeroth element is the begin state (don't parse back to it)
    // this also reflects linkages as positions, though their parsing is irrelevant
    // the final link will be ignored no matter what
    long parseLen = (endIndex - beginIndex) * 2;

    // convention in this function: 0 => low coverage, 1 => high coverage
    int* priorState[2];
    float* probability[2];
    for (int n = 0; n < 2; ++n){
      priorState[n] = new int[ parseLen+2 ];
      probability[n] = new float[ parseLen+2 ];
    }

    // for fast access of probabilities
    // 0+0=0, 0+1=1, 1+0=1, 1+1=2 <= adding the two states gives the correct probability
    float stayTheSameProb = float(1) - changeProb;
    // these must be logs for adding
    float transitionProb[3] = {log(stayTheSameProb), log(changeProb), log(stayTheSameProb)};

    // set up the initial cells
    switch (beginState){
    case LOWCOVERAGE: probability[0][0] = 0; probability[1][0] = -1000000; break;
    case HIGHCOVERAGE: probability[0][0] = -1000000; probability[1][0] = 0; break;
    case AMBIGUOUSCOVERAGE: probability[0][0] = log(0.5); probability[1][0] = log(0.5); break;
    default: throw AssemblyException::ArgError("RD:CWS this isnt going to happen");
    }

    // get arrays of the scores for easy access
    float* scoresAndLinks[2];
    scoresAndLinks[0] = _contig->getSubScores(beginIndex, endIndex-beginIndex, '+');
    if (endIndex == _contig->size()){
      scoresAndLinks[1] = _contig->getSubLinks(beginIndex, endIndex-beginIndex-1, '+');
      --parseLen;
    } else {
      scoresAndLinks[1] = _contig->getSubLinks(beginIndex, endIndex-beginIndex, '+');
    }

    // an array for the two emission models so that they can be accessed immediately based on remainder
    EmissionCarrier* emissions[2];
    emissions[0] = scoreEmissions;
    emissions[1] = linkEmissions;

    // remember that pN is one behind the current HMM array position
    for (long pN = 0; pN < parseLen; ++pN){
      long pNp1 = pN + 1;
      int emType;
      long emIndex = pN % 2;
      if (scoresAndLinks[emIndex][pN / 2] > emissions[emIndex]->_scoreCutoff){ emType = 1; }
      else { emType = 0; }
      float emissionProb[2] = { emissions[emIndex]->_lowStateEmission[emType], emissions[emIndex]->_highStateEmission[emType] };

      for (int state = 0; state < 2; ++state){
	// the two cells are the two possible prior states
	float newProb[] = { probability[0][pN] + transitionProb[0+state] + emissionProb[state],
			    probability[1][pN] + transitionProb[1+state] + emissionProb[state] };
	int betterProbIndex;
	if (newProb[0] < newProb[1]){ betterProbIndex = 1; }
	else { betterProbIndex = 0; }
	priorState[state][pNp1] = betterProbIndex;
	probability[state][pNp1] = newProb[ betterProbIndex ];
      }
    }
    delete [] scoresAndLinks[0];
    delete [] scoresAndLinks[1];

    // figure out the final state for the traceback
    int stateIndex;
    switch (endState){
    case LOWCOVERAGE: stateIndex = 0; break;
    case HIGHCOVERAGE: stateIndex = 1; break;
    case AMBIGUOUSCOVERAGE:
      if (probability[0][parseLen] < probability[1][parseLen]){ stateIndex = 1; }
      else { stateIndex = 0; }
      break;
    default: throw AssemblyException::ArgError("RD:CWS this isnt going to happen");
    }

    // these three variables keep track of when an HC block is encountered
    bool hcBlockOpen = false;
    long hcBlockEnd;
    if (stateIndex == 1){
      hcBlockOpen = true;
      hcBlockEnd = (parseLen+1) / 2 + beginIndex;
    }

    // do the traceback
    for (long pN = parseLen-1; pN > 0; --pN){
      if (stateIndex != priorState[stateIndex][pN]){
	long hcBlockStart;
	switch (stateIndex){
	case 0:
	  hcBlockOpen = true;
	  hcBlockEnd = pN / 2 + beginIndex;
	  break;
	case 1:
	  hcBlockOpen = false;
	  hcBlockStart = pN / 2 + beginIndex;
	  outputFragments->push_back( new ContigFragmentSimple(_contig, hcBlockStart, hcBlockEnd - hcBlockStart) );
	  break;
	default:
	  throw AssemblyException::LogicError("RD::hmmParse state should be 0 or 1");
	}
	stateIndex = priorState[stateIndex][pN];
      }
    }

    // if an hc block is still open, close it and add it to the collection
    if (hcBlockOpen){
      outputFragments->push_back( new ContigFragmentSimple(_contig, beginIndex, hcBlockEnd - beginIndex) );
    }

    // clean up
    for (int n = 0; n < 2; ++n){
      delete [] priorState[n];
      delete [] probability[n];
    }
  }
}

RepeatDetector::ContigFragment::~ContigFragment(){}


RepeatDetector::ContigFragmentWindow::ContigFragmentWindow(ContigWithStats* cWithStats, long windowIndex) :
  _contig(cWithStats),
  _windowIndex(windowIndex)
{}
RepeatDetector::ContigFragmentWindow::~ContigFragmentWindow(){}
float RepeatDetector::ContigFragmentWindow::score(){ return _contig->_scoreMeans[_windowIndex]; }
float RepeatDetector::ContigFragmentWindow::normScore(){
  return _contig->_scoreMeans[_windowIndex] / _contig->_windowNorm[_windowIndex];
}
float RepeatDetector::ContigFragmentWindow::normCount(){ return _contig->_windowNorm[_windowIndex]; }
long RepeatDetector::ContigFragmentWindow::startPos(){ return _contig->_windowStarts[_windowIndex]; }
long RepeatDetector::ContigFragmentWindow::size(){ return _contig->_windowSizes[_windowIndex]; }
ScoredSeq* RepeatDetector::ContigFragmentWindow::sourceSeq(){ return _contig->_contig; }


RepeatDetector::ContigFragmentSimple::ContigFragmentSimple(ScoredSeq* sourceSeq, long startPos, long size):
  _contig(sourceSeq),
  _startPos(startPos),
  _size(size)
{}
RepeatDetector::ContigFragmentSimple::~ContigFragmentSimple(){}
long RepeatDetector::ContigFragmentSimple::startPos(){ return _startPos; }
long RepeatDetector::ContigFragmentSimple::size(){ return _size; }
ScoredSeq* RepeatDetector::ContigFragmentSimple::sourceSeq(){ return _contig; }


#endif
