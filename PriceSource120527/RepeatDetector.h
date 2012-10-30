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


#ifndef REPEATDETECTOR_H
#define REPEATDETECTOR_H

#include <set>
#include "ScoredSeq.h"
#include "ParamsAlignment.h"
using namespace::std;

class RepeatDetector {
 public:
  RepeatDetector(long minRepeatSize, ParamsAlignment* alignParams);
  ~RepeatDetector();

  enum RepeatDetectorThreadedness{ NOTTHREADED=0, THREADED=1 };

  ScoredSeq** findRepeats(set<ScoredSeq*>* contigs, float numStdDevs, float minFoldIncrease,
			  RepeatDetectorThreadedness threadedness);


 private:
  enum CoverageState{ HIGHCOVERAGE=0, LOWCOVERAGE=1, AMBIGUOUSCOVERAGE=2 };
  long _windowSize;
  long _minRepeatSize;
  ParamsAlignment* _alignParams;

  class ContigFragment;
  class ContigFragmentWindow;
  class EmissionCarrier;

  // a collection of the data structures associated with a contig
  // that are necessary for analysis of sub-sequence scores
  class ContigWithStats {
  public:
    ContigWithStats(ScoredSeq* seq, long windowSize, RepeatDetectorThreadedness threadedness);
    ~ContigWithStats();
    ContigFragmentWindow* getWindow(long windowNum);
    // must end with a NULL element
    ScoredSeq** getHighCoverageFragments(float windowScoreCutoff,
					 EmissionCarrier* scoreEmissions,
					 EmissionCarrier* linkEmissions);

    long _numWindows;
    ScoredSeq* _contig;
    long* _windowStarts;
    long* _windowSizes;

    float _windowCount;
    float _scoreSum;

    float* _scoreMeans;
    float* _windowNorm;
  private:
    void windowStatHelper(long winNum, long windowSize, long cSize, float* scores, float* links, float maxLengthDenom);
    void parseCoverageHmm(vector<ContigFragment*>* outputFragments, long beginIndex, long endIndex,
			  CoverageState beginState, CoverageState endState,
			  float changeProb, EmissionCarrier* scoreEmissions, EmissionCarrier* linkEmissions);
  };

  struct SortFragmentsByStart {
    bool operator() (ContigFragment* fragA, ContigFragment* fragB);
  };

  class EmissionCarrier {
  public:
    EmissionCarrier(float scoreCutoff, float hcProbLcRegion, float hcProbHcRegion);
    ~EmissionCarrier();
    // publicly-accessible fields
    float _scoreCutoff;
    // these values are the logs (will be added)
    float _lowStateEmission[2];
    float _highStateEmission[2];
  };

  // an abstract class that represents a subset of a contig (implementations will be
  // optimized based on context)
  class ContigFragment {
  public:
    virtual ~ContigFragment();
    //virtual float score() = 0;
    //virtual float normCount() = 0;
    virtual long startPos() = 0;
    virtual long size() = 0;
    virtual ScoredSeq* sourceSeq() = 0;
  };

  class ContigFragmentSimple : public ContigFragment {
  public:
    ContigFragmentSimple(ScoredSeq* sourceSeq, long startPos, long size);
    ~ContigFragmentSimple();
    long startPos();
    long size();
    ScoredSeq* sourceSeq();
  private:
    ScoredSeq* _contig;
    long _startPos;
    long _size;
  };

  class ContigFragmentWindow : public ContigFragment {
  public:
    ContigFragmentWindow(ContigWithStats* cWithStats, long windowIndex);
    ~ContigFragmentWindow();
    long startPos();
    long size();
    ScoredSeq* sourceSeq();
    // class-specific
    float score();
    float normScore();
    float normCount();
  private:
    ContigWithStats* _contig;
    long _windowIndex;
    long _windowStart;
    long _windowSize;
  };


  float getWeightedPercentileScore(float percentile, ContigFragmentWindow** sortedWindows, long numWindows, float normNumWindows);
  void getHcCounts(long* nucThenLink, ContigFragmentWindow* window, float scoreCutoff);
  //float getWeightedMeanScore(ContigWindow** sortedWindows, long numWindows);



  void sortWindowArray(ContigFragmentWindow** windowArray, long wCount);
  struct SortWindowsByCount {
    bool operator() (ContigFragmentWindow* wA, ContigFragmentWindow* wB);
  };

};

#endif
