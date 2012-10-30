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

/* Makes alignments to edges; terminates early if no alignment is found.
*/



#ifndef DYNAMICPROGRAMMINGALIGNER_H
#define DYNAMICPROGRAMMINGALIGNER_H

#include <set>
#include <map>
#include "AlignmentScoreMatrix.h"
#include "AlignmentGapped.h"
#include "ScoredSeq.h"
using namespace::std;

class DynamicProgrammingAligner {


 public:
  DynamicProgrammingAligner();
  // default min score is zero
  DynamicProgrammingAligner(float fractId, long minOverlap, AlignmentScoreMatrix * scoreMatrix);
  DynamicProgrammingAligner(float fractId, long minOverlap, long maxOverlap, AlignmentScoreMatrix * scoreMatrix);
  DynamicProgrammingAligner(float fractId, long minOverlap, AlignmentScoreMatrix * scoreMatrix, long minScore);
  DynamicProgrammingAligner(float fractId, long minOverlap, long maxOverlap, AlignmentScoreMatrix * scoreMatrix, long minScore);
  DynamicProgrammingAligner(DynamicProgrammingAligner* dpa);
  ~DynamicProgrammingAligner();

  // get the input info (useful for the ScoredSeqCollection)
  float getFractId();
  long getMinOverlap();
  long getMinScore();
  AlignmentScoreMatrix * getScoreMatrix();
  DynamicProgrammingAligner * copy();

  // input max overlap constrains, but does not weaken, global max overlap
  Alignment * align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation);
  Alignment * align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long maxOverlap);
  Alignment * align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long minOffset, long maxOffset);
  Alignment * align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, set<long>* inputOffsets);
  Alignment * align(ScoredSeq* extSeqA, ScoredSeq* extSeqB, char orientation, long* offsetArray, long numOffsets);

  // make a copy that changes one or more of the input params; copy is shallow
  DynamicProgrammingAligner * minScoreCopy(long minScore);

  void OK();


 private:
  AlignmentScoreMatrix * _scoreMatrix;
  float _fractId;
  long _minOverlap;
  bool _usingMaxOverlap;
  long _maxOverlap; // less than 0 means no max overlap
  // this is an externally-imposed minScore that is zero by default
  long _minScore;


  class Square {
  public:
    Square(long posA, long posB, long score);
    Square(long posA, long posB, ScoredSeq* seqA, ScoredSeq* seqB, char senseA, char senseB, Square* origin, AlignmentScoreMatrix* scoreMatrix);
    Square(long posA, long posB, long extraNucPos, ScoredSeq* extraNucSeq, char extraNucSense, Square * origin, AlignmentScoreMatrix * scoreMatrix);
    ~Square();
    //static long _count;
    bool isMatch();
    bool isGap();
    bool okToDelete(); // no children, not best, not in grid
    long _score;
    bool _isMatch;
    long _posA;
    long _posB;
    bool _isBest;
    bool _isInGrid;
    long _sumMatch;
    long _sumSkip;
    // the gaps and mismatches are never considered independently
    long _sumGapAndMis;
    // sum of all the gaps and aligned positions so far in the alignment
    long _denominator;
    ScoredSeq * _extraNucSeq;
    Square * _origin;
    bool _visitedTargetOffset;
    // this one allows recursive deletion to be run iteratively
    bool _deleteDeep;
  private:
    int _hasChildSum;
  };

  class EdgeBoundaryCarrier {
  public:
    EdgeBoundaryCarrier(long minOffset, long maxOffset);
    EdgeBoundaryCarrier(long minOffset, long maxOffset, long targetOffset);
    ~EdgeBoundaryCarrier();
    // i am leaving these exposed for expediency
    long _minOffset;
    long _maxOffset;
    long* _distToOkOffset;
    // these methods are for manipulating and combining Ebc's
    bool isContinuous(EdgeBoundaryCarrier* otherEbc);
    EdgeBoundaryCarrier* combine(EdgeBoundaryCarrier* otherEbc);
    bool overlaps(EdgeBoundaryCarrier* otherEbc);
    EdgeBoundaryCarrier* getOverlap(EdgeBoundaryCarrier* otherEbc);
  };

  class EdgeBoundaryCarrierCollection {
  public:
    EdgeBoundaryCarrierCollection();
    EdgeBoundaryCarrierCollection(long size);
    ~EdgeBoundaryCarrierCollection();
    EdgeBoundaryCarrier** _ebcArray;
    long _size;
  };

  // this class is designed to optimize the functions of both finding
  // a square in a line and iterating throught the non-null squares from a line
  // of the dynamic programming grid.
  class DpaLine {
  public:
    //DpaLine();
    DpaLine(long minIndex, long maxIndex, long extraCells);
    ~DpaLine();
    // instead of creating a new one
    void recycle(long minIndex, long maxIndex);

    // for getting or modifying squares
    bool hasSquare(long posB);
    void addSquare(Square* square, long posB);
    Square* getSquare(long posB);
    // number of active elements
    long size();

    Square** _squareArray;
    long _rowSize;
    long* _activeCells;
    long _numActive;
    // these define the limits of the row's arrays
    long _minIndex;
    long _maxIndex;
    // these will define limits for the next row's arrays
    long _currentMin;
    long _currentMax;
    bool _currentSet;

    class Iterator{
    public:
      //Iterator();
      Iterator(DpaLine* source, long itemNum);
      Iterator(const Iterator& it);
      Iterator operator++();
      bool operator==(const Iterator& it);
      bool operator!=(const Iterator& it);
      Square* operator*();
      long _itemNum;
      DpaLine* _source;
    };

    // for iterating through the non-null elements
    DpaLine::Iterator begin();
    DpaLine::Iterator end();
  };


  // private helper functions
  //Alignment * alignPrivate(ScoredSeq * extSeqA, ScoredSeq * extSeqB, char orientation, int maxOverlap, bool useMinScore, int minScore);

  // REQUIRES: seqA is as long or longer than seqB
  EdgeBoundaryCarrierCollection* getOverlapLimits(ScoredSeq* seqA, ScoredSeq* seqB, long minOverlap, long maxOverlap);
  // WORKS BEST IF: seqA is as long or longer than seqB
  Square * alignHelp(ScoredSeq * seqA, ScoredSeq * seqB, char senseA, char senseB, EdgeBoundaryCarrier* ebc);
  // helps align; delete oldEbcc and its redundant content
  EdgeBoundaryCarrierCollection* makeEbccNr(EdgeBoundaryCarrierCollection* oldEbcc);
  void sortEbccByOffsets(EdgeBoundaryCarrierCollection* ebcc, long* sortedOffsets, long numOffsets);

  static void alignmentFromSquare(Square* bestSquare, AlignmentGapped* alignment, ScoredSeq* seqA, ScoredSeq* seqB);
  // the one on top is the old one and should be removed when the new one works
  Square* evaluateNewPosition(Square* candNewSquare, DpaLine* nextLine, long seqAsizeM1, long seqBsizeM1,
			      Square* bestSquare, EdgeBoundaryCarrier* ebc, long maxNumMismatch);
  Square* getNewBestRef(Square* bestSquare, Square* bestCand);
  bool isSquareLegit(Square* newSquare, long seqAsizeM1, long seqBsizeM1, EdgeBoundaryCarrier* ebc, long maxNumMismatch);

  // helper to simplify code dealing with overlap lengths
  long getOverlapLength(long seqLenA, long seqLenB, long offset);
  // helper for defining the longest possible overlap for the window; useful
  // for defining the max number of allowed mismatches for the corridor defined
  // by the provided EBC
  long getMaxOverlapLength(long seqLenA, long seqLenB, EdgeBoundaryCarrier* ebc);
  // this is used at multiple points (corridor width, maxMismatch)
  long adjustMaxMismatch(long maxMismatch);

};

#endif
