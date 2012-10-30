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
This interface provides the opportunity to filter away potentially problematic input data.  It provides a test for
for both an individual sequence AND a pair of sequences (it is applied to paired-end sequences), and informs the 
user whether that sequence (or pair) passes or fails the test.  Generally, are used to eliminate reads or input contigs
with very low sequence complexity, overall low quality scores, etc.  The provision of a single test for a pair of sequences
allows for the boolean combination of the single-sequence test to be adjusted according to the test (i.e. does the pair
pass if either sequence passes, or only if both sequences pass?) and also permits the implementation of tests that apply
to the relationship of reads to one another (for example, keeping or eliminating pairs for which the two sequences extensively
overlap, indicating a short underlying amplicon).</li>
*/

#ifndef READPAIRFILTER_H
#define READPAIRFILTER_H

#include <set>
#include "ScoredSeq.h"
#include "ParameterizedInitialFile.h"
#include "ScoredSeqCollectionBwt.h"
using namespace::std;

// this is a pure virtual class that will allow the contigs generated
// at the end of each cycle to be filtered according to some arbitrary
// criteria set by the implementing classes, which are below.

class ReadPairFilter {
 public:
  virtual ~ReadPairFilter();
  virtual ReadPairFilter* copy() = 0;
  // allows the filters to be parsed before being called (important because they will
  // be called such a huge number of times).  cycles are indexed from zero
  virtual bool useThisCycle(int cycleNum) = 0;
  // can decide whether or not a test will be applied to a sequence based on that seq's length;
  // this is not considered by isReadOk!!!!
  virtual bool applyToLength(long seqSize) = 0;
  virtual bool testKillsPair() = 0;
  virtual bool isReadOk(ScoredSeq* read) = 0;
};


// i am including the implementing classes in this file

// this one always returns true for if the sequence is OK
// and always false for if it should be run
class ReadPairFilterNull : public ReadPairFilter {
 public:
  ReadPairFilterNull();
  ~ReadPairFilterNull();
  ReadPairFilterNull* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
};

// prevents long strings of single nucleotides
class ReadPairFilterHomoPolymer : public ReadPairFilter {
 public:
  ReadPairFilterHomoPolymer(long maxHomoPolymer);
  ReadPairFilterHomoPolymer(long maxHomoPolymer, long maxSeqLength);
  ~ReadPairFilterHomoPolymer();
  ReadPairFilterHomoPolymer* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  long _maxHomoPolymer;
  bool _allLengths;
  // only has a value if _allLengths==false
  long _maxLength;
};

// prevents long strings of a dinucleotide repeat
class ReadPairFilterDinucleotide : public ReadPairFilter {
 public:
  // max length can be an even or odd number; the edge of the repeat can be at either
  // nucleotide at either end
  ReadPairFilterDinucleotide(long maxHomoPolymer);
  ReadPairFilterDinucleotide(long maxHomoPolymer, long maxSeqLength);
  ~ReadPairFilterDinucleotide();
  ReadPairFilterDinucleotide* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  long _maxDinucleotide;
  bool _allLengths;
  // only has a value if _allLengths==false
  long _maxLength;
};

// gets rid of sequences that match some target
class ReadPairFilterAvoidSeqs : public ReadPairFilter {
 public:
  ReadPairFilterAvoidSeqs(set<ScoredSeq*>* seqsToAvoid, float minFractId);
  ReadPairFilterAvoidSeqs(set<ScoredSeq*>* seqsToAvoid, float minFractId, long maxSeqLength);
  ~ReadPairFilterAvoidSeqs();
  ReadPairFilterAvoidSeqs* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  ScoredSeqCollectionBwt* _seqsToAvoid;
  bool _allLengths;
  // only has a value if _allLengths==false
  long _maxLength;
};

// remove low-quality sequences
class ReadPairFilterQuality : public ReadPairFilter {
 public:
  ReadPairFilterQuality(float maxFractBad, float scoreThreshold);
  ReadPairFilterQuality(float maxFractBad, float scoreThreshold, int numSkipCycles, int numRunCycles);
  ReadPairFilterQuality(float maxFractBad, float scoreThreshold, long maxSeqLength);
  ReadPairFilterQuality(float maxFractBad, float scoreThreshold, int numSkipCycles, int numRunCycles, long maxSeqLength);
  ~ReadPairFilterQuality();
  ReadPairFilterQuality* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  float _maxFractBad;
  float _scoreThreshold;
  // the cycle fields only have values if the bool is 'false'
  bool _useEveryCycle;
  int _firstCycle;
  int _lastCycleP1;
  bool _allLengths;
  // only has a value if _allLengths==false
  long _maxLength;
};

// remove low-quality sequences based on the number of 
class ReadPairFilterUncalledBases : public ReadPairFilter {
 public:
  ReadPairFilterUncalledBases(float maxFractBad);
  ReadPairFilterUncalledBases(float maxFractBad, int numSkipCycles, int numRunCycles);
  ReadPairFilterUncalledBases(float maxFractBad, long maxSeqLength);
  ReadPairFilterUncalledBases(float maxFractBad, int numSkipCycles, int numRunCycles, long maxSeqLength);
  ~ReadPairFilterUncalledBases();
  ReadPairFilterUncalledBases* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  float _maxFractBad;
  // the cycle fields only have values if the bool is 'false'
  bool _useEveryCycle;
  int _firstCycle;
  int _lastCycleP1;
  bool _allLengths;
  // only has a value if _allLengths==false
  long _maxLength;
};

// imposes multiple tests
class ReadPairFilterMulti : public ReadPairFilter {
 public:
  ReadPairFilterMulti(ReadPairFilter** testArray, int numTests);
  ~ReadPairFilterMulti();
  ReadPairFilterMulti* copy();
  bool useThisCycle(int cycleNum);
  bool applyToLength(long seqSize);
  bool testKillsPair();
  bool isReadOk(ScoredSeq* read);
 private:
  ReadPairFilter** _testArray;
  int _numTests;
  bool _testKillsPair;
};


#endif


