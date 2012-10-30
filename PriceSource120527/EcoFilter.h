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
This interface represents a test for contigs (ScoredSeqs), and implementing classes take in a set of sequences and
sort them into two categories: "passing" or "not passing".  Generally, these tests are used to identify contigs
that are not of interest and therefore should be deleted at the end of an assembly cycle.  Current classes implement
tests for minimum length, minimum coverage levels, matches to reference sequences, etc.</li>
*/

#ifndef ECOFILTER_H
#define ECOFILTER_H

#include <set>
#include "ScoredSeq.h"
#include "ParameterizedInitialFile.h"
#include "AlignmentScoreMatrix.h"
using namespace::std;

// this is a pure virtual class that will allow the contigs generated
// at the end of each cycle to be filtered according to some arbitrary
// criteria set by the implementing classes, which are below.

class EcoFilter {
 public:
  virtual ~EcoFilter();
  virtual EcoFilter* copy() = 0;
  enum FilterThreadedness{ notThreaded=0, threaded=1 };
  virtual void filterSeqs(set<ScoredSeq*>* fullSet,
			  set<ScoredSeq*>* retainedSet,
			  set<ScoredSeq*>* removedSet,
			  int cycleNum, FilterThreadedness threadedness) = 0;
};


// i am including the implementing classes in this file

// this one always returns true
class EcoFilterNull : public EcoFilter {
 public:
  EcoFilterNull();
  ~EcoFilterNull();
  EcoFilterNull* copy();
  void filterSeqs(set<ScoredSeq*>* fullSet,
		  set<ScoredSeq*>* retainedSet,
		  set<ScoredSeq*>* removedSet,
		  int cycleNum, FilterThreadedness threadedness);
};

// makes sure a seq has enough score counts
class EcoFilterMinCount : public EcoFilter {
 public:
  EcoFilterMinCount();
  EcoFilterMinCount(float minCount);
  ~EcoFilterMinCount();
  EcoFilterMinCount* copy();
  void filterSeqs(set<ScoredSeq*>* fullSet,
		  set<ScoredSeq*>* retainedSet,
		  set<ScoredSeq*>* removedSet,
		  int cycleNum, FilterThreadedness threadedness);
  static EcoFilterMinCount* defaultFilter();
 private:
  bool isOk(ScoredSeq* seq);
  float _minCount;
};


// makes sure a seq is long enough
class EcoFilterMinLength : public EcoFilter {
 public:
  EcoFilterMinLength();
  EcoFilterMinLength(long minLength, int skipCycles);
  ~EcoFilterMinLength();
  // modify it to add a new min length - won't modify specified
  // cycles after it - will replace redundant entries
  void addMinLength(long minLength, int skipCycles);
  EcoFilterMinLength* copy();
  void filterSeqs(set<ScoredSeq*>* fullSet,
		  set<ScoredSeq*>* retainedSet,
		  set<ScoredSeq*>* removedSet,
		  int cycleNum, FilterThreadedness threadedness);
 private:
  // at cycles beyond this one, the final value is used
  long _maxCycleNum;
  // zero or one can signify no filtering
  long* _minLenByCycle;
  // marks cycles at which the length filter changes
  bool* _isTransition;
  // private constructor for the copy method's use
  EcoFilterMinLength(EcoFilterMinLength* efml);
};


// seq must match among the initial contigs that have been included so far
class EcoFilterInitialContigMatch : public EcoFilter {
 public:
  EcoFilterInitialContigMatch();
  EcoFilterInitialContigMatch(float minFractId, int cyclesToSkip, bool fullFile = false);
  EcoFilterInitialContigMatch(float minFractId, int cyclesToSkip,
			      int numFilterCycles, int numSkipCycles, bool fullFile = false);
  ~EcoFilterInitialContigMatch();
  EcoFilterInitialContigMatch* copy();
  void setAlignmentScoreMatrix(AlignmentScoreMatrix* alSM);
  void filterSeqs(set<ScoredSeq*>* fullSet,
		  set<ScoredSeq*>* retainedSet,
		  set<ScoredSeq*>* removedSet,
		  int cycleNum, FilterThreadedness threadedness);
  // class-specific
  void addFiles(set<ParameterizedInitialFile*>* initialContigFiles);
 private:
  set<ParameterizedInitialFile*> _initialContigFiles;
  float _minFractId;
  int _cyclesToSkip;
  int _numFilterCycles;
  int _numSkipCycles;
  bool _fullFile;
  AlignmentScoreMatrix* _asMatrix;
};

// imposes multiple tests
class EcoFilterMulti : public EcoFilter {
 public:
  EcoFilterMulti();
  EcoFilterMulti(EcoFilter** testArray, int numTests);
  ~EcoFilterMulti();
  EcoFilterMulti* copy();
  void filterSeqs(set<ScoredSeq*>* fullSet,
		  set<ScoredSeq*>* retainedSet,
		  set<ScoredSeq*>* removedSet,
		  int cycleNum, FilterThreadedness threadedness);
 private:
  EcoFilter** _testArray;
  int _numTests;
};


#endif


