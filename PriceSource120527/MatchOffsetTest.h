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

#ifndef MATCHOFFSETTEST_H
#define MATCHOFFSETTEST_H

#include "ScoredSeq.h"
using namespace::std;

// this is a pure virtual class for tests that can be passed to
// ScoredSeqCollectionBwt to evaluate sequences' appropriateness for
// return based on properties other than the potential alignments they
// might form

// these classes are going to accelerate the assembler by performing tests
// on potential match sequences prior to generating an alignment rather than
// waiting until afterwards.  provided implementations are below.
class MatchOffsetTest {
 public:
  virtual ~MatchOffsetTest();
  virtual bool offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB) = 0;
};


// i am including the implementing classes in this file

// this one always returns true
class MatchOffsetTestNull : public MatchOffsetTest {
 public:
  MatchOffsetTestNull();
  ~MatchOffsetTestNull();
  bool offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB);
};

// makes sure a seq is long enough
class MatchOffsetTestFullOverlapA : public MatchOffsetTest {
 public:
  MatchOffsetTestFullOverlapA();
  ~MatchOffsetTestFullOverlapA();
  bool offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB);
};

// seq MUST be in the set
class MatchOffsetTestMinOverlap : public MatchOffsetTest {
 public:
  MatchOffsetTestMinOverlap();
  MatchOffsetTestMinOverlap(long minOverlap);
  ~MatchOffsetTestMinOverlap();
  bool offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB);
 private:
  long _minOverlap;
};

// imposes multiple tests
class MatchOffsetTestMulti : public MatchOffsetTest {
 public:
  MatchOffsetTestMulti();
  MatchOffsetTestMulti(MatchOffsetTest** testArray, int numTests);
  ~MatchOffsetTestMulti();
  bool offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB);
 private:
  MatchOffsetTest** _testArray;
  int _numTests;
};


#endif


