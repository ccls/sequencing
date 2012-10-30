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

#ifndef MATCHSEQTEST_H
#define MATCHSEQTEST_H

#include <set>
#include "ScoredSeq.h"
using namespace::std;

// this is a pure virtual class for tests that can be passed to
// ScoredSeqCollectionBwt to evaluate sequences' appropriateness for
// return based on properties other than the potential alignments they
// might form

// these classes are going to accelerate the assembler by performing tests
// on potential match sequences prior to generating an alignment rather than
// waiting until afterwards.  provided implementations are below.
class MatchSeqTest {
 public:
  static MatchSeqTest* getNullTest();
  virtual ~MatchSeqTest();
  virtual bool seqIsOk(ScoredSeq* seq) = 0;
};


// i am including the implementing classes in this file

// this one always returns true
class MatchSeqTestNull : public MatchSeqTest {
 public:
  MatchSeqTestNull();
  ~MatchSeqTestNull();
  bool seqIsOk(ScoredSeq* seq);
};

// makes sure a seq is long enough
class MatchSeqTestMinLength : public MatchSeqTest {
 public:
  MatchSeqTestMinLength(long minLength);
  ~MatchSeqTestMinLength();
  bool seqIsOk(ScoredSeq* seq);
 private:
  long _minLength;
};

// looks for an exact length match
class MatchSeqTestExactLength : public MatchSeqTest {
 public:
  MatchSeqTestExactLength(long length);
  ~MatchSeqTestExactLength();
  bool seqIsOk(ScoredSeq* seq);
 private:
  long _length;
};

// seq MUST be in the set
class MatchSeqTestInSet : public MatchSeqTest {
 public:
  MatchSeqTestInSet(set<ScoredSeq*>* seqSet);
  ~MatchSeqTestInSet();
  bool seqIsOk(ScoredSeq* seq);
 private:
  set<ScoredSeq*> _seqSet;
};

// seq MUST NOT be in the set
class MatchSeqTestNotInSet : public MatchSeqTest {
 public:
  MatchSeqTestNotInSet();
  MatchSeqTestNotInSet(set<ScoredSeq*>* seqSet);
  ~MatchSeqTestNotInSet();
  bool seqIsOk(ScoredSeq* seq);
 private:
  set<ScoredSeq*> _seqSet;
};


// imposes multiple tests
class MatchSeqTestMulti : public MatchSeqTest {
 public:
  MatchSeqTestMulti();
  MatchSeqTestMulti(MatchSeqTest** testArray, int numTests);
  ~MatchSeqTestMulti();
  bool seqIsOk(ScoredSeq* seq);
 private:
  MatchSeqTest** _testArray;
  int _numTests;
};


#endif


