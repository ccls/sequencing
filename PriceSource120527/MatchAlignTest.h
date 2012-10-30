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

#ifndef MATCHALIGNTEST_H
#define MATCHALIGNTEST_H

#include <set>
#include "Alignment.h"
using namespace::std;

// this is a pure virtual class for tests that can be passed to
// AlignmentCollectionBwt to evaluate aluences' appropriateness for
// return based on properties other than the potential alignments they
// might form

// these classes are going to accelerate the assembler by performing tests
// on potential match aluences prior to generating an alignment rather than
// waiting until afterwards.  provided implementations are below.
class MatchAlignTest {
 public:
  static MatchAlignTest* getNullTest();
  virtual ~MatchAlignTest();
  virtual bool alignIsOk(Alignment* al) = 0;
};


// i am including the implementing classes in this file

// this one always returns true
class MatchAlignTestNull : public MatchAlignTest {
 public:
  MatchAlignTestNull();
  ~MatchAlignTestNull();
  bool alignIsOk(Alignment* al);
};

// makes sure that seqA does not overhang seqB
class MatchAlignTestNoOverhangA : public MatchAlignTest {
 public:
  MatchAlignTestNoOverhangA();
  ~MatchAlignTestNoOverhangA();
  bool alignIsOk(Alignment* al);
};

// imposes multiple tests
class MatchAlignTestMulti : public MatchAlignTest {
 public:
  MatchAlignTestMulti();
  MatchAlignTestMulti(MatchAlignTest** testArray, int numTests);
  ~MatchAlignTestMulti();
  bool alignIsOk(Alignment* al);
 private:
  MatchAlignTest** _testArray;
  int _numTests;
};


#endif


