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

#ifndef MATCHOFFSETTEST_CPP
#define MATCHOFFSETTEST_CPP

#include "MatchOffsetTest.h"


MatchOffsetTest::~MatchOffsetTest(){}


MatchOffsetTestNull::MatchOffsetTestNull(){}
MatchOffsetTestNull::~MatchOffsetTestNull(){}
bool MatchOffsetTestNull::offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB){ return true; }


MatchOffsetTestFullOverlapA::MatchOffsetTestFullOverlapA(){}
MatchOffsetTestFullOverlapA::~MatchOffsetTestFullOverlapA(){}
bool MatchOffsetTestFullOverlapA::offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB){
  return offset >= 0 and seqA->size() <= seqB->size() and seqB->size() + offset <= seqA->size();
}



MatchOffsetTestMinOverlap::MatchOffsetTestMinOverlap(){}
MatchOffsetTestMinOverlap::MatchOffsetTestMinOverlap(long minOverlap) : _minOverlap(minOverlap) {}
MatchOffsetTestMinOverlap::~MatchOffsetTestMinOverlap(){}
bool MatchOffsetTestMinOverlap::offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB){
  // if the first bool is in danger of being violated, then the offset will be positive
  // if the second bool is in danger of being violated, then the offset will be negative
  return seqA->size() - offset >= _minOverlap and seqB->size() + offset >= _minOverlap;
}



MatchOffsetTestMulti::MatchOffsetTestMulti(){}
MatchOffsetTestMulti::MatchOffsetTestMulti(MatchOffsetTest** testArray, int numTests) : _numTests(numTests){
  _testArray = new MatchOffsetTest*[ _numTests ];
  for (int n = 0; n < _numTests; ++n){ _testArray[n] = testArray[n]; }
}
MatchOffsetTestMulti::~MatchOffsetTestMulti(){
  delete [] _testArray;
}
bool MatchOffsetTestMulti::offsetIsOk(long offset, ScoredSeq* seqA, ScoredSeq* seqB){
  int testNum = 0;
  while ( testNum < _numTests and _testArray[testNum]->offsetIsOk(offset,seqA,seqB) ){ ++testNum; }
  return testNum == _numTests;
}


#endif

