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

#ifndef MATCHALIGNTEST_CPP
#define MATCHALIGNTEST_CPP

#include "MatchAlignTest.h"


MatchAlignTest::~MatchAlignTest(){}
MatchAlignTest* MatchAlignTest::getNullTest(){ return new MatchAlignTestNull(); }


MatchAlignTestNull::MatchAlignTestNull(){}
MatchAlignTestNull::~MatchAlignTestNull(){}
bool MatchAlignTestNull::alignIsOk(Alignment* al){ return true; }


MatchAlignTestNoOverhangA::MatchAlignTestNoOverhangA(){}
MatchAlignTestNoOverhangA::~MatchAlignTestNoOverhangA(){}
bool MatchAlignTestNoOverhangA::alignIsOk(Alignment* al){
  ScoredSeq* seqA = al->seqA();
  long lastPos = seqA->size() - 1;
  return ( (! al->isGapped(0, seqA) ) and (! al->isGapped(lastPos, seqA) ) );
}



MatchAlignTestMulti::MatchAlignTestMulti(){}
MatchAlignTestMulti::MatchAlignTestMulti(MatchAlignTest** testArray, int numTests) : _numTests(numTests){
  _testArray = new MatchAlignTest*[ _numTests ];
  for (int n = 0; n < _numTests; ++n){ _testArray[n] = testArray[n]; }
}
MatchAlignTestMulti::~MatchAlignTestMulti(){
  delete [] _testArray;
}
bool MatchAlignTestMulti::alignIsOk(Alignment* al){
  int testNum = 0;
  while ( testNum < _numTests and _testArray[testNum]->alignIsOk(al) ){ ++testNum; }
  return testNum == _numTests;
}


#endif

