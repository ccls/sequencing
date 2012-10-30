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

#ifndef MATCHSEQTEST_CPP
#define MATCHSEQTEST_CPP

#include "MatchSeqTest.h"


MatchSeqTest::~MatchSeqTest(){}
MatchSeqTest* MatchSeqTest::getNullTest(){ return new MatchSeqTestNull(); }


MatchSeqTestNull::MatchSeqTestNull(){}
MatchSeqTestNull::~MatchSeqTestNull(){}
bool MatchSeqTestNull::seqIsOk(ScoredSeq* seq){ return true; }


MatchSeqTestMinLength::MatchSeqTestMinLength(long minLength) : _minLength(minLength){}
MatchSeqTestMinLength::~MatchSeqTestMinLength(){}
bool MatchSeqTestMinLength::seqIsOk(ScoredSeq* seq){ return seq->size() >= _minLength; }


MatchSeqTestExactLength::MatchSeqTestExactLength(long length) : _length(length){}
MatchSeqTestExactLength::~MatchSeqTestExactLength(){}
bool MatchSeqTestExactLength::seqIsOk(ScoredSeq* seq){ return seq->size() == _length; }


MatchSeqTestInSet::MatchSeqTestInSet(set<ScoredSeq*>* seqSet){
  _seqSet.insert(seqSet->begin(), seqSet->end());
}
MatchSeqTestInSet::~MatchSeqTestInSet(){}
bool MatchSeqTestInSet::seqIsOk(ScoredSeq* seq){ return _seqSet.count(seq) > 0; }



MatchSeqTestNotInSet::MatchSeqTestNotInSet(){}
MatchSeqTestNotInSet::MatchSeqTestNotInSet(set<ScoredSeq*>* seqSet){
  _seqSet.insert(seqSet->begin(), seqSet->end());
}
MatchSeqTestNotInSet::~MatchSeqTestNotInSet(){}
bool MatchSeqTestNotInSet::seqIsOk(ScoredSeq* seq){ return _seqSet.count(seq) == 0; }



MatchSeqTestMulti::MatchSeqTestMulti(){}
MatchSeqTestMulti::MatchSeqTestMulti(MatchSeqTest** testArray, int numTests) : _numTests(numTests){
  _testArray = new MatchSeqTest*[ _numTests ];
  for (int n = 0; n < _numTests; ++n){ _testArray[n] = testArray[n]; }
}
MatchSeqTestMulti::~MatchSeqTestMulti(){
  delete [] _testArray;
}
bool MatchSeqTestMulti::seqIsOk(ScoredSeq* seq){
  int testNum = 0;
  while ( testNum < _numTests and _testArray[testNum]->seqIsOk(seq) ){ ++testNum; }
  return testNum == _numTests;
}


#endif

