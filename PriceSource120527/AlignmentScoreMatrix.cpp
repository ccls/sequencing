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


#ifndef ALIGNMENTSCOREMATRIX_CPP
#define ALIGNMENTSCOREMATRIX_CPP

#include "AlignmentScoreMatrix.h"
#include "AssemblyException.h"
#include <iostream>
using namespace::std;

AlignmentScoreMatrix* AlignmentScoreMatrix::getDefault(){ return new AlignmentScoreMatrix(1,-2,-5,-2); }

AlignmentScoreMatrix::AlignmentScoreMatrix(){}
AlignmentScoreMatrix::AlignmentScoreMatrix(long match, long mismatch, long newGap, long extendGap) :
  _match(match),
  _mismatch(mismatch),
  _newGap(newGap),
  _extendGap(extendGap)
{ }

AlignmentScoreMatrix::AlignmentScoreMatrix(AlignmentScoreMatrix * alsm){
  _match = alsm->_match;
  _mismatch = alsm->_mismatch;
  _newGap = alsm->_newGap;
  _extendGap = alsm->_extendGap;
}

/*
long AlignmentScoreMatrix::newGap(){
  //OK();
  return _newGap;
}
long AlignmentScoreMatrix::extendGap(){
  //OK();
  return _extendGap;
}
long AlignmentScoreMatrix::match(){
  //OK();
  return _match;
}
long AlignmentScoreMatrix::mismatch(){
  //OK();
  return _mismatch;
}
*/
long AlignmentScoreMatrix::match(char nucA, char nucB){
  //OK();
  if (nucA != nucB or nucA == 'N'){ return _mismatch; }
  else { return _match; }
}
bool AlignmentScoreMatrix::isMatch(char nucA, char nucB){
  //OK();
  return (nucA == nucB and nucA != 'N');
}

#endif
