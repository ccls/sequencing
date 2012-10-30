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



#ifndef ALIGNMENT_CPP
#define ALIGNMENT_CPP

#include "Alignment.h"
#include <limits.h>
using namespace::std;

Alignment::~Alignment(){}


float Alignment::_SCORE_MAX = float(INT_MAX / 2);



Alignment::ScoreCalculator::~ScoreCalculator(){}

Alignment::ScoreCalcAdd::ScoreCalcAdd(){}
Alignment::ScoreCalcAdd::~ScoreCalcAdd(){}


float Alignment::ScoreCalcAdd::adjustA(float a){ return a; }
float Alignment::ScoreCalcAdd::adjustB(float b){ return b; }
float Alignment::ScoreCalcAdd::aPlusB(float a, float b){
  if ( _SCORE_MAX - a < b ){ return _SCORE_MAX; }
  else { return a + b; }
}
float Alignment::ScoreCalcAdd::aMinusB(float a, float b){
  if (a < b){ return 0; }
  else { return a - b; }
}
float Alignment::ScoreCalcAdd::bMinusA(float b, float a){
  if (b < a){ return 0; }
  else { return b - a; }
}



Alignment::ScoreCalcNormalizeA::ScoreCalcNormalizeA(){};
Alignment::ScoreCalcNormalizeA::~ScoreCalcNormalizeA(){};

Alignment::ScoreCalcNormalizeA::ScoreCalcNormalizeA(int denominator) :
  _denominator(denominator){}
float Alignment::ScoreCalcNormalizeA::adjustA(float a){ return a / _denominator; }
float Alignment::ScoreCalcNormalizeA::adjustB(float b){ return b; }
float Alignment::ScoreCalcNormalizeA::aPlusB(float a, float b){
  a = adjustA(a);
  if ( _SCORE_MAX - a < b ){ return _SCORE_MAX; }
  else { return a + b; }
}
float Alignment::ScoreCalcNormalizeA::aMinusB(float a, float b){
  a = adjustA(a);
  if (a < b){ return 0; }
  else { return a - b; }
}
float Alignment::ScoreCalcNormalizeA::bMinusA(float b, float a){
  a = adjustA(a);
  if (b < a){ return 0; }
  else { return b - a; }
}


Alignment::ScoreCalcXy::~ScoreCalcXy(){};
Alignment::ScoreCalcXy* Alignment::ScoreCalcXy::getScoreCalcXy(ScoreCalculator* calc, bool xIsA){
  if (xIsA){ return new ScoreCalcXyA(calc); }
  else { return new ScoreCalcXyB(calc); }
}

Alignment::ScoreCalcXyA::ScoreCalcXyA(ScoreCalculator* calc) :  _calc(calc){}
Alignment::ScoreCalcXyB::ScoreCalcXyB(ScoreCalculator* calc) :  _calc(calc){}
Alignment::ScoreCalcXyA::~ScoreCalcXyA(){};
Alignment::ScoreCalcXyB::~ScoreCalcXyB(){};

float Alignment::ScoreCalcXyA::adjustX(float x){ return _calc->adjustA(x); }
float Alignment::ScoreCalcXyB::adjustX(float x){ return _calc->adjustB(x); }

float Alignment::ScoreCalcXyA::adjustY(float y){ return _calc->adjustB(y); }
float Alignment::ScoreCalcXyB::adjustY(float y){ return _calc->adjustA(y); }

float Alignment::ScoreCalcXyA::xPlusY(float x, float y){ return _calc->aPlusB(x,y); }
float Alignment::ScoreCalcXyB::xPlusY(float x, float y){ return _calc->aPlusB(y,x); }

float Alignment::ScoreCalcXyA::xMinusY(float x, float y){ return _calc->aMinusB(x,y); }
float Alignment::ScoreCalcXyB::xMinusY(float x, float y){ return _calc->aMinusB(y,x); }

float Alignment::ScoreCalcXyA::yMinusX(float y, float x){ return _calc->bMinusA(y,x); }
float Alignment::ScoreCalcXyB::yMinusX(float y, float x){ return _calc->bMinusA(x,y); }


#endif
