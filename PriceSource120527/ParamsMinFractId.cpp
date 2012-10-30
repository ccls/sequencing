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

#ifndef PARAMSMINFRACTID_CPP
#define PARAMSMINFRACTID_CPP

#include "ParamsMinFractId.h"
#include "AssemblyException.h"
#include <cmath>
#include <iostream>
using namespace::std;


ParamsMinFractId* ParamsMinFractId::getDefault(){ return new ParamsMinFractId(0.85,20,25); }

// constructors
ParamsMinFractId::ParamsMinFractId(){}
ParamsMinFractId::ParamsMinFractId(float globalFractId, long fractIdThreshold, long denominatorLogStep) :
  _globalFractId(globalFractId),
  _fractIdThreshold(fractIdThreshold),
  _denominatorLogStep(denominatorLogStep){
  if (_globalFractId < 0 or _globalFractId > 1){
    throw AssemblyException::ArgError("for minFractId, globalFractId must be between (including) 0 and 1.");
  }
  calculateDerivedValues();
}
ParamsMinFractId::ParamsMinFractId(ParamsMinFractId * pmf){
  _globalFractId = pmf->globalFractId();
  _fractIdThreshold = pmf->fractIdThreshold();
  _denominatorLogStep = pmf->denominatorLogStep();
  if (_globalFractId < 0 or _globalFractId > 1){
    throw AssemblyException::ArgError("for minFractId, globalFractId must be between (including) 0 and 1.");
  }
  calculateDerivedValues();
}

void ParamsMinFractId::calculateDerivedValues(){
  // derived
  _globalFractMis = 1.0 - _globalFractId;
  _fidLogConst = log10( float(_denominatorLogStep) );
  if (_fractIdThreshold > 1){
    _thresholdLogConst = log10( _fractIdThreshold ) / _fidLogConst;
  } else { _thresholdLogConst = 0; }
}

// get the values
float ParamsMinFractId::globalFractId(){ return _globalFractId; }
long ParamsMinFractId::fractIdThreshold(){ return _fractIdThreshold; }
long ParamsMinFractId::denominatorLogStep(){ return _denominatorLogStep; }

float ParamsMinFractId::calculateMinFractId(long readNum){
  // CHANGE THIS!!! it would be better for the denominator to define a constant decay rate (half-life)
  // OR MAYBE NOT.  as written, this will cut the min % mis quickly at first but will approach the
  // asypmtote of 100% ID slowly, which would be good.
  if ( readNum <= _fractIdThreshold ){ return _globalFractId; }
  else {
    float exponent = log10( readNum ) / _fidLogConst - _thresholdLogConst;
    return 1.0 - _globalFractMis * pow ( float(0.5), exponent );
  }
}

float ParamsMinFractId::calculateMinFractId(set<ScoredSeq*>* reads){
  return calculateMinFractId(reads, 1.0);
}

float ParamsMinFractId::calculateMinFractId(set<ScoredSeq*>* reads, float readCountFactor){
  long readNum = 0;
  //for (set<ScoredSeq*>::iterator it = reads->begin(); it != reads->end(); ++it){ readNum += (*it)->size(); }
  readNum = reads->size();
  readNum = long(float(readNum) * readCountFactor);

  if ( readNum <= _fractIdThreshold ){ return _globalFractId; }
  else {
    float exponent = log10( readNum ) / _fidLogConst - _thresholdLogConst;
    return 1.0 - _globalFractMis * pow ( float(0.5), exponent );
  }
}


#endif
