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

#ifndef PARAMSALIGNMENT_CPP
#define PARAMSALIGNMENT_CPP

#include "ParamsAlignment.h"

// constructors
ParamsAlignment::ParamsAlignment(){}
ParamsAlignment::ParamsAlignment(ParamsMinOverlap* pmo, ParamsMinFractId* pmf, AlignmentScoreMatrix * sm){
  _alScoreMatrix = new AlignmentScoreMatrix( sm );
  _paramsMinOverlap = new ParamsMinOverlap( pmo );
  _paramsMinFractId = new ParamsMinFractId( pmf );
}
ParamsAlignment::ParamsAlignment(ParamsAlignment * pa){
  _alScoreMatrix = pa->getAsm();
  _paramsMinOverlap = new ParamsMinOverlap( dynamic_cast<ParamsMinOverlap*>(pa) );
  _paramsMinFractId = new ParamsMinFractId( dynamic_cast<ParamsMinFractId*>(pa) );
}
ParamsAlignment::~ParamsAlignment(){
  delete _alScoreMatrix;
  delete _paramsMinOverlap;
  delete _paramsMinFractId;
}

// get the values
long ParamsAlignment::globalMinOverlap(){ return _paramsMinOverlap->globalMinOverlap(); }
long ParamsAlignment::minOvlThreshold(){ return _paramsMinOverlap->minOvlThreshold(); }
float ParamsAlignment::minOvlSlope(){ return _paramsMinOverlap->minOvlSlope(); }
long ParamsAlignment::calculateMinOverlap(long readNum){ return _paramsMinOverlap->calculateMinOverlap(readNum); }
long ParamsAlignment::altMinOvlBaseline(){ return _paramsMinOverlap->altMinOvlBaseline(); }

float ParamsAlignment::globalFractId(){ return _paramsMinFractId->globalFractId(); }
long ParamsAlignment::fractIdThreshold(){ return _paramsMinFractId->fractIdThreshold(); }
long ParamsAlignment::denominatorLogStep(){ return _paramsMinFractId->denominatorLogStep(); }
// the actual calculation
float ParamsAlignment::calculateMinFractId(long readNum){ return _paramsMinFractId->calculateMinFractId(readNum); }
float ParamsAlignment::calculateMinFractId(set<ScoredSeq*>* reads){ return _paramsMinFractId->calculateMinFractId(reads); }
float ParamsAlignment::calculateMinFractId(set<ScoredSeq*>* reads, float readCountFactor){ return _paramsMinFractId->calculateMinFractId(reads,readCountFactor); }

AlignmentScoreMatrix* ParamsAlignment::getAsm(){ return new AlignmentScoreMatrix(_alScoreMatrix); }

#endif
