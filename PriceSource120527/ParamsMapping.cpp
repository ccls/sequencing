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

#ifndef PARAMSMAPPING_CPP
#define PARAMSMAPPING_CPP

#include "ParamsMapping.h"
#include "AssemblyException.h"


ParamsMapping* ParamsMapping::getDefault(){ return new ParamsMapping(25000, 0.9); }

// constructors
ParamsMapping::ParamsMapping(){}
ParamsMapping::ParamsMapping(long stepSize, float fractId) :
  _stepSize(stepSize),
  _fractId(fractId) {
  // make sure that stepSizeAddendum doesn't get zeroed.
  _stepSizeAddendum = stepSize / 3;
  if (_stepSizeAddendum == 0){ _stepSizeAddendum = stepSize; }
  if (_stepSize < 1){ throw AssemblyException::ArgError("ParamsMapping: mapping must occur in steps of size > 0."); }
  if (_fractId < 0 or _fractId > 1){
    throw AssemblyException::ArgError("for mapping, fractId must be between (including) 0 and 1.");
  }
}
ParamsMapping::ParamsMapping(long stepSize, long stepSizeAddendum, float fractId) :
  _stepSize(stepSize),
  _stepSizeAddendum(stepSizeAddendum),
  _fractId(fractId) {
  if (_stepSize < 1){ throw AssemblyException::ArgError("ParamsMapping: mapping must occur in steps of size > 0."); }
  if (_fractId < 0 or _fractId > 1){
    throw AssemblyException::ArgError("for mapping, fractId must be between (including) 0 and 1.");
  }
}
ParamsMapping::ParamsMapping(ParamsMapping * pmap){
  _stepSize = pmap->stepSize();
  _stepSizeAddendum = pmap->stepSizeAddendum();
  _fractId = pmap->fractId();
  if (_fractId < 0 or _fractId > 1){
    throw AssemblyException::ArgError("for mapping, fractId must be between (including) 0 and 1.");
  }
}

// get the values
long ParamsMapping::stepSize(){ return _stepSize; }
long ParamsMapping::stepSizeAddendum(){ return _stepSizeAddendum; }

float ParamsMapping::fractId(){ return _fractId; }


#endif
