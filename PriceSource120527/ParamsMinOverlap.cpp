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

#ifndef PARAMSMINOVERLAP_CPP
#define PARAMSMINOVERLAP_CPP

#include "ParamsMinOverlap.h"
#include <cmath>

// constructors
ParamsMinOverlap::ParamsMinOverlap(){}
ParamsMinOverlap::ParamsMinOverlap(long globalMinOverlap, long minOvlThreshold, float minOvlSlope) :
  _globalMinOverlap(globalMinOverlap),
  _minOvlThreshold(minOvlThreshold),
  _minOvlSlope(minOvlSlope),
  _altBaseline(0){
}
ParamsMinOverlap::ParamsMinOverlap(long globalMinOverlap, long minOvlThreshold, float minOvlSlope, long altBaseline) :
  _globalMinOverlap(globalMinOverlap),
  _minOvlThreshold(minOvlThreshold),
  _minOvlSlope(minOvlSlope),
  _altBaseline(altBaseline){
}
ParamsMinOverlap::ParamsMinOverlap(ParamsMinOverlap * pmo){
  _globalMinOverlap = pmo->globalMinOverlap();
  _minOvlThreshold = pmo->minOvlThreshold();
  _minOvlSlope = pmo->minOvlSlope();
  _altBaseline = pmo->altMinOvlBaseline();
}
ParamsMinOverlap::ParamsMinOverlap(ParamsMinOverlap * pmo, long altBaseline) :
  _altBaseline(altBaseline){
  _globalMinOverlap = pmo->globalMinOverlap();
  _minOvlThreshold = pmo->minOvlThreshold();
  _minOvlSlope = pmo->minOvlSlope();
}

// get the values
long ParamsMinOverlap::globalMinOverlap(){ return _globalMinOverlap; }
long ParamsMinOverlap::minOvlThreshold(){ return _minOvlThreshold; }
float ParamsMinOverlap::minOvlSlope(){ return _minOvlSlope; }
long ParamsMinOverlap::altMinOvlBaseline(){ return _altBaseline; }

// calculate the real value
long ParamsMinOverlap::calculateMinOverlap(long readNum){
  long tempMin;
  if ( readNum <= _minOvlThreshold ){ tempMin = _globalMinOverlap; }
  else { tempMin = long(log10( readNum - _minOvlThreshold ) * _minOvlSlope) + _globalMinOverlap; }
  if (tempMin < _altBaseline){ tempMin = _altBaseline; }
  return tempMin;
}


#endif
