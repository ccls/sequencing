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

/* packages paramaters that will be used by AssemblyJobFactory
 * to calculate an appropriate minOverlap.  it is super simple,
 * and exists just to consolidate a series of values and make
 * them clearly distinct from parameters for other things (like,
 * say, min fract ID).  there are no default values here; those are
 * determined by the AssemblyJobFactory.
 */


#ifndef PARAMSMINOVERLAP_H
#define PARAMSMINOVERLAP_H


class ParamsMinOverlap {

 public:
  ParamsMinOverlap();
  // minOvlSlope is the num of nucs that will be added per order-of-magnitude increase above the minOvlThreshold
  ParamsMinOverlap(long globalMinOverlap, long minOvlThreshold, float minOvlSlope);
  ParamsMinOverlap(ParamsMinOverlap * pmo);
  // altBaseline allows scaling to remain the same but with a new bottom min overlap
  ParamsMinOverlap(long globalMinOverlap, long minOvlThreshold, float minOvlSlope, long altBaseline);
  ParamsMinOverlap(ParamsMinOverlap * pmo, long altBaseline);

  // get the values
  virtual long globalMinOverlap();
  virtual long minOvlThreshold();
  virtual float minOvlSlope();
  virtual long altMinOvlBaseline(); // 0 if not specified

  // the actual calculation
  virtual long calculateMinOverlap(long readNum);

 private:

  // three variables for minFractId
  long _globalMinOverlap;
  long _minOvlThreshold;
  float _minOvlSlope; // the number of nucs that the minOvl will increse per log10 increase above the minOvlThreshold
  long _altBaseline;
};

#endif
