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
 * to calculate an appropriate minFractId.  it is super simple,
 * and exists just to consolidate a series of values and make
 * them clearly distinct from parameters for other things (like,
 * say, min overlap).  there are no default values here; those are
 * determined by the AssemblyJobFactory.
 */


#ifndef PARAMSMINFRACTID_H
#define PARAMSMINFRACTID_H

#include <set>
#include "ScoredSeq.h"
using namespace::std;

class ParamsMinFractId {

 public:
  ParamsMinFractId();
  ParamsMinFractId(float globalFractId, long fractIdThreshold, long denominatorLogStep);
  ParamsMinFractId(ParamsMinFractId * pmf);

  static ParamsMinFractId* getDefault();

  // get the values
  virtual float globalFractId();
  virtual long fractIdThreshold();
  virtual long denominatorLogStep();

  // the actual calculation
  virtual float calculateMinFractId(long readNum);
  virtual float calculateMinFractId(set<ScoredSeq*>* reads);
  virtual float calculateMinFractId(set<ScoredSeq*>* reads, float readCountFactor);

 private:
  // a helper for the constructor
  void calculateDerivedValues();

  // three variables for minFractId
  float _globalFractId;
  long _fractIdThreshold;
  long _denominatorLogStep; // same as fidLogBase;
  // defines the log increment along which the fraction mismatch allowed will decrement

  // derived values:
  float _globalFractMis;
  float _fidLogConst;
  float _thresholdLogConst;

};

#endif
