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

/* packages all of the parameters that are needed for generating/scoring
 * the alignments that are generated during an assembly process.  this
 * class simplifies argument lists since it can be cast as any of the 
 * three Params types whose info it carries.
 */


#ifndef PARAMSALIGNMENT_H
#define PARAMSALIGNMENT_H

#include "AlignmentScoreMatrix.h"
#include "ParamsMinOverlap.h"
#include "ParamsMinFractId.h"

class ParamsAlignment : public ParamsMinOverlap, public ParamsMinFractId {

 public:
  ParamsAlignment();
  ParamsAlignment(ParamsMinOverlap* pmo, ParamsMinFractId* pmf, AlignmentScoreMatrix * sm);
  ParamsAlignment(ParamsAlignment * pa);
  ~ParamsAlignment();

  // implementations of the methods from the three parent classes
  long globalMinOverlap();
  long minOvlThreshold();
  float minOvlSlope();
  long altMinOvlBaseline();
  long calculateMinOverlap(long readNum);

  float globalFractId();
  long fractIdThreshold();
  long denominatorLogStep();
  // the actual calculation
  float calculateMinFractId(long readNum);
  float calculateMinFractId(set<ScoredSeq*>* reads);
  float calculateMinFractId(set<ScoredSeq*>* reads, float readCountFactor);



  AlignmentScoreMatrix* getAsm();

 private:

  // the three represented parent classes
  AlignmentScoreMatrix * _alScoreMatrix;
  ParamsMinOverlap * _paramsMinOverlap;
  ParamsMinFractId * _paramsMinFractId;

};

#endif
