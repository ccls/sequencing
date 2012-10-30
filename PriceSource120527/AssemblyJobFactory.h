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

/* the front-end for generating assemblies 
*/


#ifndef ASSEMBLYJOBFACTORY_H
#define ASSEMBLYJOBFACTORY_H

#include <set>
using namespace::std;

#include "ScoredSeq.h"
#include "AssemblyJob.h"
#include "AlignmentScoreMatrix.h"
#include "ParamsMinOverlap.h"
#include "ParamsMinFractId.h"
#include "ParamsAlignment.h"
#include "ParamsDeBruijn.h"
#include "AssemblyJobNull.h"
#include "AssemblyJobHierarchy.h"
#include "AssemblerListener.h"

class AssemblyJobFactory {

 public:
  AssemblyJobFactory();
  AssemblyJobFactory(ParamsAlignment* pa,
		     ParamsDeBruijn* pdb,
		     AssemblyJob::AssemblyStrandedness strandedness,
		     AssemblerListener* listener);
  AssemblyJobFactory(ParamsAlignment* pa,
		     ParamsDeBruijn* pdb,
		     long maxOverlap,
		     AssemblyJob::AssemblyStrandedness strandedness,
		     AssemblerListener* listener);
  AssemblyJobFactory(ParamsMinOverlap * pmo,
		     ParamsMinFractId * pmf,
		     AlignmentScoreMatrix * alScoreMatrix,
		     ParamsDeBruijn* pdb,
		     AssemblyJob::AssemblyStrandedness strandedness,
		     AssemblerListener* listener);
  AssemblyJobFactory(ParamsMinOverlap * pmo,
		     ParamsMinFractId * pmf,
		     AlignmentScoreMatrix * alScoreMatrix,
		     ParamsDeBruijn* pdb,
		     long maxOverlap,
		     AssemblyJob::AssemblyStrandedness strandedness,
		     AssemblerListener* listener);

  // for replicating the factory, maybe with a change
  AssemblyJobFactory(AssemblyJobFactory* ajf);
  AssemblyJobFactory(AssemblyJobFactory* ajf, ParamsMinOverlap* pmo);
  AssemblyJobFactory(AssemblyJobFactory* ajf, ParamsMinFractId* pmf);
  AssemblyJobFactory(AssemblyJobFactory* ajf, AlignmentScoreMatrix* alScoreMatrix);
  AssemblyJobFactory(AssemblyJobFactory* ajf, AssemblyJob::AssemblyStrandedness strandedness);
  AssemblyJobFactory(AssemblyJobFactory* ajf, long maxOverlap);

  ~AssemblyJobFactory();

  // get back the input parameters
  AlignmentScoreMatrix * getScoreMatrix();
  ParamsMinOverlap * getParamsMinOverlap();
  ParamsMinFractId * getParamsMinFractId();
  ParamsDeBruijn * getParamsDeBruijn();
  long getMaxOverlap();
  AssemblyJob::AssemblyStrandedness strandedness();
  AssemblerListener* getListener();

  // only assembles using sense (i.e. (+)strand) overlaps of the input sequences
  // OR
  // input sequences can be assembled using both (+) and (-) -orientation alignments
  // DEPENDING ON OBJECT FIELD
  AssemblyJob * assemblyJob( set<ScoredSeq*>* inputSeqs);
  AssemblyJob * assemblyJob( set<ScoredSeq*>* inputSeqs, float contigFactor);

  // only removes highly redundant sequences
  AssemblyJob * redundancyJob( set<ScoredSeq*>* inputSeqs);
  AssemblyJob * redundancyJob( set<ScoredSeq*>* inputSeqs, float contigFactor);


 private:
  // helper methods to consolidate the implementations of the two assembly job methods above
  AssemblyJob * generalAssemblyJob(set<ScoredSeq*>* inputSeqs, long inputSeqNum);
  AssemblyJob * generalRedundancyJob(set<ScoredSeq*>* inputSeqs, long inputSeqNum);

  // can be assigned to the default or an input
  AlignmentScoreMatrix * _alScoreMatrix;
  ParamsMinOverlap * _paramsMinOverlap;
  ParamsMinFractId * _paramsMinFractId;
  ParamsDeBruijn* _paramsDeBruijn;
  long _maxOverlap;
  AssemblyJob::AssemblyStrandedness _strandedness;
  AssemblerListener* _listener;

};

#endif
