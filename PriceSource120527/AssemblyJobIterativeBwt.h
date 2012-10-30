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

/* lower fidelity, higher speed */


#ifndef ASSEMBLYJOBITERATIVEBWT_H
#define ASSEMBLYJOBITERATIVEBWT_H

#include <set>
#include <map>
#include <queue>
#include "AssemblyJob.h"
#include "Alignment.h"
#include "ScoredSeq.h"
#include "ScoredSeqNormalized.h"
#include "ScoredSeqFlip.h"
#include "DynamicProgrammingAligner.h"
#include "ScoredSeqCollectionBwt.h"
#include "MatchSeqTest.h"
using namespace::std;

class AssemblyJobIterativeBwt: public AssemblyJob {

 public:
  AssemblyJobIterativeBwt();
  // HERE, NON-THREADED IS DEFAULT
  AssemblyJobIterativeBwt(set<ScoredSeq*>* seqs,
			  AssemblyJob::AssemblyStrandedness strandedness,
			  AssemblyJob::GapStatus gapStatus,
			  DynamicProgrammingAligner * als,
			  AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);
  AssemblyJobIterativeBwt(set<ScoredSeq*>* seqs,
			  AssemblyJob::AssemblyStrandedness strandedness,
			  AssemblyJob::GapStatus gapStatus,
			  DynamicProgrammingAligner * als,
			  long maxOverlap,
			  AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);

  // deletes the DynamicProgrammingAligner that was input
  ~AssemblyJobIterativeBwt();

  // indicates whether the original input min ovl should always be used
  // (static) or if it should be changed to match the query size (dynamic)
  enum MinOvlMode { STATICMINOVL=0, DYNAMICMINOVL=1 };

  void OK();

  /* see specs from parent class */
  void makeShallow(bool shallowDelete);
  void runJob(AssemblyJob::AssemblyThreadedness threadedness);
  bool wasRun();

  void allInputSeqs( set<ScoredSeq*>* );
  void remainingInputSeqs( set<ScoredSeq*>* );
  void discardedInputSeqs( set<ScoredSeq*>* );
  void newlyAssembledSeqs( set<ScoredSeq*>* );
  void allAssembledSeqs( set<ScoredSeq*>* );

 private:
  bool _gapped;
  bool _usingMaxOverlap;

  void constructorHelper(set<ScoredSeq*>* inputSeqs);
  //void subAssembly(AssemblyJob::AssemblyThreadedness threadedness);
  DynamicProgrammingAligner * _als;
  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _retainedSeqs;
  set<ScoredSeq*> _novelSeqs;
  set<ScoredSeq*> _usedSeqs;
  set<ScoredSeq*> _tooShortSeqs; // special for this class; do not participate in the assembly

  bool _wasRun;
  bool _inputRemoved;
  bool _novelRemoved;
  AssemblyJob::AssemblyStrandedness _strandedness;
  long _maxOverlap;
  long _minOverlap;
  AssemblyJob::AssemblyMode _assemblyMode;

  // helper methods (make threading easier)
  void runJobGetMatchesHelper(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, ScoredSeqCollectionBwt* seqCollectionBwt,
			      MatchSeqTest** seqTestArray, AssemblyThreadedness threadedness);
};

#endif
