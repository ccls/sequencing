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




#ifndef ASSEMBLYJOBHIERARCHY_H
#define ASSEMBLYJOBHIERARCHY_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "AssemblyJob.h"
#include "DynamicProgrammingAligner.h"
#include "ParamsMinOverlap.h"
#include "ParamsMinFractId.h"
#include "AssemblerListener.h"
#include "ParamsDeBruijn.h"
using namespace::std;

class AssemblyJobHierarchy: public AssemblyJob {

 public:
  AssemblyJobHierarchy();
  AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
		       ParamsMinOverlap * pmo,
		       ParamsMinFractId * pmf,
		       AlignmentScoreMatrix * sm,
		       ParamsDeBruijn* pmdb,
		       AssemblyJob::AssemblyStrandedness strandedness,
		       AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);
  AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
		       ParamsMinOverlap * pmo,
		       ParamsMinFractId * pmf,
		       AlignmentScoreMatrix * sm,
		       ParamsDeBruijn* pmdb,
		       AssemblyJob::AssemblyStrandedness strandedness,
		       AssemblerListener* listener,
		       AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);

  // multiplies the true number of input reads by this factor
  // for calculating assembly requirements; above, default is 1.0
  AssemblyJobHierarchy(float readCountFactor,
		       set<ScoredSeq*>* seqs,
		       ParamsMinOverlap * pmo,
		       ParamsMinFractId * pmf,
		       AlignmentScoreMatrix * sm,
		       ParamsDeBruijn* pmdb,
		       AssemblyJob::AssemblyStrandedness strandedness,
		       AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);

  AssemblyJobHierarchy(set<ScoredSeq*>* seqs,
		       ParamsMinOverlap * pmo,
		       ParamsMinFractId * pmf,
		       AlignmentScoreMatrix * sm,
		       ParamsDeBruijn* pmdb,
		       AssemblyJob::AssemblyStrandedness strandedness,
		       long maxOverlap,
		       AssemblerListener* listener,
		       AssemblyJob::AssemblyMode assemblyMode = AssemblyJob::FULLASSEMBLY);




  // deletes the DynamicProgrammingAligner that was input
  ~AssemblyJobHierarchy();

  /* see specs from parent class */
  void makeShallow(bool shallowDelete);
  void runJob(AssemblyJob::AssemblyThreadedness threadedness);
  bool wasRun();

  void allInputSeqs( set<ScoredSeq*>* );
  void remainingInputSeqs( set<ScoredSeq*>* );
  void discardedInputSeqs( set<ScoredSeq*>* );
  void newlyAssembledSeqs( set<ScoredSeq*>* );
  void allAssembledSeqs( set<ScoredSeq*>* );

  void OK();

 private:

  // normal or redundancy-only
  AssemblyJob::AssemblyMode _assemblyMode;
  long getN50(set<ScoredSeq*>* seqs);
  ParamsDeBruijn* _paramsDeBruijn;

  // a helper for runJob to shorten the code
  void executeAndCleanUpJob(AssemblyJob* aj, AssemblyJob::AssemblyThreadedness threadedness);
  void runJobHelpPerfect(AssemblyJob::AssemblyThreadedness threadedness);
  void runJobHelpDeBruijn(AssemblyJob::AssemblyThreadedness threadedness);
  void runJobHelpIterativeBwt(AssemblyJob::GapStatus gapStatus, AssemblyJob::AssemblyThreadedness threadedness);

  // decides whether or not another cycle of assembly will be performed by
  // setting a threshold for the fraction of sequences that must have been
  // collapsed during the prior threshold, with the prior set of overlap/ID
  // requirements.  this value is actually the maximum fraction of reads that
  // are allowed to be left over versus the beginning of the last round.
  float _maxRepeatFract; // default is 0.9

  // to maintain the edge-ness
  bool _usingMaxOverlap;
  long _maxOverlap;

  // these must be deleted because
  // they are copies of the input variables
  ParamsMinOverlap* _pmo;
  ParamsMinFractId* _pmf;
  float _readCountFactor;

  // these should be empty at the beginning+end of runJob()
  set<ScoredSeq*> _seqCarrier;
  set<ScoredSeq*> _localUsedSeqs;

  AlignmentScoreMatrix * _sm;
  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _allAssembledSeqs;
  set<ScoredSeq*> _novelSeqs;
  set<ScoredSeq*> _remainingSeqs;
  set<ScoredSeq*> _usedSeqs;
  bool _wasRun;
  bool _inputRemoved;
  bool _novelRemoved;
  AssemblyJob::AssemblyStrandedness _strandedness;

  // this is so that i can choose to give the object a listener
  // or not.  the first variable will always be bound to a null
  // listener that is created by the constructor and deleted by
  // the destructor.  if no listener is provided, the second 
  // variable will also be bound to this temporary object, but
  // alternatively a provided listener will be bound to that
  // variable.  whichever listener is bound to the second variable
  // will be used, and delete will not be called on that variable
  AssemblerListener* _nullListener;
  AssemblerListener* _listener;

};

#endif
