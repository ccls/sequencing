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

/* This class is not meant to perform a near-complete assembly;
 * its goal is to simply reduce the complexity of downstream
 * *real* assemblies by getting rid of sequence entries that are
 * almost but not quite identical.  It will therefore favor
 * efficiency over accuracy/completeness.
 */


#ifndef ASSEMBLYJOBNEARPERFECT_H
#define ASSEMBLYJOBNEARPERFECT_H

#include <set>
#include <map>
#include <queue>
#include "ScoredSeq.h"
#include "ScoredSeqNormalized.h"
#include "AlignmentFull.h"
#include "AssemblyJob.h"
#include "AlignmentScoreMatrix.h"
using namespace::std;

class AssemblyJobNearPerfect: public AssemblyJob {

 public:
  AssemblyJobNearPerfect();
 AssemblyJobNearPerfect(set<ScoredSeq*>* seqs,
			float minFractId,
			AssemblyJob::AssemblyStrandedness strandedness);


  // deletes the AlignmentStrategy that was input
  ~AssemblyJobNearPerfect();

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

  // just to make passing params more efficient
  class WindowStartAndLenCarrier{
  public:
    WindowStartAndLenCarrier();
    WindowStartAndLenCarrier(long seqLen, long maxMisNum);
    ~WindowStartAndLenCarrier();
    long* _winStart;
    long* _winLen;
    long _numWindows;
  };

  class SeqCarrier{
  public:
    SeqCarrier(ScoredSeq* seq, long numWindows, bool wasInput);
    ~SeqCarrier();
    ScoredSeq* _seq;
    long* _stringPos;
    long* _seqPos;
    bool _wasInput;
  };

  // THE NEW
  void subAssembly2(long seqLen, vector<ScoredSeq*>* monoLenSeqSet);


  // THE OLD
  void subAssembly(long seqLen, vector<ScoredSeq*>* monoLenSeqSet);
  void subAssembly(long seqLen, ScoredSeqNormalized** monoLenSeqSet);
  void insertSeqHelper(ScoredSeqNormalized* newContig,
		       WindowStartAndLenCarrier* wslc,
		       map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs);
  void removeSeqHelper(ScoredSeqNormalized* oldSeq,
		       WindowStartAndLenCarrier* wslc,
		       map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs);
  void getLegitMatchesHelper(ScoredSeq* inputSeq, char sense,
			     WindowStartAndLenCarrier* wslc,  long maxMisNum,
			     map<string, set<ScoredSeqNormalized*>* >* windowToSubstrToSeqs,
			     map<ScoredSeqNormalized*,AlignmentFull*>* legitMatches);
  inline long compareSeqHelper(char* seqA, char* seqB, long seqSize, char sense, long maxMisNum);

  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _usedSeqs;
  set<ScoredSeq*> _retainedSeqs;
  set<ScoredSeq*> _novelSeqs;

  bool _wasRun;
  bool _inputRemoved;
  bool _novelRemoved;
  AssemblyJob::AssemblyStrandedness _strandedness;

  float _maxFractMis; // keep track of this to make calculating the num allowed mismatches faster
  AlignmentScoreMatrix * _mismatchCountAsm;
  AlignmentScoreMatrix * _matchCountAsm;
};

#endif
