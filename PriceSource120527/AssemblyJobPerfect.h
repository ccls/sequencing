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




#ifndef ASSEMBLYJOBPERFECT_H
#define ASSEMBLYJOBPERFECT_H

#include <set>
#include <map>
#include <queue>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "AssemblyJob.h"
using namespace::std;

class AssemblyJobPerfect: public AssemblyJob {

 public:
  AssemblyJobPerfect();
  AssemblyJobPerfect(set<ScoredSeq*>* seqs,
		     AssemblyJob::AssemblyStrandedness strandedness);


  // deletes the DynamicProgrammingAligner that was input
  ~AssemblyJobPerfect();

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

  // REQUIRES: the vector is set-ordered
  void subAssembly(vector<ScoredSeq*>* monoLenSeqSet, set<ScoredSeq*>* retainedSeqs,
		   set<ScoredSeq*>* usedSeqs, set<ScoredSeq*>* novelSeqs);
  //ScoredSeq* subAssemblyHelper(ScoredSeq* newInput, ScoredSeq* matchSeq, char sense);

  class SameSizeSeqCombiner {
  public:
    SameSizeSeqCombiner(ScoredSeq* original);
    ~SameSizeSeqCombiner();
    ScoredSeq* getSeq();
    void addSeq(ScoredSeq* seq, char sense);
    bool isOriginal();
    ScoredSeq* getOriginal();
  private:
    ScoredSeq* _original;
    bool _isOriginal;
    // only created when _isOriginal becomes false; NOT deleted by the destructor because they
    // will be input into the seq returned by getSeq directly (rep exposed)
    float* _scores;
    float* _links;
    // to check that the above requirement happened
    bool _seqGotten;
  };

  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _usedSeqs;
  set<ScoredSeq*> _retainedSeqs;
  set<ScoredSeq*> _novelSeqs;

  set<ScoredSeq*> _newInputSet;

  set<ScoredSeq*> _allAssembledSeqs;
  set<ScoredSeq*> _remainingSeqs;
  bool _wasRun;
  bool _inputRemoved;
  bool _novelRemoved;
  AssemblyJob::AssemblyStrandedness _strandedness;

};

#endif
