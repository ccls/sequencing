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
 * efficiency over accuracy/completeness.  This class is for
 * combining sequences that are perfect subsets of other sequences,
 * only for those sequences that are below a threshold length.
 */


#ifndef ASSEMBLYJOBSUBSET_H
#define ASSEMBLYJOBSUBSET_H

#include <set>
#include <map>
#include <queue>
#include "ScoredSeq.h"
#include "ScoredSeqNormalized.h"
#include "ScoredSeqCollectionBwt.h"
#include "Alignment.h"
#include "AssemblyJob.h"
using namespace::std;

class AssemblyJobSubset: public AssemblyJob {

 public:
  AssemblyJobSubset();
  AssemblyJobSubset(set<ScoredSeq*>* seqs,
		    float minFractId,
		    AssemblyJob::AssemblyStrandedness strandedness);


  // deletes the AlignmentStrategy that was input
  ~AssemblyJobSubset();

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

  // this is only public so that I can run tests on it
  ScoredSeqNormalized* combineHelper(vector<Alignment*>* alSet, ScoredSeqNormalized* oldSeq);

 private:

  // temporary, for troubleshooting
  void printAlignNw(Alignment * al);

  struct SortSsnByLength {
    bool operator() (ScoredSeqNormalized* ssnA, ScoredSeqNormalized* ssnB);
  };

  set<ScoredSeq*> _inputSeqs;
  set<ScoredSeq*> _remainSeqs;
  set<ScoredSeq*> _usedSeqs;
  set<ScoredSeq*> _novelSeqs;

  set<ScoredSeq*> _allAssembledSeqs;

  bool _wasRun;
  bool _inputRemoved;
  bool _novelRemoved;
  AssemblyJob::AssemblyStrandedness _strandedness;

  float _minFractId;

  void runJobMapping(AssemblyJob::AssemblyThreadedness threadedness);
  void runJobCombine();

  // these will be created in runJobMapping and used/destroyed in runJobCombine
  long _numInput;
  ScoredSeqNormalized** _shortSeqArray;
  vector<Alignment*>** _shortIndexedAlSets;
  // indexed from zero; these indexes are in alignment with the order in _initialLongs
  long* _longToShortIndexes;
  vector<Alignment*>** _longAlignArray;
  ScoredSeqNormalized** _initialLongs;
  long _longCount;

  // all sequences that will be combined exist as keys in at least one
  // of these two maps as keys.
  //map<ScoredSeqNormalized*, vector<Alignment*>* > _shortToAlignments; // keys are seqA
  // entries of this table will be removed as longs are combined
  //map<ScoredSeqNormalized*, vector<Alignment*>* > _longToAlignments; // keys are seqB
  // this list will be sorted by length, short->long
  // don't delete these, they are found in the set below
  //vector<ScoredSeqNormalized*> _initialLongs;
  // this set will keep track of which normalized sequences were
  // part of the input and thus don't need to be separately dealt with.
  //vector<ScoredSeqNormalized*> _initialCarriers; // for non-redundant outer deleteion
  /*
  // things that are created as intermediates
  set<ScoredSeqNormalized*> _transientCarriers; // for non-redundant deep deleteion
  */
};

#endif
