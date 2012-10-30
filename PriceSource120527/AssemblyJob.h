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

/* A modularized assembly job, involving a set of ScoredSeqs and
 * an AlignmentStrategy (i.e. a set of allowed criteria for
 * combining seqs into contigs).
 */


#ifndef ASSEMBLYJOB_H
#define ASSEMBLYJOB_H

#include <set>
#include <map>
class ScoredSeq;
using namespace::std;

class AssemblyJob{

 public:
  virtual ~AssemblyJob();

  // used to make the assembly jobs thread-safe; must be run before runJob
  virtual void makeShallow(bool shallowDelete) = 0;

  /* access sequences (or subsets of sequences)
   * REQUIRES: wasRun() -> true (except allInputSeqs)
   * REQUIRES: those sequences have not been removed */
  virtual void allInputSeqs( set<ScoredSeq*>* ) = 0;
  virtual void remainingInputSeqs( set<ScoredSeq*>* ) = 0;
  virtual void discardedInputSeqs( set<ScoredSeq*>* ) = 0;
  virtual void newlyAssembledSeqs( set<ScoredSeq*>* ) = 0;
  virtual void allAssembledSeqs( set<ScoredSeq*>* ) = 0;

  virtual void OK() = 0;

  // used to determine AssemblyJob strandedness
  enum AssemblyStrandedness{ SINGLESTRANDED=0, DOUBLESTRANDED=1 };
  enum GapStatus{ UNGAPPED=0, GAPPED=1 };
  enum AssemblyMode{ FULLASSEMBLY=0, REDUNDANTASSEMBLY=1 };
  enum AssemblyThreadedness{ NOTTHREADED=0, THREADED=1 };

  // doesn't re-run job, just runs it the first time
  virtual void runJob(AssemblyThreadedness threadedness) = 0;
  // false unless job is complete
  virtual bool wasRun() = 0;
};

#endif
