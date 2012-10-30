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


/* this is for when there are few/no input seqs to be assembled, so
 * the input should just be returned.
 */

#ifndef ASSEMBLYJOBNULL_H
#define ASSEMBLYJOBNULL_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "AssemblyJob.h"
using namespace::std;

class AssemblyJobNull: public AssemblyJob {

 public:
  AssemblyJobNull();
  AssemblyJobNull(set<ScoredSeq*>* seqs);
  ~AssemblyJobNull();

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
  set<ScoredSeq*> _inputSeqs;
  bool _inputRemoved;
  // even though nothing happens when this job is run,
  // having a consistent definition of wasRun with respect
  // to whether runJob was called across all classes is
  // important for maintaining consistent expectations of
  // behavior and finding errors.
  bool _wasRun;

};

#endif
