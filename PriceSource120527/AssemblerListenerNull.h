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

/* this implementation prints out nothing, thereby preventing
 * stdout.
 */

#ifndef ASSEMBLERLISTENERNULL_H
#define ASSEMBLERLISTENERNULL_H

#include "AssemblerListener.h"
#include "ScoredSeq.h"
#include <string>
using namespace::std;


class AssemblerListenerNull : public AssemblerListener {

 public:

  AssemblerListenerNull();
  ~AssemblerListenerNull();

  void showParameters();
  void assemblyJobBeginStats(AssemblyJob* aj);
  void assemblyJobEndStats(AssemblyJob* aj);

  void beginAssembly(string message);
  void endAssembly(string message);

  void beginCycle(int cycleNum);
  void writeOutfile(string outfileName);
  void endCycle(set<ScoredSeq*>* contigs);

  void collectedInitialContigs(set<ScoredSeq*>* contigs);
  void filteredInitialContigs(set<ScoredSeq*>* contigs);
  void beginAddInitialContigs(set<ScoredSeq*>* contigs);
  void endAddInitialContigs(set<ScoredSeq*>* contigs);

  void beginMapping(long totalNumReads);
  void updateMapping(long addReadsToTally);
  void endMapping(long readsFiltered);

  void beginEdgeAssembly(long totalNumJobs);
  void updateEdgeAssembly(long addJobsToTally);
  void endEdgeAssembly();

  void beginMetaAssembly(long totalNumContigs);
  void endMetaAssembly(long finalNumContigs);

  void beginRepeatDetection();
  void endRepeatDetection(long numRepeats, long totalLength);

  void message(string message);
  void verboseMessage(string message);
  void errorMessage(string message);

  void timeStamp();
};

#endif
