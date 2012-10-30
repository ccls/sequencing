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

/* implementations of this abstract class will monitor the
 * progress of Assembler and thereby handle things like
 * stdout, GUIs, or just take the place of the annoing cout
 * streams that I would otherwise have in place.
 */

#ifndef ASSEMBLERLISTENER_H
#define ASSEMBLERLISTENER_H
#include "ScoredSeq.h"
#include "AssemblyJob.h"
#include <string>
#include <set>
using namespace::std;


class AssemblerListener {

 public:

  virtual ~AssemblerListener();

  virtual void showParameters() = 0;
  virtual void assemblyJobBeginStats(AssemblyJob* aj) = 0;
  virtual void assemblyJobEndStats(AssemblyJob* aj) = 0;

  virtual void beginAssembly(string message) = 0;
  virtual void endAssembly(string message) = 0;

  virtual void beginCycle(int cycleNum) = 0;
  virtual void writeOutfile(string outfileName) = 0;
  virtual void endCycle(set<ScoredSeq*>* contigs) = 0;

  virtual void collectedInitialContigs(set<ScoredSeq*>* contigs) = 0;
  virtual void filteredInitialContigs(set<ScoredSeq*>* contigs) = 0;
  virtual void beginAddInitialContigs(set<ScoredSeq*>* contigs) = 0;
  virtual void endAddInitialContigs(set<ScoredSeq*>* contigs) = 0;

  virtual void beginMapping(long totalNumReads) = 0;
  virtual void updateMapping(long addReadsToTally) = 0;
  virtual void endMapping(long readsFiltered) = 0;

  virtual void beginEdgeAssembly(long totalNumJobs) = 0;
  virtual void updateEdgeAssembly(long addJobsToTally) = 0;
  virtual void endEdgeAssembly() = 0;

  virtual void beginMetaAssembly(long totalNumContigs) = 0;
  virtual void endMetaAssembly(long finalNumContigs) = 0;

  virtual void beginRepeatDetection() = 0;
  virtual void endRepeatDetection(long numRepeats, long totalLength) = 0;

  virtual void message(string message) = 0;
  virtual void verboseMessage(string message) = 0;
  virtual void errorMessage(string message) = 0;

  virtual void timeStamp() = 0;

  static long getTotalSize(set<ScoredSeq*>* contigs);
  static long getMaxSize(set<ScoredSeq*>* contigs);
  static long getN50(set<ScoredSeq*>* contigs);

  // used by other classes to make sure that the output file is
  // in a real directory - does so by creating a temporary version
  // of the file that will then be deleted (if the file does not
  // already exist).
  static void ensureWritability(const char* filename);
};

#endif
