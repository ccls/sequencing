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

/* packages paramaters that will be used by Assembler to decide
 * how many cycles will be run and
 * how quickly initial contigs will be added to the assembly process.
 * maxContigLinkages is the maximum number of existing contigs that
 * a read may be mapped to in order for it to be replaced by those contigs
 * in an edge assembly job.
 * NOTE: externally, cycles are indexed from 1
 */


#ifndef PARAMSCYCLEMANAGEMENT_H
#define PARAMSCYCLEMANAGEMENT_H

#include "ReadFile.h"
#include "ScoredSeq.h"
#include "ExtendCycle.h"
#include "EcoFilter.h"
#include <set>

class ParamsCycleManagement{

 public:
  ParamsCycleManagement();
  ParamsCycleManagement(int totalCycles);
  //ParamsCycleManagement(int totalCycles, long maxContigLinkages);
  ParamsCycleManagement(ParamsCycleManagement* pcm);
  ~ParamsCycleManagement();

  enum WhenInCycle{ BEFORECYCLE=0, AFTERCYCLE=1 };


  // returns a min count = 1.5 eco filter
  static EcoFilter* defaultCountEcoFilter();
  // the maximum number of existing contigs that a read may be mapped to 
  // in order for it to be replaced by those contigs in an edge assembly job.
  static long defaultMaxContigLinkages();

  // get the only immutable value
  int totalCycles();

  // get/set this value
  long maxContigLinkages();
  void setMaxContigLinkages(long mcl);

  // from EcoFilter
  void addEcoFilter(EcoFilter* sizeFilter);
  EcoFilter* getFullEcoFilter();

  // there is a separate access point for this portion of the filter since
  // it will be used in other contexts.  size filters should not be added
  // to the full filter set because they will then only be applied post-meta-assembly.

  // adds size filtering limits (note: use 0 as min size for "null")
  void addSizeFilter(long minLength, long skipCycles);
  // erases all size-filtering information that has been added up to now
  void resetSizeFilter();
  // replaces all size-filtering information that has been added up to now
  // with that in the new object (which is copied)
  void setSizeFilter(EcoFilterMinLength* sizeFilter);
  // returns a copy of the size filter
  EcoFilterMinLength* getSizeFilter();

  ParamsCycleManagement* copy();

  // adds the indicated trim to only the indicated cycle, zero-indexed; not erased by
  // basal trim, but will replace basal trim value
  void addTrim(int cycleNum, float coverageLevel);
  void addTrim(int cycleNum, float coverageLevel, long minLength);
  // numPriorCycles => number of cycles skipped before imposing trim.  remember
  // that trim is applied at the END of the cycle.  this works like the length
  // filter: a given trim will be applied until a cycle where another is specified
  void addBasalTrim(int numPriorCycles, float coverageLevel);
  void addBasalTrim(int numPriorCycles, float coverageLevel, long minLength);
  // these are applied only to input contigs prior to assembly
  void addInitialTrim(float coverageLevel);
  void addInitialTrim(float coverageLevel, long minLength);
  // remember: cycleNum is indexed from 1
  void executeTrim(set<ScoredSeq*>* contigSet, int cycleNum); // calls deepDelete on discarded
  void executeTrim(set<ScoredSeq*>* contigSet, int cycleNum, set<ScoredSeq*>* discardSet);
  void executeInitialTrim(set<ScoredSeq*>* contigSet); // calls deepDelete on discarded
  void executeInitialTrim(set<ScoredSeq*>* contigSet, set<ScoredSeq*>* discardSet);

  // regulates how often an output file is writting; must be >= 1
  int getOutputInterval();
  void setOutputInterval(int interval);
  // first cycle is 1
  bool outputThisCycle(int cycleNum);

  // if the first function returns True, then all of the contigs should be
  // re-set to "unchanged" and included as non-leftover input contigs for
  // that cycle.  the param carrier can be modified using the second function
  // to add points where that will happen.  cycles indexed from 1.
  bool shouldContigsReset(int cycleNum);
  void addContigReset(int cycleNum);

  void addRepeatDetection(int cycleNum, float numStdDevs, float minFoldIncrease, WhenInCycle whenInCycle,
			  long minRepeatSize, float matchFractId);
  void addRepeatDetection(int cycleNum, float numStdDevs, float minFoldIncrease, WhenInCycle whenInCycle,
			  long minRepeatSize, float matchFractId, string outfileName);
  bool seekDynamicRepeats(int cycleNum, WhenInCycle whenInCycle);
  // REQUIRES: seekDynamicRepeats is TRUE
  long repeatMinSize(int cycleNum, WhenInCycle whenInCycle);
  float repeatStandDevs(int cycleNum, WhenInCycle whenInCycle);
  float repeatMinFoldIncrease(int cycleNum, WhenInCycle whenInCycle);
  float repeatMatchFractId(int cycleNum, WhenInCycle whenInCycle);
  bool makeRepeatOutfile(int cycleNum, WhenInCycle whenInCycle);
  // REQUIRES: makeRepeatOutfile is TRUE
  string repeatOutfileName(int cycleNum, WhenInCycle whenInCycle);

 private:
  int _totalCycles;
  long _maxContigLinkages;
  EcoFilterMinLength* _sizeFilter;
  EcoFilter** _otherFilters;
  int _numOtherFilters;
  set<int> _resetCycles;

  // the sets have n+1 (one-indexed) cycle counts, while the arrays use
  // zero-indexed cycle counts
  set<int> _seekRepsBeforeCycles;
  long* _repMinSizeBefore;
  float* _repStdDevsBefore;
  float* _repMinFoldIncBefore;
  float* _repMatchFractIdBefore;
  bool* _repMakeFileBefore;
  string* _repFileNameBefore;
  set<int> _seekRepsAfterCycles;
  long* _repMinSizeAfter;
  float* _repStdDevsAfter;
  float* _repMinFoldIncAfter;
  float* _repMatchFractIdAfter;
  bool* _repMakeFileAfter;
  string* _repFileNameAfter;

  // constructor helpers
  void setupFilters();
  void setupTrimming();
  void setupTrimming(ParamsCycleManagement* pcm);
  void setupRepeatDetection();

  // for trimming
  void executeTrimHelper(set<ScoredSeq*>* contigSet, set<ScoredSeq*>* discardSet, float minScore, long minLength);
  ScoredSeq* trimThisContig(ScoredSeq* old, float minScore, long minLength);

  float* _cycleToTrimMinScore; // zero signifies no trimming
  long* _cycleToTrimMinLength; // 0 for no min length
  // helps with filling in the intermediate values: marks the beginning
  // of a block set by the basal trim method
  bool* _isTrimTransition;
  // keeps specifically-set values from being over-written
  bool* _wasTrimSetSpecifically;
  /*
  bool* _cycleWillTrim;
  float* _cycleToMinScore;
  long* _cycleToMinLength; // 0 for no min length
  */

  // for trimming of initial contigs
  //bool _initialWillTrim;
  float _initTrimMinScore;
  long _initTrimMinLength; // 0 for no min length

  // for outfile writing
  int _outputInterval;
};

#endif
