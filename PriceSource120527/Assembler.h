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

/* This is the programmatic interface.
 */


#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <set>
#include <algorithm>
#include <vector>
#include "ScoredSeq.h"
#include "ReadFile.h"
#include "ParamsMapping.h"
#include "ParamsCycleManagement.h"
#include "ParamsAlignment.h"
#include "AssemblyJob.h"
#include "AssemblyJobFactory.h"
#include "ExtendCycle.h"
#include "OutputFile.h"
#include "ParameterizedReadFile.h"
#include "ParameterizedInitialFile.h"
#include "EcoFilter.h"
#include "AssemblerListener.h"
#include "ParamsDeBruijn.h"
#include "ReadPairFilter.h"
using namespace::std;

class Assembler {

 public:
  Assembler(int totalCycles);
  ~Assembler();

  void runAssembly();
  void getContigs(set<ScoredSeq*>* outputContigs);

  long getTotalSize();
  long getMaxSize();
  long getN50();

  static ParamsMinOverlap* getDefaultParamsMinOverlap();
  static ParamsMinOverlap* getNullParamsMinOverlap();

  void addInitialFile(string filename, int numSteps, int cyclesPerStep, float countFactor);
  void addInitialFile(long numSeqsToUse, string filename, int numSteps, int cyclesPerStep, float countFactor);
  void addInitialFileNoTarget(string filename, int numSteps, int cyclesPerStep, float countFactor);
  void addInitialFileNoTarget(long numSeqsToUse, string filename, int numSteps, int cyclesPerStep, float countFactor);

  // PAIRED ENDS
  // single files, alternating pairs
  void addReadFile(string filename, long ampSize);
  void addReadFile(string filename, long ampSize, float fractId);
  void addReadFileLimitUse(int startCycle, int numCycles, string filename, long ampSize);
  void addReadFileLimitUse(int startCycle, int numCycles, string filename, long ampSize, float fractId);
  void addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize);
  void addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize, float fractId);
  // two files
  void addReadFile(string filenameA, string filenameB, long ampSize);
  void addReadFile(string filenameA, string filenameB, long ampSize, float fractId);
  void addReadFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize);
  void addReadFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize, float fractId);
  void addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize);
  void addReadFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize, float fractId);

  // MATE PAIRS
  // single files, alternating pairs
  void addMatePairFile(string filename, long ampSize);
  void addMatePairFile(string filename, long ampSize, float fractId);
  void addMatePairFileLimitUse(int startCycle, int numCycles, string filename, long ampSize);
  void addMatePairFileLimitUse(int startCycle, int numCycles, string filename, long ampSize, float fractId);
  void addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize);
  void addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long ampSize, float fractId);
  // two files
  void addMatePairFile(string filenameA, string filenameB, long ampSize);
  void addMatePairFile(string filenameA, string filenameB, long ampSize, float fractId);
  void addMatePairFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize);
  void addMatePairFileLimitUse(int startCycle, int numCycles, string filenameA, string filenameB, long ampSize, float fractId);
  void addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize);
  void addMatePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filenameA, string filenameB, long ampSize, float fractId);

  // FALSE PAIRS
  void addFalsePairFile(string filename, long readLength, long ampSize);
  void addFalsePairFile(string filename, long readLength, long ampSize, float fractId);
  void addFalsePairFileLimitUse(int startCycle, int numCycles, string filename, long readLength, long ampSize);
  void addFalsePairFileLimitUse(int startCycle, int numCycles, string filename, long readLength, long ampSize, float fractId);
  void addFalsePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long readLength, long ampSize);
  void addFalsePairFileLimitUse(int startCycle, int onCycles, int skipCycles, string filename, long readLength, long ampSize, float fractId);



  void addOutputFile(string filename);

  // FILTER READS
  void addBadSequenceFilter(string badSeqFile, float fractId);
  void addHomopolymerFilter(long maxLength);
  void addDinucRepeatFilter(long maxLength);
  void addReadQualityFilter(float minFractGood, float minProbCorrect);
  void addReadQualityFilter(float minFractGood, float minProbCorrect, int numSkipCycles, int numRunCycles);
  void addReadCalledBasesFilter(float minFractGood);
  void addReadCalledBasesFilter(float minFractGood, int numSkipCycles, int numRunCycles);

  void findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			   long minRepeatSize, float matchFractId);
  void findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			   long minRepeatSize, float matchFractId, string outfileName);
  void findRepsAfterCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			  long minRepeatSize, float matchFractId);
  void findRepsAfterCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			  long minRepeatSize, float matchFractId, string outfileName);

  // FILTER INITIAL CONTIGS
  void icBadSequenceFilter(string badSeqFile, float fractId);
  void icBadSequenceFilter(string badSeqFile, float fractId, long maxSeqLength);
  void icHomopolymerFilter(long maxLength);
  void icHomopolymerFilter(long maxLength, long maxSeqLength);
  void icDinucRepeatFilter(long maxLength);
  void icDinucRepeatFilter(long maxLength, long maxSeqLength);
  void icReadQualityFilter(float minFractGood, float minProbCorrect);
  void icReadQualityFilter(float minFractGood, float minProbCorrect, long maxSeqLength);
  void icReadCalledBasesFilter(float minFractGood);
  void icReadCalledBasesFilter(float minFractGood, long maxSeqLength);


  // FILTER CONTIGS
  void addLengthFilter(long minLength, int cyclesToSkip);
  void addTargetMode(float minFractId, int cyclesToSkip, bool fullFile = false);
  void addTargetMode(float minFractId, int cyclesToSkip,
		     int numFilterCycles, int numSkipCycles, bool fullFile = false);

  void verboseLogFile(string filename);
  void makeLogNull();
  void makeLogVerbose();

  void setMaxThreadsPerFile(int maxThreadsPerFile);

  void setDeBruijnKmerSize(int kmerSize);
  void setDeBruijnMaxSeqLength(long maxSeqLength);
  void setDeBruijnMinSeqNum(long minSeqNum);

  long getLinkMax();
  void setLinkMax(long maxNumLinks);

  void getResetCycles(set<int>* resetCycles);
  void addResetCycle(int cycleNum);

  int getOutputInterval();
  void setOutputInterval(int interval);

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

  // adjust the min overlap requirements for mini assemblies...
  void setMiniMinOverlap(long minOvl);
  void setMiniMinOvlContigThreshold(long contigThreshold);
  void setMiniMinOvlScalingSlope(float slope);

  // ...and min fract ID for mini assemblies...
  void setMiniMinFractId(float minFractId);
  void setMiniFractIdContigThreshold(long contigThreshold);
  void setMiniFractIdScalingSlope(float slope);

  // ...and min fract ID for meta assemblies...
  void setMetaMinFractId(float minFractId);
  void setMetaFractIdContigThreshold(long contigThreshold);
  void setMetaFractIdScalingSlope(float slope);

  // ...and alignment score matrix values
  void setNucMatchScore(long score);
  void setNucMismatchPenalty(long penalty);
  void setOpenGapPenalty(long penalty);
  void setExtendGapPenalty(long penalty);
  //void setAlignmentScoreMatrix(AlignmentScoreMatrix* params);

 private:

  // consolidates all of the various inputs and checks for inconsistencies
  void runSetup();
  // cleans up at the end of a run; different from the destructor, it will leave params intact
  void runCleanup();

  // a helper method for finding repeats
  void detectRepeatsDynamically(int cycleNumP1, ParamsCycleManagement::WhenInCycle whenInCycle);

  void icBadSequenceFilterHelper(string badSeqFile, float fractId, bool useMaxLength, long maxSeqLength);

  // output files as they are being collected
  set<OutputFile*> _outfileSet;
  // the one that will be used; should be generated when assembly is run and deleted at the end of the run
  OutputFile* _outfile;

  // these two are used for set-up
  AssemblerListener* _stdOutListener;
  AssemblerListener* _fileListener;
  // this one is used by the job
  AssemblerListener* _listener;


  // If numInitialContigs==0, then all of the contigs in all the files will be used.
  // Same goes for if numInitialContigs > the num of sequences in the files.
  void generateInitialContigs(int cycleNum);
  void extendOneCycle(int cycleNum);
  void writeOutfile(string addendum);
  bool reverseN50Helper(long i,long j);
  // helps sort sequences - only modifies leftover by adding components of current
  void findLeftoverContigs(set<ScoredSeq*>* current, set<ScoredSeq*>* prior, set<ScoredSeq*>* leftover);

  // temp storage
  set<ParameterizedReadFile*> _inputPairedReadFiles;
  set<ParameterizedInitialFile*> _inputInitialContigFiles;
  // will not be targeted
  set<ParameterizedInitialFile*> _inputInitialContigFilesNoTarget;
  // possibly modified, used for actual run
  set<ParameterizedReadFile*> _pairedReadFiles;
  set<ParameterizedInitialFile*> _initialContigFiles;

  // FILTERING OF READS
  // collects the input
  vector<ReadPairFilter*> _waitingRpf;
  // NULL if no repeats have been searched for (or none found)
  ReadPairFilter* _dynamicRepeatFilter;

  // FILTERING OF INITIAL CONTIGS
  vector<ReadPairFilter*> _waitingIcf;


  // so that contigs can be added at the last minute
  EcoFilterInitialContigMatch* _waitingTargetMode;
  // is created immediately; EcoFilters are added to it directly
  ParamsCycleManagement* _cycleInfo;

  set<ScoredSeq*> _priorCycleContigs;
  set<ScoredSeq*> _currentCycleContigs;
  set<ScoredSeq*> _priorCycleLeftovers;
  set<ScoredSeq*> _currentCycleLeftovers;

  // the alignment score parameters are the same for both assembly steps
  AlignmentScoreMatrix* _alScoreMatrix;

  // these three fields will be stored...
  ParamsMinOverlap* _miniMinOverlapParams;
  ParamsMinFractId* _miniMinFractIdParams;
  // ...and used to create the transient field below
  ParamsAlignment* _miniAlignParams;

  // these two fields will be stored...
  ParamsMinFractId* _metaMinFractIdParams;
  // ...and used to create the transient field below (OVL is null)
  ParamsAlignment* _metaAlignParams;

  ParamsDeBruijn* _paramsDeBruijn;

  int _maxThreadsPerFile;
};

#endif
