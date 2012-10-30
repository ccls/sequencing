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

#ifndef EXTENDCYCLE_H
#define EXTENDCYCLE_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "ReadFile.h"
#include "ReadPairFilter.h"
#include "DynamicProgrammingAligner.h"
#include "ScoredSeqCollectionBwt.h"
#include "AssemblyJobFactory.h"
#include "AssemblyJob.h"
#include "ScoredSeqNormalized.h"
#include "ScoredSeqWithCarriers.h"
#include "ScoredSeqFlip.h"
#include "ParamsMapping.h"
#include "ParamsMinOverlap.h"
#include "ParamsMinFractId.h"
#include "ParamsAlignment.h"
#include "ParamsDeBruijn.h"
#include "ParameterizedReadFile.h"
#include "AssemblerListener.h"
#include "ExtendJobMapper.h"
#include "EcoFilter.h"
#include "MapReadSetupHelper.h"
using namespace::std;
class ExtendJobMapper;
class MapReadSetupHelper;

class ExtendCycle {

 public:
  ExtendCycle(int cycleNum, set<ParameterizedReadFile*>* readFileSet, ReadPairFilter* rpfSingle, ReadPairFilter* rpfDouble,
	      set<ScoredSeq*>* inputContigs, AssemblerListener* listener);
  ExtendCycle(int cycleNum, set<ParameterizedReadFile*>* readFileSet, ReadPairFilter* rpfSingle, ReadPairFilter* rpfDouble,
	      set<ScoredSeq*>* inputContigs, set<ScoredSeq*>* doNotExtend, AssemblerListener* listener);
  ~ExtendCycle();

  // make these accessible but not mutable
  static int plusIndex();
  static int minusIndex();
  static int frontIndex();
  static int backIndex();

  // in case there is nothing to filter

  void runCycle(ParamsAlignment * paMini,
		ParamsAlignment * paMeta,
		ParamsDeBruijn* pdb,
		long maxContigLinkages,
		EcoFilter* preMetaFilter,
		EcoFilter* fullFilter);
  bool cycleHasRun();
  void getContigs(set<ScoredSeq*>* outputContigs);

 private:
  int _cycleNum;

  /*
  // hitReadsByFile can be set to NULL if no reads are to be collected
  void mapReads(MapReadSetupHelper* mrsh, vector<ScoredSeqWithCarriers*>** hitReadsByFile);
  */
  // the last three items can be set to NULL if no reads are to be collected
  // REQUIRES that each item in the read array is matched by its actual paired-end in the pair array
  void mapReads(MapReadSetupHelper* mrsh, vector<ScoredSeqWithCarriers*>** hitReadsByFile,
		vector<ScoredSeqWithCarriers*>** hitPairsByFile);
  void runMiniAssemblies(ParamsAlignment* pa, ParamsDeBruijn* pdb, long maxContigLinkages, EcoFilter* preMetaFilter);
  void runMetaAssembly(ParamsAlignment* pa, ParamsDeBruijn* pdb);

  // a simple helper method
  void clearLinkages(vector<ScoredSeqWithCarriers*>* linkedSeqs);

  AssemblerListener* _listener;
  ReadPairFilter* _rpfSingle;
  ReadPairFilter* _rpfDouble;

  bool _cycleHasRun;
  void constructorHelper(set<ParameterizedReadFile*>* readFileSet,
			 set<ScoredSeq*>* inputContigs,
			 set<ScoredSeq*>* doNotExtend);

  set<ReadFile*> _readFiles; // needs to be provided initially

  set<ScoredSeqNormalized*> _inputContigs;
  set<ScoredSeqNormalized*> _doNotExtendContigs;
  ParameterizedReadFile** _readFileArray;
  int _numReadFiles;

  int _numThreads;

  // for buffering and unbuffering
  map<ScoredSeq*,long> _seqToBufferUseCount;
  // these do the buffering/unbuffering
  void addBufferedSeqs(set<ScoredSeq*>* inputSet);
  void removeBufferedSeqs(set<ScoredSeq*>* inputSet);
  // and for overall usage (to delete)
  map<ScoredSeq*,long> _seqToTotalUseCount;
  void deleteReads(set<ScoredSeq*>* inputSet);

  // for general reference for the SSWC access to distinct collections
  static int _frontIndex; // 1
  static int _backIndex;  // 0
  static int _plusIndex;  // 0
  static int _minusIndex; // 1

  // for read mapping
  long _mappingStepSize;
  long _totalQueryNum;

  // WC encases Normalized encases Shallow
  //   carried: WC encasing Normalized encasing Read
  //              carried: WC encasing Normalized Shallow
  // 0->frontSet, 1->backSet
  set<ScoredSeqWithCarriers*> _mapSets;

  set<ScoredSeq*> _metaInputContigs;
  set<ScoredSeq*> _assembledContigs;

  void runOneMiniAssembly(AssemblyJob* aj, set<ScoredSeq*>* metaInputByThread, EcoFilter* preMetaFilter,
			  AssemblyJob::AssemblyThreadedness threadedness);

};

#endif
