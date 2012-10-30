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


#ifndef MAPREADSETUPHELPER_H
#define MAPREADSETUPHELPER_H

#include <set>
#include <map>
#include <vector>
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
using namespace::std;

class ExtendJobMapper;

// used by mapReads to separate behavior between initial mapping and 
class MapReadSetupHelper{
 public:
  virtual ~MapReadSetupHelper();
  virtual int numFilesTotal() = 0;
  virtual int numFilesActive() = 0;
  virtual long totalReadCount() = 0;
  virtual bool isFileActive(int fileNum) = 0;
  virtual ExtendJobMapper* getNewEjm(set<ScoredSeqWithCarriers*>* parallelMapSet, int fileNum) = 0;
  virtual ScoredSeqWithCarriers*** getNewReadArrays(int* fileNumArray, int numFiles) = 0;
  virtual void getPairedReadArrays(ScoredSeqWithCarriers*** readArrays, ScoredSeqWithCarriers*** pairArrays,
				   long* numPairArray, int* fileNumArray, int numFiles) = 0;
  // returns the number of READS deemed un-mappable
  virtual long decideMappability(bool* shouldReadsMap, bool* shouldPairsMap, int fileNum,
				 ScoredSeqWithCarriers** reads, ScoredSeqWithCarriers** pairs, long rpCount) = 0;
  virtual void triggerFiles() = 0;
  virtual void updateActiveFiles() = 0;
  virtual bool includeContigCarrier(ScoredSeqWithCarriers* carrier) = 0;
};


class MrsMainHelper : public MapReadSetupHelper {
 public:
  MrsMainHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf,
		set<ScoredSeqWithCarriers*>* doNotExtendContigs);
  ~MrsMainHelper();
  int numFilesTotal();
  int numFilesActive();
  long totalReadCount();
  bool isFileActive(int fileNum);
  ExtendJobMapper* getNewEjm(set<ScoredSeqWithCarriers*>* parallelMapSet, int fileNum);
  ScoredSeqWithCarriers*** getNewReadArrays(int* fileNumArray, int numFiles);
  void getPairedReadArrays(ScoredSeqWithCarriers*** readArrays, ScoredSeqWithCarriers*** pairArrays,
			   long* numPairArray, int* fileNumArray, int numFiles);
  long decideMappability(bool* shouldReadsMap, bool* shouldPairsMap, int fileNum,
			 ScoredSeqWithCarriers** reads, ScoredSeqWithCarriers** pairs, long rpCount);
  void triggerFiles();
  void updateActiveFiles();
  bool includeContigCarrier(ScoredSeqWithCarriers* carrier);
 private:
  ReadPairFilter* _rpf;
  int _numReadFiles;
  int _numActiveFiles;
  long _totalFileSize;
  bool* _fileIsActive;
  ParamsMapping** _fileToPmap;
  ParameterizedReadFile** _readFileArray;
  set<ScoredSeqWithCarriers*> _doNotExtendContigs;
};

class MrsAddendumHelper : public MapReadSetupHelper {
 public:
  MrsAddendumHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf,
		    vector<ScoredSeqWithCarriers*>** hitReadsByFile);
  MrsAddendumHelper(ParameterizedReadFile** readFileArray, int numReadFiles, ReadPairFilter* rpf,
		    vector<ScoredSeqWithCarriers*>** hitReadsByFile,
		    vector<ScoredSeqWithCarriers*>** hitPairsByFile);
  ~MrsAddendumHelper();
  int numFilesTotal();
  int numFilesActive();
  long totalReadCount();
  bool isFileActive(int fileNum);
  ExtendJobMapper* getNewEjm(set<ScoredSeqWithCarriers*>* parallelMapSet, int fileNum);
  ScoredSeqWithCarriers*** getNewReadArrays(int* fileNumArray, int numFiles);
  void getPairedReadArrays(ScoredSeqWithCarriers*** readArrays, ScoredSeqWithCarriers*** pairArrays,
			   long* numPairArray, int* fileNumArray, int numFiles);
  long decideMappability(bool* shouldReadsMap, bool* shouldPairsMap, int fileNum,
			 ScoredSeqWithCarriers** reads, ScoredSeqWithCarriers** pairs, long rpCount);
  void triggerFiles();
  void updateActiveFiles();
  bool includeContigCarrier(ScoredSeqWithCarriers* carrier);
 private:
  ReadPairFilter* _rpf;
  int _numReadFiles;
  int _numActiveFiles;
  long _totalNumQueries;
  bool* _fileIsActive;
  ParamsMapping** _fileToPmap;
  ParameterizedReadFile** _readFileArray;

  // i am replacing the single data structure above with the parallel ones below
  ScoredSeqWithCarriers*** _readQueries;
  ScoredSeqWithCarriers*** _pairQueries;
  long* _numPairs;
  // keeps track of the current position in the arrays for extracting read pairs
  long* _readPairCounter;
};

#endif
