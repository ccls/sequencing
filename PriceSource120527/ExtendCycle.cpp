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

#ifndef EXTENDCYCLE_CPP
#define EXTENDCYCLE_CPP

#include "ExtendCycle.h"
#include "ExtendJobCreator.h"
#include "Read.h"
#include "ScoredSeqSubseq.h"
#include "ScoredSeqNested.h"
#include "ScoredSeqMonoScore.h"
#include <typeinfo>
#include <sstream>
#include <omp.h>
#include <limits.h>
using namespace::std;

int ExtendCycle::_plusIndex = 0;
int ExtendCycle::_minusIndex = 1;
int ExtendCycle::_frontIndex = 1;
int ExtendCycle::_backIndex = 0;

int ExtendCycle::plusIndex(){ return _plusIndex; }
int ExtendCycle::minusIndex(){ return _minusIndex; }
int ExtendCycle::frontIndex(){ return _frontIndex; }
int ExtendCycle::backIndex(){ return _backIndex; }


void ExtendCycle::constructorHelper(set<ParameterizedReadFile*>* readFileSet,
				    set<ScoredSeq*>* inputContigs,
				    set<ScoredSeq*>* doNotExtend){

  _cycleHasRun = false;

  _numReadFiles = readFileSet->size();
  _readFileArray = new ParameterizedReadFile*[ _numReadFiles ];
  int rfCounter = 0;
  for (set<ParameterizedReadFile*>::iterator rfIt = readFileSet->begin(); rfIt != readFileSet->end(); ++rfIt){
    _readFileArray[rfCounter] = (*rfIt);
    rfCounter++;
  }

  // make shallow copies of the input contigs so that they can be deleted
  for (set<ScoredSeq*>::iterator icIt = inputContigs->begin(); icIt != inputContigs->end(); ++icIt){
    //ScoredSeqNormalized* inSeq = new ScoredSeqNormalized( (*icIt)->shallowCopy() );
    ScoredSeqNormalized* inSeq = new ScoredSeqNormalized( *icIt );
    _inputContigs.insert( inSeq );
    if ( doNotExtend->count( (*icIt) ) != 0 ){ _doNotExtendContigs.insert( inSeq ); }
  }

  _numThreads = omp_get_max_threads();
}

ExtendCycle::ExtendCycle(int cycleNum, set<ParameterizedReadFile*>* readFileSet, ReadPairFilter* rpfSingle, ReadPairFilter* rpfDouble, set<ScoredSeq*>* inputContigs, AssemblerListener* listener) :
  _cycleNum(cycleNum),
  _listener(listener) {
  _rpfSingle = rpfSingle->copy();
  _rpfDouble = rpfDouble->copy();
  set<ScoredSeq*> emptyDoNotExtend;
  constructorHelper(readFileSet, inputContigs, &emptyDoNotExtend);
}

ExtendCycle::ExtendCycle(int cycleNum, set<ParameterizedReadFile*>* readFileSet, ReadPairFilter* rpfSingle, ReadPairFilter* rpfDouble, set<ScoredSeq*>* inputContigs, set<ScoredSeq*>* doNotExtend, AssemblerListener* listener) :
  _cycleNum(cycleNum),
  _listener(listener) {
  _rpfSingle = rpfSingle->copy();
  _rpfDouble = rpfDouble->copy();
  constructorHelper(readFileSet, inputContigs, doNotExtend);
}

ExtendCycle::~ExtendCycle(){
  delete _rpfSingle;
  delete _rpfDouble;
  delete [] _readFileArray;
  // get rid of the wrappers of the contigs (the map set outermost wrappers)
  if (_mapSets.size() != 0){
    throw AssemblyException::ImplementationError("EC:destructor, mapSets still exist; these should have already been deleted."); 
  }
}

void ExtendCycle::runCycle(ParamsAlignment * paMini,
			   ParamsAlignment * paMeta,
			   ParamsDeBruijn* pdb,
			   long maxContigLinkages,
			   EcoFilter* preMetaFilter,
			   EcoFilter* fullFilter){

  // set up the carriers that will hold the read-contig mapping linkages
  set<ScoredSeqWithCarriers*> dneWithCarriers;
  for (set<ScoredSeqNormalized*>::iterator it = _inputContigs.begin(); it != _inputContigs.end(); ++it){
    ScoredSeqNormalized* contig = (*it); // the contig will be normalized
    ScoredSeqWithCarriers* carrier = new ScoredSeqWithCarriers(contig, 2); // 0->frontSet, 1->backSet
    _mapSets.insert( carrier );
    if ( _doNotExtendContigs.count(contig) != 0){ dneWithCarriers.insert(carrier); }
  }

  // first, map to the contig edges of "alive" contigs to establish the
  // set of relevant reads
  MapReadSetupHelper* mrsh;
  ParameterizedReadFile** nullRfa = new ParameterizedReadFile*[1];
  nullRfa[0] = NULL;
  if (dneWithCarriers.size() == _inputContigs.size()){
    mrsh = new MrsMainHelper(nullRfa, 0, _rpfDouble, &dneWithCarriers);
  } else {
    mrsh = new MrsMainHelper(_readFileArray, _numReadFiles, _rpfDouble, &dneWithCarriers);
  }

  // this is used to collect reads that mapped from the first job and pass them to the second
  vector<ScoredSeqWithCarriers*>** hitReadsByFile = new vector<ScoredSeqWithCarriers*>*[ _numReadFiles ];
  vector<ScoredSeqWithCarriers*>** hitPairsByFile = new vector<ScoredSeqWithCarriers*>*[ _numReadFiles ];
  for (int n = 0; n < _numReadFiles; ++n){ hitReadsByFile[n] = new vector<ScoredSeqWithCarriers*>; }
  for (int n = 0; n < _numReadFiles; ++n){ hitPairsByFile[n] = new vector<ScoredSeqWithCarriers*>; }

  // the FIRST mapping
  mapReads(mrsh, hitReadsByFile, hitPairsByFile);
  delete mrsh;

  // second, map the relevant reads to ALL of the contigs.  "relevant" includes
  // the reads that matched in the first step plus their paired ends.  first,
  // the matches from the first step must be erased; they will either be re-discovered
  // or replaced by better-scoring matches to other parts of the assembly
  for (int n = 0; n < _numReadFiles; ++n){ clearLinkages(hitReadsByFile[n]); }
  for (int n = 0; n < _numReadFiles; ++n){ clearLinkages(hitPairsByFile[n]); }
  mrsh = new MrsAddendumHelper(_readFileArray, _numReadFiles, _rpfSingle, hitReadsByFile, hitPairsByFile);

  // i can now get rid of these
  for (int n = 0; n < _numReadFiles; ++n){
    delete hitReadsByFile[n];
    delete hitPairsByFile[n];
  }
  delete [] hitReadsByFile;
  delete [] hitPairsByFile;

  // the SECOND mapping - no reads need to be collected this time, so NULL is used
  mapReads(mrsh, NULL, NULL);
  delete mrsh;

  // now, perform the assemblies
  runMiniAssemblies(paMini, pdb, maxContigLinkages, preMetaFilter);
  runMetaAssembly(paMeta, pdb);
  // use the full filter to filter the potential output
  _metaInputContigs.clear();
  set<ScoredSeq*> failingContigs;
  fullFilter->filterSeqs(&_assembledContigs, &_metaInputContigs, &failingContigs, _cycleNum, EcoFilter::threaded);
  // i don't need to do anything if there were no sequences filtered out
  if (failingContigs.size() != 0){
    _assembledContigs.clear();
    for (set<ScoredSeq*>::iterator it = failingContigs.begin(); it != failingContigs.end(); ++it){ (*it)->deepDelete(); }
    runMetaAssembly(paMeta, pdb);
  }

  delete [] nullRfa;

  _cycleHasRun = true;
}



void ExtendCycle::clearLinkages(vector<ScoredSeqWithCarriers*>* linkedSeqs){
  for (vector<ScoredSeqWithCarriers*>::iterator it = linkedSeqs->begin(); it != linkedSeqs->end(); ++it){
    (*it)->clearCarriedSeqs(0);
    (*it)->clearCarriedSeqs(1);
  }
}








void ExtendCycle::mapReads(MapReadSetupHelper* mrsh, vector<ScoredSeqWithCarriers*>** hitReadsByFile,
			   vector<ScoredSeqWithCarriers*>** hitPairsByFile){

  bool returnMatchPairs;
  if (hitReadsByFile==NULL or hitPairsByFile==NULL){
    if (hitReadsByFile!=NULL or hitPairsByFile!=NULL){
      throw AssemblyException::ArgError("EC::mapReads, the return-seq arrays should be either both NULL or neither NULL");
    }
    returnMatchPairs = false;
  } else { returnMatchPairs = true; }

  // no need to do anything if there are no active files
  if (mrsh->numFilesActive() == 0){
    _listener->beginMapping(0);
    _listener->endMapping(0);
  } else {
    // for building arrays
    int numThreadsP1 = _numThreads + 1;

    // set up for mapping - one EJM for each file, for each thread
    _listener->verboseMessage("Creating mapper objects...");
    ExtendJobMapper** fileToEjm = new ExtendJobMapper*[ _numReadFiles+1];
    #pragma omp parallel for schedule(dynamic)
    for (int fileN = 0; fileN < mrsh->numFilesTotal(); ++fileN){
      if ( mrsh->isFileActive(fileN) ){
	fileToEjm[ fileN ] = mrsh->getNewEjm(&_mapSets, fileN);
      } else { fileToEjm[ fileN ] = 0; }
    }

    _listener->beginMapping(mrsh->totalReadCount());
    // keeps track of the number of reads that are filtered out by RPFs
    long readsFiltered = 0;

    while (mrsh->numFilesActive() > 0){
      // prepare a subset of reads for mapping - fill these in below
      ScoredSeqWithCarriers*** readsToMap = new ScoredSeqWithCarriers**[ numThreadsP1 ];
      ScoredSeqWithCarriers*** pairsToMap = new ScoredSeqWithCarriers**[ numThreadsP1 ];
      long numPairs[ numThreadsP1 ];
      long seqsWithHitsCount[ numThreadsP1 ];

      // get the files to use this round, the max num of files allowed by the num threads
      set<int> thisRoundFileNums;
      int fn = 0;
      while (fn < mrsh->numFilesTotal() and thisRoundFileNums.size() < _numThreads){
	if ( mrsh->isFileActive(fn) ){ thisRoundFileNums.insert(fn); }
	++fn;
      }

      // allocate threads to each file
      int numRoundFiles = thisRoundFileNums.size();
      int numThreadsUsed = 0;
      int threadToFileNum[numThreadsP1];

      // set up an array of targets
      ExtendJobMapper* ejmArray[ numThreadsP1 ];

      _listener->verboseMessage("getting and buffering reads...");
      while (numThreadsUsed < _numThreads and thisRoundFileNums.size() > 0){
        // get the files ready for read extraction by creating local thread arrays
	int numLocalThreads = thisRoundFileNums.size();
	if (numLocalThreads + numThreadsUsed > _numThreads){ numLocalThreads = _numThreads - numThreadsUsed; }
	int localIndex = 0;
	for (set<int>::iterator it = thisRoundFileNums.begin(); it != thisRoundFileNums.end(); ++it){
	  if (localIndex < numLocalThreads){
	    threadToFileNum[localIndex + numThreadsUsed] = *it;
	    localIndex++;
	  }
	}

        // get the queries for each thread (threaded by file) - do not thread this, the mrsh will thread it
	long numLocalThreadsP1 = numLocalThreads + 1;
	int* localFileN = new int[ numLocalThreadsP1 ];
	int* localThreadN = new int[ numLocalThreadsP1 ];
	for (int localN = 0; localN < numLocalThreads; ++localN){
  	  // get all of the necessary local variables
	  int threadN = localN + numThreadsUsed;
	  int fileN = threadToFileNum[threadN];
	  ejmArray[threadN] = fileToEjm[fileN];
	  localFileN[localN] = fileN;
	  localThreadN[localN] = threadN;
	}

	// create the local versions of these data structures and fill using the mrsh
	ScoredSeqWithCarriers*** localReadsToMap = new ScoredSeqWithCarriers**[ numLocalThreadsP1 ];
	ScoredSeqWithCarriers*** localPairsToMap = new ScoredSeqWithCarriers**[ numLocalThreadsP1 ];
	long* localNumPairs = new long[ numLocalThreadsP1 ];
	mrsh->getPairedReadArrays(localReadsToMap, localPairsToMap, localNumPairs, localFileN, numLocalThreads);
	for (int n = 0; n < numLocalThreads; ++n){
	  readsToMap[localThreadN[n]] = localReadsToMap[n];
	  pairsToMap[localThreadN[n]] = localPairsToMap[n];
	  numPairs[ localThreadN[n] ] = localNumPairs[n];
	}
	delete [] localReadsToMap;
	delete [] localPairsToMap;
	delete [] localNumPairs;
	delete [] localFileN;
	delete [] localThreadN;
	numThreadsUsed += numLocalThreads;

        // get rid of the finished file nums
	mrsh->updateActiveFiles();
	set<int> fileNumsToRemove;
	for (set<int>::iterator it = thisRoundFileNums.begin(); it != thisRoundFileNums.end(); ++it){
	  if (! mrsh->isFileActive(*it)){ fileNumsToRemove.insert( *it ); }
	}
	for (set<int>::iterator it = fileNumsToRemove.begin(); it != fileNumsToRemove.end(); ++it){ thisRoundFileNums.erase( *it ); }
      }

      // trigger the buffering of sequences OUTSIDE the parallel block
      _listener->verboseMessage("triggering files...");
      mrsh->triggerFiles();

      // do the mapping
      long readsBeingMapped = 0;
      for (int n = 0; n < numThreadsUsed; ++n){ readsBeingMapped += numPairs[n] * 2; }

      _listener->verboseMessage("mapping reads...");
      // if either read or pair has a match, both will be output
      bool** eitherHasMatch = new bool*[ numThreadsUsed+1 ];

      // do the mapping, and decide what needs to be deleted, in a threaded block
      long* numReadsFiltered = new long[numThreadsUsed+1];

      #pragma omp parallel for schedule(dynamic)
      for (int threadN = 0; threadN < numThreadsUsed; ++threadN){
	bool* shouldReadsMap = new bool[ numPairs[threadN]+1 ];
	bool* shouldPairsMap = new bool[ numPairs[threadN]+1 ];
	numReadsFiltered[threadN] = mrsh->decideMappability(shouldReadsMap, shouldPairsMap, threadToFileNum[threadN],
							    readsToMap[threadN], pairsToMap[threadN], numPairs[threadN]);
	bool* readHasMatch = new bool[ numPairs[threadN]+1 ];
	bool* pairHasMatch = new bool[ numPairs[threadN]+1 ];
	eitherHasMatch[threadN] = new bool[ numPairs[threadN]+1 ];
	seqsWithHitsCount[threadN] = ejmArray[threadN]->mapPairedQueries(readsToMap[threadN], pairsToMap[threadN],
									 shouldReadsMap, shouldPairsMap,
									 readHasMatch, pairHasMatch, numPairs[threadN]);
	// decide if the unmapped read pairs should be deleted
	for (long n = 0; n < numPairs[ threadN ]; ++n){
	  eitherHasMatch[threadN][n] = (readHasMatch[n] or pairHasMatch[n]);
	}
	delete [] readHasMatch;
	delete [] pairHasMatch;
	delete [] shouldReadsMap;
	delete [] shouldPairsMap;
      }

      // account for the filtered reads
      for (int n = 0; n < numThreadsUsed; ++n){ readsFiltered += numReadsFiltered[n]; }
      delete [] numReadsFiltered;

      // the following operations can be threaded by file - so I need an array of file nums
      set<int> usedFileSet;
      for (int n = 0; n < numThreadsUsed; ++n){ usedFileSet.insert( threadToFileNum[n] ); }
      int numUsedFiles = usedFileSet.size();
      int* usedFiles = new int[numUsedFiles+1];
      int ufCount = 0;
      for (set<int>::iterator it = usedFileSet.begin(); it != usedFileSet.end(); ++it){
	usedFiles[ufCount] = *it;
	++ufCount;
      }

      #pragma omp parallel for schedule(dynamic)
      for (int ufN = 0; ufN < numUsedFiles; ++ufN){
	int fileN = usedFiles[ufN];
	// iterate through all of the threads using this file
	for (int threadN = 0; threadN < numThreadsUsed; ++threadN){
	  if (fileN == threadToFileNum[threadN]){
	    // either just deletion, or deletion + return
	    if (returnMatchPairs){
	      for (long n = 0; n < numPairs[ threadN ]; ++n){
		if (eitherHasMatch[threadN][n]){
		  hitReadsByFile[fileN]->push_back( readsToMap[threadN][n] );
		  hitPairsByFile[fileN]->push_back( pairsToMap[threadN][n] );
		} else {
		  readsToMap[threadN][n]->deepDelete();
		  pairsToMap[threadN][n]->deepDelete();
		}
	      }
	    } else {
	      for (long n = 0; n < numPairs[ threadN ]; ++n){
		if (! eitherHasMatch[threadN][n]){
		  readsToMap[threadN][n]->deepDelete();
		  pairsToMap[threadN][n]->deepDelete();
		}
	      }
	    }
	  }
	}
      }
      delete [] usedFiles;

      // clean up the EJMapper objects for non-active files.  I will set them
      // to zero so that I can still call delete on them in later mapping rounds.
      _listener->verboseMessage("cleaning up...");
      mrsh->updateActiveFiles();
      for (int n = 0; n < mrsh->numFilesTotal(); ++n){
	if (! mrsh->isFileActive(n)){
	  delete fileToEjm[n];
	  fileToEjm[n] = NULL;
	}
      }
      for (int n = 0; n < numThreadsUsed; ++n){
	delete [] readsToMap[n];
	delete [] pairsToMap[n];
	delete [] eitherHasMatch[n];
      }
      delete [] readsToMap;
      delete [] pairsToMap;
      delete [] eitherHasMatch;

      _listener->updateMapping(readsBeingMapped);
    }
    delete [] fileToEjm;

    _listener->endMapping(readsFiltered);
  }
}



void ExtendCycle::addBufferedSeqs(set<ScoredSeq*>* inputSet){
  #pragma omp critical (ECManageBuffer)
  {
    for (set<ScoredSeq*>::iterator seqIt = inputSet->begin(); seqIt != inputSet->end(); ++seqIt){
      map<ScoredSeq*,long>::iterator foundIt = _seqToBufferUseCount.find( *seqIt );
      if (foundIt == _seqToBufferUseCount.end()){
	_seqToBufferUseCount.insert( pair<ScoredSeq*,long>(*seqIt,1) );
	(*seqIt)->bottomBuffer();
      } else { foundIt->second++; }
    }
  }
}
void ExtendCycle::removeBufferedSeqs(set<ScoredSeq*>* inputSet){
  #pragma omp critical (ECManageBuffer)
  {
    for (set<ScoredSeq*>::iterator seqIt = inputSet->begin(); seqIt != inputSet->end(); ++seqIt){
      map<ScoredSeq*,long>::iterator foundIt = _seqToBufferUseCount.find( *seqIt );
      if (foundIt == _seqToBufferUseCount.end()){
	throw AssemblyException::LogicError("EC::removeBufferedSeqs cannot remove a non-existant buffered seq");
      }
      if (foundIt->second == 0){
	throw AssemblyException::LogicError("EC::removeBufferedSeqs cannot remove a zero-count buffered seq");
      }
      foundIt->second--;
      if (foundIt->second == 0){
	_seqToBufferUseCount.erase( foundIt );
	(*seqIt)->deepUnbuffer();
      }
    }
  }
}



void ExtendCycle::runOneMiniAssembly(AssemblyJob* aj, set<ScoredSeq*>* metaInputByThread, EcoFilter* preMetaFilter,
				     AssemblyJob::AssemblyThreadedness threadedness){
  // the input params of this factory need to be further utilized
  int threadNum = omp_get_thread_num();

  ostringstream outMessage1;
  outMessage1 << "Thread " << threadNum << ": Job Start";
  _listener->verboseMessage(outMessage1.str());

  // generate the buffer-threaded enum data type
  ReadFile::BufferThreadedness isBufferThreaded;
  if (threadedness == AssemblyJob::THREADED){ isBufferThreaded = ReadFile::bufferThreaded; }
  else { isBufferThreaded = ReadFile::bufferNotThreaded; }

  // make shallow copies of everything so that they can be deleted without worry
  set<ScoredSeq*> rawInput;
  aj->allInputSeqs( &rawInput );
  addBufferedSeqs( &rawInput );
  for (int n = 0; n < _numReadFiles; ++n){ ReadFileCommunicator::bufferFill(_readFileArray[n], isBufferThreaded); }
  aj->makeShallow(false);
  _listener->assemblyJobBeginStats(aj);
  removeBufferedSeqs( &rawInput );
  deleteReads( &rawInput );

  aj->runJob(threadedness);
  set<ScoredSeq*> unfilteredOutput;
  aj->allAssembledSeqs( &unfilteredOutput );

  // get rid of the too-$hort sequences
  set<ScoredSeq*> filteredMetaInput;
  set<ScoredSeq*> unsuitableMetaInput;

  // get the input data that was combined into the output contigs and delete it
  set<ScoredSeq*> discardedInput;
  aj->discardedInputSeqs( &discardedInput );
  for (set<ScoredSeq*>::iterator disIt = discardedInput.begin(); disIt != discardedInput.end(); ++disIt){ (*disIt)->deepDelete(); }

  _listener->assemblyJobEndStats(aj);
  ostringstream outMessage2;
  outMessage2 << "Thread " << threadNum << ": Job End";
  _listener->verboseMessage(outMessage2.str());
  _listener->updateEdgeAssembly(1);

  // deal with output
  if (isBufferThreaded){
    preMetaFilter->filterSeqs(&unfilteredOutput, metaInputByThread, &unsuitableMetaInput, _cycleNum, EcoFilter::threaded);
  } else {
    preMetaFilter->filterSeqs(&unfilteredOutput, metaInputByThread, &unsuitableMetaInput, _cycleNum, EcoFilter::notThreaded);
  }
  for (set<ScoredSeq*>::iterator cit = unsuitableMetaInput.begin(); cit != unsuitableMetaInput.end(); ++cit){ (*cit)->deepDelete(); }
}

void ExtendCycle::deleteReads(set<ScoredSeq*>* inputSet){
  #pragma omp critical (ECDeleteReads)
  {
    for (set<ScoredSeq*>::iterator seqIt = inputSet->begin(); seqIt != inputSet->end(); ++seqIt){
      // now delete, if used
      long totalCount = dynamic_cast<ScoredSeqNormalized*>( *seqIt )->getCount();
      map<ScoredSeq*,long>::iterator foundIt = _seqToTotalUseCount.find( *seqIt );
      if (foundIt != _seqToTotalUseCount.end()){
	foundIt->second++;
	if (foundIt->second == totalCount){
	  (*seqIt)->deepDelete();
	  _seqToTotalUseCount.erase( foundIt );
	}
      }
    }
  }
}


void ExtendCycle::runMiniAssemblies(ParamsAlignment* pa, ParamsDeBruijn* pdb, long maxContigLinkages, EcoFilter* preMetaFilter){

  AssemblyJobFactory * overhangFactory = new AssemblyJobFactory(pa, pdb, AssemblyJob::SINGLESTRANDED, _listener); //single-stranded

  // CREATE THE JOBS (this function contains threaded sections)
  set<ScoredSeqWithCarriers*> garbageDump;
  vector<AssemblyJob*>* jobSet = ExtendJobCreator::createJobs(&_mapSets, &_doNotExtendContigs,
							      overhangFactory, maxContigLinkages,
							      &garbageDump);

  set<AssemblyJob*> tempSet;
  for (vector<AssemblyJob*>::iterator it = jobSet->begin(); it != jobSet->end(); ++it){ tempSet.insert(*it); }
  jobSet->clear();
  for (set<AssemblyJob*>::iterator it = tempSet.begin(); it != tempSet.end(); ++it){ jobSet->push_back(*it); }

  delete overhangFactory;
  // get rid of the wrappers of the contigs (the map set outermost wrappers)
  for (set<ScoredSeqWithCarriers*>::iterator inIt = _mapSets.begin(); inIt != _mapSets.end(); ++inIt){
    delete (*inIt); // shallow delete since inner contigs were already deleted
  }
  _mapSets.clear();

  // assembly jobs have been constructed, so carrier objects can be deleted and
  // nested objects stopred for future deletion
  set<ScoredSeq*> nestedGarbageDump;
  for (set<ScoredSeqWithCarriers*>::iterator it = garbageDump.begin(); it != garbageDump.end(); ++it){
    nestedGarbageDump.insert( (*it)->getNested() );
    delete *it;
  }

  // make an initial unsorted version of the array
  AssemblyJob** jobArray = new AssemblyJob*[ jobSet->size() + 1 ];
  long totalNumJobs = 0;
  for (vector<AssemblyJob*>::iterator jobIt = jobSet->begin(); jobIt != jobSet->end(); ++jobIt){
    jobArray[totalNumJobs] = *jobIt;
    totalNumJobs++;
  }
  delete jobSet;

  // i want the job array to be sorted by size, so first i need to know
  // how many jobs of each size there are.  while i'm at it, i am also
  // adding normalization counts to the SSNormalized job components.
  // NOTE: because of buffering, this loop is not thread-able, despite how it may first appear
  map<long, long> sizeToJobCount;
  long* arrayOfJobSizes = new long[ totalNumJobs ];
  long totalOfJobSizes = 0;
  for (long oc = 0; oc < totalNumJobs; ++oc){
    set<ScoredSeq*> inputSeqs;
    jobArray[oc]->allInputSeqs( &inputSeqs );

    for (set<ScoredSeq*>::iterator seqIt = inputSeqs.begin(); seqIt != inputSeqs.end(); ++seqIt){
      dynamic_cast<ScoredSeqNormalized*>( *seqIt )->addCount();
      // this sequence will be deleted separately
      nestedGarbageDump.erase( *seqIt );
      // only delete it locally if it is not an input sequence
      if (_inputContigs.count(dynamic_cast<ScoredSeqNormalized*>(*seqIt))==0){
	_seqToTotalUseCount.insert( pair<ScoredSeq*,long>( *seqIt, 0 ) );
      }
    }

    // DEFINE the job size as a long
    // i will sort by the product of the biggest contig times the summed lengths of all the rest
    // jobs with just one contig are null, and they will be appropriately sized as zero by this metric
    // note: getting sizes will require buffering

    long jobSize = 0;
    if (inputSeqs.size() > 1){
      addBufferedSeqs( &inputSeqs );
      for (int n = 0; n < _numReadFiles; ++n){ ReadFileCommunicator::bufferFill(_readFileArray[n], ReadFile::bufferThreaded); }
      long totalContigSize = 0;
      long maxContigSize = 0;
      for (set<ScoredSeq*>::iterator seqIt = inputSeqs.begin(); seqIt != inputSeqs.end(); ++seqIt){
	long contigSize = (*seqIt)->size();
	totalContigSize += contigSize;
	if (contigSize > maxContigSize){ maxContigSize = contigSize; }
      }
      removeBufferedSeqs( &inputSeqs );
      if (maxContigSize != 0){
	long upperBound = LONG_MAX / maxContigSize;
	long remainingLength = totalContigSize - maxContigSize;
	if (remainingLength < upperBound){
	  jobSize = maxContigSize * (totalContigSize - maxContigSize);
	} else { jobSize = LONG_MAX; }
      }
    }

    arrayOfJobSizes[oc] = jobSize;
    totalOfJobSizes += jobSize;
    map<long,long>::iterator stjcIt = sizeToJobCount.find( jobSize );
    if (stjcIt == sizeToJobCount.end()){ sizeToJobCount.insert( pair<long,long>(jobSize,1) ); }
    else { stjcIt->second++; }
  }

  // delete the nested reads that are not part of jobs
  for (set<ScoredSeq*>::iterator it = nestedGarbageDump.begin(); it != nestedGarbageDump.end(); ++it){
    (*it)->deepDelete();
  }


  // SORT THE JOBS

  // create another data structure to make sorting of the list easier
  map<long, long> sizeToStartIndex;
  map<long, long>::iterator stsiIt = sizeToStartIndex.begin();
  long currentStartIndex = totalNumJobs;
  for (map<long,long>::iterator it = sizeToJobCount.begin(); it != sizeToJobCount.end(); ++it){
    currentStartIndex -= it->second;
    stsiIt = sizeToStartIndex.insert( stsiIt, pair<long,long>(it->first, currentStartIndex) );
  }

  // now i can use the information that i gathered above to sort the array
  AssemblyJob** sortedJobArray = new AssemblyJob*[ totalNumJobs ];
  long* arrayOfSortedJobSizes = new long[ totalNumJobs ];
  for (long oc = 0; oc < totalNumJobs; ++oc){
    map<long,long>::iterator indexIt = sizeToStartIndex.find( arrayOfJobSizes[oc] );
    arrayOfSortedJobSizes[ indexIt->second ] = arrayOfJobSizes[oc];
    sortedJobArray[ indexIt->second ] = jobArray[oc];
    indexIt->second++;
  }

  // replace the unsorted with the sorted versions
  delete [] jobArray;
  jobArray = sortedJobArray;
  delete [] arrayOfJobSizes;
  arrayOfJobSizes = arrayOfSortedJobSizes;

  // organize the jobs into those to be threaded and those to not be threaded
  long remainingJobSizeTotal = totalOfJobSizes;
  long startParallelIndex = 0;
  int numThreads = omp_get_max_threads();
  // if there is only one thread, no jobs will be run in parallel
  if (numThreads == 1){ startParallelIndex = totalNumJobs; }
  else if (totalNumJobs > 0){
    remainingJobSizeTotal -= arrayOfJobSizes[startParallelIndex];
    while (startParallelIndex < totalNumJobs and arrayOfJobSizes[startParallelIndex] * numThreads >= remainingJobSizeTotal){
      startParallelIndex++;
      if (startParallelIndex < totalNumJobs){ remainingJobSizeTotal -= arrayOfJobSizes[startParallelIndex]; }
    }
  }
  if (totalNumJobs - startParallelIndex < numThreads){ startParallelIndex = totalNumJobs; }

  _listener->beginEdgeAssembly(totalNumJobs);

  // an output collection for each thread; this will be consolidated after
  set<ScoredSeq*>* metaInputByThread[ _numThreads ];
  for (int tn = 0; tn < _numThreads; ++tn){ metaInputByThread[tn] = new set<ScoredSeq*>; }


  // RUN THE JOBS

  // first run the long jobs individually threaded
  for (long oc = 0; oc < startParallelIndex; ++oc){
    runOneMiniAssembly(jobArray[oc], metaInputByThread[0], preMetaFilter, AssemblyJob::THREADED);
  }

  // now run the smaller jobs in parallel (biggest launched first)
  #pragma omp parallel for schedule(dynamic)
  for (long oc = startParallelIndex; oc < totalNumJobs; ++oc){
    int threadNum = omp_get_thread_num();
    runOneMiniAssembly(jobArray[oc], metaInputByThread[threadNum], preMetaFilter, AssemblyJob::NOTTHREADED);
  }

  // clean up array out here
  for (long oc = 0; oc < totalNumJobs; ++oc){ delete jobArray[oc]; }
  delete [] jobArray;
  delete [] arrayOfJobSizes;

  // consolidate the output and clean up
  for (int n = 0; n < _numThreads; ++n){
    _metaInputContigs.insert( metaInputByThread[n]->begin(), metaInputByThread[n]->end() );
    delete metaInputByThread[n];
  }

  _listener->endEdgeAssembly();
}





void ExtendCycle::runMetaAssembly(ParamsAlignment* pa, ParamsDeBruijn* pdb){
  _listener->beginMetaAssembly( _metaInputContigs.size() );
  AssemblyJobFactory * contigFactory = new AssemblyJobFactory(pa, pdb, AssemblyJob::DOUBLESTRANDED, _listener);
  AssemblyJob* aj = contigFactory->redundancyJob( &_metaInputContigs );

  //singleAssemblyHelper(aj, &_assembledContigs, AssemblyJob::threaded);

  aj->runJob(AssemblyJob::THREADED);
  aj->allAssembledSeqs( &_assembledContigs );

  // get the input data that was combined into the output contigs and delete it
  set<ScoredSeq*> discardedInput;
  aj->discardedInputSeqs( &discardedInput );
  for (set<ScoredSeq*>::iterator disIt = discardedInput.begin(); disIt != discardedInput.end(); ++disIt){
    (*disIt)->deepDelete();
  }

  delete aj;
  delete contigFactory;

  // delete all of the input contigs - this will include the _doNotExtendContigs
  // shallow delete since I am only deleting the shells
  for (set<ScoredSeqNormalized*>::iterator icIt = _inputContigs.begin(); icIt != _inputContigs.end(); ++icIt){ delete *icIt; }
  _inputContigs.clear();

  _listener->endMetaAssembly( _assembledContigs.size() );
}



bool ExtendCycle::cycleHasRun(){ return _cycleHasRun; }
void ExtendCycle::getContigs(set<ScoredSeq*>* outputContigs){
  outputContigs->insert(_assembledContigs.begin(), _assembledContigs.end());
}


#endif
