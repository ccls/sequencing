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

#ifndef PARAMSCYCLEMANAGEMENT_CPP
#define PARAMSCYCLEMANAGEMENT_CPP

#include "ParamsCycleManagement.h"
#include "AssemblyException.h"
#include "ScoredSeqSubseq.h"
#include <iostream>
using namespace::std;



EcoFilter* ParamsCycleManagement::defaultCountEcoFilter(){ return new EcoFilterMinCount(1.5); }
long ParamsCycleManagement::defaultMaxContigLinkages(){ return 2; }

// constructors
ParamsCycleManagement::ParamsCycleManagement(){}

ParamsCycleManagement::ParamsCycleManagement(int totalCycles) :
  _totalCycles(totalCycles),
  _outputInterval(1)
{
  _maxContigLinkages = defaultMaxContigLinkages();
  setupFilters();
  setupTrimming();
  setupRepeatDetection();
}
/*
ParamsCycleManagement::ParamsCycleManagement(int totalCycles, long maxContigLinkages) :
  _totalCycles(totalCycles),
  _maxContigLinkages(maxContigLinkages),
  _outputInterval(1)
{
  setupFilters();
  setupTrimming();
  setupRepeatDetection();
}
*/


void ParamsCycleManagement::setupFilters(){
  _sizeFilter = new EcoFilterMinLength();
  _otherFilters = new EcoFilter*[1];
  _numOtherFilters = 0;
}
void ParamsCycleManagement::setupRepeatDetection(){
  // before-cycle data structures
  _repMinSizeBefore = new long[ _totalCycles+1 ];
  _repStdDevsBefore = new float[ _totalCycles+1 ];
  _repMinFoldIncBefore = new float[ _totalCycles+1 ];
  _repMatchFractIdBefore = new float[ _totalCycles+1 ];
  _repMakeFileBefore = new bool[ _totalCycles+1 ];
  _repFileNameBefore = new string[ _totalCycles+1 ];
  // after-cycle data structures
  _repMinSizeAfter = new long[ _totalCycles+1 ];
  _repStdDevsAfter = new float[ _totalCycles+1 ];
  _repMinFoldIncAfter = new float[ _totalCycles+1 ];
  _repMatchFractIdAfter = new float[ _totalCycles+1 ];
  _repMakeFileAfter = new bool[ _totalCycles+1 ];
  _repFileNameAfter = new string[ _totalCycles+1 ];
}
// initially null
void ParamsCycleManagement::setupTrimming(){
  _cycleToTrimMinScore = new float[_totalCycles];
  _cycleToTrimMinLength = new long[_totalCycles];
  _isTrimTransition = new bool[_totalCycles];
  _wasTrimSetSpecifically = new bool[_totalCycles];
  // give all the cycles null values
  for (long n = 0; n < _totalCycles; ++n){
    _cycleToTrimMinScore[n] = 0;
    _cycleToTrimMinLength[n] = 0;
    _isTrimTransition[n] = false;
    _wasTrimSetSpecifically[n] = false;
  }
  // ...and the initial (input) trimming values
  _initTrimMinScore = 0;
  _initTrimMinLength = 0;
}
// copied from another instance
void ParamsCycleManagement::setupTrimming(ParamsCycleManagement* pcm){
  _cycleToTrimMinScore = new float[_totalCycles];
  _cycleToTrimMinLength = new long[_totalCycles];
  _isTrimTransition = new bool[_totalCycles];
  _wasTrimSetSpecifically = new bool[_totalCycles];
  // copy the values
  for (long n = 0; n < _totalCycles; ++n){
    _cycleToTrimMinScore[n] = pcm->_cycleToTrimMinScore[n];
    _cycleToTrimMinLength[n] = pcm->_cycleToTrimMinLength[n];
    _isTrimTransition[n] = pcm->_isTrimTransition[n];
    _wasTrimSetSpecifically[n] = pcm->_wasTrimSetSpecifically[n];
  }
  // ...and the initial (input) trimming values
  _initTrimMinScore = pcm->_initTrimMinScore;
  _initTrimMinLength = pcm->_initTrimMinLength;
}


ParamsCycleManagement::~ParamsCycleManagement(){
  delete [] _isTrimTransition;
  delete [] _wasTrimSetSpecifically;
  delete [] _cycleToTrimMinScore;
  delete [] _cycleToTrimMinLength;
  delete _sizeFilter;
  for (int n = 0; n < _numOtherFilters; ++n){ delete _otherFilters[n]; }
  delete [] _otherFilters;

  delete [] _repMinSizeBefore;
  delete [] _repStdDevsBefore;
  delete [] _repMinFoldIncBefore;
  delete [] _repMatchFractIdBefore;
  delete [] _repMakeFileBefore;
  delete [] _repFileNameBefore;

  delete [] _repMinSizeAfter;
  delete [] _repStdDevsAfter;
  delete [] _repMinFoldIncAfter;
  delete [] _repMatchFractIdAfter;
  delete [] _repMakeFileAfter;
  delete [] _repFileNameAfter;
}


ParamsCycleManagement::ParamsCycleManagement(ParamsCycleManagement* pcm) :
  _totalCycles(pcm->_totalCycles),
  _maxContigLinkages(pcm->_maxContigLinkages)
{
  _sizeFilter = pcm->_sizeFilter->copy();
  _numOtherFilters = pcm->_numOtherFilters;
  _otherFilters = new EcoFilter*[ _numOtherFilters + 1 ];
  _resetCycles.insert(pcm->_resetCycles.begin(), pcm->_resetCycles.end());
  for (int n = 0; n < _numOtherFilters; ++n){
    _otherFilters[n] = pcm->_otherFilters[n]->copy();
  }
  setupTrimming(pcm);
  _outputInterval = pcm->_outputInterval;

  _seekRepsBeforeCycles.insert(pcm->_seekRepsBeforeCycles.begin(), pcm->_seekRepsBeforeCycles.end());
  _seekRepsAfterCycles.insert(pcm->_seekRepsAfterCycles.begin(), pcm->_seekRepsAfterCycles.end());
  setupRepeatDetection();
  for (set<int>::iterator it = _seekRepsBeforeCycles.begin(); it !=_seekRepsBeforeCycles.end(); ++it){
    int n = *it - 1;
    _repMinSizeBefore[n] = pcm->_repMinSizeBefore[n];
    _repStdDevsBefore[n] = pcm->_repStdDevsBefore[n];
    _repMinFoldIncBefore[n] = pcm->_repMinFoldIncBefore[n];
    _repMatchFractIdBefore[n] = pcm->_repMatchFractIdBefore[n];
    _repMakeFileBefore[n] = pcm->_repMakeFileBefore[n];
    if ( _repMakeFileBefore[n] ){ _repFileNameBefore[n] = pcm->_repFileNameBefore[n]; }
  }
  for (set<int>::iterator it = _seekRepsAfterCycles.begin(); it !=_seekRepsAfterCycles.end(); ++it){
    int n = *it - 1;
    _repMinSizeAfter[n] = pcm->_repMinSizeAfter[n];
    _repStdDevsAfter[n] = pcm->_repStdDevsAfter[n];
    _repMinFoldIncAfter[n] = pcm->_repMinFoldIncAfter[n];
    _repMatchFractIdAfter[n] = pcm->_repMatchFractIdAfter[n];
    _repMakeFileAfter[n] = pcm->_repMakeFileAfter[n];
    if ( _repMakeFileAfter[n] ){ _repFileNameAfter[n] = pcm->_repFileNameAfter[n]; }
  }
}
ParamsCycleManagement* ParamsCycleManagement::copy(){ return new ParamsCycleManagement(this); }



bool ParamsCycleManagement::shouldContigsReset(int cycleNum){
  if (cycleNum < 1){ throw AssemblyException::ArgError("PCM:contig history reset cycle num must be >0 (indexed from 1)"); }
  return (_resetCycles.count(cycleNum) != 0);
}
void ParamsCycleManagement::addContigReset(int cycleNum){
  if (cycleNum < 1){ throw AssemblyException::ArgError("PCM:contig history reset cycle num must be >0 (indexed from 1)"); }
  _resetCycles.insert(cycleNum);
}


void ParamsCycleManagement::addRepeatDetection(int cycleNum, float numStdDevs, float minFoldIncrease,
					       WhenInCycle whenInCycle, long minRepeatSize, float matchFractId){
  if (cycleNum < 1 or cycleNum > _totalCycles){
    throw AssemblyException::ArgError("PCM:seek repeats add, cycle num must be >0 and not greater than the num of cycles");
  }
  if (whenInCycle == BEFORECYCLE){
    if( _seekRepsBeforeCycles.count(cycleNum) != 0 ){
      throw AssemblyException::ArgError("PCM::addRepDetect cannot have redundant detection instances");
    }
    _seekRepsBeforeCycles.insert(cycleNum);
    _repMinSizeBefore[cycleNum-1] = minRepeatSize;
    _repStdDevsBefore[cycleNum-1] = numStdDevs;
    _repMinFoldIncBefore[cycleNum-1] = minFoldIncrease;
    _repMatchFractIdBefore[cycleNum-1] = matchFractId;
    _repMakeFileBefore[cycleNum-1] = false;
  } else if (whenInCycle == AFTERCYCLE){
    if( _seekRepsAfterCycles.count(cycleNum) != 0 ){
      throw AssemblyException::ArgError("PCM::addRepDetect cannot have redundant detection instances");
    }
    _seekRepsAfterCycles.insert(cycleNum);
    _repMinSizeAfter[cycleNum-1] = minRepeatSize;
    _repStdDevsAfter[cycleNum-1] = numStdDevs;
    _repMinFoldIncAfter[cycleNum-1] = minFoldIncrease;
    _repMatchFractIdAfter[cycleNum-1] = matchFractId;
    _repMakeFileAfter[cycleNum-1] = false;
  } else { throw AssemblyException::ArgError("PCM:addRepeatDetection bad whenInCycle arg"); }
}
void ParamsCycleManagement::addRepeatDetection(int cycleNum, float numStdDevs, float minFoldIncrease,
					       WhenInCycle whenInCycle, long minRepeatSize, float matchFractId, string outfileName){
  // this is just called as a fail-fast way of checking that the file name will be usable
  ReadFile::FileType fileType = ReadFile::getFileType(outfileName);
  if (fileType != ReadFile::FASTAFILE and fileType != ReadFile::PRICEQFILE){
    throw AssemblyException::ArgError("PCM::addRepDetect, filename was not of a usable type");
  }

  // now fill in the values
  if (cycleNum < 1 or cycleNum > _totalCycles){
    throw AssemblyException::ArgError("PCM:seek repeats add, cycle num must be >0 and not greater than the num of cycles");
  }
  if (whenInCycle == BEFORECYCLE){
    if( _seekRepsBeforeCycles.count(cycleNum) != 0 ){
      throw AssemblyException::ArgError("PCM::addRepDetect cannot have redundant detection instances");
    }
    _seekRepsBeforeCycles.insert(cycleNum);
    _repMinSizeBefore[cycleNum-1] = minRepeatSize;
    _repStdDevsBefore[cycleNum-1] = numStdDevs;
    _repMinFoldIncBefore[cycleNum-1] = minFoldIncrease;
    _repMatchFractIdBefore[cycleNum-1] = matchFractId;
    _repMakeFileBefore[cycleNum-1] = true;
    _repFileNameBefore[cycleNum-1] = outfileName;
  } else if (whenInCycle == AFTERCYCLE){
    if( _seekRepsAfterCycles.count(cycleNum) != 0 ){
      throw AssemblyException::ArgError("PCM::addRepDetect cannot have redundant detection instances");
    }
    _seekRepsAfterCycles.insert(cycleNum);
    _repMinSizeAfter[cycleNum-1] = minRepeatSize;
    _repStdDevsAfter[cycleNum-1] = numStdDevs;
    _repMinFoldIncAfter[cycleNum-1] = minFoldIncrease;
    _repMatchFractIdAfter[cycleNum-1] = matchFractId;
    _repMakeFileAfter[cycleNum-1] = true;
    _repFileNameAfter[cycleNum-1] = outfileName;
  } else { throw AssemblyException::ArgError("PCM:addRepeatDetection bad whenInCycle arg"); }
}
bool ParamsCycleManagement::seekDynamicRepeats(int cycleNum, WhenInCycle whenInCycle){
  if (cycleNum < 1 or cycleNum > _totalCycles){
    throw AssemblyException::ArgError("PCM:seek repeats, cycle num must be >0 and not greater than the num of cycles");
  }
  if (whenInCycle == BEFORECYCLE){
    return (_seekRepsBeforeCycles.count(cycleNum) != 0);
  } else if (whenInCycle == AFTERCYCLE){
    return (_seekRepsAfterCycles.count(cycleNum) != 0);
  } else {
    throw AssemblyException::ArgError("PCM:seekDynamicRepeats bad whenInCycle arg");
  }
}
long ParamsCycleManagement::repeatMinSize(int cycleNum, WhenInCycle whenInCycle){
  if (! seekDynamicRepeats(cycleNum, whenInCycle)){
    throw AssemblyException::ArgError("PCM:repMinSize requires seekDynamicRepeats is TRUE");
  }
  if (whenInCycle == BEFORECYCLE){ return _repMinSizeBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repMinSizeAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:repMinSize bad whenInCycle arg"); }
}
float ParamsCycleManagement::repeatStandDevs(int cycleNum, WhenInCycle whenInCycle){
  if (! seekDynamicRepeats(cycleNum, whenInCycle)){
    throw AssemblyException::ArgError("PCM:repStandDevs requires seekDynamicRepeats is TRUE");
  }
  if (whenInCycle == BEFORECYCLE){ return _repStdDevsBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repStdDevsAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:repStandDevs bad whenInCycle arg"); }
}
float ParamsCycleManagement::repeatMinFoldIncrease(int cycleNum, WhenInCycle whenInCycle){
  if (! seekDynamicRepeats(cycleNum, whenInCycle)){
    throw AssemblyException::ArgError("PCM:repMinFold requires seekDynamicRepeats is TRUE");
  }
  if (whenInCycle == BEFORECYCLE){ return _repMinFoldIncBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repMinFoldIncAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:repMinFold bad whenInCycle arg"); }
}
float ParamsCycleManagement::repeatMatchFractId(int cycleNum, WhenInCycle whenInCycle){
  if (! seekDynamicRepeats(cycleNum, whenInCycle)){
    throw AssemblyException::ArgError("PCM:repMinFold requires seekDynamicRepeats is TRUE");
  }
  if (whenInCycle == BEFORECYCLE){ return _repMatchFractIdBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repMatchFractIdAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:repMatchFractId bad whenInCycle arg"); }
}
bool ParamsCycleManagement::makeRepeatOutfile(int cycleNum, WhenInCycle whenInCycle){
  if (! seekDynamicRepeats(cycleNum, whenInCycle)){ return false; }
  else if (whenInCycle == BEFORECYCLE){ return _repMakeFileBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repMakeFileAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:makeRepOutfile bad whenInCycle arg"); }
}
string ParamsCycleManagement::repeatOutfileName(int cycleNum, WhenInCycle whenInCycle){
  if (! makeRepeatOutfile(cycleNum, whenInCycle)){
    throw AssemblyException::ArgError("PCM:repOutfileName, requires makeRepeatOutfile is TRUE");
  }
  if (whenInCycle == BEFORECYCLE){ return _repFileNameBefore[cycleNum-1]; }
  else if (whenInCycle == AFTERCYCLE){ return _repFileNameAfter[cycleNum-1]; }
  else { throw AssemblyException::ArgError("PCM:repOutfileName bad whenInCycle arg"); }
}


void ParamsCycleManagement::addSizeFilter(long minLength, long skipCycles){
  _sizeFilter->addMinLength(minLength,skipCycles);
}
void ParamsCycleManagement::resetSizeFilter(){
  delete _sizeFilter;
  _sizeFilter = new EcoFilterMinLength();
}
void ParamsCycleManagement::setSizeFilter(EcoFilterMinLength* sizeFilter){
  delete _sizeFilter;
  _sizeFilter = sizeFilter->copy();
}
EcoFilterMinLength* ParamsCycleManagement::getSizeFilter(){ return _sizeFilter->copy(); }



void ParamsCycleManagement::addEcoFilter(EcoFilter* ecoFilter){
  EcoFilter** filterHolder = new EcoFilter*[ _numOtherFilters + 2 ];
  for (int n = 0; n < _numOtherFilters; ++n){ filterHolder[n] = _otherFilters[n]; }
  filterHolder[ _numOtherFilters ] = ecoFilter->copy();
  delete [] _otherFilters;
  _otherFilters = filterHolder;
  _numOtherFilters++;
}
EcoFilter* ParamsCycleManagement::getFullEcoFilter(){
  _otherFilters[ _numOtherFilters ] = _sizeFilter;
  return new EcoFilterMulti(_otherFilters, _numOtherFilters + 1);
}




// get the values
int ParamsCycleManagement::totalCycles(){ return _totalCycles; }

// deal with maxContigLinkages
long ParamsCycleManagement::maxContigLinkages(){ return _maxContigLinkages; }
void ParamsCycleManagement::setMaxContigLinkages(long mcl){ _maxContigLinkages = mcl; }


void ParamsCycleManagement::addTrim(int cycleNum, float coverageLevel){
  // just add in the null value for min contig length
  addTrim(cycleNum, coverageLevel, 0);
}
void ParamsCycleManagement::addTrim(int cycleNum, float coverageLevel, long minLength){
  if (_totalCycles < cycleNum){ throw AssemblyException::ArgError("PCM cannot add trim to a cycle beyond the num of cycles (1)."); }
  if (cycleNum < 1){ throw AssemblyException::ArgError("PCM cannot add trim to a zero cycle"); }
  cycleNum = cycleNum - 1;
  _cycleToTrimMinScore[cycleNum] = coverageLevel;
  _cycleToTrimMinLength[cycleNum] = minLength;
  _wasTrimSetSpecifically[cycleNum] = true;
}

void ParamsCycleManagement::addBasalTrim(int numPriorCycles, float coverageLevel){
  // just add in the null value for min contig length
  addBasalTrim(numPriorCycles, coverageLevel, 0);
}
void ParamsCycleManagement::addBasalTrim(int numPriorCycles, float coverageLevel, long minLength){
  if (_totalCycles <= numPriorCycles){
    throw AssemblyException::ArgError("PCM cannot add trim to a cycle beyond the num of cycles (1).");
  }
  int n = numPriorCycles;
  while (n < _totalCycles and (! _isTrimTransition[n])){
    if (! _wasTrimSetSpecifically[n]){
      _cycleToTrimMinScore[n] = coverageLevel;
      _cycleToTrimMinLength[n] = minLength;
    }
    ++n;
  }
  // this is true even if there is a conflict with a specifically-set cycle
  _isTrimTransition[numPriorCycles] = true;
}
void ParamsCycleManagement::addInitialTrim(float coverageLevel){
    _initTrimMinScore = coverageLevel;
    _initTrimMinLength = 0;
}
void ParamsCycleManagement::addInitialTrim(float coverageLevel, long minLength){
    _initTrimMinScore = coverageLevel;
    _initTrimMinLength = minLength;
}

// methods for calling trim
void ParamsCycleManagement::executeTrim(set<ScoredSeq*>* contigSet, int cycleNum){
  set<ScoredSeq*> discardSet;
  executeTrim(contigSet, cycleNum, &discardSet);
  for (set<ScoredSeq*>::iterator it = discardSet.begin(); it != discardSet.end(); ++it){ (*it)->deepDelete(); }
}
void ParamsCycleManagement::executeTrim(set<ScoredSeq*>* contigSet, int cycleNum, set<ScoredSeq*>* discardSet){
  // if the cycle num is above the max, use the max's values
  if (cycleNum > _totalCycles){ throw AssemblyException::ArgError("PCM::executeTrim cannot trim for a cycle beyond the max cycle"); }
  if (cycleNum < 1){ throw AssemblyException::ArgError("PCM::executeTrim cannot trim for a zero cycle"); }
  cycleNum = cycleNum - 1;
  if ( _cycleToTrimMinScore[cycleNum] > 0 ){
    executeTrimHelper(contigSet, discardSet, _cycleToTrimMinScore[cycleNum], _cycleToTrimMinLength[cycleNum]);
  }
}
void ParamsCycleManagement::executeInitialTrim(set<ScoredSeq*>* contigSet){
  set<ScoredSeq*> discardSet;
  executeInitialTrim(contigSet, &discardSet);
  for (set<ScoredSeq*>::iterator it = discardSet.begin(); it != discardSet.end(); ++it){ (*it)->deepDelete(); }
}
void ParamsCycleManagement::executeInitialTrim(set<ScoredSeq*>* contigSet, set<ScoredSeq*>* discardSet){
  if ( _initTrimMinScore > 0 ){ executeTrimHelper(contigSet, discardSet, _initTrimMinScore, _initTrimMinLength); }
}

// support methods for executing trim
void ParamsCycleManagement::executeTrimHelper(set<ScoredSeq*>* contigSet, set<ScoredSeq*>* discardSet, float minScore, long minLength){
  set<ScoredSeq*> kept;
  for (set<ScoredSeq*>::iterator it = contigSet->begin(); it != contigSet->end(); ++it){
    ScoredSeq* input = (*it);
    ScoredSeq* trimmed = trimThisContig( input, minScore, minLength );
    if (input == trimmed){ kept.insert( input ); }
    else {
      discardSet->insert( input );
      if ( trimmed != NULL ){ kept.insert( trimmed ); }
    }
  }
  contigSet->clear();
  contigSet->insert( kept.begin(), kept.end() );
}
ScoredSeq* ParamsCycleManagement::trimThisContig(ScoredSeq* old, float minScore, long minLength){
  // a contig cannot be null
  if (minLength < 1){ minLength = 1; }
  long contigSize = old->size();
  long startPos = 0;
  long endPos = contigSize - 1;
  while (startPos < endPos and old->scoreAtPosPlus(startPos) < minScore){ startPos++; }
  while (endPos >= startPos and old->scoreAtPosPlus(endPos) < minScore){ endPos--; }
  // a null will be returned if the resulting contig isn't long enough
  if (endPos - startPos + 1 < minLength){ return NULL; }
  // no trim, so return the original
  else if (startPos == 0 and endPos == contigSize - 1){ return old; }
  else {
    ScoredSeq* tempSeq = new ScoredSeqSubseq(old, startPos, endPos - startPos + 1);
    ScoredSeq* returnSeq = tempSeq->shallowCopy();
    delete tempSeq;
    return returnSeq;
  }
}


void ParamsCycleManagement::setOutputInterval(int interval){
  if (interval < 1){ throw AssemblyException::ArgError("PCM output interval needs to be >= 1."); }
  _outputInterval = interval;
}
int ParamsCycleManagement::getOutputInterval(){ return _outputInterval; }
bool ParamsCycleManagement::outputThisCycle(int cycleNum){
  if (cycleNum < 1){ throw AssemblyException::ArgError("PCM::outputThisCycle, there is no cycle less than 1."); }
  return (cycleNum % _outputInterval == 0);
}

#endif
