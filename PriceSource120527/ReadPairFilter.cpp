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

#ifndef READPAIRFILTER_CPP
#define READPAIRFILTER_CPP

#include "ReadPairFilter.h"
#include "AssemblyException.h"
#include "Alignment.h"
#include "ScoredSeqCollectionBwt.h"
#include "MatchSeqTest.h"
#include "AssemblyJobPerfect.h"


ReadPairFilter::~ReadPairFilter(){}


ReadPairFilterNull::ReadPairFilterNull(){}
ReadPairFilterNull::~ReadPairFilterNull(){}
// this is never doing anything, so it would never want to be applied
ReadPairFilterNull* ReadPairFilterNull::copy(){ return new ReadPairFilterNull(); }
bool ReadPairFilterNull::useThisCycle(int cycleNum){ return false; }
bool ReadPairFilterNull::applyToLength(long seqSize){ return false; }
bool ReadPairFilterNull::testKillsPair(){ return false; }
bool ReadPairFilterNull::isReadOk(ScoredSeq* read){ return true; }



ReadPairFilterHomoPolymer::ReadPairFilterHomoPolymer(long maxHomoPolymer) :
  _maxHomoPolymer(maxHomoPolymer),
  _allLengths(true)
{}
ReadPairFilterHomoPolymer::ReadPairFilterHomoPolymer(long maxHomoPolymer, long maxSeqLength) :
  _maxHomoPolymer(maxHomoPolymer),
  _allLengths(false),
  _maxLength(maxSeqLength)
{}
ReadPairFilterHomoPolymer::~ReadPairFilterHomoPolymer(){}
ReadPairFilterHomoPolymer* ReadPairFilterHomoPolymer::copy(){
  if (_allLengths){ return new ReadPairFilterHomoPolymer(_maxHomoPolymer); }
  else { return new ReadPairFilterHomoPolymer(_maxHomoPolymer, _maxLength); }
}
bool ReadPairFilterHomoPolymer::useThisCycle(int cycleNum){ return true; }
bool ReadPairFilterHomoPolymer::applyToLength(long seqSize){ return _allLengths or seqSize <= _maxLength; }
bool ReadPairFilterHomoPolymer::testKillsPair(){ return true; }
bool ReadPairFilterHomoPolymer::isReadOk(ScoredSeq* read){
  long rs = read->size();
  long pos = 0;
  while (pos < rs and read->nucAtPosPlus(pos)=='N' and pos <= _maxHomoPolymer){ pos++; }
  // considers Ns to extend the prior nuc
  long currentCount = pos;
  // considers Ns to extend the following nuc
  long futureCount = pos;
  // the highest count so far (will be continuously updated)
  long maxCount = pos;
  // this is literally the prior nucleotide identity (at pos 0, there was no ID, hence N)
  char priorNuc = 'N';
  // N as the prior nuc is only valid as an initial condition for this variable
  char priorLegitNuc = 'N';

  while (pos < rs and maxCount <= _maxHomoPolymer){
    char currentNuc = read->nucAtPosPlus(pos);
    if (currentNuc == 'N'){
      currentCount++;
      futureCount++;
    } else {
      if (currentNuc == priorLegitNuc){
	currentCount++;
	futureCount++;
      } else if (priorNuc == 'N'){
	currentCount = 1;
	futureCount++;
      } else {
	currentCount = 1;
	futureCount = 0;
      }
      priorLegitNuc = currentNuc;
    }
    priorNuc = currentNuc;
    pos++;
    if (currentCount > maxCount){ maxCount = currentCount; }
    if (futureCount > maxCount){ maxCount = futureCount; }
  }
  return (maxCount <= _maxHomoPolymer);
}



ReadPairFilterDinucleotide::ReadPairFilterDinucleotide(long maxDinucleotide) :
  _maxDinucleotide(maxDinucleotide),
  _allLengths(true)
{
  if (maxDinucleotide < 2){ throw AssemblyException::ArgError("RPFDinuc, maxDi may not be less than two"); }
}
ReadPairFilterDinucleotide::ReadPairFilterDinucleotide(long maxDinucleotide, long maxSeqLength) :
  _maxDinucleotide(maxDinucleotide),
  _allLengths(false),
  _maxLength(maxSeqLength)
{
  if (maxDinucleotide < 2){ throw AssemblyException::ArgError("RPFDinuc, maxDi may not be less than two"); }
}
ReadPairFilterDinucleotide::~ReadPairFilterDinucleotide(){}
ReadPairFilterDinucleotide* ReadPairFilterDinucleotide::copy(){
  if (_allLengths){ return new ReadPairFilterDinucleotide(_maxDinucleotide); }
  else { return new ReadPairFilterDinucleotide(_maxDinucleotide, _maxLength); }
}
bool ReadPairFilterDinucleotide::useThisCycle(int cycleNum){ return true; }
bool ReadPairFilterDinucleotide::applyToLength(long seqSize){ return _allLengths or seqSize <= _maxLength; }
bool ReadPairFilterDinucleotide::testKillsPair(){ return true; }
bool ReadPairFilterDinucleotide::isReadOk(ScoredSeq* read){
  long rs = read->size();

  if (rs < 3){ return true; }
  else {
    // this is like the hp filter, except there are 'even' and 'odd' values being tracked simultaneously
    // the baseline value is always 2
    long currentCount = 2;
    long futureCount = 2;
    long maxCount = 2;

    // go ahead and fill in the first two values, then start scanning
    char priorNuc[2] = {read->nucAtPosPlus(0),read->nucAtPosPlus(1)};
    char priorLegitNuc[2] = {priorNuc[0],priorNuc[1]};

    long pos = 2;
    while (pos < rs and maxCount <= _maxDinucleotide){
      // decides if even or odd is being evaluated
      long eo = pos % 2;
      char currentNuc = read->nucAtPosPlus(pos);
      if (currentNuc == 'N'){
	currentCount++;
	futureCount++;
      } else {
	if (currentNuc == priorLegitNuc[eo]){
	  currentCount++;
	  futureCount++;
	} else if (priorNuc[eo] == 'N'){
	  currentCount = 2;
	  futureCount++;
	} else {
	  currentCount = 2;
	  futureCount = 2;
	}
	priorLegitNuc[eo] = currentNuc;
      }
      priorNuc[eo] = currentNuc;
      pos++;
      if (currentCount > maxCount){ maxCount = currentCount; }
      if (futureCount > maxCount){ maxCount = futureCount; }
    }
    return (maxCount <= _maxDinucleotide);
  }
}




ReadPairFilterAvoidSeqs::ReadPairFilterAvoidSeqs(set<ScoredSeq*>* seqsToAvoid, float minFractId) :
  _allLengths(true)
{
  set<ScoredSeq*> seqCopies;
  for (set<ScoredSeq*>::iterator it = seqsToAvoid->begin(); it != seqsToAvoid->end(); ++it){
    seqCopies.insert( (*it)->shallowCopy() );
  }
  _seqsToAvoid = new ScoredSeqCollectionBwt(&seqCopies, minFractId);
}
ReadPairFilterAvoidSeqs::ReadPairFilterAvoidSeqs(set<ScoredSeq*>* seqsToAvoid, float minFractId, long maxSeqLength) :
  _allLengths(false),
  _maxLength(maxSeqLength)
{
  set<ScoredSeq*> seqCopies;
  for (set<ScoredSeq*>::iterator it = seqsToAvoid->begin(); it != seqsToAvoid->end(); ++it){
    seqCopies.insert( (*it)->shallowCopy() );
  }
  _seqsToAvoid = new ScoredSeqCollectionBwt(&seqCopies, minFractId);
}
ReadPairFilterAvoidSeqs::~ReadPairFilterAvoidSeqs(){
  set<ScoredSeq*> seqsToDelete;
  _seqsToAvoid->getSeqs( &seqsToDelete );
  for (set<ScoredSeq*>::iterator it = seqsToDelete.begin(); it != seqsToDelete.end(); ++it){
    (*it)->deepDelete();
  }
  delete _seqsToAvoid;
}
ReadPairFilterAvoidSeqs* ReadPairFilterAvoidSeqs::copy(){
  set<ScoredSeq*> seqsToCopy;
  _seqsToAvoid->getSeqs( &seqsToCopy );
  if (_allLengths){ return new ReadPairFilterAvoidSeqs(&seqsToCopy, _seqsToAvoid->getFractId()); }
  else { return new ReadPairFilterAvoidSeqs(&seqsToCopy, _seqsToAvoid->getFractId(), _maxLength); }
}
bool ReadPairFilterAvoidSeqs::useThisCycle(int cycleNum){ return true; }
bool ReadPairFilterAvoidSeqs::applyToLength(long seqSize){ return _allLengths or seqSize <= _maxLength; }
bool ReadPairFilterAvoidSeqs::testKillsPair(){ return false; }
bool ReadPairFilterAvoidSeqs::isReadOk(ScoredSeq* read){
  return (! _seqsToAvoid->hasMatch(read, read->size(), ScoredSeqCollectionBwt::softMinOvl) );
}





ReadPairFilterQuality::ReadPairFilterQuality(float maxFractBad, float scoreThreshold) :
  _maxFractBad(maxFractBad),
  _scoreThreshold(scoreThreshold),
  _useEveryCycle(true),
  _allLengths(true)
{}
ReadPairFilterQuality::ReadPairFilterQuality(float maxFractBad, float scoreThreshold, int numSkipCycles, int numRunCycles) :
  _maxFractBad(maxFractBad),
  _scoreThreshold(scoreThreshold),
  _useEveryCycle(false),
  _firstCycle(numSkipCycles),
  _lastCycleP1(numSkipCycles + numRunCycles),
  _allLengths(true)
{}
ReadPairFilterQuality::ReadPairFilterQuality(float maxFractBad, float scoreThreshold, long maxSeqLength) :
  _maxFractBad(maxFractBad),
  _scoreThreshold(scoreThreshold),
  _useEveryCycle(true),
  _allLengths(false),
  _maxLength(maxSeqLength)
{}
ReadPairFilterQuality::ReadPairFilterQuality(float maxFractBad, float scoreThreshold, int numSkipCycles, int numRunCycles, long maxSeqLength) :
  _maxFractBad(maxFractBad),
  _scoreThreshold(scoreThreshold),
  _useEveryCycle(false),
  _firstCycle(numSkipCycles),
  _lastCycleP1(numSkipCycles + numRunCycles),
  _allLengths(false),
  _maxLength(maxSeqLength)
{}
ReadPairFilterQuality::~ReadPairFilterQuality(){}
ReadPairFilterQuality* ReadPairFilterQuality::copy(){
  if (_useEveryCycle){
    if (_allLengths){ return new ReadPairFilterQuality(_maxFractBad, _scoreThreshold); }
    else { return new ReadPairFilterQuality(_maxFractBad, _scoreThreshold, _maxLength); }
  } else {
    if (_allLengths){ return new ReadPairFilterQuality(_maxFractBad, _scoreThreshold, _firstCycle, _lastCycleP1 - _firstCycle); }
    else { return new ReadPairFilterQuality(_maxFractBad, _scoreThreshold, _firstCycle, _lastCycleP1 - _firstCycle, _maxLength); }
  }
}
bool ReadPairFilterQuality::useThisCycle(int cycleNum){
  return ( _useEveryCycle or (cycleNum >= _firstCycle and cycleNum < _lastCycleP1));
}
bool ReadPairFilterQuality::applyToLength(long seqSize){ return _allLengths or seqSize <= _maxLength; }
bool ReadPairFilterQuality::testKillsPair(){ return true; }
bool ReadPairFilterQuality::isReadOk(ScoredSeq* read){
  long seqPos = read->size();
  long maxNumMis = long(_maxFractBad * float(seqPos));
  bool seqOk = true;
  while (seqPos > 0){
    --seqPos;
    if ( _scoreThreshold > read->scoreAtPosPlus(seqPos) ){
      if (maxNumMis == 0){
	seqOk = false;
	seqPos = 0;
      } else { --maxNumMis; }
    }
  }
  return seqOk;
}



ReadPairFilterUncalledBases::ReadPairFilterUncalledBases(float maxFractBad) :
  _maxFractBad(maxFractBad),
  _useEveryCycle(true),
  _allLengths(true)
{}
ReadPairFilterUncalledBases::ReadPairFilterUncalledBases(float maxFractBad, int numSkipCycles, int numRunCycles) :
  _maxFractBad(maxFractBad),
  _useEveryCycle(false),
  _firstCycle(numSkipCycles),
  _lastCycleP1(numSkipCycles + numRunCycles),
  _allLengths(true)
{}
ReadPairFilterUncalledBases::ReadPairFilterUncalledBases(float maxFractBad, long maxSeqLength) :
  _maxFractBad(maxFractBad),
  _useEveryCycle(true),
  _allLengths(false),
  _maxLength(maxSeqLength)
{}
ReadPairFilterUncalledBases::ReadPairFilterUncalledBases(float maxFractBad, int numSkipCycles, int numRunCycles, long maxSeqLength) :
  _maxFractBad(maxFractBad),
  _useEveryCycle(false),
  _firstCycle(numSkipCycles),
  _lastCycleP1(numSkipCycles + numRunCycles),
  _allLengths(false),
  _maxLength(maxSeqLength)
{}
ReadPairFilterUncalledBases::~ReadPairFilterUncalledBases(){}
ReadPairFilterUncalledBases* ReadPairFilterUncalledBases::copy(){
  if (_useEveryCycle){
    if (_allLengths){ return new ReadPairFilterUncalledBases(_maxFractBad); }
    else { return new ReadPairFilterUncalledBases(_maxFractBad, _maxLength); }
  } else {
    if (_allLengths){ return new ReadPairFilterUncalledBases(_maxFractBad, _firstCycle, _lastCycleP1 - _firstCycle); }
    else { return new ReadPairFilterUncalledBases(_maxFractBad, _firstCycle, _lastCycleP1 - _firstCycle, _maxLength); }
  }
}
bool ReadPairFilterUncalledBases::useThisCycle(int cycleNum){
  return ( _useEveryCycle or (cycleNum >= _firstCycle and cycleNum < _lastCycleP1));
}
bool ReadPairFilterUncalledBases::applyToLength(long seqSize){ return _allLengths or seqSize <= _maxLength; }
bool ReadPairFilterUncalledBases::testKillsPair(){ return true; }
bool ReadPairFilterUncalledBases::isReadOk(ScoredSeq* read){
  long seqPos = read->size();
  long maxNumMis = long(_maxFractBad * float(seqPos));
  bool seqOk = true;
  while (seqPos > 0){
    --seqPos;
    if ( 'N' == read->nucAtPosPlus(seqPos) ){
      if (maxNumMis == 0){
	seqOk = false;
	seqPos = 0;
      } else { --maxNumMis; }
    }
  }
  return seqOk;
}




ReadPairFilterMulti::ReadPairFilterMulti(ReadPairFilter** testArray, int numTests) : _numTests(numTests){
  _testArray = new ReadPairFilter*[ _numTests+1 ];
  for (int n = 0; n < _numTests; ++n){ _testArray[n] = testArray[n]->copy(); }
}
ReadPairFilterMulti::~ReadPairFilterMulti(){
  for (int testNum = 0; testNum < _numTests; ++testNum){ delete _testArray[testNum]; }
  delete [] _testArray;
}
ReadPairFilterMulti* ReadPairFilterMulti::copy(){ return new ReadPairFilterMulti(_testArray,_numTests); }
bool ReadPairFilterMulti::useThisCycle(int cycleNum){
  throw AssemblyException::CallingError("RPFMulti::useThisCycle, calling is inappropriate, could hide the preferences of contained RPFs.");
}
bool ReadPairFilterMulti::applyToLength(long seqSize){
  throw AssemblyException::CallingError("RPFMulti::applyToLength, calling is inappropriate, could hide the preferences of contained RPFs.");
}
bool ReadPairFilterMulti::testKillsPair(){
  throw AssemblyException::CallingError("RPFMulti::testKillsPair, calling is inappropriate, could hide the preferences of contained RPFs.");
}
bool ReadPairFilterMulti::isReadOk(ScoredSeq* read){
  bool testsPass = true;
  int testNum = 0;
  while (testsPass and testNum < _numTests){
    testsPass = _testArray[testNum]->isReadOk(read);
    testNum++;
  }
  return testsPass;
}


#endif

