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



#ifndef ASSEMBLYJOBFACTORY_CPP
#define ASSEMBLYJOBFACTORY_CPP

#include "AssemblyJobFactory.h"

#include "AssemblyException.h"
#include <cmath>
#include <iostream>
using namespace::std;



AssemblyJobFactory::AssemblyJobFactory(){}
AssemblyJobFactory::AssemblyJobFactory(ParamsAlignment* pa,
				       ParamsDeBruijn* pdb,
				       AssemblyJob::AssemblyStrandedness strandedness,
				       AssemblerListener* listener) :
  _listener(listener),
  _strandedness(strandedness),
  _maxOverlap(-1){
  _alScoreMatrix = pa->getAsm();
  _paramsMinOverlap = new ParamsMinOverlap( dynamic_cast<ParamsMinOverlap*>(pa) );
  _paramsMinFractId = new ParamsMinFractId( dynamic_cast<ParamsMinFractId*>(pa) );
  _paramsDeBruijn = new ParamsDeBruijn(pdb);
}
AssemblyJobFactory::AssemblyJobFactory(ParamsAlignment* pa,
				       ParamsDeBruijn* pdb,
				       long maxOverlap,
				       AssemblyJob::AssemblyStrandedness strandedness,
				       AssemblerListener* listener) :
  _listener(listener),
  _strandedness(strandedness),
  _maxOverlap(maxOverlap){
  _alScoreMatrix = pa->getAsm();
  _paramsMinOverlap = new ParamsMinOverlap( dynamic_cast<ParamsMinOverlap*>(pa) );
  _paramsMinFractId = new ParamsMinFractId( dynamic_cast<ParamsMinFractId*>(pa) );
  _paramsDeBruijn = new ParamsDeBruijn(pdb);
}
AssemblyJobFactory::AssemblyJobFactory(ParamsMinOverlap * pmo,
				       ParamsMinFractId * pmf,
				       AlignmentScoreMatrix * sm,
				       ParamsDeBruijn* pdb,
				       AssemblyJob::AssemblyStrandedness strandedness,
				       AssemblerListener* listener) :
  _listener(listener),
  _strandedness(strandedness),
  _maxOverlap(-1){
  _alScoreMatrix = new AlignmentScoreMatrix( sm );
  _paramsMinOverlap = new ParamsMinOverlap( pmo );
  _paramsMinFractId = new ParamsMinFractId( pmf );
  _paramsDeBruijn = new ParamsDeBruijn(pdb);
}
AssemblyJobFactory::AssemblyJobFactory(ParamsMinOverlap * pmo,
				       ParamsMinFractId * pmf,
				       AlignmentScoreMatrix * sm,
				       ParamsDeBruijn* pdb,
				       long maxOverlap,
				       AssemblyJob::AssemblyStrandedness strandedness,
				       AssemblerListener* listener) :
  _listener(listener),
  _strandedness(strandedness),
  _maxOverlap(maxOverlap){
  _alScoreMatrix = new AlignmentScoreMatrix( sm );
  _paramsMinOverlap = new ParamsMinOverlap( pmo );
  _paramsMinFractId = new ParamsMinFractId( pmf );
  _paramsDeBruijn = new ParamsDeBruijn(pdb);
}

AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf){
  _listener = ajf->_listener;
  _strandedness = ajf->_strandedness;
  _maxOverlap = ajf->_maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( ajf->_alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( ajf->_paramsMinOverlap );
  _paramsMinFractId = new ParamsMinFractId( ajf->_paramsMinFractId );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}
AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf, ParamsMinOverlap* pmo){
  _listener = ajf->_listener;
  _strandedness = ajf->_strandedness;
  _maxOverlap = ajf->_maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( ajf->_alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( pmo );
  _paramsMinFractId = new ParamsMinFractId( ajf->_paramsMinFractId );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}
AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf, ParamsMinFractId* pmf){
  _listener = ajf->_listener;
  _strandedness = ajf->_strandedness;
  _maxOverlap = ajf->_maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( ajf->_alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( ajf->_paramsMinOverlap );
  _paramsMinFractId = new ParamsMinFractId( pmf );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}
AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf, AlignmentScoreMatrix* alScoreMatrix){
  _listener = ajf->_listener;
  _strandedness = ajf->_strandedness;
  _maxOverlap = ajf->_maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( ajf->_paramsMinOverlap );
  _paramsMinFractId = new ParamsMinFractId( ajf->_paramsMinFractId );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}
AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf, AssemblyJob::AssemblyStrandedness strandedness){
  _listener = ajf->_listener;
  _strandedness = strandedness;
  _maxOverlap = ajf->_maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( ajf->_alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( ajf->_paramsMinOverlap );
  _paramsMinFractId = new ParamsMinFractId( ajf->_paramsMinFractId );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}
AssemblyJobFactory::AssemblyJobFactory(AssemblyJobFactory* ajf, long maxOverlap){
  _listener = ajf->_listener;
  _strandedness = ajf->_strandedness;
  _maxOverlap = maxOverlap;
  _alScoreMatrix = new AlignmentScoreMatrix( ajf->_alScoreMatrix );
  _paramsMinOverlap = new ParamsMinOverlap( ajf->_paramsMinOverlap );
  _paramsMinFractId = new ParamsMinFractId( ajf->_paramsMinFractId );
  _paramsDeBruijn = new ParamsDeBruijn( ajf->_paramsDeBruijn );
}


AssemblyJobFactory::~AssemblyJobFactory(){
  delete _alScoreMatrix;
  delete _paramsMinOverlap;
  delete _paramsMinFractId;
  delete _paramsDeBruijn;
}



AlignmentScoreMatrix * AssemblyJobFactory::getScoreMatrix(){
  return new AlignmentScoreMatrix(_alScoreMatrix); }
ParamsMinOverlap * AssemblyJobFactory::getParamsMinOverlap(){
  return new ParamsMinOverlap(_paramsMinOverlap); }
ParamsMinFractId * AssemblyJobFactory::getParamsMinFractId(){
  return new ParamsMinFractId(_paramsMinFractId); }
ParamsDeBruijn * AssemblyJobFactory::getParamsDeBruijn(){
  return new ParamsDeBruijn(_paramsDeBruijn); }
long AssemblyJobFactory::getMaxOverlap(){ return _maxOverlap; }
AssemblyJob::AssemblyStrandedness AssemblyJobFactory::strandedness(){ return _strandedness; }
AssemblerListener* AssemblyJobFactory::getListener(){ return _listener; }


AssemblyJob * AssemblyJobFactory::assemblyJob(set<ScoredSeq*>* inputSeqs){
  long effectiveSeqNum;
  if (_strandedness == AssemblyJob::DOUBLESTRANDED){ effectiveSeqNum = inputSeqs->size() * 2; }
  else { effectiveSeqNum = inputSeqs->size(); }
  return generalAssemblyJob(inputSeqs,effectiveSeqNum);
}
AssemblyJob * AssemblyJobFactory::assemblyJob(set<ScoredSeq*>* inputSeqs, float contigFactor){
  long effectiveSeqNum;
  if (_strandedness == AssemblyJob::DOUBLESTRANDED){ effectiveSeqNum = long(float(inputSeqs->size() * 2) * contigFactor); }
  else { effectiveSeqNum = long(float(inputSeqs->size()) * contigFactor); }
  return generalAssemblyJob(inputSeqs,effectiveSeqNum);
}

AssemblyJob * AssemblyJobFactory::generalAssemblyJob(set<ScoredSeq*>* inputSeqs, long inputSeqNum){
  if ( inputSeqs->size() < 2 ){ return new AssemblyJobNull( inputSeqs ); }
  else if (_maxOverlap == -1) {
    return new AssemblyJobHierarchy(inputSeqs, _paramsMinOverlap, _paramsMinFractId, _alScoreMatrix, _paramsDeBruijn,
				    _strandedness, _listener);
  } else {
    return new AssemblyJobHierarchy(inputSeqs, _paramsMinOverlap, _paramsMinFractId, _alScoreMatrix, _paramsDeBruijn,
				    _strandedness, _maxOverlap, _listener);
  }
}



AssemblyJob * AssemblyJobFactory::redundancyJob(set<ScoredSeq*>* inputSeqs){
  long effectiveSeqNum;
  if (_strandedness == AssemblyJob::DOUBLESTRANDED){ effectiveSeqNum = inputSeqs->size() * 2; }
  else { effectiveSeqNum = inputSeqs->size(); }
  return generalRedundancyJob(inputSeqs,effectiveSeqNum);
}
AssemblyJob * AssemblyJobFactory::redundancyJob(set<ScoredSeq*>* inputSeqs, float contigFactor){
  long effectiveSeqNum;
  if (_strandedness == AssemblyJob::DOUBLESTRANDED){ effectiveSeqNum = long(float(inputSeqs->size() * 2) * contigFactor); }
  else { effectiveSeqNum = long(float(inputSeqs->size()) * contigFactor); }
  return generalRedundancyJob(inputSeqs,effectiveSeqNum);
}

AssemblyJob * AssemblyJobFactory::generalRedundancyJob(set<ScoredSeq*>* inputSeqs, long inputSeqNum){
  if ( inputSeqs->size() < 2 ){ return new AssemblyJobNull( inputSeqs ); }
  else if (_maxOverlap == -1) {
    return new AssemblyJobHierarchy(inputSeqs, _paramsMinOverlap, _paramsMinFractId, _alScoreMatrix, _paramsDeBruijn,
				    _strandedness, _listener, AssemblyJob::REDUNDANTASSEMBLY);
  } else {
    return new AssemblyJobHierarchy(inputSeqs, _paramsMinOverlap, _paramsMinFractId, _alScoreMatrix, _paramsDeBruijn,
				    _strandedness, _maxOverlap, _listener, AssemblyJob::REDUNDANTASSEMBLY);
  }
}



#endif
