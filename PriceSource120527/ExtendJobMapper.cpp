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

#ifndef EXTENDJOBMAPPER_CPP
#define EXTENDJOBMAPPER_CPP

#include "ExtendJobMapper.h"
#include <typeinfo>
#include <cstring>
#include <omp.h>
#include "ScoredSeqSubseq.h"
#include "ScoredSeqMonoScore.h"
using namespace::std;


ExtendJobMapper::~ExtendJobMapper(){}



//ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(){}

ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestNull();
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestNull();
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestEdgeProximity(edgeTest, fractMapId);
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers,
						   float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestEdgeProximity(edgeTest, fractMapId);
}


void ExtendJobMapperWithLinks::constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers){

  _shortArray = new ScoredSeq*[ contigsWithCarriers->size() * 2 + 1 ];
  _flipArray = new ScoredSeq*[ contigsWithCarriers->size() * 2 + 1 ];

  _plusIndex = ExtendCycle::plusIndex();
  _minusIndex = ExtendCycle::minusIndex();
  _frontIndex = ExtendCycle::frontIndex();
  _backIndex = ExtendCycle::backIndex();

  // use the contigsWithCarriers to construct the front and back BWT subseq collections
  // minContigSplitLen -> edgeSize

  // first the front edges
  set<ScoredSeq*> inputSet;
  _disposableIndex = 0;

  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeqFront;
    ScoredSeq* keySeqBack;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeqFront = new ScoredSeqSubseq( carrier, 0, _edgeSize);
      keySeqBack = new ScoredSeqSubseq( carrier, carrier->size() - _edgeSize, _edgeSize);
    } else {
      keySeqFront = new ScoredSeqSubseq( carrier, 0, carrier->size() );
      keySeqBack = new ScoredSeqSubseq( carrier, 0, carrier->size());
    }
    ScoredSeq* inSeqFront = ScoredSeqFlip::getFlip(keySeqFront,'-');
    ScoredSeq* inSeqBack = ScoredSeqFlip::getFlip(keySeqBack,'+');
    inputSet.insert(inSeqFront);
    inputSet.insert(inSeqBack);

    _shortArray[_disposableIndex] = keySeqFront;
    _flipArray[_disposableIndex] = inSeqFront;
    ++_disposableIndex;
    _shortArray[_disposableIndex] = keySeqBack;
    _flipArray[_disposableIndex] = inSeqBack;
    ++_disposableIndex;

    keySeqFront->buffer();
    keySeqBack->buffer();
  }
  _mappingCollection = new ScoredSeqCollectionBwt(&inputSet, _fractMapId, false);
  _mappingCollection->disableEdgeScaling();
}



ExtendJobMapperWithLinks::~ExtendJobMapperWithLinks(){
  // get rid of the subseq keys, not the contig values
  for (long n = 0; n < _disposableIndex; ++n){
    delete _shortArray[n];
    delete _flipArray[n];
  }
  delete [] _shortArray;
  delete [] _flipArray;

  // and delete the collections themselves
  delete _mappingCollection;
  delete _alTest;
}


long ExtendJobMapperWithLinks::mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
						bool* seekReadMatch, bool* seekPairMatch,
						bool* readMatchFound, bool* pairMatchFound, long numQueries){
  /*
  // simply forward the method - for now
  long hitCount = mapQueries(readQueries, seekReadMatch, readMatchFound, numQueries);
  hitCount += mapQueries(pairQueries, seekPairMatch, pairMatchFound, numQueries);
  return hitCount;
  */

  long numWithMatch = 0;
  // these will be filled in by the setup function
  vector<ScoredSeqWithCarriers*>** readMatchCarrier = new vector<ScoredSeqWithCarriers*>*[2];
  vector<ScoredSeqWithCarriers*>** pairMatchCarrier = new vector<ScoredSeqWithCarriers*>*[2];

  // do everything twice: once for reads, once for pairs
  // READS:
  long* qToUniqueReads = new long[numQueries+1];
  bool* uReadPasses = mapQueriesSetup(readQueries, seekReadMatch, numQueries, qToUniqueReads, readMatchCarrier);
  vector<ScoredSeqWithCarriers*>* readMatchContigsPlus = readMatchCarrier[0];
  vector<ScoredSeqWithCarriers*>* readMatchContigsMinus = readMatchCarrier[1];
  // PAIRS:
  long* qToUniquePairs = new long[numQueries+1];
  bool* uPairPasses = mapQueriesSetup(pairQueries, seekPairMatch, numQueries, qToUniquePairs, pairMatchCarrier);
  vector<ScoredSeqWithCarriers*>* pairMatchContigsPlus = pairMatchCarrier[0];
  vector<ScoredSeqWithCarriers*>* pairMatchContigsMinus = pairMatchCarrier[1];

  // now create links where appropriate
  for (long n = 0; n < numQueries; ++n){
    long urN = qToUniqueReads[n];
    long upN = qToUniquePairs[n];
    bool readPasses = seekReadMatch[n] and uReadPasses[urN];
    bool pairPasses = seekPairMatch[n] and uPairPasses[upN];
    if (readPasses or pairPasses){
      // this pair is ok, so fill in the values of reads, then pairs, as appropriate
      // READS:
      if ( seekReadMatch[n] ){
	removeReciprocalMatches(&readMatchContigsPlus[urN], &readMatchContigsMinus[urN]);
	if (readMatchContigsPlus[urN].size() + readMatchContigsMinus[urN].size() > 0){
	  addToContigMatchesHelper( readQueries[n], _plusIndex, _backIndex, &readMatchContigsPlus[urN]);
	  addToContigMatchesHelper( readQueries[n], _minusIndex, _frontIndex, &readMatchContigsMinus[urN]);
	  readMatchFound[n] = true;
	  ++numWithMatch;
	} else { readMatchFound[n] = false; }
      } else { readMatchFound[n] = false; }
      // PAIRS:
      if ( seekPairMatch[n] ){
	removeReciprocalMatches(&pairMatchContigsPlus[upN], &pairMatchContigsMinus[upN]);
	if (pairMatchContigsPlus[upN].size() + pairMatchContigsMinus[upN].size() > 0){
	  addToContigMatchesHelper( pairQueries[n], _plusIndex, _backIndex, &pairMatchContigsPlus[upN]);
	  addToContigMatchesHelper( pairQueries[n], _minusIndex, _frontIndex, &pairMatchContigsMinus[upN]);
	  pairMatchFound[n] = true;
	  ++numWithMatch;
	} else { pairMatchFound[n] = false; }
      } else { pairMatchFound[n] = false; }
    } else {
      readMatchFound[n] = false;
      pairMatchFound[n] = false;
    }
    readQueries[n]->deepUnbuffer();
    pairQueries[n]->deepUnbuffer();
  }

  delete [] readMatchCarrier;
  delete [] pairMatchCarrier;
  delete [] qToUniqueReads;
  delete [] qToUniquePairs;
  delete [] uReadPasses;
  delete [] uPairPasses;
  delete [] readMatchContigsPlus;
  delete [] readMatchContigsMinus;
  delete [] pairMatchContigsPlus;
  delete [] pairMatchContigsMinus;

  return numWithMatch;
}


bool* ExtendJobMapperWithLinks::mapQueriesSetup(ScoredSeqWithCarriers** queries, bool* seekMatch,
						long numQueries, long* qToUnique,
						vector<ScoredSeqWithCarriers*>** matchCarrier){

  // applies the test and sorts out the results but doesn't create actual linkages
  ScoredSeq** monoQueries = new ScoredSeq*[numQueries+1];
  vector<Alignment*>** matches = new vector<Alignment*>*[numQueries+1];
  long numUnique = mapQueriesHelper(queries, qToUnique, seekMatch, matches, monoQueries, numQueries);
  vector<ScoredSeqWithCarriers*>* matchContigsPlus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  vector<ScoredSeqWithCarriers*>* matchContigsMinus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  bool* uPass = new bool[numUnique+1];

  // i NEED to sort the matches even if the sequence doesn't pass the test because any
  // one of the reads with that sequences' paired ends could pass
  for (long n = 0; n < numUnique; ++n){
    uPass[n] = _alTest->passes( matches[n], monoQueries[n] );
    sortMatchesHelper(matches[n], &matchContigsPlus[n], &matchContigsMinus[n]);
    for (vector<Alignment*>::iterator it = matches[n]->begin(); it != matches[n]->end(); ++it){ delete *it; }
    delete matches[n];
    delete monoQueries[n];
  }
  delete [] matches;
  delete [] monoQueries;
  matchCarrier[0] = matchContigsPlus;
  matchCarrier[1] = matchContigsMinus;
  return uPass;
}


long ExtendJobMapperWithLinks::mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchesFound, long numQueries){
  vector<Alignment*>** matchesArray = new vector<Alignment*>*[numQueries+1];
  ScoredSeq** monoQueries = new ScoredSeq*[numQueries+1];
  long* qToUnique = new long[numQueries+1];
  long numUnique = mapQueriesHelper(queries, qToUnique, seekMatch, matchesArray, monoQueries, numQueries);
  long numThatPass = 0;

  vector<ScoredSeqWithCarriers*>* matchingContigsPlus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  vector<ScoredSeqWithCarriers*>* matchingContigsMinus = new vector<ScoredSeqWithCarriers*>[numUnique+1];

  bool* uniquePasses = new bool[numUnique+1];
  for (long n = 0; n < numUnique; ++n){
    uniquePasses[n] = _alTest->passes(matchesArray[n], monoQueries[n]);
    if (uniquePasses[n]){
      sortMatchesHelper(matchesArray[n], &matchingContigsPlus[n], &matchingContigsMinus[n]);
    }
    delete monoQueries[n];
  }
  delete [] monoQueries;

  for (long n = 0; n < numQueries; ++n){
    long uN = qToUnique[n];
    if ( seekMatch[n] and uniquePasses[qToUnique[n]] ){
      removeReciprocalMatches(&matchingContigsPlus[uN], &matchingContigsMinus[uN]);
      if (matchingContigsPlus[uN].size() + matchingContigsMinus[uN].size() > 0){
	addToContigMatchesHelper( queries[n], _plusIndex, _backIndex, &matchingContigsPlus[uN]);
	addToContigMatchesHelper( queries[n], _minusIndex, _frontIndex, &matchingContigsMinus[uN]);
	++numThatPass;
	matchesFound[n] = true;
      } else { matchesFound[n] = false; }
    } else { matchesFound[n] = false; }
    queries[n]->deepUnbuffer();
  }

  for (long n = 0; n < numUnique; ++n){
    for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
    delete matchesArray[n];
  }
  delete [] matchesArray;
  delete [] matchingContigsPlus;
  delete [] matchingContigsMinus;
  delete [] uniquePasses;
  delete [] qToUnique;

  return numThatPass;
}

/*
bool ExtendJobMapperWithLinks::passesTest(vector<Alignment*>* matches){
  return matches->size() > 0;
  //return true;
}
*/

ExtendJobMapperWithLinks::AlignmentTest::~AlignmentTest(){}

// the null test always returns true
ExtendJobMapperWithLinks::AlTestNull::AlTestNull(){}
ExtendJobMapperWithLinks::AlTestNull::~AlTestNull(){}
bool ExtendJobMapperWithLinks::AlTestNull::passes(vector<Alignment*>* matches, ScoredSeq* query){
  return matches->size() > 0;
}
ExtendJobMapperWithLinks::AlTestEdgeProximity::AlTestEdgeProximity(long edgeProximity, float fractId) :
  _edgeProximity(edgeProximity),
  _fractId(fractId){}
ExtendJobMapperWithLinks::AlTestEdgeProximity::~AlTestEdgeProximity(){}

bool ExtendJobMapperWithLinks::AlTestEdgeProximity::passes(vector<Alignment*>* matches, ScoredSeq* query){
  bool foundPass = false;
  // this allows the query to overhang the window to a degree acceptable by the percent ID,
  // while also including a -1 that compensates for the zero-indexing of the final aligned pos
  long effQsizeM1 = long(float(query->size()) * _fractId) - 1;
  long aLastPos = query->size() - 1;

  vector<Alignment*>::iterator alIt = matches->begin();
  vector<Alignment*>::iterator alEnd = matches->end();
  while (alIt != alEnd){
    Alignment* al = *alIt;
    long bLastPos;
    if ( al->isLinked(aLastPos, query) ){ bLastPos = al->getLinkage(aLastPos, query); }
    else { bLastPos = al->gapPairedAfter(aLastPos, query); }

    // test if the full sequence of the query would fit in the window given where the
    // alignment ends (it is possible that gaps push the 5p end of the query beyond the
    // window, but i don't want that to penalize the match).
    foundPass = _edgeProximity >= al->seqB()->size() - bLastPos + effQsizeM1;

    // if the test passed, skip to the end
    if (foundPass){ alIt = alEnd; }
    else { ++alIt; }
  }
  return foundPass;
}


// MODIFIES both input vectors
void ExtendJobMapperWithLinks::removeReciprocalMatches(vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						       vector<ScoredSeqWithCarriers*>* matchingContigsMinus){
  // the smaller of the two vectors will be made into a set (vector A); the
  // order-of-growth is NlogN for size(contigsA) and N for size(contigsB)
  vector<ScoredSeqWithCarriers*>* contigsA;
  vector<ScoredSeqWithCarriers*>* contigsB;
  if (matchingContigsPlus->size() < matchingContigsMinus->size()){
    contigsA = matchingContigsPlus;
    contigsB = matchingContigsMinus;
  } else {
    contigsB = matchingContigsPlus;
    contigsA = matchingContigsMinus;
  }

  // i don't need to do anything if the smaller vector is size zero
  if (contigsA->size() > 0){

    // this keeps the contigs that are ok from B; if not all of the contigs are being kept,
    // then the contents of contigsB can be replaced by the contents of this vector
    vector<ScoredSeqWithCarriers*> keepB;

    // this is the painful part - both of these operations require log entry/search times
    set<ScoredSeqWithCarriers*> nrSetA;
    for (vector<ScoredSeqWithCarriers*>::iterator it = contigsA->begin(); it != contigsA->end(); ++it){ nrSetA.insert(*it); }
    for (vector<ScoredSeqWithCarriers*>::iterator it = contigsB->begin(); it != contigsB->end(); ++it){
      set<ScoredSeqWithCarriers*>::iterator findIt = nrSetA.find(*it);
      if (findIt == nrSetA.end()){ keepB.push_back(*it); }
      else { nrSetA.erase( findIt ); }
    }

    // i don't have to do anything if contigsB and keepB are the same length
    if (keepB.size() < contigsB->size()){
      // replace the contents of B
      contigsB->clear();
      contigsB->insert(contigsB->begin(), keepB.begin(), keepB.end());
      // replace the contents of A (the non-redundant set can be used)
      contigsA->clear();
      for (set<ScoredSeqWithCarriers*>::iterator it = nrSetA.begin(); it != nrSetA.end(); ++it){ contigsA->push_back(*it); }
    }
  }
}


long ExtendJobMapperWithLinks::mapQueriesHelper(ScoredSeqWithCarriers** queries, long* qToUnique, bool* seekMatch,
						vector<Alignment*>** matchesArray, ScoredSeq** monoQueries, long numQueries){

  // make a non-redundant sequence set
  //map<string, QueryCarrier*> seqToQuerySet;
  map<string, vector<long>*> seqToQuerySet;
  long numNr = 0;
  for (long n = 0; n < numQueries; ++n){
    if (seekMatch[n]){
      char* querySeq = queries[n]->getSeq('+');
      map<string, vector<long>*>::iterator foundSeqIt = seqToQuerySet.find( querySeq );
      if (foundSeqIt == seqToQuerySet.end()){
	vector<long>* newCarrier = new vector<long>;
	newCarrier->push_back(n);
	seqToQuerySet.insert( pair<string, vector<long>*> (querySeq, newCarrier) );
      } else {
	foundSeqIt->second->push_back(n);
      }
      delete [] querySeq;
      queries[n]->deepUnbuffer();
    }
  }

  // set up the arrays
  long numMonoP1 = seqToQuerySet.size() + 1;
  vector<long>** likeSequences = new vector<long>*[ numMonoP1 ];
  //ScoredSeq** monoQueries = new ScoredSeq*[ numMonoP1 ];
  long* queryLengths = new long[ numMonoP1 ];
  // now use this as an index iterator - it will return to the full value
  long numMono = 0;
  for (map<string, vector<long>*>::iterator seqIt = seqToQuerySet.begin(); seqIt != seqToQuerySet.end(); ++seqIt){
    likeSequences[numMono] = seqIt->second;
    queryLengths[numMono] = seqIt->first.size();
    monoQueries[numMono] = new ScoredSeqMonoScore(seqIt->first, 1);
    matchesArray[numMono] = new vector<Alignment*>;
    ++numMono;
  }

  _mappingCollection->getBestMatches(numMono, matchesArray, monoQueries, '+', queryLengths,
				     ScoredSeqCollectionBwt::notThreaded, ScoredSeqCollectionBwt::softMinOvl);

  for (long mqN = 0; mqN < numMono; ++mqN){
    for (vector<long>::iterator it = likeSequences[mqN]->begin(); it != likeSequences[mqN]->end(); ++it){
      queries[*it]->deepUnbuffer();
      qToUnique[*it] = mqN;
    }
    //delete monoQueries[mqN];
    delete likeSequences[mqN];
  }
  //delete [] monoQueries;
  delete [] likeSequences;
  delete [] queryLengths;

  return numMono;
}


// MODIFIES: matchingContigs AND alreadyQueriedToHits
void ExtendJobMapperWithLinks::findMatchesHelper(ScoredSeq* currentQuery,
						 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						 vector<ScoredSeqWithCarriers*>* matchingContigsMinus){
  vector<Alignment*> matchAlignments;
  _mappingCollection->getBestMatches( &matchAlignments, currentQuery, '+', currentQuery->size(), ScoredSeqCollectionBwt::softMinOvl);

  for (vector<Alignment*>::iterator alIt = matchAlignments.begin(); alIt != matchAlignments.end(); ++alIt){
    // find the full contig
    ScoredSeqFlip* flipContig = dynamic_cast<ScoredSeqFlip*>( (*alIt)->seqB() );
    ScoredSeqSubseq* partialContig = dynamic_cast<ScoredSeqSubseq*>( flipContig->getNested() );
    ScoredSeqWithCarriers* fullContig = dynamic_cast<ScoredSeqWithCarriers*>( partialContig->getNested() );
    switch ( flipContig->getSense() ){
    case '+': matchingContigsPlus->push_back( fullContig ); break;
    case '-': matchingContigsMinus->push_back( fullContig ); break;
    default: throw AssemblyException::LogicError("EJMwLinks::findMatchesHelper got a bad sense char");
    }
    delete *alIt;
  }
}

// MODIFIES: matchingContigs AND alreadyQueriedToHits
void ExtendJobMapperWithLinks::sortMatchesHelper(vector<Alignment*>* matches,
						 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						 vector<ScoredSeqWithCarriers*>* matchingContigsMinus){

  for (vector<Alignment*>::iterator alIt = matches->begin(); alIt != matches->end(); ++alIt){
    // find the full contig
    ScoredSeqFlip* flipContig = dynamic_cast<ScoredSeqFlip*>( (*alIt)->seqB() );
    ScoredSeqSubseq* partialContig = dynamic_cast<ScoredSeqSubseq*>( flipContig->getNested() );
    ScoredSeqWithCarriers* fullContig = dynamic_cast<ScoredSeqWithCarriers*>( partialContig->getNested() );
    switch ( flipContig->getSense() ){
    case '+': matchingContigsPlus->push_back( fullContig ); break;
    case '-': matchingContigsMinus->push_back( fullContig ); break;
    default: throw AssemblyException::LogicError("EJMwLinks::findMatchesHelper got a bad sense char");
    }
  }
}


void ExtendJobMapperWithLinks::addToContigMatchesHelper(ScoredSeqWithCarriers* seqToAdd, int side, int frontOrBack,
							vector<ScoredSeqWithCarriers*>* matchingContigs){
  #pragma omp critical (EJMlink)
  {
    for (vector<ScoredSeqWithCarriers*>::iterator mcIt = matchingContigs->begin(); mcIt != matchingContigs->end(); ++mcIt){
      (*mcIt)->addCarriedSeq( seqToAdd, frontOrBack );
      seqToAdd->addCarriedSeq( (*mcIt), side );
    }
  }
}




ExtendJobMapperNoLinks::ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
}
ExtendJobMapperNoLinks::ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
}
void ExtendJobMapperNoLinks::constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers){

  // use the contigsWithCarriers to construct the front and back BWT subseq collection (one collection)
  set<ScoredSeq*> shortSet;

  // first the front edges
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeq;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeq = new ScoredSeqSubseq( carrier, 0, _edgeSize);
    } else {
      keySeq = new ScoredSeqSubseq( carrier, 0, carrier->size() );
    }
    shortSet.insert(keySeq);
    keySeq->buffer();
  }

  // now the back edges
  set<ScoredSeq*> shortBackSet;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeq;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeq = new ScoredSeqSubseq( carrier, carrier->size() - _edgeSize, _edgeSize);
    } else {
      keySeq = new ScoredSeqSubseq( carrier, 0, carrier->size());
    }
    shortSet.insert(keySeq);
    keySeq->buffer();
  }

  _mappingCollection = new ScoredSeqCollectionBwt(&shortSet, _fractMapId, false);
  _mappingCollection->disableEdgeScaling();

  // MORE HERE???
}

ExtendJobMapperNoLinks::~ExtendJobMapperNoLinks(){
  // get rid of the subseq keys, not the contig values
  set<ScoredSeq*> outerToDelete;
  _mappingCollection->getSeqs( &outerToDelete );
  for (set<ScoredSeq*>::iterator it = outerToDelete.begin(); it != outerToDelete.end(); ++it){ delete (*it); }
  delete _mappingCollection;
}


ExtendJobMapper::QueryCarrier::QueryCarrier(char* query, long size) : _query(query), _size(size) {}
ExtendJobMapper::QueryCarrier::~QueryCarrier(){}

/*
long ExtendJobMapperNoLinks::mapQueries(ScoredSeqWithCarriers** queries, long numQueries){
  bool* seekMatch = new bool[numQueries+1];
  for (long n = 0; n < numQueries; ++n){ seekMatch[n] = true; }
  bool* matchFound = new bool[numQueries+1];
  long numFound = mapQueries(queries, seekMatch, matchFound, numQueries);
  delete [] seekMatch;
  ScoredSeqWithCarriers** justFoundCopy = new ScoredSeqWithCarriers*[numFound+1];
  // copy the ones with hits to the new array, delete the rest
  long qN = 0;
  long fN = 0;
  while (fN < numFound){
    if (matchFound[qN]) { justFoundCopy[fN] = queries[qN]; ++fN; }
    else { queries[qN]->deepDelete(); }
    ++qN;
  }
  // delete the rest
  while (qN < numQueries){ queries[qN]->deepDelete(); ++qN; }
  // now copy back to the original array
  for (long n = 0; n < numFound; ++n){ queries[n] = justFoundCopy[n]; }
  delete [] matchFound;
  delete [] justFoundCopy;
  return numFound;
}
*/

long ExtendJobMapperNoLinks::mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
					      bool* seekReadMatch, bool* seekPairMatch,
					      bool* readMatchFound, bool* pairMatchFound, long numQueries){
  // simply forward the method - for now
  long hitCount = mapQueries(readQueries, seekReadMatch, readMatchFound, numQueries);
  hitCount += mapQueries(pairQueries, seekPairMatch, pairMatchFound, numQueries);
  return hitCount;
}
long ExtendJobMapperNoLinks::mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchFound, long numQueries){

  // save the ones with hits in this array; at the end, place them in the input array
  //ScoredSeqWithCarriers** queriesWithHits = new ScoredSeqWithCarriers*[ numQueries + 1 ];
  long countWithHits = 0;

  // make a non-redundant sequence set
  map<string, vector<long>*> seqToQuerySet;
  long numNr = 0;
  for (long n = 0; n < numQueries; ++n){
    // these false values will be changed to true if a match is found
    matchFound[n] = false;
    if (seekMatch[n]){
      char* querySeq = queries[n]->getSeq('+');
      map<string, vector<long>*>::iterator foundSeqIt = seqToQuerySet.find( querySeq );
      if (foundSeqIt == seqToQuerySet.end()){
	vector<long>* newCarrier = new vector<long>;
	newCarrier->push_back(n);
	seqToQuerySet.insert( pair<string, vector<long>*> (querySeq, newCarrier) );
      } else {
	foundSeqIt->second->push_back(n);
      }
      delete [] querySeq;
    } else {
      queries[n]->deepUnbuffer();
    }
  }

  // set up the arrays
  long numMonoP1 = seqToQuerySet.size() + 1;
  vector<long>** likeSequences = new vector<long>*[ numMonoP1 ];
  ScoredSeq** monoQueries = new ScoredSeq*[ numMonoP1 ];
  long* queryLengths = new long[ numMonoP1 ];
  // now use this as an index iterator - it will return to the full value
  long numMono = 0;
  for (map<string, vector<long>*>::iterator seqIt = seqToQuerySet.begin(); seqIt != seqToQuerySet.end(); ++seqIt){
    likeSequences[numMono] = seqIt->second;
    queryLengths[numMono] = seqIt->first.size();
    monoQueries[numMono] = new ScoredSeqMonoScore(seqIt->first, 1);
    ++numMono;
  }

  bool* qHasMatch = _mappingCollection->hasMatch(numMono, monoQueries, queryLengths,
						 ScoredSeqCollectionBwt::notThreaded,
						 ScoredSeqCollectionBwt::softMinOvl );
  for (long mqN = 0; mqN < numMono; ++mqN){
    if ( qHasMatch[mqN] ){
      for (vector<long>::iterator it = likeSequences[mqN]->begin(); it != likeSequences[mqN]->end(); ++it){
	matchFound[*it] = true;
	queries[*it]->deepUnbuffer();
      }
      countWithHits += likeSequences[mqN]->size();
    } else {
      for (vector<long>::iterator it = likeSequences[mqN]->begin(); it != likeSequences[mqN]->end(); ++it){
	queries[*it]->deepUnbuffer();
      }
    }
    delete monoQueries[mqN];
    delete likeSequences[mqN];
  }
  delete [] monoQueries;
  delete [] likeSequences;
  delete [] queryLengths;
  delete [] qHasMatch;

  return countWithHits;
}







#endif
