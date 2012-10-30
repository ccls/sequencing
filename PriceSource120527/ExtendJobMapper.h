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

/* This class separates out the task of mapping reads to existing
 * contigs from the rest of the business of running an ExtendCycle.
 */
#ifndef EXTENDJOBMAPPER_H
#define EXTENDJOBMAPPER_H

#include <set>
#include <map>
#include "ScoredSeq.h"
#include "ScoredSeqWithCarriers.h"
#include "ScoredSeqFlip.h"
#include "ScoredSeqCollectionBwt.h"
#include "ParamsMapping.h"
#include "ExtendCycle.h"
using namespace::std;


class ExtendJobMapper {
 public:
  virtual ~ExtendJobMapper();
  //virtual long mapQueries(ScoredSeqWithCarriers** queries, long numQueries) = 0;
  // the returned long is the number of queries for which the matchFound value is "true"
  virtual long mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchFound, long numQueries) = 0;
  // the returned long is the total number of queries for which the matchFound value is "true"
  virtual long mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
				bool* seekReadMatch, bool* seekPairMatch,
				bool* readMatchFound, bool* pairMatchFound, long numQueries) = 0;

  class QueryCarrier {
  public:
    QueryCarrier(char* query, long size);
    ~QueryCarrier();
    char* _query;
    long _size;
    vector<ScoredSeqWithCarriers*> _seqs;
  };
};


class ExtendJobMapperWithLinks : public ExtendJobMapper {
 public:
  //ExtendJobMapperWithLinks();
  // no max value
  ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId);
  // edgeSize < 0 for no max value
  ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize);

  ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId);
  ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize);

  ~ExtendJobMapperWithLinks();

  // deep deletes any queries without matches; modifies the query array by
  // placing those that mapped to something at the front of the array and
  // deleting everything with no hits; returns the number of returned queries
  // with hits
  //long mapQueries(ScoredSeqWithCarriers** queries, long numQueries);
  // does NOT delete reads without hits
  long mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchFound, long numQueries);
  long mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
			bool* seekReadMatch, bool* seekPairMatch,
			bool* readMatchFound, bool* pairMatchFound, long numQueries);

 private:
  void constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers);

  long _edgeSize;
  int _plusIndex;
  int _minusIndex;
  int _frontIndex;
  int _backIndex;

  // for read mapping
  float _fractMapId;

  // references are Subseq encasing WC encasing Normalized encasing Shallow
  ScoredSeqCollectionBwt* _mappingCollection;
  ScoredSeq** _shortArray;
  ScoredSeq** _flipArray;
  long _disposableIndex;
  // RETURNS array of unique hits to answer about passing the test
  // MODIFIES qToUnique, matchCarrier
  // for matchCarrier, arrays will be created; 0=>plus, 1=>minus
  bool* mapQueriesSetup(ScoredSeqWithCarriers** queries, bool* seekMatch,
			long numQueries, long* qToUnique,
			vector<ScoredSeqWithCarriers*>** matchCarrier);
  // the returned long is the number of queries for which the matchFound value is "true"
  long mapQueriesHelper(ScoredSeqWithCarriers** queries, long* qToUnique, bool* seekMatch,
			vector<Alignment*>** matchesArray, ScoredSeq** monoQueries, long numQueries);
  //bool passesTest(vector<Alignment*>* matches);
  void sortMatchesHelper(vector<Alignment*>* matches,
			 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
			 vector<ScoredSeqWithCarriers*>* matchingContigsMinus);
  void findMatchesHelper(ScoredSeq* currentQuery,
			 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
			 vector<ScoredSeqWithCarriers*>* matchingContigsMinus);
  // here, "side" should be 0 for plus and 1 for minus
  void addToContigMatchesHelper(ScoredSeqWithCarriers* seqToAdd, int side, int frontOrBack, vector<ScoredSeqWithCarriers*>* matchingContigs);

  // this method gets rid of links to both strands of the hit contig
  // MODIFIES both input vectors
  void removeReciprocalMatches(vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
			       vector<ScoredSeqWithCarriers*>* matchingContigsMinus);

  // used to ensure that at least one read still has a best match to the end of a contig
  class AlignmentTest {
  public:
    virtual ~AlignmentTest();
    // REQUIRES: query is seqA in alignment
    virtual bool passes(vector<Alignment*>* matches, ScoredSeq* query) = 0;
  };
  class AlTestNull : public AlignmentTest {
  public:
    AlTestNull();
    ~AlTestNull();
    bool passes(vector<Alignment*>* matches, ScoredSeq* query);
  };
  class AlTestEdgeProximity : public AlignmentTest {
  public:
    AlTestEdgeProximity(long edgeProximity, float fractId);
    ~AlTestEdgeProximity();
    bool passes(vector<Alignment*>* matches, ScoredSeq* query);
  private:
    long _edgeProximity;
    float _fractId;
  };
  AlignmentTest* _alTest;
};



class ExtendJobMapperNoLinks : public ExtendJobMapper {
 public:
  // no max value
  ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId);
  // edgeSize < 0 for no max value
  ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize);
  ~ExtendJobMapperNoLinks();
  // deep deletes any queries without matches; modifies the query array by
  // placing those that mapped to something at the front of the array and
  // deleting everything with no hits; returns the number of returned queries
  // with hits
  //long mapQueries(ScoredSeqWithCarriers** queries, long numQueries);
  long mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchFound, long numQueries);
  long mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
			bool* seekReadMatch, bool* seekPairMatch,
			bool* readMatchFound, bool* pairMatchFound, long numQueries);

 private:
  void constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers);

  long _edgeSize;
  // for read mapping
  float _fractMapId;
  // references are Subseq encasing WC encasing Normalized encasing Shallow
  ScoredSeqCollectionBwt* _mappingCollection;
};

#endif
