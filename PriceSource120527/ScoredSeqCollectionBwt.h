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

/* A collection of ScoredSeq objects indexed for quick retrieval

This version is optimized for high speed when high fractId is required.
The set of ScoredSeqs that it contains cannot be modified.

The contents of this collection are typed.
 */


#ifndef SCOREDSEQCOLLECTIONBWT_H
#define SCOREDSEQCOLLECTIONBWT_H

#include <omp.h>
#include <set>
#include <map>
#include <vector>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "DynamicProgrammingAligner.h"
#include "AlignmentScoreMatrix.h"
#include "burrowsWheelerUtilities.h"
#include "MatchSeqTest.h"
#include "MatchOffsetTest.h"
using namespace::std;


class ScoredSeqCollectionBwt{

 public:
  // for constructors, note that there would be no advantage to passing in an array

  // these two will be used if ungapped alignments are to be returned
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, bool penalizeEdgeGaps=true);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap, long maxOverlap);
  // these two will be used if gapped alignments are to be returned
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DynamicProgrammingAligner * asi);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DynamicProgrammingAligner * asi, long globalMaxOverlap);
  ~ScoredSeqCollectionBwt();

 
  enum GapMode { gapped=0, ungapped=1 };
  enum MinOvlStringency { hardMinOvl=0, softMinOvl=1 };
  enum AlignmentType { semiGlobal=0, fullAlignment=1 };
  enum AlignmentThreadedness{ notThreaded=0, threaded=1 };
  enum MatchesSought{ allMatches=0, bestMatches=1, firstMatch=2 };

  // DEFAULT: edge scaling is enabled
  // that means that matches to smaller words from the edges of query contigs
  // will be sought in order to find matches with small overlaps.
  void enableEdgeScaling();
  void disableEdgeScaling();

  // returns a deep copy of the collection that uses only the
  // shallow versions of the contained objects so that there are
  // no pointers that could be used in some other part of the
  // software (excepting singletons).
  ScoredSeqCollectionBwt * copy();

  // methods for dealing with contents
  bool contains(ScoredSeq* s);
  void getSeqs( set<ScoredSeq*>* );

  // learn about the collection
  float getFractId();
  long getMinOverlap();

  // In the methods below, a pre-created set of alignments is passed in as a reference.
  // It will be modified by having the new matches added to it.

  // these methods return just the best matches
  // (according to score; all equal matches are returned)
  //void getBestFullMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense = '.', MatchSeqTest* seqTest = NULL);

  // these allow a new minOverlap to be set
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'

  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* seqTest,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* seqTest,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'

  // THREADED
  /*
  void getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
			  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  void getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
			  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  */
  // these allow a new minOverlap to be set
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** seqTestArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);


  // these just determine if there is a full match (either full query or full target)
  bool hasFullMatch(ScoredSeq* seq, char sense = '.', MatchSeqTest* seqTest = NULL);
  // minOverlap MUST be set
  bool hasMatch(ScoredSeq* seq, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, char sense, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  // THREADED

  bool* hasFullMatch(long numQueries, ScoredSeq** seqArray,
		     AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  bool* hasFullMatch(long numQueries, ScoredSeq** seqArray, char sense,
		     AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);

  bool* hasMatch(long numQueries, ScoredSeq** seqArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, char sense, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, MatchSeqTest** seqTestArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);




  // GETMATCHES

  // and i have equivalents that return all matches
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* matchTest,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* matchTest,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  // these will be threaded
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** matchTestArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** matchTestArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'


  // these two methods deal with the buffering for all the member seqs.  they 
  // don't change the abstract state; they are provided as tools for saving memory
  // without destroying the collection; their converse operations are private and 
  // will be called based on demand (i.e. when matches are requested).
  void unbufferSeqs();
  void unindexSeqs();

  long size(); // number of member ScoredSeqs
  void OK();

  void bufferSeqs();

  // this is only public so that I can test it
  //Alignment* reverseAlignmentSeqA(Alignment* oldAl, ScoredSeq* newSeqA);



 private:
  void constructorHelper(set<ScoredSeq*>* inputSeqs);

  bool _penalizeEdgeGaps;
  static long _defaultMinScore;

  // this should be true during assembly and false during read mapping
  bool _edgeScalingEnabled;

  // this is used in several places
  static long getOverlapFromOffset(long offset, long seqSizeA, long seqSizeB);

  // private versions of these methods; the above methods will forward to them
  void runSearchesConverter(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
			    MatchSeqTest* seqTest, long minOverlap,
			    AlignmentType alType, MatchesSought matchesSought,
			    MinOvlStringency ovlStringency, AlignmentThreadedness threadedness);
  void runSearchesPrivate(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
			  MatchSeqTest** seqTestArray, long* minOvlArray,
			  AlignmentType alType, MatchesSought matchesSought,
			  MinOvlStringency ovlStringency, AlignmentThreadedness threadedness);

  // these methods create an internal representation with a bigger memory
  // footprint and so are called based on demand (i.e. when matches are 
  // requested).
  void indexSeqs();
  void indexSeq(ScoredSeq * seq); // some functional abstraction

  float _fractId;
  float _maxFractMis;
  long _maxOverlap;
  long _minOverlap;
  bool _usingMaxOverlap;

  GapMode _gapMode;
  // a possible alternative to _minOverlap for full alignments; will be
  // set to the length of the shortest sequence or _minOverlap, whichever
  // is longer
  long _minFragmentSize;

  // these are for the gapped alignment mode
  DynamicProgrammingAligner * _asi;


  // BWT stuff
  long* _sortedSuffixes; // also an array of text length
  long* _bwTransform; // also an array of text length
  long* _occCounts; // also an array of text length
  long* _tableC; // an array of alphabet length
  long _textSize;

  // I NEED A CLASS OF SCOREDSEQ WRAPPERS THAT DOES THE OVERLAP FUNCTION AND KEEPS TRACK
  // OF COORD POSITION IN THE META-TEXT USED FOR BWT MAPPING.
  // used to store the ScoredSeqs so they can be queried
  class ScoredSeqLocus {
  public:
    ScoredSeqLocus();
    ScoredSeqLocus(ScoredSeq* seq, long start, long end, long index);
    bool overlaps(long otherStart, long otherEnd);
    bool contains(long otherStart, long otherEnd);
    ScoredSeq* _seq;
    long _start;
    long _end;
    long _index;
  private:
  };



  // helps enforce the limits of offsets that may be generated from
  // a given query's word match
  class WordQuery{
  public:
    virtual ~WordQuery();
    virtual long size() = 0;
    virtual long start() = 0;
    // end is pos+1, good for iterating
    virtual long end() = 0;
    virtual bool acceptableOverlap(long targetStart, ScoredSeq* target) = 0;
    virtual bool notContinuous(WordQuery* otherWq) = 0;
  };


  class WordQuerySimple : public WordQuery {
  public:
    WordQuerySimple(long start, long size);
    WordQuerySimple(WordQuery* wq);
    ~WordQuerySimple();
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    // class-specific; REQUIRES continuity; modifies this
    void merge(WordQuery* wq);
  private:
    // define the souce seq and position on it of the query word
    long _start;
    long _end;
    long _size;
  };

  class WordQueryNoEdgeLimit : public WordQuery {
  public:
    WordQueryNoEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense);
    ~WordQueryNoEdgeLimit();
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
  private:
    // define the souce seq and position on it of the query word
    ScoredSeq* _sourceSeq;
    long _start;
    long _end;
    long _size;
    char _sense;
  };

  class WordQueryFiveEdgeLimit : public WordQuery {
  public:
    WordQueryFiveEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense, long edgeLimit);
    ~WordQueryFiveEdgeLimit();
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
  private:
    // define the souce seq and position on it of the query word
    ScoredSeq* _sourceSeq;
    long _start;
    long _end;
    long _size;
    char _sense;
    // define the distance to the edge of a contig within
    // which the word match must be found (or its relevance)
    long _edgeLimit;
  };
  class WordQueryThreeEdgeLimit : public WordQuery {
  public:
    WordQueryThreeEdgeLimit(ScoredSeq* sourceSeq, long start, long size, char sense, long edgeLimit);
    ~WordQueryThreeEdgeLimit();
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
  private:
    // define the souce seq and position on it of the query word
    ScoredSeq* _sourceSeq;
    long _start;
    long _end;
    long _size;
    char _sense;
    // define the distance to the edge of a contig within
    // which the word match must be found (or its relevance)
    long _edgeLimit;
  };

  class WordQueryFullTarget : public WordQuery {
  public:
    WordQueryFullTarget(ScoredSeq* sourceSeq, long start, long size, char sense, float fractId, long maxTargetSize);
    ~WordQueryFullTarget();
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
  private:
    // define the souce seq and position on it of the query word
    ScoredSeq* _sourceSeq;
    long _start;
    long _end;
    long _size;
    char _sense;
    float _fractId;
    long _maxTargetSize;
  };

  // used to deal with redundant word queries
  class WordQueryBinOrganizer {
  public:
    WordQueryBinOrganizer(vector<WordQuery*>* wordQueries);
    ~WordQueryBinOrganizer();

    long numQueries();
    bool queryStillValid(long qIndex);
    WordQuery* getQuery(long qIndex);

    void removeQuery(long qIndex);
    void removeQuery(long qIndex, long currentBlock);
    void removeSupersetQueries(long qIndex);

    void getQueries(vector<WordQuery*>* wqSet);
    long numBlocks();
    bool blockHasQueries(long block);

    long blockStart(long block);
    long blockEnd(long block);
    bool blockHasQuery(long block, long qIndex);

    //void getBlockQueryIndexes(long block, set<long>* indexes);
    // these two methods replace the one above
    long getBlockQueryIndexCount(long block);
    long* getBlockQueryIndexes(long block);

  private:
    WordQuery** _queries;
    bool* _queryIsValid;
    long _numQueries;
    long _numBlocks;
    // for these, the N value in End equals the N+1 value in start; pre-computing
    // is faster than adding one every time the boundaries are checked
    long* _blockToStart;
    long* _blockToEnd;
    long* _startToBlock;
    // organizes the data by block
    long* _blockToBinCount;
    // has the qIndexes of the queries in that bin in an arbitrary order;
    // # useful elements == _blockToBinCOunt value
    long** _blockToQisInBin;
    // sparse: for all of the qIndexes in the bin, there is a value showing
    // where they appear in the _blockToQisInBin array
    long** _blockToQiBinIndex;
  };


  class OffsetAndQueryCarrier {
  public:
    OffsetAndQueryCarrier(long offset, long numQueries);
    ~OffsetAndQueryCarrier();
    long _offset;
    long _numQueries;
    WordQuery** _queries;
  };

  // for keeping track of offsets and passing them around
  class OffsetTracker {
  public:
    OffsetTracker(ScoredSeq* target);
    ~OffsetTracker();
    void addOffset(long offset, WordQuery* seedingQuery, long normFactor=1);
    long maxOffsetCount();
    long maxCountOffset();
    long numOffsets();
    //long* getOffsets();
    OffsetAndQueryCarrier** getOffsets(bool includeQueries);
    ScoredSeq* targetContig();
  private:
    class OffsetFields {
    public:
      OffsetFields(float count, WordQuery* maxWord);
      ~OffsetFields();
      void addOffset(float localAddition, WordQuery* seedingQuery);
      float _count;
      WordQuerySimple* _maxWord;
      vector<WordQuerySimple*> _otherWords;
    };
    map<long,OffsetFields*> _offsetToFields;
    long _maxCountOffset;
    float _maxOffsetCount;
    ScoredSeq* _targetContig;
    // for sorting word queries by start coord
    struct SortWqByStart {
      bool operator() (WordQuery* wqA, WordQuery* wqB);
    };

  };




  // pure abstract class that will replace the one above; I will change
  // the name when I am done with it
  class AlignmentMaker {
  public:
    virtual ~AlignmentMaker();
    virtual Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier) = 0;
    virtual long getScore(Alignment* al) = 0;
    virtual long getMinScore() = 0;
    virtual void setMinScore(long minScore) = 0;
    virtual AlignmentScoreMatrix* getScoreMatrix() = 0;
    //static bool makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis, long alLength);
    static bool makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis,
					long exBlockCount, long* examineStarts, long* examineEnds);
    //static bool makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis, long alLength);
    static bool makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis,
					 long exBlockCount, long* examineStarts, long* examineEnds);
  };
  // implementing classes
  class AlMakerNoEdge : public AlignmentMaker {
  public:
    AlMakerNoEdge();
    AlMakerNoEdge(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis);
    ~AlMakerNoEdge();
    Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier);
    long getScore(Alignment* al);
    long getMinScore();
    void setMinScore(long minScore);
    AlignmentScoreMatrix* getScoreMatrix();
  private:
    AlignmentScoreMatrix* _scoreMatrix;
    long _minScore;
    float _maxFractMis;
  };
  class AlMakerPenalizeEdgeA : public AlignmentMaker {
  public:
    AlMakerPenalizeEdgeA();
    AlMakerPenalizeEdgeA(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis);
    ~AlMakerPenalizeEdgeA();
    Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier);
    long getScore(Alignment* al);
    long getMinScore();
    void setMinScore(long minScore);
    AlignmentScoreMatrix* getScoreMatrix();
  private:
    AlignmentScoreMatrix* _scoreMatrix;
    long _minScore;
    float _maxFractMis;
  };



  // I will use this to find the ScoredSeqs to which queries were mapped quickly
  // indexes are 0 -> total text length / _binSize (continuous)
  // this will be an array of ScoredSeqLocus* sets.
  vector<ScoredSeqLocus*>* _binToSeqs;

  // an array of loci
  ScoredSeqLocus** _orderedContigLoci;
  long _numContigLoci;

  long _binSize;
  long _maxBin; // for safety checks
  long getOverlapContigIndex(long start, long end);
  AlignmentScoreMatrix* _matchCountAsm;
  AlignmentScoreMatrix* _misCountAsm; // gaps count as mismatches here (terminal gaps will exist)
  AlignmentScoreMatrix* _scoreAsm;

  // helpers for functions above that use private datatypes
  void runAlignmentsHelper(long numQueries, OffsetTracker*** contigToOffsetsArray,
			   vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, GapMode gapMode,
			   MatchSeqTest** seqTestArray, long* minOvlArray, long* bestScoreArray,
			   MatchesSought matchesSought, AlignmentThreadedness threadedness);
  AlignmentMaker* makeAlMaker(GapMode gapMode, long bestScoreSoFar);
  void dealWithMatches(vector<Alignment*>* tempMatches, vector<Alignment*>* realMatches,
		       bool getAllMatches, bool justFirstMatch,
		       long index, long* bestScoreArray, long newBestScore, bool* matchSoughtArray);


  // if no offsets need to be generated
  void setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
			       AlignmentMaker** alMakerArray, ScoredSeq** seqFlipArray);
  void setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
			       ScoredSeq** seqFlipArray);
  // if offsets need to be generated
  long setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
			       long* minOvlArray, MatchSeqTest** matchTestArray,
			       AlignmentMaker** alMakerArray, ScoredSeq** seqFlipArray,
			       OffsetTracker*** conToOffByQuery, long* otCountByQuery, AlignmentType alType);
  long setupPrivateMultiHelper(long index, ScoredSeq** seqArray, char sense,
			       long* minOvlArray, MatchSeqTest** matchTestArray,
			       ScoredSeq** seqFlipArray,
			       OffsetTracker*** conToOffByQuery, long* otCountByQuery, AlignmentType alType);

  // this is used in many contexts; its calling is controlled by the hasMatch, getBestMatches, and getMatches
  // wrapper methods.  note that the default min score so far is zero, but i am going to force that to be
  // input so that i don't fuck up and forget to include it somewhere
  // NOTE: if onlyBestAlignment==false, then the input bestScoreSoFar will just be spit back without being updated
  long makeAlignmentHelper(OffsetTracker* tracker, ScoredSeq* seq, vector<Alignment*>* matches, long minOverlap,
			   AlignmentMaker* alMaker, ScoredSeq* seqFlip, char sense, long bestScoreSoFar,
			   bool onlyBestAlignment, GapMode gapMode, MatchesSought matchesSought);


  // used to define the set of subtext strings for which perfect matches will be sought
  static int _specialCaseLenDenom[6];
  static int _specialCaseStepDenom[6];
  static int _specialCaseWindowNum[6];
  void getSubtextWordQueries(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet, char sense, long minOverlap);
  void getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
				   long blockStartPos, long blockEndPos, bool useLimits);
  void subtextWordQueriesFilterNs(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet);

  void getSubtextWordQueriesFullAlignment(ScoredSeq* querySeq, char* seqSeq, vector<WordQuery*>* wordQuerySet,
					  char sense, long minTargetLength);
  void getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
				   long fragmentLength);

  long getNumWindowsHelper(long numAllowedMis);

  // a helper for getMatches/getBestMatches
  OffsetTracker** getOffsets(ScoredSeq* seq, char sense, long minOverlap, MatchSeqTest* matchTest, AlignmentType alType);

  // for sorting word queries by length, shortest to longest
  struct SortWqByLength {
    bool operator() (WordQuery* wqA, WordQuery* wqB);
  };

  // these keep track of the total content, binned by whether or not the
  // work that needs to be done to get the collection ready to use has
  // already been performed.
  set<ScoredSeq*> _allSeqs;

};

#endif
