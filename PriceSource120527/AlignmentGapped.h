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

/* an alignment between two ScoredSeqs 

Conceptually, this is an immutable class.  However, realizing that
its construction is going to generally be a little bit complicated,
I moved the establishment of linkages between nucleotide positions
to a post-constructor method that adds linkages one at a time in a
manner that should be easy to understand for the user and therefore
difficult to make a careless mistake.

In order to prevent clients from inappropriately mutating alignments,
I included a "lock" function that disables all methods that mutate the
abstract state.  Those methods will raise an exception if they are
called after the alignment has been locked, and an "isLocked" method
is provided to determine if this is the case before calling the mutator
function.

I am also going to specify that linkages cannot be over-written, i.e.
once a position is linked in either sequence it cannot be re-linked.
This will prevent errors during creation, before "lock" is called.
Calling "lock" will prevent any future addition of linkages.  In no
situation can a linkage be removed.

*/

#ifndef ALIGNMENTGAPPED_H
#define ALIGNMENTGAPPED_H

#include <map>
#include "Alignment.h"
#include "ScoredSeq.h"
#include "ScoredSeqShallow.h"
#include "AlignmentScoreMatrix.h"
using namespace std;


class AlignmentGapped : public Alignment {

 public:
  ~AlignmentGapped();
  /* The orientation of the two input ScoredSeqs also indicates how
     the position integers should be interpreted.  It will be assumed
     that position integers refer to seqA in the (+) orientation and
     to seqB in the orientation specified in the constructor (also
     returned by the "orientation" method).
     REQUIRES: seqA != seqB
     REQUIRES: sense is '+' or '-'
  */
  AlignmentGapped(ScoredSeq * seqA, ScoredSeq* seqB, char orientation);

  /* allows the same alignment to be applied to two new seqeunces.
   * new sequences must be the same length as the old sequences.
   */
  AlignmentGapped* copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB);
  void seqReplace(ScoredSeq* seqA, ScoredSeq* seqB);
  AlignmentGapped* copyRcSeqA(ScoredSeq* newSeqA);

  ScoredSeq * seqA();
  ScoredSeq* seqB();
  char orientation();
  bool isNull();
  bool isLinked(long pos, ScoredSeq * seq);
  bool isGapped(long pos, ScoredSeq * seq);
  long getLinkage(long pos, ScoredSeq * seq);
  long gapPairedAfter(long pos, ScoredSeq * seq);
  long score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps);
  long scoreOverhangA(AlignmentScoreMatrix* scoreMatrix);
  bool isLocked();
  ScoredSeq* combine();
  ScoredSeq* combine(int denominator);
  bool hasConstantOffset();
  long getConstantOffset();

  void OK(); // checkRep

  // ############ MUTATOR METHODS ############ (should be specific to Valid class?)

  /* makes the Alignment immutable from the point lock is called forward */
  void lock();

  /* Creates a linkage between the indicated position on seqA (posSeqA) and
     seqB (posSeqB).  For posSeqA, numbering refers to the (+) strand position; 
     for posSeqB, numbering refers to the position on the strand indicated by
     the "orientation" function.  In the overloaded version, the strands of 
     the positions can be specified.
     THROWS: ??? if a linkage already exists for either seq/position
     THROWS: ??? if posSeq >= length of seq
     THROWS: ??? if this is locked
     REQUIRES: sense is '+' or '-'
  */
  void addLinkage(ScoredSeq * seq, long posSeq, long posPair);

  /* for each position on seqA starting at posSeqA and proceeding 5p->3p for
     blockLength nucleotides, a linkage is created for each position starting
     at posSeqB on seqB proceeding 5p-3p.  For posSeqA, numbering refers to the
     (+) strand position; for posSeqB, numbering refers to the position on the 
     strand indicated by the "orientation" function.  In the overloaded version,
     the strands of the positions can be specified.
     THROWS: ??? if a linkage already exists for any seq/position
     THROWS: ??? if posSeqX + blockLength > length of seqX
     THROWS: ??? if this is locked
     REQUIRES: sense is '+' or '-'
  */
  void addLinkageBlock(ScoredSeq * seq, long posSeq, long posPair, long blockLength);

  /* posGap is the position before the ins; if the gap is at the
   * beginning of the sequence, posGap = -1
   * NOTE: posGap is from the perspective of the sense in which
   * the sequence is stored in the alignment.
   */
  void addGap(ScoredSeq* inSeq, long posIns, long posGap);

  /* posGap is the position before the ins; if the gap is at the
   * beginning of the sequence, posGap = -1.  block will extend
   * in the 5p->3p direction of inSeq in its orientation in the alignment.
   */
  void addGapBlock(ScoredSeq* inSeq, long posIns, long posGap, long blockLength);


 private:
  AlignmentGapped(AlignmentGapped* oldAlignment, ScoredSeq* newSeqA, ScoredSeq* newSeqB);

  bool _locked;
  ScoredSeq * _seqA;
  ScoredSeq* _seqB;
  char _orientation;
  bool _isNull;  // only gains meaning when alignment is locked

  // -5 signals that the space is empty
  // size is the same as the key sequence
  long* _aToB; // A coord (key) is (+), B coord is (orientation)
  long* _bToA; // A coord is (+), B coord (key) is (orientation)

  // keys must be in range [0, seq length) - any gapped position
  // values must be in range [-1, seq length) - these are the position 
  //   of the nucleotide BEFORE the gap in the aligned sequence.
  // -5 signals that the space is empty
  long* _aGaps;
  long* _bGaps;

  // public fields; this is just to carry around sequences and their info
  // in order to minimize method calling
  class SeqInfoCarrier{
  public:
    SeqInfoCarrier(ScoredSeq* seq, char sense);
    ~SeqInfoCarrier();
    ScoredSeq* _seq;
    char* _nucs;
    float* _scores;
    float* _links;
    char _sense;
  };

  // these will be kept track of as the alignment is being filled in
  // so that the time for locking scales with the alignment length
  // rather than the sum of the sequence lengths.
  long _aMin;
  long _aMax;
  long _bMin;
  long _bMax;
  long _aMinLink;
  long _aMaxLink;

  // helpers for the score method
  long scoreGapBlock(long aStartPos, long aEndPos,
		     long bStartPos, long bEndPos,
		     AlignmentScoreMatrix* scoreMatrix,
		     bool penalizeTerminalGaps);
  // to use in place of a map for keeping track of/communicating gap positions/sizes;
  // when it is used, the assumption is that the positions were added in ascending order
  class GapToSize{
  public:
    GapToSize(long maxLen);
    ~GapToSize();
    void addPosition(long pos);
    void addCountToCurrent();
    long _gCount;
    long* _gStarts;
    long* _gLengths;
  private:
    long _gIndex;
  };
  GapToSize* getSegmentsFromBlockSeqA(long startPos, long endPos,
				      long maxOpposingPosition, bool penalizeTerminalGaps);
  GapToSize* getSegmentsFromBlockSeqB(long startPos, long endPos,
				      long maxOpposingPosition, bool penalizeTerminalGaps);

  // these are implemented as inline methods
  bool isGappedSeqA(long pos);
  bool isGappedSeqB(long pos);
  bool isLinkedSeqA(long pos);
  bool isLinkedSeqB(long pos);
  long getLinkageSeqA(long pos);
  long getLinkageSeqB(long pos);
  long gapPairedAfterSeqA(long pos);
  long gapPairedAfterSeqB(long pos);

  long reversePosCoord(long pos, long seqLenM1);
  long reversePosCoordGap(long pos, long seqLenM1);

  void addGapSeqA(long posIns, long posGap);
  void addGapSeqB(long posIns, long posGap);
  void addGapNoUpdateSeqA(long posIns, long posGap);
  void addGapNoUpdateSeqB(long posIns, long posGap);

  void helpAddLinkSeqA(long posSeq, long posPair);
  void helpAddLinkSeqB(long posSeq, long posPair);
  void helpAddGapBlockSeqA(long posIns, long posGap, long blockLength);
  void helpAddGapBlockSeqB(long posIns, long posGap, long blockLength);

  void updateMinMaxA(long pos);
  void updateMinMaxB(long pos);

  // these are for optimization of the inline methods
  long _sizeA;
  long _sizeB;
  long _sizeAm1;
  long _sizeBm1;
  long _sizeAp5;
  long _sizeBp5;


  // this is just a vehicle for carrying around data; its data fields are public
  class CombineDataCarrier{
  public:
    CombineDataCarrier();
    CombineDataCarrier(long maxSize);
    ~CombineDataCarrier();
    void addNuc(char nuc, float score);
    void addNuc(char nuc, float score, float link);
    void addLink(float link);
    void decrementLink();
    ScoredSeq* getSeq();
  private:
    long _size;
    char* _seq;
    float* _scores;
    float* _links;
    long _linkOffset;
    long _maxSize;
  };


  static bool isPositionFinal(ScoredSeq * seq, long pos);
  static bool isPositionFinal(ScoredSeq * seqX, long posX, ScoredSeq * seqY, long posY);

  // I will eventually replace the method above with the one below
  ScoredSeq* combine(ScoreCalculator* calc);


  void combineGapBlock(SeqInfoCarrier* infoA, SeqInfoCarrier* infoB,
		       long aStartPos, long aEndPos,
		       long bStartPos, long bEndPos,
		       CombineDataCarrier* carrier, ScoreCalculator* outsideCalc);


  void getSegmentsFromBlock(map<long,float>* links, long startPos, long endPos, char sense,
			    ScoredSeq* seq, ScoreCalcXy* calc);
  static void fillInGaps(CombineDataCarrier* carrier,
			 SeqInfoCarrier* seqInfo, long startPos, long endPos,
			 float penaltyConstant, map<long,float>* links, ScoreCalcXy* calc);


};


#endif

