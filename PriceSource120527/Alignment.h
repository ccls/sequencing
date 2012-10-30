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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <map>
#include "ScoredSeq.h"
#include "ScoredSeqShallow.h"
#include "AlignmentScoreMatrix.h"
using namespace std;


class Alignment {
 public:
  virtual ~Alignment();

  virtual Alignment* copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB) = 0;
  virtual void seqReplace(ScoredSeq* seqA, ScoredSeq* seqB) = 0;
  // the replacement for seqA (newSeqA) is its revers complement, so the alignment
  // orientation needs to be reversed
  virtual Alignment* copyRcSeqA(ScoredSeq* newSeqA) = 0;

  /* Get a reference to the first ScoredSeq (the one whose coords are
     always reported for the (+) strand). 
  */
  virtual ScoredSeq * seqA() = 0;

  /* Get a reference to the second ScoredSeq (the one whose coords are
     reported for the strand indicated by the "orientation" method).
  */
  virtual ScoredSeq* seqB() = 0;

  /* The orientation of seqB in the alignment; coords reported for seqB
     are numbered with respect to this orientation (see ScoredSeq specs 
     for more description).
     RETURNS: '+' or '-' char
  */
  virtual char orientation() = 0;

  /* If there are no aligned nucleotides in the alignment, then the
   * alignment is "null" (i.e. empty)
   * REQUIRES: alignment isLocked
   */
  virtual bool isNull() = 0;

  /* Indicates whether or not the nucleotide at the indicated position
     of the indicated ScoredSeq is aligned.  The overloaded version
     permits specification of the strand on seq to which the provided
     pos refers.
     THROWS: ??? if pos >= length of seq
     REQUIRES: sense is '+' or '-'
  */
  virtual bool isLinked(long pos, ScoredSeq * seq) = 0;

  /* like isLinked, but tells whether or not a nucleotide is part
     of a gap
  */
  virtual bool isGapped(long pos, ScoredSeq * seq) = 0;

  /* For the nucleotide at the indicated position of the indicated ScoredSeq,
     returns the position in the other sequence to which it is aligned.  The
     sense of the returned position is determined by the convention of the 
     class: for seqA, it is the (+) strand coord; for seqB, it is the coord
     for the strand indicated by the "orientation" function.  The overloaded
     version permits specification of the strand on seq to which the provided
     posSeq refers (senseProvided) and for the returned position's seq (senseOther).
     THROWS: ??? if pos of seq is unlinked (isLinked can provide that info
             ahead of time.
     THROWS: ??? if pos >= length of seq
     REQUIRES: sense is '+' or '-'
  */
  virtual long getLinkage(long pos, ScoredSeq * seq) = 0;

  /* Like getLinkage, but gets the position in the other sequence that the
   * pos in seq is gapped AFTER
   * if the pos in seq is paired with gaps before the entire other sequence,
   * then -1 is returned.
   * REQUIRES: pos in seq isGapped
   */
  virtual long gapPairedAfter(long pos, ScoredSeq * seq) = 0;

  /* calculates a score using the given scoring matrix;
   * note that this alignment may not be the optimal alignment
   * for the given scoring matrix. 
   * REQUIRES: alignment is locked.*/
  virtual long score(AlignmentScoreMatrix* scoreMatrix, bool penalizeTerminalGaps) = 0;
  virtual long scoreOverhangA(AlignmentScoreMatrix* scoreMatrix) = 0;

  /* determines if alignment can be mutated. */
  virtual bool isLocked() = 0;

  virtual void OK() = 0; // checkRep

  // I am going to replace the static methods above with these, eventually
  virtual ScoredSeq* combine() = 0;
  // seqA will be divided by denominator for normalization
  virtual ScoredSeq* combine(int denominator) = 0;

  // used for multiple simultaneous collapses, which are only possible for
  // ungapped alignments (for which this value is true)
  // REQUIRES: alignment is not null
  virtual bool hasConstantOffset() = 0;
  // offset is the position on seqA where seqB begins (negative if seqB overhangs in front)
  // REQUIRES: above function is true
  virtual long getConstantOffset() = 0;

  // requires that the alignment is ungapped and that A is a subset of B; the bool method
  // will say if this is legit or not
  /*
  virtual bool canCollapseIntoB() = 0;
  virtual void collapseIntoB() = 0;
  virtual void collapseIntoB(int denominator) = 0;
  */

 protected:

  // THIS STUFF IS HERE SO THAT I DON'T NEED TO RE-WRITE THESE CLASSES FOR EVERY IMPLEMENTING CLASS

  static float _SCORE_MAX;

  // this is a way of providing altenative methods for calculating scores
  class ScoreCalculator {
  public:
    virtual ~ScoreCalculator();
    virtual float adjustA(float a) = 0;
    virtual float adjustB(float b) = 0;
    virtual float aPlusB(float a, float b) = 0;
    virtual float aMinusB(float a, float b) = 0;
    virtual float bMinusA(float b, float a) = 0;
  };
  class ScoreCalcAdd : public ScoreCalculator {
  public:
    ScoreCalcAdd();
    ~ScoreCalcAdd();
    float adjustA(float a);
    float adjustB(float b);
    float aPlusB(float a, float b);
    float aMinusB(float a, float b);
    float bMinusA(float b, float a);
  };
  class ScoreCalcNormalizeA : public ScoreCalculator {
  public:
    ScoreCalcNormalizeA();
    ScoreCalcNormalizeA(int denominator); // divides seqA by denominator
    ~ScoreCalcNormalizeA();
    float adjustA(float a);
    float adjustB(float b);
    float aPlusB(float a, float b);
    float aMinusB(float a, float b);
    float bMinusA(float b, float a);
  private:
    int _denominator;
  };

  class ScoreCalcXy { // allows easy swapping of A vs B identity without ambiguous labelling
  public:
    static ScoreCalcXy* getScoreCalcXy(ScoreCalculator* calc, bool xIsA);
    virtual ~ScoreCalcXy();
    virtual float adjustX(float x) = 0;
    virtual float adjustY(float y) = 0;
    virtual float xPlusY(float x, float y) = 0;
    virtual float xMinusY(float x, float y) = 0;
    virtual float yMinusX(float y, float x) = 0;
  };
  class ScoreCalcXyA : public ScoreCalcXy {
  public:
    ScoreCalcXyA(ScoreCalculator* calc);
    ~ScoreCalcXyA();
    float adjustX(float x);
    float adjustY(float y);
    float xPlusY(float x, float y);
    float xMinusY(float x, float y);
    float yMinusX(float y, float x);
  private:
    ScoreCalculator* _calc;
  };
  class ScoreCalcXyB : public ScoreCalcXy {
  public:
    ScoreCalcXyB(ScoreCalculator* calc);
    ~ScoreCalcXyB();
    float adjustX(float x);
    float adjustY(float y);
    float xPlusY(float x, float y);
    float xMinusY(float x, float y);
    float yMinusX(float y, float x);
  private:
    ScoreCalculator* _calc;
  };



};



#endif

