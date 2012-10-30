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

#ifndef ALIGNMENTNULL_H
#define ALIGNMENTNULL_H

#include <map>
#include "Alignment.h"
#include "ScoredSeq.h"
#include "ScoredSeqShallow.h"
#include "AlignmentScoreMatrix.h"
using namespace std;

class AlignmentNull : public Alignment {

 public:
  ~AlignmentNull();
  /* The orientation of the two input ScoredSeqs also indicates how
     the position integers should be interpreted.  It will be assumed
     that position integers refer to seqA in the (+) orientation and
     to seqB in the orientation specified in the constructor (also
     returned by the "orientation" method).
     REQUIRES: seqA != seqB
     REQUIRES: sense is '+' or '-'
  */
  AlignmentNull(ScoredSeq * seqA, ScoredSeq* seqB, char orientation);

  AlignmentNull* copySeqReplace(ScoredSeq* seqA, ScoredSeq* seqB);
  void seqReplace(ScoredSeq* seqA, ScoredSeq* seqB);
  AlignmentNull* copyRcSeqA(ScoredSeq* newSeqA);

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

 private:
  ScoredSeq* _seqA;
  ScoredSeq* _seqB;
  char _orientation;

};


#endif

