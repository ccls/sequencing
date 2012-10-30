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

/*
 * This is a wrapper for a ScoredSeq that allows its scores to
 * be automatically normalized (divided by the normalization count).
 * It allows that count to be added to in order to facilitate easy
 * and fast counting/normalizing.  The initial count is zero, and the
 * returned link scores are not normalized in that case (i.e. the
 * same result as when the count is 1).
 */

#ifndef SCOREDSEQNORMALIZED_H
#define SCOREDSEQNORMALIZED_H

#include "ScoredSeq.h"
#include "ScoredSeqNested.h"
#include "ScoredSeqPaired.h"
#include "AssemblyException.h"
#include "AssemblyException.h"


class ScoredSeqNormalized : public ScoredSeqNested, public ScoredSeqPaired {
 public:
  ScoredSeqNormalized(); //default constructor
  ScoredSeqNormalized(ScoredSeq* innerSeq);
  ~ScoredSeqNormalized();

  static long _readCount;

  /* see parent abstract class ScoredSeq */
  float scoreAtPosition(long position, char sense);
  float scoreAtPosPlus(long position);
  float scoreAtPosMinus(long position);
  float linkAfterPosition(long position, char sense);
  float linkAfterPosPlus(long position);
  float linkAfterPosMinus(long position);
  char nucAtPosition(long position, char sense);
  char nucAtPosPlus(long position);
  char nucAtPosMinus(long position);
  char* getSeq(char sense);
  float* getScores(char sense);
  float* getLinks(char sense);
  char* getSubseq(long position, long length, char sense);
  float* getSubScores(long position, long length, char sense);
  float* getSubLinks(long position, long length, char sense);
  void buffer();
  void unbuffer();
  void bottomBuffer();
  long size();
  bool hasPairedEnd();
  ScoredSeq * getPairedEnd();
  ScoredSeq * getTempPairedEnd();
  ScoredSeq * shallowCopy();
  ScoredSeq * flipCopy();
  bool isNested();
  /* see parent abstract class ScoredSeqNested */
  ScoredSeq* getNested();
  void deepDelete();
  void deepUnbuffer();

  /* see parent abstract class ScoredSeqPaired */
  void makeAlive();
  bool isAlive();

  // class-specific
  void addCount();
  void addCounts(long counts);
  long getCount();
  void resetCount();

  void OK();

 protected:
  friend class ReadFileCommunicator;
  void addPairedEnd(ScoredSeqNormalized * ss);

 private:

  /* informs this that its pair is dead; authenticity is verified by
   * passing the dying Read (oldPair) and its destructable replacement. 
   * note that passing out the paired end again should restore its
   * indestructable status */
  void makeDead();
  void replaceDeadPe(ScoredSeqNormalized * deadReplacement);

  void shallowDelete(); // does not delete a dead PE
  long _normCount;
  ScoredSeq* _innerSeq;
  ScoredSeq* _bufferSeq;
  ScoredSeqNormalized* _pe;
  bool _thisIsAlive;
  bool _okToDeletePe;
};

#endif
