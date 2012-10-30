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
 * This is a wrapper for a ScoredSeq that allows it to carry other
 * ScoredSeqs with it (like, say, matches from a mapping exercise).
 * It avoids the use of a lookup table in that case.
 * It allows as many distinct sets as the user wants, each with a
 * unique numerical index.
 */

#ifndef SCOREDSEQWITHCARRIERS_H
#define SCOREDSEQWITHCARRIERS_H

#include "ScoredSeqNested.h"
#include "ScoredSeqPaired.h"
#include "AssemblyException.h"
#include <set>

using namespace::std;

class ScoredSeqWithCarriers : public ScoredSeqNested, public ScoredSeqPaired {
 public:
  ScoredSeqWithCarriers(); //default constructor
  ScoredSeqWithCarriers(ScoredSeq* innerSeq, int numOfSets);
  ~ScoredSeqWithCarriers();

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
  /* see parent abstract class ScoredSeqNested */
  bool isNested();
  ScoredSeq* getNested();
  void deepDelete();
  void deepUnbuffer();
  /* see parent abstract class ScoredSeqPaired */
  void makeAlive();
  bool isAlive();

  // class-specific
  void addCarriedSeq(ScoredSeq* seq, int setNum);
  void addCarriedSeqs(set<ScoredSeq*>* seqs, int setNum);
  void removeCarriedSeq(ScoredSeq* seq, int setNum);
  void getCarriedSeqs(set<ScoredSeq*>* seqs, int setNum);
  void clearCarriedSeqs(int setNum);
  long numCarriedSeqs(int setNum);
  int numCarriers();
  bool hasCarriedSeq(ScoredSeq* seq, int setNum);


 protected:
  friend class ReadFileCommunicator;
  void addPairedEnd(ScoredSeqWithCarriers * ss);

 private:
  /* informs this that its pair is dead; authenticity is verified by
   * passing the dying Read (oldPair) and its destructable replacement. 
   * note that passing out the paired end again should restore its
   * indestructable status */
  void makeDead();
  void replaceDeadPe(ScoredSeqWithCarriers * deadReplacement);

  void shallowDelete(); // does not delete a dead PE
  ScoredSeq* _innerSeq;
  int _numOfSets;
  set<ScoredSeq*>** _carriedSeqSets;
  ScoredSeqWithCarriers* _pe;
  bool _thisIsAlive;
  bool _okToDeletePe;

};

#endif
