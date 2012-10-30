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

/* For when the scores at all positions are the same.  Does not support paired-ends right now.
 */

#ifndef SCOREDSEQMONOSCORE_H
#define SCOREDSEQMONOSCORE_H

#include "ScoredSeqPaired.h"


class ScoredSeqMonoScore : public ScoredSeqPaired {
 public:
  ScoredSeqMonoScore(string seq, float score);
  ScoredSeqMonoScore(char* seq, float score, long size, bool expRep=false); 
  ScoredSeqMonoScore(ScoredSeq* seq, float score);
 ~ScoredSeqMonoScore();

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
  void deepDelete();
  void deepUnbuffer();

  void OK();

  void makeAlive();
  bool isAlive();


 protected:
  friend class ReadFileCommunicator;
  void addPairedEnd(ScoredSeqMonoScore * ss);


 private:
  ScoredSeqMonoScore(ScoredSeqMonoScore* seq, char sense);  //copy constructor
  void constructorSeqCopyHelper(char* seq);


  /* informs this that its pair is dead; authenticity is verified by
   * passing the dying ScoredSeq (oldPair) and its destructable replacement. 
   * note that passing out the paired end again should restore its
   * indestructable status */
  void makeDead();
  void replaceDeadPe(ScoredSeqMonoScore * deadReplacement);

  char* _seq;
  long _size;
  char* _rcSeq;
  float _score;
  ScoredSeqMonoScore* _pe;
  bool _thisIsAlive;

  void constructorSeqHelper(string seq);

};

#endif
