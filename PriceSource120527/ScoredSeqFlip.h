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
 * This is a wrapper for a ScoredSeq that just flips the sequence
 * over.  It is just supposed to speed up some operations.  It can
 * also wrap a sequence without flipping it.
 */

#ifndef SCOREDSEQFLIP_H
#define SCOREDSEQFLIP_H

#include "ScoredSeqNested.h"
#include "AssemblyException.h"
#include "AssemblyException.h"

class ScoredSeqFlip : public ScoredSeqNested {
 public:
  static ScoredSeqFlip* getFlip(ScoredSeq* innerSeq, char sense);
  virtual ~ScoredSeqFlip();
  // class-specific; other methods inherited from implemented classes
  virtual char getSense() = 0;
};


class ScoredSeqFlipPlus : public ScoredSeqFlip {
 public:
  ScoredSeqFlipPlus(ScoredSeq* innerSeq);
  ~ScoredSeqFlipPlus();

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

  // class-specific
  void OK();
  char getSense();

 private:
  ScoredSeq* _innerSeq;

};


class ScoredSeqFlipMinus : public ScoredSeqFlip {
 public:
  ScoredSeqFlipMinus(ScoredSeq* innerSeq);
  ~ScoredSeqFlipMinus();

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

  // class-specific
  char getSense();

 private:
  ScoredSeq* _innerSeq;
};

#endif
