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
 */

#ifndef SCOREDSEQSHALLOW_H
#define SCOREDSEQSHALLOW_H

#include "ScoredSeq.h"
#include "AssemblyException.h"
#include "AssemblyException.h"

class ScoredSeqShallow : public ScoredSeq {
 public:
  ScoredSeqShallow(); //default constructor
  ScoredSeqShallow(string seq, vector<float>* scores); // sets all linkages to 1
  ScoredSeqShallow(string seq, vector<float>* scores, vector<float>* links);
  // REQUIRES that len(scores)==len(seq) AND len(links)==len(seq)-1
  ScoredSeqShallow(string seq, float* scores); // sets all linkages to 1
  ScoredSeqShallow(string seq, float* scores, float* links);
  // REPLACES strings with char* arrays
  ScoredSeqShallow(char* seq, float* scores, long size);
  ScoredSeqShallow(char* seq, float* scores, float* links, long size);
  // ALLOWS rep to be exposed (i.e. arrays stored are the actual ones input)
  ScoredSeqShallow(bool expRep, char* seq, float* scores, float* links, long size);


  ScoredSeqShallow(ScoredSeq *readToCopy, char sense); //constructor
  void OK();
  ~ScoredSeqShallow();

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
  void deepDelete();
  void deepUnbuffer();

 private:
  char* _seq;
  long _size;
  char* _rcSeq;
  float* _scores;
  float* _links;
  ScoredSeqShallow * _pe;

  void constructorSeqHelper(string seq);

};

#endif
