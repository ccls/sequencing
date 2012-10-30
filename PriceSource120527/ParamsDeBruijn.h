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

/* packages paramaters that will be used by Assembler to decide
 * how to execute de Bruijn graph assemblies
 */


#ifndef PARAMSDEBRUIJN_H
#define PARAMSDEBRUIJN_H

#include "ReadFile.h"
#include "ScoredSeq.h"
#include <set>

class ParamsDeBruijn {

 public:
  ParamsDeBruijn();
  ParamsDeBruijn(int kmerSize);
  ParamsDeBruijn(int kmerSize, float minCount);
  ParamsDeBruijn(int kmerSize, float minCount, long maxSeqLength);
  ParamsDeBruijn(int kmerSize, float minCount, long maxSeqLength, long minNumSeqs);
  ParamsDeBruijn(ParamsDeBruijn* pcm);
  ~ParamsDeBruijn();

  static float defaultMinCount();
  static int defaultKmerSize();
  static long defaultMaxSeqLength();
  static long defaultMinNumSeqs();
  static ParamsDeBruijn* getDefault();

  // get the values
  float minCount();
  int kmerSize();
  long maxSeqLength();
  long minNumSeqs();

  // set the values
  void setMinCount(float minCount);
  void setKmerSize(int kmerSize);
  void setMaxSeqLength(long maxSeqLength);
  void setMinNumSeqs(long minNumSeqs);

  // check to see if a sequence should go into DeBruijn assembly
  bool inputToDeBruijn(ScoredSeq* seq);
  // check to see if an output seq should be kept
  bool retainDeBruijnOutput(ScoredSeq* seq);

  ParamsDeBruijn* copy();

 private:
  int _kmerSize;
  float _minCount;
  long _maxSeqLength;
  // the number of input seqs required for de Bruijn to be invoked
  long _minNumSeqs;
};

#endif
