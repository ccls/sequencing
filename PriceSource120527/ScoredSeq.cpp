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


#ifndef SCOREDSEQ_CPP
#define SCOREDSEQ_CPP

#include "ScoredSeq.h"
#include "ScoredSeqMonoScore.h"
#include "ScoredSeqShallow.h"
#include "AssemblyException.h"
#include <iostream>
using namespace::std;

ScoredSeq::~ScoredSeq() {}


// static method
ScoredSeq* ScoredSeq::getScoredSeq(char* seq, float score, long size){
  return new ScoredSeqMonoScore(seq,score,size);
}
ScoredSeq* ScoredSeq::getScoredSeq(char* seq, float* scores, long size){
  return new ScoredSeqShallow(seq,scores,size);
}
ScoredSeq* ScoredSeq::getScoredSeq(char* seq, float* scores, float* links, long size){
  return new ScoredSeqShallow(seq,scores,links,size);
}

ScoredSeq* ScoredSeq::repExposedSeq(char* seq, float score, long size){
  return new ScoredSeqMonoScore(seq,score,size,true);
}
ScoredSeq* ScoredSeq::repExposedSeq(char* seq, float* scores, long size){
  float* links = new float[size+1];
  for (long n = 0; n < size; ++n){ links[n] = 1; }
  return new ScoredSeqShallow(true,seq,scores,links,size);
}
ScoredSeq* ScoredSeq::repExposedSeq(char* seq, float* scores, float* links, long size){
  return new ScoredSeqShallow(true,seq,scores,links,size);
}

ScoredSeq* ScoredSeq::copyShallowSeq(ScoredSeq* seq, char sense){
  return new ScoredSeqShallow(seq,sense);
}



char* ScoredSeq::reverseComplement(char* seq, long length){
  char* rcSeq = new char[length+1];
  for (long n = 0; n < length; n++){
    switch(seq[n]) {
    case 'A': rcSeq[length - n - 1] = 'T'; break;
    case 'C': rcSeq[length - n - 1] = 'G'; break;
    case 'G': rcSeq[length - n - 1] = 'C'; break;
    case 'T': rcSeq[length - n - 1] = 'A'; break;
    case 'N': rcSeq[length - n - 1] = 'N'; break;
    default:
      /*
      char errorMessage[100];
      char part1[] = "invalid nucleotide:";
      sprintf(errorMessage, "%s %c", part1, n);
      throw AssemblyException::ArgError(errorMessage);
      */
      throw AssemblyException::ArgError("invalid nucleotide in ScoredSeq::reverseComplement");
    }
  }
  rcSeq[length] = '\0'; // the null ending that makes it a cstring
  return rcSeq;
}


#endif
