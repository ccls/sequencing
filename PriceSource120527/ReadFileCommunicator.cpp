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

#ifndef READFILECOMMUNICATOR_CPP
#define READFILECOMMUNICATOR_CPP

#include "ReadFileCommunicator.h"
#include "AssemblyException.h"
#include <omp.h>


Read* ReadFileCommunicator::makeRead(ReadFile *rs, ReadFileIndex *rfi){
  return new Read(rs, rfi);
}


// sort of a checkRep for input args
void ReadFileCommunicator::checkFileHasRead(ReadFile * rf, Read * r){
  if ( rf != r->getReadFile() ){
    throw AssemblyException::ArgError("Read does not derive from provided ReadFile.");
  }
}


// pairs ends
void ReadFileCommunicator::pairEnds(Read * rA, Read * rB){
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
}
void ReadFileCommunicator::pairEnds(ScoredSeqWithCarriers * rA, ScoredSeqWithCarriers * rB){
  if ( rA->getNested()->hasPairedEnd() ){
    if (! rB->getNested()->hasPairedEnd() ){
      throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
    } else if ( rA->getNested()->getTempPairedEnd() != rB->getNested() ){
      throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired to each other: SSWC");
    }
  } else if ( rB->getNested()->hasPairedEnd() ){
    throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
  } else {
    throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired at all: SSWC");
  }
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
}
void ReadFileCommunicator::pairEnds(ScoredSeqNormalized * rA, ScoredSeqNormalized * rB){
  if ( rA->getNested()->hasPairedEnd() ){
    if (! rB->getNested()->hasPairedEnd() ){
      throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
    } else if ( rA->getNested()->getTempPairedEnd() != rB->getNested() ){
      throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired to each other: SSWC");
    }
  } else if ( rB->getNested()->hasPairedEnd() ){
    throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
  } else {
    throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired at all: SSWC");
  }
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
}
// this is to help out with tests
void ReadFileCommunicator::pairEnds(ScoredSeqMonoScore * rA, ScoredSeqMonoScore * rB){
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
}
/*
void ReadFileCommunicator::pairEnds(ScoredSeqAlignCarriers * rA, ScoredSeqAlignCarriers * rB){
  if ( rA->getNested()->hasPairedEnd() ){
    if (! rB->getNested()->hasPairedEnd() ){
      throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
    } else if ( rA->getNested()->getTempPairedEnd() != rB->getNested() ){
      throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired to each other: SSWC");
    }
  } else if ( rB->getNested()->hasPairedEnd() ){
    throw AssemblyException::LogicError("Can't pair in RFC; one nested seq is paired, one isn't: SSWC");
  } else {
    throw AssemblyException::LogicError("Can't pair in RFC; nested seqs aren't paired at all: SSWC");
  }
  rA->addPairedEnd(rB);
  rB->addPairedEnd(rA);
}
*/

// used by Read to communicate with ReadFile

void ReadFileCommunicator::bufferQueue(ReadFile * rf, Read * r){
  checkFileHasRead(rf,r);
  rf->bufferQueue(r);
}
void ReadFileCommunicator::bufferUnqueue(ReadFile * rf, Read * r){
  checkFileHasRead(rf,r);
  rf->bufferUnqueue(r);
}

void ReadFileCommunicator::bufferFill(ReadFile * rf){
  rf->bufferFill();
}

void ReadFileCommunicator::bufferFill(ReadFile * rf, ReadFile::BufferThreadedness threadedness){
  rf->bufferFill(threadedness);
}



// used by ReadFile to communicate with Read
void ReadFileCommunicator::addScoredSeqImp(Read * r, ScoredSeq* ssi){
  r->addScoredSeqImp(ssi);
}

ReadFileIndex * ReadFileCommunicator::getRfiRef(Read * r){
  return r->getRfiRef();
}





#endif
