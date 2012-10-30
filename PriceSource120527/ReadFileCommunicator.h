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
Permits privelidged communication between Read and all of ReadFile's subclasses.
This class should never be used outside of the ReadFactory facade.  It breaks the
modularity of the Read class from the module of the assembler.  In any case, the
external client should not be able to use this method on reads returned by the 
ReadFactory because those are cast as ScoredSeqs, and the methods here require
that they be cast as the ScoredSeq subclass Read.

In order for this class to succeed in simplifying the implementation of read-relevant
classes, any subclass/implementation of Read or ReadFile should friend this class.
If that is done, neither this class nor any of the other Read/ReadFile classes need
to be modified when a new class is added since all of the communication betweeen 
them is handled by this class' public methods.

Just static methods.  No need to instantiate.
 */


#ifndef READFILECOMMUNICATOR_H
#define READFILECOMMUNICATOR_H

#include <vector>
#include <string>
#include "Read.h"
#include "ScoredSeqMonoScore.h"
#include "ScoredSeqShallow.h"
#include "ScoredSeqNormalized.h"
#include "ScoredSeqWithCarriers.h"
//#include "ScoredSeqAlignCarriers.h"
class ReadFile;
class ReadFileIndex;
using std::string;
using std::vector;


class ReadFileCommunicator {
 public:


  /* Adds r to the queue of Reads that will be seq/score buffered by 
     ReadFile rf as soon as the ScoredSeq methods are called on any
     of the Reads in the buffer queue.  This method must be called 
     prior to bufferFill in order Read r to be buffered.
     THROWS: ???? if r is not a member of rf
     MODIFIES: rf by adding r to the bufferQueue
  */
  static void bufferQueue(ReadFile * rf, Read * r);
  static void bufferUnqueue(ReadFile * rf, Read * r);

  // factory method that will allow the Read class to be 
  // decoupled from the read file classes (leaving only
  // a weak dependency)
  static Read* makeRead(ReadFile *rs, ReadFileIndex *rfi); //constructor

  /* Actually gets the sequences/scores for the buffered Reads and
     adds them to the Reads' internal memory.
     MODIFIES: rf by removing all queued Reads from the bufferQueue
     MODIFIES: all Reads in the rf bufferQueue by calling addScoredSeqImp
     on them.
  */
  static void bufferFill(ReadFile * rf);
  // WARNING: only call a threaded buffer from a non-threaded block;
  // buffering is non-critical in that case
  static void bufferFill(ReadFile * rf, ReadFile::BufferThreadedness threadedness);


  /* Places the actual sequence/score info into a Read so that it will
     then be accessible by Read's ScoredSeq methods.  This step is a 
     requirement for Read to succeed in fulfilling its specifications
     but is an internal, implementation-specific requirement; the abstract
     state of Read does not change, but this function will need to be
     called in order for Read to meet its specifications.
     MODIFIES: Abstract State: None
               Concrete State: Adds seq/scores to Read r
     REQUIRES: ssi is the correct seq/scores for r (Read has no way of 
               checking this)
  */
  static void addScoredSeqImp(Read* r, ScoredSeq* ssi);


  /* Accesses a private piece of data from the Read that can
     be used by the ReadFile to gather the sequence/scores for
     the read.  Should only be called by a ReadFile on a member
     read.
     REQUIRES: r is a member of the calling ReadFile (this can be
     checked ahead of time with the doesFileHaveRead method).
  */
  static ReadFileIndex * getRfiRef(Read * r);


  /* Pairs paired-end reads.  This must be done after the reads are created
     but should not be done by clients, which is why the function is here.
     MODIFIES: rA and rB by pairing them with each other (hasPairedEnd->true)
     THROWS: ???? if rA and rB do not derive from the same ReadFile
  */
  static void pairEnds(Read * rA, Read * rB);
  static void pairEnds(ScoredSeqWithCarriers * rA, ScoredSeqWithCarriers * rB);
  static void pairEnds(ScoredSeqNormalized * rA, ScoredSeqNormalized * rB);
  static void pairEnds(ScoredSeqMonoScore * rA, ScoredSeqMonoScore * rB);
  //static void pairEnds(ScoredSeqAlignCarriers * rA, ScoredSeqAlignCarriers * rB);


  /* sort of a checkRep (OK) for making sure that the Read and
     its ReadFile are correctly matched. */
  static void checkFileHasRead(ReadFile * rf, Read * r);

};

#endif
