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
This is a read from a ReadSet.  This base class will be used longernally
by the ReadSet and will be extended to include the features of a ScoredSeq.

SPEC FIELDS:
readSet: the ReadSet from which this BareRead derives
specifier: an identifier to distinguish this BareRead from other
   BareReads from the same ReadSet

ABSTRACTION FUNCTION:
ReadSet readSetVar -> a reference to the instance of ReadSet from which 
   this BareRead derives
long specifierVar -> the unique id for this BaseRead within readSet

REP. INVARIANT:
-this BareRead must be a member of readSet

-because this is carrying a polonger, the copy constructor, copy assignment
 operator, and destructor must all be implemented.
 */

#ifndef READ_H
#define READ_H

#include "ReadFileIndex.h"
#include "ReadFile.h"
#include "ScoredSeqPaired.h"
#include "ScoredSeqShallow.h"
#include "ReadFileCommunicator.h"

class Read : public ScoredSeqPaired {
 public:
  Read(); //default constructor
  Read(ReadFile *rs, ReadFileIndex *rfi); //constructor
  ~Read();

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

  /* see parent abstract class ScoredSeq 
     NON-SPECIFICATION NOTE: in this class, buffer is used;
     to reduce memory, sequences/scores are not stored and
     therefore need to be obtained from the source ReadFile
     when needed; to expedite that process, the retrieval
     tasks are bundled by having called 'buffer' ahead of
     time.  This is not a requirement, getSeq etc should still
     work even if buffer has not been called. */
  void buffer();
  void unbuffer();
  void bottomBuffer();

  /* see parent abstract class ScoredSeq */
  long size();
  bool hasPairedEnd();
  ScoredSeq * getPairedEnd();
  ScoredSeq * getTempPairedEnd();
  ScoredSeq * shallowCopy();
  ScoredSeq * flipCopy();
  bool isNested();
  void deepDelete();
  void deepUnbuffer();


  /* see parent abstract class ScoredSeqPaired */
  void makeAlive();
  bool isAlive();

 protected:
  friend class ReadFileCommunicator;
  // this method can only be called by a ReadFile and
  // is only called by a ReadFile upon request from a
  // Read, and those requests are all found in critical blocks.
  void addScoredSeqImp(ScoredSeq* ssi);
  void addPairedEnd(Read * ss);
  ReadFile * getReadFile();
  ReadFileIndex * getRfiRef();


 private:
  // these are used internally to manage state and are therefore
  // not critical, though they should only be called from within
  // critical blocks
  void bufferPrivate();
  void unbufferPrivate();
  // these do not have public equivalents but are labelled as
  // private because they also should not be called outside a
  // critical block
  void fillBufferPrivate();
  void unqueuePrivate();

  /* the copy constructor is private because it is a shallow copy (returns
   * a copy with copies of the references, not the referenced objects), so
   * copies will not work after original has been deleted.  ditto for the 
   * copy assignment operator. both raise exceptions if called. */
  Read(const Read& r);  //copy constructor
  Read& operator=(const Read& that); //copy assignment operator

  /* informs this that its pair is dead; authenticity is verified by
   * passing the dying Read (oldPair) and its destructable replacement. 
   * note that passing out the paired end again should restore its
   * indestructable status */
  void makeDead();
  void replaceDeadPe(Read * deadReplacement);

  void shallowDelete(); // does not delete a dead PE
  ReadFile* _readFile;
  ReadFileIndex* _rfi;
  ScoredSeq* _buffer;
  bool _requestedBuffer;
  Read* _pe;
  bool _thisIsAlive;
  bool _okToDeletePe;

};

#endif
