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
This is the implementation for Fasta-format read files:
"""
>SCS:2:1:8:926#0/1_1
GAAACACAGTCAAAATAAATGAAGAG
GCGCA
"""
Sequences may contain A, T, C, and G nucleotides, or Ns.  They may also contain
other nucleotide IUPAC characters.  U will be interpreted as T, all other letters
will count as N.  Gap characters or any other character will simply be eliminated
without an error being raised.  ">" is not a valid character in the sequences.
Uppercase and lowercase are both valid.
 */


#ifndef READFILEFASTASINGLE_H
#define READFILEFASTASINGLE_H

#include <vector>
#include <string>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include "ReadFile.h"
#include "ReadFileCommunicator.h"
#include "Read.h"
#include "ScoredSeq.h"
#include "ReadFileIndex.h"
using namespace std;


class ReadFileFastaSingle : public ReadFile {
 public:

  /* Non-paired end reads will be returned; every nucleotide will have the given score;
   * default score for first constructor is 1.0.
   */
  ReadFileFastaSingle(string filename);
  ReadFileFastaSingle(string filename, float nucScore, bool invert=false);
  // for false-paired reads (invert must be specified since one of the two file objects
  // must be inverted)
  ReadFileFastaSingle(string filename, float nucScore, long readSize, bool invert);
  ~ReadFileFastaSingle();

  /* see parent class ReadFile */
  ReadFileFastaSingle* copy();
  void splitFile(set<ReadFile*>* putFilesHere, int numFiles);
  ReadFile* subFile(long numReadsSkip, long numReadsKeep);
  string getName();
  void open();
  void open(bool getBothEnds);
  void openWithCarriers(int numSets);
  void openWithCarriers(int numSets, bool getBothEnds);
  void openNormalized();
  void openNormalized(bool getBothEnds);
  void openNormWithCarriers(int numSets);
  void openNormWithCarriers(int numSets, bool getBothEnds);
  bool hasRead();
  ScoredSeq* getRead();
  void skipRead();
  void skipReads(long numToSkip);
  void close();
  long ampliconSize();
  long numReads();

 protected:
  friend class ReadFileCommunicator;
  void bufferQueue(Read * r);
  void bufferUnqueue(Read * r);
  void bufferFill();
  void bufferFill(BufferThreadedness threadedness);

 private:
  void constructorHelp();

  long bufferFillSetupHelper(char** rawSeqArray, long* rfiReadSize, Read** readArray);
  void bufferFillInterpretHelper(long seqN, char** rawSeqArray, long* rfiReadSize,
				 Read** readArray, ScoredSeq** bufferArray);
  void bufferFillInterpretSplitHelper(long seqN, char** rawSeqArray, long* rfiReadSize,
				      Read** readArray, ScoredSeq** bufferArray);

  // the nucleotide score
  float _nucScore;
  // if true, sequences are RC'ed
  bool _invert;
  // if true, the second value is used to split each single read into
  // two "paired-end" reads of the specified size
  bool _splitMode;
  long _maxReadSize;
  // for string interpretation
  set<char> _okToSkip;


  // private constructor for making subset files
  ReadFileFastaSingle(ReadFileFastaSingle* precursor, int intialBlock, long initialPos, long numReads, bool numReadsDetermined);

  void updateBuffer();
  void emptyBuffer();

  bool _openWithCarriers;
  int _numCarrierSets;
  bool _openNormalized;

  ScoredSeq* getReadHelper();
  void readFromFileHelper();
  bool _isOpen;
  ifstream _getReadsFile; // the file object for obtaining reads
  queue<ScoredSeq*> _getReadsQueue;


  char* _filename;
  int _minCountedScore;

  set<Read*> _bufferSet;

  long _numReads;
  // this will be true once the number of reads has been counter
  // or if it was provided to the private constructor
  bool _numReadsDetermined;

  // these are used during the process of iterating through the file
  bool _hasReadWaiting;
  long _seqStart;
  long _seqLen;
  int _currentBlock;
  long _currentPos;

  class BlockPosFileCarrier {
  public:
    BlockPosFileCarrier();
    BlockPosFileCarrier(ReadFileFastaSingle* host);
    ~BlockPosFileCarrier();
    void seekBlockPos(int newBlock, long newPos);
    void getString(char* cstringToFill, streamsize length);
  private:
    ReadFileFastaSingle* _host;
    ifstream _file;
    int _block;
    long _pos;
  };

  // stuff for dealing with the differnt long int-limited sections 
  // of the file.  REP INVARIANT: _numPosBlocks
  long* _posBlocks;
  int _numPosBlocks;
  // these define any initially skipped part of the file
  // and are by default set to zero by all the public constructors
  int _initialBlock;
  long _initialPos;
  // identifies the last read in a file if _numReadsDetermined==true
  long _numReadsRead;
  // defines the max _pos value before a new block is created
  // default is MAX LONG
  static long _maxBlockSize;


  class ReadFileIndexJustSeq : public ReadFileIndex {
  public:
    ReadFileIndexJustSeq(); //default constructor
    ReadFileIndexJustSeq(long readStart, long readSize, int blockNum); //constructor
    ~ReadFileIndexJustSeq();
    long readStart();
    long readSize();
    long scoreStart(); // meaningless but must be implemented; throws an exception if called
    long scoreSize(); // meaningless but must be implemented; throws an exception if called
    long linkStart(); // meaningless but must be implemented; throws an exception if called
    long linkSize(); // meaningless but must be implemented; throws an exception if called
    int blockNum();
  private:
    long _readStart;
    long _readSize;
    int _blockNum;
  };

  struct SortReadByRfiStart {
    bool operator() (Read* readA, Read* readB);
  };

};

#endif
