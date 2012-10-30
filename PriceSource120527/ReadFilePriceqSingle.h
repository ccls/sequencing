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
 * special format for PRICE input/ouput
 *
 * Six-line convention:
 * Line 1: "@" then sequence name
 * Line 2: DNA sequence
 * Line 3: "+" then optionally sequence name
 * Line 4: ASCII-encoded nucleotide scores (length == length of DNA sequence)
 * Line 5: "~" then optionally sequence name
 * Line 6: ASCII-encoded linkage scores (length == length of DNA sequence minus one)
 *
 * Scores are supposed to reflect the PRICE internally-tracked sequence scores (which
 * are supposed to reflect supporting read counts), but are lower-resolution for compactness.
 *
 * AD = decimal value of ASCII character
 *     _AD_      _Score_
 *      33         0.0
 *      34         0.3
 *     >=35     1.5 ^ (AD - 35)
 *     >=126    1.5 ^ (126 - 35)
 * NOTE: in the score->AD conversion, any score >= 0.9 will be upgraded to 35
 * max score (calculated by python): >>> 1.5**(126-35) = 658887371.45190775
 */


#ifndef READFILEPRICEQSINGLE_H
#define READFILEPRICEQSINGLE_H

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


class ReadFilePriceqSingle : public ReadFile {
 public:

  /* Non-paired end reads will be returned; the whole line for each entry
     will be treated as a read.
  */
  ReadFilePriceqSingle(string filename);
  ReadFilePriceqSingle(bool invert, string filename);

  /* Non-paired end reads will be returned.  Only gets the portion of the read
     that is described by the specified readStart and readLength.  readStart
     is zero-indexed from the 5p end of the read; readLength is in units of
     nucleotides.
     THROWS: ???? if readStart + readLength > length of any read in the file.
     * NOTE: zero read length means the whole read length, whatever that is
  */
  ReadFilePriceqSingle(string filename, int readStart, int readLength);
  ReadFilePriceqSingle(bool invert, string filename, int readStart, int readLength);
  ~ReadFilePriceqSingle();

  /* see parent class ReadFile */
  ReadFilePriceqSingle* copy();
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
  void constructorHelper(string filename);
  float* cstringToScores(char* scoreCstring, long readSize);

  long bufferFillSetupHelper(char** rawSeqArray, char** rawScoresArray, char** rawLinksArray,
			     long* rfiReadSize, Read** readArray);
  void bufferFillInterpretHelper(long seqN, char** rawSeqArray, char** rawScoresArray, char** rawLinksArray, 
				 long* rfiReadSize, Read** readArray, ScoredSeq** bufferArray);

  // private constructor for making subset files
  ReadFilePriceqSingle(ReadFilePriceqSingle* precursor, int intialBlock, long initialPos, long numReads, bool numReadsDetermined);

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
  int _readStart; // zero-indexed
  int _readLength; // -1 if read length is not fixed; >=0 otherwise
  bool _invert;

  // for string interpretation
  set<char> _okToSkip;

  set<Read*> _bufferSet;

  long _numReads;
  // this will be true once the number of reads has been counter
  // or if it was provided to the private constructor
  bool _numReadsDetermined;

  // these are used during the process of iterating through the file
  long _seqStart;
  long _seqLen;
  long _scoreStart;
  long _scoreLen;
  long _linkStart;
  long _linkLen;
  bool _hasReadWaiting;
  int _currentBlock;
  long _currentPos;

  class BlockPosFileCarrier {
  public:
    BlockPosFileCarrier();
    BlockPosFileCarrier(ReadFilePriceqSingle* host);
    ~BlockPosFileCarrier();
    void seekBlockPos(int newBlock, long newPos);
    void getString(char* cstringToFill, streamsize length);
  private:
    ReadFilePriceqSingle* _host;
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



  class ReadFileIndex3start1len : public ReadFileIndex {
  public:
    ReadFileIndex3start1len(); //default constructor
    ReadFileIndex3start1len(long readStart, long scoreStart, long linkStart, long readSize, int blockNum);
    ~ReadFileIndex3start1len();
    long readStart();
    long readSize();
    long scoreStart();
    long scoreSize();
    long linkStart();
    long linkSize();
    int blockNum();
  private:
    long _readStart;
    long _readSize;
    long _scoreStart;
    long _linkStart;
    int _blockNum;
  };

  struct SortReadByRfiStart {
    bool operator() (Read* readA, Read* readB);
  };

};

#endif
