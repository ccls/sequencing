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
This is the implementation for Fastq-format read files.  Really, it was
implemented for Illumina's _sequence.txt files, which follow the fastq
format.  Each read in these files is given four newline-separated lines
in the file:
"""
@SCS:2:1:8:926#0/1_1
GAAACACAGTCAAAATAAATGAAGAG
+SCS:2:1:8:926#0/1_1
^\aa]V`W^\`TQ`\]`bab^aa`aa
"""
The first and third lines of each quatrain are the name of the read; names
should be the same between those two postions and unique within a file.
In the first line, the name is preceded by "@"; in the third, it is preceded
by "+".  The second line contains the nucleotide sequence in capital letters:
A, T, C, G, and N (meaning 'uncalled') are supported.  The fourth line contains
a string of ascii characters whose 8-bit integer values are confidence scores
for the bases at the same position in the second line (the second and fourth lines
must be of equal length).  Higher scores have higher confidence (lower probability
of a miscall).  Scores are assumed to have a negative linear relationship with the
log of the probability that the base is miscalled.  Ns should have a low score that
reflects the zero probability that 'N' is the correctly-identified nucleotide.

This class has two modes for interpreting scores: fastqSystem and illuminaSystem.
Both interpret the ASCII characters as Phred scores, with the fastqSystem zero
Phred score set to 33 and the illuminaSystem zero Phred score set to 64.  For
both systems, scores of 2 or less are interpreted as a zero score in PRICE.
this is the specified interpretation of illumina 1.5+ and a heuristic measure
for earlier illumina and sanger formats.
 */


#ifndef READFILEFASTQSINGLE_H
#define READFILEFASTQSINGLE_H

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


class ReadFileFastqSingle : public ReadFile {
 public:

  enum Encoding{ fastqSystem=0, illuminaSystem=1 };

  ReadFileFastqSingle();

  /* Non-paired end reads will be returned; the whole line for each entry
     will be treated as a read.
  */
  ReadFileFastqSingle(string filename, Encoding encoding);
  ReadFileFastqSingle(string filename, Encoding encoding, float countFactor);
  ReadFileFastqSingle(bool invert, string filename, Encoding encoding, float countFactor);

  /* Non-paired end reads will be returned.  Only gets the portion of the read
     that is described by the specified readStart and readLength.  readStart
     is zero-indexed from the 5p end of the read; readLength is in units of
     nucleotides.
     THROWS: ???? if readStart + readLength > length of any read in the file.
     * NOTE: zero read length means the whole read length, whatever that is
  */
  ReadFileFastqSingle(string filename, Encoding encoding, int readStart, int readLength);
  ReadFileFastqSingle(string filename, Encoding encoding, int readStart, int readLength, float countFactor);
  ReadFileFastqSingle(bool invert, string filename, Encoding encoding, int readStart, int readLength, float countFactor);

  // this is the false-pair constructor (since invert is required)
  ReadFileFastqSingle(bool invert, long readLength, string filename, Encoding encoding, float countFactor);

  ~ReadFileFastqSingle();

  /* see parent class ReadFile */
  ReadFileFastqSingle* copy();
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
  void constructorConstHelper();


  long bufferFillSetupHelper(char** rawSeqArray, char** rawScoresArray, long* rfiReadSize, Read** readArray);
  void bufferFillInterpretHelper(long seqN, char** rawSeqArray, char** rawScoresArray, long* rfiReadSize,
				 Read** readArray, ScoredSeq** bufferArray);
  void bufferFillInterpretSplitHelper(long seqN, char** rawSeqArray, char** rawScoresArray, long* rfiReadSize,
				      Read** readArray, ScoredSeq** bufferArray);

  // private constructor for making subset files
  ReadFileFastqSingle(ReadFileFastqSingle* precursor, int intialBlock, long initialPos, long numReads, bool numReadsDetermined);
  void setEncodingHelper(Encoding encoding, float countFactor);



  // SCORE CONVERSION ENCODING
  // these will contain accelerated methods for converting score chars into floats; 
  // i will pre-compute the values and use switch statements
  class ScoreConverter {
  public:
    virtual ~ScoreConverter();
    // positions in the seqString with close-to-zero scores will be converted to 'N' chars
    // I will cut it off at less-than-one-bit-of-information
    virtual float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString) = 0;
    virtual Encoding getEncoding() = 0;
  };
  // Phred+33 encoding
  class ScoreConverterFastq : public ScoreConverter {
  public:
    ScoreConverterFastq(float countFactor);
    ~ScoreConverterFastq();
    float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString);
    Encoding getEncoding();
  private:
    float charToScore(char nuc, char* seqString, long seqPos);
    float _scores[100];
  };
  // Phred+64 encoding
  class ScoreConverterIllumina : public ScoreConverter {
  public:
    ScoreConverterIllumina(float countFactor);
    ~ScoreConverterIllumina();
    float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString);
    Encoding getEncoding();
  private:
    float charToScore(char nuc, char* seqString, long seqPos);
    float _scores[100];
  };

  // encoding for Illumina or Sanger fastq scores
  ScoreConverter* _encodingConverter;


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

  // these are mode controllers
  char* _filename;
  int _readStart; // zero-indexed
  int _readLength; // -1 if read length is not fixed; >=0 otherwise
  float _countFactor; // every read count is multiplied by this
  bool _invert;
  // if true, the second value is used to split each single read into
  // two "paired-end" reads of the specified size
  bool _splitMode;
  long _maxReadSize;

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
  bool _hasReadWaiting;
  int _currentBlock;
  long _currentPos;

  class BlockPosFileCarrier {
  public:
    BlockPosFileCarrier();
    BlockPosFileCarrier(ReadFileFastqSingle* host);
    ~BlockPosFileCarrier();
    void seekBlockPos(int newBlock, long newPos);
    void getString(char* cstringToFill, streamsize length);
  private:
    ReadFileFastqSingle* _host;
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


  // the RFI class for this file type
  class ReadFileIndex2start1len : public ReadFileIndex {
  public:
    ReadFileIndex2start1len(); //default constructor
    ReadFileIndex2start1len(long readStart, long scoreStart, long readSize, int blockNum); //constructor
    ~ReadFileIndex2start1len(); //default constructor
    long readStart();
    long readSize();
    long scoreStart();
    long scoreSize();
    long linkStart(); // meaningless but must be implemented; throws an exception if called
    long linkSize(); // meaningless but must be implemented; throws an exception if called
    int blockNum();
  private:
    long _readStart;
    long _readSize;
    long _scoreStart;
    int _blockNum;
  };

  struct SortReadByRfiStart {
    bool operator() (Read* readA, Read* readB);
  };

};

#endif
