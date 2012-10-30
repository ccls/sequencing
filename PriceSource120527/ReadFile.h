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
This is an abstract class for a file of reads.  I will create an implementation
for each format of input file, paired-end or not, etc.
 */


#ifndef READFILE_H
#define READFILE_H

# include <vector>
# include <string>
# include <set>
class Read;
class ScoredSeq;
class ReadFileIndex;
using namespace::std;


class ReadFile {
 public:


  // returns a string with the file type based on the append of the filename
  // NONEFILE: no described filetype
  // FASTAFILE: .fa or .fasta
  // FASTQFILE: .fq or .fastq or _sequence.txt
  enum FileType { NONEFILE, FASTAFILE, FASTQFILE, PRICEQFILE, ILLUMINAFILE };
  static FileType getFileType(string filename);

  // factory methods
  static ReadFile* makeBasicReadFile(string filename, float countFactor = 1);
  static ReadFile* makeInvertedReadFile(string filename, float countFactor = 1);
  // assumes that reads face towards one another
  static ReadFile* makePairedReadFile(string filename, long ampliconSize, float countFactor = 1);
  static ReadFile* makePairedReadFile(string filenameA, string filenameB, long ampliconSize, float countFactor = 1);
  // assumes that reads face away from one another in the file; files return RC seqs
  static ReadFile* makeMatePairFile(string filename, long ampliconSize, float countFactor = 1);
  static ReadFile* makeMatePairFile(string filenameA, string filenameB, long ampliconSize, float countFactor = 1);
  // treats each read as a pair of reads by taking fragments that face one another and halving the
  // scores of overlapping nucleotides
  static ReadFile* makeFalsePairFile(string filename, long readSize, long ampliconSize, float countFactor = 1);

  // these define an interface for use and re-use of read files
  virtual ~ReadFile();

  // the name should allow the user to know where the file(s) is(are),
  // but it need not be directly interpretable by the operating system
  // (that would be impossible for ReadFile objects that include multiple
  // actual files.
  virtual string getName() = 0;
  virtual ReadFile* copy() = 0;

  // will split up the file contents as evenly as possible between the different
  // files.  up to numFiles will be created, but could be fewer; empty files
  // will not neccessarily be added to the collection (but could be).
  virtual void splitFile(set<ReadFile*>* putFilesHere, int numFiles) = 0;
  // REQUIRES do not try to read beyond the last read (i.e. check that
  // numReadsSkip + numReadsKeep <= numReads).
  virtual ReadFile* subFile(long numReadsSkip, long numReadsKeep) = 0;
  // by default, both paired ends will be returned individually
  virtual void open() = 0;
  virtual void open(bool getBothEnds) = 0;
  // these perform the same function as open but specify a
  // ScoredSeq derived class that will be returned by getRead
  virtual void openWithCarriers(int numSets) = 0;
  virtual void openWithCarriers(int numSets, bool getBothEnds) = 0;
  virtual void openNormalized() = 0;
  virtual void openNormalized(bool getBothEnds) = 0;
  // WithCarriers encapsulates Norm
  virtual void openNormWithCarriers(int numSets) = 0;
  virtual void openNormWithCarriers(int numSets, bool getBothEnds) = 0;

  virtual bool hasRead() = 0;
  virtual ScoredSeq* getRead() = 0;
  virtual void skipRead() = 0;
  virtual void skipReads(long numToSkip) = 0;
  virtual void close() = 0;
  virtual long numReads() = 0;

  // expected distance between outside edges of paired ends.
  // treat as a maximum.  if reads are not paired-end, default
  // return value is zero (for paried-end data, zero should be
  // treated as infinity).
  virtual long ampliconSize() = 0;

  // old; get rid of this
  //virtual void getReads(vector<ScoredSeq*>* readCarrier) = 0;
  enum BufferThreadedness{ bufferNotThreaded=0, bufferThreaded=1 };


  // basic string interpretation
  class StringInterpreter{
  public:
    StringInterpreter(char* rawSeq, long rawLen, long maxLen, set<char>* okToSkip, bool revComp);
    ~StringInterpreter();
    // does NOT get deleted, just created (values are publicly accessible)
    char* _seqString;
    long _seqLen;
  };


 protected:
  friend class ReadFileCommunicator;
  friend class OutputFile;
  virtual void bufferQueue(Read * r) = 0;
  virtual void bufferUnqueue(Read * r) = 0;
  virtual void bufferFill() = 0; // default = not threaded
  virtual void bufferFill(BufferThreadedness threadedness) = 0;

 private:
  virtual void updateBuffer() = 0;
  virtual void emptyBuffer() = 0;
   

};

#endif
