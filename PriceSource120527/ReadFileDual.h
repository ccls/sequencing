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
This is for paired-end reads where the 1st and 2nd reads are found in
separate files.  The two files must contain the same number of reads.
 */


#ifndef READFILEDUAL_H
#define READFILEDUAL_H

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


class ReadFileDual : public ReadFile {
 public:

  ReadFileDual();
  // copies the files, so these can be deleted if they are not needed independently
  ReadFileDual(ReadFile* fileA, ReadFile* fileB, long ampliconSize);
  ~ReadFileDual();

  /* see parent class ReadFile */
  ReadFileDual* copy();
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
  long numReads();
  long ampliconSize();

  // old; get rid of this
  //void getReads(vector<ScoredSeq*>* readCarrier);


 protected:
  friend class ReadFileCommunicator;
  void bufferQueue(Read * r);
  void bufferUnqueue(Read * r);
  void bufferFill();
  void bufferFill(BufferThreadedness threadedness);

 private:
  // used privately for creating the split files
  ReadFileDual(ReadFileDual* parent, ReadFile* rfA, ReadFile* rfB);

  bool _getBothEnds;

  void updateBuffer();
  void emptyBuffer();

  bool _openWithCarriers;
  int _numCarrierSets;
  bool _openNormalized;
  bool _isOpen;

  long _ampliconSize;

  void getReadPairHelper();
  queue<ScoredSeq*> _getReadsQueue;


  ReadFile* _readFileA;
  ReadFile* _readFileB;

};

#endif
