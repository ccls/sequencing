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
 * when paired-end read sets will be used for the assembly process.
 * will delete the read file when it is deleted
 */


#ifndef PARAMETERIZEDREADFILE_H
#define PARAMETERIZEDREADFILE_H

#include "ReadFile.h"
#include "ScoredSeq.h"
#include "ParamsMapping.h"
#include <set>

class ParameterizedReadFile : public ReadFile{

 public:
  ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm);
  ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle);
  ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle, int numCycles);
  ParameterizedReadFile(ReadFile* readFile, ParamsMapping* pm, int initialCycle, int cyclesOn, int cyclesSkipped);
  ~ParameterizedReadFile();

  bool useThisCycle(int cycleNum);
  ParamsMapping* getMapParams(); // returns a copy
  // returns the same as class function below, but prevents the need for recasting
  void splitFile(set<ParameterizedReadFile*>* putFilesHere, int numFiles);
  ParameterizedReadFile* copy();

  // inherited interface from ReadFile; these are forwarded to the carried file
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

 protected:
  void bufferQueue(Read * r);
  void bufferUnqueue(Read * r);
  void bufferFill();
  void bufferFill(BufferThreadedness threadedness);

 private:
  ParameterizedReadFile(ParameterizedReadFile* readFile, ReadFile* replacementFile);

  void updateBuffer();
  void emptyBuffer();

  int _initialCycle;

  // there are two use-limiting modes, either finite or cyclic.  each has
  // variables associated with it that will not be used if the other mode
  // is in play.  note that it is also possible for neither of these limiting
  // modes
  bool _isCyclic;
  int _cyclesOn;
  int _cyclesSkipped;

  bool _isFinite;
  int _totalCycles;

  ReadFile* _readFile;
  ParamsMapping* _pm;

};

#endif
