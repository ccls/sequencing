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
 * how quickly initial contigs will be added to the assembly process.
 * will delete the read file when it is deleted
 */


#ifndef PARAMETERIZEDINITIALFILE_H
#define PARAMETERIZEDINITIALFILE_H

#include "ReadFile.h"
#include "ScoredSeq.h"
#include <set>

class ParameterizedInitialFile {

 public:
  ParameterizedInitialFile();
  ParameterizedInitialFile(ReadFile* readFile);
  ParameterizedInitialFile(ReadFile* readFile, int numSteps, int cyclesPerStep);
  ParameterizedInitialFile(ReadFile* readFile, int numSteps, int cyclesPerStep, int initialCycle);
  ParameterizedInitialFile(long totalInitial, ReadFile* readFile, int numSteps, int cyclesPerStep);
  ParameterizedInitialFile(long totalInitial, ReadFile* readFile, int numSteps, int cyclesPerStep, int initialCycle);
  ~ParameterizedInitialFile();

  // the min number of cycles required for this file to deploy its
  // contigs over the provided number of steps (all steps complete)
  int totalCycles();
  void getContigs(set<ScoredSeq*>* initialContigs, int cycleNum);

 private:
  void makeInitDeploymentArray(int numSteps, int cyclesPerStep);
  void fillCycleSelectionArray(bool* selectionArray, int cycle);

  // keys are zero-indexed cycle nums
  long* _cycleToInputNum;
  long* _cycleToOutputNum;
  long* _cycleToPastInputCount;

  long _totalInput;
  long _totalOutput;

  int _initialCycle;
  int _totalCycles;
  ReadFile* _readFile;

};

#endif
