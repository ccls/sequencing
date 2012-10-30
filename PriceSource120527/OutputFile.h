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
 * This is an abstract class for an output file of contigs.  Having it as
 * an interface will allow the client to select from multiple formats to
 * feed to the Assembler instance, which can treat all different output
 * formats uniformly.
 */


#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#include <set>
#include "ScoredSeq.h"
using namespace::std;


class OutputFile {
 public:
  static OutputFile* makeOutputFile(string filename);
  // the output array ends with a NULL
  static ScoredSeq** makeSortedContigArray(set<ScoredSeq*>* contigs);

  virtual ~OutputFile();

  // these open the filename specified at creation; contents
  // of an existing file will be deleted.  the addendum allows
  // multiple output files with the same core name to be created.
  virtual void open() = 0;
  virtual void open(string addendum) = 0;
  virtual bool isOpen() = 0;
  // seqs will appear in the output file in the order in which they are delivered in the array;
  // the array must end with a NULL
  virtual void writeContigs(ScoredSeq** contigs) = 0;
  virtual void writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchangedContigs) = 0;
  virtual void close() = 0;

  // find out the name, with or without addendum
  virtual string getName() = 0;
  virtual string getName(string addendum) = 0;

  // used by other classes to make sure that the output file is
  // in a real directory - does so by creating a temporary version
  // of the file that will then be deleted (if the file does not
  // already exist).
  static void ensureWritability(const char* filename);

 private:
  struct SortSsByLength {
    bool operator() (ScoredSeq* ssnA, ScoredSeq* ssnB);
  };

};

#endif
