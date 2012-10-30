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
 * allows an output file to be provided without needing to write
 * to an actual file
 */


#ifndef OUTPUTFILENULL_H
#define OUTPUTFILENULL_H

#include <set>
#include <fstream>
#include <string>
#include "ScoredSeq.h"
#include "OutputFile.h"
using namespace::std;


class OutputFileNull : public OutputFile {
 public:
  OutputFileNull();
  ~OutputFileNull();

  /* see parent class for doc */
  void open();
  void open(string addendum);
  bool isOpen();
  void writeContigs(ScoredSeq** contigs);
  void writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchangedContigs);
  void close();
  string getName();
  string getName(string addendum);

 private:
  bool _isOpen;

};

#endif
