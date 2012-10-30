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


#ifndef OUTPUTFILEPRICEQ_H
#define OUTPUTFILEPRICEQ_H

#include <set>
#include <fstream>
#include <string>
#include "ScoredSeq.h"
#include "OutputFile.h"
using namespace::std;


class OutputFilePriceq : public OutputFile {
 public:
  OutputFilePriceq();
  OutputFilePriceq(string filename);
  ~OutputFilePriceq();

  /* see parent class for doc */
  void open();
  void open(string addendum);
  bool isOpen();
  void writeContigs(ScoredSeq** contigs);
  void writeContigs(ScoredSeq** contigs, set<ScoredSeq*>* unchangedContigs);
  void close();
  string getName();
  string getName(string addendum);

  // these are public so that I can test them and also so that I
  // can have the two methods in the same place and still be accessible
  // to other classes that may use them (like the input file class).
  static char scoreToChar(float score);
  static float charToScore(char scoreChar);

 private:
  string _filenameBegin;
  string _filenameEnd;
  int _charPerLine;

  bool _isOpen;
  ofstream _outputFileObj;
  long _contigCounter;
};

#endif
