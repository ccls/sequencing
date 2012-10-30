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

#ifndef BURROWSWHEELERUTILITIES_H
#define BURROWSWHEELERUTILITIES_H

#include <stdexcept>
#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <set>
#include <map>
#include "AssemblyException.h"
using namespace::std;

namespace burrowsWheelerUtilities {

  // function defs

  // for sorting suffix arrays
  inline bool leq(long a1, long a2, long b1, long b2);
  inline bool leq(long a1, long a2, long a3,  long b1, long b2, long b3);
  inline void radixPass(long* a, long* b, long* r, long n, long K);
  void suffixArrayHelper(long* s, long* SA, long n, long K);
  void suffixArray(long* s, long* SA, long n, long K);

  // sorts the suffixes of the text lexicographically
  // MODIFIES: sortedSuffixes by adding the output
  // REQUIRES: text and sortedSuffixes are at least as long as textLen
  // REQUIRES: text has all characters less than alphaSize and >=0
  void sortSuffixes(long* text, long* sortedSuffixes, long textLen, long alphaSize);


  // for converting contig seqs to number arrays
  long getAlphabetSize(); // size of the alphabet used by makeNumString

  // translates a DNA sequence to numbers.  valid letters are 'A', 'T', 'C', 'G',
  // '&' (used to separate distinct seqs), '#' (marks the end of a text)
  // MODIFIES: longSeq by adding the numbers from the seq
  // REQUIRES: longSeq must have room for all of the elements in seq after they are
  // offset by the specified amount
  void makeNumString(char* seq, long seqLen, long offset, long* longSeq);


  // set up for bwt searches
  // MODIFIES: columnL according to burrows-wheeler (assumed empty when input)
  void bwTransform(long* text, long* sortedSuffixes, long textLen, long* columnL);
  // MODIFIES: occCounts
  // occCounts is the number of counts for the letter at the given position.
  //   for other letters, one will need to backtrack until the letter is found
  //   or the beginning of the text is reached.  this will not be a problem due
  //   to the finite size of the alphabet and the anticipated sizes of the text
  //   relative to it.
  // REQUIRES: column L is properly encoded with num representations for
  //   letters as in tableC
  // REQUIRES: columnL len == textLen
  // REQUIRES: occCounts len >= textLen (should be ==, but > is ok)
  // RETURNS: tableC has as many indexes as letters from the alphabet PLUS ONE
  long* fillBwtHelperArrays(long* columnL, long textLen, long alphaSize, long* occCounts);


  // perform bwt searches
  long findOccCount(long pos, long c, long* text, long textSize, long* occCounts);
  // this is the algorithm from fig 2 of:
  // Ferrangina + Manzini, "An experimental study of an opportunistic index" in
  // Proceedings of the twelfth annual ACM-SIAM symposium on discrete algorithms
  // 2001: 269-278.
  // the algorithm was modified for zero-indexing
  // REQUIRES: qSize > 0
  // RETURNS: first index in return value is the number of additional indexes; those
  //    indexes are filled with start positions.
  //long* bwtCountAndFind(long* query, long qSize, long* columnL, long* sortedSuffixes, long textSize, long* cCounts, long* occCounts);
  long* bwtCountAndFind(long* query, long qStart, long qSize, long* columnL, long* sortedSuffixes, long textSize, long* cCounts, long* occCounts);

}


#endif
