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

#ifndef BURROWSWHEELERUTILITIES_CPP
#define BURROWSWHEELERUTILITIES_CPP

#include "burrowsWheelerUtilities.h"

#include <stdexcept>
#include <string>
#include <map>
#include <stdio.h>
#include <iostream>
#include "AssemblyException.h"
using namespace::std;



/* These first functions are for sorting suffixes in constant time.  They were lifted from
 * the text of the cited paper. They assume that the input text is formatted as an array of
 * longs.  The original had them as arrays of ints, so I changed them to longs because I am
 * going to have very large texts to index.
 * NOTE: I modified the functions below so that K is the size of the alphabet rather than the max char
 */

// SOURCE FROM: "simple linear work suffix array construction", karkkainen and sanders, ICALP, 203, p 943-955.
inline bool burrowsWheelerUtilities::leq(long a1, long a2, long b1, long b2) { // lexic. order for pairs
  return(a1 < b1 || a1 == b1 && a2 <= b2);
}                                             // and triples
inline bool burrowsWheelerUtilities::leq(long a1, long a2, long a3, long b1, long b2, long b3) {
  return(a1 < b1 || a1 == b1 && leq(a2,a3, b2,b3));
}

// SOURCE FROM: "simple linear work suffix array construction", karkkainen and sanders, ICALP, 203, p 943-955.
// stably sort a[0..n-1] to b[0..n-1] with keys in [0..K) from r
inline void burrowsWheelerUtilities::radixPass(long* a, long* b, long* r, long n, long K){ 
  // count occurrences
  long* c = new long[K];                          // counter array
  for (long i = 0;  i < K;  i++) c[i] = 0;         // reset counters
  for (long i = 0;  i < n;  i++) c[r[a[i]]]++;    // count occurences
  for (long i = 0, sum = 0;  i < K;  i++) { // exclusive prefix sums
    long t = c[i];  c[i] = sum;  sum += t;
  }
  for (long i = 0;  i < n;  i++) b[c[r[a[i]]]++] = a[i];      // sort
  delete [] c;
}


// just a translator so that suffixArray's input will be intuitive to use, but
// the function "suffixArray" can keep the naming conventions from "simple linear
// work suffix array construction", karkkainen and sanders, ICALP, 203, p 943-955.
// MODIFIES: sortedSuffixes by adding the output
// REQUIRES: text and sortedSuffixes are at least as long as textLen
// REQUIRES: text has all characters less than alphaSize and >=0
void burrowsWheelerUtilities::sortSuffixes(long* text, long* sortedSuffixes, long textLen, long alphaSize) {
  burrowsWheelerUtilities::suffixArray(text, sortedSuffixes, textLen, alphaSize);
}

// this allows a series of checks and safety procedures to be performed 
// initially and not repeated during recursion of the helper function
void burrowsWheelerUtilities::suffixArray(long* s, long* SA, long saLen, long K) {
  long sLen = saLen + 3; // three zeros will be added to the end
  long* newS = new long[sLen]; // create a new array that will have the appended zeros
  for (long n = 0; n < sLen; ++n){
    if (n < saLen){
      if (s[n] >= K){ throw AssemblyException::ArgError("in suffixArray, letter num too high for alphabet"); }
      newS[n] = s[n];
    }
    else { newS[n] = 0; }
  }
  suffixArrayHelper(newS, SA, saLen, K);
  delete [] newS;
}


// SOURCE FROM: "simple linear work suffix array construction", karkkainen and sanders, ICALP, 203, p 943-955.
// find the suffix array SA of s[0..n-1] in {1..K}^n
// require s[n]=s[n+1]=s[n+2]=0, n>=2
void burrowsWheelerUtilities::suffixArrayHelper(long* s, long* SA, long n, long K) {
  long n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2; 
  long* s12  = new long[n02 + 3];  s12[n02]= s12[n02+1]= s12[n02+2]=0; 
  long* SA12 = new long[n02 + 3]; SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;
  long* s0   = new long[n0];
  long* SA0  = new long[n0];
 
  // generate positions of mod 1 and mod  2 suffixes
  // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
  for (long i=0, j=0;  i < n+(n0-n1);  i++) if (i%3 != 0) s12[j++] = i;

  // lsb radix sort the mod 1 and mod 2 triples
  radixPass(s12 , SA12, s+2, n02, K);
  radixPass(SA12, s12 , s+1, n02, K);  
  radixPass(s12 , SA12, s  , n02, K);

  // find lexicographic names of triples
  long name = 0, c0 = -1, c1 = -1, c2 = -1;
  for (long i = 0;  i < n02;  i++) {
    if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) { 
      name++;  c0 = s[SA12[i]];  c1 = s[SA12[i]+1];  c2 = s[SA12[i]+2];
    }
    if (SA12[i] % 3 == 1) { s12[SA12[i]/3]      = name; } // left half
    else                  { s12[SA12[i]/3 + n0] = name; } // right half
  }

  // recurse if names are not yet unique
  if (name < n02) {
    suffixArrayHelper(s12, SA12, n02, name+1);
    // store unique names in s12 using the suffix array 
    for (long i = 0;  i < n02;  i++) s12[SA12[i]] = i + 1;
  } else // generate the suffix array of s12 directly
    for (long i = 0;  i < n02;  i++) SA12[s12[i] - 1] = i; 

  // stably sort the mod 0 suffixes from SA12 by their first character
  for (long i=0, j=0;  i < n02;  i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
  radixPass(s0, SA0, s, n0, K);

  // merge sorted SA0 suffixes and sorted SA12 suffixes
  for (long p=0,  t=n0-n1,  k=0;  k < n;  k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
    long i = GetI(); // pos of current offset 12 suffix
    long j = SA0[p]; // pos of current offset 0  suffix
    if (SA12[t] < n0 ? 
        leq(s[i],       s12[SA12[t] + n0], s[j],       s12[j/3]) :
        leq(s[i],s[i+1],s12[SA12[t]-n0+1], s[j],s[j+1],s12[j/3+n0]))
      { // suffix from SA12 is smaller
	SA[k] = i;  t++;
	if (t == n02) { // done --- only SA0 suffixes left
	  for (k++;  p < n0;  p++, k++) SA[k] = SA0[p];
	}
      } else { 
      SA[k] = j;  p++; 
      if (p == n0)  { // done --- only SA12 suffixes left
        for (k++;  t < n02;  t++, k++) SA[k] = GetI(); 
      }
    }  
  } 
  delete [] s12; delete [] SA12; delete [] SA0; delete [] s0; 
}



/* This next section is for converting texts to the number-array format used by the
 * functions above.
 */

// returns the max char number
long burrowsWheelerUtilities::getAlphabetSize(){ return 7; }

/* REQUIRES: longSeq must have room for all of the elements in seq after they are
 * offset by the specified amount
 * FYI, the encoding: '#'=0, '&'=1, 'A'=2, 'C'=3, 'G'=4, 'N'=5, 'T'=6
 */
void burrowsWheelerUtilities::makeNumString(char* seq, long seqLen, long offset, long* longSeq){
  for (long n = 0; n < seqLen; ++n){
    switch( seq[n] ) {
    case '#': longSeq[offset] = 0; break;
    case '&': longSeq[offset] = 1; break;
    case 'A': longSeq[offset] = 2; break;
    case 'C': longSeq[offset] = 3; break;
    case 'G': longSeq[offset] = 4; break;
    case 'N': longSeq[offset] = 5; break;
    case 'T': longSeq[offset] = 6; break;
    default:
      char errorMessage[100];
      char part1[] = "invalid nucleotide:";
      sprintf(errorMessage, "%s %c", part1, seq[n]);
      throw AssemblyException::ArgError(errorMessage);
    }
    ++offset;
  }
}



/* This next section is functions for doing the burrows-wheeler transform and search.
 */

// MODIFIES: columnL according to burrows-wheeler (assumed empty when input)
void burrowsWheelerUtilities::bwTransform(long* text, long* sortedSuffixes, long textLen, long* columnL){
  for (long n = 0; n < textLen; ++n){
    // note: suffixStart = sortedSuffixes[n];
    if (sortedSuffixes[n] == 0){ columnL[n] = text[textLen-1]; }
    else { columnL[n] = text[sortedSuffixes[n]-1]; }
  }
}


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
long* burrowsWheelerUtilities::fillBwtHelperArrays(long* columnL, long textLen, long alphaSize, long* occCounts){
  long alphaSizeP1 = alphaSize + 1;
  long* tableC = new long[alphaSizeP1];
  for (long c = 0; c < alphaSizeP1; ++c){ tableC[c] = 0; }
  for (long n = 0; n < textLen; ++n){
    tableC[columnL[n]+1]++;
    occCounts[n] = tableC[columnL[n]+1];
  }
  for (long c = 1; c < alphaSizeP1; ++c){ tableC[c] += tableC[c-1]; }
  if (tableC[alphaSize] != textLen){
    throw AssemblyException::LogicError("BWT: the sum of counts in tableC does not equal the text length.");
  }
  return tableC;
}


// this is the algorithm from fig 2 of:
// Ferrangina + Manzini, "An experimental study of an opportunistic index" in
// Proceedings of the twelfth annual ACM-SIAM symposium on discrete algorithms
// 2001: 269-278.
// the algorithm was modified for zero-indexing
// REQUIRES: qSize > 0
// RETURNS: first index in return value is the number of additional indexes; those
//    indexes are filled with start positions.
long* burrowsWheelerUtilities::bwtCountAndFind(long* query, long qStart, long qSize,
					       long* columnL, long* sortedSuffixes, long textSize,
					       long* cCounts, long* occCounts){

  long i = qSize + qStart - 2;
  long c = query[i+1];
  long sp = cCounts[c];
  long ep = cCounts[c+1] - 1;
  while ( sp <= ep and i >= qStart ){
    c = query[i];
    sp = cCounts[c] + findOccCount(sp-1, c, columnL, textSize, occCounts);
    ep = cCounts[c] + findOccCount(ep, c, columnL, textSize, occCounts) - 1;
    i--;
  }

  long count = ep - sp + 1;
  if (count < 0){ count = 0; }
  long* indexes = new long[count+1]; // first index is the count; remaining indexes are the indexes of the matches
  indexes[0] = count;
  for (long n = 0; n < count; ++n){
    indexes[n+1] = sortedSuffixes[sp+n];
  }
  return indexes;
}

inline long burrowsWheelerUtilities::findOccCount(long pos, long c, long* text, long textSize, long* occCounts){
  if (pos >= textSize){ pos = textSize - 1; }
  while (pos >= 0 and text[pos] != c){ --pos; }
  if (pos == -1){ return 0; }
  else { return occCounts[pos]; }
}


#endif
