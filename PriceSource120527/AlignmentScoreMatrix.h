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
 */

#ifndef ALIGNMENTSCOREMATRIX_H
#define ALIGNMENTSCOREMATRIX_H


class AlignmentScoreMatrix {
 public:
  static AlignmentScoreMatrix* getDefault();

  AlignmentScoreMatrix();
  AlignmentScoreMatrix(long match, long mismatch, long newGap, long extendGap);
  AlignmentScoreMatrix(AlignmentScoreMatrix * alsm);
  long match(char nucA, char nucB);
  bool isMatch(char nucA, char nucB);

  // the match score should have a positive value (but doesn't have to)
  long _match;
  // these penalties should have negative values (but don't have to)
  long _mismatch;
  long _newGap;
  long _extendGap;
 private:

};


#endif
