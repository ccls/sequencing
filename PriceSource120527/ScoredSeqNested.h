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

/* This subclass of ScoredSeq is for those ScoredSeq types that wrap other
 * ScoredSeq objects and those objects are/can be accessible directly to
 * the external environment.
 */

#ifndef SCOREDSEQNESTED_H
#define SCOREDSEQNESTED_H

#include "ScoredSeq.h"

class ScoredSeqNested : virtual public ScoredSeq {

 public:
  virtual ~ScoredSeqNested();
  virtual ScoredSeq * getNested() = 0;

};

#endif
