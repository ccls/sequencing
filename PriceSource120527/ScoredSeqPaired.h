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

/* This subclass of ScoredSeq is for those ScoredSeq types that have paired
 * ends and need these methods to maintain there paired ends when deleting.
 */

#ifndef SCOREDSEQPAIRED_H
#define SCOREDSEQPAIRED_H

#include "ScoredSeq.h"

class ScoredSeqPaired : virtual public ScoredSeq {

 public:
  virtual ~ScoredSeqPaired();
  /* checks if the paired-end seq is ok to delete (i.e. is not alive
   * somewhere).  also, the power is provided publicly to enliven a
   * seq.  The power to de-enliven a seq is protected by the individual
   * seqs themselves.
   */
  virtual void makeAlive() = 0;
  virtual bool isAlive() = 0;
};

#endif
