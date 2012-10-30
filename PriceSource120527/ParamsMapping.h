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

/* packages paramaters that will be used by ExtendCycle when
 * mapping reads to contigs.
 */


#ifndef PARAMSMAPPING_H
#define PARAMSMAPPING_H


class ParamsMapping {

 public:
  static ParamsMapping* getDefault();

  ParamsMapping();
  ParamsMapping(long stepSize, float fractId);
  ParamsMapping(long stepSize, long stepSizeAddendum, float fractId);
  ParamsMapping(ParamsMapping * pmap);

  // get the values
  long stepSize();
  long stepSizeAddendum();
  float fractId();

 private:

  long _stepSize;
  long _stepSizeAddendum;
  float _fractId;

};

#endif
