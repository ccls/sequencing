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

/* this is an abstract class; multiple implementations of it will
   be created for return by a factory*/

#ifndef READFILEINDEX_H
#define READFILEINDEX_H
#include <iostream>
using namespace std;


class ReadFileIndex {
 public:
  virtual ~ReadFileIndex();
  /* NEED METHOD SPECS */
  virtual long readStart() = 0;
  virtual long readSize() = 0;
  virtual long scoreStart() = 0;
  virtual long scoreSize() = 0;
  virtual long linkStart() = 0;
  virtual long linkSize() = 0;
  virtual int blockNum() = 0;
};





#endif
