/*
  © Copyright 2008 Randal A. Koene <randalk@netmorph.org>
  
  With design assistance from J. van Pelt & A. van Ooyen, and support
  from the Netherlands Organization for Scientific Research (NWO)
  Program Computational Life Sciences grant CLS2003 (635.100.005) and
  from the EC Marie Curie Research and Training Network (RTN)
  NEURoVERS-it 019247.

  This file is part of NETMORPH.

  NETMORPH is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  NETMORPH is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.
*/
// state_storable.hh
// Randal A. Koene, 20050827
//
// Class definitions that obligate state storable objects to provide specific
// functions that enable the storage and retrieval of simulation state.

#ifndef __STATE_STORABLE_HH
#define __STATE_STORABLE_HH

#include "BigString.hh"

/* Documentation:
   The state_storable class is inherited by all objects for which state needs
   to be stored and retrieved so that a simulation may be processed further at
   another time.
   The specific functions that are defined in the derived classes should
   call on a derivative object of state_file_format to store and retrieve
   such elemental things as numerical data, list syntax, etc.
 */

class state_storable {
protected:
  state_storable() {}
  virtual ~state_storable() {}
public:
  virtual String state_output();
  virtual void state_input(String s);
};

#endif
