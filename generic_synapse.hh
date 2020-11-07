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
// generic_synapse.hh (nibr)
// Randal A. Koene, 20041111

#ifndef __GENERIC_SYNAPSE_HH
#define __GENERIC_SYNAPSE_HH

// Required class headers
#include "synapse.hh"

// Class pointer forward declarations (to avoid recursive inclusion)
#ifndef __CONNECTION_HH
class connection;
#endif

// functional requirements met:
//   synapse formed by presynaptic and postsynaptic structures
// physiology implemented:
//   presynaptic and postsynaptic structures have specific arborization

class generic_synapse: public synapse {
protected:
  /// presynaptic_node prenode; // terminal
  /// postsynaptic_node postnode; // receptive zone
public:
  generic_synapse(connection & conn); // (see generic_synapse.cc)
  //~generic_synapse();
};

#endif

