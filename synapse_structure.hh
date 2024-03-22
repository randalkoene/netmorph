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
// synapse_structure.hh
// Randal A. Koene, 20050215
//
// Classes that define spatial structure detail of a synapse.

#ifndef __SYNAPSE_STRUCTURE_HH
#define __SYNAPSE_STRCUTURE_HH

// Required class headers

#include "spatial.hh"
#include "state_storable.hh"

// Class pointer forward declarations (to avoid recursive inclusion)
#ifndef __FIG_OBJECT_HH
class Fig_Object;
#endif
#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
#endif

class synapse_structure: public Segment, public state_storable {
  // The Segment is initialized without specification for the axonsegment or
  // dendritesegment spatial information.
  // The Segment is properly defined by calling dist3D_Segment_to_Segment with
  // the axon and dendrite segments, and this object as the third parameter,
  // see synapse_formation_model::evaluate_possible_connection().
protected:
  fibre_segment * axonsegment, * dendritesegment; // *** perhaps remove if unused
  int usedbyNsynapses; // counts how many synapses use this structure definition
public:
  synapse_structure(fibre_segment & a, fibre_segment & d): axonsegment(&a), dendritesegment(&d), usedbyNsynapses(0) {}
  void move_add(spatial & addvec);
  void fanin_rot();
  Fig_Object * net_Fig();
  void require() { usedbyNsynapses++; }
  void release() { usedbyNsynapses--; if (usedbyNsynapses<=0) delete this; }
  fibre_segment * AxonSegment() { return axonsegment; }
  fibre_segment * DendriteSegment() { return dendritesegment; }
};

#endif
