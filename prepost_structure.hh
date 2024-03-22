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
// prepost_structure.hh
// Randal A. Koene, 20041118
//
// Classes for presynaptic and postsynaptic structure objects.

#ifndef __PREPOST_STRUCTURE_HH
#define __PREPOST_STRUCTURE_HH

// Required class headers

#include "fibre_structure.hh"
#include "spatial.hh"
#include "Network_Generated_Statistics.hh"
#include "state_storable.hh"

// Class pointer forward declarations (to avoid recursive inclusion)
#ifndef __NEURON_HH
class neuron;
#endif

class presynaptic_structure: public fibre_structure {
  // axonal morphology
public:
  presynaptic_structure(neuron & _n, Segment & initseg, spatial & acoords, int _typeid = -1);
  virtual terminal_segment_ptr * array_of_terminal_segments();
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(network_statistics_base & nsb);
#endif
};

class postsynaptic_structure: public fibre_structure {
  // dendritic morphology
public:
  postsynaptic_structure(neuron & _n, Segment & initseg, spatial & acoords, int _typeid = -1);
  virtual terminal_segment_ptr * array_of_terminal_segments();
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(network_statistics_base & nsb);
#endif
};

extern natural_schema_parent_set axons_most_specific_natural_set[];
extern natural_schema_parent_set dendrites_most_specific_natural_set[];

#endif
