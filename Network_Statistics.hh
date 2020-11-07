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
// Network_Statistics.hh
// Randal A. Koene, 20041118
//
// Classes that provide statistics for network development.

#ifndef __NETWORK_STATISTICS_HH
#define __NETWORK_STATISTICS_HH

#include "templates.hh"
#include "BigString.hh"
#include "neuron.hh"

class Network_Statistics: public PLLHandle<Network_Statistics> {
  // provides statistics for network generation
protected:
  neuron_type ntype; // network statistics for this type
  double ratio; // ratio of cells in network that are of type ntype
  int populationsize; // number of cells in network that are of type ntype
public:
  Network_Statistics(neuron_type nt, double r): ntype(nt), ratio(r), populationsize(-1) {}
  Network_Statistics(neuron_type nt, int n): ntype(nt), ratio(0.0), populationsize(n) {}
  neuron_type TypeID() { return ntype; }
  double Ratio() { return ratio; }
  int PopulationSize() { return populationsize; }
  String str(); // (see nibr.cc)
  String all_str(char sep = '\n'); // (see nibr.cc)
};

class Network_Statistics_Root: public PLLRoot<Network_Statistics> {
  bool exactpopulationsizes;
public:
  Network_Statistics_Root(bool epsizes = false): exactpopulationsizes(epsizes)  {}
  bool UseExactPopulationSizes() const { return exactpopulationsizes; }
};

#endif
