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
// Network_Statistics.cc
// Randal A. Koene, 20041118

#include "Network_Statistics.hh"

String Network_Statistics::str() {
  String s(neuron_type_name[ntype]);
  if (populationsize<0) s += " cell proportion = " + String(ratio,"%.3f");
  else s += " cells = " + String((long) populationsize);
  return s;
}

String Network_Statistics::all_str(char sep) {
  String s = str();
  PLL_LOOP_FORWARD(Network_Statistics,Next(),1) s += sep + e->str();
  return s;
}

