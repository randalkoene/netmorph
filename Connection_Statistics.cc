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
// Connection_Statistics.cc
// Randal A. Koene, 20041118

#include "Connection_Statistics.hh"

String Connection_Statistics::str() {
  String s(neuron_type_name[sourcetype]);
  s += " to ";
  s += neuron_type_name[targettype];
  s += ": #connections [" + String((long) minconn) + ',' + String((long) maxconn) + "], range " + String(range,"%.3f") + " micrometers";
  return s;
}

String Connection_Statistics::all_str(char sep) {
  String s = str();
  PLL_LOOP_FORWARD(Connection_Statistics,Next(),1) s += sep + e->str();
  return s;
}

Connection_Statistics * Connection_Statistics_Root::get_Connection_Statistics(neuron_type stype, neuron_type ttype) {
  // [***NOTE] This is slow, an alternative is to set up a look-up table as
  // elements are linked to the root and to then simply return pointers from
  // the look-up table in this function.
  PLL_LOOP_FORWARD(Connection_Statistics,head(),1) if ((e->SourceType()==stype) && (e->TargetType()==ttype)) return e;
  return NULL;
}
