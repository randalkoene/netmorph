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
// Connection_Statistics.hh
// Randal A. Koene, 20041118
//
// Classes for the specification of network connection statistics.

#ifndef __CONNECTION_STATISTICS_HH
#define __CONNECTION_STATISTICS_HH

#include "BigString.hh"
#include "templates.hh"
#include "neuron.hh"

class Connection_Statistics: public PLLHandle<Connection_Statistics> {
protected:
  neuron_type sourcetype;
  neuron_type targettype;
  int minconn, maxconn;
  double range;
public:
  Connection_Statistics(neuron_type stype, neuron_type ttype, int minc, int maxc, double r): sourcetype(stype), targettype(ttype), minconn(minc), maxconn(maxc), range(r) {}
  neuron_type SourceType() { return sourcetype; }
  neuron_type TargetType() { return targettype; }
  int MinConn() { return minconn; }
  int MaxConn() { return maxconn; }
  double Range() { return range; }
  String str();
  String all_str(char sep = '\n');
};

class Connection_Statistics_Root: public PLLRoot<Connection_Statistics> {
public:
  Connection_Statistics_Root() {}
  Connection_Statistics * get_Connection_Statistics(neuron_type stype, neuron_type ttype);
};

#endif
