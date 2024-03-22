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
// synapse_structure.cc
// Randal A. Koene, 20050215

#include "synapse_structure.hh"
#include "Fig_Object.hh"
#include "Sampled_Output.hh"
#include "global.hh"

void synapse_structure::move_add(spatial & addvec) {
  P0 += addvec;
  P1 += addvec;
}

void synapse_structure::fanin_rot() {
  P0.convert_to_spherical();
  P0.set_Y(0.0);
  P0.convert_from_spherical();
  P1.convert_to_spherical();
  P1.set_Y(0.0);
  P1.convert_from_spherical();
}

Fig_Object * synapse_structure::net_Fig() {
  if ((!figattr_fibres_nobox) && ((!fig_in_zoom(P0)) || (!fig_in_zoom(P1)))) return NULL;
  double x1, y1, x2, y2;
  P0.plane_mapped(x1,y1);
  P1.plane_mapped(x2,y2);
  return new Fig_Line(0,1,_figvar_connection_color,7,_figvar_connection_depth,0,0.0,0,0,FIGSCALED(x1),FIGSCALED(y1),FIGSCALED(x2),FIGSCALED(y2));
}
