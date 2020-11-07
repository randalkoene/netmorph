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
// nibr.hh
// Randal A. Koene, 20040830

#ifndef __NIBR_HH
#define __NIBR_HH

#include "Command_Line_Parameters.hh"
#include "network.hh"
#include "spatial.hh"

// definitions

#ifdef VECTOR3D
#define RCFILE "./.nibrrc"
#define LASTCLPFILE ".nibr.clp"
#endif
#ifdef VECTOR2D
#define RCFILE "./.nibr2Drc"
#define LASTCLPFILE ".nibr2D.clp"
#endif

// function declarations

electrodearray * make_electrodes_hexagon(double centerx = 0.0, double centery = 0.0);

void running_instance_preop();
void running_instance_postop();
void single_instance_restriction();

void report_form_aware(Command_Line_Parameters & clp, CLP_Modifiable * clpm);
void report_form_aware(Command_Line_Parameters & clp, String s);

#endif
