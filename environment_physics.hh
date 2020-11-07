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
// environment_physics.hh
// Randal A. Koene, 20070710

/* Documentation:
   This implements environmental effects that represent physical influences
   on network development, such as environmental constraint within a sphere
   (the physical boundaries of the brain) and consequent effects on
   neurite direction and neurite elongation. Other models that can be
   included here may represent physical barriers, channels and obstacles.

   At this time, these classes are applied equally to all neurons in a
   network and therefore do not participate in the NETMORPH interface
   protocol for model and parameter selection according to the smallest
   matching set.
 */

#ifndef __ENVIRONMENT_PHYSICS_HH
#define __ENVIRONMENT_PHYSICS_HH

#include "templates.hh"
//#include "fibre_structure.hh"
#include "spatial.hh"
#include "Command_Line_Parameters.hh"
//#include "global.hh"

enum physical_boundary_types { spherical_pb, point_attractor_pb, NUM_pb };

extern const char physical_boundary_label[NUM_pb][20];

#ifdef DEBUG_SPHERICAL_BOUNDARY
extern unsigned long dsb_zerolength;
extern unsigned long dsb_in_after_adjust;
extern unsigned long dsb_safe;
extern unsigned long dsb_from_beyond;
extern unsigned long dsb_no_move_to_border;
extern unsigned long dsb_adjusted;
#endif

class physical_boundary: public PLLHandle<physical_boundary>, public CLP_Modifiable {
  // All the boundaries are considered solid, so each chained
  // boundary effect must be applied in full the result of the
  // preceding boundary effect.
  // Also, all boundaries are considered independent actors, so
  // that chaining is really just a convenient way to create an
  // array. Consequently, all labels are at the top level, all
  // boundaries must have a label, and setting of a boundary
  // begins with providing a new label.
protected:
  String thislabel;
  //  physical_boundary * contributing;
  //double contributingweight;
public:
  physical_boundary(String & label): thislabel(label) {} //, contributing(pbcontrib), contributingweight(pbweight) {}
  physical_boundary(String & label, Command_Line_Parameters & clp): thislabel(label) {} //, contributing(NULL), contributingweight(0.0) {}
  virtual ~physical_boundary() {}
  virtual void direction_boundary_effect(spatial & P, spatial & angularoffset) = 0;
};

class spherical_boundary: public physical_boundary {
protected:
  spatial center;
  double radius;
  double maxrange; // effects are applied between radius-maxrange and radius
  double c_repulse; // coefficient of the repulsive force that corrects neurite growth for the density in the boundary field
  bool direction_adjust(spatial & P, spatial & angularoffset, int passnum);
public:
  spherical_boundary(String & label): physical_boundary(label), radius(600.0), maxrange(50.0), c_repulse(0.5) {}
  spherical_boundary(String & label, Command_Line_Parameters & clp): physical_boundary(label,clp), radius(600.0), maxrange(50.0), c_repulse(0.5) { parse_CLP(clp); }
#ifdef DEBUG_SPHERICAL_BOUNDARY
  virtual ~spherical_boundary();
#endif
  virtual void direction_boundary_effect(spatial & P, spatial & angularoffset);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class point_attractor: public physical_boundary {
protected:
  spatial attractor;
  double maxrange; // effects are applied up to a radius of maxrange around the attractor
  double c_attract; // coefficient of the attractive force that affects neurite growth
  bool direction_adjust(spatial & P, spatial & angularoffset); //, int passnum);
public:
  point_attractor(String & label): physical_boundary(label), maxrange(100.0), c_attract(0.5) {}
  point_attractor(String & label, Command_Line_Parameters & clp): physical_boundary(label,clp), maxrange(100.0), c_attract(0.5) { parse_CLP(clp); }
  virtual void direction_boundary_effect(spatial & P, spatial & angularoffset);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class environment_physics: public PLLRoot<physical_boundary>, public CLP_Modifiable {
public:
  environment_physics() {}
  int select_physical_boundary(String prefixlabel, Command_Line_Parameters & clp);
  void direction_boundary_effect(spatial & P, spatial & angularoffset) {
    PLL_LOOP_FORWARD(physical_boundary,head(),1) e->direction_boundary_effect(P,angularoffset);
  }
  virtual void parse_CLP(Command_Line_Parameters & clp) { select_physical_boundary("",clp); }
  virtual String report_parameters() {
    String res;
    PLL_LOOP_FORWARD(physical_boundary,head(),1) res += e->report_parameters();
    return res;
  }
};

extern environment_physics pb;

#endif
