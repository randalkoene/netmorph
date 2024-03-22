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
// environment_physics.cc
// Randal A. Koene, 20070710

#include <math.h>
#include "environment_physics.hh"
#include "event.hh"
#include "neuron.hh"
#include "network.hh"
#include "spatial.hh"
#include "global.hh"

// constants

const char physical_boundary_label[NUM_pb][20] = {
  "spherical",
  "point_attractor"
};

#define MAX_PB_ADJUSTS 100

// variables

environment_physics pb;

// classes
  // Note that handling labels is implemented slightly differently for the chains created here
  // than in the standard chaining method of models inherited from schemas and parents.

void spherical_boundary::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_radius"))>=0) radius = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_maxrange"))>=0) maxrange = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_c_repulse"))>=0) {
    double new_c_repulse = atof(clp.ParValue(n));
    if ((new_c_repulse>=0.0) && (new_c_repulse<=1.0)) c_repulse = new_c_repulse;
    else warning("Warning: The repulsion coefficient of a spherical boundary ("+thislabel+".spherical_c_repulse) must lie in [0,1].\n");
  }
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_center"))>=0) {
    if (!eq) warning("Warning: In order to set "+thislabel+".spherical_center to a neuron index the Event_Queue must be initialized first.\n");
    else if (!eq->Net()) warning("Warning: In order to set "+thislabel+".spherical_center to a neuron index the network must be initialized first.\n");
    else {
      neuron * nptr = static_cast<PLLRoot<neuron> *>(eq->Net())->el(atoi(clp.ParValue(n)));
      if (nptr) center = nptr->Pos();
      else warning("Parameter ERROR, "+thislabel+".spherical_center: Neuron index "+clp.ParValue(n)+" does not exist in network.\n");
    }
  }
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_centerX"))>=0) center.set_X(atof(clp.ParValue(n)));
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_centerY"))>=0) center.set_Y(atof(clp.ParValue(n)));
#ifdef VECTOR3D
  if ((n=clp.Specifies_Parameter(thislabel+".spherical_centerZ"))>=0) center.set_Z(atof(clp.ParValue(n)));
#endif
  if (Root()) static_cast<environment_physics *>(Root())->select_physical_boundary(thislabel,clp); // check for more with thislabel as prefix
}

String spherical_boundary::report_parameters() {
  String res("Spherical Boundary: center=(");
  res += String(center.X(),"%.3f,");
  res += String(center.Y(),"%.3f");
#ifdef VECTOR3D
  res += String(center.Z(),",%.3f");
#endif
  res += ") radius=";
  res += String(radius,"%.3f");
  res += " maxrange=";
  res += String(maxrange,"%.3f");
  res += " c_repulse=";
  res += String(c_repulse,"%.3f");
  return res;
}

#ifdef DEBUG_SPHERICAL_BOUNDARY
unsigned long dsb_zerolength = 0;
unsigned long dsb_in_after_adjust = 0;
unsigned long dsb_safe = 0;
unsigned long dsb_from_beyond = 0;
unsigned long dsb_no_move_to_border = 0;
unsigned long dsb_adjusted = 0;
#define RETURN_DSB_ZEROLENGTH { dsb_zerolength++; return false; }
#define RETURN_DSB_IN_AFTER_ADJUST { dsb_in_after_adjust++; return false; }
#define RETURN_DSB_SAFE { dsb_safe++; return false; }
#define RETURN_DSB_FROM_BEYOND { dsb_from_beyond++; return false; }
#define RETURN_DSB_NO_MOVE_TO_BORDER { dsb_no_move_to_border++; return false; }
#define RETURN_DSB_ADJUSTED { dsb_adjusted++; return true; }
#else 
#define RETURN_DSB_ZEROLENGTH return false
#define RETURN_DSB_IN_AFTER_ADJUST return false
#define RETURN_DSB_SAFE return false
#define RETURN_DSB_FROM_BEYOND return false
#define RETURN_DSB_NO_MOVE_TO_BORDER return false
#define RETURN_DSB_ADJUSTED return true
#endif

#ifdef DEBUG_SPHERICAL_BOUNDARY
spherical_boundary::~spherical_boundary() {
  cout << "DEBUG_SPHERICAL_BOUNDARY, number of returns per situation:\n";
  cout << "  zero length terminal segment              : " << dsb_zerolength << '\n';
  cout << "  within sphere after adjustment            : " << dsb_in_after_adjust << '\n';
  cout << "  safely within sphere before maxrange field: " << dsb_safe << '\n';
  cout << "  originates and remains beyond sphere      : " << dsb_from_beyond << '\n';
  cout << "  grows toward center, not border           : " << dsb_no_move_to_border << '\n';
  cout << "  ADJUSTED                                  : " << dsb_adjusted << '\n';
}
#endif

bool spherical_boundary::direction_adjust(spatial & P, spatial & angularoffset, int passnum) {
  // For growth cones within the sphere, this function insures repulsion from
  // the spherical boundary.
  // Note the case when P lies outside the spherical boundary: When implemented
  // without a special test case for d_P < 0, this function would leave growth
  // cones already heading toward the sphere unmodified, while for those cases
  // where d_P - d_Pnew > 0, it sets up a force (repulse = -d_Pnew + d_p*c_repulse > 0),
  // which is aimed toward the center of the sphere, thereby turning those growth
  // cones around to head toward the sphere.
  // It may be more desirable to handle that as a separate intended use, that may
  // distriminate between randomly directed growth cones, ignoring or repulsing
  // those headed in a different direction, while attracting those headed in the
  // direction of the sphere. I will implement that possibility as:
  // stochastic_discriminating_sphere.
  // Returns false if no further tests are required and true if another pass
  // should be made to verify if Pnew is now within the sphere (note the special
  // test in case both P and Pnew are outside the sphere).
  if (angularoffset.X()==0.0) RETURN_DSB_ZEROLENGTH; // catch this special case faster than after computing d_approach
  spatial a(angularoffset); // local working copy
  a.convert_from_spherical();
  spatial Pnew(P);
  Pnew += a;
  double d_Pnew = radius - Pnew.distance(center);
  if ((passnum>0) && (d_Pnew>=0.0)) RETURN_DSB_IN_AFTER_ADJUST; // after adjustment, within sphere, no more adjustments needed
  if (d_Pnew > maxrange) RETURN_DSB_SAFE;
  double d_P = radius - P.distance(center); // distances from border
  if ((passnum>0) && (d_P<0.0)) RETURN_DSB_FROM_BEYOND; // after adjustment of growth cone originating outside sphere, no more than one adjustment
  double d_approach = d_P - d_Pnew; // if positive, then got closer
  if (d_approach <= 0.0) RETURN_DSB_NO_MOVE_TO_BORDER;
  // the following is applied to growth cones that have been heading toward the boundary, and are within the maxrange field to the boundary
  double repulse;
  if ((d_Pnew >= 0.0) && (maxrange>0.0)) {
    double fieldratio = d_Pnew / maxrange;
    repulse = d_approach*c_repulse*(1 - (fieldratio*fieldratio)); // repulsion to correct growth for density in the boundary field
  } else repulse = (-d_Pnew) + (d_P*c_repulse); // if passed beyond boundary, or if maxrange==0.0 ([***NOTE] What does maxrange<0 mean?)
  spatial repulsevec(center);
  repulsevec -= Pnew;
  repulsevec.normalize();
  repulsevec *= repulse;
  a += repulsevec;
  a.convert_to_spherical();
  a.set_X(angularoffset.X()); // retain length
  angularoffset = a; // modify
  RETURN_DSB_ADJUSTED;
}

void spherical_boundary::direction_boundary_effect(spatial & P, spatial & angularoffset) {
  // In the case of the spherical boundary, boundary effects begin even before a segment
  // actually intersects the boundary. This is based on the hypothesis that the accumulation
  // of constrained fibers near a boundary leads to an increasing density toward the
  // boundary asymptote, so that there is increasing pressure to move away from that
  // physical boundary.
  // As also noted in the paper, if a spherical boundary of this type is placed entirely
  // outside of a set of cell bodies, it acts as an attractor in the environment for the
  // growing neurites of those cells.
  int i = 0;
  for (; ((direction_adjust(P,angularoffset,i)) && (i<MAX_PB_ADJUSTS)); i++) ; // iterative for strictness
  if (i>=MAX_PB_ADJUSTS) warning("Warning: Reached maximum number of attempts to adjust to physical boundary effects.\n");
  // [***NOTE] We currently deal with other environment physics in the linked list at the root level.
}

void point_attractor::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_maxrange"))>=0) maxrange = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_c_attract"))>=0) {
    double new_c_attract = atof(clp.ParValue(n));
    if ((new_c_attract>=-1.0) && (new_c_attract<=1.0)) c_attract = new_c_attract;
    else warning("Warning: The attraction coefficient of a point attractor ("+thislabel+".attractor_c_attract) must lie in [-1,1].\n");
  }
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_point"))>=0) {
    if (!eq) warning("Warning: In order to set "+thislabel+".attractor_point to a neuron index the Event_Queue must be initialized first.\n");
    else if (!eq->Net()) warning("Warning: In order to set "+thislabel+".attractor_point to a neuron index the network must be initialized first.\n");
    else {
      neuron * nptr = static_cast<PLLRoot<neuron> *>(eq->Net())->el(atoi(clp.ParValue(n)));
      if (nptr) attractor = nptr->Pos();
      else warning("Parameter ERROR, "+thislabel+".attractor_point: Neuron index "+clp.ParValue(n)+" does not exist in network.\n");
    }
  }
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_pointX"))>=0) attractor.set_X(atof(clp.ParValue(n)));
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_pointY"))>=0) attractor.set_Y(atof(clp.ParValue(n)));
#ifdef VECTOR3D
  if ((n=clp.Specifies_Parameter(thislabel+".attractor_pointZ"))>=0) attractor.set_Z(atof(clp.ParValue(n)));
#endif
  if (Root()) static_cast<environment_physics *>(Root())->select_physical_boundary(thislabel,clp); // check for more with thislabel as prefix
}

String point_attractor::report_parameters() {
  String res("Point Attractor: attractor=(");
  res += String(attractor.X(),"%.3f,");
  res += String(attractor.Y(),"%.3f");
#ifdef VECTOR3D
  res += String(attractor.Z(),",%.3f");
#endif
  res += ") maxrange=";
  res += String(maxrange,"%.3f");
  res += " c_attract=";
  res += String(c_attract,"%.3f");
  return res;
}

bool point_attractor::direction_adjust(spatial & P, spatial & angularoffset) { //, int passnum) {
  // The point attractor is fairly simple, but it does discriminate between
  // "interested" and "uninterested" growth cones. Only those already getting closer
  // to a point attractor are considered interested and are further encouraged.
  if (angularoffset.X()==0.0) return false; // catch this special case faster than after computing d_approach
  spatial a(angularoffset); // local working copy
  a.convert_from_spherical();
  spatial Pnew(P);
  Pnew += a;
  double d_Pnew = Pnew.distance(attractor);
  //if ((passnum>0) && (d_Pnew>=0.0)) RETURN_DSB_IN_AFTER_ADJUST; // after adjustment, within sphere, no more adjustments needed
  if (d_Pnew > maxrange) return false; //RETURN_DSB_SAFE;
  double d_P = P.distance(attractor);
  //if ((passnum>0) && (d_P<0.0)) RETURN_DSB_FROM_BEYOND; // after adjustment of growth cone originating outside sphere, no more than one adjustment
  double d_approach = d_P - d_Pnew; // if positive, then got closer
  if (d_approach <= 0.0) return false; //RETURN_DSB_NO_MOVE_TO_BORDER; "uninterested" growth cone
  // the following is applied to growth cones that have been heading toward the attractor, and are within the attractor maxrange field
  double attract;
  double fieldratio = d_Pnew / maxrange;
  attract = d_approach*c_attract*(1 - (fieldratio*fieldratio));
  spatial attractvec(attractor);
  attractvec -= Pnew;
  attractvec.normalize();
  attractvec *= attract;
  a += attractvec;
  a.convert_to_spherical();
  a.set_X(angularoffset.X()); // retain length
  angularoffset = a; // modify
  return true; //RETURN_DSB_ADJUSTED; 
}

void point_attractor::direction_boundary_effect(spatial & P, spatial & angularoffset) {
  direction_adjust(P,angularoffset);
  /*int i = 0;
  for (; ((direction_adjust(P,angularoffset,i)) && (i<MAX_PB_ADJUSTS)); i++) ; // iterative for strictness
  if (i>=MAX_PB_ADJUSTS) warning("Warning: Reached maximum number of attempts to adjust to point attractor effects.\n");*/
  // [***NOTE] We currently deal with other environment physics in the linked list at the root level.
}

int environment_physics::select_physical_boundary(String prefixlabel, Command_Line_Parameters & clp) {
  // Note that handling labels is implemented slightly differently for the chains created here
  // than in the standard chaining method of models inherited from schemas and parents.
  // Note that valid calls to this function are the first call with prefixlabel="", subsequent
  // calls from similar existing objects looking for more objects to add to the linked-list.
  // Those calls can specify one or many labels (separated by spaces).
  // Returns the number of environment physics objects added to the linked list.
  int n;
  if (!prefixlabel.empty()) if (prefixlabel.lastchar()!='.') prefixlabel += '.';
  if ((n=clp.Specifies_Parameter(prefixlabel+"environment_physics"))>=0) {
    String epsstr(clp.URI_unescape_ParValue(n));
    int added = 0;
    while (!epsstr.empty()) {
      String epname(epsstr.before(' '));
      if (epname.empty()) {
	epname = epsstr;
	epsstr = "";
      } else epsstr = epsstr.after(' ');
      if ((n=clp.Specifies_Parameter(epname+".physical_boundary"))>=0) {
	String downcaselabel(downcase(clp.ParValue(n)));
	bool recognized = false;
	for (physical_boundary_types pb = (physical_boundary_types) 0; pb < NUM_pb; pb = (physical_boundary_types) (pb + 1)) {
	  if (downcaselabel.matches(physical_boundary_label[pb])) {
	    switch (pb) {
	    case spherical_pb:
	      link_before(new spherical_boundary(epname,clp));
	      recognized = true;
	      added++;
	      break;;
	    case point_attractor_pb:
	      link_before(new point_attractor(epname,clp));
	      recognized = true;
	      added++;
	      break;;
	    case NUM_pb: break;;
	      //default: ;
	    }
	    if (recognized) break;
	  }
	}
	if (!recognized) warning("Warning: Unknown physical boundary with label '"+downcaselabel+"', not applied\n");
      } else warning("Warning: Environment physics label "+epname+" unused\n");
    }
    return added;
  }
  return 0;
}
