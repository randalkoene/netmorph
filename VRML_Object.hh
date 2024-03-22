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
// VRML_Object.hh
// Randal A. Koene, 20080226
//
// Classes that define VRML_Object and objects that inherit it.

// This format was developed for DIL#20050216080412.1. Instead
// of catering specifically to AC3D, as in the original task
// description, this code now aims at X3D output, the
// successor to VRML.

#ifndef __VRML_OBJECT_HH
#define __VRML_OBJECT_HH

#include <limits.h>
#include "BigString.hh"
#include "templates.hh"

// definitions

// functions

// classes

class VRML_Object: public PLLHandle<VRML_Object> {
public:
  VRML_Object() {}
  virtual String str() { return String(""); }
};

class VRML_Header: public VRML_Object {
protected:
  String text;
public:
  VRML_Header();
  virtual String str() { return text; }
};

class VRML_Group: public VRML_Object {
protected:
  String comment;
  PLLRoot<VRML_Object> vrmlobjects;
public:
  VRML_Group(String c): comment(c) {}
  VRML_Group() {}
  virtual String str();
  void Add_VRML_Object(VRML_Object * vrmlobject);
};

void X3D_cylinder_center(spatial & p0, spatial & p1, spatial & center);
bool X3D_SFrotation_to_vector(spatial & p0, spatial & p1, spatial & rotaxis, double & theta, double & len);

extern String * VRML_neuronlist;
extern String * VRML_synapselist;
extern String * VRML_fiberlist;
extern String * VRML_continuationnodelist;
extern String * VRML_bifurcationnodelist;
extern String * VRML_terminalgrowthconelist;
extern long VRML_neuronindex;
extern long VRML_synapseindex;
extern long VRML_nodeindex;

#endif
