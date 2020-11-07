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
// Spatial_Presentation.hh
// Randal A. Koene, 20050218
//
// Classes that provide control over the spatial presentation of output.

#ifndef __SPATIAL_PRESENTATION_HH
#define __SPATIAL_PRESENTATION_HH

#include <math.h>
#include "Command_Line_Parameters.hh"
#include "spatial.hh"

#ifdef VECTOR3D
class Spatial_Presentation: public Rotation, public Projection, public CLP_Modifiable {
protected:
  spatial rotateinterval, rotatestepsize;
  bool useprojection;
public:
  Spatial_Presentation(const spatial & a, const spatial & o): Rotation(a,o), rotateinterval(2.0*M_PI,2.0*M_PI,2.0*M_PI), rotatestepsize(1.0,1.0,1.0), useprojection(false) {}
  //Spatial_Presentation(const spatial & a): Rotation(a), rotateinterval(2.0*M_PI,2.0*M_PI,2.0*M_PI), rotatestepsize(1.0,1.0,1.0), useprojection(false) {}
  Spatial_Presentation(): rotateinterval(2.0*M_PI,2.0*M_PI,2.0*M_PI), rotatestepsize(1.0,1.0,1.0), useprojection(false) {
    spatial a(M_PI/4.0,M_PI/4.0,M_PI/4.0);
    set_rotation(a);
  }
  Spatial_Presentation(const spatial & o): rotateinterval(2.0*M_PI,2.0*M_PI,2.0*M_PI), rotatestepsize(1.0,1.0,1.0), useprojection(false) {
    spatial a(M_PI/4.0,M_PI/4.0,M_PI/4.0);
    set_rotation(a,o);
  }
  virtual ~Spatial_Presentation() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  void set_by_angles(double ax, double ay, double az, const spatial & o);
  void set_by_angles(const spatial & a, const spatial & o) { set_rotation(a,o); }
  void add_angular(double ax, double ay, double az);
  void add_angular(const spatial & a);
  void set_rotate_speed(int numsteps);
  void update_rotation() { add_angular(rotatestepsize); }
  void extract_view(const spatial & s, double & x, double & y) { if (useprojection) project_extract_view(s,x,y); else rotate_extract_view(s,x,y); }
};

bool get_spatial(Command_Line_Parameters & clp, String prefix, spatial & s);

extern Spatial_Presentation * spat_pres; // global presentation of the entire spatial environment
#endif

#endif

