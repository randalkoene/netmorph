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
// Spatial_Presentation.cc
// Randal A. Koene, 20050218

#include "Spatial_Presentation.hh"
#include "global.hh"

#ifdef VECTOR3D
Spatial_Presentation * spat_pres = NULL;

bool get_spatial(Command_Line_Parameters & clp, String prefix, spatial & s) {
  int n = -1;
  if ((n=clp.Specifies_Parameter(prefix))>=0) { 
    int dlen;
    double * dptr = str2vec(clp.ParValue(n),&dlen);
    if (dlen<3) {
      warning("Warning: "+prefix+" location incomplete, ignoring parameters.\n");
      return false;
    } else {
      s.set_all(dptr[0],dptr[1],dptr[2]);
      delete[] dptr;
    }
  } else {
    if ((n=clp.Specifies_Parameter(prefix+"_x"))>=0) {
      s.set_X(clp.get_double(n));
      if ((n=clp.Specifies_Parameter(prefix+"_y"))>=0) {
	s.set_Y(clp.get_double(n));
	if ((n=clp.Specifies_Parameter(prefix+"_z"))>=0) s.set_Z(clp.get_double(n));
      }
    }
  }
  return (n>=0);
}

void Spatial_Presentation::parse_CLP(Command_Line_Parameters & clp) {
  spatial s;
  useprojection = get_spatial(clp,"camera",s);
  if (useprojection) {
    cameralocation = s;
    report("Using PROJECTION viewing. Any Rotation view parameters are ignored.\n");
    if (get_spatial(clp,"camera_ROT",s)) camerarotation = s;
    if (get_spatial(clp,"viewer",s)) viewerposition = s;
    set_projection();
  } else {
    report("Using ROTATION viewing. Any Projection view parameters are ignored.\n");
    if (get_spatial(clp,"ROT_origin",s)) origin = s;
    if (get_spatial(clp,"ROT",s)) angles = s;
    if (get_spatial(clp,"ROT_interval",s)) rotateinterval = s;
    set_rotation(Angles());
  }
}

String Spatial_Presentation::report_parameters() {
  if (useprojection) {
    String res("Presentation Projection:\n  camera location = (");
    res += String(cameralocation.X(),"%.2f,")+String(cameralocation.Y(),"%.2f,")+String(cameralocation.Z(),"%.2f)\n");
    res += "  camera rotation = (" + String(camerarotation.X(),"%.2f,")+String(camerarotation.Y(),"%.2f,")+String(camerarotation.Z(),"%.2f)\n");
    res += "  viewer position = (" + String(viewerposition.X(),"%.2f,")+String(viewerposition.Y(),"%.2f,")+String(viewerposition.Z(),"%.2f)\n");
    return res;
  } else {
    double ox, oy, oz, ax, ay, az, rix, riy, riz;
    String res("Presentation origin of rotation: ROT_origin_<x,y,z>=(");
    Origin().get_all(ox,oy,oz);
    res += String(ox,"%.2f,") + String(oy,"%.2f,") + String(oz,"%.2f)\n");
    res += "Presentation initial axes rotations: ROT_<x,y,z>=(";
    Angles().get_all(ax,ay,az);
    res += String(ax,"%.2f,") + String(ay,"%.2f,") + String(az,"%.2f)\n");
    res += "Presentation axes rotation intervals: ROT_interval_<x,y,z>=(";
    rotateinterval.get_all(rix,riy,riz);
    res += String(rix,"%.2f,") + String(riy,"%.2f,") + String(riz,"%.2f)\n");
    return res;
  }
}

void Spatial_Presentation::set_by_angles(double ax, double ay, double az, const spatial & o) {
  spatial a(ax,ay,az);
  set_rotation(a,o);
}

void Spatial_Presentation::add_angular(const spatial & a) {
  spatial A(Angles());
  A += a;
  set_rotation(A);
}

void Spatial_Presentation::add_angular(double ax, double ay, double az) {
  spatial a(ax,ay,az);
  add_angular(a);
}

void Spatial_Presentation::set_rotate_speed(int numsteps) {
  rotatestepsize = rotateinterval;
  rotatestepsize /= ((double) numsteps);
}
#endif

