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
// vector.hh
// Randal A. Koene, 20050112
// Vector and derived  classes that select 2D or 3D vectors and points
// according to the compiler directive VECTOR

#ifndef __SPATIAL_HH
#define __SPATIAL_HH

// Please note:
// Two different optimization goals are strived for with the use of these
// classes:
// (1) Two maintain only one source file for each set of classes and
//     functions, while compiling both nibr and nibr2D progreams.
// (2) To use inline compilation where possible to improve efficiency.
//
// Goal (1) is the reason why this header file and those that include
// it attempt to avoid compiler directives and the inclusion of vector
// header code in other header files, as well as the reason why
// compiler directives are used in those source files where direct access
// to spatial variables or functions is needed and depends on the number
// of dimensions.
//
// Goal (2) is the reason why this header file cannot simply provide a
// single interface (i.e. reference different vector headers only in
// the .cc file). Only by including the compiler directives and code
// differences in this header file to some degree is it possible to
// use inline compilation for some of the vector code.
//
// The result: Maintain combined source files and header files. Produce
// distinct 3D and 2D object files where necessary for linkage in
// nibr and nibr2D.

// 3D Display Options:
// It is possible to use gnuplot (with all its conversion capabilities,
// as shown in the scripts when images are clicked on at
// http://ayapin.film.s.dendai.ac.jp/~matuda/Gnuplot/pm3d.html)
// to draw the network in 3D by using spheres and pipes.
//
// AC3D format is also quite simple, so I can at least export in that format
// as well, even if actual manual manipulations in X are needed to use it
// (since AC3D does not have a command line mode).
//
// One of the most useful things is a simple 3D interface to
// Fig_Objects, including some rotation and translation matrices. That may
// require much smaller files than the GNUplot approach and can produce
// high quality figures for publication.

#ifdef VECTOR3D
#include "Vector3D.hh"
#include "Point3D.hh"
#define Vector Vector3D
#define Point Point3D
#define dot(u,v)   ((u).X() * (v).X() + (u).Y() * (v).Y() + (u).Z() * (v).Z())

class Spatial_Presentation;
#endif
#ifdef VECTOR2D
#include "Vector2D.hh"
#include "Point2D.hh"
#define Vector Vector2D
#define Point Point2D
#define dot(u,v)   ((u).X() * (v).X() + (u).Y() * (v).Y())
#endif

#define SMALL_NUM  0.00000001 // anything that avoids division overflow
//#define Pabs(x)     ((x) >= 0 ? (x) : -(x))   // absolute value
#define norm(v)    sqrt(dot(v,v))  // norm = length of vector
//#define Pmax(a,b) ((a) >= (b) ? (a) : (b)) // maximum of two values

// VECTOR2D or VECTOR3D should be supplied by a compiler directive

extern unsigned long int zerocoordinateconversions;

class spatial: public Vector {
public:
  spatial() {}
  spatial(const spatial & s): Vector(s) {}
#ifdef VECTOR3D
  spatial(double X, double Y, double Z = 0.0): Vector(X,Y,Z) {}
  void set_all(double X, double Y, double Z = 0.0) { x=X; y=Y; z=Z; }
  void get_all(double & X, double & Y, double & Z) { X=x; Y=y; Z=z; }
  void plane_mapped(double & px, double & py);
#endif
#ifdef VECTOR2D
  spatial(double X, double Y, double Z = 0.0): Vector(X,Y) {}
  void set_all(double X, double Y, double Z = 0.0) { x=X; y=Y; }
  void get_all(double & X, double & Y, double & Z) { X=x; Y=y; Z=0.0; }
  void plane_mapped(double & px, double & py) { px = x; py = y; }
#endif
  //operator Vector () { return ::Vector(); }
  //spatial operator-( Vector3D);        // vector subtract
  void add_to_all(double m) { // scalar addition to all elements
    x += m;
    y += m;
#ifdef VECTOR3D
    z += m;
#endif
  }
  void convert_to_spherical();
  void convert_from_spherical();
  void RandomAngular(double minlen, double maxlen, double minangle, double maxangle);
  void AddRandomAngular(double minlen, double maxlen, double minangle, double maxangle);
  void AddRandomAngularAllConstrained(double minlen, double maxlen, double minangle, double maxangle);
  void AngularOpposite();
  void Negate();
};

// [***Note] For the angle between two vectors see the turn() function in
// dendritic_growth_model.cc.

class Segment {
public:
  spatial P0, P1;
  Segment() {}
  //Segment(spatial p0, spatial p1): P0(p0), P1(p1) {}
  Segment(spatial & p0, spatial & p1): P0(p0), P1(p1) {}
  bool iszerolength() { return (P0==P1); }
  void set_by_angular_coordinates(spatial & acoords);
  void get_angles_and_length(spatial & acoords);
  double Length();
};

class Plane {
public:
  spatial V0, n;
  Plane() {};
  Plane(spatial & v0, spatial & N): V0(v0), n(N) {}
  Plane(double vx, double nx, double vy, double ny, double vz = 0.0, double nz = 0.0): V0(vx,vy,vz), n(nx,ny,nz) {}
};

class Rotation {
  // Rotations are achieved by a rotation matrix that is the result of
  // multiplication of rotation matrices Rx Ry Rz, each of which applies
  // rotation about one the of x, y or z axes.
  // Full rotation around an origin is accomplished by translating to the
  // origin, then rotating and translating back.
  // The rotate_extract_view() function only translated to the origin,
  // then applies rotation to the x and y coordinates so that a rotated
  // spatial coordinates is mapped to the x-y plane for display.
protected:
  spatial xr, yr, zr; // row vectors of the rotation matrix
  spatial origin;
  spatial angles;
public:
  Rotation() { set_rotation(angles); }
  Rotation(const spatial & a, const spatial & o) { set_rotation(a,o); }
  Rotation(const spatial & a) { set_rotation(a); }
  void set_rotation(const spatial & a, const spatial & o);
  void set_rotation(const spatial & a);
  spatial & Angles() { return angles; }
  spatial & Origin() { return origin; }
  spatial & XR() { return xr; }
  spatial & YR() { return yr; }
  spatial & ZR() { return zr; }
  void rotate(const spatial & s, spatial & rotated);
  void rotate_extract_view(const spatial & s, double & x, double & y);
};

class Projection {
  // This object contains parameters, matrices of the parameters, and functions
  // used to create 3D perspective projections (see ./wikipedia.3D-projection.pdf).
protected:
  spatial xp, yp, zp; // row vectors of the projection matrix
  spatial cameralocation;
  spatial camerarotation;
  spatial viewerposition;
public:
  Projection() { set_projection(); }
  Projection(const spatial & cl, const spatial & cr, const spatial & vp) { set_projection(cl,cr,vp); }
  void set_projection();
  void set_projection(const spatial & cl, const spatial & cr, const spatial & vp);
  void set_camera_location(const spatial & cl) { cameralocation = cl; set_projection(); }
  void set_camera_rotation(const spatial & cr) { camerarotation = cr; set_projection(); }
  void set_viewer_position(const spatial & vp) { viewerposition = vp; set_projection(); }
  spatial & CameraLocation() { return cameralocation; }
  spatial & CameraRotation() { return camerarotation; }
  spatial & ViewerPosition() { return viewerposition; }
  //  void project(const spatial & s, spatial & projected);
  bool project_extract_view(const spatial & s, double & x, double & y);
};

double dist3D_Segment_to_Segment( Segment S1, Segment S2, Segment * distbar = NULL);

#ifdef VECTOR3D
// This is normally equal to spat_pres (see Spatial_Presentation.hh).
extern Spatial_Presentation * spat_map; // global mapping of the spatial environment
#endif

#endif
