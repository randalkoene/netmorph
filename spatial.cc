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
// spatial.cc
// Randal A. Koene, 20050112

#include <math.h>
//#include <errno.h>
#include "spatial.hh"
#include "Spatial_Presentation.hh"
#include "global.hh"

// special interfaces (e.g. calls by Fig_Object functions)

unsigned long int zerocoordinateconversions = 0;

#ifdef VECTOR3D
Spatial_Presentation * spat_map = NULL;

void spatial::plane_mapped(double & px, double & py) {
  // Project onto a plain, taking into acount rotations and such.
  // Adapted from the intersect3D_SegmentPlane() example in
  // ~/download/manuals/geometry-algorithms.softsurfer.com.tar.gz.
  // given: point in view plane Pn.V0, normal vector of view plane, Pn.n
  //        and of course a point (this)
  //spatial w(*this); w -= viewplane.V0;
  //double D = viewplane.n.len2(); // *** dot(viewplane.n,u);
  //double N = -dot(viewplane.n,w);
  //double sl = N / D;
  //px = x + (sl * viewplane.n.X());
  //py = y + (sl * viewplane.n.Y());
  spat_map->extract_view(*this,px,py);
}
#endif

void spatial::convert_to_spherical() {
  // 3D spherical: x is rho, y is theta, z is phi
  // 2D polar: x is r, y is theta
  // If the spatial vector is a zero vector, then a zero angular vector is returned.
  double S = (x*x) + (y*y);
#ifdef TEST_FOR_NAN
  bool showtest = false;
  if ((x>-0.00001) && (x<0.00001)) {
    showtest = true;
    cout << "Converting ZERO: x=" << x << ", y=" << y << ", S=" << S << '\n'; cout.flush();
  }
#endif
#ifdef VECTOR3D
  double rho = S + (z*z);
#endif
  S = sqrt(S);
#ifdef VECTOR2D
  y = atan2(y,x);
  x = S;
#endif
#ifdef TEST_FOR_NAN
  if (showtest) { cout << "After some steps: x=" << x << ", y=" << y << ", S=" << S << '\n'; cout.flush(); }
#endif
#ifdef VECTOR3D
#ifdef TEST_FOR_NAN
  double rhosquared = rho;
#endif
  if ((rho>-DBL_EPSILON) && (rho<DBL_EPSILON)) {
    z = 0.0;
    y = 0.0;
    x = 0.0;
    zerocoordinateconversions++;
  } else {
    rho = sqrt(rho);
#ifdef TEST_FOR_NAN
    String rhonanstr(rho,"%f");
    if (rhonanstr==String("nan")) { cout << "convert_to_spherical rho==NAN (" << rho << "), rhosquared=" << rhosquared << ", S=" << S << ", z=" << z << ", x=" << x << ", y=" << y << '\n'; cout.flush(); }
#endif
    z = acos(z/rho);
    if (x>=0) y = asin(y/S);
    else y = M_PI - asin(y/S);
    x = rho;
  }
#endif
#ifdef TEST_FOR_NAN
  String nanstr(x,"%f");
  if (nanstr==String("nan")) { cout << "convert_to_spherical resulted in NAN!\n"; cout.flush(); }
#endif
}

void spatial::convert_from_spherical() {
  // In 3D, theta (y) is assumed to go fro 0 to 2*M_PI (it is actually
  // unrestricted). It is the angle between the x axis and the projection
  // of the spatial vector on the x-y plane.
  // Phi (z) is assumed to go from 0 to M_PI. It is the angle between the
  // z axis and the spatial vector. (This applies only to 3D.)
  double lensinphi = x;
#ifdef VECTOR3D
  lensinphi *= sin(z);
  z = x*cos(z);
#endif
  x = lensinphi*cos(y);
  y = lensinphi*sin(y);
}

void spatial::RandomAngular(double minlen, double maxlen, double minangle, double maxangle) {
  // This sets a random length between minlen and maxlen, and random
  // angles between minangle and maxangle.
  x = X_misc.get_rand_range_real1(minlen,maxlen); // [***NOTE] Should the random generator be specified as a parameter to this function?
  y = X_misc.get_rand_range_real1(minangle,maxangle);
#ifdef VECTOR3D
  if (minangle<0.0) minangle=0.0;
  if (minangle>M_PI) minangle=M_PI;
  if (maxangle<0.0) maxangle=0.0;
  if (maxangle>M_PI) maxangle=M_PI;
  z = X_misc.get_rand_range_real1(minangle,maxangle);
#endif
}

void spatial::AddRandomAngular(double minlen, double maxlen, double minangle, double maxangle) {
  // This sets a random length between minlen and maxlen, and adds random
  // angles between minangle and maxangle to angle theta and randomly
  // rotates phi.
  x = X_misc.get_rand_range_real1(minlen,maxlen); // [***NOTE] Should the random generator be specified as a parameter to this function?
  y += X_misc.get_rand_range_real1(minangle,maxangle);
#ifdef VECTOR3D
  /*if (minangle<0.0) minangle=0.0;
  if (minangle>M_PI) minangle=M_PI;
  if (maxangle<0.0) maxangle=0.0;
  if (maxangle>M_PI) maxangle=M_PI;*/
  z = X_misc.get_rand_range_real1(0.0,M_PI);
#endif
}

void spatial::AddRandomAngularAllConstrained(double minlen, double maxlen, double minangle, double maxangle) {
  // This sets a random length between minlen and maxlen, and adds random
  // angles between minangle and maxangle to angles theta and phi. This
  // constrains the deviation of phi just as the deviation of theta is
  // constrained. In 2D this function is identical to AddRandomAngular().
  x = X_misc.get_rand_range_real1(minlen,maxlen); // [***NOTE] Should the random generator be specified as a parameter to this function?
  y += X_misc.get_rand_range_real1(minangle,maxangle);
#ifdef VECTOR3D
  z += X_misc.get_rand_range_real1(minangle,maxangle);
  if (z<0) {
    y += M_PI;
    z = -z;
    return;
  }
  if (z>=M_PI) {
    y += M_PI;
    z = (2*M_PI) - z;
    return;
  }
#endif
}

void spatial::AngularOpposite() {
  // This assumes that a spatial location is given in angular coordinates and
  // computes its opposite as determined along a normal vector through the
  // origin.
  y += M_PI;
#ifdef VECTOR3D
  z = M_PI - z;
#endif
}

void spatial::Negate() {
  // This function inverts the sign of all elements. When the spatial location
  // is given in Euler cartesian coordinates this computes the same opposite
  // as the AngularOpposite() function above.
  x = -x;
  y = -y;
#ifdef VECTOR3D
  z = -z;
#endif
}

void Segment::set_by_angular_coordinates(spatial & acoords) {
  // This sets up a segment for which P0 is already defined as
  // P1 = P0 + Euler(acoords), where Euler(x) is a conversion from
  // spherical/polar coordinates to Euler cartesian coordinates.
  // In 3D, theta (acoords[1], ie. acoords.Y()) is assumed to go from
  // 0 to 2*M_PI (it is actually unrestricted). It is the angle between
  // the x axis and the projection of the spatial vector on the x-y plane.
  // Phi (acoords[2], ie. acoords.Z()) is assumed to go from 0 to M_PI. It
  // is the angle between the z axis and the spatial vector.
  double lensinphi = acoords.X();
#ifdef VECTOR3D
  P1.set_Z(lensinphi*cos(acoords.Z()));
  lensinphi *= sin(acoords.Z());
#endif
  P1.set_X(lensinphi*cos(acoords.Y()));
  P1.set_Y(lensinphi*sin(acoords.Y()));
  P1 += P0;
}

void Segment::get_angles_and_length(spatial & acoords) {
  // This subtracts P0 from P1 and then converts the result into
  // angles and length to return vector data when P0 is the origin.
  acoords = P1;
  acoords -= P0;
  /* or acoords = P0; acoords.Negate(); acoords += P1; */
  acoords.convert_to_spherical();
}

double Segment::Length() {
  // Length as Euclidean distance between points P0 and P1 of the Segment.
  double dx = P0.X() - P1.X();
  double dy = P0.Y() - P1.Y();
#ifdef VECTOR3D
  double dz = P0.Z() - P1.Z();
  return sqrt(dx*dx + dy*dy + dz*dz);
#endif
#ifdef VECTOR2D
  return sqrt(dx*dx + dy*dy);
#endif
}

void Rotation::set_rotation(const spatial & a) {
  angles = a;
  double sinax = sin(a.X()), cosax = cos(a.X());
  double sinay = sin(a.Y()), cosay = cos(a.Y());
#ifdef VECTOR3D
  double sinaz = sin(a.Z()), cosaz = cos(a.Z());
#endif
#ifdef VECTOR2D
  double sinaz = 0.0, cosaz = 1.0;
#endif
  spatial // Rx in rows
    rxx(1.0,0.0,0.0), 
    rxy(0.0,cosax,sinax), 
    rxz(0.0,-sinax,cosax);
  spatial //Ry in rows
    ryx(cosay,0.0,-sinay), 
    ryy(0.0,1.0,0.0), 
    ryz(sinay,0.0,cosay);
  spatial // Rz in columns
    czx(cosaz,
	-sinaz,
	0.0), 
    czy(sinaz,
	cosaz,
	0.0), 
    czz(0.0,
	0.0,
	1.0);
  spatial // Ry Rz in columns
    cyzx(dot(ryx,czx),
	 dot(ryy,czx),
	 dot(ryz,czx)), // Ry czx
    cyzy(dot(ryx,czy),
	 dot(ryy,czy),
	 dot(ryz,czy)), // Ry czy
    cyzz(dot(ryx,czz),
	 dot(ryy,czz),
	 dot(ryz,czz)); // Ry czz
  // Rx (Ry Rz) in rows
  xr.set_all(dot(rxx,cyzx),dot(rxx,cyzy),dot(rxx,cyzz));
  yr.set_all(dot(rxy,cyzx),dot(rxy,cyzy),dot(rxy,cyzz));
  zr.set_all(dot(rxz,cyzx),dot(rxz,cyzy),dot(rxz,cyzz));
}

void Rotation::set_rotation(const spatial & a, const spatial & o) {
  origin = o;
  set_rotation(a);
}

void Rotation::rotate(const spatial & s, spatial & rotated) {
  spatial translated(s); translated -= origin;
  double x = dot(xr,translated);
  double y = dot(yr,translated);
  double z = dot(zr,translated);
  rotated.set_all(x,y,z);
  rotated += origin;
}

void Rotation::rotate_extract_view(const spatial & s, double & x, double & y) {
  //cout << origin.X() << ',' << origin.Y() << '\n';
  spatial translated(s); translated -= origin;
  x = dot(xr,translated);
  y = dot(yr,translated);
}

double dist3D_Segment_to_Segment( Segment S1, Segment S2, Segment * distbar) {
  // Input:  two 3D line segments S1 and S2
  // Return: the shortest distance between S1 and S2
  // distbar can return a segment for the nearest points between S1 and S2
  spatial u(S1.P1); u -= S1.P0;
  spatial v(S2.P1); v -= S2.P0;
  spatial w(S1.P0); w -= S2.P0;
  double    a = dot(u,u);        // always >= 0
  double    b = dot(u,v);
  double    c = dot(v,v);        // always >= 0
  double    d = dot(u,w);
  double    e = dot(v,w);
  double    D = a*c - b*b;       // always >= 0
  double    sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
  double    tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0
                                    
  // compute the line parameters of the two closest points
  if (D < SMALL_NUM) { // the lines are almost parallel
    sN = 0.0;        // force using point P0 on segment S1
    sD = 1.0;        // to prevent possible division by 0.0 later
    tN = e;
    tD = c;
  } else {                // get the closest points on the infinite lines
    sN = (b*e - c*d);
    tN = (a*e - b*d);
    if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
      sN = 0.0;
      tN = e;
      tD = c;
    } else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
      sN = sD;
      tN = e + b;
      tD = c;
    }
  }
  if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
    tN = 0.0;
    // recompute sc for this edge
    if (-d < 0.0) sN = 0.0;
    else if (-d > a) sN = sD;
    else {
      sN = -d;
      sD = a;
    }
  } else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
    tN = tD;
    // recompute sc for this edge
    if ((-d + b) < 0.0) sN = 0;
    else if ((-d + b) > a) sN = sD;
    else {
      sN = (-d + b);
      sD = a;
    }
  }
  // finally do the division to get sc and tc
  double abssN = Pabs(sN);
  double abstN = Pabs(tN);
  if (abssN < SMALL_NUM) sc = 0.0;
  else sc = sN / sD;
  if (abstN < SMALL_NUM) tc = 0.0;
  else tc = tN / tD;
  // get the difference of the two closest points
  u *= sc; v *= tc;
  Vector dP = w + u - v;  // = S1(sc) - S2(tc)
  if (distbar) {
    distbar->P0 = S1.P0; distbar->P0 += u;
    distbar->P1 = S2.P0; distbar->P1 += v;
  }
  return norm(dP);   // return the closest distance
}

void Projection::set_projection() {
  double coscrx = cos(camerarotation.X()), sincrx = sin(camerarotation.X());
  double coscry = cos(camerarotation.Y()), sincry = sin(camerarotation.Y());
#ifdef VECTOR3D
  double coscrz = cos(camerarotation.Z()), sincrz = sin(camerarotation.Z());
#endif
#ifdef VECTOR2D
  double coscrz = 1.0, sincrz = 1.0;
#endif
  xp.set_all(coscry*coscrz,coscry*sincrz,-sincry);
  double sxsy = sincrx*sincry, cxsy = coscrx*sincry;
  yp.set_all((sxsy*coscrz)-(coscrx*sincrz),(sxsy*sincrz)+(coscrx*coscrz),sincrx*coscry);
  zp.set_all((cxsy*coscrz)+(sincrx*sincrz),(cxsy*sincrz)-(sincrx*coscrz),coscrx*coscry);
}

void Projection::set_projection(const spatial & cl, const spatial & cr, const spatial & vp) {
  cameralocation = cl;
  camerarotation = cr;
  viewerposition = vp;
  set_projection();
}

//void Projection::project(const spatial & s, spatial & projected) {
//}

bool Projection::project_extract_view(const spatial & s, double & x, double & y) {
  spatial ac(s); ac-=cameralocation;
  double dx = dot(xp,ac);
  double dy = dot(yp,ac);
#ifdef VECTOR3D
  double dz = dot(zp,ac);
  if (dz==0.0) return false;
  double ezdivdz = viewerposition.Z()/dz;
#endif
#ifdef VECTOR2D
  double ezdivdz = 1.0;
#endif
  x = (dx - viewerposition.X())*ezdivdz;
  y = (dy - viewerposition.Y())*ezdivdz;
  return true;
}
