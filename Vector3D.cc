//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050112

#include "Point3D.hh"
#include "Vector3D.hh"

Vector3D Vector3D::operator-() {
  Vector3D v;
  v.x = -x; v.y = -y; v.z = -z;
  return v;
}

// Unary 2D perp operator
Vector3D Vector3D::operator~() {
  // continue even if not 2D
  Vector3D v;
  v.x = -y; v.y = x; v.z = z;
  return v;
}

Vector3D operator*( int c, Vector3D w ) {
  Vector3D v;
  v.x = c * w.x;
  v.y = c * w.y;
  v.z = c * w.z;
  return v;
}

Vector3D operator*( double c, Vector3D w ) {
  Vector3D v;
  v.x = c * w.x;
  v.y = c * w.y;
  v.z = c * w.z;
  return v;
}

Vector3D operator*( Vector3D w, int c ) {
  Vector3D v;
  v.x = c * w.x;
  v.y = c * w.y;
  v.z = c * w.z;
  return v;
}

Vector3D operator*( Vector3D w, double c ) {
  Vector3D v;
  v.x = c * w.x;
  v.y = c * w.y;
  v.z = c * w.z;
  return v;
}

// Scalar division
Vector3D operator/( Vector3D w, int c ) {
  Vector3D v;
  v.x = w.x / c;
  v.y = w.y / c;
  v.z = w.z / c;
  return v;
}

Vector3D operator/( Vector3D w, double c ) {
  Vector3D v;
  v.x = w.x / c;
  v.y = w.y / c;
  v.z = w.z / c;
  return v;
}

Vector3D Vector3D::operator+( Vector3D w ) {
  Vector3D v;
  v.x = x + w.x;
  v.y = y + w.y;
  v.z = z + w.z;
  return v;
}

Vector3D Vector3D::operator-( Vector3D w ) {
  Vector3D v;
  v.x = x - w.x;
  v.y = y - w.y;
  v.z = z - w.z;
  return v;
}

// 3D Exterior Cross Product
Vector3D Vector3D::operator^( Vector3D w ) {
  Vector3D v;
  v.x = y * w.z - z * w.y;
  v.y = z * w.x - x * w.z;
  v.z = x * w.y - y * w.x;
  return v;
}

Vector3D& Vector3D::operator^=( Vector3D w ) { // 3D exterior cross product
  double ox=x, oy=y, oz=z;
  x = oy * w.z - oz * w.y;
  y = oz * w.x - ox * w.z;
  z = ox * w.y - oy * w.x;
  return *this;
}

double Vector3D::normalize() { 
  double ln = sqrt( x*x + y*y + z*z );
  if (ln == 0.0) return 0.0;
  x /= ln;
  y /= ln;
  z /= ln;
  return ln;
}

Vector3D sum( int n, int c[], Vector3D w[] ) {     // vector sum
  Vector3D  v;
  for (int i=0; i<n; i++) {
    v.x += c[i] * w[i].x;
    v.y += c[i] * w[i].y;
    v.z += c[i] * w[i].z;
  }
  return v;
}

Vector3D sum( int n, double c[], Vector3D w[] ) {  // vector sum
  Vector3D  v;
  for (int i=0; i<n; i++) {
    v.x += c[i] * w[i].x;
    v.y += c[i] * w[i].y;
    v.z += c[i] * w[i].z;
  }
  return v;
}
