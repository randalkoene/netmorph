//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050112

#include "Point2D.hh"
#include "Vector2D.hh"

Vector2D Vector2D::operator-() {
  Vector2D v;
  v.x = -x; v.y = -y;
  return v;
}

// Unary 2D perp operator
Vector2D Vector2D::operator~() {
  // continue even if not 2D
  Vector2D v;
  v.x = -y; v.y = x;
  return v;
}

Vector2D operator*( int c, Vector2D w ) {
  Vector2D v;
  v.x = c * w.x;
  v.y = c * w.y;
  return v;
}

Vector2D operator*( double c, Vector2D w ) {
  Vector2D v;
  v.x = c * w.x;
  v.y = c * w.y;
  return v;
}

Vector2D operator*( Vector2D w, int c ) {
  Vector2D v;
  v.x = c * w.x;
  v.y = c * w.y;
  return v;
}

Vector2D operator*( Vector2D w, double c ) {
  Vector2D v;
  v.x = c * w.x;
  v.y = c * w.y;
  return v;
}

// Scalar division
Vector2D operator/( Vector2D w, int c ) {
  Vector2D v;
  v.x = w.x / c;
  v.y = w.y / c;
  return v;
}

Vector2D operator/( Vector2D w, double c ) {
  Vector2D v;
  v.x = w.x / c;
  v.y = w.y / c;
   return v;
}

Vector2D Vector2D::operator+( Vector2D w ) {
  Vector2D v;
  v.x = x + w.x;
  v.y = y + w.y;
  return v;
}

Vector2D Vector2D::operator-( Vector2D w ) {
  Vector2D v;
  v.x = x - w.x;
  v.y = y - w.y;
  return v;
}

double Vector2D::normalize() { 
  double ln = sqrt(x*x + y*y);
  if (ln == 0.0) return 0.0;
  x /= ln;
  y /= ln;
  return ln;
}

Vector2D sum( int n, int c[], Vector2D w[] ) {     // vector sum
  Vector2D  v;
  for (int i=0; i<n; i++) {
    v.x += c[i] * w[i].x;
    v.y += c[i] * w[i].y;
  }
  return v;
}

Vector2D sum( int n, double c[], Vector2D w[] ) {  // vector sum
  Vector2D  v;
  for (int i=0; i<n; i++) {
    v.x += c[i] * w[i].x;
    v.y += c[i] * w[i].y;
  }
  return v;
}
