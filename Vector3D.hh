//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#ifndef SS_Vector3D_H
#include <cmath>
#include "Point3D.hh"
#define SS_Vector3D_H

class Vector3D : public Point3D {
public:
  Vector3D(): Point3D() {};
  Vector3D(double X, double Y, double Z) : Point3D(X,Y,Z) {};
  //Vector3D(Vector3D & v): Point3D(v.x,v.y,v.z) {};
  //*** Using implicit copy constructor
  ~Vector3D() {};

  Vector3D operator-();                // unary minus
  Vector3D operator~();                // unary 2D perp operator

  // Scalar Multiplication
  friend Vector3D operator*( int, Vector3D);
  friend Vector3D operator*( double, Vector3D);
  friend Vector3D operator*( Vector3D, int);
  friend Vector3D operator*( Vector3D, double);
  // Scalar Division
  friend Vector3D operator/( Vector3D, int);
  friend Vector3D operator/( Vector3D, double);

  // Vector3D Arithmetic Operations
  Vector3D operator+( Vector3D);        // vector add
  Vector3D operator-( Vector3D);        // vector subtract
  double operator*( Vector3D w) { return (x * w.x + y * w.y + z * w.z); } // inner dot product
  double operator|( Vector3D w) { return (x * w.y - y * w.x); } // 2D exterior perp product
  Vector3D operator^( Vector3D);        // 3D exterior cross product
  
  Vector3D& operator*=(double c) { // vector scalar mult
    x *= c;
    y *= c;
    z *= c;
    return *this;
  }
  Vector3D& operator/=(double c) { // vector scalar div
    x /= c;
    y /= c;
    z /= c;
    return *this;
  }
  Vector3D& operator+=( Vector3D w) { // vector increment
    x += w.x;
    y += w.y;
    z += w.z;
    return *this;
  }
  Vector3D& operator-=( Vector3D w) { // vector decrement
    x -= w.x;
    y -= w.y;
    z -= w.z;
    return *this;
  }
  Vector3D& operator^=( Vector3D);      // 3D exterior cross product
  // operator== inherited from Point

  double len() { return sqrt(x*x + y*y + z*z); } // vector length
  double len2() { return (x*x + y*y + z*z); } // vector length squared (faster)
  // Special Operations
  double normalize();                 // convert vector to unit length
  friend Vector3D sum( int, int[], Vector3D[]);     // vector sum
  friend Vector3D sum( int, double[], Vector3D[]);  // vector sum
};

#endif
