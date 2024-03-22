//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#ifndef SS_Vector2D_H
#include <cmath>
#include "Point2D.hh"
#define SS_Vector2D_H

class Vector2D : public Point2D {
public:
  Vector2D(): Point2D() {};
  Vector2D(double X, double Y) : Point2D(X,Y) {};
  //Vector2D(Vector2D & v): Point2D(v.x,v.y) {};
  //*** Using implicit copy constructor
  ~Vector2D() {};

  Vector2D operator-();                // unary minus
  Vector2D operator~();                // unary 2D perp operator

  // Scalar Multiplication
  friend Vector2D operator*( int, Vector2D);
  friend Vector2D operator*( double, Vector2D);
  friend Vector2D operator*( Vector2D, int);
  friend Vector2D operator*( Vector2D, double);
  // Scalar Division
  friend Vector2D operator/( Vector2D, int);
  friend Vector2D operator/( Vector2D, double);

  // Vector2D Arithmetic Operations
  Vector2D operator+( Vector2D);        // vector add
  Vector2D operator-( Vector2D);        // vector subtract
  double operator*( Vector2D w) { return (x * w.x + y * w.y); } // inner dot product
  double operator|( Vector2D w) { return (x * w.y - y * w.x); } // 2D exterior perp product
  
  Vector2D& operator*=(double c) { // vector scalar mult
    x *= c;
    y *= c;
    return *this;
  }
  Vector2D& operator/=(double c) { // vector scalar div
    x /= c;
    y /= c;
    return *this;
  }
  Vector2D& operator+=(Vector2D w) { // vector increment
    x += w.x;
    y += w.y;
    return *this;
  }
  Vector2D& operator-=(Vector2D w) { // vector decrement
    x -= w.x;
    y -= w.y;
    return *this;
  }
  // operator== inherited from Point
  
  double len() { return sqrt(x*x + y*y); } // vector length
  double len2() { return (x*x + y*y); } // vector length squared (faster)
  // Special Operations
  double normalize();                 // convert vector to unit length
  friend Vector2D sum( int, int[], Vector2D[]);     // vector sum
  friend Vector2D sum( int, double[], Vector2D[]);  // vector sum
};

#endif
