//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#ifndef SS_Point3D_H
#define SS_Point3D_H

#include <iostream>
using namespace std;

#ifndef SS_Vector3D_H
class Vector3D;
#endif

// utility macros
#define	Pabs(x)		((x) >= 0 ? x : -(x))
#define	Pmin(x,y)	((x) < (y) ? (x) : (y))
#define	Pmax(x,y)	((x) > (y) ? (x) : (y))

// Now:
// (2) Update the Makefile to get compiler directives in line with those
//     that I use.
// (3) Make the basic point classes very lean, add derived class in a
//     separate object file that has all the operators a type usually has.
// (4) Possibly turn these into templates, so that I can choose different
//     types other than double as the member variables and thereby alter
//     precision, memory consumption and more.

#define PDIMSIZE 3

class Point3D {
  friend class Vector3D;
 protected:
  double x, y, z;      // z=0 for 2D, y=z=0 for 1D
 public:
  Point3D(): x(0.0), y(0.0), z(0.0) {}
  Point3D(double X, double Y, double Z = 0.0): x(X), y(Y), z(Z) {}
  ~Point3D() {};

  double X() const { return x; }
  double Y() const { return y; }
  double Z() const { return z; }
  void set_X(double X) { x = X; }
  void set_Y(double Y) { y = Y; }
  void set_Z(double Z) { z = Z; }
  void add_X(double X) { x += X; }
  void add_Y(double Y) { y += Y; }
  void add_Z(double Z) { z += Z; }

  double & operator [] (unsigned int n) {
    if (n==0) return x;
    else if (n==1) return y;
    else return z;
  }

  friend istream& operator>>( istream & input, Point3D & P);
  friend ostream& operator<<( ostream & output, Point3D P);

  int operator==( Point3D Q) { return (x==Q.x && y==Q.y && z==Q.z); }
  int operator!=( Point3D Q) { return (x!=Q.x || y!=Q.y || z!=Q.z); }

  Vector3D operator-( Point3D);       // Vector3D difference
  Point3D  operator+( Vector3D);      // +translate
  Point3D  operator-( Vector3D);      // -translate
  Point3D& operator+=( Vector3D);     // inc translate
  Point3D& operator-=( Vector3D);     // dec translate
  
  //----------------------------------------------------------
  // Point3D Scalar Operations (convenient but often illegal)
  // using any type of scalar (int, float, or double)
  //    are not valid for points in general,
  //    unless they are 'affine' as coeffs of 
  //    a sum in which all the coeffs add to 1,
  //    such as: the sum (a*P + b*Q) with (a+b == 1).
  //    The programmer must enforce this (if they want to).
  
  // Scalar Multiplication
  friend Point3D operator*( int, Point3D);
  friend Point3D operator*( double, Point3D);
  friend Point3D operator*( Point3D, int);
  friend Point3D operator*( Point3D, double);
  // Scalar Division
  friend Point3D operator/( Point3D, int);
  friend Point3D operator/( Point3D, double);

  //----------------------------------------------------------
  // Point3D Addition (also convenient but often illegal)
  //    is not valid unless part of an affine sum.
  //    The programmer must enforce this (if they want to).
  friend Point3D operator+( Point3D, Point3D);     // add points

  // Affine Sum
  // Returns weighted sum, even when not affine, but...
  // Tests if coeffs add to 1.  If not, sets: err = Esum.
  friend Point3D asum( int, int[], Point3D[], int * err = 0);
  friend Point3D asum( int, double[], Point3D[], int * err = 0);

  //----------------------------------------------------------
  // Point3D Relations
  friend double distance( Point3D, Point3D);         // Distance
  friend double distance2( Point3D, Point3D);        // Distance^2
  double distance(Point3D & P); // Distance
  double isLeft( Point3D, Point3D, int * err = 0); // 2D only
  double Area( Point3D, Point3D); 		// any dim for triangle PPP
  
  // Collinearity Conditions (any dim n)
  bool isOnLine( Point3D, Point3D, char);  // is On line (char= flag)
  bool isOnLine( Point3D, Point3D);        // is On line (flag= all)
  bool isBefore( Point3D, Point3D);        // is On line (flag= before)
  bool isBetween( Point3D, Point3D);       // is On line (flag= between)
  bool isAfter( Point3D, Point3D);         // is On line (flag= after)
  bool isOnRay( Point3D, Point3D);         // is On line (flag= between|after)
	
  bool isZero() { return ((x==0.0) && (y==0.0) && (z==0.0)); }
	
};

#endif
