//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#ifndef SS_Point2D_H
#define SS_Point2D_H

#include <iostream>
using namespace std;

#ifndef SS_Vector2D_H
class Vector2D;
#endif

// utility macros
#define	Pabs(x)		((x) >= 0 ? x : -(x))
#define	Pmin(x,y)	((x) < (y) ? (x) : (y))
#define	Pmax(x,y)	((x) > (y) ? (x) : (y))

#define PDIMSIZE 2

class Point2D {
  friend class Vector2D;
 protected:
  double x, y;      // y=0 for 1D
 public:
  Point2D(): x(0.0), y(0.0) {}
  Point2D(double X, double Y): x(X), y(Y) {}
  ~Point2D() {};

  double X() const { return x; }
  double Y() const { return y; }
  void set_X(double X) { x = X; }
  void set_Y(double Y) { y = Y; }
  void add_X(double X) { x += X; }
  void add_Y(double Y) { y += Y; }

  double & operator [] (unsigned int n) {
    if (n==0) return x;
    else return y;
  }

  friend istream& operator>>( istream & input, Point2D & P);
  friend ostream& operator<<( ostream & output, Point2D P);

  int operator==( Point2D Q) { return (x==Q.x && y==Q.y); }
  int operator!=( Point2D Q) { return (x!=Q.x || y!=Q.y); }

  Vector2D operator-( Point2D);       // Vector2D difference
  Point2D  operator+( Vector2D);      // +translate
  Point2D  operator-( Vector2D);      // -translate
  Point2D& operator+=( Vector2D);     // inc translate
  Point2D& operator-=( Vector2D);     // dec translate
  
  //----------------------------------------------------------
  // Point2D Scalar Operations (convenient but often illegal)
  // using any type of scalar (int, float, or double)
  //    are not valid for points in general,
  //    unless they are 'affine' as coeffs of 
  //    a sum in which all the coeffs add to 1,
  //    such as: the sum (a*P + b*Q) with (a+b == 1).
  //    The programmer must enforce this (if they want to).
  
  // Scalar Multiplication
  friend Point2D operator*( int, Point2D);
  friend Point2D operator*( double, Point2D);
  friend Point2D operator*( Point2D, int);
  friend Point2D operator*( Point2D, double);
  // Scalar Division
  friend Point2D operator/( Point2D, int);
  friend Point2D operator/( Point2D, double);

  //----------------------------------------------------------
  // Point2D Addition (also convenient but often illegal)
  //    is not valid unless part of an affine sum.
  //    The programmer must enforce this (if they want to).
  friend Point2D operator+( Point2D, Point2D);     // add points

  // Affine Sum
  // Returns weighted sum, even when not affine, but...
  // Tests if coeffs add to 1.  If not, sets: err = Esum.
  friend Point2D asum( int, int[], Point2D[], int * err = 0);
  friend Point2D asum( int, double[], Point2D[], int * err = 0);

  //----------------------------------------------------------
  // Point2D Relations
  friend double distance( Point2D, Point2D);         // Distance
  friend double distance2( Point2D, Point2D);        // Distance^2
  double distance(Point2D & P); // Distance
  double isLeft( Point2D, Point2D, int * err = 0); // 2D only
  double Area( Point2D, Point2D); 		// any dim for triangle PPP
  
  // Collinearity Conditions (any dim n)
  bool isOnLine( Point2D, Point2D, char);  // is On line (char= flag)
  bool isOnLine( Point2D, Point2D);        // is On line (flag= all)
  bool isBefore( Point2D, Point2D);        // is On line (flag= before)
  bool isBetween( Point2D, Point2D);       // is On line (flag= between)
  bool isAfter( Point2D, Point2D);         // is On line (flag= after)
  bool isOnRay( Point2D, Point2D);         // is On line (flag= between|after)

  bool isZero() { return ((x==0.0) && (y==0.0)); }
	
};

#endif
