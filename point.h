//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#ifndef SS_Point_H
#define SS_Point_H

#include "common.h"

//==================================================================
//  Point Class Definition
//==================================================================

class Point {
  friend class Vector;
 protected:
  double x, y, z;      // z=0 for 2D, y=z=0 for 1D
 public:
  Point() { x=y=z=0; err=Enot; }
	// 1D Point
	Point( int a) {
		dimn=1; x=a; y=z=0; err=Enot; }
	Point( double a) {
		dimn=1; x=a; y=z=0; err=Enot; }
	// 2D Point
	Point( int a, int b) {
		dimn=2; x=a; y=b; z=0; err=Enot; }
	Point( double a, double b) {
		dimn=2; x=a; y=b; z=0; err=Enot; }
	// 3D Point
	Point( int a, int b, int c) {
		dimn=3; x=a; y=b; z=c; err=Enot; }
	Point( double a, double b, double c) {
		dimn=3; x=a; y=b; z=c; err=Enot; }
	// n-dim Point
	Point( int n, int a[]);
	Point( int n, double a[]);
	// Destructor
	~Point() {};

	//----------------------------------------------------------
	// Input/Output streams
	friend istream& operator>>( istream&, Point&);
	friend ostream& operator<<( ostream&, Point);

	//----------------------------------------------------------
	// Assignment "=": use the default to copy all members
	int dim() { return dimn; }      // get dimension
	int setdim( int);               // set new dimension

	//----------------------------------------------------------
	// Comparison (dimension must match, or not)
	int operator==( Point);
	int operator!=( Point);

	//----------------------------------------------------------
	// Point and Vector Operations (always valid) 
	Vector operator-( Point);       // Vector difference
	Point  operator+( Vector);      // +translate
	Point  operator-( Vector);      // -translate
	Point& operator+=( Vector);     // inc translate
	Point& operator-=( Vector);     // dec translate

	//----------------------------------------------------------
	// Point Scalar Operations (convenient but often illegal)
	// using any type of scalar (int, float, or double)
	//    are not valid for points in general,
	//    unless they are 'affine' as coeffs of 
	//    a sum in which all the coeffs add to 1,
	//    such as: the sum (a*P + b*Q) with (a+b == 1).
	//    The programmer must enforce this (if they want to).

	// Scalar Multiplication
	friend Point operator*( int, Point);
	friend Point operator*( double, Point);
	friend Point operator*( Point, int);
	friend Point operator*( Point, double);
	// Scalar Division
	friend Point operator/( Point, int);
	friend Point operator/( Point, double);

	//----------------------------------------------------------
	// Point Addition (also convenient but often illegal)
	//    is not valid unless part of an affine sum.
	//    The programmer must enforce this (if they want to).
	friend Point operator+( Point, Point);     // add points

	// Affine Sum
	// Returns weighted sum, even when not affine, but...
	// Tests if coeffs add to 1.  If not, sets: err = Esum.
	friend Point asum( int, int[], Point[]);
	friend Point asum( int, double[], Point[]);

	//----------------------------------------------------------
	// Point Relations
	friend double d( Point, Point);         // Distance
	friend double d2( Point, Point);        // Distance^2
	double isLeft( Point, Point);           // 2D only
	double Area( Point, Point); 		// any dim for triangle PPP

	// Collinearity Conditions (any dim n)
	boolean isOnLine( Point, Point, char);  // is On line (char= flag)
	boolean isOnLine( Point, Point);        // is On line (flag= all)
	boolean isBefore( Point, Point);        // is On line (flag= before)
	boolean isBetween( Point, Point);       // is On line (flag= between)
	boolean isAfter( Point, Point);         // is On line (flag= after)
	boolean isOnRay( Point, Point);         // is On line (flag= between|after)

	//----------------------------------------------------------
	// Error Handling
	void  clerr() { err = Enot;}            // clear error
	int   geterr() { return err;}           // get error
	char* errstr();                         // return error string
};
#endif SS_Point_H
