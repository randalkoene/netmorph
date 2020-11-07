//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#include "Point2D.hh"
#include "Vector2D.hh"

// Read input Point2D format: "(%f)" or "(%f, %f)"
istream& operator>>( istream& input, Point2D& P) {
	char c;
	input >> c;                // skip '('
	input >> P.x;
	input >> c;                
	if (c == ')') {
		P.y=0.0;           // 1D coord
		return input;
	}
	input >> P.y;
	input >> c;                // skip ')'
	return input;
}

// Write output Point2D in format: "(%f, %f)"
ostream& operator<<( ostream& output, Point2D P) {
  output << "(" << P.x << ", " << P.y << ")";
  return output;
}

Vector2D Point2D::operator-( Point2D Q)        // Vector2D diff of Point2Ds
{
	Vector2D v;
	v.x = x - Q.x;
	v.y = y - Q.y;
	return v;
}

Point2D Point2D::operator+( Vector2D v)        // +ve translation
{
	Point2D P;
	P.x = x + v.x;
	P.y = y + v.y;
	return P;
}

Point2D Point2D::operator-( Vector2D v)        // -ve translation
{
	Point2D P;
	P.x = x - v.x;
	P.y = y - v.y;
	return P;
}

Point2D& Point2D::operator+=( Vector2D v)        // +ve translation
{
	x += v.x;
	y += v.y;
	return *this;
}

Point2D& Point2D::operator-=( Vector2D v)        // -ve translation
{
	x -= v.x;
	y -= v.y;
	return *this;
}

//------------------------------------------------------------------
// Point2D Scalar Operations (convenient but often illegal)
//        are not valid for points in general,
//        unless they are 'affine' as coeffs of 
//        a sum in which all the coeffs add to 1,
//        such as: the sum (a*P + b*Q) with (a+b == 1).
//        The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point2D operator*( int c, Point2D Q) {
	Point2D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	return P;
}

Point2D operator*( double c, Point2D Q) {
	Point2D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	return P;
}

Point2D operator*( Point2D Q, int c) {
	Point2D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	return P;
}

Point2D operator*( Point2D Q, double c) {
	Point2D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	return P;
}

Point2D operator/( Point2D Q, int c) {
	Point2D P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	return P;
}

Point2D operator/( Point2D Q, double c) {
	Point2D P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	return P;
}

//------------------------------------------------------------------
// Point2D Addition (also convenient but often illegal)
//    is not valid unless part of an affine sum.
//    The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point2D operator+( Point2D Q, Point2D R)
{
	Point2D P;
	P.x = Q.x + R.x;
	P.y = Q.y + R.y;
	return P;
}

//------------------------------------------------------------------
// Affine Sums
// Returns weighted sum, even when not affine, but...
// Tests if coeffs add to 1.  If not, sets: err = Esum.
//------------------------------------------------------------------

Point2D asum( int n, int c[], Point2D Q[], int * err) {
  int        cs = 0;
  Point2D      P;
  
  if (err) {
    *err = 0;
    for (int i=0; i<n; i++) cs += c[i];
    if (cs != 1)     // not an affine sum
      *err = -1;     // flag error, but compute sum anyway
  }

  for (int i=0; i<n; i++) {
    P.x += c[i] * Q[i].x;
    P.y += c[i] * Q[i].y;
  }
  return P;
}

Point2D asum( int n, double c[], Point2D Q[], int * err) {
  double     cs = 0.0;
  Point2D      P;

  if (err) {
    *err = 0;
    for (int i=0; i<n; i++) cs += c[i];
    if (cs != 1)      // not an affine sum
      *err = -1;      // flag error, but compute sum anyway
  }
  
  for (int i=0; i<n; i++) {
    P.x += c[i] * Q[i].x;
    P.y += c[i] * Q[i].y;
  }
  return P;
}

//------------------------------------------------------------------
// Distance between Point2Ds
//------------------------------------------------------------------

double Point2D::distance( Point2D & P ) {      // Euclidean distance
  double dx = P.x - x;
  double dy = P.y - y;
  return sqrt(dx*dx + dy*dy);
}

double distance( Point2D P, Point2D Q) {      // Euclidean distance
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	return sqrt(dx*dx + dy*dy);
}

double distance2( Point2D P, Point2D Q) {     // squared distance (more efficient)
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	return (dx*dx + dy*dy);
}

//------------------------------------------------------------------
// Sidedness of a Point2D wrt a directed line P1->P2
//        - makes sense in 2D only
//------------------------------------------------------------------

double Point2D::isLeft( Point2D P1, Point2D P2, int * err) {
  if (err) *err = 0;
  return ((P1.x - x) * (P2.y - y) - (P2.x - x) * (P1.y - y));
}
