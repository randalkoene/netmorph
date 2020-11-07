//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

// Heavily modified by Randal A. Koene, 20050111

#include "Point3D.hh"
#include "Vector3D.hh"

// Read input Point3D format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
istream& operator>>( istream& input, Point3D& P) {
	char c;
	input >> c;                // skip '('
	input >> P.x;
	input >> c;                
	if (c == ')') {
		P.y=0.0; P.z=0.0;  // 1D coord
		return input;
	}
	input >> P.y;
	input >> c;
	if (c == ')') {
		P.z=0.0;           // 2D coord
		return input;
	}
	input >> P.z;
	input >> c;                // skip ')'
	return input;
}

// Write output Point3D in format: "(%f)", "(%f, %f)", or "(%f, %f, %f)"
ostream& operator<<( ostream& output, Point3D P) {
  output << "(" << P.x << ", " << P.y << ", " << P.z << ")";
  return output;
}

Vector3D Point3D::operator-( Point3D Q)        // Vector3D diff of Point3Ds
{
	Vector3D v;
	v.x = x - Q.x;
	v.y = y - Q.y;
	v.z = z - Q.z;
	return v;
}

Point3D Point3D::operator+( Vector3D v)        // +ve translation
{
	Point3D P;
	P.x = x + v.x;
	P.y = y + v.y;
	P.z = z + v.z;
	return P;
}

Point3D Point3D::operator-( Vector3D v)        // -ve translation
{
	Point3D P;
	P.x = x - v.x;
	P.y = y - v.y;
	P.z = z - v.z;
	return P;
}

Point3D& Point3D::operator+=( Vector3D v)        // +ve translation
{
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

Point3D& Point3D::operator-=( Vector3D v)        // -ve translation
{
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

//------------------------------------------------------------------
// Point3D Scalar Operations (convenient but often illegal)
//        are not valid for points in general,
//        unless they are 'affine' as coeffs of 
//        a sum in which all the coeffs add to 1,
//        such as: the sum (a*P + b*Q) with (a+b == 1).
//        The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point3D operator*( int c, Point3D Q) {
	Point3D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	return P;
}

Point3D operator*( double c, Point3D Q) {
	Point3D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	return P;
}

Point3D operator*( Point3D Q, int c) {
	Point3D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	return P;
}

Point3D operator*( Point3D Q, double c) {
	Point3D P;
	P.x = c * Q.x;
	P.y = c * Q.y;
	P.z = c * Q.z;
	return P;
}

Point3D operator/( Point3D Q, int c) {
	Point3D P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	P.z = Q.z / c;
	return P;
}

Point3D operator/( Point3D Q, double c) {
	Point3D P;
	P.x = Q.x / c;
	P.y = Q.y / c;
	P.z = Q.z / c;
	return P;
}

//------------------------------------------------------------------
// Point3D Addition (also convenient but often illegal)
//    is not valid unless part of an affine sum.
//    The programmer must enforce this (if they want to).
//------------------------------------------------------------------

Point3D operator+( Point3D Q, Point3D R)
{
	Point3D P;
	P.x = Q.x + R.x;
	P.y = Q.y + R.y;
	P.z = Q.z + R.z;
	return P;
}

//------------------------------------------------------------------
// Affine Sums
// Returns weighted sum, even when not affine, but...
// Tests if coeffs add to 1.  If not, sets: err = Esum.
//------------------------------------------------------------------

Point3D asum( int n, int c[], Point3D Q[], int * err) {
  int        cs = 0;
  Point3D      P;
  
  if (err) {
    *err = 0;
    for (int i=0; i<n; i++) cs += c[i];
    if (cs != 1)     // not an affine sum
      *err = -1;     // flag error, but compute sum anyway
  }

  for (int i=0; i<n; i++) {
    P.x += c[i] * Q[i].x;
    P.y += c[i] * Q[i].y;
    P.z += c[i] * Q[i].z;
  }
  return P;
}

Point3D asum( int n, double c[], Point3D Q[], int * err) {
  double     cs = 0.0;
  Point3D      P;
  
  if (err) {
    *err = 0;
    for (int i=0; i<n; i++) cs += c[i];
    if (cs != 1)      // not an affine sum
      *err = -1;      // flag error, but compute sum anyway
  }    
  
  for (int i=0; i<n; i++) {
    P.x += c[i] * Q[i].x;
    P.y += c[i] * Q[i].y;
    P.z += c[i] * Q[i].z;
  }
  return P;
}

//------------------------------------------------------------------
// Distance between Point3Ds
//------------------------------------------------------------------

double Point3D::distance( Point3D & P ) {      // Euclidean distance
  double dx = P.x - x;
  double dy = P.y - y;
  double dz = P.z - z;
  return sqrt(dx*dx + dy*dy + dz*dz); 
}

double distance( Point3D P, Point3D Q) {      // Euclidean distance
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	double dz = P.z - Q.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double distance2( Point3D P, Point3D Q) {     // squared distance (more efficient)
	double dx = P.x - Q.x;
	double dy = P.y - Q.y;
	double dz = P.z - Q.z;
	return (dx*dx + dy*dy + dz*dz);
}

//------------------------------------------------------------------
// Sidedness of a Point3D wrt a directed line P1->P2
//        - makes sense in 2D only
//------------------------------------------------------------------

double Point3D::isLeft( Point3D P1, Point3D P2, int * err) {
	if (err) {
	  if (z != 0.0 || P1.Z() != 0.0 || P2.Z() != 0.0) *err = -1; // flag error, but compute anyway
	  else *err = 0;
	}
	return ((P1.x - x) * (P2.y - y) - (P2.x - x) * (P1.y - y));
}
