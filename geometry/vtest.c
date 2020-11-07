//==================================================================
// Copyright 2002, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from it's use.
// Users of this code must verify correctness for their application.
//==================================================================

#include "point.h"
#include "vector.h"

//------------------------------------------------------------------
// vtest - test Vector class methods
//------------------------------------------------------------------

main()
{
	double	a, c[3];
	Vector	u, v(1), w(0,1);
	Vector	vx[3];


	// output preset coords of v and w
	cout << "Initial v = " << v << endl;
	cout << "Initial w = " << w << endl;
	cout << endl;

	// input new values for vectors v and w
	// to test ">>" and "<<" Vector operators
	cout << "Input v = '(x,y,z)': ";
	cin  >> v;
	cout << "v = " << v << endl;
	cout << "Input w = '(x,y,z)': ";
	cin  >> w;
	cout << "w = " << w << endl;
	cout << endl;

	// test operators
	u = -v;			cout << "-v = " << u << endl;
	u = ~v;			cout << "~v = " << u << endl;
	u = 2*v;		cout << "2*v = " << u << endl;
	u = v + w;		cout << "v + w = " << u << endl;
	u = (v + w)/2;		cout << "(v + w)/2 = " << u << endl;
	u = v - w;		cout << "v - w = " << u << endl;
	u = (v - w)*0.5;	cout << "(v - w)*0.5 = " << u << endl;
	u = v ^ w;		cout << "v ^ w = " << u << endl;
	a = v * w;		cout << "v * w = " << a << endl;
	a = v | w;		cout << "v | w = " << a << endl;
	a = ~v * w;		cout << "~v * w = " << a << endl;
	cout << endl;

	// test functions
	c[0]= 0.5; c[1]= 0.5;	cout << "c[]= {" << c[0] << ", " << c[1] << "}\n";
	vx[0]= v;  vx[1]= w;	cout << "vx[]= {" << vx[0] << ", " << vx[1] << "}\n";
	u = sum(2,c,vx);	cout << "sum(2,c,vx) = " << u << endl;

	u = v.len2();		cout << "v.len2() = " << u << endl;
	u = v.len();		cout << "v.len()  = " << u << endl;
	v.normalize();		cout << "v.normalize() = " << v << endl;
	w.normalize();		cout << "w.normalize() = " << w << endl;
	cout << endl;

	// check that a few identities are always true
	int cnt=0, errs=0;
	Vector u0, uinc(2,0,-2);
	Vector v0, vinc(1,1,1);
	Vector w0, winc(0.2, 0.2, 0.0);
	Vector x, x0(-1,-2,-3), xinc(0.5,-0.5,0.1);
	double dl, dr;		// left and right side doubles
	Vector vl, vr;		// left and right side vectors
	double er, e=0.000000000001; // allow small float error e

	for (u=u0; u.x < 7; u+=uinc) 
	for (v=v0; v.x < 11; v+=vinc) 
	for (w=w0; w.x < 1.1; w+=winc) {
		++cnt;
		// Triple product
		dl = u * (v ^ w);
		dr = (u ^ v) * w;
		er = abs(dl-dr);
		if (er > e) {	// identity fails > e
			++errs;	// count errors
			cout << "Error: err= " << er << endl;
		}
		++cnt;
		// Left Association of Cross Product
		vl = (u ^ v) ^ w;
		vr = (u * w) * v - (v * w) * u;
		er = (vl-vr).len();
		if (er > e) {	// identity fails > e
			++errs;	// count errors
			cout << "Error: err= " << er << endl;
		}
		++cnt;
		// Lagrange Identity
		for (x=x0; x.x < 1.1; x+=xinc) {
			dl = (u ^ v) * (w ^ x);
			dr = (u * w) * (v * x) - (v * w) * (u * x);
			er = abs(dl-dr);
			if (er > e) {	// identity fails > e
				++errs;	// count errors
				cout << "Error: err= " << er << endl;
			}
			++cnt;
		}
	}
	cout << "Did " << cnt << " identity chks. Had " << errs << " errors\n";
}
