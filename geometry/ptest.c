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
// ptest - test Point class methods
//------------------------------------------------------------------

main()
{
	double	a, c[3];
	Point	P1(1), P2(0,1);
	Point	Q[3], R;
	Vector	v, w;

	// output preset coords of P1 and P2
	cout << "Initial P1 = " << P1 << endl;
	cout << "Initial P2 = " << P2 << endl;
	cout << endl;

	// input new values for Points P1 and P2
	// to test ">>" and "<<" Point operators
	cout << "Input new P1 = '(x,y,z)': ";
	cin  >> P1;
	cout << "P1 = " << P1 << endl;

	cout << "Input new P2 = '(x,y,z)': ";
	cin  >> P2;
	cout << "P2 = " << P2 << endl;
	cout << endl;

	// test all operators
	v = P1 - P2;		cout << "P1 - P2 = " << v << endl;
	R = 2*P1;		cout << "2*P1 = " << R << endl;
	R = P1/2;		cout << "P1/2 = " << R << endl;
	R = P1 + P2;		cout << "P1 + P2 = " << R << endl;
	R = 2*P1 + P2;		cout << "2*P1 + P2 = " << R << endl;
	R = (P1 + P2)/2;	cout << "(P1 + P2)/2 = " << R << endl;
	R = (P1 + P2)*0.5;	cout << "(P1 + P2)*0.5 = " << R << endl;
	cout << endl;

	// test all functions
	a = d(P1,P2);		cout << "d(P1,P2) = " << a << endl;
	a = d2(P1,P2);		cout << "d2(P1,P2) = " << a << endl;

	c[0]= 0.5; c[1]= 0.5;	cout << "c[]= {" << c[0] << ", " << c[1] << "}" << endl;
	Q[0]= P1;  Q[1]= P2;	cout << "Q[]= {" << Q[0] << ", " << Q[1] << "}" << endl;
	R = asum(2,c,Q);		cout << "asum(2,c,Q) = " << R << endl;
}
