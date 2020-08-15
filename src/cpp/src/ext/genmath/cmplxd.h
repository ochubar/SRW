/************************************************************************//**
 * File: cmplxd.h
 * Description: Definition of the simplest complex<double> structure
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __CMPLXD_H
#define __CMPLXD_H

//#include <math.h>

//-------------------------------------------------------------------------

struct TComplexD {
	double x, y;

	TComplexD(double xx =0, double yy =0) { x=xx; y=yy;}
	TComplexD(double* dArray) { x=dArray[0]; y=dArray[1];}

	double AbsE2() { return (x*x + y*y);}
	double Abs() { return sqrt(x*x + y*y);}
	void Conjugate() { y = -y;}

	TComplexD& operator +=(const TComplexD& c)
	{
		x += c.x; y += c.y; return *this;
	}
	TComplexD& operator *=(const TComplexD& c)
	{
		double ax = x*c.x - y*c.y;
		double ay = x*c.y + y*c.x;
		x = ax; y = ay; 
		return *this;
	}

	inline friend TComplexD operator +(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator +(const double, const TComplexD&);
	inline friend TComplexD operator +(const TComplexD&, const double);
	inline friend TComplexD operator -(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator -(const double, const TComplexD&);
	inline friend TComplexD operator -(const TComplexD&, const double);
	inline friend TComplexD operator *(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator *(const double, const TComplexD&);
	inline friend TComplexD operator *(const TComplexD&, const double);
	inline friend TComplexD operator /(const TComplexD&, const TComplexD&);
	inline friend TComplexD operator /(const double, const TComplexD&);
	inline friend TComplexD operator /(const TComplexD&, const double);
};

//-------------------------------------------------------------------------

inline TComplexD operator +(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x + P2.x, P1.y + P2.y);
	return Res;
}
inline TComplexD operator +(const double a, const TComplexD& P)
{
	TComplexD Res(a + P.x, P.y);
	return Res;
}
inline TComplexD operator +(const TComplexD& P, const double a)
{
	TComplexD Res(a + P.x, P.y);
	return Res;
}
inline TComplexD operator -(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x - P2.x, P1.y - P2.y);
	return Res;
}
inline TComplexD operator -(const double a, const TComplexD& P)
{
	TComplexD Res(a - P.x, -P.y);
	return Res;
}
inline TComplexD operator -(const TComplexD& P, const double a)
{
	TComplexD Res(P.x - a, P.y);
	return Res;
}
inline TComplexD operator *(const TComplexD& P1, const TComplexD& P2)
{
	TComplexD Res(P1.x*P2.x - P1.y*P2.y, P1.x*P2.y + P1.y*P2.x);
	return Res;
}
inline TComplexD operator *(const double a, const TComplexD& P)
{
	TComplexD Res(a*P.x, a*P.y);
	return Res;
}
inline TComplexD operator *(const TComplexD& P, const double a)
{
	TComplexD Res(a*P.x, a*P.y);
	return Res;
}
inline TComplexD operator /(const TComplexD& P1, const TComplexD& P2)
{
	double InvAbsP2E2 = 1./(P2.x*P2.x + P2.y*P2.y);
	TComplexD Res((P1.x*P2.x + P1.y*P2.y)*InvAbsP2E2, (P1.y*P2.x - P1.x*P2.y)*InvAbsP2E2);
	return Res;
}
inline TComplexD operator /(const double a, const TComplexD& P)
{
	double InvAbsP2E2 = 1./(P.x*P.x + P.y*P.y);
	TComplexD Res(a*P.x*InvAbsP2E2, -a*P.y*InvAbsP2E2);
	return Res;
}
inline TComplexD operator /(const TComplexD& P, const double a)
{
	double Inv_a = 1./a;
	TComplexD Res(Inv_a*P.x, -Inv_a*P.y);
	return Res;
}

//-------------------------------------------------------------------------

#endif
