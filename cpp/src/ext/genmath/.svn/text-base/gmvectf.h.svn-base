/*-------------------------------------------------------------------------
*
* File name:      gmvectf.h
*
* Project:        RADIA, ...
*
* Description:    Definition of simple vector(3x1) and matrix(3x3) structures
*
* Author(s):      Oleg Chubar
*
* First release:  1997
* 
* Copyright (C):  1997 by European Synchrotron Radiation Facility, France
*                 All Rights Reserved
*
-------------------------------------------------------------------------*/

#ifndef __GMVECTF_H
#define __GMVECTF_H

#include "gmvect.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TMatrix3df;

//-------------------------------------------------------------------------

struct TVector3df {
	float x,y,z;

	TVector3df(float xx =0, float yy =0, float zz =0) { x=xx; y=yy; z=zz;}
	TVector3df(float* dArray) { x=dArray[0]; y=dArray[1]; z=dArray[2];}

	TVector3df& operator +=(const TVector3df& AnotherVect)
	{
		x+=AnotherVect.x; y+=AnotherVect.y; z+=AnotherVect.z; return *this;
	}

	TVector3df& operator =(const TVector3d& V)
	{
		x = float(V.x); y = float(V.y); z = float(V.z); return *this;
	}

	friend TVector3df operator +(const TVector3df&, const TVector3df&);
	friend TVector3df operator -(const TVector3df&, const TVector3df&);
	friend TVector3df operator *(const float, const TVector3df&);
	friend float operator *(const TVector3df&, const TVector3df&); // Scalar product
	inline friend TVector3df operator ^(const TVector3df&, const TVector3df&); // Vector product
	friend TVector3df operator *(const TMatrix3df&, const TVector3df&);
};

//-------------------------------------------------------------------------

inline TVector3df operator +(const TVector3df& P1, const TVector3df& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3df(P1.x+P2.x, P1.y+P2.y, P1.z+P2.z);
}

//-------------------------------------------------------------------------

inline TVector3df operator -(const TVector3df& P1, const TVector3df& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3df(P1.x-P2.x, P1.y-P2.y, P1.z-P2.z);
}

//-------------------------------------------------------------------------

inline TVector3df operator *(const float D, const TVector3df& P)
{
	// The following can cause problems with Code Warrior
	return TVector3df(D*P.x, D*P.y, D*P.z);
}

//-------------------------------------------------------------------------

inline float operator *(const TVector3df& P1, const TVector3df& P2)
{
	return P1.x*P2.x+P1.y*P2.y+P1.z*P2.z;
}

//-------------------------------------------------------------------------

inline TVector3df operator ^(const TVector3df& v1, const TVector3df& v2)
{
	// The following can cause problems with Code Warrior
	return TVector3df(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TMatrix3df {
	TVector3df Str0, Str1, Str2;

	TMatrix3df(const TVector3df& InpStr0, const TVector3df& InpStr1, const TVector3df& InpStr2)
	{
		Str0=InpStr0; Str1=InpStr1; Str2=InpStr2;
	}
	TMatrix3df() {}

	TMatrix3df& operator +=(const TMatrix3df& AnotherMatrix)
	{
		Str0+=AnotherMatrix.Str0; Str1+=AnotherMatrix.Str1; Str2+=AnotherMatrix.Str2; return *this;
	}

	TMatrix3df& operator =(const TMatrix3d& M)
	{
		Str0 = M.Str0; Str1 = M.Str1; Str2 = M.Str2; return *this;
	}

	friend float detMatrix3d(const TMatrix3df&);
	friend TMatrix3df Matrix3d_inv(const TMatrix3df&);
	friend void Matrix3d_inv(const TMatrix3df&, TMatrix3df&);

	friend TMatrix3df operator +(const TMatrix3df&, const TMatrix3df&);
	friend TMatrix3df operator -(const TMatrix3df&, const TMatrix3df&);
	friend TMatrix3df operator *(float d, const TMatrix3df&);
	friend TVector3df operator *(const TMatrix3df&, const TVector3df&);
	friend TVector3d operator *(const TMatrix3df&, const TVector3d&);
	friend TMatrix3df operator *(const TMatrix3df&, const TMatrix3df&);
	friend TMatrix3d operator *(const TMatrix3d&, const TMatrix3df&);
};

//-------------------------------------------------------------------------

inline float detMatrix3d(const TMatrix3df& M)
{
	return M.Str0.x*M.Str1.y*M.Str2.z + M.Str0.y*M.Str1.z*M.Str2.x + M.Str0.z*M.Str1.x*M.Str2.y
		  -M.Str0.z*M.Str1.y*M.Str2.x - M.Str0.x*M.Str1.z*M.Str2.y - M.Str0.y*M.Str1.x*M.Str2.z;
}

//-------------------------------------------------------------------------

inline TMatrix3df Matrix3d_inv(const TMatrix3df& M)
{
	TVector3df St0( M.Str1.y*M.Str2.z-M.Str1.z*M.Str2.y,-M.Str0.y*M.Str2.z+M.Str0.z*M.Str2.y, M.Str0.y*M.Str1.z-M.Str0.z*M.Str1.y);
	TVector3df St1(-M.Str1.x*M.Str2.z+M.Str1.z*M.Str2.x, M.Str0.x*M.Str2.z-M.Str0.z*M.Str2.x,-M.Str0.x*M.Str1.z+M.Str0.z*M.Str1.x);
	TVector3df St2( M.Str1.x*M.Str2.y-M.Str1.y*M.Str2.x,-M.Str0.x*M.Str2.y+M.Str0.y*M.Str2.x, M.Str0.x*M.Str1.y-M.Str0.y*M.Str1.x);
	float invDet = float(1./detMatrix3d(M));
	// The following can cause problems with Code Warrior
	return TMatrix3df(invDet*St0, invDet*St1, invDet*St2);
}

//-------------------------------------------------------------------------

inline void Matrix3d_inv(const TMatrix3df& M, TMatrix3df& OutInvMatr)
{
	if(M.Str0.y==0. && (M.Str0.z==0. && (M.Str1.x==0. && (M.Str1.z==0. && (M.Str2.x==0. && M.Str2.y==0.)))))
	{
		OutInvMatr = M;
		float d0 = M.Str0.x;
		float d1 = M.Str1.y;
		float d2 = M.Str2.z;
		if(d0!=0.) OutInvMatr.Str0.x = float(1./d0);
		if(d1!=0.) OutInvMatr.Str1.y = float(1./d1);
		if(d2!=0.) OutInvMatr.Str2.z = float(1./d2);
		return;
	}

	const TVector3df& InSt0 = M.Str0;
	const TVector3df& InSt1 = M.Str1;
	const TVector3df& InSt2 = M.Str2;
	TVector3df St0( InSt1.y*InSt2.z-InSt1.z*InSt2.y,-InSt0.y*InSt2.z+InSt0.z*InSt2.y, InSt0.y*InSt1.z-InSt0.z*InSt1.y);
	TVector3df St1(-InSt1.x*InSt2.z+InSt1.z*InSt2.x, InSt0.x*InSt2.z-InSt0.z*InSt2.x,-InSt0.x*InSt1.z+InSt0.z*InSt1.x);
	TVector3df St2( InSt1.x*InSt2.y-InSt1.y*InSt2.x,-InSt0.x*InSt2.y+InSt0.y*InSt2.x, InSt0.x*InSt1.y-InSt0.y*InSt1.x);
	float invDet = float(1./detMatrix3d(M));	// Do something if Det=0 ?

	OutInvMatr.Str0 = invDet*St0;
	OutInvMatr.Str1 = invDet*St1;
	OutInvMatr.Str2 = invDet*St2;
}
//-------------------------------------------------------------------------

inline TMatrix3df operator +(const TMatrix3df& M1, const TMatrix3df& M2)
{
	// The following can cause problems with Code Warrior
	return TMatrix3df(M1.Str0+M2.Str0, M1.Str1+M2.Str1, M1.Str2+M2.Str2);
}

//-------------------------------------------------------------------------

inline TMatrix3df operator -(const TMatrix3df& M1, const TMatrix3df& M2)
{
	// The following can cause problems with Code Warrior
	return TMatrix3df(M1.Str0-M2.Str0, M1.Str1-M2.Str1, M1.Str2-M2.Str2);
}

//-------------------------------------------------------------------------

inline TMatrix3df operator *(float d, const TMatrix3df& M)
{
	// The following can cause problems with Code Warrior
	return TMatrix3df(d*M.Str0, d*M.Str1, d*M.Str2);
}

//-------------------------------------------------------------------------

inline TVector3df operator *(const TMatrix3df& M, const TVector3df& P)
{
	// The following can cause problems with Code Warrior
	return TVector3df(M.Str0*P, M.Str1*P, M.Str2*P);
}

//-------------------------------------------------------------------------

inline TVector3d operator *(const TMatrix3df& M, const TVector3d& InP)
{
	TVector3df P;
	P = InP;
	// The following can cause problems with Code Warrior
	return TVector3d(M.Str0*P, M.Str1*P, M.Str2*P);
}

//-------------------------------------------------------------------------

inline TMatrix3df operator *(const TMatrix3df& M1, const TMatrix3df& M2)
{
	const TVector3df& M1Str0 = M1.Str0;
	const TVector3df& M1Str1 = M1.Str1;
	const TVector3df& M1Str2 = M1.Str2;
	const TVector3df& M2Str0 = M2.Str0;
	const TVector3df& M2Str1 = M2.Str1;
	const TVector3df& M2Str2 = M2.Str2;
	TVector3df St0(M1Str0.x*M2Str0.x+M1Str0.y*M2Str1.x+M1Str0.z*M2Str2.x, M1Str0.x*M2Str0.y+M1Str0.y*M2Str1.y+M1Str0.z*M2Str2.y, M1Str0.x*M2Str0.z+M1Str0.y*M2Str1.z+M1Str0.z*M2Str2.z);
	TVector3df St1(M1Str1.x*M2Str0.x+M1Str1.y*M2Str1.x+M1Str1.z*M2Str2.x, M1Str1.x*M2Str0.y+M1Str1.y*M2Str1.y+M1Str1.z*M2Str2.y, M1Str1.x*M2Str0.z+M1Str1.y*M2Str1.z+M1Str1.z*M2Str2.z);
	TVector3df St2(M1Str2.x*M2Str0.x+M1Str2.y*M2Str1.x+M1Str2.z*M2Str2.x, M1Str2.x*M2Str0.y+M1Str2.y*M2Str1.y+M1Str2.z*M2Str2.y, M1Str2.x*M2Str0.z+M1Str2.y*M2Str1.z+M1Str2.z*M2Str2.z);
	// The following can cause problems with Code Warrior
	return TMatrix3df(St0, St1, St2);
}

//-------------------------------------------------------------------------

inline TMatrix3d operator *(const TMatrix3d& M1, const TMatrix3df& M2)
{
	const TVector3d& M1Str0 = M1.Str0;
	const TVector3d& M1Str1 = M1.Str1;
	const TVector3d& M1Str2 = M1.Str2;
	const TVector3df& M2Str0 = M2.Str0;
	const TVector3df& M2Str1 = M2.Str1;
	const TVector3df& M2Str2 = M2.Str2;
	TVector3d St0(M1Str0.x*M2Str0.x+M1Str0.y*M2Str1.x+M1Str0.z*M2Str2.x, M1Str0.x*M2Str0.y+M1Str0.y*M2Str1.y+M1Str0.z*M2Str2.y, M1Str0.x*M2Str0.z+M1Str0.y*M2Str1.z+M1Str0.z*M2Str2.z);
	TVector3d St1(M1Str1.x*M2Str0.x+M1Str1.y*M2Str1.x+M1Str1.z*M2Str2.x, M1Str1.x*M2Str0.y+M1Str1.y*M2Str1.y+M1Str1.z*M2Str2.y, M1Str1.x*M2Str0.z+M1Str1.y*M2Str1.z+M1Str1.z*M2Str2.z);
	TVector3d St2(M1Str2.x*M2Str0.x+M1Str2.y*M2Str1.x+M1Str2.z*M2Str2.x, M1Str2.x*M2Str0.y+M1Str2.y*M2Str1.y+M1Str2.z*M2Str2.y, M1Str2.x*M2Str0.z+M1Str2.y*M2Str1.z+M1Str2.z*M2Str2.z);
	// The following can cause problems with Code Warrior
	return TMatrix3d(St0, St1, St2);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

inline TVector3d::TVector3d(const TVector3df& Vf)
{
	x = Vf.x; y = Vf.y; z = Vf.z;
}

//-------------------------------------------------------------------------

inline TMatrix3d::TMatrix3d(const TMatrix3df& Mf)
{
	Str0 = TVector3d(Mf.Str0);
	Str1 = TVector3d(Mf.Str1);
	Str2 = TVector3d(Mf.Str2);
}

//-------------------------------------------------------------------------

#endif
