/************************************************************************//**
 * File: gmvect.h
 * Description: Definition of simple vector(3x1) and matrix(3x3) structures
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMVECT_H
#define __GMVECT_H

#include <math.h>

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TMatrix3d;
struct TMatrix3df;
struct TVector3df;

//-------------------------------------------------------------------------

struct TVector3d {
	double x,y,z;

	TVector3d(double xx =0, double yy =0, double zz =0) { x=xx; y=yy; z=zz;}
	TVector3d(double* dArray) { x=dArray[0]; y=dArray[1]; z=dArray[2];}
	TVector3d(const TVector3df&);
	TVector3d(const char* XorYorZ)
	{
		x = y = z = 0;
		if(XorYorZ != 0)
		{
			if((*XorYorZ == 'x') || (*XorYorZ == 'X')) x = 1.;
			else if((*XorYorZ == 'y') || (*XorYorZ == 'Y')) y = 1.;
			else if((*XorYorZ == 'z') || (*XorYorZ == 'Z')) z = 1.;
		}
	}

	TVector3d& operator +=(const TVector3d& AnotherVect)
	{
		x+=AnotherVect.x; y+=AnotherVect.y; z+=AnotherVect.z; return *this;
	}
	TVector3d& operator -=(const TVector3d& AnotherVect)
	{
		x-=AnotherVect.x; y-=AnotherVect.y; z-=AnotherVect.z; return *this;
	}
	TVector3d& operator *=(double a)
	{
		x *= a; y *= a; z *= a; return *this;
	}
	TVector3d& operator /=(double a)
	{
		x /= a; y /= a; z /= a; return *this;
	}

	inline double Abs() { return sqrt(x*x + y*y + z*z);}
	inline double AmpE2() { return x*x + y*y + z*z;}
	inline void Zero() { x = y = z = 0;}
	inline bool isZero() { return (x == 0) && (y == 0) && (z == 0);}
	inline double absMaxElem()
	{
		double ax = ::fabs(x), ay = ::fabs(y), az = ::fabs(z);
		double res = ax;
		if(res < ay) res = ay;
		if(res < az) res = az;
		return res;
	}
	inline void Normalize()
	{
		if((x == 0) && (y == 0) && (z == 0)) return;
		double invNorm = 1./sqrt(x*x + y*y + z*z);
		x *= invNorm; y *= invNorm; z *= invNorm; 
	}

	inline friend TVector3d operator +(const TVector3d&, const TVector3d&);
	inline friend TVector3d operator -(const TVector3d&, const TVector3d&);
	inline friend TVector3d operator *(const double, const TVector3d&);
	inline friend TVector3d operator *(const TVector3d&, const double);
	inline friend TVector3d operator /(const TVector3d&, const double);
	inline friend double operator *(const TVector3d&, const TVector3d&); // Scalar product
	inline friend TVector3d operator ^(const TVector3d&, const TVector3d&); // Vector product
	inline friend TVector3d operator *(const TMatrix3d&, const TVector3d&);
	friend TVector3d operator *(const TMatrix3df&, const TVector3d&);
	inline friend double NormAbs(const TVector3d&);

	inline friend int operator <(const TVector3d&, const TVector3d&);
	inline friend int operator ==(const TVector3d&, const TVector3d&);
	inline friend bool operator !=(const TVector3d&, const TVector3d&);
	inline friend bool operator >(const TVector3d&, const TVector3d&);

	inline friend int PracticallyEqual(const TVector3d&, const TVector3d&, double);
	inline friend int PracticallyEqualOrAnti(const TVector3d& v1, const TVector3d& v2, double Tol);
	inline friend int PracticallyParallel(const TVector3d&, const TVector3d&, double);

	inline friend double AngleBwUnitVectors(const TVector3d& v1, const TVector3d& v2, const TVector3d* p_vRefNorm);
};

//-------------------------------------------------------------------------

inline TVector3d operator +(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x+P2.x, P1.y+P2.y, P1.z+P2.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator -(const TVector3d& P1, const TVector3d& P2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P1.x-P2.x, P1.y-P2.y, P1.z-P2.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator *(const double D, const TVector3d& P)
{
	// The following can cause problems with Code Warrior
	return TVector3d(D*P.x, D*P.y, D*P.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator *(const TVector3d& P, const double D)
{
	// The following can cause problems with Code Warrior
	return TVector3d(D*P.x, D*P.y, D*P.z);
}

//-------------------------------------------------------------------------

inline TVector3d operator /(const TVector3d& P, const double D)
{
	// The following can cause problems with Code Warrior
	return TVector3d(P.x/D, P.y/D, P.z/D);
}

//-------------------------------------------------------------------------

inline double operator *(const TVector3d& P1, const TVector3d& P2)
{
	return P1.x*P2.x+P1.y*P2.y+P1.z*P2.z;
}

//-------------------------------------------------------------------------

inline TVector3d operator ^(const TVector3d& v1, const TVector3d& v2)
{
	// The following can cause problems with Code Warrior
	return TVector3d(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

//-------------------------------------------------------------------------

inline int operator <(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x*P1.x + P1.y*P1.y + P1.z*P1.z) < (P2.x*P2.x + P2.y*P2.y + P2.z*P2.z));
}

//-------------------------------------------------------------------------

inline int operator ==(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x == P2.x) && (P1.y == P2.y) && (P1.z == P2.z));
}

//-------------------------------------------------------------------------

inline bool operator !=(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x != P2.x) || (P1.y != P2.y) || (P1.z != P2.z));
}

//-------------------------------------------------------------------------

inline bool operator >(const TVector3d& P1, const TVector3d& P2)
{
	return ((P1.x*P1.x + P1.y*P1.y + P1.z*P1.z) > (P2.x*P2.x + P2.y*P2.y + P2.z*P2.z));
}

//-------------------------------------------------------------------------

inline double NormAbs(const TVector3d& v)
{
	double vx = v.x, vy = v.y, vz = v.z;
	double Abs_vx = (vx>=0.)? vx : -vx;
	double Abs_vy = (vy>=0.)? vy : -vy;
	double Abs_vz = (vz>=0.)? vz : -vz;
	return (Abs_vx > Abs_vy)? ((Abs_vx > Abs_vz)? Abs_vx : Abs_vz) : ((Abs_vy > Abs_vz)? Abs_vy : Abs_vz);
}

//-------------------------------------------------------------------------

inline int PracticallyEqual(const TVector3d& v1, const TVector3d& v2, double Tol)
{
	return (fabs(v2.x-v1.x) < Tol) && (fabs(v2.y-v1.y) < Tol) && (fabs(v2.z-v1.z) < Tol);
}

//-------------------------------------------------------------------------

inline int PracticallyEqualOrAnti(const TVector3d& v1, const TVector3d& v2, double Tol)
{
	return ((fabs(v2.x-v1.x) < Tol) && (fabs(v2.y-v1.y) < Tol) && (fabs(v2.z-v1.z) < Tol))
		|| ((fabs(v2.x+v1.x) < Tol) && (fabs(v2.y+v1.y) < Tol) && (fabs(v2.z+v1.z) < Tol));
}

//-------------------------------------------------------------------------

inline int PracticallyParallel(const TVector3d& v1, const TVector3d& v2, double Tol)
{
	TVector3d v1U = (1./sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z))*v1;
	TVector3d v2U = (1./sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z))*v2;
	TVector3d VectProd = v1U^v2U;
	return (NormAbs(VectProd) < Tol);
}

//-------------------------------------------------------------------------

inline double AngleBwUnitVectors(const TVector3d& v1, const TVector3d& v2, const TVector3d* p_vRefNorm)
{//returns 0 <= angle < 2*Pi
	//TVector3d vE1=v1, vE2=v2;
	//vE1.Normalize(); vE2.Normalize();
	const double TwoPi = 2.*3.141592653589793238;
	double cosAng = v1*v2;

	if(cosAng < -1.) cosAng = -1. + 1.E-13; //OC160209
	else if(cosAng > 1.) cosAng = 1. - 1.E-13;

	double angle = acos(cosAng); // 0 <= angle < Pi
	
	if(p_vRefNorm != 0)
	{
		double testScalProd = (v1^v2)*(*p_vRefNorm);
		if(testScalProd < 0.) angle = TwoPi - angle;
	}
	return angle;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TMatrix3d {
	TVector3d Str0, Str1, Str2;

	TMatrix3d(const TVector3d& InpStr0, const TVector3d& InpStr1, const TVector3d& InpStr2)
	{
		Str0=InpStr0; Str1=InpStr1; Str2=InpStr2;
	}
	TMatrix3d(const double* arElem)
	{//arElem is assumed to contain matrix elements string by string!
		if(arElem == 0) return;
		Str0.x = arElem[0]; Str0.y = arElem[1]; Str0.z = arElem[2]; 
		Str1.x = arElem[3]; Str1.y = arElem[4]; Str1.z = arElem[5]; 
		Str2.x = arElem[6]; Str2.y = arElem[7]; Str2.z = arElem[8]; 
	}
	//TMatrix3d(const TMatrix3d& InM) //??
	//{
	//	*this = InM;
	//}
	TMatrix3d() 
	{
		Str0.x = Str0.y = Str0.z = 0.;
		Str1.x = Str1.y = Str1.z = 0.;
		Str2.x = Str2.y = Str2.z = 0.;
	}

	TMatrix3d(const TMatrix3df&);

	TMatrix3d& operator +=(const TMatrix3d& AnotherMatrix)
	{
		Str0+=AnotherMatrix.Str0; Str1+=AnotherMatrix.Str1; Str2+=AnotherMatrix.Str2; return *this;
	}

	TMatrix3d& operator -=(const TMatrix3d& AnotherMatrix)
	{
		Str0-=AnotherMatrix.Str0; Str1-=AnotherMatrix.Str1; Str2-=AnotherMatrix.Str2; return *this;
	}
	TMatrix3d& operator *=(double a)
	{
		Str0 *= a; Str1 *= a; Str2 *= a; return *this;
	}
	TMatrix3d& operator /=(double a)
	{
		Str0 /= a; Str1 /= a; Str2 /= a; return *this;
	}

	double det()
	{
		return Str0.x*Str1.y*Str2.z + Str0.y*Str1.z*Str2.x + Str0.z*Str1.x*Str2.y
			-Str0.z*Str1.y*Str2.x - Str0.x*Str1.z*Str2.y - Str0.y*Str1.x*Str2.z;
	}
	double absMaxElem()
	{
		double rStr0 = Str0.absMaxElem(), rStr1 = Str1.absMaxElem(), rStr2 = Str2.absMaxElem();
		double res = rStr0;
		if(res < rStr1) res = rStr1;
		if(res < rStr2) res = rStr2;
		return res;
	}
	bool isZero()
	{
		return Str0.isZero() && Str1.isZero() && Str2.isZero();
	}

	friend double detMatrix3d(const TMatrix3d&);
	friend TMatrix3d Matrix3d_inv(const TMatrix3d&);
	friend void Matrix3d_inv(const TMatrix3d&, TMatrix3d&);

	friend TMatrix3d operator +(const TMatrix3d&, const TMatrix3d&);
	friend TMatrix3d operator -(const TMatrix3d&, const TMatrix3d&);
	friend TMatrix3d operator *(double d, const TMatrix3d&);
	friend TVector3d operator *(const TMatrix3d&, const TVector3d&);
	friend TMatrix3d operator *(const TMatrix3d&, const TMatrix3d&);
	friend TMatrix3d operator *(const TMatrix3d&, const TMatrix3df&);
};

//-------------------------------------------------------------------------

inline double detMatrix3d(const TMatrix3d& M)
{
	return M.Str0.x*M.Str1.y*M.Str2.z + M.Str0.y*M.Str1.z*M.Str2.x + M.Str0.z*M.Str1.x*M.Str2.y
		  -M.Str0.z*M.Str1.y*M.Str2.x - M.Str0.x*M.Str1.z*M.Str2.y - M.Str0.y*M.Str1.x*M.Str2.z;
}

//-------------------------------------------------------------------------

inline TMatrix3d Matrix3d_inv(const TMatrix3d& M)
{
	TVector3d St0( M.Str1.y*M.Str2.z-M.Str1.z*M.Str2.y,-M.Str0.y*M.Str2.z+M.Str0.z*M.Str2.y, M.Str0.y*M.Str1.z-M.Str0.z*M.Str1.y);
	TVector3d St1(-M.Str1.x*M.Str2.z+M.Str1.z*M.Str2.x, M.Str0.x*M.Str2.z-M.Str0.z*M.Str2.x,-M.Str0.x*M.Str1.z+M.Str0.z*M.Str1.x);
	TVector3d St2( M.Str1.x*M.Str2.y-M.Str1.y*M.Str2.x,-M.Str0.x*M.Str2.y+M.Str0.y*M.Str2.x, M.Str0.x*M.Str1.y-M.Str0.y*M.Str1.x);
	double invDet = 1./detMatrix3d(M);
	// The following can cause problems with Code Warrior
	return TMatrix3d(invDet*St0, invDet*St1, invDet*St2);
}

//-------------------------------------------------------------------------

inline void Matrix3d_inv(const TMatrix3d& M, TMatrix3d& OutInvMatr)
{
	if(M.Str0.y==0. && (M.Str0.z==0. && (M.Str1.x==0. && (M.Str1.z==0. && (M.Str2.x==0. && M.Str2.y==0.)))))
	{
		OutInvMatr = M;
		double d0 = M.Str0.x;
		double d1 = M.Str1.y;
		double d2 = M.Str2.z;
		if(d0!=0.) OutInvMatr.Str0.x = 1./d0;
		if(d1!=0.) OutInvMatr.Str1.y = 1./d1;
		if(d2!=0.) OutInvMatr.Str2.z = 1./d2;
		return;
	}

	const TVector3d& InSt0 = M.Str0;
	const TVector3d& InSt1 = M.Str1;
	const TVector3d& InSt2 = M.Str2;
	TVector3d St0( InSt1.y*InSt2.z-InSt1.z*InSt2.y,-InSt0.y*InSt2.z+InSt0.z*InSt2.y, InSt0.y*InSt1.z-InSt0.z*InSt1.y);
	TVector3d St1(-InSt1.x*InSt2.z+InSt1.z*InSt2.x, InSt0.x*InSt2.z-InSt0.z*InSt2.x,-InSt0.x*InSt1.z+InSt0.z*InSt1.x);
	TVector3d St2( InSt1.x*InSt2.y-InSt1.y*InSt2.x,-InSt0.x*InSt2.y+InSt0.y*InSt2.x, InSt0.x*InSt1.y-InSt0.y*InSt1.x);
	double invDet = 1./detMatrix3d(M);	// Do something if Det=0 ?

	OutInvMatr.Str0 = invDet*St0;
	OutInvMatr.Str1 = invDet*St1;
	OutInvMatr.Str2 = invDet*St2;
}
//-------------------------------------------------------------------------

inline TMatrix3d operator +(const TMatrix3d& M1, const TMatrix3d& M2)
{
	// The following can cause problems with Code Warrior
	return TMatrix3d(M1.Str0+M2.Str0, M1.Str1+M2.Str1, M1.Str2+M2.Str2);
}

//-------------------------------------------------------------------------

inline TMatrix3d operator -(const TMatrix3d& M1, const TMatrix3d& M2)
{
	// The following can cause problems with Code Warrior
	return TMatrix3d(M1.Str0-M2.Str0, M1.Str1-M2.Str1, M1.Str2-M2.Str2);
}

//-------------------------------------------------------------------------

inline TMatrix3d operator *(double d, const TMatrix3d& M)
{
	// The following can cause problems with Code Warrior
	return TMatrix3d(d*M.Str0, d*M.Str1, d*M.Str2);
}

//-------------------------------------------------------------------------

inline TVector3d operator *(const TMatrix3d& M, const TVector3d& P)
{
	// The following can cause problems with Code Warrior
	return TVector3d(M.Str0*P, M.Str1*P, M.Str2*P);
}

//-------------------------------------------------------------------------

inline TMatrix3d operator *(const TMatrix3d& M1, const TMatrix3d& M2)
{
	const TVector3d& M1Str0 = M1.Str0;
	const TVector3d& M1Str1 = M1.Str1;
	const TVector3d& M1Str2 = M1.Str2;
	const TVector3d& M2Str0 = M2.Str0;
	const TVector3d& M2Str1 = M2.Str1;
	const TVector3d& M2Str2 = M2.Str2;
	TVector3d St0(M1Str0.x*M2Str0.x+M1Str0.y*M2Str1.x+M1Str0.z*M2Str2.x, M1Str0.x*M2Str0.y+M1Str0.y*M2Str1.y+M1Str0.z*M2Str2.y, M1Str0.x*M2Str0.z+M1Str0.y*M2Str1.z+M1Str0.z*M2Str2.z);
	TVector3d St1(M1Str1.x*M2Str0.x+M1Str1.y*M2Str1.x+M1Str1.z*M2Str2.x, M1Str1.x*M2Str0.y+M1Str1.y*M2Str1.y+M1Str1.z*M2Str2.y, M1Str1.x*M2Str0.z+M1Str1.y*M2Str1.z+M1Str1.z*M2Str2.z);
	TVector3d St2(M1Str2.x*M2Str0.x+M1Str2.y*M2Str1.x+M1Str2.z*M2Str2.x, M1Str2.x*M2Str0.y+M1Str2.y*M2Str1.y+M1Str2.z*M2Str2.y, M1Str2.x*M2Str0.z+M1Str2.y*M2Str1.z+M1Str2.z*M2Str2.z);
	// The following can cause problems with Code Warrior
	return TMatrix3d(St0, St1, St2);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TVector2d {
	double x, y;
	TVector2d(double In_x =0., double In_y =0.)
	{
		x = In_x; y = In_y;
	}

	TVector2d& operator +=(const TVector2d& AnotherVect)
	{
		x+=AnotherVect.x; y+=AnotherVect.y; return *this;
	}
	TVector2d& operator -=(const TVector2d& AnotherVect)
	{
		x-=AnotherVect.x; y-=AnotherVect.y; return *this;
	}
	TVector2d& operator *=(double a)
	{
		x *= a; y *= a; return *this;
	}

	inline void Normalize()
	{
		if((x == 0) && (y == 0)) return;
		double invNorm = 1./sqrt(x*x + y*y);
		x *= invNorm; y *= invNorm; 
	}
	inline double Abs() { return sqrt(x*x + y*y);}

	inline friend TVector2d operator +(const TVector2d&, const TVector2d&);
	inline friend TVector2d operator -(const TVector2d&, const TVector2d&);

	inline friend TVector2d operator *(const double, const TVector2d&);
	inline friend double operator *(const TVector2d&, const TVector2d&); // Scalar product

	inline friend int operator <(const TVector2d&, const TVector2d&);
	inline friend int operator ==(const TVector2d&, const TVector2d&);
};

//-------------------------------------------------------------------------

inline double operator *(const TVector2d& P1, const TVector2d& P2)
{
	return P1.x*P2.x+P1.y*P2.y;
}

//-------------------------------------------------------------------------

inline TVector2d operator +(const TVector2d& P1, const TVector2d& P2)
{
	TVector2d OutV(P1.x+P2.x, P1.y+P2.y);
	return OutV;
}

//-------------------------------------------------------------------------

inline TVector2d operator -(const TVector2d& P1, const TVector2d& P2)
{
	TVector2d OutV(P1.x-P2.x, P1.y-P2.y);
	return OutV;
}

//-------------------------------------------------------------------------

inline TVector2d operator *(const double D, const TVector2d& P)
{
	TVector2d OutV(D*P.x, D*P.y);
	return OutV;
}

//-------------------------------------------------------------------------

inline int operator <(const TVector2d& P1, const TVector2d& P2)
{
	return (P1.x*P1.x + P1.y*P1.y < P2.x*P2.x + P2.y*P2.y);
}

//-------------------------------------------------------------------------

inline int operator ==(const TVector2d& P1, const TVector2d& P2)
{
	return ((P1.x == P2.x) && (P1.y == P2.y));
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

struct TMatrix2d {
	TVector2d Str0, Str1;
	TMatrix2d(const TVector2d& InpStr0, const TVector2d& InpStr1)
	{
		Str0=InpStr0; Str1=InpStr1;
	}
	TMatrix2d(double a00, double a01, double a10, double a11)
	{
		Str0.x = a00; Str0.y = a01;
		Str1.x = a10; Str1.y = a11;
	}
	TMatrix2d() {}

	TMatrix2d& operator +=(const TMatrix2d& AnotherMatrix)
	{
		Str0+=AnotherMatrix.Str0; Str1+=AnotherMatrix.Str1; return *this;
	}

	TMatrix2d& operator *=(const TMatrix2d& M2) // right-multiply
	{
        const TVector2d& M2Str0 = M2.Str0;
        const TVector2d& M2Str1 = M2.Str1;
        TVector2d St0(Str0.x*M2Str0.x+Str0.y*M2Str1.x, Str0.x*M2Str0.y+Str0.y*M2Str1.y);
        TVector2d St1(Str1.x*M2Str0.x+Str1.y*M2Str1.x, Str1.x*M2Str0.y+Str1.y*M2Str1.y);
		Str0 = St0; Str1 = St1;
		return *this;
	}

//dddddddddddddddddd
	//friend double detMatrix3d(const TMatrix3d&);
	//friend TMatrix3d Matrix3d_inv(const TMatrix3d&);
	//friend void Matrix3d_inv(const TMatrix3d&, TMatrix3d&);

	//friend TMatrix3d operator +(const TMatrix3d&, const TMatrix3d&);
	//friend TMatrix3d operator -(const TMatrix3d&, const TMatrix3d&);
	//friend TMatrix3d operator *(double d, const TMatrix3d&);
	//friend TVector3d operator *(const TMatrix3d&, const TVector3d&);

    friend TMatrix2d operator *(const TMatrix2d&, const TMatrix2d&);
};

//-------------------------------------------------------------------------

//inline double detMatrix3d(const TMatrix3d& M)
//{
//	return M.Str0.x*M.Str1.y*M.Str2.z + M.Str0.y*M.Str1.z*M.Str2.x + M.Str0.z*M.Str1.x*M.Str2.y
//		  -M.Str0.z*M.Str1.y*M.Str2.x - M.Str0.x*M.Str1.z*M.Str2.y - M.Str0.y*M.Str1.x*M.Str2.z;
//}

//-------------------------------------------------------------------------

//inline TMatrix3d Matrix3d_inv(const TMatrix3d& M)
//{
//	double invDet = 1./detMatrix3d(M);
//	TVector3d St0(invDet*(M.Str1.y*M.Str2.z-M.Str1.z*M.Str2.y), invDet*(-M.Str0.y*M.Str2.z+M.Str0.z*M.Str2.y), invDet*(M.Str0.y*M.Str1.z-M.Str0.z*M.Str1.y));
//	TVector3d St1(invDet*(-M.Str1.x*M.Str2.z+M.Str1.z*M.Str2.x), invDet*(M.Str0.x*M.Str2.z-M.Str0.z*M.Str2.x), invDet*(-M.Str0.x*M.Str1.z+M.Str0.z*M.Str1.x));
//	TVector3d St2(invDet*(M.Str1.x*M.Str2.y-M.Str1.y*M.Str2.x), invDet*(-M.Str0.x*M.Str2.y+M.Str0.y*M.Str2.x), invDet*(M.Str0.x*M.Str1.y-M.Str0.y*M.Str1.x));
//	// The following can cause problems with Code Warrior
//	return TMatrix3d(St0, St1, St2);
//}

//-------------------------------------------------------------------------

//inline void Matrix3d_inv(const TMatrix3d& M, TMatrix3d& OutInvMatr)
//{
//	if(M.Str0.y==0. && (M.Str0.z==0. && (M.Str1.x==0. && (M.Str1.z==0. && (M.Str2.x==0. && M.Str2.y==0.)))))
//	{
//		OutInvMatr = M;
//		double d0 = M.Str0.x;
//		double d1 = M.Str1.y;
//		double d2 = M.Str2.z;
//		if(d0!=0.) OutInvMatr.Str0.x = 1./d0;
//		if(d1!=0.) OutInvMatr.Str1.y = 1./d1;
//		if(d2!=0.) OutInvMatr.Str2.z = 1./d2;
//		return;
//	}
//
//	const TVector3d& InSt0 = M.Str0;
//	const TVector3d& InSt1 = M.Str1;
//	const TVector3d& InSt2 = M.Str2;
//	TVector3d St0( InSt1.y*InSt2.z-InSt1.z*InSt2.y,-InSt0.y*InSt2.z+InSt0.z*InSt2.y, InSt0.y*InSt1.z-InSt0.z*InSt1.y);
//	TVector3d St1(-InSt1.x*InSt2.z+InSt1.z*InSt2.x, InSt0.x*InSt2.z-InSt0.z*InSt2.x,-InSt0.x*InSt1.z+InSt0.z*InSt1.x);
//	TVector3d St2( InSt1.x*InSt2.y-InSt1.y*InSt2.x,-InSt0.x*InSt2.y+InSt0.y*InSt2.x, InSt0.x*InSt1.y-InSt0.y*InSt1.x);
//	double invDet = 1./detMatrix3d(M);	// Do something if Det=0 ?
//
//	OutInvMatr.Str0 = Mult(invDet, St0);
//	OutInvMatr.Str1 = Mult(invDet, St1);
//	OutInvMatr.Str2 = Mult(invDet, St2);
//}

//-------------------------------------------------------------------------

//inline TMatrix3d operator +(const TMatrix3d& M1, const TMatrix3d& M2)
//{
//	// The following can cause problems with Code Warrior
//	return TMatrix3d(Plus(M1.Str0, M2.Str0), Plus(M1.Str1, M2.Str1), Plus(M1.Str2, M2.Str2));
//}

//-------------------------------------------------------------------------

//inline TMatrix3d operator -(const TMatrix3d& M1, const TMatrix3d& M2)
//{
//	// The following can cause problems with Code Warrior
//	return TMatrix3d(Minus(M1.Str0, M2.Str0), Minus(M1.Str1, M2.Str1), Minus(M1.Str2, M2.Str2));
//}

//-------------------------------------------------------------------------

//inline TMatrix3d operator *(double d, const TMatrix3d& M)
//{
//	// The following can cause problems with Code Warrior
//	return TMatrix3d(Mult(d, M.Str0), Mult(d, M.Str1), Mult(d, M.Str2));
//}

//-------------------------------------------------------------------------

//inline TVector3d operator *(const TMatrix3d& M, const TVector3d& P)
//{
//	// The following can cause problems with Code Warrior
//	return TVector3d(Mult(M.Str0, P), Mult(M.Str1, P), Mult(M.Str2, P));
//}

//-------------------------------------------------------------------------

inline TMatrix2d operator *(const TMatrix2d& M1, const TMatrix2d& M2)
{
	const TVector2d& M1Str0 = M1.Str0;
	const TVector2d& M1Str1 = M1.Str1;
	const TVector2d& M2Str0 = M2.Str0;
	const TVector2d& M2Str1 = M2.Str1;
	TVector2d St0(M1Str0.x*M2Str0.x+M1Str0.y*M2Str1.x, M1Str0.x*M2Str0.y+M1Str0.y*M2Str1.y);
	TVector2d St1(M1Str1.x*M2Str0.x+M1Str1.y*M2Str1.x, M1Str1.x*M2Str0.y+M1Str1.y*M2Str1.y);
	// The following can cause problems with Code Warrior
	return TMatrix2d(St0, St1);
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

enum TAxisOrient { ParallelToX, ParallelToY, ParallelToZ };

//-------------------------------------------------------------------------

#endif
