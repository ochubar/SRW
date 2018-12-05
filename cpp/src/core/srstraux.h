/************************************************************************//**
 * File: srstraux.h
 * Description: Auxiliary structures for various SR calculation methods (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRSTRAUX_H
#define __SRSTRAUX_H

#include "srradstr.h"
#include "srwfrsmp.h"
#include "gmvect.h"
#include "srebmdat.h"
#include "srradind.h"
#include "srinterf.h"
#include "cmplxd.h"
#include "srercode.h"
#include "srwlib.h"

#include <vector>

using namespace std;

//*************************************************************************

struct srLambXYZ {
	double Lamb, x, y, z;
};

//*************************************************************************

struct srTEXZ {
	double e, x, z;
	char VsXorZ;
	//long aux_offset;
	long long aux_offset;
};

//*************************************************************************

struct srTEXZY {
	double e, x, z, y;
	double dx, dz;
	srTEXZY() { e = x = z = y = dx = dz = 0;}
};

//*************************************************************************

struct srTStNoFiNo {
	//int StNo, FiNo, StOffset;
	long long StNo, FiNo, StOffset;
	//srTStNoFiNo(int InStNo, int InFiNo) { StNo = InStNo; FiNo = InFiNo; StOffset = 0;}
	srTStNoFiNo(long long InStNo, long long InFiNo) { StNo = InStNo; FiNo = InFiNo; StOffset = 0;}
	//srTStNoFiNo(int InStNo, int InFiNo, int InStOffset) 
	srTStNoFiNo(long long InStNo, long long InFiNo, long long InStOffset) 
	{ 
		StNo = InStNo; FiNo = InFiNo; StOffset = InStNo; StOffset = InStOffset;
	}
	srTStNoFiNo() { StNo = FiNo = StOffset = 0;}

	inline friend int operator <(const srTStNoFiNo&, const srTStNoFiNo&);
	inline friend int operator ==(const srTStNoFiNo&, const srTStNoFiNo&);
};

inline int operator <(const srTStNoFiNo& P1, const srTStNoFiNo& P2)
{
	return P1.StNo < P2.StNo;
}

inline int operator ==(const srTStNoFiNo& P1, const srTStNoFiNo& P2)
{
	return (P1.StNo == P2.StNo) && (P1.FiNo == P2.FiNo) && (P1.StOffset == P2.StOffset);
}

//*************************************************************************

struct srTEFourier {
	double EwX_Re, EwX_Im, EwZ_Re, EwZ_Im;
	srTEFourier(double InEwX_Re =0., double InEwX_Im =0., double InEwZ_Re =0., double InEwZ_Im =0.) 
	{ 
		EwX_Re = InEwX_Re; EwX_Im = InEwX_Im; EwZ_Re = InEwZ_Re; EwZ_Im = InEwZ_Im;
	}

	void ZeroAll()
	{
		EwX_Re = EwX_Im = EwZ_Re = EwZ_Im = 0.;
	}

	srTEFourier& operator +=(const srTEFourier& AnotherE)
	{
		EwX_Re += AnotherE.EwX_Re; EwX_Im += AnotherE.EwX_Im; EwZ_Re += AnotherE.EwZ_Re; EwZ_Im += AnotherE.EwZ_Im; 
		return *this;
	}
	srTEFourier& operator *=(const double D)
	{
		EwX_Re *= D; EwX_Im *= D; EwZ_Re *= D; EwZ_Im *= D; 
		return *this;
	}

	inline friend srTEFourier operator +(const srTEFourier&, const srTEFourier&);
	inline friend srTEFourier operator -(const srTEFourier&, const srTEFourier&);
	inline friend srTEFourier operator *(const double, const srTEFourier&);
	inline friend int operator <(const srTEFourier&, const srTEFourier&);
	inline friend int operator ==(const srTEFourier&, const srTEFourier&);
};

inline srTEFourier operator +(const srTEFourier& E1, const srTEFourier& E2)
{
	return srTEFourier(E1.EwX_Re + E2.EwX_Re, E1.EwX_Im + E2.EwX_Im, E1.EwZ_Re + E2.EwZ_Re, E1.EwZ_Im + E2.EwZ_Im);
}
inline srTEFourier operator -(const srTEFourier& E1, const srTEFourier& E2)
{
	return srTEFourier(E1.EwX_Re - E2.EwX_Re, E1.EwX_Im - E2.EwX_Im, E1.EwZ_Re - E2.EwZ_Re, E1.EwZ_Im - E2.EwZ_Im);
}
inline srTEFourier operator *(const double d, const srTEFourier& E1)
{
	return srTEFourier(d*E1.EwX_Re, d*E1.EwX_Im, d*E1.EwZ_Re, d*E1.EwZ_Im);
}

inline int operator <(const srTEFourier& P1, const srTEFourier& P2)
{
	return P1.EwX_Re < P2.EwX_Re;
}
inline int operator ==(const srTEFourier& P1, const srTEFourier& P2)
{
	return (P1.EwX_Re == P2.EwX_Re) && (P1.EwX_Im == P2.EwX_Im) && (P1.EwZ_Re == P2.EwZ_Re) && (P1.EwZ_Im == P2.EwZ_Im);
}

//*************************************************************************

struct srTStokes {
	double s0, s1, s2, s3;
	srTStokes(double In_s0 =0., double In_s1 =0., double In_s2 =0., double In_s3 =0.) 
	{ 
		s0 = In_s0; s1 = In_s1; s2 = In_s2; s3 = In_s3;
	}

	void ZeroAll()
	{
		s0 = s1 = s2 = s3 = 0.;
	}

	srTStokes& operator +=(const srTStokes& AnotherStokes)
	{
		s0 += AnotherStokes.s0; s1 += AnotherStokes.s1; s2 += AnotherStokes.s2; s3 += AnotherStokes.s3; 
		return *this;
	}
	srTStokes& operator *=(const double D)
	{
		s0 *= D; s1 *= D; s2 *= D; s3 *= D; 
		return *this;
	}

	inline friend srTStokes operator +(const srTStokes&, const srTStokes&);
	inline friend srTStokes operator -(const srTStokes&, const srTStokes&);
	inline friend srTStokes operator *(double, const srTStokes&);
	inline friend srTStokes operator *(const srTStokes&, double);

	inline friend int operator <(const srTStokes&, const srTStokes&);
	inline friend int operator ==(const srTStokes&, const srTStokes&);
};

inline srTStokes operator +(const srTStokes& E1, const srTStokes& E2)
{
	return srTStokes(E1.s0 + E2.s0, E1.s1 + E2.s1, E1.s2 + E2.s2, E1.s3 + E2.s3);
}
inline srTStokes operator -(const srTStokes& E1, const srTStokes& E2)
{
	return srTStokes(E1.s0 - E2.s0, E1.s1 - E2.s1, E1.s2 - E2.s2, E1.s3 - E2.s3);
}
inline srTStokes operator *(double a, const srTStokes& E)
{
	return srTStokes(a*E.s0, a*E.s1, a*E.s2, a*E.s3);
}
inline srTStokes operator *(const srTStokes& E, double a)
{
	return srTStokes(a*E.s0, a*E.s1, a*E.s2, a*E.s3);
}

inline int operator <(const srTStokes& P1, const srTStokes& P2)
{
	return P1.s0 < P2.s0;
}
inline int operator ==(const srTStokes& P1, const srTStokes& P2)
{
	return (P1.s0 == P2.s0) && (P1.s1 == P2.s1) && (P1.s2 == P2.s2) && (P1.s3 == P2.s3);
}

//*************************************************************************

struct srTStokesC {
	TComplexD s0, s1, s2, s3;
};

//*************************************************************************

struct srTEFieldPtrs {
	float *pExRe, *pExIm, *pEzRe, *pEzIm;
	srTEFieldPtrs(float* In_pExRe =0, float* In_pExIm =0, float* In_pEzRe =0, float* In_pEzIm =0) 
	{ 
		pExRe = In_pExRe; pExIm = In_pExIm; pEzRe = In_pEzRe; pEzIm = In_pEzIm;
	}
};

//*************************************************************************

struct srTEFieldPtrsX {
	float *pExReHorL, *pExImHorL, *pEyReHorL, *pEyImHorL;
	float *pExReHorR, *pExImHorR, *pEyReHorR, *pEyImHorR;
	float *pExReVerD, *pExImVerD, *pEyReVerD, *pEyImVerD;
	float *pExReVerU, *pExImVerU, *pEyReVerU, *pEyImVerU;
	float *pExReC, *pExImC, *pEyReC, *pEyImC;
	double x, y, ePh, RobsX, RobsY;
	double xStep, yStep;
	char WaveFrontTermCanBeTreated;

	srTEFieldPtrsX() 
	{ 
		pExReHorL = pExImHorL = pEyReHorL = pEyImHorL = 0;
		pExReHorR = pExImHorR = pEyReHorR = pEyImHorR = 0;
		pExReVerD = pExImVerD = pEyReVerD = pEyImVerD = 0;
		pExReVerU = pExImVerU = pEyReVerU = pEyImVerU = 0;
		pExReC = pExImC = pEyReC = pEyImC = 0;
		WaveFrontTermCanBeTreated = 0;
	}
	void Initialize()
	{
		pExReHorL = pExImHorL = pEyReHorL = pEyImHorL = 0;
		pExReHorR = pExImHorR = pEyReHorR = pEyImHorR = 0;
		pExReVerD = pExImVerD = pEyReVerD = pEyImVerD = 0;
		pExReVerU = pExImVerU = pEyReVerU = pEyImVerU = 0;
		pExReC = pExImC = pEyReC = pEyImC = 0;
		WaveFrontTermCanBeTreated = 0;
	}
};

//*************************************************************************

struct srTSegmOfMesh2D {
	double VarSt, VarFi, Con;
	char NeedsFurtherSubd;
	srTSegmOfMesh2D() { NeedsFurtherSubd = 1;}
};

//*************************************************************************

struct srTCellOfMesh2D {
	srTSegmOfMesh2D Side1, Side2, Side3, Side4;
	srTEFourier f1, f2, f3, f4;
};

//*************************************************************************

typedef vector<srTStNoFiNo, allocator<srTStNoFiNo> > srTStNoFiNoVect;
typedef vector<srTCellOfMesh2D, allocator<srTCellOfMesh2D> > srTMesh2DVect;
typedef vector<srTEFourier, allocator<srTEFourier> > srTEFourierVect;
typedef vector<char*, allocator<char*> > srTStringVect;
typedef vector<int, allocator<int> > srTIntVect;
typedef vector<double, allocator<double> > srTDoubleVect;

//*************************************************************************

#ifdef __MWERKS__
/*
null_template
struct iterator_traits <srTEFourier*> {
     typedef ptrdiff_t difference_type;
     typedef srTEFourier value_type;
     typedef srTEFourier* pointer;
     typedef srTEFourier& reference;
     typedef random_access_iterator_tag iterator_category;
};

null_template
struct iterator_traits <const srTEFourier*> {
     typedef ptrdiff_t                    difference_type;
     typedef const srTEFourier            value_type;
     typedef const srTEFourier*           pointer;
     typedef const srTEFourier&           reference;
     typedef random_access_iterator_tag   iterator_category;
};

null_template
struct iterator_traits <srTStNoFiNo*> {
     typedef ptrdiff_t difference_type;
     typedef srTStNoFiNo value_type;
     typedef srTStNoFiNo* pointer;
     typedef srTStNoFiNo& reference;
     typedef random_access_iterator_tag   iterator_category;
};

null_template
struct iterator_traits <char**> {
     typedef ptrdiff_t difference_type;
     typedef char* value_type;
     typedef char** pointer;
     typedef char*& reference;
     typedef random_access_iterator_tag iterator_category;
};
*/
#endif

//*************************************************************************

struct srTRadIntervVal {
	//int InitLevelNo, SecondInt;
	long long InitLevelNo, SecondInt;
	double sStart, sEnd;
};

//*************************************************************************

struct srTSMeshData {
	double s;
	//int LevelNo, NoOfPoOnLev;
	long long LevelNo, NoOfPoOnLev;
};

//*************************************************************************

struct srTMomentsPtrs {
	//float *pTotPhot, *pX, *pXP, *pZ, *pZP, *pXX, *pXXP, *pXPXP, *pZZ, *pZZP, *pZPZP;
	//float Mxx, Mxxp, Mxpxp, Mzz, Mzzp, Mzpzp;
	//float SqrtMxx, SqrtMxpxp, SqrtMzz, SqrtMzpzp;
	//OC130311
	double *pTotPhot, *pX, *pXP, *pZ, *pZP, *pXX, *pXXP, *pXPXP, *pZZ, *pZZP, *pZPZP;
	double Mxx, Mxxp, Mxpxp, Mzz, Mzzp, Mzpzp;
	double SqrtMxx, SqrtMxpxp, SqrtMzz, SqrtMzpzp;
	bool precCenMomIsOK;

	//srTMomentsPtrs(float* tMom, int ie=0)
	srTMomentsPtrs(double* tMom, int ie=0) //OC130311
	{
		const int numStatMom = 11;
		if(tMom != 0)
		{
			tMom += numStatMom*ie;

			pTotPhot = tMom++;

			pX = tMom++; pXP = tMom++; 
			pZ = tMom++; pZP = tMom++; 
			pXX = tMom++; pXXP = tMom++; pXPXP = tMom++; 
			pZZ = tMom++; pZZP = tMom++; pZPZP = tMom++;

			ComputeCentralMoments();
		}
	}
	srTMomentsPtrs() {}

	void ComputeCentralMoments()
	{
		Mxx = *pXX - (*pX)*(*pX);
		Mxxp = *pXXP - (*pX)*(*pXP);
		Mxpxp = *pXPXP - (*pXP)*(*pXP);
		Mzz = *pZZ - (*pZ)*(*pZ);
		Mzzp = *pZZP - (*pZ)*(*pZP);
		Mzpzp = *pZPZP - (*pZP)*(*pZP);

		precCenMomIsOK = true;
		if((Mxx < 0) || (Mxpxp < 0) || (Mzz < 0) || (Mzpzp < 0)) precCenMomIsOK = false;

		//OC test (trying to walk-around problem with economic wavefront propag.) 280809
		double SafeOffsetCoef = 0.1;
		//OC130311 float->double
		if(Mxx < SafeOffsetCoef*(*pXX)) Mxx = SafeOffsetCoef*(*pXX);
		if(Mxxp < SafeOffsetCoef*(*pXXP)) Mxxp = SafeOffsetCoef*(*pXXP);
		if(Mxpxp < SafeOffsetCoef*(*pXPXP)) Mxpxp = SafeOffsetCoef*(*pXPXP);
		if(Mzz < SafeOffsetCoef*(*pZZ)) Mzz = SafeOffsetCoef*(*pZZ);
		if(Mzzp < SafeOffsetCoef*(*pZZP)) Mzzp = SafeOffsetCoef*(*pZZP);
		if(Mzpzp < SafeOffsetCoef*(*pZPZP)) Mzpzp = SafeOffsetCoef*(*pZPZP);

		SqrtMxx = sqrt(::fabs(Mxx)); SqrtMxpxp = sqrt(::fabs(Mxpxp));
		SqrtMzz = sqrt(::fabs(Mzz)); SqrtMzpzp = sqrt(::fabs(Mzpzp));
	}
};

//*************************************************************************

struct srTMomentsVals {
	//float TotPhot, X, XP, Z, ZP, XX, XXP, XPXP, ZZ, ZZP, ZPZP;
	//float Mxx, Mxxp, Mxpxp, Mzz, Mzzp, Mzpzp;
	//float SqrtMxx, SqrtMxpxp, SqrtMzz, SqrtMzpzp;
	//OC130311
	double TotPhot, X, XP, Z, ZP, XX, XXP, XPXP, ZZ, ZZP, ZPZP;
	double Mxx, Mxxp, Mxpxp, Mzz, Mzzp, Mzpzp;
	double SqrtMxx, SqrtMxpxp, SqrtMzz, SqrtMzpzp;

	srTMomentsVals(srTMomentsPtrs& MomentsPtrs)
	{
		TotPhot = *(MomentsPtrs.pTotPhot);
		X = *(MomentsPtrs.pX); XP = *(MomentsPtrs.pXP); 
		Z = *(MomentsPtrs.pZ); ZP = *(MomentsPtrs.pZP); 
		XX = *(MomentsPtrs.pXX); XXP = *(MomentsPtrs.pXXP); XPXP = *(MomentsPtrs.pXPXP); 
		ZZ = *(MomentsPtrs.pZZ); ZZP = *(MomentsPtrs.pZZP); ZPZP = *(MomentsPtrs.pZPZP);

		Mxx = MomentsPtrs.Mxx; Mxxp = MomentsPtrs.Mxxp; Mxpxp = MomentsPtrs.Mxpxp;
		Mzz = MomentsPtrs.Mzz; Mzzp = MomentsPtrs.Mzzp; Mzpzp = MomentsPtrs.Mzpzp;
		SqrtMxx = MomentsPtrs.SqrtMxx; SqrtMxpxp = MomentsPtrs.SqrtMxpxp;
		SqrtMzz = MomentsPtrs.SqrtMzz; SqrtMzpzp = MomentsPtrs.SqrtMzpzp;
	}
	srTMomentsVals() {}
};

//*************************************************************************

struct srTMomentsRatios {
	//float RxxMomX, RxpxpMomX, RzzMomX, RzpzpMomX;
	//float RxxMomZ, RxpxpMomZ, RzzMomZ, RzpzpMomZ;
	double RxxMomX, RxpxpMomX, RzzMomX, RzpzpMomX; //OC130311
	double RxxMomZ, RxpxpMomZ, RzzMomZ, RzpzpMomZ;

	void GiveMaxRatios(double& RxxMax, double& RxpxpMax, double& RzzMax, double& RzpzpMax)
	{
		RxxMax = (RxxMomX > RxxMomZ)? RxxMomX : RxxMomZ;
		RxpxpMax = (RxpxpMomX > RxpxpMomZ)? RxpxpMomX : RxpxpMomZ;
		RzzMax = (RzzMomX > RzzMomZ)? RzzMomX : RzzMomZ;
		RzpzpMax = (RzpzpMomX > RzpzpMomZ)? RzpzpMomX : RzpzpMomZ;
	}
};

//*************************************************************************

struct srTMinMaxEParam {

	float MaxReEx, MaxImEx, MaxReEz, MaxImEz, MinReEx, MinImEx, MinReEz, MinImEz;
	long xIndMaxReEx, xIndMaxImEx, xIndMaxReEz, xIndMaxImEz, xIndMinReEx, xIndMinImEx, xIndMinReEz, xIndMinImEz;
	long zIndMaxReEx, zIndMaxImEx, zIndMaxReEz, zIndMaxImEz, zIndMinReEx, zIndMinImEx, zIndMinReEz, zIndMinImEz;

	srTMinMaxEParam()
	{
		xIndMaxReEx = xIndMaxImEx = xIndMaxReEz = xIndMaxImEz = xIndMinReEx = xIndMinImEx = xIndMinReEz = xIndMinImEz = 0;
		zIndMaxReEx = zIndMaxImEx = zIndMaxReEz = zIndMaxImEz = zIndMinReEx = zIndMinImEx = zIndMinReEz = zIndMinImEz = 0;
	}

	//void FindAbsMaxAmongReAndIm(float& MaxAbsEx, long& xIndMaxAbsEx, long& zIndMaxAbsEx, float& MaxAbsEz, long& xIndMaxAbsEz, long& zIndMaxAbsEz)
	void FindAbsMaxAmongReAndIm(float& MaxAbsEx, long long& xIndMaxAbsEx, long long& zIndMaxAbsEx, float& MaxAbsEz, long long& xIndMaxAbsEz, long long& zIndMaxAbsEz)
	{
		float AbsMaxReEx, AbsMaxImEx;
		//long xIndAbsMaxReEx, zIndAbsMaxReEx, xIndAbsMaxImEx, zIndAbsMaxImEx;
		long long xIndAbsMaxReEx, zIndAbsMaxReEx, xIndAbsMaxImEx, zIndAbsMaxImEx;
		if(::fabs(MaxReEx) > ::fabs(MinReEx))
		{
			AbsMaxReEx = (float)::fabs(MaxReEx); xIndAbsMaxReEx = xIndMaxReEx; zIndAbsMaxReEx = zIndMaxReEx; 
		}
		else
		{
			AbsMaxReEx = (float)::fabs(MinReEx); xIndAbsMaxReEx = xIndMinReEx; zIndAbsMaxReEx = zIndMinReEx; 
		}

		if(::fabs(MaxImEx) > ::fabs(MinImEx))
		{
			AbsMaxImEx = (float)::fabs(MaxImEx); xIndAbsMaxImEx = xIndMaxImEx; zIndAbsMaxImEx = zIndMaxImEx; 
		}
		else
		{
			AbsMaxImEx = (float)::fabs(MinImEx); xIndAbsMaxImEx = xIndMinImEx; zIndAbsMaxImEx = zIndMinImEx; 
		}

		if(AbsMaxReEx > AbsMaxImEx)
		{
			MaxAbsEx = AbsMaxReEx; xIndMaxAbsEx = xIndAbsMaxReEx; zIndMaxAbsEx = zIndAbsMaxReEx;
		}
		else
		{
			MaxAbsEx = AbsMaxImEx; xIndMaxAbsEx = xIndAbsMaxImEx; zIndMaxAbsEx = zIndAbsMaxImEx;
		}

		float AbsMaxReEz, AbsMaxImEz;
		//long xIndAbsMaxReEz, zIndAbsMaxReEz, xIndAbsMaxImEz, zIndAbsMaxImEz;
		long long xIndAbsMaxReEz, zIndAbsMaxReEz, xIndAbsMaxImEz, zIndAbsMaxImEz;

		if(::fabs(MaxReEz) > ::fabs(MinReEz))
		{
			AbsMaxReEz = (float)::fabs(MaxReEz); xIndAbsMaxReEz = xIndMaxReEz; zIndAbsMaxReEz = zIndMaxReEz; 
		}
		else
		{
			AbsMaxReEz = (float)::fabs(MinReEz); xIndAbsMaxReEz = xIndMinReEz; zIndAbsMaxReEz = zIndMinReEz; 
		}

		if(::fabs(MaxImEz) > ::fabs(MinImEz))
		{
			AbsMaxImEz = (float)::fabs(MaxImEz); xIndAbsMaxImEz = xIndMaxImEz; zIndAbsMaxImEz = zIndMaxImEz; 
		}
		else
		{
			AbsMaxImEz = (float)::fabs(MinImEz); xIndAbsMaxImEz = xIndMinImEz; zIndAbsMaxImEz = zIndMinImEz; 
		}

		if(AbsMaxReEz > AbsMaxImEz)
		{
			MaxAbsEz = AbsMaxReEz; xIndMaxAbsEz = xIndAbsMaxReEz; zIndMaxAbsEz = zIndAbsMaxReEz;
		}
		else
		{
			MaxAbsEz = AbsMaxImEz; xIndMaxAbsEz = xIndAbsMaxImEz; zIndMaxAbsEz = zIndAbsMaxImEz;
		}
	}

	//void FindGenMaxAbsE(float& MaxAbsE, long& xIndMaxAbsE, long& zIndMaxAbsE)
	void FindGenMaxAbsE(float& MaxAbsE, long long& xIndMaxAbsE, long long& zIndMaxAbsE)
	{
		float MaxAbsEx, MaxAbsEz;
		//long xIndMaxAbsEx, zIndMaxAbsEx, xIndMaxAbsEz, zIndMaxAbsEz;
		long long xIndMaxAbsEx, zIndMaxAbsEx, xIndMaxAbsEz, zIndMaxAbsEz;
		FindAbsMaxAmongReAndIm(MaxAbsEx, xIndMaxAbsEx, zIndMaxAbsEx, MaxAbsEz, xIndMaxAbsEz, zIndMaxAbsEz);
		if(MaxAbsEx > MaxAbsEz)
		{
			MaxAbsE = MaxAbsEx; xIndMaxAbsE = xIndMaxAbsEx; zIndMaxAbsE = zIndMaxAbsEx;
		}
		else
		{
			MaxAbsE = MaxAbsEz; xIndMaxAbsE = xIndMaxAbsEz; zIndMaxAbsE = zIndMaxAbsEz;
		}
	}
};

//*************************************************************************

struct srTRadResize1D {
	double pm, pd;

	double RelCenPos; // 0 <= RelCenPos <= 1
	double RelCenPosTol; // Never change this in code
	char UseOtherSideFFT;
	char DoNotTreatSpherTerm;
	srTRadResize1D() 
	{
		pm = pd = 1.;
		RelCenPos = 0.5;
		RelCenPosTol = 0.1; // To steer only here		
		UseOtherSideFFT = 0;
		DoNotTreatSpherTerm = 0;
	}
	void SetRelCenPosTolToEnforceCenterCorr()
	{
		double CurRelCenPosDif = ::fabs(RelCenPos - 0.5);
		if(CurRelCenPosDif < 1.E-06) return;
		RelCenPosTol = 0.5*CurRelCenPosDif;
	}
};

//*************************************************************************

struct srTRadResize {
	//double pxm, pxd, pzm, pzd;
	double pem, ped, pxm, pxd, pzm, pzd; //OC111111

	//double RelCenPosX, RelCenPosZ; // 0 <= RelCenPos <= 1
	double RelCenPosE, RelCenPosX, RelCenPosZ; //OC111111
	double RelCenPosTol; // Never change this in code

	char ShiftTypeBeforeRes; //OC090311
	double eCenShift, xCenShift, zCenShift;
	
	//char UseOtherSideFFT;
	//char DoNotTreatSpherTerm;
	char ModeBits; //OC090311
	double PropAutoPrec; //OC090311

	//OC011213
	double vLxOut, vLyOut, vLzOut; //Coordinates of the output Optical Axis vector
	double vHxOut, vHyOut; //Coordinates of the Horizontal Base vector of the output frame

	srTRadResize() 
	{
		pem = ped = pxm = pxd = pzm = pzd = 1.;

		RelCenPosE = RelCenPosX = RelCenPosZ = 0.5;
		RelCenPosTol = 1.e-06; //0.001; //OC141014 //0.1; // To steer only here	
		
		PropAutoPrec = 1.;
		ShiftTypeBeforeRes = 0;
		eCenShift = xCenShift = zCenShift = 0.;

		//UseOtherSideFFT = 0;
		//DoNotTreatSpherTerm = 0;
		ModeBits = 0;
		//bits: 
		//#0- UseOtherSideFFT
		//#1- DoNotTreatSpherTerm
		//#2- Propagation: Auto-Resize Before Propagation
		//#3- Propagation: Auto-Resize After Propagation
		//#4- Propagation: Allow Under-Sampling Mode

		//OC011213
		vLxOut = vLyOut = vLzOut = 0; //Default coordinates of the output Optical Axis vector
		vHxOut = vHyOut = 0; //Default coordinates of the Horizontal Base vector of the output frame
	}

	char useOtherSideFFT(int in=-1) 
	{
		if(in == 0)     ModeBits &= 126; //30;//1111110
		else if(in > 0) ModeBits |= 1; //0000001
		return (1 & ModeBits);
	}
	char doNotTreatSpherTerm(int in=-1) 
	{ 
		if(in == 0)     ModeBits &= 125; //29;//1111101
		else if(in > 0) ModeBits |= 2; //0000010
		return (2 & ModeBits) >> 1;
	}
	char propAutoResizeBefore(int in=-1) 
	{ 
		if(in == 0)     ModeBits &= 123; //27;//1111011
		else if(in > 0) ModeBits |= 4; //0000100
		return (4 & ModeBits) >> 2;
	}
	char propAutoResizeAfter(int in=-1) 
	{
		if(in == 0)     ModeBits &= 119; //23;//1110111
		else if(in > 0) ModeBits |= 8; //0001000
		return (8 & ModeBits) >> 3;
	}
	char propAllowUnderSamp(int in=-1) 
	{//uses bytes 4, 5 for encoding type of analytical treatment at propagation with resizing
		if(in == 0)     ModeBits &= 15;//0001111
		//else if(in > 0) ModeBits |= 16;//10000
		//return (16 & ModeBits) >> 4;
		else if(in > 0) 
		{
			in <<= 4;
			ModeBits |= in;//XXX0000
		}
		return (112 & ModeBits) >> 4;
	}
};

typedef vector<srTRadResize, allocator<srTRadResize> > srTRadResizeVect;

//*************************************************************************

struct srTSRWRadStructWaveNames {

	char NameRad[MAX_OBJ_NAME+1];
	char NameRadX[MAX_OBJ_NAME+1], NameRadZ[MAX_OBJ_NAME+1];
	char NameElecBeam[MAX_OBJ_NAME+1], NameTrj[MAX_OBJ_NAME+1];
	char Name4x4PropMatr[MAX_OBJ_NAME+1];
	char NameMomX[MAX_OBJ_NAME+1], NameMomZ[MAX_OBJ_NAME+1];
	char NameWfrAuxData[MAX_OBJ_NAME+1];

	srTSRWRadStructWaveNames()
	{
		*NameRad = '\0';
		*NameRadX = '\0'; *NameRadZ = '\0';
		*NameElecBeam = '\0'; *NameTrj = '\0';
		*Name4x4PropMatr = '\0';
		*NameMomX = '\0'; *NameMomZ = '\0';
		*NameWfrAuxData = '\0';
	}
};

//*************************************************************************

struct srTSRWRadStructWaveKeys {

	char wRad_;
	char wRadX_, wRadZ_;
	char wElecBeam_, wTrj_;
	char w4x4PropMatr_;
	char wMomX_, wMomZ_;
	char wWfrAuxData_;

	srTSRWRadStructWaveKeys()
	{
		wRad_= wRadX_= wRadZ_= wElecBeam_= wTrj_= w4x4PropMatr_= wMomX_= wMomZ_= wWfrAuxData_=0;
	}
};

//*************************************************************************

struct srTRadSect1D {

	float *pEx, *pEz;
	double ArgStep, ArgStart, ArgStartTr;
	long np;

	double eVal, OtherCoordVal;
	char VsXorZ;
	long icOtherCoord;

	double Robs;
	double RobsAbsErr;
	double cArg;

	double WfrMin, WfrMax; // Exact borders of the Wavefront
	char WfrEdgeCorrShouldBeDone; // To switch off/on manually

	char ThereIsUnderSamplingIn2D;

	char Pres; // 0- Coord, 1- Ang.
	char LengthUnit; // 0- m; 1- mm; 
	char PhotEnergyUnit; // 0- eV; 1- keV;

	long AuxLong1, AuxLong2, AuxLong3, AuxLong4;

	char NameOfWave[50];
	char DeleteArraysAtDestruction;

	srTRadSect1D(float* In_pEx, float* In_pEz, char vs_x_or_z, long In_icOtherCoord, srTSRWRadStructAccessData* pRadAccessData)
	{
		pEx = In_pEx; pEz = In_pEz;

		icOtherCoord = In_icOtherCoord;
		VsXorZ = vs_x_or_z;

		if(vs_x_or_z == 'x')
		{
			ArgStep = pRadAccessData->xStep;
			ArgStart = pRadAccessData->xStart;
			np = pRadAccessData->nx;

			Robs = pRadAccessData->RobsX;
			RobsAbsErr = pRadAccessData->RobsXAbsErr;
			cArg = pRadAccessData->xc;

			WfrMin = pRadAccessData->xWfrMin;
			WfrMax = pRadAccessData->xWfrMax;
			WfrEdgeCorrShouldBeDone = pRadAccessData->WfrEdgeCorrShouldBeDone;

			Pres = pRadAccessData->Pres;
			LengthUnit = pRadAccessData->LengthUnit;
			PhotEnergyUnit = pRadAccessData->PhotEnergyUnit;

			OtherCoordVal = icOtherCoord*pRadAccessData->zStep + pRadAccessData->zStart;
		}
		else
		{
			ArgStep = pRadAccessData->zStep;
			ArgStart = pRadAccessData->zStart;
			np = pRadAccessData->nz;

			Robs = pRadAccessData->RobsZ;
			RobsAbsErr = pRadAccessData->RobsZAbsErr;
			cArg = pRadAccessData->zc;

			WfrMin = pRadAccessData->zWfrMin;
			WfrMax = pRadAccessData->zWfrMax;
			WfrEdgeCorrShouldBeDone = pRadAccessData->WfrEdgeCorrShouldBeDone;

			Pres = pRadAccessData->Pres;
			LengthUnit = pRadAccessData->LengthUnit;
			PhotEnergyUnit = pRadAccessData->PhotEnergyUnit;

			OtherCoordVal = icOtherCoord*pRadAccessData->xStep + pRadAccessData->xStart;
		}

		eVal = pRadAccessData->eStart;
		ThereIsUnderSamplingIn2D = 0;

		DeleteArraysAtDestruction = 0;

		// Add new srTSRWRadStructAccessData elements, if any !!!
	}

	srTRadSect1D() 
	{
		pEx = pEz = 0;
		ThereIsUnderSamplingIn2D = 0;
		DeleteArraysAtDestruction = 0;
	}

	~srTRadSect1D()
	{
		if(DeleteArraysAtDestruction) DeleteArrays();
	}

	int SetupDupl(srTRadSect1D& RadSect1D)
	{
		if(RadSect1D.DeleteArraysAtDestruction) RadSect1D.DeleteArrays();

		RadSect1D = *this;

		long Two_np = np << 1;
		RadSect1D.pEx = new float[Two_np];
		if(RadSect1D.pEx == 0) return MEMORY_ALLOCATION_FAILURE;
		RadSect1D.pEz = new float[Two_np];
		if(RadSect1D.pEz == 0) return MEMORY_ALLOCATION_FAILURE;

		float *tEx = pEx, *tExD = RadSect1D.pEx;
		float *tEz = pEz, *tEzD = RadSect1D.pEz;
		for(int i=0; i<Two_np; i++)
		{
			*(tExD++) = *(tEx++); *(tEzD++) = *(tEz++);
		}

		RadSect1D.DeleteArraysAtDestruction = 1;
		return 0;
	}

	void DeleteArrays()
	{
		if(pEx != 0) delete[] pEx; pEx = 0;
		if(pEz != 0) delete[] pEz; pEz = 0;
	}

	void SetNonZeroWavefrontLimitsToFullRange()
	{// This switches off the sharp Wfr edges treatment
		WfrMin = ArgStart;
		WfrMax = ArgStart + np*ArgStep;
	}
	void PreserveLogicsOfWfrLimitsAtRangeResizing(srTRadSect1D* pOldRadData)
	{
		double DistTol = 0.01*ArgStep;
		if((::fabs(pOldRadData->WfrMin - pOldRadData->ArgStart) < DistTol) && (::fabs(pOldRadData->ArgStart + pOldRadData->np*pOldRadData->ArgStep - pOldRadData->WfrMax) < DistTol))
		{
			WfrMin = ArgStart; WfrMax = ArgStart + np*ArgStep;
		}
		else
		{
			WfrMin = pOldRadData->WfrMin; WfrMax = pOldRadData->WfrMax;
		}
	}
	void FindMinMaxReE(srTMinMaxEParam& a)
	{
		float &MaxReEx = a.MaxReEx, &MaxImEx = a.MaxImEx, &MaxReEz = a.MaxReEz, &MaxImEz = a.MaxImEz, &MinReEx = a.MinReEx, &MinImEx = a.MinImEx, &MinReEz = a.MinReEz, &MinImEz = a.MinImEz;
		long &xIndMaxReEx = a.xIndMaxReEx, &xIndMaxImEx = a.xIndMaxImEx, &xIndMaxReEz = a.xIndMaxReEz, &xIndMaxImEz = a.xIndMaxImEz, &xIndMinReEx = a.xIndMinReEx, &xIndMinImEx = a.xIndMinImEx, &xIndMinReEz = a.xIndMinReEz, &xIndMinImEz = a.xIndMinImEz;
		MaxReEx = MaxImEx = MaxReEz = MaxImEz = (float)(-1.E+23);
		MinReEx = MinImEx = MinReEz = MinImEz = (float)(1.E+23);
		float *tReEx = pEx, *tImEx = pEx + 1, *tReEz = pEz, *tImEz = pEz + 1;
		for(long i=0; i<np; i++)
		{
			if(*tReEx > MaxReEx)
			{
				MaxReEx = *tReEx; xIndMaxReEx = i;
			}
			if(*tImEx > MaxImEx)
			{
				MaxImEx = *tImEx; xIndMaxImEx = i;
			}
			if(*tReEz > MaxReEz)
			{
				MaxReEz = *tReEz; xIndMaxReEz = i;
			}
			if(*tImEz > MaxImEz)
			{
				MaxImEz = *tImEz; xIndMaxImEz = i;
			}
			if(*tReEx < MinReEx)
			{
				MinReEx = *tReEx; xIndMinReEx = i;
			}
			if(*tImEx < MinImEx)
			{
				MinImEx = *tImEx; xIndMinImEx = i;
			}
			if(*tReEz < MinReEz)
			{
				MinReEz = *tReEz; xIndMinReEz = i;
			}
			if(*tImEz < MinImEz)
			{
				MinImEz = *tImEz; xIndMinImEz = i;
			}
			tReEx += 2; tImEx += 2; tReEz += 2; tImEz += 2;
		}
	}
	void FindWfrEdgeIndexes(long& iWfrMin, long& iWfrMax)
	{
		long fWfrMin = (long)((WfrMin - ArgStart)/ArgStep);
		long Abs_fWfrMin = (long)::fabs((long double)fWfrMin);
		iWfrMin = ((fWfrMin - Abs_fWfrMin) > 1.E-03*ArgStep)? (Abs_fWfrMin + 1) : Abs_fWfrMin;

		long fWfrMax = (long)((WfrMax - ArgStart)/ArgStep);
		long Abs_fWfrMax = (long)::fabs((long double)fWfrMax);
		iWfrMax = ((fWfrMax - Abs_fWfrMax) > 1.E-03*ArgStep)? Abs_fWfrMax : (Abs_fWfrMax + 1);
	}
	char DetermineMoreIntensivePolarComp()
	{
		if(pEx == 0)
		{
			if(pEz != 0) return 'z';
			else return 0;
		}
		if(pEz == 0)
		{
			if(pEx != 0) return 'x';
			else return 0;
		}

		char cOut = 'x';
		float Imax = 0.;
		float *tEx = pEx, *tEz = pEz;
		for(long i=0; i<np; i++)
		{
			float Ix = (*tEx)*(*tEx) + (*(tEx+1))*(*(tEx+1)); tEx += 2;
			float Iz = (*tEz)*(*tEz) + (*(tEz+1))*(*(tEz+1)); tEz += 2;
			if(Ix > Iz) if(Ix > Imax) { Imax = Ix; cOut = 'x';}
			else if(Iz > Ix) if(Iz > Imax) { Imax = Iz; cOut = 'z';}
		}
		return cOut;
	}
};

//*************************************************************************

class srTStokesStructAccessData : public CGenObject {
public:

	float *pBaseSto;
	float *pS0, *pS1, *pS2, *pS3; //SRWLib
	waveHndl wSto;
	int hStateSto;
	bool MemoryWasAllocatedInternally;

	double eStep, eStart, xStep, xStart, zStep, zStart, yStep, yStart;
	long ne, nx, nz, ny;

	double dx, dz; //for flux calculations, for internal use only

	long AuxLong1, AuxLong2, AuxLong3, AuxLong4;

	srTStokesStructAccessData()
	{
		Initialize();
	}
	srTStokesStructAccessData(double In_eStart, double In_eStep, int In_ne, double In_xStart, double In_xStep, int In_nx, double In_zStart, double In_zStep, int In_nz, bool MemoryShouldBeAllocated)
	{
		Initialize();

		eStep = In_eStep; eStart = In_eStart; xStep = In_xStep; xStart = In_xStart; zStep = In_zStep; zStart = In_zStart;
		ne = In_ne; nx = In_nx; nz = In_nz;
		ny = 1; yStep = 0; yStart = 0;

		if(MemoryShouldBeAllocated)
		{
			//long Np = In_ne*In_nx*In_nz;
			long long Np = ((long long)In_ne)*((long long)In_nx)*((long long)In_nz);
			if(Np <= 0) return;
			
			pBaseSto = new float[Np << 2];
			MemoryWasAllocatedInternally = true;
		}
	}
	srTStokesStructAccessData(double In_eStart, double In_eStep, int In_ne, double In_xStart, double In_xStep, int In_nx, double In_zStart, double In_zStep, int In_nz, double In_yStart, double In_yStep, int In_ny, bool MemoryShouldBeAllocated)
	{
		Initialize();

		eStep = In_eStep; eStart = In_eStart; xStep = In_xStep; xStart = In_xStart; zStep = In_zStep; zStart = In_zStart; yStep = In_yStep; yStart = In_yStart;
		ne = In_ne; nx = In_nx; nz = In_nz; ny = In_ny;

		if(MemoryShouldBeAllocated)
		{
			//long Np = In_ne*In_nx*In_nz*In_ny;
			long long Np = ((long long)In_ne)*((long long)In_nx)*((long long)In_nz)*((long long)In_ny);
			if(Np <= 0) return;
			
			pBaseSto = new float[Np << 2];
			MemoryWasAllocatedInternally = true;
		}
	}
    srTStokesStructAccessData(srTWfrSmp* pWfrSmp)
	{
		Initialize();

		if(pWfrSmp != 0)
		{
			eStart = pWfrSmp->LambStart;
            ne = pWfrSmp->nLamb;
			eStep = (ne > 1)? (pWfrSmp->LambEnd - eStart)/(ne - 1) : 0;
	
			xStart = pWfrSmp->xStart; 
			nx = pWfrSmp->nx; 
			xStep = (nx > 1)? (pWfrSmp->xEnd - xStart)/(nx - 1) : 0;
				
			zStart = pWfrSmp->zStart; 
			nz = pWfrSmp->nz; 
			zStep = (nz > 1)? (pWfrSmp->zEnd - zStart)/(nz - 1) : 0;

			yStart = pWfrSmp->yStart; 
			ny = pWfrSmp->ny; 
			yStep = (ny > 1)? (pWfrSmp->yEnd - yStart)/(ny - 1) : 0;

			//long Np = ne*nx*nz*ny;
			long long Np = ((long long)ne)*((long long)nx)*((long long)nz)*((long long)ny);
			if(Np <= 0) throw INCORRECT_GRID_FOR_WAVEFRONT;
			pBaseSto = new float[Np << 2];
			MemoryWasAllocatedInternally = true;
		}
	}
    srTStokesStructAccessData(srTSRWStokesInData* pStokesInData)
	{
		Initialize();

		if(pStokesInData != 0)
		{
			eStart = pStokesInData->eStart;
			eStep = pStokesInData->eStep;
            ne = pStokesInData->ne;

			xStart = pStokesInData->xStart;
			xStep = pStokesInData->xStep;
			nx = pStokesInData->nx;
				
			zStart = pStokesInData->zStart;
			zStep = pStokesInData->xStep;
			nz = pStokesInData->nz;
			
 			yStart = pStokesInData->yStart;
			yStep = pStokesInData->yStep;
			ny = pStokesInData->ny;

            pBaseSto = pStokesInData->pBaseSto;
            wSto = pStokesInData->wSto;
            hStateSto = pStokesInData->hStateSto;
			MemoryWasAllocatedInternally = false;
		}
	}

	~srTStokesStructAccessData()
	{
		if(MemoryWasAllocatedInternally)
		{
			if(pBaseSto != 0) delete[] pBaseSto;
		}
	}

	void Initialize()
	{
		pBaseSto = 0;
		pS0 = pS1 = pS2 = pS3 = 0;
		wSto = 0;
		hStateSto = -1;
		AuxLong1 = AuxLong2 = AuxLong3 = AuxLong4 = 0;

		eStep = eStart = xStep = xStart = zStep = zStart = yStep = yStart = 0.;
		ne = nx = nz = ny = 0;
		dx = dz = 0;

		MemoryWasAllocatedInternally = false;
	}

	void ZeroStokesData()
	{
		if(pBaseSto == 0) return;
		//long LenData = (ne << 2)*nx*nz;
		//long LenData = (ne << 2)*nx*nz*ny;
		long long LenData = (ne << 2)*((long long)nx)*((long long)nz)*((long long)ny);
		if(LenData <= 0) return;
		float *tData = pBaseSto;
		//for(long i=0; i<LenData; i++) *(tData++) = 0.;
		for(long long i=0; i<LenData; i++) *(tData++) = 0.;
	}

	void OutDataPtrs(srTSRWStokesInData* pStokesInData)
	{
		if(pStokesInData == 0) return;

		pStokesInData->eStart = eStart;
        pStokesInData->eStep = eStep;
        pStokesInData->ne = ne;

		pStokesInData->xStart = xStart;
		pStokesInData->xStep = xStep;
        pStokesInData->nx = nx;
				
		pStokesInData->zStart = zStart;
		pStokesInData->xStep = zStep;
        pStokesInData->nz = nz;
		
		pStokesInData->yStart = yStart;
		pStokesInData->yStep = yStep;
        pStokesInData->ny = ny;

        pStokesInData->pBaseSto = pBaseSto;
        pStokesInData->wSto = wSto;
        pStokesInData->hStateSto = hStateSto;
	}

	void UpdateLongitudinalGridParams(double In_yStart, double In_yStep, int In_ny)
	{
        yStep = In_yStep; yStart = In_yStart;
        ny = In_ny;
	}
	void UpdateSlitDimensions(double In_dx, double In_dz)
	{
        dx = In_dx; dz = In_dz;
	}

	double xFinMin()
	{
		return xStart + (nx - 1)*xStep;
	}
	double zFinMin()
	{
		return zStart + (nz - 1)*zStep;
	}
	double eFinMin()
	{
		return eStart + (ne - 1)*eStep;
	}
};

//*************************************************************************

class srTPowDensStructAccessData : public CGenObject {
public:

	float *pBasePowDens;
	waveHndl wPowDens;
	int hStatePowDens;
	bool BasePowDensWasEmulated;

	double xStep, xStart, zStep, zStart;
	long nx, nz;

	CSmartPtr<double> m_spObSurfData;

	srTPowDensStructAccessData(srTWfrSmp* pWfrSmp)
	{
		if(pWfrSmp == 0) throw INCORRECT_PARAMS_SR_COMP;

		nx = pWfrSmp->nx;
		nz = pWfrSmp->nz;
		xStart = pWfrSmp->xStart;
		zStart = pWfrSmp->zStart;
		xStep = zStep = 0.;
		if(pWfrSmp->nx > 1) xStep = (pWfrSmp->xEnd - pWfrSmp->xStart)/(pWfrSmp->nx - 1);
		if(pWfrSmp->nz > 1) zStep = (pWfrSmp->zEnd - pWfrSmp->zStart)/(pWfrSmp->nz - 1);

		//long LenData = nx*nz;
		long long LenData = ((long long)nx)*((long long)nz);
		if(LenData <= 0) throw INCORRECT_PARAMS_SR_COMP;
		pBasePowDens = 0;
		pBasePowDens = new float[LenData];
		if(pBasePowDens == 0) throw MEMORY_ALLOCATION_FAILURE;
		BasePowDensWasEmulated = true;
	}
	srTPowDensStructAccessData(srTSRWPowDensInData* pPowDensIn)
	{
		if(pPowDensIn == 0) throw INCORRECT_PARAMS_SR_COMP;

		nx = pPowDensIn->nx;
		nz = pPowDensIn->nz;
		xStart = pPowDensIn->xStart;
		zStart = pPowDensIn->zStart;
		xStep = pPowDensIn->xStep;
        zStep = pPowDensIn->zStep;

		pBasePowDens = pPowDensIn->pBasePowDens;
        wPowDens = pPowDensIn->wPowDens;
        hStatePowDens = pPowDensIn->hStatePowDens;
        BasePowDensWasEmulated = false;
	}
	srTPowDensStructAccessData() 
	{
		pBasePowDens = 0;
		hStatePowDens = 0;
		BasePowDensWasEmulated = false;
	}
	~srTPowDensStructAccessData()
	{
		if(BasePowDensWasEmulated && (pBasePowDens != 0)) delete[] pBasePowDens;
	}

	void EnsureNonNegativeValues()
	{
		if(pBasePowDens == 0) return;
        float *tData = pBasePowDens;
		for(long i=0; i<nz; i++)
		{
			for(long j=0; j<nx; j++)
			{
				if(*tData < 0) *tData = 0;
				tData++;
			}
		}
	}

	void SetRadSamplingFromObs(srTWfrSmp& DistrInfoDat)
	{
		xStep = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
		nx = DistrInfoDat.nx;
		
		zStep = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
		nz = DistrInfoDat.nz;
	
		if(DistrInfoDat.ObsPlaneIsTransverse())
		{
			xStart = DistrInfoDat.xStart;
			zStart = DistrInfoDat.zStart;
		}
		else
		{
			xStart = -0.5*xStep*(nx - 1); //local mesh with respect to center point in arbitrarily-oriented observation plane
			zStart = -0.5*zStep*(nz - 1);
		}

		// To walk around a bug in Igor
		if(xStep == 0.) { xStep = (xStart != 0.)? (1.e-08)*(::fabs(xStart)) : 1.e-10;}
		if(zStep == 0.) { zStep = (zStart != 0.)? (1.e-08)*(::fabs(zStart)) : 1.e-10;}

		m_spObSurfData = DistrInfoDat.m_spSurfData;
	}
};

//*************************************************************************

struct srTElecBeamMoments {
	DOUBLE E, Mx, Mxp, Mz, Mzp;
	DOUBLE Mee;
	DOUBLE Mxx, Mxxp, Mxpxp, Mzz, Mzzp, Mzpzp; // Central Moments
	DOUBLE Mxz, Mxpz, Mxzp, Mxpzp;

	srTElecBeamMoments(DOUBLE* pElecBeamData)
	{
		if(pElecBeamData != 0)
		{
			E = *pElecBeamData;
			Mx = *(pElecBeamData + 2);
			Mxp = *(pElecBeamData + 3);
			Mz = *(pElecBeamData + 4);
			Mzp = *(pElecBeamData + 5);
			
			Mee	= *(pElecBeamData + 13);
			
			Mxx = *(pElecBeamData + 20);
			Mxxp = *(pElecBeamData + 21);
			Mxpxp = *(pElecBeamData + 22);
			Mzz = *(pElecBeamData + 23);
			Mzzp = *(pElecBeamData + 24);
			Mzpzp = *(pElecBeamData + 25);
			
			Mxz = *(pElecBeamData + 26);
			Mxpz = *(pElecBeamData + 27);
			Mxzp = *(pElecBeamData + 28);
			Mxpzp = *(pElecBeamData + 29);
		}
		else ZeroAll();
	}

	srTElecBeamMoments(srTEbmDat* pEbm)
	{
		ZeroAll();
		if(pEbm == 0) return;

		E = pEbm->Energy;

		Mx = pEbm->x0;
        Mxp = pEbm->dxds0;
        Mz = pEbm->z0;
		Mzp = pEbm->dzds0;

		Mee	= pEbm->Mee;
        Mxx = pEbm->Mxx;
		Mxxp = pEbm->Mxxp;
		Mxpxp = pEbm->Mxpxp;
		Mzz = pEbm->Mzz;
		Mzzp = pEbm->Mzzp;
		Mzpzp = pEbm->Mzpzp;
	}

	srTElecBeamMoments() { ZeroAll();}

	void ZeroAll()
	{
		E = Mx = Mxp = Mz = Mzp = 0.;
		Mee = 0.;
		Mxx = Mxxp = Mxpxp = Mzz = Mzzp = Mzpzp = 0.; // Central Moments
		Mxz = Mxpz = Mxzp = Mxpzp = 0.;
	}

	void OutData(srTEbmDat* pEbm)
	{
		if(pEbm == 0) return;
		
		pEbm->Energy = E;

		pEbm->x0 = Mx;
        pEbm->dxds0 = Mxp;
        pEbm->z0 = Mz;
		pEbm->dzds0 = Mzp;

		pEbm->Mee = Mee;
        pEbm->Mxx = Mxx;
		pEbm->Mxxp = Mxxp;
		pEbm->Mxpxp = Mxpxp;
		pEbm->Mzz = Mzz;
		pEbm->Mzzp = Mzzp;
		pEbm->Mzpzp = Mzpzp;
	}
};

//*************************************************************************

struct srTWaveAccessData {
	char* pWaveData;
	char WaveType[2]; // 'f'|'d'|'cf'|'cd'
	long AmOfDims;
	long DimSizes[10];
	double DimStartValues[10];
	double DimSteps[10];
	char DimUnits[10][255];
	char DataUnits[255];

	waveHndl wHndl;
	int hState;

	char NameOfWave[50];

	srTWaveAccessData()
	{
		Init(); //OC13112018

		//pWaveData = 0;
		//*WaveType = '\0';
		//AmOfDims = -1;
		//for(int i=0; i<10; i++) 
		//{
		//	DimSizes[i] = -1;
		//	DimStartValues[i] = 1.E+23;
		//	DimSteps[i] = 1.E+23;
		//	*(DimUnits[i]) = '\0';
		//}
		//*NameOfWave = '\0';
		//*DataUnits = '\0';
	}

	srTWaveAccessData(char* pcData, char typeData, SRWLRadMesh* pMesh)
	{//OC13112018
		Init();

		pWaveData = pcData;
		WaveType[0] = typeData; WaveType[1] = '\0';
	
		int nDims = 0;
		int n1 = 0, n2 = 0, n3 = 0;
		double start1 = 0, start2 = 0, start3 = 0;
		double step1 = 0, step2 = 0, step3 = 0;
		if(pMesh->ne > 1) 
		{
			nDims++; 
			n1 = pMesh->ne;
			start1 = pMesh->eStart;
			step1 = (pMesh->eFin - start1)/(n1 - 1);
		}
		if(pMesh->nx > 1) 
		{
			nDims++;
			if(n1 == 0) 
			{
				n1 = pMesh->nx;
				start1 = pMesh->xStart;
				step1 = (pMesh->xFin - start1)/(n1 - 1);
			}
			else 
			{
				n2 = pMesh->nx;
				start2 = pMesh->xStart;
				step2 = (pMesh->xFin - start2)/(n2 - 1);
			}
		}
		if(pMesh->ny > 1) 
		{
			nDims++;
			if(n1 == 0) 
			{
				n1 = pMesh->ny;
				start1 = pMesh->yStart;
				step1 = (pMesh->yFin - start1)/(n1 - 1);
			}
			else if(n2 == 0) 
			{
				n2 = pMesh->ny;
				start2 = pMesh->yStart;
				step2 = (pMesh->yFin - start2)/(n2 - 1);
			}
			else 
			{
				n3 = pMesh->ny;
				start3 = pMesh->yStart;
				step3 = (pMesh->yFin - start3)/(n3 - 1);
			}
		}

		AmOfDims = nDims;

		DimSizes[0] = n1;
		DimSizes[1] = n2;
		DimSizes[2] = n3;
		DimStartValues[0] = start1;
		DimStartValues[1] = start2;
		DimStartValues[2] = start3;
		DimSteps[0] = step1;
		DimSteps[1] = step2;
		DimSteps[2] = step3;

		//To process Mutual Intensity case: pMesh->type == 'm' !
	}

	void Init()
	{//OC13112018
		pWaveData = 0;
		*WaveType = '\0';
		AmOfDims = -1;
		for(int i=0; i<10; i++) 
		{
			DimSizes[i] = -1;
			DimStartValues[i] = 1.E+23;
			DimSteps[i] = 1.E+23;
			*(DimUnits[i]) = '\0';
		}
		*NameOfWave = '\0';
		*DataUnits = '\0';
	}

	void OutRealData(double* ArrToFill, long MaxLen)
	{
		if((ArrToFill == 0) || (pWaveData == 0) || (AmOfDims == 0)) return;

		long ActLen = 1;
		for(int i=0; i<AmOfDims; i++) ActLen *= DimSizes[i];
		if(ActLen < MaxLen) MaxLen = ActLen;
		
		float* fpData = 0;
		double* dpData = 0;
		int Per = 1;
		if(WaveType[0] == 'f') fpData = (float*)pWaveData;
		else if(WaveType[0] == 'd') dpData = (double*)pWaveData;
		else if(WaveType[0] == 'c') 
		{ 
			Per = 2;
			if(WaveType[1] == 'f') fpData = (float*)pWaveData; 
            else if(WaveType[1] == 'd') dpData = (double*)pWaveData;
		}
		else return;

		long IndCount = 0;
		//double CurVal;
		for(int j=0; j<AmOfDims; j++) 
		{
			if(fpData != 0) 
			{
				for(int k=0; k<DimSizes[j]; k++) { ArrToFill[IndCount++] = *fpData; fpData += Per;}
			}
			else if(dpData != 0) 
			{
				for(int k=0; k<DimSizes[j]; k++) { ArrToFill[IndCount++] = *dpData; dpData += Per;}
			}
		}
	}
};

//*************************************************************************

struct srTWaveAccessDataD1D {

	DOUBLE* pData;
	//long np;
	long long np;
	double Start;
	double Step;
	char ArgUnits[255];
	char DataUnits[255];
	char NameOfWave[50];

	waveHndl wHndl;
	int hState;

// Buf Vars
	double InvStep;

	srTWaveAccessDataD1D()
	{
		np = -1;
		Start = 1.E+23;
		Step = 1.E+23;
		*ArgUnits = '\0';
		*DataUnits = '\0';
		*NameOfWave = '\0';

		pData = 0; // This is checked !
		wHndl = NIL; // This is checked at Finishing Working with ...
	}

	void SetupBufVars()
	{
		InvStep = 1./Step;
	}
};

//*************************************************************************

struct srTRadExtract {

	int PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
	int Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E); 4: 1-e Flux; 5: Multi-e Flux; 6- Im(E); 7- Time or Photon Energy Integrated Intensity
	int PlotType; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
	int TransvPres; // 0: Spatial; 1: Angular

	double ePh, x, z;

	float* pExtractedData;
	DOUBLE* pExtractedDataD;

	waveHndl wExtractedData;
	int hStateExtractedData;

	srTRadExtract(int In_PolarizCompon, int In_Int_or_Phase, int In_SectID, int In_TransvPres, double In_e, double In_x, double In_z, char* In_pData)
	{
        PolarizCompon = In_PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
		Int_or_Phase = In_Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E)
		PlotType = In_SectID; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
		TransvPres = In_TransvPres; // 0: Spatial; 1: Angular
		ePh = In_e; x = In_x; z = In_z;

		pExtractedData = 0; pExtractedDataD = 0;

		if(In_Int_or_Phase != 2) pExtractedData = (float*)In_pData;
		else pExtractedDataD = (DOUBLE*)In_pData;
	}
	srTRadExtract() {};

	void SetupExtractedWaveAccessData(srTWaveAccessData* pWaveAccessData)
	{
		pWaveAccessData->pWaveData = (char*)(pExtractedData);
		pWaveAccessData->wHndl = wExtractedData;
		pWaveAccessData->hState = hStateExtractedData;
	}
};

//*************************************************************************

struct srTInterpolAux01 {

	double cAx0z1, cAx0z2, cAx0z3, cAx1z0, cAx1z1, cAx1z2, cAx1z3;
	double cAx2z0, cAx2z1, cAx2z2, cAx2z3, cAx3z0, cAx3z1, cAx3z2, cAx3z3;
	double cLAx1z0, cLAx0z1, cLAx1z1;

	srTInterpolAux01()
	{
		cAx0z1 = 0.1666666667;
		cAx0z2 = 0.5;
		cAx0z3 = 0.1666666667;
		cAx1z0 = 0.1666666667;
		cAx1z1 = 0.027777777778;
		cAx1z2 = 0.083333333333;
		cAx1z3 = 0.027777777778;
		cAx2z0 = 0.5;
		cAx2z1 = 0.083333333333;
		cAx2z2 = 0.25;
		cAx2z3 = 0.083333333333;
		cAx3z0 = 0.16666666667;
		cAx3z1 = 0.027777777778;
		cAx3z2 = 0.083333333333;
		cAx3z3 = 0.027777777778;

		cLAx1z0 = 1.;
		cLAx0z1 = 1.;
		cLAx1z1 = 1.;
	}
};

//*************************************************************************

struct srTInterpolAux01_1D {

	double cA1, cA2, cA3;
	double cLA1;

	srTInterpolAux01_1D()
	{
		cA1 = 0.1666666667;
		cA2 = 0.5;
		cA3 = 0.1666666667;

		cLA1 = 1.;
	}
};

//*************************************************************************

struct srTInterpolAux02 {

	double Ax0z0, Ax0z1, Ax0z2, Ax0z3, Ax1z0, Ax1z1, Ax1z2, Ax1z3;
	double Ax2z0, Ax2z1, Ax2z2, Ax2z3, Ax3z0, Ax3z1, Ax3z2, Ax3z3;
};

//*************************************************************************

struct srTInterpolAux02_1D {

	double A0, A1, A2, A3;
};

//*************************************************************************

struct srTInterpolAuxF {

	float f00, f10, f20, f30;
	float f01, f11, f21, f31;
	float f02, f12, f22, f32;
	float f03, f13, f23, f33;

	float fAvg, fNorm;
	void SetUpAvg()
	{
		fAvg = (float)(0.0625*(f00 + f10 + f20 + f30 + f01 + f11 + f21 + f31 + f02 + f12 + f22 + f32 + f03 + f13 + f23 + f33));
	}
	void NormalizeByAvg()
	{
		const float CritNorm = 1.;
		if(::fabs(fAvg) > CritNorm)
		{
			float a = (float)(1./fAvg);
			f00 *= a; f10 *= a; f20 *= a; f30 *= a;
			f01 *= a; f11 *= a; f21 *= a; f31 *= a;
			f02 *= a; f12 *= a; f22 *= a; f32 *= a;
			f03 *= a; f13 *= a; f23 *= a; f33 *= a;
			fNorm = fAvg;
		}
		else fNorm = 1.;
	}
};

//*************************************************************************

struct srTInterpolAuxF_1D {

	float f0, f1, f2, f3;
	float fAvg, fNorm;

	void SetUpAvg()
	{
		fAvg = (float)(0.25*(f0 + f1 + f2 + f3));
	}
	void NormalizeByAvg()
	{
		const float CritNorm = 1.;
		if(::fabs(fAvg) > CritNorm)
		{
			float a = (float)(1./fAvg);
			f0 *= a; f1 *= a; f2 *= a; f3 *= a;

			fNorm = fAvg;
		}
		else fNorm = 1.;
	}
};

//*************************************************************************

struct srTMagFieldAccessData {

	waveHndl wBx, wBz;
	int hStateBx, hStateBz;
};

//*************************************************************************

struct srTDataPtrsForWfrEdgeCorr {

	float *ExpArrXSt, *ExpArrXFi;
	float *ExpArrZSt, *ExpArrZFi;

	float *FFTArrXStEx, *FFTArrXFiEx;
	float *FFTArrZStEx, *FFTArrZFiEx;
	float *FFTArrXStEz, *FFTArrXFiEz;
	float *FFTArrZStEz, *FFTArrZFiEz;

	float fxStzSt[4], fxFizSt[4], fxStzFi[4], fxFizFi[4];

	double dxSt, dxFi, dzSt, dzFi, dx, dz;
	char WasSetup;

	srTDataPtrsForWfrEdgeCorr()
	{
		InitializeAll();
	}

	void InitializeAll()
	{
		ExpArrXSt = ExpArrXFi = 0;
		ExpArrZSt = ExpArrZFi = 0;

		FFTArrXStEx = FFTArrXFiEx = 0;
		FFTArrZStEx = FFTArrZFiEx = 0;
		FFTArrXStEz = FFTArrXFiEz = 0;
		FFTArrZStEz = FFTArrZFiEz = 0;

		dxSt = dxFi = dzSt = dzFi = dx = dz = 0.;

		for(int i=0; i<4; i++)
		{
			fxStzSt[i] = 0.; fxFizSt[i] = 0.; fxStzFi[i] = 0.; fxFizFi[i] = 0.;
		}
		WasSetup = 0;
	}
	void DisposeData()
	{
		if(ExpArrXSt != 0) delete[] ExpArrXSt;
		if(ExpArrXFi != 0) delete[] ExpArrXFi;
		if(ExpArrZSt != 0) delete[] ExpArrZSt;
		if(ExpArrZFi != 0) delete[] ExpArrZFi;
		if(FFTArrXStEx != 0) delete[] FFTArrXStEx;
		if(FFTArrXFiEx != 0) delete[] FFTArrXFiEx;
		if(FFTArrZStEx != 0) delete[] FFTArrZStEx;
		if(FFTArrZFiEx != 0) delete[] FFTArrZFiEx;
		InitializeAll();
	}
};

//*************************************************************************

struct srTDataPtrsForWfrEdgeCorr1D {

	float *ExpArrSt, *ExpArrFi;
	float fSt[4], fFi[4];
	double dSt, dFi, d;

	char WasSetup;

	srTDataPtrsForWfrEdgeCorr1D()
	{
		InitializeAll();
	}

	void InitializeAll()
	{
		ExpArrSt = ExpArrFi = 0;
		dSt = dFi = d = 0.;
		for(int i=0; i<4; i++)
		{
			fSt[i] = 0.; fFi[i] = 0.;
		}
		WasSetup = 0;
	}

	void DisposeData()
	{
		if(ExpArrSt != 0) delete[] ExpArrSt;
		if(ExpArrFi != 0) delete[] ExpArrFi;

		InitializeAll();
	}
};

//*************************************************************************

struct srTFringeInfo {
	long AmOfFringes;
	double LeftPointsPerFringe, RightPointsPerFringe;
};

//*************************************************************************

struct srTPropagScenario1D {

	srTRadResize1D ResizeBefore, ResizeAfter;
	srTFringeInfo CurFringeInfo;
	srTPropagScenario1D()
	{
		ResizeBefore.pm = ResizeBefore.pd = ResizeAfter.pm = ResizeAfter.pd = 1.;
	}
};

//*************************************************************************

struct srTPredictedPropagData1D {

	double CoordMomRat, AngMomRat;
	char CoordMomRadCanBeModified, AngMomRatCanBeModified;

	double SigmaCoord, SigmaAng;

	srTPredictedPropagData1D()
	{
		CoordMomRadCanBeModified = AngMomRatCanBeModified = 1;
	}
};

//*************************************************************************

struct srTRelAndAbsTolerance {

	double RelTol, AbsTol;
	double AuxCoef1, AuxCoef2;

	srTRelAndAbsTolerance(double InRelTol =0., double InAbsTol =0.)
	{
		RelTol = InRelTol; AbsTol = InAbsTol;
	}
};

//*************************************************************************

struct srTAuxTestIntegValues {

	float DecrResVal, NormVal, IncrResVal, IncrRes2Val;
	float HalfDecrResVal, HalfNormVal, HalfIncrResVal, HalfIncrRes2Val;
	char DecrResValIsNeeded, NormValIsNeeded, IncrResValIsNeeded, IncrRes2ValIsNeeded;
	char TreatExOrEz;

	srTAuxTestIntegValues()
	{
		DecrResValIsNeeded = NormValIsNeeded = IncrResValIsNeeded = IncrRes2ValIsNeeded = 1;
		TreatExOrEz = 'x';
	}
};

//*************************************************************************

struct srTAuxTestValues {

	srTAuxTestIntegValues PrevTestIntegVal, CurrentTestIntegVal;

	double PrevParamValue, CurrentParamValue, SuggestedParamValue;
	char PrevParamValueIsGood, CurrentParamValueIsGood;
	char PrevResolCheckResult, CurrentResolCheckResult;
	char ExitLoop, ExitSuccessful;
	char NeedToChangeStrategy;

	srTAuxTestValues()
	{
		PrevTestIntegVal.DecrResVal = PrevTestIntegVal.NormVal = PrevTestIntegVal.IncrResVal = PrevTestIntegVal.IncrRes2Val = 0.;
		CurrentTestIntegVal.DecrResVal = CurrentTestIntegVal.NormVal = CurrentTestIntegVal.IncrResVal = CurrentTestIntegVal.IncrRes2Val = 0.;

		PrevParamValueIsGood = CurrentParamValueIsGood = 0;
		ExitLoop = 0; ExitSuccessful = 0;
		NeedToChangeStrategy = 0;
	}
};

//*************************************************************************

class srTCosAndSinComp {

	double HalfPI, PI, TwoPI, ThreePIdTwo, One_dTwoPI, One_dHalfPI; // Constants
	double a2c, a4c, a6c, a8c, a10c, a12c, a3s, a5s, a7s, a9s, a11s, a13s;

public:

	srTCosAndSinComp()
	{
		Initialize();
	}

	void Initialize()
	{
		HalfPI = 1.5707963267949; PI = 3.141592653590; TwoPI = 6.2831853071796; ThreePIdTwo = 4.7123889803847;
		One_dTwoPI = 0.1591549430919; One_dHalfPI = 0.636619772367581;
		a2c = -0.5; a4c = 0.041666666666667; a6c = -0.0013888888888889; a8c = 0.000024801587301587; a10c = -2.755731922E-07;
		a3s = -0.16666666666667; a5s = 0.0083333333333333; a7s = -0.0001984126984127; a9s = 2.755731922E-06; a11s = -2.505210839E-08;
	}
	void CosAndSin(double x, double& Cos, double& Sin)
	{
		x -= TwoPI*int(x*One_dTwoPI);
		if(x < 0.) x += TwoPI;

		char ChangeSign=0;
		if(x > ThreePIdTwo) x -= TwoPI;
		else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

		double xe2 = x*x;
		Cos = 1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c))));
		Sin = x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s)))));
		if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
	}
};

//*************************************************************************

struct srTransvLimits {
	double Xmin, Xmax, Ymin, Ymax;
	bool LimitsAreDefined;
	srTransvLimits() { LimitsAreDefined = false;}
	
	void Setup(double InXmin, double InDelX, double InYmin, double InDelY)
	{
		Xmin = InXmin; Xmax = Xmin + InDelX; Ymin = InYmin; Ymax = Ymin + InDelY;
		LimitsAreDefined = true;
	}
};

//*************************************************************************

struct srTRadSASE {
	double Power; // [W]
	double WaistDiam; // [m]
	double WaistLongPos; // [m]
	double PhotonEnergySim; // [eV]
};

//*************************************************************************

struct srTPrecSASE {
	long npart;
	double rmax0;
	long ncar;
	long nptr;
	long nscr;
	long nscz;
	char lbc;
	double delz;
	double zstop;
	char iorb;
	char UseElecDistr;
	char CreateElecDistr;

	char itdp;
	long nslice;
	double zsep;
	long ntail;

	double photEn_xlamds; //photon energy in [eV] to calc. xlamds

	int alignradf; // "ALIGNRADF" <>0 imported radfile is aligned to electron beam
	int offsetradf; 	

	char AllowAutoChoiceOfNxNzForPropagat;
	double NxNzOversamplingParam;

	//int iPower;
	//int iRadHorSize;
	//int iRadVertSize;
};

//*************************************************************************

struct srTControlAccessSASE {

	static const int cMaxNumHarm = 7;

	//float *pBasePower_vs_s;
	float *ar_pBasePower_vs_s[cMaxNumHarm];

	float *pBaseRadPhase_vs_s; 
	float *pBaseRadSize_vs_s; 
	float *pBaseBmSizeX_vs_s;
	float *pBaseBmSizeZ_vs_s;

	//float *pBasePower_vs_t;
	//float *pBaseEnergy_vs_s;
	//float *pBasePeakPower_vs_s;
	float *ar_pBasePower_vs_t[cMaxNumHarm];
	float *ar_pBaseEnergy_vs_s[cMaxNumHarm];
	float *ar_pBasePeakPower_vs_s[cMaxNumHarm];

	//float *pLocEnergy_vs_s; //for internal use (in C part) only
	float *ar_pLocEnergy_vs_s[cMaxNumHarm]; //for internal use (in C part) only

	//float *pBaseBunchFact1_vs_s;
	//float *pBaseBunchFact2_vs_s;
	//float *pBaseBunchFact3_vs_s;
	//float *pBaseBunchFact4_vs_s;
	//float *pBaseBunchFact5_vs_s;
	float *ar_pBaseBunchFact_vs_s[cMaxNumHarm];

	waveHndl wControl;
	//waveHndl wPower_vs_s;
	waveHndl ar_wPower_vs_s[cMaxNumHarm];
	waveHndl wRadPhase_vs_s; 
	waveHndl wRadSize_vs_s; 
	waveHndl wBmSizeX_vs_s;
	waveHndl wBmSizeZ_vs_s;
	//waveHndl wPower_vs_t;
	//waveHndl wPeakPower_vs_s;
	//waveHndl wEnergy_vs_s;
	waveHndl ar_wPower_vs_t[cMaxNumHarm];
	waveHndl ar_wPeakPower_vs_s[cMaxNumHarm];
	waveHndl ar_wEnergy_vs_s[cMaxNumHarm];

	//waveHndl wBunchFact1_vs_s;
	//waveHndl wBunchFact2_vs_s;
	//waveHndl wBunchFact3_vs_s;
	//waveHndl wBunchFact4_vs_s;
	//waveHndl wBunchFact5_vs_s;
	waveHndl ar_wBunchFact_vs_s[cMaxNumHarm];

	//int hStatePower_vs_s;
	int ar_hStatePower_vs_s[cMaxNumHarm];
	int hStateRadPhase_vs_s; 
	int hStateRadSize_vs_s; 
	int hStateBmSizeX_vs_s;
	int hStateBmSizeZ_vs_s;
	//int hStatePower_vs_t;
	//int hStateEnergy_vs_s;
	//int hStatePeakPower_vs_s;
	int ar_hStatePower_vs_t[cMaxNumHarm];
	int ar_hStateEnergy_vs_s[cMaxNumHarm];
	int ar_hStatePeakPower_vs_s[cMaxNumHarm];

	//int hStateBunchFact1_vs_s;
	//int hStateBunchFact2_vs_s;
	//int hStateBunchFact3_vs_s;
	//int hStateBunchFact4_vs_s;
	//int hStateBunchFact5_vs_s;
	int ar_hStateBunchFact_vs_s[cMaxNumHarm];

	long ns, nt;
	double sStart, sStep, tStart, tStep;

	srTControlAccessSASE()
	{
		Initialize();
	}

	~srTControlAccessSASE()
	{
		DeallocLocalCont();
	}

	void Initialize()
	{
		wControl = NIL;
		wRadPhase_vs_s = wRadSize_vs_s = wBmSizeX_vs_s = wBmSizeZ_vs_s = NIL;

		// This is checked at Finish... !!!
		pBaseRadPhase_vs_s = pBaseRadSize_vs_s = pBaseBmSizeX_vs_s = pBaseBmSizeZ_vs_s = 0;
		//pLocEnergy_vs_s = 0;
		//pBaseBunchFact_vs_s = 0;
		//pBaseBunchFact1_vs_s = pBaseBunchFact2_vs_s = pBaseBunchFact3_vs_s = pBaseBunchFact4_vs_s = pBaseBunchFact5_vs_s = 0;

		for(int i=0; i<cMaxNumHarm; i++)
		{
			ar_wPower_vs_s[i] = NIL;
			ar_wPower_vs_t[i] = NIL;
			ar_wEnergy_vs_s[i] = NIL;
			ar_wPeakPower_vs_s[i] = NIL;
			ar_wBunchFact_vs_s[i] = NIL;

			ar_pBasePower_vs_s[i] = 0;
			ar_pBasePower_vs_t[i] = 0; 
			ar_pBaseEnergy_vs_s[i] = 0;
			ar_pBasePeakPower_vs_s[i] = 0;
			ar_pBaseBunchFact_vs_s[i] = 0;
			ar_pLocEnergy_vs_s[i] = 0;
		}
	}
	void AllocLocalCont(int numHarm)
	{
		for(int i=0; i<numHarm; i++)
		{
			if(ar_pBaseEnergy_vs_s[i] != 0)
			{
				if(ar_pLocEnergy_vs_s[i] == 0) ar_pLocEnergy_vs_s[i] = new float[ns];
			}
		}
	}
	void DeallocLocalCont()
	{
		for(int i=0; i<cMaxNumHarm; i++)
		{
			if(ar_pLocEnergy_vs_s[i] != 0) { delete[] ar_pLocEnergy_vs_s[i]; ar_pLocEnergy_vs_s[i] = 0;}
		}
	}
};

//*************************************************************************

struct srTWigComSASE {

//typedef long int f2c_integer;
//typedef float f2c_real;
//typedef double f2c_doublereal;
//typedef struct { f2c_real r, i; } f2c_complex;
//typedef struct { f2c_doublereal r, i; } f2c_doublecomplex;
//	struct {
//		f2c_doublereal aw0, delaw, xlamd, fbess, wcoefz[3], awz[2502], awd, awdr, 
//			quadf, quadd, fbess0, qfdx, qfdy, fl, dl, drl, f1st, qx[2502], qy[2502], awerx[2502], awery[2502], xkx, xky;
//		f2c_integer iseed, iwityp, iertyp, iseed0;
//	} wigcom_;

	double aw0, delaw, xlamd, fbess, wcoefz[3], awd, awdr, quadf, quadd, fbess0, qfdx, qfdy, fl, dl, drl, f1st, xkx, xky;
	long iseed, iwityp, iertyp, iseed0;
	long simcom_nwig, simcom_nsec;

	double chic_bfield, chic_magl, chic_dril; //for HGHG computation/setup

	char UndulatorJustPassed; // Aux. hack, used at setting-up of this structure
	double GapLen; // Longitudinal gap bw sections 

	srTWigComSASE()
	{
		Initialize();
	}

	void Initialize()
	{
		xlamd = 0.; //this is checked; this means that the parameters were not set up
		fl = dl = drl = f1st = 0;
		quadf = quadd = qfdx = qfdy = 0;
		UndulatorJustPassed = 0;
		chic_bfield = chic_magl = chic_dril = 0;
	}

};

//*************************************************************************

class srTUtiDataMD {
public:

/**
	char* pData;
	char DataType[2]; // 'f'|'d'|'cf'|'cd'
	long AmOfDims;
	long DimSizes[10];
	double DimStartValues[10];
	double DimSteps[10];
	char DimUnits[10][255];
	char DataUnits[255];
	char DataName[255];
**/

	static double ExtractValue(srTDataMD* pDataMD, long FlatInd, double* ImPart)
	{
		if((pDataMD == 0) || (FlatInd < 0)) return 0;
		if(pDataMD->pData == 0) return 0;
		if(ImPart != 0) *ImPart = 0;

		if(pDataMD->DataType[0] == 'f')
		{
			float *p0 = (float*)(pDataMD->pData);
			return (double)(*(p0 + FlatInd));
		}
		else if(pDataMD->DataType[0] == 'd')
		{
			double *p0 = (double*)(pDataMD->pData);
			return *(p0 + FlatInd);
		}
		else if(pDataMD->DataType[0] == 'c')
		{
			if(pDataMD->DataType[1] == 'f')
			{
				float *p0 = (float*)(pDataMD->pData);
				if(ImPart != 0) *ImPart = (double)(*(p0 + FlatInd + 1));
				return (double)(*(p0 + FlatInd));
			}
			else if(pDataMD->DataType[1] == 'd')
			{
				double *p0 = (double*)(pDataMD->pData);
				if(ImPart != 0) *ImPart = *(p0 + FlatInd + 1);
				return *(p0 + FlatInd);
			}
		}
		return 0;
	}

	static int ExtractDimNum(srTDataMD* pDataMD)
	{
		if(pDataMD == 0) return 0;
		if(pDataMD->pData == 0) return 0;
		int OutAmOfDims = (int)(pDataMD->AmOfDims);
		if(OutAmOfDims > 10) OutAmOfDims = 10;
		return OutAmOfDims;
	}

	static int ExtractDimSizes(srTDataMD* pDataMD, long* ArDimSizes)
	{
		if((pDataMD == 0) || (ArDimSizes == 0)) return 0;
		if(pDataMD->pData == 0) return 0;
		int DimNum = (int)(pDataMD->AmOfDims);
		if(DimNum > 10) DimNum = 10;

		long *tArDimSizes = ArDimSizes, *tDimSizes = pDataMD->DimSizes;
		for(int i=0; i<DimNum; i++)
		{
			*(tArDimSizes++) = *(tDimSizes++);
		}
		return DimNum;
	}

	static void ExtractDataPointer(srTDataMD* pDataMD, float*& fpData, double*& dpData)
	{
		if(pDataMD == 0) return;
		if(pDataMD->pData == 0) return;

		fpData = 0; dpData = 0;

		if(pDataMD->DataType[0] == 'f')
		{
			fpData = (float*)(pDataMD->pData);
		}
		else if(pDataMD->DataType[0] == 'd')
		{
			dpData = (double*)(pDataMD->pData);
		}
		else if(pDataMD->DataType[0] == 'c')
		{
			if(pDataMD->DataType[1] == 'f')
			{
                fpData = (float*)(pDataMD->pData);
			}
			else if(pDataMD->DataType[1] == 'd')
			{
				dpData = (double*)(pDataMD->pData);
            }
		}
	}

	static double ExtractDimStartValue(srTDataMD* pDataMD, int DimNo)
	{
		if((pDataMD == 0) || (DimNo < 0) || (DimNo > 4)) return 0;
		return pDataMD->DimStartValues[DimNo];
	}

	static double ExtractDimStep(srTDataMD* pDataMD, int DimNo)
	{
		if((pDataMD == 0) || (DimNo < 0) || (DimNo > 4)) return 0;
		return pDataMD->DimSteps[DimNo];
	}

};

//*************************************************************************

#endif
