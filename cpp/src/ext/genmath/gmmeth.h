/************************************************************************//**
 * File: gmmeth.h
 * Description: Auxiliary mathematical utilities (header)
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMMETH_H
#define __GMMETH_H

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

#include "gmvect.h"
#include <math.h>
#include <complex>
#include "stlstart.h"

//*************************************************************************

inline double Max(double x1, double x2) { return (x1<x2)? x2 : x1;}

//*************************************************************************

inline double Abs(double x) { return (x<0)? -x : x;}

//*************************************************************************

inline double Norm(double x) { return (x<0)? -x : x;}

//*************************************************************************

inline double Norm(complex<double> z) { return Max( Abs(z.real()), Abs(z.imag()));}

//*************************************************************************

inline void MakeItLarge(double& x) { x = 1.E+23;}

//*************************************************************************

inline void MakeItLarge(complex<double>& z) { z = complex<double>(1.E+23, 1.E+23);}

//*************************************************************************

inline void MakeItZero(double& x) { x = 0.;}

//*************************************************************************

inline void MakeItZero(complex<double>& z) { z = complex<double>(0., 0.);}

//*************************************************************************

template <class T, class ValClass> void FormalOneFoldInteg(
	T* PtrT, void (T::*FunPtr)(double,ValClass*), 
	int LenVal, double* RelPrecAndLimitsArray, short* ElemCompNotFinished, ValClass** IntegVal)
{// This uses Newton method (n=3)
 // This function requies memory allocation from outside:
 //
 // double* RelPrecAndLimitsArray = new double[LenVal+2];
 //	ValClass* IntegVal[6];
 //	for(int j=0; j<6; j++) IntegVal[j] = new ValClass[LenVal];
 // short* ElemCompNotFinished = new short[LenVal];

	const double IntegWeight[] = {3./8., 9./8., 9./8., 3./4.};
	ValClass LargeEntity, ZeroEntity; 
	MakeItLarge(LargeEntity); MakeItZero(ZeroEntity);

	ValClass* LocSumArray = IntegVal[1];
	ValClass* GenSumArray = IntegVal[2];
	ValClass* InstIntegVal = IntegVal[3];
	ValClass* PrInstIntegVal = IntegVal[4];
	ValClass* ValF = IntegVal[5];

	for(int k=0; k<LenVal; k++)
	{
		PrInstIntegVal[k] = LargeEntity; 
		LocSumArray[k] = GenSumArray[k] = (IntegVal[0])[k] = ZeroEntity;
		ElemCompNotFinished[k] = 1;
	}

	double ArgStrt = RelPrecAndLimitsArray[LenVal];
	double ArgFin = RelPrecAndLimitsArray[LenVal+1];

	double StepArg, Arg;
	double IntervLength = ArgFin - ArgStrt;

	short IndForWeight, IndForPass;
	short NotFirstPass = 0;
	short MorePassesNeeded = 1;

	int AmOfPoi = 4;
	int AmOfPoi_mi_1;

	while(MorePassesNeeded)
	{
		AmOfPoi_mi_1 = AmOfPoi - 1; StepArg = IntervLength/AmOfPoi_mi_1;
		Arg = ArgStrt;
		IndForWeight = IndForPass = 0;

		for(int i=0; i<AmOfPoi; i++)
		{
			if(IndForPass==2) IndForPass = 0;
			if(IndForWeight==4) IndForWeight = 1;
			if(!(NotFirstPass && (IndForPass==0)))
			{
				if(i==AmOfPoi_mi_1) IndForWeight = 0;

				(PtrT->*FunPtr)(Arg, ValF);

				for(int kk=0; kk<LenVal; kk++)
					if(ElemCompNotFinished[kk]) LocSumArray[kk] += IntegWeight[IndForWeight]*ValF[kk];
			}
			IndForPass++; IndForWeight++; Arg += StepArg;
		}
		short MorePassesNeededProbe = 0;
		for(int kkk=0; kkk<LenVal; kkk++)
			if(ElemCompNotFinished[kkk])
			{
				GenSumArray[kkk] += LocSumArray[kkk];

				ValClass& InstVal = InstIntegVal[kkk]; 
				InstVal = StepArg*GenSumArray[kkk];
				ValClass& PrVal = PrInstIntegVal[kkk];

				if(Norm(InstVal - PrVal)/Norm(InstVal) < RelPrecAndLimitsArray[kkk]) ElemCompNotFinished[kkk] = 0;
				else
				{
					MorePassesNeededProbe = 1;
					PrVal = InstVal; LocSumArray[kkk] = ZeroEntity;
				}
			}
		MorePassesNeeded = MorePassesNeededProbe;

		AmOfPoi = AmOfPoi_mi_1*2 + 1;
		NotFirstPass = 1;
	}
	for(int kkkk=0; kkkk<LenVal; kkkk++) (IntegVal[0])[kkkk] = InstIntegVal[kkkk];
}

//*************************************************************************

class CGenMathMeth
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{

public:

	static bool VectCheckIfCollinear(double xV1, double yV1, double zV1, double xV2, double yV2, double zV2, double RelPrec);

	//static double Integ1D_FuncWithEdgeDer(double (*pF)(double), double x1, double x2, double dFdx1, double dFdx2, double RelPrec);
	static double Integ1D_FuncWithEdgeDer(double (*pF)(double, void*), double x1, double x2, double dFdx1, double dFdx2, double RelPrec, void* pv=0); //OC20112018
	static double Integ1D_Func(double (*pF)(double, void*), double x1, double x2, double RelPrec, void* pv=0); //OC02122018

	//static double Integ1D_FuncDefByArray(double* FuncArr, long Np, double Step);
	//static double Integ1D_FuncDefByArray(float* FuncArr, long Np, double Step);
	//template <class T> static double Integ1D_FuncDefByArray(T* FuncArr, long Np, double Step)
	template <class T> static double Integ1D_FuncDefByArray(T* FuncArr, long long Np, double Step)
	{
		if((FuncArr == 0) || (Np < 2) || (Step == 0)) return 0;
		//if(Np == 2) return (double)(0.5*(FuncArr[0] + FuncArr[1]));
		if(Np == 2) return (double)(0.5*Step*(FuncArr[0] + FuncArr[1])); //02122018

		T *tFuncArr = FuncArr + 1;
		bool NpIsEven = (Np == ((Np >> 1) << 1));

		//long NpSim = Np;
		long long NpSim = Np;
		if(NpIsEven) NpSim--;

		//Simpson part
		double Sum1 = 0, Sum2 = 0;
		//for(long i = 1; i < ((NpSim - 3) >> 1); i++)
		for(long long i = 1; i < ((NpSim - 3) >> 1); i++)
		{
			Sum1 += *(tFuncArr++);
			Sum2 += *(tFuncArr++);
		}
		Sum1 += *(tFuncArr++);
		double res = (Step/3.)*(FuncArr[0] + 4.*Sum1 + 2.*Sum2 + (*tFuncArr));

		//if(NpIsEven) res += (double)(0.5*(FuncArr[NpSim - 1] + FuncArr[NpSim])); //Last step is "trapethoidal"
		//OC291214
		if(NpIsEven) res += (double)(0.5*Step*(FuncArr[NpSim - 1] + FuncArr[NpSim])); //Last step is "trapethoidal"

		return res;
	}

	//template <class T> static double Integ2D_FuncDefByArray(T* arFlatFunc, long np1, long np2, double step1, double step2)
	template <class T> static double Integ2D_FuncDefByArray(T* arFlatFunc, long long np1, long long np2, double step1, double step2)
	{
		if((arFlatFunc == 0) || (np1 < 2) || (step1 == 0) || (np2 < 2) || (step2 == 0)) return 0;

		double *arResInt1D = new double[np2];
		double *t_arResInt1D = arResInt1D;
		T *p_arFlatFunc = arFlatFunc;
		//for(long i2 = 0; i2 < np2; i2++)
		for(long long i2 = 0; i2 < np2; i2++)
		{
			*(t_arResInt1D++) = Integ1D_FuncDefByArray(p_arFlatFunc, np1, step1);
			p_arFlatFunc += np1;
		}
		double res = Integ1D_FuncDefByArray(arResInt1D, np2, step2);

		if(arResInt1D != 0) delete[] arResInt1D;
		return res;
	}

	//template <class T> static T tabFunc2D(int ix, int iy, int nx, T* pF)
	template <class T> static T tabFunc2D(long long ix, long long iy, long long nx, T* pF)
	{//just function value
		if(pF == 0) return (T)0.;
		//long ofst = iy*nx + ix;
		long long ofst = iy*nx + ix;
		return *(pF + ofst);
	}

	template <class T> static T interpFunc2D(double x, double y, double xStart, double xEnd, int nx, double yStart, double yEnd, int ny, T* pF)
	{//uses bilinear interpolation
		if((pF == 0) || (nx <= 0) || (ny <= 0)) return (T)0.;

		double xStep = 0, yStep = 0;
		if(nx > 1) xStep = (xEnd - xStart)/(nx - 1);
		if(ny > 1) yStep = (yEnd - yStart)/(ny - 1);

		long ix1, iy1;
		if(x < xStart) ix1 = 0; 
		else if(x <= xEnd) ix1 = (int)((x - xStart)/xStep + 0.001*xStep); 
		else ix1 = nx - 1;

		if(y < yStart) iy1 = 0;
		else if(y <= yEnd) iy1 = (int)((y - yStart)/yStep + 0.001*yStep);
		else iy1 = ny - 1;

		long ix2 = ix1, iy2 = iy1;
		int nx_mi_1 = nx - 1;
		if(ix1 < nx_mi_1) ix2 = ix1 + 1;
		else if(ix1 >= nx_mi_1) { ix1 = nx - 2; ix2 = ix1 + 1;}

		int ny_mi_1 = ny - 1;
		if(iy1 < ny_mi_1) iy2 = iy1 + 1;
		else if(iy1 >= ny_mi_1) { iy1 = ny - 2; iy2 = iy1 + 1;}

		double xr = x - (xStart + ix1*xStep);
		double yr = y - (yStart + iy1*yStep);

		//long ofst11 = nx*iy1 + ix1;
		//long ofst12 = nx*iy2 + ix1;
		//long ofst21 = nx*iy1 + ix2;
		//long ofst22 = nx*iy2 + ix2;

		long long ofst11 = nx*iy1 + ix1;
		long long ofst12 = nx*iy2 + ix1;
		long long ofst21 = nx*iy1 + ix2;
		long long ofst22 = nx*iy2 + ix2;

		if(nx == 1)
		{
			if(ny == 1) return *(pF + ofst11);
			else return (*(pF + ofst11)) + ((*(pF + ofst12)) - (*(pF + ofst11)))*yr/yStep;
		}
		else
		{
			if(ny == 1) return (*(pF + ofst11)) + ((*(pF + ofst21)) - (*(pF + ofst11)))*xr/xStep;
			else 
			{
				double xt = xr/xStep, yt = yr/yStep;
				return (*(pF + ofst11))*(1 - xt)*(1 - yt) + (*(pF + ofst21))*xt*(1 - yt) + (*(pF + ofst12))*(1 - xt)*yt + (*(pF + ofst22))*xt*yt;
			}
		}
	}

	template <class T> static T tabTangOrtsToSurf2D(TVector3d& vHorRes, TVector3d& vVertRes, int ix, int iz, int nx, int nz, double xStep, double zStep, T* pF)
	{//calculates  "horizontal" and "vertical" tangential vector and inner normal (assuming y - longitudinal coordinate)
	 //returns function value, same as tabFunc2D
		if((pF == 0) || (nx <= 1) || (nz <= 1)) return (T)0.;

		double hx = xStep, hz = zStep;
		//long ofst0 = iz*nx + ix;
		//long ofstX1, ofstX2, ofstZ1, ofstZ2;
		long long ofst0 = ((long long)iz)*((long long)nx) + ix;
		long long ofstX1, ofstX2, ofstZ1, ofstZ2;
		if(ix == 0)
		{
			ofstX1 = ofst0; ofstX2 = ofst0 + 1;
		}
		else if(ix == (nx - 1))
		{
			ofstX1 = ofst0 - 1; ofstX2 = ofst0;
		}
		else
		{
			ofstX1 = ofst0 - 1; ofstX2 = ofst0 + 1; hx *= 2;
		}
		if(iz == 0)
		{
			ofstZ1 = ofst0; ofstZ2 = ofst0 + nx;
		}
		else if(iz == (nz - 1))
		{
			ofstZ1 = ofst0 - nx; ofstZ2 = ofst0;
		}
		else
		{
			ofstZ1 = ofst0 - nx; ofstZ2 = ofst0 + nx; hz *= 2;
		}

		vHorRes.x = hx; vHorRes.y = (double)(*(pF + ofstX2) - *(pF + ofstX1)); vHorRes.z = 0;
		vVertRes.x = 0; vVertRes.y = (double)(*(pF + ofstZ2) - *(pF + ofstZ1)); vVertRes.z = hz;
		vHorRes.Normalize(); vVertRes.Normalize();
		return *(pF + ofst0);
	}

	template <class T> static T differentiate(T* f, double h, int PoIndx, int AmOfPo)
	{//differentiate function of real argument
		if(AmOfPo==5)
		{
			if(PoIndx==2) return 0.08333333333333*(f[0] - 8.*f[1] + 8.*f[3] - f[4])/h;
			else if(PoIndx==1) return 0.08333333333333*(-3.*f[0] - 10.*f[1] + 18.*f[2] - 6.*f[3] + f[4])/h;
			else if(PoIndx==3) return 0.08333333333333*(-f[0] + 6.*f[1] - 18.*f[2] + 10.*f[3] + 3.*f[4])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
			else return 0.5*(f[2] - 4.*f[3] + 3.*f[4])/h;
			//else if(PoIndx==4) return 0.5*(f[2] - 4.*f[3] + 3.*f[4])/h;
			//else return 1.E+23;
		}
		else if(AmOfPo==4)
		{
			if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
			else if(PoIndx==2) return 0.5*(-f[1] + f[3])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
			else return 0.5*(f[1] - 4.*f[2] + 3.*f[3])/h;
			//else if(PoIndx==3) return 0.5*(f[1] - 4.*f[2] + 3.*f[3])/h;
			//else return 1.E+23;
		}
		else if(AmOfPo==3)
		{
			if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
			else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
			else return 0.5*(f[0] - 4.*f[1] + 3.*f[2])/h;
			//else if(PoIndx==2) return 0.5*(f[0] - 4.*f[1] + 3.*f[2])/h;
			//else return 1.E+23;
		}
		else return (-f[0] + f[1])/h;
		//else if(AmOfPo==2) return (-f[0] + f[1])/h;
		//else return 1.E+23;
	}

	static double radicalOnePlusSmall(double x)
	{
		if(::fabs(x) > 0.01) return sqrt(1. + x);
		return 1. + x*(0.5 + x*(-0.125 + x*(0.0625 + x*(-0.0390625 + x*(0.02734375 + x*(-0.0205078125 + x*0.01611328125))))));
	}

	static double radicalOnePlusSmallMinusOne(double x)
	{
		if(::fabs(x) > 0.01) return sqrt(1. + x) - 1.;
		return x*(0.5 + x*(-0.125 + x*(0.0625 + x*(-0.0390625 + x*(0.02734375 + x*(-0.0205078125 + x*0.01611328125))))));
	}

	static double azimAngFrXY(double x, double y)
	{//Returns azimuth angle within limits [-pi, pi) based on Cartesian coordinates (x, y)
		const double pi = 3.141592653589793;
		const double halfPi = 0.5*pi;
		if(x == 0.)
		{
			if(y < 0.) return -halfPi;
			else if(y > 0.) return halfPi;
			else return 0.;
		}

		double absAng = atan(fabs(y/x));
		if(x < 0.)
		{
			if(y < 0.) return -pi + absAng;
			else return pi - absAng;
		}
		else
		{
			if(y < 0.) return -absAng;
			else return absAng;
		}
	}
};

//-------------------------------------------------------------------------

#endif


