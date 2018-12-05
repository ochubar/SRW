/************************************************************************//**
 * File: gminterp.h
 * Description: Interpolation routines header
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMINTERP_H
#define __GMINTERP_H

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

#ifndef MATH_INTERP_STRUCT_WAS_NOT_SETUP
#define MATH_INTERP_STRUCT_WAS_NOT_SETUP 132 + 10000 //in line with SRW
#endif
#ifndef MEMORY_ALLOCATION_FAILURE
#define MEMORY_ALLOCATION_FAILURE 8 + 10000 //in line with SRW
#endif
#ifndef CAN_NOT_FIND_IND_FOR_INTERP
#define CAN_NOT_FIND_IND_FOR_INTERP 187 + 10000 //in line with SRW
#endif

//-------------------------------------------------------------------------

class CGenMathInterp
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{
	double *mSplineY2Arr, *mSplineArgTabArr, *mSplineValTabArr;
	double mArgStep, mArgStart;
	int mSplineTabNp;
	int mMethNo;

	//from former srTMathInterpol1D:
	double* AllCf;
	double** PlnCf;
	double Orig_sStart, Orig_sStep, InvStep;
	int OrigNp;

public:

	CGenMathInterp(int MethNo, double *x, double *y, int np);
	CGenMathInterp(int MethNo, double xStart, double xStep, double *y, int np);
	CGenMathInterp(double* OrigF, int InOrigNp, double InOrig_sStart, double InOrig_sStep)
	{
		if((OrigF == 0) || (InOrigNp == 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;

		mSplineY2Arr = mSplineArgTabArr = mSplineValTabArr = 0; //OC011213

        AllCf = 0;
		PlnCf = 0;
        OrigNp = InOrigNp;
        Orig_sStart = InOrig_sStart;
		Orig_sStep = InOrig_sStep;
        InvStep = 1./Orig_sStep;

		SetupPolinomCfs(OrigF);
	}

	CGenMathInterp()
	{
		mMethNo = 0;
        mSplineY2Arr = mSplineArgTabArr = mSplineValTabArr = 0;
		mArgStep = 0; mArgStart = 0;
		mSplineTabNp = 0;

		AllCf = 0;
		PlnCf = 0;
		OrigNp = 0;
	}
	~CGenMathInterp()
	{
		if(mSplineY2Arr != 0) { delete[] mSplineY2Arr; mSplineY2Arr = 0;}
		if(mSplineArgTabArr != 0) { delete[] mSplineArgTabArr; mSplineArgTabArr = 0;}
		if(mSplineValTabArr != 0) { delete[] mSplineValTabArr; mSplineValTabArr = 0;}
		
		DeallocateMemoryForCfs();
	}

	//static void CubSplinePrep(double *x, double *y, int n, double yp1, double ypn, double *y2);
	static void InterpCubicSplinePrep(double *x, double *y, int n, double *y2);
	static void InterpCubicSplinePrepU(double xStart, double xStep, double *y, int n, double *y2);
	static double InterpCubicSpline(double *xa, double *ya, double *y2a, int n, double x);
	static double Deriv1(double* f, double h, int PoIndx, int AmOfPo);
	static double Derivative(double* f, double h, int PoIndx, int AmOfPo=5);
	//static bool TryToFindRectMesh(double* arX, double* arY, int nPtTot, double* arResMeshX, double* arResMeshY, int& nMeshX, int& nMeshY, int* arPtMeshInd, double arMinMaxX[2], double arMinMaxY[2]);
	static bool TryToFindRectMesh(double* arX, double* arY, int nPtTot, double* arResMeshX, double* arResMeshY, int& nMeshX, int& nMeshY, double arMinMaxX[2], double arMinMaxY[2]);
	//static void SelectPointsForInterp2d(double x, double y, double* arX, double* arY, int nPtTot, int ord, int* arResInd, int& nResInd);
	static int SelectPointsForInterp2d(double x, double y, double* arX, double* arY, int nPtTot, int& ord, int* arResInd, int& nResInd, bool& resMeshIsRect);
	static void SelectPointsForInterp1d2d(double* arGaps, double* arPhases, int nVals, int* arResInds, int& nResInds, double arPrecPar[5]);
	static int TryToFindMeshPointForPars(double* arPar1, double* arPar2, int nVals, double* arPrecPar);

	void Interpolate(double sSt, double sStp, int Np, double* pInterpData);
	void CompDerivForOrigData(double* OrigF, double* DerF);

	void InitCubicSpline(double *x, double *y, int np);
	void InitCubicSplineU(double xStart, double xStep, double *y, int np);

	double Interp1D(double x)
	{
		if((mMethNo == 1) && (mSplineY2Arr != 0) && (mSplineArgTabArr != 0) && (mSplineValTabArr != 0) && (mSplineTabNp > 0))
		{
			return InterpCubicSpline(mSplineArgTabArr, mSplineValTabArr, mSplineY2Arr, mSplineTabNp, x);
		}
		return 0;
	}
	
	double InterpRelCubicSpline(double b, int i0, double curArgStep)
	{
		if((mSplineY2Arr == 0) || (mSplineValTabArr == 0)) return 0;

		//b = (x - xa[klo])/h; 
		double a = 1. - b; //(xa[khi] - x)/h; //Cubic spline polynomial is now evaluated.
		double ya_lo = mSplineValTabArr[i0], ya_hi = mSplineValTabArr[i0 + 1];
		double y2a_lo = mSplineY2Arr[i0], y2a_hi = mSplineY2Arr[i0 + 1];
		return a*ya_lo + b*ya_hi + ((a*a*a - a)*y2a_lo + (b*b*b - b)*y2a_hi)*(curArgStep*curArgStep)/6.0;
	}
	double InterpRelCubicSplineU(double b, int i0)
	{
		if((mSplineY2Arr == 0) || (mSplineValTabArr == 0) || (mArgStep == 0)) return 0;

		//b = (x - xa[klo])/h; 
		double a = 1. - b; //(xa[khi] - x)/h; //Cubic spline polynomial is now evaluated.
		double ya_lo = mSplineValTabArr[i0], ya_hi = mSplineValTabArr[i0 + 1];
		double y2a_lo = mSplineY2Arr[i0], y2a_hi = mSplineY2Arr[i0 + 1];
		return a*ya_lo + b*ya_hi + ((a*a*a - a)*y2a_lo + (b*b*b - b)*y2a_hi)*(mArgStep*mArgStep)/6.0;
	}

	void SetupPolinomCfs(double* OrigF)
	{
		if((OrigF == 0) || (OrigNp <= 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
		AllocateMemoryForCfs();

		double* DerF = new double[OrigNp];
		CompDerivForOrigData(OrigF, DerF);
		CalcPlnCfs(OrigF, DerF);
		if(DerF != 0) delete[] DerF;
	}
	void CalcPlnCfs(double* OrigF, double* DerF)
	{
		if((OrigF == 0) || (DerF == 0) || (OrigNp <= 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
        //int LenFieldData_m_1 = OrigNp - 1; //OC020110
        double f1 = OrigF[0], f2;
        double fpr1 = DerF[0], fpr2;
        for(int is=1; is<OrigNp; is++)
        {
            f2 = OrigF[is];
            fpr2 = DerF[is];
            CubicPlnCfs(f1, f2, fpr1, fpr2, Orig_sStep, PlnCf[is-1]);
            f1 = f2; fpr1 = fpr2;
		}
	}
    void AllocateMemoryForCfs()
    {
		if(OrigNp <= 0) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
        DeallocateMemoryForCfs();
        int LenFieldData_m_1 = OrigNp - 1;

        PlnCf = new double*[LenFieldData_m_1];
        if(PlnCf == 0) throw MEMORY_ALLOCATION_FAILURE;
        AllCf = new double[LenFieldData_m_1*4];
        if(AllCf == 0) { delete[] PlnCf; throw MEMORY_ALLOCATION_FAILURE;}
        double* tAllCf = AllCf;
        for(int i=0; i<LenFieldData_m_1; i++) { PlnCf[i] = tAllCf; tAllCf += 4;}
    }
    void DeallocateMemoryForCfs()
    {
        if(AllCf != 0) { delete[] AllCf; AllCf = 0;}
        if(PlnCf != 0) { delete[] PlnCf; PlnCf = 0;}
    }

	static int FindIndOfPointOnIrregMesh2d(double x, double y, double* arX, double* arY, int nPtTot, double relTol=1.e-09)
	{
		//if((nPtTot <= 0) || (arX == 0) || (arY == 0)) return -1;
		if((nPtTot <= 0) || (arX == 0) || (arY == 0)) throw CAN_NOT_FIND_IND_FOR_INTERP;

		if(nPtTot == 1)
		{
			double dXi = x - (*arX);
			double curTolX = (*arX)*relTol;
			if((-curTolX <= dXi) && (dXi <= curTolX))
			{
				double dYi = y - (*arY);
				double curTolY = (*arY)*relTol;
				if((-curTolY <= dYi) && (dYi <= curTolY))
				{
					return 0;
				}
			}
		}

		//double curTolX = (arX[1] - arX[0])*relTol;
		//double curTolY = (arY[1] - arY[0])*relTol;
		int nPtTot_mi_1 = nPtTot - 1; //OC31102017
		double curTolX = (arX[nPtTot_mi_1] - arX[0])*relTol;
		double curTolY = (arY[nPtTot_mi_1] - arY[0])*relTol;
		for(int i=0; i<nPtTot; i++)
		{
			double dXi = x - arX[i];
			if((-curTolX <= dXi) && (dXi <= curTolX))
			{
				double dYi = y - arY[i];
				if((-curTolY <= dYi) && (dYi <= curTolY)) return i;
			}

			//OC31102017 (commented-out)
			//if(i > 1)
			//{
			//	curTolX = (arX[i] - arX[i - 1])*relTol;
			//	curTolY = (arY[i] - arY[i - 1])*relTol;
			//}
		}
		throw CAN_NOT_FIND_IND_FOR_INTERP;
		return -1;
	}

	static double Interp3dBilin(double* inP, double* inArrArgBounds, double* inArrFunc)
	{
	// The approach taken from Numerical Recipes (p. 104), extended to 3d
	// assumes:
	// inP[] = {x,y,z}; //point were the function should be computed
	// inArrArgBounds[] = {x0,x1,y0,y1,z0,z1}; //argument values at the corners of the cube
	// inArrFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
		double xt = 0, yt = 0, zt = 0;
		double &x0 = inArrArgBounds[0], &x1 = inArrArgBounds[1], &y0 = inArrArgBounds[2], &y1 = inArrArgBounds[3], &z0 = inArrArgBounds[4], &z1 = inArrArgBounds[5];
		if(x1 != x0) xt = (inP[0] - x0)/(x1 - x0);
		if(y1 != y0) yt = (inP[1] - y0)/(y1 - y0);
		if(z1 != z0) zt = (inP[2] - z0)/(z1 - z0);
		double one_mi_xt = 1 - xt, one_mi_yt = 1 - yt, one_mi_zt = 1 - zt;
		return inArrFunc[0]*one_mi_xt*one_mi_yt*one_mi_zt
			+ inArrFunc[1]*xt*one_mi_yt*one_mi_zt
			+ inArrFunc[2]*one_mi_xt*yt*one_mi_zt
			+ inArrFunc[3]*one_mi_xt*one_mi_yt*zt
			+ inArrFunc[4]*xt*yt*one_mi_zt
			+ inArrFunc[5]*xt*one_mi_yt*zt
			+ inArrFunc[6]*one_mi_xt*yt*zt
			+ inArrFunc[7]*xt*yt*zt;
	}

	static double Interp3dBilinRel(double xt, double yt, double zt, double* inArFunc)
	{// The approach taken from Numerical Recipes (p. 104), extended to 3d
	 // assumes:
	 // inArFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
		double one_mi_xt = 1.- xt, one_mi_yt = 1.- yt, one_mi_zt = 1.- zt;
		return inArFunc[0]*one_mi_xt*one_mi_yt*one_mi_zt
			+ inArFunc[1]*xt*one_mi_yt*one_mi_zt
			+ inArFunc[2]*one_mi_xt*yt*one_mi_zt
			+ inArFunc[3]*one_mi_xt*one_mi_yt*zt
			+ inArFunc[4]*xt*yt*one_mi_zt
			+ inArFunc[5]*xt*one_mi_yt*zt
			+ inArFunc[6]*one_mi_xt*yt*zt
			+ inArFunc[7]*xt*yt*zt;
	}

	static double Interp3dQuadRel(double xt, double yt, double zt, double* arF)
	{//assumes:
	 //double arF[] = {f(x0,y0,zm1),f(x0,ym1,z0),f(xm1,y0,z0),f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x1,y1,z0),f(x0,y0,z1),f(x1,y0,z1),f(x0,y1,z1)};
		double ax2 = 0.5*(-2.*arF[3] + arF[4] + arF[2]); //1/2 (-2 f000 + f100 + fm100)
		double ay2 = 0.5*(-2.*arF[3] + arF[5] + arF[1]); //1/2 (-2 f000 + f010 + f0m10)
		double az2 = 0.5*(-2.*arF[3] + arF[7] + arF[0]); //1/2 (-2 f000 + f001 + f00m1)
		double axy = arF[3] - arF[5] - arF[4] + arF[6]; //f000 - f010 - f100 + f110
		double axz = arF[3] - arF[7] - arF[4] + arF[8]; //f000 - f001 - f100 + f101
		double ayz = arF[3] - arF[7] - arF[5] + arF[9]; //f000 - f001 - f010 + f011
		double ax = 0.5*(arF[4] - arF[2]); //(f100 - fm100)/2
		double ay = 0.5*(arF[5] - arF[1]); //(f010 - f0m10)/2
		double az = 0.5*(arF[7] - arF[0]); //(f001 - f00m1)/2
		return arF[3] + xt*(ax + xt*ax2 + yt*axy + zt*axz) + yt*(ay + yt*ay2 + zt*ayz) + zt*(az + zt*az2);
	}

	static double Interp3dCubicRel(double xt, double yt, double zt, double* arF)
	{//assumes:
	 //double arF[] = {f00m1,f0m10,fm100,f000,f100,f200,f010,f110,f210,f020,f120,f001,f101,f201,f011,f111,f021,f002,f102,f012};
		const double i6 = 1./6.;
		double ax3 = i6*(3.*(arF[3] - arF[4]) + arF[5] - arF[2]); //1/6 (3 f000 - 3 f100 + f200 - fm100)
		double ay3 = i6*(3.*(arF[3] - arF[6]) + arF[9] - arF[1]); //1/6 (3 f000 - 3 f010 + f020 - f0m10)
		double az3 = i6*(3.*(arF[3] - arF[11]) + arF[17] - arF[0]); //1/6 (3 f000 - 3 f001 + f002 - f00m1)
		double ax2y = 0.5*(-arF[3] + arF[6] + 2.*(arF[4] - arF[7]) - arF[5] + arF[8]); //1/2 (-f000 + f010 + 2 f100 - 2 f110 - f200 + f210)
		double ax2z = 0.5*(-arF[3] + arF[11] + 2.*(arF[4] - arF[12]) - arF[5] + arF[13]); //1/2 (-f000 + f001 + 2 f100 - 2 f101 - f200 + f201)
		double axy2 = 0.5*(-arF[3] + 2.*(arF[6] - arF[7]) - arF[9] + arF[4] + arF[10]); //1/2 (-f000 + 2 f010 - f020 + f100 - 2 f110 + f120)
		double axz2 = 0.5*(-arF[3] + 2.*arF[11] - arF[17] + arF[4] - 2.*arF[12] + arF[18]); //1/2 (-f000 + 2 f001 - f002 + f100 - 2 f101 + f102)
		double ay2z = 0.5*(-arF[3] + arF[11] + 2.*arF[6] - 2.*arF[14] - arF[9] + arF[16]); //1/2 (-f000 + f001 + 2 f010 - 2 f011 - f020 + f021)
		double ayz2 = 0.5*(-arF[3] + 2.*arF[11] - arF[17] + arF[6] - 2.*arF[14] + arF[19]); //1/2 (-f000 + 2 f001 - f002 + f010 - 2 f011 + f012)
		double axyz = -arF[3] + arF[11] + arF[6] - arF[14] + arF[4] - arF[12] - arF[7] + arF[15]; //-f000 + f001 + f010 - f011 + f100 - f101 - f110 + f111
		double ax2 = 0.5*(-2.*arF[3] + arF[4] + arF[2]); //1/2 (-2 f000 + f100 + fm100)
		double ay2 = 0.5*(-2.*arF[3] + arF[6] + arF[1]); //1/2 (-2 f000 + f010 + f0m10)
		double az2 = 0.5*(-2.*arF[3] + arF[11] + arF[0]); //1/2 (-2 f000 + f001 + f00m1)
		double axy = 0.5*(4.*arF[3] - 5.*arF[6] + arF[9] - 5.*arF[4] + 6.*arF[7] - arF[10] + arF[5] - arF[8]); //1/2 (4 f000 - 5 f010 + f020 - 5 f100 + 6 f110 - f120 + f200 - f210)
		double axz = 0.5*(4.*arF[3] - 5.*arF[11] + arF[17] - 5.*arF[4] + 6.*arF[12] - arF[18] + arF[5] - arF[13]); //1/2 (4 f000 - 5 f001 + f002 - 5 f100 + 6 f101 - f102 + f200 - f201)
		double ayz = 0.5*(4.*arF[3] - 5.*arF[11] + arF[17] - 5.*arF[6] + 6.*arF[14] - arF[19] + arF[9] - arF[16]); //1/2 (4 f000 - 5 f001 + f002 - 5 f010 + 6 f011 - f012 + f020 - f021)
		double ax = i6*(-3.*arF[3] + 6.*arF[4] - arF[5] - 2.*arF[2]); //1/6 (-3 f000 + 6 f100 - f200 - 2 fm100)
		double ay = i6*(-3.*arF[3] + 6.*arF[6] - arF[9] - 2.*arF[1]); //1/6 (-3 f000 + 6 f010 - f020 - 2 f0m10)
		double az = i6*(-3.*arF[3] + 6.*arF[11] - arF[17] - 2.*arF[0]); //1/6 (-3 f000 + 6 f001 - f002 - 2 f00m1)
		return arF[3] + xt*(ax + xt*(ax2 + xt*ax3 + yt*ax2y + zt*ax2z) + yt*(axy + yt*axy2 + zt*axyz) + zt*(axz + zt*axz2)) +
						yt*(ay + yt*(ay2 + yt*ay3 + xt*axy2 + zt*ay2z) + zt*(ayz + zt*ayz2)) + zt*(az + zt*(az2 + zt*az3));
	}

	static double Interp3dBiCubic32pRel(double xt, double yt, double zt, double* arF)
	{
		double *p = arF;
		double f00m1,f10m1,f01m1,f11m1,f0m10,f1m10,fm100,f000,f100,f200,fm110,f010,f110,f210,f020,f120,f0m11,f1m11,fm101,f001,f101,f201,fm111,f011,f111,f211,f021,f121,f002,f102,f012,f112;
		f00m1=*(p++); f10m1=*(p++); f01m1=*(p++); f11m1=*(p++);
		f0m10=*(p++); f1m10=*(p++); fm100=*(p++); f000=*(p++); f100=*(p++); f200=*(p++); fm110=*(p++); f010=*(p++); f110=*(p++); f210=*(p++); f020=*(p++); f120=*(p++); 
		f0m11=*(p++); f1m11=*(p++); fm101=*(p++); f001=*(p++); f101=*(p++); f201=*(p++); fm111=*(p++); f011=*(p++); f111=*(p++); f211=*(p++); f021=*(p++); f121=*(p++); 
		f002=*(p++); f102=*(p++); f012=*(p++); f112=*(p++);

		const double i6 = 1./6.;
		double b5 = 3.*(f000 - f001 - f010 + f011 - f100 + f101 + f110 - f111);
		double a311 = i6*(b5 + f200 - f201 - f210 + f211 - fm100 + fm101 + fm110 - fm111);
		double a131 = i6*(b5 + f020 - f021 - f0m10 + f0m11 - f120 + f121 + f1m10 - f1m11);
		double a113 = i6*(b5 + f002 - f00m1 - f102 - f012 + f01m1 + f10m1 + f112 - f11m1);
		double b4xy = 3.*(f010 - f000 + f100 - f110);
		double a310 = i6*(b4xy - f200 + f210 + fm100 - fm110);
		double a130 = i6*(b4xy - f020 + f0m10 + f120 - f1m10);
		double b4xz = 3.*(f001 - f000 + f100 - f101);
		double a301 = i6*(b4xz - f200 + f201 + fm100 - fm101);
		double a103 = i6*(b4xz - f002 + f00m1 + f102 - f10m1);
		double b4yz = 3.*(f001 - f000 + f010 - f011);
		double a031 = i6*(b4yz - f020 + f021 + f0m10 - f0m11);
		double a013 = i6*(b4yz - f002 + f00m1 + f012 - f01m1);
		double a211 = 0.5*(2.*(f001 - f000 + f010 - f011) + f100 - f101 - f110 + f111 + fm100 - fm101 - fm110 + fm111);
		double a121 = 0.5*(2.*(f001 - f000 + f100 - f101) + f010 - f011 + f0m10 - f0m11 - f110 + f111 - f1m10 + f1m11);
		double a112 = 0.5*(2.*(f010 - f000 + f100 - f110) + f001 + f00m1 - f011 - f01m1 - f101 - f10m1 + f111 + f11m1);
		double a300 = i6*(3.*(f000 - f100) + f200 - fm100);
		double a030 = i6*(3.*(f000 - f010) + f020 - f0m10);
		double a003 = i6*(3.*(f000 - f001) + f002 - f00m1);
		double a210 = 0.5*(2.*(f000 - f010) - f100 + f110 - fm100 + fm110); 
		double a201 = 0.5*(2.*(f000 - f001) - f100 + f101 - fm100 + fm101);
		double a120 = 0.5*(2.*(f000 - f100) - f010 - f0m10 + f110 + f1m10);
		double a021 = 0.5*(2.*(f000 - f001) - f010 + f011 - f0m10 + f0m11);
		double a102 = 0.5*(2.*(f000 - f100) - f001 - f00m1 + f101 + f10m1);
		double a012 = 0.5*(2.*(f000 - f010) - f001 - f00m1 + f011 + f01m1);
		double a111 = f111 + i6*(3.*(f000 - f011 - f101 - f110) + 2.*(f01m1 - f00m1 - f0m10 + f0m11 + f10m1 - f11m1 + f1m10 - f1m11 - fm100 + fm101 + fm110 - fm111) - f002 + f012 - f020 + f021 + f102 - f112 + f120 - f121 - f200 + f201 + f210 - f211);
		double a200 = 0.5*(f100 + fm100) - f000;
		double a020 = 0.5*(f010 + f0m10) - f000;
		double a002 = 0.5*(f001 + f00m1) - f000;
		double a110 = f110 + i6*(-3.*(f010 + f100) + 2.*(f0m10 - f1m10 + fm100 - fm110) + f020 - f120 + f200 - f210);
		double a101 = f101 + i6*(-3.*(f001 + f100) + 2.*(f00m1 - f10m1 + fm100 - fm101) + f002 - f102 + f200 - f201);
		double a011 = f011 + i6*(-3.*(f001 + f010) + 2.*(f00m1 - f01m1 + f0m10 - f0m11) + f002 - f012 + f020 - f021);
		double a100 = f100 + i6*(-3.*f000 - 2.*fm100 - f200);
		double a010 = f010 + i6*(-3.*f000 - 2.*f0m10 - f020);
		double a001 = f001 + i6*(-3.*f000 - 2.*f00m1 - f002);

		return f000 + xt*(a100 + xt*(a200 + xt*(a300 + yt*(a310 + a311*zt) + a301*zt) + yt*(a210 + a211*zt) + a201*zt) 
							   + yt*(a110 + yt*(a120 + yt*(a130 + a131*zt) + a121*zt) + zt*(a111 + zt*(a112 + a113*zt))) 
							   + zt*(a101 + zt*(a102 + a103*zt))) 
					+ yt*(a010 + yt*(a020 + yt*(a030 + a031*zt) + a021*zt) + zt*(a011 + zt*(a012 + a013*zt)))
					+ zt*(a001 + zt*(a002 + a003*zt));
	}

	static double Interp2dBiCubic12pRel(double xt, double yt, double* arF)
	{
		double *p = arF;
		double f0m1,f1m1,fm10,f00,f10,f20,fm11,f01,f11,f21,f02,f12;
		f0m1=*(p++); f1m1=*(p++); fm10=*(p++); f00=*(p++); f10=*(p++); f20=*(p++); fm11=*(p++); f01=*(p++); f11=*(p++); f21=*(p++); f02=*(p++); f12=*(p++); 

		const double i6 = 1./6.;
		double b4xy = 3.*(f01 - f00 + f10 - f11);
		double a31 = i6*(b4xy - f20 + f21 + fm10 - fm11);
		double a13 = i6*(b4xy - f02 + f0m1 + f12 - f1m1);
		double a30 = i6*(3.*(f00 - f10) + f20 - fm10);
		double a03 = i6*(3.*(f00 - f01) + f02 - f0m1);
		double a21 = 0.5*(2.*(f00 - f01) - f10 + f11 - fm10 + fm11); 
		double a12 = 0.5*(2.*(f00 - f10) - f01 - f0m1 + f11 + f1m1);
		double a20 = 0.5*(f10 + fm10) - f00;
		double a02 = 0.5*(f01 + f0m1) - f00;
		double a11 = f11 + i6*(-3.*(f01 + f10) + 2.*(f0m1 - f1m1 + fm10 - fm11) + f02 - f12 + f20 - f21);
		double a10 = f10 + i6*(-3.*f00 - 2.*fm10 - f20);
		double a01 = f01 + i6*(-3.*f00 - 2.*f0m1 - f02);

		return f00 + xt*(a10 + xt*(a20 + xt*(a30 + yt*a31) + yt*a21) + yt*(a11 + yt*(a12 + yt*a13))) + yt*(a01 + yt*(a02 + yt*a03));
	}

	static double Interp2dBiCubic12pRecVar(double x, double y, double* arXY, double* arF)
	{//bi-cubic interpolation on rectangular, but variable step-size mesh (12 points), for relative arguments, "central" point (that correspopnds to f00) is x = 0, y = 0
		double *p = arF;
		double f0m1,f1m1,fm10,f00,f10,f20,fm11,f01,f11,f21,f02,f12;
		f0m1=*(p++); f1m1=*(p++); fm10=*(p++); f00=*(p++); f10=*(p++); f20=*(p++); fm11=*(p++); f01=*(p++); f11=*(p++); f21=*(p++); f02=*(p++); f12=*(p++); 

		double xm1 = *(arXY++), ym1 = *(arXY++);
		//double x0 = 0, y0 = 0;
		double x1 = *(arXY++), y1 = *(arXY++);
		double x2 = *(arXY++), y2 = *(arXY++);

		double x1_mi_xm1 = x1 - xm1, x1_mi_x2 = x1 - x2, x2_mi_xm1 = x2 - xm1;
		double y1_mi_ym1 = y1 - ym1, y2_mi_ym1 = y2 - ym1, y1_mi_y2 = y1 - y2;

		double fm10_d_x1_mi_xm1_x2_mi_xm1_xm1 = fm10/(x1_mi_xm1*x2_mi_xm1*xm1);
		double x1_mi_x2_x2_x2_mi_xm1 = x1_mi_x2*x2*x2_mi_xm1;
		double x1_x1_mi_x2_x1_mi_xm1 = x1*x1_mi_x2*x1_mi_xm1;
		double f10_d_x1_x1_mi_x2_x1_mi_xm1 = f10/x1_x1_mi_x2_x1_mi_xm1;
		double f20_d_x1_mi_x2_x2_x2_mi_xm1 = f20/x1_mi_x2_x2_x2_mi_xm1;
		double x1_mi_xm1_xm1_x2_mi_xm1 = x1_mi_xm1*xm1*x2_mi_xm1;
		double x1_mi_xm1_xm1_x2_mi_xm1_y1 = x1_mi_xm1_xm1_x2_mi_xm1*y1;
		double y1_mi_ym1_y2_mi_ym1_ym1 = y1_mi_ym1*y2_mi_ym1*ym1;
		double f0m1_d_y1_mi_ym1_y2_mi_ym1_ym1 = f0m1/y1_mi_ym1_y2_mi_ym1_ym1;
		double x1_p_x2 = x1 + x2, x1_p_xm1 = x1 + xm1, x2_p_xm1 = x2 + xm1;
		double x1_x2 = x1*x2, x1_xm1 = x1*xm1, x2_xm1 = x2*xm1;
		double x1_x2_xm1 = x1_x2*xm1;
		double y1_mi_ym1_ym1_y2_mi_ym1 = y1_mi_ym1*ym1*y2_mi_ym1;
		double x1_y1_mi_ym1_ym1_y2_mi_ym1 = x1*y1_mi_ym1_ym1_y2_mi_ym1;
		double f0m1_mi_f1m1_d_x1_y1_mi_ym1_ym1_y2_mi_ym1 = (f0m1 - f1m1)/x1_y1_mi_ym1_ym1_y2_mi_ym1;
		double y1_mi_y2_y2_y2_mi_ym1 = y1_mi_y2*y2*y2_mi_ym1;
		double f02_d_y1_mi_y2_y2_y2_mi_ym1 = f02/y1_mi_y2_y2_y2_mi_ym1, x1_y1_mi_y2_y2_y2_mi_ym1 = x1*y1_mi_y2_y2_y2_mi_ym1;
		double y1_mi_y2_y1_mi_ym1 = y1_mi_y2*y1_mi_ym1;
		double y1_mi_y2_y1_mi_ym1_x1 = y1_mi_y2_y1_mi_ym1*x1;
		double y1_y1_mi_y2_y1_mi_ym1 = y1*y1_mi_y2_y1_mi_ym1;
		double f01_d_y1_y1_mi_y2_y1_mi_ym1 = f01/y1_y1_mi_y2_y1_mi_ym1;
		double f20_mi_f21_d_x1_mi_x2_x2_x2_mi_xm1_y1 = (f20 - f21)/(x1_mi_x2_x2_x2_mi_xm1*y1);
		double x1_x1_mi_x2_x1_mi_xm1_y1 = x1_x1_mi_x2_x1_mi_xm1*y1;
		double f11_mi_f10_d_x1_x1_mi_x2_x1_mi_xm1_y1 = (f11 - f10)/x1_x1_mi_x2_x1_mi_xm1_y1;
		double x1_mi_x2_x2_mi_xm1_x2 = x1_mi_x2*x2_mi_xm1*x2;
		double fm10_mi_fm11_d_x1_mi_xm1_xm1_x2_mi_xm1_y1 = (fm10 - fm11)/x1_mi_xm1_xm1_x2_mi_xm1_y1;
		double y2_y1_mi_y2_y2_mi_ym1 = y2*y1_mi_y2*y2_mi_ym1;
		double y1_p_y2 = y1 + y2, y1_p_ym1 = y1 + ym1, y2_p_ym1 = y2 + ym1;
		double y1_y2 = y1*y2, y1_ym1 = y1*ym1, y2_ym1 = y2*ym1;
		double y1_y2_ym1 = y1_y2*ym1;
		double f00_mi_f10_d_x1_y1_y2_ym1 = (f00 - f10)/(x1*y1_y2_ym1);
		double f00_mi_f01_d_x1_x2_xm1_y1 = (f00 - f01)/(x1_x2_xm1*y1);
		double f02_mi_f12_d_x1_y1_mi_y2_y2_y2_mi_ym1 = (f02 - f12)/x1_y1_mi_y2_y2_y2_mi_ym1;
		double f11_mi_f01_d_x1_y1_y1_mi_y2_y1_mi_ym1 = (f11 - f01)/(x1*y1_y1_mi_y2_y1_mi_ym1);
		double x2_xm1_d_x1_x1_mi_x2_x1_mi_xm1_y1 = x2_xm1/x1_x1_mi_x2_x1_mi_xm1_y1;

		double a10 = x1_x2*fm10_d_x1_mi_xm1_x2_mi_xm1_xm1 - x1_xm1*f20_d_x1_mi_x2_x2_x2_mi_xm1 + x2_xm1*f10_d_x1_x1_mi_x2_x1_mi_xm1 - f00*(1/x1 + 1/x2 + 1/xm1);
		double a20 = -x1_p_x2*fm10_d_x1_mi_xm1_x2_mi_xm1_xm1 + x1_p_xm1*f20_d_x1_mi_x2_x2_x2_mi_xm1 - x2_p_xm1*f10_d_x1_x1_mi_x2_x1_mi_xm1 + f00*(1/x1_x2 + 1/x1_xm1 + 1/x2_xm1);
		double a30 = fm10_d_x1_mi_xm1_x2_mi_xm1_xm1 + f10_d_x1_x1_mi_x2_x1_mi_xm1 - f20_d_x1_mi_x2_x2_x2_mi_xm1 - f00/x1_x2_xm1;
		double a01 = y1_y2*f0m1_d_y1_mi_ym1_y2_mi_ym1_ym1 - y1_ym1*f02_d_y1_mi_y2_y2_y2_mi_ym1 + y2_ym1*f01_d_y1_y1_mi_y2_y1_mi_ym1 - f00*(1/y1 + 1/y2 + 1/ym1);
		double a11 = -y1_y2*f0m1_mi_f1m1_d_x1_y1_mi_ym1_ym1_y2_mi_ym1 + y1_ym1*f02_mi_f12_d_x1_y1_mi_y2_y2_y2_mi_ym1 - x1_x2*fm10_mi_fm11_d_x1_mi_xm1_xm1_x2_mi_xm1_y1 + x1_xm1*f20_mi_f21_d_x1_mi_x2_x2_x2_mi_xm1_y1 - f10*(x2_xm1_d_x1_x1_mi_x2_x1_mi_xm1_y1 + 1/(x1*y2) + 1/(x1*ym1)) - f01*(1/(x2*y1) + 1/(xm1*y1) + y2_ym1/(y1*y1_mi_y2_y1_mi_ym1_x1)) + f11*(x2_xm1_d_x1_x1_mi_x2_x1_mi_xm1_y1 + (ym1 - y1_mi_y2)/y1_mi_y2_y1_mi_ym1_x1) + f00*(1/(x1*y1) + 1/(x2*y1) + 1/(xm1*y1) + 1/(x1*y2) + 1/(x1*ym1));
		double a21 = x1_p_x2*fm10_mi_fm11_d_x1_mi_xm1_xm1_x2_mi_xm1_y1 - x1_p_xm1*f20_mi_f21_d_x1_mi_x2_x2_x2_mi_xm1_y1 - x2_p_xm1*f11_mi_f10_d_x1_x1_mi_x2_x1_mi_xm1_y1 - (x1_p_x2 + xm1)*f00_mi_f01_d_x1_x2_xm1_y1;
		double a31 = -fm10_mi_fm11_d_x1_mi_xm1_xm1_x2_mi_xm1_y1 + f00_mi_f01_d_x1_x2_xm1_y1 + f11_mi_f10_d_x1_x1_mi_x2_x1_mi_xm1_y1 + f20_mi_f21_d_x1_mi_x2_x2_x2_mi_xm1_y1;
		double a02 = -y1_p_y2*f0m1_d_y1_mi_ym1_y2_mi_ym1_ym1 + y1_p_ym1*f02_d_y1_mi_y2_y2_y2_mi_ym1 + f00*(1/y1_y2 + 1/y1_ym1 + 1/y2_ym1) - y2_p_ym1*f01_d_y1_y1_mi_y2_y1_mi_ym1; 
		double a12 = y1_p_y2*f0m1_mi_f1m1_d_x1_y1_mi_ym1_ym1_y2_mi_ym1 - y1_p_ym1*f02_mi_f12_d_x1_y1_mi_y2_y2_y2_mi_ym1 - y2_p_ym1*f11_mi_f01_d_x1_y1_y1_mi_y2_y1_mi_ym1 - (y1_p_y2 + ym1)*f00_mi_f10_d_x1_y1_y2_ym1;
		double a03 = f0m1_d_y1_mi_ym1_y2_mi_ym1_ym1 + f01_d_y1_y1_mi_y2_y1_mi_ym1 - f02_d_y1_mi_y2_y2_y2_mi_ym1 - f00/y1_y2_ym1;
		double a13 = f11_mi_f01_d_x1_y1_y1_mi_y2_y1_mi_ym1 - f0m1_mi_f1m1_d_x1_y1_mi_ym1_ym1_y2_mi_ym1 + f00_mi_f10_d_x1_y1_y2_ym1 + f02_mi_f12_d_x1_y1_mi_y2_y2_y2_mi_ym1;

		return f00 + x*(a10 + x*(a20 + a21*y + x*(a30 + a31*y)) + y*(a11 + y*(a12 + a13*y))) + y*(a01 + y*(a02 + a03*y));
	}

	static double Interp2dBiLinRec(double xt, double yt, double* arF)
	{//bilinear interpolation on rectangular mesh, for normalized arguments (0 <= xt <= 1, 0 <= yt <= 1)
		//double *p = arF;
		double f00 = *(arF++);
		double f10 = *(arF++);
		double f01 = *(arF++);
		double f11 = *arF;
		double a01 = f01 - f00;
		double a10 = f10 - f00;
		double a11 = f00 - f01 - f10 + f11;
		return xt*(a10 + a11*yt) + a01*yt + f00;
	}

	static double Interp2dBiLinVar(double x, double y, double* arXY, double* arF)
	{//bilinear interpolation on irregular mesh, for relative arguments, first point is x = 0, y = 0
	 //arXY is flat array of coordinates of 3 other points {x10, y10, x01, y01, x11, y11}
		double f00 = *(arF++);
		double f10 = *(arF++);
		double f01 = *(arF++);
		double f11 = *(arF++);
		double x10 = *(arXY++);
		double y10 = *(arXY++);
		double x01 = *(arXY++);
		double y01 = *(arXY++);
		double x11 = *(arXY++);
		double y11 = *arXY;
		double D = 1./(x10*x11*y01*(y10 - y11) + x01*(x10*(y01 - y10)*y11 + x11*y10*(y11 - y01)));
		double a01 = D*(f11*x01*x10*(y01 - y10) + x11*(f01*x10*(y10 - y11) + f10*x01*(y11 - y01)) + f00*(x10*x11*(y11 - y10) + x01*((x11 - x10)*y01 + x10*y10 - x11*y11)));
		double a10 = D*((f00 - f11)*(x01 - x10)*y01*y10 + ((f10 - f00)*(x01 - x11)*y01 + (f00 - f01)*(x10 - x11)*y10)*y11);
		double a11 = D*(f10*(x11*y01 - x01*y11) + f11*(x01*y10 - x10*y01) + f01*(x10*y11 - x11*y10) + f00*(x11*(y10 - y01) + x01*(y11 - y10) + x10*(y01 - y11)));
		return x*(a10 + a11*y) + a01*y + f00;
	}

	static double Interp2dBiQuad5Rec(double xt, double yt, double* arF)
	{//bi-quadratic (5-point) interpolation on rectangular regular mesh, for normalized arguments (-1 <= xt <= 1, -1 <= yt <= 1)
		//double *p = arF;
		double f0m1 = *(arF++);
		double fm10 = *(arF++);
		double f00 = *(arF++);
		double f10 = *(arF++);
		double f01 = *arF;
		double a10 = 0.5*(f10 - fm10);
		double a01 = 0.5*(f01 - f0m1);
		double a20 = 0.5*(f10 + fm10) - f00;
		double a02 = 0.5*(f01 + f0m1) - f00;
		return xt*(xt*a20 + a10) + yt*(yt*a02 + a01) + f00;
	}

	static double Interp2dBiQuad5RecVar(double x, double y, double* arXY, double* arF)
	{//bi-quadratic interpolation on rectangular, but variable step-size mesh, for relative arguments, central (3rd) point is x = 0, y = 0
		//double *p = arF;
		double f0m1 = *(arF++);
		double fm10 = *(arF++);
		double f00 = *(arF++);
		double f10 = *(arF++);
		double f01 = *arF;
		double xm1 = *(arXY++);
		double ym1 = *(arXY++);
		double x1 = *(arXY++);
		double y1 = *arXY;
		double DX = 1./((x1 - xm1)*x1*xm1);
		double DY = 1./((y1 - ym1)*y1*ym1);
		double f00_mi_fm10_x1 = (f00 - fm10)*x1;
		double f00_mi_f0m1_y1 = (f00 - f0m1)*y1;
		double f00_mi_f10_xm1 = (f00 - f10)*xm1;
		double f00_mi_f01_ym1 = (f00 - f01)*ym1;
		double a10 = DX*(f00_mi_f10_xm1*xm1 - f00_mi_fm10_x1*x1);
		double a01 = DY*(f00_mi_f01_ym1*ym1 - f00_mi_f0m1_y1*y1);
		double a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1);
		double a02 = DY*(f00_mi_f0m1_y1 - f00_mi_f01_ym1);
		return x*(x*a20 + a10) + y*(y*a02 + a01) + f00;
	}

	static double Interp2dBiQuad5Var(double x, double y, double* arXY, double* arF)
	{//bi-quadratic interpolation on irregular mesh, for relative arguments, central (3rd) point is x = 0, y = 0
		//double *p = arF;
		double f0m1 = *(arF++);
		double fm10 = *(arF++);
		double f00 = *(arF++);
		double f10 = *(arF++);
		double f01 = *arF;
		double x0m1 = *(arXY++);
		double y0m1 = *(arXY++);
		double xm10 = *(arXY++);
		double ym10 = *(arXY++);
		double x10 = *(arXY++);
		double y10 = *(arXY++);
		double x01 = *(arXY++);
		double y01 = *arXY;

		double x0m1e2 = x0m1*x0m1;
		double y0m1e2 = y0m1*y0m1;
		double xm10e2 = xm10*xm10;
		double ym10e2 = ym10*ym10;
		double x10e2 = x10*x10;
		double y10e2 = y10*y10;
		double x01e2 = x01*x01;
		double y01e2 = y01*y01;

		double y01_y10_mi_y01 = y01*(y10 - y01);
		double y01_y01_mi_y0m1 = y01*(y01 - y0m1);
		double y0m1_y0m1_mi_y01 = y0m1*(y0m1 - y01);
		double y01_y01_mi_y10 = y01*(y01 - y10);
		double x01_x01_mi_xm10 = x01*(x01 - xm10);
		double y0m1_y0m1_mi_y10 = y0m1*(y0m1 - y10);
		double x01_x01_mi_x10 = x01*(x01 - x10);
		double x0m1_x0m1_mi_x01 = x0m1*(x0m1 - x01);
		double x0m1_x0m1_mi_x10 = x0m1*(x0m1 - x10);
		double x01_x10_mi_x01 = x01*(x10 - x01);
		double x01_x01_mi_x0m1 = x01*(x01 - x0m1);
		double xm10e2_y01_mi_y0m1 = xm10e2*(y01 - y0m1);
		double y0m1_y10_mi_y0m1 = y0m1*(y10 - y0m1);
		double y01_mi_y10_y0m1_mi_y10 = (y01 - y10)*(y0m1 - y10);
		double y01_mi_y0m1_y01_mi_ym10_y0m1_mi_ym10 = (y01 - y0m1)*(y01 - ym10)*(y0m1 - ym10);
		double y01_mi_y10_y01_mi_ym10 = (y01 - y10)*(y01 - ym10);
		double y0m1_mi_y10_y0m1_mi_ym10 = (y0m1 - y10)*(y0m1 - ym10);
		double x10_x10_mi_xm10 = x10*(x10 - xm10);
		double x0m1_x0m1_mi_xm10 = x0m1*(x0m1 - xm10);
		double x0m1_x10_mi_x0m1 = x0m1*(x10 - x0m1);
		double x10e2_y01_mi_ym10 = x10e2*(y01 - ym10);

		double DD = 1./((xm10*(x10*xm10*y01_y01_mi_y0m1*y0m1 + x10e2*y01*y0m1_y0m1_mi_y01 + y10*(x0m1e2*y01_y01_mi_y10 - x01_x01_mi_xm10*y0m1_y0m1_mi_y10 + x0m1*xm10*y01_y10_mi_y01)) + (x0m1*x10*(x10 - x0m1)*y01e2 + x01_x01_mi_x10*x10*y0m1e2 + x01*x0m1_x0m1_mi_x01*y10e2)*ym10 + (x10*(x0m1_x0m1_mi_x10*y01 + x01_x10_mi_x01*y0m1) + x01_x01_mi_x0m1*x0m1*y10)*ym10e2));
		double a10 = DD*(xm10e2*(f10*y01_y01_mi_y0m1*y0m1 + y10*(f01*y0m1_y0m1_mi_y10 + f0m1*y01_y10_mi_y01)) + fm10*(x10e2*y01*y0m1_y0m1_mi_y01 + y10*(x0m1e2*y01_y01_mi_y10 + x01e2*y0m1_y10_mi_y0m1)) + f00*(-xm10e2_y01_mi_y0m1*y01_mi_y10_y0m1_mi_y10 + x10e2*y01_mi_y0m1_y01_mi_ym10_y0m1_mi_ym10 - (x0m1e2*y01_mi_y10_y01_mi_ym10 - x01e2*y0m1_mi_y10_y0m1_mi_ym10)*(y10 - ym10)) + (-f10*x0m1e2*y01e2 + f0m1*x10e2*y01e2 + f10*x01e2*y0m1e2 - f01*x10e2*y0m1e2 - f0m1*x01e2*y10e2 + f01*x0m1e2*y10e2)*ym10 + (f10*x0m1e2*y01 - f0m1*x10e2*y01 - f10*x01e2*y0m1 + f01*x10e2*y0m1 + f0m1*x01e2*y10 - f01*x0m1e2*y10)*ym10e2);
		double a01 = (fm10*x10*(x0m1_x0m1_mi_x10*y01e2 + x01_x10_mi_x01*y0m1e2) + fm10*x01_x01_mi_x0m1*x0m1*y10e2 + xm10*(f10*x0m1*(xm10 - x0m1)*y01e2 + f10*x01_x01_mi_xm10*y0m1e2 + x10_x10_mi_xm10*(f0m1*y01e2 - f01*y0m1e2) + (f01*x0m1_x0m1_mi_xm10 + f0m1*x01*(xm10 - x01))*y10e2) + (f10*x01*x0m1_x0m1_mi_x01 + x10*(f0m1*x01_x01_mi_x10 + f01*x0m1_x10_mi_x0m1))*ym10e2 + f00*(x01_x01_mi_xm10*xm10*(-y0m1e2 + y10e2) + x10e2*(-(xm10*y01e2) - x01*y0m1e2 + xm10*y0m1e2 + x01*ym10e2) + x10*(xm10e2_y01_mi_y0m1*(y01 + y0m1) + x01e2*(y0m1e2 - ym10e2)) + x0m1e2*(xm10*(y01e2 - y10e2) + x01*(y10e2 - ym10e2) + x10*(ym10e2 - y01e2)) + x0m1*(xm10e2*(y10e2 - y01e2) + x10e2_y01_mi_ym10*(y01 + ym10) + x01e2*(ym10e2 - y10e2))))/
			(y10*(xm10*(x0m1*xm10*y01_y01_mi_y10 + x01_x01_mi_xm10*y0m1_y0m1_mi_y10 + x0m1e2*y01_y10_mi_y01) + x01_x01_mi_x0m1*x0m1*y10*ym10 + x01*x0m1_x0m1_mi_x01*ym10e2) + x10e2*(xm10*y01_y01_mi_y0m1*y0m1 + ym10*(x01*y0m1*(y0m1 - ym10) + x0m1*y01*(ym10 - y01))) + x10*(xm10e2*y01*y0m1_y0m1_mi_y01 + ym10*(x0m1e2*y01*(y01 - ym10) + x01e2*y0m1*(ym10 - y0m1))));
		double a20 = DD*(fm10*(x10*y01_y01_mi_y0m1*y0m1 + y10*(x01*y0m1_y0m1_mi_y10 + x0m1*y01_y10_mi_y01)) + xm10*(f10*y01*y0m1_y0m1_mi_y01 + y10*(f0m1*y01_y01_mi_y10 + f01*y0m1_y10_mi_y0m1)) + f00*(xm10*(y01 - y0m1)*y01_mi_y10_y0m1_mi_y10 - x10*y01_mi_y0m1_y01_mi_ym10_y0m1_mi_ym10 + (x0m1*y01_mi_y10_y01_mi_ym10 - x01*y0m1_mi_y10_y0m1_mi_ym10)*(y10 - ym10)) + (f10*x0m1*y01e2 - f0m1*x10*y01e2 - f10*x01*y0m1e2 + f01*x10*y0m1e2 + (f0m1*x01 - f01*x0m1)*y10e2)*ym10 + (f0m1*x10*y01 - f10*x0m1*y01 + f10*x01*y0m1 - f01*x10*y0m1 - f0m1*x01*y10 + f01*x0m1*y10)*ym10e2);
		double a02 = DD*(fm10*x10*(x0m1_x0m1_mi_x10*y01 + x01_x10_mi_x01*y0m1) + fm10*x01_x01_mi_x0m1*x0m1*y10 + xm10*(x10_x10_mi_xm10*(f0m1*y01 - f01*y0m1) + f10*(-x0m1e2*y01 + x0m1*xm10*y01 + x01_x01_mi_xm10*y0m1) + (f01*x0m1_x0m1_mi_xm10 + f0m1*x01*(xm10 - x01))*y10) + (f10*x01*x0m1_x0m1_mi_x01 + x10*(f0m1*x01_x01_mi_x10 + f01*x0m1_x10_mi_x0m1))*ym10 + f00*(-(x01_x01_mi_xm10*xm10*(y0m1 - y10)) + x10*(xm10e2_y01_mi_y0m1 + x01e2*(y0m1 - ym10)) + x10e2*(xm10*y0m1 - xm10*y01 - x01*y0m1 + x01*ym10) + x0m1e2*(xm10*y01 + x01*y10 - xm10*y10 - x01*ym10 + x10*(ym10 - y01)) + x0m1*(xm10e2*(y10 - y01) + x10e2_y01_mi_ym10 + x01e2*(ym10 - y10))));
		return x*(x*a20 + a10) + y*(y*a02 + a01) + f00;
	}

	static double Interp1dLinRel(double xr, double f0, double f1)
	{//simple linear interpolation for xr being relative argument (f = f0 at xr = 0, f = f1 at xr = 1)
		return f0 + xr*(f1 - f0);
	}

	static double Interp1dQuadVarRel(double xr, double hmdhp, double fm1, double f0, double fp1)
	{//simple quadratic interpolation for xr being relative argument (f = fm1 at xr = -hmdhp, f = f0 at xr = 0, f = fp1 at xr = 1)
		double buf = 1./((hmdhp + 1.)*hmdhp);
		double a1 = (f0 - fm1 + (fp1 - f0)*hmdhp*hmdhp)*buf;
		double a2 = (fm1 - f0 + (fp1 - f0)*hmdhp)*buf;
		return f0 + (a1 + a2*xr)*xr;
	}

	static double Interp1dCubVarRel(double xr, double hmdhp1, double hp2dhp1, double fm1, double f0, double fp1, double fp2)
	{//simple cubic interpolation for xr being relative argument (f = fm1 at xr = -hmdhp1, f = f0 at xr = 0, f = fp1 at xr = 1, f = fp2 at xr = hp2dhp1)
		double buf = 1./(hmdhp1*(1 + hmdhp1)*(hp2dhp1 - 1)*hp2dhp1*(hmdhp1 + hp2dhp1));
		double hmdhp1E2 = hmdhp1*hmdhp1;
		double hmdhp1E3 = hmdhp1E2*hmdhp1;
		double hp2dhp1E2 = hp2dhp1*hp2dhp1;
		double hp2dhp1E3 = hp2dhp1E2*hp2dhp1;
		double a1 = -((fp2-f0)*hmdhp1E2 + (fp2-f0)*hmdhp1E3 + (f0-fm1)*hp2dhp1E2 + (f0-fp1)*hmdhp1E3*hp2dhp1E2 + (fm1-f0)*hp2dhp1E3 + (f0-fp1)*hmdhp1E2*hp2dhp1E3)*buf;
		double buf1 = f0/(hmdhp1*hp2dhp1);
		double a2 = buf1*(hmdhp1 - hp2dhp1 - 1.) - (fp2*hmdhp1 - fp2*hmdhp1E3 + fm1*hp2dhp1 + fp1*hmdhp1E3*hp2dhp1 - fm1*hp2dhp1E3 - fp1*hmdhp1*hp2dhp1E3)*buf;
		double a3 = buf1 + (fp2*hmdhp1 + fp2*hmdhp1E2 + fm1*hp2dhp1 - fp1*hmdhp1E2*hp2dhp1 - fm1*hp2dhp1E2 - fp1*hmdhp1*hp2dhp1E2)*buf;
		return f0 + (a1 + (a2 + a3*xr)*xr)*xr;
	}

	static double InterpCubHalfStep(double* f, int i)
	{//interpolation by cubic poligon for a point in the middle between equidistant points for which the function f is defined 
		if(i < 0) return 0.0625*(5.*f[0] + 15.*f[1] - 5.*f[2] + f[3]);
		else if(i == 0) return 0.0625*(-f[0] + 9.*f[1] + 9.*f[2] - f[3]);
		else return 0.0625*(f[0] - 5.*f[1] + 15.*f[2] + 5.*f[3]);
	}

	static void CubicPlnCfs(double f1, double f2, double fpr1, double fpr2, double sStep, double* aa)
	{
		double f1mf2_d_s1ms2 = (f2 - f1)/sStep;
		aa[0] = f1;
		aa[1] = fpr1;
		aa[2] = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)/sStep;
		aa[3] = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)/(sStep*sStep);
	}

	static double Interp2D4pRel(double f00, double f10, double f01, double f11, double xr, double yr)
	{// assumes 0 <= xr <= 1; 0 <= yr <= 1; f00 := f|xr=0,yr=0; f10 := f|xr=1,yr=0;
		return (f11 + f00 - f10 - f01)*xr*yr + (f10 - f00)*xr + (f01 - f00)*xr + f00;
	}

	template<class T> static double InterpOnRegMesh2d(double x, double y, double x_min, double x_step, long nx, double y_min, double y_step, long ny, T* ar_f, char ord=3, long ix_per=1, long ix_ofst=0)
	{//OC20112018: "copied" from uti_math.py: uti_math.interp_2d

		if((x_step == 0) || (y_step == 0) || (ar_f == 0) || (ord < 1) || (ord > 3)) throw CAN_NOT_FIND_IND_FOR_INTERP;
		const double truncTol = 1.e-12; //to steer

		if(ord == 1) //bi-linear interpolation based on 4 points
		{
			long ix0 = (long)((x - x_min)/x_step + truncTol);
	        if(ix0 < 0) ix0 = 0;
			else if(ix0 >= nx - 1) ix0 = nx - 2;

			long ix1 = ix0 + 1;
			double tx = (x - (x_min + x_step*ix0))/x_step;

	        long iy0 = (long)((y - y_min)/y_step + truncTol);
			if(iy0 < 0) iy0 = 0;
			else if(iy0 >= ny - 1) iy0 = ny - 2;

			long iy1 = iy0 + 1;
			double ty = (y - (y_min + y_step*iy0))/y_step;

			long long nx_ix_per = nx*ix_per;
			long long iy0_nx_ix_per = iy0*nx_ix_per;
			long long iy1_nx_ix_per = iy1*nx_ix_per;
			long long ix0_ix_per_p_ix_ofst = ix0*ix_per + ix_ofst;
			long long ix1_ix_per_p_ix_ofst = ix1*ix_per + ix_ofst;

			double a00 = *(ar_f + (iy0_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f10 = *(ar_f + (iy0_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double f01 = *(ar_f + (iy1_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f11 = *(ar_f + (iy1_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double a10 = f10 - a00;
			double a01 = f01 - a00;
			double a11 = a00 - f01 - f10 + f11;
			return a00 + tx*(a10 + ty*a11) + ty*a01;
		}
		else if(ord == 2) //bi-quadratic interpolation based on 6 points
		{
			long ix0 = (long)((x - x_min)/x_step + truncTol);
			if(ix0 < 1) ix0 = 1;
			else if(ix0 >= nx - 1) ix0 = nx - 2;

			long ixm1 = ix0 - 1;
			long ix1 = ix0 + 1;
			double tx = (x - (x_min + x_step*ix0))/x_step;

			long iy0 = (long)((y - y_min)/y_step + truncTol);
			if(iy0 < 1) iy0 = 1;
			else if(iy0 >= ny - 1) iy0 = ny - 2;

			long iym1 = iy0 - 1;
			long iy1 = iy0 + 1;
			double ty = (y - (y_min + y_step*iy0))/y_step;

			long long nx_ix_per = nx*ix_per;
			long long iym1_nx_ix_per = iym1*nx_ix_per;
			long long iy0_nx_ix_per = iy0*nx_ix_per;
			long long iy1_nx_ix_per = iy1*nx_ix_per;
			long long ixm1_ix_per_p_ix_ofst = ixm1*ix_per + ix_ofst;
			long long ix0_ix_per_p_ix_ofst = ix0*ix_per + ix_ofst;
			long long ix1_ix_per_p_ix_ofst = ix1*ix_per + ix_ofst;

			double fm10 = *(ar_f + (iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst));
			double a00 = *(ar_f + (iy0_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f10 = *(ar_f + (iy0_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double f0m1 = *(ar_f + (iym1_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f01 = *(ar_f + (iy1_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f11 = *(ar_f + (iy1_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double a10 = 0.5*(f10 - fm10);
			double a01 = 0.5*(f01 - f0m1);
			double a11 = a00 - f01 - f10 + f11;
			double a20 = 0.5*(f10 + fm10) - a00;
			double a02 = 0.5*(f01 + f0m1) - a00;
			return a00 + tx*(a10 + tx*a20 + ty*a11) + ty*(a01 + ty*a02);
		}
		else if(ord == 3) //bi-cubic interpolation based on 12 points
		{
			long ix0 = (long)((x - x_min)/x_step + truncTol);
			if(ix0 < 1) ix0 = 1;
			else if(ix0 >= nx - 2) ix0 = nx - 3;

			long ixm1 = ix0 - 1;
			long ix1 = ix0 + 1;
			long ix2 = ix0 + 2;
			double tx = (x - (x_min + x_step*ix0))/x_step;

			long iy0 = (long)((y - y_min)/y_step + truncTol);
			if(iy0 < 1) iy0 = 1;
			else if(iy0 >= ny - 2) iy0 = ny - 3;

			long iym1 = iy0 - 1;
			long iy1 = iy0 + 1;
			long iy2 = iy0 + 2;
			double ty = (y - (y_min + y_step*iy0))/y_step;

			long long nx_ix_per = nx*ix_per;
			long long iym1_nx_ix_per = iym1*nx_ix_per;
			long long iy0_nx_ix_per = iy0*nx_ix_per;
			long long iy1_nx_ix_per = iy1*nx_ix_per;
			long long iy2_nx_ix_per = iy2*nx_ix_per;
			long long ixm1_ix_per_p_ix_ofst = ixm1*ix_per + ix_ofst;
			long long ix0_ix_per_p_ix_ofst = ix0*ix_per + ix_ofst;
			long long ix1_ix_per_p_ix_ofst = ix1*ix_per + ix_ofst;
			long long ix2_ix_per_p_ix_ofst = ix2*ix_per + ix_ofst;

			double f0m1 = *(ar_f + (iym1_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f1m1 = *(ar_f + (iym1_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double fm10 = *(ar_f + (iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst));
			double a00 = *(ar_f + (iy0_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f10 = *(ar_f + (iy0_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double f20 = *(ar_f + (iy0_nx_ix_per + ix2_ix_per_p_ix_ofst));
			double fm11 = *(ar_f + (iy1_nx_ix_per + ixm1_ix_per_p_ix_ofst));
			double f01 = *(ar_f + (iy1_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f11 = *(ar_f + (iy1_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double f21 = *(ar_f + (iy1_nx_ix_per + ix2_ix_per_p_ix_ofst));
			double f02 = *(ar_f + (iy2_nx_ix_per + ix0_ix_per_p_ix_ofst));
			double f12 = *(ar_f + (iy2_nx_ix_per + ix1_ix_per_p_ix_ofst));
			double a10 = -0.5*a00 + f10 - f20/6. - fm10/3.;
			double a01 = -0.5*a00 + f01 - f02/6. - f0m1/3.;
			double a11 = -0.5*(f01 + f10) + (f02 - f12 + f20 - f21)/6. + (f0m1 - f1m1 + fm10 - fm11)/3. + f11;
			double a20 = -a00 + 0.5*(f10 + fm10);
			double a02 = -a00 + 0.5*(f01 + f0m1);
			double a21 = a00 - f01 + 0.5*(f11 - f10 - fm10 + fm11);
			double a12 = a00 - f10 + 0.5*(f11 - f01 - f0m1 + f1m1);
			double a30 = 0.5*(a00 - f10) + (f20 - fm10)/6.;
			double a03 = 0.5*(a00 - f01) + (f02 - f0m1)/6.;
			double a31 = 0.5*(f01 + f10 - f11 - a00) + (f21 + fm10 - f20 - fm11)/6.;
			double a13 = 0.5*(f10 - f11 - a00 + f01) + (f0m1 + f12 - f02 - f1m1)/6.;
			return a00 + tx*(a10 + tx*(a20 + tx*(a30 + ty*a31) + ty*a21) + ty*a11) + ty*(a01 + ty*(a02 + ty*(a03 + tx*a13) + tx*a12));
		}
		else return 0;
	}
};

//-------------------------------------------------------------------------

#endif
