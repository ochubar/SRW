/************************************************************************//**
 * File: srfft.h
 * Description: Auxiliary utilities to work with FFTW library (header)
 * Project: 
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMFFT_H
#define __GMFFT_H

#include "fftw.h"

//#ifdef __IGOR_PRO__
//#include "srigintr.h"
//#endif
//#include "srercode.h"

//#include <cmath>
#include <math.h>

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

#ifndef MEMORY_ALLOCATION_FAILURE
#define MEMORY_ALLOCATION_FAILURE 8 + 10000 //in line with SRW
#endif
#ifndef ERROR_IN_FFT
#define ERROR_IN_FFT 40 + 10000
#endif

//*************************************************************************

class CGenMathFFT //{
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{//OC01052013
	double a2c, a4c, a6c, a8c, a10c, a12c;
	double a3s, a5s, a7s, a9s, a11s, a13s;

protected:

	static long GoodNumbers[];
	static long LenGoodNumbers;
	static long GoodNum100s[];
	static long LenGoodNum100s;
	static long GoodNum1000s[];
	static long LenGoodNum1000s;
	static long GoodNum10000s[];
	static long LenGoodNum10000s;

public:

	double HalfPI, PI, TwoPI, ThreePIdTwo, One_dTwoPI; // Constants

	CGenMathFFT()
	{
		HalfPI = 1.5707963267949;
		PI = 3.141592653590;
		TwoPI = 6.2831853071796;
		ThreePIdTwo = 4.7123889803847;
		One_dTwoPI = 0.1591549430919;
		a2c = -0.5; a4c = 0.041666666666667; a6c = -0.0013888888888889; a8c = 0.000024801587301587; a10c = -2.755731922E-07;
		a3s = -0.16666666666667; a5s = 0.0083333333333333; a7s = -0.0001984126984127; a9s = 2.755731922E-06; a11s = -2.505210839E-08;
	}

	void CosAndSin(double x, float& Cos, float& Sin)
	{
		x -= TwoPI*int(x*One_dTwoPI);
		if(x < 0.) x += TwoPI;

		char ChangeSign=0;
		if(x > ThreePIdTwo) x -= TwoPI;
		else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

		double xe2 = x*x;
		Cos = float(1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c)))));
		Sin = float(x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s))))));
		if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
	}
	void NextCorrectNumberForFFT(long&);
};

//*************************************************************************

struct CGenMathFFT2DInfo {
	float* pData;
	char Dir; // >0: forward; <0: backward
	double xStep, yStep, xStart, yStart;
	double xStepTr, yStepTr, xStartTr, yStartTr;
	long Nx, Ny;

	long howMany; //OC151014
	long iStride, iDist; //OC151014
	//From FFTW 2.1.5 Tutorial
	//iStride and iDist describe the input array(s). 
	//There are howMany multi-dimensional input arrays; the first one is pointed to by in (= pData), 
	//the second one is pointed to by in + iDist, and so on, up to in + (howMany - 1) * iDist. 
	//Each multi-dimensional input array consists of complex numbers (see Section Data Types), 
	//stored in row-major format (see Section Multi-dimensional Array Format), which are not necessarily contiguous in memory. 
	//Specifically, in[0] is the first element of the first array, in[istride] is the second element of the first array, and so on. 
	//In general, the i-th element of the j-th input array will be in position in[i * istride + j * idist]. 
	//Note that, here, i refers to an index into the row-major format for the multi-dimensional array, rather than an index in any particular dimension. 
	//In-place transforms:  For plans created with the FFTW_IN_PLACE option, the transform is computed in-place--the output is returned in the in array, 
	//using the same strides, etcetera, as were used in the input. 

	char UseGivenStartTrValues;
	double ExtraMult; //OC20112017

	CGenMathFFT2DInfo() 
	{ 
		howMany = 1; iStride = 1; iDist = 0; //OC151014
		UseGivenStartTrValues = 0;
		ExtraMult = 1.; //OC20112017
	}
};

//*************************************************************************

class CGenMathFFT2D : public CGenMathFFT {

	long Nx, Ny;
	long HalfNx, HalfNy;
	char NeedsShiftBeforeX, NeedsShiftBeforeY, NeedsShiftAfterX, NeedsShiftAfterY;
	float *ArrayShiftX, *ArrayShiftY;

public:
	CGenMathFFT2D()
	{
		NeedsShiftBeforeX = NeedsShiftBeforeY = NeedsShiftAfterX = NeedsShiftAfterY = 0;
	}

	int Make2DFFT(CGenMathFFT2DInfo&);
	int AuxDebug_TestFFT_Plans();

	void SetupLimitsTr(CGenMathFFT2DInfo& FFT2DInfo)
	{// Modify this if Make2DFFT is modified !
		Nx = FFT2DInfo.Nx; Ny = FFT2DInfo.Ny; 
		HalfNx = (Nx >> 1); HalfNy = (Ny >> 1);

		double xStartTr = -0.5/FFT2DInfo.xStep;
		FFT2DInfo.xStepTr = -xStartTr/HalfNx;

		double yStartTr = -0.5/FFT2DInfo.yStep;
		FFT2DInfo.yStepTr = -yStartTr/HalfNy;

		if(!FFT2DInfo.UseGivenStartTrValues)
		{
			FFT2DInfo.xStartTr = xStartTr;
			FFT2DInfo.yStartTr = yStartTr;
		}
	}

	void FillArrayShift(char x_or_y, double t0, double tStep)
	{
		float* tArrayShift;
		long N;
		if(x_or_y == 'x') { tArrayShift = ArrayShiftX; N = Nx;}
		else { tArrayShift = ArrayShiftY; N = Ny;}

		float *tp = tArrayShift + N;
		*tp = 1.; *(tp+1) = 0.; tp += 2;
		float *tm = tp - 4;
		
		double t0TwoPI = t0*TwoPI;
		double q = tStep;
		long HalfN = N >> 1;
		for(int i=0; i<HalfN - 1; i++)
		{
			CosAndSin(q*t0TwoPI, *tp, *(tp+1));
			*tm = *tp; *(tm+1) = -(*(tp+1));
			tp += 2; tm -= 2; q += tStep;
		}
		CosAndSin(-q*t0TwoPI, *tm, *(tm+1));
	}

	void RotateDataAfter2DFFT(FFTW_COMPLEX* pAfterFFT)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
	 //OC281117: Consider combining RotateDataAfter2DFFT, RepairSignAfter2DFFT, NormalizeDataAfter2DFFT
		//long HalfNyNx = HalfNy*Nx;
		long long HalfNyNx = ((long long)HalfNy)*((long long)Nx);
		FFTW_COMPLEX *t1 = pAfterFFT, *t2 = pAfterFFT + (HalfNyNx + HalfNx);
	    FFTW_COMPLEX *t3 = pAfterFFT + HalfNx, *t4 = pAfterFFT + HalfNyNx;
		FFTW_COMPLEX Buf;
		for(long jj=0; jj<HalfNy; jj++)
		{
			for(long ii=0; ii<HalfNx; ii++)
			{
				Buf = *t1; *(t1++) = *t2; *(t2++) = Buf;
				Buf = *t3; *(t3++) = *t4; *(t4++) = Buf;
			}
			t1 += HalfNx; t2 += HalfNx; t3 += HalfNx; t4 += HalfNx;
		}
	}

	void RepairSignAfter2DFFT(FFTW_COMPLEX* pAfterFFT)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
		FFTW_COMPLEX *t = pAfterFFT;
		FFTW_REAL sx0 = 1., sy0 = 1., s;
		for(long iy=0; iy<Ny; iy++)
		{
			s = sy0*sx0;
			for(long ix=0; ix<Nx; ix++)
			{
				t->re *= s; (t++)->im *= s; s = -s;
			}
			sy0 = -sy0;
		}
	}

	void NormalizeDataAfter2DFFT(FFTW_COMPLEX* pAfterFFT, double Mult)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
		//long NxNy = Nx*Ny;
		long long NxNy = ((long long)Nx)*((long long)Ny);
		FFTW_COMPLEX *t = pAfterFFT;
		//for(long i=0; i<NxNy; i++)
		for(long long i=0; i<NxNy; i++)
		{
			t->re *= (FFTW_REAL)Mult; (t++)->im *= (FFTW_REAL)Mult;
		}
	}

	void TreatShifts(FFTW_COMPLEX* pData)
	{
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		char NeedsShiftY = NeedsShiftBeforeY || NeedsShiftAfterY;

		FFTW_COMPLEX *t = pData;
		float *tShiftY = ArrayShiftY;
		float MultY_Re = 1., MultY_Im = 0., MultX_Re = 1., MultX_Im = 0.;
		float MultRe, MultIm;

		for(long iy=0; iy<Ny; iy++)
		{
			if(NeedsShiftY) { MultY_Re = *(tShiftY++); MultY_Im = *(tShiftY++);}
			float *tShiftX = ArrayShiftX;
			for(long ix=0; ix<Nx; ix++)
			{
				if(NeedsShiftX) 
				{ 
					MultX_Re = *(tShiftX++); MultX_Im = *(tShiftX++);
					if(NeedsShiftY)
					{
						MultRe = MultX_Re*MultY_Re - MultX_Im*MultY_Im;
						MultIm = MultX_Re*MultY_Im + MultX_Im*MultY_Re;
					}
					else
					{
						MultRe = MultX_Re; MultIm = MultX_Im;
					}
				}
				else
				{
					MultRe = MultY_Re; MultIm = MultY_Im;
				}

				float NewRe = t->re*MultRe - t->im*MultIm;
				float NewIm = t->re*MultIm + t->im*MultRe;

				t->re = NewRe;
				(t++)->im = NewIm;
			}
		}
	}
};

//*************************************************************************

struct CGenMathFFT1DInfo {
	float *pInData, *pOutData;
	char Dir; // >0: forward; <0: backward
	double xStep, xStart;
	double xStepTr, xStartTr;
	long Nx;
	long HowMany;
	char UseGivenStartTrValue;
	double MultExtra;

	char TreatSharpEdges;
	double LeftSharpEdge, RightSharpEdge;
	char ApplyAutoShiftAfter;

	CGenMathFFT1DInfo() 
	{ 
		HowMany = 1; UseGivenStartTrValue = 0;
		TreatSharpEdges = 0;
		MultExtra = 1.;
		ApplyAutoShiftAfter = 1;
	}
};

//*************************************************************************

struct CGenMathAuxDataForSharpEdgeCorr1D {

	float *ExpArrSt, *ExpArrFi;
	double dSt, dFi, d;
	long iSt, iFi;

	char WasSetUp;

	CGenMathAuxDataForSharpEdgeCorr1D()
	{
		Initialize();
	}

	void Initialize()
	{
		ExpArrSt = ExpArrFi = 0;
		dSt = dFi = d = 0.;
		iSt = iFi = 0;
		WasSetUp = 0;
	}

	void Dispose()
	{
		if(ExpArrSt != 0) delete[] ExpArrSt;
		if(ExpArrFi != 0) delete[] ExpArrFi;

		Initialize();
	}
};

//*************************************************************************

class CGenMathFFT1D : public CGenMathFFT {

	long Nx;
	long HalfNx;
	char NeedsShiftBeforeX, NeedsShiftAfterX;
	float *m_ArrayShiftX;

public:
	CGenMathFFT1D()
	{
		NeedsShiftBeforeX = NeedsShiftAfterX = 0;
	}

	int Make1DFFT(CGenMathFFT1DInfo&);
	int Make1DFFT_InPlace(CGenMathFFT1DInfo& FFT1DInfo);

	void SetupLimitsTr(CGenMathFFT1DInfo& FFT1DInfo)
	{ // Modify this if Make1DFFT is modified !
		Nx = FFT1DInfo.Nx;
		HalfNx = (Nx >> 1);

		double xStartTr = -0.5/FFT1DInfo.xStep;
		FFT1DInfo.xStepTr = -xStartTr/HalfNx;

		if(!FFT1DInfo.UseGivenStartTrValue)
		{
			FFT1DInfo.xStartTr = xStartTr;
		}
	}

	void FillArrayShift(double t0, double tStep)
	{
		float *tArrayShift = m_ArrayShiftX;
		long N = Nx;

		float *tp = tArrayShift + N;
		*tp = 1.; *(tp+1) = 0.; tp += 2;
		float *tm = tp - 4;

		double t0TwoPI = t0*TwoPI;
		double q = tStep;
		long HalfN = N >> 1;

		for(int i=0; i<HalfN - 1; i++)
		{
			CosAndSin(q*t0TwoPI, *tp, *(tp+1));
			*tm = *tp; *(tm+1) = -(*(tp+1));
			tp += 2; tm -= 2; q += tStep;
		}
		CosAndSin(-q*t0TwoPI, *tm, *(tm+1));
	}

	void TreatShift(FFTW_COMPLEX* pData, long HowMany)
	{
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		if(!NeedsShiftX) return;

		FFTW_COMPLEX *t = pData;
		float *tShiftX = m_ArrayShiftX;

		for(long ix=0; ix<Nx; ix++)
		{
			float MultX_Re = *(tShiftX++), MultX_Im = *(tShiftX++);

			FFTW_COMPLEX *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				float NewRe = tMany->re*MultX_Re - tMany->im*MultX_Im;
				float NewIm = tMany->re*MultX_Im + tMany->im*MultX_Re;

				tMany->re = NewRe; tMany->im = NewIm;
				tMany += Nx;
			}
		}
	}

	void RepairSignAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany)
	{// Assumes Nx even !
		FFTW_COMPLEX *t = pAfterFFT;
		int s = 1;
		for(long ix=0; ix<Nx; ix++)
		{
			if(s < 0)
			{
				FFTW_COMPLEX *tMany = t;
				for(long k=0; k<HowMany; k++)
				{
					tMany->re = -tMany->re; tMany->im = -tMany->im;
					tMany += Nx;
				}
			}
			t++; s = -s;
		}
	}

	void RotateDataAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany)
	{// Assumes Nx even !
		FFTW_COMPLEX *t1 = pAfterFFT, *t2 = pAfterFFT + HalfNx;
		FFTW_COMPLEX Buf;
		for(long ix=0; ix<HalfNx; ix++)
		{
			FFTW_COMPLEX *t1Many = t1++, *t2Many = t2++;
			for(long k=0; k<HowMany; k++)
			{
				Buf = *t1Many; *t1Many = *t2Many; *t2Many = Buf;
				t1Many += Nx; t2Many += Nx; 
			}
		}
	}

	void NormalizeDataAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany, double Mult)
	{// Assumes Nx even !
		FFTW_COMPLEX *t = pAfterFFT;
		for(long ix=0; ix<Nx; ix++)
		{
			FFTW_COMPLEX *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				tMany->re *= (FFTW_REAL)Mult; tMany->im *= (FFTW_REAL)Mult;
				tMany += Nx;
			}
		}
	}

	int SetupAuxDataForSharpEdgeCorr(CGenMathFFT1DInfo&, CGenMathAuxDataForSharpEdgeCorr1D&);
	void MakeSharpEdgeCorr(CGenMathFFT1DInfo&, CGenMathAuxDataForSharpEdgeCorr1D&);

	void SetupSharpEdgeExpCorrArray(float* pCmpData, long AmOfPt, double x, double qStart, double qStep)
	{
		const double TwoPi = 6.28318530717959;
		double TwoPiX = TwoPi*x;
		double q = qStart;
		float *tCmpData = pCmpData;
		for(long i=0; i<AmOfPt; i++)
		{
			double Arg = TwoPiX*q;
			float Co, Si;
			CosAndSin(Arg, Co, Si);
			*(tCmpData++) = Co; *(tCmpData++) = -Si;
			q += qStep;
		}
	}
	int ProcessSharpEdges(CGenMathFFT1DInfo& FFT1DInfo)
	{
		CGenMathAuxDataForSharpEdgeCorr1D AuxDataForSharpEdgeCorr;
		int result = SetupAuxDataForSharpEdgeCorr(FFT1DInfo, AuxDataForSharpEdgeCorr);
		if(result) return result;
		MakeSharpEdgeCorr(FFT1DInfo, AuxDataForSharpEdgeCorr);
		AuxDataForSharpEdgeCorr.Dispose();
		return 0;
	}
};

//*************************************************************************

#endif
