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

#ifdef _FFTW3 //OC28012019
#include "fftw3.h"
#else
#include "fftw.h"
#endif

//#ifdef __IGOR_PRO__
//#include "srigintr.h"
//#endif
//#include "srercode.h"

//#include <cmath>
#include <math.h>

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

#ifdef _WITH_OMP //OC31102018: Pre-processor definition for compiling SRW with OpenMP library
#include "omp.h"
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
	void CosAndSin(double x, double& Cos, double& Sin) //OC02022019
	{
		//x -= TwoPI*int(x*One_dTwoPI);
		x -= TwoPI*((long long)(x*One_dTwoPI));

		if(x < 0.) x += TwoPI;

		char ChangeSign=0;
		if(x > ThreePIdTwo) x -= TwoPI;
		else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

		double xe2 = x*x;
		Cos = 1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c))));
		Sin = x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s)))));
		if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
	}

	//void NextCorrectNumberForFFT(long long&); //OC26042019
	void NextCorrectNumberForFFT(long&);
};

//*************************************************************************

struct CGenMathFFT2DInfo {
	float* pData;
	double* pdData; //OC31012019

	char Dir; // >0: forward; <0: backward
	double xStep, yStep, xStart, yStart;
	double xStepTr, yStepTr, xStartTr, yStartTr;
	long Nx, Ny;
	//long long Nx, Ny;

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

		pData = 0; //OC31012019
		pdData = 0;
	}
};

//*************************************************************************

class CGenMathFFT2D : public CGenMathFFT {

	long Nx, Ny;
	long HalfNx, HalfNy;
	//long long Nx, Ny;
	//long long HalfNx, HalfNy;
	char NeedsShiftBeforeX, NeedsShiftBeforeY, NeedsShiftAfterX, NeedsShiftAfterY;
	//float *ArrayShiftX, *ArrayShiftY;
	float *m_ArrayShiftX, *m_ArrayShiftY; //OC02022019
	double *m_dArrayShiftX, *m_dArrayShiftY; 

public:
	CGenMathFFT2D()
	{
		NeedsShiftBeforeX = NeedsShiftBeforeY = NeedsShiftAfterX = NeedsShiftAfterY = 0;
	}

	//int Make2DFFT(CGenMathFFT2DInfo&);
	//Modification by S.Yakubov for parallelizing SRW via OpenMP:
#ifdef _FFTW3 //28012019
	int Make2DFFT(CGenMathFFT2DInfo&, fftwf_plan* pPrecreatedPlan2DFFT=0, fftw_plan* pdPrecreatedPlan2DFFT=0); //OC02022019
	//int Make2DFFT(CGenMathFFT2DInfo&, fftwf_plan* pPrecreatedPlan2DFFT=0);
#else
	int Make2DFFT(CGenMathFFT2DInfo&, fftwnd_plan* pPrecreatedPlan2DFFT=0); //OC27102018
#endif

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

	template <class T> void FillArrayShift(char x_or_y, double t0, double tStep, T* arShift) //OC02022019
	//void FillArrayShift(char x_or_y, double t0, double tStep)
	{
		T* tArrayShift = arShift;
		//float* tArrayShift;
		//long N;
		long N = (x_or_y == 'x')? Nx : Ny;
		//if(x_or_y == 'x') { tArrayShift = m_ArrayShiftX; N = Nx;}
		//else { tArrayShift = m_ArrayShiftY; N = Ny;}

		T *tp = tArrayShift + N;
		//float *tp = tArrayShift + N;
		*tp = 1.; *(tp+1) = 0.; tp += 2;
		T *tm = tp - 4;
		//float *tm = tp - 4;
		
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

#ifdef _FFTW3 //OC29012019
	template <class T> void RotateDataAfter2DFFT(T* pAfterFFT)
	//void RotateDataAfter2DFFT(fftwf_complex* pAfterFFT)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
	 //OC281117: Consider combining RotateDataAfter2DFFT, RepairSignAfter2DFFT, NormalizeDataAfter2DFFT
		//long HalfNyNx = HalfNy*Nx;
		long long HalfNyNx = ((long long)HalfNy)*((long long)Nx);

		//fftwf_complex *t1 = pAfterFFT, *t2 = pAfterFFT + (HalfNyNx + HalfNx);
		//fftwf_complex *t3 = pAfterFFT + HalfNx, *t4 = pAfterFFT + HalfNyNx;
		//fftwf_complex Buf;
		T *t1 = pAfterFFT, *t2 = pAfterFFT + (HalfNyNx + HalfNx);
	    T *t3 = pAfterFFT + HalfNx, *t4 = pAfterFFT + HalfNyNx;
		T Buf;

		for(long jj=0; jj<HalfNy; jj++)
		{
			for(long ii=0; ii<HalfNx; ii++)
			{
				Buf[0] = (*t1)[0]; Buf[1] = (*t1)[1]; //Buf = *t1;
				(*t1)[0] = (*t2)[0]; (*(t1++))[1] = (*t2)[1]; //*(t1++) = *t2;
				(*t2)[0] = Buf[0]; (*(t2++))[1] = Buf[1]; //*(t2++) = Buf;

				Buf[0] = (*t3)[0]; Buf[1] = (*t3)[1]; //Buf = *t3; 
				(*t3)[0] = (*t4)[0]; (*(t3++))[1] = (*t4)[1]; //*(t3++) = *t4; 
				(*t4)[0] = Buf[0]; (*(t4++))[1] = Buf[1]; //*(t4++) = Buf;
			}
			t1 += HalfNx; t2 += HalfNx; t3 += HalfNx; t4 += HalfNx;
		}
	}
#else
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
#endif

#ifdef _FFTW3 //OC29012019
	void RepairSignAfter2DFFT(fftwf_complex* pAfterFFT)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
		fftwf_complex *t = pAfterFFT;
		float sx0 = 1., sy0 = 1., s;
		for(long iy=0; iy<Ny; iy++)
		{
			s = sy0*sx0;
			for(long ix=0; ix<Nx; ix++)
			{
				(*t)[0] *= s; (*(t++))[1] *= s; s = -s;
			}
			sy0 = -sy0;
		}
	}
	void RepairSignAfter2DFFT(fftw_complex* pAfterFFT)
	{// Assumes Nx, Ny even !
	 //OC281117: Make it work for odd Nx, Ny as well!
		fftw_complex *t = pAfterFFT;
		double sx0 = 1., sy0 = 1., s;
		for(long iy=0; iy<Ny; iy++)
		{
			s = sy0*sx0;
			for(long ix=0; ix<Nx; ix++)
			{
				(*t)[0] *= s; (*(t++))[1] *= s; s = -s;
			}
			sy0 = -sy0;
		}
	}
#else
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
#endif

#ifdef _FFTW3 //OC29012019
	void NormalizeDataAfter2DFFT(fftwf_complex* pAfterFFT, double Mult)
	{// Assumes Nx, Ny even !
	 //OC281117: To make it work for odd Nx, Ny as well in the future!
		float fMult = (float)Mult;
		long long NxNy = ((long long)Nx)*((long long)Ny);
		fftwf_complex *t = pAfterFFT;
		for(long long i=0; i<NxNy; i++)
		{
			(*t)[0] *= fMult; (*(t++))[1] *= fMult; 
		}
	}
	void NormalizeDataAfter2DFFT(fftw_complex* pAfterFFT, double Mult)
	{// Assumes Nx, Ny even !
	 //OC281117: To make it work for odd Nx, Ny as well in the future!
		long long NxNy = ((long long)Nx)*((long long)Ny);
		fftw_complex *t = pAfterFFT;
		for(long long i=0; i<NxNy; i++)
		{
			(*t)[0] *= Mult; (*(t++))[1] *= Mult; 
		}
	}
#else
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
#endif

#ifdef _FFTW3 //OC29012019
	void TreatShifts(fftwf_complex* pData)
	{
		fftwf_complex *t = pData;
#else
	void TreatShifts(FFTW_COMPLEX* pData)
	{
		FFTW_COMPLEX *t = pData;
#endif
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		char NeedsShiftY = NeedsShiftBeforeY || NeedsShiftAfterY;

		float *tShiftY = m_ArrayShiftY;
		float MultY_Re = 1., MultY_Im = 0., MultX_Re = 1., MultX_Im = 0.;
		float MultRe, MultIm;

		for(long iy=0; iy<Ny; iy++)
		{
			if(NeedsShiftY) { MultY_Re = *(tShiftY++); MultY_Im = *(tShiftY++);}
			float *tShiftX = m_ArrayShiftX;
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

#ifdef _FFTW3 //OC29012019
				float NewRe = (*t)[0]*MultRe - (*t)[1]*MultIm;
				float NewIm = (*t)[0]*MultIm + (*t)[1]*MultRe;
				(*t)[0] = NewRe;
				(*(t++))[1] = NewIm;
#else
				float NewRe = t->re*MultRe - t->im*MultIm;
				float NewIm = t->re*MultIm + t->im*MultRe;
				t->re = NewRe;
				(t++)->im = NewIm;
#endif
			}
		}
	}
#ifdef _FFTW3 //OC02022019
	void TreatShifts(fftw_complex* pData)
	{
		fftw_complex *t = pData;
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		char NeedsShiftY = NeedsShiftBeforeY || NeedsShiftAfterY;

		double *tShiftY = m_dArrayShiftY;
		double MultY_Re = 1., MultY_Im = 0., MultX_Re = 1., MultX_Im = 0.;
		double MultRe, MultIm;

		for(long iy=0; iy<Ny; iy++)
		{
			if(NeedsShiftY) { MultY_Re = *(tShiftY++); MultY_Im = *(tShiftY++);}
			double *tShiftX = m_dArrayShiftX;
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

				double NewRe = (*t)[0]*MultRe - (*t)[1]*MultIm;
				double NewIm = (*t)[0]*MultIm + (*t)[1]*MultRe;
				(*t)[0] = NewRe;
				(*(t++))[1] = NewIm;
			}
		}
	}
#endif
};

//*************************************************************************

struct CGenMathFFT1DInfo {
	float *pInData, *pOutData;
	double *pdInData, *pdOutData; //OC31012019

	char Dir; // >0: forward; <0: backward
	double xStep, xStart;
	double xStepTr, xStartTr;
	long Nx;
	//long long Nx;
	long HowMany;
	//long long HowMany;
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

		pInData = 0; //OC31012019
		pOutData = 0;
		pdInData = 0;
		pdOutData = 0;
	}
};

//*************************************************************************

struct CGenMathAuxDataForSharpEdgeCorr1D {

	float *ExpArrSt, *ExpArrFi;
	double *dExpArrSt, *dExpArrFi;

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
		dExpArrSt = dExpArrFi = 0;

		dSt = dFi = d = 0.;
		iSt = iFi = 0;
		WasSetUp = 0;
	}

	void Dispose()
	{
		if(ExpArrSt != 0) delete[] ExpArrSt;
		if(ExpArrFi != 0) delete[] ExpArrFi;

		if(dExpArrSt != 0) delete[] dExpArrSt;
		if(dExpArrFi != 0) delete[] dExpArrFi;

		Initialize();
	}
};

//*************************************************************************

class CGenMathFFT1D : public CGenMathFFT {

	long Nx;
	long HalfNx;
	//long long Nx;
	//long long HalfNx;
	char NeedsShiftBeforeX, NeedsShiftAfterX;
	float *m_ArrayShiftX;
	double *m_dArrayShiftX; //OC02022019

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

	template <class T> void FillArrayShift(double t0, double tStep, T* arShiftX) //OC02022019
	//void FillArrayShift(double t0, double tStep)
	{
		//float *tArrayShift = m_ArrayShiftX;
		T *tArrayShift = arShiftX; //OC02022019
		long N = Nx;

		//float *tp = tArrayShift + N;
		T *tp = tArrayShift + N; //OC02022019
		*tp = 1.; *(tp+1) = 0.; tp += 2;
		//float *tm = tp - 4;
		T *tm = tp - 4;

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

#ifdef _FFTW3 //OC29012019
	void TreatShift(fftwf_complex* pData, long HowMany)
	{
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		if(!NeedsShiftX) return;

		fftwf_complex *t = pData;
		float *tShiftX = m_ArrayShiftX;

		for(long ix=0; ix<Nx; ix++)
		{
			float MultX_Re = *(tShiftX++), MultX_Im = *(tShiftX++);

			fftwf_complex *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				float NewRe = (*tMany)[0]*MultX_Re - (*tMany)[1]*MultX_Im;
				float NewIm = (*tMany)[0]*MultX_Im + (*tMany)[1]*MultX_Re;
				(*tMany)[0] = NewRe; (*tMany)[1] = NewIm;
				tMany += Nx;
			}
		}
	}
	void TreatShift(fftw_complex* pData, long HowMany) //OC02022019
	{
		char NeedsShiftX = NeedsShiftBeforeX || NeedsShiftAfterX;
		if(!NeedsShiftX) return;

		fftw_complex *t = pData;
		double *tShiftX = m_dArrayShiftX;

		for(long ix=0; ix<Nx; ix++)
		{
			double MultX_Re = *(tShiftX++), MultX_Im = *(tShiftX++);

			fftw_complex *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				double NewRe = (*tMany)[0]*MultX_Re - (*tMany)[1]*MultX_Im;
				double NewIm = (*tMany)[0]*MultX_Im + (*tMany)[1]*MultX_Re;
				(*tMany)[0] = NewRe; (*tMany)[1] = NewIm;
				tMany += Nx;
			}
		}
	}

#else
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
#endif

#ifdef _FFTW3 //OC29012019
	template <class T> void RepairSignAfter1DFFT(T* pAfterFFT, long HowMany) //OC02022019
	//void RepairSignAfter1DFFT(fftwf_complex* pAfterFFT, long HowMany)
	{// Assumes Nx even ! - to be improved
		//OC27102018
		//SY: optimized, adopt for OpenMP
#ifdef _WITH_OMP
		#pragma omp parallel for
#endif
		for(long ix=1; ix<Nx; ix+=2)
		{
			//FFTW_COMPLEX *t = pAfterFFT + ix;
			//FFTW_COMPLEX *tMany = t;
			//OC27102018
			//fftwf_complex *tMany = pAfterFFT + ix;
			T *tMany = pAfterFFT + ix;
			for(long k=0; k<HowMany; k++)
			{
				(*tMany)[0] = -(*tMany)[0]; (*tMany)[1] = -(*tMany)[1];
				tMany += Nx;
			}
		}
	}
#else
	void RepairSignAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany)
	{// Assumes Nx even !
		//FFTW_COMPLEX *t = pAfterFFT;
		//int s = 1;
		//for(long ix=0; ix<Nx; ix++)
		//{
		//	if(s < 0)
		//	{
		//		FFTW_COMPLEX *tMany = t;
		//		for(long k=0; k<HowMany; k++)
		//		{
		//			tMany->re = -tMany->re; tMany->im = -tMany->im;
		//			tMany += Nx;
		//		}
		//	}
		//	t++; s = -s;
		//}
		//OC27102018
		//SY: optimized, adopt for OpenMP
#ifdef _WITH_OMP
		#pragma omp parallel for
#endif
		for(long ix=1; ix<Nx; ix+=2)
		{
			//FFTW_COMPLEX *t = pAfterFFT + ix;
			//FFTW_COMPLEX *tMany = t;
			//OC27102018
			FFTW_COMPLEX *tMany = pAfterFFT + ix;
			for(long k=0; k<HowMany; k++)
			{
				tMany->re = -tMany->re; tMany->im = -tMany->im;
				tMany += Nx;
			}
		}
	}
#endif

#ifdef _FFTW3 //OC29012019
	template <class T> void RotateDataAfter1DFFT(T* pAfterFFT, long HowMany) //OC02022019
	//void RotateDataAfter1DFFT(fftwf_complex* pAfterFFT, long HowMany)
	{// Assumes Nx even !
#ifndef _WITH_OMP //OC27102018
		//fftwf_complex *t1 = pAfterFFT, *t2 = pAfterFFT + HalfNx;
		//fftwf_complex Buf;
		T *t1 = pAfterFFT, *t2 = pAfterFFT + HalfNx, Buf;
		for(long ix=0; ix<HalfNx; ix++)
		{
			//fftwf_complex *t1Many = t1++, *t2Many = t2++;
			T *t1Many = t1++, *t2Many = t2++;
			for(long k=0; k<HowMany; k++)
			{
				Buf[0] = (*t1Many)[0]; Buf[1] = (*t1Many)[1];
				(*t1Many)[0] = (*t2Many)[0]; (*t1Many)[1] = (*t2Many)[1];
				(*t2Many)[0] = Buf[0]; (*t2Many)[1] = Buf[1]; 
				t1Many += Nx; t2Many += Nx; 
			}
		}
#else //OC27102018
		//SY: adopted for OpenMP
		#pragma omp parallel for
		for(long ix=0; ix<HalfNx; ix++)
		{
			//fftwf_complex *t1Many = pAfterFFT + ix;
			//fftwf_complex *t2Many = pAfterFFT + HalfNx + ix;
			//fftwf_complex Buf;
			T *t1Many = pAfterFFT + ix; //OC02022019
			T *t2Many = pAfterFFT + HalfNx + ix;
			T Buf;
			for(long k=0; k<HowMany; k++)
			{
				Buf[0] = (*t1Many)[0]; Buf[1] = (*t1Many)[1]; 
				(*t1Many)[0] = (*t2Many)[0]; (*t1Many)[1] = (*t2Many)[1];
				(*t2Many)[0] = Buf[0]; (*t2Many)[1] = Buf[1];
				t1Many += Nx; t2Many += Nx; 
			}
		}
#endif
	}
#else
	void RotateDataAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany)
	{// Assumes Nx even !
#ifndef _WITH_OMP //OC27102018
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
#else //OC27102018
		//SY: adopted for OpenMP
		#pragma omp parallel for
		for(long ix=0; ix<HalfNx; ix++)
		{
			FFTW_COMPLEX *t1Many = pAfterFFT + ix;
			FFTW_COMPLEX *t2Many = pAfterFFT + HalfNx + ix;
			FFTW_COMPLEX Buf;
			for(long k=0; k<HowMany; k++)
			{
				Buf = *t1Many; *t1Many = *t2Many; *t2Many = Buf;
				t1Many += Nx; t2Many += Nx; 
			}
		}
#endif
	}
#endif

#ifdef _FFTW3 //OC29012019
	void NormalizeDataAfter1DFFT(fftwf_complex* pAfterFFT, long HowMany, double Mult)
	{// Assumes Nx even !
		float fMult = (float)Mult;
#ifndef _WITH_OMP //OC27102018
		fftwf_complex *t = pAfterFFT;
		for(long ix=0; ix<Nx; ix++)
		{
			fftwf_complex *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				(*tMany)[0] *= fMult; (*tMany)[1] *= fMult;
				tMany += Nx;
			}
		}
#else //OC27102018
		//SY: adopted for OpenMP
		#pragma omp parallel for
		for(long ix=0; ix<Nx; ix++)
		{
			fftwf_complex *tMany = pAfterFFT + ix;
			for(long k=0; k<HowMany; k++)
			{
				(*tMany)[0] *= fMult; (*tMany)[1] *= fMult;
				tMany += Nx;
			}
		}
#endif
	}
	void NormalizeDataAfter1DFFT(fftw_complex* pAfterFFT, long HowMany, double Mult)
	{// Assumes Nx even !
#ifndef _WITH_OMP //OC27102018
		fftw_complex *t = pAfterFFT;
		for(long ix=0; ix<Nx; ix++)
		{
			fftw_complex *tMany = t++;
			for(long k=0; k<HowMany; k++)
			{
				(*tMany)[0] *= Mult; (*tMany)[1] *= Mult;
				tMany += Nx;
			}
		}
#else //OC27102018
		//SY: adopted for OpenMP
		#pragma omp parallel for
		for(long ix=0; ix<Nx; ix++)
		{
			fftw_complex *tMany = pAfterFFT + ix;
			for(long k=0; k<HowMany; k++)
			{
				(*tMany)[0] *= Mult; (*tMany)[1] *= Mult;
				tMany += Nx;
			}
		}
#endif
	}

#else //ifndef _FFTW3
	void NormalizeDataAfter1DFFT(FFTW_COMPLEX* pAfterFFT, long HowMany, double Mult)
	{// Assumes Nx even !
#ifndef _WITH_OMP //OC27102018
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
#else //OC27102018
		//SY: adopted for OpenMP
		#pragma omp parallel for
		for(long ix=0; ix<Nx; ix++)
		{
			FFTW_COMPLEX *tMany = pAfterFFT + ix;
			for(long k=0; k<HowMany; k++)
			{
				tMany->re *= (FFTW_REAL)Mult; tMany->im *= (FFTW_REAL)Mult;
				tMany += Nx;
			}
		}
#endif
	}
#endif

	int SetupAuxDataForSharpEdgeCorr(CGenMathFFT1DInfo&, CGenMathAuxDataForSharpEdgeCorr1D&, char dataType='f'); //OC02022019
	//int SetupAuxDataForSharpEdgeCorr(CGenMathFFT1DInfo&, CGenMathAuxDataForSharpEdgeCorr1D&);
	void MakeSharpEdgeCorr(CGenMathFFT1DInfo&, CGenMathAuxDataForSharpEdgeCorr1D&);

	template <class T> void SetupSharpEdgeExpCorrArray(T* pCmpData, long AmOfPt, double x, double qStart, double qStep) //OC02022019
	//void SetupSharpEdgeExpCorrArray(float* pCmpData, long AmOfPt, double x, double qStart, double qStep)
	{
		const double TwoPi = 6.28318530717959;
		double TwoPiX = TwoPi*x;
		double q = qStart;
		//float *tCmpData = pCmpData;
		T *tCmpData = pCmpData;
		for(long i=0; i<AmOfPt; i++)
		{
			double Arg = TwoPiX*q;
			//float Co, Si;
			T Co, Si;
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
