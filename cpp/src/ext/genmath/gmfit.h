/************************************************************************//**
 * File: gmfit.h
 * Description: Fitting utilities (taken from "Numerical Recipes in C") header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMFIT_H
#define __GMFIT_H

#ifndef NR_END
#define NR_END 1
#endif
#ifndef FREE_ARG
#define FREE_ARG char*
#endif

//#include <cmath>
#include <math.h>

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

#ifndef MEMORY_ALLOCATION_FAILURE
#define MEMORY_ALLOCATION_FAILURE 8 + 10000 //in line with SRW
#endif
#ifndef ERROR_IN_FFT
#define INTERNAL_ERROR_AT_LINEAR_FITTING 53 + 10000 //in line with SRW
#endif

//*************************************************************************

class CGenMathFit
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{
public:

	int lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[], int ma, float **covar, float *chisq, void (CGenMathFit::*FunPtr)(float, float[], int));
	void covsrt(float **covar, int ma, int ia[], int mfit);
	int gaussj(float **a, int n, float **b, int m);

	int AllocateMatrix(long nrl, long nrh, long ncl, long nch, float**& matrix);
	int AllocateVector(long nl, long nh, float*& v);
	int AllocateVector(long nl, long nh, int*& v);

	void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
	void free_vector(float *v, long nl, long nh);
	void free_vector(int *v, long nl, long nh);

	void gser(float *gamser, float a, float x, float *gln);
	void gcf(float *gammcf, float a, float x, float *gln);

	void fpoly(float x, float p[], int np)
	{
		p[1]=1.0;
		for (int j=2;j<=np;j++) p[j]=p[j-1]*x;
	}
	int FitPolynomial(float x[], float y[], float sig[], int ndat, float a[], int ia[], int ma, float *chisq, float *q)
	{
		int result;
		float** covar;
		if(result = AllocateMatrix(1, ma, 1, ma, covar)) return result;
		if(result = lfit(x, y, sig, ndat, a, ia, ma, covar, chisq, &CGenMathFit::fpoly)) 
		{
			free_matrix(covar, 1, ma, 1, ma); //OC291009
			return result;
		}
		*q = gammq(float(0.5*(ndat - ma)), float(0.5*(*chisq)));

		free_matrix(covar, 1, ma, 1, ma); //OC291009
		return 0;
	}

	float gammq(float a, float x)
	{
		float gamser,gammcf,gln;
		if(x < 0.0 || a <= 0.0) return float(-1.E+23); // Error
		if(x < (a + 1.0)) { gser(&gamser,a,x,&gln); return float(1.0-gamser);}
		else { gcf(&gammcf,a,x,&gln); return gammcf;}
	}
	float gammln(float xx)
	{
		double x,y,tmp,ser;
		static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
		y=x=xx;
		tmp=x+5.5;
		tmp -= (x+0.5)*log(tmp);
		ser=1.000000000190015;
		for(int j=0;j<=5;j++) ser += cof[j]/++y;
		return float(-tmp+log(2.5066282746310005*ser/x));
	}
};

//*************************************************************************

#endif
