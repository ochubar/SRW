/************************************************************************//**
 * File: gmfit.cpp
 * Description: Function fitting utilities (some taken from "Numerical Recipes in C")
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "gmfit.h"
//#include "srercode.h"

#ifdef __IGOR_PRO__
#include "srigintr.h"
#endif

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

#include <stdlib.h>

//*************************************************************************

int CGenMathFit::AllocateMatrix(long nrl, long nrh, long ncl, long nch, float**& m)
// allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;

// allocate pointers to rows
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) return MEMORY_ALLOCATION_FAILURE;
	m += NR_END;
	m -= nrl;

// allocate rows and set pointers to them
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) return MEMORY_ALLOCATION_FAILURE;
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return 0;
}

//*************************************************************************

int CGenMathFit::AllocateVector(long nl, long nh, float*& v)
/* allocate a float vector with subscript range v[nl..nh] */
{
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) return MEMORY_ALLOCATION_FAILURE;
	v -= nl; v += NR_END;
	return 0;
}

//*************************************************************************

int CGenMathFit::AllocateVector(long nl, long nh, int*& v)
/* allocate an int vector with subscript range v[nl..nh] */
{
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) return MEMORY_ALLOCATION_FAILURE;
	v -= nl; v += NR_END;
	return 0;
}

//*************************************************************************

void CGenMathFit::free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

//*************************************************************************

void CGenMathFit::free_vector(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

//*************************************************************************

void CGenMathFit::free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

//*************************************************************************

int CGenMathFit::lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[], int ma, float **covar, float *chisq, void (CGenMathFit::*FunPtr)(float, float[], int))
{
	int i,j,k,l,m,mfit=0;
	float ym,wt,sum,sig2i,**beta,*afunc;

	int result;
	if(result = AllocateMatrix(1, ma, 1, 1, beta)) return result;
	if(result = AllocateVector(1, ma, afunc)) return result;

	for (j=1;j<=ma;j++) if (ia[j]) mfit++;

	for (j=1;j<=mfit;j++) 
	{
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}

	for (i=1;i<=ndat;i++) 
	{
		(this->*FunPtr)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) 
		{
			for (j=1;j<=ma;j++) if (!ia[j]) ym -= a[j]*afunc[j];
		}

		sig2i=(float)(1.0/(sig[i]*sig[i]));
		for (j=0,l=1;l<=ma;l++) 
		{
			if (ia[l]) 
			{
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++) if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) covar[k][j]=covar[j][k];
	if(result = gaussj(covar,mfit,beta,1)) return result;
	for (j=0,l=1;l<=ma;l++) if (ia[l]) a[l]=beta[++j][1];
	*chisq=0.0;
	for (i=1;i<=ndat;i++) 
	{
		(this->*FunPtr)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];

		float BufVar = (y[i]-sum)/sig[i];
		*chisq += BufVar*BufVar;
	}

	covsrt(covar,ma,ia,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
	return 0;
}

//*************************************************************************

void CGenMathFit::covsrt(float **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	float swap;
   
	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) 
	{
		if (ia[j]) 
		{
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}

//*************************************************************************

int CGenMathFit::gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,swap;

	int result;
	if(result = AllocateVector(1, n, indxc)) return result;
	if(result = AllocateVector(1, n, indxr)) return result;
	if(result = AllocateVector(1, n, ipiv)) return result;

	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) 
	{
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) 
				{
					if (ipiv[k] == 0) 
					{
						if (::fabs(a[j][k]) >= big) 
						{
							big=(float)::fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} 
					else if (ipiv[k] > 1) return INTERNAL_ERROR_AT_LINEAR_FITTING;
				}
		++(ipiv[icol]);
		if (irow != icol) 
		{
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) return INTERNAL_ERROR_AT_LINEAR_FITTING;
		pivinv=(float)(1.0/a[icol][icol]);
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) 
			{
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) 
	{
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++) SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}

	free_vector(ipiv,1,n);
	free_vector(indxr,1,n);
	free_vector(indxc,1,n);
	return 0;
}

//*************************************************************************

void CGenMathFit::gser(float *gamser, float a, float x, float *gln)
{
	const int ITMAX = 1000;
	const double EPS = 3.0e-7;

	float sum,del,ap;
	*gln=gammln(a);
   
	if(x <= 0.0) 
	{
		*gamser=0.0; return;
	} 
	else 
	{
		ap=a;
		del=sum=(float)(1.0/a);
		for(int n=1; n<=ITMAX; n++) 
		{
			++ap;
			del *= x/ap;
			sum += del;
			if(::fabs(del) < ::fabs(sum)*EPS) 
			{
				*gamser=(float)(sum*exp(-x+a*log(x)-(*gln))); return;
			}
		}
	}
}

//*************************************************************************

void CGenMathFit::gcf(float *gammcf, float a, float x, float *gln)
{
	const int ITMAX = 1000;
	const double EPS = 3.0e-7;
	const double FPMIN = 1.0e-30;

	*gln=gammln(a);
	float b=(float)(x+1.0-a);
	float c=(float)(1.0/FPMIN);
	float d=(float)(1.0/b);
	float h=d;
	float an,del;

	for(int i=1;i<=ITMAX;i++) 
	{
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if(::fabs(d) < FPMIN) d=(float)FPMIN;
		c=b+an/c;
		if(::fabs(c) < FPMIN) c=(float)FPMIN;
		d=(float)(1.0/d);
		del=d*c;
		h *= del;
		if(::fabs(del-1.0) < EPS) break;
	}
	double ExpArg = -x+a*log(x)-(*gln);
	*gammcf=(float)(exp(ExpArg)*h);
}

//*************************************************************************

#undef SWAP
