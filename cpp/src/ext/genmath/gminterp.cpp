/************************************************************************//**
 * File: gminterp.cpp
 * Description: Interpolation routines (some taken from "Numerical Recipes in C")
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "gminterp.h"

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------

CGenMathInterp::CGenMathInterp(int MethNo, double *x, double *y, int np)
{
	mSplineY2Arr = mSplineArgTabArr = mSplineValTabArr = 0;
	mSplineTabNp = 0;

	mMethNo = MethNo;
	if(MethNo == 1) // Cubic Spline
	{
		InitCubicSpline(x, y, np);
	}
}

//-------------------------------------------------------------------------
// For uniform mesh
//-------------------------------------------------------------------------

CGenMathInterp::CGenMathInterp(int MethNo, double xStart, double xStep, double *y, int np)
{
	//mSplineY2Arr = mSplineValTabArr = 0;
	mSplineY2Arr = mSplineArgTabArr = mSplineValTabArr = 0; //OC011213
	mArgStart = mArgStep = 0;
	mSplineTabNp = 0;

	mMethNo = MethNo;
	if(MethNo == 1) // Cubic Spline
	{
		InitCubicSplineU(xStart, xStep, y, np);
	}
}

//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------

void CGenMathInterp::InitCubicSpline(double *x, double *y, int np)
{
	mSplineY2Arr = new double[np];
    InterpCubicSplinePrep(x, y, np, mSplineY2Arr);

    mSplineArgTabArr = new double[np];
    mSplineValTabArr = new double[np];
	double *tx = x, *ty = y, *tSplineArgTabArr = mSplineArgTabArr, *tSplineValTabArr = mSplineValTabArr;
	for(int i=0; i<np; i++)
	{
        *(tSplineArgTabArr++) = *(tx++);
        *(tSplineValTabArr++) = *(ty++);
	}
	mSplineTabNp = np;
}

//-------------------------------------------------------------------------
// For uniform mesh spline
//-------------------------------------------------------------------------

void CGenMathInterp::InitCubicSplineU(double xStart, double xStep, double *y, int np)
{
	mArgStart = xStart;
	mArgStep = xStep;
	mSplineY2Arr = new double[np];

	InterpCubicSplinePrepU(xStart, xStep, y, np, mSplineY2Arr);
	//mSplineArgTabArr = new double[np];
	mSplineValTabArr = new double[np];

	//double *tx = x, *ty = y, *tSplineArgTabArr = mSplineArgTabArr, *tSplineValTabArr = mSplineValTabArr;
	double *ty = y, *tSplineValTabArr = mSplineValTabArr;
	for(int i=0; i<np; i++)
	{
        //*(tSplineArgTabArr++) = *(tx++);
        *(tSplineValTabArr++) = *(ty++);
	}
	mSplineTabNp = np;
	mArgStart = xStart;
	mArgStep = xStep;
}

//-------------------------------------------------------------------------
// From "Numerical Recipes in C"
//-------------------------------------------------------------------------

//void CGenMathInterp::CubSplinePrep(double *x, double *y, int n, double yp1, double ypn, double *y2)
void CGenMathInterp::InterpCubicSplinePrep(double *x, double *y, int n, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

	if(n < 2) throw "ERROR: A function tabulated at more than one point is required";

	int NpEdgeDer = 2;
    double yp1 = Deriv1(y, x[1] - x[0], 0, NpEdgeDer);
	double ypn = Deriv1(y + (n-2), x[n-1] - x[n-2], 1, NpEdgeDer);

	//u=vector(1,n-1);
	u = new double[n - 1];

	//if(yp1 > 0.99e30) y2[1]=u[1]=0.0; //The lower boundary condition is set either to be “natural”
	if(yp1 > 0.99e30) y2[0]=u[0]=0.0; //The lower boundary condition is set either to be “natural”
	else 
	{ //or else to have a specified first derivative.
		//y2[1] = -0.5;
		y2[0] = -0.5;
		//u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}

	//for(i=2; i<=n-1; i++) 
	for(i=1; i<=n-2; i++) 
	{	
		//This is the decomposition loop of the tridiagonal algorithm.
		//y2 and u are used for temporary
		//storage of the decomposed
		//factors.
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1] + 2.0;
		y2[i] = (sig - 1.0)/p;
		u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
	}

	if(ypn > 0.99e30) qn=un=0.0; //The upper boundary condition is set either to be “natural” 
	else 
	{ //or else to have a specified first derivative.
		qn = 0.5;
		//un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
		un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
	}
		
	//y2[n] = (un - qn*u[n-1])/(qn*y2[n-1] + 1.0);
	y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);

	//for(k=n-1; k>=1; k--) y2[k] = y2[k]*y2[k+1]+u[k]; //This is the backsubstitution loop of the tridiagonal algorithm. 
	for(k=n-2; k>=0; k--) y2[k] = y2[k]*y2[k+1] + u[k]; //This is the backsubstitution loop of the tridiagonal algorithm. 
		
	//free_vector(u,1,n-1);
	delete[] u;
}

//-------------------------------------------------------------------------
// From "Numerical Recipes in C", for uniform mesh
//-------------------------------------------------------------------------

void CGenMathInterp::InterpCubicSplinePrepU(double xStart, double xStep, double *y, int n, double *y2)
{
	int i,k;
	double p,qn,un;

	if(n < 2) throw "ERROR: A function tabulated at more than one point is required";

	int NpEdgeDer = 2;
    double yp1 = Deriv1(y, xStep, 0, NpEdgeDer);
	double ypn = Deriv1(y + (n-2), xStep, 1, NpEdgeDer);

	double *u = new double[n - 1];

	if(yp1 > 0.99e30) y2[0]=u[0]=0.0; //The lower boundary condition is set either to be “natural”
	else 
	{ //or else to have a specified first derivative.
		//y2[1] = -0.5;
		y2[0] = -0.5;
		//u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		u[0] = (3.0/xStep)*((y[1] - y[0])/xStep - yp1);
	}

	//for(i=2; i<=n-1; i++) 
	for(i=1; i<=n-2; i++) 
	{	
		//This is the decomposition loop of the tridiagonal algorithm.
		//y2 and u are used for temporary
		//storage of the decomposed factors.

		//sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		//p = sig*y2[i-1] + 2.0;
		p = 0.5*y2[i-1] + 2.;

		//y2[i] = (sig - 1.0)/p;
		y2[i] = -0.5/p;
		//u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
		u[i] = (y[i+1] - 2*y[i] + y[i-1])/xStep;
		//u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
		u[i] = (3.*u[i]/xStep - 0.5*u[i-1])/p;
	}

	if(ypn > 0.99e+30) qn=un=0.; //The upper boundary condition is set either to be “natural” 
	else 
	{ //or else to have a specified first derivative.
		qn = 0.5;
		//un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
		un = (3.0/xStep)*(ypn - (y[n-1] - y[n-2])/xStep);
	}

	y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.);

	for(k=n-2; k>=0; k--) y2[k] = y2[k]*y2[k+1] + u[k]; //This is the backsubstitution loop of the tridiagonal algorithm. 

	delete[] u;
}

//-------------------------------------------------------------------------
// From "Numerical Recipes in C"
//-------------------------------------------------------------------------

//void CGenMathInterp::CubSplineInterp(double *xa, double *ya, double *y2a, int n, double x, double *y)
double CGenMathInterp::InterpCubicSpline(double *xa, double *ya, double *y2a, int n, double x)
{
	//void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;
	//We will find the right place in the table by means of
	//bisection. This is optimal if sequential calls to this
	//routine are at random values of x. If sequential calls
	//are in order, and closely spaced, one would do better
	//to store previous values of klo and khi and test if
	//they remain appropriate on the next call.
	//klo = 1; 
	klo = 0; 
	//khi = n;
	khi = n - 1;
	while(khi - klo > 1) 
	{
		k = (khi + klo) >> 1;
		if(xa[k] > x) khi = k;
		else klo = k;
	} //klo and khi now bracket the input value of x.

	h = xa[khi] - xa[klo];
	//if(h == 0.0) nrerror("Bad xa input to routine splint"); //The xa’s must be distinct.
	//if(h == 0.0) { *y = *ya; return;}
	if(h == 0.0) return *ya;

	a = (xa[khi] - x)/h;
	b = (x - xa[klo])/h; //Cubic spline polynomial is now evaluated.
	//*y = a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
	return a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}

//-------------------------------------------------------------------------
// Calculates first derivative using from 2 to 5 points
//-------------------------------------------------------------------------

double CGenMathInterp::Deriv1(double* f, double h, int PoIndx, int AmOfPo)
{
	if(AmOfPo==5)
	{
		if(PoIndx==2) return 0.08333333333333*(f[0] - 8.*f[1] + 8.*f[3] - f[4])/h;
		else if(PoIndx==1) return 0.08333333333333*(-3.*f[0] - 10.*f[1] + 18.*f[2] - 6.*f[3] + f[4])/h;
		else if(PoIndx==3) return 0.08333333333333*(-f[0] + 6.*f[1] - 18.*f[2] + 10.*f[3] + 3.*f[4])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==4) return 0.5*(f[2] - 4.*f[3] + 3.*f[4])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==4)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==2) return 0.5*(-f[1] + f[3])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==3) return 0.5*(f[1] - 4.*f[2] + 3.*f[3])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==3)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==2) return 0.5*(f[0] - 4.*f[1] + 3.*f[2])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==2) return (-f[0] + f[1])/h;
	else return 0.;
}

//-------------------------------------------------------------------------

void CGenMathInterp::CompDerivForOrigData(double* OrigF, double* DerF)
{
	if((OrigF == 0) || (DerF == 0) || (OrigNp <= 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;

	double SubArray[5];

	for(int k=0; k<5; k++) 
	{
		SubArray[k] = OrigF[k];
	}
	DerF[0] = Derivative(SubArray, Orig_sStep, 0);
	DerF[1] = Derivative(SubArray, Orig_sStep, 1);
	DerF[2] = Derivative(SubArray, Orig_sStep, 2);

	int LenFieldData_m_2 = OrigNp - 2;
	for(int i=3; i<LenFieldData_m_2; i++)
	{
		int i_m_2 = i - 2;
		for(int k=0; k<5; k++) 
		{
			SubArray[k] = OrigF[i_m_2 + k];
		}
		DerF[i] = Derivative(SubArray, Orig_sStep, 2);
	}

	DerF[LenFieldData_m_2] = Derivative(SubArray, Orig_sStep, 3);
	DerF[OrigNp - 1] = Derivative(SubArray, Orig_sStep, 4);
}

//-------------------------------------------------------------------------

double CGenMathInterp::Derivative(double* f, double h, int PoIndx, int AmOfPo)
{
	if(AmOfPo==5)
	{
		if(PoIndx==2) return 0.08333333333333*(f[0] - 8.*f[1] + 8.*f[3] - f[4])/h;
		else if(PoIndx==1) return 0.08333333333333*(-3.*f[0] - 10.*f[1] + 18.*f[2] - 6.*f[3] + f[4])/h;
		else if(PoIndx==3) return 0.08333333333333*(-f[0] + 6.*f[1] - 18.*f[2] + 10.*f[3] + 3.*f[4])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==4) return 0.5*(f[2] - 4.*f[3] + 3.*f[4])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==4)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==2) return 0.5*(-f[1] + f[3])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==3) return 0.5*(f[1] - 4.*f[2] + 3.*f[3])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==3)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==2) return 0.5*(f[0] - 4.*f[1] + 3.*f[2])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==2) return (-f[0] + f[1])/h;
	else return 1.E+23;
}

//-------------------------------------------------------------------------

void CGenMathInterp::Interpolate(double sSt, double sStp, int Np, double* pInterpData)
{
	double s = sSt;
	for(long i=0; i<Np; i++)
	{
		int Indx = int((s - Orig_sStart)/Orig_sStep); 
		if(Indx >= OrigNp - 1) Indx = OrigNp - 2;

		double sb = Orig_sStart + Indx*Orig_sStep;
		double smsb = s - sb;
		double *B_CfP = PlnCf[Indx];
        *(pInterpData++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));

		s += sStp;
	}
}

//-------------------------------------------------------------------------
