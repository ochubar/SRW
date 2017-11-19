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
#include "auxparse.h"
#include <algorithm>
#include <math.h>

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

//bool CGenMathInterp::TryToFindRectMesh(double* arX, double* arY, int nPtTot, double* arResMeshX, double* arResMeshY, int& nMeshX, int& nMeshY, int* arPtMeshInd, double arMinMaxX[2], double arMinMaxY[2])
bool CGenMathInterp::TryToFindRectMesh(double* arX, double* arY, int nPtTot, double* arResMeshX, double* arResMeshY, int& nMeshX, int& nMeshY, double arMinMaxX[2], double arMinMaxY[2])
{
	if((arX == 0) || (arY == 0) || (nPtTot < 2)) return false;

	//Find different grid points vs X and Y
	double relTolSameGridPoint = 1.e-10; //to tune
	double xMin=0, xMax=0;
	long ixMin=-1, ixMax=-1;
	CAuxParse::FindMinMax(arX, (long)nPtTot, xMin, ixMin, xMax, ixMax);
	double absTolSameGridPointX = relTolSameGridPoint*(xMax - xMin);
	double yMin=0, yMax=0;
	long iyMin=-1, iyMax=-1;
	CAuxParse::FindMinMax(arY, (long)nPtTot, yMin, iyMin, yMax, iyMax);
	double absTolSameGridPointY = relTolSameGridPoint*(yMax - yMin);

	vector<double> vMeshX, vMeshY;
	for(int i=0; i<nPtTot; i++)
	{
		double curX = arX[i];
		if(CAuxParse::FindIndVectElemEqualToWithinTol(curX, vMeshX, absTolSameGridPointX) < 0) vMeshX.push_back(curX);
		double curY = arY[i];
		if(CAuxParse::FindIndVectElemEqualToWithinTol(curY, vMeshY, absTolSameGridPointY) < 0) vMeshY.push_back(curY);
	}
	sort(vMeshX.begin(), vMeshX.end());
	sort(vMeshY.begin(), vMeshY.end());

	int npX = (int)vMeshX.size();
	int npY = (int)vMeshY.size();

	arMinMaxX[0] = vMeshX[0]; arMinMaxX[1] = vMeshX[npX - 1];
	arMinMaxY[0] = vMeshY[0]; arMinMaxY[1] = vMeshY[npY - 1];

	if((npX < 2) && (npY >= 2))
	{
		arResMeshX[0] = vMeshX[0];
		nMeshX = 1;
		for(int iy=0; iy<npY; iy++) arResMeshY[iy] = vMeshY[iy];
		nMeshY = npY;
		return true;
	}
	else if((npY < 2) && (npX >= 2))
	{
		arResMeshY[0] = vMeshY[0];
		nMeshY = 1;
		for(int ix=0; ix<npX; ix++) arResMeshX[ix] = vMeshX[ix];
		nMeshX = npX;
		return true;
	}

	//For every "point" of grid vx Y, make sure that the grid points over X are the same
	const double maxPartExceptFromRectMesh = 0.3; //0.2;
	//const double badArgVal = 1.e+23;

	//double dMaxNumOfExceptFromRectMeshX = npX*maxPartExceptFromRectMesh;
	//int maxNumOfExceptFromRectMeshX = int(dMaxNumOfExceptFromRectMeshX);
	//if((dMaxNumOfExceptFromRectMeshX - maxNumOfExceptFromRectMeshX) >= 0.5) maxNumOfExceptFromRectMeshX++;
	//
	//double dMaxNumOfExceptFromRectMeshY = npY*maxPartExceptFromRectMesh;
	//int maxNumOfExceptFromRectMeshY = int(dMaxNumOfExceptFromRectMeshY);
	//if((dMaxNumOfExceptFromRectMeshY - maxNumOfExceptFromRectMeshY) >= 0.5) maxNumOfExceptFromRectMeshY++;

	int *arNumBadX = new int[npY];
	int *arNumBadY = new int[npX];
	for (int iix = 0; iix < npX; iix++) arNumBadY[iix] = 0;
	
	//int exceptCountY = 0;
	//int *t_arPtMeshInd = arPtMeshInd;
	for(int iy=0; iy<npY; iy++)
	{
		//int exceptCountX = 0;
		arNumBadX[iy] = 0;

		double curY = vMeshY[iy];
		for(int ix=0; ix<npX; ix++)
		{
			double curX = vMeshX[ix];
			bool thisPointExists = false;
			//*t_arPtMeshInd = -1;
			for(int i=0; i<nPtTot; i++)
			{
				double curDifX = arX[i] - curX;
				if((-absTolSameGridPointX <= curDifX) && (curDifX <= absTolSameGridPointX))
				{
					double curDifY = arY[i] - curY;
					if((-absTolSameGridPointY <= curDifY) && (curDifY <= absTolSameGridPointY))
					{
						thisPointExists = true;
						//*t_arPtMeshInd = i;
						break;
					}
				}
			}
			if(!thisPointExists)
			{
				arNumBadY[ix] += 1;
				arNumBadX[iy] += 1;
				//return false;
				//exceptCountX++;
				//if(exceptCountX > maxNumOfExceptFromRectMeshX) return false;
				//else
				//{
				//	arResMeshX[ix] = badArgVal;
				//	break;
				//}
			}
			if(iy == 0) arResMeshX[ix] = curX;
			//t_arPtMeshInd++;
		}
		arResMeshY[iy] = curY;
	}
	nMeshX = npX;
	nMeshY = npY;

	//Analyze "bad" (/missing) nodes
	int maxNumBadNodesX = 0, maxNumBadNodesY = 0, i, curNumBad;
	int numBadRowsX = 0, numBadRowsY = 0;
	int *pArNumBad = arNumBadX;
	int *pMaxNumBadNodes = &maxNumBadNodesX;
	int *pNumBadRows = &numBadRowsX;
	int *pNp = &npY;
	for(i=0; i<(*pNp); i++)
	{
		curNumBad = pArNumBad[i];
		if((*pMaxNumBadNodes) < curNumBad) *pMaxNumBadNodes = curNumBad;
		if(curNumBad > 0) (*pNumBadRows)++;
	}

	pArNumBad = arNumBadY;
	pMaxNumBadNodes = &maxNumBadNodesY;
	pNumBadRows = &numBadRowsY;
	pNp = &npX;
	for(i=0; i<(*pNp); i++)
	{
		curNumBad = pArNumBad[i];
		if((*pMaxNumBadNodes) < curNumBad) *pMaxNumBadNodes = curNumBad;
		if(curNumBad > 0) (*pNumBadRows)++;
	}

	int maxNumBadNodes = maxNumBadNodesX;
	int numBadRows = numBadRowsX;
	int totNumRows = npY; //npX;
	pArNumBad = arNumBadX;
	double *pArResMesh = arResMeshY; //arResMeshX;
	int *pnMesh = &nMeshY; //&nMeshX;
	if(maxNumBadNodesX < maxNumBadNodesY)
	{
		maxNumBadNodes = maxNumBadNodesY;
		numBadRows = numBadRowsY;

		totNumRows = npX; //npY;
		pArNumBad = arNumBadY;
		pArResMesh = arResMeshX; //arResMeshY;
		pnMesh = &nMeshX; //&nMeshY;
	}

	if(maxNumBadNodes <= 0)
	{
		delete[] arNumBadX; delete[] arNumBadY;
		return true;
	}

	if(numBadRows > maxPartExceptFromRectMesh*totNumRows)
	{
		delete[] arNumBadX; delete[] arNumBadY;
		return false;
	}

	vector<double> vAuxResMesh;
	for(i=0; i<totNumRows; i++)
	{
		if(pArNumBad[i] == 0) vAuxResMesh.push_back(pArResMesh[i]);
	}
	*pnMesh = (int)(vAuxResMesh.size());

	for(i=0; i<(*pnMesh); i++) pArResMesh[i] = vAuxResMesh[i];

	delete[] arNumBadX; delete[] arNumBadY;
	return true;
}

//-------------------------------------------------------------------------

//void CGenMathInterp::SelectPointsForInterp2d(double x, double y, double* arX, double* arY, int nPtTot, int ord, int* arResInd, int& nResInd)
int CGenMathInterp::SelectPointsForInterp2d(double x, double y, double* arX, double* arY, int nPtTot, int& ord, int* arResInd, int& nResInd, bool& resMeshIsRect)
{//Selects points for 2D interpolation on irregular mesh
 //Returns effective number of dimensions suggested for the interpolation
	nResInd = 0;
	if((arX == 0) || (arY == 0)) return 0; //throw exception?

	//if((ord < 1) && (ord > 2)) return; //throw exception?
	if(ord < 1) return 0; //throw exception?
	if(ord > 2) ord = 2; //to update when/if necessary

	double *arMeshX = new double[nPtTot];
	double *arMeshY = new double[nPtTot];
	int nMeshX=0, nMeshY=0;
	//int *arPtMeshInd = new int[nPtTot*nPtTot]; //excessively large allocation?
	double arMinMaxX[2], arMinMaxY[2];
	//bool meshIsRect = TryToFindRectMesh(arX, arY, nPtTot, arMeshX, arMeshY, nMeshX, nMeshY, arPtMeshInd, arMinMaxX, arMinMaxY);
	resMeshIsRect = TryToFindRectMesh(arX, arY, nPtTot, arMeshX, arMeshY, nMeshX, nMeshY, arMinMaxX, arMinMaxY);

	int nMeshXmi1 = nMeshX - 1, nMeshYmi1 = nMeshY - 1;

	if(resMeshIsRect)
	{
		if((x < arMeshX[0]) || (x > arMeshX[nMeshXmi1]) || (y < arMeshY[0]) || (y > arMeshY[nMeshYmi1]))
		{
			if((arMinMaxX[0] <= x) && (x <= arMinMaxX[1]) && (arMinMaxY[0] <= y) && (y <= arMinMaxY[1]))
			{
				resMeshIsRect = false; //to try calculation on non-rect. mesh
			}
			else
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
				delete[] arMeshX; delete[] arMeshY; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
			}
		}
	}

	if(resMeshIsRect)
	{
		if((nMeshX < 1) || (nMeshY < 1) || ((nMeshY < 2) && (y != arMeshY[0])) || ((nMeshX < 2) && (x != arMeshX[0]))) 
		{
			//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
			delete[] arMeshX; delete[] arMeshY; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
		}

		if((nMeshY < 2) && (y == arMeshY[0]) && (nMeshX < 2) && (x == arMeshX[0]))
		{
			//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return 0; //point (x,y) is exactly on a single mesh point
			delete[] arMeshX; delete[] arMeshY; return 0; //point (x,y) is exactly on a single mesh point
		}

		//Treating effective 1D mesh here because formally it is a "rectangular" mesh
		double *arMesh1D = 0, arg1D; 
		int nMesh1D = 0;
		if((nMeshY < 2) && (y == arMeshY[0]) && (nMeshX >= 2))
		{//Treat as 1D interpolation vs X
			arMesh1D = arMeshX; nMesh1D = nMeshX; arg1D = x;
		}
		else if((nMeshX < 2) && (x == arMeshX[0]) && (nMeshY >= 2))
		{//Treat as 1D interpolation vs Y
			arMesh1D = arMeshY; nMesh1D = nMeshY; arg1D = y;
		}
		if((arMesh1D != 0) && (nMesh1D > 0))
		{
			if((arg1D < arMesh1D[0]) || (arg1D > arMesh1D[nMesh1D - 1]))
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
				delete[] arMeshX; delete[] arMeshY; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
			}

			int i0_1D=-1;
			for(int i=1; i<nMesh1D; i++)
			{
				if((arMesh1D[i-1] <= arg1D) && (arg1D < arMesh1D[i]))
				{
					i0_1D = i; break;
				}
			}
			if(i0_1D < 0)
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
				delete[] arMeshX; delete[] arMeshY; return -1; //point (x,y) is beyond the mesh ranges, or no mesh was found return error
			}

			int nMesh1D_mi_1 = nMesh1D - 1, nMesh1D_mi_2 = nMesh1D - 2, nMesh1D_mi_3 = nMesh1D - 3;

			if(nMesh1D == 2)
			{
				if(ord > 1) ord = 1;
			}
			else if(nMesh1D == 3)
			{
				if(ord > 2) ord = 2;
			}

			nResInd = -1;
			arResInd[0] = -1; arResInd[1] = -1; arResInd[2] = -1; arResInd[3] = -1;

			if(ord == 1)
			{
				if(i0_1D >= nMesh1D_mi_1) i0_1D = nMesh1D_mi_2;

				nResInd = 2;
				if(arMesh1D == arMeshX)
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D], y, arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D + 1], y, arX, arY, nPtTot);
				}
				else
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D], arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D + 1], arX, arY, nPtTot);
				}
			}
			else if(ord == 2)
			{
				if(i0_1D <= 0) i0_1D = 1;
				else if(i0_1D >= nMesh1D_mi_1) i0_1D = nMesh1D_mi_2;
			
				nResInd = 3;
				if(arMesh1D == arMeshX)
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D - 1], y, arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D], y, arX, arY, nPtTot);
					arResInd[2] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D + 1], y, arX, arY, nPtTot);
				}
				else
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D - 1], arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D], arX, arY, nPtTot);
					arResInd[2] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D + 1], arX, arY, nPtTot);
				}
			}
			else if(ord == 3)
			{
				if(i0_1D <= 0) i0_1D = 1;
				else if(i0_1D >= nMesh1D_mi_2) i0_1D = nMesh1D_mi_3;
			
				nResInd = 4;
				if(arMesh1D == arMeshX)
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D - 1], y, arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D], y, arX, arY, nPtTot);
					arResInd[2] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D + 1], y, arX, arY, nPtTot);
					arResInd[3] = FindIndOfPointOnIrregMesh2d(arMesh1D[i0_1D + 2], y, arX, arY, nPtTot);
				}
				else
				{
					arResInd[0] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D - 1], arX, arY, nPtTot);
					arResInd[1] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D], arX, arY, nPtTot);
					arResInd[2] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D + 1], arX, arY, nPtTot);
					arResInd[3] = FindIndOfPointOnIrregMesh2d(x, arMesh1D[i0_1D + 2], arX, arY, nPtTot);
				}
			}

			bool meshDataIsConsitent = true;
			if(nResInd <= 0) meshDataIsConsitent = false;
			else
			{
				for(int iii=0; iii<nResInd; iii++)
				{
					if(arResInd[iii] < 0)
					{
						meshDataIsConsitent = false; break;
					}
				}
			}
			
			delete[] arMeshX; delete[] arMeshY; //delete[] arPtMeshInd;
			
			if(meshDataIsConsitent) return -1;
			else return 1;
		}//End of treating effective 1D mesh
		
		//Start treating true 2D case
		if((nMeshX < 3) || (nMeshY < 3)) ord = 1; //Reduce polynomial interpolation order if not enaugh points
		else if((nMeshX < 4) || (nMeshY < 4)) 
		{
			if(ord == 3) ord = 2;
		}

		//Finding "central" mesh point
		int ixm1=-1, ix0=0, ixp1=-1, ixp2=-1;
		if(nMeshX > 1) ix0 = CAuxParse::FindLargestElemSmallerThanOrEqualTo(x, arMeshX, nMeshX);
		if((ix0 < 0) || (ix0 > nMeshXmi1)) 
		{
			//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1;
			delete[] arMeshX; delete[] arMeshY; return -1;
		}
		if((ix0 == nMeshXmi1) && (nMeshX > 1)) ix0--;

		int iym1=-1, iy0=0, iyp1=-1, iyp2=-1;
		if(nMeshY > 1) iy0 = CAuxParse::FindLargestElemSmallerThanOrEqualTo(y, arMeshY, nMeshY);
		if((iy0 < 0) || (iy0 > nMeshYmi1)) 
		{
			//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1;
			delete[] arMeshX; delete[] arMeshY; return -1;
		}
		if((iy0 == nMeshYmi1) && (nMeshY > 1)) iy0--;

		ixp1 = ix0 + 1;
		iyp1 = iy0 + 1;

		if(ord == 1)
		{
			for(int ii=0; ii<4; ii++) arResInd[ii] = -1; 

			double x0 = arMeshX[ix0], x1 = arMeshX[ixp1];
			double y0 = arMeshX[iy0], y1 = arMeshX[iyp1];

			nResInd = 4;
			arResInd[0] = FindIndOfPointOnIrregMesh2d(x0, y0, arX, arY, nPtTot);
			arResInd[1] = FindIndOfPointOnIrregMesh2d(x1, y0, arX, arY, nPtTot);
			arResInd[2] = FindIndOfPointOnIrregMesh2d(x0, y1, arX, arY, nPtTot);
			arResInd[3] = FindIndOfPointOnIrregMesh2d(x1, y1, arX, arY, nPtTot);

			//long iy0_nMeshX = iy0*nMeshX;
			//long iyp1_nMeshX = iyp1*nMeshX;
			//arResInd[0] = arPtMeshInd[iy0_nMeshX + ix0];
			//arResInd[1] = arPtMeshInd[iy0_nMeshX + ixp1];
			//arResInd[2] = arPtMeshInd[iyp1_nMeshX + ix0];
			//arResInd[3] = arPtMeshInd[iyp1_nMeshX + ixp1];
		}
		else if(ord == 2) //More tests for the case of 5-point "cross" interpolation
		{
			for(int ii=0; ii<5; ii++) arResInd[ii] = -1; 

			if(ix0 <= 0)
			{
				ixm1 = 0; ix0 = 1; ixp1 = 2;
			}
			else
			{
				ixm1 = ix0 - 1;
				double abs_x_mi_xm1 = fabs(x - arMeshX[ixm1]);
				double abs_x_mi_x0 = fabs(x - arMeshX[ix0]);
				double abs_x_mi_xp1 = fabs(x - arMeshX[ixp1]);
				if((abs_x_mi_xm1 < abs_x_mi_x0) && (ixm1 > 0)) { ixm1--; ix0--; ixp1--;}
				else if((abs_x_mi_xp1 < abs_x_mi_x0) && (ixp1 < nMeshXmi1)) { ixm1++; ix0++; ixp1++;}
			}

			if(iy0 <= 0)
			{
				iym1 = 0; iy0 = 1; iyp1 = 2;
			}
			else
			{
				iym1 = iy0 - 1;
				double abs_y_mi_ym1 = fabs(y - arMeshY[iym1]);
				double abs_y_mi_y0 = fabs(y - arMeshY[iy0]);
				double abs_y_mi_yp1 = fabs(y - arMeshY[iyp1]);
				if((abs_y_mi_ym1 < abs_y_mi_y0) && (iym1 > 0)) { iym1--; iy0--; iyp1--;}
				else if((abs_y_mi_yp1 < abs_y_mi_y0) && (iyp1 < nMeshYmi1)) { iym1++; iy0++; iyp1++;}
			}

			double xm1 = arMeshX[ixm1], x0 = arMeshX[ix0], x1 = arMeshX[ixp1];
			double ym1 = arMeshY[iym1], y0 = arMeshY[iy0], y1 = arMeshY[iyp1];

			nResInd = 5;
			arResInd[0] = FindIndOfPointOnIrregMesh2d(x0, ym1, arX, arY, nPtTot);
			arResInd[1] = FindIndOfPointOnIrregMesh2d(xm1, y0, arX, arY, nPtTot);
			arResInd[2] = FindIndOfPointOnIrregMesh2d(x0, y0, arX, arY, nPtTot);
			arResInd[3] = FindIndOfPointOnIrregMesh2d(x1, y0, arX, arY, nPtTot);
			arResInd[4] = FindIndOfPointOnIrregMesh2d(x0, y1, arX, arY, nPtTot);

			//long iy0_nMeshX = iy0*nMeshX;
			//arResInd[0] = arPtMeshInd[iym1*nMeshX + ix0];
			//arResInd[1] = arPtMeshInd[iy0_nMeshX + ixm1];
			//arResInd[2] = arPtMeshInd[iy0_nMeshX + ix0];
			//arResInd[3] = arPtMeshInd[iy0_nMeshX + ixp1];
			//arResInd[4] = arPtMeshInd[iyp1*nMeshX + ix0];
		}
		else if(ord == 3) //To test
		{
			for(int ii=0; ii<12; ii++) arResInd[ii] = -1; 

			int nMeshXmi3 = nMeshX - 3;
			if(ix0 < 1) ix0 = 1;
			else if(ix0 > nMeshXmi3) ix0 = nMeshXmi3;

			ixm1 = ix0 - 1; ixp1 = ix0 + 1; ixp2 = ix0 + 2;

			int nMeshYmi3 = nMeshY - 3;
			if(iy0 < 1) iy0 = 1;
			else if(iy0 > nMeshYmi3) iy0 = nMeshYmi3;

			iym1 = iy0 - 1; iyp1 = iy0 + 1; iyp2 = iy0 + 2;
		
			double xm1 = arMeshX[ixm1], x0 = arMeshX[ix0], x1 = arMeshX[ixp1], x2 = arMeshX[ixp2];
			double ym1 = arMeshY[iym1], y0 = arMeshY[iy0], y1 = arMeshY[iyp1], y2 = arMeshY[iyp2];

			nResInd = 12;
			arResInd[0] = FindIndOfPointOnIrregMesh2d(x1, ym1, arX, arY, nPtTot);
			arResInd[1] = FindIndOfPointOnIrregMesh2d(x0, ym1, arX, arY, nPtTot);
			arResInd[2] = FindIndOfPointOnIrregMesh2d(xm1, y0, arX, arY, nPtTot);
			arResInd[3] = FindIndOfPointOnIrregMesh2d(x0, y0, arX, arY, nPtTot);
			arResInd[4] = FindIndOfPointOnIrregMesh2d(x1, y0, arX, arY, nPtTot);
			arResInd[5] = FindIndOfPointOnIrregMesh2d(x2, y0, arX, arY, nPtTot);
			arResInd[6] = FindIndOfPointOnIrregMesh2d(xm1, y1, arX, arY, nPtTot);
			arResInd[7] = FindIndOfPointOnIrregMesh2d(x0, y1, arX, arY, nPtTot);
			arResInd[8] = FindIndOfPointOnIrregMesh2d(x1, y1, arX, arY, nPtTot);
			arResInd[9] = FindIndOfPointOnIrregMesh2d(x2, y1, arX, arY, nPtTot);
			arResInd[10] = FindIndOfPointOnIrregMesh2d(x0, y2, arX, arY, nPtTot);
			arResInd[11] = FindIndOfPointOnIrregMesh2d(x1, y2, arX, arY, nPtTot);
			//f0m1=*(p++); f1m1=*(p++); fm10=*(p++); f00=*(p++); f10=*(p++); f20=*(p++); fm11=*(p++); f01=*(p++); f11=*(p++); f21=*(p++); f02=*(p++); f12=*(p++); 
		}
	}
	else
	{//General case of a non-rectangular mesh...

		if(ord > 2) ord = 2; //to update when/if necessary
		
		//Find closest poont of mesh based on distances from a given point to all available points
		vector<pair<int, double> > vDist;
		double dx, dy;
		double xMin=arX[0], xMax=arX[0], yMin=arY[0], yMax=arY[0];
		for(int i=0; i<nPtTot; i++)
		{
			double curX = arX[i], curY = arY[i];

			if(xMin > curX) xMin = curX;
			else if(xMax < curX) xMax = curX;

			if(yMin > curY) yMin = curY;
			else if(yMax < curY) yMax = curY;

			dx = x - curX; dy = y - curY;
			vDist.push_back(pair<int, double>(i, sqrt(dx*dx + dy*dy)));
		}
		//Sort all points based on distance from (x,y), lower distances first
		sort(vDist.begin(), vDist.end(), CAuxParse::LessInPairBasedOnSecond<int, double>);

		int im11=-1, i01=-1, i11=-1; //To be used at (ord == 1)
		int im10=-1, i00=-1, i10=-1;
		int im1m1=-1, i0m1=-1, i1m1=-1;

		double abs_dxm11=1.e+23, abs_dx01=1.e+23, abs_dx11=1.e+23; //To be used at (ord == 1)
		double abs_dxm10=1.e+23, abs_dx00=1.e+23, abs_dx10=1.e+23;
		double abs_dxm1m1=1.e+23, abs_dx0m1=1.e+23, abs_dx1m1=1.e+23;

		double abs_dym11=1.e+23, abs_dy01=1.e+23, abs_dy11=1.e+23; //To be used at (ord == 1)
		double abs_dym10=1.e+23, abs_dy00=1.e+23, abs_dy10=1.e+23;
		double abs_dym1m1=1.e+23, abs_dy0m1=1.e+23, abs_dy1m1=1.e+23;

		int i0 = vDist[0].first; //To be used at (ord == 2)

		double relTol = 1.e-06; //to tune
		double absTolX = (xMax - xMin)*relTol;
		double absTolY = (yMax - yMin)*relTol;

		//Find candidates for other mesh points
		int ixm1=-1, ixm2=-1, ixp1=-1, ixp2=-1;
		int iym1=-1, iym2=-1, iyp1=-1, iyp2=-1;
		double abs_dxm1=1.e+23, abs_dxm2=1.e+23, abs_dxp1=1.e+23, abs_dxp2=1.e+23;
		double abs_dym1=1.e+23, abs_dym2=1.e+23, abs_dyp1=1.e+23, abs_dyp2=1.e+23;

		int iStart = 0;
		if(ord == 2) iStart = 1;

		for(int i=iStart; i<nPtTot; i++)
		{
			int ii = vDist[i].first;
			dx = x - arX[ii]; dy = y - arY[ii];

			double abs_rx = fabs(dx/absTolX);
			double abs_ry = fabs(dy/absTolY);

			double abs_dx = fabs(dx);
			double abs_dy = fabs(dy);

			if(ord == 1)
			{
				if(dx < 0)
				{
					if(dy < 0)
					{
						if(im1m1 < 0) { im1m1 = ii; abs_dxm1m1 = abs_dx; abs_dym1m1 = abs_dy;} 
					}
					else //if(dy >= 0)
					{
						if(im10 < 0) { im10 = ii; abs_dxm10 = abs_dx; abs_dym10 = abs_dy;} 
						else if((im11 < 0) && (abs_dxm10 < abs_dx)) { im11 = ii; abs_dxm11 = abs_dx; abs_dym11 = abs_dy;}
					}
				}
				else //if(dx >= 0)
				{
					if(dy < 0)
					{
						if(i0m1 < 0) { i0m1 = ii; abs_dx0m1 = abs_dx; abs_dy0m1 = abs_dy;}
						else if((i1m1 < 0) && (abs_dx0m1 < abs_dx)) { i1m1 = ii; abs_dx1m1 = abs_dx; abs_dy1m1 = abs_dy;}
					}
					else //if(dy >= 0)
					{
						if(i00 < 0) { i00 = ii; abs_dx00 = abs_dx; abs_dy00 = abs_dy;}
						else if((i10 < 0) && (abs_dx00 < abs_dx)) { i10 = ii; abs_dx10 = abs_dx; abs_dy10 = abs_dy;}
						else if((i01 < 0) && (abs_dy00 < abs_dy)) { i01 = ii; abs_dx01 = abs_dx; abs_dy01 = abs_dy;}
						else if((i11 < 0) && (abs_dx10 < abs_dx) && (abs_dy01 < abs_dy)) { i11 = ii; abs_dx11 = abs_dx; abs_dy11 = abs_dy;}
					}
				}
			}
			else if(ord == 2)
			{
				if((dx < -absTolX) && ((ixm1 < 0) || (ixm2 < 0)))
				{
					if((dy < -absTolY) && ((iym1 < 0) || (iym2 < 0)))
					{
						if(abs_rx >= abs_ry)
						{
							if((ixm1 < 0) && (abs_dx < abs_dxm1)) { ixm1 = ii; abs_dxm1 = abs_dx; }
							else if((ixm2 < 0) && (abs_dxm1 < abs_dx) && (abs_dx < abs_dxm2)) { ixm2 = ii; abs_dxm2 = abs_dx; }
							else if((iym1 < 0) && (abs_dy < abs_dym1)) { iym1 = ii; abs_dym1 = abs_dy; }
							else if((iym2 < 0) && (abs_dym1 < abs_dy) && (abs_dy < abs_dym2)) { iym2 = ii; abs_dym2 = abs_dy; }
						}
						else
						{
							if((iym1 < 0) && (abs_dy < abs_dym1)) { iym1 = ii; abs_dym1 = abs_dy; }
							else if((iym2 < 0) && (abs_dym1 < abs_dy) && (abs_dy < abs_dym2)) { iym2 = ii; abs_dym2 = abs_dy; }
							else if((ixm1 < 0) && (abs_dx < abs_dxm1)) { ixm1 = ii; abs_dxm1 = abs_dx; }
							else if((ixm2 < 0) && (abs_dxm1 < abs_dx) && (abs_dx < abs_dxm2)) { ixm2 = ii; abs_dxm2 = abs_dx; }
						}
					}
					else if((dy > absTolY) && ((iyp1 < 0) || (iyp2 < 0)))
					{
						if(abs_rx >= abs_ry)
						{
							if((ixm1 < 0) && (abs_dx < abs_dxm1)) { ixm1 = ii; abs_dxm1 = abs_dx; }
							else if((ixm2 < 0) && (abs_dxm1 < abs_dx) && (abs_dx < abs_dxm2)) { ixm2 = ii; abs_dxm2 = abs_dx; }
							else if((iyp1 < 0) && (abs_dy < abs_dyp1)) { iyp1 = ii; abs_dyp1 = abs_dy; }
							else if((iyp2 < 0) && (abs_dxp1 < abs_dy) && (abs_dy < abs_dyp2)) { iyp2 = ii; abs_dyp2 = abs_dy; }
						}
						else
						{
							if((iyp1 < 0) && (abs_dy < abs_dyp1)) { iyp1 = ii; abs_dyp1 = abs_dy; }
							else if((iyp2 < 0) && (abs_dyp1 < abs_dy) && (abs_dy < abs_dyp2)) { iyp2 = ii; abs_dyp2 = abs_dy; }
							else if((ixm1 < 0) && (abs_dx < abs_dxm1)) { ixm1 = ii; abs_dxm1 = abs_dx; }
							else if((ixm2 < 0) && (abs_dxm1 < abs_dx) && (abs_dx < abs_dxm2)) { ixm2 = ii; abs_dxm2 = abs_dx; }
						}
					}
				}
				else if((dx > absTolX) && ((ixp1 < 0) || (ixp2 < 0)))
				{
					if((dy < -absTolY) && ((iym1 < 0) || (iym2 < 0)))
					{
						if(abs_rx >= abs_ry)
						{
							if((ixp1 < 0) && (abs_dx < abs_dxp1)) { ixp1 = ii; abs_dxp1 = abs_dx; }
							else if((ixp2 < 0) && (abs_dxp1 < abs_dx) && (abs_dx < abs_dxp2)) { ixp2 = ii; abs_dxp2 = abs_dx; }
							else if((iym1 < 0) && (abs_dy < abs_dym1)) { iym1 = ii; abs_dym1 = abs_dy; }
							else if((iym2 < 0) && (abs_dym1 < abs_dy) && (abs_dy < abs_dym2)) { iym2 = ii; abs_dym2 = abs_dy; }
						}
						else
						{
							if((iym1 < 0) && (abs_dy < abs_dym1)) { iym1 = ii; abs_dym1 = abs_dy; }
							else if((iym2 < 0) && (abs_dym1 < abs_dy) && (abs_dy < abs_dym2)) { iym2 = ii; abs_dym2 = abs_dy; }
							else if((ixp1 < 0) && (abs_dx < abs_dxp1)) { ixp1 = ii; abs_dxp1 = abs_dx; }
							else if((ixp2 < 0) && (abs_dxp1 < abs_dx) && (abs_dx < abs_dxp2)) { ixp2 = ii; abs_dxp2 = abs_dx; }
						}
					}
					else if((dy > absTolY) && ((iyp1 < 0) || (iyp2 < 0)))
					{
						if(abs_rx >= abs_ry)
						{
							if((ixp1 < 0) && (abs_dx < abs_dxp1)) { ixp1 = ii; abs_dxp1 = abs_dx; }
							else if((ixp2 < 0) && (abs_dxp1 < abs_dx) && (abs_dx < abs_dxp2)) { ixp2 = ii; abs_dxp2 = abs_dx; }
							else if((iyp1 < 0) && (abs_dy < abs_dyp1)) { iyp1 = ii; abs_dyp1 = abs_dy; }
							else if((iyp2 < 0) && (abs_dyp1 < abs_dy) && (abs_dy < abs_dyp2)) { iyp2 = ii; abs_dyp2 = abs_dy; }
						}
						else
						{
							if((iyp1 < 0) && (abs_dy < abs_dyp1)) { iyp1 = ii; abs_dyp1 = abs_dy; }
							else if((iyp2 < 0) && (abs_dyp1 < abs_dy) && (abs_dy < abs_dyp2)) { iyp2 = ii; abs_dyp2 = abs_dy; }
							else if((ixp1 < 0) && (abs_dx < abs_dxp1)) { ixp1 = ii; abs_dxp1 = abs_dx; }
							else if((ixp2 < 0) && (abs_dxp1 < abs_dx) && (abs_dx < abs_dxp2)) { ixp2 = ii; abs_dxp2 = abs_dx; }
						}
					}
				}
				if((ixm1 >= 0) && (ixm2 >= 0) && (iym1 >= 0) && (iym2 >= 0)) break;
			}
		}

		if(ord == 1)
		{
			for(int ii=0; ii<4; ii++) arResInd[ii] = -1; 

			if((i00 < 0) && (i10 < 0) && (i01 < 0) && (i11 < 0))
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1;
				delete[] arMeshX; delete[] arMeshY; return -1;
			}
			nResInd = 4;
			//double f00 = *(arF++);
			//double f10 = *(arF++);
			//double f01 = *(arF++);
			//double f11 = *(arF++);

			if((i00 >= 0) && (i10 >= 0) && (i01 >= 0) && (i11 >= 0))
			{
				arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
			}
			else if(i11 < 0)
			{
				if((im10 > 0) && (i00 > 0) && (im11 > 0) && (i01 > 0))
				{
					arResInd[0] = im10; arResInd[1] = i00; arResInd[2] = im11; arResInd[3] = i01;
				}
				else if((i0m1 > 0) && (i1m1 > 0) && (i00 > 0) && (i10 > 0))
				{
					arResInd[0] = i0m1; arResInd[1] = i1m1; arResInd[2] = i00; arResInd[3] = i10;
				}
				else if((im1m1 > 0) && (i0m1 > 0) && (im10 > 0) && (i00 > 0))
				{
					arResInd[0] = im1m1; arResInd[1] = i0m1; arResInd[2] = im10; arResInd[3] = i00;
				}
			}
			else if(i10 < 0)
			{
				if((im10 > 0) && (i00 > 0) && (im11 > 0) && (i01 > 0))
				{
					arResInd[0] = im10; arResInd[1] = i00; arResInd[2] = im11; arResInd[3] = i01;
				}
				else if((im1m1 > 0) && (i0m1 > 0) && (im10 > 0) && (i00 > 0))
				{
					arResInd[0] = im1m1; arResInd[1] = i0m1; arResInd[2] = im10; arResInd[3] = i00;
				}
			}
			else if(i1m1 < 0)
			{
				if((im1m1 > 0) && (i0m1 > 0) && (im10 > 0) && (i00 > 0))
				{
					arResInd[0] = im1m1; arResInd[1] = i0m1; arResInd[2] = im10; arResInd[3] = i00;
				}
				else if((i00 > 0) && (i10 > 0) && (i01 > 0) && (i11 > 0))
				{
					arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
				}
				else if((im10 > 0) && (i00 > 0) && (im11 > 0) && (i01 > 0))
				{
					arResInd[0] = im10; arResInd[1] = i00; arResInd[2] = im11; arResInd[3] = i01;
				}
			}
			else if(i0m1 < 0)
			{
				if((im10 > 0) && (i00 > 0) && (im11 > 0) && (i01 > 0))
				{
					arResInd[0] = im10; arResInd[1] = i00; arResInd[2] = im11; arResInd[3] = i01;
				}
				else if((i00 > 0) && (i10 > 0) && (i01 > 0) && (i11 > 0))
				{
					arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
				}
			}
			else if(im1m1 < 0)
			{
				if((im10 > 0) && (i00 > 0) && (im11 > 0) && (i01 > 0))
				{
					arResInd[0] = im10; arResInd[1] = i00; arResInd[2] = im11; arResInd[3] = i01;
				}
				else if((i0m1 > 0) && (i1m1 > 0) && (i00 > 0) && (i10 > 0))
				{
					arResInd[0] = i0m1; arResInd[1] = i1m1; arResInd[2] = i00; arResInd[3] = i10;
				}
				else if((i00 > 0) && (i10 > 0) && (i01 > 0) && (i11 > 0))
				{
					arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
				}
			}
			else if(im10 < 0)
			{
				if((i0m1 > 0) && (i1m1 > 0) && (i00 > 0) && (i10 > 0))
				{
					arResInd[0] = i0m1; arResInd[1] = i1m1; arResInd[2] = i00; arResInd[3] = i10;
				}
				else if((i00 > 0) && (i10 > 0) && (i01 > 0) && (i11 > 0))
				{
					arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
				}
			}
			else if(im11 < 0)
			{
				if((im1m1 > 0) && (i0m1 > 0) && (im10 > 0) && (i00 > 0))
				{
					arResInd[0] = im1m1; arResInd[1] = i0m1; arResInd[2] = im10; arResInd[3] = i00;
				}
				else if((i00 > 0) && (i10 > 0) && (i01 > 0) && (i11 > 0))
				{
					arResInd[0] = i00; arResInd[1] = i10; arResInd[2] = i01; arResInd[3] = i11;
				}
				else if((i0m1 > 0) && (i1m1 > 0) && (i00 > 0) && (i10 > 0))
				{
					arResInd[0] = i0m1; arResInd[1] = i1m1; arResInd[2] = i00; arResInd[3] = i10;
				}
			}
			else if(i01 < 0)
			{
				if((im1m1 > 0) && (i0m1 > 0) && (im10 > 0) && (i00 > 0))
				{
					arResInd[0] = im1m1; arResInd[1] = i0m1; arResInd[2] = im10; arResInd[3] = i00;
				}
				else if((i0m1 > 0) && (i1m1 > 0) && (i00 > 0) && (i10 > 0))
				{
					arResInd[0] = i0m1; arResInd[1] = i1m1; arResInd[2] = i00; arResInd[3] = i10;
				}
			}

			for(int iii=0; iii<4; iii++)
			{
				if(arResInd[iii] < 0)
				{
					delete[] arMeshX; delete[] arMeshY; return -1;
				}
			}
		}
		else if(ord == 2)
		{
			for(int ii=0; ii<5; ii++) arResInd[ii] = -1; 

			if((ixm1 < 0) && (ixm2 < 0) && (iym1 < 0) && (iym2 < 0))
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1;
				delete[] arMeshX; delete[] arMeshY; return -1;
			}
			nResInd = 5;
			//double f0m1 = *(arF++);
			//double fm10 = *(arF++);
			//double f00 = *(arF++);
			//double f10 = *(arF++);
			//double f01 = *arF;

			if((ixm1 >= 0) && (ixp1 >= 0) && (iym1 >= 0) && (iyp1 >= 0))
			{
				arResInd[0] = iym1; arResInd[1] = ixm1; arResInd[2] = i0; arResInd[3] = ixp1; arResInd[4] = iyp1;
			}
			else if((iym1 < 0) && (iyp1 >= 0) && (iyp2 >= 0) && (ixm1 >= 0) && (ixp1 >= 0))
			{
				arResInd[0] = i0; arResInd[1] = ixm1; arResInd[2] = iyp1; arResInd[3] = ixp1; arResInd[4] = iyp2;
			}
			else if((iyp1 < 0) && (iym1 >= 0) && (iym2 >= 0) && (ixm1 >= 0) && (ixp1 >= 0))
			{
				arResInd[0] = iym2; arResInd[1] = ixm1; arResInd[2] = iym1; arResInd[3] = ixp1; arResInd[4] = i0;
			}
			else if((ixm1 < 0) && (ixp1 >= 0) && (ixp2 >= 0) && (iym1 >= 0) && (iyp1 >= 0))
			{
				arResInd[0] = iym1; arResInd[1] = i0; arResInd[2] = ixp1; arResInd[3] = ixp2; arResInd[4] = iyp1;
			}
			else if((ixp1 < 0) && (ixm1 >= 0) && (ixm2 >= 0) && (iym1 >= 0) && (iyp1 >= 0))
			{
				arResInd[0] = iym1; arResInd[1] = ixm2; arResInd[2] = ixm1; arResInd[3] = i0; arResInd[4] = iyp1;
			}
			else if((ixm1 < 0) && (iym1 < 0) && (ixp1 >= 0) && (ixp2 >= 0) && (iyp1 >= 0) && (iyp2 >= 0))
			{
				arResInd[0] = iyp1; arResInd[1] = ixp1; arResInd[2] = i0; arResInd[3] = ixp2; arResInd[4] = iyp2; //?
			}
			else if((ixm1 < 0) && (iyp1 < 0) && (ixp1 >= 0) && (ixp2 >= 0) && (iym1 >= 0) && (iym2 >= 0))
			{
				arResInd[0] = iym2; arResInd[1] = ixp1; arResInd[2] = i0; arResInd[3] = ixp2; arResInd[4] = iym1; //?
			}
			else if((ixp1 < 0) && (iym1 < 0) && (ixm1 >= 0) && (ixm2 >= 0) && (iyp1 >= 0) && (iyp2 >= 0))
			{
				arResInd[0] = iyp1; arResInd[1] = ixm2; arResInd[2] = i0; arResInd[3] = ixm1; arResInd[4] = iyp2; //?
			}
			else if((ixp1 < 0) && (iyp1 < 0) && (ixm1 >= 0) && (ixm2 >= 0) && (iym1 >= 0) && (iym2 >= 0))
			{
				arResInd[0] = iym2; arResInd[1] = ixm2; arResInd[2] = i0; arResInd[3] = ixm1; arResInd[4] = iym1; //?
			}
			else
			{
				//delete[] arMeshX; delete[] arMeshY; delete[] arPtMeshInd; return -1;
				delete[] arMeshX; delete[] arMeshY; return -1;
			}

			for(int iii=0; iii<5; iii++)
			{
				if(arResInd[iii] < 0)
				{
					delete[] arMeshX; delete[] arMeshY; return -1;
				}
			}

		}
	}

	delete[] arMeshX;
	delete[] arMeshY;
	//delete[] arPtMeshInd;
	return 2; //Normal return indicating 2d interpolation
}

//-------------------------------------------------------------------------

int CGenMathInterp::TryToFindMeshPointForPars(double* arPar1, double* arPar2, int nVals, double* arPrecPar)
{
	if((arPar1 == 0) && (arPar2 == 0)) throw CAN_NOT_FIND_IND_FOR_INTERP; //return; //throw exception?
	if(nVals < 1) throw CAN_NOT_FIND_IND_FOR_INTERP; //return;

	int dimInterp = (int)arPrecPar[0];
	double par1 = arPrecPar[1];
	double par2 = arPrecPar[2];

	if(arPar1 == 0)
	{
		arPar1 = arPar2;
		par1 = par2;
		dimInterp = 1;
	}
	if(arPar2 == 0) dimInterp = 1;

	const double relTolPar = 1.e-09;
	int nVals_mi_1 = nVals - 1;
	double tolPar1 = relTolPar*(::abs(arPar1[nVals_mi_1] - arPar1[0]));
	double tolPar2 = (dimInterp == 2)? relTolPar*(::abs(arPar2[nVals_mi_1] - arPar2[0])) : 0.;

	for(int i=0; i<nVals; i++)
	{
		double curMeshPar = arPar1[i];
		if(((curMeshPar - tolPar1) <= par1) && (par1 <= (curMeshPar + tolPar1)))
		{
			if(dimInterp == 2)
			{
				curMeshPar = arPar2[i];
				if(((curMeshPar - tolPar2) <= par2) && (par2 <= (curMeshPar + tolPar2))) return i;
			}
			else return i;
		}
	}
	return -1;
}

//-------------------------------------------------------------------------

void CGenMathInterp::SelectPointsForInterp1d2d(double* arPar1, double* arPar2, int nVals, int* arResInds, int& nResInds, double arPrecPar[5])
{//This calculates arResInds, nResInds, and can eventually modify arPrecPar
	nResInds = 0;

	if((arPar1 == 0) && (arPar2 == 0)) throw CAN_NOT_FIND_IND_FOR_INTERP; //return; //throw exception?
	if(arResInds == 0) throw CAN_NOT_FIND_IND_FOR_INTERP; //return; //throw exception?
	//if(nVals <= 1) throw CAN_NOT_FIND_IND_FOR_INTERP; //return;
	if(nVals < 1) throw CAN_NOT_FIND_IND_FOR_INTERP; //return;

	int indMeshPoint = TryToFindMeshPointForPars(arPar1, arPar2, nVals, arPrecPar);
	if((indMeshPoint >= 0) && (indMeshPoint < nVals))
	{
		arResInds[0] = indMeshPoint; nResInds = 1; 
		arPrecPar[3] = 0; //?
		return;
	}

	int dimInterp = (int)arPrecPar[0];
	double par1 = arPrecPar[1];
	double par2 = arPrecPar[2];
	int ordInterp = (int)arPrecPar[3];
	int meshIsRect = (int)arPrecPar[4];

	if(arPar1 == 0)
	{
		arPar1 = arPar2;
		par1 = par2;
		dimInterp = 1;
	}
	if(arPar2 == 0) dimInterp = 1;

	if(dimInterp == 2)
	{
		bool resMeshIsRect = false;
		int dimInterpNew = SelectPointsForInterp2d(par1, par2, arPar1, arPar2, nVals, ordInterp, arResInds, nResInds, resMeshIsRect);
		if(dimInterpNew <= 0) return;

		dimInterp = dimInterpNew;
		meshIsRect = resMeshIsRect? 1 : 0;
	}

	if(dimInterp == 1)
	{//This assumes that arPar1 values are in the increasing order
		meshIsRect = 1;

		int i0 = -1;
		for(int i=0; i<nVals; i++)
		{
			if(par1 < arPar1[i]) { i0 = i - 1; break;}
		}

		if(i0 < 0)
		{
			const double relTolPar1 = 1.e-09;
			int nVals_mi_1 = nVals - 1;
			double tolPar1 = relTolPar1*(::abs(arPar1[nVals_mi_1] - arPar1[0])); //::abs() is required for compilation with GCC (OC16112017)

			if((par1 <= arPar1[0]) && (par1 >= (arPar1[0] - tolPar1))) i0 = 0;
			else if((par1 >= arPar1[nVals_mi_1]) && (par1 <= (arPar1[nVals_mi_1] + tolPar1))) i0 = nVals - 2;
		}

		if(ordInterp == 2)
		{
			if(nVals < 3) ordInterp = 1;
		}
		else if(ordInterp == 3)
		{
			if(nVals < 4) ordInterp = 2;
			else if(nVals < 3) ordInterp = 1;
		}
		else if(ordInterp > 3)
		{
			if(nVals >= 4) ordInterp = 3;
			else if(nVals >= 3) ordInterp = 2;
			else ordInterp = 1;
		}

		int nVals_mi_2 = nVals - 2;
		if(ordInterp == 1)
		{
			if(i0 < 0) i0 = 0;
			else if(i0 > nVals_mi_2) i0 = nVals_mi_2;

			arResInds[0] = i0;
			arResInds[1] = i0 + 1;
			nResInds = 2;
		}
		else if(ordInterp == 2)
		{
			if(i0 < 1) i0 = 1;
			else if(i0 > nVals_mi_2) i0 = nVals_mi_2;

			arResInds[0] = i0 - 1;
			arResInds[1] = i0;
			arResInds[2] = i0 + 1;
			nResInds = 3;
		}
		else //if(ordInterp == 3)
		{
			ordInterp = 3;
			int nVals_mi_3 = nVals - 3;
			if(i0 < 1) i0 = 1;
			else if(i0 > nVals_mi_3) i0 = nVals_mi_3;

			arResInds[0] = i0 - 1;
			arResInds[1] = i0;
			arResInds[2] = i0 + 1;
			arResInds[3] = i0 + 2;
			nResInds = 4;
		}
	}

	arPrecPar[0] = (double)dimInterp;
	arPrecPar[3] = (double)ordInterp;
	arPrecPar[4] = (double)meshIsRect;
}

//-------------------------------------------------------------------------
