/************************************************************************//**
 * File: gmmeth.cpp
 * Description: Auxiliary mathematical utilities
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "gmmeth.h"
//#include <math.h>

//-------------------------------------------------------------------------
// General mathematical methods
//-------------------------------------------------------------------------

bool CGenMathMeth::VectCheckIfCollinear(double xV1, double yV1, double zV1, double xV2, double yV2, double zV2, double RelPrec)
{
	double V1scV2 = xV1*xV2 + yV1*yV2 + zV1*zV2;
	double V1scV2e2 = V1scV2*V1scV2;

	double V1e2 = xV1*xV1 + yV1*yV1 + zV1*zV1;
	double V2e2 = xV2*xV2 + yV2*yV2 + zV2*zV2;
	double V1e2V2e2 = V1e2*V2e2;

	if(::fabs(::fabs(V1scV2e2) - ::fabs(V1e2V2e2)) < RelPrec*V1e2V2e2) return true;
	return false;
}

//-------------------------------------------------------------------------

double CGenMathMeth::Integ1D_FuncWithEdgeDer(double (*pF)(double, void*), double x1, double x2, double dFdx1, double dFdx2, double RelPrec, void* pv) //OC20112018
//double CGenMathMeth::Integ1D_FuncWithEdgeDer(double (*pF)(double), double x1, double x2, double dFdx1, double dFdx2, double RelPrec)
{
	//if((pF == 0) || (RelPrec <= 0)) return 0;
	if((pF == 0) || (RelPrec <= 0) || (x1 == x2)) return 0;

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	//long NpOnLevel = 5; // Must be non-even!
	long long NpOnLevel = 5; // Must be non-even!
	double xStep = (x2 - x1)/(NpOnLevel - 1);

	double Sum1=0., Sum2=0.;
	double wF;
	int LevelNo=0;

	wF = (*pF)(x1, pv); //OC20112018
	//wF = (*pF)(x1);
    double x = x1 + xStep;
	//int AmOfPass = (NpOnLevel - 3) >> 1;
	long long AmOfPass = (NpOnLevel - 3) >> 1;
	for(int i=0; i<AmOfPass; i++)
	{
        Sum1 += (*pF)(x, pv); x += xStep; //OC20112018
		//Sum1 += (*pF)(x); x += xStep;
		Sum2 += (*pF)(x, pv); x += xStep; //OC20112018
		//Sum2 += (*pF)(x); x += xStep;
	}
	Sum1 += (*pF)(x, pv); x += xStep; //OC20112018
	//Sum1 += (*pF)(x); x += xStep;
	wF += (*pF)(x, pv); //OC20112018
	//wF += (*pF)(x);
	wF *= wfe;

	double DifDer = dFdx1 - dFdx2;
	double wDifDer = wd*DifDer;
	double Int = xStep*(wF + wf1*Sum1 + wf2*Sum2 + xStep*wDifDer);

	NpOnLevel--;
	char NotFinishedYet = 1;
	char ExtraPassForAnyCase = 0; 
	while(NotFinishedYet)
	{
		Sum2 += Sum1;
		Sum1 = 0.;

		//char ThisMayBeTheLastLoop = 0;
		LevelNo++;

		double HalfStep = 0.5*xStep;
		x = x1 + HalfStep;

		//long NpOnLevel_mi_1 = NpOnLevel - 1;
		//for(int i=0; i<NpOnLevel; i++)
		for(long long i=0; i<NpOnLevel; i++)
		{
			Sum1 += (*pF)(x, pv); x += xStep; //OC20112018
			//Sum1 += (*pF)(x); x += xStep;
		}

		double LocInt = HalfStep*(wF + wf1*Sum1 + wf2*Sum2 + HalfStep*wDifDer);

        double TestVal = ::fabs(LocInt - Int);
        char SharplyGoesDown = (::fabs(LocInt) < 0.2*::fabs(Int));
        char NotFinishedYetFirstTest = (TestVal > RelPrec*(::fabs(LocInt)));

		if(!NotFinishedYetFirstTest)
		{
            if(ExtraPassForAnyCase || SharplyGoesDown) NotFinishedYet = 0;
            else ExtraPassForAnyCase = 1;
		}

		Int = LocInt;
		xStep = HalfStep; NpOnLevel *= 2;
	}
	return Int;
}

//-------------------------------------------------------------------------

double CGenMathMeth::Integ1D_Func(double (*pF)(double, void*), double x1, double x2, double RelPrec, void* pv) //OC02122018
{
	if((pF == 0) || (RelPrec <= 0) || (x1 == x2)) return 0;

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	//long NpOnLevel = 5; // Must be non-even!
	long long NpOnLevel = 5; // Must be non-even!
	double xStep = (x2 - x1)/(NpOnLevel - 1);

	double dF1, dF2;
	double Sum1=0., Sum2=0.;
	double wF;
	int LevelNo=0;

	wF = (*pF)(x1, pv); //OC20112018
	double Fedge1 = wF;

	//wF = (*pF)(x1);
    double x = x1 + xStep;
	//int AmOfPass = (NpOnLevel - 3) >> 1;
	long long AmOfPass = (NpOnLevel - 3) >> 1;
	for(int i=0; i<AmOfPass; i++)
	{
        Sum1 += (*pF)(x, pv); x += xStep; //OC20112018
		//Sum1 += (*pF)(x); x += xStep;
		if(i == 0) dF1 = Sum1 - Fedge1;

		Sum2 += (*pF)(x, pv); x += xStep; //OC20112018
		//Sum2 += (*pF)(x); x += xStep;
	}
	double auxF = (*pF)(x, pv);
	Sum1 += auxF; x += xStep; //OC02122018
	//Sum1 += (*pF)(x, pv); x += xStep; //OC20112018
	//Sum1 += (*pF)(x); x += xStep;
	dF2 = -auxF;

	double Fedge2 = (*pF)(x, pv);
	wF += Fedge2; //OC02122018
	//wF += (*pF)(x, pv); //OC20112018
	//wF += (*pF)(x);
	dF2 += Fedge2;
	wF *= wfe;

	double wDifFedge = wd*(dF1 - dF2); 
	double Int = xStep*(wF + wf1*Sum1 + wf2*Sum2 + wDifFedge);
	//double DifDer = dFdx1 - dFdx2;
	//double wDifDer = wd*DifDer;
	//double Int = xStep*(wF + wf1*Sum1 + wf2*Sum2 + xStep*wDifDer);

	NpOnLevel--;
	char NotFinishedYet = 1;
	char ExtraPassForAnyCase = 0; 
	while(NotFinishedYet)
	{
		Sum2 += Sum1;
		//Sum1 = 0.;

		//char ThisMayBeTheLastLoop = 0;
		LevelNo++;

		double HalfStep = 0.5*xStep;
		x = x1 + HalfStep;

		Sum1 = (*pF)(x, pv); x += xStep; //OC02122018
		dF1 = Sum1 - Fedge1;

		for(long long i=1; i<(NpOnLevel-1); i++) //OC02122018
		//for(long long i=0; i<NpOnLevel; i++)
		{
			Sum1 += (*pF)(x, pv); x += xStep;
		}

		auxF = (*pF)(x, pv);
		Sum1 += auxF; //OC02122018
		dF2 = Fedge2 - auxF;
		wDifFedge = wd*(dF1 - dF2); 
		double LocInt = HalfStep*(wF + wf1*Sum1 + wf2*Sum2 + wDifFedge);
		//double LocInt = HalfStep*(wF + wf1*Sum1 + wf2*Sum2 + HalfStep*wDifDer);

        double TestVal = ::fabs(LocInt - Int);
        char SharplyGoesDown = (::fabs(LocInt) < 0.2*::fabs(Int));
        char NotFinishedYetFirstTest = (TestVal > RelPrec*(::fabs(LocInt)));

		if(!NotFinishedYetFirstTest)
		{
            if(ExtraPassForAnyCase || SharplyGoesDown) NotFinishedYet = 0;
            else ExtraPassForAnyCase = 1;
		}

		Int = LocInt;
		xStep = HalfStep; NpOnLevel *= 2;
	}
	return Int;
}

//-------------------------------------------------------------------------
//
//double CGenMathMeth::Integ1D_FuncDefByArray(double* FuncArr, long Np, double Step)
//{
//	if((FuncArr == 0) || (Np < 2) || (Step == 0)) return 0;
//
//	double *tFuncArr = FuncArr + 1;
//	bool NpIsEven = (Np == ((Np >> 1) << 1));
//	if(NpIsEven)
//	{//apply method of trapeth
//		double Sum = 0;
//		for(long i=1; i<(Np - 1); i++)
//		{
//            Sum += *(tFuncArr++);
//		}
//		Sum += 0.5*((*tFuncArr) + FuncArr[0]);
//		return Sum*Step;
//	}
//	else
//	{//apply Simplson method, to check!!!
//		double Sum1 = 0, Sum2 = 0;
//		for(long i=1; i<((Np - 3) >> 1); i++)
//		{
//            Sum1 += *(tFuncArr++);
//            Sum2 += *(tFuncArr++);
//		}
//		Sum1 += *(tFuncArr++);
//		return (Step/3.)*(FuncArr[0] + 4.*Sum1 + 2.*Sum2 + (*tFuncArr));
//	}
//}
//
//-------------------------------------------------------------------------
//
//double CGenMathMeth::Integ1D_FuncDefByArray(float* FuncArr, long Np, double Step)
//{
//	if((FuncArr == 0) || (Np < 2) || (Step == 0)) return 0;
//
//	float *tFuncArr = FuncArr + 1;
//	bool NpIsEven = (Np == ((Np >> 1) << 1));
//	if(NpIsEven)
//	{//apply method of trapeth
//		double Sum = 0;
//		for(long i=1; i<(Np - 1); i++)
//		{
//            Sum += *(tFuncArr++);
//		}
//		Sum += 0.5*((*tFuncArr) + FuncArr[0]);
//		return Sum*Step;
//	}
//	else
//	{//apply Simplson method, to check!!!
//		double Sum1 = 0, Sum2 = 0;
//		for(long i=1; i<((Np - 3) >> 1); i++)
//		{
//            Sum1 += *(tFuncArr++);
//            Sum2 += *(tFuncArr++);
//		}
//		Sum1 += *(tFuncArr++);
//		//return (Step/3.)*(FuncArr[0] + 4.*Sum1 + 2.*Sum1 + (*tFuncArr));
//		return (Step/3.)*(FuncArr[0] + 4.*Sum1 + 2.*Sum2 + (*tFuncArr));
//	}
//}
//
//-------------------------------------------------------------------------


