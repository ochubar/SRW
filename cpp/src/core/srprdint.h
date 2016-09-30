/************************************************************************//**
 * File: srprdint.h
 * Description: Auxiliary (obsolete?) SR calculation method from ~Arbitrary Transversely-Uniform Magnetic Field, in Near-Field observation conditions (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRPRTINT_H
#define __SRPRTINT_H

#include "srstraux.h"
#include "srtrjdat.h"

//*************************************************************************

class srTPartAutoRadInt {

	double AuxDataArray[1000];
	double* AuxDataPtr;

	double *BtxArrP[50], *XArrP[50], *IntBtxE2ArrP[50], *BxArrP[50];
	double *BtzArrP[50], *ZArrP[50], *IntBtzE2ArrP[50], *BzArrP[50];
	//int AmOfPointsOnLevel[50];
	long long AmOfPointsOnLevel[50];
	//int NumberOfLevelsFilled, NpOnZeroLevel;
	long long NumberOfLevelsFilled, NpOnZeroLevel;

	double sIntegStart, sIntegFin, sIntegStep, sIntegRelPrec;
	double MaxFluxDensVal, CurrentAbsPrec;

	srTTrjDat* TrjDatPtr;
	srTWfrSmp* DistrInfoDatPtr;
	srLambXYZ* ObsCoorPtr;

	//srTSend Send;

	double wfe, wf1, wf2, wd;
	double PI, TwoPI, ThreePIdTwo, FivePIdFour, HalfPI, One_dTwoPI, PIm10e6, PIm10e6dEnCon, TenEm6dTwoPI; // Constants
	double a2c, a4c, a6c, a8c, a10c, a12c;
	double a3s, a5s, a7s, a9s, a11s, a13s;

public:

	char SomethingIsWrong;

	//srTPartAutoRadInt(srTTrjDat* InTrjDatPtr, srTWfrSmp* InDistrInfoDatPtr, srLambXYZ* InObsCoorPtr, double In_sStart, double In_sFin, double In_sRelPrec, int InNpOnZeroLevel);
	srTPartAutoRadInt(srTTrjDat* InTrjDatPtr, srTWfrSmp* InDistrInfoDatPtr, srLambXYZ* InObsCoorPtr, double In_sStart, double In_sFin, double In_sRelPrec, long long InNpOnZeroLevel);

	~srTPartAutoRadInt() 
	{
		if(AuxDataPtr != 0) delete[] AuxDataPtr;
		for(int i=3; i<NumberOfLevelsFilled; i++) delete[] (BtxArrP[i]);
	}

	int Integrate(double& IntXRe, double& IntXIm, double& IntZRe, double& IntZIm, double* InEdgeFunDer =0);

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

	//int FillNextLevel(int LevelNo, double sStart, double sEnd, long Np)
	int FillNextLevel(int LevelNo, double sStart, double sEnd, long long Np)
	{
		double* BasePtr;
		if(LevelNo == 2) BasePtr = (AuxDataPtr == 0)? &(AuxDataArray[(2*NpOnZeroLevel - 1)*8]) : &(AuxDataPtr[(2*NpOnZeroLevel - 1)*8]);
		else
		{
			BasePtr = new double[Np*8];
			//if(BasePtr == 0) { Send.ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
			if(BasePtr == 0) { return MEMORY_ALLOCATION_FAILURE;}
		}
		BtxArrP[LevelNo] = BasePtr;
		BasePtr += Np; XArrP[LevelNo] = BasePtr;
		BasePtr += Np; IntBtxE2ArrP[LevelNo] = BasePtr;
		BasePtr += Np; BxArrP[LevelNo] = BasePtr;
		BasePtr += Np; BtzArrP[LevelNo] = BasePtr;
		BasePtr += Np; ZArrP[LevelNo] = BasePtr;
		BasePtr += Np; IntBtzE2ArrP[LevelNo] = BasePtr;
		BasePtr += Np; BzArrP[LevelNo] = BasePtr;
		TrjDatPtr->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevelNo], BtzArrP[LevelNo], XArrP[LevelNo], ZArrP[LevelNo], IntBtxE2ArrP[LevelNo], IntBtzE2ArrP[LevelNo], BxArrP[LevelNo], BzArrP[LevelNo]);
		AmOfPointsOnLevel[LevelNo] = Np;
		NumberOfLevelsFilled++;
		return 0;
	}
};

//*************************************************************************

#endif
