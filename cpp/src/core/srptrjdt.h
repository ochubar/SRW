
#ifndef __SRPTRJDT_H
#define __SRPTRJDT_H
/************************************************************************//**
 * File: srptrjdt.h
 * Description: Calculation of Electron Trajectory in Periodic Magnetic Field of an ID (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srmagfld.h"
#include "srebmdat.h"
#include "srgtrjdt.h"

//*************************************************************************

class srTPerTrjDat : public srTGenTrjDat {
public:

	srTMagFieldPeriodic MagPer;

	srTPerTrjDat() {}

	void CompTrjDataDerivedAtPointPowDens(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz);

	void CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays);
	//void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);
	void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);
	//void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);
	void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);

	int SetUpFieldBasedArraysAtOnePeriod(srTFieldBasedArrayKeys&, srTFieldBasedArrays&);
	int SetUpFieldBasedArraysTotal(srTFieldBasedArrayKeys&, srTFieldBasedArrays&);

	int ConvertToArbTrjDat(char, srTWfrSmp&, srTGenTrjHndl&);
	//int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long& Ns, int& NperTot, int& NperLeft);
	int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long long& Ns, int& NperTot, int& NperLeft);

	int MagFieldPeriodicity() { return 2;} // Periodic
	char MagFieldIsConstant() { return 0;}

	void AnalizeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ);
	void ShowFullLimits(double& sIntegStart, double& sIntegFin)
	{
		int Nper = int(MagPer.TotLength/MagPer.PerLength);
		if(((Nper >> 1) << 1) != Nper) // Odd
		{
			int NperLeftRight = (Nper - 1) >> 1;
			sIntegStart = -(0.5 + NperLeftRight)*MagPer.PerLength;
		}
		else // Even
		{
			int NperLeft = Nper >> 1;
			sIntegStart = -NperLeft*MagPer.PerLength;
		}
		sIntegFin = -sIntegStart;
	}

	int ShowLimitsAndInitInteg(srTWfrSmp&, char LongIntType, double& sIntegStart, double& sIntegFin, int& AmOfPer, bool doInit = true) 
	{
		if(LongIntType == 1) // Total range
		{
			AmOfPer = int(MagPer.TotLength/MagPer.PerLength);
			double ActualTotLength = AmOfPer*MagPer.PerLength;
			sIntegStart = -(AmOfPer >> 1)*MagPer.PerLength;
			sIntegFin = sIntegStart + ActualTotLength;
			AmOfPer = 1;
		}
		else if(LongIntType == 2) // One period
		{
			sIntegStart = 0.;
			sIntegFin = MagPer.PerLength;
			AmOfPer = int(MagPer.TotLength/MagPer.PerLength);;
		}
		return 0;
	}
	int InitTrjComp() 
	{ 
		CheckIfHorOrVertFieldIsZero();
		CompBetaNormConst();
		return 0;
	}
	void CheckIfHorOrVertFieldIsZero()
	{
		char LocHorFieldIsNotZero = 0, LocVerFieldIsNotZero = 0;
		for(int i=0; i<MagPer.AmOfHarm; i++)
		{
			srTMagHarm& Harm = (MagPer.HarmVect)[i];
			if(Harm.XorZ == 'x') LocHorFieldIsNotZero = 1;
			if(Harm.XorZ == 'z') LocVerFieldIsNotZero = 1;
		}
		HorFieldIsNotZero = LocHorFieldIsNotZero;
		VerFieldIsNotZero = LocVerFieldIsNotZero;
	}
	double MaxAbsField()
	{
		double Bcon = 0.010709839006/MagPer.PerLength;
		double MaxB = 0.;
		for(int i=0; i<MagPer.AmOfHarm; i++)
		{
			srTMagHarm& Harm = (MagPer.HarmVect)[i];
			double B0 = ::fabs(Harm.K*Bcon);
			if(MaxB < B0) MaxB = B0;
		}
		return MaxB;
	}
	void CountFieldExtremums()
	{// Improve sometime
		if(HorFieldIsNotZero) AmOfExtremInBx = 2;
		else AmOfExtremInBx = 0;
		if(VerFieldIsNotZero) AmOfExtremInBz = 2;
		else AmOfExtremInBz = 0;
	}

	//long EstimMinNpForRadInteg(char typeInt) //	virtual 
	long long EstimMinNpForRadInteg(char typeInt) //	virtual 
	{//treat eventually monochromatic case (typeInt == 1)
		const int nPtPerPeriod = 4; //5;

		int nPerEstim = 1;
		if((MagPer.PerLength > 0) && (MagPer.TotLength > 0)) nPerEstim = (int)(MagPer.TotLength/MagPer.PerLength) + 1;

		return nPtPerPeriod*nPerEstim*(MagPer.AmOfHarm);
	}
};

//*************************************************************************

#endif
