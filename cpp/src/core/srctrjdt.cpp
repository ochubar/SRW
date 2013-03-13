/************************************************************************//**
 * File: srctrjdt.cpp
 * Description: Electron trajectory calculation, constant magnetic field case
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srtrjdat.h"
#include "gmtrans.h"
#include "srctrjdt.h"

//*************************************************************************

int srTConstTrjDat::ConvertToArbTrjDat(char Periodicity, srTWfrSmp& DistrInfoDat, srTGenTrjHndl& TrjHndl)
{// This performs only the actions that are normally performed at setting up srTTrjDat via srTSend !!!
 // if(Periodicity == 1) - Arb. field
 // if(Periodicity == 2) - Periodic field
	//int result;
	const int AmOfFieldPoints = 300; // To steer

	srTTrjDat* pNewTrjDat = new srTTrjDat();
	if(pNewTrjDat == 0) return MEMORY_ALLOCATION_FAILURE;

//El. Beam
	pNewTrjDat->EbmDat = EbmDat;

//Mag. Field
	pNewTrjDat->HorFieldIsNotZero = HorFieldIsNotZero;
	pNewTrjDat->VerFieldIsNotZero = VerFieldIsNotZero;

	pNewTrjDat->LenFieldData = AmOfFieldPoints;
	pNewTrjDat->AuxBxLen = AmOfFieldPoints;
	pNewTrjDat->AuxBzLen = AmOfFieldPoints;

	pNewTrjDat->BxInData = new srTFunDer[AmOfFieldPoints];
	if(pNewTrjDat->BxInData == 0) return MEMORY_ALLOCATION_FAILURE;
	pNewTrjDat->BzInData = new srTFunDer[AmOfFieldPoints];
	if(pNewTrjDat->BzInData == 0) return MEMORY_ALLOCATION_FAILURE;

	pNewTrjDat->NperTot = 1;
	pNewTrjDat->NperLeft = 0;

	double sStart, sEnd;
	DetermineIntegLimitsForArbTrj(DistrInfoDat, sStart, sEnd);

	pNewTrjDat->sStep = (sEnd - sStart)/(AmOfFieldPoints - 1);
	pNewTrjDat->Inv_Step = 1./pNewTrjDat->sStep;
	pNewTrjDat->sStart = sStart;

	//double FieldZeroTol = pNewTrjDat->FieldZeroTolerance;
	srTFunDer *tFunDerX = pNewTrjDat->BxInData;
	srTFunDer *tFunDerZ = pNewTrjDat->BzInData;
	double Bx = MagConst.Bx, Bz = MagConst.Bz;
	for(long i=0; i<AmOfFieldPoints; i++)
	{
		tFunDerX->f = Bx;
		(tFunDerX++)->dfds = 0.;
		tFunDerZ->f = Bz;
		(tFunDerZ++)->dfds = 0.;
	}

	TrjHndl = srTGenTrjHndl(pNewTrjDat);
	return 0;
}

//*************************************************************************

void srTConstTrjDat::DetermineIntegLimitsForArbTrj(srTWfrSmp& DistrInfoDat, double& sStart, double& sEnd)
{
	const int AmOfSafetyInterv = 10; // To steer
	const double RelTolNotRotate = 1.E-05; // To steer

	TVector3d uLab(MagConst.Bx, 0., MagConst.Bz);
	double Bcon = sqrt(uLab.x*uLab.x + uLab.z*uLab.z);
	double Rmag = 3.33564076253*(EbmDat.Energy)/Bcon;
	double ds = 0;

	if(DistrInfoDat.AssumeAllPhotonEnergies)
	{
        ds = Rmag*sqrt(EbmDat.GammaEm2);
	}
	else
	{
		double LambMax_m;
		if(DistrInfoDat.TreatLambdaAsEnergyIn_eV)
		{
			double eMin = DistrInfoDat.LambStart;
			LambMax_m = 1.24E-06/eMin;
		}
		else
		{
			LambMax_m = 1.E-09*((DistrInfoDat.nLamb > 1)? DistrInfoDat.LambEnd : DistrInfoDat.LambStart);
		}
        double ds1 = pow(LambMax_m*Rmag*Rmag, 1./3.);
        double ds2 = LambMax_m/EbmDat.GammaEm2;

		//ds = (ds1 > ds2)? ds1 : ds2;
		ds = (ds1 < ds2)? ds1 : ds2;
	}

	uLab = (1./Bcon)*uLab;
	TVector3d uLoc(0.,0.,1.), Zero(0.,0.,0.);
	TVector3d uLab_mi_uLoc = uLab - uLoc;
	//double NormEst = Abs(uLab_mi_uLoc);
	double NormEst = uLab_mi_uLoc.Abs();

	gmTrans Rot;
	if(NormEst < RelTolNotRotate)
	{
		Rot.SetupIdent();
	}
	else
	{
		TVector3d AxVect = uLab^uLoc;
		double uLab_uLoc = uLab*uLoc;
		double Angle = acos(uLab_uLoc);
		
		Rot.SetupRotation(Zero, AxVect, Angle);
		TVector3d TestVect = Rot.TrBiPoint(uLab);
		double TestScal = TestVect*uLab;
		if(TestScal < 0.)
		{
			Rot.SetupRotation(Zero, AxVect, -Angle);
		}
	}

	TVector3d P1(DistrInfoDat.xStart, 0., DistrInfoDat.zStart);
	P1 = Rot.TrPoint(P1);
	double xLocMin = P1.x, xLocMax = P1.x;
	if(DistrInfoDat.nx > 1)
	{
		TVector3d P2(DistrInfoDat.xEnd, 0., DistrInfoDat.zStart);
		P2 = Rot.TrPoint(P2);
		if(xLocMin > P2.x) xLocMin = P2.x;
		if(xLocMax < P2.x) xLocMax = P2.x;
	}
	if(DistrInfoDat.nz > 1)
	{
		TVector3d P3(DistrInfoDat.xStart, 0., DistrInfoDat.zEnd);
		P3 = Rot.TrPoint(P3);
		if(xLocMin > P3.x) xLocMin = P3.x;
		if(xLocMax < P3.x) xLocMax = P3.x;

		if(DistrInfoDat.nx > 1)
		{
			TVector3d P4(DistrInfoDat.xEnd, 0., DistrInfoDat.zEnd);
			P4 = Rot.TrPoint(P4);
			if(xLocMin > P4.x) xLocMin = P4.x;
			if(xLocMax < P4.x) xLocMax = P4.x;
		}
	}

	TVector3d InitCoord(EbmDat.x0, 0., EbmDat.z0);
	InitCoord = Rot.TrPoint(InitCoord);
	TVector3d InitAng(EbmDat.dxds0, 0., EbmDat.dzds0);
	InitAng = Rot.TrPoint(InitAng);

	double y1Inv = 1./(DistrInfoDat.yStart - EbmDat.s0);
	double xAngMin = (xLocMin - InitCoord.x)*y1Inv - InitAng.x;
	double xAngMax = (xLocMax - InitCoord.x)*y1Inv - InitAng.x;

	sStart = EbmDat.s0 + xAngMin*Rmag - AmOfSafetyInterv*ds;
	sEnd = EbmDat.s0 + xAngMax*Rmag + AmOfSafetyInterv*ds;
}

//*************************************************************************

srTGenTrjDat* srTMagFieldConstant::CreateAndSetupNewTrjDat(srTEbmDat* pEbmDat) //virtual
{
	srTConstTrjDat* pOut = new srTConstTrjDat();
	pOut->MagConst = *this;
	if(pEbmDat != 0) pOut->EbmDat = *pEbmDat;
	pOut->CheckIfHorOrVertFieldIsZero();
	return pOut;
}

//*************************************************************************
