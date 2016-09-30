/************************************************************************//**
 * File: srthckbm2.cpp
 * Description: SR Stokes parameters calculation method for ~Arbitrary magnetic field (used rarely, under-programmed)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srthckbm.h"
#include "srmagfld.h"
#include "srradmnp.h"
#include "gmrand.h"
#include "srradint.h"

//*************************************************************************

void srTRadIntThickBeam::EstimateExtraObservRangesToIncludeEmittance(srTStokesStructAccessData* pStokes, srTEbmDat* pElecBeam, double* ArrExtraRanges)
{
	ArrExtraRanges[0] = ArrExtraRanges[1] = 0;
	double Length = pStokes->yStart - pElecBeam->s0;
	double p4x4PropMatr[] = {
        1., Length, 0., 0.,
		0., 1., 0., 0.,
        0., 0., 1., Length,
		0., 0., 0., 1.,
	};
	double p4Vect[] = { 0., 0., 0., 0.};
	srTElecBeamMoments ElecBeamMom(pElecBeam);
	srTRadGenManip::PropagateElecBeamMoments(ElecBeamMom, p4x4PropMatr, p4Vect);

	const double SigmaExtentFactor = 3; // to tune
	ArrExtraRanges[0] = SigmaExtentFactor*(::sqrt(ElecBeamMom.Mxx)); // extension for one side
    ArrExtraRanges[1] = SigmaExtentFactor*(::sqrt(ElecBeamMom.Mzz));
}

//*************************************************************************

srTSRWRadStructAccessData* srTRadIntThickBeam::CreateNewRadStructWithConstParams(srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTStokesStructAccessData* pStokes, srTWfrSmp*& pWfrSmp)
{
	double ArrExtraRanges[2];
	EstimateExtraObservRangesToIncludeEmittance(pStokes, pElecBeam, ArrExtraRanges);
	double &HalfExtraRangeX = ArrExtraRanges[0], &HalfExtraRangeZ = ArrExtraRanges[1];

	int NpMin = 31;

	double xStep = pStokes->xStep;
	if(xStep <= 0) xStep = HalfExtraRangeX/NpMin;
	double xRangeStokes = (pStokes->nx - 1)*xStep;
	int NumStepsInHalfExtraRangeX = (int)(HalfExtraRangeX/xStep) + 1;
	double ActualHalfExtraRangeX = NumStepsInHalfExtraRangeX*xStep;

	double hSt = pStokes->xStart - ActualHalfExtraRangeX;
	int hN = pStokes->nx + 2*NumStepsInHalfExtraRangeX;
    double hFi = pStokes->xStart + xRangeStokes + ActualHalfExtraRangeX;
    //double hSt = pStokes->xStart;
    //int hN = pStokes->nx;
    //double hFi = hSt + xRangeStokes;
	
	double zStep = pStokes->zStep;
	if(zStep <= 0) zStep = HalfExtraRangeZ/NpMin;
	double zRangeStokes = (pStokes->nz - 1)*zStep;
	int NumStepsInHalfExtraRangeZ = (int)(HalfExtraRangeZ/zStep) + 1;
	double ActualHalfExtraRangeZ = NumStepsInHalfExtraRangeZ*zStep;

	double vSt = pStokes->zStart - ActualHalfExtraRangeZ;
	int vN = pStokes->nz + 2*NumStepsInHalfExtraRangeZ;
    double vFi = pStokes->zStart + zRangeStokes + ActualHalfExtraRangeZ;
	//double vSt = pStokes->zStart;
	//int vN = pStokes->nz;
	//double vFi = vSt + zRangeStokes;

	double eSt = pStokes->eStart;
	int eN = pStokes->ne;
    double eFi = eSt + (eN - 1)*(pStokes->eStep);

	//pWfrSmp = new srTWfrSmp(pStokes->yStart, hSt, hFi, hN, vSt, vFi, vN, eSt, eFi, eN, "EV");
	pWfrSmp = new srTWfrSmp(pStokes->yStart, hSt, hFi, hN, vSt, vFi, vN, 0, eSt, eFi, eN, "EV");
	return new srTSRWRadStructAccessData(pElecBeam, pTrjDat, pWfrSmp, 0);
}

//*************************************************************************

double srTRadIntThickBeam::GetNextElecEnergyFromGausDistrib(srTEbmDat& OrigElecBeam, CGenMathRand& RandGen)
{
	if(m_SpareElecEnergyVal > 0) 
	{
		double OutVal = m_SpareElecEnergyVal;
		m_SpareElecEnergyVal = 0;
		return OutVal;
	}

	double newE1=0, newE2=0;
	//RandGen.NextTwoGaussRand1D(OrigElecBeam.Energy, OrigElecBeam.GetSigmaE_GeV(), newE1, newE2);
	RandGen.NextRandGauss2(OrigElecBeam.Energy, OrigElecBeam.GetSigmaE_GeV(), newE1, newE2);
	m_SpareElecEnergyVal = newE2;
	return newE1;

	//return RandGen.NextGaussRandP(OrigElecBeam.Energy, OrigElecBeam.GetSigmaE_GeV());
}

//*************************************************************************

//double srTRadIntThickBeam::UpdateResultStokesData(float* ArrAuxDataS0, float* ArrAuxDataS1, float* ArrAuxDataS2, float* ArrAuxDataS3, srTWfrSmp* pWfrSmp, int iElecEn, srTStokesStructAccessData* pStokes)
double srTRadIntThickBeam::UpdateResultStokesData(float* ArrAuxDataS0, float* ArrAuxDataS1, float* ArrAuxDataS2, float* ArrAuxDataS3, srTWfrSmp* pWfrSmp, long long iElecEn, srTStokesStructAccessData* pStokes)
{//this assumes that Stokes S0 is always defined
	double RelPrec = 1E+23;

	int xAmStepsLeft = (int)((pStokes->xStart - pWfrSmp->xStart)/(pStokes->xStep) + 0.000001);
	int xAmStepsRight = (int)((pWfrSmp->xEnd - (pStokes->xStart + (pStokes->nx - 1)*(pStokes->xStep)))/(pStokes->xStep) + 0.000001);
	int zAmStepsLeft = (int)((pStokes->zStart - pWfrSmp->zStart)/(pStokes->zStep) + 0.000001);
	//int zAmStepsRight = (int)((pWfrSmp->zEnd - (pStokes->zStart + (pStokes->nz - 1)*(pStokes->zStep)))/(pStokes->zStep) + 0.000001);

	//int Offset0 = (zAmStepsLeft*(pWfrSmp->nx) + xAmStepsLeft)*(pWfrSmp->nLamb);
	long long Offset0 = (((long long)zAmStepsLeft)*((long long)(pWfrSmp->nx)) + (long long)xAmStepsLeft)*((long long)(pWfrSmp->nLamb));
	int ShiftBetweenRowsX = (xAmStepsRight + xAmStepsLeft)*(pWfrSmp->nLamb);

	float *tBaseSto = pStokes->pBaseSto;
	float *tAuxDataS0 = 0, *tAuxDataS1 = 0, *tAuxDataS2 = 0, *tAuxDataS3 = 0;
	bool S0_IsDefined = false, S1_IsDefined = false, S2_IsDefined = false, S3_IsDefined = false;

	tAuxDataS0 = ArrAuxDataS0 + Offset0; S0_IsDefined = true;
	if(ArrAuxDataS1 != 0) { tAuxDataS1 = ArrAuxDataS1 + Offset0; S1_IsDefined = true;}
	if(ArrAuxDataS2 != 0) { tAuxDataS2 = ArrAuxDataS2 + Offset0; S2_IsDefined = true;}
	if(ArrAuxDataS3 != 0) { tAuxDataS3 = ArrAuxDataS3 + Offset0; S3_IsDefined = true;}

	double Inv_ip1 = 1./double(iElecEn + 1);
	double i_d_ip1 = iElecEn*Inv_ip1;
	double SumDifSq = 0, SumNewS0 = 0;

	for(int iz=0; iz<(pStokes->nz); iz++)
	{
		for(int ix=0; ix<(pStokes->nx); ix++)
		{
			for(int ie=0; ie<(pStokes->ne); ie++)
			{
				double NewS0 = (*tBaseSto)*i_d_ip1 + (*(tAuxDataS0++))*Inv_ip1;
				double DifNewOld = NewS0 - (*tBaseSto);
				SumDifSq += DifNewOld*DifNewOld;
				SumNewS0 += NewS0;
				*(tBaseSto++) = (float)NewS0;

				*(tBaseSto++) = (float)(S1_IsDefined? ((*tBaseSto)*i_d_ip1 + (*(tAuxDataS1++))*Inv_ip1) : 0);
				*(tBaseSto++) = (float)(S2_IsDefined? ((*tBaseSto)*i_d_ip1 + (*(tAuxDataS2++))*Inv_ip1) : 0);
				*(tBaseSto++) = (float)(S3_IsDefined? ((*tBaseSto)*i_d_ip1 + (*(tAuxDataS3++))*Inv_ip1) : 0);
			}
		}
		if(S0_IsDefined) tAuxDataS0 += ShiftBetweenRowsX;
		if(S1_IsDefined) tAuxDataS1 += ShiftBetweenRowsX;
		if(S2_IsDefined) tAuxDataS2 += ShiftBetweenRowsX;
		if(S3_IsDefined) tAuxDataS3 += ShiftBetweenRowsX;
	}
	//long TotNp = (pStokes->nz)*(pStokes->nx)*(pStokes->ne);
	long long TotNp = ((long long)(pStokes->nz))*((long long)(pStokes->nx))*((long long)(pStokes->ne));
	double AbsPrecSigma = ::sqrt(SumDifSq/double(TotNp));
	double AvgS0 = SumNewS0/double(TotNp);
	if(AvgS0 == 0) AvgS0 = 1.E-14;
	RelPrec = AbsPrecSigma/AvgS0;

	return RelPrec;
}

//*************************************************************************

void srTRadIntThickBeam::ComputeTotalStokesDistrViaSingleElec(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes)
{
	if(pMagFldTrUnif == 0) throw NO_MAG_FIELD_DEFINED;
    if((pElecBeam == 0) || (pPrcPar == 0)) throw INCORRECT_PARAMS_SR_COMP;
	if(pStokes == 0) throw NO_STOKES_STRUCTURE_SUPPLIED;
	
	srTEbmDat& OrigElecBeam = *pElecBeam;
	srTEbmDat LocElecBeam(OrigElecBeam);
	//CRandGen RandGen;
	CGenMathRand RandGen;

	LocElecBeam.SetNewEnergy(GetNextElecEnergyFromGausDistrib(OrigElecBeam, RandGen));
	srTTrjDat* pTrjDat = (srTTrjDat*)(pMagFldTrUnif->CreateAndSetupNewTrjDat(&LocElecBeam));
	//gAuxPar.Setup(*pElecBeam);
	
	srTWfrSmp* pWfrSmp = 0;
	srTSRWRadStructAccessData* pRad = CreateNewRadStructWithConstParams(pElecBeam, pTrjDat, pStokes, pWfrSmp);
	CHGenObj hRad(pRad);
	srTRadGenManip RadGenManip(hRad);
	srTRadInt RadInt;

	int IntegMethNo = 2; //Auto Wiggler
	double sEndInt = pTrjDat->sStart + (pTrjDat->LenFieldData - 1)*(pTrjDat->sStep);
	srTParPrecElecFld PrecElecFld(IntegMethNo, pPrcPar->RelPrecOrStep, pTrjDat->sStart, sEndInt, 0);

	int ExtractIntType = 1; // Multi-e intensity
	int ExtractPT = 3; // vs x & z
	if(pStokes->ne > 1) ExtractPT = 6;
	double xc = pStokes->xStart + 0.5*(pStokes->nx - 1)*(pStokes->xStep);
	double zc = pStokes->zStart + 0.5*(pStokes->nz - 1)*(pStokes->zStep);
	double ec = pStokes->eStart;

	//long AmOfAuxData = (pWfrSmp->nLamb)*(pWfrSmp->nx)*(pWfrSmp->nz);
	long long AmOfAuxData = ((long long)(pWfrSmp->nLamb))*((long long)(pWfrSmp->nx))*((long long)(pWfrSmp->nz));
	float* ArrAuxDataS0 = new float[AmOfAuxData];
	float* ArrAuxDataS1 = new float[AmOfAuxData];
	float* ArrAuxDataS2 = new float[AmOfAuxData];
	float* ArrAuxDataS3 = new float[AmOfAuxData];

	pStokes->ZeroStokesData();

	double RelPrecElecEn = pPrcPar->RelPrecOrStep;
	double ActualRelPrecIntens = 1E+23;
	//int iElecEn = 0;
	long long iElecEn = 0;

	//int MaxNumIter = 1000000;
	long long MaxNumIter = 1000000;
	if(pPrcPar->NumIter > 0)
	{
        MaxNumIter = pPrcPar->NumIter;
	}

	while(ActualRelPrecIntens > RelPrecElecEn)
	{
		if(ActualRelPrecIntens != 1E+23)
		{
            LocElecBeam.SetNewEnergy(GetNextElecEnergyFromGausDistrib(OrigElecBeam, RandGen));
			if(pTrjDat != 0) { delete pTrjDat; pTrjDat = 0;}
			pTrjDat = (srTTrjDat*)(pMagFldTrUnif->CreateAndSetupNewTrjDat(&LocElecBeam));
			pRad->EmulateElectronBeamStruct(LocElecBeam);
		}

        RadInt.ComputeElectricFieldFreqDomain(pTrjDat, pWfrSmp, &PrecElecFld, pRad);
        RadGenManip.ExtractRadiation(-1, ExtractIntType, ExtractPT, 0, ec, xc, zc, (char*)ArrAuxDataS0);
        RadGenManip.ExtractRadiation(-2, ExtractIntType, ExtractPT, 0, ec, xc, zc, (char*)ArrAuxDataS1);
        RadGenManip.ExtractRadiation(-3, ExtractIntType, ExtractPT, 0, ec, xc, zc, (char*)ArrAuxDataS2);
        RadGenManip.ExtractRadiation(-4, ExtractIntType, ExtractPT, 0, ec, xc, zc, (char*)ArrAuxDataS3);

		ActualRelPrecIntens = UpdateResultStokesData(ArrAuxDataS0, ArrAuxDataS1, ArrAuxDataS2, ArrAuxDataS3, pWfrSmp, iElecEn, pStokes);
		iElecEn++;
		if(iElecEn >= MaxNumIter) break;
	}

	if(pTrjDat != 0) { delete pTrjDat; pTrjDat = 0;}
	if(pWfrSmp != 0) { delete pWfrSmp; pWfrSmp = 0;}
	//if(pRad != 0) { delete pRad; pRad = 0;}

	if(ArrAuxDataS0 != 0) { delete ArrAuxDataS0; ArrAuxDataS0 = 0;}
	if(ArrAuxDataS1 != 0) { delete ArrAuxDataS1; ArrAuxDataS1 = 0;}
	if(ArrAuxDataS2 != 0) { delete ArrAuxDataS2; ArrAuxDataS2 = 0;}
	if(ArrAuxDataS3 != 0) { delete ArrAuxDataS3; ArrAuxDataS3 = 0;}
}

//*************************************************************************
