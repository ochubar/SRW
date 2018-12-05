/************************************************************************//**
 * File: srradstr.cpp
 * Description: Auxiliary structures for various SR calculation methods
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srradstr.h"
#include "srradind.h"
#include "srstraux.h"
#include "srmagfld.h"
#include "srgsnbm.h"
#include "sroptelm.h"
#include "gmfft.h"
#include "gmmeth.h"
#include "gminterp.h"

#ifdef __IGOR_PRO__
#ifndef __SRSEND_H
#include "srsend.h"
#endif
#endif

//*************************************************************************

extern int (*pgWfrExtModifFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp);
extern int (*gpWfrModifFunc)(int action, SRWLWfr* pWfrIn, char pol); //SRWLIB
extern char* (*gpAllocArrayFunc)(char type, long long len); //OC15082018

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTEbmDat* pEbmDat, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{
	if((pEbmDat == 0) || (pTrjDat == 0) || (pWfrSmp == 0)) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
	EmulateElectronBeamStruct(*pEbmDat);
	double From_s0ToObsPoint = pWfrSmp->yStart - pEbmDat->s0;
	yStart = pWfrSmp->yStart;

	AuxSetupActions(pTrjDat, pWfrSmp, From_s0ToObsPoint, NxNzOversamplingFactor);
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTGsnBeam* pGsnBeam, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{
	if((pGsnBeam == 0) || (pWfrSmp == 0)) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
	EmulateElectronBeamStruct(*pGsnBeam);
	double From_s0ToObsPoint = pWfrSmp->yStart - (pGsnBeam->EbmDat).s0;
	yStart = pWfrSmp->yStart;

	AuxSetupActions(0, pWfrSmp, From_s0ToObsPoint, NxNzOversamplingFactor);
}

//*************************************************************************

void srTSRWRadStructAccessData::AuxSetupActions(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, double NxNzOversamplingFactor)
{
	Alloc4x4PropMatr();
    AllocWfrAux();

 	SetRadSamplingFromObs(*pWfrSmp);
	InitialSetupOf4x4PropMatr(From_s0ToObsPoint);
	RobsX = RobsZ = From_s0ToObsPoint;
	yStart = pWfrSmp->yStart;

	if(pTrjDat != 0)
	{
		double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
		int res = 0;
		if(res = FindAverageDistanceToSource(*pTrjDat, *pWfrSmp, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) throw res;
		RobsX = RobsZ = Robs;
		RobsXAbsErr = RobsZAbsErr = RobsAbsErr;
	}
	else
	{
        RobsXAbsErr = RobsZAbsErr = 0.01*From_s0ToObsPoint;
	}


	//if(pMagElem != 0)
	//{
 //       RobsXAbsErr = RobsZAbsErr = 0.5*(pMagElem->GetLongExtent());
	//}
	//else
	//{
 //       RobsXAbsErr = RobsZAbsErr = 0.01*From_s0ToObsPoint;
	//}

	SetupXcZcFromElecData();
	ProcessNxNzForPropag(pWfrSmp, NxNzOversamplingFactor);
	AllocBaseRadAccordingToNeNxNz();
  
	AllocStatMom();
	SetupNonZeroWavefrontLimitsAtCreation();
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTEbmDat* pEbmDat, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{
	if((pRadInData == 0) || (pEbmDat == 0) || (pWfrSmp == 0)) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
	InSRWRadPtrs(pRadInData);
	double From_s0ToObsPoint = pWfrSmp->yStart - pEbmDat->s0;
	yStart = pWfrSmp->yStart;

	AuxSetupActions2(pRadInData, pTrjDat, pWfrSmp, From_s0ToObsPoint, NxNzOversamplingFactor);
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTEbmDat* pEbmDat, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecElecFld* pPrecElecFld)
{
	if((pRadInData == 0) || (pEbmDat == 0) || (pWfrSmp == 0)) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
	InSRWRadPtrs(pRadInData);
	double From_s0ToObsPoint = pWfrSmp->yStart - pEbmDat->s0;
	yStart = pWfrSmp->yStart;

	AuxSetupActions2SR(pRadInData, pTrjDat, pWfrSmp, From_s0ToObsPoint, pPrecElecFld);
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTGsnBeam* pGsnBeam, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{
	if((pRadInData == 0) || (pGsnBeam == 0) || (pWfrSmp == 0)) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
	InSRWRadPtrs(pRadInData);
	double From_s0ToObsPoint = pWfrSmp->yStart - (pGsnBeam->EbmDat).s0;
	yStart = pWfrSmp->yStart;

	AuxSetupActions2(pRadInData, 0, pWfrSmp, From_s0ToObsPoint, NxNzOversamplingFactor);
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTWfrSmp* pWfrSmp, bool AllocateData)
{//used for extracting transmission characteristics of optical elements
	if(pWfrSmp == 0) throw INCORRECT_PARAMS_SR_COMP;

	Initialize();
 	SetRadSamplingFromObs(*pWfrSmp);
	yStart = pWfrSmp->yStart;

	if(AllocateData)
	{
		AllocBaseRadAccordingToNeNxNz('x');
	}

	Pres = 0; // 0- Coord, 1- Ang.
	PresT = 0; // 0- Frequency (Photon Energy), 1- Time
	LengthUnit = 0; // 0- m; 1- mm; 
	PhotEnergyUnit = 0; // 0- eV; 1- keV; 
	ElecFldUnit = 0; // 0- Arb. Units, 1- sqrt(Phot/s/0.1%bw/mm^2)
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(SRWLWfr* pWfr, srTTrjDat* pTrjDat, double* precPar)
{
	if(pWfr == 0) throw SRWL_INCORRECT_WFR_STRUCT;
	
	Initialize();
	InSRWRadPtrs(*pWfr);
	m_newExtWfrCreateNotAllowed = true; //since new wavefront creation is not implemented in SRWLIB

	//if(pTrjDat != 0) AuxSetupActions2SR(*pWfr, *pTrjDat, precPar);
	if(pTrjDat != 0) 
	{
		double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
		int res = 0;
		if(res = FindAverageDistanceToSource(*pTrjDat, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, precPar)) throw res;
		
		//arPrecPar = array('d', [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp])
		double NxNzOversamplingFactor = precPar[6];
		AuxSetupActionsArbSrc(*pWfr, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, NxNzOversamplingFactor);
	}
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(SRWLWfr* pWfr, srTGsnBeam* pGsnBm, double* precPar)
{
	if(pWfr == 0) throw SRWL_INCORRECT_WFR_STRUCT;
	
	Initialize();
	InSRWRadPtrs(*pWfr);
	m_newExtWfrCreateNotAllowed = true; //since new wavefront creation is not implemented in SRWLIB

	if(pGsnBm != 0)
	{
		//yStart = pWfr->zStart;
		yStart = pWfr->mesh.zStart;
		//double Robs = pWfr->zStart - (pGsnBm->EbmDat).s0;
		double Robs = pWfr->mesh.zStart - (pGsnBm->EbmDat).s0;
		double RobsAbsErr = 0.01*(::fabs(Robs));
		double xElAtYsrc = (pGsnBm->EbmDat).x0, zElAtYsrc = (pGsnBm->EbmDat).z0;
		double NxNzOversamplingFactor = precPar[0];
		AuxSetupActionsArbSrc(*pWfr, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, NxNzOversamplingFactor);

		if(pWfr->presFT == 1) //OC171215
		{
			avgPhotEn = pGsnBm->m_AvgPhotEn;
		}
	}
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(SRWLWfr* pWfr, double longPosSrc, double* precPar)
{//Used for setting up spherical wave
	if(pWfr == 0) throw SRWL_INCORRECT_WFR_STRUCT;
	
	Initialize();
	InSRWRadPtrs(*pWfr);
	m_newExtWfrCreateNotAllowed = true; //since new wavefront creation is not implemented in SRWLIB

	yStart = pWfr->mesh.zStart;
	double Robs = pWfr->mesh.zStart - longPosSrc;
	double RobsAbsErr = 0.01*(::fabs(Robs));
	double xElAtYsrc = pWfr->xc, zElAtYsrc = pWfr->yc;
	double NxNzOversamplingFactor = precPar[0];
	AuxSetupActionsArbSrc(*pWfr, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, NxNzOversamplingFactor);

	if(pWfr->presFT == 1)
	{
		avgPhotEn = 0.5*(pWfr->mesh.eStart + pWfr->mesh.eFin);
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::AuxSetupActions2(srTSRWRadInData* pRadInData, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, double NxNzOversamplingFactor)
{
 	SetRadSamplingFromObs(*pWfrSmp);
	InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	RobsX = RobsZ = From_s0ToObsPoint;
	yStart = pWfrSmp->yStart;

	if(pTrjDat != 0)
	{
		double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
		int res = 0;
		if(res = FindAverageDistanceToSource(*pTrjDat, *pWfrSmp, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) throw res;
		RobsX = RobsZ = Robs;
		RobsXAbsErr = RobsZAbsErr = RobsAbsErr;

		xc = xElAtYsrc;
		zc = zElAtYsrc;
	}
	else
	{
        RobsXAbsErr = RobsZAbsErr = 0.01*From_s0ToObsPoint;
        SetupXcZcFromElecData(); //this does not seem correct in general case
	}

	//SetupXcZcFromElecData(); //this does not seem correct in general case

	ProcessNxNzForPropag(pWfrSmp, NxNzOversamplingFactor);

	if(NxNzOversamplingFactor > 0.)
	{
        OutSRWRadPtrs(pRadInData);
		(*pgWfrExtModifFunc)(2, pRadInData, 0);
        InSRWRadPtrs(pRadInData);
	}
	SetupNonZeroWavefrontLimitsAtCreation();
}

//*************************************************************************

void srTSRWRadStructAccessData::AuxSetupActions2SR(srTSRWRadInData* pRadInData, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, srTParPrecElecFld* pPrecElecFld)
{
 	SetRadSamplingFromObs(*pWfrSmp);
	InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	RobsX = RobsZ = From_s0ToObsPoint;
	yStart = pWfrSmp->yStart;

	if(pTrjDat != 0)
	{
		double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
		int res = 0;
		if(res = FindAverageDistanceToSource(*pTrjDat, *pWfrSmp, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, pPrecElecFld)) throw res;
		RobsX = RobsZ = Robs;
		RobsXAbsErr = RobsZAbsErr = RobsAbsErr;

		xc = xElAtYsrc;
		zc = zElAtYsrc;
	}
	else
	{
        RobsXAbsErr = RobsZAbsErr = 0.01*From_s0ToObsPoint;
        SetupXcZcFromElecData(); //this does not seem correct in general case
	}

	ProcessNxNzForPropag(pWfrSmp, pPrecElecFld->NxNzOversamplingFactor);

	if(pPrecElecFld->NxNzOversamplingFactor > 0.)
	{
        OutSRWRadPtrs(pRadInData);
		(*pgWfrExtModifFunc)(2, pRadInData, 0);
        InSRWRadPtrs(pRadInData);
	}
	SetupNonZeroWavefrontLimitsAtCreation();
}

//*************************************************************************

void srTSRWRadStructAccessData::AuxSetupActionsArbSrc(SRWLWfr& srwlWfr, double Robs, double RobsAbsErr, double xElAtYsrc, double zElAtYsrc, double NxNzOversamplingFactor)
{//ATTENTION: this may request for changing numbers of points in the wavefront mesh
 	//SetRadSamplingFromObs(*pWfrSmp);
	//yStart = srwlWfr.zStart;
	yStart = srwlWfr.mesh.zStart;
	//double from_s0ToObsPoint = srwlWfr.zStart - srwlWfr.partBeam.partStatMom1.z;
	double from_s0ToObsPoint = srwlWfr.mesh.zStart - srwlWfr.partBeam.partStatMom1.z;
	InitialSetupOf4x4PropMatr(from_s0ToObsPoint);
	
	//double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	//int res = 0;
	//if(res = FindAverageDistanceToSource(trjDat, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, precPar)) throw res;

	RobsX = RobsZ = Robs;
	RobsXAbsErr = RobsZAbsErr = RobsAbsErr;
	xc = xElAtYsrc;
	zc = zElAtYsrc;

	//to conserve these data in view of eventual resizing:
	srwlWfr.Rx = RobsX; 
	srwlWfr.dRx = RobsXAbsErr;
	srwlWfr.xc = xc;
	srwlWfr.Ry = RobsZ; 
	srwlWfr.dRy = RobsZAbsErr;
	srwlWfr.yc = zc;

	//arPrecPar = array('d', [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp])
	//double NxNzOversamplingFactor = precPar[6];
	if(NxNzOversamplingFactor > 0.)
	{
		long nxNew=-1, nzNew=-1; //?
		UnderSamplingX = UnderSamplingZ = 1.;
		CheckNxNzForSR(NxNzOversamplingFactor, nxNew, nzNew);

		if((nx != nxNew) || (nz != nzNew))
		{
			//srwlWfr.nx = nxNew; srwlWfr.ny = nzNew;
			srwlWfr.mesh.nx = nxNew; srwlWfr.mesh.ny = nzNew;
			//OutSRWRadPtrs(srwlWfr);
			//(*pgWfrExtModifFunc)(2, pRadInData, 0);
			if(gpWfrModifFunc != 0) 
			{//wavefront resizing from external application!
				if((*gpWfrModifFunc)(2, &srwlWfr, 0)) throw SRWL_WFR_EXT_MODIF_FAILED; 
			}
			else throw SRWL_WFR_EXT_FUNC_NOT_DEFINED; 

			InSRWRadPtrs(srwlWfr);
		}
	}
	SetupNonZeroWavefrontLimitsAtCreation();
}

/**
void srTSRWRadStructAccessData::AuxSetupActions2SR(SRWLWfr& srwlWfr, srTTrjDat& trjDat, double* precPar)
{//ATTENTION: this may request for changing numbers of points in the wavefront mesh
 	//SetRadSamplingFromObs(*pWfrSmp);
	yStart = srwlWfr.zStart;
	double from_s0ToObsPoint = srwlWfr.zStart - srwlWfr.partBeam.partStatMom1.z;
	InitialSetupOf4x4PropMatr(from_s0ToObsPoint);
	
	double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	int res = 0;
	if(res = FindAverageDistanceToSource(trjDat, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc, precPar)) throw res;

	RobsX = RobsZ = Robs;
	RobsXAbsErr = RobsZAbsErr = RobsAbsErr;
	xc = xElAtYsrc;
	zc = zElAtYsrc;

	//to conserve these data in view of eventual resizing:
	srwlWfr.Rx = RobsX; 
	srwlWfr.dRx = RobsXAbsErr;
	srwlWfr.xc = xc;
	srwlWfr.Ry = RobsZ; 
	srwlWfr.dRy = RobsZAbsErr;
	srwlWfr.yc = zc;

	//arPrecPar = array('d', [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp])
	double NxNzOversamplingFactor = precPar[6];
	if(NxNzOversamplingFactor > 0.)
	{
		long nxNew=-1, nzNew=-1; //?
		UnderSamplingX = UnderSamplingZ = 1.;
		CheckNxNzForSR(NxNzOversamplingFactor, nxNew, nzNew);

		if((nx != nxNew) || (nz != nzNew))
		{
			srwlWfr.nx = nxNew; srwlWfr.ny = nzNew;
			//OutSRWRadPtrs(srwlWfr);
			//(*pgWfrExtModifFunc)(2, pRadInData, 0);
			if(gpWfrModifFunc != 0) 
			{//wavefront resizing from external application!
				if((*gpWfrModifFunc)(2, &srwlWfr, 0)) throw SRWL_WFR_EXT_MODIF_FAILED; 
			}
			else throw SRWL_WFR_EXT_FUNC_NOT_DEFINED; 

			InSRWRadPtrs(srwlWfr);
		}
	}
	SetupNonZeroWavefrontLimitsAtCreation();
}
**/
//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, int CopyData)
{
	if(pRadInData == 0) throw INCORRECT_ARGUMENTS;

	Initialize();
	bool DataShouldBeCopied = false;
	if(CopyData != 0) DataShouldBeCopied = true;
	InSRWRadPtrs(pRadInData, DataShouldBeCopied);
}

//*************************************************************************

void srTSRWRadStructAccessData::InSRWRadPtrs(srTSRWRadInData* p, bool DataShouldBeCopied)
{
	if(p == 0) throw INCORRECT_PARAMS_SR_COMP;

	pBaseRadX = p->pBaseRadX; pBaseRadZ = p->pBaseRadZ;
	wRad = p->wRad; wRadX = p->wRadX; wRadZ = p->wRadZ;
	hStateRadX = p->hStateRadX; hStateRadZ = p->hStateRadZ;
	eStep = p->eStep; eStart = p->eStart; 
	xStep = p->xStep; xStart = p->xStart; 
	zStep = p->zStep; zStart = p->zStart;
	ne = p->ne; nx = p->nx; nz = p->nz;

	RobsX = p->RobsX; RobsZ = p->RobsZ;
	RobsXAbsErr = p->RobsXAbsErr; RobsZAbsErr = p->RobsZAbsErr;
	xc = p->xc; zc = p->zc;
	xWfrMin = p->xWfrMin; xWfrMax = p->xWfrMax;
	zWfrMin = p->zWfrMin; zWfrMax = p->zWfrMax;

	UnderSamplingX = p->UnderSamplingX; UnderSamplingZ = p->UnderSamplingZ;
	AllowAutoSwitchToPropInUnderSamplingMode = p->AllowAutoSwitchToPropInUnderSamplingMode;
	InvUnderSamplingThreshold = p->InvUnderSamplingThreshold;

	Pres = p->Pres;
	PresT = p->PresT;
	LengthUnit = p->LengthUnit;
	PhotEnergyUnit = p->PhotEnergyUnit;
	ElecFldUnit = p->ElecFldUnit;
	avgPhotEn = p->avgPhotEn;

	pElecBeam = p->pElecBeam;
	wElecBeam = p->wElecBeam;
	hStateElecBeam = p->hStateElecBeam;

	wTrj = p->wTrj; // Can be this or Electron Beam
	hStateTrj = p->hStateTrj;

	p4x4PropMatr = p->p4x4PropMatr;
	w4x4PropMatr = p->w4x4PropMatr;
	hState4x4PropMatr = p->hState4x4PropMatr;

	pMomX = p->pMomX; pMomZ = p->pMomZ;
	wMomX = p->wMomX; wMomZ = p->wMomZ;
	hStateMomX = p->hStateMomX; hStateMomZ = p->hStateMomZ;

	pWfrAuxData = p->pWfrAuxData;
	wWfrAuxData = p->wWfrAuxData;
	hStateWfrAuxData = p->hStateWfrAuxData;

	if(DataShouldBeCopied)
	{
		AllocBaseRadAccordingToNeNxNz(); CopyBaseRadData(p->pBaseRadX, p->pBaseRadZ);
		AllocStatMom(); CopyStatMomData(p->pMomX, p->pMomZ);
		AllocElectronBeam(); CopyElectronBeamData(p->pElecBeam);
		Alloc4x4PropMatr(); Copy4x4PropMatrData(p->p4x4PropMatr);
        AllocWfrAux(); CopyWfrAuxData(p->pWfrAuxData);
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::InSRWRadPtrs(SRWLWfr& srwlWfr)
{
	pBaseRadX = (float*)srwlWfr.arEx; pBaseRadZ = (float*)srwlWfr.arEy;
	pBaseRadXaux = (float*)srwlWfr.arExAux; pBaseRadZaux = (float*)srwlWfr.arEyAux; //OC151115
	wRad = 0; wRadX = 0; wRadZ = 0;
	hStateRadX = 0; hStateRadZ = 0;

	//eStart = srwlWfr.eStart;
	//eStep = (srwlWfr.ne <= 1)? 0 : (srwlWfr.eFin - srwlWfr.eStart)/(srwlWfr.ne - 1); 
	//xStart = srwlWfr.xStart; 
	//xStep = (srwlWfr.nx <= 1)? 0 : (srwlWfr.xFin - srwlWfr.xStart)/(srwlWfr.nx - 1); 
	//zStart = srwlWfr.yStart;
	//zStep = (srwlWfr.ny <= 1)? 0 : (srwlWfr.yFin - srwlWfr.yStart)/(srwlWfr.ny - 1); 
	//ne = srwlWfr.ne; nx = srwlWfr.nx; nz = srwlWfr.ny;
	//yStart = srwlWfr.zStart; //OC21092011

	SRWLStructRadMesh &mesh = srwlWfr.mesh;
	eStart = mesh.eStart;
	eStep = (mesh.ne <= 1)? 0 : (mesh.eFin - mesh.eStart)/(mesh.ne - 1); 

	//xStart = mesh.xStart; 
	//xStep = (mesh.nx <= 1)? 0 : (mesh.xFin - mesh.xStart)/(mesh.nx - 1); 
	if(mesh.nx <= 1) //OC170615
	{
		xStart = 0.5*(mesh.xStart + mesh.xFin);
		xStep = 0.;
	}
	else
	{
		xStart = mesh.xStart; 
		xStep = (mesh.xFin - mesh.xStart)/(mesh.nx - 1); 
	}

	//zStart = mesh.yStart;
	//zStep = (mesh.ny <= 1)? 0 : (mesh.yFin - mesh.yStart)/(mesh.ny - 1); 
	if(mesh.ny <= 1) //OC170615
	{
		zStart = 0.5*(mesh.yStart + mesh.yFin);
		zStep = 0.;
	}
	else
	{
		zStart = mesh.yStart;
		zStep = (mesh.yFin - mesh.yStart)/(mesh.ny - 1); 
	}

	ne = mesh.ne; nx = mesh.nx; nz = mesh.ny;
	yStart = mesh.zStart; //OC21092011

	RobsX = srwlWfr.Rx; RobsZ = srwlWfr.Ry;
	RobsXAbsErr = srwlWfr.dRx; RobsZAbsErr = srwlWfr.dRy;
	xc = srwlWfr.xc; zc = srwlWfr.yc;

	//xWfrMin = srwlWfr.xStart; xWfrMax = srwlWfr.xFin;
	//zWfrMin = srwlWfr.yStart; zWfrMax = srwlWfr.yFin;
	xWfrMin = mesh.xStart; xWfrMax = mesh.xFin;
	zWfrMin = mesh.yStart; zWfrMax = mesh.yFin;

	UnderSamplingX = 1.; UnderSamplingZ = 1.;
	AllowAutoSwitchToPropInUnderSamplingMode = 0; //?
	InvUnderSamplingThreshold = 0.01; //?

	Pres = srwlWfr.presCA;
	PresT = srwlWfr.presFT;
	ElecFldUnit = srwlWfr.unitElFld;
	ElecFldAngUnit = srwlWfr.unitElFldAng; //OC20112017

	avgPhotEn = srwlWfr.avgPhotEn;
	LengthUnit = 0; // 0- m; 1- mm; 
	PhotEnergyUnit = 0; // 0- eV; 1- keV; 

	avgT = (PresT == 1)? (eStart + 0.5*eStep*(ne - 1)) : 0; //OC101015 //???

	EmulateElectronBeamStruct(srwlWfr.partBeam);
	wElecBeam = 0;
	hStateElecBeam = 0;

	wTrj = 0; //?
	hStateTrj = 0;

	p4x4PropMatr = srwlWfr.arElecPropMatr;
	w4x4PropMatr = 0;
	hState4x4PropMatr = 0;

	pWfrAuxData = srwlWfr.arWfrAuxData;
	wWfrAuxData = 0;
	hStateWfrAuxData = 0;

	wMomX = 0; wMomZ = 0;
	hStateMomX = 0; hStateMomZ = 0;
	pMomX = srwlWfr.arMomX; //OC130311
	pMomZ = srwlWfr.arMomY;

	m_pExtWfr = &srwlWfr; //OC130311

	//int LenMomComp = ne*11; // to steer
	//if(srwlWfr.arMomX != 0)
	//{
	//	if(MomWereEmulated) { if(pMomX != 0) delete[] pMomX;}
	//	//pMomX = new float[LenMomComp];
	//	//float *tMomX = pMomX;
	//	pMomX = new double[LenMomComp]; //OC130311
	//	double *tMomX = pMomX;

	//	double *tInMomX = srwlWfr.arMomX;
	//	//for(int i=0; i<LenMomComp; i++) *(tMomX++) = (float)(*(tInMomX++));
	//	for(int i=0; i<LenMomComp; i++) *(tMomX++) = *(tInMomX++); //OC130311
	//	MomWereEmulated = true;
	//}
	//if(srwlWfr.arMomY != 0)
	//{
	//	if(MomWereEmulated) { if(pMomZ != 0) delete[] pMomZ;}
	//	//pMomZ = new float[LenMomComp];
	//	//float *tMomZ = pMomZ;
	//	pMomZ = new double[LenMomComp];
	//	double *tMomZ = pMomZ;

	//	double *tInMomY = srwlWfr.arMomY;
	//	//for(int i=0; i<LenMomComp; i++) *(tMomZ++) = (float)(*(tInMomY++));
	//	for(int i=0; i<LenMomComp; i++) *(tMomZ++) = *(tInMomY++); //OC130311
	//	MomWereEmulated = true;
	//}

	//if(DataShouldBeCopied)
	//{
	//	AllocBaseRadAccordingToNeNxNz(); CopyBaseRadData(p->pBaseRadX, p->pBaseRadZ);
	//	AllocStatMom(); CopyStatMomData(p->pMomX, p->pMomZ);
	//	AllocElectronBeam(); CopyElectronBeamData(p->pElecBeam);
	//	Alloc4x4PropMatr(); Copy4x4PropMatrData(p->p4x4PropMatr);
	//	AllocWfrAux(); CopyWfrAuxData(p->pWfrAuxData);
	//}
}

//*************************************************************************

void srTSRWRadStructAccessData::OutSRWRadPtrs(srTSRWRadInData* p)
{
	if(p == 0) throw INCORRECT_PARAMS_SR_COMP;

	p->pBaseRadX = pBaseRadX; p->pBaseRadZ = pBaseRadZ;
	p->wRad = wRad; p->wRadX = wRadX; p->wRadZ = wRadZ;
	p->hStateRadX = hStateRadX; p->hStateRadZ = hStateRadZ;
	p->eStep = eStep; p->eStart = eStart; 
	p->xStep = xStep; p->xStart = xStart; 
	p->zStep = zStep; p->zStart = zStart;
	p->ne = ne; p->nx = nx; p->nz = nz;

	p->RobsX = RobsX; p->RobsZ = RobsZ;
	p->RobsXAbsErr = RobsXAbsErr; p->RobsZAbsErr = RobsZAbsErr;
	p->xc = xc; p->zc = zc;
	p->xWfrMin = xWfrMin; p->xWfrMax = xWfrMax;
	p->zWfrMin = zWfrMin; p->zWfrMax = zWfrMax;

	p->UnderSamplingX = UnderSamplingX; p->UnderSamplingZ = UnderSamplingZ;
	p->AllowAutoSwitchToPropInUnderSamplingMode = AllowAutoSwitchToPropInUnderSamplingMode;
	p->InvUnderSamplingThreshold = InvUnderSamplingThreshold;

	p->Pres = Pres;
	p->PresT = PresT;
	p->LengthUnit = LengthUnit;
	p->PhotEnergyUnit = PhotEnergyUnit;
	p->ElecFldUnit = ElecFldUnit;
	p->avgPhotEn = avgPhotEn;

	p->pElecBeam = pElecBeam;
	p->wElecBeam = wElecBeam;
	p->hStateElecBeam = hStateElecBeam;

	p->wTrj = wTrj; // Can be this or Electron Beam
	p->hStateTrj = hStateTrj;

	p->p4x4PropMatr = p4x4PropMatr;
	p->w4x4PropMatr = w4x4PropMatr;
	p->hState4x4PropMatr = hState4x4PropMatr;

	p->pMomX = pMomX; p->pMomZ = pMomZ;
	p->wMomX = wMomX; p->wMomZ = wMomZ;
	p->hStateMomX = hStateMomX; p->hStateMomZ = hStateMomZ;

	p->pWfrAuxData = pWfrAuxData;
	p->wWfrAuxData = wWfrAuxData;
	p->hStateWfrAuxData = hStateWfrAuxData;
}

//*************************************************************************

void srTSRWRadStructAccessData::OutSRWRadPtrs(SRWLWfr& srwlWfr)
{
	srwlWfr.arEx = (char*)pBaseRadX; srwlWfr.arEy = (char*)pBaseRadZ;
	srwlWfr.arExAux = (char*)pBaseRadXaux; srwlWfr.arEyAux = (char*)pBaseRadZaux; //OC151115
	//p->wRad = wRad; p->wRadX = wRadX; p->wRadZ = wRadZ;
	//p->hStateRadX = hStateRadX; p->hStateRadZ = hStateRadZ;
	
	//srwlWfr.eStart = eStart; srwlWfr.eFin = eStart + eStep*(ne - 1); 
	//srwlWfr.xStart = xStart; srwlWfr.xFin = xStart + xStep*(nx - 1); 
	//srwlWfr.yStart = zStart; srwlWfr.yFin = zStart + zStep*(nz - 1); 
	//srwlWfr.ne = ne; srwlWfr.nx = nx; srwlWfr.ny = nz;

	SRWLStructRadMesh &mesh = srwlWfr.mesh;
	mesh.eStart = eStart; mesh.eFin = eStart + eStep*(ne - 1); 
	mesh.xStart = xStart; mesh.xFin = xStart + xStep*(nx - 1); 
	mesh.yStart = zStart; mesh.yFin = zStart + zStep*(nz - 1); 
	mesh.ne = ne; mesh.nx = nx; mesh.ny = nz;

	srwlWfr.Rx = RobsX; srwlWfr.Ry = RobsZ;
	srwlWfr.dRx = RobsXAbsErr; srwlWfr.dRy = RobsZAbsErr;
	srwlWfr.xc = xc; srwlWfr.yc = zc;

	//p->xWfrMin = xWfrMin; p->xWfrMax = xWfrMax;
	//p->zWfrMin = zWfrMin; p->zWfrMax = zWfrMax;
	//p->UnderSamplingX = UnderSamplingX; p->UnderSamplingZ = UnderSamplingZ;
	//p->AllowAutoSwitchToPropInUnderSamplingMode = AllowAutoSwitchToPropInUnderSamplingMode;
	//p->InvUnderSamplingThreshold = InvUnderSamplingThreshold;

	srwlWfr.presCA = Pres;
	srwlWfr.presFT = PresT;

	//p->LengthUnit = LengthUnit;
	//p->PhotEnergyUnit = PhotEnergyUnit;

	srwlWfr.unitElFld = ElecFldUnit;
	srwlWfr.avgPhotEn = avgPhotEn;
	
	OutElectronBeamStruct(srwlWfr.partBeam);
	//p->wElecBeam = wElecBeam;
	//p->hStateElecBeam = hStateElecBeam;
	//p->wTrj = wTrj; // Can be this or Electron Beam
	//p->hStateTrj = hStateTrj;

	srwlWfr.arElecPropMatr = p4x4PropMatr;
	//p->w4x4PropMatr = w4x4PropMatr;
	//p->hState4x4PropMatr = hState4x4PropMatr;

	//p->pMomX = pMomX; p->pMomZ = pMomZ;
	//p->wMomX = wMomX; p->wMomZ = wMomZ;
	//p->hStateMomX = hStateMomX; p->hStateMomZ = hStateMomZ;

	//int LenMomComp = ne*11; // to steer
	//float *tMomX = pMomX, *tMomZ = pMomZ;
	//double *tOutMomX = srwlWfr.arMomX, *tOutMomY = srwlWfr.arMomY;
	//bool momXareDefined = (srwlWfr.arMomX != 0);
	//bool momYareDefined = (srwlWfr.arMomY != 0);
	//for(int i=0; i<LenMomComp; i++) 
	//{
	//	if(momXareDefined) *(tOutMomX++) = *(tMomX++); 
	//	if(momYareDefined) *(tOutMomY++) = *(tMomZ++);
	//}
	//OC130311
	srwlWfr.arMomX = pMomX;
	srwlWfr.arMomY = pMomZ;

	srwlWfr.arWfrAuxData = pWfrAuxData;
	//p->wWfrAuxData = wWfrAuxData;
	//p->hStateWfrAuxData = hStateWfrAuxData;
}

//*************************************************************************

//int srTSRWRadStructAccessData::ModifyWfrNeNxNz(char PolarizComp)
int srTSRWRadStructAccessData::ModifyWfrNeNxNz(char PolarizComp, bool backupIsReq) 
{//OC131115
#if defined(_SRWDLL) || defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	int res = 0;
	if(BaseRadWasEmulated) return ReAllocBaseRadAccordingToNeNxNz(PolarizComp);
	
	if(pgWfrExtModifFunc != 0)
	{//to be removed!!
		srTSRWRadInData AuxRadInData;
		OutSRWRadPtrs(&AuxRadInData);
		//if((*pgWfrExtModifFunc)(2, &AuxRadInData, PolarizComp)) return SRWL_WFR_EXT_MODIF_FAILED;
		if((*pgWfrExtModifFunc)(2, &AuxRadInData, PolarizComp)) return SRWL_WFR_EXT_MODIF_FAILED;

		InSRWRadPtrs(&AuxRadInData);
	}
	else if((gpWfrModifFunc != 0) && (m_pExtWfr != 0))
	{
		SRWLWfr *pExtWfr = (SRWLWfr*)m_pExtWfr;
		OutSRWRadPtrs(*pExtWfr);

		int actNum = 2;
		if(backupIsReq) actNum = 12; //OC131115
		if((*gpWfrModifFunc)(actNum, pExtWfr, PolarizComp)) return SRWL_WFR_EXT_MODIF_FAILED;

		InSRWRadPtrs(*pExtWfr);
	}
	else if(gpWfrModifFunc == 0) return SRWL_WFR_EXT_FUNC_NOT_DEFINED;
	else if(pgWfrExtModifFunc == 0) return SRWL_WFR_EXT_FUNC_NOT_DEFINED;

	return 0;
#endif
#ifdef __IGOR_PRO__
	srTSend Send;
	return Send.ModifyRadNeNxNz(*this, PolarizComp);
#endif
}

//*************************************************************************

int srTSRWRadStructAccessData::AllocExtIntArray(char type, char dep, char*& pcAlloc) 
{//OC18082018
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	
	pcAlloc = 0;
	if(gpAllocArrayFunc != 0)
	{
		char typeAr = 'f';
		if(type == 4) typeAr = 'd'; //single-e rad. phase
	
		long long np = GetIntNumPts(dep);
		if(np > 0)
		{
			pcAlloc = (*gpAllocArrayFunc)(typeAr, np);
			if(pcAlloc == 0) return SRWL_EXT_ARRAY_ALLOC_FAILED;
		}
	}

#endif
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::DeleteWfrBackupData(char PolarizComp)
{//OC131115
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 

	if((gpWfrModifFunc != 0) && (m_pExtWfr != 0))
	{
		SRWLWfr *pExtWfr = (SRWLWfr*)m_pExtWfr;
		OutSRWRadPtrs(*pExtWfr);

		int actNum = 20; //Delete backup data
		if((*gpWfrModifFunc)(actNum, pExtWfr, PolarizComp)) return SRWL_WFR_EXT_MODIF_FAILED;

		InSRWRadPtrs(*pExtWfr);
	}

#endif
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::GetWfrStructNames(srTSRWRadStructWaveNames& Names)
{
//#ifdef _SRWDLL
#if defined(_SRWDLL) || defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	if(BaseRadWasEmulated) return 0;
	if(pgWfrExtModifFunc == 0) return SRWL_WFR_EXT_FUNC_NOT_DEFINED;

	int res = 0;
	srTSRWRadInData AuxRadInData;
	OutSRWRadPtrs(&AuxRadInData);
	if((*pgWfrExtModifFunc)(4, &AuxRadInData, 0)) return SRWL_WFR_EXT_MODIF_FAILED;

	strcpy(Names.NameRad, AuxRadInData.NameRad);
	strcpy(Names.NameRadX, AuxRadInData.NameRadX); 
	strcpy(Names.NameRadZ, AuxRadInData.NameRadZ);
	strcpy(Names.NameElecBeam, AuxRadInData.NameElecBeam);
	strcpy(Names.NameTrj, AuxRadInData.NameTrj);
	strcpy(Names.Name4x4PropMatr, AuxRadInData.Name4x4PropMatr);
	strcpy(Names.NameMomX, AuxRadInData.NameMomX); 
	strcpy(Names.NameMomZ, AuxRadInData.NameMomZ);
	strcpy(Names.NameWfrAuxData, AuxRadInData.NameWfrAuxData);
	return 0;
#endif
#ifdef __IGOR_PRO__
	srTSend Send;
    return Send.GetRadStructNames(*this, Names);
#endif
}

//*************************************************************************

int srTSRWRadStructAccessData::DeleteWfrStructWaves(srTSRWRadStructWaveKeys& RadKeys)
{
//#ifdef _SRWDLL
#if defined(_SRWDLL) || defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	if(BaseRadWasEmulated) return 0;
	if(pgWfrExtModifFunc == 0) return SRWL_WFR_EXT_FUNC_NOT_DEFINED;

	int res = 0;
	srTSRWRadInData AuxRadInData;
	OutSRWRadPtrs(&AuxRadInData);

	AuxRadInData.wRad_ = RadKeys.wRad_;
	AuxRadInData.wRadX_ = RadKeys.wRadX_;
	AuxRadInData.wRadZ_ = RadKeys.wRadZ_;
	AuxRadInData.wElecBeam_ = RadKeys.wElecBeam_;
	AuxRadInData.wTrj_ = RadKeys.wTrj_;
	AuxRadInData.w4x4PropMatr_ = RadKeys.w4x4PropMatr_;
	AuxRadInData.wMomX_ = RadKeys.wMomX_;
	AuxRadInData.wMomZ_ = RadKeys.wMomZ_;
	AuxRadInData.wWfrAuxData_ = RadKeys.wWfrAuxData_;

	if((*pgWfrExtModifFunc)(0, &AuxRadInData, 0)) return SRWL_WFR_EXT_MODIF_FAILED;
    
	InSRWRadPtrs(&AuxRadInData);
	return 0;
#endif
#ifdef __IGOR_PRO__
	srTSend Send;
	return Send.DeleteRadStructWaves(*this, RadKeys);
#endif
}

//*************************************************************************

int srTSRWRadStructAccessData::RenameWfrStruct(srTSRWRadStructWaveNames& Names)
{
//#ifdef _SRWDLL
#if defined(_SRWDLL) || defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	if(BaseRadWasEmulated) return 0;
	if(pgWfrExtModifFunc == 0) return SRWL_WFR_EXT_FUNC_NOT_DEFINED;

	int res = 0;
	srTSRWRadInData AuxRadInData;
	OutSRWRadPtrs(&AuxRadInData);
	if(res = (*pgWfrExtModifFunc)(4, &AuxRadInData, 0)) return SRWL_WFR_EXT_MODIF_FAILED;

	strcpy(AuxRadInData.NameRad, Names.NameRad);
	strcpy(AuxRadInData.NameRadX, Names.NameRadX); 
	strcpy(AuxRadInData.NameRadZ, Names.NameRadZ);
	strcpy(AuxRadInData.NameElecBeam, Names.NameElecBeam);
	strcpy(AuxRadInData.NameTrj, Names.NameTrj);
	strcpy(AuxRadInData.Name4x4PropMatr, Names.Name4x4PropMatr);
	strcpy(AuxRadInData.NameMomX, Names.NameMomX); 
	strcpy(AuxRadInData.NameMomZ, Names.NameMomZ);
	strcpy(AuxRadInData.NameWfrAuxData, Names.NameWfrAuxData);

	if((*pgWfrExtModifFunc)(3, &AuxRadInData, 0)) return SRWL_WFR_EXT_MODIF_FAILED;

    InSRWRadPtrs(&AuxRadInData);
	return 0;
#endif
#ifdef __IGOR_PRO__
	srTSend Send;
	return Send.RenameRadStruct(*this, Names);
#endif
}

//*************************************************************************

int srTSRWRadStructAccessData::CreateNewWfrStruct(srTSRWRadStructWaveNames& Names)
{
#if defined(_SRWDLL) || defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) 
	if(BaseRadWasEmulated) return 0;
	int res = 0;

	if(pgWfrExtModifFunc != 0)
	{
		srTSRWRadInData AuxRadInData;
		OutSRWRadPtrs(&AuxRadInData);

		strcpy(AuxRadInData.NameRad, Names.NameRad);
		strcpy(AuxRadInData.NameRadX, Names.NameRadX); 
		strcpy(AuxRadInData.NameRadZ, Names.NameRadZ);
		strcpy(AuxRadInData.NameElecBeam, Names.NameElecBeam);
		strcpy(AuxRadInData.NameTrj, Names.NameTrj);
		strcpy(AuxRadInData.Name4x4PropMatr, Names.Name4x4PropMatr);
		strcpy(AuxRadInData.NameMomX, Names.NameMomX); 
		strcpy(AuxRadInData.NameMomZ, Names.NameMomZ);
		strcpy(AuxRadInData.NameWfrAuxData, Names.NameWfrAuxData);

		//if(res = (*pgWfrExtModifFunc)(1, &AuxRadInData, 0)) return res;
		if((*pgWfrExtModifFunc)(1, &AuxRadInData, 0)) return SRWL_WFR_EXT_MODIF_FAILED;

		InSRWRadPtrs(&AuxRadInData);
	}
	else return SRWL_WFR_EXT_FUNC_NOT_DEFINED;
	//SRWLIB version not implemented: all modifications go through srTSRWRadStructAccessData::ModifyWfrNeNxNz(char PolarizComp)

	//else if(gpWfrModifFunc != 0)
	//{
	//	SRWLWfr auxWfr;
	//	OutSRWRadPtrs(auxWfr);
	//	if(res = (*gpWfrModifFunc)(1, &auxWfr, 0)) return res;
	//	InSRWRadPtrs(auxWfr);
	//}
	return 0;
#endif
#ifdef __IGOR_PRO__
	srTSend Send;
	return Send.CreateNewRadStruct(*this, Names);
#endif
}

//*************************************************************************

void srTSRWRadStructAccessData::ProcessNxNzForPropag(srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{// Assumes that Robs, xc, zc were already set !
	if(NxNzOversamplingFactor <= 0.) return;

	UnderSamplingX = UnderSamplingZ = 1;

	pWfrSmp->nx = pWfrSmp->nz = -1;
	CheckNxNzForSR(pWfrSmp, NxNzOversamplingFactor);

	if(pWfrSmp->DimensionsWereSetAuto)
	{
        UpdateObsParam(*pWfrSmp);
		//if(result = Send.ModifyRadNeNxNz(SRWRadStructAccessData)) return result;
		pWfrSmp->DimensionsWereSetAuto = 0;
	}
	pWfrSmp->AllowAutoChoiceOfNxNzForPropagat = 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::CheckNxNzForSR(srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor)
{// Assumes that Robs, xc, zc were already set !
	long Nx = pWfrSmp->nx, Nz = pWfrSmp->nz;
	if((Nx > 0) && (Nz > 0)) return;
	if(NxNzOversamplingFactor <= 0.) return;

	const int SmallestN = 8;
	//double WavelengthIn_m = (pWfrSmp->TreatLambdaAsEnergyIn_eV)? 1.239854E-06/(pWfrSmp->LambStart) : 1.E-06*(pWfrSmp->LambEnd);
	double WavelengthIn_m = (pWfrSmp->TreatLambdaAsEnergyIn_eV)? 1.239842E-06/(pWfrSmp->LambStart) : 1.E-06*(pWfrSmp->LambEnd);

	CGenMathFFT2D FFT;
	if(pWfrSmp->CoordOrAngPresentation == CoordPres)
	{
		if(Nx < 1)
		{
            double HalfLambR = 0.5*WavelengthIn_m*RobsX;
			double xStartRel = pWfrSmp->xStart - xc;
			double xEndRel = pWfrSmp->xEnd - xc;
			double dxStart = ::fabs(HalfLambR/xStartRel);
			double dxEnd = ::fabs(HalfLambR/xEndRel);
			double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/NxNzOversamplingFactor;
			dx /= 1.2;

			Nx = long(::fabs(xEndRel - xStartRel)/dx) + 1;
			if(((Nx>>1)<<1) != Nx) Nx++;
			FFT.NextCorrectNumberForFFT(Nx);
			if(Nx < SmallestN) Nx = SmallestN;
		}
		if(Nz < 1)
		{
            double HalfLambR = 0.5*WavelengthIn_m*RobsZ;
			double zStartRel = pWfrSmp->zStart - zc;
			double zEndRel = pWfrSmp->zEnd - zc;
			double dzStart = ::fabs(HalfLambR/zStartRel);
			double dzEnd = ::fabs(HalfLambR/zEndRel);
			double dz = ((dzStart < dzEnd)? dzStart : dzEnd)/NxNzOversamplingFactor;
			dz /= 1.2;

			Nz = long(::fabs(zEndRel - zStartRel)/dz) + 1;
			if(((Nz>>1)<<1) != Nz) Nz++;
			FFT.NextCorrectNumberForFFT(Nz);
			if(Nz < SmallestN) Nz = SmallestN;
		}
	}
// Continue here if(DistrInfoDat.CoordOrAngPresentation == AngPres)

	pWfrSmp->nx = long(Nx); pWfrSmp->nz = long(Nz);
	pWfrSmp->DimensionsWereSetAuto = 1;
}

//*************************************************************************

void srTSRWRadStructAccessData::CheckNxNzForSR(double NxNzOversamplingFactor, long& _nx, long& _nz)
{// Assumes that Robs, xc, zc and all other main members were already set in this!

	//long Nx = pWfrSmp->nx, Nz = pWfrSmp->nz;
	if((_nx > 0) && (_nz > 0)) return;
	if(NxNzOversamplingFactor <= 0.) return;

	const int SmallestN = 8;
	//double WavelengthIn_m = 1.239842E-06/eStart;
	double ePhRef = (PresT == 0)? eStart : avgPhotEn; //OC041115 (making sure this works in time domain as well)
	double WavelengthIn_m = 1.239842E-06/ePhRef; //OC041115
	//(pWfrSmp->TreatLambdaAsEnergyIn_eV)? 1.239842E-06/(pWfrSmp->LambStart) : 1.E-06*(pWfrSmp->LambEnd);

	CGenMathFFT2D FFT;
	//if(pWfrSmp->CoordOrAngPresentation == CoordPres)
	if(Pres == 0)
	{//Coordinate presentation
		if(_nx < 1)
		{
            double HalfLambR = 0.5*WavelengthIn_m*RobsX;
			//double xStartRel = pWfrSmp->xStart - xc;
			//double xEndRel = pWfrSmp->xEnd - xc;
			double xStartRel = xStart - xc;
			double xEnd = xStart + (nx - 1)*xStep; //?
			double xEndRel = xEnd - xc;
			double dxStart = ::fabs(HalfLambR/xStartRel);
			double dxEnd = ::fabs(HalfLambR/xEndRel);
			double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/(NxNzOversamplingFactor*1.2);

			_nx = long(::fabs(xEndRel - xStartRel)/dx) + 1;
			if(((_nx>>1)<<1) != _nx) _nx++;
			FFT.NextCorrectNumberForFFT(_nx);
			if(_nx < SmallestN) _nx = SmallestN;
		}
		if(_nz < 1)
		{
            double HalfLambR = 0.5*WavelengthIn_m*RobsZ;
			//double zStartRel = pWfrSmp->zStart - zc;
			//double zEndRel = pWfrSmp->zEnd - zc;
			double zStartRel = zStart - zc;
			double zEnd = zStart + (nz - 1)*zStep; //?
			double zEndRel = zEnd - zc;
			double dzStart = ::fabs(HalfLambR/zStartRel);
			double dzEnd = ::fabs(HalfLambR/zEndRel);
			double dz = ((dzStart < dzEnd)? dzStart : dzEnd)/(NxNzOversamplingFactor*1.2);

			_nz = long(::fabs(zEndRel - zStartRel)/dz) + 1;
			if(((_nz>>1)<<1) != _nz) _nz++;
			FFT.NextCorrectNumberForFFT(_nz);
			if(_nz < SmallestN) _nz = SmallestN;
		}
	}
// Continue here if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	//pWfrSmp->nx = long(Nx); pWfrSmp->nz = long(Nz);
	//pWfrSmp->DimensionsWereSetAuto = 1;
}

//*************************************************************************

void srTSRWRadStructAccessData::EstimateOversamplingFactors(double& estimOverSampX, double& estimOverSampZ)
{
	//const int SmallestN = 8;
	//double WavelengthIn_m = (pWfrSmp->TreatLambdaAsEnergyIn_eV)? 1.239854E-06/(pWfrSmp->LambStart) : 1.E-06*(pWfrSmp->LambEnd);
	
	//double WavelengthIn_m = 1.239854E-06/eStart;
	double WavelengthIn_m = 1.239842E-06/eStart;
	CGenMathFFT2D FFT;

	if(Pres == 0) //Coord.
	{
		double HalfLambRx = 0.5*WavelengthIn_m*RobsX;
		double xStartRel = xStart - xc;
		double xEndAbs = xStart + (nx - 1)*xStep;
		double xEndRel = xEndAbs - xc;
		double dxStart = ::fabs(HalfLambRx/xStartRel);
		double dxEnd = ::fabs(HalfLambRx/xEndRel);
		double dx = ((dxStart < dxEnd)? dxStart : dxEnd); // /NxNzOversamplingFactor;
		dx /= 1.2;
		long Nx_nom = long(::fabs(xEndRel - xStartRel)/dx) + 1;
		if(((Nx_nom>>1)<<1) != Nx_nom) Nx_nom++;
		FFT.NextCorrectNumberForFFT(Nx_nom);
		//if(Nx_nom < SmallestN) Nx_nom = SmallestN;
		estimOverSampX = double(nx)/double(Nx_nom);

		double HalfLambRz = 0.5*WavelengthIn_m*RobsZ;
		double zStartRel = zStart - zc;
		double zEndAbs = zStart + (nz - 1)*zStep;
		double zEndRel = zEndAbs - zc;
		double dzStart = ::fabs(HalfLambRz/zStartRel);
		double dzEnd = ::fabs(HalfLambRz/zEndRel);
		double dz = ((dzStart < dzEnd)? dzStart : dzEnd); // /NxNzOversamplingFactor;
		dz /= 1.2;
		long Nz_nom = long(::fabs(zEndRel - zStartRel)/dz) + 1;
		if(((Nz_nom>>1)<<1) != Nz_nom) Nz_nom++;
		FFT.NextCorrectNumberForFFT(Nz_nom);
		//if(Nz_nom < SmallestN) Nz_nom = SmallestN;
		estimOverSampZ = double(nz)/double(Nz_nom);
	}
// Continue here if(DistrInfoDat.CoordOrAngPresentation == AngPres)
}

//*************************************************************************

void srTSRWRadStructAccessData::CopyBaseRadData(float* pInBaseRadX, float* pInBaseRadZ)
{
	//long LenRadData = (ne << 1)*nx*nz;
	long long LenRadData = (ne << 1)*((long long)nx)*((long long)nz);
	bool NeedRadX = (LenRadData > 0) && (pInBaseRadX != 0) && (pBaseRadX != 0);
	bool NeedRadZ = (LenRadData > 0) && (pInBaseRadZ != 0) && (pBaseRadZ != 0);

	if(NeedRadX)
	{
		float *tBaseRadX = pBaseRadX;
		float *tInBaseRadX = pInBaseRadX;
		//for(long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
		for(long long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
		BaseRadWasEmulated = true;
	}
	if(NeedRadZ)
	{
		float *tBaseRadZ = pBaseRadZ;
		float *tInBaseRadZ = pInBaseRadZ;
		//for(long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
		for(long long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
		BaseRadWasEmulated = true;
	}
}

//*************************************************************************

//void srTSRWRadStructAccessData::CopyStatMomData(float* pInMomX, float* pInMomZ)
void srTSRWRadStructAccessData::CopyStatMomData(double* pInMomX, double* pInMomZ) //OC130311
{
	int LenMom = 11*ne; // to steer
	if((pInMomX != 0) && (pMomX != 0))
	{
		//float *tMomX = pMomX, *tInMomX = pInMomX;
		double *tMomX = pMomX, *tInMomX = pInMomX;
		for(int i=0; i<LenMom; i++) *(tMomX++) = *(tInMomX++);
		MomWereEmulated = true;
	}
	if((pInMomZ != 0) && (pMomZ != 0))
	{
		//float *tMomZ = pMomZ, *tInMomZ = pInMomZ;
		double *tMomZ = pMomZ, *tInMomZ = pInMomZ;
		for(int i=0; i<LenMom; i++) *(tMomZ++) = *(tInMomZ++);
		MomWereEmulated = true;
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::CopyElectronBeamData(DOUBLE* pInElecBeam)
{
	if((pInElecBeam != 0) && (pElecBeam != 0)) 
	{
		const int LenElecData = 30; // to steer
		DOUBLE *tElecBeam = pElecBeam, *tInElecBeam = pInElecBeam;
		for(int i=0; i<LenElecData; i++) *(tElecBeam++) = *(tInElecBeam++);
		ElectronBeamEmulated = true;
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::Copy4x4PropMatrData(DOUBLE* pIn4x4PropMatr)
{
	if((pIn4x4PropMatr != 0) && (p4x4PropMatr != 0))
	{
		const int Len4x4PropMatr = 20; //OC fix 16082004 //16; // to steer 
		DOUBLE *t4x4PropMatr = p4x4PropMatr, *tIn4x4PropMatr = pIn4x4PropMatr;
		for(int i=0; i<Len4x4PropMatr; i++) *(t4x4PropMatr++) = *(tIn4x4PropMatr++);
		PropMatrWasEmulated = true;
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::CopyWfrAuxData(DOUBLE* pInWfrAuxData)
{
	if((pInWfrAuxData != 0) && (pWfrAuxData != 0))
	{
		const int LenWfrAuxData = 12; // to steer
		DOUBLE *tWfrAuxData = pWfrAuxData, *tInWfrAuxData = pInWfrAuxData;
		for(int i=0; i<LenWfrAuxData; i++) *(tWfrAuxData++) = *(tInWfrAuxData++);
		WfrAuxDataWasEmulated = true;
	}
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(srTSRWRadStructAccessData* pInRadStruct)
{// copies all data/substructures for internal use only
	Initialize();
	
	if(pInRadStruct == 0) return;
	srTSRWRadStructAccessData& InRadStruct = *pInRadStruct;

	//long LenRadData = (InRadStruct.ne << 1)*(InRadStruct.nx)*(InRadStruct.nz);
	long long LenRadData = (InRadStruct.ne << 1)*((long long)InRadStruct.nx)*((long long)InRadStruct.nz);
	bool NeedRadX = (LenRadData > 0) && (InRadStruct.pBaseRadX != 0);
	bool NeedRadZ = (LenRadData > 0) && (InRadStruct.pBaseRadZ != 0);

	if(NeedRadX)
	{
		pBaseRadX = new float[LenRadData];
		float *tBaseRadX = pBaseRadX;
		float *tInBaseRadX = InRadStruct.pBaseRadX;
		//for(long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
		for(long long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
		BaseRadWasEmulated = true;
	}
	if(NeedRadZ)
	{
		pBaseRadZ = new float[LenRadData];
		float *tBaseRadZ = pBaseRadZ;
		float *tInBaseRadZ = InRadStruct.pBaseRadZ;
		//for(long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
		for(long long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
		BaseRadWasEmulated = true;
	}

	eStep = InRadStruct.eStep;
	eStart = InRadStruct.eStart;
	xStep = InRadStruct.xStep;
	xStart = InRadStruct.xStart;
	zStep = InRadStruct.zStep;
	zStart = InRadStruct.zStart;
	ne = InRadStruct.ne;
	nx = InRadStruct.nx;
	nz = InRadStruct.nz;
	xStartTr = InRadStruct.xStartTr;
	zStartTr = InRadStruct.zStartTr;

	RobsX = InRadStruct.RobsX;
	RobsZ = InRadStruct.RobsZ;
	RobsXAbsErr = InRadStruct.RobsXAbsErr;
	RobsZAbsErr = InRadStruct.RobsZAbsErr;
	xc = InRadStruct.xc;
	zc = InRadStruct.zc;
	xWfrMin = InRadStruct.xWfrMin;
	xWfrMax = InRadStruct.xWfrMax;
	zWfrMin = InRadStruct.zWfrMin;
	zWfrMax = InRadStruct.zWfrMax;

	WfrEdgeCorrShouldBeDone = InRadStruct.WfrEdgeCorrShouldBeDone;
	UnderSamplingX = InRadStruct.UnderSamplingX;
	UnderSamplingZ = InRadStruct.UnderSamplingZ;
	AllowAutoSwitchToPropInUnderSamplingMode = InRadStruct.AllowAutoSwitchToPropInUnderSamplingMode;
	InvUnderSamplingThreshold = InRadStruct.InvUnderSamplingThreshold;

	if(InRadStruct.pResAfter != 0) 
	{
		pResAfter = new srTRadResize(*InRadStruct.pResAfter);
		ResAfterWasEmulated = true;
	}

	Pres = InRadStruct.Pres;
	PresT = InRadStruct.PresT;
	LengthUnit = InRadStruct.LengthUnit;
	PhotEnergyUnit = InRadStruct.PhotEnergyUnit;

	if(InRadStruct.pElecBeam != 0) 
	{
		const int LenElecData = 30; // to steer
		pElecBeam = new DOUBLE[LenElecData << 1];
		DOUBLE *tElecBeam = pElecBeam, *tInElecBeam = InRadStruct.pElecBeam;
		for(int i=0; i<LenElecData; i++) *(tElecBeam++) = *(tInElecBeam++);
		ElectronBeamEmulated = true;
	}
	if(InRadStruct.p4x4PropMatr != 0)
	{
		const int Len4x4PropMatr = 16; // to steer
		p4x4PropMatr = new DOUBLE[Len4x4PropMatr << 1];
		DOUBLE *t4x4PropMatr = p4x4PropMatr, *tIn4x4PropMatr = InRadStruct.p4x4PropMatr;
		for(int i=0; i<Len4x4PropMatr; i++) *(t4x4PropMatr++) = *(tIn4x4PropMatr++);
		PropMatrWasEmulated = true;
	}
	if(InRadStruct.pMomX != 0)
	{
		//const int LenMomX = 11; // to steer
		//pMomX = new float[LenMomX << 1];
		int LenMomX = ne*11; // to steer
		//pMomX = new float[LenMomX];
		//float *tMomX = pMomX, *tInMomX = InRadStruct.pMomX;
		pMomX = new double[LenMomX]; //OC130311
		double *tMomX = pMomX, *tInMomX = InRadStruct.pMomX;

		for(int i=0; i<LenMomX; i++) *(tMomX++) = *(tInMomX++);
		MomWereEmulated = true;
	}
	if(InRadStruct.pMomZ != 0)
	{
		const int LenMomZ = 11; // to steer
		//pMomZ = new float[LenMomZ << 1];
		//float *tMomZ = pMomZ, *tInMomZ = InRadStruct.pMomZ;
		pMomZ = new double[LenMomZ << 1]; //OC130311
		double *tMomZ = pMomZ, *tInMomZ = InRadStruct.pMomZ;

		for(int i=0; i<LenMomZ; i++) *(tMomZ++) = *(tInMomZ++);
		MomWereEmulated = true;
	}
	if(InRadStruct.pWfrAuxData != 0)
	{
		const int LenWfrAuxData = 12; // to steer
		pWfrAuxData = new DOUBLE[LenWfrAuxData];
		DOUBLE *tWfrAuxData = pWfrAuxData, *tInWfrAuxData = InRadStruct.pWfrAuxData;
		for(int i=0; i<LenWfrAuxData; i++) *(tWfrAuxData++) = *(tInWfrAuxData++);
		WfrAuxDataWasEmulated = true;
	}

	AuxLong1 = InRadStruct.AuxLong1;
	AuxLong2 = InRadStruct.AuxLong2;
	AuxLong3 = InRadStruct.AuxLong3;
	AuxLong4 = InRadStruct.AuxLong4;
	
	WfrQuadTermCanBeTreatedAtResizeX = InRadStruct.WfrQuadTermCanBeTreatedAtResizeX; // is used at the time of one resize only
	WfrQuadTermCanBeTreatedAtResizeZ = InRadStruct.WfrQuadTermCanBeTreatedAtResizeZ;

	//add/copy more members here, if they appear
}

//*************************************************************************

srTSRWRadStructAccessData::srTSRWRadStructAccessData(const srTSRWRadStructAccessData& inRad, bool createNewEmulStruct)
{
	*this = inRad;

	if(createNewEmulStruct)
	{
		if(BaseRadWasEmulated)
		{
			//long LenRadData = (inRad.ne << 1)*(inRad.nx)*(inRad.nz);
			long long LenRadData = (inRad.ne << 1)*((long long)inRad.nx)*((long long)inRad.nz);
			bool NeedRadX = (LenRadData > 0) && (inRad.pBaseRadX != 0);
			bool NeedRadZ = (LenRadData > 0) && (inRad.pBaseRadZ != 0);
			if(NeedRadX)
			{
				pBaseRadX = new float[LenRadData];
				float *tBaseRadX = pBaseRadX;
				float *tInBaseRadX = inRad.pBaseRadX;
				//for(long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
				for(long long i=0; i<LenRadData; i++) *(tBaseRadX++) = *(tInBaseRadX++);
			}
			if(NeedRadZ)
			{
				pBaseRadZ = new float[LenRadData];
				float *tBaseRadZ = pBaseRadZ;
				float *tInBaseRadZ = inRad.pBaseRadZ;
				//for(long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
				for(long long i=0; i<LenRadData; i++) *(tBaseRadZ++) = *(tInBaseRadZ++);
			}
		}
		if(ResAfterWasEmulated && (inRad.pResAfter != 0))
		{
			pResAfter = new srTRadResize(*inRad.pResAfter);
		}
		if(ElectronBeamEmulated && (inRad.pElecBeam != 0))
		{
			const int LenElecData = 30; // to steer
			pElecBeam = new DOUBLE[LenElecData << 1];
			DOUBLE *tElecBeam = pElecBeam, *tInElecBeam = inRad.pElecBeam;
			for(int i=0; i<LenElecData; i++) *(tElecBeam++) = *(tInElecBeam++);
		}
		if(PropMatrWasEmulated && (inRad.p4x4PropMatr != 0))
		{
			const int Len4x4PropMatr = 16; // to steer
			p4x4PropMatr = new DOUBLE[Len4x4PropMatr << 1];
			DOUBLE *t4x4PropMatr = p4x4PropMatr, *tIn4x4PropMatr = inRad.p4x4PropMatr;
			for(int i=0; i<Len4x4PropMatr; i++) *(t4x4PropMatr++) = *(tIn4x4PropMatr++);
		}
		if(MomWereEmulated)
		{
			if(inRad.pMomX != 0)
			{
				int LenMomX = ne*11; // to steer
				pMomX = new double[LenMomX]; //OC130311
				double *tMomX = pMomX, *tInMomX = inRad.pMomX;

				for(int i=0; i<LenMomX; i++) *(tMomX++) = *(tInMomX++);
			}
			if(inRad.pMomZ != 0)
			{
				const int LenMomZ = 11; // to steer
				pMomZ = new double[LenMomZ << 1]; //OC130311
				double *tMomZ = pMomZ, *tInMomZ = inRad.pMomZ;

				for(int i=0; i<LenMomZ; i++) *(tMomZ++) = *(tInMomZ++);
			}
		}
		if(WfrAuxDataWasEmulated && (inRad.pWfrAuxData != 0))
		{
			const int LenWfrAuxData = 12; // to steer
			pWfrAuxData = new DOUBLE[LenWfrAuxData];
			DOUBLE *tWfrAuxData = pWfrAuxData, *tInWfrAuxData = inRad.pWfrAuxData;
			for(int i=0; i<LenWfrAuxData; i++) *(tWfrAuxData++) = *(tInWfrAuxData++);
		}
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::Initialize()
{
	wRad = NIL;
	pBaseRadX = 0; pBaseRadZ = 0; // This is checked in FinishWorkingWithSRWRadStruct !!!
	pBaseRadXaux = 0; pBaseRadZaux = 0; //OC151115

	wRadX = wRadZ = NIL;
	BaseRadWasEmulated = false;
	
	UseStartTrToShiftAtChangingRepresToCoord = false;
	//UseStartTrToShiftAtChangingRepresToTime = false; //OC091115
	avgT = 0; //OC101115
	
	DoNotResizeAfter = false;
	ResAfterWasEmulated = false;
	
	pElecBeam = 0; wElecBeam = NIL;
	wTrj = NIL;
	ElectronBeamEmulated = 0;
	
	p4x4PropMatr = 0; w4x4PropMatr = NIL;
	PropMatrWasEmulated = false;
	
	pMomX = pMomZ = 0;
	wMomX = wMomZ = NIL;
	MomWereEmulated = false;
	MomWereCalcNum = false;
	
	pWfrAuxData = 0; wWfrAuxData = NIL;
	WfrAuxDataWasEmulated = false;
	
	WfrEdgeCorrShouldBeDone = 0; //1; // Turn it on/off manually
	pResAfter = 0;
	
	UnderSamplingX = UnderSamplingZ = 1.;
	AllowAutoSwitchToPropInUnderSamplingMode = 0;
	
	ElecFldUnit = 1;
	
	WfrQuadTermCanBeTreatedAtResizeX = false; // is used at the time of one resize only
	WfrQuadTermCanBeTreatedAtResizeZ = false;

	m_xQuadPhaseTermWasSubtracted = false;
	m_zQuadPhaseTermWasSubtracted = false; 
	m_xLinOnlyPhaseTermWasSubtracted = false;
	m_zLinOnlyPhaseTermWasSubtracted = false;
	m_dxcSub = m_dzcSub = 0;

	m_pExtWfr = 0;
	m_newExtWfrCreateNotAllowed = false;
}

//*************************************************************************

void srTSRWRadStructAccessData::DisposeEmulatedStructs()
{
	if(BaseRadWasEmulated)
	{
		if(pBaseRadX != 0) delete[] pBaseRadX;
		if(pBaseRadZ != 0) delete[] pBaseRadZ;
		pBaseRadX = pBaseRadZ = 0;
		BaseRadWasEmulated = false;
	}
	if(ResAfterWasEmulated)
	{
		if(pResAfter != 0) delete pResAfter;
		pResAfter = 0;
		ResAfterWasEmulated = false;
	}
	
	if(ElectronBeamEmulated && (pElecBeam != 0)) delete[] pElecBeam;
	pElecBeam = 0;
	ElectronBeamEmulated = 0;
	
	if(PropMatrWasEmulated)
	{
		if(p4x4PropMatr != 0) delete[] p4x4PropMatr;
		p4x4PropMatr = 0;
		PropMatrWasEmulated = false;
	}
	if(MomWereEmulated)
	{
		if(pMomX != 0) delete[] pMomX;
		if(pMomZ != 0) delete[] pMomZ;
		pMomX = pMomZ = 0;
		MomWereEmulated = false;
	}
	if(WfrAuxDataWasEmulated)
	{
		if(pWfrAuxData != 0) delete[] pWfrAuxData;
		pWfrAuxData = 0;
		WfrAuxDataWasEmulated = false;
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::PreserveLogicsOfWfrLimitsAtRangeResizing(srTSRWRadStructAccessData* pOldRadData, char x_or_z)
{
	if(x_or_z == 'x')
	{
		double xDistTol = 0.01*xStep;
		if((::fabs(pOldRadData->xWfrMin - pOldRadData->xStart) < xDistTol) && (::fabs(pOldRadData->xStart + pOldRadData->nx*pOldRadData->xStep - pOldRadData->xWfrMax) < xDistTol))
		{
			xWfrMin = xStart; xWfrMax = xStart + nx*xStep;
		}
		else
		{
			xWfrMin = pOldRadData->xWfrMin; xWfrMax = pOldRadData->xWfrMax;
		}
	}
	else //if((x_or_z == 'z') || (x_or_z == 'y'))
	{
		double zDistTol = 0.01*zStep;
		if((::fabs(pOldRadData->zWfrMin - pOldRadData->zStart) < zDistTol) && (::fabs(pOldRadData->zStart + pOldRadData->nz*pOldRadData->zStep - pOldRadData->zWfrMax) < zDistTol))
		{
			zWfrMin = zStart; zWfrMax = zStart + nz*zStep;
		}
		else
		{
			zWfrMin = pOldRadData->zWfrMin; zWfrMax = pOldRadData->zWfrMax;
		}
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::FindMinMaxReE(srTMinMaxEParam& a)
{
	float &MaxReEx = a.MaxReEx, &MaxImEx = a.MaxImEx, &MaxReEz = a.MaxReEz, &MaxImEz = a.MaxImEz, &MinReEx = a.MinReEx, &MinImEx = a.MinImEx, &MinReEz = a.MinReEz, &MinImEz = a.MinImEz;
	long &xIndMaxReEx = a.xIndMaxReEx, &xIndMaxImEx = a.xIndMaxImEx, &xIndMaxReEz = a.xIndMaxReEz, &xIndMaxImEz = a.xIndMaxImEz, &xIndMinReEx = a.xIndMinReEx, &xIndMinImEx = a.xIndMinImEx, &xIndMinReEz = a.xIndMinReEz, &xIndMinImEz = a.xIndMinImEz;
	long &zIndMaxReEx = a.zIndMaxReEx, &zIndMaxImEx = a.zIndMaxImEx, &zIndMaxReEz = a.zIndMaxReEz, &zIndMaxImEz = a.zIndMaxImEz, &zIndMinReEx = a.zIndMinReEx, &zIndMinImEx = a.zIndMinImEx, &zIndMinReEz = a.zIndMinReEz, &zIndMinImEz = a.zIndMinImEz;
	MaxReEx = MaxImEx = MaxReEz = MaxImEz = (float)(-1.E+23);
	MinReEx = MinImEx = MinReEz = MinImEz = (float)(1.E+23);
	float *tReEx = pBaseRadX, *tImEx = pBaseRadX + 1, *tReEz = pBaseRadZ, *tImEz = pBaseRadZ + 1;
	for(long iz=0; iz<nz; iz++)
	{
		for(long ix=0; ix<nx; ix++)
		{
			if(*tReEx > MaxReEx)
			{
				MaxReEx = *tReEx; xIndMaxReEx = ix; zIndMaxReEx = iz;
			}
			if(*tImEx > MaxImEx)
			{
				MaxImEx = *tImEx; xIndMaxImEx = ix; zIndMaxImEx = iz;
			}
			if(*tReEz > MaxReEz)
			{
				MaxReEz = *tReEz; xIndMaxReEz = ix; zIndMaxReEz = iz;
			}
			if(*tImEz > MaxImEz)
			{
				MaxImEz = *tImEz; xIndMaxImEz = ix; zIndMaxImEz = iz;
			}
			if(*tReEx < MinReEx)
			{
				MinReEx = *tReEx; xIndMinReEx = ix; zIndMinReEx = iz;
			}
			if(*tImEx < MinImEx)
			{
				MinImEx = *tImEx; xIndMinImEx = ix; zIndMinImEx = iz;
			}
			if(*tReEz < MinReEz)
			{
				MinReEz = *tReEz; xIndMinReEz = ix; zIndMinReEz = iz;
			}
			if(*tImEz < MinImEz)
			{
				MinImEz = *tImEz; xIndMinImEz = ix; zIndMinImEz = iz;
			}
			tReEx += 2; tImEx += 2; tReEz += 2; tImEz += 2;
		}
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::AllocElectronBeam()
{
	int MaxLenElecBeam = 50;
	pElecBeam = new DOUBLE[MaxLenElecBeam];
	if(pElecBeam == 0) throw MEMORY_ALLOCATION_FAILURE;
	ElectronBeamEmulated = 1;
	
	DOUBLE *tElecBeam = pElecBeam;
	for(int i=0; i<MaxLenElecBeam; i++) *(tElecBeam++) = 0.;
}

//*************************************************************************

int srTSRWRadStructAccessData::EmulateElectronBeamStruct(srTEbmDat& EbmDat)
{
	if(pElecBeam == 0)
	{
		pElecBeam = new DOUBLE[50];
		if(pElecBeam == 0) return MEMORY_ALLOCATION_FAILURE;
		ElectronBeamEmulated = 1;
	}
	
	DOUBLE *tElecBeam = pElecBeam;
	for(int i=0; i<50; i++) *(tElecBeam++) = 0.;
	
	*pElecBeam = EbmDat.Energy;
	*(pElecBeam + 1) = EbmDat.Current;
	*(pElecBeam + 2) = EbmDat.x0;
	*(pElecBeam + 3) = EbmDat.dxds0;
	*(pElecBeam + 4) = EbmDat.z0;
	*(pElecBeam + 5) = EbmDat.dzds0;
	*(pElecBeam + 6) = EbmDat.s0;
	
	*(pElecBeam + 13) = EbmDat.SigmaRelE;

	*(pElecBeam + 20) = EbmDat.Mxx;
	*(pElecBeam + 21) = EbmDat.Mxxp;
	*(pElecBeam + 22) = EbmDat.Mxpxp;
	*(pElecBeam + 23) = EbmDat.Mzz;
	*(pElecBeam + 24) = EbmDat.Mzzp;
	*(pElecBeam + 25) = EbmDat.Mzpzp;
	*(pElecBeam + 26) = EbmDat.Mxz;
	*(pElecBeam + 27) = EbmDat.Mxpz;
	*(pElecBeam + 28) = EbmDat.Mxzp;
	*(pElecBeam + 29) = EbmDat.Mxpzp;
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::EmulateElectronBeamStruct(const SRWLPartBeam& srwlPartBeam)
{
	const int nValElecBeam = 60;
	if(pElecBeam == 0)
	{
		pElecBeam = new DOUBLE[nValElecBeam];
		if(pElecBeam == 0) return MEMORY_ALLOCATION_FAILURE;
		ElectronBeamEmulated = 1;
	}
	DOUBLE *tElecBeam = pElecBeam;
	for(int i=0; i<nValElecBeam; i++) *(tElecBeam++) = 0.;

	const double elecEn0 = 0.51099890221e-03; //[GeV]
	*pElecBeam = srwlPartBeam.partStatMom1.gamma*srwlPartBeam.partStatMom1.relE0*elecEn0;
	*(pElecBeam + 1) = srwlPartBeam.Iavg;
	*(pElecBeam + 2) = srwlPartBeam.partStatMom1.x;
	*(pElecBeam + 3) = srwlPartBeam.partStatMom1.xp;
	*(pElecBeam + 4) = srwlPartBeam.partStatMom1.y;
	*(pElecBeam + 5) = srwlPartBeam.partStatMom1.yp;
	*(pElecBeam + 6) = srwlPartBeam.partStatMom1.z;

	*(pElecBeam + 13) = ::sqrt(srwlPartBeam.arStatMom2[10]); //EbmDat.SigmaRelE;

	*(pElecBeam + 20) = srwlPartBeam.arStatMom2[0]; //EbmDat.Mxx;
	*(pElecBeam + 21) = srwlPartBeam.arStatMom2[1]; //EbmDat.Mxxp;
	*(pElecBeam + 22) = srwlPartBeam.arStatMom2[2]; //EbmDat.Mxpxp;
	*(pElecBeam + 23) = srwlPartBeam.arStatMom2[3]; //EbmDat.Mzz;
	*(pElecBeam + 24) = srwlPartBeam.arStatMom2[4]; //EbmDat.Mzzp;
	*(pElecBeam + 25) = srwlPartBeam.arStatMom2[5]; //EbmDat.Mzpzp;
	*(pElecBeam + 26) = srwlPartBeam.arStatMom2[6]; //EbmDat.Mxz;
	*(pElecBeam + 27) = srwlPartBeam.arStatMom2[7]; //EbmDat.Mxpz;
	*(pElecBeam + 28) = srwlPartBeam.arStatMom2[8]; //EbmDat.Mxzp;
	*(pElecBeam + 29) = srwlPartBeam.arStatMom2[9]; //EbmDat.Mxpzp;
	*(pElecBeam + 33) = srwlPartBeam.arStatMom2[11]; //EbmDat.Mss;
	*(pElecBeam + 34) = srwlPartBeam.arStatMom2[12]; //EbmDat.Mse;
	*(pElecBeam + 35) = srwlPartBeam.arStatMom2[13]; //EbmDat.Mxe;
	*(pElecBeam + 36) = srwlPartBeam.arStatMom2[14]; //EbmDat.Mxpe;
	*(pElecBeam + 37) = srwlPartBeam.arStatMom2[15]; //EbmDat.Mze;
	*(pElecBeam + 38) = srwlPartBeam.arStatMom2[16]; //EbmDat.Mzpe;
	*(pElecBeam + 39) = srwlPartBeam.arStatMom2[17]; //EbmDat.Mxs;
	*(pElecBeam + 40) = srwlPartBeam.arStatMom2[18]; //EbmDat.Mxps;
	*(pElecBeam + 41) = srwlPartBeam.arStatMom2[19]; //EbmDat.Mzs;
	*(pElecBeam + 42) = srwlPartBeam.arStatMom2[20]; //EbmDat.Mzps;
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::EmulateElectronBeamStruct(srTGsnBeam& GsnBeam)
{
	if(pElecBeam == 0)
	{
		pElecBeam = new DOUBLE[50];
		if(pElecBeam == 0) return MEMORY_ALLOCATION_FAILURE;
		ElectronBeamEmulated = 1;
	}
	DOUBLE *tElecBeam = pElecBeam;
	for(int i=0; i<50; i++) *(tElecBeam++) = 0.;
	
    srTEbmDat& EbmDat = GsnBeam.EbmDat;

	*pElecBeam = 1.;
	*(pElecBeam + 1) = 1.;
	*(pElecBeam + 2) = EbmDat.x0;
	*(pElecBeam + 3) = EbmDat.dxds0;
	*(pElecBeam + 4) = EbmDat.z0;
	*(pElecBeam + 5) = EbmDat.dzds0;
	*(pElecBeam + 6) = EbmDat.s0;

	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::OutElectronBeamStruct(srTEbmDat& EbmDat)
{
	if(pElecBeam == 0) return ELECTRON_BEAM_WAS_NOT_SET_UP;

	EbmDat.Energy = *pElecBeam;
	EbmDat.Current = *(pElecBeam + 1);
	EbmDat.x0 = *(pElecBeam + 2);
	EbmDat.dxds0 = *(pElecBeam + 3);
	EbmDat.z0 = *(pElecBeam + 4);
	EbmDat.dzds0 = *(pElecBeam + 5);
	EbmDat.s0 = *(pElecBeam + 6);
	
	EbmDat.SigmaRelE = *(pElecBeam + 13);
	EbmDat.Mee = (EbmDat.SigmaRelE)*(EbmDat.SigmaRelE);

	EbmDat.Mxx = *(pElecBeam + 20);
	EbmDat.Mxxp = *(pElecBeam + 21);
	EbmDat.Mxpxp = *(pElecBeam + 22);
	EbmDat.Mzz = *(pElecBeam + 23);
	EbmDat.Mzzp = *(pElecBeam + 24);
	EbmDat.Mzpzp = *(pElecBeam + 25);
	EbmDat.Mxz = *(pElecBeam + 26);
	EbmDat.Mxpz = *(pElecBeam + 27);
	EbmDat.Mxzp = *(pElecBeam + 28);
	EbmDat.Mxpzp = *(pElecBeam + 29);

	EbmDat.TypeDistrTransverse = 2; // 1- uniform; 2- Gaussian; 3- parabolic
	EbmDat.TypeDistrLongitudinal = 2; // 1- infinite uniform; 2- gaussian
	//consider setting this from data

	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::OutElectronBeamStruct(SRWLPartBeam& srwlPartBeam)
{
	if(pElecBeam == 0) return ELECTRON_BEAM_WAS_NOT_SET_UP;

	const double elecEn0 = 0.51099890221e-03; //[GeV]
	//*pElecBeam = srwlPartBeam.partStatMom1.gamma*srwlPartBeam.partStatMom1.relE0*elecEn0;
	srwlPartBeam.partStatMom1.gamma = pElecBeam[0]/(srwlPartBeam.partStatMom1.relE0*elecEn0);

	srwlPartBeam.Iavg = pElecBeam[1];
	srwlPartBeam.partStatMom1.x = pElecBeam[2];
	srwlPartBeam.partStatMom1.xp = pElecBeam[3];
	srwlPartBeam.partStatMom1.y = pElecBeam[4];
	srwlPartBeam.partStatMom1.yp = pElecBeam[5];
	srwlPartBeam.partStatMom1.z = pElecBeam[6];

	srwlPartBeam.arStatMom2[0] = pElecBeam[20]; //EbmDat.Mxx;
	srwlPartBeam.arStatMom2[1] = pElecBeam[21]; //EbmDat.Mxxp;
	srwlPartBeam.arStatMom2[2] = pElecBeam[22]; //EbmDat.Mxpxp;
	srwlPartBeam.arStatMom2[3] = pElecBeam[23]; //EbmDat.Mzz;
	srwlPartBeam.arStatMom2[4] = pElecBeam[24]; //EbmDat.Mzzp;
	srwlPartBeam.arStatMom2[5] = pElecBeam[25]; //EbmDat.Mzpzp;
	srwlPartBeam.arStatMom2[6] = pElecBeam[26]; //EbmDat.Mxz;
	srwlPartBeam.arStatMom2[7] = pElecBeam[27]; //EbmDat.Mxpz;
	srwlPartBeam.arStatMom2[8] = pElecBeam[28]; //EbmDat.Mxzp;
	srwlPartBeam.arStatMom2[9] = pElecBeam[29]; //EbmDat.Mxpzp;
	srwlPartBeam.arStatMom2[10] = pElecBeam[13]*pElecBeam[13];
	//*(pElecBeam + 13) = ::sqrt(srwlPartBeam.arStatMom2[10]); //EbmDat.SigmaRelE;

	srwlPartBeam.arStatMom2[11] = pElecBeam[33]; //EbmDat.Mss;
	srwlPartBeam.arStatMom2[12] = pElecBeam[34]; //EbmDat.Mse;
	srwlPartBeam.arStatMom2[13] = pElecBeam[35]; //EbmDat.Mxe;
	srwlPartBeam.arStatMom2[14] = pElecBeam[36]; //EbmDat.Mxpe;
	srwlPartBeam.arStatMom2[15] = pElecBeam[37]; //EbmDat.Mze;
	srwlPartBeam.arStatMom2[16] = pElecBeam[38]; //EbmDat.Mzpe;
	srwlPartBeam.arStatMom2[17] = pElecBeam[39]; //EbmDat.Mxs;
	srwlPartBeam.arStatMom2[18] = pElecBeam[40]; //EbmDat.Mxps;
	srwlPartBeam.arStatMom2[19] = pElecBeam[41]; //EbmDat.Mzs;
	srwlPartBeam.arStatMom2[20] = pElecBeam[42]; //EbmDat.Mzps;
	return 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::EstimateAndSetUnderSampling()
{
	//double HalfWavelength_m = 0.5*1.239854E-06/eStart; // Assumes eStart in eV
	double HalfWavelength_m = 0.5*1.239842E-06/eStart; // Assumes eStart in eV
	double HalfLambRx = HalfWavelength_m*RobsX;
	double HalfLambRz = HalfWavelength_m*RobsZ;
	
	double xStartRel = xStart - xc;
	double xEndRel = xStartRel + xStep*(nx - 1);
	double dxStart = ::fabs(HalfLambRx/xStartRel);
	double dxEnd = ::fabs(HalfLambRx/xEndRel);
	double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/1.2;
	double Nx = ::fabs(xEndRel - xStartRel)/dx + 1.;
	double TestInvUnderSampX = double(nx)/Nx;
	if(TestInvUnderSampX <= InvUnderSamplingThreshold) UnderSamplingX = 1./TestInvUnderSampX;
	
	double zStartRel = zStart - zc;
	double zEndRel = zStartRel + zStep*(nz - 1);
	double dzStart = ::fabs(HalfLambRz/zStartRel);
	double dzEnd = ::fabs(HalfLambRz/zEndRel);
	double dz = ((dzStart < dzEnd)? dzStart : dzEnd)/1.2;
	double Nz = ::fabs(zEndRel - zStartRel)/dz + 1.;
	double TestInvUnderSampZ = double(nz)/Nz;
	if(TestInvUnderSampZ <= InvUnderSamplingThreshold) UnderSamplingZ = 1./TestInvUnderSampZ;
}

//*************************************************************************

void srTSRWRadStructAccessData::ZeroPtrs()
{
	pBaseRadX = pBaseRadZ = 0;
	pResAfter = 0;
	pElecBeam = 0;
	p4x4PropMatr = 0;
	pMomX = pMomZ = 0;
	pWfrAuxData = 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::ReAllocBaseRadAccordingToNeNxNz(char PolarizComp)
{
	//long LenRadData = (ne << 1)*nx*nz;
	long long LenRadData = (ne << 1)*((long long)nx)*((long long)nz);
	bool TreatPolCompX = ((PolarizComp == 0) || (PolarizComp == 'x')) && (LenRadData > 0);
	bool TreatPolCompZ = ((PolarizComp == 0) || (PolarizComp == 'z')) && (LenRadData > 0);

	if(TreatPolCompX)
	{
		if(pBaseRadX != 0) { delete[] pBaseRadX; pBaseRadX = 0;}
		pBaseRadX = new float[LenRadData];
		if(pBaseRadX == 0) return MEMORY_ALLOCATION_FAILURE;
		BaseRadWasEmulated = true;
	}
	if(TreatPolCompZ)
	{
		if(pBaseRadZ != 0) { delete[] pBaseRadZ; pBaseRadZ = 0;}
		pBaseRadZ = new float[LenRadData];
		if(pBaseRadZ == 0) return MEMORY_ALLOCATION_FAILURE;
		BaseRadWasEmulated = true;
	}
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::AllocBaseRadAccordingToNeNxNz(char PolarizComp)
{
	//long LenRadData = (ne << 1)*nx*nz;
	long long LenRadData = (ne << 1)*((long long)nx)*((long long)nz);
	bool TreatPolCompX = ((PolarizComp == 0) || (PolarizComp == 'x')) && (LenRadData > 0);
	bool TreatPolCompZ = ((PolarizComp == 0) || (PolarizComp == 'z')) && (LenRadData > 0);

	if(TreatPolCompX)
	{
		pBaseRadX = 0;
		pBaseRadX = new float[LenRadData];
		if(pBaseRadX == 0) return MEMORY_ALLOCATION_FAILURE;
		BaseRadWasEmulated = true;
	}
	if(TreatPolCompZ)
	{
		pBaseRadZ = 0;
		pBaseRadZ = new float[LenRadData];
		if(pBaseRadZ == 0) return MEMORY_ALLOCATION_FAILURE;
		BaseRadWasEmulated = true;
	}
	return 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::DeAllocBaseRadAccordingToNeNxNz(char PolarizComp)
{
	//long LenRadData = (ne << 1)*nx*nz;
	long long LenRadData = (ne << 1)*((long long)nx)*((long long)nz);
	bool TreatPolCompX = ((PolarizComp == 0) || (PolarizComp == 'x')) && (LenRadData > 0);
	bool TreatPolCompZ = ((PolarizComp == 0) || (PolarizComp == 'z')) && (LenRadData > 0);

	if(TreatPolCompX)
	{
		if(pBaseRadX != 0) { delete[] pBaseRadX; pBaseRadX = 0;}
	}
	if(TreatPolCompZ)
	{
		if(pBaseRadZ != 0) { delete[] pBaseRadZ; pBaseRadZ = 0;}
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::UpdateObsParam(srTWfrSmp& DistrInfoDat)
{
	xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
	xStart = DistrInfoDat.xStart;
	nx = DistrInfoDat.nx;
	zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
	zStart = DistrInfoDat.zStart;
	nz = DistrInfoDat.nz;
}

//*************************************************************************

void srTSRWRadStructAccessData::SetObsParamFromWfr(srTWfrSmp& smp)
{
	smp.Initialize();

	//if(Pres == 0)
	if(PresT == 0) //OC071112
	{
		smp.LambStart = eStart;
		smp.LambEnd = eStart + (ne - 1)*eStep;
		smp.nLamb = ne;
		smp.tStart = smp.tEnd = 0;
		//smp.nt = 0;
		smp.nt = 1; //OC191215
	}
	else
	{
		smp.tStart = eStart;
		smp.tEnd = eStart + (ne - 1)*eStep;
		smp.nt = ne;
		smp.LambStart = smp.LambEnd = avgPhotEn;
		//smp.nLamb = 0;
		smp.nLamb = 1; //OC191215
	}

	smp.nx = nx;
	smp.xStart = xStart;
	smp.xEnd = xStart + (nx - 1)*xStep;

	smp.nz = nz;
	smp.zStart = zStart;
	smp.zEnd = zStart + (nz - 1)*zStep;

	smp.yStart = smp.yEnd = yStart; //longitudinal position of observation point

	smp.AllowAutoChoiceOfNxNzForPropagat = 0;
	smp.NxNzOversamplingParam = 0;
	smp.DimensionsWereSetAuto = 0;

	smp.obsPlaneIsTransv = true;

	smp.PresT = PresT;
	smp.PhotonEnergyWavelengthUnits = 1; //eV
	smp.BandwidthUnits = 0;
	smp.CoordUnits = 0; // 0- m, 1- mm
	smp.FluxComp = 0; // if !=0, Integrated flux is computed

	smp.TreatLambdaAsEnergyIn_eV = 1;
	smp.ShowPhaseOnly = 0;

	smp.InputWasModified = 1;
	smp.RadDistrDataContShouldBeRebuild = 0;
	smp.OnlyOnePoint = 0;
	smp.AssumeAllPhotonEnergies = 0;
	smp.AngPresToSpeedUpCoordPres = 0;

	strcpy(smp.LoopOrder, "yzxw");
	smp.DistrValType = StokesParam;
	smp.DistrPolariz = HorAndVer;
	smp.CoordOrAngPresentation = CoordPres;
}

//*************************************************************************

void srTSRWRadStructAccessData::SetupRadMomentsPtrs(srTMomentsPtrs& MomPtrsX, srTMomentsPtrs& MomPtrsZ)
{
	MomPtrsX = OneSetOfMomentsPtrs(pMomX); 
	MomPtrsZ = OneSetOfMomentsPtrs(pMomZ); 
}

//*************************************************************************

//srTMomentsPtrs srTSRWRadStructAccessData::OneSetOfMomentsPtrs(float* tMom)
srTMomentsPtrs srTSRWRadStructAccessData::OneSetOfMomentsPtrs(double* tMom) //OC130311
{
	return srTMomentsPtrs(tMom);
}

//*************************************************************************

void srTSRWRadStructAccessData::SetRadSamplingFromObs(srTWfrSmp& DistrInfoDat)
{
	eStart = DistrInfoDat.LambStart;
	eStep = (DistrInfoDat.nLamb > 1)? (DistrInfoDat.LambEnd - DistrInfoDat.LambStart)/(DistrInfoDat.nLamb - 1) : 0.;
	ne = DistrInfoDat.nLamb;
	
	xStart = DistrInfoDat.xStart;
	xStep = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
	nx = DistrInfoDat.nx;
	
	zStart = DistrInfoDat.zStart;
	zStep = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
	nz = DistrInfoDat.nz;

	PresT = DistrInfoDat.PresT; //OC191215

	//if(DistrInfoDat.nt > 1)
	if(DistrInfoDat.PresT == 1) //OC191215
	{
		eStart = DistrInfoDat.tStart;
		eStep = (DistrInfoDat.tEnd - DistrInfoDat.tStart)/(DistrInfoDat.nt - 1);
		ne = DistrInfoDat.nt;
		//PresT = 1; //OC191215
	}
	//else PresT = 0;
	
	// To walk around a bug in Igor
	if(eStep == 0.) { eStep = (eStart != 0.)? (1.e-08)*(::fabs(eStart)) : 1.e-10;}
	if(xStep == 0.) { xStep = (xStart != 0.)? (1.e-08)*(::fabs(xStart)) : 1.e-10;}
	if(zStep == 0.) { zStep = (zStart != 0.)? (1.e-08)*(::fabs(zStart)) : 1.e-10;}

	Pres = DistrInfoDat.CoordOrAngPresentation; //0- coord., 1- ang.
	LengthUnit = DistrInfoDat.CoordUnits; // 0- m; 1- mm; 
	PhotEnergyUnit = 0; // 0- eV; 1- keV; 

	//DistrInfoDat.PhotonEnergyWavelengthUnits; // 0- keV, 1- eV, 2- Ang, 3- nm, 4- micron
}

//*************************************************************************

void srTSRWRadStructAccessData::Alloc4x4PropMatr()
{
	const int Len4x4PropMatr = 20; //OC modified 16082004 //16; // to steer
	p4x4PropMatr = new DOUBLE[Len4x4PropMatr << 1];
	PropMatrWasEmulated = true;
}

//*************************************************************************

void srTSRWRadStructAccessData::AllocStatMom()
{
	const int LenMomX = 11; // to steer
	//pMomX = new float[LenMomX*ne];
	pMomX = new double[LenMomX*ne]; //OC130311
	const int LenMomZ = 11; // to steer
	//pMomZ = new float[LenMomZ*ne];
	pMomZ = new double[LenMomZ*ne];
	MomWereEmulated = true;
}

//*************************************************************************

void srTSRWRadStructAccessData::AllocWfrAux()
{
	const int LenWfrAuxData = 12; // to steer
	pWfrAuxData = new DOUBLE[LenWfrAuxData];
	WfrAuxDataWasEmulated = true;
}

//*************************************************************************

int srTSRWRadStructAccessData::FindAverageDistanceToSource(srTTrjDat& TrjDat, srTWfrSmp& DistrInfoDat, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc, srTParPrecElecFld* pPrecElecFld)
{// Should be called after the trajectory is already computed!!!
	
	double sStart = TrjDat.sStart;
	double sRange = (TrjDat.LenFieldData - 1)*TrjDat.sStep;
	double sEnd = sStart + sRange;
	//int NpVsS = TrjDat.LenFieldData;
	long long NpVsS = TrjDat.LenFieldData;

	if(pPrecElecFld != 0)
	{
		double sStartIntPrec = pPrecElecFld->sStartInt;
		double sEndIntPrec = pPrecElecFld->sEndInt;
		double sStart0 = sStart, sEnd0 = sEnd;
		bool SpecLimitsMayBeDefined = (sStartIntPrec < sEndIntPrec);
		if(SpecLimitsMayBeDefined && (sStartIntPrec > sStart) && (sStartIntPrec < sEnd)) 
		{
			sStart = sStartIntPrec;
		}
		if(SpecLimitsMayBeDefined && (sEndIntPrec > sStart) && (sEndIntPrec < sEnd)) 
		{
			sEnd = sEndIntPrec;
		}
		if((sStart != sStart0) || (sEnd != sEnd0))
		{
			NpVsS = (int)((sEnd - sStart)/TrjDat.sStep + 0.00001) + 1;
			sEnd = sStart + (NpVsS - 1)*TrjDat.sStep;
			sRange = sEnd - sStart;
		}
	}

	double Ysrc = 0.;

	char DistUnits = 0; // X, Z in mm
	double *TmpDataStorage = new double[TrjDat.LenFieldData << 2];
	if(TmpDataStorage == 0) return MEMORY_ALLOCATION_FAILURE;
	double *BtxArr = TmpDataStorage;
	double *BtzArr = TmpDataStorage + TrjDat.LenFieldData;
	double *xArr = TmpDataStorage + (TrjDat.LenFieldData << 1);
	double *zArr = TmpDataStorage + (TrjDat.LenFieldData*3);
	//TrjDat.CompTotalTrjDataTrjDisp(TrjDat.sStart, sEnd, TrjDat.LenFieldData, BtxArr, BtzArr, xArr, zArr, DistUnits);
	TrjDat.CompTotalTrjDataTrjDisp(sStart, sEnd, NpVsS, BtxArr, BtzArr, xArr, zArr, DistUnits);

	//int Len_mi_1 = TrjDat.LenFieldData - 1;
	// Len_mi_1 = NpVsS - 1;
	long long Len_mi_1 = NpVsS - 1;
	double *pBtx = BtxArr + Len_mi_1, *pBtz = BtzArr + Len_mi_1, *pX = xArr + Len_mi_1, *pZ = zArr + Len_mi_1;
	double RobsLoc = DistrInfoDat.yStart - sEnd;
	double InvRobsLoc = 1./RobsLoc;
	double MisFitXStLast = (DistrInfoDat.xStart - *pX)*InvRobsLoc - *pBtx;
	double MisFitXFiLast = (DistrInfoDat.xEnd - *pX)*InvRobsLoc - *pBtx;
	double MisFitZStLast = (DistrInfoDat.zStart - *pZ)*InvRobsLoc - *pBtz;
	double MisFitZFiLast = (DistrInfoDat.zEnd - *pZ)*InvRobsLoc - *pBtz;

	double midObsX = 0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd);
	double midObsZ = 0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd);
	double MisFitXMidLast = (midObsX - *pX)*InvRobsLoc - *pBtx;
	double MisFitZMidLast = (midObsZ - *pZ)*InvRobsLoc - *pBtz;

	const double VeryLarge = 1.E+23;
	double RobsXSt = VeryLarge, RobsXFi = VeryLarge, RobsZSt = VeryLarge, RobsZFi = VeryLarge;
	double RobsXMid = VeryLarge, RobsZMid = VeryLarge;

	//for(int is=1; is<TrjDat.LenFieldData; is++)
	//for(int is=1; is<NpVsS; is++)
	for(long long is=1; is<NpVsS; is++)
	{
		RobsLoc += TrjDat.sStep;
		InvRobsLoc = 1./RobsLoc;
		pBtx--; pBtz--; pX--; pZ--;

		if(RobsXSt == VeryLarge)
		{
			double MisFitXSt = (DistrInfoDat.xStart - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXSt*MisFitXStLast < 0.) RobsXSt = RobsLoc;
		}
		if(RobsXFi == VeryLarge)
		{
			double MisFitXFi = (DistrInfoDat.xEnd - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXFi*MisFitXFiLast < 0.) RobsXFi = RobsLoc;
		}
		if(RobsXMid == VeryLarge)
		{
			double MisFitXMid = (midObsX - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXMid*MisFitXMidLast < 0.) RobsXMid = RobsLoc;
		}

		if(RobsZSt == VeryLarge)
		{
			double MisFitZSt = (DistrInfoDat.zStart - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZSt*MisFitZStLast < 0.) RobsZSt = RobsLoc;
		}
		if(RobsZFi == VeryLarge)
		{
			double MisFitZFi = (DistrInfoDat.zEnd - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZFi*MisFitZFiLast < 0.) RobsZFi = RobsLoc;
		}
		if(RobsZMid == VeryLarge)
		{
			double MisFitZMid = (midObsZ - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZMid*MisFitZMidLast < 0.) RobsZMid = RobsLoc;
		}
	}
	double MinRobsX = (RobsXSt < RobsXFi)? RobsXSt : RobsXFi; //OC test roll back 130208
	//double MinRobsX = RobsXMid; //???
	double MinRobsZ = (RobsZSt < RobsZFi)? RobsZSt : RobsZFi;
	//double MinRobsZ = RobsZMid; //???

	double MinRobs = (MinRobsX < MinRobsZ)? MinRobsX : MinRobsZ;
	if(MinRobs != VeryLarge)
	{
		Robs = MinRobs;
		RobsAbsErr = 0.25*sRange;
        Ysrc = DistrInfoDat.yStart - Robs;
	}
	else
	{
		//if((TrjDat.sStart < 0.) && (sEnd > 0.)) Ysrc = 0.35*sRange;
		//else Ysrc = TrjDat.sStart + 0.75*sRange;

		if((sStart < 0.) && (sEnd > 0.)) Ysrc = 0.35*sRange;
		else Ysrc = sStart + 0.75*sRange;

		Robs = DistrInfoDat.yStart - Ysrc;
		RobsAbsErr = 0.25*sRange;
	}

	// Make more smart estimation of the "Source Point" later !!!
	//xElAtYsrc = TrjDat.EbmDat.x0; zElAtYsrc = TrjDat.EbmDat.z0;

	//int YsrcIndNo = (Ysrc - TrjDat.sStart)/TrjDat.sStep;
	//int YsrcIndNo = (int)((Ysrc - sStart)/TrjDat.sStep + 0.00001);
	long long YsrcIndNo = (long long)((Ysrc - sStart)/TrjDat.sStep + 0.00001);
	if(YsrcIndNo < 0) YsrcIndNo = 0;
	//if(YsrcIndNo >= TrjDat.LenFieldData) YsrcIndNo = TrjDat.LenFieldData - 1;
	if(YsrcIndNo >= NpVsS) YsrcIndNo = NpVsS - 1;

	xElAtYsrc = xArr[YsrcIndNo];
	zElAtYsrc = zArr[YsrcIndNo];
		
	if(TmpDataStorage != 0) delete[] TmpDataStorage;
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::FindAverageDistanceToSource(srTTrjDat& TrjDat, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc, double* precPar)
{// Should be called after the trajectory is already computed!!!
	
	double sStart = TrjDat.sStart;
	double sRange = (TrjDat.LenFieldData - 1)*TrjDat.sStep;
	double sEnd = sStart + sRange;
	//int NpVsS = TrjDat.LenFieldData;
	long long NpVsS = TrjDat.LenFieldData;

	if(precPar != 0)
	{
		//arPrecPar = array('d', [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, tuneNxNyForProp, sampFactNxNyForProp])
		double sStartIntPrec = precPar[2]; //pPrecElecFld->sStartInt;
		double sEndIntPrec = precPar[3]; //pPrecElecFld->sEndInt;

		double sStart0 = sStart, sEnd0 = sEnd;
		bool SpecLimitsMayBeDefined = (sStartIntPrec < sEndIntPrec);

		if(SpecLimitsMayBeDefined && (sStartIntPrec > sStart) && (sStartIntPrec < sEnd)) 
		{
			sStart = sStartIntPrec;
		}
		if(SpecLimitsMayBeDefined && (sEndIntPrec > sStart) && (sEndIntPrec < sEnd)) 
		{
			sEnd = sEndIntPrec;
		}
		if((sStart != sStart0) || (sEnd != sEnd0))
		{
			NpVsS = (int)((sEnd - sStart)/TrjDat.sStep + 0.00001) + 1;
			sEnd = sStart + (NpVsS - 1)*TrjDat.sStep;
			sRange = sEnd - sStart;
		}
	}

	double Ysrc = 0.;
	char DistUnits = 0; // X, Z in mm ??? - to check
	double *TmpDataStorage = new double[TrjDat.LenFieldData << 2];
	if(TmpDataStorage == 0) return MEMORY_ALLOCATION_FAILURE;
	double *BtxArr = TmpDataStorage;
	double *BtzArr = TmpDataStorage + TrjDat.LenFieldData;
	double *xArr = TmpDataStorage + (TrjDat.LenFieldData << 1);
	double *zArr = TmpDataStorage + (TrjDat.LenFieldData*3);
	TrjDat.CompTotalTrjDataTrjDisp(sStart, sEnd, NpVsS, BtxArr, BtzArr, xArr, zArr, DistUnits);

	//int Len_mi_1 = NpVsS - 1;
	long long Len_mi_1 = NpVsS - 1;
	double *pBtx = BtxArr + Len_mi_1, *pBtz = BtzArr + Len_mi_1, *pX = xArr + Len_mi_1, *pZ = zArr + Len_mi_1;

	double RobsLoc = yStart - sEnd;
	double InvRobsLoc = 1./RobsLoc;
	double MisFitXStLast = (xStart - *pX)*InvRobsLoc - *pBtx;
	double xEnd = xStart + xStep*(nx - 1);
	double zEnd = zStart + zStep*(nz - 1);
	double MisFitXFiLast = (xEnd - *pX)*InvRobsLoc - *pBtx;
	double MisFitZStLast = (zStart - *pZ)*InvRobsLoc - *pBtz;
	double MisFitZFiLast = (zEnd - *pZ)*InvRobsLoc - *pBtz;
	double midObsX = 0.5*(xStart + xEnd);
	double midObsZ = 0.5*(zStart + zEnd);
	//double MisFitXMidLast = (midObsX - *pX)*InvRobsLoc - *pBtx;
	//double MisFitZMidLast = (midObsZ - *pZ)*InvRobsLoc - *pBtz;

	const double VeryLarge = 1.E+23;
	double RobsXSt = VeryLarge, RobsXFi = VeryLarge, RobsZSt = VeryLarge, RobsZFi = VeryLarge;
	//double RobsXMid = VeryLarge, RobsZMid = VeryLarge;

	//for(int is=1; is<NpVsS; is++)
	for(long long is=1; is<NpVsS; is++)
	{
		RobsLoc += TrjDat.sStep;
		InvRobsLoc = 1./RobsLoc;
		pBtx--; pBtz--; pX--; pZ--;

		if(RobsXSt == VeryLarge)
		{
			double MisFitXSt = (xStart - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXSt*MisFitXStLast < 0.) RobsXSt = RobsLoc;
		}
		if(RobsXFi == VeryLarge)
		{
			double MisFitXFi = (xEnd - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXFi*MisFitXFiLast < 0.) RobsXFi = RobsLoc;
		}
		//if(RobsXMid == VeryLarge)
		//{
		//	double MisFitXMid = (midObsX - *pX)*InvRobsLoc - *pBtx;
		//	if(MisFitXMid*MisFitXMidLast < 0.) RobsXMid = RobsLoc;
		//}
		if(RobsZSt == VeryLarge)
		{
			double MisFitZSt = (zStart - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZSt*MisFitZStLast < 0.) RobsZSt = RobsLoc;
		}
		if(RobsZFi == VeryLarge)
		{
			double MisFitZFi = (zEnd - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZFi*MisFitZFiLast < 0.) RobsZFi = RobsLoc;
		}
		//if(RobsZMid == VeryLarge)
		//{
		//	double MisFitZMid = (midObsZ - *pZ)*InvRobsLoc - *pBtz;
		//	if(MisFitZMid*MisFitZMidLast < 0.) RobsZMid = RobsLoc;
		//}

		if((RobsXSt != VeryLarge) && (RobsXFi != VeryLarge) && (RobsZSt != VeryLarge) && (RobsZFi != VeryLarge)) break; //OC190414
	}

	double MinRobsX = (RobsXSt < RobsXFi)? RobsXSt : RobsXFi; //OC test roll back 130208
	double MinRobsZ = (RobsZSt < RobsZFi)? RobsZSt : RobsZFi;
	double MinRobs = (MinRobsX < MinRobsZ)? MinRobsX : MinRobsZ;

	//Estimating MaxRobs OC190414
	pBtx = BtxArr; pBtz = BtzArr; pX = xArr; pZ = zArr;
	RobsLoc = yStart - sStart;
	InvRobsLoc = 1./RobsLoc;
	MisFitXStLast = (xStart - *pX)*InvRobsLoc - *pBtx;
	MisFitXFiLast = (xEnd - *pX)*InvRobsLoc - *pBtx;
	MisFitZStLast = (zStart - *pZ)*InvRobsLoc - *pBtz;
	MisFitZFiLast = (zEnd - *pZ)*InvRobsLoc - *pBtz;
	RobsXSt = VeryLarge; RobsXFi = VeryLarge; RobsZSt = VeryLarge; RobsZFi = VeryLarge;

	//for(int is=1; is<NpVsS; is++)
	for(long long is=1; is<NpVsS; is++)
	{
		RobsLoc -= TrjDat.sStep;
		InvRobsLoc = 1./RobsLoc;
		pBtx++; pBtz++; pX++; pZ++;

		if(RobsXSt == VeryLarge)
		{
			double MisFitXSt = (xStart - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXSt*MisFitXStLast < 0.) RobsXSt = RobsLoc;
		}
		if(RobsXFi == VeryLarge)
		{
			double MisFitXFi = (xEnd - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXFi*MisFitXFiLast < 0.) RobsXFi = RobsLoc;
		}
		if(RobsZSt == VeryLarge)
		{
			double MisFitZSt = (zStart - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZSt*MisFitZStLast < 0.) RobsZSt = RobsLoc;
		}
		if(RobsZFi == VeryLarge)
		{
			double MisFitZFi = (zEnd - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZFi*MisFitZFiLast < 0.) RobsZFi = RobsLoc;
		}
	}

	double MaxRobsX = (RobsXSt < RobsXFi)? RobsXSt : RobsXFi; 
	double MaxRobsZ = (RobsZSt < RobsZFi)? RobsZSt : RobsZFi;
	double MaxRobs = (MinRobsX < MinRobsZ)? MaxRobsX : MaxRobsZ;

	double Ravg = 0.5*(MinRobs + MaxRobs); //OC190414
	if(MinRobs == VeryLarge)
	{
		if(MaxRobs == VeryLarge) Ravg = VeryLarge;
		else Ravg = MaxRobs;
	}
	else
	{
		if(MaxRobs == VeryLarge) Ravg = MinRobs;
	}

	//if(MinRobs != VeryLarge)
	if(Ravg != VeryLarge)
	{
		Robs = Ravg;
		RobsAbsErr = 0.25*sRange;
        Ysrc = yStart - Robs;
	}
	else
	{
		if((sStart < 0.) && (sEnd > 0.)) Ysrc = 0.35*sRange;
		else Ysrc = sStart + 0.75*sRange;

		Robs = yStart - Ysrc;
		RobsAbsErr = 0.25*sRange;
	}

	//int YsrcIndNo = (int)((Ysrc - sStart)/TrjDat.sStep + 0.00001);
	long long YsrcIndNo = (long long)((Ysrc - sStart)/TrjDat.sStep + 0.00001);
	if(YsrcIndNo < 0) YsrcIndNo = 0;

	if(YsrcIndNo >= NpVsS) YsrcIndNo = NpVsS - 1;

	xElAtYsrc = xArr[YsrcIndNo];
	zElAtYsrc = zArr[YsrcIndNo];
		
	if(TmpDataStorage != 0) delete[] TmpDataStorage;
	return 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::AddStokesAtPoint(srTEXZ& EXZ, float* pStokesVal)
{
	double x = EXZ.x, z = EXZ.z, e = EXZ.e;
	double eFinMin = eStart + (ne - 1)*eStep;
	double xFinMin = xStart + (nx - 1)*xStep;
	double zFinMin = zStart + (nz - 1)*zStep;

	const double RelEqStepTol = 0.1;
	double AbsEqStepTolE = RelEqStepTol*eStep;
	double AbsEqStepTolX = RelEqStepTol*xStep;
	double AbsEqStepTolZ = RelEqStepTol*zStep;

	if((e < eStart - AbsEqStepTolE) || (e > eFinMin + AbsEqStepTolE)) return;
	if((x < xStart - AbsEqStepTolX) || (x > xFinMin + AbsEqStepTolX)) return;
	if((z < zStart - AbsEqStepTolZ) || (z > zFinMin + AbsEqStepTolZ)) return;

	long ie = 0;
	if(ne > 1)
	{
		ie = (long)((e - eStart)/eStep);
		if(ie < 0) ie = 0;
		else if(ie >= ne) ie = ne - 1;
	}
	long ixMin = 0, ixMax = 0;
	if(nx > 1)
	{
		ixMin = (long)((x - xStart)/xStep);
        ixMax = ixMin + 1;
		long nx_mi_1 = nx - 1;
		if(ixMin < 0) { ixMin = ixMax = 0;}
		else if(ixMin >= nx_mi_1) { ixMin = nx_mi_1; ixMax = ixMin;}
	}
	long izMin = 0, izMax = 0;
	if(nz > 1)
	{
		izMin = (long)((z - zStart)/zStep);
        izMax = izMin + 1;
		long nz_mi_1 = nz - 1;
		if(izMin < 0) { izMin = izMax = 0;}
		else if(izMin >= nz_mi_1) { izMin = nz_mi_1; izMax = izMin;}
	}

	double xr = 0, zr = 0;
	if(ixMax != ixMin) xr = (x - (xStart + xStep*ixMin))/xStep;
	if(xr < 0) xr = 0;
	else if(xr > 1) xr = 1;
	if(izMax != izMin) zr = (z - (zStart + zStep*izMin))/zStep;
	if(zr < 0) zr = 0;
	else if(zr > 1) zr = 1;

	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	//long TwoIe = ie << 1;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long TwoIe = ie << 1;

	//long ixMinPerX = ixMin*PerX, ixMaxPerX = ixMax*PerX;
	//long izMinPerZ = izMin*PerZ, izMaxPerZ = izMax*PerZ;
	long long ixMinPerX = ixMin*PerX, ixMaxPerX = ixMax*PerX;
	long long izMinPerZ = izMin*PerZ, izMaxPerZ = izMax*PerZ;

	//long Offset00 = izMinPerZ + ixMinPerX + TwoIe;
	//long Offset10 = izMinPerZ + ixMaxPerX + TwoIe;
	//long Offset01 = izMaxPerZ + ixMinPerX + TwoIe;
	//long Offset11 = izMaxPerZ + ixMaxPerX + TwoIe;
	long long Offset00 = izMinPerZ + ixMinPerX + TwoIe;
	long long Offset10 = izMinPerZ + ixMaxPerX + TwoIe;
	long long Offset01 = izMaxPerZ + ixMinPerX + TwoIe;
	long long Offset11 = izMaxPerZ + ixMaxPerX + TwoIe;

	float ZeroArr[] = {0, 0};
	float *pEx00, *pEz00, *pEx10, *pEz10, *pEx01, *pEz01, *pEx11, *pEz11;
	if(pBaseRadX != 0)
	{
        pEx00 = pBaseRadX + Offset00; pEx10 = pBaseRadX + Offset10; pEx01 = pBaseRadX + Offset01; pEx11 = pBaseRadX + Offset11;
	}
	else pEx00 = pEx10 = pEx01 = pEx11 = ZeroArr;
	if(pBaseRadZ != 0)
	{
        pEz00 = pBaseRadZ + Offset00; pEz10 = pBaseRadZ + Offset10; pEz01 = pBaseRadZ + Offset01; pEz11 = pBaseRadZ + Offset11;
	}
	else pEz00 = pEz10 = pEz01 = pEz11 = ZeroArr;

    float pSt00[4], pSt10[4], pSt01[4], pSt11[4];
	CalcStokesFromE(pEx00, pEz00, pSt00);
	CalcStokesFromE(pEx10, pEz10, pSt10);
	CalcStokesFromE(pEx01, pEz01, pSt01);
	CalcStokesFromE(pEx11, pEz11, pSt11);
    float *tSt00 = pSt00, *tSt10 = pSt10, *tSt01 = pSt01, *tSt11 = pSt11, *tStokesVal = pStokesVal;
	//for(int i=0; i<4; i++) *(tStokesVal++) += (float)TMathInterpolMultiD::Interp2D4pRel(*(tSt00++), *(tSt10++), *(tSt01++), *(tSt11++), xr, zr);
	for(int i=0; i<4; i++) *(tStokesVal++) += (float)CGenMathInterp::Interp2D4pRel(*(tSt00++), *(tSt10++), *(tSt01++), *(tSt11++), xr, zr);
}

//*************************************************************************

void srTSRWRadStructAccessData::CheckAndSubtractPhaseTermsLin(double newXc, double newZc)
{
	const double ratAllowSubtract = 0.2;
	const double minNumOptCycles = 10;
	double lambda_m = 3.1415926535898/(eStart*2.53384080189E+06);
	const double twoPi = 2.*3.1415926535898;
	
	bool xLinPhaseTermCanBeSubtracted = false;
	double dxcSubNew = newXc - xc;
	if(RobsX != 0)
	{
		double xWfrRange = xStep*(nx - 1);
		double xNumOptCycles = 0.25*xWfrRange*xWfrRange/(lambda_m*RobsX);
		xLinPhaseTermCanBeSubtracted = (fabs(RobsX) > fabs(RobsXAbsErr)) && (xNumOptCycles > minNumOptCycles) && (fabs(dxcSubNew)/xWfrRange > ratAllowSubtract);
	}
	if(!xLinPhaseTermCanBeSubtracted) dxcSubNew = 0;
	
	bool zLinPhaseTermCanBeSubtracted = false;
	double dzcSubNew = newZc - zc;
	if(RobsZ != 0)
	{
		double zWfrRange = zStep*(nz - 1);
		double zNumOptCycles = 0.25*zWfrRange*zWfrRange/(lambda_m*RobsZ);
		zLinPhaseTermCanBeSubtracted = (fabs(RobsZ) > fabs(RobsZAbsErr)) && (zNumOptCycles > minNumOptCycles) && (fabs(dzcSubNew)/zWfrRange > ratAllowSubtract);
	}
	if(!zLinPhaseTermCanBeSubtracted) dzcSubNew = 0;

	if((!xLinPhaseTermCanBeSubtracted) && (!zLinPhaseTermCanBeSubtracted)) return;

	double dxc = dxcSubNew;
	if(m_xLinOnlyPhaseTermWasSubtracted)
	{
		dxc -= m_dxcSub;
		if(fabs(dxc)/fabs(m_dxcSub) < ratAllowSubtract) 
		{
			dxc = 0; dxcSubNew = m_dxcSub;
		}
	}
	m_dxcSub = dxcSubNew;
	if(dxcSubNew == 0) m_xLinOnlyPhaseTermWasSubtracted = false;
	else m_xLinOnlyPhaseTermWasSubtracted = true;
	double xMult = -twoPi*dxc/(lambda_m*RobsX);

	double dzc = dzcSubNew;
	if(zLinPhaseTermCanBeSubtracted && m_zLinOnlyPhaseTermWasSubtracted)
	{
		dzc -= m_dzcSub;
		if(fabs(dzc)/fabs(m_dzcSub) < ratAllowSubtract) 
		{
			dzc = 0; dzcSubNew = m_dzcSub;
		}
	}
	m_dzcSub = dzcSubNew;
	if(dzcSubNew == 0) m_zLinOnlyPhaseTermWasSubtracted = false;
	else m_zLinOnlyPhaseTermWasSubtracted = true;
	double zMult = -twoPi*dzc/(lambda_m*RobsZ);

	if((xMult == 0) && (zMult == 0)) return;
	
	MultiplyElFieldByPhaseLin(xMult, zMult);
}

//*************************************************************************

void srTSRWRadStructAccessData::CheckAndResetPhaseTermsLin()
{
	if((!m_xLinOnlyPhaseTermWasSubtracted) && (!m_zLinOnlyPhaseTermWasSubtracted)) return;
	
	double lambda_m = 3.1415926535898/(eStart*2.53384080189E+06);
	const double twoPi = 2.*3.1415926535898;
	
	double xMult = 0, zMult = 0;
	if(m_xLinOnlyPhaseTermWasSubtracted && (m_dxcSub != 0) && (RobsX != 0))
	{
		xMult = twoPi*m_dxcSub/(lambda_m*RobsX);
	}
	if(m_zLinOnlyPhaseTermWasSubtracted && (m_dzcSub != 0) && (RobsZ != 0))
	{
		zMult = twoPi*m_dzcSub/(lambda_m*RobsZ);
	}
	
	m_xLinOnlyPhaseTermWasSubtracted = false;
	m_zLinOnlyPhaseTermWasSubtracted = false;
	m_dxcSub = 0;
	m_dzcSub = 0;
	
	if((xMult == 0) && (zMult == 0)) return;
	
	MultiplyElFieldByPhaseLin(xMult, zMult);
}

//*************************************************************************

void srTSRWRadStructAccessData::MirrorFieldData(int sx, int sz)
{// sx < 0 means mirroring should be done vs x 
 // sz < 0 means mirroring should be done vs z 
	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	float buf;
	float *pEX0 = pBaseRadX;
	float *pEZ0 = pBaseRadZ;

	if((sx > 0) && (sz > 0)) return; //no mirroring is necessary 
	else if((sx < 0) && (sz > 0)) //mirroring with respect to x
	{
		long nx_mi_1 = nx - 1;
		for(long ie=0; ie<ne; ie++)
		{
			long Two_ie = ie << 1;
			for(long iz=0; iz<nz; iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				float *pEX_StartForX = pEX0 + izPerZ;
				float *pEZ_StartForX = pEZ0 + izPerZ;

				for(long ix=0; ix<(nx >> 1); ix++)
				{
					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					float *pEX = pEX_StartForX + ixPerX_p_Two_ie;
					float *pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

					//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					float *rev_pEX = pEX_StartForX + rev_ixPerX_p_Two_ie;
					float *rev_pEZ = pEZ_StartForX + rev_ixPerX_p_Two_ie;

					if(pEX0 != 0)
					{
						buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
						buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
					}
					if(pEZ0 != 0)
					{
						buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
						buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
					}
				}
			}
		}
	}
	else if((sx > 0) && (sz < 0))
	{
		long nz_mi_1 = nz - 1;
		for(long ie=0; ie<ne; ie++)
		{
			long Two_ie = ie << 1;
			for(long iz=0; iz<(nz >> 1); iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				float *pEX_StartForX = pEX0 + izPerZ;
				float *pEZ_StartForX = pEZ0 + izPerZ;

				//long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				long long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				float *rev_pEX_StartForX = pEX0 + rev_izPerZ;
				float *rev_pEZ_StartForX = pEZ0 + rev_izPerZ;

				for(long ix=0; ix<nx; ix++)
				{
					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					float *pEX = pEX_StartForX + ixPerX_p_Two_ie;
					float *pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

					float *rev_pEX = rev_pEX_StartForX + ixPerX_p_Two_ie;
					float *rev_pEZ = rev_pEZ_StartForX + ixPerX_p_Two_ie;

					if(pEX0 != 0)
					{
						buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
						buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
					}
					if(pEZ0 != 0)
					{
						buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
						buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
					}
				}
			}
		}
	}
	else
	{
		long nx_mi_1 = nx - 1;
		long nz_mi_1 = nz - 1;
		for(long ie=0; ie<ne; ie++)
		{
			long Two_ie = ie << 1;
			for(long iz=0; iz<(nz >> 1); iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				float *pEX_StartForX = pEX0 + izPerZ;
				float *pEZ_StartForX = pEZ0 + izPerZ;

				//long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				long long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				float *rev_pEX_StartForX = pEX0 + rev_izPerZ;
				float *rev_pEZ_StartForX = pEZ0 + rev_izPerZ;

				for(long ix=0; ix<nx; ix++)
				{
					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					float *pEX = pEX_StartForX + ixPerX_p_Two_ie;
					float *pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

					//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					float *rev_pEX = rev_pEX_StartForX + rev_ixPerX_p_Two_ie;
					float *rev_pEZ = rev_pEZ_StartForX + rev_ixPerX_p_Two_ie;

					if(pEX0 != 0)
					{
						buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
						buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
					}
					if(pEZ0 != 0)
					{
						buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
						buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
					}
				}
			}
			if(((nz >> 1) << 1) != nz)
			{
				//long izPerZ = ((nz >> 1) + 1)*PerZ;
				long long izPerZ = ((nz >> 1) + 1)*PerZ;
				float *pEX_StartForX = pEX0 + izPerZ;
				float *pEZ_StartForX = pEZ0 + izPerZ;

				for(long ix=0; ix<(nx >> 1); ix++)
				{
					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					float *pEX = pEX_StartForX + ixPerX_p_Two_ie;
					float *pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

					//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					float *rev_pEX = pEX_StartForX + rev_ixPerX_p_Two_ie;
					float *rev_pEZ = pEZ_StartForX + rev_ixPerX_p_Two_ie;

					if(pEX0 != 0)
					{
						buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
						buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
					}
					if(pEZ0 != 0)
					{
						buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
						buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
					}
				}
			}
		}
	}
}
//*************************************************************************

int srTSRWRadStructAccessData::ExtractSliceConstEorT(long ie, float*& pOutEx, float*& pOutEz)
{// ATTENTION: In the case of single energy, it simply returns pointers to pBaseRadX, pBaseRadZ!!!

	float *pEx0 = pBaseRadX;
	float *pEz0 = pBaseRadZ;

	if(ne == 1)
	{
		pOutEx = pEx0; pOutEz = pEz0;
		return 0;
	}

	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;

	//long izPerZ = 0;
	//long iePerE = ie << 1;
	long long izPerZ = 0;
	long long iePerE = ie << 1;

	float *tOutEx = pOutEx, *tOutEz = pOutEz;
	for(int iz=0; iz<nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pEx = pEx_StartForX + ixPerX_p_iePerE;
			float *pEz = pEz_StartForX + ixPerX_p_iePerE;

			*(tOutEx++) = *(pEx++); *(tOutEx++) = *pEx;
			*(tOutEz++) = *(pEz++); *(tOutEz++) = *pEz;

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::SetupSliceConstEorT(long ie, float* pInEx, float* pInEz)
{
	float *pEx0 = pBaseRadX;
	float *pEz0 = pBaseRadZ;
	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;

	//long izPerZ = 0;
	//long iePerE = ie << 1;
	long long izPerZ = 0;
	long long iePerE = ie << 1;

	float *tInEx = pInEx, *tInEz = pInEz;
	for(int iz=0; iz<nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pEx = pEx_StartForX + ixPerX_p_iePerE;
			float *pEz = pEz_StartForX + ixPerX_p_iePerE;

			*(pEx++) = *(tInEx++); *pEx = *(tInEx++);
			*(pEz++) = *(tInEz++); *pEz = *(tInEz++);

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
	return 0;
}

//*************************************************************************

int srTSRWRadStructAccessData::SetupWfrEdgeCorrData(float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrsForWfrEdgeCorr)
{//Copied from srTGenOptElem::
	int result;

	double xAbsTol = 0.05*xStep;
	double zAbsTol = 0.05*zStep;
//--X
	double xWfrMinOffsetFromStart = xWfrMin - xStart;
	long ixWfrMinLower = long(xWfrMinOffsetFromStart/xStep + 1.E-13);
	double xWfrMinLowerMisfit = xWfrMinOffsetFromStart - ixWfrMinLower*xStep;

	double xWfrMaxOffsetFromStart = xWfrMax - xStart;
	long ixWfrMaxLower = long(xWfrMaxOffsetFromStart/xStep + 1.E-13);
	double xWfrMaxLowerMisfit = xWfrMaxOffsetFromStart - ixWfrMaxLower*xStep;
	
	char xWfrMinIsBetweenMeshPoints = (xWfrMinLowerMisfit > xAbsTol);
	char xWfrMaxIsBetweenMeshPoints = (xWfrMaxLowerMisfit > xAbsTol);
	char xWfrMaxIsSmallerThanDataEnd = (::fabs((xStart + nx*xStep) - xWfrMax) > xAbsTol);
	char xWfrCorrNeeded = (xWfrMinIsBetweenMeshPoints || xWfrMaxIsBetweenMeshPoints || xWfrMaxIsSmallerThanDataEnd);

	float dxSt = 0.;
	if(xWfrMinIsBetweenMeshPoints) dxSt = (float)(xStep - xWfrMinLowerMisfit);

	float dxFi = 0.;
	if(xWfrMaxIsBetweenMeshPoints) dxFi = (float)(xStep - xWfrMaxLowerMisfit);
	else if(xWfrMaxIsSmallerThanDataEnd) dxFi = (float)(xStep);
//--Z
	double zWfrMinOffsetFromStart = zWfrMin - zStart;

	long izWfrMinLower = long(zWfrMinOffsetFromStart/zStep + 1.E-13);
	double zWfrMinLowerMisfit = zWfrMinOffsetFromStart - izWfrMinLower*zStep;

	double zWfrMaxOffsetFromStart = zWfrMax - zStart;
	long izWfrMaxLower = long(zWfrMaxOffsetFromStart/zStep + 1.E-13);
	double zWfrMaxLowerMisfit = zWfrMaxOffsetFromStart - izWfrMaxLower*zStep;
	
	char zWfrMinIsBetweenMeshPoints = (zWfrMinLowerMisfit > zAbsTol);
	char zWfrMaxIsBetweenMeshPoints = (zWfrMaxLowerMisfit > zAbsTol);
	char zWfrMaxIsSmallerThanDataEnd = (::fabs((zStart + nz*zStep) - zWfrMax) > zAbsTol);
	char zWfrCorrNeeded = (zWfrMinIsBetweenMeshPoints || zWfrMaxIsBetweenMeshPoints || zWfrMaxIsSmallerThanDataEnd);

	float dzSt = 0.;
	if(zWfrMinIsBetweenMeshPoints) dzSt = (float)(zStep - zWfrMinLowerMisfit);

	float dzFi = 0.;
	if(zWfrMaxIsBetweenMeshPoints) dzFi = (float)(zStep - zWfrMaxLowerMisfit);
	else if(zWfrMaxIsSmallerThanDataEnd) dzFi = (float)(zStep);

//--Gen
	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = xStep;
	FFT2DInfo.yStep = zStep;
	FFT2DInfo.xStart = xStart;
	FFT2DInfo.yStart = zStart;
	FFT2DInfo.Nx = nx;
	FFT2DInfo.Ny = nz;
	FFT2DInfo.UseGivenStartTrValues = 0;
	CGenMathFFT2D FFT2D;
	FFT2D.SetupLimitsTr(FFT2DInfo);

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.HowMany = 2;
	FFT1DInfo.UseGivenStartTrValue = 0;

	if(xWfrCorrNeeded || zWfrCorrNeeded)
	{
		DataPtrsForWfrEdgeCorr.dx = xStep;
		DataPtrsForWfrEdgeCorr.dz = zStep;

		long TwoNx = nx << 1;
		long TwoNz = nz << 1;

		if(dxSt != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrXSt = new float[TwoNx];
			if(DataPtrsForWfrEdgeCorr.ExpArrXSt == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrXStEx = new float[TwoNz << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrXStEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrXStEz = DataPtrsForWfrEdgeCorr.FFTArrXStEx + TwoNz;
			DataPtrsForWfrEdgeCorr.dxSt = dxSt;

			long jxSt = ixWfrMinLower + 1;
			double xjSt = xStart + jxSt*xStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrXSt, nx, xjSt, FFT2DInfo.xStartTr, FFT2DInfo.xStepTr);
			SetupRadXorZSectFromSliceConstEorT(pDataEx, pDataEz, nx, nz, 'z', jxSt, DataPtrsForWfrEdgeCorr.FFTArrXStEx, DataPtrsForWfrEdgeCorr.FFTArrXStEz);

			if(dzSt != 0.)
			{
				long jzSt2 = (izWfrMinLower + 1) << 1;
				DataPtrsForWfrEdgeCorr.fxStzSt[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzSt2);
				DataPtrsForWfrEdgeCorr.fxStzSt[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzSt2 + 1);
				DataPtrsForWfrEdgeCorr.fxStzSt[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzSt2);
				DataPtrsForWfrEdgeCorr.fxStzSt[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzSt2 + 1);
			}
			if(dzFi != 0.)
			{
				long jzFi2 = izWfrMaxLower << 1;
				DataPtrsForWfrEdgeCorr.fxStzFi[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzFi2);
				DataPtrsForWfrEdgeCorr.fxStzFi[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzFi2 + 1);
				DataPtrsForWfrEdgeCorr.fxStzFi[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzFi2);
				DataPtrsForWfrEdgeCorr.fxStzFi[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzFi2 + 1);
			}

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrXStEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = zStep;
			FFT1DInfo.xStart = zStart;
			FFT1DInfo.Nx = nz;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dxFi != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrXFi = new float[TwoNx];
			if(DataPtrsForWfrEdgeCorr.ExpArrXFi == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrXFiEx = new float[TwoNz << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrXFiEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrXFiEz = DataPtrsForWfrEdgeCorr.FFTArrXFiEx + TwoNz;
			DataPtrsForWfrEdgeCorr.dxFi = dxFi;

			double xjFi = xStart + ixWfrMaxLower*xStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrXFi, nx, xjFi, FFT2DInfo.xStartTr, FFT2DInfo.xStepTr);
			SetupRadXorZSectFromSliceConstEorT(pDataEx, pDataEz, nx, nz, 'z', ixWfrMaxLower, DataPtrsForWfrEdgeCorr.FFTArrXFiEx, DataPtrsForWfrEdgeCorr.FFTArrXFiEz);

			if(dzSt != 0.)
			{
				long jzSt2 = (izWfrMinLower + 1) << 1;
				DataPtrsForWfrEdgeCorr.fxFizSt[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzSt2);
				DataPtrsForWfrEdgeCorr.fxFizSt[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzSt2 + 1);
				DataPtrsForWfrEdgeCorr.fxFizSt[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzSt2);
				DataPtrsForWfrEdgeCorr.fxFizSt[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzSt2 + 1);
			}
			if(dzFi != 0.)
			{
				long jzFi2 = izWfrMaxLower << 1;
				DataPtrsForWfrEdgeCorr.fxFizFi[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzFi2);
				DataPtrsForWfrEdgeCorr.fxFizFi[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzFi2 + 1);
				DataPtrsForWfrEdgeCorr.fxFizFi[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzFi2);
				DataPtrsForWfrEdgeCorr.fxFizFi[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzFi2 + 1);
			}

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrXFiEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = zStep;
			FFT1DInfo.xStart = zStart;
			FFT1DInfo.Nx = nz;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dzSt != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrZSt = new float[TwoNz];
			if(DataPtrsForWfrEdgeCorr.ExpArrZSt == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrZStEx = new float[TwoNx << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrZStEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrZStEz = DataPtrsForWfrEdgeCorr.FFTArrZStEx + TwoNx;
			DataPtrsForWfrEdgeCorr.dzSt = dzSt;

			long jzSt = izWfrMinLower + 1;
			double zjSt = zStart + jzSt*zStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrZSt, nz, zjSt, FFT2DInfo.yStartTr, FFT2DInfo.yStepTr);
			SetupRadXorZSectFromSliceConstEorT(pDataEx, pDataEz, nx, nz, 'x', jzSt, DataPtrsForWfrEdgeCorr.FFTArrZStEx, DataPtrsForWfrEdgeCorr.FFTArrZStEz);

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrZStEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = xStep;
			FFT1DInfo.xStart = xStart;
			FFT1DInfo.Nx = nx;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dzFi != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrZFi = new float[TwoNz];
			if(DataPtrsForWfrEdgeCorr.ExpArrZFi == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrZFiEx = new float[TwoNx << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrZFiEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrZFiEz = DataPtrsForWfrEdgeCorr.FFTArrZFiEx + TwoNx;
			DataPtrsForWfrEdgeCorr.dzFi = dzFi;

			double zjFi = zStart + izWfrMaxLower*zStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrZFi, nz, zjFi, FFT2DInfo.yStartTr, FFT2DInfo.yStepTr);
			SetupRadXorZSectFromSliceConstEorT(pDataEx, pDataEz, nx, nz, 'x', izWfrMaxLower, DataPtrsForWfrEdgeCorr.FFTArrZFiEx, DataPtrsForWfrEdgeCorr.FFTArrZFiEz);

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrZFiEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = xStep;
			FFT1DInfo.xStart = xStart;
			FFT1DInfo.Nx = nx;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		DataPtrsForWfrEdgeCorr.WasSetup = 1;
	}
	return 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::MakeWfrEdgeCorrection(float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrs)
{//Copied from srTGenOptElem::
	float *tEx = pDataEx, *tEz = pDataEz;

	double dxSt_dzSt = DataPtrs.dxSt*DataPtrs.dzSt;
	double dxSt_dzFi = DataPtrs.dxSt*DataPtrs.dzFi;
	double dxFi_dzSt = DataPtrs.dxFi*DataPtrs.dzSt;
	double dxFi_dzFi = DataPtrs.dxFi*DataPtrs.dzFi;
	float fSSExRe = DataPtrs.fxStzSt[0], fSSExIm = DataPtrs.fxStzSt[1], fSSEzRe = DataPtrs.fxStzSt[2], fSSEzIm = DataPtrs.fxStzSt[3];
	float fFSExRe = DataPtrs.fxFizSt[0], fFSExIm = DataPtrs.fxFizSt[1], fFSEzRe = DataPtrs.fxFizSt[2], fFSEzIm = DataPtrs.fxFizSt[3];
	float fSFExRe = DataPtrs.fxStzFi[0], fSFExIm = DataPtrs.fxStzFi[1], fSFEzRe = DataPtrs.fxStzFi[2], fSFEzIm = DataPtrs.fxStzFi[3];
	float fFFExRe = DataPtrs.fxFizFi[0], fFFExIm = DataPtrs.fxFizFi[1], fFFEzRe = DataPtrs.fxFizFi[2], fFFEzIm = DataPtrs.fxFizFi[3];

	float bRe, bIm, cRe, cIm;
	for(long iz=0; iz<nz; iz++)
	{
		long Two_iz = iz << 1;
		long Two_iz_p_1 = Two_iz + 1;

		for(long ix=0; ix<nx; ix++)
		{
			long Two_ix = ix << 1;
			long Two_ix_p_1 = Two_ix + 1;

			float ExRe = *tEx, ExIm = *(tEx+1);
			float EzRe = *tEz, EzIm = *(tEz+1);

			if(DataPtrs.dxSt != 0.)
			{
				float ExpXStRe = DataPtrs.ExpArrXSt[Two_ix], ExpXStIm = DataPtrs.ExpArrXSt[Two_ix_p_1];

				bRe = DataPtrs.FFTArrXStEx[Two_iz]; bIm = DataPtrs.FFTArrXStEx[Two_iz_p_1];
				ExRe += (float)(DataPtrs.dxSt*(ExpXStRe*bRe - ExpXStIm*bIm));
				ExIm += (float)(DataPtrs.dxSt*(ExpXStRe*bIm + ExpXStIm*bRe));

				bRe = DataPtrs.FFTArrXStEz[Two_iz]; bIm = DataPtrs.FFTArrXStEz[Two_iz_p_1];
				EzRe += (float)(DataPtrs.dxSt*(ExpXStRe*bRe - ExpXStIm*bIm));
				EzIm += (float)(DataPtrs.dxSt*(ExpXStRe*bIm + ExpXStIm*bRe));

				if(DataPtrs.dzSt != 0.)
				{
					bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
					cRe = ExpXStRe*bRe - ExpXStIm*bIm; cIm = ExpXStRe*bIm + ExpXStIm*bRe;

					ExRe += (float)(dxSt_dzSt*(fSSExRe*cRe - fSSExIm*cIm));
					ExIm += (float)(dxSt_dzSt*(fSSExRe*cIm + fSSExIm*cRe));
					EzRe += (float)(dxSt_dzSt*(fSSEzRe*cRe - fSSEzIm*cIm));
					EzIm += (float)(dxSt_dzSt*(fSSEzRe*cIm + fSSEzIm*cRe));
				}
				if(DataPtrs.dzFi != 0.)
				{
					bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
					cRe = ExpXStRe*bRe - ExpXStIm*bIm; cIm = ExpXStRe*bIm + ExpXStIm*bRe;

					ExRe -= (float)(dxSt_dzFi*(fSFExRe*cRe - fSFExIm*cIm));
					ExIm -= (float)(dxSt_dzFi*(fSFExRe*cIm + fSFExIm*cRe));
					EzRe -= (float)(dxSt_dzFi*(fSFEzRe*cRe - fSFEzIm*cIm));
					EzIm -= (float)(dxSt_dzFi*(fSFEzRe*cIm + fSFEzIm*cRe));
				}
			}
			if(DataPtrs.dxFi != 0.)
			{
				float ExpXFiRe = DataPtrs.ExpArrXFi[Two_ix], ExpXFiIm = DataPtrs.ExpArrXFi[Two_ix_p_1];

				bRe = DataPtrs.FFTArrXFiEx[Two_iz]; bIm = DataPtrs.FFTArrXFiEx[Two_iz_p_1];
				ExRe -= (float)(DataPtrs.dxFi*(ExpXFiRe*bRe - ExpXFiIm*bIm));
				ExIm -= (float)(DataPtrs.dxFi*(ExpXFiRe*bIm + ExpXFiIm*bRe));

				bRe = DataPtrs.FFTArrXFiEz[Two_iz]; bIm = DataPtrs.FFTArrXFiEz[Two_iz_p_1];
				EzRe -= (float)(DataPtrs.dxFi*(ExpXFiRe*bRe - ExpXFiIm*bIm));
				EzIm -= (float)(DataPtrs.dxFi*(ExpXFiRe*bIm + ExpXFiIm*bRe));

				if(DataPtrs.dzSt != 0.)
				{
					bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
					cRe = ExpXFiRe*bRe - ExpXFiIm*bIm; cIm = ExpXFiRe*bIm + ExpXFiIm*bRe;

					ExRe -= (float)(dxFi_dzSt*(fFSExRe*cRe - fFSExIm*cIm));
					ExIm -= (float)(dxFi_dzSt*(fFSExRe*cIm + fFSExIm*cRe));
					EzRe -= (float)(dxFi_dzSt*(fFSEzRe*cRe - fFSEzIm*cIm));
					EzIm -= (float)(dxFi_dzSt*(fFSEzRe*cIm + fFSEzIm*cRe));
				}
				if(DataPtrs.dzFi != 0.)
				{
					bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
					cRe = ExpXFiRe*bRe - ExpXFiIm*bIm; cIm = ExpXFiRe*bIm + ExpXFiIm*bRe;

					ExRe += (float)(dxFi_dzFi*(fFFExRe*cRe - fFFExIm*cIm));
					ExIm += (float)(dxFi_dzFi*(fFFExRe*cIm + fFFExIm*cRe));
					EzRe += (float)(dxFi_dzFi*(fFFEzRe*cRe - fFFEzIm*cIm));
					EzIm += (float)(dxFi_dzFi*(fFFEzRe*cIm + fFFEzIm*cRe));
				}
			}
			if(DataPtrs.dzSt != 0.)
			{
				float ExpZStRe = DataPtrs.ExpArrZSt[Two_iz], ExpZStIm = DataPtrs.ExpArrZSt[Two_iz_p_1];

				bRe = DataPtrs.FFTArrZStEx[Two_ix]; bIm = DataPtrs.FFTArrZStEx[Two_ix_p_1];
				ExRe += (float)(DataPtrs.dzSt*(ExpZStRe*bRe - ExpZStIm*bIm));
				ExIm += (float)(DataPtrs.dzSt*(ExpZStRe*bIm + ExpZStIm*bRe));

				bRe = DataPtrs.FFTArrZStEz[Two_ix]; bIm = DataPtrs.FFTArrZStEz[Two_ix_p_1];
				EzRe += (float)(DataPtrs.dzSt*(ExpZStRe*bRe - ExpZStIm*bIm));
				EzIm += (float)(DataPtrs.dzSt*(ExpZStRe*bIm + ExpZStIm*bRe));
			}
			if(DataPtrs.dzFi != 0.)
			{
				float ExpZFiRe = DataPtrs.ExpArrZFi[Two_iz], ExpZFiIm = DataPtrs.ExpArrZFi[Two_iz_p_1];

				bRe = DataPtrs.FFTArrZFiEx[Two_ix]; bIm = DataPtrs.FFTArrZFiEx[Two_ix_p_1];
				ExRe -= (float)(DataPtrs.dzFi*(ExpZFiRe*bRe - ExpZFiIm*bIm));
				ExIm -= (float)(DataPtrs.dzFi*(ExpZFiRe*bIm + ExpZFiIm*bRe));

				bRe = DataPtrs.FFTArrZFiEz[Two_ix]; bIm = DataPtrs.FFTArrZFiEz[Two_ix_p_1];
				EzRe -= (float)(DataPtrs.dzFi*(ExpZFiRe*bRe - ExpZFiIm*bIm));
				EzIm -= (float)(DataPtrs.dzFi*(ExpZFiRe*bIm + ExpZFiIm*bRe));
			}

			*tEx = ExRe; *(tEx+1) = ExIm;
			*tEz = EzRe; *(tEz+1) = EzIm;
			tEx += 2; tEz += 2;
		}
	}
}

//*************************************************************************

int srTSRWRadStructAccessData::ShiftWfrByInterpolVsXZ(double shiftX, double shiftZ)
{//Shift the wavefront E-field data by interpolation, in whatever representation (coord. of ang.), keeping same mesh
 //Note: it also modifies xc, zc (requires for best treatment of quadratic phase term) !

	//long nTot = (ne << 1)*nx*nz;
	long long nTot = (ne << 1)*((long long)nx)*((long long)nz);
	float *pAuxBaseRadX = 0;
	float *pAuxBaseRadZ = 0;
	if(pBaseRadX != 0) 
	{
		pAuxBaseRadX = new float[nTot];
		float *tAuxBaseRadX = pAuxBaseRadX;
		//for(long i=0; i<nTot; i++) *(tAuxBaseRadX++) = 0;
		for(long long i=0; i<nTot; i++) *(tAuxBaseRadX++) = 0;
	}
	if(pBaseRadZ != 0)
	{
		pAuxBaseRadZ = new float[nTot];
		float *tAuxBaseRadZ = pAuxBaseRadZ;
		//for(long i=0; i<nTot; i++) *(tAuxBaseRadZ++) = 0;
		for(long long i=0; i<nTot; i++) *(tAuxBaseRadZ++) = 0;
	}

	char PolComp = -1;
	bool WaveFrontTermWasTreated = false;
	if(QuadPhaseTermCanBeTreated())
	{
		if(pBaseRadX != 0)
		{
			if(pBaseRadZ != 0) PolComp = 0;
			else PolComp = 'x';
		}
		else if(pBaseRadZ != 0) PolComp = 'z';
		
		TreatQuadPhaseTermTerm('r', PolComp);
		WaveFrontTermWasTreated = true;
	}

	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long nx_mi_1 = nx - 1;
	long nz_mi_1 = nz - 1;
	double arF[5];

	for(long ie=0; ie<ne; ie++)
	{
		long Two_ie = ie << 1;
		double z = zStart - shiftZ;

		for(long iz=0; iz<nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			float *pEX_NewStartForX = pAuxBaseRadX + izPerZ;
			float *pEZ_NewStartForX = pAuxBaseRadZ + izPerZ;

			double d_izOld = (z - zStart)/zStep;
			if((d_izOld < 0) || (d_izOld > nz_mi_1)) 
			{
				z += zStep; continue;
			}

			long izOld = (long)d_izOld;
			if((d_izOld - izOld) >= 0.5) izOld++;
			if(izOld < 0) izOld = 0;
			else if(izOld > nz_mi_1) izOld = nz_mi_1;

			long izOld_mi_1 = izOld - 1;
			if(izOld_mi_1 < 0) izOld_mi_1 = 0;
			long izOld_pl_1 = izOld + 1;
			if(izOld_pl_1 > nz_mi_1) izOld_pl_1 = nz_mi_1;

			double rz = z - (zStart + izOld*zStep);
			double zt = (zStep > 0)? rz/zStep : 0.;

			//long izOld_PerZ = izOld*PerZ;
			//long izOld_mi_1_PerZ = izOld_mi_1*PerZ;
			//long izOld_pl_1_PerZ = izOld_pl_1*PerZ;
			long long izOld_PerZ = izOld*PerZ;
			long long izOld_mi_1_PerZ = izOld_mi_1*PerZ;
			long long izOld_pl_1_PerZ = izOld_pl_1*PerZ;

			double x = xStart - shiftX;

			for(long ix=0; ix<nx; ix++)
			{
				//long ixPerX_p_Two_ie = ix*PerX + Two_ie; //offset for the new data
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie; //offset for the new data
				float *pEX_New = pEX_NewStartForX + ixPerX_p_Two_ie;
				float *pEZ_New = pEZ_NewStartForX + ixPerX_p_Two_ie;

				double d_ixOld = (x - xStart)/xStep;
				if((d_ixOld < 0) || (d_ixOld > nx_mi_1)) 
				{
					x += xStep; continue;
				}

				long ixOld = (long)d_ixOld;
				if((d_ixOld - ixOld) >= 0.5) ixOld++;
				if(ixOld < 0) ixOld = 0;
				else if(ixOld > nx_mi_1) ixOld = nx_mi_1;

				long ixOld_mi_1 = ixOld - 1;
				if(ixOld_mi_1 < 0) ixOld_mi_1 = 0;
				long ixOld_pl_1 = ixOld + 1;
				if(ixOld_pl_1 > nx_mi_1) ixOld_pl_1 = nx_mi_1;

				double rx = x - (xStart + ixOld*xStep);
				double xt = (xStep > 0)? rx/xStep : 0.;

				//long ixOld_PerX = ixOld*PerX;
				//long ixOld_mi_1_PerX = ixOld_mi_1*PerX;
				//long ixOld_pl_1_PerX = ixOld_pl_1*PerX;
				long long ixOld_PerX = ixOld*PerX;
				long long ixOld_mi_1_PerX = ixOld_mi_1*PerX;
				long long ixOld_pl_1_PerX = ixOld_pl_1*PerX;
				
				//long ofstOld_0m1 = ixOld_PerX + izOld_mi_1_PerZ + Two_ie; //offset for the new data
				//long ofstOld_m10 = ixOld_mi_1_PerX + izOld_PerZ + Two_ie;
				//long ofstOld_00 = ixOld_PerX + izOld_PerZ + Two_ie;
				//long ofstOld_10 = ixOld_pl_1_PerX + izOld_PerZ + Two_ie;
				//long ofstOld_01 = ixOld_PerX + izOld_pl_1_PerZ + Two_ie;
				long long ofstOld_0m1 = ixOld_PerX + izOld_mi_1_PerZ + Two_ie; //offset for the new data
				long long ofstOld_m10 = ixOld_mi_1_PerX + izOld_PerZ + Two_ie;
				long long ofstOld_00 = ixOld_PerX + izOld_PerZ + Two_ie;
				long long ofstOld_10 = ixOld_pl_1_PerX + izOld_PerZ + Two_ie;
				long long ofstOld_01 = ixOld_PerX + izOld_pl_1_PerZ + Two_ie;

				//long ofstOld_0m1_p1 = ofstOld_0m1 + 1; //offset for the new data
				//long ofstOld_m10_p1 = ofstOld_m10 + 1;
				//long ofstOld_00_p1 = ofstOld_00 + 1;
				//long ofstOld_10_p1 = ofstOld_10 + 1;
				//long ofstOld_01_p1 = ofstOld_01 + 1;
				long long ofstOld_0m1_p1 = ofstOld_0m1 + 1; //offset for the new data
				long long ofstOld_m10_p1 = ofstOld_m10 + 1;
				long long ofstOld_00_p1 = ofstOld_00 + 1;
				long long ofstOld_10_p1 = ofstOld_10 + 1;
				long long ofstOld_01_p1 = ofstOld_01 + 1;

				if(pBaseRadX != 0)
				{
					double *t_arF = arF;
					*(t_arF++) = pBaseRadX[ofstOld_0m1]; //double f0m1 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_m10]; //double fm10 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_00]; //double f00 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_10]; //double f10 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_01]; //double f01 = *arF;
					*pEX_New = (float)CGenMathInterp::Interp2dBiQuad5Rec(xt, zt, arF);
					//OCTEST
					//*pEX_New = (float)pBaseRadX[ofstOld_00];

					t_arF = arF;
					*(t_arF++) = pBaseRadX[ofstOld_0m1_p1]; //double f0m1 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_m10_p1]; //double fm10 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_00_p1]; //double f00 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_10_p1]; //double f10 = *(arF++);
					*(t_arF++) = pBaseRadX[ofstOld_01_p1]; //double f01 = *arF;
					*(pEX_New + 1) = (float)CGenMathInterp::Interp2dBiQuad5Rec(xt, zt, arF);
					//OCTEST
					//*(pEX_New + 1) = (float)pBaseRadX[ofstOld_00_p1];
				}
				if(pBaseRadZ != 0)
				{
					double *t_arF = arF;
					*(t_arF++) = pBaseRadZ[ofstOld_0m1]; //double f0m1 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_m10]; //double fm10 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_00]; //double f00 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_10]; //double f10 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_01]; //double f01 = *arF;
					*pEZ_New = (float)CGenMathInterp::Interp2dBiQuad5Rec(xt, zt, arF);
					//OCTEST
					//*pEZ_New = (float)pBaseRadZ[ofstOld_00];

					t_arF = arF;
					*(t_arF++) = pBaseRadZ[ofstOld_0m1_p1]; //double f0m1 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_m10_p1]; //double fm10 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_00_p1]; //double f00 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_10_p1]; //double f10 = *(arF++);
					*(t_arF++) = pBaseRadZ[ofstOld_01_p1]; //double f01 = *arF;
					*(pEZ_New + 1) = (float)CGenMathInterp::Interp2dBiQuad5Rec(xt, zt, arF);
					//OCTEST
					//*(pEZ_New + 1) = (float)pBaseRadZ[ofstOld_00_p1];
				}
				x += xStep;
			}
			z += zStep;
		}
	}

	if(pBaseRadX != 0) 
	{
		float *tAuxRadX = pAuxBaseRadX, *tRadX = pBaseRadX;
		//for(long i=0; i<nTot; i++) *(tRadX++) = *(tAuxRadX++);
		for(long long i=0; i<nTot; i++) *(tRadX++) = *(tAuxRadX++);
	}
	if(pBaseRadZ != 0) 
	{
		float *tAuxRadZ = pAuxBaseRadZ, *tRadZ = pBaseRadZ;
		//for(long i=0; i<nTot; i++) *(tRadZ++) = *(tAuxRadZ++);
		for(long long i=0; i<nTot; i++) *(tRadZ++) = *(tAuxRadZ++);
	}

	//OC180813: don't correct it here; will be corrected in sep. function
	xc += shiftX; 
	zc += shiftZ;

	if(WaveFrontTermWasTreated) TreatQuadPhaseTermTerm('a', PolComp);

	//OC180813: don't correct it here; will be corrected in sep. function
	xc -= shiftX; 
	zc -= shiftZ;

	if(pAuxBaseRadX != 0) delete[] pAuxBaseRadX;
	if(pAuxBaseRadZ != 0) delete[] pAuxBaseRadZ;
	return 0;
}

//*************************************************************************

void srTSRWRadStructAccessData::FlipFieldData(bool flipOverX, bool flipOverZ)
{
	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;

	long halfNz = nz >> 1, nz_mi_1 = nz - 1;
	long halfNx = nx >> 1, nx_mi_1 = nx - 1;

	bool treatEx = (pBaseRadX != 0);
	bool treatEz = (pBaseRadZ != 0);

	if(flipOverZ)
	{
		if(flipOverX)
		{
			for(long iz=0; iz<halfNz; iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				for(long ix=0; ix<halfNx; ix++)
				{
					//long offset = izPerZ + ix*PerX;
					long long offset = izPerZ + ix*PerX;
					float* pOrigDataEx = pBaseRadX + offset;
					float* pOrigDataEz = pBaseRadZ + offset;

					//long offsetSym = izPerZ + (nx_mi_1 - ix)*PerX;
					long long offsetSym = izPerZ + (nx_mi_1 - ix)*PerX;
					float* pSymDataEx = pBaseRadX + offsetSym;
					float* pSymDataEz = pBaseRadZ + offsetSym;
					SwapDataInEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, treatEx, treatEz);
				}
			}
		}
		for(long iz=0; iz<halfNz; iz++)
		{
			//long izPerZ = iz*PerZ, BufZ = (nz_mi_1 - iz)*PerZ;
			long long izPerZ = iz*PerZ, BufZ = (nz_mi_1 - iz)*PerZ;
			for(long ix=0; ix<nx; ix++)
			{			
				//long ixPerX = ix*PerX;
				//long offset = izPerZ + ixPerX;
				long long ixPerX = ix*PerX;
				long long offset = izPerZ + ixPerX;
				float* pOrigDataEx = pBaseRadX + offset;
				float* pOrigDataEz = pBaseRadZ + offset;

				//long offsetSym = BufZ + ixPerX;
				long long offsetSym = BufZ + ixPerX;
				float* pSymDataEx = pBaseRadX + offsetSym;
				float* pSymDataEz = pBaseRadZ + offsetSym;
				SwapDataInEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, treatEx, treatEz);
			}
		}
	}
	else if(flipOverX)
	{
		for(long iz=0; iz<nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			for(long ix=0; ix<halfNx; ix++)
			{
				//long offset = izPerZ + ix*PerX;
				long long offset = izPerZ + ix*PerX;
				float* pOrigDataEx = pBaseRadX + offset;
				float* pOrigDataEz = pBaseRadZ + offset;

				//long offsetSym = izPerZ + (nx_mi_1 - ix)*PerX;
				long long offsetSym = izPerZ + (nx_mi_1 - ix)*PerX;
				float* pSymDataEx = pBaseRadX + offsetSym;
				float* pSymDataEz = pBaseRadZ + offsetSym;
				SwapDataInEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, treatEx, treatEz);
			}
		}
	}
}

//*************************************************************************

void srTSRWRadStructAccessData::TransposeFieldData()
{//To improve (make more memory-economic version)

	double xStartOld = xStart, xStepOld = xStep;
	long nxOld = nx;
	xStart = zStart; xStep = zStep; nx = nz;
	zStart = xStartOld; zStep = xStepOld; nz = nxOld;

	//long PerX = ne << 1;
	//long PerZ = PerX*nx;
	long long PerX = ne << 1;
	long long PerZ = PerX*nx;

	//long nTot = (ne << 1)*nx*nz;
	long long nTot = (ne << 1)*((long long)nx)*((long long)nz);
	float *arE = new float[nTot];

	if(pBaseRadX != 0)
	{
		float *t_arE = arE;
		float *t_pRad = pBaseRadX;
		//for(long i=0; i<nTot; i++) *(t_arE++) = *(t_pRad++);
		for(long long i=0; i<nTot; i++) *(t_arE++) = *(t_pRad++);

		t_arE = arE;
		for(long ix = 0; ix < nx; ix++)
		{
			for(long iz = 0; iz < nz; iz++)
			{
				//long ofst0 = iz*PerZ + ix*PerX;
				long long ofst0 = iz*PerZ + ix*PerX;
				for(long ie = 0; ie < ne; ie++)
				{
					//long ofst = ofst0 + (ie << 1);
					long long ofst = ofst0 + (ie << 1);
					*(pBaseRadX + ofst) = *(t_arE++);
					*(pBaseRadX + ofst + 1) = *(t_arE++);
				}
			}
		}
	}
	if(pBaseRadZ != 0)
	{
		float *t_arE = arE;
		float *t_pRad = pBaseRadZ;
		//for(long i=0; i<nTot; i++) *(t_arE++) = *(t_pRad++);
		for(long long i=0; i<nTot; i++) *(t_arE++) = *(t_pRad++);

		t_arE = arE;
		for(long ix = 0; ix < nx; ix++)
		{
			for(long iz = 0; iz < nz; iz++)
			{
				//long ofst0 = iz*PerZ + ix*PerX;
				long long ofst0 = iz*PerZ + ix*PerX;
				for(long ie = 0; ie < ne; ie++)
				{
					//long ofst = ofst0 + (ie << 1);
					long long ofst = ofst0 + (ie << 1);
					*(pBaseRadZ + ofst) = *(t_arE++);
					*(pBaseRadZ + ofst + 1) = *(t_arE++);
				}
			}
		}
	}
	if(arE != 0) delete[] arE;
}

//*************************************************************************

int srTSRWRadStructAccessData::SetRepresCA(char CoordOrAng)
{//Copied from srTGenOptElem::SetRadRepres(srTSRWRadStructAccessData*, char);
// 'c' or 'C' or 0- to coord.; 'a' or 'A' or 1- to ang.
// 0- to coord.; 1- to ang.
	int result=0;
	if((CoordOrAng == 'c') || (CoordOrAng == 'C')) CoordOrAng = 0;
	else CoordOrAng = 1;

	//char WfrEdgeCorrShouldBeTreated = pRadAccessData->WfrEdgeCorrShouldBeDone; // Turn on/off here

	if(Pres == CoordOrAng) return 0;
	char DirFFT = (CoordOrAng == 0)? -1 : 1;

	CGenMathFFT2DInfo FFT2DInfo;

	FFT2DInfo.xStep = xStep;
	FFT2DInfo.yStep = zStep;
	FFT2DInfo.xStart = xStart;
	FFT2DInfo.yStart = zStart;
	FFT2DInfo.Nx = nx;
	FFT2DInfo.Ny = nz;
	FFT2DInfo.Dir = DirFFT;
	FFT2DInfo.UseGivenStartTrValues = 0;

//New
	if((AuxLong4 == 7777777) || ((CoordOrAng == 0) && UseStartTrToShiftAtChangingRepresToCoord))
	{
		FFT2DInfo.UseGivenStartTrValues = 1;
		FFT2DInfo.xStartTr = xStartTr;
		FFT2DInfo.yStartTr = zStartTr;
	}
//End New

	CGenMathFFT2D FFT2D;
	const double constPhotEnWavelenConv = 1.239842e-06;
	double avgWaveLength_m = 0.;

	if(ne == 1)
	{
		srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
		if(WfrEdgeCorrShouldBeDone)
		{
			if(CoordOrAng == 1)
			{
				if(result = SetupWfrEdgeCorrData(pBaseRadX, pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;
			}
		}

		if(ElecFldAngUnit == 1) //OC20112017
		{
			double avgPhotEnLoc = (PresT == 0)? eStart : avgPhotEn;
			avgWaveLength_m = constPhotEnWavelenConv/avgPhotEnLoc;
			if(CoordOrAng == 0)
			{
				FFT2DInfo.ExtraMult = avgWaveLength_m;
				double multArg = 1./avgWaveLength_m;
				FFT2DInfo.xStart *= multArg; FFT2DInfo.xStep *= multArg;
				FFT2DInfo.yStart *= multArg; FFT2DInfo.yStep *= multArg;
			}
			else FFT2DInfo.ExtraMult = 1./avgWaveLength_m;
		}

		FFT2DInfo.pData = pBaseRadX;
		if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
		FFT2DInfo.pData = pBaseRadZ;
		if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

		if(WfrEdgeCorrShouldBeDone)
		{
			if(CoordOrAng == 1)
			{
				if(DataPtrsForWfrEdgeCorr.WasSetup)
				{
					MakeWfrEdgeCorrection(pBaseRadX, pBaseRadZ, DataPtrsForWfrEdgeCorr);
					DataPtrsForWfrEdgeCorr.DisposeData();
				}
			}
		}
	}
	else
	{
		if(ElecFldAngUnit == 1) //OC20112017
		{
			double avgPhotEnLoc = (PresT == 0)? (eStart + 0.5*ne*eStep) : avgPhotEn;
			avgWaveLength_m = constPhotEnWavelenConv/avgPhotEnLoc;
			if(CoordOrAng == 0)
			{
				FFT2DInfo.ExtraMult = avgWaveLength_m; //NOTE: this may be incorrect for large bandwidth !?
				double multArg = 1./avgWaveLength_m;
				FFT2DInfo.xStart *= multArg; FFT2DInfo.xStep *= multArg;
				FFT2DInfo.yStart *= multArg; FFT2DInfo.yStep *= multArg;
			}
			else FFT2DInfo.ExtraMult = 1./avgWaveLength_m; //NOTE: this may be incorrect for large bandwidth !?
		}

		//long TwoNxNz = (nx*nz) << 1;
		long long TwoNxNz = (((long long)nx)*((long long)nz)) << 1;
		float* AuxEx = new float[TwoNxNz];
		if(AuxEx == 0) return MEMORY_ALLOCATION_FAILURE;
		float* AuxEz = new float[TwoNxNz];
		if(AuxEz == 0) return MEMORY_ALLOCATION_FAILURE;

		for(long ie=0; ie<ne; ie++)
		{
			if(result = ExtractSliceConstEorT(ie, AuxEx, AuxEz)) return result;

			srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
			if(WfrEdgeCorrShouldBeDone)
			{
				if(CoordOrAng == 1)
				{
					if(result = SetupWfrEdgeCorrData(AuxEx, AuxEz, DataPtrsForWfrEdgeCorr)) return result;
				}
			}

			FFT2DInfo.pData = AuxEx;
			if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
			FFT2DInfo.pData = AuxEz;
			if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

			if(WfrEdgeCorrShouldBeDone)
			{
				if(CoordOrAng == 1)
				{
					if(DataPtrsForWfrEdgeCorr.WasSetup)
					{
						MakeWfrEdgeCorrection(AuxEx, AuxEz, DataPtrsForWfrEdgeCorr);
						DataPtrsForWfrEdgeCorr.DisposeData();
					}
				}
			}

			if(result = SetupSliceConstEorT(ie, AuxEx, AuxEz)) return result;
		}

		if(AuxEx != 0) delete[] AuxEx;
		if(AuxEz != 0) delete[] AuxEz;
	}

	xStep = FFT2DInfo.xStepTr;
	zStep = FFT2DInfo.yStepTr;
	xStart = FFT2DInfo.xStartTr;
	zStart = FFT2DInfo.yStartTr;

	if((ElecFldAngUnit == 1) && (CoordOrAng == 1)) //OC20112017
	{
		double multArg = avgWaveLength_m;
		xStart *= multArg; xStep *= multArg;
		zStart *= multArg; zStep *= multArg;
	}

	Pres = CoordOrAng;

	SetNonZeroWavefrontLimitsToFullRange();
	return result;
}

//*************************************************************************

int srTSRWRadStructAccessData::SetRepresFT(char FreqOrTime) 
{//set Frequency or Time representation
// 'f' or 'F' or 0- to freq.; 't' or 'T' or 1- to time
//Conversions are done assuming intensity units to be:
//	- in Time domain: [W/mm^2]
//	- in Frequency domain: [Photons/0.1%bw/mm^2]

	int result=0;
	if((FreqOrTime == 'f') || (FreqOrTime == 'F')) FreqOrTime = 0; //to Frequency
	else FreqOrTime = 1; //to Time

	if(FreqOrTime == PresT) return 0;
	if(ne <= 1) return 0; //nothing to FT

	//char DirFFT = (FreqOrTime == 0)? -1 : 1;

	const double multConvHz2eV = 4.135667175e-15;
	//const double multConv_eV2PhperBW = 6.24146e+15;
	double multConv_eV2PhperBW = 6.24146e+15; 
	if(ElecFldUnit == 2) multConv_eV2PhperBW = 1; //OC170813

	//default- to time:
	double normArgConv = multConvHz2eV;
	//double multPreInt = sqrt(multConvHz2eV);
	double multPreInt = sqrt(multConvHz2eV/multConv_eV2PhperBW);

	//char DirFFT = -1; //1; //OCTEST111211
	char DirFFT = 1; //OC041215

	double shiftArgBefore = -avgPhotEn;
	//double shiftArgAfter = 0;
	double shiftArgAfter = avgT; //OC101115 (Need to define it prior to this!)

	if(FreqOrTime == 0) //to Frequency
	{
		//multConv = multConvHz2eV;
		//multPreInt = sqrt(multConvHz2eV);
		//multPreInt *= sqrt(multConv_eV2PhperBW);

		multPreInt = sqrt(multConvHz2eV*multConv_eV2PhperBW);
		//DirFFT = 1; //-1; //OCTEST111211
		DirFFT = -1; //OC041215

		avgT = eStart + 0.5*eStep*(ne - 1); //OC101115 //???

		//shiftArgBefore = 0;
		shiftArgBefore = -avgT; //OC101115
		shiftArgAfter = avgPhotEn;
	}

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = eStep/normArgConv; //[s] or [eV]
	FFT1DInfo.xStart = (eStart + shiftArgBefore)/normArgConv; //[s] or [eV]
	FFT1DInfo.Nx = ne;
	FFT1DInfo.HowMany = nx*nz; //May result in overflow?
	FFT1DInfo.Dir = DirFFT;
	FFT1DInfo.UseGivenStartTrValue = 0;
	FFT1DInfo.MultExtra = multPreInt; //1./multConv; //???
	FFT1DInfo.ApplyAutoShiftAfter = false;

	CGenMathFFT1D FFT1D;
	if(pBaseRadX != 0)
	{
		FFT1DInfo.pInData = pBaseRadX;
		FFT1DInfo.pOutData = pBaseRadX;
		if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;
	}
	if(pBaseRadZ != 0)
	{
		FFT1DInfo.pInData = pBaseRadZ;
		FFT1DInfo.pOutData = pBaseRadZ;
		if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;
	}

	eStep = FFT1DInfo.xStepTr;
	eStart = FFT1DInfo.xStartTr + shiftArgAfter;

	if(FreqOrTime == 1) //to Time
	{//OC101015 //???
		double tCenAux = eStart + 0.5*eStep*(ne - 1);
		double tShiftAux = avgT - tCenAux;
		eStart += tShiftAux;
	}

	PresT = FreqOrTime;
	return result;
}

//*************************************************************************

int srTSRWRadStructAccessData::ComputeRadMoments()
{
	srTGenOptElem GenOptElem;
	return GenOptElem.ComputeRadMoments(this);
}

//*************************************************************************

bool srTSRWRadStructAccessData::CheckIfQuadTermTreatIsBenefit(char cutX_or_Z, char fldX_or_Z)
{
	if((pBaseRadX == 0) && (pBaseRadZ == 0)) return false;

	if((Pres != 0) && (PresT != 0)) return true; //this test is currently imlemented only for the Coordinate-Frequency domain (other cases to consider)

	//long ieCen = 0;
	long long ieCen = 0;
	if(ne > 1) ieCen = ne >> 1;
	double eCen = eStart + eStep*ieCen;
	double halfWaveNum = 0.5*(5.06773065E+06)*eCen; //Pi/lambda_m

	//long ofst0, per;
	long long ofst0, per;
	double argStep, argStart, argN;
	double argCen, argR;

	if((cutX_or_Z == 'x') || (cutX_or_Z == 'X'))
	{
		per = ne << 1;

		long nOtherMi1 = nz - 1;
		double dicOther = (zc - zStart)/zStep;
		long icOther = (long)dicOther;
		if((dicOther - icOther) >= 0.5) icOther++;
		if(icOther < 0) icOther = 0;
		else if(icOther > nOtherMi1) icOther = nOtherMi1;

		ofst0 = (ieCen << 1) + icOther*(per*nx);
		argStep = xStep;
		argStart = xStart;
		argN = nx;
		argCen = xc;
		argR = RobsX;
	}
	else
	{
		per = (ne << 1)*nx;

		long nOtherMi1 = nx - 1;
		double dicOther = (xc - xStart)/xStep;
		long icOther = (long)dicOther;
		if((dicOther - icOther) >= 0.5) icOther++;
		if(icOther < 0) icOther = 0;
		else if(icOther > nOtherMi1) icOther = nOtherMi1;

		ofst0 = icOther*(ne << 1) + (ieCen << 1);
		argStep = zStep;
		argStart = zStart;
		argN = nz;
		argCen = zc;
		argR = RobsZ;
	}
	double invArgStep = 1./argStep;
	double coefPh = halfWaveNum/argR;

	bool fldX_ShouldBeTreated = (fldX_or_Z != 'z') && (fldX_or_Z != 'Z') && (fldX_or_Z != 'y') && (fldX_or_Z != 'Y') && (pBaseRadX != 0);
	bool fldZ_ShouldBeTreated = (fldX_or_Z != 'x') && (fldX_or_Z != 'X') && (pBaseRadZ != 0);

	float *pFldX = pBaseRadX + ofst0, *pFldZ = pBaseRadZ + ofst0;
	float *tFldX = pFldX, *tFldZ = pFldZ;
	double reEx=0, imEx=0, reEz=0, imEz=0;
	double maxIntX=0., maxIntZ=0.;
	for(long i=0; i<argN; i++)
	{
		if(fldX_ShouldBeTreated)
		{
			reEx = *tFldX; imEx = *(tFldX + 1); tFldX += per;
			double curIntX = reEx*reEx + imEx*imEx;
			if(maxIntX < curIntX) maxIntX = curIntX;
		}
		if(fldZ_ShouldBeTreated)
		{
			reEz = *tFldZ; imEz = *(tFldZ + 1); tFldZ += per;
			double curIntZ = reEz*reEz + imEz*imEz;
			if(maxIntZ < curIntZ) maxIntZ = curIntZ;
		}
	}

	float *tFld = pFldX;
	double maxInt = maxIntX;
	if(maxIntZ > maxIntX)
	{
		tFld = pFldZ; maxInt = maxIntZ;
	}

	double threshInt = 0.01*maxInt; //to tune!

	double re=0, im=0;
	double prevRe=0, prevIm=0;
	double prevReAfter=0, prevImAfter=0;

	double difArg = argStart - argCen;
	double phShift = coefPh*difArg*difArg;
	double cosPh = cos(phShift), sinPh = sin(phShift);
	prevRe = *tFld; prevIm = *(tFld + 1);
	prevReAfter = prevRe*cosPh + prevIm*sinPh;
	prevImAfter = prevIm*cosPh - prevRe*sinPh;

	double reAfter=0, imAfter=0;
	double derRe=0, derIm=0;
	double derReAfter=0, derImAfter=0;
	double prevDerRe=0, prevDerIm=0;
	double prevDerReAfter=0, prevDerImAfter=0;
	int numDerReSignChange=0, numDerImSignChange=0;
	int numDerReSignChangeAfter=0, numDerImSignChangeAfter=0;
	
	const double twoPi = 6.2831853;
	double phShiftPrev;

	for(long i=0; i<argN; i++)
	{
		phShift = coefPh*difArg*difArg;

		bool phShiftIsSmall = false;
		if(i > 0) //OC05012017
		{
			if(::fabs(phShift - phShiftPrev) < twoPi) phShiftIsSmall = true;
		}
		phShiftPrev = phShift;

		difArg += argStep;
		cosPh = cos(phShift); sinPh = sin(phShift);

		re = *tFld; im = *(tFld + 1); tFld += per;
		derRe = (re - prevRe)*invArgStep;
		derIm = (im - prevIm)*invArgStep;

		reAfter = re*cosPh + im*sinPh;
		imAfter = im*cosPh - re*sinPh;
		derReAfter = (reAfter - prevReAfter)*invArgStep;
		derImAfter = (imAfter - prevImAfter)*invArgStep;

		double curInt = re*re + im*im;
		//if(curInt > threshInt)
		if((curInt > threshInt) && phShiftIsSmall) //OC05012017 //take into account only "resolved" oscillations
		{
			if(derRe*prevDerRe < 0.) numDerReSignChange++;
			if(derIm*prevDerIm < 0.) numDerImSignChange++;
			if(derReAfter*prevDerReAfter < 0.) numDerReSignChangeAfter++;
			if(derImAfter*prevDerImAfter < 0.) numDerImSignChangeAfter++;
		}
		prevDerRe = derRe; prevDerIm = derIm;
		prevRe = re; prevIm = im;
		prevDerReAfter = derReAfter; prevDerImAfter = derImAfter;
		prevReAfter = reAfter; prevImAfter = imAfter;
	}

	int numDerE_SignChange = numDerReSignChange, numDerE_SignChangeAfter = numDerReSignChangeAfter;
	if(numDerImSignChange > numDerReSignChange)
	{
		numDerE_SignChange = numDerImSignChange;
		numDerE_SignChangeAfter = numDerImSignChangeAfter;
	}
	
	return (numDerE_SignChangeAfter <= numDerE_SignChange);
}

//*************************************************************************

void srTSRWRadStructAccessData::GetIntMesh(char dep, SRWLRadMesh& mesh) //OC23082018
{//This assumes center values for the intensity distribution are defined in mesh.eStart, mesh.xStart, mesh.yStart at input
	mesh.ne = mesh.nx = mesh.ny = 1;
	if(dep == 0) 
	{
		mesh.ne = ne;
		mesh.eStart = eStart; 
		mesh.eFin = eStart + eStep*(ne - 1);
		//Keep mesh.xStart, mesh.yStart as they define "central" values of the intensity distribution
	}
	else if(dep == 1) 
	{
		mesh.nx = nx;
		mesh.xStart = xStart; 
		mesh.xFin = xStart + xStep*(nx - 1);
		//Keep mesh.eStart, mesh.yStart as they define "central" values of the intensity distribution
	}
	else if(dep == 2) 
	{
		mesh.ny = nz;
		mesh.yStart = zStart; 
		mesh.yFin = zStart + zStep*(nz - 1);
		//Keep mesh.eStart, mesh.xStart as they define "central" values of the intensity distribution
	}
	else if(dep == 3) 
	{
		mesh.nx = nx;
		mesh.xStart = xStart; 
		mesh.xFin = xStart + xStep*(nx - 1);
		mesh.ny = nz;
		mesh.yStart = zStart; 
		mesh.yFin = zStart + zStep*(nz - 1);
		//Keep mesh.eStart as it defines "central" value of the intensity distribution
	}
	else if(dep == 4) 
	{
		mesh.ne = ne;
		mesh.eStart = eStart; 
		mesh.eFin = eStart + eStep*(ne - 1);
		mesh.nx = nx;
		mesh.xStart = xStart; 
		mesh.xFin = xStart + xStep*(nx - 1);
		//Keep mesh.yStart as it defines "central" value of the intensity distribution
	}
	else if(dep == 5) 
	{
		mesh.ne = ne;
		mesh.eStart = eStart; 
		mesh.eFin = eStart + eStep*(ne - 1);
		mesh.ny = nz;
		mesh.yStart = zStart; 
		mesh.yFin = zStart + zStep*(nz - 1);
		//Keep mesh.xStart as it defines "central" value of the intensity distribution
	}
	else if(dep == 6) 
	{
		mesh.ne = ne;
		mesh.eStart = eStart; 
		mesh.eFin = eStart + eStep*(ne - 1);
		mesh.nx = nx;
		mesh.xStart = xStart; 
		mesh.xFin = xStart + xStep*(nx - 1);
		mesh.ny = nz;
		mesh.yStart = zStart; 
		mesh.yFin = zStart + zStep*(nz - 1);
	}
}

//*************************************************************************
/**
void srTSRWRadStructAccessData::EstimWfrRadCen(double& resR, double& resCen, char cutX_or_Z, char fldX_or_Z, double relArgRange, double relArgCenOther)
{
	resR = 0; resCen = 0;

	if((pBaseRadX == 0) && (pBaseRadZ == 0)) return;
	if((relArgRange <= 0.) || (relArgRange > 1.)) return;

	long ieCen = 0;
	if(ne > 1) ieCen = ne >> 1;
	double eCen = eStart + eStep*ieCen;
	double waveNum = (5.06773065E+06)*eCen; //2*Pi/lambda_m

	const long minNp = 5; //to tune
	long ofst0, per, iStart, iEnd;
	double argStep, argStart;

	if((cutX_or_Z == 'x') || (cutX_or_Z == 'X'))
	{
		per = ne << 1;

		long nOtherMi1 = nz - 1;
		double dicOther = nOtherMi1*relArgCenOther;
		long icOther = (long)dicOther;
		if((dicOther - icOther) >= 0.5) icOther++;
		if(icOther < 0) icOther = 0;
		else if(icOther > nOtherMi1) icOther = nOtherMi1;

		double diStart = 0.5*(1. - relArgRange)*(nx - 1);
		iStart = (long)diStart;
		if((diStart - iStart) >= 0.5) iStart++;
		if(iStart < 0) iStart = 0;
		long actNp = nx - 2*iStart;
		if(actNp < minNp) iStart = (nx - actNp) >> 1;
		if(iStart < 0) iStart = 0;

		iEnd = nx - 1 - iStart;
		ofst0 = (ieCen << 1) + iStart*per + icOther*(per*nx);
		argStep = xStep;
		argStart = xStart;
	}
	else
	{
		per = (ne << 1)*nx;

		long nOtherMi1 = nx - 1;
		double dicOther = nOtherMi1*relArgCenOther;
		long icOther = (long)dicOther;
		if((dicOther - icOther) >= 0.5) icOther++;
		if(icOther < 0) icOther = 0;
		else if(icOther > nOtherMi1) icOther = nOtherMi1;

		double diStart = 0.5*(1. - relArgRange)*(nz - 1);
		iStart = (long)diStart;
		if((diStart - iStart) >= 0.5) iStart++;
		if(iStart < 0) iStart = 0;
		long actNp = nz - 2*iStart;
		if(actNp < minNp) iStart = (nz - actNp) >> 1;
		if(iStart < 0) iStart = 0;

		iEnd = nz - 1 - iStart;
		ofst0 = iStart*per + icOther*(ne << 1) + (ieCen << 1);
		argStep = zStep;
		argStart = zStart;
	}

	float *pFld = 0;
	if((fldX_or_Z == 'x') || (fldX_or_Z == 'X')) pFld = pBaseRadX;
	else if((fldX_or_Z == 'z') || (fldX_or_Z == 'Z') || (fldX_or_Z == 'y') || (fldX_or_Z == 'Y')) pFld = pBaseRadZ;
	else if(fldX_or_Z == 0)
	{
		if(pBaseRadX == 0) pFld = pBaseRadZ;
		else if(pBaseRadZ == 0) pFld = pBaseRadX;
		else
		{//Find max. field component
			double maxIntX = 0, maxIntZ = 0;
			float *tBaseRadX = pBaseRadX + ofst0, *tBaseRadZ = pBaseRadZ + ofst0;
			for(long i=iStart; i<=iEnd; i++)
			{
				double reEx = *tBaseRadX, imEx = *(tBaseRadX + 1);
				double curIntX = reEx*reEx + imEx*imEx;
				if(maxIntX < curIntX) maxIntX = curIntX;
				tBaseRadX += per;

				double reEz = *tBaseRadZ, imEz = *(tBaseRadZ + 1);
				double curIntZ = reEz*reEz + imEz*imEz;
				if(maxIntZ < curIntZ) maxIntZ = curIntZ;
				tBaseRadZ += per;
			}
			if(maxIntX >= maxIntZ) pFld = pBaseRadX;
			else pFld = pBaseRadZ;
		}
	}

	double invArgStep = 1./argStep;
	double invTwoArgStep = 0.5*invArgStep;

	long np = (iEnd - iStart + 1) << 1;
	double *arA1 = new double[np];
	double *arA2 = new double[np];
	double *arB = new double[np];

	float *tFld = pFld + ofst0;
	double ReF = *tFld, ImF = *(tFld + 1);
	double dReFdx = (*(tFld + 2) - ReF)*invArgStep;
	double dImFdx = (*(tFld + 3) - ImF)*invArgStep;
	tFld += 2;
	double arg = argStart + argStep*iStart;
	double *t_arA1 = arA1, *t_arA2 = arA2, *t_arB = arB;

	double waveNum_ImF = waveNum*ImF, waveNum_ReF = waveNum*ReF;
	*(t_arA1++) = -waveNum_ImF*arg;
	*(t_arA1++) = waveNum_ReF*arg;
	*(t_arA2++) = -waveNum_ImF;
	*(t_arA2++) = waveNum_ReF;
	*(t_arB++) = dReFdx;
	*(t_arB++) = dImFdx;
	arg += argStep;

	for(long i=(iStart+1); i<iEnd; i++)
	{
		ReF = *tFld; ImF = *(tFld + 1);
		dReFdx = (*(tFld + 2) - *(tFld - 2))*invTwoArgStep;
		dImFdx = (*(tFld + 3) - *(tFld - 1))*invTwoArgStep;
		tFld += 2;

		waveNum_ImF = waveNum*ImF; waveNum_ReF = waveNum*ReF;
		*(t_arA1++) = -waveNum_ImF*arg;
		*(t_arA1++) = waveNum_ReF*arg;
		*(t_arA2++) = -waveNum_ImF;
		*(t_arA2++) = waveNum_ReF;
		*(t_arB++) = dReFdx;
		*(t_arB++) = dImFdx;
		arg += argStep;
	}

	ReF = *tFld; ImF = *(tFld + 1);
	dReFdx = (ReF - *(tFld - 2))*invArgStep;
	dImFdx = (ImF - *(tFld - 1))*invArgStep;

	waveNum_ImF = waveNum*ImF; waveNum_ReF = waveNum*ReF;
	*(t_arA1++) = -waveNum_ImF*arg;
	*t_arA1 = waveNum_ReF*arg;
	*(t_arA2++) = -waveNum_ImF;
	*t_arA2 = waveNum_ReF;
	*(t_arB++) = dReFdx;
	*t_arB = dImFdx;

	//Finding solution (linear fit): (AT*A)^(-1)*AT*B
	double AtB1=0., AtB2=0., AtA11=0., AtA12=0., AtA22=0.;
	t_arA1 = arA1; t_arA2 = arA2; t_arB = arB;
	double a1, a2, b;

	for(long i=0; i<np; i++)
	{
		a1 = *(t_arA1++); a2 = *(t_arA2++); b = *(t_arB++);
		
		AtB1 += a1*b;
		AtB2 += a2*b;
		AtA11 += a1*a1; AtA12 += a1*a2;
		AtA22 += a2*a2;
	}
	double AtA21 = AtA12;

	double detAtA = AtA11*AtA22 - AtA12*AtA22;
	if(detAtA == 0.) return;

	double invDetAtA = 1./detAtA;
	double invAtA11 = invDetAtA*AtA22, invAtA12 = -invDetAtA*AtA12;
	double invAtA21 = -invDetAtA*AtA21, invAtA22 = invDetAtA*AtA11;
	double u1 = invAtA11*AtB1 + invAtA12*AtB2;
	double u2 = invAtA21*AtB1 + invAtA22*AtB2;

	resR = 1./u1; 
	resCen = -u2*resR;

	delete[] arA1;
	delete[] arA2;
	delete[] arB;
}
**/
//*************************************************************************
