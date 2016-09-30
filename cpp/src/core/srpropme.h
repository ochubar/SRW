/************************************************************************//**
 * File: srpropme.h
 * Description: Feasilibility tests for computation of emission and propagation of Partially-Coherent SR from Finite-Emittance Electron Beam (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRPROPME_H
#define __SRPROPME_H

#include "srstraux.h"
#include "sroptelm.h"
#include "gmrand.h"

//*************************************************************************

struct srTSigleElecVars {
	bool VarX, VarXp, VarZ, VarZp, VarE, VarPh;

	srTSigleElecVars() {}
	srTSigleElecVars(bool SameForAll) 
	{
		VarX = VarXp = VarZ = VarZp = VarE = VarPh = SameForAll;
	}
	srTSigleElecVars(srTEbmDat& EbmDat) 
	{
		VarX = VarXp = VarZ = VarZp = VarE = VarPh = false;

		if((EbmDat.Mxx > 0.) || (EbmDat.Mxxp > 0.)) VarX = true;
		if((EbmDat.Mxpxp > 0.) || (EbmDat.Mxxp > 0.)) VarXp = true;
		if((EbmDat.Mzz > 0.) || (EbmDat.Mzzp > 0.)) VarZ = true;
		if((EbmDat.Mzpzp > 0.) || (EbmDat.Mzzp > 0.)) VarZp = true;
		//if((EbmDat.SigmaRelE > 0.) || (EbmDat.Mee > 0.)) VarE = true;
		//add more if/when implemented
	}
	//srTSigleElecVars(double* pInData)
	//{
	//	VarX = VarXp = VarZ = VarZp = VarE = VarPh = false;
	//	if(pInData == 0) return;
	//	if(*pInData == 1) VarX = true;
	//	if(*(pInData+1) == 1) VarXp = true;
	//	if(*(pInData+2) == 1) VarZ = true;
	//	if(*(pInData+3) == 1) VarZp = true;
	//	if(*(pInData+4) == 1) VarE = true;
	//	if(*(pInData+5) == 1) VarPh = true;
	//}
};

//*************************************************************************

class srTPropagMultiE {

	static CGenMathRand gRandGen;

public:

	static int PropagateElecFieldStokes(srTEbmDat& EbmDat, srTSRWRadStructAccessData&, srTGenOptElemHndl, double* pPrecPar, srTStokesStructAccessData&);
	static int PropagateElecFieldStokesAuto(srTEbmDat& ThickEbmDat, srTSRWRadStructAccessData& InWfr, srTGenOptElemHndl OptHndl, double* pPrecPar, srTStokesStructAccessData& OutStokes);
	static int EmitPropagElecFieldStokes(srTTrjDat& TrjDat, srTGenOptElemHndl, double* pPrecPar, srTStokesStructAccessData& OutStokes);

	static int ReallocateStokesAccordingToWfr(srTSRWRadStructAccessData& LocWfr, srTStokesStructAccessData& OutStokes);
	
	//static int AddWfrToStokes(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long MacroPartCount, double& CurRelPrec);
	static int AddWfrToStokes(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long long MacroPartCount, double& CurRelPrec);
	//static int AddWfrToStokesWithInterpXZ(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long MacroPartCount);
	static int AddWfrToStokesWithInterpXZ(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long long MacroPartCount);

	static void SetupNextThinEbm(srTEbmDat& EbmDat, srTSigleElecVars&, srTEbmDat& OutThinEbmDat);
	static void SimulateWfrFromOffAxisEbm(srTEbmDat& OnAxisEbmDat, srTEbmDat& OffAxisEbmDat, srTSigleElecVars&, srTSRWRadStructAccessData& Wfr);

	static void CalcStokesFromE(float* tEx, float* tEz, float* Stokes)
	{
		float &EwX_Re = *tEx, &EwX_Im = *(tEx + 1);
		float &EwZ_Re = *tEz, &EwZ_Im = *(tEz + 1);
		double LinHor = EwX_Re*EwX_Re + EwX_Im*EwX_Im;
		double LinVer = EwZ_Re*EwZ_Re + EwZ_Im*EwZ_Im;
		Stokes[0] = (float)(LinHor + LinVer);
		Stokes[1] = (float)(LinHor - LinVer);
		Stokes[2] = (float)(-2.*(EwX_Re*EwZ_Re + EwX_Im*EwZ_Im));
		Stokes[3] = (float)(2.*(-EwX_Re*EwZ_Im + EwX_Im*EwZ_Re));
	}

	static int ReallocateStokesAccordingToWfr(srTWfrSmp& WfrSmp, srTStokesStructAccessData& OutStokes);
	//{
	//	int OldNe = OutStokes.ne, OldNx = OutStokes.nx, OldNz = OutStokes.nz;

	//	OutStokes.eStart = WfrSmp.LambStart;
	//	OutStokes.eStep = 0;
	//	if(WfrSmp.nLamb > 1) OutStokes.eStep = (WfrSmp.LambEnd - WfrSmp.LambStart)/WfrSmp.nLamb;
	//	OutStokes.ne = WfrSmp.nLamb;

	//	OutStokes.xStart = WfrSmp.xStart;
	//	OutStokes.xStep = 0;
	//	if(WfrSmp.nx > 1) OutStokes.xStep = (WfrSmp.xEnd - WfrSmp.xStart)/WfrSmp.nx;
	//	OutStokes.nx = WfrSmp.nx;

	//	OutStokes.zStart = WfrSmp.zStart;
	//	OutStokes.zStep = 0;
	//	if(WfrSmp.nz > 1) OutStokes.eStep = (WfrSmp.zEnd - WfrSmp.zStart)/WfrSmp.nz;
	//	OutStokes.nz = WfrSmp.nz;

	//	if((OldNe != OutStokes.ne) || (OldNx != OutStokes.nx) || (OldNz != OutStokes.nz)) return srTSend::ModifyStokesNeNxNz(OutStokes);
	//	else return 0;
	//	//DLL_IMPLEMENT
	//}

};

//*************************************************************************

#endif
