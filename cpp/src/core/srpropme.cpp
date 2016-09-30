/************************************************************************//**
 * File: srpropme.cpp
 * Description: Feasilibility tests for computation of emission and propagation of Partially-Coherent SR from Finite-Emittance Electron Beam 
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srpropme.h"
#include "srsend.h"

//*************************************************************************

CGenMathRand srTPropagMultiE::gRandGen;

//*************************************************************************

int srTPropagMultiE::PropagateElecFieldStokesAuto(srTEbmDat& ThickEbmDat, srTSRWRadStructAccessData& InWfr, srTGenOptElemHndl OptHndl, double* pPrecPar, srTStokesStructAccessData& OutStokes)
{
	int result = 0;
	//double RelPrec = 0.001; // corresponding to PrecParMultiE = 1
	//long MaxAmOfMacroPrt = 1000000;
	long long MaxAmOfMacroPrt = 1000000;
	//int AmOfSecurityPasses = 1;
	gRandGen.Initialize();

	double PrecParMultiE = 1.;
	if(pPrecPar != 0) PrecParMultiE = *pPrecPar;
	//if((PrecParMultiE > 0.) && (PrecParMultiE != 1.)) RelPrec /= PrecParMultiE;

	if(pPrecPar[1] >= 1.) MaxAmOfMacroPrt = (long)pPrecPar[1];

	//srTEbmDat ThickEbmDat;
	//if(result = InWfr.OutElectronBeamStruct(ThickEbmDat)) return result;
	srTSigleElecVars SigleElecVars(ThickEbmDat);

	srTRadResizeVect RadResizeVect, DummyResizeVect;
	int AutoMethNo = 2; //, ManMethNo = 0;

//make single-e propagation in autom. mode to set initial ranges

	srTSRWRadStructAccessData* pLocWfr = new srTSRWRadStructAccessData(&InWfr);
	pLocWfr->DoNotResizeAfter = true; 
	// this is a temporary measure, to prevent wfr range from shrinking
	// in the future, estimate thick beam wfr range here

	srTParPrecWfrPropag AutoParPrecWfrPropag(AutoMethNo, 1, 0, PrecParMultiE, 0.5);
	//srTParPrecWfrPropag ManParPrecWfrPropag(ManMethNo, 0, 0, 1., 0.5);

	srTGenOptElem *pOptElem = (srTGenOptElem*)(OptHndl.rep);
	if(result = pOptElem->PropagateRadiation(pLocWfr, AutoParPrecWfrPropag, RadResizeVect)) return result;

	if(result = ReallocateStokesAccordingToWfr(*pLocWfr, OutStokes)) return result;
	OutStokes.ZeroStokesData();

	//double CurRelPrec = -1.;
	if(result = AddWfrToStokesWithInterpXZ(*pLocWfr, OutStokes, 0)) return result;

	//srTGenOptElemPtrList ActOptElemsList;
	//((srTGenOptElem*)(OptHndl.rep))->AddPtrOfActualOptElem(ActOptElemsList);

	//int AmOfOptElemFromResize = ((int)(((double)RadResizeVect.size() + 0.0000001)/((double)(pLocWfr->ne)))) >> 1;
	//int AmOfOptElem = ActOptElemsList.size();
	//if(AmOfOptElemFromResize != AmOfOptElem) return NON_COMPATIBLE_WAVEFRONT_AND_STOKES_STRUCTS;

	if(pLocWfr != 0) delete pLocWfr;

		double StestXc = 0., StestMxx = 0.;
		double StestXpc = 0., StestMxpxp = 0.;

//repeat propagation in "automatic" mode

	srTEbmDat CurThinEbmDat = ThickEbmDat;
	//int SecurityPassCount = 0;
	//bool PrecLevelWasNotReached = true;
	
	//srTEbmDat SmallThickEbmDat = ThickEbmDat;
	//SmallThickEbmDat.

	//for(int i=1; i<MaxAmOfMacroPrt; i++)
	for(long long i=1; i<MaxAmOfMacroPrt; i++)
	{
		pLocWfr = new srTSRWRadStructAccessData(&InWfr);
		//pLocWfr->UseStartTrToShiftAtChangingRepresToCoord = true;
		//pLocWfr->xStartTr = OutStokes.xStart;
		//pLocWfr->zStartTr = OutStokes.zStart;
		pLocWfr->DoNotResizeAfter = true;
		pLocWfr->WfrEdgeCorrShouldBeDone = 0; // ????

		SetupNextThinEbm(ThickEbmDat, SigleElecVars, CurThinEbmDat);

			//OC debug
			//CurThinEbmDat.dxds0 = 2.3e-005;
			//CurThinEbmDat.x0 = 0;
			//CurThinEbmDat.dzds0 = 3.5e-006;
			//CurThinEbmDat.z0 = 0;

		SimulateWfrFromOffAxisEbm(ThickEbmDat, CurThinEbmDat, SigleElecVars, *pLocWfr);

        if(result = pOptElem->PropagateRadiation(pLocWfr, AutoParPrecWfrPropag, RadResizeVect)) return result;
		if(result = AddWfrToStokesWithInterpXZ(*pLocWfr, OutStokes, i)) return result;

		if(pLocWfr != 0) delete pLocWfr;

			StestXc += CurThinEbmDat.z0; StestMxx += (CurThinEbmDat.z0*CurThinEbmDat.z0);
			StestXpc += CurThinEbmDat.dzds0; StestMxpxp += (CurThinEbmDat.dzds0*CurThinEbmDat.dzds0);
	}

		//double testXc = StestXc/MaxAmOfMacroPrt, testMxx = StestMxx/MaxAmOfMacroPrt;
		//double testXpc = StestXpc/MaxAmOfMacroPrt, testMxpxp = StestMxpxp/MaxAmOfMacroPrt;

	return 0;
}

//*************************************************************************

int srTPropagMultiE::PropagateElecFieldStokes(srTEbmDat& ThickEbmDat, srTSRWRadStructAccessData& InWfr, srTGenOptElemHndl OptHndl, double* pPrecPar, srTStokesStructAccessData& OutStokes)
{
	int result = 0;
	double RelPrec = 0.001; // corresponding to PrecParMultiE = 1
	//long MaxAmOfMacroPrt = 1000000;
	long long MaxAmOfMacroPrt = 1000000;
	int AmOfSecurityPasses = 1;
	gRandGen.Initialize();

	double PrecParMultiE = 1.;
	if(pPrecPar != 0) PrecParMultiE = *pPrecPar;
	if((PrecParMultiE > 0.) && (PrecParMultiE != 1.))
	{
		RelPrec /= PrecParMultiE;
	}
	//if(pPrecPar[1] >= 1.) MaxAmOfMacroPrt = (long)pPrecPar[1];
	if(pPrecPar[1] >= 1.) MaxAmOfMacroPrt = (long long)pPrecPar[1];

	//srTEbmDat ThickEbmDat;
	//if(result = InWfr.OutElectronBeamStruct(ThickEbmDat)) return result;
	srTSigleElecVars SigleElecVars(ThickEbmDat);

	srTRadResizeVect RadResizeVect, DummyResizeVect;
	int AutoMethNo = 2, ManMethNo = 0;

//make single-e propagation in autom. mode and get all pre- and post-resizing params

	srTSRWRadStructAccessData* pLocWfr = new srTSRWRadStructAccessData(&InWfr);
	pLocWfr->DoNotResizeAfter = true; 
	// this is a temporary measure, to prevent wfr range from shrinking
	// in the future, estimate thick beam wfr range here

	//srTParPrecWfrPropag AutoParPrecWfrPropag(AutoMethNo, 1, 0, 1., 0.5); //srTParPrecWfrPropag(char In_MethNo, char In_UseResBefore, char In_UseResAfter, double In_PrecFact, double In_UnderSampThresh)
	srTParPrecWfrPropag AutoParPrecWfrPropag(AutoMethNo, 1, 0, PrecParMultiE, 0.5); //srTParPrecWfrPropag(char In_MethNo, char In_UseResBefore, char In_UseResAfter, double In_PrecFact, double In_UnderSampThresh)
	srTParPrecWfrPropag ManParPrecWfrPropag(ManMethNo, 0, 0, 1., 0.5);

	//if(result = ((srTGenOptElem*)(OptHndl.rep))->PropagateRadiation(pLocWfr, AutoMethNo, RadResizeVect)) return result;
	if(result = ((srTGenOptElem*)(OptHndl.rep))->PropagateRadiation(pLocWfr, AutoParPrecWfrPropag, RadResizeVect)) return result;

	if(result = ReallocateStokesAccordingToWfr(*pLocWfr, OutStokes)) return result;
	OutStokes.ZeroStokesData();

		//srTEbmDat CurThinEbmDat = ThickEbmDat;
		//SetupNextThinEbm(ThickEbmDat, SigleElecVars, CurThinEbmDat);
		//SimulateWfrFromOffAxisEbm(ThickEbmDat, CurThinEbmDat, SigleElecVars, InWfr);
		//InWfr.UseStartTrToShiftAtChangingRepresToCoord = true;
		//InWfr.xStartTr = OutStokes.xStart;
		//InWfr.zStartTr = OutStokes.zStart;
		//OptHndl.rep->PropagateRadiation(&InWfr, ManMethNo, DummyResizeVect);

	double CurRelPrec = -1.;
	long MacroPartCount = 0;
	if(result = AddWfrToStokes(*pLocWfr, OutStokes, ++MacroPartCount, CurRelPrec)) return result;
	//++MacroPartCount;

	srTGenOptElemPtrList ActOptElemsList;
	//OptHndl.rep->AddPtrOfActualOptElem(ActOptElemsList);
	((srTGenOptElem*)(OptHndl.rep))->AddPtrOfActualOptElem(ActOptElemsList);

	int AmOfOptElemFromResize = ((int)(((double)RadResizeVect.size() + 0.0000001)/((double)(pLocWfr->ne)))) >> 1;
	int AmOfOptElem = (int)ActOptElemsList.size();
	if(AmOfOptElemFromResize != AmOfOptElem) return NON_COMPATIBLE_WAVEFRONT_AND_STOKES_STRUCTS;

	if(pLocWfr != 0) delete pLocWfr;

		double StestXc = 0., StestMxx = 0.;
		double StestXpc = 0., StestMxpxp = 0.;

//repeat propagation in "manual" mode with resize params found for different coord/angle offsets

	srTEbmDat CurThinEbmDat = ThickEbmDat;
	int SecurityPassCount = 0;
	bool PrecLevelWasNotReached = true;
	while(PrecLevelWasNotReached)
	{
		pLocWfr = new srTSRWRadStructAccessData(&InWfr);

		pLocWfr->UseStartTrToShiftAtChangingRepresToCoord = true;
		pLocWfr->xStartTr = OutStokes.xStart;
		pLocWfr->zStartTr = OutStokes.zStart;
		pLocWfr->DoNotResizeAfter = true;

		SetupNextThinEbm(ThickEbmDat, SigleElecVars, CurThinEbmDat);
		SimulateWfrFromOffAxisEbm(ThickEbmDat, CurThinEbmDat, SigleElecVars, *pLocWfr);

		int ResCount = 0;
		for(srTGenOptElemPtrList::iterator iter = ActOptElemsList.begin(); iter != ActOptElemsList.end(); ++iter)
		{
			srTRadResize ResBefore = RadResizeVect[ResCount++];
			srTRadResize ResAfter = RadResizeVect[ResCount++];
			
			srTGenOptElem *pOptElem = *iter;

			pLocWfr->WfrEdgeCorrShouldBeDone = 0; // ????

			if(result = pOptElem->RadResizeGen(*pLocWfr, ResBefore)) return result;

			//if(result = pOptElem->PropagateRadiation(pLocWfr, ManMethNo, DummyResizeVect)) return result;
			if(result = pOptElem->PropagateRadiation(pLocWfr, ManParPrecWfrPropag, DummyResizeVect)) return result;
			
			//if(result = pOptElem->MakePostPropagationProc(pLocWfr, ResAfter)) return result;
		}
		if(result = AddWfrToStokes(*pLocWfr, OutStokes, ++MacroPartCount, CurRelPrec)) return result;
		//PrecLevelWasNotReached = (CurRelPrec > RelPrec);

		if(CurRelPrec <= RelPrec)
		{
			if(SecurityPassCount >= AmOfSecurityPasses) PrecLevelWasNotReached = false;
			else SecurityPassCount++;
		}
		else SecurityPassCount = 0;

		if(pLocWfr != 0) delete pLocWfr;

			StestXc += CurThinEbmDat.z0; StestMxx += (CurThinEbmDat.z0*CurThinEbmDat.z0);
			StestXpc += CurThinEbmDat.dzds0; StestMxpxp += (CurThinEbmDat.dzds0*CurThinEbmDat.dzds0);

		if(MacroPartCount >= MaxAmOfMacroPrt) PrecLevelWasNotReached = false;
	}

		//double testXc = StestXc/MaxAmOfMacroPrt, testMxx = StestMxx/MaxAmOfMacroPrt;
		//double testXpc = StestXpc/MaxAmOfMacroPrt, testMxpxp = StestMxpxp/MaxAmOfMacroPrt;

	return 0;
}

//*************************************************************************

int srTPropagMultiE::EmitPropagElecFieldStokes(srTTrjDat& TrjDat, srTGenOptElemHndl, double* pPrecPar, srTStokesStructAccessData& OutStokes)
{

//ddddddddddddddddddddddddddddddd


	return 0;
}

//*************************************************************************

int srTPropagMultiE::ReallocateStokesAccordingToWfr(srTSRWRadStructAccessData& LocWfr, srTStokesStructAccessData& OutStokes)
{
	OutStokes.eStart = LocWfr.eStart;
	OutStokes.eStep = LocWfr.eStep;
	OutStokes.ne = LocWfr.ne;
	OutStokes.zStart = LocWfr.zStart;
	OutStokes.zStep = LocWfr.zStep;
	OutStokes.nz = LocWfr.nz;
	OutStokes.xStart = LocWfr.xStart;
	OutStokes.xStep = LocWfr.xStep;
	OutStokes.nx = LocWfr.nx;

	srTSend Send;
	return Send.ModifyStokesNeNxNz(OutStokes);
	//DLL_IMPLEMENT
	return 0;
}

//*************************************************************************

//int srTPropagMultiE::AddWfrToStokesWithInterpXZ(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long MacroPartCountZeroBased)
int srTPropagMultiE::AddWfrToStokesWithInterpXZ(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long long MacroPartCountZeroBased)
{
	//int result = 0;
	//if(AllowResize)
	//{
	//	//check and resize stokes, if necessary
	//}

	//The following assumes that Stokes have necessary ranges

	float fInvN = (float)(1./(MacroPartCountZeroBased + 1));
	float fMacroPartCountZeroBased = (float)MacroPartCountZeroBased;

	//long PerE = 4;
	//long PerX = Stokes.ne*PerE;
	//long PerZ = Stokes.nx*PerX;
	long long PerE = 4;
	long long PerX = Stokes.ne*PerE;
	long long PerZ = Stokes.nx*PerX;

	srTEXZ EXZ;
    EXZ.z = Stokes.zStart;

	for(int iz=0; iz<Stokes.nz; iz++)
	{
		//long izPerZ = iz*PerZ;
		long long izPerZ = iz*PerZ;
		EXZ.x = Stokes.xStart;
		for(int ix=0; ix<Stokes.nx; ix++)
		{
			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;
			EXZ.e = Stokes.eStart;
			for(int ie=0; ie<Stokes.ne; ie++)
			{
				float *pSto = Stokes.pBaseSto + (izPerZ + ixPerX +ie*PerE);

				float *tSto = pSto;
				for(int k=0; k<4; k++) *(tSto++) *= fMacroPartCountZeroBased; 
				Wfr.AddStokesAtPoint(EXZ, pSto);
				tSto = pSto;
				for(int i=0; i<4; i++) *(tSto++) *= fInvN; 

				//if(result = srYield.Check()) return result;
				//if(result = CompProgressInd.UpdateIndicator(PointCount++)) return result;
				
				EXZ.e += Stokes.eStep;
			}
			EXZ.x += Stokes.xStep;
		}
		EXZ.z += Stokes.zStep;
	}
	return 0;
}

//*************************************************************************

//int srTPropagMultiE::AddWfrToStokes(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long MacroPartCount, double& CurRelPrec)
int srTPropagMultiE::AddWfrToStokes(srTSRWRadStructAccessData& Wfr, srTStokesStructAccessData& Stokes, long long MacroPartCount, double& CurRelPrec)
{
	if((Wfr.ne != Stokes.ne) || (Wfr.nx != Stokes.nx) || (Wfr.nz != Stokes.nz)) return NON_COMPATIBLE_WAVEFRONT_AND_STOKES_STRUCTS;

	//long AmOfPts = Stokes.ne*Stokes.nx*Stokes.nz;
	long long AmOfPts = ((long long)Stokes.ne)*((long long)Stokes.nx)*((long long)Stokes.nz);
	float *tSto = Stokes.pBaseSto;
	float *tEx = Wfr.pBaseRadX, *tEz = Wfr.pBaseRadZ; 

	float StokesVal[4];
	double InvN = 1./MacroPartCount;
	//long Nmi1 = MacroPartCount - 1;
	long long Nmi1 = MacroPartCount - 1;
	double SqDif = 0.;
	//for(long i=0; i<AmOfPts; i++)
	for(long long i=0; i<AmOfPts; i++)
	{
		CalcStokesFromE(tEx, tEz, StokesVal);
		
		float OldS0 = *tSto;
		float *tStoVal = StokesVal;
		*tSto *= (float)Nmi1; *tSto += (float)(*(tStoVal++)); *tSto *= (float)InvN;
		double ds = (*tSto - OldS0);
		if(*tSto == 0.) ds = 0.;
		else ds /= (*tSto);

		SqDif += ds*ds;

		tSto++;
		*tSto *= (float)Nmi1; *tSto += (float)(*(tStoVal++)); *(tSto++) *= (float)InvN;
		*tSto *= (float)Nmi1; *tSto += (float)(*(tStoVal++)); *(tSto++) *= (float)InvN;
		*tSto *= (float)Nmi1; *tSto += (float)(*tStoVal); *(tSto++) *= (float)InvN;

		//*tSto = ((*tSto)*Nmi1 + (*StokesVal))*InvN;
		tEx += 2; tEz += 2;
	}
	CurRelPrec = sqrt(SqDif/AmOfPts);
	return 0;
}

//*************************************************************************

void srTPropagMultiE::SetupNextThinEbm(srTEbmDat& EbmDat, srTSigleElecVars& SigleElecVars, srTEbmDat& OutThinEbmDat)
{
	//Take into account here:
	// - coupling x - xp
	// - distributions other then Gaussian
	// - energy spread

	bool VarZ_or_VarZp = (SigleElecVars.VarZ || SigleElecVars.VarZp);
	bool VarX_or_VarXp = (SigleElecVars.VarX || SigleElecVars.VarXp);

	if((!VarX_or_VarXp) && VarZ_or_VarZp)
	{
		//gRandGen.NextGaussRand2D(EbmDat.z0, sqrt(EbmDat.Mzz), EbmDat.dzds0, sqrt(EbmDat.Mzpzp), OutThinEbmDat.z0, OutThinEbmDat.dzds0);
		gRandGen.NextRandGauss2D(EbmDat.z0, sqrt(EbmDat.Mzz), EbmDat.dzds0, sqrt(EbmDat.Mzpzp), OutThinEbmDat.z0, OutThinEbmDat.dzds0);

		//OutThinEbmDat.z0 = 0.; OutThinEbmDat.dzds0 = 2.e-06;
	}
	else if(VarX_or_VarXp && (!VarZ_or_VarZp))
	{
		//gRandGen.NextGaussRand2D(EbmDat.x0, sqrt(EbmDat.Mxx), EbmDat.dxds0, sqrt(EbmDat.Mxpxp), OutThinEbmDat.x0, OutThinEbmDat.dxds0);
		gRandGen.NextRandGauss2D(EbmDat.x0, sqrt(EbmDat.Mxx), EbmDat.dxds0, sqrt(EbmDat.Mxpxp), OutThinEbmDat.x0, OutThinEbmDat.dxds0);
	}
	else if(VarX_or_VarXp && VarZ_or_VarZp)
	{
		double XcArr[] = {EbmDat.x0, EbmDat.dxds0, EbmDat.z0, EbmDat.dzds0};
		double SigmaXArr[] = {sqrt(EbmDat.Mxx), sqrt(EbmDat.Mxpxp), sqrt(EbmDat.Mzz), sqrt(EbmDat.Mzpzp)};

		//gRandGen.NextGaussRand4D(XcArr, SigmaXArr, OutThinEbmDat.x0, OutThinEbmDat.dxds0, OutThinEbmDat.z0, OutThinEbmDat.dzds0);
		gRandGen.NextRandGauss4D(XcArr, SigmaXArr, OutThinEbmDat.x0, OutThinEbmDat.dxds0, OutThinEbmDat.z0, OutThinEbmDat.dzds0);
		
		//int aha = 1;
		//make sure Multi-D LpTau is implemented
	}

/**
	if(SigleElecVars.VarX)
	{
		OutThinEbmDat.x0 = gRandGen.NextGaussRand(EbmDat.x0, sqrt(EbmDat.Mxx));
	}
	if(SigleElecVars.VarXp)
	{
		OutThinEbmDat.dxds0 = gRandGen.NextGaussRand(EbmDat.dxds0, sqrt(EbmDat.Mxpxp));
	}
	if(SigleElecVars.VarZ)
	{
		OutThinEbmDat.z0 = gRandGen.NextGaussRand(EbmDat.z0, sqrt(EbmDat.Mzz));
	}
	if(SigleElecVars.VarZp)
	{
		OutThinEbmDat.dzds0 = gRandGen.NextGaussRand(EbmDat.dzds0, sqrt(EbmDat.Mzpzp));
	}
**/
}

//*************************************************************************

void srTPropagMultiE::SimulateWfrFromOffAxisEbm(srTEbmDat& OnAxisEbmDat, srTEbmDat& OffAxisEbmDat, srTSigleElecVars& SigleElecVars, srTSRWRadStructAccessData& Wfr)
{
	if(SigleElecVars.VarX)
	{
		double dx = OffAxisEbmDat.x0 - OnAxisEbmDat.x0;
		Wfr.xStart += dx;
		Wfr.xc += dx;
	}
	if(SigleElecVars.VarZ)
	{
		double dz = OffAxisEbmDat.z0 - OnAxisEbmDat.z0;
		Wfr.zStart += dz;
		Wfr.zc += dz;
	}

	double AngX = OffAxisEbmDat.dxds0 - OnAxisEbmDat.dxds0;
	double AngZ = OffAxisEbmDat.dzds0 - OnAxisEbmDat.dzds0;

	if(SigleElecVars.VarXp || SigleElecVars.VarZp)
	{
		//double kMult = 5.0676816037856e+06; //(TwoPi/(1.239854e-06))
		double kMult = 5.067730652e+06; //(TwoPi/(1.239842e-06))

		float *pEx0 = Wfr.pBaseRadX;
		float *pEz0 = Wfr.pBaseRadZ;
		//long PerX = Wfr.ne << 1;
		//long PerZ = PerX*Wfr.nx;
		long long PerX = Wfr.ne << 1;
		long long PerZ = PerX*Wfr.nx;

		srTEFieldPtrs EFieldPtrs;
		srTEXZ EXZ;
		EXZ.z = Wfr.zStart;
		//long izPerZ = 0;
		long long izPerZ = 0;

		for(int iz=0; iz<Wfr.nz; iz++)
		{
			float *pEx_StartForX = pEx0 + izPerZ;
			float *pEz_StartForX = pEz0 + izPerZ;
			EXZ.x = Wfr.xStart;
			//long ixPerX = 0;
			long long ixPerX = 0;
		
			for(int ix=0; ix<Wfr.nx; ix++)
			{
				float *pEx_StartForE = pEx_StartForX + ixPerX;
				float *pEz_StartForE = pEz_StartForX + ixPerX;
				EXZ.e = Wfr.eStart;
				//long iePerE = 0;
				long long iePerE = 0;

				for(int ie=0; ie<Wfr.ne; ie++)
				{
					// e in eV; Length in m !!!
					double WaveNumb = EXZ.e*kMult;
					double Ph = WaveNumb*(AngX*EXZ.x + AngZ*EXZ.z);
					double CosPh = cos(Ph), SinPh = sin(Ph);
					//check absolute/relative angle here

					if(pEx0 != 0)
					{
						float *pExRe = pEx_StartForE + iePerE;
						float *pExIm = pExRe + 1;

						float NewExRe = (float)((*pExRe)*CosPh - (*pExIm)*SinPh);
						float NewExIm = (float)((*pExRe)*SinPh + (*pExIm)*CosPh);
						*pExRe = NewExRe; *pExIm = NewExIm;

						//if(*pExRe != 0.)
						//{
						//	int aha = 1;
						//}
					}
					if(pEz0 != 0)
					{
						float *pEzRe = pEz_StartForE + iePerE;
						float *pEzIm = pEzRe + 1;

						float NewEzRe = (float)((*pEzRe)*CosPh - (*pEzIm)*SinPh);
						float NewEzIm = (float)((*pEzRe)*SinPh + (*pEzIm)*CosPh);
						*pEzRe = NewEzRe; *pEzIm = NewEzIm;
					}

					iePerE += 2;
					EXZ.e += Wfr.eStep;
				}
				ixPerX += PerX;
				EXZ.x += Wfr.xStep;
			}
			izPerZ += PerZ;
			EXZ.z += Wfr.zStep;
		}

		if(SigleElecVars.VarXp)
		{
			double dx = Wfr.RobsX*AngX;
			Wfr.xStart += dx;
			Wfr.xWfrMin += dx; Wfr.xWfrMax += dx;
			Wfr.xc += dx;
		}
		if(SigleElecVars.VarZp)
		{
			double dz = Wfr.RobsZ*AngZ;

				//double dz = 50.e-06;

			Wfr.zStart += dz;
			Wfr.zWfrMin += dz; Wfr.zWfrMax += dz;
			Wfr.zc += dz;
		}
	}
/*
Variable WaveNumb = (2*Pi)/(6.27774e-11)
Variable dxp = 0.000004
Variable ph = WaveNumb*dxp*x
return cos(ph) + cmplx(0,1)*sin(ph)
*/
}

//*************************************************************************

int srTPropagMultiE::ReallocateStokesAccordingToWfr(srTWfrSmp& WfrSmp, srTStokesStructAccessData& OutStokes)
{
	int OldNe = OutStokes.ne, OldNx = OutStokes.nx, OldNz = OutStokes.nz;

	OutStokes.eStart = WfrSmp.LambStart;
	OutStokes.eStep = 0;
	if(WfrSmp.nLamb > 1) OutStokes.eStep = (WfrSmp.LambEnd - WfrSmp.LambStart)/WfrSmp.nLamb;
	OutStokes.ne = WfrSmp.nLamb;

	OutStokes.xStart = WfrSmp.xStart;
	OutStokes.xStep = 0;
	if(WfrSmp.nx > 1) OutStokes.xStep = (WfrSmp.xEnd - WfrSmp.xStart)/WfrSmp.nx;
	OutStokes.nx = WfrSmp.nx;

	OutStokes.zStart = WfrSmp.zStart;
	OutStokes.zStep = 0;
	if(WfrSmp.nz > 1) OutStokes.eStep = (WfrSmp.zEnd - WfrSmp.zStart)/WfrSmp.nz;
	OutStokes.nz = WfrSmp.nz;

	if((OldNe != OutStokes.ne) || (OldNx != OutStokes.nx) || (OldNz != OutStokes.nz)) return srTSend::ModifyStokesNeNxNz(OutStokes);
	else return 0;
	//DLL_IMPLEMENT
}

//*************************************************************************
