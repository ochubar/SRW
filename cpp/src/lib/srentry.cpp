/************************************************************************//**
 * File: srentry.cpp
 * Description: C/C++ API of Synchrotron Radiation Workshop (SRW) - OBSOLETE version
 * Project: Synchrotron Radiation Workshop
 * First release: October 2002
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srinterf.h"

#include "srigorre.h"
#include "stlstart.h"
#include "objcont.h"
#include "srobject.h"
#include "srgtrjdt.h"
#include "srstraux.h"
#include "srerror.h"
#include "srctrjdt.h"
#include "srgsnbm.h"
#include "srmagcnt.h"
#include "srradint.h"
#include "srpowden.h"
#include "srradmnp.h"
#include "sroptdrf.h"
//#include "srmamet.h"
#include "gmmeth.h"
#include "gmfunc.h"
#include "srcradint.h"
#include "auxparse.h"

#include "srwlib.h" //from SRWLIB, to allow for compiling

//-------------------------------------------------------------------------

#ifndef __IGOR_PRO__
#ifdef _WINDOWS
int SpinProcess() { return 0;}
Handle NewHandle(long size)
{//to check
	if(size <= 0) return 0;
	char** aHandle = new char*();
	*aHandle = new char[size];
	return aHandle;
}
#endif
#endif

//-------------------------------------------------------------------------
// Global Variables

CObjCont<CGenObject> gSRObjects;

#ifdef __IGOR_PRO__
extern srTIntVect gVectWarnNos;
#else
srTIntVect gVectWarnNos;
srTYield srYield;
int gCallSpinProcess = 1;
#endif

CGenMathRand gMathRand; //random generator object

//-------------------------------------------------------------------------

int (*pgCompProgressIndicFunc)(double CurVal);
int (*pgWfrExtModifFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp);
int (*pgOptElemGetInfByNameFunc)(const char* sNameOptElem, char** pDescrStr, int* LenDescr, void*);

int (*gpWfrModifFunc)(int action, SRWLWfr* pWfrIn, char pol) = 0; //from SRWLIB, to allow for compiling

//-------------------------------------------------------------------------

void UtiWarnCheck()
{
	if(!gVectWarnNos.empty())
	{
		int CurWarnNo = gVectWarnNos[0];
        gVectWarnNos.erase(gVectWarnNos.begin(), gVectWarnNos.end());
		throw CurWarnNo;
	}
}

//-------------------------------------------------------------------------

//void UtiAuxFillStringVect(srTStringVect& AuxStringVect, char** pDescrStr, int LenDescr)
//{//to move to a dedicated class ?
//	if((pDescrStr == 0) || (LenDescr == 0)) return;
//	for(int i=0; i<LenDescr; i++) 
//	{
//		char* CurStr = pDescrStr[i];
//		AuxStringVect.push_back(CurStr);
//	}
//}

//-------------------------------------------------------------------------

EXP int CALL srUtiVerNo(char* VerNoStr)
{//to modify at each new release!
	char CurVerStr[] = "3.965"; //"3.87";

	strcpy(VerNoStr, CurVerStr);
	return 0;
}

//-------------------------------------------------------------------------

EXP void CALL srUtiProgrIndFuncSet(int (*pExtFunc)(double CurVal))
{
	if(pExtFunc != 0) pgCompProgressIndicFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP void CALL srUtiWfrExtModifFuncSet(int (*pExtFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp))
{
	if(pExtFunc != 0) pgWfrExtModifFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP void CALL srUtiOptElemGetInfByNameFuncSet(int (*pExtFunc)(const char*, char**, int*, void*))
{
	if(pExtFunc != 0) pgOptElemGetInfByNameFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP void CALL srUtiWarnTextGet(char* t, int WarnNo)
{
	//const char* strWarn = CErrWarn::GetWarning(WarnNo);
	
	CErrWarn srwlErWar; //OC030412
	const char* strWarn = srwlErWar.GetWarning(WarnNo);
	if(strWarn != 0) strcpy(t, strWarn);
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamSet(int* i, double I, double Neb, double* pMom1, int nMom1, double* pMom2, int nMom2, double s)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTEbmDat(I, Neb, pMom1, nMom1, pMom2, nMom2, s));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamContSet(int* i, int* pElecElem, int nElecElem)
{
	try 
	{
		//dddddddddddddddddddddddd
        //*i = gSRObjects.insert(new srTEbmDat(I, Neb, pMom1, nMom1, pMom2, nMom2, s));
		//
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamPropagate(int iElecBeam, int iOptElem)
{// modify to support magnetic optical elements
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iOptElem))) return SR_OBJECT_DOES_NOT_EXIST;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;

	srTDriftSpace* pDriftSpace = dynamic_cast<srTDriftSpace*>(gSRObjects.getPtr(iOptElem));
    if(pDriftSpace == 0) return OBJECT_IS_NOT_DRIFT_SPACE;

	try 
	{ 
		//int res = 0;
		//if(res = pDriftSpace->PropagateElecBeamMoments(pElecBeam)) return res;
		int res = pDriftSpace->PropagateElecBeamMoments(pElecBeam); //OC150208
		if(res > 0) throw res;
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamTwiss2Mom(double* pMom2, int nMom2, double* pTwissHor, double* pTwissVer, double RmsEnSpr)
{
	//todo: conversion

	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamMom2EmitAndTwiss(double* pEmit, double* pBeta, double* pAlpha, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double Eta, double EtaPrime)
{
	return srTEbmDat::Mom2EmitAndTwiss(pEmit, pBeta, pAlpha, Sigma2, SigmaPrime2, MixMom, SigmaE2, Eta, EtaPrime);
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamMomAndEmit2Twiss(double* pBeta, double* pAlpha, double* pEta, double* pEmit, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double EtaPrime)
{
	return srTEbmDat::MomAndEmit2Twiss(pBeta, pAlpha, pEta, pEmit, Sigma2, SigmaPrime2, MixMom, SigmaE2, EtaPrime);
}

//-------------------------------------------------------------------------

EXP int CALL srElecBeamGet(double* pI, double* pMom1, double* pMom2, double* ps0, int iElecBeam)
{
	if(!gSRObjects.exists(iElecBeam)) return SR_OBJECT_DOES_NOT_EXIST;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;

	try 
	{ 
		if(pI != 0) *pI = pElecBeam->Current;
		if(ps0 != 0) *ps0 = pElecBeam->s0;
		if(pMom1 != 0) pElecBeam->GetMom1(pMom1);
		if(pMom2 != 0) pElecBeam->GetMom2(pMom2);
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------
/**
EXP int CALL srMagFldPerHarmSet(int* i, int n, char v_or_h, double K, double Ph)
{
	if(n <= 0) return INCORRECT_MAG_HARM_NUMBER;
	try 
	{ 
		char XorZ = ((v_or_h == 'v')? 'z' : 'x');
		*i = gSRObjects.insert(new srTMagHarm(n, XorZ, K, Ph));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}
**/
//-------------------------------------------------------------------------

EXP int CALL srMagFldPerSet(int* i, double Per, double L, double s0, int nHarm, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase, char Type, double SpecPar)
{
	if(Per <= 0) return INCORRECT_MAG_FIELD_PER;
	if(L <= 0) return INCORRECT_PER_MAG_FIELD_LENGTH;
	//if(nHarm == 0) return INCORRECT_NUMBER_OF_HARMONICS;

	try 
	{
		srTMagHarm* HarmArr = 0;
		if((nHarm > 0) && (ArrHarmNo != 0) && (ArrXorZ != 0) && (ArrK != 0))
		{
			//bool AssumeZeroPhase = (ArrPhase == 0);
			HarmArr = new srTMagHarm[nHarm];
			srTMagHarm *tHarmArr = HarmArr;

			int *tArrHarmNo = ArrHarmNo;
			char *tArrXorZ = ArrXorZ;
			double *tArrK = ArrK;
			double *tArrPhase = ArrPhase;

			for(int k=0; k<nHarm; k++)
			{
				(tHarmArr++)->InData(*(tArrHarmNo++), *(tArrXorZ++), *(tArrK++), *(tArrPhase++));
			}
		}
		*i = gSRObjects.insert(new srTMagFieldPeriodic(Per, L, s0, HarmArr, nHarm, Type, SpecPar));

		if(HarmArr != 0) delete[] HarmArr;
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------
/**
EXP int CALL srMagFldPerSet(int* i, double Per, double L, double s0, int* pHarm, int nHarm, char Type, double SpecPar)
{
	if(Per <= 0) return INCORRECT_MAG_FIELD_PER;
	if(L <= 0) return INCORRECT_PER_MAG_FIELD_LENGTH;
	//if((pHarm == 0) || (nHarm == 0)) return INCORRECT_NUMBER_OF_HARMONICS;

	try 
	{
		srTMagHarm* HarmArr = 0;

		if((pHarm != 0) && (nHarm > 0))
		{
			HarmArr = new srTMagHarm[nHarm];
			for(int k=0; k<nHarm; k++)
			{
				srTMagHarm* tHarm = dynamic_cast<srTMagHarm*>(gSRObjects.getPtr(pHarm[k]));
				if(tHarm != 0) HarmArr[k] = *tHarm;
			}
		}
		*i = gSRObjects.insert(new srTMagFieldPeriodic(Per, L, s0, HarmArr, nHarm, Type, SpecPar));

		if(HarmArr != 0) delete[] HarmArr;
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}
**/
//-------------------------------------------------------------------------

EXP int CALL srMagFldPerSetFromTrUnif(int* i, int iMagFld, double RelPrec, int MaxHarm, double MaxPerLen_m)
{
	if(RelPrec <= 0) return INCORRECT_RELATIVE_PRECISION;
	if(MaxHarm <= 0) return INCORRECT_NUMBER_OF_HARMONICS;

	if(!gSRObjects.exists(iMagFld)) return SR_OBJECT_DOES_NOT_EXIST;
	
	srTMagFldTrUnif* pMagFldTrUnif = dynamic_cast<srTMagFldTrUnif*>(gSRObjects.getPtr(iMagFld));
    if(pMagFldTrUnif == 0) return OBJECT_IS_NOT_MAG_TRANSV_UNIF;

	try 
	{
		*i = gSRObjects.insert(pMagFldTrUnif->CreateAndSetupMagFieldPeriodic(RelPrec, MaxHarm, MaxPerLen_m));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFldPerGet(int iMagFldPer, double* pPer, double* pL, double* psCen, int* pnHarm, char* pType, double* pSpecPar, double* pkeVperGeV2)
{
    if(!gSRObjects.exists(iMagFldPer)) return SR_OBJECT_DOES_NOT_EXIST;

	CHGenObj hMagFldPer = gSRObjects.get(iMagFldPer);
	srTMagFieldPeriodic* pMagFldPer = dynamic_cast<srTMagFieldPeriodic*>(hMagFldPer.ptr());
    if(pMagFldPer == 0) return OBJECT_IS_NOT_MAG_PER;

	try 
	{
		pMagFldPer->OutBasicData(pPer, pL, psCen, pnHarm, pType, pSpecPar, pkeVperGeV2);
        UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFldPerHarmGet(int iMagFldPer, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase)
{
    if(!gSRObjects.exists(iMagFldPer)) return SR_OBJECT_DOES_NOT_EXIST;

	CHGenObj hMagFldPer = gSRObjects.get(iMagFldPer);
	srTMagFieldPeriodic* pMagFldPer = dynamic_cast<srTMagFieldPeriodic*>(hMagFldPer.ptr());
    if(pMagFldPer == 0) return OBJECT_IS_NOT_MAG_PER;

	try 
	{
		pMagFldPer->OutHarmData(ArrHarmNo, ArrXorZ, ArrK, ArrPhase);
        UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFldConstSet(int* i, double BH, double BV)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTMagFieldConstant(BH, BV));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFldTrUnifSet(int* i, double sStart, double sStep, int np, double* pBh, double* pBv)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTMagFldTrUnif(sStart, sStep, np, pBh, pBv, 1));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFld3dSet(int* i, double xStart, double xStep, int nx, double yStart, double yStep, int ny, double zStart, double zStep, int nz, double* pBx, double* pBy, double* pBz)
{
	try 
	{
		//*i = gSRObjects.insert(new srTMagFld3d(xStart, xStep, nx, yStart, yStep, ny, zStart, zStep, nz, pBx, pBy, pBz, 1));
		//*i = gSRObjects.insert(new srTMagFld3d(xStart, xStep, nx, yStart, yStep, ny, zStart, zStep, nz, pBx, pBy, pBz, 1, 1));
		*i = gSRObjects.insert(new srTMagFld3d(xStart, xStep, nx, yStart, yStep, ny, zStart, zStep, nz, pBx, pBy, pBz, 1, 1, 1));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srMagFldContSet(int* i, int* pMagElem, int nMagElem)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTMagFldCont(pMagElem, nMagElem));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srGausBeamSet(int* i, double SpecFlux, int Polar, double sigX, int mx, double sigZ, int mz, double sigT, int typeT, double* pMom1, double s0, double RepRate, double PulseEn, double AvgPhotEn)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTGsnBeam(SpecFlux, Polar, sigX, mx, sigZ, mz, sigT, typeT, pMom1, s0, RepRate, PulseEn, AvgPhotEn));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrSmpSet(int* i, double s, double hSt, double hFi, int hN, double vSt, double vFi, int vN, double* pSurfData, double eSt, double eFi, int eN, char* PhotEnUnit, double tSt, double tFi, int tN, int presT, double* horOrtObsPlane, double* inNormObsPlane)
{
	try 
	{ 
		*i = gSRObjects.insert(new srTWfrSmp(s, hSt, hFi, hN, vSt, vFi, vN, pSurfData, eSt, eFi, eN, PhotEnUnit, tSt, tFi, tN, presT, horOrtObsPlane, inNormObsPlane));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrSet(int* piWfr, srTSRWRadInData* pSRWRadInData, int CopyData)
{
	try 
	{ 
		*piWfr = gSRObjects.insert(new srTSRWRadStructAccessData(pSRWRadInData, CopyData));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srOptElemSet(int* piOptElem, char** pDescrStr, int* pLenDescr, srTDataMD* pExtraData, int iWfr)
{
	if((piOptElem == 0) || (pDescrStr == 0)) return INCORRECT_ARGUMENTS;
	//ATTENTION: it may return updated strings (if they are allocated outside)

	srTSRWRadStructAccessData *pWfr = 0;
	if(iWfr != 0)
	{
		CHGenObj hWfr = gSRObjects.get(iWfr);
		pWfr = dynamic_cast<srTSRWRadStructAccessData*>(hWfr.ptr());
		if(pWfr == 0) return OBJECT_IS_NOT_WAVEFRONT;
	}

	try 
	{
		srTStringVect AuxStringVect;
		//UtiAuxFillStringVect(AuxStringVect, pDescrStr, LenDescr);
		CAuxParse::StringArr2VectCStr(pDescrStr, *pLenDescr, AuxStringVect);

		srTGenOptElemHndl OptElemHndl;
		//int res;
		//if(res = srTGenOptElem::SetupOpticalElement(&AuxStringVect, pExtraData, pWfr, OptElemHndl)) throw res;
		int res = srTGenOptElem::SetupOpticalElement(&AuxStringVect, pExtraData, pWfr, OptElemHndl);
		if(res > 0) throw res;
		//ATTENTION: it returns eventually updated strings!

		*piOptElem = gSRObjects.insert(OptElemHndl);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

//EXP int CALL srOptElemGetInf(char** pDescrStr, int* pLenDescr, int iOptElem)
//{
//	if((pDescrStr == 0) || (pLenDescr == 0)) return 0;
//	if(!gSRObjects.exists(iOptElem)) return SR_OBJECT_DOES_NOT_EXIST;
//
//	srTGenOptElem* pOptElem = dynamic_cast<srTGenOptElem*>(gSRObjects.getPtr(iOptElem));
//    if(pOptElem == 0) return OBJECT_IS_NOT_OPT_ELEM;
//
//	try 
//	{
//		pOptElem->OutOptElemInfo(pDescrStr, pLenDescr);
//		UtiWarnCheck();
//	}
//	catch(int ErrNo) { return ErrNo;}
//	return 0;
//}

//-------------------------------------------------------------------------

EXP int CALL srOptDriftSet(int* piOptElem, double Length)
{
	try 
	{ 
		*piOptElem = gSRObjects.insert(new srTDriftSpace(Length));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srOptContSet(int* piOptElem, int* piMembArr, int nMemb)
{
	try 
	{ 
		//ddddddddddddddddddddddddddddddddd

		//*piOptElem = gSRObjects.insert(new srTDriftSpace(Length));
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecTrjComp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep, int iElecBeam, int iMagFld)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld))) return SR_OBJECT_DOES_NOT_EXIST;
	if((pOutBtxData == 0) && (pOutXData == 0) && (pOutBtzData == 0) && (pOutZData == 0)) return 0;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;
	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTGenTrjDat *pTrjDat = 0;

	try 
	{
		//srTTrjDat TrjDat(pElecBeam, pMagElem);
		//TrjDat.CompTrjDataForDisp(pOutBtxData, pOutXData, pOutBtzData, pOutZData);

		pTrjDat = pMagElem->CreateAndSetupNewTrjDat(pElecBeam);
		pTrjDat->CompTrjDataForDisp(pOutBtxData, pOutXData, pOutBtyData, pOutYData, pOutBtzData, pOutZData, ns, sStart, sStep);

		delete pTrjDat; pTrjDat = 0;
		UtiWarnCheck();
	}
	catch(int ErrNo) 
	{
		if(pTrjDat != 0) delete pTrjDat;
		return ErrNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecFldComp(int* piWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pVoidPrecElecFld)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pVoidPrecElecFld == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;

	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTMagFldTrUnif* pMagFldTrUnif = dynamic_cast<srTMagFldTrUnif*>(pMagElem);
    if(pMagFldTrUnif == 0) return OBJECT_IS_NOT_MAG_TRANSV_UNIF;

	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
        srTParPrecElecFld* pPrecElecFld = (srTParPrecElecFld*)pVoidPrecElecFld;
		srTRadInt RadInt;
		//srTTrjDat TrjDat(pElecBeam, pMagElem);
		srTTrjDat TrjDat(pElecBeam, pMagFldTrUnif);
		int res = TrjDat.ComputeInterpolatingStructure();
		if(res > 0) throw res;

		srTSRWRadStructAccessData* pWfr = new srTSRWRadStructAccessData(pElecBeam, &TrjDat, pWfrSmp, pPrecElecFld->NxNzOversamplingFactor);
		*piWfr = gSRObjects.insert(pWfr);
		RadInt.ComputeElectricFieldFreqDomain(&TrjDat, pWfrSmp, pPrecElecFld, pWfr);

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecFldCompExt(void* pvWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pVoidPrecElecFld)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pVoidPrecElecFld == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;

	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTMagFldTrUnif* pMagFldTrUnif = dynamic_cast<srTMagFldTrUnif*>(pMagElem);
    if(pMagFldTrUnif == 0) return OBJECT_IS_NOT_MAG_TRANSV_UNIF;

	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
        srTSRWRadInData* pRadInData = (srTSRWRadInData*)pvWfr;
        srTParPrecElecFld* pPrecElecFld = (srTParPrecElecFld*)pVoidPrecElecFld;

        //srTTrjDat TrjDat(pElecBeam, pMagElem);
        srTTrjDat TrjDat(pElecBeam, pMagFldTrUnif);
        int res = TrjDat.ComputeInterpolatingStructure();
		if(res > 0) throw res;

		//srTSRWRadStructAccessData Wfr(pRadInData, pElecBeam, &TrjDat, pWfrSmp, pPrecElecFld->NxNzOversamplingFactor);
		srTSRWRadStructAccessData Wfr(pRadInData, pElecBeam, &TrjDat, pWfrSmp, pPrecElecFld);

        srTRadInt RadInt;
		//char showProgressInd = 1; //0; //to switch for special releases
		char showProgressInd = (char)(pPrecElecFld->ShowProgrIndic); //OC180312
		RadInt.ComputeElectricFieldFreqDomain(&TrjDat, pWfrSmp, pPrecElecFld, &Wfr, showProgressInd);

		Wfr.OutSRWRadPtrs(pRadInData);

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL csrElecFldCompExt(void* pvWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pPrcPar == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	if(pvWfr == 0) return INCORRECT_PARAMS_SR_COMP_PREC;

	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;

	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTMagFldTrUnif* pMagFldTrUnif = dynamic_cast<srTMagFldTrUnif*>(pMagElem);
    if(pMagFldTrUnif == 0) return OBJECT_IS_NOT_MAG_TRANSV_UNIF;

	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
        srTSRWRadInData* pRadInData = (srTSRWRadInData*)pvWfr;
		double* pdPrcPar = (double*)pPrcPar;

        //srTTrjDat TrjDat(pElecBeam, pMagElem);
        srTTrjDat TrjDat(pElecBeam, pMagFldTrUnif);
        int res = TrjDat.ComputeInterpolatingStructure();
		if(res > 0) throw res;

		double UseAutoSamp = pdPrcPar[4], NxNzOversampFact = pdPrcPar[5];
		if(UseAutoSamp <= 0) NxNzOversampFact = 0;
		srTSRWRadStructAccessData Wfr(pRadInData, pElecBeam, &TrjDat, pWfrSmp, NxNzOversampFact);

		srTCSR csrInt(TrjDat, pdPrcPar, Wfr);
		csrInt.computeElectricFieldFreqDomain();

		Wfr.OutSRWRadPtrs(pRadInData);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;

	//try 
	//{
	//	srTSRWStokesInData* pStokesInData = (srTSRWStokesInData*)pVoidStokesInData;
	//	srTStokesStructAccessData Stokes(pStokesInData);
	//	pMagElem->ComputeSR_Stokes(pElecBeam, pWfrSmp, pPrcPar, &Stokes);
	//	Stokes.OutDataPtrs(pStokesInData);
	//	UtiWarnCheck();
	//}
}

//-------------------------------------------------------------------------

EXP int CALL srElecFldGausBeamComp(int* piWfr, int iGsnBeam, int iWfrSmp, void* pVoidPrecElecFldGaus)
{
	if(!(gSRObjects.exists(iGsnBeam) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pVoidPrecElecFldGaus == 0) return INCORRECT_PARAMS_SR_COMP_PREC;

	srTGsnBeam* pGsnBeam = dynamic_cast<srTGsnBeam*>(gSRObjects.getPtr(iGsnBeam));
    if(pGsnBeam == 0) return OBJECT_IS_NOT_GAUSBEAM;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTParPrecElecFldGaus* pPrecElecFldGaus = (srTParPrecElecFldGaus*)pVoidPrecElecFldGaus;
        srTSRWRadStructAccessData* pWfr = new srTSRWRadStructAccessData(pGsnBeam, pWfrSmp, pPrecElecFldGaus->NxNzOversamplingFactor);
		*piWfr = gSRObjects.insert(pWfr);
		//pGsnBeam->ComputeElectricFieldFreqDomain(pWfrSmp, pWfr);
		pGsnBeam->ComputeElectricField(pWfrSmp, pWfr);

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srElecFldGausBeamCompExt(void* pvWfr, int iGsnBeam, int iWfrSmp, void* pVoidPrecElecFldGaus)
{
	if(!(gSRObjects.exists(iGsnBeam) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pVoidPrecElecFldGaus == 0) return INCORRECT_PARAMS_SR_COMP_PREC;

	srTGsnBeam* pGsnBeam = dynamic_cast<srTGsnBeam*>(gSRObjects.getPtr(iGsnBeam));
    if(pGsnBeam == 0) return OBJECT_IS_NOT_GAUSBEAM;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTParPrecElecFldGaus* pPrecElecFldGaus = (srTParPrecElecFldGaus*)pVoidPrecElecFldGaus;
        srTSRWRadInData* pRadInData = (srTSRWRadInData*)pvWfr;
        srTSRWRadStructAccessData Wfr(pRadInData, pGsnBeam, pWfrSmp, pPrecElecFldGaus->NxNzOversamplingFactor);
        //pGsnBeam->ComputeElectricFieldFreqDomain(pWfrSmp, &Wfr);
        pGsnBeam->ComputeElectricField(pWfrSmp, &Wfr);
		Wfr.OutSRWRadPtrs(pRadInData);

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srStokesComp(int* piWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar)
{
	if(piWfr == 0) return SR_REFERENCE_DOES_NOT_EXIST;

	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pPrcPar == 0) return INCORRECT_PARAMS_SR_COMP_PREC;

	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;
	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTStokesStructAccessData* pStokes = new srTStokesStructAccessData(pWfrSmp);
        *piWfr = gSRObjects.insert(pStokes);
		pMagElem->ComputeSR_Stokes(pElecBeam, pWfrSmp, pPrcPar, pStokes);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srStokesCompExt(void* pVoidStokesInData, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pPrcPar == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	if(pVoidStokesInData == 0) return INCORRECT_PARAMS_SR_COMP_PREC;

	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;
	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTSRWStokesInData* pStokesInData = (srTSRWStokesInData*)pVoidStokesInData;
		srTStokesStructAccessData Stokes(pStokesInData);
		pMagElem->ComputeSR_Stokes(pElecBeam, pWfrSmp, pPrcPar, &Stokes);
		Stokes.OutDataPtrs(pStokesInData);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srPowDensComp(int* piPow, int iElecBeam, int iMagFld, int iWfrSmp, void* pvPrecPowDens)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pvPrecPowDens == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;
	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTParPrecPowDens* pPrecPowDens = (srTParPrecPowDens*)pvPrecPowDens;
        srTRadIntPowerDensity RadIntPowDens;
        srTPowDensStructAccessData* pPow = new srTPowDensStructAccessData(pWfrSmp);
        *piPow = gSRObjects.insert(pPow);
        RadIntPowDens.ComputePowerDensity(pElecBeam, pMagElem, pWfrSmp, pPrecPowDens, pPow);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srPowDensCompExt(void* pvPow, int iElecBeam, int iMagFld, int iWfrSmp, void* pvPrecPowDens)
{
	if(!(gSRObjects.exists(iElecBeam) && gSRObjects.exists(iMagFld) && gSRObjects.exists(iWfrSmp))) return SR_OBJECT_DOES_NOT_EXIST;
	if(pvPrecPowDens == 0) return INCORRECT_PARAMS_SR_COMP_PREC;
	
	srTEbmDat* pElecBeam = dynamic_cast<srTEbmDat*>(gSRObjects.getPtr(iElecBeam));
    if(pElecBeam == 0) return OBJECT_IS_NOT_EBEAM;
	srTMagElem* pMagElem = dynamic_cast<srTMagElem*>(gSRObjects.getPtr(iMagFld));
    if(pMagElem == 0) return OBJECT_IS_NOT_MAG;
	srTWfrSmp* pWfrSmp = dynamic_cast<srTWfrSmp*>(gSRObjects.getPtr(iWfrSmp));
    if(pWfrSmp == 0) return OBJECT_IS_NOT_WFR_SMP;

	try 
	{
		srTParPrecPowDens* pPrecPowDens = (srTParPrecPowDens*)pvPrecPowDens;
        srTSRWPowDensInData* pPowDensInData = (srTSRWPowDensInData*)pvPow;
		srTRadIntPowerDensity RadIntPowDens;
        srTPowDensStructAccessData Pow(pPowDensInData);
        RadIntPowDens.ComputePowerDensity(pElecBeam, pMagElem, pWfrSmp, pPrecPowDens, &Pow);

		//OC140110
		//Since srTPowDensStructAccessData Pow can be modified by ComputePowerDensity
		pPowDensInData->xStep = Pow.xStep; pPowDensInData->xStart = Pow.xStart; pPowDensInData->nx = Pow.nx;
		pPowDensInData->zStep = Pow.zStep; pPowDensInData->zStart = Pow.zStart; pPowDensInData->nz = Pow.nz;

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrSmpGet(double* pHorCenRange, int* nHor, double* pVerCenRange, int* nVer, double* pEnStartEnd, int* nEn, int iRad)
{
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrComponGet(char* pData, int iWfr, int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z)
{
	if(!gSRObjects.exists(iWfr)) return SR_OBJECT_DOES_NOT_EXIST;
	if(pData == 0) return INCORRECT_ARGUMENTS;

	CHGenObj hWfr = gSRObjects.get(iWfr);
	srTSRWRadStructAccessData* pWfr = dynamic_cast<srTSRWRadStructAccessData*>(hWfr.ptr());
    if(pWfr == 0) return OBJECT_IS_NOT_WAVEFRONT;

	try 
	{
        srTRadGenManip RadGenManip(hWfr);
		RadGenManip.ExtractRadiation(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrComponGetExt(char* pData, void* pVoidRadInData, int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z)
{
	//if(!gSRObjects.exists(iWfr)) return SR_OBJECT_DOES_NOT_EXIST;
	if(pData == 0) return INCORRECT_ARGUMENTS;
	if(pVoidRadInData == 0) return INCORRECT_WAVEFRONT_STRUCTURE;

	//CHGenObj hWfr = gSRObjects.get(iWfr);
	//srTSRWRadStructAccessData* pWfr = dynamic_cast<srTSRWRadStructAccessData*>(hWfr.ptr());
	//if(pWfr == 0) return OBJECT_IS_NOT_WAVEFRONT;

	try 
	{
		srTSRWRadInData* pRadInData = (srTSRWRadInData*)pVoidRadInData; //OC
		srTSRWRadStructAccessData Wfr(pRadInData, 0); //don't copy electric field data //OC

        //srTRadGenManip RadGenManip(hWfr);
		CHGenObj hWfr(&Wfr, true);
		srTRadGenManip RadGenManip(hWfr);

		RadGenManip.ExtractRadiation(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);

		Wfr.OutSRWRadPtrs(pRadInData); //OC

		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrComponInteg(srTDataMD* pIntensOrigData, srTDataMD* pIntegParData, srTDataMD* pIntegResData)
{//assumes all data is allocated by calling function
	if((pIntensOrigData == 0) || (pIntegParData == 0) || (pIntegResData == 0)) return 0;
	try 
	{
		srTRadGenManip::ComponInteg(pIntensOrigData, pIntegParData, pIntegResData);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}

	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrReprSet(int iWfr, char cReprID)
{
	if(!gSRObjects.exists(iWfr)) return SR_OBJECT_DOES_NOT_EXIST;
	CHGenObj hWfr = gSRObjects.get(iWfr);
	srTSRWRadStructAccessData* pWfr = dynamic_cast<srTSRWRadStructAccessData*>(hWfr.ptr());
    if(pWfr == 0) return OBJECT_IS_NOT_WAVEFRONT;

	int res = 0;
	try 
	{
		if((cReprID == 'c') || (cReprID == 'C') || (cReprID == 'a') || (cReprID == 'A'))
		{
			//if(res = pWfr->SetRepresCA(cReprID)) throw res;
			res = pWfr->SetRepresCA(cReprID);
			if(res > 0) throw res;
		}
		else if((cReprID == 'f') || (cReprID == 'F') || (cReprID == 't') || (cReprID == 'T'))
		{
			//if(res = pWfr->SetRepresFT(cReprID)) throw res;
			res = pWfr->SetRepresFT(cReprID);
			if(res > 0) throw res;
		}
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srWfrReprSetExt(void* pVoidRadInData, char cReprID)
{
	if(pVoidRadInData == 0) return INCORRECT_WAVEFRONT_STRUCTURE;

	int res = 0;
	try 
	{
		srTSRWRadInData* pRadInData = (srTSRWRadInData*)pVoidRadInData;
		srTSRWRadStructAccessData Wfr(pRadInData, 0); //don't copy electric field data

		if((cReprID == 'c') || (cReprID == 'C') || (cReprID == 'a') || (cReprID == 'A'))
		{
			//if(res = Wfr.SetRepresCA(cReprID)) throw res;
			res = Wfr.SetRepresCA(cReprID);
			if(res > 0) throw res;
		}
		else if((cReprID == 'f') || (cReprID == 'F') || (cReprID == 't') || (cReprID == 'T'))
		{
			//if(res = Wfr.SetRepresFT(cReprID)) throw res;
			res = Wfr.SetRepresFT(cReprID);
			if(res > 0) throw res;

			//OC041108: Temp. walk-around: re-calculate stat. moments after moving to frequency domain
			//because they may not be calculated correctly in Time domain.
			//To remove after fixing the srTGenOptElem::ComputeRadMoments() function!
			if((cReprID == 'f') || (cReprID == 'F')) Wfr.ComputeRadMoments();
		}
		else throw INCORRECT_WAVEFRONT_REPRES_ID;

		Wfr.OutSRWRadPtrs(pRadInData);
		UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srUtiTrfOptMirSet(double* pM, double* pV, double* pAngles, double* pP0)
{
	if((pAngles == 0) || (pM == 0) || (pV == 0)) return INCORRECT_ARGUMENTS_ARRAY;

	gmTrans T1, T2, T3; //, TrfTot, TrfTotFin;
	TVector3d PointOnAxVect(0, 0, 0);
	TVector3d AxVect1(-1, 0, 0), AxVect2(0, 0, 1), N_Orig(0, -1, 0);

	T1.SetupRotation(PointOnAxVect, AxVect1, pAngles[0]);
	T2.SetupRotation(PointOnAxVect, AxVect2, pAngles[1]);
    TrProduct(&T2, &T1, T3);

	TVector3d AxVect3 = T3.TrBiPoint(N_Orig);
	T1.SetupRotation(PointOnAxVect, AxVect3, pAngles[2]);
    TrProduct(&T1, &T3, T2);

	if(pP0 != 0)
	{
		double x0 = pP0[0], z0 = pP0[1];
		if((x0 != 0.) || (z0 != 0.))
		{
            TVector3d CenShift(x0, z0);
            T3.SetupTranslation(CenShift);
		    TrProduct(&T3, &T2, T1);
			T2 = T1;
		}
	}

	TMatrix3d OutM = T2.GetMatrix_inv();
	TVector3d OutV = T2.GetVector_inv();

	pM[0] = OutM.Str0.x; pM[1] = OutM.Str0.y; pM[2] = OutM.Str0.z;
	pM[3] = OutM.Str1.x; pM[4] = OutM.Str1.y; pM[5] = OutM.Str1.z;
	pM[6] = OutM.Str2.x; pM[7] = OutM.Str2.y; pM[8] = OutM.Str2.z;
	pV[0] = OutV.x; pV[1] = OutV.y; pV[2] = OutV.z;

	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srUtiFuncIntKnXn(double* ResInt, int type, double x1, double x2, double prec)
{
	if(prec == 0) return 0;
	try 
	{
		//*ResInt = srTMathFunctions::IntKnXn(type, x1, x2, prec);
		*ResInt = CGenMathFunc::IntKnXn(type, x1, x2, prec);
        UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srUtiRandGsnMD(double* pGsnND, double* arXc, double* arSigma, int n, char init, char rand_mode)
{
	if((pGsnND == 0) || (n == 0)) return 0;
	try 
	{
		if(init) gMathRand.Initialize();
		gMathRand.NextRandGaussND(arXc, arSigma, n, pGsnND, rand_mode);

        UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srUtiConvolComp(srTDataMD* pSrcData, srTDataMD* pDstData)
{
	if((pSrcData == 0) || (pDstData == 0)) return 0;
	try 
	{
        //srTRadGenManip RadGenManip(hWfr);
		//RadGenManip.ExtractRadiation(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);
		//UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}

	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srDel(int iObj)
{
	if(!gSRObjects.exists(iObj)) return SR_OBJECT_DOES_NOT_EXIST;
	gSRObjects.erase(iObj);
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srOptTransmCharGet(float* pData, int iOptElem, int CharType, double xc, double xr, int nx, double zc, double zr, int nz)
{
	if(pData == 0) return INCORRECT_ARGUMENTS;
	if(!gSRObjects.exists(iOptElem)) return SR_OBJECT_DOES_NOT_EXIST;
	srTGenOptElem* pOptElem = dynamic_cast<srTGenOptElem*>(gSRObjects.getPtr(iOptElem));
    if(pOptElem == 0) return OBJECT_IS_NOT_OPT_ELEM;

	try
	{
		pOptElem->ExtractTransmCharact(CharType, xc, xr, nx, zc, zr, nz, pData);
        UtiWarnCheck();
	}
	catch(int ErrNo) { return ErrNo;}
	return 0;
}

//-------------------------------------------------------------------------
