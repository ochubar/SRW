/************************************************************************//**
 * File: srigintr.cpp
 * Description: Interface to IGOR Pro (WaveMetrics)
 * Project: Synchrotron Radiation Workshop (SRW)
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srigintr.h"
#include "srigsend.h"
#include "srradind.h"
#include "srinterf.h"
#include "srappl.h"
#include "srmlttsk.h"
//#include "srsource.h" //OC031112
#include "auxparse.h"

//*************************************************************************

// All structures are 2-byte-aligned.
//#if GENERATINGPOWERPC
//#if TARGET_CPU_PPC
//	#pragma options align=mac68k
//#endif
//#ifdef _WINDOWS_
//	#pragma pack(2)
//#endif

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

//*************************************************************************

// Global Variables
srTApplication SR;
srTYield srYield;
srTIntVect gVectWarnNos;
int gCallSpinProcess = 1;

//for compatibility with lib/dll
//int (*pgCompProgressIndicFunc)(double CurVal);
//int (*pgWfrExtModifFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp);

//*************************************************************************

static int
AddCStringToHandle(						// Concatenates C string to handle.
	char *theStr,
	Handle theHand)
{
	return PtrAndHand(theStr, theHand, strlen(theStr));
}

//*************************************************************************
// Use this to trace and repair memory leaks
#ifdef __VC__
#ifdef _DEBUG

long UseOfNew = 0;
long UseOfDelete = 0;

#endif
#endif

//*************************************************************************

int CorrectErrorCode(int ErrCode)
{
	return SR.CorrectErrorCode(ErrCode);
}

//*************************************************************************

static void ProcErr(int ErrNo)
{
	if(ErrNo == 0) return;
	else if(ErrNo < 0) //warning
	{
		char WarnTextBuf[2048];
		srUtiWarnTextGet(WarnTextBuf, ErrNo);
        srTIgorSend::WarningMessage(WarnTextBuf);
	}
	else throw CorrectErrorCode(ErrNo);
}


//*************************************************************************

static void DeleteObjects(int* ObjInds, int AmOfObj)
{
	if((ObjInds == 0) || (AmOfObj == 0)) return;
	int* tInd = ObjInds;
	for(int i=0; i<AmOfObj; i++)
	{
		if(*tInd != 0) srDel(*tInd);
		tInd++;
	}
}
static void DeleteObjects(int Obj1, int Obj2)
{
	int ObjsToDel[] = {Obj1, Obj2};
	DeleteObjects(ObjsToDel, 2);
}
static void DeleteObjects(int Obj1, int Obj2, int Obj3)
{
	int ObjsToDel[] = {Obj1, Obj2, Obj3};
	DeleteObjects(ObjsToDel, 3);
}
static void DeleteObjects(int Obj1, int Obj2, int Obj3, int Obj4)
{
	int ObjsToDel[] = {Obj1, Obj2, Obj3, Obj4};
	DeleteObjects(ObjsToDel, 4);
}
static void DeleteObjects(int Obj1, int Obj2, int Obj3, int Obj4, int Obj5)
{
	int ObjsToDel[] = {Obj1, Obj2, Obj3, Obj4, Obj5};
	DeleteObjects(ObjsToDel, 5);
}

//*************************************************************************

static int srLoop(void* p_void)
{
	//SR.Send.SetIgorRadInputStruct((srTIgorRadInputStruct*)p_void);
	//return SR.ComputeSR_DirectOut();

	srTIgorRadInputStruct* pStr = (srTIgorRadInputStruct*)p_void;

	int iElecBeam = 0, iMagFld = 0, iWfrSmp = 0;
	srTSRWRadInData SRWRadInData;

	try
	{
		ProcErr(srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData));

		double I, Neb, sElecBeam, Mom1Arr[6];
		ProcErr(srTIgorSend::GetElecBeamThin(pStr->wElectronBeam, I, Neb, Mom1Arr, sElecBeam));
		ProcErr(srElecBeamSet(&iElecBeam, I, Neb, Mom1Arr, 5, 0, 0, sElecBeam));

		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wField));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));

		bool AllowAutoChoiceOfNxNzForPropag = true;
		double NxNzOversamplingParam = 1.;
		ProcErr(srTIgorSend::GetPrecParamWfrSamplingForPropag(pStr->wAuxParam, AllowAutoChoiceOfNxNzForPropag, NxNzOversamplingParam));

		int IntegMethNo = 2; 
		double RelPrecOrStep = 0.005, sStartInt = 0., sEndInt = 0.;
		bool ShowProgrIndic = true;
		//ProcErr(srTIgorSend::GetPrecParamElectricFieldComp(pStr->wIntPar, IntegMethNo, RelPrecOrStep, sStartInt, sEndInt));
		ProcErr(srTIgorSend::GetPrecParamElectricFieldComp(pStr->wIntPar, IntegMethNo, RelPrecOrStep, sStartInt, sEndInt, ShowProgrIndic));
		//srTParPrecElecFld PrecElecFld(IntegMethNo, RelPrecOrStep, sStartInt, sEndInt, NxNzOversamplingParam);
		srTParPrecElecFld PrecElecFld(IntegMethNo, RelPrecOrStep, sStartInt, sEndInt, NxNzOversamplingParam, ShowProgrIndic); //OC180312

		//int iWfr = 0;
		//ProcErr(srElecFldComp(&iWfr, iElecBeam, iMagFld, iWfrSmp, &PrecElecFld));
		ProcErr(srElecFldCompExt(&SRWRadInData, iElecBeam, iMagFld, iWfrSmp, &PrecElecFld));

		ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData);
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srWfrFromTrj(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetIgorWfrFromTrjInputStruct((srTIgorWfrFromTrjInputStruct*)p_void);
	return SR.ComputeWfrFromTrj();
}

//*************************************************************************

struct srTIgorUtiWfrCSR {
	waveHndl wRad;
	waveHndl wIntPar;
//Integration Parameters wave:
// 0: Meth. No (=0 for simplest method)
// 1: Prec. Param. (~1)
// 2: Long. Pos. to Start integ.
// 3: Long. Pos. to Finish integ.
// 4: Use auto-sampling for propag. (=0 means don't use)
// 5: Over-sampling param.
	waveHndl wObservation;
	waveHndl wField;
	waveHndl wElectronBeam;
};
static int srWfrCSR(void* p_void)
{
	srTIgorUtiWfrCSR* pStr = (srTIgorUtiWfrCSR*)p_void;

	int iElecBeam = 0, iMagFld = 0, iWfrSmp = 0;
	srTSRWRadInData SRWRadInData;

	try
	{
		ProcErr(srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData));

		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wField));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));
		ProcErr(srTIgorSend::GetAndSetElecBeamGen(&iElecBeam, pStr->wElectronBeam));

		double arPrecParam[8];
		double *pPrecParam = arPrecParam;
		long auxNp = 0;
        ProcErr(srTIgorSend::GetArrDoubleFromNumWave1D(pStr->wIntPar, 8, pPrecParam, auxNp));

		ProcErr(csrElecFldCompExt(&SRWRadInData, iElecBeam, iMagFld, iWfrSmp, arPrecParam));

		ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData);
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srWfrIsotrSrc(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetIgorIsotrSrcInputStruct((srTIgorIsotrSrcInputStruct*)p_void);
	return SR.ComputeWfrIsotrSrc();
}

//*************************************************************************

static int srWfrGsnBeam(void* p_void)
{
	//SR.Send.SetIgorGsnBeamInputStruct((srTIgorGsnBeamInputStruct*)p_void);
	//return SR.ComputeWfrGsnBeam();

	//int res = 0;
	srTIgorGsnBeamInputStruct* pStr = (srTIgorGsnBeamInputStruct*)p_void;

	int iGsnBeam = 0, iWfrSmp = 0;
	srTSRWRadInData SRWRadInData;

	try
	{
		ProcErr(srTIgorSend::GetAndSetGaussianBeam(&iGsnBeam, pStr->wGsnBeam));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));
		ProcErr(srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData));

		bool AllowAutoChoiceOfNxNzForPropag = true;
		double NxNzOversamplingParam = -1.;
		ProcErr(srTIgorSend::GetPrecParamWfrSamplingForPropag(pStr->wAuxParam, AllowAutoChoiceOfNxNzForPropag, NxNzOversamplingParam));

		srTParPrecElecFldGaus PrecElecFldGaus(NxNzOversamplingParam);

		//int iWfr = 0;
		//ProcErr(srElecFldGausBeamComp(&iWfr, iGsnBeam, iWfrSmp, &PrecElecFldGaus));
		ProcErr(srElecFldGausBeamCompExt(&SRWRadInData, iGsnBeam, iWfrSmp, &PrecElecFldGaus));

		ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
		DeleteObjects(iWfrSmp, iGsnBeam);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData);
		DeleteObjects(iWfrSmp, iGsnBeam);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srWfrSASE(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetIgorWfrSASEInputStruct((srTIgorWfrSASEInputStruct*)p_void);
	return SR.ComputeWfrSASE();
}

//*************************************************************************

static int srTraj(void* p_void)
{
	//SR.Send.SetIgorTrjInputStruct((srTIgorTrjInputStruct*)p_void);
	//return SR.ComputeTrajectory();

	int res = 0;
	srTIgorTrjInputStruct* pStr = (srTIgorTrjInputStruct*)p_void;

	double *pOutBtxData=0, *pOutXData=0, *pOutBtyData=0, *pOutYData=0, *pOutBtzData=0, *pOutZData=0;
	int hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor;
	
	int iElecBeam = 0, iMagFld = 0;
	try
	{
        int ns;
		double sStart, sStep;
		ProcErr(srTIgorSend::GetTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, pOutBtxData, pOutXData, pOutBtzData, pOutZData, ns, sStart, sStep, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor));
		ProcErr(srTIgorSend::GetAndSetElecBeamThin(&iElecBeam, pStr->wElectronBeam));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wField));

		ProcErr(srElecTrjComp(pOutBtxData, pOutXData, pOutBtyData, pOutYData, pOutBtzData, pOutZData, ns, sStart, sStep, iElecBeam, iMagFld));

		DeleteObjects(iElecBeam, iMagFld);
		ProcErr(srTIgorSend::FinishWorkingWithTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor));
	}
	catch(int ErrNo)
	{
		if(res = srTIgorSend::FinishWorkingWithTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor)) return res;
		DeleteObjects(iElecBeam, iMagFld);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct srTIgorTrj3dInputStruct {
	waveHndl wOutVerCoor; // Vertical trajectory coordinate [m] vs s [m]
	waveHndl wOutVerAng; // Vertical trajectory angle [r] vs s [m]
	waveHndl wOutLongCoor; // Longitudinal trajectory coordinate [m] vs s [m]
	waveHndl wOutLongDer; // Horizontal trajectory angle [r] vs s [m]
	waveHndl wOutHorCoor; // Horizontal trajectory coordinate [m] vs s [m]
	waveHndl wOutHorAng; // Horizontal trajectory angle [r] vs s [m]

	waveHndl wField; // Magnetic Field structure (text wave)
	waveHndl wElectronBeam; // Electron Beam wave
};
static int srTraj3d(void* p_void)
{
	int res = 0;
	srTIgorTrj3dInputStruct* pStr = (srTIgorTrj3dInputStruct*)p_void;

	double *pOutBtxData=0, *pOutXData=0, *pOutBtyData=0, *pOutYData=0, *pOutBtzData=0, *pOutZData=0;
	int hStateHorAng, hStateHorCoor, hStateLongDer, hStateLongCoor, hStateVerAng, hStateVerCoor;
	
	int iElecBeam = 0, iMagFld = 0;
	try
	{
        int ns;
		double sStart, sStep;
		ProcErr(srTIgorSend::GetTrj3dDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutLongDer, pStr->wOutLongCoor, pStr->wOutVerAng, pStr->wOutVerCoor, pOutBtxData, pOutXData, pOutBtyData, pOutYData, pOutBtzData, pOutZData, ns, sStart, sStep, hStateHorAng, hStateHorCoor, hStateLongDer, hStateLongCoor, hStateVerAng, hStateVerCoor));
		ProcErr(srTIgorSend::GetAndSetElecBeamThin(&iElecBeam, pStr->wElectronBeam));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wField));

		ProcErr(srElecTrjComp(pOutBtxData, pOutXData, pOutBtyData, pOutYData, pOutBtzData, pOutZData, ns, sStart, sStep, iElecBeam, iMagFld));

		DeleteObjects(iElecBeam, iMagFld);
        ProcErr(srTIgorSend::FinishWorkingWithTrj3dDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutLongDer, pStr->wOutLongCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateLongDer, hStateLongCoor, hStateVerAng, hStateVerCoor));
	}
	catch(int ErrNo)
	{
		if(res = srTIgorSend::FinishWorkingWithTrj3dDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutLongDer, pStr->wOutLongCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateLongDer, hStateLongCoor, hStateVerAng, hStateVerCoor)) return res;
		DeleteObjects(iElecBeam, iMagFld);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srElecBeamPropag(void* p_void)
{
	//int res = 0;
	srTIgorElecBeamPropag* pStr = (srTIgorElecBeamPropag*)p_void;

	int iElecBeam = 0, iOptElem = 0;
	try
	{
		double I, Neb, s0_ebm, Mom1Arr[6], Mom2Arr[30];
		double HorTwissArr[4], VertTwissArr[4], HorEmit, VertEmit;
		int nMom2;
		int TypeDistrTrans, TypeDistrLong;
		double NoiseFactor;
		ProcErr(srTIgorSend::GetElecBeamThick(pStr->wElectronBeam, I, Neb, Mom1Arr, Mom2Arr, nMom2, s0_ebm, TypeDistrTrans, TypeDistrLong, NoiseFactor));
		ProcErr(srTIgorSend::GetElecBeamTwiss(pStr->wElectronBeam, HorTwissArr, VertTwissArr));
		ProcErr(srElecBeamSet(&iElecBeam, I, Neb, Mom1Arr, 5, Mom2Arr, nMom2, s0_ebm));

		double &SigmaX2 = Mom2Arr[0], &MixMomX = Mom2Arr[1], &SigmaPrimeX2 = Mom2Arr[2];
		double &SigmaZ2 = Mom2Arr[3], &MixMomZ = Mom2Arr[4], &SigmaPrimeZ2 = Mom2Arr[5];
		double &SigmaE2 = Mom2Arr[10];
        ProcErr(srElecBeamMom2EmitAndTwiss(&HorEmit, HorTwissArr, HorTwissArr + 1, SigmaX2, SigmaPrimeX2, MixMomX, SigmaE2, HorTwissArr[2], HorTwissArr[3]));
        ProcErr(srElecBeamMom2EmitAndTwiss(&VertEmit, VertTwissArr, VertTwissArr + 1, SigmaZ2, SigmaPrimeZ2, MixMomZ, SigmaE2, VertTwissArr[2], VertTwissArr[3]));

		vector<string> OptElemInfo;
		ProcErr(srTIgorSend::GetVectorOfStrings(pStr->wOptElem, &OptElemInfo));
		
		int AmOfStrings = OptElemInfo.size();
		if(AmOfStrings == 0) return 0;
        char** DescrArr = new char*[AmOfStrings];
		for(int i=0; i<AmOfStrings; i++) 
		{
			const char* CurStr = OptElemInfo[i].c_str();
			if(CurStr == 0) continue;
			int CurLen = strlen(CurStr);
			if(CurLen == 0) continue;
			DescrArr[i] = new char[CurLen + 1];
			strcpy(DescrArr[i], CurStr);
		}
		ProcErr(srOptElemSet(&iOptElem, DescrArr, &AmOfStrings, 0, 0));

		ProcErr(srElecBeamPropagate(iElecBeam, iOptElem));
        ProcErr(srElecBeamGet(&I, Mom1Arr, Mom2Arr, &s0_ebm, iElecBeam));

		//BUG fixed: propagation in presence of dispersion (??!!): Eta chnages, EtaPrime stays constant
        //ProcErr(srElecBeamMom2EmitAndTwiss(&HorEmit, HorTwissArr, HorTwissArr + 1, SigmaX2, SigmaPrimeX2, MixMomX, SigmaE2, HorTwissArr[2], HorTwissArr[3]));
        //ProcErr(srElecBeamMom2EmitAndTwiss(&VertEmit, VertTwissArr, VertTwissArr + 1, SigmaZ2, SigmaPrimeZ2, MixMomZ, SigmaE2, VertTwissArr[2], VertTwissArr[3]));

		ProcErr(srElecBeamMomAndEmit2Twiss(HorTwissArr, HorTwissArr + 1, HorTwissArr + 2, &HorEmit, SigmaX2, SigmaPrimeX2, MixMomX, SigmaE2, HorTwissArr[3])); //EtaPrime is input only
		ProcErr(srElecBeamMomAndEmit2Twiss(VertTwissArr, VertTwissArr + 1, VertTwissArr + 2, &VertEmit, SigmaZ2, SigmaPrimeZ2, MixMomZ, SigmaE2, VertTwissArr[3])); //EtaPrime is input only

		ProcErr(srTIgorSend::SetElecBeamThick(pStr->wElectronBeam, I, Mom1Arr, Mom2Arr, s0_ebm, TypeDistrTrans, TypeDistrLong, NoiseFactor));
		ProcErr(srTIgorSend::SetElecBeamEmitAndTwiss(pStr->wElectronBeam, HorEmit, VertEmit, HorTwissArr, VertTwissArr));

		DeleteObjects(iElecBeam, iOptElem);
        
		for(int k=0; k<AmOfStrings; k++)
		{
			if(DescrArr[k] != 0) delete[] (DescrArr[k]);
			DescrArr[k] = 0;
		}
		if(DescrArr != 0) delete[] DescrArr;

	}
	catch(int ErrNo)
	{
		DeleteObjects(iElecBeam, iOptElem);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srVerNo(void* p_void)
{
	char cVersionStr[200];
    srUtiVerNo(cVersionStr);

	long LenVersionStr = strlen(cVersionStr);
    Handle VersionStr = NewHandle(LenVersionStr);
	strncpy(*VersionStr, cVersionStr, LenVersionStr);
	((srTIgorVersionStruct*)p_void)->result = VersionStr;
	return 0;
}

//*************************************************************************

static int srRadMom(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetHandleOfSRWRadStruct((srTHandleOfSRWRadStruct*)p_void);
	return SR.ComputeMomentsOfSR();
}

//*************************************************************************

static int srFFT1D(void* p_void)
{//to modify: use DLL API !!!
	return SR.Make1DFFT((srTHandleOfOneWaveStruct*)p_void);
}

//*************************************************************************

static int srFFT2D(void* p_void)
{//to modify: use DLL API !!!
	return SR.Make2DFFT((srTHandleOfOneWaveStruct*)p_void);
}

//*************************************************************************

static int srRadResizeXZ(void* p_void)
{//to modify: use DLL API !!!
	return SR.RadResize((srTIgorRadResizeInputStruct*)p_void);
}

//*************************************************************************

static int srObsSetNxNz(void* p_void)
{//to modify: use DLL API !!!
	return SR.ObsSetNxNz((srTIgorObsSetNxNzInputStruct*)p_void);
}

//*************************************************************************

static int srRadPropag(void* p_void)
{//to modify: use DLL API !!!
	return SR.PropagateRadiation((srTIgorRadPropagInputStruct*)p_void);
}

//*************************************************************************

static int srWfrPropag(void* p_void)
{//to modify: use DLL API !!!
	return SR.PropagateWavefront((srTIgorWfrPropag*)p_void);
}

//*************************************************************************

//static int srWfrChangeRep(void* p_void)
//{//to modify: use DLL API !!!
// //obsolette version
//	srTIgorWfrChangeRepInputStruct *pInStr = (srTIgorWfrChangeRepInputStruct*)p_void;
//	srTHandleOfSRWRadStruct AuxHandleOfSRWRadStruct;
//	AuxHandleOfSRWRadStruct.wRad = (pInStr->wRad);
//	SR.Send.SetHandleOfSRWRadStruct(&AuxHandleOfSRWRadStruct);
//
//	return SR.WfrChangeRep((int)(pInStr->MethNo));
//}

//*************************************************************************

struct srTIgorWfrSetRepresInputStruct {
	Handle sCoordOrTime;
	waveHndl wRad;
};
static int srWfrSetRepres(void* p_void)
{
	//srTIgorWfrChangeRepInputStruct *pInStr = (srTIgorWfrChangeRepInputStruct*)p_void;
	//srTHandleOfSRWRadStruct AuxHandleOfSRWRadStruct;
	//AuxHandleOfSRWRadStruct.wRad = (pInStr->wRad);
	//SR.Send.SetHandleOfSRWRadStruct(&AuxHandleOfSRWRadStruct);
	//return SR.WfrChangeRep((int)(pInStr->MethNo));

	int res = 0;
	srTIgorWfrSetRepresInputStruct* p = (srTIgorWfrSetRepresInputStruct*)p_void;

	//int iWfr = 0;
	srTSRWRadInData SRWRadInData;
	if(res = srTIgorSend::GetSRWRadInData(p->wRad, &SRWRadInData)) return res;

	try
	{
		//ProcErr(srWfrSet(&iWfr, &SRWRadInData, 0)); //do not copy radiation data
		ProcErr(srWfrReprSetExt(&SRWRadInData, **(p->sCoordOrTime)));
		//ProcErr(srDel(iWfr));
        ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
	}
	catch(int ErrNo)
	{
        //srDel(iWfr);
		if(res = srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData)) return res;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srRadExtract(void* p_void)
{
	//return SR.ExtractRadiation((srTIgorRadExtractInputStruct*)p_void);

	int res = 0;
	srTIgorRadExtractInputStruct* pStr = (srTIgorRadExtractInputStruct*)p_void;

			//DEBUG
			//srTIgorSend::WarningMessage("srRadExtract: just in");
			//END DEBUG

	int iWfr = 0;
	srTSRWRadInData SRWRadInData;
	if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	try
	{
		int PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
		int Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E)
		int PlotType; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
		int TransvPres; // 0: Spatial; 1: Angular
		//int FreqTimePres; // 0: Frequency; 1: Time
		double e, x, z;

			//DEBUG
			//srTIgorSend::WarningMessage("srRadExtract: before GetIntensExtractParam");
			//END DEBUG

		ProcErr(srTIgorSend::GetIntensExtractParam(pStr->wExtractParam, PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z));
		char* pExtrData = 0;
		int hStateExtractedData;
		ProcErr(srTIgorSend::GetIntensExtractData(pStr->wExtractedData, hStateExtractedData, pExtrData));
		//srTParIntensExtract ParIntensExtract(PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z); //OC031208

		srTIgorWaveAccessData ExtrWaveData;
		//ProcErr(srWfrSet(&iWfr, &SRWRadInData, 1));

		//ProcErr(srWfrSet(&iWfr, &SRWRadInData, 0)); //do not copy radiation data //OC060211

			//DEBUG
			//srTIgorSend::WarningMessage("srRadExtract: before srWfrComponGet");
			//END DEBUG

		//ProcErr(srWfrComponGet(pExtrData, iWfr, PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z));
		ProcErr(srWfrComponGetExt(pExtrData, &SRWRadInData, PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z)); //OC031208

			//DEBUG
			//srTIgorSend::WarningMessage("srRadExtract: after srWfrComponGet");
			//END DEBUG

		//ProcErr(srTIgorSend::SetupExtractedWaveData(&SRWRadInData, Int_or_Phase, PlotType, TransvPres, pExtrData, pStr->wExtractedData, hStateExtractedData, &ExtrWaveData));
		ProcErr(srTIgorSend::SetupExtractedWaveData(&SRWRadInData, Int_or_Phase, PlotType, pExtrData, pStr->wExtractedData, hStateExtractedData, &ExtrWaveData));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&ExtrWaveData));
		//ProcErr(srDel(iWfr)); //OC031208
        ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
	}
	catch(int ErrNo)
	{
        srDel(iWfr);
		if(res = srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData)) return res;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srStokesUnd(void* p_void)
{
	//SR.Send.SetIgorPerStokesInputStruct((srTIgorPerStokesInputStruct*)p_void);
	//return SR.ComputePerStokes();

	//int res = 0;
	srTIgorPerStokesInputStruct* pStr = (srTIgorPerStokesInputStruct*)p_void;
	int iElecBeam=0, iMagFld=0, iWfrSmp=0;
	srTSRWStokesInData StokesInData;

	try
	{
		ProcErr(srTIgorSend::GetAndSetElecBeamThick(&iElecBeam, pStr->wElectronBeam));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wPerField));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));

		int InitHarm, FinHarm;
		double Kns, Knphi;
		char IntensityOrFlux;
		double minPhotEnExtRight = 1; //OC170713
        //ProcErr(srTIgorSend::GetPrecParamStokesPerComp(pStr->wPerIntPar, InitHarm, FinHarm, Kns, Knphi, IntensityOrFlux));
        ProcErr(srTIgorSend::GetPrecParamStokesPerComp(pStr->wPerIntPar, InitHarm, FinHarm, Kns, Knphi, IntensityOrFlux, minPhotEnExtRight)); //OC170713
		//srTParPrecStokesPer PrecStokesPer(InitHarm, FinHarm, Kns, Knphi, IntensityOrFlux);
		srTParPrecStokesPer PrecStokesPer(InitHarm, FinHarm, Kns, Knphi, IntensityOrFlux, minPhotEnExtRight); //OC170713

		ProcErr(srTIgorSend::GetSRWStokesInData(pStr->wStokes, &StokesInData));

		//int iWfr = 0;
		//ProcErr(srStokesComp(&iWfr, iElecBeam, iMagFld, iWfrSmp, &PrecStokesPer));
		ProcErr(srStokesCompExt(&StokesInData, iElecBeam, iMagFld, iWfrSmp, &PrecStokesPer));

		ProcErr(srTIgorSend::FinishWorkingWithSRWStokesStruct(&StokesInData));
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWStokesStruct(&StokesInData);
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srMagArb2Per(void* p_void)
{
	//int res = 0;
	srTIgorMagArb2Per* pStr = (srTIgorMagArb2Per*)p_void;
	int iMagArb=0, iMagPer=0;

	try
	{
		//double FieldAbsZeroTol = 1.E-12, sStartB, sStepB, *pBH, *pBV;
		//long NpB;
		//bool BHIsZero, BVIsZero;
		//ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wTrUnifMagField, FieldAbsZeroTol, sStartB, sStepB, NpB, pBH, BHIsZero, pBV, BVIsZero));
		//ProcErr(srMagFldTrUnifSet(&iMagArb, sStartB, sStepB, NpB, pBH, pBV));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagArb, pStr->wTrUnifMagField));

		double RelPrec, MaxPerLen_m;
		int MaxAmOfHarm;
		ProcErr(srTIgorSend::GetPrecParamMagArb2Per(pStr->wPrecPar, RelPrec, MaxAmOfHarm, MaxPerLen_m));

		ProcErr(srMagFldPerSetFromTrUnif(&iMagPer, iMagArb, RelPrec, MaxAmOfHarm, MaxPerLen_m));

		double PerLength, TotLength, sCen, SpecPar, TaperPar_TU=0, PhaseSh_OK=0, Fund_keV_per_GeV2;
		int AmOfHarm, TypeOfUnd; // 0- infinite; 1- normal; 2- tapered; 3- optical clystron
		char Type; //'c'- conventional, 't'- tapered, 'k'- optical klystron, 'i'- infinite

		ProcErr(srMagFldPerGet(iMagPer, &PerLength, &TotLength, &sCen, &AmOfHarm, &Type, &SpecPar, &Fund_keV_per_GeV2));

		if(Type == 'i') TypeOfUnd = 0;
		else if(Type == 'c') TypeOfUnd = 1;
		else if(Type == 't') { TypeOfUnd = 2; TaperPar_TU = SpecPar;}
		else if(Type == 'k') { TypeOfUnd = 3; PhaseSh_OK = SpecPar;}

		int* ArrHarmNo=0;
		char* ArrXorZ=0;
		double *ArrK=0, *ArrPhase=0;
		if(AmOfHarm > 0)
		{
            ArrHarmNo = new int[AmOfHarm];
            ArrXorZ = new char[AmOfHarm];
            ArrK = new double[AmOfHarm];
			ArrPhase = new double[AmOfHarm];
            ProcErr(srMagFldPerHarmGet(iMagPer, ArrHarmNo, ArrXorZ, ArrK, ArrPhase));
        }

		if(AmOfHarm > 20) //limitation for Igor Pro
		{
            AmOfHarm = 20;
			ProcErr(TOO_LARGE_MAX_MAG_FIELD_HARMONIC_NUMBER);
		}

        ProcErr(srTIgorSend::SetMagFieldPeriodic(pStr->wPerMagField, PerLength, TotLength, AmOfHarm, ArrHarmNo, ArrXorZ, ArrK, ArrPhase, TypeOfUnd, TaperPar_TU, PhaseSh_OK, Fund_keV_per_GeV2));

		if(ArrHarmNo != 0) delete[] ArrHarmNo;
		if(ArrXorZ != 0) delete[] ArrXorZ;
		if(ArrK != 0) delete[] ArrK;
		if(ArrPhase != 0) delete[] ArrPhase;

		DeleteObjects(iMagArb, iMagPer);
	}
	catch(int ErrNo)
	{
		DeleteObjects(iMagArb, iMagPer);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srStokesWig(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetIgorStokesWigInputStruct((srTIgorStokesWigInputStruct*)p_void);
	return SR.ComputeStokesWig();
}

//*************************************************************************

static int srStokesConst(void* p_void)
{//to modify: use DLL API !!!
	SR.Send.SetIgorStokesConstInputStruct((srTIgorStokesConstInputStruct*)p_void);
	return SR.ComputeStokesConst();
}

//*************************************************************************

static int srPowDens(void* p_void)
{
	//SR.Send.SetIgorPowDensInputStruct((srTIgorPowDensInputStruct*)p_void);
	//return SR.ComputePowerDensity();

	int res = 0;
	srTIgorPowDensInputStruct* pStr = (srTIgorPowDensInputStruct*)p_void;
	
	int iElecBeam=0, iMagFld=0, iWfrSmp=0;
	srTSRWPowDensInData SRWPowDensInData;
    if(res = srTIgorSend::GetSRWPowDensInData(pStr->wPowDens, &SRWPowDensInData)) return res;

	try
	{
		ProcErr(srTIgorSend::GetAndSetElecBeamThick(&iElecBeam, pStr->wElectronBeam));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wGenMagField));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));

		double PrecFact;
		int Method, UseSpecIntLim;
		double sIntStart, sIntFin;
		ProcErr(srTIgorSend::GetPrecParamPowDensComp(pStr->wPowDensIntPar, PrecFact, Method, UseSpecIntLim, sIntStart, sIntFin));
		srTParPrecPowDens PrecPowDens(Method, PrecFact, UseSpecIntLim, sIntStart, sIntFin);

		//int iWfr = 0;
		//ProcErr(srPowDensComp(&iWfr, iElecBeam, iMagFld, iWfrSmp, &PrecPowDens));
		ProcErr(srPowDensCompExt(&SRWPowDensInData, iElecBeam, iMagFld, iWfrSmp, &PrecPowDens));

		ProcErr(srTIgorSend::FinishWorkingWithSRWPowDensStruct(&SRWPowDensInData));
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWPowDensStruct(&SRWPowDensInData);
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

void srEnsureSrwInit()
{
	char TestName[] = "SrwVerProc";
	char StrDummy[50];
	int SRW_WasNotInit = FetchStrVar(TestName, StrDummy);
	if(SRW_WasNotInit)
	{
		char InitStr[] = "SrwInit(1)";
		XOPCommand(InitStr);
	}
	SetXOPType(RESIDENT);
}

//*************************************************************************

//static int srOptThinGenSetup(void* p_void)
//{//to modify: use DLL API !!!
//	return SR.OptGenTransSetup((srTIgorOptGenTransSetupInputStruct*)p_void);
//}

struct TIgor_srOptElemSetup {
	waveHndl wOptElem; //first argument
};
static int srOptElemSetup(void* p_void) //used for ThinGen, ThickMirGen
{//replaced the int srTApplication::OptGenTransSetup(srTIgorOptGenTransSetupInputStruct* pInputStruct)

	TIgor_srOptElemSetup *p = (TIgor_srOptElemSetup*)p_void;
	int iOptElem = 0, iWfrDummy = 0;
	try
	{
		ProcErr(srTIgorSend::GetAndSetOptElem(&iOptElem, &iWfrDummy, p->wOptElem, NULL));
		//the above makes setup and update of "p->wOptElem" eventually

		srDel(iOptElem);
	}
	catch(int ErrNo)
	{
		srDel(iOptElem);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srOptMatConst(void* p_void)
{//to modify: use DLL API !!!
	return SR.OptMatConst((srTIgorOptMatConstInputStruct*)p_void);
}

//*************************************************************************

static int srOptZonePlateSpecSetup(void* p_void)
{//to modify: use DLL API !!!
	return SR.OptZonePlateSpecSetup((srTIgorOptZonePlateSpecSetupInputStruct*)p_void);
}

//*************************************************************************

static int srUtiInterTime(void* p_void)
{
	srYield.Init(((srTIgorUtiInterTimeInputStruct*)p_void)->Delta);
	return 0;
}

//*************************************************************************

static int srUtiSpotInfo(void* p_void)
{//to modify: use DLL API !!!
//Calculates the following params of the Spot:
//[0]: Max. value //*(pOut++) = (float)MaxVal;
//[1]: Hor. pos. of max. //*(pOut++) = (float)xMax;
//[2]: Vert. pos. of max. //*(pOut++) = (float)yMax;
//[3]: Hor. FWHM size //*(pOut++) = (float)(xRight - xLeft);
//[4]: Vert. FWHM size //*(pOut++) = (float)(yRight - yLeft);
	return SR.UtiSpotInfo((srTIgorUtiSpotInfoInputStruct*)p_void);
}

//*************************************************************************

static int srUtiWfrLimits(void* p_void)
{//to modify: use DLL API !!!
	return SR.UtiWfrLimits((srTIgorUtiWrfLimitsInputStruct*)p_void);
}

//*************************************************************************

static int srUtiRemoveFlips(void* p_void)
{
	return SR.UtiRemoveFlips((srTIgorUtiRemoveFlipsInputStruct*)p_void);
}

//*************************************************************************

static int srUtiMagRad(void* p_void)
{//to modify: use DLL API !!!
	srTIgorUtiMagRadInputStruct* p = (srTIgorUtiMagRadInputStruct*)p_void;
	double Rmag = SR.UtiMagRad(p->Bconst, p->ElecEnergy, p->OutUnit);
	((srTIgorUtiMagRadInputStruct*)p_void)->Rmag = Rmag;
	return 0;
}

//*************************************************************************

static int srUtiMagCritPhotEn(void* p_void)
{//to modify: use DLL API !!!
	srTIgorUtiMagCritPhotEnInputStruct* p = (srTIgorUtiMagCritPhotEnInputStruct*)p_void;
	double CritPhotEn = SR.UtiMagCritPhotEn(p->Bconst, p->ElecEnergy, p->OutUnit);
	((srTIgorUtiMagCritPhotEnInputStruct*)p_void)->CritPhotEn = CritPhotEn;
	return 0;
}

//*************************************************************************

static int srUtiUndK(void* p_void)
{//to modify: use DLL API !!!
	srTIgorUtiUndKInputStruct* p = (srTIgorUtiUndKInputStruct*)p_void;
	double K = SR.UtiUndK(p->Bpeak, p->Period);
	((srTIgorUtiUndKInputStruct*)p_void)->K = K;
	return 0;
}

//*************************************************************************

static int srUtiUndFundPhotEn(void* p_void)
{//to modify: use DLL API !!!
	srTIgorUtiUndFundPhotEnInputStruct* p = (srTIgorUtiUndFundPhotEnInputStruct*)p_void;

	//int aha = 1;
	double FundPhotEn = SR.UtiUndFundPhotEn(p->Bpeak, p->Period, p->ElecEnergy, p->OutUnit);
	((srTIgorUtiUndFundPhotEnInputStruct*)p_void)->FundPhotEn = FundPhotEn;
	return 0;
}

//*************************************************************************

struct srTIgorUtiRandGsnInputStruct {
    waveHndl wGsn;
	DOUBLE meth;
	DOUBLE init;
	DOUBLE n;
    waveHndl wXcSig;

	DOUBLE res; //first gsn number
};
static int srUtiRandGsn(void* p_void)
{
	double *pGsn = 0, *pXcSig = 0;
	try
	{
		srTIgorUtiRandGsnInputStruct* p = (srTIgorUtiRandGsnInputStruct*)p_void;

		p->res = 0;
		if(p->n <= 0) return 0;

        long nGsn = 0;
		ProcErr(srTIgorSend::GetArrDoubleFromNumWave1D(p->wGsn, -1, pGsn, nGsn));
		if((pGsn == 0) || (nGsn < p->n)) 
		{
			ProcErr(srTIgorSend::ReDimNumWave1D(p->wGsn, (int)(p->n)));
			ProcErr(srTIgorSend::GetArrDoubleFromNumWave1D(p->wGsn, -1, pGsn, nGsn));
			nGsn = (int)(p->n);
		}

        int NumDims = 0, NumCols = 0;
		ProcErr(srTIgorSend::GetArrDoubleFromNumWave2D(p->wXcSig, NumDims, NumCols, pXcSig));
		if((pXcSig == 0) || (NumDims <= 0) || (NumCols < 2)) throw INCORRECT_GAUSS_CENTERS_SIGMAS;
		double *pXc = pXcSig, *pSigma = pXcSig + NumDims;
		ProcErr(srUtiRandGsnMD(pGsn, pXc, pSigma, (int)(p->n), (char)(p->init), (char)(p->meth)));

		ProcErr(srTIgorSend::SetDataInNumWave(p->wGsn, pGsn, nGsn));
		p->res = *pGsn;

		if(pGsn != 0) delete[] pGsn;
		if(pXcSig != 0) delete[] pXcSig;
	}
	catch(int ErrNo)
	{
		if(pGsn != 0) delete[] pGsn;
		if(pXcSig != 0) delete[] pXcSig;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srKn(void* p_void)
{//to modify: use DLL API !!!
	return SR.Kmu((srTIgorKmuInputStruct*)p_void);
}

//*************************************************************************

struct srTIgorUtiIntKnXnInputStruct {
	DOUBLE prec;
	DOUBLE x2;
	DOUBLE x1;
	DOUBLE type;

	DOUBLE res;
};
static int srUtiIntKnXn(void* p_void)
{
	try
	{
		double ResInt;
        ProcErr(srUtiFuncIntKnXn(&ResInt, (int)(((srTIgorUtiIntKnXnInputStruct*)p_void)->type), (double)(((srTIgorUtiIntKnXnInputStruct*)p_void)->x1), (double)(((srTIgorUtiIntKnXnInputStruct*)p_void)->x2), (double)(((srTIgorUtiIntKnXnInputStruct*)p_void)->prec)));
		((srTIgorUtiIntKnXnInputStruct*)p_void)->res = ResInt;
	}
	catch(int ErrNo)
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srStokesArb(void* p_void)
{
	srTIgorStokesArbInputStruct* pStr = (srTIgorStokesArbInputStruct*)p_void;

	int iElecBeam = 0, iMagFld = 0, iWfrSmp = 0;
	srTParPrecStokesArb PrecStokesArb;
	srTSRWStokesInData StokesInData;

	try
	{
		ProcErr(srTIgorSend::GetAndSetElecBeamThick(&iElecBeam, pStr->wElectronBeam));
		ProcErr(srTIgorSend::GetAndSetMagFieldGen(&iMagFld, pStr->wGenMagField));
		ProcErr(srTIgorSend::GetAndSetWfrSampling(&iWfrSmp, pStr->wObservation));
		ProcErr(srTIgorSend::GetPrecParamStokesArbComp(pStr->wPrecPar, &PrecStokesArb));
		ProcErr(srTIgorSend::GetSRWStokesInData(pStr->wStokes, &StokesInData));

		ProcErr(srStokesCompExt(&StokesInData, iElecBeam, iMagFld, iWfrSmp, &PrecStokesArb));

		ProcErr(srTIgorSend::FinishWorkingWithSRWStokesStruct(&StokesInData));
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWStokesStruct(&StokesInData);
		DeleteObjects(iElecBeam, iMagFld, iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct srTIgorUtiTrfOptMir {
    waveHndl wVect;
    waveHndl wMatr;
	waveHndl wPointAndAng;
	DOUBLE d; //return
};
static int srUtiTrfOptMir(void* p_void)
{
	try
	{
        srTIgorUtiTrfOptMir *p = (srTIgorUtiTrfOptMir*)p_void;

        double *pPointAndAng=0;
        long nPointAndAng = 0;
		ProcErr(srTIgorSend::GetArrDoubleFromNumWave1D(p->wPointAndAng, -1, pPointAndAng, nPointAndAng));
		if((nPointAndAng < 5) || (pPointAndAng == 0)) return 0;

		double Angles[] = {pPointAndAng[0], pPointAndAng[1], pPointAndAng[2]};
		double P0[] = {pPointAndAng[3], pPointAndAng[4]};

		double M[9], V[3];
		ProcErr(srUtiTrfOptMirSet(M, V, Angles, P0));
		p->d = 0;

        ProcErr(srTIgorSend::ReDimNumWave2D(p->wMatr, 3, 3));
        ProcErr(srTIgorSend::ReDimNumWave2D(p->wVect, 3, 0));

        ProcErr(srTIgorSend::SetDataInNumWave(p->wMatr, M, 9));
        ProcErr(srTIgorSend::SetDataInNumWave(p->wVect, V, 3));
	}
	catch(int ErrNo) 
	{
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

struct srTIgorRadIntensInteg {
    waveHndl wIntegRes;
	waveHndl wIntegPar;
	waveHndl wIntensOrig;
};
static int srRadIntensInteg(void* p_void)
{
	srTIgorRadIntensInteg *p = (srTIgorRadIntensInteg*)p_void;
	//int hStateIntensOrig, hStateIntegPar, hStateIntegRes;
	srTDataMD IntensOrigData, IntegParData, IntegResData;
	try
	{
		ProcErr(srTIgorSend::GetNumWaveData(p->wIntensOrig, &IntensOrigData)); //, hStateIntensOrig));
		ProcErr(srTIgorSend::GetNumWaveData(p->wIntegPar, &IntegParData)); //, hStateIntegPar));
		ProcErr(srTIgorSend::GetNumWaveData(p->wIntegRes, &IntegResData)); //, hStateIntegRes));

		ProcErr(srWfrComponInteg(&IntensOrigData, &IntegParData, &IntegResData));

		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntensOrigData, p->wIntensOrig)); //, hStateIntensOrig));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntegParData, p->wIntegPar)); //, hStateIntegPar));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntegResData, p->wIntegRes)); //, hStateIntegRes));
	}
	catch(int ErrNo) 
	{
		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntensOrigData, p->wIntensOrig)); //, hStateIntensOrig));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntegParData, p->wIntegPar)); //, hStateIntegPar));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&IntegResData, p->wIntegRes)); //, hStateIntegRes));
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srElecEnergyModul(void* p_void)
{
	int res = 0;
//ddddddddddddddddddddddddddddddddd
	//srTIgorElecEnergyModul* pStr = (srTIgorElecEnergyModul*)p_void;

	//int iElecBeam = 0, iMagFld = 0, iWfrSmp = 0;
	//srTSRWRadInData SRWRadInData;
	//if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	return res;
}

//*************************************************************************
//TEST
static int srUti3DView(void* p_void)
{//to allow 2D drawing and viewing beamline feometries

	//SR.Send.SetIgor3DViewInputStruct((srTIgor3DViewInputStruct*)p_void);
	int result = 0;
	
	CHGenObj hSrc;
	//if(result = SR.Send.GetRadSourceData(p_void, "srTIgor3DViewInputStruct", hSrc)) return result;
	char nmStr[] = "srTIgor3DViewInputStruct\0";
	if(result = SR.Send.GetRadSourceData(p_void, nmStr, hSrc)) return result;

	srTEbmDat ElBeam;
	//if(result = SR.Send.GetTotalElectronBeamDataFormat3(ElBeam)) return result;
	srTMagElem MagField;


//ddddddddddddddddddddddddddddddddd
	return SR.Uti3DView(ElBeam, MagField);
}

//*************************************************************************
//TEST
static int srRadPropagStokesMultiE(void* p_void)
{
	srTIgorRadPropagStokesMultiElecInputStruct* pRadPropagInputStruct = (srTIgorRadPropagStokesMultiElecInputStruct*)p_void;
	SR.Send.SetIgorRadPropagStokesMultiElecInputStruct(pRadPropagInputStruct);
	return SR.PropagateRadiationStokesMultiElec(pRadPropagInputStruct);
}

//*************************************************************************
//TEST
static int srWfrEmitPropagStokesMultiE(void* p_void)
{
	srTIgorRadEmitPropagStokesMultiElecInputStruct* pRadPropagInputStruct = (srTIgorRadEmitPropagStokesMultiElecInputStruct*)p_void;
	SR.Send.SetIgorRadEmitPropagStokesMultiElecInputStruct(pRadPropagInputStruct);
	return SR.EmitPropagRadStokesMultiElec(pRadPropagInputStruct);
}

//*************************************************************************
//TEST
static int srWfrEmitPropag(void* p_void)
{
	SR.Send.SetIgorWfrEmitPropagInputStruct((srTIgorWfrEmitPropagInputStruct*)p_void);
	return SR.ComputeWfrEmitPropag();
}

//*************************************************************************
//TEST
static int srWfrReflect(void* p_void)
{
	SR.Send.SetIgorWfrReflectInputStruct((srTIgorWfrReflectInputStruct*)p_void);
	return SR.ReflectWavefront();
}

//*************************************************************************

struct TIgor_srOptThinTransmCharExtract {
	waveHndl wToExtract; //last argument
	waveHndl wExtractPar;
	waveHndl wOptElem; //first argument
};
static int srOptThinTransmCharExtract(void* p_void)
{
	TIgor_srOptThinTransmCharExtract *p = (TIgor_srOptThinTransmCharExtract*)p_void;
	//int hStateToExtract;
    srTDataMD mdToExtract;
	int iOptElem = 0, iWfrAux = 0;
	try
	{
		ProcErr(srTIgorSend::GetNumWaveData(p->wToExtract, &mdToExtract)); //, hStateToExtract));
		if(mdToExtract.DataType[0] != 'f') throw NT_FP32_WAVE_REQUIRED;

		double ParData[10];
		double *pParData = ParData;
		long LenParData = 7;
        ProcErr(srTIgorSend::GetArrDoubleFromNumWave1D(p->wExtractPar, 10, pParData, LenParData));

		ProcErr(srTIgorSend::GetAndSetOptElem(&iOptElem, &iWfrAux, p->wOptElem, NULL));

		//double xr = mdToExtract.DimSteps[0]*(mdToExtract.DimSizes[0] - 1);
		//double zr = mdToExtract.DimSteps[1]*(mdToExtract.DimSizes[1] - 1);
		int CharType = (int)ParData[0];
		double xc = ParData[1], xr = ParData[2], zc = ParData[4], zr = ParData[5];
		int nx = (int)ParData[3], nz = (int)ParData[6];
		if((nx > 0) && (nz == 0)) nz = 1;
		if((nx == 0) && (nz > 0)) nx = 1;

		if(nx <= 1) xr = 0;
		if(nz <= 1) zr = 0;

		ProcErr(srOptTransmCharGet((float*)mdToExtract.pData, iOptElem, CharType, xc, xr, nx, zc, zr, nz));

		ProcErr(srTIgorSend::FinishWorkingWithWave(&mdToExtract, p->wToExtract)); //, hStateToExtract));
		srDel(iOptElem);
	}
	catch(int ErrNo)
	{
        ProcErr(srTIgorSend::FinishWorkingWithWave(&mdToExtract, p->wToExtract)); //, hStateToExtract));
		srDel(iOptElem);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************
/*	RegisterFunction()
	Igor calls this at startup time to find the address of the
	XFUNCs added by this XOP. See XOP manual regarding "Direct XFUNCs".
*/
static long RegisterFunction()
{
	int funcIndex;

	funcIndex = GetXOPItem(0);		// Which function is Igor asking about?
	switch (funcIndex) {
		case 0:
			return((long)srLoop);
			break;
		case 1:
			return((long)srWfrFromTrj);
			break;
		case 2:
			return((long)srWfrIsotrSrc);
			break;
		case 3:
			return((long)srWfrGsnBeam);
			break;
		case 4:
			return((long)srWfrSASE);
			break;
		case 5:
			return((long)srTraj);
			break;
		case 6:
			return((long)srVerNo);
			break;
		case 7:
			return((long)srRadMom);
			break;
		case 8:
			return((long)srFFT1D);
			break;
		case 9:
			return((long)srFFT2D);
			break;
		case 10:
			return((long)srRadResizeXZ);
			break;
		case 11:
			return((long)srObsSetNxNz);
			break;
		case 12:
			return((long)srRadPropag);
			break;
		case 13:
			return((long)srWfrPropag);
			break;
		case 14:
			return((long)srRadExtract);
			break;
		case 15:
			return((long)srStokesUnd);
			break;
		case 16:
			return((long)srStokesWig);
			break;
		case 17:
			return((long)srStokesConst);
			break;
		case 18:
			return((long)srPowDens);
			break;
		case 19:
			return((long)srUtiInterTime);
			break;
		case 20:
			//return((long)srOptThinGenSetup);
			return((long)srOptElemSetup);
			break;
		case 21:
			return((long)srOptMatConst);
			break;
		case 22:
			return((long)srOptZonePlateSpecSetup);
			break;
		case 23:
			return((long)srUtiSpotInfo);
			break;
		case 24:
			return((long)srUtiWfrLimits);
			break;
		case 25:
			return((long)srUtiRemoveFlips);
			break;
		case 26:
			return((long)srUtiMagRad);
			break;
		case 27:
			return((long)srUtiMagCritPhotEn);
			break;
		case 28:
			return((long)srUtiUndK);
			break;
		case 29:
			return((long)srUtiUndFundPhotEn);
			break;
		case 30:
			return((long)srKn);
			break;
		case 31:
			return((long)srElecBeamPropag);
			break;
		case 32:
			return((long)srMagArb2Per);
			break;
		case 33:
			return((long)srStokesArb);
			break;
		case 34:
			return((long)srUtiTrfOptMir);
			break;
		case 35:
			return((long)srUti3DView);
			break;
		case 36:
			return((long)srRadPropagStokesMultiE);
			break;
		//case 37:
		//	return((long)srWfrEmitPropagStokesMultiE);
		//	break;
		case 37:
			return((long)srWfrEmitPropag);
			break;
		case 38:
			return((long)srWfrReflect);
			break;
		case 39:
			return((long)srWfrSetRepres);
			break;
		case 40:
			return((long)srUtiIntKnXn);
			break;
		case 41:
			return((long)srRadIntensInteg);
			break;
		case 42:
			return((long)srOptThinTransmCharExtract);
			break;
		case 43:
			return((long)srWfrCSR);
			break;
		case 44:
			return((long)srTraj3d);
			break;
		case 45:
			return((long)srUtiRandGsn);
			break;
		//case 45:
		//	return((long)srOptThickMirGenSetup);
		//	break;
	}
	return(NIL);
}

//*************************************************************************
/*	DoFunction()
	
	Igor calls this when the user invokes one if the XOP's XFUNCs
	if we returned NIL for the XFUNC from RegisterFunction. In this
	XOP, we always use the direct XFUNC method, so Igor will never call
	this function. See XOP manual regarding "Direct XFUNCs".
*/
static int DoFunction()
{
	int funcIndex;
	void *p;				// Pointer to structure containing function parameters and result.
	int err;

	funcIndex = GetXOPItem(0);		// Which function is being invoked ?
	p = (void *)GetXOPItem(1);		// Get pointer to params/result.
	switch (funcIndex) {
		case 0:
			err = srLoop(p);
			break;
		case 1:
			err = srWfrFromTrj(p);
			break;
		case 2:
			err = srWfrIsotrSrc(p);
			break;
		case 3:
			err = srWfrGsnBeam(p);
			break;
		case 4:
			err = srWfrSASE(p);
			break;
		case 5:
			err = srTraj(p);
			break;
		case 6:
			err = srVerNo(p);
			break;
		case 7:
			err = srRadMom(p);
			break;
		case 8:
			err = srFFT1D(p);
			break;
		case 9:
			err = srFFT2D(p);
			break;
		case 10:
			err = srRadResizeXZ(p);
			break;
		case 11:
			err = srObsSetNxNz(p);
			break;
		case 12:
			err = srRadPropag(p);
			break;
		case 13:
			err = srWfrPropag(p);
			break;
		case 14:
			err = srRadExtract(p);
			break;
		case 15:
			err = srStokesUnd(p);
			break;
		case 16:
			err = srStokesWig(p);
			break;
		case 17:
			err = srStokesConst(p);
			break;
		case 18:
			err = srPowDens(p);
			break;
		case 19:
			err = srUtiInterTime(p);
			break;
		case 20:
			//err = srOptThinGenSetup(p);
			err = srOptElemSetup(p);
			break;
		case 21:
			err = srOptMatConst(p);
			break;
		case 22:
			err = srOptZonePlateSpecSetup(p);
			break;
		case 23:
			err = srUtiSpotInfo(p);
			break;
		case 24:
			err = srUtiWfrLimits(p);
			break;
		case 25:
			err = srUtiRemoveFlips(p);
			break;
		case 26:
			err = srUtiMagRad(p);
			break;
		case 27:
			err = srUtiMagCritPhotEn(p);
			break;
		case 28:
			err = srUtiUndK(p);
			break;
		case 29:
			err = srUtiUndFundPhotEn(p);
			break;
		case 30:
			err = srKn(p);
			break;
		case 31:
			err = srElecBeamPropag(p);
			break;
		case 32:
			err = srMagArb2Per(p);
			break;
		case 33:
			err = srStokesArb(p);
			break;
		case 34:
			err = srUtiTrfOptMir(p);
			break;
		case 35:
			err = srUti3DView(p);
			break;
		case 36:
			err = srRadPropagStokesMultiE(p);
			break;
		//case 37:
		//	err = srWfrEmitPropagStokesMultiE(p);
		//	break;
		case 37:
			err = srWfrEmitPropag(p);
			break;
		case 38:
			err = srWfrReflect(p);
			break;
		case 39:
			err = srWfrSetRepres(p);
			break;
		case 40:
			err = srUtiIntKnXn(p);
			break;
		case 41:
			err = srRadIntensInteg(p);
			break;
		case 42:
			err = srOptThinTransmCharExtract(p);
			break;
		case 43:
			err = srWfrCSR(p);
			break;
		case 44:
			err = srTraj3d(p);
			break;
		case 45:
			err = srUtiRandGsn(p);
			break;
		//case 45:
		//	err = srOptThickMirGenSetup(p);
		//	break;
	}
	return(err);
}

//*************************************************************************
/*	XOPEntry()

	This is the entry point from the host application to the XOP for all messages after the
	INIT message.
*/
static void
XOPEntry(void)
{	
	long result = 0;

	switch (GetXOPMessage()) {
		case FUNCTION:						// Our external function being invoked ?
			result = DoFunction();
			break;

		case FUNCADDRS:
			result = RegisterFunction();
			break;

		case IDLE:
			srEnsureSrwInit();
			break;
	}
	SetXOPResult(result);
}

//*************************************************************************
/*	main(ioRecHandle)

	This is the initial entry point at which the host application calls XOP.
	The message sent by the host must be INIT.
	main() does any necessary initialization and then sets the XOPEntry field of the
	ioRecHandle to the address to be called for future messages.
*/
//HOST_IMPORT void main(IORecHandle ioRecHandle)
HOST_IMPORT int main(IORecHandle ioRecHandle) //OC030110, required for GCC 4.2 Mac OSX
{
	//#ifdef applec					// For MPW C for 68K only.
	//	void _DATAINIT(void);
	//	_DATAINIT();				// For MPW C only.
	//	UnloadSeg(_DATAINIT);
	//#endif
	
	//#ifdef XOP_GLOBALS_ARE_A4_BASED
	//	#ifdef __MWERKS__
	//		SetCurrentA4();							// Set up correct A4. This allows globals to work.
	//		SendXOPA4ToIgor(ioRecHandle, GetA4());	// And communicate it to Igor.
	//	#endif
	//#endif

	//LoadXOPSegs();
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	if(igorVersion < 200) SetXOPResult(OLD_IGOR);
	else SetXOPResult(0L);
	
	SetXOPType(RESIDENT | IDLES);
	SR.Initialize();
	srYield.Init(0.5);

	srUtiProgrIndFuncSet(&(srTIgorSend::ShowCompProgress));
	srUtiWfrExtModifFuncSet(&(srTIgorSend::WfrModify));
	srUtiOptElemGetInfByNameFuncSet(&(srTIgorSend::GetOptElemInfByName));

	return 0;
}

//*************************************************************************

// All structures are 2-byte-aligned.
//#if GENERATINGPOWERPC
//#if TARGET_CPU_PPC
//	#pragma options align=reset
//#endif
//#ifdef _WINDOWS_
//	#pragma pack()
//#endif

#include "XOPStructureAlignmentReset.h"	// All structures passed between Igor and XOP are two-byte aligned.

//*************************************************************************
