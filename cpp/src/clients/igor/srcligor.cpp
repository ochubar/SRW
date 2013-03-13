
#ifndef __SRIGINTR_H
#include "srigintr.h"
#endif
#ifndef __SRIGSEND_H
#include "srigsend.h"
#endif

#ifndef __SRRADIND_H
#include "srradind.h"
#endif
#ifndef __SRINTERF_H
#include "srinterf.h"
#endif

//#ifndef __STLSTART_H
//#include "stlstart.h"
//#endif
#include <vector>

// All structures are 2-byte-aligned.
#if GENERATINGPOWERPC
	#pragma options align=mac68k
#endif
#ifdef _WINDOWS_
	#pragma pack(2)
#endif

//*************************************************************************

// Global Variables
//int gCallSpinProcess = 1;
//static char gErrWarnBuffer[2048];

//*************************************************************************

static int
AddCStringToHandle(						// Concatenates C string to handle.
	char *theStr,
	Handle theHand)
{
	return PtrAndHand(theStr, theHand, strlen(theStr));
}

//*************************************************************************

static int CorrectErrorCode(int ErrCode)
{
	ErrCode += (FIRST_XOP_ERR - FIRST_ERR_SRWDLL);
	if(ErrCode == 1) return SR_COMP_PROC_ABORTED; // 1 is returned when user presses Abort
	return ErrCode;
}

//*************************************************************************

static void ProcErr(int ErrNo)
{
	if(ErrNo == 0) return;
	else if(ErrNo < 0) //warning
	{
		char WarnTextBuf[2048];
		srGetWarnText(WarnTextBuf, ErrNo);
        srTIgorSend::WarningMessage(WarnTextBuf);
	}
	else throw CorrectErrorCode(ErrNo);
}

//*************************************************************************

static int srVerNo(void* p_void)
{
	char cVersionStr[] = "3.9";
	long LenVersionStr = 10;

	Handle VersionStr = NewHandle(LenVersionStr);
	strncpy(*VersionStr, cVersionStr, LenVersionStr);

	((srTIgorVersionStruct*)p_void)->result = VersionStr;
	return 0;
}

//*************************************************************************

static int srLoop(void* p_void)
{
	int res = 0;
	srTIgorRadInputStruct* pStr = (srTIgorRadInputStruct*)p_void;

	int iElecBeam = 0, iMagFld = 0, iWfrSmp = 0;
	srTSRWRadInData SRWRadInData;
	//if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	try
	{
		ProcErr(srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData));

		double I, sElecBeam, Mom1Arr[6];
		ProcErr(srTIgorSend::GetElecBeamThin(pStr->wElectronBeam, I, Mom1Arr, sElecBeam));

		//double sStartBH, sStepBH, sStartBV, sStepBV;
		double sStartB, sStepB;
		double *pBH = 0, *pBV = 0;
		//long NpBH = 0, NpBV = 0;
		long NpB = 0;
		double FieldAbsZeroTol = 1.E-12; //[T]
		bool BHIsZero = true, BVIsZero = true;
		//ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wField, FieldAbsZeroTol, sStartBH, sStepBH, NpBH, pBH, BHIsZero, sStartBV, sStepBV, NpBV, pBV, BVIsZero));
		ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wField, FieldAbsZeroTol, sStartB, sStepB, NpB, pBH, BHIsZero, pBV, BVIsZero));

		double yObs, zSt, zFi, xSt, xFi, eSt, eFi;
		int nz, nx, ne;
		char PhotEnUnits[10];
		ProcErr(srTIgorSend::GetWfrSampling(pStr->wObservation, yObs, zSt, zFi, nz, xSt, xFi, nx, eSt, eFi, ne, PhotEnUnits));

		bool AllowAutoChoiceOfNxNzForPropag = true;
		double NxNzOversamplingParam = 1.;
		ProcErr(srTIgorSend::GetPrecParamWfrSamplingForPropag(pStr->wAuxParam, AllowAutoChoiceOfNxNzForPropag, NxNzOversamplingParam));

		int IntegMethNo = 2; 
		double RelPrecOrStep = 0.005, sStartInt = 0., sEndInt = 0.;
		ProcErr(srTIgorSend::GetPrecParamElectricFieldComp(pStr->wIntPar, IntegMethNo, RelPrecOrStep, sStartInt, sEndInt));
		srTParPrecElecFld PrecElecFld(0, IntegMethNo, RelPrecOrStep, sStartInt, sEndInt, NxNzOversamplingParam);

		ProcErr(srSetElecBeam(&iElecBeam, I, Mom1Arr, 5, 0, 0, sElecBeam));
		ProcErr(srSetMagFldTrUnif(&iMagFld, sStartB, sStepB, NpB, pBH, pBV));
		ProcErr(srSetWfrSmp(&iWfrSmp, yObs, xSt, xFi, nx, zSt, zFi, nz, eSt, eFi, ne, PhotEnUnits));

		//int iWfr = 0;
		//ProcErr(srCompElecFld(&iWfr, iElecBeam, iMagFld, iWfrSmp, &PrecElecFld));
		ProcErr(srCompElecFld(&SRWRadInData, iElecBeam, iMagFld, iWfrSmp, &PrecElecFld));

		ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
		ProcErr(srDel(iElecBeam));
		ProcErr(srDel(iMagFld));
		ProcErr(srDel(iWfrSmp));
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData);
		if(iElecBeam != 0) srDel(iElecBeam);
		if(iMagFld != 0) srDel(iMagFld);
		if(iWfrSmp != 0) srDel(iWfrSmp);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srWfrFromTrj(void* p_void)
{
//	SR.Send.SetIgorWfrFromTrjInputStruct((srTIgorWfrFromTrjInputStruct*)p_void);
//	return SR.ComputeWfrFromTrj();
	return 0;
}

//*************************************************************************

static int srWfrIsotrSrc(void* p_void)
{
//	SR.Send.SetIgorIsotrSrcInputStruct((srTIgorIsotrSrcInputStruct*)p_void);
//	return SR.ComputeWfrIsotrSrc();
	return 0;
}

//*************************************************************************

static int srWfrGsnBeam(void* p_void)
{
//	SR.Send.SetIgorGsnBeamInputStruct((srTIgorGsnBeamInputStruct*)p_void);
//	return SR.ComputeWfrGsnBeam();
	return 0;
}

//*************************************************************************

static int srWfrSASE(void* p_void)
{
//	SR.Send.SetIgorWfrSASEInputStruct((srTIgorWfrSASEInputStruct*)p_void);
//	return SR.ComputeWfrSASE();
	return 0;
}

//*************************************************************************

static int srTraj(void* p_void)
{
	int res = 0;
	srTIgorTrjInputStruct* pStr = (srTIgorTrjInputStruct*)p_void;

	double *pOutBtxData=0, *pOutXData=0, *pOutBtzData=0, *pOutZData=0;
	int hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor;
	
	int iElecBeam = 0, iMagFld = 0;
	try
	{
		ProcErr(srTIgorSend::GetTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, pOutBtxData, pOutXData, pOutBtzData, pOutZData, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor));

		double I, sElecBeam, Mom1Arr[6];
		ProcErr(srTIgorSend::GetElecBeamThin(pStr->wElectronBeam, I, Mom1Arr, sElecBeam));

		//double sStartBH, sStepBH, sStartBV, sStepBV;
		double sStartB, sStepB;
		double *pBH = 0, *pBV = 0;
		//long NpBH = 0, NpBV = 0;
		long NpB = 0;
		double FieldAbsZeroTol = 1.E-12; //[T]
		bool BHIsZero = true, BVIsZero = true;
		//ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wField, FieldAbsZeroTol, sStartBH, sStepBH, NpBH, pBH, BHIsZero, sStartBV, sStepBV, NpBV, pBV, BVIsZero));
		ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wField, FieldAbsZeroTol, sStartB, sStepB, NpB, pBH, BHIsZero, pBV, BVIsZero));

		ProcErr(srSetElecBeam(&iElecBeam, I, Mom1Arr, 5, 0, 0, sElecBeam));
		ProcErr(srSetMagFldTrUnif(&iMagFld, sStartB, sStepB, NpB, pBH, pBV));

		ProcErr(srCompTrj(pOutBtxData, pOutXData, pOutBtzData, pOutZData, iElecBeam, iMagFld));

		ProcErr(srDel(iElecBeam));
		ProcErr(srDel(iMagFld));
		ProcErr(srTIgorSend::FinishWorkingWithTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor));
	}
	catch(int ErrNo)
	{
		if(res = srTIgorSend::FinishWorkingWithTrjDataPointers(pStr->wOutHorAng, pStr->wOutHorCoor, pStr->wOutVerAng, pStr->wOutVerCoor, hStateHorAng, hStateHorCoor, hStateVerAng, hStateVerCoor)) return res;
		if(iElecBeam != 0) srDel(iElecBeam);
		if(iMagFld != 0) srDel(iMagFld);
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srRadMom(void* p_void)
{
	//SR.Send.SetHandleOfSRWRadStruct((srTHandleOfSRWRadStruct*)p_void);
	//return SR.ComputeMomentsOfSR();
	return 0;
}

//*************************************************************************

static int srFFT1D(void* p_void)
{
	//return SR.Make1DFFT((srTHandleOfOneWaveStruct*)p_void);
	return 0;
}

//*************************************************************************

static int srFFT2D(void* p_void)
{
	//return SR.Make2DFFT((srTHandleOfOneWaveStruct*)p_void);
	return 0;
}

//*************************************************************************

static int srRadResizeXZ(void* p_void)
{
	//return SR.RadResize((srTIgorRadResizeInputStruct*)p_void);
	return 0;
}

//*************************************************************************

static int srObsSetNxNz(void* p_void)
{
	//return SR.ObsSetNxNz((srTIgorObsSetNxNzInputStruct*)p_void);
	return 0;
}

//************************************************************************

static int SetupOpticalElement(int* piOptElem, vector<string>* pOptElemInfo, srTSRWRadInData* pRadInData)
{
	const char* ElemID = ((*pOptElemInfo)[0]).c_str();

	if(strcmp(ElemID, "Drift") == 0)
	{
		double Length = atof(((*pOptElemInfo)[1]).c_str());
		ProcErr(srSetOptDrift(piOptElem, Length));
	}
	else if(strcmp(ElemID, "Container") == 0)
	{
		//int AmOfElem = atoi((*pOptElemInfo)[1]);
		//OptElemHndl = srTGenOptElemHndl(new srTCompositeOptElem(pOptElemInfo, AmOfElem, pRad));
	}

/*

	else if(!strcmp(ElemID, "RectAperture"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTRectAperture(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "CircAperture"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTCircAperture(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "ThinLens"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTThinLens(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "SpherMirror"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTSpherMirror(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "WaveguideRect"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTWaveguideRect(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "ZonePlate"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTZonePlate(pOptElemInfo));
	}

	else if(!strcmp(ElemID, "ThinGen"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTGenTransmission(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "PhaseShift"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTPhaseShift(pOptElemInfo, pRad));
	}
	else if(!strcmp(ElemID, "ZonePlateSpec"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTZonePlateSpec(pOptElemInfo, pRad));
	}

// Continue for new Optical Elements

	else return UNKNOWN_OPTICAL_ELEMENT;

	if(OptElemHndl.rep->ErrorCode != 0) return OptElemHndl.rep->ErrorCode;

*/
	return 0;
}

//*************************************************************************

static int SetupOpticalElementContainer(int* piOptElem, vector<string>* pElemInfo, srTSRWRadInData* pRadInData)
{
	int res = 0;
	int AmOfMembers = pElemInfo->size() - 1;
	int* piContElemArr = 0;
	if(AmOfMembers > 0) piContElemArr = new int[AmOfMembers];

	for(int i=1; i<=AmOfMembers; i++)
	{
		const char* MemberID = ((*pElemInfo)[i]).c_str();
		vector<string> MemberInfo;

		ProcErr(srTIgorSend::GetVectorOfStrings(MemberID, &MemberInfo));
		ProcErr(SetupOpticalElement((piContElemArr + (i - 1)), &MemberInfo, pRadInData));
	}
	ProcErr(srSetOptCont(piOptElem, piContElemArr, AmOfMembers));
	if(piContElemArr != 0) delete[] piContElemArr;
	return 0;
}

//*************************************************************************

static int srRadPropag(void* p_void)
{
	int res = 0;
	srTIgorRadPropagInputStruct* pStr = (srTIgorRadPropagInputStruct*)p_void;

	srTSRWRadInData SRWRadInData;
	if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	int iOptElem = 0;

	try
	{
		vector<string> OptElemInfo;
		ProcErr(srTIgorSend::GetVectorOfStrings(pStr->wOptElem, &OptElemInfo));
		ProcErr(SetupOpticalElement(&iOptElem, &OptElemInfo, &SRWRadInData));



		ProcErr(srDel(iOptElem));
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData);
		if(iOptElem != 0) srDel(iOptElem);
		return ErrNo;
	}

	//dddddddddddddddddddddddddddddd
/*
	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTStringVect OptElemInfo;
	int result;
	if(result = Send.GetSRWRadStructAndOptElemNames(pRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return result;

	srTGenOptElemHndl OptElemHndl;
	srTOptElemSummary OptElemSummary;
	result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	if(result != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return result;}

	if(result = PerformMethodDependentPrePropagationActions((int)(pRadPropagInputStruct->MethNo), (int)(pRadPropagInputStruct->AuxPar1), OptElemHndl, SRWRadStructAccessData)) return result;
	if(result = OptElemHndl.rep->CheckRadStructForPropagation(&SRWRadStructAccessData)) return result;

	srTRadResizeVect AuxResizeVect;
	result = OptElemHndl.rep->PropagateRadiation(&SRWRadStructAccessData, (int)(pRadPropagInputStruct->MethNo), AuxResizeVect);

//Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);
	return CorrectErrorCode(result);
*/
	//return SR.PropagateRadiation((srTIgorRadPropagInputStruct*)p_void);
	return 0;
}

//*************************************************************************
//
//static int srWfrChangeRep(void* p_void)
//{
//	srTIgorWfrChangeRepInputStruct *pInStr = (srTIgorWfrChangeRepInputStruct*)p_void;
//	srTHandleOfSRWRadStruct AuxHandleOfSRWRadStruct;
//	AuxHandleOfSRWRadStruct.wRad = (pInStr->wRad);
//	SR.Send.SetHandleOfSRWRadStruct(&AuxHandleOfSRWRadStruct);
//
//	return SR.WfrChangeRep((int)(pInStr->MethNo));
//}
//
//*************************************************************************

static int srRadExtract(void* p_void)
{
	int res = 0;
	srTIgorRadExtractInputStruct* pStr = (srTIgorRadExtractInputStruct*)p_void;

	srTSRWRadInData SRWRadInData;
	if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	try
	{
		int PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
		int Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E)
		int PlotType; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
		int TransvPres; // 0: Spatial; 1: Angular
		double e, x, z;
		ProcErr(srTIgorSend::GetIntensExtractParam(pStr->wExtractParam, PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z));
		char* pExtrData = 0;
		int hStateExtractedData;
		ProcErr(srTIgorSend::GetIntensExtractData(pStr->wExtractedData, hStateExtractedData, pExtrData));
		srTParIntensExtract ParIntensExtract(PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z);

		int iWfr = 0;
		ProcErr(srSetWfr(&iWfr, &SRWRadInData, 1));
		ProcErr(srGetWfrCompon(pExtrData, iWfr, PolarizCompon, Int_or_Phase, PlotType, TransvPres, e, x, z));

		srTIgorWaveAccessData ExtrWaveData;
		ProcErr(srTIgorSend::SetupExtractedWaveData(&SRWRadInData, Int_or_Phase, PlotType, TransvPres, pExtrData, pStr->wExtractedData, hStateExtractedData, &ExtrWaveData));
		ProcErr(srTIgorSend::FinishWorkingWithWave(&ExtrWaveData));
		ProcErr(srDel(iWfr));
        ProcErr(srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData));
	}
	catch(int ErrNo)
	{
		if(res = srTIgorSend::FinishWorkingWithSRWRadStruct(&SRWRadInData)) return res;
		return ErrNo;
	}
	return 0;
}

//*************************************************************************

static int srStokesUnd(void* p_void)
{
//	SR.Send.SetIgorPerStokesInputStruct((srTIgorPerStokesInputStruct*)p_void);
//	return SR.ComputePerStokes();
	return 0;
}

//*************************************************************************

static int srStokesWig(void* p_void)
{
//	SR.Send.SetIgorStokesWigInputStruct((srTIgorStokesWigInputStruct*)p_void);
//	return SR.ComputeStokesWig();
	return 0;
}

//*************************************************************************

static int srStokesConst(void* p_void)
{
	//SR.Send.SetIgorStokesConstInputStruct((srTIgorStokesConstInputStruct*)p_void);
	//return SR.ComputeStokesConst();
	return 0;
}

//*************************************************************************

static int SetMagFldPerHarmonics(int*& iHarmArr, int AmOfHarm, srTMagFldHarm* MagFldHarmArr)
{
	if((AmOfHarm <= 0) || (MagFldHarmArr == 0)) return 0;

	iHarmArr = new int[AmOfHarm];
	int* tiHarm = iHarmArr;
	srTMagFldHarm* tMagFldHarm = MagFldHarmArr;

	for(int i=0; i<AmOfHarm; i++)
	{
		char v_or_h = ((tMagFldHarm->XorZ == 'x') || (tMagFldHarm->XorZ == 'X'))? 'h' : 'v';
		ProcErr(srSetMagFldPerHarm(tiHarm, tMagFldHarm->n, v_or_h, tMagFldHarm->K, tMagFldHarm->Phase));
		tMagFldHarm++;
		tiHarm++;
	}
	return 0;
}

//*************************************************************************

static int srPowDens(void* p_void)
{
	int res = 0;
	srTIgorPowDensInputStruct* pStr = (srTIgorPowDensInputStruct*)p_void;
	
	int iElecBeam=0, iMagFld=0, iWfrSmp=0;
	srTSRWPowDensInData SRWPowDensInData;
    if(res = srTIgorSend::GetSRWPowDensInData(pStr->wPowDens, &SRWPowDensInData)) return res;

	try
	{
		double I, s0_ebm, Mom1Arr[6], Mom2Arr[30];
		int nMom2;
		int TypeDistrTrans, TypeDistrLong;
		double NoiseFactor;
		ProcErr(srTIgorSend::GetElecBeamThick(pStr->wElectronBeam, I, Mom1Arr, Mom2Arr, nMom2, s0_ebm, TypeDistrTrans, TypeDistrLong, NoiseFactor));
		ProcErr(srSetElecBeam(&iElecBeam, I, Mom1Arr, 5, Mom2Arr, nMom2, s0_ebm));

        int MagFieldType = 0;
		ProcErr(srTIgorSend::IdentifyMagFieldType(pStr->wGenMagField, MagFieldType));
		if(MagFieldType == 1) 
		{
            double FieldAbsZeroTol = 1.E-12, sStartB, sStepB, *pBH, *pBV;
            long NpB;
            bool BHIsZero, BVIsZero;
			ProcErr(srTIgorSend::GetMagFieldTransvUnif(pStr->wGenMagField, FieldAbsZeroTol, sStartB, sStepB, NpB, pBH, BHIsZero, pBV, BVIsZero));
            ProcErr(srSetMagFldTrUnif(&iMagFld, sStartB, sStepB, NpB, pBH, pBV));
		}
		else if(MagFieldType == 2) 
		{
			double PerLength, TotLength, TaperPar_TU, PhaseSh_OK, FldErrRMS, NatFocNxSASE, NatFocNySASE, TaperStartSASE, TaperRelFldChgSASE;
			int AmOfHarm, TypeOfUnd, FldErrTypeSASE, TaperTypeSASE;
			// 0- infinite; 1- normal; 2- tapered; 3- optical clystron
			srTMagFldHarm* MagFldHarmArr = 0;
            ProcErr(srTIgorSend::GetMagFieldPeriodic(pStr->wGenMagField, PerLength, TotLength, AmOfHarm, MagFldHarmArr, TypeOfUnd, TaperPar_TU, PhaseSh_OK, FldErrTypeSASE, FldErrRMS, NatFocNxSASE, NatFocNySASE, TaperTypeSASE, TaperStartSASE, TaperRelFldChgSASE));
			int* iHarmArr = 0;
			ProcErr(SetMagFldPerHarmonics(iHarmArr, AmOfHarm, MagFldHarmArr));
			double SpecPar = 0.;
			char Type = 'c'; //'c'- conventional, 't'- tapered, 'k'- optical klystron, 'i'- infinite
			if(TypeOfUnd == 0) Type = 'i';
			else if(TypeOfUnd == 1) Type = 'c';
			else if(TypeOfUnd == 2) { Type = 't'; SpecPar = TaperPar_TU;}
			else if(TypeOfUnd == 3) { Type = 'k'; SpecPar = PhaseSh_OK;}
			ProcErr(srSetMagFldPer(&iMagFld, PerLength, TotLength, 0., iHarmArr, AmOfHarm, Type, SpecPar));
			if(MagFldHarmArr != 0) delete[] MagFldHarmArr;
		}
		else if(MagFieldType == 3) 
		{
			double BH, BV;
            ProcErr(srTIgorSend::GetMagFieldConstant(pStr->wGenMagField, BH, BV));
			ProcErr(srSetMagFldConst(&iMagFld, BH, BV));
		}

		double yObs, zSt, zFi, xSt, xFi, eSt, eFi;
		int nz, nx, ne;
		char PhotEnUnits[10];
		ProcErr(srTIgorSend::GetWfrSampling(pStr->wObservation, yObs, zSt, zFi, nz, xSt, xFi, nx, eSt, eFi, ne, PhotEnUnits));
		ProcErr(srSetWfrSmp(&iWfrSmp, yObs, xSt, xFi, nx, zSt, zFi, nz, eSt, eFi, ne, PhotEnUnits));

		double PrecFact;
		int Method;
		ProcErr(srTIgorSend::GetPrecParamPowDensComp(pStr->wPowDensIntPar, PrecFact, Method));
		srTParPrecPowDens PrecPowDens(0, Method, PrecFact);

		//int iWfr = 0;
		//ProcErr(srCompPowDens(&iWfr, iElecBeam, iMagFld, iWfrSmp, &PrecPowDens));
		ProcErr(srCompPowDens(&SRWPowDensInData, iElecBeam, iMagFld, iWfrSmp, &PrecPowDens));

		ProcErr(srTIgorSend::FinishWorkingWithSRWPowDensStruct(&SRWPowDensInData));
		ProcErr(srDel(iElecBeam));
		ProcErr(srDel(iMagFld));
		ProcErr(srDel(iWfrSmp));
	}
	catch(int ErrNo)
	{
		srTIgorSend::FinishWorkingWithSRWPowDensStruct(&SRWPowDensInData);
		if(iElecBeam != 0) srDel(iElecBeam);
		if(iMagFld != 0) srDel(iMagFld);
		if(iWfrSmp != 0) srDel(iWfrSmp);
		return ErrNo;
	}
	return 0;
}

////*************************************************************************
//
//void srEnsureSrwInit()
//{
//	char TestName[] = "SrwVerProc";
//	char StrDummy[50];
//	int SRW_WasNotInit = FetchStrVar(TestName, StrDummy);
//	if(SRW_WasNotInit)
//	{
//		char InitStr[] = "SrwInit(1)";
//		XOPCommand(InitStr);
//	}
//	SetXOPType(RESIDENT);
//}
//
////*************************************************************************
//
//static int srOptThinGenSetup(void* p_void)
//{
//	return SR.OptGenTransSetup((srTIgorOptGenTransSetupInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srOptMatConst(void* p_void)
//{
//	return SR.OptMatConst((srTIgorOptMatConstInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srOptZonePlateSpecSetup(void* p_void)
//{
//	return SR.OptZonePlateSpecSetup((srTIgorOptZonePlateSpecSetupInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srUtiInterTime(void* p_void)
//{
//	srYield.Init(((srTIgorUtiInterTimeInputStruct*)p_void)->Delta);
//	return 0;
//}
//
////*************************************************************************
//
//static int srUtiSpotInfo(void* p_void)
//{
//	return SR.UtiSpotInfo((srTIgorUtiSpotInfoInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srUtiWfrLimits(void* p_void)
//{
//	return SR.UtiWfrLimits((srTIgorUtiWrfLimitsInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srUtiRemoveFlips(void* p_void)
//{
//	return SR.UtiRemoveFlips((srTIgorUtiRemoveFlipsInputStruct*)p_void);
//}
//
////*************************************************************************
//
//static int srUtiMagRad(void* p_void)
//{
//	srTIgorUtiMagRadInputStruct* p = (srTIgorUtiMagRadInputStruct*)p_void;
//	double Rmag = SR.UtiMagRad(p->Bconst, p->ElecEnergy, p->OutUnit);
//	((srTIgorUtiMagRadInputStruct*)p_void)->Rmag = Rmag;
//	return 0;
//}
//
////*************************************************************************
//
//static int srUtiMagCritPhotEn(void* p_void)
//{
//	srTIgorUtiMagCritPhotEnInputStruct* p = (srTIgorUtiMagCritPhotEnInputStruct*)p_void;
//	double CritPhotEn = SR.UtiMagCritPhotEn(p->Bconst, p->ElecEnergy, p->OutUnit);
//	((srTIgorUtiMagCritPhotEnInputStruct*)p_void)->CritPhotEn = CritPhotEn;
//	return 0;
//}
//
////*************************************************************************
//
//static int srUtiUndK(void* p_void)
//{
//	srTIgorUtiUndKInputStruct* p = (srTIgorUtiUndKInputStruct*)p_void;
//	double K = SR.UtiUndK(p->Bpeak, p->Period);
//	((srTIgorUtiUndKInputStruct*)p_void)->K = K;
//	return 0;
//}
//
////*************************************************************************
//
//static int srUtiUndFundPhotEn(void* p_void)
//{
//	srTIgorUtiUndFundPhotEnInputStruct* p = (srTIgorUtiUndFundPhotEnInputStruct*)p_void;
//
//	int aha = 1;
//	double FundPhotEn = SR.UtiUndFundPhotEn(p->Bpeak, p->Period, p->ElecEnergy, p->OutUnit);
//	((srTIgorUtiUndFundPhotEnInputStruct*)p_void)->FundPhotEn = FundPhotEn;
//	return 0;
//}
//
////*************************************************************************
//
//static int srKn(void* p_void)
//{
//	return SR.Kmu((srTIgorKmuInputStruct*)p_void);
//}
//
////*************************************************************************
////TEST
//static int srUti3DView(void* p_void)
//{//to allow 2D drawing and viewing beamline feometries
//
//	//SR.Send.SetIgor3DViewInputStruct((srTIgor3DViewInputStruct*)p_void);
//	int result = 0;
//	
//	CHGenObj hSrc;
//	if(result = SR.Send.GetRadSourceData(p_void, "srTIgor3DViewInputStruct", hSrc)) return result;
//
//	srTEbmDat ElBeam;
//	//if(result = SR.Send.GetTotalElectronBeamDataFormat3(ElBeam)) return result;
//	srTMagElem MagField;
//
//
////ddddddddddddddddddddddddddddddddd
//	return SR.Uti3DView(ElBeam, MagField);
//}
//
////*************************************************************************
////TEST
//static int srRadPropagStokesMultiE(void* p_void)
//{
//	srTIgorRadPropagStokesMultiElecInputStruct* pRadPropagInputStruct = (srTIgorRadPropagStokesMultiElecInputStruct*)p_void;
//	SR.Send.SetIgorRadPropagStokesMultiElecInputStruct(pRadPropagInputStruct);
//	return SR.PropagateRadiationStokesMultiElec(pRadPropagInputStruct);
//}
//
////*************************************************************************
////TEST
//static int srWfrEmitPropag(void* p_void)
//{
//	SR.Send.SetIgorWfrEmitPropagInputStruct((srTIgorWfrEmitPropagInputStruct*)p_void);
//	return SR.ComputeWfrEmitPropag();
//}
//
////*************************************************************************
//
//static int srWfrReflect(void* p_void)
//{
//	SR.Send.SetIgorWfrReflectInputStruct((srTIgorWfrReflectInputStruct*)p_void);
//	return SR.ReflectWavefront();
//}
//
//*************************************************************************

static int srWfrAdd(void* p_void)
{
	int res = 0;
	srTIgorRadExtractInputStruct* pStr = (srTIgorRadExtractInputStruct*)p_void;

	srTSRWRadInData SRWRadInData;
	if(res = srTIgorSend::GetSRWRadInData(pStr->wRad, &SRWRadInData)) return res;

	//ddddddddddddddddddddddddddddddd

//	SR.Send.SetIgorWfrReflectInputStruct((srTIgorWfrReflectInputStruct*)p_void);
//	return SR.ReflectWavefront();
	return 0;
}

////*************************************************************************
///*OC
//static int srStokesArb(void* p_void)
//{
//	SR.Send.SetIgorStokesArbInputStruct((srTIgorStokesArbInputStruct*)p_void);
//	return SR.ComputeStokesArb();
//}
//OC*/

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
			return((long)srRadExtract);
			break;
		case 14:
			return((long)srStokesUnd);
			break;
		case 15:
			return((long)srStokesWig);
			break;
		case 16:
			return((long)srStokesConst);
			break;
		case 17:
			return((long)srPowDens);
			break;
		//case 18:
		//	return((long)srUtiInterTime);
		//	break;
		//case 19:
		//	return((long)srOptThinGenSetup);
		//	break;
		//case 20:
		//	return((long)srOptMatConst);
		//	break;
		//case 21:
		//	return((long)srOptZonePlateSpecSetup);
		//	break;
		//case 22:
		//	return((long)srUtiSpotInfo);
		//	break;
		//case 23:
		//	return((long)srUtiWfrLimits);
		//	break;
		//case 24:
		//	return((long)srUtiRemoveFlips);
		//	break;
		//case 25:
		//	return((long)srUtiMagRad);
		//	break;
		//case 26:
		//	return((long)srUtiMagCritPhotEn);
		//	break;
		//case 27:
		//	return((long)srUtiUndK);
		//	break;
		//case 28:
		//	return((long)srUtiUndFundPhotEn);
		//	break;
		//case 29:
		//	return((long)srKn);
		//	break;
		//case 30:
		//	return((long)srUti3DView);
		//	break;
		//case 31:
		//	return((long)srRadPropagStokesMultiE);
		//	break;
		//case 32:
		//	return((long)srWfrEmitPropag);
		//	break;
		//case 33:
		//	return((long)srWfrReflect);
		//	break;
		//case 34:
		//	return((long)srWfrChangeRep);
		//	break;
		case 35:
			return((long)srWfrAdd);
			break;
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
	int err = 0;

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
			err = srRadExtract(p);
			break;
		case 14:
			err = srStokesUnd(p);
			break;
		case 15:
			err = srStokesWig(p);
			break;
		case 16:
			err = srStokesConst(p);
			break;
		case 17:
			err = srPowDens(p);
			break;
		//case 18:
		//	err = srUtiInterTime(p);
		//	break;
		//case 19:
		//	err = srOptThinGenSetup(p);
		//	break;
		//case 20:
		//	err = srOptMatConst(p);
		//	break;
		//case 21:
		//	err = srOptZonePlateSpecSetup(p);
		//	break;
		//case 22:
		//	err = srUtiSpotInfo(p);
		//	break;
		//case 23:
		//	err = srUtiWfrLimits(p);
		//	break;
		//case 24:
		//	err = srUtiRemoveFlips(p);
		//	break;
		//case 25:
		//	err = srUtiMagRad(p);
		//	break;
		//case 26:
		//	err = srUtiMagCritPhotEn(p);
		//	break;
		//case 27:
		//	err = srUtiUndK(p);
		//	break;
		//case 28:
		//	err = srUtiUndFundPhotEn(p);
		//	break;
		//case 29:
		//	err = srKn(p);
		//	break;
		//case 30:
		//	err = srUti3DView(p);
		//	break;

		//case 31:
		//	err = srRadPropagStokesMultiE(p);
		//	break;
		//case 32:
		//	err = srWfrEmitPropag(p);
		//	break;
		//case 33:
		//	err = srWfrReflect(p);
		//	break;
		//case 34:
		//	err = srWfrChangeRep(p);
		//	break;
		case 35:
			err = srWfrAdd(p);
			break;
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

		//case IDLE:
		//	srEnsureSrwInit();
		//	break;
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
HOST_IMPORT void main(IORecHandle ioRecHandle)
{
	//#ifdef applec					// For MPW C for 68K only.
	//	void _DATAINIT(void);
	//	_DATAINIT();				// For MPW C only.
	//	UnloadSeg(_DATAINIT);
	//#endif
	
	#ifdef XOP_GLOBALS_ARE_A4_BASED
		#ifdef __MWERKS__
			SetCurrentA4();							// Set up correct A4. This allows globals to work.
			SendXOPA4ToIgor(ioRecHandle, GetA4());	// And communicate it to Igor.
		#endif
	#endif

	LoadXOPSegs();
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	if(igorVersion < 200) SetXOPResult(OLD_IGOR);
	else SetXOPResult(0L);
	
	SetXOPType(RESIDENT | IDLES);
	//SR.Initialize();
	//srYield.Init(0.5);

	srSetCompProgressIndicFunc(&(srTIgorSend::ShowCompProgress));
	srSetWfrExtModifFunc(&(srTIgorSend::WfrModify));
}

//*************************************************************************

// All structures are 2-byte-aligned.
#if GENERATINGPOWERPC
	#pragma options align=reset
#endif
#ifdef _WINDOWS_
	#pragma pack()
#endif
