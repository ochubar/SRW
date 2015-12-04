/************************************************************************//**
 * File: srsend.h
 * Description: Interface (input-output) functions essentially for IGOR Pro XOP (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRSEND_H
#define __SRSEND_H

#include "gmvect.h"

#ifdef __IGOR_PRO__
#include "srigintr.h"
#endif

#include <complex>
#include <time.h>
#include <iostream>
#include <sstream>

#include "stlstart.h" //OC030411

//#ifdef __VC__
//using namespace std;
//#endif
//#ifdef __MWERKS__
////using namespace std;
//#endif
//#ifdef __GCC__
//#define std
//#endif

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include "srstraux.h"
#include "srgtrjdt.h"

//*************************************************************************

extern int gCallSpinProcess;

//*************************************************************************

class srTTrjDat;
class srTRadInt;
class srTEbmDat;
class srTMagFieldPeriodic;
class srTMagFieldConstant;
class srTMagGroup;
class srTIsotrSrc;
class srTGsnBeam;
struct srTFunDer;
struct srTRadIntPerStoPrec;
struct srTRadIntWigPrec;
struct srTRadIntConstPrec;
struct srTZonePlateSpecNumAccessData;
struct srTRadIntPowDenPrec;
struct srTRadSASE;

//*************************************************************************

class srTSend {

#ifdef __IGOR_PRO__

	srTIgorTrjInputStruct* pIgorTrjInputStruct;
	srTIgorRadInputStruct* pIgorRadInputStruct;
	srTIgorRadFocusInputStruct* pIgorRadFocusInputStruct;
	srTIgorPerStokesInputStruct* pIgorPerStokesInputStruct;
	srTIgorStokesWigInputStruct* pIgorStokesWigInputStruct;
	srTIgorStokesConstInputStruct* pIgorStokesConstInputStruct;
	srTIgorStokesArbInputStruct* pIgorStokesArbInputStruct;
	srTIgorPowDensInputStruct* pIgorPowDensInputStruct;
	srTIgorIsotrSrcInputStruct* pIgorIsotrSrcInputStruct;
	srTIgorGsnBeamInputStruct* pIgorGsnBeamInputStruct;
	srTIgorWfrFromTrjInputStruct* pIgorWfrFromTrjInputStruct;
	srTIgorWfrSASEInputStruct* pIgorWfrSASEInputStruct;
	srTIgorRadPropagStokesMultiElecInputStruct* pIgorRadPropagStokesMultiElecInputStruct;
    srTIgorRadEmitPropagStokesMultiElecInputStruct* pIgorRadEmitPropagStokesMultiElecInputStruct;
	srTIgorWfrEmitPropagInputStruct* pIgorWfrEmitPropagInputStruct;
	srTIgorWfrReflectInputStruct* pIgorWfrReflectInputStruct;
	//srTIgor3DViewInputStruct* pIgor3DViewInputStruct;
	
	srTHandleOfSRWRadStruct HandleOfSRWRadStruct;
	
	DOUBLE *pOutBtxData, *pOutXData, *pOutBtzData, *pOutZData;
	int hStateOutBtxData, hStateOutXData, hStateOutBtzData, hStateOutZData;

	waveHndl wavH_OutBtxData, wavH_OutXData, wavH_OutBtzData, wavH_OutZData;

	float *pOutExData, *pOutEzData, *tOutExData, *tOutEzData;
	int hStateOutExData, hStateOutEzData;
	waveHndl wavH_OutExData, wavH_OutEzData;

	waveHndl wavH_AuxElectronBeam, wavH_AuxTrajectory;

	char ProgressIndicatorWinName[256];
	long TotalAmOfOutPoints;
	clock_t StartLoopClock, PrevUpdateClock;
#endif

	char ProgressIndicatorShouldBeUsed;
	char ProgressIndicatorIsUsed;

public:
	srTSend() 
	{
		ProgressIndicatorIsUsed = 0;
		ProgressIndicatorShouldBeUsed = 0;

#ifdef __VC__
		ProgressIndicatorShouldBeUsed = 1;
#endif
	}

	void ErrorMessage(const char*);
	void OrdinaryMessage(const char*);
	void WarningMessage(const char*);
	void WarningMessages(srTIntVect*);

	void AddWarningMessage(srTIntVect*, int);

	void String(char*);
	void Long(long);
	void Int(int);
	void IntList(int*, int);
	void Double(double);
	void DoubleList(double*, int);
	void DoubleListWithArg(double, double, double*, int);

	void InitOutList(int);

	void ArbNestedArrays(double*, int*, int);
	void SubArbNestedArrays(double*, int*, int, int&);

	void MyMLPutDouble(double);

	int GetTotalElectronData(srTTrjDat&);

	int OutRadDistrFormat1(srTRadInt&);
	int OutRadDistrFormat2(srTSRWRadStructAccessData&, srTRadInt&);

	void OutPhaseFormat1(srTRadInt&);

	int InitTrjOutFormat1(srTTrjDat&);
	inline void FinishTrjOutFormat1();

	int InitRadDistrOutFormat1(srTWfrSmp&);
	int InitRadDistrOutFormat2(srTWfrSmp&);
	int InitRadDistrOutFormat3(srTSRWRadStructAccessData&, srTWfrSmp&);

	inline void FinishRadDistrOutFormat1();
	inline int RadValDirectOut(complex<double>*);

	inline int RadValDirectOut(srTEFourier&);
	inline int RadValDirectOut_FluxDens(srTEFourier&);

	int GetTotalElectronBeamDataFormat2(srTTrjDat&);
	int GetTotalElectronBeamDataFormat3(srTEbmDat&);
	int GetTotalElectronDistribDataFormat1(srTEbmDat&);

	int GetIsotrSrcExtraDataFormat1(srTIsotrSrc&);
	int GetTotalGsnBeamDataFormat1(srTGsnBeam&);

	int GetTotalFieldDataFormat1(srTMagFieldAccessData&, srTTrjDat&);
	int GetTotalFieldDataFormat2(srTTrjDat&);

	int GetPeriodicFieldDataFormat1(srTMagFieldPeriodic&, waveHndl h=0);
	int GetConstantFieldDataFormat1(srTMagFieldConstant&);
	int GetGenMagFieldDataFormat1(srTGenTrjHndl&);
	int IdentifyMagFieldType(char&);
	int IdentifyMagFieldTypeFromName(char*);

	int GetGeneralMagneticFieldDataFormat1(srTMagGroup&);

	int GetTotalTrajectoryDataFormat1(srTTrjDat&);
	int GetTrajectoryComponentDataFormat1(srTWaveAccessDataD1D&);

	int GetSASEInputRadDataFormat1(srTRadSASE&, srTSRWRadStructAccessData&);
	int GetSASEPrecDataFormat1(srTPrecSASE&);
	int GetSRWRadStructArray(srTSRWRadStructAccessData*& arSRWRadStructAccessData, int& numRadStruct);
	int SetupSASEControlStruct(srTControlAccessSASE&, int numHarm);
	int SetupControlStruct1D(waveHndl wText, long* Indices, long NewNp, DOUBLE NewStart, DOUBLE NewStep, waveHndl& wHndl, char*& pBase, int& hState);
	int SetupControlStructMD(waveHndl wText, long* Indices, long* NewNpAr, DOUBLE* NewStartAr, DOUBLE* NewStepAr, int numDim, waveHndl& wHndl, char*& pBase, int& hState);

	int GetPropagRadStokesMultiElecDataFormat1(double* pPrecPar);
	int GetWfrEmitPropagPrec(double* pPrecPar);

	int FinishWorkingWithControlSASEStruct(srTControlAccessSASE& ControlAccessSASE);

	int GetTotalObservationDataFormat1(srTWfrSmp&);
	int GetTotalObservationDataFormat2(srTWfrSmp&); // For transformed Light; May be changed
	int GetTotalObservationDataFormat3(srTWfrSmp&); // To read directly from SRW Observation wave

	int GetTotalRadIntegrationParamDataFormat1(srTRadInt&);
	int GetRadIntPeriodicParamDataFormat1(srTRadIntPerStoPrec&);
	int GetRadIntWigglerParamDataFormat1(srTRadIntWigPrec&);
	int GetRadIntConstParamDataFormat1(srTRadIntConstPrec&);
	int GetRadIntPowDensParamDataFormat1(srTRadIntPowDenPrec&);

	int GetAuxObsTreatParamFormat1(srTWfrSmp&);

	int GetOptElemNames(srTStringVect&);

	//int GetTotalOpticalElemDataFormat1(void*); //Obsolete
	int OutOpticsIncRadDistrFormat1(double*, int, double, double); // Temporary

	int GetWavelengthAndCoordData(srTWfrSmp&, char);
	int GetTotalOutTrjParam(short*, double*, short&);
	int GetCoordData(double&, double&, int&);
	int GetPrtclTrjInitData(double&, double&, double&, double&, double&);
	int GetString(char*&, int);
	int GetDouble(double&);

	int GetVectorOfStrings(char*, srTStringVect*);
	int GetVectorOfStrings(waveHndl&, srTStringVect*);
	int FinishWorkingWithWave(srTWaveAccessData*);
	int FetchNumWave(char*, srTWaveAccessData*);
	int ModifyRadNeNxNz(srTSRWRadStructAccessData&, char =0);
	static int ModifyStokesNeNxNz(srTStokesStructAccessData&);
	int GetRadStructNames(srTSRWRadStructAccessData&, srTSRWRadStructWaveNames&);
	int CreateNewRadStruct(srTSRWRadStructAccessData&, srTSRWRadStructWaveNames&);
	int DeleteRadStructWaves(srTSRWRadStructAccessData&, srTSRWRadStructWaveKeys&);
	int RenameRadStruct(srTSRWRadStructAccessData&, char*);
	int RenameRadStruct(srTSRWRadStructAccessData&, srTSRWRadStructWaveNames&);
	int MakeWaveAccordingToWaveAccessData(srTWaveAccessData&);

	void DeleteStringVect(srTStringVect& OptElemInfo)
	{
		for(int k=0; k<(int)(OptElemInfo.size()); k++)
		{
			char* aStr = OptElemInfo[k]; 
			if(aStr != 0) delete[] aStr;
		}
	}

	int SetUpPhaseShiftWave(srTWaveAccessData& PhShData, int FunNo)
	{
#ifdef __IGOR_PRO__
		ostringstream OutStream;
		OutStream << "SrwOptPhaseShiftLoop(" << FunNo << ", \"" << PhShData.NameOfWave << "\");" << "\r" << ends;
		basic_string<char> BufString = OutStream.str();
		const char* OutString = BufString.c_str();

		HSetState((Handle)(PhShData.wHndl), PhShData.hState);
		WaveHandleModified(PhShData.wHndl);

		int result;
		if(result = XOPSilentCommand(OutString)) return result;

		long dataOffset;
		if(result = MDAccessNumericWaveData(PhShData.wHndl, kMDWaveAccessMode0, &dataOffset)) return result;
		PhShData.hState = 0; //MoveLockHandle(PhShData.wHndl);
		PhShData.pWaveData = (char*)(*(PhShData.wHndl)) + dataOffset;
		return 0;
#else
		//todo
		return 0;
#endif
	}
	void ShowOutTrjDataPointers(DOUBLE*& Out_pOutBtxData, DOUBLE*& Out_pOutXData, DOUBLE*& Out_pOutBtzData, DOUBLE*& Out_pOutZData)
	{
#ifdef __IGOR_PRO__
		Out_pOutBtxData = pOutBtxData; Out_pOutXData = pOutXData; Out_pOutBtzData = pOutBtzData; Out_pOutZData = pOutZData;
#endif
	}
	int GetSRWRadStructAccessData(srTSRWRadStructAccessData*);

#ifdef __IGOR_PRO__
	int GetSRWRadStructAndResizeData(srTIgorRadResizeInputStruct*, srTSRWRadStructAccessData*, srTRadResize*);
#endif
	
	int FinishWorkingWithSRWRadStruct(srTSRWRadStructAccessData*);

#ifdef __IGOR_PRO__
	srTIgorRadInputStruct* GetIgorRadInputStructPtr() 
	{ 
		return pIgorRadInputStruct;
//#else
		//todo
		//return 0;
	}
	void SetIgorRadInputStruct(srTIgorRadInputStruct* In_pIgorRadInputStruct)
	{
//#ifdef __IGOR_PRO__
		ZeroAllInputStructPtrs();
		pIgorRadInputStruct = In_pIgorRadInputStruct;
		HandleOfSRWRadStruct.wRad = pIgorRadInputStruct->wRad;
//#else
		//todo
	}
	int ChangeObservationData(srTHandleOfOneWaveStruct*, srTWfrSmp&);
	int GetSRWRadStructAndExtractData(srTIgorRadExtractInputStruct*, srTSRWRadStructAccessData*, srTRadExtract*);
	int GetWaveAccessData(srTHandleOfOneWaveStruct*, srTWaveAccessData*);
	int GetSRWRadStructAndOptElemNames(srTIgorRadPropagInputStruct*, srTSRWRadStructAccessData*, srTStringVect*);
	int GetSRWRadStructAndOptElemNames(srTIgorRadPropagStokesMultiElecInputStruct*, srTSRWRadStructAccessData*, srTStringVect*);

#endif

	int Finish2DFFT(srTWaveAccessData* pWaveAccessData)
	{
#ifdef __IGOR_PRO__
		int result;
		DOUBLE xStart = pWaveAccessData->DimStartValues[0];
		DOUBLE xStep = pWaveAccessData->DimSteps[0];
		if(result = MDSetWaveScaling(pWaveAccessData->wHndl, ROWS, &xStep, &xStart)) return result;

		DOUBLE yStart = pWaveAccessData->DimStartValues[1];
		DOUBLE yStep = pWaveAccessData->DimSteps[1];
		if(result = MDSetWaveScaling(pWaveAccessData->wHndl, COLUMNS, &yStep, &yStart)) return result;

		HSetState((Handle)(pWaveAccessData->wHndl), pWaveAccessData->hState);
		WaveHandleModified(pWaveAccessData->wHndl);
		return 0;
#else
		//todo
		return 0;
#endif
	}
	int FinishWorkingWithTrajectoryData(srTTrjDat&);
	
	int GetStokesStructAccessData(srTStokesStructAccessData*);
	int FinishWorkingWithStokesStruct(srTStokesStructAccessData* pStokesAccessData)
	{
#ifdef __IGOR_PRO__
		HSetState((Handle)(pStokesAccessData->wSto), pStokesAccessData->hStateSto);
		WaveHandleModified(pStokesAccessData->wSto);
		return 0;
#else
		//todo
		return 0;
#endif
	}
	int GetPowDensStructAccessData(srTPowDensStructAccessData*);
	int FinishWorkingWithPowDensStruct(srTPowDensStructAccessData* pPowDensAccessData)
	{
#ifdef __IGOR_PRO__
		HSetState((Handle)(pPowDensAccessData->wPowDens), pPowDensAccessData->hStatePowDens);
		WaveHandleModified(pPowDensAccessData->wPowDens);
		return 0;
#else
		//todo
		return 0;
#endif
	}
	int FinishWorkingWithElecDistrStruct(srTEbmDat& EbmDat)
	{
#ifdef __IGOR_PRO__
		return ReleaseWave(EbmDat.wElecDistr, EbmDat.hStateElecDistr);
#else
		return 0;
#endif
	}
	int UpdateTextWave(waveHndl& wavH, srTStringVect* pStringVect);
	int DeleteWave(srTWaveAccessData& WaveAccessData)
	{
#ifdef __IGOR_PRO__
		int result;
		HSetState((Handle)(WaveAccessData.wHndl), WaveAccessData.hState);
		WaveHandleModified(WaveAccessData.wHndl);
		if(result = KillWave(WaveAccessData.wHndl)) return result;
		WaveAccessData.pWaveData = 0;
		return 0;
#else
		//todo
		return 0;
#endif
	}

#ifdef __IGOR_PRO__
	void SetHandleOfSRWRadStruct(srTHandleOfSRWRadStruct* In_pHandleOfSRWRadStruct)
	{
		HandleOfSRWRadStruct = *In_pHandleOfSRWRadStruct;
//#else
		//todo
	}
	void ZeroAllInputStructPtrs()
	{
		HandleOfSRWRadStruct.wRad = NIL;

		pIgorTrjInputStruct = 0;
		pIgorRadInputStruct = 0;
		pIgorRadFocusInputStruct = 0;
		pIgorPerStokesInputStruct = 0;
		pIgorStokesWigInputStruct = 0;
		pIgorStokesConstInputStruct = 0;
		pIgorStokesArbInputStruct = 0;
		pIgorPowDensInputStruct = 0;
		pIgorIsotrSrcInputStruct = 0;
		pIgorGsnBeamInputStruct = 0;
		pIgorWfrFromTrjInputStruct = 0;
		pIgorWfrSASEInputStruct = 0;
		pIgorRadPropagStokesMultiElecInputStruct = 0;
		pIgorWfrEmitPropagInputStruct = 0;
		pIgorWfrReflectInputStruct = 0;
		//pIgor3DViewInputStruct = 0;
	}
#endif

	void UpdateInterface();


#ifdef __IGOR_PRO__

	void SetIgorTrjInputStruct(srTIgorTrjInputStruct* In_pIgorTrjInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorTrjInputStruct = In_pIgorTrjInputStruct;
	}
	void SetIgorRadFocusInputStruct(srTIgorRadFocusInputStruct* In_pIgorRadFocusInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorRadFocusInputStruct = In_pIgorRadFocusInputStruct;
	}
	void SetIgorPerStokesInputStruct(srTIgorPerStokesInputStruct* In_pIgorPerStokesInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorPerStokesInputStruct = In_pIgorPerStokesInputStruct;
	}
	void SetIgorStokesWigInputStruct(srTIgorStokesWigInputStruct* In_pIgorStokesWigInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorStokesWigInputStruct = In_pIgorStokesWigInputStruct;
	}
	void SetIgorStokesConstInputStruct(srTIgorStokesConstInputStruct* In_pIgorStokesConstInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorStokesConstInputStruct = In_pIgorStokesConstInputStruct;
	}
	void SetIgorStokesArbInputStruct(srTIgorStokesArbInputStruct* In_pIgorStokesArbInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorStokesArbInputStruct = In_pIgorStokesArbInputStruct;
	}
	void SetIgorPowDensInputStruct(srTIgorPowDensInputStruct* In_pIgorPowDensInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorPowDensInputStruct = In_pIgorPowDensInputStruct;
	}
	void SetIgorIsotrSrcInputStruct(srTIgorIsotrSrcInputStruct* In_pIgorIsotrSrcInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorIsotrSrcInputStruct = In_pIgorIsotrSrcInputStruct;
		HandleOfSRWRadStruct.wRad = pIgorIsotrSrcInputStruct->wRad;
	}
	void SetIgorGsnBeamInputStruct(srTIgorGsnBeamInputStruct* In_pIgorGsnBeamInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorGsnBeamInputStruct = In_pIgorGsnBeamInputStruct;
		HandleOfSRWRadStruct.wRad = pIgorGsnBeamInputStruct->wRad;
	}
	void SetIgorWfrFromTrjInputStruct(srTIgorWfrFromTrjInputStruct* In_pIgorWfrFromTrjInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorWfrFromTrjInputStruct = In_pIgorWfrFromTrjInputStruct;
		HandleOfSRWRadStruct.wRad = pIgorWfrFromTrjInputStruct->wRad;
	}
	void SetIgorWfrSASEInputStruct(srTIgorWfrSASEInputStruct* In_pIgorWfrSASEInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorWfrSASEInputStruct = In_pIgorWfrSASEInputStruct;
		HandleOfSRWRadStruct.wRad = pIgorWfrSASEInputStruct->wRad;
	}
	void SetIgorRadPropagStokesMultiElecInputStruct(srTIgorRadPropagStokesMultiElecInputStruct* In_pIgorRadPropagStokesMultiElecInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorRadPropagStokesMultiElecInputStruct = In_pIgorRadPropagStokesMultiElecInputStruct;
	}
	void SetIgorRadEmitPropagStokesMultiElecInputStruct(srTIgorRadEmitPropagStokesMultiElecInputStruct* In_pIgorRadEmitPropagStokesMultiElecInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorRadEmitPropagStokesMultiElecInputStruct = In_pIgorRadEmitPropagStokesMultiElecInputStruct;
	}
	void SetIgorWfrEmitPropagInputStruct(srTIgorWfrEmitPropagInputStruct* In_pIgorWfrEmitPropagInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorWfrEmitPropagInputStruct = In_pIgorWfrEmitPropagInputStruct;
	}
	void SetIgorWfrReflectInputStruct(srTIgorWfrReflectInputStruct* In_pIgorWfrReflectInputStruct)
	{
		ZeroAllInputStructPtrs();
		pIgorWfrReflectInputStruct = In_pIgorWfrReflectInputStruct;
	}
	//void SetIgor3DViewInputStruct(srTIgor3DViewInputStruct* In_pIgor3DViewInputStruct)
	//{
	//	ZeroAllInputStructPtrs();
	//	pIgor3DViewInputStruct = In_pIgor3DViewInputStruct;
	//}

	srTIgorRadFocusInputStruct* GetIgorRadFocusInputStructPtr() { return pIgorRadFocusInputStruct;}

	int UpdateNumberPositionInSRWRad(srTSRWRadStructAccessData* pSRWRadStructAccessData, long* RadIndices, char* CharBuf, int AmOfBytes)
	{
		char* tCharBuf = CharBuf;

		Handle textH = NewHandle(AmOfBytes);
		strncpy(*textH, CharBuf, AmOfBytes);
		int result;
		if(result = MDSetTextWavePointValue(pSRWRadStructAccessData->wRad, RadIndices, textH)) return result;
		DisposeHandle(textH);

		for(int i=0; i<15; i++) *(tCharBuf++) = 0;
		return 0;
	}
	int UpdateTextPositionInSRWRad(srTSRWRadStructAccessData* pRad, int Index, char* Text)
	{
		int result;
		long RadIndices[MAX_DIMENSIONS]; RadIndices[0] = Index;
		int AmOfBytes = (int)strlen(Text);
		Handle textH = NewHandle(AmOfBytes);
		strncpy(*textH, Text, AmOfBytes);
		if(result = MDSetTextWavePointValue(pRad->wRad, RadIndices, textH)) return result;
		DisposeHandle(textH);
		return 0;
	}
	int UpdateTextWave1DPointValue(waveHndl& wavH, int Index, char* Text)
	{
		int result;
		long Indices[MAX_DIMENSIONS]; Indices[0] = Index;
		int AmOfBytes = (int)strlen(Text);
		Handle textH = NewHandle(AmOfBytes);
		strncpy(*textH, Text, AmOfBytes);
		if(result = MDSetTextWavePointValue(wavH, Indices, textH)) return result;
		DisposeHandle(textH);
		return 0;
	}

	int GetDoubleFromTextWave1D(waveHndl& wavH, int IndNo, double& Value)
	{
		long RadIndices[MAX_DIMENSIONS];
		RadIndices[0] = IndNo;
		Handle textH = NewHandle(0L);
		int result;
		if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
		if(textH != 0)
		{
			if(*textH != 0)
			{
				Value = atof(*textH);
				if((Value == 0.) && (**textH != '0') && (**textH != ' ')) 
				{
					DisposeHandle(textH);
					return ERROR_IN_READING_NUMBER_FROM_TEXT_WAVE;
				}
			}
		}
		DisposeHandle(textH);
		return 0;
	}

	double IgorStrToDouble(char* Str)
	{
		char CharBuf[256];
		char* tCharBuf = CharBuf;
		char* pCh = Str;
		for(int i=0; i<255; i++) 
		{
			if(!(isdigit(*pCh) || (*pCh == '.') || (*pCh == '+') || (*pCh == '-') || (*pCh == 'e') || (*pCh == 'E'))) break;
			*(tCharBuf++) = *(pCh++);
		}
		*tCharBuf = '\0';
		return atof(CharBuf);
	}

	int ReleaseWave(waveHndl wHndl, int hState)
	{
		if(wHndl != NIL)
		{
			HSetState((Handle)wHndl, hState);
			WaveHandleModified(wHndl);
		}
		return 0;
	}

	int GetDataPointerAndStateAndLockNumericWave(waveHndl wavH, char*& DataPtr, int& hState)
	{
		long dataOffset, result = 0;
		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		hState = 0; //MoveLockHandle(wavH);
		DataPtr = (char*)(*wavH) + dataOffset;

		return result;
	}

	int GetTotalElectronBeamDataFormat3_FromIgorWave(srTEbmDat& EbmDat, waveHndl wavH);
	int GetRadSourceData(void* pVoidIn, char* StructNameIn, CHGenObj& hSrcOut);

	int ShowRadSect1D(srTRadSect1D&, char);
	int MakeWaveToShowData(srTWaveAccessData&);

#endif

	bool CheckIfRadStructWasSet()
	{
#ifdef __IGOR_PRO__
		if(HandleOfSRWRadStruct.wRad == NIL) return false;

		Handle textH = NewHandle(0L);
		long RadIndices[MAX_DIMENSIONS];
		for(int i=0; i<MAX_DIMENSIONS; i++) RadIndices[i] = 0;
		
		RadIndices[0] = 0;
		if(MDGetTextWavePointValue(HandleOfSRWRadStruct.wRad, RadIndices, textH)) return false;
		if(PtrAndHand("\0", textH, 1)) return false;
		int len = (int)strlen(*textH);
		if(len < 5) return false;
		if(strcmp(*textH + len - 5, "X_rae") != 0) return false;
		DisposeHandle(textH);

		textH = NewHandle(0L);
		RadIndices[0] = 1;
		if(MDGetTextWavePointValue(HandleOfSRWRadStruct.wRad, RadIndices, textH)) return false;
		if(PtrAndHand("\0", textH, 1)) return false;
		len = (int)strlen(*textH);
		if(len < 5) return false;
		if(strcmp(*textH + len - 5, "Z_rae") != 0) return false;
		DisposeHandle(textH);
#endif

		return true;
	}

};

//*************************************************************************

inline void srTSend::FinishTrjOutFormat1()
{
#if defined(__MATHEMATICA__)

#elif defined(__IGOR_PRO__)

	HSetState((Handle)wavH_OutBtxData, hStateOutBtxData);
	WaveHandleModified(wavH_OutBtxData);
	HSetState((Handle)wavH_OutXData, hStateOutXData);
	WaveHandleModified(wavH_OutXData);
	HSetState((Handle)wavH_OutBtzData, hStateOutBtzData);
	WaveHandleModified(wavH_OutBtzData);
	HSetState((Handle)wavH_OutZData, hStateOutZData);
	WaveHandleModified(wavH_OutZData);

#endif
}

//*************************************************************************

inline void srTSend::FinishRadDistrOutFormat1()
{
#if defined(__MATHEMATICA__)

#elif defined(__IGOR_PRO__)

	wavH_OutExData = NIL;
	wavH_OutEzData = NIL; // The handles are released by FinishWorkingWithSRWRadStruct

#ifdef __VC__
	if(ProgressIndicatorIsUsed)
	{
		char TotOutStr[300];
		strcpy(TotOutStr, "DoWindow/K ");
		strcat(TotOutStr, ProgressIndicatorWinName);
		XOPCommand(TotOutStr);
	}
#endif
#endif
}

//*************************************************************************

inline int srTSend::RadValDirectOut(complex<double>* RadVal)
{
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	*(tOutExData++) = float(RadVal->real());
	*(tOutExData++) = float((RadVal++)->imag());
	*(tOutEzData++) = float(RadVal->real());
	*(tOutEzData++) = float(RadVal->imag());

	if(gCallSpinProcess && SpinProcess()) return SR_COMP_PROC_ABORTED;
	return 0;

#endif
}

//*************************************************************************

inline int srTSend::RadValDirectOut(srTEFourier& RadVal)
{
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	*(tOutExData++) = float(RadVal.EwX_Re);
	*(tOutExData++) = float(RadVal.EwX_Im);
	*(tOutEzData++) = float(RadVal.EwZ_Re);
	*(tOutEzData++) = float(RadVal.EwZ_Im);

	if(gCallSpinProcess && SpinProcess()) return SR_COMP_PROC_ABORTED;
	return 0;

#endif
}

//*************************************************************************

inline int srTSend::RadValDirectOut_FluxDens(srTEFourier& RadVal)
{
#if defined(__IGOR_PRO__)

	float xRe = float(RadVal.EwX_Re);
	float xIm = float(RadVal.EwX_Im);
	*(tOutExData++) = xRe*xRe + xIm*xIm;
	float zRe = float(RadVal.EwZ_Re);
	float zIm = float(RadVal.EwZ_Im);
	*(tOutEzData++) = zRe*zRe + zIm*zIm;

	if(gCallSpinProcess && SpinProcess()) return SR_COMP_PROC_ABORTED;
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

#endif
