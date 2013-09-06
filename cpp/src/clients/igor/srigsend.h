/************************************************************************//**
 * File: srigsend.h
 * Description: Interface (header), mostly used for communications with IGOR Pro (WaveMetrics)
 * Project: Synchrotron Radiation Workshop (SRW)
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRIGSEND_H
#define __SRIGSEND_H

#include "srigintr.h"
#include "stlstart.h"
#include <vector>
#include <string>

//*************************************************************************

struct srTSRWRadInData;
struct srTSRWStokesInData;
struct srTSRWPowDensInData;
struct srTIgorWaveAccessData;
struct srTMagFldHarm;
struct srTParPrecStokesArb;
struct srTDataMD;

//*************************************************************************

class srTIgorSend {

	static char gProgressIndicatorWinName[256];
	static long gProgressIndicatorScale;

public:

	static void WarningMessage(const char* WarnText);
	static int GetDoubleFromTextWave1D(waveHndl& wavH, int IndNo, double& Value);
	static int UpdateTextWave1D(waveHndl& wavH, vector<string>* pVectStr);

	static int GetArrDoubleFromNumWave1D(waveHndl wavH, long MaxNp, double*& pData, long& Np);
	static int GetArrDoubleFromNumWave2D(waveHndl wavH, int& NumRows, int& NumCols, double*& pData);

	static int ReDimNumWave1D(waveHndl wavH, int NumRows);
    static int ReDimNumWave2D(waveHndl wavH, int NumRows, int NumCols);
    static int SetDataInNumWave(waveHndl wavH, double* pData, long Np);

	static int FetchNumWave(char* sWaveName, srTDataMD* pWaveData);
	static int GetNumWaveData(waveHndl wavH, srTDataMD* pWaveData); //, int& hState);
	static int FinishWorkingWithWave(srTDataMD* pWaveData, waveHndl wavH); //, int hState);

	//static int GetElecBeamThin(waveHndl wavH, double& I, double* Mom1, double& s);
	static int GetElecBeamThin(waveHndl wavH, double& I, double& Neb, double* Mom1, double& s);
	static int SetElecBeamThin(waveHndl wavH, double I, double* Mom1, double s);
	//static int GetElecBeamThick(waveHndl wavH, double& I, double* pMom1, double* pMom2, int& nMom2, double& s0, int& TypeDistrTransverse, int& TypeDistrLongitudinal, double& ShortNoiseFactor);
	static int GetElecBeamThick(waveHndl wavH, double& I, double& Neb, double* pMom1, double* pMom2, int& nMom2, double& s0, int& TypeDistrTransverse, int& TypeDistrLongitudinal, double& ShortNoiseFactor);
    static int SetElecBeamThick(waveHndl wavH, double I, double* pMom1, double* pMom2, double s0, int TypeDistrTransverse, int TypeDistrLongitudinal, double ShortNoiseFactor);
	static int GetElecBeamTwiss(waveHndl wavH, double* pHorTwiss, double* pVertTwiss);
	static int SetElecBeamEmitAndTwiss(waveHndl wavH, double HorEmit, double VertEmit, double* pHorTwiss, double* pVertTwiss);
	
	static int GetMagFieldTransvUnif(waveHndl wavH, double FieldAbsZeroTol, double& sStart, double& sStep, long& np, double*& pBH, bool& BHIsZero, double*& pBV, bool& BVIsZero);
	static int GetMagField3d(waveHndl wField, double FieldZeroTolerance, double& sStart, double& sStep, long& NpS, double& xStart, double& xStep, long& NpX, double& zStart, double& zStep, long& NpZ , double*& pBsData, bool& BsIsZero, double*& pBxData, bool& BxIsZero, double*& pBzData, bool& BzIsZero);
	static int GetMagField3dComponent(waveHndl wavH, double FieldZeroTolerance, double& sStart, double& sStep, long& NpS, double& xStart, double& xStep, long& NpX, double& zStart, double& zStep, long& NpZ, double*& pB_Data, bool& B_IsZero);

	static int GetMagFieldPeriodic(waveHndl wavH, double& PerLength, double& TotLength, int& AmOfHarm, srTMagFldHarm*& MagFldHarmArr, int& TypeOfUnd, double& TaperPar_TU, double& PhaseSh_OK, int& FldErrTypeSASE, double& FldErrRMS, double& NatFocNxSASE, double& NatFocNySASE, int& TaperTypeSASE, double& TaperStartSASE, double& TaperRelFldChgSASE);
	static int GetMagFieldPeriodic(waveHndl wavH, double& PerLength, double& TotLength, int& AmOfHarm, int*& ArrHarmNo, char*& ArrXorZ, double*& ArrK, double*& ArrPhase, int& TypeOfUnd, double& TaperPar_TU, double& PhaseSh_OK, int& FldErrTypeSASE, double& FldErrRMS, double& NatFocNxSASE, double& NatFocNySASE, int& TaperTypeSASE, double& TaperStartSASE, double& TaperRelFldChgSASE);
	static int SetMagFieldPeriodic(waveHndl wavH, double PerLength, double TotLength, int AmOfHarm, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase, int TypeOfUnd, double TaperPar_TU, double PhaseSh_OK, double Fund_keV_per_GeV2);
	static int GetMagFieldConstant(waveHndl wavH, double& Bx, double& Bz);
	static int IdentifyMagFieldType(waveHndl wavH, int& MagFieldType);
	static int IdentifyMagFieldTypeFromName(char* MagFieldName);
	static int GetTrjDataPointers(waveHndl wOutHorAng, waveHndl wOutHorCoor, waveHndl wOutVerAng, waveHndl wOutVerCoor, double*& pOutBtxData, double*& pOutXData, double*& pOutBtzData, double*& pOutZData, int& ns, double& sStart, double& sStep, int& hStateOutBtxData, int& hStateOutXData, int& hStateOutBtzData, int& hStateOutZData);
	static int GetTrj3dDataPointers(waveHndl wavH_OutBtxData, waveHndl wavH_OutXData, waveHndl wavH_OutBtyData, waveHndl wavH_OutYData, waveHndl wavH_OutBtzData, waveHndl wavH_OutZData, double*& pOutBtxData, double*& pOutXData, double*& pOutBtyData, double*& pOutYData, double*& pOutBtzData, double*& pOutZData, int& ns, double& sStart, double& sStep, int& hStateOutBtxData, int& hStateOutXData, int& hStateOutBtyData, int& hStateOutYData, int& hStateOutBtzData, int& hStateOutZData);
	static int GetPrecParamWfrSamplingForPropag(waveHndl wavH, bool& AllowAutoChoiceOfNxNzForPropag, double& NxNzOversamplingParam);
	static int GetWfrSampling(waveHndl wavH, double& s, double& zSt, double& zFi, int& nz, double& xSt, double& xFi, int& nx, double& eSt, double& eFi, int& ne, char* PhotEnUnits, double& tSt, double& tFi, int& nt, int& presT, double*& pSurfData, waveHndl& wSurfData, int& hStateSurfData, double* horOrtObsPlane =0, double* inNormObsPlane =0);
	//static int GetPrecParamElectricFieldComp(waveHndl wavH, int& IntegMeth, double& RelPrecOrStep, double& sStart, double& sEnd);
	static int GetPrecParamElectricFieldComp(waveHndl wavH, int& IntegMeth, double& RelPrecOrStep, double& sStart, double& sEnd, bool& showProgrIndic);
	//static int GetPrecParamStokesPerComp(waveHndl wavH, int& InitHarm, int& FinHarm, double& Kns, double& Knphi, char& IntensityOrFlux);
	static int GetPrecParamStokesPerComp(waveHndl wavH, int& InitHarm, int& FinHarm, double& Kns, double& Knphi, char& IntensityOrFlux, double& minPhotEnExtRight);
	static int GetPrecParamPowDensComp(waveHndl wavH, double& PrecFact, int& Method, int& UseSpecIntLim, double& sIntStart, double& sIntFin);
	static int GetPrecParamMagArb2Per(waveHndl wavH, double& RelPrec, int& MaxAmOfHarm, double& MaxPerLen_m);
	static int GetPrecParamStokesArbComp(waveHndl wavH, srTParPrecStokesArb* pPrecStokesArb);
	static int GetSRWRadInData(waveHndl wavH, srTSRWRadInData* pSRWRadInData);
	static int GetSRWStokesInData(waveHndl wavH, srTSRWStokesInData* pStokesAccessData);
	static int GetSRWPowDensInData(waveHndl wavH, srTSRWPowDensInData* pSRWPowDensInData);
	static int GetIntensExtractParam(waveHndl wExtractParam, int& PolarizCompon, int& Int_or_Phase, int& PlotType, int& TransvPres, double& ePh, double& x, double& z);
    static int GetIntensExtractData(waveHndl wExtractedData, int& hStateExtractedData, char*& pExtractedData);
	static int GetVectorOfStrings(waveHndl wavH, vector<string>*);
	static int GetVectorOfStrings(const char* StructName, vector<string>*);
	static int GetOptElemInfByName(const char* sNameOptElem, char** pDescrStr, int* LenDescr, void*);

	static int GetGaussianBeam(waveHndl& wavH, double& Phot, int& Polar, double& SigmaX, double& SigmaZ, int& mx, int& mz, double* Mom1Arr, double& s0, double& SigmaT, int& TypeDistrT, double& RepRate, double& PulseEn, double& AvgPhotEn);

    static int GetAndSetElecBeamThick(int* piElecBeam, waveHndl& wavH);
	static int GetAndSetElecBeamThin(int* piElecBeam, waveHndl& wavH);
	static int GetAndSetElecBeamGen(int* piElecBeam, waveHndl& wavH);
	static int GetAndSetElecBeamCont(int* piElecBeam, waveHndl& wavH);

	static int GetAndSetGaussianBeam(int* piGsnBeam, waveHndl& wavH);
    static int GetAndSetMagFieldGen(int* piMagFld, waveHndl& wavH);
	static int GetAndSetWfrSampling(int* piWfrSmp, waveHndl& wavH);
	static int GetAndSetWfr(int* piWfr, waveHndl& wavH);

	static int GetAndSetOptElem(int* piOptElem, int* piWfrAux, waveHndl& wavH, waveHndl wavH_Wfr);

	static int ShowCompProgress(double CurNormVal);
    static int InitCompProgressIndicator();
    static int UpdateCompProgressIndicator(double CurNormVal);
    static void DestroyCompProgressIndicator();

	static int WfrModify(int ActionNo, srTSRWRadInData* pNewRadStructAccessData, char PolarizComp);
	static int WfrCreateNew(srTSRWRadInData* pNewRadStructAccessData);
	static int WfrModifyNeNxNz(srTSRWRadInData* pNewRadStructAccessData, char PolarizComp);
	static int WfrDelete(srTSRWRadInData* pNewRadStructAccessData);
	static int WfrRename(srTSRWRadInData* pNewRadStructAccessData);
	static int WfrGetNames(srTSRWRadInData* pNewRadStructAccessData);

	//static int SetupExtractedWaveData(srTSRWRadInData* pSRWRadInData, int Int_or_Phase, int PlotType, int TransvPres, char* pExtrData, waveHndl wavH, int hStateExtractedData, srTIgorWaveAccessData* pExtrWaveData);
	static int SetupExtractedWaveData(srTSRWRadInData* pSRWRadInData, int Int_or_Phase, int PlotType, char* pExtrData, waveHndl wavH, int hStateExtractedData, srTIgorWaveAccessData* pExtrWaveData);
	static int FinishWorkingWithSRWRadStruct(srTSRWRadInData* pSRWRadInData);
	static int FinishWorkingWithSRWStokesStruct(srTSRWStokesInData* pStokesAccessData);
	static int FinishWorkingWithSRWPowDensStruct(srTSRWPowDensInData* pPowDensAccessData);
	static int FinishWorkingWithWave(srTIgorWaveAccessData* pExtrWaveData);
	static int FinishWorkingWithTrjDataPointers(waveHndl wOutHorAng, waveHndl wOutHorCoor, waveHndl wOutVerAng, waveHndl wOutVerCoor, int hStateHorAng, int hStateHorCoor, int hStateVerAng, int hStateVerCoor);
    static int FinishWorkingWithTrj3dDataPointers(waveHndl wOutHorAng, waveHndl wOutHorCoor, waveHndl wOutLongDer, waveHndl wOutLongCoor, waveHndl wOutVerAng, waveHndl wOutVerCoor, int hStateHorAng, int hStateHorCoor, int hStateLongDer, int hStateLongCoor, int hStateVerAng, int hStateVerCoor);

private:

	static int UpdateNumberPositionInSRWRad(srTSRWRadInData* pSRWRadStructAccessData, long* RadIndices, char* CharBuf, int AmOfBytes);
	static int UpdateSrwWfrAuxData(srTSRWRadInData* pSRWRadInData);
	static int SetupSrwWfrAuxData(srTSRWRadInData* pSRWRadInData);
	static int UpdateTextPositionInSRWRad(srTSRWRadInData* pRad, int Index, char* Text);

	static int UtiCopyStringVect2ppChar(vector<string>& StrVect, char**& StrArr);
	static int UtiDelStrArr(char**& StringArr, int LenStringArr);
};

//*************************************************************************

#endif