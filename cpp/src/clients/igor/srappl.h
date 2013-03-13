/************************************************************************//**
 * File: srappl.h
 * Description: Application class (obsolete; still used for IGOR Pro SRW version) header
 * Project: Synchrotron Radiation Workshop (SRW)
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRAPPL_H
#define __SRAPPL_H

#include "sroptelm.h"
#include "srtrjdat.h"
#include "srradint.h"
//#include "sroptics.h" //OC221012
#include "srstraux.h"
#include "srradmnp.h"
#include "srebmdat.h"
#include "srpersto.h"
#include "srstowig.h"
#include "srradinc.h"
#include "srpowden.h"
#include "srobject.h"
#include "srsend.h"
#include "srigintr.h"
#include "srigintr.h"
#include "gmfunc.h"

#include <time.h>

//*************************************************************************
//
// Units at the input and within the code:
//
// All Cartesian coordinates and Distances - m
// All Angles (incl. dxds, etc.) - r
// El. Current - A
// El. Energy - GeV
// Photon Energy - eV (Exception: at UR Spectrum through a slit computation internal units are keV)
// Wavelength = nm
// 
//-------------------------------------------------------------------------
//
// Make units transformation in accordance with version for Igor when preparing a version for Mathematica !!!
//
//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

class srTApplication {

	CMHGenObj ObjHandleMap;

	short SendingIsRequired;
	char DirectRadDistrWasComputed, FocusRadDistrWasComputed;

	clock_t ClockTime;

public:

	srTSend Send;
	srTTrjDat TrjDat;

	srTRadInt RadInt;
	srTRadIntPeriodic RadIntPer;
	srTRadIntWiggler RadIntWiggler;
	srTRadIntConst RadIntConst;
	srTRadIntPowerDensity RadIntPowDens;

	//srTOpticsHndl OpticsHndl; //OC221012

	srTApplication() 
	{
		Initialize();
	}
	~srTApplication()
	{
		RadInt.DeallocateMemForRadDistr();
	}

	void Initialize()
	{
		RadInt.Initialize();
		TrjDat.Initialize();
		RadIntPer.Initialize();

		SendingIsRequired = 1;
		RadInt.SetInputTrjData(&TrjDat);
		DirectRadDistrWasComputed = FocusRadDistrWasComputed = 0;
		ClockTime = clock();
	}

	inline int ComputeTrajectory();
	int ComputeSR_DirectOut();
	int ComputeWfrFromTrj();
	int ComputeWfrIsotrSrc();
	int ComputeWfrGsnBeam();
	int ComputeWfrSASE();
	int ComputePerStokes();
	int ComputeStokesWig();
	int ComputeStokesConst();
	int ComputeStokesArb();
	int ComputePowerDensity();
	inline int ComputeMomentsOfSR();

	int Make1DFFT(srTHandleOfOneWaveStruct*);
	int Make2DFFT(srTHandleOfOneWaveStruct*);

	int FindAverageDistanceToSource(srTTrjDat&, double&, double&, double&, double&);
	int FindAverageDistanceToSource(double&, double&, double&, double&);

	int CheckNxNzForSR(double, double, double, double, srTWfrSmp&);
	void CheckAndCorrectNxNzForUnderSamplingMode(srTWfrSmp&);

	inline int RadResize(srTIgorRadResizeInputStruct*);
	inline int ObsSetNxNz(srTIgorObsSetNxNzInputStruct*);
	int ProcessNxNzForPropag(srTWfrSmp&, srTSRWRadStructAccessData&);

	int PropagateRadiation(srTIgorRadPropagInputStruct*);
	int PropagateWavefront(srTIgorWfrPropag*);
	int PerformMethodDependentPrePropagationActions(int MethNo, double AuxPar, srTGenOptElemHndl& hOptElem, srTSRWRadStructAccessData& Rad);
	int FinalizeSetupOfSRWRadStruct(srTWfrSmp& DistrInfoDat, srTTrjDat& TrjDat, srTSRWRadStructAccessData& SRWRadStructAccessData);

	int WfrChangeRep(int MethNo);

	inline int ExtractRadiation(srTIgorRadExtractInputStruct*);
	//int OptGenTransSetup(srTIgorOptGenTransSetupInputStruct*); //obsolette; was generalized

	int OptMatConst(srTIgorOptMatConstInputStruct*);
	int UtiSpotInfo(srTIgorUtiSpotInfoInputStruct*);
	int UtiRemoveFlips(srTIgorUtiRemoveFlipsInputStruct*);
	int Uti3DView(srTEbmDat&, srTMagElem&);

	int UtiWfrLimits(srTIgorUtiWrfLimitsInputStruct*);

	//TEST
	int PropagateRadiationStokesMultiElec(srTIgorRadPropagStokesMultiElecInputStruct*);
	int EmitPropagRadStokesMultiElec(srTIgorRadEmitPropagStokesMultiElecInputStruct* pRadEmitPropag);
	//TEST
	//int ComputeWfrEmitPropag(srTIgorWfrEmitPropagInputStruct*);
	int ComputeWfrEmitPropag();
	//TEST
	int ReflectWavefront();

	double UtiMagRad(double Bconst, double ElecEnergy, double OutUnitStr);
	double UtiMagCritPhotEn(double Bconst, double ElecEnergy, double OutUnitStr);
	double UtiUndK(double Bpeak, double Period);
	double UtiUndFundPhotEn(double Bpeak, double ElecEnergy, double Period, double OutUnit);
	inline int Kmu(srTIgorKmuInputStruct*);

	int CorrectErrorCode(int ErrCode);

// Obsolete
	inline int ComputeSR();
	int OptZonePlateSpecSetup(srTIgorOptZonePlateSpecSetupInputStruct*);
	int ExtractPhase(srTIgorExtractPhaseInputStruct*);

};

//*************************************************************************

inline int srTApplication::ComputeTrajectory()
{
	TrjDat.InputWasModified = 1;
	int result;
	if(result = Send.GetTotalElectronBeamDataFormat2(TrjDat)) return result;
	if(result = Send.GetTotalFieldDataFormat2(TrjDat)) return result;
	if(result = TrjDat.ComputeInterpolatingStructure()) return result;
	if(result = Send.InitTrjOutFormat1(TrjDat)) return result;

	double sSt = TrjDat.sStart;
	double sFi = sSt + TrjDat.sStep*(TrjDat.LenFieldData - 1);
	char DistUnits = 1; // m for coord.
	DOUBLE *Loc_pOutBtxData, *Loc_pOutXData, *Loc_pOutBtzData, *Loc_pOutZData;
	Send.ShowOutTrjDataPointers(Loc_pOutBtxData, Loc_pOutXData, Loc_pOutBtzData, Loc_pOutZData);
	TrjDat.CompTotalTrjDataTrjDisp(sSt, sFi, TrjDat.LenFieldData, (double*)Loc_pOutBtxData, (double*)Loc_pOutBtzData, (double*)Loc_pOutXData, (double*)Loc_pOutZData, DistUnits);
	Send.FinishTrjOutFormat1();

	return 0;
}

//*************************************************************************

inline int srTApplication::ComputeSR()
{
	if(DirectRadDistrWasComputed) 
	{
		TrjDat.InputWasModified = 0;
		RadInt.DistrInfoDat.RadDistrDataContShouldBeRebuild = 0;
		RadInt.TrjDataContShouldBeRebuild = 0;
	}
	else
	{
		TrjDat.InputWasModified = 1;
		RadInt.DistrInfoDat.RadDistrDataContShouldBeRebuild = 1;
		RadInt.TrjDataContShouldBeRebuild = 1;
	}

	int result;
	if(result = Send.GetTotalElectronBeamDataFormat2(TrjDat)) return result;
	if(result = Send.GetTotalFieldDataFormat2(TrjDat)) return result;

	if(TrjDat.InputWasModified) if(result = TrjDat.ComputeInterpolatingStructure()) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadInt.DistrInfoDat)) return result;
	if(result = Send.GetTotalRadIntegrationParamDataFormat1(RadInt)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;

	RadInt.TrjDataContShouldBeRebuild = (RadInt.TrjDataContShouldBeRebuild || TrjDat.InputWasModified);
	RadInt.DeallocateMemForRadDistr();

	clock_t NewClockTime = clock();
	RadInt.ProbablyTheSameLoop = 0;
	if(!RadInt.TrjDataContShouldBeRebuild) 
	{
		RadInt.ProbablyTheSameLoop = (NewClockTime - ClockTime) < CLOCKS_PER_SEC;
		ClockTime = NewClockTime;
	}

	if(result = RadInt.ComputeTotalRadDistr()) return result;

	if(SendingIsRequired) if(result = Send.OutRadDistrFormat2(SRWRadStructAccessData, RadInt)) return result;
	DirectRadDistrWasComputed = 1;

	return 0;
}

//*************************************************************************

inline int srTApplication::RadResize(srTIgorRadResizeInputStruct* pRadResizeInputStruct)
{
	int result;
	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTRadResize RadResizeStruct;
	if(result = Send.GetSRWRadStructAndResizeData(pRadResizeInputStruct, &SRWRadStructAccessData, &RadResizeStruct)) return result;

	//double pxdIn = RadResizeStruct.pxd, pzdIn = RadResizeStruct.pzd;
	//char ThereWasUnderSamplingX = SRWRadStructAccessData.ThereIsUnderSamplingX();
	//char ThereWasUnderSamplingZ = SRWRadStructAccessData.ThereIsUnderSamplingZ();

	srTGenOptElem GenOptElem;
	if(result = GenOptElem.RadResizeGen(SRWRadStructAccessData, RadResizeStruct)) return result;

	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	return result;
}

//*************************************************************************

inline int srTApplication::ObsSetNxNz(srTIgorObsSetNxNzInputStruct* pObsSetNxNzInputStruct)
{
	srTIgorRadInputStruct IgorRadInputStruct;

	IgorRadInputStruct.wObservation = pObsSetNxNzInputStruct->wObservation;
	IgorRadInputStruct.wField = pObsSetNxNzInputStruct->wMagnField;
	IgorRadInputStruct.wElectronBeam = pObsSetNxNzInputStruct->wElectronBeam;

	srTIgorRadInputStruct* pOldIgorRadInputStruct = Send.GetIgorRadInputStructPtr();
	Send.SetIgorRadInputStruct(&IgorRadInputStruct);

	int result;
	if(result = Send.GetTotalElectronBeamDataFormat2(TrjDat)) return result;
	if(result = Send.GetTotalFieldDataFormat2(TrjDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadInt.DistrInfoDat)) return result;
	RadInt.DistrInfoDat.nx = RadInt.DistrInfoDat.nz = -1;

	double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	if(result = FindAverageDistanceToSource(Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) return result;

	if(result = CheckNxNzForSR(Robs, xElAtYsrc, zElAtYsrc, double(pObsSetNxNzInputStruct->p1), RadInt.DistrInfoDat)) return result;

	srTHandleOfOneWaveStruct HandleOfOneWaveStruct;
	HandleOfOneWaveStruct.WaveHndl = pObsSetNxNzInputStruct->wObservation;

	if(result = Send.ChangeObservationData(&HandleOfOneWaveStruct, RadInt.DistrInfoDat)) return result;

	Send.SetIgorRadInputStruct(pOldIgorRadInputStruct);
	return 0;
}

//*************************************************************************

inline int srTApplication::ComputeMomentsOfSR()
{// Here Lengths are in m !
	int result;
	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;

	srTGenOptElem GenOptElem;
	result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData);

	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	return result;
}

//*************************************************************************

inline int srTApplication::ExtractRadiation(srTIgorRadExtractInputStruct* pRadExtractInputStruct)
{// In Extraction, all Lengths are in m, Photon Energy - in eV !
	int result;
	//srTSRWRadStructAccessData SRWRadStructAccessData;
	srTSRWRadStructAccessData* pSRWRadStructAccessData = new srTSRWRadStructAccessData();
	srTSRWRadStructAccessData& SRWRadStructAccessData = *pSRWRadStructAccessData;

	srTRadExtract RadExtract;
	if(result = Send.GetSRWRadStructAndExtractData(pRadExtractInputStruct, &SRWRadStructAccessData, &RadExtract)) return result;

	//srTRadGenManip RadGenManip(SRWRadStructAccessData);
	CHGenObj hRad(pSRWRadStructAccessData);
	srTRadGenManip RadGenManip(hRad);

	srTWaveAccessData ExtractedWaveData;
	result = RadGenManip.ExtractRadiation(RadExtract, ExtractedWaveData);
	
	Send.FinishWorkingWithWave(&ExtractedWaveData);
	//Send.FinishWorkingWithSRWRadStruct(&(RadGenManip.RadAccessData));
	Send.FinishWorkingWithSRWRadStruct(pSRWRadStructAccessData);

	Send.WarningMessages(&gVectWarnNos);
	return result;
}

//*************************************************************************

inline int srTApplication::Kmu(srTIgorKmuInputStruct* p)
{
	//srTSpecialFunctions SpecFunc;
	//int result = srTMathFunctions::Kmu((int)(p->nInt), (double)(p->mu), (double)(p->x), p->Kmu);
	int result = CGenMathFunc::Kmu((int)(p->nInt), (double)(p->mu), (double)(p->x), p->Kmu);
	return result;
}

//*************************************************************************

#endif
