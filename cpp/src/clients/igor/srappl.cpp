/************************************************************************//**
 * File: srappl.cpp
 * Description: Application class (obsolete)
 * Project: Synchrotron Radiation Workshop (SRW)
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srappl.h"
#include "gmfft.h"
#include "sroptmat.h"
#include "srmatsta.h"
#include "srmemory.h"
#include "srisosrc.h"
#include "srgsnbm.h"
#include "srsase.h"
#include "srremflp.h"
#include "srclcuti.h"
#include "srpropme.h"
#include "sremitpr.h"
#include "sroptdrf.h"
#include "srinterf.h"
#include "gmercode.h"

#ifdef __IGOR_PRO__
#include "srigsend.h"
#endif

//*************************************************************************

int srTApplication::Make1DFFT(srTHandleOfOneWaveStruct* pInputOneWaveStruct)
{
	int result;
	srTWaveAccessData WaveAccessData;
	if(result = Send.GetWaveAccessData(pInputOneWaveStruct, &WaveAccessData)) return result;

	if(WaveAccessData.AmOfDims != 1) 
	{
		Send.FinishWorkingWithWave(&WaveAccessData); return NEEDS_1D_WAVE;
	}
	if((WaveAccessData.WaveType[0] != 'c') || (WaveAccessData.WaveType[1] != 'f')) 
	{
		Send.FinishWorkingWithWave(&WaveAccessData); return NT_FP32_COMPLEX_WAVE_REQUIRED;
	}
	long Nx = *(WaveAccessData.DimSizes);
	if(((Nx>>1)<<1) != Nx) 
	{
		Send.FinishWorkingWithWave(&WaveAccessData); return WAVE_SIZES_SHOULD_BE_EVEN;
	}

	float* AuxDataCont = new float[Nx << 1];
	if(AuxDataCont == 0) return MEMORY_ALLOCATION_FAILURE;

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.pInData = (float*)(WaveAccessData.pWaveData);
	FFT1DInfo.pOutData = AuxDataCont;
	FFT1DInfo.Dir = char(pInputOneWaveStruct->p1);
	FFT1DInfo.xStep = WaveAccessData.DimSteps[0];
	FFT1DInfo.xStart = WaveAccessData.DimStartValues[0];
	FFT1DInfo.Nx = Nx;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;

	CGenMathFFT1D FFT1D;
	FFT1D.Make1DFFT(FFT1DInfo);

	WaveAccessData.DimStartValues[0] = FFT1DInfo.xStartTr;
	WaveAccessData.DimSteps[0] = FFT1DInfo.xStepTr;

	float *tOut = FFT1DInfo.pInData, *t = AuxDataCont;
	for(int ix=0; ix<(Nx<<1); ix++) *(tOut++) = *(t++);

	Send.FinishWorkingWithWave(&WaveAccessData);

	if(AuxDataCont != 0) delete[] AuxDataCont;
	return 0;
}
//*************************************************************************

int srTApplication::Make2DFFT(srTHandleOfOneWaveStruct* pInputOneWaveStruct)
{
	int result;
	srTWaveAccessData WaveAccessData;
	if(result = Send.GetWaveAccessData(pInputOneWaveStruct, &WaveAccessData)) return result;

	if(WaveAccessData.AmOfDims != 2) 
	{
		Send.Finish2DFFT(&WaveAccessData); return NEEDS_2D_WAVE;
	}
	if((WaveAccessData.WaveType[0] != 'c') || (WaveAccessData.WaveType[1] != 'f')) 
	{
		Send.Finish2DFFT(&WaveAccessData); return NT_FP32_COMPLEX_WAVE_REQUIRED;
	}

	long Nx = *(WaveAccessData.DimSizes);
	long Ny = *(WaveAccessData.DimSizes + 1);
	if((((Nx>>1)<<1) != Nx) || (((Ny>>1)<<1) != Ny)) 
	{
		Send.Finish2DFFT(&WaveAccessData); return WAVE_SIZES_SHOULD_BE_EVEN;
	}

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.pData = (float*)(WaveAccessData.pWaveData);
	FFT2DInfo.Dir = char(pInputOneWaveStruct->p1);
	FFT2DInfo.xStep = WaveAccessData.DimSteps[0];
	FFT2DInfo.yStep = WaveAccessData.DimSteps[1];
	FFT2DInfo.xStart = WaveAccessData.DimStartValues[0];
	FFT2DInfo.yStart = WaveAccessData.DimStartValues[1];
	FFT2DInfo.Nx = Nx;
	FFT2DInfo.Ny = Ny;
	FFT2DInfo.UseGivenStartTrValues = 0;

	CGenMathFFT2D FFT2D;
	FFT2D.Make2DFFT(FFT2DInfo);

	WaveAccessData.DimStartValues[0] = FFT2DInfo.xStartTr;
	WaveAccessData.DimStartValues[1] = FFT2DInfo.yStartTr;
	WaveAccessData.DimSteps[0] = FFT2DInfo.xStepTr;
	WaveAccessData.DimSteps[1] = FFT2DInfo.yStepTr;

	Send.Finish2DFFT(&WaveAccessData);
	return 0;
}

//*************************************************************************

int srTApplication::FindAverageDistanceToSource(double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc)
{
	return FindAverageDistanceToSource(TrjDat, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc);
}

//*************************************************************************

int srTApplication::FindAverageDistanceToSource(srTTrjDat& TrjDat, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc)
{// Should be called after the trajectory is already computed!!!
	double sRange = (TrjDat.LenFieldData - 1)*TrjDat.sStep;
	double sEnd = TrjDat.sStart + sRange;
	double Ysrc = 0.;

	char DistUnits = 0; // X, Z in mm
	double *TmpDataStorage = new double[TrjDat.LenFieldData << 2];
	if(TmpDataStorage == 0) return MEMORY_ALLOCATION_FAILURE;
	double *BtxArr = TmpDataStorage;
	double *BtzArr = TmpDataStorage + TrjDat.LenFieldData;
	double *xArr = TmpDataStorage + (TrjDat.LenFieldData << 1);
	double *zArr = TmpDataStorage + (TrjDat.LenFieldData*3);
	TrjDat.CompTotalTrjDataTrjDisp(TrjDat.sStart, sEnd, TrjDat.LenFieldData, BtxArr, BtzArr, xArr, zArr, DistUnits);

	int Len_mi_1 = TrjDat.LenFieldData - 1;
	double *pBtx = BtxArr + Len_mi_1, *pBtz = BtzArr + Len_mi_1, *pX = xArr + Len_mi_1, *pZ = zArr + Len_mi_1;
	double RobsLoc = RadInt.DistrInfoDat.yStart - sEnd;
	double InvRobsLoc = 1./RobsLoc;
	double MisFitXStLast = (RadInt.DistrInfoDat.xStart - *pX)*InvRobsLoc - *pBtx;
	double MisFitXFiLast = (RadInt.DistrInfoDat.xEnd - *pX)*InvRobsLoc - *pBtx;
	double MisFitZStLast = (RadInt.DistrInfoDat.zStart - *pZ)*InvRobsLoc - *pBtz;
	double MisFitZFiLast = (RadInt.DistrInfoDat.zEnd - *pZ)*InvRobsLoc - *pBtz;

	const double VeryLarge = 1.E+23;
	double RobsXSt = VeryLarge, RobsXFi = VeryLarge, RobsZSt = VeryLarge, RobsZFi = VeryLarge;

	for(int is=1; is<TrjDat.LenFieldData; is++)
	{
		RobsLoc += TrjDat.sStep;
		InvRobsLoc = 1./RobsLoc;
		pBtx--; pBtz--; pX--; pZ--;

		if(RobsXSt == VeryLarge)
		{
			double MisFitXSt = (RadInt.DistrInfoDat.xStart - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXSt*MisFitXStLast < 0.) RobsXSt = RobsLoc;
		}
		if(RobsXFi == VeryLarge)
		{
			double MisFitXFi = (RadInt.DistrInfoDat.xEnd - *pX)*InvRobsLoc - *pBtx;
			if(MisFitXFi*MisFitXFiLast < 0.) RobsXFi = RobsLoc;
		}
		if(RobsZSt == VeryLarge)
		{
			double MisFitZSt = (RadInt.DistrInfoDat.zStart - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZSt*MisFitZStLast < 0.) RobsZSt = RobsLoc;
		}
		if(RobsZFi == VeryLarge)
		{
			double MisFitZFi = (RadInt.DistrInfoDat.zEnd - *pZ)*InvRobsLoc - *pBtz;
			if(MisFitZFi*MisFitZFiLast < 0.) RobsZFi = RobsLoc;
		}
	}
	double MinRobsX = (RobsXSt < RobsXFi)? RobsXSt : RobsXFi;
	double MinRobsZ = (RobsZSt < RobsZFi)? RobsZSt : RobsZFi;
	double MinRobs = (MinRobsX < MinRobsZ)? MinRobsX : MinRobsZ;
	if(MinRobs != VeryLarge)
	{
		Robs = MinRobs;
		RobsAbsErr = 0.25*sRange;
	}
	else
	{
		if((TrjDat.sStart < 0.) && (sEnd > 0.)) Ysrc = 0.35*sRange;
		else Ysrc = TrjDat.sStart + 0.75*sRange;

		Robs = RadInt.DistrInfoDat.yStart - Ysrc;
		RobsAbsErr = 0.25*sRange;
	}

	// Make more smart estimation of the "Source Point" later !!!
	xElAtYsrc = TrjDat.EbmDat.x0; zElAtYsrc = TrjDat.EbmDat.z0;

	if(TmpDataStorage != 0) delete[] TmpDataStorage;
	return 0;
}

//*************************************************************************

int srTApplication::CheckNxNzForSR(double Robs, double xElAtYsrc, double zElAtYsrc, double MultFactor, srTWfrSmp& DistrInfoDat)
{
	long Nx = DistrInfoDat.nx, Nz = DistrInfoDat.nz;
	if((Nx > 0) && (Nz > 0)) return 0;

	const int SmallestN = 8;
	double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239854E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;

	CGenMathFFT2D FFT;
	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		double dx, dz;
		double HalfLambR = 0.5*WavelengthIn_m*Robs;
		if(Nx < 1)
		{
			double xStartRel = DistrInfoDat.xStart - xElAtYsrc;
			double xEndRel = DistrInfoDat.xEnd - xElAtYsrc;
			double dxStart = ::fabs(HalfLambR/xStartRel);
			double dxEnd = ::fabs(HalfLambR/xEndRel);
			dx = ((dxStart < dxEnd)? dxStart : dxEnd)/MultFactor;
			dx /= 1.2; // Remove this if Prop. Meth.2 does not work!

			Nx = long(::fabs(xEndRel - xStartRel)/dx) + 1;
			if(((Nx>>1)<<1) != Nx) Nx++;

			FFT.NextCorrectNumberForFFT(Nx);
			if(Nx < SmallestN) Nx = SmallestN;
		}
		if(Nz < 1)
		{
			double zStartRel = DistrInfoDat.zStart - zElAtYsrc;
			double zEndRel = DistrInfoDat.zEnd - zElAtYsrc;
			double dzStart = ::fabs(HalfLambR/zStartRel);
			double dzEnd = ::fabs(HalfLambR/zEndRel);
			dz = ((dzStart < dzEnd)? dzStart : dzEnd)/MultFactor;
			dz /= 1.2; // Remove this if Prop. Meth.2 does not work!

			Nz = long(::fabs(zEndRel - zStartRel)/dz) + 1;
			if(((Nz>>1)<<1) != Nz) Nz++;
			FFT.NextCorrectNumberForFFT(Nz);
			if(Nz < SmallestN) Nz = SmallestN;
		}
	}
// Continue here if(DistrInfoDat.CoordOrAngPresentation == AngPres)

	DistrInfoDat.nx = long(Nx); DistrInfoDat.nz = long(Nz);
	DistrInfoDat.DimensionsWereSetAuto = 1;
	return 0;
}

//*************************************************************************

void srTApplication::CheckAndCorrectNxNzForUnderSamplingMode(srTWfrSmp& DistrInfoDat)
{//Assumes the necessary sampling for propagation was already set
	//const int MaxAmOfPointsWithNoUnderSampling = 500; // To steer

// To Do it or not to do it ?

}

//*************************************************************************

int srTApplication::ProcessNxNzForPropag(srTWfrSmp& DistrInfoDat, srTSRWRadStructAccessData& SRWRadStructAccessData)
{// Assumes that Robs, xc, zc were already set in SRWRadStructAccessData
	int result;
	SRWRadStructAccessData.UnderSamplingX = SRWRadStructAccessData.UnderSamplingZ = 1;

	if(!DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat) return 0;
	
	DistrInfoDat.nx = DistrInfoDat.nz = -1;
	if(result = CheckNxNzForSR(SRWRadStructAccessData.RobsX, SRWRadStructAccessData.xc, SRWRadStructAccessData.zc, DistrInfoDat.NxNzOversamplingParam, DistrInfoDat)) return result;
	
	//CheckAndCorrectNxNzForUnderSamplingMode(DistrInfoDat);

	DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat = 0;

	if(DistrInfoDat.DimensionsWereSetAuto)
	{
		SRWRadStructAccessData.UpdateObsParam(DistrInfoDat);
		if(result = Send.ModifyRadNeNxNz(SRWRadStructAccessData)) return result;
		DistrInfoDat.DimensionsWereSetAuto = 0;
	}
	return 0;
}

//*************************************************************************

int srTApplication::ComputeSR_DirectOut()
{
	int result = 0, result1;
	srTGenOptElem GenOptElem;

	RadInt.pSend = &Send;

	if(result = Send.GetTotalElectronBeamDataFormat2(TrjDat)) return result;
	if(result = Send.GetTotalFieldDataFormat2(TrjDat)) return result;
	if(result = Send.GetAuxObsTreatParamFormat1(RadInt.DistrInfoDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadInt.DistrInfoDat)) return result;
	if(result = Send.GetTotalRadIntegrationParamDataFormat1(RadInt)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;
	SRWRadStructAccessData.SetRadSamplingFromObs(RadInt.DistrInfoDat);

	TrjDat.InputWasModified = 1;
	RadInt.TrjDataContShouldBeRebuild = 1;
	RadInt.DistrInfoDat.RadDistrDataContShouldBeRebuild = 0;
	double From_s0ToObsPoint = RadInt.DistrInfoDat.yStart - TrjDat.EbmDat.s0;

	if(result = TrjDat.ComputeInterpolatingStructure()) goto FinishWorkingWithSRWRad;
	SRWRadStructAccessData.InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	if(result = FindAverageDistanceToSource(Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) goto FinishWorkingWithSRWRad;
	SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = Robs;
	SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = RobsAbsErr;
	SRWRadStructAccessData.SetupXcZcFromElecData();
	
	if(result = ProcessNxNzForPropag(RadInt.DistrInfoDat, SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
	SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	if(result = RadInt.CheckInputConsistency()) goto FinishWorkingWithSRWRad;
	if(result = RadInt.ComputeTotalRadDistrDirectOut(SRWRadStructAccessData)) 
	{
		Send.FinishRadDistrOutFormat1();
		goto FinishWorkingWithSRWRad;
	}
	DirectRadDistrWasComputed = 1;

	if(result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;

FinishWorkingWithSRWRad:
	if(result1 = Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData)) return result1;
	Send.WarningMessages(&gVectWarnNos);

	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeWfrFromTrj()
{
	int result = 0;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTGenOptElem GenOptElem;
	double From_s0ToObsPoint;
	RadInt.pSend = &Send;

	char OldCompFromTrj = TrjDat.CompFromTrj;
	TrjDat.CompFromTrj = 1;

	if(result = Send.GetTotalTrajectoryDataFormat1(TrjDat)) return result; // also sets s0, x0, ...
	if(result = Send.GetAuxObsTreatParamFormat1(RadInt.DistrInfoDat)) goto FinishAndReleaseAll;
	if(result = Send.GetTotalObservationDataFormat3(RadInt.DistrInfoDat)) goto FinishAndReleaseAll;
	if(result = Send.GetTotalRadIntegrationParamDataFormat1(RadInt)) goto FinishAndReleaseAll;

	if(result = SRWRadStructAccessData.EmulateElectronBeamStruct(TrjDat.EbmDat)) goto FinishAndReleaseAll; // call before Send.GetSRWRadStructAccessData
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) goto FinishAndReleaseAll;
	SRWRadStructAccessData.SetRadSamplingFromObs(RadInt.DistrInfoDat);

	TrjDat.InputWasModified = 1;
	RadInt.TrjDataContShouldBeRebuild = 1;
	RadInt.DistrInfoDat.RadDistrDataContShouldBeRebuild = 0;

	if(result = TrjDat.ComputeInterpolatingStructure_FromTrj()) goto FinishAndReleaseAll;

	From_s0ToObsPoint = RadInt.DistrInfoDat.yStart - TrjDat.EbmDat.s0;
	SRWRadStructAccessData.InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	if(result = FindAverageDistanceToSource(Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) goto FinishAndReleaseAll;
	SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = Robs;
	SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = RobsAbsErr;
	SRWRadStructAccessData.xc = xElAtYsrc;
	SRWRadStructAccessData.zc = zElAtYsrc;

	if(result = ProcessNxNzForPropag(RadInt.DistrInfoDat, SRWRadStructAccessData)) goto FinishAndReleaseAll;
	SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	if(result = RadInt.CheckInputConsistency()) goto FinishAndReleaseAll;
	if(result = RadInt.ComputeTotalRadDistrDirectOut(SRWRadStructAccessData)) 
	{
		Send.FinishRadDistrOutFormat1(); goto FinishAndReleaseAll;
	}
	DirectRadDistrWasComputed = 1;

	if(result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData)) goto FinishAndReleaseAll;

FinishAndReleaseAll:
	Send.FinishWorkingWithTrajectoryData(TrjDat);
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);

	TrjDat.CompFromTrj = OldCompFromTrj;
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeWfrIsotrSrc()
{
	int result = 0;
	srTGenOptElem GenOptElem;
	srTIsotrSrc IsotrSrc;

	if(result = Send.GetTotalElectronBeamDataFormat3(IsotrSrc.EbmDat)) return result;
	if(result = Send.GetIsotrSrcExtraDataFormat1(IsotrSrc)) return result;
	if(result = Send.GetAuxObsTreatParamFormat1(IsotrSrc.DistrInfoDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(IsotrSrc.DistrInfoDat)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;
	SRWRadStructAccessData.SetRadSamplingFromObs(IsotrSrc.DistrInfoDat);

	double From_s0ToObsPoint = IsotrSrc.DistrInfoDat.yStart - IsotrSrc.EbmDat.s0;
	SRWRadStructAccessData.InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = From_s0ToObsPoint;
	SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = 0.001*From_s0ToObsPoint; // To steer
	SRWRadStructAccessData.xc = IsotrSrc.EbmDat.x0;
	SRWRadStructAccessData.zc = IsotrSrc.EbmDat.z0;

	if(result = ProcessNxNzForPropag(IsotrSrc.DistrInfoDat, SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
	SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	if(result = IsotrSrc.CreateWavefrontElField(SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
	if(result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;

FinishWorkingWithSRWRad:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeWfrGsnBeam()
{
	int result = 0;
	srTGenOptElem GenOptElem;
	srTGsnBeam GsnBeam;

	if(result = Send.GetTotalGsnBeamDataFormat1(GsnBeam)) return result;
	if(result = Send.GetAuxObsTreatParamFormat1(GsnBeam.DistrInfoDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(GsnBeam.DistrInfoDat)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;
	SRWRadStructAccessData.SetRadSamplingFromObs(GsnBeam.DistrInfoDat);

	double From_s0ToObsPoint = GsnBeam.DistrInfoDat.yStart - GsnBeam.EbmDat.s0;
	SRWRadStructAccessData.InitialSetupOf4x4PropMatr(From_s0ToObsPoint);

	SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = From_s0ToObsPoint;
	SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = 0.01*From_s0ToObsPoint; // To steer
	SRWRadStructAccessData.xc = GsnBeam.EbmDat.x0;
	SRWRadStructAccessData.zc = GsnBeam.EbmDat.z0;

	if(result = ProcessNxNzForPropag(GsnBeam.DistrInfoDat, SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
	SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	if(result = GsnBeam.CreateWavefrontElFieldFreqDomain(SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
	if(result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;

FinishWorkingWithSRWRad:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeWfrSASE()
{
	int result = 0;
	double From_s0ToObsPoint = 0.;
	//srTGenOptElem GenOptElem;
	//srTStringVect SRWRadStructVect;
	
	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTSRWRadStructAccessData *arSRWRadStructAccessData = 0;
	int numRadStruct = 0;
	//srTGenOptElem GenOptElem;
	srTSASE SASE;
	SASE.pSend = &Send;

	if(result = Send.GetTotalElectronBeamDataFormat3(SASE.EbmDat)) return result;
	//{
	//	if(Send.GetTotalElectronDistribDataFormat1(SASE.EbmDat)) return result;
	//}

	if(result = Send.GetGeneralMagneticFieldDataFormat1(SASE.MagDat)) return result;
	if(result = Send.GetSASEInputRadDataFormat1(SASE.InRad, SASE.SeedRad)) return result;

	SASE.DistrInfoDat.Initialize(); //Assuming that observation parameters may not be set
	if(result = Send.GetTotalObservationDataFormat3(SASE.DistrInfoDat)) return result;
	if(result = Send.GetSASEPrecDataFormat1(SASE.PrecDat)) return result;
	
	bool radStructArWasSetUp = false;
	if(Send.CheckIfRadStructWasSet())
	{
		if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;
		SRWRadStructAccessData.SetRadSamplingFromObs(SASE.DistrInfoDat);

		From_s0ToObsPoint = SASE.DistrInfoDat.yStart;
		SRWRadStructAccessData.InitialSetupOf4x4PropMatr(From_s0ToObsPoint);
		SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = From_s0ToObsPoint;
		SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = 0.05*From_s0ToObsPoint; // To steer
		SRWRadStructAccessData.xc = 0.; //SASE.EbmDat.x0;
		SRWRadStructAccessData.zc = 0.; //SASE.EbmDat.z0;
		if(result = ProcessNxNzForPropag(SASE.DistrInfoDat, SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
		SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

		if(result = SASE.CreateWavefrontElField(SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
		//if(result = GenOptElem.ComputeRadMoments(&SRWRadStructAccessData)) goto FinishWorkingWithSRWRad;
		//if(result = SRWRadStructAccessData.ComputeRadMoments()) goto FinishWorkingWithSRWRad;
		//moved to SASE.CreateWavefrontElField
	}
	else
	{
		if(result = Send.GetSRWRadStructArray(arSRWRadStructAccessData, numRadStruct)) return result;

		From_s0ToObsPoint = SASE.DistrInfoDat.yStart;
		srTSRWRadStructAccessData *t_arSRWRadStruct = arSRWRadStructAccessData;
		for(int i=0; i<numRadStruct; i++)
		{
			//srTSRWRadStructAccessData &curSRWRadStruct = arSRWRadStructAccessData[i];
			t_arSRWRadStruct->SetRadSamplingFromObs(SASE.DistrInfoDat);
			t_arSRWRadStruct->InitialSetupOf4x4PropMatr(From_s0ToObsPoint);
			t_arSRWRadStruct->RobsX = t_arSRWRadStruct->RobsZ = From_s0ToObsPoint;
			t_arSRWRadStruct->RobsXAbsErr = t_arSRWRadStruct->RobsZAbsErr = 0.05*From_s0ToObsPoint; // To steer
			t_arSRWRadStruct->xc = 0.; //SASE.EbmDat.x0;
			t_arSRWRadStruct->zc = 0.; //SASE.EbmDat.z0;
			if(result = ProcessNxNzForPropag(SASE.DistrInfoDat, *t_arSRWRadStruct)) goto FinishWorkingWithSRWRad;
			t_arSRWRadStruct->SetupNonZeroWavefrontLimitsAtCreation();
			t_arSRWRadStruct++;
		}

		if(result = SASE.CreateWavefrontElField(arSRWRadStructAccessData, numRadStruct)) goto FinishWorkingWithSRWRad;

		//t_arSRWRadStruct = arSRWRadStructAccessData;
		//for(int i=0; i<numRadStruct; i++) (t_arSRWRadStruct++)->ComputeRadMoments();
		//moved to SASE.CreateWavefrontElField

		radStructArWasSetUp = true;
	}

FinishWorkingWithSRWRad:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.FinishWorkingWithControlSASEStruct(SASE.ControlSASE);
	Send.FinishWorkingWithElecDistrStruct(SASE.EbmDat);
	
	if(radStructArWasSetUp) 
	{
		for(int i=0; i<numRadStruct; i++) Send.FinishWorkingWithSRWRadStruct(arSRWRadStructAccessData + i);
		delete[] arSRWRadStructAccessData;
	}

	Send.WarningMessages(&gVectWarnNos);
	return CorrectErrorCode(result);
}

//*************************************************************************
//Works only with IGOR!
int srTApplication::PropagateRadiation(srTIgorRadPropagInputStruct* pRadPropagInputStruct)
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTStringVect OptElemInfo;
	int result;
	if(result = Send.GetSRWRadStructAndOptElemNames(pRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return result;

	//OC140208
	//in some optical elements, StringArr[1] contains name of numerical wave describing optical element
	srTDataMD OptElemNumData;
	OptElemNumData.pData = 0;
	if(OptElemInfo[1] != 0)
	{
		if(*OptElemInfo[1] != '\0')
		{
#ifdef __IGOR_PRO__
			srTIgorSend::FetchNumWave(OptElemInfo[1], &OptElemNumData);
			//if fetch doesn't succeed, this is normal here
#endif
		}
	}

	srTGenOptElemHndl OptElemHndl;
	//srTOptElemSummary OptElemSummary;
	//result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	//result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, &OptElemNumData, &SRWRadStructAccessData, OptElemHndl);

	if(result != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return result;}

	if(result = PerformMethodDependentPrePropagationActions((int)(pRadPropagInputStruct->MethNo), (int)(pRadPropagInputStruct->AuxPar1), OptElemHndl, SRWRadStructAccessData)) return result;
	if(result = ((srTGenOptElem*)(OptElemHndl.rep))->CheckRadStructForPropagation(&SRWRadStructAccessData)) return result;

	srTRadResizeVect AuxResizeVect;
	
	char MethNo = (char)(pRadPropagInputStruct->MethNo);
	char UseResizeBefore = 0, UseResizeAfter = 0;
	if(MethNo > 0) { UseResizeBefore = UseResizeAfter = 1;}
	srTParPrecWfrPropag ParPrecWfrPropag(MethNo, UseResizeBefore, UseResizeAfter, 1., (double)(pRadPropagInputStruct->AuxPar1));

	//result = ((srTGenOptElem*)(OptElemHndl.rep))->PropagateRadiation(&SRWRadStructAccessData, (int)(pRadPropagInputStruct->MethNo), AuxResizeVect);
	result = ((srTGenOptElem*)(OptElemHndl.rep))->PropagateRadiation(&SRWRadStructAccessData, ParPrecWfrPropag, AuxResizeVect);

//Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);
	return CorrectErrorCode(result);
}

//*************************************************************************
//Works only with IGOR!
int srTApplication::PropagateWavefront(srTIgorWfrPropag* pWfrPropagInputStruct)
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTStringVect OptElemInfo;

	srTIgorRadPropagInputStruct AuxRadPropagInputStruct;
    AuxRadPropagInputStruct.wRad = pWfrPropagInputStruct->wRad;
    AuxRadPropagInputStruct.wOptElem = pWfrPropagInputStruct->wOptElem;

	int result;
	if(result = Send.GetSRWRadStructAndOptElemNames(&AuxRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return result;

	//OC140208
	//in some optical elements, StringArr[1] contains name of numerical wave describing optical element
	srTDataMD OptElemNumData;
	OptElemNumData.pData = 0;
	if(OptElemInfo.size() > 1)
	{
		if(OptElemInfo[1] != 0)
		{
			if(*OptElemInfo[1] != '\0')
			{
#ifdef __IGOR_PRO__
				srTIgorSend::FetchNumWave(OptElemInfo[1], &OptElemNumData);
				//if fetch doesn't succeed - this is normal here
#endif
			}
		}
	}

	srTGenOptElemHndl OptElemHndl;
	//result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, &OptElemNumData, &SRWRadStructAccessData, OptElemHndl);

	if(result != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return result;}

	srTHandleOfOneWaveStruct AuxOneWaveStruct;
	AuxOneWaveStruct.WaveHndl = pWfrPropagInputStruct->wPrecPar;
	srTWaveAccessData AuxWaveAccessData;
    if(result = Send.GetWaveAccessData(&AuxOneWaveStruct, &AuxWaveAccessData)) return result;
	
	//double AuxPrecParArr[4];
	//AuxWaveAccessData.OutRealData(AuxPrecParArr, 4);
	double AuxPrecParArr[5];
	AuxWaveAccessData.OutRealData(AuxPrecParArr, 5);
	Send.FinishWorkingWithWave(&AuxWaveAccessData);

	char UseResizeBefore = (char)AuxPrecParArr[0];
	char UseResizeAfter = (char)AuxPrecParArr[1];
	int MethNo = 0; 
	if(UseResizeBefore || UseResizeAfter) MethNo = 2;

	//srTParPrecWfrPropag ParPrecWfrPropag(MethNo, UseResizeBefore, UseResizeAfter, AuxPrecParArr[2], AuxPrecParArr[3]);
	//srTParPrecWfrPropag ParPrecWfrPropag(MethNo, UseResizeBefore, UseResizeAfter, AuxPrecParArr[2], AuxPrecParArr[3], AuxPrecParArr[4]);
	srTParPrecWfrPropag ParPrecWfrPropag(MethNo, UseResizeBefore, UseResizeAfter, AuxPrecParArr[2], AuxPrecParArr[3], (char)AuxPrecParArr[4]); //OC230810

	if(result = PerformMethodDependentPrePropagationActions(ParPrecWfrPropag.MethNo, ParPrecWfrPropag.UnderSampThresh, OptElemHndl, SRWRadStructAccessData)) return result;
	if(result = ((srTGenOptElem*)(OptElemHndl.rep))->CheckRadStructForPropagation(&SRWRadStructAccessData)) return result;

	if(result = SRWRadStructAccessData.EnsureFreqRepres()) return result;

	srTRadResizeVect AuxResizeVect;
	result = ((srTGenOptElem*)(OptElemHndl.rep))->PropagateRadiation(&SRWRadStructAccessData, ParPrecWfrPropag, AuxResizeVect);

//Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);

	return CorrectErrorCode(result);
}

//*************************************************************************
//obsolete
int srTApplication::WfrChangeRep(int MethNo)
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	srTSRWRadStructAccessData SRWRadStructAccessData;
	int res = 0;

	if(res = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return res;

	srTGenOptElemHndl OptElemHndl = srTGenOptElemHndl(new srTDriftSpace(1000.));
	res = ((srTGenOptElem*)(OptElemHndl.rep))->CheckRadStructForPropagation(&SRWRadStructAccessData);
	if(res != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return res;}
	res = ((srTGenOptElem*)(OptElemHndl.rep))->ChangeWfrRepres(&SRWRadStructAccessData, MethNo);

//Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	return CorrectErrorCode(res);
}

//*************************************************************************
//TEST
int srTApplication::PropagateRadiationStokesMultiElec(srTIgorRadPropagStokesMultiElecInputStruct* pRadPropagInputStruct)
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	int result = 0;
	//Send.SetIgorRadPropagStokesMultiElecInputStruct(pRadPropagInputStruct);

	srTEbmDat EbmDat;
	if(result = Send.GetTotalElectronBeamDataFormat3(EbmDat)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTStringVect OptElemInfo;
	if(result = Send.GetSRWRadStructAndOptElemNames(pRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return CorrectErrorCode(result);

	srTStokesStructAccessData StokesStructAccessData;
	double PrecPar[2];
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return CorrectErrorCode(result);
	if(result = Send.GetPropagRadStokesMultiElecDataFormat1(PrecPar)) return CorrectErrorCode(result);

	srTGenOptElemHndl OptElemHndl;
	//srTOptElemSummary OptElemSummary;
	//result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	//result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, 0, &SRWRadStructAccessData, OptElemHndl);

	if(result != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return CorrectErrorCode(result);}

	int MethNoSingleE = 2;
	if(result = PerformMethodDependentPrePropagationActions(MethNoSingleE, 0, OptElemHndl, SRWRadStructAccessData)) goto Termination;
	if(result = ((srTGenOptElem*)(OptElemHndl.rep))->CheckRadStructForPropagation(&SRWRadStructAccessData)) goto Termination;

	//result = srTPropagMultiE::PropagateElecFieldStokes(EbmDat, SRWRadStructAccessData, OptElemHndl, PrecPar, StokesStructAccessData);
	result = srTPropagMultiE::PropagateElecFieldStokesAuto(EbmDat, SRWRadStructAccessData, OptElemHndl, PrecPar, StokesStructAccessData);

Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);

	return CorrectErrorCode(result);
}

//*************************************************************************
//TEST
int srTApplication::EmitPropagRadStokesMultiElec(srTIgorRadEmitPropagStokesMultiElecInputStruct* pRadEmitPropag)
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	int result = 0;
	//Send.SetIgorRadPropagStokesMultiElecInputStruct(pRadPropagInputStruct);

	srTTrjDat LocTrjDat;
    srTStringVect OptElemInfo;
	srTGenOptElemHndl OptElemHndl;
    srTWfrSmp DistrInfoDat;
    double PrecPar[2];
	srTStokesStructAccessData StokesStructAccessData;

  	if(result = Send.GetTotalElectronBeamDataFormat2(LocTrjDat)) return CorrectErrorCode(result);
	if(result = Send.GetTotalFieldDataFormat2(LocTrjDat)) return CorrectErrorCode(result);
	if(result = Send.GetOptElemNames(OptElemInfo)) return CorrectErrorCode(result);
	//if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, NULL)) return CorrectErrorCode(result);
	if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, NULL, NULL, OptElemHndl)) return CorrectErrorCode(result);

	if(result = Send.GetTotalObservationDataFormat3(DistrInfoDat)) return CorrectErrorCode(result);
	if(result = Send.GetPropagRadStokesMultiElecDataFormat1(PrecPar)) return CorrectErrorCode(result);
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return CorrectErrorCode(result);

	if(result = srTPropagMultiE::ReallocateStokesAccordingToWfr(DistrInfoDat, StokesStructAccessData))
	{
		Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
		return CorrectErrorCode(result);
	}
	result = srTPropagMultiE::EmitPropagElecFieldStokes(LocTrjDat, OptElemHndl, PrecPar, StokesStructAccessData);

//Termination:
	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);
	return CorrectErrorCode(result);
}

//*************************************************************************
//TEST
int srTApplication::ComputeWfrEmitPropag()
{// In Propagation, all Lengths are in m, Photon Energy - in eV !
	int result = 0;
	
	srTTrjDat LocTrjDat;
	srTStringVect OptElemInfo;
	//srTOptElemSummary OptElemSummary;
	srTGenOptElemHndl OptElemHndl;
	srTWfrSmp DistrInfoDat;
	double PrecPar[2];
	srTSRWRadStructAccessData SRWRadStructAccessData;

  	if(result = Send.GetTotalElectronBeamDataFormat2(LocTrjDat)) return CorrectErrorCode(result);
	if(result = Send.GetTotalFieldDataFormat2(LocTrjDat)) return CorrectErrorCode(result);
	if(result = Send.GetOptElemNames(OptElemInfo)) return CorrectErrorCode(result);
	if(result = Send.GetTotalObservationDataFormat3(DistrInfoDat)) return CorrectErrorCode(result);
	if(result = Send.GetWfrEmitPropagPrec(PrecPar)) return CorrectErrorCode(result);
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return CorrectErrorCode(result);
	
	if(result = TrjDat.ComputeInterpolatingStructure()) goto Termination;
	if(result = FinalizeSetupOfSRWRadStruct(DistrInfoDat, LocTrjDat, SRWRadStructAccessData)) goto Termination;
	//if(result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData)) goto Termination;
	//if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData)) goto Termination;
	if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, 0, &SRWRadStructAccessData, OptElemHndl)) goto Termination;

	if(result = PerformMethodDependentPrePropagationActions(2, 0, OptElemHndl, SRWRadStructAccessData)) goto Termination;
	
	if(result = srTEmitPropag::ComputeRadiation(LocTrjDat, OptElemHndl, SRWRadStructAccessData, PrecPar)) goto Termination;
	
Termination:

	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.DeleteStringVect(OptElemInfo);
	Send.WarningMessages(&gVectWarnNos);

	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::FinalizeSetupOfSRWRadStruct(srTWfrSmp& DistrInfoDat, srTTrjDat& LocTrjDat, srTSRWRadStructAccessData& SRWRadStructAccessData)
{
	int result = 0;

	SRWRadStructAccessData.SetRadSamplingFromObs(DistrInfoDat);
	SRWRadStructAccessData.InitialSetupOf4x4PropMatr(DistrInfoDat.yStart - TrjDat.EbmDat.s0);

	double Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc;
	if(result = FindAverageDistanceToSource(LocTrjDat, Robs, RobsAbsErr, xElAtYsrc, zElAtYsrc)) return result;

	SRWRadStructAccessData.RobsX = SRWRadStructAccessData.RobsZ = Robs;
	SRWRadStructAccessData.RobsXAbsErr = SRWRadStructAccessData.RobsZAbsErr = RobsAbsErr;
	SRWRadStructAccessData.SetupXcZcFromElecData();

	if(result = ProcessNxNzForPropag(DistrInfoDat, SRWRadStructAccessData)) return result;
	SRWRadStructAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	return 0;
}

//*************************************************************************

int srTApplication::PerformMethodDependentPrePropagationActions(int MethNo, double AuxPar1, srTGenOptElemHndl& hOptElem, srTSRWRadStructAccessData& Rad)
{
	int result;

	//if(pRadPropagInputStruct->MethNo == 0) // Manual method
	if(MethNo == 0) // Manual method
	{
	}

	//if(pRadPropagInputStruct->MethNo == 1) // Obsolete Automatic method
	if(MethNo == 1) // Obsolete Automatic method
	{
		char SimpleOversamplingTestAndCorrectionShouldBeDone = 1; // To turn on/off manually
		if(SimpleOversamplingTestAndCorrectionShouldBeDone)
		{
			if(result = ((srTGenOptElem*)(hOptElem.rep))->MakeSimpleOversamplingTestAndCorrection(&Rad)) return result;
		}
	}

	//if(pRadPropagInputStruct->MethNo == 2) // Automatic method
	if(MethNo == 2) // Automatic method
	{
	}

	//if(pRadPropagInputStruct->MethNo == 3) // Automatic method with auto-switching to propag. with undersampling
	if(MethNo == 3) // Automatic method with auto-switching to propag. with undersampling
	{
		//pRadPropagInputStruct->MethNo = 2;

		Rad.AllowAutoSwitchToPropInUnderSamplingMode = 1;
		//Rad.InvUnderSamplingThreshold = pRadPropagInputStruct->AuxPar1;
		Rad.InvUnderSamplingThreshold = AuxPar1;
		Rad.EstimateAndSetUnderSampling();
	}
	return 0;
}

//*************************************************************************

int srTApplication::ComputePerStokes()
{
	srTIntVect Warnings;

	int result = 0;
	if(result = Send.GetTotalElectronBeamDataFormat3(RadIntPer.EbmDat)) return result;
	if(result = Send.GetPeriodicFieldDataFormat1(RadIntPer.MagPer)) return result;

	RadIntPer.DistrInfoDat.FluxComp = 1;
	if(result = Send.GetTotalObservationDataFormat3(RadIntPer.DistrInfoDat)) return result;
	if(result = Send.GetRadIntPeriodicParamDataFormat1(RadIntPer.IntPerStoPrec)) return result;

	RadIntPer.pWarningsGen = &Warnings;
	if(result = RadIntPer.CheckInputConsistency()) return result;

	srTStokesStructAccessData StokesStructAccessData;
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return result;

// WARNING: In this particular case the Photon Energy is treated internally in keV
// But this is an exception (normally, everywhere it is in eV)
	RadIntPer.DistrInfoDat.LambStart *= 0.001;
	RadIntPer.DistrInfoDat.LambEnd *= 0.001;

	//result = RadIntPer.ComputeTotalStokesDistr(StokesStructAccessData);
	result = RadIntPer.ComputeTotalStokesDistr(&StokesStructAccessData); //OC020112

	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&Warnings);

	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeStokesWig()
{
	int result = 0;
	srTGenTrjHndl hGenTrj;
	if(result = Send.GetGenMagFieldDataFormat1(hGenTrj)) return result;
	if(result = Send.GetTotalElectronBeamDataFormat3(hGenTrj.rep->EbmDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadIntWiggler.DistrInfoDat)) return result;
	if(result = Send.GetRadIntWigglerParamDataFormat1(RadIntWiggler.IntWigPrec)) return result;

	srTStokesStructAccessData StokesStructAccessData;
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return result;

	char LongIntTypeG = hGenTrj.rep->MagFieldPeriodicity();
	if(LongIntTypeG == 2) RadIntWiggler.CheckPossibilityOfFarFieldPeriodicComp(hGenTrj, LongIntTypeG);
	RadIntWiggler.LongIntTypeG = LongIntTypeG;

	srTGenTrjHndl LocTrjHndl = hGenTrj;
	if(result = LocTrjHndl.rep->ConvertToArbTrjDat(LongIntTypeG, RadIntWiggler.DistrInfoDat, hGenTrj)) return result;
	RadIntWiggler.TrjDatPtr = (srTTrjDat*)(hGenTrj.rep);

	result = RadIntWiggler.ComputeTotalStokesDistr(StokesStructAccessData);

	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&gVectWarnNos);

	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeStokesConst()
{
	int result = 0;
	srTGenTrjHndl hGenTrj;
	if(result = Send.GetGenMagFieldDataFormat1(hGenTrj)) return result;
	if(result = Send.GetTotalElectronBeamDataFormat3(hGenTrj.rep->EbmDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadIntConst.DistrInfoDat)) return result;
	if(result = Send.GetRadIntConstParamDataFormat1(RadIntConst.IntConstPrec)) return result;

	srTStokesStructAccessData StokesStructAccessData;
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return result;

	RadIntConst.ConstTrjDatPtr = (srTConstTrjDat*)(hGenTrj.rep);
	result = RadIntConst.ComputeTotalStokesDistr(StokesStructAccessData);

	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputeStokesArb()
{
	int result = 0;

	srTGenTrjHndl hGenTrj;
	if(result = Send.GetGenMagFieldDataFormat1(hGenTrj)) return result;
	if(result = Send.GetTotalElectronBeamDataFormat3(hGenTrj.rep->EbmDat)) return result;
	if(result = Send.GetTotalObservationDataFormat3(RadIntConst.DistrInfoDat)) return result;
/*
	if(result = Send.GetRadIntConstParamDataFormat1(RadIntConst.IntConstPrec)) return result;

	srTStokesStructAccessData StokesStructAccessData;
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return result;

	RadIntConst.ConstTrjDatPtr = (srTConstTrjDat*)(hGenTrj.rep);
	result = RadIntConst.ComputeTotalStokesDistr(StokesStructAccessData);

	Send.FinishWorkingWithStokesStruct(&StokesStructAccessData);
	Send.WarningMessages(&gVectWarnNos);
*/
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ComputePowerDensity()
{
	srTIntVect Warnings;
	int result = 0;
	srTGenTrjHndl hGenTrj;
	if(result = Send.GetGenMagFieldDataFormat1(hGenTrj)) return result;
	if(result = Send.GetTotalElectronBeamDataFormat3(hGenTrj.rep->EbmDat)) return result;
	RadIntPowDens.TrjHndl = hGenTrj;

	if(result = Send.GetTotalObservationDataFormat3(RadIntPowDens.DistrInfoDat)) return result;
	if(result = Send.GetRadIntPowDensParamDataFormat1(RadIntPowDens.IntPowDenPrec)) return result;

	srTPowDensStructAccessData PowDensStructAccessData;
	if(result = Send.GetPowDensStructAccessData(&PowDensStructAccessData)) return result;
	
	RadIntPowDens.pWarningsGen = &Warnings;
	result = RadIntPowDens.ComputeTotalPowerDensityDistr(PowDensStructAccessData);

	Send.FinishWorkingWithPowDensStruct(&PowDensStructAccessData);
	Send.WarningMessages(&Warnings);
	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::ExtractPhase(srTIgorExtractPhaseInputStruct* pInputStruct)
{
	int result = 0;
	srTHandleOfOneWaveStruct InputOneWaveStruct;

	InputOneWaveStruct.WaveHndl = pInputStruct->InWaveHndl;
	srTWaveAccessData InWaveAccessData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &InWaveAccessData)) return result;
	if((InWaveAccessData.WaveType[0] != 'c') || (InWaveAccessData.WaveType[1] != 'f')) 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); return NT_FP32_COMPLEX_WAVE_REQUIRED;
	}
	if(InWaveAccessData.AmOfDims > 2) 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); return NEEDS_1D_OR_2D_WAVE;
	}

	InputOneWaveStruct.WaveHndl = pInputStruct->Out1WaveHndl;
	srTWaveAccessData Out1WaveAccessData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &Out1WaveAccessData)) return result;
	if(Out1WaveAccessData.WaveType[0] != 'd') 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); 
		Send.FinishWorkingWithWave(&Out1WaveAccessData); 
		return NT_FP64_WAVE_REQUIRED;
	}
	if(Out1WaveAccessData.AmOfDims > 2) 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); 
		Send.FinishWorkingWithWave(&Out1WaveAccessData); 
		return NEEDS_1D_OR_2D_WAVE;
	}

	InputOneWaveStruct.WaveHndl = pInputStruct->Out2WaveHndl;
	srTWaveAccessData Out2WaveAccessData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &Out2WaveAccessData)) return result;
	if(Out2WaveAccessData.WaveType[0] != 'd')
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); 
		Send.FinishWorkingWithWave(&Out1WaveAccessData); 
		Send.FinishWorkingWithWave(&Out2WaveAccessData); 
		return NT_FP64_WAVE_REQUIRED;
	}
	if(Out2WaveAccessData.AmOfDims > 2) 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); 
		Send.FinishWorkingWithWave(&Out1WaveAccessData); 
		Send.FinishWorkingWithWave(&Out2WaveAccessData); 
		return NEEDS_1D_OR_2D_WAVE;
	}

	if((InWaveAccessData.AmOfDims != Out1WaveAccessData.AmOfDims) || (InWaveAccessData.AmOfDims != Out2WaveAccessData.AmOfDims)) 
	{
		Send.FinishWorkingWithWave(&InWaveAccessData); 
		Send.FinishWorkingWithWave(&Out1WaveAccessData); 
		Send.FinishWorkingWithWave(&Out2WaveAccessData); 
		return NEEDS_WAVE_OF_EQUEAL_SIZES;
	}
	for(int i=0; i<InWaveAccessData.AmOfDims; i++)
		if((InWaveAccessData.DimSizes[i] != Out1WaveAccessData.DimSizes[i]) || (InWaveAccessData.DimSizes[i] != Out2WaveAccessData.DimSizes[i]))
		{
			Send.FinishWorkingWithWave(&InWaveAccessData); 
			Send.FinishWorkingWithWave(&Out1WaveAccessData); 
			Send.FinishWorkingWithWave(&Out2WaveAccessData); 
			return NEEDS_WAVE_OF_EQUEAL_SIZES;
		}

	double* pOutMag = (double*)(Out1WaveAccessData.pWaveData);
	double* pOutPhase = (double*)(Out2WaveAccessData.pWaveData);
	int xToRemove = 0, zToRemove = 0;
	srTGenOptElem GenOptElem;
	result = GenOptElem.GenExtractPhase(InWaveAccessData, pOutMag, pOutPhase, xToRemove, zToRemove);

	Send.FinishWorkingWithWave(&InWaveAccessData);
	Send.FinishWorkingWithWave(&Out1WaveAccessData);
	Send.FinishWorkingWithWave(&Out2WaveAccessData);
	return result;
}

//*************************************************************************
/**
int srTApplication::OptGenTransSetup(srTIgorOptGenTransSetupInputStruct* pInputStruct)
{//Obsolette
	int result = 0;
	srTStringVect OptElemInfo;
	if(result = Send.GetVectorOfStrings(pInputStruct->wOptElem, &OptElemInfo)) return result;
	
	srTSRWRadStructAccessData DummyRadAccessData;
	srTGenOptElemHndl OptElemHndl;
	//srTOptElemSummary OptElemSummary;
	//if(result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &DummyRadAccessData)) return result;
	//if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &DummyRadAccessData)) return result;
	if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, 0, &DummyRadAccessData, OptElemHndl)) return result;

	result = Send.UpdateTextWave(pInputStruct->wOptElem, &OptElemInfo);
	Send.WarningMessages(&gVectWarnNos);
	Send.DeleteStringVect(OptElemInfo);
	return result;
}
**/
//*************************************************************************

int srTApplication::ReflectWavefront()
{// To implement: reflection from multiple infinite planes
	int result = 0;

	srTEbmDat EbmDat;
	srTSRWRadStructAccessData SRWRadStructAccessData;
	if(result = Send.GetTotalElectronBeamDataFormat3(EbmDat)) return result;
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;


//read in mirrors - opt. components
//ddddddddddddddddddddddddddddd

/**

	srTStringVect OptElemInfo;
	if(result = Send.GetSRWRadStructAndOptElemNames(pRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return CorrectErrorCode(result);

	srTStokesStructAccessData StokesStructAccessData;
	double PrecPar[2];
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return CorrectErrorCode(result);
	if(result = Send.GetPropagRadStokesMultiElecDataFormat1(PrecPar)) return CorrectErrorCode(result);

	srTGenOptElemHndl OptElemHndl;
	srTOptElemSummary OptElemSummary;
	result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData);
	if(result != 0) { Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData); return CorrectErrorCode(result);}

	int MethNoSingleE = 2;
	if(result = PerformMethodDependentPrePropagationActions(MethNoSingleE, 0, OptElemHndl, SRWRadStructAccessData)) goto Termination;
	if(result = OptElemHndl.rep->CheckRadStructForPropagation(&SRWRadStructAccessData)) goto Termination;

	result = srTPropagMultiE::PropagateElecFieldStokes(EbmDat, SRWRadStructAccessData, OptElemHndl, PrecPar, StokesStructAccessData);
**/

//Termination:
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.WarningMessages(&gVectWarnNos);

	return CorrectErrorCode(result);
}

//*************************************************************************

int srTApplication::OptMatConst(srTIgorOptMatConstInputStruct* pInputStruct)
{
	int MatNo = int(pInputStruct->MatNo);
	int Z = int(pInputStruct->AtomNum);
	double Density = pInputStruct->Density;
	double PhotEn = pInputStruct->PhotEn;

	int result;
	srTOptMater OptMat;
	double RefrDelta;
	if(result = OptMat.CompRefractDelta(MatNo, Z, Density, PhotEn, RefrDelta)) return result;
	pInputStruct->Delta = RefrDelta;
	pInputStruct->AttenLength = OptMat.AttenLength(MatNo, PhotEn);

	Send.WarningMessages(&gVectWarnNos);
	return 0;
}

//*************************************************************************

int srTApplication::UtiSpotInfo(srTIgorUtiSpotInfoInputStruct* pInputStruct)
{
	int result;
	srTHandleOfOneWaveStruct InputOneWaveStruct;
	srTWaveAccessData InWaveData;

	InputOneWaveStruct.WaveHndl = pInputStruct->wData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &InWaveData)) return result;

	srTWaveAccessData OutSpotInfoData;

	char InfoSuf[] = "_inf"; // To steer

	(OutSpotInfoData.WaveType)[0] = 'f';
	OutSpotInfoData.AmOfDims = 1;
	(OutSpotInfoData.DimSizes)[0] = 5;
	(OutSpotInfoData.DimSizes)[1] = 0;
	(OutSpotInfoData.DimStartValues)[0] = 0;
	(OutSpotInfoData.DimSteps)[0] = 1;

	strcpy(OutSpotInfoData.NameOfWave, InWaveData.NameOfWave);
	strcat(OutSpotInfoData.NameOfWave, InfoSuf);

	if(result = Send.MakeWaveAccordingToWaveAccessData(OutSpotInfoData)) return result;

	InputOneWaveStruct.WaveHndl = OutSpotInfoData.wHndl;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &OutSpotInfoData)) return result;

	long LenSpotInfoWaveName = (long)strlen(OutSpotInfoData.NameOfWave);

	Handle SpotInfoWaveNameStr = NewHandle(LenSpotInfoWaveName);

	strncpy(*SpotInfoWaveNameStr, OutSpotInfoData.NameOfWave, LenSpotInfoWaveName);
	pInputStruct->SpotInfoName = SpotInfoWaveNameStr;

	srTAuxMatStat AuxMatStat;
	if(result = AuxMatStat.FindSimplestStat(InWaveData, OutSpotInfoData))
	{
		Send.DeleteWave(OutSpotInfoData);
		Send.FinishWorkingWithWave(&InWaveData);
		return result;
	}

	Send.FinishWorkingWithWave(&InWaveData);
	Send.FinishWorkingWithWave(&OutSpotInfoData);
	return 0;
}

//*************************************************************************

int srTApplication::UtiWfrLimits(srTIgorUtiWrfLimitsInputStruct* pInputStruct)
{
	int result;
	srTHandleOfOneWaveStruct InputOneWaveStruct;
	srTWaveAccessData InWaveData;

	InputOneWaveStruct.WaveHndl = pInputStruct->wData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &InWaveData)) return result;

	srTWaveAccessData OutInfoData;
	char InfoSuf[] = "_inf"; // To steer

	(OutInfoData.WaveType)[0] = 'f';
	OutInfoData.AmOfDims = 1;
	(OutInfoData.DimSizes)[0] = 5;
	(OutInfoData.DimSizes)[1] = 0;
	(OutInfoData.DimStartValues)[0] = 0;
	(OutInfoData.DimSteps)[0] = 1;

	strcpy(OutInfoData.NameOfWave, InWaveData.NameOfWave);
	strcat(OutInfoData.NameOfWave, InfoSuf);
	if(result = Send.MakeWaveAccordingToWaveAccessData(OutInfoData)) return result;
	InputOneWaveStruct.WaveHndl = OutInfoData.wHndl;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &OutInfoData)) return result;
	long LenInfoWaveName = (long)strlen(OutInfoData.NameOfWave);
	Handle InfoWaveNameStr = NewHandle(LenInfoWaveName);

	strncpy(*InfoWaveNameStr, OutInfoData.NameOfWave, LenInfoWaveName);
	pInputStruct->InfoName = InfoWaveNameStr;

	srTAuxMatStat AuxMatStat;

	if(result = AuxMatStat.FindIntensityLimits(InWaveData, (double)(pInputStruct->PowLevel), OutInfoData))
	{
		Send.DeleteWave(OutInfoData);
		Send.FinishWorkingWithWave(&InWaveData);
		return result;
	}

	Send.FinishWorkingWithWave(&InWaveData);
	Send.FinishWorkingWithWave(&OutInfoData);
	return 0;
}

//*************************************************************************

int srTApplication::UtiRemoveFlips(srTIgorUtiRemoveFlipsInputStruct* pInputStruct)
{
	int result;
	srTHandleOfOneWaveStruct InputOneWaveStruct;
	srTWaveAccessData InWaveData;
	InputOneWaveStruct.WaveHndl = pInputStruct->wData;
	if(result = Send.GetWaveAccessData(&InputOneWaveStruct, &InWaveData)) return result;

	srTAuxRemoveFlips RemoveFlips;
	if(result = RemoveFlips.GenRemoveFlips(InWaveData)) return result;

	Send.FinishWorkingWithWave(&InWaveData);
	return 0;
}

//*************************************************************************

int srTApplication::OptZonePlateSpecSetup(srTIgorOptZonePlateSpecSetupInputStruct* pInputStruct)
{
	int result;
	srTStringVect OptElemInfo;
	if(result = Send.GetVectorOfStrings(pInputStruct->wOptElem, &OptElemInfo)) return result;

	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTHandleOfSRWRadStruct HandleOfSRWRadStruct;
	HandleOfSRWRadStruct.wRad = pInputStruct->wRad;
	Send.SetHandleOfSRWRadStruct(&HandleOfSRWRadStruct);
	if(result = Send.GetSRWRadStructAccessData(&SRWRadStructAccessData)) return result;

	srTGenOptElemHndl OptElemHndl;
	//srTOptElemSummary OptElemSummary;
	//if(result = OptElemSummary.SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData)) return result;
	//if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, OptElemHndl, &SRWRadStructAccessData)) return result;
	if(result = srTGenOptElem::SetupOpticalElement(&OptElemInfo, 0, &SRWRadStructAccessData, OptElemHndl)) return result;

	result = Send.UpdateTextWave(pInputStruct->wOptElem, &OptElemInfo);
	Send.FinishWorkingWithSRWRadStruct(&SRWRadStructAccessData);
	Send.DeleteStringVect(OptElemInfo);

	return result;
}

//*************************************************************************

double srTApplication::UtiMagRad(double Bconst, double ElecEnergy, double OutUnit)
{
	return srTCalcUtils::MagnetRadius(Bconst, ElecEnergy, (int)OutUnit);
}

//*************************************************************************

double srTApplication::UtiMagCritPhotEn(double Bconst, double ElecEnergy, double OutUnit)
{
	return srTCalcUtils::MagnetCritPhotEn(Bconst, ElecEnergy, (int)OutUnit);
}

//*************************************************************************

double srTApplication::UtiUndK(double Bpeak, double Period)
{
	return srTCalcUtils::UndK(Bpeak, Period);
}

//*************************************************************************

double srTApplication::UtiUndFundPhotEn(double Bpeak, double Period, double ElecEnergy, double OutUnit)
{
	return srTCalcUtils::UndFundPhotEn(Bpeak, Period, ElecEnergy, (int)OutUnit);
}

//*************************************************************************

int srTApplication::CorrectErrorCode(int ErrCode)
{
	if(ErrCode == 1) return SR_COMP_PROC_ABORTED; // 1 is returned when user presses Abort
	else if(ErrCode == GM_RK_MAX_NUM_STEPS_REACHED) return SRW_GM_RK_MAX_NUM_STEPS_REACHED;
	else if(ErrCode == GM_RK_STEP_SIZE_TOO_SMALL) return SRW_GM_RK_STEP_SIZE_TOO_SMALL;
	else return ErrCode;
}

//*************************************************************************
//TEST
int srTApplication::Uti3DView(srTEbmDat& ElBeam, srTMagElem& MagField)
{
	int result = 0;
/*
	Send.SetIgor3DViewInputStruct(p3DViewInputStruct);

	srTEbmDat EbmDat;
	if(result = Send.GetTotalElectronBeamDataFormat3(EbmDat)) return result;
*/

/*	

	srTSRWRadStructAccessData SRWRadStructAccessData;
	srTStringVect OptElemInfo;
	if(result = Send.GetSRWRadStructAndOptElemNames(pRadPropagInputStruct, &SRWRadStructAccessData, &OptElemInfo)) return CorrectErrorCode(result);

	srTStokesStructAccessData StokesStructAccessData;
	double PrecPar[2];
	if(result = Send.GetStokesStructAccessData(&StokesStructAccessData)) return CorrectErrorCode(result);
*/

	return result;
}

//*************************************************************************

