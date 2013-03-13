/************************************************************************//**
 * File: sroptzps.h
 * Description: Optical element: "Special" Zone Plate header (used to simulate "phase corrections") - Obsolete? 
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTZPS_H
#define __SROPTZPS_H

#include "sroptfoc.h"

//*************************************************************************

class srTZonePlateSpec : public srTFocusingElem {

	srTWaveAccessData ZonePlateNumData;
	char IsFocusing;

public:

	srTZonePlateSpec(srTStringVect* pElemInfo, srTSRWRadStructAccessData* pRad);
	~srTZonePlateSpec()
	{
		if(ZonePlateNumData.pWaveData != 0)
		{
			//srTSend Send; Send.FinishWorkingWithWave(&ZonePlateNumData);
			//???????????????
		}
	}

	int SetUpNumStruct(srTSRWRadStructAccessData* pRad, char* NumStructName, char TreatExOrEz, double yDist, double xPos, double zPos);
	void SetUpZonePlateNumAccessData(srTSRWRadStructAccessData& Rad, char* NumStructName);
	int ComputeOptPath(srTSRWRadStructAccessData& Rad, char TreatExOrEz, double yDist, double xPos, double zPos);

	int EstimateFocalDist(srTSRWRadStructAccessData* pRad, double yDist)
	{// To improve?
		const double MinFocFact = 1.E-03; // To steer
		const double MaxFocFact = 1.E+10; // To steer
		FocDistX = pRad->RobsX*yDist/(pRad->RobsX + yDist);
		FocDistZ = pRad->RobsZ*yDist/(pRad->RobsZ + yDist);
		double AbsFx = ::fabs(FocDistX), AbsFz = ::fabs(FocDistZ);
		//double AbsRx = ::fabs(pRad->RobsX), AbsRz = ::fabs(pRad->RobsZ);
		if((AbsFx > MinFocFact*AbsFx) && (AbsFx < MaxFocFact*AbsFx) && (AbsFz > MinFocFact*AbsFz) && (AbsFz < MaxFocFact*AbsFz))
			IsFocusing = 1;
		return 0;
	}
	int FetchNumStruct(char* NumStructName);
	//{
	//	srTSend Send;
	//	return Send.FetchNumWave(NumStructName, &ZonePlateNumData);
	//	//DLL_IMPLEMENT
	//}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		if(!IsFocusing) return 0;
		return srTFocusingElem::PropagateRadMoments(pRadAccessData, MomRatArray);
	}
	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(!IsFocusing) return 0;
		return srTFocusingElem::PropagateWaveFrontRadius(pRadAccessData);
	}
	int PropagateWaveFrontRadius1D(srTRadSect1D* pSect1D)
	{
		if(!IsFocusing) return 0;
		return srTFocusingElem::PropagateWaveFrontRadius1D(pSect1D);
	}
	int Propagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData) 
	{
		if(!IsFocusing) return 0;
		return srTFocusingElem::Propagate4x4PropMatr(pRadAccessData);
	}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		char &MethNo = ParPrecWfrPropag.MethNo;
		//if(MethNo == 2) return PropagateRadiationMeth_2(pRadAccessData, ResBeforeAndAfterVect);
		if(MethNo == 2) return PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
		return 0;
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		return TraverseRadZXE(pRadAccessData);
	}
  	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		int result;
		if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
		return TraverseRad1D(pSect1D);
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
		long Nx = (ZonePlateNumData.DimSizes)[0], Nz = (ZonePlateNumData.DimSizes)[1];
		long Nxmi2 = Nx - 2, Nzmi2 = Nz - 2;
		double xStart = (ZonePlateNumData.DimStartValues)[0], zStart = (ZonePlateNumData.DimStartValues)[1];
		double xStep = (ZonePlateNumData.DimSteps)[0], zStep = (ZonePlateNumData.DimSteps)[1];
		double xEnd = xStart + (Nx - 1)*xStep, zEnd = zStart + (Nz - 1)*zStep;
		if((xRel < xStart) || (xRel > xEnd) || (zRel < zStart) || (zRel > zEnd))
		{
			if(EPtrs.pExRe != 0) { *(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;}
			if(EPtrs.pEzRe != 0) { *(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;}
			return;
		}
		long ix = (long)((xRel - xStart)/xStep), iz = (long)((zRel - zStart)/zStep);
		if(ix > Nxmi2) ix = Nxmi2;
		if(iz > Nzmi2) iz = Nzmi2;
		double xr = (xRel - (ix*xStep + xStart))/xStep, zr = (zRel - (iz*zStep + zStart))/zStep;
		float *p00 = (float*)(ZonePlateNumData.pWaveData) + iz*Nx + ix;
		float *p10 = p00 + 1, *p01 = p00 + Nx;
		float *p11 = p01 + 1;
		double TwoPi_d_Lambda_m = EXZ.e*5.0676816042E+06;
		double Ph = TwoPi_d_Lambda_m*((*p00 - *p01 - *p10 + *p11)*xr*zr + (*p10 - *p00)*xr + (*p01 - *p00)*zr + *p00);
		float CosPh, SinPh; CosAndSin(Ph, CosPh, SinPh);
		if(EPtrs.pExRe != 0)
		{
			float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
			float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
			*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		}
		if(EPtrs.pEzRe != 0)
		{
			float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
			float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
			*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
		}
	}
  	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
		long Nx = (ZonePlateNumData.DimSizes)[0], Nz = (ZonePlateNumData.DimSizes)[1];
		long Nxmi2 = Nx - 2, Nzmi2 = Nz - 2;
		double xStart = (ZonePlateNumData.DimStartValues)[0], zStart = (ZonePlateNumData.DimStartValues)[1];
		double xStep = (ZonePlateNumData.DimSteps)[0], zStep = (ZonePlateNumData.DimSteps)[1];
		double xEnd = xStart + (Nx - 1)*xStep, zEnd = zStart + (Nz - 1)*zStep;
		if((xRel < xStart) || (xRel > xEnd) || (zRel < zStart) || (zRel > zEnd))
		{
			if(EPtrs.pExRe != 0) { *(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;}
			if(EPtrs.pEzRe != 0) { *(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;}
			return;
		}
		long ix = (long)((xRel - xStart)/xStep), iz = (long)((zRel - zStart)/zStep);
		if(ix > Nxmi2) ix = Nxmi2;
		if(iz > Nzmi2) iz = Nzmi2;
		double xr = (xRel - (ix*xStep + xStart))/xStep, zr = (zRel - (iz*zStep + zStart))/zStep;
		float *p0 = (float*)(ZonePlateNumData.pWaveData) + iz*Nx + ix;
		double TwoPi_d_Lambda_m = EXZ.e*5.0676816042E+06, Ph;
		if(EXZ.VsXorZ == 'x') Ph = TwoPi_d_Lambda_m*((*(p0 + 1) - *p0)*xr + *p0);
		else Ph = TwoPi_d_Lambda_m*((*(p0 + Nx) - *p0)*zr + *p0);
		float CosPh, SinPh; CosAndSin(Ph, CosPh, SinPh);
		if(EPtrs.pExRe != 0)
		{
			float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
			float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
			*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		}
		if(EPtrs.pEzRe != 0)
		{
			float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
			float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
			*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm;
		}
	}

	int RangeShouldBeAdjustedAtPropag() { return 0;} // Or switch it On
	int ResolutionShouldBeAdjustedAtPropag() { return 1;}
};

//*************************************************************************

#endif
