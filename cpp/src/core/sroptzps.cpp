/************************************************************************//**
 * File: sroptzps.cpp
 * Description: Optical element: "Special" Zone Plate (used to simulate "phase corrections") - Obsolete? 
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include <cstdio> //??

#include "sroptzps.h"
#include "srsend.h"

//*************************************************************************

srTZonePlateSpec::srTZonePlateSpec(srTStringVect* pElemInfo, srTSRWRadStructAccessData* pRad) 
{
	IsFocusing = 0;
	ZonePlateNumData.pWaveData = 0;

	char NumStructName[256];
	strcpy(NumStructName, (*pElemInfo)[1]);
	char FxStr[50], FzStr[50];
	strcpy(FxStr, (*pElemInfo)[2]);
	strcpy(FzStr, (*pElemInfo)[3]);
	if(!((*FxStr == '\0') || (*FzStr == '\0')))
	{
		FocDistX = atof(FxStr);
		FocDistZ = atof(FzStr);
		if((FocDistX == 0.) || (FocDistZ == 0.))
		{
			ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;
		}

		const double MinFocFact = 1.E-03; // To steer
		const double MaxFocFact = 1.E+10; // To steer
		double AbsFx = ::fabs(FocDistX), AbsFz = ::fabs(FocDistZ);
		//double AbsRx = ::fabs(pRad->RobsX), AbsRz = ::fabs(pRad->RobsZ);
		if((AbsFx > MinFocFact*AbsFx) && (AbsFx < MaxFocFact*AbsFx) && (AbsFz > MinFocFact*AbsFz) && (AbsFz < MaxFocFact*AbsFz))
			IsFocusing = 1;

		ErrorCode = FetchNumStruct(NumStructName);
		//DLL_IMPLEMENT

		return;
	}
	char TreatExOrEz = atoi((*pElemInfo)[4]);
	if(TreatExOrEz == 0) TreatExOrEz = 'x';
	else TreatExOrEz = 'z';
	double yDist, xPos, zPos;
	yDist = atof((*pElemInfo)[5]);
	xPos = atof((*pElemInfo)[6]);
	zPos = atof((*pElemInfo)[7]);

	if(ErrorCode = SetUpNumStruct(pRad, NumStructName, TreatExOrEz, yDist, xPos, zPos)) return;

	char* aStr;
	for(int k=2; k<(int)(pElemInfo->size()); k++)
	{
		aStr = (*pElemInfo)[k];
		if(aStr != 0) delete[] aStr;
	}
	pElemInfo->erase(pElemInfo->begin() + 2, pElemInfo->end());

	aStr = new char[256]; sprintf(aStr, "%f", FocDistX);
	pElemInfo->push_back(aStr);
	aStr = new char[256]; sprintf(aStr, "%f", FocDistZ);
	pElemInfo->push_back(aStr);
}

//*************************************************************************

int srTZonePlateSpec::SetUpNumStruct(srTSRWRadStructAccessData* pRad, char* NumStructName, char TreatExOrEz, double yDist, double xPos, double zPos)
{
	int result;
	if(result = EstimateFocalDist(pRad, yDist)) return result;
	SetUpZonePlateNumAccessData(*pRad, NumStructName);
	srTSend Send;
	if(result = Send.MakeWaveAccordingToWaveAccessData(ZonePlateNumData)) return result;
	//DLL_IMPLEMENT

	if(result = ComputeOptPath(*pRad, TreatExOrEz, yDist, xPos, zPos)) return result;
	return 0;
}

//*************************************************************************

void srTZonePlateSpec::SetUpZonePlateNumAccessData(srTSRWRadStructAccessData& Rad, char* NumStructName)
{
	ZonePlateNumData.pWaveData = 0;
	ZonePlateNumData.WaveType[0] = 'f'; // 'f'|'d'|'cf'|'cd'
	ZonePlateNumData.WaveType[1] = '\0'; 
	ZonePlateNumData.AmOfDims = 2;
	ZonePlateNumData.DimSizes[0] = Rad.nx;
	ZonePlateNumData.DimSizes[1] = Rad.nz;
	ZonePlateNumData.DimStartValues[0] = Rad.xStart;
	ZonePlateNumData.DimStartValues[1] = Rad.zStart;
	ZonePlateNumData.DimSteps[0] = Rad.xStep;
	ZonePlateNumData.DimSteps[1] = Rad.zStep;
	strcpy(ZonePlateNumData.DimUnits[0], "m");
	strcpy(ZonePlateNumData.DimUnits[1], "m");
	strcpy(ZonePlateNumData.NameOfWave, NumStructName);
}

//*************************************************************************

int srTZonePlateSpec::ComputeOptPath(srTSRWRadStructAccessData& Rad, char TreatExOrEz, double yDist, double xPos, double zPos)
{
	if(ZonePlateNumData.pWaveData == 0) return 0;

	float *tOpt = (float*)(ZonePlateNumData.pWaveData);
	float *tRad = (TreatExOrEz == 'x')? Rad.pBaseRadX : Rad.pBaseRadZ;
	if(tRad == 0) return IMPROPER_RADIATION_STRUCTURE;

	const double PI = 3.141592653590;
	const double TwoPI = 6.2831853071796;
	double ePh_keV = Rad.eStart;
	if(Rad.PhotEnergyUnit == 0) ePh_keV *= 0.001;
	double WavelengthIn_m = 1.239854*1.E-09/ePh_keV;
	double WaveNum = TwoPI/WavelengthIn_m;
	double InvWaveNum = 1./WaveNum;

	double BufR = 0.5/(::fabs(yDist));

	//long xPer = Rad.ne << 1;
	long long xPer = Rad.ne << 1;
	double zz = Rad.zStart;
	for(long iz=0; iz<Rad.nz; iz++)
	{
		double zr = zz - zPos;
		double zre2 = zr*zr;

		double xx = Rad.xStart;
		for(long ix=0; ix<Rad.nx; ix++)
		{
			float ReRad = *tRad, ImRad = *(tRad + 1);

			double re2RadInv = 1./double(ReRad*ReRad + ImRad*ImRad);
			double ReRadInv = ReRad*re2RadInv, ImRadInv = -ImRad*re2RadInv;

			double xr = xx - xPos;
			double xre2 = xr*xr;
			double re2 = xre2 + zre2;
			double BufRre2 = BufR*re2;
			double BufRe2re2 = BufRre2*BufR;
			double PhFoc = -WaveNum*BufRre2*(1. - BufRe2re2*(1. - BufRe2re2*2.));
			float CosPh, SinPh; CosAndSin(PhFoc, CosPh, SinPh);

			float ReF = (float)(ReRadInv*CosPh - ImRadInv*SinPh);
			float ImF = (float)(ReRadInv*SinPh + ImRadInv*CosPh);
			*(tOpt++) = (float)(InvWaveNum*(FormalPhase(ReF, ImF) + PI));

			tRad += xPer;
			xx += Rad.xStep;
		}
		zz += Rad.zStep;
	}
	return 0;
}

//*************************************************************************

int srTZonePlateSpec::FetchNumStruct(char* NumStructName)
{
	srTSend Send;
	return Send.FetchNumWave(NumStructName, &ZonePlateNumData);
	//DLL_IMPLEMENT
}

//*************************************************************************
