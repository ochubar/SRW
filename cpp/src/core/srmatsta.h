/************************************************************************//**
 * File: srmatsta.h
 * Description: Basic statistical characteristics of intensity distributions (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMATSTA_H
#define __SRMATSTA_H

#include "srobject.h"

#ifdef __IGOR_PRO__
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#else
#ifndef __SRIGORRE_H
#include "srigorre.h"
#endif
#endif

#include <complex>

//*************************************************************************

struct srTWaveAccessData;

//*************************************************************************

class srTAuxMatStat {
public:

	int FindSimplestStat(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	int FindSimplestStat1D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	int FindSimplestStat2D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	int ValidateSpotData(srTWaveAccessData& InWaveData);
	//void FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long& iMax, long& jMax);
	void FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long long& iMax, long long& jMax);
	
	int FindIntensityLimits(srTWaveAccessData& InWaveData, double PowLevel, srTWaveAccessData& OutSpotInfoData);
	int FindIntensityLimits1D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo);
	int FindIntensityLimits2D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo);
	//int FindIntensityLimitsInds(srTSRWRadStructAccessData&, int ie, double RelPow, int* IndLims);
	int FindIntensityLimitsInds(CHGenObj&, int ie, double RelPow, int* IndLims);

	//void FindMax1D(float* p0, long LenArr, double& MaxVal, long& iMax);
	//void FindMax1D(DOUBLE* p0, long LenArr, double& MaxVal, long& iMax);
	void FindMax1D(float* p0, long long LenArr, double& MaxVal, long long& iMax);
	void FindMax1D(DOUBLE* p0, long long LenArr, double& MaxVal, long long& iMax);
	//void FindIndHalfMaxLeftRight1D(float* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight);
	//void FindIndHalfMaxLeftRight1D(DOUBLE* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight);
	void FindIndHalfMaxLeftRight1D(float* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight);
	void FindIndHalfMaxLeftRight1D(DOUBLE* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight);
	double IntegrateSimple(srTWaveAccessData& InWaveData);

	//template <class T> long FindLimit1DLeft(T* p0, long n, double IntToStop);
	//template <class T> long FindLimit1DRight(T* p0, long n, double IntToStop);
	template <class T> long long FindLimit1DLeft(T* p0, long long n, double IntToStop);
	template <class T> long long FindLimit1DRight(T* p0, long long n, double IntToStop);
	//template <class T> int IntegrateOverX(T* p0, long ixStart, long ixEnd, double xStep, long Nx, long Ny, double* AuxArrIntOverX);
	//template <class T> int IntegrateOverY(T* p0, long iyStart, long iyEnd, double yStep, long Nx, double* AuxArrIntOverY);
	template <class T> int IntegrateOverX(T* p0, long long ixStart, long long ixEnd, double xStep, long long Nx, long long Ny, double* AuxArrIntOverX);
	template <class T> int IntegrateOverY(T* p0, long long iyStart, long long iyEnd, double yStep, long long Nx, double* AuxArrIntOverY);

	//template <class T> double SumUpArray(T* p0, long iStart, long iEnd, long iStep)
	template <class T> double SumUpArray(T* p0, long long iStart, long long iEnd, long long iStep)
	{
		if((p0 == 0) || (iStart >= iEnd) || (iEnd <= 0)) return 0.;
		if(iStart < 0) iStart = 0;
		double Sum = 0.;
		//for(long i=iStart; i<=iEnd; i+=iStep) Sum += *(p0 + i);
		for(long long i=iStart; i<=iEnd; i+=iStep) Sum += *(p0 + i);
		return Sum;
	}
};

//*************************************************************************

#endif
