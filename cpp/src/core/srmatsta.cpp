/************************************************************************//**
 * File: srmatsta.cpp
 * Description: Basic statistical characteristics of intensity distributions
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 0.06
 ***************************************************************************/

#include "srmatsta.h"
#include "srradmnp.h"

//*************************************************************************

int srTAuxMatStat::FindSimplestStat(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo)
{
	int result;
	if(result = ValidateSpotData(InWaveData)) return result;

	if(InWaveData.AmOfDims == 1) return FindSimplestStat1D(InWaveData, OutSpotInfo);
	if(InWaveData.AmOfDims == 2) return FindSimplestStat2D(InWaveData, OutSpotInfo);
	return 0;
};

//*************************************************************************

int srTAuxMatStat::FindSimplestStat1D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo)
{
	double xStart = (InWaveData.DimStartValues)[0], xStep = (InWaveData.DimSteps)[0];
	//long LenArr = (InWaveData.DimSizes)[0];
	long long LenArr = (InWaveData.DimSizes)[0];

	double MaxVal;
	//long iMax;
	long long iMax;
	if(*(InWaveData.WaveType) == 'f')
	{
		float* p0 = (float*)(InWaveData.pWaveData);
		FindMax1D(p0, LenArr, MaxVal, iMax);
	}
	else
	{
		DOUBLE* p0 = (DOUBLE*)(InWaveData.pWaveData);
		FindMax1D(p0, (InWaveData.DimSizes)[0], MaxVal, iMax);
	}
	double xMax = xStart + iMax*xStep;

	double HalfMaxVal = 0.5*MaxVal;
	//long iHalfMaxLeft, iHalfMaxRight;
	long long iHalfMaxLeft, iHalfMaxRight;
	double yLeft1, yLeft2, yRight1, yRight2;
	if(*(InWaveData.WaveType) == 'f')
	{
		float* p0 = (float*)(InWaveData.pWaveData);
		FindIndHalfMaxLeftRight1D(p0, LenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);

		yLeft1 = *(p0 + iHalfMaxLeft); yLeft2 = *(p0 + iHalfMaxLeft + 1);
		yRight1 = *(p0 + iHalfMaxRight); yRight2 = *(p0 + iHalfMaxRight + 1);
	}
	else
	{
		DOUBLE* p0 = (DOUBLE*)(InWaveData.pWaveData);
		FindIndHalfMaxLeftRight1D(p0, LenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);

		yLeft1 = *(p0 + iHalfMaxLeft); yLeft2 = *(p0 + iHalfMaxLeft + 1);
		yRight1 = *(p0 + iHalfMaxRight); yRight2 = *(p0 + iHalfMaxRight + 1);
	}

	double xLeft1 = xStart + iHalfMaxLeft*xStep;
	double xLeft = xLeft1 + ((yLeft1 - HalfMaxVal)/(yLeft1 - yLeft2))*xStep;

	double xRight1 = xStart + iHalfMaxRight*xStep;
	double xRight = xRight1 + ((yRight1 - HalfMaxVal)/(yRight1 - yRight2))*xStep;

	float* pOut = (float*)(OutSpotInfo.pWaveData);
	*(pOut++) = (float)MaxVal;
	*(pOut++) = (float)xMax;
	*(pOut++) = 0;
	*(pOut++) = (float)(xRight - xLeft);
	*pOut = 0;
	return 0;
};

//*************************************************************************

int srTAuxMatStat::FindSimplestStat2D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo)
{
	double xStart = (InWaveData.DimStartValues)[0], xStep = (InWaveData.DimSteps)[0];
	//long xLenArr = (InWaveData.DimSizes)[0];
	long long xLenArr = (InWaveData.DimSizes)[0];
	double yStart = (InWaveData.DimStartValues)[1], yStep = (InWaveData.DimSteps)[1];
	//long yLenArr = (InWaveData.DimSizes)[1];
	long long yLenArr = (InWaveData.DimSizes)[1];

	double MaxVal;
	//long iMax, jMax;
	long long iMax, jMax;
	FindMax2D(InWaveData, MaxVal, iMax, jMax);
	double xMax = xStart + iMax*xStep;
	double yMax = yStart + jMax*yStep;

	double HalfMaxVal = 0.5*MaxVal;

	//long iHalfMaxLeft, iHalfMaxRight, jHalfMaxLeft, jHalfMaxRight;
	long long iHalfMaxLeft, iHalfMaxRight, jHalfMaxLeft, jHalfMaxRight;
	double VxLeft1, VxLeft2, VxRight1, VxRight2, VyLeft1, VyLeft2, VyRight1, VyRight2;

	DOUBLE* AuxCont = new DOUBLE[yLenArr];
	if(AuxCont == 0) return MEMORY_ALLOCATION_FAILURE;

	if(*(InWaveData.WaveType) == 'f')
	{
		float* p0 = (float*)(InWaveData.pWaveData) + jMax*xLenArr;
		FindIndHalfMaxLeftRight1D(p0, xLenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);
		VxLeft1 = *(p0 + iHalfMaxLeft); VxLeft2 = *(p0 + iHalfMaxLeft + 1);
		VxRight1 = *(p0 + iHalfMaxRight); VxRight2 = *(p0 + iHalfMaxRight + 1);
		
		DOUBLE *tAux = AuxCont;
		float *t = (float*)(InWaveData.pWaveData) + iMax;
		//for(long k=0; k<yLenArr; k++) { *(tAux++) = *t; t += xLenArr;}
		for(long long k=0; k<yLenArr; k++) { *(tAux++) = *t; t += xLenArr;}
	}
	else
	{
		DOUBLE* p0 = (DOUBLE*)(InWaveData.pWaveData) + jMax*xLenArr;
		FindIndHalfMaxLeftRight1D(p0, xLenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);
		VxLeft1 = *(p0 + iHalfMaxLeft); VxLeft2 = *(p0 + iHalfMaxLeft + 1);
		VxRight1 = *(p0 + iHalfMaxRight); VxRight2 = *(p0 + iHalfMaxRight + 1);

		DOUBLE *tAux = AuxCont;
		DOUBLE *t = (DOUBLE*)(InWaveData.pWaveData) + iMax;
		//for(long k=0; k<yLenArr; k++) { *(tAux++) = *t; t += xLenArr;}
		for(long long k=0; k<yLenArr; k++) { *(tAux++) = *t; t += xLenArr;}
	}
	FindIndHalfMaxLeftRight1D(AuxCont, yLenArr, HalfMaxVal, jHalfMaxLeft, jHalfMaxRight);
	VyLeft1 = *(AuxCont + jHalfMaxLeft); VyLeft2 = *(AuxCont + jHalfMaxLeft + 1);
	VyRight1 = *(AuxCont + jHalfMaxRight); VyRight2 = *(AuxCont + jHalfMaxRight + 1);

	double xLeft1 = xStart + iHalfMaxLeft*xStep;
	double xLeft = xLeft1 + ((VxLeft1 - HalfMaxVal)/(VxLeft1 - VxLeft2))*xStep;
	double xRight1 = xStart + iHalfMaxRight*xStep;
	double xRight = xRight1 + ((VxRight1 - HalfMaxVal)/(VxRight1 - VxRight2))*xStep;
	double yLeft1 = yStart + jHalfMaxLeft*yStep;
	double yLeft = yLeft1 + ((VyLeft1 - HalfMaxVal)/(VyLeft1 - VyLeft2))*yStep;
	double yRight1 = yStart + jHalfMaxRight*yStep;
	double yRight = yRight1 + ((VyRight1 - HalfMaxVal)/(VyRight1 - VyRight2))*yStep;

	delete[] AuxCont;

	float* pOut = (float*)(OutSpotInfo.pWaveData);
	*(pOut++) = (float)MaxVal;
	*(pOut++) = (float)xMax;
	*(pOut++) = (float)yMax;
	*(pOut++) = (float)(xRight - xLeft);
	*pOut = (float)(yRight - yLeft);
	return 0;
};

//*************************************************************************
//Searches for intensity limits: returns a rectangle within which the power is located
int srTAuxMatStat::FindIntensityLimits(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo)
{
	int result;
	if(result = ValidateSpotData(InWaveData)) return result;

	float* pOut = (float*)(OutSpotInfo.pWaveData);
	for(int i=0; i<5; i++) pOut[i] = 0.;
	*pOut = (float)IntegrateSimple(InWaveData);

	if(InWaveData.AmOfDims == 1) return FindIntensityLimits1D(InWaveData, RelPowLevel, OutSpotInfo);
	if(InWaveData.AmOfDims == 2) return FindIntensityLimits2D(InWaveData, RelPowLevel, OutSpotInfo);

	return 0;
};

//*************************************************************************

//void srTAuxMatStat::FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long& iMax, long& jMax)
void srTAuxMatStat::FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long long& iMax, long long& jMax)
{
	//long xLenArr = (InWaveData.DimSizes)[0], yLenArr = (InWaveData.DimSizes)[1];
	long long xLenArr = (InWaveData.DimSizes)[0], yLenArr = (InWaveData.DimSizes)[1];
	MaxVal = -1E+23;

	float* pf = 0;
	DOUBLE* pd = 0;
	if(*(InWaveData.WaveType) == 'f')
	{
		pf = (float*)(InWaveData.pWaveData);
		//for(long j=0; j<yLenArr; j++)
		for(long long j=0; j<yLenArr; j++)
		{
			double MaxValLoc = -1E+23;
			//long iMaxLoc;
			long long iMaxLoc;
			FindMax1D(pf, xLenArr, MaxValLoc, iMaxLoc);
			if(MaxVal < MaxValLoc)
			{
				MaxVal = MaxValLoc; iMax = iMaxLoc; jMax = j;
			}
			pf += xLenArr;
		}
	}
	else
	{
		pd = (DOUBLE*)(InWaveData.pWaveData);
		//for(long j=0; j<yLenArr; j++)
		for(long long j=0; j<yLenArr; j++)
		{
			double MaxValLoc = -1E+23;
			//long iMaxLoc;
			long long iMaxLoc;
			FindMax1D(pd, xLenArr, MaxValLoc, iMaxLoc);
			if(MaxVal < MaxValLoc)
			{
				MaxVal = MaxValLoc; iMax = iMaxLoc; jMax = j;
			}
			pd += xLenArr;
		}
	}
}

//*************************************************************************

//void srTAuxMatStat::FindMax1D(float* p0, long LenArr, double& MaxVal, long& iMax)
void srTAuxMatStat::FindMax1D(float* p0, long long LenArr, double& MaxVal, long long& iMax)
{
	float vMax = (float)(-1.E+23);
	float* t = p0;
	//for(long i=0; i<LenArr; i++)
	for(long long i=0; i<LenArr; i++)
	{
		if(vMax < *t) { vMax = *t; iMax = i;}
		t++;
	}
	MaxVal = double(vMax);
}

//*************************************************************************

//void srTAuxMatStat::FindMax1D(DOUBLE* p0, long LenArr, double& MaxVal, long& iMax)
void srTAuxMatStat::FindMax1D(DOUBLE* p0, long long LenArr, double& MaxVal, long long& iMax)
{
	double vMax = -1.E+23;
	DOUBLE* t = p0;
	//for(long i=0; i<LenArr; i++)
	for(long long i=0; i<LenArr; i++)
	{
		if(vMax < *t) { vMax = *t; iMax = i;}
		t++;
	}
	MaxVal = double(vMax);
}

//*************************************************************************

//void srTAuxMatStat::FindIndHalfMaxLeftRight1D(float* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight)
void srTAuxMatStat::FindIndHalfMaxLeftRight1D(float* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight)
{
	iHalfMaxLeft = 0; iHalfMaxRight = LenArr - 1;

	float* t = p0;
	float HalfMax = (float)HalfMaxVal;
	//for(long i=0; i<LenArr; i++)
	for(long long i=0; i<LenArr; i++)
	{
		if(*(t++) < HalfMax) iHalfMaxLeft = i;
		else break;
	}
	//long LenArr_mi_1 = LenArr - 1;
	long long LenArr_mi_1 = LenArr - 1;
	t = p0 + LenArr_mi_1;
	//for(long j=LenArr_mi_1; j>=0; j--)
	for(long long j=LenArr_mi_1; j>=0; j--)
	{
		iHalfMaxRight = j;
		if(*(t--) >= HalfMax) break;
	}
}

//*************************************************************************

//void srTAuxMatStat::FindIndHalfMaxLeftRight1D(DOUBLE* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight)
void srTAuxMatStat::FindIndHalfMaxLeftRight1D(DOUBLE* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight)
{
	iHalfMaxLeft = 0; iHalfMaxRight = LenArr - 1;

	DOUBLE* t = p0;
	double HalfMax = HalfMaxVal;
	//for(long i=0; i<LenArr; i++)
	for(long long i=0; i<LenArr; i++)
	{
		if(*(t++) < HalfMax) iHalfMaxLeft = i;
		else break;
	}
	//long LenArr_mi_1 = LenArr - 1;
	long long LenArr_mi_1 = LenArr - 1;
	t = p0 + LenArr_mi_1;
	//for(long j=LenArr_mi_1; j>=0; j--)
	for(long long j=LenArr_mi_1; j>=0; j--)
	{
		iHalfMaxRight = j;
		if(*(t--) >= HalfMax) break;
	}
}

//*************************************************************************

int srTAuxMatStat::ValidateSpotData(srTWaveAccessData& InWaveData)
{
	if((*(InWaveData.WaveType) != 'f') && (*(InWaveData.WaveType) != 'd')) return NT_FP32_OR_NT_FP64_WAVE_REQUIRED;
	if(InWaveData.AmOfDims > 2) return WAVE_1D_OR_2D_REQUIRED;

	return 0;
};

//*************************************************************************

double srTAuxMatStat::IntegrateSimple(srTWaveAccessData& InWaveData)
{
	bool Is2D = (InWaveData.AmOfDims == 2);

	//long AmOfVals = (InWaveData.DimSizes)[0];
	long long AmOfVals = (InWaveData.DimSizes)[0];
	if(Is2D) AmOfVals *= (InWaveData.DimSizes)[1];

	double Sum = 0.;
	if(*(InWaveData.WaveType) == 'f') Sum = SumUpArray((float*)(InWaveData.pWaveData), 0, AmOfVals - 1, 1);
	else Sum = SumUpArray((DOUBLE*)(InWaveData.pWaveData), 0, AmOfVals - 1, 1);

	double OutVal = Sum*((InWaveData.DimSteps)[0]);
	if(Is2D) OutVal *= (InWaveData.DimSteps)[1];

	return OutVal;
}

//*************************************************************************

int srTAuxMatStat::FindIntensityLimits1D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo)
{
	//long AmOfVals = (InWaveData.DimSizes)[0];
	long long AmOfVals = (InWaveData.DimSizes)[0];
	if(AmOfVals <= 0) return 0;
	double xStart = (InWaveData.DimStartValues)[0], xStep = (InWaveData.DimSteps)[0];
	
	double TotPower = (*((float*)(OutSpotInfo.pWaveData)))/xStep;
	double AbsPowerToStopOn = TotPower*(1. - RelPowLevel)*0.5;

	float* pOut = (float*)(OutSpotInfo.pWaveData);
	if(*(InWaveData.WaveType) == 'f') 
	{
		*(pOut + 1) = (float)(xStart + xStep*FindLimit1DLeft((float*)(InWaveData.pWaveData), AmOfVals, AbsPowerToStopOn));
		*(pOut + 2) = (float)(xStart + xStep*FindLimit1DRight((float*)(InWaveData.pWaveData), AmOfVals, AbsPowerToStopOn));
	}
	else
	{
		*(pOut + 1) = (float)(xStart + xStep*FindLimit1DLeft((DOUBLE*)(InWaveData.pWaveData), AmOfVals, AbsPowerToStopOn));
		*(pOut + 2) = (float)(xStart + xStep*FindLimit1DRight((DOUBLE*)(InWaveData.pWaveData), AmOfVals, AbsPowerToStopOn));
	}
	return 0;
}

//*************************************************************************

//template <class T> long srTAuxMatStat::FindLimit1DLeft(T* p0, long n, double IntToStop)
template <class T> long long srTAuxMatStat::FindLimit1DLeft(T* p0, long long n, double IntToStop)
{
	if((p0 == 0) || (n == 0)) return 0;

	double Sum = 0.;
	//long iStop = -1;
	long long iStop = -1;
	T* t = p0;
	//for(long i=0; i<n; i++) 
	for(long long i=0; i<n; i++) 
	{
		Sum += (*(t++));
		if(Sum > IntToStop)
		{
			iStop = i; break;
		}
	}
	if((iStop < 0) || (iStop >= n)) return 0;
	//if(iStop > 0) iStop--;

	return iStop;
}

//*************************************************************************

//template <class T> long srTAuxMatStat::FindLimit1DRight(T* p0, long n, double IntToStop)
template <class T> long long srTAuxMatStat::FindLimit1DRight(T* p0, long long n, double IntToStop)
{
	if((p0 == 0) || (n == 0)) return (n - 1);

	double Sum = 0.;
	//long iStop = -1;
	long long iStop = -1;
	T* t = p0 + (n - 1);

	//for(long i=0; i<n; i++) 
	for(long long i=0; i<n; i++) 
	{
		Sum += (*(t--));
		if(Sum > IntToStop)
		{
			iStop = i; break;
		}
	}
	if((iStop < 0) || (iStop >= n)) return (n - 1);
	//if(iStop > 0) iStop--;

	return ((n - 1) - iStop);
}

//*************************************************************************

int srTAuxMatStat::FindIntensityLimits2D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo)
{
	long Nx = (InWaveData.DimSizes)[0], Ny = (InWaveData.DimSizes)[1];
	if((Nx <= 0) || (Ny <= 0)) return INCORRECT_ARGUMENTS;

	double xStart = (InWaveData.DimStartValues)[0], xStep = (InWaveData.DimSteps)[0];
	double yStart = (InWaveData.DimStartValues)[1], yStep = (InWaveData.DimSteps)[1];

	double TotPower = (*((float*)(OutSpotInfo.pWaveData)));
	double AbsPowerToStopOn = TotPower*(1. - RelPowLevel)*0.25;

	float *pf0 = 0;
	DOUBLE *pD0 = 0;
	if(*(InWaveData.WaveType) == 'f') pf0 = (float*)(InWaveData.pWaveData);
	else pD0 = (DOUBLE*)(InWaveData.pWaveData);

	int res = 0;
	double* AuxArrIntOverX = new double[Ny];

	if(pf0 != 0) 
	{ 
		if(res = IntegrateOverX(pf0, 0, Nx - 1, xStep, Nx, Ny, AuxArrIntOverX)) 
		{
			delete[] AuxArrIntOverX; return res;
		}
	}
	else 
	{ 
		if(res = IntegrateOverX(pD0, 0, Nx - 1, xStep, Nx, Ny, AuxArrIntOverX)) 
		{
			delete[] AuxArrIntOverX; return res;
		}
	}
	double CurAbsPowerToStopOn = AbsPowerToStopOn/yStep;
	//long iyStart = FindLimit1DLeft(AuxArrIntOverX, Ny, CurAbsPowerToStopOn);
	//long iyEnd = FindLimit1DRight(AuxArrIntOverX, Ny, CurAbsPowerToStopOn);
	long long iyStart = FindLimit1DLeft(AuxArrIntOverX, Ny, CurAbsPowerToStopOn);
	long long iyEnd = FindLimit1DRight(AuxArrIntOverX, Ny, CurAbsPowerToStopOn);
	if(iyEnd <= iyStart) 
	{
		delete[] AuxArrIntOverX; return INCORRECT_ARGUMENTS;
	}

	float* pOut = (float*)(OutSpotInfo.pWaveData);
	*(pOut + 3) = (float)(yStart + yStep*iyStart);
	*(pOut + 4) = (float)(yStart + yStep*iyEnd);

	double* AuxArrIntOverY = new double[Nx];

	if(pf0 != 0) 
	{ 
		if(res = IntegrateOverY(pf0, iyStart, iyEnd, yStep, Nx, AuxArrIntOverY)) 
		{
			delete[] AuxArrIntOverY; return res;
		}
	}
	else 
	{ 
		if(res = IntegrateOverY(pD0, iyStart, iyEnd, yStep, Nx, AuxArrIntOverY)) 
		{
			delete[] AuxArrIntOverY; return res;
		}
	}
	CurAbsPowerToStopOn = AbsPowerToStopOn/xStep;
	//long ixStart = FindLimit1DLeft(AuxArrIntOverY, Nx, CurAbsPowerToStopOn);
	//long ixEnd = FindLimit1DRight(AuxArrIntOverY, Nx, CurAbsPowerToStopOn);
	long long ixStart = FindLimit1DLeft(AuxArrIntOverY, Nx, CurAbsPowerToStopOn);
	long long ixEnd = FindLimit1DRight(AuxArrIntOverY, Nx, CurAbsPowerToStopOn);
	*(pOut + 1) = (float)(xStart + xStep*ixStart);
	*(pOut + 2) = (float)(xStart + xStep*ixEnd);

	if(AuxArrIntOverX != 0) delete[] AuxArrIntOverX;
	if(AuxArrIntOverY != 0) delete[] AuxArrIntOverY;
	return 0;
}

//*************************************************************************

//template <class T> int srTAuxMatStat::IntegrateOverX(T* p0, long ixStart, long ixEnd, double xStep, long Nx, long Ny, double* AuxArrIntOverX)
template <class T> int srTAuxMatStat::IntegrateOverX(T* p0, long long ixStart, long long ixEnd, double xStep, long long Nx, long long Ny, double* AuxArrIntOverX)
{
	if((p0 == 0) || (AuxArrIntOverX == 0) || (Nx <= 0) || (Ny <= 0) || (ixStart >= ixEnd) || (ixEnd <= 0)) return INCORRECT_ARGUMENTS;
	if(ixStart < 0) ixStart = 0;
	
	T* t = p0;
	//for(int i=0; i<Ny; i++)
	for(long long i=0; i<Ny; i++)
	{
		AuxArrIntOverX[i] = xStep*SumUpArray(t, ixStart, ixEnd, 1);
		t += Nx;
	}
	return 0;
}

//*************************************************************************

//template <class T> int srTAuxMatStat::IntegrateOverY(T* p0, long iyStart, long iyEnd, double yStep, long Nx, double* AuxArrIntOverY)
template <class T> int srTAuxMatStat::IntegrateOverY(T* p0, long long iyStart, long long iyEnd, double yStep, long long Nx, double* AuxArrIntOverY)
{
	if((p0 == 0) || (AuxArrIntOverY == 0) || (Nx <= 0) || (iyStart >= iyEnd) || (iyEnd <= 0)) return INCORRECT_ARGUMENTS;
	if(iyStart < 0) iyStart = 0;

	T* t = p0;
	//long iyStartNx = iyStart*Nx, iyEndNx = iyEnd*Nx;
	long long iyStartNx = iyStart*Nx, iyEndNx = iyEnd*Nx;
	//for(int i=0; i<Nx; i++)
	for(long long i=0; i<Nx; i++)
	{
		AuxArrIntOverY[i] = yStep*SumUpArray(t, iyStartNx, iyEndNx, Nx);
		t++;
	}
	return 0;
}

//*************************************************************************

//int srTAuxMatStat::FindIntensityLimitsInds(srTSRWRadStructAccessData& Rad, int ie, double RelPow, int* IndLims)
int srTAuxMatStat::FindIntensityLimitsInds(CHGenObj& hRad, int ie, double RelPow, int* IndLims)
{
	srTSRWRadStructAccessData& Rad = *((srTSRWRadStructAccessData*)(hRad.ptr()));

	IndLims[0] = 0;
	IndLims[1] = Rad.nx - 1;
	IndLims[2] = 0;
	IndLims[3] = Rad.nz - 1;

	srTRadExtract RadExtract;
	RadExtract.PolarizCompon = 6;
	RadExtract.Int_or_Phase = 0;
	RadExtract.PlotType = 3;
	RadExtract.TransvPres = Rad.Pres;

	RadExtract.ePh = Rad.eStart + ie*Rad.eStep;
	RadExtract.pExtractedData = new float[Rad.nx*Rad.nz];

	//srTRadGenManip RadGenManip(Rad);
	srTRadGenManip RadGenManip(hRad);
	srTWaveAccessData ExtractedWaveData;
	int res = 0;
	if(res = RadGenManip.ExtractRadiation(RadExtract, ExtractedWaveData)) 
	{
		delete[] RadExtract.pExtractedData; return res;
	}

	float AuxArrF[5];
	srTWaveAccessData OutInfoData;
	//char InfoSuf[] = "_inf"; // To steer //OC030110
	(OutInfoData.WaveType)[0] = 'f';
	OutInfoData.AmOfDims = 1;
	(OutInfoData.DimSizes)[0] = 5;
	(OutInfoData.DimSizes)[1] = 0;
	(OutInfoData.DimStartValues)[0] = 0;
	(OutInfoData.DimSteps)[0] = 1;
	OutInfoData.pWaveData = (char*)AuxArrF;
	for(int i=0; i<5; i++) AuxArrF[i] = 0.;
	AuxArrF[0] = (float)IntegrateSimple(ExtractedWaveData);
	if(res = FindIntensityLimits2D(ExtractedWaveData, RelPow, OutInfoData)) 
	{
		delete[] RadExtract.pExtractedData; return res;
	}

	IndLims[0] = (int)((AuxArrF[1] - Rad.xStart)*1.0000001/Rad.xStep);
	if(IndLims[0] < 0) IndLims[0] = 0;
	IndLims[1] = (int)((AuxArrF[2] - Rad.xStart)*1.0000001/Rad.xStep);
	if(IndLims[1] >= Rad.nx) IndLims[1] = Rad.nx - 1;
	IndLims[2] = (int)((AuxArrF[3] - Rad.zStart)*1.0000001/Rad.zStep);
	if(IndLims[2] < 0) IndLims[2] = 0;
	IndLims[3] = (int)((AuxArrF[4] - Rad.zStart)*1.0000001/Rad.zStep);
	if(IndLims[3] >= Rad.nz) IndLims[3] = Rad.nz - 1;

	delete[] RadExtract.pExtractedData;
	return 0;
}

//*************************************************************************
