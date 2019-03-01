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
//#include "srercode.h"
#include "srstraux.h" //OC07012019 (added to make it compile on GCC)

#ifdef __IGOR_PRO__
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#else
#ifndef __SRIGORRE_H
#include "srigorre.h"
#endif
#endif

#include <complex>

//*************************************************************************

//struct srTWaveAccessData; //OC07012019 (commented-out to make it compile on GCC)

//*************************************************************************

class srTAuxMatStat {
public:

	int FindSimplestStat(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo, double* arPar=0, int nPar=0); //OC29122018
	//int FindSimplestStat1D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo, double* arPar=0);
	//int FindSimplestStat2D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo, double* arPar=0);
	//int FindSimplestStat(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	//int FindSimplestStat1D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	//int FindSimplestStat2D(srTWaveAccessData& InWaveData, srTWaveAccessData& OutSpotInfo);
	int ValidateSpotData(srTWaveAccessData& InWaveData);
	//void FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long& iMax, long& jMax);
	//void FindMax2D(srTWaveAccessData& InWaveData, double& MaxVal, long long& iMax, long long& jMax);
	
	int FindIntensityLimits(srTWaveAccessData& InWaveData, double PowLevel, srTWaveAccessData& OutSpotInfoData);
	int FindIntensityLimits1D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo);
	int FindIntensityLimits2D(srTWaveAccessData& InWaveData, double RelPowLevel, srTWaveAccessData& OutSpotInfo);
	//int FindIntensityLimitsInds(srTSRWRadStructAccessData&, int ie, double RelPow, int* IndLims);
	int FindIntensityLimitsInds(CHGenObj&, int ie, double RelPow, int* IndLims);

	//void FindMax1D(float* p0, long LenArr, double& MaxVal, long& iMax);
	//void FindMax1D(DOUBLE* p0, long LenArr, double& MaxVal, long& iMax);
	//void FindMax1D(float* p0, long long LenArr, double& MaxVal, long long& iMax);
	//void FindMax1D(DOUBLE* p0, long long LenArr, double& MaxVal, long long& iMax);

	//void FindIndHalfMaxLeftRight1D(float* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight);
	//void FindIndHalfMaxLeftRight1D(DOUBLE* p0, long LenArr, double HalfMaxVal, long& iHalfMaxLeft, long& iHalfMaxRight);
	//void FindIndHalfMaxLeftRight1D(float* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight);
	//void FindIndHalfMaxLeftRight1D(DOUBLE* p0, long long LenArr, double HalfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight);
	double IntegrateSimple(srTWaveAccessData& InWaveData);

	//template <class T> long FindLimit1DLeft(T* p0, long n, double IntToStop);
	//template <class T> long FindLimit1DRight(T* p0, long n, double IntToStop);
	template <class T> long long FindLimit1DLeft(T* p0, long long n, double IntToStop);
	template <class T> long long FindLimit1DRight(T* p0, long long n, double IntToStop);
	//template <class T> int IntegrateOverX(T* p0, long ixStart, long ixEnd, double xStep, long Nx, long Ny, double* AuxArrIntOverX);
	//template <class T> int IntegrateOverY(T* p0, long iyStart, long iyEnd, double yStep, long Nx, double* AuxArrIntOverY);
	template <class T> int IntegrateOverX(T* p0, long long ixStart, long long ixEnd, double xStep, long long Nx, long long Ny, double* AuxArrIntOverX);
	template <class T> int IntegrateOverY(T* p0, long long iyStart, long long iyEnd, double yStep, long long Nx, double* AuxArrIntOverY);

	template <class T> void FindIndHalfMaxLeftRightFromPeak1D(T* ar, long long lenArr, double halfMaxVal, long long iMax, long long& iHalfMaxLeft, long long& iHalfMaxRight)
	{//OC02012019
		iHalfMaxLeft = 0; iHalfMaxRight = lenArr - 1;

		T halfMax = (T)halfMaxVal;
		//T *t0 = p0 + iMax;
		//T *t = t0;
		for(long long i=(iMax-1); i>=0; i--)
		{
			//if(*(t--) <= halfMax) { iHalfMaxLeft = i; break;}
			if(ar[i] <= halfMax) { iHalfMaxLeft = i; break;}
		}
		//t = t0;
		for(long long i=(iMax+1); i<lenArr; i++)
		{
			//if(*(t++) < halfMax) { iHalfMaxRight = i - 1; break;}
			if(ar[i] < halfMax) { iHalfMaxRight = i - 1; break;}
		}
	}

	template <class T> void FindIndHalfMaxLeftRight1D(T* p0, long long lenArr, double halfMaxVal, long long& iHalfMaxLeft, long long& iHalfMaxRight)
	{
		iHalfMaxLeft = 0; iHalfMaxRight = lenArr - 1;

		T* t = p0;
		T halfMax = (T)halfMaxVal;
		for(long long i=0; i<lenArr; i++)
		{
			if(*(t++) < halfMax) iHalfMaxLeft = i;
			else break;
		}
		long long lenArr_mi_1 = lenArr - 1;
		t = p0 + lenArr_mi_1;
		for(long long j=lenArr_mi_1; j>=0; j--)
		{
			iHalfMaxRight = j;
			if(*(t--) >= halfMax) break;
		}
	}

	template <class T> void FindMax1D(T* p0, long long lenArr, double& maxVal, long long& iMax, long long per=1) //OC29122018
	{
		T *t = p0;
		T vMax = *t; //(float)(-1.E+23);
		long long iMaxLoc = 0;
		for(long long i=0; i<lenArr; i++)
		{
			if(vMax < *t) { vMax = *t; iMaxLoc = i;}
			//t++;
			t += per; //OC02012019
		}
		maxVal = double(vMax);
		iMax = iMaxLoc;
	}

	template <class T> void FindMaxExt1D(T* p0, long long lenArr, long long per, double& maxVal, long long& iMax) //OC02012019
	{
		T *t = p0;
		T vMax = *t; //(float)(-1.E+23);
		long long iMaxLoc = 0;
		for(long long i=0; i<lenArr; i++)
		{
			if(vMax < *t) { vMax = *t; iMaxLoc = i;}
			t++;
		}
		maxVal = double(vMax);
		iMax = iMaxLoc;
	}

	template <class T> void FillOutResData(T* pOut, double* arRes, int nRes) //OC03012019
	//template <class T> void FillOutResData(T* pOut, double valMax, double xMax, double yMax, double xFWHM, double yFWHM)
	{
		if(pOut == 0) return;

		//*(pOut++) = (T)valMax;
		//*(pOut++) = (T)xMax;
		//*(pOut++) = (T)yMax;
		//*(pOut++) = (T)(xFWHM);
		//*pOut = (T)(yFWHM);

		double *tRes = arRes;
		for(int i=0; i<nRes; i++) *(pOut++) = (T)(*(tRes++));
	}

	template <class T> void FindSimplestStat1D(T* p0, double xStart, double xStep, long long lenArr, srTWaveAccessData& OutSpotInfo, double* arPar=0, int nPar=0) //OC29122018
	{
		double MaxVal;
		long long iMax;
		FindMax1D(p0, lenArr, MaxVal, iMax);
		double xMax = xStart + iMax*xStep;

		double HalfMaxVal = 0.5*MaxVal;
		long long iHalfMaxLeft, iHalfMaxRight;

		bool calcSizeFromPeak = false;
		if(arPar != 0) { if(arPar[0] == 1) calcSizeFromPeak = true;}

		if(calcSizeFromPeak) FindIndHalfMaxLeftRightFromPeak1D(p0, lenArr, HalfMaxVal, iMax, iHalfMaxLeft, iHalfMaxRight);
		else FindIndHalfMaxLeftRight1D(p0, lenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);

		double yLeft1 = *(p0 + iHalfMaxLeft), yLeft2 = *(p0 + iHalfMaxLeft + 1);
		double yRight1 = *(p0 + iHalfMaxRight), yRight2 = *(p0 + iHalfMaxRight + 1);
		double xLeft1 = xStart + iHalfMaxLeft*xStep;
		double xLeft = xLeft1 + ((yLeft1 - HalfMaxVal)/(yLeft1 - yLeft2))*xStep;
		double xRight1 = xStart + iHalfMaxRight*xStep;
		double xRight = xRight1 + ((yRight1 - HalfMaxVal)/(yRight1 - yRight2))*xStep;
		double xFWHM = xRight - xLeft;

		double xFWFM = 0.;
		int nRes = 5;
		if(nPar > 1)
		{//OC03012019
			HalfMaxVal = arPar[1]*MaxVal;

			if(calcSizeFromPeak) FindIndHalfMaxLeftRightFromPeak1D(p0, lenArr, HalfMaxVal, iMax, iHalfMaxLeft, iHalfMaxRight);
			else FindIndHalfMaxLeftRight1D(p0, lenArr, HalfMaxVal, iHalfMaxLeft, iHalfMaxRight);

			yLeft1 = *(p0 + iHalfMaxLeft); yLeft2 = *(p0 + iHalfMaxLeft + 1);
			yRight1 = *(p0 + iHalfMaxRight); yRight2 = *(p0 + iHalfMaxRight + 1);
			xLeft1 = xStart + iHalfMaxLeft*xStep;
			xLeft = xLeft1 + ((yLeft1 - HalfMaxVal)/(yLeft1 - yLeft2))*xStep;
			xRight1 = xStart + iHalfMaxRight*xStep;
			xRight = xRight1 + ((yRight1 - HalfMaxVal)/(yRight1 - yRight2))*xStep;
			xFWFM = xRight - xLeft;
			nRes = 7;
		}

		//float* pOut = (float*)(OutSpotInfo.pWaveData); //make it more general (allow double and float)
		//*(pOut++) = (float)MaxVal;
		//*(pOut++) = (float)xMax;
		//*(pOut++) = 0;
		//*(pOut++) = (float)(xRight - xLeft);
		//*pOut = 0;

		double arRes[] = {MaxVal, xMax, 0., xFWHM, 0., xFWFM, 0.};
		if(*(OutSpotInfo.WaveType) == 'f') FillOutResData((float*)(OutSpotInfo.pWaveData), arRes, nRes); //OC03012019
		else FillOutResData((double*)(OutSpotInfo.pWaveData), arRes, nRes);
		//if(*(OutSpotInfo.WaveType) == 'f') FillOutResData((float*)(OutSpotInfo.pWaveData), MaxVal, xMax, 0., xFWHM, 0., xFWFM, 0.);
		//else FillOutResData((double*)(OutSpotInfo.pWaveData), MaxVal, xMax, 0., xFWHM, 0., xFWFM, 0.);
	}

	template <class T> void FindIndHalfMaxLeftRight2D(T* p0, long long xLenArr, long long yLenArr, double HalfMaxVal, long long ixMax, long long iyMax, bool searchFromPeak, long long *arResInd, double *arResVxVy)
	{//02012019
		long long &ixHalfMaxLeft = arResInd[0], &ixHalfMaxRight = arResInd[1], &iyHalfMaxLeft = arResInd[2], &iyHalfMaxRight = arResInd[3];
		//ixHalfMaxLeft = 0; ixHalfMaxRight = xLenArr - 1;
		ixHalfMaxLeft = -1; ixHalfMaxRight = -1;
		iyHalfMaxLeft = -1; iyHalfMaxRight = -1;

		for(int ii=0; ii<8; ii++) arResVxVy[ii] = HalfMaxVal;
		double &VxLeft1 = arResVxVy[0], &VxLeft2 = arResVxVy[1], &VxRight1 = arResVxVy[2], &VxRight2 = arResVxVy[3];
		double &VyLeft1 = arResVxVy[4], &VyLeft2 = arResVxVy[5], &VyRight1 = arResVxVy[6], &VyRight2 = arResVxVy[7];

		double curTest;
		long long iyAuxMax = iyMax;
		if(searchFromPeak)
		{
			for(long long ix=(ixMax-1); ix>=0; ix--)
			{
				FindMax1D(p0 + ix, yLenArr, curTest, iyAuxMax, xLenArr);
				if(curTest <= HalfMaxVal) { ixHalfMaxLeft = ix; break;}
			}
		}
		else //serach from extremities
		{
			for(long long ix=0; ix<ixMax; ix++)
			{
				FindMax1D(p0 + ix, yLenArr, curTest, iyAuxMax, xLenArr);
				if(curTest > HalfMaxVal) { ixHalfMaxLeft = ix - 1; break;}
			}
		}
		if(ixHalfMaxLeft < 0)
		{//OC05012019
			ixHalfMaxLeft = ixMax - 1;
			if(ixHalfMaxLeft < 0) ixHalfMaxLeft = 0;
		}
		long long ofst = ixHalfMaxLeft + iyAuxMax*xLenArr;
		VxLeft1 = (double)(*(p0 + ofst));
		VxLeft2 = (double)(*(p0 + (ofst + 1)));
	
		iyAuxMax = iyMax;
		if(searchFromPeak)
		{
			for(long long ix=(ixMax+1); ix<xLenArr; ix++)
			{
				FindMax1D(p0 + ix, yLenArr, curTest, iyAuxMax, xLenArr);
				if(curTest < HalfMaxVal) { ixHalfMaxRight = ix - 1; break;}
			}
		}
		else //serach from extremities
		{
			for(long long ix=(xLenArr-1); ix>=ixMax; ix--)
			{
				FindMax1D(p0 + ix, yLenArr, curTest, iyAuxMax, xLenArr);
				if(curTest >= HalfMaxVal) { ixHalfMaxRight = ix; break;}
			}
		}
		if(ixHalfMaxRight < 0)
		{//OC05012019
			ixHalfMaxRight = ixMax + 1;
			if(ixHalfMaxRight >= xLenArr) ixHalfMaxRight = xLenArr - 1;
		}
		ofst = ixHalfMaxRight + iyAuxMax*xLenArr;
		VxRight1 = (double)(*(p0 + ofst));
		VxRight2 = (double)(*(p0 + (ofst + 1)));

		long long ixAuxMax = ixMax;
		if(searchFromPeak)
		{
			for(long long iy=(iyMax-1); iy>=0; iy--)
			{
				FindMax1D(p0 + iy*xLenArr, xLenArr, curTest, ixAuxMax);
				if(curTest <= HalfMaxVal) { iyHalfMaxLeft = iy; break;}
			}
		}
		else //serach from extremities
		{
			for(long long iy=0; iy<iyMax; iy++)
			{
				FindMax1D(p0 + iy*xLenArr, xLenArr, curTest, ixAuxMax);
				if(curTest > HalfMaxVal) { iyHalfMaxLeft = iy - 1; break;}
			}
		}
		if(iyHalfMaxLeft < 0)
		{//OC05012019
			iyHalfMaxLeft = iyMax - 1;
			if(iyHalfMaxLeft < 0) iyHalfMaxLeft = 0;
		}
		ofst = ixAuxMax + iyHalfMaxLeft*xLenArr;
		VyLeft1 = (double)(*(p0 + ofst));
		VyLeft2 = (double)(*(p0 + (ofst + xLenArr)));

		ixAuxMax = ixMax;
		if(searchFromPeak)
		{
			for(long long iy=(iyMax+1); iy<yLenArr; iy++)
			{
				FindMax1D(p0 + iy*xLenArr, xLenArr, curTest, ixAuxMax);
				if(curTest < HalfMaxVal) { iyHalfMaxRight = iy - 1; break;}
			}
		}
		else //serach from extremities
		{
			for(long long iy=(yLenArr-1); iy>=iyMax; iy--)
			{
				FindMax1D(p0 + iy*xLenArr, xLenArr, curTest, ixAuxMax);
				if(curTest >= HalfMaxVal) { iyHalfMaxRight = iy; break;}
			}
		}
		if(iyHalfMaxRight < 0)
		{//OC05012019
			iyHalfMaxRight = iyMax + 1;
			if(iyHalfMaxRight >= yLenArr) iyHalfMaxRight = yLenArr - 1;
		}
		ofst = ixAuxMax + iyHalfMaxRight*xLenArr;
		VyRight1 = (double)(*(p0 + ofst));
		VyRight2 = (double)(*(p0 + (ofst + xLenArr)));
	}

	void AuxCalcWidth2D(double xStart, double xStep, double yStart, double yStep, double HalfMaxVal, long long arResInd[4], double arResVxVy[8], double& xWidth, double& yWidth)
	{
		long long ixHalfMaxLeft = arResInd[0], ixHalfMaxRight = arResInd[1], iyHalfMaxLeft = arResInd[2], iyHalfMaxRight = arResInd[3];
		double VxLeft1 = arResVxVy[0], VxLeft2 = arResVxVy[1], VxRight1 = arResVxVy[2], VxRight2 = arResVxVy[3];
		double VyLeft1 = arResVxVy[4], VyLeft2 = arResVxVy[5], VyRight1 = arResVxVy[6], VyRight2 = arResVxVy[7];

		double xLeft1 = xStart + ixHalfMaxLeft*xStep;
		double xLeft = xLeft1 + ((VxLeft1 - HalfMaxVal)/(VxLeft1 - VxLeft2))*xStep;
		double xRight1 = xStart + ixHalfMaxRight*xStep;
		double xRight = xRight1 + ((VxRight1 - HalfMaxVal)/(VxRight1 - VxRight2))*xStep;
		xWidth = xRight - xLeft;
		double yLeft1 = yStart + iyHalfMaxLeft*yStep;
		double yLeft = yLeft1 + ((VyLeft1 - HalfMaxVal)/(VyLeft1 - VyLeft2))*yStep;
		double yRight1 = yStart + iyHalfMaxRight*yStep;
		double yRight = yRight1 + ((VyRight1 - HalfMaxVal)/(VyRight1 - VyRight2))*yStep;
		yWidth = yRight - yLeft;
	}

	template <class T> void FindSimplestStat2D(T* p0, double xStart, double xStep, long long xLenArr, double yStart, double yStep, long long yLenArr, srTWaveAccessData& OutSpotInfo, double* arPar=0, int nPar=0) //OC02012019
	{
		double MaxVal;
		long long iMaxTot, lenArrTot = xLenArr*yLenArr;
		//FindMax2D(InWaveData, MaxVal, iMax, jMax);
		FindMax1D(p0, lenArrTot, MaxVal, iMaxTot);
		long long iyMax = (long long)(iMaxTot/xLenArr);
		long long ixMax = iMaxTot - iyMax*xLenArr;

		double xMax = xStart + ixMax*xStep;
		double yMax = yStart + iyMax*yStep;

		double HalfMaxVal = 0.5*MaxVal;

		bool calcSizeFromPeak = false;
		if(arPar != 0) { if(arPar[0] == 1) calcSizeFromPeak = true;}

		long long arResInd[4];
		double arResVxVy[8]; 
		FindIndHalfMaxLeftRight2D(p0, xLenArr, yLenArr, HalfMaxVal, ixMax, iyMax, calcSizeFromPeak, arResInd, arResVxVy);

		double xFWHM = 0., yFWHM = 0.;
		AuxCalcWidth2D(xStart, xStep, yStart, yStep, HalfMaxVal, arResInd, arResVxVy, xFWHM, yFWHM);
		//long long ixHalfMaxLeft = arResInd[0], ixHalfMaxRight = arResInd[1], iyHalfMaxLeft = arResInd[2], iyHalfMaxRight = arResInd[3];
		//double VxLeft1 = arResVxVy[0], VxLeft2 = arResVxVy[1], VxRight1 = arResVxVy[2], VxRight2 = arResVxVy[3];
		//double VyLeft1 = arResVxVy[4], VyLeft2 = arResVxVy[5], VyRight1 = arResVxVy[6], VyRight2 = arResVxVy[7];

		//double xLeft1 = xStart + ixHalfMaxLeft*xStep;
		//double xLeft = xLeft1 + ((VxLeft1 - HalfMaxVal)/(VxLeft1 - VxLeft2))*xStep;
		//double xRight1 = xStart + ixHalfMaxRight*xStep;
		//double xRight = xRight1 + ((VxRight1 - HalfMaxVal)/(VxRight1 - VxRight2))*xStep;
		//double xFWHM = xRight - xLeft;
		//double yLeft1 = yStart + iyHalfMaxLeft*yStep;
		//double yLeft = yLeft1 + ((VyLeft1 - HalfMaxVal)/(VyLeft1 - VyLeft2))*yStep;
		//double yRight1 = yStart + iyHalfMaxRight*yStep;
		//double yRight = yRight1 + ((VyRight1 - HalfMaxVal)/(VyRight1 - VyRight2))*yStep;
		//double yFWHM = yRight - yLeft;

		double xFWFM = 0., yFWFM = 0.;
		int nRes = 5;
		if(nPar > 2)
		{//OC03012019
			double xFract = arPar[1], yFract = arPar[2];
			if((xFract == yFract) && (xFract > 0.))
			{
				HalfMaxVal = xFract*MaxVal;
				FindIndHalfMaxLeftRight2D(p0, xLenArr, yLenArr, HalfMaxVal, ixMax, iyMax, calcSizeFromPeak, arResInd, arResVxVy);
				AuxCalcWidth2D(xStart, xStep, yStart, yStep, HalfMaxVal, arResInd, arResVxVy, xFWFM, yFWFM);
			}
			else 
			{
				if(xFract > 0.)
				{
					HalfMaxVal = xFract*MaxVal;
					FindIndHalfMaxLeftRight2D(p0, xLenArr, yLenArr, HalfMaxVal, ixMax, iyMax, calcSizeFromPeak, arResInd, arResVxVy);
					double yAuxFWFM;
					AuxCalcWidth2D(xStart, xStep, yStart, yStep, HalfMaxVal, arResInd, arResVxVy, xFWFM, yAuxFWFM);
				}
				if(yFract > 0.)
				{
					HalfMaxVal = yFract*MaxVal;
					FindIndHalfMaxLeftRight2D(p0, xLenArr, yLenArr, HalfMaxVal, ixMax, iyMax, calcSizeFromPeak, arResInd, arResVxVy);
					double xAuxFWFM;
					AuxCalcWidth2D(xStart, xStep, yStart, yStep, HalfMaxVal, arResInd, arResVxVy, xAuxFWFM, yFWFM);
				}
			}
			nRes = 7;
		}

		//float* pOut = (float*)(OutSpotInfo.pWaveData); //make it more general (allow double and float)
		//*(pOut++) = (float)MaxVal;
		//*(pOut++) = (float)xMax;
		//*(pOut++) = (float)yMax;
		//*(pOut++) = (float)(xRight - xLeft);
		//*pOut = (float)(yRight - yLeft);

		double arRes[] = {MaxVal, xMax, yMax, xFWHM, yFWHM, xFWFM, yFWFM};
		if(*(OutSpotInfo.WaveType) == 'f') FillOutResData((float*)(OutSpotInfo.pWaveData), arRes, nRes); //OC03012019
		else FillOutResData((double*)(OutSpotInfo.pWaveData), arRes, nRes);

		//if(*(OutSpotInfo.WaveType) == 'f') FillOutResData((float*)(OutSpotInfo.pWaveData), MaxVal, xMax, yMax, xRight - xLeft, yRight - yLeft);
		//else FillOutResData((double*)(OutSpotInfo.pWaveData), MaxVal, xMax, yMax, xRight - xLeft, yRight - yLeft);
	}

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
