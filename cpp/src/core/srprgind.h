/************************************************************************//**
 * File: srprgind.h
 * Description: Computation "Progress Indicator" (programmed and is proved to work only with IGOR Pro) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRPRGIND_H
#define __SRPRGIND_H

#include <time.h>

#ifdef __IGOR_PRO__
#include "srigintr.h"
#endif

#ifdef _SRWDLL
extern int (*pgCompProgressIndicFunc)(double CurVal);
#endif

//*************************************************************************

class srTCompProgressIndicator {

#ifdef __IGOR_PRO__
	char ProgressIndicatorWinName[256];
#endif

	char ProgressIndicatorIsUsed, CallsAreCountedInside;
	//long TotalAmOfOutPoints, PrevAmOfPoints, PrevAmOfPointsShown, CallsCount;
	long long TotalAmOfOutPoints, PrevAmOfPoints, PrevAmOfPointsShown, CallsCount;
	clock_t UpdateTimeInt, PrevUpdateClock, StartCompClock;

public:

	char ErrorCode;

	//srTCompProgressIndicator(long InTotalAmOfOutPoints, double UpdateTimeInt_s, char CountCallsInside=0)
	srTCompProgressIndicator(long long InTotalAmOfOutPoints, double UpdateTimeInt_s, char CountCallsInside=0)
	{
		ProgressIndicatorIsUsed = 0; ErrorCode = 0;
		if(InTotalAmOfOutPoints <= 0) return;

		ErrorCode = InitializeIndicator(InTotalAmOfOutPoints, UpdateTimeInt_s, CountCallsInside);
	}
	srTCompProgressIndicator()
	{
		ProgressIndicatorIsUsed = 0; ErrorCode = 0;
	}
	~srTCompProgressIndicator()
	{
		DestroyIndicator();
	}

	//int InitializeIndicator(long InTotalAmOfOutPoints, double UpdateTimeInt_s, char CountCallsInside=0);
	int InitializeIndicator(long long InTotalAmOfOutPoints, double UpdateTimeInt_s, char CountCallsInside=0);

	//int UpdateIndicator(long CurPoint=0)
	int UpdateIndicator(long long CurPoint=0)
	{
		int result = 0;
		if((!ProgressIndicatorIsUsed) || (ErrorCode != 0)) return 0;

		if(CallsAreCountedInside) CurPoint = CallsCount;
		
		clock_t CurrentClock = clock();
		if(CurrentClock < (UpdateTimeInt + PrevUpdateClock)) return 0;

		//CurPoint = EstimateCurrentPosition(CurrentClock, CurPoint);
		//PrevAmOfPoints = CurPoint;
		PrevUpdateClock = CurrentClock;
		
#ifdef __IGOR_PRO__
#ifdef __VC__
		char TotOutStr[200];
		char PtNoStr[12];
		sprintf(PtNoStr, "%d\n", CurPoint);
		strcpy(TotOutStr, "ValDisplay srIgorCompProgressBar value= #\"");
		strcat(TotOutStr, PtNoStr);
		strcat(TotOutStr, "\"");
		strcat(TotOutStr, ", win=");
		strcat(TotOutStr, ProgressIndicatorWinName);
		if(result = XOPCommand(TotOutStr)) return result;
#endif
		if(result = SpinProcess()) return result;
#endif

#ifdef _SRWDLL
		if(pgCompProgressIndicFunc != 0) 
		{
			double NormIndVal = ((double)CurPoint)/((double)TotalAmOfOutPoints);
			(*pgCompProgressIndicFunc)(NormIndVal);
		}
#endif

		if(CallsAreCountedInside) CallsCount++;
		return result;
	}
	
	void DestroyIndicator()
	{
		if(!ProgressIndicatorIsUsed) return;

#ifdef __VC__
#ifdef __IGOR_PRO__
		char TotOutStr[300];
		strcpy(TotOutStr, "DoWindow/K ");
		strcat(TotOutStr, ProgressIndicatorWinName);
		XOPCommand(TotOutStr);
#endif
#endif

#ifdef _SRWDLL
		if(pgCompProgressIndicFunc != 0) (*pgCompProgressIndicFunc)(1.);
#endif
	}

	//long EstimateCurrentPosition(clock_t CurrentClock, long CurPoint)
	long long EstimateCurrentPosition(clock_t CurrentClock, long long CurPoint)
	{
		clock_t TimePassed = CurrentClock - StartCompClock;

		//double InvSpeed = double(CurrentClock - PrevUpdateClock)/double(CurPoint - PrevAmOfPoints);
		double InvSpeed = double(TimePassed)/double(CurPoint);

		clock_t EstTimeLeft = (clock_t)(InvSpeed*double(TotalAmOfOutPoints - CurPoint));

		//long EstCurPoint = long(TotalAmOfOutPoints*double(TimePassed)/double(EstTimeLeft + TimePassed));
		long long EstCurPoint = (long long)(TotalAmOfOutPoints*double(TimePassed)/double(EstTimeLeft + TimePassed));
		
		PrevAmOfPoints = CurPoint;
		PrevUpdateClock = CurrentClock;

		if(PrevAmOfPointsShown < EstCurPoint) PrevAmOfPointsShown = EstCurPoint;
		return PrevAmOfPointsShown;
	}
};

//*************************************************************************

#endif
