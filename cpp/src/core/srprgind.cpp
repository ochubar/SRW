/************************************************************************//**
 * File: srprgind.cpp
 * Description: Computation "Progress Indicator" (programmed and is proved to work only with IGOR Pro)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srprgind.h"

//*************************************************************************

//int srTCompProgressIndicator::InitializeIndicator(long InTotalAmOfOutPoints, double InUpdateTimeInt_s, char CountCallsInside)
int srTCompProgressIndicator::InitializeIndicator(long long InTotalAmOfOutPoints, double InUpdateTimeInt_s, char CountCallsInside)
{
	TotalAmOfOutPoints = InTotalAmOfOutPoints;

#ifdef __VC__
#ifdef __IGOR_PRO__

	int result;
	if(result = SetIgorIntVar("IgorTotAmOfPo", TotalAmOfOutPoints, 0)) return result;
	if(result = XOPCommand("NewPanel/W=(3.75,43.25,163,85)/K=2 as \" \"")) return result;
	//if(result = XOPCommand("String IgorIndicatorPanelWinName = WinName(0,64)")) return result;
	if(result = XOPCommand("String/G IgorIndicatorPanelWinName = WinName(0,64)")) return result;
	if(result = XOPCommand("SetDrawEnv fname= \"Helvetica\",fsize= 12; DrawText 15,18,\"Computation Progress\"; DoUpdate;")) return result;
	if(result = FetchStrVar("IgorIndicatorPanelWinName", ProgressIndicatorWinName)) return result;
	if(result = XOPCommand("ValDisplay srIgorCompProgressBar,pos={15,23},size={130,10},font=\"Helvetica\",limits={0,IgorTotAmOfPo,0},barmisc={30,0},value= #\"0\"")) return result;

#endif
#endif

#ifdef _SRWDLL
	if(pgCompProgressIndicFunc != 0) (*pgCompProgressIndicFunc)(0.);
#endif

	UpdateTimeInt = (clock_t)(CLOCKS_PER_SEC*InUpdateTimeInt_s);

	StartCompClock = clock();
	PrevUpdateClock = StartCompClock;

	ProgressIndicatorIsUsed = 1;
	PrevAmOfPoints = PrevAmOfPointsShown = 0;

	CallsCount = 0;
    CallsAreCountedInside = 0;
	if(CountCallsInside) CallsAreCountedInside = 1;

	return 0;
}

//*************************************************************************

