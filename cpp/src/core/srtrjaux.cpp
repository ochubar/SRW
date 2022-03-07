/************************************************************************//**
 * File: srtrjaux.cpp
 * Description: Auxiliary structures for Electron Trajectory calculation
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srtrjaux.h"

//*************************************************************************

int srTTrjArraysAux::AllocateArraysOnLevels(int StLevNo, int FiLevNo)
{
	if((StLevNo < 0) || (StLevNo > FiLevNo)) return 0;
	if((FiLevNo < 0) || (FiLevNo > LEN_POINTERS - 1)) return 0;
	
	//long TotAmOfPointsWithThis = TotAmOfPointsIncludingTheLevel(FiLevNo);
	long long TotAmOfPointsWithThis = TotAmOfPointsIncludingTheLevel(FiLevNo);
	char UseActualAlloc = (TotAmOfPointsWithThis*10 > LEN_SMALL_CONT);
	
	for(int iLev=StLevNo; iLev<=FiLevNo; iLev++)
	{
		DeallocateArraysOnLevel(iLev);
		//long NpOnLevel = AmOfPointsOnLevel(iLev);
		long long NpOnLevel = AmOfPointsOnLevel(iLev);
		
		double **tBxArrP = BxArrP+iLev, **tBtxArrP = BtxArrP+iLev, **tXArrP = XArrP+iLev, **tIntBtxE2ArrP = IntBtxE2ArrP+iLev, **tdBxdsArrP = dBxdsArrP+iLev;
		double **tBzArrP = BzArrP+iLev, **tBtzArrP = BtzArrP+iLev, **tZArrP = ZArrP+iLev, **tIntBtzE2ArrP = IntBtzE2ArrP+iLev, **tdBzdsArrP = dBzdsArrP+iLev;
		
		if(UseActualAlloc)
		{
			*tBxArrP = new double[NpOnLevel*10];
			if(*tBxArrP == 0) { DeallocateArrays(); return MEMORY_ALLOCATION_FAILURE;}
			LevIsActuallyAlloc[iLev] = 1;
		}
		else
		{
			//long Offset = 10*TotAmOfPointsIncludingTheLevel(iLev - 1);
			long long Offset = 10*TotAmOfPointsIncludingTheLevel(iLev - 1);
			*tBxArrP = SmallDataCont + Offset;
		}
		
		*tBtxArrP = *tBxArrP + NpOnLevel;
		*tXArrP = *tBtxArrP + NpOnLevel;
		*tIntBtxE2ArrP = *tXArrP + NpOnLevel;
		*tdBxdsArrP = *tIntBtxE2ArrP + NpOnLevel;
		
		*tBzArrP = *tdBxdsArrP + NpOnLevel;
		*tBtzArrP = *tBzArrP + NpOnLevel;
		*tZArrP = *tBtzArrP + NpOnLevel;
		*tIntBtzE2ArrP = *tZArrP + NpOnLevel;
		*tdBzdsArrP = *tIntBtzE2ArrP + NpOnLevel;
	}
	return 0;
}

//*************************************************************************

int srTTrjArraysAux::AllocateArraysOnThisLevelNeglectingOthers(int LevNo)
{
	if((LevNo < 0) || (LevNo > LEN_POINTERS - 1)) return 0;

	//long NpOnLevel = AmOfPointsOnLevel(LevNo);
	long long NpOnLevel = AmOfPointsOnLevel(LevNo);
	char UseActualAlloc = (NpOnLevel*10 > LEN_SMALL_CONT);

	DeallocateArraysOnLevel(LevNo);

	double **tBxArrP = BxArrP+LevNo, **tBtxArrP = BtxArrP+LevNo, **tXArrP = XArrP+LevNo, **tIntBtxE2ArrP = IntBtxE2ArrP+LevNo, **tdBxdsArrP = dBxdsArrP+LevNo;
	double **tBzArrP = BzArrP+LevNo, **tBtzArrP = BtzArrP+LevNo, **tZArrP = ZArrP+LevNo, **tIntBtzE2ArrP = IntBtzE2ArrP+LevNo, **tdBzdsArrP = dBzdsArrP+LevNo;

	if(UseActualAlloc)
	{
		//long TotAmDat = NpOnLevel*10;
		long long TotAmDat = NpOnLevel*10;
		*tBxArrP = new double[TotAmDat];
		if(*tBxArrP == 0) { DeallocateArrays(); return MEMORY_ALLOCATION_FAILURE;}
		LevIsActuallyAlloc[LevNo] = 1;
	}
	else
	{
		//long Offset = 0;
		long long Offset = 0;
		*tBxArrP = SmallDataCont + Offset;
	}

	*tBtxArrP = *tBxArrP + NpOnLevel;
	*tXArrP = *tBtxArrP + NpOnLevel;
	*tIntBtxE2ArrP = *tXArrP + NpOnLevel;
	*tdBxdsArrP = *tIntBtxE2ArrP + NpOnLevel;
	
	*tBzArrP = *tdBxdsArrP + NpOnLevel;
	*tBtzArrP = *tBzArrP + NpOnLevel;
	*tZArrP = *tBtzArrP + NpOnLevel;
	*tIntBtzE2ArrP = *tZArrP + NpOnLevel;
	*tdBzdsArrP = *tIntBtzE2ArrP + NpOnLevel;
	return 0;
}

//*************************************************************************
