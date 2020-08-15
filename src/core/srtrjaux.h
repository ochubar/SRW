/************************************************************************//**
 * File: srtrjaux.cpp
 * Description: Auxiliary structures for Electron Trajectory calculation (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRTRJAUX_H
#define __SRTRJAUX_H

#include "srtrjdat.h"

//*************************************************************************

#define LEN_SMALL_CONT 10000
#define LEN_POINTERS 100

//*************************************************************************

class srTTrjArraysAux {

	double SmallDataCont[LEN_SMALL_CONT];

public:

	double *BxArrP[LEN_POINTERS], *BtxArrP[LEN_POINTERS], *XArrP[LEN_POINTERS], *IntBtxE2ArrP[LEN_POINTERS], *dBxdsArrP[LEN_POINTERS];
	double *BzArrP[LEN_POINTERS], *BtzArrP[LEN_POINTERS], *ZArrP[LEN_POINTERS], *IntBtzE2ArrP[LEN_POINTERS], *dBzdsArrP[LEN_POINTERS];
	char LevIsActuallyAlloc[LEN_POINTERS], LevelIsFilled[LEN_POINTERS];

	//int AmOfPointsOnZeroLevel;
	long long AmOfPointsOnZeroLevel;

	srTTrjArraysAux(int InAmOfPointsOnZeroLevel =5)
	{
		Initialize(InAmOfPointsOnZeroLevel);
	}
	~srTTrjArraysAux()
	{
		DeallocateArrays();
	}

	//void Initialize(int InAmOfPointsOnZeroLevel)
	void Initialize(long long InAmOfPointsOnZeroLevel)
	{
		AmOfPointsOnZeroLevel = InAmOfPointsOnZeroLevel;
		for(int k=0; k<LEN_POINTERS; k++) 
		{
			LevIsActuallyAlloc[k] = 0; LevelIsFilled[k] = 0;
			BxArrP[k] = 0;
		}
	}
	//void Reset(int InAmOfPointsOnZeroLevel)
	void Reset(long long InAmOfPointsOnZeroLevel)
	{
		DeallocateArrays();
		Initialize(InAmOfPointsOnZeroLevel);
	}

	int AllocateArraysOnLevels(int StLevNo, int FiLevNo);
	int AllocateArraysOnThisLevelNeglectingOthers(int LevNo);

	void DeallocateArraysOnLevel(int LevNo)
	{
		if((LevNo < 0) || (LevNo > LEN_POINTERS - 1)) return;
		double **tBxArrP = BxArrP+LevNo;
		if(LevIsActuallyAlloc[LevNo])
		{
			if(*tBxArrP != 0) delete[] (*tBxArrP);
		}

		*tBxArrP = 0;
		LevIsActuallyAlloc[LevNo] = 0;
		LevelIsFilled[LevNo] = 0;
	}
	void DeallocateArrays()
	{
		for(int k=0; k<LEN_POINTERS; k++) DeallocateArraysOnLevel(k);
	}

	int FillLevel(int LevNo, double sStart, double sEnd, srTGenTrjHndl& TrjHndl)
	{
		if((LevNo < 0) || (LevNo > LEN_POINTERS - 1)) return MEMORY_ALLOCATION_FAILURE;
		if(LevelIsFilled[LevNo]) return 0;

		int result;
		if(result = AllocateArraysOnLevels(LevNo, LevNo)) return result;

		//long Np = AmOfPointsOnLevel(LevNo);
		long long Np = AmOfPointsOnLevel(LevNo);
		TrjHndl.rep->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevNo], BtzArrP[LevNo], XArrP[LevNo], ZArrP[LevNo], IntBtxE2ArrP[LevNo], IntBtzE2ArrP[LevNo], BxArrP[LevNo], BzArrP[LevNo], dBxdsArrP[LevNo], dBzdsArrP[LevNo]);
		LevelIsFilled[LevNo] = 1;
		return 0;
	}

	int FillNextAndDestroyPrevLevel(int LevNo, double sStart, double sEnd, srTTrjDat* TrjDatPtr)
	{
		int LevNo_mi_1 = LevNo - 1;
		DeallocateArraysOnLevel(LevNo_mi_1);

		int result;
		if(result = AllocateArraysOnThisLevelNeglectingOthers(LevNo)) return result;

		//long Np = AmOfPointsOnLevel(LevNo);
		long long Np = AmOfPointsOnLevel(LevNo);
		TrjDatPtr->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevNo], BtzArrP[LevNo], XArrP[LevNo], ZArrP[LevNo], IntBtxE2ArrP[LevNo], IntBtzE2ArrP[LevNo], BxArrP[LevNo], BzArrP[LevNo], dBxdsArrP[LevNo], dBzdsArrP[LevNo]);
		LevelIsFilled[LevNo] = 1;
		return 0;
	}

	//long AmOfPointsOnLevel(int LevNo)
	long long AmOfPointsOnLevel(int LevNo)
	{
		return (LevNo == 0)? AmOfPointsOnZeroLevel : ((AmOfPointsOnZeroLevel - 1) << (LevNo - 1));
		//return (LevNo == 0)? AmOfPointsOnZeroLevel : ((AmOfPointsOnZeroLevel - 1) << ((long long)(LevNo - 1)));
	}
	//long TotAmOfPointsIncludingTheLevel(int LevNo)
	long long TotAmOfPointsIncludingTheLevel(int LevNo)
	{
		return (LevNo < 0)? 0 : (1 + ((AmOfPointsOnZeroLevel - 1) << LevNo));
	}

	void SetupPtrs_Btx_X_IntBtxE2_Bx_Btz_Z_IntBtzE2_Bz(int LevNo, double*** TrjPtrs)
	{
		*(TrjPtrs[0]) = BtxArrP[LevNo]; *(TrjPtrs[1]) = XArrP[LevNo]; *(TrjPtrs[2]) = IntBtxE2ArrP[LevNo]; *(TrjPtrs[3]) = BxArrP[LevNo];
		*(TrjPtrs[4]) = BtzArrP[LevNo]; *(TrjPtrs[5]) = ZArrP[LevNo]; *(TrjPtrs[6]) = IntBtzE2ArrP[LevNo]; *(TrjPtrs[7]) = BzArrP[LevNo];
	}
};

//*************************************************************************

//#undef LEN_SMALL_CONT
//#undef LEN_POINTERS

//*************************************************************************

#endif
