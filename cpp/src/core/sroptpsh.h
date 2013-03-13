/************************************************************************//**
 * File: sroptpsh.h
 * Description: Optical element: Phase Shift (obsolete?) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTPSH_H
#define __SROPTPSH_H

#include "sroptfoc.h"

//*************************************************************************

class srTPhaseShift : public srTFocusingElem {

	DOUBLE* tPhaseShiftData;

public:
	int FunNo;
	char IsFocusing;

	srTPhaseShift(srTStringVect* pElemInfo, srTSRWRadStructAccessData* pRad) 
	{
		IsFocusing = 0;
		FunNo = atoi((*pElemInfo)[1]);

		tPhaseShiftData = 0;
		if(ErrorCode = TryToSetUpFocalDist(pRad)) return;
	}

	int TryToSetUpFocalDist(srTSRWRadStructAccessData* pRad)
	{// Sets Up:
		int result;
		srTRadSect1D PointSourceSect1D[2];
		if(result = SetUpPointSourceSect1D(PointSourceSect1D, pRad)) return result;
		for(int k=0; k<2; k++)
		{
			if(result = SetUpFocalDistByPropag1D(PointSourceSect1D[k])) return result;
		}
		return 0;
	}

	int SetUpPointSourceSect1D(srTRadSect1D* PointSourceSect1D, srTSRWRadStructAccessData* pRad);
	int SetUpFocalDistByPropag1D(srTRadSect1D&);
	int SetUpPhaseShiftWave1D(srTRadSect1D& Sect1D, srTWaveAccessData& PhShData1D);
	int SetUpPhaseShiftWave(srTSRWRadStructAccessData& RadData, srTWaveAccessData& PhShData);

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
		srTWaveAccessData PhShWaveAccessData;
		if(result = SetUpPhaseShiftWave(*pRadAccessData, PhShWaveAccessData)) return result;
		tPhaseShiftData = (DOUBLE*)(PhShWaveAccessData.pWaveData);

		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;

		//srTSend Send;
		//if(result = Send.FinishWorkingWithWave(&PhShWaveAccessData)) return result;
		tPhaseShiftData = 0;
		return 0;
	}
  	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		int result;
		srTWaveAccessData PhShWaveAccessData1D;
		if(result = SetUpPhaseShiftWave1D(*pSect1D, PhShWaveAccessData1D)) return result;
		tPhaseShiftData = (DOUBLE*)(PhShWaveAccessData1D.pWaveData);

		if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
		if(result = TraverseRad1D(pSect1D)) return result;

		//srTSend Send;
		//if(result = Send.FinishWorkingWithWave(&PhShWaveAccessData1D)) return result;
		tPhaseShiftData = 0;
		return 0;
	}
	
	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double TwoPi_d_Lambda_m = EXZ.e*5.0676816042E+06;
		double PhaseShift = TwoPi_d_Lambda_m*(*(tPhaseShiftData++));
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);

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
		double TwoPi_d_Lambda_m = EXZ.e*5.0676816042E+06;
		double PhaseShift = TwoPi_d_Lambda_m*(*(tPhaseShiftData++));
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);

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
