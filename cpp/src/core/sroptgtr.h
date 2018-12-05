/************************************************************************//**
 * File: sroptgtr.h
 * Description: Optical element: Transmission (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTGTR_H
#define __SROPTGTR_H

#include "sroptfoc.h"

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include <stdio.h>

struct SRWLStructOpticsTransmission;
typedef struct SRWLStructOpticsTransmission SRWLOptT;

//*************************************************************************

class srTGenTransmission : public srTFocusingElem {

	//srTWaveAccessData GenTransNumData;
	srTDataMD GenTransNumData;

	char OptPathOrPhase; // 1- Opt. Path, 2- Phase
	char OuterTransmIs; // 1- Zero, 2- Same as on border
	double eMid;
	double DxContin, DzContin; // Minimal intervals between discontinuties

public:

	srTGenTransmission(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	srTGenTransmission(const SRWLOptT& tr);
	~srTGenTransmission()
	{
		if(GenTransNumData.pData != 0)
		{
			//srTSend Send; Send.FinishWorkingWithWave(&GenTransNumData);
			//?????????????????
		}
	}

	void EnsureTransmissionForField();
	double DetermineAppropriatePhotEnergyForFocDistTest(double Rx, double Rz);
	int EstimateFocalDistancesAndCheckSampling();
	int SetUpPointSourceSect1D(char x_or_z, double R, double RelOtherCoord, srTRadSect1D& PointSourceSect1D);
	int DetermineFocalDistByPropag1D(srTRadSect1D& Sect1D, double& F);
	int DetermineFocalDistDirectAndCheckSampling(char x_or_z, double RelOtherCoord, double& Fx);
	int ExtractNumStructSect1DAndCheckSampling(char x_or_z, double RelOtherCoord, srTRadSect1D& Sect1D, double*& PhaseCont, char& PhaseIsContinuous);
	void CopyNumStructValuesToSect1DAndCheckSampling(srTRadSect1D& Sect1D, double* PhaseCont, char& PhaseIsContinuous);
	void EstimateEffPointsRange(char x_or_z, long icOtherCoord, long& iFirst, long& iLast, double& ArgFirst, double& ArgLast);
	int EstimateMinimalContinuousIntervals();

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		return srTFocusingElem::PropagateRadMoments(pRadAccessData, MomRatArray);
	}
	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{
		return srTFocusingElem::PropagateWaveFrontRadius(pRadAccessData);
	}
	int PropagateWaveFrontRadius1D(srTRadSect1D* pSect1D) 
	{
		return srTFocusingElem::PropagateWaveFrontRadius1D(pSect1D);
	}
	int Propagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData) 
	{
		return srTFocusingElem::Propagate4x4PropMatr(pRadAccessData);
	}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterArr)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterArr)
	{
		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			pRadAccessData->CheckAndSubtractPhaseTermsLin(TransvCenPoint.x, TransvCenPoint.y);
		//}

		char &MethNo = ParPrecWfrPropag.MethNo;
		
		int result = 0;

		if(MethNo == 0) result = PropagateRadiationMeth_0(pRadAccessData);
		else result = PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterArr);
		
		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			if(!ParPrecWfrPropag.DoNotResetAnalTreatTermsAfterProp) pRadAccessData->CheckAndResetPhaseTermsLin();
		//}

		return result;
	}

	//int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadDataSingleE)
	{
		int result;
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = PropagateRadiationSimple(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
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

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs);
  	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs);

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz);

	int RangeShouldBeAdjustedAtPropag() { return 0;} // Or switch it On
	int ResolutionShouldBeAdjustedAtPropag() { return 1;}
};

//*************************************************************************

#endif
