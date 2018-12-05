/************************************************************************//**
 * File: sroptel2.cpp
 * Description: Optical element (general functions)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptelm.h"
#include "srsysuti.h"
#include "srmlttsk.h"
#include "srerror.h"

//OC31102018: added by SY at parallelizing SRW via OpenMP
//#include "srwlib.h"
//#include "stdio.h"

#ifdef _WITH_OMP //OC31102018: added by SY at parallelizing SRW via OpenMP
#include "omp.h"
#endif

extern srTYield srYield;

//*************************************************************************

double srTGenOptElem::CheckMemoryAvailable()
{
	return srTSystemUtils::CheckMemoryAvailable();
}

//*************************************************************************

int srTGenOptElem::PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData) 
{//Moved from derived classes: loops over E, calls derived PropagateRadiationSingleE_Meth_0
 //This propagation method doesn't allow for true wavefront "resizing/resampling" 
 //(which results in changing numbers of points) in "slices" vs photon energy.
 //However, it may allow for changing wavefront step size and range (keeping numbers of points constant)
 //in the "slices".
 //It is virtual ("standard" processing). Known re-definitions: in srTDriftSpace

	//OC31102018: added by SY at parallelizing SRW via OpenMP
	//double start;
	//get_walltime (&start);

	int result=0;

#ifndef _WITH_OMP //OC31102018

	srTSRWRadStructAccessData *pRadDataSingleE = 0, *pPrevRadDataSingleE = 0;
	if(pRadAccessData->ne == 1)
	{
		pRadDataSingleE = pRadAccessData;
	}
	else
	{
		if(result = SetupNewRadStructFromSliceConstE(pRadAccessData, -1, pRadDataSingleE)) return result;
		//allocates new pRadDataSingleE !
	}

	if(!m_PropWfrInPlace)
	{
		if(result = SetupNewRadStructFromSliceConstE(pRadAccessData, -1, pPrevRadDataSingleE)) return result;
	}

	//separate processing of wavefront radius is necessary
	double origRobsX = pRadAccessData->RobsX, origRobsXAbsErr = pRadAccessData->RobsXAbsErr;
	double origRobsZ = pRadAccessData->RobsZ, origRobsZAbsErr = pRadAccessData->RobsZAbsErr;
	double origXc = pRadAccessData->xc; //OC180813
	double origZc = pRadAccessData->zc;

	const int AmOfMoments = 11;
	int neOrig = pRadAccessData->ne;

	vector<srTSRWRadStructAccessData> vRadSlices; //to keep grid parameters of slices for eventual change of the main 3D grid and re-interpolation at the end
	bool gridParamWereModifInSlices = false;

	for(int ie=0; ie<neOrig; ie++)
	{
		if(pRadDataSingleE != pRadAccessData)
		{
			if(result = ExtractRadSliceConstE(pRadAccessData, ie, pRadDataSingleE->pBaseRadX, pRadDataSingleE->pBaseRadZ)) return result;
			pRadDataSingleE->eStart = pRadAccessData->eStart + ie*pRadAccessData->eStep;
			long OffsetMom = AmOfMoments*ie;
			pRadDataSingleE->pMomX = pRadAccessData->pMomX + OffsetMom;
			pRadDataSingleE->pMomZ = pRadAccessData->pMomZ + OffsetMom;
			pRadDataSingleE->RobsX = origRobsX; pRadDataSingleE->RobsXAbsErr = origRobsXAbsErr;
			pRadDataSingleE->RobsZ = origRobsZ; pRadDataSingleE->RobsZAbsErr = origRobsZAbsErr;

			pRadDataSingleE->xc = origXc; //OC180813
			pRadDataSingleE->zc = origZc;

			//if(pRadAccessData->xStart != pRadDataSingleE->xStart) pRadDataSingleE->xStart = pRadAccessData->xStart;
			//if(pRadAccessData->xStep != pRadDataSingleE->xStep) pRadDataSingleE->xStep = pRadAccessData->xStep;
			//if(pRadAccessData->zStart != pRadDataSingleE->zStart) pRadDataSingleE->zStart = pRadAccessData->zStart;
			//if(pRadAccessData->zStep != pRadDataSingleE->zStep) pRadDataSingleE->zStep = pRadAccessData->zStep;
			pRadDataSingleE->xStart = pRadAccessData->xStart;
			pRadDataSingleE->xStep = pRadAccessData->xStep;
			pRadDataSingleE->zStart = pRadAccessData->zStart;
			pRadDataSingleE->zStep = pRadAccessData->zStep;
		}
		if(pPrevRadDataSingleE != 0)
		{
			if(result = ExtractRadSliceConstE(pRadAccessData, ie, pPrevRadDataSingleE->pBaseRadX, pPrevRadDataSingleE->pBaseRadZ, true)) return result; //OC120908
			pPrevRadDataSingleE->eStart = pRadDataSingleE->eStart;
			pPrevRadDataSingleE->pMomX = pRadDataSingleE->pMomX;
			pPrevRadDataSingleE->pMomZ = pRadDataSingleE->pMomZ;
		}

		if(result = PropagateRadiationSingleE_Meth_0(pRadDataSingleE, pPrevRadDataSingleE)) return result; //from derived classes

		if(pRadDataSingleE != pRadAccessData)
		{
			if(result = UpdateGenRadStructSliceConstE_Meth_0(pRadDataSingleE, ie, pRadAccessData)) return result;
			//the above doesn't change the transverse grid parameters in *pRadAccessData

			//vRadSlices.push_back(*pRadDataSingleE); //this automatically calls destructor, which can eventually delete "emulated" structs!
			//srTSRWRadStructAccessData copyRadDataSingleE(*pRadDataSingleE); //this doesn't assume to copy pBaseRadX, pBaseRadZ
			srTSRWRadStructAccessData copyRadDataSingleE(*pRadDataSingleE, false); //OC290813 fixing memory leak(?) //this doesn't assume to copy pBaseRadX, pBaseRadZ
			copyRadDataSingleE.pBaseRadX = copyRadDataSingleE.pBaseRadZ = 0; copyRadDataSingleE.BaseRadWasEmulated = false;
			vRadSlices.push_back(copyRadDataSingleE);
			
			if((pRadDataSingleE->nx != pRadAccessData->nx) || (pRadDataSingleE->xStart != pRadAccessData->xStart) || (pRadDataSingleE->xStep != pRadAccessData->xStep)) gridParamWereModifInSlices = true;
			if((pRadDataSingleE->nz != pRadAccessData->nz) || (pRadDataSingleE->zStart != pRadAccessData->zStart) || (pRadDataSingleE->zStep != pRadAccessData->zStep)) gridParamWereModifInSlices = true;
		}
	}

	if((pRadAccessData->RobsX != 0) && (pRadAccessData->RobsXAbsErr == 0)) pRadAccessData->RobsXAbsErr = ::fabs(0.1*pRadAccessData->RobsX);
	if((pRadAccessData->RobsZ != 0) && (pRadAccessData->RobsZAbsErr == 0)) pRadAccessData->RobsZAbsErr = ::fabs(0.1*pRadAccessData->RobsZ);

	if(gridParamWereModifInSlices)
	{//to test!
		if(result = ReInterpolateWfrDataOnNewTransvMesh(vRadSlices, pRadDataSingleE, pRadAccessData)) return result;
	}

	if((pRadDataSingleE != 0) && (pRadDataSingleE != pRadAccessData)) delete pRadDataSingleE;
	if((pPrevRadDataSingleE != 0) && (pPrevRadDataSingleE != pRadAccessData)) delete pPrevRadDataSingleE;

#else //OC31102018: modified by SY at parallelizing SRW via OpenMP

	if(pRadAccessData->ne == 1)
	{
		//SY: nothing more is actually needed in this case
		return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0); //from derived classes

		//OC31102018: added by SY (for profiling?) at parallelizing SRW via OpenMP
		//srwlPrintTime(": PropagateRadiationMeth_0 : PropagateRadiationSingleE_Meth_0 - single",&start);
	}

	srTSRWRadStructAccessData *pPrevRadDataSingleE = 0;

	if(!m_PropWfrInPlace)
	{
		if(result = SetupNewRadStructFromSliceConstE(pRadAccessData, -1, pPrevRadDataSingleE)) return result;
	}

	//separate processing of wavefront radius is necessary
	double origRobsX = pRadAccessData->RobsX, origRobsXAbsErr = pRadAccessData->RobsXAbsErr;
	double origRobsZ = pRadAccessData->RobsZ, origRobsZAbsErr = pRadAccessData->RobsZAbsErr;
	double origXc = pRadAccessData->xc; //OC180813
	double origZc = pRadAccessData->zc;

	const int AmOfMoments = 11;
	int neOrig = pRadAccessData->ne;

	vector<srTSRWRadStructAccessData> vRadSlices; //to keep grid parameters of slices for eventual change of the main 3D grid and re-interpolation at the end
	srTSRWRadStructAccessData** pRadSlices= new srTSRWRadStructAccessData* [neOrig];

	bool gridParamWereModifInSlices = false;
	bool* gridParamWereModif = new bool[neOrig];

	//SY: return outside of parallel regions is not allowed - we do it outside

	int* results = new int[neOrig];
	if(results == 0) return MEMORY_ALLOCATION_FAILURE;
	for(int ie = 0; ie < neOrig; ie++) results[ie] = 0;

	//OC31102018: added by SY (for profiling?) at parallelizing SRW via OpenMP
	//srwlPrintTime(": PropagateRadiationMeth_0 : before cycle",&start);

	//SY: we cannot do it in parallel if previous field is needed (which is not the case in the current version of SRW)
	#pragma omp parallel if (pPrevRadDataSingleE == 0)
	{
		int threadNum = omp_get_thread_num();
		srTSRWRadStructAccessData *pRadDataSingleE = 0;
		results[threadNum] = SetupNewRadStructFromSliceConstE(pRadAccessData, -1, pRadDataSingleE);
		//allocates new pRadDataSingleE !
		if(!results[threadNum])
		{
			#pragma omp for
			for(int ie=0; ie<neOrig; ie++)
			{
				gridParamWereModif[ie] = false;
				if(results[ie] = ExtractRadSliceConstE(pRadAccessData, ie, pRadDataSingleE->pBaseRadX, pRadDataSingleE->pBaseRadZ)) continue;
				pRadDataSingleE->eStart = pRadAccessData->eStart + ie*pRadAccessData->eStep;
				long OffsetMom = AmOfMoments*ie;
				pRadDataSingleE->pMomX = pRadAccessData->pMomX + OffsetMom;
				pRadDataSingleE->pMomZ = pRadAccessData->pMomZ + OffsetMom;
				pRadDataSingleE->RobsX = origRobsX; pRadDataSingleE->RobsXAbsErr = origRobsXAbsErr;
				pRadDataSingleE->RobsZ = origRobsZ; pRadDataSingleE->RobsZAbsErr = origRobsZAbsErr;

				pRadDataSingleE->xc = origXc; //OC180813
				pRadDataSingleE->zc = origZc;

				//if(pRadAccessData->xStart != pRadDataSingleE->xStart) pRadDataSingleE->xStart = pRadAccessData->xStart;
				//if(pRadAccessData->xStep != pRadDataSingleE->xStep) pRadDataSingleE->xStep = pRadAccessData->xStep;
				//if(pRadAccessData->zStart != pRadDataSingleE->zStart) pRadDataSingleE->zStart = pRadAccessData->zStart;
				//if(pRadAccessData->zStep != pRadDataSingleE->zStep) pRadDataSingleE->zStep = pRadAccessData->zStep;
				pRadDataSingleE->xStart = pRadAccessData->xStart;
				pRadDataSingleE->xStep = pRadAccessData->xStep;
				pRadDataSingleE->zStart = pRadAccessData->zStart;
				pRadDataSingleE->zStep = pRadAccessData->zStep;

				if(pPrevRadDataSingleE != 0)
				{
					if(results[ie] = ExtractRadSliceConstE(pRadAccessData, ie, pPrevRadDataSingleE->pBaseRadX, pPrevRadDataSingleE->pBaseRadZ, true)) continue; //OC120908
					pPrevRadDataSingleE->eStart = pRadDataSingleE->eStart;
					pPrevRadDataSingleE->pMomX = pRadDataSingleE->pMomX;
					pPrevRadDataSingleE->pMomZ = pRadDataSingleE->pMomZ;
				}
				if(results[ie] = PropagateRadiationSingleE_Meth_0(pRadDataSingleE, pPrevRadDataSingleE)) continue; //from derived classes

				if(results[ie] = UpdateGenRadStructSliceConstE_Meth_0(pRadDataSingleE, ie, pRadAccessData,1)) continue;
				//the above doesn't change the transverse grid parameters in *pRadAccessData

				//vRadSlices.push_back(*pRadDataSingleE); //this automatically calls destructor, which can eventually delete "emulated" structs!
				//srTSRWRadStructAccessData copyRadDataSingleE(*pRadDataSingleE); //this doesn't assume to copy pBaseRadX, pBaseRadZ
				//srTSRWRadStructAccessData copyRadDataSingleE(*pRadDataSingleE, false); //OC290813 fixing memory leak(?) //this doesn't assume to copy pBaseRadX, pBaseRadZ
				//SY: we cannot use push here, so use normal array instead (probably using of vector would be still possible, but just to be sure)

				pRadSlices[ie] = new srTSRWRadStructAccessData(*pRadDataSingleE, false);
				pRadSlices[ie]->pBaseRadX = pRadSlices[ie]->pBaseRadZ = 0; pRadSlices[ie]->BaseRadWasEmulated = false;

				if((pRadDataSingleE->nx != pRadAccessData->nx) || (pRadDataSingleE->xStart != pRadAccessData->xStart) || (pRadDataSingleE->xStep != pRadAccessData->xStep)) gridParamWereModif[ie] = true;
				if((pRadDataSingleE->nz != pRadAccessData->nz) || (pRadDataSingleE->zStart != pRadAccessData->zStart) || (pRadDataSingleE->zStep != pRadAccessData->zStep)) gridParamWereModif[ie] = true;
			}
		}
		if(pRadDataSingleE != 0) delete pRadDataSingleE;
	}

	for(int ie = 0; ie < neOrig; ie++) if(results[ie]) return results[ie];
	delete[]  results;

	//OC31102018: added by SY (for profiling?) at parallelizing SRW via OpenMP
	//char str[256];
	//sprintf(str,"%s %d",":PropagateRadiationMeth_0 : PropagateRadiationSingleE_Meth_0 - cycles:",neOrig);
	//srwlPrintTime(str,&start);

	//SY: update limits and Radii (cannot be done in parallel)
	for(int ie = 0; ie < neOrig; ie++) if(result = UpdateGenRadStructSliceConstE_Meth_0(pRadSlices[ie], ie, pRadAccessData,2)) return result;

	//OC31102018: added by SY (for profiling?) at parallelizing SRW via OpenMP
	//srwlPrintTime(":PropagateRadiationMeth_0 : UpdateGenRadStructSliceConstE_Meth_0:",&start);

	if((pRadAccessData->RobsX != 0) && (pRadAccessData->RobsXAbsErr == 0)) pRadAccessData->RobsXAbsErr = ::fabs(0.1*pRadAccessData->RobsX);
	if((pRadAccessData->RobsZ != 0) && (pRadAccessData->RobsZAbsErr == 0)) pRadAccessData->RobsZAbsErr = ::fabs(0.1*pRadAccessData->RobsZ);

	for(int i=0;i<neOrig;i++) if(gridParamWereModif[i]) gridParamWereModifInSlices = true;

	if(gridParamWereModifInSlices)
	{//to test!
		for(int ie = 0; ie < neOrig; ie++)
		{
			vRadSlices.push_back(*pRadSlices[ie]);
			delete pRadSlices[ie];
		}
		// SY: instead of using pRadDataSingleE (which is not accessible here) we create an auxiliary structure
		srTSRWRadStructAccessData *pAuxRadData = 0;
		if(result = SetupNewRadStructFromSliceConstE(pRadAccessData, -1, pAuxRadData)) return result;
		if(result = ReInterpolateWfrDataOnNewTransvMesh(vRadSlices, pAuxRadData, pRadAccessData)) return result;
		if(pAuxRadData != 0) delete pAuxRadData;
	}

	delete[] pRadSlices;
	delete[] gridParamWereModif;

	if(pPrevRadDataSingleE != 0) delete pPrevRadDataSingleE;

#endif

	return result;
}

//*************************************************************************

//int srTGenOptElem::PropagateRadiationMeth_2(srTSRWRadStructAccessData* pRadAccessData, srTRadResizeVect& ResizeBeforeAndAfterVect)
int srTGenOptElem::PropagateRadiationMeth_2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResizeBeforeAndAfterVect)
{
	srTSRWRadStructAccessData *pRadDataSingleE = 0; //, *RadDataStructArray = 0;
	if(pRadAccessData->ne == 1) pRadDataSingleE = pRadAccessData;
	//else
	//{
	//	RadDataStructArray = new srTSRWRadStructAccessData[pRadAccessData->ne];
	//	if(RadDataStructArray == 0) return MEMORY_ALLOCATION_FAILURE;
	//	pRadDataSingleE = RadDataStructArray;
	//}

	int result;
	for(int ie=0; ie<pRadAccessData->ne; ie++)
	{
		if(pRadDataSingleE != pRadAccessData)
		{
			if(result = SetupNewRadStructFromSliceConstE(pRadAccessData, ie, pRadDataSingleE)) return result;
			//if(result = RemoveSliceConstE_FromGenRadStruct(pRadAccessData, ie)) return result; // To save memory
		}
		
		//if(result = PropagateRadiationSingleE_Meth_2(pRadDataSingleE, ResizeBeforeAndAfterVect)) return result;
		if(result = PropagateRadiationSingleE_Meth_2(pRadDataSingleE, ParPrecWfrPropag, ResizeBeforeAndAfterVect)) return result;
		
		if(pRadDataSingleE != pRadAccessData)
		{
			//To implement !!! for structures with (eventually) different nx, nz !!!:
			if(result = UpdateGenRadStructSliceConstE_Meth_2(pRadDataSingleE, ie, pRadAccessData)) return result;
			delete pRadDataSingleE;
			pRadDataSingleE = 0;
		}
		//pRadDataSingleE++;
	}

	//if(pRadAccessData->ne > 1)
	//{
	//	if(result = UpdateGenRadStructFromSlicesConstE(pRadAccessData, RadDataStructArray)) return result;
	//}

	//if(RadDataStructArray != 0) delete[] RadDataStructArray;
	return 0;
}

//*************************************************************************

//int srTGenOptElem::PropagateRadiationSingleE_Meth_2(srTSRWRadStructAccessData* pRadAccessData, srTRadResizeVect& ActResizeBeforeAndAfterVect)
int srTGenOptElem::PropagateRadiationSingleE_Meth_2(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ActResizeBeforeAndAfterVect)
{
	int result;

// Ensure coordinate representation
	if(result = EnsureCoordRepres(pRadAccessData)) return result;

	srTPredictedPropagData1D PredictedPropagData[2];
	srTMomentsPtrs MomX(pRadAccessData->pMomX), MomZ(pRadAccessData->pMomZ);
	char TreatExOrEz = (*(MomX.pTotPhot) > *(MomZ.pTotPhot))? 'x' : 'z';
	if(TreatExOrEz == 'x')
	{
		PredictedPropagData[0].SigmaCoord = MomX.SqrtMxx;
		PredictedPropagData[0].SigmaAng = MomX.SqrtMxpxp;
		PredictedPropagData[1].SigmaCoord = MomX.SqrtMzz;
		PredictedPropagData[1].SigmaAng = MomX.SqrtMzpzp;
	}
	else
	{
		PredictedPropagData[0].SigmaCoord = MomZ.SqrtMxx;
		PredictedPropagData[0].SigmaAng = MomZ.SqrtMxpxp;
		PredictedPropagData[1].SigmaCoord = MomZ.SqrtMzz;
		PredictedPropagData[1].SigmaAng = MomZ.SqrtMzpzp;
	}

	srTMomentsRatios MomRatios;
	if(result = PropagateRadMoments(pRadAccessData, &MomRatios)) return result;
	MomRatios.GiveMaxRatios(PredictedPropagData[0].CoordMomRat, PredictedPropagData[0].AngMomRat, PredictedPropagData[1].CoordMomRat, PredictedPropagData[1].AngMomRat);

// Define Propagation Scenario
	char ThereWasUnderSamplingX = pRadAccessData->ThereIsUnderSamplingX();
	char ThereWasUnderSamplingZ = pRadAccessData->ThereIsUnderSamplingZ();

	srTPropagScenario1D PropagScenario[2];
	if(result = DefinePropagScenario(pRadAccessData, ParPrecWfrPropag, PredictedPropagData, PropagScenario)) return result;

	srTRadResize ResBefore, ResAfter;
	TransferResizeParam(PropagScenario, ResBefore, ResAfter);

		pRadAccessData->WfrEdgeCorrShouldBeDone = 0; //test !!!

// Resize before Propagation, if necessary

	if(ParPrecWfrPropag.UseResBefore)
	{
		if((ResBefore.pxd != 1.) && ThereWasUnderSamplingX)
		{
			pRadAccessData->UnderSamplingX *= ResBefore.pxd;
			ResBefore.pxd = 1.;
		}
		if((ResBefore.pzd != 1.) && ThereWasUnderSamplingZ)
		{
			pRadAccessData->UnderSamplingZ *= ResBefore.pzd;
			ResBefore.pzd = 1.;
		}
		if((ResBefore.pxm != 1.) || (ResBefore.pxd != 1.) || (ResBefore.pzm != 1.) || (ResBefore.pzd != 1.))
		{
			if(result = RadResizeGen(*pRadAccessData, ResBefore)) return result;
		}
	}
	else
	{
		ResBefore.pxd = ResBefore.pzd = ResBefore.pxm = ResBefore.pzm = 1.;
	}

	pRadAccessData->pResAfter = &ResAfter;

// Simple Propagation
	if(result = PropagateRadiationSimple(pRadAccessData)) return result;

// Propagate Wavefront Radius
	if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
	if(result = Propagate4x4PropMatr(pRadAccessData)) return result;

// Resize after Propagation, if necessary
	//if(!pRadAccessData->DoNotResizeAfter)
	if(ParPrecWfrPropag.UseResAfter)
	{
		if(result = MakePostPropagationProc(pRadAccessData, ResAfter)) return result;
	}
	else
	{
		ResAfter.pxd = ResAfter.pzd = ResAfter.pxm = ResAfter.pzm = 1.;
	}
	if(result = RecomputeRadMomentsIfPossible(pRadAccessData)) return result; // Make it at a condition

	ActResizeBeforeAndAfterVect.push_back(ResBefore);
	ActResizeBeforeAndAfterVect.push_back(ResAfter);

	return 0;
}

//*************************************************************************

void srTGenOptElem::FindWidestWfrMeshParam(vector<srTSRWRadStructAccessData>& vRadSlices, srTSRWRadStructAccessData* pRad, bool keepConstNumPoints)
{
	int numWfr = (int)vRadSlices.size();
	if((pRad == 0) || (numWfr <= 0)) return;

	srTSRWRadStructAccessData &vRadSlice_0 = vRadSlices[0];
	if(numWfr == 1)
	{
		pRad->xStart = vRadSlice_0.xStart; pRad->xStep = vRadSlice_0.xStep; pRad->nx = vRadSlice_0.nx;
		pRad->zStart = vRadSlice_0.zStart; pRad->zStep = vRadSlice_0.zStep; pRad->nz = vRadSlice_0.nz;
		return;
	}

	double res_xStart, res_xEnd, res_zStart, res_zEnd;
	int res_nx, res_nz;

	for(int i=0; i<numWfr; i++)
	{
		srTSRWRadStructAccessData &rad = vRadSlices[i];
		if(i == 0)
		{
			res_xStart = rad.xStart; res_nx = rad.nx; res_xEnd = res_xStart + rad.xStep*res_nx;
			res_zStart = rad.zStart; res_nz = rad.nz; res_zEnd = res_zStart + rad.zStep*res_nz;
			continue;
		}

		double &test_xStart = rad.xStart, &test_zStart = rad.zStart, &test_xStep = rad.xStep, &test_zStep = rad.zStep;
		long &test_nx = rad.nx, &test_nz = rad.nz;
		double test_xEnd = test_xStart + test_nx*test_xStep;
		double test_zEnd = test_zStart + test_nz*test_zStep;

		if(res_xStart > test_xStart) res_xStart = test_xStart;
		if(res_xEnd < test_xEnd) res_xEnd = test_xEnd;
		if(res_nx < test_nx) res_nx = test_nx;

		if(res_zStart > test_zStart) res_zStart = test_zStart;
		if(res_zEnd < test_zEnd) res_zEnd = test_zEnd;
		if(res_nz < test_nz) res_nz = test_nz;
	}

	if(!keepConstNumPoints) pRad->nx = res_nx;
	pRad->xStart = res_xStart;
	pRad->xStep = (pRad->nx > 0)? (res_xEnd - res_xStart)/(pRad->nx) : 0;

	if(!keepConstNumPoints) pRad->nz = res_nz;
	pRad->zStart = res_zStart;
	pRad->zStep = (pRad->nz > 0)? (res_zEnd - res_zStart)/(pRad->nz) : 0;
}

//*************************************************************************

int srTGenOptElem::ReInterpolateWfrDataOnNewTransvMesh(vector<srTSRWRadStructAccessData>& vRadSlices, srTSRWRadStructAccessData* pAuxRadSingleE, srTSRWRadStructAccessData* pRadRes)
{//this requires same nx, nz in all rad. structures
 //assumes that field data was allocated for pAuxRadSingleE, pRadRes;

	FindWidestWfrMeshParam(vRadSlices, pRadRes, true);

	int numSlicesE = (int)vRadSlices.size();
	if((numSlicesE <= 0) || (pAuxRadSingleE == 0) || (pRadRes == 0)) return 0;

	//if(pAuxRadSingleE->nx <= pRadRes->nx) return 0; //to fire error?
	if(pAuxRadSingleE->nx != pRadRes->nx) return 0; //to fire error?

	//double xAbsTol = 0.000001*pRadRes->xStep;
	//double zAbsTol = 0.000001*pRadRes->zStep;
	double xAbsTol = 0.0001*pRadRes->xStep; //OC271214
	double zAbsTol = 0.0001*pRadRes->zStep;

	int result = 0;

	float *pOrigBufEX = pAuxRadSingleE->pBaseRadX, *pOrigBufEZ = pAuxRadSingleE->pBaseRadZ;
	double origMultiRobsX = pRadRes->RobsX, origMultiErrRobsX = pRadRes->RobsXAbsErr;
	double origMultiRobsZ = pRadRes->RobsZ, origMultiErrRobsZ = pRadRes->RobsZAbsErr;

	for(int ie=0; ie<numSlicesE; ie++)
	{
		srTSRWRadStructAccessData &radMesh = vRadSlices[ie];

		if((radMesh.nx == pRadRes->nx) && (::fabs(radMesh.xStart - pRadRes->xStart) < xAbsTol) && (::fabs(radMesh.xStep - pRadRes->xStep) < xAbsTol) &&
		   (radMesh.nz == pRadRes->nz) && (::fabs(radMesh.zStart - pRadRes->zStart) < zAbsTol) && (::fabs(radMesh.zStep - pRadRes->zStep) < zAbsTol)) continue;

		if(result = ExtractRadSliceConstE(pRadRes, ie, pOrigBufEX, pOrigBufEZ, true)) return result;
		*pAuxRadSingleE = radMesh;
		pAuxRadSingleE->pBaseRadX = pOrigBufEX; pAuxRadSingleE->pBaseRadZ = pOrigBufEZ;
		//we require transverse mesh parameters and wfr radii of curvature and their errors!

		pRadRes->RobsX = radMesh.RobsX; pRadRes->RobsXAbsErr = radMesh.RobsXAbsErr;
		pRadRes->RobsZ = radMesh.RobsZ; pRadRes->RobsZAbsErr = radMesh.RobsZAbsErr;
		//to allow for removing / adding the quadratic phase term

		if(result = ReInterpolateWfrSliceSingleE(*pAuxRadSingleE, *pRadRes, ie)) return result;
	}
	pRadRes->RobsX = origMultiRobsX; pRadRes->RobsXAbsErr = origMultiErrRobsX;
	pRadRes->RobsZ = origMultiRobsZ; pRadRes->RobsZAbsErr = origMultiErrRobsZ;

	return 0;
}

//*************************************************************************

int srTGenOptElem::RecomputeRadMomentsIfPossible(srTSRWRadStructAccessData* pRadAccessData)
{// This assumes that moments were already propagated !
	int result;
	srTMomentsPtrs PropMomXPtrs(pRadAccessData->pMomX), PropMomZPtrs(pRadAccessData->pMomZ);
	srTMomentsVals PropMomX(PropMomXPtrs), PropMomZ(PropMomZPtrs);

	char BadConditionsForCoordComp = !WaveFrontTermCanBeTreated(*pRadAccessData);
	if(result = ComputeRadMoments(pRadAccessData)) return result;

	if(BadConditionsForCoordComp) // Setting back angular moments
	{
		*(PropMomXPtrs.pXP) = PropMomX.XP; *(PropMomZPtrs.pXP) = PropMomZ.XP;
		*(PropMomXPtrs.pZP) = PropMomX.ZP; *(PropMomZPtrs.pZP) = PropMomZ.ZP;
		*(PropMomXPtrs.pXXP) = PropMomX.XXP; *(PropMomZPtrs.pXXP) = PropMomZ.XXP;
		*(PropMomXPtrs.pXPXP) = PropMomX.XPXP; *(PropMomZPtrs.pXPXP) = PropMomZ.XPXP;
		*(PropMomXPtrs.pZZP) = PropMomX.ZZP; *(PropMomZPtrs.pZZP) = PropMomZ.ZZP;
		*(PropMomXPtrs.pZPZP) = PropMomX.ZPZP; *(PropMomZPtrs.pZPZP) = PropMomZ.ZPZP;
	}
	
	CheckAndCorrectSecondOrderRadAngMoments(pRadAccessData);
	return 0;
}

//*************************************************************************

void srTGenOptElem::CheckAndCorrectSecondOrderRadAngMoments(srTSRWRadStructAccessData* pRadAccessData)
{
	srTMomentsPtrs MomX(pRadAccessData->pMomX), MomZ(pRadAccessData->pMomZ);

	double Lambda_m_d_4Pi = 9.8664446e-08/pRadAccessData->eStart;
	double Lambda_m_d_4Pi_E2 = Lambda_m_d_4Pi*Lambda_m_d_4Pi;
	
	double MaxMxxpMomXe2 = MomX.Mxx*MomX.Mxpxp - Lambda_m_d_4Pi_E2;
	if(MaxMxxpMomXe2 < 0.)
	{
		*(MomX.pXXP) = (*(MomX.pX))*(*(MomX.pXP));
		*(MomX.pXPXP) = (float)(Lambda_m_d_4Pi_E2/(MomX.Mxx) + (*(MomX.pXP))*(*(MomX.pXP)));
	}
	//else //OC051206
	//{
	//	*(MomX.pXXP) = (float)((*(MomX.pX))*(*(MomX.pXP)) + sqrt(MaxMxxpMomXe2));
	//}
	double MaxMzzpMomXe2 = MomX.Mzz*MomX.Mzpzp - Lambda_m_d_4Pi_E2;
	if(MaxMzzpMomXe2 < 0.)
	{
		*(MomX.pZZP) = (*(MomX.pZ))*(*(MomX.pZP));
		*(MomX.pZPZP) = (float)(Lambda_m_d_4Pi_E2/(MomX.Mzz) + (*(MomX.pZP))*(*(MomX.pZP)));
	}
	//else
	//{
	//	*(MomX.pZZP) = (float)((*(MomX.pZ))*(*(MomX.pZP)) + sqrt(MaxMzzpMomXe2));
	//}

	double MaxMxxpMomZe2 = MomZ.Mxx*MomZ.Mxpxp - Lambda_m_d_4Pi_E2;
	if(MaxMxxpMomZe2 < 0.)
	{
		*(MomZ.pXXP) = (*(MomZ.pX))*(*(MomZ.pXP));
		*(MomZ.pXPXP) = (float)(Lambda_m_d_4Pi_E2/(MomZ.Mxx) + (*(MomZ.pXP))*(*(MomZ.pXP)));
	}
	//else
	//{
	//	*(MomZ.pXXP) = (float)((*(MomZ.pX))*(*(MomZ.pXP)) + sqrt(MaxMxxpMomZe2));
	//}
	double MaxMzzpMomZe2 = MomZ.Mzz*MomZ.Mzpzp - Lambda_m_d_4Pi_E2;
	if(MaxMzzpMomZe2 < 0.)
	{
		*(MomZ.pZZP) = (*(MomZ.pZ))*(*(MomZ.pZP));
		*(MomZ.pZPZP) = (float)(Lambda_m_d_4Pi_E2/(MomZ.Mzz) + (*(MomZ.pZP))*(*(MomZ.pZP)));
	}
	//else
	//{
	//	*(MomZ.pZZP) = (float)((*(MomZ.pZ))*(*(MomZ.pZP)) + sqrt(MaxMzzpMomZe2));
	//}
}

//*************************************************************************

void srTGenOptElem::TransferResizeParam(srTPropagScenario1D* PropagScenario, srTRadResize& ResBefore, srTRadResize& ResAfter)
{
	ResBefore.pxm = PropagScenario->ResizeBefore.pm;
	ResBefore.pxd = PropagScenario->ResizeBefore.pd;
	ResBefore.RelCenPosX = PropagScenario->ResizeBefore.RelCenPos;
	ResBefore.RelCenPosTol = PropagScenario->ResizeBefore.RelCenPosTol;

	ResAfter.pxm = PropagScenario->ResizeAfter.pm;
	ResAfter.pxd = PropagScenario->ResizeAfter.pd;
	ResAfter.RelCenPosX = PropagScenario->ResizeAfter.RelCenPos;
	ResAfter.RelCenPosTol = PropagScenario->ResizeAfter.RelCenPosTol;

	ResBefore.pzm = PropagScenario[1].ResizeBefore.pm;
	ResBefore.pzd = PropagScenario[1].ResizeBefore.pd;
	ResBefore.RelCenPosZ = PropagScenario[1].ResizeBefore.RelCenPos;
	if(ResBefore.RelCenPosTol > PropagScenario[1].ResizeBefore.RelCenPosTol) 
		ResBefore.RelCenPosTol = PropagScenario[1].ResizeBefore.RelCenPosTol;

	ResAfter.pzm = PropagScenario[1].ResizeAfter.pm;
	ResAfter.pzd = PropagScenario[1].ResizeAfter.pd;
	ResAfter.RelCenPosZ = PropagScenario[1].ResizeAfter.RelCenPos;
	if(ResAfter.RelCenPosTol > PropagScenario[1].ResizeAfter.RelCenPosTol) 
		ResAfter.RelCenPosTol = PropagScenario[1].ResizeAfter.RelCenPosTol;

	RejectSmallResize(ResBefore);
	RejectSmallResize(ResAfter);
}

//*************************************************************************

int srTGenOptElem::DefinePropagScenario(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTPredictedPropagData1D* PredictedPropagData, srTPropagScenario1D* PropagScenario)
{
	int result;

	srTRadSect1D Sect1D[2];
	if(result = SetupCharacteristicSections1D(pRadAccessData, Sect1D)) return result;
	//check/correct potential problem with off-axis aperture
	//to avoid zero intensity cuts
	//OC: seems fixed 170304

	PredictedPropagData[0].CoordMomRadCanBeModified = PredictedPropagData[0].AngMomRatCanBeModified = 1;
	PredictedPropagData[1].CoordMomRadCanBeModified = PredictedPropagData[1].AngMomRatCanBeModified = 1;
	if(result = DefinePropagScenario1D(Sect1D[0], ParPrecWfrPropag, PredictedPropagData[0], PropagScenario[0])) return result;
	if(result = DefinePropagScenario1D(Sect1D[1], ParPrecWfrPropag, PredictedPropagData[1], PropagScenario[1])) return result;

	if(result = CheckIfScenarioCanBeExecutedOrSuggestReduced(pRadAccessData, Sect1D, ParPrecWfrPropag, PredictedPropagData, PropagScenario)) return result;

	return 0;
}

//*************************************************************************

int srTGenOptElem::CheckIfScenarioCanBeExecutedOrSuggestReduced(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTPredictedPropagData1D* PredictedPropagData, srTPropagScenario1D* PropagScenario)
{// Check total memory needed for the execution and suggest a Reduced scenario if not enough memory
	int result;
	//char UnderSamplingModeTakesPlace = (pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ());
	char UnderSamplingModeWasSet = 0;

	//OC test: to uncomment if necessary
	//if(AllowAutoSwitchToUndersamplingMode() && (!UnderSamplingModeTakesPlace) && UnderSamplingModeCanBeSuggested(pRadAccessData, PropagScenario))
	//{
	//	if(result = TryToSetUnderSamplingMode(pRadAccessData, Sect1D, PropagScenario, UnderSamplingModeWasSet)) return 0;
	//	if(UnderSamplingModeWasSet) return 0;
	//}

	double ExtraSizeBeforeProp, ExtraSizeAfterProp;
	if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
		//test: checking available memory next time
		if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;

	if(AllowAutoSwitchToUndersamplingMode() && (!UnderSamplingModeWasSet))
	//if(AllowAutoSwitchToUndersamplingMode && (!UnderSamplingModeTakesPlace))
	{
		//to update for new treatment of the undersampling !!!
		//print warning ??
		if(result = TryToSetUnderSamplingMode(pRadAccessData, Sect1D, PropagScenario, UnderSamplingModeWasSet)) return 0;
		if(UnderSamplingModeWasSet) return 0;
	}

	CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);

	const double MinPoPerFrToDoResize = 1.1; // To steer
	const double MinPoPerFrToDoGradualResize = 1.3; //1.5; // To steer
	const int MaxAmOfTrendsForGradualResize = 30; //12; //4; // To steer

	srTFringeInfo &xFringeInfo = PropagScenario->CurFringeInfo, &zFringeInfo = (PropagScenario + 1)->CurFringeInfo;

	int xAmOfTrendsForGradualResize = MaxAmOfTrendsForGradualResize;
	double xMinPoPerFr = (xFringeInfo.LeftPointsPerFringe < xFringeInfo.RightPointsPerFringe)? xFringeInfo.LeftPointsPerFringe : xFringeInfo.RightPointsPerFringe;
	if(xMinPoPerFr < MinPoPerFrToDoGradualResize) 
	{
		xAmOfTrendsForGradualResize = 1;
		if(xMinPoPerFr < MinPoPerFrToDoResize) xAmOfTrendsForGradualResize = 0;
	}
	double xResCoef = (xMinPoPerFr > MinPoPerFrToDoResize)? pow(1./xMinPoPerFr, 1./xAmOfTrendsForGradualResize) : 1.;

	int zAmOfTrendsForGradualResize = MaxAmOfTrendsForGradualResize;
	double zMinPoPerFr = (zFringeInfo.LeftPointsPerFringe < zFringeInfo.RightPointsPerFringe)? zFringeInfo.LeftPointsPerFringe : zFringeInfo.RightPointsPerFringe;
	if(zMinPoPerFr < MinPoPerFrToDoGradualResize) 
	{
		zAmOfTrendsForGradualResize = 1;
		if(zMinPoPerFr < MinPoPerFrToDoResize) zAmOfTrendsForGradualResize = 0;
	}
	double zResCoef = (zMinPoPerFr > MinPoPerFrToDoResize)? pow(1./zMinPoPerFr, 1./zAmOfTrendsForGradualResize) : 1.;

	srTRadResize1D ResizeParamX, ResizeParamZ;

	while((xAmOfTrendsForGradualResize > 0) || (zAmOfTrendsForGradualResize > 0))
	{
		if(result = srYield.Check()) return result;

		if(xAmOfTrendsForGradualResize > 0)
		{
			//srTRadResize1D ResizeParam; ResizeParam.pd = xResCoef;
			ResizeParamX.pd *= xResCoef;

			srTRadSect1D SectDupl;
			if(result = Sect1D->SetupDupl(SectDupl)) return result;
			//if(result = RadResizeGen1D(SectDupl, ResizeParam)) return result;
			if(result = RadResizeGen1D(SectDupl, ResizeParamX)) return result;

			//if(result = DefinePropagScenario1D(SectDupl, PredictedPropagData[0], PropagScenario[0])) return result;
			if(result = DefinePropagScenario1D(SectDupl, ParPrecWfrPropag, PredictedPropagData[0], PropagScenario[0])) return result;
			//PropagScenario[0].ResizeBefore.pd *= ResizeParam.pd;
			PropagScenario[0].ResizeBefore.pd *= ResizeParamX.pd;

			xAmOfTrendsForGradualResize--;
		}
		if(zAmOfTrendsForGradualResize > 0)
		{
			//srTRadResize1D ResizeParam; ResizeParam.pd = zResCoef;
			ResizeParamZ.pd *= zResCoef;

			srTRadSect1D SectDupl;
			if(result = (Sect1D + 1)->SetupDupl(SectDupl)) return result;
			//if(result = RadResizeGen1D(SectDupl, ResizeParam)) return result;
			if(result = RadResizeGen1D(SectDupl, ResizeParamZ)) return result;

			//if(result = DefinePropagScenario1D(SectDupl, PredictedPropagData[1], PropagScenario[1])) return result;
			if(result = DefinePropagScenario1D(SectDupl, ParPrecWfrPropag, PredictedPropagData[1], PropagScenario[1])) return result;
			//PropagScenario[1].ResizeBefore.pd *= ResizeParam.pd;
			PropagScenario[1].ResizeBefore.pd *= ResizeParamZ.pd;

			zAmOfTrendsForGradualResize--;
		}
		if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
			//test: checking available memory next time
			if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
	}
	return SuggestScenarioThatFitsMemory(pRadAccessData, PropagScenario);
}

//*************************************************************************

int srTGenOptElem::SuggestScenarioThatFitsMemory(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario)
{
	const int MinNp = 18;

	double ExtraSizeBeforeProp, ExtraSizeAfterProp;
	if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
		//test: checking available memory next time
		if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;

	if(ExtraSizeAfterProp > ExtraSizeBeforeProp)
	{
		PropagScenario->ResizeAfter.pm = 1.;
		PropagScenario->ResizeAfter.pd = 1.;
		(PropagScenario + 1)->ResizeAfter.pm = 1.;
		(PropagScenario + 1)->ResizeAfter.pd = 1.;
		if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
			//test: checking available memory next time
			if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
	}

	const double ReduceFactor = 0.8; // To steer
	const int MaxAmOfTrends = 40; // To steer
	const double WfrEdgeSafetyCoef = 0.9;

	double InvReduceFactor = 1./ReduceFactor;

	double &pmxBefore = PropagScenario->ResizeBefore.pm;
	double &pdxBefore = PropagScenario->ResizeBefore.pd;
	double &pmzBefore = (PropagScenario + 1)->ResizeBefore.pm;
	double &pdzBefore = (PropagScenario + 1)->ResizeBefore.pd;

	double xHalfMagn = (pRadAccessData->nx >> 1)*pRadAccessData->xStep;
	double xMid = pRadAccessData->xStart + xHalfMagn;
	double zHalfMagn = (pRadAccessData->nz >> 1)*pRadAccessData->zStep;
	double zMid = pRadAccessData->zStart + zHalfMagn;

	for(int i=0; i<MaxAmOfTrends; i++)
	{
		double pxTotBeforeCurr = pdxBefore*pmxBefore;
		if(pxTotBeforeCurr > 1.)
		{
			pdxBefore *= ReduceFactor;
			double pxTot = pdxBefore*pmxBefore;
			if(pRadAccessData->nx*pxTot < MinNp) pdxBefore *= InvReduceFactor;
			
			if(pRadAccessData->WfrEdgeCorrShouldBeDone)
			{
				double xNewHalfMagnExtr = xHalfMagn*pmxBefore*ReduceFactor*WfrEdgeSafetyCoef;
				if((pRadAccessData->xWfrMin > xMid - xNewHalfMagnExtr) && (pRadAccessData->xWfrMax < xMid + xNewHalfMagnExtr))
				{
					pmxBefore *= ReduceFactor;
					pxTot = pdxBefore*pmxBefore;
					if(pRadAccessData->nx*pxTot < MinNp) pmxBefore *= InvReduceFactor;
				}
			}
			else
			{
				pmxBefore *= ReduceFactor;
				pxTot = pdxBefore*pmxBefore;
				if(pRadAccessData->nx*pxTot < MinNp) pmxBefore *= InvReduceFactor;
			}
		}

		double pzTotBeforeCurr = pdzBefore*pmzBefore;
		if(pzTotBeforeCurr > 1.)
		{
			pdzBefore *= ReduceFactor;
			double pzTot = pdzBefore*pmzBefore;
			if(pRadAccessData->nz*pzTot < MinNp) pdzBefore *= InvReduceFactor;
			
			if(pRadAccessData->WfrEdgeCorrShouldBeDone)
			{
				double zNewHalfMagnExtr = zHalfMagn*pmzBefore*ReduceFactor*WfrEdgeSafetyCoef;
				if((pRadAccessData->zWfrMin > zMid - zNewHalfMagnExtr) && (pRadAccessData->zWfrMax < zMid + zNewHalfMagnExtr))
				{
					pmzBefore *= ReduceFactor;
					pzTot = pdzBefore*pmzBefore;
					if(pRadAccessData->nz*pzTot < MinNp) pmzBefore *= InvReduceFactor;
				}
			}
			else
			{
				pmzBefore *= ReduceFactor;
				pzTot = pdzBefore*pmzBefore;
				if(pRadAccessData->nz*pzTot < MinNp) pmzBefore *= InvReduceFactor;
			}
		}

		if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
			//test: checking available memory next time
			if(MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp)) return 0;
	}

	pmxBefore = pdxBefore = pmzBefore = pdzBefore = 1.;
	PropagScenario->ResizeAfter.pm = 1.;
	PropagScenario->ResizeAfter.pd = 1.;
	(PropagScenario + 1)->ResizeAfter.pm = 1.;
	(PropagScenario + 1)->ResizeAfter.pd = 1.;
	return 0;
}

//*************************************************************************

void srTGenOptElem::EstimateMemoryNeededForPropag(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario, double& TotExtraSizeBeforeProp, double& TotExtraSizeAfterProp)
{
	double pmxBefore = PropagScenario->ResizeBefore.pm;
	double pdxBefore = PropagScenario->ResizeBefore.pd;
	double pmzBefore = (PropagScenario + 1)->ResizeBefore.pm;
	double pdzBefore = (PropagScenario + 1)->ResizeBefore.pd;

	long nxCurRad = pRadAccessData->nx, nzCurRad = pRadAccessData->nz;
	TotExtraSizeBeforeProp = ExtraMemSizeForResize(nxCurRad, nzCurRad, pmxBefore, pdxBefore, pmzBefore, pdzBefore, 0);

	double pmxAfter = PropagScenario->ResizeAfter.pm;
	double pdxAfter = PropagScenario->ResizeAfter.pd;
	double pmzAfter = (PropagScenario + 1)->ResizeAfter.pm;
	double pdzAfter = (PropagScenario + 1)->ResizeAfter.pd;

	//nxCurRad *= long(pmxBefore*pdxBefore); //?
	//nzCurRad *= long(pmzBefore*pdzBefore); //?

	double d_nxCurRad = ((double)nxCurRad)*pmxBefore*pdxBefore;
	nxCurRad = (long)d_nxCurRad;
	double d_nzCurRad = ((double)nzCurRad)*pmzBefore*pdzBefore;
	nzCurRad = (long)d_nzCurRad;

	//TotExtraSizeAfterProp = ExtraMemSizeForResize(nxCurRad, nzCurRad, pmxAfter, pdxAfter, pmzAfter, pdzAfter, 0);
	TotExtraSizeAfterProp = ExtraMemSizeForResize(nxCurRad, nzCurRad, pmxAfter, pdxAfter, pmzAfter, pdzAfter, 1);
}

//*************************************************************************

double srTGenOptElem::ExtraMemSizeForResize(long nxCurRad, long nzCurRad, double pxm, double pxd, double pzm, double pzd, char Mode)
{// Mode = 0 - economic for memory; 1 - normal (better performance)
	double TotExtraSize = 0.;
	if((pxm == 1.) && (pxd == 1.) && (pzm == 1.) && (pzd == 1.)) return TotExtraSize;

	double ResizeParamProd = pxm*pxd*pzm*pzd;

	double SizeOnePoint = 2.*sizeof(float);
	double SizeAllCurPoints = SizeOnePoint*nxCurRad*nzCurRad;

	if(Mode == 0) // Economic mode
	{
		if(ResizeParamProd > 1.) 
		{
			double Buf = 2.*(ResizeParamProd - 1)*SizeAllCurPoints;
			TotExtraSize = SizeAllCurPoints + Buf;
		}
		else
		{
			TotExtraSize = ResizeParamProd*SizeAllCurPoints;
		}
	}
	else if(Mode == 1) // Normal mode
	{
		TotExtraSize = 2.*ResizeParamProd*SizeAllCurPoints;
	}
	return TotExtraSize;
}

//*************************************************************************

double srTGenOptElem::ExtraMemSizeForResizeE(long neCurRad, long nxCurRad, long nzCurRad, double pem, double ped, char Mode)
{// Mode = 0 - economic for memory; 1 - normal (better performance)
	double TotExtraSize = 0.;
	if((pem == 1.) && (ped == 1.)) return TotExtraSize;

	double ResizeParamProd = pem*ped;

	double SizeOnePoint = 2.*sizeof(float);
	double SizeAllCurPoints = SizeOnePoint*neCurRad*nxCurRad*nzCurRad;

	if(Mode == 0) // Economic mode
	{
		if(ResizeParamProd > 1.) 
		{
			double Buf = 2.*(ResizeParamProd - 1)*SizeAllCurPoints;
			TotExtraSize = SizeAllCurPoints + Buf;
		}
		else
		{
			TotExtraSize = ResizeParamProd*SizeAllCurPoints;
		}
	}
	else if(Mode == 1) // Normal mode
	{
		TotExtraSize = 2.*ResizeParamProd*SizeAllCurPoints;
	}
	return TotExtraSize;
}

//*************************************************************************

int srTGenOptElem::SetupCharacteristicSections1D(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* Sect1DArr)
{
	int result;

	srTRadSect1D &SectVsX = Sect1DArr[0], &SectVsZ = Sect1DArr[1];
	SectVsX.VsXorZ = 'x';
	SectVsZ.VsXorZ = 'z';

	long &ixc = SectVsZ.icOtherCoord;
	long &izc = SectVsX.icOtherCoord;

	double xWfrMinLoc = pRadAccessData->xStart;
	double xWfrMaxLoc = pRadAccessData->xStart + (pRadAccessData->xStep)*(pRadAccessData->nx);
	double zWfrMinLoc = pRadAccessData->zStart;
	double zWfrMaxLoc = pRadAccessData->zStart + (pRadAccessData->zStep)*(pRadAccessData->nz);
	
	if(pRadAccessData->WfrXLimitsSeemDefined())
	{
		xWfrMinLoc = pRadAccessData->xWfrMin;
		xWfrMaxLoc = pRadAccessData->xWfrMax;
	}
	if(pRadAccessData->WfrZLimitsSeemDefined())
	{
		zWfrMinLoc = pRadAccessData->zWfrMin;
		zWfrMaxLoc = pRadAccessData->zWfrMax;
	}

	if(WfrTransmLimits.LimitsAreDefined)
	{
		if(xWfrMinLoc < WfrTransmLimits.Xmin) xWfrMinLoc = WfrTransmLimits.Xmin;
		if(xWfrMaxLoc > WfrTransmLimits.Xmax) xWfrMaxLoc = WfrTransmLimits.Xmax;
		if(zWfrMinLoc < WfrTransmLimits.Ymin) zWfrMinLoc = WfrTransmLimits.Ymin;
		if(zWfrMaxLoc > WfrTransmLimits.Ymax) zWfrMaxLoc = WfrTransmLimits.Ymax;
	}

	//long ixWFmin = (pRadAccessData->xWfrMin - pRadAccessData->xStart)/pRadAccessData->xStep;
	long ixWFmin = (long)((xWfrMinLoc - pRadAccessData->xStart)/pRadAccessData->xStep);
	if((ixWFmin < 0) || (ixWFmin >= pRadAccessData->nx)) ixWFmin = 0;
	//long ixWFmax = (pRadAccessData->xWfrMax - pRadAccessData->xStart)/pRadAccessData->xStep;
	long ixWFmax = (long)((xWfrMaxLoc - pRadAccessData->xStart)/pRadAccessData->xStep);
	if((ixWFmax < 0) || (ixWFmax >= pRadAccessData->nx)) ixWFmax = pRadAccessData->nx - 1;
	long ixRel = ((ixWFmax - ixWFmin) >> 1) - 1; if(ixRel < 0) ixRel = 0;
	ixc = ixWFmin + ixRel;

	//long izWFmin = (pRadAccessData->zWfrMin - pRadAccessData->zStart)/pRadAccessData->zStep;
	long izWFmin = (long)((zWfrMinLoc - pRadAccessData->zStart)/pRadAccessData->zStep);
	if((izWFmin < 0) || (izWFmin >= pRadAccessData->nz)) izWFmin = 0;
	//long izWFmax = (pRadAccessData->zWfrMax - pRadAccessData->zStart)/pRadAccessData->zStep;
	long izWFmax = (long)((zWfrMaxLoc - pRadAccessData->zStart)/pRadAccessData->zStep);
	if((izWFmax < 0) || (izWFmax >= pRadAccessData->nz)) izWFmax = pRadAccessData->nz - 1;
	long izRel = ((izWFmax - izWFmin) >> 1) - 1; if(izRel < 0) izRel = 0;
	izc = izWFmin + izRel;

	if(ixc < 0) ixc = 0;
	if(ixc >= pRadAccessData->nx) ixc = pRadAccessData->nx - 1;
	if(izc < 0) izc = 0;
	if(izc >= pRadAccessData->nz) izc = pRadAccessData->nz - 1;

	if(result = SetupSectionArraysVsXandZ(pRadAccessData, SectVsX, SectVsZ)) return result;

	const double RatForMinE = 0.05;
	srTMinMaxEParam MinMaxE2D;
	pRadAccessData->FindMinMaxReE(MinMaxE2D);
	float MaxAbsE2D;
	//long xIndMaxAbsE2D, zIndMaxAbsE2D;
	long long xIndMaxAbsE2D, zIndMaxAbsE2D;
	MinMaxE2D.FindGenMaxAbsE(MaxAbsE2D, xIndMaxAbsE2D, zIndMaxAbsE2D);

	float MaxAbsEVsX, MaxAbsEVsZ;
	//long IndMaxAbsEVsX, IndMaxAbsEVsZ, IndDummy;
	long long IndMaxAbsEVsX, IndMaxAbsEVsZ, IndDummy;
	srTMinMaxEParam MinMaxE1D;
	SectVsX.FindMinMaxReE(MinMaxE1D);
	MinMaxE1D.FindGenMaxAbsE(MaxAbsEVsX, IndMaxAbsEVsX, IndDummy);
	SectVsZ.FindMinMaxReE(MinMaxE1D);
	MinMaxE1D.FindGenMaxAbsE(MaxAbsEVsZ, IndMaxAbsEVsZ, IndDummy);

	//long ixcPrev = ixc, izcPrev = izc;
	long long ixcPrev = ixc, izcPrev = izc;
	float MaxAbsE2D_RatForMinE = (float)(MaxAbsE2D*RatForMinE);
	if(MaxAbsEVsX < MaxAbsE2D_RatForMinE) 
	{
		izc = (long)zIndMaxAbsE2D;
	}
	if(MaxAbsEVsZ < MaxAbsE2D_RatForMinE)
	{
		ixc = (long)xIndMaxAbsE2D;
	}
	
	//if((izc != izcPrev) || (ixc != ixcPrev))
	//OC1703
	if(((izc != izcPrev) && (izc >= izWFmin) && (izc <= izWFmax)) || ((ixc != ixcPrev) && (ixc >= ixWFmin) && (ixc <= ixWFmax)))
	{
		SectVsX.DeleteArrays(); SectVsX.DeleteArrays();
		if(result = SetupSectionArraysVsXandZ(pRadAccessData, SectVsX, SectVsZ)) return result;
	}

	if(pRadAccessData->ThereIsUnderSamplingX())
	{
		srTRadResize1D Resize1D;
		Resize1D.pm = 1.; Resize1D.pd = pRadAccessData->UnderSamplingX;
		if(result = RadResizeGen1D(SectVsX, Resize1D)) return result;
		SectVsX.ThereIsUnderSamplingIn2D = 1;
		SectVsZ.ThereIsUnderSamplingIn2D = 1;
	}
	if(pRadAccessData->ThereIsUnderSamplingZ())
	{
		srTRadResize1D Resize1D;
		Resize1D.pm = 1.; Resize1D.pd = pRadAccessData->UnderSamplingZ;
		if(result = RadResizeGen1D(SectVsZ, Resize1D)) return result;
		SectVsX.ThereIsUnderSamplingIn2D = 1;
		SectVsZ.ThereIsUnderSamplingIn2D = 1;
	}

	return 0;
}

//*************************************************************************

int srTGenOptElem::DefinePropagScenario1D(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTPredictedPropagData1D& PredictedPropagData, srTPropagScenario1D& PropagScenario)
{
	const double RelResizeTol = 0.02; //0.1; //Small tolerance is necessary for wide wavefronts
	//const int AmOfIterToAdjustRangeAndResol = 10; // To steer
	int result;

	char TreatExOrEz = ChooseTreatExOrEzBasedOnMax(Sect1D);

	srTRadResize1D ResizeParam;
	if((PredictedPropagData.CoordMomRat > 0.) && (::fabs(PredictedPropagData.CoordMomRat - 1.) > RelResizeTol)) 
	{
		ResizeParam.pm = PredictedPropagData.CoordMomRat;

		//OC
		//double CoordSize = 4.*PredictedPropagData.SigmaCoord; // To steer
		//double CoordRange = Sect1D.ArgStep*Sect1D.np;
		//if(CoordRange > CoordSize)
		//{
		//	if(ResizeParam.pm > 1.) ResizeParam.pm *= CoordSize/CoordRange;
		//}
	}
	if((PredictedPropagData.AngMomRat > 0.) && (::fabs(PredictedPropagData.AngMomRat - 1.) > RelResizeTol)) 
	{
		ResizeParam.pd = PredictedPropagData.AngMomRat;

		//OC
		//double AngSize = 4.*PredictedPropagData.SigmaAng; // To steer
		//double InvLambda_m = Sect1D.eVal*806546.577258;
		//double AngRange = 1./(Sect1D.ArgStep*InvLambda_m);
		//if(AngRange > AngSize)
		//{
		//	if(ResizeParam.pd > 1.) ResizeParam.pd *= AngSize/AngRange;
		//}
	}

	const double MinInitResizePm = 0.5;
	const double MinInitResizePd = 0.5;
	if(ResizeParam.pm < MinInitResizePm) ResizeParam.pm = MinInitResizePm;
	if(ResizeParam.pd < MinInitResizePd) ResizeParam.pd = MinInitResizePd;

	srTRadResize1D InitialResizeParam = ResizeParam;

	char RangeTuningIsNeeded = (RangeShouldBeAdjustedAtPropag() && PredictedPropagData.CoordMomRadCanBeModified);
	char ResolutionTuningIsNeeded = (ResolutionShouldBeAdjustedAtPropag() && PredictedPropagData.AngMomRatCanBeModified);

// Find Minimal pd allowing not to loose the Resolution at propagation
	srTRadResize1D ResizeParamToKeepResol = ResizeParam;
	if(ResolutionTuningIsNeeded)
	{
		//if(result = TuneAndKeepResolution1D(Sect1D, ResizeParamToKeepResol, PropagScenario.CurFringeInfo)) return result;
		if(result = TuneAndKeepResolution1D(Sect1D, ParPrecWfrPropag, ResizeParamToKeepResol, PropagScenario.CurFringeInfo)) return result;
		
		ResizeParam.pd = ResizeParamToKeepResol.pd;
	}
	else 
	{
		ResizeParam.pd = 1.;
	}

// Find Minimal pm that does not degrade the Precision at propagation
	srTRadResize1D PostResizeForRange;
	if(RangeTuningIsNeeded)
	{
		if(!Sect1D.ThereIsUnderSamplingIn2D)
		{
			srTRadResize1D ResizeParamMinRangeForPrec;
			if(ResizeParam.pm > 1. + RelResizeTol) ResizeParamMinRangeForPrec = ResizeParam;

			//if(result = TuneRangeNotDegradingPrec1D(Sect1D, ResizeParamMinRangeForPrec)) return result;
			if(result = TuneRangeNotDegradingPrec1D(Sect1D, ParPrecWfrPropag, ResizeParamMinRangeForPrec)) return result;
			
			if(ResizeParamMinRangeForPrec.pm > 1. + RelResizeTol) ResizeParam.pm = ResizeParamMinRangeForPrec.pm;
		}
		else
		{
			ResizeParam.pm = 1.;
		}

		////test
		//char ErrorMesTitle[] = "SRW Debug";
        //char ErrorStr[100];
		//int j = sprintf(ErrorStr, "pm: %e", ResizeParam.pm);
		////j += sprintf(ErrorStr + j, "  pxm: %e", ResBefore.pxm);
		////j += sprintf(ErrorStr + j, "  pzd: %e", ResBefore.pzd);
		////j += sprintf(ErrorStr + j, "  pzm: %e", ResBefore.pzm);

		//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
        //int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
		////end test

		PostResizeForRange = ResizeParam;

		//if(result = FindPostResizeForRange1D(Sect1D, PostResizeForRange)) return result;
		if(result = FindPostResizeForRange1D(Sect1D, ParPrecWfrPropag, PostResizeForRange)) return result;

		if(ResizeParam.pm < 1.)
		{
			ResizeParam.pm *= PostResizeForRange.pm;
		}

		if(result = AnalizeFringes(Sect1D, TreatExOrEz, PropagScenario.CurFringeInfo)) return result;
	}

	if(ResizeParam.pm < 1.)
	{
		PropagScenario.ResizeAfter.pm = ResizeParam.pm;
	}
	else if(ResizeParam.pm >= 1.)
	{
		PropagScenario.ResizeBefore.pm = ResizeParam.pm;
		if(PostResizeForRange.pm < 1.)
		{
			PropagScenario.ResizeAfter.pm = PostResizeForRange.pm;
		}
	}

	if(ResizeParam.pd < 1.)
	{
		PropagScenario.ResizeAfter.pd = ResizeParam.pd;
	}
	else if(ResizeParam.pd >= 1.)
	{
		PropagScenario.ResizeBefore.pd = ResizeParam.pd;
	}

	if(RangeTuningIsNeeded)
	{
		//if(result = FindPostResizeCenterCorrection(Sect1D, PropagScenario)) return result;
		if(result = FindPostResizeCenterCorrection(Sect1D, ParPrecWfrPropag, PropagScenario)) return result;
	}
	CorrectResParMinNumPo(Sect1D.np, PropagScenario.ResizeBefore, PropagScenario.ResizeAfter);

	return 0;
}

//*************************************************************************

//int srTGenOptElem::FindPostResizeCenterCorrection(srTRadSect1D& Sect1D, srTPropagScenario1D& PropagScenario)
int srTGenOptElem::FindPostResizeCenterCorrection(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTPropagScenario1D& PropagScenario)
{
	double RelRangeToleranceForIntensity = 0.001; // Steer well or set from external variable !
//OC: to scale with precision
	RelRangeToleranceForIntensity /= ParPrecWfrPropag.PrecFact; //??
	
	//const double RelTolForE = sqrt(RelRangeToleranceForIntensity); //OC030110

	const double MinPostResize = 0.1; // To steer
	const double NotResizeTol = 0.05;

	srTRadResize1D& ResizeParam = PropagScenario.ResizeBefore;
	srTRadResize1D& PostResizeParam = PropagScenario.ResizeAfter;

	int result;
	srTRadSect1D SectDpl;
	if(result = Sect1D.SetupDupl(SectDpl)) return result;
	if((ResizeParam.pm < 1.) || (ResizeParam.pd < 1.))
	{
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result;
		if(result = RadResizeGen1D(SectDpl, ResizeParam)) return result;
	}
	else
	{
		if(result = RadResizeGen1D(SectDpl, ResizeParam)) return result;
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result;
	}

	//long MaxAbsReExInd;
	long long MaxAbsReExInd;
	float MaxAbsReEx, MaxAbsEx, MaxAbsEz;
	char TreatExOrEz;
	//long IndMaxAbsEx, IndMaxAbsEz;
	long long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(SectDpl, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	if(MaxAbsEx > MaxAbsEz)
	{
		TreatExOrEz = 'x'; MaxAbsReEx = MaxAbsEx; MaxAbsReExInd = IndMaxAbsEx;
	}
	else
	{
		TreatExOrEz = 'z'; MaxAbsReEx = MaxAbsEz; MaxAbsReExInd = IndMaxAbsEz;
	}

	long iLeftThreshBorder, iRightThreshBorder;
	FindIntensityBorders1D(SectDpl, TreatExOrEz, RelRangeToleranceForIntensity, iLeftThreshBorder, iRightThreshBorder);

	bool ModifyPmEvenIfCenPosIsNotSet = false;
	CheckRelCenPosAndSetPostResizeParamPmIfNecessary(SectDpl.np, iLeftThreshBorder, iRightThreshBorder, PostResizeParam, ModifyPmEvenIfCenPosIsNotSet);
	if(PostResizeParam.pm < MinPostResize) PostResizeParam.pm = MinPostResize;
	if(::fabs(PostResizeParam.pm - 1) < NotResizeTol) PostResizeParam.pm = 1;

/**
	long ic = iLeftThreshBorder + ((iRightThreshBorder - iLeftThreshBorder) >> 1);
	double xc = SectDpl.ArgStart + ic*SectDpl.ArgStep;
	
	double ArgRange = SectDpl.np*SectDpl.ArgStep;
	double RelCenPos = (xc - SectDpl.ArgStart)/ArgRange;
	if(::fabs(RelCenPos - 0.5) > PostResizeParam.RelCenPosTol)
	{
		PostResizeParam.RelCenPos = RelCenPos;
		double NewRange = (iRightThreshBorder - iLeftThreshBorder)*SectDpl.ArgStep;
		PostResizeParam.pm = (NewRange/(SectDpl.np*SectDpl.ArgStep))*1.3; // To steer

		if(PostResizeParam.pm < MinPostResize) PostResizeParam.pm = MinPostResize;
	}
**/
	return 0;
}

//*************************************************************************

void srTGenOptElem::CheckRelCenPosAndSetPostResizeParamPmIfNecessary(long np, long iLeftThreshBorder, long iRightThreshBorder, srTRadResize1D& PostResizeParam, bool ModifyPmEvenIfCenPosIsNotSet)
{
	const long MinNp = 50;

	if(np <= MinNp) return;

	iLeftThreshBorder -= 1;
	iRightThreshBorder += 1;

	if(iLeftThreshBorder < 0) iLeftThreshBorder = 0;
	if(iRightThreshBorder >= np) iRightThreshBorder = np - 1;
	if(iLeftThreshBorder >= iRightThreshBorder) return;
	long OrigNp = np;
	long iOrigRightLeftIntDif = iRightThreshBorder - iLeftThreshBorder;
	long iOrigIntCen = iLeftThreshBorder + (iOrigRightLeftIntDif >> 1);

/**
	if(PostResizeParam.pm != 1)
	{
		long NewNp = long(np*PostResizeParam.pm);
		long ExtraNp = (NewNp - np) >> 1;
		long New_iLeftThreshBorder = iLeftThreshBorder + ExtraNp;
		long New_iRightThreshBorder = iRightThreshBorder + ExtraNp;

		if((New_iLeftThreshBorder < 0) || (New_iRightThreshBorder >= NewNp))
		{
			long iDeltaLeft = 0;
			if(New_iLeftThreshBorder < 0) iDeltaLeft = -New_iLeftThreshBorder;

			long iDeltaRight = 0;
			if(New_iRightThreshBorder >= NewNp) iDeltaRight = New_iRightThreshBorder - NewNp + 1;

			long iDeltaMax = (iDeltaLeft > iDeltaRight)? iDeltaLeft : iDeltaRight;
			NewNp += (iDeltaMax << 1);

			PostResizeParam.pm = double(NewNp)/double(np);
			ExtraNp = (NewNp - np) >> 1;
			New_iLeftThreshBorder = iLeftThreshBorder + ExtraNp;
			New_iRightThreshBorder = iRightThreshBorder + ExtraNp;

			if(New_iLeftThreshBorder < 0) New_iLeftThreshBorder = 0;
			if(New_iRightThreshBorder >= (NewNp - 1)) New_iRightThreshBorder = (NewNp - 1);
		}

		iLeftThreshBorder = New_iLeftThreshBorder;
		iRightThreshBorder = New_iRightThreshBorder;
		np = NewNp;
	}
**/

	//checking what would be the position of center if wfr is resized without center correction
	long iCurCen = (np >> 1);
	long iDeltaLeft = iCurCen - iLeftThreshBorder;
	long iDeltaRight = iRightThreshBorder - iCurCen;
	long iDeltaMax = (iDeltaLeft > iDeltaRight)? iDeltaLeft : iDeltaRight;
	//if(iDeltaMax <= 0) return;

	long TestNewNp = (iDeltaMax << 1);
	//long iTestNewCen = iDeltaMax; //(TestNewNp >> 1);
	long iSubst = iCurCen - iDeltaMax;
	if(iSubst < 0) iSubst = 0;

	long iTestLeftThreshBorder = iLeftThreshBorder - iSubst;
	long iTestRightThreshBorder = iRightThreshBorder - iSubst;
	if(iTestLeftThreshBorder < 0) iTestLeftThreshBorder = 0;
	if(iTestRightThreshBorder >= (TestNewNp - 1)) iTestRightThreshBorder = (TestNewNp - 1);

	long iTestRightLeftDif = iTestRightThreshBorder - iTestLeftThreshBorder;
	if(iTestRightLeftDif <= 0) iTestRightLeftDif = MinNp;

	long iTestCen = iTestLeftThreshBorder + (iTestRightLeftDif >> 1);
	double TestRelCenPos = double(iTestCen)/double(TestNewNp);

	//long iIntCen = iLeftThreshBorder + (iRightLeftDif >> 1);
	////double TestRelCenPos = double(iIntCen)/double(NewNp);
	//double TestRelCenPos = double(iIntCen)/double(OrigNp);

	if(::fabs(TestRelCenPos - 0.5) > PostResizeParam.RelCenPosTol)
	{//set resize parameter Pm as for modified center position

		//PostResizeParam.RelCenPos = TestRelCenPos;
		PostResizeParam.RelCenPos = double(iOrigIntCen)/double(OrigNp);
		PostResizeParam.SetRelCenPosTolToEnforceCenterCorr(); 
		PostResizeParam.pm = (double(iOrigRightLeftIntDif)/double(OrigNp))*1.3; // To steer
	}
	else
	{//calculate Pm for the case without modifying center position
		if(ModifyPmEvenIfCenPosIsNotSet)
		{
			PostResizeParam.pm = double(TestNewNp)/double(OrigNp);
			PostResizeParam.RelCenPos = 0.5;
		}

	//	double ExtraPm = double(NewNp)/double(np);
	//	PostResizeParam.pm *= ExtraPm;
	}
}

//*************************************************************************

void srTGenOptElem::CorrectResParMinNumPo(long Np, srTRadResize1D& ResizeBefore, srTRadResize1D& ResizeAfter)
{
	const long MinNpAllowed = 18; // To steer
	long NewNp = long(ResizeBefore.pm*ResizeBefore.pd*Np);
	if(NewNp < MinNpAllowed)
	{
		double Coef = double(MinNpAllowed)/double(NewNp);
		if(ResizeBefore.pm < 1.) ResizeBefore.pm *= Coef;
		else if(ResizeBefore.pd < 1.) ResizeBefore.pd *= Coef;
		NewNp = MinNpAllowed;
	}

	NewNp = long(ResizeAfter.pm*ResizeAfter.pd*NewNp);
	if(NewNp < MinNpAllowed)
	{
		double Coef = double(MinNpAllowed)/double(NewNp);
		if(ResizeAfter.pm < 1.) ResizeAfter.pm *= Coef;
		else if(ResizeAfter.pd < 1.) ResizeAfter.pd *= Coef;
	}
}

//*************************************************************************

//int srTGenOptElem::TuneAndKeepResolution1D(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParamToKeepResol, srTFringeInfo& FringeInfo)
int srTGenOptElem::TuneAndKeepResolution1D(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResize1D& ResizeParamToKeepResol, srTFringeInfo& FringeInfo)
{
	double RelPointsPerFringeMisfit = 0.03; // To steer

	//OC: to scale with precision
	RelPointsPerFringeMisfit /= ParPrecWfrPropag.PrecFact; //???

	int MaxIterNum = 80;//40; // To steer
	double MaxPointsPerFringe = 10.; //8.; //5.; // To steer !!!
	MaxPointsPerFringe *= ParPrecWfrPropag.PrecFact;

	int result;
	char TreatExOrEz = ChooseTreatExOrEzBasedOnMax(Sect1D);

//	srTDoubleVect FringeCoor;
	if(result = AnalizeFringes(Sect1D, TreatExOrEz, FringeInfo)) return result;
	if(FringeInfo.LeftPointsPerFringe == 1.) FringeInfo.LeftPointsPerFringe = 1. + 1.1*RelPointsPerFringeMisfit;
	if(FringeInfo.RightPointsPerFringe == 1.) FringeInfo.RightPointsPerFringe = 1. + 1.1*RelPointsPerFringeMisfit;

	if(FringeInfo.LeftPointsPerFringe > MaxPointsPerFringe) FringeInfo.LeftPointsPerFringe = MaxPointsPerFringe;
	if(FringeInfo.RightPointsPerFringe > MaxPointsPerFringe) FringeInfo.RightPointsPerFringe = MaxPointsPerFringe;

	double MinPoPerFrBefore = (FringeInfo.LeftPointsPerFringe < FringeInfo.RightPointsPerFringe)? FringeInfo.LeftPointsPerFringe : FringeInfo.RightPointsPerFringe;

	char PrevResolNeedsToChange;
	srTRadResize1D ResizeParam = ResizeParamToKeepResol;

	for(int i=0; i<MaxIterNum; i++)
	{
		if(result = srYield.Check()) return result;

		srTFringeInfo PropFringeInfo;
		if(result = ResizePropagateAndAnalizeFringes1D(Sect1D, ResizeParam, TreatExOrEz, PropFringeInfo)) return result;
		if(result = CheckIfOversamplingIsReal(Sect1D, ResizeParam, TreatExOrEz, PropFringeInfo)) return result;

		//OC210508
		//if(PropFringeInfo.LeftPointsPerFringe > MaxPointsPerFringe) PropFringeInfo.LeftPointsPerFringe = MaxPointsPerFringe;
		//if(PropFringeInfo.RightPointsPerFringe > MaxPointsPerFringe) PropFringeInfo.RightPointsPerFringe = MaxPointsPerFringe;

		double MinPoPerFrAfter = (PropFringeInfo.LeftPointsPerFringe < PropFringeInfo.RightPointsPerFringe)? PropFringeInfo.LeftPointsPerFringe : PropFringeInfo.RightPointsPerFringe;
		char FringesOK = (::fabs(MinPoPerFrBefore - MinPoPerFrAfter) < MinPoPerFrBefore*RelPointsPerFringeMisfit);

		double pdLoc = 1.;
		if(!FringesOK)
		{
			pdLoc = SuggestResolResizeCoef(MinPoPerFrBefore, MinPoPerFrAfter, PropFringeInfo.AmOfFringes);
			if(pdLoc == 1.) FringesOK = 1;
		}
		char ResolNeedsIncr = (pdLoc > 1.);
		char ResolNeedsDecr = (pdLoc < 1.);

		if(i > 0)
		{
			if(ResolNeedsIncr && (PrevResolNeedsToChange == -1)) break;
			if(ResolNeedsDecr && (PrevResolNeedsToChange == 1)) break;
		}
		PrevResolNeedsToChange = (ResolNeedsDecr? -1 : (ResolNeedsIncr? 1 : 0));

		if(!(ResolNeedsIncr || ResolNeedsDecr)) break;
		else
		{
			ResizeParam.pd *= pdLoc;
		}
	}

	ResizeParamToKeepResol = ResizeParam;

	return 0;
}

//*************************************************************************

int srTGenOptElem::AnalizeFringes(srTRadSect1D& Sect1D, char TreatExOrEz, srTFringeInfo& FringeInfo)
{
	int result;

	srTIntVect FringeContent;
	srTDoubleVect FringeCoor;
	if(result = CountFringes(Sect1D, FringeContent, TreatExOrEz, FringeCoor)) return result;

	const int AmOfFringesToSwitchToOneInterval = 9; // To steer
	long TotAmOfFringes = (long)FringeContent.size();

	int AmOfFringesToTreatAtEdges = (TotAmOfFringes < 25)? 6 : 10; //4 : 8; // To steer
	if(TotAmOfFringes > 50) AmOfFringesToTreatAtEdges = 15; //12;
	if(TotAmOfFringes > 500) AmOfFringesToTreatAtEdges = 80; //12;

	char TreatOnlyOneInterval = (TotAmOfFringes < AmOfFringesToSwitchToOneInterval);
	char TreatTwoEdgeIntervals = !TreatOnlyOneInterval;

	double LeftPointsPerFringe, RightPointsPerFringe;

	if(TreatOnlyOneInterval)
	{
		int iSt, iFi;
		if(TotAmOfFringes > 3)
		{
			iSt = 1; iFi = TotAmOfFringes - 2;
		}
		else
		{
			iSt = 0; iFi = TotAmOfFringes - 1; 
		}
		int TotFringes = iFi;

		//OC210508
		if(TotFringes <= 0) TotFringes = 1;

		//long SumOfPoints = 0;
		long long SumOfPoints = 0;
		for(int i=iSt; i<=iFi; i++) SumOfPoints += FringeContent[i];
		LeftPointsPerFringe = 0.8*double(SumOfPoints)/double(TotFringes); // To steer
		RightPointsPerFringe = LeftPointsPerFringe;
	}
	else if(TreatTwoEdgeIntervals)
	{
		//long SumOfPoints = 0;
		long long SumOfPoints = 0;
		for(int i1=1; i1<(AmOfFringesToTreatAtEdges + 1); i1++) SumOfPoints += FringeContent[i1];
		LeftPointsPerFringe = double(SumOfPoints)/double(AmOfFringesToTreatAtEdges);

		SumOfPoints = 0;
		for(int i2=(TotAmOfFringes - AmOfFringesToTreatAtEdges - 1); i2<(TotAmOfFringes - 1); i2++) SumOfPoints += FringeContent[i2];
		RightPointsPerFringe = double(SumOfPoints)/double(AmOfFringesToTreatAtEdges);
	}

	FringeInfo.AmOfFringes = TotAmOfFringes;
	FringeInfo.LeftPointsPerFringe = LeftPointsPerFringe;
	FringeInfo.RightPointsPerFringe = RightPointsPerFringe;
	return 0;
}

//*************************************************************************

int srTGenOptElem::AnalizeFringes2D(srTSRWRadStructAccessData* pRadAccessData, srTFringeInfo* FringeInfoArr)
{
	int result;
	srTRadSect1D Sect1DArr[2];
	if(result = SetupCharacteristicSections1D(pRadAccessData, Sect1DArr)) return result;

	char pol0 = Sect1DArr->DetermineMoreIntensivePolarComp();
	char pol1 = (Sect1DArr + 1)->DetermineMoreIntensivePolarComp();
	if(result = AnalizeFringes(Sect1DArr[0], pol0, FringeInfoArr[0])) return result;
	if(result = AnalizeFringes(Sect1DArr[1], pol1, FringeInfoArr[1])) return result;
	return 0;
}

//*************************************************************************

int srTGenOptElem::CountFringes(srTRadSect1D& Sect1D, srTIntVect& FringeContent, char TreatExOrEz, srTDoubleVect& FringeCoor)
{
/**
	const double RelZeroTolForIntensBorders = 0.005; //0.002; // To steer
	const double RelZeroTolForFringeCount = 0.01; //0.5*RelZeroTolForIntensBorders; // To steer
	long iFirst, iLast;
	FindIntensityBorders1D(Sect1D, TreatExOrEz, RelZeroTolForIntensBorders, iFirst, iLast);

	float *tE = ((TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz) + ((iFirst + 1) << 1);

	char DerivSign = (*tE >= *(tE-2))? 1 : -1;
	int PointsInFringe = 1;
	long AmOfFringes = 0;

	char E_IsZero = (*tE == 0.);
	char LastFringeCounted = 0;

	for(long i=iFirst + 2; i<=iLast; i++)
	{
		tE += 2;

		double CurAbsE = ::fabs(*tE);
		double PrevAbsE = ::fabs(*(tE - 2));
		double MinAbsE = CurAbsE;
		if(PrevAbsE < MinAbsE) MinAbsE = PrevAbsE;
		double AbsTolE = RelZeroTolForFringeCount*MinAbsE;

		char NewE_IsZero = (*tE == 0.);
		if(!NewE_IsZero)
		{
			LastFringeCounted = 0;
			if(E_IsZero)
			{
				DerivSign = (*tE >= *(tE-2))? 1 : -1;
			}
			else
			{
				//char NewDerivSign = (*tE >= *(tE-2))? 1 : -1;
				char NewDerivSign = (*tE >= (*(tE-2) + AbsTolE))? 1 : ((*tE < (*(tE-2) - AbsTolE))? -1 : 0);
				if(NewDerivSign == 0) NewDerivSign = DerivSign;

				if(NewDerivSign != DerivSign)
				{
					FringeContent.push_back(PointsInFringe); AmOfFringes++;

					double FrCoor = Sect1D.ArgStart + i*Sect1D.ArgStep;
					FringeCoor.push_back(FrCoor);

					PointsInFringe = 1;
				}
				else
				{
					PointsInFringe++;
				}
				DerivSign = NewDerivSign;
			}
		}
		else
		{
			if(!E_IsZero)
			{
				char NewDerivSign = (*tE >= *(tE-2))? 1 : -1;

				if(NewDerivSign != DerivSign)
				{
					FringeContent.push_back(PointsInFringe); AmOfFringes++;

					double FrCoor = Sect1D.ArgStart + i*Sect1D.ArgStep;
					FringeCoor.push_back(FrCoor);

					PointsInFringe = 1;
					LastFringeCounted = 1;
				}
				else
				{
					PointsInFringe++;
				}
				DerivSign = NewDerivSign;
			}
		}
		E_IsZero = NewE_IsZero;
	}
	if(!LastFringeCounted)
	{
		double FrCoor = Sect1D.ArgStart + iLast*Sect1D.ArgStep;

		int AmOfFringeCoor = FringeCoor.size();
		double LastFrCoor = -1E+23;
		if(AmOfFringeCoor > 0)
		{
			LastFrCoor = FringeCoor[AmOfFringeCoor - 1];
		}
		if(FrCoor != LastFrCoor) 
		{
			FringeContent.push_back(PointsInFringe); AmOfFringes++;

			FringeCoor.push_back(FrCoor);
		}
	}
**/

//this counts number of times of crossing 0

	const double RelZeroTolForIntensBorders = 0.005; //0.002; // To steer
	long iFirst, iLast;
	FindIntensityBorders1D(Sect1D, TreatExOrEz, RelZeroTolForIntensBorders, iFirst, iLast);

	float *tE = ((TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz) + ((iFirst + 1) << 1);

	char CurSign = (*tE >= 0.)? 1 : -1;
	int PointsInFringe = 1;
	long AmOfFringes = 0;

	for(long i=iFirst + 2; i<=iLast; i++)
	{
		tE += 2;

		char NewSign = (*tE >= 0.)? 1 : -1;
		if(NewSign != CurSign)
		{
            FringeContent.push_back(PointsInFringe); AmOfFringes++;

			double FrCoor = Sect1D.ArgStart + i*Sect1D.ArgStep;
			FringeCoor.push_back(FrCoor);

			PointsInFringe = 1;
		}
		else
		{
			PointsInFringe++;
		}
		CurSign = NewSign;
	}
	int AmOfFringeCoor = (int)FringeCoor.size();
	if(PointsInFringe > 1)
	{
		double FrCoor = Sect1D.ArgStart + iLast*Sect1D.ArgStep;

		double LastFrCoor = -1E+23;
		if(AmOfFringeCoor > 0)
		{
			LastFrCoor = FringeCoor[AmOfFringeCoor - 1];
		}
        if(FrCoor != LastFrCoor) 
        {
			FringeContent.push_back(PointsInFringe); AmOfFringes++;
			FringeCoor.push_back(FrCoor);
		}
	}

	AmOfFringeCoor = (int)FringeCoor.size();
	int AmOfFringeContent = (int)FringeContent.size();
	if((AmOfFringeContent == 1) && (AmOfFringeCoor <= 1))
	{
		FringeCoor.erase(FringeCoor.begin(), FringeCoor.end());
		if(iFirst == iLast)
		{
			if(iFirst > 0) iFirst--;
			else if(iLast < (Sect1D.np - 1))
			{
				iLast++;
			}
		}
		FringeCoor.push_back(Sect1D.ArgStart + iFirst*Sect1D.ArgStep);
		FringeCoor.push_back(Sect1D.ArgStart + iLast*Sect1D.ArgStep);
	}

	return 0;
}

//*************************************************************************

int srTGenOptElem::TuneRangeNotDegradingPrec1D(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResize1D& ResizeParam)
{
	double RelPrecTolForRangeSteering = 0.0007; //0.0014; //0.002; // Set or deduce from an external variable !
//OC: to scale with precision
    RelPrecTolForRangeSteering /= ParPrecWfrPropag.PrecFact; //??

	double IncrDecrCoef = 1.07; // To steer
	double CheckCoef = 1.2; // To steer
	int MaxIterNum = 100; //40; // To steer
	int AmOfTrendsToEnsureBadExit = 10; // To steer

	double InvIncrDecrCoef = 1./IncrDecrCoef;

	char TreatExOrEz = ChooseTreatExOrEzBasedOnMax(Sect1D);

	float PrevEstPrec = (float)(1.E+23), EstPrec = (float)(1.E+23);
	double pmPrev = ResizeParam.pm;
	char ProbablyBadExit = 0;
	int BadExitExtraTrendsCount = 0;
	int result = 0;

	for(int i=0; i<MaxIterNum; i++)
	{
		if(result = srYield.Check()) return result;

		PrevEstPrec = EstPrec;
		//if(result = FindRelPrecForRangeOverWfr1D(Sect1D, ResizeParam, CheckCoef, TreatExOrEz, EstPrec)) return result;
		if(result = FindRelPrecForRangeOverWfr1D(Sect1D, ParPrecWfrPropag, ResizeParam, CheckCoef, TreatExOrEz, EstPrec)) return result;

		if(EstPrec < RelPrecTolForRangeSteering)
		{
			ProbablyBadExit = 0; BadExitExtraTrendsCount = 0;
			if((PrevEstPrec > RelPrecTolForRangeSteering) && (i > 0)) break;

			pmPrev = ResizeParam.pm;
			ResizeParam.pm *= InvIncrDecrCoef;
		}
		else
		{
			if(PrevEstPrec < RelPrecTolForRangeSteering)
			{
				ProbablyBadExit = 0; BadExitExtraTrendsCount = 0;
				ResizeParam.pm = pmPrev; break;
			}
			else
			{
				if(EstPrec <= PrevEstPrec)
				{
					ProbablyBadExit = 0; BadExitExtraTrendsCount = 0;
					double pmPrevOld = pmPrev;
					pmPrev = ResizeParam.pm;
					if(ResizeParam.pm >= pmPrevOld)
					{
						ResizeParam.pm *= IncrDecrCoef;
					}
					else
					{
						ResizeParam.pm *= InvIncrDecrCoef;
					}
				}
				else
				{
					if(!ProbablyBadExit)
					{
						double pmPrevOld = pmPrev;
						pmPrev = ResizeParam.pm;
						if(ResizeParam.pm >= pmPrevOld)
						{
							ResizeParam.pm *= IncrDecrCoef;
						}
						else
						{
							ResizeParam.pm *= InvIncrDecrCoef;
						}
						ProbablyBadExit = 1; BadExitExtraTrendsCount++;
					}
					else
					{
						ProbablyBadExit = 1; BadExitExtraTrendsCount++;

						if(BadExitExtraTrendsCount > AmOfTrendsToEnsureBadExit)
						{
							ResizeParam.pm = pmPrev;
							break;
						}
					}
				}
			}
		}
	}
	if(ProbablyBadExit)
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_MAY_BE_LOW);
	}

	return 0;
}

//*************************************************************************

//int srTGenOptElem::FindRelPrecForRangeOverWfr1D(srTRadSect1D& Sect1D, srTRadResize1D& BasicResPar, double IncrDecrCoef, char TreatExOrEz, float& OutRelPrec)
int srTGenOptElem::FindRelPrecForRangeOverWfr1D(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResize1D& BasicResPar, double IncrDecrCoef, char TreatExOrEz, float& OutRelPrec)
{
	double RelRangeToleranceForIntensity = 0.004; //0.002; 
 	//The above number identifies the wavefront part which is taken into account 
	//at the estimation of the precision
	RelRangeToleranceForIntensity /= ParPrecWfrPropag.PrecFact; //??

	int result;
	srTRadSect1D SectDpl;
	if(result = Sect1D.SetupDupl(SectDpl)) return result;

	char ResizePropagOrder;
	if((BasicResPar.pm < 1.) && (BasicResPar.pd < 1.))
	{
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result;
		if(result = RadResizeGen1D(SectDpl, BasicResPar)) return result;
		ResizePropagOrder = 1;
	}
	else
	{
		if(result = RadResizeGen1D(SectDpl, BasicResPar)) return result;
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result;
		ResizePropagOrder = 2;
	}

	srTRadSect1D SectDiffRange;
	if(result = Sect1D.SetupDupl(SectDiffRange)) return result;

	if(ResizePropagOrder == 1)
	{
		srTRadResize1D ResParDiffRange;
		ResParDiffRange.pm = IncrDecrCoef;
		if(result = RadResizeGen1D(SectDiffRange, ResParDiffRange)) return result;
		if(result = PropagateRadiationSimple1D(&SectDiffRange)) return result;
		if(result = RadResizeGen1D(SectDiffRange, BasicResPar)) return result;
	}
	else
	{
		srTRadResize1D ResParDiffRange = BasicResPar;
		ResParDiffRange.pm *= IncrDecrCoef;
		if(result = RadResizeGen1D(SectDiffRange, ResParDiffRange)) return result;
		if(result = PropagateRadiationSimple1D(&SectDiffRange)) return result;
	}

	double ArgStartMax = (SectDpl.ArgStart > SectDiffRange.ArgStart)? SectDpl.ArgStart : SectDiffRange.ArgStart;
	long npMin = (SectDpl.np < SectDiffRange.np)? SectDpl.np : SectDiffRange.np;

	long iStart = IntegerOffsetCoord(SectDpl.ArgStart, SectDpl.ArgStep, ArgStartMax);
	long iStartDif = IntegerOffsetCoord(SectDiffRange.ArgStart, SectDiffRange.ArgStep, ArgStartMax);

	float *tNorm, *tDif;
	if(TreatExOrEz == 'x')
	{
		tNorm = SectDpl.pEx + (iStart << 1); 
		tDif = SectDiffRange.pEx + (iStartDif << 1); 
	}
	else
	{
		tNorm = SectDpl.pEz + (iStart << 1); 
		tDif = SectDiffRange.pEz + (iStartDif << 1); 
	}

	long iLeftThresholdBorder, iRightThresholdBorder;
	if(SectDpl.np < SectDiffRange.np)
		FindIntensityBorders1D(SectDpl, TreatExOrEz, RelRangeToleranceForIntensity, iLeftThresholdBorder, iRightThresholdBorder);
	else 
		FindIntensityBorders1D(SectDiffRange, TreatExOrEz, RelRangeToleranceForIntensity, iLeftThresholdBorder, iRightThresholdBorder);

	long AmOfWorkPo = iRightThresholdBorder - iLeftThresholdBorder + 1;

	double SumRe = 0., SumInt = 0.; //, fMax = -1.E+23;
	for(long i=0; i<npMin; i++)
	{
		if((i >= iLeftThresholdBorder) && (i <= iRightThresholdBorder)) 
		{
			double ReNorm = *tNorm, ImNorm = *(tNorm+1);
			double IntNorm = ReNorm*ReNorm + ImNorm*ImNorm;
			double ReDif = *tDif, ImDif = *(tDif+1);
			double IntDif = ReDif*ReDif + ImDif*ImDif;
			double DeltaRe = IntNorm - IntDif;
			
			SumRe += DeltaRe*DeltaRe;
			SumInt += IntNorm;
		}
		tNorm += 2; tDif += 2;
	}
	double Favg = SumInt/AmOfWorkPo;

	//OC170304
	if(Favg == 0)
	{
		OutRelPrec = (float)0.;
	}
	else
	{
		OutRelPrec = float(sqrt(SumRe/AmOfWorkPo)/Favg);
	}
	return 0;
}

//*************************************************************************

//int srTGenOptElem::FindPostResizeForRange1D(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam)
int srTGenOptElem::FindPostResizeForRange1D(srTRadSect1D& Sect1D, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResize1D& ResizeParam)
{
	double RelRangeToleranceForIntensity = 0.0015;//OC//0.007; // Steer well or set from external variable !
//OC: to scale with precision
	RelRangeToleranceForIntensity /= ParPrecWfrPropag.PrecFact; //??
	
	double MinPostResize = 0.09; // To steer
	if(Sect1D.ThereIsUnderSamplingIn2D) MinPostResize = 0.03; // To steer (to allow micro-focusing)

	//const double RelTolForE = sqrt(RelRangeToleranceForIntensity);

	int result;
	srTRadSect1D SectDpl;
	if(result = Sect1D.SetupDupl(SectDpl)) return result;
	if((ResizeParam.pm < 1.) || (ResizeParam.pd < 1.))
	{
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result; // To implement in all elements
		if(result = RadResizeGen1D(SectDpl, ResizeParam)) return result;
	}
	else
	{
		if(result = RadResizeGen1D(SectDpl, ResizeParam)) return result;
		if(result = PropagateRadiationSimple1D(&SectDpl)) return result; // To implement in all elements
	}

	//long MaxAbsReExInd;
	long long MaxAbsReExInd;
	float MaxAbsReEx;
	char TreatExOrEz;

	float MaxAbsEx, MaxAbsEz;
	//long IndMaxAbsEx, IndMaxAbsEz;
	long long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(SectDpl, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	if(MaxAbsEx > MaxAbsEz)
	{
		TreatExOrEz = 'x'; MaxAbsReEx = MaxAbsEx; MaxAbsReExInd = IndMaxAbsEx;
	}
	else
	{
		TreatExOrEz = 'z'; MaxAbsReEx = MaxAbsEz; MaxAbsReExInd = IndMaxAbsEz;
	}

	long iLeftThresholdBorder = 0, iRightThresholdBorder = SectDpl.np - 1;
	FindIntensityBorders1D(SectDpl, TreatExOrEz, RelRangeToleranceForIntensity, iLeftThresholdBorder, iRightThresholdBorder);

	srTRadResize1D PostResize;
	if((iLeftThresholdBorder > 3) && (iRightThresholdBorder < SectDpl.np - 4))
	{
		long iMid = (SectDpl.np >> 1);
		double pmLeft = double(iMid - iLeftThresholdBorder)/double(iMid);
		double pmRight = double(iRightThresholdBorder - iMid)/double(iMid);
		double pmMax = (pmLeft > pmRight)? pmLeft : pmRight;
		if(pmMax < 0.98) PostResize.pm = pmMax;
	}
	if(PostResize.pm < MinPostResize) PostResize.pm = MinPostResize;

	ResizeParam = PostResize;
	return 0;
}

//*************************************************************************

void srTGenOptElem::FindThresholdBorders(srTRadSect1D& Sect1D, double AbsThresh, char TreatExOrEz, long& iLeftThreshBorder, long& iRightThreshBorder)
{
	long np_mi_1 = Sect1D.np - 1;
	iLeftThreshBorder = -1; iRightThreshBorder = Sect1D.np;

	float *tDir = (TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz;
	float *tInv = tDir + (np_mi_1 << 1);

	for(int i=0; i<Sect1D.np; i++)
	{
		if(iLeftThreshBorder == -1)
		{
			if(::fabs(*tDir) > AbsThresh)
			{
				iLeftThreshBorder = i;
			}
		}
		if(iRightThreshBorder == Sect1D.np)
		{
			if(::fabs(*tInv) > AbsThresh)
			{
				iRightThreshBorder = np_mi_1 - i;
			}
		}
		if((iLeftThreshBorder != -1) && (iRightThreshBorder != Sect1D.np)) break;
		tDir += 2; tInv -= 2;
	}
	if(iLeftThreshBorder == -1) iLeftThreshBorder = 0;
	if(iRightThreshBorder == Sect1D.np) iRightThreshBorder = np_mi_1;
}

//*************************************************************************

//void srTGenOptElem::FindMaximumAbsReE(srTRadSect1D& Sect1D, float& MaxAbsEx, long& IndMaxAbsEx, float& MaxAbsEz, long& IndMaxAbsEz)
void srTGenOptElem::FindMaximumAbsReE(srTRadSect1D& Sect1D, float& MaxAbsEx, long long& IndMaxAbsEx, float& MaxAbsEz, long long& IndMaxAbsEz)
{
	float MaxValEx = (float)(-1.E+23), MaxValEz = (float)(-1.E+23);
	float *tEx = Sect1D.pEx, *tEz = Sect1D.pEz;

	//for(long i=0; i<Sect1D.np; i++)
	for(long long i=0; i<Sect1D.np; i++)
	{
		float TestMaxEx = (float)::fabs(*tEx); tEx += 2;
		if(TestMaxEx > MaxValEx) 
		{
			MaxValEx = TestMaxEx; IndMaxAbsEx = i;
		}
		float TestMaxEz = (float)::fabs(*tEz); tEz += 2;
		if(TestMaxEz > MaxValEz) 
		{
			MaxValEz = TestMaxEz; IndMaxAbsEz = i;
		}
	}
	MaxAbsEx = MaxValEx; MaxAbsEz = MaxValEz;
}

//*************************************************************************
/**
int srTGenOptElem::EstimateMinimalRangeAllowed1D(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam)
{
	const int AmOfIterToAdjustRangeOrResol = 10; // To steer
	int result;

	double &ParamToTune = ResizeParam.pm;
	double PrevParamValue = ParamToTune;

	char ParamValueIsGood = 0, PrevParamValueIsGood;
	for(int i=0; i<AmOfIterToAdjustRangeOrResol; i++)
	{
		double NextParamValue;

		srTRadSect1D SectDupl;
		if(result = Sect1D.SetupDupl(SectDupl)) return result;
		if(result = RadResizeGen1D(SectDupl, ResizeParam)) return result;

		PrevParamValueIsGood = ParamValueIsGood;

		ParamValueIsGood = 1;
		if(i > 0)
		{
			ParamValueIsGood = CheckRangeAndSuggestNextValue1D(SectDupl, ParamToTune, NextParamValue);
		}
		if(ParamValueIsGood)
		{
			if(result = PropagateRadiationSimple1D(&SectDupl)) return result; // To implement in all elements
			ParamValueIsGood = CheckRangeAndSuggestNextValue1D(SectDupl, ParamToTune, NextParamValue);

			if(!ParamValueIsGood) 
			{
				if(PrevParamValueIsGood) 
				{
					double AuxParamToTune = ParamToTune;
					ParamToTune = PrevParamValue;
					PrevParamValue = AuxParamToTune;
					break;
				}
				PrevParamValue = ParamToTune;
				ParamToTune = NextParamValue;
			}
			else
			{
				if((!PrevParamValueIsGood) && (i>0)) break;
				PrevParamValue = ParamToTune;
				ParamToTune = NextParamValue;
			}
		}
		else
		{
			double AuxParamToTune = ParamToTune;
			ParamToTune = PrevParamValue;
			PrevParamValue = AuxParamToTune;

			if(PrevParamValueIsGood) break;
		}
	}

	if(!ParamValueIsGood) {} // Throw Warning message

	return 0;
}
**/
//*************************************************************************
/**
char srTGenOptElem::CheckRangeAndSuggestNextValue1D(srTRadSect1D& Sect1D, double CurrentPm, double& NextPm)
{
	double RelRangeToleranceForIntensity = 0.01; // Set from external variable !
	long MaxAbsReExInd;
	float MaxAbsReEx;

	char TreatExOrEz;
	float MaxAbsEx, MaxAbsEz;
	long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(Sect1D, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	if(MaxAbsEx > MaxAbsEz)
	{
		TreatExOrEz = 'x';
		MaxAbsReEx = MaxAbsEx;
		MaxAbsReExInd = IndMaxAbsEx;
	}
	else
	{
		TreatExOrEz = 'z';
		MaxAbsReEx = MaxAbsEz;
		MaxAbsReExInd = IndMaxAbsEz;
	}

	int AmOfValuesToTreat = 3; // Steer this

	float LeftSumOfLocMaxAbsReEx = 0., RightSumOfLocMaxAbsReEx = 0.;
	long LeftLocMaxInd = -1, RightLocMaxInd = Sect1D.np;
	for(int k=0; k<AmOfValuesToTreat; k++) 
	{
		LeftSumOfLocMaxAbsReEx += ClosestLocMaxAbsReE_FromRight(Sect1D, LeftLocMaxInd, TreatExOrEz);
		RightSumOfLocMaxAbsReEx += ClosestLocMaxAbsReE_FromLeft(Sect1D, RightLocMaxInd, TreatExOrEz);
	}
	float LeftAvgMax = LeftSumOfLocMaxAbsReEx/AmOfValuesToTreat;
	float RightAvgMax = RightSumOfLocMaxAbsReEx/AmOfValuesToTreat;
	float EdgeMaxVal = (LeftAvgMax > RightAvgMax)? LeftAvgMax : RightAvgMax;

	double RatE = EdgeMaxVal/MaxAbsReEx;
	char EdgeValuesAreGood = (RatE*RatE < RelRangeToleranceForIntensity);

	const double EnlargeCoef = 1.2; // Steer this
	const double InvEnlargeCoef = 1./EnlargeCoef;

	NextPm = EdgeValuesAreGood? CurrentPm*InvEnlargeCoef : CurrentPm*EnlargeCoef;

	return EdgeValuesAreGood;
}
**/
//*************************************************************************
/**
float srTGenOptElem::MaximumAbsReEx(srTRadSect1D& Sect1D, long& MaxInd)
{
	float MaxVal = (float)(-1.E+23);
	float *tEx = Sect1D.pEx;
	for(long i=0; i<Sect1D.np; i++)
	{
		float TestMax = (float)::fabs(*tEx); tEx += 2;
		if(TestMax > MaxVal) 
		{
			MaxVal = TestMax; MaxInd = i;
		}
	}
	return MaxVal;
}
**/
//*************************************************************************
/**
float srTGenOptElem::ClosestLocMaxAbsReE_FromRight(srTRadSect1D& Sect1D, long& LocMaxInd, char TreatExOrEz)
{
	float MaxVal = (float)(-1.E+23), MinVal = (float)(1.E+23), PrevVal = (float)(-1.E+23);
	long CurLocMax_i, CurLocMin_i;

	float *tEx = (TreatExOrEz == 'x')? (Sect1D.pEx + ((++LocMaxInd) << 1)) : (Sect1D.pEz + ((++LocMaxInd) << 1));

	for(long i=0; i<(Sect1D.np - LocMaxInd); i++)
	{
		float TestVal = *tEx; tEx += 2;
		if(TestVal > MaxVal) 
		{
			MaxVal = TestVal; CurLocMax_i = i;
		}
		if(TestVal < MinVal) 
		{
			MinVal = TestVal; CurLocMin_i = i;
		}

		long i_mi_1 = i - 1;
		char NormalBreak = (((MinVal < 0.) && (MaxVal > 0.)) && ((CurLocMax_i == i_mi_1) || (CurLocMin_i == i_mi_1)));
		char BreakDueToZero = ((TestVal == PrevVal) && (MaxVal == 0.) && (MinVal == 0.));
		if(NormalBreak || BreakDueToZero) break;

		PrevVal = TestVal;
	}

	float RetVal;
	float AbsMaxVal = (float)::fabs(MaxVal), AbsMinVal = (float)::fabs(MinVal);
	if(AbsMaxVal > AbsMinVal)
	{
		RetVal = AbsMaxVal; 
		LocMaxInd += CurLocMax_i;
	}
	else
	{
		RetVal = AbsMinVal;
		LocMaxInd += CurLocMin_i;
	}
	return RetVal;
}
**/
//*************************************************************************
/**
float srTGenOptElem::ClosestLocMaxAbsReE_FromLeft(srTRadSect1D& Sect1D, long& LocMaxInd, char TreatExOrEz)
{
	float MaxVal = (float)(-1.E+23), MinVal = (float)(1.E+23), PrevVal = (float)(-1.E+23);
	long CurLocMax_i, CurLocMin_i;

	float *tEx = (TreatExOrEz == 'x')? (Sect1D.pEx + ((--LocMaxInd) << 1)) : (Sect1D.pEz + ((--LocMaxInd) << 1));

	for(long i=0; i<=LocMaxInd; i++)
	{
		float TestVal = *tEx; tEx -= 2;
		if(TestVal > MaxVal) 
		{
			MaxVal = TestVal; CurLocMax_i = i;
		}
		if(TestVal < MinVal) 
		{
			MinVal = TestVal; CurLocMin_i = i;
		}

		long i_mi_1 = i - 1;
		char NormalBreak = (((MinVal < 0.) && (MaxVal > 0.)) && ((CurLocMax_i == i_mi_1) || (CurLocMin_i == i_mi_1)));
		char BreakDueToZero = ((TestVal == PrevVal) && (MaxVal == 0.) && (MinVal == 0.));
		if(NormalBreak || BreakDueToZero) break;

		PrevVal = TestVal;
	}

	float RetVal;
	float AbsMaxVal = (float)::fabs(MaxVal), AbsMinVal = (float)::fabs(MinVal);
	if(AbsMaxVal > AbsMinVal)
	{
		RetVal = AbsMaxVal; 
		LocMaxInd += CurLocMax_i;
	}
	else
	{
		RetVal = AbsMinVal;
		LocMaxInd += CurLocMin_i;
	}
	return RetVal;
}
**/
//*************************************************************************
/**
int srTGenOptElem::CheckResolution1D(srTRadSect1D& Sect1D, char& ResolCheckResult, srTAuxTestIntegValues& TestIntegVal, srTRelAndAbsTolerance& RelAndAbsTol)
{
// ResolCheckResult: -1- Error/warning, 0- Undersampling, 1- Good Sampling, 2- Oversampling, 
	double RelToleranceForIntensity = RelAndAbsTol.RelTol;
	double AbsZeroToleranceForIntensity = RelAndAbsTol.AbsTol;
	double IncrDecrCoef = 2.;
	double CoefFromIntToField = 2.; // To steer
	double RelToleranceForField = CoefFromIntToField*RelToleranceForIntensity;
	double AbsZeroToleranceForField = CoefFromIntToField*AbsZeroToleranceForIntensity;
	int result;

	srTRadResize1D DecrResResize;
	DecrResResize.pm = 1.; DecrResResize.pd = 1./IncrDecrCoef;
	srTRadResize1D IncrResResize;
	IncrResResize.pm = 1.; IncrResResize.pd = IncrDecrCoef;
	srTRadResize1D IncrResResize2;
	IncrResResize.pm = 1.; IncrResResize2.pd = IncrDecrCoef*IncrDecrCoef;

	srTRadSect1D SectDecrRes;
	if(result = Sect1D.SetupDupl(SectDecrRes)) return result;
	if(result = RadResizeGen1D(SectDecrRes, DecrResResize)) return result;
	srTRadSect1D SectIncrRes;
	if(result = Sect1D.SetupDupl(SectIncrRes)) return result;
	if(result = RadResizeGen1D(SectIncrRes, IncrResResize)) return result;
	srTRadSect1D SectIncrRes2;
	if(result = Sect1D.SetupDupl(SectIncrRes2)) return result;
	if(result = RadResizeGen1D(SectIncrRes2, IncrResResize2)) return result;

	float MaxAbsEx, MaxAbsEz;
	long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(Sect1D, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	char TreatExOrEz = (MaxAbsEx > MaxAbsEz)? 'x' : 'z';
	float MaxAbsE = (MaxAbsEx > MaxAbsEz)? MaxAbsEx : MaxAbsEz;

	float HalfDecrResIntegVal, HalfNormIntegVal, HalfIncrResIntegVal, HalfIncrResIntegVal2;

	float DecrResIntegVal = IntegrateElField1D(SectDecrRes, TreatExOrEz, HalfDecrResIntegVal);
	float NormIntegVal = IntegrateElField1D(Sect1D, TreatExOrEz, HalfNormIntegVal);
	float IncrResIntegVal = IntegrateElField1D(SectIncrRes, TreatExOrEz, HalfIncrResIntegVal);
	float IncrResIntegVal2 = IntegrateElField1D(SectIncrRes2, TreatExOrEz, HalfIncrResIntegVal2);

	float AbsDecrResIntegVal = (float)::fabs(DecrResIntegVal), AbsNormIntegVal = (float)::fabs(NormIntegVal);
	float AbsIncrResIntegVal = (float)::fabs(IncrResIntegVal), AbsIncrResIntegVal2 = (float)::fabs(IncrResIntegVal2);

	float AbsDifDecr = (float)::fabs(DecrResIntegVal - NormIntegVal);
	float AbsDifIncr = (float)::fabs(IncrResIntegVal - NormIntegVal);
	float AbsDifIncr2 = (float)::fabs(IncrResIntegVal2 - IncrResIntegVal);

	float NormIntegVal_RelTol = (float)(AbsNormIntegVal*RelToleranceForField);

// Logics
	char ClearOverSampling = ((AbsDifIncr2 < AbsDifIncr) && (AbsDifIncr < AbsDifDecr) && (AbsDifDecr < NormIntegVal_RelTol));

	char OverOrNormSamplingInCaseZero = ((AbsIncrResIntegVal2 < AbsIncrResIntegVal) && (AbsIncrResIntegVal < AbsNormIntegVal));
	OverOrNormSamplingInCaseZero = OverOrNormSamplingInCaseZero && (AbsNormIntegVal < AbsZeroToleranceForField*(::fabs(HalfNormIntegVal)));

	char ClearNormSampling = ((AbsDifIncr2 < NormIntegVal_RelTol) && (AbsDifIncr < NormIntegVal_RelTol));

	char ClearUnderSampling = ((AbsDifIncr2 < NormIntegVal_RelTol) && (AbsDifIncr > NormIntegVal_RelTol));

	char HorribleUnderSamplingButStillHopeToRecover = (AbsDifIncr2 < 3.*NormIntegVal_RelTol);

	if(ClearOverSampling) ResolCheckResult = 2;
	else if(ClearNormSampling || OverOrNormSamplingInCaseZero) ResolCheckResult = 1;
	else if(ClearUnderSampling || HorribleUnderSamplingButStillHopeToRecover) ResolCheckResult = 0;
	else ResolCheckResult = -1;

	TestIntegVal.DecrResVal = DecrResIntegVal;
	TestIntegVal.NormVal = NormIntegVal;
	TestIntegVal.IncrResVal = IncrResIntegVal;
	TestIntegVal.IncrRes2Val = IncrResIntegVal2;

	return 0;
}
**/
//*************************************************************************
/**
int srTGenOptElem::CheckAndSuggestNextValue1D(srTRadSect1D& Sect1D, char RangeOrResol, srTAuxTestValues& TestValues, srTRelAndAbsTolerance& RelAndAbsTol)
{
	double RelToleranceForIntensity = RelAndAbsTol.RelTol;
	double AbsZeroToleranceForIntensity = RelAndAbsTol.AbsTol;
	double IncrDecrCoef = RelAndAbsTol.AuxCoef1;
	
	int result;

	char FirstPass = (TestValues.CurrentTestIntegVal.NormVal == 0.)? 1 : 0;
	if(FirstPass)
	{
		char ResolCheckResult = 0;
		srTAuxTestIntegValues TestIntegVal;

		if(result = CheckResolution1D(Sect1D, ResolCheckResult, TestIntegVal, RelAndAbsTol)) return result;

		TestValues.CurrentTestIntegVal = TestIntegVal;
		TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
		TestValues.SuggestedParamValue = (ResolCheckResult < 1)? TestValues.CurrentParamValue*IncrDecrCoef : TestValues.CurrentParamValue/IncrDecrCoef;
		TestValues.CurrentParamValueIsGood = (ResolCheckResult > 0);
		TestValues.CurrentResolCheckResult = ResolCheckResult;
		TestValues.ExitLoop = 0;
		return 0;
	}

	char ResolCheckResult = 0;
	srTAuxTestIntegValues TestIntegVal;

	if(result = CheckResolution1D(Sect1D, ResolCheckResult, TestIntegVal, RelAndAbsTol)) return result;

	if(ResolCheckResult <= 0)
	{
		TestValues.PrevParamValueIsGood = TestValues.CurrentParamValueIsGood;
		TestValues.CurrentParamValueIsGood = 0;
				
		TestValues.PrevResolCheckResult = TestValues.CurrentResolCheckResult;
		TestValues.CurrentResolCheckResult = ResolCheckResult;

		if(TestValues.CurrentResolCheckResult <= 0)
		{
			if(!TestValues.NeedToChangeStrategy)
			{
				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue *= IncrDecrCoef;

				TestValues.ExitLoop = 0;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;

				TestValues.ExitLoop = 1;
			}
		}
		else
		{
			double AuxParVal = TestValues.CurrentParamValue;
			TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
			TestValues.SuggestedParamValue = TestValues.PrevParamValue;
			TestValues.PrevParamValue = AuxParVal;

			TestValues.ExitLoop = 1;
		}
	}
	else
	{
		TestValues.PrevParamValueIsGood = TestValues.CurrentParamValueIsGood;
		TestValues.CurrentParamValueIsGood = 1;

		TestValues.PrevResolCheckResult = TestValues.CurrentResolCheckResult;
		TestValues.CurrentResolCheckResult = ResolCheckResult;

		if(TestValues.CurrentResolCheckResult <= 0)
		{
			if(!TestValues.NeedToChangeStrategy)
			{
				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;
			}

			TestValues.ExitLoop = 1;
		}
		else
		{
			if(!TestValues.NeedToChangeStrategy)
			{
				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue /= IncrDecrCoef;

				TestValues.ExitLoop = 0;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;

				TestValues.ExitLoop = 1;
			}
		}
	}

	TestValues.PrevTestIntegVal = TestValues.CurrentTestIntegVal; TestValues.CurrentTestIntegVal = TestIntegVal;

	return 0;
}
**/
//*************************************************************************
/**
int srTGenOptElem::CheckAndSuggestNextValueWithProp1D(srTRadSect1D& Sect1D, char RangeOrResol, srTAuxTestValues& TestValues, srTRelAndAbsTolerance& RelAndAbsTol)
{
	int result;
	double IncrDecrCoef = RelAndAbsTol.AuxCoef1;

	char ResolCheckResult;
	srTAuxTestIntegValues TestIntegVal;

	char FirstPass = (TestValues.CurrentTestIntegVal.NormVal == 0.)? 1 : 0;
	if(result = CheckPrecAtPropag1D(Sect1D, RangeOrResol, ResolCheckResult, TestIntegVal, RelAndAbsTol)) return result;

	if(FirstPass)
	{
		TestValues.CurrentTestIntegVal = TestIntegVal;
		TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
		TestValues.SuggestedParamValue = (ResolCheckResult < 1)? TestValues.CurrentParamValue*IncrDecrCoef : TestValues.CurrentParamValue/IncrDecrCoef;
		TestValues.CurrentParamValueIsGood = (ResolCheckResult > 0);
		TestValues.CurrentResolCheckResult = ResolCheckResult;
		TestValues.ExitLoop = 0;
		return 0;
	}

	if(ResolCheckResult <= 0)
	{
		TestValues.PrevParamValueIsGood = TestValues.CurrentParamValueIsGood;
		TestValues.CurrentParamValueIsGood = 0;
				
		TestValues.PrevResolCheckResult = TestValues.CurrentResolCheckResult;
		TestValues.CurrentResolCheckResult = ResolCheckResult;

		if(TestValues.CurrentResolCheckResult <= 0)
		{
			if(!TestValues.NeedToChangeStrategy)
			{
				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue *= IncrDecrCoef;

				TestValues.ExitLoop = 0;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;

				TestValues.ExitLoop = 1;
				TestValues.ExitSuccessful = 1;
			}
		}
		else
		{
			double AuxParVal = TestValues.CurrentParamValue;
			TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
			TestValues.SuggestedParamValue = TestValues.PrevParamValue;
			TestValues.PrevParamValue = AuxParVal;

			TestValues.ExitLoop = 1;
		}
	}
	else
	{
		if(TestValues.CurrentResolCheckResult <= 0)
		{
			TestValues.PrevParamValueIsGood = TestValues.CurrentParamValueIsGood;
			TestValues.CurrentParamValueIsGood = 1;
	
			TestValues.PrevResolCheckResult = TestValues.CurrentResolCheckResult;
			TestValues.CurrentResolCheckResult = ResolCheckResult;

			if(!TestValues.NeedToChangeStrategy)
			{
				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;
			}

			TestValues.ExitLoop = 1;
			TestValues.ExitSuccessful = 1;
		}
		else
		{
			TestValues.PrevParamValueIsGood = TestValues.CurrentParamValueIsGood;
			TestValues.CurrentParamValueIsGood = 1;

			TestValues.PrevResolCheckResult = TestValues.CurrentResolCheckResult;
			TestValues.CurrentResolCheckResult = ResolCheckResult;

			if(!TestValues.NeedToChangeStrategy)
			{
				if(TestValues.PrevParamValue < TestValues.CurrentParamValue)
				{
				}

				TestValues.PrevParamValue = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue /= IncrDecrCoef;

				TestValues.ExitLoop = 0;
			}
			else
			{
				double AuxParVal = TestValues.CurrentParamValue;
				TestValues.CurrentParamValue = TestValues.SuggestedParamValue;
				TestValues.SuggestedParamValue = TestValues.PrevParamValue;
				TestValues.PrevParamValue = AuxParVal;

				TestValues.ExitLoop = 1;
				TestValues.ExitSuccessful = 1;
			}
		}
	}
	return 0;
}
**/
//*************************************************************************
/**
int srTGenOptElem::CheckPrecAtPropag1D(srTRadSect1D& Sect1D, char RangeOrResol, char& ResolCheckResult, srTAuxTestIntegValues& TestIntegVal, srTRelAndAbsTolerance& RelAndAbsTol)
{
// ResolCheckResult: -1- Error/warning, 0- Undersampling, 1- Good Sampling, 2- Oversampling, 
	double RelToleranceForIntensity = RelAndAbsTol.RelTol;
	double AbsZeroToleranceForIntensity = RelAndAbsTol.AbsTol;
	double IncrDecrCoef = 2.;
	double CoefFromIntToField = 2.; // To steer

	double RelToleranceForField = CoefFromIntToField*RelToleranceForIntensity;
	double AbsZeroToleranceForField = CoefFromIntToField*AbsZeroToleranceForIntensity;

	float MaxAbsEx, MaxAbsEz;
	long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(Sect1D, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	char TreatExOrEz = (MaxAbsEx > MaxAbsEz)? 'x' : 'z';
	float MaxAbsE = (MaxAbsEx > MaxAbsEz)? MaxAbsEx : MaxAbsEz;

	float DecrResIntegVal, NormIntegVal, IncrResIntegVal, IncrResIntegVal2;
	float HalfDecrResIntegVal, HalfNormIntegVal, HalfIncrResIntegVal, HalfIncrResIntegVal2;
	int result = 0;

	if(TestIntegVal.DecrResValIsNeeded)
	{
		if(result = IntegAfterPropag1D(Sect1D, RangeOrResol, -1, TreatExOrEz, DecrResIntegVal, HalfDecrResIntegVal, RelAndAbsTol)) return result;
	}
	else
	{
		DecrResIntegVal = TestIntegVal.DecrResVal; HalfDecrResIntegVal = TestIntegVal.HalfDecrResVal;
	}
	if(TestIntegVal.NormValIsNeeded)
	{
		if(result = IntegAfterPropag1D(Sect1D, RangeOrResol, 0, TreatExOrEz, NormIntegVal, HalfNormIntegVal, RelAndAbsTol)) return result;
	}
	else
	{
		NormIntegVal = TestIntegVal.DecrResVal; HalfNormIntegVal = TestIntegVal.HalfNormVal;
	}
	if(TestIntegVal.IncrResValIsNeeded)
	{
		if(result = IntegAfterPropag1D(Sect1D, RangeOrResol, 1, TreatExOrEz, IncrResIntegVal, HalfIncrResIntegVal, RelAndAbsTol)) return result;
	}
	else
	{
		IncrResIntegVal = TestIntegVal.IncrResVal; HalfIncrResIntegVal = TestIntegVal.HalfIncrResVal;
	}
	if(TestIntegVal.IncrRes2ValIsNeeded)
	{
		if(result = IntegAfterPropag1D(Sect1D, RangeOrResol, 2, TreatExOrEz, IncrResIntegVal2, HalfIncrResIntegVal2, RelAndAbsTol)) return result;
	}
	else
	{
		IncrResIntegVal2 = TestIntegVal.IncrRes2Val; HalfIncrResIntegVal2 = TestIntegVal.HalfIncrRes2Val;
	}

	float AbsDecrResIntegVal = (float)::fabs(DecrResIntegVal), AbsNormIntegVal = (float)::fabs(NormIntegVal);
	float AbsIncrResIntegVal = (float)::fabs(IncrResIntegVal), AbsIncrResIntegVal2 = (float)::fabs(IncrResIntegVal2);

	float AbsDifDecr = (float)::fabs(DecrResIntegVal - NormIntegVal);
	float AbsDifIncr = (float)::fabs(IncrResIntegVal - NormIntegVal);
	float AbsDifIncr2 = (float)::fabs(IncrResIntegVal2 - IncrResIntegVal);

	float NormIntegVal_RelTol = (float)(AbsNormIntegVal*RelToleranceForField);

// Logics
	char ClearOverSampling = ((AbsDifIncr2 < AbsDifIncr) && (AbsDifIncr < AbsDifDecr) && (AbsDifDecr < NormIntegVal_RelTol));

	char OverOrNormSamplingInCaseZero = ((AbsIncrResIntegVal2 < AbsIncrResIntegVal) && (AbsIncrResIntegVal < AbsNormIntegVal));
	OverOrNormSamplingInCaseZero = OverOrNormSamplingInCaseZero && (AbsNormIntegVal < AbsZeroToleranceForField*(::fabs(HalfNormIntegVal)));

	char ClearNormSampling = ((AbsDifIncr2 < NormIntegVal_RelTol) && (AbsDifIncr < NormIntegVal_RelTol));
	char ClearUnderSampling = ((AbsDifIncr2 < NormIntegVal_RelTol) && (AbsDifIncr > NormIntegVal_RelTol));
	char HorribleUnderSamplingButStillHopeToRecover = (AbsDifIncr2 < 3.*NormIntegVal_RelTol);

	if(ClearOverSampling) ResolCheckResult = 2;
	else if(ClearNormSampling || OverOrNormSamplingInCaseZero) ResolCheckResult = 1;
	else if(ClearUnderSampling || HorribleUnderSamplingButStillHopeToRecover) ResolCheckResult = 0;
	else ResolCheckResult = -1;

	TestIntegVal.DecrResVal = DecrResIntegVal; TestIntegVal.HalfDecrResVal = HalfDecrResIntegVal; 
	TestIntegVal.NormVal = NormIntegVal; TestIntegVal.HalfNormVal = HalfNormIntegVal; 
	TestIntegVal.IncrResVal = IncrResIntegVal; TestIntegVal.HalfIncrResVal = HalfIncrResIntegVal; 
	TestIntegVal.IncrRes2Val = IncrResIntegVal2; TestIntegVal.HalfIncrRes2Val = HalfIncrResIntegVal2; 

	return 0;
}
**/
//*************************************************************************
/**
int srTGenOptElem::IntegAfterPropag1D(srTRadSect1D& Sect1D, char RangeOrResol, char DecrOrIncr, char TreatExOrEz, float& ResIntegVal, float& HalfIntegVal, srTRelAndAbsTolerance& RelAndAbsTol)
{
// DecrOrIncr: -1 - decr, 0 - norm, 1 - incr, 2 - incr2
	double IncrDecrCoef = RelAndAbsTol.AuxCoef1;
	int result;

	srTRadResize1D ResResize;
	if(RangeOrResol == 1) // Range
	{
		ResResize.pd = 1.;
		if(DecrOrIncr == -1)
		{
			ResResize.pm = 1./IncrDecrCoef; 
		}
		else if(DecrOrIncr == 0)
		{
			ResResize.pm = 1.; 
		}
		else if(DecrOrIncr == 1)
		{
			ResResize.pm = IncrDecrCoef; 
		}
		else if(DecrOrIncr == 2)
		{
			ResResize.pm = IncrDecrCoef*IncrDecrCoef; 
		}
	}
	else
	{
		ResResize.pm = 1.; 
		if(DecrOrIncr == -1)
		{
			ResResize.pd = 1./IncrDecrCoef; 
		}
		else if(DecrOrIncr == 0)
		{
			ResResize.pd = 1.; 
		}
		else if(DecrOrIncr == 1)
		{
			ResResize.pd = IncrDecrCoef; 
		}
		else if(DecrOrIncr == 2)
		{
			ResResize.pd = IncrDecrCoef*IncrDecrCoef; 
		}
	}
	srTRadSect1D SectDiffRes;
	if(result = Sect1D.SetupDupl(SectDiffRes)) return result;
	if((ResResize.pm != 1.) || (ResResize.pd != 1.))
	{
		if(result = RadResizeGen1D(SectDiffRes, ResResize)) return result;
	}
	if(result = PropagateRadiationSimple1D(&SectDiffRes)) return result; // To implement in all elements
	ResIntegVal = IntegrateElField1D(SectDiffRes, TreatExOrEz, HalfIntegVal);
	return 0;
}
**/
//*************************************************************************
