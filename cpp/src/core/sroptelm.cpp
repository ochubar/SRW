/************************************************************************//**
 * File: sroptelm.cpp
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
#include "srsend.h"

#include "sroptapt.h"
#include "sroptdrf.h"
#include "sroptwgr.h"
#include "sroptfoc.h"
#include "sroptsmr.h"
#include "sroptpsh.h"
#include "gmfft.h"
#include "sroptgtr.h"
#include "sroptzps.h"
#include "sroptzp.h"
#include "srmatsta.h"
#include "srmlttsk.h"
#include "srinterf.h"
#include "sropthck.h"
#include "sroptgrat.h"

//*************************************************************************

extern srTYield srYield;

//*************************************************************************

//int srTOptElemSummary::SetupOpticalElement(srTStringVect* pOptElemInfo, srTGenOptElemHndl& OptElemHndl, srTSRWRadStructAccessData* pRad)
//int srTGenOptElem::SetupOpticalElement(srTStringVect* pOptElemInfo, srTGenOptElemHndl& OptElemHndl, srTSRWRadStructAccessData* pRad)
int srTGenOptElem::SetupOpticalElement(srTStringVect* pOptElemInfo, srTDataMD* pExtraData, srTSRWRadStructAccessData* pRad, srTGenOptElemHndl& OptElemHndl)
{
	char* ElemID = (*pOptElemInfo)[0];

	if(!strcmp(ElemID, "Container"))
	{
		//int AmOfElem = atoi((*pOptElemInfo)[1]);
		//OptElemHndl = srTGenOptElemHndl(new srTCompositeOptElem(pOptElemInfo, AmOfElem, pRad));
		OptElemHndl = srTGenOptElemHndl(new srTCompositeOptElem(pOptElemInfo, pRad));
	}

	else if(!strcmp(ElemID, "Drift"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTDriftSpace(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "RectAperture"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTRectAperture(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "CircAperture"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTCircAperture(pOptElemInfo));
	}

	else if(!strcmp(ElemID, "RectObstacle"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTRectObstacle(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "CircObstacle"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTCircObstacle(pOptElemInfo));
	}

	else if(!strcmp(ElemID, "ThinLens"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTThinLens(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "WaveguideRect"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTWaveguideRect(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "ZonePlate"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTZonePlate(pOptElemInfo));
	}
	else if(!strcmp(ElemID, "PlaneGrating"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTGrating(pOptElemInfo));
	}

	else if(!strcmp(ElemID, "ThinGen"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTGenTransmission(pOptElemInfo, pExtraData));
	}

	else if(!strcmp(ElemID, "PhaseShift"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTPhaseShift(pOptElemInfo, pRad));
	}
	else if(!strcmp(ElemID, "ZonePlateSpec"))
	{
		OptElemHndl = srTGenOptElemHndl(new srTZonePlateSpec(pOptElemInfo, pRad));
	}

	else if(!strcmp(ElemID, "Mirror"))
	{
		OptElemHndl = srTGenOptElemHndl(srTMirror::DefineMirror(pOptElemInfo, pExtraData));
	}

	//else if(!strcmp(ElemID, "ThickMirrorGen"))
	//{
	//	OptElemHndl = srTGenOptElemHndl(new srTThickMirrorGen(pOptElemInfo, pExtraData));
	//}
	//else if(!strcmp(ElemID, "ThickMirrorToroid"))
	//{
	//	OptElemHndl = srTGenOptElemHndl(new srTThickMirrorToroid(pOptElemInfo));
	//}
	//else if(!strcmp(ElemID, "SpherMirror"))
	//{
	//	OptElemHndl = srTGenOptElemHndl(new srTSpherMirror(pOptElemInfo));
	//}

// Continue for new Optical Elements

	else return UNKNOWN_OPTICAL_ELEMENT;

	//if(OptElemHndl.rep->ErrorCode != 0) return OptElemHndl.rep->ErrorCode;
	if(((srTGenOptElem*)(OptElemHndl.rep))->ErrorCode != 0) return ((srTGenOptElem*)(OptElemHndl.rep))->ErrorCode;

	return 0;
}

//*************************************************************************

int srTGenOptElem::TraverseRadZXE(srTSRWRadStructAccessData* pRadAccessData)
{
	float *pEx0 = pRadAccessData->pBaseRadX;
	float *pEz0 = pRadAccessData->pBaseRadZ;
	//long PerX = pRadAccessData->ne << 1;
	//long PerZ = PerX*pRadAccessData->nx;
	long long PerX = pRadAccessData->ne << 1;
	long long PerZ = PerX*pRadAccessData->nx;

	srTEFieldPtrs EFieldPtrs;
	srTEXZ EXZ;
	EXZ.z = pRadAccessData->zStart;
	//long izPerZ = 0;
	//long iTotTest = 0; //OCTEST
	long long izPerZ = 0;
	long long iTotTest = 0; //OCTEST

	int result = 0;
	for(int iz=0; iz<pRadAccessData->nz; iz++)
	{
		if(result = srYield.Check()) return result;

		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		EXZ.x = pRadAccessData->xStart;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<pRadAccessData->nx; ix++)
		{
			float *pEx_StartForE = pEx_StartForX + ixPerX;
			float *pEz_StartForE = pEz_StartForX + ixPerX;
			EXZ.e = pRadAccessData->eStart;
			//long iePerE = 0;
			long long iePerE = 0;

			//DEBUG
			//if((iz == ((pRadAccessData->nz)>>1)) && (ix == ((pRadAccessData->nx)>>1)))
			//{
			//	int aha = 1;
			//}
			//END DEBUG

			for(int ie=0; ie<pRadAccessData->ne; ie++)
			{
				if(pEx0 != 0)
				{
					EFieldPtrs.pExRe = pEx_StartForE + iePerE;
					EFieldPtrs.pExIm = EFieldPtrs.pExRe + 1;
				}
				else
				{
					EFieldPtrs.pExRe = 0;
					EFieldPtrs.pExIm = 0;
				}
				if(pEz0 != 0)
				{
					EFieldPtrs.pEzRe = pEz_StartForE + iePerE;
					EFieldPtrs.pEzIm = EFieldPtrs.pEzRe + 1;
				}
				else
				{
					EFieldPtrs.pEzRe = 0;
					EFieldPtrs.pEzIm = 0;
				}

				EXZ.aux_offset = izPerZ + ixPerX + iePerE;
				RadPointModifier(EXZ, EFieldPtrs);

				iTotTest++; //OCTEST

				iePerE += 2;
				EXZ.e += pRadAccessData->eStep;
			}
			ixPerX += PerX;
			EXZ.x += pRadAccessData->xStep;
		}
		izPerZ += PerZ;
		EXZ.z += pRadAccessData->zStep;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::TraverseRad1D(srTRadSect1D* pSect1D)
{
	float *tEx = pSect1D->pEx;
	float *tEz = pSect1D->pEz;

	srTEFieldPtrs EFieldPtrs;
	srTEXZ EXZ;

	EXZ.VsXorZ = pSect1D->VsXorZ;
	EXZ.e = pSect1D->eVal;
	EXZ.x = (pSect1D->VsXorZ == 'x')? pSect1D->ArgStart : pSect1D->OtherCoordVal;
	EXZ.z = (pSect1D->VsXorZ == 'x')? pSect1D->OtherCoordVal : pSect1D->ArgStart;
	double &Arg = (pSect1D->VsXorZ == 'x')? EXZ.x : EXZ.z;

	for(int i=0; i<pSect1D->np; i++)
	{
		EFieldPtrs.pExRe = tEx;
		EFieldPtrs.pExIm = tEx + 1;
		EFieldPtrs.pEzRe = tEz;
		EFieldPtrs.pEzIm = tEz + 1;

		RadPointModifier1D(EXZ, EFieldPtrs);

		if(tEx != 0) tEx += 2; 
		if(tEz != 0) tEz += 2;
		Arg += pSect1D->ArgStep;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::ExtractRadSliceConstE(srTSRWRadStructAccessData* pRadAccessData, long ie, float*& pOutEx, float*& pOutEz, bool forceCopyField)
{// ATTENTION: In the case of single energy, it may simply return pointers from pRadAccessData!!!

	float *pEx0 = pRadAccessData->pBaseRadX;
	float *pEz0 = pRadAccessData->pBaseRadZ;

	if(!forceCopyField)
	{
		if(pRadAccessData->ne == 1)
		{
			if(pOutEx == 0) pOutEx = pEx0; 
			if(pOutEz == 0) pOutEz = pEz0;
			return 0;
		}
	}
	if((ie < 0) || (ie >= pRadAccessData->ne)) return 0;

	//long PerX = pRadAccessData->ne << 1;
	//long PerZ = PerX*pRadAccessData->nx;
	//long izPerZ = 0;
	//long iePerE = ie << 1;

	long long  PerX = pRadAccessData->ne << 1;
	long long PerZ = PerX*pRadAccessData->nx;
	long long izPerZ = 0;
	long long iePerE = ie << 1;

	float *tOutEx = pOutEx, *tOutEz = pOutEz;
	for(int iz=0; iz<pRadAccessData->nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<pRadAccessData->nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pEx = pEx_StartForX + ixPerX_p_iePerE;
			float *pEz = pEz_StartForX + ixPerX_p_iePerE;

			*(tOutEx++) = *(pEx++); 
			*(tOutEx++) = *pEx;

			*(tOutEz++) = *(pEz++); 
			*(tOutEz++) = *pEz;

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::SetupRadSliceConstE(srTSRWRadStructAccessData* pRadAccessData, long ie, float* pInEx, float* pInEz)
{
	float *pEx0 = pRadAccessData->pBaseRadX;
	float *pEz0 = pRadAccessData->pBaseRadZ;
	//long PerX = pRadAccessData->ne << 1;
	//long PerZ = PerX*pRadAccessData->nx;
	//long izPerZ = 0;
	//long iePerE = ie << 1;

	long long PerX = pRadAccessData->ne << 1;
	long long PerZ = PerX*pRadAccessData->nx;
	long long izPerZ = 0;
	long long iePerE = ie << 1;

	float *tInEx = pInEx, *tInEz = pInEz;
	for(int iz=0; iz<pRadAccessData->nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<pRadAccessData->nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pEx = pEx_StartForX + ixPerX_p_iePerE;
			float *pEz = pEz_StartForX + ixPerX_p_iePerE;

			*(pEx++) = *(tInEx++); *pEx = *(tInEx++);
			*(pEz++) = *(tInEz++); *pEz = *(tInEz++);

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::ExtractRadSectVsXorZ(srTSRWRadStructAccessData* pRadAccessData, long ie, long ix_or_iz, char Vs_x_or_z, float* pOutEx, float* pOutEz)
{
	//long Period, InitialOffset, Np;
	//long Two_ne = pRadAccessData->ne << 1, Two_ie = ie << 1;
	long long Period, InitialOffset, Np;
	long long Two_ne = pRadAccessData->ne << 1, Two_ie = ie << 1;

	if(Vs_x_or_z == 'x')
	{
		Period = Two_ne;
		InitialOffset = ix_or_iz*(pRadAccessData->nx*Two_ne) + Two_ie;
		Np = pRadAccessData->nx;
	}
	else
	{
		Period = pRadAccessData->nx*Two_ne;
		InitialOffset = ix_or_iz*Two_ne + Two_ie;
		Np = pRadAccessData->nz;
	}

	float *tEx = pRadAccessData->pBaseRadX + InitialOffset;
	float *tEz = pRadAccessData->pBaseRadZ + InitialOffset;

	float *tOutEx = pOutEx, *tOutEz = pOutEz;

	for(int i=0; i<Np; i++)
	{
		*(tOutEx++) = *tEx; 
		*(tOutEx++) = *(tEx+1); 
		*(tOutEz++) = *tEz; 
		*(tOutEz++) = *(tEz+1); 

		tEx += Period; tEz += Period;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::SetupSectionArraysVsXandZ(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D& SectVsX, srTRadSect1D& SectVsZ)
{
	long LenArrX = pRadAccessData->nx << 1;
	float* SectVsX_pEx = new float[LenArrX];
	if(SectVsX_pEx == 0) return MEMORY_ALLOCATION_FAILURE;
	float* SectVsX_pEz = new float[LenArrX];
	if(SectVsX_pEz == 0) return MEMORY_ALLOCATION_FAILURE;

	long LenArrZ = pRadAccessData->nz << 1;
	float* SectVsZ_pEx = new float[LenArrZ];
	if(SectVsZ_pEx == 0) return MEMORY_ALLOCATION_FAILURE;
	float* SectVsZ_pEz = new float[LenArrZ];
	if(SectVsZ_pEz == 0) return MEMORY_ALLOCATION_FAILURE;

	long iec = 0;
	int result;
	if(result = ExtractRadSectVsXorZ(pRadAccessData, iec, SectVsX.icOtherCoord, 'x', SectVsX_pEx, SectVsX_pEz)) return result;
	if(result = ExtractRadSectVsXorZ(pRadAccessData, iec, SectVsZ.icOtherCoord, 'z', SectVsZ_pEx, SectVsZ_pEz)) return result;

	SectVsX = srTRadSect1D(SectVsX_pEx, SectVsX_pEz, 'x', SectVsX.icOtherCoord, pRadAccessData);
	SectVsX.DeleteArraysAtDestruction = 1;

	SectVsZ = srTRadSect1D(SectVsZ_pEx, SectVsZ_pEz, 'z', SectVsZ.icOtherCoord, pRadAccessData);
	SectVsZ.DeleteArraysAtDestruction = 1;

	return 0;
}

//*************************************************************************

int srTGenOptElem::SetupNewRadStructFromSliceConstE(srTSRWRadStructAccessData* pRadAccessData, long ie, srTSRWRadStructAccessData*& pRadDataSingleE)
{// Only new Electric Field may be allocated (all the rest just points to old data !!!)
	//if(pRadAccessData->ne == 1)
	//{
	//	pRadDataSingleE = pRadAccessData; return 0;
	//}

	pRadDataSingleE = new srTSRWRadStructAccessData();
	if(pRadDataSingleE == 0) return MEMORY_ALLOCATION_FAILURE;

	*pRadDataSingleE = *pRadAccessData;

	//long LenFloatArr = (pRadAccessData->nx*pRadAccessData->nz) << 1;
	long long LenFloatArr = (pRadAccessData->nx*pRadAccessData->nz) << 1;
	pRadDataSingleE->pBaseRadX = new float[LenFloatArr];
	if(pRadDataSingleE->pBaseRadX == 0) return MEMORY_ALLOCATION_FAILURE;
	pRadDataSingleE->pBaseRadZ = new float[LenFloatArr];
	if(pRadDataSingleE->pBaseRadZ == 0) return MEMORY_ALLOCATION_FAILURE;

	pRadDataSingleE->BaseRadWasEmulated = true; //to ensure that the above arrays are deleted by destructor

	const int AmOfMoments = 11;
	//long OffsetMom = 0;
	long long OffsetMom = 0;
	double OffsetPhotEn = 0;

	int result;
	if((ie >= 0) && (ie < pRadAccessData->ne))
	{
		if(result = ExtractRadSliceConstE(pRadAccessData, ie, pRadDataSingleE->pBaseRadX, pRadDataSingleE->pBaseRadZ)) return result;
		OffsetPhotEn = ie*pRadAccessData->eStep;
		//pRadDataSingleE->eStart = pRadAccessData->eStart + ie*pRadAccessData->eStep;
		OffsetMom = AmOfMoments*ie;
	}

	pRadDataSingleE->eStep = 0.;
	pRadDataSingleE->eStart = pRadAccessData->eStart + OffsetPhotEn;
	pRadDataSingleE->ne = 1;

	pRadDataSingleE->pMomX = pRadAccessData->pMomX + OffsetMom;
	pRadDataSingleE->pMomZ = pRadAccessData->pMomZ + OffsetMom;

	//Note that all arrays except pBaseRadX, pBaseRadZ were not allocated in *pRadDataSingleE!:
	pRadDataSingleE->ResAfterWasEmulated = false; //OC13112011
	pRadDataSingleE->ElectronBeamEmulated = 0;
	pRadDataSingleE->PropMatrWasEmulated = false;
	pRadDataSingleE->MomWereEmulated = false;
	pRadDataSingleE->WfrAuxDataWasEmulated = false;

	return 0;
}

//*************************************************************************
/**
int srTGenOptElem::UpdateGenRadStructFromSlicesConstE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* arRadDataSlicesConstE)
{//Compose the Electric Field of the pRadAccessData from the slices ConstE.
 //The slices are assumed to have same dimensions over nx and nz. 
	if((pRadAccessData == 0) || (arRadDataSlicesConstE == 0)) return 0;
	if(pRadAccessData->ne <= 0) return 0;

	float *pEx0 = pRadAccessData->pBaseRadX;
	float *pEz0 = pRadAccessData->pBaseRadZ;
	
	int neCom = pRadAccessData->ne;
	int nxCom = pRadAccessData->nx;
	int nzCom = pRadAccessData->nz;

	long PerX = neCom << 1;
	long PerZ = PerX*nxCom;

	srTSRWRadStructAccessData *t_arRadDataSlicesConstE = arRadDataSlicesConstE;
	for(int ie=0; ie<neCom; ie++)
	{
		long izPerZ = 0;
		long iePerE = ie << 1;
		float *tSliceEx = t_arRadDataSlicesConstE->pBaseRadX;
		float *tSliceEz = (t_arRadDataSlicesConstE++)->pBaseRadZ;

		for(int iz=0; iz<nzCom; iz++)
		{
			float *pEx_StartForX = pEx0 + izPerZ;
			float *pEz_StartForX = pEz0 + izPerZ;
			long ixPerX = 0;

			for(int ix=0; ix<pRadAccessData->nx; ix++)
			{
				long ixPerX_p_iePerE = ixPerX + iePerE;
				float *pEx = pEx_StartForX + ixPerX_p_iePerE;
				float *pEz = pEz_StartForX + ixPerX_p_iePerE;

				*(pEx++) = *(tSliceEx++); *pEx = *(tSliceEx++);
				*(pEz++) = *(tSliceEz++); *pEz = *(tSliceEz++);

				ixPerX += PerX;
			}
			izPerZ += PerZ;
		}
	}
	return 0;
}
**/
//*************************************************************************

int srTGenOptElem::UpdateGenRadStructSliceConstE_Meth_0(srTSRWRadStructAccessData* pRadDataSliceConstE, int ie, srTSRWRadStructAccessData* pRadAccessData)
{//Compose the Electric Field of the pRadAccessData from the slices ConstE.
 //The slices are assumed to have same dimensions over nx and nz. 
	if((pRadAccessData == 0) || (pRadDataSliceConstE == 0) || (ie < 0)) return 0;

	float *pEx0 = pRadAccessData->pBaseRadX;
	float *pEz0 = pRadAccessData->pBaseRadZ;
	
	int neCom = pRadAccessData->ne;
	int nxCom = pRadAccessData->nx;
	int nzCom = pRadAccessData->nz;
	if(neCom <= 0) return 0;

	//long PerX = neCom << 1;
	//long PerZ = PerX*nxCom;
	long long PerX = neCom << 1;
	long long PerZ = PerX*nxCom;

	//long izPerZ = 0;
	//long iePerE = ie << 1;
	long long izPerZ = 0;
	long long iePerE = ie << 1;
	float *tSliceEx = pRadDataSliceConstE->pBaseRadX;
	float *tSliceEz = pRadDataSliceConstE->pBaseRadZ;

	for(int iz=0; iz<nzCom; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<pRadAccessData->nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pEx = pEx_StartForX + ixPerX_p_iePerE;
			float *pEz = pEz_StartForX + ixPerX_p_iePerE;

			*(pEx++) = *(tSliceEx++); *pEx = *(tSliceEx++);
			*(pEz++) = *(tSliceEz++); *pEz = *(tSliceEz++);

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}

	//Update wavefront limits in the main rad. structure:
	if(pRadAccessData->xWfrMin > pRadDataSliceConstE->xWfrMin) pRadAccessData->xWfrMin = pRadDataSliceConstE->xWfrMin;
	if(pRadAccessData->xWfrMax < pRadDataSliceConstE->xWfrMax) pRadAccessData->xWfrMax = pRadDataSliceConstE->xWfrMax;
	if(pRadAccessData->zWfrMin > pRadDataSliceConstE->zWfrMin) pRadAccessData->zWfrMin = pRadDataSliceConstE->zWfrMin;
	if(pRadAccessData->zWfrMax < pRadDataSliceConstE->zWfrMax) pRadAccessData->zWfrMax = pRadDataSliceConstE->zWfrMax;

	//Update wavefront radii in the main rad. structure (making average?):
	double inv_ie_p_1 = 1./(ie + 1);
	pRadAccessData->RobsX = (ie*(pRadAccessData->RobsX) + pRadDataSliceConstE->RobsX)*inv_ie_p_1;
	pRadAccessData->RobsZ = (ie*(pRadAccessData->RobsZ) + pRadDataSliceConstE->RobsZ)*inv_ie_p_1;
	pRadAccessData->RobsXAbsErr = (ie*(pRadAccessData->RobsXAbsErr) + pRadDataSliceConstE->RobsXAbsErr)*inv_ie_p_1;
	pRadAccessData->RobsZAbsErr = (ie*(pRadAccessData->RobsZAbsErr) + pRadDataSliceConstE->RobsZAbsErr)*inv_ie_p_1;

	//assuming that the mesh in all slices is transformed the same way
	//if(ie == (pRadAccessData->ne - 1)) //OC151008 //commented-out
	//{
	//	if(pRadAccessData->xStart != pRadDataSliceConstE->xStart) pRadAccessData->xStart = pRadDataSliceConstE->xStart;
	//	if(pRadAccessData->xStep != pRadDataSliceConstE->xStep) pRadAccessData->xStep = pRadDataSliceConstE->xStep;
	//	if(pRadAccessData->zStart != pRadDataSliceConstE->zStart) pRadAccessData->zStart = pRadDataSliceConstE->zStart;
	//	if(pRadAccessData->zStep != pRadDataSliceConstE->zStep) pRadAccessData->zStep = pRadDataSliceConstE->zStep;
	//}

	return 0;
}


//*************************************************************************

int srTGenOptElem::UpdateGenRadStructSliceConstE_Meth_2(srTSRWRadStructAccessData* pRadDataSliceConstE, int ie, srTSRWRadStructAccessData* pRadAccessData)
{//To implement:
 // Undate the Electric Field of the pRadAccessData from the slice ConstE.
// The slice can have different dimensions over nx and nz, compared to general pRadAccessData. 
// In such a case, maximum nx and nz (over Range and Resolution) are taken and Resizing is performed on pRadAccessData.
// No loss of data / precision here.
// To save memory, use the "realloc" C function.

	return 0;
}

//*************************************************************************

//int srTGenOptElem::UpdateGenRadStructFromSlicesConstE(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* RadDataStructArray)
//{
//	if(pRadAccessData->ne == 1) return 0;
//
//// Compose the Electric Field of the pRadAccessData from the slices ConstE.
//// The slices can have different dimensions over nx and nz. 
//// In such a case, maximum nx and nz (over Range and Resolution) are taken and Resizing is performed on all the structures which need this.
//// No loss of data / precision here.
//// To save memory, use the "realloc" C function.
//
//// Delete all Electric Field arrays in slices here !!!
//
//	return 0;
//}

//*************************************************************************

int srTGenOptElem::RemoveSliceConstE_FromGenRadStruct(srTSRWRadStructAccessData* pRadAccessData, long ie)
{
	if(pRadAccessData->ne == 1) return 0;

// Remove an Electric Field slice ConstE from the general Electric Field arrays.
// Use the "realloc" C function to save memory !!!
// No other data of the pRadAccessData is changed !!!

// This should be used only to save memory at propagation

	return 0;
}

//*************************************************************************

int srTGenOptElem::SetupWfrEdgeCorrData(srTSRWRadStructAccessData* pRadAccessData, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrsForWfrEdgeCorr)
{
	int result;

	double xAbsTol = 0.05*pRadAccessData->xStep;
	double zAbsTol = 0.05*pRadAccessData->zStep;
//--X
	double xWfrMinOffsetFromStart = pRadAccessData->xWfrMin - pRadAccessData->xStart;
	long ixWfrMinLower = long(xWfrMinOffsetFromStart/pRadAccessData->xStep + 1.E-13);
	double xWfrMinLowerMisfit = xWfrMinOffsetFromStart - ixWfrMinLower*pRadAccessData->xStep;

	double xWfrMaxOffsetFromStart = pRadAccessData->xWfrMax - pRadAccessData->xStart;
	long ixWfrMaxLower = long(xWfrMaxOffsetFromStart/pRadAccessData->xStep + 1.E-13);
	double xWfrMaxLowerMisfit = xWfrMaxOffsetFromStart - ixWfrMaxLower*pRadAccessData->xStep;
	
	char xWfrMinIsBetweenMeshPoints = (xWfrMinLowerMisfit > xAbsTol);
	char xWfrMaxIsBetweenMeshPoints = (xWfrMaxLowerMisfit > xAbsTol);
	char xWfrMaxIsSmallerThanDataEnd = (::fabs((pRadAccessData->xStart + pRadAccessData->nx*pRadAccessData->xStep) - pRadAccessData->xWfrMax) > xAbsTol);
	char xWfrCorrNeeded = (xWfrMinIsBetweenMeshPoints || xWfrMaxIsBetweenMeshPoints || xWfrMaxIsSmallerThanDataEnd);

	float dxSt = 0.;
	if(xWfrMinIsBetweenMeshPoints) dxSt = (float)(pRadAccessData->xStep - xWfrMinLowerMisfit);

	float dxFi = 0.;
	if(xWfrMaxIsBetweenMeshPoints) dxFi = (float)(pRadAccessData->xStep - xWfrMaxLowerMisfit);
	else if(xWfrMaxIsSmallerThanDataEnd) dxFi = (float)(pRadAccessData->xStep);
//--Z
	double zWfrMinOffsetFromStart = pRadAccessData->zWfrMin - pRadAccessData->zStart;

	long izWfrMinLower = long(zWfrMinOffsetFromStart/pRadAccessData->zStep + 1.E-13);
	double zWfrMinLowerMisfit = zWfrMinOffsetFromStart - izWfrMinLower*pRadAccessData->zStep;

	double zWfrMaxOffsetFromStart = pRadAccessData->zWfrMax - pRadAccessData->zStart;
	long izWfrMaxLower = long(zWfrMaxOffsetFromStart/pRadAccessData->zStep + 1.E-13);
	double zWfrMaxLowerMisfit = zWfrMaxOffsetFromStart - izWfrMaxLower*pRadAccessData->zStep;
	
	char zWfrMinIsBetweenMeshPoints = (zWfrMinLowerMisfit > zAbsTol);
	char zWfrMaxIsBetweenMeshPoints = (zWfrMaxLowerMisfit > zAbsTol);
	char zWfrMaxIsSmallerThanDataEnd = (::fabs((pRadAccessData->zStart + pRadAccessData->nz*pRadAccessData->zStep) - pRadAccessData->zWfrMax) > zAbsTol);
	char zWfrCorrNeeded = (zWfrMinIsBetweenMeshPoints || zWfrMaxIsBetweenMeshPoints || zWfrMaxIsSmallerThanDataEnd);

	float dzSt = 0.;
	if(zWfrMinIsBetweenMeshPoints) dzSt = (float)(pRadAccessData->zStep - zWfrMinLowerMisfit);

	float dzFi = 0.;
	if(zWfrMaxIsBetweenMeshPoints) dzFi = (float)(pRadAccessData->zStep - zWfrMaxLowerMisfit);
	else if(zWfrMaxIsSmallerThanDataEnd) dzFi = (float)(pRadAccessData->zStep);

//--Gen
	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.UseGivenStartTrValues = 0;
	CGenMathFFT2D FFT2D;
	FFT2D.SetupLimitsTr(FFT2DInfo);

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.HowMany = 2;
	FFT1DInfo.UseGivenStartTrValue = 0;

	if(xWfrCorrNeeded || zWfrCorrNeeded)
	{
		DataPtrsForWfrEdgeCorr.dx = pRadAccessData->xStep;
		DataPtrsForWfrEdgeCorr.dz = pRadAccessData->zStep;

		long TwoNx = pRadAccessData->nx << 1;
		long TwoNz = pRadAccessData->nz << 1;

		if(dxSt != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrXSt = new float[TwoNx];
			if(DataPtrsForWfrEdgeCorr.ExpArrXSt == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrXStEx = new float[TwoNz << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrXStEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrXStEz = DataPtrsForWfrEdgeCorr.FFTArrXStEx + TwoNz;
			DataPtrsForWfrEdgeCorr.dxSt = dxSt;

			long jxSt = ixWfrMinLower + 1;
			double xjSt = pRadAccessData->xStart + jxSt*pRadAccessData->xStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrXSt, pRadAccessData->nx, xjSt, FFT2DInfo.xStartTr, FFT2DInfo.xStepTr);

			SetupRadXorZSectFromSliceConstE(pDataEx, pDataEz, pRadAccessData->nx, pRadAccessData->nz, 'z', jxSt, DataPtrsForWfrEdgeCorr.FFTArrXStEx, DataPtrsForWfrEdgeCorr.FFTArrXStEz);

			if(dzSt != 0.)
			{
				long jzSt2 = (izWfrMinLower + 1) << 1;
				DataPtrsForWfrEdgeCorr.fxStzSt[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzSt2);
				DataPtrsForWfrEdgeCorr.fxStzSt[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzSt2 + 1);
				DataPtrsForWfrEdgeCorr.fxStzSt[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzSt2);
				DataPtrsForWfrEdgeCorr.fxStzSt[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzSt2 + 1);
			}
			if(dzFi != 0.)
			{
				long jzFi2 = izWfrMaxLower << 1;
				DataPtrsForWfrEdgeCorr.fxStzFi[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzFi2);
				DataPtrsForWfrEdgeCorr.fxStzFi[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEx + jzFi2 + 1);
				DataPtrsForWfrEdgeCorr.fxStzFi[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzFi2);
				DataPtrsForWfrEdgeCorr.fxStzFi[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXStEz + jzFi2 + 1);
			}

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrXStEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = pRadAccessData->zStep;
			FFT1DInfo.xStart = pRadAccessData->zStart;
			FFT1DInfo.Nx = pRadAccessData->nz;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dxFi != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrXFi = new float[TwoNx];
			if(DataPtrsForWfrEdgeCorr.ExpArrXFi == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrXFiEx = new float[TwoNz << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrXFiEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrXFiEz = DataPtrsForWfrEdgeCorr.FFTArrXFiEx + TwoNz;
			DataPtrsForWfrEdgeCorr.dxFi = dxFi;

			double xjFi = pRadAccessData->xStart + ixWfrMaxLower*pRadAccessData->xStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrXFi, pRadAccessData->nx, xjFi, FFT2DInfo.xStartTr, FFT2DInfo.xStepTr);

			SetupRadXorZSectFromSliceConstE(pDataEx, pDataEz, pRadAccessData->nx, pRadAccessData->nz, 'z', ixWfrMaxLower, DataPtrsForWfrEdgeCorr.FFTArrXFiEx, DataPtrsForWfrEdgeCorr.FFTArrXFiEz);

			if(dzSt != 0.)
			{
				long jzSt2 = (izWfrMinLower + 1) << 1;
				DataPtrsForWfrEdgeCorr.fxFizSt[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzSt2);
				DataPtrsForWfrEdgeCorr.fxFizSt[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzSt2 + 1);
				DataPtrsForWfrEdgeCorr.fxFizSt[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzSt2);
				DataPtrsForWfrEdgeCorr.fxFizSt[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzSt2 + 1);
			}
			if(dzFi != 0.)
			{
				long jzFi2 = izWfrMaxLower << 1;
				DataPtrsForWfrEdgeCorr.fxFizFi[0] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzFi2);
				DataPtrsForWfrEdgeCorr.fxFizFi[1] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEx + jzFi2 + 1);
				DataPtrsForWfrEdgeCorr.fxFizFi[2] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzFi2);
				DataPtrsForWfrEdgeCorr.fxFizFi[3] = *(DataPtrsForWfrEdgeCorr.FFTArrXFiEz + jzFi2 + 1);
			}

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrXFiEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = pRadAccessData->zStep;
			FFT1DInfo.xStart = pRadAccessData->zStart;
			FFT1DInfo.Nx = pRadAccessData->nz;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dzSt != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrZSt = new float[TwoNz];
			if(DataPtrsForWfrEdgeCorr.ExpArrZSt == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrZStEx = new float[TwoNx << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrZStEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrZStEz = DataPtrsForWfrEdgeCorr.FFTArrZStEx + TwoNx;
			DataPtrsForWfrEdgeCorr.dzSt = dzSt;

			long jzSt = izWfrMinLower + 1;
			double zjSt = pRadAccessData->zStart + jzSt*pRadAccessData->zStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrZSt, pRadAccessData->nz, zjSt, FFT2DInfo.yStartTr, FFT2DInfo.yStepTr);

			SetupRadXorZSectFromSliceConstE(pDataEx, pDataEz, pRadAccessData->nx, pRadAccessData->nz, 'x', jzSt, DataPtrsForWfrEdgeCorr.FFTArrZStEx, DataPtrsForWfrEdgeCorr.FFTArrZStEz);

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrZStEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = pRadAccessData->xStep;
			FFT1DInfo.xStart = pRadAccessData->xStart;
			FFT1DInfo.Nx = pRadAccessData->nx;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		if(dzFi != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrZFi = new float[TwoNz];
			if(DataPtrsForWfrEdgeCorr.ExpArrZFi == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.FFTArrZFiEx = new float[TwoNx << 1];
			if(DataPtrsForWfrEdgeCorr.FFTArrZFiEx == 0) return MEMORY_ALLOCATION_FAILURE;
			DataPtrsForWfrEdgeCorr.FFTArrZFiEz = DataPtrsForWfrEdgeCorr.FFTArrZFiEx + TwoNx;
			DataPtrsForWfrEdgeCorr.dzFi = dzFi;

			double zjFi = pRadAccessData->zStart + izWfrMaxLower*pRadAccessData->zStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrZFi, pRadAccessData->nz, zjFi, FFT2DInfo.yStartTr, FFT2DInfo.yStepTr);

			SetupRadXorZSectFromSliceConstE(pDataEx, pDataEz, pRadAccessData->nx, pRadAccessData->nz, 'x', izWfrMaxLower, DataPtrsForWfrEdgeCorr.FFTArrZFiEx, DataPtrsForWfrEdgeCorr.FFTArrZFiEz);

			FFT1DInfo.pInData = DataPtrsForWfrEdgeCorr.FFTArrZFiEx;
			FFT1DInfo.pOutData = 0;
			FFT1DInfo.xStep = pRadAccessData->xStep;
			FFT1DInfo.xStart = pRadAccessData->xStart;
			FFT1DInfo.Nx = pRadAccessData->nx;
			CGenMathFFT1D FFT1D;
			if(result = FFT1D.Make1DFFT_InPlace(FFT1DInfo)) return result;
		}
		DataPtrsForWfrEdgeCorr.WasSetup = 1;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::SetupWfrEdgeCorrData1D(srTRadSect1D* pRadSect1D, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr1D& DataPtrsForWfrEdgeCorr)
{
	//int result;
	double AbsTol = 0.05*pRadSect1D->ArgStep;

	double WfrMinOffsetFromStart = pRadSect1D->WfrMin - pRadSect1D->ArgStart;
	long iWfrMinLower = long(WfrMinOffsetFromStart/pRadSect1D->ArgStep + 1.E-13);
	double WfrMinLowerMisfit = WfrMinOffsetFromStart - iWfrMinLower*pRadSect1D->ArgStep;

	double WfrMaxOffsetFromStart = pRadSect1D->WfrMax - pRadSect1D->ArgStart;
	long iWfrMaxLower = long(WfrMaxOffsetFromStart/pRadSect1D->ArgStep + 1.E-13);
	double WfrMaxLowerMisfit = WfrMaxOffsetFromStart - iWfrMaxLower*pRadSect1D->ArgStep;

	char WfrMinIsBetweenMeshPoints = (WfrMinLowerMisfit > AbsTol);
	char WfrMaxIsBetweenMeshPoints = (WfrMaxLowerMisfit > AbsTol);
	char WfrMaxIsSmallerThanDataEnd = (::fabs((pRadSect1D->ArgStart + pRadSect1D->np*pRadSect1D->ArgStep) - pRadSect1D->WfrMax) > AbsTol);
	char WfrCorrNeeded = (WfrMinIsBetweenMeshPoints || WfrMaxIsBetweenMeshPoints || WfrMaxIsSmallerThanDataEnd);

	float dSt = 0.;
	if(WfrMinIsBetweenMeshPoints) dSt = (float)(pRadSect1D->ArgStep - WfrMinLowerMisfit);
	float dFi = 0.;
	if(WfrMaxIsBetweenMeshPoints) dFi = (float)(pRadSect1D->ArgStep - WfrMaxLowerMisfit);
	else if(WfrMaxIsSmallerThanDataEnd) dFi = (float)(pRadSect1D->ArgStep);

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = pRadSect1D->ArgStep;
	FFT1DInfo.xStart = pRadSect1D->ArgStart;
	FFT1DInfo.Nx = pRadSect1D->np;
	FFT1DInfo.UseGivenStartTrValue = 0;
	CGenMathFFT1D FFT1D;
	FFT1D.SetupLimitsTr(FFT1DInfo);

	if(WfrCorrNeeded)
	{
		DataPtrsForWfrEdgeCorr.d = pRadSect1D->ArgStep;

		long TwoN = pRadSect1D->np << 1;

		if(dSt != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrSt = new float[TwoN];
			if(DataPtrsForWfrEdgeCorr.ExpArrSt == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.dSt = dSt;

			long jSt = iWfrMinLower + 1;
			double ArgjSt = pRadSect1D->ArgStart + jSt*pRadSect1D->ArgStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrSt, pRadSect1D->np, ArgjSt, FFT1DInfo.xStartTr, FFT1DInfo.xStepTr);
		}
		if(dFi != 0.)
		{
			DataPtrsForWfrEdgeCorr.ExpArrFi = new float[TwoN];
			if(DataPtrsForWfrEdgeCorr.ExpArrFi == 0) return MEMORY_ALLOCATION_FAILURE;

			DataPtrsForWfrEdgeCorr.dFi = dFi;

			double ArgjFi = pRadSect1D->ArgStart + iWfrMaxLower*pRadSect1D->ArgStep;
			SetupExpCorrArray(DataPtrsForWfrEdgeCorr.ExpArrFi, pRadSect1D->np, ArgjFi, FFT1DInfo.xStartTr, FFT1DInfo.xStepTr);
		}
		DataPtrsForWfrEdgeCorr.WasSetup = 1;
	}
	return 0;
}

//*************************************************************************

void srTGenOptElem::MakeWfrEdgeCorrection(srTSRWRadStructAccessData* pRadAccessData, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrs)
{
	float *tEx = pDataEx, *tEz = pDataEz;

	double dxSt_dzSt = DataPtrs.dxSt*DataPtrs.dzSt;
	double dxSt_dzFi = DataPtrs.dxSt*DataPtrs.dzFi;
	double dxFi_dzSt = DataPtrs.dxFi*DataPtrs.dzSt;
	double dxFi_dzFi = DataPtrs.dxFi*DataPtrs.dzFi;

	float fSSExRe = DataPtrs.fxStzSt[0], fSSExIm = DataPtrs.fxStzSt[1], fSSEzRe = DataPtrs.fxStzSt[2], fSSEzIm = DataPtrs.fxStzSt[3];
	float fFSExRe = DataPtrs.fxFizSt[0], fFSExIm = DataPtrs.fxFizSt[1], fFSEzRe = DataPtrs.fxFizSt[2], fFSEzIm = DataPtrs.fxFizSt[3];
	float fSFExRe = DataPtrs.fxStzFi[0], fSFExIm = DataPtrs.fxStzFi[1], fSFEzRe = DataPtrs.fxStzFi[2], fSFEzIm = DataPtrs.fxStzFi[3];
	float fFFExRe = DataPtrs.fxFizFi[0], fFFExIm = DataPtrs.fxFizFi[1], fFFEzRe = DataPtrs.fxFizFi[2], fFFEzIm = DataPtrs.fxFizFi[3];

	float bRe, bIm, cRe, cIm;

	for(long iz=0; iz<pRadAccessData->nz; iz++)
	{
		//long Two_iz = iz << 1;
		//long Two_iz_p_1 = Two_iz + 1;
		long long Two_iz = iz << 1;
		long long Two_iz_p_1 = Two_iz + 1;

		for(long ix=0; ix<pRadAccessData->nx; ix++)
		{
			//long Two_ix = ix << 1;
			//long Two_ix_p_1 = Two_ix + 1;
			long long Two_ix = ix << 1;
			long long Two_ix_p_1 = Two_ix + 1;

			float ExRe = *tEx, ExIm = *(tEx+1);
			float EzRe = *tEz, EzIm = *(tEz+1);

			if(DataPtrs.dxSt != 0.)
			{
				float ExpXStRe = DataPtrs.ExpArrXSt[Two_ix], ExpXStIm = DataPtrs.ExpArrXSt[Two_ix_p_1];

				bRe = DataPtrs.FFTArrXStEx[Two_iz]; bIm = DataPtrs.FFTArrXStEx[Two_iz_p_1];
				ExRe += (float)(DataPtrs.dxSt*(ExpXStRe*bRe - ExpXStIm*bIm));
				ExIm += (float)(DataPtrs.dxSt*(ExpXStRe*bIm + ExpXStIm*bRe));

				bRe = DataPtrs.FFTArrXStEz[Two_iz]; bIm = DataPtrs.FFTArrXStEz[Two_iz_p_1];
				EzRe += (float)(DataPtrs.dxSt*(ExpXStRe*bRe - ExpXStIm*bIm));
				EzIm += (float)(DataPtrs.dxSt*(ExpXStRe*bIm + ExpXStIm*bRe));

				if(DataPtrs.dzSt != 0.)
				{
					bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
					cRe = ExpXStRe*bRe - ExpXStIm*bIm; cIm = ExpXStRe*bIm + ExpXStIm*bRe;

					ExRe += (float)(dxSt_dzSt*(fSSExRe*cRe - fSSExIm*cIm));
					ExIm += (float)(dxSt_dzSt*(fSSExRe*cIm + fSSExIm*cRe));
					EzRe += (float)(dxSt_dzSt*(fSSEzRe*cRe - fSSEzIm*cIm));
					EzIm += (float)(dxSt_dzSt*(fSSEzRe*cIm + fSSEzIm*cRe));
				}
				if(DataPtrs.dzFi != 0.)
				{
					bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
					cRe = ExpXStRe*bRe - ExpXStIm*bIm; cIm = ExpXStRe*bIm + ExpXStIm*bRe;

					ExRe -= (float)(dxSt_dzFi*(fSFExRe*cRe - fSFExIm*cIm));
					ExIm -= (float)(dxSt_dzFi*(fSFExRe*cIm + fSFExIm*cRe));
					EzRe -= (float)(dxSt_dzFi*(fSFEzRe*cRe - fSFEzIm*cIm));
					EzIm -= (float)(dxSt_dzFi*(fSFEzRe*cIm + fSFEzIm*cRe));
				}
			}
			if(DataPtrs.dxFi != 0.)
			{
				float ExpXFiRe = DataPtrs.ExpArrXFi[Two_ix], ExpXFiIm = DataPtrs.ExpArrXFi[Two_ix_p_1];

				bRe = DataPtrs.FFTArrXFiEx[Two_iz]; bIm = DataPtrs.FFTArrXFiEx[Two_iz_p_1];
				ExRe -= (float)(DataPtrs.dxFi*(ExpXFiRe*bRe - ExpXFiIm*bIm));
				ExIm -= (float)(DataPtrs.dxFi*(ExpXFiRe*bIm + ExpXFiIm*bRe));

				bRe = DataPtrs.FFTArrXFiEz[Two_iz]; bIm = DataPtrs.FFTArrXFiEz[Two_iz_p_1];
				EzRe -= (float)(DataPtrs.dxFi*(ExpXFiRe*bRe - ExpXFiIm*bIm));
				EzIm -= (float)(DataPtrs.dxFi*(ExpXFiRe*bIm + ExpXFiIm*bRe));

				if(DataPtrs.dzSt != 0.)
				{
					bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
					cRe = ExpXFiRe*bRe - ExpXFiIm*bIm; cIm = ExpXFiRe*bIm + ExpXFiIm*bRe;

					ExRe -= (float)(dxFi_dzSt*(fFSExRe*cRe - fFSExIm*cIm));
					ExIm -= (float)(dxFi_dzSt*(fFSExRe*cIm + fFSExIm*cRe));
					EzRe -= (float)(dxFi_dzSt*(fFSEzRe*cRe - fFSEzIm*cIm));
					EzIm -= (float)(dxFi_dzSt*(fFSEzRe*cIm + fFSEzIm*cRe));
				}
				if(DataPtrs.dzFi != 0.)
				{
					bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
					cRe = ExpXFiRe*bRe - ExpXFiIm*bIm; cIm = ExpXFiRe*bIm + ExpXFiIm*bRe;

					ExRe += (float)(dxFi_dzFi*(fFFExRe*cRe - fFFExIm*cIm));
					ExIm += (float)(dxFi_dzFi*(fFFExRe*cIm + fFFExIm*cRe));
					EzRe += (float)(dxFi_dzFi*(fFFEzRe*cRe - fFFEzIm*cIm));
					EzIm += (float)(dxFi_dzFi*(fFFEzRe*cIm + fFFEzIm*cRe));
				}
			}
			if(DataPtrs.dzSt != 0.)
			{
				float ExpZStRe = DataPtrs.ExpArrZSt[Two_iz], ExpZStIm = DataPtrs.ExpArrZSt[Two_iz_p_1];

				bRe = DataPtrs.FFTArrZStEx[Two_ix]; bIm = DataPtrs.FFTArrZStEx[Two_ix_p_1];
				ExRe += (float)(DataPtrs.dzSt*(ExpZStRe*bRe - ExpZStIm*bIm));
				ExIm += (float)(DataPtrs.dzSt*(ExpZStRe*bIm + ExpZStIm*bRe));

				bRe = DataPtrs.FFTArrZStEz[Two_ix]; bIm = DataPtrs.FFTArrZStEz[Two_ix_p_1];
				EzRe += (float)(DataPtrs.dzSt*(ExpZStRe*bRe - ExpZStIm*bIm));
				EzIm += (float)(DataPtrs.dzSt*(ExpZStRe*bIm + ExpZStIm*bRe));
			}
			if(DataPtrs.dzFi != 0.)
			{
				float ExpZFiRe = DataPtrs.ExpArrZFi[Two_iz], ExpZFiIm = DataPtrs.ExpArrZFi[Two_iz_p_1];

				bRe = DataPtrs.FFTArrZFiEx[Two_ix]; bIm = DataPtrs.FFTArrZFiEx[Two_ix_p_1];
				ExRe -= (float)(DataPtrs.dzFi*(ExpZFiRe*bRe - ExpZFiIm*bIm));
				ExIm -= (float)(DataPtrs.dzFi*(ExpZFiRe*bIm + ExpZFiIm*bRe));

				bRe = DataPtrs.FFTArrZFiEz[Two_ix]; bIm = DataPtrs.FFTArrZFiEz[Two_ix_p_1];
				EzRe -= (float)(DataPtrs.dzFi*(ExpZFiRe*bRe - ExpZFiIm*bIm));
				EzIm -= (float)(DataPtrs.dzFi*(ExpZFiRe*bIm + ExpZFiIm*bRe));
			}

			*tEx = ExRe; *(tEx+1) = ExIm;
			*tEz = EzRe; *(tEz+1) = EzIm;

			tEx += 2; tEz += 2;
		}
	}
}

//*************************************************************************

void srTGenOptElem::MakeWfrEdgeCorrection1D(srTRadSect1D* pRadSect1D, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr1D& DataPtrs)
{
	float *tEx = pDataEx, *tEz = pDataEz;

	float fSExRe = DataPtrs.fSt[0], fSExIm = DataPtrs.fSt[1], fSEzRe = DataPtrs.fSt[2], fSEzIm = DataPtrs.fSt[3];
	float fFExRe = DataPtrs.fFi[0], fFExIm = DataPtrs.fFi[1], fFEzRe = DataPtrs.fFi[2], fFEzIm = DataPtrs.fFi[3];

	float bRe, bIm;

	//for(long i=0; i<pRadSect1D->np; i++)
	for(long long i=0; i<pRadSect1D->np; i++)
	{
		//long Two_i = i << 1;
		//long Two_i_p_1 = Two_i + 1;
		long long Two_i = i << 1;
		long long Two_i_p_1 = Two_i + 1;

		float ExRe = *tEx, ExIm = *(tEx+1);
		float EzRe = *tEz, EzIm = *(tEz+1);

		if(DataPtrs.dSt != 0.)
		{
			float ExpStRe = DataPtrs.ExpArrSt[Two_i], ExpStIm = DataPtrs.ExpArrSt[Two_i_p_1];

			bRe = fSExRe; bIm = fSExIm;
			ExRe += (float)(DataPtrs.dSt*(ExpStRe*bRe - ExpStIm*bIm));
			ExIm += (float)(DataPtrs.dSt*(ExpStRe*bIm + ExpStIm*bRe));

			bRe = fSEzRe; bIm = fSEzIm;
			EzRe += (float)(DataPtrs.dSt*(ExpStRe*bRe - ExpStIm*bIm));
			EzIm += (float)(DataPtrs.dSt*(ExpStRe*bIm + ExpStIm*bRe));
		}
		if(DataPtrs.dFi != 0.)
		{
			float ExpFiRe = DataPtrs.ExpArrFi[Two_i], ExpFiIm = DataPtrs.ExpArrFi[Two_i_p_1];

			bRe = fFExRe; bIm = fFExIm;
			ExRe -= (float)(DataPtrs.dFi*(ExpFiRe*bRe - ExpFiIm*bIm));
			ExIm -= (float)(DataPtrs.dFi*(ExpFiRe*bIm + ExpFiIm*bRe));

			bRe = fFEzRe; bIm = fFEzIm;
			EzRe -= (float)(DataPtrs.dFi*(ExpFiRe*bRe - ExpFiIm*bIm));
			EzIm -= (float)(DataPtrs.dFi*(ExpFiRe*bIm + ExpFiIm*bRe));
		}

		*tEx = ExRe; *(tEx+1) = ExIm;
		*tEz = EzRe; *(tEz+1) = EzIm;

		tEx += 2; tEz += 2;
	}
}

//*************************************************************************

//int srTGenOptElem::SetRadRepres(srTSRWRadStructAccessData* pRadAccessData, char CoordOrAng)
int srTGenOptElem::SetRadRepres(srTSRWRadStructAccessData* pRadAccessData, char CoordOrAng, double* ar_xStartInSlicesE, double* ar_zStartInSlicesE)
{// 0- to coord.; 1- to ang.
	int result;

	char WfrEdgeCorrShouldBeTreated = pRadAccessData->WfrEdgeCorrShouldBeDone; // Turn on/off here

	if(pRadAccessData->Pres == CoordOrAng) return 0;
	char DirFFT = (CoordOrAng == 0)? -1 : 1;

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.Dir = DirFFT;
	FFT2DInfo.UseGivenStartTrValues = 0;
//New
	if((pRadAccessData->AuxLong4 == 7777777) || ((CoordOrAng == 0) && pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord))
	{
		FFT2DInfo.UseGivenStartTrValues = 1;
		FFT2DInfo.xStartTr = pRadAccessData->xStartTr;
		FFT2DInfo.yStartTr = pRadAccessData->zStartTr;
	}
//End New

	CGenMathFFT2D FFT2D;

	if(pRadAccessData->ne == 1)
	{
		srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
		if(WfrEdgeCorrShouldBeTreated)
		{
			if(CoordOrAng == 1)
			{
				if(result = SetupWfrEdgeCorrData(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;
			}
		}

		if(ar_xStartInSlicesE != 0) FFT2DInfo.xStart = *ar_xStartInSlicesE;
		if(ar_zStartInSlicesE != 0) FFT2DInfo.yStart = *ar_zStartInSlicesE;

		FFT2DInfo.pData = pRadAccessData->pBaseRadX;
		if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
		FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
		if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

		if(WfrEdgeCorrShouldBeTreated)
		{
			if(CoordOrAng == 1)
			{
				if(DataPtrsForWfrEdgeCorr.WasSetup)
				{
					MakeWfrEdgeCorrection(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr);
					DataPtrsForWfrEdgeCorr.DisposeData();
				}
			}
		}
	}
	else
	{
		////OC151014 (test)
		//if(!WfrEdgeCorrShouldBeTreated)
		//{
		//	FFT2DInfo.howMany = pRadAccessData->ne;
		//	FFT2DInfo.iStride = pRadAccessData->ne; //periodicity for extracting next element of 2D array
		//	FFT2DInfo.iDist = 1; //periodicity for starting extraction of next 2D array

		//	FFT2DInfo.pData = pRadAccessData->pBaseRadX;
		//	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
		//	FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
		//	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
		//}
		//else
		//{
		//This is the original version; works by "slices"
		//long TwoNxNz = (pRadAccessData->nx*pRadAccessData->nz) << 1;
		long long TwoNxNz = (((long long)(pRadAccessData->nx))*((long long)(pRadAccessData->nz))) << 1;
		float* AuxEx = new float[TwoNxNz];
		if(AuxEx == 0) return MEMORY_ALLOCATION_FAILURE;
		float* AuxEz = new float[TwoNxNz];
		if(AuxEz == 0) return MEMORY_ALLOCATION_FAILURE;

		for(long ie = 0; ie < pRadAccessData->ne; ie++)
		{
			if(result = ExtractRadSliceConstE(pRadAccessData, ie, AuxEx, AuxEz)) return result;

			srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
			if(WfrEdgeCorrShouldBeTreated)
			{
				if(CoordOrAng == 1)
				{
					if(result = SetupWfrEdgeCorrData(pRadAccessData, AuxEx, AuxEz, DataPtrsForWfrEdgeCorr)) return result;
				}
			}

			//After the FFT, all slices will be authomatically brought to the same mesh
			if(ar_xStartInSlicesE != 0) FFT2DInfo.xStart = ar_xStartInSlicesE[ie];
			if(ar_zStartInSlicesE != 0) FFT2DInfo.yStart = ar_zStartInSlicesE[ie];

			FFT2DInfo.pData = AuxEx;
			if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
			FFT2DInfo.pData = AuxEz;
			if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

			if(WfrEdgeCorrShouldBeTreated)
			{
				if(CoordOrAng == 1)
				{
					if(DataPtrsForWfrEdgeCorr.WasSetup)
					{
						MakeWfrEdgeCorrection(pRadAccessData, AuxEx, AuxEz, DataPtrsForWfrEdgeCorr);
						DataPtrsForWfrEdgeCorr.DisposeData();
					}
				}
			}

			if(result = SetupRadSliceConstE(pRadAccessData, ie, AuxEx, AuxEz)) return result;
		}

		if(AuxEx != 0) delete[] AuxEx;
		if(AuxEz != 0) delete[] AuxEz;
		//}
	}

	pRadAccessData->xStep = FFT2DInfo.xStepTr;
	pRadAccessData->zStep = FFT2DInfo.yStepTr;
	pRadAccessData->xStart = FFT2DInfo.xStartTr;
	pRadAccessData->zStart = FFT2DInfo.yStartTr;
	pRadAccessData->Pres = CoordOrAng;

	pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();

	return 0;
}

//*************************************************************************

int srTGenOptElem::SetRadRepres1D(srTRadSect1D* pRadSect1D, char CoordOrAng)
{// 0- to coord.; 1- to ang.
	int result;

	char WfrEdgeCorrShouldBeTreated = pRadSect1D->WfrEdgeCorrShouldBeDone; // Turn on/off here

	if(pRadSect1D->Pres == CoordOrAng) return 0;
	char DirFFT = (CoordOrAng == 0)? -1 : 1;

	float* AuxDataCont = new float[pRadSect1D->np << 2];
	if(AuxDataCont == 0) return MEMORY_ALLOCATION_FAILURE;

	//long Two_np = pRadSect1D->np << 1;
	long long Two_np = pRadSect1D->np << 1;
	float *tEx = pRadSect1D->pEx, *tEz = pRadSect1D->pEz;
	float *tx = AuxDataCont, *tz = AuxDataCont + Two_np;
	for(int i=0; i<Two_np; i++)
	{
		*(tx++) = *(tEx++); *(tz++) = *(tEz++);
	}

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = pRadSect1D->ArgStep;
	FFT1DInfo.xStart = pRadSect1D->ArgStart;
	FFT1DInfo.Nx = pRadSect1D->np;
	FFT1DInfo.Dir = DirFFT;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
//New
	if(pRadSect1D->AuxLong4 == 7777777)
	{
		FFT1DInfo.UseGivenStartTrValue = 1;
		FFT1DInfo.xStartTr = pRadSect1D->ArgStartTr;
	}
//End New
	CGenMathFFT1D FFT1D;

	srTDataPtrsForWfrEdgeCorr1D DataPtrsForWfrEdgeCorr;
	if(WfrEdgeCorrShouldBeTreated)
	{
		if(CoordOrAng == 1)
		{
			if(result = SetupWfrEdgeCorrData1D(pRadSect1D, pRadSect1D->pEx, pRadSect1D->pEz, DataPtrsForWfrEdgeCorr)) return result;
		}
	}

	FFT1DInfo.pInData = AuxDataCont;
	FFT1DInfo.pOutData = pRadSect1D->pEx;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	FFT1DInfo.pInData = AuxDataCont + Two_np;
	FFT1DInfo.pOutData = pRadSect1D->pEz;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	if(WfrEdgeCorrShouldBeTreated)
	{
		if(CoordOrAng == 1)
		{
			if(DataPtrsForWfrEdgeCorr.WasSetup)
			{
				MakeWfrEdgeCorrection1D(pRadSect1D, pRadSect1D->pEx, pRadSect1D->pEz, DataPtrsForWfrEdgeCorr);
				DataPtrsForWfrEdgeCorr.DisposeData();
			}
		}
	}

	pRadSect1D->ArgStep = FFT1DInfo.xStepTr;
	pRadSect1D->ArgStart = FFT1DInfo.xStartTr;
	pRadSect1D->Pres = CoordOrAng;

	pRadSect1D->SetNonZeroWavefrontLimitsToFullRange();

	if(AuxDataCont != 0) delete[] AuxDataCont;
	return 0;
}

//*************************************************************************

int srTGenOptElem::ComputeRadMoments(srTSRWRadStructAccessData* pSRWRadStructAccessData)
{// Here Lengths are in m and Phot energy in eV!
 // This function seems to work correctly only in Frequency-Coordinate representation

	if((pSRWRadStructAccessData->nx <= 1) && (pSRWRadStructAccessData->nz <= 1)) return 0; //OC090112 or return error?

	const double TwoPi = 3.141592653590*2.;
	const double FourPi = TwoPi*2.;
	//const double Inv_eV_In_m = 1.239854E-06;
	const double Inv_eV_In_m = 1.239842E-06;

	srTAuxMatStat AuxMatStat; //OC13112010 (uncommented)
	int IndLims[4];
	const double RelPowForLimits = 0.9; //to steer
	
	CHGenObj hRad(pSRWRadStructAccessData, true); //OC13112010
	int result = 0;

	//if(pSRWRadStructAccessData->Pres != 0)
	//	if(result = SetRadRepres(pSRWRadStructAccessData, 0)) return result;
	bool IsCoordRepres = (pSRWRadStructAccessData->Pres == 0);
	bool IsFreqRepres = (pSRWRadStructAccessData->PresT == 0);

	char WaveFrontTermWasTreated = 0;
	//if(IsFreqRepres && IsCoordRepres && WaveFrontTermCanBeTreated(pSRWRadStructAccessData))
	if(IsFreqRepres && IsCoordRepres && WaveFrontTermCanBeTreated(*pSRWRadStructAccessData))
	{
		pSRWRadStructAccessData->wfrReffX = -1.E+023; //sign to use only RobsX, RobsZ for Quad. Term treatment
		pSRWRadStructAccessData->wfrReffZ = pSRWRadStructAccessData->wfrReffX;

		TreatStronglyOscillatingTerm(*pSRWRadStructAccessData, 'r');
		WaveFrontTermWasTreated = 1;
	}

	double SumsX[22], SumsZ[22], ff[22];

	float *fpX0 = pSRWRadStructAccessData->pBaseRadX;
	float *fpZ0 = pSRWRadStructAccessData->pBaseRadZ;
	bool ExIsOK = fpX0 != 0; //13112011
	bool EzIsOK = fpZ0 != 0;

	//long PerX = pSRWRadStructAccessData->ne << 1;
	//long PerZ = PerX*pSRWRadStructAccessData->nx;
	long long PerX = pSRWRadStructAccessData->ne << 1;
	long long PerZ = PerX*pSRWRadStructAccessData->nx;

	int nx_mi_1 = pSRWRadStructAccessData->nx - 1;
	int nz_mi_1 = pSRWRadStructAccessData->nz - 1;

	double ePh = pSRWRadStructAccessData->eStart; //This assumes wavefront in Time domain; Photon Energy in eV !
	//If wavefront is in Time-domain representation, average photon energy should be used

	//float *fpMomX = pSRWRadStructAccessData->pMomX, *fpMomZ = pSRWRadStructAccessData->pMomZ;
	double *fpMomX = pSRWRadStructAccessData->pMomX, *fpMomZ = pSRWRadStructAccessData->pMomZ; //130311
	int AmOfMom = 11;

	double InvStepX = 1./pSRWRadStructAccessData->xStep;
	double InvStepZ = 1./pSRWRadStructAccessData->zStep;
	double InvStepXe2 = InvStepX*InvStepX, InvStepZe2 = InvStepZ*InvStepZ;

	char ActualScanX = (pSRWRadStructAccessData->nx > 1);
	char ActualScanZ = (pSRWRadStructAccessData->nz > 1);
	char ActualScansXZ = (ActualScanX && ActualScanZ);

	for(int ie=0; ie<pSRWRadStructAccessData->ne; ie++)
	{
		if(!IsFreqRepres)
		{
			ePh = pSRWRadStructAccessData->avgPhotEn; //?? OC041108
		}

		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;
		for(int k=0; k<22; k++) SumsZ[k] = 0.;

		double Lamb_d_FourPi = Inv_eV_In_m/(FourPi*ePh);
		double Lamb_d_TwoPi = 2.*Lamb_d_FourPi;
		double Lamb_d_TwoPiE2 = Lamb_d_TwoPi*Lamb_d_TwoPi;

		double Lamb_m = Lamb_d_FourPi*FourPi;

		double FourPi_d_Lamb = 1./Lamb_d_FourPi;
		//double FourPi_d_Lamb_d_Rx = FourPi_d_Lamb/pSRWRadStructAccessData->RobsX;
		//double FourPi_d_Lamb_d_Rz = FourPi_d_Lamb/pSRWRadStructAccessData->RobsZ;

		double LocRobsX = pSRWRadStructAccessData->RobsX; //OC030409
		if(LocRobsX == 0.) LocRobsX = 100.*Lamb_m;
		double LocRobsZ = pSRWRadStructAccessData->RobsZ;
		if(LocRobsZ == 0.) LocRobsZ = 100.*Lamb_m;

		double FourPi_d_Lamb_d_Rx = FourPi_d_Lamb/LocRobsX;
		double FourPi_d_Lamb_d_Rz = FourPi_d_Lamb/LocRobsZ;

		double FourPi_d_Lamb_d_Rx_xStep = pSRWRadStructAccessData->xStep*FourPi_d_Lamb_d_Rx;
		double FourPi_d_Lamb_d_Rz_zStep = pSRWRadStructAccessData->zStep*FourPi_d_Lamb_d_Rz;
		double TwoPi_d_Lamb_d_Rx_xStep = 0.5*FourPi_d_Lamb_d_Rx_xStep;
		double TwoPi_d_Lamb_d_Rz_zStep = 0.5*FourPi_d_Lamb_d_Rz_zStep;

		double TwoPi_d_Lamb_d_Rx_xStepE2 = TwoPi_d_Lamb_d_Rx_xStep*TwoPi_d_Lamb_d_Rx_xStep;
		double TwoPi_d_Lamb_d_Rz_zStepE2 = TwoPi_d_Lamb_d_Rz_zStep*TwoPi_d_Lamb_d_Rz_zStep;

		srTMomentsPtrs MomXPtrs(fpMomX), MomZPtrs(fpMomZ);

		AuxMatStat.FindIntensityLimitsInds(hRad, ie, RelPowForLimits, IndLims);

		//AuxMatStat.FindIntensityLimitsInds(*pSRWRadStructAccessData, ie, RelPowForLimits, IndLims);
		//not good for computing precisely intensity
		//make decision
		for(int iz=0; iz<pSRWRadStructAccessData->nz; iz++)
		//for(int iz=IndLims[2]; iz<=IndLims[3]; iz++)
		{
			if(result = srYield.Check()) return result;

			bool vertCoordInsidePowLim = ((iz >= IndLims[2]) && (iz <= IndLims[3]));

			for(int k=0; k<22; k++) SumsX[k] = 0.;

			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			float *fpX_StartForX = fpX0 + izPerZ;
			float *fpZ_StartForX = fpZ0 + izPerZ;

			double z = pSRWRadStructAccessData->zStart + iz*pSRWRadStructAccessData->zStep;

			for(int ix=0; ix<pSRWRadStructAccessData->nx; ix++)
			//for(int ix=IndLims[0]; ix<=IndLims[1]; ix++)
			{
				bool horCoordInsidePowLim = ((ix >= IndLims[0]) && (ix <= IndLims[1]));
				bool coordInsidePowLim = vertCoordInsidePowLim && horCoordInsidePowLim;

				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				float *fpX = fpX_StartForX + ixPerX_p_Two_ie;
				float *fpZ = fpZ_StartForX + ixPerX_p_Two_ie;

				double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
				if(ExIsOK)
				{
					ExRe = *fpX;
					ExIm = *(fpX+1);
				}
				if(EzIsOK)
				{
					EzRe = *fpZ;
					EzIm = *(fpZ+1);
				}

				double x = pSRWRadStructAccessData->xStart + ix*pSRWRadStructAccessData->xStep;
				ff[0] = ExRe*ExRe + ExIm*ExIm; // NormX
				ff[11] = EzRe*EzRe + EzIm*EzIm; // NormZ

				ff[1] = x*ff[0]; // <x>
				ff[3] = z*ff[0]; // <z>
				ff[12] = x*ff[11]; // <x>
				ff[14] = z*ff[11]; // <z>

				if(coordInsidePowLim) //OC13112010
				{
					ff[5] = x*ff[1]; // <xx>
					ff[8] = z*ff[3]; // <zz>
					ff[16] = x*ff[12]; // <xx>
					ff[19] = z*ff[14]; // <zz>
				}
				else
				{
					ff[5] = 0.; // <xx>
					ff[8] = 0.; // <zz>
					ff[16] = 0.; // <xx>
					ff[19] = 0.; // <zz>
				}

				if(IsCoordRepres && (ix > 0))
				{
					float *fpX_Prev = fpX - PerX;
					float *fpZ_Prev = fpZ - PerX;

					double ExReM = 0., ExImM = 0., EzReM = 0., EzImM = 0.;
					if(ExIsOK)
					{
						ExReM = *fpX_Prev; ExImM = *(fpX_Prev+1);
					}
					if(EzIsOK)
					{
						EzReM = *fpZ_Prev; EzImM = *(fpZ_Prev+1);
					}

					double ExReP_mi_ExReM = ExRe - ExReM;
					double ExImP_mi_ExImM = ExIm - ExImM;
					double EzReP_mi_EzReM = EzRe - EzReM;
					double EzImP_mi_EzImM = EzIm - EzImM;

					double ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm = ExImP_mi_ExImM*ExRe - ExReP_mi_ExReM*ExIm;
					ff[2] = ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm + TwoPi_d_Lamb_d_Rx_xStep*x*ff[0]; // <x'>

					double EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm = EzImP_mi_EzImM*EzRe - EzReP_mi_EzReM*EzIm;
					ff[13] = EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm + TwoPi_d_Lamb_d_Rx_xStep*x*ff[11]; // <x'>

					if(coordInsidePowLim) //OC13112010
					{
						ff[6] = x*ff[2]; // <xx'>
						ff[7] = (ExReP_mi_ExReM*ExReP_mi_ExReM + ExImP_mi_ExImM*ExImP_mi_ExImM) 
								+ ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm*TwoPi_d_Lamb_d_Rx_xStep*x
								+ TwoPi_d_Lamb_d_Rx_xStepE2*x*x*ff[0]; // <x'x'>
						ff[17] = x*ff[13]; // <xx'>
						ff[18] = EzReP_mi_EzReM*EzReP_mi_EzReM + EzImP_mi_EzImM*EzImP_mi_EzImM
								+ EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm*TwoPi_d_Lamb_d_Rx_xStep*x
								+ TwoPi_d_Lamb_d_Rx_xStepE2*x*x*ff[11]; // <x'x'>
					}
					else
					{
						ff[6] = 0.; // <xx'>
						ff[7] = 0.; // <x'x'>
						ff[17] = 0.; // <xx'>
						ff[18] = 0.; // <x'x'>
					}
				}
				else
				{
					ff[2] = 0.; // <x'>
					ff[6] = 0.; // <xx'>
					ff[7] = 0.; // <x'x'>
					ff[13] = 0.; // <x'>
					ff[17] = 0.; // <xx'>
					ff[18] = 0.; // <x'x'>
				}

				if(IsCoordRepres && (iz > 0))
				{
					float *fpX_Prev = fpX - PerZ;
					float *fpZ_Prev = fpZ - PerZ;

					double ExReM = 0., ExImM = 0, EzReM = 0., EzImM = 0.;
					if(ExIsOK)
					{
						ExReM = *fpX_Prev; ExImM = *(fpX_Prev+1);
					}
					if(EzIsOK)
					{
						EzReM = *fpZ_Prev; EzImM = *(fpZ_Prev+1);
					}

					double ExReP_mi_ExReM = ExRe - ExReM;
					double ExImP_mi_ExImM = ExIm - ExImM;
					double EzReP_mi_EzReM = EzRe - EzReM;
					double EzImP_mi_EzImM = EzIm - EzImM;

					double ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm = ExImP_mi_ExImM*ExRe - ExReP_mi_ExReM*ExIm;
					ff[4] = ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm + TwoPi_d_Lamb_d_Rz_zStep*z*ff[0]; // <z'>

					double EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm = EzImP_mi_EzImM*EzRe - EzReP_mi_EzReM*EzIm;
					ff[15] = EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm + TwoPi_d_Lamb_d_Rz_zStep*z*ff[11]; // <z'>

					if(coordInsidePowLim) //OC13112010
					{
						ff[9] = z*ff[4]; // <zz'>
						ff[10] = ExReP_mi_ExReM*ExReP_mi_ExReM + ExImP_mi_ExImM*ExImP_mi_ExImM
								+ ExImP_mi_ExImM_ExRe_mi_ExReP_mi_ExReM_ExIm*TwoPi_d_Lamb_d_Rz_zStep*z
								+ TwoPi_d_Lamb_d_Rz_zStepE2*z*z*ff[0]; // <z'z'>
						ff[20] = z*ff[15]; // <zz'>
						ff[21] = EzReP_mi_EzReM*EzReP_mi_EzReM + EzImP_mi_EzImM*EzImP_mi_EzImM
								+ EzImP_mi_EzImM_EzRe_mi_EzReP_mi_EzReM_EzIm*TwoPi_d_Lamb_d_Rz_zStep*z
								+ TwoPi_d_Lamb_d_Rz_zStepE2*z*z*ff[11]; // <z'z'>
					}
					else
					{
						ff[9] = 0.; // <zz'>
						ff[10] = 0.; // <z'z'>
						ff[20] = 0.; // <zz'>
						ff[21] = 0.;
					}
				}
				else
				{
					ff[4] = 0.; // <z'>
					ff[9] = 0.; // <zz'>
					ff[10] = 0.; // <z'z'>
					ff[15] = 0.; // <z'>
					ff[20] = 0.; // <zz'>
					ff[21] = 0.; // <z'z'>
				}

				if((ix == 0) || (ix == nx_mi_1)) for(int k=0; k<22; k++) ff[k] *= 0.5;
				if(ix == 1)
				{
					ff[2] *= 0.5; // <x'>
					ff[6] *= 0.5; // <xx'>
					ff[7] *= 0.5; // <x'x'>
					ff[13] *= 0.5; // <x'>
					ff[17] *= 0.5; // <xx'>
					ff[18] *= 0.5; // <x'x'>
				}
				for(int k1=0; k1<22; k1++) SumsX[k1] += ff[k1];
			}
			
			if((iz == 0) || (iz == nz_mi_1)) for(int k2=0; k2<22; k2++) SumsX[k2] *= 0.5;
			if(iz == 1)
			{
				SumsX[4] *= 0.5; // <z'>
				SumsX[9] *= 0.5; // <zz'>
				SumsX[10] *= 0.5; // <z'z'>
				SumsX[15] *= 0.5; // <z'>
				SumsX[20] *= 0.5; // <zz'>
				SumsX[21] *= 0.5; // <z'z'>
			}

			for(int kk=0; kk<22; kk++) SumsZ[kk] += SumsX[kk];
		}

		double xStep_zStep_mm2 = (pSRWRadStructAccessData->xStep)*(pSRWRadStructAccessData->zStep)*1.E+06;

		if(ExIsOK) //OC13112011
		{
			if(SumsZ[0] != 0.)
			{
				double InvNormX = 1./SumsZ[0];
				double InvNormXLamb_d_TwoPi = InvNormX*Lamb_d_TwoPi;
				double InvNormXLamb_d_TwoPiInvStepX = InvNormXLamb_d_TwoPi*InvStepX;
				double InvNormXLamb_d_TwoPiInvStepZ = InvNormXLamb_d_TwoPi*InvStepZ;
				double InvNormXLamb_d_TwoPiE2 = InvNormX*Lamb_d_TwoPiE2;
				double InvNormXLamb_d_TwoPiE2InvStepXe2 = InvNormXLamb_d_TwoPiE2*InvStepXe2;
				double InvNormXLamb_d_TwoPiE2InvStepZe2 = InvNormXLamb_d_TwoPiE2*InvStepZe2;

				if(IsCoordRepres)
				{
					*(MomXPtrs.pTotPhot) = (float)(ActualScansXZ? SumsZ[0]*xStep_zStep_mm2 : 0.); // Since we have Length in m here and Phot./s/mm^2/(0.1%BW) for Intensity
					*(MomXPtrs.pX) = (float)(ActualScanX? SumsZ[1]*InvNormX : pSRWRadStructAccessData->xStart); // <x>

					*(MomXPtrs.pXP) = (float)(ActualScanX? SumsZ[2]*InvNormXLamb_d_TwoPiInvStepX : 0.); // <x'>
					//*(MomXPtrs.pXP) = (float)(ActualScanX? -SumsZ[2]*InvNormXLamb_d_TwoPiInvStepX : 0.); // <x'> //OC041206

					*(MomXPtrs.pZ) = (float)(ActualScanZ? SumsZ[3]*InvNormX : pSRWRadStructAccessData->zStart); // <z>

					*(MomXPtrs.pZP) = (float)(ActualScanZ? SumsZ[4]*InvNormXLamb_d_TwoPiInvStepZ : 0.); // <z'>
					//*(MomXPtrs.pZP) = (float)(ActualScanZ? -SumsZ[4]*InvNormXLamb_d_TwoPiInvStepZ : 0.); // <z'> //OC041206

					*(MomXPtrs.pXX) = (float)(ActualScanX? SumsZ[5]*InvNormX : (pSRWRadStructAccessData->xStart)*(pSRWRadStructAccessData->xStart)); // <xx>

					*(MomXPtrs.pXXP) = (float)(ActualScanX? SumsZ[6]*InvNormXLamb_d_TwoPiInvStepX : 0.); // <xx'>
					//*(MomXPtrs.pXXP) = (float)(ActualScanX? -SumsZ[6]*InvNormXLamb_d_TwoPiInvStepX : 0.); // <xx'> //OC041206

					*(MomXPtrs.pXPXP) = (float)(ActualScanX? SumsZ[7]*InvNormXLamb_d_TwoPiE2InvStepXe2 : 0.); // <x'x'>
					*(MomXPtrs.pZZ) = (float)(ActualScanZ? SumsZ[8]*InvNormX : (pSRWRadStructAccessData->zStart)*(pSRWRadStructAccessData->zStart)); // <zz>

					*(MomXPtrs.pZZP) = (float)(ActualScanZ? SumsZ[9]*InvNormXLamb_d_TwoPiInvStepZ : 0.); // <zz'>
					//*(MomXPtrs.pZZP) = (float)(ActualScanZ? -SumsZ[9]*InvNormXLamb_d_TwoPiInvStepZ : 0.); // <zz'> //OC041206

					*(MomXPtrs.pZPZP) = (float)(ActualScanZ? SumsZ[10]*InvNormXLamb_d_TwoPiE2InvStepZe2 : 0.); // <z'z'>
				}
				else
				{
					*(MomXPtrs.pXP) = (float)(Lamb_m*(ActualScanX? SumsZ[1]*InvNormX : pSRWRadStructAccessData->xStart)); // <x'>
					*(MomXPtrs.pZP) = (float)(Lamb_m*(ActualScanZ? SumsZ[3]*InvNormX : pSRWRadStructAccessData->zStart)); // <z'>
					*(MomXPtrs.pXPXP) = (float)(Lamb_m*Lamb_m*(ActualScanX? SumsZ[5]*InvNormX : (pSRWRadStructAccessData->xStart)*(pSRWRadStructAccessData->xStart))); // <x'x'>
					*(MomXPtrs.pZPZP) = (float)(Lamb_m*Lamb_m*(ActualScanZ? SumsZ[8]*InvNormX : (pSRWRadStructAccessData->zStart)*(pSRWRadStructAccessData->zStart))); // <z'z'>
				}
			}
			else
			{
				*(MomXPtrs.pTotPhot) = 0.;
				*(MomXPtrs.pX) = 0.; // <x>
				*(MomXPtrs.pXP) = 0.; // <x'>
				*(MomXPtrs.pZ) = 0.; // <z>
				*(MomXPtrs.pZP) = 0.; // <z'>
				*(MomXPtrs.pXX) = 0.; // <xx>
				*(MomXPtrs.pXXP) = 0.; // <xx'>
				*(MomXPtrs.pXPXP) = 0.; // <x'x'>
				*(MomXPtrs.pZZ) = 0.; // <zz>
				*(MomXPtrs.pZZP) = 0.; // <zz'>
				*(MomXPtrs.pZPZP) = 0.; // <z'z'>
			}
		}

		if(EzIsOK) //OC13112011
		{
			if(SumsZ[11] != 0.)
			{
				double InvNormZ = 1./SumsZ[11];
				double InvNormZLamb_d_TwoPi = InvNormZ*Lamb_d_TwoPi;
				double InvNormZLamb_d_TwoPiInvStepX = InvNormZLamb_d_TwoPi*InvStepX;
				double InvNormZLamb_d_TwoPiInvStepZ = InvNormZLamb_d_TwoPi*InvStepZ;
				double InvNormZLamb_d_TwoPiE2 = InvNormZ*Lamb_d_TwoPiE2;
				double InvNormZLamb_d_TwoPiE2InvStepXe2 = InvNormZLamb_d_TwoPiE2*InvStepXe2;
				double InvNormZLamb_d_TwoPiE2InvStepZe2 = InvNormZLamb_d_TwoPiE2*InvStepZe2;

				if(IsCoordRepres)
				{
					*(MomZPtrs.pTotPhot) = (float)(ActualScansXZ? SumsZ[11]*xStep_zStep_mm2 : 0.); // Since we have Length in m here and Phot./s/mm^2/(0.1%BW) for Intensity
					*(MomZPtrs.pX) = (float)(ActualScanX? SumsZ[12]*InvNormZ : pSRWRadStructAccessData->xStart); // <x>

					*(MomZPtrs.pXP) = (float)(ActualScanX? SumsZ[13]*InvNormZLamb_d_TwoPiInvStepX : 0.); // <x'>
					//*(MomZPtrs.pXP) = (float)(ActualScanX? -SumsZ[13]*InvNormZLamb_d_TwoPiInvStepX : 0.); // <x'> //OC041206

					*(MomZPtrs.pZ) = (float)(ActualScanZ? SumsZ[14]*InvNormZ : pSRWRadStructAccessData->zStart); // <z>

					*(MomZPtrs.pZP) = (float)(ActualScanZ? SumsZ[15]*InvNormZLamb_d_TwoPiInvStepZ : 0.); // <z'>
					//*(MomZPtrs.pZP) = (float)(ActualScanZ? -SumsZ[15]*InvNormZLamb_d_TwoPiInvStepZ : 0.); // <z'> //OC041206

					*(MomZPtrs.pXX) = (float)(ActualScanX? SumsZ[16]*InvNormZ : (pSRWRadStructAccessData->xStart)*(pSRWRadStructAccessData->xStart)); // <xx>

					*(MomZPtrs.pXXP) = (float)(ActualScanX? SumsZ[17]*InvNormZLamb_d_TwoPiInvStepX : 0.); // <xx'>
					//*(MomZPtrs.pXXP) = (float)(ActualScanX? -SumsZ[17]*InvNormZLamb_d_TwoPiInvStepX : 0.); // <xx'> //OC041206

					*(MomZPtrs.pXPXP) = (float)(ActualScanX? SumsZ[18]*InvNormZLamb_d_TwoPiE2InvStepXe2 : 0.); // <x'x'>
					*(MomZPtrs.pZZ) = (float)(ActualScanZ? SumsZ[19]*InvNormZ : (pSRWRadStructAccessData->zStart)*(pSRWRadStructAccessData->zStart)); // <zz>

					*(MomZPtrs.pZZP) = (float)(ActualScanZ? SumsZ[20]*InvNormZLamb_d_TwoPiInvStepZ : 0.); // <zz'>
					//*(MomZPtrs.pZZP) = (float)(ActualScanZ? -SumsZ[20]*InvNormZLamb_d_TwoPiInvStepZ : 0.); // <zz'> //OC041206

					*(MomZPtrs.pZPZP) = (float)(ActualScanZ? SumsZ[21]*InvNormZLamb_d_TwoPiE2InvStepZe2 : 0.); // <z'z'>
				}
				else
				{
					*(MomZPtrs.pXP) = (float)(Lamb_m*(ActualScanX? SumsZ[12]*InvNormZ : pSRWRadStructAccessData->xStart)); // <x'>
					*(MomZPtrs.pZP) = (float)(Lamb_m*(ActualScanZ? SumsZ[14]*InvNormZ : pSRWRadStructAccessData->zStart)); // <z'>
					*(MomZPtrs.pXPXP) = (float)(Lamb_m*Lamb_m*(ActualScanX? SumsZ[16]*InvNormZ : (pSRWRadStructAccessData->xStart)*(pSRWRadStructAccessData->xStart))); // <x'x'>
					*(MomZPtrs.pZPZP) = (float)(Lamb_m*Lamb_m*(ActualScanZ? SumsZ[19]*InvNormZ : (pSRWRadStructAccessData->zStart)*(pSRWRadStructAccessData->zStart))); // <z'z'>
				}
			}
			else
			{
				*(MomZPtrs.pTotPhot) = 0.; // Since we have Length in m here and Phot./s/mm^2/(0.1%BW) for Intensity
				*(MomZPtrs.pX) = 0.; // <x>
				*(MomZPtrs.pXP) = 0.; // <x'>
				*(MomZPtrs.pZ) = 0.; // <z>
				*(MomZPtrs.pZP) = 0.; // <z'>
				*(MomZPtrs.pXX) = 0.; // <xx>
				*(MomZPtrs.pXXP) = 0.; // <xx'>
				*(MomZPtrs.pXPXP) = 0.; // <x'x'>
				*(MomZPtrs.pZZ) = 0.; // <zz>
				*(MomZPtrs.pZZP) = 0.; // <zz'>
				*(MomZPtrs.pZPZP) = 0.; // <z'z'>
			}
		}
		
		ePh += pSRWRadStructAccessData->eStep;
		fpMomX += AmOfMom; fpMomZ += AmOfMom;
	}

	if(WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(*pSRWRadStructAccessData, 'a');
	pSRWRadStructAccessData->MomWereCalcNum = true;
	return 0;
}

//*************************************************************************

//int srTGenOptElem::GenAuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, float** ax, float** az, srTMomentsRatios* MomRatArray)
int srTGenOptElem::GenAuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, double** ax, double** az, srTMomentsRatios* MomRatArray) //OC130311
{// Drift Space has its own realization of this function
	//float bxStr0[] = { (float)(ax[0][0]*ax[0][0]), (float)(2.*ax[0][0]*ax[0][1]), (float)(ax[0][1]*ax[0][1]) };
	//float bxStr1[] = { (float)(ax[0][0]*ax[1][0]), (float)(ax[0][1]*ax[1][0] + ax[0][0]*ax[1][1]), (float)(ax[0][1]*ax[1][1]) };
	//float bxStr2[] = { (float)(ax[1][0]*ax[1][0]), (float)(2.*ax[1][0]*ax[1][1]), (float)(ax[1][1]*ax[1][1]) };
	//float* bx[] = { bxStr0, bxStr1, bxStr2 };

	//float bzStr0[] = { (float)(az[0][0]*az[0][0]), (float)(2.*az[0][0]*az[0][1]), (float)(az[0][1]*az[0][1]) };
	//float bzStr1[] = { (float)(az[0][0]*az[1][0]), (float)(az[0][1]*az[1][0] + az[0][0]*az[1][1]), (float)(az[0][1]*az[1][1]) };
	//float bzStr2[] = { (float)(az[1][0]*az[1][0]), (float)(2.*az[1][0]*az[1][1]), (float)(az[1][1]*az[1][1]) };
	//float* bz[] = { bzStr0, bzStr1, bzStr2 };

	//float v3dAux1[3], v3dAux2[3];
	//float v2dAux1[2], v2dAux2[2];

	//OC130311
	double bxStr0[] = { (ax[0][0]*ax[0][0]), (2.*ax[0][0]*ax[0][1]), (ax[0][1]*ax[0][1]) };
	double bxStr1[] = { (ax[0][0]*ax[1][0]), (ax[0][1]*ax[1][0] + ax[0][0]*ax[1][1]), (ax[0][1]*ax[1][1]) };
	double bxStr2[] = { (ax[1][0]*ax[1][0]), (2.*ax[1][0]*ax[1][1]), (ax[1][1]*ax[1][1]) };
	double* bx[] = { bxStr0, bxStr1, bxStr2 };
	double bzStr0[] = { (az[0][0]*az[0][0]), (2.*az[0][0]*az[0][1]), (az[0][1]*az[0][1]) };
	double bzStr1[] = { (az[0][0]*az[1][0]), (az[0][1]*az[1][0] + az[0][0]*az[1][1]), (az[0][1]*az[1][1]) };
	double bzStr2[] = { (az[1][0]*az[1][0]), (2.*az[1][0]*az[1][1]), (az[1][1]*az[1][1]) };
	double* bz[] = { bzStr0, bzStr1, bzStr2 };
	double v3dAux1[3], v3dAux2[3];
	double v2dAux1[2], v2dAux2[2];

	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	int OffsetE = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + OffsetE);

		//float OldMomX_xx = *(MomX.pXX), OldMomX_xpxp = *(MomX.pXPXP);
		//float OldMomX_zz = *(MomX.pZZ), OldMomX_zpzp = *(MomX.pZPZP);
		double OldMomX_xx = *(MomX.pXX), OldMomX_xpxp = *(MomX.pXPXP); //OC130311
		double OldMomX_zz = *(MomX.pZZ), OldMomX_zpzp = *(MomX.pZPZP);

		*v2dAux1 = *(MomX.pX); *(v2dAux1+1) = *(MomX.pXP);
		MultSquareMatrByVect(ax, v2dAux1, 2, v2dAux2);
		*(MomX.pX) = *v2dAux2; *(MomX.pXP) = *(v2dAux2+1);

		*v2dAux1 = *(MomX.pZ); *(v2dAux1+1) = *(MomX.pZP);
		MultSquareMatrByVect(az, v2dAux1, 2, v2dAux2);
		*(MomX.pZ) = *v2dAux2; *(MomX.pZP) = *(v2dAux2+1);

		*v3dAux1 = *(MomX.pXX); *(v3dAux1+1) = *(MomX.pXXP); *(v3dAux1+2) = *(MomX.pXPXP);
		MultSquareMatrByVect(bx, v3dAux1, 3, v3dAux2);
		*(MomX.pXX) = *v3dAux2; *(MomX.pXXP) = *(v3dAux2+1); *(MomX.pXPXP) = *(v3dAux2+2);

		*v3dAux1 = *(MomX.pZZ); *(v3dAux1+1) = *(MomX.pZZP); *(v3dAux1+2) = *(MomX.pZPZP);
		MultSquareMatrByVect(bz, v3dAux1, 3, v3dAux2);
		*(MomX.pZZ) = *v3dAux2; *(MomX.pZZP) = *(v3dAux2+1); *(MomX.pZPZP) = *(v3dAux2+2);

		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + OffsetE);

		//float OldMomZ_xx = *(MomZ.pXX), OldMomZ_xpxp = *(MomZ.pXPXP);
		//float OldMomZ_zz = *(MomZ.pZZ), OldMomZ_zpzp = *(MomZ.pZPZP);
		double OldMomZ_xx = *(MomZ.pXX), OldMomZ_xpxp = *(MomZ.pXPXP); //OC130311
		double OldMomZ_zz = *(MomZ.pZZ), OldMomZ_zpzp = *(MomZ.pZPZP);

		*v2dAux1 = *(MomZ.pX); *(v2dAux1+1) = *(MomZ.pXP);
		MultSquareMatrByVect(ax, v2dAux1, 2, v2dAux2);
		*(MomZ.pX) = *v2dAux2; *(MomZ.pXP) = *(v2dAux2+1);

		*v2dAux1 = *(MomZ.pZ); *(v2dAux1+1) = *(MomZ.pZP);
		MultSquareMatrByVect(az, v2dAux1, 2, v2dAux2);
		*(MomZ.pZ) = *v2dAux2; *(MomZ.pZP) = *(v2dAux2+1);

		*v3dAux1 = *(MomZ.pXX); *(v3dAux1+1) = *(MomZ.pXXP); *(v3dAux1+2) = *(MomZ.pXPXP);
		MultSquareMatrByVect(bx, v3dAux1, 3, v3dAux2);
		*(MomZ.pXX) = *v3dAux2; *(MomZ.pXXP) = *(v3dAux2+1); *(MomZ.pXPXP) = *(v3dAux2+2);

		*v3dAux1 = *(MomZ.pZZ); *(v3dAux1+1) = *(MomZ.pZZP); *(v3dAux1+2) = *(MomZ.pZPZP);
		MultSquareMatrByVect(bz, v3dAux1, 3, v3dAux2);
		*(MomZ.pZZ) = *v3dAux2; *(MomZ.pZZP) = *(v3dAux2+1); *(MomZ.pZPZP) = *(v3dAux2+2);

		if(FillMomRatios)
		{
			//tMomRatArray->RxxMomX = (float)((*(MomX.pXX) > 0.)? sqrt(*(MomX.pXX)/OldMomX_xx) : -1.);
			//tMomRatArray->RxpxpMomX = (float)((*(MomX.pXPXP) > 0.)? sqrt(*(MomX.pXPXP)/OldMomX_xpxp) : -1.);
			//tMomRatArray->RzzMomX = (float)((*(MomX.pZZ) > 0.)? sqrt(*(MomX.pZZ)/OldMomX_zz) : -1.);
			//tMomRatArray->RzpzpMomX = (float)((*(MomX.pZPZP) > 0.)? sqrt(*(MomX.pZPZP)/OldMomX_zpzp) : -1.);
			//tMomRatArray->RxxMomZ = (float)((*(MomZ.pXX) > 0.)? sqrt(*(MomZ.pXX)/OldMomZ_xx) : -1.);
			//tMomRatArray->RxpxpMomZ = (float)((*(MomZ.pXPXP) > 0.)? sqrt(*(MomZ.pXPXP)/OldMomZ_xpxp) : -1.);
			//tMomRatArray->RzzMomZ = (float)((*(MomZ.pZZ) > 0.)? sqrt(*(MomZ.pZZ)/OldMomZ_zz) : -1.);
			//tMomRatArray->RzpzpMomZ = (float)((*(MomZ.pZPZP) > 0.)? sqrt(*(MomZ.pZPZP)/OldMomZ_zpzp) : -1.);
			//OC130311
			tMomRatArray->RxxMomX = ((*(MomX.pXX) > 0.)? sqrt(*(MomX.pXX)/OldMomX_xx) : -1.);
			tMomRatArray->RxpxpMomX = ((*(MomX.pXPXP) > 0.)? sqrt(*(MomX.pXPXP)/OldMomX_xpxp) : -1.);
			tMomRatArray->RzzMomX = ((*(MomX.pZZ) > 0.)? sqrt(*(MomX.pZZ)/OldMomX_zz) : -1.);
			tMomRatArray->RzpzpMomX = ((*(MomX.pZPZP) > 0.)? sqrt(*(MomX.pZPZP)/OldMomX_zpzp) : -1.);
			tMomRatArray->RxxMomZ = ((*(MomZ.pXX) > 0.)? sqrt(*(MomZ.pXX)/OldMomZ_xx) : -1.);
			tMomRatArray->RxpxpMomZ = ((*(MomZ.pXPXP) > 0.)? sqrt(*(MomZ.pXPXP)/OldMomZ_xpxp) : -1.);
			tMomRatArray->RzzMomZ = ((*(MomZ.pZZ) > 0.)? sqrt(*(MomZ.pZZ)/OldMomZ_zz) : -1.);
			tMomRatArray->RzpzpMomZ = ((*(MomZ.pZPZP) > 0.)? sqrt(*(MomZ.pZPZP)/OldMomZ_zpzp) : -1.);
			tMomRatArray++;
		}

		OffsetE += 11;
	}
	pRadAccessData->MomWereCalcNum = false; //OC13112010

	return 0;
}

//*************************************************************************

//void srTGenOptElem::SetupMxxMzzArr(srTSRWRadStructAccessData* pRadAccessData, float* MxxArr, float* MzzArr)
void srTGenOptElem::SetupMxxMzzArr(srTSRWRadStructAccessData* pRadAccessData, double* MxxArr, double* MzzArr) //OC130311
{
	int OffsetE = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + OffsetE);
		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + OffsetE);

		MxxArr[ie] = (*(MomX.pXX) > *(MomZ.pXX))? *(MomX.pXX) : *(MomZ.pXX);
		MzzArr[ie] = (*(MomX.pZZ) > *(MomZ.pZZ))? *(MomX.pZZ) : *(MomZ.pZZ);

		OffsetE += 11;
	}
}

//*************************************************************************

//void srTGenOptElem::FindMinMaxRatio(float* Arr1, float* Arr2, int n, float& MinRat2to1, float& MaxRat2to1)
void srTGenOptElem::FindMinMaxRatio(double* Arr1, double* Arr2, int n, double& MinRat2to1, double& MaxRat2to1) //OC130311
{
	//float Min = (float)(1.E+23), Max = (float)(1.E-23);
	double Min = (1.E+23), Max = (1.E-23);
	for(int i=0; i<n; i++)
	{
		//float Rat = Arr2[i]/Arr1[i];
		double Rat = Arr2[i]/Arr1[i];
		if(Rat > Max) Max = Rat;
		if(Rat < Min) Min = Rat;
	}
	MinRat2to1 = Min;
	MaxRat2to1 = Max;
}

//*************************************************************************

int srTGenOptElem::RadResizeGen(srTSRWRadStructAccessData& SRWRadStructAccessData, srTRadResize& RadResizeStruct)
{
	if((RadResizeStruct.pxm == 1.) && (RadResizeStruct.pxd == 1.) && (RadResizeStruct.pzm == 1.) && (RadResizeStruct.pzd == 1.)) return 0;
	int result = 0;

	bool ExIsOK = SRWRadStructAccessData.pBaseRadX != 0;
	bool EzIsOK = SRWRadStructAccessData.pBaseRadZ != 0;

//New
	long nxOld = SRWRadStructAccessData.nx;
	long nzOld = SRWRadStructAccessData.nz;
	double xStartOld = SRWRadStructAccessData.xStart, zStartOld = SRWRadStructAccessData.zStart;
	double xStepOld = SRWRadStructAccessData.xStep;
	double zStepOld = SRWRadStructAccessData.zStep;
	double xRangeOld = xStepOld*nxOld;
	double zRangeOld = zStepOld*nzOld;
	char UseStartTrWithOtherSide = (::fabs(xStartOld + 0.5*xRangeOld) > 0.5*SRWRadStructAccessData.xStep) ||
								   (::fabs(zStartOld + 0.5*zRangeOld) > 0.5*SRWRadStructAccessData.zStep);
	double pxmIn = RadResizeStruct.pxm, pxdIn = RadResizeStruct.pxd;
	double pzmIn = RadResizeStruct.pzm, pzdIn = RadResizeStruct.pzd;
//End New

	//if(RadResizeStruct.UseOtherSideFFT)
	if(RadResizeStruct.useOtherSideFFT()) //OC090311
	{
		char ToRepres = (SRWRadStructAccessData.Pres == 0)? 1 : 0;

		if(UseStartTrWithOtherSide)
		{
			SRWRadStructAccessData.xStart = -(SRWRadStructAccessData.nx >> 1)*SRWRadStructAccessData.xStep;
			SRWRadStructAccessData.zStart = -(SRWRadStructAccessData.nz >> 1)*SRWRadStructAccessData.zStep;
			
			double xShift = SRWRadStructAccessData.xStart - xStartOld;
			double zShift = SRWRadStructAccessData.zStart - zStartOld;

			SRWRadStructAccessData.xWfrMin += xShift; SRWRadStructAccessData.xWfrMax += xShift;
			SRWRadStructAccessData.zWfrMin += zShift; SRWRadStructAccessData.zWfrMax += zShift;
		}

		if(result = SetRadRepres(&SRWRadStructAccessData, ToRepres)) return result;

		double pxmNew = RadResizeStruct.pxd, pxdNew = RadResizeStruct.pxm;
		double pzmNew = RadResizeStruct.pzd, pzdNew = RadResizeStruct.pzm;

		RadResizeStruct.pxm = pxmNew; RadResizeStruct.pxd = pxdNew;
		RadResizeStruct.pzm = pzmNew; RadResizeStruct.pzd = pzdNew;
	}

	//srTSRWRadStructAccessData NewSRWRadStructAccessData = SRWRadStructAccessData;
	srTSRWRadStructAccessData NewSRWRadStructAccessData(SRWRadStructAccessData); //OC140411

	double pxTot = RadResizeStruct.pxm*RadResizeStruct.pxd;
	double pzTot = RadResizeStruct.pzm*RadResizeStruct.pzd;

	long NewNx = SRWRadStructAccessData.nx; 
	long NewNz = SRWRadStructAccessData.nz; 

	CGenMathFFT2D FFT;
	char RadShouldBeChanged = ((pxTot != 1.) || (pzTot != 1.));
	if(RadShouldBeChanged)
	{
		double pxTotNewNx = pxTot*NewNx;
		long Long_pxTotNewNx = long(pxTotNewNx);
		NewNx = ((pxTotNewNx - Long_pxTotNewNx) >= 0.5)? Long_pxTotNewNx + 1 : Long_pxTotNewNx;

		double pzTotNewNz = pzTot*NewNz;
		long Long_pzTotNewNz = long(pzTotNewNz);
		NewNz = ((pzTotNewNz - Long_pzTotNewNz) >= 0.5)? Long_pzTotNewNz + 1 : Long_pzTotNewNz;

		FFT.NextCorrectNumberForFFT(NewNx);
		FFT.NextCorrectNumberForFFT(NewNz);

		NewSRWRadStructAccessData.nx = NewNx;
		NewSRWRadStructAccessData.nz = NewNz;
	}

	char CenterIsOffset = 0;
	long iXc = (SRWRadStructAccessData.nx >> 1);
	
	if(::fabs(RadResizeStruct.RelCenPosX - 0.5) > RadResizeStruct.RelCenPosTol)
	{
		double xRangeOld = SRWRadStructAccessData.nx*SRWRadStructAccessData.xStep;
		double xcAppr = SRWRadStructAccessData.xStart + xRangeOld*RadResizeStruct.RelCenPosX;
		iXc = IntegerOffsetCoord(SRWRadStructAccessData.xStart, SRWRadStructAccessData.xStep, xcAppr);
		CenterIsOffset = 1;
	}
	double Xc = SRWRadStructAccessData.xStart + SRWRadStructAccessData.xStep*iXc;

	long iZc = (SRWRadStructAccessData.nz >> 1);
	if(::fabs(RadResizeStruct.RelCenPosZ - 0.5) > RadResizeStruct.RelCenPosTol)
	{
		double zRangeOld = SRWRadStructAccessData.nz*SRWRadStructAccessData.zStep;
		double zcAppr = SRWRadStructAccessData.zStart + zRangeOld*RadResizeStruct.RelCenPosZ;
		iZc = IntegerOffsetCoord(SRWRadStructAccessData.zStart, SRWRadStructAccessData.zStep, zcAppr);
		CenterIsOffset = 1;
	}
	double Zc = SRWRadStructAccessData.zStart + SRWRadStructAccessData.zStep*iZc;

	if(CenterIsOffset)
	{
//x--
		long NewNxLeftHalf = NewNx >> 1;
		NewSRWRadStructAccessData.xStep = SRWRadStructAccessData.xStep/RadResizeStruct.pxd;
		NewSRWRadStructAccessData.xStart = Xc - NewNxLeftHalf*NewSRWRadStructAccessData.xStep;
		double NewNxEndInterp = NewSRWRadStructAccessData.xStart + (NewNx - 1)*NewSRWRadStructAccessData.xStep;
		double OldNxEndInterp = SRWRadStructAccessData.xStart + (SRWRadStructAccessData.nx - 1)*SRWRadStructAccessData.xStep;

		TuneStepToKeepInterpLimitsTheSameAtResize(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'x', iXc);

		if(NewSRWRadStructAccessData.xStart < SRWRadStructAccessData.xStart) // ixStart
			NewSRWRadStructAccessData.AuxLong1 = IntegerOffsetCoord(NewSRWRadStructAccessData.xStart, NewSRWRadStructAccessData.xStep, SRWRadStructAccessData.xStart);
		else NewSRWRadStructAccessData.AuxLong1 = 0;

		if(NewNxEndInterp > OldNxEndInterp) // ixEnd
			NewSRWRadStructAccessData.AuxLong2 = IntegerOffsetCoord(NewSRWRadStructAccessData.xStart, NewSRWRadStructAccessData.xStep, OldNxEndInterp);
		else NewSRWRadStructAccessData.AuxLong2 = NewNx - 1;

		if(RadResizeStruct.pxm > 1.) NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'x');
		else
		{
			if(NewSRWRadStructAccessData.xStart > NewSRWRadStructAccessData.xWfrMin) NewSRWRadStructAccessData.xWfrMin = NewSRWRadStructAccessData.xStart;
			double xEndNew = NewSRWRadStructAccessData.xStart + NewSRWRadStructAccessData.nx*NewSRWRadStructAccessData.xStep;
			if(xEndNew < NewSRWRadStructAccessData.xWfrMax) NewSRWRadStructAccessData.xWfrMax = xEndNew;
		}
//z--
		long NewNzLeftHalf = NewNz >> 1;
		NewSRWRadStructAccessData.zStep = SRWRadStructAccessData.zStep/RadResizeStruct.pzd;
		NewSRWRadStructAccessData.zStart = Zc - NewNzLeftHalf*NewSRWRadStructAccessData.zStep;

		TuneStepToKeepInterpLimitsTheSameAtResize(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'z', iZc);

		if(NewSRWRadStructAccessData.zStart < SRWRadStructAccessData.zStart) // izStart
			NewSRWRadStructAccessData.AuxLong3 = 1 + IntegerOffsetCoord(NewSRWRadStructAccessData.zStart, NewSRWRadStructAccessData.zStep, SRWRadStructAccessData.zStart);
		else NewSRWRadStructAccessData.AuxLong3 = 0;

		double NewNzEndInterp = NewSRWRadStructAccessData.zStart + (NewNz - 1)*NewSRWRadStructAccessData.zStep;
		double OldNzEndInterp = SRWRadStructAccessData.zStart + (SRWRadStructAccessData.nz - 1)*SRWRadStructAccessData.zStep;
		if(NewNzEndInterp > OldNzEndInterp) // izEnd
			NewSRWRadStructAccessData.AuxLong4 = IntegerOffsetCoord(NewSRWRadStructAccessData.zStart, NewSRWRadStructAccessData.zStep, OldNzEndInterp);
		else NewSRWRadStructAccessData.AuxLong4 = NewNz - 1;

		if(RadResizeStruct.pzm > 1.) NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'z');
		else
		{
			if(NewSRWRadStructAccessData.zStart > NewSRWRadStructAccessData.zWfrMin) NewSRWRadStructAccessData.zWfrMin = NewSRWRadStructAccessData.zStart;
			double zEndNew = NewSRWRadStructAccessData.zStart + NewSRWRadStructAccessData.nz*NewSRWRadStructAccessData.zStep;
			if(zEndNew < NewSRWRadStructAccessData.zWfrMax) NewSRWRadStructAccessData.zWfrMax = zEndNew;
		}
	}
	else
	{
//x--
		double xMagOld = (SRWRadStructAccessData.nx - 1)*SRWRadStructAccessData.xStep;
		if(RadResizeStruct.pxm > 1.)
		{
			long NxInner = (RadResizeStruct.pxd != 1.)? long(NewNx/RadResizeStruct.pxm) : SRWRadStructAccessData.nx;

			long NewNx_mi_NxInner = NewNx - NxInner;
			long NxOuterLeft = NewNx_mi_NxInner >> 1;
			if((NxOuterLeft << 1) != NewNx_mi_NxInner) NxOuterLeft++;
		
			NewSRWRadStructAccessData.xStep = xMagOld/(NxInner - 1);
			NewSRWRadStructAccessData.xStart = SRWRadStructAccessData.xStart - NxOuterLeft*NewSRWRadStructAccessData.xStep;
			NewSRWRadStructAccessData.AuxLong1 = NxOuterLeft; // ixStart
			NewSRWRadStructAccessData.AuxLong2 = NxOuterLeft + NxInner - 1; // ixEnd

			NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'x');
		}
		else
		{
			double xEndOld = SRWRadStructAccessData.xStart + (SRWRadStructAccessData.nx - 1)*SRWRadStructAccessData.xStep;
			double xMid = 0.5*(SRWRadStructAccessData.xStart + xEndOld);
			double xMagNew = xMagOld*RadResizeStruct.pxm;
			NewSRWRadStructAccessData.xStep = xMagNew/(NewNx - 1);

			if(RadResizeStruct.pxd == 1)
			{
				NewSRWRadStructAccessData.xStep = SRWRadStructAccessData.xStep;
				xMagNew = (NewNx - 1)*NewSRWRadStructAccessData.xStep;
			}

			NewSRWRadStructAccessData.xStart = xMid - 0.5*xMagNew;
			NewSRWRadStructAccessData.AuxLong1 = 0; // ixStart
			NewSRWRadStructAccessData.AuxLong2 = NewNx - 1; // ixEnd

			if(NewSRWRadStructAccessData.xStart > NewSRWRadStructAccessData.xWfrMin) NewSRWRadStructAccessData.xWfrMin = NewSRWRadStructAccessData.xStart;
			if((NewSRWRadStructAccessData.xStart + xMagNew) < NewSRWRadStructAccessData.xWfrMax) NewSRWRadStructAccessData.xWfrMax = NewSRWRadStructAccessData.xStart + xMagNew;
		}
//z--
		double zMagOld = (SRWRadStructAccessData.nz - 1)*SRWRadStructAccessData.zStep;
		if(RadResizeStruct.pzm > 1.)
		{
			long NzInner = (RadResizeStruct.pzd != 1.)? long(NewNz/RadResizeStruct.pzm) : SRWRadStructAccessData.nz;

			long NewNz_mi_NzInner = NewNz - NzInner;
			long NzOuterLeft = NewNz_mi_NzInner >> 1;
			if((NzOuterLeft << 1) != NewNz_mi_NzInner) NzOuterLeft++;
		
			NewSRWRadStructAccessData.zStep = zMagOld/(NzInner - 1);
			NewSRWRadStructAccessData.zStart = SRWRadStructAccessData.zStart - NzOuterLeft*NewSRWRadStructAccessData.zStep;
			NewSRWRadStructAccessData.AuxLong3 = NzOuterLeft; // izStart
			NewSRWRadStructAccessData.AuxLong4 = NzOuterLeft + NzInner - 1; // izEnd

			NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'z');
		}
		else
		{
			double zEndOld = SRWRadStructAccessData.zStart + (SRWRadStructAccessData.nz - 1)*SRWRadStructAccessData.zStep;
			double zMid = 0.5*(SRWRadStructAccessData.zStart + zEndOld);
			double zMagNew = zMagOld*RadResizeStruct.pzm;
			NewSRWRadStructAccessData.zStep = zMagNew/(NewNz - 1);

			if(RadResizeStruct.pzd == 1)
			{
				NewSRWRadStructAccessData.zStep = SRWRadStructAccessData.zStep;
				zMagNew = (NewNz - 1)*NewSRWRadStructAccessData.zStep;
			}

			NewSRWRadStructAccessData.zStart = zMid - 0.5*zMagNew;
			NewSRWRadStructAccessData.AuxLong3 = 0; // izStart
			NewSRWRadStructAccessData.AuxLong4 = NewNz - 1; // izEnd

			if(NewSRWRadStructAccessData.zStart > NewSRWRadStructAccessData.zWfrMin) NewSRWRadStructAccessData.zWfrMin = NewSRWRadStructAccessData.zStart;
			if((NewSRWRadStructAccessData.zStart + zMagNew) < NewSRWRadStructAccessData.zWfrMax) NewSRWRadStructAccessData.zWfrMax = NewSRWRadStructAccessData.zStart + zMagNew;
		}
	}

	//long TotAmOfOldData = (SRWRadStructAccessData.ne*SRWRadStructAccessData.nx*SRWRadStructAccessData.nz) << 1;
	//long TotAmOfNewData = (NewSRWRadStructAccessData.ne*NewSRWRadStructAccessData.nx*NewSRWRadStructAccessData.nz) << 1;
	long long TotAmOfOldData = (((long long)SRWRadStructAccessData.ne)*((long long)SRWRadStructAccessData.nx)*((long long)SRWRadStructAccessData.nz)) << 1;
	long long TotAmOfNewData = (((long long)NewSRWRadStructAccessData.ne)*((long long)NewSRWRadStructAccessData.nx)*((long long)NewSRWRadStructAccessData.nz)) << 1;

	//char TreatPolarizSepar = 0;
	char TreatPolarizSepar = (!ExIsOK) || (!EzIsOK); //OC13112011
	if(!TreatPolarizSepar)
	{
		long nxCurRad = SRWRadStructAccessData.nx, nzCurRad = SRWRadStructAccessData.nz;
		double ExtraMemSize = ExtraMemSizeForResize(nxCurRad, nzCurRad, pxmIn, pxdIn, pzmIn, pzdIn, 1);
		double MemoryAvail = CheckMemoryAvailable();

		double SecurityMemResizeCoef = 0.95; //to tune
		TreatPolarizSepar = (ExtraMemSize > SecurityMemResizeCoef*MemoryAvail)? 1 : 0;
	}

	float *OldRadXCopy = 0, *OldRadZCopy = 0;
	if(!TreatPolarizSepar)
	{
		//if(pxmIn*pxdIn*pzmIn*pzdIn >= 1.)
		//if((pxmIn*pxdIn*pzmIn*pzdIn >= 1.) || (SRWRadStructAccessData.m_newExtWfrCreateNotAllowed)) //OC140311
		if(pxmIn*pxdIn*pzmIn*pzdIn >= 1.) //OC161115
		{//Is this part really necessary?
			OldRadXCopy = new float[TotAmOfOldData];
			if(OldRadXCopy == 0) return MEMORY_ALLOCATION_FAILURE;
			OldRadZCopy = new float[TotAmOfOldData];
			if(OldRadZCopy == 0) return MEMORY_ALLOCATION_FAILURE;
			
			float *tOldRadXCopy = OldRadXCopy, *tOldRadZCopy = OldRadZCopy;
			float *tBaseRadX = SRWRadStructAccessData.pBaseRadX, *tBaseRadZ = SRWRadStructAccessData.pBaseRadZ;
			//for(long i=0; i<TotAmOfOldData; i++) 
			for(long long i=0; i<TotAmOfOldData; i++) 
			{
				*(tOldRadXCopy++) = *(tBaseRadX++);
				*(tOldRadZCopy++) = *(tBaseRadZ++);
			}
			
			if(RadShouldBeChanged)
			{
				if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
				{
					if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz()) return result;
				}
				//else if(result = Send.ModifyRadNeNxNz(NewSRWRadStructAccessData)) return result;
				else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz()) return result;
			}
			
			tBaseRadX = NewSRWRadStructAccessData.pBaseRadX;
			tBaseRadZ = NewSRWRadStructAccessData.pBaseRadZ;
			//for(long j=0; j<TotAmOfNewData; j++)
			for(long long j=0; j<TotAmOfNewData; j++)
			{
				*(tBaseRadX++) = 0.; *(tBaseRadZ++) = 0.; 
			}
			
			SRWRadStructAccessData.pBaseRadX = OldRadXCopy;
			SRWRadStructAccessData.pBaseRadZ = OldRadZCopy;
			
			if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct)) return result;
			
			if(OldRadXCopy != 0) delete[] OldRadXCopy;
			if(OldRadZCopy != 0) delete[] OldRadZCopy;
		}
		else
		{
			srTSRWRadStructWaveNames RadStructNames, OldRadStructNames;

			if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
			{
				if(result = NewSRWRadStructAccessData.AllocBaseRadAccordingToNeNxNz()) return result;
			}
			else
			{
#ifdef __IGOR_PRO__
				//if(result = Send.GetRadStructNames(SRWRadStructAccessData, RadStructNames)) return result;
				if(result = SRWRadStructAccessData.GetWfrStructNames(RadStructNames)) return result;
				OldRadStructNames = RadStructNames;
				
				char AuxRadName[] = "SrwWfrAux_rad\0", AuxExName[] = "SrwWfrAuxX_rae\0", AuxEzName[] = "SrwWfrAuxZ_rae\0";
				strcpy(RadStructNames.NameRad, AuxRadName);
				strcpy(RadStructNames.NameRadX, AuxExName);
				strcpy(RadStructNames.NameRadZ, AuxEzName);
				
				NewSRWRadStructAccessData.wRad = NIL;
				NewSRWRadStructAccessData.pBaseRadX = 0; NewSRWRadStructAccessData.wRadX = NIL;
				NewSRWRadStructAccessData.pBaseRadZ = 0; NewSRWRadStructAccessData.wRadZ = NIL;
				//if(result = Send.CreateNewRadStruct(NewSRWRadStructAccessData, RadStructNames)) return result;
				if(result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;		
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz(0, true)) return result;
				//OCTEST
				//SRWRadStructAccessData.pBaseRadX = NewSRWRadStructAccessData.pBaseRadXaux;
				//SRWRadStructAccessData.pBaseRadZ = NewSRWRadStructAccessData.pBaseRadZaux;
#endif
			}
			
			float *tRadX = NewSRWRadStructAccessData.pBaseRadX, *tRadZ = NewSRWRadStructAccessData.pBaseRadZ;
			//for(long j=0; j<TotAmOfNewData; j++)
			for(long long j=0; j<TotAmOfNewData; j++)
			{
				*(tRadX++) = 0.; *(tRadZ++) = 0.; 
			}

			if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct)) return result;

			if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
			{
				SRWRadStructAccessData.DeAllocBaseRadAccordingToNeNxNz();
			}
			else
			{
#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys RadKeys;
				RadKeys.wRad_= RadKeys.wRadX_= RadKeys.wRadZ_= 1;
				//if(result = Send.DeleteRadStructWaves(SRWRadStructAccessData, RadKeys)) return result;
				if(result = SRWRadStructAccessData.DeleteWfrStructWaves(RadKeys)) return result;
				
				//if(result = Send.RenameRadStruct(NewSRWRadStructAccessData, OldRadStructNames)) return result;
				if(result = NewSRWRadStructAccessData.RenameWfrStruct(OldRadStructNames)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if(result = NewSRWRadStructAccessData.DeleteWfrBackupData()) return result;
#endif
			}
		}
	}
	else //TreatPolarizSepar
	{
		//if(pxmIn*pxdIn*pzmIn*pzdIn >= 1.)
		//if((pxmIn*pxdIn*pzmIn*pzdIn >= 1.) || (SRWRadStructAccessData.m_newExtWfrCreateNotAllowed)) //OC140311
		if(pxmIn*pxdIn*pzmIn*pzdIn >= 1.) //OC161115
		{//Is this part necessary at all?
			if(ExIsOK) //OC13112011
			{
				OldRadXCopy = new float[TotAmOfOldData];
				if(OldRadXCopy == 0) return MEMORY_ALLOCATION_FAILURE;

				float *tOldRadXCopy = OldRadXCopy;
				float *tBaseRadX = SRWRadStructAccessData.pBaseRadX;
				//for(long i=0; i<TotAmOfOldData; i++) 
				for(long long i=0; i<TotAmOfOldData; i++) 
				{
					*(tOldRadXCopy++) = *(tBaseRadX++);
				}
				if(RadShouldBeChanged)
				{
					if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
					{
						if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz('x')) return result;
					}
					//else if(result = Send.ModifyRadNeNxNz(NewSRWRadStructAccessData, 'x')) return result;
					else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('x')) return result;
				}
				tBaseRadX = NewSRWRadStructAccessData.pBaseRadX;
				//for(long j=0; j<TotAmOfNewData; j++) 
				for(long long j=0; j<TotAmOfNewData; j++) 
				{
					*(tBaseRadX++) = 0.;
				}
				SRWRadStructAccessData.pBaseRadX = OldRadXCopy;
				if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'x')) return result;
				if(OldRadXCopy != 0) delete[] OldRadXCopy;
			}
			if(EzIsOK)
			{
				OldRadZCopy = new float[TotAmOfOldData];
				if(OldRadZCopy == 0) return MEMORY_ALLOCATION_FAILURE;

				float *tOldRadZCopy = OldRadZCopy;
				float *tBaseRadZ = SRWRadStructAccessData.pBaseRadZ;
				//for(long i=0; i<TotAmOfOldData; i++) 
				for(long long i=0; i<TotAmOfOldData; i++) 
				{
					float testVal = *(tBaseRadZ++);
					*(tOldRadZCopy++) = testVal;
					//*(tOldRadZCopy++) = *(tBaseRadZ++);
				}
				if(RadShouldBeChanged)
				{
					if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
					{
						if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz('z')) return result;
					}
					//else if(result = Send.ModifyRadNeNxNz(NewSRWRadStructAccessData, 'z')) return result;
					else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('z')) return result;
				}
				tBaseRadZ = NewSRWRadStructAccessData.pBaseRadZ;
				//for(long j=0; j<TotAmOfNewData; j++) 
				for(long long j=0; j<TotAmOfNewData; j++) 
				{
					*(tBaseRadZ++) = 0.;
				}
				SRWRadStructAccessData.pBaseRadZ = OldRadZCopy;
				if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'z')) return result;
				if(OldRadZCopy != 0) delete[] OldRadZCopy;
			}
		}
		else
		{
#ifdef __IGOR_PRO__
			srTSRWRadStructWaveNames RadStructNames;
			//if(result = Send.GetRadStructNames(SRWRadStructAccessData, RadStructNames)) return result;
			if(result = SRWRadStructAccessData.GetWfrStructNames(RadStructNames)) return result;			
			srTSRWRadStructWaveNames OldRadStructNames = RadStructNames;

			char AuxRadName[] = "SrwWfrAux_rad\0", AuxExName[] = "SrwWfrAuxX_rae\0", AuxEzName[] = "SrwWfrAuxZ_rae\0";
			strcpy(RadStructNames.NameRad, AuxRadName);
			strcpy(RadStructNames.NameRadZ, AuxEzName);
			strcpy(RadStructNames.NameRadX, AuxExName);

			NewSRWRadStructAccessData.wRad = NIL;
#endif
//Ex
			if(ExIsOK)
			{
#ifdef __IGOR_PRO__
				NewSRWRadStructAccessData.pBaseRadX = 0; NewSRWRadStructAccessData.wRadX = NIL;
				//if(result = Send.CreateNewRadStruct(NewSRWRadStructAccessData, RadStructNames)) return result;
				if(result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;		
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if (result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('x', true)) return result;
#endif
				float *tRadX = NewSRWRadStructAccessData.pBaseRadX;
				//for(long j=0; j<TotAmOfNewData; j++) *(tRadX++) = 0.;
				for(long long j=0; j<TotAmOfNewData; j++) *(tRadX++) = 0.;

				if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'x')) return result;

#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys Keys;
				//Keys.wRadX_= 1;
				Keys.wRad_ = Keys.wRadX_= 1; Keys.wRadZ_ = 0; //OC161115
				//if(result = Send.DeleteRadStructWaves(SRWRadStructAccessData, Keys)) return result;
				if(result = SRWRadStructAccessData.DeleteWfrStructWaves(Keys)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if (result = NewSRWRadStructAccessData.DeleteWfrBackupData('x')) return result;
#endif
			}
//Ez
			if(EzIsOK)
			{
#ifdef __IGOR_PRO__
				NewSRWRadStructAccessData.pBaseRadZ = 0; NewSRWRadStructAccessData.wRadZ = NIL;
				//if(result = Send.CreateNewRadStruct(NewSRWRadStructAccessData, RadStructNames)) return result;
				if (result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('z', true)) return result;
#endif
				float *tRadZ = NewSRWRadStructAccessData.pBaseRadZ;
				//for(long i = 0; i < TotAmOfNewData; i++) *(tRadZ++) = 0.;
				for(long long i = 0; i < TotAmOfNewData; i++) *(tRadZ++) = 0.;

				if(result = RadResizeCore(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'z')) return result;

#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys Keys;
				Keys.wRadX_ = 0; Keys.wRad_ = Keys.wRadZ_ = 1;
				//if(result = Send.DeleteRadStructWaves(SRWRadStructAccessData, Keys)) return result;
				if (result = SRWRadStructAccessData.DeleteWfrStructWaves(Keys)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if (result = NewSRWRadStructAccessData.DeleteWfrBackupData('z')) return result;
#endif
			}
#ifdef __IGOR_PRO__
			//if(result = Send.RenameRadStruct(NewSRWRadStructAccessData, OldRadStructNames)) return result;
			if(result = NewSRWRadStructAccessData.RenameWfrStruct(OldRadStructNames)) return result;
#endif
		}
	}

	SRWRadStructAccessData = NewSRWRadStructAccessData;
	NewSRWRadStructAccessData.ZeroPtrs();

	//if(RadResizeStruct.UseOtherSideFFT)
	if(RadResizeStruct.useOtherSideFFT()) //OC090311
	{
		char ToRepres = (SRWRadStructAccessData.Pres == 0)? 1 : 0;
		if(result = SetRadRepres(&SRWRadStructAccessData, ToRepres)) return result;

		if(UseStartTrWithOtherSide)
		{
			double xcOld = xStartOld + (nxOld >> 1)*xStepOld;
			SRWRadStructAccessData.xStart = xcOld - (SRWRadStructAccessData.nx >> 1)*SRWRadStructAccessData.xStep;
			double zcOld = zStartOld + (nzOld >> 1)*zStepOld;
			SRWRadStructAccessData.zStart = zcOld - (SRWRadStructAccessData.nz >> 1)*SRWRadStructAccessData.zStep;
			SRWRadStructAccessData.SetNonZeroWavefrontLimitsToFullRange();
		}

		double pxmNew = RadResizeStruct.pxd, pxdNew = RadResizeStruct.pxm;
		double pzmNew = RadResizeStruct.pzd, pzdNew = RadResizeStruct.pzm;

		RadResizeStruct.pxm = pxmNew; RadResizeStruct.pxd = pxdNew;
		RadResizeStruct.pzm = pzmNew; RadResizeStruct.pzd = pzdNew;
	}

	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeCore(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp)
{
	const double RelStepTol = 1.e-05; // To steer
	char OnlyMakeLargerRange = ((RadResizeStruct.pxd == 1.) && (RadResizeStruct.pzd == 1.) && (RadResizeStruct.pxm >= 1.) && (RadResizeStruct.pzm >= 1.)) 
							&& ((::fabs(OldRadAccessData.xStep - NewRadAccessData.xStep) < RelStepTol) && (::fabs(OldRadAccessData.zStep - NewRadAccessData.zStep) < RelStepTol));

	if(OnlyMakeLargerRange) return RadResizeCore_OnlyLargerRange(OldRadAccessData, NewRadAccessData, RadResizeStruct, PolComp);

	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x'));
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z'));

	const double DistAbsTol = 1.E-10;

	int ixStart = int(NewRadAccessData.AuxLong1);
	int ixEnd = int(NewRadAccessData.AuxLong2);
	int izStart = int(NewRadAccessData.AuxLong3);
	int izEnd = int(NewRadAccessData.AuxLong4);

	char WaveFrontTermWasTreated = 0;
	bool OrigWfrQuadTermCanBeTreatedAtResizeX = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX;
	bool OrigWfrQuadTermCanBeTreatedAtResizeZ = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ;

	//if((!RadResizeStruct.DoNotTreatSpherTerm) && WaveFrontTermCanBeTreated(OldRadAccessData))
	if((!RadResizeStruct.doNotTreatSpherTerm()) && WaveFrontTermCanBeTreated(OldRadAccessData)) //OC090311
	{
		NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX;
		NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ;

		TreatStronglyOscillatingTerm(OldRadAccessData, 'r', PolComp);

		NewRadAccessData.wfrReffX = OldRadAccessData.wfrReffX;
		NewRadAccessData.wfrReffZ = OldRadAccessData.wfrReffZ;

		WaveFrontTermWasTreated = 1;
	}

	double xStepInvOld = 1./OldRadAccessData.xStep;
	double zStepInvOld = 1./OldRadAccessData.zStep;
	int nx_mi_1Old = OldRadAccessData.nx - 1;
	int nz_mi_1Old = OldRadAccessData.nz - 1;
	int nx_mi_2Old = nx_mi_1Old - 1;
	int nz_mi_2Old = nz_mi_1Old - 1;

	srTInterpolAux01 InterpolAux01;

	srTInterpolAux02 InterpolAux02[4], InterpolAux02I[2];
	srTInterpolAuxF AuxF[4], AuxFI[2];
	int ixStOld, izStOld, ixStOldPrev = -1000, izStOldPrev = -1000;

	float *pEX0_New = 0, *pEZ0_New = 0;
	if(TreatPolCompX) pEX0_New = NewRadAccessData.pBaseRadX;
	if(TreatPolCompZ) pEZ0_New = NewRadAccessData.pBaseRadZ;

	//long PerX_New = NewRadAccessData.ne << 1;
	//long PerZ_New = PerX_New*NewRadAccessData.nx;
	long long PerX_New = NewRadAccessData.ne << 1;
	long long PerZ_New = PerX_New*NewRadAccessData.nx;

	//long PerX_Old = PerX_New;
	//long PerZ_Old = PerX_Old*OldRadAccessData.nx;
	long long PerX_Old = PerX_New;
	long long PerZ_Old = PerX_Old*OldRadAccessData.nx;

	float BufF[4], BufFI[2];
	char UseLowOrderInterp_PolCompX, UseLowOrderInterp_PolCompZ;

	int result = 0;
	for(int ie=0; ie<NewRadAccessData.ne; ie++)
	{
		ixStOldPrev = -1000; izStOldPrev = -1000;

		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;
		for(int iz=izStart; iz<=izEnd; iz++)
		{
			if(result = srYield.Check()) return result;

			double zAbs = NewRadAccessData.zStart + iz*NewRadAccessData.zStep;

			char FieldShouldBeZeroedDueToZ = 0;
			if(NewRadAccessData.WfrEdgeCorrShouldBeDone)
			{
				if((zAbs < NewRadAccessData.zWfrMin - DistAbsTol) || (zAbs > NewRadAccessData.zWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToZ = 1;
			}

			int izcOld = int((zAbs - OldRadAccessData.zStart)*zStepInvOld + 1.E-06);

			double zRel = zAbs - (OldRadAccessData.zStart + izcOld*OldRadAccessData.zStep);

			if(izcOld == nz_mi_1Old) { izStOld = izcOld - 3; zRel += 2.*OldRadAccessData.zStep;}
			else if(izcOld == nz_mi_2Old) { izStOld = izcOld - 2; zRel += OldRadAccessData.zStep;}
			else if(izcOld == 0) { izStOld = izcOld; zRel -= OldRadAccessData.zStep;}
			else izStOld = izcOld - 1;

			zRel *= zStepInvOld;

			int izcOld_mi_izStOld = izcOld - izStOld;
			//long izPerZ_New = iz*PerZ_New;
			long long izPerZ_New = iz*PerZ_New;

			float *pEX_StartForX_New = 0, *pEZ_StartForX_New = 0;
			if(TreatPolCompX) pEX_StartForX_New = pEX0_New + izPerZ_New;
			if(TreatPolCompZ) pEZ_StartForX_New = pEZ0_New + izPerZ_New;

			for(int ix=ixStart; ix<=ixEnd; ix++)
			{
				//long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
				long long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
				float *pEX_New = 0, *pEZ_New = 0;
				if(TreatPolCompX) pEX_New = pEX_StartForX_New + ixPerX_New_p_Two_ie;
				if(TreatPolCompZ) pEZ_New = pEZ_StartForX_New + ixPerX_New_p_Two_ie;

				double xAbs = NewRadAccessData.xStart + ix*NewRadAccessData.xStep;

				char FieldShouldBeZeroedDueToX = 0;
				if(NewRadAccessData.WfrEdgeCorrShouldBeDone)
				{
					if((xAbs < NewRadAccessData.xWfrMin - DistAbsTol) || (xAbs > NewRadAccessData.xWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToX = 1;
				}
				char FieldShouldBeZeroed = (FieldShouldBeZeroedDueToX || FieldShouldBeZeroedDueToZ);

				int ixcOld = int((xAbs - OldRadAccessData.xStart)*xStepInvOld + 1.E-06);
				double xRel = xAbs - (OldRadAccessData.xStart + ixcOld*OldRadAccessData.xStep);

				if(ixcOld == nx_mi_1Old) { ixStOld = ixcOld - 3; xRel += 2.*OldRadAccessData.xStep;}
				else if(ixcOld == nx_mi_2Old) { ixStOld = ixcOld - 2; xRel += OldRadAccessData.xStep;}
				else if(ixcOld == 0) { ixStOld = ixcOld; xRel -= OldRadAccessData.xStep;}
				else ixStOld = ixcOld - 1;

				xRel *= xStepInvOld;

				int ixcOld_mi_ixStOld = ixcOld - ixStOld;

				if((izStOld != izStOldPrev) || (ixStOld != ixStOldPrev))
				{
					UseLowOrderInterp_PolCompX = 0, UseLowOrderInterp_PolCompZ = 0;

					//long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old + Two_ie;
					long long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old + Two_ie;

					if(TreatPolCompX)
					{
						float* pExSt_Old = OldRadAccessData.pBaseRadX + TotOffsetOld;
						GetCellDataForInterpol(pExSt_Old, PerX_Old, PerZ_Old, AuxF);

						SetupCellDataI(AuxF, AuxFI);
						UseLowOrderInterp_PolCompX = CheckForLowOrderInterp(AuxF, AuxFI, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02, InterpolAux02I);

						if(!UseLowOrderInterp_PolCompX)
						{
							for(int i=0; i<2; i++) 
							{
								SetupInterpolAux02(AuxF + i, &InterpolAux01, InterpolAux02 + i);
							}
							SetupInterpolAux02(AuxFI, &InterpolAux01, InterpolAux02I);
						}
					}
					if(TreatPolCompZ)
					{
						float* pEzSt_Old = OldRadAccessData.pBaseRadZ + TotOffsetOld;
						GetCellDataForInterpol(pEzSt_Old, PerX_Old, PerZ_Old, AuxF+2);

						SetupCellDataI(AuxF+2, AuxFI+1);
						UseLowOrderInterp_PolCompZ = CheckForLowOrderInterp(AuxF+2, AuxFI+1, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02+2, InterpolAux02I+1);

						if(!UseLowOrderInterp_PolCompZ)
						{
							for(int i=0; i<2; i++) 
							{
								SetupInterpolAux02(AuxF+2+i, &InterpolAux01, InterpolAux02+2+i);
							}
							SetupInterpolAux02(AuxFI+1, &InterpolAux01, InterpolAux02I+1);
						}
					}

					ixStOldPrev = ixStOld; izStOldPrev = izStOld;
				}

				if(TreatPolCompX)
				{
					if(UseLowOrderInterp_PolCompX) 
					{
						InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 0);
						InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 0);
					}
					else
					{
						InterpolF(InterpolAux02, xRel, zRel, BufF, 0);
						InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 0);
					}

					(*BufFI) *= AuxFI->fNorm;
					ImproveReAndIm(BufF, BufFI);

					if(FieldShouldBeZeroed)
					{
						*BufF = 0.; *(BufF+1) = 0.;
					}

					*pEX_New = *BufF;
					*(pEX_New+1) = *(BufF+1);
				}
				if(TreatPolCompZ)
				{
					if(UseLowOrderInterp_PolCompZ) 
					{
						InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 2);
						InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 1);
					}
					else
					{
						InterpolF(InterpolAux02, xRel, zRel, BufF, 2);
						InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 1);
					}

					(*(BufFI+1)) *= (AuxFI+1)->fNorm;
					ImproveReAndIm(BufF+2, BufFI+1);

					if(FieldShouldBeZeroed)
					{
						*(BufF+2) = 0.; *(BufF+3) = 0.;
					}

					*pEZ_New = *(BufF+2);
					*(pEZ_New+1) = *(BufF+3);
				}
			}
		}
	}
	if(WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(NewRadAccessData, 'a', PolComp);

	OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;
	NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;

	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeGenE(srTSRWRadStructAccessData& SRWRadStructAccessData, srTRadResize& RadResizeStruct)
{
	if((RadResizeStruct.pem == 1.) && (RadResizeStruct.ped == 1.)) return 0;
	int result = 0;

	bool ExIsOK = SRWRadStructAccessData.pBaseRadX != 0;
	bool EzIsOK = SRWRadStructAccessData.pBaseRadZ != 0;

	long neOld = SRWRadStructAccessData.ne;
	double eStartOld = SRWRadStructAccessData.eStart;
	double eStepOld = SRWRadStructAccessData.eStep;
	double eRangeOld = eStepOld*neOld;
	char UseStartTrWithOtherSide = (::fabs(eStartOld + 0.5*eRangeOld) > 0.5*SRWRadStructAccessData.eStep);
	double pemIn = RadResizeStruct.pem, pedIn = RadResizeStruct.ped;

	if(RadResizeStruct.useOtherSideFFT())
	{
		char ToRepres = (SRWRadStructAccessData.PresT == 0)? 1 : 0;
		if(UseStartTrWithOtherSide)
		{
			SRWRadStructAccessData.eStart = -(SRWRadStructAccessData.ne >> 1)*SRWRadStructAccessData.eStep;
			//double eShift = SRWRadStructAccessData.eStart - eStartOld;
			//SRWRadStructAccessData.eWfrMin += eShift; SRWRadStructAccessData.eWfrMax += eShift;
		}
		//if(result = SetRadRepres(&SRWRadStructAccessData, ToRepres)) return result;
		if(result = SRWRadStructAccessData.SetRepresFT(ToRepres)) return result;

		double pemNew = RadResizeStruct.ped, pedNew = RadResizeStruct.pem;
		RadResizeStruct.pem = pemNew; RadResizeStruct.ped = pedNew;
	}

	srTSRWRadStructAccessData NewSRWRadStructAccessData(SRWRadStructAccessData);

	double peTot = RadResizeStruct.pem*RadResizeStruct.ped;
	long NewNe = SRWRadStructAccessData.ne; 

	CGenMathFFT1D FFT;
	char RadShouldBeChanged = (peTot != 1.);
	if(RadShouldBeChanged)
	{
		double peTotNewNe = peTot*NewNe;
		long Long_peTotNewNe = long(peTotNewNe);
		NewNe = ((peTotNewNe - Long_peTotNewNe) >= 0.5)? Long_peTotNewNe + 1 : Long_peTotNewNe;

		FFT.NextCorrectNumberForFFT(NewNe);
		NewSRWRadStructAccessData.ne = NewNe;
	}

	char CenterIsOffset = 0;
	long iEc = (SRWRadStructAccessData.ne >> 1);

	if(::fabs(RadResizeStruct.RelCenPosE - 0.5) > RadResizeStruct.RelCenPosTol)
	{
		double eRangeOld = SRWRadStructAccessData.ne*SRWRadStructAccessData.eStep;
		double ecAppr = SRWRadStructAccessData.eStart + eRangeOld*RadResizeStruct.RelCenPosE;
		iEc = IntegerOffsetCoord(SRWRadStructAccessData.eStart, SRWRadStructAccessData.eStep, ecAppr);
		CenterIsOffset = 1;
	}
	double Ec = SRWRadStructAccessData.eStart + SRWRadStructAccessData.eStep*iEc;

	if(CenterIsOffset)
	{
		long NewNeLeftHalf = NewNe >> 1;
		NewSRWRadStructAccessData.eStep = SRWRadStructAccessData.eStep/RadResizeStruct.ped;
		NewSRWRadStructAccessData.eStart = Ec - NewNeLeftHalf*NewSRWRadStructAccessData.eStep;

		double NewNeEndInterp = NewSRWRadStructAccessData.eStart + (NewNe - 1)*NewSRWRadStructAccessData.eStep;
		double OldNeEndInterp = SRWRadStructAccessData.eStart + (SRWRadStructAccessData.ne - 1)*SRWRadStructAccessData.eStep;

		TuneStepToKeepInterpLimitsTheSameAtResize(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'e', iEc);

		if(NewSRWRadStructAccessData.eStart < SRWRadStructAccessData.eStart) // ieStart
			NewSRWRadStructAccessData.AuxLong1 = IntegerOffsetCoord(NewSRWRadStructAccessData.eStart, NewSRWRadStructAccessData.eStep, SRWRadStructAccessData.eStart);
		else NewSRWRadStructAccessData.AuxLong1 = 0;

		if(NewNeEndInterp > OldNeEndInterp) // ieEnd
			NewSRWRadStructAccessData.AuxLong2 = IntegerOffsetCoord(NewSRWRadStructAccessData.eStart, NewSRWRadStructAccessData.eStep, OldNeEndInterp);
		else NewSRWRadStructAccessData.AuxLong2 = NewNe - 1;

		//if(RadResizeStruct.pem > 1.) NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'e'); //not necessary for e
		//else
		//{
		//	if(NewSRWRadStructAccessData.eStart > NewSRWRadStructAccessData.eWfrMin) NewSRWRadStructAccessData.eWfrMin = NewSRWRadStructAccessData.eStart;
		//	double eEndNew = NewSRWRadStructAccessData.eStart + NewSRWRadStructAccessData.ne*NewSRWRadStructAccessData.eStep;
		//	if(eEndNew < NewSRWRadStructAccessData.eWfrMax) NewSRWRadStructAccessData.eWfrMax = eEndNew;
		//}
	}
	else
	{
		double eMagOld = (SRWRadStructAccessData.ne - 1)*SRWRadStructAccessData.eStep;
		if(RadResizeStruct.pem > 1.)
		{
			long NeInner = (RadResizeStruct.ped != 1.)? long(NewNe/RadResizeStruct.pem) : SRWRadStructAccessData.ne;
			long NewNe_mi_NeInner = NewNe - NeInner;
			long NeOuterLeft = NewNe_mi_NeInner >> 1;
			if((NeOuterLeft << 1) != NewNe_mi_NeInner) NeOuterLeft++;

			NewSRWRadStructAccessData.eStep = eMagOld/(NeInner - 1);
			NewSRWRadStructAccessData.eStart = SRWRadStructAccessData.eStart - NeOuterLeft*NewSRWRadStructAccessData.eStep;
			NewSRWRadStructAccessData.AuxLong1 = NeOuterLeft; // ieStart
			NewSRWRadStructAccessData.AuxLong2 = NeOuterLeft + NeInner - 1; // ixEnd

			//NewSRWRadStructAccessData.PreserveLogicsOfWfrLimitsAtRangeResizing(&SRWRadStructAccessData, 'e'); //not necessary for e
		}
		else
		{
			double eEndOld = SRWRadStructAccessData.eStart + (SRWRadStructAccessData.ne - 1)*SRWRadStructAccessData.eStep;
			double eMid = 0.5*(SRWRadStructAccessData.eStart + eEndOld);
			double eMagNew = eMagOld*RadResizeStruct.pem;
			NewSRWRadStructAccessData.eStep = eMagNew/(NewNe - 1);
			if(RadResizeStruct.ped == 1)
			{
				NewSRWRadStructAccessData.eStep = SRWRadStructAccessData.eStep;
				eMagNew = (NewNe - 1)*NewSRWRadStructAccessData.eStep;
			}
			NewSRWRadStructAccessData.eStart = eMid - 0.5*eMagNew;
			NewSRWRadStructAccessData.AuxLong1 = 0; // ieStart
			NewSRWRadStructAccessData.AuxLong2 = NewNe - 1; // ieEnd

			//if(NewSRWRadStructAccessData.eStart > NewSRWRadStructAccessData.eWfrMin) NewSRWRadStructAccessData.eWfrMin = NewSRWRadStructAccessData.eStart; //not necessary for e
			//if((NewSRWRadStructAccessData.eStart + eMagNew) < NewSRWRadStructAccessData.eWfrMax) NewSRWRadStructAccessData.eWfrMax = NewSRWRadStructAccessData.eStart + eMagNew;
		}
	}

	//long TotAmOfOldData = (SRWRadStructAccessData.ne*SRWRadStructAccessData.nx*SRWRadStructAccessData.nz) << 1;
	//long TotAmOfNewData = (NewSRWRadStructAccessData.ne*NewSRWRadStructAccessData.nx*NewSRWRadStructAccessData.nz) << 1;
	long long TotAmOfOldData = (((long long)SRWRadStructAccessData.ne)*((long long)SRWRadStructAccessData.nx)*((long long)SRWRadStructAccessData.nz)) << 1;
	long long TotAmOfNewData = (((long long)NewSRWRadStructAccessData.ne)*((long long)NewSRWRadStructAccessData.nx)*((long long)NewSRWRadStructAccessData.nz)) << 1;

	char TreatPolarizSepar = (!ExIsOK) || (!EzIsOK);
	if(!TreatPolarizSepar)
	{
		long neCurRad = SRWRadStructAccessData.ne;
		long nxCurRad = SRWRadStructAccessData.nx;
		long nzCurRad = SRWRadStructAccessData.nz;
		double ExtraMemSize = ExtraMemSizeForResizeE(neCurRad, nxCurRad, nzCurRad, pemIn, pedIn, 1);
		double MemoryAvail = CheckMemoryAvailable();
		double SecurityMemResizeCoef = 0.95; //to tune
		TreatPolarizSepar = (ExtraMemSize > SecurityMemResizeCoef*MemoryAvail)? 1 : 0;
	}

	float *OldRadXCopy = 0, *OldRadZCopy = 0;
	if(!TreatPolarizSepar)
	{
		//if((pemIn*pedIn >= 1.) || SRWRadStructAccessData.m_newExtWfrCreateNotAllowed)
		if(pemIn*pedIn >= 1.) //OC161115
		{//Is this part necessary at all?
			OldRadXCopy = new float[TotAmOfOldData];
			if(OldRadXCopy == 0) return MEMORY_ALLOCATION_FAILURE;
			OldRadZCopy = new float[TotAmOfOldData];
			if(OldRadZCopy == 0) return MEMORY_ALLOCATION_FAILURE;

			float *tOldRadXCopy = OldRadXCopy, *tOldRadZCopy = OldRadZCopy;
			float *tBaseRadX = SRWRadStructAccessData.pBaseRadX, *tBaseRadZ = SRWRadStructAccessData.pBaseRadZ;
			//for(long i=0; i<TotAmOfOldData; i++) 
			for(long long i=0; i<TotAmOfOldData; i++) 
			{
				*(tOldRadXCopy++) = *(tBaseRadX++);
				*(tOldRadZCopy++) = *(tBaseRadZ++);
			}

			if(RadShouldBeChanged)
			{
				if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
				{
					if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz()) return result;
				}
				else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz()) return result;
			}

			tBaseRadX = NewSRWRadStructAccessData.pBaseRadX;
			tBaseRadZ = NewSRWRadStructAccessData.pBaseRadZ;
			//for(long j=0; j<TotAmOfNewData; j++)
			for(long long j=0; j<TotAmOfNewData; j++)
			{
				*(tBaseRadX++) = 0.; *(tBaseRadZ++) = 0.; 
			}

			SRWRadStructAccessData.pBaseRadX = OldRadXCopy;
			SRWRadStructAccessData.pBaseRadZ = OldRadZCopy;

			if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct)) return result;
			
			if(OldRadXCopy != 0) delete[] OldRadXCopy;
			if(OldRadZCopy != 0) delete[] OldRadZCopy;
		}
		else
		{//This processing saves more memory
			srTSRWRadStructWaveNames RadStructNames, OldRadStructNames;
			if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
			{
				if(result = NewSRWRadStructAccessData.AllocBaseRadAccordingToNeNxNz()) return result;
			}
			else
			{
#ifdef __IGOR_PRO__
				if(result = SRWRadStructAccessData.GetWfrStructNames(RadStructNames)) return result;
				OldRadStructNames = RadStructNames;

				char AuxRadName[] = "SrwWfrAux_rad\0", AuxExName[] = "SrwWfrAuxX_rae\0", AuxEzName[] = "SrwWfrAuxZ_rae\0";
				strcpy(RadStructNames.NameRad, AuxRadName);
				strcpy(RadStructNames.NameRadX, AuxExName);
				strcpy(RadStructNames.NameRadZ, AuxEzName);

				NewSRWRadStructAccessData.wRad = NIL;
				NewSRWRadStructAccessData.pBaseRadX = 0; NewSRWRadStructAccessData.wRadX = NIL;
				NewSRWRadStructAccessData.pBaseRadZ = 0; NewSRWRadStructAccessData.wRadZ = NIL;
				if(result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC141115
				if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz(0, true)) return result;
#endif
			}

			float *tRadX = NewSRWRadStructAccessData.pBaseRadX, *tRadZ = NewSRWRadStructAccessData.pBaseRadZ;
			//for(long j=0; j<TotAmOfNewData; j++)
			for(long long j=0; j<TotAmOfNewData; j++)
			{
				*(tRadX++) = 0.; *(tRadZ++) = 0.; 
			}

			if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct)) return result;

			if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
			{
				SRWRadStructAccessData.DeAllocBaseRadAccordingToNeNxNz();
			}
			else
			{
#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys RadKeys;
				RadKeys.wRad_= RadKeys.wRadX_= RadKeys.wRadZ_= 1;
				if(result = SRWRadStructAccessData.DeleteWfrStructWaves(RadKeys)) return result;
				
				if(result = NewSRWRadStructAccessData.RenameWfrStruct(OldRadStructNames)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED)
				if(result = NewSRWRadStructAccessData.DeleteWfrBackupData()) return result;
#endif
			}
		}
	}
	else //TreatPolarizSepar
	{
		//if((pemIn*pedIn >= 1.) || (SRWRadStructAccessData.m_newExtWfrCreateNotAllowed))
		if(pemIn*pedIn >= 1.) //OC161115
		{//Is this part necessary at all?
			if(ExIsOK)
			{
				OldRadXCopy = new float[TotAmOfOldData];
				if(OldRadXCopy == 0) return MEMORY_ALLOCATION_FAILURE;

				float *tOldRadXCopy = OldRadXCopy;
				float *tBaseRadX = SRWRadStructAccessData.pBaseRadX;
				//for(long i=0; i<TotAmOfOldData; i++) 
				for(long long i=0; i<TotAmOfOldData; i++) 
				{
					*(tOldRadXCopy++) = *(tBaseRadX++);
				}

				if(RadShouldBeChanged)
				{
					if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
					{
						if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz('x')) return result;
					}
					else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('x')) return result;
				}
				
				tBaseRadX = NewSRWRadStructAccessData.pBaseRadX;
				//for(long j=0; j<TotAmOfNewData; j++) *(tBaseRadX++) = 0.;
				for(long long j=0; j<TotAmOfNewData; j++) *(tBaseRadX++) = 0.;

				SRWRadStructAccessData.pBaseRadX = OldRadXCopy;
				if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'x')) return result;
				if(OldRadXCopy != 0) delete[] OldRadXCopy;
			}
			if(EzIsOK)
			{
				OldRadZCopy = new float[TotAmOfOldData];
				if(OldRadZCopy == 0) return MEMORY_ALLOCATION_FAILURE;

				float *tOldRadZCopy = OldRadZCopy;
				float *tBaseRadZ = SRWRadStructAccessData.pBaseRadZ;
				//for(long i=0; i<TotAmOfOldData; i++) 
				for(long long i=0; i<TotAmOfOldData; i++) 
				{
					//float testVal = *(tBaseRadZ++);
					//*(tOldRadZCopy++) = testVal;
					*(tOldRadZCopy++) = *(tBaseRadZ++);
				}

				if(RadShouldBeChanged)
				{
					if(NewSRWRadStructAccessData.BaseRadWasEmulated) 
					{
						if(result = NewSRWRadStructAccessData.ReAllocBaseRadAccordingToNeNxNz('z')) return result;
					}
					else if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('z')) return result;
				}

				tBaseRadZ = NewSRWRadStructAccessData.pBaseRadZ;
				//for(long j=0; j<TotAmOfNewData; j++) *(tBaseRadZ++) = 0.;
				for(long long j=0; j<TotAmOfNewData; j++) *(tBaseRadZ++) = 0.;

				SRWRadStructAccessData.pBaseRadZ = OldRadZCopy;
				if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'z')) return result;
				if(OldRadZCopy != 0) delete[] OldRadZCopy;
			}
		}
		else
		{
#ifdef __IGOR_PRO__
			srTSRWRadStructWaveNames RadStructNames;
			if(result = SRWRadStructAccessData.GetWfrStructNames(RadStructNames)) return result;			
			srTSRWRadStructWaveNames OldRadStructNames = RadStructNames;

			char AuxRadName[] = "SrwWfrAux_rad\0", AuxExName[] = "SrwWfrAuxX_rae\0", AuxEzName[] = "SrwWfrAuxZ_rae\0";
			strcpy(RadStructNames.NameRad, AuxRadName);
			strcpy(RadStructNames.NameRadZ, AuxEzName);
			strcpy(RadStructNames.NameRadX, AuxExName);

			NewSRWRadStructAccessData.wRad = NIL;
#endif
			if(ExIsOK)
			{
#ifdef __IGOR_PRO__
				NewSRWRadStructAccessData.pBaseRadX = 0; NewSRWRadStructAccessData.wRadX = NIL;
				if(result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;		
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if(result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('x', true)) return result;
#endif
				float *tRadX = NewSRWRadStructAccessData.pBaseRadX;
				//for(long j=0; j<TotAmOfNewData; j++) *(tRadX++) = 0.;
				for(long long j=0; j<TotAmOfNewData; j++) *(tRadX++) = 0.;

				if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'x')) return result;

#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys Keys; 
				Keys.wRad_= Keys.wRadX_= 1;
				if(result = SRWRadStructAccessData.DeleteWfrStructWaves(Keys)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if(result = NewSRWRadStructAccessData.DeleteWfrBackupData('x')) return result;
#endif
			}
			//if(ExIsOK)
			if(EzIsOK) //OC161115
			{
#ifdef __IGOR_PRO__
				NewSRWRadStructAccessData.pBaseRadZ = 0; NewSRWRadStructAccessData.wRadZ = NIL;
				if(result = NewSRWRadStructAccessData.CreateNewWfrStruct(RadStructNames)) return result;		
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if (result = NewSRWRadStructAccessData.ModifyWfrNeNxNz('z', true)) return result;
#endif
				float *tRadZ = NewSRWRadStructAccessData.pBaseRadZ;
				//for(long i=0; i<TotAmOfNewData; i++) *(tRadZ++) = 0.;
				for(long long i=0; i<TotAmOfNewData; i++) *(tRadZ++) = 0.;

				if(result = RadResizeCoreE(SRWRadStructAccessData, NewSRWRadStructAccessData, RadResizeStruct, 'z')) return result;

#ifdef __IGOR_PRO__
				srTSRWRadStructWaveKeys Keys; 
				Keys.wRad_= Keys.wRadZ_= 1;
				if(result = SRWRadStructAccessData.DeleteWfrStructWaves(Keys)) return result;
#endif
#if defined(SRWLIB_STATIC) || defined(SRWLIB_SHARED) //OC161115
				if (result = NewSRWRadStructAccessData.DeleteWfrBackupData('z')) return result;
#endif
			}
#ifdef __IGOR_PRO__
			if(result = NewSRWRadStructAccessData.RenameWfrStruct(OldRadStructNames)) return result;
#endif
		}
	}

	SRWRadStructAccessData = NewSRWRadStructAccessData;
	NewSRWRadStructAccessData.ZeroPtrs();

	if(RadResizeStruct.useOtherSideFFT())
	{
		char ToRepres = (SRWRadStructAccessData.PresT == 0)? 1 : 0;
		//if(result = SetRadRepres(&SRWRadStructAccessData, ToRepres)) return result;
		if(result = SRWRadStructAccessData.SetRepresFT(ToRepres)) return result;

		if(UseStartTrWithOtherSide)
		{
			double ecOld = eStartOld + (neOld >> 1)*eStepOld;
			SRWRadStructAccessData.eStart = ecOld - (SRWRadStructAccessData.ne >> 1)*SRWRadStructAccessData.eStep;
			SRWRadStructAccessData.SetNonZeroWavefrontLimitsToFullRange();
		}

		double pemNew = RadResizeStruct.ped, pedNew = RadResizeStruct.pem;
		RadResizeStruct.pem = pemNew; RadResizeStruct.ped = pedNew;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeCoreE(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp)
{
	const double RelStepTol = 1.e-05; // To steer
	char OnlyMakeLargerRange = ((RadResizeStruct.ped == 1.) && (RadResizeStruct.pem >= 1.)) 
							&& ((::fabs(OldRadAccessData.eStep - NewRadAccessData.eStep) < RelStepTol));

	if(OnlyMakeLargerRange) return RadResizeCore_OnlyLargerRangeE(OldRadAccessData, NewRadAccessData, RadResizeStruct, PolComp);

	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x')) && (OldRadAccessData.pBaseRadX != 0);
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z')) && (OldRadAccessData.pBaseRadZ != 0);

	const double DistAbsTol = 1.E-10;

	int ieStart = int(NewRadAccessData.AuxLong1);
	int ieEnd = int(NewRadAccessData.AuxLong2);

	char WaveFrontTermWasTreated = 0;
	bool OrigWfrQuadTermCanBeTreatedAtResizeX = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX;
	bool OrigWfrQuadTermCanBeTreatedAtResizeZ = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ;

	if((!RadResizeStruct.doNotTreatSpherTerm()) && WaveFrontTermCanBeTreated(OldRadAccessData)) //OC090311
	{
		NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX;
		NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ;

		TreatStronglyOscillatingTerm(OldRadAccessData, 'r', PolComp);

		NewRadAccessData.wfrReffX = OldRadAccessData.wfrReffX;
		NewRadAccessData.wfrReffZ = OldRadAccessData.wfrReffZ;

		WaveFrontTermWasTreated = 1;
	}

	double eStepInvOld = 1./OldRadAccessData.eStep;
	int ne_mi_1Old = OldRadAccessData.ne - 1;
	int ne_mi_2Old = ne_mi_1Old - 1;

	srTInterpolAux01_1D InterpolAux01;
	srTInterpolAux02_1D InterpolAux02[4], InterpolAux02I[2];
	srTInterpolAuxF_1D AuxF[4], AuxFI[2];

	int ieStOld, ieStOldPrev = -1000;

	float *pEX0_New = 0, *pEZ0_New = 0;
	if(TreatPolCompX) pEX0_New = NewRadAccessData.pBaseRadX;
	if(TreatPolCompZ) pEZ0_New = NewRadAccessData.pBaseRadZ;

	//long PerX_New = NewRadAccessData.ne << 1;
	//long PerZ_New = PerX_New*NewRadAccessData.nx;
	long long PerX_New = NewRadAccessData.ne << 1;
	long long PerZ_New = PerX_New*NewRadAccessData.nx;

	//long PerX_Old = OldRadAccessData.ne << 1;
	//long PerZ_Old = PerX_Old*OldRadAccessData.nx;
	long long PerX_Old = OldRadAccessData.ne << 1;
	long long PerZ_Old = PerX_Old*OldRadAccessData.nx;

	float BufF[4], BufFI[2];
	char UseLowOrderInterp_PolCompX, UseLowOrderInterp_PolCompZ;

	int result = 0;
	for(int iz=0; iz<NewRadAccessData.nz; iz++)
	{
		//long iz_PerZ_New = iz*PerZ_New;
		//long iz_PerZ_Old = iz*PerZ_Old;
		long long iz_PerZ_New = iz*PerZ_New;
		long long iz_PerZ_Old = iz*PerZ_Old;

		for(int ix=0; ix<NewRadAccessData.nx; ix++)
		{
			if(result = srYield.Check()) return result;

			//long iz_PerZ_New_p_ix_PerX_New = iz_PerZ_New + ix*PerX_New;
			//long iz_PerZ_Old_p_ix_PerX_Old = iz_PerZ_Old + ix*PerX_Old;
			long long iz_PerZ_New_p_ix_PerX_New = iz_PerZ_New + ix*PerX_New;
			long long iz_PerZ_Old_p_ix_PerX_Old = iz_PerZ_Old + ix*PerX_Old;

			ieStOldPrev = -1000;

			for(int ie=ieStart; ie<=ieEnd; ie++)
			{
				//long ofstNew = iz_PerZ_New_p_ix_PerX_New + (ie << 1);
				long long ofstNew = iz_PerZ_New_p_ix_PerX_New + (ie << 1);
				float *pEX_New = pEX0_New + ofstNew;
				float *pEZ_New = pEZ0_New + ofstNew;

				double eAbs = NewRadAccessData.eStart + ie*NewRadAccessData.eStep;

				//char FieldShouldBeZeroedDueToE = 0;
				//if(NewRadAccessData.WfrEdgeCorrShouldBeDone)
				//{
				//	if((eAbs < NewRadAccessData.eWfrMin - DistAbsTol) || (eAbs > NewRadAccessData.eWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToE = 1;
				//}
				//char FieldShouldBeZeroed = FieldShouldBeZeroedDueToE;

				int iecOld = int((eAbs - OldRadAccessData.eStart)*eStepInvOld + 1.E-08);
				double eRel = eAbs - (OldRadAccessData.eStart + iecOld*OldRadAccessData.eStep);

				if(iecOld == ne_mi_1Old) { ieStOld = iecOld - 3; eRel += 2.*OldRadAccessData.eStep;}
				else if(iecOld == ne_mi_2Old) { ieStOld = iecOld - 2; eRel += OldRadAccessData.eStep;}
				else if(iecOld == 0) { ieStOld = iecOld; eRel -= OldRadAccessData.eStep;}
				else ieStOld = iecOld - 1;

				eRel *= eStepInvOld;
				int iecOld_mi_ieStOld = iecOld - ieStOld;

				if(ieStOld != ieStOldPrev)
				{
					UseLowOrderInterp_PolCompX = 0, UseLowOrderInterp_PolCompZ = 0;
					//long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old + Two_ie;
					//long TotOffsetOld = iz_PerZ_Old_p_ix_PerX_Old + (ieStOld << 1);
					long long TotOffsetOld = iz_PerZ_Old_p_ix_PerX_Old + (ieStOld << 1);

					if(TreatPolCompX)
					{
						float *pExSt_Old = OldRadAccessData.pBaseRadX + TotOffsetOld;

						GetCellDataForInterpol1D(pExSt_Old, 2, AuxF);
						SetupCellDataI1D(AuxF, AuxFI);
						UseLowOrderInterp_PolCompX = CheckForLowOrderInterp1D(AuxF, AuxFI, iecOld_mi_ieStOld, &InterpolAux01, InterpolAux02, InterpolAux02I);
						if(!UseLowOrderInterp_PolCompX)
						{
							for(int i=0; i<2; i++) 
							{
								SetupInterpolAux02_1D(AuxF + i, &InterpolAux01, InterpolAux02 + i);
							}
							SetupInterpolAux02_1D(AuxFI, &InterpolAux01, InterpolAux02I);
						}
					}
					if(TreatPolCompZ)
					{
						float *pEzSt_Old = OldRadAccessData.pBaseRadZ + TotOffsetOld;

						GetCellDataForInterpol1D(pEzSt_Old, 2, AuxF + 2);
						SetupCellDataI1D(AuxF + 2, AuxFI + 1);
						UseLowOrderInterp_PolCompZ = CheckForLowOrderInterp1D(AuxF + 2, AuxFI + 1, iecOld_mi_ieStOld, &InterpolAux01, InterpolAux02 + 2, InterpolAux02I + 1);
						if(!UseLowOrderInterp_PolCompZ)
						{
							for(int i=0; i<2; i++) 
							{
								SetupInterpolAux02_1D(AuxF + 2 + i, &InterpolAux01, InterpolAux02 + 2 + i);
							}
							SetupInterpolAux02_1D(AuxFI + 1, &InterpolAux01, InterpolAux02I + 1);
						}
					}
					ieStOldPrev = ieStOld;
				}

				if(TreatPolCompX)
				{
					if(UseLowOrderInterp_PolCompX) 
					{
						InterpolF_LowOrder1D(InterpolAux02, eRel, BufF, 0);
						InterpolFI_LowOrder1D(InterpolAux02I, eRel, BufFI, 0);
					}
					else
					{
						InterpolF1D(InterpolAux02, eRel, BufF, 0);
						InterpolFI1D(InterpolAux02I, eRel, BufFI, 0);
					}

					(*BufFI) *= AuxFI->fNorm;
					ImproveReAndIm(BufF, BufFI);

					//if(FieldShouldBeZeroed)
					//{
					//	*BufF = 0.; *(BufF+1) = 0.;
					//}

					*pEX_New = *BufF;
					*(pEX_New + 1) = *(BufF + 1);
				}
				if(TreatPolCompZ)
				{
					if(UseLowOrderInterp_PolCompZ) 
					{
						InterpolF_LowOrder1D(InterpolAux02, eRel, BufF, 2);
						InterpolFI_LowOrder1D(InterpolAux02I, eRel, BufFI, 1);
					}
					else
					{
						InterpolF1D(InterpolAux02, eRel, BufF, 2);
						InterpolFI1D(InterpolAux02I, eRel, BufFI, 1);
					}

					*(BufFI + 1) *= (AuxFI + 1)->fNorm;
					ImproveReAndIm(BufF + 2, BufFI + 1);

					//if(FieldShouldBeZeroed)
					//{
					//	*(BufF+2) = 0.; *(BufF+3) = 0.;
					//}

					*pEZ_New = *(BufF + 2);
					*(pEZ_New + 1) = *(BufF + 3);
				}
			}
		}
	}
	if(WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(NewRadAccessData, 'a', PolComp);

	OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	OldRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;
	NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	NewRadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;
	return 0;
}

//*************************************************************************

int srTGenOptElem::ReInterpolateWfrSliceSingleE(srTSRWRadStructAccessData& oldRadSingleE, srTSRWRadStructAccessData& newRadMultiE, int ie)
{//similar to "RadResizeCore"; eventually used for propagation at different photon energies
	const double DistAbsTol = 1.E-10;
	bool TreatPolCompX=true, TreatPolCompZ=true;

	int ixStart = 0; //int(NewRadAccessData.AuxLong1);
	int ixEnd = newRadMultiE.nx - 1; //int(NewRadAccessData.AuxLong2);
	int izStart = 0; //int(NewRadAccessData.AuxLong3);
	int izEnd = newRadMultiE.nz - 1; //int(NewRadAccessData.AuxLong4);

	bool WaveFrontTermWasTreated = 0;
	bool OrigWfrQuadTermCanBeTreatedAtResizeX = oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeX;
	bool OrigWfrQuadTermCanBeTreatedAtResizeZ = oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeZ;

	//if(WaveFrontTermCanBeTreated(oldRadSingleE))
	if(WaveFrontTermCanBeTreated(oldRadSingleE, true)) //OC251214
	{
		newRadMultiE.WfrQuadTermCanBeTreatedAtResizeX = oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeX;
		newRadMultiE.WfrQuadTermCanBeTreatedAtResizeZ = oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeZ;
		TreatStronglyOscillatingTerm(oldRadSingleE, 'r', 0);
		WaveFrontTermWasTreated = true;
	}

	double xStepInvOld = 1./oldRadSingleE.xStep;
	double zStepInvOld = 1./oldRadSingleE.zStep;
	int nx_mi_1Old = oldRadSingleE.nx - 1;
	int nz_mi_1Old = oldRadSingleE.nz - 1;
	int nx_mi_2Old = nx_mi_1Old - 1;
	int nz_mi_2Old = nz_mi_1Old - 1;

	srTInterpolAux01 InterpolAux01;
	srTInterpolAux02 InterpolAux02[4], InterpolAux02I[2];
	srTInterpolAuxF AuxF[4], AuxFI[2];
	int ixStOld, izStOld, ixStOldPrev = -1000, izStOldPrev = -1000;

	float *pEX0_New = 0, *pEZ0_New = 0;
	if(TreatPolCompX) pEX0_New = newRadMultiE.pBaseRadX;
	if(TreatPolCompZ) pEZ0_New = newRadMultiE.pBaseRadZ;

	//long PerX_New = newRadMultiE.ne << 1;
	//long PerZ_New = PerX_New*newRadMultiE.nx;
	long long PerX_New = newRadMultiE.ne << 1;
	long long PerZ_New = PerX_New*newRadMultiE.nx;

	//long PerX_Old = 2; //PerX_New;
	//long PerZ_Old = PerX_Old*oldRadSingleE.nx;
	long long PerX_Old = 2; //PerX_New;
	long long PerZ_Old = PerX_Old*oldRadSingleE.nx;

	float BufF[4], BufFI[2];
	int UseLowOrderInterp_PolCompX, UseLowOrderInterp_PolCompZ;
	int result = 0;

	//for(int ie=0; ie<NewRadAccessData.ne; ie++)
	//{
	//ixStOldPrev = -1000; izStOldPrev = -1000;

	//long Two_ie = ie << 1;
	long long Two_ie = ie << 1;
	for(int iz=izStart; iz<=izEnd; iz++)
	{
		if(result = srYield.Check()) return result;

		double zAbs = newRadMultiE.zStart + iz*newRadMultiE.zStep;
		char FieldShouldBeZeroedDueToZ = 0;
		if(newRadMultiE.WfrEdgeCorrShouldBeDone)
		{
			if((zAbs < newRadMultiE.zWfrMin - DistAbsTol) || (zAbs > newRadMultiE.zWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToZ = 1;
		}
		int izcOld = int((zAbs - oldRadSingleE.zStart)*zStepInvOld + 1.E-06);
		if((izcOld < 0) || (izcOld > nz_mi_1Old))
		{
			//set El. field to 0 for all ix
			FieldShouldBeZeroedDueToZ = 1;
		}

		double zRel = zAbs - (oldRadSingleE.zStart + izcOld*oldRadSingleE.zStep);

		if(izcOld == nz_mi_1Old) { izStOld = izcOld - 3; zRel += 2.*oldRadSingleE.zStep;}
		else if(izcOld == nz_mi_2Old) { izStOld = izcOld - 2; zRel += oldRadSingleE.zStep;}
		else if(izcOld == 0) { izStOld = izcOld; zRel -= oldRadSingleE.zStep;}
		else izStOld = izcOld - 1;

		zRel *= zStepInvOld;
		int izcOld_mi_izStOld = izcOld - izStOld;
		//long izPerZ_New = iz*PerZ_New;
		long long izPerZ_New = iz*PerZ_New;

		float *pEX_StartForX_New = 0, *pEZ_StartForX_New = 0;
		if(TreatPolCompX) pEX_StartForX_New = pEX0_New + izPerZ_New;
		if(TreatPolCompZ) pEZ_StartForX_New = pEZ0_New + izPerZ_New;

		for(int ix=ixStart; ix<=ixEnd; ix++)
		{
			//long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
			long long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
			float *pEX_New = 0, *pEZ_New = 0;
			if(TreatPolCompX) pEX_New = pEX_StartForX_New + ixPerX_New_p_Two_ie;
			if(TreatPolCompZ) pEZ_New = pEZ_StartForX_New + ixPerX_New_p_Two_ie;

			double xAbs = newRadMultiE.xStart + ix*newRadMultiE.xStep;
			char FieldShouldBeZeroedDueToX = 0;
			if(newRadMultiE.WfrEdgeCorrShouldBeDone)
			{
				if((xAbs < newRadMultiE.xWfrMin - DistAbsTol) || (xAbs > newRadMultiE.xWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToX = 1;
			}

			int ixcOld = int((xAbs - oldRadSingleE.xStart)*xStepInvOld + 1.E-06);
			if((ixcOld < 0) || (ixcOld > nx_mi_1Old))
			{
				FieldShouldBeZeroedDueToX = 1;
			}
			char FieldShouldBeZeroed = (FieldShouldBeZeroedDueToX || FieldShouldBeZeroedDueToZ);

			if(FieldShouldBeZeroed)
			{
				//*BufF = 0.; *(BufF+1) = 0.;
				if(TreatPolCompX)
				{
					*pEX_New = 0.;
					*(pEX_New+1) = 0.;
				}
				if(TreatPolCompZ)
				{
					*pEZ_New = 0.;
					*(pEZ_New+1) = 0.;
				}
				continue;
			}

			double xRel = xAbs - (oldRadSingleE.xStart + ixcOld*oldRadSingleE.xStep);

			if(ixcOld == nx_mi_1Old) { ixStOld = ixcOld - 3; xRel += 2.*oldRadSingleE.xStep;}
			else if(ixcOld == nx_mi_2Old) { ixStOld = ixcOld - 2; xRel += oldRadSingleE.xStep;}
			else if(ixcOld == 0) { ixStOld = ixcOld; xRel -= oldRadSingleE.xStep;}
			else ixStOld = ixcOld - 1;

			xRel *= xStepInvOld;
			int ixcOld_mi_ixStOld = ixcOld - ixStOld;

			if((izStOld != izStOldPrev) || (ixStOld != ixStOldPrev))
			{
				UseLowOrderInterp_PolCompX = 0, UseLowOrderInterp_PolCompZ = 0;
				//long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old + Two_ie;
				//long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old; //old is single slice
				long long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old; //old is single slice

				if(TreatPolCompX)
				{
					float* pExSt_Old = oldRadSingleE.pBaseRadX + TotOffsetOld;
					GetCellDataForInterpol(pExSt_Old, PerX_Old, PerZ_Old, AuxF);
					SetupCellDataI(AuxF, AuxFI);
					UseLowOrderInterp_PolCompX = CheckForLowOrderInterp(AuxF, AuxFI, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02, InterpolAux02I);

					if(!UseLowOrderInterp_PolCompX)
					{
						for(int i=0; i<2; i++) 
						{
							SetupInterpolAux02(AuxF + i, &InterpolAux01, InterpolAux02 + i);
						}
						SetupInterpolAux02(AuxFI, &InterpolAux01, InterpolAux02I);
					}
				}
				if(TreatPolCompZ)
				{
					float* pEzSt_Old = oldRadSingleE.pBaseRadZ + TotOffsetOld;
					GetCellDataForInterpol(pEzSt_Old, PerX_Old, PerZ_Old, AuxF+2);
					SetupCellDataI(AuxF+2, AuxFI+1);
					UseLowOrderInterp_PolCompZ = CheckForLowOrderInterp(AuxF+2, AuxFI+1, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02+2, InterpolAux02I+1);

					if(!UseLowOrderInterp_PolCompZ)
					{
						for(int i=0; i<2; i++) 
						{
							SetupInterpolAux02(AuxF+2+i, &InterpolAux01, InterpolAux02+2+i);
						}
						SetupInterpolAux02(AuxFI+1, &InterpolAux01, InterpolAux02I+1);
					}
				}
				ixStOldPrev = ixStOld; izStOldPrev = izStOld;
			}

			if(TreatPolCompX)
			{
				if(UseLowOrderInterp_PolCompX) 
				{
					InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 0);
					InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 0);
				}
				else
				{
					InterpolF(InterpolAux02, xRel, zRel, BufF, 0);
					InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 0);
				}

				(*BufFI) *= AuxFI->fNorm;
				ImproveReAndIm(BufF, BufFI);

				//if(FieldShouldBeZeroed)
				//{
				//	*BufF = 0.; *(BufF+1) = 0.;
				//}

				*pEX_New = *BufF;
				*(pEX_New+1) = *(BufF+1);
			}
			if(TreatPolCompZ)
			{
				if(UseLowOrderInterp_PolCompZ) 
				{
					InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 2);
					InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 1);
				}
				else
				{
					InterpolF(InterpolAux02, xRel, zRel, BufF, 2);
					InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 1);
				}

				(*(BufFI+1)) *= (AuxFI+1)->fNorm;
				ImproveReAndIm(BufF+2, BufFI+1);

				//if(FieldShouldBeZeroed)
				//{
				//	*(BufF+2) = 0.; *(BufF+3) = 0.;
				//}

				*pEZ_New = *(BufF+2);
				*(pEZ_New+1) = *(BufF+3);
			}
		}
	}
	//}
	if(WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(newRadMultiE, 'a', 0, ie);

	oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	oldRadSingleE.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;
	newRadMultiE.WfrQuadTermCanBeTreatedAtResizeX = OrigWfrQuadTermCanBeTreatedAtResizeX;
	newRadMultiE.WfrQuadTermCanBeTreatedAtResizeZ = OrigWfrQuadTermCanBeTreatedAtResizeZ;
	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeCore_OnlyLargerRange(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp)
{
	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x'));
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z'));

	float *pEX0_New = NewRadAccessData.pBaseRadX;
	float *pEZ0_New = NewRadAccessData.pBaseRadZ;

	float *pEX0_Old = OldRadAccessData.pBaseRadX;
	float *pEZ0_Old = OldRadAccessData.pBaseRadZ;

	//long PerX_New = NewRadAccessData.ne << 1;
	//long PerZ_New = PerX_New*NewRadAccessData.nx;
	long long PerX_New = NewRadAccessData.ne << 1;
	long long PerZ_New = PerX_New*NewRadAccessData.nx;

	//long PerX_Old = PerX_New;
	//long PerZ_Old = PerX_Old*OldRadAccessData.nx;
	long long PerX_Old = PerX_New;
	long long PerZ_Old = PerX_Old*OldRadAccessData.nx;

	int ixStart = int(NewRadAccessData.AuxLong1);
	int ixEnd = int(NewRadAccessData.AuxLong2);
	int izStart = int(NewRadAccessData.AuxLong3);
	int izEnd = int(NewRadAccessData.AuxLong4);

	double xStepInvOld = 1./OldRadAccessData.xStep;
	double zStepInvOld = 1./OldRadAccessData.zStep;

	for(long ie=0; ie<NewRadAccessData.ne; ie++)
	{
		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;
		for(long iz=izStart; iz<=izEnd; iz++)
		{
			//long izPerZ_New = iz*PerZ_New;
			long long izPerZ_New = iz*PerZ_New;
			float *pEX_StartForX_New = pEX0_New + izPerZ_New;
			float *pEZ_StartForX_New = pEZ0_New + izPerZ_New;

			//long izPerZ_Old = (iz - izStart)*PerZ_Old;

			double zAbs = NewRadAccessData.zStart + iz*NewRadAccessData.zStep;
			long izOld = long((zAbs - OldRadAccessData.zStart)*zStepInvOld + 1.E-08);
			//long izPerZ_Old = izOld*PerZ_Old;
			long long izPerZ_Old = izOld*PerZ_Old;

			float *pEX_StartForX_Old = pEX0_Old + izPerZ_Old;
			float *pEZ_StartForX_Old = pEZ0_Old + izPerZ_Old;

			for(long ix=ixStart; ix<=ixEnd; ix++)
			{
				//long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
				long long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
				float *pEX_New = pEX_StartForX_New + ixPerX_New_p_Two_ie;
				float *pEZ_New = pEZ_StartForX_New + ixPerX_New_p_Two_ie;

				//long ixPerX_Old_p_Two_ie = (ix - ixStart)*PerX_Old + Two_ie;

				double xAbs = NewRadAccessData.xStart + ix*NewRadAccessData.xStep;
				long ixOld = long((xAbs - OldRadAccessData.xStart)*xStepInvOld + 1.E-08);
				//long ixPerX_Old_p_Two_ie = ixOld*PerX_Old + Two_ie;
				long long ixPerX_Old_p_Two_ie = ixOld*PerX_Old + Two_ie;

				float *pEX_Old = pEX_StartForX_Old + ixPerX_Old_p_Two_ie;
				float *pEZ_Old = pEZ_StartForX_Old + ixPerX_Old_p_Two_ie;

				if(TreatPolCompX) { *pEX_New = *pEX_Old; *(pEX_New + 1) = *(pEX_Old + 1);}
				if(TreatPolCompZ) { *pEZ_New = *pEZ_Old; *(pEZ_New + 1) = *(pEZ_Old + 1);}
			}
		}
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeCore_OnlyLargerRangeE(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp)
{
	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x')) && (OldRadAccessData.pBaseRadX != 0);
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z')) && (OldRadAccessData.pBaseRadZ != 0);

	float *pEX0_New = NewRadAccessData.pBaseRadX;
	float *pEZ0_New = NewRadAccessData.pBaseRadZ;

	float *pEX0_Old = OldRadAccessData.pBaseRadX;
	float *pEZ0_Old = OldRadAccessData.pBaseRadZ;

	//long PerX_New = NewRadAccessData.ne << 1;
	//long PerZ_New = PerX_New*NewRadAccessData.nx;
	long long PerX_New = NewRadAccessData.ne << 1;
	long long PerZ_New = PerX_New*NewRadAccessData.nx;

	//long PerX_Old = OldRadAccessData.ne << 1;
	//long PerZ_Old = PerX_Old*OldRadAccessData.nx;
	long long PerX_Old = OldRadAccessData.ne << 1;
	long long PerZ_Old = PerX_Old*OldRadAccessData.nx;

	int ieStart = int(NewRadAccessData.AuxLong1);
	int ieEnd = int(NewRadAccessData.AuxLong2);

	double eStepInvOld = 1./OldRadAccessData.eStep;
	
	for(long iz=0; iz<NewRadAccessData.nz; iz++)
	{
		//long iz_PerZ_New = iz*PerZ_New;
		//long iz_PerZ_Old = iz*PerZ_Old;
		long long iz_PerZ_New = iz*PerZ_New;
		long long iz_PerZ_Old = iz*PerZ_Old;

		for(long ix=0; ix<NewRadAccessData.nx; ix++)
		{
			//long iz_PerZ_New_p_ix_PerX_New = iz_PerZ_New + ix*PerX_New;
			//long iz_PerZ_Old_p_ix_PerX_Old = iz_PerZ_Old + ix*PerX_Old;
			long long iz_PerZ_New_p_ix_PerX_New = iz_PerZ_New + ix*PerX_New;
			long long iz_PerZ_Old_p_ix_PerX_Old = iz_PerZ_Old + ix*PerX_Old;

			for(long ie=ieStart; ie<=ieEnd; ie++)
			{
				//long ofstNew = iz_PerZ_New_p_ix_PerX_New + (ie << 1);
				long long ofstNew = iz_PerZ_New_p_ix_PerX_New + (ie << 1);
				float *pEX_New = pEX0_New + ofstNew;
				float *pEZ_New = pEZ0_New + ofstNew;

				double eAbs = NewRadAccessData.eStart + ie*NewRadAccessData.eStep;
				long ieOld = long((eAbs - OldRadAccessData.eStart)*eStepInvOld + 1.E-08);

				//long ofstOld = iz_PerZ_Old_p_ix_PerX_Old + (ieOld << 1);
				long long ofstOld = iz_PerZ_Old_p_ix_PerX_Old + (ieOld << 1);
				float *pEX_Old = pEX0_Old + ofstOld;
				float *pEZ_Old = pEZ0_Old + ofstOld;

				if(TreatPolCompX) { *pEX_New = *pEX_Old; *(pEX_New + 1) = *(pEX_Old + 1);}
				if(TreatPolCompZ) { *pEZ_New = *pEZ_Old; *(pEZ_New + 1) = *(pEZ_Old + 1);}
			}
		}
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeGen1D(srTRadSect1D& RadSect1D, srTRadResize1D& RadResizeStruct1D)
{
  	if((RadResizeStruct1D.pm == 1.) && (RadResizeStruct1D.pd == 1.)) return 0;
	int result;

	double pTot = RadResizeStruct1D.pm*RadResizeStruct1D.pd;
	char RadShouldBeChanged = (pTot != 1.);
	long NewN = RadSect1D.np; 
	CGenMathFFT1D FFT;
	if(RadShouldBeChanged)
	{
		double pTotNewN = pTot*NewN;
		long Long_pTotNewN = long(pTotNewN);
		NewN = ((pTotNewN - Long_pTotNewN) >= 0.5)? (Long_pTotNewN + 1) : Long_pTotNewN;
		FFT.NextCorrectNumberForFFT(NewN);
	}

	if((NewN == RadSect1D.np) && (RadResizeStruct1D.pd == 1.)) return 0; //OC fixed here: October 2003

	if(RadResizeStruct1D.UseOtherSideFFT)
	{
		char ToRepres = (RadSect1D.Pres == 0)? 1 : 0;

		if(result = SetRadRepres1D(&RadSect1D, ToRepres)) return result;

		double pmNew = RadResizeStruct1D.pd, pdNew = RadResizeStruct1D.pm;
		RadResizeStruct1D.pm = pmNew; RadResizeStruct1D.pd = pdNew;
	}

	srTRadSect1D NewRadSect1D = RadSect1D;
	NewRadSect1D.DeleteArraysAtDestruction = 0;
	
	if(RadShouldBeChanged)
	{
		NewRadSect1D.np = NewN;
	}

	double MagOld = (RadSect1D.np - 1)*RadSect1D.ArgStep;
	if(RadResizeStruct1D.pm > 1.)
	{
		long NInner = (RadResizeStruct1D.pd != 1.)? long(NewN/RadResizeStruct1D.pm) : RadSect1D.np;

		long NewN_mi_NInner = NewN - NInner;
		long NOuterLeft = NewN_mi_NInner >> 1;

		if((NOuterLeft << 1) != NewN_mi_NInner) NOuterLeft++;

		NewRadSect1D.ArgStep = MagOld/(NInner - 1);
		NewRadSect1D.ArgStart = RadSect1D.ArgStart - NOuterLeft*NewRadSect1D.ArgStep;
		NewRadSect1D.AuxLong1 = NOuterLeft; // iStart
		NewRadSect1D.AuxLong2 = NOuterLeft + NInner - 1; // iEnd
		NewRadSect1D.PreserveLogicsOfWfrLimitsAtRangeResizing(&RadSect1D);
	}
	else
	{//bug fixed here: May 2003

		double EndOld = RadSect1D.ArgStart + (RadSect1D.np - 1)*RadSect1D.ArgStep;
		//double Mid = RadSect1D.ArgStart + (RadSect1D.np >> 1)*RadSect1D.ArgStep;
		double Mid = 0.5*(RadSect1D.ArgStart + EndOld); //OC
		double MagNew = MagOld*RadResizeStruct1D.pm;

		if(RadResizeStruct1D.pd != 1.)
		{
			NewRadSect1D.ArgStep = MagNew/(NewN - 1);
		}
		else
		{
			NewRadSect1D.ArgStep = RadSect1D.ArgStep;
		}
		//NewRadSect1D.ArgStart = Mid - (NewN >> 1)*NewRadSect1D.ArgStep;
		NewRadSect1D.ArgStart = Mid - 0.5*MagNew; //OC
		NewRadSect1D.AuxLong1 = 0; // iStart
		NewRadSect1D.AuxLong2 = NewN - 1; // iEnd

		if(NewRadSect1D.ArgStart > NewRadSect1D.WfrMin) NewRadSect1D.WfrMin = NewRadSect1D.ArgStart;
		if((NewRadSect1D.ArgStart + MagNew) < NewRadSect1D.WfrMax) NewRadSect1D.WfrMax = NewRadSect1D.ArgStart + MagNew;
	}

	//long TotAmOfOldData = RadSect1D.np << 1;
	long long TotAmOfOldData = RadSect1D.np << 1;
	float *OldRadXCopy = new float[TotAmOfOldData];
	if(OldRadXCopy == 0) return MEMORY_ALLOCATION_FAILURE;
	float *OldRadZCopy = new float[TotAmOfOldData];
	if(OldRadZCopy == 0) return MEMORY_ALLOCATION_FAILURE;

	float *tOldRadXCopy = OldRadXCopy, *tOldRadZCopy = OldRadZCopy;
	float *tBaseRadX = RadSect1D.pEx, *tBaseRadZ = RadSect1D.pEz;
	//for(long i=0; i<TotAmOfOldData; i++) 
	for(long long i=0; i<TotAmOfOldData; i++) 
	{
		*(tOldRadXCopy++) = *(tBaseRadX++);
		*(tOldRadZCopy++) = *(tBaseRadZ++);
	}

	if(RadShouldBeChanged)
	{
		//long Two_np = NewRadSect1D.np << 1;
		long long Two_np = NewRadSect1D.np << 1;
		NewRadSect1D.pEx = new float[Two_np];
		if(NewRadSect1D.pEx == 0) return MEMORY_ALLOCATION_FAILURE;
		NewRadSect1D.pEz = new float[Two_np];
		if(NewRadSect1D.pEz == 0) return MEMORY_ALLOCATION_FAILURE;
	}

	//long TotAmOfNewData = NewN << 1;
	long long TotAmOfNewData = NewN << 1;
	tBaseRadX = NewRadSect1D.pEx;
	tBaseRadZ = NewRadSect1D.pEz;
	//for(long j=0; j<TotAmOfNewData; j++)
	for(long long j=0; j<TotAmOfNewData; j++)
	{
		*(tBaseRadX++) = 0.; *(tBaseRadZ++) = 0.; 
	}

	if(RadShouldBeChanged)
	{
		if(RadSect1D.DeleteArraysAtDestruction) RadSect1D.DeleteArrays(); //MEM LEAK REPAIR
	}

	RadSect1D.pEx = OldRadXCopy;
	RadSect1D.pEz = OldRadZCopy;

	if(result = RadResizeCore1D(RadSect1D, NewRadSect1D, RadResizeStruct1D)) return result;

	if(OldRadXCopy != 0) delete[] OldRadXCopy;
	if(OldRadZCopy != 0) delete[] OldRadZCopy;

	char OldDeleteArraysAtDestruction = RadSect1D.DeleteArraysAtDestruction;

	RadSect1D = NewRadSect1D;

	RadSect1D.DeleteArraysAtDestruction = OldDeleteArraysAtDestruction;
	NewRadSect1D.DeleteArraysAtDestruction = 0;

	if(RadResizeStruct1D.UseOtherSideFFT)
	{
		char ToRepres = (RadSect1D.Pres == 0)? 1 : 0;
		if(result = SetRadRepres1D(&RadSect1D, ToRepres)) return result;

		double pmNew = RadResizeStruct1D.pd, pdNew = RadResizeStruct1D.pm;
		RadResizeStruct1D.pm = pmNew; RadResizeStruct1D.pd = pdNew;
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::RadResizeCore1D(srTRadSect1D& OldRadSect1D, srTRadSect1D& NewRadSect1D, srTRadResize1D& RadResizeStruct1D)
{
	const double DistAbsTol = 1.E-10;
	const double RelStepTol = 1.e-05; // To steer

	int iStart = int(NewRadSect1D.AuxLong1);
	int iEnd = int(NewRadSect1D.AuxLong2);

	char OnlyMakeLargerRange = ((RadResizeStruct1D.pd == 1.) && (RadResizeStruct1D.pm > 1.))
							&& (::fabs(OldRadSect1D.ArgStep - NewRadSect1D.ArgStep) < RelStepTol);

	char WaveFrontTermWasTreated = 0;
	if((!RadResizeStruct1D.DoNotTreatSpherTerm) && WaveFrontTermCanBeTreated1D(OldRadSect1D))
	{
		if(!OnlyMakeLargerRange)
		{
			TreatStronglyOscillatingTerm1D(OldRadSect1D, 'r');
			WaveFrontTermWasTreated = 1;
		}
	}

	double StepInvOld = 1./OldRadSect1D.ArgStep;
	int np_mi_1Old = OldRadSect1D.np - 1;
	int np_mi_2Old = np_mi_1Old - 1;

	double argStartOld = OldRadSect1D.ArgStart; //OC091207
	double argEndOld = OldRadSect1D.ArgStart + np_mi_1Old*OldRadSect1D.ArgStep;

	srTInterpolAux01_1D InterpolAux01;
	srTInterpolAux02_1D InterpolAux02[4], InterpolAux02I[2];
	srTInterpolAuxF_1D AuxF[4], AuxFI[2];

	int iStOld, iStOldPrev = -1000;
	float *pEX0_New = NewRadSect1D.pEx, *pEZ0_New = NewRadSect1D.pEz;

	//long Per_New = 2;
	//long Per_Old = 2;
	long long Per_New = 2;
	long long Per_Old = 2;

	float BufF[4], BufFI[2];
	char UseLowOrderInterp_PolCompX, UseLowOrderInterp_PolCompZ;

	for(int i=iStart; i<=iEnd; i++)
	{
		//long iPer_New = i*Per_New;
		long long iPer_New = i*Per_New;
		float *pEX_New = pEX0_New + iPer_New;
		float *pEZ_New = pEZ0_New + iPer_New;

		if(OnlyMakeLargerRange)
		{
			//long TotOffsetOld = (i - iStart)*Per_Old;
			long long TotOffsetOld = (i - iStart)*Per_Old;
			float* pExSt_Old = OldRadSect1D.pEx + TotOffsetOld;
			float* pEzSt_Old = OldRadSect1D.pEz + TotOffsetOld;

			*pEX_New = *pExSt_Old; *(pEX_New + 1) = *(pExSt_Old + 1);
			*pEZ_New = *pEzSt_Old; *(pEZ_New + 1) = *(pEzSt_Old + 1);
			continue;
		}

		double ArgAbs = NewRadSect1D.ArgStart + i*NewRadSect1D.ArgStep;
		if((ArgAbs < argStartOld) || (ArgAbs > argEndOld)) //OC091207
		{
			*pEX_New = 0.; *(pEX_New + 1) = 0.;
			*pEZ_New = 0.; *(pEZ_New + 1) = 0.;
			continue;
		}

		char FieldShouldBeZeroed = 0;
		if(NewRadSect1D.WfrEdgeCorrShouldBeDone)
		{
			if((ArgAbs < NewRadSect1D.WfrMin - DistAbsTol) || (ArgAbs > NewRadSect1D.WfrMax + DistAbsTol)) FieldShouldBeZeroed = 1;
		}

		int icOld = int((ArgAbs - OldRadSect1D.ArgStart)*StepInvOld + 1.E-06);
		double ArgRel = ArgAbs - (OldRadSect1D.ArgStart + icOld*OldRadSect1D.ArgStep);

		if(icOld == np_mi_1Old) { iStOld = icOld - 3; ArgRel += 2.*OldRadSect1D.ArgStep;}
		else if(icOld == np_mi_2Old) { iStOld = icOld - 2; ArgRel += OldRadSect1D.ArgStep;}
		else if(icOld == 0) { iStOld = icOld; ArgRel -= OldRadSect1D.ArgStep;}
		else iStOld = icOld - 1;

		ArgRel *= StepInvOld;
		int icOld_mi_iStOld = icOld - iStOld;

		if(iStOld != iStOldPrev)
		{
			UseLowOrderInterp_PolCompX = 0, UseLowOrderInterp_PolCompZ = 0;

			//long TotOffsetOld = iStOld*Per_Old;
			long long TotOffsetOld = iStOld*Per_Old;

			float* pExSt_Old = OldRadSect1D.pEx + TotOffsetOld;
			GetCellDataForInterpol1D(pExSt_Old, Per_Old, AuxF);
			SetupCellDataI1D(AuxF, AuxFI);

			UseLowOrderInterp_PolCompX = CheckForLowOrderInterp1D(AuxF, AuxFI, icOld_mi_iStOld, &InterpolAux01, InterpolAux02, InterpolAux02I);
			if(!UseLowOrderInterp_PolCompX)
			{
				for(int i=0; i<2; i++) 
				{
					SetupInterpolAux02_1D(AuxF + i, &InterpolAux01, InterpolAux02 + i);
				}
				SetupInterpolAux02_1D(AuxFI, &InterpolAux01, InterpolAux02I);
			}

			float* pEzSt_Old = OldRadSect1D.pEz + TotOffsetOld;
			GetCellDataForInterpol1D(pEzSt_Old, Per_Old, AuxF+2);
			SetupCellDataI1D(AuxF+2, AuxFI+1);

			UseLowOrderInterp_PolCompZ = CheckForLowOrderInterp1D(AuxF+2, AuxFI+1, icOld_mi_iStOld, &InterpolAux01, InterpolAux02+2, InterpolAux02I+1);
			if(!UseLowOrderInterp_PolCompZ)
			{
				for(int i=0; i<2; i++) 
				{
					SetupInterpolAux02_1D(AuxF+2+i, &InterpolAux01, InterpolAux02+2+i);
				}
				SetupInterpolAux02_1D(AuxFI+1, &InterpolAux01, InterpolAux02I+1);
			}

			iStOldPrev = iStOld;
		}

		if(UseLowOrderInterp_PolCompX) 
		{
			InterpolF_LowOrder1D(InterpolAux02, ArgRel, BufF, 0);
			InterpolFI_LowOrder1D(InterpolAux02I, ArgRel, BufFI, 0);
		}
		else
		{
			InterpolF1D(InterpolAux02, ArgRel, BufF, 0);
			InterpolFI1D(InterpolAux02I, ArgRel, BufFI, 0);
		}

		(*BufFI) *= AuxFI->fNorm;
		ImproveReAndIm(BufF, BufFI);
		if(FieldShouldBeZeroed)
		{
			*BufF = 0.; *(BufF+1) = 0.;
		}
		*pEX_New = *BufF;
		*(pEX_New+1) = *(BufF+1);

		if(UseLowOrderInterp_PolCompZ) 
		{
			InterpolF_LowOrder1D(InterpolAux02, ArgRel, BufF, 2);
			InterpolFI_LowOrder1D(InterpolAux02I, ArgRel, BufFI, 1);
		}
		else
		{
			InterpolF1D(InterpolAux02, ArgRel, BufF, 2);
			InterpolFI1D(InterpolAux02I, ArgRel, BufFI, 1);
		}

		(*(BufFI+1)) *= (AuxFI+1)->fNorm;
		ImproveReAndIm(BufF+2, BufFI+1);
		if(FieldShouldBeZeroed)
		{
			*(BufF+2) = 0.; *(BufF+3) = 0.;
		}
		*pEZ_New = *(BufF+2);
		*(pEZ_New+1) = *(BufF+3);
	}
	if(WaveFrontTermWasTreated) TreatStronglyOscillatingTerm1D(NewRadSect1D, 'a');

	return 0;
}

//*************************************************************************

char srTGenOptElem::WaveFrontTermCanBeTreated(srTSRWRadStructAccessData& RadAccessData, bool checkBenefit)
{
	//Later treat X and Z fully separately here and at removing the corresponding terms from Phase !!!

	const double CritRatTransvLong = 0.1;
	//const double CritRelRobsErr = 0.2; //0.1; //0.2;
	const double CritRelRobsErr = 0.4; //0.1; //0.2; //OC17032016 ?
	const double Pi = 3.14159265358979;

	bool RobsXErrIsSmall = ::fabs(RadAccessData.RobsXAbsErr) < CritRelRobsErr*(::fabs(RadAccessData.RobsX));
	bool RobsZErrIsSmall = ::fabs(RadAccessData.RobsZAbsErr) < CritRelRobsErr*(::fabs(RadAccessData.RobsZ));

	if(RadAccessData.Pres == 0) // Coord
	{
		/**
		double xMagn = ::fabs((RadAccessData.nx - 1)*RadAccessData.xStep);
		double zMagn = ::fabs((RadAccessData.nz - 1)*RadAccessData.zStep);

		char AnglesXAreSmall = (xMagn < CritRatTransvLong*(::fabs(RadAccessData.RobsX)));
		char AnglesZAreSmall = (zMagn < CritRatTransvLong*(::fabs(RadAccessData.RobsZ)));
		**/

		//RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = (AnglesXAreSmall && RobsXErrIsSmall);
		//RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = (AnglesZAreSmall && RobsZErrIsSmall);
		//OC260114
		//Removed requirement of small angles, because this was producing a buggy behavior for grazing-incidence mirrors
		//Actually, the quadratic term can be subtracted and added even if angles are large(?) 
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = RobsXErrIsSmall; 
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = RobsZErrIsSmall;

		//return (RadAccessData.WfrQuadTermCanBeTreatedAtResizeX || RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ);
	}
	else // Ang
	{// Not sure about this...
		/**
		double ePh = RadAccessData.eStart;
		double xRatMax = 1.E-23, zRatMax = 1.E-23, MaxPhaseChange = 1.E-23;
		for(int ie=0; ie<RadAccessData.ne; ie++)
		{
			//double Lambda_m = 1.239854e-06/ePh;
			double Lambda_m = 1.239842e-06/ePh;
			double xTetMagn = ::fabs((RadAccessData.nx - 1)*RadAccessData.xStep*Lambda_m);
			double zTetMagn = ::fabs((RadAccessData.nz - 1)*RadAccessData.zStep*Lambda_m);
			double PhaseChange = ::fabs((Pi/Lambda_m)*(RadAccessData.RobsX*xTetMagn*xTetMagn + RadAccessData.RobsZ*zTetMagn*zTetMagn));

			if(xTetMagn > xRatMax) xRatMax = xTetMagn;
			if(zTetMagn > zRatMax) zRatMax = zTetMagn;
			if(PhaseChange > MaxPhaseChange) MaxPhaseChange = PhaseChange;

			ePh += RadAccessData.eStep;
		}

		char AnglesXAreSmall = (xRatMax < CritRatTransvLong);
		char AnglesZAreSmall = (zRatMax < CritRatTransvLong);
		char PhaseChangeIsLarge = (MaxPhaseChange > 2.*Pi);
		**/

		//RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = (AnglesXAreSmall && RobsXErrIsSmall);
		//RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = (AnglesZAreSmall && RobsZErrIsSmall);
		//OC260114
		//Removed requirement of small angles, because this was producing a buggy behavior for grazing-incidence mirrors
		//Actually, the quadratic term can be subtracted and added even angles are large(?) 
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = RobsXErrIsSmall; 
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = RobsZErrIsSmall;

		//return ((RadAccessData.WfrQuadTermCanBeTreatedAtResizeX || RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) && PhaseChangeIsLarge);
		//return (RadAccessData.WfrQuadTermCanBeTreatedAtResizeX || RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ); //OC260114
	}

	//OC290914
	if(checkBenefit)
	{
		if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeX) RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = RadAccessData.CheckIfQuadTermTreatIsBenefit('x');
		if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = RadAccessData.CheckIfQuadTermTreatIsBenefit('z');
	}

	return (RadAccessData.WfrQuadTermCanBeTreatedAtResizeX || RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ);
}

//*************************************************************************

void srTGenOptElem::TreatStronglyOscillatingTerm(srTSRWRadStructAccessData& RadAccessData, char AddOrRem, char PolComp, int ieOnly)
{
	//Later treat X and Z coordinates separately here!!!

	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x')) && (RadAccessData.pBaseRadX != 0); //OC13112011
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z')) && (RadAccessData.pBaseRadZ != 0);

	double Rx = RadAccessData.RobsX;
	double Rz = RadAccessData.RobsZ;

	const double Pi = 3.14159265358979;
	//double Const = Pi*1.E+06/1.239854; // Assumes m and eV
	double Const = Pi*1.E+06/1.23984186; // Assumes m and eV

/**
	//Correcting Effective Wavefront Radius taking into account Rayleigh Length (relies on Stat. Moments...)
	//This was introduced in attempt to improve the quality of resizing (before sample) at the CHX (NSLS-II) CDI experiment simulation;
	//however later was found to introduce problems in simulations with grazing-incidence mirrors (SRWLIB_Example10.py) and with VLS Grating (SRWLIB_Example12.py)
	//therefore "rolled back"
	if(AddOrRem == 'r')
	{
		bool doTuneR = true;
		if((::fabs(RadAccessData.wfrReffX + 1.E+023) < 1.E+18) && (::fabs(RadAccessData.wfrReffZ + 1.E+023) < 1.E+18)) doTuneR = false;

		if(doTuneR)
		{
			int ieCen = 0;
			double photEn0 = RadAccessData.eStart;
			if ((RadAccessData.ne > 1) && (RadAccessData.eStep > 0))
			{
				photEn0 = RadAccessData.avgPhotEn;

				double d_ieCen = (photEn0 - RadAccessData.eStart) / (RadAccessData.eStep);
				ieCen = (int)d_ieCen;
				if ((d_ieCen - ieCen) > 0.5) ieCen++;
				if (ieCen < 0) ieCen = 0;
				else if (ieCen >= RadAccessData.ne) ieCen = RadAccessData.ne - 1;
			}

			srTMomentsPtrs MomX(RadAccessData.pMomX, ieCen), MomZ(RadAccessData.pMomZ, ieCen);
			char TreatExOrEz = (*(MomX.pTotPhot) >= *(MomZ.pTotPhot)) ? 'x' : 'z';
			//double SigX=0, SigZ=0;
			double SigXp = 0, SigZp = 0;
			if (TreatExOrEz == 'x')
			{
				//if((!MomX.precCenMomIsOK) || (MomX.SqrtMxx == 0) || (MomX.SqrtMzz == 0))
				////if((!MomX.precCenMomIsOK) || (MomX.SqrtMxpxp == 0) || (MomX.SqrtMzpzp == 0))
				//{
				//	ComputeRadMoments(&RadAccessData);
				//	MomX.ComputeCentralMoments();
				//}//OC150914: eliminated because ComputeRadMoments calls this function
				//SigX = MomX.SqrtMxx;
				SigXp = MomX.SqrtMxpxp;
				//SigZ = MomX.SqrtMzz;
				SigZp = MomX.SqrtMzpzp;
			}
			else
			{
				//if((!MomZ.precCenMomIsOK) || (MomZ.SqrtMxx == 0) || (MomZ.SqrtMzz == 0))
				////if((!MomZ.precCenMomIsOK) || (MomZ.SqrtMxpxp == 0) || (MomZ.SqrtMzpzp == 0))
				//{//OC13112010: uncommented
				//	ComputeRadMoments(&RadAccessData);
				//	MomZ.ComputeCentralMoments();
				//}//OC150914: eliminated because ComputeRadMoments calls this function
				//SigX = MomZ.SqrtMxx;
				SigXp = MomZ.SqrtMxpxp;
				//SigZ = MomZ.SqrtMzz;
				SigZp = MomZ.SqrtMzpzp;
			}

			const double coefAngRange = 0.35; //0.2; //0.5; //0.1; //to tune
			double apSigXp = coefAngRange*SigXp; //apparent RMS angular divergence
			double apSigZp = coefAngRange*SigZp; //apparent RMS angular divergence

			const double fourPi = 4.*Pi;
			double lambM = 1.23984186E-06 / photEn0;
			double absLenRayleighX = 0, absLenRayleighZ = 0;
			if (apSigXp > 0.)
			{
				absLenRayleighX = lambM / (fourPi*apSigXp*apSigXp);
				double ratRayleighX = absLenRayleighX / RadAccessData.RobsX;
				Rx = RadAccessData.RobsX*(1. + ratRayleighX*ratRayleighX);
			}
			if (apSigZp > 0.)
			{
				absLenRayleighZ = lambM / (fourPi*apSigZp*apSigZp);
				double ratRayleighZ = absLenRayleighZ / RadAccessData.RobsZ;
				Rz = RadAccessData.RobsZ*(1. + ratRayleighZ*ratRayleighZ);
			}
		}
		RadAccessData.wfrReffX = Rx;
		RadAccessData.wfrReffZ = Rz;
	}
	else if(AddOrRem == 'a')
	{
		Rx = RadAccessData.wfrReffX;
		Rz = RadAccessData.wfrReffZ;
	}
**/

	//double RxEst=1.E+23, xcEst=0., RzEst=1.E+23, zcEst=0.;
	//RadAccessData.EstimWfrRadCen(RxEst, xcEst, 'x', 0, 0.1);
	//RadAccessData.EstimWfrRadCen(RzEst, zcEst, 'z', 0, 0.1);

	double ConstRx = (RadAccessData.Pres == 0)? Const/Rx : -Const*Rx;
	double ConstRz = (RadAccessData.Pres == 0)? Const/Rz : -Const*Rz;

	if(AddOrRem == 'r') { ConstRx = -ConstRx; ConstRz = -ConstRz;}

	double ConstRxE, ConstRzE;
	double ePh = RadAccessData.eStart, x, z, zE2;
	double Phase;
	float CosPh, SinPh;

	float *pEX0 = 0, *pEZ0 = 0;
	if(TreatPolCompX) pEX0 = RadAccessData.pBaseRadX;
	if(TreatPolCompZ) pEZ0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//srTFFT2D AuxFFT2D;
	int ieStart=0, ieBefEnd=RadAccessData.ne;
	if((ieOnly >= 0) && (ieOnly < RadAccessData.ne))
	{
		ieStart = ieOnly; ieBefEnd = ieOnly + 1;
	}

	//for(int ie=0; ie<RadAccessData.ne; ie++)
	for(int ie=ieStart; ie<ieBefEnd; ie++) //OC161008
	{
		if(RadAccessData.PresT == 1)
		{
			ePh = RadAccessData.avgPhotEn; //?? OC041108
		}

		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;

		ConstRxE = ConstRx*ePh;
		ConstRzE = ConstRz*ePh;

		if(RadAccessData.Pres == 1)
		{
			//double Lambda_m = 1.239854e-06/ePh;
			double Lambda_m = 1.239842e-06/ePh;
			if(RadAccessData.PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

			double Lambda_me2 = Lambda_m*Lambda_m;
			ConstRxE *= Lambda_me2;
			ConstRzE *= Lambda_me2;
		}

		z = RadAccessData.zStart - RadAccessData.zc; //To check: this is probably not correct in Angular representation?

		zE2 = z*z;
		double PhaseAddZ = 0.;
		if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE*zE2;

		for(int iz=0; iz<RadAccessData.nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			float *pEX_StartForX = pEX0 + izPerZ;
			float *pEZ_StartForX = pEZ0 + izPerZ;

			x = RadAccessData.xStart - RadAccessData.xc; //To check: this is probably not correct in Angular representation?

			for(int ix=0; ix<RadAccessData.nx; ix++)
			{
				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie;

				//Phase = ConstRxE*x*x + ConstRzE*zE2;
				Phase = PhaseAddZ;
				if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeX) Phase += ConstRxE*x*x;

				//AuxFFT2D.CosAndSin(Phase, CosPh, SinPh);
				CosAndSin(Phase, CosPh, SinPh);

				if(TreatPolCompX)
				{
					float *pExRe = pEX_StartForX + ixPerX_p_Two_ie;
					float *pExIm = pExRe + 1;
					double ExReNew = (*pExRe)*CosPh - (*pExIm)*SinPh;
					double ExImNew = (*pExRe)*SinPh + (*pExIm)*CosPh;
					*pExRe = (float)ExReNew; *pExIm = (float)ExImNew;
				}
				if(TreatPolCompZ)
				{
					float *pEzRe = pEZ_StartForX + ixPerX_p_Two_ie;
					float *pEzIm = pEzRe + 1;
					double EzReNew = (*pEzRe)*CosPh - (*pEzIm)*SinPh;
					double EzImNew = (*pEzRe)*SinPh + (*pEzIm)*CosPh;
					*pEzRe = (float)EzReNew; *pEzIm = (float)EzImNew;
				}

				x += RadAccessData.xStep;
			}
			z += RadAccessData.zStep;
			zE2 = z*z;
			PhaseAddZ = 0.;
			if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE*zE2;
		}
		ePh += RadAccessData.eStep;
	}
}

//*************************************************************************

//void srTGenOptElem::TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData& RadAccessData, float* arRayTrCoord, float xMin, float xMax, float zMin, float zMax, char AddOrRem, char PolComp, int ieOnly)
//void srTGenOptElem::TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData& RadAccessData, double* arRayTrCoord, double xMin, double xMax, double zMin, double zMax, char AddOrRem, char PolComp, int ieOnly, double anamorphMagnX, double anamorphMagnZ)
void srTGenOptElem::TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData& RadAccessData, double* arRayTrCoord, double xMin, double xMax, double zMin, double zMax, char AddOrRem, char PolComp, int ieOnly)
{
	//Later treat X and Z coordinates separately here!!!

	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x')) && (RadAccessData.pBaseRadX != 0); //OC13112011
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z')) && (RadAccessData.pBaseRadZ != 0);

	double Rx = RadAccessData.RobsX;
	double Rz = RadAccessData.RobsZ;

	const double Pi = 3.14159265358979;
	double Const = Pi*1.E+06/1.239854; // Assumes m and eV

	double ConstRx = (RadAccessData.Pres == 0)? Const/Rx : -Const*Rx;
	double ConstRz = (RadAccessData.Pres == 0)? Const/Rz : -Const*Rz;

	if(AddOrRem == 'r') { ConstRx = -ConstRx; ConstRz = -ConstRz;}

	double ConstRxE, ConstRzE;
	double ePh = RadAccessData.eStart, x, z; //, zE2;
	double Phase;
	//float CosPh, SinPh;
	double CosPh, SinPh;

	float *pEX0 = 0, *pEZ0 = 0;
	if(TreatPolCompX) pEX0 = RadAccessData.pBaseRadX;
	if(TreatPolCompZ) pEZ0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	int ieStart=0, ieBefEnd=RadAccessData.ne;
	if((ieOnly >= 0) && (ieOnly < RadAccessData.ne))
	{
		ieStart = ieOnly; ieBefEnd = ieOnly + 1;
	}

	//const float relStepTol = (float)1.E-03;
	const double relStepTol = 1.E-03;

	//float xStepTol = (float)(relStepTol*fabs(RadAccessData.xStep));
	//float zStepTol = (float)(relStepTol*fabs(RadAccessData.zStep));
	//float xMin_mi_xStepTol = xMin - xStepTol;
	//float xMax_pl_xStepTol = xMax + xStepTol;
	//float zMin_mi_zStepTol = zMin - zStepTol;
	//float zMax_pl_zStepTol = zMax + zStepTol;

	double xStepTol = relStepTol*fabs(RadAccessData.xStep);
	double zStepTol = relStepTol*fabs(RadAccessData.zStep);
	double xMin_mi_xStepTol = xMin - xStepTol;
	double xMax_pl_xStepTol = xMax + xStepTol;
	double zMin_mi_zStepTol = zMin - zStepTol;
	double zMax_pl_zStepTol = zMax + zStepTol;

	//double invAnamorphMagnXe2 = 1./(anamorphMagnX*anamorphMagnX);
	//double invAnamorphMagnZe2 = 1./(anamorphMagnZ*anamorphMagnZ);

	//for(int ie=0; ie<RadAccessData.ne; ie++)
	for(int ie=ieStart; ie<ieBefEnd; ie++) //OC161008
	{
		if(RadAccessData.PresT == 1)
		{
			ePh = RadAccessData.avgPhotEn; //?? OC041108
		}

		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;

		ConstRxE = ConstRx*ePh;
		ConstRzE = ConstRz*ePh;

		if(RadAccessData.Pres == 1)
		{
			//double Lambda_m = 1.239854e-06/ePh;
			double Lambda_m = 1.239842e-06/ePh;
			if(RadAccessData.PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

			double Lambda_me2 = Lambda_m*Lambda_m;
			ConstRxE *= Lambda_me2;
			ConstRzE *= Lambda_me2;
		}

		//z = RadAccessData.zStart - RadAccessData.zc;
		//zE2 = z*z;
		//double PhaseAddZ = 0.;
		//if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE*zE2;

		for(int iz=0; iz<RadAccessData.nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			float *pEX_StartForX = pEX0 + izPerZ;
			float *pEZ_StartForX = pEZ0 + izPerZ;

			//float *pRayTrCrd_StartForX = arRayTrCoord + izPerZ;
			double *pRayTrCrd_StartForX = arRayTrCoord + izPerZ;

			//x = RadAccessData.xStart - RadAccessData.xc;

			for(int ix=0; ix<RadAccessData.nx; ix++)
			{
				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				//long ixPerX_p_Two_ie_p_1 = ixPerX_p_Two_ie + 1;
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie_p_1 = ixPerX_p_Two_ie + 1;

				x = *(pRayTrCrd_StartForX + ixPerX_p_Two_ie);
				if((xMin_mi_xStepTol <= x) && (x <= xMax_pl_xStepTol))
				{
					x -= RadAccessData.xc;
					z = *(pRayTrCrd_StartForX + ixPerX_p_Two_ie_p_1);
					if((zMin_mi_zStepTol <= z) && (z <= zMax_pl_zStepTol))
					{
						z -= RadAccessData.zc;
						//Phase = ConstRxE*x*x + ConstRzE*zE2;
						//Phase = PhaseAddZ;
						Phase = 0;
						if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeX) Phase += ConstRxE*x*x;
						if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) Phase += ConstRzE*z*z;
						//if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeX) Phase += ConstRxE*x*x*invAnamorphMagnXe2; //OC220214
						//if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) Phase += ConstRzE*z*z*invAnamorphMagnZe2;

						//CosAndSin(Phase, CosPh, SinPh);
						CosPh = cos(Phase); SinPh = sin(Phase); //OC260114

						if(TreatPolCompX)
						{
							float *pExRe = pEX_StartForX + ixPerX_p_Two_ie;
							float *pExIm = pExRe + 1;
							double ExReNew = (*pExRe)*CosPh - (*pExIm)*SinPh;
							double ExImNew = (*pExRe)*SinPh + (*pExIm)*CosPh;
							*pExRe = (float)ExReNew; *pExIm = (float)ExImNew;
							//test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
							//*pExRe = (float)z; *pExIm = (float)z;
						}
						if(TreatPolCompZ)
						{
							float *pEzRe = pEZ_StartForX + ixPerX_p_Two_ie;
							float *pEzIm = pEzRe + 1;
							double EzReNew = (*pEzRe)*CosPh - (*pEzIm)*SinPh;
							double EzImNew = (*pEzRe)*SinPh + (*pEzIm)*CosPh;
							*pEzRe = (float)EzReNew; *pEzIm = (float)EzImNew;
						}
					}
				}
			}
		}
		ePh += RadAccessData.eStep;
	}
}

//*************************************************************************

char srTGenOptElem::WaveFrontTermCanBeTreated1D(srTRadSect1D& RadSect1D)
{
	const double CritRatTransvLong = 0.1;
	const double CritRelRobsErr = 0.1; //0.2;

	const double Pi = 3.14159265358979;

	char RobsErrIsSmall = ::fabs(RadSect1D.RobsAbsErr) < CritRelRobsErr*(::fabs(RadSect1D.Robs));
	if(RadSect1D.Pres == 0) // Coord
	{
		double Magn = ::fabs((RadSect1D.np - 1)*RadSect1D.ArgStep);
		char AnglesAreSmall = (Magn < CritRatTransvLong*(::fabs(RadSect1D.Robs)));
		return (AnglesAreSmall && RobsErrIsSmall);
	}
	else // Ang
	{// Not sure for this...
		double ePh = RadSect1D.eVal;

		//double Lambda_m = 1.239854e-06/ePh;
		double Lambda_m = 1.239842e-06/ePh;
		double RatMax = ::fabs((RadSect1D.np - 1)*RadSect1D.ArgStep*Lambda_m);
		double MaxPhaseChange = ::fabs((Pi/Lambda_m)*(RadSect1D.Robs*RatMax*RatMax));

		char AnglesAreSmall = (RatMax < CritRatTransvLong);
		char PhaseChangeIsLarge = (MaxPhaseChange > 2.*Pi);

		return (AnglesAreSmall && PhaseChangeIsLarge && RobsErrIsSmall);
	}
}

//*************************************************************************

void srTGenOptElem::TreatStronglyOscillatingTerm1D(srTRadSect1D& RadSect1D, char AddOrRem)
{
	const double Pi = 3.14159265358979;
	double Const = Pi*1.E+06/1.239854; // Assumes m and eV

	double ConstR = (RadSect1D.Pres == 0)? Const/RadSect1D.Robs : -Const*RadSect1D.Robs;
	if(AddOrRem == 'r') ConstR = -ConstR;

	double ePh = RadSect1D.eVal, Arg, Phase;
	double ConstRE = ConstR*ePh;
	float CosPh, SinPh;

	float *pEX0 = RadSect1D.pEx, *pEZ0 = RadSect1D.pEz;
	//long Per = 2;
	long long Per = 2;

	if(RadSect1D.Pres == 1)
	{
		//double Lambda_m = 1.239854e-06/ePh;
		double Lambda_m = 1.239842e-06/ePh;
		if(RadSect1D.PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

		double Lambda_me2 = Lambda_m*Lambda_m;
		ConstRE *= Lambda_me2;
	}

	Arg = RadSect1D.ArgStart - RadSect1D.cArg;

	for(int i=0; i<RadSect1D.np; i++)
	{
		//long iPer = i*Per;
		long long iPer = i*Per;
		Phase = ConstRE*Arg*Arg;
		CosAndSin(Phase, CosPh, SinPh);

		float *pExRe = pEX0 + iPer;
		float *pExIm = pExRe + 1;
		double ExReNew = (*pExRe)*CosPh - (*pExIm)*SinPh;
		double ExImNew = (*pExRe)*SinPh + (*pExIm)*CosPh;
		*pExRe = (float)ExReNew; *pExIm = (float)ExImNew;

		float *pEzRe = pEZ0 + iPer;
		float *pEzIm = pEzRe + 1;
		double EzReNew = (*pEzRe)*CosPh - (*pEzIm)*SinPh;
		double EzImNew = (*pEzRe)*SinPh + (*pEzIm)*CosPh;
		*pEzRe = (float)EzReNew; *pEzIm = (float)EzImNew;

		Arg += RadSect1D.ArgStep;
	}
}

//*************************************************************************

int srTGenOptElem::RadRearrangeOnRegularMesh(srTSRWRadStructAccessData* pRadAccessData, float* CoordX, float* CoordZ)
{
	// Continue here

	return 0;
}

//*************************************************************************

int srTGenOptElem::GenExtractPhase(srTWaveAccessData& InWaveAccessData, double* pOutMag, double* pOutPhase, int xToRemove, int zToRemove)
{
	if(InWaveAccessData.AmOfDims == 1) 
	{
		double* AuxMagArray = new double[InWaveAccessData.DimSizes[0]];
		if(AuxMagArray == 0) return MEMORY_ALLOCATION_FAILURE;
		double* AuxPhaseArray = new double[InWaveAccessData.DimSizes[0]];
		if(AuxPhaseArray == 0) return MEMORY_ALLOCATION_FAILURE;
		
		int result;
		if(result = ExtractPhase1D(InWaveAccessData, AuxMagArray, AuxPhaseArray, -1, 0.)) return result;

		int RemoveCount = 0;
		double *tOutMag = pOutMag, *tOutPhase = pOutPhase;
		double *tAuxMagArray = AuxMagArray, *tAuxPhaseArray = AuxPhaseArray;
		for(long i=0; i<InWaveAccessData.DimSizes[0]; i++)
		{
			if(++RemoveCount > xToRemove)
			{
				*(tOutMag++) = *tAuxMagArray; 
				*(tOutPhase++) = *tAuxPhaseArray; 
				RemoveCount = 0;
			}
			tAuxMagArray++; tAuxPhaseArray++;
		}
		if(AuxMagArray != 0) delete[] AuxMagArray;
		if(AuxPhaseArray != 0) delete[] AuxPhaseArray;
	}
	else
	{
		long nx = InWaveAccessData.DimSizes[0], nz = InWaveAccessData.DimSizes[1];

		double* AuxMagArrayX = new double[nx];
		if(AuxMagArrayX == 0) return MEMORY_ALLOCATION_FAILURE;
		double* AuxPhaseArrayX = new double[nx];
		if(AuxPhaseArrayX == 0) return MEMORY_ALLOCATION_FAILURE;
		double* AuxMagArrayZ = new double[nz];
		if(AuxMagArrayZ == 0) return MEMORY_ALLOCATION_FAILURE;
		double* AuxPhaseArrayZ = new double[nz];
		if(AuxPhaseArrayZ == 0) return MEMORY_ALLOCATION_FAILURE;
		float* AuxContainer = new float[nz << 1];
		if(AuxContainer == 0) return MEMORY_ALLOCATION_FAILURE;

		//long iz, ixContrPass = 0;
		//long PerZ = nx << 1;
		long long iz, ixContrPass = 0;
		long long PerZ = nx << 1;

		float *tAuxContainer = AuxContainer, *tInWaveData = (float*)InWaveAccessData.pWaveData + (ixContrPass << 1);
		for(iz=0; iz<nz; iz++)
		{
			*(tAuxContainer++) = *tInWaveData; *(tAuxContainer++) = *(tInWaveData+1);
			tInWaveData += PerZ;
		}

		srTWaveAccessData AuxWaveAccessData;
		AuxWaveAccessData.pWaveData = (char*)AuxContainer;
		AuxWaveAccessData.WaveType[0] = 'c';
		AuxWaveAccessData.WaveType[1] = 'f';
		AuxWaveAccessData.AmOfDims = 1;
		AuxWaveAccessData.DimSizes[0] = nz;
		AuxWaveAccessData.DimStartValues[0] = InWaveAccessData.DimStartValues[1];
		AuxWaveAccessData.DimSteps[0] = InWaveAccessData.DimSteps[1];
		int result;
		if(result = ExtractPhase1D(AuxWaveAccessData, AuxMagArrayZ, AuxPhaseArrayZ, -1, 0.)) return result;

		AuxWaveAccessData.DimSizes[0] = nx;
		AuxWaveAccessData.DimStartValues[0] = InWaveAccessData.DimStartValues[0];
		AuxWaveAccessData.DimSteps[0] = InWaveAccessData.DimSteps[0];

		//long OffsetForX = 0;
		long long OffsetForX = 0;
		double *tOutMag = pOutMag, *tOutPhase = pOutPhase;
		int zRemoveCount = 0;
		for(iz=0; iz<nz; iz++)
		{
			if(++zRemoveCount > zToRemove)
			{
				AuxWaveAccessData.pWaveData = (char*)((float*)InWaveAccessData.pWaveData + OffsetForX);
				if(result = ExtractPhase1D(AuxWaveAccessData, AuxMagArrayX, AuxPhaseArrayX, ixContrPass, AuxPhaseArrayZ[iz])) return result;

				double *tAuxMagArrayX = AuxMagArrayX, *tAuxPhaseArrayX = AuxPhaseArrayX;
				int xRemoveCount = 0;
				//for(long ix=0; ix<nx; ix++)
				for(long long ix=0; ix<nx; ix++)
				{
					if(++xRemoveCount > xToRemove)
					{
						*(tOutMag++) = *tAuxMagArrayX; 
						*(tOutPhase++) = *tAuxPhaseArrayX;
						xRemoveCount = 0;
					}
					tAuxMagArrayX++; tAuxPhaseArrayX++;
				}
				zRemoveCount = 0;
			}
			OffsetForX += PerZ;
		}

		if(AuxMagArrayX != 0) delete[] AuxMagArrayX;
		if(AuxPhaseArrayX != 0) delete[] AuxPhaseArrayX;
		if(AuxMagArrayZ != 0) delete[] AuxMagArrayZ;
		if(AuxPhaseArrayZ != 0) delete[] AuxPhaseArrayZ;
		if(AuxContainer != 0) delete[] AuxContainer;
	}

	return 0;
}

//*************************************************************************

//int srTGenOptElem::ExtractPhase1D(srTWaveAccessData& InWaveAccessData, double* pOutMag, double* pOutPhase, long i0, double Phi0)
int srTGenOptElem::ExtractPhase1D(srTWaveAccessData& InWaveAccessData, double* pOutMag, double* pOutPhase, long long i0, double Phi0)
{
	const double Pi = 3.14159265358979;
	//const double TwoPi = 6.2831853071796;
	const int MaxAmOfPer = 16;

	float* pStartData = (float*)(InWaveAccessData.pWaveData);
	//long iTot = InWaveAccessData.DimSizes[0];
	long long iTot = InWaveAccessData.DimSizes[0];

	double PhPrev = FormalPhase(*pStartData, *(pStartData+1));
	*pOutPhase = PhPrev;
	*pOutMag = FormalMag(*pStartData, *(pStartData+1), PhPrev);
	//long kPrev = 0;
	long long kPrev = 0;

	float* tData = pStartData + 2;
	double* tMag = pOutMag + 1;
	double* tPhase = pOutPhase + 1;
	//long i;
	long long i;
	for(i=1; i<iTot; i++)
	{
		double PredictedPhase = PredictPhase(*(tPhase-1), *(tData-2), *(tData-1), *tData, *(tData+1));
		//long kStart = kPrev - (MaxAmOfPer >> 1);
		long long kStart = kPrev - (MaxAmOfPer >> 1);
		double AuxPhase = FormalPhase(*tData, *(tData+1)) + kStart*Pi;

		double MinPhaseDiff = 1.E+23;
		double GoodPhase = AuxPhase;

		//long kNewPrev;
		long long kNewPrev;
		//for(long k=0; k<MaxAmOfPer; k++)
		for(long long k=0; k<MaxAmOfPer; k++)
		{
			double PhaseDiff = ::fabs(AuxPhase - PredictedPhase);
			if(PhaseDiff < MinPhaseDiff) { MinPhaseDiff = PhaseDiff; kNewPrev = k; GoodPhase = AuxPhase;}
			AuxPhase += Pi;
		}
		kPrev = kStart + kNewPrev;
		*(tPhase++) = GoodPhase;
		*(tMag++) = FormalMag(*tData, *(tData+1), GoodPhase);
		tData += 2;
	}

	if(i0 >= 0)
	{
		double PhaseCorr = Phi0 - pOutPhase[i0];
		tPhase = pOutPhase;
		for(i=0; i<iTot; i++) *(tPhase++) += PhaseCorr;
	}

	return 0;
}

//*************************************************************************

int srTGenOptElem::GenAuxPropagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData, double* OptElem4x4Matr, double* OptElem4Vect)
{
	int i, j, k;
	DOUBLE *OldMatrStrPtrs[4];
	DOUBLE *OldVect = pRadAccessData->p4x4PropMatr + 16;
	double ResMatr[16], ResVect[4];
	double *ResMatrStrPtrs[4], *OptElemMatrStrPtrs[4];
	for(i=0; i<4; i++)
	{
		int i4 = i << 2;
		ResMatrStrPtrs[i] = ResMatr + i4;
		OldMatrStrPtrs[i] = pRadAccessData->p4x4PropMatr + i4;
		OptElemMatrStrPtrs[i] = OptElem4x4Matr + i4;
	}

	for(i=0; i<4; i++) // String No
	{
		double& ResVect_i = ResVect[i];
		ResVect_i = 0.;

		double* OptElemMatrStrPtrs_i = OptElemMatrStrPtrs[i];
		double* ResMatrStrPtrs_i = ResMatrStrPtrs[i];

		for(j=0; j<4; j++) // Column No
		{
			double Res_ij = 0.;
			for(k=0; k<4; k++) 
			{
				Res_ij += (*(OptElemMatrStrPtrs_i + k))*(*(OldMatrStrPtrs[k] + j));
			}
			*(ResMatrStrPtrs_i + j) = Res_ij;
			ResVect_i += (*(OptElemMatrStrPtrs_i + j))*OldVect[j];
		}
	}

	DOUBLE* t4x4PropMatr = pRadAccessData->p4x4PropMatr;
	double* tResMatr = ResMatr;
	for(i=0; i<16; i++) *(t4x4PropMatr++) = *(tResMatr++);

	for(i=0; i<4; i++) OldVect[i] = ResVect[i] + OptElem4Vect[i];

	return 0;
}

//*************************************************************************

int srTGenOptElem::FillOutRadFromInRad(srTSRWRadStructAccessData* pInRadAccessData, srTSRWRadStructAccessData* pOutRadAccessData)
{
	// Fill here

	return 0;
}

//*************************************************************************

int srTGenOptElem::MakeSimpleOversamplingTestAndCorrection(srTSRWRadStructAccessData* pRadAccessData)
{
	if(pRadAccessData->Pres == 1) return 0; // Do nothing if Angular repr.

	const double AcceptedMisfitRatio = 1.6;
	const double pmTolRatio = 1.1;

	double pxm = 1., pzm = 1.;

	double LambdaStart_m;
	if(pRadAccessData->PhotEnergyUnit == 0) LambdaStart_m = 1.239854*1.E-06/pRadAccessData->eStart; // eV
	else if(pRadAccessData->PhotEnergyUnit == 1) LambdaStart_m = 1.239854*1.E-09/pRadAccessData->eStart; // keV
	double HalfLambdaStart_m = 0.5*LambdaStart_m;

	double Abs_xWfrMin = ::fabs(pRadAccessData->xWfrMin), Abs_xWfrMax = ::fabs(pRadAccessData->xWfrMax);
	double xWfrM = (Abs_xWfrMin > Abs_xWfrMax)? Abs_xWfrMin : Abs_xWfrMax;

	char BetterDoNotTouchX = ((Abs_xWfrMin < pRadAccessData->xStep) || (Abs_xWfrMax < pRadAccessData->xStep));
	if(!BetterDoNotTouchX)
	{
		double RobsX_m = pRadAccessData->RobsX;
		double dxRequired = HalfLambdaStart_m*RobsX_m/xWfrM;
		double StepX_m = pRadAccessData->xStep;
		double xRatioRequired = dxRequired/StepX_m;

		double xMinAbs = ::fabs(pRadAccessData->xStart), xMaxAbs = ::fabs(pRadAccessData->xStart + pRadAccessData->xStep*(pRadAccessData->nx - 1));
		double CurRatioX = (Abs_xWfrMin > Abs_xWfrMax)? ::fabs(xMinAbs/Abs_xWfrMin) : ::fabs(xMaxAbs/Abs_xWfrMax);

		double xRatioStillNeeded = xRatioRequired/CurRatioX;
		if(xRatioStillNeeded > AcceptedMisfitRatio*pmTolRatio) pxm = xRatioStillNeeded/AcceptedMisfitRatio;
	}

	double Abs_zWfrMin = ::fabs(pRadAccessData->zWfrMin), Abs_zWfrMax = ::fabs(pRadAccessData->zWfrMax);
	double zWfrM = (Abs_zWfrMin > Abs_zWfrMax)? Abs_zWfrMin : Abs_zWfrMax;

	char BetterDoNotTouchZ = ((Abs_zWfrMin < pRadAccessData->zStep) || (Abs_zWfrMax < pRadAccessData->zStep));
	if(!BetterDoNotTouchZ)
	{
		double RobsZ_m = pRadAccessData->RobsZ;
		double dzRequired = HalfLambdaStart_m*RobsZ_m/zWfrM;
		double StepZ_m = pRadAccessData->zStep;
		double zRatioRequired = dzRequired/StepZ_m;

		double zMinAbs = ::fabs(pRadAccessData->zStart), zMaxAbs = ::fabs(pRadAccessData->zStart + pRadAccessData->zStep*(pRadAccessData->nz - 1));
		double CurRatioZ = (Abs_zWfrMin > Abs_zWfrMax)? ::fabs(zMinAbs/Abs_zWfrMin) : ::fabs(zMaxAbs/Abs_zWfrMax);

		double zRatioStillNeeded = zRatioRequired/CurRatioZ;
		if(zRatioStillNeeded > AcceptedMisfitRatio*pmTolRatio) pzm = zRatioStillNeeded/AcceptedMisfitRatio;
	}

	if((pxm != 1.) || (pzm != 1.))
	{
		srTRadResize RadResizeStruct;
		RadResizeStruct.pxm = pxm; RadResizeStruct.pxd = 1.;
		RadResizeStruct.pzm = pzm; RadResizeStruct.pzd = 1.;

		int result;
		if(result = RadResizeGen(*pRadAccessData, RadResizeStruct)) return result;
	}

	return 0;
}

//*************************************************************************

int srTGenOptElem::PropagateRadiationTest(srTSRWRadStructAccessData* pInRadAccessData, srTSRWRadStructAccessData* pOutRadAccessData) 
{
	srTRadResizeVect AuxResizeVect;
	srTParPrecWfrPropag ParPrecWfrPropag;
	ParPrecWfrPropag.MethNo = 1;
	return PropagateRadiation(pInRadAccessData, ParPrecWfrPropag, AuxResizeVect);
}

//*************************************************************************

void srTGenOptElem::ExtractTransmCharact(int CharType, double xc, double xr, int nx, double zc, double zr, int nz, float* pData)
{//virtual; those optical elements which don't conform with this treatment should implement another version of this function
	if((CharType < 1) || (CharType > 3)) throw TRANSM_CHARACT_NOT_DEFINED;
	if((nx <= 0) || (nz <= 0)) throw NUM_POINTS_SHOULD_BE_POSITIVE;
	if(pData == 0) throw TRANSM_CHARACT_NOT_DEFINED;

	if((CharType == 1) || (CharType == 2)) //amplitude or intensity transmission
	{
		double sAux = 10., eSt = 1000., eFi = 1000.;
		int eN = 1;
		//srTWfrSmp AuxWfrSmp(sAux, xc - 0.5*xr, xc + 0.5*xr, nx, zc - 0.5*zr, zc + 0.5*zr, nz, eSt, eFi, eN, "EV");
		srTWfrSmp AuxWfrSmp(sAux, xc - 0.5*xr, xc + 0.5*xr, nx, zc - 0.5*zr, zc + 0.5*zr, nz, 0, eSt, eFi, eN, "EV");
		srTSRWRadStructAccessData AuxRadStruct(&AuxWfrSmp, true);
        AuxRadStruct.AssignElFieldToConst((float)1., (float)0.);

		PropagateRadiationSimple(&AuxRadStruct);
		AuxRadStruct.ExtractElFieldAmplitude(CharType, pData);
	}
	else if(CharType == 3) //optical path difference
	{
		float *tData = pData;
		srTEXZ LocEXZ;
		LocEXZ.VsXorZ = 'x';
		LocEXZ.e = 1000.; //dummy

		double zStep = 0, xStep = 0;
		if(nz > 1) zStep = zr/(nz - 1);
		if(nx > 1) xStep = xr/(nx - 1);

		LocEXZ.z = zc - 0.5*zr;
		for(int iz=0; iz<nz; iz++)
		{
            LocEXZ.x = xc - 0.5*xr;
			for(int ix=0; ix<nx; ix++)
			{
				*(tData++) = (float)RadOptPathDiff(LocEXZ);
				LocEXZ.x += xStep;
			}
			LocEXZ.z += zStep;
		}
	}
}

//*************************************************************************
