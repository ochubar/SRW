/************************************************************************//**
 * File: sroptdrf.cpp
 * Description: Optical element: Drift space
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptdrf.h"
#include "gmfft.h"
#include "srradmnp.h"
#include "gmmeth.h"

//*************************************************************************

//int srTDriftSpace::AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, float** ax, float** az, srTMomentsRatios* MomRatArray)
int srTDriftSpace::AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, double** ax, double** az, srTMomentsRatios* MomRatArray) //OC130311
{// There is a general function like this
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

	double SigXe2, MinSigX, MinSigXe2, SigZe2, MinSigZ, MinSigZe2;

	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	int OffsetE = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + OffsetE);
		if(*(MomX.pTotPhot) != 0.)
		{
			//float OldMomX_xx = *(MomX.pXX), OldMomX_xpxp = *(MomX.pXPXP);
			//float OldMomX_zz = *(MomX.pZZ), OldMomX_zpzp = *(MomX.pZPZP);
			double OldMomX_xx = *(MomX.pXX), OldMomX_xpxp = *(MomX.pXPXP);
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
			
			//Protection against zero or negative spot sizes 
			SigXe2 = (*(MomX.pXX)) - (*(MomX.pX))*(*(MomX.pX));
			MinSigX = DiffractionLimitedPropagatedSpotSize('x', pRadAccessData, ie);
			MinSigXe2 = MinSigX*MinSigX;
			if(SigXe2 < MinSigXe2) *(MomX.pXX) = (MinSigXe2 + (*(MomX.pX))*(*(MomX.pX)));
			
			SigZe2 = (*(MomX.pZZ)) - (*(MomX.pZ))*(*(MomX.pZ));
			MinSigZ = DiffractionLimitedPropagatedSpotSize('z', pRadAccessData, ie);
			MinSigZe2 = MinSigZ*MinSigZ;
			if(SigZe2 < MinSigZe2) *(MomX.pZZ) = (MinSigZe2 + (*(MomX.pZ))*(*(MomX.pZ)));
			
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomX = ((*(MomX.pXX) > 0.)? sqrt(*(MomX.pXX)/OldMomX_xx) : -1.);
				tMomRatArray->RxpxpMomX = ((*(MomX.pXPXP) > 0.)? sqrt(*(MomX.pXPXP)/OldMomX_xpxp) : -1.);
				tMomRatArray->RzzMomX = ((*(MomX.pZZ) > 0.)? sqrt(*(MomX.pZZ)/OldMomX_zz) : -1.);
				tMomRatArray->RzpzpMomX = ((*(MomX.pZPZP) > 0.)? sqrt(*(MomX.pZPZP)/OldMomX_zpzp) : -1.);
			}
		}
		else
		{
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomX = 1;
				tMomRatArray->RxpxpMomX = 1;
				tMomRatArray->RzzMomX = 1;
				tMomRatArray->RzpzpMomX = 1;
			}
		}

		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + OffsetE);
		if(*(MomZ.pTotPhot) != 0.)
		{
			//float OldMomZ_xx = *(MomZ.pXX), OldMomZ_xpxp = *(MomZ.pXPXP);
			//float OldMomZ_zz = *(MomZ.pZZ), OldMomZ_zpzp = *(MomZ.pZPZP);
			double OldMomZ_xx = *(MomZ.pXX), OldMomZ_xpxp = *(MomZ.pXPXP);
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
			
			//Protection against zero or negative spot sizes 
			SigXe2 = (*(MomZ.pXX)) - (*(MomZ.pX))*(*(MomZ.pX));
			if(SigXe2 < MinSigXe2) *(MomZ.pXX) = (MinSigXe2 + (*(MomZ.pX))*(*(MomZ.pX)));
			
			SigZe2 = (*(MomZ.pZZ)) - (*(MomZ.pZ))*(*(MomZ.pZ));
			if(SigZe2 < MinSigZe2) *(MomZ.pZZ) = (MinSigZe2 + (*(MomZ.pZ))*(*(MomZ.pZ)));
			
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomZ = ((*(MomZ.pXX) > 0.)? sqrt(*(MomZ.pXX)/OldMomZ_xx) : -1.);
				tMomRatArray->RxpxpMomZ = ((*(MomZ.pXPXP) > 0.)? sqrt(*(MomZ.pXPXP)/OldMomZ_xpxp) : -1.);
				tMomRatArray->RzzMomZ = ((*(MomZ.pZZ) > 0.)? sqrt(*(MomZ.pZZ)/OldMomZ_zz) : -1.);
				tMomRatArray->RzpzpMomZ = ((*(MomZ.pZPZP) > 0.)? sqrt(*(MomZ.pZPZP)/OldMomZ_zpzp) : -1.);
			}
		}
		else
		{
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomZ = 1;
				tMomRatArray->RxpxpMomZ = 1;
				tMomRatArray->RzzMomZ = 1;
				tMomRatArray->RzpzpMomZ = 1;
			}
		}
		
		if(FillMomRatios)
		{
			tMomRatArray++;
		}

		OffsetE += 11;
	}
	return 0;
}

//*************************************************************************

int srTDriftSpace::PropagateElecBeamMoments(srTEbmDat* pEbm)
{
	if(pEbm == 0) return 0;

    //float aStr0[] = { 1., (float)Length };
    //float aStr1[] = { 0., 1. };
    //float* a[] = { aStr0, aStr1 };
	double p4x4PropMatr[] = {
        1., Length, 0., 0.,
		0., 1., 0., 0.,
        0., 0., 1., Length,
		0., 0., 0., 1.,
	};
	double p4Vect[] = { 0., 0., 0., 0.};

	srTElecBeamMoments EbmMom(pEbm);
    srTRadGenManip::PropagateElecBeamMoments(EbmMom, p4x4PropMatr, p4Vect);
	EbmMom.OutData(pEbm);
    pEbm->s0 += Length;

	return 0;
}

//*************************************************************************

int srTDriftSpace::TuneRadForPropMeth_1(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& PostResize)
{
	srTMomentsRatios* MomRatArray = new srTMomentsRatios[pRadAccessData->ne];
	if(MomRatArray == 0) return MEMORY_ALLOCATION_FAILURE;

	int result;
	if(pRadAccessData->Pres != 0) // Go to spatial...
		if(result = SetRadRepres(pRadAccessData, 0)) return result;

	if(result = PropagateRadMoments(pRadAccessData, MomRatArray)) return result;
	
	srTMomentsRatios* tMomRatArray = MomRatArray;

	//float pxMaxMomX = (float)(1.e-23), pxMinMomX = (float)(1.e+23), pzMaxMomX = (float)(1.e-23), pzMinMomX = (float)(1.e+23);
	//float pxMaxMomZ = (float)(1.e-23), pxMinMomZ = (float)(1.e+23), pzMaxMomZ = (float)(1.e-23), pzMinMomZ = (float)(1.e+23);
	double pxMaxMomX = (1.e-23), pxMinMomX = (1.e+23), pzMaxMomX = (1.e-23), pzMinMomX = (1.e+23); //OC130311
	double pxMaxMomZ = (1.e-23), pxMinMomZ = (1.e+23), pzMaxMomZ = (1.e-23), pzMinMomZ = (1.e+23);

	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		if(tMomRatArray->RxxMomX > pxMaxMomX) pxMaxMomX = tMomRatArray->RxxMomX;
		if(tMomRatArray->RxxMomZ > pxMaxMomZ) pxMaxMomZ = tMomRatArray->RxxMomZ;
		if(tMomRatArray->RxxMomX < pxMinMomX) pxMinMomX = tMomRatArray->RxxMomX;
		if(tMomRatArray->RxxMomZ < pxMinMomZ) pxMinMomZ = tMomRatArray->RxxMomZ;

		if(tMomRatArray->RzzMomX > pzMaxMomX) pzMaxMomX = tMomRatArray->RzzMomX;
		if(tMomRatArray->RzzMomZ > pzMaxMomZ) pzMaxMomZ = tMomRatArray->RzzMomZ;
		if(tMomRatArray->RzzMomX < pzMinMomX) pzMinMomX = tMomRatArray->RzzMomX;
		if(tMomRatArray->RzzMomZ < pzMinMomZ) pzMinMomZ = tMomRatArray->RzzMomZ;
	
		tMomRatArray++;
	}

	//float pxMax = (pxMaxMomX > pxMaxMomZ)? pxMaxMomX : pxMaxMomZ;
	//float pxMin = (pxMinMomX < pxMinMomZ)? pxMinMomX : pxMinMomZ;
	//float pzMax = (pzMaxMomX > pzMaxMomZ)? pzMaxMomX : pzMaxMomZ;
	//float pzMin = (pzMinMomX < pzMinMomZ)? pzMinMomX : pzMinMomZ;
	double pxMax = (pxMaxMomX > pxMaxMomZ)? pxMaxMomX : pxMaxMomZ; //OC130311
	double pxMin = (pxMinMomX < pxMinMomZ)? pxMinMomX : pxMinMomZ;
	double pzMax = (pzMaxMomX > pzMaxMomZ)? pzMaxMomX : pzMaxMomZ;
	double pzMin = (pzMinMomX < pzMinMomZ)? pzMinMomX : pzMinMomZ;

	char xPostResizeUndefined = 0, zPostResizeUndefined = 0;
	if((pxMax < 0.) || (pxMin < 0.)) xPostResizeUndefined = 1;
	if((pzMax < 0.) || (pzMin < 0.)) zPostResizeUndefined = 1;

	srTRadResize RadResize;
	RadResize.pxm = RadResize.pxd = RadResize.pzm = RadResize.pzd = 1.;

	const double ResizeTol = 0.15;
	const double DiffractionFactor = 1.1;
	char xResizeNeeded = (pxMax - 1. > ResizeTol);
	char zResizeNeeded = (pzMax - 1. > ResizeTol);
	if(xResizeNeeded) RadResize.pxm = DiffractionFactor*pxMax;
	if(zResizeNeeded) RadResize.pzm = DiffractionFactor*pzMax;
	if(xResizeNeeded || zResizeNeeded) if(result = RadResizeGen(*pRadAccessData, RadResize)) return result;

	PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;

	if(!xPostResizeUndefined)
	{
		char xPostResizeNeeded = (1.- ResizeTol > pxMax);
		if(xPostResizeNeeded) 
		{
			PostResize.pxm = pxMax;
		}
	}
	else PostResize.pxm = -1.;

	if(!zPostResizeUndefined)
	{
		char zPostResizeNeeded = (1.- ResizeTol > pzMax);
		if(zPostResizeNeeded) 
		{
			PostResize.pzm = pzMax;
		}
	}
	else PostResize.pzm = -1.;

	if(MomRatArray != 0) delete[] MomRatArray;
	return 0;
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationMeth_1(srTSRWRadStructAccessData* pRadAccessData)
{
	int result;
	srTRadResize PostResize;
	PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;

	//float* OldMxxArr = new float[pRadAccessData->ne];
	double* OldMxxArr = new double[pRadAccessData->ne]; //OC130311
	if(OldMxxArr == 0) return MEMORY_ALLOCATION_FAILURE;
	//float* OldMzzArr = new float[pRadAccessData->ne];
	double* OldMzzArr = new double[pRadAccessData->ne]; //OC130311
	if(OldMzzArr == 0) return MEMORY_ALLOCATION_FAILURE;

	SetupMxxMzzArr(pRadAccessData, OldMxxArr, OldMzzArr);

	//float *NewMxxArr = 0, *NewMzzArr = 0;
	double *NewMxxArr = 0, *NewMzzArr = 0;

	if(result = TuneRadForPropMeth_1(pRadAccessData, PostResize)) return result;
	if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;

	if(pRadAccessData->Pres != 1) if(result = SetRadRepres(pRadAccessData, 1)) return result;
	if(result = TraverseRadZXE(pRadAccessData)) return result;
	if(result = SetRadRepres(pRadAccessData, 0)) return result;

	//const double ResizeTol = 0.15;
	if((PostResize.pxm != -1) && (PostResize.pzm != -1))
	{
		char PostResizeNeeded = (::fabs(PostResize.pxm - 1.) || ::fabs(PostResize.pzm - 1.));
		if(PostResizeNeeded) if(result = RadResizeGen(*pRadAccessData, PostResize)) return result;
	}
	else
	{
		if(result = ComputeRadMoments(pRadAccessData)) return result;

		//NewMxxArr = new float[pRadAccessData->ne];
		NewMxxArr = new double[pRadAccessData->ne]; //OC130311
		if(NewMxxArr == 0) return MEMORY_ALLOCATION_FAILURE;
		//NewMzzArr = new float[pRadAccessData->ne];
		NewMzzArr = new double[pRadAccessData->ne];
		if(NewMzzArr == 0) return MEMORY_ALLOCATION_FAILURE;
		SetupMxxMzzArr(pRadAccessData, NewMxxArr, NewMzzArr);

		//float pxmMinE2, pxmMaxE2, pzmMinE2, pzmMaxE2;
		double pxmMinE2, pxmMaxE2, pzmMinE2, pzmMaxE2; //OC130311
		FindMinMaxRatio(OldMxxArr, NewMxxArr, pRadAccessData->ne, pxmMinE2, pxmMaxE2);
		FindMinMaxRatio(OldMzzArr, NewMzzArr, pRadAccessData->ne, pzmMinE2, pzmMaxE2);

		PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;
		PostResize.pxm = sqrt(pxmMaxE2);
		PostResize.pzm = sqrt(pzmMaxE2);
		char PostResizeNeeded = (::fabs(PostResize.pxm - 1.) || ::fabs(PostResize.pzm - 1.));

		if(PostResizeNeeded) if(result = RadResizeGen(*pRadAccessData, PostResize)) return result;
	}

	if(result = Propagate4x4PropMatr(pRadAccessData)) return result;

	pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();

	if(OldMxxArr != 0) delete[] OldMxxArr;
	if(OldMzzArr != 0) delete[] OldMzzArr;
	if(NewMxxArr != 0) delete[] NewMxxArr;
	if(NewMzzArr != 0) delete[] NewMzzArr;
	return 0;
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationSimple_PropToWaist(srTSRWRadStructAccessData* pRadAccessData)
{// e in eV; Length in m !!!
	int result;

	pRadAccessData->SetNonZeroWavefrontLimitsToFullRange(); //OCtest271214

	SetupPropBufVars_PropToWaist(pRadAccessData);
	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;

	PropBufVars.PassNo = 1;
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	//OC240114 (commented-out)
	//if(result = ResizeBeforePropToWaistIfNecessary(pRadAccessData)) return result;

	//DEBUG
	//PropBufVars.PassNo = 3;
	//if(result = TraverseRadZXE(pRadAccessData)) return result;
	//END DEBUG

	double InvLambda_m = pRadAccessData->eStart*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.UseGivenStartTrValues = 0;

	CGenMathFFT2D FFT2D;

	//To remove this?
	srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
	if(result = SetupWfrEdgeCorrData(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;

	FFT2DInfo.pData = pRadAccessData->pBaseRadX;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

	//To remove this?
	if(DataPtrsForWfrEdgeCorr.WasSetup)
	{
		MakeWfrEdgeCorrection(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr);
		DataPtrsForWfrEdgeCorr.DisposeData();
	}

	double qxCen = -pRadAccessData->xc*InvLambda_m/pRadAccessData->RobsX; 
	double qzCen = -pRadAccessData->zc*InvLambda_m/pRadAccessData->RobsZ;

// Puting back some parameters
	pRadAccessData->xStart = (FFT2DInfo.xStartTr + qxCen)*LambdaM_Length;
	pRadAccessData->zStart = (FFT2DInfo.yStartTr + qzCen)*LambdaM_Length;
	pRadAccessData->xStep = FFT2DInfo.xStepTr*LambdaM_Length;
	pRadAccessData->zStep = FFT2DInfo.yStepTr*LambdaM_Length;

	PropBufVars.PassNo = 2;
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	pRadAccessData->UnderSamplingX = 1; // Assuming successful propagation to waist
	pRadAccessData->UnderSamplingZ = 1;

	return 0;
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationSimple_PropFromWaist(srTSRWRadStructAccessData* pRadAccessData)
{//Should be very similar to PropagateRadiationSimple_PropToWaist, consider merging
	int result = 0;
	
	SetupPropBufVars_PropFromWaist(pRadAccessData);
	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;

	LocalPropMode = 2; // prop. from waist
	PropBufVars.PassNo = 1;
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	double InvLambda_m = pRadAccessData->eStart*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.UseGivenStartTrValues = 0;

	//OCTEST (commented-out "edge correction")
	//srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
	//if(result = SetupWfrEdgeCorrData(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;

	CGenMathFFT2D FFT2D;
	FFT2DInfo.pData = pRadAccessData->pBaseRadX;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

	//OCTEST (commented-out "edge correction")
	//if(DataPtrsForWfrEdgeCorr.WasSetup)
	//{
	//	MakeWfrEdgeCorrection(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr);
	//	DataPtrsForWfrEdgeCorr.DisposeData();
	//}

// Re-scaling
	pRadAccessData->xStart = (FFT2DInfo.xStartTr)*LambdaM_Length;
	pRadAccessData->zStart = (FFT2DInfo.yStartTr)*LambdaM_Length;
	pRadAccessData->xStep = FFT2DInfo.xStepTr*LambdaM_Length;
	pRadAccessData->zStep = FFT2DInfo.yStepTr*LambdaM_Length;

	PropBufVars.PassNo = 2;
	if(result = TraverseRadZXE(pRadAccessData)) return result;
	return result;
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationSimple_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData)
{// e in eV; Length in m !!!
	int result = 0;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime(&start);

	SetupPropBufVars_AnalytTreatQuadPhaseTerm(pRadAccessData);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:SetupPropBufVars_AnalytTreatQuadPhaseTerm",&start);

	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:SetRadRepres 1",&start);

	PropBufVars.PassNo = 1; //Remove quadratic term from the Phase in coord. repres.
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:TraverseRadZXE 1",&start);

		//testOC09302011
		//if(Length == 34.63) return 0;

	double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
	pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
	pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
	double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;

	pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
	pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;

		pRadAccessData->WfrEdgeCorrShouldBeDone = 0;

	if(result = SetRadRepres(pRadAccessData, 1)) return result; //To angular repres.

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:SetRadRepres 2",&start);

	PropBufVars.PassNo = 2; //Loop in angular repres.
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:TraverseRadZXE 2",&start);

		if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		{
			pRadAccessData->xStartTr += xShift;
			pRadAccessData->zStartTr += zShift;
		}

	if(result = SetRadRepres(pRadAccessData, 0)) return result; //Back to coord. repres.

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:SetRadRepres 3",&start);

	pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;
		if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		{
			pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
			pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
		}

	//Change scale
	//double kx = (pRadAccessData->RobsX + Length)/pRadAccessData->RobsX;
	//pRadAccessData->xStart = kx*(pRadAccessData->xStart) - (Length/pRadAccessData->RobsX)*(pRadAccessData->xc);
	//pRadAccessData->xStep *= kx;
	pRadAccessData->xStart = PropBufVars.kx_AnalytTreatQuadPhaseTerm*(pRadAccessData->xStart) - PropBufVars.kxc_AnalytTreatQuadPhaseTerm*(pRadAccessData->xc);
	pRadAccessData->xStep *= PropBufVars.kx_AnalytTreatQuadPhaseTerm;

	//double kz = (pRadAccessData->RobsZ + Length)/pRadAccessData->RobsZ;
	//pRadAccessData->zStart = kz*(pRadAccessData->zStart) - (Length/pRadAccessData->RobsZ)*(pRadAccessData->zc);
	//pRadAccessData->zStep *= kz;
	pRadAccessData->zStart = PropBufVars.kz_AnalytTreatQuadPhaseTerm*(pRadAccessData->zStart) - PropBufVars.kzc_AnalytTreatQuadPhaseTerm*(pRadAccessData->zc);
	pRadAccessData->zStep *= PropBufVars.kz_AnalytTreatQuadPhaseTerm;

	PropBufVars.PassNo = 3; //Add new quadratic term to the Phase in coord. repres.
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:TraverseRadZXE 3",&start);

	//pRadAccessData->MirrorFieldData(sign(kx), sign(kz));
	pRadAccessData->MirrorFieldData((int)sign(PropBufVars.kx_AnalytTreatQuadPhaseTerm), (int)sign(PropBufVars.kz_AnalytTreatQuadPhaseTerm));

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":PropagateRadiationSimple_AnalytTreatQuadPhaseTerm:MirrorFieldData",&start);

	//if(kx < 0)
	if(PropBufVars.kx_AnalytTreatQuadPhaseTerm < 0)
	{
		double xEnd = pRadAccessData->xStart + (pRadAccessData->nx - 1)*pRadAccessData->xStep; //or nx ???
		pRadAccessData->xStart = xEnd;
		pRadAccessData->xStep *= -1;
	}
	//if(kz < 0)
	if(PropBufVars.kz_AnalytTreatQuadPhaseTerm < 0)
	{
		double zEnd = pRadAccessData->zStart + (pRadAccessData->nz - 1)*pRadAccessData->zStep; //or nz ???
		pRadAccessData->zStart = zEnd;
		pRadAccessData->zStep *= -1;
	}

	pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
	return result;
}

//*************************************************************************

void srTDriftSpace::SetupPropBufVars_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData)
{// Compute any necessary buf. vars

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime(&start);

	PropBufVars.xc = pRadAccessData->xc;
	PropBufVars.zc = pRadAccessData->zc;
	PropBufVars.invRx = PropBufVars.invRz = 0;
	PropBufVars.Pi_d_LambdaM_d_Rx = PropBufVars.Pi_d_LambdaM_d_Rz = 0;
	PropBufVars.Lx = PropBufVars.Lz = 0;
	PropBufVars.invRxL = PropBufVars.invRzL = 1./Length;

	const double infLarge = 1E+23;
	double Pi_d_LambdaM = pRadAccessData->eStart*2.53384080189E+06;

	PropBufVars.kx_AnalytTreatQuadPhaseTerm = infLarge;
	PropBufVars.kxc_AnalytTreatQuadPhaseTerm = infLarge;
	PropBufVars.kz_AnalytTreatQuadPhaseTerm = infLarge;
	PropBufVars.kzc_AnalytTreatQuadPhaseTerm = infLarge;

	//double Lx_eff_max, Lz_eff_max, trueRx, trueRz;
	//EstimateTrueWfrRadAndMaxLeff_AnalytTreatQuadPhaseTerm(pRadAccessData, trueRx, trueRz, Lx_eff_max, Lz_eff_max);
	double trueRx = pRadAccessData->RobsX;
	double trueRz = pRadAccessData->RobsZ;

	//if(pRadAccessData->ne > 1) PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;
	//OC120412 (commented-out)
	//OC180813 (uncommented)
	//OC151014 (commented-out for complicance with steady-state simulations for IXS at EXFEL)

	//testOC30092011
	if(!PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm)
	{
		if(PropBufVars.AnalytTreatSubType == 1) 
		{
			EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm(pRadAccessData, trueRx, trueRz);

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":SetupPropBufVars_AnalytTreatQuadPhaseTerm:EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm",&start);
		}
		//OC15102011 -- under testing (disadvantage of the previous version is the dependence of "trueR" on statistical moments)
		else if(PropBufVars.AnalytTreatSubType == 2) 
		{
			EstimateWfrRadToSub2_AnalytTreatQuadPhaseTerm(pRadAccessData, trueRx, trueRz); //OC22042013 (uncommented)
		
			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":SetupPropBufVars_AnalytTreatQuadPhaseTerm:EstimateWfrRadToSub2_AnalytTreatQuadPhaseTerm",&start);
		}
	}

	//if(pRadAccessData->RobsX != 0) 
	if(trueRx != 0) 
	{
		//PropBufVars.invRx = 1./pRadAccessData->RobsX;
		PropBufVars.invRx = 1./trueRx;
		PropBufVars.Pi_d_LambdaM_d_Rx = Pi_d_LambdaM*PropBufVars.invRx;

		//PropBufVars.kx_AnalytTreatQuadPhaseTerm = (pRadAccessData->RobsX + Length)/pRadAccessData->RobsX;
		//PropBufVars.kxc_AnalytTreatQuadPhaseTerm = Length/pRadAccessData->RobsX;
		PropBufVars.kx_AnalytTreatQuadPhaseTerm = (trueRx + Length)/trueRx;
		PropBufVars.kxc_AnalytTreatQuadPhaseTerm = Length/trueRx;

		//if(-Length != pRadAccessData->RobsX)
		if(-Length != trueRx)
		{
			PropBufVars.Lx = Length/(1. + Length*PropBufVars.invRx);
			//PropBufVars.invRxL = 1./(Length + pRadAccessData->RobsX);
			PropBufVars.invRxL = 1./(Length + trueRx);
		}
		else 
		{
			PropBufVars.Lx = infLarge;
			PropBufVars.invRxL = infLarge;
		}
	}

	//if(pRadAccessData->RobsZ != 0) 
	if(trueRz != 0) 
	{
		//PropBufVars.invRz = 1./pRadAccessData->RobsZ;
		PropBufVars.invRz = 1./trueRz;
		PropBufVars.Pi_d_LambdaM_d_Rz = Pi_d_LambdaM*PropBufVars.invRz;

		//PropBufVars.kz_AnalytTreatQuadPhaseTerm = (pRadAccessData->RobsZ + Length)/pRadAccessData->RobsZ;
		//PropBufVars.kzc_AnalytTreatQuadPhaseTerm = Length/pRadAccessData->RobsZ;
		PropBufVars.kz_AnalytTreatQuadPhaseTerm = (trueRz + Length)/trueRz;
		PropBufVars.kzc_AnalytTreatQuadPhaseTerm = Length/trueRz;

		//if(-Length != pRadAccessData->RobsZ)
		if(-Length != trueRz)
		{
			PropBufVars.Lz = Length/(1. + Length*PropBufVars.invRz);
			//PropBufVars.invRzL = 1./(Length + pRadAccessData->RobsZ);
			PropBufVars.invRzL = 1./(Length + trueRz);
		}
		else 
		{
			PropBufVars.Lz = infLarge;
			PropBufVars.invRzL = infLarge;
		}
	}

	//double Lx_eff_max, Lz_eff_max;
	//EstimateMaxLeff_AnalytTreatQuadPhaseTerm(pRadAccessData, Lx_eff_max, Lz_eff_max);

/**test OC061108
	if(::fabs(PropBufVars.Lx) > Lx_eff_max)
	{
		int signLx = sign(PropBufVars.Lx);
		PropBufVars.Lx = signLx*Lx_eff_max;

		//double Rx_eff = sign(pRadAccessData->RobsX)*infLarge;
		double Rx_eff = sign(trueRx)*infLarge;
		if(PropBufVars.Lx != Length) Rx_eff = Length*PropBufVars.Lx/(Length - PropBufVars.Lx);
		PropBufVars.invRx = 1./Rx_eff;

		double Rx_eff_p_Length = Rx_eff + Length;
		PropBufVars.invRxL = (Rx_eff_p_Length == 0.)? infLarge : 1./Rx_eff_p_Length;
		PropBufVars.kx_AnalytTreatQuadPhaseTerm = Rx_eff_p_Length/Rx_eff;
		PropBufVars.kxc_AnalytTreatQuadPhaseTerm = Length/Rx_eff;
	}

	if(::fabs(PropBufVars.Lz) > Lz_eff_max)
	{
		int signLz = sign(PropBufVars.Lz);
		PropBufVars.Lz = signLz*Lz_eff_max;

		//double Rz_eff = sign(pRadAccessData->RobsZ)*infLarge;
		double Rz_eff = sign(trueRz)*infLarge;
		if(PropBufVars.Lz != Length) Rz_eff = Length*PropBufVars.Lz/(Length - PropBufVars.Lz);
		PropBufVars.invRz = 1./Rz_eff;

		double Rz_eff_p_Length = Rz_eff + Length;
		PropBufVars.invRzL = (Rz_eff_p_Length == 0.)? infLarge : 1./Rz_eff_p_Length;
		PropBufVars.kz_AnalytTreatQuadPhaseTerm = Rz_eff_p_Length/Rz_eff;
		PropBufVars.kzc_AnalytTreatQuadPhaseTerm = Length/Rz_eff;
	}
**/

	PropBufVars.sqrt_LxLz_d_L = ::sqrt(::fabs(PropBufVars.Lx*PropBufVars.Lz))/Length;
	PropBufVars.phase_term_signLxLz = 0.25*3.141592653589793*(2. - sign(PropBufVars.Lx) - sign(PropBufVars.Lz));

	// Continue for more buf vars
}

//*************************************************************************

void srTDriftSpace::EstimateWfrRadToSub2_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& effRx, double& effRz)
{//this simple version is not based on statistical moments
	if(pRadAccessData == 0) return;

	effRx = pRadAccessData->RobsX;
	effRz = pRadAccessData->RobsZ;

	const double factRadEr = 3.; //to steer
	double minAbsEffRx = factRadEr*(pRadAccessData->RobsXAbsErr);
	double minAbsEffRz = factRadEr*(pRadAccessData->RobsZAbsErr);

	double effRx_p_L = effRx + Length;
	if(::fabs(effRx) < minAbsEffRx)
	{
		double sgnRx = (effRx >= 0)? 1. : -1.;
		effRx = minAbsEffRx*sgnRx;
	}
	else if(::fabs(effRx_p_L) < minAbsEffRx)
	{
		double sgnRx_p_L = (effRx_p_L >= 0)? 1. : -1.;
		effRx = minAbsEffRx*sgnRx_p_L - Length;
	}

	double effRz_p_L = effRz + Length;
	if(::fabs(effRz) < minAbsEffRz)
	{
		double sgnRz = (effRz >= 0)? 1. : -1.;
		effRz = minAbsEffRz*sgnRz;
	}
	else if(::fabs(effRz_p_L) < minAbsEffRz)
	{
		double sgnRz_p_L = (effRz_p_L >= 0)? 1. : -1.;
		effRz = minAbsEffRz*sgnRz_p_L - Length;
	}
}

//*************************************************************************

void srTDriftSpace::EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& effRx, double& effRz)
{
	if(pRadAccessData == 0) return;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime(&start);

	const double infLarge = 1E+23;
	const double coefAngRange = 0.2; //0.5; //0.6; //0.1; //to tune
	//const double coefCoordRange = 0.1;
	const double pi = 3.14159265359;
	const double fourPi = 4.*pi;

	double photEn0 = pRadAccessData->eStart;
	//long offsetMom = 0;
	//const int numStatMom = 11;
	int ieCen = 0;
	if((pRadAccessData->ne > 1) && (pRadAccessData->eStep > 0)) 
	{
		photEn0 = pRadAccessData->avgPhotEn;

		double d_ieCen = (photEn0 - pRadAccessData->eStart)/(pRadAccessData->eStep);
		ieCen = (int)d_ieCen;
		if((d_ieCen - ieCen) > 0.5) ieCen++;
		if(ieCen < 0) ieCen = 0;
		else if(ieCen >= pRadAccessData->ne) ieCen = pRadAccessData->ne - 1;

		//offsetMom = numStatMom*ieCen;
	}

	//srTMomentsPtrs MomX(pRadAccessData->pMomX + offsetMom), MomZ(pRadAccessData->pMomZ + offsetMom);
	srTMomentsPtrs MomX(pRadAccessData->pMomX, ieCen), MomZ(pRadAccessData->pMomZ, ieCen);

	double s1X = pRadAccessData->RobsX, s1Z = pRadAccessData->RobsZ;
	double s2X = s1X + Length, s2Z = s1Z + Length;
	double abs_s1X = ::fabs(s1X), abs_s1Z = ::fabs(s1Z);
	double abs_s2X = ::fabs(s2X), abs_s2Z = ::fabs(s2Z);

	//double SigX=0, SigXp=0, SigZ=0, SigZp=0;
	double SigXp=0, SigZp=0;
	//double SigX=0, SigZ=0; //OC13112010

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:setup",&start);

	if((*(MomX.pTotPhot) == 0) && (*(MomZ.pTotPhot) == 0)) ComputeRadMoments(pRadAccessData); //OC14092011

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:ComputeRadMoments 1",&start);

	char TreatExOrEz = (*(MomX.pTotPhot) >= *(MomZ.pTotPhot))? 'x' : 'z';

	if(TreatExOrEz == 'x')
	{
		if((!MomX.precCenMomIsOK) || (MomX.SqrtMxpxp == 0) || (MomX.SqrtMzpzp == 0))
		//if((!MomX.precCenMomIsOK) || (MomX.SqrtMxpxp == 0) || (MomX.SqrtMzpzp == 0) || ((abs_s1X <= pRadAccessData->RobsXAbsErr) && (!pRadAccessData->MomWereCalcNum)))
		{//OC13112010: uncommented
			ComputeRadMoments(pRadAccessData);

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:ComputeRadMoments 2",&start);

			MomX.ComputeCentralMoments();

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:ComputeCentralMoments 1",&start);
		}
		//SigX = MomX.SqrtMxx;
		SigXp = MomX.SqrtMxpxp;
		//SigZ = MomX.SqrtMzz;
		SigZp = MomX.SqrtMzpzp;
	}
	else
	{
		if((!MomZ.precCenMomIsOK) || (MomZ.SqrtMxpxp == 0) || (MomZ.SqrtMzpzp == 0))
		//if((!MomZ.precCenMomIsOK) || (MomZ.SqrtMxpxp == 0) || (MomZ.SqrtMzpzp == 0) || ((abs_s1Z <= pRadAccessData->RobsZAbsErr) && (!pRadAccessData->MomWereCalcNum)))
		{//OC13112010: uncommented
			ComputeRadMoments(pRadAccessData);

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:ComputeRadMoments 3",&start);

			MomZ.ComputeCentralMoments();

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:ComputeCentralMoments 2",&start);
		}
		//SigX = MomZ.SqrtMxx;
		SigXp = MomZ.SqrtMxpxp;
		//SigZ = MomZ.SqrtMzz;
		SigZp = MomZ.SqrtMzpzp;
	}

	if((SigXp == 0) || (SigZp == 0))
	//if((SigX == 0) || (SigZ == 0))
	{
		effRx = pRadAccessData->RobsX;
		if(effRx == 0) effRx = infLarge;
		effRz = pRadAccessData->RobsZ;
		if(effRz == 0) effRz = infLarge;
		return;
	}

	double lambM = 1.239842E-06/photEn0;

	double apSigXp = coefAngRange*SigXp; //apparent RMS angular divergence
	double apSigZp = coefAngRange*SigZp; //apparent RMS angular divergence
	double sminX = lambM/(fourPi*apSigXp*apSigXp);
	double sminZ = lambM/(fourPi*apSigZp*apSigZp);

	//double apSigX = SigX*coefCoordRange; //apparent RMS angular divergence
	//double apSigZ = SigZ*coefCoordRange; //apparent RMS angular divergence
	////double sminX = fourPi*apSigX*apSigX/lambM;
	////double sminZ = fourPi*apSigZ*apSigZ/lambM;
	//double sminX = apSigX*apSigX/lambM;
	//double sminZ = apSigZ*apSigZ/lambM;

	//double half_abs_s1X = 0.5*abs_s1X;
	//if(abs_s1X <= pRadAccessData->RobsXAbsErr)
	//{
	//	if(sminX < half_abs_s1X) sminX = half_abs_s1X;
	//}
	//else sminX = half_abs_s1X;

	//double half_abs_s1Z = 0.5*abs_s1Z;
	//if(abs_s1Z <= pRadAccessData->RobsZAbsErr)
	//{
	//	if(sminZ < half_abs_s1Z) sminZ = half_abs_s1Z;
	//}
	//else sminZ = half_abs_s1Z;

	double RminX = 2*sminX, RminZ = 2*sminZ;

	if(abs_s1X < sminX)
	{
		if(abs_s2X < sminX)
		{
			if(abs_s1X < abs_s2X)
			{
				if(Length > 0)
				{
					effRx = RminX;
				}
				else
				{
					effRx = -RminX;
				}
			}
			else
			{
				if(Length > 0)
				{
					effRx = -RminX - Length;
				}
				else
				{
					effRx = RminX - Length;
				}
			}
		}
		else
		{
			if(Length > 0)
			{
				effRx = RminX;
			}
			else
			{
				effRx = -RminX;
			}
		}
	}
	else
	{
		effRx = s1X + sminX*sminX/s1X;
		//effRx = s1X;
		double effRx_p_Length = effRx + Length;

		//double minRminX = 0.003*(::fabs(effRx));
		double minRminX = 0.007*(::fabs(effRx)); //OC260709
		if(RminX < minRminX) RminX = minRminX;

		if(::fabs(effRx_p_Length) < RminX)
		{
			if(effRx_p_Length < 0)
			{
				effRx = -RminX - Length;
			}
			else
			{
				effRx = RminX - Length;
			}
		}
	}

	if(abs_s1Z < sminZ)
	{
		if(abs_s2Z < sminZ)
		{
			if(abs_s1Z < abs_s2Z)
			{
				if(Length > 0)
				{
					effRz = RminZ;
				}
				else
				{
					effRz = -RminZ;
				}
			}
			else
			{
				if(Length > 0)
				{
					effRz = -RminZ - Length;
				}
				else
				{
					effRz = RminZ - Length;
				}
			}
		}
		else
		{
			if(Length > 0)
			{
				effRz = RminZ;
			}
			else
			{
				effRz = -RminZ;
			}
		}
	}
	else
	{
		effRz = s1Z + sminZ*sminZ/s1Z;
		//effRz = s1Z;
		double effRz_p_Length = effRz + Length;

		//double minRminZ = 0.003*(::fabs(effRz));
		double minRminZ = 0.007*(::fabs(effRz)); //OC260709
		if(RminZ < minRminZ) RminZ = minRminZ;

		if(::fabs(effRz_p_Length) < RminZ)
		{
			if(effRz_p_Length < 0)
			{
				effRz = -RminZ - Length;
			}
			else
			{
				effRz = RminZ - Length;
			}
		}
	}
	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm:rest",&start);
}

//*************************************************************************

void srTDriftSpace::EstimateTrueWfrRadAndMaxLeff_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& trueRx, double& trueRz, double& Lx_eff_max, double& Lz_eff_max)
{//not used anymore
	if(pRadAccessData == 0) return;

	const double infLarge = 1E+23;
	const double coefAngRange = 0.4; //0.6; //0.1; //to tune

	double SigX=0, SigXp=0, SigZ=0, SigZp=0;
	srTMomentsPtrs MomX(pRadAccessData->pMomX), MomZ(pRadAccessData->pMomZ);
	char TreatExOrEz = (*(MomX.pTotPhot) > *(MomZ.pTotPhot))? 'x' : 'z';
	if(TreatExOrEz == 'x')
	{
		SigX = MomX.SqrtMxx;
		SigXp = MomX.SqrtMxpxp;
		SigZ = MomX.SqrtMzz;
		SigZp = MomX.SqrtMzpzp;
	}
	else
	{
		SigX = MomZ.SqrtMxx;
		SigXp = MomZ.SqrtMxpxp;
		SigZ = MomZ.SqrtMzz;
		SigZp = MomZ.SqrtMzpzp;
	}

	trueRx = infLarge;
	trueRz = infLarge;

	double photEn0 = pRadAccessData->eStart;
	if(pRadAccessData->ne > 1) photEn0 = pRadAccessData->avgPhotEn;
	double Pi_d_LambdaM = photEn0*2.53384080189E+06;
	
	if((pRadAccessData->RobsX != 0) && (SigXp != 0))
	{
		//double w0x = 1./(Pi_d_LambdaM*SigXp);
		//double ax = Pi_d_LambdaM*w0x*w0x;
		double ax = 1./(Pi_d_LambdaM*SigXp*SigXp*coefAngRange*coefAngRange);
		trueRx = pRadAccessData->RobsX + ax*ax/(pRadAccessData->RobsX);
	}
	if(::fabs(trueRx) > infLarge) trueRx = infLarge;

	if((pRadAccessData->RobsZ != 0) && (SigZp != 0))
	{
		//double w0z = 1./(Pi_d_LambdaM*SigZp);
		//double az = Pi_d_LambdaM*w0z*w0z;
		double az = 1./(Pi_d_LambdaM*SigZp*SigZp*coefAngRange*coefAngRange);
		trueRz = pRadAccessData->RobsZ + az*az/(pRadAccessData->RobsZ);
	}
	if(::fabs(trueRz) > infLarge) trueRz = infLarge;

	double Xm = 0.5*pRadAccessData->xStep*(pRadAccessData->nx - 1);
	double Zm = 0.5*pRadAccessData->zStep*(pRadAccessData->nz - 1);

	Lx_eff_max = Xm/(coefAngRange*SigXp);
	Lz_eff_max = Zm/(coefAngRange*SigZp);
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationSimple_NumIntFresnel(srTSRWRadStructAccessData* pRadAccessData) 
{//OC100914 Aux. method for testing / benchmarking
//This method attempts to calculate 2D Fresnel integral using the standard numerical integration
//(e.g. for testing accuracy of the FTT-based integration)

	double *arReExVsX = new double[pRadAccessData->nx];
	double *arImExVsX = new double[pRadAccessData->nx];
	double *arReEzVsX = new double[pRadAccessData->nx];
	double *arImEzVsX = new double[pRadAccessData->nx];

	double *arReFxVsZ = new double[pRadAccessData->nz];
	double *arImFxVsZ = new double[pRadAccessData->nz];
	double *arReFzVsZ = new double[pRadAccessData->nz];
	double *arImFzVsZ = new double[pRadAccessData->nz];

	//long nTot = 2*(pRadAccessData->ne)*(pRadAccessData->nx)*(pRadAccessData->nz);
	long long nTot = (((long long)(pRadAccessData->ne))*((long long)(pRadAccessData->nx))*((long long)(pRadAccessData->nz))) << 1;
	float *arAuxEx=0, *arAuxEz=0;
	if(pRadAccessData->pBaseRadX != 0)
	{
		arAuxEx = new float[nTot];
		float *t = pRadAccessData->pBaseRadX, *tAux = arAuxEx;
		//for(long i = 0; i < nTot; i++) *(tAux++) = *(t++);
		for(long long i = 0; i < nTot; i++) *(tAux++) = *(t++);
	}
	if(pRadAccessData->pBaseRadZ != 0)
	{
		arAuxEz = new float[nTot];
		float *t = pRadAccessData->pBaseRadZ, *tAux = arAuxEz;
		//for(long i = 0; i < nTot; i++) *(tAux++) = *(t++);
		for(long long i = 0; i < nTot; i++) *(tAux++) = *(t++);
	}

	//long perX = pRadAccessData->ne << 1;
	//long perZ = perX*pRadAccessData->nx;
	long long perX = pRadAccessData->ne << 1;
	long long perZ = perX*pRadAccessData->nx;

	double invL = 1./Length;

	double ePh = pRadAccessData->eStart;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		double pi_d_lambda_m = ePh*2.533865612E+06;
		double invPiLambda_m = ePh*invL*0.80655447456E+06;

		//long iePerE = ie << 1;
		long long iePerE = ie << 1;
		double z = pRadAccessData->zStart;
		for(long iz=0; iz<pRadAccessData->nz; iz++)
		{
			//long izPerZ = iz*perZ;
			long long izPerZ = iz*perZ;
			double x = pRadAccessData->xStart;
			for(long ix=0; ix<pRadAccessData->nx; ix++)
			{			
				double z1 = pRadAccessData->zStart;
				for(long iz1=0; iz1<pRadAccessData->nz; iz1++)
				{
					double dz = (z1 - z);
					double dze2 = dz*dz;

					//long iz1PerZ = iz1*perZ;
					long long iz1PerZ = iz1*perZ;
					double x1 = pRadAccessData->xStart;
					for(long ix1=0; ix1<pRadAccessData->nx; ix1++)
					{
						double dx = (x1 - x);
						double dxe2 = dx*dx;
						double quadTerm = (dxe2 + dze2)*invL;
						double a = quadTerm*invL;
						double phase = pi_d_lambda_m*quadTerm*(1. - 0.25*a + 0.125*a*a);
						double cosPh = cos(phase), sinPh = sin(phase);

						//long ix1PerX = ix1*perX;
						//long ofst1 = iz1PerZ + ix1PerX + iePerE;
						long long ix1PerX = ix1*perX;
						long long ofst1 = iz1PerZ + ix1PerX + iePerE;
						float *pEx1 = arAuxEx + ofst1;
						float *pEz1 = arAuxEz + ofst1;

						if(arAuxEx != 0)
						{
							float ReEx1 = *pEx1, ImEx1 = *(pEx1 + 1);
							arReExVsX[ix1] = invPiLambda_m*(ReEx1*sinPh + ImEx1*cosPh);
							arImExVsX[ix1] = invPiLambda_m*(ImEx1*sinPh - ReEx1*cosPh);
						}
						if(arAuxEz != 0)
						{
							float ReEz1 = *pEz1, ImEz1 = *(pEz1 + 1);
							arReEzVsX[ix1] = invPiLambda_m*(ReEz1*sinPh + ImEz1*cosPh);
							arImEzVsX[ix1] = invPiLambda_m*(ImEz1*sinPh - ReEz1*cosPh);
						}
						x1 += pRadAccessData->xStep;
					}
					if(arAuxEx != 0)
					{
						arReFxVsZ[iz1] = CGenMathMeth::Integ1D_FuncDefByArray(arReExVsX, pRadAccessData->nx, pRadAccessData->xStep);
						arImFxVsZ[iz1] = CGenMathMeth::Integ1D_FuncDefByArray(arImExVsX, pRadAccessData->nx, pRadAccessData->xStep);
					}
					if(arAuxEz != 0)
					{
						arReFzVsZ[iz1] = CGenMathMeth::Integ1D_FuncDefByArray(arReEzVsX, pRadAccessData->nx, pRadAccessData->xStep);
						arImFzVsZ[iz1] = CGenMathMeth::Integ1D_FuncDefByArray(arImEzVsX, pRadAccessData->nx, pRadAccessData->xStep);
					}
					z1 += pRadAccessData->zStep;
				}

				//long ixPerX = ix*perX;
				//long ofst = izPerZ + ixPerX + iePerE;
				long long ixPerX = ix*perX;
				long long ofst = izPerZ + ixPerX + iePerE;
				float *pEx = pRadAccessData->pBaseRadX + ofst;
				float *pEz = pRadAccessData->pBaseRadZ + ofst;

				if(arAuxEx != 0)
				{
					*pEx = (float)CGenMathMeth::Integ1D_FuncDefByArray(arReFxVsZ, pRadAccessData->nz, pRadAccessData->zStep);
					*(pEx + 1) = (float)CGenMathMeth::Integ1D_FuncDefByArray(arImFxVsZ, pRadAccessData->nz, pRadAccessData->zStep);
				}
				if(arAuxEz != 0)
				{
					*pEz = (float)CGenMathMeth::Integ1D_FuncDefByArray(arReFzVsZ, pRadAccessData->nz, pRadAccessData->zStep);
					*(pEz + 1) = (float)CGenMathMeth::Integ1D_FuncDefByArray(arImFzVsZ, pRadAccessData->nz, pRadAccessData->zStep);
				}
				x += pRadAccessData->xStep;
			}
			z += pRadAccessData->zStep;
		}
		ePh += pRadAccessData->eStep;
	}

	if(arReExVsX != 0) delete[] arReExVsX;
	if(arImExVsX != 0) delete[] arImExVsX;
	if(arReEzVsX != 0) delete[] arReEzVsX;
	if(arImEzVsX != 0) delete[] arImEzVsX;

	if(arReFxVsZ != 0) delete[] arReFxVsZ;
	if(arImFxVsZ != 0) delete[] arImFxVsZ;
	if(arReFzVsZ != 0) delete[] arReFzVsZ;
	if(arImFzVsZ != 0) delete[] arImFzVsZ;

	if(arAuxEx != 0) delete[] arAuxEx;
	if(arAuxEz != 0) delete[] arAuxEz;
	return 0;
}

//*************************************************************************

int srTDriftSpace::PropagateRadiationSimple1D_PropToWaist(srTRadSect1D* pSect1D)
{
	int result;
	SetupPropBufVars_PropToWaist(pSect1D);
	if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;

	PropBufVars.PassNo = 1;
	if(result = TraverseRad1D(pSect1D)) return result;

	if(result = ResizeBeforePropToWaistIfNecessary1D(pSect1D)) return result;

	double InvLambda_m = pSect1D->eVal*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	//long TwoNp = pSect1D->np << 1;
	long long TwoNp = pSect1D->np << 1;
	float *AuxCont = new float[TwoNp << 1];
	if(AuxCont == 0) return MEMORY_ALLOCATION_FAILURE;

	float *pAuxX = AuxCont, *pAuxZ = AuxCont + TwoNp;
	float *tAuxX = pAuxX, *tAuxZ = pAuxZ;
	float *tEx = pSect1D->pEx, *tEz = pSect1D->pEz;
	//for(long i=0; i<TwoNp; i++)
	for(long long i=0; i<TwoNp; i++)
	{
		*(tAuxX++) = *(tEx++); 
		*(tAuxZ++) = *(tEz++);
	}

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = pSect1D->ArgStep;
	FFT1DInfo.xStart = pSect1D->ArgStart;
	FFT1DInfo.Nx = pSect1D->np;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
	CGenMathFFT1D FFT1D;

	srTDataPtrsForWfrEdgeCorr1D DataPtrsForWfrEdgeCorr1D;
	if(result = SetupWfrEdgeCorrData1D(pSect1D, pAuxX, pAuxZ, DataPtrsForWfrEdgeCorr1D)) return result;

	FFT1DInfo.pInData = pAuxX;
	FFT1DInfo.pOutData = pSect1D->pEx;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	FFT1DInfo.pInData = pAuxZ;
	FFT1DInfo.pOutData = pSect1D->pEz;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	if(DataPtrsForWfrEdgeCorr1D.WasSetup)
	{
		MakeWfrEdgeCorrection1D(pSect1D, pSect1D->pEx, pSect1D->pEz, DataPtrsForWfrEdgeCorr1D);
		DataPtrsForWfrEdgeCorr1D.DisposeData();
	}

	double qCen = -pSect1D->cArg*InvLambda_m/pSect1D->Robs;

	pSect1D->ArgStart = (FFT1DInfo.xStartTr + qCen)*LambdaM_Length;
	pSect1D->ArgStep = FFT1DInfo.xStepTr*LambdaM_Length;

	PropBufVars.PassNo = 2;
	if(result = TraverseRad1D(pSect1D)) return result;

	if(AuxCont != 0) delete[] AuxCont;
	return 0;
}

//*************************************************************************

int srTDriftSpace::ResizeBeforePropToWaistIfNecessary(srTSRWRadStructAccessData* pRadAccessData)
{
	int result;
	const int MinNafter = 60; // To steer
	const double DiffAllowResize = 0.05; // To steer

	double InvLambda_m = pRadAccessData->eStart*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	double xRange = pRadAccessData->nx*pRadAccessData->xStep;
	double zRange = pRadAccessData->nz*pRadAccessData->zStep;

	double xStartAbs = ::fabs(pRadAccessData->xStart);
	double xAbsMax = ::fabs(pRadAccessData->xStart + xRange);
	if(xAbsMax < xStartAbs) xAbsMax = xStartAbs;

	double zStartAbs = ::fabs(pRadAccessData->zStart);
	double zAbsMax = ::fabs(pRadAccessData->zStart + zRange);
	if(zAbsMax < zStartAbs) zAbsMax = zStartAbs;

	double pxm = 1.4*LambdaM_Length*pRadAccessData->UnderSamplingX/(xRange*pRadAccessData->xStep); // To steer
	if(pxm < 1.) pxm = 1.;
	if(::fabs(pxm  - 1.) < DiffAllowResize) pxm = 1.;

	double pzm = 1.4*LambdaM_Length*pRadAccessData->UnderSamplingZ/(zRange*pRadAccessData->zStep); // To steer
	if(pzm < 1.) pzm = 1.;
	if(::fabs(pzm  - 1.) < DiffAllowResize) pzm = 1.;

	int NxResWell, NzResWell;
	EstimateMinNxNzBeforePropToWaist(pRadAccessData, NxResWell, NzResWell);
	double MinNxInRange = (NxResWell > MinNafter)? NxResWell : MinNafter;
	double MinNzInRange = (NzResWell > MinNafter)? NzResWell : MinNafter;

	double ResAfter_pxm = 1;
	if(pRadAccessData->pResAfter != 0) ResAfter_pxm = pRadAccessData->pResAfter->pxm;

	double pxd = ::fabs(2.3*ResAfter_pxm*xRange*pRadAccessData->xStep/(LambdaM_Length)); 
	double pxdOutOfSpot = ::fabs(2*(pRadAccessData->RobsX + Length)*xAbsMax*pRadAccessData->xStep/((pRadAccessData->RobsX)*LambdaM_Length)); 
	if(pxd < pxdOutOfSpot) pxd = pxdOutOfSpot;

	//to improve !!!

	if(::fabs(pxd - 1.) < DiffAllowResize) pxd = 1.;

	if((pRadAccessData->nx)*pxd < MinNxInRange)
	{
		pxd = MinNxInRange/double(pRadAccessData->nx);
		if(::fabs(pxd - 1.) < DiffAllowResize) pxd = 1.;
		double xRangeNew = pxd*LambdaM_Length/pRadAccessData->xStep;

		//double xRangeShouldBe = xRange*pRadAccessData->pResAfter->pxm;
		double xRangeShouldBe = xRange*ResAfter_pxm;

		if(pRadAccessData->pResAfter != 0)
		{
			pRadAccessData->pResAfter->pxm = xRangeShouldBe/xRangeNew;
			if(pRadAccessData->pResAfter->pxm > 1.) pRadAccessData->pResAfter->pxm = 1.;
		}
	}
	else
	{
		if(pRadAccessData->pResAfter != 0) pRadAccessData->pResAfter->pxm = 1.;
	}

	double ResAfter_pzm = 1;
	if(pRadAccessData->pResAfter != 0) ResAfter_pzm = pRadAccessData->pResAfter->pzm;

	double pzd = ::fabs(1.*ResAfter_pzm*zRange*pRadAccessData->zStep/(LambdaM_Length));
	double pzdOutOfSpot = ::fabs(2*(pRadAccessData->RobsZ + Length)*zAbsMax*pRadAccessData->zStep/((pRadAccessData->RobsZ)*LambdaM_Length)); 
	if(pzd < pzdOutOfSpot) pzd = pzdOutOfSpot;

	if(::fabs(pzd - 1.) < DiffAllowResize) pzd = 1.;

	if((pRadAccessData->nz)*pzd < MinNzInRange)
	{
		pzd = MinNzInRange/double(pRadAccessData->nz);
		if(::fabs(pzd - 1.) < DiffAllowResize) pzd = 1.;

		double zRangeNew = pzd*LambdaM_Length/pRadAccessData->zStep;
		//double zRangeShouldBe = zRange*pRadAccessData->pResAfter->pzm;
		double zRangeShouldBe = zRange*ResAfter_pzm;

		if(pRadAccessData->pResAfter != 0)
		{
			pRadAccessData->pResAfter->pzm = zRangeShouldBe/zRangeNew;
			if(pRadAccessData->pResAfter->pzm > 1.) pRadAccessData->pResAfter->pzm = 1.;
		}
	}
	else
	{
		if(pRadAccessData->pResAfter != 0) pRadAccessData->pResAfter->pzm = 1.;
	}

	//if((pRadAccessData->pResAfter != 0) && (::fabs(pxm - 1.) > DiffAllowResize) || (::fabs(pxd - 1.) > DiffAllowResize) || (::fabs(pzm - 1.) > DiffAllowResize) || (::fabs(pzd - 1.) > DiffAllowResize))
	if((pRadAccessData->pResAfter != 0) && ((::fabs(pxm - 1.) > DiffAllowResize) || (::fabs(pxd - 1.) > DiffAllowResize) || (::fabs(pzm - 1.) > DiffAllowResize) || (::fabs(pzd - 1.) > DiffAllowResize)))
	{
		srTRadResize Resize;
		Resize.pxd = pxd; Resize.pxm = pxm; Resize.pzd = pzd; Resize.pzm = pzm;
		//Resize.DoNotTreatSpherTerm = 1;
		Resize.doNotTreatSpherTerm(1); //OC090311

		const double PdReduceCoef = 0.95; // To steer
		long nxCurRad = pRadAccessData->nx, nzCurRad = pRadAccessData->nz;
		double MemForResize = ExtraMemSizeForResize(nxCurRad, nzCurRad, Resize.pxm, Resize.pxd, Resize.pzm, Resize.pzd, 0);
		double MemAvail = CheckMemoryAvailable(), PrevMemForResize;
		char ResolutionWasReduced = 0;
		while(MemAvail < MemForResize)
		{
			Resize.pxd *= PdReduceCoef; 
			double xRangeNew = Resize.pxd*LambdaM_Length*pRadAccessData->UnderSamplingX/pRadAccessData->xStep;

			double xRangeShouldBe = xRange*pRadAccessData->pResAfter->pxm;
			pRadAccessData->pResAfter->pxm = xRangeShouldBe/xRangeNew;
			if(pRadAccessData->pResAfter->pxm > 1.) pRadAccessData->pResAfter->pxm = 1.;
			
			Resize.pzd *= PdReduceCoef; 
			double zRangeNew = Resize.pzd*LambdaM_Length*pRadAccessData->UnderSamplingZ/pRadAccessData->zStep;

			double zRangeShouldBe = zRange*pRadAccessData->pResAfter->pzm;
			pRadAccessData->pResAfter->pzm = zRangeShouldBe/zRangeNew;
			if(pRadAccessData->pResAfter->pzm > 1.) pRadAccessData->pResAfter->pzm = 1.;

			if(Resize.pxm > 1.) Resize.pxm *= PdReduceCoef;
			if(Resize.pzm > 1.) Resize.pzm *= PdReduceCoef;

			ResolutionWasReduced = 1;
			PrevMemForResize = MemForResize;
			long nxCurRad = pRadAccessData->nx, nzCurRad = pRadAccessData->nz;
			MemForResize = ExtraMemSizeForResize(nxCurRad, nzCurRad, Resize.pxm, Resize.pxd, Resize.pzm, Resize.pzd, 0);

			if(MemForResize >= PrevMemForResize) break; // Dangerous, yet necessary to break infinite loop
		}

		double SmallLength = 0.001*Length; // To steer
		double RxOld = pRadAccessData->RobsX, RzOld = pRadAccessData->RobsZ;
		//if(::fabs(Length + RxOld) > SmallLength) Resize.DoNotTreatSpherTerm = 0;
		//if(::fabs(Length + RzOld) > SmallLength) Resize.DoNotTreatSpherTerm = 0;
		if(::fabs(Length + RxOld) > SmallLength) Resize.doNotTreatSpherTerm(0); //OC090311
		if(::fabs(Length + RzOld) > SmallLength) Resize.doNotTreatSpherTerm(0);

		//if(!(Resize.DoNotTreatSpherTerm))
		if(!(Resize.doNotTreatSpherTerm())) //OC090311
		{
			pRadAccessData->RobsX = Length*RxOld/(Length + RxOld);
			pRadAccessData->RobsZ = Length*RzOld/(Length + RzOld);
		}

				//Resize.pxm = Resize.pzm = 1;//OC

		if(result = RadResizeGen(*pRadAccessData, Resize)) return result;

		pRadAccessData->RobsX = RxOld;
		pRadAccessData->RobsZ = RzOld;

		if(ResolutionWasReduced)
		{
			CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
		}
	}

	if(pRadAccessData->pResAfter != 0)
	{
		pRadAccessData->pResAfter->RelCenPosX = 0.5; // To check
		pRadAccessData->pResAfter->RelCenPosZ = 0.5;
	}

	return 0;
}

//*************************************************************************

void srTDriftSpace::EstimateMinNxNzBeforePropToWaist(srTSRWRadStructAccessData* pRadAccessData, int& Nx, int& Nz)
{
	//double WavelengthIn_m = 1.239854E-06/pRadAccessData->eStart;
	double WavelengthIn_m = 1.239842E-06/pRadAccessData->eStart;
	double SmallLength = 0.001*Length;

	double Length_p_Rx = Length + pRadAccessData->RobsX;
	if(::fabs(Length_p_Rx) < SmallLength) Length_p_Rx = SmallLength;
	double Rx = ::fabs(Length*pRadAccessData->RobsX/Length_p_Rx);
	double HalfLambRx = 0.5*WavelengthIn_m*Rx;
	double xStartRel = pRadAccessData->xStart - pRadAccessData->xc;
	double xEndRel = pRadAccessData->xStart + pRadAccessData->xStep*pRadAccessData->nx - pRadAccessData->xc;
	double dxStart = ::fabs(HalfLambRx/xStartRel);
	double dxEnd = ::fabs(HalfLambRx/xEndRel);
	double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/1.4; // To steer
	Nx = int(::fabs(xEndRel - xStartRel)/dx) + 1;
	if(((Nx >> 1) << 1) != Nx) Nx++;

	double Length_p_Rz = Length + pRadAccessData->RobsZ;
	if(::fabs(Length_p_Rz) < SmallLength) Length_p_Rz = SmallLength;
	double Rz = ::fabs(Length*pRadAccessData->RobsZ/Length_p_Rz);
	double HalfLambRz = 0.5*WavelengthIn_m*Rz;
	double zStartRel = pRadAccessData->zStart - pRadAccessData->zc;
	double zEndRel = pRadAccessData->zStart + pRadAccessData->zStep*pRadAccessData->nz - pRadAccessData->zc;
	double dzStart = ::fabs(HalfLambRz/zStartRel);
	double dzEnd = ::fabs(HalfLambRz/zEndRel);
	double dz = ((dzStart < dzEnd)? dzStart : dzEnd)/1.4; // To steer
	Nz = int(::fabs(zEndRel - zStartRel)/dz) + 1;
	if(((Nz >> 1) << 1) != Nz) Nz++;
}

//*************************************************************************

void srTDriftSpace::EstimateMinNpBeforePropToWaist1D(srTRadSect1D* pSect1D, int& Np)
{
	//double WavelengthIn_m = 1.239854E-06/pSect1D->eVal;
	double WavelengthIn_m = 1.239842E-06/pSect1D->eVal;
	double SmallLength = 0.001*Length;

	double Length_p_R = Length + pSect1D->Robs;
	if(::fabs(Length_p_R) < SmallLength) Length_p_R = SmallLength;
	double R = ::fabs(Length*pSect1D->Robs/Length_p_R);
	double HalfLambR = 0.5*WavelengthIn_m*R;
	double xStartRel = pSect1D->ArgStart - pSect1D->cArg;
	double xEndRel = pSect1D->ArgStart + pSect1D->ArgStep*pSect1D->np - pSect1D->cArg;
	double dxStart = ::fabs(HalfLambR/xStartRel);
	double dxEnd = ::fabs(HalfLambR/xEndRel);
	double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/1.4;  //1.3; // To steer
	Np = int(::fabs(xEndRel - xStartRel)/dx) + 1;
	if(((Np >> 1) << 1) != Np) Np++;
}

//*************************************************************************

int srTDriftSpace::ResizeBeforePropToWaistIfNecessary1D(srTRadSect1D* pSect1D)
{// ATTENTION: Mesh may be shifted (not good for steering resize parameters !)
	int result;
	const int MinNafter = 60; // To steer
	const double DiffAllowResize = 0.05; // To steer

	double InvLambda_m = pSect1D->eVal*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	double Range = pSect1D->np*pSect1D->ArgStep;

	double pm = LambdaM_Length/(Range*pSect1D->ArgStep);
	if(pm < 1.) pm = 1.;
	if(::fabs(pm  - 1.) < DiffAllowResize) pm = 1.;

	int NpResWell;
	EstimateMinNpBeforePropToWaist1D(pSect1D, NpResWell);
	double MinNpInRange = (NpResWell > MinNafter)? NpResWell : MinNafter;

	double pd = Range*pSect1D->ArgStep/LambdaM_Length;
	if(::fabs(pd - 1.) < DiffAllowResize) pd = 1.;

	if(double(pSect1D->np)*pd < MinNafter)
	{
		pd = MinNpInRange/double(pSect1D->np);
		if(::fabs(pd - 1.) < DiffAllowResize) pd = 1.;
	}

	if((::fabs(pm - 1.) > DiffAllowResize) || (::fabs(pd - 1.) > DiffAllowResize))
	{
		srTRadResize1D Resize1D;
		Resize1D.pd = pd; Resize1D.pm = pm;
		Resize1D.DoNotTreatSpherTerm = 1;

		double SmallLength = 0.001*Length; // To steer
		double R_Old = pSect1D->Robs;
		if(::fabs(Length + R_Old) > SmallLength) Resize1D.DoNotTreatSpherTerm = 0;

		if(!(Resize1D.DoNotTreatSpherTerm)) pSect1D->Robs = Length*R_Old/(Length + R_Old);
		if(result = RadResizeGen1D(*pSect1D, Resize1D)) return result;
		pSect1D->Robs = R_Old;
	}
	return 0;
}

//*************************************************************************

//void srTDriftSpace::CheckAndSubtractPhaseTermsLin(srTSRWRadStructAccessData* pRadAccessData)
//{
//	if(pRadAccessData == 0) return;
//}

//*************************************************************************

//void srTDriftSpace::CheckAndResetPhaseTermsLin(srTSRWRadStructAccessData* pRadAccessData)
//{
//	if(pRadAccessData == 0) return;
//}

//*************************************************************************
