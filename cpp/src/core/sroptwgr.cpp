/************************************************************************//**
 * File: sroptwgr.cpp
 * Description: Optical element: Rectangular Waveguide with perfectly conducting walls
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptwgr.h"
#include "gmfft.h"

//*************************************************************************

int srTWaveguideRect::PrepareWavefrontForPropagation(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pOutWfr)
{
	int result;

	long NewNx = pRadAccessData->nx;
	long NewNz = pRadAccessData->nz;
	double NewStepX = pRadAccessData->xStep;
	double NewStepZ = pRadAccessData->zStep;
	double NewStartX = pRadAccessData->xStart;
	double NewStartZ = pRadAccessData->zStart;

	CGenMathFFT aFFT;

	if(!BufVars.PropagInFreeSpaceHoriz)
	{
		NewNx = (long)(Dx/pRadAccessData->xStep);
		aFFT.NextCorrectNumberForFFT(NewNx);
		NewNx = ((NewNx >> 1) << 1); //ensure even
		NewNx <<= 1;
		NewStepX = 2.*Dx/NewNx;
		NewStartX = TransvCenPoint.x - Dx;
	}
	if(!BufVars.PropagInFreeSpaceVert)
	{
		NewNz = (long)(Dz/pRadAccessData->zStep);
		aFFT.NextCorrectNumberForFFT(NewNz);
		NewNz = ((NewNz >> 1) << 1); //ensure even
		NewNz <<= 1;
		NewStepZ = 2.*Dz/NewNz;
		NewStartZ= TransvCenPoint.y - Dz;
	}

	pOutWfr->eStep = pRadAccessData->eStep;
	pOutWfr->eStart = pRadAccessData->eStart;
	pOutWfr->ne = pRadAccessData->ne;
	pOutWfr->xStep = NewStepX;
	pOutWfr->xStart = NewStartX;
	pOutWfr->nx = NewNx;
	pOutWfr->zStep = NewStepZ;
	pOutWfr->zStart = NewStartZ;
	pOutWfr->nz = NewNz;
	pOutWfr->pBaseRadX = 0;
	pOutWfr->pBaseRadZ = 0;

	//long TotData = (pRadAccessData->ne)*NewNx*(NewNz<<1);
	long long TotData = (pRadAccessData->ne)*NewNx*(NewNz<<1);
	if(pRadAccessData->pBaseRadX != 0)
	{
		pOutWfr->pBaseRadX = new float[TotData];
		if(pOutWfr->pBaseRadX == 0) return MEMORY_ALLOCATION_FAILURE;
		float *tBaseRadX = pOutWfr->pBaseRadX;
		for(long j=0; j<TotData; j++) *(tBaseRadX++) = 0.;
	}
	if(pRadAccessData->pBaseRadZ != 0)
	{
		pOutWfr->pBaseRadZ = new float[TotData];
		if(pOutWfr->pBaseRadZ == 0) return MEMORY_ALLOCATION_FAILURE;
		float *tBaseRadZ = pOutWfr->pBaseRadZ;
		for(long j=0; j<TotData; j++) *(tBaseRadZ++) = 0.;
	}

	srTRadResize AuxResize; // approximate values
	AuxResize.pxd = (pOutWfr->xStep)/(pRadAccessData->xStep);
	AuxResize.pzd = (pOutWfr->zStep)/(pRadAccessData->zStep);
	AuxResize.pxm = ((pOutWfr->xStep)*(pOutWfr->nx))/((pRadAccessData->xStep)*(pRadAccessData->nx));
	AuxResize.pzm = ((pOutWfr->zStep)*(pOutWfr->nz))/((pRadAccessData->zStep)*(pRadAccessData->nz));

//set limits for interpolation
	pOutWfr->AuxLong1 = (pOutWfr->nx) >> 2;
	pOutWfr->AuxLong2 = pOutWfr->AuxLong1 + ((pOutWfr->nx) >> 1);
	pOutWfr->AuxLong3 = (pOutWfr->nz) >> 2;
	pOutWfr->AuxLong4 = pOutWfr->AuxLong3 + ((pOutWfr->nz) >> 1);

	double xInterpStart = pOutWfr->xStart + (pOutWfr->AuxLong1)*(pOutWfr->xStep);
	if(xInterpStart <= pRadAccessData->xStart)
	{
		xInterpStart = pRadAccessData->xStart + 0.0001*pRadAccessData->xStep;
		pOutWfr->AuxLong1 = (long)((xInterpStart - pOutWfr->xStart)/(pOutWfr->xStep));
	}

	double xInterpEnd = pOutWfr->xStart + (pOutWfr->AuxLong2)*(pOutWfr->xStep);
	double xEndOld = pRadAccessData->xStart + (pRadAccessData->xStep)*(pRadAccessData->nx);
	if(xInterpEnd >= xEndOld)
	{
		xInterpEnd = xEndOld - 0.0001*pRadAccessData->xStep;
		pOutWfr->AuxLong2 = (long)((xInterpEnd - pOutWfr->xStart)/(pOutWfr->xStep));
	}

	double zInterpStart = pOutWfr->zStart + (pOutWfr->AuxLong3)*(pOutWfr->zStep);
	if(zInterpStart <= pRadAccessData->zStart)
	{
		zInterpStart = pRadAccessData->zStart + 0.0001*pRadAccessData->zStep;
		pOutWfr->AuxLong3 = (long)((zInterpStart - pOutWfr->zStart)/(pOutWfr->zStep));
	}

	double zInterpEnd = pOutWfr->zStart + (pOutWfr->AuxLong4)*(pOutWfr->zStep);
	double zEndOld = pRadAccessData->zStart + (pRadAccessData->zStep)*(pRadAccessData->nz);
	if(zInterpEnd >= zEndOld)
	{
		zInterpEnd = zEndOld - 0.0001*pRadAccessData->zStep;
		pOutWfr->AuxLong4 = (long)((zInterpEnd - pOutWfr->zStart)/(pOutWfr->zStep));
	}

	//this resizes and fills in central part of data (inside the waveguide)
	if(result = RadResizeCore(*pRadAccessData, *pOutWfr, AuxResize, 0)) return result;

	//this fills in external parts of data (outside the waveguide)
	if(result = FillInSymmetricalPartsOutsideWaveguide(*pOutWfr)) return result;

	//todo: treat large aperture size as propagation in free space
		//PropagInFreeSpaceHoriz = false;
		//PropagInFreeSpaceVert = false;

	return 0;
}

//*************************************************************************

int srTWaveguideRect::FillInSymmetricalPartsOutsideWaveguide(srTSRWRadStructAccessData& Wfr)
{
	//long PerX = Wfr.ne << 1;
	//long PerZ = PerX*Wfr.nx;
	long long PerX = Wfr.ne << 1;
	long long PerZ = PerX*Wfr.nx;

	long QuarterNz = Wfr.nz >> 2;
	long QuarterNx = Wfr.nx >> 2;
	long HalfNz = Wfr.nz >> 1;
	long HalfNx = Wfr.nx >> 1;
	long ThreeQuarterNz = 3*QuarterNz;
	long ThreeQuarterNx = 3*QuarterNx;

	long iz, ix;//, OffsetAux;
	long long OffsetAux;

	float *pEx0 = Wfr.pBaseRadX;
	float *pEz0 = Wfr.pBaseRadZ;
	float *pEx, *pEz;

	for(long ie=0; ie<Wfr.ne; ie++)
	{
		long Two_ie = ie << 1;

		for(iz=0; iz<=QuarterNz; iz++)
		{
			//long izPerZ = iz*PerZ;
			//long izOrigPerZ = (HalfNz - iz)*PerZ;
			//long izPerZ_p_Two_ie = izPerZ + Two_ie;
			long long izPerZ = iz*PerZ;
			long long izOrigPerZ = (HalfNz - iz)*PerZ;
			long long izPerZ_p_Two_ie = izPerZ + Two_ie;

			pEx = pEx0 + izPerZ_p_Two_ie;
			pEz = pEz0 + izPerZ_p_Two_ie;
			OffsetAux = izOrigPerZ + Two_ie + HalfNx*PerX;
			float *pExOrig = pEx0 + OffsetAux;
			float *pEzOrig = pEz0 + OffsetAux;
			for(ix=0; ix<=QuarterNx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}

			//ix = QuarterNx // set zero
			OffsetAux = izPerZ + Two_ie + QuarterNx*PerX;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			OffsetAux = izPerZ_p_Two_ie + (QuarterNx + 1)*PerX;
			pEx = pEx0 + OffsetAux;
			pEz = pEz0 + OffsetAux;
			OffsetAux = izOrigPerZ + Two_ie + (QuarterNx + 1)*PerX;
			pExOrig = pEx0 + OffsetAux;
			pEzOrig = pEz0 + OffsetAux;
			for(ix=(QuarterNx + 1); ix<=ThreeQuarterNx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = (*pEzOrig); *(pEz + 1) = (*(pEzOrig + 1));
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig += PerX; pEzOrig += PerX;
			}

			//ix = ThreeQuarterNx // set zero
			OffsetAux = izPerZ + ThreeQuarterNx*PerX + Two_ie;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			OffsetAux = izPerZ_p_Two_ie + (ThreeQuarterNx + 1)*PerX;
			pEx = pEx0 + OffsetAux;
			pEz = pEz0 + OffsetAux;
			OffsetAux = izOrigPerZ + Two_ie + (ThreeQuarterNx - 1)*PerX;
			pExOrig = pEx0 + OffsetAux;
			pEzOrig = pEz0 + OffsetAux;
			for(ix=(ThreeQuarterNx + 1); ix<Wfr.nx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}
		}

		//iz = QuarterNz // set zero
		OffsetAux = QuarterNz*PerZ + Two_ie;
		pEx = pEx0 + OffsetAux;
		for(ix=0; ix<Wfr.nx; ix++)
		{
			*pEx = 0.; *(pEx + 1) = 0.;
			pEx += PerX;
		}

		for(iz=(QuarterNz + 1); iz<=ThreeQuarterNz; iz++)
		{
			//long izPerZ = iz*PerZ;
			//long izPerZ_p_Two_ie = izPerZ + Two_ie;
			long long izPerZ = iz*PerZ;
			long long izPerZ_p_Two_ie = izPerZ + Two_ie;

			pEx = pEx0 + izPerZ_p_Two_ie;
			pEz = pEz0 + izPerZ_p_Two_ie;
			OffsetAux = izPerZ_p_Two_ie + HalfNx*PerX;
			float *pExOrig = pEx0 + OffsetAux;
			float *pEzOrig = pEz0 + OffsetAux;
			for(ix=0; ix<=QuarterNx; ix++)
			{
				*pEx = (*pExOrig); *(pEx + 1) = (*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}

			//ix = QuarterNx // set zero
			OffsetAux = izPerZ + Two_ie + QuarterNx*PerX;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			//ix = 3*QuarterNx // set zero
			OffsetAux = izPerZ + Two_ie + ThreeQuarterNx*PerX;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			OffsetAux = izPerZ_p_Two_ie + (ThreeQuarterNx + 1)*PerX;
			pEx = pEx0 + OffsetAux;
			pEz = pEz0 + OffsetAux;
			OffsetAux = izPerZ_p_Two_ie + (ThreeQuarterNx - 1)*PerX;
			pExOrig = pEx0 + OffsetAux;
			pEzOrig = pEz0 + OffsetAux;
			for(ix=(ThreeQuarterNx + 1); ix<Wfr.nx; ix++)
			{
				*pEx = (*pExOrig); *(pEx + 1) = (*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}
		}

		//iz = 3*QuarterNz // set zero
		OffsetAux = ThreeQuarterNz*PerZ + Two_ie;
		pEx = pEx0 + OffsetAux;
		for(ix=0; ix<Wfr.nx; ix++)
		{
			*pEx = 0.; *(pEx + 1) = 0.;
			pEx += PerX;
		}

		for(iz=(ThreeQuarterNz + 1); iz<Wfr.nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			//long izOrigPerZ = ((ThreeQuarterNz << 1) - iz)*PerZ;
			//long izPerZ_p_Two_ie = izPerZ + Two_ie;
			long long izPerZ = iz*PerZ;
			long long izOrigPerZ = ((ThreeQuarterNz << 1) - iz)*PerZ;
			long long izPerZ_p_Two_ie = izPerZ + Two_ie;

			pEx = pEx0 + izPerZ_p_Two_ie;
			pEz = pEz0 + izPerZ_p_Two_ie;
			OffsetAux = izOrigPerZ + Two_ie + HalfNx*PerX;
			float *pExOrig = pEx0 + OffsetAux;
			float *pEzOrig = pEz0 + OffsetAux;
			for(ix=0; ix<=QuarterNx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}

			//ix = QuarterNx // set zero
			OffsetAux = izPerZ + Two_ie + QuarterNx*PerX;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			OffsetAux = izPerZ_p_Two_ie + (QuarterNx + 1)*PerX;
			pEx = pEx0 + OffsetAux;
			pEz = pEz0 + OffsetAux;
			OffsetAux = izOrigPerZ + Two_ie + (QuarterNx + 1)*PerX;
			pExOrig = pEx0 + OffsetAux;
			pEzOrig = pEz0 + OffsetAux;
			for(ix=(QuarterNx + 1); ix<=ThreeQuarterNx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = (*pEzOrig); *(pEz + 1) = (*(pEzOrig + 1));
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig += PerX; pEzOrig += PerX;
			}

			//ix = ThreeQuarterNx // set zero
			OffsetAux = izPerZ + ThreeQuarterNx*PerX + Two_ie;
			*(pEz0 + OffsetAux) = 0.; *(pEz0 + OffsetAux + 1) = 0.;

			OffsetAux = izPerZ_p_Two_ie + (ThreeQuarterNx + 1)*PerX;
			pEx = pEx0 + OffsetAux;
			pEz = pEz0 + OffsetAux;
			OffsetAux = izOrigPerZ + Two_ie + (ThreeQuarterNx - 1)*PerX;
			pExOrig = pEx0 + OffsetAux;
			pEzOrig = pEz0 + OffsetAux;
			for(ix=(ThreeQuarterNx + 1); ix<Wfr.nx; ix++)
			{
				*pEx = -(*pExOrig); *(pEx + 1) = -(*(pExOrig + 1)); 
				*pEz = -(*pEzOrig); *(pEz + 1) = -(*(pEzOrig + 1)); 
				pEx += PerX; pEz += PerX; // it's not an error !
				pExOrig -= PerX; pEzOrig -= PerX;
			}
		}
	}
	return 0;
}

//*************************************************************************

int srTWaveguideRect::CopyElecFieldDataForOut(srTSRWRadStructAccessData& WfrIn, srTSRWRadStructAccessData& WfrOut)
{
	//int result = 0;
	WfrOut.nx = WfrIn.nx;
	WfrOut.xStart = WfrIn.xStart;
	WfrOut.xStep = WfrIn.xStep;
	WfrOut.nz = WfrIn.nz;
	WfrOut.zStart = WfrIn.zStart;
	WfrOut.zStep = WfrIn.zStep;
	//srTSend LocSend;
	//LocSend.ModifyRadNeNxNz(WfrOut);
	WfrOut.ModifyWfrNeNxNz();

	float *tx = WfrOut.pBaseRadX, *tz = WfrOut.pBaseRadZ;
	float *t1x = WfrIn.pBaseRadX, *t1z = WfrIn.pBaseRadZ;
	//for(long j=0; j<((WfrOut.nx*WfrOut.nz) << 1); j++)
	for(long long j=0; j<((WfrOut.nx*WfrOut.nz) << 1); j++)
	{
		*(tx++) = *(t1x++); *(tz++) = *(t1z++); 
	}
	//LocSend.FinishWorkingWithSRWRadStruct(&AuxStruct);
	return 0;
}

//*************************************************************************
/*
int srTDriftSpace::AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, float** ax, float** az, srTMomentsRatios* MomRatArray)
{// There is a general function like this
	float bxStr0[] = { ax[0][0]*ax[0][0], 2.*ax[0][0]*ax[0][1], ax[0][1]*ax[0][1] };
	float bxStr1[] = { ax[0][0]*ax[1][0], ax[0][1]*ax[1][0] + ax[0][0]*ax[1][1], ax[0][1]*ax[1][1] };
	float bxStr2[] = { ax[1][0]*ax[1][0], 2.*ax[1][0]*ax[1][1], ax[1][1]*ax[1][1] };
	float* bx[] = { bxStr0, bxStr1, bxStr2 };

	float bzStr0[] = { az[0][0]*az[0][0], 2.*az[0][0]*az[0][1], az[0][1]*az[0][1] };
	float bzStr1[] = { az[0][0]*az[1][0], az[0][1]*az[1][0] + az[0][0]*az[1][1], az[0][1]*az[1][1] };
	float bzStr2[] = { az[1][0]*az[1][0], 2.*az[1][0]*az[1][1], az[1][1]*az[1][1] };
	float* bz[] = { bzStr0, bzStr1, bzStr2 };

	float v3dAux1[3], v3dAux2[3];
	float v2dAux1[2], v2dAux2[2];

	double SigXe2, MinSigX, MinSigXe2, SigZe2, MinSigZ, MinSigZe2;

	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	int OffsetE = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + OffsetE);
		if(*(MomX.pTotPhot) != 0.)
		{
			float OldMomX_xx = *(MomX.pXX), OldMomX_xpxp = *(MomX.pXPXP);
			float OldMomX_zz = *(MomX.pZZ), OldMomX_zpzp = *(MomX.pZPZP);
			
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
			if(SigXe2 < MinSigXe2) *(MomX.pXX) = MinSigXe2 + (*(MomX.pX))*(*(MomX.pX));
			
			SigZe2 = (*(MomX.pZZ)) - (*(MomX.pZ))*(*(MomX.pZ));
			MinSigZ = DiffractionLimitedPropagatedSpotSize('z', pRadAccessData, ie);
			MinSigZe2 = MinSigZ*MinSigZ;
			if(SigZe2 < MinSigZe2) *(MomX.pZZ) = MinSigZe2 + (*(MomX.pZ))*(*(MomX.pZ));
			
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomX = (*(MomX.pXX) > 0.)? sqrt(*(MomX.pXX)/OldMomX_xx) : -1.;
				tMomRatArray->RxpxpMomX = (*(MomX.pXPXP) > 0.)? sqrt(*(MomX.pXPXP)/OldMomX_xpxp) : -1.;
				tMomRatArray->RzzMomX = (*(MomX.pZZ) > 0.)? sqrt(*(MomX.pZZ)/OldMomX_zz) : -1.;
				tMomRatArray->RzpzpMomX = (*(MomX.pZPZP) > 0.)? sqrt(*(MomX.pZPZP)/OldMomX_zpzp) : -1.;
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
			float OldMomZ_xx = *(MomZ.pXX), OldMomZ_xpxp = *(MomZ.pXPXP);
			float OldMomZ_zz = *(MomZ.pZZ), OldMomZ_zpzp = *(MomZ.pZPZP);
			
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
			if(SigXe2 < MinSigXe2) *(MomZ.pXX) = MinSigXe2 + (*(MomZ.pX))*(*(MomZ.pX));
			
			SigZe2 = (*(MomZ.pZZ)) - (*(MomZ.pZ))*(*(MomZ.pZ));
			if(SigZe2 < MinSigZe2) *(MomZ.pZZ) = MinSigZe2 + (*(MomZ.pZ))*(*(MomZ.pZ));
			
			if(FillMomRatios)
			{
				tMomRatArray->RxxMomZ = (*(MomZ.pXX) > 0.)? sqrt(*(MomZ.pXX)/OldMomZ_xx) : -1.;
				tMomRatArray->RxpxpMomZ = (*(MomZ.pXPXP) > 0.)? sqrt(*(MomZ.pXPXP)/OldMomZ_xpxp) : -1.;
				tMomRatArray->RzzMomZ = (*(MomZ.pZZ) > 0.)? sqrt(*(MomZ.pZZ)/OldMomZ_zz) : -1.;
				tMomRatArray->RzpzpMomZ = (*(MomZ.pZPZP) > 0.)? sqrt(*(MomZ.pZPZP)/OldMomZ_zpzp) : -1.;
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
*/
//*************************************************************************
/*
int srTDriftSpace::TuneRadForPropMeth_1(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& PostResize)
{
	srTMomentsRatios* MomRatArray = new srTMomentsRatios[pRadAccessData->ne];
	if(MomRatArray == 0) return MEMORY_ALLOCATION_FAILURE;

	int result;
	if(pRadAccessData->Pres != 0) // Go to spatial...
		if(result = SetRadRepres(pRadAccessData, 0)) return result;

	if(result = PropagateRadMoments(pRadAccessData, MomRatArray)) return result;
	
	srTMomentsRatios* tMomRatArray = MomRatArray;

	float pxMaxMomX = (float)(1.e-23), pxMinMomX = (float)(1.e+23), pzMaxMomX = (float)(1.e-23), pzMinMomX = (float)(1.e+23);
	float pxMaxMomZ = (float)(1.e-23), pxMinMomZ = (float)(1.e+23), pzMaxMomZ = (float)(1.e-23), pzMinMomZ = (float)(1.e+23);

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

	float pxMax = (pxMaxMomX > pxMaxMomZ)? pxMaxMomX : pxMaxMomZ;
	float pxMin = (pxMinMomX < pxMinMomZ)? pxMinMomX : pxMinMomZ;
	float pzMax = (pzMaxMomX > pzMaxMomZ)? pzMaxMomX : pzMaxMomZ;
	float pzMin = (pzMinMomX < pzMinMomZ)? pzMinMomX : pzMinMomZ;

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
*/
//*************************************************************************
/*
int srTDriftSpace::PropagateRadiationMeth_1(srTSRWRadStructAccessData* pRadAccessData)
{
	int result;
	srTRadResize PostResize;
	PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;

	float* OldMxxArr = new float[pRadAccessData->ne];
	if(OldMxxArr == 0) return MEMORY_ALLOCATION_FAILURE;
	float* OldMzzArr = new float[pRadAccessData->ne];
	if(OldMzzArr == 0) return MEMORY_ALLOCATION_FAILURE;
	SetupMxxMzzArr(pRadAccessData, OldMxxArr, OldMzzArr);

	float *NewMxxArr = 0, *NewMzzArr = 0;

	if(result = TuneRadForPropMeth_1(pRadAccessData, PostResize)) return result;
	if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;

	if(pRadAccessData->Pres != 1) if(result = SetRadRepres(pRadAccessData, 1)) return result;
	if(result = TraverseRadZXE(pRadAccessData)) return result;
	if(result = SetRadRepres(pRadAccessData, 0)) return result;

	const double ResizeTol = 0.15;
	if((PostResize.pxm != -1) && (PostResize.pzm != -1))
	{
		char PostResizeNeeded = (::fabs(PostResize.pxm - 1.) || ::fabs(PostResize.pzm - 1.));
		if(PostResizeNeeded) if(result = RadResizeGen(*pRadAccessData, PostResize)) return result;
	}
	else
	{
		if(result = ComputeRadMoments(pRadAccessData)) return result;

		NewMxxArr = new float[pRadAccessData->ne];
		if(NewMxxArr == 0) return MEMORY_ALLOCATION_FAILURE;
		NewMzzArr = new float[pRadAccessData->ne];
		if(NewMzzArr == 0) return MEMORY_ALLOCATION_FAILURE;
		SetupMxxMzzArr(pRadAccessData, NewMxxArr, NewMzzArr);

		float pxmMinE2, pxmMaxE2, pzmMinE2, pzmMaxE2;
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
*/
//*************************************************************************
/*
int srTDriftSpace::PropagateRadiationSimple_PropToWaist(srTSRWRadStructAccessData* pRadAccessData)
{// e in eV; Length in m !!!
	int result;

	SetupPropBufVars_PropToWaist(pRadAccessData);
	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;

	PropBufVars.PassNo = 1;
	if(result = TraverseRadZXE(pRadAccessData)) return result;

	if(result = ResizeBeforePropToWaistIfNecessary(pRadAccessData)) return result;

	double InvLambda_m = pRadAccessData->eStart*806546.577258;
	double InvLambda_m_d_Length = InvLambda_m/Length;
	double LambdaM_Length = 1./InvLambda_m_d_Length;

	srTFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.UseGivenStartTrValues = 0;

	srTFFT2D FFT2D;

	srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
	if(result = SetupWfrEdgeCorrData(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;

	FFT2DInfo.pData = pRadAccessData->pBaseRadX;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	
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
*/
//*************************************************************************
/*
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

	srTFFT2DInfo FFT2DInfo;
	FFT2DInfo.xStep = pRadAccessData->xStep;
	FFT2DInfo.yStep = pRadAccessData->zStep;
	FFT2DInfo.xStart = pRadAccessData->xStart;
	FFT2DInfo.yStart = pRadAccessData->zStart;
	FFT2DInfo.Nx = pRadAccessData->nx;
	FFT2DInfo.Ny = pRadAccessData->nz;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.UseGivenStartTrValues = 0;

	srTDataPtrsForWfrEdgeCorr DataPtrsForWfrEdgeCorr;
	if(result = SetupWfrEdgeCorrData(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr)) return result;

	srTFFT2D FFT2D;
	FFT2DInfo.pData = pRadAccessData->pBaseRadX;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	FFT2DInfo.pData = pRadAccessData->pBaseRadZ;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

	if(DataPtrsForWfrEdgeCorr.WasSetup)
	{
		MakeWfrEdgeCorrection(pRadAccessData, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, DataPtrsForWfrEdgeCorr);
		DataPtrsForWfrEdgeCorr.DisposeData();
	}

// Re-scaling
	pRadAccessData->xStart = (FFT2DInfo.xStartTr)*LambdaM_Length;
	pRadAccessData->zStart = (FFT2DInfo.yStartTr)*LambdaM_Length;
	pRadAccessData->xStep = FFT2DInfo.xStepTr*LambdaM_Length;
	pRadAccessData->zStep = FFT2DInfo.yStepTr*LambdaM_Length;

	PropBufVars.PassNo = 2;
	if(result = TraverseRadZXE(pRadAccessData)) return result;
	return result;
}
/*
//*************************************************************************
/*
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

	long TwoNp = pSect1D->np << 1;
	float *AuxCont = new float[TwoNp << 1];
	if(AuxCont == 0) return MEMORY_ALLOCATION_FAILURE;

	float *pAuxX = AuxCont, *pAuxZ = AuxCont + TwoNp;
	float *tAuxX = pAuxX, *tAuxZ = pAuxZ;
	float *tEx = pSect1D->pEx, *tEz = pSect1D->pEz;
	for(long i=0; i<TwoNp; i++)
	{
		*(tAuxX++) = *(tEx++); 
		*(tAuxZ++) = *(tEz++);
	}

	srTFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = pSect1D->ArgStep;
	FFT1DInfo.xStart = pSect1D->ArgStart;
	FFT1DInfo.Nx = pSect1D->np;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
	srTFFT1D FFT1D;

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
*/
//*************************************************************************
/*
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

	double pxd = pRadAccessData->pResAfter->pxm*xRange*pRadAccessData->xStep/(LambdaM_Length);
	if(::fabs(pxd - 1.) < DiffAllowResize) pxd = 1.;

	if((pRadAccessData->nx)*pxd < MinNxInRange)
	{
		pxd = MinNxInRange/double(pRadAccessData->nx);
		if(::fabs(pxd - 1.) < DiffAllowResize) pxd = 1.;
		double xRangeNew = pxd*LambdaM_Length/pRadAccessData->xStep;

		double xRangeShouldBe = xRange*pRadAccessData->pResAfter->pxm;
		pRadAccessData->pResAfter->pxm = xRangeShouldBe/xRangeNew;
		if(pRadAccessData->pResAfter->pxm > 1.) pRadAccessData->pResAfter->pxm = 1.;
	}
	else
	{
		pRadAccessData->pResAfter->pxm = 1.;
	}

	double pzd = pRadAccessData->pResAfter->pzm*zRange*pRadAccessData->zStep/(LambdaM_Length);
	if(::fabs(pzd - 1.) < DiffAllowResize) pzd = 1.;

	if((pRadAccessData->nz)*pzd < MinNzInRange)
	{
		pzd = MinNzInRange/double(pRadAccessData->nz);
		if(::fabs(pzd - 1.) < DiffAllowResize) pzd = 1.;

		double zRangeNew = pzd*LambdaM_Length/pRadAccessData->zStep;
		double zRangeShouldBe = zRange*pRadAccessData->pResAfter->pzm;
		pRadAccessData->pResAfter->pzm = zRangeShouldBe/zRangeNew;
		if(pRadAccessData->pResAfter->pzm > 1.) pRadAccessData->pResAfter->pzm = 1.;
	}
	else
	{
		pRadAccessData->pResAfter->pzm = 1.;
	}

	if((::fabs(pxm - 1.) > DiffAllowResize) || (::fabs(pxd - 1.) > DiffAllowResize) || (::fabs(pzm - 1.) > DiffAllowResize) || (::fabs(pzd - 1.) > DiffAllowResize))
	{
		srTRadResize Resize;
		Resize.pxd = pxd; Resize.pxm = pxm; Resize.pzd = pzd; Resize.pzm = pzm;
		Resize.DoNotTreatSpherTerm = 1;

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
		if(::fabs(Length + RxOld) > SmallLength) Resize.DoNotTreatSpherTerm = 0;
		if(::fabs(Length + RzOld) > SmallLength) Resize.DoNotTreatSpherTerm = 0;

		if(!(Resize.DoNotTreatSpherTerm))
		{
			pRadAccessData->RobsX = Length*RxOld/(Length + RxOld);
			pRadAccessData->RobsZ = Length*RzOld/(Length + RzOld);
		}

		if(result = RadResizeGen(*pRadAccessData, Resize)) return result;

		pRadAccessData->RobsX = RxOld;
		pRadAccessData->RobsZ = RzOld;

		if(ResolutionWasReduced)
		{
			srTSend Send; Send.AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
		}
	}

	pRadAccessData->pResAfter->RelCenPosX = 0.5; // To check
	pRadAccessData->pResAfter->RelCenPosZ = 0.5;

	return 0;
}
*/
//*************************************************************************
/*
void srTDriftSpace::EstimateMinNxNzBeforePropToWaist(srTSRWRadStructAccessData* pRadAccessData, int& Nx, int& Nz)
{
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
*/
//*************************************************************************
/*
void srTDriftSpace::EstimateMinNpBeforePropToWaist1D(srTRadSect1D* pSect1D, int& Np)
{
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
	double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/1.3; // To steer
	Np = int(::fabs(xEndRel - xStartRel)/dx) + 1;
	if(((Np >> 1) << 1) != Np) Np++;
}
*/
//*************************************************************************
/*
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
*/
//*************************************************************************
