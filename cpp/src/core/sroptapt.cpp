/************************************************************************//**
 * File: sroptapt.cpp
 * Description: Optical element: Aperture
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptapt.h"

//*************************************************************************

int srTAperture::TraverseRadZXE_TracingZeroField(srTSRWRadStructAccessData* pRadAccessData)
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
	long long izPerZ = 0;

	int ixStartNonZeroOld = pRadAccessData->nx + 1, izStartNonZeroOld = -1;
	int ixStartZeroAgainOld = -1, izStartZeroAgainOld = -1;
	int ixStartNonZeroNew = pRadAccessData->nx + 1, izStartNonZeroNew = -1;
	int ixStartZeroAgainNew = -1, izStartZeroAgainNew = -1;

	char zPrevOldWasZero = 1, zPrevNewWasZero = 1;

	for(int iz=0; iz<pRadAccessData->nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		EXZ.x = pRadAccessData->xStart;
		//long ixPerX = 0;
		long long ixPerX = 0;

		char xPrevOldWasZero = 1, xPrevNewWasZero = 1;
		int ixStartNonZeroOldLoc = -1, ixStartZeroAgainOldLoc = -1;
		int ixStartNonZeroNewLoc = -1, ixStartZeroAgainNewLoc = -1;

		for(int ix=0; ix<pRadAccessData->nx; ix++)
		{
			float *pEx_StartForE = pEx_StartForX + ixPerX;
			float *pEz_StartForE = pEz_StartForX + ixPerX;
			EXZ.e = pRadAccessData->eStart;
			//long iePerE = 0;
			long long iePerE = 0;

			char EFieldOldWasNonZero = 0, EFieldNewIsNonZero = 0;
			for(int ie=0; ie<pRadAccessData->ne; ie++)
			{
				EFieldPtrs.pExRe = pEx_StartForE + iePerE;
				EFieldPtrs.pExIm = EFieldPtrs.pExRe + 1;
				EFieldPtrs.pEzRe = pEz_StartForE + iePerE;
				EFieldPtrs.pEzIm = EFieldPtrs.pEzRe + 1;

				char EFieldOldWasNonZeroLoc = !((*(EFieldPtrs.pExRe) == 0.) && (*(EFieldPtrs.pExIm) == 0.) && (*(EFieldPtrs.pEzRe) == 0.) && (*(EFieldPtrs.pEzIm) == 0.));
				if(EFieldOldWasNonZeroLoc)
				{
					RadPointModifier(EXZ, EFieldPtrs);

					EFieldOldWasNonZero = 1;
					EFieldNewIsNonZero = !((*(EFieldPtrs.pExRe) == 0.) && (*(EFieldPtrs.pExIm) == 0.) && (*(EFieldPtrs.pEzRe) == 0.) && (*(EFieldPtrs.pEzIm) == 0.));
				}

				iePerE += 2;
				EXZ.e += pRadAccessData->eStep;
			}

			if(EFieldOldWasNonZero && xPrevOldWasZero)
			{
				ixStartNonZeroOldLoc = ix; xPrevOldWasZero = 0;
			}
			if((!EFieldOldWasNonZero) && (!xPrevOldWasZero))
			{
				ixStartZeroAgainOldLoc = ix; xPrevOldWasZero = 1;
			}

			if(EFieldNewIsNonZero && xPrevNewWasZero)
			{
				ixStartNonZeroNewLoc = ix; xPrevNewWasZero = 0;
			}
			if((!EFieldNewIsNonZero) && (!xPrevNewWasZero))
			{
				ixStartZeroAgainNewLoc = ix; xPrevNewWasZero = 1;
			}

			ixPerX += PerX;
			EXZ.x += pRadAccessData->xStep;
		}

		if(ixStartZeroAgainOldLoc == -1) ixStartZeroAgainOldLoc = pRadAccessData->nx;
		if(ixStartZeroAgainNewLoc == -1) ixStartZeroAgainNewLoc = pRadAccessData->nx;

		if((ixStartNonZeroOldLoc != -1) && (ixStartNonZeroOldLoc < ixStartNonZeroOld)) ixStartNonZeroOld = ixStartNonZeroOldLoc;
		if((ixStartZeroAgainOldLoc != -1) && (ixStartZeroAgainOldLoc > ixStartZeroAgainOld)) ixStartZeroAgainOld = ixStartZeroAgainOldLoc;
		if((ixStartNonZeroOldLoc != -1) && zPrevOldWasZero)
		{
			izStartNonZeroOld = iz; zPrevOldWasZero = 0;
		}
		if((ixStartNonZeroOldLoc == -1) && (!zPrevOldWasZero))
		{
			izStartZeroAgainOld = iz; zPrevOldWasZero = 1;
		}

		if((ixStartNonZeroNewLoc != -1) && (ixStartNonZeroNewLoc < ixStartNonZeroNew)) ixStartNonZeroNew = ixStartNonZeroNewLoc;
		if((ixStartZeroAgainNewLoc != -1) && (ixStartZeroAgainNewLoc > ixStartZeroAgainNew)) ixStartZeroAgainNew = ixStartZeroAgainNewLoc;
		if((ixStartNonZeroNewLoc != -1) && zPrevNewWasZero)
		{
			izStartNonZeroNew = iz; zPrevNewWasZero = 0;
		}
		if((ixStartNonZeroNewLoc == -1) && (!zPrevNewWasZero))
		{
			izStartZeroAgainNew = iz; zPrevNewWasZero = 1;
		}

		izPerZ += PerZ;
		EXZ.z += pRadAccessData->zStep;
	}

	if(izStartZeroAgainOld == -1) izStartZeroAgainOld = pRadAccessData->nz;
	if(izStartZeroAgainNew == -1) izStartZeroAgainNew = pRadAccessData->nz;

	if((ixStartNonZeroOld < 0) && (ixStartNonZeroOld >= pRadAccessData->nx)) ixStartNonZeroOld = 0;
	if((ixStartZeroAgainOld < 0) && (ixStartZeroAgainOld >= pRadAccessData->nx)) ixStartZeroAgainOld = 0;
	if((ixStartNonZeroNew < 0) && (ixStartNonZeroNew >= pRadAccessData->nx)) ixStartNonZeroNew = 0;
	if((ixStartZeroAgainNew < 0) && (ixStartZeroAgainNew >= pRadAccessData->nx)) ixStartZeroAgainNew = 0;
	if((izStartNonZeroOld < 0) && (izStartNonZeroOld >= pRadAccessData->nz)) izStartNonZeroOld = 0;
	if((izStartZeroAgainOld < 0) && (izStartZeroAgainOld >= pRadAccessData->nz)) izStartZeroAgainOld = 0;
	if((izStartNonZeroNew < 0) && (izStartNonZeroNew >= pRadAccessData->nz)) izStartNonZeroNew = 0;
	if((izStartZeroAgainNew < 0) && (izStartZeroAgainNew >= pRadAccessData->nz)) izStartZeroAgainNew = 0;

	xStartNonZeroOld = pRadAccessData->xStart + ixStartNonZeroOld*pRadAccessData->xStep;
	zStartNonZeroOld = pRadAccessData->zStart + izStartNonZeroOld*pRadAccessData->zStep;
	xStartZeroAgainOld = pRadAccessData->xStart + ixStartZeroAgainOld*pRadAccessData->xStep;
	zStartZeroAgainOld = pRadAccessData->zStart + izStartZeroAgainOld*pRadAccessData->zStep;
	xStartNonZeroNew = pRadAccessData->xStart + ixStartNonZeroNew*pRadAccessData->xStep;
	zStartNonZeroNew = pRadAccessData->zStart + izStartNonZeroNew*pRadAccessData->zStep;
	xStartZeroAgainNew = pRadAccessData->xStart + ixStartZeroAgainNew*pRadAccessData->xStep;
	zStartZeroAgainNew = pRadAccessData->zStart + izStartZeroAgainNew*pRadAccessData->zStep;
	return 0;
}

//*************************************************************************

int srTAperture::TuneRadAfterPropMeth_1(srTSRWRadStructAccessData* pRadAccessData)
{
	if((xStartNonZeroOld == xStartZeroAgainOld) || (xStartNonZeroNew == xStartZeroAgainNew)) return 0;
	if((zStartNonZeroOld == zStartZeroAgainOld) || (zStartNonZeroNew == zStartZeroAgainNew)) return 0;

	double XrNonZeroLeftOld = xStartNonZeroOld - pRadAccessData->xc;
	double XrNonZeroLeftNew = xStartNonZeroNew - pRadAccessData->xc;
	double XrNonZeroRightOld = xStartZeroAgainOld - pRadAccessData->xc;
	double XrNonZeroRightNew = xStartZeroAgainNew - pRadAccessData->xc;

	double ZrNonZeroLeftOld = zStartNonZeroOld - pRadAccessData->zc;
	double ZrNonZeroLeftNew = zStartNonZeroNew - pRadAccessData->zc;
	double ZrNonZeroRightOld = zStartZeroAgainOld - pRadAccessData->zc;
	double ZrNonZeroRightNew = zStartZeroAgainNew - pRadAccessData->zc;

	//double Lambda_m = 1.239854e-06/pRadAccessData->eStart;
	double Lambda_m = 1.239842e-06/pRadAccessData->eStart;
	if(pRadAccessData->PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

	double XcNonZeroNew = 0.5*(XrNonZeroLeftNew + XrNonZeroRightNew);
	double XcNonZeroNewE2 = XcNonZeroNew*XcNonZeroNew;

	double HorSizeEst1 = (XrNonZeroRightNew*XrNonZeroRightNew - XcNonZeroNewE2)/(pRadAccessData->RobsX*Lambda_m);
	double HorSizeEst2 = (XcNonZeroNewE2 - XrNonZeroLeftNew*XrNonZeroLeftNew)/(pRadAccessData->RobsX*Lambda_m);
	char HorSizeIsVerySmall = ((::fabs(HorSizeEst1) < 1.) && (::fabs(HorSizeEst2) < 1.));
	if(HorSizeIsVerySmall)
	{
		pRadAccessData->RobsX = 0.; pRadAccessData->xc = 0.5*(xStartNonZeroNew + xStartZeroAgainNew);
	}

	double ZcNonZeroNew = 0.5*(ZrNonZeroLeftNew + ZrNonZeroRightNew);
	double ZcNonZeroNewE2 = ZcNonZeroNew*ZcNonZeroNew;

	double VerSizeEst1 = (ZrNonZeroRightNew*ZrNonZeroRightNew - ZcNonZeroNewE2)/(pRadAccessData->RobsZ*Lambda_m);
	double VerSizeEst2 = (ZcNonZeroNewE2 - ZrNonZeroLeftNew*ZrNonZeroLeftNew)/(pRadAccessData->RobsZ*Lambda_m);
	char VertSizeIsVerySmall = ((::fabs(VerSizeEst1) < 1.) && (::fabs(VerSizeEst2) < 1.));
	if(VertSizeIsVerySmall)
	{
		pRadAccessData->RobsZ = 0.; pRadAccessData->zc = 0.5*(zStartNonZeroNew + zStartZeroAgainNew);
	}

	srTRadResize RadResize;
	RadResize.pxm = RadResize.pxd = RadResize.pzm = RadResize.pzd = 1.;
	double PxLeft, PxRight, PzLeft, PzRight;

	if((XrNonZeroLeftNew >= 0.) || (XrNonZeroLeftOld >= 0.)) PxLeft = -1.;
	else PxLeft = XrNonZeroLeftNew/XrNonZeroLeftOld;
	if((XrNonZeroRightNew <= 0.) || (XrNonZeroRightOld <= 0.)) PxRight = -1.;
	else PxRight = XrNonZeroRightNew/XrNonZeroRightOld;
	if((ZrNonZeroLeftNew >= 0.) || (ZrNonZeroLeftOld >= 0.)) PzLeft = -1.;
	else PzLeft = ZrNonZeroLeftNew/ZrNonZeroLeftOld;
	if((ZrNonZeroRightNew <= 0.) || (ZrNonZeroRightOld <= 0.)) PzRight = -1.;
	else PzRight = ZrNonZeroRightNew/ZrNonZeroRightOld;

	int result;
	double pxMax = (PxLeft > PxRight)? PxLeft : PxRight;
	double pzMax = (PzLeft > PzRight)? PzLeft : PzRight;
	if((pxMax != -1.) && (pzMax != -1.))
	{
		const double ResizeTol = 0.15;
		const int MinAmOfPoints = 10;

		char xResizeNeeded = (1.- pxMax > ResizeTol);
		char zResizeNeeded = (1.- pzMax > ResizeTol);
		if(xResizeNeeded) 
		{
			int NewAmOfPo = (int)(pxMax*pRadAccessData->nx);
			if(NewAmOfPo < MinAmOfPoints) pxMax = MinAmOfPoints/double(pRadAccessData->nx);
			RadResize.pxm = RadResize.pxd = pxMax;
		}
		if(zResizeNeeded) 
		{ 
			int NewAmOfPo = (int)(pzMax*pRadAccessData->nz);
			if(NewAmOfPo < MinAmOfPoints) pzMax = MinAmOfPoints/double(pRadAccessData->nz);
			RadResize.pzm = RadResize.pzd = pzMax;
		}
		if(xResizeNeeded || zResizeNeeded) if(result = RadResizeGen(*pRadAccessData, RadResize)) return result;
	}
	return 0;
}

//*************************************************************************

int srTAperture::PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
{
	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX = (1./MomX.SqrtMxx);
			tMomRatArray->RxpxpMomX = (1./MomX.SqrtMxpxp);
			tMomRatArray->RzzMomX = (1./MomX.SqrtMzz);
			tMomRatArray->RzpzpMomX = (1./MomX.SqrtMzpzp);
			tMomRatArray->RxxMomZ = (1./MomZ.SqrtMxx);
			tMomRatArray->RxpxpMomZ = (1./MomZ.SqrtMxpxp);
			tMomRatArray->RzzMomZ = (1./MomZ.SqrtMzz);
			tMomRatArray->RzpzpMomZ = (1./MomZ.SqrtMzpzp);

			//Offset += 10;
			Offset += 11; //OC071108
			tMomRatArray++;
		}
		tMomRatArray = MomRatArray;
	}

	int Offset = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

		//float MomX_X = *(MomX.pX), MomX_Z = *(MomX.pZ), MomZ_X = *(MomZ.pX), MomZ_Z = *(MomZ.pZ);
		//float MomX_SqrtMxx = MomX.SqrtMxx, MomX_SqrtMzz = MomX.SqrtMzz, MomZ_SqrtMxx = MomZ.SqrtMxx, MomZ_SqrtMzz = MomZ.SqrtMzz;
		double MomX_X = *(MomX.pX), MomX_Z = *(MomX.pZ), MomZ_X = *(MomZ.pX), MomZ_Z = *(MomZ.pZ); //OC130311
		double MomX_SqrtMxx = MomX.SqrtMxx, MomX_SqrtMzz = MomX.SqrtMzz, MomZ_SqrtMxx = MomZ.SqrtMxx, MomZ_SqrtMzz = MomZ.SqrtMzz;
		if(TransHndl.rep != 0)
		{
			TVector3d PointMomX(MomX_X, 0., MomX_Z);
			FromLabToLocFrame_Point(PointMomX);
			MomX_X = (PointMomX.x + TransvCenPoint.x);
			MomX_Z = (PointMomX.z + TransvCenPoint.y);

			TVector3d VectMomX(MomX_SqrtMxx, 0., MomX_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomX);
			MomX_SqrtMxx = VectMomX.x;
			MomX_SqrtMzz = VectMomX.z;

			TVector3d PointMomZ(MomZ_X, 0., MomZ_Z);
			FromLabToLocFrame_Point(PointMomZ);
			MomZ_X = (PointMomZ.x + TransvCenPoint.x);
			MomZ_Z = (PointMomZ.z + TransvCenPoint.y);

			TVector3d VectMomZ(MomZ_SqrtMxx, 0., MomZ_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomZ);
			MomZ_SqrtMxx = VectMomZ.x;
			MomZ_SqrtMzz = VectMomZ.z;
		}

		//float MomX_SqrtMxx_Mult = (float)(1.7*MomX_SqrtMxx), MomX_SqrtMzz_Mult = (float)(1.7*MomX_SqrtMzz), MomZ_SqrtMxx_Mult = (float)(1.7*MomZ_SqrtMxx), MomZ_SqrtMzz_Mult = (float)(1.7*MomZ_SqrtMzz);
		double MomX_SqrtMxx_Mult = (1.7*MomX_SqrtMxx), MomX_SqrtMzz_Mult = (1.7*MomX_SqrtMzz), MomZ_SqrtMxx_Mult = (1.7*MomZ_SqrtMxx), MomZ_SqrtMzz_Mult = (1.7*MomZ_SqrtMzz); //OC130311

		if(CheckIfMomentsShouldBeRecomputed(MomX_X, MomX_Z, MomZ_X, MomZ_Z, MomX_SqrtMxx_Mult, MomX_SqrtMzz_Mult, MomZ_SqrtMxx_Mult, MomZ_SqrtMzz_Mult))
		{
			int result;
			if(result = ComputeRadMoments(pRadAccessData)) return result;
			break;
		}

		//Offset += 10;
		Offset += 11; //OC071108
	}

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs NewMomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs NewMomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX *= NewMomX.SqrtMxx;
			tMomRatArray->RxpxpMomX *= NewMomX.SqrtMxpxp;
			tMomRatArray->RzzMomX *= NewMomX.SqrtMzz;
			tMomRatArray->RzpzpMomX *= NewMomX.SqrtMzpzp;
			tMomRatArray->RxxMomZ *= NewMomZ.SqrtMxx;
			tMomRatArray->RxpxpMomZ *= NewMomZ.SqrtMxpxp;
			tMomRatArray->RzzMomZ *= NewMomZ.SqrtMzz;
			tMomRatArray->RzpzpMomZ *= NewMomZ.SqrtMzpzp;

			//Offset += 10;
			Offset += 11; //OC071108
			tMomRatArray++;
		}
	}
	return 0;
}

//*************************************************************************
/**
int srTRectAperture::PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
{
	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX = (float)(1./MomX.SqrtMxx);
			tMomRatArray->RxpxpMomX = (float)(1./MomX.SqrtMxpxp);
			tMomRatArray->RzzMomX = (float)(1./MomX.SqrtMzz);
			tMomRatArray->RzpzpMomX = (float)(1./MomX.SqrtMzpzp);
			tMomRatArray->RxxMomZ = (float)(1./MomZ.SqrtMxx);
			tMomRatArray->RxpxpMomZ = (float)(1./MomZ.SqrtMxpxp);
			tMomRatArray->RzzMomZ = (float)(1./MomZ.SqrtMzz);
			tMomRatArray->RzpzpMomZ = (float)(1./MomZ.SqrtMzpzp);

			Offset += 10;
			tMomRatArray++;
		}
		tMomRatArray = MomRatArray;
	}

	int Offset = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

		float MomX_X = *(MomX.pX), MomX_Z = *(MomX.pZ), MomZ_X = *(MomZ.pX), MomZ_Z = *(MomZ.pZ);
		float MomX_SqrtMxx = MomX.SqrtMxx, MomX_SqrtMzz = MomX.SqrtMzz, MomZ_SqrtMxx = MomZ.SqrtMxx, MomZ_SqrtMzz = MomZ.SqrtMzz;
		if(TransHndl.rep != 0)
		{
			TVector3d PointMomX(MomX_X, 0., MomX_Z);
			FromLabToLocFrame_Point(PointMomX);
			MomX_X = (float)(PointMomX.x + TransvCenPoint.x);
			MomX_Z = (float)(PointMomX.z + TransvCenPoint.y);

			TVector3d VectMomX(MomX_SqrtMxx, 0., MomX_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomX);
			MomX_SqrtMxx = (float)VectMomX.x;
			MomX_SqrtMzz = (float)VectMomX.z;

			TVector3d PointMomZ(MomZ_X, 0., MomZ_Z);
			FromLabToLocFrame_Point(PointMomZ);
			MomZ_X = (float)(PointMomZ.x + TransvCenPoint.x);
			MomZ_Z = (float)(PointMomZ.z + TransvCenPoint.y);

			TVector3d VectMomZ(MomZ_SqrtMxx, 0., MomZ_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomZ);
			MomZ_SqrtMxx = (float)VectMomZ.x;
			MomZ_SqrtMzz = (float)VectMomZ.z;
		}

		float MomX_SqrtMxx_Mult = (float)(1.7*MomX_SqrtMxx), MomX_SqrtMzz_Mult = (float)(1.7*MomX_SqrtMzz), MomZ_SqrtMxx_Mult = (float)(1.7*MomZ_SqrtMxx), MomZ_SqrtMzz_Mult = (float)(1.7*MomZ_SqrtMzz);

		if((((MomX_X - MomX_SqrtMxx_Mult) < (TransvCenPoint.x - HalfDx)) || ((MomX_X + MomX_SqrtMxx_Mult) > (TransvCenPoint.x + HalfDx))) ||
		   (((MomX_Z - MomX_SqrtMzz_Mult) < (TransvCenPoint.y - HalfDz)) || ((MomX_Z + MomX_SqrtMzz_Mult) > (TransvCenPoint.y + HalfDz))) ||
		   (((MomZ_X - MomZ_SqrtMxx_Mult) < (TransvCenPoint.x - HalfDx)) || ((MomZ_X + MomZ_SqrtMxx_Mult) > (TransvCenPoint.x + HalfDx))) ||
		   (((MomZ_Z - MomZ_SqrtMzz_Mult) < (TransvCenPoint.y - HalfDz)) || ((MomZ_Z + MomZ_SqrtMzz_Mult) > (TransvCenPoint.y + HalfDz))))
		{
			int result;
			if(result = ComputeRadMoments(pRadAccessData)) return result;
			break;
		}
		Offset += 10;
	}

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs NewMomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs NewMomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX *= NewMomX.SqrtMxx;
			tMomRatArray->RxpxpMomX *= NewMomX.SqrtMxpxp;
			tMomRatArray->RzzMomX *= NewMomX.SqrtMzz;
			tMomRatArray->RzpzpMomX *= NewMomX.SqrtMzpzp;
			tMomRatArray->RxxMomZ *= NewMomZ.SqrtMxx;
			tMomRatArray->RxpxpMomZ *= NewMomZ.SqrtMxpxp;
			tMomRatArray->RzzMomZ *= NewMomZ.SqrtMzz;
			tMomRatArray->RzpzpMomZ *= NewMomZ.SqrtMzpzp;

			Offset += 10;
			tMomRatArray++;
		}
	}
	return 0;
}
**/
//*************************************************************************
/**
int srTCircAperture::PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
{
	char FillMomRatios = (MomRatArray != 0);
	srTMomentsRatios* tMomRatArray = MomRatArray;

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX = (float)(1./MomX.SqrtMxx);
			tMomRatArray->RxpxpMomX = (float)(1./MomX.SqrtMxpxp);
			tMomRatArray->RzzMomX = (float)(1./MomX.SqrtMzz);
			tMomRatArray->RzpzpMomX = (float)(1./MomX.SqrtMzpzp);
			tMomRatArray->RxxMomZ = (float)(1./MomZ.SqrtMxx);
			tMomRatArray->RxpxpMomZ = (float)(1./MomZ.SqrtMxpxp);
			tMomRatArray->RzzMomZ = (float)(1./MomZ.SqrtMzz);
			tMomRatArray->RzpzpMomZ = (float)(1./MomZ.SqrtMzpzp);

			Offset += 10;
			tMomRatArray++;
		}
		tMomRatArray = MomRatArray;
	}

	int Offset = 0;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		srTMomentsPtrs MomX(pRadAccessData->pMomX + Offset);
		srTMomentsPtrs MomZ(pRadAccessData->pMomZ + Offset);

		float MomX_X = *(MomX.pX), MomX_Z = *(MomX.pZ), MomZ_X = *(MomZ.pX), MomZ_Z = *(MomZ.pZ);
		float MomX_SqrtMxx = MomX.SqrtMxx, MomX_SqrtMzz = MomX.SqrtMzz, MomZ_SqrtMxx = MomZ.SqrtMxx, MomZ_SqrtMzz = MomZ.SqrtMzz;

		if(TransHndl.rep != 0)
		{
			TVector3d PointMomX(MomX_X, 0., MomX_Z);
			FromLabToLocFrame_Point(PointMomX);
			MomX_X = (float)(PointMomX.x + TransvCenPoint.x);
			MomX_Z = (float)(PointMomX.z + TransvCenPoint.y);

			TVector3d VectMomX(MomX_SqrtMxx, 0., MomX_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomX);
			MomX_SqrtMxx = (float)VectMomX.x;
			MomX_SqrtMzz = (float)VectMomX.z;

			TVector3d PointMomZ(MomZ_X, 0., MomZ_Z);
			FromLabToLocFrame_Point(PointMomZ);
			MomZ_X = (float)(PointMomZ.x + TransvCenPoint.x);
			MomZ_Z = (float)(PointMomZ.z + TransvCenPoint.y);

			TVector3d VectMomZ(MomZ_SqrtMxx, 0., MomZ_SqrtMzz);
			FromLabToLocFrame_Vector(VectMomZ);
			MomZ_SqrtMxx = (float)VectMomZ.x;
			MomZ_SqrtMzz = (float)VectMomZ.z;
		}

		float MomX_SqrtMxx_Mult = (float)(1.7*MomX_SqrtMxx), MomX_SqrtMzz_Mult = (float)(1.7*MomX_SqrtMzz), MomZ_SqrtMxx_Mult = (float)(1.7*MomZ_SqrtMxx), MomZ_SqrtMzz_Mult = (float)(1.7*MomZ_SqrtMzz);

		double CPx_mi_R = TransvCenPoint.x - R, CPx_pl_R = TransvCenPoint.x + R;
		double CPy_mi_R = TransvCenPoint.y - R, CPy_pl_R = TransvCenPoint.y + R;

		if((((MomX_X - MomX_SqrtMxx_Mult) < CPx_mi_R) || ((MomX_X + MomX_SqrtMxx_Mult) > CPx_pl_R)) ||
		   (((MomX_Z - MomX_SqrtMzz_Mult) < CPy_mi_R) || ((MomX_Z + MomX_SqrtMzz_Mult) > CPy_pl_R)) ||
		   (((MomZ_X - MomZ_SqrtMxx_Mult) < CPx_mi_R) || ((MomZ_X + MomZ_SqrtMxx_Mult) > CPx_pl_R)) ||
		   (((MomZ_Z - MomZ_SqrtMzz_Mult) < CPy_mi_R) || ((MomZ_Z + MomZ_SqrtMzz_Mult) > CPy_pl_R)))
		{
			int result;
			if(result = ComputeRadMoments(pRadAccessData)) return result;
			break;
		}
		Offset += 10;
	}

	if(FillMomRatios)
	{
		int Offset = 0;
		for(long ie=0; ie<pRadAccessData->ne; ie++)
		{
			srTMomentsPtrs NewMomX(pRadAccessData->pMomX + Offset);
			srTMomentsPtrs NewMomZ(pRadAccessData->pMomZ + Offset);

			tMomRatArray->RxxMomX *= NewMomX.SqrtMxx;
			tMomRatArray->RxpxpMomX *= NewMomX.SqrtMxpxp;
			tMomRatArray->RzzMomX *= NewMomX.SqrtMzz;
			tMomRatArray->RzpzpMomX *= NewMomX.SqrtMzpzp;
			tMomRatArray->RxxMomZ *= NewMomZ.SqrtMxx;
			tMomRatArray->RxpxpMomZ *= NewMomZ.SqrtMxpxp;
			tMomRatArray->RzzMomZ *= NewMomZ.SqrtMzz;
			tMomRatArray->RzpzpMomZ *= NewMomZ.SqrtMzpzp;

			Offset += 10;
			tMomRatArray++;
		}
	}
	return 0;
}
**/
//*************************************************************************

