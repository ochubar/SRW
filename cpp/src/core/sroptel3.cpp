/************************************************************************//**
 * File: sroptel3.cpp
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
#include "srmlttsk.h"
#include "srerror.h"

extern srTYield srYield;

//*************************************************************************

int srTGenOptElem::CheckRadStructForPropagation(srTSRWRadStructAccessData* pRad)
{
	if((pRad->nx <= 1) || (pRad->nz <= 1)) return CAN_NOT_PROPAGATE_1D_RAD;
	//if(pRad->ne > 1) return CAN_NOT_PROPAGATE_SEVERAL_PHOTON_ENERGIES; // To remove

	return 0;
}

//*************************************************************************

int srTGenOptElem::MakePostPropagationProc(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter)
{// ??? This default proc. should be applied to everything for the exception of Drifts and Rectangular Apertures (they have their own MakePostPropagationProc)
	int result;
	//OC//if(result = CheckPostResizeCenterCorrection(pRadAccessData, ResAfter)) return 0;
	//if(result = PostResizeAndTryToImproveResolInSmallSpot(pRadAccessData, ResAfter)) return 0;
	if(result = PostResizeAndTryToImproveResolInSmallSpot(pRadAccessData, ResAfter)) return result;
	return 0;
}

//*************************************************************************

int srTGenOptElem::PostResizeAndTryToImproveResolInSmallSpot(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter)
{
	int result;
	srTRadResize ImproveRes;
	char SpotShouldBeResized;
 	if(result = CheckIfSpotShouldBeResized(pRadAccessData, SpotShouldBeResized, ImproveRes)) return result;
	//ImproveRes.UseOtherSideFFT = 1;
	ImproveRes.useOtherSideFFT(1); //OC090311

	SteerPostResizeParam(pRadAccessData, ResAfter);

	char SpotWasResized = 0;
	if(SpotShouldBeResized)
	{
		int ImproveResolNow = MemoryIsSufficientForTwoResize(pRadAccessData, ImproveRes, ResAfter);
		if(ImproveResolNow)
		{
			if(result = RadResizeGen(*pRadAccessData, ImproveRes)) return result;
			SpotWasResized = 1;
		}
	}
	if(result = RadResizeGen(*pRadAccessData, ResAfter)) return result;

	if(SpotShouldBeResized && (!SpotWasResized))
	{
		int ImproveResolNow = MemoryIsSufficientForResize(pRadAccessData, ImproveRes);
		if(ImproveResolNow)
		{
			if(result = RadResizeGen(*pRadAccessData, ImproveRes)) return result;
		}
		else
		{
			if(ImproveRes.pxd > 1.) ImproveRes.pxd *= 0.8;
			if(ImproveRes.pzd > 1.) ImproveRes.pzd *= 0.8;

			while((ImproveRes.pxd > 1.1) || (ImproveRes.pzd > 1.1))
			{
				if(result = srYield.Check()) return result;

				int ImproveResolNow = MemoryIsSufficientForResize(pRadAccessData, ImproveRes);
				if(ImproveResolNow)
				{
					if(result = RadResizeGen(*pRadAccessData, ImproveRes)) return result;
					break;
				}
				if(ImproveRes.pxd > 1.) ImproveRes.pxd *= 0.8;
				if(ImproveRes.pzd > 1.) ImproveRes.pzd *= 0.8;
			}
		}
	}

	return 0;
}

//*************************************************************************

void srTGenOptElem::SteerPostResizeParam(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter)
{
	const double IntRelThresh = 0.0008; //3.e-04; // To steer
	const double NotResizeTol = 0.1; //0.05; // To steer
	int NpMin = 40; // To steer

	if((pRadAccessData->nx < NpMin) || (pRadAccessData->nz < NpMin)) return;
	if(!RangeShouldBeAdjustedAtPropag()) return; //OC ???

	double IMax = 0.;
	long ixc, izc, ix, iz;
	float *tEx = pRadAccessData->pBaseRadX, *tEz = pRadAccessData->pBaseRadZ;
	for(iz=0; iz<pRadAccessData->nz; iz++)
	{
		for(ix=0; ix<pRadAccessData->nx; ix++)
		{
			double ReEx = *(tEx++), ReEz = *(tEz++);
			double ImEx = *(tEx++), ImEz = *(tEz++);
			double I = ReEx*ReEx + ImEx*ImEx + ReEz*ReEz + ImEz*ImEz;
			if(IMax < I) { IMax = I; ixc = ix; izc = iz;}
		}
	}

	double IntThresh = IntRelThresh*IMax;
	long izMin = 0, izMax = pRadAccessData->nz - 1, ixMin = 0, ixMax = pRadAccessData->nx - 1;

	for(iz=0; iz<pRadAccessData->nz; iz++)
	{
		double IMaxInStr = MaxIntInHorString(iz, pRadAccessData);
		if(IMaxInStr >= IntThresh) { izMin = iz; break;}
	}
	for(iz=(pRadAccessData->nz - 1); iz>=0; iz--)
	{
		double IMaxInStr = MaxIntInHorString(iz, pRadAccessData);
		if(IMaxInStr >= IntThresh) { izMax = iz; break;}
	}
	for(ix=0; ix<pRadAccessData->nx; ix++)
	{
		double IMaxInStr = MaxIntInVertString(ix, izMin, izMax, pRadAccessData);
		if(IMaxInStr >= IntThresh) { ixMin = ix; break;}
	}
	for(ix=(pRadAccessData->nx - 1); ix>=0; ix--)
	{
		double IMaxInStr = MaxIntInVertString(ix, izMin, izMax, pRadAccessData);
		if(IMaxInStr >= IntThresh) { ixMax = ix; break;}
	}

	if((ixMax - ixMin - 1) < NpMin)
	{
		long HalfExtra = (NpMin - (ixMax - ixMin)) >> 1;
		ixMin -= HalfExtra;
		if(ixMin < 0) ixMin = 0;
		ixMax += HalfExtra;
		if(ixMax >= pRadAccessData->nx) ixMax = pRadAccessData->nx - 1;
	}
	if((izMax - izMin - 1) < NpMin)
	{
		long HalfExtra = (NpMin - (izMax - izMin)) >> 1;
		izMin -= HalfExtra;
		if(izMin < 0) izMin = 0;
		izMax += HalfExtra;
		if(izMax >= pRadAccessData->nz) izMax = pRadAccessData->nz - 1;
	}

	srTRadResize1D ResizeAfterX, ResizeAfterZ;
	ResizeAfterX.pm = ResAfter.pxm; ResizeAfterZ.pm = ResAfter.pzm;

	bool ModifyPmEvenIfCenPosIsNotSet = true;
	CheckRelCenPosAndSetPostResizeParamPmIfNecessary(pRadAccessData->nx, ixMin, ixMax, ResizeAfterX, ModifyPmEvenIfCenPosIsNotSet);
	CheckRelCenPosAndSetPostResizeParamPmIfNecessary(pRadAccessData->nz, izMin, izMax, ResizeAfterZ, ModifyPmEvenIfCenPosIsNotSet);

	ResAfter.pxm = ResizeAfterX.pm; ResAfter.RelCenPosX = ResizeAfterX.RelCenPos; 
	if(ResAfter.RelCenPosTol > ResizeAfterX.RelCenPosTol) ResAfter.RelCenPosTol = ResizeAfterX.RelCenPosTol;

	ResAfter.pzm = ResizeAfterZ.pm; ResAfter.RelCenPosZ = ResizeAfterZ.RelCenPos;
	if(ResAfter.RelCenPosTol > ResizeAfterZ.RelCenPosTol) ResAfter.RelCenPosTol = ResizeAfterZ.RelCenPosTol;

	if(::fabs(ResAfter.pxm - 1) < NotResizeTol) ResAfter.pxm = 1;
	if(::fabs(ResAfter.pzm - 1) < NotResizeTol) ResAfter.pzm = 1;


/**
	long ixCurCenter = (pRadAccessData->nx) >> 1;
	long ixDeltaLeft = ixCurCenter - ixMin;
	long ixDeltaRight = ixMax - ixCurCenter;
	long ixDeltaMax = (ixDeltaLeft > ixDeltaRight)? ixDeltaLeft : ixDeltaRight;
	ResAfter.pxm = double(ixDeltaMax)/double(ixCurCenter);
	if(::fabs(ResAfter.pxm - 1) < NotResizeTol) ResAfter.pxm = 1;

	long izCurCenter = (pRadAccessData->nz) >> 1;
	long izDeltaLeft = izCurCenter - izMin;
	long izDeltaRight = izMax - izCurCenter;
	long izDeltaMax = (izDeltaLeft > izDeltaRight)? izDeltaLeft : izDeltaRight;
	ResAfter.pzm = double(izDeltaMax)/double(izCurCenter);
	if(::fabs(ResAfter.pzm - 1) < NotResizeTol) ResAfter.pzm = 1;

	long ixcI = (ixMax + ixMin) >> 1, izcI = (izMax + izMin) >> 1;

	double RelCenPosX = double(ixcI)/double(pRadAccessData->nx);
	if(::fabs(RelCenPosX - 0.5) < ResAfter.RelCenPosTol) RelCenPosX = 0.5;
	double RelCenPosZ = double(izcI)/double(pRadAccessData->nz);
	if(::fabs(RelCenPosZ - 0.5) < ResAfter.RelCenPosTol) RelCenPosZ = 0.5;

	ResAfter.RelCenPosX = RelCenPosX;
	ResAfter.RelCenPosZ = RelCenPosZ;
**/

	//ResAfter.pxm = double(ixMax - ixMin)/double(pRadAccessData->nx);
	//if(::fabs(ResAfter.pxm - 1) < NotResizeTol) ResAfter.pxm = 1;
	//ResAfter.pzm = double(izMax - izMin)/double(pRadAccessData->nz);
	//if(::fabs(ResAfter.pzm - 1) < NotResizeTol) ResAfter.pzm = 1;
}

//*************************************************************************

int srTGenOptElem::MemoryIsSufficientForResize(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& Resize)
{
	double MemAvail = CheckMemoryAvailable();
	double ExtraMemForResize = ExtraMemSizeForResize(pRadAccessData->nx, pRadAccessData->nz, Resize.pxm, Resize.pxd, Resize.pzm, Resize.pzd, 0);
	return (ExtraMemForResize < MemAvail);
}

//*************************************************************************

int srTGenOptElem::MemoryIsSufficientForTwoResize(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& Resize1, srTRadResize& Resize2)
{
	double MemAvail = CheckMemoryAvailable();
	double ExtraSizeForResize1 = ExtraMemSizeForResize(pRadAccessData->nx, pRadAccessData->nz, Resize1.pxm, Resize1.pxd, Resize1.pzm, Resize1.pzd, 0);
	if(ExtraSizeForResize1 > MemAvail) return 0;
	long nxNext = (pRadAccessData->nx)*long(Resize1.pxm*Resize1.pxd);
	long nzNext = (pRadAccessData->nz)*long(Resize1.pzm*Resize1.pzd);
	double ExtraSizeForResize2 = ExtraMemSizeForResize(nxNext, nzNext, Resize2.pxm, Resize2.pxd, Resize2.pzm, Resize2.pzd, 0);
	return (ExtraSizeForResize2 < MemAvail);
}

//*************************************************************************

int srTGenOptElem::CheckIfSpotShouldBeResized(srTSRWRadStructAccessData* pRadAccessData, char& ShouldBeResized, srTRadResize& ImproveRes)
{
	int result;
	const double SpotSizeToRangeRatioToTreatSmall = 0.4; // To steer
	const double PoPerSpotImproved = 9.; // To steer
	const double RelZeroTolForIntens = 0.005; // To steer
	const double IntSizeToFringeSizeCritRatio = 4.; // To steer
	ShouldBeResized = 0;

	//double xEnd = pRadAccessData->xStart + pRadAccessData->nx*pRadAccessData->xStep;
	//double zEnd = pRadAccessData->zStart + pRadAccessData->nz*pRadAccessData->zStep;
	//char xWfrLimitsAreSmallerThanRange = (pRadAccessData->xWfrMin > pRadAccessData->xStart) || (pRadAccessData->xWfrMax < xEnd);
	//char zWfrLimitsAreSmallerThanRange = (pRadAccessData->zWfrMin > pRadAccessData->zStart) || (pRadAccessData->zWfrMax < zEnd);
	//if(xWfrLimitsAreSmallerThanRange || zWfrLimitsAreSmallerThanRange) { ShouldBeResized = 0; return 0;}
	//OC190411

	char TreatExOrEz;
	srTMomentsPtrs MomPtrsX(pRadAccessData->pMomX), MomPtrsZ(pRadAccessData->pMomZ);
	srTMomentsPtrs* pMomPtrs = 0;
	if(*(MomPtrsX.pTotPhot) > *(MomPtrsZ.pTotPhot))
	{
		pMomPtrs = &MomPtrsX; TreatExOrEz = 'x';
	}
	else
	{
		pMomPtrs = &MomPtrsZ; TreatExOrEz = 'z';
	}

	srTRadSect1D SectVsX, SectVsZ;
	SectVsX.VsXorZ = 'x'; SectVsZ.VsXorZ = 'z';

	double Xc = *(MomPtrsX.pX), Zc = *(MomPtrsX.pZ);
	long &ixc = SectVsZ.icOtherCoord;
	long &izc = SectVsX.icOtherCoord;
	ixc = IntegerOffsetCoord(pRadAccessData->xStart, pRadAccessData->xStep, Xc);
	if(ixc < 0) ixc = 0;
	if(ixc >= pRadAccessData->nx - 1) ixc = pRadAccessData->nx - 1;
	izc = IntegerOffsetCoord(pRadAccessData->zStart, pRadAccessData->zStep, Zc);
	if(izc < 0) izc = 0;
	if(izc >= pRadAccessData->nz - 1) izc = pRadAccessData->nz - 1;

	if(result = SetupSectionArraysVsXandZ(pRadAccessData, SectVsX, SectVsZ)) return result;

	long CritOffsetFromEdgeX = long(SectVsX.np*0.2), CritOffsetFromEdgeZ = long(SectVsZ.np*0.2); // To steer
	if(CritOffsetFromEdgeX == 0) CritOffsetFromEdgeX = 2;
	if(CritOffsetFromEdgeZ == 0) CritOffsetFromEdgeZ = 2;

	long iFirstX, iLastX, iFirstZ, iLastZ;
	FindIntensityBorders1D(SectVsX, TreatExOrEz, RelZeroTolForIntens, iFirstX, iLastX);
	char SmallIntAtEdgesX = ((iFirstX > CritOffsetFromEdgeX) && (iLastX < SectVsX.np - 1 - CritOffsetFromEdgeX));
	FindIntensityBorders1D(SectVsZ, TreatExOrEz, RelZeroTolForIntens, iFirstZ, iLastZ);
	char SmallIntAtEdgesZ = ((iFirstZ > CritOffsetFromEdgeZ) && (iLastZ < SectVsZ.np - 1 - CritOffsetFromEdgeZ));

	if(!(SmallIntAtEdgesX && SmallIntAtEdgesZ)) return 0;

	long iFirstSizeX, iLastSizeX;
	FindIntensityBorders1D(SectVsX, TreatExOrEz, 0.5, iFirstSizeX, iLastSizeX);
	double xSize = (iLastSizeX - iFirstSizeX)*SectVsX.ArgStep;
	long iFirstSizeZ, iLastSizeZ;
	FindIntensityBorders1D(SectVsZ, TreatExOrEz, 0.5, iFirstSizeZ, iLastSizeZ);
	double zSize = (iLastSizeZ - iFirstSizeZ)*SectVsZ.ArgStep;
	double xRange = pRadAccessData->xStep*pRadAccessData->nx;
	double zRange = pRadAccessData->zStep*pRadAccessData->nz;
	char xSpotIsSmallerThanRange = (xSize/xRange < SpotSizeToRangeRatioToTreatSmall);
	char zSpotIsSmallerThanRange = (zSize/zRange < SpotSizeToRangeRatioToTreatSmall);

	if(!(xSpotIsSmallerThanRange || zSpotIsSmallerThanRange)) return 0;

	if(xSpotIsSmallerThanRange)
	{
		double Xc, DeltaX;
		if(result = CheckWidthMax1D(SectVsX, TreatExOrEz, Xc, DeltaX)) return result;

		double xPoPerFr, xFrSize;
		if(result = AnalizeFringesAroundPoint(SectVsX, TreatExOrEz, ixc, xPoPerFr, xFrSize)) return result;
		char xFringeSizeIsLargeComparedToIntSize = (xFrSize*IntSizeToFringeSizeCritRatio > DeltaX);
		if(xFringeSizeIsLargeComparedToIntSize)
		{
			double PoPerFringeImproved = PoPerSpotImproved*xFrSize/DeltaX;
			double pd = PoPerFringeImproved/xPoPerFr;
			if(pd > 1.1) // To steer
			{
				ImproveRes.pxd = pd; ShouldBeResized = 1;
			}
		}
	}
	if(zSpotIsSmallerThanRange)
	{
		double Zc, DeltaZ;
		if(result = CheckWidthMax1D(SectVsZ, TreatExOrEz, Zc, DeltaZ)) return result;

		double zPoPerFr, zFrSize;
		if(result = AnalizeFringesAroundPoint(SectVsZ, TreatExOrEz, izc, zPoPerFr, zFrSize)) return result;
		char zFringeSizeIsLargeComparedToIntSize = (zFrSize*IntSizeToFringeSizeCritRatio > DeltaZ);
		if(zFringeSizeIsLargeComparedToIntSize)
		{
			double PoPerFringeImproved = PoPerSpotImproved*zFrSize/DeltaZ;
			double pd = PoPerFringeImproved/zPoPerFr;
			if(pd > 1.1) // To steer
			{
				ImproveRes.pzd = pd; ShouldBeResized = 1;
			}
		}
	}
	return 0;
}

//*************************************************************************

int srTGenOptElem::CheckWidthMax1D(srTRadSect1D& Sect1D, char TreatExOrEz, double& Xc, double& DeltaX)
{
	long iMax, iHalfLeft = -1, iHalfRight = -1;
	double ValMax = 0.;
	float* tE = (TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz;
	for(long i=0; i<Sect1D.np; i++)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;
		if(ValMax < Val) { ValMax = Val; iMax = i;}
		tE += 2;
	}
	float* pMaxE = (TreatExOrEz == 'x')? Sect1D.pEx + (iMax << 1) : Sect1D.pEz + (iMax << 1);

	double HalfValMax = 0.5*ValMax;
	tE = pMaxE;
	for(long im=iMax; im>=0; im--)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;
		if(Val < HalfValMax) { iHalfLeft = im; break;}
		tE -= 2;
	}
	if(iHalfLeft == -1) iHalfLeft = 0;

	tE = pMaxE;
	for(long ip=iMax; ip<Sect1D.np; ip++)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;
		if(Val < HalfValMax) { iHalfRight = ip; break;}
		tE += 2;
	}
	if(iHalfRight == -1) iHalfRight = Sect1D.np - 1;

	Xc = Sect1D.ArgStart + iMax*Sect1D.ArgStep;
	DeltaX = Sect1D.ArgStep*(iHalfRight - iHalfLeft);
	return 0;
}

//*************************************************************************

int srTGenOptElem::CheckRMS_Sizes1D(srTRadSect1D& Sect1D, char TreatExOrEz, double& Xc, double& SigmaX)
{
	double Sum0 = 0., SumX = 0., SumXX = 0.;
	float* tE = (TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz;
	double x = Sect1D.ArgStart;
	for(long i=0; i<Sect1D.np; i++)
	{
		float ReE = *tE, ImE = *(tE + 1);
		float f0 = ReE*ReE + ImE*ImE;
		float fx = (float)(f0*x);
		float fxx = (float)(fx*x);
		Sum0 += f0; SumX += fx; SumXX += fxx;

		x += Sect1D.ArgStep;
		tE += 2;
	}
	Xc = SumX/Sum0;
	double Mxx = SumXX/Sum0;
	SigmaX = sqrt(::fabs(Mxx - Xc*Xc));
	return 0;
}

//*************************************************************************

int srTGenOptElem::AnalizeFringesAroundPoint(srTRadSect1D& Sect1D, char TreatExOrEz, long ic, double& PoPerFr, double& FringeSize)
{
	int AmOfFrAvg = 3; // To steer
	int result;
	srTIntVect FringeContent;
	srTDoubleVect FringeCoor;
	if(result = CountFringes(Sect1D, FringeContent, TreatExOrEz, FringeCoor)) return result;
	long TotAmOfFringes = (long)FringeContent.size();

	double cArg = Sect1D.ArgStart + ic*Sect1D.ArgStep;
	long iFringe = 0;
	for(srTDoubleVect::iterator iter = FringeCoor.begin(); iter != FringeCoor.end(); ++iter)
	{
		if(*iter > cArg) break;
		iFringe++;
	}

	long iFrStart = iFringe - (AmOfFrAvg >> 1);
	if(iFrStart < 0) iFrStart = 0;
	long iFrEnd = iFrStart + AmOfFrAvg;
	if(iFrEnd > TotAmOfFringes - 1)
	{
		AmOfFrAvg = TotAmOfFringes - iFrStart;
	}

	double SumPoPerFr = 0., SumFrSize = 0.;
	int size_FringeCoor = (int)FringeCoor.size();
	for(int i=iFrStart; i<(iFrStart + AmOfFrAvg); i++)
	{
		SumPoPerFr += FringeContent[i];
		//SumFrSize += (i > 0)? (FringeCoor[i] - FringeCoor[i - 1]) : (FringeCoor[i + 1] - FringeCoor[i]);
		if(i > 0)
		{
			SumFrSize += FringeCoor[i] - FringeCoor[i - 1];
		}
		else
		{
			if(size_FringeCoor > 1)
			{
				SumFrSize += FringeCoor[i + 1] - FringeCoor[i];
			}
		}
	}
	PoPerFr = SumPoPerFr/double(AmOfFrAvg);
	FringeSize = SumFrSize/double(AmOfFrAvg);
	return 0;
}

//*************************************************************************

int srTGenOptElem::CheckPostResizeCenterCorrection(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter)
{
	srTMomentsPtrs MomPtrsX(pRadAccessData->pMomX), MomPtrsZ(pRadAccessData->pMomZ);

	double ReprXc, ReprZc;
	if(*(MomPtrsX.pTotPhot) > *(MomPtrsZ.pTotPhot))
	{
		ReprXc = *(MomPtrsX.pX); ReprZc = *(MomPtrsX.pZ);
	}
	else
	{
		ReprXc = *(MomPtrsZ.pX); ReprZc = *(MomPtrsZ.pZ);
	}
	
	double xRange = pRadAccessData->xStep*pRadAccessData->nx;
	double zRange = pRadAccessData->zStep*pRadAccessData->nz;

	double xcForm = pRadAccessData->xStart + 0.5*xRange;
	double zcForm = pRadAccessData->zStart + 0.5*zRange;
	double xcRes = pRadAccessData->xStart + ResAfter.RelCenPosX*xRange;
	double zcRes = pRadAccessData->zStart + ResAfter.RelCenPosZ*zRange;

	double AbsCenPosTolX = ResAfter.RelCenPosTol*xRange;
	double AbsCenPosTolZ = ResAfter.RelCenPosTol*zRange;

	if((::fabs(ReprXc - xcForm) < AbsCenPosTolX) || (::fabs(ReprXc - xcRes) > AbsCenPosTolX)) ResAfter.RelCenPosX = 0.5;
	if((::fabs(ReprZc - zcForm) < AbsCenPosTolZ) || (::fabs(ReprZc - zcRes) > AbsCenPosTolZ)) ResAfter.RelCenPosZ = 0.5;

	return 0;
}

//*************************************************************************

void srTGenOptElem::FindIntensityBorders1D(srTRadSect1D& Sect1D, char TreatExOrEz, double RelZeroTolForIntens, long& iFirst, long& iLast)
{
	long iMax = 0;
	iFirst = -1; iLast = -1;

	double ValMax = 0.;
	float* tE0 = (TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz;
	float* tE = tE0;
	for(long i=0; i<Sect1D.np; i++)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;
		if(ValMax < Val) { ValMax = Val; iMax = i;}
		tE += 2;
	}
	//float* pMaxE = (TreatExOrEz == 'x')? Sect1D.pEx + (iMax << 1) : Sect1D.pEz + (iMax << 1);

	double ValThresh = RelZeroTolForIntens*ValMax;

	tE = tE0;
	for(long im=0; im<Sect1D.np; im++)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;
		if(Val > ValThresh) { iFirst = im - 1; break;}
		tE += 2;
	}
	if(iFirst < 0) iFirst = 0;

	tE = ((TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz) + ((Sect1D.np - 1) << 1);
	for(long ip=(Sect1D.np - 1); ip>=0; ip--)
	{
		float ReE = *tE, ImE = *(tE + 1);
		double Val = ReE*ReE + ImE*ImE;

		if(Val > ValThresh) { iLast = ip + 1; break;}
		tE -= 2;
	}
	if((iLast > (Sect1D.np - 1)) || (iLast < 0)) iLast = Sect1D.np - 1;
}

//*************************************************************************

int srTGenOptElem::TuneStepToKeepInterpLimitsTheSameAtResize(srTSRWRadStructAccessData& SRWRadStructAccessData, srTSRWRadStructAccessData& NewSRWRadStructAccessData, srTRadResize& RadResize, char XorZ, long ic)
{// This can modify NewStep and NewStart !!!
	char LeftBorderNeedsCorr=0, RightBorderNeedsCorr=0;
	double NewEndInterp, OldEndInterp;
	double *pNewStart, *pOldStart, *pNewStep, *pOldStep;
	long *pNewNp;
	double pm;

	if(XorZ == 'x')
	{
		NewEndInterp = NewSRWRadStructAccessData.xStart + (NewSRWRadStructAccessData.nx - 1)*NewSRWRadStructAccessData.xStep;
		OldEndInterp = SRWRadStructAccessData.xStart + (SRWRadStructAccessData.nx - 1)*SRWRadStructAccessData.xStep;

		pNewStart = &(NewSRWRadStructAccessData.xStart);
		pNewStep = &(NewSRWRadStructAccessData.xStep);
		pOldStart = &(SRWRadStructAccessData.xStart);
		pOldStep = &(SRWRadStructAccessData.xStep);
		pNewNp = &(NewSRWRadStructAccessData.nx);
		pm = RadResize.pxm;
	}
	else if((XorZ == 'z') || (XorZ == 'y'))
	{
		NewEndInterp = NewSRWRadStructAccessData.zStart + (NewSRWRadStructAccessData.nz - 1)*NewSRWRadStructAccessData.zStep;
		OldEndInterp = SRWRadStructAccessData.zStart + (SRWRadStructAccessData.nz - 1)*SRWRadStructAccessData.zStep;

		pNewStart = &(NewSRWRadStructAccessData.zStart);
		pNewStep = &(NewSRWRadStructAccessData.zStep);
		pOldStart = &(SRWRadStructAccessData.zStart);
		pOldStep = &(SRWRadStructAccessData.zStep);
		pNewNp = &(NewSRWRadStructAccessData.nz);
		pm = RadResize.pzm;
	}
	else //if(XorZ == 'e') //OC111111
	{
		NewEndInterp = NewSRWRadStructAccessData.eStart + (NewSRWRadStructAccessData.ne - 1)*NewSRWRadStructAccessData.eStep;
		OldEndInterp = SRWRadStructAccessData.eStart + (SRWRadStructAccessData.ne - 1)*SRWRadStructAccessData.eStep;
	
		pNewStart = &(NewSRWRadStructAccessData.eStart);
		pNewStep = &(NewSRWRadStructAccessData.eStep);
		pOldStart = &(SRWRadStructAccessData.eStart);
		pOldStep = &(SRWRadStructAccessData.eStep);
		pNewNp = &(NewSRWRadStructAccessData.ne);
		pm = RadResize.pem;
	}

	LeftBorderNeedsCorr = (*pNewStart <= *pOldStart);
	RightBorderNeedsCorr = (NewEndInterp >= OldEndInterp);

	if(!(LeftBorderNeedsCorr || RightBorderNeedsCorr)) return 0;

	double cArg = *pOldStart + ic*(*pOldStep);
	long HalfNewNp = *pNewNp >> 1;
	if(LeftBorderNeedsCorr && (!RightBorderNeedsCorr))
	{
		double DistLeft = cArg - *pNewStart;
		long AmOfNewStepsLeft = long(DistLeft/(*pNewStep) + 1.E-10);

		*pNewStep = DistLeft/double(AmOfNewStepsLeft);
		*pNewStart = cArg - HalfNewNp*(*pNewStep);
	}
	else if((!LeftBorderNeedsCorr) && RightBorderNeedsCorr)
	{
		double DistRight = OldEndInterp - cArg;
		long AmOfNewStepsRight = long(DistRight/(*pNewStep) + 1.E-10);

		*pNewStep = DistRight/double(AmOfNewStepsRight);
		*pNewStart = cArg - HalfNewNp*(*pNewStep);
	}
	else if(LeftBorderNeedsCorr && RightBorderNeedsCorr)
	{
		double OldRange = OldEndInterp - *pOldStart;
		if(pm != 1.)
		{
			long AmOfNewStepsInsideOldRange = long(OldRange/(*pNewStep) + 1.E-10);
			long AmOfNewStepsLeftFromOldRange = long((*pOldStart - *pNewStart)/(*pNewStep) + 1.E-10);

			*pNewStep = OldRange/double(AmOfNewStepsInsideOldRange);
			*pNewStart = *pOldStart - AmOfNewStepsLeftFromOldRange*(*pNewStep);
		}
		else
		{
			*pNewStep = OldRange/double(*pNewNp - 1);
			*pNewStart = *pOldStart;
		}
	}
	return 0;
}

//*************************************************************************

char srTGenOptElem::UnderSamplingModeCanBeSuggested(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario)
{
	const double MinPd = 2.; // To steer
	const int MaxNpResized = 800; // To steer

	if(!SuitableConditionsForUnderSamplingMode(pRadAccessData, PropagScenario)) return 0;

	srTPropagScenario1D *pScenX = PropagScenario, *pScenZ = PropagScenario + 1;
	double &pxd = pScenX->ResizeBefore.pd; //, &pxm = pScenX->ResizeBefore.pm;
	double &pzd = pScenZ->ResizeBefore.pd; //, &pzm = pScenZ->ResizeBefore.pm;

	//double pdMin = (pxd < pzd)? pxd : pzd;
	double pdMax = (pxd > pzd)? pxd : pzd;

	if(pdMax < MinPd) return 0;
	if((pRadAccessData->nx*pxd < MaxNpResized) && (pRadAccessData->nz*pzd < MaxNpResized)) return 0;

	return 1;
}


//*************************************************************************

char srTGenOptElem::SuitableConditionsForUnderSamplingMode(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario)
{
	//const double MinPd = 1.1; //1.5; // To steer
	//const double MinPm = 1.1; //1.5; // To steer

	const double MinPd = 1.5; // To steer
	const double MinDifPm = 0.5; // To steer

	srTPropagScenario1D *pScenX = PropagScenario, *pScenZ = PropagScenario + 1;
	double &pxd = pScenX->ResizeBefore.pd, &pxm = pScenX->ResizeBefore.pm;
	double &pzd = pScenZ->ResizeBefore.pd, &pzm = pScenZ->ResizeBefore.pm;

	char ConditionsForUnderSamplingMode = 0;

	char PxdIsLarge = ((pxd > MinPd) && (::fabs(pxm - 1.) < MinDifPm));
	char PzdIsLarge = ((pzd > MinPd) && (::fabs(pzm - 1.) < MinDifPm));

	//char PxdIsLarge = (pxd > MinPd); //OC291206
	//char PzdIsLarge = (pzd > MinPd);
	//char PxmIsLarge = (pxm > MinPm);
	//char PzmIsLarge = (pzm > MinPm);

	char RxIsGood = (::fabs(pRadAccessData->RobsX) > 3.*pRadAccessData->RobsXAbsErr); // To steer
	char RzIsGood = (::fabs(pRadAccessData->RobsZ) > 3.*pRadAccessData->RobsZAbsErr); // To steer

	if((PxdIsLarge || PzdIsLarge) && RxIsGood && RzIsGood)
	//if((PxdIsLarge || PzdIsLarge || PxmIsLarge || PzmIsLarge) && RxIsGood && RzIsGood) //OC291206
	{
		ConditionsForUnderSamplingMode = 1;
	}
	// Make more tests (?)
	// Make checking for particular optical component

	return ConditionsForUnderSamplingMode;
}

//*************************************************************************

void srTGenOptElem::ShowCurrentOverSamplingFactors(srTPropagScenario1D* PropagScenario, double& Fx, double& Fz)
{
	const double NominalPointsPerFringe = 1.2;

	srTFringeInfo &FringeInfoX = PropagScenario[0].CurFringeInfo, &FringeInfoZ = PropagScenario[1].CurFringeInfo;
	double PoPerFringeX = (FringeInfoX.LeftPointsPerFringe < FringeInfoX.RightPointsPerFringe)? FringeInfoX.LeftPointsPerFringe : FringeInfoX.RightPointsPerFringe;
	double PoPerFringeZ = (FringeInfoZ.LeftPointsPerFringe < FringeInfoZ.RightPointsPerFringe)? FringeInfoZ.LeftPointsPerFringe : FringeInfoZ.RightPointsPerFringe;

	Fx = PoPerFringeX/NominalPointsPerFringe;
	Fz = PoPerFringeZ/NominalPointsPerFringe;
}

//*************************************************************************

long srTGenOptElem::EstimateMinNpForQuadTerm(double en, double R, double xStart, double xEnd)
{
	//double WavelengthIn_m = 1.239854E-06/en;
	double WavelengthIn_m = 1.239842E-06/en;

	if(R == 0.) return 1000000000;
	R = ::fabs(R);
	if(::fabs(xStart) < WavelengthIn_m) xStart = WavelengthIn_m;
	if(::fabs(xEnd) < WavelengthIn_m) xEnd = WavelengthIn_m;

	double HalfLambR = 0.5*WavelengthIn_m*R;
	double dxStart = ::fabs(HalfLambR/xStart);
	double dxEnd = ::fabs(HalfLambR/xEnd);
	double dx = ((dxStart < dxEnd)? dxStart : dxEnd)/1.2; // To steer

	long Np = int(::fabs(xEnd - xStart)/dx) + 1;
	if(((Np >> 1) << 1) != Np) Np++;
	return Np;
}

//*************************************************************************

int srTGenOptElem::EstimateNominalNpForUnderSampling(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* Sect1D, double& NomNxForUnderSampl, double& NomNzForUnderSampl)
{
	const double MinNominalNpForUnderSampling = 150; // To steer: this corresponds to OverSampling = 1
	NomNxForUnderSampl = MinNominalNpForUnderSampling;
	NomNzForUnderSampl = MinNominalNpForUnderSampling;

	double AbsLostRx1 = ::fabs(-(pRadAccessData->RobsX)*(pRadAccessData->RobsX + pRadAccessData->RobsXAbsErr)/(pRadAccessData->RobsXAbsErr));
	double AbsLostRx2 = ::fabs((pRadAccessData->RobsX)*(pRadAccessData->RobsX - pRadAccessData->RobsXAbsErr)/(pRadAccessData->RobsXAbsErr));
	double MinAbsLostRx = (AbsLostRx1 < AbsLostRx2)? AbsLostRx1 : AbsLostRx2;
	double xStartRel = pRadAccessData->xStart - pRadAccessData->xc;
	double xEndRel = xStartRel + (pRadAccessData->nx - 1)*pRadAccessData->xStep;
	long NxLostQuadTerm = EstimateMinNpForQuadTerm(pRadAccessData->eStart, MinAbsLostRx, xStartRel, xEndRel);
	if(NomNxForUnderSampl < double(NxLostQuadTerm)) NomNxForUnderSampl = NxLostQuadTerm;

	double AbsLostRz1 = ::fabs(-(pRadAccessData->RobsZ)*(pRadAccessData->RobsZ + pRadAccessData->RobsZAbsErr)/(pRadAccessData->RobsZAbsErr));
	double AbsLostRz2 = ::fabs((pRadAccessData->RobsZ)*(pRadAccessData->RobsZ - pRadAccessData->RobsZAbsErr)/(pRadAccessData->RobsZAbsErr));
	double MinAbsLostRz = (AbsLostRz1 < AbsLostRz2)? AbsLostRz1 : AbsLostRz2;
	double zStartRel = pRadAccessData->zStart - pRadAccessData->zc;
	double zEndRel = zStartRel + (pRadAccessData->nz - 1)*pRadAccessData->zStep;
	long NzLostQuadTerm = EstimateMinNpForQuadTerm(pRadAccessData->eStart, MinAbsLostRz, zStartRel, zEndRel);
	if(NomNzForUnderSampl < double(NzLostQuadTerm)) NomNzForUnderSampl = NzLostQuadTerm;

	double MinNxToResolveOptElm, MinNzToResolveOptElm;
	EstimateMinNpToResolveOptElem(pRadAccessData, MinNxToResolveOptElm, MinNzToResolveOptElm);
	if(NomNxForUnderSampl < MinNxToResolveOptElm) NomNxForUnderSampl = MinNxToResolveOptElm;
	if(NomNzForUnderSampl < MinNzToResolveOptElm) NomNzForUnderSampl = MinNzToResolveOptElm;
	return 0;
}

//*************************************************************************

int srTGenOptElem::TryToSetUnderSamplingMode(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* Sect1D, srTPropagScenario1D* PropagScenario, char& UnderSamplingModeWasSet)
{
	//const int MinNpAllowed = 50; // To steer
	//const double ReduceFactor = 0.8; // To steer
	const double MinDifPd_ToAllowResize = 0.15; // To steer
	const double MinDifPd_ToAllowUnderSampling = 0.25; // To steer

	//const double InvReduceFactor = 1./ReduceFactor;

	if(!SuitableConditionsForUnderSamplingMode(pRadAccessData, PropagScenario))
	{
		UnderSamplingModeWasSet = 0; return 0;
	}

	double NomNxForUnderSampling, NomNzForUnderSampling;
	EstimateNominalNpForUnderSampling(pRadAccessData, Sect1D, NomNxForUnderSampling, NomNzForUnderSampling);


//ddddddddddddddddddddddddddddddddddd


	//double Fx, Fz;
	//ShowCurrentOverSamplingFactors(PropagScenario, Fx, Fz);
	//double dNxUnderSampl = NomNxForUnderSampling*Fx;
	//double dNzUnderSampl = NomNzForUnderSampling*Fz;

	srTRadResize1D &ResBeforeX = PropagScenario->ResizeBefore, &ResBeforeZ = (PropagScenario + 1)->ResizeBefore;
	double &pxd = ResBeforeX.pd; //, &pxm = ResBeforeX.pm;
	double &pzd = ResBeforeZ.pd; //, &pzm = ResBeforeZ.pm;

	//double dNxRequired = pRadAccessData->nx*pxd;
	//double dNzRequired = pRadAccessData->nz*pzd;

	//double CurUnderSamplingX = dNxRequired/dNxUnderSampl;
	//double CurUnderSamplingZ = dNzRequired/dNzUnderSampl;

	double CurUnderSamplingX = 1;
	double CurUnderSamplingZ = 1;

	if(((pxd - 1.) > MinDifPd_ToAllowUnderSampling))
	{
		CurUnderSamplingX = pxd; pxd = 1.;
		
		double TestRatio = double(NomNxForUnderSampling/(pRadAccessData->nx));
		if((TestRatio - 1.) > MinDifPd_ToAllowResize)
		{
			pxd = TestRatio; CurUnderSamplingX /= TestRatio;
			if(CurUnderSamplingX < 1.) 
			{
				pxd *= CurUnderSamplingX; CurUnderSamplingX = 1;
			}
		}
	}
	if(((pzd - 1.) > MinDifPd_ToAllowUnderSampling))
	{
		CurUnderSamplingZ = pzd; pzd = 1.;
		
		double TestRatio = double(NomNzForUnderSampling/(pRadAccessData->nz));
		if((TestRatio - 1.) > MinDifPd_ToAllowResize)
		{
			pzd = TestRatio; CurUnderSamplingZ /= TestRatio;
			if(CurUnderSamplingZ < 1.) 
			{
				pzd *= CurUnderSamplingZ; CurUnderSamplingZ = 1;
			}
		}
	}


/**
	double pxdNew = dNxUnderSampl/double(pRadAccessData->nx);
	double pzdNew = dNzUnderSampl/double(pRadAccessData->nz);

	if(dNxUnderSampl < dNxRequired) pxd = pxdNew;
	else CurUnderSamplingX = 1.;

	if(dNzUnderSampl < dNzRequired) pzd = pzdNew;
	else CurUnderSamplingZ = 1.;

	double ExtraSizeBeforeProp, PrevExtraSizeBeforeProp, DummyExtraSizeAfterProp;
	char FitMemOK = MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, DummyExtraSizeAfterProp);
	char PrecWasReduced = 0;

	while(!FitMemOK)
	{
		pxd *= ReduceFactor; pzd *= ReduceFactor; 
		CurUnderSamplingX *= InvReduceFactor; CurUnderSamplingZ *= InvReduceFactor;
		PrecWasReduced = 1;

		PrevExtraSizeBeforeProp = ExtraSizeBeforeProp;
		FitMemOK = MemoryIsSufficientForProp(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, DummyExtraSizeAfterProp);
		if(PrevExtraSizeBeforeProp >= ExtraSizeBeforeProp) return PROP_CAN_NOT_BE_DONE_DUE_TO_MEMORY_LIMIT;
	}
	if(PrecWasReduced)
	{
		srTSend Send; Send.AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
	}

	if((::fabs(pxd - 1.) < MinDifPd_ToAllowResize) && (CurUnderSamplingX > 1.) && (pRadAccessData->nx > MinNpAllowed))
	{
		CurUnderSamplingX *= pxd; pxd = 1.;
	}
	if((::fabs(pzd - 1.) < MinDifPd_ToAllowResize) && (CurUnderSamplingZ > 1.) && (pRadAccessData->nz > MinNpAllowed))
	{
		CurUnderSamplingZ *= pzd; pzd = 1.;
	}
**/

	if((pxd < 1.) && (CurUnderSamplingX > 1.))// OC??
	{
		CurUnderSamplingX *= pxd; pxd = 1.;
	}
	if((pzd < 1.) && (CurUnderSamplingZ > 1.))
	{
		CurUnderSamplingZ *= pzd; pzd = 1.;
	}

	if((::fabs(CurUnderSamplingX - 1.) < MinDifPd_ToAllowUnderSampling))
	{
		pxd *= CurUnderSamplingX; CurUnderSamplingX = 1; 
	}
	if((::fabs(CurUnderSamplingZ - 1.) < MinDifPd_ToAllowUnderSampling))
	{
		pzd *= CurUnderSamplingZ; CurUnderSamplingZ = 1; 
	}

	pRadAccessData->UnderSamplingX *= CurUnderSamplingX;
	pRadAccessData->UnderSamplingZ *= CurUnderSamplingZ;

	UnderSamplingModeWasSet = (pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ())? 1 : 0;
	return 0;
}

//*************************************************************************

int srTGenOptElem::TryToRemoveUndersamplingByResizing(srTSRWRadStructAccessData* pRadAccessData)
{
	int result;
	const double pxdReduceFact = 0.8; //To steer

	srTRadResize Resize;
	double &pxd = Resize.pxd, &pxm = Resize.pxm;
	double &pzd = Resize.pzd, &pzm = Resize.pzm;

	pxd = pRadAccessData->UnderSamplingX; pxm = 1.;
	pzd = pRadAccessData->UnderSamplingZ; pzm = 1.;

	if((pxd == 1.) && (pzd == 1.)) return 0;

	double MemForResize = 1.E+23, PrevMemForResize;
	char ResolutionWasReduced = 0;
	while((pxd > 1.) && (pzd > 1.))
	{
		PrevMemForResize = MemForResize;
		long nxCurRad = pRadAccessData->nx, nzCurRad = pRadAccessData->nz;
		MemForResize = ExtraMemSizeForResize(nxCurRad, nzCurRad, pxm, pxd, pzm, pzd, 0);

		if(MemForResize >= PrevMemForResize) break;

		double MemAvail = CheckMemoryAvailable();
		if(MemAvail >= MemForResize)
		{
			if(result = RadResizeGen(*pRadAccessData, Resize)) return result;
			pRadAccessData->UnderSamplingX = 1.;
			pRadAccessData->UnderSamplingZ = 1.;
			return 0;
		}

		if(pxd > 1.) 
		{
			pxd *= pxdReduceFact; ResolutionWasReduced = 1;
		}
		if(pzd > 1.)
		{
			pzd *= pxdReduceFact; ResolutionWasReduced = 1;
		}
	}

	return PROP_CAN_NOT_BE_DONE_DUE_TO_MEMORY_LIMIT;
}

//*************************************************************************

int srTGenOptElem::RemoveUndersamplingByResizingOrStop(srTSRWRadStructAccessData* pRadAccessData)
{
	int result = 0;
	
	double GoodNx = pRadAccessData->nx*pRadAccessData->UnderSamplingX;
	double GoodNz = pRadAccessData->nz*pRadAccessData->UnderSamplingZ;

	if(((int)(GoodNx + 1.e-12) == pRadAccessData->nx) && ((int)(GoodNz + 1.e-12) == pRadAccessData->nz)) return 0;
	
	if(result = TryToRemoveUndersamplingByResizing(pRadAccessData)) return result;
	
	if((GoodNx*0.7 > double(pRadAccessData->nx)) || (GoodNz*0.7 > double(pRadAccessData->nz)))
	{// To steer
		CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
	}
	return result;
}

//*************************************************************************

int srTGenOptElem::ReduceBiggerResizeParamUntilFitMemory(srTSRWRadStructAccessData& RadAccessData, srTRadResize& Resize, double& UnderSamplingX, double& UnderSamplingZ)
{
	int result = 0;
	const double ReductionCoef = 0.8;
	if(MemoryIsSufficientForResize(&RadAccessData, Resize)) return result;

	double pxd0 = Resize.pxd, pxm0 = Resize.pxm, pzd0 = Resize.pzd, pzm0 = Resize.pzm;
	double* PtrsToCoef[] = {&(Resize.pxd), &(Resize.pxm), &(Resize.pzd), &(Resize.pzm)};
	char ReductionNeeded = 1;

	while(ReductionNeeded)
	{
		ReductionNeeded = 0;
		for(int j=0; j<4; j++) if(*(PtrsToCoef[j]) > 1.) { *(PtrsToCoef[j]) *= ReductionCoef; ReductionNeeded = 1;}
		if(MemoryIsSufficientForResize(&RadAccessData, Resize)) break;
	}
	UnderSamplingX = (pxd0/(*(PtrsToCoef[0])))*(pxm0/(*(PtrsToCoef[1])));
	UnderSamplingZ = (pzd0/(*(PtrsToCoef[2])))*(pzm0/(*(PtrsToCoef[3])));
	return result;
}

//*************************************************************************

int srTGenOptElem::ChangeWfrRepres(srTSRWRadStructAccessData* pRad, int MethNo)
{
	int res = 0;
	double GoodNx = pRad->nx*pRad->UnderSamplingX;
	double GoodNz = pRad->nz*pRad->UnderSamplingZ;

	if(res = TryToRemoveUndersamplingByResizing(pRad)) return res;
	if(pRad->ThereIsUnderSamplingX() || pRad->ThereIsUnderSamplingZ()) return PROP_CAN_NOT_BE_DONE_DUE_TO_MEMORY_LIMIT;

	if((GoodNx*0.7 > double(pRad->nx)) || (GoodNz*0.7 > double(pRad->nz)))
	{// To steer
		CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
	}

	if(MethNo == 0) return ChangeWfrRepresMeth_0(pRad);
	else return ChangeWfrRepresMeth_1(pRad);

	return 0;
} 

//*************************************************************************

int srTGenOptElem::ChangeWfrRepresMeth_0(srTSRWRadStructAccessData* pRad)
{
	int res = 0;
	if(pRad == 0) return 0;

	int DirFromCoordToAng = 0;
	if(pRad->Pres != 1) DirFromCoordToAng = 1;

	if(res = ComputeRadMoments(pRad)) return res;
	if(res = SetRadRepres(pRad, DirFromCoordToAng)) return res;
	pRad->SetNonZeroWavefrontLimitsToFullRange();
	if(res = ComputeRadMoments(pRad)) return res;

	return 0;
}

//*************************************************************************

int srTGenOptElem::ChangeWfrRepresMeth_1(srTSRWRadStructAccessData* pRad)
{
//dddddddddddddddddddddddddddddddd
	return 0;
}

//*************************************************************************

