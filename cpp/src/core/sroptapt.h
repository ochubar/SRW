/************************************************************************//**
 * File: sroptapt.h
 * Description: Optical element: Aperture (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTAPT_H
#define __SROPTAPT_H

#include "sroptshp.h"

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include <stdio.h>
//#include "srwlib.h"

//*************************************************************************

class srTAperture : public srTShapedOptElem {

	double xStartNonZeroOld, zStartNonZeroOld;
	double xStartZeroAgainOld, zStartZeroAgainOld;
	double xStartNonZeroNew, zStartNonZeroNew;
	double xStartZeroAgainNew, zStartZeroAgainNew;

public:
	srTAperture () {}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		char &MethNo = ParPrecWfrPropag.MethNo;

		if(MethNo == 0) return PropagateRadiationMeth_0(pRadAccessData);
		else if(MethNo == 1) return PropagateRadiationMeth_1(pRadAccessData);
		//else if(MethNo == 2) return PropagateRadiationMeth_2(pRadAccessData, ResBeforeAndAfterVect);
		else if(MethNo == 2) return PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);

		return 0;
	}
	//int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;

		SetNewNonZeroWfrLimits(pRadAccessData);
		return 0;
	}
	int PropagateRadiationMeth_1(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE_TracingZeroField(pRadAccessData)) return result;
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = TuneRadAfterPropMeth_1(pRadAccessData)) return result;

		SetNewNonZeroWfrLimits(pRadAccessData);
		return 0;
	}

	int TraverseRadZXE_TracingZeroField(srTSRWRadStructAccessData*);
	int TuneRadAfterPropMeth_1(srTSRWRadStructAccessData*);

	virtual void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData*) {}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;

		SetNewNonZeroWfrLimits(pRadAccessData);
		return 0;
	}
	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		int result;
		if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
		if(result = TraverseRad1D(pSect1D)) return result;
		return 0;
	}

	//int RangeShouldBeAdjustedAtPropag() { return 0;}
	int RangeShouldBeAdjustedAtPropag() { return 1;}
	int ResolutionShouldBeAdjustedAtPropag() { return 0;}

	int PropagateRadMoments(srTSRWRadStructAccessData*, srTMomentsRatios*);
	//virtual int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) { return 1;}
	virtual int CheckIfMomentsShouldBeRecomputed(double MomX_X, double MomX_Z, double MomZ_X, double MomZ_Z, double MomX_SqrtMxx_Mult, double MomX_SqrtMzz_Mult, double MomZ_SqrtMxx_Mult, double MomZ_SqrtMzz_Mult) { return 1;} //OC130311
};

//*************************************************************************

class srTRectAperture : public srTAperture {
protected:
	double HalfDx, HalfDz;
public:
	double Dx, Dz;

	srTRectAperture(srTStringVect* pElemInfo)
	{
		Dx = atof((*pElemInfo)[1]); // input in m
		Dz = atof((*pElemInfo)[2]); // input in m
		HalfDx = 0.5*Dx; HalfDz = 0.5*Dz;

		if(pElemInfo->size() > 3)
		{
			TransvCenPoint.x = atof((*pElemInfo)[3]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[4]); // input in m
		}

		WfrTransmLimits.Setup(TransvCenPoint.x - 0.5*Dx, Dx, TransvCenPoint.y - 0.5*Dz, Dz);
	}
	srTRectAperture(double InDx, double InDz, double InCx, double InCz)
	{
		Dx = InDx; HalfDx = 0.5*Dx;
		Dz = InDz; HalfDz = 0.5*Dz;
		TransvCenPoint.x = InCx;
		TransvCenPoint.y = InCz;

		WfrTransmLimits.Setup(TransvCenPoint.x - 0.5*Dx, Dx, TransvCenPoint.y - 0.5*Dz, Dz);
	}

	srTRectAperture() {}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		const double SmallOffset = 1.E-10;
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;

		if(TransHndl.rep != 0)
		{
			TVector3d Point(EXZ.x, 0., EXZ.z);
			FromLabToLocFrame_Point(Point);
			xRel = Point.x; zRel = Point.z;
		}

		//if((xRel < -HalfDx) || (xRel >= HalfDx) || (zRel < -HalfDz) || (zRel >= HalfDz))
		double EffHalfDx = HalfDx + SmallOffset, EffHalfDz = HalfDz + SmallOffset;
		if((xRel < -EffHalfDx) || (xRel > EffHalfDx) || (zRel < -EffHalfDz) || (zRel > EffHalfDz))
		{
			*(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;
			*(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;
		}
	}
	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		RadPointModifier(EXZ, EPtrs); // Check this
	}

	void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData* pRadAccessData) 
	{
		double xWfrMinNew = TransvCenPoint.x - HalfDx;
		if(pRadAccessData->xWfrMin < xWfrMinNew) pRadAccessData->xWfrMin = xWfrMinNew;
		
		double xWfrMaxNew = TransvCenPoint.x + HalfDx;
		if(pRadAccessData->xWfrMax > xWfrMaxNew) pRadAccessData->xWfrMax = xWfrMaxNew;

		double zWfrMinNew = TransvCenPoint.y - HalfDz;
		if(pRadAccessData->zWfrMin < zWfrMinNew) pRadAccessData->zWfrMin = zWfrMinNew;

		double zWfrMaxNew = TransvCenPoint.y + HalfDz;
		if(pRadAccessData->zWfrMax > zWfrMaxNew) pRadAccessData->zWfrMax = zWfrMaxNew;
	}

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
	{
		const double MinPoInAperture = 6; // To steer
		double xWfrSize = pRadAccessData->xWfrMax - pRadAccessData->xWfrMin;
		MinNx = (xWfrSize/Dx)*MinPoInAperture;
		double zWfrSize = pRadAccessData->zWfrMax - pRadAccessData->zWfrMin;
		MinNz = (zWfrSize/Dz)*MinPoInAperture;
		return 0;
	}

	//int PropagateRadMoments(srTSRWRadStructAccessData*, srTMomentsRatios*);
	//int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) 
	int CheckIfMomentsShouldBeRecomputed(double MomX_X, double MomX_Z, double MomZ_X, double MomZ_Z, double MomX_SqrtMxx_Mult, double MomX_SqrtMzz_Mult, double MomZ_SqrtMxx_Mult, double MomZ_SqrtMzz_Mult) //OC130311
	{
        if((((MomX_X - MomX_SqrtMxx_Mult) < (TransvCenPoint.x - HalfDx)) || ((MomX_X + MomX_SqrtMxx_Mult) > (TransvCenPoint.x + HalfDx))) ||
		   (((MomX_Z - MomX_SqrtMzz_Mult) < (TransvCenPoint.y - HalfDz)) || ((MomX_Z + MomX_SqrtMzz_Mult) > (TransvCenPoint.y + HalfDz))) ||
		   (((MomZ_X - MomZ_SqrtMxx_Mult) < (TransvCenPoint.x - HalfDx)) || ((MomZ_X + MomZ_SqrtMxx_Mult) > (TransvCenPoint.x + HalfDx))) ||
		   (((MomZ_Z - MomZ_SqrtMzz_Mult) < (TransvCenPoint.y - HalfDz)) || ((MomZ_Z + MomZ_SqrtMzz_Mult) > (TransvCenPoint.y + HalfDz)))) return 1;
		else return 0;
	}
};

//*************************************************************************

class srTRectObstacle : public srTRectAperture {
public:

	srTRectObstacle(srTStringVect* pElemInfo)
	{
		Dx = atof((*pElemInfo)[1]); // input in m
		Dz = atof((*pElemInfo)[2]); // input in m
		HalfDx = 0.5*Dx; HalfDz = 0.5*Dz;

		if(pElemInfo->size() > 3)
		{
			TransvCenPoint.x = atof((*pElemInfo)[3]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[4]); // input in m
		}
		//WfrTransmLimits.Setup(TransvCenPoint.x - 0.5*Dx, Dx, TransvCenPoint.y - 0.5*Dz, Dz);
	}
	srTRectObstacle(double InDx, double InDz, double InCx, double InCz)
	{
		Dx = InDx; HalfDx = 0.5*Dx;
		Dz = InDz; HalfDz = 0.5*Dz;
		TransvCenPoint.x = InCx;
		TransvCenPoint.y = InCz;
		//WfrTransmLimits.Setup(TransvCenPoint.x - 0.5*Dx, Dx, TransvCenPoint.y - 0.5*Dz, Dz);
	}
	srTRectObstacle() {}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		const double SmallOffset = 1.E-10;
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;

		if(TransHndl.rep != 0)
		{
			TVector3d Point(EXZ.x, 0., EXZ.z);
			FromLabToLocFrame_Point(Point);
			xRel = Point.x; zRel = Point.z;
		}

		double EffHalfDx = HalfDx + SmallOffset, EffHalfDz = HalfDz + SmallOffset;
		if((xRel >= -EffHalfDx) && (xRel <= EffHalfDx) && (zRel >= -EffHalfDz) && (zRel <= EffHalfDz))
		{
			*(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;
			*(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;
		}
	}

	void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData* pRadAccessData) 
	{
		double xWfrMinEdge = TransvCenPoint.x - HalfDx;
		double xWfrMaxEdge = TransvCenPoint.x + HalfDx;
		double zWfrMinEdge = TransvCenPoint.y - HalfDz;
		double zWfrMaxEdge = TransvCenPoint.y + HalfDz;

		//if((pRadAccessData->xWfrMin < xWfrMaxEdge) && (xWfrMaxEdge < pRadAccessData->xWfrMax)) pRadAccessData->xWfrMin = xWfrMaxEdge;
		//if((pRadAccessData->xWfrMin < xWfrMinEdge) && (xWfrMinEdge < pRadAccessData->xWfrMax)) pRadAccessData->xWfrMax = xWfrMinEdge;
		//if((pRadAccessData->zWfrMin < zWfrMaxEdge) && (zWfrMaxEdge < pRadAccessData->zWfrMax)) pRadAccessData->zWfrMin = zWfrMaxEdge;
		//if((pRadAccessData->zWfrMin < zWfrMinEdge) && (zWfrMinEdge < pRadAccessData->zWfrMax)) pRadAccessData->zWfrMax = zWfrMinEdge;

		//OC fix 050107
		if((xWfrMinEdge < pRadAccessData->xWfrMin) && (pRadAccessData->xWfrMin < xWfrMaxEdge) &&
		   (zWfrMinEdge < pRadAccessData->zWfrMin) && (pRadAccessData->zWfrMin < zWfrMaxEdge) &&
		   (zWfrMinEdge < pRadAccessData->zWfrMax) && (pRadAccessData->zWfrMax < zWfrMaxEdge)) pRadAccessData->xWfrMin = xWfrMaxEdge;
		if((xWfrMinEdge < pRadAccessData->xWfrMax) && (pRadAccessData->xWfrMax < xWfrMaxEdge) &&
		   (zWfrMinEdge < pRadAccessData->zWfrMin) && (pRadAccessData->zWfrMin < zWfrMaxEdge) &&
		   (zWfrMinEdge < pRadAccessData->zWfrMax) && (pRadAccessData->zWfrMax < zWfrMaxEdge)) pRadAccessData->xWfrMax = xWfrMinEdge;

		if((zWfrMinEdge < pRadAccessData->zWfrMin) && (pRadAccessData->zWfrMin < zWfrMaxEdge) &&
		   (xWfrMinEdge < pRadAccessData->xWfrMin) && (pRadAccessData->xWfrMin < xWfrMaxEdge) &&
		   (xWfrMinEdge < pRadAccessData->xWfrMax) && (pRadAccessData->xWfrMax < xWfrMaxEdge)) pRadAccessData->zWfrMin = zWfrMaxEdge;
		if((zWfrMinEdge < pRadAccessData->zWfrMax) && (pRadAccessData->zWfrMax < zWfrMaxEdge) &&
		   (xWfrMinEdge < pRadAccessData->xWfrMin) && (pRadAccessData->xWfrMin < xWfrMaxEdge) &&
		   (xWfrMinEdge < pRadAccessData->xWfrMax) && (pRadAccessData->xWfrMax < xWfrMaxEdge)) pRadAccessData->zWfrMax = zWfrMinEdge;
	}

	//int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) { return 1;}
	int CheckIfMomentsShouldBeRecomputed(double MomX_X, double MomX_Z, double MomZ_X, double MomZ_Z, double MomX_SqrtMxx_Mult, double MomX_SqrtMzz_Mult, double MomZ_SqrtMxx_Mult, double MomZ_SqrtMzz_Mult) { return 1;} //OC130311
};

//*************************************************************************

class srTCircAperture : public srTAperture {
protected:
	double R, Re2;
public:
	double D; // diameter

	srTCircAperture(srTStringVect* pElemInfo)
	{
		D = atof((*pElemInfo)[1]); // input in m
		R = 0.5*D; Re2 = R*R;

		if(pElemInfo->size() > 2)
		{
			TransvCenPoint.x = atof((*pElemInfo)[2]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[3]); // input in m
		}

		WfrTransmLimits.Setup(TransvCenPoint.x - R, D, TransvCenPoint.y - R, D);
	}
	srTCircAperture(double InD, double InCx, double InCz)
	{
		D = InD;
		R = 0.5*D; Re2 = R*R;
		TransvCenPoint.x = InCx;
		TransvCenPoint.y = InCz;

		WfrTransmLimits.Setup(TransvCenPoint.x - R, D, TransvCenPoint.y - R, D);
	}

	srTCircAperture() {}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;

		if(TransHndl.rep != 0)
		{
			TVector3d Point(EXZ.x, 0., EXZ.z);
			FromLabToLocFrame_Point(Point);
			xRel = Point.x; zRel = Point.z;
		}
		if(xRel*xRel + zRel*zRel > Re2)
		{
			*(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;
			*(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;
		}
	}
	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		RadPointModifier(EXZ, EPtrs); // Check this
	}

	void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData* pRadAccessData) 
	{
		double WfrMinNew = TransvCenPoint.x - R;
		if(pRadAccessData->xWfrMin < WfrMinNew) pRadAccessData->xWfrMin = WfrMinNew;
		
		double WfrMaxNew = TransvCenPoint.x + R;
		if(pRadAccessData->xWfrMax > WfrMaxNew) pRadAccessData->xWfrMax = WfrMaxNew;

		WfrMinNew = TransvCenPoint.y - R;
		if(pRadAccessData->zWfrMin < WfrMinNew) pRadAccessData->zWfrMin = WfrMinNew;

		WfrMaxNew = TransvCenPoint.y + R;
		if(pRadAccessData->zWfrMax > WfrMaxNew) pRadAccessData->zWfrMax = WfrMaxNew;
	}

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
	{
		const double MinPoInAperture = 6; // To steer
		double xWfrSize = pRadAccessData->xWfrMax - pRadAccessData->xWfrMin;
		MinNx = (xWfrSize/D)*MinPoInAperture;
		double zWfrSize = pRadAccessData->zWfrMax - pRadAccessData->zWfrMin;
		MinNz = (zWfrSize/D)*MinPoInAperture;
		return 0;
	}

	//int PropagateRadMoments(srTSRWRadStructAccessData*, srTMomentsRatios*);
	//int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) 
	int CheckIfMomentsShouldBeRecomputed(double MomX_X, double MomX_Z, double MomZ_X, double MomZ_Z, double MomX_SqrtMxx_Mult, double MomX_SqrtMzz_Mult, double MomZ_SqrtMxx_Mult, double MomZ_SqrtMzz_Mult) //OC130311
	{
        double CPx_mi_R = TransvCenPoint.x - R, CPx_pl_R = TransvCenPoint.x + R;
		double CPy_mi_R = TransvCenPoint.y - R, CPy_pl_R = TransvCenPoint.y + R;

		if((((MomX_X - MomX_SqrtMxx_Mult) < CPx_mi_R) || ((MomX_X + MomX_SqrtMxx_Mult) > CPx_pl_R)) ||
		   (((MomX_Z - MomX_SqrtMzz_Mult) < CPy_mi_R) || ((MomX_Z + MomX_SqrtMzz_Mult) > CPy_pl_R)) ||
		   (((MomZ_X - MomZ_SqrtMxx_Mult) < CPx_mi_R) || ((MomZ_X + MomZ_SqrtMxx_Mult) > CPx_pl_R)) ||
		   (((MomZ_Z - MomZ_SqrtMzz_Mult) < CPy_mi_R) || ((MomZ_Z + MomZ_SqrtMzz_Mult) > CPy_pl_R))) return 1;
		else return 0;
	}
};

//*************************************************************************

class srTCircObstacle : public srTCircAperture {
public:

	srTCircObstacle(srTStringVect* pElemInfo)
	{
		D = atof((*pElemInfo)[1]); // input in m
		R = 0.5*D; Re2 = R*R;

		if(pElemInfo->size() > 2)
		{
			TransvCenPoint.x = atof((*pElemInfo)[2]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[3]); // input in m
		}

		//WfrTransmLimits.Setup(TransvCenPoint.x - R, D, TransvCenPoint.y - R, D);
	}
	srTCircObstacle(double InD, double InCx, double InCz)
	{
		D = InD;
		R = 0.5*D; Re2 = R*R;
		TransvCenPoint.x = InCx;
		TransvCenPoint.y = InCz;

		//WfrTransmLimits.Setup(TransvCenPoint.x - R, D, TransvCenPoint.y - R, D);
	}

	srTCircObstacle() {}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;

		if(TransHndl.rep != 0)
		{
			TVector3d Point(EXZ.x, 0., EXZ.z);
			FromLabToLocFrame_Point(Point);
			xRel = Point.x; zRel = Point.z;
		}
		if(xRel*xRel + zRel*zRel <= Re2)
		{
			*(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;
			*(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;
		}
	}

	void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData* pRadAccessData) {}
	//int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) { return 1;}
	int CheckIfMomentsShouldBeRecomputed(double MomX_X, double MomX_Z, double MomZ_X, double MomZ_Z, double MomX_SqrtMxx_Mult, double MomX_SqrtMzz_Mult, double MomZ_SqrtMxx_Mult, double MomZ_SqrtMzz_Mult) { return 1;} //OC130311
};

//*************************************************************************

#endif

