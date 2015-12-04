/************************************************************************//**
 * File: sroptfoc.h
 * Description: Optical element: Focusing (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTFOC_H
#define __SROPTFOC_H

#include "sroptshp.h"

//*************************************************************************

class srTFocusingElem : public srTShapedOptElem {

protected:
	bool m_wfrRadWasProp;

public:
	double FocDistX, FocDistZ;
	srTFocusingElem() {}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		//float aStr0[] = { 1., 0. };
		//float axStr1[] = { (float)(-1./FocDistX), 1. };
		//float azStr1[] = { (float)(-1./FocDistZ), 1. };
		//float* ax[] = { aStr0, axStr1 };
		//float* az[] = { aStr0, azStr1 };
		double aStr0[] = { 1., 0. }; //OC130311
		double axStr1[] = { (-1./FocDistX), 1. };
		double azStr1[] = { (-1./FocDistZ), 1. };
		double* ax[] = { aStr0, axStr1 };
		double* az[] = { aStr0, azStr1 };
		return GenAuxPropagateRadMoments(pRadAccessData, ax, az, MomRatArray);
	}
	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{
		//double MagnX = FocDistX/(FocDistX - pRadAccessData->RobsX);
		//double MagnZ = FocDistZ/(FocDistZ - pRadAccessData->RobsZ);

		const double dSmall = 1E-23;
		double difX = (FocDistX == pRadAccessData->RobsX)? dSmall : FocDistX - pRadAccessData->RobsX;
		double difZ = (FocDistZ == pRadAccessData->RobsZ)? dSmall : FocDistZ - pRadAccessData->RobsZ;

		double MagnX = FocDistX/difX;
		double MagnZ = FocDistZ/difZ;

		//pRadAccessData->RobsX = pRadAccessData->RobsX*MagnX;
		//pRadAccessData->RobsZ = pRadAccessData->RobsZ*MagnZ;
		pRadAccessData->RobsX *= MagnX;
		pRadAccessData->RobsZ *= MagnZ;

		pRadAccessData->RobsXAbsErr *= (MagnX*MagnX);
		pRadAccessData->RobsZAbsErr *= (MagnZ*MagnZ);

		pRadAccessData->xc = (pRadAccessData->xc - TransvCenPoint.x)*MagnX;
		pRadAccessData->zc = (pRadAccessData->zc - TransvCenPoint.y)*MagnZ;

		m_wfrRadWasProp = true;
		return 0;
	}
	int PropagateWaveFrontRadius1D(srTRadSect1D* pSect1D) 
	{
		const double dSmall = 1E-23;
		if(pSect1D->VsXorZ == 'x')
		{
			//double Magn = FocDistX/(FocDistX - pSect1D->Robs);

			double dif = (FocDistX == pSect1D->Robs)? dSmall : FocDistX - pSect1D->Robs;
			double Magn = FocDistX/dif;

			pSect1D->Robs = pSect1D->Robs*Magn;
			pSect1D->RobsAbsErr *= (Magn*Magn);
			pSect1D->cArg = (pSect1D->cArg - TransvCenPoint.x)*Magn;
		}
		else
		{
			//double Magn = FocDistZ/(FocDistZ - pSect1D->Robs);

			double dif = (FocDistZ == pSect1D->Robs)? dSmall : FocDistZ - pSect1D->Robs;
			double Magn = FocDistZ/dif;

			pSect1D->Robs = pSect1D->Robs*Magn;
			pSect1D->RobsAbsErr *= (Magn*Magn);
			pSect1D->cArg = (pSect1D->cArg - TransvCenPoint.y)*Magn;
		}
		return 0;
	}

	int Propagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData) 
	{
		double Foc4x4Matr[] = { 1., 0., 0., 0.,
							   -1./FocDistX, 1., 0., 0.,
								0., 0., 1., 0.,
								0., 0., -1./FocDistZ, 1.};
		double Foc4Vect[] = { 0., TransvCenPoint.x/FocDistX, 0., TransvCenPoint.y/FocDistZ };

		return GenAuxPropagate4x4PropMatr(pRadAccessData, Foc4x4Matr, Foc4Vect);
	}

	//int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadDataSingleE)
	{
		int result = 0;
		//if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		//if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		//if(result = PropagateRadiationSimple(pRadAccessData)) return result;
		//if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		if(result = PropagateRadiationSimple(pRadAccessData)) return result; //in first place because previous wavefront radius may be required for some derived classes
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int TuneRadForPropMeth_1(srTSRWRadStructAccessData*, srTRadResize&);
};

//*************************************************************************

class srTThinLens : public srTFocusingElem {
public:

	srTThinLens(srTStringVect* pElemInfo) 
	{
		FocDistX = atof((*pElemInfo)[1]); // input in m
		FocDistZ = atof((*pElemInfo)[2]); // input in m

		if(pElemInfo->size() > 3)
		{
			TransvCenPoint.x = atof((*pElemInfo)[3]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[4]); // input in m
		}
	}
	srTThinLens(double InFocDistX, double InFocDistZ, double InCx, double InCz)
	{
		FocDistX = InFocDistX; FocDistZ = InFocDistZ; 
		TransvCenPoint.x = InCx; TransvCenPoint.y = InCz;
	}
	srTThinLens() {}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			pRadAccessData->CheckAndSubtractPhaseTermsLin(TransvCenPoint.x, TransvCenPoint.y);
		//}

		char &MethNo = ParPrecWfrPropag.MethNo;
		int result = 0;

		if(MethNo == 0) result = PropagateRadiationMeth_0(pRadAccessData);
		else if(MethNo == 1) result = PropagateRadiationMeth_1(pRadAccessData);
		else if(MethNo == 2) result = PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);

		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			if(!ParPrecWfrPropag.DoNotResetAnalTreatTermsAfterProp) pRadAccessData->CheckAndResetPhaseTermsLin();
		//}

		return result;
	}
	int PropagateRadiationMeth_1(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		srTRadResize PostResize;
		PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;

		if(result = TuneRadForPropMeth_1(pRadAccessData, PostResize)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;

		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;

		//const double ResizeTol = 0.15;
		char PostResizeNeeded = (::fabs(PostResize.pxm - 1.) || ::fabs(PostResize.pzm - 1.) || ::fabs(PostResize.pxd - 1.) || ::fabs(PostResize.pzd - 1.));
		if(PostResizeNeeded) if(result = RadResizeGen(*pRadAccessData, PostResize)) return result;

		if(result = ComputeRadMoments(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;
		return 0;
	}
  	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		int result;
		if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
		if(result = TraverseRad1D(pSect1D)) return result;
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		//double Pi_d_Lambda_m = EXZ.e*2.533840802E+06;
		double Pi_d_Lambda_m = EXZ.e*2.533865612E+06;
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;

		double PhaseShift = -Pi_d_Lambda_m*(xRel*xRel/FocDistX + zRel*zRel/FocDistZ);
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}

  	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double Pi_d_Lambda_m = EXZ.e*2.533840802E+06;
		double ArgRel = (EXZ.VsXorZ == 'x')? (EXZ.x - TransvCenPoint.x) : (EXZ.z - TransvCenPoint.y);
		double FocDist = (EXZ.VsXorZ == 'x')? FocDistX : FocDistZ;

		double PhaseShift = -Pi_d_Lambda_m*(ArgRel*ArgRel/FocDist);
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}

	double RadOptPathDiff(srTEXZ& EXZ) //virtual
	{
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
		return -0.5*(xRel*xRel/FocDistX + zRel*zRel/FocDistZ);
	}

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
	{
		const double MinPo = 40; // To steer
		MinNx = MinNz = MinPo;
		return 0;
	}

	int RangeShouldBeAdjustedAtPropag() { return 0;} // Or switch it On
	int ResolutionShouldBeAdjustedAtPropag() { return 1;}
};

//*************************************************************************

#endif
