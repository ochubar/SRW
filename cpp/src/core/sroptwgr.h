/************************************************************************//**
 * File: sroptwgr.h
 * Description: Optical element: Rectangular Waveguide with perfectly conducting walls (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTWGR_H
#define __SROPTWGR_H

#include "sroptshp.h"
#include "sroptapt.h"

//*************************************************************************

struct srTWaveguideRectPropBufVars {

	bool PropagInFreeSpaceHoriz;
	bool PropagInFreeSpaceVert;

	double InvDxe2, InvDze2;
	double xStartTr, zStartTr;
	double xStepTr, zStepTr;
	double xFractStepTr, zFractStepTr;
	long HalfNx, HalfNz;
	//double Inv_xRangeTr, Inv_zRangeTr;
	double Inv_xStepTr, Inv_zStepTr;

	srTWaveguideRectPropBufVars()
	{
		PropagInFreeSpaceHoriz = PropagInFreeSpaceVert = false;
	}
	void SetLimitsVars(srTSRWRadStructAccessData& Wfr)
	{
		xStartTr = Wfr.xStart;
		zStartTr = Wfr.zStart;
		xStepTr = Wfr.xStep;
		zStepTr = Wfr.zStep;
		xFractStepTr = 0.1*xStepTr;
		zFractStepTr = 0.1*zStepTr;
		//Inv_xRangeTr = 1./((Wfr.nx)*(Wfr.xStep));
		//Inv_zRangeTr = 1./((Wfr.nz)*(Wfr.zStep));
		Inv_xStepTr = 1./(Wfr.xStep);
		Inv_zStepTr = 1./(Wfr.zStep);

		HalfNx = Wfr.nx >> 1;
		HalfNz = Wfr.nz >> 1;
	}
};

//*************************************************************************

class srTWaveguideRect : public srTShapedOptElem {

	srTWaveguideRectPropBufVars BufVars;

public:

	double Length;
	double Dx, Dz;

	srTWaveguideRect(srTStringVect* pElemInfo) 
	{
		Length = atof((*pElemInfo)[1]);
		Dx = atof((*pElemInfo)[2]); // input in m
		Dz = atof((*pElemInfo)[3]); // input in m

		if(pElemInfo->size() > 4)
		{
			TransvCenPoint.x = atof((*pElemInfo)[4]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[5]); // input in m
		}
		
		BufVars.PropagInFreeSpaceHoriz = false;
		BufVars.PropagInFreeSpaceVert = false;
		if(Dx <= 0) 
		{
			BufVars.PropagInFreeSpaceHoriz = true;
			BufVars.InvDxe2 = 1.E-23;
		}
		else
		{
			BufVars.InvDxe2 = 1./(Dx*Dx);
		}

		if(Dz <= 0) 
		{
			BufVars.PropagInFreeSpaceVert = true;
			BufVars.InvDze2 = 1.E-23;
		}
		else
		{
			BufVars.InvDze2 = 1./(Dz*Dz);
		}
	}

	//srTWaveguideRect(double _L, double _Dx, double _Dy, double _x =0, double _y =0) 
	srTWaveguideRect(double inL, double inDx, double inDy, double in_x=0, double in_y=0) 
	{
		Length = inL;
		Dx = inDx;
		Dz = inDy;

		TransvCenPoint.x = in_x; // input in m
		TransvCenPoint.y = in_y; // input in m

		BufVars.PropagInFreeSpaceHoriz = false;
		BufVars.PropagInFreeSpaceVert = false;
		if(Dx <= 0) 
		{
			BufVars.PropagInFreeSpaceHoriz = true;
			BufVars.InvDxe2 = 1.E-23;
		}
		else
		{
			BufVars.InvDxe2 = 1./(Dx*Dx);
		}

		if(Dz <= 0) 
		{
			BufVars.PropagInFreeSpaceVert = true;
			BufVars.InvDze2 = 1.E-23;
		}
		else
		{
			BufVars.InvDze2 = 1./(Dz*Dz);
		}
	}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResizeBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResizeBeforeAndAfterVect)
	{
		//Checks current sampling "resolution" in hor. and vert. directions
		//Makes necessary sampling for propag. through the waveguide (fit the waveguide with approx. the same resolution, include all harmonics until the cut-off)

		int result = 0;
		if(result = RemoveUndersamplingByResizingOrStop(pRadAccessData)) return result;

		//srTFringeInfo TwoSectFringeInfo[2];
		//if(result = AnalizeFringes2D(pRadAccessData, TwoSectFringeInfo)) return result;

		srTSRWRadStructAccessData AuxWfrData(*pRadAccessData);
		if(result = PrepareWavefrontForPropagation(pRadAccessData, &AuxWfrData)) return result;

		if(result = PropagateRadiationSimple_AngRepres(&AuxWfrData)) return result;

		srTRectAperture RectAp(Dx, Dz, TransvCenPoint.x, TransvCenPoint.y);
		if(result = RectAp.TraverseRadZXE(&AuxWfrData)) return result;

		if(result = CopyElecFieldDataForOut(AuxWfrData, *pRadAccessData)) return result;
		AuxWfrData.DeleteElecFieldArrays(); //deletes Ex, Ez only

		RectAp.SetNewNonZeroWfrLimits(pRadAccessData);
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		//if(result = Propagate4x4PropMatr(pRadAccessData)) return result;

		//if(result = ComputeRadMoments(pRadAccessData)) return result;
		if(result = RecomputeRadMomentsIfPossible(pRadAccessData)) return result;

		//if(MethNo > 0)
		//{
		//Consider doing: After propagation: if automatic mode, check resolution and make it the same as it was before (i.e., reduce sanmpling, if necessary)
		//}

		return 0;
	}

	int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
		pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
		pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
		double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;

		pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
		pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;

		pRadAccessData->WfrEdgeCorrShouldBeDone = 0;	
		if(pRadAccessData->Pres != 1) 
		{
			if(result = SetRadRepres(pRadAccessData, 1)) return result;
		}

			BufVars.SetLimitsVars(*pRadAccessData);

		if(result = TraverseRadZXE(pRadAccessData)) return result;

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStartTr += xShift;
				pRadAccessData->zStartTr += zShift;
			}

		if(result = SetRadRepres(pRadAccessData, 0)) return result;

		pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
				pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
			}

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!
		//double Lambda_m = 1.239854e-06/EXZ.e;
		double Lambda_m = 1.239842e-06/EXZ.e;
		if(!HarmExists(EXZ.x, EXZ.z, Lambda_m))
		{//filter-out harmonics taking into account cut-off
			*(EPtrs.pExRe) = 0; *(EPtrs.pExIm) = 0; 
			*(EPtrs.pEzRe) = 0; *(EPtrs.pEzIm) = 0; 
			return;
		}
		
		double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
		double qx2_p_qz2 = EXZ.x*EXZ.x + EXZ.z*EXZ.z;
		double c1qx2_p_qz2 = c1*qx2_p_qz2;
		double PhaseShift = c0*qx2_p_qz2*(1. + c1qx2_p_qz2 + c1qx2_p_qz2*c1qx2_p_qz2);
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}

	int HarmExists(double x, double z, double Lambda_m)
	{
		//long ix = (long)(x + BufVars.xFractStepTr - BufVars.xStartTr)*BufVars.Inv_xRangeTr;
		long ix = (long)((x + BufVars.xFractStepTr - BufVars.xStartTr)*BufVars.Inv_xStepTr);

		long xHarmNum = BufVars.HalfNx - ix;
		if(BufVars.PropagInFreeSpaceHoriz) xHarmNum = 0;

		//long iz = (long)(z + BufVars.zFractStepTr - BufVars.zStartTr)*BufVars.Inv_zRangeTr;
		long iz = (long)((z + BufVars.zFractStepTr - BufVars.zStartTr)*BufVars.Inv_zStepTr);

		long zHarmNum = BufVars.HalfNz - iz;
		if(BufVars.PropagInFreeSpaceVert) zHarmNum = 0;

		return (Lambda_m*Lambda_m*(xHarmNum*xHarmNum*BufVars.InvDxe2 + zHarmNum*zHarmNum*BufVars.InvDze2) < 4.)? 1 : 0;
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{// This is not fully correct! Try to improve (use fitting ?)

		if(BufVars.PropagInFreeSpaceHoriz)
		{
			pRadAccessData->RobsX += Length;
		}
		else
		{
			pRadAccessData->RobsX = 0.5*Length;
			pRadAccessData->RobsXAbsErr = 0.5*(pRadAccessData->RobsX);

			pRadAccessData->xc = TransvCenPoint.x;
		}
		if(BufVars.PropagInFreeSpaceVert)
		{
			pRadAccessData->RobsZ += Length;
		}
		else
		{
			pRadAccessData->RobsZ = 0.5*Length;
			pRadAccessData->RobsZAbsErr = 0.5*(pRadAccessData->RobsZ);

			pRadAccessData->zc = TransvCenPoint.y;
		}
		// RobsAbsErr not changed;
		return 0;
	}

	int PrepareWavefrontForPropagation(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pOutWfr);
	int FillInSymmetricalPartsOutsideWaveguide(srTSRWRadStructAccessData& Wfr);
	int CopyElecFieldDataForOut(srTSRWRadStructAccessData& WfrIn, srTSRWRadStructAccessData& WfrOut);


/*

	int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(result = PropagateRadiationSimple(pRadAccessData)) return result;

		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}
	int PropagateRadiationMeth_1(srTSRWRadStructAccessData*);

	void ChooseLocalPropMode(srTSRWRadStructAccessData* pRadAccessData)
	{
		int LocPropToWaistCanBeApplied = PropToWaistCanBeApplied(pRadAccessData);
		if(LocPropToWaistCanBeApplied)
		{
			if(pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ())
			{// Check necessity of this mode (if pRadAccessData is undersampled ...)
				LocalPropMode = 1;
				PropBufVars.UnderSamplingX = pRadAccessData->UnderSamplingX;
				PropBufVars.UnderSamplingZ = pRadAccessData->UnderSamplingZ;

				return;
			}
		}
		if((!LocPropToWaistCanBeApplied) && (pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ()))
		{
			LocalPropMode = -1;
			return;
		}

		LocalPropMode = 0; // Normal, through ang. repres.
	}

	int PropToWaistCanBeApplied(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(!AllowPropToWaist) return 0;

		const double DistWavelengthFactor = 5; // To steer
		double LambdaM = 3.1415926535898/(pRadAccessData->eStart*2.53384080189E+06);
		double CritDist = LambdaM*DistWavelengthFactor;
		if((::fabs(pRadAccessData->RobsX) < CritDist) || (::fabs(pRadAccessData->RobsZ) < CritDist)) return 0;

		double NewRobsX = pRadAccessData->RobsX + Length;
		double NewRobsZ = pRadAccessData->RobsZ + Length;
		char ToWaistX = (::fabs(NewRobsX) < 0.6*(::fabs(pRadAccessData->RobsX))); // To steer
		char ToWaistZ = (::fabs(NewRobsZ) < 0.6*(::fabs(pRadAccessData->RobsZ))); // To steer
		if(!(ToWaistX || ToWaistZ)) return 0;

		// Check more...

		 return 1;
	}
	int PropStatPhaseCanBeApplied(srTSRWRadStructAccessData* pRadAccessData)
	{
		const double GoodRelPrecR = 0.15; // To steer

		char HorRadOK = (pRadAccessData->RobsXAbsErr < ::fabs(pRadAccessData->RobsX));
		char VertRadOK = (pRadAccessData->RobsZAbsErr < ::fabs(pRadAccessData->RobsZ));
		return (HorRadOK && VertRadOK);
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(LocalPropMode == 0) return PropagateRadiationSimple_AngRepres(pRadAccessData);
		else if(LocalPropMode == 1) return PropagateRadiationSimple_PropToWaist(pRadAccessData);
		else return 0;
	}
	int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
		pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
		pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
		double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;

		pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
		pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;

			pRadAccessData->WfrEdgeCorrShouldBeDone = 0;	

		if(pRadAccessData->Pres != 1) 
		{
			if(result = SetRadRepres(pRadAccessData, 1)) return result;
		}
		if(result = TraverseRadZXE(pRadAccessData)) return result;

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStartTr += xShift;
				pRadAccessData->zStartTr += zShift;
			}

		if(result = SetRadRepres(pRadAccessData, 0)) return result;

		pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
				pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
			}

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
		return 0;
	}
	int PropagateRadiationSimple_PropToWaist(srTSRWRadStructAccessData* pRadAccessData);
	int PropagateRadiationSimple_PropFromWaist(srTSRWRadStructAccessData* pRadAccessData);

	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		if(LocalPropMode == 0) return PropagateRadiationSimple1D_AngRepres(pSect1D);
		else if(LocalPropMode == 1) return PropagateRadiationSimple1D_PropToWaist(pSect1D);
		else return 0;
	}
	int PropagateRadiationSimple1D_AngRepres(srTRadSect1D* pSect1D)
	{
		int result;
		double OldStart = pSect1D->ArgStart;
		pSect1D->ArgStart = -(pSect1D->np >> 1)*pSect1D->ArgStep;
		double Shift = pSect1D->ArgStart - OldStart;

		pSect1D->WfrMin += Shift;
		pSect1D->WfrMax += Shift;

		if(pSect1D->Pres != 1) 
		{
			if(result = SetRadRepres1D(pSect1D, 1)) return result;
		}
		if(result = TraverseRad1D(pSect1D)) return result;
		if(result = SetRadRepres1D(pSect1D, 0)) return result;

		pSect1D->ArgStart = OldStart;
		pSect1D->SetNonZeroWavefrontLimitsToFullRange();
		return 0;
	}
	int PropagateRadiationSimple1D_PropToWaist(srTRadSect1D* pSect1D);

	int ResizeBeforePropToWaistIfNecessary(srTSRWRadStructAccessData* pRadAccessData);
	void EstimateMinNxNzBeforePropToWaist(srTSRWRadStructAccessData* pRadAccessData, int& Nx, int& Nz);
	int ResizeBeforePropToWaistIfNecessary1D(srTRadSect1D* pSect1D);
	void EstimateMinNpBeforePropToWaist1D(srTRadSect1D* pSect1D, int& Np);

	void SetupPropBufVars_PropToWaist(srTSRWRadStructAccessData* pRadAccessData)
	{// Compute any buf. vars for Stat. Phase propagation in PropStatPhaseBufVars
		PropBufVars.InvLength = 1./Length;
		double Pi_d_LambdaM = pRadAccessData->eStart*2.53384080189E+06;
		PropBufVars.Pi_d_LambdaM_d_Length = Pi_d_LambdaM*(PropBufVars.InvLength);
		PropBufVars.InvLength_d_Lambda = PropBufVars.InvLength*pRadAccessData->eStart*806546.577258;
		PropBufVars.xc = pRadAccessData->xc;
		PropBufVars.zc = pRadAccessData->zc;
		PropBufVars.ExtraConstPhase = Pi_d_LambdaM*(pRadAccessData->xc*pRadAccessData->xc/pRadAccessData->RobsX 
												+ pRadAccessData->zc*pRadAccessData->zc/pRadAccessData->RobsZ);

		double TwoPi_d_LambdaM = 2.*Pi_d_LambdaM;
		PropBufVars.TwoPiXc_d_LambdaMRx = TwoPi_d_LambdaM*pRadAccessData->xc/pRadAccessData->RobsX;
		PropBufVars.TwoPiZc_d_LambdaMRz = TwoPi_d_LambdaM*pRadAccessData->zc/pRadAccessData->RobsZ;

		// Continue for more buf vars
	}

	void SetupPropBufVars_PropFromWaist(srTSRWRadStructAccessData* pRadAccessData)
	{// Compute any buf. vars for Stat. Phase propagation in PropStatPhaseBufVars
		PropBufVars.InvLength = 1./Length;
		double Pi_d_LambdaM = pRadAccessData->eStart*2.53384080189E+06;
		PropBufVars.Pi_d_LambdaM_d_Length = Pi_d_LambdaM*(PropBufVars.InvLength);
		PropBufVars.InvLength_d_Lambda = PropBufVars.InvLength*pRadAccessData->eStart*806546.577258;

		//PropBufVars.xc = pRadAccessData->xc;
		//PropBufVars.zc = pRadAccessData->zc;
		//PropBufVars.ExtraConstPhase = Pi_d_LambdaM*(pRadAccessData->xc*pRadAccessData->xc/pRadAccessData->RobsX 
		//										+ pRadAccessData->zc*pRadAccessData->zc/pRadAccessData->RobsZ);

		//double TwoPi_d_LambdaM = 2.*Pi_d_LambdaM;
		//PropBufVars.TwoPiXc_d_LambdaMRx = TwoPi_d_LambdaM*pRadAccessData->xc/pRadAccessData->RobsX;
		//PropBufVars.TwoPiZc_d_LambdaMRz = TwoPi_d_LambdaM*pRadAccessData->zc/pRadAccessData->RobsZ;

		// Continue for more buf vars
	}

	void SetupPropBufVars_PropToWaist(srTRadSect1D* pSect1D)
	{// Compute any buf. vars for Stat. Phase propagation in PropStatPhaseBufVars
		PropBufVars.InvLength = 1./Length;
		double Pi_d_LambdaM = pSect1D->eVal*2.53384080189E+06;
		PropBufVars.Pi_d_LambdaM_d_Length = Pi_d_LambdaM*(PropBufVars.InvLength);
		PropBufVars.InvLength_d_Lambda = PropBufVars.InvLength*pSect1D->eVal*806546.577258;
		double TwoPi_d_LambdaM = 2.*Pi_d_LambdaM;
		PropBufVars.TwoPiXc_d_LambdaMRx = TwoPi_d_LambdaM*pSect1D->cArg/pSect1D->Robs;

		// Continue for more buf vars
	}

	int PropagateRadiationTest(srTSRWRadStructAccessData* pInRadAccessData, srTSRWRadStructAccessData* pOutRadAccessData)
	{
		// Develop here: Propagate by computing the integral directly
		pOutRadAccessData->AuxLong1 = 1;
		return 0;
	}
	int TuneRadForPropMeth_1(srTSRWRadStructAccessData*, srTRadResize&);

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		if(LocalPropMode == 0) { RadPointModifier_AngRepres(EXZ, EPtrs); return;}
		else if(LocalPropMode == 1) { RadPointModifier_PropToWaist(EXZ, EPtrs); return;}
		else if(LocalPropMode == 2) { RadPointModifier_PropFromWaist(EXZ, EPtrs); return;}
	}
	void RadPointModifier_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!
		double Lambda_m = 1.239842e-06/EXZ.e;
		double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
		double qx2_p_qz2 = EXZ.x*EXZ.x + EXZ.z*EXZ.z;
		double c1qx2_p_qz2 = c1*qx2_p_qz2;
		double PhaseShift = c0*qx2_p_qz2*(1. + c1qx2_p_qz2 + c1qx2_p_qz2*c1qx2_p_qz2);
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}
	void RadPointModifier_PropToWaist(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		double rx = EXZ.x, rz = EXZ.z;
		double PhaseShift = PropBufVars.Pi_d_LambdaM_d_Length*(rx*rx + rz*rz);
		if(PropBufVars.PassNo == 1) PhaseShift += PropBufVars.TwoPiXc_d_LambdaMRx*rx + PropBufVars.TwoPiZc_d_LambdaMRz*rz;
		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);

		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
		
		if(PropBufVars.PassNo == 2)
		{
			NewExRe = (*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda;
			NewExIm = (*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda;
			NewEzRe = (*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda; 
			NewEzIm = (*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda; 
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}
	void RadPointModifier_PropFromWaist(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		double rx = EXZ.x, rz = EXZ.z;
		double PhaseShift = PropBufVars.Pi_d_LambdaM_d_Length*(rx*rx + rz*rz);
		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);

		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 

		if(PropBufVars.PassNo == 2)
		{
			NewExRe = (*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda;
			NewExIm = (*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda;
			NewEzRe = (*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda; 
			NewEzIm = (*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda; 
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}

	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		if(LocalPropMode == 0) { RadPointModifier1D_AngRepres(EXZ, EPtrs); return; }
		else if(LocalPropMode == 1) { RadPointModifier1D_PropToWaist(EXZ, EPtrs); return; }
	}
	void RadPointModifier1D_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!
		double Lambda_m = 1.239842e-06/EXZ.e;
		double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
		double Arg = (EXZ.VsXorZ == 'x')? EXZ.x : EXZ.z;
		double ArgE2 = Arg*Arg;
		double c1ArgE2 = c1*ArgE2;
		double PhaseShift = c0*ArgE2*(1. + c1ArgE2 + c1ArgE2*c1ArgE2);
		float CosPh, SinPh;
		CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}
	void RadPointModifier1D_PropToWaist(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
		double rx = EXZ.x, rz = EXZ.z;
		double PhaseShift = PropBufVars.Pi_d_LambdaM_d_Length*(rx*rx + rz*rz);

		if(PropBufVars.PassNo == 1)
		{
			if(EXZ.VsXorZ == 'x') PhaseShift += PropBufVars.TwoPiXc_d_LambdaMRx*rx;
			else PhaseShift += PropBufVars.TwoPiXc_d_LambdaMRx*rz;
		}
		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
		
		if(PropBufVars.PassNo == 2)
		{
			NewExRe = (*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda; // Change Norm for 1D !!!
			NewExIm = (*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda; // Change Norm for 1D
			NewEzRe = (*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda; // Change Norm for 1D
			NewEzIm = (*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda; // Change Norm for 1D
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		float aStr0[] = { 1., Length };
		float aStr1[] = { 0., 1. };
		float* a[] = { aStr0, aStr1 };
		return AuxPropagateRadMoments(pRadAccessData, a, a, MomRatArray);
	}
	int AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, float** ax, float** az, srTMomentsRatios* MomRatArray);

	int PropagateWaveFrontRadius1D(srTRadSect1D* pSect1D)
	{// This is not fully correct: It does not take into account diffraction...
		pSect1D->Robs += Length;
		pSect1D->SetNonZeroWavefrontLimitsToFullRange();
		// RobsAbsErr not changed;
		// xc, zc not changed
		return 0;
	}

	int Propagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData) 
	{
		double Drift4x4Matr[] = { 1., Length, 0., 0.,
								  0., 1., 0., 0.,
								  0., 0., 1., Length,
								  0., 0., 0., 1.};
		double Drift4Vect[] = { 0., 0., 0., 0.};

		return GenAuxPropagate4x4PropMatr(pRadAccessData, Drift4x4Matr, Drift4Vect);
	}

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
	{
		const double MinPo = 40; // To steer
		MinNx = MinNz = MinPo;
		return 0;
	}

	double DiffractionLimitedPropagatedSpotSize(char x_or_z, srTSRWRadStructAccessData* pRadAccessData, long ie)
	{
		double PhotEn = pRadAccessData->eStart + ie*pRadAccessData->eStep;
		double Lambda_m = 1.239842e-06/PhotEn;
		double dLim = (x_or_z == 'x')? (pRadAccessData->nx - 1)*pRadAccessData->xStep : (pRadAccessData->nz - 1)*pRadAccessData->zStep;
		double dWfr = (x_or_z == 'x')? (pRadAccessData->xWfrMax - pRadAccessData->xWfrMin) : (pRadAccessData->zWfrMax - pRadAccessData->zWfrMin);
		double d = ((dWfr > 0.) && (dWfr < dLim))? dWfr : dLim;
		return Lambda_m*Length/d;
	}

	int RangeShouldBeAdjustedAtPropag() { return 1;}
	int ResolutionShouldBeAdjustedAtPropag() { return 0;} // Or switch it Off
*/
};

//*************************************************************************

#endif
