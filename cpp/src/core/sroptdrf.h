/************************************************************************//**
 * File: sroptdrf.h
 * Description: Optical element: Drift space (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTDRF_H
#define __SROPTDRF_H

#ifndef __SROPTELM_H
#include "sroptelm.h"
#endif
#ifndef __SRERROR_H
#include "srerror.h"
#endif

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include <stdio.h>
//#include "srwlib.h"

//*************************************************************************

struct srTDriftPropBufVars {
	int PassNo;

	double Pi_d_LambdaM_d_Length;
	double InvLength;
	double InvLength_d_Lambda;
	double xc, zc;
	double ExtraConstPhase;
	double invRx, invRz;
	double Lx, Lz, invRxL, invRzL, sqrt_LxLz_d_L, phase_term_signLxLz;
	double Pi_d_LambdaM_d_Rx, Pi_d_LambdaM_d_Rz;
	double kx_AnalytTreatQuadPhaseTerm, kxc_AnalytTreatQuadPhaseTerm, kz_AnalytTreatQuadPhaseTerm, kzc_AnalytTreatQuadPhaseTerm;

	double TwoPiXc_d_LambdaMRx, TwoPiZc_d_LambdaMRz;
	double UnderSamplingX, UnderSamplingZ;
	bool UseExactRxRzForAnalytTreatQuadPhaseTerm;
	char AnalytTreatSubType; //OC24042013
	srTDriftPropBufVars()
	{
		UnderSamplingX = UnderSamplingZ = 1.;
		UseExactRxRzForAnalytTreatQuadPhaseTerm = false;
		AnalytTreatSubType = 0;
	}
};

//*************************************************************************

class srTDriftSpace : public srTGenOptElem {

	char AllowPropToWaist; // To remove

	char LocalPropMode; // -1- abort; 
						// 0- normal (through angular repres.);
						// 1- prop. to waist;
						// 2- prop. from waist
						// 3- prop. without quad. phase term
	//srTDriftPropBufVars PropBufVars;
	char TreatPath; // switch specifying whether the absolute optical path should be taken into account in radiation phase (=1) or not (=0, default)

public:
	double Length;
	srTDriftPropBufVars PropBufVars;

	srTDriftSpace(double InLength =0., char InTreatPath =0) 
	{ 
		Length = InLength;
		TreatPath = InTreatPath; //OC010813
		AllowPropToWaist = 1; // To switch
	}
	srTDriftSpace(srTStringVect* pElemInfo) 
	{ 
		Length = atof((*pElemInfo)[1]);
		AllowPropToWaist = 1; // To switch
	}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResizeBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResizeBeforeAndAfterVect)
	{
		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//double start;
		//get_walltime(&start);

		int result = 0;
		
		ChooseLocalPropMode(pRadAccessData, ParPrecWfrPropag);
		if(LocalPropMode == -1)
		{
			double GoodNx = pRadAccessData->nx*pRadAccessData->UnderSamplingX;
			double GoodNz = pRadAccessData->nz*pRadAccessData->UnderSamplingZ;

			if(result = TryToRemoveUndersamplingByResizing(pRadAccessData)) return result;
			if(pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ()) return PROP_CAN_NOT_BE_DONE_DUE_TO_MEMORY_LIMIT;
			else LocalPropMode = 0;

			if((GoodNx*0.7 > double(pRadAccessData->nx)) || (GoodNz*0.7 > double(pRadAccessData->nz)))
			{// To steer
				CErrWarn::AddWarningMessage(&gVectWarnNos, PROPAG_PREC_REDUCED_DUE_TO_MEMORY_LIMIT);
			}
		}

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiation : LocalPropMode == -1",&start);

		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
		//OC25102010: commented-out because of errors in case of partially-coherent emission and B fiber
		//	pRadAccessData->CheckAndSubtractPhaseTermsLin(pRadAccessData->GetWfrMiddleHor(), pRadAccessData->GetWfrMiddleVer());
		//}
		//return result; //test

		char &MethNo = ParPrecWfrPropag.MethNo;

		if(MethNo == 0) result = PropagateRadiationMeth_0(pRadAccessData);
		else if(MethNo == 1) result = PropagateRadiationMeth_1(pRadAccessData);
		else if(MethNo == 2) result = PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResizeBeforeAndAfterVect);
		
		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
		//OC25102010: commented-out because of errors in case of partially-coherent emission and B fiber
		//	if(!ParPrecWfrPropag.DoNotResetAnalTreatTermsAfterProp) pRadAccessData->CheckAndResetPhaseTermsLin();
		//}

		return result;
	}

	//int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{//it works for many photon energies too!
		int result;
		if(result = PropagateRadiationSimple(pRadAccessData)) return result;
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData) //virtual in srTGenOptElem
	{//because for the Drift, the following works for many photon energies too!
		//return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0);
		//OC251214
		if((LocalPropMode == 0) || (LocalPropMode == 3) || (pRadAccessData->ne == 1)) return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0);
		else
		{
			pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
			return srTGenOptElem::PropagateRadiationMeth_0(pRadAccessData); //since (LocalPropMode == 1) and (LocalPropMode == 2) - propagation to/from waist introduces some dispersion and can potentially modify mesh
		}
	}

	int PropagateRadiationMeth_1(srTSRWRadStructAccessData*);

	void ChooseLocalPropMode(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag)
	{
		//OC test
		//	LocalPropMode = 3; return;
		//end OC test

/**
//OC240114 (commented-out)
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
			else LocPropToWaistCanBeApplied = 0; //test 03012007
		}

		PropBufVars.AnalytTreatSubType = ParPrecWfrPropag.AnalTreatment; //OC24042013

		int GoodCondForPropWithoutQuadTerm = ((ParPrecWfrPropag.AnalTreatment != 0) && (!LocPropToWaistCanBeApplied)); // && //test 03012007
			//(pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ() || 
			//(!pRadAccessData->CheckIfSamplingResolvesElecFieldWithQuadPhaseTerm()));
		//if(GoodCondForPropWithoutQuadTerm && (!ParPrecWfrPropag.UseResBefore) && (!ParPrecWfrPropag.UseResAfter))
		if(GoodCondForPropWithoutQuadTerm && (!ParPrecWfrPropag.UseResBefore))
		{//propagation without quad. term is temporary not allowed in Automatic Mode
			LocalPropMode = 3; return;
		}

		if((!LocPropToWaistCanBeApplied) && (!GoodCondForPropWithoutQuadTerm) && (pRadAccessData->ThereIsUnderSamplingX() || pRadAccessData->ThereIsUnderSamplingZ()))
		{
			LocalPropMode = -1; return;
		}
**/

		LocalPropMode = 0; // Normal, through ang. repres.

		//OC240114
		if((ParPrecWfrPropag.AnalTreatment == 1) || (ParPrecWfrPropag.AnalTreatment == 2)) 
		{
			LocalPropMode = 3; PropBufVars.AnalytTreatSubType = ParPrecWfrPropag.AnalTreatment;
		}
		else if(ParPrecWfrPropag.AnalTreatment == 3) LocalPropMode = 2; //Propagation From Waist
		else if(ParPrecWfrPropag.AnalTreatment == 4) LocalPropMode = 1; //Propagation To Waist

		//OC100914 Aux. methods for testing / benchmarking
		else if(ParPrecWfrPropag.AnalTreatment >= 100) LocalPropMode = 100; //Propagation To Waist
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
		char ToWaistX = (::fabs(NewRobsX) < 0.3*(::fabs(pRadAccessData->RobsX))); // To steer //0.6
		char ToWaistZ = (::fabs(NewRobsZ) < 0.3*(::fabs(pRadAccessData->RobsZ))); // To steer //0.6
		//if(!(ToWaistX || ToWaistZ)) return 0;
		if((!ToWaistX) || (!ToWaistZ)) return 0; //OCfix 170307

		// Check more...

		 return 1;
	}
	int PropStatPhaseCanBeApplied(srTSRWRadStructAccessData* pRadAccessData)
	{
		//const double GoodRelPrecR = 0.15; // To steer

		char HorRadOK = (pRadAccessData->RobsXAbsErr < ::fabs(pRadAccessData->RobsX));
		char VertRadOK = (pRadAccessData->RobsZAbsErr < ::fabs(pRadAccessData->RobsZ));
		return (HorRadOK && VertRadOK);
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(LocalPropMode == 0) return PropagateRadiationSimple_AngRepres(pRadAccessData);
		else if(LocalPropMode == 1) return PropagateRadiationSimple_PropToWaist(pRadAccessData);
		else if(LocalPropMode == 2) return PropagateRadiationSimple_PropFromWaist(pRadAccessData); //OC240114 (added)
		else if(LocalPropMode == 3) return PropagateRadiationSimple_AnalytTreatQuadPhaseTerm(pRadAccessData);

		//OC100914 Aux. methods for testing / benchmarking
		else if(LocalPropMode == 100) return PropagateRadiationSimple_NumIntFresnel(pRadAccessData);

		else return 0;
	}
	int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRadAccessData)
	{
		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//double start;
		//get_walltime(&start);

		int result;
		double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
		pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
		pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
		double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;

		pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
		pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;

			pRadAccessData->WfrEdgeCorrShouldBeDone = 0;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiationSimple_AngRepres:setup",&start);

		if(pRadAccessData->Pres != 1) 
		{
			if(result = SetRadRepres(pRadAccessData, 1)) return result;
		}

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiationSimple_AngRepres:SetRadRepres 1",&start);

		if(result = TraverseRadZXE(pRadAccessData)) return result;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiationSimple_AngRepres:TraverseRadZXE",&start);

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStartTr += xShift;
				pRadAccessData->zStartTr += zShift;
			}

		if(result = SetRadRepres(pRadAccessData, 0)) return result;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiationSimple_AngRepres:SetRadRepres 2",&start);

		pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;

			if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
			{
				pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
				pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
			}

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":PropagateRadiationSimple_AngRepres:SetNonZeroWavefrontLimitsToFullRange 2",&start);

		return 0;
	}
	int PropagateRadiationSimple_PropToWaist(srTSRWRadStructAccessData* pRadAccessData);
	int PropagateRadiationSimple_PropFromWaist(srTSRWRadStructAccessData* pRadAccessData);
	int PropagateRadiationSimple_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData);
	int PropagateRadiationSimple_NumIntFresnel(srTSRWRadStructAccessData* pRadAccessData); //OC100914 Aux. method for testing / benchmarking

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
	
	//void CheckAndSubtractPhaseTermsLin(srTSRWRadStructAccessData* pRadAccessData);
	//void CheckAndResetPhaseTermsLin(srTSRWRadStructAccessData* pRadAccessData);
	
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

	void SetupPropBufVars_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData);
	void EstimateTrueWfrRadAndMaxLeff_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& trueRx, double& trueRz, double& Lx_eff_max, double& Lz_eff_max);
	void EstimateWfrRadToSub_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& effRx, double& effRz);
	void EstimateWfrRadToSub2_AnalytTreatQuadPhaseTerm(srTSRWRadStructAccessData* pRadAccessData, double& effRx, double& effRz);

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
		else if(LocalPropMode == 3) { RadPointModifier_AnalytTreatQuadPhaseTerm(EXZ, EPtrs); return;}
	}
	void RadPointModifier_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!
		//double Lambda_m = 1.239854e-06/EXZ.e;
		double Lambda_m = 1.239842e-06/EXZ.e;
		double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
		double qx2_p_qz2 = EXZ.x*EXZ.x + EXZ.z*EXZ.z;
		double c1qx2_p_qz2 = c1*qx2_p_qz2;
		double PhaseShift = c0*qx2_p_qz2*(1. + c1qx2_p_qz2 + c1qx2_p_qz2*c1qx2_p_qz2);

		////double ConstPhaseShift = 6.2831853071795*Length/Lambda_m;
		////PhaseShift += ConstPhaseShift; //may be necessary for 3D wavefronts?
		if(TreatPath == 1) //OC010813
		{//??
			PhaseShift += (5.067730652e+06)*Length*EXZ.e;
			//+= 6.2831853072*Length/Lambda_m;
		}

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

		//if(PropBufVars.PassNo == 3) 
		//{
		//	//if(*(EPtrs.pExRe) != 0)
		//	//{
		//	//	int aha = 1;
		//	//}
		//	return;
		//}

		if(PropBufVars.PassNo == 1) 
		{
			PhaseShift += PropBufVars.TwoPiXc_d_LambdaMRx*rx + PropBufVars.TwoPiZc_d_LambdaMRz*rz;

			if(TreatPath == 1) //OC010813
			{//??
				PhaseShift += (5.067730652e+06)*Length*EXZ.e;
				//+= 6.2831853072*Length/Lambda_m;
			}
		}
		////else if(PropBufVars.PassNo == 2)
		////{
		////	double Lambda_m = 1.239842e-06/EXZ.e;
		////	double ConstPhaseShift = 6.2831853071795*Length/Lambda_m;
		////	PhaseShift += ConstPhaseShift; //may be necessary for 3D wavefronts?
		////}
		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);

		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
		
		if(PropBufVars.PassNo == 2)
		{
			NewExRe = (float)((*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda);
			NewExIm = (float)((*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda);
			NewEzRe = (float)((*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda); 
			NewEzIm = (float)((*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda); 
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}
	void RadPointModifier_PropFromWaist(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		double rx = EXZ.x, rz = EXZ.z;
		double PhaseShift = PropBufVars.Pi_d_LambdaM_d_Length*(rx*rx + rz*rz);

		////double Lambda_m = 1.239842e-06/EXZ.e;
		////double ConstPhaseShift = 6.2831853071795*Length/Lambda_m;
		////PhaseShift += ConstPhaseShift; //may be necessary for 3D wavefronts?

		if(TreatPath == 1) //OC010813
		{
			if(PropBufVars.PassNo == 2)//??
			{
				PhaseShift += (5.067730652e+06)*Length*EXZ.e;
				//+= 6.2831853072*Length/Lambda_m;
			}
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
			NewExRe = (float)((*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda);
			NewExIm = (float)((*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda);
			NewEzRe = (float)((*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda); 
			NewEzIm = (float)((*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda); 
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}
	void RadPointModifier_AnalytTreatQuadPhaseTerm(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{//don't use RobsX, RobsZ directly here!

		double PhaseShift = 0;
		//double Lambda_m = 1.239854e-06/EXZ.e;
		double Lambda_m = 1.239842e-06/EXZ.e;
		if(PropBufVars.PassNo == 1) //removing quad. term from the phase on coordinate side
		{
			double rx = EXZ.x - PropBufVars.xc, rz = EXZ.z - PropBufVars.zc;
			double Pi_d_Lambda_m = 3.1415926536/Lambda_m;
			PhaseShift = -Pi_d_Lambda_m*(PropBufVars.invRx*rx*rx + PropBufVars.invRz*rz*rz);
		}
		else if(PropBufVars.PassNo == 2) //loop on angular side
		{// e in eV; Length in m !!! Operates on Angles side !!!
			double Pi_Lambda_m = 3.1415926536*Lambda_m;
			double c0x = -Pi_Lambda_m*PropBufVars.Lx, c0z = -Pi_Lambda_m*PropBufVars.Lz; //c1 = 0.25*Lambda_m*Lambda_m;
			
			//double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
			//double qx2_p_qz2 = EXZ.x*EXZ.x + EXZ.z*EXZ.z;
			//double c1qx2_p_qz2 = c1*qx2_p_qz2;
			//double PhaseShift = c0*qx2_p_qz2*(1. + c1qx2_p_qz2 + c1qx2_p_qz2*c1qx2_p_qz2);
			PhaseShift = c0x*EXZ.x*EXZ.x + c0z*EXZ.z*EXZ.z + PropBufVars.phase_term_signLxLz;
		}
		else if(PropBufVars.PassNo == 3) //adding new quad. term to the phase on coordinate side
		{
			double rx = EXZ.x - PropBufVars.xc, rz = EXZ.z - PropBufVars.zc;
			//double rx = EXZ.x /*- PropBufVars.xc*/, rz = EXZ.z; /*- PropBufVars.zc;*/
			double Pi_d_Lambda_m = 3.1415926536/Lambda_m;
			PhaseShift = Pi_d_Lambda_m*(PropBufVars.invRxL*rx*rx + PropBufVars.invRzL*rz*rz);

			//OCTEST
			//PhaseShift += 2*Pi_d_Lambda_m*Length; //may be necessary for 3D wavefronts?
			//END OCTEST
			if(TreatPath == 1) //OC010813
			{
				PhaseShift += 2*Pi_d_Lambda_m*Length;
			}
		}

		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);
		float NewExRe = (*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh;
		float NewExIm = (*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh;
		float NewEzRe = (*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh;

		if(PropBufVars.PassNo == 2)
		{
			double& multL = PropBufVars.sqrt_LxLz_d_L;
			//NewExRe *= multL;
			//NewExIm *= multL;
			//NewEzRe *= multL;
			//NewEzIm *= multL;
			NewExRe = (float)(NewExRe*multL);
			NewExIm = (float)(NewExIm*multL);
			NewEzRe = (float)(NewEzRe*multL);
			NewEzIm = (float)(NewEzIm*multL);
		}

		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}

	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		if(LocalPropMode == 0) { RadPointModifier1D_AngRepres(EXZ, EPtrs); return; }
		else if(LocalPropMode == 1) { RadPointModifier1D_PropToWaist(EXZ, EPtrs); return; }
	}
	void RadPointModifier1D_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!
		//double Lambda_m = 1.239854e-06/EXZ.e;
		double Lambda_m = 1.239842e-06/EXZ.e;
		double c0 = -3.1415926536*Length*Lambda_m, c1 = 0.25*Lambda_m*Lambda_m;
		double Arg = (EXZ.VsXorZ == 'x')? EXZ.x : EXZ.z;
		double ArgE2 = Arg*Arg;
		double c1ArgE2 = c1*ArgE2;
		double PhaseShift = c0*ArgE2*(1. + c1ArgE2 + c1ArgE2*c1ArgE2);
		if(TreatPath == 1) //OC010813
		{//??
			PhaseShift += (5.067730652e+06)*Length*EXZ.e;
			//+= 6.2831853072*Length/Lambda_m;
		}

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

			if(TreatPath == 1) //OC010813
			{//??
				PhaseShift += (5.067730652e+06)*Length*EXZ.e;
				//+= 6.2831853072*Length/Lambda_m;
			}
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
			NewExRe = (float)((*(EPtrs.pExRe))*PropBufVars.InvLength_d_Lambda); // Change Norm for 1D !!!
			NewExIm = (float)((*(EPtrs.pExIm))*PropBufVars.InvLength_d_Lambda); // Change Norm for 1D
			NewEzRe = (float)((*(EPtrs.pEzRe))*PropBufVars.InvLength_d_Lambda); // Change Norm for 1D
			NewEzIm = (float)((*(EPtrs.pEzIm))*PropBufVars.InvLength_d_Lambda); // Change Norm for 1D
			
			*(EPtrs.pExRe) = NewExIm; *(EPtrs.pExIm) = -NewExRe;
			*(EPtrs.pEzRe) = NewEzIm; *(EPtrs.pEzIm) = -NewEzRe;
		}
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		//float aStr0[] = { 1., (float)Length };
		//float aStr1[] = { 0., 1. };
		//float* a[] = { aStr0, aStr1 };
		double aStr0[] = { 1., Length };
		double aStr1[] = { 0., 1. };
		double* a[] = { aStr0, aStr1 };
		return AuxPropagateRadMoments(pRadAccessData, a, a, MomRatArray);
	}
	//int AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, float** ax, float** az, srTMomentsRatios* MomRatArray);
	int AuxPropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, double** ax, double** az, srTMomentsRatios* MomRatArray); //OC130311
	int PropagateElecBeamMoments(srTEbmDat* pEbm);

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{// This is not fully correct: It does not take into account diffraction...
	 // RobsX, RobsZ are treated here as distances to waists 

		pRadAccessData->RobsX += Length;
		pRadAccessData->RobsZ += Length;

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
		// RobsAbsErr not changed;
		// xc, zc not changed
		return 0;
	}
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
		//double Lambda_m = 1.239854e-06/PhotEn;
		double Lambda_m = 1.239842e-06/PhotEn;
		double dLim = (x_or_z == 'x')? (pRadAccessData->nx - 1)*pRadAccessData->xStep : (pRadAccessData->nz - 1)*pRadAccessData->zStep;
		double dWfr = (x_or_z == 'x')? (pRadAccessData->xWfrMax - pRadAccessData->xWfrMin) : (pRadAccessData->zWfrMax - pRadAccessData->zWfrMin);
		double d = ((dWfr > 0.) && (dWfr < dLim))? dWfr : dLim;
		return Lambda_m*Length/d;
	}

	int RangeShouldBeAdjustedAtPropag() { return 1;}
	int ResolutionShouldBeAdjustedAtPropag() { return 0;} // Or switch it Off
	int AllowAutoSwitchToUndersamplingMode() { return 1;} //0;}
};

//*************************************************************************

#endif
