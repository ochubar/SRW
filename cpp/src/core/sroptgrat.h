/************************************************************************//**
 * File: sroptgrat.h
 * Description: Optical element: Grating (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTGRAT_H
#define __SROPTGRAT_H

#include "sroptshp.h"

//*************************************************************************

struct srTGratingPropBufVars {
	double L2, ThetaM0, ThetaM, ThetaIt, SinThetaI;
	double Sin_ThetaM0_p_ThetaIt, Cos_ThetaM0_p_ThetaIt, Tg_ThetaIt;
	double td_NominTermL2, td_NominMultX2, td_MultInvDenom;
	double OptPathCorTiltMultX2;
	double Sin_ThetaM_p_ThetaIt, Cos_ThetaM_p_ThetaIt;
	double CurPhotEn, CurWaveNumb, ReflectAmp;
	double wfrR, Lambda;
	bool WaveFrontTermWasTreated;
	//double ConstRxE, ConstRzE;
	double AnamorphMagn, PowerConservMultE;

	srTGratingPropBufVars()
	{
		AnamorphMagn = 1.;
	}
};

//*************************************************************************

class srTGrating : public srTShapedOptElem {
//For the moment, plane grating only (to be extended?)

	double m_Period; //Grating Period [m]
	int m_Order; 
	double m_ReflectAvgInt;

	//Polynomial coefficients for groove density, so that the Groove Density in [lines/m] (units to check/refine!?) equals to:
	//m_grDen + m_grDen1*s + m_grDen2*s^2 + m_grDen3*s^3 + m_grDen4*s^4
	//where s is longitudinal position along grating in [m]
	double m_grDen, m_grDen1, m_grDen2, m_grDen3, m_grDen4;  

	srTGratingPropBufVars m_PropBufVars;

public:

	srTGrating(srTStringVect* pElemInfo)
	{
		//m_PropWfrInPlace = false; //previous electric field is necessary for the propagation
		m_PropWfrInPlace = true; //OC151008 //previous electric field is NOT necessary for the propagation

		if(pElemInfo == 0) return;

		double GrooveDens = atof((*pElemInfo)[1]); //expected in [lines/mm]
		m_Period = 1.e-03/GrooveDens;
		
		RotPlane = atoi((*pElemInfo)[2]); //treat as Dispersion (Deflection) Plane; can be 'h' or 'v'
		if(RotPlane == 1) RotPlane = 'h';
		else RotPlane = 'v';

		Theta = HalfPI - (atof((*pElemInfo)[3]))*PI/180.; //in line with definition in base class
		m_Order = atoi((*pElemInfo)[4]); 
		m_ReflectAvgInt = atof((*pElemInfo)[5]); 
	}
	//srTGrating(double _grDen, char _disPl, double _ang, int _m, double _refl)
	srTGrating(double _grDen, char _disPl, double _ang, int _m, double _refl, double _grDen1=0, double _grDen2=0, double _grDen3=0, double _grDen4=0)
	{
		m_Period = 1.e-03/_grDen;

		if((_disPl == 'x') || (_disPl == 'h')) RotPlane = 'h';
		else if((_disPl == 'y') || (_disPl == 'v')) RotPlane = 'v';

		Theta = HalfPI - _ang; //in line with definition in base class
		m_Order = _m;
		m_ReflectAvgInt = _refl;

		m_grDen = _grDen*1e+03; //[lines/m]?
		m_grDen1 = _grDen1*1e+06; //[lines/m^2]?
		m_grDen2 = _grDen2*1e+09; //[lines/m^3]?
		m_grDen3 = _grDen3*1e+12; //[lines/m^4]?
		m_grDen4 = _grDen4*1e+15; //[lines/m^5]?
	}
	srTGrating() 
	{
		//m_PropWfrInPlace = false; //previous electric field is necessary for the propagation
		m_PropWfrInPlace = true; //OC151008 //previous electric field is NOT necessary for the propagation
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		//char &MethNo = ParPrecWfrPropag.MethNo;
		SetupPropBufVars_Gen(pRadAccessData);

		return PropagateRadiationMeth_0(pRadAccessData);

		//to program all methods later!

		//if(MethNo == 0) return PropagateRadiationMeth_0(pRadAccessData);
		//else if(MethNo == 1) return PropagateRadiationMeth_1(pRadAccessData);
		//else if(MethNo == 2) return PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
		//return 0;
	}
	int PropagateRadiationSingleE_Meth_0_old(srTSRWRadStructAccessData* pWfr, srTSRWRadStructAccessData* pPrevWfr)
	{
		int result = 0;
		if(pPrevWfr == 0) return 0; //to fire error?

		m_pPrevWfr = pPrevWfr;
		if(m_pPrevWfr->Pres != 0) if(result = SetRadRepres(m_pPrevWfr, 0)) return result;
		if(pWfr->Pres != 0) if(result = SetRadRepres(pWfr, 0)) return result;

		SetupPropBufVars_SingleE(pWfr->eStart);

		if(result = TraverseRadZXE(pWfr)) return result;

		//To simulate anamorphic magnification in the functions:
		//if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pWfr)) return result;
		//if(result = Propagate4x4PropMatr(pRadAccessData)) return result;

		pWfr->SetNonZeroWavefrontLimitsToFullRange();
		if(m_PropBufVars.WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(*m_pPrevWfr, 'a');
		m_pPrevWfr = 0;
		return 0;
	}
	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pWfr, srTSRWRadStructAccessData* pPrevWfr)
	{//this version doesn't use pPrevWfr
	 //however, it may modify me

		int result = 0;
		//if(pPrevWfr == 0) return 0;
		//m_pPrevWfr = pPrevWfr;
		//if(m_pPrevWfr->Pres != 0) if(result = SetRadRepres(m_pPrevWfr, 0)) return result;

		m_pPrevWfr = 0;
		if(pWfr->Pres != 0) if(result = SetRadRepres(pWfr, 0)) return result; //ensure coordinate repres.

		SetupPropBufVars_SingleE(pWfr->eStart); //check if all these buf. var. are necessary
		AdjustWfrMeshParamToTreatAnamorphMagn(pWfr);

		if(result = TraverseRadZXE(pWfr)) return result;

		//To simulate anamorphic magnification in the functions:
		//if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pWfr)) return result;
		//if(result = Propagate4x4PropMatr(pRadAccessData)) return result;

		pWfr->SetNonZeroWavefrontLimitsToFullRange();
		//if(m_PropBufVars.WaveFrontTermWasTreated) TreatStronglyOscillatingTerm(*m_pPrevWfr, 'a');
		return 0;
	}

	void SetupPropBufVars_Gen(srTSRWRadStructAccessData* pWfr);
	void SetupPropBufVars_SingleE(double PhotEn);
	
	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{//Photon energy in eV; Length in m !!! Operates on Coordinate side !!!
		//This version takes into account varying incidence angle over the input wavefront
		//and it treats the amplitude change due to anamorphous magnification
		//It is assumed that the wavefront was already "stretched" (mesh changed) due to the anamorphous magnification
		//"x->-x" is not taken into account

		if(EXZ.e != m_PropBufVars.CurPhotEn) SetupPropBufVars_SingleE(EXZ.e);

		double x2 = (RotPlane == 'h')? EXZ.x : EXZ.z;

		double x1 = -x2/m_PropBufVars.AnamorphMagn;
		double instThetaI = Theta;
		if(m_PropBufVars.wfrR != 0) instThetaI += x1/m_PropBufVars.wfrR;

		double instThetaM = asin(m_Order*m_PropBufVars.Lambda/m_Period - sin(instThetaI));

		//Simplified version: just an angle (dependent on x1!)
		//double angDisp = instThetaM - m_PropBufVars.ThetaM0;
		double angDisp = (instThetaM - m_PropBufVars.ThetaM0) + (instThetaI - Theta); //to check!!!
		//double angDisp = -(instThetaM - m_PropBufVars.ThetaM0) - (instThetaI - Theta); //to check!!!

		double phaseShift = m_PropBufVars.CurWaveNumb*angDisp*x2;
		float cosPh, sinPh;
		CosAndSin(phaseShift, cosPh, sinPh);
		float NewExRe = (float)(((*(EPtrs.pExRe))*cosPh - (*(EPtrs.pExIm))*sinPh)*m_PropBufVars.PowerConservMultE);
		float NewExIm = (float)(((*(EPtrs.pExRe))*sinPh + (*(EPtrs.pExIm))*cosPh)*m_PropBufVars.PowerConservMultE);
		float NewEzRe = (float)(((*(EPtrs.pEzRe))*cosPh - (*(EPtrs.pEzIm))*sinPh)*m_PropBufVars.PowerConservMultE);
		float NewEzIm = (float)(((*(EPtrs.pEzRe))*sinPh + (*(EPtrs.pEzIm))*cosPh)*m_PropBufVars.PowerConservMultE);
		//to check/correct aventual modification of different component phases at reflection!

		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{//This is not fully correct: It does not take into account diffraction...

		if((m_PropBufVars.AnamorphMagn != 0) && (m_PropBufVars.AnamorphMagn != 1))
		{
			if(RotPlane == 'h') pRadAccessData->RobsX *= m_PropBufVars.AnamorphMagn;
			else if(RotPlane == 'v') pRadAccessData->RobsZ *= m_PropBufVars.AnamorphMagn;
		}
		// RobsAbsErr not changed;
		// xc, zc not changed
		return 0;
	}

	void AdjustWfrMeshParamToTreatAnamorphMagn(srTSRWRadStructAccessData* pWfr)
	{//To call only after SetupPropBufVars_SingleE
		if(pWfr == 0) return;
		double magnTol = 1.e-5;
		if(::fabs(m_PropBufVars.AnamorphMagn - 1.) < magnTol) return;

		double *pStart = &(pWfr->xStart), *pStep = &(pWfr->xStep);
		long *pN = &(pWfr->nx);
		if(RotPlane == 'v')
		{
			pStart = &(pWfr->zStart); pStep = &(pWfr->zStep);
			pN = &(pWfr->nz);
		}

		//almost same actions as in Resize
		double oldRange = ((*pN) - 1)*(*pStep);
		double oldEnd = (*pStart) + oldRange;
		double mid = 0.5*((*pStart) + oldEnd);
		double newRange = m_PropBufVars.AnamorphMagn*oldRange;
		*pStep = (*pN > 1)? newRange/((*pN) - 1) : 0;
		*pStart = mid - 0.5*newRange;
		//Amplitude correction due to the anamorph. magn. should be made in RadPointModifier
	}

	//int PropagateRadMoments(srTSRWRadStructAccessData*, srTMomentsRatios*);

	//virtual void SetNewNonZeroWfrLimits(srTSRWRadStructAccessData*) {}

	//int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	//{
	//	int result;
	//	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
	//	if(result = TraverseRadZXE(pRadAccessData)) return result;

	//	SetNewNonZeroWfrLimits(pRadAccessData);
	//	return 0;
	//}
	//int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	//{
	//	int result;
	//	if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
	//	if(result = TraverseRad1D(pSect1D)) return result;
	//	return 0;
	//}

	//int RangeShouldBeAdjustedAtPropag() { return 1;}
	//int ResolutionShouldBeAdjustedAtPropag() { return 0;}

	//virtual int CheckIfMomentsShouldBeRecomputed(float MomX_X, float MomX_Z, float MomZ_X, float MomZ_Z, float MomX_SqrtMxx_Mult, float MomX_SqrtMzz_Mult, float MomZ_SqrtMxx_Mult, float MomZ_SqrtMzz_Mult) { return 1;}
};

//*************************************************************************

#endif
