/************************************************************************//**
 * File: srgsnbm.h
 * Description: Gaussian beam (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRGSNBM_H
#define __SRGSNBM_H

#include "srobject.h"
#include "srebmdat.h"
#include "srstraux.h"

//*************************************************************************

class srTGsnBeam : public CGenObject {

	double LongDist;
	double NormConstElField;
	double InvTwoSigXe2, InvTwoSigZe2, m_InvTwoSigTe2, m_InvTwoSigPhotEnE2;
	double PropagMultX, PropagMultZ;
	double x0Prop, z0Prop;

	double m_NormConstElecFldT; //pre-exponential factor for time-domain electric field
	//which ensures intensity in [W/mm^2] given by squared amplitude of the electric field

	double m_AvgPhotEn_enMult_PropInvRx_T, m_AvgPhotEn_enMult_PropInvRz_T, m_PhaseLongDelay_T;
	double m_PropInvTwoSigXe2_T, m_PropInvTwoSigZe2_T;
	double m_PropInvSigX_T, m_PropInvSigZ_T, m_NormConstElecFldT_Prop;

public:

	srTEbmDat EbmDat;

	double SigmaX, SigmaZ, SigmaT, SigmaPhotEn; //these are RMS values with respect to Electric Field, and not the Intensity!
	int mx, mz, Polar, TypeT;
	double PhotPerBW;
	double RefractIndex;
	double m_RepRate, m_PulseEn, m_AvgPhotEn;

	srTWfrSmp DistrInfoDat;

	srTGsnBeam(double SpecFlux, int Polar, double sigX, int mx, double sigZ, int mz, double sigT, int typeT, double* pMom1, double s0, double RepRate=0, double PulseEn=0, double AvgPhotEn=0);
	srTGsnBeam()
	{
		Initialize();
	}
	void Initialize()
	{
		RefractIndex = 1.;
		m_RepRate = m_PulseEn = m_AvgPhotEn = 0;
		m_NormConstElecFldT = 0;
	}

	int CheckInputConsistency();
	void SetupSourceConstantsFreqDomain();
	void SetupSourceConstantsTimeDomain();
	int CreateWavefrontElFieldFreqDomain(srTSRWRadStructAccessData&);
	int CreateWavefrontElFieldTimeDomain(srTSRWRadStructAccessData&);

	//void ComputeElectricFieldFreqDomain(srTWfrSmp* pWfrSmp, srTSRWRadStructAccessData* pWfr);
	void ComputeElectricField(srTWfrSmp* pWfrSmp, srTSRWRadStructAccessData* pWfr);

	double Factorial(long n);
	double HermitePolynomial(int n, double x);

	void SetupProperPolariz(double ReA, double ImA, float* tRadX, float* tRadZ)
	{// Same for Isotropic Source and Gaussian Beam
		const double c = 0.70710678118655;
		switch(Polar) 
		{
		case 1: // Lin. Hor.
			*tRadX = (float)ReA; *(tRadX+1) = (float)ImA; *tRadZ = 0.; *(tRadZ+1) = 0.;
			break;
		case 2: // Lin. Vert.
			*tRadX = 0.; *(tRadX+1) = 0.; *tRadZ = (float)ReA; *(tRadZ+1) = (float)ImA;
			break;
		case 3: // Lin. 45
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(c*ReA); *(tRadZ+1) = (float)(c*ImA);
			break;
		case 4: // Lin. 135
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(-c*ReA); *(tRadZ+1) = (float)(-c*ImA);
			break;
		case 5: // Circ. Right
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(-c*ImA); *(tRadZ+1) = (float)(c*ReA);
			break;
		case 6: // Circ. Left
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(c*ImA); *(tRadZ+1) = (float)(-c*ReA);
			break;
		}
	}
};

//*************************************************************************

#endif
