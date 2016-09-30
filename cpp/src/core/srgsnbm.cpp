/************************************************************************//**
 * File: srgsnbm.cpp
 * Description: Gaussian beam
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srgsnbm.h"
#include "sroptelm.h"

//*************************************************************************

srTGsnBeam::srTGsnBeam(double SpecFlux, int InPolar, double InSigX_intens_m, int In_mx, double InSigZ_intens_m, int In_mz, double InSigT_intens_s, int InTypeT, double* pMom1, double s0, double RepRate_Hz, double PulseEn_J, double AvgPhotEn_eV)
{
	PhotPerBW = SpecFlux;
    Polar = InPolar;

	const double sqrt_2 = sqrt(2.);
	const double Pi = 3.1415926535898;
	const double sqrt_Pi = sqrt(Pi);
	const double sqrt_TwoPi = sqrt(2*Pi);

	SigmaX = InSigX_intens_m*sqrt_2;
	mx = In_mx;
    SigmaZ = InSigZ_intens_m*sqrt_2;
	mz = In_mz;
    SigmaT = InSigT_intens_s*sqrt_2;
	TypeT = InTypeT;

	SigmaPhotEn = 1.e+23;
	if(InSigT_intens_s > 0)
	{
		SigmaPhotEn = (0.5*6.58211889E-16/InSigT_intens_s)*sqrt_2;
	}

	m_RepRate = RepRate_Hz; 
	m_PulseEn = PulseEn_J;
	m_AvgPhotEn = AvgPhotEn_eV;

	m_NormConstElecFldT = 0;
	if((PulseEn_J > 0) && (InSigT_intens_s > 0) && (InSigX_intens_m > 0) && (InSigZ_intens_m > 0))
	{
		const double TwoPi_pow3d2 = sqrt_TwoPi*sqrt_TwoPi*sqrt_TwoPi;
		m_NormConstElecFldT = 0.001*sqrt(PulseEn_J/(TwoPi_pow3d2*InSigX_intens_m*InSigZ_intens_m*InSigT_intens_s));
	}
	
	if(SpecFlux <= 0)
	{//Calculate spectral flux from pulse energy
		if((PulseEn_J > 0) && (InSigT_intens_s > 0))
		{
			PhotPerBW = 6.2415096E+15*PulseEn_J/(sqrt_Pi*SigmaPhotEn); //to get fluence in [Ph/mm^2/0.1%bw]
			if(RepRate_Hz > 0) PhotPerBW *= RepRate_Hz; //to get intensity in [Ph/s/mm^2/0.1%bw]
		}
	}

	double AuxConst = double(1 << (mx + mz))*Factorial(mx)*Factorial(mz); //moved here from SetupSourceConstantsFreqDomain
	//double ConstSpecFluxPerUnSurf = PhotPerBW/(AuxConst*((1.E+06)*Pi));
	NormConstElField = sqrt(PhotPerBW/(AuxConst*((1.E+06)*Pi))); //moved here from SetupSourceConstantsFreqDomain

	double ElecMom1[] = {1., 0., 0., 0., 0.};
	if(pMom1 != 0) 
	{
		for(int i=0; i<4; i++) ElecMom1[i + 1] = pMom1[i];
	}

	double InNeb = 0;
    srTEbmDat LocEbmDat(1., InNeb, ElecMom1, 5, 0, 0, s0);
	EbmDat = LocEbmDat;

	RefractIndex = 1;
}

//*************************************************************************

int srTGsnBeam::CheckInputConsistency()
{
	if((DistrInfoDat.LambStart == 0.) || (DistrInfoDat.LambEnd == 0.)) return PHOT_EN_SHOULD_BE_POSITIVE;
// Fill-in more
	return 0;
}

//*************************************************************************

double srTGsnBeam::Factorial(long n)
{
	if(n == 0) return 1.;
	else return n*Factorial(n - 1);
}

//*************************************************************************

double srTGsnBeam::HermitePolynomial(int n, double x)
{
	if(n == 0) return 1.;
	if(n == 1) return 2.*x;
	return 2.*(x*HermitePolynomial(n - 1, x) - (n - 1)*HermitePolynomial(n - 2, x));
}

//*************************************************************************

void srTGsnBeam::SetupSourceConstantsFreqDomain()
{
	const double Pi = 3.1415926535898;
	//const double PhEnConv = 1.239854E-06;
	const double PhEnConv = 1.239842E-06;

	LongDist = DistrInfoDat.yStart - EbmDat.s0;
	if(LongDist == 0.)
	{// Artificial shift of observation position to avoid singularity
		double Wavelength = PhEnConv/DistrInfoDat.LambStart;
		LongDist = 0.01*Wavelength;
		DistrInfoDat.yStart = EbmDat.s0 + LongDist;
	}

// To give Phot/s/0.1%bw/mm^2
// assuming EbmDat.Current in A, PhotPerBW in Phot/0.1%bw, SigmaX in m
	//double AuxConst = double(1 << (mx + mz))*Factorial(mx)*Factorial(mz); //moved to constructor
	//double ConstSpecFluxPerUnSurf = EbmDat.Current*PhotPerBW/(AuxConst*((1.602189246E-13)*Pi));
	
	//double ConstSpecFluxPerUnSurf = PhotPerBW/(AuxConst*((1.E+06)*Pi)); //moved to constructor
	//NormConstElField = sqrt(ConstSpecFluxPerUnSurf); //moved to constructor

	InvTwoSigXe2 = 0.5/(SigmaX*SigmaX);
	InvTwoSigZe2 = 0.5/(SigmaZ*SigmaZ);
	m_InvTwoSigPhotEnE2 = 0.5/(SigmaPhotEn*SigmaPhotEn);

	double Buf = (Pi*RefractIndex)/(PhEnConv*LongDist);
	PropagMultX = Buf/InvTwoSigXe2;
	PropagMultZ = Buf/InvTwoSigZe2;

	x0Prop = EbmDat.x0 + EbmDat.dxds0*LongDist;
	z0Prop = EbmDat.z0 + EbmDat.dzds0*LongDist;
}

//*************************************************************************

void srTGsnBeam::SetupSourceConstantsTimeDomain()
{
	const double Pi = 3.1415926535898;
	//const double PhEnConv = 1.239854E-06;
	const double PhEnConv = 1.239842E-06;

	LongDist = DistrInfoDat.yStart - EbmDat.s0;
	if(LongDist == 0.)
	{// Artificial shift of observation position to avoid singularity
		double Wavelength = PhEnConv/DistrInfoDat.LambStart;
		LongDist = 0.01*Wavelength;
		DistrInfoDat.yStart = EbmDat.s0 + LongDist;
	}

// To give W/mm^2
// assuming EbmDat.Current in A, PhotPerBW in Phot/0.1%bw, SigmaX in m
	double AuxConstHighModes = double(1 << (mx + mz))*Factorial(mx)*Factorial(mz);
	////double ConstSpecFluxPerUnSurf = EbmDat.Current*PhotPerBW/(AuxConst*((1.602189246E-13)*Pi));
	//double ConstSpecFluxPerUnSurf = PhotPerBW/(AuxConst*((1.E+06)*Pi));
	//NormConstElField = sqrt(ConstSpecFluxPerUnSurf);

	InvTwoSigXe2 = 0.5/(SigmaX*SigmaX);
	InvTwoSigZe2 = 0.5/(SigmaZ*SigmaZ);
	m_InvTwoSigTe2 = (SigmaT == 0)? 0 : 0.5/(SigmaT*SigmaT);

	double Buf = (Pi*RefractIndex)/(PhEnConv*LongDist);
	PropagMultX = Buf/InvTwoSigXe2;
	PropagMultZ = Buf/InvTwoSigZe2;

	const double enMult = 2.53384080189E+06;

	double InvLongDist = 1./LongDist;

	double PropRatX_T = m_AvgPhotEn*PropagMultX, PropRatZ_T = m_AvgPhotEn*PropagMultZ;
	double PropRatXe2_T = PropRatX_T*PropRatX_T, PropRatZe2_T = PropRatZ_T*PropRatZ_T;
	double DisMultX_T = 1. + PropRatXe2_T, DisMultZ_T = 1. + PropRatZe2_T;
	double PropInvRx_T = InvLongDist/DisMultX_T, PropInvRz_T = InvLongDist/DisMultZ_T;
	m_AvgPhotEn_enMult_PropInvRx_T = m_AvgPhotEn*enMult*PropInvRx_T;
	m_AvgPhotEn_enMult_PropInvRz_T = m_AvgPhotEn*enMult*PropInvRz_T;

	double InvPropRatX_T = 1./PropRatX_T, InvPropRatZ_T = 1./PropRatZ_T;
	double SigMultX_T = 1. + InvPropRatX_T*InvPropRatX_T, SigMultZ_T = 1. + InvPropRatZ_T*InvPropRatZ_T;
	m_PropInvTwoSigXe2_T = InvTwoSigXe2/SigMultX_T;
	m_PropInvTwoSigZe2_T = InvTwoSigZe2/SigMultZ_T;
	m_PropInvSigX_T = sqrt(2.*m_PropInvTwoSigXe2_T); 
	m_PropInvSigZ_T = sqrt(2.*m_PropInvTwoSigZe2_T);
	m_NormConstElecFldT_Prop = m_NormConstElecFldT/sqrt(sqrt(SigMultX_T*SigMultZ_T)*AuxConstHighModes);

	double NuX_T = atan(InvPropRatX_T), NuZ_T = atan(InvPropRatZ_T);
	m_PhaseLongDelay_T = (mx + 0.5)*NuX_T + (mz + 0.5)*NuZ_T;

	x0Prop = EbmDat.x0 + EbmDat.dxds0*LongDist;
	z0Prop = EbmDat.z0 + EbmDat.dzds0*LongDist;
}

//*************************************************************************

//void srTGsnBeam::ComputeElectricFieldFreqDomain(srTWfrSmp* pWfrSmp, srTSRWRadStructAccessData* pWfr)
void srTGsnBeam::ComputeElectricField(srTWfrSmp* pWfrSmp, srTSRWRadStructAccessData* pWfr)
{// m, eV !
	if((pWfrSmp == 0) || (pWfr == 0)) throw INCORRECT_PARAMS_SR_COMP;

	DistrInfoDat = *pWfrSmp;
	DistrInfoDat.EnsureZeroTransverseRangesForSinglePoints();
	pWfr->SetRadSamplingFromObs(DistrInfoDat);

	int res = 0;
	if(pWfr->PresT)
	{
		if(res = CreateWavefrontElFieldTimeDomain(*pWfr)) throw res;
	}
	else
	{
		if(res = CreateWavefrontElFieldFreqDomain(*pWfr)) throw res;
	}

	pWfr->SetNonZeroWavefrontLimitsToFullRange();

	srTGenOptElem GenOptElem;
	if(res = GenOptElem.ComputeRadMoments(pWfr)) throw res;
}

//*************************************************************************

int srTGsnBeam::CreateWavefrontElFieldFreqDomain(srTSRWRadStructAccessData& RadAccessData)
{// m, eV or s!
//can create frequency- or time-domain electric field
	int result;
	const double TwoPI = 6.28318530717959;
	const double InvTwoPI = 1./TwoPI;
	const double enMult = 2.53384080189E+06;

	RadAccessData.SetAvgPhotEnergyFromLimits(); //OC180314

	if(result = CheckInputConsistency()) return result;
	SetupSourceConstantsFreqDomain();

	//x0Prop = EbmDat.x0 + EbmDat.dxds0*LongDist;
	//z0Prop = EbmDat.z0 + EbmDat.dzds0*LongDist;

	double ActNormConstElField = NormConstElField; //OC081014
	if(RadAccessData.ElecFldUnit == 2) ActNormConstElField *= sqrt(1.602176462e-16); //case of field units: sqrt(J/eV/mm^2) or sqrt(W/mm^2)

	double InvLongDist = 1./LongDist;

	float* tRadX = RadAccessData.pBaseRadX;
	float* tRadZ = RadAccessData.pBaseRadZ;

	double z = RadAccessData.zStart - z0Prop;
	double z_mi_z0 = RadAccessData.zStart - EbmDat.z0; //OC210413

	for(int iz=0; iz<RadAccessData.nz; iz++)
	{
		double ze2 = z*z;
		double z_mi_z0_e2 = z_mi_z0*z_mi_z0; //OC210413
		double zz = z_mi_z0 + EbmDat.z0; //OC210413

		double x = RadAccessData.xStart - x0Prop;
		double x_mi_x0 = RadAccessData.xStart - EbmDat.x0; //OC210413

		for(int ix=0; ix<RadAccessData.nx; ix++)
		{
			double xe2 = x*x;
			double x_mi_x0_e2 = x_mi_x0*x_mi_x0; //OC210413
			double xx = x_mi_x0 + EbmDat.x0; //OC210413

			double en = RadAccessData.eStart;
			for(int ie=0; ie<RadAccessData.ne; ie++)
			{
				double PropRatX = en*PropagMultX, PropRatZ = en*PropagMultZ;
				double PropRatXe2 = PropRatX*PropRatX, PropRatZe2 = PropRatZ*PropRatZ;
				double InvPropRatX = 1./PropRatX, InvPropRatZ = 1./PropRatZ;

				double SigMultX = 1. + InvPropRatX*InvPropRatX, SigMultZ = 1. + InvPropRatZ*InvPropRatZ;
				double DisMultX = 1. + PropRatXe2, DisMultZ = 1. + PropRatZe2;

				double PropInvTwoSigXe2 = InvTwoSigXe2/SigMultX, PropInvTwoSigZe2 = InvTwoSigZe2/SigMultZ;
				double PropInvSigX = sqrt(2.*PropInvTwoSigXe2), PropInvSigZ = sqrt(2.*PropInvTwoSigZe2);
				double PropInvRx = InvLongDist/DisMultX, PropInvRz = InvLongDist/DisMultZ;

				double xpMult = (1./PropInvRx - LongDist)*EbmDat.dxds0, zpMult = (1./PropInvRz - LongDist)*EbmDat.dzds0; //OC210413

				double NuX = atan(InvPropRatX), NuZ = atan(InvPropRatZ);

				//double Phase = en*enMult*(xe2*PropInvRx + ze2*PropInvRz) + (mx + 0.5)*NuX + (mz + 0.5)*NuZ;
				double Phase = en*enMult*((x_mi_x0_e2 + 2*LongDist*EbmDat.x0*EbmDat.dxds0 + xpMult*(2*xx - LongDist*EbmDat.dxds0))*PropInvRx 
										+ (z_mi_z0_e2 + 2*LongDist*EbmDat.z0*EbmDat.dzds0 + zpMult*(2*zz - LongDist*EbmDat.dzds0))*PropInvRz) 
							   + (mx + 0.5)*NuX + (mz + 0.5)*NuZ; //OC210413 //??

				Phase -= TwoPI*long(Phase*InvTwoPI);
				double CosPh = cos(Phase), SinPh = sin(Phase);

				//double ExpAr = exp(-xe2*PropInvTwoSigXe2 - ze2*PropInvTwoSigZe2);
				double argForExp = -xe2*PropInvTwoSigXe2 - ze2*PropInvTwoSigZe2;
				if(m_AvgPhotEn > 0) 
				{
					double dPhotEn = en - m_AvgPhotEn;
					argForExp -= dPhotEn*dPhotEn*m_InvTwoSigPhotEnE2;
				}

				double ExpAr = exp(argForExp);
				double HermX = HermitePolynomial(mx, x*PropInvSigX), HermZ = HermitePolynomial(mz, z*PropInvSigZ);

				//double BufA = NormConstElField*sqrt(PropInvSigX*PropInvSigZ)*ExpAr*HermX*HermZ;
				double BufA = ActNormConstElField*sqrt(PropInvSigX*PropInvSigZ)*ExpAr*HermX*HermZ; //OC081014
				double ReA = BufA*CosPh, ImA = BufA*SinPh;

				SetupProperPolariz(ReA, ImA, tRadX, tRadZ);

				tRadX += 2; tRadZ += 2;
				en += RadAccessData.eStep;
			}
			x += RadAccessData.xStep;
			x_mi_x0 += RadAccessData.xStep; //OC210413
		}
		z += RadAccessData.zStep;
		z_mi_z0 += RadAccessData.zStep; //OC210413
	}
	RadAccessData.Pres = 0;
	RadAccessData.PresT = 0;
	return 0;
}

//*************************************************************************

int srTGsnBeam::CreateWavefrontElFieldTimeDomain(srTSRWRadStructAccessData& RadAccessData)
{// m, eV or s!
//can create frequency- or time-domain electric field
	int result;

	RadAccessData.avgPhotEn = m_AvgPhotEn; 

	const double TwoPI = 6.28318530717959;
	const double InvTwoPI = 1./TwoPI;
	//const double enMult = 2.53384080189E+06;

	if(result = CheckInputConsistency()) return result;
	SetupSourceConstantsTimeDomain();

	//double InvLongDist = 1./LongDist;

	float* tRadX = RadAccessData.pBaseRadX;
	float* tRadZ = RadAccessData.pBaseRadZ;

	double z = RadAccessData.zStart - z0Prop;
	for(int iz=0; iz<RadAccessData.nz; iz++)
	{
		double ze2 = z*z;

		double x = RadAccessData.xStart - x0Prop;
		for(int ix=0; ix<RadAccessData.nx; ix++)
		{
			double xe2 = x*x;

			double Phase_T = m_AvgPhotEn_enMult_PropInvRx_T*xe2 + m_AvgPhotEn_enMult_PropInvRz_T*ze2 + m_PhaseLongDelay_T;
			//Phase_T -= TwoPI*long(Phase_T*InvTwoPI);
			Phase_T -= TwoPI*((long long)(Phase_T*InvTwoPI));
			double CosPh_T = cos(Phase_T), SinPh_T = sin(Phase_T);
			double argForExp0_T = -xe2*m_PropInvTwoSigXe2_T - ze2*m_PropInvTwoSigZe2_T;

			double HermX_T = HermitePolynomial(mx, x*m_PropInvSigX_T);
			double HermZ_T = HermitePolynomial(mz, z*m_PropInvSigZ_T);

			double t = RadAccessData.eStart; //en is actually time here
			for(int it=0; it<RadAccessData.ne; it++)
			{
				double argForExp_T = argForExp0_T - t*t*m_InvTwoSigTe2;
				double ExpAr_T = exp(argForExp_T);
				double BufA_T = m_NormConstElecFldT_Prop*ExpAr_T*HermX_T*HermZ_T;
				double ReA_T = BufA_T*CosPh_T, ImA_T = BufA_T*SinPh_T;
				SetupProperPolariz(ReA_T, ImA_T, tRadX, tRadZ);

				tRadX += 2; tRadZ += 2;
				t += RadAccessData.eStep;
			}
			x += RadAccessData.xStep;
		}
		z += RadAccessData.zStep;
	}
	RadAccessData.Pres = 0;
	RadAccessData.PresT = 1;
	return 0;
}

//*************************************************************************
