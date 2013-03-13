/************************************************************************//**
 * File: srebmdat.h
 * Description: Electron beam class (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SREBMDAT_H
#define __SREBMDAT_H

#include <math.h>
#include "srobject.h"

//*************************************************************************

class srTEbmDat : public CGenObject {
public:
// Use normal units in this structure: GeV, A, m, r !!!

	double Energy, Current;
	double Neb, CurrentPeak; //number of electrons in bunch and peak current - two parameters which are not independent
	double s0, x0, dxds0, z0, dzds0, sc;
	double Mxx, Mxxp, Mxpxp;
	double Mzz, Mzzp, Mzpzp;
	double Mxz, Mxpz, Mxzp, Mxpzp;
	double Mee; //(SigmaE/E)^2
	double Mss;
	double Mse;
	double Mxe, Mxpe, Mze, Mzpe;
	double Mxs, Mxps, Mzs, Mzps;

	double SigmaRelE;
	double Gamma, GammaEm2;
	int nQ; //SRWL: number and sign of charges in units of electron charge; "-1" for electron, "1" for proton

	int TypeDistrTransverse; // 1- uniform; 2- Gaussian; 3- parabolic
	int TypeDistrLongitudinal; // 1- infinite uniform; 2- gaussian
	double ShotNoiseFactor;

	float *pElecDistr;
	//waveHndl wElecDistr;
	char** wElecDistr;
	int hStateElecDistr;
	long nTotMacroPart; //(number of macro-particles) x (number of slices vs time)

	char InputWasModified;

	srTEbmDat()
	{
		PresetAll();
	}
	srTEbmDat(double I, double InNeb, double* pMom1, int nMom1, double* pMom2, int nMom2, double s, int nq =-1)
	{
		PresetAll();

		Current = I;
		Neb = InNeb;
		Energy = x0 = dxds0 = z0 = dzds0 = sc = 0.;
		if(pMom1 != 0) 
		{
			if(nMom1 > 5) { Energy = pMom1[0]; x0 = pMom1[1]; dxds0 = pMom1[2]; z0 = pMom1[3]; dzds0 = pMom1[4]; sc = pMom1[5];}
			else if(nMom1 > 4) { Energy = pMom1[0]; x0 = pMom1[1]; dxds0 = pMom1[2]; z0 = pMom1[3]; dzds0 = pMom1[4];}
			else if(nMom1 > 3) { Energy = pMom1[0]; x0 = pMom1[1]; dxds0 = pMom1[2]; z0 = pMom1[3];}
			else if(nMom1 > 2) { Energy = pMom1[0]; x0 = pMom1[1]; dxds0 = pMom1[2];}
			else if(nMom1 > 1) { Energy = pMom1[0]; x0 = pMom1[1];}
			else if(nMom1 > 0) { Energy = pMom1[0];}
		}
		Mxx = Mxxp = Mxpxp = Mzz = Mzzp = Mzpzp = Mxz = Mxpz = Mxzp = Mxpzp = Mee = Mss = 0.;
		Mse = Mxe = Mxpe = Mze = Mzpe = Mxs = Mxps = Mzs = Mzps = 0.;
		if(pMom2 != 0) 
		{
			if(nMom2 > 20) 
			{ 
				Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10]; Mss = pMom2[11];
				Mse = pMom2[12]; Mxe = pMom2[13]; Mxpe = pMom2[14]; Mze = pMom2[15]; Mzpe = pMom2[16]; Mxs = pMom2[17]; Mxps = pMom2[18]; Mzs = pMom2[19]; Mzps = pMom2[20];
			}
			else if(nMom2 > 11) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10]; Mss = pMom2[11];}
			else if(nMom2 > 10) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10];}
			else if(nMom2 > 9) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9];}
			else if(nMom2 > 8) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8];}
			else if(nMom2 > 7) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7];}
			else if(nMom2 > 6) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6];}
			else if(nMom2 > 5) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5];}
			else if(nMom2 > 4) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4];}
			else if(nMom2 > 3) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3];}
			else if(nMom2 > 2) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2];}
			else if(nMom2 > 1) { Mxx = pMom2[0]; Mxxp = pMom2[1];}
			else if(nMom2 > 0) { Mxx = pMom2[0];}
		}
		s0 = s; 
		SetupGamma(Energy);
		SigmaRelE = 0.;
		if(Mee > 0.) SigmaRelE = sqrt(Mee);

		nQ = nq;

		pElecDistr = 0;
		nTotMacroPart = 0;
	}

	void PresetAll()
	{
		InputWasModified = 1;
		Energy = Current = Gamma = GammaEm2 = SigmaRelE = 0.;
		Neb = CurrentPeak = 0;

        s0 = sc = x0 = dxds0 = z0 = dzds0 = 0.;
		Mxx = Mxxp = Mxpxp = 0.;
		Mzz = Mzzp = Mzpzp = 0.;
		Mxz = Mxpz = Mxzp = Mxpzp = 0.;
		Mee = 0; //(SigmaE/E)^2
		Mse = Mxe = Mxpe = Mze = Mzpe = 0.;
		Mxs = Mxps = Mzs = Mzps = 0.;
		Mss = 1.E+23;

		nQ = -1; //SRWL
		
		TypeDistrTransverse = 2;
		TypeDistrLongitudinal = 1;
		ShotNoiseFactor = 1;

		pElecDistr = 0;
		wElecDistr = 0;
		hStateElecDistr = 0;
		nTotMacroPart = 0;
	}

	void SetupGamma(double ElectronEnergyInGeV)
	{
		Energy = ElectronEnergyInGeV;
		Gamma = ElectronEnergyInGeV*1000./0.511003414;
		if(Gamma != 0.) GammaEm2 = 1./Gamma/Gamma;
	}

	void SetNewEnergy(double NewEnergy_GeV)
	{
		Energy = NewEnergy_GeV;
		SetupGamma(NewEnergy_GeV);
	}

	void SetCurrentAndMom2(double inI, double* pMom2, int nMom2)
	{
		Current = inI;

		Mxx = Mxxp = Mxpxp = Mzz = Mzzp = Mzpzp = Mxz = Mxpz = Mxzp = Mxpzp = Mee = Mss = 0.;
		Mse = Mxe = Mxpe = Mze = Mzpe = Mxs = Mxps = Mzs = Mzps = 0.;
		if(pMom2 != 0)
		{
			if(nMom2 > 20)
			{
				Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10]; Mss = pMom2[11];
				Mse = pMom2[12]; Mxe = pMom2[13]; Mxpe = pMom2[14]; Mze = pMom2[15]; Mzpe = pMom2[16]; Mxs = pMom2[17]; Mxps = pMom2[18]; Mzs = pMom2[19]; Mzps = pMom2[20];
			}
			else if(nMom2 > 11) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10]; Mss = pMom2[11];}
			else if(nMom2 > 10) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9]; Mee = pMom2[10];}
			else if(nMom2 > 9) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8]; Mxpzp = pMom2[9];}
			else if(nMom2 > 8) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7]; Mxzp = pMom2[8];}
			else if(nMom2 > 7) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6]; Mxpz = pMom2[7];}
			else if(nMom2 > 6) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5]; Mxz = pMom2[6];}
			else if(nMom2 > 5) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4]; Mzpzp = pMom2[5];}
			else if(nMom2 > 4) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3]; Mzzp = pMom2[4];}
			else if(nMom2 > 3) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2]; Mzz = pMom2[3];}
			else if(nMom2 > 2) { Mxx = pMom2[0]; Mxxp = pMom2[1]; Mxpxp = pMom2[2];}
			else if(nMom2 > 1) { Mxx = pMom2[0]; Mxxp = pMom2[1];}
			else if(nMom2 > 0) { Mxx = pMom2[0];}
		}
	}

	void GetMom1(double* pMom1)
	{
		if(pMom1 == 0) return;
        pMom1[0] = Energy; 
		pMom1[1] = x0; 
		pMom1[2] = dxds0; 
		pMom1[3] = z0; 
		pMom1[4] = dzds0;
	}
	void GetMom2(double* pMom2)
	{
		if(pMom2 == 0) return;
        pMom2[0] = Mxx; 
		pMom2[1] = Mxxp; 
		pMom2[2] = Mxpxp; 
		pMom2[3] = Mzz; 
		pMom2[4] = Mzzp; 
		pMom2[5] = Mzpzp; 
		pMom2[6] = Mxz; 
		pMom2[7] = Mxpz; 
		pMom2[8] = Mxzp; 
		pMom2[9] = Mxpzp; 
		pMom2[10] = Mee; 
		pMom2[11] = Mss;
	}
	double GetSigmaE_GeV()
	{
		SigmaRelE = ::sqrt(::fabs(Mee));
		return SigmaRelE*Energy;
	}
	double GetSigma_x() { return ::sqrt(::fabs(Mxx));}
	double GetSigma_dxds() { return ::sqrt(::fabs(Mxpxp));}
	double GetSigma_z() { return ::sqrt(::fabs(Mzz));}
	double GetSigma_dzds() { return ::sqrt(::fabs(Mzpzp));}
	double GetSigma_s() { return ::sqrt(::fabs(Mss));}

	double EmittanceX()
	{
		return sqrt(::fabs(Mxx*Mxpxp - Mxxp*Mxxp));
	}
	double EmittanceZ()
	{
		return sqrt(::fabs(Mzz*Mzpzp - Mzzp*Mzzp));
	}
	double NormalizedEmittanceX()
	{
		return EmittanceX()*Gamma; // Check this !!!
		//return EmittanceX()*Gamma/3.1415926; // Check this !!!
	}
	double NormalizedEmittanceZ()
	{
		return EmittanceZ()*Gamma; // Check this !!!
		//return EmittanceZ()*Gamma/3.1415926; // Check this !!!
	}
	double AlphaX()
	{
		return -Mxxp/EmittanceX();
	}
	double AlphaZ()
	{
		return -Mzzp/EmittanceZ();
	}

	double MultForPartDistrFunc(char x_or_z)
	{
		if(x_or_z == 'x') return 0.5/(Mxx*Mxpxp - Mxxp*Mxxp);
		else return 0.5/(Mzz*Mzpzp - Mzzp*Mzzp);
	}
	double BetaPartDistr(char x_or_z) { return ((x_or_z == 'x')? Mxx : Mzz)*MultForPartDistrFunc(x_or_z);}
	double GammaPartDistr(char x_or_z) { return ((x_or_z == 'x')? Mxpxp : Mzpzp)*MultForPartDistrFunc(x_or_z);}
	double AlphaPartDistr(char x_or_z) { return ((x_or_z == 'x')? -Mxxp : -Mzzp)*MultForPartDistrFunc(x_or_z);}
	double ThetaPartDistr(char x_or_xp, char z_or_zp) { return 0;} //xz cross-terms: program when necessary
	//double BgamPartDistr() { return (Mee > 0)? sqrt(Mee) : 0;}
	double BgamPartDistr() { return (Mee > 0)? 0.5/Mee : 1e+50;}
	double PggPartDistr_NoCross() { return (Mee > 0)? 0.5/Mee : 1e+50;}
	double PssPartDistr_NoCross() { return (Mss > 0)? 0.5/Mss : 1e+50;}

	double Inv_B_rho() { return 0.299792458*nQ/Energy;} //1/(B*rho) in [1/(T*m)]

	void analyseDistribSymmetry(char& ElecBeamIsSymOverX, char& ElecBeamIsSymOverZ)
	{
		ElecBeamIsSymOverX = ElecBeamIsSymOverZ = 1;
		if((Mxz != 0) || (Mxpz != 0) || (Mxzp != 0) || (Mxpzp != 0))
		{
            ElecBeamIsSymOverX = ElecBeamIsSymOverZ = 0;
		}
		if((Mxe != 0) || (Mxpe != 0) || (Mxs != 0) || (Mxps != 0)) ElecBeamIsSymOverX = 0;
		if((Mze != 0) || (Mzpe != 0) || (Mzs != 0) || (Mzps != 0)) ElecBeamIsSymOverZ = 0;
	}

	static int Mom2EmitAndTwiss(double* pEmit, double* pBeta, double* pAlpha, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double Eta, double EtaPrime)
	{
		double MixMom2 = MixMom*MixMom, Eta2 = Eta*Eta, EtaPrime2 = EtaPrime*EtaPrime;
		double Emit2 = Sigma2*SigmaPrime2 - MixMom2 - SigmaE2*(EtaPrime2*Sigma2 + Eta2*SigmaPrime2 - 2*MixMom*Eta*EtaPrime);
		if(Emit2 <= 0.) 
		{
			*pEmit = 0.; *pBeta = 0.; *pAlpha = 0.;
			return 0;
		}

		*pEmit = sqrt(Emit2);
		*pBeta = (Sigma2 - SigmaE2*Eta2)/(*pEmit);
		*pAlpha = (SigmaE2*Eta*EtaPrime - MixMom)/(*pEmit);
		return 0;
	}

	static int MomAndEmit2Twiss(double* pBeta, double* pAlpha, double* pEta, double* pEmit, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double EtaPrime)
	{
		double MixMom2 = MixMom*MixMom, EtaPrime2 = EtaPrime*EtaPrime;
		//double Emit2 = Sigma2*SigmaPrime2 - MixMom2 - SigmaE2*(EtaPrime2*Sigma2 + Eta2*SigmaPrime2 - 2*MixMom*Eta*EtaPrime);
		
		if((Sigma2 == 0) || (SigmaPrime2 == 0))
		{
			*pBeta = 0.; *pAlpha = 0.; *pEta = 0.; *pEmit = 0.;
			return 0;
		}
		else if(SigmaE2 <= 0.) //Modify *pEmit
		{
			double Emit2 = Sigma2*SigmaPrime2 - MixMom2;
			if(Emit2 <= 0)
			{
				*pBeta = 0.; *pAlpha = 0.; *pEta = 0.; *pEmit = 0.;
				return 0;
			}
			*pEmit = sqrt(Emit2);
		}
		else //Assume *pEmit is defined at input
		{
			double C = Sigma2*EtaPrime2 + ((*pEmit)*(*pEmit) - Sigma2*SigmaPrime2 + MixMom2)/SigmaE2;
			double argRoot = MixMom2*EtaPrime2 - SigmaPrime2*C;
			if(argRoot < 0.)
			{
				*pEta = 0.; //?
			}
			else *pEta = (MixMom*EtaPrime + sqrt(argRoot))/SigmaPrime2;
		}

		*pBeta = (Sigma2 - SigmaE2*(*pEta)*(*pEta))/(*pEmit);
		*pAlpha = (SigmaE2*(*pEta)*EtaPrime - MixMom)/(*pEmit);
		return 0;
	}

};

//*************************************************************************

//typedef CHandle<srTEbmDat> CHElBeam;
typedef CSmartPtr<srTEbmDat> CHElBeam;

//*************************************************************************

#endif
