/************************************************************************//**
 * File: srptrjdt.cpp
 * Description: Calculation of Electron Trajectory in Periodic Magnetic Field of an ID
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srtrjdat.h"
#include "srptrjdt.h"

//*************************************************************************

//int srTMagFieldPeriodic::SetupFieldBasedArrays(srTEbmDat& EbmDat, int Ns, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtE2)
int srTMagFieldPeriodic::SetupFieldBasedArrays(srTEbmDat& EbmDat, long long Ns, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtE2)
{
	const double pi = 3.14159265358979;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
	double InvGamma = 1./EbmDat.Gamma;
	double Btx0 = Btx0Gam*InvGamma, Btz0 = Btz0Gam*InvGamma, x0 = x0Gam*InvGamma, z0 = z0Gam*InvGamma;

	double *IntBtE2Arr = pIntBtE2;

	double pi_d_lu = pi/PerLength;
	double sStep = (Ns > 1)? PerLength/(Ns - 1) : PerLength;
	double HalfStep = 0.5*sStep;

	double s = 0.;
	double PrevBtE2 = 0.;

	//for(int is=0; is<Ns; is++)
	for(long long is=0; is<Ns; is++)
	{
		*pBtx = Btx0; *pBtz = Btz0; *pX = x0 + Btx0*s; *pZ = z0 + Btz0*s; *pIntBtE2 = 0.;
		for(int i=0; i<AmOfHarm; i++)
		{
			srTMagHarm& Harm = HarmVect[i];

			double a = pi_d_lu*Harm.HarmNo;
			double as = a*s;
			double aspphi = as + Harm.Phase;
			double Sin_as = sin(as), Sin_aspphi = sin(aspphi), Cos_aspphi = cos(aspphi), Sin_phi = sin(Harm.Phase);
			//double K_d_gamma_d_j = Harm.K/(EbmDat.Gamma*Harm.HarmNo);
			double K_d_gamma_d_j = Harm.K/EbmDat.Gamma; //OC221102

			double AddBt = 2.*K_d_gamma_d_j*Sin_as*Cos_aspphi;
			double AddCrd = K_d_gamma_d_j*(Sin_as*Sin_aspphi/a - s*Sin_phi);

			if(Harm.XorZ == 'z')
			{
				*pBtx -= AddBt; *pX -= AddCrd;
			}
			else if(Harm.XorZ == 'x')
			{
				*pBtz += AddBt; *pZ += AddCrd;
			}
		}

		double BtE2 = (*pBtx)*(*pBtx) + (*pBtz)*(*pBtz);
		if(is > 0)
		{
			*pIntBtE2 = *(pIntBtE2 - 1) + HalfStep*(BtE2 + PrevBtE2);
		}
		PrevBtE2 = BtE2;

		s += sStep;
		pBtx++; pBtz++; pX++; pZ++; pIntBtE2++;
	}

	HalfKxE2pKzE2 = IntBtE2Arr[Ns - 1]/(EbmDat.GammaEm2*PerLength);
	return 0;
}

//*************************************************************************

int srTPerTrjDat::SetUpFieldBasedArraysAtOnePeriod(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{
	const int NsPerPeriod = 7;

	FieldBasedArrays.Ns = NsPerPeriod;
	FieldBasedArrays.Nper = int(MagPer.TotLength/MagPer.PerLength);
	FieldBasedArrays.sStart = 0.;
	FieldBasedArrays.sStep = MagPer.PerLength/(NsPerPeriod - 1);

	int result;
	if(result = FieldBasedArrays.AllocateArrays(NsPerPeriod, Keys)) return result;
	CompTotalTrjData(Keys, FieldBasedArrays);
	return 0;
}

//*************************************************************************

int srTPerTrjDat::SetUpFieldBasedArraysTotal(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{
	const int NsPerPeriod = 7;

	int AmOfPer = int(MagPer.TotLength/MagPer.PerLength);
	double ActualTotLength = AmOfPer*MagPer.PerLength;

	FieldBasedArrays.Ns = NsPerPeriod*AmOfPer;
	FieldBasedArrays.Nper = 1;
	FieldBasedArrays.sStart = -(AmOfPer >> 1)*MagPer.PerLength;
	FieldBasedArrays.sStep = ActualTotLength/(FieldBasedArrays.Ns - 1);

	int result;
	if(result = FieldBasedArrays.AllocateArrays(FieldBasedArrays.Ns, Keys)) return result;
	CompTotalTrjData(Keys, FieldBasedArrays);
	return 0;
}

//*************************************************************************

void srTPerTrjDat::CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{// Length units are m here !
	const double pi = 3.14159265358979;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	MagPer.FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
	double InvGamma = 1./EbmDat.Gamma;

	double Btx0 = Btx0Gam*InvGamma, Btz0 = Btz0Gam*InvGamma, x0 = x0Gam*InvGamma, z0 = z0Gam*InvGamma;
	double pi_d_lu = pi/MagPer.PerLength;
	double Bcon = 0.010709839006/MagPer.PerLength;

	double sStep = FieldBasedArrays.sStep;
	double HalfStep = 0.5*sStep;
	double s = FieldBasedArrays.sStart;

	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	double *pX = FieldBasedArrays.XArr, *pZ = FieldBasedArrays.ZArr;
	double *pIntBtxE2 = FieldBasedArrays.IntBtxE2Arr, *pIntBtzE2 = FieldBasedArrays.IntBtzE2Arr;
	double *pdBxds = FieldBasedArrays.dBxdsArr, *pdBzds = FieldBasedArrays.dBzdsArr;

	double PrevBtxE2 = 0., PrevBtzE2 = 0.;
	double Btx, Btz;

	//for(int is=0; is<FieldBasedArrays.Ns; is++)
	for(long long is=0; is<FieldBasedArrays.Ns; is++)
	{
		Btx = Btx0; Btz = Btz0;

		if(Keys.dBxds_) *pdBxds = 0.;
		if(Keys.dBzds_) *pdBzds = 0.;
		if(Keys.Bx_) *pBx = 0.;
		if(Keys.Bz_) *pBz = 0.;
		if(Keys.Btx_) *pBtx = Btx;
		if(Keys.Btz_) *pBtz = Btz;
		if(Keys.X_) *pX = x0 + Btx0*s;
		if(Keys.Z_) *pZ = z0 + Btz0*s;
		if(Keys.IntBtxE2_) *pIntBtxE2 = 0.;
		if(Keys.IntBtzE2_) *pIntBtzE2 = 0.;

		for(int i=0; i<MagPer.AmOfHarm; i++)
		{
			srTMagHarm& Harm = (MagPer.HarmVect)[i];
			double a = pi_d_lu*Harm.HarmNo;
			double as = a*s;
			double aspphi = as + Harm.Phase;
			double Sin_as = sin(as), Sin_aspphi = sin(aspphi), Cos_aspphi = cos(aspphi), Sin_phi = sin(Harm.Phase);
			double K_d_gamma_d_j = Harm.K/(EbmDat.Gamma*Harm.HarmNo);
			double AddBt = 2.*K_d_gamma_d_j*Sin_as*Cos_aspphi;
			double AddCrd = K_d_gamma_d_j*(Sin_as*Sin_aspphi/a - s*Sin_phi);
			double B0 = Harm.K*Bcon;
			double dBds0 = -2.*a*B0;

			if(Harm.XorZ == 'z')
			{
				Btx -= AddBt;

				if(Keys.dBzds_) *pdBzds += dBds0*sin(2.*as + Harm.Phase);
				if(Keys.Bz_) *pBz += B0*cos(2.*as + Harm.Phase);
				if(Keys.Btx_) *pBtx = Btx; 
				if(Keys.X_) *pX -= AddCrd;
			}
			else if(Harm.XorZ == 'x')
			{
				Btz += AddBt;

				if(Keys.dBxds_) *pdBxds += dBds0*sin(2.*as + Harm.Phase);
				if(Keys.Bx_) *pBx += B0*cos(2.*as + Harm.Phase);
				if(Keys.Btz_) *pBtz = Btz;
				if(Keys.Z_) *pZ += AddCrd;
			}
		}

		if(Keys.IntBtxE2_)
		{
			double BtxE2 = Btx*Btx;
			if(is > 0)
			{
				*pIntBtxE2 = *(pIntBtxE2 - 1) + HalfStep*(BtxE2 + PrevBtxE2);
			}
			PrevBtxE2 = BtxE2;
		}
		if(Keys.IntBtzE2_)
		{
			double BtzE2 = Btz*Btz;
			if(is > 0)
			{
				*pIntBtzE2 = *(pIntBtzE2 - 1) + HalfStep*(BtzE2 + PrevBtzE2);
			}
			PrevBtzE2 = BtzE2;
		}

		s += sStep;

		if(Keys.dBxds_) pdBxds++;
		if(Keys.dBzds_) pdBzds++;
		if(Keys.Bx_) pBx++;
		if(Keys.Bz_) pBz++;
		if(Keys.Btx_) pBtx++;
		if(Keys.Btz_) pBtz++;
		if(Keys.X_) pX++;
		if(Keys.Z_) pZ++;
		if(Keys.IntBtxE2_) pIntBtxE2++;
		if(Keys.IntBtzE2_) pIntBtzE2++;
	}
}

//*************************************************************************

//void srTPerTrjDat::CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
void srTPerTrjDat::CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
{// Length units are m here !
	const double pi = 3.14159265358979;
	//const double HalfPi = 0.5*pi;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	MagPer.FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
	double InvGamma = 1./EbmDat.Gamma;

	double Btx0 = Btx0Gam*InvGamma, Btz0 = Btz0Gam*InvGamma, x0 = x0Gam*InvGamma, z0 = z0Gam*InvGamma;
	double pi_d_lu = pi/MagPer.PerLength;
	double Bcon = 0.010709839006/MagPer.PerLength;

	double sStep = (sEn - sSt)/(Np - 1);
	//double HalfStep = 0.5*sStep;
	double s = sSt;

	//for(long is=0; is<Np; is++)
	for(long long is=0; is<Np; is++)
	{
		*pBx = 0.; *pBz = 0.;
		*pBtx = Btx0; *pBtz = Btz0;
		*pX = x0 + Btx0*s; *pZ = z0 + Btz0*s;

		for(int i=0; i<MagPer.AmOfHarm; i++)
		{
			srTMagHarm& Harm = (MagPer.HarmVect)[i];
			double HarmPhase = Harm.Phase;

			double a = pi_d_lu*Harm.HarmNo;
			double as = a*s;
			double aspphi = as + HarmPhase;
			double Sin_as = sin(as), Sin_aspphi = sin(aspphi), Cos_aspphi = cos(aspphi), Sin_phi = sin(HarmPhase);
			double K_d_gamma_d_j = Harm.K/(EbmDat.Gamma*Harm.HarmNo);
			double AddBt = 2.*K_d_gamma_d_j*Sin_as*Cos_aspphi;
			double AddCrd = K_d_gamma_d_j*(Sin_as*Sin_aspphi/a - s*Sin_phi);
			double B0 = Harm.K*Bcon;

			if(Harm.XorZ == 'z')
			{
				*pBz += B0*cos(2.*as + HarmPhase);
				*pBtx -= AddBt; 
				*pX -= AddCrd;
			}
			else if(Harm.XorZ == 'x')
			{
				*pBx += B0*cos(2.*as + HarmPhase);
				*pBtz += AddBt;
				*pZ += AddCrd;
			}
		}

		s += sStep;
		pBx++; pBz++; pBtx++; pBtz++; pX++; pZ++;
	}
}

//*************************************************************************

//void srTPerTrjDat::CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
void srTPerTrjDat::CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
{// This does not compute correctly pIntBtxE2, pIntBtzE2.

	const double pi = 3.14159265358979;
	//const double HalfPi = 0.5*pi;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	MagPer.FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
	double InvGamma = 1./EbmDat.Gamma;

	double Btx0 = Btx0Gam*InvGamma, Btz0 = Btz0Gam*InvGamma, x0 = x0Gam*InvGamma, z0 = z0Gam*InvGamma;
	double pi_d_lu = pi/MagPer.PerLength;
	double Bcon = 0.010709839006/MagPer.PerLength;

	double sStep = (sEn - sSt)/(Np - 1);
	double HalfStep = 0.5*sStep;
	double s = sSt;
	double PrevBtxE2 = 0., PrevBtzE2 = 0.;

	//for(long is=0; is<Np; is++)
	for(long long is=0; is<Np; is++)
	{
		*pdBxds = 0.; *pdBzds = 0.;
		*pBx = 0.; *pBz = 0.;
		*pBtx = Btx0; *pBtz = Btz0;
		*pX = x0 + Btx0*s; *pZ = z0 + Btz0*s;
		*pIntBtxE2 = 0.; *pIntBtzE2 = 0.;

		for(int i=0; i<MagPer.AmOfHarm; i++)
		{
			srTMagHarm& Harm = (MagPer.HarmVect)[i];
			double HarmPhase = Harm.Phase;

			double a = pi_d_lu*Harm.HarmNo;
			double as = a*s;
			double aspphi = as + HarmPhase;
			double Sin_as = sin(as), Sin_aspphi = sin(aspphi), Cos_aspphi = cos(aspphi), Sin_phi = sin(HarmPhase);
			double K_d_gamma_d_j = Harm.K/(EbmDat.Gamma*Harm.HarmNo);
			double AddBt = 2.*K_d_gamma_d_j*Sin_as*Cos_aspphi;
			double AddCrd = K_d_gamma_d_j*(Sin_as*Sin_aspphi/a - s*Sin_phi);
			double B0 = Harm.K*Bcon;
			double dBds0 = -2.*a*B0;

			if(Harm.XorZ == 'z')
			{
				*pdBzds += dBds0*sin(2.*as + Harm.Phase);
				*pBz += B0*cos(2.*as + HarmPhase);
				*pBtx -= AddBt; 
				*pX -= AddCrd;
			}
			else if(Harm.XorZ == 'x')
			{
				*pdBxds += dBds0*sin(2.*as + Harm.Phase);
				*pBx += B0*cos(2.*as + HarmPhase);
				*pBtz += AddBt;
				*pZ += AddCrd;
			}
		}

		double BtxE2 = (*pBtx)*(*pBtx);
		double BtzE2 = (*pBtz)*(*pBtz);
		if(is > 0)
		{
			*pIntBtxE2 = *(pIntBtxE2 - 1) + HalfStep*(BtxE2 + PrevBtxE2);
			*pIntBtzE2 = *(pIntBtzE2 - 1) + HalfStep*(BtzE2 + PrevBtzE2);
		}
		PrevBtxE2 = BtxE2;
		PrevBtzE2 = BtzE2;
		
		s += sStep;
		pdBxds++; pdBzds++; pBx++; pBz++; pBtx++; pBtz++; pX++; pZ++; pIntBtxE2++; pIntBtzE2++;
	}
}

//*************************************************************************

void srTPerTrjDat::CompTrjDataDerivedAtPointPowDens(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz)
{// Length units are m here !
	const double pi = 3.14159265358979;
	//const double HalfPi = 0.5*pi;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	MagPer.FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
	double InvGamma = 1./EbmDat.Gamma;

	double Btx0 = Btx0Gam*InvGamma, Btz0 = Btz0Gam*InvGamma, x0 = x0Gam*InvGamma, z0 = z0Gam*InvGamma;
	double pi_d_lu = pi/MagPer.PerLength;
	double Bcon = 0.010709839006/MagPer.PerLength;

	Bx = 0.; Bz = 0.;
	Btx = Btx0; Btz = Btz0;
	X = x0 + Btx0*s; Z = z0 + Btz0*s;

	for(int i=0; i<MagPer.AmOfHarm; i++)
	{
		srTMagHarm& Harm = (MagPer.HarmVect)[i];
		double HarmPhase = Harm.Phase;

		double a = pi_d_lu*Harm.HarmNo;
		double as = a*s;
		double aspphi = as + HarmPhase;
		double Sin_as = sin(as), Sin_aspphi = sin(aspphi), Cos_aspphi = cos(aspphi), Sin_phi = sin(HarmPhase);
		double K_d_gamma_d_j = Harm.K/(EbmDat.Gamma*Harm.HarmNo);
		double AddBt = 2.*K_d_gamma_d_j*Sin_as*Cos_aspphi;
		double AddCrd = K_d_gamma_d_j*(Sin_as*Sin_aspphi/a - s*Sin_phi);
		double B0 = Harm.K*Bcon;

		if(Harm.XorZ == 'z')
		{
			Bz += B0*cos(2.*as + HarmPhase);
			Btx -= AddBt; 
			X -= AddCrd;
		}
		else if(Harm.XorZ == 'x')
		{
			Bx += B0*cos(2.*as + HarmPhase);
			Btz += AddBt;
			Z += AddCrd;
		}
	}
}

//*************************************************************************

int srTPerTrjDat::ConvertToArbTrjDat(char Periodicity, srTWfrSmp&, srTGenTrjHndl& TrjHndl)
{// This performs only the actions that are normally performed at setting up srTTrjDat via srTSend !!!
 // if(Periodicity == 1) - Arb. field
 // if(Periodicity == 2) - Periodic field
	const int AmOfPoPerPeriod = 400; // To steer
	int result;

	srTTrjDat* pNewTrjDat = new srTTrjDat();
	if(pNewTrjDat == 0) return MEMORY_ALLOCATION_FAILURE;

//El. Beam
	pNewTrjDat->EbmDat = EbmDat;

	double Btx0Gam, Btz0Gam, x0Gam, z0Gam;
	MagPer.FindInitialTrjAngles(Btx0Gam, Btz0Gam, x0Gam, z0Gam);
//--
	double InvGam = 1./EbmDat.Gamma;
	pNewTrjDat->EbmDat.x0 = -x0Gam*InvGam;
	pNewTrjDat->EbmDat.dxds0 = -Btx0Gam*InvGam;
	pNewTrjDat->EbmDat.z0 = -z0Gam*InvGam;
	pNewTrjDat->EbmDat.dzds0 = -Btz0Gam*InvGam; //Not sure ...

//--
//Mag. Field
	CheckIfHorOrVertFieldIsZero();
	pNewTrjDat->HorFieldIsNotZero = HorFieldIsNotZero;
	pNewTrjDat->VerFieldIsNotZero = VerFieldIsNotZero;

	int Nper = int(MagPer.TotLength/MagPer.PerLength);
	//long TotAmOfPo = (Periodicity == 2)? ((AmOfPoPerPeriod << 1) + 1) : AmOfPoPerPeriod*Nper;
	long long TotAmOfPo = (Periodicity == 2)? ((AmOfPoPerPeriod << 1) + 1) : AmOfPoPerPeriod*Nper;

	pNewTrjDat->LenFieldData = TotAmOfPo;
	pNewTrjDat->AuxBxLen = TotAmOfPo;
	pNewTrjDat->AuxBzLen = TotAmOfPo;

	pNewTrjDat->BxInData = new srTFunDer[TotAmOfPo];
	if(pNewTrjDat->BxInData == 0) return MEMORY_ALLOCATION_FAILURE;

	pNewTrjDat->BzInData = new srTFunDer[TotAmOfPo];
	if(pNewTrjDat->BzInData == 0) return MEMORY_ALLOCATION_FAILURE;

	pNewTrjDat->NperTot = (Periodicity == 2)? Nper : 1;
	pNewTrjDat->sStep = (Periodicity == 2)? 2.*MagPer.PerLength/(TotAmOfPo - 1) : MagPer.PerLength*Nper/(TotAmOfPo - 1);
	pNewTrjDat->Inv_Step = 1./pNewTrjDat->sStep;

	if(((Nper >> 1) << 1) != Nper) // Odd
	{
		if(Periodicity == 2)
		{
			pNewTrjDat->NperLeft = (Nper - 1) >> 1;
			pNewTrjDat->sStart = -0.5*MagPer.PerLength;
		}
		else
		{
			pNewTrjDat->NperLeft = 0;
			int NperLeftRight = (Nper - 1) >> 1;
			pNewTrjDat->sStart = -(0.5 + NperLeftRight)*MagPer.PerLength;
		}
	}
	else // Even
	{
		if(Periodicity == 2)
		{
			pNewTrjDat->NperLeft = Nper >> 1;
			pNewTrjDat->sStart = 0.;
		}
		else
		{
			pNewTrjDat->NperLeft = 0;
			int NperLeft = Nper >> 1;
			pNewTrjDat->sStart = -NperLeft*MagPer.PerLength;
		}
	}
	pNewTrjDat->EbmDat.s0 = pNewTrjDat->sStart;

	srTFieldBasedArrayKeys Keys;
	Keys.Bx_ = Keys.dBxds_ = Keys.Bz_ = Keys.dBzds_ = 1;
	srTFieldBasedArrays AuxFieldArrays;
	if(result = AuxFieldArrays.AllocateArrays(TotAmOfPo, Keys)) return result;

	AuxFieldArrays.sStart = pNewTrjDat->sStart;
	AuxFieldArrays.sStep = pNewTrjDat->sStep;
	AuxFieldArrays.Ns = TotAmOfPo;
	AuxFieldArrays.Nper = 1;

	CompTotalTrjData(Keys, AuxFieldArrays);

	//double FieldZeroTol = pNewTrjDat->FieldZeroTolerance;
	srTFunDer *tFunDerX = pNewTrjDat->BxInData;
	srTFunDer *tFunDerZ = pNewTrjDat->BzInData;
	double *tBx = AuxFieldArrays.BxArr, *tdBxds = AuxFieldArrays.dBxdsArr;
	double *tBz = AuxFieldArrays.BzArr, *tdBzds = AuxFieldArrays.dBzdsArr;
	//for(long i=0; i<TotAmOfPo; i++)
	for(long long i=0; i<TotAmOfPo; i++)
	{
		tFunDerX->f = *(tBx++);
		(tFunDerX++)->dfds = *(tdBxds++);
		tFunDerZ->f = *(tBz++);
		(tFunDerZ++)->dfds = *(tdBzds++);
	}

	srTGenTrjHndl LocTrjHndl(pNewTrjDat);
	TrjHndl = LocTrjHndl;
	return 0;
}

//*************************************************************************

//int srTPerTrjDat::SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long& Ns, int& NperTot, int& NperLeft) 
int srTPerTrjDat::SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long long& Ns, int& NperTot, int& NperLeft) 
{
	const double MinStepSteerCoef = 0.4;
	const int NsMin = 10;
	
	double AbsMax = MaxAbsField();
	double Rmin = 3.3*EbmDat.Energy/AbsMax;
	double dsMin = MinStepSteerCoef*Rmin/EbmDat.Gamma;
	
	NperTot = int(MagPer.TotLength/MagPer.PerLength);

	if(LongIntType == 1) // Total range
	{
		double ActualTotLength = NperTot*MagPer.PerLength;

		sStart = -0.5*ActualTotLength;
		//Ns = int(ActualTotLength/dsMin);
		Ns = (long long)(ActualTotLength/dsMin);
		if(Ns < NsMin) Ns = NsMin;

		sStep = ActualTotLength/double(Ns - 1);

		NperTot = 1;
		NperLeft = 0;
	}
	else if(LongIntType == 2) // One period
	{
		char NperTotIsEven = (((NperTot >> 1) << 1) == NperTot);
		if(NperTotIsEven)
		{
			sStart = 0.;
			NperLeft = NperTot >> 1;
		}
		else
		{
			sStart = -0.5*MagPer.PerLength;
			NperLeft = (NperTot - 1) >> 1;
		}
		//Ns = int(MagPer.PerLength/dsMin);
		Ns = (long long)(MagPer.PerLength/dsMin);
		if(Ns < NsMin) Ns = NsMin;

		sStep = MagPer.PerLength/double(Ns - 1);
	}
	return 0;
}

//*************************************************************************

void srTPerTrjDat::AnalizeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ)
{
	FieldIsSymOverX = FieldIsSymOverZ = 0;
	MagPer.AnalyzeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
}

//*************************************************************************
