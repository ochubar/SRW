/************************************************************************//**
 * File: srtrjdat.h
 * Description: Electron Trajectory (and relevant characteristics) calculation for different types of Magnetic Fields (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRTRJDAT_H
#define __SRTRJDAT_H

//#ifndef __SRSEND_H
//#include "srsend.h"
//#endif

//#ifdef __IGOR_PRO__
//#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
//#else
//#ifndef __SRIGORRE_H
//#include "srigorre.h"
//#endif
//#endif

#include "srstraux.h"
#include "srgtrjdt.h"
#include "srsysuti.h"
#include "srwlib.h"

//*************************************************************************

class srTMagFldTrUnif;
//struct SRWLStructParticleTrajectory;
//typedef struct SRWLStructParticleTrajectory SRWLPrtTrj;

//*************************************************************************

struct srTFunDer {
	double f, dfds;
	srTFunDer(double f_In =0, double dfds_In =0) { f=f_In; dfds=dfds_In;}
};

//*************************************************************************

class srTTrjDat : public srTGenTrjDat {

	double xCorr, BtxCorr, zCorr, BtzCorr, IntBtxE2Corr, IntBtzE2Corr;
	double BtxCorrForX, BtzCorrForZ, BtxCorrForXe2, BtzCorrForZe2;
	//srTSend Send;

	double* AllCf;
	double** BxPlnCf;
	double** BzPlnCf;
	double** BtxPlnCf;
	double** BtzPlnCf;
	double** xPlnCf;
	double** zPlnCf;
	double** IntBtx2PlnCf;
	double** IntBtz2PlnCf;

	DOUBLE *IntBtxE2Arr, *IntBtzE2Arr;

	//int m_estimMinNpForRadInteg; //OC220112
	long long m_estimMinNpForRadInteg; //OC220112

public:

	//int AuxBxLen, AuxBzLen;
	long long AuxBxLen, AuxBzLen;

	srTFunDer* BxInData;
	srTFunDer* BzInData;
	//int LenFieldData;
	long long LenFieldData;
	double sStart, sStep, Inv_Step;

	srTWaveAccessDataD1D xTrjInData, zTrjInData; // for computing SR from Trajectory
	char CompFromTrj;

	int NperTot, NperLeft;

	double FieldZeroTolerance;
	double MaxMemAvail;
	char InputWasModified, LastCompWasNotOK;
	bool m_doNotDeleteData;

	//srTTrjDat(srTEbmDat* pEbmDat, srTMagElem* pMagElem);
	srTTrjDat(srTEbmDat* pEbmDat, srTMagFldTrUnif* pMagElem);
	srTTrjDat(SRWLPrtTrj* pTrj);
	//srTTrjDat(srTTrjDat& anotherTrjDat);
	srTTrjDat()
	{
		Initialize();
	}
	~srTTrjDat()
	{
		if(!m_doNotDeleteData)
		{
			DeleteInitialFieldData();
			DeallocateMemoryForCfs();
			DeallocateQuadPhaseTermsArr();
		}
	}

	void Initialize()
	{
		BxInData = 0; BzInData = 0; 
		BxPlnCf = BzPlnCf = 0;
		BtxPlnCf = BtzPlnCf = 0;
		xPlnCf = zPlnCf = 0;
		IntBtx2PlnCf = IntBtz2PlnCf = 0;
		AllCf = 0;

		FieldZeroTolerance = 1.E-10; // To steer. Do not change in the code.

		InputWasModified = 1;
		AuxBxLen = AuxBzLen = 0;

		MaxMemAvail = 8.;
		LastCompWasNotOK = 1;

		NperTot = 1;
		NperLeft = 0;

		IntBtxE2Arr = IntBtzE2Arr = 0;
		CompFromTrj = 0;

        xCorr = BtxCorr = zCorr = BtzCorr = IntBtxE2Corr = IntBtzE2Corr = 0;
        BtxCorrForX = BtzCorrForZ = BtxCorrForXe2 = BtzCorrForZe2 = 0;

		m_doNotDeleteData = false;
		m_estimMinNpForRadInteg = -1;
	}

	int DeleteInitialFieldData() 
	{
		if(BxInData != 0) delete[] BxInData;
		BxInData = 0;
		if(BzInData != 0) delete[] BzInData;
		BzInData = 0;
		return 0;
	}

	int CompDerivForFieldData(srTFunDer*);
	inline void SetupFldPlnCf(srTFunDer*, double**);

	int AllocateMemoryForCfs();
	int AllocateMemoryForCfs_FromTrj();
	//int AllocateMemoryForCfsFromTrj(int np);
	int AllocateMemoryForCfsFromTrj(long long np);

	int DeallocateMemoryForCfs();

	inline double Derivative(double*, double, int, int =5);
	inline void CubPln(double, double, double, double, double*);
	inline double Field(char, double);

	inline void CompTrjDataDerivedAtPointPowDens(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz);
	inline void CompTrjDataDerivedAtPointPowDens_FromTrj(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz);
	inline void CompTrjDataAndFieldWithDerAtPoint(char, double, double&, double&, double&, double&, double&);
	inline void CompTrjDataAndFieldWithDerAtPoint_FromTrj(char x_or_z, double s, double& dBds, double& B, double& Bt, double& Crd, double& IntBtE2);
	inline void CompTrjDataAndFieldWithDerAtPoint_FromTrjInitial(char x_or_z, double s, double& dBds, double& B, double& Bt, double& Crd, double& IntBtE2);
	inline void CompTrjDataDerivedAtPoint(char, double, double&, double&, double&);
	inline void CompTrjDataDerivedAtPoint_FromTrj(char x_or_z, double s, double& Bt, double& Crd, double& IntBtE2);
	inline void CompTrjDataDerivedAtPoint(double sArg, double& Btx, double& Crdx, double& IntBtE2x, double& Btz, double& Crdz, double& IntBtE2z);
	inline void CompTrjDataDerivedAtPoint_FromTrj(double s, double& Btx, double& Crdx, double& IntBtxE2, double& Btz, double& Crdz, double& IntBtzE2);
	//void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz);
	//void CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz);
	//void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);
	//void CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);
	//void CompTotalTrjDataTrjDisp(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn);
	//void CompTotalTrjDataTrjDisp_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn);
	//void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);
	//void CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);
	void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz);
	void CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz);
	void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);
	void CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds);
	void CompTotalTrjDataTrjDisp(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn);
	void CompTotalTrjDataTrjDisp_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn);
	void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);
	void CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz);

	void CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays);
	void CompTotalTrjData_FromTrj(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays);

	inline double Pol04(double, double*);
	inline double Pol05(double, double*);
	inline double Pol09(double, double*);

	//inline void ComponentArray(int, char, double, double, int, double*);
	inline void ComponentArray(int, char, double, double, long long, double*);

	void SetupIntegrPlnCfs(char);
	inline int ComputeInterpolatingStructure();
	inline int ComputeInterpolatingStructure_FromTrj();
	int ComputeInterpolatingStructureFromTrj(SRWLPrtTrj* pTrj);
	int ComputeInterpolatingStructure_FromTrj1D(char x_or_z);
	int ComputeInterpolatingStructureFromTrj1D(char x_or_z, const SRWLPrtTrj& trj);
	inline void CoefsPol5thOrder(double h, DOUBLE* pf0, double* pCoef);

	inline void CompCorrectionsForTrjDataDerived();

	inline void MakeZeroInputFieldDer();
	inline int EvaluateMemAvailBeforeTrjComp();

	void RecomputeHorOrVertField(int AmOfBxValues, int AmOfBzValues) 
	{ /* To implement: make the initial arrays equal */ }

	int MagFieldPeriodicity() 
	{
		if(NperTot == 1) return 1; // Arbitrary
		else return 2; // Periodic
	}
	char MagFieldIsConstant() { return 0;}

	int SetUpFieldBasedArraysAtOnePeriod(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays) 
	{ 
		return SetUpFieldBasedArraysTotal(Keys, FieldBasedArrays);
	}
	inline int SetUpFieldBasedArraysTotal(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays);
	void AnalizeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ);
	inline int ShowLimitsAndInitInteg(srTWfrSmp& DistrInfoDat, char LongIntType, double& sIntegStart, double& sIntegFin, int& AmOfPer, bool doInit = true);
	int InitTrjComp() { return ComputeInterpolatingStructure();}
	//int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long& Ns, int& OutNperTot, int& OutNperLeft);
	int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long long& Ns, long long& OutNperTot, long long& OutNperLeft);
	void CountFieldExtrem(int& AmOfExtrem, double& AbsMax);
	void FindFieldLimitsBasedOnTolValue(double AbsTolForField, double& sStartLoc, double& sFinLoc);
	void CheckFromTrjIfFieldCompAreZero(SRWLPrtTrj& trj, short& horFieldIsNotZero, short& verFieldIsNotZero);

	int ConvertToArbTrjDat(char, srTWfrSmp&, srTGenTrjHndl&) { return 0;}
	void ShowFullLimits(double& sIntegStart, double& sIntegFin) 
	{
		double sRange = sStep*(LenFieldData - 1);
		sIntegStart = sStart;
		if(NperLeft > 0) sIntegStart -= NperLeft*sRange;
		sIntegFin = sStart + sRange;
		int NperRight = NperTot - NperLeft - 1;
		if(NperRight > 0) sIntegFin += NperRight*sRange;
	}

	int CheckAndSetupTrajectoryLimits();
	int SetupSourcePointFromTrajectory();
	int ComputeOneQuadPhaseTermFromTrj(char x_or_z);
	int ComputeQuadPhaseTermsFromTrj(const SRWLPrtTrj& trj);

	inline int ComputeQuadPhaseTermsFromTrj();
	inline void TrjCoordAngField(double s, char x_or_z, double& x, double& dxds, double& B);
	//inline void FindOffestAndRelArg(double Arg, srTWaveAccessDataD1D& TrjInData, long& Offset, double& RelArg);
	inline void FindOffestAndRelArg(double Arg, srTWaveAccessDataD1D& TrjInData, long long& Offset, double& RelArg);
	//inline void FindOffestAndRelArgFromTrj(double Arg, const SRWLPrtTrj& trj, long& Offset, double& RelArg);
	inline void FindOffestAndRelArgFromTrj(double Arg, const SRWLPrtTrj& trj, long long& Offset, double& RelArg);
	inline void InterpFuncAndDerivs(double h, double x, DOUBLE* pf0, double& f, double& dfdx, double& d2fdx2, double& d3fdx3);
	inline int AllocateQuadPhaseTermsArrFromTrj();
	inline int AllocateQuadPhaseTermsArrFromTrj(const SRWLPrtTrj& trj);
	inline void DeallocateQuadPhaseTermsArr();
	inline int CheckIfFieldIsZeroFromTrj();
	int FieldComponIsZero_FromTrj(char x_or_z);

	void CountFieldExtremums();
	//void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep); //virtual
	void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, long long ns, double sStart, double sStep); //virtual

	//double EstimateIntermedVal(double* pFunc, int Ind, double sStp, double ds)
	double EstimateIntermedVal(double* pFunc, long long Ind, double sStp, double ds)
	{
		if((pFunc == 0) || (Ind < 0)) return 0;
		double f1 = pFunc[Ind];
		if((ds <= 0.) || (sStp <= 0.)) return f1;
		double f2 = pFunc[Ind + 1];
		return f1 + (ds*(f2 - f1))/sStp;
	}

	//long EstimMinNpForRadInteg(char typeInt) //virtual 
	long long EstimMinNpForRadInteg(char typeInt) //virtual 
	{//typeInt == 1: monochromatic emission in frequency domain
	 //typeInt == 2: power density (integral over all photon energies)
		int npMin = 5; //to tune
		if(m_estimMinNpForRadInteg > npMin) return m_estimMinNpForRadInteg;
		else return npMin; //to be re-defined in derived classes
	}
};

//*************************************************************************

inline double srTTrjDat::Derivative(double* f, double h, int PoIndx, int AmOfPo)
{
	if(AmOfPo==5)
	{
		if(PoIndx==2) return 0.08333333333333*(f[0] - 8.*f[1] + 8.*f[3] - f[4])/h;
		else if(PoIndx==1) return 0.08333333333333*(-3.*f[0] - 10.*f[1] + 18.*f[2] - 6.*f[3] + f[4])/h;
		else if(PoIndx==3) return 0.08333333333333*(-f[0] + 6.*f[1] - 18.*f[2] + 10.*f[3] + 3.*f[4])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==4) return 0.5*(f[2] - 4.*f[3] + 3.*f[4])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==4)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==2) return 0.5*(-f[1] + f[3])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==3) return 0.5*(f[1] - 4.*f[2] + 3.*f[3])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==3)
	{
		if(PoIndx==1) return 0.5*(-f[0] + f[2])/h;
		else if(PoIndx==0) return 0.5*(-3.*f[0] + 4.*f[1] - f[2])/h;
		else if(PoIndx==2) return 0.5*(f[0] - 4.*f[1] + 3.*f[2])/h;
		else return 1.E+23;
	}
	else if(AmOfPo==2) return (-f[0] + f[1])/h;
	else return 1.E+23;
}

//*************************************************************************

inline void srTTrjDat::CubPln(double f1, double f2, double fpr1, double fpr2, double* aa)
{
	double f1mf2_d_s1ms2 = (f2 - f1)/sStep;
	aa[0] = f1;
	aa[1] = fpr1;
	aa[2] = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)/sStep;
	aa[3] = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)/(sStep*sStep);
}

//*************************************************************************

inline double srTTrjDat::Field(char FieldID, double s)
{
	int Indx = int((s - sStart)/sStep);
	double sb = sStart + Indx*sStep;
	double smsb = s - sb;
	double smsbe2 = smsb*smsb;
	double* cPt = (FieldID=='x')? BxPlnCf[Indx] : ((FieldID=='z')? BzPlnCf[Indx] : 0);

	return (cPt!=0)? (cPt[0] + (cPt[1])*smsb + (cPt[2])*smsbe2 + (cPt[3])*smsbe2*smsb) : 0.;
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPoint(char x_or_z, double s, double& Bt, double& Crd, double& IntBtE2)
{
	if(CompFromTrj) { CompTrjDataDerivedAtPoint_FromTrj(x_or_z, s, Bt, Crd, IntBtE2); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	double sb = sStart + Indx*sStep;
	double smsb = s - sb;

	double *Bt_CfP, *C_CfP, *IntBt2_CfP;
	if(x_or_z=='x')
	{
		if(VerFieldIsNotZero)
		{
			Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];
			Bt = BtxCorr + BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
			double BufCrd = BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
			Crd = xCorr + BtxCorrForX*s + BufCrd;
			IntBtE2 = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]);
		}
		else { Bt = dxds0; Crd = x0 + dxds0*(s - s0); IntBtE2 = dxds0*dxds0*(s - s0);}
	}
	else if(x_or_z=='z')
	{
		if(HorFieldIsNotZero)
		{
			Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];
			Bt = BtzCorr - BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
			double BufCrd = -BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
			Crd = zCorr + BtzCorrForZ*s + BufCrd;
			IntBtE2 = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]);
		}
		else { Bt = dzds0; Crd = z0 + dzds0*(s - s0); IntBtE2 = dzds0*dzds0*(s - s0);}
	}
	else return;
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPoint_FromTrj(char x_or_z, double s, double& Bt, double& Crd, double& IntBtE2)
{
	double dBdsDummy, BDummy;
	CompTrjDataAndFieldWithDerAtPoint_FromTrj(x_or_z, s, dBdsDummy, BDummy, Bt, Crd, IntBtE2);
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPoint(double s, double& Btx, double& Crdx, double& IntBtxE2, double& Btz, double& Crdz, double& IntBtzE2)
{
	if(CompFromTrj) { CompTrjDataDerivedAtPoint_FromTrj(s, Btx, Crdx, IntBtxE2, Btz, Crdz, IntBtzE2); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//int Indx = int((s - sStart)*Inv_Step); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	long long Indx = (long long)((s - sStart)*Inv_Step); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	double sb = sStart + Indx*sStep;
	double smsb = s - sb;
	double *Bt_CfP, *C_CfP, *IntBt2_CfP;
	if(VerFieldIsNotZero)
	{
		Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];
		Btx = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
		double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
		Crdx = xCorr + BtxCorrForX*s + BufCrd;
		IntBtxE2 = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
	}
	else { Btx = dxds0; Crdx = x0 + dxds0*(s - s0); IntBtxE2 = dxds0*dxds0*(s - s0);}
	if(HorFieldIsNotZero)
	{
		Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];
		Btz = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
		double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
		Crdz = zCorr + BtzCorrForZ*s + BufCrd;
		IntBtzE2 = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
	}
	else { Btz = dzds0; Crdz = z0 + dzds0*(s - s0); IntBtzE2 = dzds0*dzds0*(s - s0);}
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPoint_FromTrj(double s, double& Btx, double& Crdx, double& IntBtxE2, double& Btz, double& Crdz, double& IntBtzE2)
{
	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	double sr, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//int Indx;
	long long Indx;

	if(VerFieldIsNotZero)
	{
		//int Indx = int((s - xTrjInData.Start)*xTrjInData.InvStep); 
		Indx = (long long)((s - xTrjInData.Start)*xTrjInData.InvStep); 
		if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
		if(Indx < 2) sr += (Indx - 2)*xTrjInData.Step;
		else if(Indx < xTrjInData.np - 3) ;
		else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
		else sr += (xTrjInData.Step + xTrjInData.Step);
		
		pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pIntBtE2_Cf = *(IntBtx2PlnCf+Indx);
		IntBtxE2 = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
		Crdx = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		Btx = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
	}
	else { double Buf = dxds0*(s - s0); Btx = dxds0; Crdx = x0 + Buf; IntBtxE2 = dxds0*Buf;}
	if(HorFieldIsNotZero)
	{
		//Indx = int((s - zTrjInData.Start)*zTrjInData.InvStep); 
		Indx = (long long)((s - zTrjInData.Start)*zTrjInData.InvStep); 
		if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
		if(Indx < 2) sr += (Indx - 2)*zTrjInData.Step;
		else if(Indx < zTrjInData.np - 3) ;
		else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
		else sr += (zTrjInData.Step + zTrjInData.Step);
		
		pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pIntBtE2_Cf = *(IntBtz2PlnCf+Indx);
		IntBtzE2 = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
		Crdz = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		Btz = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
	}
	else { double Buf = dzds0*(s - s0); Btz = dzds0; Crdz = z0 + Buf; IntBtzE2 = dzds0*Buf;}
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPointPowDens(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz)
{
	if(CompFromTrj) { CompTrjDataDerivedAtPointPowDens_FromTrj(s, Btx, Btz, X, Z, Bx, Bz); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//int Indx = int((s - sStart)*Inv_Step); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	long long Indx = (long long)((s - sStart)*Inv_Step); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	double sb = sStart + Indx*sStep;
	double smsb = s - sb;

	double *B_CfP, *Bt_CfP, *C_CfP;
	if(VerFieldIsNotZero)
	{
		B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx];
		Bz = *B_CfP + smsb*(B_CfP[1] + smsb*(B_CfP[2] + smsb*(B_CfP[3])));
		Btx = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
		double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
		X = (xCorr + BtxCorrForX*s + BufCrd);
	}
	else 
	{ 
		Bz = 0.; Btx = dxds0; 
		X = (x0 + dxds0*(s - s0));
	}
	if(HorFieldIsNotZero)
	{
		B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx];
		Bx = *B_CfP + smsb*(B_CfP[1] + smsb*(B_CfP[2] + smsb*(B_CfP[3])));
		Btz = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
		double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
		Z = (zCorr + BtzCorrForZ*s + BufCrd);
	}
	else 
	{ 
		Bx = 0.; Btz = dzds0; 
		Z = (z0 + dzds0*(s - s0));
	}
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataDerivedAtPointPowDens_FromTrj(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz)
{
	double *pB_Cf, *pBt_Cf, *pCrd_Cf;

	//int Indx = int((s - xTrjInData.Start)*xTrjInData.InvStep); 
	long long Indx = (long long)((s - xTrjInData.Start)*xTrjInData.InvStep); 
	if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
	if(Indx < 0) Indx = 0;
	double sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
	if(Indx < 2) sr += (Indx - 2)*xTrjInData.Step;
	else if(Indx < xTrjInData.np - 3) ;
	else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
	else sr += (xTrjInData.Step + xTrjInData.Step);

	pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pB_Cf = *(BzPlnCf+Indx);
	X = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
	Btx = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
	Bz = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));

	//Indx = int((s - zTrjInData.Start)*zTrjInData.InvStep); 
	Indx = (long long)((s - zTrjInData.Start)*zTrjInData.InvStep); 
	if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
	if(Indx < 0) Indx = 0;
	sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
	if(Indx < 2) sr += (Indx - 2)*zTrjInData.Step;
	else if(Indx < zTrjInData.np - 3) ;
	else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
	else sr += (zTrjInData.Step + zTrjInData.Step);

	pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pB_Cf = *(BxPlnCf+Indx);
	Z = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
	Btz = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
	Bx = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataAndFieldWithDerAtPoint(char x_or_z, double s, double& dBds, double& B, double& Bt, double& Crd, double& IntBtE2)
{
	if(CompFromTrj) { CompTrjDataAndFieldWithDerAtPoint_FromTrj(x_or_z, s, dBds, B, Bt, Crd, IntBtE2); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	double sb = sStart + Indx*sStep;
	double smsb = s - sb;
	double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;
	if(x_or_z=='x')
	{
		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];
			dBds = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
			B = B_CfP[0] + smsb*(B_CfP[1] + smsb*(B_CfP[2] + smsb*(B_CfP[3])));
			Bt = BtxCorr + BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
			double BufCrd = BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
			Crd = xCorr + BtxCorrForX*s + BufCrd;
			IntBtE2 = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]);
		}
		else { dBds = 0.; B = 0.; Bt = dxds0; Crd = x0 + dxds0*(s - s0); IntBtE2 = dxds0*dxds0*(s - s0);}
	}
	else if(x_or_z=='z')
	{
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];
			dBds = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
			B = B_CfP[0] + smsb*(B_CfP[1] + smsb*(B_CfP[2] + smsb*(B_CfP[3])));
			Bt = BtzCorr - BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
			double BufCrd = -BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
			Crd = zCorr + BtzCorrForZ*s + BufCrd;
			IntBtE2 = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]);
		}
		else { dBds = 0.; B = 0.; Bt = dzds0; Crd = z0 + dzds0*(s - s0); IntBtE2 = dzds0*dzds0*(s - s0);}
	}
	else return;
}

//*************************************************************************

//inline void srTTrjDat::ComponentArray(int ComponNo, char x_or_z, double sSt, double sEn, int Ns, double* CompArray)
inline void srTTrjDat::ComponentArray(int ComponNo, char x_or_z, double sSt, double sEn, long long Ns, double* CompArray)
{
	double h = (sEn - sSt)/(Ns - 1);
	double s = sSt;
	//for(int i=0; i<Ns; i++)
	for(long long i=0; i<Ns; i++)
	{
		if(ComponNo==0)
		{
			CompArray[i] = Field(x_or_z, s);
		}
		else
		{
			double Bt, Crd, IntBtE2;
			CompTrjDataDerivedAtPoint(x_or_z, s, Bt, Crd, IntBtE2);

			CompArray[i] = (ComponNo==1)? Bt : ((ComponNo==2)? Crd : ((ComponNo==3)? IntBtE2 : 0.));
		}
		s += h;
	}
}

//*************************************************************************

inline double srTTrjDat::Pol04(double s, double* c)
{
	return s*(*c + s*(*(c+1) + s*(*(c+2) + s*(*(c+3)))));
}

//*************************************************************************

inline double srTTrjDat::Pol05(double s, double* c)
{
	return s*(*c + s*(*(c+1) + s*(*(c+2) + s*(*(c+3) + s*(*(c+4))))));
}

//*************************************************************************

inline double srTTrjDat::Pol09(double s, double* c)
{
	return s*(*c + s*(*(c+1) + s*(*(c+2) + s*(*(c+3) + s*(*(c+4) + s*(*(c+5) + s*(*(c+6) + s*(*(c+7) + s*c[8]))))))));
}

//*************************************************************************

inline void srTTrjDat::CompCorrectionsForTrjDataDerived()
{
	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//int Indx = int((s0 - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	long long Indx = (long long)((s0 - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
	double sb = sStart + Indx*sStep;
	double smsb = s0 - sb;
	double *Bt_CfP, *C_CfP, *IntBt2_CfP;

	xCorr = BtxCorr = zCorr = BtzCorr = IntBtxE2Corr = IntBtzE2Corr = BtxCorrForX = BtzCorrForZ = 0.;
	if(VerFieldIsNotZero)
	{
		Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];

		BtxCorrForX = dxds0 - BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
		BtxCorrForXe2 = BtxCorrForX*BtxCorrForX;
		BtxCorr = BtxCorrForX;

		double BufX = BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
		xCorr = x0 - (BtxCorrForX*s0 + BufX);
		IntBtxE2Corr = -(BtxCorrForXe2*s0 + 2.*BtxCorrForX*BufX + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]));
	}
	else
	{
		BtxCorr = dxds0;
		xCorr = x0;
	}

	if(HorFieldIsNotZero)
	{
		Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];
		
		BtzCorrForZ = dzds0 + BetaNormConst*(Pol04(smsb, &(Bt_CfP[1])) + Bt_CfP[0]);
		BtzCorrForZe2 = BtzCorrForZ*BtzCorrForZ;
		BtzCorr = BtzCorrForZ;

		double BufZ = -BetaNormConst*(Pol05(smsb, &(C_CfP[1])) + C_CfP[0]);
		zCorr = z0 - (BtzCorrForZ*s0 + BufZ);
		IntBtzE2Corr = -(BtzCorrForZe2*s0 + 2.*BtzCorrForZ*BufZ + BetaNormConstE2*(Pol09(smsb, &(IntBt2_CfP[1])) + IntBt2_CfP[0]));
	}
	else
	{
		BtzCorr = dzds0;
		zCorr = z0;
	}
}

//*************************************************************************

inline void srTTrjDat::SetupFldPlnCf(srTFunDer* InitialFieldData, double** FldPlnCf)
{
	//int LenFieldData_m_1 = LenFieldData - 1; //OC020110

	double f1 = InitialFieldData[0].f, f2;
	double fpr1 = InitialFieldData[0].dfds, fpr2;

	//for(int is=1; is<LenFieldData; is++)
	for(long long is=1; is<LenFieldData; is++)
	{
		f2 = InitialFieldData[is].f;
		fpr2 = InitialFieldData[is].dfds;
		CubPln(f1, f2, fpr1, fpr2, FldPlnCf[is-1]);
		f1 = f2; fpr1 = fpr2;
	}
 }

//*************************************************************************

inline int srTTrjDat::ComputeInterpolatingStructure()
{
	//double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;
	double &s0 = EbmDat.s0; //OC020110

	if((s0 < sStart) || (s0 > sStart + sStep*(LenFieldData - 1))) return S0_OUT_OF_FIELD_DEFINITION_LIMITS;

	int result;
	if(result = EvaluateMemAvailBeforeTrjComp()) return result;

	Inv_Step = 1./sStep;

	m_estimMinNpForRadInteg = -1; //OC220112

	MakeZeroInputFieldDer();
	if(result = AllocateMemoryForCfs()) return result;
	if(HorFieldIsNotZero) 
	{
		CompDerivForFieldData(BxInData);
		SetupFldPlnCf(BxInData, BxPlnCf);
		SetupIntegrPlnCfs('z'); // This is not an error: Horizontal field is responsible for Vertical dynamics
	}
	if(VerFieldIsNotZero)
	{
		CompDerivForFieldData(BzInData);
		SetupFldPlnCf(BzInData, BzPlnCf);
		SetupIntegrPlnCfs('x');
	}
	CompBetaNormConst();
	CompCorrectionsForTrjDataDerived();

	CountFieldExtremums();

	LastCompWasNotOK = 0; 
	return 0;
}

//*************************************************************************

inline int srTTrjDat::ComputeInterpolatingStructure_FromTrj()
{
	int result;
	if((xTrjInData.pData == 0) || (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;

	if(result = ComputeQuadPhaseTermsFromTrj()) return result;
	if(result = AllocateMemoryForCfs_FromTrj()) return result;

	xTrjInData.SetupBufVars(); zTrjInData.SetupBufVars();

	CompBetaNormConst();
	m_estimMinNpForRadInteg = -1; //OC220112

	if(result = ComputeInterpolatingStructure_FromTrj1D('x')) return result;
	if(result = ComputeInterpolatingStructure_FromTrj1D('z')) return result;

	DeallocateQuadPhaseTermsArr();
	LastCompWasNotOK = 0; 
	return 0;
}

//*************************************************************************

inline void srTTrjDat::MakeZeroInputFieldDer()
{
	srTFunDer *tBxInData = BxInData, *tBzInData = BzInData;
	//for(int k=0; k<LenFieldData; k++)
	for(long long k=0; k<LenFieldData; k++)
	{
		(tBxInData++)->dfds = 0.;
		(tBzInData++)->dfds = 0.;
	}
}

//*************************************************************************

inline int srTTrjDat::EvaluateMemAvailBeforeTrjComp()
{
	//double CurrMemAvail = MaxMemAvail*1.E+06 - LenFieldData*450. - 0.8*1.E+06;
	//if(CurrMemAvail < 0.) { LastCompWasNotOK = 1; return NOT_ENOUGH_MEMORY_FOR_SR_COMP;}
	//else return 0;

#ifdef WIN32
	return 0;
#endif
#ifdef __MAC__
	double CurrMemAvail = srTSystemUtils::CheckMemoryAvailable();
	double RequiredMem = LenFieldData*450. + 0.8*1.E+06;
	if(CurrMemAvail*2. < RequiredMem) { LastCompWasNotOK = 1; return NOT_ENOUGH_MEMORY_FOR_SR_COMP;}
	else return 0;
#endif

}

//*************************************************************************

inline int srTTrjDat::SetUpFieldBasedArraysTotal(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{
	//long Ns = LenFieldData >> 1; // To steer
	long long Ns = LenFieldData >> 1; // To steer
	FieldBasedArrays.sStart = sStart;
	FieldBasedArrays.sStep = sStep;
	FieldBasedArrays.Ns = Ns;
	FieldBasedArrays.Nper = 1;

	int result;
	if(result = FieldBasedArrays.AllocateArrays(Ns, Keys)) return result;
	CompTotalTrjData(Keys, FieldBasedArrays);
	return 0;
}

//*************************************************************************

inline int srTTrjDat::ShowLimitsAndInitInteg(srTWfrSmp&, char LongIntType, double& sIntegStart, double& sIntegFin, int& AmOfPer, bool doInit) 
{
	sIntegStart = sStart;
	sIntegFin = sStart + (LenFieldData - 1)*sStep;
	AmOfPer = 1;

	if(CompFromTrj == 1) return 0; //OC220112: for SRWLib (since it calculates interpolating structure from trajectory only)

	if(doInit) return ComputeInterpolatingStructure(); //OC030412
	else return 0;
}

//*************************************************************************

inline void srTTrjDat::TrjCoordAngField(double s, char x_or_z, double& x, double& dxds, double& B)
{// for computing SR from Trajectory
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData; 
	//long StartInd;
	long long StartInd;
	double sr, d2xds2, d3xds3;
	FindOffestAndRelArg(s, TrjInData, StartInd, sr);

	DOUBLE *pfArr = TrjInData.pData + StartInd;
	InterpFuncAndDerivs(TrjInData.Step, sr, pfArr, x, dxds, d2xds2, d3xds3);

	B = d2xds2*InvBetaNormConst;
	if(x_or_z != 'x') B = -B;
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataAndFieldWithDerAtPoint_FromTrj(char x_or_z, double s, double& dBds, double& B, double& Bt, double& Crd, double& IntBtE2)
{// for computing SR from Trajectory
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData; 

	//int Indx = int((s - TrjInData.Start)/TrjInData.Step); 
	long long Indx = (long long)((s - TrjInData.Start)/TrjInData.Step); 
	if(Indx >= TrjInData.np - 1) Indx = TrjInData.np - 2;
	if(Indx < 0) Indx = 0;

	double sr = s - (TrjInData.Start + TrjInData.Step*Indx);
	if(Indx < 2) sr -= (2 - Indx)*TrjInData.Step;
	else if(Indx < TrjInData.np - 3) ;
	else if(Indx < TrjInData.np - 2) sr += TrjInData.Step;
	else sr += 2*TrjInData.Step;

	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	if(x_or_z == 'x')
	{
		pB_Cf = *(BzPlnCf+Indx); pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pIntBtE2_Cf = *(IntBtx2PlnCf+Indx);
	}
	else
	{
		pB_Cf = *(BxPlnCf+Indx); pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pIntBtE2_Cf = *(IntBtz2PlnCf+Indx);
	}

	IntBtE2 = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
	Crd = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
	Bt = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
	B = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
	dBds = *(pB_Cf+1) + sr*(*(pB_Cf+2)*2 + sr*(*(pB_Cf+3)*3));
}

//*************************************************************************

inline void srTTrjDat::CompTrjDataAndFieldWithDerAtPoint_FromTrjInitial(char x_or_z, double s, double& dBds, double& B, double& Bt, double& Crd, double& IntBtE2)
{// for computing SR from Trajectory
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData; 
	DOUBLE *pIntBtE2Arr = (x_or_z == 'x')? IntBtxE2Arr : IntBtzE2Arr; 

	//long StartInd;
	long long StartInd;
	double sr, d2xds2, d3xds3;
	FindOffestAndRelArg(s, TrjInData, StartInd, sr);

	DOUBLE *pfArr = TrjInData.pData + StartInd;
	InterpFuncAndDerivs(TrjInData.Step, sr, pfArr, Crd, Bt, d2xds2, d3xds3);

	B = d2xds2*InvBetaNormConst;
	dBds = d3xds3*InvBetaNormConst;
	if(x_or_z != 'x') { B = -B; dBds = -dBds;}

	double DummyDer1, DummyDer2, DummyDer3;
	pfArr = pIntBtE2Arr + StartInd;
	InterpFuncAndDerivs(TrjInData.Step, sr, pfArr, IntBtE2, DummyDer1, DummyDer2, DummyDer3);
}

//*************************************************************************

//inline void srTTrjDat::FindOffestAndRelArg(double Arg, srTWaveAccessDataD1D& TrjInData, long& Offset, double& RelArg)
inline void srTTrjDat::FindOffestAndRelArg(double Arg, srTWaveAccessDataD1D& TrjInData, long long& Offset, double& RelArg)
{
	//int Indx = int((Arg - TrjInData.Start)/TrjInData.Step); 
	long long Indx = (long long)((Arg - TrjInData.Start)/TrjInData.Step); 
	if(Indx >= TrjInData.np - 1) Indx = TrjInData.np - 2;
	if(Indx < 0) Indx = 0;

	RelArg = Arg - (TrjInData.Start + TrjInData.Step*Indx);
	if(Indx < 2)
	{
		Offset = 0;
		RelArg -= (2 - Indx)*TrjInData.Step;
	}
	else if(Indx < TrjInData.np - 3)
	{
		Offset = Indx - 2;
	}
	else if(Indx < TrjInData.np - 2)
	{
		Offset = Indx - 3;
		RelArg += TrjInData.Step;
	}
	else
	{
		Offset = Indx - 4;
		RelArg += 2*TrjInData.Step;
	}
}

//*************************************************************************

//inline void srTTrjDat::FindOffestAndRelArgFromTrj(double Arg, const SRWLPrtTrj& trj, long& Offset, double& RelArg)
inline void srTTrjDat::FindOffestAndRelArgFromTrj(double Arg, const SRWLPrtTrj& trj, long long& Offset, double& RelArg)
{//assumes that sStart, sStep were already setup
	Offset = 0; RelArg = 0;

	//int np = trj.np;
	long long np = trj.np;
	//sStart = trj.ctStart + trj.partInitCond.z;
	//sStep = (trj.ctEnd - trj.ctStart)/(np - 1);
	//if(trj.arZ != 0)
	//{//use tabulated longitudinal position
	//	sStart = trj.arZ[0];
	//	sStep = (trj.arZ[trj.np - 1] - sStart)/(np - 1);
	//}

	//int Indx = int((Arg - sStart)/sStep); 
	long long Indx = (long long)((Arg - sStart)/sStep); 
	if(Indx >= np - 1) Indx = np - 2;
	else if(Indx < 0) Indx = 0;
	RelArg = Arg - (sStart + sStep*Indx);

	//OC150815 (commented-out)
	//if(trj.arZ != 0)
	//{//use tabulated longitudinal position
	//	if(trj.arZ[Indx] > Arg)
	//	{
	//		for(int ii=Indx; ii>=0; ii--)
	//		{
	//			if(trj.arZ[ii] <= Arg) { Indx = ii; break;}
	//		}
	//	}
	//	else if(trj.arZ[Indx + 1] <= Arg)
	//	{
	//		for(int ii=Indx; ii<np; ii++)
	//		{
	//			if(trj.arZ[ii] > Arg) { Indx = ii - 1; break;}
	//		}
	//	}
	//	if(Indx >= np - 1) Indx = np - 2;
	//	else if(Indx < 0) Indx = 0;
	//	double zs = trj.arZ[Indx];
	//	RelArg = Arg - zs; sStep = trj.arZ[Indx + 1] - zs;
	//}

	if(Indx < 2) { Offset = 0; RelArg -= (2 - Indx)*sStep;}
	else if(Indx < np - 3) { Offset = Indx - 2;}
	else if(Indx < np - 2) { Offset = Indx - 3; RelArg += sStep;}
	else { Offset = Indx - 4; RelArg += 2*sStep;}
}

//*************************************************************************

inline void srTTrjDat::InterpFuncAndDerivs(double h, double x, DOUBLE* pf0, double& f, double& dfdx, double& d2fdx2, double& d3fdx3)
{
	double f0 = *pf0, f1 = *(pf0+1), f2 = *(pf0+2), f3 = *(pf0+3), f4 = *(pf0+4), f5 = *(pf0+5);

	double hInv = 1./h;
	double hInvE2 = hInv*hInv;
	double d1 = 0.016666666666667*hInv;
	double d2 = 0.041666666666667*hInvE2;
	double d3 = d2*hInv;
	double d4 = d3*hInv;
	double d5 = 0.2*d4*hInv;

	double a0 = f2;
	double a1 = (3*f0 - 30*f1 - 20*f2 + 60*f3 - 15*f4 + 2*f5)*d1;
	double a2 = -(f0 - 16*f1 + 30*f2 - 16*f3 + f4)*d2;
	double a3 = -(f0 + f1 - 10*f2 + 14*f3 - 7*f4 + f5)*d3;
	double a4 = (f0 - 4*(f1 + f3) + 6*f2 + f4)*d4;
	double a5 = (-f0 + 5*(f1 - f4) + 10*(f3 - f2) + f5)*d5;

	double b0 = a1, b1 = 2*a2, b2 = 3*a3, b3 = 4*a4, b4 = 5*a5; 
	double c0 = b1, c1 = 2*b2, c2 = 3*b3, c3 = 4*b4; 

	f = a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5))));
	dfdx = b0 + x*(b1 + x*(b2 + x*(b3 + x*b4)));
	d2fdx2 = c0 + x*(c1 + x*(c2 + x*c3));
	d3fdx3 = c1 + x*(2*c2 + x*3*c3);
}

//*************************************************************************

inline void srTTrjDat::CoefsPol5thOrder(double h, DOUBLE* pf0, double* pCoef)
{
	double f0 = *pf0, f1 = *(pf0+1), f2 = *(pf0+2), f3 = *(pf0+3), f4 = *(pf0+4), f5 = *(pf0+5);

	double hInv = 1./h;
	double hInvE2 = hInv*hInv;
	double d1 = 0.016666666666667*hInv;
	double d2 = 0.041666666666667*hInvE2;
	double d3 = d2*hInv;
	double d4 = d3*hInv;
	double d5 = 0.2*d4*hInv;

	*pCoef = f2;
	pCoef[1] = (3*f0 - 30*f1 - 20*f2 + 60*f3 - 15*f4 + 2*f5)*d1;
	pCoef[2] = -(f0 - 16*f1 + 30*f2 - 16*f3 + f4)*d2;
	pCoef[3] = -(f0 + f1 - 10*f2 + 14*f3 - 7*f4 + f5)*d3;
	pCoef[4] = (f0 - 4*(f1 + f3) + 6*f2 + f4)*d4;
	pCoef[5] = (-f0 + 5*(f1 - f4) + 10*(f3 - f2) + f5)*d5;
}

//*************************************************************************

inline int srTTrjDat::AllocateQuadPhaseTermsArrFromTrj()
{
	DeallocateQuadPhaseTermsArr();
	if((xTrjInData.pData == 0) || (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;

	IntBtxE2Arr = new DOUBLE[xTrjInData.np];
	if(IntBtxE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;
	IntBtzE2Arr = new DOUBLE[zTrjInData.np];
	if(IntBtzE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;
	return 0;
}

//*************************************************************************

inline int srTTrjDat::AllocateQuadPhaseTermsArrFromTrj(const SRWLPrtTrj& trj)
{
	DeallocateQuadPhaseTermsArr();

	IntBtxE2Arr = new DOUBLE[trj.np];
	if(IntBtxE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;
	IntBtzE2Arr = new DOUBLE[trj.np];
	if(IntBtzE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;
	return 0;
}

//*************************************************************************

inline void srTTrjDat::DeallocateQuadPhaseTermsArr()
{
	if(IntBtxE2Arr != 0) delete[] IntBtxE2Arr; IntBtxE2Arr = 0;
	if(IntBtzE2Arr != 0) delete[] IntBtzE2Arr; IntBtzE2Arr = 0;
}

//*************************************************************************

inline int srTTrjDat::ComputeQuadPhaseTermsFromTrj()
{
	int result;
	if((xTrjInData.pData == 0) || (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;
	if(result = AllocateQuadPhaseTermsArrFromTrj()) return result;
	if(result = ComputeOneQuadPhaseTermFromTrj('x')) return result;
	if(result = ComputeOneQuadPhaseTermFromTrj('z')) return result;
	return 0;
}

//*************************************************************************

inline int srTTrjDat::CheckIfFieldIsZeroFromTrj()
{
	if((xTrjInData.pData == 0) || (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;
	VerFieldIsNotZero = FieldComponIsZero_FromTrj('x')? 0 : 1; // x_or_z defines Trajectory component !!
	HorFieldIsNotZero = FieldComponIsZero_FromTrj('z')? 0 : 1;
	return 0;
}

//*************************************************************************

#endif
