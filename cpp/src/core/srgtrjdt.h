/************************************************************************//**
 * File: srgtrjdt.h
 * Description: Electron trajectory calculation (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * Portions Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRGTRJDT_H
#define __SRGTRJDT_H

#include "srebmdat.h"
#include "srercode.h"
#include "srmagelem.h"
#include "gmmeth.h"
#include "gmintrk.h"
#include <vector>

//*************************************************************************

struct srTFieldBasedArrayKeys {
	char Bx_, Bz_, Btx_, Btz_, X_, Z_, IntBtxE2_, IntBtzE2_, dBxds_, dBzds_;
	char X1p_, Z1p_, X2p_, Z2p_, X1_, Z1_, X2_, Z2_;
    char IntX01_, IntX02_, IntX11_, IntX12_, IntX22_, IntZ01_, IntZ02_, IntZ11_, IntZ12_, IntZ22_;
	char IntX01toS0_, IntX02toS0_, IntX11toS0_, IntX12toS0_, IntX22toS0_;
    char IntZ01toS0_, IntZ02toS0_, IntZ11toS0_, IntZ12toS0_, IntZ22toS0_;

	srTFieldBasedArrayKeys()
	{
		ZeroAllKeys();
	}
	void ZeroAllKeys()
	{
		Bx_ = Bz_ = Btx_ = Btz_ = X_ = Z_ = IntBtxE2_ = IntBtzE2_ = dBxds_ = dBzds_ = 0;
        X1p_ = Z1p_ = X2p_ = Z2p_ = X1_ = Z1_ = X2_ = Z2_ = 0;
        IntX01_ = IntX02_ = IntX11_ = IntX12_ = IntX22_ = IntZ01_ = IntZ02_ = IntZ11_ = IntZ12_ = IntZ22_ = 0;
	}

    char AllEqData_ShouldBeSet()
	{
		return (Bx_ & Bz_ & Btx_ & Btz_ & X_ & Z_ & IntBtxE2_ & IntBtzE2_ & dBxds_ & dBzds_);
	}
	char AllNonEqTrajData_ShouldBeSet()
	{
		return (X1p_ & Z1p_ & X2p_ & Z2p_ & X1_ & Z1_ & X2_ & Z2_);
	}
	char AllNonEqIntegData_ShouldBeSet()
	{
		return (IntX01_ & IntX02_ & IntX11_ & IntX12_ & IntX22_ & IntZ01_ & IntZ02_ & IntZ11_ & IntZ12_ & IntZ22_);
	}
    char AllNonEqData_ShouldBeSet()
	{
		return (AllNonEqTrajData_ShouldBeSet() & AllNonEqIntegData_ShouldBeSet());
	}
	char CheckIfCalcIsForThickBeam()
	{
		return (X1p_ | Z1p_ | X2p_ | Z2p_ | X1_ | Z1_ | X2_ | Z2_ | IntX01_ | IntX02_ | IntZ01_ | IntZ02_);
	}
	char CheckIfCalcCorrIntIsNecessary()
	{
		return (IntX01toS0_ | IntX02toS0_ | IntX11toS0_ | IntX12toS0_ | IntX22toS0_ | IntZ01toS0_ | IntZ02toS0_ | IntZ11toS0_ | IntZ12toS0_ | IntZ22toS0_);
	}
	char CheckIfNonEqNumIntegIsNecessary()
	{
        if((!X1p_) && (IntX01_ || IntX11_ || IntX12_)) return 1;
        else if((!X2p_) && (IntX02_ || IntX22_)) return 1;
        else if((!Z1p_) && (IntZ01_ || IntZ11_ || IntZ12_)) return 1;
        else if((!Z2p_) && (IntZ02_ || IntZ22_)) return 1;
		else return 0;
	}
};

//*************************************************************************

struct srTFieldBasedArrays {

	double *BxArr, *BzArr, *BtxArr, *BtzArr, *XArr, *ZArr, *IntBtxE2Arr, *IntBtzE2Arr, *dBxdsArr, *dBzdsArr;
	double *X1pArr, *Z1pArr, *X2pArr, *Z2pArr, *X1Arr, *Z1Arr, *X2Arr, *Z2Arr;
    double *IntX01Arr, *IntX02Arr, *IntX11Arr, *IntX12Arr, *IntX22Arr;
    double *IntZ01Arr, *IntZ02Arr, *IntZ11Arr, *IntZ12Arr, *IntZ22Arr;

	double absBxMax, absBzMax;

	double IntX01toS0, IntX02toS0, IntX11toS0, IntX12toS0, IntX22toS0;
    double IntZ01toS0, IntZ02toS0, IntZ11toS0, IntZ12toS0, IntZ22toS0;

	double sStart, sStep;
	//long Ns;
	//int Nper;
	//int NperLeft;
	long long Ns;
	long long Nper;
	long long NperLeft;

	srTFieldBasedArrays()
	{
        ZeroPtrs();

		Ns = 0;
		Nper = 1;
		NperLeft = 0;
		absBxMax = absBzMax = 0;
	}
	~srTFieldBasedArrays()
	{
		DisposeArrays();
	}

	void ZeroPtrs()
	{
		BxArr = BzArr = BtxArr = BtzArr = XArr = ZArr = IntBtxE2Arr = IntBtzE2Arr = dBxdsArr = dBzdsArr = 0;
        X1pArr = Z1pArr = X2pArr = Z2pArr = X1Arr = Z1Arr = X2Arr = Z2Arr = 0;
        IntX01Arr = IntX02Arr = IntX11Arr = IntX12Arr = IntX22Arr = 0;
        IntZ01Arr = IntZ02Arr = IntZ11Arr = IntZ12Arr = IntZ22Arr = 0;

        IntX01toS0 = IntX02toS0 = IntX11toS0 = IntX12toS0 = IntX22toS0 = 0;
        IntZ01toS0 = IntZ02toS0 = IntZ11toS0 = IntZ12toS0 = IntZ22toS0 = 0;
	}

	void DisposeArrays()
	{
		if(BxArr != 0) delete[] BxArr;
		if(BzArr != 0) delete[] BzArr;
		if(BtxArr != 0) delete[] BtxArr;
		if(BtzArr != 0) delete[] BtzArr;
		if(XArr != 0) delete[] XArr;
		if(ZArr != 0) delete[] ZArr;
		if(IntBtxE2Arr != 0) delete[] IntBtxE2Arr;
		if(IntBtzE2Arr != 0) delete[] IntBtzE2Arr;
		if(dBxdsArr != 0) delete[] dBxdsArr;
		if(dBzdsArr != 0) delete[] dBzdsArr;

		if(X1pArr != 0) delete[] X1pArr;
		if(Z1pArr != 0) delete[] Z1pArr;
		if(X2pArr != 0) delete[] X2pArr;
		if(Z2pArr != 0) delete[] Z2pArr;
		if(X1Arr != 0) delete[] X1Arr;
		if(Z1Arr != 0) delete[] Z1Arr;
		if(X2Arr != 0) delete[] X2Arr;
		if(Z2Arr != 0) delete[] Z2Arr;
		if(IntX01Arr != 0) delete[] IntX01Arr;
		if(IntX02Arr != 0) delete[] IntX02Arr;
		if(IntX11Arr != 0) delete[] IntX11Arr;
		if(IntX12Arr != 0) delete[] IntX12Arr;
		if(IntX22Arr != 0) delete[] IntX22Arr;
		if(IntZ01Arr != 0) delete[] IntZ01Arr;
		if(IntZ02Arr != 0) delete[] IntZ02Arr;
		if(IntZ11Arr != 0) delete[] IntZ11Arr;
		if(IntZ12Arr != 0) delete[] IntZ12Arr;
		if(IntZ22Arr != 0) delete[] IntZ22Arr;

		ZeroPtrs();
	}
	//int AllocateArrays(long InNs, srTFieldBasedArrayKeys& Keys);
	int AllocateArrays(long long InNs, srTFieldBasedArrayKeys& Keys);
};

//*************************************************************************

class srTGenTrjDat;
class srTMagElem;
class srTWfrSmp;
struct SRWLStructKickMatrix;
typedef struct SRWLStructKickMatrix SRWLKickM;
typedef CSmartPtr<srTGenTrjDat> srTGenTrjHndl;

//*************************************************************************

class srTGenTrjDat {

	double m_Mult2ndDerRK;
	CSmartPtr<srTMagElem> m_hMagElem; //SRWLIB

public:

	srTEbmDat EbmDat;

	short HorFieldIsNotZero, VerFieldIsNotZero;
	double BetaNormConst, BetaNormConstE2, InvBetaNormConst;

	//long AmOfExtremInBx, AmOfExtremInBz;
	long long AmOfExtremInBx, AmOfExtremInBz;

	srTGenTrjDat(srTEbmDat* pEbmDat)
	{
		HorFieldIsNotZero = VerFieldIsNotZero = 0;
		AmOfExtremInBx = AmOfExtremInBz = -1;
		if(pEbmDat != 0) EbmDat = *pEbmDat;
	}
	srTGenTrjDat(srTEbmDat* pEbmDat, CSmartPtr<srTMagElem>& hMagElem) //SRWLIB
	{
		HorFieldIsNotZero = VerFieldIsNotZero = 0;
		AmOfExtremInBx = AmOfExtremInBz = -1;
		
		if(pEbmDat != 0) EbmDat = *pEbmDat;
		m_hMagElem = hMagElem;
	}

	srTGenTrjDat() { AmOfExtremInBx = AmOfExtremInBz = -1;}
	virtual ~srTGenTrjDat() {}

	static srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat* pEbmDat, srTMagElem* pMagElem);

	virtual void CompTrjDataDerivedAtPointPowDens(double s, double& Btx, double& Btz, double& X, double& Z, double& Bx, double& Bz) {}
	virtual void CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays) {}
	//virtual void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz) {}
	virtual void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz) {}
	//virtual void CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds) {}
	virtual void CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds) {}

	virtual int SetUpFieldBasedArraysAtOnePeriod(srTFieldBasedArrayKeys&, srTFieldBasedArrays&) { return 0;}
	virtual int SetUpFieldBasedArraysTotal(srTFieldBasedArrayKeys&, srTFieldBasedArrays&) { return 0;}
	virtual int MagFieldPeriodicity() { return 1;}
	virtual char MagFieldIsConstant() { return 0;}
	virtual void AnalizeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ) { FieldIsSymOverX = FieldIsSymOverZ = 0;}

	virtual void ShowFullLimits(double& sIntegStart, double& sIntegFin) {}
	virtual int ShowLimitsAndInitInteg(srTWfrSmp& DistrInfoDat, char LongIntType, double& sIntegStart, double& sIntegFin, int& AmOfPer, bool doInit = true) { return 0;}
	virtual int InitTrjComp() { return 0;}
	//virtual int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long& Ns, int& NperTot, int& NperLeft) { return 0;}
	virtual int SetupLimitsByAnalizingField(char LongIntType, double& sStart, double& sStep, long long& Ns, int& NperTot, int& NperLeft) { return 0;}
	virtual int ConvertToArbTrjDat(char, srTWfrSmp&, srTGenTrjHndl&) { return 0;}

	//virtual void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep) {}; //virtual
	virtual void CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, long long ns, double sStart, double sStep) {}; //virtual

	//void CompTrjCrdVelRK(double sSt, double sEn, long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData);
	void CompTrjCrdVelRK(double sSt, double sEn, long long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData);
	//void CompTrjKickMatr(SRWLKickM* arKickM, int nKickM, double sSt, double sEn, long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData);
	void CompTrjKickMatr(SRWLKickM* arKickM, int nKickM, double sSt, double sEn, long long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData);
	//void IntegrateKicks(SRWLKickM* arKickM, vector<vector<int> >& vIndNonOverlapKickGroups, vector<pair<double, double> >& vIndNonOverlapKickGroupRanges, double inv_B_pho, double* initCond, double sStart, double sEnd, int ns, double* pTrjRes);
	void IntegrateKicks(SRWLKickM* arKickM, vector<vector<int> >& vIndNonOverlapKickGroups, vector<pair<double, double> >& vIndNonOverlapKickGroupRanges, double inv_B_pho, double* initCond, double sStart, double sEnd, long long ns, double* pTrjRes);

	//virtual long EstimMinNpForRadInteg(char typeInt)
	virtual long long EstimMinNpForRadInteg(char typeInt)
	{//typeInt == 1: monochromatic emission in frequency domain
	 //typeInt == 2: power density (integral over all photon energies)
		return 5; //to be re-defined in derived classes
	}

	void CompBetaNormConst()
	{// Assume lengths in m and field in Tesla at the input
	 // then output beta is in radians
	
		const double e = 1.60217646263E-19; //1.602189246E-19; // Charge of electron in Coulomb
		const double m = 9.1093818872E-31; //9.10953447E-31; // Mass of electron in kg
		const double c = 2.99792458E+08; // Speed of light in m/s

		BetaNormConst = -e/(EbmDat.Gamma*m*c);
		BetaNormConstE2 = BetaNormConst*BetaNormConst;
		InvBetaNormConst = 1./BetaNormConst;
	}

	//void CompTrjCrdVel(double sSt, double sEn, long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData)
	void CompTrjCrdVel(double sSt, double sEn, long long np, double* pInPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData)
	{
		int methNo = 1; //RK4
		//methNo = 2; //RK5
		double *pPrec=0;
		if(pInPrecPar != 0)
		{
			if(pInPrecPar[0] > 0)
			{
				methNo = (int)pInPrecPar[1];
				pPrec = pInPrecPar;
			}
		}

		if((methNo == 1) || (methNo == 2)) CompTrjCrdVelRK(sSt, sEn, np, pPrec, pOutBtxData, pOutXData, pOutBtyData, pOutYData, pOutBtzData, pOutZData, pOutBxData, pOutByData, pOutBzData);
		//else ...
	}

	void funcDerivRK(double s, double* arr_F, double* arr_dFds)
	{
		double xd = arr_F[1], yd = arr_F[3]; //, zd = arr_F[5];
		double zd = CGenMathMeth::radicalOnePlusSmall(-(EbmDat.GammaEm2 + xd*xd + yd*yd));
		TVector3d P(arr_F[0], arr_F[2], arr_F[4]), B;
		//m_Fld3d.compB(P, B);
		m_hMagElem.rep->compB(P, B);

		arr_dFds[0] = xd;
		arr_dFds[1] = m_Mult2ndDerRK*(yd*B.z - zd*B.y);
		arr_dFds[2] = yd;
		arr_dFds[3] = m_Mult2ndDerRK*(zd*B.x - xd*B.z);
		arr_dFds[4] = zd;
		//arr_dFds[5] = m_Mult2ndDer*(xd*B.y - yd*B.x);
	}

	//template<class T1, class T2> bool auxLessInPairBasedOnSecond(pair<T1, T2> p1, pair<T1, T2> p2)
	//{
	//	return p1.second < p2.second;
	//}

};

//*************************************************************************

#endif
