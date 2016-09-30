/************************************************************************//**
 * File: srcradint.h
 * Description: CSR calculation (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2006
 *
 * Copyright (C) Synchrotron SOLEIL, Gif-sur-Yvette, France
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SRCRADINT_H
#define __SRCRADINT_H

//#ifdef __IGOR_PRO__
//#ifndef __SRSEND_H
//#include "srsend.h"
//#endif
//#endif

#include "srtrjdat.h"
#include "srmagcnt.h"
#include "srradint.h"
#include "gmrand.h"

#include <complex>

//*************************************************************************

//typedef CSmartPtr<srTPartAutoRadInt> srTPartAutoRadIntHndl;
//
//#ifdef __GCC__
//typedef vector<srTPartAutoRadIntHndl> srTPartAutoRadIntHndlVect;
//#else
//typedef vector<srTPartAutoRadIntHndl, allocator<srTPartAutoRadIntHndl> > srTPartAutoRadIntHndlVect;
//#endif
//
//struct srTParPrecElecFld;

//*************************************************************************

struct TPrecParamCSR {

	int m_MethNo;
	double m_PrecPar;
	double m_sStep, m_sIntegStart, m_sIntegEnd;
	double m_NxNzOversampFact;
	int m_MC_NumMacroPart;

	void setupFromArray(double* _pdPar)
	{
		m_NxNzOversampFact = 0;
		if(_pdPar != 0)
		{
			m_MethNo = (int)_pdPar[0];
			if(m_MethNo == 0)
			{
                m_PrecPar = 0;
				m_sStep = _pdPar[1];
			}
			else
			{
                m_PrecPar = _pdPar[1];
				m_sStep = 0;
			}
			m_sIntegStart = _pdPar[2];
			m_sIntegEnd = _pdPar[3];

			if((_pdPar[4] != 0) && (_pdPar[5] > 0))
			{
				m_NxNzOversampFact = _pdPar[5];
			}

			if(_pdPar[6] > 0)
			{
				//m_MethNo = 10; // Monte-Carlo
				m_MC_NumMacroPart = (int)_pdPar[7];
			}
			else
			{
                m_MC_NumMacroPart = 0;
			}

		//Integration Parameters wave:
		// 0: Meth. No (=0 for simplest method)
		// 1: Prec. Param. (~1) or Step
		// 2: Long. Pos. to Start integ.
		// 3: Long. Pos. to Finish integ.
		// 4: Use auto-sampling for propag. (=0 means don't use)
		// 5: Over-sampling param.
		// 6: Use Monte-Carlo multi-particle integration
		// 7: Number of macro-particles for Monte-Carlo integration
		}
	}
};

//*************************************************************************

struct TAuxParamForIntegCSR {

	srTEFourier* arrEw;

	double pi, half_pi, sqrt_pi, half_sqrt_pi, two_sqrt_pi, piE2;
	double half_k_d_e, k_d_e;
	
	double pxx, pzz, pxx1, pzz1, px1x1, pz1z1;
	double pxz, px1z, pxz1, px1z1;
	double pxg, pzg, px1g, pz1g;
	double pxs, pzs, px1s, pz1s;
	double psg, pss, pgg;

	double constW1, constW2, constW3a, constW3b, constW4b, constW5, constQ6;
	double cf, cn;

	TAuxParamForIntegCSR()
	{
		arrEw = 0;
	}
	~TAuxParamForIntegCSR()
	{
		disposeArrays();
	}

	//void allocateArrays(int meth, long ns)
	void allocateArrays(int meth, long long ns)
	{
		disposeArrays();

		if(meth == 0)
		{
			if(ns <= 0) return;
            arrEw = new srTEFourier[ns];
		}
		else if(meth == 1) {}
		else if(meth == 2) {}
	}
	void disposeArrays()
	{
		if(arrEw != 0) { delete[] arrEw; arrEw = 0;}
	}
	void setupConstParams(srTEbmDat& eBeam);
};

//*************************************************************************

class srTCSR {

	srTTrjDat& m_TrjDat;
	srTSRWRadStructAccessData& m_Wfr;
	TPrecParamCSR m_PrecParams;
	srTFieldBasedArrays m_FldArr;
    TAuxParamForIntegCSR m_AuxIntPar;

	srTRadInt* m_pRadInt;
	srTWfrSmp* m_pWfrSmpAux;
	srTSRWRadStructAccessData* m_pWfrAux;
	srTParPrecElecFld* m_pParPrecElecFldSingle;
	srTTrjDat* m_pTrjDatAux;
	CGenMathRand m_gmRand;
	double m_xcArr[6], m_sigArr[6];

public:

	srTCSR(srTTrjDat& _TrjDat, double* _pdPrcPar, srTSRWRadStructAccessData& _Wfr) : m_TrjDat(_TrjDat), m_Wfr(_Wfr)
	{
		m_PrecParams.setupFromArray(_pdPrcPar);
		checkAndCorrectIntegLimits();
		m_pRadInt = 0; m_pWfrSmpAux = 0; m_pWfrAux = 0; m_pParPrecElecFldSingle = 0; m_pTrjDatAux = 0;
	}

	~srTCSR()
	{
		//DeallocateMemForRadDistr();	
	}

	void checkAndCorrectIntegLimits() 
	{
		double FieldDataFin = m_TrjDat.sStart + m_TrjDat.sStep*(m_TrjDat.LenFieldData - 1);
		if(m_PrecParams.m_sIntegStart < m_TrjDat.sStart) m_PrecParams.m_sIntegStart = m_TrjDat.sStart;
		if(m_PrecParams.m_sIntegEnd > FieldDataFin) m_PrecParams.m_sIntegEnd = FieldDataFin;
	}
	void genRadIntegration(srTEXZY& exzy, srTEFourier& Ew)
	{// Put here more functionality (switching to different methods) later

		if(m_PrecParams.m_MC_NumMacroPart > 0) { radIntegrationMonteCarlo(exzy, Ew); return;}

		srTEFourier EwResid, dEwdsAtEdges[2];

		radIntegrationResiduals(exzy, EwResid, dEwdsAtEdges);
		Ew.EwX_Re = Ew.EwX_Im = Ew.EwZ_Re = Ew.EwZ_Im = 0;

		if(m_PrecParams.m_MethNo == 0) radIntegrationManual(exzy, dEwdsAtEdges, Ew);
		else if(m_PrecParams.m_MethNo == 1) radIntegrationAutoUnd(exzy, dEwdsAtEdges, Ew);
		else if(m_PrecParams.m_MethNo == 2) radIntegrationAutoWig(exzy, dEwdsAtEdges, Ew);
		
		Ew += EwResid;
		Ew *= m_AuxIntPar.cn*exzy.e;
	}
	complex<double> computeResidualRow(complex<double>* dAds, complex<double>* dBds, int numTerms)
	{//assumes dAds[0] = A, dAds[1] = dA/ds; dBds[0] = B, dBds[1] = dB/ds;
	 //add analysis of accuracy here
		complex<double> inv_dBds = 1./(*(dBds + 1));
		complex<double> multC = exp(*dBds)*inv_dBds;

		complex<double> rowC = *dAds;
		if(numTerms > 1) 
		{
			rowC += (-dAds[1] + ((*dAds)*dBds[2]*inv_dBds))*inv_dBds;
		}
		if(numTerms > 2) 
		{
			rowC += (dAds[2] + (((3.*(((*dAds)*dBds[2]*dBds[2]*inv_dBds) - (dAds[1]*dBds[2]))) - ((*dAds)*dBds[3]))*inv_dBds))*inv_dBds*inv_dBds;
		}
		return rowC*multC;
	}
	//complex<double> sqrtC(complex<double>& c)
	complex<double> sqrtC(complex<double> c) //OC030110 required for GCC 4.2
	{
		double rec = c.real(), imc = c.imag();
		double amp = sqrt(sqrt(rec*rec + imc*imc));
		double arg = (rec == 0)? 0 : 0.5*atan(imc/rec);
		complex<double> res(amp*cos(arg), amp*sin(arg));
		return res;
	}
	complex<double> invSqrtC(complex<double>& c)
	{
		double rec = c.real(), imc = c.imag();
		double amp = 1./sqrt(sqrt(rec*rec + imc*imc));
		double arg = (rec == 0)? 0 : -0.5*atan(imc/rec);
		complex<double> res(amp*cos(arg), amp*sin(arg));
		return res;
	}

	void checkInputConsistency();
	void estimateAbsoluteTolerance();
	void performMethodDependentSetupActions();
	void performMethodDependentFinishActions();
	//void createAuxWfrSmp();
	void analyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ);
	void fillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ);
    void copySymEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, char SymWithRespectToXax, char SymWithRespectToZax);
	void radIntegrationManual(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew);
	void radIntegrationAutoUnd(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew);
	void radIntegrationAutoWig(srTEXZY& exzy, srTEFourier* arr_dEwds, srTEFourier& Ew);
    void radIntegrationMonteCarlo(srTEXZY& exzy, srTEFourier& Ew);

	void radIntegrationResiduals(srTEXZY& exzy, srTEFourier& Ew, srTEFourier* dEwds);
	void computeTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont);
	void computeOffAxisTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont);
	void setupInitialTrajArrays(srTMagFldCont* pMagLensCont);
	//void integrateSimpleEwArr(srTEFourier* arrEw, long np, double h, srTEFourier* pDer, srTEFourier& resEw);
	void integrateSimpleEwArr(srTEFourier* arrEw, long long np, double h, srTEFourier* pDer, srTEFourier& resEw);
	//void computeFuncToIntegAtOnePointOnTrj(long i, srTEXZY exzy, srTEFourier& Ew, complex<double>& ampX, complex<double>& ampZ, complex<double>& arg);
	void computeFuncToIntegAtOnePointOnTrj(long long i, srTEXZY exzy, srTEFourier& Ew, complex<double>& ampX, complex<double>& ampZ, complex<double>& arg);

    void computeElectricFieldFreqDomain();
};

//*************************************************************************

#endif
