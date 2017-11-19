/************************************************************************//**
 * File: srmagfld.h
 * Description: Magnetic elements (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * Portions Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMAGFLD_H
#define __SRMAGFLD_H

#include <stdlib.h>

#include "srmagelem.h"
#include "gmvect.h"
#include <vector>

#include "srobject.h"
#include "srtrjdat.h"
#include "srercode.h"
#include "srstraux.h"
#include "objcont.h"
#include "gminterp.h"

//*************************************************************************

class srTGenTrjDat;
class srTTrjDat;
class srTTrjDat3d;
class srTMagFieldPeriodic;
class srTMagFldCont;
struct srTWigComSASE;

struct SRWLStructMagneticFieldUndulator;
typedef struct SRWLStructMagneticFieldUndulator SRWLMagFldU;

//*************************************************************************

class srTMagHarm : public srTMagElem {
public:

	int HarmNo;
	char XorZ; //'x' or 'z'
	double K;
	double Phase;
	int s; //symmetry vs longitudinal position: 1 - symmetric, -1 - anti-symmetric
	double TrA; //coefficient for transverse dependence

	srTMagHarm(int In_HarmNo, char In_XorZ, double In_K, double In_Phase, int In_s=1, double In_TrA=1)
	{
		HarmNo = In_HarmNo; XorZ = In_XorZ; K = In_K; Phase = In_Phase; s = In_s; TrA = In_TrA;
	}
	srTMagHarm() {}

	void InData(int InHarmNo, char InXorZ, double InK, double InPhase)
	{//used by DLL only
		HarmNo = InHarmNo; K = InK; Phase = InPhase;
		if((InXorZ == 'x') || (InXorZ == 'h')) XorZ = 'x';
		else if((InXorZ == 'z') || (InXorZ == 'v')) XorZ = 'z';
	}

	void OutData(int* pHarmNo, char* pXorZ, double* pK, double* pPhase)
	{//used by DLL only
		*pHarmNo = HarmNo; *pK = K; *pPhase = Phase;
		if(XorZ == 'x') *pXorZ = 'h';
		else if(XorZ == 'z') *pXorZ = 'v';
	}

	inline friend bool operator <(const srTMagHarm&, const srTMagHarm&);
	inline friend bool operator ==(const srTMagHarm&, const srTMagHarm&);
	inline friend bool operator >(const srTMagHarm&, const srTMagHarm&); //OC071214
};

inline bool operator <(const srTMagHarm& P1, const srTMagHarm& P2)
{
	return P1.K < P2.K;
}

inline bool operator ==(const srTMagHarm& P1, const srTMagHarm& P2)
{
	return (P1.HarmNo < P2.HarmNo) && (P1.XorZ < P2.XorZ) && (P1.K < P2.K) && (P1.Phase < P2.Phase);
}

inline bool operator >(const srTMagHarm& P1, const srTMagHarm& P2) //OC071214
{
	return P1.K > P2.K;
}

//*************************************************************************

typedef vector<srTMagHarm, allocator<srTMagHarm> > srTMagHarmVect;

//*************************************************************************

class srTMagFldArb1D : public srTMagElem {

	double *pB;
	int Np;
	double sStart, sStep;
	char v_or_h;

public:

	srTMagFldArb1D(char* Name, char in_v_or_h, double in_s0, double in_ds, int in_np, double* in_pB) : srTMagElem(Name)
	{
        pB = in_pB; Np = in_np;
        sStart = in_s0; sStep = in_ds;
        v_or_h = in_v_or_h;

		gsStart = in_s0; 
		gsEnd = in_s0 + (in_np - 1)*in_ds;
	}
	srTMagFldArb1D() 
	{
		pB = 0; Np = 0;
		gsStart = gsEnd = 0.;
	}
	//~srTMagFldArb1D() 
	//{
	//	DeleteArray();
	//}
	void DeleteArray()
	{
		if(pB != 0) { delete[] pB; pB = 0;}
	}
};

//*************************************************************************

class srTMagFldTrUnif : public srTMagElem {

	double *BxArr, *BzArr;
	int Np;
	double sStart, sStep;
	char ArraysWereAllocated;

public:

	srTMagFldTrUnif(double In_sStart, double In_sStep, int In_Np, double* In_BxArr, double* In_BzArr, char ArraysShouldBeAllocated)
	{
		sStart = In_sStart; sStep = In_sStep;
		Np = In_Np;
		if(ArraysShouldBeAllocated)
		{
			AllocateArrays();
			CopyArrays(In_BxArr, In_BzArr);
		}
		else
		{
            BxArr = In_BxArr; BzArr = In_BzArr;
			ArraysWereAllocated = 0;
		}
		gsStart = In_sStart;
        gsEnd = In_sStart + (In_Np - 1)*In_sStep;
	}
	srTMagFldTrUnif()
	{
		BxArr = BzArr = 0;
		Np = 0;
		sStart = sStep = 0;
		ArraysWereAllocated = 0;
		gsStart = gsEnd = 0.;
	}
	~srTMagFldTrUnif()
	{
		if(ArraysWereAllocated) DeleteArrays();
	}

	double GetLongExtent() //virtual
	{
		return sStep*(Np - 1);
	}
	srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat*); //virtual
	void SetupTrjDat(srTTrjDat*);

	srTMagFieldPeriodic* CreateAndSetupMagFieldPeriodic(double RelPrec, int MaxHarm, double MaxPerLen_m);
	srTMagFieldPeriodic* CreateAndSetupMagFieldPeriodicOld(double RelPrec, int MaxHarm, double MaxPerLen_m);
    void FindBasicFieldPeriodicParamAr(double* pB, int nB, double sInit, double sDelta, double absTolB, double& Per, double& L, double& sCen, double*& ar_sStartOnePer, int& nStartPer);
    void FindBasicFieldPeriodicParam(double* pB, int nB, double sInit, double sDelta, double absTolB, double& Per, double& L, double& sCen, double& sStartOnePer);
	void FindFieldHarmonicsAr(double* pB, int nB, double sInit, double sDelta, double Per, double* ar_sStartOnePer, int nPer, double RelPrec, char XorZ, int& NumHarm, srTMagHarm*& MarHarmArr);
	void FindFieldHarmonics(double* pB, int nB, double sInit, double sDelta, double Per, double sStartOnePer, double RelPrec, char XorZ, int& NumHarm, srTMagHarm*& MarHarmArr);
	double FindMaxAbsVal(double* Arr, int np);
    void FindFieldZeros(double* pB, int nB, double sStart, double sStep, double absTolB, double* ArgFldZerosIncr, double* ArgFldZerosDecr, int& AmOfZeros);
    void FindOnePeriodAr(double* ArgFldZeros, int AmOfZeros, double& Per, double* ar_sStartOnePer, int& nStartPer);
    void FindOnePeriod(double* ArgFldZeros, int AmOfZeros, double& sStartOnePer, double& Per);
    void InterpolateOnePeriodData(double* pB, int nB, double sInit, double sDelta, double sStartOnePer, double Per, double* InterpB, int AmOfInterpPts);
	void RotateOnePeriodData(double* InterpB, int AmOfInterpPts);
	void AnalyzeForHarmonics(double* pB, int AmOfPts, double Per, double RelPrec, char XorZ, int& AmOfHarm, srTMagHarm*& MarHarmArr);
	void AnalyzeForHarmonics_DeleteAuxArrays(float*& AuxDataContIn, float*& AuxDataContOut, double*& CkArr, double*& PhikArr, int*& HarmNoArr);
	void ChooseDominantBasicFieldPeriodicParamAr(
		double Per_HorFld, double L_HorFld, double sCen_HorFld, double* ar_sStartPer_HorFld, int nPer_HorFld, double MaxAbsHorFld,
		double Per_VertFld, double L_VertFld, double sCen_VertFld, double* ar_sStartPer_VertFld, int nPer_VertFld, double MaxAbsVertFld, 
		double& Per, double& L, double& sCen, double*& ar_sStartPer, int& nPer);
	void ChooseDominantBasicFieldPeriodicParam(
		double Per_HorFld, double L_HorFld, double sCen_HorFld, double sStartPer_HorFld, double MaxAbsHorFld,
		double Per_VertFld, double L_VertFld, double sCen_VertFld, double sStartPer_VertFld, double MaxAbsVertFld, 
		double& Per, double& L, double& sCen, double& sStartPer);

	void SumUpFieldHarmonics(srTMagHarm*& MagHarmArr_HorFld, int NumHarm_HorFld, srTMagHarm*& MagHarmArr_VertFld, int NumHarm_VertFld, srTMagHarm*& TotHarmArr, int& TotAmOfHarm);

	static srTMagFldTrUnif* SumUpSeveralFldTrUnif(srTMagFldCont* pMagTrUnifCont, srTMagFldCont* pMagContOpt);

	void ComputeSR_Stokes(srTEbmDat* pElecBeam, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes); //virtual

	double IntersectLineWithZero(double x1, double x2, double y1, double y2)
    {
        if(y1 == y2) return x1;
        return (x2*y1 - x1*y2)/(y1 - y2);
    }

	void AllocateArrays()
	{
        BxArr = new double[Np];
		if(BxArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		BzArr = new double[Np];
		if(BzArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		ArraysWereAllocated = 1;
	}
	void CopyArrays(double* In_BxArr, double* In_BzArr)
	{
		double *tBx = BxArr, *tBz = BzArr;
		double *tBxIn = In_BxArr, *tBzIn = In_BzArr;
		for(int i=0; i<Np; i++)
		{
			*(tBx++) = *(tBxIn++);
			*(tBz++) = *(tBzIn++);
		}
	}
	void DeleteArrays()
	{
		if(BxArr != 0) { delete[] BxArr; BxArr = 0;}
		if(BzArr != 0) { delete[] BzArr; BzArr = 0;}
		ArraysWereAllocated = 0;
	}
};

//*************************************************************************

class srTMagFieldPeriodic : public srTMagElem {
	
	//char FieldLongSym;

public:
// Use normal SI units in this structure: GeV, A, m, r !!!
	double PerLength;
	double TotLength;
	double Fund_keV_per_GeV2;
	int AmOfHarm;
	char TypeOfUnd; // 0- infinite, 1- normal, 2- tapered, 3- optical klystron
	double TaperPar_TU, PhaseSh_OK;
	double HalfKxE2pKzE2;

	srTMagHarmVect HarmVect;
	char FieldSymmetryInd; // 0- planar, 1- circular, 2- arbitrary
	char FldPolar; // 0, 'h', 'v'

	double sCen; // longitudinal coordinate of the undulator center

	//int NatFocTypeSASE; // 0- natural, 1- other
	double NatFocNxSASE, NatFocNySASE;
	int FldErrTypeSASE; // 0- No errors; 1- Uniform uncorrelated; 2- Uniform correlated; 3- Gaussian uncorrelated; 4- Gaussian correlated
	double FldErrRMS;
	int TaperTypeSASE; // 0- No taper; 1- Linear; 2- Quadratic
	double TaperStartSASE;
	double TaperRelFldChgSASE;

	srTMagFieldPeriodic(double Per, double L, double In_sCen, srTMagHarm* pHarm, int nHarm, char Type, double SpecPar);
	//srTMagFieldPeriodic(const SRWLMagFldU& inUnd, const TVector3d& inCenP); //SRWLIB
	srTMagFieldPeriodic(const SRWLMagFldU& inUnd, const TVector3d& inCenP, const TVector3d& inAxV, double inAng=0); //OC170615 //SRWLIB
	srTMagFieldPeriodic()
	{
		AmOfHarm = 0; TypeOfUnd = 1; 
		TaperPar_TU = PhaseSh_OK = 0.;
		HalfKxE2pKzE2 = 0.;
		sCen = 0.;

		TaperTypeSASE = 0; TaperStartSASE = 0; TaperRelFldChgSASE = 0;
		//NatFocTypeSASE = 0; NatFocNxSASE = 0;
		NatFocNxSASE = 0; NatFocNySASE = 0;
		FldErrTypeSASE = 0; FldErrRMS = 0;
		FldPolar = 0;
		gsStart = gsEnd = 0.;
	}
	~srTMagFieldPeriodic()
	{
		HarmVect.erase(HarmVect.begin(), HarmVect.end());
	}

	void OutBasicData(double* pPer, double* pL, double* psCen, int* pnHarm, char* pType, double* pSpecPar, double* pkeVperGeV2)
	{//used by DLL only
        *pPer = PerLength;
		*pL = TotLength;
		*psCen = sCen;
		*pnHarm = AmOfHarm;
        *pType = 'c'; //'c'- conventional, 't'- tapered, 'k'- optical klystron, 'i'- infinite
		*pSpecPar = 0.;
        if(TypeOfUnd == 0) *pType = 'i';
        else if(TypeOfUnd == 1) *pType = 'c';
        else if(TypeOfUnd == 2) { *pType = 't'; *pSpecPar = TaperPar_TU;}
		else if(TypeOfUnd == 3) { *pType = 'k'; *pSpecPar = PhaseSh_OK;}
		*pkeVperGeV2 = 0.; // To program eventually
	}
	void OutHarmData(int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase)
	{//used by DLL only
        int* tArrHarmNo = ArrHarmNo;
		char* tArrXorZ = ArrXorZ;
		double* tArrK = ArrK;
		double* tArrPhase = ArrPhase;
		for(int i=0; i<AmOfHarm; i++)
		{
			srTMagHarm& Harm = HarmVect[i];
			Harm.OutData(tArrHarmNo, tArrXorZ, tArrK, tArrPhase);
			tArrHarmNo++, tArrXorZ++, tArrK++, tArrPhase++;
		}
	}

	//int SetupFieldBasedArrays(srTEbmDat& EbmDat, int Ns, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtE2);
	int SetupFieldBasedArrays(srTEbmDat& EbmDat, long long Ns, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtE2);
	
	void AnalyzeFieldSymmetry()
	{
		char zFieldPresent = 0, xFieldPresent = 0;
        FldPolar = 0;
        for(int i=0; i<AmOfHarm; i++)
		{
			char XorZ = HarmVect[i].XorZ;
			if(XorZ == 'x') xFieldPresent = 1;
			if(XorZ == 'z') zFieldPresent = 1;
		}
		
		if(zFieldPresent && xFieldPresent) 
		{
			FieldSymmetryInd = 2; // arbitrary
			FldPolar = 0;
		}
		else 
		{
			if(zFieldPresent) FldPolar = 'v';
			else if(xFieldPresent) FldPolar = 'h';

			FieldSymmetryInd = 0; // planar
		}
	}
	void AnalyzeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ)
	{
		char zFieldPresent = 0, xFieldPresent = 0;
		for(int i=0; i<AmOfHarm; i++)
		{
			char XorZ = HarmVect[i].XorZ;
			if(XorZ == 'x') xFieldPresent = 1;
			if(XorZ == 'z') zFieldPresent = 1;
		}
		if(zFieldPresent && xFieldPresent) 
		{
			FieldSymmetryInd = 2; // arbitrary
		}
		else 
		{
			if(zFieldPresent) FldPolar = 'v';
			else if(xFieldPresent) FldPolar = 'h';

			FieldSymmetryInd = 0; // planar
		}

		FieldIsSymOverX = FieldIsSymOverZ = 0;
		if(!xFieldPresent) FieldIsSymOverZ = 1;
		if(!zFieldPresent) FieldIsSymOverX = 1;
	}

	void FindInitialTrjAngles(double& Btx0Gam, double& Btz0Gam, double& x0Gam, double& z0Gam)
	{
		Btx0Gam = Btz0Gam = x0Gam = z0Gam = 0.;
		for(int i=0; i<AmOfHarm; i++)
		{
			srTMagHarm& Harm = HarmVect[i];
			double K_d_n = Harm.K/Harm.HarmNo;
			double Buf = K_d_n*sin(Harm.Phase);
			double BufC = K_d_n*cos(Harm.Phase)/Harm.HarmNo;

			if(Harm.XorZ == 'z') 
			{
				Btx0Gam -= Buf; x0Gam += BufC;
			}
			else
			{
				Btz0Gam += Buf; z0Gam -= BufC;
			}
		}
		double lu_d_TwoPi = 0.1591549430919*PerLength;
		x0Gam *= lu_d_TwoPi; z0Gam *= lu_d_TwoPi;
	}

	void compB(TVector3d& inP, TVector3d& outB) //virtual
	{//this adds field to any previous value already in outB (as in Radia)
	 //"3/4 - 1/4" terminations are assumed
	 //Z is longitudinal coordinate
		
		//double xr = inP.x - mCenP.x, yr = inP.y - mCenP.y, zr = inP.z - mCenP.z;
		TVector3d Bloc = mTrans.TrVectField_inv(outB); //OC160615
		TVector3d Ploc = mTrans.TrPoint_inv(inP);
		double xr = Ploc.x, yr = Ploc.y, zr = Ploc.z;

		const double pi = 3.1415926535898;
		const double twoPi = 6.2831853072;
		const double piE2 = pi*pi;
		const double piE3 = piE2*pi;
		const double factK2B = 93.372904175766;

		double HalfTotLength = 0.5*TotLength, HalfPerLength = 0.5*PerLength;
		//double HalfTotLenWithTerm = HalfTotLength + PerLength;
		double HalfTotLenWithTerm = HalfTotLength + 4*PerLength;
		if((zr < -HalfTotLenWithTerm) || (zr > HalfTotLenWithTerm)) return;

		for(int i=0; i<AmOfHarm; i++)
		{
			srTMagHarm& Harm = HarmVect[i];
			double curPhase = Harm.Phase;
			//long nTwoPi = (long)(fabs(curPhase)/twoPi);
			long long nTwoPi = (long long)(fabs(curPhase)/twoPi);
			if(curPhase > 0) curPhase -= twoPi*nTwoPi;
			else if(curPhase < 0) curPhase += twoPi*nTwoPi;

			double PerLength_d_HarmNo = PerLength/Harm.HarmNo;
			double HalfPerLength_d_HarmNo = 0.5*PerLength_d_HarmNo;

			//double s0 = -curPhase*PerLength/(twoPi*Harm.HarmNo);
			double s0 = -curPhase*PerLength_d_HarmNo/twoPi;
			double zrr = zr - s0;
			if((zrr < -HalfTotLenWithTerm) || (zrr > HalfTotLenWithTerm)) continue;

			double UndWaveNum = twoPi/PerLength_d_HarmNo;
			//double Bm = (Harm.K/(PerLength_d_HarmNo*factK2B))*cosh(UndWaveNum*yr);
			double Bm = (Harm.K/(PerLength*factK2B))*cosh(UndWaveNum*yr); //OC261012  
			//OC261012: corrected since K is defined from harmonic field amplitude as: 93.37290417576577*PerLength*t_FldH->B
			//and it is used assuming this definition for firld/trajectory calculations in SRW for IGOR 

			double arg = UndWaveNum*zrr;
			double dB = 0; 

			double HalfTotLengthThisHarm = HalfTotLength;
			if(Harm.s == 1) HalfTotLengthThisHarm += 0.25*PerLength; //in case of symmetric field

			//if(zrr < -HalfTotLength)
			if(zrr < -HalfTotLengthThisHarm)
			{//1st termination
				//double d_zr = zrr + HalfTotLenWithTerm;
				//Bm *= ((d_zr < HalfPerLength)? 0.25 : 0.75);

				double perE2 = PerLength_d_HarmNo*PerLength_d_HarmNo;
				double p = 2.*piE2/(3.*perE2);
				//double zArg = zrr + HalfTotLength;
				double zArg = zrr + HalfTotLengthThisHarm;
				double zArgE2 = zArg*zArg;
				dB = (twoPi*Bm*zArg/PerLength_d_HarmNo)*(1. - 4.*piE2*zArgE2/(9.*perE2))*exp(-p*zArgE2);

				//double zTestSign = -HalfTotLength + 0.5*HalfPerLength_d_HarmNo;
				double zTestSign = -HalfTotLengthThisHarm + 0.5*HalfPerLength_d_HarmNo;
				//double testTrig = ((FieldLongSym == 1)? cos(UndWaveNum*zTestSign) : sin(UndWaveNum*zTestSign));
				double testTrig = ((Harm.s == 1)? cos(UndWaveNum*zTestSign) : sin(UndWaveNum*zTestSign));

				if(testTrig < 0) dB = -dB;
			}
			//else if(zrr > HalfTotLength)
			else if(zrr > HalfTotLengthThisHarm)
			{//2nd termination
				//double d_zr = zrr - HalfTotLength;
				//Bm *= ((d_zr < HalfPerLength)? 0.75 : 0.25);

				double perE2 = PerLength_d_HarmNo*PerLength_d_HarmNo;
				double p = 2.*piE2/(3.*perE2);
				//double zArg = zrr - HalfTotLength;
				double zArg = zrr - HalfTotLengthThisHarm;
				double zArgE2 = zArg*zArg;
				dB = (twoPi*Bm*zArg/PerLength_d_HarmNo)*(1. - 4.*piE2*zArgE2/(9.*perE2))*exp(-p*zArgE2);

				//double zTestSign = HalfTotLength - 0.5*HalfPerLength_d_HarmNo;
				double zTestSign = HalfTotLengthThisHarm - 0.5*HalfPerLength_d_HarmNo;
				//double testTrig = ((FieldLongSym == 1)? cos(UndWaveNum*zTestSign) : sin(UndWaveNum*zTestSign));
				double testTrig = ((Harm.s == 1)? cos(UndWaveNum*zTestSign) : sin(UndWaveNum*zTestSign));
				if(testTrig > 0) dB = -dB;
			}
			else
			{
				//dB = Bm*((FieldLongSym == 1)? cos(arg) : sin(arg)); //symmetric or anti-symmetric field
				dB = Bm*((Harm.s == 1)? cos(arg) : sin(arg)); //symmetric or anti-symmetric field
			}

			//if(Harm.XorZ == 'x') outB.x += dB;
			//else outB.y += dB;
			if(Harm.XorZ == 'x') Bloc.x += dB; //OC160615
			else Bloc.y += dB;
		}
		outB = mTrans.TrVectField(Bloc); //OC160615
	}

	srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat*); //virtual

	double GetLongExtent() { return TotLength;} //virtual
	void SetupWigSASE(srTWigComSASE&); //sets up SASE wiggler for Genesis
	void SetupExtMagFldU(SRWLMagFldU&, double&);

	void ComputeSR_Stokes(srTEbmDat* pElecBeam, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes); //virtual
	srTMagFldTrUnif* CreateAndSetupMagFldTrUnif();
};

//*************************************************************************

class srTMagFld3d : public srTMagElem {

	double *BxArr, *ByArr, *BzArr;
	int nx, ny, nz;
	double xStart, xStep, yStart, yStep, zStart, zStep;
	int nRep; //SRWLIB
	char ArraysWereAllocated;
	double xEnd, yEnd, zEnd;
	double *xArr, *yArr, *zArr; //SRWLIB
	int mInterp; //SRWLIB

	int m_nx_mi_2, m_ny_mi_2, m_nz_mi_2;

	map<pair<int, int>, CGenMathInterp*> mAuxSplineDataB;

public:

	//srTMagFld3d(double _xStart, double _xStep, int _nx, double _yStart, double _yStep, int _ny, double _zStart, double _zStep, int _nz, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated)
	srTMagFld3d(double _xStart, double _xStep, int _nx, double _yStart, double _yStep, int _ny, double _zStart, double _zStep, int _nz, double* _pBx, double* _pBy, double* _pBz, int _nRep, int _interp, char _arraysShouldBeAllocated)
	{
		SetupGrid(_xStart, _xStep, _nx, _yStart, _yStep, _ny, _zStart, _zStep, _nz, _pBx, _pBy, _pBz, _nRep, _arraysShouldBeAllocated);
		mInterp = _interp;
	}
	//srTMagFld3d(double _xStart, double _xStep, int _nx, double _yStart, double _yStep, int _ny, double _zStart, double _zStep, int _nz, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated, const TVector3d& inCenP) : srTMagElem(inCenP)
	//srTMagFld3d(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated, const TVector3d& inCenP) : srTMagElem(inCenP)
	//srTMagFld3d(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pX, double* _pY, double* _pZ, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated, const TVector3d& inCenP) : srTMagElem(inCenP)
	//srTMagFld3d(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pX, double* _pY, double* _pZ, double* _pBx, double* _pBy, double* _pBz, int _nRep, int _interp, char _arraysShouldBeAllocated, const TVector3d& inCenP) : srTMagElem(inCenP)
	srTMagFld3d(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pX, double* _pY, double* _pZ, double* _pBx, double* _pBy, double* _pBz, int _nRep, int _interp, char _arraysShouldBeAllocated, const TVector3d& inCenP, const TVector3d& inAxV, double inAng=0) : srTMagElem(inCenP, inAxV, inAng)
	{
		//SetupGridFromRange(_xRange, _nx, _yRange, _ny, _zRange, _nz, _pX, _pY, _pZ, _pBx, _pBy, _pBz, _nRep, _arraysShouldBeAllocated, inCenP);
		SetupGridFromRange(_xRange, _nx, _yRange, _ny, _zRange, _nz, _pX, _pY, _pZ, _pBx, _pBy, _pBz, _nRep, _arraysShouldBeAllocated);
		mInterp = _interp;
	}

	srTMagFld3d()
	{
		BxArr = ByArr = BzArr = 0;
		xArr = yArr = zArr = 0;
		nx = ny = nz = 0;
		m_nx_mi_2 = m_ny_mi_2 = m_nz_mi_2 = 0;
		mInterp = 1;

		xStart = xStep = yStart = yStep = zStart = zStep = 0;
		ArraysWereAllocated = 0;
		gsStart = gsEnd = 0.;
	}
	~srTMagFld3d()
	{
		if(ArraysWereAllocated) DeleteArrays();
		DeleteAuxSplineData();
		//DeallocAuxData(); //virtual in srTMagElem
	}

	void SetupGrid(double _xStart, double _xStep, int _nx, double _yStart, double _yStep, int _ny, double _zStart, double _zStep, int _nz, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated)
	{
		xArr = yArr = zArr = 0;
		xStart = _xStart; xStep = _xStep; nx = _nx;
		yStart = _yStart; yStep = _yStep; ny = _ny;
		zStart = _zStart; zStep = _zStep; nz = _nz;

		if(_nx <= 1) xStep = 0;
		if(_ny <= 1) yStep = 0;
		if(_nz <= 1) zStep = 0;

		m_nx_mi_2 = nx - 2; m_ny_mi_2 = ny - 2; m_nz_mi_2 = nz - 2;

		if(_arraysShouldBeAllocated)
		{
			AllocateArrays();
			CopyArrays(_pBx, _pBy, _pBz);
		}
		else
		{
            BxArr = _pBx; ByArr = _pBy; BzArr = _pBz;
			ArraysWereAllocated = 0;
		}
        //gsStart = _yStart;
		//gsEnd = _yStart + (_ny - 1)*_yStep;

		xEnd = _xStart + (_nx - 1)*_xStep; //SRWL
		yEnd = _yStart + (_ny - 1)*_yStep;
		zEnd = _zStart + (_nz - 1)*_zStep;
		//zEnd = _zStart + (_nz - 1)*_zStep*nRep;
		
		nRep = _nRep;

        gsStart = zStart;
		gsEnd = zEnd;
	}

	//void SetupGridFromRange(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pX, double* _pY, double* _pZ, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated, const TVector3d& inCenP)
	void SetupGridFromRange(double _xRange, int _nx, double _yRange, int _ny, double _zRange, int _nz, double* _pX, double* _pY, double* _pZ, double* _pBx, double* _pBy, double* _pBz, int _nRep, char _arraysShouldBeAllocated) //OC150815
	{
		xStart = -0.5*_xRange; xEnd = 0.5*_xRange; nx = _nx; xStep = (_nx <= 1)? 0 : _xRange/(_nx - 1);
		yStart = -0.5*_yRange; yEnd = 0.5*_yRange; ny = _ny; yStep = (_ny <= 1)? 0 : _yRange/(_ny - 1);
		zStart = -0.5*_zRange; zEnd = 0.5*_zRange; nz = _nz; zStep = (_nz <= 1)? 0 : _zRange/(_nz - 1);

		m_nx_mi_2 = nx - 2; m_ny_mi_2 = ny - 2; m_nz_mi_2 = nz - 2;

		xArr = yArr = zArr = 0;
		if(_arraysShouldBeAllocated)
		{
			AllocateArrays();
			CopyArrays(_pBx, _pBy, _pBz, _pX, _pY, _pZ);
		}
		else
		{
            BxArr = _pBx; ByArr = _pBy; BzArr = _pBz;
            xArr = _pX; yArr = _pY; zArr = _pZ;
			ArraysWereAllocated = 0;
		}

		if(_pX != 0) { xStart = _pX[0]; xEnd = _pX[_nx - 1]; xStep = (_nx <= 1)? 0 : (xEnd - xStart)/(_nx - 1);}
		if(_pY != 0) { yStart = _pY[0]; yEnd = _pY[_ny - 1]; yStep = (_ny <= 1)? 0 : (yEnd - yStart)/(_ny - 1);}
		if(_pZ != 0) { zStart = _pZ[0]; zEnd = _pZ[_nz - 1]; zStep = (_nz <= 1)? 0 : (zEnd - zStart)/(_nz - 1);}

		nRep = _nRep;

		gsStart = zStart;
		gsEnd = zEnd;

		//OC150815: commented-out (because CenP was moved to up to srTMagElem)
		//double halfTotLength = 0.5*_nRep*(zEnd - zStart); //OC01302011
		//gsStart = inCenP.z - halfTotLength;
		//gsEnd = inCenP.z + halfTotLength;
	}

	double GetLongExtent() //virtual
	{
		return yStep*(ny - 1);
	}

	srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat*); //virtual
    //void SetupTrjDat(srTTrjDat3d*);
	
	void compB(TVector3d& inP, TVector3d& outB) //virtual
	{//this adds field to any previous value already in outB (as in Radia)
		const double smallRelConst = 1.e-12;

		//double xr = inP.x - mCenP.x, yr = inP.y - mCenP.y, zr = inP.z - mCenP.z;
		TVector3d Bloc = mTrans.TrVectField_inv(outB); //OC160615
		TVector3d Ploc = mTrans.TrPoint_inv(inP);
		double xr = Ploc.x, yr = Ploc.y, zr = Ploc.z;

		if(((xr < xStart) || (xr >= xEnd)) && (xStart < xEnd)) return;
		if(((yr < yStart) || (yr >= yEnd)) && (yStart < yEnd)) return;
		//return immediately if point is outside the field definition regions

		if(nRep > 1)
		{
			double perLen = zEnd - zStart;
			if(perLen <= 0) return;
			int nRep_mi_1 = nRep - 1;
			zr += 0.5*perLen*nRep_mi_1;

			int iPer = (int)((zr - zStart)/perLen);
			if(iPer > nRep_mi_1) iPer = nRep_mi_1;
			if(iPer > 0)
			{
				zr -= iPer*perLen;
			}
		}
		if(((zr < zStart) || (zr >= zEnd)) && (zStart < zEnd)) return;
		//boundary point is only included at one side: [zStart, zEnd)
		//periods (nRep) are taken into account

		int ix = 0, ix1 = 0; 
		double x0 = xStart, xt = 0.;
		double xStepLocVar = xStep;
		if(nx > 1) 
		{
			if(nx < 3) ix = 0; 
			else 
			{
				ix = (int)((xr - xStart)/xStep + smallRelConst);
				if(ix < 0) ix = 0;
				else if(ix > m_nx_mi_2) ix = m_nx_mi_2;

				x0 = xStart + ix*xStep;
				if(xArr != 0)
				{
					if(xArr[ix] > xr)
					{
						int ii_start = ix - 1;
						ix = 0;
						for(int ii = ii_start; ii >= 0; ii--)
						{
							if(xArr[ii] <= xr) { ix = ii; break;}
						}
					}
					else
					{
						int ii_start = ix + 1;
						ix = m_nx_mi_2;
						for(int ii = ii_start; ii < nx; ii++)
						{
							if(xArr[ii] > xr) { ix = ii - 1; break;}
						}
						if(ix < 0) ix = 0;
						else if(ix > m_nx_mi_2) ix = m_nx_mi_2;
					}
					x0 = xArr[ix];
					xStepLocVar = xArr[ix + 1] - x0;
				}
			}
			ix1 = ix + 1;
			xt = (xr - x0)/xStepLocVar;
		}

		int iy = 0, iy1 = 0;
		double y0 = yStart, yt = 0.;
		double yStepLocVar = yStep;
		if(ny > 1) 
		{
			if(ny < 3) iy = 0; 
			else 
			{
				iy = (int)((yr - yStart)/yStep + smallRelConst);
				if(iy < 0) iy = 0;
				else if(iy > m_ny_mi_2) iy = m_ny_mi_2;

				y0 = yStart + iy*yStep;
				if(yArr != 0)
				{
					if(yArr[iy] > yr)
					{
						int ii_start = iy - 1;
						iy = 0;
						for(int ii = ii_start; ii >= 0; ii--)
						{
							if(yArr[ii] <= yr) { iy = ii; break;}
						}
					}
					else
					{
						int ii_start = iy + 1;
						iy = m_ny_mi_2;
						for(int ii = ii_start; ii < ny; ii++)
						{
							if(yArr[ii] > yr) { iy = ii - 1; break;}
						}
						if(iy < 0) iy = 0;
						else if(iy > m_ny_mi_2) iy = m_ny_mi_2;
					}
					y0 = yArr[iy];
					yStepLocVar = yArr[iy + 1] - y0;
				}
			}
			iy1 = iy + 1;
			yt = (yr - y0)/yStepLocVar;
		}

		int iz = 0, iz1 = 0;
		double z0 = zStart, zt = 0.;
		double zStepLocVar = zStep;
		if(nz > 1) 
		{
			if(nz < 3) iz = 0; 
			else 
			{
				iz = (int)((zr - zStart)/zStep + smallRelConst);
				if(iz < 0) iz = 0;
				else if(iz > m_nz_mi_2) iz = m_nz_mi_2;

				z0 = zStart + iz*zStep;
				if(zArr != 0)
				{
					if(zArr[iz] > zr)
					{
						int ii_start = iz - 1;
						iz = 0;
						for(int ii = ii_start; ii >= 0; ii--)
						{
							if(zArr[ii] <= zr) { iz = ii; break;}
						}
					}
					else
					{
						int ii_start = iz + 1;
						iz = m_nz_mi_2;
						for(int ii = ii_start; ii < nz; ii++)
						{
							if(zArr[ii] > zr) { iz = ii - 1; break;}
						}
						if(iz < 0) iz = 0;
						else if(iz > m_nz_mi_2) iz = m_nz_mi_2;
					}
					z0 = zArr[iz];
					zStepLocVar = zArr[iz + 1] - z0;
				}
			}
			iz1 = iz + 1;
			zt = (zr - z0)/zStepLocVar;
		}

		//long perY = nx;
		//long perZ = perY*ny;
		long long perY = nx;
		long long perZ = perY*ny;

		if(mInterp <= 1)
		{
			//long perX = ny;
			//long perZ = perX*nx;
			//long ofst000 = ix*perX + iy + iz*perZ, ofst100 = ix1*perX + iy + iz*perZ, ofst010 = ix*perX + iy1 + iz*perZ, ofst001 = ix*perX + iy + iz1*perZ;
			//long ofst110 = ix1*perX + iy1 + iz*perZ, ofst101 = ix1*perX + iy + iz1*perZ, ofst011 = ix*perX + iy1 + iz1*perZ, ofst111 = ix1*perX + iy1 + iz1*perZ;
			//inArrFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
			//long ofst000 = ix + iy*perY + iz*perZ;
			//long ofst100 = ix1 + iy*perY + iz*perZ;
			//long ofst010 = ix + iy1*perY + iz*perZ;
			//long ofst001 = ix + iy*perY + iz1*perZ;
			//long ofst110 = ix1 + iy1*perY + iz*perZ;
			//long ofst101 = ix1 + iy*perY + iz1*perZ;
			//long ofst011 = ix + iy1*perY + iz1*perZ;
			//long ofst111 = ix1 + iy1*perY + iz1*perZ;
			long long ofst000 = ix + iy*perY + iz*perZ;
			long long ofst100 = ix1 + iy*perY + iz*perZ;
			long long ofst010 = ix + iy1*perY + iz*perZ;
			long long ofst001 = ix + iy*perY + iz1*perZ;
			long long ofst110 = ix1 + iy1*perY + iz*perZ;
			long long ofst101 = ix1 + iy*perY + iz1*perZ;
			long long ofst011 = ix + iy1*perY + iz1*perZ;
			long long ofst111 = ix1 + iy1*perY + iz1*perZ;

			if(BxArr != 0)
			{
				double arrFunc[] = {BxArr[ofst000], BxArr[ofst100], BxArr[ofst010], BxArr[ofst001], BxArr[ofst110], BxArr[ofst101], BxArr[ofst011], BxArr[ofst111]};
				//outB.x += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc);
				Bloc.x += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc); //OC160615
			}
			if(ByArr != 0)
			{
				double arrFunc[] = {ByArr[ofst000], ByArr[ofst100], ByArr[ofst010], ByArr[ofst001], ByArr[ofst110], ByArr[ofst101], ByArr[ofst011], ByArr[ofst111]};
				//outB.y += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc);
				Bloc.y += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc); //OC160615
			}
			if(BzArr != 0)
			{
				double arrFunc[] = {BzArr[ofst000], BzArr[ofst100], BzArr[ofst010], BzArr[ofst001], BzArr[ofst110], BzArr[ofst101], BzArr[ofst011], BzArr[ofst111]};
				//outB.z += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc);
				Bloc.z += CGenMathInterp::Interp3dBilinRel(xt, yt, zt, arrFunc); //OC160615
			}
		}
		else if(mInterp == 2)
		{
			int ix0 = ix, iy0 = iy, iz0 = iz;
			if((xt >= 0.5) && (ix0 < m_nx_mi_2)) { ix0++; xt -= 1.; ix1++;}
			if((yt >= 0.5) && (iy0 < m_ny_mi_2)) { iy0++; yt -= 1.; iy1++;}
			if((zt >= 0.5) && (iz0 < m_nz_mi_2)) { iz0++; zt -= 1.; iz1++;}

			int ixm1 = ix0 - 1, iym1 = iy0 - 1, izm1 = iz0 - 1;
			if(ixm1 < 0) ixm1 = 0;
			if(iym1 < 0) iym1 = 0;
			if(izm1 < 0) izm1 = 0;

			//long ofst_00m1 = ix0 + iy0*perY + izm1*perZ;
			//long ofst_0m10 = ix0 + iym1*perY + iz0*perZ;
			//long ofst_m100 = ixm1 + iy0*perY + iz0*perZ;
			//long ofst_000 = ix0 + iy0*perY + iz0*perZ;
			//long ofst_100 = ix1 + iy0*perY + iz0*perZ;
			//long ofst_010 = ix0 + iy1*perY + iz0*perZ;
			//long ofst_110 = ix1 + iy1*perY + iz0*perZ;
			//long ofst_001 = ix0 + iy0*perY + iz1*perZ;
			//long ofst_101 = ix1 + iy0*perY + iz1*perZ;
			//long ofst_011 = ix0 + iy1*perY + iz1*perZ;

			long long ofst_00m1 = ix0 + iy0*perY + izm1*perZ;
			long long ofst_0m10 = ix0 + iym1*perY + iz0*perZ;
			long long ofst_m100 = ixm1 + iy0*perY + iz0*perZ;
			long long ofst_000 = ix0 + iy0*perY + iz0*perZ;
			long long ofst_100 = ix1 + iy0*perY + iz0*perZ;
			long long ofst_010 = ix0 + iy1*perY + iz0*perZ;
			long long ofst_110 = ix1 + iy1*perY + iz0*perZ;
			long long ofst_001 = ix0 + iy0*perY + iz1*perZ;
			long long ofst_101 = ix1 + iy0*perY + iz1*perZ;
			long long ofst_011 = ix0 + iy1*perY + iz1*perZ;

			if(BxArr != 0)
			{
				double arF[] = {BxArr[ofst_00m1], BxArr[ofst_0m10], BxArr[ofst_m100], BxArr[ofst_000], BxArr[ofst_100], BxArr[ofst_010], BxArr[ofst_110], BxArr[ofst_001], BxArr[ofst_101], BxArr[ofst_011]};
				//outB.x += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF);
				Bloc.x += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF); //OC160615
			}
			if(ByArr != 0)
			{
				double arF[] = {ByArr[ofst_00m1], ByArr[ofst_0m10], ByArr[ofst_m100], ByArr[ofst_000], ByArr[ofst_100], ByArr[ofst_010], ByArr[ofst_110], ByArr[ofst_001], ByArr[ofst_101], ByArr[ofst_011]};
				//outB.y += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF);
				Bloc.y += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF); //OC160615
			}
			if(BzArr != 0)
			{
				double arF[] = {BzArr[ofst_00m1], BzArr[ofst_0m10], BzArr[ofst_m100], BzArr[ofst_000], BzArr[ofst_100], BzArr[ofst_010], BzArr[ofst_110], BzArr[ofst_001], BzArr[ofst_101], BzArr[ofst_011]};
				//outB.z += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF);
				Bloc.z += CGenMathInterp::Interp3dQuadRel(xt, yt, zt, arF); //OC160615
			}
		}
		else if(mInterp == 3)
		{
			int ix0 = ix, iy0 = iy, iz0 = iz;
			int ixm1 = ix0 - 1, iym1 = iy0 - 1, izm1 = iz0 - 1;
			int ix2 = ix1 + 1, iy2 = iy1 + 1, iz2 = iz1 + 1;

			if(ixm1 < 0) 
			{
				ixm1 = 0;
				//if(nx > 3) { ix0++; ix1++; ix2++; xt -= 1.;}
			}
			if(iym1 < 0) 
			{
				iym1 = 0;
				//if(ny > 3) { iy0++; iy1++; iy2++; yt -= 1.;}
			}
			if(izm1 < 0) 
			{
				izm1 = 0;
				//if(nz > 3) { iz0++; iz1++; iz2++; zt -= 1.;}
			}

			if(ix2 >= nx) 
			{
				ix2 = ix1;
				//if(nx > 3) { ixm1--; ix0--; ix1--; xt += 1.;}
			}
			if(iy2 >= ny) 
			{
				iy2 = iy1;
				//if(ny > 3) { iym1--; iy0--; iy1--; yt += 1.;}
			}
			if(iz2 >= nz) 
			{
				iz2 = iz1;
				//if(nz > 3) { izm1--; iz0--; iz1--; zt += 1.;}
			}

			//long ofst_00m1 = ix0 + iy0*perY + izm1*perZ;
			//long ofst_10m1 = ix1 + iy0*perY + izm1*perZ;
			//long ofst_01m1 = ix0 + iy1*perY + izm1*perZ;
			//long ofst_11m1 = ix1 + iy1*perY + izm1*perZ;

			//long ofst_0m10 = ix0 + iym1*perY + iz0*perZ;
			//long ofst_1m10 = ix1 + iym1*perY + iz0*perZ;
			//long ofst_m100 = ixm1 + iy0*perY + iz0*perZ;
			//long ofst_000 = ix0 + iy0*perY + iz0*perZ;
			//long ofst_100 = ix1 + iy0*perY + iz0*perZ;
			//long ofst_200 = ix2 + iy0*perY + iz0*perZ;
			//long ofst_m110 = ixm1 + iy1*perY + iz0*perZ;
			//long ofst_010 = ix0 + iy1*perY + iz0*perZ;
			//long ofst_110 = ix1 + iy1*perY + iz0*perZ;
			//long ofst_210 = ix2 + iy1*perY + iz0*perZ;
			//long ofst_020 = ix0 + iy2*perY + iz0*perZ;
			//long ofst_120 = ix1 + iy2*perY + iz0*perZ;

			//long ofst_0m11 = ix0 + iym1*perY + iz1*perZ;
			//long ofst_1m11 = ix1 + iym1*perY + iz1*perZ;
			//long ofst_m101 = ixm1 + iy0*perY + iz1*perZ;
			//long ofst_001 = ix0 + iy0*perY + iz1*perZ;
			//long ofst_101 = ix1 + iy0*perY + iz1*perZ;
			//long ofst_201 = ix2 + iy0*perY + iz1*perZ;
			//long ofst_m111 = ixm1 + iy1*perY + iz1*perZ;
			//long ofst_011 = ix0 + iy1*perY + iz1*perZ;
			//long ofst_111 = ix1 + iy1*perY + iz1*perZ;
			//long ofst_211 = ix2 + iy1*perY + iz1*perZ;
			//long ofst_021 = ix0 + iy2*perY + iz1*perZ;
			//long ofst_121 = ix1 + iy2*perY + iz1*perZ;

			//long ofst_002 = ix0 + iy0*perY + iz2*perZ;
			//long ofst_102 = ix1 + iy0*perY + iz2*perZ;
			//long ofst_012 = ix0 + iy1*perY + iz2*perZ;
			//long ofst_112 = ix1 + iy1*perY + iz2*perZ;

			long long ofst_00m1 = ix0 + iy0*perY + izm1*perZ;
			long long ofst_10m1 = ix1 + iy0*perY + izm1*perZ;
			long long ofst_01m1 = ix0 + iy1*perY + izm1*perZ;
			long long ofst_11m1 = ix1 + iy1*perY + izm1*perZ;

			long long ofst_0m10 = ix0 + iym1*perY + iz0*perZ;
			long long ofst_1m10 = ix1 + iym1*perY + iz0*perZ;
			long long ofst_m100 = ixm1 + iy0*perY + iz0*perZ;
			long long ofst_000 = ix0 + iy0*perY + iz0*perZ;
			long long ofst_100 = ix1 + iy0*perY + iz0*perZ;
			long long ofst_200 = ix2 + iy0*perY + iz0*perZ;
			long long ofst_m110 = ixm1 + iy1*perY + iz0*perZ;
			long long ofst_010 = ix0 + iy1*perY + iz0*perZ;
			long long ofst_110 = ix1 + iy1*perY + iz0*perZ;
			long long ofst_210 = ix2 + iy1*perY + iz0*perZ;
			long long ofst_020 = ix0 + iy2*perY + iz0*perZ;
			long long ofst_120 = ix1 + iy2*perY + iz0*perZ;

			long long ofst_0m11 = ix0 + iym1*perY + iz1*perZ;
			long long ofst_1m11 = ix1 + iym1*perY + iz1*perZ;
			long long ofst_m101 = ixm1 + iy0*perY + iz1*perZ;
			long long ofst_001 = ix0 + iy0*perY + iz1*perZ;
			long long ofst_101 = ix1 + iy0*perY + iz1*perZ;
			long long ofst_201 = ix2 + iy0*perY + iz1*perZ;
			long long ofst_m111 = ixm1 + iy1*perY + iz1*perZ;
			long long ofst_011 = ix0 + iy1*perY + iz1*perZ;
			long long ofst_111 = ix1 + iy1*perY + iz1*perZ;
			long long ofst_211 = ix2 + iy1*perY + iz1*perZ;
			long long ofst_021 = ix0 + iy2*perY + iz1*perZ;
			long long ofst_121 = ix1 + iy2*perY + iz1*perZ;

			long long ofst_002 = ix0 + iy0*perY + iz2*perZ;
			long long ofst_102 = ix1 + iy0*perY + iz2*perZ;
			long long ofst_012 = ix0 + iy1*perY + iz2*perZ;
			long long ofst_112 = ix1 + iy1*perY + iz2*perZ;

			if(BxArr != 0)
			{
				double arF[] = {
					BxArr[ofst_00m1],BxArr[ofst_10m1],BxArr[ofst_01m1],BxArr[ofst_11m1],
					BxArr[ofst_0m10],BxArr[ofst_1m10],BxArr[ofst_m100],BxArr[ofst_000],BxArr[ofst_100],BxArr[ofst_200],BxArr[ofst_m110],BxArr[ofst_010],BxArr[ofst_110],BxArr[ofst_210],BxArr[ofst_020],BxArr[ofst_120],
					BxArr[ofst_0m11],BxArr[ofst_1m11],BxArr[ofst_m101],BxArr[ofst_001],BxArr[ofst_101],BxArr[ofst_201],BxArr[ofst_m111],BxArr[ofst_011],BxArr[ofst_111],BxArr[ofst_211],BxArr[ofst_021],BxArr[ofst_121],
					BxArr[ofst_002],BxArr[ofst_102],BxArr[ofst_012],BxArr[ofst_112]
				};
				//outB.x += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF);
				Bloc.x += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF); //OC160615
			}
			if(ByArr != 0)
			{
				double arF[] = {
					ByArr[ofst_00m1],ByArr[ofst_10m1],ByArr[ofst_01m1],ByArr[ofst_11m1],
					ByArr[ofst_0m10],ByArr[ofst_1m10],ByArr[ofst_m100],ByArr[ofst_000],ByArr[ofst_100],ByArr[ofst_200],ByArr[ofst_m110],ByArr[ofst_010],ByArr[ofst_110],ByArr[ofst_210],ByArr[ofst_020],ByArr[ofst_120],
					ByArr[ofst_0m11],ByArr[ofst_1m11],ByArr[ofst_m101],ByArr[ofst_001],ByArr[ofst_101],ByArr[ofst_201],ByArr[ofst_m111],ByArr[ofst_011],ByArr[ofst_111],ByArr[ofst_211],ByArr[ofst_021],ByArr[ofst_121],
					ByArr[ofst_002],ByArr[ofst_102],ByArr[ofst_012],ByArr[ofst_112]
				};
				//outB.y += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF);
				Bloc.y += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF); //OC160615
			}
			if(BzArr != 0)
			{
				double arF[] = {
					BzArr[ofst_00m1],BzArr[ofst_10m1],BzArr[ofst_01m1],BzArr[ofst_11m1],
					BzArr[ofst_0m10],BzArr[ofst_1m10],BzArr[ofst_m100],BzArr[ofst_000],BzArr[ofst_100],BzArr[ofst_200],BzArr[ofst_m110],BzArr[ofst_010],BzArr[ofst_110],BzArr[ofst_210],BzArr[ofst_020],BzArr[ofst_120],
					BzArr[ofst_0m11],BzArr[ofst_1m11],BzArr[ofst_m101],BzArr[ofst_001],BzArr[ofst_101],BzArr[ofst_201],BzArr[ofst_m111],BzArr[ofst_011],BzArr[ofst_111],BzArr[ofst_211],BzArr[ofst_021],BzArr[ofst_121],
					BzArr[ofst_002],BzArr[ofst_102],BzArr[ofst_012],BzArr[ofst_112]
				};
				//outB.z += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF);
				Bloc.z += CGenMathInterp::Interp3dBiCubic32pRel(xt, yt, zt, arF); //OC160615
			}
/**20-points version:
			long ofst_00m1 = ix0 + iy0*perY + izm1*perZ;
			long ofst_0m10 = ix0 + iym1*perY + iz0*perZ;
			long ofst_m100 = ixm1 + iy0*perY + iz0*perZ;
			long ofst_000 = ix0 + iy0*perY + iz0*perZ;
			long ofst_100 = ix1 + iy0*perY + iz0*perZ;
			long ofst_200 = ix2 + iy0*perY + iz0*perZ;
			long ofst_010 = ix0 + iy1*perY + iz0*perZ;
			long ofst_110 = ix1 + iy1*perY + iz0*perZ;
			long ofst_210 = ix2 + iy1*perY + iz0*perZ;
			long ofst_020 = ix0 + iy2*perY + iz0*perZ;
			long ofst_120 = ix1 + iy2*perY + iz0*perZ;
			long ofst_001 = ix0 + iy0*perY + iz1*perZ;
			long ofst_101 = ix1 + iy0*perY + iz1*perZ;
			long ofst_201 = ix2 + iy0*perY + iz1*perZ;
			long ofst_011 = ix0 + iy1*perY + iz1*perZ;
			long ofst_111 = ix1 + iy1*perY + iz1*perZ;
			long ofst_021 = ix0 + iy2*perY + iz1*perZ;
			long ofst_002 = ix0 + iy0*perY + iz2*perZ;
			long ofst_102 = ix1 + iy0*perY + iz2*perZ;
			long ofst_012 = ix0 + iy1*perY + iz2*perZ;
			if(BxArr != 0)
			{
				double arF[] = {BxArr[ofst_00m1],BxArr[ofst_0m10],BxArr[ofst_m100],BxArr[ofst_000],BxArr[ofst_100],BxArr[ofst_200],BxArr[ofst_010],BxArr[ofst_110],BxArr[ofst_210],BxArr[ofst_020],BxArr[ofst_120],
								BxArr[ofst_001],BxArr[ofst_101],BxArr[ofst_201],BxArr[ofst_011],BxArr[ofst_111],BxArr[ofst_021],BxArr[ofst_002],BxArr[ofst_102],BxArr[ofst_012]};
				//outB.x += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF);
				Bloc.x += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF); //OC160615
			}
			if(ByArr != 0)
			{
				double arF[] = {ByArr[ofst_00m1],ByArr[ofst_0m10],ByArr[ofst_m100],ByArr[ofst_000],ByArr[ofst_100],ByArr[ofst_200],ByArr[ofst_010],ByArr[ofst_110],ByArr[ofst_210],ByArr[ofst_020],ByArr[ofst_120],
								ByArr[ofst_001],ByArr[ofst_101],ByArr[ofst_201],ByArr[ofst_011],ByArr[ofst_111],ByArr[ofst_021],ByArr[ofst_002],ByArr[ofst_102],ByArr[ofst_012]};
				//outB.y += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF);
				Bloc.y += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF); //OC160615
			}
			if(BzArr != 0)
			{
				double arF[] = {BzArr[ofst_00m1],BzArr[ofst_0m10],BzArr[ofst_m100],BzArr[ofst_000],BzArr[ofst_100],BzArr[ofst_200],BzArr[ofst_010],BzArr[ofst_110],BzArr[ofst_210],BzArr[ofst_020],BzArr[ofst_120],
								BzArr[ofst_001],BzArr[ofst_101],BzArr[ofst_201],BzArr[ofst_011],BzArr[ofst_111],BzArr[ofst_021],BzArr[ofst_002],BzArr[ofst_102],BzArr[ofst_012]};
				//outB.z += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF);
				Bloc.z += CGenMathInterp::Interp3dCubicRel(xt, yt, zt, arF); //OC160615
			}
**/
		}
		else if(mInterp == 4)
		{
			double *arAuxBx_vs_Z=0, *arAuxBy_vs_Z=0, *arAuxBz_vs_Z=0; //for spline interpolations
			
			int ix0 = ix, iy0 = iy, iz0 = iz;
			int ixm1 = ix0 - 1, iym1 = iy0 - 1, izm1 = iz0 - 1;
			int ix2 = ix1 + 1, iy2 = iy1 + 1, iz2 = iz1 + 1;
			if(ixm1 < 0) ixm1 = 0;
			if(iym1 < 0) iym1 = 0;
			if(izm1 < 0) izm1 = 0;
			if(ix2 >= nx) ix2 = ix1;
			if(iy2 >= ny) iy2 = iy1;
			if(iz2 >= nz) iz2 = iz1;
			
			pair<int,int> arPairInd[] = {
				pair<int,int>(ix0,iym1), pair<int,int>(ix1,iym1),
				pair<int,int>(ixm1,iy0), pair<int,int>(ix0,iy0), pair<int,int>(ix1,iy0), pair<int,int>(ix2,iy0),
				pair<int,int>(ixm1,iy1), pair<int,int>(ix0,iy1), pair<int,int>(ix1,iy1), pair<int,int>(ix2,iy1),
				pair<int,int>(ix0,iy2), pair<int,int>(ix1,iy2)
			};

			double arCellBx[12], arCellBy[12], arCellBz[12];
			map<pair<int, int>, CGenMathInterp*>::const_iterator it;
			for(int i=0; i<12; i++)
			{
				pair<int,int> *pCurPairInd = arPairInd + i;
				//long ofst = pCurPairInd->first + (pCurPairInd->second)*perY + izm1*perZ;
				//long ofst0 = pCurPairInd->first + (pCurPairInd->second)*perY;
				long long ofst0 = pCurPairInd->first + (pCurPairInd->second)*perY;

				CGenMathInterp *curSplineDataB = 0;
				it = mAuxSplineDataB.find(*pCurPairInd);
				if(it == mAuxSplineDataB.end())
				{
					curSplineDataB = new CGenMathInterp[3];
					if(BxArr != 0)
					{
						if(arAuxBx_vs_Z == 0) arAuxBx_vs_Z = new double[nz];
						double *tBx = arAuxBx_vs_Z, *tBxOrig = BxArr + ofst0;
						for(int iz=0; iz<nz; iz++) { *(tBx++) = *tBxOrig; tBxOrig += perZ;}

						if(zArr == 0) curSplineDataB->InitCubicSplineU(zStart, zStep, arAuxBx_vs_Z, nz);
						else curSplineDataB->InitCubicSpline(zArr, arAuxBx_vs_Z, nz);
					}
					if(ByArr != 0)
					{
						if(arAuxBy_vs_Z == 0) arAuxBy_vs_Z = new double[nz];
						double *tBy = arAuxBy_vs_Z, *tByOrig = ByArr + ofst0;
						for(int iz=0; iz<nz; iz++) { *(tBy++) = *tByOrig; tByOrig += perZ;}

						if(zArr == 0) (curSplineDataB + 1)->InitCubicSplineU(zStart, zStep, arAuxBy_vs_Z, nz);
						else (curSplineDataB + 1)->InitCubicSpline(zArr, arAuxBy_vs_Z, nz);
					}
					if(BzArr != 0)
					{
						if(arAuxBz_vs_Z == 0) arAuxBz_vs_Z = new double[nz];
						double *tBz = arAuxBz_vs_Z, *tBzOrig = BzArr + ofst0;
						for(int iz=0; iz<nz; iz++) { *(tBz++) = *tBzOrig; tBzOrig += perZ;}

						if(zArr == 0) (curSplineDataB + 2)->InitCubicSplineU(zStart, zStep, arAuxBz_vs_Z, nz);
						else (curSplineDataB + 2)->InitCubicSpline(zArr, arAuxBy_vs_Z, nz);
					}
					mAuxSplineDataB[*pCurPairInd] = curSplineDataB;
				}
				else curSplineDataB = it->second;

				if(BxArr != 0)
				{
					arCellBx[i] = (zArr == 0)? curSplineDataB->InterpRelCubicSplineU(zt, iz0) : curSplineDataB->InterpRelCubicSpline(zt, iz0, zStepLocVar);
				}
				if(ByArr != 0)
				{
					arCellBy[i] = (zArr == 0)? (curSplineDataB + 1)->InterpRelCubicSplineU(zt, iz0) : (curSplineDataB + 1)->InterpRelCubicSpline(zt, iz0, zStepLocVar);
				}
				if(BzArr != 0)
				{
					arCellBz[i] = (zArr == 0)? (curSplineDataB + 2)->InterpRelCubicSplineU(zt, iz0) : (curSplineDataB + 2)->InterpRelCubicSpline(zt, iz0, zStepLocVar);
				}
			}

			//use 2D cubic interpolation (based on 12 points) to find actual B:
			//if(BxArr != 0) outB.x += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBx);
			//if(ByArr != 0) outB.y += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBy);
			//if(BzArr != 0) outB.z += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBz);
			if(BxArr != 0) Bloc.x += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBx); //OC160615
			if(ByArr != 0) Bloc.y += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBy);
			if(BzArr != 0) Bloc.z += CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arCellBz);

			if(arAuxBx_vs_Z != 0) { delete[] arAuxBx_vs_Z; arAuxBx_vs_Z=0;}
			if(arAuxBy_vs_Z != 0) { delete[] arAuxBy_vs_Z; arAuxBy_vs_Z=0;}
			if(arAuxBz_vs_Z != 0) { delete[] arAuxBz_vs_Z; arAuxBz_vs_Z=0;}
		}
		outB = mTrans.TrVectField(Bloc); //OC160615
	}

	void tabulateB(srTMagElem* pMagElem)
	{
		if(pMagElem == 0) return;

		TVector3d vP, vB;
		double *tBx = BxArr, *tBy = ByArr, *tBz = BzArr;
		double z = zStart; //+ mCenP.z; //OC160615
		for(int iz=0; iz<nz; iz++)
		{
			if(zArr != 0) z = zArr[iz]; //+ mCenP.z; //OC160615
			//vP.z = z;
			double y = yStart; //+ mCenP.y; //OC160615
			for(int iy=0; iy<ny; iy++)
			{
				if(yArr != 0) y = yArr[iy]; //+ mCenP.y; //OC160615
				//vP.y = y;
				double x = xStart; //+ mCenP.x; //OC160615
				for(int ix=0; ix<nx; ix++)
				{
					if(xArr != 0) x = xArr[ix]; //+ mCenP.x; //OC160615
					
					vP.x = x; vP.y = y; vP.z = z; //OC160615
					vP = mTrans.TrPoint(vP);
					vB.x = vB.y = vB.z = 0.;
					pMagElem->compB(vP, vB);
					vB = mTrans.TrVectField_inv(vB); //OC170615??

					if(BxArr != 0) *(tBx++) = vB.x;
					if(ByArr != 0) *(tBy++) = vB.y;
					if(BzArr != 0) *(tBz++) = vB.z;

					x += xStep;
				}
				y += yStep;
			}
			z += zStep;
		}
	}

	void tabInterpB(srTMagFldCont& magCont, double arPrecPar[6], double* arPar1, double* arPar2, double* arCoefBx, double* arCoefBy); //OC02112017
	//void tabInterpB(srTMagFldCont& magCont, double* arPrecPar, double* arPar1, double* arPar2, double* arCoefBx, double* arCoefBy);

	//void DeallocAuxData() //virtual in srTMagElem
	//{
	//	if(ArraysWereAllocated) DeleteArrays();
	//	DeleteAuxSplineData();
	//}

	void AllocateArrays()
	{
		//long Np = nx*ny*nz;
		long long Np = ((long long)nx)*((long long)ny)*((long long)nz);
        BxArr = new double[Np];
		if(BxArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		ByArr = new double[Np];
		if(ByArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		BzArr = new double[Np];
		if(BzArr == 0) throw MEMORY_ALLOCATION_FAILURE;

        xArr = new double[nx];
		if(xArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		yArr = new double[ny];
		if(yArr == 0) throw MEMORY_ALLOCATION_FAILURE;
		zArr = new double[nz];
		if(zArr == 0) throw MEMORY_ALLOCATION_FAILURE;

		ArraysWereAllocated = 1;
	}
	void CopyArrays(double* In_BxArr, double* In_ByArr, double* In_BzArr, double* In_xArr=0, double* In_yArr=0, double* In_zArr=0)
	{
		double *tBx = BxArr, *tBy = ByArr, *tBz = BzArr, *tx = xArr, *ty = yArr, *tz = zArr;
		double *tBxIn = In_BxArr, *tByIn = In_ByArr, *tBzIn = In_BzArr, *txIn = In_xArr, *tyIn = In_yArr, *tzIn = In_zArr;
		for(int iz=0; iz<nz; iz++)
		{
			for(int ix=0; ix<nx; ix++)
			{
				for(int iy=0; iy<ny; iy++)
				{
					if(In_BxArr != 0) *(tBx++) = *(tBxIn++);
					if(In_ByArr != 0) *(tBy++) = *(tByIn++);
					if(In_BzArr != 0) *(tBz++) = *(tBzIn++);
				}
			}
			if(In_zArr != 0) *(tz++) = *(tzIn++);
		}
		if(In_xArr != 0)
		{
			for(int ix=0; ix<nx; ix++) *(tx++) = *(txIn++);
		}
		if(In_yArr != 0)
		{
			for(int iy=0; iy<ny; iy++) *(ty++) = *(tyIn++);
		}
	}
	void DeleteArrays()
	{
		if(BxArr != 0) { delete[] BxArr; BxArr = 0;}
		if(ByArr != 0) { delete[] ByArr; ByArr = 0;}
		if(BzArr != 0) { delete[] BzArr; BzArr = 0;}

		if(xArr != 0) { delete[] xArr; xArr = 0;}
		if(yArr != 0) { delete[] yArr; yArr = 0;}
		if(zArr != 0) { delete[] zArr; zArr = 0;}

		ArraysWereAllocated = 0;
	}
	void DeleteAuxSplineData()
	{
		if(mAuxSplineDataB.empty()) return;

		for(map<pair<int, int>, CGenMathInterp*>::iterator it = mAuxSplineDataB.begin(); it != mAuxSplineDataB.end(); ++it)
		{
			delete[] it->second;
			it->second = 0;
		}
		mAuxSplineDataB.erase(mAuxSplineDataB.begin(), mAuxSplineDataB.end());
	}
};

//*************************************************************************

//class srTMagQuad : public srTMagElem {
class srTMagMult : public srTMagElem {
//Multipole magnet, replaced former srTMagQuad (Quadrupoole)

	double m_HalfLenModConst, m_LenModEdge; // Constants for analytically integrable soft-edge quadrupole model
	//validation is required!

	double Strength; // [T] for dipole, [T/m] for quadrupole (negative means defocusing for x), [T/m^2] for sextupole, [T/m^3] for octupole
	char m; // Order: 1 for dipole, 2 for quadrupole, 3 for sextupole, 4 for octupole
	char n_or_s; // Normal ('n') or skew ('s')
	double Length; // Effective length [m]
	double LenEdge; // Edge length for field variation from 10% to 90% [m]
	double m_R; // Radius of curvature of central trajectory [m] (for simulating e.g. quadrupole component integrated to a bending magnet; effective if > 0)

public:
	
	//TVector2d TransvCenPoint;
	//double sCen; // longitudinal coordinate of the center

	//srTMagQuad(double InLength=0., double InStrength=0., double In_sCen=0.)
	//srTMagQuad(double InStrength=0., double InLength=0., double InLenEdge=0., double In_sCen=0.)
	srTMagMult(double InStrength=0., char In_m=2, double InLength=0., double InLenEdge=0., double In_sCen=0.)
	{
		//TransvCenPoint.x = TransvCenPoint.y = 0.;
		//sCen = In_sCen;
		//mCenP = TVector3d(0, 0, In_sCen);
		TVector3d cenP(0,0,In_sCen), axV(0,0,1); //OC160615
		SetupOrient(cenP, axV);

		Length = InLength; Strength = InStrength; LenEdge = InLenEdge;
		gsStart = In_sCen - 0.5*Length;
		gsEnd = In_sCen + 0.5*Length;

		m_LenModEdge = InLenEdge/1.23789045853;
		m_HalfLenModConst = 0.5*(InLength - 1.2689299897*InLenEdge); //validation is required!
		
		m = In_m;
		n_or_s = 'n';
		m_R = 0;
	}
	//srTMagQuad(double InStrength, char In_n_or_s, double InLength, double InLenEdge, const TVector3d& inCenP) : srTMagElem(inCenP)
	//srTMagMult(double InStrength, char In_m, char In_n_or_s, double InLength, double InLenEdge, const TVector3d& inCenP, const TVector3d& inAxV, double inAng=0) : srTMagElem(inCenP, inAxV, inAng)
	srTMagMult(double InStrength, char In_m, char In_n_or_s, double InLength, double InLenEdge, double InR, const TVector3d& inCenP, const TVector3d& inAxV, double inAng=0) : srTMagElem(inCenP, inAxV, inAng)
	{
		Length = InLength; 
		m = In_m;
		Strength = InStrength; 
		LenEdge = InLenEdge;
		n_or_s = In_n_or_s;

		m_R = InR; //OC310715

		//TransvCenPoint.x = inCenP.x;
		//TransvCenPoint.y = inCenP.y;
		//sCen = inCenP.z;
		//gsStart = inCenP.z - 0.5*Length;
		//gsEnd = inCenP.z + 0.5*Length;
		gsStart = -0.5*Length; //OC170615
		gsEnd = 0.5*Length;

		m_LenModEdge = InLenEdge/1.23789045853;
		m_HalfLenModConst = 0.5*(InLength - 1.2689299897*InLenEdge); //validation is required!
	}
	//srTMagQuad(srTStringVect* pElemInfo)
	srTMagMult(srTStringVect* pElemInfo)
	{
		char* ElemID = (*pElemInfo)[0];
		if((!strcmp(ElemID, "Dipole")) || (!strcmp(ElemID, "dipole")) || (!strcmp(ElemID, "DIPOLE"))) m = 1;
		else if((!strcmp(ElemID, "Quadrupole")) || (!strcmp(ElemID, "quadrupole")) || (!strcmp(ElemID, "QUADRUPOLE"))) m = 2;
		else if((!strcmp(ElemID, "Sextupole")) || (!strcmp(ElemID, "sextupole")) || (!strcmp(ElemID, "SEXTUPOLE"))) m = 3;
		else if((!strcmp(ElemID, "Octupole")) || (!strcmp(ElemID, "octupole")) || (!strcmp(ElemID, "OCTUPOLE"))) m = 4;

		Strength = atof((*pElemInfo)[1]); // [T/m]
		Length = atof((*pElemInfo)[2]); // [m]

		m_R = 0; // [m]

		//sCen = 0;
		//mCenP = TVector3d(0,0,0);

		TVector3d cenP(0,0,0), vZero(0,0,0); //OC170615
		if(pElemInfo->size() > 3)
		{
			//TransvCenPoint.x = atof((*pElemInfo)[3]); // [m]
			//TransvCenPoint.y = atof((*pElemInfo)[4]); // [m]
			//mCenP.x = atof((*pElemInfo)[3]); // [m]
			//mCenP.y = atof((*pElemInfo)[4]); // [m]
			cenP.x = atof((*pElemInfo)[3]); // [m] //OC170615
			cenP.y = atof((*pElemInfo)[4]); // [m]
		}

		SetupOrient(cenP, vZero); //OC170615

		//gsStart = mCenP.z - 0.5*Length;
		//gsEnd = mCenP.z + 0.5*Length;
		gsStart = -0.5*Length; //?? OC170615 
		gsEnd = 0.5*Length;
	}

	void SetupWigSASE(srTWigComSASE&); //sets up SASE wiggler for Genesis

	void ComputeParticlePropagMatrix(double s, TMatrix2d& Mx, TMatrix2d& Mz) // virtual
	{
		//to implement !!!
	}

	void compB(TVector3d& inP, TVector3d& outB) //virtual, used by SRWLIB
	{//this adds field to any previous value already in outB (as in Radia)
		//Z is longitudinal coord. here
		//double xr = inP.x - mCenP.x, yr = inP.y - mCenP.y, zr = inP.z - mCenP.z;
		TVector3d Bloc = mTrans.TrVectField_inv(outB); //OC170615
		TVector3d Ploc = mTrans.TrPoint_inv(inP);
		double xr = Ploc.x, yr = Ploc.y, zr = Ploc.z;

		if(m_R != 0)
		{	//OC310715
			//double Rmix0 = m_R - xr;
			//xr = m_R - sqrt(Rmix0*Rmix0 + zr*zr);

			//OC23102016
			if(m_R < 0)
			{
				double XmiR = xr - m_R;
				xr = sqrt(XmiR*XmiR + zr*zr) + m_R;
				double alp = atan(zr/XmiR);
				zr = -alp*m_R;
			}
			else
			{
				double RmiX = m_R - xr;
				xr = m_R - sqrt(RmiX*RmiX + zr*zr);
				double alp = atan(zr/RmiX);
				zr = alp*m_R;
			}
			//double RpX = m_R + xr;
			//xr = sqrt(RpX*RpX + zr*zr) - m_R;
			//double alp = atan(zr/RpX);
			//zr = m_R*alp;
		}

		double Strength_u = 0.;

		if(m_LenModEdge <= 0)
		{//hard-edge case
			if((zr < -m_HalfLenModConst) || (zr > m_HalfLenModConst)) return;
			Strength_u = Strength;
		}
		else
		{
			double extEdges = 15*m_LenModEdge;
			if((zr < -m_HalfLenModConst - extEdges) || (zr > m_HalfLenModConst + extEdges)) return;
			
			Strength_u = Strength;
			if(zr < -m_HalfLenModConst)
			{
				double r = (zr + m_HalfLenModConst)/m_LenModEdge;
				double sqr_u = 1./(1. + r*r);
				Strength_u *= sqr_u*sqr_u;
			}
			else if(zr > m_HalfLenModConst)
			{
				double r = (zr - m_HalfLenModConst)/m_LenModEdge;
				double sqr_u = 1./(1. + r*r);
				Strength_u *= sqr_u*sqr_u;
			}
		}

		if(n_or_s == 'n')
		{//normal (to check!)
			if(m == 1)
			{//dipole
				//outB.x -= Strength_u;
				//outB.y -= Strength_u; //Strength = -By = B1
				Bloc.y -= Strength_u; //OC170615 //Strength = -By = B1
			}
			else if(m == 2)
			{//quadrupole
				//outB.x -= Strength_u*yr; //Strength = -dBx/dy = -dBy/dx = 2*B2
				//outB.y -= Strength_u*xr;
				Bloc.x -= Strength_u*yr; //OC170715 //Strength = -dBx/dy = -dBy/dx = 2*B2
				Bloc.y -= Strength_u*xr;
			}
			else if(m == 3)
			{//sextupole
				//outB.x -= Strength_u*xr*yr; //Strength = -d2By/dx2 = d2By/dy2 = 6*B3
				//outB.y -= 0.5*Strength_u*(xr*xr - yr*yr);
				Bloc.x -= Strength_u*xr*yr; //OC170615 //Strength = -d2By/dx2 = d2By/dy2 = 6*B3
				Bloc.y -= 0.5*Strength_u*(xr*xr - yr*yr);
			}
			else if(m == 4)
			{//octupole
				//outB.x += Strength_u*(-0.5*xr*xr*yr + yr*yr*yr/6.); //Strength = -d3By/dx3 = d3Bx/dy3 = 24*B4
				//outB.y += Strength_u*(-xr*xr*xr/6. + 0.5*xr*yr*yr);
				Bloc.x += Strength_u*(-0.5*xr*xr*yr + yr*yr*yr/6.); //OC170615 //Strength = -d3By/dx3 = d3Bx/dy3 = 24*B4
				Bloc.y += Strength_u*(-xr*xr*xr/6. + 0.5*xr*yr*yr);
			}
		}
		else if(n_or_s == 's')
		{//skew (to check!)
			if(m == 1)
			{//dipole
				//outB.x += Strength_u; //Strength = Bx = -A1
				Bloc.x += Strength_u; //OC170615 //Strength = Bx = -A1
				//outB.y -= Strength_u;
			}
			else if(m == 2)
			{//quadrupole
				//outB.x += Strength_u*xr; //Strength = dBx/dx = -dBy/dy = -2*A2
				//outB.y -= Strength_u*yr;
				Bloc.x += Strength_u*xr; //OC170615 //Strength = dBx/dx = -dBy/dy = -2*A2
				Bloc.y -= Strength_u*yr;
			}
			else if(m == 3)
			{//sextupole
				//outB.x += 0.5*Strength_u*(xr*xr - yr*yr); //Strength = d2Bx/dx2 = -d2Bx/dy2 = -6*A3
				//outB.y -= Strength_u*xr*yr;
				Bloc.x += 0.5*Strength_u*(xr*xr - yr*yr); //OC170615 //Strength = d2Bx/dx2 = -d2Bx/dy2 = -6*A3
				Bloc.y -= Strength_u*xr*yr;
			}
			else if(m == 4)
			{//octupole
				//outB.x += Strength_u*(xr*xr*xr/6. - 0.5*xr*yr*yr); //Strength = d3Bx/dx3 = d3By/dy3 = -24*A4
				//outB.y += Strength_u*(yr*yr*yr/6. - 0.5*xr*xr*yr);
				Bloc.x += Strength_u*(xr*xr*xr/6. - 0.5*xr*yr*yr); //OC170615 //Strength = d3Bx/dx3 = d3By/dy3 = -24*A4
				Bloc.y += Strength_u*(yr*yr*yr/6. - 0.5*xr*xr*yr);
			}
		}
		outB = mTrans.TrVectField(Bloc); //OC170615
	}
};

//*************************************************************************

class srTMagSol : public srTMagElem {

	double m_HalfLenModConst, m_LenModEdge; // Constants for analytically inttegratable soft-edge quadrupole model

	double B; //Magnetic field [T]
	double Length; //Effective length [m]
	double LenEdge; //Edge length for field variation from 10% to 90% [m]

public:

	//srTMagSol(double InB, double InLength, double InLenEdge, const TVector3d& inCenP) : srTMagElem(inCenP)
	srTMagSol(double InB, double InLength, double InLenEdge, const TVector3d& inCenP, const TVector3d& inAxV, double inAng=0) : srTMagElem(inCenP, inAxV, inAng)
	{
		Length = InLength; B = InB; LenEdge = InLenEdge;

		gsStart = inCenP.z - 0.5*Length;
		gsEnd = inCenP.z + 0.5*Length;

		m_HalfLenModConst = 0.5*InLength;
		m_LenModEdge = InLenEdge;
	}

	void compB(TVector3d& inP, TVector3d& outB) //virtual, used by SRWLIB
	{//this adds field to any previous value already in outB (as in Radia)
		//Z is longitudinal coord. here
		//double zr = inP.z - mCenP.z;
		TVector3d Bloc = mTrans.TrVectField_inv(outB); //OC170615
		TVector3d Ploc = mTrans.TrPoint_inv(inP);
		//double xr = Ploc.x, yr = Ploc.y, zr = Ploc.z;
		double zr = Ploc.z;

		if((zr > -m_HalfLenModConst) && (zr < m_HalfLenModConst))
		{
			//outB.z += B;	
			Bloc.z += B;	
		}
		outB = mTrans.TrVectField(Bloc); //OC170615
	}
};

//*************************************************************************

class srTMagDrift : public srTMagElem {

	double Length; // Effective length, [m]

public:
	
	srTMagDrift(double InLength=0., double In_sStart=0.)
	{
		Length = InLength;
		gsStart = In_sStart;
        gsEnd = In_sStart + InLength;
	}

	void ComputeParticlePropagMatrix(double s, TMatrix2d& Mx, TMatrix2d& Mz) // virtual
	{
		//to implement !!!
		//take into account partial propagation, depending on s !
	}
};

//*************************************************************************

class srTMagChicane : public srTMagElem {
public:

	double m_DipoleB, m_DipoleLength;
	int m_DipoleNum;
	double m_TotDriftLength;
	double m_sCen; // longitudinal coordinate of the center

	srTMagChicane(double DipoleB=0., double DipoleLength=0., int DipoleNum=0, double TotDriftLength=0., double sCen=0.)
	{
		m_DipoleB = DipoleB;
		m_DipoleLength = DipoleLength;
		m_DipoleNum = DipoleNum;
		m_TotDriftLength = TotDriftLength;
		m_sCen = sCen;

		double TotLength = CalcTotLength();
		gsStart = m_sCen - 0.5*TotLength;
		gsEnd = m_sCen + 0.5*TotLength;
	}

	srTMagChicane(srTStringVect* pElemInfo)
	{
		m_DipoleB = atof((*pElemInfo)[1]);
		m_DipoleLength = atof((*pElemInfo)[2]);
		m_DipoleNum = atoi((*pElemInfo)[3]);
		m_TotDriftLength = atof((*pElemInfo)[4]);
		m_sCen = 0;

		double TotLength = CalcTotLength();
		gsStart = m_sCen - 0.5*TotLength;
		gsEnd = m_sCen + 0.5*TotLength;
	}

	double CalcTotLength()
	{
		return m_TotDriftLength + m_DipoleNum*m_DipoleLength;
	}

	void SetupWigSASE(srTWigComSASE& WigCom)//sets up SASE wiggler for Genesis
	{
		WigCom.chic_bfield = m_DipoleB;
		WigCom.chic_magl = m_DipoleLength;
		WigCom.chic_dril = m_TotDriftLength;
	}

	void ComputeParticlePropagMatrix(double s, TMatrix2d& Mx, TMatrix2d& Mz) // virtual
	{
		//to implement !!!
	}
};

//*************************************************************************

typedef CSmartPtr<srTMagElem> CHMagFld;

//*************************************************************************

struct srTMagPosAndElem {
	double s;
	CHMagFld MagHndl;
};

typedef vector<srTMagPosAndElem, allocator<srTMagPosAndElem> > srTMagPosAndElemVect;

//*************************************************************************

class srTMagElemSummary {
public:
	int SetupMagElement(srTStringVect*, CHMagFld&);
};

//*************************************************************************

class srTMagGroup : public srTMagElem {
public:
	srTMagPosAndElemVect PosAndElemVect; // vector of srTMagPosAndElem

	srTMagGroup() {}

	void SetupWigSASE(srTWigComSASE&);
};

//*************************************************************************

#endif
