/************************************************************************//**
 * File: srthckbm.h
 * Description: SR Stokes parameters calculation method for ~Arbitrary magnetic field (used rarely, under-programmed) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRTHCKBM_H
#define __SRTHCKBM_H

#include "srgtrjdt.h"
#include "srstraux.h"
#include "srinterf.h"

//*************************************************************************

struct srTParPrecStokesArb;
class srTMagFldTrUnif;
class srTMagFldCont;
class srTTrjDat;
class srTStokesStructAccessData;
class CGenMathRand;

//*************************************************************************

struct srTExternIntKey2D {
	char k00, k10, k20, k21, k22;
	srTExternIntKey2D(char In_k00=0, char In_k10=0, char In_k20=0, char In_k21=0, char In_k22=0)
	{
        k00 = In_k00; k10 = In_k10; k20 = In_k20; k21 = In_k21; k22 = In_k22;
	}
};

//*************************************************************************

struct srTRadIntThickBeamAuxParams {
	double PI;
	double Half_k_d_e, k_d_e;

	double xc, xpc, zc, zpc;
	double Fc, Fn, C1, Cm;
	double BetX, BetZ, AlpX, AlpZ, GamX, GamZ, ThetXZ, ThetXpZ, ThetXZp, ThetXpZp, Bgam, GammaEm2;
	double InvBgam, HalfInvBgam;

	TComplexD Aux_deltz1000, Aux_alpz1000, Aux_deltx1000, Aux_alpx1000;

	srTEXZY EXZY;
	double xObsE2, zObsE2, xzObs;

	void Setup(srTEbmDat& ElecBeam);
};

//*************************************************************************

class srTRadIntThickBeam {

	srTFieldBasedArrays gFldArr;
	srTRadIntThickBeamAuxParams gAuxPar;
	//srTEXZY gEXZY;

	TComplexD *gCoefA, *gCoefB, *gBottomArrA, *gBottomArrB, *gRightArrA, *gRightArrB;
	const static int gNumCoefForOnePointA = 6; //Aij
	const static int gNumCoefForOnePointB = 7; //Bij, C00

	double m_SpareElecEnergyVal;

public:

	//srTRadIntThickBeam(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTMagFldCont* pMagLensCont, void* pPrcPar)
	srTRadIntThickBeam()
    {
		Initialize();

		//gPrec = *((srTParPrecStokesArb*)pPrcPar);
        //gpGenTrjDat = 0;
		//if((pElecBeam != 0) && (pMagFldTrUnif != 0))
		//{
        //	gpGenTrjDat = pMagFldTrUnif->CreateAndSetupNewTrjDat(pElecBeam);
		//}
	}
	~srTRadIntThickBeam()
	{
		DeleteAuxStruct();	
	}

	void Initialize()
	{
		gCoefA = 0;
		gCoefB = 0;
		gBottomArrA = gBottomArrB = gRightArrA = gRightArrB = 0;
		m_SpareElecEnergyVal = 0;
	}
	void DeleteAuxStruct()
	{
		if(gCoefA != 0) { delete gCoefA; gCoefA = 0;}
		if(gCoefB != 0) { delete gCoefB; gCoefB = 0;}

		if(gBottomArrA != 0) { delete gBottomArrA; gBottomArrA = 0;}
		if(gBottomArrB != 0) { delete gBottomArrB; gBottomArrB = 0;}
		if(gRightArrA != 0) { delete gRightArrA; gRightArrA = 0;}
		if(gRightArrB != 0) { delete gRightArrB; gRightArrB = 0;}
	}

	static void ComputeStokes(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes);
    void ComputeTotalStokesDistr(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes);
    void ComputeTotalStokesDistrViaSingleElec(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes);
	void SetupInitialTrajArrays(srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar);
    void ComputeTrajArrays(srTFieldBasedArrays& FldArr, srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont);
	void ComputeOffAxisTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont);
    void DetermineLongPosGridLimits(srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont, double& sStart, double& sEnd);
    void AnalyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ, srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTMagFldCont* pMagLensCont, srTStokesStructAccessData* pStokes);
    void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData* pStokes);
    void CopySymEnergySlice(float* pOrigData, float* pSymData, long Ne, char SymWithRespectToXax, char SymWithRespectToZax);
    //long FindTotalAmOfPointsToCalc(srTStokesStructAccessData* pStokes, char FinalResAreSymOverX, char FinalResAreSymOverZ);
    long long FindTotalAmOfPointsToCalc(srTStokesStructAccessData* pStokes, char FinalResAreSymOverX, char FinalResAreSymOverZ);
	void ComputeStokesAtOneObsPoint(srTEXZY EXZY, srTParPrecStokesArb& Prec, srTStokes& CurSt);
    void ComputeStokesAtOneObsPoint_InternIntens_EvenMesh(srTFieldBasedArrays& FldArr, srTStokes& CurSt);
	
	void ComputeStokesAtOneObsPoint_ExternIntens_AddBotLeft(srTFieldBasedArrays& FldArr, srTStokes& CurSt);
    void ComputeStokesAtOneObsPoint_ExternIntens_AddCen(srTFieldBasedArrays& FldArr, char b_or_r, srTStokes& CurSt);
    void ComputeStokesAtOneObsPoint_ExternIntens_AddBotRight(srTFieldBasedArrays& FldArr, srTStokes& CurSt);
    void ComputeStokesAtOneObsPoint_ExternIntens_AddTopRight(srTFieldBasedArrays& FldArr, srTStokes& CurSt);
    //void ComputeStokesAtOneObsPoint_Intens_PrepAandB(srTFieldBasedArrays& FldArr, int CurNs, int CurNst, int iStart, int itStart, TComplexD* ArrA, TComplexD* ArrB);
    void ComputeStokesAtOneObsPoint_Intens_PrepAandB(srTFieldBasedArrays& FldArr, long long CurNs, long long CurNst, long long iStart, long long itStart, TComplexD* ArrA, TComplexD* ArrB);
	
	//void Integrate1DStokesArr(srTStokes* StokesArr, int Np, double Step, srTStokes* pInitDer, srTStokes* pFinDer, srTStokes& ResSt);
	void Integrate1DStokesArr(srTStokes* StokesArr, long long Np, double Step, srTStokes* pInitDer, srTStokes* pFinDer, srTStokes& ResSt);

	void EstimateExtraObservRangesToIncludeEmittance(srTStokesStructAccessData* pStokes, srTEbmDat* pElecBeam, double* ArrExtraRanges);
    srTSRWRadStructAccessData* CreateNewRadStructWithConstParams(srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTStokesStructAccessData* pStokes, srTWfrSmp*& pWfrSmp);
    double GetNextElecEnergyFromGausDistrib(srTEbmDat& OrigElecBeam, CGenMathRand& RandGen);
	//double UpdateResultStokesData(float* ArrAuxDataS0, float* ArrAuxDataS1, float* ArrAuxDataS2, float* ArrAuxDataS3, srTWfrSmp* pWfrSmp, int iElecEn, srTStokesStructAccessData* pStokes);
	double UpdateResultStokesData(float* ArrAuxDataS0, float* ArrAuxDataS1, float* ArrAuxDataS2, float* ArrAuxDataS3, srTWfrSmp* pWfrSmp, long long iElecEn, srTStokesStructAccessData* pStokes);

    void ComputeStokesAtOneObsPoint_ExternIntens(srTFieldBasedArrays& FldArr, srTStokes& CurSt)
	{
        CurSt.ZeroAll();
		ComputeStokesAtOneObsPoint_ExternIntens_AddBotLeft(FldArr, CurSt);
		ComputeStokesAtOneObsPoint_ExternIntens_AddCen(FldArr, 'b', CurSt);
		ComputeStokesAtOneObsPoint_ExternIntens_AddBotRight(FldArr, CurSt);
		ComputeStokesAtOneObsPoint_ExternIntens_AddCen(FldArr, 'r', CurSt);
		ComputeStokesAtOneObsPoint_ExternIntens_AddTopRight(FldArr, CurSt);

		CurSt *= gAuxPar.Cm;
	}

	void ComputeIntensFuncPartsForInteg2D(double s, double us, TComplexD* ArrA, TComplexD* ArrB, TComplexD* pA_Stokes, TComplexD& B)
	{
		//double xObsE2 = xObs*xObs, zObsE2 = zObs*zObs, xzObs = xObs*zObs;
		double &xObs = gAuxPar.EXZY.x, &yObs = gAuxPar.EXZY.y, &zObs = gAuxPar.EXZY.z;
		double &xObsE2 = gAuxPar.xObsE2, &zObsE2 = gAuxPar.zObsE2, &xzObs = gAuxPar.xzObs;

		B = ArrB[0] + (xObs*ArrB[1]) + (zObs*ArrB[2]) + (xObsE2*ArrB[3]) + (xzObs*ArrB[4]) + (zObsE2*ArrB[5]);
		TComplexD Mult = (1./((yObs - s)*(yObs - us)))*ArrB[6];
		int IndA = 0;
		for(int i=0; i<4; i++)
		{
            pA_Stokes[i] = Mult*(ArrA[IndA] + (xObs*ArrA[IndA + 1]) + (zObs*ArrA[IndA + 2]) + (xObsE2*ArrA[IndA + 3]) + (xzObs*ArrA[IndA + 4]) + (zObsE2*ArrA[IndA + 5]));
            IndA += gNumCoefForOnePointA;
		}
	}

	//void ComputeStokesAtOneObsPoint_FuncForInteg2D_AB(srTFieldBasedArrays& FldArr, int i, int it, TComplexD* pA_Stokes, TComplexD& B)
	void ComputeStokesAtOneObsPoint_FuncForInteg2D_AB(srTFieldBasedArrays& FldArr, long long i, long long it, TComplexD* pA_Stokes, TComplexD& B)
	{
		TComplexD* tA_Stokes = pA_Stokes;
		bool AB_ShouldBeComputed = true;

		if(it < 4)
		{
			if((gBottomArrA != 0) && (gBottomArrB != 0))
			{
                //long OffsB = it*FldArr.Ns + i;
                long long OffsB = it*FldArr.Ns + i;
                B = gBottomArrB[OffsB];
				TComplexD* tBottomArrA = gBottomArrA + (OffsB<<2);
                for(int k=0; k<4; k++) *(tA_Stokes++) = *(tBottomArrA++);
				AB_ShouldBeComputed = false;
			}
		}
		else if(i >= (FldArr.Ns - 4))
		{
			if((gRightArrA != 0) && (gRightArrB != 0))
			{
                //long OffsB = ((it - 4) << 2) + (i - (FldArr.Ns - 4));
                long long OffsB = ((it - 4) << 2) + (i - (FldArr.Ns - 4));
                B = gRightArrB[OffsB];
				TComplexD* tBottomArrA = gBottomArrA + (OffsB<<2);
                for(int k=0; k<4; k++) *(tA_Stokes++) = *(tBottomArrA++);
				AB_ShouldBeComputed = false;
			}
		}
		if(AB_ShouldBeComputed)
		{
			double s = FldArr.sStart + i*FldArr.sStep;
			double st = FldArr.sStart + it*FldArr.sStep;

			//long Offset1 = it*((((FldArr.Ns) << 1) - 1 - it) >> 1) + i;
			long long Offset1 = it*((((FldArr.Ns) << 1) - 1 - it) >> 1) + i;
			//long OffsetCoefA = Offset1*(gNumCoefForOnePointA << 2);
			//long OffsetCoefB = Offset1*gNumCoefForOnePointB;
			long long OffsetCoefA = Offset1*(gNumCoefForOnePointA << 2);
			long long OffsetCoefB = Offset1*gNumCoefForOnePointB;
			ComputeIntensFuncPartsForInteg2D(s, st, gCoefA + OffsetCoefA, gCoefB + OffsetCoefB, pA_Stokes, B);
		}
	}

	void ComputeStokesFromAB(TComplexD* A_Stokes, TComplexD& B, srTStokes& ResSt)
	{
        double ExpReB = exp(B.x);
        TComplexD ExpB(ExpReB*cos(B.y), ExpReB*sin(B.y));

		TComplexD LocA_Stokes[4];
		TComplexD *tA_Stokes = LocA_Stokes;
		for(int k=0; k<4; k++) *(tA_Stokes++) = A_Stokes[k]*ExpB;
		ResSt.s0 = 2*(LocA_Stokes->x);
		ResSt.s1 = 2*((LocA_Stokes+1)->x);
		ResSt.s2 = 2*((LocA_Stokes+2)->x);
		ResSt.s3 = 2*((LocA_Stokes+3)->x);
	}

	//void ComputeStokesAtOneObsPoint_FuncForInteg2D(srTFieldBasedArrays& FldArr, int i, int it, srTStokes& ResSt)
	void ComputeStokesAtOneObsPoint_FuncForInteg2D(srTFieldBasedArrays& FldArr, long long i, long long it, srTStokes& ResSt)
	{
		TComplexD A_Stokes[4], B;
        ComputeStokesAtOneObsPoint_FuncForInteg2D_AB(FldArr, i, it, A_Stokes, B);
		ComputeStokesFromAB(A_Stokes, B, ResSt);
        //double ExpReB = exp(B.x);
        //TComplexD ExpB(ExpReB*cos(B.y), ExpReB*sin(B.y));
		//TComplexD *tA_Stokes = A_Stokes;
		//for(int k=0; k<4; k++) *(tA_Stokes++) *= ExpB;
		//ResSt.s0 = 2*(A_Stokes->x);
		//ResSt.s1 = 2*((A_Stokes+1)->x);
		//ResSt.s2 = 2*((A_Stokes+2)->x);
		//ResSt.s3 = 2*((A_Stokes+3)->x);
	}

    //void Integrate1DStokesFunc_EvenMesh(srTFieldBasedArrays& FldArr, int it, srTStokes& ResSt)
    void Integrate1DStokesFunc_EvenMesh(srTFieldBasedArrays& FldArr, long long it, srTStokes& ResSt)
    {
		//int Ns = FldArr.Ns;
		//int Ns_mi_4 = Ns - 4;
		long long Ns = FldArr.Ns;
		long long Ns_mi_4 = Ns - 4;
        if(it >= Ns_mi_4) 
		{
            srTStokes F;
			if(it == Ns_mi_4) Integrate1DStokesFunc_EvenMesh_4p(FldArr, it, ResSt, F);
			else if(it == (Ns - 3)) Integrate1DStokesFunc_EvenMesh_3p(FldArr, it, ResSt, F);
			else if(it == (Ns - 2)) Integrate1DStokesFunc_EvenMesh_2p(FldArr, it, ResSt, F);
		}
		else if(it != ((it>>1)<<1)) //it is odd, np is even
		{
			srTStokes ExtraSt(0,0,0,0), F2(0,0,0,0);
			Integrate1DStokesFunc_EvenMesh_2p(FldArr, it, ExtraSt, F2);
			Integrate1DStokesFunc_EvenMesh_OddNp(FldArr, it, 1, &F2, ResSt);
			ResSt += ExtraSt;
		}
		else Integrate1DStokesFunc_EvenMesh_OddNp(FldArr, it, 0, 0, ResSt);
		ResSt *= gAuxPar.Cm;
	}

    //void Integrate1DStokesFunc_EvenMesh_2p(srTFieldBasedArrays& FldArr, int it, srTStokes& ResSt, srTStokes& F2)
    void Integrate1DStokesFunc_EvenMesh_2p(srTFieldBasedArrays& FldArr, long long it, srTStokes& ResSt, srTStokes& F2)
	{
        srTStokes F1;
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it, it, F1);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 1, it, F2);
		ResSt = (0.5*FldArr.sStep)*(F1 + F2);
	}
    //void Integrate1DStokesFunc_EvenMesh_3p(srTFieldBasedArrays& FldArr, int it, srTStokes& ResSt, srTStokes& F3)
    void Integrate1DStokesFunc_EvenMesh_3p(srTFieldBasedArrays& FldArr, long long it, srTStokes& ResSt, srTStokes& F3)
	{
        srTStokes F1, F2;
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it, it, F1);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 1, it, F2);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 2, it, F3);
		ResSt = (FldArr.sStep/3.)*(F1 + (4*F2) + F3);
	}
    //void Integrate1DStokesFunc_EvenMesh_4p(srTFieldBasedArrays& FldArr, int it, srTStokes& ResSt, srTStokes& F4)
    void Integrate1DStokesFunc_EvenMesh_4p(srTFieldBasedArrays& FldArr, long long it, srTStokes& ResSt, srTStokes& F4)
	{
        srTStokes F1, F2, F3;
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it, it, F1);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 1, it, F2);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 2, it, F3);
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, it + 3, it, F4);
		ResSt = (0.375*FldArr.sStep)*(F1 + (3*(F2 + F3)) + F4);
	}
    //void Integrate1DStokesFunc_EvenMesh_OddNp(srTFieldBasedArrays& FldArr, int it, int i_Offset, srTStokes* pFe1, srTStokes& ResSt);
    void Integrate1DStokesFunc_EvenMesh_OddNp(srTFieldBasedArrays& FldArr, long long it, long long i_Offset, srTStokes* pFe1, srTStokes& ResSt);

	//void ComputeStokesIntensFuncAtPoint(long is, long ist, srTFieldBasedArrays& FldArr, srTStokesC& A, srTStokesC& B, srTStokesC& C);
	//void ComputeStokesFluxFuncAtPoint(long is, long ist, srTFieldBasedArrays& FldArr, srTStokesC& A, srTStokesC& B, srTStokesC& C);
	//void ComputeStokesFuncAtPoint(long is, long ist, char IntOrFlux, srTFieldBasedArrays& FldArr, srTStokesC& A, srTStokesC& B, srTStokesC& C)
	//{
	//	if(IntOrFlux == 'i') ComputeStokesIntensFuncAtPoint(is, ist, FldArr, A, B, C);
	//	else if(IntOrFlux == 'f') ComputeStokesFluxFuncAtPoint(is, ist, FldArr, A, B, C);
	//}
	//void ComputeExpCoefForOneObsPoint(long is, long ist, double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB);
	void ComputeExpCoefForOneObsPoint(long long is, long long ist, double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB);
	void ComputeExpCoefXZArraysForInteg2D_EvenMesh(double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB);

	void ComputeExpCoefXZArraysForInteg2D(double yObs, double eObs, srTParPrecStokesArb& PrecPar)
    {
        if(PrecPar.MethNo == 1)
        {
            ComputeExpCoefXZArraysForInteg2D_EvenMesh(yObs, eObs, gFldArr, gCoefA, gCoefB);
        }
		//else if
	}

	srTStokes DerivRight3p(srTStokes* f, double dx, int per=1)
	{//to calc. derivative at the start of array: [0]
		return (0.5/dx)*(((-3.)*f[0]) + (4.*f[per]) - f[per << 1]);
	}
	TComplexD DerivRight3p(TComplexD* f, double dx, int per=1)
	{//to calc. derivative at the start of array: [0]
		return (0.5/dx)*(((-3.)*f[0]) + (4.*f[per]) - f[per << 1]);
	}
	srTStokes DerivLeft3p(srTStokes* f, double dx, int per=1)
	{//to calc. derivative at the end of array: [2]
		return (0.5/dx)*(f[0] - (4.*f[per]) + (3.*f[per << 1]));
	}
	TComplexD DerivLeft3p(TComplexD* f, double dx, int per=1)
	{//to calc. derivative at the end of array: [2]
		return (0.5/dx)*(f[0] - (4.*f[per]) + (3.*f[per << 1]));
	}
	TComplexD DerivCen3p(TComplexD* f, double dx, int per=1)
	{
		return (0.5/dx)*(f[per << 1] - f[0]);
	}
	TComplexD Deriv3p(TComplexD* f, int i0, double dx, int per=1)
	{
		if(i0 == 0) return (0.5/dx)*(((-3.)*f[0]) + (4.*f[per]) - f[per << 1]);
		else if(i0 == 1) return (0.5/dx)*(f[per << 1] - f[0]);
		else return (0.5/dx)*(f[0] - (4.*f[per]) + (3.*f[per << 1])); //if(i0 == 2) 
	}

	TComplexD Deriv4p(TComplexD* f, int i0, double dx, int per=1)
	{//to calc. derivative at start or at the end of array (i0 = 0 or i0 = 3 respectively)
		if(i0 == 0) return (((-11)*f[0]) + (18*f[per]) - (9*f[per << 1]) + (2*f[3*per]))/(6*dx);
		else if(i0 == 1) return (((-2)*f[0]) - (3*f[per]) + (6*f[per << 1]) - f[3*per])/(6*dx);
		else if(i0 == 2) return (f[0] - (6*f[per]) + (3*f[per << 1]) + (2*f[3*per]))/(6*dx);
		else return (((-2)*f[0]) + (9*f[per]) - (18*f[per << 1]) + (11*f[3*per]))/(6*dx); //if(i0 == 3) 
	}

	TComplexD Deriv2Cen3p(TComplexD* f, double dx, int per=1)
	{//2nd derivative
		return (f[0] - (2.*f[per]) + f[per << 1])/(dx*dx);
	}
	TComplexD Deriv2_4p(TComplexD* f, int i0, double dx, int per=1)
	{//2nd derivative at start or at the end of array (i0 = 0 or i0 = 3 respectively)
		if(i0 == 0) return ((2.*f[0]) - (5.*f[per]) + (4.*f[per << 1]) - f[3*per])/(dx*dx);
		else if(i0 == 1) return (f[0] - (2.*f[per]) + f[per << 1])/(dx*dx);
		else if(i0 == 2) return (f[per] - (2.*f[per << 1]) + f[3*per])/(dx*dx);
		else return ((4.*f[per]) - (5.*f[per << 1]) + (2.*f[3*per]) - f[0])/(dx*dx);
	}

	void ComputeResidTermA_Stokes(TComplexD* A, TComplexD* B, int i0, double dx, TComplexD* pA_Stokes)
	{//assumes A[4*4] and B[4]; i0 = 0 (for left residual) or 3 (for right residual); pA_Stokes[4]
		TComplexD B0 = B[i0];
		//TComplexD dBdx = Deriv4p(B, i0, dx);
		TComplexD dBdx = Deriv4p(B, i0, dx);
		TComplexD d2Bdx2 = Deriv2_4p(B, i0, dx);

        TComplexD Inv_dBdx = 1./dBdx;
        TComplexD Inv_dBdxE3 = Inv_dBdx/(dBdx*dBdx);
		TComplexD A0, dAdx;
		TComplexD At1, At2;
		TComplexD d2Bdx2_d_dBdxE3;

		for(int i=0; i<4; i++)
		{
            A0 = A[i + (i0<<2)];
            dAdx = Deriv4p(A + i, i0, dx, 4);

			At1 = A0*Inv_dBdx;

			d2Bdx2_d_dBdxE3 = d2Bdx2*Inv_dBdxE3;

			At2 = ((A0*d2Bdx2) - (dAdx*dBdx))*Inv_dBdxE3;

            pA_Stokes[i] = (At1 + At2);
		}
		//double ExpReB0 = exp(B0.x);
        //TComplexD ExpB(ExpReB0*cos(B0.y), ExpReB0*sin(B0.y));
		//return ExpB*((A0/dBdx) + (((A0*d2Bdx2) - (dAdx*dBdx))/(dBdx*dBdx*dBdx)));
	}

	//srTStokes DerivWithExpStokes3p(TComplexD* A, TComplexD* B, double sStep, int ic)
	srTStokes DerivWithExpStokes3p(TComplexD* A, TComplexD* B, double sStep, long long ic)
	{// assumes A[12] and B[3]; ic=0 for start of array, ic=2 for end of array
        double ExpReB0 = exp(B[ic].x);
        TComplexD ExpB(ExpReB0*cos(B[ic].y), ExpReB0*sin(B[ic].y));
        TComplexD dBds;
		if(ic == 0) dBds = DerivRight3p(B, sStep);
        else dBds = DerivLeft3p(B, sStep);

        TComplexD dAds, AuxStokes[4];
        for(int k=0; k<4; k++) 
        {
            if(ic == 0) dAds = DerivRight3p(A + k, sStep, 4);
			else dAds = DerivLeft3p(A + k, sStep, 4);

            AuxStokes[k] = (dAds + (A[(ic<<2) + k]*dBds))*ExpB;
        }
        srTStokes OutStokes;
        OutStokes.s0 = 2*(AuxStokes[0].x);
        OutStokes.s1 = 2*(AuxStokes[1].x);
        OutStokes.s2 = 2*(AuxStokes[2].x);
        OutStokes.s3 = 2*(AuxStokes[3].x);
		return OutStokes;
	}

	//void ComputeStokesAtOneObsPoint_DerivOverS_FromAB(srTFieldBasedArrays& FldArr, int i, int it, int ic, srTStokes& StokesDerF, srTStokes* ArrStokesF)
	void ComputeStokesAtOneObsPoint_DerivOverS_FromAB(srTFieldBasedArrays& FldArr, long long i, long long it, long long ic, srTStokes& StokesDerF, srTStokes* ArrStokesF)
	{//assumes ic=0 for start of array, ic=2 for end of array
		//int iStart = i - ic;
		long long iStart = i - ic;
        TComplexD A_Stokes[12], B[3];
		TComplexD *tA_Stokes = A_Stokes;
		//srTStokes *tStokesF = ArrStokesF;
		for(int k=0; k<3; k++)
		{
            ComputeStokesAtOneObsPoint_FuncForInteg2D_AB(FldArr, iStart, it, tA_Stokes, B[k]);
			ComputeStokesFromAB(tA_Stokes, B[k], ArrStokesF[k]);
			iStart++; tA_Stokes += 4;
		}
        StokesDerF = DerivWithExpStokes3p(A_Stokes, B, FldArr.sStep, ic);
	}

	//void AllocateCoefArraysForInteg2D(long Ns)
	void AllocateCoefArraysForInteg2D(long long Ns)
	{
		if(gCoefA != 0) { delete gCoefA; gCoefA = 0;}
		if(gCoefB != 0) { delete gCoefB; gCoefB = 0;}

		//long NsTriang = ((Ns + 1)*Ns) >> 1;
		//long NpA = (4*gNumCoefForOnePointA)*NsTriang;
		//long NpB = gNumCoefForOnePointB*NsTriang;
		long long NsTriang = ((Ns + 1)*Ns) >> 1;
		long long NpA = (4*gNumCoefForOnePointA)*NsTriang;
		long long NpB = gNumCoefForOnePointB*NsTriang;

		gCoefA = new TComplexD[NpA];
		if(gCoefA == 0) throw MEMORY_ALLOCATION_FAILURE;
		gCoefB = new TComplexD[NpB];
		if(gCoefB == 0) throw MEMORY_ALLOCATION_FAILURE;
	}
	//void AllocateFuncArraysForExternInteg(long Ns)
	void AllocateFuncArraysForExternInteg(long long Ns)
	{
		if(gBottomArrA != 0) { delete gBottomArrA; gBottomArrA = 0;}
		if(gBottomArrB != 0) { delete gBottomArrB; gBottomArrB = 0;}
		if(gRightArrA != 0) { delete gRightArrA; gRightArrA = 0;}
		if(gRightArrB != 0) { delete gRightArrB; gRightArrB = 0;}

		//long NpBottomB = Ns<<2;
		long long NpBottomB = Ns<<2;
		gBottomArrA = new TComplexD[NpBottomB<<2];
		gBottomArrB = new TComplexD[NpBottomB];

		//long NpRightB = (Ns-4)<<2;
		long long NpRightB = (Ns-4)<<2;
		gRightArrA = new TComplexD[NpRightB<<2];
		gRightArrB = new TComplexD[NpRightB];
	}
};

//*************************************************************************

#endif
