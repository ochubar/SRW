/************************************************************************//**
 * File: srstowig.h
 * Description: Calculation of Stokes parameters of Wiggler Radiation (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRSTOWIG_H
#define __SRSTOWIG_H

#include "srtrjdat.h"
#include "srtrjaux.h"
#include "srmlttsk.h"

//*************************************************************************

extern srTYield srYield;
extern srTIntVect gVectWarnNos;

//*************************************************************************

struct srTRadIntWigPrec {
	double PrecFact;
	char TreatInterf;

	srTRadIntWigPrec()
	{
        PrecFact = 1;
        TreatInterf = 0;
	}
};

//*************************************************************************

struct srTRadFuncsAndDers {
	double s;
	double Ph, dPhds, d2Phds2;
	double Ax, dAxds, Az, dAzds;
	double CosPh, SinPh;

	srTRadFuncsAndDers()
	{
		Ph = dPhds = d2Phds2 = 0.;
		Ax = dAxds = Az = dAzds = 0.;
		CosPh = SinPh = 0.;
	}
};

//*************************************************************************

struct srTRadContribInterval {

	double Bx_St, Bz_St, Btx_St, Btz_St, X_St, Z_St, IntBtxE2_St, IntBtzE2_St;
	double Bx_Fi, Bz_Fi, Btx_Fi, Btz_Fi, X_Fi, Z_Fi, IntBtxE2_Fi, IntBtzE2_Fi;
	//int IndSt, IndFi;
	long long IndSt, IndFi;

	double Ph_St, dPhds_St, d2Phds2_St, Ph_Fi, dPhds_Fi, d2Phds2_Fi;
	double Ax_St, dAxds_St, Az_St, dAzds_St, Ax_Fi, dAxds_Fi, Az_Fi, dAzds_Fi;
	double CosPh_St, SinPh_St, CosPh_Fi, SinPh_Fi;
	char PoIsGoodForAn_St, PoIsGoodForAn_Fi;

	double SumReEx, SumImEx, SumReEz, SumImEz;
	srTEFourier E;

	double sc, Xc, Zc;

	srTRadContribInterval() 
	{ 
		IndSt = IndFi = -1;
		Bx_St = Bz_St = Btx_St = Btz_St = X_St = Z_St = IntBtxE2_St = IntBtzE2_St = 0.;
		Bx_Fi = Bz_Fi = Btx_Fi = Btz_Fi = X_Fi = Z_Fi = IntBtxE2_Fi = IntBtzE2_Fi = 0.;

		Ph_St = dPhds_St = d2Phds2_St = Ph_Fi = dPhds_Fi = d2Phds2_Fi = 0.;
		Ax_St = dAxds_St = Az_St = dAzds_St = Ax_Fi = dAxds_Fi = Az_Fi = dAzds_Fi = 0.;
		CosPh_St = SinPh_St = CosPh_Fi = SinPh_Fi = 0.;

		SumReEx = SumImEx = SumReEz = SumImEz = 0.;
	}

	void OutLeftFuncs(srTRadFuncsAndDers& Funcs)
	{
		Funcs.Ph = Ph_St; Funcs.dPhds = dPhds_St; Funcs.d2Phds2 = d2Phds2_St;
		Funcs.Ax = Ax_St; Funcs.dAxds = dAxds_St; 
		Funcs.Az = Az_St; Funcs.dAzds = dAzds_St;
		Funcs.CosPh = CosPh_St; Funcs.SinPh = SinPh_St;
	}
	void OutRightFuncs(srTRadFuncsAndDers& Funcs)
	{
		Funcs.Ph = Ph_Fi; Funcs.dPhds = dPhds_Fi; Funcs.d2Phds2 = d2Phds2_Fi;
		Funcs.Ax = Ax_Fi; Funcs.dAxds = dAxds_Fi; 
		Funcs.Az = Az_Fi; Funcs.dAzds = dAzds_Fi;
		Funcs.CosPh = CosPh_Fi; Funcs.SinPh = SinPh_Fi;
	}
};

//*************************************************************************

class srTRadIntWiggler {

	double PI, TwoPI, ThreePIdTwo, FivePIdFour, HalfPI, One_dTwoPI, PIm10e6, PIm10e6dEnCon, TenEm6dTwoPI; // Constants
	srTCosAndSinComp CosAndSinComp;

	double sIntegRelPrecG, RelTolForAnTermsG;
	char NumIntegG;
	double NormalizingConst;

	srTEXZ EXZ;
	srTFieldBasedArrays FieldBasedArrays;
	srTTrjArraysAux TrjArraysAux;
	srTRadContribInterval *RadContribIntervals, *RadNumIntervals;
	int *MaxPoIndArray, *MinPoIndArray;
	double MaxIntValWithinPointG, MinEstIntValWithinPointG;

	float *CrossTermsContribArray, *tCrossTermsContribG;
	char IncludeCrossTermsG;
	double Gx, Gz;

public:

	srTTrjDat* TrjDatPtr;

	srTWfrSmp DistrInfoDat;
	srTRadIntWigPrec IntWigPrec;

	char LongIntTypeG; // 1- all field range; 2- one period;

	srTRadIntWiggler()
	{
		Initialize();
	}
	void Initialize() 
	{
		DistrInfoDat.CoordUnits = 0; // To ensure m for coord.
		RadContribIntervals = RadNumIntervals = 0;
		MaxPoIndArray = MinPoIndArray = 0;

		CrossTermsContribArray = 0;
		tCrossTermsContribG = 0;

		IncludeCrossTermsG = IntWigPrec.TreatInterf; //this will only work if Initialize() is called after input is finished
		//IncludeCrossTermsG = 0; //1; // Make Switch/Analyzer later
		//OC modified 160304
		//OC modified 130105

		PI = 3.141592653590;
		TwoPI = 2.*PI;
		ThreePIdTwo = 1.5*PI;
		FivePIdFour = 1.25*PI;
		HalfPI = 0.5*PI;
		One_dTwoPI = 1./TwoPI;
		PIm10e6 = PI*1.E+06;
		PIm10e6dEnCon = PIm10e6*0.80654658;
		TenEm6dTwoPI = 1.E-06/TwoPI;

		NumIntegG = 1; // Make Switch/Analizer latrer

		CosAndSinComp.Initialize();
		FieldBasedArrays.ZeroPtrs();
	}

	int ComputeTotalStokesDistr(srTStokesStructAccessData&);
	int ComputeStokesAtPoint(float* pStokes);
	void DetermineRadIntervals(int& AmOfInterv);
	void AnalizeAllRadPoints(int& AmOfMax, int& AmOfMin, int& AmOfNumInterv);
	void MergeAdjacentIntervals(int& AmOfInterv);
	void SetupIntervalsForMcDonaldEstimation(int AmOfMax, int AmOfMin, int& AmOfInterv);
	int ComputeAndAddOneMainTermNum(int TermNo);
	int ComputeAndAddTerminations(int TermNo);
	int ComputeNormalResidual(srTRadFuncsAndDers& Funcs, int NumberOfTerms, srTEFourier& Res);
	void ChooseNumIntervals(int AmOfInterv, int& StartIntervNo, int& FinIntervNo);
	int ComputeCrossTermsAtPointArb(int AmOfInterv, srTStokes& Stokes);
	int ComputeCrossTermsAtPointPer(int StartIntervNo, int FinIntervNo, srTStokes& Stokes);
	int TreatFiniteElecBeamEmittance(srTStokesStructAccessData&, char MainOrCrossTerms, double ElBeamMomFact);
	void ExtractStokesSliceConstE(srTStokesStructAccessData& StokesAccessData, long ie, int StokesNo, char MainOrCrossTerms, float* pOutS);
	void UpdateStokesSliceConstE(float* StokesCmpnArr, long ie, int is, char MainOrCrossTerms, srTStokesStructAccessData& StokesAccessData);
	int TreatFiniteElecBeamEmittanceOneComp1D(float* CmpnArr, double ElBeamMomFact);
	int TreatFiniteElecBeamEmittanceOneComp2D(float* CmpnArr, double ElBeamMomFact);
	void DetermineSingleElecDistrEffSizes2D(float* CmpnArr, double& Mxx, double& Mzz);
	void DetermineSingleElecDistrEffSizes1D(float* CmpnArr, char VsXorZ, double& M_Cen);
	void DetermineResizeBeforeConv2D(double MxxElecEff, double MzzElecEff, double MxxPowSingleE, double MzzPowSingleE, srTRadResize& Resize);
	void DetermineResizeBeforeConv1D(double M_ElecEff, double M_DistrSingleE, char VsXorZ, srTRadResize1D& Resize);
	void ConstructDataForConv2D(float* CmpnArr, float* NewData, long NewNx, long NewNz);
	void ConstructDataForConv1D(float* CmpnArr, float* AuxConvData, long NpOld, long NpNew);
	int PerformConvolutionWithGaussian2D(float* ConvData, long NewNx, long NewNz, double MxxElecEff, double MzzElecEff);
	int PerformConvolutionWithGaussian1D(float* AuxConvData, long NpAux, double M_ElecEff, char VsXorZ);
	void ExtractDataAfterConv2D(float* AuxConvData, long NxAux, long NzAux, float* CmpnArr);
	//void ExtractDataAfterConv1D(float* AuxConvData, long NpAux, long Np, float* CmpnArr);
	void ExtractDataAfterConv1D(float* AuxConvData, long long NpAux, long long Np, float* CmpnArr);
	void SuppressNegativeValues(float* StokesCmpnArr);
	
	int CheckInputConsistency();
	void CheckPossibilityOfFarFieldPeriodicComp(srTGenTrjHndl& TrjHndl, char& LongIntType);
	int SetUpFieldBasedArrays();
	void AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ);
	void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData& StokesAccessData);

	char FindMaxFieldPlane(srTFieldBasedArrays& FieldBasedArrays);
	void SetupCentralValuesInNumInterv(int AmOfInterv);

	void CopySymEnergySlice(float* pOrigData, float* pSymData, char SymWithRespectToXax, char SymWithRespectToZax)
	{
		char ChangeSignS2 = !(SymWithRespectToXax && SymWithRespectToZax);
		char ChangeSignS3 = SymWithRespectToXax;
		float *tOrig = pOrigData, *tSym = pSymData;
		for(int ie=0; ie<DistrInfoDat.nLamb; ie++)
		{
			*(tSym++) = *(tOrig++); *(tSym++) = *(tOrig++);
			*(tSym++) = ChangeSignS2? -(*(tOrig++)) : *(tOrig++);
			*(tSym++) = ChangeSignS3? -(*(tOrig++)) : *(tOrig++);
		}
	}
	//int AllocateIntervalsArray(long Ns)
	int AllocateIntervalsArray(long long Ns)
	{
		DeallocateIntervalsArray();
		RadContribIntervals = new srTRadContribInterval[Ns];
		if(RadContribIntervals == 0) goto SlivaiVodu;
		RadNumIntervals = new srTRadContribInterval[Ns];
		if(RadNumIntervals == 0) goto SlivaiVodu;

		MaxPoIndArray = new int[Ns];
		if(MaxPoIndArray == 0) goto SlivaiVodu;

		MinPoIndArray = new int[Ns];
		if(MinPoIndArray == 0) goto SlivaiVodu;
		return 0;

	SlivaiVodu:
		DeallocateIntervalsArray(); return MEMORY_ALLOCATION_FAILURE;
	}
	void DeallocateIntervalsArray()
	{
		if(RadContribIntervals != 0) delete[] RadContribIntervals; RadContribIntervals = 0;
		if(RadNumIntervals != 0) delete[] RadNumIntervals; RadNumIntervals = 0;

		if(MaxPoIndArray != 0) delete[] MaxPoIndArray; MaxPoIndArray = 0;
		if(MinPoIndArray != 0) delete[] MinPoIndArray; MinPoIndArray = 0;
	}
	void SetIntegPrecLevel()
	{
		const double NominalPrec = 0.02; // To steer
		sIntegRelPrecG = NominalPrec/IntWigPrec.PrecFact;
		//RelTolForAnTermsG = sIntegRelPrecG*0.4; // To steer
		RelTolForAnTermsG = sIntegRelPrecG*1; // To steer

		const double NominalIntVal = 1.E+12;
		const double MinEstCoef = 5.E-03;
		const double NominalCur = 0.2;
		const double NominalE = 6.;
		const double NominalR = 10.;
		double RatI = TrjDatPtr->EbmDat.Current/NominalCur;
		double RatE = TrjDatPtr->EbmDat.Energy/NominalE;
		double RatR = NominalR/DistrInfoDat.yStart;
		MinEstIntValWithinPointG = NominalIntVal*MinEstCoef*RatI*(RatE*RatE)*(RatR*RatR);
	}
	void SetupNormalizingConst()
	{//Assume Spectral flux density is in Photons/(s*mm^2*(0.1%BW))
		const double e_coulomb = 1.602189246E-19;
		const double Alpha = 1./137.0360411; // Fine-structure constant
		NormalizingConst = sqrt(Alpha*(TrjDatPtr->EbmDat.Current)*1.E+09/e_coulomb);
		DistrInfoDat.NormalizingConst = NormalizingConst;
	}
	void E2Stokes(srTEFourier& E, srTStokes& Stokes)
	{
		double LinHor = E.EwX_Re*E.EwX_Re + E.EwX_Im*E.EwX_Im;
		double LinVer = E.EwZ_Re*E.EwZ_Re + E.EwZ_Im*E.EwZ_Im;
		Stokes.s0 = LinHor + LinVer;
		Stokes.s1 = LinHor - LinVer;
		Stokes.s2 = -2.*(E.EwX_Re*E.EwZ_Re + E.EwX_Im*E.EwZ_Im);
		Stokes.s3 = 2.*(-E.EwX_Re*E.EwZ_Im + E.EwX_Im*E.EwZ_Re);
	}

	int SetUpCrossTermsContribArray()
	{
		DeallocateCrossTermsContribArray();
		//long TotPo = (DistrInfoDat.nz*DistrInfoDat.nx*DistrInfoDat.nLamb) << 2;
		long long TotPo = (((long long)DistrInfoDat.nz)*((long long)DistrInfoDat.nx)*((long long)DistrInfoDat.nLamb)) << 2;
		CrossTermsContribArray = new float[TotPo];
		if(CrossTermsContribArray == 0) return MEMORY_ALLOCATION_FAILURE;

		tCrossTermsContribG = CrossTermsContribArray;
		return 0;
	}
	void DeallocateCrossTermsContribArray()
	{
		if(CrossTermsContribArray != 0) delete[] CrossTermsContribArray; CrossTermsContribArray = 0;
	}
	void AddCrossTermsContrib(srTStokesStructAccessData& StokesAccessData)
	{
		float *tGen = StokesAccessData.pBaseSto, *tExtra = CrossTermsContribArray;
		for(long iz=0; iz<StokesAccessData.nz; iz++)
			for(long ix=0; ix<StokesAccessData.nx; ix++)
				for(long ie=0; ie<StokesAccessData.ne; ie++)
				{
					*tGen += *(tExtra++);
					if(*tGen < 0.) *tGen = 0.;
					tGen++;
					for(int is=1; is<4; is++) *(tGen++) += *(tExtra++);
				}
	}
	int SetupThickBeamConsts()
	{
		srTEbmDat& Ebm = TrjDatPtr->EbmDat;
		double y = DistrInfoDat.yStart - Ebm.s0;
		double ye2 = y*y;

		double BufX = Ebm.Mxx + (Ebm.Mxpxp)*ye2 + 2.*(Ebm.Mxxp)*y;
		double BufZ = Ebm.Mzz + (Ebm.Mzpzp)*ye2 + 2.*(Ebm.Mzzp)*y;

		if((BufX <= 0.) || (BufZ <= 0.)) return THICK_EL_BEAM_WAS_NOT_SET_UP;

		Gx = 0.5/BufX;
		Gz = 0.5/BufZ;
		return 0;
	}
	double MaxPhaseDifForThisObsPoint(char NearField)
	{
		double PIdLamb_Inv_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 2.533840802E+06*EXZ.e : 3.14159265359E+09/EXZ.e;

		//int Ns_mi_1 = FieldBasedArrays.Ns - 1;
		long long Ns_mi_1 = FieldBasedArrays.Ns - 1;
		double PerLen = FieldBasedArrays.sStep*Ns_mi_1;
		double yObs = DistrInfoDat.yStart;
		double BufMag = ((FieldBasedArrays.IntBtxE2Arr)[Ns_mi_1] - *(FieldBasedArrays.IntBtxE2Arr))
					   +((FieldBasedArrays.IntBtzE2Arr)[Ns_mi_1] - *(FieldBasedArrays.IntBtzE2Arr));
		if(NearField)
		{
			double BufDist = 1./(yObs - FieldBasedArrays.sStart);
			double BufDist_End = 1./(yObs - (FieldBasedArrays.sStart + PerLen));

			double BufX_End = EXZ.x - (FieldBasedArrays.XArr)[Ns_mi_1], BufX = EXZ.x - *(FieldBasedArrays.XArr);
			double BufZ_End = EXZ.z - (FieldBasedArrays.ZArr)[Ns_mi_1], BufZ = EXZ.z - *(FieldBasedArrays.ZArr);
			
			return PIdLamb_Inv_m*((TrjDatPtr->EbmDat.GammaEm2)*PerLen + BufMag + (BufX_End*BufX_End + BufZ_End*BufZ_End)*BufDist_End - (BufX*BufX + BufZ*BufZ)*BufDist);
		}
		else
		{
			double xAng = EXZ.x/yObs, zAng = EXZ.z/yObs;
			double xDif = (FieldBasedArrays.XArr)[Ns_mi_1] - *(FieldBasedArrays.XArr);
			double zDif = (FieldBasedArrays.ZArr)[Ns_mi_1] - *(FieldBasedArrays.ZArr);
			return PIdLamb_Inv_m*((TrjDatPtr->EbmDat.GammaEm2 + xAng*xAng + zAng*zAng)*PerLen + BufMag - 2.*(xAng*xDif + zAng*zDif));
		}
	}
};

//*************************************************************************

#endif
