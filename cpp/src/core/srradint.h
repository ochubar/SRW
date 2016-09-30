/************************************************************************//**
 * File: srradint.h
 * Description: SR calculation method from ~Arbitrary Transversely-Uniform Magnetic Field, in Near-Field observation conditions (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRRADINT_H
#define __SRRADINT_H

//#ifdef __IGOR_PRO__
#include "srsend.h"
//#endif

#include "srstraux.h"
#include "srmagfld.h"
#include "srtrjdat.h"
//#include "srmamet.h"
#include "gmmeth.h"
#include "srprdint.h"
#include "smartptr.h"
#include "gmvect.h"

//*************************************************************************

typedef CSmartPtr<srTPartAutoRadInt> srTPartAutoRadIntHndl;

#ifdef __GCC__
typedef vector<srTPartAutoRadIntHndl> srTPartAutoRadIntHndlVect;
#else
typedef vector<srTPartAutoRadIntHndl, allocator<srTPartAutoRadIntHndl> > srTPartAutoRadIntHndlVect;
#endif

struct srTParPrecElecFld;

//*************************************************************************

class srTRadInt {

	double* IntegSubdArray;
	int LenIntegSubdArray;

	srLambXYZ ObsCoor;
	char x_or_z;
	double NormalizingConst;

	complex<double> InitDerMan[2], FinDerMan[2]; // To change
	char ComputeDer;

	double *BtxArr, *XArr, *IntBtxE2Arr, *BtzArr, *ZArr, *IntBtzE2Arr, *BxArr, *BzArr;

	double PI, TwoPI, ThreePIdTwo, FivePIdFour, HalfPI, One_dTwoPI, PIm10e6, PIm10e6dEnCon, TenEm6dTwoPI; // Constants

	double a2c, a4c, a6c, a8c, a10c, a12c;
	double a3s, a5s, a7s, a9s, a11s, a13s;

	double *BtxArrP[50], *XArrP[50], *IntBtxE2ArrP[50], *BxArrP[50];
	double *BtzArrP[50], *ZArrP[50], *IntBtzE2ArrP[50], *BzArrP[50];
	//int AmOfPointsOnLevel[50];
	//int NumberOfLevelsFilled;
	long long AmOfPointsOnLevel[50];
	long long NumberOfLevelsFilled;
	double MaxFluxDensVal, CurrentAbsPrec;
	int MaxLevelForMeth_10_11;

	srTStNoFiNoVect StNoFiNoVectArr[50];

	srTRadIntervVal ComputedRadIntervAuto[512], NumRadIntervAuto[512];
	double IntXReArrAuto[513], IntXImArrAuto[513], IntZReArrAuto[513], IntZImArrAuto[513];
	double PhArrAuto[513], dPhdsArrAuto[513], d2Phds2ArrAuto[513], AxArrAuto[513], dAxdsArrAuto[513], AzArrAuto[513], dAzdsArrAuto[513];
	char PointIsGoodForAnAuto[513];
	double d2Phds2ddPhdsE2ArrAuto[513];
	
	srTPartAutoRadIntHndlVect PartAutoRadIntHndlVect;

	complex<double> SmallContComplex[4];
	double SmallContDouble[4];

	char ComputeNormalDerivative; // For Diffraction computation in Near Field
	TVector3d SurfNorm;
	double SmallContDoubleForNormDer[4];
	double InitDerXReND, InitDerXImND, InitDerZReND, InitDerZImND, FinDerXReND, FinDerXImND, FinDerZReND, FinDerZImND;

	double EstimatedAbsoluteTolerance;
	char m_CalcResidTerminTerms;

public:

	srTSend* pSend;
	//DLL_IMPLEMENT

	srTTrjDat* TrjDatPtr;
	srTWfrSmp DistrInfoDat;

	double *RadDistrPhotonFluxHorPol, *RadDistrPhotonFluxVerPol;
	double *RadDistrPhotonFluxHorPolTravers, *RadDistrPhotonFluxVerPolTravers;
	complex<double> *RadDistrFieldFourierHorPol, *RadDistrFieldFourierVerPol;
	complex<double> *RadDistrFieldFourierHorPolTravers, *RadDistrFieldFourierVerPolTravers;
	double *IntegratedPhotonFlux;
	double *AuxPhaseArray;

	double *dExdlRe, *dExdlIm, *dEzdlRe, *dEzdlIm;
	double *dExdlReTravers, *dExdlImTravers, *dEzdlReTravers, *dEzdlImTravers;

	double sIntegStart, sIntegFin, sIntegStep, sIntegRelPrec, sIntegStep_Input;
	double Inv_sIntegStep;
	//int AmOfPointsForManIntegr;
	long long AmOfPointsForManIntegr;
	char sIntegMethod, UseManualSlower;
	double MaxMemAvail, CurrMemAvail;
	//long MaxNumPoToSave;
	long long MaxNumPoToSave;
	
	char TryToApplyNearFieldResidual;
	char TryToApplyNearFieldResidual_AtRight;

	char TrjDataContShouldBeRebuild, ProbablyTheSameLoop;

	srTRadInt()
	{
		Initialize();	
	}
	~srTRadInt()
	{
		DeallocateMemForRadDistr();	
	}

	void Initialize(); // Same as constructor (to solve the problem with CW)

	inline void SetInputTrjData(srTTrjDat*);
	inline void SetupNormalizingConst();

	int AllocateMemForRadDistr();
	inline void SetIntegSubdArray(double* =0, int =0);

	inline double FindCentrPointForRadInteg(double s1, double s2, double sacc);
	inline double FindCandidateToCentrPointNewt(double s1, double s2, double sacc, int JMAX, short& Inacc);
	inline void FuncForCenPoint(double, double*, double*);

	inline int ComputeTotalRadDistr();
	int ComputeTotalRadDistrLoops();
	int ComputeTotalRadDistrDirectOut(srTSRWRadStructAccessData&, char showProgressInd = 1);
	inline int GenRadIntegration(complex<double>*, srTEFourier*);
	inline int RadIntegrationAutoByPieces(complex<double>*);
	inline int RadIntegrationResiduals(complex<double>*, srTEFourier*);
	int ComputeNormalResidual(double, int, complex<double>*, srTEFourier*);

	inline double AbsE2OfComplex(complex<double>&);
	inline void InitializeTraverses();

	inline void FunForRadInt(double, complex<double>*);

	inline void ReformatRadDistrCompRes();
	inline double IntegrateDistrLine(int);

	inline void DeallocateMemForRadDistr();
	inline void AnglesFromMrToRad();
	inline void AnglesFromRadToMr();

	inline int ScanPhase();
	inline double PhaseFun(double, int);
	inline int SetupParamForManualIntegr();

	int RadIntegrationManualSlower(double&, double&, double&, double&, srTEFourier*);

	int RadIntegrationManualFaster0(double&, double&, double&, double&, srTEFourier*);
	int RadIntegrationManualFaster1(double&, double&, double&, double&, srTEFourier*);
	int RadIntegrationManualFaster2(double&, double&, double&, double&, srTEFourier*);

	inline void FunForRadIntWithDer(double, complex<double>*, complex<double>*);
	
	int RadIntegrationAuto1(double&, double&, double&, double&, srTEFourier*);
	//int RadIntegrationAuto1M(double sStart, double sEnd, double* FunArr, double* EdgeDerArr, int AmOfInitPo, int NextLevNo, double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm);
	int RadIntegrationAuto1M(double sStart, double sEnd, double* FunArr, double* EdgeDerArr, long long AmOfInitPo, int NextLevNo, double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm);
	int RadIntegrationAuto2(double&, double&, double&, double&, srTEFourier*);

	inline void CosAndSin(double x, double& Cos, double& Sin);

	inline void AxAzPhFarField(double s, double& Ax, double& Az, double& Ph);
	inline void AxAzPhNearField(double s, double& Ax, double& Az, double& Ph);
	inline void AxAzPhFarField2(int LevelNo, int IndxOnLevel, double s, double& Ax, double& Az, double& Ph);
	inline void AxAzPhNearField2(int LevelNo, int IndxOnLevel, double s, double& Ax, double& Az, double& Ph);

	//inline int FillNextLevel(int LevelNo, double sStart, double sEnd, long Np);
	inline int FillNextLevel(int LevelNo, double sStart, double sEnd, long long Np);
	//int FillNextLevelPart(int LevelNo, double sStart, double sEnd, long Np, double*** TrjPtrs);
	int FillNextLevelPart(int LevelNo, double sStart, double sEnd, long long Np, double*** TrjPtrs);

	inline int SetupRadCompStructures();
	inline int SetupRadCompStructMethAuto2();

	inline void EnableNormDerComp(TVector3d&);
	inline void DisableNormDerComp();

	inline int EvaluateMemAvailAfterTrjComp();

	srTEFourier ComputeRightPartOfRadIntegralInNearField(double sSt, double Btx, double xSt, double IntBtE2xSt, double Btz, double zSt, double IntBtE2zSt);
	inline srTEFourier RadFuncInDriftSpace(double s, double sSt, double Btx, double xSt, double Btz, double zSt, double IntBtE2xzSt, double& Phase);
	void DetermineIntegIntervalsForRightResidual(double sStGen, const int AmOfParts, double* EdgePoints);

	int ComputePreResid(double sStNew, double sStOld, double Btx, double Crdx, double IntBtE2x, double Btz, double Crdz, double IntBtE2z, srTEFourier*, srTEFourier*, char LeftOrRight);
	double X_FreeSpace(double s, double s0, double Bt0, double x0) { return x0 + Bt0*(s-s0);}

	int CheckInputConsistency();
	inline void EstimateAbsoluteTolerance();

	void AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ);
	void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTSRWRadStructAccessData&);
	inline void CopySymEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, char SymWithRespectToXax, char SymWithRespectToZax);

	inline void CheckFurtherSubdNeed(srTEFourierVect*, long, long, char*, char*);
	int CheckFurtherSubdNeedForOneCoord(srTEFourier**);

	int RadInterpolationOnePointXZ(srTEFourierVect*, int, int, double, double, srTEFourier*);

    void SetPrecParams(srTParPrecElecFld*);
    void ComputeElectricFieldFreqDomain(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecElecFld* pPrecElecFld, srTSRWRadStructAccessData* pWfr, char showProgressInd = 1);
};

//*************************************************************************

inline double srTRadInt::IntegrateDistrLine(int kk)
{
	double Step;
	if(DistrInfoDat.IntegrDistrMap == MapVsHor)
		Step = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
	else Step = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;

	int nn = (DistrInfoDat.IntegrDistrMap == MapVsHor)? DistrInfoDat.nz : DistrInfoDat.nx;
	double* LocDistrPtr = &(RadDistrPhotonFluxHorPol[kk*nn]);
	double LocSum = 0.5*(*(LocDistrPtr++));
	for(int k=1; k<(nn - 1); k++) LocSum += *(LocDistrPtr++);
	LocSum += 0.5*(*(LocDistrPtr));
	return LocSum * Step;
}

//*************************************************************************

inline void srTRadInt::SetInputTrjData(srTTrjDat* InTrjDatPtr) 
{ 
	TrjDatPtr = InTrjDatPtr;

	double FieldDataFin = TrjDatPtr->sStart + TrjDatPtr->sStep*(TrjDatPtr->LenFieldData - 1);
	if(sIntegStart < TrjDatPtr->sStart) sIntegStart = TrjDatPtr->sStart;
	if(sIntegFin > FieldDataFin) sIntegFin = FieldDataFin;
}

//*************************************************************************

inline void srTRadInt::InitializeTraverses()
{
	RadDistrPhotonFluxHorPolTravers = RadDistrPhotonFluxHorPol;
	RadDistrPhotonFluxVerPolTravers = RadDistrPhotonFluxVerPol;
	RadDistrFieldFourierHorPolTravers = RadDistrFieldFourierHorPol;
	RadDistrFieldFourierVerPolTravers = RadDistrFieldFourierVerPol;

	if(ComputeNormalDerivative)
	{
		dExdlReTravers = dExdlRe;
		dExdlImTravers = dExdlIm;
		dEzdlReTravers = dEzdlRe;
		dEzdlImTravers = dEzdlIm;
	}
}

//*************************************************************************

inline double srTRadInt::AbsE2OfComplex(complex<double>& ComplexNumber)
{
	double Re = ComplexNumber.real();
	double Im = ComplexNumber.imag();
	return Re*Re + Im*Im;
}

//*************************************************************************

inline void srTRadInt::SetupNormalizingConst()
{// We assume: 
//Wavelength in nm;
//Length in m;
//Then Electric field Fourier component is in ESU units (following Landau presentation),
//Spectral flux density is in Photons/(s*mr^2*(0.1%BW)) or in Photons/(s*mm^2*(0.1%BW))

	//const double e_esu = 4.80324214E-10; // Charge of electron in esu
	const double e_coulomb = 1.602189246E-19;
	//const double c = 2.99792458E+10; // Speed of light in cm/s
	const double Alpha = 1./137.0360411; // Fine-structure constant
	//const double PI = 3.141592653590;
	//const double TwoPI = 2.*PI;

	if((DistrInfoDat.DistrValType == StokesParam) || (DistrInfoDat.DistrValType == FieldFourier))
	{
		if(DistrInfoDat.CoordOrAngPresentation == CoordPres) 
			NormalizingConst = sqrt(Alpha*(TrjDatPtr->EbmDat.Current)*(1.E+09)/e_coulomb);
		else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
			NormalizingConst = sqrt(Alpha*(TrjDatPtr->EbmDat.Current)*(1.E+03)/e_coulomb);
	}

	DistrInfoDat.NormalizingConst = NormalizingConst;
}

//*************************************************************************

inline int srTRadInt::ComputeTotalRadDistr()
{
	int result;
	if(result = SetupRadCompStructures()) return result;
	if(DistrInfoDat.ShowPhaseOnly) return ScanPhase();

	if(DistrInfoDat.OnlyOnePoint)
	{
		ObsCoor.Lamb = DistrInfoDat.LambStart;
		ObsCoor.x = DistrInfoDat.xStart; ObsCoor.y = DistrInfoDat.yStart; ObsCoor.z = DistrInfoDat.zStart;

		complex<double> RadIntegValues[2];
		srTEFourier EwNormDer;

		if(result = GenRadIntegration(RadIntegValues, &EwNormDer)) return result;

		*RadDistrFieldFourierHorPol = *RadIntegValues;
		*RadDistrFieldFourierVerPol = *(RadIntegValues + 1);

		if(ComputeNormalDerivative)
		{
			*dExdlRe = EwNormDer.EwX_Re;
			*dExdlIm = EwNormDer.EwX_Im;
			*dEzdlRe = EwNormDer.EwZ_Re;
			*dEzdlIm = EwNormDer.EwZ_Im;
		}

		return 0;
	}
	else return ComputeTotalRadDistrLoops();
}

//*************************************************************************

inline int srTRadInt::GenRadIntegration(complex<double>* RadIntegValues, srTEFourier* pEwNormDer)
{// Put here more functionality (switching to different methods) later
	int result;
	complex<double> ResidVal[2];
	srTEFourier EwNormDer;

	double IntXRe = 0., IntXIm = 0.;
	double IntZRe = 0., IntZIm = 0.;

	if(m_CalcResidTerminTerms > 0)
	{
		if(result = RadIntegrationResiduals(ResidVal, &EwNormDer)) return result;
		IntXRe = (*ResidVal).real(); IntXIm = (*ResidVal).imag();
		IntZRe = (*(ResidVal+1)).real(); IntZIm = (*(ResidVal+1)).imag();
	}
//TEST!
//		IntXRe = IntXIm = 0.;
//		IntZRe = IntZIm = 0.;
//END TEST

	if(sIntegMethod == 0) 
	{
		if(UseManualSlower) { if(result = RadIntegrationManualSlower(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
		else { if(result = RadIntegrationManualFaster0(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
	}
	else if(sIntegMethod == 1) 
	{ 
		if(UseManualSlower) { if(result = RadIntegrationManualSlower(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
		else { if(result = RadIntegrationManualFaster1(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
	}
	else if(sIntegMethod == 2) 
	{ 
		if(UseManualSlower) { if(result = RadIntegrationManualSlower(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
		else { if(result = RadIntegrationManualFaster2(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
	}
	else if(sIntegMethod == 10) { if(result = RadIntegrationAuto1(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}
	else if(sIntegMethod == 11) { if(result = RadIntegrationAuto2(IntXRe, IntXIm, IntZRe, IntZIm, &EwNormDer)) return result;}

	complex<double> IntX(IntXRe, IntXIm), IntZ(IntZRe, IntZIm);
	*(RadIntegValues++) = IntX;
	*RadIntegValues = IntZ;

	if(ComputeNormalDerivative)
	{
		pEwNormDer->EwX_Re = EwNormDer.EwX_Re;
		pEwNormDer->EwX_Im = EwNormDer.EwX_Im;
		pEwNormDer->EwZ_Re = EwNormDer.EwZ_Re;
		pEwNormDer->EwZ_Im = EwNormDer.EwZ_Im;
	}
	return 0;
}

//*************************************************************************

inline int srTRadInt::RadIntegrationAutoByPieces(complex<double>* RadIntegValues)
{
	int LenVal = ((DistrInfoDat.DistrPolariz == HorOnly) || (DistrInfoDat.DistrPolariz == VerOnly))? 1 : 2;

	double RelPrecAndLimitsArray[4];
	complex<double> ActualIntegVal[12];
	complex<double>* IntegVal[6];
 	for(int j=0; j<6; j++) IntegVal[j] = &(ActualIntegVal[j*2]);
	short ElemCompNotFinished[2];

	complex<double> &xPartVal = ActualIntegVal[0], &zPartVal = ActualIntegVal[1];
	complex<double> &xTotVal = *RadIntegValues, &zTotVal = *(RadIntegValues + 1);
	double *RelPrecAndLimitsArrayPtr = RelPrecAndLimitsArray;
	*(RelPrecAndLimitsArrayPtr++) = sIntegRelPrec;
	if(LenVal == 2) *(RelPrecAndLimitsArrayPtr++) = sIntegRelPrec;

	double *IntegSubdArrayPtr = IntegSubdArray;
	for(int i=0; i<LenIntegSubdArray-1; i++)
	{
		*RelPrecAndLimitsArrayPtr = *IntegSubdArrayPtr;
		*(RelPrecAndLimitsArrayPtr+1) = *(++IntegSubdArrayPtr);
		FormalOneFoldInteg(this, &srTRadInt::FunForRadInt, LenVal, RelPrecAndLimitsArray, ElemCompNotFinished, IntegVal);
		xTotVal += xPartVal; zTotVal += zPartVal;
	}
	return 1;
}

//*************************************************************************

inline int srTRadInt::RadIntegrationResiduals(complex<double>* RadIntegValues, srTEFourier* pEwNormDer)
{
	complex<double> RightResidValues[2], LeftResidValues[2];
	complex<double> cZero(0., 0.);
	*RightResidValues = cZero; *(RightResidValues + 1) = cZero;
	*LeftResidValues = cZero; *(LeftResidValues + 1) = cZero;

	srTEFourier EwNormDerRightResid, EwNormDerLeftResid;

	int NumberOfTerms = 3;
	int result;

	if((m_CalcResidTerminTerms == 1) || (m_CalcResidTerminTerms == 2))
	{
		ComputeDer = 1;
		if(result = ComputeNormalResidual(sIntegStart, NumberOfTerms, LeftResidValues, &EwNormDerLeftResid)) return result;
	}

	if((m_CalcResidTerminTerms == 1) || (m_CalcResidTerminTerms == 3))
	{
		ComputeDer = 2;
		if(result = ComputeNormalResidual(sIntegFin, NumberOfTerms, RightResidValues, &EwNormDerRightResid)) return result;
	}

	for(int i=0; i<2; i++) *(RadIntegValues++) = (RightResidValues[i] - LeftResidValues[i]);

	pEwNormDer->EwX_Re = EwNormDerRightResid.EwX_Re - EwNormDerLeftResid.EwX_Re;
	pEwNormDer->EwX_Im = EwNormDerRightResid.EwX_Im - EwNormDerLeftResid.EwX_Im;
	pEwNormDer->EwZ_Re = EwNormDerRightResid.EwZ_Re - EwNormDerLeftResid.EwZ_Re;
	pEwNormDer->EwZ_Im = EwNormDerRightResid.EwZ_Im - EwNormDerLeftResid.EwZ_Im;
	return 0;
}

//*************************************************************************

inline double srTRadInt::FindCentrPointForRadInteg(double s1, double s2, double sPrec)
{
	int J_Max = 20;
	short Inacc = 0;
	double CandidateToCentrPoint = FindCandidateToCentrPointNewt(s1, s2, sPrec, J_Max, Inacc);
	if(!Inacc) return CandidateToCentrPoint;
	else
	{
		double F1, F2, dF;
		FuncForCenPoint(s1, &F1, &dF); F1 = ::fabs(F1);
		FuncForCenPoint(s2, &F2, &dF); F2 = ::fabs(F2);
		return (F1<F2)? s1 : s2;
	}
}

//*************************************************************************

inline double srTRadInt::FindCandidateToCentrPointNewt(double x1, double x2, double xacc, int JMAX, short& Inaccurate)
{
	int j;
	double df,dx,f,rtn;
   
	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) 
	{
		FuncForCenPoint(rtn,&f,&df);
		dx=f/df;
 		rtn -= dx;
		if(x1-rtn > 0.) { Inaccurate = 1; return x1;}
		if(x2-rtn < 0.) { Inaccurate = 1; return x2;}
		if(::fabs(dx) < xacc) return rtn;
	}
	Inaccurate = 1;
	return rtn;
}

//*************************************************************************

inline void srTRadInt::FuncForCenPoint(double s, double* F, double* dF)
{// Not used currently
	double Btx, Btz, X, Z, Bx, Bz;
	TrjDatPtr->CompTrjDataDerivedAtPointPowDens(s, Btx, Btz, X, Z, Bx, Bz);
	
	double Fld, Bt, Crd;
	if(x_or_z == 'x')
	{
		Fld = Bz; Bt = Btx; Crd = X;
	}
	else
	{
		Fld = Bx; Bt = Btz; Crd = Z;
	}

	int Sign = -1;
	double ObCrd = ObsCoor.x;
	if(x_or_z == 'z')
	{
		Sign = 1; ObCrd = ObsCoor.z;
	}
	double ObCrd_mi_Crd = ObCrd - Crd;
	double ObsY_mi_s = ObsCoor.y - s;

	*F = Bt - ObCrd_mi_Crd/ObsY_mi_s;
	*dF = Sign*TrjDatPtr->BetaNormConst*Fld - (ObCrd_mi_Crd - Bt*ObsY_mi_s)/(ObsY_mi_s*ObsY_mi_s);
}

//*************************************************************************

inline void srTRadInt::SetIntegSubdArray(double* InSubdArray, int InLenSubdArray)
{
	int DefaultLenIntegSubdArray = 19; // To steer!

	if((InSubdArray != 0) && (InLenSubdArray != 0))
	{
		LenIntegSubdArray = InLenSubdArray;
		IntegSubdArray = new double[LenIntegSubdArray];
		for(int i=0; i<LenIntegSubdArray; i++) IntegSubdArray[i] = InSubdArray[i];
	}
	else
	{
		LenIntegSubdArray = DefaultLenIntegSubdArray;
		IntegSubdArray = new double[LenIntegSubdArray];
		double SubdStep = (sIntegFin - sIntegStart)/(LenIntegSubdArray - 1);
		double BufVal = IntegSubdArray[0] = sIntegStart;
		for(int i=1; i<LenIntegSubdArray; i++)
		{
			BufVal += SubdStep;
			IntegSubdArray[i] = BufVal;
		}
	}
}

//*************************************************************************

inline void srTRadInt::FunForRadIntWithDer(double s, complex<double>* FunPtr, complex<double>* DerPtr)
{
	double dBxds=0., dBzds=0., Bx=0., Bz=0., Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;

	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('x', s, dBzds, Bz, Btx, Crdx, IntBtE2x);
	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('z', s, dBxds, Bx, Btz, Crdz, IntBtE2z);

	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double Nx, Nz, Ax, dAxds, dAzds, Az, Phase, dPhds, Btx_mi_Nx, Btz_mi_Nz;
	complex<double> *FunTravers = FunPtr, *DerTravers = DerPtr;

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		double One_d_ymis = 1./(ObsCoor.y - s);
		double xObs_mi_x = xObs - Crdx, zObs_mi_z = zObs - Crdz;
		Phase = PIm10e9_d_Lamb*(s*TrjDatPtr->EbmDat.GammaEm2 + IntBtE2x + IntBtE2z + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		Phase -= TwoPI*int(Phase/TwoPI);
		complex<double> ExpIPhase(cos(Phase), sin(Phase));

		Nx = xObs_mi_x*One_d_ymis; Btx_mi_Nx = Btx - Nx; Ax = Btx_mi_Nx*One_d_ymis;
		Nz = zObs_mi_z*One_d_ymis; Btz_mi_Nz = Btz - Nz; Az = Btz_mi_Nz*One_d_ymis;

		*(FunTravers++) = Ax*ExpIPhase;
		*FunTravers = Az*ExpIPhase;

		dPhds = PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
		dAxds = (2.*Ax + (TrjDatPtr->BetaNormConst)*Bz)*One_d_ymis;
		dAzds = (2.*Az + (-TrjDatPtr->BetaNormConst)*Bx)*One_d_ymis;

		complex<double> DerPreExpX(dAxds, Ax*dPhds), DerPreExpZ(dAzds, Az*dPhds);
		*(DerTravers++) = DerPreExpX*ExpIPhase;
		*DerTravers = DerPreExpZ*ExpIPhase;
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		Ax = Btx - xObs; Az = Btz - zObs;
		Phase = PIm10e9_d_Lamb*(s*(TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs) + IntBtE2x + IntBtE2z - 2.*(xObs*Crdx + zObs*Crdz));
		Phase -= TwoPI*int(Phase/TwoPI);
		complex<double> ExpIPhase(cos(Phase), sin(Phase));

		*(FunTravers++) = Ax*ExpIPhase;
		*FunTravers = Az*ExpIPhase;

		dPhds = PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Ax*Ax + Az*Az);
		dAxds = (TrjDatPtr->BetaNormConst)*Bz;
		dAzds = (-TrjDatPtr->BetaNormConst)*Bx;

		complex<double> DerPreExpX(dAxds, Ax*dPhds), DerPreExpZ(dAzds, Az*dPhds);
		*(DerTravers++) = DerPreExpX*ExpIPhase;
		*DerTravers = DerPreExpZ*ExpIPhase;
	}
}

//*************************************************************************

inline void srTRadInt::DeallocateMemForRadDistr()
{
	if(DistrInfoDat.RadDistrDataContShouldBeRebuild)
	{
		if(RadDistrPhotonFluxHorPol != SmallContDouble)
		{
			if(RadDistrPhotonFluxHorPol != 0) { delete[] RadDistrPhotonFluxHorPol; RadDistrPhotonFluxHorPol = 0;}
			if(RadDistrPhotonFluxVerPol != 0) { delete[] RadDistrPhotonFluxVerPol; RadDistrPhotonFluxVerPol = 0;}
		}
		if(RadDistrFieldFourierHorPol != SmallContComplex)
		{
			if(RadDistrFieldFourierHorPol != 0) { delete[] RadDistrFieldFourierHorPol; RadDistrFieldFourierHorPol = 0;}
			if(RadDistrFieldFourierVerPol != 0) { delete[] RadDistrFieldFourierVerPol; RadDistrFieldFourierVerPol = 0;}
		}

		if(ComputeNormalDerivative)
		{
			if(dExdlRe != SmallContDoubleForNormDer)
			{
				if(dExdlRe != 0) { delete[] dExdlRe; dExdlRe = dExdlReTravers = 0;}
				if(dExdlIm != 0) { delete[] dExdlIm; dExdlIm = dExdlImTravers = 0;}
				if(dEzdlRe != 0) { delete[] dEzdlRe; dEzdlRe = dEzdlReTravers = 0;}
				if(dEzdlIm != 0) { delete[] dEzdlIm; dEzdlIm = dEzdlImTravers = 0;}
			}
		}
	}
	if(AuxPhaseArray != 0) { delete[] AuxPhaseArray; AuxPhaseArray = 0;}

	if(TrjDataContShouldBeRebuild)
	{
		if(BtxArr != 0) { delete[] BtxArr; BtxArr = 0;}
		if(XArr != 0) { delete[] XArr; XArr = 0;}
		if(IntBtxE2Arr != 0) { delete[] IntBtxE2Arr; IntBtxE2Arr = 0;}
		if(BtzArr != 0) { delete[] BtzArr; BtzArr = 0;}
		if(ZArr != 0) { delete[] ZArr; ZArr = 0;}
		if(IntBtzE2Arr != 0) { delete[] IntBtzE2Arr; IntBtzE2Arr = 0;}
		if(BxArr != 0) { delete[] BxArr; BxArr = 0;}
		if(BzArr != 0) { delete[] BzArr; BzArr = 0;}

		for(int i=0; i<NumberOfLevelsFilled; i++)
		{
			double*& BtxArrP_i = BtxArrP[i];
			if(BtxArrP_i != 0) 
			{ 
				delete[] BtxArrP_i; 
				BtxArrP_i = 0;
			}
			AmOfPointsOnLevel[i] = 0;
		}

		for(int k=0; k<50; k++)
		{
			srTStNoFiNoVect& StNoFiNoVect = StNoFiNoVectArr[k];
			if(!StNoFiNoVect.empty())
			{
				double*& BtxArrP_k = BtxArrP[k];
				if(BtxArrP_k != 0) { delete[] BtxArrP_k; BtxArrP_k = 0;}
				AmOfPointsOnLevel[k] = 0;
				StNoFiNoVect.erase(StNoFiNoVect.begin(), StNoFiNoVect.end());
			}
		}
		NumberOfLevelsFilled = 0;
		PartAutoRadIntHndlVect.erase(PartAutoRadIntHndlVect.begin(), PartAutoRadIntHndlVect.end());
	}
}

//*************************************************************************

inline void srTRadInt::ReformatRadDistrCompRes()
{
	if(DistrInfoDat.ShowPhaseOnly) return;

	if(DistrInfoDat.DistrValType == StokesParam)
	{
		int TotAmOfData = DistrInfoDat.nLamb * DistrInfoDat.ny * DistrInfoDat.nz * DistrInfoDat.nx;
		double s0, s1, s2, s3;

		complex<double>* LocHorPolTravers = RadDistrFieldFourierHorPol;
		complex<double>* LocVerPolTravers = RadDistrFieldFourierVerPol;

		for(int k=0; k<TotAmOfData; k++)
		{
			double HorRe = (*LocHorPolTravers).real();
			double HorIm = (*LocHorPolTravers).imag();
			double VerRe = (*LocVerPolTravers).real();
			double VerIm = (*LocVerPolTravers).imag();

			double Hor = HorRe*HorRe + HorIm*HorIm;
			double Ver = VerRe*VerRe + VerIm*VerIm;
			s0 = Hor + Ver;
			s1 = Hor - Ver;
			s2 = -2.*(HorRe*VerRe + HorIm*VerIm);
			s3 = 2.*(-HorRe*VerIm + HorIm*VerRe);

			complex<double> BufS0S1(s0, s1);
			*(LocHorPolTravers++) = BufS0S1;
			complex<double> BufS2S3(s2, s3);
			*(LocVerPolTravers++) = BufS2S3;
		}
	}
}

//*************************************************************************

inline void srTRadInt::AnglesFromMrToRad()
{
	DistrInfoDat.xStart *= 1.E-03;
	DistrInfoDat.xEnd *= 1.E-03;
	DistrInfoDat.zStart *= 1.E-03;
	DistrInfoDat.zEnd *= 1.E-03;
}

//*************************************************************************

inline void srTRadInt::AnglesFromRadToMr()
{
	DistrInfoDat.xStart *= 1.E+03;
	DistrInfoDat.xEnd *= 1.E+03;
	DistrInfoDat.zStart *= 1.E+03;
	DistrInfoDat.zEnd *= 1.E+03;
}

//*************************************************************************

inline double srTRadInt::PhaseFun(double sArg, int DerOrder)
{
	const double PI = 3.141592653590;
	//const double TwoPI = 2.*PI;
	const double PIm10e6 = PI*1.E+06;
	const double PIm10e6dEnCon = PIm10e6*0.80654658;
	//const double TenEm6dTwoPI = 1.E-06/TwoPI;

	double dBxds=0., dBzds=0., Bx=0., Bz=0., Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('x', sArg, dBzds, Bz, Btx, Crdx, IntBtE2x);
	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('z', sArg, dBxds, Bx, Btz, Crdz, IntBtE2z);

	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		double One_d_ymis = 1./(ObsCoor.y - sArg);
		double xObs_mi_x = xObs - Crdx, zObs_mi_z = zObs - Crdz;
		double Nx = xObs_mi_x*One_d_ymis, Btx_mi_Nx = Btx - Nx;
		double Nz = zObs_mi_z*One_d_ymis, Btz_mi_Nz = Btz - Nz;

		if(DerOrder == 0) return PIm10e9_d_Lamb*(sArg*TrjDatPtr->EbmDat.GammaEm2 + IntBtE2x + IntBtE2z + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		else if(DerOrder == 1) 
			return PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
		else if(DerOrder == 2)
		{
			double ConBtxBzpAx = (TrjDatPtr->BetaNormConst)*Bz + Btx_mi_Nx*One_d_ymis;
			double ConBtzBxpAz = (-TrjDatPtr->BetaNormConst)*Bx + Btz_mi_Nz*One_d_ymis;
			return (2.*PIm10e9_d_Lamb)*(Btx_mi_Nx*ConBtxBzpAx + Btz_mi_Nz*ConBtzBxpAz);
		}
		else return 0.;
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		double Ax = Btx - xObs, Az = Btz - zObs;
		if(DerOrder == 0) return PIm10e9_d_Lamb*(sArg*(TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs) + IntBtE2x + IntBtE2z - 2.*(xObs*Crdx + zObs*Crdz));
		if(DerOrder == 1) return PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Ax*Ax + Az*Az);
		if(DerOrder == 2) return (2.*PIm10e9_d_Lamb)*(Ax*(TrjDatPtr->BetaNormConst)*Bz + Az*(-TrjDatPtr->BetaNormConst)*Bx);
	}
	return 0.;
}

//*************************************************************************

inline int srTRadInt::ScanPhase()
{
	if(AuxPhaseArray != 0) { delete[] AuxPhaseArray; AuxPhaseArray = 0;}
	int AmOfP = int((sIntegFin - sIntegStart)/sIntegStep) + 1;
	AuxPhaseArray = new double[AmOfP];
	//if(AuxPhaseArray == 0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
	if(AuxPhaseArray == 0) { return MEMORY_ALLOCATION_FAILURE;}

	ObsCoor.Lamb = DistrInfoDat.LambStart; 
	ObsCoor.x = DistrInfoDat.xStart; 
	ObsCoor.y = DistrInfoDat.yStart; 
	ObsCoor.z = DistrInfoDat.zStart; 

	double sArg = sIntegStart;
	double* tAuxPhaseArray = AuxPhaseArray;
	for(int i=0; i<AmOfP; i++)
	{
		//if(sArg > -500.)
		//{
		//	int StopHere = 1;
		//}

		double Phase = PhaseFun(sArg, DistrInfoDat.PhaseDerOrder);
		*(tAuxPhaseArray++) = Phase; sArg += sIntegStep;
	}
	return 0;
}

//*************************************************************************

inline int srTRadInt::SetupParamForManualIntegr()
{
	if(sIntegMethod < 2)
	{
		int BufInt = int((sIntegFin - sIntegStart)/sIntegStep) >> 1;
		AmOfPointsForManIntegr = (BufInt << 1) + 1;
	}
	else if(sIntegMethod == 2)
	{
		int BufInt = int((sIntegFin - sIntegStart)/(sIntegStep*3.));
		AmOfPointsForManIntegr = BufInt*3 + 1;
	}
	sIntegStep = (sIntegFin - sIntegStart)/(AmOfPointsForManIntegr - 1);

	if(!UseManualSlower)
	{
		if(BtxArr!=0) { delete[] BtxArr; BtxArr = 0;}
		//BtxArr = new double[AmOfPointsForManIntegr]; if(BtxArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		BtxArr = new double[AmOfPointsForManIntegr]; if(BtxArr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(XArr!=0) { delete[] XArr; XArr = 0;}
		//XArr = new double[AmOfPointsForManIntegr]; if(XArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		XArr = new double[AmOfPointsForManIntegr]; if(XArr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(IntBtxE2Arr!=0) { delete[] IntBtxE2Arr; IntBtxE2Arr = 0;}
		//IntBtxE2Arr = new double[AmOfPointsForManIntegr]; if(IntBtxE2Arr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		IntBtxE2Arr = new double[AmOfPointsForManIntegr]; if(IntBtxE2Arr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(BxArr!=0) { delete[] BxArr; BxArr = 0;}
		//BxArr = new double[AmOfPointsForManIntegr]; if(BxArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		BxArr = new double[AmOfPointsForManIntegr]; if(BxArr==0) { return MEMORY_ALLOCATION_FAILURE;}

		if(BtzArr!=0) { delete[] BtzArr; BtzArr = 0;}
		//BtzArr = new double[AmOfPointsForManIntegr]; if(BtzArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		BtzArr = new double[AmOfPointsForManIntegr]; if(BtzArr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(ZArr!=0) { delete[] ZArr; ZArr = 0;}
		//ZArr = new double[AmOfPointsForManIntegr]; if(ZArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		ZArr = new double[AmOfPointsForManIntegr]; if(ZArr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(IntBtzE2Arr!=0) { delete[] IntBtzE2Arr; IntBtzE2Arr = 0;}
		//IntBtzE2Arr = new double[AmOfPointsForManIntegr]; if(IntBtzE2Arr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		IntBtzE2Arr = new double[AmOfPointsForManIntegr]; if(IntBtzE2Arr==0) { return MEMORY_ALLOCATION_FAILURE;}
		
		if(BzArr!=0) { delete[] BzArr; BzArr = 0;}
		//BzArr = new double[AmOfPointsForManIntegr]; if(BzArr==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		BzArr = new double[AmOfPointsForManIntegr]; if(BzArr==0) { return MEMORY_ALLOCATION_FAILURE;}
	}
	return 0;
}

//*************************************************************************

inline void srTRadInt::CosAndSin(double x, double& Cos, double& Sin)
{
	x -= TwoPI*int(x*One_dTwoPI);
	if(x < 0.) x += TwoPI;

	char ChangeSign=0;
	if(x > ThreePIdTwo) x -= TwoPI;
	else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

	double xe2 = x*x;
	Cos = 1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c))));
	Sin = x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s)))));
	if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
}

//*************************************************************************

inline void srTRadInt::FunForRadInt(double sArg, complex<double>* FunPtr)
{
	double Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataDerivedAtPoint(sArg, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z);
	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	complex<double>* FunTravers = FunPtr;
	double CosPhase, SinPhase;

// Angles assumed in mr!
	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		double One_d_ymis = 1./(ObsCoor.y - sArg);
		double xObs_mi_x = xObs - Crdx, zObs_mi_z = zObs - Crdz;
		double Phase = PIm10e9_d_Lamb*(sArg*TrjDatPtr->EbmDat.GammaEm2 + IntBtE2x + IntBtE2z + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		CosAndSin(Phase, CosPhase, SinPhase);

		double Nx = xObs_mi_x*One_d_ymis, BufPreExpReX = One_d_ymis*(Btx - Nx);
		complex<double> BufTotComplX(BufPreExpReX*CosPhase, BufPreExpReX*SinPhase);
		*(FunTravers++) = BufTotComplX;
		double Nz = zObs_mi_z*One_d_ymis, BufPreExpReZ = One_d_ymis*(Btz - Nz);
		complex<double> BufTotComplZ(BufPreExpReZ*CosPhase, BufPreExpReZ*SinPhase);
		*FunTravers = BufTotComplZ;
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		double Btx_mi_HorAng = Btx - xObs, Btz_mi_HorAng = Btz - zObs;
		double Phase = PIm10e9_d_Lamb*(sArg*(TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs) + IntBtE2x + IntBtE2z - 2.*(xObs*Crdx + zObs*Crdz));
		CosAndSin(Phase, CosPhase, SinPhase);

		complex<double> BufTotComplX(Btx_mi_HorAng*CosPhase, Btx_mi_HorAng*SinPhase);
		*(FunTravers++) = BufTotComplX;
		complex<double> BufTotComplZ(Btz_mi_HorAng*CosPhase, Btz_mi_HorAng*SinPhase);
		*FunTravers = BufTotComplZ;
	}
}

//*************************************************************************

inline void srTRadInt::AxAzPhNearField(double sArg, double& Ax, double& Az, double& Ph)
{
	double Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataDerivedAtPoint(sArg, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z);
	double xObs = ObsCoor.x, zObs = ObsCoor.z;

	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	double One_d_ymis = 1./(ObsCoor.y - sArg);
	double xObs_mi_x = xObs - Crdx, zObs_mi_z = zObs - Crdz;
	Ph = PIm10e9_d_Lamb*(sArg*TrjDatPtr->EbmDat.GammaEm2 + IntBtE2x + IntBtE2z + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
	Ax = (Btx - xObs_mi_x*One_d_ymis)*One_d_ymis;
	Az = (Btz - zObs_mi_z*One_d_ymis)*One_d_ymis;
}

//*************************************************************************

inline void srTRadInt::AxAzPhFarField(double sArg, double& Ax, double& Az, double& Ph)
{
	double Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataDerivedAtPoint(sArg, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z);
	double xObs = ObsCoor.x, zObs = ObsCoor.z;

	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	
	double AngPhConst = TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs;
	double Two_xObs = 2.*xObs, Two_zObs = 2.*zObs;
	Ph = PIm10e9_d_Lamb*(sArg*AngPhConst + IntBtE2x + IntBtE2z - (Two_xObs*Crdx + Two_zObs*Crdz));
	Ax = Btx - xObs;
	Az = Btz - zObs;
}

//*************************************************************************

inline void srTRadInt::AxAzPhNearField2(int LevelNo, int IndxOnLevel, double s, double& Ax, double& Az, double& Ph)
{
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	double *pBtx = BtxArrP[LevelNo], *pBtz=BtzArrP[LevelNo], *pX=XArrP[LevelNo], *pZ=ZArrP[LevelNo], *pIntBtxE2=IntBtxE2ArrP[LevelNo], *pIntBtzE2=IntBtzE2ArrP[LevelNo];
	double One_d_ymis = 1./(ObsCoor.y - s);
	double xObs_mi_x = ObsCoor.x - *(pX+IndxOnLevel);
	double zObs_mi_z = ObsCoor.z - *(pZ+IndxOnLevel);
	double Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;

	Ph = PIm10e9_d_Lamb*(s*TrjDatPtr->EbmDat.GammaEm2 + *(pIntBtxE2+IndxOnLevel) + *(pIntBtzE2+IndxOnLevel) + xObs_mi_x*Nx + zObs_mi_z*Nz);
	Ax = (*(pBtx+IndxOnLevel) - Nx)*One_d_ymis;
	Az = (*(pBtz+IndxOnLevel) - Nz)*One_d_ymis;
}

//*************************************************************************

inline void srTRadInt::AxAzPhFarField2(int LevelNo, int IndxOnLevel, double s, double& Ax, double& Az, double& Ph)
{
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	double *pBtx = BtxArrP[LevelNo], *pBtz=BtzArrP[LevelNo], *pX=XArrP[LevelNo], *pZ=ZArrP[LevelNo], *pIntBtxE2=IntBtxE2ArrP[LevelNo], *pIntBtzE2=IntBtzE2ArrP[LevelNo];
	double xObs = ObsCoor.x, zObs = ObsCoor.z;

	double AngPhConst = TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs;
	double Two_xObs = 2.*xObs, Two_zObs = 2.*zObs;
	Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2+IndxOnLevel) + *(pIntBtzE2+IndxOnLevel) - (Two_xObs*(*(pX+IndxOnLevel)) + Two_zObs*(*(pZ+IndxOnLevel))));
	Ax = *(pBtx+IndxOnLevel) - xObs;
	Az = *(pBtz+IndxOnLevel) - zObs;
}

//*************************************************************************

//inline int srTRadInt::FillNextLevel(int LevelNo, double sStart, double sEnd, long Np)
inline int srTRadInt::FillNextLevel(int LevelNo, double sStart, double sEnd, long long Np)
{
	double* BasePtr = new double[Np*8];
	//if(BasePtr == 0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
	if(BasePtr == 0) { return MEMORY_ALLOCATION_FAILURE;}

	BtxArrP[LevelNo] = BasePtr;
	BasePtr += Np; XArrP[LevelNo] = BasePtr;
	BasePtr += Np; IntBtxE2ArrP[LevelNo] = BasePtr;
	BasePtr += Np; BxArrP[LevelNo] = BasePtr;

	BasePtr += Np; BtzArrP[LevelNo] = BasePtr;
	BasePtr += Np; ZArrP[LevelNo] = BasePtr;
	BasePtr += Np; IntBtzE2ArrP[LevelNo] = BasePtr;
	BasePtr += Np; BzArrP[LevelNo] = BasePtr;

	TrjDatPtr->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevelNo], BtzArrP[LevelNo], XArrP[LevelNo], ZArrP[LevelNo], IntBtxE2ArrP[LevelNo], IntBtzE2ArrP[LevelNo], BxArrP[LevelNo], BzArrP[LevelNo]);

	AmOfPointsOnLevel[LevelNo] = Np;
	NumberOfLevelsFilled++;
	return 0;
}

//*************************************************************************

inline int srTRadInt::SetupRadCompStructures()
{
	int result;
	if(TrjDataContShouldBeRebuild)
	{
		EvaluateMemAvailAfterTrjComp();
		if(sIntegMethod < 10) 
		{
			if(result = SetupParamForManualIntegr()) return result;
			if(!UseManualSlower) 
			{
				TrjDatPtr->CompTotalTrjData(sIntegStart, sIntegFin, AmOfPointsForManIntegr, BtxArr, BtzArr, XArr, ZArr, IntBtxE2Arr, IntBtzE2Arr, BxArr, BzArr);
			}
		}
		if(sIntegMethod == 10)
		{
			MaxFluxDensVal = CurrentAbsPrec = 0.;
		}
	}

	if(DistrInfoDat.RadDistrDataContShouldBeRebuild)
	{
		if(result = AllocateMemForRadDistr()) return result;
	}

	SetupNormalizingConst();
	return 0;
}

//*************************************************************************

inline int srTRadInt::SetupRadCompStructMethAuto2()
{// Old Auto2 - int by pieces - not used currently
	int NumberOfParts = 5; // To steer!

	double SubIntervLength = (sIntegFin - sIntegStart)/NumberOfParts;
	double sStart = sIntegStart;

	for(int k=0; k<NumberOfParts; k++)
	{
		double sFin = sStart + SubIntervLength;
		srTPartAutoRadInt* pPartAutoRadInt = new srTPartAutoRadInt(TrjDatPtr, &DistrInfoDat, &ObsCoor, sStart, sFin, sIntegRelPrec, 5);
		//if(pPartAutoRadInt == 0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		if(pPartAutoRadInt == 0) { return MEMORY_ALLOCATION_FAILURE;}

		srTPartAutoRadIntHndl PartAutoRadIntHndl(pPartAutoRadInt);
		PartAutoRadIntHndlVect.push_back(PartAutoRadIntHndl);
		sStart = sFin;
	}
	return 0;
}

//*************************************************************************

inline void srTRadInt::EnableNormDerComp(TVector3d& InSurfNorm)
{
	ComputeNormalDerivative = 1; // For Diffraction computation in Near Field
	SurfNorm = InSurfNorm;
}

//*************************************************************************

inline void srTRadInt::DisableNormDerComp()
{
	ComputeNormalDerivative = 0; // For Diffraction computation in Near Field
}

//*************************************************************************

inline int srTRadInt::EvaluateMemAvailAfterTrjComp()
{
	if(sIntegMethod < 10)
	{
		int BufInt = int((sIntegFin - sIntegStart)/sIntegStep) >> 1;
		int LocAmOfPointsForManIntegr = (BufInt << 1) + 1;

		if(MaxNumPoToSave <= LocAmOfPointsForManIntegr) UseManualSlower = 1;
		else UseManualSlower = 0;
	}
	else if((sIntegMethod == 10) || (sIntegMethod == 11))
	{
		const int NpOnZeroLev = 5; // Change here if it is changed in Auto methods
		double BufDouble = (double)((MaxNumPoToSave - 1)/(NpOnZeroLev - 1));

		MaxLevelForMeth_10_11 = int(log(BufDouble)*1.443);
	}
	return (CurrMemAvail > 0.)? 0 : -1;
}

//*************************************************************************

inline void srTRadInt::EstimateAbsoluteTolerance()
{//To check this !!! 

	double GammaE2 = (TrjDatPtr->EbmDat.Gamma)*(TrjDatPtr->EbmDat.Gamma);
	double AvgObsDist = (DistrInfoDat.yStart + DistrInfoDat.yEnd)*0.5;
	double AvgObsWaveLength = (DistrInfoDat.LambStart + DistrInfoDat.LambEnd)*0.5;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*AvgObsWaveLength : PIm10e6*1000./AvgObsWaveLength;

	double c10e9_d_Lamb = PIm10e9_d_Lamb/PI;
	double InvDiffractAngleE2 = AvgObsDist*c10e9_d_Lamb;
	double Min_GammaE2_InvDiffractAngleE2 = (GammaE2 < InvDiffractAngleE2)? GammaE2 : InvDiffractAngleE2;
	double InvAvgObsDist = 1./AvgObsDist;
	double ExpectedIntesityValues = (3.47E+12)*(TrjDatPtr->EbmDat.Current)*Min_GammaE2_InvDiffractAngleE2*InvAvgObsDist*InvAvgObsDist;

	EstimatedAbsoluteTolerance = (1.E-05)*ExpectedIntesityValues*sIntegRelPrec;
}

//*************************************************************************

inline void srTRadInt::CheckFurtherSubdNeed(srTEFourierVect* pEwVect, long StartInd, long LenSideX, char* xSubdNeedArr, char* zSubdNeedArr)
{
	srTEFourier* EwPtrsArr[5];
	int SectNo = StartInd;
	for(int iz=0; iz<5; iz++)
	{
		for(int ix=0; ix<5; ix++)
		{
			EwPtrsArr[ix] = &((*pEwVect)[SectNo+ix]);
		}
		xSubdNeedArr[iz] = CheckFurtherSubdNeedForOneCoord(EwPtrsArr);
		SectNo += LenSideX;
	}
	SectNo = StartInd;
	for(int ix=0; ix<5; ix++)
	{
		int BufNo = SectNo;
		for(int iz=0; iz<5; iz++)
		{
			EwPtrsArr[iz] = &((*pEwVect)[BufNo]);
			BufNo += LenSideX;
		}
		zSubdNeedArr[ix] = CheckFurtherSubdNeedForOneCoord(EwPtrsArr);
		SectNo++;
	}
}

//*************************************************************************

inline void srTRadInt::CopySymEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, char SymWithRespectToXax, char SymWithRespectToZax)
{
	char ChangeSignEx = SymWithRespectToZax;
	char ChangeSignEz = SymWithRespectToXax;

	float *tOrigEx = pOrigDataEx, *tSymEx = pSymDataEx;
	float *tOrigEz = pOrigDataEz, *tSymEz = pSymDataEz;
	for(int ie=0; ie<DistrInfoDat.nLamb; ie++)
	{
		*tSymEx = *(tOrigEx++); *(tSymEx + 1) = *(tOrigEx++);
		if(ChangeSignEx) { *tSymEx = -(*tSymEx); *(tSymEx + 1) *= -1;}
		tSymEx += 2;

		*tSymEz = *(tOrigEz++); *(tSymEz + 1) = *(tOrigEz++);
		if(ChangeSignEz) { *tSymEz = -(*tSymEz); *(tSymEz + 1) *= -1;}
		tSymEz += 2;
	}
}

//*************************************************************************

inline srTEFourier srTRadInt::RadFuncInDriftSpace(double s, double sSt, double Btx, double xSt, double Btz, double zSt, double IntBtE2xzSt, double& Phase)
{
	const double CritRatioForRootArg = 1.E-03;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double sRel = s - sSt;
	double x = xSt + sRel*Btx, z = zSt + sRel*Btz;
	double dx = ObsCoor.x - x, dz = ObsCoor.z - z;
	double dyLarge = ObsCoor.y - s, dySmall = 0.5*(IntBtE2xzSt + (Btx*Btx + Btz*Btz)*(s - sSt));
	double RootArgSmall = dx*dx + dz*dz + dySmall*dySmall + 2*dySmall*dyLarge;
	double RootArgLarge = dyLarge*dyLarge;
	double RootArgRatio = RootArgSmall/RootArgLarge;
	double Ph, InvR, CosPh, SinPh;
	if((RootArgRatio > CritRatioForRootArg) || (dyLarge < 0.))
	{
		double R = sqrt(RootArgSmall + RootArgLarge);
		Ph = (2*PIm10e9_d_Lamb)*(R - dyLarge + 0.5*(TrjDatPtr->EbmDat.GammaEm2)*s);
		InvR = 1./R;
	}
	else
	{
		Ph = PIm10e9_d_Lamb*((TrjDatPtr->EbmDat.GammaEm2)*s + dyLarge*RootArgRatio*(1 - 0.25*RootArgRatio*(1 - 0.5*RootArgRatio)));
		InvR = 1./dyLarge;
	}
	double Nx = dx*InvR, Nz = dz*InvR, InvRe2_d_k = InvR*InvR/(2*PIm10e9_d_Lamb);
	double AxRe = (Btx - Nx)*InvR, AxIm = -Nx*InvRe2_d_k;
	double AzRe = (Btz - Nz)*InvR, AzIm = -Nz*InvRe2_d_k;
	CosAndSin(Ph, CosPh, SinPh);
	Phase = Ph;
	srTEFourier Ew(AxRe*CosPh - AxIm*SinPh, AxRe*SinPh + AxIm*CosPh, AzRe*CosPh - AzIm*SinPh, AzRe*SinPh + AzIm*CosPh);
	return 	Ew;
}

//*************************************************************************

#endif
