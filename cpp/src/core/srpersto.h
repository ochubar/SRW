/************************************************************************//**
 * File: srpersto.h
 * Description: Calculation of Stokes parameters / Spectral Flux of Undulator Radiation from a Finite-Emittance Electron Beam, collected by a finite-size rectangular aperture (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRPERSTO_H
#define __SRPERSTO_H

#include "srptrjdt.h"
#include "srebmdat.h"
#include "gmfft.h"
#include "srmlttsk.h"
#include "srerror.h"

//*************************************************************************

extern srTYield srYield;

//*************************************************************************

struct srTRadIntPerStoPrec {
	int InitHarm, FinHarm;
	double Kns, Knphi;
	char IntensityOrFlux;
	double MinPhotEnExtRight; //OC170713

	//srTRadIntPerStoPrec(int InInitHarm, int InFinHarm, double InKns, double InKnphi, char InIntensityOrFlux)
	//{
    //       InitHarm = InInitHarm; FinHarm = InFinHarm;
    //       Kns = InKns; Knphi = InKnphi;
    //       IntensityOrFlux = InIntensityOrFlux;
	//}
};

//*************************************************************************

struct srTEnergyAdjustFinGridPar {
	double keStep, eShift;
	long AmOfExtraStepsFromLeftAndRightOnTheOtherSide;
	int Multip;
	long AmStExtraLeft, AmStExtraRight;
	int AmOfInitPoToIgnore;

	srTEnergyAdjustFinGridPar()
	{
		eShift = 0.;
		AmOfExtraStepsFromLeftAndRightOnTheOtherSide = 0;
		Multip = 1;
		AmStExtraLeft = 0;
		AmStExtraRight = 0;
		AmOfInitPoToIgnore = 0;
	}
};

//*************************************************************************

struct srTEnergyAzimuthGrid {
	long Ne;
	double eStart, eFin;
	double eMinEff, eMaxEff;
	double eCritForHarm;
	long NeExtraLeft, NeExtraRight;

	srTEnergyAdjustFinGridPar EnergyAdjustFinGridPar;

	double PhiMin, PhiMax;
	int* AmOfAzPoints;
	char AmsOfPointsOverPhiAreConstant;

	double eStepToResolve, PhiLenToResolve;

	double eStepToResolveObsPixels, eStepToResolveFinNperAndEnSpr;
	char EnsureEnResolvingObsPixels;
	srTEFourier IntOfInfNperData;

	double *CosLookUpArray, *SinLookUpArray;
	srTCosAndSinComp CosAndSinComp;

	srTEnergyAzimuthGrid(double& eStartIn, double& eFinIn, long& NeIn)
	{
		eStart = eStartIn; eFin = eFinIn; Ne = NeIn;
		AmOfAzPoints = 0;
		if(Ne > 0)
		{
			AmOfAzPoints = new int[Ne];
			MakeAllAzPointsConstant(1);
		}
		PhiMin = 0.; PhiMax = 6.28318530717959;
		AmsOfPointsOverPhiAreConstant = 0;
		CosLookUpArray = SinLookUpArray = 0;

		EnsureEnResolvingObsPixels = 1;
	}
	srTEnergyAzimuthGrid()
	{
		AmOfAzPoints = 0;
		PhiMin = 0.; PhiMax = 6.28318530717959;
		AmsOfPointsOverPhiAreConstant = 0;
		CosLookUpArray = SinLookUpArray = 0;

		EnsureEnResolvingObsPixels = 1;
	}
	~srTEnergyAzimuthGrid()
	{
		if(AmOfAzPoints != 0) delete[] AmOfAzPoints;

		if(CosLookUpArray != 0) delete[] CosLookUpArray;
		if(SinLookUpArray != 0) delete[] SinLookUpArray;
	}

	int AllocateAmOfPointsArray()
	{
		if(Ne == 0) return 0;
		if(AmOfAzPoints != 0)
		{
			delete[] AmOfAzPoints; AmOfAzPoints = 0;
		}
		AmOfAzPoints = new int[Ne];
		if(AmOfAzPoints == 0) return MEMORY_ALLOCATION_FAILURE;
		MakeAllAzPointsConstant(1);
		return 0;
	}

	void MakeAllAzPointsConstant(int ConstAmOfPoints)
	{
		if(AmOfAzPoints == 0) return;
		int *t = AmOfAzPoints;
		for(long i=0; i<Ne; i++) *(t++) = ConstAmOfPoints;
		AmsOfPointsOverPhiAreConstant = 1;
	}

	int SetUpCosAndSinLookUpArrays();
};

//*************************************************************************

class srTRadIntPeriodic {

	double HalfPI, PI, TwoPI, ThreePIdTwo, One_dTwoPI, One_dHalfPI; // Constants
	double a2c, a4c, a6c, a8c, a10c, a12c, a3s, a5s, a7s, a9s, a11s, a13s;

	double One_d_SqrtPi, Two_d_SqrtPi;

	double txE2ptzE2pgEm2_Fs, InvLambn_Fs, PiInvLambn_Fs, tx_Fs, tz_Fs, Twotx_Fs, Twotz_Fs;
	
	double dxHalfSlit, dzHalfSlit, Inv_SigmaXpSqrt2, Inv_SigmaZpSqrt2;
	double xPartFactRangeGen, zPartFactRangeGen;
	double SigFactGen;

	double Inv_dsSlit_mmEm2;
	double SigmaXpGen, SigmaZpGen;

	double NperGen, Con_PhPerSecPer10em3bw;
	char EnergySpreadShouldBeTreated;

	double EnAvg_FinNper, eExtra_FinNper, GridMisfit_FinNper;

	char s0ShouldBeTreatedGen, s1ShouldBeTreatedGen, s2ShouldBeTreatedGen, s3ShouldBeTreatedGen;
	double PartDistrNormGen;

	double Tet_Fphi, n_Fphi;
	CGenMathAuxDataForSharpEdgeCorr1D AuxDataForSharpEdgeCorrGen;

	//float* pA_Fphi;
	double* pA_Fphi; //OC020112
	double txMin_Fphi, txMax_Fphi, tzMin_Fphi, tzMax_Fphi;
	srTEnergyAzimuthGrid* pEnAzGrid_Fphi;

	int NsGen, NphiGen;

	srTEXZ EXZ;
	double *BtxArr, *BtzArr, *XArr, *ZArr, *IntBtE2Arr;

public:

	srTEbmDat EbmDat;
	srTMagFieldPeriodic MagPer;
	srTWfrSmp DistrInfoDat;
	srTRadIntPerStoPrec IntPerStoPrec;

	srTIntVect* pWarningsGen;

	srTRadIntPeriodic()
	{
		Initialize();
	}
	srTRadIntPeriodic(srTEbmDat* pElecBeam, srTMagFieldPeriodic* pMagFldPer, srTWfrSmp* pWfrSmp, void* pPrcPar);
	~srTRadIntPeriodic()
	{
		DisposeFieldBasedArrays();
	}

	void Initialize();
    static void ComputeStokes(srTEbmDat* pElecBeam, srTMagFieldPeriodic* pMagFldPer, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes);

	int CheckInputConsistency();
	//int ComputeTotalStokesDistr(srTStokesStructAccessData&);
	int ComputeTotalStokesDistr(srTStokesStructAccessData* pStokesAccessData, SRWLStructStokes* pStokesSRWL=0);

	int DeduceGridOverPhotonEnergyAndAzimuth(int n, double& eStart, double& eFin, long& ne, srTEnergyAzimuthGrid& EnAzGrid);
	void CorrectGridForPassingThroughCritEnergy(int n, double& eStart, double& eStep, long& ne);
	void CorrectGridToAllowRangeResizeOnTheOtherSide(srTEnergyAzimuthGrid& EnAzGrid);
	void CorrectGridForOnePoint(srTEnergyAzimuthGrid& EnAzGrid);
	int SetUpVariableGridOverAzimuth(int n, srTEnergyAzimuthGrid& EnAzGrid);
	void FindAngularObsGrid(double& xAngStart, double& xAngStep, double& zAngStart, double& zAngStep);

	//int ComputeLongIntForEnAndAz(int, srTEnergyAzimuthGrid&, float**&, int**&);
	int ComputeLongIntForEnAndAz(int, srTEnergyAzimuthGrid&, double**&, int**&); //OC020112
	//int RestoreLongIntArray(long ie, srTEnergyAzimuthGrid& EnAzGrid, float** LongIntArrays, int** LongIntArrInfo, float*& pArr);
	int RestoreLongIntArray(long ie, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, double*& pArr); //OC020112
	int AllocateLongIntArraysForEnAndAz(srTEnergyAzimuthGrid&, float**&);
	//void DisposeLongIntArraysForEnAndAz(srTEnergyAzimuthGrid&, float**&, int**&);
	void DisposeLongIntArraysForEnAndAz(srTEnergyAzimuthGrid&, double**&, int**&); //OC020112
	void FindObservationLimits(double&, double&, double&, double&);
	//void FindMostOffsetPixelCenter(double&, double&);
	void FindLeastAndMostOffsetPixelCenters(double& tx0, double& tz0, double& tx1, double& tz1);
	void FindTetMinMaxE2_FromTetxTetz(double txMin, double txMax, double tzMin, double tzMax, double& TetMinE2, double& TetMaxE2);

	void EstimateEnergyStepAndPhiLenToResolveObsPixels(int, double&, double&);
	void EstimateEnergyStepToResolveFinNper(int, srTEnergyAzimuthGrid&);
	void EstimateEnergyAndPhiObsLimits(int n, double& eMinEff, double& eMaxEff, double& PhiMinEff, double& PhiMaxEff);
	void FindPhiIntervalForVectors(TVector2d* Vectors, int LenVectors, double& PhiMin, double& PhiMax);
	double PhiIntToResolveBox(double x1, double x2, double y1, double y2, double rAvg);

	//int ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, float** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice);
	//int ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt); //OC020112
	int ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt); //OC020112
	void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData*, SRWLStructStokes*); //080612
	//void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData&);

	//int TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice);
	//int TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt); //OC020112
	int TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt); //OC020112
	//int FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice);
	//int FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt); //OC020112
	int FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt); //OC020112

	//int ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice);
	//int ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt); //OC020112
	int ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt); //OC020112
	void FindIntegralOfInfNperData(int n, srTEnergyAzimuthGrid& EnAzGrid, float* InfNperHarmData);
	double EstimateTaperResCurveWidth(int n);

	int Int1D_Simpson(double xSt, double xFi, long Nx, char VsSorPhi, srTEFourier&);

	int SetupConvolutionData_Normal(int n, float* ConvFactorData, double eStart, double eFin, long Ne);
	int SetupConvolutionData_Tapered(int n, float* ConvFactorData, double eStart, double eFin, long Ne);
	int SetupConvolutionData_OpticalKlystron(int n, float* ConvFactorData, double eStart, double eFin, long Ne);
	int SetupConvolutionData(int n, float* ConvFactorData, double eStart, double eFin, long Ne)
	{
		if(MagPer.TypeOfUnd == 1) return SetupConvolutionData_Normal(n, ConvFactorData, eStart, eFin, Ne);
		else if(MagPer.TypeOfUnd == 2) return SetupConvolutionData_Tapered(n, ConvFactorData, eStart, eFin, Ne);
		else if(MagPer.TypeOfUnd == 3) return SetupConvolutionData_OpticalKlystron(n, ConvFactorData, eStart, eFin, Ne);
		else return 0;
	}

	void AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
	{
		FinalResAreSymOverX = FinalResAreSymOverZ = 0;
		if(MagPer.FieldSymmetryInd != 0) return;
		if(MagPer.AmOfHarm > 1) return;

		if(DistrInfoDat.nx > 1)
		{
			double xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
			double xTol = xStep*0.001; // To steer
			FinalResAreSymOverX = (::fabs(DistrInfoDat.xStart + DistrInfoDat.xEnd - EbmDat.dxds0) < xTol)? 1 : 0;
		}
		if(DistrInfoDat.nz > 1)
		{
			double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
			double zTol = zStep*0.001; // To steer
			FinalResAreSymOverZ = (::fabs(DistrInfoDat.zStart + DistrInfoDat.zEnd - EbmDat.dzds0) < zTol)? 1 : 0;
		}
	}
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
	//void CopySymEnergySliceSRWL(SRWLStructStokes& stk, long ofstOrigData, long ofstSymData, char SymWithRespectToXax, char SymWithRespectToZax)
	void CopySymEnergySliceSRWL(SRWLStructStokes& stk, long long ofstOrigData, long long ofstSymData, char SymWithRespectToXax, char SymWithRespectToZax)
	{
		char ChangeSignS2 = !(SymWithRespectToXax && SymWithRespectToZax);
		char ChangeSignS3 = SymWithRespectToXax;
		float *tOrigS0 = ((float*)(stk.arS0)) + ofstOrigData, *tSymS0 = ((float*)(stk.arS0)) + ofstSymData;
		float *tOrigS1 = ((float*)(stk.arS1)) + ofstOrigData, *tSymS1 = ((float*)(stk.arS1)) + ofstSymData;
		float *tOrigS2 = ((float*)(stk.arS2)) + ofstOrigData, *tSymS2 = ((float*)(stk.arS2)) + ofstSymData;
		float *tOrigS3 = ((float*)(stk.arS3)) + ofstOrigData, *tSymS3 = ((float*)(stk.arS3)) + ofstSymData;
		for(int ie=0; ie<DistrInfoDat.nLamb; ie++)
		{
			*(tSymS0++) = *(tOrigS0++);
			*(tSymS1++) = *(tOrigS1++);
			*(tSymS2++) = ChangeSignS2? -(*(tOrigS2++)) : *(tOrigS2++);
			*(tSymS3++) = ChangeSignS3? -(*(tOrigS3++)) : *(tOrigS3++);
		}
	}

	int A(int n, double tx, double tz, srTEFourier& E)
	{// Returns Stokes
		int result;
		txE2ptzE2pgEm2_Fs = tx*tx + tz*tz + EbmDat.GammaEm2;
		InvLambn_Fs = (n << 1)/(MagPer.PerLength*(EbmDat.GammaEm2*MagPer.HalfKxE2pKzE2 + txE2ptzE2pgEm2_Fs));
		PiInvLambn_Fs = PI*InvLambn_Fs;
		tx_Fs = tx; tz_Fs = tz;
		Twotx_Fs = tx*2.; Twotz_Fs = tz*2.;

		if(((NsGen >> 1) << 1) == NsGen) NsGen++;
		if(result = Int1D_Simpson(0., MagPer.PerLength, NsGen, 's', E)) return result;
		E.EwX_Re *= InvLambn_Fs; E.EwX_Im *= InvLambn_Fs; E.EwZ_Re *= InvLambn_Fs; E.EwZ_Im *= InvLambn_Fs;

		double LinHor = E.EwX_Re*E.EwX_Re + E.EwX_Im*E.EwX_Im;
		double LinVer = E.EwZ_Re*E.EwZ_Re + E.EwZ_Im*E.EwZ_Im;
		double s0 = LinHor + LinVer, s1 = LinHor - LinVer;
		double s2 = -2.*(E.EwX_Re*E.EwZ_Re + E.EwX_Im*E.EwZ_Im);
		double s3 = 2.*(-E.EwX_Re*E.EwZ_Im + E.EwX_Im*E.EwZ_Re);
		E.EwX_Re = s0; E.EwX_Im = s1; E.EwZ_Re = s2; E.EwZ_Im = s3;
		return 0;
	}
	void Fs(double s, int is, srTEFourier& E)
	{
		double Phase = PiInvLambn_Fs*(txE2ptzE2pgEm2_Fs*s + IntBtE2Arr[is] - (Twotx_Fs*XArr[is] + Twotz_Fs*ZArr[is]));
		double CosPh, SinPh;
		CosAndSin(Phase, CosPh, SinPh);
		double Ax = BtxArr[is] - tx_Fs, Az = BtzArr[is] - tz_Fs;
		E.EwX_Re = Ax*CosPh; E.EwX_Im = Ax*SinPh;
		E.EwZ_Re = Az*CosPh; E.EwZ_Im = Az*SinPh;
	}

	void Fphi(double phi, int iphi, srTEFourier& Stokes)
	{
		double CosPhi, SinPhi;
		if(pEnAzGrid_Fphi->AmsOfPointsOverPhiAreConstant)
		{
			CosPhi = *(pEnAzGrid_Fphi->CosLookUpArray + iphi); SinPhi = *(pEnAzGrid_Fphi->SinLookUpArray + iphi);
		}
		else CosAndSin(phi, CosPhi, SinPhi);
		double TetCosPhi = Tet_Fphi*CosPhi, TetSinPhi = Tet_Fphi*SinPhi;

		double BufX = EXZ.x - EbmDat.dxds0 - TetCosPhi, BufZ = EXZ.z - EbmDat.dzds0 - TetSinPhi;
		if((BufX < -xPartFactRangeGen) || (BufX > xPartFactRangeGen) || (BufZ < -zPartFactRangeGen) || (BufZ > zPartFactRangeGen))
		{
			Stokes.EwX_Re = 0.; Stokes.EwX_Im = 0.; Stokes.EwZ_Re = 0.; Stokes.EwZ_Im = 0.; return;
		}
		double Fact = 0.25*(Erf(Inv_SigmaXpSqrt2*(BufX + dxHalfSlit)) - Erf(Inv_SigmaXpSqrt2*(BufX - dxHalfSlit)))
							*(Erf(Inv_SigmaZpSqrt2*(BufZ + dzHalfSlit)) - Erf(Inv_SigmaZpSqrt2*(BufZ - dzHalfSlit)));

		//float* tA = pA_Fphi + (iphi << 2);
		double* tA = pA_Fphi + (iphi << 2); //OC020112
		Stokes.EwX_Re = (*(tA++))*Fact;  Stokes.EwX_Im = (*(tA++))*Fact; 
		Stokes.EwZ_Re = (*(tA++))*Fact;  Stokes.EwZ_Im = (*(tA++))*Fact;
	}

	char FilamentTreatmentIsPossible(srTEnergyAzimuthGrid& EnAzGrid)
	{
		//const double FactorToAllowFilamentTreatment = 100.;
		//OC111110: =100.; was producing buggy on- and off-axis spectra through a small aperture
		const double FactorToAllowFilamentTreatment = 100000.; //-to revise this criterion
		double eResolRat = EnAzGrid.eStepToResolveFinNperAndEnSpr/(EnAzGrid.eStepToResolveObsPixels + 1.E-20);

		return ((FactorToAllowFilamentTreatment < eResolRat) && 
				(EnAzGrid.EnergyAdjustFinGridPar.Multip >= 2));
	}

	void FindPartDistrNorm()
	{
		double tx = EXZ.x - EbmDat.dxds0, tz = EXZ.z - EbmDat.dzds0;
		PartDistrNormGen = 1./F_PartDistr(tx, tz);
	}
	double F_PartDistr(double tx, double tz)
	{
		double BufX = EXZ.x - EbmDat.dxds0 - tx;
		double BufZ = EXZ.z - EbmDat.dzds0 - tz;
		if((BufX < -xPartFactRangeGen) || (BufX > xPartFactRangeGen) || (BufZ < -zPartFactRangeGen) || (BufZ > zPartFactRangeGen)) return 0.;

		return 0.25*(Erf(Inv_SigmaXpSqrt2*(BufX + dxHalfSlit)) - Erf(Inv_SigmaXpSqrt2*(BufX - dxHalfSlit)))
				   *(Erf(Inv_SigmaZpSqrt2*(BufZ + dzHalfSlit)) - Erf(Inv_SigmaZpSqrt2*(BufZ - dzHalfSlit)));
	}

	void SetUpAvgEnergy(int n)
	{
		const double c0 = 1.239854E-09;
		if(DistrInfoDat.nLamb == 1) EnAvg_FinNper = DistrInfoDat.LambStart; // "global" in this class
		else
		{
			double CentralAngleX = EXZ.x - EbmDat.dxds0, CentralAngleZ = EXZ.z - EbmDat.dzds0;
			double CentralAngleE2 = CentralAngleX*CentralAngleX + CentralAngleZ*CentralAngleZ;
			EnAvg_FinNper = c0*(n << 1)/(MagPer.PerLength*(EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2) + CentralAngleE2)); // "global" in this class
		}
	}

	void DeduceSlitAngDims()
	{
		if(DistrInfoDat.nx == 1)
		{
			dxHalfSlit = 0.5*(DistrInfoDat.xEnd - DistrInfoDat.xStart)/DistrInfoDat.yStart;
		}
		else
		{
			dxHalfSlit = 0.5*(DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.yStart*(DistrInfoDat.nx - 1));
		}
		if(DistrInfoDat.nz == 1)
		{
			dzHalfSlit = 0.5*(DistrInfoDat.zEnd - DistrInfoDat.zStart)/DistrInfoDat.yStart;
		}
		else
		{
			dzHalfSlit = 0.5*(DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.yStart*(DistrInfoDat.nz - 1));
		}
		Inv_dsSlit_mmEm2 = 0.25E-06/(dxHalfSlit*dzHalfSlit*DistrInfoDat.yStart*DistrInfoDat.yStart);
	}
	void DeduceE_BeamConstants()
	{// This should be called after DeduceSlitAngDims()
		double ye2 = DistrInfoDat.yStart*DistrInfoDat.yStart;
		double SigXpE2 = EbmDat.Mxpxp + EbmDat.Mxx/ye2 + 2.*EbmDat.Mxxp/DistrInfoDat.yStart;
		double SigZpE2 = EbmDat.Mzpzp + EbmDat.Mzz/ye2 + 2.*EbmDat.Mzzp/DistrInfoDat.yStart;

		SigmaXpGen = sqrt(SigXpE2); SigmaZpGen = sqrt(SigZpE2);
		Inv_SigmaXpSqrt2 = 0.707106781186547/SigmaXpGen;
		Inv_SigmaZpSqrt2 = 0.707106781186547/SigmaZpGen;

		xPartFactRangeGen = dxHalfSlit + SigFactGen*SigmaXpGen;
		zPartFactRangeGen = dzHalfSlit + SigFactGen*SigmaZpGen;
	}
	void DeduceNormConstants()
	{
		NperGen = int(MagPer.TotLength/MagPer.PerLength);
		if(NperGen < 7) // To steer
		{
			//srTSend Send; Send.AddWarningMessage(pWarningsGen, TOO_FEW_PERIODS);
			CErrWarn::AddWarningMessage(pWarningsGen, TOO_FEW_PERIODS);
		}

		Con_PhPerSecPer10em3bw = EbmDat.Current*NperGen/(137.0360411*1000.*1.602189246E-19*MagPer.PerLength);
		EnergySpreadShouldBeTreated = (MagPer.TypeOfUnd > 0); // Change if needed
	}

	void DeduceNsForOnePeriod(int n)
	{
		const int PointsPerFringe = 13; // To steer

		double dNs = n*PointsPerFringe*IntPerStoPrec.Kns;
		NsGen = int(dNs);
		if(dNs - NsGen > 0.01) NsGen++;
		if(((NsGen >> 1) << 1) == NsGen) NsGen++;
	}

	void DeduceNphi(int n)
	{
		const int PointsPerTwoPi_Even = 67; // To steer
		const int PointsPerTwoPi_Odd = 55; // To steer

		int PointsPerTwoPi = (((n >> 1) << 1) == n)? PointsPerTwoPi_Even : PointsPerTwoPi_Odd;
		double dNphiGen = PointsPerTwoPi*IntPerStoPrec.Knphi;
		NphiGen = int(dNphiGen);
		if(dNphiGen - NphiGen > 0.01) NphiGen++;
		if(((NphiGen >> 1) << 1) == NphiGen) NphiGen++;
	}
	char HarmonicIsNeeded(int n)
	{
		return ((n >= IntPerStoPrec.InitHarm) && (n <= IntPerStoPrec.FinHarm));
	}
	void AddHarmonicData(float* HarmonicData, float* TotData, long ne)
	{
		float *tIn = HarmonicData, *tOut = TotData;
		for(int ie=0; ie<ne; ie++)
		{
			char BadCaseNoticed = 0;
			float s0 = *(tIn++), s1 = *(tIn++), s2 = *(tIn++), s3 = *(tIn++);
			float LinHor = (float)(0.5*(s0 + s1)); if(LinHor < 0.) { LinHor = 0.; BadCaseNoticed = 1;}
			float LinVer = (float)(0.5*(s0 - s1)); if(LinVer < 0.) { LinVer = 0.; BadCaseNoticed = 1;}

			if(BadCaseNoticed)
			{
				s0 = LinHor + LinVer; s1 = LinHor - LinVer; s2 = s3 = 0.;
			}

			*(tOut++) += s0; *(tOut++) += s1; *(tOut++) += s2; *(tOut++) += s3;
		}
	}

	int AllocateFieldBasedArrays()
	{
		DisposeFieldBasedArrays();

		BtxArr = new double[NsGen*5];
		if(BtxArr == 0) return MEMORY_ALLOCATION_FAILURE;
		BtzArr = BtxArr + NsGen;
		XArr = BtzArr + NsGen;
		ZArr = XArr + NsGen;
		IntBtE2Arr = ZArr + NsGen;

		return 0;
	}
	void DisposeFieldBasedArrays()
	{
		if(BtxArr != 0) delete[] BtxArr;
		BtxArr = BtzArr = XArr = ZArr = IntBtE2Arr = 0;
	}
	//void ZeroOutData(srTStokesStructAccessData& StokesAccessData)
	void ZeroOutData(srTStokesStructAccessData* pStokesAccessData, SRWLStructStokes* pStokesSRWL)
	{
		long Ne, Nx, Nz;
		if(pStokesAccessData != 0)
		{
			Ne = pStokesAccessData->ne;
			Nx = pStokesAccessData->nx;
			Nz = pStokesAccessData->nz;
			float *tStokes = pStokesAccessData->pBaseSto;
			for(long iz=0; iz<Nz; iz++)
			{
				for(long ix=0; ix<Nx; ix++)
				{
					for(long ie=0; ie<Ne; ie++)
					{
						for(int i=0; i<4; i++) *(tStokes++) = 0.;
					}
				}
			}
		}
		else if(pStokesSRWL != 0)
		{
			//Ne = pStokesSRWL->ne;
			//Nx = pStokesSRWL->nx;
			//Nz = pStokesSRWL->ny;
			Ne = pStokesSRWL->mesh.ne;
			Nx = pStokesSRWL->mesh.nx;
			Nz = pStokesSRWL->mesh.ny;

			float *afS0 = 0, *afS1 = 0, *afS2 = 0, *afS3 = 0, *tfS0 = 0, *tfS1 = 0, *tfS2 = 0, *tfS3 = 0;
			double *adS0 = 0, *adS1 = 0, *adS2 = 0, *adS3 = 0, *tdS0 = 0, *tdS1 = 0, *tdS2 = 0, *tdS3 = 0;
			if(pStokesSRWL->numTypeStokes == 'f') 
			{
				if(pStokesSRWL->arS0 != 0) { afS0 = (float*)(pStokesSRWL->arS0); tfS0 = afS0;}
				if(pStokesSRWL->arS1 != 0) { afS1 = (float*)(pStokesSRWL->arS1); tfS1 = afS1;}
				if(pStokesSRWL->arS2 != 0) { afS2 = (float*)(pStokesSRWL->arS2); tfS2 = afS2;}
				if(pStokesSRWL->arS3 != 0) { afS3 = (float*)(pStokesSRWL->arS3); tfS3 = afS3;}
				for(long iz=0; iz<Nz; iz++)
				{
					for(long ix=0; ix<Nx; ix++)
					{
						for(long ie=0; ie<Ne; ie++)
						{
							if(afS0 != 0) *(tfS0++) = 0.;
							if(afS1 != 0) *(tfS1++) = 0.;
							if(afS2 != 0) *(tfS2++) = 0.;
							if(afS3 != 0) *(tfS3++) = 0.;
						}
					}
				}
			}
			else if(pStokesSRWL->numTypeStokes == 'd') 
			{
				if(pStokesSRWL->arS0 != 0) { adS0 = (double*)(pStokesSRWL->arS0); tdS0 = adS0;}
				if(pStokesSRWL->arS1 != 0) { adS1 = (double*)(pStokesSRWL->arS1); tdS1 = adS1;}
				if(pStokesSRWL->arS2 != 0) { adS2 = (double*)(pStokesSRWL->arS2); tdS2 = adS2;}
				if(pStokesSRWL->arS3 != 0) { adS3 = (double*)(pStokesSRWL->arS3); tdS3 = adS3;}
				for(long iz=0; iz<Nz; iz++)
				{
					for(long ix=0; ix<Nx; ix++)
					{
						for(long ie=0; ie<Ne; ie++)
						{
							if(adS0 != 0) *(tdS0++) = 0.;
							if(adS1 != 0) *(tdS1++) = 0.;
							if(adS2 != 0) *(tdS2++) = 0.;
							if(adS3 != 0) *(tdS3++) = 0.;
						}
					}
				}
			}
		}
	}

	double AzimuthOfArcAndHorLineIntersect(double y0, double r, double x0)
	{
		if(::fabs(r) < ::fabs(y0)) return -1.E+23;
		if(x0 < 0.)
		{
			if(y0 < 0.) return -PI - asin(y0/r);
			else return PI - asin(y0/r);
		}
		else return asin(y0/r);
	}
	double AzimuthOfArcAndVerLineIntersect(double x0, double r, double y0)
	{
		if(::fabs(r) < ::fabs(x0)) return -1.E+23;
		if(y0 < 0.) return acos(-x0/r) - PI;
		else return acos(x0/r);
	}

	void EstimateRadAvgOfBox(double x1, double x2, double y1, double y2, double& r)
	{
		double x1x1 = x1*x1, x2x2 = x2*x2, y1y1 = y1*y1, y2y2 = y2*y2;
		double R11 = sqrt(x1x1 + y1y1), R12 = sqrt(x1x1 + y2y2), R21 = sqrt(x2x2 + y1y1), R22 = sqrt(x2x2 + y2y2);
		double Buf = 2*(x1*(y1*R11 - y2*R12) + x2*(y2*R22 - y1*R21)) 
				   + y1y1*y1*log((x1 + R11)/(x2 + R21)) + y2y2*y2*log((x2 + R22)/(x1 + R12)) 
				   + x1x1*x1*log((y1 + R11)/(y2 + R12)) + x2x2*x2*log((y2 + R22)/(y1 + R21));
		r = Buf/(6.*(x2 - x1)*(y2 - y1));
	}

	double FreqWindowFunc(double q, double a)
	{
		if((q < -a) || (q >= a)) return 0.;
		else if(q < 0.) return (a + q)/a;
		else return (a - q)/a;
	}

	double FormalPhase(float Re, float Im)
	{// Gives -Pi <= Phase < Pi
		const double HalhPi = 1.5707963267949;
		const double Pi = 3.1415926535898;
		if(Re != 0.) 
		{
			if(Im <= 0.)
			{
				if(Re < 0.) return atan(double(Im/Re)) - Pi;
				else return atan(double(Im/Re));
			}
			else
			{
				if(Re < 0.) return atan(double(Im/Re)) + Pi;
				else return atan(double(Im/Re));
			}
		}
		else
		{
			if(Im == 0.) return  0.;
			else return (Im > 0.)? HalhPi : -HalhPi;
		}
	}

	void CosAndSin(double x, double& Cos, double& Sin)
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
	double Erf(double x)
	{
		const double SmallLargeBorder = 2.3;
		const double RelPrecForSmall = 1.E-07;
		//const double RelPrecForSmall = 1.E-09;

		const int MaxAmOfTermsForSmall = 25;
		double xe2 = x*x;
		if(::fabs(x) < SmallLargeBorder)
		{
			double an = x, Sum = x;
			for(int i=1; i<MaxAmOfTermsForSmall; i++)
			{
				int Twoi = i << 1;
				double cn = -xe2*double(Twoi - 1)/double(i*(Twoi + 1));
				an *= cn; Sum += an;
				if(::fabs(an) < RelPrecForSmall*::fabs(Sum)) break;
			}
			return Two_d_SqrtPi*Sum;
		}
		else
		{
			double b1 = -0.5/xe2;
			double b2 = b1*b1*3.;
			double b3 = b2*b1*5.;
			return ((x>=0.)? 1. : -1.) - One_d_SqrtPi*exp(-xe2)*(1. + b1 + b2 + b3)/x;
		}
	}
};

//*************************************************************************

#endif
