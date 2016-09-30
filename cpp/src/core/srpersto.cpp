/************************************************************************//**
 * File: srpersto.cpp
 * Description: Calculation of Stokes parameters / Spectral Flux of Undulator Radiation from a Finite-Emittance Electron Beam, collected by a finite-size rectangular aperture
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srpersto.h"
#include "srmagfld.h"
#include "srprgind.h"
#include "srinterf.h"

//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

int srTEnergyAzimuthGrid::SetUpCosAndSinLookUpArrays()
{
	if(!AmsOfPointsOverPhiAreConstant) return 0;

	if(CosLookUpArray != 0) delete[] CosLookUpArray; CosLookUpArray = 0;
	CosLookUpArray = new double[*AmOfAzPoints];
	if(CosLookUpArray == 0) return MEMORY_ALLOCATION_FAILURE;
	if(SinLookUpArray != 0) delete[] SinLookUpArray; SinLookUpArray = 0;
	SinLookUpArray = new double[*AmOfAzPoints];
	if(SinLookUpArray == 0) return MEMORY_ALLOCATION_FAILURE;

	const double HalfPI = 1.5707963267949; 
	const double TwoPI = 6.2831853071796;
	int HalfAmOfPhiPoints = (*AmOfAzPoints) >> 1;
	double phStep = TwoPI/(*AmOfAzPoints);
	double PhiTol = phStep*0.001; // To steer
	double Phi = 0.;
	double *tCos = CosLookUpArray, *tSin = SinLookUpArray;
	for(int iph=0; iph<(*AmOfAzPoints); iph++)
	{
		if(Phi > HalfPI + PhiTol) break;
		CosAndSinComp.CosAndSin(Phi, *tCos, *tSin); 

		int i2Ofst = HalfAmOfPhiPoints - iph;
		if(i2Ofst != iph) { *(CosLookUpArray + i2Ofst) = -(*tCos); *(SinLookUpArray + i2Ofst) = *tSin;}
		int i3Ofst = iph + HalfAmOfPhiPoints;
		if(i3Ofst != i2Ofst) { *(CosLookUpArray + i3Ofst) = -(*tCos); *(SinLookUpArray + i3Ofst) = -(*tSin);}
		int i4Ofst = (iph > 0)? ((*AmOfAzPoints) - iph) : 0;
		if((i4Ofst != i3Ofst) && (i4Ofst != iph)) { *(CosLookUpArray + i4Ofst) = *tCos; *(SinLookUpArray + i4Ofst) = -(*tSin);}
		tCos++; tSin++; Phi += phStep;
	}
	return 0;
}

//*************************************************************************

srTRadIntPeriodic::srTRadIntPeriodic(srTEbmDat* pElecBeam, srTMagFieldPeriodic* pMagFldPer, srTWfrSmp* pWfrSmp, void* pPrcPar)
{
	Initialize();

	if(pElecBeam != 0)
	{
        EbmDat = *pElecBeam;
	}
	if(pMagFldPer != 0)
	{
        MagPer = *pMagFldPer;
	}
	if(pWfrSmp != 0)
	{
        DistrInfoDat = *pWfrSmp;
		// WARNING: In this particular case the Photon Energy is treated internally in keV
        // But this is an exception (normally, everywhere it is in eV)
        DistrInfoDat.LambStart *= 0.001;
        DistrInfoDat.LambEnd *= 0.001;
	}
	if(pPrcPar != 0)
	{
		srTParPrecStokesPer* pParPrecStokesPer = (srTParPrecStokesPer*)pPrcPar;
		IntPerStoPrec.InitHarm = pParPrecStokesPer->InitHarm; 
		IntPerStoPrec.FinHarm = pParPrecStokesPer->FinHarm; 
        IntPerStoPrec.Kns = pParPrecStokesPer->PrecS;
		IntPerStoPrec.Knphi = pParPrecStokesPer->PrecPhi;
        IntPerStoPrec.IntensityOrFlux = pParPrecStokesPer->IntOrFlux;
		IntPerStoPrec.MinPhotEnExtRight = pParPrecStokesPer->MinPhotEnExtRight; //OC170713

		//if(IntPerStoPrec.IntensityOrFlux != 'f') DistrInfoDat.EnsureZeroTransverseRangesForSinglePoints();
		//this leads to bug
	}
    pWarningsGen = &gVectWarnNos;

	int res = CheckInputConsistency();
    if(res != 0) throw res;
}

//*************************************************************************

void srTRadIntPeriodic::Initialize()
{
	DistrInfoDat.CoordUnits = 0; // To ensure m for coord.
	DistrInfoDat.FluxComp = 1;

	HalfPI = 1.5707963267949;
	PI = 3.141592653590;
	TwoPI = 6.2831853071796;
	ThreePIdTwo = 4.7123889803847;
	One_dTwoPI = 0.1591549430919;
	One_dHalfPI = 0.636619772367581;
	a2c = -0.5; a4c = 0.041666666666667; a6c = -0.0013888888888889; a8c = 0.000024801587301587; a10c = -2.755731922E-07;
	a3s = -0.16666666666667; a5s = 0.0083333333333333; a7s = -0.0001984126984127; a9s = 2.755731922E-06; a11s = -2.505210839E-08;

	One_d_SqrtPi = 0.564189583547756;
	Two_d_SqrtPi = 2.*One_d_SqrtPi;

	BtxArr = BtzArr = XArr = ZArr = IntBtE2Arr = 0;

	EnergySpreadShouldBeTreated = 1;
	s0ShouldBeTreatedGen = s1ShouldBeTreatedGen = s2ShouldBeTreatedGen = s3ShouldBeTreatedGen = 1;
	PartDistrNormGen = 1.;

	SigFactGen = 4.2; // To steer. Very important parameter !!!

    pWarningsGen = &gVectWarnNos;
}

//*************************************************************************

int srTRadIntPeriodic::CheckInputConsistency()
{
	DeduceE_BeamConstants();

	if(DistrInfoDat.yStart <= 0.) return OBS_DIST_SHOULD_BE_POSITIVE;

	if((DistrInfoDat.xStart == DistrInfoDat.xEnd) && (SigmaXpGen == 0.)) return ZERO_DIVERG_AND_SLIT_SIZE;
	if((DistrInfoDat.zStart == DistrInfoDat.zEnd) && (SigmaZpGen == 0.)) return ZERO_DIVERG_AND_SLIT_SIZE;

	if((MagPer.TypeOfUnd == 3) && (MagPer.PhaseSh_OK > 4.5)) // optical klystron
	{
		//srTSend Send; Send.AddWarningMessage(pWarningsGen, TOO_LARGE_OK_PHASE_SHIFT_PARAM);
		CErrWarn::AddWarningMessage(pWarningsGen, TOO_LARGE_OK_PHASE_SHIFT_PARAM);
	}

	return 0;
}

//*************************************************************************

void srTRadIntPeriodic::ComputeStokes(srTEbmDat* pElecBeam, srTMagFieldPeriodic* pMagFldPer, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes)
{
	if((pElecBeam == 0) || (pMagFldPer == 0) || (pWfrSmp == 0) || (pPrcPar == 0) || (pStokes == 0)) throw INCORRECT_PARAMS_SR_COMP;
	srTRadIntPeriodic* pRadInt = new srTRadIntPeriodic(pElecBeam, pMagFldPer, pWfrSmp, pPrcPar);
	//int res = pRadInt->ComputeTotalStokesDistr(*pStokes);
	int res = pRadInt->ComputeTotalStokesDistr(pStokes); //OC020112
	delete pRadInt;
	if(res != 0) throw res;
}

//*************************************************************************

//int srTRadIntPeriodic::ComputeTotalStokesDistr(srTStokesStructAccessData& StokesAccessData)
int srTRadIntPeriodic::ComputeTotalStokesDistr(srTStokesStructAccessData* pStokesAccessData, SRWLStructStokes* pStokesSRWL) //OC020112
{
	int result;

	MagPer.AnalyzeFieldSymmetry();
	DeduceSlitAngDims();
	DeduceE_BeamConstants(); // This should be called after DeduceSlitAngDims()
	DeduceNormConstants();
	//ZeroOutData(StokesAccessData);
	ZeroOutData(pStokesAccessData, pStokesSRWL); //OC020112
	FindPartDistrNorm();

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	AnalizeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ);

	//long PerX = DistrInfoDat.nLamb << 2;
	//long PerZ = PerX*DistrInfoDat.nx;
	long long PerX = DistrInfoDat.nLamb << 2;
	long long PerZ = PerX*DistrInfoDat.nx;

	//long PerX1 = DistrInfoDat.nLamb;
	//long PerZ1 = PerX1*DistrInfoDat.nx;
	long long PerX1 = DistrInfoDat.nLamb;
	long long PerZ1 = PerX1*DistrInfoDat.nx;

	double xAngStart, xAngStep, zAngStart, zAngStep;
	FindAngularObsGrid(xAngStart, xAngStep, zAngStart, zAngStep);
	double xTol = xAngStep*0.001; // To steer
	double zTol = zAngStep*0.001; // To steer

	srTCompProgressIndicator CompProgressInd;
	//long TotalAmOfOutCounts = DistrInfoDat.nz*DistrInfoDat.nx*(IntPerStoPrec.FinHarm - IntPerStoPrec.InitHarm + 1);
	long long TotalAmOfOutCounts = DistrInfoDat.nz*DistrInfoDat.nx*(IntPerStoPrec.FinHarm - IntPerStoPrec.InitHarm + 1);
	if(FinalResAreSymOverX) TotalAmOfOutCounts >>= 1;
	if(FinalResAreSymOverZ) TotalAmOfOutCounts >>= 1;
	char ProgressIndicatorEnabled = (TotalAmOfOutCounts >= 5); // To steer
	double UpdateTimeInt_s = 0.5; // To steer
	if(ProgressIndicatorEnabled) if(result = CompProgressInd.InitializeIndicator(TotalAmOfOutCounts, UpdateTimeInt_s)) return result;
	//long ProgressCount = 0;
	long long ProgressCount = 0;

	for(int n=IntPerStoPrec.InitHarm; n<=IntPerStoPrec.FinHarm; n++)
	{
		DeduceNsForOnePeriod(n);
		if(result = AllocateFieldBasedArrays()) return result;
		if(result = MagPer.SetupFieldBasedArrays(EbmDat, NsGen, BtxArr, BtzArr, XArr, ZArr, IntBtE2Arr)) return result;

		double eStart = DistrInfoDat.LambStart, eFin = DistrInfoDat.LambEnd;
		long Ne = DistrInfoDat.nLamb;
		srTEnergyAzimuthGrid EnAzGrid; // Keep it as local
		if(result = DeduceGridOverPhotonEnergyAndAzimuth(n, eStart, eFin, Ne, EnAzGrid)) return result;
		if(result = EnAzGrid.SetUpCosAndSinLookUpArrays()) return result;

		//float** LongIntArrays = 0;
		double** LongIntArrays = 0; //OC020112
		int** LongIntArrInfo = 0;
		if(result = ComputeLongIntForEnAndAz(n, EnAzGrid, LongIntArrays, LongIntArrInfo)) return result;

		for(int iz=0; iz<DistrInfoDat.nz; iz++)
		{
			EXZ.z = zAngStart + iz*zAngStep;
			if(FinalResAreSymOverZ) 
			{ 
				if(EXZ.z - EbmDat.dzds0 > zTol) break;
			}

			for(int ix=0; ix<DistrInfoDat.nx; ix++)
			{
				EXZ.x = xAngStart + ix*xAngStep;
				if(FinalResAreSymOverX) 
				{ 
					if(EXZ.x - EbmDat.dxds0 > xTol) break;
				}
				//float* pStartEnSlice = StokesAccessData.pBaseSto + iz*PerZ + ix*PerX;
				float *pStartEnSlice = 0;
				if(pStokesAccessData != 0) pStartEnSlice = pStokesAccessData->pBaseSto + iz*PerZ + ix*PerX; //OC020112

				SetUpAvgEnergy(n); // Since it depends on EXZ.x, EXZ.z

				//if(result = ComputeHarmContribToSpecAtDir(n, EnAzGrid, LongIntArrays, LongIntArrInfo, pStartEnSlice)) return result;
				//long ofstSt = iz*PerZ1 + ix*PerX1; //OC020112
				long long ofstSt = iz*PerZ1 + ix*PerX1; //OC020112
				if(result = ComputeHarmContribToSpecAtDir(n, EnAzGrid, LongIntArrays, LongIntArrInfo, pStartEnSlice, pStokesSRWL, ofstSt)) return result;

				if(result = srYield.Check()) return result;
				if(ProgressIndicatorEnabled) if(result = CompProgressInd.UpdateIndicator(ProgressCount++)) return result;
			}
		}

		DisposeLongIntArraysForEnAndAz(EnAzGrid, LongIntArrays, LongIntArrInfo);
		DisposeFieldBasedArrays();
	}

	if(FinalResAreSymOverZ || FinalResAreSymOverX) 
		FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, pStokesAccessData, pStokesSRWL); //OC060812
		//FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, *pStokesAccessData); //OC020112
		//FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, StokesAccessData);

	return 0;
}

//*************************************************************************

//int srTRadIntPeriodic::ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, float** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice)
//int srTRadIntPeriodic::ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt)
int srTRadIntPeriodic::ComputeHarmContribToSpecAtDir(int n, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt)
{
	int result;

	float* HarmDataArray = new float[EnAzGrid.Ne << 2];
	//double* HarmDataArray = new double[EnAzGrid.Ne << 2]; //OC020112
	if(HarmDataArray == 0) return MEMORY_ALLOCATION_FAILURE;

	FindObservationLimits(txMin_Fphi, txMax_Fphi, tzMin_Fphi, tzMax_Fphi); // To limit computation of the Erf factor

	const double c0 = 1.239854E-09;
	double be0 = (n << 1)*c0/MagPer.PerLength;
	double be1 = EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2);
	double CritEnergyForHarm = be0/be1;

	double Con_PhPerSecPer10em3bwForE = Con_PhPerSecPer10em3bw*c0;

	const double OverCritFactor = 1.000001; // To steer

	//double AuxX = EXZ.x - EbmDat.dxds0; //OC
	//double AuxX1 = AuxX - xPartFactRangeGen; //OC
	//double AuxX2 = AuxX + xPartFactRangeGen; //OC
	//double AuxX1e2 = AuxX1*AuxX1, AuxX2e2 = AuxX2*AuxX2; //OC
	//double AuxZ = EXZ.z - EbmDat.dzds0; //OC
	//double AuxZ1 = AuxZ - zPartFactRangeGen; //OC
	//double AuxZ2 = AuxZ + zPartFactRangeGen; //OC
	//double AuxZ1e2 = AuxZ1*AuxZ1, AuxZ2e2 = AuxZ2*AuxZ2; //OC

	char NonZeroInfUndResult = 0;

	double eStep = (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1);
	EXZ.e = EnAzGrid.eStart;
	float *tHarmData = HarmDataArray;
	//double *tHarmData = HarmDataArray; //OC020112
	for(long ie=0; ie<EnAzGrid.Ne; ie++)
	{
		srTEFourier Stokes;
		if((EnAzGrid.eMinEff <= EXZ.e) && (EXZ.e <= EnAzGrid.eMaxEff*OverCritFactor))
		{
			if(EXZ.e > CritEnergyForHarm*OverCritFactor)
			{ 
				Stokes.EwX_Re = Stokes.EwX_Im = Stokes.EwZ_Re = Stokes.EwZ_Im = 0.;
			}
			else 
			{
				double ArgSqrt = be0/EXZ.e - be1;

				//OC
				//if((ArgSqrt < AuxX1e2) || (ArgSqrt < AuxX2e2) || (ArgSqrt < AuxZ1e2) || (ArgSqrt < AuxZ2e2))
				//{
				//	Stokes.EwX_Re = Stokes.EwX_Im = Stokes.EwZ_Re = Stokes.EwZ_Im = 0.;
				//}
				//else
				//{
				Tet_Fphi = (ArgSqrt > 0.)? sqrt(ArgSqrt) : 0.;

				pA_Fphi = 0;
				if(result = RestoreLongIntArray(ie, EnAzGrid, LongIntArrays, LongIntArrInfo, pA_Fphi)) return result;

				pEnAzGrid_Fphi = &EnAzGrid;
				if(result = Int1D_Simpson(0., TwoPI, (EnAzGrid.AmOfAzPoints)[ie] + 1, 'p', Stokes)) return result; 

				double Fact = Con_PhPerSecPer10em3bwForE/EXZ.e;
				if(IntPerStoPrec.IntensityOrFlux == 'i') Fact *= Inv_dsSlit_mmEm2;

				Stokes.EwX_Re *= Fact; Stokes.EwX_Im *= Fact; Stokes.EwZ_Re *= Fact; Stokes.EwZ_Im *= Fact;

				if(MagPer.FldPolar == 'v')
				{
					Stokes.EwX_Im = Stokes.EwX_Re;
					Stokes.EwZ_Re = Stokes.EwZ_Im = 0.;
				}
				else if(MagPer.FldPolar == 'h')
				{
					Stokes.EwX_Im = -Stokes.EwX_Re;
					Stokes.EwZ_Re = Stokes.EwZ_Im = 0.;
				}

				if((Stokes.EwX_Re != 0.) || (Stokes.EwX_Im != 0.) || (Stokes.EwZ_Re != 0.) || (Stokes.EwZ_Im != 0.))
					NonZeroInfUndResult = 1;

				if(pA_Fphi != 0) delete[] pA_Fphi; // delete restored temporary array
				//}
			}
		}
		else
		{ 
			Stokes.EwX_Re = Stokes.EwX_Im = Stokes.EwZ_Re = Stokes.EwZ_Im = 0.;
		}

		*(tHarmData++) = (float)Stokes.EwX_Re;
		*(tHarmData++) = (float)Stokes.EwX_Im;
		*(tHarmData++) = (float)Stokes.EwZ_Re;
		*(tHarmData++) = (float)Stokes.EwZ_Im;
		//*(tHarmData++) = Stokes.EwX_Re; //OC020112
		//*(tHarmData++) = Stokes.EwX_Im;
		//*(tHarmData++) = Stokes.EwZ_Re;
		//*(tHarmData++) = Stokes.EwZ_Im;
		EXZ.e += eStep;

		if(result = srYield.Check()) return result;
	}

	if(NonZeroInfUndResult)
	{
		if((!EnergySpreadShouldBeTreated) && (MagPer.TypeOfUnd == 0))
		{
			if(pOutEnSlice != 0)
			{
				float *tOut = pOutEnSlice, *t = HarmDataArray;
				for(long ie=0; ie<(EnAzGrid.Ne << 2); ie++) *(tOut++) += *(t++);
			}
			else if((pStokesSRWL != 0) && (ofstSt >= 0)) //OC020112
			{
				float *t = HarmDataArray;
				if(pStokesSRWL->numTypeStokes == 'f')
				{
					float *aS0 = (float*)(pStokesSRWL->arS0);
					float *aS1 = (float*)(pStokesSRWL->arS1);
					float *aS2 = (float*)(pStokesSRWL->arS2);
					float *aS3 = (float*)(pStokesSRWL->arS3);
					float *tS0 = aS0 + ofstSt, *tS1 = aS1 + ofstSt, *tS2 = aS2 + ofstSt, *tS3 = aS3 + ofstSt;
					for(long ie=0; ie<(EnAzGrid.Ne); ie++) 
					{
						if(aS0 != 0) *(tS0++) += *t;
						if(aS1 != 0) *(tS1++) += *(t + 1);
						if(aS2 != 0) *(tS2++) += *(t + 2);
						if(aS3 != 0) *(tS3++) += *(t + 3);
						t += 4;
					}
				}
				else if(pStokesSRWL->numTypeStokes == 'd')
				{
					double *aS0 = (double*)(pStokesSRWL->arS0);
					double *aS1 = (double*)(pStokesSRWL->arS1);
					double *aS2 = (double*)(pStokesSRWL->arS2);
					double *aS3 = (double*)(pStokesSRWL->arS3);
					double *tS0 = aS0 + ofstSt, *tS1 = aS1 + ofstSt, *tS2 = aS2 + ofstSt, *tS3 = aS3 + ofstSt;
					for(long ie=0; ie<(EnAzGrid.Ne); ie++) 
					{
						if(aS0 != 0) *(tS0++) += *t;
						if(aS1 != 0) *(tS1++) += *(t + 1);
						if(aS2 != 0) *(tS2++) += *(t + 2);
						if(aS3 != 0) *(tS3++) += *(t + 3);
						t += 4;
					}
				}
			}
			if(HarmDataArray != 0) delete[] HarmDataArray; HarmDataArray = 0;
			return 0;
		}

		//if(result = TreatEnergySpreadAndFiniteNumberOfPeriods(n, EnAzGrid, HarmDataArray, pOutEnSlice)) return result;
		if(result = TreatEnergySpreadAndFiniteNumberOfPeriods(n, EnAzGrid, HarmDataArray, pOutEnSlice, pStokesSRWL, ofstSt)) return result; //OC020112
	}

	if(HarmDataArray != 0) delete[] HarmDataArray;
	return 0;
}

//*************************************************************************

void srTRadIntPeriodic::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData* pStokesAccessData, SRWLStructStokes* pStokesSRWL)
//void srTRadIntPeriodic::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData& StokesAccessData)
{
	//long PerX = DistrInfoDat.nLamb << 2; //This is for IGOR version data alignment
	//long PerZ = PerX*DistrInfoDat.nx;
	long PerX=0;
	if(pStokesAccessData != 0)
	{//This is for IGOR version data alignment
		PerX = DistrInfoDat.nLamb << 2; 
	}
	if(pStokesSRWL != 0)
	{//This is for SRWLib version data alignment
		PerX = DistrInfoDat.nLamb; //This is for IGOR version data alignment
	}
	long PerZ = PerX*DistrInfoDat.nx;

	char SymWithRespectToXax, SymWithRespectToZax;
	int HalfNz = DistrInfoDat.nz >> 1, Nz_mi_1 = DistrInfoDat.nz - 1;
	
	if(FinalResAreSymOverZ && FinalResAreSymOverX)
	{
		if((HalfNz << 1) != DistrInfoDat.nz) HalfNz++;
	}

	int HalfNx = DistrInfoDat.nx >> 1, Nx_mi_1 = DistrInfoDat.nx - 1;
	int iz, ix;
	if(FinalResAreSymOverZ)
	{
		if(FinalResAreSymOverX)
		{
			SymWithRespectToXax = 0; SymWithRespectToZax = 1;
			for(iz=0; iz<HalfNz; iz++)
			{
				long izPerZ = iz*PerZ;
				for(ix=0; ix<HalfNx; ix++)
				{
					//float* pOrigData = StokesAccessData.pBaseSto + izPerZ + ix*PerX;
					//float* pSymData = StokesAccessData.pBaseSto + izPerZ + (Nx_mi_1 - ix)*PerX;
					//CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
					long ofstOrigData = izPerZ + ix*PerX;
					long ofstSymData = izPerZ + (Nx_mi_1 - ix)*PerX;
					if(pStokesAccessData != 0) //OC080812
					{
						float* pOrigData = pStokesAccessData->pBaseSto + ofstOrigData;
						float* pSymData = pStokesAccessData->pBaseSto + ofstSymData;
						CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
					}
					if(pStokesSRWL != 0)
					{
						CopySymEnergySliceSRWL(*pStokesSRWL, ofstOrigData, ofstSymData, SymWithRespectToXax, SymWithRespectToZax);
					}
				}
			}
		}
		SymWithRespectToXax = 1; SymWithRespectToZax = 0;
		for(iz=0; iz<HalfNz; iz++)
		{
			long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
			for(ix=0; ix<DistrInfoDat.nx; ix++)
			{
				long ixPerX = ix*PerX;
				//float* pOrigData = StokesAccessData.pBaseSto + izPerZ + ixPerX;
				//float* pSymData = StokesAccessData.pBaseSto + BufZ + ixPerX;
				//CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				long ofstOrigData = izPerZ + ixPerX;
				long ofstSymData = BufZ + ixPerX;
				if(pStokesAccessData != 0) //OC080812
				{
					float* pOrigData = pStokesAccessData->pBaseSto + ofstOrigData;
					float* pSymData = pStokesAccessData->pBaseSto + ofstSymData;
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
				if(pStokesSRWL != 0)
				{
					CopySymEnergySliceSRWL(*pStokesSRWL, ofstOrigData, ofstSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
	}
	else if(FinalResAreSymOverX)
	{
		SymWithRespectToXax = 0; SymWithRespectToZax = 1;
		for(iz=0; iz<DistrInfoDat.nz; iz++)
		{
			long izPerZ = iz*PerZ;
			for(ix=0; ix<HalfNx; ix++)
			{
				//float* pOrigData = StokesAccessData.pBaseSto + izPerZ + ix*PerX;
				//float* pSymData = StokesAccessData.pBaseSto + izPerZ + (Nx_mi_1 - ix)*PerX;
				//CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				long ofstOrigData = izPerZ + ix*PerX;
				long ofstSymData = izPerZ + (Nx_mi_1 - ix)*PerX;
				if(pStokesAccessData != 0) //OC080812
				{
					float* pOrigData = pStokesAccessData->pBaseSto + ofstOrigData;
					float* pSymData = pStokesAccessData->pBaseSto + ofstSymData;
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
				if(pStokesSRWL != 0)
				{
					CopySymEnergySliceSRWL(*pStokesSRWL, ofstOrigData, ofstSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
	}
}

//*************************************************************************

//int srTRadIntPeriodic::ComputeLongIntForEnAndAz(int n, srTEnergyAzimuthGrid& EnAzGrid, float**& LongIntArrays, int**& LongIntArrInfo)
int srTRadIntPeriodic::ComputeLongIntForEnAndAz(int n, srTEnergyAzimuthGrid& EnAzGrid, double**& LongIntArrays, int**& LongIntArrInfo)
{
	int result;

	//LongIntArrays = new float*[EnAzGrid.Ne];
	LongIntArrays = new double*[EnAzGrid.Ne]; //OC020112
	if(LongIntArrays == 0) return MEMORY_ALLOCATION_FAILURE;
	LongIntArrInfo = new int*[EnAzGrid.Ne];
	if(LongIntArrInfo == 0) return MEMORY_ALLOCATION_FAILURE;

	double txMin1, txMax1, tzMin1, tzMax1;
	FindObservationLimits(txMin1, txMax1, tzMin1, tzMax1);

	double txMin2 = -txMax1, txMax2 = -txMin1, tzMin2 = tzMin1, tzMax2 = tzMax1;
	double txMin3 = txMin1, txMax3 = txMax1, tzMin3 = -tzMax1, tzMax3 = -tzMin1;
	double txMin4 = -txMax1, txMax4 = -txMin1, tzMin4 = -tzMax1, tzMax4 = -tzMin1;

	const double c0 = 1.239854E-09;
	double b1 = (n << 1)*c0/MagPer.PerLength;
	double b2 = EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2);

	const double OverCritFactor = 1.000001; // To steer

	double eStep = (EnAzGrid.Ne > 1)? (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1) : 0.;
	double e = EnAzGrid.eStart;
	for(long ie=0; ie<EnAzGrid.Ne; ie++)
	{
		int AmOfPhiPoints = (EnAzGrid.AmOfAzPoints)[ie];

		//float *BufLongIntArr = 0;
		double *BufLongIntArr = 0; //OC020112
		int *BufLongIntArrInfo = 0;
		int NonZeroPhiPointsCount = 0;

		if((EnAzGrid.eMinEff <= e) && (e <= EnAzGrid.eMaxEff*OverCritFactor))
		{
			if(result = srYield.Check()) return result;

			//BufLongIntArr = new float[AmOfPhiPoints << 2];
			BufLongIntArr = new double[AmOfPhiPoints << 2]; //OC020112
			if(BufLongIntArr == 0) return MEMORY_ALLOCATION_FAILURE;
			BufLongIntArrInfo = new int[AmOfPhiPoints];
			if(BufLongIntArrInfo == 0) return MEMORY_ALLOCATION_FAILURE;

			double AgrSqrt = b1/e - b2;
			if(AgrSqrt < 0.) AgrSqrt = 0.;
			double Tet = (AgrSqrt > 0.)? sqrt(AgrSqrt) : 0.;

			double PhiStep = TwoPI/AmOfPhiPoints;
			double PhiTol = PhiStep*0.001; // To steer
			double Phi = 0.;
			//float *t0 = BufLongIntArr;
			//float *t = t0;
			double *t0 = BufLongIntArr; //OC020112
			double *t = t0;
			int *tInfo = BufLongIntArrInfo;

			for(int iph=0; iph<AmOfPhiPoints; iph++)
			{
				if((Phi > HalfPI + PhiTol) && (MagPer.FieldSymmetryInd == 0)) break; // planar

				double SinPh, CosPh; 
				if(EnAzGrid.AmsOfPointsOverPhiAreConstant)
				{
					CosPh = *(EnAzGrid.CosLookUpArray + iph); SinPh = *(EnAzGrid.SinLookUpArray + iph);
				}
				else CosAndSin(Phi, CosPh, SinPh);
				double tx = Tet*CosPh, tz = Tet*SinPh;

				char InsideSomething = (txMin1 < tx) && (tx < txMax1) && (tzMin1 < tz) && (tz < tzMax1);

				if((!InsideSomething) && (MagPer.FieldSymmetryInd == 0)) 
				{
					InsideSomething = (txMin2 < tx) && (tx < txMax2) && (tzMin2 < tz) && (tz < tzMax2);
					if(!InsideSomething)
					{
						InsideSomething = (txMin3 < tx) && (tx < txMax3) && (tzMin3 < tz) && (tz < tzMax3);
						if(!InsideSomething) InsideSomething = (txMin4 < tx) && (tx < txMax4) && (tzMin4 < tz) && (tz < tzMax4);
					}
				}

				if(InsideSomething)
				{
					srTEFourier Stokes;
					if(result = A(n, tx, tz, Stokes)) return result;
					//*t = (float)Stokes.EwX_Re; *(t + 1) = (float)Stokes.EwX_Im; *(t + 2) = (float)Stokes.EwZ_Re; *(t + 3) = (float)Stokes.EwZ_Im;
					*t = Stokes.EwX_Re; *(t + 1) = Stokes.EwX_Im; *(t + 2) = Stokes.EwZ_Re; *(t + 3) = Stokes.EwZ_Im;

					*(tInfo++) = iph; NonZeroPhiPointsCount++;

					//float s0 = *t, s1 = *(t + 1), s2 = *(t + 2), s3 = *(t + 3);
					double s0 = *t, s1 = *(t + 1), s2 = *(t + 2), s3 = *(t + 3); //OC020112
					t += 4;
					if(MagPer.FieldSymmetryInd == 0) // planar
					{
						int HalfAmOfPhiPoints = AmOfPhiPoints >> 1;

						int t2Ofst = HalfAmOfPhiPoints - iph;
						if(t2Ofst != iph)
						{
							*(tInfo++) = t2Ofst; NonZeroPhiPointsCount++;
							*(t++) = s0; *(t++) = s1; *(t++) = -s2; *(t++) = s3;
						}

						int t3Ofst = iph + HalfAmOfPhiPoints;
						if(t3Ofst != t2Ofst)
						{
							*(tInfo++) = t3Ofst; NonZeroPhiPointsCount++;
							*(t++) = s0; *(t++) = s1; *(t++) = s2; *(t++) = -s3; 
						}

						int t4Ofst = (iph > 0)? (AmOfPhiPoints - iph) : 0;
						if((t4Ofst != t3Ofst) && (t4Ofst != iph))
						{
							*(tInfo++) = t4Ofst; NonZeroPhiPointsCount++;
							*(t++) = s0; *(t++) = s1; *(t++) = -s2; *(t++) = -s3; 
						}
					}
				}
				Phi += PhiStep;
			}
		}

		if(NonZeroPhiPointsCount == 0)
		{
			LongIntArrays[ie] = 0;
			LongIntArrInfo[ie] = 0;
		}
		else
		{
			//LongIntArrays[ie] = new float[NonZeroPhiPointsCount << 2];
			LongIntArrays[ie] = new double[NonZeroPhiPointsCount << 2]; //OC020112
			if(LongIntArrays[ie] == 0) return MEMORY_ALLOCATION_FAILURE;

			LongIntArrInfo[ie] = new int[NonZeroPhiPointsCount + 1];
			if(LongIntArrInfo[ie] == 0) return MEMORY_ALLOCATION_FAILURE;

			//float *tBuf = BufLongIntArr, *tNew = LongIntArrays[ie];
			double *tBuf = BufLongIntArr, *tNew = LongIntArrays[ie]; //OC020112
			int *tBufInfo = BufLongIntArrInfo, *tNewInfo = LongIntArrInfo[ie];

			*(tNewInfo++) = NonZeroPhiPointsCount;

			for(int i=0; i<NonZeroPhiPointsCount; i++)
			{
				for(int k=0; k<4; k++) *(tNew++) = *(tBuf++);
				*(tNewInfo++) = *(tBufInfo++);
			}
		}

		e += eStep;

		if(BufLongIntArr != 0) delete[] BufLongIntArr;
		if(BufLongIntArrInfo != 0) delete[] BufLongIntArrInfo;
	}
	return 0;
}

//*************************************************************************

//int srTRadIntPeriodic::RestoreLongIntArray(long ie, srTEnergyAzimuthGrid& EnAzGrid, float** LongIntArrays, int** LongIntArrInfo, float*& pArr)
int srTRadIntPeriodic::RestoreLongIntArray(long ie, srTEnergyAzimuthGrid& EnAzGrid, double** LongIntArrays, int** LongIntArrInfo, double*& pArr)
{
	int NphiLoc = *(EnAzGrid.AmOfAzPoints + ie);
	if(NphiLoc == 0) { pArr = 0; return 0;}

	int FourNphiLoc = NphiLoc << 2;
	//pArr = new float[FourNphiLoc];
	pArr = new double[FourNphiLoc]; //OC020112
	if(pArr == 0) return MEMORY_ALLOCATION_FAILURE;

	//float *tArr = pArr;
	double *tArr = pArr; //OC020112
	for(int i=0; i<FourNphiLoc; i++) *(tArr++) = 0.;

	int *tEncInfo = LongIntArrInfo[ie];
	if(tEncInfo == 0) return 0;

	int AmOfNonZeroPhiPoints = *(tEncInfo++);

	//float *tEnc = LongIntArrays[ie];
	double *tEnc = LongIntArrays[ie]; //OC020112
	for(int k=0; k<AmOfNonZeroPhiPoints; k++)
	{
		//float *pArrEl = pArr + ((*(tEncInfo++)) << 2);
		double *pArrEl = pArr + ((*(tEncInfo++)) << 2);
		*(pArrEl++) = *(tEnc++); *(pArrEl++) = *(tEnc++); *(pArrEl++) = *(tEnc++); *pArrEl = *(tEnc++); 
	}
	return 0;
}

//*************************************************************************

void srTRadIntPeriodic::FindAngularObsGrid(double& xAngStart, double& xAngStep, double& zAngStart, double& zAngStep)
{
	if(DistrInfoDat.nx > 1)
	{
		xAngStart = DistrInfoDat.xStart/DistrInfoDat.yStart;
		double xAngEnd = DistrInfoDat.xEnd/DistrInfoDat.yStart;
		xAngStep = (xAngEnd - xAngStart)/(DistrInfoDat.nx - 1);
	}
	else
	{
		xAngStart = 0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd)/DistrInfoDat.yStart;
		xAngStep = 0.;
	}
	if(DistrInfoDat.nz > 1)
	{
		zAngStart = DistrInfoDat.zStart/DistrInfoDat.yStart;
		double zAngEnd = DistrInfoDat.zEnd/DistrInfoDat.yStart;
		zAngStep = (zAngEnd - zAngStart)/(DistrInfoDat.nz - 1);
	}
	else
	{
		zAngStart = 0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd)/DistrInfoDat.yStart;
		zAngStep = 0.;
	}
}

//*************************************************************************

void srTRadIntPeriodic::FindObservationLimits(double& txMin, double& txMax, double& tzMin, double& tzMax)
{
	double xAngStart, xAngEnd, zAngStart, zAngEnd;
	if(DistrInfoDat.nx > 1)
	{
		xAngStart = DistrInfoDat.xStart/DistrInfoDat.yStart;
		xAngEnd = DistrInfoDat.xEnd/DistrInfoDat.yStart;
	}
	else
	{
		xAngStart = 0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd)/DistrInfoDat.yStart;
		xAngEnd = xAngStart;
	}
	if(DistrInfoDat.nz > 1)
	{
		zAngStart = DistrInfoDat.zStart/DistrInfoDat.yStart;
		zAngEnd = DistrInfoDat.zEnd/DistrInfoDat.yStart;
	}
	else
	{
		zAngStart = 0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd)/DistrInfoDat.yStart;
		zAngEnd = zAngStart;
	}

	txMin = xAngStart - EbmDat.dxds0 - xPartFactRangeGen;
	txMax = xAngEnd - EbmDat.dxds0 + xPartFactRangeGen;
	tzMin = zAngStart - EbmDat.dzds0 - zPartFactRangeGen;
	tzMax = zAngEnd - EbmDat.dzds0 + zPartFactRangeGen;
}

//*************************************************************************

void srTRadIntPeriodic::FindLeastAndMostOffsetPixelCenters(double& tx0, double& tz0, double& tx1, double& tz1)
{// tx0, tz0 is the closest center

	double xAngStart, xAngEnd, zAngStart, zAngEnd;
	if(DistrInfoDat.nx > 1)
	{
		xAngStart = DistrInfoDat.xStart/DistrInfoDat.yStart;
		xAngEnd = DistrInfoDat.xEnd/DistrInfoDat.yStart;
	}
	else
	{
		xAngStart = 0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd)/DistrInfoDat.yStart;
		xAngEnd = xAngStart;
	}
	if(DistrInfoDat.nz > 1)
	{
		zAngStart = DistrInfoDat.zStart/DistrInfoDat.yStart;
		zAngEnd = DistrInfoDat.zEnd/DistrInfoDat.yStart;
	}
	else
	{
		zAngStart = 0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd)/DistrInfoDat.yStart;
		zAngEnd = zAngStart;
	}

	double txMin = xAngStart - EbmDat.dxds0, txMax = xAngEnd - EbmDat.dxds0;
	double tzMin = zAngStart - EbmDat.dzds0, tzMax = zAngEnd - EbmDat.dzds0;

	double txMinE2 = txMin*txMin, txMaxE2 = txMax*txMax, tzMinE2 = tzMin*tzMin, tzMaxE2 = tzMax*tzMax;
	double Tet00e2 = txMinE2 + tzMinE2,  Tet01e2 = txMinE2 + tzMaxE2, Tet10e2 = txMaxE2 + tzMinE2, Tet11e2 = txMaxE2 + tzMaxE2;
	
	double TetMinE2 = Tet00e2;
	tx0 = txMin; tz0 = tzMin;
	if(Tet01e2 < TetMinE2) { TetMinE2 = Tet01e2; tx0 = txMin; tz0 = tzMax;}
	if(Tet10e2 < TetMinE2) { TetMinE2 = Tet10e2; tx0 = txMax; tz0 = tzMin;}
	if(Tet11e2 < TetMinE2) { TetMinE2 = Tet11e2; tx0 = txMax; tz0 = tzMax;} 
	
	double TetMaxE2 = Tet00e2;
	tx1 = txMin; tz1 = tzMin; 
	if(Tet01e2 > TetMaxE2) { TetMaxE2 = Tet01e2; tx1 = txMin; tz1 = tzMax;}
	if(Tet10e2 > TetMaxE2) { TetMaxE2 = Tet10e2; tx1 = txMax; tz0 = tzMin;}
	if(Tet11e2 > TetMaxE2) { TetMaxE2 = Tet11e2; tx1 = txMax; tz1 = tzMax;} 

	if((txMin < 0) && (txMax > 0)) // to improve
	{
		tx0 = 0.; 
	}
	if((tzMin < 0) && (tzMax > 0))
	{
		tz0 = 0.; 
	}
}

//*************************************************************************

void srTRadIntPeriodic::FindTetMinMaxE2_FromTetxTetz(double txMin, double txMax, double tzMin, double tzMax, double& TetMinE2, double& TetMaxE2)
{
	double txMinE2 = txMin*txMin, txMaxE2 = txMax*txMax, tzMinE2 = tzMin*tzMin, tzMaxE2 = tzMax*tzMax;
	double Tet00e2 = txMinE2 + tzMinE2, Tet01e2 = txMinE2 + tzMaxE2, Tet10e2 = txMaxE2 + tzMinE2, Tet11e2 = txMaxE2 + tzMaxE2;
	
	TetMaxE2 = Tet00e2;
	if(Tet01e2 > TetMaxE2) TetMaxE2 = Tet01e2;
	if(Tet10e2 > TetMaxE2) TetMaxE2 = Tet10e2;
	if(Tet11e2 > TetMaxE2) TetMaxE2 = Tet11e2;

	char PixelCoversZeroPoint = ((txMin < 0.) && (txMax > 0.) && (tzMin < 0.) && (tzMax > 0.));

	TetMinE2 = Tet00e2;
	if(PixelCoversZeroPoint) TetMinE2 = 0.;
	else
	{
		if(Tet01e2 < TetMinE2) TetMinE2 = Tet01e2;
		if(Tet10e2 < TetMinE2) TetMinE2 = Tet10e2;
		if(Tet11e2 < TetMinE2) TetMinE2 = Tet11e2;
	}
}

//*************************************************************************

int srTRadIntPeriodic::AllocateLongIntArraysForEnAndAz(srTEnergyAzimuthGrid& EnAzGrid, float**& LongIntArrays)
{
	LongIntArrays = new float*[EnAzGrid.Ne];
	if(LongIntArrays == 0) return MEMORY_ALLOCATION_FAILURE;

	for(long i=0; i<EnAzGrid.Ne; i++)
	{
		LongIntArrays[i] = new float[(*(EnAzGrid.AmOfAzPoints + i)) << 2];
		if(LongIntArrays[i] == 0)
		{
			for(long k=0; k<i; k++) if(LongIntArrays[k] != 0) delete[] (LongIntArrays[k]);
			delete[] LongIntArrays; LongIntArrays = 0;
			return MEMORY_ALLOCATION_FAILURE;
		}
	}
	return 0;
}

//*************************************************************************

//void srTRadIntPeriodic::DisposeLongIntArraysForEnAndAz(srTEnergyAzimuthGrid& EnAzGrid, float**& LongIntArrays, int**& LongIntArrInfo)
void srTRadIntPeriodic::DisposeLongIntArraysForEnAndAz(srTEnergyAzimuthGrid& EnAzGrid, double**& LongIntArrays, int**& LongIntArrInfo)
{
	if(LongIntArrays == 0) return;

	for(long i=0; i<EnAzGrid.Ne; i++)
	{
		if(LongIntArrays[i] != 0) delete[] (LongIntArrays[i]);
		if(LongIntArrInfo[i] != 0) delete[] (LongIntArrInfo[i]);
	}
	delete[] LongIntArrays; LongIntArrays = 0;
	delete[] LongIntArrInfo; LongIntArrInfo = 0;
}

//*************************************************************************

void srTRadIntPeriodic::FindIntegralOfInfNperData(int n, srTEnergyAzimuthGrid& EnAzGrid, float* InfNperHarmData)
{
	int Np = EnAzGrid.Ne;
	double eStart = EnAzGrid.eStart;
	double eFin = EnAzGrid.eFin;
	double eStep = (eFin - eStart)/(Np - 1);

	float F[4], Edges[]={0.,0.,0.,0.}, Sum1[]={0.,0.,0.,0.}, Sum2[]={0.,0.,0.,0.}, ExtraInt[]={0.,0.,0.,0.};
	float *pStartDat = InfNperHarmData;

	int s;
	for(s=0; s<4; s++) F[s] = *(pStartDat + s);

	float *tDat = pStartDat + 4;

	if(((Np >> 1) << 1) == Np)
	{
		float HalfStep = (float)(0.5*eStep);

		for(s=0; s<4; s++) 
		{
			float Buf = *(tDat + s);
			ExtraInt[s] = HalfStep*(F[s] + Buf);
			F[s] = Buf;
		}
		Np--; tDat += 4;
	}
	
	float *pLast = pStartDat + ((EnAzGrid.Ne - 1) << 2);
	for(s=0; s<4; s++) Edges[s] = F[s] + *(pLast + s);

	for(long i=0; i<((Np-3)>>1); i++)
	{
		float *tDatp4 = tDat + 4;
		for(s=0; s<4; s++) 
		{
			Sum1[s] += *(tDat++); Sum2[s] += *(tDatp4++);
		}
		tDat = tDatp4;
	}
	for(s=0; s<4; s++) Sum1[s] += *(tDat++);

	double Mult = 0.3333333333*eStep;
	for(s=0; s<4; s++) F[s] = (float)(Mult*(Edges[s] + 4.*Sum1[s] + 2.*Sum2[s]) + ExtraInt[s]);
	srTEFourier& Res = EnAzGrid.IntOfInfNperData;
	Res.EwX_Re = F[0]; Res.EwX_Im = F[1]; Res.EwZ_Re = F[2]; Res.EwZ_Im = F[3];
}

//*************************************************************************

//int srTRadIntPeriodic::FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* InfNperHarmData, float* pOutEnSlice)
//int srTRadIntPeriodic::FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* InfNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt) //OC020112
int srTRadIntPeriodic::FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* InfNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt) //OC020112
{
	int result;

	double eStart = DistrInfoDat.LambStart, eFin = DistrInfoDat.LambEnd;
	long Ne = DistrInfoDat.nLamb;
	srTEnergyAzimuthGrid EnAzGridLoc; // Keep it as local
	EnAzGridLoc.EnsureEnResolvingObsPixels = 0;
	if(result = DeduceGridOverPhotonEnergyAndAzimuth(n, eStart, eFin, Ne, EnAzGridLoc)) return result;

	FindIntegralOfInfNperData(n, EnAzGrid, InfNperHarmData);
	EnAzGridLoc.IntOfInfNperData = EnAzGrid.IntOfInfNperData;

	float DummyF;
	//return TreatEnergySpreadAndFiniteNumberOfPeriods(n, EnAzGridLoc, &DummyF, pOutEnSlice);
	return TreatEnergySpreadAndFiniteNumberOfPeriods(n, EnAzGridLoc, &DummyF, pOutEnSlice, pStokesSRWL, ofstSt);
}

//*************************************************************************

//int srTRadIntPeriodic::TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice)
//int srTRadIntPeriodic::TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt) //OC020112
int srTRadIntPeriodic::TreatEnergySpreadAndFiniteNumberOfPeriods(int n, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt) //OC020112
{
	if(FilamentTreatmentIsPossible(EnAzGrid)) return FilamentTreatEnergySpreadAndFiniteNumberOfPeriods(n, EnAzGrid, FinNperHarmData, pOutEnSlice, pStokesSRWL, ofstSt);

	double eStep = (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1);
	double eStart = EnAzGrid.eStart - eStep*EnAzGrid.NeExtraLeft;
	long Ne = EnAzGrid.Ne + EnAzGrid.NeExtraLeft + EnAzGrid.NeExtraRight;
	double eFin = eStart + eStep*(Ne - 1);

	float *ConvFactorData = new float[Ne << 1];
	if(ConvFactorData == 0) return MEMORY_ALLOCATION_FAILURE;

	int result;
	if(result = SetupConvolutionData(n, ConvFactorData, eStart, eFin, Ne))
	{
		if(ConvFactorData != 0) delete[] ConvFactorData; return result;
	}

	if(result = srYield.Check()) return result;

	double EnergyCut = 1.239854E-09*(n << 1)/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	if((eStart < EnergyCut) && (EnergyCut <= eFin))
	{
		double eStartFict = -eStep*(Ne >> 1);
  		CGenMathFFT1DInfo FFT1DInfo;
		FFT1DInfo.xStep = eStep;
		FFT1DInfo.xStart = eStartFict;
		FFT1DInfo.Nx = Ne;
		FFT1DInfo.UseGivenStartTrValue = 0;
		FFT1DInfo.Dir = 1;
		FFT1DInfo.TreatSharpEdges = 1;
		FFT1DInfo.LeftSharpEdge = eStartFict;
		FFT1DInfo.RightSharpEdge = (EnergyCut - eStart) + eStartFict;

		CGenMathFFT1D FFT1D;
		if(result = FFT1D.SetupAuxDataForSharpEdgeCorr(FFT1DInfo, AuxDataForSharpEdgeCorrGen)) return result;
	}

	for(int iSto=0; iSto<4; iSto++)
	{
		if((pOutEnSlice == 0) && (pStokesSRWL != 0) && (ofstSt > 0))
		{//OC020112
			if(((iSto == 0) && (pStokesSRWL->arS0 == 0)) ||
			   ((iSto == 1) && (pStokesSRWL->arS1 == 0)) ||
			   ((iSto == 2) && (pStokesSRWL->arS2 == 0)) ||
			   ((iSto == 3) && (pStokesSRWL->arS3 == 0))) continue;
		}

		//if(result = ConvStokesCompon(iSto, EnAzGrid, FinNperHarmData, ConvFactorData, pOutEnSlice)) return result;
		if(result = ConvStokesCompon(iSto, EnAzGrid, FinNperHarmData, ConvFactorData, pOutEnSlice, pStokesSRWL, ofstSt)) return result; //OC020112

		if(result = srYield.Check()) return result;
	}

	AuxDataForSharpEdgeCorrGen.Dispose();
	if(ConvFactorData != 0) delete[] ConvFactorData;

	return 0;
}

//*************************************************************************

//int srTRadIntPeriodic::ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice)
//int srTRadIntPeriodic::ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long ofstSt) //OC020112
int srTRadIntPeriodic::ConvStokesCompon(int StokesNo, srTEnergyAzimuthGrid& EnAzGrid, float* FinNperHarmData, float* ConvFactorData, float* pOutEnSlice, SRWLStructStokes* pStokesSRWL, long long ofstSt) //OC020112
{
	int result;

	long Ne = EnAzGrid.Ne + EnAzGrid.NeExtraLeft + EnAzGrid.NeExtraRight;
	double eStep = (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1);
	double eStart = EnAzGrid.eStart - eStep*EnAzGrid.NeExtraLeft;
	//double eFin = eStart + eStep*(Ne - 1);

	long PointsToAddOnTheOtherSide = EnAzGrid.EnergyAdjustFinGridPar.AmOfExtraStepsFromLeftAndRightOnTheOtherSide << 1;
	long NewNe = Ne + PointsToAddOnTheOtherSide;

  	float *AuxCont = new float[NewNe << 2];
	if(AuxCont == 0) return MEMORY_ALLOCATION_FAILURE;

	long TwoNe = Ne << 1, i;
	double eStartFict = -eStep*(Ne >> 1);
  	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = eStep;
	FFT1DInfo.xStart = eStartFict;
	FFT1DInfo.Nx = Ne;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
	FFT1DInfo.pInData = AuxCont;
	FFT1DInfo.pOutData = AuxCont + TwoNe;
	CGenMathFFT1D FFT1D;

	float *t;
	if(EnAzGrid.EnsureEnResolvingObsPixels)
	{
		t = AuxCont;
		for(i=0; i<TwoNe; i++) *(t++) = 0.;

		float *tIn = FinNperHarmData + StokesNo; 
		t = AuxCont + (EnAzGrid.NeExtraLeft << 1);
		for(i=0; i<EnAzGrid.Ne; i++) { *t = *tIn; t += 2; tIn += 4;}

		if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

		if(AuxDataForSharpEdgeCorrGen.WasSetUp) FFT1D.MakeSharpEdgeCorr(FFT1DInfo, AuxDataForSharpEdgeCorrGen);
	}
	else
	{
		FFT1D.SetupLimitsTr(FFT1DInfo);

		float NormConst;
		if(StokesNo == 0) NormConst = (float)EnAzGrid.IntOfInfNperData.EwX_Re;
		else if(StokesNo == 1) NormConst = (float)EnAzGrid.IntOfInfNperData.EwX_Im;
		else if(StokesNo == 2) NormConst = (float)EnAzGrid.IntOfInfNperData.EwZ_Re;
		else NormConst = (float)EnAzGrid.IntOfInfNperData.EwZ_Im;

		t = FFT1DInfo.pOutData;
		for(i=0; i<Ne; i++) { *(t++) = NormConst; *(t++) = 0.;}

		if(StokesNo == 0) EnAzGrid.EnergyAdjustFinGridPar.eShift += ((eStart - eStartFict) - EnAvg_FinNper);
	}

	float *tGen = FFT1DInfo.pOutData, *tCon = ConvFactorData;
	for(long j=0; j<Ne; j++)
	{
		float &ReGen = *tGen, &ImGen = *(tGen + 1);
		float Re = *tCon, Im = *(tCon + 1);

		float ReRes = ReGen*Re - ImGen*Im;
		float ImRes = ImGen*Re + ReGen*Im;

		ReGen = ReRes; ImGen = ImRes;

		tGen += 2; tCon += 2;
	}

	FFT1DInfo.xStep = FFT1DInfo.xStepTr;
	FFT1DInfo.xStart = FFT1DInfo.xStartTr;
	FFT1DInfo.Dir = -1;
	FFT1DInfo.pInData = AuxCont + TwoNe;
	FFT1DInfo.pOutData = AuxCont;

	if(::fabs(EnAzGrid.EnergyAdjustFinGridPar.eShift) > 1.E-06*eStep) // if there is a shift
	{
		FFT1DInfo.UseGivenStartTrValue = 1;
		FFT1DInfo.xStartTr = eStartFict + EnAzGrid.EnergyAdjustFinGridPar.eShift;
	}

	if(PointsToAddOnTheOtherSide != 0)
	{
		long Offset = PointsToAddOnTheOtherSide*3;
		float *tOld = AuxCont + (TwoNe << 1) - 1, *tNew = tOld + Offset;
		for(i=0; i<TwoNe; i++) *(tNew--) = *(tOld--);

		float *pNewStartData = AuxCont + TwoNe + (PointsToAddOnTheOtherSide << 1);
		float *t1 = pNewStartData, *t2 = pNewStartData + PointsToAddOnTheOtherSide + TwoNe;
		for(i=0; i<PointsToAddOnTheOtherSide; i++) { *(t1++) = 0.; *(t2++) = 0.;}

		FFT1DInfo.Nx += PointsToAddOnTheOtherSide;
		FFT1DInfo.xStart -= FFT1DInfo.xStepTr*EnAzGrid.EnergyAdjustFinGridPar.AmOfExtraStepsFromLeftAndRightOnTheOtherSide;
		FFT1DInfo.pInData = pNewStartData;
	}

	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	int Multip = EnAzGrid.EnergyAdjustFinGridPar.Multip;
	long AmStExtraLeft = EnAzGrid.EnergyAdjustFinGridPar.AmStExtraLeft;
	long AmStExtraRight = EnAzGrid.EnergyAdjustFinGridPar.AmStExtraRight;
	long AmOfInitPoToIgnore = EnAzGrid.EnergyAdjustFinGridPar.AmOfInitPoToIgnore;

	t = FFT1DInfo.pOutData + (AmOfInitPoToIgnore << 1);
	long AmOfPoints = NewNe - AmOfInitPoToIgnore;
	int CountMult = 0;

	if(pOutEnSlice != 0)
	{
		float *tOut = pOutEnSlice + (AmStExtraLeft*4) + StokesNo;
		if(AmStExtraLeft < 0)
		{
			tOut = pOutEnSlice + StokesNo;
			t = FFT1DInfo.pOutData + 2*(-AmStExtraLeft)*Multip;
			AmOfPoints += AmStExtraLeft*Multip;
		}
		if(AmStExtraRight < 0) AmOfPoints += AmStExtraRight*Multip;

		for(long ie=0; ie<AmOfPoints; ie++)
		{
			if(CountMult == 0) { *tOut += *t; tOut += 4;}
			t += 2; CountMult++;
			if(CountMult == Multip) CountMult = 0;
		}
	}
	else if((pStokesSRWL != 0) && (ofstSt >= 0)) //OC020112
	{
		if(pStokesSRWL->numTypeStokes == 'f')
		{
			//float *tOut = pOutEnSlice + (AmStExtraLeft*4) + StokesNo;
			float *tOut = 0;
			if(StokesNo == 0) tOut = (float*)(pStokesSRWL->arS0) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 1) tOut = (float*)(pStokesSRWL->arS1) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 2) tOut = (float*)(pStokesSRWL->arS2) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 3) tOut = (float*)(pStokesSRWL->arS3) + ofstSt + AmStExtraLeft;
			
			if(AmStExtraLeft < 0)
			{
				//tOut = pOutEnSlice + StokesNo;
				if(StokesNo == 0) tOut = (float*)(pStokesSRWL->arS0) + ofstSt;
				else if(StokesNo == 1) tOut = (float*)(pStokesSRWL->arS1) + ofstSt;
				else if(StokesNo == 2) tOut = (float*)(pStokesSRWL->arS2) + ofstSt;
				else if(StokesNo == 3) tOut = (float*)(pStokesSRWL->arS3) + ofstSt;
			
				t = FFT1DInfo.pOutData + 2*(-AmStExtraLeft)*Multip;
				AmOfPoints += AmStExtraLeft*Multip;
			}

			if(AmStExtraRight < 0) AmOfPoints += AmStExtraRight*Multip;
	
			for(long ie=0; ie<AmOfPoints; ie++)
			{
				if(CountMult == 0) *(tOut++) += *t;
				t += 2; CountMult++;
				if(CountMult == Multip) CountMult = 0;
			}
		}
		else if(pStokesSRWL->numTypeStokes == 'd')
		{
			//float *tOut = pOutEnSlice + (AmStExtraLeft*4) + StokesNo;
			double *tOut = 0;
			if(StokesNo == 0) tOut = (double*)(pStokesSRWL->arS0) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 1) tOut = (double*)(pStokesSRWL->arS1) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 2) tOut = (double*)(pStokesSRWL->arS2) + ofstSt + AmStExtraLeft;
			else if(StokesNo == 3) tOut = (double*)(pStokesSRWL->arS3) + ofstSt + AmStExtraLeft;

			if(AmStExtraLeft < 0)
			{
				//tOut = pOutEnSlice + StokesNo;
				if(StokesNo == 0) tOut = (double*)(pStokesSRWL->arS0) + ofstSt;
				else if(StokesNo == 1) tOut = (double*)(pStokesSRWL->arS1) + ofstSt;
				else if(StokesNo == 2) tOut = (double*)(pStokesSRWL->arS2) + ofstSt;
				else if(StokesNo == 3) tOut = (double*)(pStokesSRWL->arS3) + ofstSt;
			
				t = FFT1DInfo.pOutData + 2*(-AmStExtraLeft)*Multip;
				AmOfPoints += AmStExtraLeft*Multip;
			}

			if(AmStExtraRight < 0) AmOfPoints += AmStExtraRight*Multip;

			for(long ie=0; ie<AmOfPoints; ie++)
			{
				if(CountMult == 0) *(tOut++) += *t;
				t += 2; CountMult++;
				if(CountMult == Multip) CountMult = 0;
			}
		}
	}

// Remove extra data if Multip != 1
// Now pad "0" on this side according to 
// EnAzGrid.EnergyAdjustFinGridPar.AmStExtraLeft
// EnAzGrid.EnergyAdjustFinGridPar.AmStExtraRight
// Copy or add to the real (!) pOutEnSlice

	if(AuxCont != 0) delete[] AuxCont;
	return 0;
}

//*************************************************************************

void srTRadIntPeriodic::CorrectGridToAllowRangeResizeOnTheOtherSide(srTEnergyAzimuthGrid& EnAzGrid)
{
	const double RelTol = 1.E-06; // To steer

	double eStepInput = (DistrInfoDat.LambEnd - DistrInfoDat.LambStart)/(DistrInfoDat.nLamb - 1);
	double eStepOld = (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1);
	long Ne = EnAzGrid.NeExtraLeft + EnAzGrid.Ne + EnAzGrid.NeExtraRight;
	long HalfNe = Ne >> 1;
	//double emOld = eStepOld*HalfNe;

	long AmOfStepsLeftFromCritEn = long((EnAzGrid.eCritForHarm - EnAzGrid.eStart)/eStepOld + 1.E-06);
	long AmOfStepsRightFromCritEn = long((EnAzGrid.eFin - EnAzGrid.eCritForHarm)/eStepOld + 1.E-06);

	double StepRatOld = eStepOld/eStepInput;
	double StepRatOld_mi_1;
	double eStepNew, StepRatNew;
	long Nq;
	int &Multip = EnAzGrid.EnergyAdjustFinGridPar.Multip;
	if(::fabs(StepRatOld - 1.) < RelTol)
	{
		EnAzGrid.EnergyAdjustFinGridPar.AmOfExtraStepsFromLeftAndRightOnTheOtherSide = 0;
		Multip = 1;

		eStepNew = eStepOld;
		Nq = Ne;
	}
	else
	{
		Multip = int(::fabs(eStepInput/eStepOld - RelTol)) + 1;

		StepRatOld *= Multip;
		StepRatOld_mi_1 = StepRatOld - 1;

		long k = long(StepRatOld_mi_1*HalfNe + 1.E-06);
		long NqOld = Ne + (k << 1);
		Nq = NqOld;
		CGenMathFFT1D FFT1D;
		FFT1D.NextCorrectNumberForFFT(Nq);
		k += ((Nq - NqOld) >> 1);

		StepRatNew = 1. + double(k)/double(HalfNe);
		eStepNew = StepRatNew*eStepInput/double(Multip);

		EnAzGrid.eStart = EnAzGrid.eCritForHarm - AmOfStepsLeftFromCritEn*eStepNew;
		EnAzGrid.eFin = EnAzGrid.eCritForHarm + AmOfStepsRightFromCritEn*eStepNew;

		EnAzGrid.EnergyAdjustFinGridPar.AmOfExtraStepsFromLeftAndRightOnTheOtherSide = k;
	}

	double &eShift = EnAzGrid.EnergyAdjustFinGridPar.eShift;
	int &AmOfInitPoToIgnore = EnAzGrid.EnergyAdjustFinGridPar.AmOfInitPoToIgnore;
	double eStartNew = EnAzGrid.eStart - eStepNew*EnAzGrid.NeExtraLeft;
	double eStepCorr = eStepInput/double(Multip);
	double Ofst = eStartNew - DistrInfoDat.LambStart;

	double eStartUse;
	if(Ofst > 0.)
	{
		long ieOfst = long(::fabs(Ofst)/eStepInput + 1.E-06);
		eStartUse = DistrInfoDat.LambStart + (ieOfst + 1)*eStepInput;
	}
	else
	{
		eStartUse = DistrInfoDat.LambStart;
	}
	double eNotUsed = eStartUse - eStartNew;
	AmOfInitPoToIgnore = int(eNotUsed/eStepCorr + 1.E-06);
	eShift = eNotUsed - AmOfInitPoToIgnore*eStepCorr;

	long &AmStExtraLeft = EnAzGrid.EnergyAdjustFinGridPar.AmStExtraLeft;
	long &AmStExtraRight = EnAzGrid.EnergyAdjustFinGridPar.AmStExtraRight;

	double eStartAfterShift = eStartNew + AmOfInitPoToIgnore*eStepCorr + eShift; // Modify here
	AmStExtraLeft = long(::fabs(eStartAfterShift - DistrInfoDat.LambStart)/eStepInput + 1.E-06);

	long NeAfterBackFFT = long(double(Nq - AmOfInitPoToIgnore - 1)/double(Multip) + 1.E-06) + 1;
	double eFinAfterShift = eStartAfterShift + eStepInput*(NeAfterBackFFT - 1);
	char SigRight = (eFinAfterShift <= DistrInfoDat.LambEnd)? 1 : -1;
	AmStExtraRight = SigRight*long(::fabs(DistrInfoDat.LambEnd - eFinAfterShift)/eStepInput + 1.E-06);
}

//*************************************************************************

int srTRadIntPeriodic::SetupConvolutionData_Normal(int n, float* ConvFactorData, double eStart, double eFin, long Ne)
{
	const double MinArgExp = -20.; // To steer

	double eStep = (eFin - eStart)/(Ne - 1);
	double eStartFict = -eStep*(Ne >> 1);

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = eStep;
	FFT1DInfo.xStart = eStartFict;
	FFT1DInfo.Nx = Ne;
	FFT1DInfo.UseGivenStartTrValue = 0;
	CGenMathFFT1D FFT1D;
	FFT1D.SetupLimitsTr(FFT1DInfo);
	double q = FFT1DInfo.xStartTr;
	double qStep = FFT1DInfo.xStepTr;

	double a = n*NperGen/EnAvg_FinNper;

	const double c0 = 1.239854E-09;
	const double TwoPIe2 = TwoPI*PI;
	double SigmaE_Gamma = 4*n*c0*EbmDat.SigmaRelE/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double MultArgEnSpr = -TwoPIe2*SigmaE_Gamma*SigmaE_Gamma;

	float *t = ConvFactorData;
	for(long i=0; i<Ne; i++)
	{
		float Fact = float(FreqWindowFunc(q, a));
		if((Fact > 0.) && EnergySpreadShouldBeTreated)
		{
			double ArgExp = MultArgEnSpr*q*q;
			if(ArgExp > MinArgExp) Fact *= (float)exp(ArgExp);
			else Fact = 0.;
		}
		*(t++) = Fact; *(t++) = 0.;
		q += qStep;
	}

	return 0;
}

//*************************************************************************

int srTRadIntPeriodic::SetupConvolutionData_OpticalKlystron(int n, float* ConvFactorData, double eStart, double eFin, long Ne)
{
	const double MinArgExp = -20.; // To steer

	double eStep = (eFin - eStart)/(Ne - 1);
	double eStartFict = -eStep*(Ne >> 1);

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = eStep;
	FFT1DInfo.xStart = eStartFict;
	FFT1DInfo.Nx = Ne;
	FFT1DInfo.UseGivenStartTrValue = 0;
	CGenMathFFT1D FFT1D;
	FFT1D.SetupLimitsTr(FFT1DInfo);
	double q = FFT1DInfo.xStartTr;
	double qStep = FFT1DInfo.xStepTr;

	const double c0 = 1.239854E-09;
	double Buf0 = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double CritEnergyForHarm = (n << 1)*Buf0;

	double AlphaOK = MagPer.PhaseSh_OK*EnAvg_FinNper/CritEnergyForHarm;

	double a0 = 0.5*n*NperGen/EnAvg_FinNper;
	double ac = a0*(1. + AlphaOK);
	double a1 = ac - a0, a2 = ac + a0;

	double PiNnAl = PI*NperGen*n*AlphaOK;
	double CosPiNnAl, SinPiNnAl; CosAndSin(PiNnAl, CosPiNnAl, SinPiNnAl);

	const double TwoPIe2 = TwoPI*PI;
	double SigmaE_Gamma = 4*n*c0*EbmDat.SigmaRelE/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double MultArgEnSpr = -TwoPIe2*SigmaE_Gamma*SigmaE_Gamma;

	float *t = ConvFactorData;
	for(long i=0; i<Ne; i++)
	{
		char InInt1=0, InInt2=0, InInt3=0, InInt4=0, InInt5=0, InInt6=0;

		InInt1 = ((-a2 < q) && (q <= -ac))? 1 : 0;
		if(!InInt1) InInt2 = ((-ac < q) && (q <= -a1))? 1 : 0;

		InInt3 = ((-a0 < q) && (q <= 0))? 1 : 0;
		if(!InInt3) InInt4 = ((0. < q) && (q <= a0))? 1 : 0;

		InInt5 = ((a1 < q) && (q <= ac))? 1 : 0;
		if(!InInt5) InInt6 = ((ac < q) && (q <= a2))? 1 : 0;

		if(InInt1 || InInt2 || InInt3 || InInt4 || InInt5 || InInt6)
		{
			float Fact = 0., FactRe = 0., FactIm = 0.;
			if(InInt1) 
			{
				FactRe = (float)(0.5*(a2 + q)/(a2 - ac));
				FactIm = -FactRe;
			}
			else if(InInt2) 
			{
				FactRe = (float)(0.5*(a1 + q)/(a1 - ac));
				FactIm = -FactRe;
			}
			else if(InInt5) 
			{
				FactRe = (float)(0.5*(a1 - q)/(a1 - ac));
				FactIm = FactRe;
			}
			else if(InInt6) 
			{
				FactRe = (float)(0.5*(a2 - q)/(a2 - ac));
				FactIm = FactRe;
			}

			if(InInt3) Fact += (float)(1. + q/a0);
			else if(InInt4) Fact += (float)(1. - q/a0);

			Fact += (float)(CosPiNnAl*FactRe);
			FactIm *= (float)SinPiNnAl;

			if(EnergySpreadShouldBeTreated)
			{
				double ArgExp = MultArgEnSpr*q*q;

				if(ArgExp > MinArgExp) 
				{
					double EnSprF = exp(ArgExp);
					Fact *= (float)EnSprF;
					if(FactIm != 0.) FactIm *= (float)EnSprF;
				}
				else
				{
					Fact = 0.; FactIm = 0.; 
				}
			}
			*(t++) = Fact;
			*(t++) = FactIm;
		}
		else
		{
			*(t++) = 0.; *(t++) = 0.;
		}

		q += qStep;
	}
	return 0;
}

//*************************************************************************

int srTRadIntPeriodic::SetupConvolutionData_Tapered(int n, float* ConvFactorData, double eStart, double eFin, long Ne)
{
	const double MinArgExp = -20.; // To steer

	double eStep = (eFin - eStart)/(Ne - 1);
	double eStartFict = -eStep*(Ne >> 1);

	float *AuxCont = new float[Ne << 1];
	if(AuxCont == 0) return MEMORY_ALLOCATION_FAILURE;

	double Buf1 = TwoPI*n/EnAvg_FinNper;
	double Buf2 = 0.5*MagPer.TaperPar_TU/(NperGen*NperGen);
	double Buf3 = n/(NperGen*EnAvg_FinNper);

	const double c0 = 1.239854E-09;
	double Buf0 = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double HalfFundamentalE = Buf0;

	double EnAvg_mi_HalfFun = EnAvg_FinNper - HalfFundamentalE;
	double EnAvg_p_HalfFun = EnAvg_FinNper + HalfFundamentalE;

	double e = eStartFict + EnAvg_FinNper;
	float *t = AuxCont;
	for(long i=0; i<Ne; i++)
	{
		if((e < EnAvg_mi_HalfFun) || (e > EnAvg_p_HalfFun))
		{
			*(t++) = 0.;
		}
		else
		{
			double SumC = 0., SumS = 0.;
			double Buf1e = Buf1*e;
			for(int k=0; k<NperGen; k++)
			{
				double Arg = Buf1e*k*(1 + Buf2*(k - NperGen + 1));
				double CosArg, SinArg;
				CosAndSin(Arg, CosArg, SinArg);
				SumC += CosArg; SumS += SinArg;
			}

			*(t++) = (float)(Buf3*(SumC*SumC + SumS*SumS));
		}
		*(t++) = 0.;
		e += eStep;
	}

  	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.xStep = eStep;
	FFT1DInfo.xStart = eStartFict;
	FFT1DInfo.Nx = Ne;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
	FFT1DInfo.pInData = AuxCont;
	FFT1DInfo.pOutData = ConvFactorData;
	CGenMathFFT1D FFT1D;

	int result;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;
	if(AuxCont != 0) delete AuxCont;

	const double TwoPIe2 = TwoPI*PI;
	double SigmaE_Gamma = 4*n*c0*EbmDat.SigmaRelE/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double MultArgEnSpr = -TwoPIe2*SigmaE_Gamma*SigmaE_Gamma;

	double q = FFT1DInfo.xStartTr;
	double qStep = FFT1DInfo.xStepTr;

	t = ConvFactorData;
	for(long j=0; j<Ne; j++)
	{
		double ArgExp = MultArgEnSpr*q*q;
		if(ArgExp > MinArgExp)
		{
			float Fact = (float)exp(ArgExp);
			*t *= Fact;
			*(t + 1) *= Fact;
		}
		else
		{
			*t = 0.;
			*(t + 1) = 0.;
		}
		t += 2;
		q += qStep;
	}
	return 0;
}

//*************************************************************************

double srTRadIntPeriodic::EstimateTaperResCurveWidth(int n)
{
	const double c0 = 1.239854E-09;
	double HalfEnergyFund = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double CritEnergyForHarm = (n << 1)*HalfEnergyFund;

	double Buf1 = TwoPI*n/CritEnergyForHarm;
	double Buf2 = 0.5*MagPer.TaperPar_TU/(NperGen*NperGen);

	const int AmOfTrends = 10; // To steer
	double RelTol = 0.1; // To steer

	double eMin = CritEnergyForHarm;
	//double eMax = CritEnergyForHarm + HalfEnergyFund;
	double e = eMin, eRange = HalfEnergyFund;
	double Fc;

	for(int i=0; i<AmOfTrends; i++)
	{
		double SumC = 0., SumS = 0.;
		double Buf1e = Buf1*e;
		for(int k=0; k<NperGen; k++)
		{
			double Arg = Buf1e*k*(1 + Buf2*(k - NperGen + 1));
			double CosArg, SinArg;
			CosAndSin(Arg, CosArg, SinArg);
			SumC += CosArg; SumS += SinArg;
		}

		double F = SumC*SumC + SumS*SumS;
		if(i == 0) 
		{
			Fc = F;
			e += eRange;
		}
		else
		{
			double Fn = F/Fc;
			if(Fn > RelTol)
			{
				if(i == 1) return 2.*HalfEnergyFund;
				else
				{
					e += eRange;
				}
			}
			else
			{
				e -= eRange;
			}
		}
		eRange *= 0.5;
	}
	return 2.*(e - eMin);
}

//*************************************************************************

int srTRadIntPeriodic::DeduceGridOverPhotonEnergyAndAzimuth(int n, double& eStart, double& eFin, long& ne, srTEnergyAzimuthGrid& EnAzGrid)
{
	double eMinEff, eMaxEff, PhiMinEff, PhiMaxEff;
	EstimateEnergyAndPhiObsLimits(n, eMinEff, eMaxEff, PhiMinEff, PhiMaxEff);

	//OCTEST 150713
	//double eStartTot = eStart, eFinTot = eFin;
	//if(eMinEff > eStartTot) eMinEff = eStartTot; //??
	//END OCTEST

	EnAzGrid.eMinEff = eMinEff;
	EnAzGrid.eMaxEff = eMaxEff;
	EnAzGrid.PhiMin = PhiMinEff;
	EnAzGrid.PhiMax = PhiMaxEff;

	const double c0 = 1.239854E-09;
	double Buf0 = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double CritEnergyForHarm = (n << 1)*Buf0;
	double FundE = CritEnergyForHarm/double(n); //OC170713
	EnAzGrid.eCritForHarm = CritEnergyForHarm;

	double &PhiLenToResolve = EnAzGrid.PhiLenToResolve, &eStepToResolve = EnAzGrid.eStepToResolve;
	eStepToResolve = 1.E+23;
	EnAzGrid.eStepToResolveObsPixels = 1.E+23;

	double eStepToResolveObsPixels;
	EstimateEnergyStepAndPhiLenToResolveObsPixels(n, eStepToResolveObsPixels, PhiLenToResolve);
	if(EnAzGrid.EnsureEnResolvingObsPixels) 
	{
		eStepToResolve = eStepToResolveObsPixels;
		EnAzGrid.eStepToResolveObsPixels = eStepToResolveObsPixels;
	}

	if((MagPer.TypeOfUnd > 0) || (EnergySpreadShouldBeTreated))
	{
		EstimateEnergyStepToResolveFinNper(n, EnAzGrid); // Don't call it at any diff. place, only after EstimateEnergyStepAndPhiLenToResolveObsPixels !!!
	}

	if((!EnergySpreadShouldBeTreated) && (MagPer.TypeOfUnd == 0)) 
	{
		EnAzGrid.Ne = ne; EnAzGrid.eStart = eStart; EnAzGrid.eFin = eFin;
		return SetUpVariableGridOverAzimuth(n, EnAzGrid);
	}

// In the case of finite Nper and Energy Spread:
// - Make such a grid that resolves well all the observation pixels (by the proper Teta angles).
// - Take care of the edges of the energy interval.

	char VeryLargeEnergySpreadAndFinNperContr = 0;
	double EnergySpreadContr = 0., FiniteNperContr = 0., EnergySpreadAndFinNperContr = 0.;
	if(EnergySpreadShouldBeTreated)
	{
		double SigmaEnGam = (n << 2)*EbmDat.SigmaRelE*Buf0;
		EnergySpreadContr = 2.35*SigmaEnGam;
	}
	if(MagPer.TypeOfUnd == 1) // Conventional undulator
	{
		//FiniteNperContr = eMaxEff/(n*NperGen); //Modify
		FiniteNperContr = CritEnergyForHarm/(n*NperGen); //OC
	}
	else if(MagPer.TypeOfUnd == 2) // Tapered
	{
		FiniteNperContr = EstimateTaperResCurveWidth(n);

		//double FundE = CritEnergyForHarm/double(n); //OC170713
		if(FiniteNperContr > 0.7*FundE) // To steer
		{
			VeryLargeEnergySpreadAndFinNperContr = 1;
		}

		if(FiniteNperContr > 0.8*FundE) // To steer
		{
			//srTSend Send; Send.AddWarningMessage(pWarningsGen, TOO_LARGE_TAPER_PARAM);
			CErrWarn::AddWarningMessage(pWarningsGen, TOO_LARGE_TAPER_PARAM);
		}
	}
	else if(MagPer.TypeOfUnd == 3) // Optical Klystron
	{
		FiniteNperContr = (1. + 0.5*MagPer.PhaseSh_OK)*eMaxEff/(n*NperGen); // Modify

		if(FiniteNperContr > 0.7*CritEnergyForHarm/double(n)) // To steer
		{
			VeryLargeEnergySpreadAndFinNperContr = 1;
		}
	}

	if(EnergySpreadShouldBeTreated && (MagPer.TypeOfUnd > 0))
	{
		EnergySpreadAndFinNperContr = sqrt(EnergySpreadContr*EnergySpreadContr + FiniteNperContr*FiniteNperContr);
	}
	else if(EnergySpreadShouldBeTreated) EnergySpreadAndFinNperContr = EnergySpreadContr;
	else EnergySpreadAndFinNperContr = FiniteNperContr;

	int AmOfExtraResWidths = 4; // To steer !
	//int AmOfExtraResWidths = 6; //OCTEST To steer !

	int MinAmOfPoints = 20; // To steer !

	if(VeryLargeEnergySpreadAndFinNperContr)
	{
		AmOfExtraResWidths = 1;
	}

	char ActuallyOnlyOnePointIsNeeded = 0;
	int ExtraOnePoLeft = 0, ExtraOnePoRight = 0;
	double eStep, eExtraLoc;
	if(ne == 1)
	{
		ActuallyOnlyOnePointIsNeeded = 1;

		//MinAmOfPoints = 20; // To steer !
		const int PointsPerFringe = 4; // To steer !
		ne = (AmOfExtraResWidths*PointsPerFringe) << 1;
		if(ne < MinAmOfPoints) ne = MinAmOfPoints;

		double EnAvgLoc = DistrInfoDat.LambStart;
		eExtraLoc = AmOfExtraResWidths*EnergySpreadAndFinNperContr;

		int neToResolve = int(2*eExtraLoc/eStepToResolve) + 1;
		if(ne < neToResolve) 
		{
			ne = neToResolve;
		}
		eStart = EnAvgLoc - eExtraLoc;
		eStep = eExtraLoc/(ne >> 1);

		long ir = long(::fabs(CritEnergyForHarm - EnAvgLoc)/eStep);
		if(ir > 0)
		{
			eStep = ::fabs(CritEnergyForHarm - EnAvgLoc)/double(ir);
		}
		else
		{
			if(::fabs(CritEnergyForHarm - EnAvgLoc) < 0.25*eStep) // To steer
			{
				EnAvgLoc = CritEnergyForHarm;
			}
			else
			{
				double eStepOld = eStep;
				eStep = ::fabs(CritEnergyForHarm - EnAvgLoc);
				ne *= long(eStepOld/eStep);
			}
		}

		eExtraLoc = eStep*(ne >> 1);
		eStart = EnAvgLoc - eExtraLoc;

		ExtraOnePoLeft = ne >> 1;
		ExtraOnePoRight = (ne >> 1) - 1;
	}
	else
	{
		eExtraLoc = AmOfExtraResWidths*EnergySpreadAndFinNperContr;
		eStep = (eFin - eStart)/(ne - 1);
		if(eStepToResolve < eStep) 
		{
			ne = int((eFin - eStart)/eStepToResolve) + 1;
			eStep = eStepToResolve;
		}

		double dAmOfExtraPointsFromLeft = eExtraLoc/eStep;
		long AmOfExtraPointsFromLeft = long(dAmOfExtraPointsFromLeft);
		if((dAmOfExtraPointsFromLeft - AmOfExtraPointsFromLeft) > 0.5*eStep) AmOfExtraPointsFromLeft++;
		if(AmOfExtraPointsFromLeft == 0) AmOfExtraPointsFromLeft = 1;

		long AmOfExtraPointsFromRight = AmOfExtraPointsFromLeft;
		
		if((eStart - eExtraLoc) > eMinEff) 
		{
			double dAmOfExtraPointsFromLeft = (eStart - eMinEff)/eStep;
			AmOfExtraPointsFromLeft = long(dAmOfExtraPointsFromLeft);
			if((dAmOfExtraPointsFromLeft - AmOfExtraPointsFromLeft) > 0.5*eStep) AmOfExtraPointsFromLeft++;
		}

		//double NextHarmPeak = ((n + 1) << 1)*Buf0; //OC
		//double eExtraLocRight = NextHarmPeak - eMaxEff;
		//eExtraLocRight = (eExtraLocRight > eExtraLoc)? eExtraLocRight : eExtraLoc;
		//double dAmOfExtraPointsFromRight = eExtraLocRight/eStep;
		//AmOfExtraPointsFromRight = long(dAmOfExtraPointsFromRight);
		//if((dAmOfExtraPointsFromRight - AmOfExtraPointsFromRight) > 0.5*eStep) AmOfExtraPointsFromRight++;
		
		//if((eFin + eExtraLocRight) < eMaxEff) 
		if((eFin + eExtraLoc) < eMaxEff) 
		{
			double dAmOfExtraPointsFromRight = (eMaxEff - eFin)/eStep;
			AmOfExtraPointsFromRight = long(dAmOfExtraPointsFromRight);
			if((dAmOfExtraPointsFromRight - AmOfExtraPointsFromRight) > 0.5*eStep) AmOfExtraPointsFromRight++;
		}

		eStart = DistrInfoDat.LambStart - AmOfExtraPointsFromLeft*eStep;
		ne += (AmOfExtraPointsFromLeft + AmOfExtraPointsFromRight);

		CorrectGridForPassingThroughCritEnergy(n, eStart, eStep, ne);
	}

	double eStartOld = eStart;
	double eFinOld = eStart + eStep*(ne - 1);
	//long neOld = ne;

	double eEffAvg = 0.5*(eMinEff + eMaxEff);
	double eFinApprox = eStart + (ne - 1)*eStep;
	if((eStart < eMinEff) && (eFinApprox > eEffAvg))
	{
		long NeLoc = ne;
		for(long i=0; i<NeLoc; i++)
		{
			double eStartBuf = eStart + eStep;
			if((eStartBuf < eMinEff) && (ne > MinAmOfPoints))
			{
				eStart = eStartBuf; ne--;
			}
			else break;
		}
	}
	double FractionOfStep = 0.01*eStep; // To steer
	eFinApprox = eStart + (ne - 1)*eStep;
	if(eFinApprox > eMaxEff + FractionOfStep)
	{
		long NeLoc = ne;
		for(long i=0; i<NeLoc; i++)
		{
			double eFinBuf = eFinApprox - eStep;
			if((eFinBuf > eMaxEff + FractionOfStep) && (ne > MinAmOfPoints))
			{
				eFinApprox = eFinBuf; ne--;
			}
			else break;
		}

		if(ne > MinAmOfPoints)
		{
			ne--; eFinApprox -= eStep;
		}
	}

	eFin = eStart + (ne - 1)*eStep;

	EnAzGrid.Ne = ne;
	EnAzGrid.eStart = eStart;
	EnAzGrid.eFin = eFin;

	long neExtraLeft = 0, neExtraRight = 0;
	double eExtraLocNew = 2*eExtraLoc; // To steer

	long nLeft = long(0.5*eStart/Buf0); //OC
	double PrevHarmPeak = (nLeft >= 1)? (nLeft << 1)*Buf0 : Buf0;
	//double PrevHarmPeak = (n > 1)? ((n - 1) << 1)*Buf0 : Buf0;
	double NextHarmPeak = ((n + 1) << 1)*Buf0;

	double eExtraLocNewLeft = (PrevHarmPeak < eStart)? (eStart - PrevHarmPeak) : 0;
	//if(eExtraLocNewLeft < eExtraLocNew) eExtraLocNewLeft = eExtraLocNew;

	double eExtraLocNewRight = (eFin < NextHarmPeak)? (NextHarmPeak - eFin) : 0;
	if(eExtraLocNewRight < eExtraLocNew) eExtraLocNewRight = eExtraLocNew;

	if(eStartOld < eStart)
	{
		double Buf = eStart - eStartOld;
		//if(Buf > eExtraLocNew) neExtraLeft = long(eExtraLocNew/eStep); //OC
		if(Buf > eExtraLocNewLeft) neExtraLeft = long(eExtraLocNewLeft/eStep);
		else neExtraLeft = long(Buf/eStep);
	}
	if(eFin < eFinOld)
	{
		double Buf = eFinOld - eFin;
		//if(Buf > eExtraLocNew) neExtraRight = long(eExtraLocNew/eStep); //OC
		if(Buf > eExtraLocNewRight) neExtraRight = long(eExtraLocNewRight/eStep);
		else neExtraRight = long(Buf/eStep);
	}

	EnAzGrid.NeExtraLeft = neExtraLeft;

	//if((EnAzGrid.eMinEff < EnAzGrid.eCritForHarm) && (EnAzGrid.eMaxEff*1.00001 >= EnAzGrid.eCritForHarm))
	//{// To make the sharp edge approximately in the middle
	//	long BufExtraRight = long((EnAzGrid.Ne + EnAzGrid.NeExtraLeft)*0.6);
	//	if(neExtraRight < BufExtraRight) neExtraRight = BufExtraRight;
	//}
	//OC
	//if((EnAzGrid.eMaxEff*1.00001 - EnAzGrid.eCritForHarm)*4. < (EnAzGrid.eCritForHarm - EnAzGrid.eMinEff)) //To steer
	//{// To make the sharp edge approximately in the middle
	//	double AuxEmax = (2.*EnAzGrid.eCritForHarm - (EnAzGrid.eMinEff + neExtraLeft*eStep));
	//	if(AuxEmax > EnAzGrid.eMaxEff)
	//	{
	//		double AuxRat = (AuxEmax - EnAzGrid.eMaxEff)/(EnAzGrid.eMaxEff - EnAzGrid.eMinEff);
	//		long BufExtraRight = long(AuxRat*EnAzGrid.Ne);
	//		if(neExtraRight < BufExtraRight) neExtraRight = BufExtraRight;
	//	}
	//}

	//OCTEST 150713 (commented-out the section below)
	//if((eFin*1.00001 - EnAzGrid.eCritForHarm)*4. < (EnAzGrid.eCritForHarm - eStart)) //To steer
	//{// To make the sharp edge approximately in the middle
	//	double AuxEmax = (2.*EnAzGrid.eCritForHarm - (eStart + neExtraLeft*eStep));
	//	if(AuxEmax > eFin)
	//	{
	//		double AuxRat = (AuxEmax - eFin)/(eFin - eStart);
	//		long BufExtraRight = long(AuxRat*EnAzGrid.Ne);
	//		if(neExtraRight < BufExtraRight) neExtraRight = BufExtraRight;
	//	}
	//}

	//OC170713
	//Attempt to move artificial discontuinities (because of convolution of harmonics 
	//for taking into account finite und. length and Nper) to higher photon energies
	double eMaxExtentRight = CritEnergyForHarm + FundE*IntPerStoPrec.MinPhotEnExtRight;
	if(EnAzGrid.eFin + neExtraRight*eStep < eMaxExtentRight)
	{
		neExtraRight = long((eMaxExtentRight - EnAzGrid.eFin)/eStep);
	}

	//OCTEST 150713
	//if(EnAzGrid.eFin + neExtraRight*eStep < 3*eFinTot)
	//{//To make sure that "sharp edge" (harmonic calculation cut-off) never happens
	//	neExtraRight = long((3*eFinTot - EnAzGrid.eFin)/eStep);
	//}
	//END OCTEST

	EnAzGrid.NeExtraRight = neExtraRight;

	long Ne = EnAzGrid.NeExtraLeft + EnAzGrid.Ne + EnAzGrid.NeExtraRight;
	long NeOld = Ne;
	CGenMathFFT1D FFT1D;
	FFT1D.NextCorrectNumberForFFT(Ne);
	long ExtraNe = Ne - NeOld;
	EnAzGrid.NeExtraRight += ExtraNe;

	if(!ActuallyOnlyOnePointIsNeeded) CorrectGridToAllowRangeResizeOnTheOtherSide(EnAzGrid);
	else CorrectGridForOnePoint(EnAzGrid);

	return SetUpVariableGridOverAzimuth(n, EnAzGrid);
}

//*************************************************************************

void srTRadIntPeriodic::CorrectGridForOnePoint(srTEnergyAzimuthGrid& EnAzGrid)
{
		double EnAvgLoc = DistrInfoDat.LambStart;
		double eStep = (EnAzGrid.eFin - EnAzGrid.eStart)/(EnAzGrid.Ne - 1);

		char SigRight = (EnAvgLoc > EnAzGrid.eStart)? 1 : -1;
		long AmPoLeft = EnAzGrid.NeExtraLeft + SigRight*long(::fabs(EnAvgLoc - EnAzGrid.eStart)/eStep + 1.E-06);

		char SigLeft = (EnAzGrid.eFin > EnAvgLoc)? 1 : -1;
		long AmPoRight = EnAzGrid.NeExtraRight + SigLeft*long(::fabs(EnAzGrid.eFin - EnAvgLoc)/eStep + 1.E-06);

		EnAzGrid.EnergyAdjustFinGridPar.eShift = 0.;
		EnAzGrid.EnergyAdjustFinGridPar.AmOfExtraStepsFromLeftAndRightOnTheOtherSide = 0;
		EnAzGrid.EnergyAdjustFinGridPar.Multip = 1;
		EnAzGrid.EnergyAdjustFinGridPar.AmStExtraLeft = -AmPoLeft;
		EnAzGrid.EnergyAdjustFinGridPar.AmStExtraRight = -AmPoRight;
		EnAzGrid.EnergyAdjustFinGridPar.AmOfInitPoToIgnore = 0;
}

//*************************************************************************

int srTRadIntPeriodic::SetUpVariableGridOverAzimuth(int n, srTEnergyAzimuthGrid& EnAzGrid)
{
	int result = EnAzGrid.AllocateAmOfPointsArray();
	if(result) return result;
	DeduceNphi(n); NphiGen++;

	const double c0 = 1.239854E-09;
	double be0 = (n << 1)*c0/MagPer.PerLength;
	double be1 = EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2);
	double MinEnergy = EnAzGrid.eMinEff;
	double ArgSqrt = be0/MinEnergy - be1;
	double MaxTet = (ArgSqrt > 0.)? sqrt(ArgSqrt) : 0.;
	int NphiNonUniform = int(MaxTet*TwoPI/EnAzGrid.PhiLenToResolve);
	if(((NphiNonUniform >> 1) << 1) != NphiNonUniform) NphiNonUniform++;
	int NphiLoc = (NphiNonUniform > NphiGen)? NphiNonUniform : NphiGen;

	NphiGen = NphiLoc;
	EnAzGrid.MakeAllAzPointsConstant(NphiGen);

	return 0;
}

//*************************************************************************

void srTRadIntPeriodic::EstimateEnergyStepToResolveFinNper(int n, srTEnergyAzimuthGrid& EnAzGrid)
{
	double qMax = 0.;

	double EnVal = (DistrInfoDat.nLamb > 1)? EnAzGrid.eMinEff : DistrInfoDat.LambStart;

	if(MagPer.TypeOfUnd == 1) // normal
	{
		const double OverResolveFactor = 1.05; // To steer
		qMax = OverResolveFactor*n*NperGen/EnVal;
	}
	else if(MagPer.TypeOfUnd == 2) // tapered
	{
		//Almost same as normal:
		const double OverResolveFactor = 0.5; // To steer
		qMax = OverResolveFactor*n*NperGen/EnVal;
	}
	else if(MagPer.TypeOfUnd == 3) // optical klystron
	{
		const double OverResolveFactor = 1.05; // To steer
		qMax = (OverResolveFactor*n*NperGen/EnVal)*(1 + 0.5*MagPer.PhaseSh_OK);
	}

	const double c0 = 1.239854E-09;
	if(EnergySpreadShouldBeTreated)
	{
		const double LocSigmaFact = 3.; //To steer
		double Buf0 = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
		double SigmaEnGam = (n << 2)*EbmDat.SigmaRelE*Buf0;

		double qMaxEnSpr = LocSigmaFact/(TwoPI*SigmaEnGam);
		if(qMaxEnSpr < qMax) qMax = qMaxEnSpr; // This can reduce fr. range, i.e. make e-step larger
	}

	if(EnAzGrid.eMaxEff < EnAzGrid.eCritForHarm) // Off-axis case: one can make the step larger
	{
		const double LocSigmaFact = 4.; //To steer
		double RangeE = 0.5*(EnAzGrid.eMaxEff - EnAzGrid.eMinEff);

		double qMaxInfOffAx = LocSigmaFact*LocSigmaFact/(TwoPI*RangeE);
		if(qMaxInfOffAx < qMax) qMax = qMaxInfOffAx; // This can reduce fr. range, i.e. make e-step larger
	}

	double eStepToResolveLoc = 0.5/qMax;
	EnAzGrid.eStepToResolveFinNperAndEnSpr = eStepToResolveLoc;
	if(eStepToResolveLoc < EnAzGrid.eStepToResolve) EnAzGrid.eStepToResolve = eStepToResolveLoc;
}

//*************************************************************************

void srTRadIntPeriodic::EstimateEnergyStepAndPhiLenToResolveObsPixels(int n, double& eStepToResolve, double& PhiLenToResolve)
{
	//double AmOfPointsPerPixelE = 7; // To steer
	double AmOfPointsPerPixelE = 14; // To steer

	double AmOfPointsPerPixelPhi = 6; // To steer
	//const double SigFact = 1.4; // To steer
	const double SigFact = 0.5; // To steer

	double xEffHalfPixelSize = dxHalfSlit + SigFact*SigmaXpGen;
	double zEffHalfPixelSize = dzHalfSlit + SigFact*SigmaZpGen;

	double tx1, tz1, tx0, tz0;
	//FindMostOffsetPixelCenter(tx, tz);
	FindLeastAndMostOffsetPixelCenters(tx0, tz0, tx1, tz1);
	double txMin = tx1 - xEffHalfPixelSize, txMax = tx1 + xEffHalfPixelSize, tzMin = tz1 - zEffHalfPixelSize, tzMax = tz1 + zEffHalfPixelSize;
	double TetMinE2_MostOffset, TetMaxE2_MostOffset;
	FindTetMinMaxE2_FromTetxTetz(txMin, txMax, tzMin, tzMax, TetMinE2_MostOffset, TetMaxE2_MostOffset);

	double txMin0 = tx0 - xEffHalfPixelSize, txMax0 = tx0 + xEffHalfPixelSize, tzMin0 = tz0 - zEffHalfPixelSize, tzMax0 = tz0 + zEffHalfPixelSize;
	double TetMinE2_LeastOffset, TetMaxE2_LeastOffset;
	FindTetMinMaxE2_FromTetxTetz(txMin0, txMax0, tzMin0, tzMax0, TetMinE2_LeastOffset, TetMaxE2_LeastOffset);

	double kSharpMax = 9.; //To steer
	double kSharpSigmaTrans = 0.08; //To steer
	double kSharpSigmaTransE2 = kSharpSigmaTrans*kSharpSigmaTrans;

	double SharpRatX = 1.15*SigmaXpGen/dxHalfSlit;
	double SharpRatZ = 1.15*SigmaZpGen/dzHalfSlit;
	double SharpRatMin = (SharpRatX < SharpRatZ)? SharpRatX : SharpRatZ;
	double SharpArgExp = -0.5*SharpRatMin*SharpRatMin/kSharpSigmaTransE2;
	double kSharpFactor = 1. + (kSharpMax - 1.)*exp(SharpArgExp);
	AmOfPointsPerPixelPhi *= kSharpFactor;

	AmOfPointsPerPixelE *= (1. + 1.*(kSharpFactor - 1.)); // ??

	const double c0 = 1.239854E-09;
	double Buf1 = (n << 1)*c0;
	double Buf2 = EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2);

	double eMin0 = Buf1/(MagPer.PerLength*(TetMaxE2_LeastOffset + Buf2));
	double eMax0 = Buf1/(MagPer.PerLength*(TetMinE2_LeastOffset + Buf2));
	double eStepToResolve0 = (eMax0 - eMin0)/(AmOfPointsPerPixelE - 1);

	double eMin1 = Buf1/(MagPer.PerLength*(TetMaxE2_MostOffset + Buf2));
	double eMax1 = Buf1/(MagPer.PerLength*(TetMinE2_MostOffset + Buf2));
	double eStepToResolve1 = (eMax1 - eMin1)/(AmOfPointsPerPixelE - 1);

	eStepToResolve = eStepToResolve0;
	if(eStepToResolve > eStepToResolve1) eStepToResolve = eStepToResolve1;

	double TetAvg, PhiInterval = 0.;
	EstimateRadAvgOfBox(txMin, txMax, tzMin, tzMax, TetAvg);

	const double ShiftFactor = 1.E-09; // To steer
	double SmallShiftX = ShiftFactor*(txMax - txMin);
	double SmallShiftZ = ShiftFactor*(tzMax - tzMin);

	if(txMin == 0.) txMin += SmallShiftX;
	if(txMax == 0.) txMax -= SmallShiftX;
	if(tzMin == 0.) tzMin += SmallShiftZ;
	if(tzMax == 0.) tzMax -= SmallShiftZ;

	double txAm, txAp, tzAm, tzAp;
	char SplitByVerLine = 0, SplitByHorLine = 0;
	if(txMin*txMax < 0.)
	{
		txAm = -SmallShiftX; txAp = SmallShiftX; SplitByVerLine = 1;
	}
	if(tzMin*tzMax < 0.)
	{
		tzAm = -SmallShiftZ; tzAp = SmallShiftZ; SplitByHorLine = 1;
	}

	if(SplitByVerLine)
	{
		if(SplitByHorLine)
		{
			PhiInterval += PhiIntToResolveBox(txMin, txAm, tzMin, tzAm, TetAvg);
			PhiInterval += PhiIntToResolveBox(txAp, txMax, tzMin, tzAm, TetAvg);
			PhiInterval += PhiIntToResolveBox(txMin, txAm, tzAp, tzMax, TetAvg);
			PhiInterval += PhiIntToResolveBox(txAp, txMax, tzAp, tzMax, TetAvg);
		}
		else
		{
			PhiInterval += PhiIntToResolveBox(txMin, txAm, tzMin, tzMax, TetAvg);
			PhiInterval += PhiIntToResolveBox(txAp, txMax, tzMin, tzMax, TetAvg);
		}
	}
	else
	{
		if(SplitByHorLine)
		{
			PhiInterval += PhiIntToResolveBox(txMin, txMax, tzMin, tzAm, TetAvg);
			PhiInterval += PhiIntToResolveBox(txMin, txMax, tzAp, tzMin, TetAvg);
		}
		else
		{
			PhiInterval += PhiIntToResolveBox(txMin, txMax, tzMin, tzMax, TetAvg);
		}
	}
	PhiLenToResolve = TetAvg*PhiInterval/(IntPerStoPrec.Knphi*(AmOfPointsPerPixelPhi - 1));
}

//*************************************************************************

double srTRadIntPeriodic::PhiIntToResolveBox(double x1, double x2, double y1, double y2, double rAvg)
{//Assumes rectangle not to cross the axes; So |PhiInt|<HalfPI.

	double x1x1 = x1*x1, x2x2 = x2*x2, y1y1 = y1*y1, y2y2 = y2*y2;
	double R11e2 = x1x1 + y1y1, R12e2 = x1x1 + y2y2, R21e2 = x2x2 + y1y1, R22e2 = x2x2 + y2y2;
	double re2 = rAvg*rAvg;

	double BufHorLine1 = (R11e2 - re2)*(R21e2 - re2);
	char InrsctHorLine1 = (BufHorLine1 <= 0.);

	double BufHorLine2 = (R12e2 - re2)*(R22e2 - re2);
	char InrsctHorLine2 = (BufHorLine2 <= 0.);

	double BufVerLine1 = (R11e2 - re2)*(R12e2 - re2);
	char InrsctVerLine1 = (BufVerLine1 <= 0.);

	double BufVerLine2 = (R21e2 - re2)*(R22e2 - re2);
	char InrsctVerLine2 = (BufVerLine2 <= 0.);

	if(!(InrsctHorLine1 || InrsctHorLine2 || InrsctVerLine1 || InrsctVerLine2)) return 0.;

	double Phi1 = 1.E+23, Phi2 = 1.E+23;
	double *pPhi = &Phi1;
	if(InrsctHorLine1) *pPhi = AzimuthOfArcAndHorLineIntersect(y1, rAvg, 0.5*(x1 + x2));
	
	pPhi = (Phi1 == 1.E+23)? &Phi1 : &Phi2;
	if(InrsctVerLine1) *pPhi = AzimuthOfArcAndVerLineIntersect(x1, rAvg, 0.5*(y1 + y2));

	pPhi = (Phi1 == 1.E+23)? &Phi1 : &Phi2;
	if(InrsctHorLine2) *pPhi = AzimuthOfArcAndHorLineIntersect(y2, rAvg, 0.5*(x1 + x2));

	pPhi = (Phi1 == 1.E+23)? &Phi1 : &Phi2;
	if(InrsctVerLine2) *pPhi = AzimuthOfArcAndVerLineIntersect(x2, rAvg, 0.5*(y1 + y2));

	if((Phi1 == 1.E+23) || (Phi2 == 1.E+23)) return 0.;

	double PhiInt = ::fabs(Phi2 - Phi1);
	if(PhiInt > HalfPI) PhiInt = TwoPI - PhiInt;
	return PhiInt;
}

//*************************************************************************

void srTRadIntPeriodic::EstimateEnergyAndPhiObsLimits(int n, double& eMinEff, double& eMaxEff, double& PhiMinEff, double& PhiMaxEff)
{
	double txMin, txMax, tzMin, tzMax;
	FindObservationLimits(txMin, txMax, tzMin, tzMax);
	char PixelsCoverZeroPoint = ((txMin < 0.) && (txMax > 0.) && (tzMin < 0.) && (tzMax > 0.));

	double txMinE2 = txMin*txMin, txMaxE2 = txMax*txMax, tzMinE2 = tzMin*tzMin, tzMaxE2 = tzMax*tzMax;
	double Tet00e2 = txMinE2 + tzMinE2, Tet01e2 = txMinE2 + tzMaxE2, Tet10e2 = txMaxE2 + tzMinE2, Tet11e2 = txMaxE2 + tzMaxE2;
	double TetMaxE2 = Tet00e2;
	if(Tet01e2 > TetMaxE2) TetMaxE2 = Tet01e2;
	if(Tet10e2 > TetMaxE2) TetMaxE2 = Tet10e2;
	if(Tet11e2 > TetMaxE2) TetMaxE2 = Tet11e2;

	double TetMinE2 = Tet00e2;
	if(PixelsCoverZeroPoint) TetMinE2 = 0.;
	else if((txMin < 0.) && (txMax > 0.))
	{
		if(tzMinE2 < TetMinE2) TetMinE2 = tzMinE2;
		if(tzMaxE2 < TetMinE2) TetMinE2 = tzMaxE2;
	}
	else if((tzMin < 0.) && (tzMax > 0.))
	{
		if(txMinE2 < TetMinE2) TetMinE2 = txMinE2;
		if(txMaxE2 < TetMinE2) TetMinE2 = txMaxE2;
	}
	else
	{
		if(Tet01e2 < TetMinE2) TetMinE2 = Tet01e2;
		if(Tet10e2 < TetMinE2) TetMinE2 = Tet10e2;
		if(Tet11e2 < TetMinE2) TetMinE2 = Tet11e2;
	}

	const double c0 = 1.239854E-09;
	double Buf1 = (n << 1)*c0;
	double Buf2 = EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2);
	eMinEff = Buf1/(MagPer.PerLength*(TetMaxE2 + Buf2));
	eMaxEff = Buf1/(MagPer.PerLength*(TetMinE2 + Buf2));

	//double Buf3 = c0/(MagPer.PerLength*Buf2); //OC
	//if(n > 1)
	//{
	//	double ePeakPrevHarm = ((n - 1) << 1)*Buf3;
	//	if(eMinEff > ePeakPrevHarm) eMinEff = ePeakPrevHarm;
	//}

	//double ePeakNextHarm = ((n + 0.6)*2.)*Buf3;
	//if(eMaxEff < ePeakNextHarm) eMaxEff = ePeakNextHarm;

	if(PixelsCoverZeroPoint)
	{
		PhiMinEff = 0.; PhiMaxEff = TwoPI;
	}
	else
	{
		TVector2d p0(txMin, tzMin), p1(txMin, tzMax), p2(txMax, tzMin), p3(txMax, tzMax);
		TVector2d Vectors[] = {p0, p1, p2, p3};
		FindPhiIntervalForVectors(Vectors, 4, PhiMinEff, PhiMaxEff);
	}
}

//*************************************************************************

void srTRadIntPeriodic::FindPhiIntervalForVectors(TVector2d* Vectors, int LenVectors, double& PhiMin, double& PhiMax)
{
	TVector2d *pMin = Vectors, *pMax = Vectors;
	for(int i=1; i<LenVectors; i++)
	{
		TVector2d *p = Vectors + i;
		double TestMin = pMin->x*p->y - p->x*pMin->y;
		if(TestMin < 0.) pMin = p;
		double TestMax = pMax->x*p->y - p->x*pMax->y;
		if(TestMax > 0.) pMax = p;
	}

	PhiMin = FormalPhase((float)(pMin->x), (float)(pMin->y)); 
	if(PhiMin < 0.) PhiMin += TwoPI;
	PhiMax = FormalPhase((float)(pMax->x), (float)(pMax->y)); 
	if(PhiMax < 0.) PhiMax += TwoPI;
	if(PhiMax < PhiMin) PhiMax += TwoPI;
}

//*************************************************************************

void srTRadIntPeriodic::CorrectGridForPassingThroughCritEnergy(int n, double& eStart, double& eStep, long& ne)
{
	if(ne < 20) return; // Steer this to avoid fluctuations

	const double c0 = 1.239854E-09;
	double Buf0 = c0/(MagPer.PerLength*EbmDat.GammaEm2*(1. + MagPer.HalfKxE2pKzE2));
	double CritEnergyForHarm = (n << 1)*Buf0;

	double eFin = eStart + (ne - 1)*eStep;
	if(!((eStart < CritEnergyForHarm) && (CritEnergyForHarm < eFin))) return;

	long ir = (long)((CritEnergyForHarm - eStart)/eStep);
	double GridMisfit = CritEnergyForHarm - (eStep*ir + eStart);
	if(GridMisfit > 0.5*eStep) GridMisfit = -(eStep - GridMisfit);

	double GridShiftTol = 1.E-06*eStep; // To steer
	if(!(::fabs(GridMisfit) > GridShiftTol)) return;

	eStart += GridMisfit;
}

//*************************************************************************

int srTRadIntPeriodic::Int1D_Simpson(double xSt, double xFi, long Nx, char VsSorPhi, srTEFourier& Res)
{
	int result;
	if(((Nx >> 1) << 1) == Nx) Nx++;
	double Step = (xFi - xSt)/(Nx - 1);

	srTEFourier F, Edges, Sum1, Sum2;

	if(VsSorPhi == 's')
	{
		Fs(xSt, 0, Edges); Fs(xFi, Nx - 1, F);
	}
	else if(VsSorPhi == 'p')
	{
		Fphi(xSt, 0, Edges); 
		if(::fabs(xFi - TwoPI) < 1.E-03*Step) F = Edges;
		else Fphi(xFi, Nx - 1, F);
	}
	Edges.EwX_Re += F.EwX_Re; Edges.EwX_Im += F.EwX_Im; Edges.EwZ_Re += F.EwZ_Re; Edges.EwZ_Im += F.EwZ_Im; 

	double x = xSt + Step;
	long i_gen = 1;

	if(VsSorPhi == 's')
	{
		for(long i=0; i<((Nx-3)>>1); i++)
		{
			Fs(x, i_gen++, F); x += Step;
			Sum1.EwX_Re += F.EwX_Re; Sum1.EwX_Im += F.EwX_Im; Sum1.EwZ_Re += F.EwZ_Re; Sum1.EwZ_Im += F.EwZ_Im; 
			Fs(x, i_gen++, F); x += Step;
			Sum2.EwX_Re += F.EwX_Re; Sum2.EwX_Im += F.EwX_Im; Sum2.EwZ_Re += F.EwZ_Re; Sum2.EwZ_Im += F.EwZ_Im;
		}
		Fs(x, i_gen, F);
	}
	else
	{
		for(long i=0; i<((Nx-3)>>1); i++)
		{
			Fphi(x, i_gen++, F); x += Step;
			Sum1.EwX_Re += F.EwX_Re; Sum1.EwX_Im += F.EwX_Im; Sum1.EwZ_Re += F.EwZ_Re; Sum1.EwZ_Im += F.EwZ_Im; 
			Fphi(x, i_gen++, F); x += Step;
			Sum2.EwX_Re += F.EwX_Re; Sum2.EwX_Im += F.EwX_Im; Sum2.EwZ_Re += F.EwZ_Re; Sum2.EwZ_Im += F.EwZ_Im;
		}
		Fphi(x, i_gen, F);
	}

	Sum1.EwX_Re += F.EwX_Re; Sum1.EwX_Im += F.EwX_Im; Sum1.EwZ_Re += F.EwZ_Re; Sum1.EwZ_Im += F.EwZ_Im;

	double Mult = 0.3333333333*Step;
	Res.EwX_Re = Mult*(Edges.EwX_Re + 4.*Sum1.EwX_Re + 2.*Sum2.EwX_Re);
	Res.EwX_Im = Mult*(Edges.EwX_Im + 4.*Sum1.EwX_Im + 2.*Sum2.EwX_Im);
	Res.EwZ_Re = Mult*(Edges.EwZ_Re + 4.*Sum1.EwZ_Re + 2.*Sum2.EwZ_Re);
	Res.EwZ_Im = Mult*(Edges.EwZ_Im + 4.*Sum1.EwZ_Im + 2.*Sum2.EwZ_Im);

	result = srYield.Check();
	return result;
}

//*************************************************************************

