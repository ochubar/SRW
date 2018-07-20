/************************************************************************//**
 * File: srradint.cpp
 * Description: SR calculation method from ~Arbitrary Transversely-Uniform Magnetic Field, in Near-Field observation conditions
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifdef __IGOR_PRO__
#include "srigintr.h"
#endif

#include "srradint.h"
#include "srprgind.h"
#include "srinterf.h"
#include "srmlttsk.h"
#include "sroptelm.h"
#include "srerror.h"

//*************************************************************************

extern srTYield srYield;

extern srTIntVect gVectWarnNos;
extern char* srWarningDynamic;

//*************************************************************************

void srTRadInt::Initialize()
{
	DistrInfoDat.Initialize();

	TrjDatPtr = 0;
	IntegSubdArray = 0; LenIntegSubdArray = 0;
	RadDistrPhotonFluxHorPol = RadDistrPhotonFluxVerPol = 0; 
	RadDistrFieldFourierHorPol = RadDistrFieldFourierVerPol = 0;
	IntegratedPhotonFlux = 0;
	AuxPhaseArray = 0;
	ComputeDer = 0;

	BtxArr = XArr = IntBtxE2Arr = BtzArr = ZArr = IntBtzE2Arr = BxArr = BzArr = 0;

	PI = 3.141592653590;
	TwoPI = 2.*PI;
	ThreePIdTwo = 1.5*PI;
	FivePIdFour = 1.25*PI;

	HalfPI = 0.5* PI;
	One_dTwoPI = 1./TwoPI;
	PIm10e6 = PI*1.E+06;
	PIm10e6dEnCon = PIm10e6*0.80654658;
	TenEm6dTwoPI = 1.E-06/TwoPI;

	a2c = -0.5; a4c = 1./24.; a6c = -1./720.; a8c = 1./40320.; a10c = -1./3628800.; a12c = 1./479001600.;
	a3s = -1./6.; a5s = 1./120.; a7s = -1./5040.; a9s = 1./362880.; a11s = -1./39916800.; a13s = 1./6227020800.;

	NumberOfLevelsFilled = 0;
	for(int k=0; k<50; k++) 
	{
		AmOfPointsOnLevel[k] = 0; BtxArrP[k] = 0;
	}
	MaxFluxDensVal = CurrentAbsPrec = 0.;

	TrjDataContShouldBeRebuild = 1;
	ProbablyTheSameLoop = 0;
	MaxMemAvail = 8.;
	MaxNumPoToSave = 200000;

	TryToApplyNearFieldResidual = 0;
    TryToApplyNearFieldResidual_AtRight = 0;

	MaxLevelForMeth_10_11 = 18; //15; //OC170815
	UseManualSlower = 0;

	ComputeNormalDerivative = 0;
	dExdlRe = dExdlIm = dEzdlRe = dEzdlIm = 0;
	dExdlReTravers = dExdlImTravers = dEzdlReTravers = dEzdlImTravers = 0;

	EstimatedAbsoluteTolerance = 0.;
	DistrInfoDat.CoordUnits = 1; // To ensure mm for coord.

	m_CalcResidTerminTerms = 1; // Do calculate residual terminating terms by default
}

//*************************************************************************

int srTRadInt::CheckInputConsistency()
{
	//if(TrjDatPtr->EbmDat.Gamma < 10.) return ELEC_BEAM_IS_NOT_ULTRARELATIVISTIC;
	if(TrjDatPtr->EbmDat.Gamma <= 1.) return ELEC_BEAM_IS_NOT_ULTRARELATIVISTIC;
	if(TrjDatPtr->EbmDat.Gamma <= 10.)
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, WARN_ELEC_BEAM_IS_NOT_ULTRARELATIVISTIC);
	}

	if((TrjDatPtr->EbmDat.s0 < TrjDatPtr->sStart) || (TrjDatPtr->EbmDat.s0 > TrjDatPtr->sStart + TrjDatPtr->sStep*(TrjDatPtr->LenFieldData - 1))) return S0_OUT_OF_FIELD_DEFINITION_LIMITS;

	if(DistrInfoDat.yStart < sIntegFin) return LONGIT_COORD_OF_OBS_POINT_WITHIN_RAD_INTEGR_LIMITS;

	//double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239854E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239842E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	if((DistrInfoDat.yStart - sIntegFin) < 3.*WavelengthIn_m) return WAVELENGTH_COMPARABLE_WITH_OBS_DISTANCE;

	double Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataDerivedAtPoint(sIntegFin, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z);
	//double Dx1 = ::fabs(DistrInfoDat.xStart - Crdx), Dx2 = ::fabs(DistrInfoDat.xEnd - Crdx);
	//double MaxDx = (Dx1 > Dx2)? Dx1 : Dx2;
	//double MaxNx = MaxDx/(DistrInfoDat.yStart - sIntegFin);
	//double Dz1 = ::fabs(DistrInfoDat.zStart - Crdz), Dz2 = ::fabs(DistrInfoDat.zEnd - Crdz);
	//double MaxDz = (Dz1 > Dz2)? Dz1 : Dz2;
	//double MaxNz = MaxDz/(DistrInfoDat.yStart - sIntegFin);
	//if((MaxNx > 0.4) || (MaxNz > 0.4)) return INSTANT_OBS_ANGLES_TOO_LARGE;
	//removed 31/08/2003

// Make more tests here

	return 0;
}

//*************************************************************************

int srTRadInt::AllocateMemForRadDistr()
{
	if(DistrInfoDat.OnlyOnePoint)
	{
		if((DistrInfoDat.DistrValType == StokesParam) || (DistrInfoDat.DistrValType == FieldFourier))
		{
			RadDistrFieldFourierHorPol = SmallContComplex;
			RadDistrFieldFourierHorPolTravers = RadDistrFieldFourierHorPol;

			RadDistrFieldFourierVerPol = SmallContComplex + 2;
			RadDistrFieldFourierVerPolTravers = RadDistrFieldFourierVerPol;

			if(ComputeNormalDerivative)
			{
				dExdlRe = dExdlReTravers = SmallContDoubleForNormDer;
				dExdlIm = dExdlImTravers = SmallContDoubleForNormDer + 1;
				dEzdlRe = dEzdlReTravers = SmallContDoubleForNormDer + 2;
				dEzdlIm = dEzdlImTravers = SmallContDoubleForNormDer + 3;
			}
			return 0;
		}
	}

	int TotAmOfData = DistrInfoDat.nLamb * DistrInfoDat.ny * DistrInfoDat.nz * DistrInfoDat.nx;

	if((DistrInfoDat.DistrValType == StokesParam) || (DistrInfoDat.DistrValType == FieldFourier))
	{
		char DeletionMayBeNeeded = ((RadDistrFieldFourierHorPol != 0) && (RadDistrFieldFourierHorPol != SmallContComplex));

		if(DeletionMayBeNeeded) delete[] RadDistrFieldFourierHorPol;
		RadDistrFieldFourierHorPol = new complex<double>[TotAmOfData];
		//if(RadDistrFieldFourierHorPol==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		if(RadDistrFieldFourierHorPol==0) { return MEMORY_ALLOCATION_FAILURE;}
		RadDistrFieldFourierHorPolTravers = RadDistrFieldFourierHorPol;

		if(RadDistrFieldFourierVerPol != 0) if(DeletionMayBeNeeded) delete[] RadDistrFieldFourierVerPol;
		RadDistrFieldFourierVerPol = new complex<double>[TotAmOfData];
		//if(RadDistrFieldFourierVerPol==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
		if(RadDistrFieldFourierVerPol==0) { return MEMORY_ALLOCATION_FAILURE;}
		RadDistrFieldFourierVerPolTravers = RadDistrFieldFourierVerPol;

		if(ComputeNormalDerivative)
		{
			if(dExdlRe != 0) if(DeletionMayBeNeeded) delete[] dExdlRe;
			dExdlRe = new double[TotAmOfData];
			//if(dExdlRe==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
			if(dExdlRe==0) { return MEMORY_ALLOCATION_FAILURE;}
			dExdlReTravers = dExdlRe;

			if(dExdlIm != 0) if(DeletionMayBeNeeded) delete[] dExdlIm;
			dExdlIm = new double[TotAmOfData];
			//if(dExdlIm==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
			if(dExdlIm==0) { return MEMORY_ALLOCATION_FAILURE;}
			dExdlImTravers = dExdlIm;

			if(dEzdlRe != 0) if(DeletionMayBeNeeded) delete[] dEzdlRe;
			dEzdlRe = new double[TotAmOfData];
			//if(dEzdlRe==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
			if(dEzdlRe==0) { return MEMORY_ALLOCATION_FAILURE;}
			dEzdlReTravers = dEzdlRe;

			if(dEzdlIm != 0) if(DeletionMayBeNeeded) delete[] dEzdlIm;
			dEzdlIm = new double[TotAmOfData];
			//if(dEzdlIm==0) { pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;}
			if(dEzdlIm==0) { return MEMORY_ALLOCATION_FAILURE;}
			dEzdlImTravers = dEzdlIm;
		}
		return 0;
	}
	return -1;
}

//*************************************************************************

int srTRadInt::ComputeTotalRadDistrLoops()
{
	int result;
	double StepLambda = (DistrInfoDat.nLamb > 1)? (DistrInfoDat.LambEnd - DistrInfoDat.LambStart)/(DistrInfoDat.nLamb - 1) : 0.;
	double StepX = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
	double StepY = (DistrInfoDat.ny > 1)? (DistrInfoDat.yEnd - DistrInfoDat.yStart)/(DistrInfoDat.ny - 1) : 0.;
	double StepZ = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
	InitializeTraverses();

	double *pStartArgLoop[4], *pArgLoop[4], *pStepArgLoop[4];
	int *pNArgLoop[4];

	for(int i=0; i<4; i++)
	{
		char LoopID = DistrInfoDat.LoopOrder[i];
		if(LoopID == 'y') 
		{ 
			pStartArgLoop[i] = &(DistrInfoDat.yStart);
			pStepArgLoop[i] = &StepY;
			pArgLoop[i] = &(ObsCoor.y);
			pNArgLoop[i] = &(DistrInfoDat.ny);
		}
		else if(LoopID == 'w') 
		{ 
			pStartArgLoop[i] = &(DistrInfoDat.LambStart);
			pStepArgLoop[i] = &StepLambda;
			pArgLoop[i] = &(ObsCoor.Lamb);
			pNArgLoop[i] = &(DistrInfoDat.nLamb);
		}
		else if(LoopID == 'x') 
		{ 
			pStartArgLoop[i] = &(DistrInfoDat.xStart);
			pStepArgLoop[i] = &StepX;
			pArgLoop[i] = &(ObsCoor.x);
			pNArgLoop[i] = &(DistrInfoDat.nx);
		}
		else if(LoopID == 'z') 
		{ 
			pStartArgLoop[i] = &(DistrInfoDat.zStart);
			pStepArgLoop[i] = &StepZ;
			pArgLoop[i] = &(ObsCoor.z);
			pNArgLoop[i] = &(DistrInfoDat.nz);
		}
	}
	*(*pArgLoop) = *(*pStartArgLoop);
	for(int iLoop0 = 0; iLoop0 < *(pNArgLoop[0]); iLoop0++)
	{
		*(pArgLoop[1]) = *(pStartArgLoop[1]);
		for(int iLoop1 = 0; iLoop1 < *(pNArgLoop[1]); iLoop1++)
		{
			*(pArgLoop[2]) = *(pStartArgLoop[2]);
			for(int iLoop2 = 0; iLoop2 < *(pNArgLoop[2]); iLoop2++)
			{
				*(pArgLoop[3]) = *(pStartArgLoop[3]);
				for(int iLoop3 = 0; iLoop3 < *(pNArgLoop[3]); iLoop3++)
				{
					complex<double> RadIntegValues[2];
					srTEFourier EwNormDer;

					if(result = GenRadIntegration(RadIntegValues, &EwNormDer)) return result;

					*(RadDistrFieldFourierHorPolTravers++) = *RadIntegValues;
					*(RadDistrFieldFourierVerPolTravers++) = *(RadIntegValues + 1);

					if(ComputeNormalDerivative)
					{
						*(dExdlReTravers++) = EwNormDer.EwX_Re;
						*(dExdlImTravers++) = EwNormDer.EwX_Im;
						*(dEzdlReTravers++) = EwNormDer.EwZ_Re;
						*(dEzdlImTravers++) = EwNormDer.EwZ_Im;
					}

					*(pArgLoop[3]) += *(pStepArgLoop[3]);
				}
				*(pArgLoop[2]) += *(pStepArgLoop[2]);
			}
			*(pArgLoop[1]) += *(pStepArgLoop[1]);
		}
		*(*pArgLoop) += *(*pStepArgLoop);
	}
	return 0;
}

//*************************************************************************

int srTRadInt::ComputeTotalRadDistrDirectOut(srTSRWRadStructAccessData& SRWRadStructAccessData, char showProgressInd)
{
	int result = 0;
	EstimateAbsoluteTolerance();
	ProbablyTheSameLoop = 1;

	SRWRadStructAccessData.UnderSamplingX = SRWRadStructAccessData.UnderSamplingZ = 1; //OC290805
	DeallocateMemForRadDistr(); // Do not remove this

	//if(result = pSend->InitRadDistrOutFormat3(SRWRadStructAccessData, DistrInfoDat)) return result;
	//?????

	if(result = SetupRadCompStructures()) return result;
	if(DistrInfoDat.ShowPhaseOnly) return ScanPhase();

	double StepLambda = (DistrInfoDat.nLamb > 1)? (DistrInfoDat.LambEnd - DistrInfoDat.LambStart)/(DistrInfoDat.nLamb - 1) : 0.;
	double StepX = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
	double StepZ = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;

	//long PerX = DistrInfoDat.nLamb << 1;
	//long PerZ = DistrInfoDat.nx*PerX;
	long long PerX = DistrInfoDat.nLamb << 1;
	long long PerZ = DistrInfoDat.nx*PerX;

	float *pEx0 = SRWRadStructAccessData.pBaseRadX;
	float *pEz0 = SRWRadStructAccessData.pBaseRadZ;

	complex<double> RadIntegValues[2];
	srTEFourier EwNormDer;

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	AnalizeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ);

	double xc = TrjDatPtr->EbmDat.x0;
	double zc = TrjDatPtr->EbmDat.z0;
	double xTol = StepX*0.001, zTol = StepZ*0.001; // To steer

	//long TotalAmOfOutPointsForInd = DistrInfoDat.nz*DistrInfoDat.nx*DistrInfoDat.nLamb;
	long TotalAmOfOutPointsForInd = ((long long)DistrInfoDat.nz)*((long long)DistrInfoDat.nx)*((long long)DistrInfoDat.nLamb);
	if(FinalResAreSymOverX) TotalAmOfOutPointsForInd >>= 1;
	if(FinalResAreSymOverZ) TotalAmOfOutPointsForInd >>= 1;
	//long PointCount = 0;
	long long PointCount = 0;
	double UpdateTimeInt_s = 0.5;
	//srTCompProgressIndicator* pCompProgressInd = 0;
	//if(showProgressInd) pCompProgressInd = new srTCompProgressIndicator(TotalAmOfOutPoints, UpdateTimeInt_s);

	if(!showProgressInd) TotalAmOfOutPointsForInd = 0;
	srTCompProgressIndicator compProgressInd(TotalAmOfOutPointsForInd, UpdateTimeInt_s);

	//long AbsPtCount = 0;
	ObsCoor.y = DistrInfoDat.yStart;
	ObsCoor.z = DistrInfoDat.zStart;
	for(int iz=0; iz<DistrInfoDat.nz; iz++)
	{
		if(FinalResAreSymOverZ) { if((ObsCoor.z - zc) > zTol) break;}

		//long izPerZ = iz*PerZ;
		long long izPerZ = iz*PerZ;
		ObsCoor.x = DistrInfoDat.xStart;
		for(int ix=0; ix<DistrInfoDat.nx; ix++)
		{
			if(FinalResAreSymOverX) { if((ObsCoor.x - xc) > xTol) break;}

			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;
			ObsCoor.Lamb = DistrInfoDat.LambStart;
			for(int iLamb=0; iLamb<DistrInfoDat.nLamb; iLamb++)
			{
				if(result = GenRadIntegration(RadIntegValues, &EwNormDer)) return result;

				//long Offset = izPerZ + ixPerX + (iLamb << 1);
				long long Offset = izPerZ + ixPerX + (iLamb << 1);
				float *pEx = pEx0 + Offset, *pEz = pEz0 + Offset;

				*pEx = float(RadIntegValues->real());
				*(pEx+1) = float(RadIntegValues->imag());
				*pEz = float(RadIntegValues[1].real());
				*(pEz+1) = float(RadIntegValues[1].imag());

				if(showProgressInd) 
				{
                    //if(result = pCompProgressInd->UpdateIndicator(PointCount++)) return result;
                    if(result = compProgressInd.UpdateIndicator(PointCount++)) return result;
				}
				if(result = srYield.Check()) return result;

				ObsCoor.Lamb += StepLambda;
			}
			ObsCoor.x += StepX;
		}
		ObsCoor.z += StepZ;
	}

	if(FinalResAreSymOverZ || FinalResAreSymOverX) 
		FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, SRWRadStructAccessData);

	//if(showProgressInd && (pCompProgressInd != 0)) 
	//{
	//	delete pCompProgressInd; pCompProgressInd = 0;
	//}

	//pSend->FinishRadDistrOutFormat1();
	//????
	return result;
}

//*************************************************************************

int srTRadInt::ComputeNormalResidual(double s, int NumberOfTerms, complex<double>* ResidValues, srTEFourier* pEwNormDerResid)
{
// Steerable parameters
	const double CritFarFieldTermRatio = 0.35;
	//const double CritNearFieldTermRatio = 0.6;

// End Steerable parameters

	if((NumberOfTerms < 1) || (NumberOfTerms > 3)) return -1;
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;

	double dBxds=0., dBzds=0., Bx=0., Bz=0., Btx=0., Btz=0., Crdx=0., Crdz=0., IntBtE2x=0., IntBtE2z=0.;
	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('x', s, dBzds, Bz, Btx, Crdx, IntBtE2x);
	TrjDatPtr->CompTrjDataAndFieldWithDerAtPoint('z', s, dBxds, Bx, Btz, Crdz, IntBtE2z);

	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double Ax, Az, dAxds, dAzds, d2Axds2, d2Azds2, Ph, dPhds, d2Phds2, d3Phds3, ConBtxBz, ConBtzBx, CosPh, SinPh;
	double ConBtxBzpAx, ConBtzBxpAz, TwoPIm10e9_d_Lamb;
	complex<double> gExpIPhase;

	double AxNormDer, AzNormDer, dAxdsNormDer, dAzdsNormDer;
	double Ny, L_x_N, t2xNormDer = 0., t2zNormDer = 0.;

	double t1xd =0., t2xd =0., t3xd =0., t1zd =0., t2zd =0., t3zd =0.;
	double d2Phds2_d_dPhdsE2_ForExpansionTest;

	complex<double> Zero(0., 0.);

// Angles assumed in mr!
	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		double One_d_ymis = 1./(ObsCoor.y - s);
		double xObs_mi_x = xObs - Crdx, zObs_mi_z = zObs - Crdz;

		double LongTerm = IntBtE2x + IntBtE2z;
		//double a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis;
		//double a = a0*One_d_ymis;
		//Ph = PIm10e9_d_Lamb*(s*TrjDatPtr->EbmDat.GammaEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
		//OC_test
        	double a0 = LongTerm + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*TrjDatPtr->EbmDat.GammaEm2 + a0);
		//end OC_TEST

		CosAndSin(Ph, CosPh, SinPh);

		complex<double> ExpIPhase(CosPh, SinPh);
		gExpIPhase = ExpIPhase;

		double Nx = xObs_mi_x*One_d_ymis, Btx_mi_Nx = Btx - Nx;
		double Nz = zObs_mi_z*One_d_ymis, Btz_mi_Nz = Btz - Nz;
		dPhds = PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
		complex<double> IdPhds(0., dPhds);
		double One_d_dPhds = 1./dPhds;
		complex<double> mOne_d_IdPhds(0., One_d_dPhds);

		Ax = Btx_mi_Nx*One_d_ymis;
		Az = Btz_mi_Nz*One_d_ymis;
		//complex<double> PreExpX = Ax;
		//complex<double> PreExpZ = Az;
		complex<double> PreExpX = Ax, PreExpX01 = Ax, PreExpX02 = Ax;
		complex<double> PreExpZ = Az, PreExpZ01 = Az, PreExpZ02 = Az;

		t1xd = Ax; t1zd = Az;

		if(ComputeNormalDerivative)
		{
			Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
			L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
			AxNormDer = Ax*L_x_N;
			AzNormDer = Az*L_x_N;
		}

		if(NumberOfTerms > 1)
		{
			ConBtxBz = (TrjDatPtr->BetaNormConst)*Bz;
			ConBtzBx = (-TrjDatPtr->BetaNormConst)*Bx;
			ConBtxBzpAx = ConBtxBz + Ax;
			ConBtzBxpAz = ConBtzBx + Az;
			TwoPIm10e9_d_Lamb = 2.*PIm10e9_d_Lamb;

			d2Phds2 = TwoPIm10e9_d_Lamb*(Btx_mi_Nx*ConBtxBzpAx + Btz_mi_Nz*ConBtzBxpAz);
			dAxds = (2.*Ax + ConBtxBz)*One_d_ymis;
			dAzds = (2.*Az + ConBtzBx)*One_d_ymis;

			double d2Phds2_d_dPhds = d2Phds2*One_d_dPhds;
			t2xd = (dAxds - Ax*d2Phds2_d_dPhds)*One_d_dPhds; 
			t2zd = (dAzds - Az*d2Phds2_d_dPhds)*One_d_dPhds;

			d2Phds2_d_dPhdsE2_ForExpansionTest = d2Phds2_d_dPhds*One_d_dPhds;

			complex<double> t2x(0., t2xd);
			PreExpX += t2x;
			PreExpX02 = PreExpX;

			complex<double> t2z(0., t2zd);
			PreExpZ += t2z;
			PreExpZ02 = PreExpZ;

			if(ComputeNormalDerivative)
			{
				double dL_x_Nds = -Ax*(SurfNorm.x - SurfNorm.y*(Nx/Ny)) - Az*(SurfNorm.z - SurfNorm.y*(Nz/Ny));
				dAxdsNormDer = L_x_N*dAxds + Ax*dL_x_Nds;
				dAzdsNormDer = L_x_N*dAzds + Az*dL_x_Nds;
				t2xNormDer = (dAxdsNormDer - AxNormDer*d2Phds2_d_dPhds)*One_d_dPhds;
				t2zNormDer = (dAzdsNormDer - AzNormDer*d2Phds2_d_dPhds)*One_d_dPhds;
			}
		}
		if(NumberOfTerms > 2)
		{
			//double dPhdsE2 = dPhds*dPhds;
			double ConBtxdBzds = (TrjDatPtr->BetaNormConst)*dBzds;
			double ConBtzdBxds = (-TrjDatPtr->BetaNormConst)*dBxds;

			d3Phds3 = TwoPIm10e9_d_Lamb*(ConBtxBzpAx*ConBtxBzpAx + ConBtzBxpAz*ConBtzBxpAz 
				+ Btx_mi_Nx*(One_d_ymis*(2.*Btx_mi_Nx*One_d_ymis + ConBtxBz) + ConBtxdBzds) 
				+ Btz_mi_Nz*(One_d_ymis*(2.*Btz_mi_Nz*One_d_ymis + ConBtzBx) + ConBtzdBxds));

			d2Axds2 = One_d_ymis*(3.*dAxds + ConBtxdBzds);
			d2Azds2 = One_d_ymis*(3.*dAzds + ConBtzdBxds);

			double One_d_dPhdsE2 = One_d_dPhds*One_d_dPhds;
			double d2Phds2_mu_d2Phds2_d_dPhdsE2 = d2Phds2*d2Phds2*One_d_dPhdsE2;

			t3xd = (-d2Axds2 + (3.*dAxds*d2Phds2 + Ax*d3Phds3)*One_d_dPhds - 3.*Ax*d2Phds2_mu_d2Phds2_d_dPhdsE2)*One_d_dPhdsE2; 
			t3zd = (-d2Azds2 + (3.*dAzds*d2Phds2 + Az*d3Phds3)*One_d_dPhds - 3.*Az*d2Phds2_mu_d2Phds2_d_dPhdsE2)*One_d_dPhdsE2;

			complex<double> t3x(t3xd, 0.);
			PreExpX += t3x;
			complex<double> t3z(t3zd, 0.);
			PreExpZ += t3z;
		}

		bool FarFieldExpansionWithTwoTermsIsPracticallyOK = false, FarFieldExpansionWithThreeTermsIsPracticallyOK = false;
		if(NumberOfTerms > 1)
		{
            FarFieldExpansionWithTwoTermsIsPracticallyOK = 
				((t2xd == 0.) && (t1xd == 0.)) || ((t2xd != 0.) && (::fabs(t2xd) < 0.7*(::fabs(t1xd))));
		}
		if(NumberOfTerms > 2)
		{
            FarFieldExpansionWithThreeTermsIsPracticallyOK = FarFieldExpansionWithTwoTermsIsPracticallyOK &&
				(((t3xd == 0.) && (t2xd == 0.)) || ((t3xd != 0.) && (::fabs(t3xd) < 0.7*(::fabs(t2xd)))));
		}

		//char PracticalProblemWithFarFieldExpansion = (NumberOfTerms > 2) && 
		//	((::fabs(t3xd) > ::fabs(t2xd)) && (::fabs(t2xd) > ::fabs(t1xd)));
		//char FarFieldExpansionIsPracticallyOK = (NumberOfTerms > 2) && 
		//	((::fabs(t3xd) < 0.7*(::fabs(t2xd))) && (::fabs(t2xd) < 0.7*(::fabs(t1xd))));

		//char PhaseDerivCriterionIsSatisfied = ::fabs(d2Phds2_d_dPhdsE2_ForExpansionTest) < CritFarFieldTermRatio;

		char FarFieldExpansDoesNotWork = false;
		//if(FarFieldExpansionIsPracticallyOK) FarFieldExpansDoesNotWork = false;
		//else
		//{
		//	if(PracticalProblemWithFarFieldExpansion) FarFieldExpansDoesNotWork = true;
		//	else FarFieldExpansDoesNotWork = !PhaseDerivCriterionIsSatisfied;
		//}

		if(NumberOfTerms > 2) //OC180105
		{
			if(!FarFieldExpansionWithThreeTermsIsPracticallyOK)
			{
				if(!FarFieldExpansionWithTwoTermsIsPracticallyOK)
				{
                    FarFieldExpansDoesNotWork = true;
                    PreExpX = Zero;
                    PreExpZ = Zero;
				}
				else
				{
                    FarFieldExpansDoesNotWork = false;
                    PreExpX = PreExpX02;
                    PreExpZ = PreExpZ02;
				}
			}
			else FarFieldExpansDoesNotWork = false;
		}
		else if(NumberOfTerms > 1)
		{
            if(!FarFieldExpansionWithTwoTermsIsPracticallyOK)
			{
                FarFieldExpansDoesNotWork = true;
                PreExpX = Zero;
                PreExpZ = Zero;
			}
			else FarFieldExpansDoesNotWork = false;
		}

		if(FarFieldExpansDoesNotWork && TryToApplyNearFieldResidual && (::fabs(Bx) < TrjDatPtr->FieldZeroTolerance) && (::fabs(Bz) < TrjDatPtr->FieldZeroTolerance))
		{
			if(ComputeDer == 1) // Left residual
			{
				double sStNew1 = s - 2.*(PI/(PIm10e9_d_Lamb*TrjDatPtr->EbmDat.GammaEm2)); // Steer 3 or ...
				double sStNew2 = 0.5*(ObsCoor.y - sqrt(ObsCoor.y*ObsCoor.y + 4.*(xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)/(TrjDatPtr->EbmDat.GammaEm2 + Btx*Btx + Btz*Btz))); // Steer ...
				double sStNew = (sStNew1 < sStNew2)? sStNew1 : sStNew2;

				double One_d_ymisNew = 1./(ObsCoor.y - sStNew);
				double xObs_mi_xNew = xObs - X_FreeSpace(sStNew, s, Btx, Crdx);
				double zObs_mi_zNew = zObs - X_FreeSpace(sStNew, s, Btz, Crdz);
				double Btx_mi_NxNew = Btx - xObs_mi_xNew*One_d_ymisNew;
				double Btz_mi_NzNew = Btz - zObs_mi_zNew*One_d_ymisNew;
				double dPhdsNew = PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Btx_mi_NxNew*Btx_mi_NxNew + Btz_mi_NzNew*Btz_mi_NzNew);
				double Two_d_ymisNewE2 = 2.*One_d_ymisNew*One_d_ymisNew;
				double AxNew = Btx_mi_NxNew*One_d_ymisNew, AzNew = Btz_mi_NzNew*One_d_ymisNew;
				double dAxdsNew = Btx_mi_NxNew*Two_d_ymisNewE2;
				double dAzdsNew = Btz_mi_NzNew*Two_d_ymisNewE2;
				double AxdPhdsNew = AxNew*dPhdsNew, AzdPhdsNew = AzNew*dPhdsNew;
				double d2Phds2New = TwoPIm10e9_d_Lamb*(Btx_mi_NxNew*AxNew + Btz_mi_NzNew*AzNew);
				double One_d_dPhdsNew = 1./dPhdsNew;
				double d2Phds2_d_dPhdsNew = d2Phds2New*One_d_dPhdsNew;

				double PhNew = PIm10e9_d_Lamb*(sStNew*TrjDatPtr->EbmDat.GammaEm2 + (IntBtE2x + IntBtE2z) + (Btx*Btx + Btz*Btz)*(sStNew - s) + (xObs_mi_xNew*xObs_mi_xNew + zObs_mi_zNew*zObs_mi_zNew)*One_d_ymisNew);
				double CosPhNew, SinPhNew;
				CosAndSin(PhNew, CosPhNew, SinPhNew);

				double DiffReX = (dAxdsNew*CosPhNew - AxdPhdsNew*SinPhNew) - (dAxds*CosPh - Ax*dPhds*SinPh);
				double DiffImX = (AxdPhdsNew*CosPhNew + dAxdsNew*SinPhNew) - (Ax*dPhds*CosPh + dAxds*SinPh);
				double DiffReZ = (dAzdsNew*CosPhNew - AzdPhdsNew*SinPhNew) - (dAzds*CosPh - Az*dPhds*SinPh);
				double DiffImZ = (AzdPhdsNew*CosPhNew + dAzdsNew*SinPhNew) - (Az*dPhds*CosPh + dAzds*SinPh);

				double PreExpXNewRe = (AxNew*d2Phds2_d_dPhdsNew - dAxdsNew)*One_d_dPhdsNew*One_d_dPhdsNew;
				double PreExpXNewIm = AxNew*One_d_dPhdsNew;
				double PreExpZNewRe = (AzNew*d2Phds2_d_dPhdsNew - dAzdsNew)*One_d_dPhdsNew*One_d_dPhdsNew;
				double PreExpZNewIm = AzNew*One_d_dPhdsNew;

				if((::fabs(PreExpXNewIm) < ::fabs(PreExpXNewRe)) || (::fabs(PreExpZNewIm) < ::fabs(PreExpZNewRe))) 
				{
					//return BAD_RAD_INT_TERMINATIONS;

					complex<double> Zero(0., 0.);
					*ResidValues = Zero; *(ResidValues+1) = Zero;

					CErrWarn::AddWarningMessage(&gVectWarnNos, COMPUTATION_WAS_NOT_DONE_BAD_RAD_INT_TERMINATIONS);
					//OC_hack 15_03_2005
					return 0;
				}

				double ActNormConstCosPhNew = ActNormConst*CosPhNew;
				double ActNormConstSinPhNew = ActNormConst*SinPhNew;
					
				srTEFourier EwLeftResid(PreExpXNewRe*ActNormConstCosPhNew - PreExpXNewIm*ActNormConstSinPhNew,
										PreExpXNewIm*ActNormConstCosPhNew + PreExpXNewRe*ActNormConstSinPhNew,
										PreExpZNewRe*ActNormConstCosPhNew - PreExpZNewIm*ActNormConstSinPhNew,
										PreExpZNewIm*ActNormConstCosPhNew + PreExpZNewRe*ActNormConstSinPhNew);
				srTEFourier DiffEdgeDer(DiffReX, DiffImX, DiffReZ, DiffImZ);
				ComputePreResid(sStNew, s, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z, &DiffEdgeDer, &EwLeftResid, 0);

				complex<double> ResidX(EwLeftResid.EwX_Re, EwLeftResid.EwX_Im), ResidZ(EwLeftResid.EwZ_Re, EwLeftResid.EwZ_Im);
				*ResidValues = ResidX;
				*(ResidValues+1) = ResidZ;
				return 0;
			}
			else // Right residual
			{
				if(TryToApplyNearFieldResidual_AtRight)
				{
                    srTEFourier EwRightResid = ComputeRightPartOfRadIntegralInNearField(s, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z);

// Revision is needed here !!!
/**
				double xEff = xObs - (Crdx + Btx*(ObsCoor.y - s));
				double zEff = zObs - (Crdz + Btz*(ObsCoor.y - s));

				double TwoLambda = TwoPI/PIm10e9_d_Lamb;
				if((::fabs(xEff) < TwoLambda) && (::fabs(zEff) < TwoLambda)) 
				{
					complex<double> Zero(0., 0.);
					*ResidValues = Zero; *(ResidValues+1) = Zero;

					pSend->AddWarningMessage(&gVectWarnNos, COMPUTATION_WAS_NOT_DONE_EL_TRAJ_TOO_CLOSE_TO_OBSERV);
					return 0;
				}

				char PreResidNeeded = 1;
				double sFiNew = (ObsCoor.y - PIm10e9_d_Lamb*(xEff*xEff + zEff*zEff));
				if(sFiNew < s) { sFiNew = s; PreResidNeeded = 0;}

				double yObs_mi_sNew = ObsCoor.y - sFiNew;
				double yObsmisNewE2 = yObs_mi_sNew*yObs_mi_sNew;
				double xObs_mi_xNew = xObs - (Crdx + Btx*(sFiNew - s));
				double zObs_mi_zNew = zObs - (Crdz + Btz*(sFiNew - s));
				double dPhdt = PIm10e9_d_Lamb*(xObs_mi_xNew*xObs_mi_xNew + zObs_mi_zNew*zObs_mi_zNew - 2.*yObs_mi_sNew*(Btx*xObs_mi_xNew + Btz*zObs_mi_zNew) + yObsmisNewE2*(TrjDatPtr->EbmDat.GammaEm2 + Btx*Btx + Btz*Btz));
				double d2Phdt2 = PIm10e9_d_Lamb*2.*(TrjDatPtr->EbmDat.GammaEm2)*yObsmisNewE2*yObs_mi_sNew;

				double Axt = yObs_mi_sNew*Btx - xObs_mi_xNew;
				double Azt = yObs_mi_sNew*Btz - zObs_mi_zNew;

				double dPhdtE2 = dPhdt*dPhdt;
				if((::fabs(d2Phdt2) >= CritNearFieldTermRatio*dPhdtE2) || (dPhdt == 0.)) 
				{
					//return BAD_RAD_INT_TERMINATIONS;

					complex<double> Zero(0., 0.);
					*ResidValues = Zero; *(ResidValues+1) = Zero;

					pSend->AddWarningMessage(&gVectWarnNos, COMPUTATION_WAS_NOT_DONE_BAD_RAD_INT_TERMINATIONS);
					return 0;
				}

				double One_d_dPhdt = 1./dPhdt;
				double d2Phdt2_d_dPhdtE2 = d2Phdt2*One_d_dPhdt*One_d_dPhdt;

				double One_d_yObsmisNew = 1./yObs_mi_sNew;
				double PhNew = PIm10e9_d_Lamb*(sFiNew*TrjDatPtr->EbmDat.GammaEm2 + (IntBtE2x + IntBtE2z) + (Btx*Btx + Btz*Btz)*(sFiNew - s) + (xObs_mi_xNew*xObs_mi_xNew + zObs_mi_zNew*zObs_mi_zNew)*One_d_yObsmisNew);
				double CosPhNew, SinPhNew;
				CosAndSin(PhNew, CosPhNew, SinPhNew);
				double BufRe = -SinPhNew + d2Phdt2_d_dPhdtE2*CosPhNew;
				double BufIm = CosPhNew + d2Phdt2_d_dPhdtE2*SinPhNew;

				double ActNormConst_d_dPhdt = ActNormConst*One_d_dPhdt;
				double ActNormConstAxt_d_dPhdt = Axt*ActNormConst_d_dPhdt;
				double ActNormConstAzt_d_dPhdt = Azt*ActNormConst_d_dPhdt;

				srTEFourier EwRightResid(ActNormConstAxt_d_dPhdt*BufRe, ActNormConstAxt_d_dPhdt*BufIm,
										 ActNormConstAzt_d_dPhdt*BufRe, ActNormConstAzt_d_dPhdt*BufIm);

				if(PreResidNeeded)
				{
					double One_d_yObsmisNewE2 = One_d_yObsmisNew*One_d_yObsmisNew;
					double dPhdsNew = dPhdt*One_d_yObsmisNewE2;
					double AxNew = Axt*One_d_yObsmisNewE2;
					double AzNew = Azt*One_d_yObsmisNewE2;
					double AxdPhdsNew = AxNew*dPhdsNew, AzdPhdsNew = AzNew*dPhdsNew;

					double Two_d_ymisNewE2 = 2.*One_d_yObsmisNewE2;
					double dAxdsNew = (Btx - xObs_mi_xNew*One_d_yObsmisNew)*Two_d_ymisNewE2;
					double dAzdsNew = (Btz - zObs_mi_zNew*One_d_yObsmisNew)*Two_d_ymisNewE2;
					
					double DiffReX = (dAxdsNew*CosPhNew - AxdPhdsNew*SinPhNew) - (dAxds*CosPh - Ax*dPhds*SinPh);
					double DiffImX = (AxdPhdsNew*CosPhNew + dAxdsNew*SinPhNew) - (Ax*dPhds*CosPh + dAxds*SinPh);
					double DiffReZ = (dAzdsNew*CosPhNew - AzdPhdsNew*SinPhNew) - (dAzds*CosPh - Az*dPhds*SinPh);
					double DiffImZ = (AzdPhdsNew*CosPhNew + dAzdsNew*SinPhNew) - (Az*dPhds*CosPh + dAzds*SinPh);
					srTEFourier DiffEdgeDer(DiffReX, DiffImX, DiffReZ, DiffImZ);
					ComputePreResid(s, sFiNew, Btx, Crdx, IntBtE2x, Btz, Crdz, IntBtE2z, &DiffEdgeDer, &EwRightResid, 1);
				}
**/
					complex<double> ResidX(EwRightResid.EwX_Re, EwRightResid.EwX_Im), ResidZ(EwRightResid.EwZ_Re, EwRightResid.EwZ_Im);
					*ResidValues = ResidX;
					*(ResidValues+1) = ResidZ;
					return 0;
				}
			}
		}

		if(FarFieldExpansDoesNotWork) 
		{
			//return BAD_RAD_INT_TERMINATIONS;

			*ResidValues = Zero; *(ResidValues+1) = Zero;

			CErrWarn::AddWarningMessage(&gVectWarnNos, COMPUTATION_WAS_NOT_DONE_BAD_RAD_INT_TERMINATIONS);
			//OC_hack 15_03_2005
			return 0;
		}

		PreExpX *= mOne_d_IdPhds;
		PreExpZ *= mOne_d_IdPhds;

		*(ResidValues++) = ActNormConst*PreExpX*ExpIPhase;
		*ResidValues = ActNormConst*PreExpZ*ExpIPhase;

		if(ComputeNormalDerivative)
		{
			CosPh = ExpIPhase.real();
			SinPh = ExpIPhase.imag();

			double BufLoc = ActNormConst*One_d_dPhds;
			pEwNormDerResid->EwX_Re = -BufLoc*(t2xNormDer*CosPh + AxNormDer*SinPh);
			pEwNormDerResid->EwX_Im = BufLoc*(AxNormDer*CosPh - t2xNormDer*SinPh);
			pEwNormDerResid->EwZ_Re = -BufLoc*(t2zNormDer*CosPh + AzNormDer*SinPh);
			pEwNormDerResid->EwZ_Im = BufLoc*(AzNormDer*CosPh - t2zNormDer*SinPh);
		}
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		Ax = Btx - xObs;
		Az = Btz - zObs;

		Ph = PIm10e9_d_Lamb*(s*(TrjDatPtr->EbmDat.GammaEm2 + xObs*xObs + zObs*zObs) + IntBtE2x + IntBtE2z - 2.*(xObs*Crdx + zObs*Crdz));
		Ph -= TwoPI*int(Ph/TwoPI);

		double cosPh = cos(Ph), sinPh = sin(Ph);
		complex<double> ExpIPhase(cosPh, sinPh);
		gExpIPhase = ExpIPhase;

		dPhds = PIm10e9_d_Lamb*(TrjDatPtr->EbmDat.GammaEm2 + Ax*Ax + Az*Az);
		complex<double> IdPhds(0., dPhds);
		double One_d_dPhds = 1./dPhds;
		complex<double> mOne_d_IdPhds(0., One_d_dPhds);

		complex<double> PreExpX = Ax;
		complex<double> PreExpZ = Az;

		t1xd = Ax; t1zd = Az;

		if(NumberOfTerms > 1)
		{
			ConBtxBz = (TrjDatPtr->BetaNormConst)*Bz;
			ConBtzBx = (-TrjDatPtr->BetaNormConst)*Bx;
			TwoPIm10e9_d_Lamb = 2.*PIm10e9_d_Lamb;

			d2Phds2 = TwoPIm10e9_d_Lamb*(Ax*ConBtxBz + Az*ConBtzBx);
			dAxds = ConBtxBz;
			dAzds = ConBtzBx;

			double d2Phds2_d_dPhds = d2Phds2*One_d_dPhds;

			t2xd = (dAxds - Ax*d2Phds2_d_dPhds)*One_d_dPhds; 
			t2zd = (dAzds - Az*d2Phds2_d_dPhds)*One_d_dPhds;

			d2Phds2_d_dPhdsE2_ForExpansionTest = d2Phds2_d_dPhds*One_d_dPhds;

			complex<double> t2x(0., t2xd);
			PreExpX += t2x;
			complex<double> t2z(0., t2zd);
			PreExpZ += t2z;
		}

		char FarFieldExpansDoesNotWork = (::fabs(d2Phds2_d_dPhdsE2_ForExpansionTest) > CritFarFieldTermRatio);
		if(FarFieldExpansDoesNotWork) 
		{
			//return BAD_RAD_INT_TERMINATIONS;

			//complex<double> Zero(0., 0.);
			*ResidValues = Zero; *(ResidValues+1) = Zero;

			CErrWarn::AddWarningMessage(&gVectWarnNos, COMPUTATION_WAS_NOT_DONE_BAD_RAD_INT_TERMINATIONS);
			//OC_hack 15_03_2005
			return 0;
		}

		PreExpX *= mOne_d_IdPhds;
		PreExpZ *= mOne_d_IdPhds;

		*(ResidValues++) = ActNormConst*PreExpX*ExpIPhase;
		*ResidValues = ActNormConst*PreExpZ*ExpIPhase;
	}

	if(NumberOfTerms > 1)
	{
		if((ComputeDer == 1) || (ComputeDer == 2))
		{
			complex<double>* pDer = (ComputeDer == 1)? InitDerMan : ((ComputeDer == 2)? FinDerMan : 0);
			complex<double> pex(dAxds, Ax*dPhds);
			*(pDer++) = pex*gExpIPhase;
			complex<double> pez(dAzds, Az*dPhds);
			*pDer = pez*gExpIPhase;

			if(ComputeNormalDerivative)
			{
				double AxNormDer_mu_dPhds = AxNormDer*dPhds;
				double AzNormDer_mu_dPhds = AzNormDer*dPhds;
				if(ComputeDer == 1)
				{
					InitDerXReND = dAxdsNormDer*CosPh - AxNormDer_mu_dPhds*SinPh;
					InitDerXImND = dAxdsNormDer*SinPh + AxNormDer_mu_dPhds*CosPh;
					InitDerZReND = dAzdsNormDer*CosPh - AzNormDer_mu_dPhds*SinPh;
					InitDerZImND = dAzdsNormDer*SinPh + AzNormDer_mu_dPhds*CosPh;
				}
				else if(ComputeDer == 2)
				{
					FinDerXReND = dAxdsNormDer*CosPh - AxNormDer_mu_dPhds*SinPh;
					FinDerXImND = dAxdsNormDer*SinPh + AxNormDer_mu_dPhds*CosPh;
					FinDerZReND = dAzdsNormDer*CosPh - AzNormDer_mu_dPhds*SinPh;
					FinDerZImND = dAzdsNormDer*SinPh + AzNormDer_mu_dPhds*CosPh;
				}
			}
		}
	}
	return 0;
}

//*************************************************************************

void srTRadInt::DetermineIntegIntervalsForRightResidual(double sStGen, const int AmOfParts, double* EdgePoints)
{
	//double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239854E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239842E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	double sEndGen = ObsCoor.y + 70.5*WavelengthIn_m;

	double AmOfWavelengths = 15.;

	EdgePoints[0] = sStGen;
	EdgePoints[AmOfParts] = sEndGen;

	if(sEndGen - sStGen < 8*AmOfWavelengths*WavelengthIn_m)
	{
		double dsAux = (sEndGen - sStGen)/double(AmOfParts);
		for(int i=1; i<AmOfParts; i++) EdgePoints[i] = EdgePoints[i - 1] + dsAux;
	}
	else
	{
		int AmOfPartsLeft = AmOfParts - 1;
		EdgePoints[AmOfPartsLeft] = sEndGen - AmOfWavelengths*WavelengthIn_m;
		double dsAux = 0.;
		if(AmOfPartsLeft > 1)
		{
			EdgePoints[AmOfPartsLeft - 1] = EdgePoints[AmOfPartsLeft] - AmOfWavelengths*WavelengthIn_m;
			AmOfPartsLeft--;
		}
		if(AmOfPartsLeft > 1)
		{
			EdgePoints[AmOfPartsLeft - 1] = EdgePoints[AmOfPartsLeft] - AmOfWavelengths*WavelengthIn_m;
			AmOfPartsLeft--;
		}
		if(AmOfPartsLeft > 1)
		{
			EdgePoints[AmOfPartsLeft - 1] = EdgePoints[AmOfPartsLeft] - AmOfWavelengths*WavelengthIn_m;
			AmOfPartsLeft--;
		}
		if(AmOfPartsLeft > 1)
		{
			EdgePoints[AmOfPartsLeft - 1] = EdgePoints[AmOfPartsLeft] - AmOfWavelengths*WavelengthIn_m;
			AmOfPartsLeft--;
		}
		if(AmOfPartsLeft > 1)
		{
			EdgePoints[AmOfPartsLeft - 1] = EdgePoints[AmOfPartsLeft] - AmOfWavelengths*WavelengthIn_m;
			AmOfPartsLeft--;
		}
		if(AmOfPartsLeft > 1)
		{
			dsAux = (EdgePoints[AmOfPartsLeft] - EdgePoints[0])/double(AmOfPartsLeft);
			for(int i=1; i<AmOfPartsLeft; i++) EdgePoints[i] = EdgePoints[i - 1] + dsAux;
		}
	}
}

//*************************************************************************

srTEFourier srTRadInt::ComputeRightPartOfRadIntegralInNearField(double sStGen, double Btx, double xSt, double IntBtE2xSt, double Btz, double zSt, double IntBtE2zSt)
{
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;
	//double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239854E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	//double WavelengthIn_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239842E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	//double sEndGen = ObsCoor.y - 0.1*WavelengthIn_m;

// Steerable parameters
	double Loc_sIntegRelPrec = (sIntegMethod < 10)? 0.001 : 0.1*sIntegRelPrec;
	const int AmOfParts = 15;
// End Steerable parameters

	double EdgePoints[AmOfParts + 1];
	DetermineIntegIntervalsForRightResidual(sStGen, AmOfParts, EdgePoints);

	const double wfe = 1./3.;
	const double wf1 = 4./3.;
	const double wf2 = 2./3.;

	//double sInterv = (sEndGen - sStGen)/double(AmOfParts);
  	//double sStart = sStGen;

	double sStart = EdgePoints[0], sEnd = EdgePoints[1];

	double IntBtE2xzSt = IntBtE2xSt + IntBtE2zSt;
	double Ph, PhPrev, PhInit;

	srTEFourier EwZero(0,0,0,0);
	srTEFourier OutEwInt(0,0,0,0);

	for(int j=0; j<AmOfParts; j++)
	{
		//int NpOnLevel = 5;
		long long NpOnLevel = 5;
		srTEFourier EwSum1(0,0,0,0), EwSum2(0,0,0,0);
		
		double sStep = (sEnd - sStart)/(NpOnLevel - 1);
		double s = sStart;
		srTEFourier EwSt = RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);
		PhInit = Ph;
		s += sStep;
		//int AmOfLoops = (NpOnLevel - 3) >> 1;
		long long AmOfLoops = (NpOnLevel - 3) >> 1;
		//for(int i=0; i<AmOfLoops; i++)
		for(long long i=0; i<AmOfLoops; i++)
		{
			EwSum1 += RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);
			s += sStep;
			EwSum2 += RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);
			s += sStep;
		}
		EwSum1 += RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);
		s += sStep;
		srTEFourier EwEdge = EwSt + RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);

		double wFxRe = wfe*EwEdge.EwX_Re, wFxIm = wfe*EwEdge.EwX_Im;
		double wFzRe = wfe*EwEdge.EwZ_Re, wFzIm = wfe*EwEdge.EwZ_Im;

		double ActNormConst_sStep = ActNormConst*sStep;

		srTEFourier EwInt;
		EwInt.EwX_Re = OutEwInt.EwX_Re + ActNormConst_sStep*(wFxRe + wf1*EwSum1.EwX_Re + wf2*EwSum2.EwX_Re);
		EwInt.EwX_Im = OutEwInt.EwX_Im + ActNormConst_sStep*(wFxIm + wf1*EwSum1.EwX_Im + wf2*EwSum2.EwX_Im);
		EwInt.EwZ_Re = OutEwInt.EwZ_Re + ActNormConst_sStep*(wFzRe + wf1*EwSum1.EwZ_Re + wf2*EwSum2.EwZ_Re);
		EwInt.EwZ_Im = OutEwInt.EwZ_Im + ActNormConst_sStep*(wFzIm + wf1*EwSum1.EwZ_Im + wf2*EwSum2.EwZ_Im);
		double SqNorm = EwInt.EwX_Re*EwInt.EwX_Re + EwInt.EwX_Im*EwInt.EwX_Im + EwInt.EwZ_Re*EwInt.EwZ_Re + EwInt.EwZ_Im*EwInt.EwZ_Im;

		NpOnLevel--;
		char NotFinishedYet = 1;
		while(NotFinishedYet)
		{
			EwSum2 += EwSum1;
			EwSum1 = EwZero;

			char ThisMayBeTheLastLoop = 1;
			PhPrev = PhInit;
			double DPhMax = 0.;

			double HalfStep = 0.5*sStep;
			s = sStart + HalfStep;
		
			//for(int i=0; i<NpOnLevel; i++)
			for(long long i=0; i<NpOnLevel; i++)
			{
				EwSum1 += RadFuncInDriftSpace(s, sStGen, Btx, xSt, Btz, zSt, IntBtE2xzSt, Ph);
				s += sStep;

				if(Ph - PhPrev > PI) ThisMayBeTheLastLoop = 0;
				double dPh = Ph - PhPrev;

				if(dPh > DPhMax) DPhMax = dPh;
				PhPrev = Ph;
			}
			double ActNormConstHalfStep = ActNormConst*HalfStep;
			srTEFourier LocEwInt;
			LocEwInt.EwX_Re = OutEwInt.EwX_Re + ActNormConstHalfStep*(wFxRe + wf1*EwSum1.EwX_Re + wf2*EwSum2.EwX_Re);
			LocEwInt.EwX_Im = OutEwInt.EwX_Im + ActNormConstHalfStep*(wFxIm + wf1*EwSum1.EwX_Im + wf2*EwSum2.EwX_Im);
			LocEwInt.EwZ_Re = OutEwInt.EwZ_Re + ActNormConstHalfStep*(wFzRe + wf1*EwSum1.EwZ_Re + wf2*EwSum2.EwZ_Re);
			LocEwInt.EwZ_Im = OutEwInt.EwZ_Im + ActNormConstHalfStep*(wFzIm + wf1*EwSum1.EwZ_Im + wf2*EwSum2.EwZ_Im);
			double LocSqNorm = LocEwInt.EwX_Re*LocEwInt.EwX_Re + LocEwInt.EwX_Im*LocEwInt.EwX_Im + LocEwInt.EwZ_Re*LocEwInt.EwZ_Re + LocEwInt.EwZ_Im*LocEwInt.EwZ_Im;
		
			if(ThisMayBeTheLastLoop)
			{
				double TestVal = ::fabs(LocSqNorm - SqNorm);
				NotFinishedYet = (TestVal > Loc_sIntegRelPrec*LocSqNorm);
			}
			EwInt = LocEwInt;
			SqNorm = LocSqNorm;
			sStep = HalfStep; NpOnLevel <<= 1;
		}
		OutEwInt = EwInt;

		sStart = sEnd;
		sEnd = EdgePoints[j + 2];
	}
	return OutEwInt;
}

//*************************************************************************

int srTRadInt::ComputePreResid(double sSt, double sFi, double Btx, double x0, double IntBtE2x0, double Btz, double z0, double IntBtE2z0, srTEFourier* pDiffDer, srTEFourier* pRes, char LeftOrRight) // 0- Left, 1- Right
{
// Steerable parameters
	double Loc_sIntegRelPrec = (sIntegMethod < 10)? 0.0001 : 0.01*sIntegRelPrec;
// End Steerable parameters

	int IntSign = (LeftOrRight == 0)? -1 : 1;
	double Sc = (LeftOrRight == 0)? sFi : sSt;

	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	const double wfe = 0.46666666666667;
	const double wf1 = 1.0666666666667;
	const double wf2 = 0.93333333333333;
	const double wd = 0.066666666666667;
	//int NpOnLevel = 5;
	long long NpOnLevel = 5;

	double sStart = sSt, sEnd = sFi;
	double sStep = (sEnd - sStart)/(NpOnLevel - 1);
	double xObs = ObsCoor.x, yObs = ObsCoor.y, zObs = ObsCoor.z, GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double Ax, Az, Ph, CosPh, SinPh, PhPrev, PhInit;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz;
	double IntBtE20 = IntBtE2x0 + IntBtE2z0;
	double BtE2 = Btx*Btx + Btz*Btz;

	double Sum1XRe=0., Sum1XIm=0., Sum1ZRe=0., Sum1ZIm=0., Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;
	double wFxRe, wFxIm, wFzRe, wFzIm;

	One_d_ymis = 1./(yObs - sStart);
	xObs_mi_x = xObs - X_FreeSpace(sStart, Sc, Btx, x0);
	zObs_mi_z = zObs - X_FreeSpace(sStart, Sc, Btz, z0);
	Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
	Ph = PIm10e9_d_Lamb*(sStart*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(sStart - Sc) + IntBtE20);
	Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
	CosAndSin(Ph, CosPh, SinPh);
	wFxRe = Ax*CosPh; wFxIm = Ax*SinPh; wFzRe = Az*CosPh; wFzIm = Az*SinPh;
	PhInit = Ph;

	double s = sStart + sStep;
	//int AmOfLoops = (NpOnLevel - 3) >> 1;
	long long AmOfLoops = (NpOnLevel - 3) >> 1;

	//for(int i=0; i<AmOfLoops; i++)
	for(long long i=0; i<AmOfLoops; i++)
	{
		One_d_ymis = 1./(yObs - s); 
		xObs_mi_x = xObs - X_FreeSpace(s, Sc, Btx, x0);
		zObs_mi_z = zObs - X_FreeSpace(s, Sc, Btz, z0);
		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
		Ph = PIm10e9_d_Lamb*(s*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(s - Sc) + IntBtE20);
		Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
		CosAndSin(Ph, CosPh, SinPh);
		Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

		One_d_ymis = 1./(yObs - s); 
		xObs_mi_x = xObs - X_FreeSpace(s, Sc, Btx, x0); 
		zObs_mi_z = zObs - X_FreeSpace(s, Sc, Btz, z0);
		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
		Ph = PIm10e9_d_Lamb*(s*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(s - Sc) + IntBtE20);
		Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
		CosAndSin(Ph, CosPh, SinPh);
		Sum2XRe += Ax*CosPh; Sum2XIm += Ax*SinPh; Sum2ZRe += Az*CosPh; Sum2ZIm += Az*SinPh; s += sStep;
	}

	One_d_ymis = 1./(yObs - s); 
	xObs_mi_x = xObs - X_FreeSpace(s, Sc, Btx, x0); 
	zObs_mi_z = zObs - X_FreeSpace(s, Sc, Btz, z0);
	Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
	Ph = PIm10e9_d_Lamb*(s*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(s - Sc) + IntBtE20);
	Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
	CosAndSin(Ph, CosPh, SinPh);
	Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

	One_d_ymis = 1./(yObs - s);
	xObs_mi_x = xObs - X_FreeSpace(s, Sc, Btx, x0); 
	zObs_mi_z = zObs - X_FreeSpace(s, Sc, Btz, z0);
	Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
	Ph = PIm10e9_d_Lamb*(s*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(s - Sc) + IntBtE20);
	Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
	CosAndSin(Ph, CosPh, SinPh);
	wFxRe += Ax*CosPh; wFxIm += Ax*SinPh; wFzRe += Az*CosPh; wFzIm += Az*SinPh;
	wFxRe *= wfe; wFxIm *= wfe; wFzRe *= wfe; wFzIm *= wfe; 

	double wDifDerXRe = wd*pDiffDer->EwX_Re, wDifDerXIm = wd*pDiffDer->EwX_Im;
	double wDifDerZRe = wd*pDiffDer->EwZ_Re, wDifDerZIm = wd*pDiffDer->EwZ_Im;

	double &OutIntXRe = pRes->EwX_Re, &OutIntXIm = pRes->EwX_Im, &OutIntZRe = pRes->EwZ_Re, &OutIntZIm = pRes->EwZ_Im;

	double ActNormConst_sStep = IntSign*ActNormConst*sStep;
	double IntXRe = OutIntXRe + ActNormConst_sStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + sStep*wDifDerXRe);
	double IntXIm = OutIntXIm + ActNormConst_sStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + sStep*wDifDerXIm);
	double IntZRe = OutIntZRe + ActNormConst_sStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + sStep*wDifDerZRe);
	double IntZIm = OutIntZIm + ActNormConst_sStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + sStep*wDifDerZIm);
	double SqNorm = IntXRe*IntXRe + IntXIm*IntXIm + IntZRe*IntZRe + IntZIm*IntZIm;

	NpOnLevel--;
	char NotFinishedYet = 1;
	while(NotFinishedYet)
	{
		Sum2XRe += Sum1XRe; Sum2XIm += Sum1XIm; Sum2ZRe += Sum1ZRe; Sum2ZIm += Sum1ZIm; 
		Sum1XRe = Sum1XIm = Sum1ZRe = Sum1ZIm = 0.;
		char ThisMayBeTheLastLoop = 1;
		PhPrev = PhInit;

		double HalfStep = 0.5*sStep;
		s = sStart + HalfStep;

		double DPhMax = 0.;

		//for(int i=0; i<NpOnLevel; i++)
		for(long long i=0; i<NpOnLevel; i++)
		{
			One_d_ymis = 1./(yObs - s);
			xObs_mi_x = xObs - X_FreeSpace(s, Sc, Btx, x0); 
			zObs_mi_z = zObs - X_FreeSpace(s, Sc, Btz, z0);
			Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*GmEm2 + xObs_mi_x*Nx + zObs_mi_z*Nz + BtE2*(s - Sc) + IntBtE20);
			Ax = (Btx - Nx)*One_d_ymis; Az = (Btz - Nz)*One_d_ymis;
			CosAndSin(Ph, CosPh, SinPh);
			Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; 
			s += sStep;

			if(Ph - PhPrev > FivePIdFour) ThisMayBeTheLastLoop = 0;
			double dPh = Ph - PhPrev;
			if(dPh > DPhMax) DPhMax = dPh;

			PhPrev = Ph;
		}
		double ActNormConstHalfStep = IntSign*ActNormConst*HalfStep;
		double LocIntXRe = OutIntXRe + ActNormConstHalfStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + HalfStep*wDifDerXRe);
		double LocIntXIm = OutIntXIm + ActNormConstHalfStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + HalfStep*wDifDerXIm);
		double LocIntZRe = OutIntZRe + ActNormConstHalfStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + HalfStep*wDifDerZRe);
		double LocIntZIm = OutIntZIm + ActNormConstHalfStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + HalfStep*wDifDerZIm);
		double LocSqNorm = LocIntXRe*LocIntXRe + LocIntXIm*LocIntXIm + LocIntZRe*LocIntZRe + LocIntZIm*LocIntZIm;

		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			NotFinishedYet = (TestVal > Loc_sIntegRelPrec*LocSqNorm);

			//if(ProbablyTheSameLoop && (MaxFluxDensVal > 0.)) NotFinishedYet = (TestVal > CurrentAbsPrec);
			//else NotFinishedYet = (TestVal > sIntegRelPrec*LocSqNorm);
		}
		IntXRe = LocIntXRe; IntXIm = LocIntXIm; IntZRe = LocIntZRe; IntZIm = LocIntZIm; 
		SqNorm = LocSqNorm;
		sStep = HalfStep; NpOnLevel *= 2;
	}

	OutIntXRe = IntXRe; OutIntXIm = IntXIm; OutIntZRe = IntZRe; OutIntZIm = IntZIm;
	//if((ProbablyTheSameLoop && (MaxFluxDensVal < SqNorm)) || !ProbablyTheSameLoop) 
	//{
	//	MaxFluxDensVal = SqNorm; CurrentAbsPrec = sIntegRelPrec*MaxFluxDensVal; 
	//	ProbablyTheSameLoop = 1;
	//}

	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationManualSlower(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;

	complex<double> Zero(0.,0.);
	complex<double> Sum1[2], Sum2[2], F1[2], F2[2], Fb[2], Fe[2];
	Sum1[0] = Sum1[1] = Sum2[0] = Sum2[1] = Zero; 

	double s = sIntegStart;
	FunForRadInt(s, Fb); 

	//int AmOfLoops = (AmOfPointsForManIntegr - 3) >> 1; // AmOfPointsForManIntegr assumed non-even
	long long AmOfLoops = (AmOfPointsForManIntegr - 3) >> 1; // AmOfPointsForManIntegr assumed non-even
	//for(int i=1; i<=AmOfLoops; i++)
	for(long long i=1; i<=AmOfLoops; i++)
	{
		s += sIntegStep; 
		FunForRadInt(s, F1);
		*Sum1 += *F1; Sum1[1] += F1[1];

		s += sIntegStep; 
		FunForRadInt(s, F2); 
		*Sum2 += *F2; Sum2[1] += F2[1];
	}
	s += sIntegStep; 
	FunForRadInt(s, F1); 
	*Sum1 += *F1; Sum1[1] += F1[1];
	FunForRadInt(sIntegFin, Fe);

	if(sIntegMethod == 0)
	{
		double ActNormConst_sIntegStep_d_3 = ActNormConst*sIntegStep*0.333333333333;
		*Sum1 = ActNormConst_sIntegStep_d_3*(*Fb + *Fe + 4.*(*Sum1) + 2.*(*Sum2));
		*(Sum1+1) = ActNormConst_sIntegStep_d_3*(Fb[1] + Fe[1] + 4.*Sum1[1] + 2.*Sum2[1]);
	}
	else if(sIntegMethod == 1)
	{
		double ActNormConst_sIntegStep_d_15 = ActNormConst*sIntegStep*0.0666666666667;
		*Sum1 = ActNormConst_sIntegStep_d_15*(7.*(*Fb + *Fe) + 16.*(*Sum1) + 14.*(*Sum2) + sIntegStep*(*InitDerMan - *FinDerMan));
		*(Sum1+1) = ActNormConst_sIntegStep_d_15*(7.*(Fb[1] + Fe[1]) + 16.*Sum1[1] + 14.*Sum2[1] + sIntegStep*(InitDerMan[1] - FinDerMan[1]));
	}
	OutIntXRe += Sum1->real(); OutIntXIm += Sum1->imag();
	OutIntZRe += (Sum1+1)->real(); OutIntZIm += (Sum1+1)->imag();
	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationManualFaster0(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;

	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	const double wf[] = {1./3., 4./3., 2./3.};

	double SumXRe=0., SumXIm=0., SumZRe=0., SumZIm=0.;
	double FxRe, FxIm, FzRe, FzIm;
	double sArg = sIntegStart;
	double *pBtx=BtxArr, *pBtz=BtzArr, *pX=XArr, *pZ=ZArr, *pIntBtxE2=IntBtxE2Arr, *pIntBtzE2=IntBtzE2Arr;

	double One_d_ymis, xObs_mi_x, zObs_mi_z, Phase, Ax, Az, Nx, Nz, W;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	double GmEm2 = TrjDatPtr->EbmDat.GammaEm2;

	double SumXReND=0., SumXImND=0., SumZReND=0., SumZImND=0.;
	double FxReND, FxImND, FzReND, FzImND;
	double Ny, L_x_N;

	char CountTo3 = 0;
	//int AmOfPointsForManIntegr_mi_1 = AmOfPointsForManIntegr - 1;
	long long AmOfPointsForManIntegr_mi_1 = AmOfPointsForManIntegr - 1;

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		//for(int i=0; i<AmOfPointsForManIntegr; i++)
		for(long long i=0; i<AmOfPointsForManIntegr; i++)
		{
			if(i==AmOfPointsForManIntegr_mi_1) CountTo3 = 0;
			if(CountTo3==3) CountTo3 = 1;

			One_d_ymis = 1./(ObsCoor.y - sArg);
			xObs_mi_x = xObs - *(pX++);
			zObs_mi_z = zObs - *(pZ++);

			double LongTerm = *(pIntBtxE2++) + *(pIntBtzE2++);
			//double a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis;
			//double a = a0*One_d_ymis;
			//Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));

			//OC_test
				double a0 = LongTerm + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis;
                Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + a0);
			//end OC_test
			
			Phase -= TwoPI*int(Phase/TwoPI);
			double CosPhase = cos(Phase), SinPhase = sin(Phase);

			Nx = xObs_mi_x*One_d_ymis; 
			Nz = zObs_mi_z*One_d_ymis;

			Ax = (*(pBtx++) - Nx)*One_d_ymis;
			FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
			Az = (*(pBtz++) - Nz)*One_d_ymis;
			FzRe = Az*CosPhase; FzIm = Az*SinPhase;

			W = (*(wf+(CountTo3++)));
			SumXRe += W*FxRe; SumXIm += W*FxIm; 
			SumZRe += W*FzRe; SumZIm += W*FzIm; 

			if(ComputeNormalDerivative)
			{
				Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
				L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
				FxReND = FxRe*L_x_N; FxImND = FxIm*L_x_N;
				FzReND = FzRe*L_x_N; FzImND = FzIm*L_x_N;

				SumXReND += W*FxReND; SumXImND += W*FxImND;
				SumZReND += W*FzReND; SumZImND += W*FzImND; 
			}

			sArg += sIntegStep;
		}
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		double AngPhConst = GmEm2 + xObs*xObs + zObs*zObs;
		//for(int i=0; i<AmOfPointsForManIntegr; i++)
		for(long long i=0; i<AmOfPointsForManIntegr; i++)
		{
			if(CountTo3==3) CountTo3 = 1;
			if(i==AmOfPointsForManIntegr_mi_1) CountTo3 = 0;

			double Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - 2.*(xObs*(*(pX++)) + zObs*(*(pZ++))));
			Phase -= TwoPI*int(Phase/TwoPI);
			double CosPhase = cos(Phase), SinPhase = sin(Phase);
			Ax = *(pBtx++) - xObs;
			FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
			Az = *(pBtz++) - zObs;
			FzRe = Az*CosPhase; FzIm = Az*SinPhase;

			W = (*(wf+(CountTo3++)));
			SumXRe += W*FxRe; SumXIm += W*FxIm; 
			SumZRe += W*FzRe; SumZIm += W*FzIm; 
			sArg += sIntegStep;
		}
	}
	double ActNormConstIntegStep = sIntegStep*ActNormConst;
	OutIntXRe += ActNormConstIntegStep*SumXRe;
	OutIntXIm += ActNormConstIntegStep*SumXIm;
	OutIntZRe += ActNormConstIntegStep*SumZRe;
	OutIntZIm += ActNormConstIntegStep*SumZIm;

	if(ComputeNormalDerivative)
	{
		pEwNormDer->EwX_Re += ActNormConstIntegStep*SumXReND;
		pEwNormDer->EwX_Im += ActNormConstIntegStep*SumXImND;
		pEwNormDer->EwZ_Re += ActNormConstIntegStep*SumZReND;
		pEwNormDer->EwZ_Im += ActNormConstIntegStep*SumZImND;
	}

	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationManualFaster1(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;

	double xObs = ObsCoor.x, zObs = ObsCoor.z, yObs = ObsCoor.y;
	double Sum1XRe=0., Sum1XIm=0., Sum1ZRe=0., Sum1ZIm=0., Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;
	double FxRe, FxIm, FzRe, FzIm;

	double sArg = sIntegStart;
	double *pBtx=BtxArr, *pBtz=BtzArr, *pX=XArr, *pZ=ZArr, *pIntBtxE2=IntBtxE2Arr, *pIntBtzE2=IntBtzE2Arr;

	double One_d_ymis, xObs_mi_x, zObs_mi_z, Phase, Ax, Az, CosPhase, SinPhase, Nx, Nz;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double GmEm2 = TrjDatPtr->EbmDat.GammaEm2;

	double Sum1XReND=0., Sum1XImND=0., Sum1ZReND=0., Sum1ZImND=0., Sum2XReND=0., Sum2XImND=0., Sum2ZReND=0., Sum2ZImND=0.;
	double FxReND, FxImND, FzReND, FzImND;
	double Ny, L_x_N, AxND, AzND;

	//char CountTo3 = 0;
	//int AmOfPointsForManIntegr_mi_1 = AmOfPointsForManIntegr - 1;

	//int AmOfLoops = (AmOfPointsForManIntegr - 3) >> 1;
	long long AmOfLoops = (AmOfPointsForManIntegr - 3) >> 1;

	const double One_d_15 = 1./15.;
	const double Seven_d_15 = 7.*One_d_15;
	double sIntegStep_d_15 = sIntegStep*One_d_15;
	double sIntegStep_d_15_mu_7 = sIntegStep*Seven_d_15;

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		One_d_ymis = 1./(yObs - sArg);
		xObs_mi_x = xObs - *(pX++);
		zObs_mi_z = zObs - *(pZ++);
		Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		CosAndSin(Phase, CosPhase, SinPhase);

		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		Ax = (*(pBtx++) - Nx)*One_d_ymis;
		FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
		Az = (*(pBtz++) - Nz)*One_d_ymis;
		FzRe = Az*CosPhase; FzIm = Az*SinPhase;

		if(ComputeNormalDerivative)
		{
			Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
			L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
			FxReND = FxRe*L_x_N; FxImND = FxIm*L_x_N;
			FzReND = FzRe*L_x_N; FzImND = FzIm*L_x_N;
		}

		sArg += sIntegStep;

		//for(int i=1; i<=AmOfLoops; i++)
		for(long long i=1; i<=AmOfLoops; i++)
		{
			One_d_ymis = 1./(yObs - sArg);
			xObs_mi_x = xObs - *(pX++);
			zObs_mi_z = zObs - *(pZ++);
			Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
			CosAndSin(Phase, CosPhase, SinPhase);

			Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

			Ax = (*(pBtx++) - Nx)*One_d_ymis;
			Sum1XRe += Ax*CosPhase; Sum1XIm += Ax*SinPhase;
			Az = (*(pBtz++) - Nz)*One_d_ymis;
			Sum1ZRe += Az*CosPhase; Sum1ZIm += Az*SinPhase;

			if(ComputeNormalDerivative)
			{
				Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
				L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
				AxND = Ax*L_x_N; AzND = Az*L_x_N;
				Sum1XReND += AxND*CosPhase; Sum1XImND += AxND*SinPhase;
				Sum1ZReND += AzND*CosPhase; Sum1ZImND += AzND*SinPhase;
			}

			sArg += sIntegStep;

			One_d_ymis = 1./(yObs - sArg);
			xObs_mi_x = xObs - *(pX++);
			zObs_mi_z = zObs - *(pZ++);
			Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
			CosAndSin(Phase, CosPhase, SinPhase);

			Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

			Ax = (*(pBtx++) - Nx)*One_d_ymis;
			Sum2XRe += Ax*CosPhase; Sum2XIm += Ax*SinPhase;
			Az = (*(pBtz++) - Nz)*One_d_ymis;
			Sum2ZRe += Az*CosPhase; Sum2ZIm += Az*SinPhase;

			if(ComputeNormalDerivative)
			{
				Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
				L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
				AxND = Ax*L_x_N; AzND = Az*L_x_N;
				Sum2XReND += AxND*CosPhase; Sum2XImND += AxND*SinPhase;
				Sum2ZReND += AzND*CosPhase; Sum2ZImND += AzND*SinPhase;
			}

			sArg += sIntegStep;
		}
		One_d_ymis = 1./(yObs - sArg);
		xObs_mi_x = xObs - *(pX++);
		zObs_mi_z = zObs - *(pZ++);
		Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		CosAndSin(Phase, CosPhase, SinPhase);

		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		Ax = (*(pBtx++) - Nx)*One_d_ymis;
		Sum1XRe += Ax*CosPhase; Sum1XIm += Ax*SinPhase;
		Az = (*(pBtz++) - Nz)*One_d_ymis;
		Sum1ZRe += Az*CosPhase; Sum1ZIm += Az*SinPhase;

		if(ComputeNormalDerivative)
		{
			Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
			L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
			AxND = Ax*L_x_N; AzND = Az*L_x_N;
			Sum1XReND += AxND*CosPhase; Sum1XImND += AxND*SinPhase;
			Sum1ZReND += AzND*CosPhase; Sum1ZImND += AzND*SinPhase;
		}

		sArg = sIntegFin;

		One_d_ymis = 1./(yObs - sArg);
		xObs_mi_x = xObs - *(pX++);
		zObs_mi_z = zObs - *(pZ++);
		Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
		CosAndSin(Phase, CosPhase, SinPhase);

		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		Ax = (*(pBtx++) - Nx)*One_d_ymis;
		FxRe += Ax*CosPhase; FxIm += Ax*SinPhase;
		Az = (*(pBtz++) - Nz)*One_d_ymis;
		FzRe += Az*CosPhase; FzIm += Az*SinPhase;

		if(ComputeNormalDerivative)
		{
			Ny = 1. - 0.5*(Nx*Nx + Nz*Nz);
			L_x_N = SurfNorm.x*Nx + SurfNorm.y*Ny + SurfNorm.z*Nz;
			AxND = Ax*L_x_N; AzND = Az*L_x_N;

			FxReND += AxND*CosPhase; FxImND += AxND*SinPhase;
			FzReND += AzND*CosPhase; FzImND += AzND*SinPhase;
		}
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		double AngPhConst = GmEm2 + xObs*xObs + zObs*zObs;
		double Two_xObs = 2.*xObs, Two_zObs = 2.*zObs;

		Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		CosAndSin(Phase, CosPhase, SinPhase);

		Ax = *(pBtx++) - xObs; FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
		Az = *(pBtz++) - zObs; FzRe = Az*CosPhase; FzIm = Az*SinPhase;
		sArg += sIntegStep;

		//for(int i=1; i<=AmOfLoops; i++)
		for(long long i=1; i<=AmOfLoops; i++)
		{
			Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			CosAndSin(Phase, CosPhase, SinPhase);

			Ax = *(pBtx++) - xObs; Sum1XRe += Ax*CosPhase; Sum1XIm += Ax*SinPhase;
			Az = *(pBtz++) - zObs; Sum1ZRe += Az*CosPhase; Sum1ZIm += Az*SinPhase;
			sArg += sIntegStep;

			Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			CosAndSin(Phase, CosPhase, SinPhase);

			Ax = *(pBtx++) - xObs; Sum2XRe += Ax*CosPhase; Sum2XIm += Ax*SinPhase;
			Az = *(pBtz++) - zObs; Sum2ZRe += Az*CosPhase; Sum2ZIm += Az*SinPhase;
			sArg += sIntegStep;
		}
		Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		CosAndSin(Phase, CosPhase, SinPhase);

		Ax = *(pBtx++) - xObs; Sum1XRe += Ax*CosPhase; Sum1XIm += Ax*SinPhase;
		Az = *(pBtz++) - zObs; Sum1ZRe += Az*CosPhase; Sum1ZIm += Az*SinPhase;
		sArg = sIntegFin;

		Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		CosAndSin(Phase, CosPhase, SinPhase);

		Ax = *(pBtx++) - xObs; FxRe += Ax*CosPhase; FxIm += Ax*SinPhase;
		Az = *(pBtz++) - zObs; FzRe += Az*CosPhase; FzIm += Az*SinPhase;
	}

	Sum1XRe = sIntegStep_d_15_mu_7*(FxRe + 2.2857142857143*Sum1XRe + 2.*Sum2XRe);
	Sum1XIm = sIntegStep_d_15_mu_7*(FxIm + 2.2857142857143*Sum1XIm + 2.*Sum2XIm);
	Sum1ZRe = sIntegStep_d_15_mu_7*(FzRe + 2.2857142857143*Sum1ZRe + 2.*Sum2ZRe);
	Sum1ZIm = sIntegStep_d_15_mu_7*(FzIm + 2.2857142857143*Sum1ZIm + 2.*Sum2ZIm);

	double wBufDer = sIntegStep*sIntegStep_d_15;

	complex<double> wDifDerX = wBufDer*(*InitDerMan - *FinDerMan);
	OutIntXRe += ActNormConst*(Sum1XRe + wDifDerX.real());
	OutIntXIm += ActNormConst*(Sum1XIm + wDifDerX.imag());
	complex<double> wDifDerZ = wBufDer*(*(InitDerMan+1) - *(FinDerMan+1));
	OutIntZRe += ActNormConst*(Sum1ZRe + wDifDerZ.real());
	OutIntZIm += ActNormConst*(Sum1ZIm + wDifDerZ.imag());

	if(ComputeNormalDerivative)
	{
		Sum1XReND = sIntegStep_d_15_mu_7*(FxReND + 2.2857142857143*Sum1XReND + 2.*Sum2XReND);
		Sum1XImND = sIntegStep_d_15_mu_7*(FxImND + 2.2857142857143*Sum1XImND + 2.*Sum2XImND);
		Sum1ZReND = sIntegStep_d_15_mu_7*(FzReND + 2.2857142857143*Sum1ZReND + 2.*Sum2ZReND);
		Sum1ZImND = sIntegStep_d_15_mu_7*(FzImND + 2.2857142857143*Sum1ZImND + 2.*Sum2ZImND);

		pEwNormDer->EwX_Re += ActNormConst*(Sum1XReND + wBufDer*(InitDerXReND - FinDerXReND));
		pEwNormDer->EwX_Im += ActNormConst*(Sum1XImND + wBufDer*(InitDerXImND - FinDerXImND));
		pEwNormDer->EwZ_Re += ActNormConst*(Sum1ZReND + wBufDer*(InitDerZReND - FinDerZReND));
		pEwNormDer->EwZ_Im += ActNormConst*(Sum1ZImND + wBufDer*(InitDerZImND - FinDerZImND));
	}

	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationManualFaster2(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;

	double xObs = ObsCoor.x, zObs = ObsCoor.z;
	const double wf[] = {3.*31./224., 3.*81./224., 3.*81./224., 6.*31./224.};
	const double wd[] = {3.*19./1120., -3.*27./1120., 3.*27./1120., 0.};

	double SumXRe=0., SumXIm=0., SumZRe=0., SumZIm=0., dSumXRe=0., dSumXIm=0., dSumZRe=0., dSumZIm=0.;
	double FxRe, FxIm, FzRe, FzIm, dFxRe, dFxIm, dFzRe, dFzIm;
	double sArg = sIntegStart;
	double *pBtx=BtxArr, *pBtz=BtzArr, *pX=XArr, *pZ=ZArr, *pIntBtxE2=IntBtxE2Arr, *pIntBtzE2=IntBtzE2Arr, *pBx = BxArr, *pBz = BzArr;

	double One_d_ymis, xObs_mi_x, zObs_mi_z, Phase, Ax, Az, Wf, Wd, dPhds, dAxds, dAzds;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double GmEm2 = TrjDatPtr->EbmDat.GammaEm2;

	char CountTo4 = 0;
	//int AmOfPointsForManIntegr_mi_1 = AmOfPointsForManIntegr - 1;
	long long AmOfPointsForManIntegr_mi_1 = AmOfPointsForManIntegr - 1;

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		//for(int i=0; i<AmOfPointsForManIntegr; i++)
		for(long long i=0; i<AmOfPointsForManIntegr; i++)
		{
			if(i==AmOfPointsForManIntegr_mi_1) CountTo4 = 0;
			if(CountTo4==4) CountTo4 = 1;

			One_d_ymis = 1./(ObsCoor.y - sArg);
			xObs_mi_x = xObs - *(pX++);
			zObs_mi_z = zObs - *(pZ++);
			Phase = PIm10e9_d_Lamb*(sArg*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + (xObs_mi_x*xObs_mi_x + zObs_mi_z*zObs_mi_z)*One_d_ymis);
			Phase -= TwoPI*int(Phase/TwoPI);
			double CosPhase = cos(Phase), SinPhase = sin(Phase);

			double Btx_mi_Nx = *(pBtx++) - xObs_mi_x*One_d_ymis;
			double Btz_mi_Nz = *(pBtz++) - zObs_mi_z*One_d_ymis;
			Ax = Btx_mi_Nx*One_d_ymis;
			FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
			Az = Btz_mi_Nz*One_d_ymis;
			FzRe = Az*CosPhase; FzIm = Az*SinPhase;

			if(CountTo4 < 3)
			{
				dPhds = PIm10e9_d_Lamb*(GmEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
				dAxds = (2.*Ax + (TrjDatPtr->BetaNormConst)*(*pBz))*One_d_ymis;
				dAzds = (2.*Az + (-TrjDatPtr->BetaNormConst)*(*pBx))*One_d_ymis;

				dFxRe = dAxds*CosPhase - Ax*dPhds*SinPhase;
				dFxIm = dAxds*SinPhase + Ax*dPhds*CosPhase;
				dFzRe = dAzds*CosPhase - Az*dPhds*SinPhase;
				dFzIm = dAzds*SinPhase + Az*dPhds*CosPhase;
				Wd = (*(wd+CountTo4));
				if(i==AmOfPointsForManIntegr_mi_1) Wd = -Wd;

				dSumXRe += Wd*dFxRe; dSumXIm += Wd*dFxIm; 
				dSumZRe += Wd*dFzRe; dSumZIm += Wd*dFzIm;
			}
			pBz++; pBx++;

			Wf = (*(wf+(CountTo4++)));
			SumXRe += Wf*FxRe; SumXIm += Wf*FxIm;
			SumZRe += Wf*FzRe; SumZIm += Wf*FzIm;
			sArg += sIntegStep;
		}
	}
	else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
	{
		double AngPhConst = GmEm2 + xObs*xObs + zObs*zObs;
		//for(int i=0; i<AmOfPointsForManIntegr; i++)
		for(long long i=0; i<AmOfPointsForManIntegr; i++)
		{
			if(CountTo4==4) CountTo4 = 1;
			if(i==AmOfPointsForManIntegr_mi_1) CountTo4 = 0;

			double Phase = PIm10e9_d_Lamb*(sArg*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - 2.*(xObs*(*(pX++)) + zObs*(*(pZ++))));
			Phase -= TwoPI*int(Phase/TwoPI);
			double CosPhase = cos(Phase), SinPhase = sin(Phase);
			Ax = *(pBtx++) - xObs;

			FxRe = Ax*CosPhase; FxIm = Ax*SinPhase;
			Az = *(pBtz++) - zObs;

			FzRe = Az*CosPhase; FzIm = Az*SinPhase;

			if(CountTo4 < 3)
			{
				dPhds = PIm10e9_d_Lamb*(GmEm2 + Ax*Ax + Az*Az);
				dAxds = (TrjDatPtr->BetaNormConst)*(*pBz);
				dAzds = (-TrjDatPtr->BetaNormConst)*(*pBx);

				dFxRe = dAxds*CosPhase - Ax*dPhds*SinPhase;
				dFxIm = dAxds*SinPhase + Ax*dPhds*CosPhase;
				dFzRe = dAzds*CosPhase - Az*dPhds*SinPhase;
				dFzIm = dAzds*SinPhase + Az*dPhds*CosPhase;

				Wd = (*(wd+CountTo4));
				if(i==AmOfPointsForManIntegr_mi_1) Wd = -Wd;

				dSumXRe += Wd*dFxRe; dSumXIm += Wd*dFxIm; 
				dSumZRe += Wd*dFzRe; dSumZIm += Wd*dFzIm;
			}
			pBz++; pBx++;

			Wf = (*(wf+(CountTo4++)));
			SumXRe += Wf*FxRe; SumXIm += Wf*FxIm; 
			SumZRe += Wf*FzRe; SumZIm += Wf*FzIm; 
			sArg += sIntegStep;
		}
	}

	double ActNormConstIntegStep = ActNormConst*sIntegStep;
	OutIntXRe += ActNormConstIntegStep*(SumXRe + sIntegStep*dSumXRe);
	OutIntXIm += ActNormConstIntegStep*(SumXIm + sIntegStep*dSumXIm);
	OutIntZRe += ActNormConstIntegStep*(SumZRe + sIntegStep*dSumZRe);
	OutIntZIm += ActNormConstIntegStep*(SumZIm + sIntegStep*dSumZIm);
	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationAuto1(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
	//const long NpOnLevelMaxNoResult = 800000000; //5000000; //2000000; // To steer; to stop computation as unsuccessful
	const long long NpOnLevelMaxNoResult = 800000000; //5000000; //2000000; // To steer; to stop computation as unsuccessful

	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	double sStart = sIntegStart;
	double sEnd = sIntegFin;

	char NearField = (DistrInfoDat.CoordOrAngPresentation == CoordPres);

	//long NpOnLevel = 5; // Must be non-even!
	long long NpOnLevel = 5; // Must be non-even!

	int result;
	if(NumberOfLevelsFilled == 0) if(result = FillNextLevel(0, sStart, sEnd, NpOnLevel)) return result;

	double sStep = (sEnd - sStart)/(NpOnLevel - 1);
	double Ax, Az, Ph, CosPh, SinPh, PhPrev, PhInit;
	double LongTerm, a0; //, a;

	double xObs = ObsCoor.x, yObs = ObsCoor.y, zObs = ObsCoor.z, GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz;
	double AngPhConst, Two_xObs, Two_zObs;

	double Sum1XRe=0., Sum1XIm=0., Sum1ZRe=0., Sum1ZIm=0., Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;
	double wFxRe, wFxIm, wFzRe, wFzIm;
	int LevelNo = 0, IndxOnLevel = 0;

	double *pBtx = *BtxArrP, *pBtz = *BtzArrP, *pX = *XArrP, *pZ = *ZArrP, *pIntBtxE2 = *IntBtxE2ArrP, *pIntBtzE2 = *IntBtzE2ArrP;
	double BtxLoc, xLoc, IntBtxE2Loc, BtzLoc, zLoc, IntBtzE2Loc;

	if(NearField) 
	{
		One_d_ymis = 1./(yObs - sStart); 
		xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		LongTerm = *(pIntBtxE2++) + *(pIntBtzE2++);
		//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
		//a = a0*One_d_ymis;
		//Ph = PIm10e9_d_Lamb*(sStart*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
		//OC_test
			a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
            Ph = PIm10e9_d_Lamb*(sStart*GmEm2 + a0);
		//end OC_test

		Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
	}
	else
	{
		AngPhConst = GmEm2 + xObs*xObs + zObs*zObs; Two_xObs = 2.*xObs; Two_zObs = 2.*zObs;
		Ph = PIm10e9_d_Lamb*(sStart*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
	}
	CosAndSin(Ph, CosPh, SinPh);
	wFxRe = Ax*CosPh; wFxIm = Ax*SinPh; wFzRe = Az*CosPh; wFzIm = Az*SinPh;
	PhInit = Ph;

	double s = sStart + sStep;
	IndxOnLevel = 1;
	//int AmOfLoops = (NpOnLevel - 3) >> 1;
	long long AmOfLoops = (NpOnLevel - 3) >> 1;

	//for(int i=0; i<AmOfLoops; i++)
	for(long long i=0; i<AmOfLoops; i++)
	{
		if(NearField)
		{
			One_d_ymis = 1./(yObs - s); 
			xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
			Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
			
			LongTerm = *(pIntBtxE2++) + *(pIntBtzE2++);
			//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
			//a = a0*One_d_ymis;
			//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
			//OC_test
				a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
                Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
			//end OC_test

			Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
		}
		CosAndSin(Ph, CosPh, SinPh);
		Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

		if(NearField) 
		{
			One_d_ymis = 1./(yObs - s); 
			xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
			Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

			LongTerm = *(pIntBtxE2++) + *(pIntBtzE2++);
			//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
			//a = a0*One_d_ymis;
			//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
			//OC_test
				a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
                Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
			//end OC_test

			Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
		}
		CosAndSin(Ph, CosPh, SinPh);
		Sum2XRe += Ax*CosPh; Sum2XIm += Ax*SinPh; Sum2ZRe += Az*CosPh; Sum2ZIm += Az*SinPh; s += sStep;
	}
	if(NearField) 
	{
		One_d_ymis = 1./(yObs - s); 
		xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		LongTerm = *(pIntBtxE2++) + *(pIntBtzE2++);
		//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
		//a = a0*One_d_ymis;
		//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
		//OC_test
			a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
            Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
		//end OC_test

		Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
	}
	else
	{
		Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
	}
	CosAndSin(Ph, CosPh, SinPh);
	Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

	if(NearField)
	{
		One_d_ymis = 1./(yObs - s); 
		xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
		Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;

		LongTerm = *pIntBtxE2 + *pIntBtzE2;
		//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
		//a = a0*One_d_ymis;
		//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
		//OC_test
			a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
            Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
		//end OC_test

		Ax = (*pBtx - Nx)*One_d_ymis; Az = (*pBtz - Nz)*One_d_ymis;
	}
	else
	{
		Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_xObs*(*pX) + Two_zObs*(*pZ)));
		Ax = *pBtx - xObs; Az = *pBtz - zObs;
	}
	CosAndSin(Ph, CosPh, SinPh);
	wFxRe += Ax*CosPh; wFxIm += Ax*SinPh; wFzRe += Az*CosPh; wFzIm += Az*SinPh;
	wFxRe *= wfe; wFxIm *= wfe; wFzRe *= wfe; wFzIm *= wfe; 

	complex<double> DifDerX = *InitDerMan - *FinDerMan;
	double wDifDerXRe = wd*DifDerX.real(), wDifDerXIm = wd*DifDerX.imag();
	complex<double> DifDerZ = *(InitDerMan+1) - *(FinDerMan+1);
	double wDifDerZRe = wd*DifDerZ.real(), wDifDerZIm = wd*DifDerZ.imag();

	double ActNormConst_sStep = ActNormConst*sStep;
	double IntXRe = OutIntXRe + ActNormConst_sStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + sStep*wDifDerXRe);
	double IntXIm = OutIntXIm + ActNormConst_sStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + sStep*wDifDerXIm);
	double IntZRe = OutIntZRe + ActNormConst_sStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + sStep*wDifDerZRe);
	double IntZIm = OutIntZIm + ActNormConst_sStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + sStep*wDifDerZIm);
	double SqNorm = IntXRe*IntXRe + IntXIm*IntXIm + IntZRe*IntZRe + IntZIm*IntZIm;

	//char ExtraPassForAnyCase = 0; 

	NpOnLevel--;
	char NotFinishedYet = 1;
	while(NotFinishedYet)
	{
		Sum2XRe += Sum1XRe; Sum2XIm += Sum1XIm; Sum2ZRe += Sum1ZRe; Sum2ZIm += Sum1ZIm; 
		Sum1XRe = Sum1XIm = Sum1ZRe = Sum1ZIm = 0.;
		char ThisMayBeTheLastLoop = 1;
		PhPrev = PhInit;
		LevelNo++;

		double HalfStep = 0.5*sStep;
		s = sStart + HalfStep;

		if(LevelNo <= MaxLevelForMeth_10_11)
		{
			if(NumberOfLevelsFilled <= LevelNo) if(result = FillNextLevel(LevelNo, s, sEnd - HalfStep, NpOnLevel)) return result;
			pBtx = BtxArrP[LevelNo]; pBtz = BtzArrP[LevelNo]; pX = XArrP[LevelNo]; pZ = ZArrP[LevelNo]; pIntBtxE2 = IntBtxE2ArrP[LevelNo]; pIntBtzE2 = IntBtzE2ArrP[LevelNo];
		}

			double DPhMax = 0.;

		//for(long i=0; i<NpOnLevel; i++)
		for(long long i=0; i<NpOnLevel; i++)
		{
			if(LevelNo > MaxLevelForMeth_10_11)
			{
				pBtx = &BtxLoc; pX = &xLoc; pIntBtxE2 = &IntBtxE2Loc;
				pBtz = &BtzLoc; pZ = &zLoc; pIntBtzE2 = &IntBtzE2Loc;
				TrjDatPtr->CompTrjDataDerivedAtPoint(s, *pBtx, *pX, *pIntBtxE2, *pBtz, *pZ, *pIntBtzE2);
			}

			if(NearField)
			{
				One_d_ymis = 1./(yObs - s);
				xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
				Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;

				LongTerm = *pIntBtxE2 + *pIntBtzE2;
				//a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
				//a = a0*One_d_ymis;
				//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
				//OC_test
					a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
                    Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
				//end OC_test

				Ax = (*pBtx - Nx)*One_d_ymis; Az = (*pBtz - Nz)*One_d_ymis;
			}
			else
			{
				Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_xObs*(*pX) + Two_zObs*(*pZ)));
				Ax = *pBtx - xObs; Az = *pBtz - zObs;
			}

			if(LevelNo <= MaxLevelForMeth_10_11)
			{
				pBtx++; pX++; pIntBtxE2++; pBtz++; pZ++; pIntBtzE2++;
			}

			CosAndSin(Ph, CosPh, SinPh);
			Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; 

			//DEBUG
			//if(::fabs(Ax) > 0.1)
			//{
			//	int aha = 1;
			//}
			//END DEBUG

			s += sStep;

			if(Ph - PhPrev > PI) ThisMayBeTheLastLoop = 0;

				double dPh = Ph - PhPrev;
				if(dPh > DPhMax) DPhMax = dPh;

			PhPrev = Ph;

				//if(i > NpOnLevel - 30)
				//{
				//	//int aha = 1;
				//	Ph *= 1.;
				//}
		}
		double ActNormConstHalfStep = ActNormConst*HalfStep;
		double LocIntXRe = OutIntXRe + ActNormConstHalfStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + HalfStep*wDifDerXRe);
		double LocIntXIm = OutIntXIm + ActNormConstHalfStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + HalfStep*wDifDerXIm);
		double LocIntZRe = OutIntZRe + ActNormConstHalfStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + HalfStep*wDifDerZRe);
		double LocIntZIm = OutIntZIm + ActNormConstHalfStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + HalfStep*wDifDerZIm);
		double LocSqNorm = LocIntXRe*LocIntXRe + LocIntXIm*LocIntXIm + LocIntZRe*LocIntZRe + LocIntZIm*LocIntZIm;

		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			//char SharplyGoesDown = (LocSqNorm < 0.1*SqNorm);

			char NotFinishedYetFirstTest;
			if(ProbablyTheSameLoop && (MaxFluxDensVal > 0.)) NotFinishedYetFirstTest = (TestVal > CurrentAbsPrec);
			else NotFinishedYetFirstTest = (TestVal > sIntegRelPrec*LocSqNorm);

			if(!NotFinishedYetFirstTest)
			{
				NotFinishedYet = 0;
			}
		}

		if(NotFinishedYet)
		{
			if(NpOnLevel > NpOnLevelMaxNoResult) 
			{
				return CAN_NOT_COMPUTE_RADIATION_INTEGRAL;
			}
		}

		IntXRe = LocIntXRe; IntXIm = LocIntXIm; IntZRe = LocIntZRe; IntZIm = LocIntZIm; 
		SqNorm = LocSqNorm;
		sStep = HalfStep; NpOnLevel *= 2;
	}

	OutIntXRe = IntXRe; OutIntXIm = IntXIm; OutIntZRe = IntZRe; OutIntZIm = IntZIm; 
	if((ProbablyTheSameLoop && (MaxFluxDensVal < SqNorm)) || !ProbablyTheSameLoop) 
	{
		MaxFluxDensVal = SqNorm; CurrentAbsPrec = sIntegRelPrec*MaxFluxDensVal; 
		ProbablyTheSameLoop = 1;
	}
	return 0;
}

//*************************************************************************

//int srTRadInt::RadIntegrationAuto1M(double sStart, double sEnd, double* FunArr, double* EdgeDerArr, int AmOfInitPo, int ThisLevNo, double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm)
int srTRadInt::RadIntegrationAuto1M(double sStart, double sEnd, double* FunArr, double* EdgeDerArr, long long AmOfInitPo, int ThisLevNo, double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm)
{
	//const long NpOnLevelMaxNoResult = 800000000; //5000000; // To steer; to stop computation as unsuccessful
	const long long NpOnLevelMaxNoResult = 800000000; //5000000; // To steer; to stop computation as unsuccessful
	//const int NpOnZeroLev = 5;
	//const int NpOnZeroLev_mi_1 = NpOnZeroLev - 1;
	//double GenInterv = sIntegFin - sIntegStart;
	//double NpOnZeroLevmi1_d_GenInterv = NpOnZeroLev_mi_1/GenInterv;

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	//int AmOfInitPo_mi_1 = AmOfInitPo - 1;
	//int AmOfInitPoMi1_mu_4 = AmOfInitPo_mi_1*4;
	long long AmOfInitPo_mi_1 = AmOfInitPo - 1;
	long long AmOfInitPoMi1_mu_4 = AmOfInitPo_mi_1*4;

	double *t1 = FunArr, *t2 = FunArr + AmOfInitPoMi1_mu_4;
	double wFxRe = wfe*(*(t1++) + *(t2++));
	double wFxIm = wfe*(*(t1++) + *(t2++));
	double wFzRe = wfe*(*(t1++) + *(t2++));
	double wFzIm = wfe*(*t1 + *(t2++));

	double PhInit = *t2, PhPrev;

	t1 = EdgeDerArr; t2 = EdgeDerArr + 4;
	double wdFxRe = wd*(*(t1++) - *(t2++));
	double wdFxIm = wd*(*(t1++) - *(t2++));
	double wdFzRe = wd*(*(t1++) - *(t2++));
	double wdFzIm = wd*(*t1 - *t2);

	double Sum1XRe=0., Sum1XIm=0., Sum1ZRe=0., Sum1ZIm=0., Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;

	double *t = FunArr + 4; 
	if(AmOfInitPo > 3)
	{
		//int AmOfLoops = (AmOfInitPo - 3) >> 1;
		long long AmOfLoops = (AmOfInitPo - 3) >> 1;
		//for(int i=0; i<AmOfLoops; i++)
		for(long long i=0; i<AmOfLoops; i++)
		{
			Sum1XRe += *(t++); Sum1XIm += *(t++); Sum1ZRe += *(t++); Sum1ZIm += *(t++);
			Sum2XRe += *(t++); Sum2XIm += *(t++); Sum2ZRe += *(t++); Sum2ZIm += *(t++);
		}
		Sum1XRe += *(t++); Sum1XIm += *(t++); Sum1ZRe += *(t++); Sum1ZIm += *t;
	}
	else
	{
		Sum1XRe += *(t++); Sum1XIm += *(t++); Sum1ZRe += *(t++); Sum1ZIm += *t;
	}

	//long Np = AmOfInitPo - 1;
	long long Np = AmOfInitPo - 1;
	double sStep = (sEnd - sStart)/double(Np), s;

	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;
	double ActNormConstsStep = ActNormConst*sStep;
	double IntXRe = ActNormConstsStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + sStep*wdFxRe);
	double IntXIm = ActNormConstsStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + sStep*wdFxIm);
	double IntZRe = ActNormConstsStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + sStep*wdFzRe);
	double IntZIm = ActNormConstsStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + sStep*wdFzIm);
	double SqNorm = IntXRe*IntXRe + IntXIm*IntXIm + IntZRe*IntZRe + IntZIm*IntZIm;

	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;
	char NearField = (DistrInfoDat.CoordOrAngPresentation == CoordPres);

	double xObs = ObsCoor.x, yObs = ObsCoor.y, zObs = ObsCoor.z, GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz, Ph, Ax, Az, CosPh, SinPh;
	double AngPhConst, Two_xObs, Two_zObs;
	AngPhConst = GmEm2 + xObs*xObs + zObs*zObs; Two_xObs = 2.*xObs; Two_zObs = 2.*zObs;

	double *pBtx, *pX, *pIntBtxE2, *pBx, *pBtz, *pZ, *pIntBtzE2, *pBz;
	double **TrjPtrs[] = {&pBtx, &pX, &pIntBtxE2, &pBx, &pBtz, &pZ, &pIntBtzE2, &pBz};

	double BtxLoc, xLoc, IntBtxE2Loc, BxLoc, BtzLoc, zLoc, IntBtzE2Loc, BzLoc;
	int result;

	//char ExtraPassForAnyCase = 0;

	int LevelNo = ThisLevNo;
	char NotFinishedYet = 1;
	while(NotFinishedYet)
	{
		Sum2XRe += Sum1XRe; Sum2XIm += Sum1XIm; Sum2ZRe += Sum1ZRe; Sum2ZIm += Sum1ZIm; 
		Sum1XRe = Sum1XIm = Sum1ZRe = Sum1ZIm = 0.;
		char ThisMayBeTheLastLoop = 1;

		PhPrev = PhInit; LevelNo++;
		double HalfStep = 0.5*sStep;
		s = sStart + HalfStep;

		if(LevelNo <= MaxLevelForMeth_10_11)
		{
			if(result = FillNextLevelPart(LevelNo, s, sEnd - HalfStep, Np, TrjPtrs)) return result;
		}
				double DPhMax = 0.;

		//for(long i=0; i<Np; i++)
		for(long long i=0; i<Np; i++)
		{
			if(LevelNo > MaxLevelForMeth_10_11)
			{
				pBtx = &BtxLoc; pX = &xLoc; pIntBtxE2 = &IntBtxE2Loc; pBx = &BxLoc;
				pBtz = &BtzLoc; pZ = &zLoc; pIntBtzE2 = &IntBtzE2Loc; pBz = &BzLoc;

				TrjDatPtr->CompTrjDataDerivedAtPoint(s, *pBtx, *pX, *pIntBtxE2, *pBtz, *pZ, *pIntBtzE2);
			}

			if(NearField)
			{
				One_d_ymis = 1./(yObs - s);
				xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
				Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;

				double LongTerm = *pIntBtxE2 + *pIntBtzE2;
				//double a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
				//double a = a0*One_d_ymis;
				//Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
				//OC_test
					double a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
					Ph = PIm10e9_d_Lamb*(s*GmEm2 + a0);
				//end OC_test

				Ax = (*pBtx - Nx)*One_d_ymis; Az = (*pBtz - Nz)*One_d_ymis;
			}
			else
			{
				Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_xObs*(*pX) + Two_zObs*(*pZ)));
				Ax = *pBtx - xObs; Az = *pBtz - zObs;
			}

			if(LevelNo <= MaxLevelForMeth_10_11)
			{
				pBtx++; pX++; pIntBtxE2++; pBtz++; pZ++; pIntBtzE2++;
			}

			//CosAndSin(Ph, CosPh, SinPh);
			CosPh = cos(Ph); SinPh = sin(Ph);

			Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; 
			s += sStep;

			if(Ph - PhPrev > PI) ThisMayBeTheLastLoop = 0;

				double dPh = Ph - PhPrev;
				if(dPh > DPhMax) DPhMax = dPh;

			PhPrev = Ph;
		}
		double ActNormConstHalfStep = ActNormConst*HalfStep;

		double LocIntXRe = ActNormConstHalfStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + HalfStep*wdFxRe);
		double LocIntXIm = ActNormConstHalfStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + HalfStep*wdFxIm);
		double LocIntZRe = ActNormConstHalfStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + HalfStep*wdFzRe);
		double LocIntZIm = ActNormConstHalfStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + HalfStep*wdFzIm);
		double LocSqNorm = LocIntXRe*LocIntXRe + LocIntXIm*LocIntXIm + LocIntZRe*LocIntZRe + LocIntZIm*LocIntZIm;

		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			//char SharplyGoesDown = (LocSqNorm < 0.2*SqNorm);
			char NotFinishedYetFirstTest = (TestVal > sIntegRelPrec*LocSqNorm);

			if(!NotFinishedYetFirstTest)
			{
				//if(ExtraPassForAnyCase || SharplyGoesDown) NotFinishedYet = 0;
				//else ExtraPassForAnyCase = 1;

				NotFinishedYet = 0;
			}
		}

		if(NotFinishedYet)
		{
			if(Np > NpOnLevelMaxNoResult) return CAN_NOT_COMPUTE_RADIATION_INTEGRAL;
		}

		IntXRe = LocIntXRe; IntXIm = LocIntXIm; IntZRe = LocIntZRe; IntZIm = LocIntZIm; SqNorm = LocSqNorm;
		sStep = HalfStep; Np *= 2;
	}
	OutIntXRe = IntXRe; OutIntXIm = IntXIm; OutIntZRe = IntZRe; OutIntZIm = IntZIm; 
	return 0;
}

//*************************************************************************

int srTRadInt::RadIntegrationAuto2(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, srTEFourier* pEwNormDer)
{
// Streering Parameters Section
	const int NpOnZeroLev = 5; // Must be non-even!
	
	//const double DelPhiSwitchToNum = 100.; // To switch to purely numerical method
	int AmOfLoops = 7;
	//int AmOfLoops = 8; //OCTEST //7;

	double RelTolForAn = sIntegRelPrec*0.45;

	double DeltaPhCrit = 0.5; // To prefer Num
	double GamMult = 1.; // To steer
// End Streering Parameters Section

	int TotAmOfInitPo = NpOnZeroLev + (NpOnZeroLev - 1)*((1 << (AmOfLoops - 1)) - 1);

	const int NpOnZeroLev_mi_1 = NpOnZeroLev - 1;
	double GenInterv = sIntegFin - sIntegStart;
	//int NpOnLevel = NpOnZeroLev;
	long long NpOnLevel = NpOnZeroLev;

	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*ObsCoor.Lamb*0.80654658E-03 : NormalizingConst/ObsCoor.Lamb;
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoor.Lamb : PIm10e6*1000./ObsCoor.Lamb;

	double xObs = ObsCoor.x, yObs = ObsCoor.y, zObs = ObsCoor.z, GmEm2 = TrjDatPtr->EbmDat.GammaEm2, TwoPIm10e9_d_Lamb = 2.*PIm10e9_d_Lamb;
	double GamFact = GamMult*GmEm2; // To steer
	double ConBtx = TrjDatPtr->BetaNormConst, ConBtz = -TrjDatPtr->BetaNormConst;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz, Btx_mi_Nx, Btz_mi_Nz, BetaFact;
	double ConBtxBz, ConBtzBx, ConBtxBzpAx, ConBtzBxpAz, SinPh, CosPh;
	double AngPhConst, Two_xObs, Two_zObs;
	double d2Phds2PartX, d2Phds2PartZ;

	char NearField = (DistrInfoDat.CoordOrAngPresentation == CoordPres);
	if(NearField)
	{
		AngPhConst = GmEm2 + xObs*xObs + zObs*zObs; Two_xObs = 2.*xObs; Two_zObs = 2.*zObs;
	}

	double *pBtx, *pX, *pIntBtxE2, *pBx, *pBtz, *pZ, *pIntBtzE2, *pBz;
	double **TrjPtrs[] = {&pBtx, &pX, &pIntBtxE2, &pBx, &pBtz, &pZ, &pIntBtzE2, &pBz};

	int result;
	if(result = FillNextLevelPart(0, sIntegStart, sIntegFin, NpOnLevel, TrjPtrs)) return result;

	//srTRadIntervVal *tComputedRadIntervAuto = ComputedRadIntervAuto;
	//int AmOfComputedRadIntervAuto = 0;

	double SumIntXRe=0, SumIntXIm=0, SumIntZRe=0, SumIntZIm=0;

	double *pIntXRe = IntXReArrAuto, *pIntXIm = IntXImArrAuto, *pIntZRe = IntZReArrAuto, *pIntZIm = IntZImArrAuto;
	double *pPh = PhArrAuto, *pdPhds = dPhdsArrAuto, *pd2Phds2 = d2Phds2ArrAuto;
	double *pAx = AxArrAuto, *pdAxds = dAxdsArrAuto, *pAz = AzArrAuto, *pdAzds = dAzdsArrAuto;
	char *pPointIsGoodForAnAuto = PointIsGoodForAnAuto;

	double *pd2Phds2ddPhdsE2 = d2Phds2ddPhdsE2ArrAuto;

	for(int kk=0; kk<=TotAmOfInitPo; kk++)
	{
		*(pPointIsGoodForAnAuto++) = 0;
	}
	pPointIsGoodForAnAuto = PointIsGoodForAnAuto;

	double s = sIntegStart;
	double sStep = GenInterv/NpOnZeroLev_mi_1;
	double HalfStep = 0.5*sStep;
	int LevelNo = 0;

	//int TotNp = NpOnLevel;
	long long TotNp = NpOnLevel;

	double t1xd = 0., t2xd = 0.,t1zd = 0., t2zd = 0., d2Phds2_d_dPhds = 0.; //OC190815

	int AmOfLoops_mi_1 = AmOfLoops - 1;
	for(int k=0; k<AmOfLoops; k++)
	{
		//for(int i=0; i<NpOnLevel; i++)
		for(long long i=0; i<NpOnLevel; i++)
		{
			if(NearField) 
			{
				One_d_ymis = 1./(yObs - s); 
				xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
				Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;

				double LongTerm = *pIntBtxE2 + *pIntBtzE2;
				//double a0 = LongTerm*(1. + 0.25*LongTerm*One_d_ymis) + xObs_mi_x*Nx + zObs_mi_z*Nz;
				//double a = a0*One_d_ymis;
				//*pPh = PIm10e9_d_Lamb*(s*GmEm2 + a0*(1 + a*(-0.25 + a*(0.125 - 0.078125*a))));
				//OC_test
					double a0 = LongTerm + xObs_mi_x*Nx + zObs_mi_z*Nz;
					*pPh = PIm10e9_d_Lamb*(s*GmEm2 + a0);
				//end OC_test

				Btx_mi_Nx = *pBtx - Nx; Btz_mi_Nz = *pBtz - Nz;

				BetaFact = Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz;

				*pAx = Btx_mi_Nx*One_d_ymis; *pAz = Btz_mi_Nz*One_d_ymis;
				*pdPhds = PIm10e9_d_Lamb*(GmEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
				ConBtxBz = ConBtx*(*pBz); ConBtzBx = ConBtz*(*pBx);
				ConBtxBzpAx = ConBtxBz + *pAx; ConBtzBxpAz = ConBtzBx + *pAz;

				d2Phds2PartX = TwoPIm10e9_d_Lamb*Btx_mi_Nx*ConBtxBzpAx;
				d2Phds2PartZ = TwoPIm10e9_d_Lamb*Btz_mi_Nz*ConBtzBxpAz;
				*pd2Phds2 = d2Phds2PartX + d2Phds2PartZ;

				*pdAxds = (2.*(*pAx) + ConBtxBz)*One_d_ymis; 
				*pdAzds = (2.*(*pAz) + ConBtzBx)*One_d_ymis;
			}
			else
			{
				*pPh = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_xObs*(*pX) + Two_zObs*(*pZ)));
				*pAx = *pBtx - xObs; *pAz = *pBtz - zObs;

				BetaFact = (*pAx)*(*pAx) + (*pAz)*(*pAz);

				*pdPhds = PIm10e9_d_Lamb*(GmEm2 + (*pAx)*(*pAx) + (*pAz)*(*pAz));
				ConBtxBz = ConBtx*(*pBz); ConBtzBx = ConBtz*(*pBx);

				d2Phds2PartX = TwoPIm10e9_d_Lamb*(*pAx)*ConBtxBz;
				d2Phds2PartZ = TwoPIm10e9_d_Lamb*(*pAz)*ConBtzBx;
				*pd2Phds2 = d2Phds2PartX + d2Phds2PartZ;

				*pdAxds = ConBtxBz; *pdAzds = ConBtzBx;
			}

			double One_d_dPhds = 1./(*pdPhds);
			double One_d_dPhdsE2 = One_d_dPhds*One_d_dPhds;

			double d2Phds2ddPhdsE2 = (*pd2Phds2)*One_d_dPhdsE2;
			*pd2Phds2ddPhdsE2 = d2Phds2ddPhdsE2;

			double d2Phds2PartXddPhdsE2 = d2Phds2PartX*One_d_dPhdsE2;
			double d2Phds2PartZddPhdsE2 = d2Phds2PartZ*One_d_dPhdsE2;

			char p2OK = (RelTolForAn > ::fabs(d2Phds2PartXddPhdsE2)) && (RelTolForAn > ::fabs(d2Phds2PartZddPhdsE2));

			char BetterUseNum = 0;
			if(i!=0)
			{
				if((*pPh - *(pPh-1)) < DeltaPhCrit) BetterUseNum = 1;
			}

			//*pPointIsGoodForAnAuto = (p2OK && (!BetterUseNum) && ((BetaFact > GamFact) || ((*pBx == 0.) && (*pBz == 0.))));
			//OC190815
			t1xd = *pAx; t1zd = *pAz;
			d2Phds2_d_dPhds = (*pd2Phds2)*One_d_dPhds;
			t2xd = ((*pdAxds) - t1xd*d2Phds2_d_dPhds)*One_d_dPhds;
			t2zd = ((*pdAzds) - t1zd*d2Phds2_d_dPhds)*One_d_dPhds;
			*pPointIsGoodForAnAuto = (p2OK && (!BetterUseNum) && (::fabs(t2xd) <= RelTolForAn*(::fabs(t1xd))) && (::fabs(t2zd) <= RelTolForAn*(::fabs(t1zd))));
			//OC190815 End 

			s += sStep;

			if(k==0)
			{
				pPh++; pdPhds++; pd2Phds2++; pAx++; pdAxds++; pAz++; pdAzds++; 
				pIntXRe++; pIntXIm++; pIntZRe++; pIntZIm++;
				pPointIsGoodForAnAuto++;
				pd2Phds2ddPhdsE2++;
			}
			else
			{
				pPh+=2; pdPhds+=2; pd2Phds2+=2; pAx+=2; pdAxds+=2; pAz+=2; pdAzds+=2;
				pIntXRe+=2; pIntXIm+=2; pIntZRe+=2; pIntZIm+=2;
				pPointIsGoodForAnAuto+=2;
				pd2Phds2ddPhdsE2+=2;
			}
			pBtx++, pBtz++, pX++, pZ++, pIntBtxE2++, pIntBtzE2++; pBx++; pBz++;
		}

		if(k==0)
		{
			NpOnLevel--;
		}
		else
		{
			sStep = HalfStep; HalfStep *= 0.5; NpOnLevel *= 2;
		}

		if(k != AmOfLoops_mi_1)
		{
			//int TotNp_mi_1 = TotNp-1;
			//int ia = TotNp_mi_1;
			long long TotNp_mi_1 = TotNp-1;
			long long ia = TotNp_mi_1;
			//for(int ii=0; ii<TotNp_mi_1; ii++)
			for(long long ii=0; ii<TotNp_mi_1; ii++)
			{
				//int ib = 2*ia;
				long long ib = 2*ia;
				IntXReArrAuto[ib] = IntXReArrAuto[ia];
				IntXImArrAuto[ib] = IntXImArrAuto[ia];
				IntZReArrAuto[ib] = IntZReArrAuto[ia];
				IntZImArrAuto[ib] = IntZImArrAuto[ia];
				PhArrAuto[ib] = PhArrAuto[ia];
				dPhdsArrAuto[ib] = dPhdsArrAuto[ia];
				d2Phds2ArrAuto[ib] = d2Phds2ArrAuto[ia];
				AxArrAuto[ib] = AxArrAuto[ia];
				dAxdsArrAuto[ib] = dAxdsArrAuto[ia];
				AzArrAuto[ib] = AzArrAuto[ia];
				dAzdsArrAuto[ib] = dAzdsArrAuto[ia];

				d2Phds2ddPhdsE2ArrAuto[ib] = d2Phds2ddPhdsE2ArrAuto[ia];

				char BufChar = PointIsGoodForAnAuto[ia];
				PointIsGoodForAnAuto[ib] = BufChar;
				if(BufChar==7)
				{
					PointIsGoodForAnAuto[ib - 1] = PointIsGoodForAnAuto[ib + 1] = 7;
				}
				else PointIsGoodForAnAuto[ib - 1] = 0;
				ia--;
			}

			s = sIntegStart + HalfStep; 
			LevelNo++; 

			if(result = FillNextLevelPart(LevelNo, sIntegStart + HalfStep, sIntegFin - HalfStep, NpOnLevel, TrjPtrs)) return result;

			pIntXRe = IntXReArrAuto+1; pIntXIm = IntXImArrAuto+1; pIntZRe = IntZReArrAuto+1; pIntZIm = IntZImArrAuto+1; 
			pPh = PhArrAuto+1; pdPhds = dPhdsArrAuto+1; pd2Phds2 = d2Phds2ArrAuto+1;
			pAx = AxArrAuto+1; pdAxds = dAxdsArrAuto+1; pAz = AzArrAuto+1; pdAzds = dAzdsArrAuto+1;
			pPointIsGoodForAnAuto = PointIsGoodForAnAuto+1;
			pd2Phds2ddPhdsE2 = d2Phds2ddPhdsE2ArrAuto + 1;

			TotNp += NpOnLevel;
		}
	}

	double Stp = (sIntegFin - sIntegStart)/(TotAmOfInitPo - 1);

	int jCurrentStartAn = 0, jCurrentStartNum = 0;

	if(!(*PointIsGoodForAnAuto)) *(PointIsGoodForAnAuto+1) = 0;
	if(!PointIsGoodForAnAuto[TotAmOfInitPo-1]) PointIsGoodForAnAuto[TotAmOfInitPo-2] = 0;

	pd2Phds2 = d2Phds2ArrAuto + 1;
	pPointIsGoodForAnAuto = PointIsGoodForAnAuto + 1;

	srTRadIntervVal *pAnRadInt = ComputedRadIntervAuto, *pNumRadInt = NumRadIntervAuto;
	int AmOfNumArrays = 0, AmOfAnArrays = 0;
	int TotAmOfInitPo_mi_1 = TotAmOfInitPo - 1, TotAmOfInitPo_mi_2 = TotAmOfInitPo - 2;

	char AllPointsAreGoodForAn = 1;

	for(int j=1; j<TotAmOfInitPo; j++)
	{
		if(*pPointIsGoodForAnAuto)
		{
			if(!(*(pPointIsGoodForAnAuto-1)))
			{
				if(((j < TotAmOfInitPo_mi_2) && *(pPointIsGoodForAnAuto+1) && *(pPointIsGoodForAnAuto+2))
					|| (j == TotAmOfInitPo_mi_2))
				{
					int AmOfPoInNumArray = j - jCurrentStartNum + 1;
					double HalfAmOfPoInNumArray = 0.5*AmOfPoInNumArray;
					char AmOfPoInNumArrayIsEven = ((HalfAmOfPoInNumArray - int(HalfAmOfPoInNumArray + 1.E-06)) < 0.4);
					int jFi = AmOfPoInNumArrayIsEven? (j + 1) : j;

					pNumRadInt->InitLevelNo = jCurrentStartNum;
					pNumRadInt->SecondInt = jFi;
					pNumRadInt->sStart = sIntegStart + Stp*jCurrentStartNum;
					pNumRadInt->sEnd = sIntegStart + Stp*jFi;
					pNumRadInt++; AmOfNumArrays++;
					jCurrentStartAn = jFi;
				}
			}
		}
		else
		{
			if(((j >= 3) && (*(pPointIsGoodForAnAuto-1) && *(pPointIsGoodForAnAuto-2) && *(pPointIsGoodForAnAuto-3)))
				|| ((j < 3) && *(pPointIsGoodForAnAuto-1) && *(pPointIsGoodForAnAuto-2)))
			{
				int jmi1 = j - 1;
				if((AmOfAnArrays == 0) || ((AmOfAnArrays > 0) && ((jmi1 - (pAnRadInt-1)->SecondInt) > 1)))
				{
					if(((jCurrentStartAn > 0) && ((jmi1 - jCurrentStartAn) > 1))
						|| ((jCurrentStartAn == 0) && (jmi1 >= 1)))
					{
						pAnRadInt->InitLevelNo = jCurrentStartAn;
						pAnRadInt->SecondInt = jmi1;
						pAnRadInt->sStart = sIntegStart + Stp*jCurrentStartAn;
						pAnRadInt->sEnd = sIntegStart + Stp*jmi1;
						pAnRadInt++; AmOfAnArrays++;
						jCurrentStartNum = jmi1;
					}
					else
					{
						pNumRadInt--; AmOfNumArrays--;
					}
				}
			}
			AllPointsAreGoodForAn = 0;
		}
		pPointIsGoodForAnAuto++;
		pd2Phds2++;
	}
	if(PointIsGoodForAnAuto[TotAmOfInitPo_mi_1])
	{
		if((AmOfNumArrays > 0) && (PointIsGoodForAnAuto[TotAmOfInitPo_mi_2]))
		{
			if((pNumRadInt-1)->SecondInt != TotAmOfInitPo_mi_1)
			{
				pAnRadInt->InitLevelNo = (pNumRadInt-1)->SecondInt;
				pAnRadInt->SecondInt = TotAmOfInitPo_mi_1;
				pAnRadInt->sStart = sIntegStart + Stp*(pAnRadInt->InitLevelNo);
				pAnRadInt->sEnd = sIntegFin;
				AmOfAnArrays++;
			}
		}
		else if((!PointIsGoodForAnAuto[TotAmOfInitPo_mi_2]) && (AmOfAnArrays > 0))
		{
			int AmOfPoInNumArray = TotAmOfInitPo_mi_1 - jCurrentStartNum + 1;
			double HalfAmOfPoInNumArray = 0.5*AmOfPoInNumArray;
			char AmOfPoInNumArrayIsEven = ((HalfAmOfPoInNumArray - int(HalfAmOfPoInNumArray + 1.E-06)) < 0.4);
			int jFi = TotAmOfInitPo_mi_1;

			if(AmOfPoInNumArrayIsEven)
			{
				//int PrevAnLen = (pAnRadInt-1)->SecondInt - (pAnRadInt-1)->InitLevelNo;
				long long PrevAnLen = (pAnRadInt-1)->SecondInt - (pAnRadInt-1)->InitLevelNo;
				if(PrevAnLen == 1)
				{
					AmOfAnArrays--;
					pNumRadInt->InitLevelNo = (pAnRadInt-1)->InitLevelNo;
				}
				else
				{
					((pAnRadInt-1)->SecondInt)--;
					(pAnRadInt-1)->sEnd -= Stp;

					pNumRadInt->InitLevelNo = (pAnRadInt-1)->SecondInt;
				}
			}
			else
			{
				pNumRadInt->InitLevelNo = (pAnRadInt-1)->SecondInt;
			}
			pNumRadInt->SecondInt = jFi;
			pNumRadInt->sStart = sIntegStart + Stp*(pNumRadInt->InitLevelNo);
			pNumRadInt->sEnd = sIntegStart + Stp*jFi;
			AmOfNumArrays++;
		}
	}
	else
	{
		int AmOfPoInNumArray = TotAmOfInitPo - jCurrentStartNum;
		double HalfAmOfPoInNumArray = 0.5*AmOfPoInNumArray;
		char AmOfPoInNumArrayIsEven = ((HalfAmOfPoInNumArray - int(HalfAmOfPoInNumArray + 1.E-06)) < 0.4);

		if(!AmOfPoInNumArrayIsEven)
		{
			pNumRadInt->InitLevelNo = jCurrentStartNum;
			pNumRadInt->SecondInt = TotAmOfInitPo_mi_1;
			pNumRadInt->sStart = sIntegStart + Stp*jCurrentStartNum;
			pNumRadInt->sEnd = sIntegFin;
			AmOfNumArrays++;
		}
		else
		{
			int jCurrentStartNum_mi_1 = jCurrentStartNum - 1;

			pNumRadInt->InitLevelNo = jCurrentStartNum_mi_1;
			pNumRadInt->SecondInt = TotAmOfInitPo_mi_1;
			pNumRadInt->sStart = sIntegStart + Stp*jCurrentStartNum_mi_1;
			pNumRadInt->sEnd = sIntegFin;
			AmOfNumArrays++;

			srTRadIntervVal *pAnRadInt_mi_1 = pAnRadInt - 1;
			if((pAnRadInt_mi_1->SecondInt - pAnRadInt_mi_1->InitLevelNo) > 1)
			{
				(pAnRadInt_mi_1->SecondInt)--;
				pAnRadInt_mi_1->sEnd = sIntegStart + Stp*(pAnRadInt_mi_1->SecondInt);
			}
			else
			{
				AmOfAnArrays--;
			}
		}
		AllPointsAreGoodForAn = 0;
	}

	if(AllPointsAreGoodForAn)
	{
		pd2Phds2ddPhdsE2 = d2Phds2ddPhdsE2ArrAuto;
		double Max = 0.;
		int kMax = -1;
		for(int k=0; k<TotAmOfInitPo; k++)
		{
			double TestMax = ::fabs(*(pd2Phds2ddPhdsE2++));
			if(TestMax > Max) { Max = TestMax; kMax = k;}
		}
		if(kMax <= 1)
		{
			pNumRadInt->InitLevelNo = 0;
			pNumRadInt->SecondInt = 2;
			pNumRadInt->sStart = sIntegStart;
			pNumRadInt->sEnd = sIntegStart + 2.*Stp;
			
			pAnRadInt->InitLevelNo = 2;
			pAnRadInt->SecondInt = TotAmOfInitPo_mi_1;
			pAnRadInt->sStart = sIntegStart + 2.*Stp;
			pAnRadInt->sEnd = sIntegFin;
			AmOfAnArrays = 1;
		}
		else if(kMax >= TotAmOfInitPo_mi_2)
		{
			int TotAmOfInitPo_mi_3 = TotAmOfInitPo - 3;
			pAnRadInt->InitLevelNo = 0;
			pAnRadInt->SecondInt = TotAmOfInitPo_mi_3;
			pAnRadInt->sStart = sIntegStart;
			pAnRadInt->sEnd = sIntegStart + TotAmOfInitPo_mi_3*Stp;
			AmOfAnArrays = 1;
 
			pNumRadInt->InitLevelNo = TotAmOfInitPo_mi_3;
			pNumRadInt->SecondInt = TotAmOfInitPo_mi_1;
			pNumRadInt->sStart = pAnRadInt->sEnd;
			pNumRadInt->sEnd = sIntegFin;
		}
		else
		{
			pAnRadInt->InitLevelNo = 0;
			pAnRadInt->SecondInt = kMax - 1;
			pAnRadInt->sStart = sIntegStart;
			pAnRadInt->sEnd = sIntegStart + Stp*(pAnRadInt->SecondInt);

			pNumRadInt->InitLevelNo = pAnRadInt->SecondInt;
			pNumRadInt->SecondInt = kMax + 1;
			pNumRadInt->sStart = pAnRadInt->sEnd;
			pNumRadInt->sEnd = sIntegStart + Stp*(pNumRadInt->SecondInt);

			++pAnRadInt;
			pAnRadInt->InitLevelNo = pNumRadInt->SecondInt;
			pAnRadInt->SecondInt = TotAmOfInitPo_mi_1;
			pAnRadInt->sStart = pNumRadInt->sEnd;
			pAnRadInt->sEnd = sIntegFin;
			AmOfAnArrays = 2;
		}
		AmOfNumArrays = 1;
	}

	double BufIntXRe, BufIntXIm, BufIntZRe, BufIntZIm;
	double FunArr[513*4], EdgeDerArr[8];
	if(AmOfAnArrays == 0)
	{
		double *t = FunArr, *td = EdgeDerArr;
		pPh = PhArrAuto; pdPhds = dPhdsArrAuto; pAx = AxArrAuto; pdAxds = dAxdsArrAuto; pAz = AzArrAuto; pdAzds = dAzdsArrAuto;
		double InitPh = *pPh;
		//for(int i=0; i<TotNp; i++)
		for(long long i=0; i<TotNp; i++)
		{
			CosAndSin(*pPh, CosPh, SinPh);
			if((i==0) || (i==(TotNp-1)))
			{
				double dPhdsSinPh = (*pdPhds)*SinPh, dPhdsCosPh = (*pdPhds)*CosPh;
				*(td++) = (*pdAxds)*CosPh - (*pAx)*dPhdsSinPh;
				*(td++) = (*pdAxds)*SinPh + (*pAx)*dPhdsCosPh;
				*(td++) = (*pdAzds)*CosPh - (*pAz)*dPhdsSinPh;
				*(td++) = (*pdAzds)*SinPh + (*pAz)*dPhdsCosPh;
			}
			*(t++) = (*pAx)*CosPh; *(t++) = (*pAx)*SinPh;
			*(t++) = (*pAz)*CosPh; *(t++) = (*pAz)*SinPh;

			pPh++; pdPhds++; pAx++; pdAxds++; pAz++; pdAzds++;
		}
		*t = InitPh;
		if(result = RadIntegrationAuto1M(sIntegStart, sIntegFin, FunArr, EdgeDerArr, TotNp, LevelNo, BufIntXRe, BufIntXIm, BufIntZRe, BufIntZIm)) return result;
		OutIntXRe += BufIntXRe; OutIntXIm += BufIntXIm; OutIntZRe += BufIntZRe; OutIntZIm += BufIntZIm;
		return 0;
	}

	pAnRadInt = ComputedRadIntervAuto;

	for(int kAn=0; kAn<AmOfAnArrays; kAn++)
	{
		//int IndSt = pAnRadInt->InitLevelNo;
		//int IndFi = pAnRadInt->SecondInt;
		long long IndSt = pAnRadInt->InitLevelNo;
		long long IndFi = pAnRadInt->SecondInt;

		double CosPh, SinPh;
		double Ph = PhArrAuto[IndFi], dPhds = dPhdsArrAuto[IndFi], d2Phds2 = d2Phds2ArrAuto[IndFi];
		double Ax = AxArrAuto[IndFi], dAxds = dAxdsArrAuto[IndFi], Az = AzArrAuto[IndFi], dAzds = dAzdsArrAuto[IndFi];
		double One_d_dPhds = 1./dPhds;
		double One_d_dPhdsE2 = One_d_dPhds*One_d_dPhds;
		double dAxdsddPhds = dAxds*One_d_dPhds;
		double dAzdsddPhds = dAzds*One_d_dPhds;
		double d2Phds2ddPhdsE2 = d2Phds2*One_d_dPhdsE2;
		double t1x = dAxdsddPhds - Ax*d2Phds2ddPhdsE2;
		double t1z = dAzdsddPhds - Az*d2Phds2ddPhdsE2;
		CosAndSin(Ph, CosPh, SinPh);
		double IntXRe = One_d_dPhds*(t1x*CosPh + Ax*SinPh);
		double IntXIm = One_d_dPhds*(t1x*SinPh - Ax*CosPh);
		double IntZRe = One_d_dPhds*(t1z*CosPh + Az*SinPh);
		double IntZIm = One_d_dPhds*(t1z*SinPh - Az*CosPh);

		Ph = PhArrAuto[IndSt]; dPhds = dPhdsArrAuto[IndSt]; d2Phds2 = d2Phds2ArrAuto[IndSt];
		Ax = AxArrAuto[IndSt]; dAxds = dAxdsArrAuto[IndSt]; Az = AzArrAuto[IndSt]; dAzds = dAzdsArrAuto[IndSt];
		One_d_dPhds = 1./dPhds;
		One_d_dPhdsE2 = One_d_dPhds*One_d_dPhds;
		dAxdsddPhds = dAxds*One_d_dPhds;
		dAzdsddPhds = dAzds*One_d_dPhds;
		d2Phds2ddPhdsE2 = d2Phds2*One_d_dPhdsE2;
		t1x = dAxdsddPhds - Ax*d2Phds2ddPhdsE2;
		t1z = dAzdsddPhds - Az*d2Phds2ddPhdsE2;
		CosAndSin(Ph, CosPh, SinPh);
		IntXRe -= One_d_dPhds*(t1x*CosPh + Ax*SinPh);
		IntXIm -= One_d_dPhds*(t1x*SinPh - Ax*CosPh);
		IntZRe -= One_d_dPhds*(t1z*CosPh + Az*SinPh);
		IntZIm -= One_d_dPhds*(t1z*SinPh - Az*CosPh);

		SumIntXRe += IntXRe; SumIntXIm += IntXIm; SumIntZRe += IntZRe; SumIntZIm += IntZIm;
		pAnRadInt++;
	}
	SumIntXRe *= ActNormConst; SumIntXIm *= ActNormConst; 
	SumIntZRe *= ActNormConst; SumIntZIm *= ActNormConst;

	pNumRadInt = NumRadIntervAuto;
	for(int n=0; n<AmOfNumArrays; n++)
	{
		//int TotAmOfPo = pNumRadInt->SecondInt - pNumRadInt->InitLevelNo + 1;
		long long TotAmOfPo = pNumRadInt->SecondInt - pNumRadInt->InitLevelNo + 1;

		double *t = FunArr, *td = EdgeDerArr;
		//int StNo = pNumRadInt->InitLevelNo, FiNo = pNumRadInt->SecondInt;
		long long StNo = pNumRadInt->InitLevelNo, FiNo = pNumRadInt->SecondInt;
		//int FiNo_mi_1 = FiNo - 1;

		pPh = PhArrAuto+StNo; pdPhds = dPhdsArrAuto+StNo; 
		pAx = AxArrAuto+StNo; pdAxds = dAxdsArrAuto+StNo; 
		pAz = AzArrAuto+StNo; pdAzds = dAzdsArrAuto+StNo;
		double InitPh = *pPh;

		//for(int jj=StNo; jj<=FiNo; jj++)
		for(long long jj=StNo; jj<=FiNo; jj++)
		{
			CosAndSin(*pPh, CosPh, SinPh);
			if((jj==StNo) || (jj==FiNo))
			{
				double dPhdsSinPh = (*pdPhds)*SinPh, dPhdsCosPh = (*pdPhds)*CosPh;
				*(td++) = (*pdAxds)*CosPh - (*pAx)*dPhdsSinPh;
				*(td++) = (*pdAxds)*SinPh + (*pAx)*dPhdsCosPh;
				*(td++) = (*pdAzds)*CosPh - (*pAz)*dPhdsSinPh;
				*(td++) = (*pdAzds)*SinPh + (*pAz)*dPhdsCosPh;
			}
			*(t++) = (*pAx)*CosPh; *(t++) = (*pAx)*SinPh;
			*(t++) = (*pAz)*CosPh; *(t++) = (*pAz)*SinPh;

			pPh++; pdPhds++; pAx++; pdAxds++; pAz++; pdAzds++;
		}
		*t = InitPh;
		if(result = RadIntegrationAuto1M(pNumRadInt->sStart, pNumRadInt->sEnd, FunArr, EdgeDerArr, TotAmOfPo, LevelNo, BufIntXRe, BufIntXIm, BufIntZRe, BufIntZIm)) return result;
		SumIntXRe += BufIntXRe; SumIntXIm += BufIntXIm; SumIntZRe += BufIntZRe; SumIntZIm += BufIntZIm;
		pNumRadInt++;
	}
	OutIntXRe += SumIntXRe; OutIntXIm += SumIntXIm; OutIntZRe += SumIntZRe; OutIntZIm += SumIntZIm;
	return 0;
}

//*************************************************************************

//int srTRadInt::FillNextLevelPart(int LevelNo, double sStart, double sEnd, long Np, double*** TrjPtrs)
int srTRadInt::FillNextLevelPart(int LevelNo, double sStart, double sEnd, long long Np, double*** TrjPtrs)
{
	srTStNoFiNoVect& StNoFiNoVect = StNoFiNoVectArr[LevelNo];
	//int TotalNpOnLevel = (LevelNo > 0)? (*AmOfPointsOnLevel - 1)*(1 << (LevelNo - 1)) : Np;
	long long TotalNpOnLevel = (LevelNo > 0)? (*AmOfPointsOnLevel - 1)*(((long long)1) << ((long long)(LevelNo - 1))) : Np;

	//int GenStartOffset = 0;
	long long GenStartOffset = 0;

	if(TotalNpOnLevel != Np)
	{
		double Two_h = (sIntegFin - sIntegStart)/TotalNpOnLevel;
		double h = 0.5*Two_h;
		//int iSt = int((sStart - sIntegStart - h)/Two_h + 1.E-08);
		//int iFi = iSt + Np - 1;
		long long iSt = (long long)((sStart - sIntegStart - h)/Two_h + 1.E-08);
		long long iFi = iSt + Np - 1;

		if(NumberOfLevelsFilled > LevelNo) 
		{
			GenStartOffset = iSt; goto SettingPointers;
		}

		if(StNoFiNoVect.empty())
		{
			double* BasePtr = new double[Np*8];
			if(BasePtr == 0) 
			{ 
				//pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;
				return MEMORY_ALLOCATION_FAILURE;
			}

			BtxArrP[LevelNo] = BasePtr;
			BasePtr += Np; XArrP[LevelNo] = BasePtr;
			BasePtr += Np; IntBtxE2ArrP[LevelNo] = BasePtr;
			BasePtr += Np; BxArrP[LevelNo] = BasePtr;
			BasePtr += Np; BtzArrP[LevelNo] = BasePtr;
			BasePtr += Np; ZArrP[LevelNo] = BasePtr;
			BasePtr += Np; IntBtzE2ArrP[LevelNo] = BasePtr;
			BasePtr += Np; BzArrP[LevelNo] = BasePtr;

			TrjDatPtr->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevelNo], BtzArrP[LevelNo], XArrP[LevelNo], ZArrP[LevelNo], IntBtxE2ArrP[LevelNo], IntBtzE2ArrP[LevelNo], BxArrP[LevelNo], BzArrP[LevelNo]);

			srTStNoFiNo StNoFiNo(iSt, iFi, 0);
			StNoFiNoVect.push_back(StNoFiNo);
		}
		else
		{
			int NoOfPart_iStIntersect = -1, NoOfPart_iFiIntersect = -1;
			int NoOfPartLeftFrom_iSt = -1, NoOfPartRightFrom_iFi = -1;
			vector<int, allocator<int> > NosOfPartsTotallyInside;

			char MergePartFromLeft = 0, MergePartFromRight = 0;
			int OldVectSize = (int)StNoFiNoVect.size();

			for(int i=0; i<OldVectSize; i++)
			{
				srTStNoFiNo StNoFiNo = StNoFiNoVect[i];

				if((iSt >= StNoFiNo.StNo) && (iFi <= StNoFiNo.FiNo)) 
				{
					GenStartOffset = StNoFiNo.StOffset + (iSt - StNoFiNo.StNo);
					goto SettingPointers;
				}

				if(StNoFiNo.StNo <= iSt) 
				{
					if(StNoFiNo.FiNo <= iSt) NoOfPartLeftFrom_iSt = i;
					else NoOfPart_iStIntersect = i;
				}
				if(StNoFiNo.StNo < iFi)
				{
					if(iFi <= StNoFiNo.FiNo) NoOfPart_iFiIntersect = i;
				}
				else if(NoOfPartRightFrom_iFi == -1) NoOfPartRightFrom_iFi = i;

				if((iSt <= StNoFiNo.StNo) && (iFi >= StNoFiNo.FiNo))
				{
					NosOfPartsTotallyInside.push_back(i);
				}
			}

			// Recent part
			int DelPartCount = 0;
			for(vector<int, allocator<int> >::iterator Iter_NosTotInside = NosOfPartsTotallyInside.begin(); Iter_NosTotInside != NosOfPartsTotallyInside.end(); ++Iter_NosTotInside)
			{
				int NoOfPartTotallyInside = *Iter_NosTotInside - DelPartCount;
				if(NoOfPartTotallyInside < (OldVectSize - 1))
				{
					srTStNoFiNo StNoFiNo_TotInside = StNoFiNoVect[NoOfPartTotallyInside];
					//int LocOffs = StNoFiNo_TotInside.StOffset;
					long long LocOffs = StNoFiNo_TotInside.StOffset;
					//int AmOfPo_TotInside = StNoFiNo_TotInside.FiNo - StNoFiNo_TotInside.StNo + 1;
					long long AmOfPo_TotInside = StNoFiNo_TotInside.FiNo - StNoFiNo_TotInside.StNo + 1;

					double *pBtx = BtxArrP[LevelNo], *pX = XArrP[LevelNo], *pIntBtxE2 = IntBtxE2ArrP[LevelNo], *pBx = BxArrP[LevelNo];
					double *pBtz = BtzArrP[LevelNo], *pZ = ZArrP[LevelNo], *pIntBtzE2 = IntBtzE2ArrP[LevelNo], *pBz = BzArrP[LevelNo];

					double *tBtx = pBtx+LocOffs, *tX = pX+LocOffs, *tIntBtxE2 = pIntBtxE2+LocOffs, *tBx = pBx+LocOffs;
					double *tBtz = pBtz+LocOffs, *tZ = pZ+LocOffs, *tIntBtzE2 = pIntBtzE2+LocOffs, *tBz = pBz+LocOffs;

					for(int k=(NoOfPartTotallyInside+1); k<OldVectSize; k++)
					{
						srTStNoFiNo& StNoFiNo_Next = StNoFiNoVect[k];
						//int OffsNext = StNoFiNo_Next.StOffset;
						long long OffsNext = StNoFiNo_Next.StOffset;
						double *tBtx_Nx = pBtx+OffsNext, *tX_Nx = pX+OffsNext, *tIntBtxE2_Nx = pIntBtxE2+OffsNext, *tBx_Nx = pBx+OffsNext;
						double *tBtz_Nx = pBtz+OffsNext, *tZ_Nx = pZ+OffsNext, *tIntBtzE2_Nx = pIntBtzE2+OffsNext, *tBz_Nx = pBz+OffsNext;

						//for(int m=StNoFiNo_Next.StNo; m<=StNoFiNo_Next.FiNo; m++)
						for(long long m=StNoFiNo_Next.StNo; m<=StNoFiNo_Next.FiNo; m++)
						{
							*(tBtx++) = *(tBtx_Nx++); *(tX++) = *(tX_Nx++); *(tIntBtxE2++) = *(tIntBtxE2_Nx++); *(tBx++) = *(tBx_Nx++);
							*(tBtz++) = *(tBtz_Nx++); *(tZ++) = *(tZ_Nx++), *(tIntBtzE2++) = *(tIntBtzE2_Nx++); *(tBz++) = *(tBz_Nx++);
						}

						StNoFiNo_Next.StOffset -= AmOfPo_TotInside;
					}
				}
				srTStNoFiNoVect::iterator Iter_TotInside = StNoFiNoVect.begin();
				for(int p=0; p<NoOfPartTotallyInside; p++) ++Iter_TotInside;

				StNoFiNoVect.erase(Iter_TotInside);
				DelPartCount++;
				OldVectSize--;

				if(OldVectSize == 0)
				{
					delete[] (BtxArrP[LevelNo]);
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

					srTStNoFiNo StNoFiNo(iSt, iFi, 0);
					StNoFiNoVect.push_back(StNoFiNo);

					GenStartOffset = 0;
					goto SettingPointers;
				}

				if(NoOfPart_iStIntersect == NoOfPartTotallyInside) NoOfPart_iStIntersect = -1;
				if(NoOfPart_iFiIntersect > -1)
				{
					if(NoOfPart_iFiIntersect == NoOfPartTotallyInside) NoOfPart_iFiIntersect = -1;
					else NoOfPart_iFiIntersect--;
				}
				if(NoOfPartRightFrom_iFi > -1) NoOfPartRightFrom_iFi--;
			}
			// End Recent part

			//int SizeToAllocate = iFi - iSt + 1;
			long long SizeToAllocate = iFi - iSt + 1;
			if(NoOfPart_iStIntersect > -1) 
			{
				srTStNoFiNo StNoFiNo = StNoFiNoVect[NoOfPart_iStIntersect];
				SizeToAllocate += iSt - StNoFiNo.StNo;

				//int DifInt = StNoFiNo.FiNo - iSt + 1;
				long long DifInt = StNoFiNo.FiNo - iSt + 1;
				sStart += Two_h*DifInt; 
				Np -= DifInt;
			}
			if(NoOfPart_iFiIntersect > -1)
			{
				srTStNoFiNo StNoFiNo = StNoFiNoVect[NoOfPart_iFiIntersect];
				SizeToAllocate += StNoFiNo.FiNo - iFi;

				//int DifInt = iFi - StNoFiNo.StNo + 1;
				long long DifInt = iFi - StNoFiNo.StNo + 1;
				sEnd -= Two_h*DifInt; 
				Np -= DifInt;
			}

			if(NoOfPartLeftFrom_iSt > -1)
			{
				for(int i=0; i<=NoOfPartLeftFrom_iSt; i++)
				{
					srTStNoFiNo StNoFiNo = StNoFiNoVect[i];
					if(iSt == StNoFiNo.FiNo)
					{
						SizeToAllocate += StNoFiNo.FiNo - StNoFiNo.StNo;
						sStart += Two_h; 
						Np--;
						MergePartFromLeft = 1;
					}
					else SizeToAllocate += StNoFiNo.FiNo - StNoFiNo.StNo + 1;
				}
			}
			if(NoOfPartRightFrom_iFi > -1)
			{
				for(int k = NoOfPartRightFrom_iFi; k < OldVectSize; k++)
				{
					srTStNoFiNo StNoFiNo = StNoFiNoVect[k];
					if(iFi == StNoFiNo.StNo)
					{
						SizeToAllocate += StNoFiNo.FiNo - StNoFiNo.StNo;
						sEnd -= Two_h; 
						Np--;
						MergePartFromRight = 1;
					}
					else SizeToAllocate += StNoFiNo.FiNo - StNoFiNo.StNo + 1;
				}
			}

			double* BasePtr = new double[SizeToAllocate*8];
			if(BasePtr == 0) 
			{ 
				//pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;
				return MEMORY_ALLOCATION_FAILURE;
			}

			double* BufBtxArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufXArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufIntBtxE2ArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufBxArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufBtzArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufZArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufIntBtzE2ArrP = BasePtr; BasePtr += SizeToAllocate; 
			double* BufBzArrP = BasePtr;

			double *pBtx = BtxArrP[LevelNo], *pX = XArrP[LevelNo], *pIntBtxE2 = IntBtxE2ArrP[LevelNo], *pBx = BxArrP[LevelNo];
			double *pBtz = BtzArrP[LevelNo], *pZ = ZArrP[LevelNo], *pIntBtzE2 = IntBtzE2ArrP[LevelNo], *pBz = BzArrP[LevelNo];

			double *TrBufBtxArrP = BufBtxArrP, *TrBufXArrP = BufXArrP, *TrBufIntBtxE2ArrP = BufIntBtxE2ArrP, *TrBufBxArrP = BufBxArrP;
			double *TrBufBtzArrP = BufBtzArrP, *TrBufZArrP = BufZArrP, *TrBufIntBtzE2ArrP = BufIntBtzE2ArrP, *TrBufBzArrP = BufBzArrP;

			if(NoOfPartLeftFrom_iSt > -1)
			{
				for(int ii = 0; ii <= NoOfPartLeftFrom_iSt; ii++)
				{
					srTStNoFiNo StNoFiNo = StNoFiNoVect[ii];
					//for(int jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
					for(long long jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
					{
						*(TrBufBtxArrP++) = *(pBtx++); *(TrBufXArrP++) = *(pX++); *(TrBufIntBtxE2ArrP++) = *(pIntBtxE2++); *(TrBufBxArrP++) = *(pBx++);
						*(TrBufBtzArrP++) = *(pBtz++); *(TrBufZArrP++) = *(pZ++); *(TrBufIntBtzE2ArrP++) = *(pIntBtzE2++); *(TrBufBzArrP++) = *(pBz++);
					}
				}
			}
			if(NoOfPart_iStIntersect > -1)
			{
				srTStNoFiNo StNoFiNo = StNoFiNoVect[NoOfPart_iStIntersect];

				//for(int jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
				for(long long jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
				{
					*(TrBufBtxArrP++) = *(pBtx++); *(TrBufXArrP++) = *(pX++); *(TrBufIntBtxE2ArrP++) = *(pIntBtxE2++); *(TrBufBxArrP++) = *(pBx++);
					*(TrBufBtzArrP++) = *(pBtz++); *(TrBufZArrP++) = *(pZ++); *(TrBufIntBtzE2ArrP++) = *(pIntBtzE2++); *(TrBufBzArrP++) = *(pBz++);
				}
			}

			TrjDatPtr->CompTotalTrjData(sStart, sEnd, Np, TrBufBtxArrP, TrBufBtzArrP, TrBufXArrP, TrBufZArrP, TrBufIntBtxE2ArrP, TrBufIntBtzE2ArrP, TrBufBxArrP, TrBufBzArrP);

			TrBufBtxArrP += Np; TrBufXArrP += Np; TrBufIntBtxE2ArrP += Np; TrBufBxArrP += Np;
			TrBufBtzArrP += Np; TrBufZArrP += Np; TrBufIntBtzE2ArrP += Np; TrBufBzArrP += Np;

			if(NoOfPart_iFiIntersect > -1)
			{
				srTStNoFiNo StNoFiNo = StNoFiNoVect[NoOfPart_iFiIntersect];

				pBtx = BtxArrP[LevelNo], pX = XArrP[LevelNo], pIntBtxE2 = IntBtxE2ArrP[LevelNo], pBx = BxArrP[LevelNo];
				pBtz = BtzArrP[LevelNo], pZ = ZArrP[LevelNo], pIntBtzE2 = IntBtzE2ArrP[LevelNo], pBz = BzArrP[LevelNo];
				//int Offset = StNoFiNo.StOffset;
				long long Offset = StNoFiNo.StOffset;
				pBtx += Offset; pX += Offset; pIntBtxE2 += Offset; pBx += Offset;
				pBtz += Offset; pZ += Offset; pIntBtzE2 += Offset; pBz += Offset;

				//for(int jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
				for(long long jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
				{
					*(TrBufBtxArrP++) = *(pBtx++); *(TrBufXArrP++) = *(pX++); *(TrBufIntBtxE2ArrP++) = *(pIntBtxE2++); *(TrBufBxArrP++) = *(pBx++);
					*(TrBufBtzArrP++) = *(pBtz++); *(TrBufZArrP++) = *(pZ++); *(TrBufIntBtzE2ArrP++) = *(pIntBtzE2++); *(TrBufBzArrP++) = *(pBz++);
				}
			}
			if(NoOfPartRightFrom_iFi > -1)
			{
				for(int ii = NoOfPartRightFrom_iFi; ii < OldVectSize; ii++)
				{
					srTStNoFiNo StNoFiNo = StNoFiNoVect[ii];
					//for(int jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
					for(long long jj = StNoFiNo.StNo; jj <= StNoFiNo.FiNo; jj++)
					{
						*(TrBufBtxArrP++) = *(pBtx++); *(TrBufXArrP++) = *(pX++); *(TrBufIntBtxE2ArrP++) = *(pIntBtxE2++); *(TrBufBxArrP++) = *(pBx++);
						*(TrBufBtzArrP++) = *(pBtz++); *(TrBufZArrP++) = *(pZ++); *(TrBufIntBtzE2ArrP++) = *(pIntBtzE2++); *(TrBufBzArrP++) = *(pBz++);
					}
				}
			}
			delete[] (BtxArrP[LevelNo]);
			BtxArrP[LevelNo] = BufBtxArrP;
			XArrP[LevelNo] = BufXArrP;
			IntBtxE2ArrP[LevelNo] = BufIntBtxE2ArrP;
			BxArrP[LevelNo] = BufBxArrP;
			BtzArrP[LevelNo] = BufBtzArrP;
			ZArrP[LevelNo] = BufZArrP;
			IntBtzE2ArrP[LevelNo] = BufIntBtzE2ArrP;
			BzArrP[LevelNo] = BufBzArrP;

			srTStNoFiNoVect LocStNoFiNoVect;
			//int LocStNo = iSt, LocFiNo = iFi, LocStOffset = -1;
			long long LocStNo = iSt, LocFiNo = iFi, LocStOffset = -1;

			if(NoOfPartLeftFrom_iSt > -1)
			{
				int ActNoOfPartLeft = MergePartFromLeft? (NoOfPartLeftFrom_iSt - 1) : NoOfPartLeftFrom_iSt;
				
				srTStNoFiNo StNoFiNoOldBuf;
				for(int ii = 0; ii <= ActNoOfPartLeft; ii++)
				{
					StNoFiNoOldBuf = StNoFiNoVect[ii];
					LocStNoFiNoVect.push_back(StNoFiNoOldBuf);
					LocStOffset = StNoFiNoOldBuf.StOffset;
				}
				if(MergePartFromLeft)
				{
					srTStNoFiNo StNoFiNoOld = StNoFiNoVect[NoOfPartLeftFrom_iSt];
					LocStNo = StNoFiNoOld.StNo;
					LocStOffset = StNoFiNoOld.StOffset;
				}
				else
				{
					LocStOffset += StNoFiNoOldBuf.FiNo - StNoFiNoOldBuf.StNo + 1;
				}
			}
			if(NoOfPart_iStIntersect > -1)
			{
				srTStNoFiNo StNoFiNoOld = StNoFiNoVect[NoOfPart_iStIntersect];
				LocStNo = StNoFiNoOld.StNo;
				LocStOffset = StNoFiNoOld.StOffset;
			}

			if(LocStOffset == -1) LocStOffset = 0;
			if(NoOfPart_iFiIntersect > -1) LocFiNo = StNoFiNoVect[NoOfPart_iFiIntersect].FiNo;
			if(MergePartFromRight) LocFiNo = StNoFiNoVect[NoOfPartRightFrom_iFi].FiNo;

			srTStNoFiNo StNoFiNoNew(LocStNo, LocFiNo, LocStOffset);
			LocStNoFiNoVect.push_back(StNoFiNoNew);

			GenStartOffset = LocStOffset;
			if(NoOfPart_iStIntersect > -1)
			{
				GenStartOffset += iSt - StNoFiNoNew.StNo;
			}

			if(NoOfPartRightFrom_iFi > -1)
			{
				int ActNoOfPartRight = MergePartFromRight? (NoOfPartRightFrom_iFi + 1) : NoOfPartRightFrom_iFi;
				//int CorrectedOffset = LocStOffset + (LocFiNo - LocStNo + 1);
				long long CorrectedOffset = LocStOffset + (LocFiNo - LocStNo + 1);
				for(int ii = ActNoOfPartRight; ii < OldVectSize; ii++)
				{
					srTStNoFiNo StNoFiNo = StNoFiNoVect[ii];
					StNoFiNo.StOffset = CorrectedOffset;
					LocStNoFiNoVect.push_back(StNoFiNo);
					CorrectedOffset += StNoFiNo.FiNo - StNoFiNo.StNo + 1;
				}
			}

			StNoFiNoVect.erase(StNoFiNoVect.begin(), StNoFiNoVect.end());  // Maybe not necessary
			StNoFiNoVect = LocStNoFiNoVect;

			if((StNoFiNoNew.StNo == 0) && (StNoFiNoNew.FiNo == TotalNpOnLevel)) 
			{
				AmOfPointsOnLevel[LevelNo] = TotalNpOnLevel;
				NumberOfLevelsFilled++;
			}
		}
	}
	else
	{
		if(NumberOfLevelsFilled > LevelNo) goto SettingPointers;

		if(!StNoFiNoVect.empty())
		{
			StNoFiNoVect.erase(StNoFiNoVect.begin(), StNoFiNoVect.end());
			delete[] (BtxArrP[LevelNo]);
		}

		double* BasePtr = new double[Np*8];
		if(BasePtr == 0) 
		{ 
			//pSend->ErrorMessage("SR::Error900"); return MEMORY_ALLOCATION_FAILURE;
			return MEMORY_ALLOCATION_FAILURE;
		}

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
		srTStNoFiNo StNoFiNo(0, Np-1);
		StNoFiNoVect.push_back(StNoFiNo);
		NumberOfLevelsFilled++;
	}

SettingPointers:

	double *pBtx = BtxArrP[LevelNo], *pX = XArrP[LevelNo], *pIntBtxE2 = IntBtxE2ArrP[LevelNo], *pBx = BxArrP[LevelNo];
	double *pBtz = BtzArrP[LevelNo], *pZ = ZArrP[LevelNo], *pIntBtzE2 = IntBtzE2ArrP[LevelNo], *pBz = BzArrP[LevelNo];

	if(GenStartOffset > 0)
	{
		pBtx += GenStartOffset; pX += GenStartOffset; pIntBtxE2 += GenStartOffset; pBx += GenStartOffset;
		pBtz += GenStartOffset; pZ += GenStartOffset; pIntBtzE2 += GenStartOffset; pBz += GenStartOffset;
	}

	**TrjPtrs = pBtx; *(TrjPtrs[1]) = pX; *(TrjPtrs[2]) = pIntBtxE2; *(TrjPtrs[3]) = pBx;
	*(TrjPtrs[4]) = pBtz; *(TrjPtrs[5]) = pZ; *(TrjPtrs[6]) = pIntBtzE2; *(TrjPtrs[7]) = pBz;
	return 0;
}

//*************************************************************************

void srTRadInt::AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;
	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;

	//TrjDatPtr->AnalizeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
	// CW 3,4 fails to make this call. So the function is here:
	if(!TrjDatPtr->HorFieldIsNotZero) FieldIsSymOverZ = 1;
	if(!TrjDatPtr->VerFieldIsNotZero) FieldIsSymOverX = 1;

	if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;

	char ObsIsSymOverX = 0, ObsIsSymOverZ = 0;
	if(FieldIsSymOverX && (DistrInfoDat.nx > 1))
	{
		double xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
		double xTol = xStep*0.01; // To steer
		char TrjAngIsSmall = (::fabs(TrjDatPtr->EbmDat.dxds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0)) < xTol);
		ObsIsSymOverX = TrjAngIsSmall && (::fabs(0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd) - TrjDatPtr->EbmDat.x0) < xTol);
	}
	if(FieldIsSymOverZ && (DistrInfoDat.nz > 1))
	{
		double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
		double zTol = zStep*0.01; // To steer
		char TrjAngIsSmall = (::fabs(TrjDatPtr->EbmDat.dzds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0)) < zTol);
		ObsIsSymOverZ = TrjAngIsSmall && (::fabs(0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd) - TrjDatPtr->EbmDat.z0) < zTol);
	}
	FinalResAreSymOverX = (FieldIsSymOverX && ObsIsSymOverX);
	FinalResAreSymOverZ = (FieldIsSymOverZ && ObsIsSymOverZ);
}

//*************************************************************************

void srTRadInt::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTSRWRadStructAccessData& Rad)
{
	//long PerX = DistrInfoDat.nLamb << 1;
	//long PerZ = PerX*DistrInfoDat.nx;
	long long PerX = DistrInfoDat.nLamb << 1;
	long long PerZ = PerX*DistrInfoDat.nx;

	char SymWithRespectToXax, SymWithRespectToZax;
	int HalfNz = DistrInfoDat.nz >> 1, Nz_mi_1 = DistrInfoDat.nz - 1;
	//int izStart = ((HalfNz << 1) == DistrInfoDat.nz)? HalfNz : (HalfNz + 1);
	int HalfNx = DistrInfoDat.nx >> 1, Nx_mi_1 = DistrInfoDat.nx - 1;
	//int ixStart = ((HalfNx << 1) == DistrInfoDat.nx)? HalfNx : (HalfNx + 1);
	int iz, ix;

	if(FinalResAreSymOverZ)
	{
		if(FinalResAreSymOverX)
		{
			SymWithRespectToXax = 0; SymWithRespectToZax = 1;
			for(iz=0; iz<HalfNz; iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				for(ix=0; ix<HalfNx; ix++)
				{
					//long Offset = izPerZ + ix*PerX;
					long long Offset = izPerZ + ix*PerX;
					float* pOrigDataEx = Rad.pBaseRadX + Offset;
					float* pOrigDataEz = Rad.pBaseRadZ + Offset;
					//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					float* pSymDataEx = Rad.pBaseRadX + OffsetSym;
					float* pSymDataEz = Rad.pBaseRadZ + OffsetSym;
					CopySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
		SymWithRespectToXax = 1; SymWithRespectToZax = 0;
		for(iz=0; iz<HalfNz; iz++)
		{
			//long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
			long long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
			for(ix=0; ix<DistrInfoDat.nx; ix++)
			{
				//long ixPerX = ix*PerX;
				//long Offset = izPerZ + ixPerX;
				long long ixPerX = ix*PerX;
				long long Offset = izPerZ + ixPerX;
				float* pOrigDataEx = Rad.pBaseRadX + Offset;
				float* pOrigDataEz = Rad.pBaseRadZ + Offset;
				//long OffsetSym = BufZ + ixPerX;
				long long OffsetSym = BufZ + ixPerX;
				float* pSymDataEx = Rad.pBaseRadX + OffsetSym;
				float* pSymDataEz = Rad.pBaseRadZ + OffsetSym;
				CopySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
			}
		}
	}
	else if(FinalResAreSymOverX)
	{
		SymWithRespectToXax = 0; SymWithRespectToZax = 1;
		for(iz=0; iz<DistrInfoDat.nz; iz++)
		{
			//long izPerZ = iz*PerZ;
			long long izPerZ = iz*PerZ;
			for(ix=0; ix<HalfNx; ix++)
			{
				//long Offset = izPerZ + ix*PerX;
				long long Offset = izPerZ + ix*PerX;
				float* pOrigDataEx = Rad.pBaseRadX + Offset;
				float* pOrigDataEz = Rad.pBaseRadZ + Offset;
				//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				float* pSymDataEx = Rad.pBaseRadX + OffsetSym;
				float* pSymDataEz = Rad.pBaseRadZ + OffsetSym;
				CopySymEnergySlice(pOrigDataEx, pOrigDataEz, pSymDataEx, pSymDataEz, SymWithRespectToXax, SymWithRespectToZax);
			}
		}
	}
}

//*************************************************************************

int srTRadInt::CheckFurtherSubdNeedForOneCoord(srTEFourier** EwPtrsArr)
{
	srTEFourier* pEw = *EwPtrsArr;
	double EwXRePr = pEw->EwX_Re;
	double EwXImPr = pEw->EwX_Im;
	double EwZRePr = pEw->EwZ_Re;
	double EwZImPr = pEw->EwZ_Im;

	char XReGoesUp=-1, XImGoesUp=-1, ZReGoesUp=-1, ZImGoesUp=-1;
	char XReGoesUpPr, XImGoesUpPr, ZReGoesUpPr, ZImGoesUpPr;
	char XReOnePeakFound=0, XImOnePeakFound=0, ZReOnePeakFound=0, ZImOnePeakFound=0;
	char XReNeedsSubd=0, XImNeedsSubd=0, ZReNeedsSubd=0, ZImNeedsSubd=0;

	for(int i=0; i<4; i++)
	{
		++pEw;
		XReGoesUpPr = XReGoesUp;
		XReGoesUp = (pEw->EwX_Re > EwXRePr)? 1 : 0;
		XImGoesUpPr = XImGoesUp;
		XImGoesUp = (pEw->EwX_Im > EwXImPr)? 1 : 0;
		ZReGoesUpPr = ZReGoesUp;
		ZReGoesUp = (pEw->EwZ_Re > EwZRePr)? 1 : 0;
		ZImGoesUpPr = ZImGoesUp;
		ZImGoesUp = (pEw->EwZ_Im > EwZImPr)? 1 : 0;
		if(i != 0)
		{
			if((XReGoesUpPr && (!XReGoesUp)) || ((!XReGoesUpPr) && XReGoesUp))
			{
				if(XReOnePeakFound) XReNeedsSubd = 1;
				else XReOnePeakFound = 1;
			}
			if((XImGoesUpPr && (!XImGoesUp)) || ((!XImGoesUpPr) && XImGoesUp))
			{
				if(XImOnePeakFound) XImNeedsSubd = 1;
				else XImOnePeakFound = 1;
			}
			if((ZReGoesUpPr && (!ZReGoesUp)) || ((!ZReGoesUpPr) && ZReGoesUp))
			{
				if(ZReOnePeakFound) ZReNeedsSubd = 1;
				else ZReOnePeakFound = 1;
			}
			if((ZImGoesUpPr && (!ZImGoesUp)) || ((!ZImGoesUpPr) && ZImGoesUp))
			{
				if(ZImOnePeakFound) ZImNeedsSubd = 1;
				else ZImOnePeakFound = 1;
			}
		}
		EwXRePr = pEw->EwX_Re;
		EwXImPr = pEw->EwX_Im;
		EwZRePr = pEw->EwZ_Re;
		EwZImPr = pEw->EwZ_Im;
	}
	return (XReNeedsSubd || XImNeedsSubd || ZReNeedsSubd || ZImNeedsSubd);
}

//*************************************************************************

int srTRadInt::RadInterpolationOnePointXZ(srTEFourierVect* pEwVect, int ixOffset, int izOffset, double xStep, double zStep, srTEFourier* pEw)
{
	double xStepInv = 1./xStep;
	double zStepInv = 1./zStep;
	int nx_mi_1 = int((DistrInfoDat.xEnd - DistrInfoDat.xStart)*xStepInv + 1.E-06);
	int nz_mi_1 = int((DistrInfoDat.zEnd - DistrInfoDat.zStart)*zStepInv + 1.E-06);
	int nx = nx_mi_1 + 1; //, nz = nz_mi_1 + 1;

	int ixc = int((ObsCoor.x - DistrInfoDat.xStart)*xStepInv + 1.E-06);
	int izc = int((ObsCoor.z - DistrInfoDat.zStart)*zStepInv + 1.E-06);
	int ixSt, izSt;
	if(ixc == nx_mi_1) ixSt = ixc - 3;
	else ixSt = ixc - ixOffset;
	if(izc == nz_mi_1) izSt = izc - 3;
	else izSt = izc - izOffset;

	double azE2 = zStepInv*zStepInv;
	double azE3 = azE2*zStepInv;
	double axaz = xStepInv*zStepInv;
	double axazE2 = axaz*zStepInv;
	double axazE3 = axazE2*zStepInv;
	double axE2 = xStepInv*xStepInv;
	double axE2az = axE2*zStepInv;
	double axE2azE2 = axE2az*zStepInv;
	double axE2azE3 = axE2azE2*zStepInv;
	double axE3 = axE2*xStepInv;
	double axE3az = axE3*zStepInv;
	double axE3azE2 = axE3az*zStepInv;
	double axE3azE3 = axE3azE2*zStepInv;

	double cAx0z1 = 0.1666666667*zStepInv;
	double cAx0z2 = 0.5*azE2;
	double cAx0z3 = 0.1666666667*azE3;
	double cAx1z0 = 0.1666666667*xStepInv;
	double cAx1z1 = 0.027777777778*axaz;
	double cAx1z2 = 0.083333333333*axazE2;
	double cAx1z3 = 0.027777777778*axazE3;
	double cAx2z0 = 0.5*axE2;
	double cAx2z1 = 0.083333333333*axE2az;
	double cAx2z2 = 0.25*axE2azE2;
	double cAx2z3 = 0.083333333333*axE2azE3;
	double cAx3z0 = 0.16666666667*axE3;
	double cAx3z1 = 0.027777777778*axE3az;
	double cAx3z2 = 0.083333333333*axE3azE2;
	double cAx3z3 = 0.027777777778*axE3azE3;

	double x = ObsCoor.x - (DistrInfoDat.xStart + (ixSt + 1)*xStep);
	double z = ObsCoor.z - (DistrInfoDat.zStart + (izSt + 1)*zStep);
	double xE2 = x*x, xz = x*z, zE2 = z*z;
	double xE3 = xE2*x, xE2z = xE2*z, xzE2 = x*zE2, zE3 = zE2*z, xE2zE2 = xE2*zE2;
	double xE3z = xE3*z, xE3zE2 = xE3*zE2, xE3zE3 = xE3*zE3, xE2zE3 = xE2*zE3, xzE3 = x*zE3;

	//long AbsIndSt = izSt*nx + ixSt;
	long long AbsIndSt = izSt*nx + ixSt;
	srTEFourierVect::iterator ItBegin = pEwVect->begin();
	srTEFourierVect::iterator It = ItBegin + AbsIndSt;
	srTEFourier &f00 = *(It++), &f10 = *(It++), &f20 = *(It++), &f30 = *(It++);
	AbsIndSt += nx;
	It = ItBegin + AbsIndSt;
	srTEFourier &f01 = *(It++), &f11 = *(It++), &f21 = *(It++), &f31 = *(It++);
	AbsIndSt += nx;
	It = ItBegin + AbsIndSt;
	srTEFourier &f02 = *(It++), &f12 = *(It++), &f22 = *(It++), &f32 = *(It++);
	AbsIndSt += nx;
	It = ItBegin + AbsIndSt;
	srTEFourier &f03 = *(It++), &f13 = *(It++), &f23 = *(It++), &f33 = *(It++);

	double Ax0z0 = f11.EwX_Re;
	double Ax0z1 = (-2*f10.EwX_Re - 3*f11.EwX_Re + 6*f12.EwX_Re - f13.EwX_Re)*cAx0z1;
	double Ax0z2 = (f10.EwX_Re + f12.EwX_Re - 2*f11.EwX_Re)*cAx0z2;
	double Ax0z3 = (f13.EwX_Re - f10.EwX_Re + 3*(f11.EwX_Re - f12.EwX_Re))*cAx0z3;
	double Ax1z0 = (-2*f01.EwX_Re - 3*f11.EwX_Re + 6*f21.EwX_Re - f31.EwX_Re)*cAx1z0;
	double Ax1z1 = (4*f00.EwX_Re + 6*(f01.EwX_Re + f10.EwX_Re - f23.EwX_Re - f32.EwX_Re) - 12*(f02.EwX_Re + f20.EwX_Re) + 2*(f03.EwX_Re + f30.EwX_Re) + 9*f11.EwX_Re - 18*(f12.EwX_Re + f21.EwX_Re) + 3*(f13.EwX_Re + f31.EwX_Re) + 36*f22.EwX_Re + f33.EwX_Re)*cAx1z1;
	double Ax1z2 = (-2*(f00.EwX_Re + f02.EwX_Re - f31.EwX_Re) + 4*f01.EwX_Re - 3*(f10.EwX_Re + f12.EwX_Re) + 6*(f11.EwX_Re + f20.EwX_Re + f22.EwX_Re) - 12*f21.EwX_Re - f30.EwX_Re - f32.EwX_Re)*cAx1z2;
	double Ax1z3 = (2*(f00.EwX_Re - f03.EwX_Re) + 6*(-f01.EwX_Re + f02.EwX_Re - f20.EwX_Re + f23.EwX_Re) + 3*(f10.EwX_Re - f13.EwX_Re - f31.EwX_Re + f32.EwX_Re) + 9*(f12.EwX_Re - f11.EwX_Re) + 18*(f21.EwX_Re - f22.EwX_Re) + f30.EwX_Re - f33.EwX_Re)*cAx1z3;
	double Ax2z0 = (f01.EwX_Re + f21.EwX_Re - 2*f11.EwX_Re)*cAx2z0;
	double Ax2z1 = (2*(-f00.EwX_Re + f13.EwX_Re - f20.EwX_Re) - 3*(f21.EwX_Re + f01.EwX_Re) + 6*(f02.EwX_Re + f11.EwX_Re + f22.EwX_Re) + 4*f10.EwX_Re - 12*f12.EwX_Re - f23.EwX_Re - f03.EwX_Re)*cAx2z1;
	double Ax2z2 = (f00.EwX_Re + f02.EwX_Re + f22.EwX_Re + f20.EwX_Re - 2*(f01.EwX_Re + f10.EwX_Re + f12.EwX_Re + f21.EwX_Re) + 4*f11.EwX_Re)*cAx2z2;
	double Ax2z3 = (-f00.EwX_Re + f03.EwX_Re - f20.EwX_Re + f23.EwX_Re + 3*(f01.EwX_Re - f02.EwX_Re + f21.EwX_Re - f22.EwX_Re) + 2*(f10.EwX_Re - f13.EwX_Re) + 6*(f12.EwX_Re - f11.EwX_Re))*cAx2z3;
	double Ax3z0 = (f31.EwX_Re - f01.EwX_Re + 3*(f11.EwX_Re - f21.EwX_Re))*cAx3z0;
	double Ax3z1 = (2*(f00.EwX_Re - f30.EwX_Re) + 3*(f01.EwX_Re - f13.EwX_Re + f23.EwX_Re - f31.EwX_Re) + 6*(-f02.EwX_Re - f10.EwX_Re + f20.EwX_Re + f32.EwX_Re) + 9*(f21.EwX_Re - f11.EwX_Re) + 18*(f12.EwX_Re - f22.EwX_Re) + f03.EwX_Re - f33.EwX_Re)*cAx3z1;
	double Ax3z2 = (f30.EwX_Re + f32.EwX_Re - f00.EwX_Re - f02.EwX_Re + 2*(f01.EwX_Re - f31.EwX_Re) + 3*(f10.EwX_Re + f12.EwX_Re - f20.EwX_Re - f22.EwX_Re) + 6*(f21.EwX_Re - f11.EwX_Re))*cAx3z2;
	double Ax3z3 = (f00.EwX_Re - f03.EwX_Re - f30.EwX_Re + f33.EwX_Re + 3*(-f01.EwX_Re + f02.EwX_Re - f10.EwX_Re + f13.EwX_Re + f20.EwX_Re - f23.EwX_Re + f31.EwX_Re - f32.EwX_Re) + 9*(f11.EwX_Re - f12.EwX_Re - f21.EwX_Re + f22.EwX_Re))*cAx3z3;
	pEw->EwX_Re = Ax3z3*xE3zE3 + Ax3z2*xE3zE2 + Ax3z1*xE3z + Ax3z0*xE3 + Ax2z3*xE2zE3 + Ax2z2*xE2zE2 + Ax2z1*xE2z + Ax2z0*xE2 + Ax1z3*xzE3 + Ax1z2*xzE2 + Ax1z1*xz + Ax1z0*x + Ax0z3*zE3 + Ax0z2*zE2 + Ax0z1*z + Ax0z0;

	Ax0z0 = f11.EwX_Im;
	Ax0z1 = (-2*f10.EwX_Im - 3*f11.EwX_Im + 6*f12.EwX_Im - f13.EwX_Im)*cAx0z1;
	Ax0z2 = (f10.EwX_Im + f12.EwX_Im - 2*f11.EwX_Im)*cAx0z2;
	Ax0z3 = (f13.EwX_Im - f10.EwX_Im + 3*(f11.EwX_Im - f12.EwX_Im))*cAx0z3;
	Ax1z0 = (-2*f01.EwX_Im - 3*f11.EwX_Im + 6*f21.EwX_Im - f31.EwX_Im)*cAx1z0;
	Ax1z1 = (4*f00.EwX_Im + 6*(f01.EwX_Im + f10.EwX_Im - f23.EwX_Im - f32.EwX_Im) - 12*(f02.EwX_Im + f20.EwX_Im) + 2*(f03.EwX_Im + f30.EwX_Im) + 9*f11.EwX_Im - 18*(f12.EwX_Im + f21.EwX_Im) + 3*(f13.EwX_Im + f31.EwX_Im) + 36*f22.EwX_Im + f33.EwX_Im)*cAx1z1;
	Ax1z2 = (-2*(f00.EwX_Im + f02.EwX_Im - f31.EwX_Im) + 4*f01.EwX_Im - 3*(f10.EwX_Im + f12.EwX_Im) + 6*(f11.EwX_Im + f20.EwX_Im + f22.EwX_Im) - 12*f21.EwX_Im - f30.EwX_Im - f32.EwX_Im)*cAx1z2;
	Ax1z3 = (2*(f00.EwX_Im - f03.EwX_Im) + 6*(-f01.EwX_Im + f02.EwX_Im - f20.EwX_Im + f23.EwX_Im) + 3*(f10.EwX_Im - f13.EwX_Im - f31.EwX_Im + f32.EwX_Im) + 9*(f12.EwX_Im - f11.EwX_Im) + 18*(f21.EwX_Im - f22.EwX_Im) + f30.EwX_Im - f33.EwX_Im)*cAx1z3;
	Ax2z0 = (f01.EwX_Im + f21.EwX_Im - 2*f11.EwX_Im)*cAx2z0;
	Ax2z1 = (2*(-f00.EwX_Im + f13.EwX_Im - f20.EwX_Im) - 3*(f21.EwX_Im + f01.EwX_Im) + 6*(f02.EwX_Im + f11.EwX_Im + f22.EwX_Im) + 4*f10.EwX_Im - 12*f12.EwX_Im - f23.EwX_Im - f03.EwX_Im)*cAx2z1;
	Ax2z2 = (f00.EwX_Im + f02.EwX_Im + f22.EwX_Im + f20.EwX_Im - 2*(f01.EwX_Im + f10.EwX_Im + f12.EwX_Im + f21.EwX_Im) + 4*f11.EwX_Im)*cAx2z2;
	Ax2z3 = (-f00.EwX_Im + f03.EwX_Im - f20.EwX_Im + f23.EwX_Im + 3*(f01.EwX_Im - f02.EwX_Im + f21.EwX_Im - f22.EwX_Im) + 2*(f10.EwX_Im - f13.EwX_Im) + 6*(f12.EwX_Im - f11.EwX_Im))*cAx2z3;
	Ax3z0 = (f31.EwX_Im - f01.EwX_Im + 3*(f11.EwX_Im - f21.EwX_Im))*cAx3z0;
	Ax3z1 = (2*(f00.EwX_Im - f30.EwX_Im) + 3*(f01.EwX_Im - f13.EwX_Im + f23.EwX_Im - f31.EwX_Im) + 6*(-f02.EwX_Im - f10.EwX_Im + f20.EwX_Im + f32.EwX_Im) + 9*(f21.EwX_Im - f11.EwX_Im) + 18*(f12.EwX_Im - f22.EwX_Im) + f03.EwX_Im - f33.EwX_Im)*cAx3z1;
	Ax3z2 = (f30.EwX_Im + f32.EwX_Im - f00.EwX_Im - f02.EwX_Im + 2*(f01.EwX_Im - f31.EwX_Im) + 3*(f10.EwX_Im + f12.EwX_Im - f20.EwX_Im - f22.EwX_Im) + 6*(f21.EwX_Im - f11.EwX_Im))*cAx3z2;
	Ax3z3 = (f00.EwX_Im - f03.EwX_Im - f30.EwX_Im + f33.EwX_Im + 3*(-f01.EwX_Im + f02.EwX_Im - f10.EwX_Im + f13.EwX_Im + f20.EwX_Im - f23.EwX_Im + f31.EwX_Im - f32.EwX_Im) + 9*(f11.EwX_Im - f12.EwX_Im - f21.EwX_Im + f22.EwX_Im))*cAx3z3;
	pEw->EwX_Im = Ax3z3*xE3zE3 + Ax3z2*xE3zE2 + Ax3z1*xE3z + Ax3z0*xE3 + Ax2z3*xE2zE3 + Ax2z2*xE2zE2 + Ax2z1*xE2z + Ax2z0*xE2 + Ax1z3*xzE3 + Ax1z2*xzE2 + Ax1z1*xz + Ax1z0*x + Ax0z3*zE3 + Ax0z2*zE2 + Ax0z1*z + Ax0z0;

	Ax0z0 = f11.EwZ_Re;
	Ax0z1 = (-2*f10.EwZ_Re - 3*f11.EwZ_Re + 6*f12.EwZ_Re - f13.EwZ_Re)*cAx0z1;
	Ax0z2 = (f10.EwZ_Re + f12.EwZ_Re - 2*f11.EwZ_Re)*cAx0z2;
	Ax0z3 = (f13.EwZ_Re - f10.EwZ_Re + 3*(f11.EwZ_Re - f12.EwZ_Re))*cAx0z3;
	Ax1z0 = (-2*f01.EwZ_Re - 3*f11.EwZ_Re + 6*f21.EwZ_Re - f31.EwZ_Re)*cAx1z0;
	Ax1z1 = (4*f00.EwZ_Re + 6*(f01.EwZ_Re + f10.EwZ_Re - f23.EwZ_Re - f32.EwZ_Re) - 12*(f02.EwZ_Re + f20.EwZ_Re) + 2*(f03.EwZ_Re + f30.EwZ_Re) + 9*f11.EwZ_Re - 18*(f12.EwZ_Re + f21.EwZ_Re) + 3*(f13.EwZ_Re + f31.EwZ_Re) + 36*f22.EwZ_Re + f33.EwZ_Re)*cAx1z1;
	Ax1z2 = (-2*(f00.EwZ_Re + f02.EwZ_Re - f31.EwZ_Re) + 4*f01.EwZ_Re - 3*(f10.EwZ_Re + f12.EwZ_Re) + 6*(f11.EwZ_Re + f20.EwZ_Re + f22.EwZ_Re) - 12*f21.EwZ_Re - f30.EwZ_Re - f32.EwZ_Re)*cAx1z2;
	Ax1z3 = (2*(f00.EwZ_Re - f03.EwZ_Re) + 6*(-f01.EwZ_Re + f02.EwZ_Re - f20.EwZ_Re + f23.EwZ_Re) + 3*(f10.EwZ_Re - f13.EwZ_Re - f31.EwZ_Re + f32.EwZ_Re) + 9*(f12.EwZ_Re - f11.EwZ_Re) + 18*(f21.EwZ_Re - f22.EwZ_Re) + f30.EwZ_Re - f33.EwZ_Re)*cAx1z3;
	Ax2z0 = (f01.EwZ_Re + f21.EwZ_Re - 2*f11.EwZ_Re)*cAx2z0;
	Ax2z1 = (2*(-f00.EwZ_Re + f13.EwZ_Re - f20.EwZ_Re) - 3*(f21.EwZ_Re + f01.EwZ_Re) + 6*(f02.EwZ_Re + f11.EwZ_Re + f22.EwZ_Re) + 4*f10.EwZ_Re - 12*f12.EwZ_Re - f23.EwZ_Re - f03.EwZ_Re)*cAx2z1;
	Ax2z2 = (f00.EwZ_Re + f02.EwZ_Re + f22.EwZ_Re + f20.EwZ_Re - 2*(f01.EwZ_Re + f10.EwZ_Re + f12.EwZ_Re + f21.EwZ_Re) + 4*f11.EwZ_Re)*cAx2z2;
	Ax2z3 = (-f00.EwZ_Re + f03.EwZ_Re - f20.EwZ_Re + f23.EwZ_Re + 3*(f01.EwZ_Re - f02.EwZ_Re + f21.EwZ_Re - f22.EwZ_Re) + 2*(f10.EwZ_Re - f13.EwZ_Re) + 6*(f12.EwZ_Re - f11.EwZ_Re))*cAx2z3;
	Ax3z0 = (f31.EwZ_Re - f01.EwZ_Re + 3*(f11.EwZ_Re - f21.EwZ_Re))*cAx3z0;
	Ax3z1 = (2*(f00.EwZ_Re - f30.EwZ_Re) + 3*(f01.EwZ_Re - f13.EwZ_Re + f23.EwZ_Re - f31.EwZ_Re) + 6*(-f02.EwZ_Re - f10.EwZ_Re + f20.EwZ_Re + f32.EwZ_Re) + 9*(f21.EwZ_Re - f11.EwZ_Re) + 18*(f12.EwZ_Re - f22.EwZ_Re) + f03.EwZ_Re - f33.EwZ_Re)*cAx3z1;
	Ax3z2 = (f30.EwZ_Re + f32.EwZ_Re - f00.EwZ_Re - f02.EwZ_Re + 2*(f01.EwZ_Re - f31.EwZ_Re) + 3*(f10.EwZ_Re + f12.EwZ_Re - f20.EwZ_Re - f22.EwZ_Re) + 6*(f21.EwZ_Re - f11.EwZ_Re))*cAx3z2;
	Ax3z3 = (f00.EwZ_Re - f03.EwZ_Re - f30.EwZ_Re + f33.EwZ_Re + 3*(-f01.EwZ_Re + f02.EwZ_Re - f10.EwZ_Re + f13.EwZ_Re + f20.EwZ_Re - f23.EwZ_Re + f31.EwZ_Re - f32.EwZ_Re) + 9*(f11.EwZ_Re - f12.EwZ_Re - f21.EwZ_Re + f22.EwZ_Re))*cAx3z3;
	pEw->EwZ_Re = Ax3z3*xE3zE3 + Ax3z2*xE3zE2 + Ax3z1*xE3z + Ax3z0*xE3 + Ax2z3*xE2zE3 + Ax2z2*xE2zE2 + Ax2z1*xE2z + Ax2z0*xE2 + Ax1z3*xzE3 + Ax1z2*xzE2 + Ax1z1*xz + Ax1z0*x + Ax0z3*zE3 + Ax0z2*zE2 + Ax0z1*z + Ax0z0;

	Ax0z0 = f11.EwZ_Im;
	Ax0z1 = (-2*f10.EwZ_Im - 3*f11.EwZ_Im + 6*f12.EwZ_Im - f13.EwZ_Im)*cAx0z1;
	Ax0z2 = (f10.EwZ_Im + f12.EwZ_Im - 2*f11.EwZ_Im)*cAx0z2;
	Ax0z3 = (f13.EwZ_Im - f10.EwZ_Im + 3*(f11.EwZ_Im - f12.EwZ_Im))*cAx0z3;
	Ax1z0 = (-2*f01.EwZ_Im - 3*f11.EwZ_Im + 6*f21.EwZ_Im - f31.EwZ_Im)*cAx1z0;
	Ax1z1 = (4*f00.EwZ_Im + 6*(f01.EwZ_Im + f10.EwZ_Im - f23.EwZ_Im - f32.EwZ_Im) - 12*(f02.EwZ_Im + f20.EwZ_Im) + 2*(f03.EwZ_Im + f30.EwZ_Im) + 9*f11.EwZ_Im - 18*(f12.EwZ_Im + f21.EwZ_Im) + 3*(f13.EwZ_Im + f31.EwZ_Im) + 36*f22.EwZ_Im + f33.EwZ_Im)*cAx1z1;
	Ax1z2 = (-2*(f00.EwZ_Im + f02.EwZ_Im - f31.EwZ_Im) + 4*f01.EwZ_Im - 3*(f10.EwZ_Im + f12.EwZ_Im) + 6*(f11.EwZ_Im + f20.EwZ_Im + f22.EwZ_Im) - 12*f21.EwZ_Im - f30.EwZ_Im - f32.EwZ_Im)*cAx1z2;
	Ax1z3 = (2*(f00.EwZ_Im - f03.EwZ_Im) + 6*(-f01.EwZ_Im + f02.EwZ_Im - f20.EwZ_Im + f23.EwZ_Im) + 3*(f10.EwZ_Im - f13.EwZ_Im - f31.EwZ_Im + f32.EwZ_Im) + 9*(f12.EwZ_Im - f11.EwZ_Im) + 18*(f21.EwZ_Im - f22.EwZ_Im) + f30.EwZ_Im - f33.EwZ_Im)*cAx1z3;
	Ax2z0 = (f01.EwZ_Im + f21.EwZ_Im - 2*f11.EwZ_Im)*cAx2z0;
	Ax2z1 = (2*(-f00.EwZ_Im + f13.EwZ_Im - f20.EwZ_Im) - 3*(f21.EwZ_Im + f01.EwZ_Im) + 6*(f02.EwZ_Im + f11.EwZ_Im + f22.EwZ_Im) + 4*f10.EwZ_Im - 12*f12.EwZ_Im - f23.EwZ_Im - f03.EwZ_Im)*cAx2z1;
	Ax2z2 = (f00.EwZ_Im + f02.EwZ_Im + f22.EwZ_Im + f20.EwZ_Im - 2*(f01.EwZ_Im + f10.EwZ_Im + f12.EwZ_Im + f21.EwZ_Im) + 4*f11.EwZ_Im)*cAx2z2;
	Ax2z3 = (-f00.EwZ_Im + f03.EwZ_Im - f20.EwZ_Im + f23.EwZ_Im + 3*(f01.EwZ_Im - f02.EwZ_Im + f21.EwZ_Im - f22.EwZ_Im) + 2*(f10.EwZ_Im - f13.EwZ_Im) + 6*(f12.EwZ_Im - f11.EwZ_Im))*cAx2z3;
	Ax3z0 = (f31.EwZ_Im - f01.EwZ_Im + 3*(f11.EwZ_Im - f21.EwZ_Im))*cAx3z0;
	Ax3z1 = (2*(f00.EwZ_Im - f30.EwZ_Im) + 3*(f01.EwZ_Im - f13.EwZ_Im + f23.EwZ_Im - f31.EwZ_Im) + 6*(-f02.EwZ_Im - f10.EwZ_Im + f20.EwZ_Im + f32.EwZ_Im) + 9*(f21.EwZ_Im - f11.EwZ_Im) + 18*(f12.EwZ_Im - f22.EwZ_Im) + f03.EwZ_Im - f33.EwZ_Im)*cAx3z1;
	Ax3z2 = (f30.EwZ_Im + f32.EwZ_Im - f00.EwZ_Im - f02.EwZ_Im + 2*(f01.EwZ_Im - f31.EwZ_Im) + 3*(f10.EwZ_Im + f12.EwZ_Im - f20.EwZ_Im - f22.EwZ_Im) + 6*(f21.EwZ_Im - f11.EwZ_Im))*cAx3z2;
	Ax3z3 = (f00.EwZ_Im - f03.EwZ_Im - f30.EwZ_Im + f33.EwZ_Im + 3*(-f01.EwZ_Im + f02.EwZ_Im - f10.EwZ_Im + f13.EwZ_Im + f20.EwZ_Im - f23.EwZ_Im + f31.EwZ_Im - f32.EwZ_Im) + 9*(f11.EwZ_Im - f12.EwZ_Im - f21.EwZ_Im + f22.EwZ_Im))*cAx3z3;
	pEw->EwZ_Im = Ax3z3*xE3zE3 + Ax3z2*xE3zE2 + Ax3z1*xE3z + Ax3z0*xE3 + Ax2z3*xE2zE3 + Ax2z2*xE2zE2 + Ax2z1*xE2z + Ax2z0*xE2 + Ax1z3*xzE3 + Ax1z2*xzE2 + Ax1z1*xz + Ax1z0*x + Ax0z3*zE3 + Ax0z2*zE2 + Ax0z1*z + Ax0z0;

	return 0;
}

//*************************************************************************

void srTRadInt::ComputeElectricFieldFreqDomain(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecElecFld* pPrecElecFld, srTSRWRadStructAccessData* pWfr, char showProgressInd)
{
	if((pTrjDat == 0) || (pWfrSmp == 0) || (pPrecElecFld == 0) || (pWfr == 0)) throw INCORRECT_PARAMS_SR_COMP;
	int res = 0;

	Initialize();

	DistrInfoDat = *pWfrSmp;
	DistrInfoDat.EnsureZeroTransverseRangesForSinglePoints();
	pWfr->SetRadSamplingFromObs(DistrInfoDat);

	TrjDataContShouldBeRebuild = 1;
	SetInputTrjData(pTrjDat);

	SetPrecParams(pPrecElecFld);

	if(res = CheckInputConsistency()) throw res;
	if(res = ComputeTotalRadDistrDirectOut(*pWfr, showProgressInd)) throw res;

	srTGenOptElem GenOptElem;
	if(res = GenOptElem.ComputeRadMoments(pWfr)) throw res;

	//setting the average photon energy:
	pWfr->SetAvgPhotEnergyFromLimits();
}

//*************************************************************************

void srTRadInt::SetPrecParams(srTParPrecElecFld* pPrecElecFld)
{
	if(pPrecElecFld == 0) return;

	int MethNo = pPrecElecFld->IntegMethNo;
	if(MethNo == 0)
	{
        sIntegStep = sIntegStep_Input = pPrecElecFld->RelPrecOrStep;
        sIntegMethod = 1;
	}
	else if(MethNo == 1)
	{
        sIntegMethod = 10;
        sIntegRelPrec = pPrecElecFld->RelPrecOrStep;
	}
	else if(MethNo == 2)
	{
        sIntegMethod = 11;
        sIntegRelPrec = pPrecElecFld->RelPrecOrStep;
	}

	sIntegStart = TrjDatPtr->sStart;
	sIntegFin = TrjDatPtr->sStart + TrjDatPtr->sStep*(TrjDatPtr->LenFieldData - 1);

	double sStartIntPrec = pPrecElecFld->sStartInt;
	double sEndIntPrec = pPrecElecFld->sEndInt;
	bool SpecLimitsMayBeDefined = (sStartIntPrec < sEndIntPrec);
	if(SpecLimitsMayBeDefined && (sStartIntPrec > sIntegStart) && (sStartIntPrec < sIntegFin)) sIntegStart = sStartIntPrec;
	if(SpecLimitsMayBeDefined && (sEndIntPrec > sIntegStart) && (sEndIntPrec < sIntegFin)) sIntegFin = sEndIntPrec;

	MaxNumPoToSave = 10000;
	TryToApplyNearFieldResidual = 1;
	TryToApplyNearFieldResidual_AtRight = 0; // because it's buggy

	m_CalcResidTerminTerms = pPrecElecFld->CalcTerminTerms;
}

//*************************************************************************

