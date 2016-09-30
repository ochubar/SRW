/************************************************************************//**
 * File: srradinc.cpp
 * Description: SR calculation from Constant Magnetic Field
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srradinc.h"
#include "srprgind.h"
#include "gmfft.h"
//#include "srmamet.h"
#include "gmmeth.h"
#include "gmfunc.h"

//*************************************************************************

int srTRadIntConst::CheckInputConsistency()
{
// Checking El. Beam
	srTEbmDat& Ebm = ConstTrjDatPtr->EbmDat;
	char ThinElBeamIsNotSetUp = (Ebm.Energy <= 0) || (Ebm.Current == 0) || (Ebm.Gamma == 0) || (Ebm.GammaEm2 == 0);
	if(ThinElBeamIsNotSetUp) return THIN_EL_BEAM_WAS_NOT_SET_UP;

	char ThickElBeamIsNotSetUp = (Ebm.Mxx < 0) || (Ebm.Mxpxp < 0) || (Ebm.Mzz < 0) || (Ebm.Mzpzp < 0);
	if(ThickElBeamIsNotSetUp) return THICK_EL_BEAM_WAS_NOT_SET_UP;

	if(!(ConstTrjDatPtr->HorFieldIsNotZero || ConstTrjDatPtr->VerFieldIsNotZero)) return MAG_FIELD_CAN_NOT_BE_ZERO;

	return 0;
}

//*************************************************************************

int srTRadIntConst::ComputeTotalStokesDistr(srTStokesStructAccessData& StokesAccessData)
{// m, eV here !!!
	int result;
	Initialize();
	if(result = CheckInputConsistency()) return result;

	SetupNativeRotation();
	SetIntegPrecLevel();
	SetupNormalizingConst();
	if(result = SetupThickBeamConsts()) return result;

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	AnalyzeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ);
	double xc = ConstTrjDatPtr->EbmDat.dxds0*(DistrInfoDat.yStart - ConstTrjDatPtr->EbmDat.s0) + ConstTrjDatPtr->EbmDat.x0;
	double zc = ConstTrjDatPtr->EbmDat.dzds0*(DistrInfoDat.yStart - ConstTrjDatPtr->EbmDat.s0) + ConstTrjDatPtr->EbmDat.z0;
	double xTol = StokesAccessData.xStep*0.01, zTol = StokesAccessData.zStep*0.01; // To steer

	//long PerE = 4;
	//long PerX = StokesAccessData.ne*PerE;
	//long PerZ = StokesAccessData.nx*PerX;
	long long PerE = 4;
	long long PerX = StokesAccessData.ne*PerE;
	long long PerZ = StokesAccessData.nx*PerX;

	EXZ.z = StokesAccessData.zStart;

	//long TotalAmOfOutPoints = DistrInfoDat.nz*DistrInfoDat.nx*DistrInfoDat.nLamb;
	long long TotalAmOfOutPoints = ((long long)DistrInfoDat.nz)*((long long)DistrInfoDat.nx)*((long long)DistrInfoDat.nLamb);
	if(FinalResAreSymOverX) TotalAmOfOutPoints >>= 1;
	if(FinalResAreSymOverZ) TotalAmOfOutPoints >>= 1;

	//long PointCount = 0;
	long long PointCount = 0;
	double UpdateTimeInt_s = 0.5;
	srTCompProgressIndicator CompProgressInd(TotalAmOfOutPoints, UpdateTimeInt_s);

	for(int iz=0; iz<DistrInfoDat.nz; iz++)
	{
		if(FinalResAreSymOverZ) { if((EXZ.z - zc) > zTol) break;}

		//long izPerZ = iz*PerZ;
		long long izPerZ = iz*PerZ;
		EXZ.x = StokesAccessData.xStart;
		for(int ix=0; ix<DistrInfoDat.nx; ix++)
		{
			if(FinalResAreSymOverX) { if((EXZ.x - xc) > xTol) break;}

			PobsLocG.x = EXZ.x; PobsLocG.y = 0.; PobsLocG.z = EXZ.z;
			PobsLocG = TrLab2Loc.TrPoint(PobsLocG);
			
			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;
			EXZ.e = StokesAccessData.eStart;
			for(int ie=0; ie<DistrInfoDat.nLamb; ie++)
			{
				float* pStokes = StokesAccessData.pBaseSto + (izPerZ + ixPerX +ie*PerE);
				if(result = ComputeStokesAtPoint(pStokes)) return result;

				if(result = srYield.Check()) return result;
				if(result = CompProgressInd.UpdateIndicator(PointCount++)) return result;
				
				EXZ.e += StokesAccessData.eStep;
			}
			EXZ.x += StokesAccessData.xStep;
		}
		EXZ.z += StokesAccessData.zStep;
	}

	if(FinalResAreSymOverZ || FinalResAreSymOverX) 
		FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, StokesAccessData);

	double ElBeamMomFact = 1.;
	if(result = TreatFiniteElecBeamEmittance(StokesAccessData, ElBeamMomFact)) return result;

	return 0;
}

//*************************************************************************

void srTRadIntConst::SetupNativeRotation()
{
	TVector3d uLab(ConstTrjDatPtr->MagConst.Bx, 0., ConstTrjDatPtr->MagConst.Bz);
	double Bcon = sqrt(uLab.x*uLab.x + uLab.z*uLab.z);
	uLab = (1./Bcon)*uLab;
	TVector3d uLoc(0.,0.,1.), Zero(0.,0.,0.);
	BconG = Bcon;

	TVector3d uLab_mi_uLoc = uLab - uLoc;
	double NormEst = uLab_mi_uLoc.Abs();
	const double RelTolNotRotate = 1.E-05; // To steer
	if(NormEst < RelTolNotRotate)
	{
		//TrLab2Loc.SetupUnit();
		TrLab2Loc.SetupIdent();
	}
	else
	{
		TVector3d AxVect = uLab^uLoc;
		double uLab_uLoc = uLab*uLoc;
		double Angle = acos(uLab_uLoc);
		
		TrLab2Loc.SetupRotation(Zero, AxVect, Angle);
		TVector3d TestVect = TrLab2Loc.TrBiPoint(uLab);
		double TestScal = TestVect*uLab;
		if(TestScal < 0.)
		{
			TrLab2Loc.SetupRotation(Zero, AxVect, -Angle);
		}
	}

	TVector3d InitCoord(ConstTrjDatPtr->EbmDat.x0, 0., ConstTrjDatPtr->EbmDat.z0);
	InitCoord = TrLab2Loc.TrPoint(InitCoord);
	
	TVector3d InitAng(ConstTrjDatPtr->EbmDat.dxds0, 0., ConstTrjDatPtr->EbmDat.dzds0);
	InitAng = TrLab2Loc.TrPoint(InitAng);
	
	LocInitCoordAng.x0 = InitCoord.x; LocInitCoordAng.z0 = InitCoord.z;
	LocInitCoordAng.dxds0 = InitAng.x; LocInitCoordAng.dzds0 = InitAng.z;
}

//*************************************************************************

void srTRadIntConst::AnalyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;
	
	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;
	ConstTrjDatPtr->MagConst.AnalyzeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
	if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;

	char ObsIsSymOverX = 0, ObsIsSymOverZ = 0;
	if(FieldIsSymOverX && (DistrInfoDat.nx > 1))
	{
		double xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
		double xTol = xStep*0.01; // To steer
		double xc = ConstTrjDatPtr->EbmDat.dxds0*(DistrInfoDat.yStart - ConstTrjDatPtr->EbmDat.s0) + ConstTrjDatPtr->EbmDat.x0;
		ObsIsSymOverX = (::fabs(0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd) - xc) < xTol);
	}
	if(FieldIsSymOverZ && (DistrInfoDat.nz > 1))
	{
		double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
		double zTol = zStep*0.01; // To steer
		double zc = ConstTrjDatPtr->EbmDat.dzds0*(DistrInfoDat.yStart - ConstTrjDatPtr->EbmDat.s0) + ConstTrjDatPtr->EbmDat.z0;
		ObsIsSymOverZ = (::fabs(0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd) - zc) < zTol);
	}
	FinalResAreSymOverX = (FieldIsSymOverX && ObsIsSymOverX);
	FinalResAreSymOverZ = (FieldIsSymOverZ && ObsIsSymOverZ);
}

//*************************************************************************

void srTRadIntConst::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData& StokesAccessData)
{
	//long PerX = StokesAccessData.ne << 2;
	//long PerZ = PerX*DistrInfoDat.nx;
	long long PerX = StokesAccessData.ne << 2;
	long long PerZ = PerX*DistrInfoDat.nx;

	char SymWithRespectToXax, SymWithRespectToZax;
	
	int HalfNz = StokesAccessData.nz >> 1, Nz_mi_1 = StokesAccessData.nz - 1;
	//int izStart = ((HalfNz << 1) == StokesAccessData.nz)? HalfNz : (HalfNz + 1);
	int HalfNx = StokesAccessData.nx >> 1, Nx_mi_1 = StokesAccessData.nx - 1;
	//int ixStart = ((HalfNx << 1) == StokesAccessData.nx)? HalfNx : (HalfNx + 1);
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
					float* pOrigData = StokesAccessData.pBaseSto + (izPerZ + ix*PerX);
					float* pSymData = StokesAccessData.pBaseSto + (izPerZ + (Nx_mi_1 - ix)*PerX);
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
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
				long long ixPerX = ix*PerX;
				float* pOrigData = StokesAccessData.pBaseSto + (izPerZ + ixPerX);
				float* pSymData = StokesAccessData.pBaseSto + (BufZ + ixPerX);
				CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
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
				float* pOrigData = StokesAccessData.pBaseSto + (izPerZ + ix*PerX);
				float* pSymData = StokesAccessData.pBaseSto + (izPerZ + (Nx_mi_1 - ix)*PerX);
				CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
			}
		}
	}
}

//*************************************************************************

int srTRadIntConst::ComputeStokesAtPoint(float* pStokes)
{
	int result;
	srTEFourier E_Loc;
	if(result = ComputeElFieldAtPointLocFrame(E_Loc)) return result;

	TVector3d ReE(E_Loc.EwX_Re, 0., E_Loc.EwZ_Re), ImE(E_Loc.EwX_Im, 0., E_Loc.EwZ_Im);
	ReE = TrLab2Loc.TrBiPoint_inv(ReE);
	ImE = TrLab2Loc.TrBiPoint_inv(ImE);

	srTEFourier E_Lab(ReE.x, ImE.x, ReE.z, ImE.z);

	srTStokes StokesLab;
	E2Stokes(E_Lab, StokesLab);
	*pStokes = (float)StokesLab.s0; 
	*(pStokes+1) = (float)StokesLab.s1; 
	*(pStokes+2) = (float)StokesLab.s2; 
	*(pStokes+3) = (float)StokesLab.s3;
	return 0;
}

//*************************************************************************

int srTRadIntConst::ComputeElFieldAtPointLocFrame(srTEFourier& E)
{
	const double Pi = 3.1415926535898;

	double PIdLamb_Invm, ActNormConst;
	if(DistrInfoDat.TreatLambdaAsEnergyIn_eV)
	{
		//PIdLamb_Invm = (Pi/1.239854E-06)*EXZ.e;
		PIdLamb_Invm = (Pi/1.239842E-06)*EXZ.e;
		ActNormConst = NormalizingConstG*EXZ.e*0.80654658E-03;
	}
	else
	{
		PIdLamb_Invm = Pi*1.E+09/EXZ.e;
		ActNormConst = NormalizingConstG/EXZ.e;
	}

	double ConArg = 2.*PIdLamb_Invm*RmaG/3.;
	double y1Inv = 1./(DistrInfoDat.yStart - ConstTrjDatPtr->EbmDat.s0);
	double Xob_mi_x0 = PobsLocG.x - LocInitCoordAng.x0;
	double Zob_mi_z0 = PobsLocG.z - LocInitCoordAng.z0;
	double Xang = Xob_mi_x0*y1Inv - LocInitCoordAng.dxds0;
	double Zang = Zob_mi_z0*y1Inv - LocInitCoordAng.dzds0;
	double GamEm2pZangE2 = ConstTrjDatPtr->EbmDat.GammaEm2 + Zang*Zang;

	double Phase = PIdLamb_Invm*((Xob_mi_x0*Xob_mi_x0 + Zob_mi_z0*Zob_mi_z0)*y1Inv + RmaG*Xang*(Xang*Xang/3. + GamEm2pZangE2));
	double Arg = ConArg*pow(GamEm2pZangE2, 1.5);

	//srTSpecialFunctions SpecFunc;
	double K2d3, K1d3;
	//srTMathFunctions::Kmu(0, 2./3., Arg, K2d3);
	CGenMathFunc::Kmu(0, 2./3., Arg, K2d3);
	//srTMathFunctions::Kmu(0, 1./3., Arg, K1d3);
	CGenMathFunc::Kmu(0, 1./3., Arg, K1d3);

	double ConA = -1.15470053838*RmaG*y1Inv*ActNormConst;
	double Ax = ConA*GamEm2pZangE2*K2d3;
	//double Az = ConA*Zang*sqrt(GamEm2pZangE2)*K1d3;
	double Az = -ConA*Zang*sqrt(GamEm2pZangE2)*K1d3;

	double CosPh, SinPh; CosAndSinComp.CosAndSin(Phase, CosPh, SinPh);

	E.EwX_Re = Ax*CosPh; E.EwX_Im = Ax*SinPh;
	//E.EwZ_Re = Az*CosPh; E.EwZ_Im = Az*SinPh;
	E.EwZ_Re = -Az*SinPh; E.EwZ_Im = Az*CosPh;
	return 0;
}

//*************************************************************************

int srTRadIntConst::TreatFiniteElecBeamEmittance(srTStokesStructAccessData& StokesAccessData, double ElBeamMomFact)
{
	int result;
	if((StokesAccessData.nx == 1) && (StokesAccessData.nz == 1)) return 0;

	//long LenFloatArr = (StokesAccessData.nx*StokesAccessData.nz);
	long long LenFloatArr = ((long long)StokesAccessData.nx)*((long long)StokesAccessData.nz);
	float *StokesCmpnArr = new float[LenFloatArr];
	if(StokesCmpnArr == 0) return MEMORY_ALLOCATION_FAILURE;

	for(int ie=0; ie<StokesAccessData.ne; ie++)
	{
		for(int is=0; is<4; is++)
		{
			ExtractStokesSliceConstE(StokesAccessData, ie, is, StokesCmpnArr);

			if((StokesAccessData.nx == 1) || (StokesAccessData.nz == 1))
			{
				if(result = TreatFiniteElecBeamEmittanceOneComp1D(StokesCmpnArr, ElBeamMomFact)) return result;
			}
			else
			{
				if(result = TreatFiniteElecBeamEmittanceOneComp2D(StokesCmpnArr, ElBeamMomFact)) return result;
			}
			if(is == 0) SuppressNegativeValues(StokesCmpnArr);

			UpdateStokesSliceConstE(StokesCmpnArr, ie, is, StokesAccessData);
		}
	}
	delete[] StokesCmpnArr;
	return 0;
}

//*************************************************************************

void srTRadIntConst::ExtractStokesSliceConstE(srTStokesStructAccessData& StokesAccessData, long ie, int StokesNo, float* pOutS)
{
	float *pS0 = StokesAccessData.pBaseSto + StokesNo;
	//long PerX = StokesAccessData.ne << 2;
	//long PerZ = PerX*StokesAccessData.nx;
	long long PerX = StokesAccessData.ne << 2;
	long long PerZ = PerX*StokesAccessData.nx;

	//long izPerZ = 0;
	//long iePerE = ie << 2;
	long long izPerZ = 0;
	long long iePerE = ie << 2;

	float *tOutS = pOutS;
	for(int iz=0; iz<StokesAccessData.nz; iz++)
	{
		float *pS_StartForX = pS0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<StokesAccessData.nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pS = pS_StartForX + ixPerX_p_iePerE;
			*(tOutS++) = *pS;

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
}

//*************************************************************************

void srTRadIntConst::UpdateStokesSliceConstE(float* StokesCmpnArr, long ie, int StokesNo, srTStokesStructAccessData& StokesAccessData)
{
	float *pS0 = StokesAccessData.pBaseSto + StokesNo;

	//long PerX = StokesAccessData.ne << 2;
	//long PerZ = PerX*StokesAccessData.nx;
	long long PerX = StokesAccessData.ne << 2;
	long long PerZ = PerX*StokesAccessData.nx;

	//long izPerZ = 0;
	//long iePerE = ie << 2;
	long long izPerZ = 0;
	long long iePerE = ie << 2;

	float *tCmpn = StokesCmpnArr;
	for(int iz=0; iz<StokesAccessData.nz; iz++)
	{
		float *pS_StartForX = pS0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(int ix=0; ix<StokesAccessData.nx; ix++)
		{
			//long ixPerX_p_iePerE = ixPerX + iePerE;
			long long ixPerX_p_iePerE = ixPerX + iePerE;
			float *pS = pS_StartForX + ixPerX_p_iePerE;
			*pS = *(tCmpn++);

			ixPerX += PerX;
		}
		izPerZ += PerZ;
	}
}

//*************************************************************************

int srTRadIntConst::TreatFiniteElecBeamEmittanceOneComp1D(float* CmpnArr, double ElBeamMomFact)
{
	int result;
	char VsXorZ = 0;
	if(DistrInfoDat.nx > 1) VsXorZ = 'x';
	else if(DistrInfoDat.nz > 1) VsXorZ = 'z';
	else return 0;

	double M_ElecEff = ElBeamMomFact*((VsXorZ == 'x')? 0.5/Gx : 0.5/Gz);
	double M_DistrSingleE;
	DetermineSingleElecDistrEffSizes1D(CmpnArr, VsXorZ, M_DistrSingleE);

	srTRadResize1D Resize;
	DetermineResizeBeforeConv1D(M_ElecEff, M_DistrSingleE, VsXorZ, Resize);

	long Np = (VsXorZ == 'x')? DistrInfoDat.nx : DistrInfoDat.nz;
	long NpAux = (long)(Resize.pm*Np);

	CGenMathFFT1D FFT;
	FFT.NextCorrectNumberForFFT(NpAux);
	float* AuxConvData = new float[NpAux << 1];
	if(AuxConvData == 0) return MEMORY_ALLOCATION_FAILURE;

	ConstructDataForConv1D(CmpnArr, AuxConvData, Np, NpAux);
	if(result = PerformConvolutionWithGaussian1D(AuxConvData, NpAux, M_ElecEff, VsXorZ)) return result;

	ExtractDataAfterConv1D(AuxConvData, NpAux, Np, CmpnArr);

	if(AuxConvData != 0) delete[] AuxConvData;
	return 0;
}

//*************************************************************************

void srTRadIntConst::DetermineSingleElecDistrEffSizes1D(float* CmpnArr, char VsXorZ, double& M_Cen)
{
	long Np;
	double Step, Arg;
	if(VsXorZ == 'x') 
	{
		Np = DistrInfoDat.nx;
		Step = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
		Arg = DistrInfoDat.xStart;
	}
	else
	{
		Np = DistrInfoDat.nz;
		Step = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
		Arg = DistrInfoDat.zStart;
	}
	long Np_mi_1 = Np - 1;
	float Sm0 = 0., Sm1 = 0., Sm2 = 0.;
	float *tDistr = CmpnArr;
	//float ArgE2 = (float)(Arg*Arg);
	float w = 0.5;

	for(int i=0; i<Np; i++)
	{
		if(i == Np_mi_1) w = 0.5;

		float Distr = (float)::fabs(w*(*(tDistr++)));
		Sm0 += Distr;
		Sm1 += (float)(Arg*Distr);
		Sm2 += (float)(Arg*Arg*Distr);

		w = 1.;
		Arg += Step;
	}
	float Norm = (float)(1./Sm0);
	float M1 = Sm1*Norm, M2 = Sm2*Norm;
	M_Cen = M2 - M1*M1;
}

//*************************************************************************

void srTRadIntConst::DetermineResizeBeforeConv1D(double M_ElecEff, double M_DistrSingleE, char VsXorZ, srTRadResize1D& Resize)
{
	const double DoNotResizeRatio = 3.; // To steer
	const double AmOfExtraSig = 2.5;

	if(M_DistrSingleE*DoNotResizeRatio > M_ElecEff)
	{
		double ExtraRange = 2.*AmOfExtraSig*sqrt(M_ElecEff);
		double CurrentRange = (VsXorZ == 'x')? (DistrInfoDat.xEnd - DistrInfoDat.xStart) : (DistrInfoDat.zEnd - DistrInfoDat.zStart);
		Resize.pm = (CurrentRange + ExtraRange)/CurrentRange;
	}
}

//*************************************************************************

void srTRadIntConst::ConstructDataForConv1D(float* CmpnArr, float* AuxConvData, long NpOld, long NpNew)
{
	long iDat = (NpNew - NpOld) >> 1;

	long i;
	float V = *CmpnArr;
	float *tNew = AuxConvData;
	for(i=0; i<iDat; i++)
	{
		*(tNew++) = V; *(tNew++) = 0;
	}
	float *tOld = CmpnArr;
	for(i=iDat; i<(iDat + NpOld); i++)
	{
		*(tNew++) = *(tOld++); *(tNew++) = 0;
	}
	V = *(CmpnArr + NpOld - 1);
	for(i=(iDat + NpOld); i<NpNew; i++)
	{
		*(tNew++) = V; *(tNew++) = 0;
	}
}

//*************************************************************************

int srTRadIntConst::PerformConvolutionWithGaussian1D(float* ConvData, long NewNp, double M_ElecEff, char VsXorZ)
{
	int result;
	const double Pi = 3.14159265358979;
	const double TwoPi = Pi*2.;
	const double TwoPiE2 = TwoPi*Pi;

	double Step = (VsXorZ == 'x')? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
	double StartFict = -Step*(NewNp >> 1);

	long TwoNewNp = NewNp << 1;
	float* AuxData = new float[TwoNewNp];
	if(AuxData == 0) return MEMORY_ALLOCATION_FAILURE;

	CGenMathFFT1DInfo FFT1DInfo;
	FFT1DInfo.pInData = ConvData;
	FFT1DInfo.pOutData = AuxData;
	FFT1DInfo.Dir = 1;
	FFT1DInfo.xStep = Step;
	FFT1DInfo.xStart = StartFict;
	FFT1DInfo.Nx = NewNp;
	FFT1DInfo.HowMany = 1;
	FFT1DInfo.UseGivenStartTrValue = 0;
	FFT1DInfo.TreatSharpEdges = 0;
	CGenMathFFT1D FFT1D;
	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;

	double C2 = TwoPiE2*M_ElecEff;
	float* tData = AuxData;
	double q = FFT1DInfo.xStartTr;
	for(long j=0; j<NewNp; j++)
	{
		//double C2qE2 = C2*q*q;
		double Magn = exp(-C2*q*q);
		*(tData++) *= (float)Magn; // Re
		*(tData++) *= (float)Magn; // Im

		q += FFT1DInfo.xStepTr;
	}

	FFT1DInfo.pInData = AuxData;
	FFT1DInfo.pOutData = ConvData;
	FFT1DInfo.Dir = -1;
	FFT1DInfo.xStep = FFT1DInfo.xStepTr; FFT1DInfo.xStepTr = Step;
	FFT1DInfo.xStart = FFT1DInfo.xStartTr; FFT1DInfo.xStartTr = StartFict;
	FFT1DInfo.UseGivenStartTrValue = 1;

	if(result = FFT1D.Make1DFFT(FFT1DInfo)) return result;
	delete[] AuxData;
	return 0;
}

//*************************************************************************

void srTRadIntConst::ExtractDataAfterConv1D(float* AuxConvData, long NpAux, long Np, float* CmpnArr)
{
	long iDat = (NpAux - Np) >> 1;

	float *tCmpn = CmpnArr;
	float *tAux = AuxConvData + (iDat << 1);
	for(long i=0; i<Np; i++)
	{
		*(tCmpn++) = *tAux; tAux += 2;
	}
}

//*************************************************************************

int srTRadIntConst::TreatFiniteElecBeamEmittanceOneComp2D(float* CmpnArr, double ElBeamMomFact)
{
	int result;
	double MxxElecEff = ElBeamMomFact*0.5/Gx, MzzElecEff = ElBeamMomFact*0.5/Gz;

	double MxxPowSingleE, MzzPowSingleE;
	DetermineSingleElecDistrEffSizes2D(CmpnArr, MxxPowSingleE, MzzPowSingleE);

	srTRadResize Resize;
	DetermineResizeBeforeConv2D(MxxElecEff, MzzElecEff, MxxPowSingleE, MzzPowSingleE, Resize);

	long NxAux = (long)(Resize.pxm*DistrInfoDat.nx);
	long NzAux = (long)(Resize.pzm*DistrInfoDat.nz);
	CGenMathFFT2D FFT;
	FFT.NextCorrectNumberForFFT(NxAux);
	FFT.NextCorrectNumberForFFT(NzAux);
	float* AuxConvData = new float[NxAux*NzAux << 1];
	if(AuxConvData == 0) return MEMORY_ALLOCATION_FAILURE;

	ConstructDataForConv2D(CmpnArr, AuxConvData, NxAux, NzAux);
	if(result = PerformConvolutionWithGaussian2D(AuxConvData, NxAux, NzAux, MxxElecEff, MzzElecEff)) return result;

	ExtractDataAfterConv2D(AuxConvData, NxAux, NzAux, CmpnArr);

	if(AuxConvData != 0) delete[] AuxConvData;
	return 0;
}

//*************************************************************************

void srTRadIntConst::DetermineSingleElecDistrEffSizes2D(float* CmpnArr, double& MxxCen, double& MzzCen)
{// Compute central moments
	float Sm0 = 0., SmX = 0., SmZ = 0., SmXX = 0., SmZZ = 0.;
	float xStep = (float)((DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.);
	float zStep = (float)((DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.);
	long nz_mi_1 = DistrInfoDat.nz - 1, nx_mi_1 = DistrInfoDat.nx - 1;

	float *tDistr = CmpnArr;
	float z = (float)DistrInfoDat.zStart;
	float ze2 = z*z;
	float wz = 0.5;

	for(int iz=0; iz<DistrInfoDat.nz; iz++)
	{
		if(iz == nz_mi_1) wz = 0.5;
		
		float x = (float)DistrInfoDat.xStart;
		float xe2 = x*x;
		for(int ix=0; ix<DistrInfoDat.nx; ix++)
		{
			float Distr = (float)::fabs(wz*(*(tDistr++)));
			if((ix == nx_mi_1) || (iz == nz_mi_1)) Distr *= 0.5;

			Sm0 += Distr;
			SmX += x*Distr;
			SmXX += xe2*Distr;
			SmZ += z*Distr;
			SmZZ += ze2*Distr;

			x += xStep;
			xe2 = x*x;
		}
		wz = 1.;
		z += zStep;
		ze2 = z*z;
	}
	float Norm = (float)(1./Sm0);
	float Mx = SmX*Norm, Mxx = SmXX*Norm;
	float Mz = SmZ*Norm, Mzz = SmZZ*Norm;

	MxxCen = Mxx - Mx*Mx;
	MzzCen = Mzz - Mz*Mz;
}

//*************************************************************************

void srTRadIntConst::DetermineResizeBeforeConv2D(double MxxElecEff, double MzzElecEff, double MxxPowSingleE, double MzzPowSingleE, srTRadResize& Resize)
{
	const double DoNotResizeRatio = 3.; // To steer
	const double AmOfExtraSig = 2.5;

	if(MxxPowSingleE*DoNotResizeRatio > MxxElecEff)
	{
		double ExtraRange = 2.*AmOfExtraSig*sqrt(MxxElecEff);
		double CurrentRange = DistrInfoDat.xEnd - DistrInfoDat.xStart;
		Resize.pxm = (CurrentRange + ExtraRange)/CurrentRange;
	}
	if(MzzPowSingleE*DoNotResizeRatio > MzzElecEff)
	{
		double ExtraRange = 2.*AmOfExtraSig*sqrt(MzzElecEff);
		double CurrentRange = DistrInfoDat.zEnd - DistrInfoDat.zStart;
		Resize.pzm = (CurrentRange + ExtraRange)/CurrentRange;
	}
}

//*************************************************************************

void srTRadIntConst::ConstructDataForConv2D(float* CmpnArr, float* NewData, long NewNx, long NewNz)
{
	long ixDat = (NewNx - DistrInfoDat.nx) >> 1;
	long izDat = (NewNz - DistrInfoDat.nz) >> 1;
	long OffsetXp = ixDat + DistrInfoDat.nx;
	long OffsetZp = izDat + DistrInfoDat.nz;

	//long NewPerZ = NewNx << 1;
	long long NewPerZ = NewNx << 1;
	long ix, iz;
	float V = *CmpnArr;
	for(iz=0; iz<izDat; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ);
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(CmpnArr + DistrInfoDat.nx - 1);
	for(iz=0; iz<izDat; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(CmpnArr + (DistrInfoDat.nz - 1)*DistrInfoDat.nx);
	for(iz=OffsetZp; iz<NewNz; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ);
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(CmpnArr + DistrInfoDat.nz*DistrInfoDat.nx - 1);
	for(iz=OffsetZp; iz<NewNz; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}

	float *tOldL = CmpnArr;
	float *tOldR = CmpnArr + DistrInfoDat.nx - 1;
	for(iz=izDat; iz<(izDat + DistrInfoDat.nz); iz++)
	{
		V = *tOldL;
		//long izNewPerZ = iz*NewPerZ;
		long long izNewPerZ = iz*NewPerZ;
		float *tNew = NewData + izNewPerZ;
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
		tOldL += DistrInfoDat.nx;

		V = *tOldR;
		tNew = NewData + (izNewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
		tOldR += DistrInfoDat.nx;
	}

	float *tOldD = CmpnArr;
	//float *tOldU = CmpnArr + (DistrInfoDat.nz - 1)*DistrInfoDat.nx;
	float *tOldU = CmpnArr + (((long long)(DistrInfoDat.nz - 1))*((long long)DistrInfoDat.nx));
	for(ix=ixDat; ix<(ixDat + DistrInfoDat.nx); ix++)
	{
		V = *(tOldD++);
		float *tNew = NewData + (ix << 1);
		for(iz=0; iz<izDat; iz++) { *tNew = V; *(tNew+1) = 0; tNew += NewPerZ;}

		V = *(tOldU++);
		tNew = NewData + (OffsetZp*NewPerZ + (ix << 1));
		for(iz=OffsetZp; iz<NewNz; iz++) { *tNew = V; *(tNew+1) = 0; tNew += NewPerZ;}
	}

	float *tOld = CmpnArr;
	for(iz=izDat; iz<(izDat + DistrInfoDat.nz); iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (ixDat << 1));
		for(ix=ixDat; ix<(ixDat + DistrInfoDat.nx); ix++)
		{
			*(tNew++) = *(tOld++); *(tNew++) = 0.;
		}
	}
}

//*************************************************************************

int srTRadIntConst::PerformConvolutionWithGaussian2D(float* ConvData, long NewNx, long NewNz, double MxxElecEff, double MzzElecEff)
{
	int result;
	double xStep = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
	double zStep = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;

	double xStartFict = -xStep*(NewNx >> 1);
	double zStartFict = -zStep*(NewNz >> 1);

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.pData = ConvData;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.xStep = xStep;
	FFT2DInfo.yStep = zStep;
	FFT2DInfo.xStart = xStartFict;
	FFT2DInfo.yStart = zStartFict;
	FFT2DInfo.Nx = NewNx;
	FFT2DInfo.Ny = NewNz;
	FFT2DInfo.UseGivenStartTrValues = 0;

	CGenMathFFT2D FFT2D;
	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;

	const double Pi = 3.14159265358979;
	const double TwoPi = Pi*2.;
	const double TwoPiE2 = TwoPi*Pi;

	double C2x = TwoPiE2*MxxElecEff;
	double C2z = TwoPiE2*MzzElecEff;

	float* tData = ConvData;
	double qz = FFT2DInfo.yStartTr;
	for(long iz=0; iz<NewNz; iz++)
	{
		double C2zqzE2 = C2z*qz*qz;
		double qx = FFT2DInfo.xStartTr;

		for(long ix=0; ix<NewNx; ix++)
		{
			double Magn = exp(-C2x*qx*qx - C2zqzE2);
			*(tData++) *= (float)Magn; // Re
			*(tData++) *= (float)Magn; // Im

			qx += FFT2DInfo.xStepTr;
		}
		qz += FFT2DInfo.yStepTr;
	}

	FFT2DInfo.pData = ConvData;
	FFT2DInfo.Dir = -1;
	FFT2DInfo.xStep = FFT2DInfo.xStepTr; FFT2DInfo.xStepTr = xStep;
	FFT2DInfo.yStep = FFT2DInfo.yStepTr; FFT2DInfo.yStepTr = zStep;
	FFT2DInfo.xStart = FFT2DInfo.xStartTr; FFT2DInfo.xStartTr = xStartFict;
	FFT2DInfo.yStart = FFT2DInfo.yStartTr; FFT2DInfo.yStartTr = zStartFict;
	FFT2DInfo.UseGivenStartTrValues = 1;

	return FFT2D.Make2DFFT(FFT2DInfo);
	return 0;
}

//*************************************************************************

void srTRadIntConst::ExtractDataAfterConv2D(float* AuxConvData, long NxAux, long NzAux, float* CmpnArr)
{
	long ixDat = (NxAux - DistrInfoDat.nx) >> 1;
	long izDat = (NzAux - DistrInfoDat.nz) >> 1;
	long Two_ixDat = ixDat << 1;
	//long AuxPerZ = NxAux << 1;
	long long AuxPerZ = NxAux << 1;

	float *tCmpn = CmpnArr;
	for(long iz=0; iz<DistrInfoDat.nz; iz++)
	{
		//float *tAux = AuxConvData + (izDat + iz)*AuxPerZ + Two_ixDat;
		float *tAux = AuxConvData + ((izDat + iz)*AuxPerZ + Two_ixDat);
		for(long ix=0; ix<DistrInfoDat.nx; ix++)
		{
			*(tCmpn++) = *tAux; tAux += 2;
		}
	}
}

//*************************************************************************

void srTRadIntConst::SuppressNegativeValues(float* StokesCmpnArr)
{
	float *t = StokesCmpnArr;
	for(int iz=0; iz<DistrInfoDat.nz; iz++)
		for(int ix=0; ix<DistrInfoDat.nx; ix++)
		{
			if(*t < 0.) *t = 0.;
			t++;
		}
}

//*************************************************************************
