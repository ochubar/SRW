/************************************************************************//**
 * File: srstowig.cpp
 * Description: Calculation of Stokes parameters of Wiggler Radiation
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srstowig.h"
#include "srprgind.h"
#include "gmfft.h"
#include "srerror.h"

//*************************************************************************

int srTRadIntWiggler::CheckInputConsistency()
{
// Checking El. Beam
	srTEbmDat& Ebm = TrjDatPtr->EbmDat;
	char ThinElBeamIsNotSetUp = (Ebm.Energy <= 0) || (Ebm.Current == 0) || (Ebm.Gamma == 0) || (Ebm.GammaEm2 == 0);
	if(ThinElBeamIsNotSetUp) return THIN_EL_BEAM_WAS_NOT_SET_UP;

	char ThickElBeamIsNotSetUp = (Ebm.Mxx <= 0) || (Ebm.Mxpxp == 0) || (Ebm.Mzz <= 0) || (Ebm.Mzpzp == 0);
	if(ThickElBeamIsNotSetUp) return THICK_EL_BEAM_WAS_NOT_SET_UP;

	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeTotalStokesDistr(srTStokesStructAccessData& StokesAccessData)
{// m, eV here !!!
	int result;
	Initialize();

	if(result = CheckInputConsistency()) return result;
	if(result = TrjDatPtr->ComputeInterpolatingStructure()) return result;

	if(result = SetUpFieldBasedArrays()) return result;
	if(result = AllocateIntervalsArray(FieldBasedArrays.Ns)) return result;

	SetIntegPrecLevel();
	SetupNormalizingConst();

	if(IncludeCrossTermsG)
	{
		if(result = SetUpCrossTermsContribArray()) return result;
	}
	if(result = SetupThickBeamConsts()) return result;

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	AnalizeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ);
	double xc = TrjDatPtr->EbmDat.dxds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0) + TrjDatPtr->EbmDat.x0;
	double zc = TrjDatPtr->EbmDat.dzds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0) + TrjDatPtr->EbmDat.z0;
	double xTol = StokesAccessData.xStep*0.01, zTol = StokesAccessData.zStep*0.01; // To steer
	
	MaxIntValWithinPointG = MinEstIntValWithinPointG;

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
			
			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;
			EXZ.e = StokesAccessData.eStart;
			for(int ie=0; ie<DistrInfoDat.nLamb; ie++)
			{
				float* pStokes = StokesAccessData.pBaseSto + izPerZ + ixPerX +ie*PerE;
				if(result = ComputeStokesAtPoint(pStokes)) return result;

				if(result = srYield.Check()) return result;
				if(result = CompProgressInd.UpdateIndicator(PointCount++)) return result;
				
				EXZ.e += StokesAccessData.eStep;
			}
			EXZ.x += StokesAccessData.xStep;
		}
		EXZ.z += StokesAccessData.zStep;
	}
	DeallocateIntervalsArray();
	FieldBasedArrays.DisposeArrays();

	if(FinalResAreSymOverZ || FinalResAreSymOverX) 
		FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, StokesAccessData);

	char MainOrCrossTerms = 'm';
	double ElBeamMomFact = 1.;
	if(result = TreatFiniteElecBeamEmittance(StokesAccessData, MainOrCrossTerms, ElBeamMomFact)) return result;

	if(IncludeCrossTermsG)
	{
		// Smoothing Cross Terms data
		ElBeamMomFact = 0.4; // To steer
		MainOrCrossTerms = 'c';
		if(result = TreatFiniteElecBeamEmittance(StokesAccessData, MainOrCrossTerms, ElBeamMomFact)) return result;

		AddCrossTermsContrib(StokesAccessData);
		DeallocateCrossTermsContribArray();
	}

	return 0;
}

//*************************************************************************

int srTRadIntWiggler::SetUpFieldBasedArrays()
{
	int result;
	double sStepLoc;
	//long NsLoc;
	long long NsLoc;
	if(result = TrjDatPtr->SetupLimitsByAnalizingField(LongIntTypeG, FieldBasedArrays.sStart, sStepLoc, NsLoc, FieldBasedArrays.Nper, FieldBasedArrays.NperLeft)) return result;
	double sEndLoc = FieldBasedArrays.sStart + sStepLoc*(NsLoc - 1);

	FieldBasedArrays.Ns = NsLoc;
	FieldBasedArrays.sStep = (sEndLoc - FieldBasedArrays.sStart)/double(FieldBasedArrays.Ns - 1);
	
	srTFieldBasedArrayKeys Keys;
	Keys.Bx_ = Keys.Bz_ = Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = Keys.IntBtxE2_ = Keys.IntBtzE2_ = Keys.dBxds_ = Keys.dBzds_ = 1;
	if(result = FieldBasedArrays.AllocateArrays(FieldBasedArrays.Ns, Keys)) return result;
	TrjDatPtr->CompTotalTrjData(Keys, FieldBasedArrays);
	return 0;
}

//*************************************************************************

void srTRadIntWiggler::CheckPossibilityOfFarFieldPeriodicComp(srTGenTrjHndl& TrjHndl, char& LongIntType)
{
	const double CritMagFieldLenToObsDistRatio = 0.05; // To steer
	
	double sIntegStart, sIntegFin;
	TrjHndl.rep->ShowFullLimits(sIntegStart, sIntegFin);
	
	if(DistrInfoDat.yStart < CritMagFieldLenToObsDistRatio*(sIntegFin - sIntegStart)) LongIntType = 1; // Int over total range (near field)
	else LongIntType = 2; // Int over one period (far field)
}

//*************************************************************************

void srTRadIntWiggler::AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;
	
	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;
	TrjDatPtr->AnalizeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
	if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;
	
	char ObsIsSymOverX = 0, ObsIsSymOverZ = 0;
	if(FieldIsSymOverX && (DistrInfoDat.nx > 1))
	{
		double xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
		double xTol = xStep*0.01; // To steer
		double xc = TrjDatPtr->EbmDat.dxds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0) + TrjDatPtr->EbmDat.x0;
		ObsIsSymOverX = (::fabs(0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd) - xc) < xTol);
	}
	if(FieldIsSymOverZ && (DistrInfoDat.nz > 1))
	{
		double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
		double zTol = zStep*0.01; // To steer
		double zc = TrjDatPtr->EbmDat.dzds0*(DistrInfoDat.yStart - TrjDatPtr->EbmDat.s0) + TrjDatPtr->EbmDat.z0;
		ObsIsSymOverZ = (::fabs(0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd) - zc) < zTol);
	}
	
	FinalResAreSymOverX = (FieldIsSymOverX && ObsIsSymOverX);
	FinalResAreSymOverZ = (FieldIsSymOverZ && ObsIsSymOverZ);
}

//*************************************************************************

void srTRadIntWiggler::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData& StokesAccessData)
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
					//long OffsetOrig = izPerZ + ix*PerX;
					//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					long long OffsetOrig = izPerZ + ix*PerX;
					long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
					float* pOrigData = StokesAccessData.pBaseSto + OffsetOrig;
					float* pSymData = StokesAccessData.pBaseSto + OffsetSym;
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);

					if(IncludeCrossTermsG)
					{
						pOrigData = CrossTermsContribArray + OffsetOrig;
						pSymData = CrossTermsContribArray + OffsetSym;
						CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
					}
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
				float* pOrigData = StokesAccessData.pBaseSto + izPerZ + ixPerX;
				float* pSymData = StokesAccessData.pBaseSto + BufZ + ixPerX;
				CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				
				if(IncludeCrossTermsG)
				{
					pOrigData = CrossTermsContribArray + izPerZ + ixPerX;
					pSymData = CrossTermsContribArray + BufZ + ixPerX;
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
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
				//long OffsetOrig = izPerZ + ix*PerX;
				//long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				long long OffsetOrig = izPerZ + ix*PerX;
				long long OffsetSym = izPerZ + (Nx_mi_1 - ix)*PerX;
				float* pOrigData = StokesAccessData.pBaseSto + OffsetOrig;
				float* pSymData = StokesAccessData.pBaseSto + OffsetSym;
				CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);

				if(IncludeCrossTermsG)
				{
					pOrigData = CrossTermsContribArray + OffsetOrig;
					pSymData = CrossTermsContribArray + OffsetSym;
					CopySymEnergySlice(pOrigData, pSymData, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
	}
}

//*************************************************************************

void srTRadIntWiggler::ChooseNumIntervals(int AmOfInterv, int& StartIntervNo, int& FinIntervNo)
{
	if(AmOfInterv <= 1) { StartIntervNo = FinIntervNo = 0; return;}

	StartIntervNo = 0; FinIntervNo = AmOfInterv - 1;
	if(LongIntTypeG == 1) return; // Arbitrary field

	double sStartTrue = FieldBasedArrays.sStart;
	double sStepTrue = FieldBasedArrays.sStep;
	double sPeriod = 0.5*sStepTrue*(FieldBasedArrays.Ns - 1);

	double sIntervStart = sStartTrue;
	char ZeroIntervIsBad = ((RadNumIntervals->IndSt == 0) && (AmOfInterv > 1));

	if(ZeroIntervIsBad)
	{
		StartIntervNo = 1; 
		sIntervStart = (RadNumIntervals[1].IndSt)*sStepTrue + sStartTrue;
	}
	double sEndCrit = sIntervStart + sPeriod*0.999999; // To steer

	for(int k=1; k<AmOfInterv; k++)
	{
		double sStartCurInt = (RadNumIntervals[k].IndSt)*sStepTrue + sStartTrue;
		if(sStartCurInt >= sEndCrit)
		{
			FinIntervNo = k - 1; break; // To check
		}
	}
}

//*************************************************************************

int srTRadIntWiggler::ComputeStokesAtPoint(float* pStokes)
{
	int result;
	*pStokes = 0.; *(pStokes+1) = 0.; *(pStokes+2) = 0.; *(pStokes+3) = 0.;

	int AmOfInterv;
	DetermineRadIntervals(AmOfInterv);

	char AllAsOneNumInterv = ((AmOfInterv == 1) && (RadNumIntervals->IndSt == 0) && (RadNumIntervals->IndFi == (FieldBasedArrays.Ns - 1)));
	if(AmOfInterv == 0) 
	{
		if(IncludeCrossTermsG)
		{
			*(tCrossTermsContribG++) = 0.; 
			*(tCrossTermsContribG++) = 0.; 
			*(tCrossTermsContribG++) = 0.; 
			*(tCrossTermsContribG++) = 0.;
			// Nowhere else tCrossTermsContribG should be changed !
		}
		return 0;
	}

	char Periodic = (LongIntTypeG == 2);
	if(Periodic && AllAsOneNumInterv)
	{// Since two periods are analysed, take only the first
		RadNumIntervals->IndFi = ((FieldBasedArrays.Ns - 1) >> 1);

		RadNumIntervals->SumReEx *= 0.5;
		RadNumIntervals->SumImEx *= 0.5;
		RadNumIntervals->SumReEz *= 0.5;
		RadNumIntervals->SumImEz *= 0.5;

		RadNumIntervals->Ph_Fi = RadNumIntervals->Ph_St + 0.5*(RadNumIntervals->Ph_Fi - RadNumIntervals->Ph_St);
		CosAndSinComp.CosAndSin(RadNumIntervals->Ph_Fi, RadNumIntervals->CosPh_Fi, RadNumIntervals->SinPh_Fi);
	}

	int StartIntervNo, FinIntervNo;
	ChooseNumIntervals(AmOfInterv, StartIntervNo, FinIntervNo);

	//char MergeEdges = 0;
	srTEFourier AuxE;
	srTStokes LocStokesCen;
	for(int i=StartIntervNo; i<=FinIntervNo; i++)
	{
		srTStokes LocStokes;
		RadNumIntervals[i].E.ZeroAll();
		if(NumIntegG)
		{
			if(result = ComputeAndAddTerminations(i)) return result;
			if(result = ComputeAndAddOneMainTermNum(i)) return result;
		}
		//else ...

		E2Stokes(RadNumIntervals[i].E, LocStokes);
		LocStokesCen += LocStokes;
	}
	if(LongIntTypeG == 2) // Periodic
	{
		double NperMult = double(FieldBasedArrays.Nper);
		LocStokesCen *= NperMult;
	}
	
	*pStokes = (float)(LocStokesCen.s0); 
	*(pStokes+1) = (float)(LocStokesCen.s1); 
	*(pStokes+2) = (float)(LocStokesCen.s2); 
	*(pStokes+3) = (float)(LocStokesCen.s3);

	if(IncludeCrossTermsG)
	{
		srTStokes LocStokes;
		if(LongIntTypeG == 1) // Arbitrary field
		{
			if(result = ComputeCrossTermsAtPointArb(AmOfInterv, LocStokes)) return result;
		}
		else if(LongIntTypeG == 2) // Periodic
		{
			if(result = ComputeCrossTermsAtPointPer(StartIntervNo, FinIntervNo, LocStokes)) return result;
		}

		*(tCrossTermsContribG++) = (float)(LocStokes.s0); 
		*(tCrossTermsContribG++) = (float)(LocStokes.s1); 
		*(tCrossTermsContribG++) = (float)(LocStokes.s2); 
		*(tCrossTermsContribG++) = (float)(LocStokes.s3);
		// Nowhere else tCrossTermsContribG should be changed !
	}
	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeCrossTermsAtPointArb(int AmOfInterv, srTStokes& Stokes)
{
	//const double RelTolIgnoreComp = 1.E-07; // To steer
	const double MinArgIgnoreComp = -15.; // To steer

	double PIdLamb_Inv_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 2.533840802E+06*EXZ.e : 3.14159265359E+09/EXZ.e;
	double TwoPIdLamb_Inv_m = 2.*PIdLamb_Inv_m;

	srTEbmDat& Ebm = TrjDatPtr->EbmDat;
	double EnSprMult = 2.*PIdLamb_Inv_m*PIdLamb_Inv_m*(Ebm.GammaEm2)*(Ebm.GammaEm2)*(Ebm.SigmaRelE)*(Ebm.SigmaRelE);

	double yObs = DistrInfoDat.yStart;

	srTRadContribInterval *ti = RadNumIntervals;
	for(int i=0; i<AmOfInterv; i++)
	{
		double si = ti->sc;
		double yObs_mi_sI_Inv = 1./(yObs - si);
		double AngXI = (EXZ.x - ti->Xc)*yObs_mi_sI_Inv;
		double AngZI = (EXZ.z - ti->Zc)*yObs_mi_sI_Inv;

		srTEFourier& EI = ti->E;
		double EIxRe = EI.EwX_Re, EIxIm = EI.EwX_Im;
		double EIzRe = EI.EwZ_Re, EIzIm = EI.EwZ_Im;

		srTRadContribInterval *tj = RadNumIntervals + i + 1;
		for(int j=(i + 1); j<AmOfInterv; j++)
		{
			double sj = tj->sc;
			double si_mi_sj = si - sj;

			double yObs_mi_sJ_Inv = 1./(yObs - sj);
			double AngXJ = (EXZ.x - tj->Xc)*yObs_mi_sJ_Inv;
			double AngZJ = (EXZ.z - tj->Zc)*yObs_mi_sJ_Inv;

			srTEFourier& EJ = tj->E;
			double EJxRe = EJ.EwX_Re, EJxIm = EJ.EwX_Im;
			double EJzRe = EJ.EwZ_Re, EJzIm = EJ.EwZ_Im;

			double T2 = PIdLamb_Inv_m*si_mi_sj*yObs_mi_sI_Inv*yObs_mi_sJ_Inv;
			double T2e2 = T2*T2, T2d4 = 0.25*T2;
			double T2e2pGxE2_Inv = 1./(T2e2 + Gx*Gx);
			double T2e2pGzE2_Inv = 1./(T2e2 + Gz*Gz);

			double T1x = TwoPIdLamb_Inv_m*(-AngXI + AngXJ);
			double T1xe2 = T1x*T1x;
			double T1xe2T2e2pGxE2_Inv = T1xe2*T2e2pGxE2_Inv;
			double T1z = TwoPIdLamb_Inv_m*(-AngZI + AngZJ);
			double T1ze2 = T1z*T1z;
			double T1ze2T2e2pGzE2_Inv = T1ze2*T2e2pGzE2_Inv;

			double ArgTot = -0.25*(Gx*T1xe2T2e2pGxE2_Inv + Gz*T1ze2T2e2pGzE2_Inv) - EnSprMult*si_mi_sj*si_mi_sj;
			if(ArgTot > MinArgIgnoreComp)
			{
				double PhTot = 0.5*(atan(T2/Gx) + atan(T2/Gz)) - T2d4*(T1xe2T2e2pGxE2_Inv + T1ze2T2e2pGzE2_Inv);
				double PhTotRe, PhTotIm; CosAndSinComp.CosAndSin(PhTot, PhTotRe, PhTotIm);
				
				double nCon = 2.*sqrt(Gx*Gz*sqrt(T2e2pGxE2_Inv*T2e2pGzE2_Inv));
				double ReMult = nCon*exp(ArgTot);
				
				double EIxReEJxIm = EIxRe*EJxIm, EIxImEJxRe = EIxIm*EJxRe, EIzReEJzIm = EIzRe*EJzIm, EIzImEJzRe = EIzIm*EJzRe;
				double EIxImEJxIm = EIxIm*EJxIm, EIxReEJxRe = EIxRe*EJxRe, EIzImEJzIm = EIzIm*EJzIm, EIzReEJzRe = EIzRe*EJzRe;
				double s0 = ReMult*((EIxReEJxIm - EIxImEJxRe + EIzReEJzIm - EIzImEJzRe)*PhTotIm 
							+(EIxImEJxIm + EIxReEJxRe + EIzImEJzIm + EIzReEJzRe)*PhTotRe);
				double s1 = ReMult*((EIxReEJxIm - EIxImEJxRe - EIzReEJzIm + EIzImEJzRe)*PhTotIm 
							+(EIxImEJxIm + EIxReEJxRe - EIzImEJzIm - EIzReEJzRe)*PhTotRe);
				double EIzReEJxIm = EIzRe*EJxIm, EIzImEJxRe = EIzIm*EJxRe, EIxReEJzIm = EIxRe*EJzIm, EIxImEJzRe = EIxIm*EJzRe;
				double EIzImEJxIm = EIzIm*EJxIm, EIzReEJxRe = EIzRe*EJxRe, EIxImEJzIm = EIxIm*EJzIm, EIxReEJzRe = EIxRe*EJzRe;
				double s2 = ReMult*((-EIzReEJxIm + EIzImEJxRe - EIxReEJzIm + EIxImEJzRe)*PhTotIm 
							+(-EIzImEJxIm - EIzReEJxRe - EIxImEJzIm - EIxReEJzRe)*PhTotRe);
				double s3 = ReMult*((-EIzImEJxIm - EIzReEJxRe + EIxImEJzIm + EIxReEJzRe)*PhTotIm 
							+(EIzReEJxIm - EIzImEJxRe - EIxReEJzIm + EIxImEJzRe)*PhTotRe);
				
				Stokes.s0 += s0; Stokes.s1 += s1; Stokes.s2 += s2; Stokes.s3 += s3;
			}
			tj++;
		}
		ti++;
	}
	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeCrossTermsAtPointPer(int StartIntervNo, int FinIntervNo, srTStokes& Stokes)
{
	//const double RelTolIgnoreComp = 1.E-07; // To steer
	const double MinArgIgnoreComp = -15.; // To steer

	int AmOfInterv = FinIntervNo - StartIntervNo + 1;
	srTRadContribInterval *pStartInterv = RadNumIntervals + StartIntervNo;

	double PIdLamb_Inv_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 2.533840802E+06*EXZ.e : 3.14159265359E+09/EXZ.e;
	double TwoPIdLamb_Inv_m = 2.*PIdLamb_Inv_m;
	double yObs = DistrInfoDat.yStart;

	double yObsInv = 1./yObs;
	double yObsInvE2 = yObsInv*yObsInv;

	srTEbmDat& Ebm = TrjDatPtr->EbmDat;
	double EnSprMult = 2.*PIdLamb_Inv_m*PIdLamb_Inv_m*(Ebm.GammaEm2)*(Ebm.GammaEm2)*(Ebm.SigmaRelE)*(Ebm.SigmaRelE);
	double PerLen = 0.5*FieldBasedArrays.sStep*(FieldBasedArrays.Ns - 1);

	//int TotAmOfInt = FieldBasedArrays.Nper*AmOfInterv;
	long long TotAmOfInt = FieldBasedArrays.Nper*AmOfInterv;
	char NearField = 0;
	double DelPhi = 0.5*MaxPhaseDifForThisObsPoint(NearField);

	int IntCountI = 0, PerCountI = 0;
	srTRadContribInterval *ti = pStartInterv;
	//for(int i=0; i<TotAmOfInt; i++)
	for(long long i=0; i<TotAmOfInt; i++)
	{
		if(IntCountI == AmOfInterv) { IntCountI = 0; ti = pStartInterv; PerCountI++;}
		double si = ti->sc + (PerCountI - FieldBasedArrays.NperLeft)*PerLen;

		srTEFourier& EI = ti->E;
		double EIxRe = EI.EwX_Re, EIxIm = EI.EwX_Im;
		double EIzRe = EI.EwZ_Re, EIzIm = EI.EwZ_Im;

		int IntCountJ = IntCountI + 1, PerCountJ = PerCountI;
		srTRadContribInterval *tj = ti + 1;
		//for(int j=(i + 1); j<TotAmOfInt; j++)
		for(long long j=(i + 1); j<TotAmOfInt; j++)
		{
			if(IntCountJ >= AmOfInterv) { IntCountJ = 0; tj = pStartInterv; PerCountJ++;}

			double sj = tj->sc + (PerCountJ - FieldBasedArrays.NperLeft)*PerLen;
			double si_mi_sj = si - sj;

			srTEFourier& EJ = tj->E;
			double EJxRe = EJ.EwX_Re, EJxIm = EJ.EwX_Im;
			double EJzRe = EJ.EwZ_Re, EJzIm = EJ.EwZ_Im;

			double DelPhiPer = DelPhi*(PerCountI - PerCountJ);

			double T2 = PIdLamb_Inv_m*si_mi_sj*yObsInvE2;
			double T2e2 = T2*T2, T2d4 = 0.25*T2;
			double T2e2pGxE2_Inv = 1./(T2e2 + Gx*Gx);
			double T2e2pGzE2_Inv = 1./(T2e2 + Gz*Gz);

			double T1x = TwoPIdLamb_Inv_m*yObsInvE2*(-EXZ.x*si_mi_sj + ti->Xc*(yObs - sj) - tj->Xc*(yObs - si));
			double T1xe2 = T1x*T1x;
			double T1xe2T2e2pGxE2_Inv = T1xe2*T2e2pGxE2_Inv;

			double T1z = TwoPIdLamb_Inv_m*yObsInvE2*(-EXZ.z*si_mi_sj + ti->Zc*(yObs - sj) - tj->Zc*(yObs - si));
			double T1ze2 = T1z*T1z;
			double T1ze2T2e2pGzE2_Inv = T1ze2*T2e2pGzE2_Inv;

			double ArgTot = -0.25*(Gx*T1xe2T2e2pGxE2_Inv + Gz*T1ze2T2e2pGzE2_Inv) - EnSprMult*si_mi_sj*si_mi_sj;
			if(ArgTot > MinArgIgnoreComp)
			{
				double PhTot = 0.5*(atan(T2/Gx) + atan(T2/Gz)) - T2d4*(T1xe2T2e2pGxE2_Inv + T1ze2T2e2pGzE2_Inv) + DelPhiPer;
				double PhTotRe, PhTotIm; CosAndSinComp.CosAndSin(PhTot, PhTotRe, PhTotIm);
				
				double nCon = 2.*sqrt(Gx*Gz*sqrt(T2e2pGxE2_Inv*T2e2pGzE2_Inv));
				double ReMult = nCon*exp(ArgTot);
				
				double EIxReEJxIm = EIxRe*EJxIm, EIxImEJxRe = EIxIm*EJxRe, EIzReEJzIm = EIzRe*EJzIm, EIzImEJzRe = EIzIm*EJzRe;
				double EIxImEJxIm = EIxIm*EJxIm, EIxReEJxRe = EIxRe*EJxRe, EIzImEJzIm = EIzIm*EJzIm, EIzReEJzRe = EIzRe*EJzRe;
				double s0 = ReMult*((EIxReEJxIm - EIxImEJxRe + EIzReEJzIm - EIzImEJzRe)*PhTotIm 
							+(EIxImEJxIm + EIxReEJxRe + EIzImEJzIm + EIzReEJzRe)*PhTotRe);
				double s1 = ReMult*((EIxReEJxIm - EIxImEJxRe - EIzReEJzIm + EIzImEJzRe)*PhTotIm 
							+(EIxImEJxIm + EIxReEJxRe - EIzImEJzIm - EIzReEJzRe)*PhTotRe);
				double EIzReEJxIm = EIzRe*EJxIm, EIzImEJxRe = EIzIm*EJxRe, EIxReEJzIm = EIxRe*EJzIm, EIxImEJzRe = EIxIm*EJzRe;
				double EIzImEJxIm = EIzIm*EJxIm, EIzReEJxRe = EIzRe*EJxRe, EIxImEJzIm = EIxIm*EJzIm, EIxReEJzRe = EIxRe*EJzRe;
				double s2 = ReMult*((-EIzReEJxIm + EIzImEJxRe - EIxReEJzIm + EIxImEJzRe)*PhTotIm 
							+(-EIzImEJxIm - EIzReEJxRe - EIxImEJzIm - EIxReEJzRe)*PhTotRe);
				double s3 = ReMult*((-EIzImEJxIm - EIzReEJxRe + EIxImEJzIm + EIxReEJzRe)*PhTotIm 
							+(EIzReEJxIm - EIzImEJxRe - EIxReEJzIm + EIxImEJzRe)*PhTotRe);
				
				Stokes.s0 += s0; Stokes.s1 += s1; Stokes.s2 += s2; Stokes.s3 += s3;
			}
			tj++; IntCountJ++;
		}
		ti++; IntCountI++;
	}
	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeAndAddOneMainTermNum(int TermNo)
{
	int result;
	const long NpOnLevelMaxNoResult = 5000000; // To steer; to stop computation as unsuccessful

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	srTRadContribInterval& ThisInterv = RadNumIntervals[TermNo];
	//int IndSt = ThisInterv.IndSt, IndFi = ThisInterv.IndFi;
	long long IndSt = ThisInterv.IndSt, IndFi = ThisInterv.IndFi;
	double sStart = FieldBasedArrays.sStart + IndSt*FieldBasedArrays.sStep;
	double sEnd = FieldBasedArrays.sStart + IndFi*FieldBasedArrays.sStep;

	srTEFourier& E = ThisInterv.E;

	//int AmOfPoOnZeroLev = IndFi - IndSt + 1;
	long long AmOfPoOnZeroLev = IndFi - IndSt + 1;
	TrjArraysAux.Reset(AmOfPoOnZeroLev);

	char NearField = (LongIntTypeG == 1);

	//int Np = AmOfPoOnZeroLev - 1;
	long long Np = AmOfPoOnZeroLev - 1;
	double sStep = (sEnd - sStart)/double(Np), s;

	double CosPh, SinPh; CosAndSinComp.CosAndSin(ThisInterv.Ph_St, CosPh, SinPh);
	double ReExSt = ThisInterv.Ax_St*CosPh, ImExSt = ThisInterv.Ax_St*SinPh;
	double ReEzSt = ThisInterv.Az_St*CosPh, ImEzSt = ThisInterv.Az_St*SinPh;

	double AxdPhds = ThisInterv.Ax_St*ThisInterv.dPhds_St;
	double DReExSt = ThisInterv.dAxds_St*CosPh - AxdPhds*SinPh;
	double DImExSt = AxdPhds*CosPh + ThisInterv.dAxds_St*SinPh;
	double AzdPhds = ThisInterv.Az_St*ThisInterv.dPhds_St;
	double DReEzSt = ThisInterv.dAzds_St*CosPh - AzdPhds*SinPh;
	double DImEzSt = AzdPhds*CosPh + ThisInterv.dAzds_St*SinPh;

	CosAndSinComp.CosAndSin(ThisInterv.Ph_Fi, CosPh, SinPh);
	double ReExFi = ThisInterv.Ax_Fi*CosPh, ImExFi = ThisInterv.Ax_Fi*SinPh;
	double ReEzFi = ThisInterv.Az_Fi*CosPh, ImEzFi = ThisInterv.Az_Fi*SinPh;

	AxdPhds = ThisInterv.Ax_Fi*ThisInterv.dPhds_Fi;
	double DReExFi = ThisInterv.dAxds_Fi*CosPh - AxdPhds*SinPh;
	double DImExFi = AxdPhds*CosPh + ThisInterv.dAxds_Fi*SinPh;

	AzdPhds = ThisInterv.Az_Fi*ThisInterv.dPhds_Fi;
	double DReEzFi = ThisInterv.dAzds_Fi*CosPh - AzdPhds*SinPh;
	double DImEzFi = AzdPhds*CosPh + ThisInterv.dAzds_Fi*SinPh;

	double wFxRe = wfe*(ReExSt + ReExFi);
	double wFxIm = wfe*(ImExSt + ImExFi);
	double wFzRe = wfe*(ReEzSt + ReEzFi);
	double wFzIm = wfe*(ImEzSt + ImEzFi);
	double wdFxRe = wd*(DReExSt - DReExFi);
	double wdFxIm = wd*(DImExSt - DImExFi);
	double wdFzRe = wd*(DReEzSt - DReEzFi);
	double wdFzIm = wd*(DImEzSt - DImEzFi);

	double PhInit = ThisInterv.Ph_St, PhPrev;

	double Sum1XRe = ThisInterv.SumReEx, Sum1XIm = ThisInterv.SumImEx;
	double Sum1ZRe = ThisInterv.SumReEz, Sum1ZIm = ThisInterv.SumImEz;
	double Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;

	double ActNormConst, PIm10e9_d_Lamb;
	if(DistrInfoDat.TreatLambdaAsEnergyIn_eV)
	{
		ActNormConst = NormalizingConst*EXZ.e*0.80654658E-03; // Assumes EXZ.e in eV
		PIm10e9_d_Lamb = PIm10e6dEnCon*EXZ.e;
	}
	else
	{
		ActNormConst = NormalizingConst/EXZ.e; // Assumes EXZ.e in nm
		PIm10e9_d_Lamb = PIm10e6*1000./EXZ.e;
	}
	double ActNormConstsStep = ActNormConst*sStep;

	const double wfeTrap = 0.5/wfe;
	double IntXRe = E.EwX_Re + ActNormConstsStep*(wfeTrap*wFxRe + Sum1XRe);
	double IntXIm = E.EwX_Im + ActNormConstsStep*(wfeTrap*wFxIm + Sum1XIm);
	double IntZRe = E.EwZ_Re + ActNormConstsStep*(wfeTrap*wFzRe + Sum1ZRe);
	double IntZIm = E.EwZ_Im + ActNormConstsStep*(wfeTrap*wFzIm + Sum1ZIm);
	double SqNorm = IntXRe*IntXRe + IntXIm*IntXIm + IntZRe*IntZRe + IntZIm*IntZIm;

	double xObs = EXZ.x, zObs = EXZ.z, yObs = DistrInfoDat.yStart;
	double GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz, Ph, Ax, Az;

	double AngPhConst, Inv_yObs, Two_Nx, Two_Nz;
	if(!NearField)
	{
		Inv_yObs = 1./yObs;
		Nx = xObs*Inv_yObs; Nz = zObs*Inv_yObs;
		Two_Nx = 2.*Nx; Two_Nz = 2.*Nz;
		AngPhConst = GmEm2 + Nx*Nx + Nz*Nz;
	}

	double *pBtx, *pX, *pIntBtxE2, *pBx, *pBtz, *pZ, *pIntBtzE2, *pBz;
	double **TrjPtrs[] = {&pBtx, &pX, &pIntBtxE2, &pBx, &pBtz, &pZ, &pIntBtzE2, &pBz};

	//char ExtraPassForAnyCase = 0, NotFinishedYet = 1;
	char NotFinishedYet = 1;
	int LevelNo = 0;

	while(NotFinishedYet)
	{
		Sum2XRe += Sum1XRe; Sum2XIm += Sum1XIm; Sum2ZRe += Sum1ZRe; Sum2ZIm += Sum1ZIm; 
		Sum1XRe = Sum1XIm = Sum1ZRe = Sum1ZIm = 0.;
		char ThisMayBeTheLastLoop = 1;

		PhPrev = PhInit; LevelNo++;
		double HalfStep = 0.5*sStep;
		s = sStart + HalfStep;

		if(result = TrjArraysAux.FillNextAndDestroyPrevLevel(LevelNo, s, sEnd - HalfStep, TrjDatPtr)) return result;
		TrjArraysAux.SetupPtrs_Btx_X_IntBtxE2_Bx_Btz_Z_IntBtzE2_Bz(LevelNo, TrjPtrs);

			double DPhMax = 0.;

		//for(int i=0; i<Np; i++)
		for(long long i=0; i<Np; i++)
		{
			if(NearField)
			{
				One_d_ymis = 1./(yObs - s);
				xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
				Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
				Ph = PIm10e9_d_Lamb*(s*GmEm2 + *pIntBtxE2 + *pIntBtzE2 + xObs_mi_x*Nx + zObs_mi_z*Nz);
				Ax = (*pBtx - Nx)*One_d_ymis; Az = (*pBtz - Nz)*One_d_ymis;
			}
			else
			{
				Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_Nx*(*pX) + Two_Nz*(*pZ)));
				Ax = (*pBtx - Nx)*Inv_yObs; Az = (*pBtz - Nz)*Inv_yObs;
			}

			pBtx++; pX++; pIntBtxE2++; pBtz++; pZ++; pIntBtzE2++;

			CosAndSinComp.CosAndSin(Ph, CosPh, SinPh);

			Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh;
			s += sStep;

			double dPh = ::fabs(Ph - PhPrev);
			if(dPh > FivePIdFour) ThisMayBeTheLastLoop = 0;

				if(dPh > DPhMax) DPhMax = dPh;

			PhPrev = Ph;
		}
		double ActNormConstHalfStep = ActNormConst*HalfStep;
		double LocIntXRe = E.EwX_Re + ActNormConstHalfStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + HalfStep*wdFxRe);
		double LocIntXIm = E.EwX_Im + ActNormConstHalfStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + HalfStep*wdFxIm);
		double LocIntZRe = E.EwZ_Re + ActNormConstHalfStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + HalfStep*wdFzRe);
		double LocIntZIm = E.EwZ_Im + ActNormConstHalfStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + HalfStep*wdFzIm);
		double LocSqNorm = LocIntXRe*LocIntXRe + LocIntXIm*LocIntXIm + LocIntZRe*LocIntZRe + LocIntZIm*LocIntZIm;

		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			//char SharplyGoesDown = (LocSqNorm < 0.2*SqNorm); // To steer

			char NotFinishedYetFirstTest;
			if(MaxIntValWithinPointG > 0.) NotFinishedYetFirstTest = (TestVal > sIntegRelPrecG*MaxIntValWithinPointG);
			else NotFinishedYetFirstTest = (TestVal > sIntegRelPrecG*LocSqNorm);

			if(!NotFinishedYetFirstTest)
			{
				NotFinishedYet = 0;
			}
		}

		TrjArraysAux.DeallocateArraysOnLevel(LevelNo);

		IntXRe = LocIntXRe; IntXIm = LocIntXIm; IntZRe = LocIntZRe; IntZIm = LocIntZIm; SqNorm = LocSqNorm;
		sStep = HalfStep; Np <<= 1;

		if(NotFinishedYet)
		{
			if(Np > NpOnLevelMaxNoResult) return CAN_NOT_COMPUTE_RADIATION_INTEGRAL;
		}
	}
	E.EwX_Re = IntXRe; E.EwX_Im = IntXIm; E.EwZ_Re = IntZRe; E.EwZ_Im = IntZIm;

	if(MaxIntValWithinPointG < SqNorm) 
	{
		MaxIntValWithinPointG = SqNorm;
	}
	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeAndAddTerminations(int TermNo)
{
	srTRadContribInterval& Interv = RadNumIntervals[TermNo];
	srTEFourier LeftResE, RightResE;

	srTRadFuncsAndDers LeftFuncs, RightFuncs;
	Interv.OutLeftFuncs(LeftFuncs);
	LeftFuncs.s = FieldBasedArrays.sStart + Interv.IndSt*FieldBasedArrays.sStep;
	Interv.OutRightFuncs(RightFuncs);
	RightFuncs.s = FieldBasedArrays.sStart + Interv.IndFi*FieldBasedArrays.sStep;

	int result;
	int NumberOfTerms = 2;

	char LeftResidShouldNotBeComp = !(Interv.PoIsGoodForAn_St);
	if(!LeftResidShouldNotBeComp)
	{
		if(result = ComputeNormalResidual(LeftFuncs, NumberOfTerms, LeftResE)) return result;
	}

	char RightResidShouldNotBeComp = !(Interv.PoIsGoodForAn_Fi);
	if(!RightResidShouldNotBeComp)
	{
		if(result = ComputeNormalResidual(RightFuncs, NumberOfTerms, RightResE)) return result;
	}
	Interv.E += (RightResE - LeftResE);
	return 0;
}

//*************************************************************************

int srTRadIntWiggler::ComputeNormalResidual(srTRadFuncsAndDers& Funcs, int NumberOfTerms, srTEFourier& Res)
{
	if((NumberOfTerms < 1) || (NumberOfTerms > 2)) NumberOfTerms = 2;

	double ActNormConst = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? NormalizingConst*EXZ.e*0.80654658E-03 : NormalizingConst/EXZ.e;

	double PreExpXRe = Funcs.Ax, PreExpXIm = 0.;
	double PreExpZRe = Funcs.Az, PreExpZIm = 0.;

	double One_d_dPhds = 1./Funcs.dPhds;
	if(NumberOfTerms > 1)
	{
		double d2Phds2_d_dPhds = Funcs.d2Phds2*One_d_dPhds;
		
		double d2Phds2_d_dPhdsE2_ForExpansionTest = d2Phds2_d_dPhds*One_d_dPhds;
		char FarFieldExpansDoesNotWork = (::fabs(d2Phds2_d_dPhdsE2_ForExpansionTest) > RelTolForAnTermsG);

		if(!FarFieldExpansDoesNotWork)
		{
			double t2xd = (Funcs.dAxds - Funcs.Ax*d2Phds2_d_dPhds)*One_d_dPhds;
			double t2zd = (Funcs.dAzds - Funcs.Az*d2Phds2_d_dPhds)*One_d_dPhds;

			PreExpXIm += t2xd;
			PreExpZIm += t2zd;
		}
	}

	double BufPreExpXRe = -PreExpXIm*One_d_dPhds;
	PreExpXIm = PreExpXRe*One_d_dPhds;
	PreExpXRe = BufPreExpXRe;
	double BufPreExpZRe = -PreExpZIm*One_d_dPhds;
	PreExpZIm = PreExpZRe*One_d_dPhds;
	PreExpZRe = BufPreExpZRe;

	Res.EwX_Re = ActNormConst*(PreExpXRe*Funcs.CosPh - PreExpXIm*Funcs.SinPh);
	Res.EwX_Im = ActNormConst*(PreExpXRe*Funcs.SinPh + PreExpXIm*Funcs.CosPh);
	Res.EwZ_Re = ActNormConst*(PreExpZRe*Funcs.CosPh - PreExpZIm*Funcs.SinPh);
	Res.EwZ_Im = ActNormConst*(PreExpZRe*Funcs.SinPh + PreExpZIm*Funcs.CosPh);
	return 0;
}

//*************************************************************************

void srTRadIntWiggler::DetermineRadIntervals(int& AmOfInterv)
{
	srTRadContribInterval *tNumInt = RadNumIntervals;
	for(int i=0; i<FieldBasedArrays.Ns; i++)
	{
		tNumInt->IndSt = tNumInt->IndFi = -1;
		tNumInt->PoIsGoodForAn_St = tNumInt->PoIsGoodForAn_Fi = 0;
		
		tNumInt++;
	}

	int AmOfMax, AmOfMin, AmOfNumInterv;
	AnalizeAllRadPoints(AmOfMax, AmOfMin, AmOfNumInterv);

	//OC140210
	RadNumIntervals[0].PoIsGoodForAn_St = 1;
	if(AmOfNumInterv > 0) RadNumIntervals[AmOfNumInterv - 1].PoIsGoodForAn_Fi = 1; //OC290410
	else RadNumIntervals[0].PoIsGoodForAn_Fi = 1;

	if(NumIntegG) 
	{ 
		AmOfInterv = AmOfNumInterv;
		SetupCentralValuesInNumInterv(AmOfInterv);
		return;
	}
	else SetupIntervalsForMcDonaldEstimation(AmOfMax, AmOfMin, AmOfInterv);
}

//*************************************************************************

void srTRadIntWiggler::AnalizeAllRadPoints(int& AmOfMax, int& AmOfMin, int& AmOfNumIntervals)
{
	double PIm10e9_d_Lamb = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*EXZ.e : PIm10e6*1000./EXZ.e;
	double TwoPIm10e9_d_Lamb = 2.*PIm10e9_d_Lamb;

	double yObs = DistrInfoDat.yStart, GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double ConBtx = TrjDatPtr->BetaNormConst, ConBtz = -TrjDatPtr->BetaNormConst;
	double Btx_mi_Nx, Btz_mi_Nz, One_d_ymis;
	double Btx_mi_Nx_Pr, Btz_mi_Nz_Pr;
	double BetaFact;

	double Nx, Nz;

	char NearField = (LongIntTypeG == 1);
	double AngPhConst, Inv_yObs, Two_Nx, Two_Nz;
	if(!NearField)
	{
		Inv_yObs = 1./yObs;
		Nx = EXZ.x*Inv_yObs; Nz = EXZ.z*Inv_yObs;
		Two_Nx = 2.*Nx; Two_Nz = 2.*Nz;
		AngPhConst = GmEm2 + Nx*Nx + Nz*Nz;
	}

	double Ph, dPhds, d2Phds2, Ax, Az, dAxds, dAzds;
	double Ph_Pr, dPhds_Pr, d2Phds2_Pr, Ax_Pr, Az_Pr, dAxds_Pr, dAzds_Pr;

	double GamFact = 1.*GmEm2; // To steer

	double s = FieldBasedArrays.sStart;
	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	double *pX = FieldBasedArrays.XArr, *pZ = FieldBasedArrays.ZArr;
	double *pIntBtxE2 = FieldBasedArrays.IntBtxE2Arr, *pIntBtzE2 = FieldBasedArrays.IntBtzE2Arr;
	int *tMaxPoInd = MaxPoIndArray, *tMinPoInd = MinPoIndArray;
	char DerSign;
	double Prev_d2Phds2ddPhdsE2 = 0.;
	AmOfMax = AmOfMin = 0;

	srTRadContribInterval *tNumInterval = RadNumIntervals;
	AmOfNumIntervals = 0;
	char PrevPointIsGoodForAn, PointIsGoodForAn;

	char LastIntervalCutArtificially = 1, ArtificialCutPassedButNotApplied = 0;
	char MainFieldPlane = FindMaxFieldPlane(FieldBasedArrays);
	//srTSend Send;

	double &rAbsBxMax = FieldBasedArrays.absBxMax;
	double &rAbsBzMax = FieldBasedArrays.absBzMax;

	for(int i=0; i<FieldBasedArrays.Ns; i++)
	{
		double ConBtxBz = ConBtx*(*pBz), ConBtzBx = ConBtz*(*pBx);

		if(NearField)
		{
			One_d_ymis = 1./(yObs - s); 
			double xObs_mi_x = EXZ.x - *pX, zObs_mi_z = EXZ.z - *pZ;
			Nx = xObs_mi_x*One_d_ymis; Nz = zObs_mi_z*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*GmEm2 + *pIntBtxE2 + *pIntBtzE2 + xObs_mi_x*Nx + zObs_mi_z*Nz);

			Btx_mi_Nx = *pBtx - Nx; Btz_mi_Nz = *pBtz - Nz;
			Ax = Btx_mi_Nx*One_d_ymis; Az = Btz_mi_Nz*One_d_ymis;
			dAxds = (2.*Ax + ConBtxBz)*One_d_ymis; dAzds = (2.*Az + ConBtzBx)*One_d_ymis;
			
			BetaFact = Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz;
			dPhds = PIm10e9_d_Lamb*(GmEm2 + BetaFact);
			d2Phds2 = TwoPIm10e9_d_Lamb*(Btx_mi_Nx*(ConBtxBz + Ax) + Btz_mi_Nz*(ConBtzBx + Az));
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_Nx*(*pX) + Two_Nz*(*pZ)));
			Btx_mi_Nx = *pBtx - Nx; Btz_mi_Nz = *pBtz - Nz;
			Ax = Btx_mi_Nx*Inv_yObs; Az = Btz_mi_Nz*Inv_yObs;
			dAxds = ConBtxBz*Inv_yObs; dAzds = ConBtzBx*Inv_yObs;
			
			BetaFact = Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz;
			dPhds = PIm10e9_d_Lamb*(GmEm2 + BetaFact);
			d2Phds2 = TwoPIm10e9_d_Lamb*(Btx_mi_Nx*ConBtxBz + Btz_mi_Nz*ConBtzBx);
		}

		double d2Phds2ddPhdsE2 = d2Phds2/(dPhds*dPhds);
		//PointIsGoodForAn = (::fabs(d2Phds2ddPhdsE2) < RelTolForAnTermsG) 
		//				&& ((BetaFact > GamFact) || ((*pBx == 0.) && (*pBz == 0.)));
		PointIsGoodForAn = (::fabs(d2Phds2ddPhdsE2) < RelTolForAnTermsG) 
			&& ((BetaFact > GamFact) || ((::fabs(*pBx) <= rAbsBxMax*sIntegRelPrecG) && (::fabs(*pBz) <= rAbsBzMax*sIntegRelPrecG)));

		//if((!(TrjDatPtr->HorFieldIsNotZero)) && (Btx_mi_Nx*Btx_mi_Nx < GamFact)) PointIsGoodForAn = 0; //OC140809 commented-out
		//if((!(TrjDatPtr->VerFieldIsNotZero)) && (Btz_mi_Nz*Btz_mi_Nz < GamFact)) PointIsGoodForAn = 0; //OC140809 commented-out
		
		char NewDerSign = (d2Phds2ddPhdsE2 < Prev_d2Phds2ddPhdsE2)? -1 : 1;
		if(i > 0)
		{
			char SumsAdded = 0;

			if(NewDerSign < DerSign)
			{
				*(tMaxPoInd++) = i - 1; AmOfMax++;
			}
			else if(NewDerSign > DerSign)
			{
				*(tMinPoInd++) = i - 1; AmOfMin++;
			}

			if(PrevPointIsGoodForAn && (!PointIsGoodForAn))
			{
				tNumInterval->IndSt = i - 1;
				tNumInterval->Ph_St = Ph_Pr; tNumInterval->dPhds_St = dPhds_Pr; tNumInterval->d2Phds2_St = d2Phds2_Pr;
				tNumInterval->Ax_St = Ax_Pr; tNumInterval->Az_St = Az_Pr;
				tNumInterval->dAxds_St = dAxds_Pr; tNumInterval->dAzds_St = dAzds_Pr;
				CosAndSinComp.CosAndSin(tNumInterval->Ph_St, tNumInterval->CosPh_St, tNumInterval->SinPh_St);

				tNumInterval->PoIsGoodForAn_St = 1;

				double CosPh, SinPh; CosAndSinComp.CosAndSin(Ph, CosPh, SinPh);
				tNumInterval->SumReEx = Ax*CosPh; tNumInterval->SumImEx = Ax*SinPh;
				tNumInterval->SumReEz = Az*CosPh; tNumInterval->SumImEz = Az*SinPh;

				tNumInterval->Bx_St = *(pBx - 1); tNumInterval->Bz_St = *(pBz - 1);
				tNumInterval->Btx_St = *(pBtx - 1); tNumInterval->Btz_St = *(pBtz - 1);
				tNumInterval->X_St = *(pX - 1); tNumInterval->Z_St = *(pZ - 1);
				tNumInterval->IntBtxE2_St = *(pIntBtxE2 - 1); tNumInterval->IntBtzE2_St = *(pIntBtzE2 - 1);

				LastIntervalCutArtificially = 0; ArtificialCutPassedButNotApplied = 0;
			}
			else if((!PrevPointIsGoodForAn) && PointIsGoodForAn)
			{
				tNumInterval->IndFi = i; 
				tNumInterval->Ph_Fi = Ph; tNumInterval->dPhds_Fi = dPhds; tNumInterval->d2Phds2_Fi = d2Phds2;
				tNumInterval->Ax_Fi = Ax; tNumInterval->Az_Fi = Az;
				tNumInterval->dAxds_Fi = dAxds; tNumInterval->dAzds_Fi = dAzds;
				CosAndSinComp.CosAndSin(tNumInterval->Ph_Fi, tNumInterval->CosPh_Fi, tNumInterval->SinPh_Fi);
				
				tNumInterval->Bx_Fi = *pBx; tNumInterval->Bz_Fi = *pBz;
				tNumInterval->Btx_Fi = *pBtx; tNumInterval->Btz_Fi = *pBtz;
				tNumInterval->X_Fi = *pX; tNumInterval->Z_Fi = *pZ;
				tNumInterval->IntBtxE2_Fi = *pIntBtxE2; tNumInterval->IntBtzE2_Fi = *pIntBtzE2;
				
				tNumInterval->PoIsGoodForAn_Fi = 1;
				tNumInterval++; AmOfNumIntervals++;
				
				LastIntervalCutArtificially = 0; ArtificialCutPassedButNotApplied = 0;
			}
			else if(!PointIsGoodForAn)
			{
				double CosPh, SinPh; CosAndSinComp.CosAndSin(Ph, CosPh, SinPh);
				tNumInterval->SumReEx += Ax*CosPh; tNumInterval->SumImEx += Ax*SinPh;
				tNumInterval->SumReEz += Az*CosPh; tNumInterval->SumImEz += Az*SinPh;

				SumsAdded = 1;
			}

			if((i < FieldBasedArrays.Ns - 2) && (!PointIsGoodForAn))
			{
				char FieldPassedUp=0;

				if(MainFieldPlane == 'z') FieldPassedUp = (*(pBz-1))*(*pBz) < 0.;
				else FieldPassedUp = (*(pBx-1))*(*pBx) < 0.;

				if(FieldPassedUp)
				{
					if(LastIntervalCutArtificially || ArtificialCutPassedButNotApplied)
					{// Cut interval artificially
						double CosPh, SinPh; CosAndSinComp.CosAndSin(Ph, CosPh, SinPh);
						if(SumsAdded)
						{
							tNumInterval->SumReEx -= Ax*CosPh; tNumInterval->SumImEx -= Ax*SinPh;
							tNumInterval->SumReEz -= Az*CosPh; tNumInterval->SumImEz -= Az*SinPh;
						}
						
						tNumInterval->IndFi = i; 
						tNumInterval->Ph_Fi = Ph; tNumInterval->dPhds_Fi = dPhds; tNumInterval->d2Phds2_Fi = d2Phds2;
						tNumInterval->Ax_Fi = Ax; tNumInterval->Az_Fi = Az;
						tNumInterval->dAxds_Fi = dAxds; tNumInterval->dAzds_Fi = dAzds;
						tNumInterval->CosPh_Fi = CosPh; tNumInterval->SinPh_Fi = SinPh;
						
						tNumInterval->Bx_Fi = *pBx; tNumInterval->Bz_Fi = *pBz;
						tNumInterval->Btx_Fi = *pBtx; tNumInterval->Btz_Fi = *pBtz;
						tNumInterval->X_Fi = *pX; tNumInterval->Z_Fi = *pZ;
						tNumInterval->IntBtxE2_Fi = *pIntBtxE2; tNumInterval->IntBtzE2_Fi = *pIntBtzE2;
						
						tNumInterval->PoIsGoodForAn_Fi = 0;
						tNumInterval++; AmOfNumIntervals++;
						
						tNumInterval->IndSt = i;
						tNumInterval->Ph_St = Ph; tNumInterval->dPhds_St = dPhds; tNumInterval->d2Phds2_St = d2Phds2;
						tNumInterval->Ax_St = Ax; tNumInterval->Az_St = Az;
						tNumInterval->dAxds_St = dAxds; tNumInterval->dAzds_St = dAzds;
						tNumInterval->CosPh_St = CosPh; tNumInterval->SinPh_St = SinPh;
						
						tNumInterval->SumReEx = 0.; tNumInterval->SumImEx = 0.;
						tNumInterval->SumReEz = 0.; tNumInterval->SumImEz = 0.;
						
						tNumInterval->Bx_St = *pBx; tNumInterval->Bz_St = *pBz;
						tNumInterval->Btx_St = *pBtx; tNumInterval->Btz_St = *pBtz;
						tNumInterval->X_St = *pX; tNumInterval->Z_St = *pZ;
						tNumInterval->IntBtxE2_St = *pIntBtxE2; tNumInterval->IntBtzE2_St = *pIntBtzE2;
						
						tNumInterval->PoIsGoodForAn_St = 0;
						
						LastIntervalCutArtificially = 1; ArtificialCutPassedButNotApplied = 0;
						//Send.AddWarningMessage(&gVectWarnNos, NOT_A_WIGGLER_CASE);
						//CErrWarn::AddWarningMessage(&gVectWarnNos, NOT_A_WIGGLER_CASE);
						//OC commented out 160304
					}
					else 
					{
						ArtificialCutPassedButNotApplied = 1;
					}
				}
			}
		}
		else
		{
			DerSign = NewDerSign;

			if(!PointIsGoodForAn)
			{
				tNumInterval->IndSt = 0;

				tNumInterval->Ph_St = Ph; tNumInterval->dPhds_St = dPhds; tNumInterval->d2Phds2_St = d2Phds2;
				tNumInterval->Ax_St = Ax; tNumInterval->Az_St = Az;
				tNumInterval->dAxds_St = dAxds; tNumInterval->dAzds_St = dAzds;
				CosAndSinComp.CosAndSin(tNumInterval->Ph_St, tNumInterval->CosPh_St, tNumInterval->SinPh_St);

				tNumInterval->PoIsGoodForAn_St = 0;

				tNumInterval->SumReEx = 0.; tNumInterval->SumImEx = 0.;
				tNumInterval->SumReEz = 0.; tNumInterval->SumImEz = 0.;

				tNumInterval->Bx_St = *pBx; tNumInterval->Bz_St = *pBz;
				tNumInterval->Btx_St = *pBtx; tNumInterval->Btz_St = *pBtz;
				tNumInterval->X_St = *pX; tNumInterval->Z_St = *pZ;
				tNumInterval->IntBtxE2_St = *pIntBtxE2; tNumInterval->IntBtzE2_St = *pIntBtzE2;
			}
		}

		PrevPointIsGoodForAn = PointIsGoodForAn;
		Ph_Pr = Ph; dPhds_Pr = dPhds; d2Phds2_Pr = d2Phds2; 
		Ax_Pr = Ax; Az_Pr = Az; dAxds_Pr = dAxds; dAzds_Pr = dAzds;

		Btx_mi_Nx_Pr = Btx_mi_Nx; Btz_mi_Nz_Pr = Btz_mi_Nz;

		DerSign = NewDerSign;
		Prev_d2Phds2ddPhdsE2 = d2Phds2ddPhdsE2;

		pBx++; pBz++; pBtx++, pBtz++, pX++, pZ++, pIntBtxE2++, pIntBtzE2++;
		s += FieldBasedArrays.sStep;
	}

	if((tNumInterval->IndSt >= 0) && (tNumInterval->IndFi < 0))
	{
		tNumInterval->IndFi = FieldBasedArrays.Ns - 1;

		tNumInterval->Ph_Fi = Ph; tNumInterval->dPhds_Fi = dPhds; tNumInterval->d2Phds2_Fi = d2Phds2;
		tNumInterval->Ax_Fi = Ax; tNumInterval->Az_Fi = Az;
		tNumInterval->dAxds_Fi = dAxds; tNumInterval->dAzds_Fi = dAzds;
		CosAndSinComp.CosAndSin(tNumInterval->Ph_Fi, tNumInterval->CosPh_Fi, tNumInterval->SinPh_Fi);

		tNumInterval->Bx_Fi = *(pBx - 1); tNumInterval->Bz_Fi = *(pBz - 1);
		tNumInterval->Btx_Fi = *(pBtx - 1); tNumInterval->Btz_Fi = *(pBtz - 1);
		tNumInterval->X_Fi = *(pX - 1); tNumInterval->Z_Fi = *(pZ - 1);
		tNumInterval->IntBtxE2_Fi = *(pIntBtxE2 - 1); tNumInterval->IntBtzE2_Fi = *(pIntBtzE2 - 1);

		tNumInterval->PoIsGoodForAn_Fi = 0;
		AmOfNumIntervals++;
	}
}

//*************************************************************************

void srTRadIntWiggler::SetupCentralValuesInNumInterv(int AmOfInterv)
{//Fills double sc, Xc, Zc;
	srTRadContribInterval *ti = RadNumIntervals;
	double Btx, IntBtE2x, Btz, IntBtE2z;
	for(int i=0; i<AmOfInterv; i++)
	{
		ti->sc = FieldBasedArrays.sStart + 0.5*(ti->IndSt + ti->IndFi)*FieldBasedArrays.sStep;
		TrjDatPtr->CompTrjDataDerivedAtPoint(ti->sc, Btx, ti->Xc, IntBtE2x, Btz, ti->Zc, IntBtE2z);
		ti++;
	}
}

//*************************************************************************

char srTRadIntWiggler::FindMaxFieldPlane(srTFieldBasedArrays& FieldBasedArrays)
{
	double MaxAbsBx = 0., MaxAbsBz = 0.;
	double *tBx = FieldBasedArrays.BxArr, *tBz = FieldBasedArrays.BzArr;
	//for(int i=0; i<FieldBasedArrays.Ns; i++)
	for(long long i=0; i<FieldBasedArrays.Ns; i++)
	{
		double AbsBx = ::fabs(*(tBx++)), AbsBz = ::fabs(*(tBz++));
		if(MaxAbsBx < AbsBx) MaxAbsBx = AbsBx;
		if(MaxAbsBz < AbsBz) MaxAbsBz = AbsBz;
	}
	return (MaxAbsBx > MaxAbsBz)? 'x' : 'z';
}

//*************************************************************************

void srTRadIntWiggler::MergeAdjacentIntervals(int& AmOfInterv)
{
	srTRadContribInterval *t = RadNumIntervals;
	srTRadContribInterval *tOld = RadNumIntervals + 1;
	int AmOfIntervOld = AmOfInterv;

	for(int i=1; i<AmOfIntervOld; i++)
	{
		if(t->IndFi == tOld->IndSt)
		{
			t->IndFi = tOld->IndFi;

			t->SumReEx += (t->Ax_Fi)*(t->CosPh_Fi) + tOld->SumReEx;
			t->SumImEx += (t->Ax_Fi)*(t->SinPh_Fi) + tOld->SumImEx;
			t->SumReEz += (t->Az_Fi)*(t->CosPh_Fi) + tOld->SumReEz;
			t->SumImEz += (t->Az_Fi)*(t->SinPh_Fi) + tOld->SumImEz;

			t->Ph_Fi = tOld->Ph_Fi; t->dPhds_Fi = tOld->dPhds_Fi; t->d2Phds2_Fi = tOld->d2Phds2_Fi;
			t->Ax_Fi = tOld->Ax_Fi; t->dAxds_Fi = tOld->dAxds_Fi; t->Az_Fi = tOld->Az_Fi; t->dAzds_Fi = tOld->dAzds_Fi;
			
			t->CosPh_Fi = tOld->CosPh_Fi; t->SinPh_Fi = tOld->SinPh_Fi;
			t->PoIsGoodForAn_Fi = tOld->PoIsGoodForAn_Fi;

			t->Bx_Fi = tOld->Bx_Fi; t->Bz_Fi = tOld->Bz_Fi;
			t->Btx_Fi = tOld->Btx_Fi; t->Btz_Fi = tOld->Btz_Fi;
			t->X_Fi = tOld->X_Fi; t->Z_Fi = tOld->Z_Fi;
			t->IntBtxE2_Fi = tOld->IntBtxE2_Fi; t->IntBtzE2_Fi = tOld->IntBtzE2_Fi;

			// Copy any new member of srTRadContribInterval !!!

			AmOfInterv--;
		}
		else
		{
			t++;
			if(AmOfInterv != AmOfIntervOld) *t = *tOld;
		}
		tOld++;
	}
}

//*************************************************************************

void srTRadIntWiggler::SetupIntervalsForMcDonaldEstimation(int AmOfMax, int AmOfMin, int& AmOfInterv)
{
	int *tSmallerExtr = 0, *tLargerExtr = 0, AmOfSmallerExtr, AmOfLargerExtr;
	if(AmOfMax < AmOfMin)
	{
		AmOfSmallerExtr = AmOfMax; AmOfLargerExtr = AmOfMin;
		tSmallerExtr = MaxPoIndArray; tLargerExtr = MinPoIndArray;
	}
	else
	{
		AmOfSmallerExtr = AmOfMin; AmOfLargerExtr = AmOfMax;
		tSmallerExtr = MinPoIndArray; tLargerExtr = MaxPoIndArray;
	}

	srTRadContribInterval *tInterval = RadContribIntervals;
	AmOfInterv = 0;

	double sStart = FieldBasedArrays.sStart, sStep = FieldBasedArrays.sStep;
	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	short LocHorFieldIsNotZero = TrjDatPtr->HorFieldIsNotZero, LocVerFieldIsNotZero = TrjDatPtr->VerFieldIsNotZero;

	double Dif1 = ::fabs((double)(*tLargerExtr - *tSmallerExtr));
	double Dif2 = ::fabs((double)(*(tLargerExtr + 1) - *tSmallerExtr));
	if(Dif1 > Dif2) 
	{
		int IndCen = *tLargerExtr;
		tInterval->sc = sStart + IndCen*sStep;
		if(LocHorFieldIsNotZero)
		{
			tInterval->Bx_St = *(pBx + IndCen);
			tInterval->Btz_St = *(pBtz + IndCen);
		}
		if(LocVerFieldIsNotZero)
		{
			tInterval->Bz_St = *(pBz + IndCen);
			tInterval->Btx_St = *(pBtx + IndCen);
		}
		tInterval++; AmOfInterv++;
		tLargerExtr++; AmOfLargerExtr--;
	}

	double GamFact = 3./(TrjDatPtr->EbmDat.Gamma);
	for(int i=0; i<AmOfSmallerExtr; i++)
	{
		char ConsiderLargerExtr = (i < AmOfLargerExtr);
		if(ConsiderLargerExtr)
		{
			char DoNotMergePoints = 0;
			char DifExtrLargerThanOneStep = (::fabs((double)(*tSmallerExtr - *tLargerExtr)) > 1);
			if(LocHorFieldIsNotZero)
			{
				double Btz1 = *(pBtz + (*tSmallerExtr));
				double Btz2 = *(pBtz + (*tLargerExtr));
				DoNotMergePoints = (DoNotMergePoints || ((::fabs(Btz1 - Btz2) > GamFact) && DifExtrLargerThanOneStep));
			}
			if(LocVerFieldIsNotZero)
			{
				double Btx1 = *(pBtx + (*tSmallerExtr));
				double Btx2 = *(pBtx + (*tLargerExtr));
				DoNotMergePoints = (DoNotMergePoints || ((::fabs(Btx1 - Btx2) > GamFact) && DifExtrLargerThanOneStep));
			}
			
			if(DoNotMergePoints)
			{
				int IndCen1 = *tSmallerExtr, IndCen2 = *tLargerExtr;
				srTRadContribInterval *tInterval1 = tInterval, *tInterval2 = tInterval + 1;
				
				tInterval1->sc = sStart + IndCen1*sStep;
				tInterval2->sc = sStart + IndCen2*sStep;
				if(LocHorFieldIsNotZero)
				{
					tInterval1->Bx_St = *(pBx + IndCen1);
					tInterval1->Btz_St = *(pBtz + IndCen1);
					tInterval2->Bx_St = *(pBx + IndCen2);
					tInterval2->Btz_St = *(pBtz + IndCen2);
				}
				if(LocVerFieldIsNotZero)
				{
					tInterval1->Bz_St = *(pBz + IndCen1);
					tInterval1->Btx_St = *(pBtx + IndCen1);
					tInterval2->Bz_St = *(pBz + IndCen2);
					tInterval2->Btx_St = *(pBtx + IndCen2);
				}
				tInterval += 2; AmOfInterv += 2;
			}
			else
			{
				int IndCen = (*tSmallerExtr + *tLargerExtr) >> 1;
				
				tInterval->sc = sStart + IndCen*sStep;
				if(LocHorFieldIsNotZero)
				{
					tInterval->Bx_St = *(pBx + IndCen);
					tInterval->Btz_St = *(pBtz + IndCen);
				}
				if(LocVerFieldIsNotZero)
				{
					tInterval->Bz_St = *(pBz + IndCen);
					tInterval->Btx_St = *(pBtx + IndCen);
				}
				tInterval++; AmOfInterv++;
			}
			tSmallerExtr++; tLargerExtr++;
		}
		else
		{
			int IndCen = *tSmallerExtr;

			tInterval->sc = sStart + IndCen*sStep;
			if(LocHorFieldIsNotZero)
			{
				tInterval->Bx_St = *(pBx + IndCen);
				tInterval->Btz_St = *(pBtz + IndCen);
			}
			if(LocVerFieldIsNotZero)
			{
				tInterval->Bz_St = *(pBz + IndCen);
				tInterval->Btx_St = *(pBtx + IndCen);
			}
			tInterval++; AmOfInterv++;
			tSmallerExtr++;
		}
	}
}

//*************************************************************************

int srTRadIntWiggler::TreatFiniteElecBeamEmittance(srTStokesStructAccessData& StokesAccessData, char MainOrCrossTerms, double ElBeamMomFact)
{
	int result;
	if((StokesAccessData.nx == 1) && (StokesAccessData.nz == 1)) return 0;

	//long LenFloatArr = (StokesAccessData.nx*StokesAccessData.nz);
	long long LenFloatArr = (((long long)StokesAccessData.nx)*((long long)StokesAccessData.nz));
	float *StokesCmpnArr = new float[LenFloatArr];
	if(StokesCmpnArr == 0) return MEMORY_ALLOCATION_FAILURE;

	for(int ie=0; ie<StokesAccessData.ne; ie++)
	{
		for(int is=0; is<4; is++)
		{
			ExtractStokesSliceConstE(StokesAccessData, ie, is, MainOrCrossTerms, StokesCmpnArr);

			if((StokesAccessData.nx == 1) || (StokesAccessData.nz == 1))
			{
				if(result = TreatFiniteElecBeamEmittanceOneComp1D(StokesCmpnArr, ElBeamMomFact)) return result;
			}
			else
			{
				if(result = TreatFiniteElecBeamEmittanceOneComp2D(StokesCmpnArr, ElBeamMomFact)) return result;
			}

			if((is == 0) && (MainOrCrossTerms == 'm')) SuppressNegativeValues(StokesCmpnArr);

			UpdateStokesSliceConstE(StokesCmpnArr, ie, is, MainOrCrossTerms, StokesAccessData);
		}
	}
	delete[] StokesCmpnArr;
	return 0;
}

//*************************************************************************

void srTRadIntWiggler::ExtractStokesSliceConstE(srTStokesStructAccessData& StokesAccessData, long ie, int StokesNo, char MainOrCrossTerms, float* pOutS)
{
	float *pS0;
	if(MainOrCrossTerms == 'm') pS0	= StokesAccessData.pBaseSto + StokesNo;
	else if(MainOrCrossTerms == 'c') pS0 = CrossTermsContribArray + StokesNo;

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

void srTRadIntWiggler::UpdateStokesSliceConstE(float* StokesCmpnArr, long ie, int StokesNo, char MainOrCrossTerms, srTStokesStructAccessData& StokesAccessData)
{
	float *pS0;
	if(MainOrCrossTerms == 'm') pS0 = StokesAccessData.pBaseSto + StokesNo;
	else if(MainOrCrossTerms == 'c') pS0 = CrossTermsContribArray + StokesNo;

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

int srTRadIntWiggler::TreatFiniteElecBeamEmittanceOneComp1D(float* CmpnArr, double ElBeamMomFact)
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

void srTRadIntWiggler::DetermineSingleElecDistrEffSizes1D(float* CmpnArr, char VsXorZ, double& M_Cen)
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

void srTRadIntWiggler::DetermineResizeBeforeConv1D(double M_ElecEff, double M_DistrSingleE, char VsXorZ, srTRadResize1D& Resize)
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

void srTRadIntWiggler::ConstructDataForConv1D(float* CmpnArr, float* AuxConvData, long NpOld, long NpNew)
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

int srTRadIntWiggler::PerformConvolutionWithGaussian1D(float* ConvData, long NewNp, double M_ElecEff, char VsXorZ)
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

//void srTRadIntWiggler::ExtractDataAfterConv1D(float* AuxConvData, long NpAux, long Np, float* CmpnArr)
void srTRadIntWiggler::ExtractDataAfterConv1D(float* AuxConvData, long long NpAux, long long Np, float* CmpnArr)
{
	//long iDat = (NpAux - Np) >> 1;
	long long iDat = (NpAux - Np) >> 1;

	float *tCmpn = CmpnArr;
	float *tAux = AuxConvData + (iDat << 1);
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		*(tCmpn++) = *tAux; tAux += 2;
	}
}

//*************************************************************************

int srTRadIntWiggler::TreatFiniteElecBeamEmittanceOneComp2D(float* CmpnArr, double ElBeamMomFact)
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

void srTRadIntWiggler::DetermineSingleElecDistrEffSizes2D(float* CmpnArr, double& MxxCen, double& MzzCen)
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

void srTRadIntWiggler::DetermineResizeBeforeConv2D(double MxxElecEff, double MzzElecEff, double MxxPowSingleE, double MzzPowSingleE, srTRadResize& Resize)
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

void srTRadIntWiggler::ConstructDataForConv2D(float* CmpnArr, float* NewData, long NewNx, long NewNz)
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
	float *tOldU = CmpnArr + (DistrInfoDat.nz - 1)*DistrInfoDat.nx;
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

int srTRadIntWiggler::PerformConvolutionWithGaussian2D(float* ConvData, long NewNx, long NewNz, double MxxElecEff, double MzzElecEff)
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
}

//*************************************************************************

void srTRadIntWiggler::ExtractDataAfterConv2D(float* AuxConvData, long NxAux, long NzAux, float* CmpnArr)
{
	long ixDat = (NxAux - DistrInfoDat.nx) >> 1;
	long izDat = (NzAux - DistrInfoDat.nz) >> 1;
	long Two_ixDat = ixDat << 1;
	long AuxPerZ = NxAux << 1;

	float *tCmpn = CmpnArr;
	for(long iz=0; iz<DistrInfoDat.nz; iz++)
	{
		float *tAux = AuxConvData + (izDat + iz)*AuxPerZ + Two_ixDat;
		for(long ix=0; ix<DistrInfoDat.nx; ix++)
		{
			*(tCmpn++) = *tAux; tAux += 2;
		}
	}
}

//*************************************************************************

void srTRadIntWiggler::SuppressNegativeValues(float* StokesCmpnArr)
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
