/************************************************************************//**
 * File: srpowden.cpp
 * Description: Calculation of Power Density (integrated over all photon energies) of Synchrotron Radiation from ~Arbitrary Transversely-Uniform Magnetic Field
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srpowden.h"
#include "srprgind.h"
#include "srctrjdt.h"
#include "srinterf.h"
#include "gmmeth.h"

//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

void srTRadIntPowerDensity::Initialize()
{
	DistrInfoDat.CoordUnits = 0; // To ensure m for coord.

	NumberOfLevelsFilledG = 0;
	MaxLevelForMeth_01G = 12; // To steer

	ProbablyTheSameLoopG = 0;
	MaxFluxDensValG = CurrentAbsPrecG = 0.;

	pWarningsGen = &gVectWarnNos;
}

//*************************************************************************

int srTRadIntPowerDensity::CheckInputConsistency()
{// Fill this

	if(IntPowDenPrec.UseSpecIntLim)
	{
		if(IntPowDenPrec.sStart >= IntPowDenPrec.sFin) return INCORRECT_INT_LIMITS_SR_POW_COMP;
	}

	return 0;
}

//*************************************************************************

int srTRadIntPowerDensity::SetUpFieldBasedArrays()
{
	int result;
	srTFieldBasedArrayKeys Keys;
	Keys.Bx_ = Keys.Bz_ = Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = 1;

	int MagType = TrjHndl.rep->MagFieldPeriodicity();
	if((IntPowDenPrec.Method == 2) || (MagType == 2)) // far-field & periodic
	{
		LongIntTypeG = 2; // integration over one period
		if(result = TrjHndl.rep->SetUpFieldBasedArraysAtOnePeriod(Keys, FieldBasedArrays)) return result;
	}
	else // integration over total range
	{
		LongIntTypeG = 1;
		if(result = TrjHndl.rep->SetUpFieldBasedArraysTotal(Keys, FieldBasedArrays)) return result;
	}
	return 0;
}

//*************************************************************************

void srTRadIntPowerDensity::SetPrecParams(srTParPrecPowDens* pPrecPowDens)
{
	if(pPrecPowDens == 0) return;
	IntPowDenPrec.Method = pPrecPowDens->MethNo;
    IntPowDenPrec.PrecFact = pPrecPowDens->PrecFact;
    IntPowDenPrec.UseSpecIntLim = pPrecPowDens->UseSpecIntLim;
    IntPowDenPrec.sStart = pPrecPowDens->sIntStart;
    IntPowDenPrec.sFin = pPrecPowDens->sIntFin;
}

//*************************************************************************

void srTRadIntPowerDensity::ComputePowerDensity(srTEbmDat* pElecBeam, srTMagElem* pMagElem, srTWfrSmp* pWfrSmp, srTParPrecPowDens* pPrecPowDens, srTPowDensStructAccessData* pPow)
{
	if((pElecBeam == 0) || (pMagElem == 0) || (pPow == 0)) throw INCORRECT_PARAMS_SR_COMP;

	int res = 0;
	Initialize();
	DistrInfoDat = *pWfrSmp;
	DistrInfoDat.EnsureZeroTransverseRangesForSinglePoints();
	pPow->SetRadSamplingFromObs(DistrInfoDat);

	//srTTrjDat* pTrjDat = new srTTrjDat(pElecBeam, pMagElem);
	//if(res = pTrjDat->ComputeInterpolatingStructure()) throw res;

	srTGenTrjDat* pGenTrjDat = srTGenTrjDat::CreateAndSetupNewTrjDat(pElecBeam, pMagElem);
	srTGenTrjHndl hGenTrj(pGenTrjDat);
	TrjHndl = hGenTrj;

	SetPrecParams(pPrecPowDens);
	if(res = ComputeTotalPowerDensityDistr(*pPow)) throw res;
}

//*************************************************************************

void srTRadIntPowerDensity::ComputePowerDensity(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecPowDens* pPrecPowDens, srTPowDensStructAccessData* pPow)
{//SRWLib
	if((pTrjDat == 0) || (pWfrSmp == 0) || (pPow == 0)) throw INCORRECT_PARAMS_SR_COMP;

	int res = 0;
	Initialize();
	DistrInfoDat = *pWfrSmp;
	DistrInfoDat.EnsureZeroTransverseRangesForSinglePoints();
	pPow->SetRadSamplingFromObs(DistrInfoDat);

	//srTGenTrjDat* pGenTrjDat = srTGenTrjDat::CreateAndSetupNewTrjDat(pElecBeam, pMagElem);
	//srTGenTrjHndl hGenTrj(pGenTrjDat);
	srTGenTrjHndl hGenTrj(pTrjDat, true);
	TrjHndl = hGenTrj;

	SetPrecParams(pPrecPowDens);
	if(res = ComputeTotalPowerDensityDistr(*pPow)) throw res;
}

//*************************************************************************

int srTRadIntPowerDensity::ComputeTotalPowerDensityDistr(srTPowDensStructAccessData& PowDensAccessData)
{
	int result;
	if(result = CheckInputConsistency()) return result;

	MagFieldIsConstG = TrjHndl.rep->MagFieldIsConstant();

	srTGenTrjHndl hGenTrjLoc = TrjHndl;
	if(MagFieldIsConstG && (IntPowDenPrec.Method == 1)) // Const. Magnetic Field and Near Field comp. method
	{
		char Periodicity = 1; // non-periodic
		if(result = hGenTrjLoc.rep->ConvertToArbTrjDat(Periodicity, DistrInfoDat, TrjHndl)) return result;
		MagFieldIsConstG = 0;
	}
	if(MagFieldIsConstG) SetupNativeRotation();

	int MagType = TrjHndl.rep->MagFieldPeriodicity();
	if((IntPowDenPrec.Method == 2) && (MagType == 2)) // far-field & periodic
	{
		LongIntTypeG = 2; // integration over one period
	}
	else
	{
		LongIntTypeG = 1; // integration over total range
	}

	DistrInfoDat.AssumeAllPhotonEnergies = 1; // for correct determination of the integration limits
	//if(result = TrjHndl.rep->ShowLimitsAndInitInteg(DistrInfoDat, LongIntTypeG, sIntegStartG, sIntegFinG, AmOfPerG)) return result;
	if(result = TrjHndl.rep->ShowLimitsAndInitInteg(DistrInfoDat, LongIntTypeG, sIntegStartG, sIntegFinG, AmOfPerG, false)) return result; //OC030412 don't calculate interpolating structure since it was already calculated outside

	if(IntPowDenPrec.UseSpecIntLim)
	{
		if(sIntegStartG < IntPowDenPrec.sStart) sIntegStartG = IntPowDenPrec.sStart;
		if(sIntegFinG > IntPowDenPrec.sFin) sIntegFinG = IntPowDenPrec.sFin;
	}

	if(((sIntegStartG > DistrInfoDat.yStart) || (sIntegFinG > DistrInfoDat.yStart)) && (IntPowDenPrec.Method == 1))
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, OBSERV_POINT_TOO_CLOSE_TO_INTEG_LIMITS);
	}

	if(DistrInfoDat.obsPlaneIsTransv)
	{
		if(result = TryToReduceIntegLimits()) return result;
	}

	sIntegRelPrecG = 0.06/IntPowDenPrec.PrecFact; // To steer
	DisposeAuxTrjArrays();
	MaxFluxDensValG = 0.;

	if(MagFieldIsConstG) 
	{
		double ElEn = TrjHndl.rep->EbmDat.Energy;
		double Robs = DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0;
		ActNormConstG = 18.11*ElEn*ElEn*ElEn*ElEn*ElEn*(TrjHndl.rep->EbmDat.Current)/(RmaG*Robs*Robs);
	}
	else
	{
		double Gam = TrjHndl.rep->EbmDat.Gamma;
		ActNormConstG = (9.167170453E-16)*Gam*Gam*Gam*Gam*Gam*Gam*(TrjHndl.rep->EbmDat.Current); // To make result in W/mm^2
	}

	gmTrans trfObsPl;
	bool trfObsPlaneIsDefined = DistrInfoDat.SetupTrfObsPlaneIfNecessary(trfObsPl);

	double *pObSurfData = PowDensAccessData.m_spObSurfData.rep;
	bool obSurfIsDefined = (pObSurfData != 0);

	TVector3d &vExP = DistrInfoDat.vHor, &vEyP = DistrInfoDat.vLong, vEzP, vEyP0, vExP0;

	char FinalResAreSymOverX = 0, FinalResAreSymOverZ = 0;
	if((!trfObsPlaneIsDefined) && (!obSurfIsDefined)) AnalizeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ); //to make more general

	double xc = TrjHndl.rep->EbmDat.dxds0*(DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0) + TrjHndl.rep->EbmDat.x0;
	double zc = TrjHndl.rep->EbmDat.dzds0*(DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0) + TrjHndl.rep->EbmDat.z0;

	double xStep = (DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.;
	double zStep = (DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.;
	double xTol = xStep*0.001; // To steer
	double zTol = zStep*0.001; // To steer

	double xStart = DistrInfoDat.xStart;
	double zStart = DistrInfoDat.zStart;
	if(trfObsPlaneIsDefined)
	{//definition in local frame
		xStart = -0.5*(DistrInfoDat.xEnd - DistrInfoDat.xStart);
		zStart = -0.5*(DistrInfoDat.zEnd - DistrInfoDat.zStart);
	}

	//TVector3d vCenSurf(0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd), DistrInfoDat.yStart, 0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd));
	//if(obSurfIsDefined)
	//{
	//	vCenSurf.y = CGenMathMeth::interpFunc2D(vCenSurf.x, vCenSurf.z, DistrInfoDat.xStart, DistrInfoDat.xEnd, DistrInfoDat.nx, DistrInfoDat.zStart, DistrInfoDat.zEnd, DistrInfoDat.nz, pObSurfData);
	//}

	//long PerZ = DistrInfoDat.nx;
	long long PerZ = DistrInfoDat.nx;

	//long TotalAmOfOutPoints = DistrInfoDat.nz*DistrInfoDat.nx;
	long long TotalAmOfOutPoints = DistrInfoDat.nz*DistrInfoDat.nx;

	if(FinalResAreSymOverX) TotalAmOfOutPoints >>= 1;
	if(FinalResAreSymOverZ) TotalAmOfOutPoints >>= 1;

	//long PointCount = 0;
	long long PointCount = 0;
	double UpdateTimeInt_s = 0.5;
	srTCompProgressIndicator CompProgressInd(TotalAmOfOutPoints, UpdateTimeInt_s);

	TVector3d vRlab(0, 0, 0), vRloc(0, 0, 0);
	double yStartOrig = DistrInfoDat.yStart;

	for(int iz=0; iz<DistrInfoDat.nz; iz++)
	{
		//EXZ.z = DistrInfoDat.zStart + iz*zStep;
		EXZ.z = zStart + iz*zStep; //OC140110
		vRloc.z = EXZ.z;
		if(FinalResAreSymOverZ) { if((EXZ.z - zc) > zTol) break;}

		for(int ix=0; ix<DistrInfoDat.nx; ix++)
		{
			//EXZ.x = DistrInfoDat.xStart + ix*xStep;
			EXZ.x = xStart + ix*xStep; //OC140110

			if(trfObsPlaneIsDefined)
			{
				if(obSurfIsDefined)
				{
					//Deviation of Long. coord of obs. point without a space transform.:
					vRloc.y = CGenMathMeth::tabTangOrtsToSurf2D(vExP0, vEzP, ix, iz, DistrInfoDat.nx, DistrInfoDat.nz, xStep, zStep, pObSurfData);
					vEyP0 = vEzP^vExP0; //vLong before space transform.

					vEyP = trfObsPl.TrBiPoint(vEyP0); //vLong after space transform., to be used in PowDensFun
					vExP = trfObsPl.TrBiPoint(vExP0); //vHor after space transform.
				}

				vRloc.x = EXZ.x;
				//if there is no surface, vRloc.y should be 0 - is it correct?
				vRlab = trfObsPl.TrPoint(vRloc);
				EXZ.x = vRlab.x;
				EXZ.z = vRlab.z;
				DistrInfoDat.yStart = vRlab.y;
			}
			else
			{
				if(obSurfIsDefined)
				{
					//DistrInfoDat.yStart = CGenMathMeth::tabFunc2D(ix, iz, DistrInfoDat.nx, pObSurfData);
					DistrInfoDat.yStart = yStartOrig + CGenMathMeth::tabTangOrtsToSurf2D(vExP, vEzP, ix, iz, DistrInfoDat.nx, DistrInfoDat.nz, xStep, zStep, pObSurfData);
					vEyP = vEzP^vExP; //defines vLong, to be used in PowDensFun

					//TVector3d mRow1(vExP.x, vEyP.x, vEzP.x);
					//TVector3d mRow2(vExP.y, vEyP.y, vEzP.y);
					//TVector3d mRow3(vExP.z, vEyP.z, vEzP.z);
					//TMatrix3d M(mRow1, mRow2, mRow3);
					//trfObsPl.SetMatrixVector(M, vCenSurf);
				}
			}

			if(FinalResAreSymOverX) { if((EXZ.x - xc) > xTol) break;}

			if(MagFieldIsConstG)
			{
				PobsLocG.x = EXZ.x; PobsLocG.y = 0.; PobsLocG.z = EXZ.z;
				PobsLocG = TrLab2Loc.TrPoint(PobsLocG);
			}

			float* pPowDens = PowDensAccessData.pBasePowDens + (iz*PerZ + ix);
			if(result = ComputePowerDensityAtPoint(pPowDens)) return result;

			DistrInfoDat.yStart = yStartOrig;
			if(result = srYield.Check()) return result;
			if(result = CompProgressInd.UpdateIndicator(PointCount++)) return result;
		}
	}

	if(FinalResAreSymOverZ || FinalResAreSymOverX) 
		FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, PowDensAccessData);

//To make this optional?
	if(result = TreatFiniteElecBeamEmittance(PowDensAccessData, &trfObsPl)) return result; //to take into accout eventual tilt of the observation plane!
	return 0;
}

//*************************************************************************

int srTRadIntPowerDensity::TryToReduceIntegLimits()
{// call only if defined MagFieldIsConstG, LongIntTypeG, sIntegStartG, sIntegFinG
	if(LongIntTypeG == 2) return 0;
	if(MagFieldIsConstG) return 0;

	const int AmOfPoToAnalize = 500; // To steer
	double SafetyAng = 5./(TrjHndl.rep->EbmDat.Gamma); // To steer

	double xMin = DistrInfoDat.xStart, xMax = DistrInfoDat.xStart;
	if(DistrInfoDat.nx > 1) xMax = DistrInfoDat.xEnd;
	double zMin = DistrInfoDat.zStart, zMax = DistrInfoDat.zStart;
	if(DistrInfoDat.nz > 1) zMax = DistrInfoDat.zEnd;

	double sStart = sIntegStartG;
	double sEnd = sIntegFinG;
	double sStep = (sEnd - sStart)/(AmOfPoToAnalize - 1);

	double* DataCont = new double[AmOfPoToAnalize*6];
	if(DataCont == 0) return MEMORY_ALLOCATION_FAILURE;
	double* BasePtr = DataCont;

	double* LocBtxArrP = BasePtr; BasePtr += AmOfPoToAnalize; 
	double* LocXArrP = BasePtr; BasePtr += AmOfPoToAnalize; 
	double* LocBxArrP = BasePtr; BasePtr += AmOfPoToAnalize; 
	double* LocBtzArrP = BasePtr; BasePtr += AmOfPoToAnalize; 
	double* LocZArrP = BasePtr; BasePtr += AmOfPoToAnalize; 
	double* LocBzArrP = BasePtr;
	TrjHndl.rep->CompTotalTrjData(sStart, sEnd, AmOfPoToAnalize, LocBtxArrP, LocBtzArrP, LocXArrP, LocZArrP, LocBxArrP, LocBzArrP);

	int AmOfPo_mi_1 = AmOfPoToAnalize - 1;
	int isStartCor = 0, isEndCor = AmOfPo_mi_1;
	double s = sStart;
	double *tBtx = LocBtxArrP, *tBtz = LocBtzArrP, *tX = LocXArrP, *tZ = LocZArrP;
	for(int i=0; i<AmOfPoToAnalize; i++)
	{
		double y1Inv = 1./(DistrInfoDat.yStart - s);
		double NxMin = (xMin - *tX)*y1Inv, NxMax = (xMax - *tX)*y1Inv;
		double NzMin = (zMin - *tZ)*y1Inv, NzMax = (zMax - *tZ)*y1Inv;

		char xIsGood = ((*tBtx < (NxMin - SafetyAng)) || (*tBtx > (NxMax + SafetyAng)));
		char zIsGood = ((*tBtz < (NzMin - SafetyAng)) || (*tBtz > (NzMax + SafetyAng)));
		if(xIsGood && zIsGood) isStartCor = i;
		else break;

		tBtx++; tBtz++; tX++; tZ++;
		s += sStep;
	}

	s = sEnd;
	tBtx = LocBtxArrP + AmOfPo_mi_1; tBtz = LocBtzArrP + AmOfPo_mi_1; 
	tX = LocXArrP + AmOfPo_mi_1; tZ = LocZArrP + AmOfPo_mi_1;
	for(int ir=AmOfPo_mi_1; ir>=0; ir--)
	{
		double y1Inv = 1./(DistrInfoDat.yStart - s);
		double NxMin = (xMin - *tX)*y1Inv, NxMax = (xMax - *tX)*y1Inv;
		double NzMin = (zMin - *tZ)*y1Inv, NzMax = (zMax - *tZ)*y1Inv;

		char xIsGood = ((*tBtx < (NxMin - SafetyAng)) || (*tBtx > (NxMax + SafetyAng)));
		char zIsGood = ((*tBtz < (NzMin - SafetyAng)) || (*tBtz > (NzMax + SafetyAng)));
		if(xIsGood && zIsGood) isEndCor = ir;
		else break;

		tBtx--; tBtz--; tX--; tZ--;
		s -= sStep;
	}

	if(isStartCor < isEndCor)
	{//OC100214
		if(isStartCor > 0) sIntegStartG += isStartCor*sStep;
		if(isEndCor < AmOfPo_mi_1) sIntegFinG -= (AmOfPo_mi_1 - isEndCor)*sStep;
	}

	if(DataCont != 0) delete[] DataCont;
	return 0;
}

//*************************************************************************

void srTRadIntPowerDensity::SetupNativeRotation()
{
	srTConstTrjDat* ConstTrjDatPtr = (srTConstTrjDat*)(TrjHndl.rep);

	TVector3d uLab(ConstTrjDatPtr->MagConst.Bx, 0., ConstTrjDatPtr->MagConst.Bz);
	double Bcon = sqrt(uLab.x*uLab.x + uLab.z*uLab.z);
	uLab = (1./Bcon)*uLab;
	TVector3d uLoc(0.,0.,1.), Zero(0.,0.,0.);

	BconG = Bcon;
	RmaG = 3.33564076253*(ConstTrjDatPtr->EbmDat.Energy)/BconG; // Magnet radius in m assuming Energy[GeV] and BconG[T]

	TVector3d uLab_mi_uLoc = uLab - uLoc;
	//double NormEst = Abs(uLab_mi_uLoc);
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

void srTRadIntPowerDensity::AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;

	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;
	TrjHndl.rep->AnalizeFieldSymmetry(FieldIsSymOverX, FieldIsSymOverZ);
	if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;

	char ObsIsSymOverX = 0, ObsIsSymOverZ = 0;
	if(FieldIsSymOverX && (DistrInfoDat.nx > 1))
	{
		double xStep = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1);
		double xTol = xStep*0.01; // To steer
		double xc = TrjHndl.rep->EbmDat.dxds0*(DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0) + TrjHndl.rep->EbmDat.x0;
		ObsIsSymOverX = (::fabs(0.5*(DistrInfoDat.xStart + DistrInfoDat.xEnd) - xc) < xTol);
	}
	if(FieldIsSymOverZ && (DistrInfoDat.nz > 1))
	{
		double zStep = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1);
		double zTol = zStep*0.01; // To steer
		double zc = TrjHndl.rep->EbmDat.dzds0*(DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0) + TrjHndl.rep->EbmDat.z0;
		ObsIsSymOverZ = (::fabs(0.5*(DistrInfoDat.zStart + DistrInfoDat.zEnd) - zc) < zTol);
	}

	FinalResAreSymOverX = (FieldIsSymOverX && ObsIsSymOverX);
	FinalResAreSymOverZ = (FieldIsSymOverZ && ObsIsSymOverZ);
}

//*************************************************************************

int srTRadIntPowerDensity::ComputePowerDensityAtPoint(float* pPowDens)
{
	if(MagFieldIsConstG) return ComputePowerDensityAtPointConstMagField(pPowDens);

	const double wfe = 7./15.;
	const double wf1 = 16./15.;
	const double wf2 = 14./15.;
	const double wd = 1./15.;

	//long AmOfExtrInBx = TrjHndl.rep->AmOfExtremInBx, AmOfExtrInBz = TrjHndl.rep->AmOfExtremInBz;
	long long AmOfExtrInBx = TrjHndl.rep->AmOfExtremInBx, AmOfExtrInBz = TrjHndl.rep->AmOfExtremInBz;
	//long MinNpAcceptedForExtrem = (AmOfExtrInBx > AmOfExtrInBz)? AmOfExtrInBx : AmOfExtrInBz;
	long long MinNpAcceptedForExtrem = (AmOfExtrInBx > AmOfExtrInBz)? AmOfExtrInBx : AmOfExtrInBz;
	if(MinNpAcceptedForExtrem < 7) MinNpAcceptedForExtrem = 7;
	//char MultiExtremCase = (MinNpAcceptedForExtrem > 100);
	
	//long MinNpAcceptedForField = TrjHndl.rep->EstimMinNpForRadInteg(2); //OC07052010
	long long MinNpAcceptedForField = TrjHndl.rep->EstimMinNpForRadInteg(2); //OC07052010

	//const long NpOnLevelMaxNoResult = 500000000; //5000000; // To steer; to stop computation as unsuccessful
	const long long NpOnLevelMaxNoResult = 500000000; //5000000; // To steer; to stop computation as unsuccessful
	const long NpOnLevelCritResultMayBeSmall = 256; // To steer
	const int MaxAmOfNotCompInterv = 50; // To steer

	const long LevelNoCritToAnalizeNotComp = 9; // To steer
	double RelRatioNotComp = (1.E-08)*sIntegRelPrecG; //1.E-05*sIntegRelPrecG; // To steer
	//long AmOfNotCompInterv = 0;
	long long AmOfNotCompInterv = 0;

	double GamEm1 = 1./TrjHndl.rep->EbmDat.Gamma;
	//double MaxBetDifAllowed = 2.*GamEm1; // To steer
	double MaxBetDifAllowed = 2.*GamEm1/IntPowDenPrec.PrecFact; //OC180408 ??
	
	//long NpOnLevel = 5; // Must be non-even!
	long long NpOnLevel = 5; // Must be non-even!

	double sStart = sIntegStartG;
	double sEnd = sIntegFinG;
	double sStep = (sEnd - sStart)/(NpOnLevel - 1);

	int result;
	if(NumberOfLevelsFilledG == 0) if(result = FillNextLevel(0, sStart, sEnd, NpOnLevel)) return result;

	double Sum1X=0., Sum1Z=0., Sum2X=0., Sum2Z=0.;
	double wFx, wFz, Fx, Fz;
	int LevelNo = 0;

	double *pBtx = *BtxArrP, *pBtz = *BtzArrP, *pX = *XArrP, *pZ = *ZArrP, *pBx = *BxArrP, *pBz = *BzArrP;
	double BtxLoc, xLoc, BxLoc, BtzLoc, zLoc, BzLoc;

	PowDensFun(sStart, *(pBx++), *(pBtx++), *(pX++), *(pBz++), *(pBtz++), *(pZ++), wFx, wFz);

	double s = sStart + sStep;
	//int AmOfPass = (NpOnLevel - 3) >> 1;
	long long AmOfPass = (NpOnLevel - 3) >> 1;
	//for(int i=0; i<AmOfPass; i++)
	for(long long i=0; i<AmOfPass; i++)
	{
		PowDensFun(s, *(pBx++), *(pBtx++), *(pX++), *(pBz++), *(pBtz++), *(pZ++), Fx, Fz);
		Sum1X += Fx; Sum1Z += Fz; s += sStep;
		PowDensFun(s, *(pBx++), *(pBtx++), *(pX++), *(pBz++), *(pBtz++), *(pZ++), Fx, Fz);
		Sum2X += Fx; Sum2Z += Fz; s += sStep;
	}
	PowDensFun(s, *(pBx++), *(pBtx++), *(pX++), *(pBz++), *(pBtz++), *(pZ++), Fx, Fz);
	Sum1X += Fx; Sum1Z += Fz; s += sStep;

	PowDensFun(s, *pBx, *pBtx, *pX, *pBz, *pBtz, *pZ, Fx, Fz);
	wFx += Fx; wFz += Fz;
	wFx *= wfe; wFz *= wfe;

	double DifDerX = 0; // *InitDerMan - *FinDerMan;
	double DifDerZ = 0; // *InitDerMan - *FinDerMan;
	double wDifDerX = wd*DifDerX, wDifDerZ = wd*DifDerZ;

	double ActNormConst = ActNormConstG*AmOfPerG; // To make proper units: W/mm^2
	double ActNormConst_sStep = ActNormConstG*sStep;

	double IntX = ActNormConst_sStep*(wFx + wf1*Sum1X + wf2*Sum2X + sStep*wDifDerX);
	double IntZ = ActNormConst_sStep*(wFz + wf1*Sum1Z + wf2*Sum2Z + sStep*wDifDerZ);
	double SqNorm = IntX + IntZ;

	NpOnLevel--;
	char NotFinishedYet = 1;
	char ExtraPassForAnyCase = 0; 
	while(NotFinishedYet)
	{
		Sum2X += Sum1X; Sum2Z += Sum1Z;
		Sum1X = Sum1Z = 0.;

		char ThisMayBeTheLastLoop = 0;
		LevelNo++;

		double HalfStep = 0.5*sStep;
		s = sStart + HalfStep;

		if(LevelNo <= MaxLevelForMeth_01G)
		{
			if(NumberOfLevelsFilledG <= LevelNo) if(result = FillNextLevel(LevelNo, s, sEnd - HalfStep, NpOnLevel)) return result;
			pBtx = BtxArrP[LevelNo]; pBtz = BtzArrP[LevelNo]; pX = XArrP[LevelNo]; pZ = ZArrP[LevelNo]; pBx = BxArrP[LevelNo]; pBz = BzArrP[LevelNo];
		}

		double MaxF = 0.;
		double BtxPrev, BtzPrev, MaxBtxDif = 0., MaxBtzDif = 0.;
		//long NpOnLevel_mi_1 = NpOnLevel - 1;
		long long NpOnLevel_mi_1 = NpOnLevel - 1;

		char PrevPointWasComputed = 0;
		int NotCompIntervCount = 0;
		for(int i=0; i<NpOnLevel; i++)
		{
			char ComputeThisPoint = 1;
			if((LevelNo > LevelNoCritToAnalizeNotComp) && (AmOfNotCompInterv > 0) && (AmOfNotCompInterv < MaxAmOfNotCompInterv))
			{
				srTPairOfFloat *pPair = NotCompInterv + NotCompIntervCount;
				if(s > pPair->f1)
				{
					if(s <= pPair->f2)
					{
						ComputeThisPoint = 0; PrevPointWasComputed = 0;
					}
					else
					{
						if(NotCompIntervCount < (AmOfNotCompInterv - 1)) NotCompIntervCount++;
						if((s > pPair->f1) && (s <= pPair->f2)) 
						{ 
							ComputeThisPoint = 0; PrevPointWasComputed = 0;
						}
					}
				}
			}

			if(ComputeThisPoint)
			{
				if(LevelNo > MaxLevelForMeth_01G)
				{
					pBtx = &BtxLoc; pX = &xLoc; pBx = &BxLoc;
					pBtz = &BtzLoc; pZ = &zLoc; pBz = &BzLoc;
					TrjHndl.rep->CompTrjDataDerivedAtPointPowDens(s, *pBtx, *pBtz, *pX, *pZ, *pBx, *pBz);
				}

					//DEBUG
					if(i == 287)
					{
						int aha = 1;
					}

				PowDensFun(s, *pBx, *pBtx, *pX, *pBz, *pBtz, *pZ, Fx, Fz);
				Sum1X += Fx; Sum1Z += Fz;

				if(LevelNo == LevelNoCritToAnalizeNotComp)
				{
					double Ftot = Fx + Fz;
					NotCompIntervBorders[i] = (float)Ftot;
					if(Ftot > MaxF) MaxF = Ftot;
				}

				if(PrevPointWasComputed && (i < NpOnLevel_mi_1))
				{
					double BtxDif = ::fabs(*pBtx - BtxPrev), BtzDif = ::fabs(*pBtz - BtzPrev);
					if(BtxDif > MaxBtxDif) MaxBtxDif = BtxDif;
					if(BtzDif > MaxBtzDif) MaxBtzDif = BtzDif;
				}
				BtxPrev = *pBtx; BtzPrev = *pBtz;
				PrevPointWasComputed = 1;
			}

			if(LevelNo <= MaxLevelForMeth_01G)
			{
				pBtx++; pX++; pBx++; pBtz++; pZ++; pBz++;
			}
			s += sStep;
		}

		if((MaxBtxDif < MaxBetDifAllowed) && (MaxBtzDif < MaxBetDifAllowed)) ThisMayBeTheLastLoop = 1;

		double ActNormConstHalfStep = ActNormConst*HalfStep;
		double LocIntX = ActNormConstHalfStep*(wFx + wf1*Sum1X + wf2*Sum2X + HalfStep*wDifDerX);
		if(LocIntX < 0) LocIntX = 0.; //OC030412
		double LocIntZ = ActNormConstHalfStep*(wFz + wf1*Sum1Z + wf2*Sum2Z + HalfStep*wDifDerZ);
		if(LocIntZ < 0) LocIntZ = 0.; //OC030412
		double LocSqNorm = LocIntX + LocIntZ;

		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			char SharplyGoesDown = (LocSqNorm < 0.2*SqNorm);

			char NotFinishedYetFirstTest;
			if(ProbablyTheSameLoopG && (MaxFluxDensValG > 0.)) NotFinishedYetFirstTest = (TestVal > CurrentAbsPrecG);
			else NotFinishedYetFirstTest = (TestVal > sIntegRelPrecG*LocSqNorm);

			//if(NpOnLevel < MinNpAcceptedForExtrem) NotFinishedYetFirstTest = 1;
			if((NpOnLevel < MinNpAcceptedForExtrem) || (NpOnLevel < MinNpAcceptedForField)) NotFinishedYetFirstTest = 1;

			if(!NotFinishedYetFirstTest)
			{
				//if(ExtraPassForAnyCase || SharplyGoesDown || MultiExtremCase) NotFinishedYet = 0;
				if(ExtraPassForAnyCase || SharplyGoesDown) NotFinishedYet = 0;
				else ExtraPassForAnyCase = 1;
			}

			if((::fabs(LocSqNorm) < CurrentAbsPrecG) && (::fabs(SqNorm) < CurrentAbsPrecG))
			{
				if(NpOnLevel < NpOnLevelCritResultMayBeSmall) 
				{
					NotFinishedYet = 1; ExtraPassForAnyCase = 0;
				}
			}
		}

		if(NotFinishedYet)
		{
			if(LevelNo == LevelNoCritToAnalizeNotComp) SetupNotCompIntervBorders(MaxF*RelRatioNotComp, sStart, sStep, NpOnLevel, AmOfNotCompInterv);

			if(NpOnLevel > NpOnLevelMaxNoResult) 
			{
				return CAN_NOT_COMPUTE_RADIATION_INTEGRAL;
			}
		}

		IntX = LocIntX; IntZ = LocIntZ;
		SqNorm = LocSqNorm;
		sStep = HalfStep; NpOnLevel *= 2;
	}

	*pPowDens = (float)(IntX + IntZ);

	if((*pPowDens < 0)) *pPowDens = 0.;//OC

	if((ProbablyTheSameLoopG && (MaxFluxDensValG < SqNorm)) || !ProbablyTheSameLoopG) 
	{
		MaxFluxDensValG = SqNorm; CurrentAbsPrecG = sIntegRelPrecG*MaxFluxDensValG; 
		ProbablyTheSameLoopG = 1;
	}
	return 0;
}

//*************************************************************************

int srTRadIntPowerDensity::ComputePowerDensityAtPointConstMagField(float* pPowDens)
{
	//int result;
	double y1Inv = 1./(DistrInfoDat.yStart - TrjHndl.rep->EbmDat.s0);
	//double Xob_mi_x0 = PobsLocG.x - LocInitCoordAng.x0;
	double Zob_mi_z0 = PobsLocG.z - LocInitCoordAng.z0;
	double Zang = Zob_mi_z0*y1Inv - LocInitCoordAng.dzds0;
	double GamZang = (TrjHndl.rep->EbmDat.Gamma)*Zang;
	double GamZangE2 = GamZang*GamZang;
	double IntX = ActNormConstG/pow(1 + GamZangE2, 2.5);
	double IntZ = (0.714285714286)*IntX*GamZangE2/(1 + GamZangE2);
	*pPowDens = (float)(IntX + IntZ);
	return 0;
}

//*************************************************************************

//void srTRadIntPowerDensity::SetupNotCompIntervBorders(double MinValNotComp, double sStart, double sStep, long Np, long& AmOfInterv)
void srTRadIntPowerDensity::SetupNotCompIntervBorders(double MinValNotComp, double sStart, double sStep, long long Np, long long& AmOfInterv)
{
	double HalfStep = 0.5*sStep;
	double s = sStart + HalfStep;
	//long IntervCount = 0;
	long long IntervCount = 0;
	char IntervStarted = 0;

	float *tNotCompIntervBorders = NotCompIntervBorders;
	srTPairOfFloat *tNotCompInterv = NotCompInterv;
	double PrevVal = 0.;
	char DerivSign = 1;
	//for(long ip=0; ip<Np; ip++)
	for(long long ip=0; ip<Np; ip++)
	{
		if(*tNotCompIntervBorders < MinValNotComp)
		{
			if(!IntervStarted)
			{
				tNotCompInterv->f1 = (float)s; 
				IntervStarted = 1;
			}
			else
			{
				if((*tNotCompIntervBorders < PrevVal) && (DerivSign > 0))
				{
					double sPrev = s - sStep;
					if(sPrev > (tNotCompInterv->f1 + 0.1*sStep))
					{
						tNotCompInterv->f2 = (float)sPrev;
						tNotCompInterv++;
						IntervCount++;
					}
					IntervStarted = 0;
				}
			}
		}
		else
		{
			if(IntervStarted)
			{
				double sPrev = s - sStep;
				if(sPrev > (tNotCompInterv->f1 + 0.1*sStep))
				{
					tNotCompInterv->f2 = (float)sPrev;
					tNotCompInterv++;
					IntervCount++;
				}
				IntervStarted = 0;
			}
		}

		DerivSign = (*tNotCompIntervBorders > PrevVal)? 1 : -1;
		PrevVal = *tNotCompIntervBorders;

		tNotCompIntervBorders++;
		s += sStep;
	}
	if(IntervStarted)
	{
		double sEnd = s - HalfStep;
		tNotCompInterv->f2 = (float)sEnd;
		IntervStarted = 0;
		IntervCount++;
	}

	if(::fabs(NotCompInterv[0].f1 - sStart - HalfStep) < 0.5*HalfStep)
	{
		NotCompInterv[0].f1 = (float)sStart;
	}

	AmOfInterv = IntervCount;
}

//*************************************************************************

void srTRadIntPowerDensity::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTPowDensStructAccessData& PowDensAccessData)
{
	//long PerX = 1;
	//long PerZ = PerX*DistrInfoDat.nx;
	long long PerX = 1;
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
					float* pOrigData = PowDensAccessData.pBasePowDens + izPerZ + ix*PerX;
					float* pSymData = PowDensAccessData.pBasePowDens + izPerZ + (Nx_mi_1 - ix)*PerX;
					*pSymData = *pOrigData;
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
				float* pOrigData = PowDensAccessData.pBasePowDens + (izPerZ + ixPerX);
				float* pSymData = PowDensAccessData.pBasePowDens + (BufZ + ixPerX);
				*pSymData = *pOrigData;
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
				float* pOrigData = PowDensAccessData.pBasePowDens + izPerZ + ix*PerX;
				float* pSymData = PowDensAccessData.pBasePowDens + izPerZ + (Nx_mi_1 - ix)*PerX;
				*pSymData = *pOrigData;
			}
		}
	}
}

//*************************************************************************

int srTRadIntPowerDensity::TreatFiniteElecBeamEmittance(srTPowDensStructAccessData& PowDensAccessData, gmTrans* pTrfObsPl)
{
	int result;
	char Treat1D = 0;
	if(PowDensAccessData.nx == 1) 
	{
		if(PowDensAccessData.nz == 1)
		{
			//srTSend Send; Send.AddWarningMessage(pWarningsGen, BEAM_EMITTANCE_WAS_NOT_TAKEN_INTO_ACCOUNT);
			CErrWarn::AddWarningMessage(pWarningsGen, BEAM_EMITTANCE_WAS_NOT_TAKEN_INTO_ACCOUNT);
			return 0;
		}
		else Treat1D = 'x';
	}
	if(PowDensAccessData.nz == 1) Treat1D = 'z';
	if(Treat1D != 0) return TreatFiniteElecBeamEmittance1D(PowDensAccessData, Treat1D);

	double MxxElecEff, MzzElecEff;
	DetermineElecEffSizes(MxxElecEff, MzzElecEff);

	if((MxxElecEff <= 0) || (MzzElecEff <= 0)) //OC210708
	{
		CErrWarn::AddWarningMessage(pWarningsGen, BEAM_EMITTANCE_WAS_NOT_TAKEN_INTO_ACCOUNT);
		return 0;
	}

	if((!DistrInfoDat.obsPlaneIsTransv) && (pTrfObsPl != 0))
	{//to program!!
		//TVector3d vProjElecDimsX(sqrt(MxxElecEff), 0, 0), vProjElecDimsZ(0, 0, sqrt(MzzElecEff));

		//TVector3d vLocProjElecDims = pTrfObsPl->TrBiPoint_inv(vProjElecDims);
		//MxxElecEff = vLocProjElecDims.x*vLocProjElecDims.x;
		//MzzElecEff = vLocProjElecDims.z*vLocProjElecDims.z;
	}

	double MxxPowSingleE, MzzPowSingleE;
	DetermineSingleElecPowDensEffSizes(PowDensAccessData, MxxPowSingleE, MzzPowSingleE);
	srTRadResize Resize;
	DetermineResizeBeforeConv(MxxElecEff, MzzElecEff, MxxPowSingleE, MzzPowSingleE, Resize);

	long NxAux = (long)(Resize.pxm*PowDensAccessData.nx);
	long NzAux = (long)(Resize.pzm*PowDensAccessData.nz);
	CGenMathFFT2D FFT;
	FFT.NextCorrectNumberForFFT(NxAux);
	FFT.NextCorrectNumberForFFT(NzAux);

	float* AuxConvData = new float[NxAux*NzAux << 1];
	if(AuxConvData == 0) return MEMORY_ALLOCATION_FAILURE;

	ConstructDataForConv(PowDensAccessData, AuxConvData, NxAux, NzAux);
	if(result = PerformConvolutionWithGaussian(AuxConvData, NxAux, NzAux, MxxElecEff, MzzElecEff)) return result;
	ExtractFinalDataAfterConv(AuxConvData, NxAux, NzAux, PowDensAccessData);

	PowDensAccessData.EnsureNonNegativeValues(); //OC

	if(AuxConvData != 0) delete[] AuxConvData;
	return 0;
}

//*************************************************************************

int srTRadIntPowerDensity::TreatFiniteElecBeamEmittance1D(srTPowDensStructAccessData& PowDensAccessData, char Treat1D)
{
	//srTSend Send; Send.AddWarningMessage(pWarningsGen, BEAM_EMITTANCE_WAS_NOT_TAKEN_INTO_ACCOUNT);
	CErrWarn::AddWarningMessage(pWarningsGen, BEAM_EMITTANCE_WAS_NOT_TAKEN_INTO_ACCOUNT);
// Later make special 1D emittance treatment and send another warning

	return 0;
}

//*************************************************************************

void srTRadIntPowerDensity::DetermineSingleElecPowDensEffSizes(srTPowDensStructAccessData& PowDensAccessData, double& MxxPowSingleE, double& MzzPowSingleE)
{// Compute central moments
	float Sm0 = 0., SmX = 0., SmZ = 0., SmXX = 0., SmZZ = 0.;
	float xStep = (float)((DistrInfoDat.nx > 1)? (DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1) : 0.);
	float zStep = (float)((DistrInfoDat.nz > 1)? (DistrInfoDat.zEnd - DistrInfoDat.zStart)/(DistrInfoDat.nz - 1) : 0.);

	long nz_mi_1 = DistrInfoDat.nz - 1, nx_mi_1 = DistrInfoDat.nx - 1;
	float *tPowDens = PowDensAccessData.pBasePowDens;
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
			//float PowDens = wz*(*tPowDens);
			float PowDens = wz*(*(tPowDens++)); // ??
			if((ix == nx_mi_1) || (iz == nz_mi_1)) PowDens *= 0.5;

			Sm0 += PowDens;
			SmX += x*PowDens;
			SmXX += xe2*PowDens;
			SmZ += z*PowDens;
			SmZZ += ze2*PowDens;

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

	MxxPowSingleE = Mxx - Mx*Mx;
	MzzPowSingleE = Mzz - Mz*Mz;
}

//*************************************************************************

void srTRadIntPowerDensity::DetermineResizeBeforeConv(double MxxElecEff, double MzzElecEff, double MxxPowSingleE, double MzzPowSingleE, srTRadResize& Resize)
{
	const double DoNotResizeRatio = 5.; //3.; // To steer
	const double AmOfExtraSig = 6.; //2.5;

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

void srTRadIntPowerDensity::ConstructDataForConv(srTPowDensStructAccessData& PowDensAccessData, float* NewData, long NewNx, long NewNz)
{
	long ixDat = (NewNx - PowDensAccessData.nx) >> 1;
	long izDat = (NewNz - PowDensAccessData.nz) >> 1;

	//long OffsetXp = ixDat + PowDensAccessData.nx;
	//long OffsetZp = izDat + PowDensAccessData.nz;
	long long OffsetXp = ixDat + PowDensAccessData.nx;
	long long OffsetZp = izDat + PowDensAccessData.nz;

	//long NewPerZ = NewNx << 1;
	long long NewPerZ = NewNx << 1;
	//long ix, iz;
	long long ix, iz;
	float V = *(PowDensAccessData.pBasePowDens);
	for(iz=0; iz<izDat; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ);
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(PowDensAccessData.pBasePowDens + PowDensAccessData.nx - 1);
	for(iz=0; iz<izDat; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(PowDensAccessData.pBasePowDens + (PowDensAccessData.nz - 1)*PowDensAccessData.nx);
	for(iz=OffsetZp; iz<NewNz; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ);
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}
	V = *(PowDensAccessData.pBasePowDens + PowDensAccessData.nz*PowDensAccessData.nx - 1);
	for(iz=OffsetZp; iz<NewNz; iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
	}

	float *tOldL = PowDensAccessData.pBasePowDens;
	float *tOldR = PowDensAccessData.pBasePowDens + PowDensAccessData.nx - 1;
	for(iz=izDat; iz<(izDat + PowDensAccessData.nz); iz++)
	{
		V = *tOldL;
		//long izNewPerZ = iz*NewPerZ;
		long long izNewPerZ = iz*NewPerZ;
		float *tNew = NewData + izNewPerZ;
		for(ix=0; ix<ixDat; ix++) { *(tNew++) = V; *(tNew++) = 0;}
		tOldL += PowDensAccessData.nx;

		V = *tOldR;
		tNew = NewData + (izNewPerZ + (OffsetXp << 1));
		for(ix=OffsetXp; ix<NewNx; ix++) { *(tNew++) = V; *(tNew++) = 0;}
		tOldR += PowDensAccessData.nx;
	}

	float *tOldD = PowDensAccessData.pBasePowDens;
	float *tOldU = PowDensAccessData.pBasePowDens + (((long long)(PowDensAccessData.nz - 1))*((long long)PowDensAccessData.nx));
	for(ix=ixDat; ix<(ixDat + PowDensAccessData.nx); ix++)
	{
		V = *(tOldD++);
		float *tNew = NewData + (ix << 1);
		for(iz=0; iz<izDat; iz++) { *tNew = V; *(tNew+1) = 0; tNew += NewPerZ;}

		V = *(tOldU++);
		tNew = NewData + (OffsetZp*NewPerZ + (ix << 1));
		for(iz=OffsetZp; iz<NewNz; iz++) { *tNew = V; *(tNew+1) = 0; tNew += NewPerZ;}
	}

	float *tOld = PowDensAccessData.pBasePowDens;
	for(iz=izDat; iz<(izDat + PowDensAccessData.nz); iz++)
	{
		float *tNew = NewData + (iz*NewPerZ + (ixDat << 1));
		for(ix=ixDat; ix<(ixDat + PowDensAccessData.nx); ix++)
		{
			*(tNew++) = *(tOld++); *(tNew++) = 0.;
		}
	}
}

//*************************************************************************

int srTRadIntPowerDensity::PerformConvolutionWithGaussian(float* ConvData, long NewNx, long NewNz, double MxxElecEff, double MzzElecEff)
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

	if(result = FFT2D.Make2DFFT(FFT2DInfo)) return result;
	return 0;
}

//*************************************************************************

void srTRadIntPowerDensity::ExtractFinalDataAfterConv(float* AuxConvData, long NxAux, long NzAux, srTPowDensStructAccessData& PowDensAccessData)
{
	//long ixDat = (NxAux - PowDensAccessData.nx) >> 1;
	//long izDat = (NzAux - PowDensAccessData.nz) >> 1;
	long long ixDat = (NxAux - PowDensAccessData.nx) >> 1;
	long long izDat = (NzAux - PowDensAccessData.nz) >> 1;

	//long Two_ixDat = ixDat << 1;
	//long AuxPerZ = NxAux << 1;
	long long Two_ixDat = ixDat << 1;
	long long AuxPerZ = NxAux << 1;

	float *tPow = PowDensAccessData.pBasePowDens;
	//for(long iz=0; iz<PowDensAccessData.nz; iz++)
	for(long long iz=0; iz<PowDensAccessData.nz; iz++)
	{
		float *tAux = AuxConvData + ((izDat + iz)*AuxPerZ + Two_ixDat);
		//for(long ix=0; ix<PowDensAccessData.nx; ix++)
		for(long long ix=0; ix<PowDensAccessData.nx; ix++)
		{
			*(tPow++) = *tAux; tAux += 2;
		}
	}
}

//*************************************************************************
