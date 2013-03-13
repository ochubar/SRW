/************************************************************************//**
 * File: srradinc.h
 * Description: SR calculation from Constant Magnetic Field (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRRADINC_H
#define __SRRADINC_H

#include "srtrjaux.h"
#include "gmtrans.h"
#include "srctrjdt.h"
#include "srmlttsk.h"

//*************************************************************************

extern srTYield srYield;
extern srTIntVect gVectWarnNos;

//*************************************************************************

struct srTRadIntConstPrec {
	double PrecFact;
};

//*************************************************************************

struct srTInitTrjCoordAng {
	double x0, dxds0, z0, dzds0;
};

//*************************************************************************

class srTRadIntConst {

	srTCosAndSinComp CosAndSinComp;

	double sIntegRelPrecG;
	double NormalizingConstG;

	srTEXZ EXZ;
	gmTrans TrLab2Loc;
	srTInitTrjCoordAng LocInitCoordAng;
	TVector3d PobsLocG;
	double BconG, RmaG;
	double Gx, Gz;

public:

	srTConstTrjDat* ConstTrjDatPtr;
	srTWfrSmp DistrInfoDat;
	srTRadIntConstPrec IntConstPrec;

	srTRadIntConst()
	{
		Initialize();
	}
	void Initialize() 
	{
		DistrInfoDat.CoordUnits = 0; // To ensure m for coord.
		CosAndSinComp.Initialize();
	}

	int CheckInputConsistency();

	int ComputeTotalStokesDistr(srTStokesStructAccessData&);
	int ComputeStokesAtPoint(float* pStokes);
	int ComputeElFieldAtPointLocFrame(srTEFourier& E);

	void SetupNativeRotation();

	int TreatFiniteElecBeamEmittance(srTStokesStructAccessData&, double ElBeamMomFact);
	void ExtractStokesSliceConstE(srTStokesStructAccessData& StokesAccessData, long ie, int StokesNo, float* pOutS);
	void UpdateStokesSliceConstE(float* StokesCmpnArr, long ie, int is, srTStokesStructAccessData& StokesAccessData);
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
	void ExtractDataAfterConv1D(float* AuxConvData, long NpAux, long Np, float* CmpnArr);
	void SuppressNegativeValues(float* StokesCmpnArr);

	void AnalyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ);
	void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData& StokesAccessData);

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
	int SetupThickBeamConsts()
	{
		srTEbmDat& Ebm = ConstTrjDatPtr->EbmDat;
		double y = DistrInfoDat.yStart - Ebm.s0;
		double ye2 = y*y;

		const double MinM2 = 1.E-20; // To steer
		double BufX = Ebm.Mxx + (Ebm.Mxpxp)*ye2 + 2.*(Ebm.Mxxp)*y;
		double BufZ = Ebm.Mzz + (Ebm.Mzpzp)*ye2 + 2.*(Ebm.Mzzp)*y;
		if((BufX <= 0.) || (BufZ <= 0.)) return THICK_EL_BEAM_WAS_NOT_SET_UP;

		if(BufX < MinM2) BufX = MinM2;
		if(BufZ < MinM2) BufZ = MinM2;
		Gx = 0.5/BufX;
		Gz = 0.5/BufZ;

		return 0;
	}
	void SetIntegPrecLevel()
	{
		const double NominalPrec = 0.01; // To steer
		sIntegRelPrecG = NominalPrec/IntConstPrec.PrecFact;
	}
	void SetupNormalizingConst()
	{//Assume Spectral flux density is in Photons/(s*mm^2*(0.1%BW))
		const double e_coulomb = 1.602189246E-19;
		const double Alpha = 1./137.0360411; // Fine-structure constant
		NormalizingConstG = sqrt(Alpha*(ConstTrjDatPtr->EbmDat.Current)*1.E+09/e_coulomb);
		DistrInfoDat.NormalizingConst = NormalizingConstG;
		RmaG = 3.33564076253*(ConstTrjDatPtr->EbmDat.Energy)/BconG; // Magnet radius in m assuming Energy[GeV] and BconG[T]
	}
	inline void PointLab2Loc(TVector3d& P)
	{
		P = TrLab2Loc.TrPoint(P);
	}
	void VectorLoc2Lab(TVector3d& V)
	{
		TVector3d V_Lab = TrLab2Loc.TrBiPoint_inv(V);
		V = V_Lab;
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
};

//*************************************************************************

#endif
