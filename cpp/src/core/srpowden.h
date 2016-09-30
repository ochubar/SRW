/************************************************************************//**
 * File: srpowden.h
 * Description: Calculation of Power Density (integrated over all photon energies) of Synchrotron Radiation from ~Arbitrary Transversely-Uniform Magnetic Field (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRPOWDEN_H
#define __SRPOWDEN_H

#include "srtrjdat.h"
#include "gmtrans.h"
#include "srpersto.h"
#include "srmlttsk.h"
#include "srradinc.h"

//*************************************************************************

extern srTYield srYield;
struct srTParPrecPowDens;

//*************************************************************************

struct srTRadIntPowDenPrec {
	double PrecFact;
	char Method; // 1- near-field, 2- far-field
	char UseSpecIntLim; // 0- don't use; 1- use
	double sStart, sFin;
};

//*************************************************************************

struct srTPairOfFloat {
	float f1, f2;
};

//*************************************************************************

class srTRadIntPowerDensity {

	srTCosAndSinComp CosAndSinComp;

	char LongIntTypeG; // 1- all field range; 2- one period;
	double sIntegStartG, sIntegFinG;
	int AmOfPerG;

	double *BtxArrP[50], *XArrP[50], *BxArrP[50];
	double *BtzArrP[50], *ZArrP[50], *BzArrP[50];
	//long AmOfPointsOnLevel[50];
	long long AmOfPointsOnLevel[50];
	long NumberOfLevelsFilledG;
	long MaxLevelForMeth_01G;
	char ProbablyTheSameLoopG;
	double MaxFluxDensValG, CurrentAbsPrecG, sIntegRelPrecG;
	double ActNormConstG;

	float NotCompIntervBorders[1100];
	srTPairOfFloat NotCompInterv[600];

	srTEXZ EXZ;
	srTFieldBasedArrays FieldBasedArrays;

	char MagFieldIsConstG; // !0 if const field at input
	gmTrans TrLab2Loc;
	srTInitTrjCoordAng LocInitCoordAng;
	TVector3d PobsLocG;
	double BconG, RmaG;

public:

	srTGenTrjHndl TrjHndl;
	srTWfrSmp DistrInfoDat;
	srTRadIntPowDenPrec IntPowDenPrec;

	srTIntVect* pWarningsGen;

	srTRadIntPowerDensity()
	{
		Initialize();
	}
	~srTRadIntPowerDensity()
	{
		DisposeAuxTrjArrays();
	}

	void Initialize();
	int CheckInputConsistency();

	void ComputePowerDensity(srTEbmDat* pElecBeam, srTMagElem* pMagElem, srTWfrSmp* pWfrSmp, srTParPrecPowDens* pPrecPowDens, srTPowDensStructAccessData* pPow);
	void ComputePowerDensity(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecPowDens* pPrecPowDens, srTPowDensStructAccessData* pPow); //SRWLib

	int ComputeTotalPowerDensityDistr(srTPowDensStructAccessData&);
	int ComputePowerDensityAtPoint(float* pPowDens);
	int SetUpFieldBasedArrays();
	void AnalizeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ);
	void FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTPowDensStructAccessData& PowDensAccessData);
	//void SetupNotCompIntervBorders(double MinValNotComp, double sStart, double sStep, long Np, long& AmOfInterv);
	void SetupNotCompIntervBorders(double MinValNotComp, double sStart, double sStep, long long Np, long long& AmOfInterv);
	int TreatFiniteElecBeamEmittance(srTPowDensStructAccessData&, gmTrans* pTrfObsPl =0);
	int TreatFiniteElecBeamEmittance1D(srTPowDensStructAccessData&, char);
	void DetermineSingleElecPowDensEffSizes(srTPowDensStructAccessData&, double& MxxPowSingleE, double& MzzPowSingleE);
	void DetermineResizeBeforeConv(double MxxElecEff, double MzzElecEff, double MxxPowSingleE, double MzzPowSingleE, srTRadResize& Resize);
	void ConstructDataForConv(srTPowDensStructAccessData& PowDensAccessData, float* NewData, long NewNx, long NewNz);
	int PerformConvolutionWithGaussian(float* AuxConvData, long NxAux, long NzAux, double MxxElecEff, double MzzElecEff);
	void ExtractFinalDataAfterConv(float* AuxConvData, long NxAux, long NzAux, srTPowDensStructAccessData& PowDensAccessData);

	void SetPrecParams(srTParPrecPowDens* pPrecPowDens);
	int TryToReduceIntegLimits();
	void SetupNativeRotation();
	int ComputePowerDensityAtPointConstMagField(float* pPowDens);

	//int FillNextLevel(int LevelNo, double sStart, double sEnd, long Np)
	int FillNextLevel(int LevelNo, double sStart, double sEnd, long long Np)
	{
		double* BasePtr = new double[Np*6];
		if(BasePtr == 0) return MEMORY_ALLOCATION_FAILURE;

		BtxArrP[LevelNo] = BasePtr;
		BasePtr += Np; XArrP[LevelNo] = BasePtr;
		BasePtr += Np; BxArrP[LevelNo] = BasePtr;
		BasePtr += Np; BtzArrP[LevelNo] = BasePtr;
		BasePtr += Np; ZArrP[LevelNo] = BasePtr;
		BasePtr += Np; BzArrP[LevelNo] = BasePtr;

		TrjHndl.rep->CompTotalTrjData(sStart, sEnd, Np, BtxArrP[LevelNo], BtzArrP[LevelNo], XArrP[LevelNo], ZArrP[LevelNo], BxArrP[LevelNo], BzArrP[LevelNo]);
		AmOfPointsOnLevel[LevelNo] = Np;
		NumberOfLevelsFilledG++;
		return 0;
	}
	void DisposeAuxTrjArrays()
	{
		for(int i=0; i<NumberOfLevelsFilledG; i++)
		{
			if(BtxArrP[i] != 0) delete[] (BtxArrP[i]);
		}
		NumberOfLevelsFilledG = 0;
		ZeroAuxTrjPtrs();
	}
	void ZeroAuxTrjPtrs()
	{
		for(int i=0; i<50; i++)
		{
			BtxArrP[i] = 0; XArrP[i] = 0; BxArrP[i] = 0;
			BtzArrP[i] = 0; ZArrP[i] = 0; BzArrP[i] = 0;
			AmOfPointsOnLevel[i] = 0;
		}
	}
	void PowDensFun(double s, double Bx, double Btx, double X, double Bz, double Btz, double Z, double& Fx, double& Fz)
	{
		if(s >= DistrInfoDat.yStart) 
		{
			Fx = 0; Fz = 0; return;
		}

		double One_d_ymis, Nx, Nz;
		double dx = (EXZ.x - X), dy = DistrInfoDat.yStart - s, dz = (EXZ.z - Z);
		if(dy == 0.) dy = 1.e-23;

		if(IntPowDenPrec.Method == 1) // Near field
		{
			//double instDist = DistrInfoDat.yStart - s;
			double instDist = sqrt(dx*dx + dy*dy + dz*dz);
			One_d_ymis = 1./instDist;

			//One_d_ymis = 1./dy;
			//Nx = (EXZ.x - X)*One_d_ymis; Nz = (EXZ.z - Z)*One_d_ymis;
			Nx = dx*One_d_ymis; Nz = dz*One_d_ymis;
		}
		else
		{
			double instDist = DistrInfoDat.yStart - 0.5*(sIntegStartG + sIntegFinG);
			One_d_ymis = 1./instDist;
			Nx = EXZ.x*One_d_ymis; Nz = EXZ.z*One_d_ymis;
		}
		double invRinst = One_d_ymis;

		double auxFact = Nx*Nx + Nz*Nz;
		//OCTEST (commented-out)
		//if(auxFact > 0.99) 
		//{
		//	double multFact = sqrt(0.99/auxFact);
		//	Nx *= multFact; Nz *= multFact;
		//	auxFact = 0.99; //to calculate angles more accurately??
		//}

		//OC27092016
		double Ny = 0.;
		if(auxFact < 1.)
		{
			Ny = (auxFact > 0.000001)? sqrt(1. - auxFact) : 1 - 0.5*auxFact - 0.125*auxFact*auxFact;
		}

		//double Ny = (auxFact > 0.01)? sqrt(1. - auxFact) : 1 - 0.5*auxFact - 0.125*auxFact*auxFact;
		//double Ny = (auxFact > 0.001)? sqrt(1. - auxFact) : 1 - 0.5*auxFact - 0.125*auxFact*auxFact;
		//double Ny = (auxFact > 0.000001)? sqrt(1. - auxFact) : 1 - 0.5*auxFact - 0.125*auxFact*auxFact; //OC27092016 (commented-out)

		double cosFactProj = Ny; //sqrt(1. - auxFact);
		double vEyP_x = 0, vEyP_y = 1, vEyP_z = 0;
		if(!DistrInfoDat.obsPlaneIsTransv)
		{
			//TVector3d &vNorm = DistrInfoDat.vNormObsPl;
			//cosFactProj = -((vNorm.x)*Nx + (vNorm.y)*sqrt(1. - auxFact) + (vNorm.z)*Nz);
			TVector3d &vEyP = DistrInfoDat.vLong;
			vEyP_x = vEyP.x; vEyP_y = vEyP.y; vEyP_z = vEyP.z;
			//cosFactProj = vEyP_x*Nx + vEyP_y*Ny + vEyP_z*Nz;
			cosFactProj = ::fabs(vEyP_x*Nx + vEyP_y*Ny + vEyP_z*Nz);
			//OC030412: to ensure positive power density even if "external" normal points towards e-beam 
			//(e.g. because of error in sign of coordinates of vectors defining the obs. plane orientation)
		}

		double Btx_mi_Nx = Btx - Nx, Btz_mi_Nz = Btz - Nz;
		double Gamma = TrjHndl.rep->EbmDat.Gamma;
		double GamE2 = Gamma*Gamma;
		double ConInvR = 586.674067035/Gamma;
		double InvRx = ConInvR*Bz, InvRz = ConInvR*Bx;
		//if(auxFact < 0.001) //to make const?
		//{//OC TEST: to comment-out
			double GamBtx_mi_Nx = Gamma*Btx_mi_Nx, GamBtz_mi_Nz = Gamma*Btz_mi_Nz;
			double GamBtx_mi_NxE2 = GamBtx_mi_Nx*GamBtx_mi_Nx, GamBtz_mi_NzE2 = GamBtz_mi_Nz*GamBtz_mi_Nz;
			double InvOne_p_GamBtE2 = 1./(1. + GamBtx_mi_NxE2 + GamBtz_mi_NzE2);
			double ComFact = One_d_ymis*One_d_ymis*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2;
			double TwoGamBt_mi_N = 2.*GamBtx_mi_Nx*GamBtz_mi_Nz;
			double BufFx = -TwoGamBt_mi_N*InvRz - (1. - GamBtx_mi_NxE2 + GamBtz_mi_NzE2)*InvRx + 2.*Btx_mi_Nx*One_d_ymis;
			double BufFz = (1. + GamBtx_mi_NxE2 - GamBtz_mi_NzE2)*InvRz + TwoGamBt_mi_N*InvRx + 2.*Btz_mi_Nz*One_d_ymis;
			Fx = ComFact*BufFx*BufFx*cosFactProj; Fz = ComFact*BufFz*BufFz*cosFactProj;

			double ComFact2 = ComFact*One_d_ymis/(GamE2*InvOne_p_GamBtE2);
			Fx += ComFact2*vEyP_x*BufFx; //to check
			Fz += ComFact2*vEyP_z*BufFz; //to check
		//	//to suppress separate treatment of Fx, Fz?
		//}
		//else
		//{
/**
			double twoGamE2 = 2*GamE2;
			double auxBetY = (1./GamE2) + Btx*Btx + Btz*Btz;
			TVector3d vBet(Btx, 1. - 0.5*auxBetY - 0.125*auxBetY*auxBetY, Btz), vN(Nx, Ny, Nz);
			TVector3d vBetPr(InvRx, -Btx*InvRx + Btz*InvRz, -InvRz), vN_mi_Bet = vN - vBet;
			TVector3d vA = (twoGamE2*(vN^(vN_mi_Bet^vBetPr))) + ((2.*invRinst)*vN_mi_Bet);

			//TVector3d vA = (twoGamE2*(vN^(vN_mi_Bet^vBetPr))); //OC TEST*********************
			
			double InvOne_p_GamBtE2 = 1./(twoGamE2*(1. - (vN*vBet)));
			double ComFact = invRinst*invRinst*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2*InvOne_p_GamBtE2;
			double ComFact2 = ComFact*invRinst/(GamE2*InvOne_p_GamBtE2);
			double ComFact_cosFactProj = ComFact*cosFactProj;
			
			double half_vAe2 = 0.5*vA.y*vA.y;
			//OC TEST*********************
			//To uncomment here!
			Fx = ComFact_cosFactProj*(vA.x*vA.x + half_vAe2); //artificial trick: adding longitudinal component to transverse ones
			Fz = ComFact_cosFactProj*(vA.z*vA.z + half_vAe2);
			//OC TEST*********************
			//Fx = ::fabs(ComFact_cosFactProj*(vA.x*vA.x)); //artificial trick: adding longitudinal component to transverse ones
			//Fz = ::fabs(ComFact_cosFactProj*(vA.z*vA.z));
			//END TEST

			double half_scalProdY = 0.5*vEyP_y*vA.y;
			//OC TEST*********************
			//To uncomment here!
			Fx = ::fabs(Fx - ComFact2*(vEyP_x*vA.x + half_scalProdY));
			Fz = ::fabs(Fz - ComFact2*(vEyP_z*vA.z + half_scalProdY));
**/
		//}
	}

	void DetermineElecEffSizes(double& MxxElecEff, double& MzzElecEff)
	{
		srTEbmDat& EbmDat = TrjHndl.rep->EbmDat;
		double R = DistrInfoDat.yStart - EbmDat.s0;
		double Re2 = R*R;
		MxxElecEff = EbmDat.Mxx + Re2*EbmDat.Mxpxp + 2.*R*EbmDat.Mxxp;
		MzzElecEff = EbmDat.Mzz + Re2*EbmDat.Mzpzp + 2.*R*EbmDat.Mzzp;
	}

	void PointLab2Loc(TVector3d& P)
	{
		TVector3d P_Loc = TrLab2Loc.TrPoint(P);
		P = P_Loc;
	}
	void VectorLoc2Lab(TVector3d& V)
	{
		TVector3d V_Lab = TrLab2Loc.TrBiPoint_inv(V);
		V = V_Lab;
	}
};

//*************************************************************************

#endif
