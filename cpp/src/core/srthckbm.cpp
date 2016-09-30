/************************************************************************//**
 * File: srthckbm.cpp
 * Description: SR Stokes parameters calculation method for ~Arbitrary magnetic field (used rarely, under-programmed)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srthckbm.h"
#include "srmagfld.h"
#include "srmagcnt.h"
#include "srprgind.h"
//#include "srmamet.h"
#include "gmmeth.h"
#include "gmfunc.h"

//*************************************************************************

void srTRadIntThickBeamAuxParams::Setup(srTEbmDat& Elec)
{
	PI = 3.141592653590;
	//TwoPI = 2.*PI;
	//ThreePIdTwo = 1.5*PI;
	//FivePIdFour = 1.25*PI;
	//HalfPI = 0.5* PI;
	//One_dTwoPI = 1./TwoPI;
	double PIm10e6 = PI*1.E+06;
	Half_k_d_e = PIm10e6*0.80654658; // PIm10e6dEnCon
	k_d_e = 2.*Half_k_d_e;

	//const double e_esu = 4.80324214E-10; // Charge of electron in esu
	//const double c = 2.99792458E+10; // Speed of light in cm/s

	BetX = Elec.BetaPartDistr('x');
    BetZ = Elec.BetaPartDistr('z');
    AlpX = Elec.AlphaPartDistr('x');
    AlpZ = Elec.AlphaPartDistr('z');
	GamX = Elec.GammaPartDistr('x');
    GamZ = Elec.GammaPartDistr('z');
	ThetXZ = Elec.ThetaPartDistr(0, 0);
	ThetXpZ = Elec.ThetaPartDistr(1, 0);
	ThetXZp = Elec.ThetaPartDistr(0, 1);
	ThetXpZp = Elec.ThetaPartDistr(1, 1);
	Bgam = Elec.BgamPartDistr();
	GammaEm2 = Elec.GammaEm2;

	InvBgam = 1./Bgam;
	HalfInvBgam = 0.5*InvBgam;

	xc = Elec.x0;
	xpc = Elec.dxds0;
	zc = Elec.z0;
	zpc = Elec.dzds0;
	double ArgFc = -BetX*xpc*xpc - BetZ*zpc*zpc - 2*ThetXpZp*xpc*zpc - 2*AlpX*xpc*xc - 2*ThetXZp*zpc*xc 
        -GamX*xc*xc - 2*AlpZ*zpc*zc - 2*ThetXpZ*xpc*zc - 2*ThetXZ*xc*zc - GamZ*zc*zc;
	Fc = exp(ArgFc);

	double AuxFn = -2*BetX*BetZ*GamX*GamZ + (-AlpX*AlpX + BetX*GamX)*(-AlpZ*AlpZ + BetZ*GamZ) + 2*AlpZ*GamX*ThetXpZ*ThetXpZp 
		+2*AlpX*BetZ*ThetXpZ*ThetXZ + (BetX*BetZ - ThetXpZp*ThetXpZp)*(GamX*GamZ - ThetXZ*ThetXZ) + 2*AlpX*GamZ*ThetXpZp*ThetXZp 
		+2*AlpZ*BetX*ThetXZ*ThetXZp - 2*ThetXpZ*ThetXpZp*ThetXZ*ThetXZp - 2*AlpX*AlpZ*(ThetXpZp*ThetXZ + ThetXpZ*ThetXZp) 
		+(BetX*GamZ - ThetXpZ*ThetXpZ)*(BetZ*GamX - ThetXZp*ThetXZp);
	Fn = sqrt(Bgam*AuxFn);

	const double Alpha = 1./137.0360411; // Fine-structure constant
	const double e_coulomb = 1.602189246E-19; // Charge of electron in C
	const double ConvPhEn = (0.80654658E-03);
	const double AuxC1 = Alpha*ConvPhEn*ConvPhEn*(1.E+09)/e_coulomb;
	C1 = (Elec.Current)*AuxC1; // To be multiplied by squared photon energy in [eV]
	Cm = C1*Fc*Fn; 

/**
	const double e_esu = 4.80324214E-10; // Charge of electron in esu
	const double e_coulomb = 1.602189246E-19;
	const double c = 2.99792458E+10; // Speed of light in cm/s
	const double Alpha = 1./137.0360411; // Fine-structure constant
	const double PI = 3.141592653590;
	const double TwoPI = 2.*PI;

	if((DistrInfoDat.DistrValType == StokesParam) || (DistrInfoDat.DistrValType == FieldFourier))
	{
		if(DistrInfoDat.CoordOrAngPresentation == CoordPres) 
			NormalizingConst = sqrt(Alpha*(TrjDatPtr->EbmDat.Current)*(1.E+09)/e_coulomb);
		else if(DistrInfoDat.CoordOrAngPresentation == AngPres)
			NormalizingConst = sqrt(Alpha*(TrjDatPtr->EbmDat.Current)*(1.E+03)/e_coulomb);
	}

	DistrInfoDat.NormalizingConst = NormalizingConst;
**/

    TComplexD AuxV0(0., -2*(AlpZ*zc + BetZ*zpc + xpc*ThetXpZp + xc*ThetXZp));
	Aux_deltz1000 = AuxV0;

    TComplexD AuxV1(0., -2.*(AlpZ*zpc + zc*GamZ + xpc*ThetXpZ + xc*ThetXZ));
	Aux_alpz1000 = AuxV1;

    TComplexD AuxV2(0., -2.*(AlpX*xc + BetX*xpc + zc*ThetXpZ + zpc*ThetXpZp));
	Aux_deltx1000 = AuxV2;

	TComplexD AuxV3(0., -2.*(AlpX*xpc + zpc*ThetXZp + GamX*xc + ThetXZ*zc));
	Aux_alpx1000 = AuxV3;
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokes(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes)
{
	if((pMagFldTrUnif == 0) && (pMagLensCont == 0)) throw NO_MAG_FIELD_DEFINED;
    if((pElecBeam == 0) || (pPrcPar == 0) || (pStokes == 0)) throw INCORRECT_PARAMS_SR_COMP;

	srTRadIntThickBeam* pRadInt = 0;
	try
	{
        pRadInt = new srTRadIntThickBeam();
		if(pPrcPar->MethNo == 3) pRadInt->ComputeTotalStokesDistrViaSingleElec(pElecBeam, pMagFldTrUnif, pPrcPar, pStokes);
		else pRadInt->ComputeTotalStokesDistr(pElecBeam, pMagFldTrUnif, pMagLensCont, pPrcPar, pStokes);
	}
	catch(int ErrNo) 
	{ 
		if(pRadInt != 0) delete pRadInt;
		throw ErrNo;
	}
	if(pRadInt != 0) delete pRadInt;
}

//*************************************************************************

void srTRadIntThickBeam::ComputeTotalStokesDistr(srTEbmDat* pElecBeam, srTMagFldTrUnif* pMagFldTrUnif, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar, srTStokesStructAccessData* pStokes)
{
	if((pMagFldTrUnif == 0) && (pMagLensCont == 0)) throw NO_MAG_FIELD_DEFINED;
    if((pElecBeam == 0) || (pPrcPar == 0)) throw INCORRECT_PARAMS_SR_COMP;
	if(pStokes == 0) throw NO_STOKES_STRUCTURE_SUPPLIED;
	
    srTTrjDat *pTrjDat=0;
	if(pMagFldTrUnif != 0) pTrjDat = (srTTrjDat*)(pMagFldTrUnif->CreateAndSetupNewTrjDat(pElecBeam));
	gAuxPar.Setup(*pElecBeam);
	
	char FinalResAreSymOverX=0, FinalResAreSymOverZ=0;
	AnalyzeFinalResultsSymmetry(FinalResAreSymOverX, FinalResAreSymOverZ, pElecBeam, pTrjDat, pMagLensCont, pStokes);
	srTCompProgressIndicator CompProgressInd(FindTotalAmOfPointsToCalc(pStokes, FinalResAreSymOverX, FinalResAreSymOverZ), 0.5, 1);

	SetupInitialTrajArrays(pTrjDat, pMagLensCont, pPrcPar);

	//long PerE = 4;
	//long PerX = (pStokes->ne)*PerE;
	//long PerZ = (pStokes->nx)*PerX;
	//long PerY = (pStokes->nz)*PerZ;
	long long PerE = 4;
	long long PerX = (pStokes->ne)*PerE;
	long long PerZ = (pStokes->nx)*PerX;
	long long PerY = (pStokes->nz)*PerZ;

	double xc = pElecBeam->x0;
	double zc = pElecBeam->z0;
	double xTol = (pStokes->xStep)*0.001, zTol = (pStokes->zStep)*0.001; // To steer

	int res = 0;
    float* pBaseStokes = pStokes->pBaseSto;

	srTEXZY EXZY;
    EXZY.dx = pStokes->dx;
    EXZY.dz = pStokes->dz;
	//if(pPrcPar->IntOrFlux == 'f')
	//{
	//	EXZY.dx = ((pStokes->nx) - 1)*(pStokes->xStep);
	//	EXZY.dz = ((pStokes->nz) - 1)*(pStokes->zStep);
	//}

	//loop according to pStokes
	EXZY.y = pStokes->yStart;
	for(int iy=0; iy<pStokes->ny; iy++)
	{
		//long iyPerY = iy*PerY;
		long long iyPerY = iy*PerY;
		EXZY.e = pStokes->eStart;
		for(int ie=0; ie<pStokes->ne; ie++)
		{
            //long iePerE = ie*PerE;
            long long iePerE = ie*PerE;

            ComputeExpCoefXZArraysForInteg2D(EXZY.y, EXZY.e, *pPrcPar);

			EXZY.z = pStokes->zStart;
			for(int iz=0; iz<pStokes->nz; iz++)
			{
				if(FinalResAreSymOverZ) { if((EXZY.z - zc) > zTol) break;}
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;

				EXZY.x = pStokes->xStart;
				for(int ix=0; ix<pStokes->nx; ix++)
				{
					if(FinalResAreSymOverX) { if((EXZY.x - xc) > xTol) break;}
					//long ixPerX = ix*PerX;
					long long ixPerX = ix*PerX;

					srTStokes CurSt;
					ComputeStokesAtOneObsPoint(EXZY, *pPrcPar, CurSt);

					float* pSto = pBaseStokes + (iyPerY + izPerZ + ixPerX + iePerE);
					*(pSto++) = (float)CurSt.s0; *(pSto++) = (float)CurSt.s1; *(pSto++) = (float)CurSt.s2; *pSto = (float)CurSt.s3;

					if(res = CompProgressInd.UpdateIndicator()) throw res;
					EXZY.x += pStokes->xStep;
				}
				EXZY.z += pStokes->zStep;
			}
			EXZY.e += pStokes->eStep;
		}
		EXZY.y += pStokes->yStep;
	}

	if(FinalResAreSymOverZ || FinalResAreSymOverX) FillInSymPartsOfResults(FinalResAreSymOverX, FinalResAreSymOverZ, pStokes);
	if(pTrjDat != 0) delete pTrjDat;
}

//*************************************************************************

void srTRadIntThickBeam::SetupInitialTrajArrays(srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont, srTParPrecStokesArb* pPrcPar)
{
    if((pTrUnifTrjDat == 0) && (pMagLensCont == 0)) throw INCORRECT_PARAMS_SR_COMP;
    if(pPrcPar == 0) throw INCORRECT_PARAMS_SR_COMP;
	if(pPrcPar->RelPrecOrStep <= 0) throw INCORRECT_STEP_OR_RELATIVE_PRECISION;

	if(pPrcPar->MethNo == 1)
	{// constant grid, fixed step vs s
        gFldArr.sStep = pPrcPar->RelPrecOrStep;
		double sEnd;
		DetermineLongPosGridLimits(pTrUnifTrjDat, pMagLensCont, gFldArr.sStart, sEnd);
		gFldArr.Ns = int((sEnd - gFldArr.sStart)/(gFldArr.sStep));

		gFldArr.Ns = ((gFldArr.Ns >> 1) << 1) + 1; // to ensure odd
		if(gFldArr.Ns < 5) throw TOO_LARGE_LONGITUD_INTEG_STEP;
		gFldArr.sStep = (sEnd - gFldArr.sStart)/(gFldArr.Ns - 1);

		srTFieldBasedArrayKeys K;
		K.Btx_ = K.Btz_ = K.X_ = K.Z_ = K.IntBtxE2_ = K.IntBtzE2_ = 1;
        K.X1p_ = K.Z1p_ = K.X2p_ = K.Z2p_ = K.X1_ = K.Z1_ = K.X2_ = K.Z2_ = 1;
        K.IntX01_ = K.IntX02_ = K.IntX11_ = K.IntX12_ = K.IntX22_ = K.IntZ01_ = K.IntZ02_ = K.IntZ11_ = K.IntZ12_ = K.IntZ22_ = 1;
        gFldArr.AllocateArrays(gFldArr.Ns, K);
		ComputeTrajArrays(gFldArr, pTrUnifTrjDat, pMagLensCont);

		AllocateCoefArraysForInteg2D(gFldArr.Ns);
		AllocateFuncArraysForExternInteg(gFldArr.Ns);
	}
	else if(pPrcPar->MethNo == 2)
	{
        //program new method here
	}
}

//*************************************************************************

void srTRadIntThickBeam::DetermineLongPosGridLimits(srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont, double& sStart, double& sEnd)
{
    if((pTrUnifTrjDat == 0) && (pMagLensCont == 0)) throw INCORRECT_PARAMS_SR_COMP;
	if(pTrUnifTrjDat != 0)
	{
        sStart = pTrUnifTrjDat->sStart;
        sEnd = pTrUnifTrjDat->sStart + (pTrUnifTrjDat->LenFieldData - 1)*(pTrUnifTrjDat->sStep);
	}
	else
	{
        sStart = pMagLensCont->gsStart;
        sEnd = pMagLensCont->gsEnd;
	}
}

//*************************************************************************

void srTRadIntThickBeam::ComputeTrajArrays(srTFieldBasedArrays& FldArr, srTTrjDat* pTrUnifTrjDat, srTMagFldCont* pMagLensCont)
{
    srTFieldBasedArrayKeys Keys;
	if(pMagLensCont != 0)
	{
		ComputeOffAxisTrajArrays(FldArr, pMagLensCont);

		if(pTrUnifTrjDat != 0)
		{
            Keys.ZeroAllKeys();
            Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = Keys.IntBtxE2_ = Keys.IntBtzE2_ = 1;
            Keys.IntX01toS0_ = Keys.IntX02toS0_ = Keys.IntX11toS0_ = Keys.IntX12toS0_ = Keys.IntX22toS0_ = 1;
			Keys.IntZ01toS0_ = Keys.IntZ02toS0_ = Keys.IntZ11toS0_ = Keys.IntZ12toS0_ = Keys.IntZ22toS0_ = 1;
			Keys.IntX01_ = Keys.IntX02_ = Keys.IntX11_ = Keys.IntX12_ = Keys.IntX22_ = Keys.IntZ01_ = Keys.IntZ02_ = Keys.IntZ11_ = Keys.IntZ12_ = Keys.IntZ22_ = 1;
			pTrUnifTrjDat->CompTotalTrjData(Keys, FldArr);
		}
	}
	else
	{
		if(pTrUnifTrjDat != 0)
		{
            Keys.ZeroAllKeys();
            Keys.Btx_ = Keys.Btz_ = Keys.X_ = Keys.Z_ = Keys.IntBtxE2_ = Keys.IntBtzE2_ = 1;
			Keys.X1p_ = Keys.Z1p_ = Keys.X2p_ = Keys.Z2p_ = Keys.X1_ = Keys.Z1_ = Keys.X2_ = Keys.Z2_ = 1;
            Keys.IntX01toS0_ = Keys.IntX02toS0_ = Keys.IntX11toS0_ = Keys.IntX12toS0_ = Keys.IntX22toS0_ = 1;
			Keys.IntZ01toS0_ = Keys.IntZ02toS0_ = Keys.IntZ11toS0_ = Keys.IntZ12toS0_ = Keys.IntZ22toS0_ = 1;
			Keys.IntX01_ = Keys.IntX02_ = Keys.IntX11_ = Keys.IntX12_ = Keys.IntX22_ = Keys.IntZ01_ = Keys.IntZ02_ = Keys.IntZ11_ = Keys.IntZ12_ = Keys.IntZ22_ = 1;
			pTrUnifTrjDat->CompTotalTrjData(Keys, FldArr);
		}
	}
}

//*************************************************************************

void srTRadIntThickBeam::ComputeOffAxisTrajArrays(srTFieldBasedArrays& FldArr, srTMagFldCont* pMagLensCont)
{
	if(pMagLensCont == 0) return;

    double *pX1p = FldArr.X1pArr, *pZ1p = FldArr.Z1pArr;
    double *pX2p = FldArr.X2pArr, *pZ2p = FldArr.Z2pArr;
    double *pX1 = FldArr.X1Arr, *pZ1 = FldArr.Z1Arr;
	double *pX2 = FldArr.X2Arr, *pZ2 = FldArr.Z2Arr;

	double s = FldArr.sStart;
	double sStp = FldArr.sStep;
	for(int i=0; i<FldArr.Ns; i++)
	{
        TMatrix2d Mx(1,0,1,0), Mz(1,0,1,0);
		pMagLensCont->ComputeParticlePropagMatrix(s, Mx, Mz);

		*(pX1++) = Mx.Str0.x; *(pX2++) = Mx.Str0.y;
		*(pX1p++) = Mx.Str1.x; *(pX2p++) = Mx.Str1.y;
		*(pZ1++) = Mz.Str0.x; *(pZ2++) = Mz.Str0.y;
		*(pZ1p++) = Mz.Str1.x; *(pZ2p++) = Mz.Str1.y;
		s += sStp;
	}
}

//*************************************************************************

void srTRadIntThickBeam::AnalyzeFinalResultsSymmetry(char& FinalResAreSymOverX, char& FinalResAreSymOverZ, srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTMagFldCont* pMagLensCont, srTStokesStructAccessData* pStokes)
{
	FinalResAreSymOverX = FinalResAreSymOverZ = 0;
	char FieldIsSymOverX = 0, FieldIsSymOverZ = 0;

	if(pMagLensCont != 0) return;
	
	if(pTrjDat != 0)
	{
        if(!pTrjDat->HorFieldIsNotZero) FieldIsSymOverZ = 1;
        if(!pTrjDat->VerFieldIsNotZero) FieldIsSymOverX = 1;
        if((!FieldIsSymOverX) && (!FieldIsSymOverZ)) return;
        FinalResAreSymOverX = FieldIsSymOverX;
		FinalResAreSymOverZ = FieldIsSymOverZ;
	}
	if((pElecBeam != 0) && (pStokes != 0))
	{// Add more checks if more params added to ElecBeam
        FinalResAreSymOverX = FinalResAreSymOverZ = 0;
        if(FieldIsSymOverX && (pStokes->nx > 1))
        {
            double xTol = (pStokes->xStep)*0.01; // To steer
            char TrjAngIsSmall = (::fabs((pElecBeam->dxds0)*(pStokes->yStart - pElecBeam->s0)) < xTol);
            char ObsIsSymOverX = TrjAngIsSmall && (::fabs((pStokes->xStart + 0.5*(pStokes->nx - 1)*(pStokes->xStep)) - pElecBeam->x0) < xTol);
            FinalResAreSymOverX = (FinalResAreSymOverX && ObsIsSymOverX);
        }
        if(FieldIsSymOverZ && (pStokes->nz > 1))
        {
            double zTol = (pStokes->zStep)*0.01; // To steer
            char TrjAngIsSmall = (::fabs((pElecBeam->dzds0)*(pStokes->yStart - pElecBeam->s0)) < zTol);
            char ObsIsSymOverZ = TrjAngIsSmall && (::fabs((pStokes->zStart + 0.5*(pStokes->nz - 1)*(pStokes->zStep)) - pElecBeam->z0) < zTol);
            FinalResAreSymOverZ = (FinalResAreSymOverZ && ObsIsSymOverZ);
        }
	}
}

//*************************************************************************

//long srTRadIntThickBeam::FindTotalAmOfPointsToCalc(srTStokesStructAccessData* pStokes, char FinalResAreSymOverX, char FinalResAreSymOverZ)
long long srTRadIntThickBeam::FindTotalAmOfPointsToCalc(srTStokesStructAccessData* pStokes, char FinalResAreSymOverX, char FinalResAreSymOverZ)
{
	//long TotalAmOfOutPoints = (pStokes->ny)*(pStokes->nz)*(pStokes->nx)*(pStokes->ne);
	long long TotalAmOfOutPoints = ((long long)(pStokes->ny))*((long long)(pStokes->nz))*((long long)(pStokes->nx))*((long long)(pStokes->ne));
	if(FinalResAreSymOverX && (pStokes->nx > 1)) TotalAmOfOutPoints >>= 1;
	if(FinalResAreSymOverZ && (pStokes->nz > 1)) TotalAmOfOutPoints >>= 1;
	return TotalAmOfOutPoints;
}

//*************************************************************************

void srTRadIntThickBeam::FillInSymPartsOfResults(char FinalResAreSymOverX, char FinalResAreSymOverZ, srTStokesStructAccessData* pStokes)
{
    if(pStokes == 0) return;

	//long PerX = (pStokes->ne) << 2;
	//long PerZ = PerX*(pStokes->nx);
	//long PerY = PerZ*(pStokes->nz);
	long long PerX = (pStokes->ne) << 2;
	long long PerZ = PerX*(pStokes->nx);
	long long PerY = PerZ*(pStokes->nz);
	char SymWithRespectToXax, SymWithRespectToZax;

	int HalfNz = (pStokes->nz) >> 1, Nz_mi_1 = (pStokes->nz) - 1;
	if(FinalResAreSymOverZ && FinalResAreSymOverX)
	{
		if((HalfNz << 1) != (pStokes->nz)) HalfNz++;
	}
	//int izStart = ((HalfNz << 1) == (pStokes->nz))? HalfNz : (HalfNz + 1);

	int HalfNx = (pStokes->nx) >> 1, Nx_mi_1 = (pStokes->nx) - 1;
	//int ixStart = ((HalfNx << 1) == (pStokes->nx))? HalfNx : (HalfNx + 1);
	int iy, iz, ix;

	if(FinalResAreSymOverZ)
	{
		if(FinalResAreSymOverX)
		{
			SymWithRespectToXax = 0; SymWithRespectToZax = 1;
			for(iy=0; iy<pStokes->ny; iy++)
			{
				//long iyPerY = iy*PerY;
				long long iyPerY = iy*PerY;
                for(iz=0; iz<HalfNz; iz++)
                {
					//long izPerZ = iz*PerZ;
					//long iyPerY_p_izPerZ = iyPerY + izPerZ;
                    long long izPerZ = iz*PerZ;
					long long iyPerY_p_izPerZ = iyPerY + izPerZ;
                    for(ix=0; ix<HalfNx; ix++)
                    {
                        float* pOrigData = pStokes->pBaseSto + (iyPerY_p_izPerZ + ix*PerX);
                        float* pSymData = pStokes->pBaseSto + (iyPerY_p_izPerZ + (Nx_mi_1 - ix)*PerX);
                        CopySymEnergySlice(pOrigData, pSymData, pStokes->ne, SymWithRespectToXax, SymWithRespectToZax);
					}
				}
			}
		}
		SymWithRespectToXax = 1; SymWithRespectToZax = 0;

		for(iy=0; iy<pStokes->ny; iy++)
		{
			//long iyPerY = iy*PerY;
			long long iyPerY = iy*PerY;
			for(iz=0; iz<HalfNz; iz++)
			{
				//long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
				//long iyPerY_p_izPerZ = iyPerY + izPerZ;
				//long iyPerY_p_BufZ = iyPerY + BufZ;
				long long izPerZ = iz*PerZ, BufZ = (Nz_mi_1 - iz)*PerZ;
				long long iyPerY_p_izPerZ = iyPerY + izPerZ;
				long long iyPerY_p_BufZ = iyPerY + BufZ;

				for(ix=0; ix<pStokes->nx; ix++)
				{
                    //long ixPerX = ix*PerX;
                    long long ixPerX = ix*PerX;
                    float* pOrigData = pStokes->pBaseSto + (iyPerY_p_izPerZ + ixPerX);
                    float* pSymData = pStokes->pBaseSto + (iyPerY_p_BufZ + ixPerX);
                    CopySymEnergySlice(pOrigData, pSymData, pStokes->ne, SymWithRespectToXax, SymWithRespectToZax);
				}
			}
		}
	}
	else if(FinalResAreSymOverX)
	{
		SymWithRespectToXax = 0; SymWithRespectToZax = 1;

		for(iy=0; iy<pStokes->ny; iy++)
		{
			//long iyPerY = iy*PerY;
			long long iyPerY = iy*PerY;
            for(iz=0; iz<pStokes->nz; iz++)
            {
				//long izPerZ = iz*PerZ;
				//long iyPerY_p_izPerZ = iyPerY + izPerZ;
                long long izPerZ = iz*PerZ;
				long long iyPerY_p_izPerZ = iyPerY + izPerZ;

				for(ix=0; ix<HalfNx; ix++)
                {
                    float* pOrigData = pStokes->pBaseSto + (iyPerY_p_izPerZ + ix*PerX);
                    float* pSymData = pStokes->pBaseSto + (iyPerY_p_izPerZ + (Nx_mi_1 - ix)*PerX);
                    CopySymEnergySlice(pOrigData, pSymData, pStokes->ne, SymWithRespectToXax, SymWithRespectToZax);
                }
            }
		}
	}
}

//*************************************************************************

void srTRadIntThickBeam::CopySymEnergySlice(float* pOrigData, float* pSymData, long Ne, char SymWithRespectToXax, char SymWithRespectToZax)
{
	char ChangeSignS2 = !(SymWithRespectToXax && SymWithRespectToZax);
	char ChangeSignS3 = SymWithRespectToXax;
	float *tOrig = pOrigData, *tSym = pSymData;
	for(int ie=0; ie<Ne; ie++)
	{
		*(tSym++) = *(tOrig++); *(tSym++) = *(tOrig++);
		*(tSym++) = ChangeSignS2? -(*(tOrig++)) : *(tOrig++);
		*(tSym++) = ChangeSignS3? -(*(tOrig++)) : *(tOrig++);
	}
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint(srTEXZY EXZY, srTParPrecStokesArb& PrecPar, srTStokes& CurSt)
{
	gAuxPar.EXZY = EXZY;
	gAuxPar.xObsE2 = EXZY.x*EXZY.x;
    gAuxPar.zObsE2 = EXZY.z*EXZY.z;
    gAuxPar.xzObs = EXZY.x*EXZY.z;

	if(PrecPar.MethNo == 1)
	{
		srTStokes StExt;
		if(PrecPar.IntOrFlux == 'i') 
		{
            ComputeStokesAtOneObsPoint_Intens_PrepAandB(gFldArr, 0, 0, gFldArr.Ns, 4, gBottomArrA, gBottomArrB);
            ComputeStokesAtOneObsPoint_Intens_PrepAandB(gFldArr, gFldArr.Ns - 4, 4, 4, gFldArr.Ns - 4, gRightArrA, gRightArrB);

			ComputeStokesAtOneObsPoint_ExternIntens(gFldArr, StExt);
            ComputeStokesAtOneObsPoint_InternIntens_EvenMesh(gFldArr, CurSt);
		}
		//else if(PrecPar.IntOrFlux == 'f') 
		//{
        //	ComputeStokesAtOneObsPoint_ExternFlux(EXZY, gFldArr, CompKey, StExt);
		//}

		//OC//CurSt += StExt; //to uncomment
	}
	else if(PrecPar.MethNo == 2)
	{
	}
}

//*************************************************************************

//void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_Intens_PrepAandB(srTFieldBasedArrays& FldArr, int iStart, int itStart, int CurNs, int CurNst, TComplexD* ArrA, TComplexD* ArrB)
void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_Intens_PrepAandB(srTFieldBasedArrays& FldArr, long long iStart, long long itStart, long long CurNs, long long CurNst, TComplexD* ArrA, TComplexD* ArrB)
{//to modify: take into account triangular alignement !
	//long ns = FldArr.Ns;
	long long ns = FldArr.Ns;
    double sStart = FldArr.sStart;
    double sStep = FldArr.sStep;

	long TotNumCoefForOnePointA = gNumCoefForOnePointA*4; // 4 stokes components

	//for(int it=itStart; it<(itStart + CurNst); it++)
	for(long long it=itStart; it<(itStart + CurNst); it++)
	{
        double st = sStart + it*sStep;

		//long it_CurNs = (it - itStart)*CurNs;
		long long it_CurNs = (it - itStart)*CurNs;
		//for(int i=iStart; i<(iStart + CurNs); i++)
		for(long long i=iStart; i<(iStart + CurNs); i++)
		{
            double s = sStart + i*sStep;

            bool IsConjugate = false;
			//int itAct = it, iAct = i;
			long long itAct = it, iAct = i;
			if(iAct < itAct)
			{
				IsConjugate = true;
                itAct = i; iAct = it;
				double Aux_st = st;
				st = s; s = Aux_st;
			}

			//long Offset1 = it*((2*ns - 1 - it) >> 1) + i;
			//long Offset1 = itAct*(((ns << 1) - 1 - itAct) >> 1) + iAct;
			//long Offset1 = (itAct*((ns << 1) - 1 - itAct) >> 1) + iAct;
			//long OffsetCoefA = Offset1*TotNumCoefForOnePointA;
			//long OffsetCoefB = Offset1*gNumCoefForOnePointB;
 			long long Offset1 = (itAct*((ns << 1) - 1 - itAct) >> 1) + iAct;
			long long OffsetCoefA = Offset1*TotNumCoefForOnePointA;
            long long OffsetCoefB = Offset1*gNumCoefForOnePointB;
			//double s = sStart + i*sStep;
            //double s = sStart + iAct*sStep;

            //long OffsetLocB = it_CurNs + (i - iStart);
            //long OffsetLocA = OffsetLocB << 2;
            long long OffsetLocB = it_CurNs + (i - iStart);
            long long OffsetLocA = OffsetLocB << 2;

            ComputeIntensFuncPartsForInteg2D(s, st, gCoefA + OffsetCoefA, gCoefB + OffsetCoefB, ArrA + OffsetLocA, *(ArrB + OffsetLocB));
			if(IsConjugate)
			{
                (ArrB + OffsetLocB)->Conjugate();
				for(int k=0; k<4; k++) (ArrA + (OffsetLocA + k))->Conjugate();
			}
		}
	}
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_ExternIntens_AddBotLeft(srTFieldBasedArrays& FldArr, srTStokes& CurSt)
{
	//TComplexD ArrA[9*4], ArrB[9];
    //ComputeStokesAtOneObsPoint_Intens_PrepAandB(FldArr, 0, 0, 3, 3, ArrA, ArrB);

	//long Ns = FldArr.Ns;
	//long TwoNs = FldArr.Ns << 1;
	//long ThreeNs = 3*FldArr.Ns;
	//long FourNs = FldArr.Ns << 2;
	long long Ns = FldArr.Ns;
	long long TwoNs = FldArr.Ns << 1;
	long long ThreeNs = 3*FldArr.Ns;
	long long FourNs = FldArr.Ns << 2;

	TComplexD At_Stokes[16];
    ComputeResidTermA_Stokes(gBottomArrA, gBottomArrB, 0, FldArr.sStep, At_Stokes);
    ComputeResidTermA_Stokes(gBottomArrA + FourNs, gBottomArrB + Ns, 0, FldArr.sStep, At_Stokes + 4);
    ComputeResidTermA_Stokes(gBottomArrA + (FourNs << 1), gBottomArrB + TwoNs, 0, FldArr.sStep, At_Stokes + 8);
    ComputeResidTermA_Stokes(gBottomArrA + (3*FourNs), gBottomArrB + ThreeNs, 0, FldArr.sStep, At_Stokes + 12);
    //ComputeResidTermA_Stokes(ArrA + 12, ArrB + 3, FldArr.sStep, At_Stokes + 4);
    //ComputeResidTermA_Stokes(ArrA + 24, ArrB + 6, FldArr.sStep, At_Stokes + 8);
	
	TComplexD Bt[] = {gBottomArrB[0], gBottomArrB[Ns], gBottomArrB[TwoNs], gBottomArrB[ThreeNs]};
	TComplexD At_StokesFin[4];
    ComputeResidTermA_Stokes(At_Stokes, Bt, 0, FldArr.sStep, At_StokesFin);

	double ExpReB0 = exp(gBottomArrB->x);
    TComplexD ExpB(ExpReB0*cos(gBottomArrB->y), ExpReB0*sin(gBottomArrB->y));
	for(int j=0; j<4; j++) At_StokesFin[j] *= ExpB;

	CurSt.s0 += At_StokesFin[0].x;
	CurSt.s1 += At_StokesFin[1].x;
	CurSt.s2 += At_StokesFin[2].x;
	CurSt.s3 += At_StokesFin[3].x;
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_ExternIntens_AddCen(srTFieldBasedArrays& FldArr, char b_or_r, srTStokes& CurSt)
{// b_or_r means 'b' or 'r' (bottom or right respectively)
	//int Ns = FldArr.Ns;
	//int Ns_mi_3 = Ns - 3;
    long long Ns = FldArr.Ns;
	long long Ns_mi_3 = Ns - 3;

    TComplexD ArrA[16];
    TComplexD ArrB[4];
	TComplexD LeftB[3], RightB[3], LeftA[12], RightA[12];
    TComplexD At_Stokes[4];
	srTStokes *LocStokesArr = new srTStokes[Ns]; 
	srTStokes *tLocStokesArr = LocStokesArr;
	//for(int i=0; i<Ns; i++)
	for(long long i=0; i<Ns; i++)
	{
		if(b_or_r == 'b')
		{
			for(int k=0; k<4; k++)
			{
				//long OffsetB = i + Ns*k;
				long long OffsetB = i + Ns*k;
				ArrB[k] = gBottomArrB[OffsetB];

				TComplexD *tgBottomArrA = gBottomArrA + (OffsetB << 2);
				TComplexD *tArrA = ArrA + (k << 2);
				for(int m=0; m<4; m++) *(tArrA++) = *(tgBottomArrA++);
			}
		}
		else
		{
			if(i < 4)
			{
                //long OffsetB0 = (Ns - 4) + i*Ns;
                long long OffsetB0 = (Ns - 4) + i*Ns;
                for(int k=0; k<4; k++)
                {
                    //long OffsetB = OffsetB0 + k;
                    long long OffsetB = OffsetB0 + k;
                    ArrB[k] = gBottomArrB[OffsetB];
                    //ArrA[k] = gBottomArrA[OffsetB << 2];

					//long OffsetA0 = OffsetB << 2;
                    //for(int m=0; m<4; m++) ArrA[k + m] = gBottomArrA[OffsetA0 + m];

                    TComplexD *tgBottomArrA = gBottomArrA + (OffsetB << 2);
                    TComplexD *tArrA = ArrA + (k << 2);
                    for(int m=0; m<4; m++) *(tArrA++) = *(tgBottomArrA++);
                }
			}
			else
			{
                //long OffsetB0 = (i - 4) << 2;
                long long OffsetB0 = (i - 4) << 2;
                for(int k=0; k<4; k++)
                {
                    //long OffsetB = OffsetB0 + k;
                    long long OffsetB = OffsetB0 + k;
                    ArrB[k] = gRightArrB[OffsetB];
                    //ArrA[k] = gRightArrA[OffsetB << 2];

                    TComplexD *tgRightArrA = gRightArrA + (OffsetB << 2);
                    TComplexD *tArrA = ArrA + (k << 2);
                    for(int m=0; m<4; m++) *(tArrA++) = *(tgRightArrA++);
                }
			}
		}
        ComputeResidTermA_Stokes(ArrA, ArrB, 0, FldArr.sStep, At_Stokes);

		if(i < 3)
		{
            LeftB[i] = ArrB[0];
			//int i_4 = i << 2;
			long long i_4 = i << 2;
			for(int k=0; k<4; k++) LeftA[i_4 + k] = At_Stokes[k];
		}
		else if(i >= Ns_mi_3)
		{
			//int i_mi_Ns_mi_3 = i - Ns_mi_3;
			//int i_mi_Ns_mi_3_4 = i_mi_Ns_mi_3 << 2;
			long long i_mi_Ns_mi_3 = i - Ns_mi_3;
            long long i_mi_Ns_mi_3_4 = i_mi_Ns_mi_3 << 2;
            RightB[i_mi_Ns_mi_3] = ArrB[0];
			for(int k=0; k<4; k++) RightA[i_mi_Ns_mi_3_4 + k] = At_Stokes[k];
		}

		int i0 = 0;
		if(b_or_r != 'b') i0 = 3;

		double ExpReB0 = exp(ArrB[i0].x);
		TComplexD ExpB(ExpReB0*cos(ArrB[i0].y),	ExpReB0*sin(ArrB[i0].y));
		for(int j=0; j<4; j++) At_Stokes[j] *= ExpB;

		tLocStokesArr->s0 = 2*(At_Stokes[0].x);
		tLocStokesArr->s1 = 2*(At_Stokes[1].x);
		tLocStokesArr->s2 = 2*(At_Stokes[2].x);
		tLocStokesArr->s3 = 2*(At_Stokes[3].x);
		tLocStokesArr++;
	}

    srTStokes InitDerStokes = DerivWithExpStokes3p(LeftA, LeftB, FldArr.sStep, 0);
    srTStokes FinDerStokes = DerivWithExpStokes3p(RightA, RightB, FldArr.sStep, 2);
	srTStokes LocStokes(0,0,0,0);
	Integrate1DStokesArr(LocStokesArr, Ns, FldArr.sStep, &InitDerStokes, &FinDerStokes, LocStokes);
    if(b_or_r != 'b') LocStokes *= (-1.);

	CurSt += LocStokes;
	if(LocStokesArr != 0) delete[] LocStokesArr;

//to modify: calculate for 4p derivatives !!!
/*
	int LocNs_mi_3 = LocNs - 3;

	for(int i=0; i<LocNs; i++)
	{
		if(b_or_r == 'b') ComputeStokesAtOneObsPoint_Intens_PrepAandB(FldArr, i + 1, 0, 1, 3, ArrA, ArrB);
		else ComputeStokesAtOneObsPoint_Intens_PrepAandB(FldArr, Ns_mi_3, i + 1, 3, 1, ArrA, ArrB);

		ComputeResidTermA_Stokes(ArrA, ArrB, FldArr.sStep, At_Stokes);
        if(b_or_r != 'b')
		{
            TComplexD *tAt_Stokes = At_Stokes;
			for(int k=0; k<4; k++) 
			{
				tAt_Stokes->x = -(tAt_Stokes->x);
				tAt_Stokes->y = -(tAt_Stokes->y);
			}
		}

		if(i < 3)
		{
            LeftB[i] = ArrB[2];
			int i_4 = i << 2;
			for(int k=0; k<4; k++) LeftA[i_4 + k] = At_Stokes[k];
		}
		else if(i >= LocNs_mi_3)
		{
			int i_mi_LocNs_mi_3 = i - LocNs_mi_3;
			int i_mi_LocNs_mi_3_4 = i_mi_LocNs_mi_3 << 2;
            RightB[i_mi_LocNs_mi_3] = ArrB[2];
			for(int k=0; k<4; k++) RightA[i_mi_LocNs_mi_3_4 + k] = At_Stokes[k];
		}

        double ExpReB0 = exp(ArrB[2].x);
        TComplexD ExpB(ExpReB0*cos(ArrB[2].y), ExpReB0*sin(ArrB[2].y));
		At_Stokes[0] *= ExpB;
		At_Stokes[1] *= ExpB;
		At_Stokes[2] *= ExpB;
		At_Stokes[3] *= ExpB;

		tLocStokesArr->s0 = 2*(At_Stokes[0].x);
		tLocStokesArr->s1 = 2*(At_Stokes[1].x);
		tLocStokesArr->s2 = 2*(At_Stokes[2].x);
		tLocStokesArr->s3 = 2*(At_Stokes[3].x);
		tLocStokesArr++;
	}
*/
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_ExternIntens_AddBotRight(srTFieldBasedArrays& FldArr, srTStokes& CurSt)
{
	//long Ns = FldArr.Ns;
	//long Ns_mi_4 = Ns - 4;
	//long TwoNs = FldArr.Ns << 1;
	//long ThreeNs = 3*FldArr.Ns;
	//long FourNs = FldArr.Ns << 2;
	long long Ns = FldArr.Ns;
	long long Ns_mi_4 = Ns - 4;
	long long TwoNs = FldArr.Ns << 1;
	long long ThreeNs = 3*FldArr.Ns;
	long long FourNs = FldArr.Ns << 2;

	TComplexD At_Stokes[16];
	//long OffsetB = Ns_mi_4;
	long long OffsetB = Ns_mi_4;
	for(int i=0; i<4; i++)
	{
        ComputeResidTermA_Stokes(gBottomArrA + (OffsetB << 2), gBottomArrB + OffsetB, 3, FldArr.sStep, At_Stokes + (i << 2));
        OffsetB += Ns;
	}

	TComplexD Bt[] = {gBottomArrB[Ns - 1], gBottomArrB[TwoNs - 1], gBottomArrB[ThreeNs - 1], gBottomArrB[FourNs - 1]};
	TComplexD At_StokesFin[4];
    ComputeResidTermA_Stokes(At_Stokes, Bt, 0, FldArr.sStep, At_StokesFin);

	double ExpReB0 = exp(gBottomArrB[Ns - 1].x);
    TComplexD ExpB(ExpReB0*cos(gBottomArrB[Ns - 1].y), ExpReB0*sin(gBottomArrB[Ns - 1].y));
	for(int j=0; j<4; j++) At_StokesFin[j] *= ExpB;

	CurSt.s0 += (-2.)*At_StokesFin[0].x;
	CurSt.s1 += (-2.)*At_StokesFin[1].x;
	CurSt.s2 += (-2.)*At_StokesFin[2].x;
	CurSt.s3 += (-2.)*At_StokesFin[3].x;
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_ExternIntens_AddTopRight(srTFieldBasedArrays& FldArr, srTStokes& CurSt)
{
	//long Ns = FldArr.Ns;
	//long Ns_mi_4 = Ns - 4;
	//long Four_Ns_mi_4 = Ns_mi_4 << 2;
	long long Ns = FldArr.Ns;
	long long Ns_mi_4 = Ns - 4;
	long long Four_Ns_mi_4 = Ns_mi_4 << 2;

	TComplexD At_Stokes[16];
	//long OffsetB = Four_Ns_mi_4 - 16;
	long long OffsetB = Four_Ns_mi_4 - 16;
	for(int i=0; i<4; i++)
	{
        ComputeResidTermA_Stokes(gRightArrA + (OffsetB << 2), gRightArrB + OffsetB, 3, FldArr.sStep, At_Stokes + (i << 2));
        OffsetB += 4;
	}

	//long AuxOffsetB = Four_Ns_mi_4 - 13;
	long long AuxOffsetB = Four_Ns_mi_4 - 13;
	TComplexD Bt[] = {gRightArrB[AuxOffsetB], gRightArrB[AuxOffsetB + 4], gRightArrB[AuxOffsetB + 8], gRightArrB[AuxOffsetB + 12]};
	TComplexD At_StokesFin[4];
    ComputeResidTermA_Stokes(At_Stokes, Bt, 3, FldArr.sStep, At_StokesFin);

	double ExpReB0 = exp(gRightArrB[Four_Ns_mi_4 - 1].x);
    TComplexD ExpB(ExpReB0*cos(gRightArrB[Four_Ns_mi_4 - 1].y), ExpReB0*sin(gRightArrB[Four_Ns_mi_4 - 1].y));
	for(int j=0; j<4; j++) At_StokesFin[j] *= ExpB;

	CurSt.s0 += At_StokesFin[0].x;
	CurSt.s1 += At_StokesFin[1].x;
	CurSt.s2 += At_StokesFin[2].x;
	CurSt.s3 += At_StokesFin[3].x;
}

//*************************************************************************

//void srTRadIntThickBeam::ComputeStokesIntensFuncAtPoint(long is, long ist, srTFieldBasedArrays& FldArr, srTStokesC& A, srTStokesC& B, srTStokesC& C)
//{
//}

//*************************************************************************

//void srTRadIntThickBeam::ComputeStokesFluxFuncAtPoint(long is, long ist, srTFieldBasedArrays& FldArr, srTStokesC& A, srTStokesC& B, srTStokesC& C)
//{
//}

//*************************************************************************

void srTRadIntThickBeam::ComputeExpCoefXZArraysForInteg2D_EvenMesh(double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB)
{
	TComplexD *pA = ArrA, *pB = ArrB;
	//long Ns = FldArr.Ns;
	long long Ns = FldArr.Ns;
	long TotNumCoefForOnePointA = gNumCoefForOnePointA*4; //for 4 stokes components

	//long AuxCount = 0;
	long long AuxCount = 0;
	//for(long ist=0; ist<Ns; ist++)
	for(long long ist=0; ist<Ns; ist++)
	{
		//for(long is=ist; is<Ns; is++)
		for(long long is=ist; is<Ns; is++)
		{
            ComputeExpCoefForOneObsPoint(is, ist, yObs, eObs, FldArr, pA, pB);
            pA += TotNumCoefForOnePointA;
            pB += gNumCoefForOnePointB;

			AuxCount++; //for debugging
		}
	}
}

//*************************************************************************

//void srTRadIntThickBeam::ComputeExpCoefForOneObsPoint(long is, long ist, double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB)
void srTRadIntThickBeam::ComputeExpCoefForOneObsPoint(long long is, long long ist, double yObs, double eObs, srTFieldBasedArrays& FldArr, TComplexD* ArrA, TComplexD* ArrB)
{
	double s = FldArr.sStart + is*FldArr.sStep;
	double us = FldArr.sStart + ist*FldArr.sStep;
	double Inv_y_mi_s = 1./(yObs - s);
	double Inv_y_mi_us = 1./(yObs - us);

	double k = gAuxPar.k_d_e*eObs;
	double Half_k = 0.5*k;
	double k_d_y_mi_s = k*Inv_y_mi_s, k_d_y_mi_us = k*Inv_y_mi_us;

	double xeq0 = *(FldArr.XArr + is), uxeq0 = *(FldArr.XArr + ist), zeq0 = *(FldArr.ZArr + is), uzeq0 = *(FldArr.ZArr + ist);
	double xpeq0 = *(FldArr.BtxArr + is), uxpeq0 = *(FldArr.BtxArr + ist), zpeq0 = *(FldArr.BtzArr + is), uzpeq0 = *(FldArr.BtzArr + ist);
	double x1 = *(FldArr.X1Arr + is), ux1 = *(FldArr.X1Arr + ist), z1 = *(FldArr.Z1Arr + is), uz1 = *(FldArr.Z1Arr + ist);
	double x1p = *(FldArr.X1pArr + is), ux1p = *(FldArr.X1pArr + ist), z1p = *(FldArr.Z1pArr + is), uz1p = *(FldArr.Z1pArr + ist);
	double x2 = *(FldArr.X2Arr + is), ux2 = *(FldArr.X2Arr + ist), z2 = *(FldArr.Z2Arr + is), uz2 = *(FldArr.Z2Arr + ist);
	double x2p = *(FldArr.X2pArr + is), ux2p = *(FldArr.X2pArr + ist), z2p = *(FldArr.Z2pArr + is), uz2p = *(FldArr.Z2pArr + ist);

	double I00 = (*(FldArr.IntBtxE2Arr + is)) + (*(FldArr.IntBtzE2Arr + is)), uI00 = (*(FldArr.IntBtxE2Arr + ist)) + (*(FldArr.IntBtzE2Arr + ist));
	double Ix01 = *(FldArr.IntX01Arr + is), uIx01 = *(FldArr.IntX01Arr + ist), Iz01 = *(FldArr.IntZ01Arr + is), uIz01 = *(FldArr.IntZ01Arr + ist);
	double Ix02 = *(FldArr.IntX02Arr + is), uIx02 = *(FldArr.IntX02Arr + ist), Iz02 = *(FldArr.IntZ02Arr + is), uIz02 = *(FldArr.IntZ02Arr + ist);
	double Ix11 = *(FldArr.IntX11Arr + is), uIx11 = *(FldArr.IntX11Arr + ist), Iz11 = *(FldArr.IntZ11Arr + is), uIz11 = *(FldArr.IntZ11Arr + ist);
	double Ix12 = *(FldArr.IntX12Arr + is), uIx12 = *(FldArr.IntX12Arr + ist), Iz12 = *(FldArr.IntZ12Arr + is), uIz12 = *(FldArr.IntZ12Arr + ist);
	double Ix22 = *(FldArr.IntX22Arr + is), uIx22 = *(FldArr.IntX22Arr + ist), Iz22 = *(FldArr.IntZ22Arr + is), uIz22 = *(FldArr.IntZ22Arr + ist);

	double x1dr = x1*Inv_y_mi_s, ux1dr = ux1*Inv_y_mi_us, z1dr = z1*Inv_y_mi_s, uz1dr = uz1*Inv_y_mi_us;
	double x2dr = x2*Inv_y_mi_s, ux2dr = ux2*Inv_y_mi_us, z2dr = z2*Inv_y_mi_s, uz2dr = uz2*Inv_y_mi_us;

	double betx0 = xpeq0 + xeq0*Inv_y_mi_s, ubetx0 = uxpeq0 + uxeq0*Inv_y_mi_us, betz0 = zpeq0 + zeq0*Inv_y_mi_s, ubetz0 = uzpeq0 + uzeq0*Inv_y_mi_us;
	double betx1 = x1p + x1dr, ubetx1 = ux1p + ux1dr, betz1 = z1p + z1dr, ubetz1 = uz1p + uz1dr;
	double betx2 = x2p + x2dr, ubetx2 = ux2p + ux2dr, betz2 = z2p + z2dr, ubetz2 = uz2p + uz2dr;

	double dbetx00 = betx0, udbetx00 = ubetx0, dbetz00 = betz0, udbetz00 = ubetz0;
	double dbet01 = -Inv_y_mi_s, udbet01 = -Inv_y_mi_us;

	double dPhi0000 = Half_k*(I00 - uI00 + (xeq0*xeq0 + zeq0*zeq0)*Inv_y_mi_s - (uxeq0*uxeq0 + uzeq0*uzeq0)*Inv_y_mi_us + (s - us)*gAuxPar.GammaEm2);
	double dPhi0010 = uxeq0*k_d_y_mi_us - xeq0*k_d_y_mi_s;
	double dPhi0001 = uzeq0*k_d_y_mi_us - zeq0*k_d_y_mi_s;
	double dPhi0020 = 0.5*(k_d_y_mi_s - k_d_y_mi_us);
	double dPhi0002 = dPhi0020;
    double dPhi0100 = -2*dPhi0000;
	double dPhi0110 = -dPhi0010;
	double dPhi0101 = -dPhi0001;
	//double dPhi0200 = 3*dPhi0000;
	//double dPhi0210 = dPhi0010;
	//double dPhi0201 = dPhi0001;

	double dalpx0000 = k*(Ix01 + x1dr*xeq0 - (uIx01 + ux1dr*uxeq0)), dalpz0000 = k*(Iz01 + z1dr*zeq0 - (uIz01 + uz1dr*uzeq0));
	double dalpx0010 = k*(ux1dr - x1dr), dalpz0001 = k*(uz1dr - z1dr);
	double dalpx01 = -dalpx0000, dalpz01 = -dalpz0000;
	double dalpx02 = dalpx0000, dalpz02 = dalpz0000;
	double ddeltx0000 = k*(Ix02 + x2dr*xeq0 - (uIx02 + ux2dr*uxeq0)), ddeltz0000 = k*(Iz02 + z2dr*zeq0 - (uIz02 + uz2dr*uzeq0));
	double ddeltx0010 = k*(ux2dr - x2dr), ddeltz0001 = k*(uz2dr - z2dr);
	double ddeltx01 = -ddeltx0000, ddeltz01 = -ddeltz0000;
	double ddeltx02 = ddeltx0000, ddeltz02 = ddeltz0000;
	double depsx = Half_k*(Ix11 + x1*x1dr - (uIx11 + ux1*ux1dr)), depsz = Half_k*(Iz11 + z1*z1dr - (uIz11 + uz1*uz1dr));
	double dzetx = k*(Ix12 + x1dr*x2 - (uIx12 + ux1dr*ux2)), dzetz = k*(Iz12 + z1dr*z2 - (uIz12 + uz1dr*uz2));
	double detax = Half_k*(Ix22 + x2dr*x2 - (uIx22 + ux2dr*ux2)), detaz = Half_k*(Iz22 + z2dr*z2 - (uIz22 + uz2dr*uz2));

	double a0000[4], a0010[4], a0001[4], a0002[4], a0011[4], a0020[4], a0100[4], a0110[4], a0101[4], a0200[4], a0210[4], a0201[4];
	double bx0000[4], bx0010[4], bx0001[4], cx0000[4], cx0010[4], cx0001[4], bz0000[4], bz0010[4], bz0001[4], cz0000[4], cz0010[4], cz0001[4];
	double bx01[4], bx02[4], cx01[4], cx02[4], bz01[4], bz02[4], cz01[4], cz02[4];
	double dx[4], ex[4], fx[4], dz[4], ez[4], fz[4], m[4], n[4], v[4], w[4];

	double dbet01_ubetx1 = dbet01*ubetx1, betx1_udbet01 = betx1*udbet01, dbet01_ubetx2 = dbet01*ubetx2, betx2_udbet01 = betx2*udbet01;
	double dbet01_ubetz1 = dbet01*ubetz1, betz1_udbet01 = betz1*udbet01, dbet01_ubetz2 = dbet01*ubetz2, betz2_udbet01 = betz2*udbet01;
	double dbetz00_udbet01 = dbetz00*udbet01, dbet01_udbetz00 = dbet01*udbetz00, dbetx00_udbet01 = dbetx00*udbet01, dbet01_udbetx00 = dbet01*udbetx00;
	double dbet01_ubetz0 = dbet01*ubetz0, betz0_udbet01 = betz0*udbet01, dbet01_ubetx0 = dbet01*ubetx0, betx0_udbet01 = betx0*udbet01;
	double dbetx00_udbetx00 = dbetx00*udbetx00, dbetz00_udbetz00 = dbetz00*udbetz00;
	double dbetx00_ubetx0 = dbetx00*ubetx0, dbetz00_ubetz0 = dbetz00*ubetz0, betx0_udbetx00 = betx0*udbetx00, betz0_udbetz00 = betz0*udbetz00;
	double betx0_ubetx0 = betx0*ubetx0, betz0_ubetz0 = betz0*ubetz0;

	a0000[0] = dbetx00_udbetx00 + dbetz00_udbetz00;
    a0010[0] = dbetx00_udbet01 + dbet01_udbetx00;
    a0001[0] = dbetz00_udbet01 + dbet01_udbetz00;
    a0020[0] = dbet01*udbet01;
    a0011[0] = 0;
    a0002[0] = a0020[0];
    a0100[0] = -dbetx00_ubetx0 - dbetz00_ubetz0 - betx0_udbetx00 - betz0_udbetz00;
    a0110[0] = -dbet01_ubetx0 - betx0_udbet01;
    a0101[0] = -dbet01_ubetz0 - betz0_udbet01;
    a0200[0] = betx0_ubetx0 + dbetx00_ubetx0 + betz0_ubetz0 + dbetz00_ubetz0 + betx0_udbetx00 + betz0_udbetz00;
    a0210[0] = -a0110[0];
    a0201[0] = -a0101[0];
    bx0000[0] = dbetx00*ubetx1 + betx1*udbetx00;
    bx0010[0] = dbet01_ubetx1 + betx1_udbet01;
	bx0001[0] = 0;
    bx01[0] = -betx1*ubetx0 - betx0*ubetx1;
    bx02[0] = -bx01[0];
    cx0000[0] = dbetx00*ubetx2 + betx2*udbetx00;
    cx0010[0] = dbet01_ubetx2 + betx2_udbet01;
    cx0001[0] = 0;
    cx01[0] = -betx2*ubetx0 - betx0*ubetx2;
    cx02[0] = -cx01[0];
    dx[0] = betx1*ubetx1;
    ex[0] = betx2*ubetx2;
    fx[0] = betx2*ubetx1 + betx1*ubetx2;
    bz0000[0] = dbetz00*ubetz1 + betz1*udbetz00;
    bz0010[0] = 0;
    bz0001[0] = dbet01_ubetz1 + betz1_udbet01;
    bz01[0] = -betz1*ubetz0 - betz0*ubetz1;
    bz02[0] = -bz01[0];
    cz0000[0] = dbetz00*ubetz2 + betz2*udbetz00;
	cz0010[0] = 0;
    cz0001[0] = dbet01_ubetz2 + betz2_udbet01;
    cz01[0] = -betz2*ubetz0 - betz0*ubetz2;
    cz02[0] = -cz01[0];
	dz[0] = betz1*ubetz1;
    ez[0] = betz2*ubetz2;
    fz[0] = betz2*ubetz1 + betz1*ubetz2;
	m[0] = 0; n[0] = 0; v[0] = 0; w[0] = 0;

    a0000[1] = dbetx00_udbetx00 - dbetz00_udbetz00;
    a0010[1] = a0010[0];
    a0001[1] = -a0001[0];
    a0020[1] = a0020[0];
    a0011[1] = 0;
    a0002[1] = -a0020[1];
    a0100[1] = -dbetx00_ubetx0 + dbetz00_ubetz0 - betx0_udbetx00 + betz0_udbetz00;
    a0110[1] = a0110[0];
    a0101[1] = -a0101[0];
    a0200[1] = betx0_ubetx0 + dbetx00_ubetx0 - betz0_ubetz0 - dbetz00_ubetz0 + betx0_udbetx00 - betz0_udbetz00;
    a0210[1] = -a0110[1];
    a0201[1] = -a0101[1];
    bx0000[1] = bx0000[0];
    bx0010[1] = bx0010[0];
    bx0001[1] = 0;
	bx01[1] = bx01[0];
    bx02[1] = -bx01[1];
    cx0000[1] = cx0000[0];
    cx0010[1] = cx0010[0];
    cx0001[1] = 0;
	cx01[1] = cx01[0];
    cx02[1] = -cx01[1];
    dx[1] = dx[0];
    ex[1] = ex[0];
    fx[1] = fx[0];
    bz0000[1] = -bz0000[0];
	bz0010[1] = 0;
    bz0001[1] = -bz0001[0];
	bz01[1] = -bz01[0];
    bz02[1] = -bz01[1];
    cz0000[1] = -cz0000[0];
    cz0010[1] = 0;
    cz0001[1] = -cz0001[0];
    cz01[1] = -cz01[0];
    cz02[1] = -cz01[1];
	dz[1] = -dz[0];
    ez[1] = -ez[0];
    fz[1] = -fz[0];
	m[1] = 0; n[1] = 0; v[1] = 0; w[1] = 0;

	double dbetz00_udbetx00 = dbetz00*udbetx00, dbetx00_udbetz00 = dbetx00*udbetz00;
	double dbetz00_ubetx0 = dbetz00*ubetx0, dbetx00_ubetz0 = dbetx00*ubetz0, betz0_udbetx00 = betz0*udbetx00, betx0_udbetz00 = betx0*udbetz00;
	double betz0_ubetx0 = betz0*ubetx0, betx0_ubetz0 = betx0*ubetz0;
	double dbetz00_ubetx1 = dbetz00*ubetx1, betx1_udbetz00 = betx1*udbetz00;
	double betz0_ubetx1 = betz0*ubetx1, betx1_ubetz0 = betx1*ubetz0, dbetz00_ubetx2 = dbetz00*ubetx2, betx2_udbetz00 = betx2*udbetz00;
	double betz0_ubetx2 = betz0*ubetx2, betx2_ubetz0 = betx2*ubetz0, dbetx00_ubetz1 = dbetx00*ubetz1, betz1_udbetx00 = betz1*udbetx00;
	double betz1_ubetx0 = betz1*ubetx0, betx0_ubetz1 = betx0*ubetz1, dbetx00_ubetz2 = dbetx00*ubetz2, betz2_udbetx00 = betz2*udbetx00;
	double betz2_ubetx0 = betz2*ubetx0, betx0_ubetz2 = betx0*ubetz2;
	double betz1_ubetx1 = betz1*ubetx1, betx1_ubetz1 = betx1*ubetz1, betz2_ubetx1 = betz2*ubetx1, betx1_ubetz2 = betx1*ubetz2;
	double betz1_ubetx2 = betz1*ubetx2, betx2_ubetz1 = betx2*ubetz1, betz2_ubetx2 = betz2*ubetx2, betx2_ubetz2 = betx2*ubetz2;

	a0000[2] = -dbetz00_udbetx00 - dbetx00_udbetz00;
    a0010[2] = -a0001[0];
    a0001[2] = -a0010[0];
    a0020[2] = 0;
    a0011[2] = -2*a0020[0];
    a0002[2] = 0;
	a0100[2] = dbetz00_ubetx0 + dbetx00_ubetz0 + betz0_udbetx00 + betx0_udbetz00;
    a0110[2] = -a0101[0];
    a0101[2] = -a0110[0];
    a0200[2] = -betz0_ubetx0 - dbetz00_ubetx0 - betx0_ubetz0 - dbetx00_ubetz0 - betz0_udbetx00 - betx0_udbetz00;
    a0210[2] = -a0110[2];
    a0201[2] = -a0101[2];
    bx0000[2] = -dbetz00_ubetx1 - betx1_udbetz00;
    bx0010[2] = 0;
    bx0001[2] = -dbet01_ubetx1 - betx1_udbet01;
	bx01[2] = betz0_ubetx1 + betx1_ubetz0;
    bx02[2] = -bx01[2];
    cx0000[2] = -dbetz00_ubetx2 - betx2_udbetz00;
    cx0010[2] = 0;
    cx0001[2] = -dbet01_ubetx2 - betx2_udbet01;
	cx01[2] = betz0_ubetx2 + betx2_ubetz0;
    cx02[2] = -cx01[2];
    dx[2] = 0; ex[2] = 0; fx[2] = 0;
    bz0000[2] = -dbetx00_ubetz1 - betz1_udbetx00;
    bz0010[2] = -dbet01_ubetz1 - betz1_udbet01;
	bz0001[2] = 0;
    bz01[2] = betz1_ubetx0 + betx0_ubetz1;
    bz02[2] = -bz01[2];
    cz0000[2] = -dbetx00_ubetz2 - betz2_udbetx00;
    cz0010[2] = -dbet01_ubetz2 - betz2_udbet01;
	cz0001[2] = 0;
    cz01[2] = betz2_ubetx0 + betx0_ubetz2;
    cz02[2] = -cz01[2];
    dz[2] = 0; ez[2] = 0; fz[2] = 0;
    m[2] = -betz1_ubetx1 - betx1_ubetz1;
    n[2] = -betz2_ubetx1 - betx1_ubetz2;
    v[2] = -betz1_ubetx2 - betx2_ubetz1;
    w[2] = -betz2_ubetx2 - betx2_ubetz2;

	a0000[3] = dbetz00_udbetx00 - dbetx00_udbetz00; //to multiply s4 by i
    a0010[3] = dbetz00_udbet01 - dbet01_udbetz00;
    a0001[3] = -dbetx00_udbet01 + dbet01_udbetx00;
    a0020[3] = 0;
    a0011[3] = 0;
    a0002[3] = 0;
    a0100[3] = -dbetz00_ubetx0 + dbetx00_ubetz0 - betz0_udbetx00 + betx0_udbetz00;
    a0110[3] = dbet01_ubetz0 - betz0_udbet01;
    a0101[3] = -dbet01_ubetx0 + betx0_udbet01;
    a0200[3] = betz0_ubetx0 + dbetz00_ubetx0 - betx0_ubetz0 - dbetx00_ubetz0 + betz0_udbetx00 - betx0_udbetz00;
    a0210[3] = -a0110[3];
    a0201[3] = -a0101[3];
    bx0000[3] = dbetz00_ubetx1 - betx1_udbetz00;
    bx0010[3] = 0;
    bx0001[3] = dbet01_ubetx1 - betx1_udbet01;
    bx01[3] = -betz0_ubetx1 + betx1_ubetz0;
    bx02[3] = -bx01[3];
    cx0000[3] = dbetz00_ubetx2 - betx2_udbetz00;
    cx0010[3] = 0;
    cx0001[3] = dbet01_ubetx2 - betx2_udbet01;
    cx01[3] = -betz0_ubetx2 + betx2_ubetz0;
    cx02[3] = -cx01[3];
    dx[3] = 0; ex[3] = 0; fx[3] = 0;
	bz0000[3] = -dbetx00_ubetz1 + betz1_udbetx00;
    bz0010[3] = -dbet01_ubetz1 + betz1_udbet01;
    bz0001[3] = 0;
    bz01[3] = -betz1_ubetx0 + betx0_ubetz1;
    bz02[3] = -bz01[3];
    cz0000[3] = -dbetx00_ubetz2 + betz2_udbetx00;
    cz0010[3] = -dbet01_ubetz2 + betz2_udbet01;
    cz0001[3] = 0;
    cz01[3] = -betz2_ubetx0 + betx0_ubetz2;
    cz02[3] = -cz01[3];
	dz[3] = 0; ez[3] = 0; fz[3] = 0;
    m[3] = betz1_ubetx1 - betx1_ubetz1;
    n[3] = betz2_ubetx1 - betx1_ubetz2;
    v[3] = betz1_ubetx2 - betx2_ubetz1;
    w[3] = betz2_ubetx2 - betx2_ubetz2;

	TComplexD ImagI(0, 1), Half_i(0, 0.5), Zero(0, 0);

	TComplexD BetZ1(gAuxPar.BetZ, -detaz);
	TComplexD InvBetZ1 = 1./BetZ1;
    TComplexD HalfInvBetZ1 = 0.5*InvBetZ1;

	TComplexD deltz1000 = ddeltz0000 + gAuxPar.Aux_deltz1000;
	double deltz1001 = ddeltz0001;
	double deltz11 = ddeltz01;
	double deltz12 = ddeltz02;
	TComplexD AlpZ1(gAuxPar.AlpZ, -0.5*dzetz);
	TComplexD AlpZ1_d_BetZ1 = AlpZ1*InvBetZ1;
	TComplexD AlpZ1_d_BetZ1E2 = AlpZ1_d_BetZ1*InvBetZ1;
	TComplexD AlpZ1E2_d_BetZ1 = AlpZ1*AlpZ1_d_BetZ1;
	TComplexD AuxGamZ(gAuxPar.GamZ, -depsz);
	TComplexD GamZ1 = AuxGamZ - AlpZ1E2_d_BetZ1;
	TComplexD InvGamZ1 = 1./GamZ1;
	TComplexD HalfInvGamZ1 = 0.5*InvGamZ1;

	TComplexD ThetXpZ1 = gAuxPar.ThetXpZ - (gAuxPar.ThetXpZp*AlpZ1_d_BetZ1);
	TComplexD ThetXZ1 = gAuxPar.ThetXZ - (gAuxPar.ThetXZp*AlpZ1_d_BetZ1);
	TComplexD alpz1000 = dalpz0000 - ((deltz1000*AlpZ1_d_BetZ1) + gAuxPar.Aux_alpz1000);
	TComplexD alpz1001 = dalpz0001 - (deltz1001*AlpZ1_d_BetZ1);
	TComplexD alpz11 = dalpz01 - (deltz11*AlpZ1_d_BetZ1);
	TComplexD alpz12 = dalpz02 - (deltz12*AlpZ1_d_BetZ1);

	TComplexD ThetXpZ1_d_GamZ1 = ThetXpZ1*InvGamZ1;
	TComplexD ThetXpZ1E2_d_GamZ1 = ThetXpZ1*ThetXpZ1_d_GamZ1;
	TComplexD ThetXpZp_d_BetZ1 = gAuxPar.ThetXpZp*InvBetZ1;
	TComplexD ThetXpZpE2_d_BetZ1 = gAuxPar.ThetXpZp*ThetXpZp_d_BetZ1;
	TComplexD AuxBetX1(gAuxPar.BetX, -detax);
	TComplexD BetX1 = AuxBetX1 - (ThetXpZ1E2_d_GamZ1 + ThetXpZpE2_d_BetZ1);
	TComplexD InvBetX1 = 1./BetX1;
	TComplexD HalfInvBetX1 = 0.5*InvBetX1;

	TComplexD AuxAlpX1(gAuxPar.AlpX, -0.5* dzetx);
	TComplexD AlpX1 = AuxAlpX1 - ((ThetXZ1*ThetXpZ1_d_GamZ1) + (gAuxPar.ThetXZp*ThetXpZp_d_BetZ1));
	TComplexD deltx1000 = (ddeltx0000 - gAuxPar.Aux_deltx1000) - ((alpz1000*ThetXpZ1_d_GamZ1) + (deltz1000*ThetXpZp_d_BetZ1));
	double deltx1010 = ddeltx0010;
	TComplexD deltx1001 = Zero - ((alpz1001*ThetXpZ1_d_GamZ1) + (deltz1001*ThetXpZp_d_BetZ1));
	TComplexD deltx11 = ddeltx01 - ThetXpZ1_d_GamZ1*alpz11 - ThetXpZp_d_BetZ1*deltz11;
	TComplexD deltx12 = ddeltx02 - ThetXpZ1_d_GamZ1*alpz12 - ThetXpZp_d_BetZ1*deltz12;

	TComplexD ThetXpZ1E2_d_GamZ1E2 = ThetXpZ1_d_GamZ1*ThetXpZ1_d_GamZ1;
	TComplexD ThetXpZ1_d_GamZ1E2 = ThetXpZ1_d_GamZ1*InvGamZ1;
	TComplexD ThetXpZpE2_d_BetZ1E2 = ThetXpZp_d_BetZ1*ThetXpZp_d_BetZ1;
	TComplexD ThetXpZp_d_BetZ1E2 = ThetXpZp_d_BetZ1*InvBetZ1;
    TComplexD ThetXZp_d_BetZ1 = gAuxPar.ThetXZp*InvBetZ1;
    TComplexD ThetXZ1_d_GamZ1 = ThetXZ1*InvGamZ1;
	TComplexD Two_ThetXpZ1_ThetXZ1_d_GamZ1E2 = (2*ThetXZ1_d_GamZ1)*ThetXpZ1_d_GamZ1;
	TComplexD Two_ThetXpZp_ThetXZp_d_BetZ1E2 = (2*ThetXpZp_d_BetZ1)*ThetXZp_d_BetZ1;

	TComplexD alpz1000_d_GamZ1 = alpz1000*InvGamZ1;
	TComplexD alpz1000_ThetXpZ1_d_GamZ1E2 = alpz1000_d_GamZ1*ThetXpZ1_d_GamZ1;
    TComplexD deltz1000_d_BetZ1 = deltz1000*InvBetZ1;
    TComplexD deltz1000_ThetXpZp_d_BetZ1E2 = deltz1000_d_BetZ1*ThetXpZp_d_BetZ1;
	TComplexD alpz1001_d_GamZ1 = alpz1001*InvGamZ1;
	TComplexD alpz1001_ThetXpZ1_d_GamZ1E2 = alpz1001*ThetXpZ1_d_GamZ1E2;
	TComplexD deltz1001_d_BetZ1 = deltz1001*InvBetZ1;
	TComplexD deltz1001_ThetXpZp_d_BetZ1E2 = ThetXpZp_d_BetZ1*deltz1001_d_BetZ1;
	TComplexD deltz11_d_BetZ1 = deltz11*InvBetZ1;
	TComplexD deltz11_ThetXpZp_d_BetZ1E2 = ThetXpZp_d_BetZ1*deltz11_d_BetZ1;
	TComplexD alpz11_d_GamZ1 = alpz11*InvGamZ1;
	TComplexD alpz11_ThetXpZ1_d_GamZ1E2 = ThetXpZ1_d_GamZ1E2*alpz11;
	TComplexD deltz12_d_BetZ1 = deltz12*InvBetZ1;
	TComplexD deltz12_ThetXpZp_d_BetZ1E2 = ThetXpZp_d_BetZ1*deltz12_d_BetZ1;
	TComplexD alpz12_d_GamZ1 = alpz12*InvGamZ1;
	TComplexD alpz12_ThetXpZ1_d_GamZ1E2 = ThetXpZ1_d_GamZ1E2*alpz12;

	TComplexD AlpX1_d_BetX1 = AlpX1*InvBetX1;
	TComplexD AlpX1E2_d_BetX1 = AlpX1*AlpX1_d_BetX1;
	TComplexD ThetXZpE2_d_BetZ1 = gAuxPar.ThetXZp*ThetXZp_d_BetZ1;
	TComplexD ThetXZ1E2_d_GamZ1 = ThetXZ1*ThetXZ1_d_GamZ1;

	TComplexD Aux_GamX1(gAuxPar.GamX, -depsx);
    TComplexD GamX1 = Aux_GamX1 - AlpX1E2_d_BetX1 - ThetXZpE2_d_BetZ1 - ThetXZ1E2_d_GamZ1;
	TComplexD InvGamX1 = 1./GamX1;
	TComplexD HalfInvGamX1 = 0.5*InvGamX1;

	TComplexD alpx1000 = dalpx0000 + gAuxPar.Aux_alpx1000 - (ThetXZ1_d_GamZ1*alpz1000) - (AlpX1_d_BetX1*deltx1000) - (ThetXZp_d_BetZ1*deltz1000);
	TComplexD alpx1010 = dalpx0010 - (AlpX1_d_BetX1*deltx1010);
	TComplexD alpx1001 = ((ThetXZp_d_BetZ1*(-deltz1001)) - (ThetXZ1_d_GamZ1*alpz1001)) - (AlpX1_d_BetX1*deltx1001);
	TComplexD alpx11 = ((dalpx01 - (ThetXZ1_d_GamZ1*alpz11)) - (AlpX1_d_BetX1*deltx11)) - (ThetXZp_d_BetZ1*deltz11);
	TComplexD alpx12 = ((dalpx02 - (ThetXZ1_d_GamZ1*alpz12)) - (AlpX1_d_BetX1*deltx12)) - (ThetXZp_d_BetZ1*deltz12);

	TComplexD AlpX1E2_d_BetX1E2 = AlpX1_d_BetX1*AlpX1_d_BetX1;
	TComplexD ThetXZ1E2_d_GamZ1E2 = ThetXZ1_d_GamZ1*ThetXZ1_d_GamZ1;
	TComplexD ThetXZpE2_d_BetZ1E2 = ThetXZp_d_BetZ1*ThetXZp_d_BetZ1;
	TComplexD deltx1000_d_BetX1 = deltx1000*InvBetX1;
	TComplexD deltx1000_AlpX1_d_BetX1E2 = deltx1000_d_BetX1*AlpX1_d_BetX1;
    TComplexD deltz1000_ThetXZp_d_BetZ1E2 = deltz1000_d_BetZ1*ThetXZp_d_BetZ1;
    TComplexD alpz1000_ThetXZ1_d_GamZ1E2 = alpz1000_d_GamZ1*ThetXZ1_d_GamZ1;
	TComplexD deltx1010_d_BetX1 = deltx1010*InvBetX1;
    TComplexD deltx1010_AlpX1_d_BetX1E2 = deltx1010_d_BetX1*AlpX1_d_BetX1;
	TComplexD deltx1001_d_BetX1 = deltx1001*InvBetX1;
    TComplexD deltx1001_AlpX1_d_BetX1E2 = deltx1001_d_BetX1*AlpX1_d_BetX1;
	TComplexD deltz1001_ThetXZp_d_BetZ1E2 = deltz1001_d_BetZ1*ThetXZp_d_BetZ1;
	TComplexD alpz1001_ThetXZ1_d_GamZ1E2 = alpz1001_d_GamZ1*ThetXZ1_d_GamZ1;
    TComplexD deltx11_d_BetX1 = deltx11*InvBetX1;
	TComplexD deltx11_AlpX1_d_BetX1E2 = deltx11_d_BetX1*AlpX1_d_BetX1;
	TComplexD deltz11_ThetXZp_d_BetZ1E2 = deltz11_d_BetZ1*ThetXZp_d_BetZ1;
	TComplexD alpz11_ThetXZ1_d_GamZ1E2 = alpz11_d_GamZ1*ThetXZ1_d_GamZ1;
	TComplexD deltx12_d_BetX1 = deltx12*InvBetX1;
	TComplexD deltx12_AlpX1_d_BetX1E2 = deltx12_d_BetX1*AlpX1_d_BetX1;
	TComplexD deltz12_ThetXZp_d_BetZ1E2 = deltz12_d_BetZ1*ThetXZp_d_BetZ1;
    TComplexD alpz12_ThetXZ1_d_GamZ1E2 = alpz12_d_GamZ1*ThetXZ1_d_GamZ1;
	TComplexD alpx1000_d_GamX1 = alpx1000*InvGamX1;
	TComplexD alpx1000E2_d_GamX1E2 = alpx1000_d_GamX1*alpx1000_d_GamX1;
	TComplexD alpz1000E2_d_GamZ1E2 = alpz1000_d_GamZ1*alpz1000_d_GamZ1;
	TComplexD deltx1000E2_d_BetX1E2 = deltx1000_d_BetX1*deltx1000_d_BetX1;
	TComplexD deltz1000E2_d_BetZ1E2 = deltz1000_d_BetZ1*deltz1000_d_BetZ1;
	TComplexD alpx1010_d_GamX1 = alpx1010*InvGamX1;
	TComplexD alpx1000_alpx1010_d_GamX1E2 = alpx1000_d_GamX1*alpx1010_d_GamX1;
	TComplexD deltx1000_deltx1010_d_BetX1E2 = deltx1000_d_BetX1*deltx1010_d_BetX1;
	TComplexD alpx1001_d_GamX1 = alpx1001*InvGamX1;
	TComplexD alpx1000_alpx1001_d_GamX1E2 = alpx1000_d_GamX1*alpx1001_d_GamX1;
	TComplexD alpz1000_alpz1001_d_GamZ1E2 = alpz1000_d_GamZ1*alpz1001_d_GamZ1;
	TComplexD deltx1000_deltx1001_d_BetX1E2 = deltx1000_d_BetX1*deltx1001_d_BetX1;
	TComplexD deltz1000_deltz1001_d_BetZ1E2 = deltz1000_d_BetZ1*deltz1001_d_BetZ1;
	TComplexD alpx1010E2_d_GamX1E2 = alpx1010_d_GamX1*alpx1010_d_GamX1;
	TComplexD deltx1010E2_d_BetX1E2 = deltx1010_d_BetX1*deltx1010_d_BetX1;
	TComplexD alpx1001_alpx1010_d_GamX1E2 = alpx1001_d_GamX1*alpx1010_d_GamX1;
	TComplexD deltx1001_deltx1010_d_BetX1E2 = deltx1001_d_BetX1*deltx1010_d_BetX1;
	TComplexD alpx1001E2_d_GamX1E2 = alpx1001_d_GamX1*alpx1001_d_GamX1;
	TComplexD alpz1001E2_d_GamZ1E2 = alpz1001_d_GamZ1*alpz1001_d_GamZ1;
	TComplexD deltx1001E2_d_BetX1E2 = deltx1001_d_BetX1*deltx1001_d_BetX1;
	TComplexD deltz1001E2_d_BetZ1E2 = deltz1001_d_BetZ1*deltz1001_d_BetZ1;
	TComplexD alpx11_d_GamX1 = alpx11*InvGamX1;
	TComplexD alpx1000_alpx11_d_GamX1E2 = alpx1000_d_GamX1*alpx11_d_GamX1;
	TComplexD alpz1000_alpz11_d_GamZ1E2 = alpz1000_d_GamZ1*alpz11_d_GamZ1;
	TComplexD deltx1000_deltx11_d_BetX1E2 = deltx1000_d_BetX1*deltx11_d_BetX1;
	TComplexD deltz1000_deltz11_d_BetZ1E2 = deltz1000_d_BetZ1*deltz11_d_BetZ1;
	TComplexD alpx1010_alpx11_d_GamX1E2 = alpx1010_d_GamX1*alpx11_d_GamX1;
	TComplexD deltx1010_deltx11_d_BetX1E2 = deltx1010_d_BetX1*deltx11_d_BetX1;
	TComplexD alpx1001_alpx11_d_GamX1E2 = alpx1001_d_GamX1*alpx11_d_GamX1;
	TComplexD alpz1001_alpz11_d_GamZ1E2 = alpz1001_d_GamZ1*alpz11_d_GamZ1;
	TComplexD deltx1001_deltx11_d_BetX1E2 = deltx1001_d_BetX1*deltx11_d_BetX1;
	TComplexD deltz1001_deltz11_d_BetZ1E2 = deltz1001_d_BetZ1*deltz11_d_BetZ1;
	TComplexD alpx12_d_GamX1 = alpx12*InvGamX1;
	TComplexD alpx11E2_d_GamX1E2 = alpx11_d_GamX1*alpx11_d_GamX1;
	TComplexD alpx1000_alpx12_d_GamX1E2 = alpx1000_d_GamX1*alpx12_d_GamX1;
	TComplexD alpz11E2_d_GamZ1E2 = alpz11_d_GamZ1*alpz11_d_GamZ1;
	TComplexD alpz1000_alpz12_d_GamZ1E2 = alpz1000_d_GamZ1*alpz12_d_GamZ1;
	TComplexD deltx11E2_d_BetX1E2 = deltx11_d_BetX1*deltx11_d_BetX1;
	TComplexD deltx1000_deltx12_d_BetX1E2 = deltx1000_d_BetX1*deltx12_d_BetX1;
	TComplexD deltz11E2_d_BetZ1E2 = deltz11_d_BetZ1*deltz11_d_BetZ1;
	TComplexD deltz1000_deltz12_d_BetZ1E2 = deltz1000_d_BetZ1*deltz12_d_BetZ1;
	TComplexD alpx1010_alpx12_d_GamX1E2 = alpx1010_d_GamX1*alpx12_d_GamX1;
	TComplexD deltx1010_deltx12_d_BetX1E2 = deltx1010_d_BetX1*deltx12_d_BetX1;
	TComplexD alpx1001_alpx12_d_GamX1E2 = alpx1001_d_GamX1*alpx12_d_GamX1;
	TComplexD alpz1001_alpz12_d_GamZ1E2 = alpz1001_d_GamZ1*alpz12_d_GamZ1;
	TComplexD deltx1001_deltx12_d_BetX1E2 = deltx1001_d_BetX1*deltx12_d_BetX1;
	TComplexD deltz1001_deltz12_d_BetZ1E2 = deltz1001_d_BetZ1*deltz12_d_BetZ1;

	//TComplexD Q000 = (ImagI*(dPhi0000 - 0.5*(srTMathFunctions::Argument(BetX1.x, BetX1.y) + srTMathFunctions::Argument(BetZ1.x, BetZ1.y) + srTMathFunctions::Argument(GamX1.x, GamX1.y) + srTMathFunctions::Argument(GamZ1.x, GamZ1.y)))) - (0.25*((alpx1000*alpx1000_d_GamX1) + (alpz1000*alpz1000_d_GamZ1) + (deltx1000*deltx1000_d_BetX1) + (deltz1000*deltz1000_d_BetZ1)));
	TComplexD Q000 = (ImagI*(dPhi0000 - 0.5*(CGenMathFunc::Argument(BetX1.x, BetX1.y) + CGenMathFunc::Argument(BetZ1.x, BetZ1.y) + CGenMathFunc::Argument(GamX1.x, GamX1.y) + CGenMathFunc::Argument(GamZ1.x, GamZ1.y)))) - (0.25*((alpx1000*alpx1000_d_GamX1) + (alpz1000*alpz1000_d_GamZ1) + (deltx1000*deltx1000_d_BetX1) + (deltz1000*deltz1000_d_BetZ1)));
	TComplexD Q010 = (ImagI*dPhi0010) - (0.5*((alpx1000_d_GamX1*alpx1010) + (deltx1000_d_BetX1*deltx1010)));
	TComplexD Q001 = (ImagI*dPhi0001) - (0.5*((alpx1000_d_GamX1*alpx1001) + (alpz1000_d_GamZ1*alpz1001) + (deltx1000_d_BetX1*deltx1001) + (deltz1000_d_BetZ1*deltz1001)));
	TComplexD Q020 = (ImagI*dPhi0020) - (0.25*((alpx1010*alpx1010_d_GamX1) + (deltx1010*deltx1010_d_BetX1)));
	TComplexD Q011 = (-0.5)*((alpx1010_d_GamX1*alpx1001) + (deltx1010_d_BetX1*deltx1001));
	TComplexD Q002 = (ImagI*dPhi0002) - (0.25*((alpx1001*alpx1001_d_GamX1) + (alpz1001*alpz1001_d_GamZ1) + (deltx1001*deltx1001_d_BetX1) + (deltz1001*deltz1001_d_BetZ1)));
	TComplexD Q100 = (ImagI*dPhi0100) - (0.5*((alpx1000_d_GamX1*alpx11) + (alpz1000_d_GamZ1*alpz11) + (deltx1000_d_BetX1*deltx11) + (deltz1000_d_BetZ1*deltz11)));
	TComplexD Q110 = (ImagI*dPhi0110) - (0.5*((alpx1010_d_GamX1*alpx11) + (deltx1010_d_BetX1*deltx11)));
	TComplexD Q101 = (ImagI*dPhi0101) - (0.5*((alpx1001_d_GamX1*alpx11) + (alpz1001_d_GamZ1*alpz11) + (deltx1001_d_BetX1*deltx11) + (deltz1001_d_BetZ1*deltz11)));
	//TComplexD Q200 = (0,1)*d\[CapitalPhi]0200 - (alpx1000dGamX1*\[Alpha]x12)/2. - (alpz1000dGamZ1*\[Alpha]z12)/2. - \[Alpha]x11**2/(4.*\[CapitalGamma]x1) - \[Alpha]z11**2/(4.*\[CapitalGamma]z1) - \[Delta]x11**2/(4.*Bx1) - (deltx1000dBetX1*\[Delta]x12)/2. - \[Delta]z11**2/(4.*Bz1) - (deltz1000dBetZ1*\[Delta]z12)/2.
	//TComplexD Q210 = (0,1)*d\[CapitalPhi]0210 - (alpx1010dGamX1*\[Alpha]x12)/2. - (deltx1010dBetX1*\[Delta]x12)/2.
	//TComplexD Q201 = (0,1)*d\[CapitalPhi]0201 - (alpx1001dGamX1*\[Alpha]x12)/2. - (alpz1001dGamZ1*\[Alpha]z12)/2. - (deltx1001dBetX1*\[Delta]x12)/2. - (deltz1001dBetZ1*\[Delta]z12)/2.
    TComplexD Q100_d_Bgam = Q100*gAuxPar.InvBgam;
	TComplexD Q100E2_d_BgamE2 = Q100_d_Bgam*Q100_d_Bgam;
	TComplexD Q110_d_Bgam = Q110*gAuxPar.InvBgam;
	TComplexD Q100_Q110_d_BgamE2 = Q100_d_Bgam*Q110_d_Bgam;
	TComplexD Q101_d_Bgam = Q101*gAuxPar.InvBgam;
	TComplexD Q100_Q101_d_BgamE2 = Q100_d_Bgam*Q101_d_Bgam;
	TComplexD Q110E2_d_BgamE2 = Q110_d_Bgam*Q110_d_Bgam;
	TComplexD Q101_Q110_d_BgamE2 = Q101_d_Bgam*Q110_d_Bgam;
	TComplexD Q101E2_d_BgamE2 = Q101_d_Bgam*Q101_d_Bgam;

	TComplexD *pOutB = ArrB;
	*(pOutB++) = Q000 + (0.25*(Q100*Q100_d_Bgam)); //B00 //-(0,0.5)*Arg(Bgam)
    *(pOutB++) = Q010 + (0.5*(Q100_d_Bgam*Q110)); //B10
	*(pOutB++) = Q001 + (0.5*(Q100_d_Bgam*Q101)); //B01
	*(pOutB++) = Q020 + (0.25*(Q110*Q110_d_Bgam)); //B20
	*(pOutB++) = Q011 + (0.5*(Q101*Q110_d_Bgam)); //B11
	*(pOutB++) = Q002 + (0.25*(Q101*Q101_d_Bgam)); //B02
    pOutB->x = eObs*eObs/sqrt(sqrt((BetZ1.AbsE2())*(GamZ1.AbsE2())*(BetX1.AbsE2())*(GamX1.AbsE2()))*(gAuxPar.Bgam)); //C00
	pOutB->y = 0;
	
	TComplexD *pOutA = ArrA;

	for(int i=0; i<4; i++) //for 4 Stokes components
	{
        TComplexD dz1 = (dz[i] - (fz[i]*AlpZ1_d_BetZ1)) + (ez[i]*(AlpZ1_d_BetZ1*AlpZ1_d_BetZ1));
		TComplexD fz_d_BetZ1 = fz[i]*InvBetZ1;
		TComplexD Half_i_fz_d_BetZ1 = Half_i*fz_d_BetZ1;

        TComplexD i_ez(0, ez[i]);
		TComplexD i_AlpZ1_ez_d_BetZ1E2 = i_ez*AlpZ1_d_BetZ1E2;
        TComplexD bz1000 = (bz0000[i] - (deltz1000*i_AlpZ1_ez_d_BetZ1E2)) + ((deltz1000*Half_i_fz_d_BetZ1) - (cz0000[i]*AlpZ1_d_BetZ1));
		TComplexD bz1010 = bz0010[i] - (cz0010[i]*AlpZ1_d_BetZ1);
        TComplexD bz1001 = (bz0001[i] - (deltz1001*i_AlpZ1_ez_d_BetZ1E2)) + ((deltz1001*Half_i_fz_d_BetZ1) - (cz0001[i]*AlpZ1_d_BetZ1));
		TComplexD bz11 = (bz01[i] - (deltz11*i_AlpZ1_ez_d_BetZ1E2) + ((deltz11*Half_i_fz_d_BetZ1) - (cz01[i]*AlpZ1_d_BetZ1)));
		TComplexD bz12 = (bz02[i] - (deltz12*i_AlpZ1_ez_d_BetZ1E2) + ((deltz12*Half_i_fz_d_BetZ1) - (cz02[i]*AlpZ1_d_BetZ1)));
		double Two_ez = 2.*ez[i];
		TComplexD m1 = (m[i] + (Two_ez*gAuxPar.ThetXZp)*AlpZ1_d_BetZ1E2) - ((n[i]*AlpZ1_d_BetZ1) + (gAuxPar.ThetXZp*fz_d_BetZ1));
		TComplexD v1 = (v[i] + (Two_ez*gAuxPar.ThetXpZp)*AlpZ1_d_BetZ1E2) - ((w[i]*AlpZ1_d_BetZ1) + (gAuxPar.ThetXpZp*fz_d_BetZ1));
		TComplexD ex1 = ex[i] + (dz1*ThetXpZ1E2_d_GamZ1E2) + (ez[i]*ThetXpZpE2_d_BetZ1E2) - (ThetXpZ1_d_GamZ1*v1) - (ThetXpZp_d_BetZ1*w[i]);
		TComplexD fx1 = fx[i] - (m1*ThetXpZ1_d_GamZ1) - (n[i]*ThetXpZp_d_BetZ1) + (dz1*Two_ThetXpZ1_ThetXZ1_d_GamZ1E2) + (ez[i]*Two_ThetXpZp_ThetXZp_d_BetZ1E2) - (ThetXZ1_d_GamZ1*v1) - (ThetXZp_d_BetZ1*w[i]);
		TComplexD cx1000 = cx0000[i] + (ImagI*(0.5*((alpz1000_d_GamZ1*v1) + (deltz1000_d_BetZ1*w[i])) - ((alpz1000_ThetXpZ1_d_GamZ1E2*dz1) + (deltz1000_ThetXpZp_d_BetZ1E2*ez[i]))) - ((bz1000*ThetXpZ1_d_GamZ1) + (cz0000[i]*ThetXpZp_d_BetZ1)));
		TComplexD cx1010 = cx0010[i] - (bz1010*ThetXpZ1_d_GamZ1) - (cz0010[i]*ThetXpZp_d_BetZ1);
		TComplexD cx1001 = ((cx0001[i] - (bz1001*ThetXpZ1_d_GamZ1)) - (cz0001[i]*ThetXpZp_d_BetZ1)) + (ImagI*(((0.5*((v1*alpz1001_d_GamZ1) + (w[i]*deltz1001_d_BetZ1))) - (dz1*alpz1001_ThetXpZ1_d_GamZ1E2)) - (ez[i]*deltz1001_ThetXpZp_d_BetZ1E2)));
		TComplexD cx11 = cx01[i] - (bz11*ThetXpZ1_d_GamZ1) - (cz01[i]*ThetXpZp_d_BetZ1) + (ImagI*(((0.5*((v1*alpz11_d_GamZ1) + (w[i]*deltz11_d_BetZ1))) - (dz1*alpz11_ThetXpZ1_d_GamZ1E2)) - (ez[i]*deltz11_ThetXpZp_d_BetZ1E2)));
		TComplexD cx12 = cx02[i] - (bz12*ThetXpZ1_d_GamZ1) - (cz02[i]*ThetXpZp_d_BetZ1) + (ImagI*(((0.5*((v1*alpz12_d_GamZ1) + (w[i]*deltz12_d_BetZ1))) - (dz1*alpz12_ThetXpZ1_d_GamZ1E2)) - (ez[i]*deltz12_ThetXpZp_d_BetZ1E2)));
		TComplexD dx1 = dx[i] + (AlpX1E2_d_BetX1E2*ex1) - (AlpX1_d_BetX1*fx1) - (m1*ThetXZ1_d_GamZ1) + (dz1*ThetXZ1E2_d_GamZ1E2) - (n[i]*ThetXZp_d_BetZ1) + (ez[i]*ThetXZpE2_d_BetZ1E2);
		TComplexD bx1000 = (bx0000[i] - (AlpX1_d_BetX1*cx1000)) + (ImagI*((0.5*((deltx1000_d_BetX1*fx1) + (alpz1000_d_GamZ1*m1) + (deltz1000_d_BetZ1*n[i]))) - (alpz1000_ThetXZ1_d_GamZ1E2*dz1) - (deltx1000_AlpX1_d_BetX1E2*ex1) - (deltz1000_ThetXZp_d_BetZ1E2*ez[i]))) - (bz1000*ThetXZ1_d_GamZ1) - (cz0000[i]*ThetXZp_d_BetZ1);
		TComplexD bx1010 = (bx0010[i] - (AlpX1_d_BetX1*cx1010)) + (ImagI*((0.5*(deltx1010_d_BetX1*fx1)) - (deltx1010_AlpX1_d_BetX1E2*ex1))) - (bz1010*ThetXZ1_d_GamZ1) - (cz0010[i]*ThetXZp_d_BetZ1);
		TComplexD bx1001 = (bx0001[i] - (AlpX1_d_BetX1*cx1001)) + (ImagI*((0.5*((deltx1001_d_BetX1*fx1) + (alpz1001_d_GamZ1*m1) + (deltz1001_d_BetZ1*n[i]))) - (alpz1001_ThetXZ1_d_GamZ1E2*dz1) - (deltx1001_AlpX1_d_BetX1E2*ex1) - (deltz1001_ThetXZp_d_BetZ1E2*ez[i]))) - (bz1001*ThetXZ1_d_GamZ1) - (cz0001[i]*ThetXZp_d_BetZ1);
		TComplexD bx11 = (bx01[i] - (AlpX1_d_BetX1*cx11)) + (ImagI*((0.5*((deltx11_d_BetX1*fx1) + (alpz11_d_GamZ1*m1) + (deltz11_d_BetZ1*n[i]))) - (alpz11_ThetXZ1_d_GamZ1E2*dz1) - (deltx11_AlpX1_d_BetX1E2*ex1) - (deltz11_ThetXZp_d_BetZ1E2*ez[i]))) - (bz11*ThetXZ1_d_GamZ1) - (cz01[i]*ThetXZp_d_BetZ1);
		TComplexD bx12 = (bx02[i] - (AlpX1_d_BetX1*cx12)) + (ImagI*((0.5*((deltx12_d_BetX1*fx1) + (alpz12_d_GamZ1*m1) + (deltz12_d_BetZ1*n[i]))) - (alpz12_ThetXZ1_d_GamZ1E2*dz1) - (deltx12_AlpX1_d_BetX1E2*ex1) - (deltz12_ThetXZp_d_BetZ1E2*ez[i]))) - (bz12*ThetXZ1_d_GamZ1) - (cz02[i]*ThetXZp_d_BetZ1);
		TComplexD a1000 = a0000[i] + (Half_i*((alpx1000_d_GamX1*bx1000) + (alpz1000_d_GamZ1*bz1000) + (cx1000*deltx1000_d_BetX1) + (cz0000[i]*deltz1000_d_BetZ1))) - 0.25*((alpx1000E2_d_GamX1E2*dx1) + (alpz1000E2_d_GamZ1E2*dz1) + (deltx1000E2_d_BetX1E2*ex1) + (deltz1000E2_d_BetZ1E2*ez[i])) + (ex1*HalfInvBetX1) + (ez[i]*HalfInvBetZ1) + (dx1*HalfInvGamX1) + (dz1*HalfInvGamZ1);
		TComplexD a1010 = a0010[i] + (Half_i*((alpx1010_d_GamX1*bx1000) + (alpx1000_d_GamX1*bx1010) + (alpz1000_d_GamZ1*bz1010) + (cx1010*deltx1000_d_BetX1) + (cx1000*deltx1010_d_BetX1) + (cz0010[i]*deltz1000_d_BetZ1))) - (0.5*((alpx1000_alpx1010_d_GamX1E2*dx1) + (deltx1000_deltx1010_d_BetX1E2*ex1)));
		TComplexD a1001 = a0001[i] + (Half_i*((alpx1001_d_GamX1*bx1000) + (alpx1000_d_GamX1*bx1001) + (alpz1001_d_GamZ1*bz1000) + (alpz1000_d_GamZ1*bz1001) + (cx1001*deltx1000_d_BetX1) + (cx1000*deltx1001_d_BetX1) + (cz0001[i]*deltz1000_d_BetZ1) + (cz0000[i]*deltz1001_d_BetZ1))) - (0.5*((alpx1000_alpx1001_d_GamX1E2*dx1) + (alpz1000_alpz1001_d_GamZ1E2*dz1) + (deltx1000_deltx1001_d_BetX1E2*ex1) + (deltz1000_deltz1001_d_BetZ1E2*ez[i])));
		TComplexD a1020 = a0020[i] + (Half_i*((alpx1010_d_GamX1*bx1010) + (cx1010*deltx1010_d_BetX1))) - (0.25*((alpx1010E2_d_GamX1E2*dx1) + (deltx1010E2_d_BetX1E2*ex1)));
		TComplexD a1011 = a0011[i] + (Half_i*((alpx1010_d_GamX1*bx1001) + (alpx1001_d_GamX1*bx1010) + (alpz1001_d_GamZ1*bz1010) + (cx1010*deltx1001_d_BetX1) + (cx1001*deltx1010_d_BetX1) + (cz0010[i]*deltz1001_d_BetZ1))) - (0.5*((alpx1001_alpx1010_d_GamX1E2*dx1) + (deltx1001_deltx1010_d_BetX1E2*ex1)));
		TComplexD a1002 = a0002[i] + (Half_i*((alpx1001_d_GamX1*bx1001) + (alpz1001_d_GamZ1*bz1001) + (cx1001*deltx1001_d_BetX1) + (cz0001[i]*deltz1001_d_BetZ1))) - (0.25*((alpx1001E2_d_GamX1E2*dx1) + (alpz1001E2_d_GamZ1E2*dz1) + (deltx1001E2_d_BetX1E2*ex1) + (deltz1001E2_d_BetZ1E2*ez[i])));
		TComplexD a1100 = a0100[i] + (Half_i*((alpx11_d_GamX1*bx1000) + (alpx1000_d_GamX1*bx11) + (alpz11_d_GamZ1*bz1000) + (alpz1000_d_GamZ1*bz11) + (cx11*deltx1000_d_BetX1) + (cx1000*deltx11_d_BetX1) + (cz01[i]*deltz1000_d_BetZ1) + (cz0000[i]*deltz11_d_BetZ1))) - (0.5*((alpx1000_alpx11_d_GamX1E2*dx1) + (alpz1000_alpz11_d_GamZ1E2*dz1) + (deltx1000_deltx11_d_BetX1E2*ex1) + (deltz1000_deltz11_d_BetZ1E2*ez[i])));
		TComplexD a1110 = a0110[i] + (Half_i*((alpx11_d_GamX1*bx1010) + (alpx1010_d_GamX1*bx11) + (alpz11_d_GamZ1*bz1010) + (cx11*deltx1010_d_BetX1) + (cx1010*deltx11_d_BetX1) + (cz0010[i]*deltz11_d_BetZ1))) - (0.5*((alpx1010_alpx11_d_GamX1E2*dx1) + (deltx1010_deltx11_d_BetX1E2*ex1)));
		TComplexD a1101 = a0101[i] + (Half_i*((alpx11_d_GamX1*bx1001) + (alpx1001_d_GamX1*bx11) + (alpz11_d_GamZ1*bz1001) + (alpz1001_d_GamZ1*bz11) + (cx11*deltx1001_d_BetX1) + (cx1001*deltx11_d_BetX1) + (cz01[i]*deltz1001_d_BetZ1) + (cz0001[i]*deltz11_d_BetZ1))) - (0.5*((alpx1001_alpx11_d_GamX1E2*dx1) + (alpz1001_alpz11_d_GamZ1E2*dz1) + (deltx1001_deltx11_d_BetX1E2*ex1) + (deltz1001_deltz11_d_BetZ1E2*ez[i])));
		TComplexD a1200 = a0200[i] + (Half_i*((alpx12_d_GamX1*bx1000) + (alpx11_d_GamX1*bx11) + (alpx1000_d_GamX1*bx12) + (alpz12_d_GamZ1*bz1000) + (alpz11_d_GamZ1*bz11) + (alpz1000_d_GamZ1*bz12) + (cx12*deltx1000_d_BetX1) + (cx11*deltx11_d_BetX1) + (cx1000*deltx12_d_BetX1) + (cz02[i]*deltz1000_d_BetZ1) + (cz01[i]*deltz11_d_BetZ1) + (cz0000[i]*deltz12_d_BetZ1))) - (0.5*((alpx1000_alpx12_d_GamX1E2*dx1) + (alpz1000_alpz12_d_GamZ1E2*dz1) + (deltx1000_deltx12_d_BetX1E2*ex1) + (deltz1000_deltz12_d_BetZ1E2*ez[i]))) - (0.25*((alpx11E2_d_GamX1E2*dx1) + (alpz11E2_d_GamZ1E2*dz1) + (deltx11E2_d_BetX1E2*ex1) + (deltz11E2_d_BetZ1E2*ez[i])));
		TComplexD a1210 = a0210[i] + (Half_i*((alpx12_d_GamX1*bx1010) + (alpx1010_d_GamX1*bx12) + (alpz12_d_GamZ1*bz1010) + (cx12*deltx1010_d_BetX1) + (cx1010*deltx12_d_BetX1) + (cz0010[i]*deltz12_d_BetZ1))) - 0.5*((alpx1010_alpx12_d_GamX1E2*dx1) + (deltx1010_deltx12_d_BetX1E2*ex1));
		TComplexD a1201 = a0201[i] + (Half_i*((alpx12_d_GamX1*bx1001) + (alpx1001_d_GamX1*bx12) + (alpz12_d_GamZ1*bz1001) + (alpz1001_d_GamZ1*bz12) + (cx12*deltx1001_d_BetX1) + (cx1001*deltx12_d_BetX1) + (cz02[i]*deltz1001_d_BetZ1) + (cz0001[i]*deltz12_d_BetZ1))) - (0.5*((alpx1001_alpx12_d_GamX1E2*dx1) + (alpz1001_alpz12_d_GamZ1E2*dz1) + (deltx1001_deltx12_d_BetX1E2*ex1) + (deltz1001_deltz12_d_BetZ1E2*ez[i])));
		
/*A00*/	*(pOutA++) = a1000 + (a1200*gAuxPar.HalfInvBgam) + (0.5*(a1100*Q100_d_Bgam)) + (0.25*(a1200*Q100E2_d_BgamE2));
/*A10*/	*(pOutA++) = a1010 + (a1210*gAuxPar.HalfInvBgam) + (0.5*((a1110*Q100_d_Bgam) + (a1200*Q100_Q110_d_BgamE2) + (a1100*Q110_d_Bgam))) + (0.25*(a1210*Q100E2_d_BgamE2));
/*A01*/	*(pOutA++) = a1001 + (a1201*gAuxPar.HalfInvBgam) + (0.5*((a1101*Q100_d_Bgam) + (a1200*Q100_Q101_d_BgamE2) + (a1100*Q101_d_Bgam))) + (0.25*(a1201*Q100E2_d_BgamE2));
/*A20*/	*(pOutA++) = a1020 + (0.5*((a1210*Q100_Q110_d_BgamE2) + (a1110*Q110_d_Bgam))) + (0.25*(a1200*Q110E2_d_BgamE2));
/*A11*/	*(pOutA++) = a1011 + (0.5*((a1210*Q100_Q101_d_BgamE2) + (a1201*Q100_Q110_d_BgamE2) + (a1110*Q101_d_Bgam) + (a1200*Q101_Q110_d_BgamE2) + (a1101*Q110_d_Bgam)));
/*A02*/	*(pOutA++) = a1002 + (0.5*((a1201*Q100_Q101_d_BgamE2) + (a1101*Q101_d_Bgam))) + (0.25*(a1200*Q101E2_d_BgamE2));
	}
}

//*************************************************************************

//void srTRadIntThickBeam::Integrate1DStokesArr(srTStokes* StokesArr, int Np, double h, srTStokes* pInitDer, srTStokes* pFinDer, srTStokes& ResSt)
void srTRadIntThickBeam::Integrate1DStokesArr(srTStokes* StokesArr, long long Np, double h, srTStokes* pInitDer, srTStokes* pFinDer, srTStokes& ResSt)
{//Np is assumed odd and >= 5!!!
	const double we = 7./15.;
    const double w1 = 16./15.;
    const double w2 = 14./15.;
    const double wd = 1./15.; 

	srTStokes Sum1 = StokesArr[1];
	srTStokes Sum2(0,0,0,0);

	srTStokes *t = StokesArr + 2;
	//int nLoops = (Np - 3) >> 1;
	long long nLoops = (Np - 3) >> 1;
    //for(int i=0; i<nLoops; i++)
    for(long long i=0; i<nLoops; i++)
	{
		Sum2 += *(t++);
		Sum1 += *(t++);
	}
    srTStokes *pFe1 = StokesArr;
    srTStokes *pFe2 = StokesArr + (Np - 1);

	ResSt = h*((we*(*pFe1 + *pFe2)) + (w1*Sum1) + (w2*Sum2));
	double he2_wd = h*h*wd;

    if(pInitDer != 0) ResSt += he2_wd*(*pInitDer);
    if(pFinDer != 0) ResSt += he2_wd*(*pFinDer);
}

//*************************************************************************
    
//void srTRadIntThickBeam::Integrate1DStokesFunc_EvenMesh_OddNp(srTFieldBasedArrays& FldArr, int it, int i_Offset, srTStokes* pFe1, srTStokes& ResSt)
void srTRadIntThickBeam::Integrate1DStokesFunc_EvenMesh_OddNp(srTFieldBasedArrays& FldArr, long long it, long long i_Offset, srTStokes* pFe1, srTStokes& ResSt)
{//Np is assumed odd and >= 5!!!
	const double we = 7./15.;
    const double w1 = 16./15.;
    const double w2 = 14./15.;
    const double wd = 1./15.; 

	//int iStart = it + i_Offset;
	//int Np = FldArr.Ns - iStart;
	long long iStart = it + i_Offset;
	long long Np = FldArr.Ns - iStart;

	srTStokes Fe1, Fe2, F1, F2, Sum1(0,0,0,0), Sum2(0,0,0,0);
    srTStokes dFe1, dFe2;

    srTStokes AuxStokes3[3];
    ComputeStokesAtOneObsPoint_DerivOverS_FromAB(FldArr, iStart, it, 0, dFe1, AuxStokes3);
	Fe1 = AuxStokes3[0];
	Sum1 += AuxStokes3[1];
	Sum2 += AuxStokes3[2];

	ComputeStokesAtOneObsPoint_DerivOverS_FromAB(FldArr, FldArr.Ns - 1,	it,	2, dFe2, AuxStokes3);
	Fe2 = AuxStokes3[2];
	Sum1 += AuxStokes3[1];

	if(Np > 5)
	{
        Sum2 += AuxStokes3[0];
        ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, iStart + 3, it, F1);
		Sum1 += F1;

		//int nLoops = (Np - 7) >> 1;
		//int ip = iStart + 4;
        long long nLoops = (Np - 7) >> 1;
		long long ip = iStart + 4;
        //for(int k=0; k<nLoops; k++)
        for(long long k=0; k<nLoops; k++)
        {
            ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, ip++, it, F2);
            Sum2 += F2;
            ComputeStokesAtOneObsPoint_FuncForInteg2D(FldArr, ip++, it, F1);
            Sum1 += F1;
        }
	}

	double h = FldArr.sStep;
	ResSt = h*((we*(Fe1 + Fe2)) + (w1*Sum1) + (w2*Sum2) + (h*wd)*(dFe1 - dFe2));
}

//*************************************************************************

void srTRadIntThickBeam::ComputeStokesAtOneObsPoint_InternIntens_EvenMesh(srTFieldBasedArrays& FldArr, srTStokes& ResSt)
{
	const double we = 7./15.;
    const double w1 = 16./15.;
    const double w2 = 14./15.;
    const double wd = 1./15.; 

	srTStokes Fe1, Fe2(0,0,0,0), F1, F2, Sum1(0,0,0,0), Sum2(0,0,0,0);

    Integrate1DStokesFunc_EvenMesh(FldArr, 0, Fe1);
    Integrate1DStokesFunc_EvenMesh(FldArr, 1, F1);
	Sum1 += F1;
    Integrate1DStokesFunc_EvenMesh(FldArr, 2, F2);
	Sum2 += F2;
	srTStokes AuxArrStokes1[] = {Fe1, F1, F2};
	srTStokes dFe1 = DerivRight3p(AuxArrStokes1, FldArr.sStep);

    Integrate1DStokesFunc_EvenMesh(FldArr, FldArr.Ns - 2, F1);
	Sum1 += F1;
	if(FldArr.Ns > 5)
	{
        Integrate1DStokesFunc_EvenMesh(FldArr, FldArr.Ns - 3, F2);
        Sum2 += F2;
	}
	srTStokes AuxArrStokes2[] = {F2, F1, Fe2};
    srTStokes dFe2 = DerivLeft3p(AuxArrStokes2, FldArr.sStep);

	if(FldArr.Ns > 5)
	{
        Integrate1DStokesFunc_EvenMesh(FldArr, 3, F1);
        Sum1 += F1;

        //int nLoops = (FldArr.Ns - 7) >> 1;
        long long nLoops = (FldArr.Ns - 7) >> 1;
        int ip = 4;
        for(int k=0; k<nLoops; k++)
        {
			Integrate1DStokesFunc_EvenMesh(FldArr, ip++, F2);
            Sum2 += F2;
            Integrate1DStokesFunc_EvenMesh(FldArr, ip++, F1);
            Sum1 += F1;
        }
	}

	double h = FldArr.sStep;
	ResSt = h*((we*(Fe1 + Fe2)) + (w1*Sum1) + (w2*Sum2) + (h*wd)*(dFe1 - dFe2));
}

//*************************************************************************
