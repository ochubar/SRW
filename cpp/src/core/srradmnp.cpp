/************************************************************************//**
 * File: srradmnp.cpp
 * Description: Various "manipulations" with Radiation data (e.g. "extraction" of Intensity from Electric Field, etc.)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srradmnp.h"
#include "gmfft.h"
#include "gmmeth.h"
//#include "gminterp.h"
#include "gmrand.h"
#include "srtrjdat.h"
#include "srradint.h"
#include "srerror.h"

#ifdef _WITH_OMPH //Pre-processor definition for compiling with OpenMP library
#include "omp.h"
#endif

//DEBUG
//#ifndef __SRIGSEND_H
//#include "srigsend.h"
//#endif
//END DEBUG

//*************************************************************************

void srTRadGenManip::SetupIntCoord(char Cmpn, double Arg, long long& i0, long long& i1, double& InvStepRelArg) //OC26042019
//void srTRadGenManip::SetupIntCoord(char Cmpn, double Arg, long& i0, long& i1, double& InvStepRelArg)
{
	double Step, Start;
	//long N;
	long long N; //OC26042019

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	if(Cmpn == 'e')
	{
		Step = RadAccessData.eStep; Start = RadAccessData.eStart; N = RadAccessData.ne;
	}
	else if(Cmpn == 'x')
	{
		Step = RadAccessData.xStep; Start = RadAccessData.xStart; N = RadAccessData.nx;
	}
	else
	{
		Step = RadAccessData.zStep; Start = RadAccessData.zStart; N = RadAccessData.nz;
	}

	//if(N == 1)
	if(N <= 1) //OC191215
	{
		i0 = i1 = 0; InvStepRelArg = 0.; return;
	}

	double InvStep = 1./Step;
	i0 = long((Arg - Start)*InvStep);
	if(i0 < 0) 
	{
		i0 = 0; i1 = i0;
	}
	else if(i0 >= N - 1) 
	{
		i0 = N - 1; i1 = i0;
	}
	else
	{
		i1 = i0 + 1;
	}

	InvStepRelArg = InvStep*(Arg - i0*Step - Start);
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecFlux1DvsE(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
	float *pI = RadExtract.pExtractedData;

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//long Nz_mi_1 = RadAccessData.nz - 1;
	//long Nx_mi_1 = RadAccessData.nx - 1;
	long long Nz_mi_1 = RadAccessData.nz - 1; //OC26042019
	long long Nx_mi_1 = RadAccessData.nx - 1;

    //long iePerE = 0;
    long long iePerE = 0;
	for(long long ie=0; ie<RadAccessData.ne; ie++) //OC26042019
	//for(long ie=0; ie<RadAccessData.ne; ie++)
	{
        //long izPerZ = 0;
        long long izPerZ = 0;

		double Sum = 0;
		for(long long iz=0; iz<RadAccessData.nz; iz++) //OC26042019
		//for(long iz=0; iz<RadAccessData.nz; iz++)
		{
			float *pEx_StartForX = pEx0 + izPerZ;
            float *pEz_StartForX = pEz0 + izPerZ;
            //long ixPerX = 0;
            long long ixPerX = 0;

			double SumX = 0.;

            for(long long ix=0; ix<RadAccessData.nx; ix++) //OC26042019
            //for(long ix=0; ix<RadAccessData.nx; ix++)
            {
				float* pEx = pEx_StartForX + (ixPerX + iePerE);
				float* pEz = pEz_StartForX + (ixPerX + iePerE);

				double CurI = IntensityComponent(pEx, pEz, PolCom, 0);
                if((ix == 0) || (ix == Nx_mi_1)) CurI *= 0.5;
				SumX += CurI;

                ixPerX += PerX; 
			}

            if((iz == 0) || (iz == Nz_mi_1)) SumX *= 0.5;
			Sum += SumX;

			izPerZ += PerZ; 
		}
		*(pI++) = (float)(Sum*(RadAccessData.xStep)*(RadAccessData.zStep)*(1E+06)); // to ensure Photons/s/0.1%bw
		iePerE += 2; 
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractMultiElecFlux1DvsE(srTRadExtract& RadExtract)
{
	int result;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	if((RadAccessData.nx == 1) || (RadAccessData.nz == 1)) return NEEDS_MORE_THAN_ONE_HOR_AND_VERT_OBS_POINT;

	long Nx = RadAccessData.nx;
	long Nz = RadAccessData.nz;
	//long long Nx = RadAccessData.nx; //OC26042019
	//long long Nz = RadAccessData.nz;

	CGenMathFFT2D FFT2D;
	FFT2D.NextCorrectNumberForFFT(Nx);
	FFT2D.NextCorrectNumberForFFT(Nz);

	//long TotAmOfNewData = (Nx*Nz) << 1;
	long long TotAmOfNewData = (((long long)Nx)*((long long)Nz)) << 1;

	float* pTempStorage = new float[TotAmOfNewData];
	if(pTempStorage == 0) return MEMORY_ALLOCATION_FAILURE;

		//test
		float *tTempStorage = pTempStorage;
		//for(long i=0; i<TotAmOfNewData; i++) *(tTempStorage++) = 0;
		for(long long i=0; i<TotAmOfNewData; i++) *(tTempStorage++) = 0;
		//end test

	srTRadExtract OwnRadExtract = RadExtract;
	OwnRadExtract.pExtractedData = pTempStorage;
	OwnRadExtract.PlotType = 3; // vs x&z
	OwnRadExtract.Int_or_Phase = 0;
	OwnRadExtract.ePh = RadAccessData.eStart;

	//long Ne = RadAccessData.ne;
	long long Ne = RadAccessData.ne; //OC26042019
	float *pF = RadExtract.pExtractedData;

	//long Nz_mi_1 = Nz - 1;
	//long Nx_mi_1 = Nx - 1;
	long long Nz_mi_1 = Nz - 1;
	long long Nx_mi_1 = Nx - 1;

	//double xStepLoc = ((RadAccessData.xStep)*(RadAccessData.nx - 1))/Nx_mi_1;
	//double zStepLoc = ((RadAccessData.zStep)*(RadAccessData.nz - 1))/Nz_mi_1;

	for(long ie=0; ie<Ne; ie++)
	{
		if(result = ExtractSingleElecIntensity2DvsXZ(OwnRadExtract)) return result;
		float *pLocI = OwnRadExtract.pExtractedData;
		if(result = ConvoluteWithElecBeamOverTransvCoord(pLocI, Nx, Nz)) return result;

		double Sum = 0;
		for(long iz=0; iz<Nz; iz++)
		{
			double SumX = 0.;
			for(long ix=0; ix<Nx; ix++)
			{
				double CurI = *pLocI;
                if((ix == 0) || (ix == Nx_mi_1)) CurI *= 0.5;
				SumX += CurI;
				pLocI += 2;
			}
            if((iz == 0) || (iz == Nz_mi_1)) SumX *= 0.5;
			Sum += SumX;
		}

		*(pF++) = (float)(Sum*(RadAccessData.xStep)*(RadAccessData.zStep)*(1E+06)); // to ensure Photons/s/0.1%bw
		OwnRadExtract.ePh += RadAccessData.eStep;
	}

	if(pTempStorage != 0) delete[] pTempStorage;
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity1DvsE(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC16042020

	float *pI = 0, *pI1 = 0, *pI2 = 0, *pI3 = 0; //OC16042020
	double *pId = 0, *pI1d = 0, *pI2d = 0, *pI3d = 0;
	long ne = RadAccessData.ne;
	//float *pI = 0;
	//DOUBLE *pId = 0;
	//double *pId = 0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	if(Int_or_ReE != 2)
	{
		pI = RadExtract.pExtractedData;
		if(allStokesReq) //OC16042020
		{
			pI1 = pI + ne; pI2 = pI1 + ne; pI3 = pI2 + ne;
		}
	}
	else
	{
		pId = RadExtract.pExtractedDataD;
		if(allStokesReq) //OC16042020
		{
			pI1d = pId + ne; pI2d = pI1d + ne; pI3d = pI2d + ne;
		}
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//long ix0=0, ix1=0, iz0=0, iz1=0;
	long long ix0=0, ix1=0, iz0=0, iz1=0; //OC26042019
	double InvStepRelArg1, InvStepRelArg2;
	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg1);
	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);

	//long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;

	//long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	float *pEx_StartForE_ix0_iz0 = pEx_StartForX_iz0 + ix0PerX, *pEx_StartForE_ix1_iz0 = pEx_StartForX_iz0 + ix1PerX;
	float *pEx_StartForE_ix0_iz1 = pEx_StartForX_iz1 + ix0PerX, *pEx_StartForE_ix1_iz1 = pEx_StartForX_iz1 + ix1PerX;
	float *pEz_StartForE_ix0_iz0 = pEz_StartForX_iz0 + ix0PerX, *pEz_StartForE_ix1_iz0 = pEz_StartForX_iz0 + ix1PerX;
	float *pEz_StartForE_ix0_iz1 = pEz_StartForX_iz1 + ix0PerX, *pEz_StartForE_ix1_iz1 = pEz_StartForX_iz1 + ix1PerX;
	//long iePerE = 0;
	long long iePerE = 0;

	double BufI, BufI1, BufI2, BufI3; //OC16042020
	for(long ie=0; ie<ne; ie++)
	//for(long ie=0; ie<RadAccessData.ne; ie++)
	{
		float *ExPtrs[] = { (pEx_StartForE_ix0_iz0 + iePerE), (pEx_StartForE_ix1_iz0 + iePerE),
							(pEx_StartForE_ix0_iz1 + iePerE), (pEx_StartForE_ix1_iz1 + iePerE)};
		float *EzPtrs[] = { (pEz_StartForE_ix0_iz0 + iePerE), (pEz_StartForE_ix1_iz0 + iePerE),
							(pEz_StartForE_ix0_iz1 + iePerE), (pEz_StartForE_ix1_iz1 + iePerE)};

		if(!allStokesReq) //OC16042020
		{
			BufI = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
			if(RadExtract.pExtractedData != 0) *(pI++) = (float)BufI;
			if(RadExtract.pExtractedDataD != 0) *(pId++) = (double)BufI; //OC17042020
			//if(RadExtract.pExtractedDataD != 0) *(pId++) = (float)BufI;
		}
		else //OC16042020
		{
			BufI = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -1, Int_or_ReE);
			BufI1 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -2, Int_or_ReE);
			BufI2 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -3, Int_or_ReE);
			BufI3 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -4, Int_or_ReE);

			if(RadExtract.pExtractedData != 0)
			{
				*(pI++) = (float)BufI; *(pI1++) = (float)BufI1; *(pI2++) = (float)BufI2; *(pI3++) = (float)BufI3;
			}
			if(RadExtract.pExtractedDataD != 0)
			{
				*(pId++) = BufI; *(pI1d++) = BufI1; *(pI2d++) = BufI2; *(pI3d++) = BufI3;
			}
		}

		iePerE += 2; 
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity1DvsX(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	if(Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC17042020

	float *pI = 0, *pI1 = 0, *pI2 = 0, *pI3 = 0; //OC17042020
	double *pId = 0, *pI1d = 0, *pI2d = 0, *pI3d = 0;
	long ne = RadAccessData.ne, nx = RadAccessData.nx;
	//float *pI = 0;
	//DOUBLE *pId = 0;
	//double *pId = 0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	if(Int_or_ReE != 2)
	{
		pI = RadExtract.pExtractedData;
		if(allStokesReq) //OC17042020
		{
			pI1 = pI + nx; pI2 = pI1 + nx; pI3 = pI2 + nx;
		}
	}
	else
	{
		pId = RadExtract.pExtractedDataD;
		if(allStokesReq) //OC17042020
		{
			pI1d = pId + nx; pI2d = pI1d + nx; pI3d = pI2d + nx;
		}
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//long ie0=0, ie1=0, iz0=0, iz1=0;
	long long ie0=0, ie1=0, iz0=0, iz1=0; //OC26042019
	double InvStepRelArg1, InvStepRelArg2;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC140813
	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);

	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC140813
	double *arAuxInt = 0, resInt, resInt1, resInt2, resInt3; //OC17042020
	//double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC140813

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	//long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019

	//long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;
	//long ixPerX = 0;
	long long ixPerX = 0;
	long ie; //OC17042020

	for(long ix=0; ix<nx; ix++) //OC17042020
	//for(long ix=0; ix<RadAccessData.nx; ix++)
	{
		float *pEx_StartForE_iz0 = pEx_StartForX_iz0 + ixPerX, *pEx_StartForE_iz1 = pEx_StartForX_iz1 + ixPerX;
		float *pEz_StartForE_iz0 = pEz_StartForX_iz0 + ixPerX, *pEz_StartForE_iz1 = pEz_StartForX_iz1 + ixPerX;

		float *ExPtrs[] = { (pEx_StartForE_iz0 + Two_ie0), (pEx_StartForE_iz0 + Two_ie1),
							(pEx_StartForE_iz1 + Two_ie0), (pEx_StartForE_iz1 + Two_ie1)};
		float *EzPtrs[] = { (pEz_StartForE_iz0 + Two_ie0), (pEz_StartForE_iz0 + Two_ie1),
							(pEz_StartForE_iz1 + Two_ie0), (pEz_StartForE_iz1 + Two_ie1)};
		//OC150813
		//if(pI != 0) *(pI++) =   IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
		//if(pId != 0) *(pId++) = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);

		if(intOverEnIsRequired) //OC150813
		{//integrate over photon energy / time
			double *tInt = arAuxInt; 
			float *pEx_StAux = pEx_StartForE_iz0;
			float *pEx_FiAux = pEx_StartForE_iz1;
			float *pEz_StAux = pEz_StartForE_iz0;
			float *pEz_FiAux = pEz_StartForE_iz1;

			if(!allStokesReq) //OC17042020
			{
				for(ie=0; ie<ne; ie++) //OC17042020
				//for(int ie=0; ie<RadAccessData.ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, PolCom, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2;
					pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
			}
			else //OC17042020
			{
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -1, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_iz0; pEx_FiAux = pEx_StartForE_iz1; pEz_StAux = pEz_StartForE_iz0; pEz_FiAux = pEz_StartForE_iz1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -2, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt1 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_iz0; pEx_FiAux = pEx_StartForE_iz1; pEz_StAux = pEz_StartForE_iz0; pEz_FiAux = pEz_StartForE_iz1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -3, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt2 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_iz0; pEx_FiAux = pEx_StartForE_iz1; pEz_StAux = pEz_StartForE_iz0; pEz_FiAux = pEz_StartForE_iz1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -4, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt3 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
			}
		}
		else
		{
			if(!allStokesReq) //OC17042020
			{
				resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
			}
			else //OC17042020
			{
				resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -1, Int_or_ReE);
				resInt1 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -2, Int_or_ReE);
				resInt2 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -3, Int_or_ReE);
				resInt3 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -4, Int_or_ReE);
			}
		}
		//OC150813
		if(pI != 0) *(pI++) = (float)resInt;
		if(pId != 0) *(pId++) = resInt; //OC17042020
		//if(pId != 0) *(pId++) = (double)resInt;
		if(allStokesReq) //OC17042020
		{
			if(RadExtract.pExtractedData != 0)
			{
				*(pI1++) = (float)resInt1; *(pI2++) = (float)resInt2; *(pI3++) = (float)resInt3;
			}
			else
			{
				*(pI1d++) = resInt1; *(pI2d++) = resInt2; *(pI3d++) = resInt3;
			}
		}
		ixPerX += PerX;
	}
	if(arAuxInt != 0) delete[] arAuxInt; //OC150813
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity1DvsZ(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	if(Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC17042020

	float *pI = 0, *pI1 = 0, *pI2 = 0, *pI3 = 0; //OC17042020
	double *pId = 0, *pI1d = 0, *pI2d = 0, *pI3d = 0;
	long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz;
	//float *pI = 0;
	//DOUBLE *pId = 0;
	//double *pId = 0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	if(Int_or_ReE != 2)
	{
		pI = RadExtract.pExtractedData;
		if(allStokesReq) //OC17042020
		{
			pI1 = pI + nz; pI2 = pI1 + nz; pI3 = pI2 + nz;
		}
	}
	else
	{
		pId = RadExtract.pExtractedDataD;
		if(allStokesReq) //OC17042020
		{
			pI1d = pId + nz; pI2d = pI1d + nz; pI3d = pI2d + nz;
		}
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long long PerX = RadAccessData.ne << 1;
	//long long PerZ = PerX*RadAccessData.nx;
	long long PerX = ((long long)ne) << 1; //OC17042020
	long long PerZ = PerX*nx;

	//long ie0=0, ie1=0, ix0=0, ix1=0;
	long long ie0=0, ie1=0, ix0=0, ix1=0; //OC26042019
	double InvStepRelArg1, InvStepRelArg2;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC150813
	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg2);

	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC150813
	double *arAuxInt = 0, resInt, resInt1, resInt2, resInt3; //OC17042020
	//double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC150813

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	//long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
	//long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	//long izPerZ = 0;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	long long izPerZ = 0;
	long ie; //OC17042020

	for(long long iz=0; iz<nz; iz++) //OC17042020
	//for(long long iz=0; iz<RadAccessData.nz; iz++) //OC26042019
	//for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		float *pEx_StartForE_ix0 = pEx_StartForX + ix0PerX, *pEx_StartForE_ix1 = pEx_StartForX + ix1PerX;
		float *pEz_StartForE_ix0 = pEz_StartForX + ix0PerX, *pEz_StartForE_ix1 = pEz_StartForX + ix1PerX;

		float *ExPtrs[] = { (pEx_StartForE_ix0 + Two_ie0), (pEx_StartForE_ix0 + Two_ie1),
							(pEx_StartForE_ix1 + Two_ie0), (pEx_StartForE_ix1 + Two_ie1)};
		float *EzPtrs[] = { (pEz_StartForE_ix0 + Two_ie0), (pEz_StartForE_ix0 + Two_ie1),
							(pEz_StartForE_ix1 + Two_ie0), (pEz_StartForE_ix1 + Two_ie1)};
		//OC150813
		//if(pI != 0) *(pI++) = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
		//if(pId != 0) *(pId++) = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);

		if(intOverEnIsRequired) //OC150813
		{//integrate over photon energy / time
			double *tInt = arAuxInt; 
			float *pEx_StAux = pEx_StartForE_ix0;
			float *pEx_FiAux = pEx_StartForE_ix1;
			float *pEz_StAux = pEz_StartForE_ix0;
			float *pEz_FiAux = pEz_StartForE_ix1;

			if(!allStokesReq) //OC17042020
			{
				for(ie=0; ie<ne; ie++) //OC17042020
				//for(int ie=0; ie<RadAccessData.ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, PolCom, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2;
					pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep); //OC17042020
				//resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
			}
			else //OC17042020
			{
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -1, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_ix0; pEx_FiAux = pEx_StartForE_ix1; pEz_StAux = pEz_StartForE_ix0; pEz_FiAux = pEz_StartForE_ix1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -2, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt1 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_ix0; pEx_FiAux = pEx_StartForE_ix1; pEz_StAux = pEz_StartForE_ix0; pEz_FiAux = pEz_StartForE_ix1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -3, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt2 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_StartForE_ix0; pEx_FiAux = pEx_StartForE_ix1; pEz_StAux = pEz_StartForE_ix0; pEz_FiAux = pEz_StartForE_ix1;
				for(ie=0; ie<ne; ie++)
				{
					*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, -4, Int_or_ReE);
					pEx_StAux += 2; pEx_FiAux += 2; pEz_StAux += 2; pEz_FiAux += 2;
				}
				resInt3 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);
			}
		}
		else
		{
			if(!allStokesReq) //OC17042020
			{
				resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
			}
			else //OC17042020
			{
				resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -1, Int_or_ReE);
				resInt1 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -2, Int_or_ReE);
				resInt2 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -3, Int_or_ReE);
				resInt3 = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, -4, Int_or_ReE);
			}
		}
		//OC150813
		if(pI != 0) *(pI++) = (float)resInt;
		if(pId != 0) *(pId++) = resInt; //OC17042020
		//if(pId != 0) *(pId++) = (double)resInt;
		if(allStokesReq) //OC17042020
		{
			if(RadExtract.pExtractedData != 0)
			{
				*(pI1++) = (float)resInt1; *(pI2++) = (float)resInt2; *(pI3++) = (float)resInt3;
			}
			else
			{
				*(pI1d++) = resInt1; *(pI2d++) = resInt2; *(pI3d++) = resInt3;
			}
		}
		izPerZ += PerZ;
	}
	if(arAuxInt != 0) delete[] arAuxInt; //OC150813
	return 0;
}

//*************************************************************************

//int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ(srTRadExtract& RadExtract)
//int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ(srTRadExtract& RadExtract, gpuUsageArg *pGpuUsage)
int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ(srTRadExtract& RadExtract, void* pvGPU) //HG02122023
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	if(Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC18042020

	float *pI = 0, *pI1 = 0, *pI2 = 0, *pI3 = 0; //OC17042020
	double *pId = 0, *pI1d = 0, *pI2d = 0, *pI3d = 0;
	long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz;
	//long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz, nwfr = RadAccessData.nwfr;
	//float *pI = 0;
	//DOUBLE *pId = 0;
	//double *pId = 0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	long long nxnz = ((long long)nx)*((long long)nz);
	if(Int_or_ReE != 2)
	{
		pI = RadExtract.pExtractedData;
		if(allStokesReq) //OC17042020
		{
			pI1 = pI + nxnz; pI2 = pI1 + nxnz; pI3 = pI2 + nxnz;
		}
	}
	else
	{
		pId = RadExtract.pExtractedDataD;
		if(allStokesReq) //OC17042020
		{
			pI1d = pId + nxnz; pI2d = pI1d + nxnz; pI3d = pI2d + nxnz;
		}
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long long PerX = RadAccessData.ne << 1;
	//long long PerZ = PerX*RadAccessData.nx;
	long long PerX = ((long long)ne) << 1; //OC18042020
	long long PerZ = PerX*nx;
	long long PerWfr = PerZ*nz;

	//long ie0=0, ie1=0;
	long long ie0=0, ie1=0; //OC26042019
	double InvStepRelArg=0;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg); //OC140813
	//bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC140813
	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (ne > 1); //OC18042020
	double *arAuxInt = 0, resInt, resInt1, resInt2, resInt3; //OC18042020
	//double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg); //OC140813

	double iter = 0, inv_iter_p_1 = 1.; //OC08052021
	double *pMeth = RadExtract.pMeth;
	if(pMeth != 0)
	{
		if(*pMeth == 1)
		{
			iter = *(pMeth + 1);
			inv_iter_p_1 = 1./(iter + 1.);
		}
		else if(*pMeth == 2) iter = -1;
	}

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	//long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
	//long izPerZ = 0;
	long ix, ie;

	//GPU_COND(pGpuUsage,
	//{
	//	ExtractSingleElecIntensity2DvsXZ_GPU(RadExtract, arAuxInt, ie0, ie1, InvStepRelArg, pGpuUsage);
	//})
	//else
	//{
		//long long iwfrPerWfr = 0;
		//for(long long iwfr=0; iwfr<nwfr; iwfr++)
		//{
#ifdef _OFFLOAD_GPU	//HG30112023 //HG02122023
	TGPUUsageArg parGPU(pvGPU); //OC19022024
	if(CAuxGPU::GPUEnabled(&parGPU)) //OC19022024
	//if(CAuxGPU::GPUEnabled((TGPUUsageArg*)pvGPU))
	{
		ExtractSingleElecIntensity2DvsXZ_GPU(RadExtract, arAuxInt, ie0, ie1, InvStepRelArg, &parGPU); //OC19022024
		//ExtractSingleElecIntensity2DvsXZ_GPU(RadExtract, arAuxInt, ie0, ie1, InvStepRelArg, (TGPUUsageArg*)pvGPU);
	}
	else
#endif
	{
		long long izPerZ = 0;
		for(long long iz=0; iz<nz; iz++) //OC18042020
		//for(long long iz=0; iz<RadAccessData.nz; iz++) //OC26042019
		//for(long iz=0; iz<RadAccessData.nz; iz++)
		{
			float *pEx_StartForX = pEx0 + izPerZ;
			float *pEz_StartForX = pEz0 + izPerZ;
			//long ixPerX = 0;

			float *pEx_St = pEx_StartForX + Two_ie0;
			float *pEz_St = pEz_StartForX + Two_ie0;
			float *pEx_Fi = pEx_StartForX + Two_ie1;
			float *pEz_Fi = pEz_StartForX + Two_ie1;

			for(ix=0; ix<nx; ix++) //OC18042020
				//for(long ix=0; ix<RadAccessData.nx; ix++)
			{
				//float *pEx_StartForE = pEx_StartForX + ixPerX;
				//float *pEz_StartForE = pEz_StartForX + ixPerX;
				//float *pEx_St = pEx_StartForE + Two_ie0, *pEx_Fi = pEx_StartForE + Two_ie1;
				//float *pEz_St = pEz_StartForE + Two_ie0, *pEz_Fi = pEz_StartForE + Two_ie1;

				//OC140813
				//if(pI != 0) *(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
				//if(pId != 0) *(pId++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);

				if(intOverEnIsRequired) //OC140813
				{//integrate over photon energy / time
					double *tInt = arAuxInt;
					float *pEx_StAux = pEx_St;
					float *pEz_StAux = pEz_St;

					if(!allStokesReq) //OC17042020
					{
						for(ie=0; ie<ne; ie++) //OC18042020
							//for(int ie=0; ie<RadAccessData.ne; ie++)
						{
							*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, PolCom, Int_or_ReE);
							pEx_StAux += 2;
							pEz_StAux += 2;
						}
						resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep); //OC18042020
						//resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
					}
					else
					{
						for(ie=0; ie<ne; ie++)
						{
							*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, -1, Int_or_ReE);
							pEx_StAux += 2; pEz_StAux += 2;
						}
						resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

						tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
						for(ie=0; ie<ne; ie++)
						{
							*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, -2, Int_or_ReE);
							pEx_StAux += 2; pEz_StAux += 2;
						}
						resInt1 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

						tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
						for(ie=0; ie<ne; ie++)
						{
							*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, -3, Int_or_ReE);
							pEx_StAux += 2; pEz_StAux += 2;
						}
						resInt2 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

						tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
						for(ie=0; ie<ne; ie++)
						{
							*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, -4, Int_or_ReE);
							pEx_StAux += 2; pEz_StAux += 2;
						}
						resInt3 = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);
					}
				}
				else
				{
					if(!allStokesReq) //OC18042020
					{
						resInt = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
					}
					else //OC18042020
					{
						resInt = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -1, Int_or_ReE);
						resInt1 = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -2, Int_or_ReE);
						resInt2 = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -3, Int_or_ReE);
						resInt3 = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -4, Int_or_ReE);
					}
				}

				if(iter == 0) //OC08052021
				{
					//OC140813
					if(pI != 0) *(pI++) = (float)resInt;
					if(pId != 0) *(pId++) = resInt; //OC18042020
					//if(pId != 0) *(pId++) = (double)resInt;
					if(allStokesReq) //OC18042020
					{
						if(RadExtract.pExtractedData != 0)
						{
							*(pI1++) = (float)resInt1; *(pI2++) = (float)resInt2; *(pI3++) = (float)resInt3;
						}
						else
						{
							*(pI1d++) = resInt1; *(pI2d++) = resInt2; *(pI3d++) = resInt3;
						}
					}
				}
				else if(iter > 0) //OC08052021
				{
					if(pI != 0)
					{
						float newI = (float)(((*pI)*iter + resInt)*inv_iter_p_1);
						*(pI++) = newI;
					}
					if(pId != 0)
					{
						double newI = ((*pId)*iter + resInt)*inv_iter_p_1;
						*(pId++) = newI;
					}
					if(allStokesReq)
					{
						if(RadExtract.pExtractedData != 0)
						{
							float newI1 = (float)(((*pI1)*iter + resInt1)*inv_iter_p_1);
							float newI2 = (float)(((*pI2)*iter + resInt2)*inv_iter_p_1);
							float newI3 = (float)(((*pI3)*iter + resInt3)*inv_iter_p_1);
							*(pI1++) = newI1; *(pI2++) = newI2; *(pI3++) = newI3;
						}
						else
						{
							double newI1 = ((*pI1d)*iter + resInt1)*inv_iter_p_1;
							double newI2 = ((*pI2d)*iter + resInt2)*inv_iter_p_1;
							double newI3 = ((*pI3d)*iter + resInt3)*inv_iter_p_1;
							*(pI1d++) = newI1; *(pI2d++) = newI2; *(pI3d++) = newI3;
						}
					}
				}
				else //OC08052021
				{
					if(pI != 0) *(pI++) += (float)resInt;
					if(pId != 0) *(pId++) += resInt;
					if(allStokesReq)
					{
						if(RadExtract.pExtractedData != 0)
						{
							*(pI1++) += (float)resInt1; *(pI2++) += (float)resInt2; *(pI3++) += (float)resInt3;
						}
						else
						{
							*(pI1d++) += resInt1; *(pI2d++) += resInt2; *(pI3d++) += resInt3;
						}
					}
				}

				//ixPerX += PerX;
				pEx_St += PerX;
				pEz_St += PerX;
				pEx_Fi += PerX;
				pEz_Fi += PerX;
			}
			izPerZ += PerZ;
		}
		//iwfrPerWfr += PerWfr;
	}
	if(arAuxInt != 0) delete[] arAuxInt; //OC150813
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity2DvsEX(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC18042020

	float *pI = RadExtract.pExtractedData;
	float *pI1 = 0, *pI2 = 0, *pI3 = 0, *pI4 = 0; //OC18042020
	long ne = RadAccessData.ne, nx = RadAccessData.nx;
	if(allStokesReq)
	{
		long long nenx = ((long long)ne)*((long long)nx);
		pI1 = pI + nenx; pI2 = pI1 + nenx; pI3 = pI2 + nenx;
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long long PerX = RadAccessData.ne << 1;
	//long long PerZ = PerX*RadAccessData.nx;
	long long PerX = ((long long)ne) << 1; //OC18042020
	long long PerZ = PerX*nx;

	//long iz0=0, iz1=0;
	long long iz0=0, iz1=0; //OC26042019
	double InvStepRelArg;
	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg);

	//long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;

	//long ixPerX = 0;
	long long ixPerX = 0;
	long ie; //OC18042020

	for(long ix=0; ix<nx; ix++) //OC18042020
	//for(long ix=0; ix<RadAccessData.nx; ix++)
	{
		float *pEx_StartForE_iz0 = pEx_StartForX_iz0 + ixPerX, *pEx_StartForE_iz1 = pEx_StartForX_iz1 + ixPerX;
		float *pEz_StartForE_iz0 = pEz_StartForX_iz0 + ixPerX, *pEz_StartForE_iz1 = pEz_StartForX_iz1 + ixPerX;
		long iePerE = 0;

		for(ie=0; ie<ne; ie++) //OC18042020
		//for(long ie=0; ie<RadAccessData.ne; ie++)
		{
			float *pEx_St = pEx_StartForE_iz0 + iePerE, *pEx_Fi = pEx_StartForE_iz1 + iePerE;
			float *pEz_St = pEz_StartForE_iz0 + iePerE, *pEz_Fi = pEz_StartForE_iz1 + iePerE;

			if(!allStokesReq) //OC18042020
			{
				*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
			}
			else //OC18042020
			{
				*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -1, Int_or_ReE);
				*(pI1++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -2, Int_or_ReE);
				*(pI2++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -3, Int_or_ReE);
				*(pI3++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -4, Int_or_ReE);
			}
			iePerE += 2;
		}
		ixPerX += PerX;
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity2DvsEZ(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC18042020

	float *pI = RadExtract.pExtractedData;
	float *pI1 = 0, *pI2 = 0, *pI3 = 0, *pI4 = 0; //OC18042020
	long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz;
	if(allStokesReq)
	{
		long long nenz = ((long long)ne)*((long long)nz);
		pI1 = pI + nenz; pI2 = pI1 + nenz; pI3 = pI2 + nenz;
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long long PerX = RadAccessData.ne << 1;
	//long long PerZ = PerX*RadAccessData.nx;
	long long PerX = ((long long)ne) << 1; //OC18042020
	long long PerZ = PerX*nx;

	//long ix0=0, ix1=0;
	long long ix0=0, ix1=0; //OC26042019
	double InvStepRelArg;
	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg);

	//long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	//long izPerZ = 0;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	long long izPerZ = 0;
	long ie;

	for(long iz=0; iz<nz; iz++) //OC18042020
	//for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;

		float *pEx_StartForE_ix0 = pEx_StartForX + ix0PerX, *pEx_StartForE_ix1 = pEx_StartForX + ix1PerX;
		float *pEz_StartForE_ix0 = pEz_StartForX + ix0PerX, *pEz_StartForE_ix1 = pEz_StartForX + ix1PerX;
		long iePerE = 0;

		for(ie=0; ie<ne; ie++) //OC18042020
		//for(long ie=0; ie<RadAccessData.ne; ie++)
		{
			float *pEx_St = pEx_StartForE_ix0 + iePerE, *pEx_Fi = pEx_StartForE_ix1 + iePerE;
			float *pEz_St = pEz_StartForE_ix0 + iePerE, *pEz_Fi = pEz_StartForE_ix1 + iePerE;

			if(!allStokesReq) //OC18042020
			{
				*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
			}
			else //OC18042020
			{
				*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -1, Int_or_ReE);
				*(pI1++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -2, Int_or_ReE);
				*(pI2++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -3, Int_or_ReE);
				*(pI3++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -4, Int_or_ReE);
			}
			iePerE += 2;
		}
		izPerZ += PerZ;
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity3D(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	bool allStokesReq = (PolCom == -5); //OC17042020

	float *pI = RadExtract.pExtractedData;
	float *pI1 = 0, *pI2 = 0, *pI3 = 0, *pI4 = 0; //OC17042020
	long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz;
	if(allStokesReq)
	{
		long long nTot = ((long long)ne)*((long long)nx)*((long long)nz);
		pI1 = pI + nTot; pI2 = pI1 + nTot; pI3 = pI2 + nTot;
	}

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long izPerZ = 0;
	//long long PerX = RadAccessData.ne << 1;
	//long long PerZ = PerX*RadAccessData.nx;
	long long PerX = ((long long)ne) << 1; //OC17042020
	long long PerZ = PerX*nx;
	long long izPerZ = 0;

	for(long iz=0; iz<nz; iz++) //OC17042020
	//for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(long ix=0; ix<nx; ix++) //OC17042020
		//for(long ix=0; ix<RadAccessData.nx; ix++)
		{
			float *pEx_StartForE = pEx_StartForX + ixPerX;
			float *pEz_StartForE = pEz_StartForX + ixPerX;
			long iePerE = 0;

			for(long ie=0; ie<ne; ie++) //OC17042020
			//for(long ie=0; ie<RadAccessData.ne; ie++)
			{
				float* pEx = pEx_StartForE + iePerE;
				float* pEz = pEz_StartForE + iePerE;

				if(!allStokesReq) //OC16042020
				{
					*(pI++) = IntensityComponent(pEx, pEz, PolCom, Int_or_ReE);
				}
				else //OC16042020
				{
					*(pI++) = IntensityComponent(pEx, pEz, -1, Int_or_ReE);
					*(pI1++) = IntensityComponent(pEx, pEz, -2, Int_or_ReE);
					*(pI2++) = IntensityComponent(pEx, pEz, -3, Int_or_ReE);
					*(pI3++) = IntensityComponent(pEx, pEz, -3, Int_or_ReE);
				}
				iePerE += 2; 
			}
			ixPerX += PerX; 
		}
		izPerZ += PerZ; 
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecMutualIntensityVsX(srTRadExtract& RadExtract)
{//OC15122019
 //This assumes "normal" data alignment in the complex Hermitian "matrix" E(x)*E*(x'), and fills-out "upper triangular" part of it plus the diagonal
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI0 = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long nz = RadAccessData.nz, nx = RadAccessData.nx, ne = RadAccessData.ne;

	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long PerArg = nx << 1; //Row or Column length of the MI matrix represented by real numbers

	long long ie0=0, ie1=0, iz0=0, iz1=0;
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool DontNeedInterp = false;
	if(((iz0 == iz1) || (fabs(InvStepRelArg2) < tolRelArg)) && 
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) DontNeedInterp = true;

	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;

	float *pExInitE0Z0 = pEx0 + iz0PerZ + Two_ie0, *pEzInitE0Z0 = pEz0 + iz0PerZ + Two_ie0;
	float *pExT = pExInitE0Z0, *pEzT = pEzInitE0Z0;
	float *pEx = pExInitE0Z0, *pEz = pEzInitE0Z0;

	float *pExInitE1Z0, *pEzInitE1Z0, *pExInitE0Z1, *pEzInitE0Z1, *pExInitE1Z1, *pEzInitE1Z1;
	float *arExPtrsT[4], *arEzPtrsT[4], *arExPtrs[4], *arEzPtrs[4];
	if(!DontNeedInterp)
	{
		pExInitE1Z0 = pEx0 + iz0PerZ + Two_ie1; pEzInitE1Z0 = pEz0 + iz0PerZ + Two_ie1;
		pExInitE0Z1 = pEx0 + iz1PerZ + Two_ie0; pEzInitE0Z1 = pEz0 + iz1PerZ + Two_ie0;
		pExInitE1Z1 = pEx0 + iz1PerZ + Two_ie1; pEzInitE1Z1 = pEz0 + iz1PerZ + Two_ie1;
		arExPtrsT[0] = pExInitE0Z0; arExPtrsT[1] = pExInitE1Z0; arExPtrsT[2] = pExInitE0Z1; arExPtrsT[3] = pExInitE1Z1;
		arEzPtrsT[0] = pEzInitE0Z0; arEzPtrsT[1] = pEzInitE1Z0; arEzPtrsT[2] = pEzInitE0Z1; arEzPtrsT[3] = pEzInitE1Z1;
		arExPtrs[0] = pExInitE0Z0; arExPtrs[1] = pExInitE1Z0; arExPtrs[2] = pExInitE0Z1; arExPtrs[3] = pExInitE1Z1;
		arEzPtrs[0] = pEzInitE0Z0; arEzPtrs[1] = pEzInitE1Z0; arEzPtrs[2] = pEzInitE0Z1; arEzPtrs[3] = pEzInitE1Z1;
	}

	double iter = 0, *pMeth = RadExtract.pMeth;
	if(pMeth != 0)
	{
		if(*pMeth == 1) iter = *(pMeth + 1);
		else if(*pMeth == 2) iter = -1;
	}

	for(long long ixt=0; ixt<nx; ixt++)
	{
		float *pMI = pMI0 + ixt*PerArg;
		for(long long ix=0; ix<=ixt; ix++)
		{
			if(DontNeedInterp)
			{
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res;
				pEx += PerX; pEz += PerX;
			}
			else
			{
				if(res = MutualIntensityComponentSimpleInterpol2D(arExPtrs, arExPtrsT, arEzPtrs, arEzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, iter, pMI)) return res;
				for(int i=0; i<4; i++) { arExPtrs[i] += PerX; arEzPtrs[i] += PerX;}
			}
			pMI += 2;
		}

		if(DontNeedInterp)
		{
			pEx = pExInitE0Z0;
			pEz = pEzInitE0Z0;
			pExT += PerX; pEzT += PerX;
		}
		else
		{
			arExPtrs[0] = pExInitE0Z0; arExPtrs[1] = pExInitE1Z0; arExPtrs[2] = pExInitE0Z1; arExPtrs[3] = pExInitE1Z1;
			arEzPtrs[0] = pEzInitE0Z0; arEzPtrs[1] = pEzInitE1Z0; arEzPtrs[2] = pEzInitE0Z1; arEzPtrs[3] = pEzInitE1Z1;
			for(int i=0; i<4; i++) { arExPtrsT[i] += PerX; arEzPtrsT[i] += PerX;}
		}
	}

/** //OC06092018
	//Previous version, assuming special alignment (diagonal data first, etc.)
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//long ie0=0, ie1=0, iz0=0, iz1=0;
	long long ie0=0, ie1=0, iz0=0, iz1=0; //OC26042019
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool NeedInterp = true;
	if(((iz0 == iz1) || (fabs(InvStepRelArg2) < tolRelArg)) && 
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) NeedInterp = false;

	int iter = 0, *pMeth = RadExtract.pMeth; //OC14122019
	if(pMeth != 0) { if(*pMeth == 1) iter = *(pMeth + 1);}

	//long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;
	float *pEx_StartForE_iz0, *pEx_StartForE_iz1;
	float *pEz_StartForE_iz0, *pEz_StartForE_iz1;

	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	//long nx = RadAccessData.nx;
	long long nx = RadAccessData.nx; //OC26042019
	long long ixPerX = 0;
	for(long long ix=0; ix<nx; ix++) //OC26042019
	//for(long ix=0; ix<nx; ix++)
	{
		if(NeedInterp)
		{
			pEx_StartForE_iz0 = pEx_StartForX_iz0 + ixPerX; pEx_StartForE_iz1 = pEx_StartForX_iz1 + ixPerX;
			pEz_StartForE_iz0 = pEz_StartForX_iz0 + ixPerX; pEz_StartForE_iz1 = pEz_StartForX_iz1 + ixPerX;
			float *ExPtrs[] = { (pEx_StartForE_iz0 + Two_ie0), (pEx_StartForE_iz0 + Two_ie1),
								(pEx_StartForE_iz1 + Two_ie0), (pEx_StartForE_iz1 + Two_ie1)};
			float *EzPtrs[] = { (pEz_StartForE_iz0 + Two_ie0), (pEz_StartForE_iz0 + Two_ie1),
								(pEz_StartForE_iz1 + Two_ie0), (pEz_StartForE_iz1 + Two_ie1)};
			*(pMI++) = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, 0);
		}
		else
		{
			*(pMI++) = IntensityComponent(pEx_StartForX_iz0 + ixPerX + Two_ie0, pEz_StartForX_iz0 + ixPerX + Two_ie0, PolCom, 0);
		}
		ixPerX += PerX;
	}

	float *pEx_StartForE_iz0_t, *pEx_StartForE_iz1_t;
	float *pEz_StartForE_iz0_t, *pEz_StartForE_iz1_t;
	float **ExPtrs, **EzPtrs, **ExPtrsT, **EzPtrsT;
	//long nx_mi_1 = nx - 1;
	long long nx_mi_1 = nx - 1; //OC26042019
	for(long long ixt=0; ixt<nx_mi_1; ixt++)
	//for(long ixt=0; ixt<nx_mi_1; ixt++)
	{
		long long ixt_PerX = ixt*PerX;
		pEx_StartForE_iz0_t = pEx_StartForX_iz0 + ixt_PerX; 
		pEz_StartForE_iz0_t = pEz_StartForX_iz0 + ixt_PerX;

		if(NeedInterp)
		{
			pEx_StartForE_iz1_t = pEx_StartForX_iz1 + ixt_PerX;
			pEz_StartForE_iz1_t = pEz_StartForX_iz1 + ixt_PerX;
			float *arExPtrs_t[] = { (pEx_StartForE_iz0_t + Two_ie0), (pEx_StartForE_iz0_t + Two_ie1),
									(pEx_StartForE_iz1_t + Two_ie0), (pEx_StartForE_iz1_t + Two_ie1)};
			float *arEzPtrs_t[] = { (pEz_StartForE_iz0_t + Two_ie0), (pEz_StartForE_iz0_t + Two_ie1),
									(pEz_StartForE_iz1_t + Two_ie0), (pEz_StartForE_iz1_t + Two_ie1)};
			ExPtrsT = arExPtrs_t; EzPtrsT = arEzPtrs_t;
		}

		for(long long ix=ixt+1; ix<nx; ix++) //OC26042019
		//for(long ix=ixt+1; ix<nx; ix++)
		{
			long long ix_PerX = ix*PerX;
			pEx_StartForE_iz0 = pEx_StartForX_iz0 + ix_PerX;
			pEz_StartForE_iz0 = pEz_StartForX_iz0 + ix_PerX;

			if(NeedInterp)
			{
				pEx_StartForE_iz1 = pEx_StartForX_iz1 + ix_PerX;
				pEz_StartForE_iz1 = pEz_StartForX_iz1 + ix_PerX;
				float *arExPtrs[] = { (pEx_StartForE_iz0 + Two_ie0), (pEx_StartForE_iz0 + Two_ie1),
									  (pEx_StartForE_iz1 + Two_ie0), (pEx_StartForE_iz1 + Two_ie1)};
				float *arEzPtrs[] = { (pEz_StartForE_iz0 + Two_ie0), (pEz_StartForE_iz0 + Two_ie1),
									  (pEz_StartForE_iz1 + Two_ie0), (pEz_StartForE_iz1 + Two_ie1)};
				ExPtrs = arExPtrs; EzPtrs = arEzPtrs;
				if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, iter, pMI)) return res; //OC14122019
				//if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, pMI)) return res;
			}
			else
			{
				float *pEx = pEx_StartForE_iz0 + Two_ie0, *pExT = pEx_StartForE_iz0_t + Two_ie0;
				float *pEz = pEz_StartForE_iz0 + Two_ie0, *pEzT = pEz_StartForE_iz0_t + Two_ie0;
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res; //OC14122019
				//if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, pMI)) return res;
			}
			pMI += 2;
		}
	}
**/
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecMutualIntensityVsZ(srTRadExtract& RadExtract)
{//OC15122019
 //This assumes "normal" data alignment in the complex Hermitian "matrix" E(y)*E*(y'), and fills-out "upper triangular" part of it plus the diagonal
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
	
	float *pMI0 = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long nz = RadAccessData.nz, nx = RadAccessData.nx, ne = RadAccessData.ne;

	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long PerArg = nz << 1; //Row or Column length of the MI matrix represented by real numbers

	long long ie0=0, ie1=0, ix0=0, ix1=0;
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool DontNeedInterp = false;
	if(((ix0 == ix1) || (fabs(InvStepRelArg2) < tolRelArg)) &&
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) DontNeedInterp = true;

	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;

	float *pExInitE0X0 = pEx0 + ix0PerX + Two_ie0, *pEzInitE0X0 = pEz0 + ix0PerX + Two_ie0;
	float *pExT = pExInitE0X0, *pEzT = pEzInitE0X0;
	float *pEx = pExInitE0X0, *pEz = pEzInitE0X0;

	float *pExInitE1X0, *pEzInitE1X0, *pExInitE0X1, *pEzInitE0X1, *pExInitE1X1, *pEzInitE1X1;
	float *arExPtrsT[4], *arEzPtrsT[4], *arExPtrs[4], *arEzPtrs[4];
	if(!DontNeedInterp)
	{
		pExInitE1X0 = pEx0 + ix0PerX + Two_ie1; pEzInitE1X0 = pEz0 + ix0PerX + Two_ie1;
		pExInitE0X1 = pEx0 + ix1PerX + Two_ie0; pEzInitE0X1 = pEz0 + ix1PerX + Two_ie0;
		pExInitE1X1 = pEx0 + ix1PerX + Two_ie1; pEzInitE1X1 = pEz0 + ix1PerX + Two_ie1;
		arExPtrsT[0] = pExInitE0X0; arExPtrsT[1] = pExInitE1X0; arExPtrsT[2] = pExInitE0X1; arExPtrsT[3] = pExInitE1X1;
		arEzPtrsT[0] = pEzInitE0X0; arEzPtrsT[1] = pEzInitE1X0; arEzPtrsT[2] = pEzInitE0X1; arEzPtrsT[3] = pEzInitE1X1;
		arExPtrs[0] = pExInitE0X0; arExPtrs[1] = pExInitE1X0; arExPtrs[2] = pExInitE0X1; arExPtrs[3] = pExInitE1X1;
		arEzPtrs[0] = pEzInitE0X0; arEzPtrs[1] = pEzInitE1X0; arEzPtrs[2] = pEzInitE0X1; arEzPtrs[3] = pEzInitE1X1;
	}

	double iter = 0, *pMeth = RadExtract.pMeth;
	if(pMeth != 0)
	{
		if(*pMeth == 1) iter = *(pMeth + 1);
		else if(*pMeth == 2) iter = -1;
	}

	for(long long izt=0; izt<nz; izt++)
	{
		float *pMI = pMI0 + izt*PerArg;
		for(long long iz=0; iz<=izt; iz++)
		{
			if(DontNeedInterp)
			{
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res;
				pEx += PerZ; pEz += PerZ;
			}
			else
			{
				if(res = MutualIntensityComponentSimpleInterpol2D(arExPtrs, arExPtrsT, arEzPtrs, arEzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, iter, pMI)) return res;
				for(int i=0; i<4; i++) { arExPtrs[i] += PerZ; arEzPtrs[i] += PerZ;}
			}
			pMI += 2;
		}

		if(DontNeedInterp)
		{
			pEx = pExInitE0X0;
			pEz = pEzInitE0X0;
			pExT += PerZ; pEzT += PerZ;
		}
		else
		{
			arExPtrs[0] = pExInitE0X0; arExPtrs[1] = pExInitE1X0; arExPtrs[2] = pExInitE0X1; arExPtrs[3] = pExInitE1X1;
			arEzPtrs[0] = pEzInitE0X0; arEzPtrs[1] = pEzInitE1X0; arEzPtrs[2] = pEzInitE0X1; arEzPtrs[3] = pEzInitE1X1;
			for(int i=0; i<4; i++) { arExPtrsT[i] += PerZ; arEzPtrsT[i] += PerZ;}
		}
	}

/** //OC10092018
	//Previous version, assuming special alignment (diagonal data first, etc.)
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	//long ie0=0, ie1=0, ix0=0, ix1=0;
	long long ie0=0, ie1=0, ix0=0, ix1=0; //OC26042019
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool NeedInterp = true;
	if(((ix0 == ix1) || (fabs(InvStepRelArg2) < tolRelArg)) && 
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) NeedInterp = false;

	int iter = 0, *pMeth = RadExtract.pMeth; //OC14122019
	if(pMeth != 0) { if(*pMeth == 1) iter = *(pMeth + 1); }

	//long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;

	float *pEx_StartForZ_ix0 = pEx0 + ix0PerX, *pEx_StartForZ_ix1 = pEx0 + ix1PerX;
	float *pEz_StartForZ_ix0 = pEz0 + ix0PerX, *pEz_StartForZ_ix1 = pEz0 + ix1PerX;
	float *pEx_StartForE_ix0, *pEx_StartForE_ix1;
	float *pEz_StartForE_ix0, *pEz_StartForE_ix1;

	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	//long nz = RadAccessData.nz;
	long long nz = RadAccessData.nz; //OC26042019
	long long izPerZ = 0;
	for(long long iz=0; iz<nz; iz++) //OC26042019
	//for(long iz=0; iz<nz; iz++)
	{
		pEx_StartForE_ix0 = pEx_StartForZ_ix0 + izPerZ;
		pEz_StartForE_ix0 = pEz_StartForZ_ix0 + izPerZ;

		if(NeedInterp)
		{
			pEx_StartForE_ix1 = pEx_StartForZ_ix1 + izPerZ;
			pEz_StartForE_ix1 = pEz_StartForZ_ix1 + izPerZ;
			float *ExPtrs[] = { (pEx_StartForE_ix0 + Two_ie0), (pEx_StartForE_ix0 + Two_ie1),
								(pEx_StartForE_ix1 + Two_ie0), (pEx_StartForE_ix1 + Two_ie1)};
			float *EzPtrs[] = { (pEz_StartForE_ix0 + Two_ie0), (pEz_StartForE_ix0 + Two_ie1),
								(pEz_StartForE_ix1 + Two_ie0), (pEz_StartForE_ix1 + Two_ie1)};

			*(pMI++) = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, 0);
		}
		else
		{
			*(pMI++) = IntensityComponent(pEx_StartForE_ix0 + Two_ie0, pEz_StartForE_ix0 + Two_ie0, PolCom, 0);
		}
		izPerZ += PerZ;
	}

	float *pEx_StartForE_ix0_t, *pEx_StartForE_ix1_t;
	float *pEz_StartForE_ix0_t, *pEz_StartForE_ix1_t;
	float **ExPtrs, **EzPtrs, **ExPtrsT, **EzPtrsT;
	//long nz_mi_1 = nz - 1;
	//for(long izt=0; izt<nz_mi_1; izt++)
	long long nz_mi_1 = nz - 1; //OC26042019
	for(long long izt=0; izt<nz_mi_1; izt++)
	{
		long long izt_PerZ = izt*PerZ;
		pEx_StartForE_ix0_t = pEx_StartForZ_ix0 + izt_PerZ;
		pEz_StartForE_ix0_t = pEz_StartForZ_ix0 + izt_PerZ;

		if(NeedInterp)
		{
			pEx_StartForE_ix1_t = pEx_StartForZ_ix1 + izt_PerZ;
			pEz_StartForE_ix1_t = pEz_StartForZ_ix1 + izt_PerZ;
			float *arExPtrs_t[] = { (pEx_StartForE_ix0_t + Two_ie0), (pEx_StartForE_ix0_t + Two_ie1),
									(pEx_StartForE_ix1_t + Two_ie0), (pEx_StartForE_ix1_t + Two_ie1)};
			float *arEzPtrs_t[] = { (pEz_StartForE_ix0_t + Two_ie0), (pEz_StartForE_ix0_t + Two_ie1),
									(pEz_StartForE_ix1_t + Two_ie0), (pEz_StartForE_ix1_t + Two_ie1)};
			ExPtrsT = arExPtrs_t; EzPtrsT = arEzPtrs_t;
		}

		for(long long iz=izt+1; iz<nz; iz++) //OC26042019
		//for(long iz=izt+1; iz<nz; iz++)
		{
			long long iz_PerZ = iz*PerZ;
			pEx_StartForE_ix0 = pEx_StartForZ_ix0 + iz_PerZ;
			pEz_StartForE_ix0 = pEz_StartForZ_ix0 + iz_PerZ;

			if(NeedInterp)
			{
				pEx_StartForE_ix1 = pEx_StartForZ_ix1 + iz_PerZ;
				pEz_StartForE_ix1 = pEz_StartForZ_ix1 + iz_PerZ;
				float *arExPtrs[] = { (pEx_StartForE_ix0 + Two_ie0), (pEx_StartForE_ix0 + Two_ie1),
									  (pEx_StartForE_ix1 + Two_ie0), (pEx_StartForE_ix1 + Two_ie1)};
				float *arEzPtrs[] = { (pEz_StartForE_ix0 + Two_ie0), (pEz_StartForE_ix0 + Two_ie1),
									  (pEz_StartForE_ix1 + Two_ie0), (pEz_StartForE_ix1 + Two_ie1)};
				ExPtrs = arExPtrs; EzPtrs = arEzPtrs;
				if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, iter, pMI)) return res; //OC14122019
				//if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, pMI)) return res;
			}
			else
			{
				float *pEx = pEx_StartForE_ix0 + Two_ie0, *pExT = pEx_StartForE_ix0_t + Two_ie0;
				float *pEz = pEz_StartForE_ix0 + Two_ie0, *pEzT = pEz_StartForE_ix0_t + Two_ie0;
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res; //OC14122019
				//if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, pMI)) return res;
			}
			pMI += 2;
		}
	}
**/
	return 0;
}

//*************************************************************************

//int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ(srTRadExtract& RadExtract)
//int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ(srTRadExtract& RadExtract, gpuUsageArg *pGpuUsage)
int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ(srTRadExtract& RadExtract, void* pvGPU) //HG30112023
{//OC13122019
 //This assumes "normal" data alignment in the complex "matrix" E(x,y)*E*(x',y')
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI0 = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long nz = RadAccessData.nz, nx = RadAccessData.nx, ne = RadAccessData.ne;
	long long nxnz = nx*nz;

	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long PerArg = nxnz << 1; //Row or Column length of the MI matrix represented by real numbers

	long long ie0=0, ie1=0;
	double InvStepRelArg=0;

	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg);

	const double tolRelArg = 1.e-08;
	bool DontNeedInterp = false;
	if((ie0 == ie1) || (fabs(InvStepRelArg) < tolRelArg)) DontNeedInterp = true;

	long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;

	float *pExInit0 = pEx0 + Two_ie0, *pEzInit0 = pEz0 + Two_ie0;
	float *pExInit1 = pEx0 + Two_ie1, *pEzInit1 = pEz0 + Two_ie1;

	float *pExT = pExInit0, *pEzT = pEzInit0;
	float *pEx = pExInit0, *pEz = pEzInit0;

	float *arExPtrsT[] = { pExInit0, pExInit1 }, *arEzPtrsT[] = { pEzInit0, pEzInit1 };
	float *arExPtrs[] = { pExInit0, pExInit1 }, *arEzPtrs[] = { pEzInit0, pEzInit1 };

	double iter = 0, Rx = 0, Rz = 0, xc = 0, zc = 0;
	double inv_iter_p_1 = 1.; //OC08052021
	long long itStart = 0, itEnd = nxnz - 1; //OC03032021
	double *pMeth = RadExtract.pMeth;
	if(pMeth != 0) 
	{ 
		if(*pMeth == 1)
		{
			iter = *(pMeth + 1);
			inv_iter_p_1 = 1./(iter + 1); //OC08052021
		}
		else if(*pMeth == 2) iter = -1;

		Rx = pMeth[2], Rz = pMeth[3], xc = pMeth[4], zc = pMeth[5];

		if(pMeth[18] >= 0) itStart = (long long)pMeth[18]; //OC03032021
		if((pMeth[19] > 0) && (pMeth[19] >= pMeth[18])) itEnd = (long long)pMeth[19]; //OC03032021
	}
	double RobsXorig = RadAccessData.RobsX;
	double RobsZorig = RadAccessData.RobsZ;
	double xc_orig = RadAccessData.xc;
	double zc_orig = RadAccessData.zc;
	bool WfrQuadTermCanBeTreatedAtResizeXorig = RadAccessData.WfrQuadTermCanBeTreatedAtResizeX;
	bool WfrQuadTermCanBeTreatedAtResizeZorig = RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ;
	if((Rx != 0) || (Rz != 0))
	{
		RadAccessData.RobsX = Rx;
		RadAccessData.RobsZ = Rz;
		RadAccessData.xc = xc;
		RadAccessData.zc = zc;
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = (Rx == 0)? false : true;
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = (Rz == 0)? false : true;
		RadAccessData.TreatQuadPhaseTerm('r', PolCom); //, int ieOnly=-1)
	}

	if(itStart > 0) //OC03032021
	{
		//pMI0 += itStart*PerArg; //OC03032021 (to check if this may be necessary)
		long long itStart_PerX = itStart*PerX; //OC03032021 (to check)
		pExT += itStart_PerX;
		pEzT += itStart_PerX;
		arExPtrsT[0] += itStart_PerX; arExPtrsT[1] += itStart_PerX;
		arEzPtrsT[0] += itStart_PerX; arEzPtrsT[1] += itStart_PerX;
	}

	//unsigned long long nPtCSD2 = nxnz*(nxnz + 1); //OC27022021 (twice number of "points" for "triangular" 4D CSD)
	//if(m_useIndE)
	//{
	//	if(m_lenArIndEforCSD != nPtCSD2)
	//	{
	//		if(m_arIndEforCSD != 0) { delete[] m_arIndEforCSD; m_arIndEforCSD = 0;}
	//	}
	//	
	//	if(m_arIndEforCSD == 0)
	//	{
	//		m_arIndEforCSD = new unsigned long[nPtCSD2];
	//		m_lenArIndEforCSD = nPtCSD2;

	//		//Fill-in indexes of E "look-up table":
	//		unsigned long long iPtCSD2 = 0;
	//		for(unsigned long it=0; it<(unsigned long)nxnz; it++)
	//		{
	//			for(unsigned int i=0; i<=it; i++)
	//			{
	//				m_arIndEforCSD[iPtCSD2++] = i;
	//				m_arIndEforCSD[iPtCSD2++] = it;
	//			}
	//		}
	//	}
	//}

#if defined(_WITH_OMPH) //OC23022021 (with OpenMP in Hybrid configuration with MPI)

/**
		long long nxnzE2 = nxnz*nxnz, i, it, ofst, ofst_t;
		float *pMI;

		//#pragma omp parallel for

		for(long long ii=0; ii<nxnzE2; ii++)
		{
			it = (long long)(ii/nxnz);
			i = ii - it*nxnz;
			if(i <= it)
			{
				ofst = i*PerX; ofst_t = it*PerX;
				if(DontNeedInterp)
				{
					pEx = pExInit0 + ofst;
					pExT = pExInit0 + ofst_t;
					pEz = pEzInit0 + ofst;
					pEzT = pEzInit0 + ofst_t;
					pMI = pMI0 + (ii << 1);
					//if(res = MutualIntensityComponentOpt(pEx, pExT, pEz, pEzT, PolCom, iter, pMI[0], pMI[1])) return res;
					//if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res;

					double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
					double ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;
					if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1); ExReT = *pExT; ExImT = *(pExT + 1); }
					if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1); EzReT = *pEzT; EzImT = *(pEzT + 1); }
					double ReMI = 0., ImMI = 0.;

					switch(PolCom)
					{
					case 0: // Lin. Hor.
					{
						ReMI = ExRe*ExReT + ExIm*ExImT;
						ImMI = ExIm*ExReT - ExRe*ExImT;
						break;
					}
					case 1: // Lin. Vert.
					{
						ReMI = EzRe*EzReT + EzIm*EzImT;
						ImMI = EzIm*EzReT - EzRe*EzImT;
						break;
					}
					case 2: // Linear 45 deg.
					{
						double ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
						double ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
						ReMI = 0.5*(ExRe_p_EzRe*ExRe_p_EzReT + ExIm_p_EzIm*ExIm_p_EzImT);
						ImMI = 0.5*(ExIm_p_EzIm*ExRe_p_EzReT - ExRe_p_EzRe*ExIm_p_EzImT);
						break;
					}
					case 3: // Linear 135 deg.
					{
						double ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
						double ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
						ReMI = 0.5*(ExRe_mi_EzRe*ExRe_mi_EzReT + ExIm_mi_EzIm*ExIm_mi_EzImT);
						ImMI = 0.5*(ExIm_mi_EzIm*ExRe_mi_EzReT - ExRe_mi_EzRe*ExIm_mi_EzImT);
						break;
					}
					case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
					//case 4: // Circ. Right
					{
						double ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
						double ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
						ReMI = 0.5*(ExRe_mi_EzIm*ExRe_mi_EzImT + ExIm_p_EzRe*ExIm_p_EzReT);
						ImMI = 0.5*(ExIm_p_EzRe*ExRe_mi_EzImT - ExRe_mi_EzIm*ExIm_p_EzReT);
						break;
					}
					case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
					//case 5: // Circ. Left
					{
						double ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
						double ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
						ReMI = 0.5*(ExRe_p_EzIm*ExRe_p_EzImT + ExIm_mi_EzRe*ExIm_mi_EzReT);
						ImMI = 0.5*(ExIm_mi_EzRe*ExRe_p_EzImT - ExRe_p_EzIm*ExIm_mi_EzReT);
						break;
					}
					case -1: // s0
					{
						ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
						ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
						break;
					}
					case -2: // s1
					{
						ReMI = ExRe*ExReT + ExIm*ExImT - (EzRe*EzReT + EzIm*EzImT);
						ImMI = ExIm*ExReT - ExRe*ExImT - (EzIm*EzReT - EzRe*EzImT);
						break;
					}
					case -3: // s2
					{
						ReMI = ExImT*EzIm + ExIm*EzImT + ExReT*EzRe + ExRe*EzReT;
						ImMI = ExReT*EzIm - ExRe*EzImT - ExImT*EzRe + ExIm*EzReT;
						break;
					}
					case -4: // s3
					{
						ReMI = ExReT*EzIm + ExRe*EzImT - ExImT*EzRe - ExIm*EzReT;
						ImMI = ExIm*EzImT - ExImT*EzIm - ExReT*EzRe + ExRe*EzReT;
						break;
					}
					default: // total mutual intensity, same as s0
					{
						ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
						ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
						break;
						//return CAN_NOT_EXTRACT_MUT_INT;
					}
					}
					if(iter == 0)
					{
						pMI[0] = (float)ReMI;
						pMI[1] = (float)ImMI;
					}
					else if(iter > 0)
					{
						double iter_p_1 = iter + 1; //OC20012020
						//long long iter_p_1 = iter + 1;
						pMI[0] = (float)((pMI[0]*iter + ReMI)/iter_p_1);
						pMI[1] = (float)((pMI[1]*iter + ImMI)/iter_p_1);
					}
					else
					{
						pMI[0] += (float)ReMI;
						pMI[1] += (float)ImMI;
					}
				}
				else
				{
					arExPtrs[0] = pExInit0 + ofst; arExPtrs[1] = pExInit1 + ofst;
					arExPtrsT[0] = pExInit0 + ofst_t; arExPtrsT[1] = pExInit1 + ofst_t;
					arEzPtrs[0] = pEzInit0 + ofst; arEzPtrs[1] = pEzInit1 + ofst;
					arEzPtrsT[0] = pEzInit0 + ofst_t; arEzPtrsT[1] = pEzInit1 + ofst_t;
					//if(res = MutualIntensityComponentSimpleInterpol(arExPtrs, arExPtrsT, arEzPtrs, arEzPtrsT, InvStepRelArg, PolCom, iter, pMI0 + (ii << 1))) return res;
				}
			}
		}
**/

	long long it_PerX, i_PerX;
	float *pMI0_aux, *pMI;

	#pragma omp parallel for //To parallelize only the outmost loop

	for(long long it=0; it<nxnz; it++)
	{
		it_PerX = it*PerX;
		pExT = pExInit0 + it_PerX;
		pEzT = pEzInit0 + it_PerX;

		//float *pMI = pMI0 + it*PerArg;
		pMI0_aux = pMI0 + it*PerArg;
		for(long long i=0; i<=it; i++)
		{
			if(DontNeedInterp)
			{
				//if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res;

				pMI = pMI0_aux + (i << 1);

				i_PerX = i*PerX;
				pEx = pExInit0 + i_PerX;
				pEz = pEzInit0 + i_PerX;

				double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
				double ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;
				if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1); ExReT = *pExT; ExImT = *(pExT + 1); }
				if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1); EzReT = *pEzT; EzImT = *(pEzT + 1); }
				double ReMI = 0., ImMI = 0.;

				switch(PolCom)
				{
				case 0: // Lin. Hor.
				{
					ReMI = ExRe*ExReT + ExIm*ExImT;
					ImMI = ExIm*ExReT - ExRe*ExImT;
					break;
				}
				case 1: // Lin. Vert.
				{
					ReMI = EzRe*EzReT + EzIm*EzImT;
					ImMI = EzIm*EzReT - EzRe*EzImT;
					break;
				}
				case 2: // Linear 45 deg.
				{
					double ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
					double ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
					ReMI = 0.5*(ExRe_p_EzRe*ExRe_p_EzReT + ExIm_p_EzIm*ExIm_p_EzImT);
					ImMI = 0.5*(ExIm_p_EzIm*ExRe_p_EzReT - ExRe_p_EzRe*ExIm_p_EzImT);
					break;
				}
				case 3: // Linear 135 deg.
				{
					double ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
					double ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
					ReMI = 0.5*(ExRe_mi_EzRe*ExRe_mi_EzReT + ExIm_mi_EzIm*ExIm_mi_EzImT);
					ImMI = 0.5*(ExIm_mi_EzIm*ExRe_mi_EzReT - ExRe_mi_EzRe*ExIm_mi_EzImT);
					break;
				}
				case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
				//case 4: // Circ. Right
				{
					double ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
					double ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
					ReMI = 0.5*(ExRe_mi_EzIm*ExRe_mi_EzImT + ExIm_p_EzRe*ExIm_p_EzReT);
					ImMI = 0.5*(ExIm_p_EzRe*ExRe_mi_EzImT - ExRe_mi_EzIm*ExIm_p_EzReT);
					break;
				}
				case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
				//case 5: // Circ. Left
				{
					double ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
					double ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
					ReMI = 0.5*(ExRe_p_EzIm*ExRe_p_EzImT + ExIm_mi_EzRe*ExIm_mi_EzReT);
					ImMI = 0.5*(ExIm_mi_EzRe*ExRe_p_EzImT - ExRe_p_EzIm*ExIm_mi_EzReT);
					break;
				}
				case -1: // s0
				{
					ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
					ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
					break;
				}
				case -2: // s1
				{
					ReMI = ExRe*ExReT + ExIm*ExImT - (EzRe*EzReT + EzIm*EzImT);
					ImMI = ExIm*ExReT - ExRe*ExImT - (EzIm*EzReT - EzRe*EzImT);
					break;
				}
				case -3: // s2
				{
					ReMI = ExImT*EzIm + ExIm*EzImT + ExReT*EzRe + ExRe*EzReT;
					ImMI = ExReT*EzIm - ExRe*EzImT - ExImT*EzRe + ExIm*EzReT;
					break;
				}
				case -4: // s3
				{
					ReMI = ExReT*EzIm + ExRe*EzImT - ExImT*EzRe - ExIm*EzReT;
					ImMI = ExIm*EzImT - ExImT*EzIm - ExReT*EzRe + ExRe*EzReT;
					break;
				}
				default: // total mutual intensity, same as s0
				{
					ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
					ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
					break;
					//return CAN_NOT_EXTRACT_MUT_INT;
				}
				}
				if(iter == 0)
				{
					pMI[0] = (float)ReMI;
					pMI[1] = (float)ImMI;
				}
				else if(iter > 0)
				{
					double iter_p_1 = iter + 1; //OC20012020
					//long long iter_p_1 = iter + 1;
					pMI[0] = (float)((pMI[0]*iter + ReMI)/iter_p_1);
					pMI[1] = (float)((pMI[1]*iter + ImMI)/iter_p_1);
				}
				else
				{
					pMI[0] += (float)ReMI;
					pMI[1] += (float)ImMI;
				}
				//pEx += PerX; pEz += PerX;
			}
			//else
			//{
			//	if(res = MutualIntensityComponentSimpleInterpol(arExPtrs, arExPtrsT, arEzPtrs, arEzPtrsT, InvStepRelArg, PolCom, iter, pMI)) return res;
			//	arExPtrs[0] += PerX; arExPtrs[1] += PerX;
			//	arEzPtrs[0] += PerX; arEzPtrs[1] += PerX;
			//}
			//pMI += 2;
		}

		//if(DontNeedInterp)
		//{
		//	pEx = pExInit0;
		//	pEz = pEzInit0;
		//	//pExT += PerX; pEzT += PerX;
		//}
		//else
		//{
		//	arExPtrs[0] = pExInit0; arExPtrs[0] = pExInit1;
		//	arEzPtrs[0] = pEzInit0; arEzPtrs[0] = pEzInit1;
		//	arExPtrsT[0] += PerX; arExPtrsT[1] += PerX;
		//	arEzPtrsT[0] += PerX; arEzPtrsT[1] += PerX;
		//}
	}

#else

	//if(DontNeedInterp)
	//{
	//	unsigned long it, it_m_1, i;
	//	unsigned long long nPtCSD2 = nxnz*(nxnz + 1); //OC27022021 (twice number of "points" for "triangular" 4D CSD)
	//	unsigned long long nPtCSD = nPtCSD2 >> 1;
	//	for(unsigned long long ii=0; ii<nPtCSD; ii++)
	//	{
	//		unsigned long long ii2 = ii << 1;
	//		//unsigned long i = m_arIndEforCSD[ii2];
	//		//unsigned long it = m_arIndEforCSD[ii2+1];

	//		it = (unsigned long)ceil(0.5*(sqrt(9. + 8.*ii) - 3.));
	//		it_m_1 = it - 1;
	//		i = (unsigned long)(ii - (0.5*it_m_1*(it_m_1 + 3.) + 1.));

	//		pEx = pExInit0 + i*PerX; pExT = pExInit0 + it*PerX;
	//		pEz = pEzInit0 + i*PerX; pEzT = pEzInit0 + it*PerX;

	//		double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
	//		double ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;
	//		if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1); ExReT = *pExT; ExImT = *(pExT + 1); }
	//		if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1); EzReT = *pEzT; EzImT = *(pEzT + 1); }

	//		double ReMI = 0., ImMI = 0.;
	//		switch(PolCom)
	//		{
	//		case 0: // Lin. Hor.
	//		{
	//			ReMI = ExRe*ExReT + ExIm*ExImT;
	//			ImMI = ExIm*ExReT - ExRe*ExImT;
	//			break;
	//		}
	//		case 1: // Lin. Vert.
	//		{
	//			ReMI = EzRe*EzReT + EzIm*EzImT;
	//			ImMI = EzIm*EzReT - EzRe*EzImT;
	//			break;
	//		}
	//		case 2: // Linear 45 deg.
	//		{
	//			double ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
	//			double ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
	//			ReMI = 0.5*(ExRe_p_EzRe*ExRe_p_EzReT + ExIm_p_EzIm*ExIm_p_EzImT);
	//			ImMI = 0.5*(ExIm_p_EzIm*ExRe_p_EzReT - ExRe_p_EzRe*ExIm_p_EzImT);
	//			break;
	//		}
	//		case 3: // Linear 135 deg.
	//		{
	//			double ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
	//			double ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
	//			ReMI = 0.5*(ExRe_mi_EzRe*ExRe_mi_EzReT + ExIm_mi_EzIm*ExIm_mi_EzImT);
	//			ImMI = 0.5*(ExIm_mi_EzIm*ExRe_mi_EzReT - ExRe_mi_EzRe*ExIm_mi_EzImT);
	//			break;
	//		}
	//		case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
	//		//case 4: // Circ. Right
	//		{
	//			double ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
	//			double ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
	//			ReMI = 0.5*(ExRe_mi_EzIm*ExRe_mi_EzImT + ExIm_p_EzRe*ExIm_p_EzReT);
	//			ImMI = 0.5*(ExIm_p_EzRe*ExRe_mi_EzImT - ExRe_mi_EzIm*ExIm_p_EzReT);
	//			break;
	//		}
	//		case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
	//		//case 5: // Circ. Left
	//		{
	//			double ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
	//			double ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
	//			ReMI = 0.5*(ExRe_p_EzIm*ExRe_p_EzImT + ExIm_mi_EzRe*ExIm_mi_EzReT);
	//			ImMI = 0.5*(ExIm_mi_EzRe*ExRe_p_EzImT - ExRe_p_EzIm*ExIm_mi_EzReT);
	//			break;
	//		}
	//		case -1: // s0
	//		{
	//			ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
	//			ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
	//			break;
	//		}
	//		case -2: // s1
	//		{
	//			ReMI = ExRe*ExReT + ExIm*ExImT - (EzRe*EzReT + EzIm*EzImT);
	//			ImMI = ExIm*ExReT - ExRe*ExImT - (EzIm*EzReT - EzRe*EzImT);
	//			break;
	//		}
	//		case -3: // s2
	//		{
	//			ReMI = ExImT*EzIm + ExIm*EzImT + ExReT*EzRe + ExRe*EzReT;
	//			ImMI = ExReT*EzIm - ExRe*EzImT - ExImT*EzRe + ExIm*EzReT;
	//			break;
	//		}
	//		case -4: // s3
	//		{
	//			ReMI = ExReT*EzIm + ExRe*EzImT - ExImT*EzRe - ExIm*EzReT;
	//			ImMI = ExIm*EzImT - ExImT*EzIm - ExReT*EzRe + ExRe*EzReT;
	//			break;
	//		}
	//		default: // total mutual intensity, same as s0
	//		{
	//			ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
	//			ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
	//			break;
	//			//return CAN_NOT_EXTRACT_MUT_INT;
	//		}
	//		}

	//		float *pMI = pMI0 + it*PerArg + (((unsigned long long)i) << 1);
	//		if(iter == 0)
	//		{
	//			pMI[0] = (float)ReMI;
	//			pMI[1] = (float)ImMI;
	//		}
	//		else if(iter > 0)
	//		{
	//			double iter_p_1 = iter + 1; //OC20012020
	//			//long long iter_p_1 = iter + 1;
	//			pMI[0] = (float)((pMI[0]*iter + ReMI)/iter_p_1);
	//			pMI[1] = (float)((pMI[1]*iter + ImMI)/iter_p_1);
	//		}
	//		else
	//		{
	//			pMI[0] += (float)ReMI;
	//			pMI[1] += (float)ImMI;
	//		}
	//	}
	//}
	//else
	//{

	if(DontNeedInterp)
	{
		//GPU_COND(pGpuUsage,
		//	{
		//		ExtractSingleElecMutualIntensityVsXZ_GPU(pEx, pEz, pMI0, nxnz, itStart, itEnd, PerX, iter, PolCom, EhOK, EvOK, pGpuUsage);
		//	})
		//else
		//{
#ifdef _OFFLOAD_GPU //HG30112023
		TGPUUsageArg parGPU(pvGPU); //OC19022024
		if(CAuxGPU::GPUEnabled(&parGPU) && nxnz <=INT_MAX) //HG26022024 Ensure that the GPU version is only called for 32-bit signed int point counts
		//if(CAuxGPU::GPUEnabled(&parGPU)) //OC19022024
		//if(CAuxGPU::GPUEnabled((TGPUUsageArg*)pvGPU))
		{
			ExtractSingleElecMutualIntensityVsXZ_GPU(pEx, pEz, pMI0, (long)nx, (long)nz, (long)ne, (long)itStart, (long)itEnd, (long)PerX, (long)iter, PolCom, EhOK, EvOK, &parGPU); //OC19022024
			//ExtractSingleElecMutualIntensityVsXZ_GPU(pEx, pEz, pMI0, (long)nx, (long)nz, (long)ne, (long)itStart, (long)itEnd, (long)PerX, (long)iter, PolCom, EhOK, EvOK, (TGPUUsageArg*)pvGPU); //OC25012024
			//ExtractSingleElecMutualIntensityVsXZ_GPU(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, PolCom, EhOK, EvOK, (TGPUUsageArg*)pvGPU);
		}
		else
#endif
		{
			for(long long it=itStart; it<=itEnd; it++) //OC16042021 (to enable partial update of MI/CSD)
				//for(long long it=0; it<=(itEnd-itStart); it++) //OC03032021 (to enable partial update of MI/CSD)
				//for(long long it=0; it<nxnz; it++)
			{
				float *pMI = pMI0 + (it - itStart)*PerArg; //OC16042021
				//float *pMI = pMI0 + it*PerArg;
				for(long long i=0; i<=it; i++)
				{
					//if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, iter, pMI)) return res;

					double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
					double ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;
					if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1); ExReT = *pExT; ExImT = *(pExT + 1); }
					if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1); EzReT = *pEzT; EzImT = *(pEzT + 1); }
					double ReMI = 0., ImMI = 0.;

					switch(PolCom)
					{
					case 0: // Lin. Hor.
					{
						ReMI = ExRe*ExReT + ExIm*ExImT;
						ImMI = ExIm*ExReT - ExRe*ExImT;
						break;
					}
					case 1: // Lin. Vert.
					{
						ReMI = EzRe*EzReT + EzIm*EzImT;
						ImMI = EzIm*EzReT - EzRe*EzImT;
						break;
					}
					case 2: // Linear 45 deg.
					{
						double ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
						double ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
						ReMI = 0.5*(ExRe_p_EzRe*ExRe_p_EzReT + ExIm_p_EzIm*ExIm_p_EzImT);
						ImMI = 0.5*(ExIm_p_EzIm*ExRe_p_EzReT - ExRe_p_EzRe*ExIm_p_EzImT);
						break;
					}
					case 3: // Linear 135 deg.
					{
						double ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
						double ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
						ReMI = 0.5*(ExRe_mi_EzRe*ExRe_mi_EzReT + ExIm_mi_EzIm*ExIm_mi_EzImT);
						ImMI = 0.5*(ExIm_mi_EzIm*ExRe_mi_EzReT - ExRe_mi_EzRe*ExIm_mi_EzImT);
						break;
					}
					case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
					//case 4: // Circ. Right
					{
						double ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
						double ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
						ReMI = 0.5*(ExRe_mi_EzIm*ExRe_mi_EzImT + ExIm_p_EzRe*ExIm_p_EzReT);
						ImMI = 0.5*(ExIm_p_EzRe*ExRe_mi_EzImT - ExRe_mi_EzIm*ExIm_p_EzReT);
						break;
					}
					case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
					//case 5: // Circ. Left
					{
						double ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
						double ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
						ReMI = 0.5*(ExRe_p_EzIm*ExRe_p_EzImT + ExIm_mi_EzRe*ExIm_mi_EzReT);
						ImMI = 0.5*(ExIm_mi_EzRe*ExRe_p_EzImT - ExRe_p_EzIm*ExIm_mi_EzReT);
						break;
					}
					case -1: // s0
					{
						ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
						ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
						break;
					}
					case -2: // s1
					{
						ReMI = ExRe*ExReT + ExIm*ExImT - (EzRe*EzReT + EzIm*EzImT);
						ImMI = ExIm*ExReT - ExRe*ExImT - (EzIm*EzReT - EzRe*EzImT);
						break;
					}
					case -3: // s2
					{
						ReMI = ExImT*EzIm + ExIm*EzImT + ExReT*EzRe + ExRe*EzReT;
						ImMI = ExReT*EzIm - ExRe*EzImT - ExImT*EzRe + ExIm*EzReT;
						break;
					}
					case -4: // s3
					{
						ReMI = ExReT*EzIm + ExRe*EzImT - ExImT*EzRe - ExIm*EzReT;
						ImMI = ExIm*EzImT - ExImT*EzIm - ExReT*EzRe + ExRe*EzReT;
						break;
					}
					default: // total mutual intensity, same as s0
					{
						ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
						ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
						break;
						//return CAN_NOT_EXTRACT_MUT_INT;
					}
					}
					if(iter == 0)
					{
						pMI[0] = (float)ReMI;
						pMI[1] = (float)ImMI;
					}
					else if(iter > 0)
					{
						//double iter_p_1 = iter + 1; //OC20012020
						//long long iter_p_1 = iter + 1;
						pMI[0] = (float)((pMI[0]*iter + ReMI)*inv_iter_p_1); //OC08052021
						pMI[1] = (float)((pMI[1]*iter + ImMI)*inv_iter_p_1);
						//pMI[0] = (float)((pMI[0]*iter + ReMI)/iter_p_1);
						//pMI[1] = (float)((pMI[1]*iter + ImMI)/iter_p_1);
					}
					else
					{
						pMI[0] += (float)ReMI;
						pMI[1] += (float)ImMI;
					}

					pEx += PerX; pEz += PerX;
					pMI += 2;
				}

				pEx = pExInit0;
				pEz = pEzInit0;
				pExT += PerX; pEzT += PerX;
			}
		}
	}
	else
	{
		for(long long it=itStart; it<=itEnd; it++) //OC16042021 (to enable partial update of MI/CSD)
		//for(long long it=0; it<=(itEnd-itStart); it++) //OC03032021 (to enable partial update of MI/CSD)
		//for(long long it=0; it<nxnz; it++)
		{
			float *pMI = pMI0 + (it - itStart)*PerArg; //OC16042021
			//float *pMI = pMI0 + it*PerArg;
			for(long long i=0; i<=it; i++)
			{
				if(res = MutualIntensityComponentSimpleInterpol(arExPtrs, arExPtrsT, arEzPtrs, arEzPtrsT, InvStepRelArg, PolCom, iter, pMI)) return res;
				arExPtrs[0] += PerX; arExPtrs[1] += PerX;
				arEzPtrs[0] += PerX; arEzPtrs[1] += PerX;
			
				pMI += 2;
			}

			arExPtrs[0] = pExInit0; arExPtrs[0] = pExInit1;
			arEzPtrs[0] = pEzInit0; arEzPtrs[0] = pEzInit1;
			arExPtrsT[0] += PerX; arExPtrsT[1] += PerX;
			arEzPtrsT[0] += PerX; arEzPtrsT[1] += PerX;
		}
	}

#endif

	if(RadAccessData.WfrQuadTermCanBeTreatedAtResizeX || RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ)
	{
		RadAccessData.TreatQuadPhaseTerm('a', PolCom); //, int ieOnly=-1)
		RadAccessData.RobsX = RobsXorig;
		RadAccessData.RobsZ = RobsZorig;
		RadAccessData.xc = xc_orig;
		RadAccessData.zc = zc_orig;
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeX = WfrQuadTermCanBeTreatedAtResizeXorig;
		RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ = WfrQuadTermCanBeTreatedAtResizeZorig;
	}

/** //OC11092018
	//Previous version, assuming special alignment (diagonal data first, etc.)
	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	long long izPerZ = 0;
	for(long long iz=0; iz<nz; iz++) //OC26042019
	//for(long iz=0; iz<nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;

		float *pEx_St = pEx_StartForX + Two_ie0;
		float *pEz_St = pEz_StartForX + Two_ie0;
		float *pEx_Fi = pEx_StartForX + Two_ie1;
		float *pEz_Fi = pEz_StartForX + Two_ie1;

		for(long long ix=0; ix<nx; ix++) //OC26042019
		//for(long ix=0; ix<nx; ix++)
		{
			if(NeedInterp)
			{
				*(pMI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, 0);
			}
			else
			{
				*(pMI++) = IntensityComponent(pEx_St, pEz_St, PolCom, 0);
			}
			pEx_St += PerX;
			pEz_St += PerX;
			pEx_Fi += PerX;
			pEz_Fi += PerX;
		}
		izPerZ += PerZ;
	}

	float **ExPtrs, **EzPtrs, **ExPtrsT, **EzPtrsT;
	long long nxnz = nx*nz;
	long long nxnz_mi_1 = nxnz - 1;
	for(long long it=0; it<nxnz_mi_1; it++)
	{
		long long it_PerX = it*PerX;
		if(NeedInterp)
		{
			float *pEx0_p_it_PerX = pEx0 + it_PerX;
			float *pEz0_p_it_PerX = pEz0 + it_PerX;
			float *arExPtrs_t[] = { (pEx0_p_it_PerX + Two_ie0), (pEx0_p_it_PerX + Two_ie1)};
			float *arEzPtrs_t[] = { (pEz0_p_it_PerX + Two_ie0), (pEz0_p_it_PerX + Two_ie1)};
			ExPtrsT = arExPtrs_t; EzPtrsT = arEzPtrs_t;
		}

		for(long long i=it+1; i<nxnz; i++)
		{
			long long i_PerX = i*PerX;
			if(NeedInterp)
			{
				float *pEx0_p_i_PerX = pEx0 + i_PerX;
				float *pEz0_p_i_PerX = pEz0 + i_PerX;
				float *arExPtrs[] = { (pEx0_p_i_PerX + Two_ie0), (pEx0_p_i_PerX + Two_ie1)};
				float *arEzPtrs[] = { (pEz0_p_i_PerX + Two_ie0), (pEz0_p_i_PerX + Two_ie1)};
				ExPtrs = arExPtrs; EzPtrs = arEzPtrs;
				if(res = MutualIntensityComponentSimpleInterpol(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg, PolCom, pMI)) return res;
			}
			else
			{
				float *pEx = pEx0 + i_PerX + Two_ie0, *pExT = pEx0 + it_PerX + Two_ie0;
				float *pEz = pEz0 + i_PerX + Two_ie0, *pEzT = pEz0 + it_PerX + Two_ie0;
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, pMI)) return res;
			}
			pMI += 2;
		}
	}
**/
	return 0;
}

//*************************************************************************
/**
int srTRadGenManip::ComputeMultiElecMutualIntensityVsXZ(srTRadExtract& RadExtract, srTTrjDat* pTrjDat)
{//OC23022020 (under development - OC23032020)
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI0 = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long nz = RadAccessData.nz, nx = RadAccessData.nx, ne = RadAccessData.ne;
	long long nxnz = nx*nz;

	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long PerArg = nxnz << 1; //Row or Column length of the complex MI matrix represented by real numbers

	long long ie0=0, ie1=0;
	double InvStepRelArg=0;

	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg);
	const double tolRelArg = 1.e-08;
	bool DontNeedInterpE = false;
	if((ie0 == ie1) || (fabs(InvStepRelArg) < tolRelArg)) DontNeedInterpE = true;

	double Rx = 0., Rz = 0., xc = 0., zc = 0.;
	int nElecXY = 1000, nElecX = 0, nElecY = 0, nElecE = 1;
	char elecIntMeth = 0;
	double *pMeth = RadExtract.pMeth;
	if(pMeth != 0)
	{
		//if(*pMeth == 1) iter = *(pMeth + 1);
		//else if(*pMeth == 2) iter = -1;
		Rx = pMeth[2], Rz = pMeth[3], xc = pMeth[4], zc = pMeth[5];
		if(pMeth[8] <= 0.)
		{
			if(pMeth[7] > 0.) nElecXY = (int)pMeth[7];
			nElecX = 0; nElecY = 0;
		}
		else
		{
			if(pMeth[7] > 0.) nElecX = (int)pMeth[7]; 
			nElecY = (int)pMeth[8]; nElecXY = 0;
		}
		if(pMeth[9] >= 1.) nElecE = (int)pMeth[9]; 
		if(pMeth[10] >= 0.) elecIntMeth = (char)pMeth[10];
	}

	srTEbmDat EbmDat;
	if(res = RadAccessData.OutElectronBeamStruct(EbmDat)) return res;
	if((EbmDat.Mee <= 0.) || (pTrjDat == 0)) nElecE = 1;

	CGenMathRand LocRand;
	srTRadInt RadInt;
	srTWfrSmp AuxSmp;

	double elecEc = EbmDat.Energy, elecSigE = 0.;
	if(nElecE > 1)
	{
		elecSigE = sqrt(EbmDat.Mee)*elecEc;
		RadAccessData.SetObsParamFromWfr(AuxSmp);
	}


	int nElecE_mi_1 = nElecE - 1;
	for(int iElecE=0; iElecE<nElecE; iElecE++)
	{


		if((pTrjDat != 0) && (iElecE < nElecE_mi_1))
		{
			double newElecE = LocRand.NextRandGauss(elecEc, elecSigE, elecIntMeth);
			pTrjDat->EbmDat.Energy = newElecE;
			pTrjDat->EbmDat.SetupGamma(newElecE);
			if(res = pTrjDat->ComputeInterpolatingStructure_FromTrj()) return res;
			
			//RadInt.ComputeElectricFieldFreqDomain(pTrjDat, &AuxSmp, &precElecFld, &wfr, 0);

		}
			
	}


	return 0;
}
**/
//*************************************************************************

int srTRadGenManip::ComputeConvolutedIntensity(srTRadExtract& RadExtract)
{
	int result;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	if((RadAccessData.nx == 1) || (RadAccessData.nz == 1)) return NEEDS_MORE_THAN_ONE_HOR_AND_VERT_OBS_POINT;
	int PT = RadExtract.PlotType;

	srTRadExtract OwnRadExtract = RadExtract;
	float* pTempStorage = 0;

	long Nx = RadAccessData.nx;
	long Nz = RadAccessData.nz;
	//long long Nx = RadAccessData.nx; //OC26042019
	//long long Nz = RadAccessData.nz;

	CGenMathFFT2D FFT2D;
	FFT2D.NextCorrectNumberForFFT(Nx);
	FFT2D.NextCorrectNumberForFFT(Nz);

	//long TotAmOfNewData = (Nx*Nz) << 1;
	long long TotAmOfNewData = (((long long)Nx)*((long long)Nz)) << 1;

	pTempStorage = new float[TotAmOfNewData];
	if(pTempStorage == 0) return MEMORY_ALLOCATION_FAILURE;

	OwnRadExtract.pExtractedData = pTempStorage;
	OwnRadExtract.PlotType = 3; // vs x&z

	//long Ne = RadAccessData.ne;
	long long Ne = RadAccessData.ne; //OC26042019
	OwnRadExtract.ePh = RadAccessData.eStart;

	char SinglePhotonEnergy = ((PT == 1) || (PT == 2) || (PT == 3));
	if(SinglePhotonEnergy) 
	{
		Ne = 1; OwnRadExtract.ePh = RadExtract.ePh;
	}

	int nComp = 1; //OC29122023
	long long ofstCmp = 0; //OC29122023
	float *pDataOrig = RadExtract.pExtractedData; //OC29122023
	double *pdDataOrig = RadExtract.pExtractedDataD; //OC29122023
	if(OwnRadExtract.PolarizCompon == -5) //OC29122023 //Proc. all Stokes components
	{
		nComp = 4;
		OwnRadExtract.PolarizCompon = -1; //Proc. each Stokes component separately
		ofstCmp = ((long long)(RadAccessData.nx))*((long long)(RadAccessData.nz));
	}

	for(long long ie=0; ie<Ne; ie++) //OC26042019
	//for(long ie=0; ie<Ne; ie++)
	{
		if(pDataOrig != 0) RadExtract.pExtractedData = pDataOrig; //OC29122023
		else if(pdDataOrig != 0) RadExtract.pExtractedDataD = pdDataOrig; //OC29122023

		for(int polComp=0; polComp<nComp; polComp++) //OC29122023
		{
			if(result = ExtractSingleElecIntensity2DvsXZ(OwnRadExtract)) return result;
			if(result = ConvoluteWithElecBeamOverTransvCoord(OwnRadExtract.pExtractedData, Nx, Nz)) return result;

			//long ie0 = ie;
			long long ie0 = ie; //OC26042019
			if(SinglePhotonEnergy && (RadAccessData.ne > 1))
			{
				//long ie1;
				long long ie1; //OC26042019
				double RelArgE_Dummy;
				SetupIntCoord('e', OwnRadExtract.ePh, ie0, ie1, RelArgE_Dummy);
			}

			PutConstPhotEnergySliceInExtractPlace(ie0, Nx, Nz, OwnRadExtract, RadExtract);

			if(pDataOrig != 0) RadExtract.pExtractedData += ofstCmp; //OC29122023
			else if(pdDataOrig != 0) RadExtract.pExtractedDataD += ofstCmp; //OC29122023
			//OwnRadExtract.ePh += RadAccessData.eStep; //OC29122023 (moved outside of this loop)

			OwnRadExtract.PolarizCompon--; //OC29122023 (all lines below were just moved into this loop)
		}
		OwnRadExtract.ePh += RadAccessData.eStep; //OC29122023
	}

	if(pTempStorage != 0) delete[] pTempStorage;
	return 0;
}

//*************************************************************************

int srTRadGenManip::ConvoluteWithElecBeamOverTransvCoord(float* DataToConv, long Nx, long Nz)
//int srTRadGenManip::ConvoluteWithElecBeamOverTransvCoord(float* DataToConv, long long Nx, long long Nz)
{
	PadImZerosToRealData(DataToConv, Nx, Nz);

			//srTFFT2D testFFT;
			//testFFT.AuxDebug_TestFFT_Plans();

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	if(RadAccessData.pElecBeam == 0) 
	{
		//srTSend Send; Send.AddWarningMessage(&gVectWarnNos, SINGLE_E_EXTRACTED_INSTEAD_OF_MULTI);
		CErrWarn::AddWarningMessage(&gVectWarnNos, SINGLE_E_EXTRACTED_INSTEAD_OF_MULTI);
		return 0;
	}

	CGenMathFFT2DInfo FFT2DInfo;
	FFT2DInfo.pData = DataToConv;
	FFT2DInfo.Dir = 1;
	FFT2DInfo.xStep = RadAccessData.xStep;
	FFT2DInfo.yStep = RadAccessData.zStep;
	FFT2DInfo.xStart = -(Nx >> 1)*RadAccessData.xStep;
	FFT2DInfo.yStart = -(Nz >> 1)*RadAccessData.zStep;
	FFT2DInfo.Nx = Nx;
	FFT2DInfo.Ny = Nz;
	FFT2DInfo.UseGivenStartTrValues = 0;

	CGenMathFFT2D FFT2D;
	FFT2D.Make2DFFT(FFT2DInfo);

	srTElecBeamMoments ElecBeamMom(RadAccessData.pElecBeam);
	PropagateElecBeamMoments(ElecBeamMom);

	const double Pi = 3.14159265358979;
	const double TwoPi = Pi*2.;
	const double TwoPiE2 = TwoPi*Pi;

	double C2x = TwoPiE2*(ElecBeamMom.Mxx);
	double C2z = TwoPiE2*(ElecBeamMom.Mzz);
	//double C1x = TwoPi*(ElecBeamMom.Mx);
	//double C1z = TwoPi*(ElecBeamMom.Mz);

	float* tData = DataToConv;
	double qz = FFT2DInfo.yStartTr;

	for(long iz=0; iz<Nz; iz++)
	{
		double C2zqzE2 = C2z*qz*qz;
		double qx = FFT2DInfo.xStartTr;

		for(long ix=0; ix<Nx; ix++)
		{
			double Magn = exp(-C2x*qx*qx - C2zqzE2);

			*(tData++) *= (float)Magn; // Re
			*(tData++) *= (float)Magn; // Im

			qx += FFT2DInfo.xStepTr;
		}
		qz += FFT2DInfo.yStepTr;
	}

	FFT2DInfo.pData = DataToConv;
	FFT2DInfo.Dir = -1;

	FFT2DInfo.xStep = FFT2DInfo.xStepTr; FFT2DInfo.xStepTr = RadAccessData.xStep;
	FFT2DInfo.yStep = FFT2DInfo.yStepTr; FFT2DInfo.yStepTr = RadAccessData.zStep;
	FFT2DInfo.xStart = FFT2DInfo.xStartTr; //FFT2DInfo.xStartTr = RadAccessData.xStart;
	FFT2DInfo.yStart = FFT2DInfo.yStartTr; //FFT2DInfo.yStartTr = RadAccessData.zStart;
	FFT2DInfo.UseGivenStartTrValues = 0;

	FFT2D.Make2DFFT(FFT2DInfo);

	return 0;
}

//*************************************************************************

void srTRadGenManip::PadImZerosToRealData(float* pData, long long Nx, long long Nz) //OC26042019
//void srTRadGenManip::PadImZerosToRealData(float* pData, long Nx, long Nz)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	//long NxNzOld = RadAccessData.nx*RadAccessData.nz;
	long long NxNzOld = ((long long)RadAccessData.nx)*((long long)RadAccessData.nz);

	float* pOrig = pData + (NxNzOld - 1);
	float* pTot = pData + ((NxNzOld << 1) - 2);

	//for(long is=0; is<NxNzOld; is++)
	for(long long is=0; is<NxNzOld; is++)
	{
		*pTot = *(pOrig--); *(pTot + 1) = 0.; 
		pTot -= 2;
	}

	if(Nz > RadAccessData.nz)
	{
		//pTot = pData + (NxNzOld << 1); //OC bug fix?
		//for(long iz=NxNzOld; iz<Nz; iz++) //OC bug fix?

        //pTot = pData + ((Nx << 1)*(Nz - 1)); //OC another fix 18jan2005
        //pTot = pData + ((Nx << 1)*(RadAccessData.nz));
        pTot = pData + ((((long long)Nx) << 1)*((long long)RadAccessData.nz));

		for(long long iz=RadAccessData.nz; iz<Nz; iz++) //OC26042019
		//for(long iz=RadAccessData.nz; iz<Nz; iz++)
		{
			for(long long ix=0; ix<Nx; ix++) //OC26042019
			//for(long ix=0; ix<Nx; ix++)
			{
				*(pTot++) = 0.; *(pTot++) = 0.;
			}
		}
	}
	if(Nx > RadAccessData.nx)
	{
		//long TwoNxOld = RadAccessData.nx << 1;
		//long NzOld_mi_1 = RadAccessData.nz - 1;
		//long ShiftPer = (Nx - RadAccessData.nx) << 1;
		long long TwoNxOld = RadAccessData.nx << 1; //OC26042019
		long long NzOld_mi_1 = RadAccessData.nz - 1;
		long long ShiftPer = (Nx - RadAccessData.nx) << 1;

		//float* pOrig_Start = pData + (TwoNxOld*NzOld_mi_1);
		float* pOrig_Start = pData + (((long long)TwoNxOld)*((long long)NzOld_mi_1));
		//long ShiftLen = ShiftPer*NzOld_mi_1;
		long long ShiftLen = ((long long)ShiftPer)*((long long)NzOld_mi_1);

		for(long long iz=0; iz<NzOld_mi_1; iz++) //OC26042019
		//for(long iz=0; iz<NzOld_mi_1; iz++)
		{
			ShiftData(pOrig_Start, TwoNxOld, ShiftLen);
			SetDataToZero(pOrig_Start + TwoNxOld + ShiftLen, ShiftPer);

			ShiftLen -= ShiftPer;
			pOrig_Start -= TwoNxOld;
		}
		SetDataToZero(pData + TwoNxOld, ShiftPer);
	}
}

//*************************************************************************

void srTRadGenManip::PropagateElecBeamMoments(srTElecBeamMoments& ElecBeamMom)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	int i, j;
	//DOUBLE *MatrStrPtrs[4];
	double *MatrStrPtrs[4]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	for(i=0; i<4; i++)
	{
		int i4 = i << 2;
		MatrStrPtrs[i] = RadAccessData.p4x4PropMatr + i4;
	}

	//DOUBLE *Vect = RadAccessData.p4x4PropMatr + 16;
	double *Vect = RadAccessData.p4x4PropMatr + 16; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)

// First-Order Moments
	//DOUBLE InitFirstOrderMom[] = { ElecBeamMom.Mx, ElecBeamMom.Mxp, ElecBeamMom.Mz, ElecBeamMom.Mzp};
	//DOUBLE FinFirstOrderMom[4];
	double InitFirstOrderMom[] = { ElecBeamMom.Mx, ElecBeamMom.Mxp, ElecBeamMom.Mz, ElecBeamMom.Mzp};
	double FinFirstOrderMom[4]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	for(i=0; i<4; i++)
	{
		//DOUBLE Res_i = 0.;
		//DOUBLE *MatrStrPtrs_i = MatrStrPtrs[i];
		double Res_i = 0.; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
		double *MatrStrPtrs_i = MatrStrPtrs[i];
		for(j=0; j<4; j++)
		{
			Res_i += MatrStrPtrs_i[j]*InitFirstOrderMom[j];
		}
		FinFirstOrderMom[i] = Res_i + Vect[i];
	}
	ElecBeamMom.Mx = *FinFirstOrderMom;
	ElecBeamMom.Mxp = FinFirstOrderMom[1];
	ElecBeamMom.Mz = FinFirstOrderMom[2];
	ElecBeamMom.Mzp = FinFirstOrderMom[3];

// Second-Order Moments
	//DOUBLE* pStr = *MatrStrPtrs;
	//DOUBLE a00 = pStr[0], a01 = pStr[1], a02 = pStr[2], a03 = pStr[3];
	double* pStr = *MatrStrPtrs; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	double a00 = pStr[0], a01 = pStr[1], a02 = pStr[2], a03 = pStr[3];
	pStr = MatrStrPtrs[1];
	//DOUBLE a10 = pStr[0], a11 = pStr[1], a12 = pStr[2], a13 = pStr[3];
	double a10 = pStr[0], a11 = pStr[1], a12 = pStr[2], a13 = pStr[3]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	pStr = MatrStrPtrs[2];
	//DOUBLE a20 = pStr[0], a21 = pStr[1], a22 = pStr[2], a23 = pStr[3];
	double a20 = pStr[0], a21 = pStr[1], a22 = pStr[2], a23 = pStr[3]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	pStr = MatrStrPtrs[3];
	//DOUBLE a30 = pStr[0], a31 = pStr[1], a32 = pStr[2], a33 = pStr[3];
	//DOUBLE Str0[] = { a00*a00, 2*a00*a01, a01*a01, a02*a02, 2*a02*a03, a03*a03, 2*a00*a02, 2*a01*a02, 2*a00*a03, 2*a01*a03};
	//DOUBLE Str1[] = { a00*a10, a01*a10 + a00*a11, a01*a11, a02*a12, a03*a12 + a02*a13, a03*a13, a02*a10 + a00*a12, a02*a11 + a01*a12, a03*a10 + a00*a13, a03*a11 + a01*a13};
	//DOUBLE Str2[] = { a10*a10, 2*a10*a11, a11*a11, a12*a12, 2*a12*a13, a13*a13, 2*a10*a12, 2*a11*a12, 2*a10*a13, 2*a11*a13};
	//DOUBLE Str3[] = { a20*a20, 2*a20*a21, a21*a21, a22*a22, 2*a22*a23, a23*a23, 2*a20*a22, 2*a21*a22, 2*a20*a23, 2*a21*a23};
	//DOUBLE Str4[] = { a20*a30, a21*a30 + a20*a31, 21*a31, a22*a32, a23*a32 + a22*a33, a23*a33, a22*a30 + a20*a32, a22*a31 + a21*a32, a23*a30 + a20*a33, a23*a31 + a21*a33};
	//DOUBLE Str5[] = { a30*a30, 2*a30*a31, a31*a31, a32*a32, 2*a32*a33, a33*a33, 2*a30*a32, 2*a31*a32, 2*a30*a33, 2*a31*a33};
	//DOUBLE Str6[] = { a00*a20, a01*a20 + a00*a21, a01*a21, a02*a22, a03*a22 + a02*a23, a03*a23, a02*a20 + a00*a22, a02*a21 + a01*a22, a03*a20 + a00*a23, a03*a21 + a01*a23};
	//DOUBLE Str7[] = { a10*a20, a11*a20 + a10*a21, a11*a21, a12*a22, a13*a22 + a12*a23, a13*a23, a12*a20 + a10*a22, a12*a21 + a11*a22, a13*a20 + a10*a23, a13*a21 + a11*a23};
	//DOUBLE Str8[] = { a00*a30, a01*a30 + a00*a31, a01*a31, a02*a32, a03*a32 + a02*a33, a03*a33, a02*a30 + a00*a32, a02*a31 + a01*a32, a03*a30 + a00*a33, a03*a31 + a01*a33};
	//DOUBLE Str9[] = { a10*a30, a11*a30 + a10*a31, a11*a31, a12*a32, a13*a32 + a12*a33, a13*a33, a12*a30 + a10*a32, a12*a31 + a11*a32, a13*a30 + a10*a33, a13*a31 + a11*a33};
	//DOUBLE *Matr10x10[10];
	double a30 = pStr[0], a31 = pStr[1], a32 = pStr[2], a33 = pStr[3]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	double Str0[] = { a00*a00, 2*a00*a01, a01*a01, a02*a02, 2*a02*a03, a03*a03, 2*a00*a02, 2*a01*a02, 2*a00*a03, 2*a01*a03};
	double Str1[] = { a00*a10, a01*a10 + a00*a11, a01*a11, a02*a12, a03*a12 + a02*a13, a03*a13, a02*a10 + a00*a12, a02*a11 + a01*a12, a03*a10 + a00*a13, a03*a11 + a01*a13};
	double Str2[] = { a10*a10, 2*a10*a11, a11*a11, a12*a12, 2*a12*a13, a13*a13, 2*a10*a12, 2*a11*a12, 2*a10*a13, 2*a11*a13};
	double Str3[] = { a20*a20, 2*a20*a21, a21*a21, a22*a22, 2*a22*a23, a23*a23, 2*a20*a22, 2*a21*a22, 2*a20*a23, 2*a21*a23};
	double Str4[] = { a20*a30, a21*a30 + a20*a31, 21*a31, a22*a32, a23*a32 + a22*a33, a23*a33, a22*a30 + a20*a32, a22*a31 + a21*a32, a23*a30 + a20*a33, a23*a31 + a21*a33};
	double Str5[] = { a30*a30, 2*a30*a31, a31*a31, a32*a32, 2*a32*a33, a33*a33, 2*a30*a32, 2*a31*a32, 2*a30*a33, 2*a31*a33};
	double Str6[] = { a00*a20, a01*a20 + a00*a21, a01*a21, a02*a22, a03*a22 + a02*a23, a03*a23, a02*a20 + a00*a22, a02*a21 + a01*a22, a03*a20 + a00*a23, a03*a21 + a01*a23};
	double Str7[] = { a10*a20, a11*a20 + a10*a21, a11*a21, a12*a22, a13*a22 + a12*a23, a13*a23, a12*a20 + a10*a22, a12*a21 + a11*a22, a13*a20 + a10*a23, a13*a21 + a11*a23};
	double Str8[] = { a00*a30, a01*a30 + a00*a31, a01*a31, a02*a32, a03*a32 + a02*a33, a03*a33, a02*a30 + a00*a32, a02*a31 + a01*a32, a03*a30 + a00*a33, a03*a31 + a01*a33};
	double Str9[] = { a10*a30, a11*a30 + a10*a31, a11*a31, a12*a32, a13*a32 + a12*a33, a13*a33, a12*a30 + a10*a32, a12*a31 + a11*a32, a13*a30 + a10*a33, a13*a31 + a11*a33};
	double *Matr10x10[10];
	//Matr10x10[0] = (DOUBLE*)&Str0; 
	//Matr10x10[1] = (DOUBLE*)&Str1; 
	//Matr10x10[2] = (DOUBLE*)&Str2; 
	//Matr10x10[3] = (DOUBLE*)&Str3; 
	//Matr10x10[4] = (DOUBLE*)&Str4; 
	//Matr10x10[5] = (DOUBLE*)&Str5; 
	//Matr10x10[6] = (DOUBLE*)&Str6; 
	//Matr10x10[7] = (DOUBLE*)&Str7; 
	//Matr10x10[8] = (DOUBLE*)&Str8; 
	//Matr10x10[9] = (DOUBLE*)&Str9; 
	Matr10x10[0] = (double*)&Str0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	Matr10x10[1] = (double*)&Str1; 
	Matr10x10[2] = (double*)&Str2; 
	Matr10x10[3] = (double*)&Str3; 
	Matr10x10[4] = (double*)&Str4; 
	Matr10x10[5] = (double*)&Str5; 
	Matr10x10[6] = (double*)&Str6; 
	Matr10x10[7] = (double*)&Str7; 
	Matr10x10[8] = (double*)&Str8; 
	Matr10x10[9] = (double*)&Str9; 
	//DOUBLE InitSecondOrderMom[] = { ElecBeamMom.Mxx, ElecBeamMom.Mxxp, ElecBeamMom.Mxpxp, ElecBeamMom.Mzz, ElecBeamMom.Mzzp, ElecBeamMom.Mzpzp, ElecBeamMom.Mxz, ElecBeamMom.Mxpz, ElecBeamMom.Mxzp, ElecBeamMom.Mxpzp};
	//DOUBLE FinSecondOrderMom[10];
	double InitSecondOrderMom[] = { ElecBeamMom.Mxx, ElecBeamMom.Mxxp, ElecBeamMom.Mxpxp, ElecBeamMom.Mzz, ElecBeamMom.Mzzp, ElecBeamMom.Mzpzp, ElecBeamMom.Mxz, ElecBeamMom.Mxpz, ElecBeamMom.Mxzp, ElecBeamMom.Mxpzp};
	double FinSecondOrderMom[10]; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	for(i=0; i<10; i++)
	{
		//DOUBLE Res_i = 0.;
		//DOUBLE *MatrStr_i = Matr10x10[i];
		double Res_i = 0.; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
		double *MatrStr_i = Matr10x10[i];
		for(j=0; j<10; j++)
		{
			Res_i += MatrStr_i[j]*InitSecondOrderMom[j];
		}
		FinSecondOrderMom[i] = Res_i;
	}
	ElecBeamMom.Mxx = *FinSecondOrderMom;
	ElecBeamMom.Mxxp = FinSecondOrderMom[1];
	ElecBeamMom.Mxpxp = FinSecondOrderMom[2];
	ElecBeamMom.Mzz = FinSecondOrderMom[3];
	ElecBeamMom.Mzzp = FinSecondOrderMom[4];
	ElecBeamMom.Mzpzp = FinSecondOrderMom[5];
	ElecBeamMom.Mxz = FinSecondOrderMom[6];
	ElecBeamMom.Mxpz = FinSecondOrderMom[7];
	ElecBeamMom.Mxzp = FinSecondOrderMom[8];
	ElecBeamMom.Mxpzp = FinSecondOrderMom[9];
}

//*************************************************************************

void srTRadGenManip::PropagateElecBeamMoments(srTElecBeamMoments& ElecBeamMom, double* p4x4PropMatr, double* p4Vect)
{
	int i, j;
	double *MatrStrPtrs[4];
	for(i=0; i<4; i++)
	{
		int i4 = i << 2;
		MatrStrPtrs[i] = p4x4PropMatr + i4;
	}

	double *Vect = p4Vect;

// First-Order Moments
	double InitFirstOrderMom[] = { ElecBeamMom.Mx, ElecBeamMom.Mxp, ElecBeamMom.Mz, ElecBeamMom.Mzp};
	double FinFirstOrderMom[4];
	for(i=0; i<4; i++)
	{
		double Res_i = 0.;
		double *MatrStrPtrs_i = MatrStrPtrs[i];
		for(j=0; j<4; j++)
		{
			Res_i += MatrStrPtrs_i[j]*InitFirstOrderMom[j];
		}
		FinFirstOrderMom[i] = Res_i + Vect[i];
	}
	ElecBeamMom.Mx = *FinFirstOrderMom;
	ElecBeamMom.Mxp = FinFirstOrderMom[1];
	ElecBeamMom.Mz = FinFirstOrderMom[2];
	ElecBeamMom.Mzp = FinFirstOrderMom[3];

// Second-Order Moments
	double* pStr = *MatrStrPtrs;
	double a00 = pStr[0], a01 = pStr[1], a02 = pStr[2], a03 = pStr[3];
	pStr = MatrStrPtrs[1];
	double a10 = pStr[0], a11 = pStr[1], a12 = pStr[2], a13 = pStr[3];
	pStr = MatrStrPtrs[2];
	double a20 = pStr[0], a21 = pStr[1], a22 = pStr[2], a23 = pStr[3];
	pStr = MatrStrPtrs[3];
	double a30 = pStr[0], a31 = pStr[1], a32 = pStr[2], a33 = pStr[3];
	double Str0[] = { a00*a00, 2*a00*a01, a01*a01, a02*a02, 2*a02*a03, a03*a03, 2*a00*a02, 2*a01*a02, 2*a00*a03, 2*a01*a03};
	double Str1[] = { a00*a10, a01*a10 + a00*a11, a01*a11, a02*a12, a03*a12 + a02*a13, a03*a13, a02*a10 + a00*a12, a02*a11 + a01*a12, a03*a10 + a00*a13, a03*a11 + a01*a13};
	double Str2[] = { a10*a10, 2*a10*a11, a11*a11, a12*a12, 2*a12*a13, a13*a13, 2*a10*a12, 2*a11*a12, 2*a10*a13, 2*a11*a13};
	double Str3[] = { a20*a20, 2*a20*a21, a21*a21, a22*a22, 2*a22*a23, a23*a23, 2*a20*a22, 2*a21*a22, 2*a20*a23, 2*a21*a23};
	double Str4[] = { a20*a30, a21*a30 + a20*a31, 21*a31, a22*a32, a23*a32 + a22*a33, a23*a33, a22*a30 + a20*a32, a22*a31 + a21*a32, a23*a30 + a20*a33, a23*a31 + a21*a33};
	double Str5[] = { a30*a30, 2*a30*a31, a31*a31, a32*a32, 2*a32*a33, a33*a33, 2*a30*a32, 2*a31*a32, 2*a30*a33, 2*a31*a33};
	double Str6[] = { a00*a20, a01*a20 + a00*a21, a01*a21, a02*a22, a03*a22 + a02*a23, a03*a23, a02*a20 + a00*a22, a02*a21 + a01*a22, a03*a20 + a00*a23, a03*a21 + a01*a23};
	double Str7[] = { a10*a20, a11*a20 + a10*a21, a11*a21, a12*a22, a13*a22 + a12*a23, a13*a23, a12*a20 + a10*a22, a12*a21 + a11*a22, a13*a20 + a10*a23, a13*a21 + a11*a23};
	double Str8[] = { a00*a30, a01*a30 + a00*a31, a01*a31, a02*a32, a03*a32 + a02*a33, a03*a33, a02*a30 + a00*a32, a02*a31 + a01*a32, a03*a30 + a00*a33, a03*a31 + a01*a33};
	double Str9[] = { a10*a30, a11*a30 + a10*a31, a11*a31, a12*a32, a13*a32 + a12*a33, a13*a33, a12*a30 + a10*a32, a12*a31 + a11*a32, a13*a30 + a10*a33, a13*a31 + a11*a33};
	double *Matr10x10[10];
	Matr10x10[0] = (double*)&Str0; 
	Matr10x10[1] = (double*)&Str1; 
	Matr10x10[2] = (double*)&Str2; 
	Matr10x10[3] = (double*)&Str3; 
	Matr10x10[4] = (double*)&Str4; 
	Matr10x10[5] = (double*)&Str5; 
	Matr10x10[6] = (double*)&Str6; 
	Matr10x10[7] = (double*)&Str7; 
	Matr10x10[8] = (double*)&Str8; 
	Matr10x10[9] = (double*)&Str9; 
	double InitSecondOrderMom[] = { ElecBeamMom.Mxx, ElecBeamMom.Mxxp, ElecBeamMom.Mxpxp, ElecBeamMom.Mzz, ElecBeamMom.Mzzp, ElecBeamMom.Mzpzp, ElecBeamMom.Mxz, ElecBeamMom.Mxpz, ElecBeamMom.Mxzp, ElecBeamMom.Mxpzp};
	double FinSecondOrderMom[10];
	for(i=0; i<10; i++)
	{
		double Res_i = 0.;
		double *MatrStr_i = Matr10x10[i];
		for(j=0; j<10; j++)
		{
			Res_i += MatrStr_i[j]*InitSecondOrderMom[j];
		}
		FinSecondOrderMom[i] = Res_i;
	}
	ElecBeamMom.Mxx = *FinSecondOrderMom; if(ElecBeamMom.Mxx < 0) ElecBeamMom.Mxx = 0;
	ElecBeamMom.Mxxp = FinSecondOrderMom[1];
	ElecBeamMom.Mxpxp = FinSecondOrderMom[2]; if(ElecBeamMom.Mxpxp < 0) ElecBeamMom.Mxpxp = 0;
	ElecBeamMom.Mzz = FinSecondOrderMom[3]; if(ElecBeamMom.Mzz < 0) ElecBeamMom.Mzz = 0;
	ElecBeamMom.Mzzp = FinSecondOrderMom[4];
	ElecBeamMom.Mzpzp = FinSecondOrderMom[5]; if(ElecBeamMom.Mzpzp < 0) ElecBeamMom.Mzpzp = 0;
	ElecBeamMom.Mxz = FinSecondOrderMom[6];
	ElecBeamMom.Mxpz = FinSecondOrderMom[7];
	ElecBeamMom.Mxzp = FinSecondOrderMom[8];
	ElecBeamMom.Mxpzp = FinSecondOrderMom[9];
}

//*************************************************************************

void srTRadGenManip::PutConstPhotEnergySliceInExtractPlace(long long ie, long long NxSlice, long long NzSlice, srTRadExtract& LocRadExtract, srTRadExtract& RadExtract) //OC26042019
//void srTRadGenManip::PutConstPhotEnergySliceInExtractPlace(long ie, long NxSlice, long NzSlice, srTRadExtract& LocRadExtract, srTRadExtract& RadExtract)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	//float *pGen0 = RadExtract.pExtractedData, *pGenStart;
	float *pGen0 = RadExtract.pExtractedData;
	float *pGenStart = pGen0; //OC160815
	float *pSlice0 = LocRadExtract.pExtractedData;

	//long Nx, Nz;
	long long Nx, Nz; //OC26042019
	double RelArgX, RelArgZ;
	double zAbs, xAbs, xAbsStart;
	//long PerZ_Slice = NxSlice << 1;
	long long PerZ_Slice = NxSlice << 1;

	//long PerX_Gen, PerZ_Gen;
	long long PerX_Gen, PerZ_Gen;

	int PT = RadExtract.PlotType;
	if(PT == 0) // vs e
	{
		Nx = 1; Nz = 1;
		xAbsStart = RadExtract.x; zAbs = RadExtract.z;  
		PerX_Gen = 0; PerZ_Gen = 0;
		pGenStart = pGen0 + ie;
	}
	else if(PT == 1) // vs x
	{
		Nx = RadAccessData.nx; Nz = 1;
		xAbsStart = RadAccessData.xStart; zAbs = RadExtract.z;  
		PerX_Gen = 1; PerZ_Gen = 0;
		pGenStart = pGen0;
	}
	else if(PT == 2) // vs z
	{
		Nx = 1; Nz = RadAccessData.nz;
		xAbsStart = RadExtract.x; zAbs = RadAccessData.zStart;
		PerX_Gen = 0; PerZ_Gen = 1;
		pGenStart = pGen0;
	}
	else if(PT == 3) // vs x&z
	{
		Nx = RadAccessData.nx; Nz = RadAccessData.nz;
		xAbsStart = RadAccessData.xStart; zAbs = RadAccessData.zStart;
		PerX_Gen = 1; PerZ_Gen = RadAccessData.nx;
		pGenStart = pGen0;
	}
	else if(PT == 4) // vs e&x
	{
		Nx = RadAccessData.nx; Nz = 1; // Similar to 1
		xAbsStart = RadAccessData.xStart; zAbs = RadExtract.z;  
		PerX_Gen = RadAccessData.ne; PerZ_Gen = 0;
		pGenStart = pGen0 + ie;
	}
	else if(PT == 5) // vs e&z
	{
		Nx = 1; Nz = RadAccessData.nz; // Similar to 2
		xAbsStart = RadExtract.x; zAbs = RadAccessData.zStart;
		PerX_Gen = 0; PerZ_Gen = RadAccessData.ne;
		pGenStart = pGen0 + ie;
	}
	else if((PT == 6) || (PT == 8)) // vs e&x&z or vs e, integrated over x&z
	{
		Nx = RadAccessData.nx; Nz = RadAccessData.nz; // Similar to 3
		xAbsStart = RadAccessData.xStart; zAbs = RadAccessData.zStart;
		PerX_Gen = RadAccessData.ne; PerZ_Gen = RadAccessData.ne*RadAccessData.nx;
		pGenStart = pGen0 + ie;
	}

	char IntegOverXandZ_IsNeeded = (PT == 8);
	float Sum = 0.;

	//long ix0_Slice, ix1_Slice, iz0_Slice, iz1_Slice;
	long long ix0_Slice, ix1_Slice, iz0_Slice, iz1_Slice; //OC26042019
	//long izPerZ_Gen = 0;
	long long izPerZ_Gen = 0;

	//long Nz_mi_1 = Nz - 1;
	//long Nx_mi_1 = Nx - 1;
	long long Nz_mi_1 = Nz - 1; //OC26042019
	long long Nx_mi_1 = Nx - 1;
	float wx, wz; //, wtot;

	for(long long iz=0; iz<Nz; iz++) //OC26042019
	//for(long iz=0; iz<Nz; iz++)
	{
		wz = 1.;
		if((iz == 0) || (iz == Nz_mi_1)) wz = 0.5;

		SetupIntCoord('z', zAbs, iz0_Slice, iz1_Slice, RelArgZ);

		//float *pSliceStartForX_iz0 = pSlice0 + iz0_Slice*PerZ_Slice;
		//float *pSliceStartForX_iz1 = pSlice0 + iz1_Slice*PerZ_Slice;
		float *pSliceStartForX_iz0 = pSlice0 + (((long long)iz0_Slice)*PerZ_Slice);
		float *pSliceStartForX_iz1 = pSlice0 + (((long long)iz1_Slice)*PerZ_Slice);

		xAbs = xAbsStart;
		//long ixPerX_Gen = 0;
		long long ixPerX_Gen = 0;

		for(long ix=0; ix<Nx; ix++)
		{
            wx = 1.;
            if((ix == 0) || (ix == Nx_mi_1)) wx = 0.5;

			SetupIntCoord('x', xAbs, ix0_Slice, ix1_Slice, RelArgX);

			float *pSlice00 = pSliceStartForX_iz0 + (ix0_Slice << 1);
			float *pSlice10 = pSliceStartForX_iz0 + (ix1_Slice << 1);
			float *pSlice01 = pSliceStartForX_iz1 + (ix0_Slice << 1);
			float *pSlice11 = pSliceStartForX_iz1 + (ix1_Slice << 1);

			float* pGen = pGenStart + izPerZ_Gen + ixPerX_Gen;
			float CurVal = (float)((*pSlice00 - *pSlice01 - *pSlice10 + *pSlice11)*RelArgX*RelArgZ + (*pSlice10 - *pSlice00)*RelArgX + (*pSlice01 - *pSlice00)*RelArgZ + *pSlice00);

			if(IntegOverXandZ_IsNeeded) 
			{
				Sum += wx*wz*CurVal;
			}
			else *pGen = CurVal;

			xAbs += RadAccessData.xStep;
			ixPerX_Gen += PerX_Gen;
		}
		zAbs += RadAccessData.zStep;
		izPerZ_Gen += PerZ_Gen; 
	}

	if(PT == 8)
	{
		*(pGen0 + ie) = (float)((1.E+06)*Sum*(RadAccessData.xStep)*(RadAccessData.zStep));
	}
}

//*************************************************************************

int srTRadGenManip::SetupExtractedWaveData(srTRadExtract& RadExtract, srTWaveAccessData& ExtrWaveData)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	ExtrWaveData.wHndl = RadExtract.wExtractedData;
	ExtrWaveData.hState = RadExtract.hStateExtractedData;

	if(RadExtract.Int_or_Phase != 2) 
	{
		ExtrWaveData.pWaveData = (char*)(RadExtract.pExtractedData);
		*(ExtrWaveData.WaveType) = 'f';
	}
	else
	{
		ExtrWaveData.pWaveData = (char*)(RadExtract.pExtractedDataD);
		*(ExtrWaveData.WaveType) = 'd';
	}

	int PT = RadExtract.PlotType;
	ExtrWaveData.AmOfDims = ((PT >= 0) && (PT < 3))? 1 : ((PT < 6)? 2 : 3);
	char TransvUnitsChar = (RadExtract.TransvPres == 0)? 'm' : 'q';

	char* pUnit0 = *(ExtrWaveData.DimUnits);
	char* pUnit1 = *(ExtrWaveData.DimUnits + 1);
	char* pUnit2 = *(ExtrWaveData.DimUnits + 2);
	for(int i=0; i<3; i++)
	{
		*(pUnit0 + i) = '\0'; *(pUnit1 + i) = '\0'; *(pUnit2 + i) = '\0';
	}
	if(PT == 0) // vs e
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.ne;
		*(ExtrWaveData.DimStartValues) = RadAccessData.eStart;
		*(ExtrWaveData.DimSteps) = RadAccessData.eStep;
		*pUnit0 = 'e'; *(pUnit0 + 1) = 'V';
	}
	else if(PT == 1) // vs x
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.nx;
		*(ExtrWaveData.DimStartValues) = RadAccessData.xStart;
		*(ExtrWaveData.DimSteps) = RadAccessData.xStep;
		*pUnit0 = TransvUnitsChar;
	}
	else if(PT == 2) // vs z
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.nz;
		*(ExtrWaveData.DimStartValues) = RadAccessData.zStart;
		*(ExtrWaveData.DimSteps) = RadAccessData.zStep;
		*pUnit0 = TransvUnitsChar;
	}
	else if(PT == 3) // vs x&z
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.nx; *(ExtrWaveData.DimSizes + 1) = RadAccessData.nz;
		*(ExtrWaveData.DimStartValues) = RadAccessData.xStart; *(ExtrWaveData.DimStartValues + 1) = RadAccessData.zStart; 
		*(ExtrWaveData.DimSteps) = RadAccessData.xStep; *(ExtrWaveData.DimSteps + 1) = RadAccessData.zStep; 
		*pUnit0 = TransvUnitsChar;
		*pUnit1 = TransvUnitsChar;
	}
	else if(PT == 4) // vs e&x
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.ne; *(ExtrWaveData.DimSizes + 1) = RadAccessData.nx;
		*(ExtrWaveData.DimStartValues) = RadAccessData.eStart; *(ExtrWaveData.DimStartValues + 1) = RadAccessData.xStart; 
		*(ExtrWaveData.DimSteps) = RadAccessData.eStep; *(ExtrWaveData.DimSteps + 1) = RadAccessData.xStep; 
		*pUnit0 = 'e'; *(pUnit0 + 1) = 'V';
		*pUnit1 = TransvUnitsChar;
	}
	else if(PT == 5) // vs e&z
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.ne; *(ExtrWaveData.DimSizes + 1) = RadAccessData.nz;
		*(ExtrWaveData.DimStartValues) = RadAccessData.eStart; *(ExtrWaveData.DimStartValues + 1) = RadAccessData.zStart; 
		*(ExtrWaveData.DimSteps) = RadAccessData.eStep; *(ExtrWaveData.DimSteps + 1) = RadAccessData.zStep; 
		*pUnit0 = 'e'; *(pUnit0 + 1) = 'V';
		*pUnit1 = TransvUnitsChar;
	}
	else if(PT == 6) // vs e&x&z
	{
		*(ExtrWaveData.DimSizes) = RadAccessData.ne; *(ExtrWaveData.DimSizes + 1) = RadAccessData.nx; *(ExtrWaveData.DimSizes + 2) = RadAccessData.nz;
		*(ExtrWaveData.DimStartValues) = RadAccessData.eStart; *(ExtrWaveData.DimStartValues + 1) = RadAccessData.xStart; *(ExtrWaveData.DimStartValues + 2) = RadAccessData.zStart;
		*(ExtrWaveData.DimSteps) = RadAccessData.eStep; *(ExtrWaveData.DimSteps + 1) = RadAccessData.xStep; *(ExtrWaveData.DimSteps + 2) = RadAccessData.zStep; 
		*pUnit0 = 'e'; *(pUnit0 + 1) = 'V';
		*pUnit1 = TransvUnitsChar;
		*pUnit2 = TransvUnitsChar;
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::TryToMakePhaseContinuous(srTWaveAccessData& WaveData)
{
	//long Nx = WaveData.DimSizes[0];
	//long Nz = WaveData.DimSizes[1];
	long long Nx = WaveData.DimSizes[0]; //OC26042019
	long long Nz = WaveData.DimSizes[1];

	double* CenterSlice = new double[Nx];
	if(CenterSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	//long ixMid = Nx >> 1, izMid = Nz >> 1;
	//long izMid = Nz >> 1;
	long long izMid = Nz >> 1; //OC26042019

	//long ix, iz;
	long long ix, iz; //OC26042019
	//DOUBLE *pData = (DOUBLE*)(WaveData.pWaveData);
	//DOUBLE *tm = pData + izMid*Nx;
	double *pData = (double*)(WaveData.pWaveData); //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
	double *tm = pData + izMid*Nx;
	double *t = CenterSlice;
	for(ix=0; ix<Nx; ix++) *(t++) = *(tm++);
	TryToMakePhaseContinuous1D(CenterSlice, Nx, -1, 0.);

	double* AuxSlice = new double[Nz];
	if(AuxSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	for(ix=0; ix<Nx; ix++)
	{
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *(t++) = *tm; tm += Nx;}
		TryToMakePhaseContinuous1D(AuxSlice, Nz, izMid, (float)CenterSlice[ix]);
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *tm = *(t++); tm += Nx;}
	}

	if(AuxSlice != 0) delete[] AuxSlice;
	if(CenterSlice != 0) delete[] CenterSlice;
	return 0;
}

//*************************************************************************

void srTRadGenManip::TryToMakePhaseContinuous1D(double* Slice, long long Np, long long i0, float Phi0) //OC26042019
//void srTRadGenManip::TryToMakePhaseContinuous1D(double* Slice, long Np, long i0, float Phi0)
{
	const double TwoPi = 6.2831853071796;
	//OCTEST17082019
	const double cFlip = TwoPi - 2.5; // To steer
	//const double cFlip = TwoPi - 3.; //2.5; // To steer

	float PhToAdd0 = (float)((i0 != -1)? (Phi0 - Slice[i0]) : 0.);

	//long HalfNp = Np >> 1;
	//long OtherHalfNp = Np - HalfNp;
	long long HalfNp = Np >> 1; //OC26042019
	long long OtherHalfNp = Np - HalfNp;

	double PhToAdd = PhToAdd0;
	double *t = Slice + HalfNp - 1; 
	*t += PhToAdd;
	double PrevPh = *(t--);
	for(long long i=0; i<(HalfNp - 1); i++) //OC26042019
	//for(long i=0; i<(HalfNp - 1); i++)
	{
		*t += PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t--);
	}

	PhToAdd = PhToAdd0;
	t = Slice + HalfNp - 1;
	PrevPh = *(t++);
	for(long j=0; j<OtherHalfNp; j++)
	{
		*t += PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t++);
	}
}

//*************************************************************************
//OC13122019 (moved to .h)
//void srTRadGenManip::ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData)
//{
//	if(pData == 0) throw INCORRECT_PARAMS_WFR_COMPON_EXTRACT;
//
//	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
//	srTGenOptElem GenOptElem;
//	srTRadExtract RadExtract(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);
//
//			//DEBUG
//			//srTIgorSend::WarningMessage("ExtractRadiation #1");
//			//END DEBUG
//
//	int res;
//	if(TransvPres != RadAccessData.Pres)
//		if(res = GenOptElem.SetRadRepres(&RadAccessData, char(TransvPres))) throw res;
//	
//	if(RadExtract.Int_or_Phase == 1)
//	{//1- Multi-Elec Intensity
//		if(res = ComputeConvolutedIntensity(RadExtract)) throw res;
//	}
//	else if(RadExtract.Int_or_Phase == 4) //OC
//	{//4- Sigle-Elec. Flux
//		if(res = ExtractFluxFromWfr(RadExtract, 's')) throw res;
//	}
//	else if(RadExtract.Int_or_Phase == 5) //OC
//	{//5- Multi-Elec. Flux
//		if(res = ExtractFluxFromWfr(RadExtract, 'm')) throw res;
//	}
//	else if(RadExtract.Int_or_Phase == 8) //OC06092018
//	{
//		if(res = ExtractSingleElecMutualIntensity(RadExtract)) throw res;
//	}
//	else
//	{
//		if(res = ExtractSingleElecIntensity(RadExtract)) throw res;
//	}
//
//	//OCTEST17082019
//	//if((RadExtract.Int_or_Phase == 2) && (RadExtract.PlotType == 3))
//	//{
//	//	srTWaveAccessData ExtrWaveData;
//	//	if(res = SetupExtractedWaveData(RadExtract, ExtrWaveData)) throw res;
//	//	if(res = TryToMakePhaseContinuous(ExtrWaveData)) throw res;
//	//}
//
//			//DEBUG
//			//srTIgorSend::WarningMessage("ExtractRadiation: ready out");
//			//END DEBUG
//
//	//if(res = SetupExtractedWaveData(RadExtract, ExtractedWaveData)) return res;
//}

//*************************************************************************

void srTRadGenManip::IntProc(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar, int nPar) //OC09032019
//void srTRadGenManip::IntProc(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar)
{//Added to move some simple "loops" from slow Py to C++
	if((pwI1 == 0) || (pwI2 == 0) || (arPar == 0)) throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	char typeProc = (char)(arPar[0]);

	if(typeProc == 1)
	{//Sum-up intensities on same mesh

		long iterAvg = -1; //Default case meaning the simple summation has to be done
		if(nPar > 1) iterAvg = (long)arPar[1]; //Iteration number (>=0) for the case when averaging has to be done

		if((pwI1->DataType == 'm') && (pwI2->DataType == 'm')) srTRadGenManip::MutualIntSumPart(pwI1, pwI2, iterAvg); //OC25042021
		//if((pwI1->DataType == 'm') && (pwI2->DataType == 'm')) srTRadGenManip::MutualIntSumPart(pwI1, pwI2); //OC20042021
		//else if((pwI1->DataType == 0) && (pwI2->DataType == 0)) srTRadGenManip::IntSumPart(pwI1, pwI2); //OC20042021
	}
	else if(typeProc == 2)
	{//Add I2 to I1 with interpolation to the mesh of I1
	
	}
	else if(typeProc == 3) srTRadGenManip::Int2DIntegOverAzim(pwI1, pwI2, arPar + 1, nPar - 1); //OC09032019
	//else if(typeProc == 3) srTRadGenManip::Int2DIntegOverAzim(pwI1, pwI2, arPar + 1);
	else if(typeProc == 4) srTRadGenManip::MutualIntFillHalfHermit(pwI1); //OC05022021
	else if(typeProc == 5) srTRadGenManip::MutualIntTreatComQuadPhTerm(pwI1, arPar + 1, nPar - 1); //OC22062021
	else if(typeProc == 6) srTRadGenManip::CohModesTreatComQuadPhTerm(pwI1, arPar + 1, nPar - 1); //OC26062021
}

//*************************************************************************

void srTRadGenManip::Int2DIntegOverAzim(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar, int nPar) //OC09032019
//void srTRadGenManip::Int2DIntegOverAzim(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar)
{//Inegrate or Average intensity over azimuth in cylindrical coordinates (take input 2D intensity distribution from *pwI2 and save resulting 1D distribution vs r to *pwI1)
	if((pwI1 == 0) || (pwI2 == 0) || (arPar == 0)) throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;
	if((pwI1->AmOfDims != 1) || (pwI1->pWaveData == 0) || (pwI2->AmOfDims != 2) || (pwI2->pWaveData == 0)) throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	bool calcAvg = ((char)arPar[0] == 1);
	char meth = (char)arPar[1];

	int nPhMax = 0;
	double relPrec = 1.e-03;
	if(meth == 1) 
	{
		nPhMax = (int)arPar[2];
		if(nPhMax <= 0) throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;
	}
	else if(meth == 2) relPrec = arPar[2];
	else throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	char ordInterp = (char)arPar[3];
	double phMin = arPar[4];
	double phMax = arPar[5];
	double x0 = arPar[6];
	double y0 = arPar[7];

	int nParMin = 8; //OC09032019 (number of input params without optional ones)
	int nRectToSkip = (int)(0.25*(nPar - nParMin) + 1e-07);
	double *arRectToSkipPar = 0;
	if(nRectToSkip > 0)
	{
		arRectToSkipPar = arPar + nParMin;
		//converting widths to half-widths for convenience of treatment
		for(int i=0; i<(nRectToSkip << 1); i++) arRectToSkipPar[2*i + 1] *= 0.5; 
	}

	const double twoPi = 2*3.141592653589793;
	if(phMin == phMax) phMax = phMin + twoPi;
	double phRange = phMax - phMin;

	float *pfI1=0, *pfI2=0; 
	double *pdI1=0, *pdI2=0; 
	char typeI1 = *(pwI1->WaveType), typeI2 = *(pwI2->WaveType);

	if(typeI1 == 'f') pfI1 = (float*)(pwI1->pWaveData);
	else if(typeI1 == 'd') pdI1 = (double*)(pwI1->pWaveData);
	else throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	float *tfI1 = pfI1;
	double *tdI1 = pdI1;

	if(typeI2 == 'f') pfI2 = (float*)(pwI2->pWaveData);
	else if(typeI2 == 'd') pdI2 = (double*)(pwI2->pWaveData);
	else throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	double xMin = *(pwI2->DimStartValues);
	double xStep = *(pwI2->DimSteps);
	//long nx = *(pwI2->DimSizes);
	long long nx = *(pwI2->DimSizes); //OC26042019
	
	double yMin = *((pwI2->DimStartValues) + 1);
	double yStep = *((pwI2->DimSteps) + 1);
	//long ny = *((pwI2->DimSizes) + 1);
	long long ny = *((pwI2->DimSizes) + 1); //OC26042019

	double *arFuncForAzimInt=0;
	if(nPhMax > 0) arFuncForAzimInt = new double[nPhMax];

	//long nr = *(pwI1->DimSizes);
	long long nr = *(pwI1->DimSizes); //OC26042019
	double rStep = *(pwI1->DimSteps);
	double r = *(pwI1->DimStartValues);
	double rMax = r + rStep*(nr - 1);

	double tolPhMeth2=0;
	double dFdPh1=0, dFdPh2=0;
	srTAuxInt2DIntegOverAzim auxStruct;
	if(meth == 2)
	{
		auxStruct.x0 = x0;
		auxStruct.y0 = y0;
		auxStruct.xMin = xMin;
		auxStruct.xStep = xStep;
		auxStruct.nx = nx;
		auxStruct.yMin = yMin;
		auxStruct.yStep = yStep;
		auxStruct.ny = ny;
		auxStruct.pfI2D = pfI2;
		auxStruct.pdI2D = pdI2;
		auxStruct.ordInterp = ordInterp;
		auxStruct.arRectToSkipPar = arRectToSkipPar;
		auxStruct.nRectToSkip = nRectToSkip;

		tolPhMeth2 = ::fabs(phRange*relPrec);
		if(::fabs(::fabs(phRange) - twoPi) < tolPhMeth2) tolPhMeth2 = 0; //meaning that calculation of edge derivatives vs phi is not necessary
	}
	bool pointToBeUsed = true; //OC09032019
	bool isMoreThanOnePoint = true;
	double curIntOverPh;
	double xx, yy;

	for(int ir=0; ir<nr; ir++)
	{
		curIntOverPh = 0.;
		isMoreThanOnePoint = true; //OC09032019
		if(meth == 1)
		{
			long nPhCur = (long)round((r/rMax)*nPhMax);
			if(nPhCur <= 1)
			{
				if(calcAvg) 
				{
					if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(x0 + r, y0, arRectToSkipPar, nRectToSkip); //OC09032019
					if(pointToBeUsed)
					{
						if(pfI2) curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						else curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
					}
					isMoreThanOnePoint = false; //OC09032019
				}
			}
			else
			{
				long nPhCur_mi_1 = nPhCur - 1;
				double phStep = phRange/nPhCur_mi_1;
				double ph = phMin;
				double *t_ar = arFuncForAzimInt;
				if(pfI2)
				{
					for(int iph=0; iph<nPhCur_mi_1; iph++)
					{
						//x = r*cos(ph); y = r*sin(ph);
						xx = x0 + r*cos(ph); yy = y0 + r*sin(ph); //OC09032019
						if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(xx, yy, arRectToSkipPar, nRectToSkip); //OC09032019
						if(pointToBeUsed) *(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(xx, yy, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						else *(t_ar++) = 0.;
						//*(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						ph += phStep;
					}
					if(::fabs(::fabs(phRange) - twoPi) < fabs(1.e-03*phStep)) *t_ar = *arFuncForAzimInt; //full circle
					else 
					{
						xx = x0 + r*cos(ph); yy = y0 + r*sin(ph); //OC09032019
						if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(xx, yy, arRectToSkipPar, nRectToSkip); //OC09032019
						if(pointToBeUsed) *t_ar = CGenMathInterp::InterpOnRegMesh2d(xx, yy, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						else *t_ar = 0;
						//*t_ar = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
					}
				}
				else
				{
					for(int iph=0; iph<nPhCur_mi_1; iph++)
					{
						//x = r*cos(ph); y = r*sin(ph);
						xx = x0 + r*cos(ph); yy = y0 + r*sin(ph); //OC09032019
						if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(xx, yy, arRectToSkipPar, nRectToSkip); //OC09032019
						if(pointToBeUsed) *(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(xx, yy, xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
						else *(t_ar++) = 0.;
						//*(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
						ph += phStep;
					}
					if(::fabs(::fabs(phRange) - twoPi) < ::fabs(1.e-03*phStep)) *t_ar = *arFuncForAzimInt; //full circle
					else 
					{
						xx = x0 + r*cos(ph); yy = y0 + r*sin(ph); //OC09032019
						if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(xx, yy, arRectToSkipPar, nRectToSkip); //OC09032019
						if(pointToBeUsed) *t_ar = CGenMathInterp::InterpOnRegMesh2d(xx, yy, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						else *t_ar = 0;
						//*t_ar = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
					}
				}
				curIntOverPh = CGenMathMeth::Integ1D_FuncDefByArray(arFuncForAzimInt, nPhCur, phStep);
			}
		}
		else if(meth == 2)
		{
			if(r == 0)
			{
				if(calcAvg) 
				{
					if(arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(x0 + r, y0, arRectToSkipPar, nRectToSkip); //OC09032019
					if(pointToBeUsed)
					{
						if(pfI2) curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						else curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
					}
					isMoreThanOnePoint = false; //OC09032019
				}
			}
			else
			{
				auxStruct.r = r;
				if(tolPhMeth2 > 0) curIntOverPh = CGenMathMeth::Integ1D_Func(&(srTRadGenManip::IntCylCrd), phMin, phMax, relPrec, (void*)(&auxStruct));
				else curIntOverPh = CGenMathMeth::Integ1D_FuncWithEdgeDer(&(srTRadGenManip::IntCylCrd), phMin, phMax, dFdPh1, dFdPh2, relPrec, (void*)(&auxStruct));
			}
		}
		
		if(calcAvg && isMoreThanOnePoint) curIntOverPh /= phRange; //OC09032019
		//if(calcAvg) curIntOverPh /= phRange;

		if(pfI1) *(tfI1++) = (float)curIntOverPh;
		else *(tdI1++) = curIntOverPh;

		r += rStep;
	}

	if(arFuncForAzimInt != 0) delete[] arFuncForAzimInt;
}

//*************************************************************************

void srTRadGenManip::MutualIntSumPart(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, long iterAvg) //OC25042021
//void srTRadGenManip::MutualIntSumPart(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2)
{//OC20042021

	//std::cout << "In srTRadGenManip::MutualIntSumPart: iterAvg=" << iterAvg << "\n"; //DEBUG
	//std::cout.flush(); //DEBUG

	float *pData1F=0;
	double *pData1D=0;
	char cType1 = pwI1->WaveType[0];
	if(cType1 == 'f') pData1F = (float*)(pwI1->pWaveData);
	else if(cType1 == 'd') pData1D = (double*)(pwI1->pWaveData);

	float *pData2F=0;
	double *pData2D=0;
	char cType2 = pwI2->WaveType[0];
	if(cType2 == 'f') pData2F = (float*)(pwI2->pWaveData);
	else if(cType2 == 'd') pData2D = (double*)(pwI2->pWaveData);

	long long nx = pwI1->DimSizes[0], nz = pwI1->DimSizes[1];
	long long nxnz = nx*nz;
	long long perArg = nxnz << 1;

	long long nx2 = pwI2->DimSizes[0], nz2 = pwI2->DimSizes[1];
	if((nx2 != nx) || (nz2 != nz)) throw INCONSISTENT_PARAMS_MI_PROC;

	long long itStart = pwI2->itStart;
	if(itStart < 0) itStart = 0;
	long long itFin = pwI2->itFin;
	if(itFin < 0) itFin = nxnz;

	double aux; //OC27042021

	if((pData1F != 0) && (pData2F != 0))
	{
		if(iterAvg < 0) //OC25042021
		{
			for(long long it=itStart; it<=itFin; it++)
			{
				float *pMI1 = pData1F + it*perArg;
				float *pMI2 = pData2F + (it - itStart)*perArg; //OC23042021
				//float *pMI2 = pData2F + it*perArg;
				for(long long i=0; i<=it; i++)
				{
					*(pMI1++) += *(pMI2++); //Re
					*(pMI1++) += *(pMI2++); //Im
				}
			}
		}
		else //OC25042021
		{
			double invIter_p_1 = 1./(iterAvg + 1.);
			for(long long it=itStart; it<=itFin; it++)
			{
				float *pMI1 = pData1F + it*perArg;
				float *pMI2 = pData2F + (it - itStart)*perArg;
				for(long long i=0; i<=it; i++)
				{
					//else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1

					//OC27042021
					aux = (((double)(*pMI1))*iterAvg + (*(pMI2++)))*invIter_p_1; //re
					*(pMI1++) = (float)aux;
					aux = (((double)(*pMI1))*iterAvg + (*(pMI2++)))*invIter_p_1; //im
					*(pMI1++) = (float)aux;
					//*(pMI1++) = (float)(((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1); //Re
					//*(pMI1++) = (float)(((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1); //Im
				}
			}
		}
	}
	else if((pData1D != 0) && (pData2D != 0))
	{
		if(iterAvg < 0) //OC25042021
		{
			for(long long it=itStart; it<=itFin; it++)
			{
				double *pMI1 = pData1D + it*perArg;
				double *pMI2 = pData2D + (it - itStart)*perArg; //OC23042021
				//double *pMI2 = pData2D + it*perArg;
				for(long long i=0; i<=it; i++)
				{
					*(pMI1++) += *(pMI2++);
					*(pMI1++) += *(pMI2++);
				}
			}
		}
		else //OC25042021
		{
			double invIter_p_1 = 1./(iterAvg + 1.);
			for(long long it=itStart; it<=itFin; it++)
			{
				double *pMI1 = pData1D + it*perArg;
				double *pMI2 = pData2D + (it - itStart)*perArg; //OC23042021
				for(long long i=0; i<=it; i++)
				{
					//else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1

					//OC27042021
					aux = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //re
					*(pMI1++) = aux;
					aux = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //im
					*(pMI1++) = aux;
					//*(pMI1++) = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //Re
					//*(pMI1++) = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //Im
				}
			}
		}
	}
	else if((pData1D != 0) && (pData2F != 0))
	{
		if(iterAvg < 0) //OC25042021
		{
			for(long long it=itStart; it<=itFin; it++)
			{
				double *pMI1 = pData1D + it*perArg;
				float *pMI2 = pData2F + (it - itStart)*perArg; //OC23042021
				//float *pMI2 = pData2F + it*perArg;
				for(long long i=0; i<=it; i++)
				{
					*(pMI1++) += *(pMI2++);
					*(pMI1++) += *(pMI2++);
				}
			}
		}
		else
		{
			double invIter_p_1 = 1./(iterAvg + 1.);
			for(long long it=itStart; it<=itFin; it++)
			{
				double *pMI1 = pData1D + it*perArg;
				float *pMI2 = pData2F + (it - itStart)*perArg; //OC23042021
				for(long long i=0; i<=it; i++)
				{
					//else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1

					//OC27042021
					aux = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //re
					*(pMI1++) = aux;
					aux = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //im
					*(pMI1++) = aux;
					//*(pMI1++) = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //Re
					//*(pMI1++) = ((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1; //Im
				}
			}
		}
	}
	else if((pData2D != 0) && (pData1F != 0))
	{
		if(iterAvg < 0) //OC25042021
		{
			for(long long it=itStart; it<=itFin; it++)
			{
				float *pMI1 = pData1F + it*perArg;
				double *pMI2 = pData2D + (it - itStart)*perArg; //OC23042021
				//double *pMI2 = pData2D + it*perArg;
				for(long long i=0; i<=it; i++)
				{
					*(pMI1++) += (float)(*(pMI2++));
					*(pMI1++) += (float)(*(pMI2++));
				}
			}
		}
		else
		{
			double invIter_p_1 = 1./(iterAvg + 1.);
			for(long long it=itStart; it<=itFin; it++)
			{
				float *pMI1 = pData1F + it*perArg;
				double *pMI2 = pData2D + (it - itStart)*perArg; //OC23042021
				for(long long i=0; i<=it; i++)
				{
					//else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1

					//OC27042021
					aux = (((double)(*pMI1))*iterAvg + (*(pMI2++)))*invIter_p_1; //re
					*(pMI1++) = (float)aux;
					aux = (((double)(*pMI1))*iterAvg + (*(pMI2++)))*invIter_p_1; //im
					*(pMI1++) = (float)aux;
					//*(pMI1++) = (float)(((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1); //Re
					//*(pMI1++) = (float)(((*pMI1)*iterAvg + (*(pMI2++)))*invIter_p_1); //Im
				}
			}
		}
	}
}

//*************************************************************************

void srTRadGenManip::MutualIntFillHalfHermit(srTWaveAccessData* pwI)
{//OC06022021
 //This assumes "normal" data alignment in the complex Hermitian "matrix" E(x,y)*E*(x',y') and fills out half of it above the diagonal, using complex conjugation
	float *pDataF=0;
	double *pDataD=0;
	char cType = pwI->WaveType[0];
	if(cType == 'f') pDataF = (float*)(pwI->pWaveData);
	else if(cType == 'd') pDataD = (double*)(pwI->pWaveData);

	long long nx = pwI->DimSizes[0], nz = pwI->DimSizes[1];
	long long nxnz = nx*nz;
	long long perArg = nxnz << 1;

	if(pDataF != 0)
	{
		for(long long it=1; it<nxnz; it++)
		{
			float *pMI = pDataF + it*perArg;
			for(long long i=0; i<it; i++)
			{
				float reMI = *(pMI++);
				float imMI = *(pMI++);
				float *pMIt = pDataF + (it*2 + i*perArg);
				*(pMIt++) = reMI;
				*pMIt = -imMI; //Hermitian matrix property
			}
		}
	}
	else if(pDataD != 0)
	{
		for(long long it=1; it<nxnz; it++)
		{
			double *pMI = pDataD + it*perArg;
			for(long long i=0; i<it; i++)
			{
				double reMI = *(pMI++);
				double imMI = *(pMI++);
				double *pMIt = pDataD + (it*2 + i*perArg);
				*(pMIt++) = reMI;
				*pMIt = -imMI; //Hermitian matrix property
			}
		}
	}
}

//*************************************************************************

void srTRadGenManip::MutualIntTreatComQuadPhTerm(srTWaveAccessData* pwI, double* arPar, int nPar) 
{//OC22062021
	float *pDataF=0;
	double *pDataD=0;
	char cType = pwI->WaveType[0];
	if(cType == 'f') pDataF = (float*)(pwI->pWaveData);
	else if(cType == 'd') pDataD = (double*)(pwI->pWaveData);

	//Call from Py and params: srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [5, -1, wfrA.Rx, wfrA.Ry, wfrA.xc, wfrA.yc])
	char AddOrRem = 'r';
	if(arPar[0] > 0.) AddOrRem = 'a';
	double Rx = arPar[1], Rz = arPar[2];
	double xc = 0., zc = 0.;
	if(nPar >= 5) { xc = arPar[3]; zc = arPar[4];}

	if((Rx == 0.) || (Rz == 0.))
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, ZERO_WFR_RAD_CURV_PH_TERM_NOT_TREATED);
	}
	if((Rx == 0.) && (Rz == 0.)) return;

	long long nx = pwI->DimSizes[0], nz = pwI->DimSizes[1];

	double xStart = pwI->DimStartValues[0];
	double zStart = pwI->DimStartValues[1]; //OC22062021: here 'z' corresponds to vertical direction, "in line" with the convention in other similar functions
	double xStep = pwI->DimSteps[0];
	double zStep = pwI->DimSteps[1];
	double phEnEv = pwI->DimStartValues[2]; //OC22062021: this processing assumes monochromatic MI/CSD, "for the moment"

	const double Pi = 3.14159265358979;
	const double Const = Pi*1.E+06/1.23984186; // Assumes m and eV
	double ConstRxE = 0., ConstRzE = 0.;
	if(Rx != 0.) ConstRxE = Const*phEnEv/Rx;
	if(Rz != 0.) ConstRzE = Const*phEnEv/Rz;

	if(AddOrRem == 'r') { ConstRxE = -ConstRxE; ConstRzE = -ConstRzE;}

	double phTermZ, phTermZt, phTerm, phTermT;

	if(pDataF != 0)
	{
		float *pMI = pDataF;
		float reMI, imMI;
		double zt = zStart - zc;
		float cosPh, sinPh;
		for(long long izt=0; izt<nz; izt++)
		{
			double xt = xStart - xc;
			phTermZt = -ConstRzE*zt*zt;
			for(long long ixt=0; ixt<nx; ixt++)
			{
				double z = zStart - zc;
				phTermT = phTermZt - ConstRxE*xt*xt;
				for(long long iz=0; iz<nz; iz++)
				{
					double x = xStart - xc;
					phTermZ = ConstRzE*z*z;
					for(long long ix=0; ix<nx; ix++)
					{
						phTerm = phTermZ + ConstRxE*x*x + phTermT;
						CosAndSin(phTerm, cosPh, sinPh);

						reMI = *pMI;
						imMI = *(pMI+1);
						*(pMI++) = reMI*cosPh - imMI*sinPh; //new reMI
						*(pMI++) = reMI*sinPh + imMI*cosPh; //new imMI

						x += xStep;
					}
					z += zStep;
				}
				xt += xStep;
			}
			zt += zStep;
		}
	}
	if(pDataD != 0)
	{
		double *pMI = pDataD;
		double reMI, imMI;
		double zt = zStart - zc;
		double cosPh, sinPh;
		for(long long izt=0; izt<nz; izt++)
		{
			double xt = xStart - xc;
			phTermZt = -ConstRzE*zt*zt;
			for(long long ixt=0; ixt<nx; ixt++)
			{
				double z = zStart - zc;
				phTermT = phTermZt - ConstRxE*xt*xt;
				for(long long iz=0; iz<nz; iz++)
				{
					double x = xStart - xc;
					phTermZ = ConstRzE*z*z;
					for(long long ix=0; ix<nx; ix++)
					{
						phTerm = phTermZ + ConstRxE*x*x + phTermT;
						cosPh = cos(phTerm);
						sinPh = sin(phTerm);
						//CosAndSin(phTerm, cosPh, sinPh);

						reMI = *pMI;
						imMI = *(pMI+1);
						*(pMI++) = reMI*cosPh - imMI*sinPh; //new reMI
						*(pMI++) = reMI*sinPh + imMI*cosPh; //new imMI

						x += xStep;
					}
					z += zStep;
				}
				xt += xStep;
			}
			zt += zStep;
		}
	}
}

//*************************************************************************

void srTRadGenManip::CohModesTreatComQuadPhTerm(srTWaveAccessData* pwI, double* arPar, int nPar)
{//OC26062021
	float *pDataF=0;
	double *pDataD=0;
	char cType = pwI->WaveType[0];
	if(cType == 'f') pDataF = (float*)(pwI->pWaveData);
	else if(cType == 'd') pDataD = (double*)(pwI->pWaveData);

	//Call from Py and params: srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [6, n_modes, -1, wfrA.Rx, wfrA.Ry, wfrA.xc, wfrA.yc])
	int nModes = (int)arPar[0];

	char AddOrRem = 'r';
	if(arPar[1] > 0.) AddOrRem = 'a';
	double Rx = arPar[2], Rz = arPar[3];
	double xc = 0., zc = 0.;
	if(nPar >= 6) { xc = arPar[4]; zc = arPar[5]; }

	if((Rx == 0.) || (Rz == 0.))
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, ZERO_WFR_RAD_CURV_PH_TERM_NOT_TREATED);
	}
	if((Rx == 0.) && (Rz == 0.)) return;

	long long ne=1, nx, nz;
	double eStart, eStep=0., xStart, xStep, zStart, zStep;
	if(pwI->AmOfDims == 2)
	{
		nx = pwI->DimSizes[0]; nz = pwI->DimSizes[1];
		eStart = pwI->DimStartValues[2];
		xStart = pwI->DimStartValues[0]; zStart = pwI->DimStartValues[1]; //OC22062021: here 'z' corresponds to vertical direction, "in line" with the convention in other similar functions
		xStep = pwI->DimSteps[0]; zStep = pwI->DimSteps[1];
	}
	else if(pwI->AmOfDims == 3)
	{
		ne = pwI->DimSizes[0]; nx = pwI->DimSizes[1]; nz = pwI->DimSizes[2];
		eStart = pwI->DimStartValues[0]; xStart = pwI->DimStartValues[1]; zStart = pwI->DimStartValues[2];
		eStep = pwI->DimSteps[0]; xStep = pwI->DimSteps[1]; zStep = pwI->DimSteps[2];
	}
	else throw SRWL_INCORRECT_PARAM_FOR_WFR_PROC;

	const double Pi = 3.14159265358979;
	const double Const = Pi*1.E+06/1.23984186; // Assumes m and eV
	double ConstRx = 0., ConstRz = 0.;
	if(Rx != 0.) ConstRx = Const/Rx;
	if(Rz != 0.) ConstRz = Const/Rz;

	if(AddOrRem == 'r') { ConstRx = -ConstRx; ConstRz = -ConstRz;}

	long long PerX = ne << 1;
	long long PerZ = PerX*nx;
	long long PerM = PerZ*nz;

	double phTermZ, phTerm;
	double x, z, ConstRxE, ConstRzE;
	for(int im=0; im<nModes; im++)
	{
		if(pDataF != 0)
		{
			float cosPh, sinPh;
			float *pE0 = pDataF + im*PerM;
			double ePh = eStart;
			for(long long ie=0; ie<ne; ie++)
			{
				long long Two_ie = ie << 1;
				ConstRxE = ConstRx*ePh;
				ConstRzE = ConstRz*ePh;
				z = zStart - zc;
				for(long long iz=0; iz<nz; iz++)
				{
					long long izPerZ = iz*PerZ;
					float *pE_StartForX = pE0 + izPerZ;
					phTermZ = ConstRzE*z*z;
					x = xStart - xc;
					for(long long ix=0; ix<nx; ix++)
					{
						long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
						phTerm = phTermZ + ConstRxE*x*x;

						CosAndSin(phTerm, cosPh, sinPh);

						float *pE_Re = pE_StartForX + ixPerX_p_Two_ie;
						float *pE_Im = pE_Re + 1;
						double E_ReNew = (*pE_Re)*cosPh - (*pE_Im)*sinPh;
						double E_ImNew = (*pE_Re)*sinPh + (*pE_Im)*cosPh;
						*pE_Re = (float)E_ReNew; *pE_Im = (float)E_ImNew;

						x += xStep;
					}
					z += zStep;
				}
				ePh += eStep;
			}
		}
		else if(pDataD != 0)
		{
			double cosPh, sinPh;
			double *pE0 = pDataD + im*PerM;
			double ePh = eStart;
			for(long long ie=0; ie<ne; ie++)
			{
				long long Two_ie = ie << 1;
				ConstRxE = ConstRx*ePh;
				ConstRzE = ConstRz*ePh;
				z = zStart - zc;
				for(long long iz=0; iz<nz; iz++)
				{
					long long izPerZ = iz*PerZ;
					double *pE_StartForX = pE0 + izPerZ;
					phTermZ = ConstRzE*z*z;
					x = xStart - xc;
					for(long long ix=0; ix<nx; ix++)
					{
						long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
						phTerm = phTermZ + ConstRxE*x*x;
						cosPh = cos(phTerm);
						sinPh = sin(phTerm);

						double *pE_Re = pE_StartForX + ixPerX_p_Two_ie;
						double *pE_Im = pE_Re + 1;
						double E_ReNew = (*pE_Re)*cosPh - (*pE_Im)*sinPh;
						double E_ImNew = (*pE_Re)*sinPh + (*pE_Im)*cosPh;
						*pE_Re = E_ReNew; *pE_Im = E_ImNew;

						x += xStep;
					}
					z += zStep;
				}
				ePh += eStep;
			}
		}
	}
}

//*************************************************************************

void srTRadGenManip::ComponInteg(srTDataMD* pIntensOrigData, srTDataMD* pIntegParData, srTDataMD* pIntegResData)
{ 
	if((pIntensOrigData == 0) || (pIntegParData == 0) || (pIntegResData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	//{IntegType, eMin, eMax, xMin, xMax, zMin, zMax}
	double DummyImg;
	int IntegType = (int)srTUtiDataMD::ExtractValue(pIntegParData, 0, &DummyImg);
	double eMin = srTUtiDataMD::ExtractValue(pIntegParData, 1, &DummyImg);
	double eMax = srTUtiDataMD::ExtractValue(pIntegParData, 2, &DummyImg);
	double xMin = srTUtiDataMD::ExtractValue(pIntegParData, 3, &DummyImg);
	double xMax = srTUtiDataMD::ExtractValue(pIntegParData, 4, &DummyImg);
	double zMin = srTUtiDataMD::ExtractValue(pIntegParData, 5, &DummyImg);
	double zMax = srTUtiDataMD::ExtractValue(pIntegParData, 6, &DummyImg);

	//"Photon Energy;Horizontal Position;Vertical Position;Hor. + Vert. Pos.;En. + Hor. Pos.;En. + Vert. Pos.;En. + Hor. + Vert. Pos."
	if(IntegType == 1) 
	{
        ComponIntegVsPhotEn(pIntensOrigData, eMin, eMax, pIntegResData);
	}
	else if(IntegType == 2)
	{
        ComponIntegVsHorOrVertPos(pIntensOrigData, 'x', xMin, xMax, pIntegResData);
	}
	else if(IntegType == 3)
	{
        ComponIntegVsHorOrVertPos(pIntensOrigData, 'z', zMin, zMax, pIntegResData);
	}
	else if(IntegType == 4)
	{
        ComponIntegVsHorAndVertPos(pIntensOrigData, xMin, xMax, zMin, zMax, pIntegResData);
	}
	else if(IntegType == 5)
	{
        ComponIntegVsPhotEnAndPos(pIntensOrigData, 'x', xMin, xMax, pIntegResData);
	}
	else if(IntegType == 6)
	{
        ComponIntegVsPhotEnAndPos(pIntensOrigData, 'z', zMin, zMax, pIntegResData);
	}
	else if(IntegType == 7)
	{
        ComponIntegVsPhotEnAndHorAndVertPos(pIntensOrigData, eMin, eMax, xMin, xMax, zMin, zMax, pIntegResData);
	}
}

//*************************************************************************

void srTRadGenManip::ComponIntegVsPhotEn(srTDataMD* pIntensOrigData, double eMin, double eMax, srTDataMD* pIntegResData)
{
	if((pIntensOrigData == 0) || (pIntegResData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	if(pIntensOrigData->DataUnits[0] != '\0')
	{
		if(strcmp(pIntensOrigData->DataUnits, "W/eV") == 0)
		{
			ConstPhotEnInteg = 1;
            strcpy(pIntegResData->DataUnits, "W");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "W/eV/mm^2") == 0) || (strcmp(pIntensOrigData->DataUnits, "W/eV/mm2") == 0))
		{
			ConstPhotEnInteg = 1;
            strcpy(pIntegResData->DataUnits, "W/mm^2");
		}
		else if(strcmp(pIntensOrigData->DataUnits, "W/eV/mm") == 0)
		{
			ConstPhotEnInteg = 1;
            strcpy(pIntegResData->DataUnits, "W/mm");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw") == 0) || (strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw") == 0))
		{
            strcpy(pIntegResData->DataUnits, "W");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm^2") == 0) || 
				(strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm^2") == 0) ||
				(strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm2") == 0) || 
				(strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm2") == 0))
		{
            strcpy(pIntegResData->DataUnits, "W/mm^2");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm") == 0) || (strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm") == 0))
		{
            strcpy(pIntegResData->DataUnits, "W/mm");
		}
	}

	//long DimSizesOrig[10];
	long long DimSizesOrig[10]; //OC26042019
	int AmOfDimOrig = srTUtiDataMD::ExtractDimSizes(pIntensOrigData, DimSizesOrig);
	if(AmOfDimOrig <= 0) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	//long DimSizesFin[10];
	long long DimSizesFin[10]; //OC26042019
	int AmOfDimFin = srTUtiDataMD::ExtractDimSizes(pIntegResData, DimSizesFin);
	if((AmOfDimFin <= 0) || (AmOfDimFin < (AmOfDimOrig - 1))) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	//long Ne = pIntensOrigData->DimSizes[0];
	long long Ne = pIntensOrigData->DimSizes[0]; //OC26042019
	double StartE = srTUtiDataMD::ExtractDimStartValue(pIntensOrigData, 0);
	double StepE = srTUtiDataMD::ExtractDimStep(pIntensOrigData, 0);

	double InvStepE = 0;
	if(StepE != 0) InvStepE = 1./StepE;

	double dIndStartE = (eMin - StartE)*InvStepE;
	//long IndStartE = (long)dIndStartE;
	long long IndStartE = (long long)dIndStartE; //OC26042019
	if((dIndStartE - IndStartE) >= 0.5) IndStartE++; //to improve
	if(IndStartE < 0) IndStartE = 0;
	if(IndStartE >= Ne) IndStartE = Ne - 1;

	double dIndEndE = (eMax - StartE)*InvStepE;
	//long IndEndE = (long)dIndEndE;
	long long IndEndE = (long long)dIndEndE; //OC26042019
	if((dIndEndE - IndEndE) >= 0.5) IndEndE++; //to improve
	if(IndEndE < 0) IndEndE = 0;
	if(IndEndE >= Ne) IndEndE = Ne - 1;

	//long EffNe = IndEndE - IndStartE;
	long long EffNe = IndEndE - IndStartE; //OC26042019
	if(EffNe < 0) EffNe = 0;
	else if(EffNe > Ne) EffNe = Ne;

	if(AmOfDimOrig > 1)
	{
		for(int k=0; k<AmOfDimFin; k++)
		{
			//int FinNp = DimSizesFin[AmOfDimFin - k - 1];
			long long FinNp = DimSizesFin[AmOfDimFin - k - 1]; //OC26042019
			if(DimSizesOrig[AmOfDimOrig - k - 1] != FinNp) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
		}
	}

	if(AmOfDimFin == 0)
	{
		AmOfDimFin = 2; DimSizesFin[0] = DimSizesFin[1] = 1;
	}
	else if(AmOfDimFin == 1)
	{
		AmOfDimFin = 2; DimSizesFin[1] = 1;
	}

	float *ftFinData = 0, *ftOrigData = 0;
	double *dtFinData = 0, *dtOrigData = 0;
	srTUtiDataMD::ExtractDataPointer(pIntegResData, ftFinData, dtFinData);
	if((ftFinData == 0) && (dtFinData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
	srTUtiDataMD::ExtractDataPointer(pIntensOrigData, ftOrigData, dtOrigData);
	if((ftOrigData == 0) && (dtOrigData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	for(int j=0; j<DimSizesFin[1]; j++)
	{
		for(int i=0; i<DimSizesFin[0]; i++)
		{
			if(ftFinData != 0) 
			{
				*(ftFinData++) = (float)(ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(ftOrigData + IndStartE, EffNe, StepE));
				ftOrigData += Ne;
			}
			else
			{
				*(dtFinData++) = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(dtOrigData + IndStartE, EffNe, StepE);
				dtOrigData += Ne;
			}
		}
	}
}

//*************************************************************************

void srTRadGenManip::ComponIntegVsHorOrVertPos(srTDataMD* pIntensOrigData, char x_or_z, double xMin, double xMax, srTDataMD* pIntegResData)
{
	if((pIntensOrigData == 0) || (pIntegResData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	double ConstPosInteg = 1000; //assuming data in [*/mm] and position in [m]

	if(pIntensOrigData->DataUnits[0] != '\0')
	{
		if((strcmp(pIntensOrigData->DataUnits, "W/eV/mm^2") == 0) || (strcmp(pIntensOrigData->DataUnits, "W/eV/mm2") == 0))
		{
            strcpy(pIntegResData->DataUnits, "W/eV/mm");
		}
		else if(strcmp(pIntensOrigData->DataUnits, "W/eV/mm") == 0)
		{
            strcpy(pIntegResData->DataUnits, "W/eV");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm^2") == 0) || 
				(strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm^2") == 0) ||
				(strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm2") == 0) || 
				(strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm2") == 0))
		{
            strcpy(pIntegResData->DataUnits, "Ph/s/0.1%bw/mm");
		}
		else if((strcmp(pIntensOrigData->DataUnits, "Ph/s/0.1%bw/mm") == 0) || (strcmp(pIntensOrigData->DataUnits, "Phot/s/0.1%bw/mm") == 0))
		{
            strcpy(pIntegResData->DataUnits, "Ph/s/0.1%bw");
		}
	}

	//long DimSizesOrig[10];
	long long DimSizesOrig[10]; //OC26042019
	int AmOfDimOrig = srTUtiDataMD::ExtractDimSizes(pIntensOrigData, DimSizesOrig);
	if(AmOfDimOrig <= 0) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	//long DimSizesFin[10];
	long long DimSizesFin[10]; //OC26042019
	int AmOfDimFin = srTUtiDataMD::ExtractDimSizes(pIntegResData, DimSizesFin);
	if((AmOfDimFin <= 0) || (AmOfDimFin < (AmOfDimOrig - 1))) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	if((x_or_z == 'x') || (x_or_z == 'X')) x_or_z = 'x';
	else if((x_or_z == 'z') || (x_or_z == 'Z')) x_or_z = 'z';
	else throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
	
	int IndDimToInteg = -1;
	char UnitFirstDim[256];
	strcpy(UnitFirstDim, "m");
	if(strlen(pIntensOrigData->DimUnits[0]) > 0)
	{
		strcpy(UnitFirstDim, pIntensOrigData->DimUnits[0]);
	}

	//long Per1 = 0, Per2 = 0;
	if(x_or_z == 'x') 
	{
		if(AmOfDimOrig == 1)
		{
			IndDimToInteg = 0;
		}
		else if((AmOfDimOrig > 1) && (AmOfDimOrig < 4))
		{
			if(strcmp(UnitFirstDim, "m") == 0) { IndDimToInteg = 0;}
			else { IndDimToInteg = 1;}
		}
		else throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
	}
	else if(x_or_z == 'z')
	{
		if(AmOfDimOrig == 1)
		{
			IndDimToInteg = 0;
		}
		else if(AmOfDimOrig == 2)
		{
			IndDimToInteg = 1;
		}
		else if(AmOfDimOrig == 3)
		{
			IndDimToInteg = 2;
		}
		else throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
	}

	//long Np = pIntensOrigData->DimSizes[IndDimToInteg];
	long long Np = pIntensOrigData->DimSizes[IndDimToInteg]; //OC26042019
	double StartP = srTUtiDataMD::ExtractDimStartValue(pIntensOrigData, IndDimToInteg);
	double StepP = srTUtiDataMD::ExtractDimStep(pIntensOrigData, IndDimToInteg);

	double InvStepP = 0;
	if(StepP != 0) InvStepP = 1./StepP;

	double dIndStartP = (xMin - StartP)*InvStepP;
	//long IndStartP = (long)dIndStartP;
	long long IndStartP = (long long)dIndStartP; //OC26042019
	if((dIndStartP - IndStartP) >= 0.5) IndStartP++; //to improve
	if(IndStartP < 0) IndStartP = 0;
	if(IndStartP >= Np) IndStartP = Np - 1;

	double dIndEndP = (xMax - StartP)*InvStepP;
	//long IndEndP = (long)dIndEndP;
	long long IndEndP = (long long)dIndEndP; //OC26042019
	if((dIndEndP - IndEndP) >= 0.5) IndEndP++;
	if(IndEndP < 0) IndEndP = 0;
	if(IndEndP >= Np) IndEndP = Np - 1;

	//long EffNp = IndEndP - IndStartP;
	long long EffNp = IndEndP - IndStartP; //OC26042019
	if(EffNp < 0) EffNp = 0;
	else if(EffNp > Np) EffNp = Np;

	float *ftFinData = 0, *ftOrigData = 0;
	double *dtFinData = 0, *dtOrigData = 0;
	srTUtiDataMD::ExtractDataPointer(pIntegResData, ftFinData, dtFinData);
	if((ftFinData == 0) && (dtFinData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
	srTUtiDataMD::ExtractDataPointer(pIntensOrigData, ftOrigData, dtOrigData);
	if((ftOrigData == 0) && (dtOrigData == 0)) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	if(IndDimToInteg == 0)
	{
        if(AmOfDimOrig > 1)
		{
			for(int k=0; k<AmOfDimFin; k++)
			{
				if(DimSizesOrig[AmOfDimOrig - k - 1] != DimSizesFin[AmOfDimFin - k - 1]) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
			}
		}
		if(AmOfDimFin == 0)
		{
			AmOfDimFin = 2; DimSizesFin[0] = DimSizesFin[1] = 1;
		}
		else if(AmOfDimFin == 1)
		{
			AmOfDimFin = 2; DimSizesFin[1] = 1;
		}

		for(int j=0; j<DimSizesFin[1]; j++)
		{
			for(int i=0; i<DimSizesFin[0]; i++)
			{
				if(ftFinData != 0) 
				{
					*(ftFinData++) = (float)(ConstPosInteg*CGenMathMeth::Integ1D_FuncDefByArray(ftOrigData + IndStartP, EffNp, StepP));
					ftOrigData += Np;
				}
				else
				{
					*(dtFinData++) = ConstPosInteg*CGenMathMeth::Integ1D_FuncDefByArray(dtOrigData + IndStartP, EffNp, StepP);
					dtOrigData += Np;
				}
			}
		}
	}
	else if(IndDimToInteg == 1)
	{
		if(DimSizesOrig[0] != DimSizesFin[0]) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
		if(AmOfDimOrig > 2)
		{
            if(DimSizesOrig[2] != DimSizesFin[1]) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
		}
		if(AmOfDimFin == 1)
		{
			AmOfDimFin = 2; DimSizesFin[1] = 1;
		}

		//long NperOrig0 = DimSizesOrig[0];
		long long NperOrig0 = DimSizesOrig[0];
		float *fAuxArr = 0;
		double *dAuxArr = 0;
		if(ftFinData != 0) 
		{
			ftFinData += IndStartP*NperOrig0;
			fAuxArr = new float[EffNp];
		}
		else 
		{
			dtFinData += IndStartP*NperOrig0;
			dAuxArr = new double[EffNp];
		}

		for(long long j=0; j<DimSizesFin[1]; j++) //OC26042019
		//for(int j=0; j<DimSizesFin[1]; j++)
		{

			//for(int i=0; i<DimSizesFin[0]; i++)
			//{
			//	if(ftFinData != 0) 
			//	{
			//		//*(ftFinData++) = ConstPosInteg*CGenMathMeth::Integ1D_FuncDefByArray(ftOrigData + IndStartP, EffNp, StepP);
			//		//ftOrigData += Np;
			//	}
			//	else
			//	{
			//		//*(dtFinData++) = ConstPosInteg*CGenMathMeth::Integ1D_FuncDefByArray(dtOrigData + IndStartP, EffNp, StepP);
			//		//dtOrigData += Np;
			//	}
			//}
		}

		if(fAuxArr != 0) delete[] fAuxArr;
		if(dAuxArr != 0) delete[] dAuxArr;

	}
	else if(IndDimToInteg == 2)
	{
		if(DimSizesOrig[0] != DimSizesFin[0]) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;
		if(DimSizesOrig[1] != DimSizesFin[1]) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	}

	//ddddddddddddddddddddddddddddddd
}

//*************************************************************************

void srTRadGenManip::ComponIntegVsHorAndVertPos(srTDataMD* pIntensOrigData, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData)
{
	//ddddddddddddddddddddddddddddddd
}

//*************************************************************************

void srTRadGenManip::ComponIntegVsPhotEnAndPos(srTDataMD* pIntensOrigData, char x_or_z, double zMin, double zMax, srTDataMD* pIntegResData)
{
	//ddddddddddddddddddddddddddddddd
}

//*************************************************************************

void srTRadGenManip::ComponIntegVsPhotEnAndHorAndVertPos(srTDataMD* pIntensOrigData, double eMin, double eMax, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData)
{
	//ddddddddddddddddddddddddddddddd
}

//*************************************************************************
