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
#include "srerror.h"

//DEBUG
//#ifndef __SRIGSEND_H
//#include "srigsend.h"
//#endif
//END DEBUG

//*************************************************************************

void srTRadGenManip::SetupIntCoord(char Cmpn, double Arg, long& i0, long& i1, double& InvStepRelArg)
{
	double Step, Start;
	long N;

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

	long Nz_mi_1 = RadAccessData.nz - 1;
	long Nx_mi_1 = RadAccessData.nx - 1;

    //long iePerE = 0;
    long long iePerE = 0;
	for(long ie=0; ie<RadAccessData.ne; ie++)
	{
        //long izPerZ = 0;
        long long izPerZ = 0;

		double Sum = 0;
		for(long iz=0; iz<RadAccessData.nz; iz++)
		{
			float *pEx_StartForX = pEx0 + izPerZ;
            float *pEz_StartForX = pEz0 + izPerZ;
            //long ixPerX = 0;
            long long ixPerX = 0;

			double SumX = 0.;

            for(long ix=0; ix<RadAccessData.nx; ix++)
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

	long Ne = RadAccessData.ne;
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

	float *pI = 0;
	DOUBLE *pId = 0;
	if(Int_or_ReE != 2) pI = RadExtract.pExtractedData;
	else pId = RadExtract.pExtractedDataD;

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ix0=0, ix1=0, iz0=0, iz1=0;
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

	for(long ie=0; ie<RadAccessData.ne; ie++)
	{
		float *ExPtrs[] = { (pEx_StartForE_ix0_iz0 + iePerE), (pEx_StartForE_ix1_iz0 + iePerE),
							(pEx_StartForE_ix0_iz1 + iePerE), (pEx_StartForE_ix1_iz1 + iePerE)};
		float *EzPtrs[] = { (pEz_StartForE_ix0_iz0 + iePerE), (pEz_StartForE_ix1_iz0 + iePerE),
							(pEz_StartForE_ix0_iz1 + iePerE), (pEz_StartForE_ix1_iz1 + iePerE)};
		double BufI = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
		if(RadExtract.pExtractedData != 0) *(pI++) = (float)BufI;
		if(RadExtract.pExtractedDataD != 0) *(pId++) = (float)BufI;

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

	float *pI = 0;
	DOUBLE *pId = 0;
	if(Int_or_ReE != 2) pI = RadExtract.pExtractedData;
	else pId = RadExtract.pExtractedDataD;

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0, iz0=0, iz1=0;
	double InvStepRelArg1, InvStepRelArg2;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC140813
	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);

	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC140813
	double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC140813

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;

	//long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;
	//long ixPerX = 0;
	long long ixPerX = 0;

	for(long ix=0; ix<RadAccessData.nx; ix++)
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

			for(int ie=0; ie<RadAccessData.ne; ie++)
			{
				*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, PolCom, Int_or_ReE);
				pEx_StAux += 2; pEx_FiAux += 2;
				pEz_StAux += 2; pEz_FiAux += 2;
			}
			resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
		}
		else resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
		//OC150813
		if(pI != 0) *(pI++) = (float)resInt;
		if(pId != 0) *(pId++) = (double)resInt;

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

	float *pI = 0;
	DOUBLE *pId = 0;
	if(Int_or_ReE != 2) pI = RadExtract.pExtractedData;
	else pId = RadExtract.pExtractedDataD;

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0, ix0=0, ix1=0;
	double InvStepRelArg1, InvStepRelArg2;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC150813
	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg2);

	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC150813
	double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1); //OC150813

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	//long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	//long izPerZ = 0;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	long long izPerZ = 0;

	for(long iz=0; iz<RadAccessData.nz; iz++)
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

			for(int ie=0; ie<RadAccessData.ne; ie++)
			{
				*(tInt++) = IntensityComponentSimpleInterpol(pEx_StAux, pEx_FiAux, pEz_StAux, pEz_FiAux, InvStepRelArg2, PolCom, Int_or_ReE);
				pEx_StAux += 2; pEx_FiAux += 2;
				pEz_StAux += 2; pEz_FiAux += 2;
			}
			resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
		}
		else resInt = IntensityComponentSimpleInterpol2D(ExPtrs, EzPtrs, InvStepRelArg1, InvStepRelArg2, PolCom, Int_or_ReE);
		//OC150813
		if(pI != 0) *(pI++) = (float)resInt;
		if(pId != 0) *(pId++) = (double)resInt;

		izPerZ += PerZ;
	}
	if(arAuxInt != 0) delete[] arAuxInt; //OC150813
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ(srTRadExtract& RadExtract)
{
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;

	if(Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pI = 0;
	DOUBLE *pId = 0;
	if(Int_or_ReE != 2) pI = RadExtract.pExtractedData;
	else pId = RadExtract.pExtractedDataD;

	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0;
	double InvStepRelArg=0;
	//SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg); //OC140813
	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1); //OC140813
	double *arAuxInt = 0, resInt;
	if(intOverEnIsRequired)
	{
		arAuxInt = new double[RadAccessData.ne];
	}
	else SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg); //OC140813

	//double ConstPhotEnInteg = 1.60219e-16; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV
	//if(RadAccessData.ElecFldUnit != 1) ConstPhotEnInteg = 1; //(?) 
	double ConstPhotEnInteg = 1.; //1 Phot/s/.1%bw correspond(s) to : 1.60219e-16 W/eV

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	//long izPerZ = 0;
	long long izPerZ = 0;

	for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;

		float *pEx_St = pEx_StartForX + Two_ie0;
		float *pEz_St = pEz_StartForX + Two_ie0;
		float *pEx_Fi = pEx_StartForX + Two_ie1;
		float *pEz_Fi = pEz_StartForX + Two_ie1;

		for(long ix=0; ix<RadAccessData.nx; ix++)
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
				for(int ie=0; ie<RadAccessData.ne; ie++)
				{
					*(tInt++) = IntensityComponent(pEx_StAux, pEz_StAux, PolCom, Int_or_ReE);
					pEx_StAux += 2;
					pEz_StAux += 2;
				}
				resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
			}
			else resInt = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
			//OC140813
			if(pI != 0) *(pI++) = (float)resInt;
			if(pId != 0) *(pId++) = (double)resInt;

			//ixPerX += PerX;
			pEx_St += PerX;
			pEz_St += PerX;
			pEx_Fi += PerX;
			pEz_Fi += PerX;
		}
		izPerZ += PerZ;
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

	float *pI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long iz0=0, iz1=0;
	double InvStepRelArg;
	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg);

	//long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;

	//long ixPerX = 0;
	long long ixPerX = 0;

	for(long ix=0; ix<RadAccessData.nx; ix++)
	{
		float *pEx_StartForE_iz0 = pEx_StartForX_iz0 + ixPerX, *pEx_StartForE_iz1 = pEx_StartForX_iz1 + ixPerX;
		float *pEz_StartForE_iz0 = pEz_StartForX_iz0 + ixPerX, *pEz_StartForE_iz1 = pEz_StartForX_iz1 + ixPerX;
		long iePerE = 0;

		for(long ie=0; ie<RadAccessData.ne; ie++)
		{
			float *pEx_St = pEx_StartForE_iz0 + iePerE, *pEx_Fi = pEx_StartForE_iz1 + iePerE;
			float *pEz_St = pEz_StartForE_iz0 + iePerE, *pEz_Fi = pEz_StartForE_iz1 + iePerE;

			*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);

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

	float *pI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ix0=0, ix1=0;
	double InvStepRelArg;
	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg);

	//long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	//long izPerZ = 0;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;
	long long izPerZ = 0;

	for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;

		float *pEx_StartForE_ix0 = pEx_StartForX + ix0PerX, *pEx_StartForE_ix1 = pEx_StartForX + ix1PerX;
		float *pEz_StartForE_ix0 = pEz_StartForX + ix0PerX, *pEz_StartForE_ix1 = pEz_StartForX + ix1PerX;
		long iePerE = 0;

		for(long ie=0; ie<RadAccessData.ne; ie++)
		{
			float *pEx_St = pEx_StartForE_ix0 + iePerE, *pEx_Fi = pEx_StartForE_ix1 + iePerE;
			float *pEz_St = pEz_StartForE_ix0 + iePerE, *pEz_Fi = pEz_StartForE_ix1 + iePerE;

			*(pI++) = IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
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

	float *pI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	//long PerX = RadAccessData.ne << 1;
	//long PerZ = PerX*RadAccessData.nx;
	//long izPerZ = 0;
	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;
	long long izPerZ = 0;

	for(long iz=0; iz<RadAccessData.nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;
		//long ixPerX = 0;
		long long ixPerX = 0;

		for(long ix=0; ix<RadAccessData.nx; ix++)
		{
			float *pEx_StartForE = pEx_StartForX + ixPerX;
			float *pEz_StartForE = pEz_StartForX + ixPerX;
			long iePerE = 0;

			for(long ie=0; ie<RadAccessData.ne; ie++)
			{
				float* pEx = pEx_StartForE + iePerE;
				float* pEz = pEz_StartForE + iePerE;

				*(pI++) = IntensityComponent(pEx, pEz, PolCom, Int_or_ReE);

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
{//OC06092018
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0, iz0=0, iz1=0;
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('z', RadExtract.z, iz0, iz1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool NeedInterp = true;
	if(((iz0 == iz1) || (fabs(InvStepRelArg2) < tolRelArg)) && 
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) NeedInterp = false;

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long iz0PerZ = iz0*PerZ, iz1PerZ = iz1*PerZ;
	float *pEx_StartForX_iz0 = pEx0 + iz0PerZ, *pEx_StartForX_iz1 = pEx0 + iz1PerZ;
	float *pEz_StartForX_iz0 = pEz0 + iz0PerZ, *pEz_StartForX_iz1 = pEz0 + iz1PerZ;
	float *pEx_StartForE_iz0, *pEx_StartForE_iz1;
	float *pEz_StartForE_iz0, *pEz_StartForE_iz1;

	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	long nx = RadAccessData.nx;
	long long ixPerX = 0;
	for(long ix=0; ix<nx; ix++)
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
	long nx_mi_1 = nx - 1;
	for(long ixt=0; ixt<nx_mi_1; ixt++)
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

		for(long ix=ixt+1; ix<nx; ix++)
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
				if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, pMI)) return res;
			}
			else
			{
				float *pEx = pEx_StartForE_iz0 + Two_ie0, *pExT = pEx_StartForE_iz0_t + Two_ie0;
				float *pEz = pEz_StartForE_iz0 + Two_ie0, *pEzT = pEz_StartForE_iz0_t + Two_ie0;
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, pMI)) return res;
			}
			pMI += 2;
		}
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecMutualIntensityVsZ(srTRadExtract& RadExtract)
{//OC10092018
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0, ix0=0, ix1=0;
	double InvStepRelArg1, InvStepRelArg2;

	SetupIntCoord('x', RadExtract.x, ix0, ix1, InvStepRelArg2);
	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg1);

	const double tolRelArg = 1.e-08;
	bool NeedInterp = true;
	if(((ix0 == ix1) || (fabs(InvStepRelArg2) < tolRelArg)) && 
	   ((ie0 == ie1) || (fabs(InvStepRelArg1) < tolRelArg))) NeedInterp = false;

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long long ix0PerX = ix0*PerX, ix1PerX = ix1*PerX;

	float *pEx_StartForZ_ix0 = pEx0 + ix0PerX, *pEx_StartForZ_ix1 = pEx0 + ix1PerX;
	float *pEz_StartForZ_ix0 = pEz0 + ix0PerX, *pEz_StartForZ_ix1 = pEz0 + ix1PerX;
	float *pEx_StartForE_ix0, *pEx_StartForE_ix1;
	float *pEz_StartForE_ix0, *pEz_StartForE_ix1;

	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	long nz = RadAccessData.nz;
	long long izPerZ = 0;
	for(long iz=0; iz<nz; iz++)
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
	long nz_mi_1 = nz - 1;
	for(long izt=0; izt<nz_mi_1; izt++)
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

		for(long iz=izt+1; iz<nz; iz++)
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
				if(res = MutualIntensityComponentSimpleInterpol2D(ExPtrs, ExPtrsT, EzPtrs, EzPtrsT, InvStepRelArg1, InvStepRelArg2, PolCom, pMI)) return res;
			}
			else
			{
				float *pEx = pEx_StartForE_ix0 + Two_ie0, *pExT = pEx_StartForE_ix0_t + Two_ie0;
				float *pEz = pEz_StartForE_ix0 + Two_ie0, *pEzT = pEz_StartForE_ix0_t + Two_ie0;
				if(res = MutualIntensityComponent(pEx, pExT, pEz, pEzT, PolCom, pMI)) return res;
			}
			pMI += 2;
		}
	}
	return 0;
}

//*************************************************************************

int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ(srTRadExtract& RadExtract)
{//OC11092018
	int res = 0;
	int PolCom = RadExtract.PolarizCompon;
	int Int_or_ReE = RadExtract.Int_or_Phase;
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	float *pMI = RadExtract.pExtractedData;
	float *pEx0 = RadAccessData.pBaseRadX;
	float *pEz0 = RadAccessData.pBaseRadZ;

	long long PerX = RadAccessData.ne << 1;
	long long PerZ = PerX*RadAccessData.nx;

	long ie0=0, ie1=0;
	double InvStepRelArg=0;

	SetupIntCoord('e', RadExtract.ePh, ie0, ie1, InvStepRelArg);

	const double tolRelArg = 1.e-08;
	bool NeedInterp = true;
	if((ie0 == ie1) || (fabs(InvStepRelArg) < tolRelArg)) NeedInterp = false;

	long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1;
	long nz = RadAccessData.nz, nx = RadAccessData.nx;

	//First RadAccessData.nx values are the "normal" intensity (real diagonal elements of mutual intensity)
	long long izPerZ = 0;
	for(long iz=0; iz<nz; iz++)
	{
		float *pEx_StartForX = pEx0 + izPerZ;
		float *pEz_StartForX = pEz0 + izPerZ;

		float *pEx_St = pEx_StartForX + Two_ie0;
		float *pEz_St = pEz_StartForX + Two_ie0;
		float *pEx_Fi = pEx_StartForX + Two_ie1;
		float *pEz_Fi = pEz_StartForX + Two_ie1;

		for(long ix=0; ix<nx; ix++)
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
	return 0;
}

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

	CGenMathFFT2D FFT2D;
	FFT2D.NextCorrectNumberForFFT(Nx);
	FFT2D.NextCorrectNumberForFFT(Nz);

	//long TotAmOfNewData = (Nx*Nz) << 1;
	long long TotAmOfNewData = (((long long)Nx)*((long long)Nz)) << 1;

	pTempStorage = new float[TotAmOfNewData];
	if(pTempStorage == 0) return MEMORY_ALLOCATION_FAILURE;

	OwnRadExtract.pExtractedData = pTempStorage;
	OwnRadExtract.PlotType = 3; // vs x&z

	long Ne = RadAccessData.ne;
	OwnRadExtract.ePh = RadAccessData.eStart;

	char SinglePhotonEnergy = ((PT == 1) || (PT == 2) || (PT == 3));
	if(SinglePhotonEnergy) 
	{
		Ne = 1; OwnRadExtract.ePh = RadExtract.ePh;
	}

	for(long ie=0; ie<Ne; ie++)
	{
		if(result = ExtractSingleElecIntensity2DvsXZ(OwnRadExtract)) return result;
		if(result = ConvoluteWithElecBeamOverTransvCoord(OwnRadExtract.pExtractedData, Nx, Nz)) return result;

		long ie0 = ie;
		if(SinglePhotonEnergy && (RadAccessData.ne > 1))
		{
			long ie1;
			double RelArgE_Dummy;
			SetupIntCoord('e', OwnRadExtract.ePh, ie0, ie1, RelArgE_Dummy);
		}

		PutConstPhotEnergySliceInExtractPlace(ie0, Nx, Nz, OwnRadExtract, RadExtract);
		OwnRadExtract.ePh += RadAccessData.eStep;
	}

	if(pTempStorage != 0) delete[] pTempStorage;
	return 0;
}

//*************************************************************************

int srTRadGenManip::ConvoluteWithElecBeamOverTransvCoord(float* DataToConv, long Nx, long Nz)
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

void srTRadGenManip::PadImZerosToRealData(float* pData, long Nx, long Nz)
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

		for(long iz=RadAccessData.nz; iz<Nz; iz++)
		{
			for(long ix=0; ix<Nx; ix++)
			{
				*(pTot++) = 0.; *(pTot++) = 0.;
			}
		}
	}
	if(Nx > RadAccessData.nx)
	{
		long TwoNxOld = RadAccessData.nx << 1;
		long NzOld_mi_1 = RadAccessData.nz - 1;
		long ShiftPer = (Nx - RadAccessData.nx) << 1;

		//float* pOrig_Start = pData + (TwoNxOld*NzOld_mi_1);
		float* pOrig_Start = pData + (((long long)TwoNxOld)*((long long)NzOld_mi_1));
		//long ShiftLen = ShiftPer*NzOld_mi_1;
		long long ShiftLen = ((long long)ShiftPer)*((long long)NzOld_mi_1);

		for(long iz=0; iz<NzOld_mi_1; iz++)
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
	DOUBLE *MatrStrPtrs[4];
	for(i=0; i<4; i++)
	{
		int i4 = i << 2;
		MatrStrPtrs[i] = RadAccessData.p4x4PropMatr + i4;
	}

	DOUBLE *Vect = RadAccessData.p4x4PropMatr + 16;

// First-Order Moments
	DOUBLE InitFirstOrderMom[] = { ElecBeamMom.Mx, ElecBeamMom.Mxp, ElecBeamMom.Mz, ElecBeamMom.Mzp};
	DOUBLE FinFirstOrderMom[4];
	for(i=0; i<4; i++)
	{
		DOUBLE Res_i = 0.;
		DOUBLE *MatrStrPtrs_i = MatrStrPtrs[i];
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
	DOUBLE* pStr = *MatrStrPtrs;
	DOUBLE a00 = pStr[0], a01 = pStr[1], a02 = pStr[2], a03 = pStr[3];
	pStr = MatrStrPtrs[1];
	DOUBLE a10 = pStr[0], a11 = pStr[1], a12 = pStr[2], a13 = pStr[3];
	pStr = MatrStrPtrs[2];
	DOUBLE a20 = pStr[0], a21 = pStr[1], a22 = pStr[2], a23 = pStr[3];
	pStr = MatrStrPtrs[3];
	DOUBLE a30 = pStr[0], a31 = pStr[1], a32 = pStr[2], a33 = pStr[3];
	DOUBLE Str0[] = { a00*a00, 2*a00*a01, a01*a01, a02*a02, 2*a02*a03, a03*a03, 2*a00*a02, 2*a01*a02, 2*a00*a03, 2*a01*a03};
	DOUBLE Str1[] = { a00*a10, a01*a10 + a00*a11, a01*a11, a02*a12, a03*a12 + a02*a13, a03*a13, a02*a10 + a00*a12, a02*a11 + a01*a12, a03*a10 + a00*a13, a03*a11 + a01*a13};
	DOUBLE Str2[] = { a10*a10, 2*a10*a11, a11*a11, a12*a12, 2*a12*a13, a13*a13, 2*a10*a12, 2*a11*a12, 2*a10*a13, 2*a11*a13};
	DOUBLE Str3[] = { a20*a20, 2*a20*a21, a21*a21, a22*a22, 2*a22*a23, a23*a23, 2*a20*a22, 2*a21*a22, 2*a20*a23, 2*a21*a23};
	DOUBLE Str4[] = { a20*a30, a21*a30 + a20*a31, 21*a31, a22*a32, a23*a32 + a22*a33, a23*a33, a22*a30 + a20*a32, a22*a31 + a21*a32, a23*a30 + a20*a33, a23*a31 + a21*a33};
	DOUBLE Str5[] = { a30*a30, 2*a30*a31, a31*a31, a32*a32, 2*a32*a33, a33*a33, 2*a30*a32, 2*a31*a32, 2*a30*a33, 2*a31*a33};
	DOUBLE Str6[] = { a00*a20, a01*a20 + a00*a21, a01*a21, a02*a22, a03*a22 + a02*a23, a03*a23, a02*a20 + a00*a22, a02*a21 + a01*a22, a03*a20 + a00*a23, a03*a21 + a01*a23};
	DOUBLE Str7[] = { a10*a20, a11*a20 + a10*a21, a11*a21, a12*a22, a13*a22 + a12*a23, a13*a23, a12*a20 + a10*a22, a12*a21 + a11*a22, a13*a20 + a10*a23, a13*a21 + a11*a23};
	DOUBLE Str8[] = { a00*a30, a01*a30 + a00*a31, a01*a31, a02*a32, a03*a32 + a02*a33, a03*a33, a02*a30 + a00*a32, a02*a31 + a01*a32, a03*a30 + a00*a33, a03*a31 + a01*a33};
	DOUBLE Str9[] = { a10*a30, a11*a30 + a10*a31, a11*a31, a12*a32, a13*a32 + a12*a33, a13*a33, a12*a30 + a10*a32, a12*a31 + a11*a32, a13*a30 + a10*a33, a13*a31 + a11*a33};
	DOUBLE *Matr10x10[10];
	Matr10x10[0] = (DOUBLE*)&Str0; 
	Matr10x10[1] = (DOUBLE*)&Str1; 
	Matr10x10[2] = (DOUBLE*)&Str2; 
	Matr10x10[3] = (DOUBLE*)&Str3; 
	Matr10x10[4] = (DOUBLE*)&Str4; 
	Matr10x10[5] = (DOUBLE*)&Str5; 
	Matr10x10[6] = (DOUBLE*)&Str6; 
	Matr10x10[7] = (DOUBLE*)&Str7; 
	Matr10x10[8] = (DOUBLE*)&Str8; 
	Matr10x10[9] = (DOUBLE*)&Str9; 
	DOUBLE InitSecondOrderMom[] = { ElecBeamMom.Mxx, ElecBeamMom.Mxxp, ElecBeamMom.Mxpxp, ElecBeamMom.Mzz, ElecBeamMom.Mzzp, ElecBeamMom.Mzpzp, ElecBeamMom.Mxz, ElecBeamMom.Mxpz, ElecBeamMom.Mxzp, ElecBeamMom.Mxpzp};
	DOUBLE FinSecondOrderMom[10];
	for(i=0; i<10; i++)
	{
		DOUBLE Res_i = 0.;
		DOUBLE *MatrStr_i = Matr10x10[i];
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

void srTRadGenManip::PutConstPhotEnergySliceInExtractPlace(long ie, long NxSlice, long NzSlice, srTRadExtract& LocRadExtract, srTRadExtract& RadExtract)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

	//float *pGen0 = RadExtract.pExtractedData, *pGenStart;
	float *pGen0 = RadExtract.pExtractedData;
	float *pGenStart = pGen0; //OC160815
	float *pSlice0 = LocRadExtract.pExtractedData;

	long Nx, Nz;
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

	long ix0_Slice, ix1_Slice, iz0_Slice, iz1_Slice;
	//long izPerZ_Gen = 0;
	long long izPerZ_Gen = 0;

	long Nz_mi_1 = Nz - 1;
	long Nx_mi_1 = Nx - 1;
	float wx, wz; //, wtot;

	for(long iz=0; iz<Nz; iz++)
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
	long Nx = WaveData.DimSizes[0];
	long Nz = WaveData.DimSizes[1];

	double* CenterSlice = new double[Nx];
	if(CenterSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	//long ixMid = Nx >> 1, izMid = Nz >> 1;
	long izMid = Nz >> 1;

	long ix, iz;
	DOUBLE *pData = (DOUBLE*)(WaveData.pWaveData);
	DOUBLE *tm = pData + izMid*Nx;
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

void srTRadGenManip::TryToMakePhaseContinuous1D(double* Slice, long Np, long i0, float Phi0)
{
	const double TwoPi = 6.2831853071796;
	const double cFlip = TwoPi - 2.5; // To steer

	float PhToAdd0 = (float)((i0 != -1)? (Phi0 - Slice[i0]) : 0.);

	long HalfNp = Np >> 1;
	long OtherHalfNp = Np - HalfNp;

	double PhToAdd = PhToAdd0;
	double *t = Slice + HalfNp - 1; 
	*t += PhToAdd;
	double PrevPh = *(t--);
	for(long i=0; i<(HalfNp - 1); i++)
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

void srTRadGenManip::ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData)
{
	if(pData == 0) throw INCORRECT_PARAMS_WFR_COMPON_EXTRACT;

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
	srTGenOptElem GenOptElem;
	srTRadExtract RadExtract(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);

			//DEBUG
			//srTIgorSend::WarningMessage("ExtractRadiation #1");
			//END DEBUG

	int res;
	if(TransvPres != RadAccessData.Pres)
		if(res = GenOptElem.SetRadRepres(&RadAccessData, char(TransvPres))) throw res;
	
	if(RadExtract.Int_or_Phase == 1)
	{//1- Multi-Elec Intensity
		if(res = ComputeConvolutedIntensity(RadExtract)) throw res;
	}
	else if(RadExtract.Int_or_Phase == 4) //OC
	{//4- Sigle-Elec. Flux
		if(res = ExtractFluxFromWfr(RadExtract, 's')) throw res;
	}
	else if(RadExtract.Int_or_Phase == 5) //OC
	{//5- Multi-Elec. Flux
		if(res = ExtractFluxFromWfr(RadExtract, 'm')) throw res;
	}
	else if(RadExtract.Int_or_Phase == 8) //OC06092018
	{
		if(res = ExtractSingleElecMutualIntensity(RadExtract)) throw res;
	}
	else
	{
		if(res = ExtractSingleElecIntensity(RadExtract)) throw res;
	}

			//DEBUG
			//srTIgorSend::WarningMessage("ExtractRadiation: ready out");
			//END DEBUG

	//if(res = SetupExtractedWaveData(RadExtract, ExtractedWaveData)) return res;
}

//*************************************************************************

void srTRadGenManip::IntProc(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar)
{//Added to move some simple "loops" from slow Py to C++
	if((pwI1 == 0) || (pwI2 == 0) || (arPar == 0)) throw SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	char typeProc = (char)(arPar[0]);

	if(typeProc == 1)
	{//Sum-up intensities on same mesh
	
	}
	else if(typeProc == 2)
	{//Add I2 to I1 with interpolation to the mesh of I1
	
	}
	else if(typeProc == 3) srTRadGenManip::Int2DIntegOverAzim(pwI1, pwI2, arPar + 1);
}

//*************************************************************************

void srTRadGenManip::Int2DIntegOverAzim(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar)
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
	long nx = *(pwI2->DimSizes);
	
	double yMin = *((pwI2->DimStartValues) + 1);
	double yStep = *((pwI2->DimSteps) + 1);
	long ny = *((pwI2->DimSizes) + 1);

	double *arFuncForAzimInt=0;
	if(nPhMax > 0) arFuncForAzimInt = new double[nPhMax];

	long nr = *(pwI1->DimSizes);
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

		tolPhMeth2 = ::fabs(phRange*relPrec);
		if(::fabs(::fabs(phRange) - twoPi) < tolPhMeth2) tolPhMeth2 = 0; //meaning that calculation of edge derivatives vs phi is not necessary
	}

	for(int ir=0; ir<nr; ir++)
	{
		double curIntOverPh = 0.;
		if(meth == 1)
		{
			long nPhCur = (long)round((r/rMax)*nPhMax);
			if(nPhCur <= 1)
			{
				if(calcAvg) 
				{
					if(pfI2) curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
					else curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
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
						*(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
						ph += phStep;
					}
					if(::fabs(::fabs(phRange) - twoPi) < fabs(1.e-03*phStep)) *t_ar = *arFuncForAzimInt;
					else *t_ar = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
				}
				else
				{
					for(int iph=0; iph<nPhCur_mi_1; iph++)
					{
						//x = r*cos(ph); y = r*sin(ph);
						*(t_ar++) = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
						ph += phStep;
					}
					if(::fabs(::fabs(phRange) - twoPi) < ::fabs(1.e-03*phStep)) *t_ar = *arFuncForAzimInt;
					else *t_ar = CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
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
					if(pfI2) curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pfI2, ordInterp);
					else curIntOverPh = CGenMathInterp::InterpOnRegMesh2d(x0 + r, y0, xMin, xStep, nx, yMin, yStep, ny, pdI2, ordInterp);
				}
			}
			else
			{
				auxStruct.r = r;
				if(tolPhMeth2 > 0) curIntOverPh = CGenMathMeth::Integ1D_Func(&(srTRadGenManip::IntCylCrd), phMin, phMax, relPrec, (void*)(&auxStruct));
				else curIntOverPh = CGenMathMeth::Integ1D_FuncWithEdgeDer(&(srTRadGenManip::IntCylCrd), phMin, phMax, dFdPh1, dFdPh2, relPrec, (void*)(&auxStruct));
			}
		}
		
		if(calcAvg) curIntOverPh /= phRange;

		if(pfI1) *(tfI1++) = (float)curIntOverPh;
		else *(tdI1++) = curIntOverPh;

		r += rStep;
	}

	if(arFuncForAzimInt != 0) delete[] arFuncForAzimInt;
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

	long DimSizesOrig[10];
	int AmOfDimOrig = srTUtiDataMD::ExtractDimSizes(pIntensOrigData, DimSizesOrig);
	if(AmOfDimOrig <= 0) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	long DimSizesFin[10];
	int AmOfDimFin = srTUtiDataMD::ExtractDimSizes(pIntegResData, DimSizesFin);
	if((AmOfDimFin <= 0) || (AmOfDimFin < (AmOfDimOrig - 1))) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	long Ne = pIntensOrigData->DimSizes[0];
	double StartE = srTUtiDataMD::ExtractDimStartValue(pIntensOrigData, 0);
	double StepE = srTUtiDataMD::ExtractDimStep(pIntensOrigData, 0);

	double InvStepE = 0;
	if(StepE != 0) InvStepE = 1./StepE;

	double dIndStartE = (eMin - StartE)*InvStepE;
	long IndStartE = (long)dIndStartE;
	if((dIndStartE - IndStartE) >= 0.5) IndStartE++; //to improve
	if(IndStartE < 0) IndStartE = 0;
	if(IndStartE >= Ne) IndStartE = Ne - 1;

	double dIndEndE = (eMax - StartE)*InvStepE;
	long IndEndE = (long)dIndEndE;
	if((dIndEndE - IndEndE) >= 0.5) IndEndE++; //to improve
	if(IndEndE < 0) IndEndE = 0;
	if(IndEndE >= Ne) IndEndE = Ne - 1;

	long EffNe = IndEndE - IndStartE;
	if(EffNe < 0) EffNe = 0;
	else if(EffNe > Ne) EffNe = Ne;

	if(AmOfDimOrig > 1)
	{
		for(int k=0; k<AmOfDimFin; k++)
		{
			int FinNp = DimSizesFin[AmOfDimFin - k - 1];
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

	long DimSizesOrig[10];
	int AmOfDimOrig = srTUtiDataMD::ExtractDimSizes(pIntensOrigData, DimSizesOrig);
	if(AmOfDimOrig <= 0) throw INCORRECT_PARAMS_WFR_COMPON_INTEG;

	long DimSizesFin[10];
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

	long Np = pIntensOrigData->DimSizes[IndDimToInteg];
	double StartP = srTUtiDataMD::ExtractDimStartValue(pIntensOrigData, IndDimToInteg);
	double StepP = srTUtiDataMD::ExtractDimStep(pIntensOrigData, IndDimToInteg);

	double InvStepP = 0;
	if(StepP != 0) InvStepP = 1./StepP;

	double dIndStartP = (xMin - StartP)*InvStepP;
	long IndStartP = (long)dIndStartP;
	if((dIndStartP - IndStartP) >= 0.5) IndStartP++; //to improve
	if(IndStartP < 0) IndStartP = 0;
	if(IndStartP >= Np) IndStartP = Np - 1;

	double dIndEndP = (xMax - StartP)*InvStepP;
	long IndEndP = (long)dIndEndP;
	if((dIndEndP - IndEndP) >= 0.5) IndEndP++;
	if(IndEndP < 0) IndEndP = 0;
	if(IndEndP >= Np) IndEndP = Np - 1;

	long EffNp = IndEndP - IndStartP;
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

		long NperOrig0 = DimSizesOrig[0];
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

		for(int j=0; j<DimSizesFin[1]; j++)
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
