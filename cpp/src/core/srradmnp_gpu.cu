/************************************************************************//**
 * File: srradmnp_gpu.cu
 * Description: Various "manipulations" with Radiation data (e.g. "extraction" of Intensity from Electric Field, etc.) (CUDA implementation)
 * Project: Synchrotron Radiation Workshop
 * First release: 2023
 *
 * Copyright (C) Brookhaven National Laboratory
 * All Rights Reserved
 *
 * @author H.Goel
 * @version 1.0
 ***************************************************************************/

#ifdef _OFFLOAD_GPU
#include <stdio.h>
#include <stdlib.h>
#include <array>
#include <string>
#include <assert.h>
#include <math.h>
#include "srradmnp.h"
#include "gmmeth.h"

template<int PolCom, bool NpIsEven>
__device__ double Integ_Intensity(srTRadGenManip *obj, float* pEx, float* pEz, int Int_or_ReE, int ne, double eStep) //HG31072024
{
	double s0 = 0., s1 = 0., s2 = 0., s3 = 0.;

	if (ne == 2)
	{
		s0 = obj->IntensityComponent(pEx, pEz, PolCom, Int_or_ReE);
		s3 = obj->IntensityComponent(pEx + 2, pEz + 2, PolCom, Int_or_ReE);
		return (s0 + s3) * 0.5 * eStep;
	}

	long long NpSim = ne;
	if (NpIsEven) NpSim--;

	float* tpEx = pEx;
	float* tpEz = pEz;

	s0 = obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
	tpEx += 2; tpEz += 2;
	for (long long i = 1; i < ((NpSim - 3) >> 1); i++)
	{
		s1 += obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
		tpEx += 2; tpEz += 2;
		s2 += obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
		tpEx += 2; tpEz += 2;
	}
	s1 += obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
	tpEx += 2; tpEz += 2;

	s3 = obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
	tpEx += 2; tpEz += 2;

	double res = (eStep/3.)*(s0 + 4.*s1 + 2.*s2 + s3);

	if (!NpIsEven)
	{
		double s4 = obj->IntensityComponent(tpEx, tpEz, PolCom, Int_or_ReE);
		
		res += (double)(0.5 * eStep * (s3 + s4));
	}
	return res;
}

//__global__ void ExtractSingleElecIntensity2DvsXZ_Kernel(srTRadExtract RadExtract, srTSRWRadStructAccessData RadAccessData, srTRadGenManip *obj, double* arAuxInt, long long ie0, long long ie1, double InvStepRelArg, int Int_or_ReE)
template <bool allStokesReq, bool intOverEnIsRequired, int PolCom, bool NpIsEven>
__global__ void ExtractSingleElecIntensity2DvsXZ_Kernel(srTRadExtract RadExtract, srTSRWRadStructAccessData* pRadAccessData, srTRadGenManip *obj, long long ie0, long long ie1, double InvStepRelArg, int Int_or_ReE) //HG31072024 Redesigned to handle ne > 1, still needs debugging
{
	int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range
    
	if (ix < pRadAccessData->nx && iz < pRadAccessData->nz) 
    {
		//int PolCom = RadExtract.PolarizCompon;
			
		//bool allStokesReq = (PolCom == -5); //OC18042020

		float* pI = 0, * pI1 = 0, * pI2 = 0, * pI3 = 0; //OC17042020
		double* pId = 0, * pI1d = 0, * pI2d = 0, * pI3d = 0;
		long ne = pRadAccessData->ne, nx = pRadAccessData->nx, nz = pRadAccessData->nz;
		//float *pI = 0;
		//DOUBLE *pId = 0;
		//double *pId = 0; //OC26112019 (related to SRW port to IGOR XOP8 on Mac)
		long long nxnz = ((long long)nx) * ((long long)nz);
		if (Int_or_ReE != 2)
		{
			pI = RadExtract.pExtractedData;
			if (allStokesReq) //OC17042020
			{
				pI1 = pI + nxnz; pI2 = pI1 + nxnz; pI3 = pI2 + nxnz;
			}
		}
		else
		{
			pId = RadExtract.pExtractedDataD;
			if (allStokesReq) //OC17042020
			{
				pI1d = pId + nxnz; pI2d = pI1d + nxnz; pI3d = pI2d + nxnz;
			}
		}

		float* pEx0 = pRadAccessData->pBaseRadX;
		float* pEz0 = pRadAccessData->pBaseRadZ;

		//long PerX = pRadAccessData->ne << 1;
		//long PerZ = PerX*pRadAccessData->nx;
		//long long PerX = pRadAccessData->ne << 1;
		//long long PerZ = PerX*pRadAccessData->nx;
		long long PerX = ((long long)ne) << 1; //OC18042020
		long long PerZ = PerX * nx;

		//bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (ne > 1); //OC18042020
		double resInt, resInt1, resInt2, resInt3;
		double ConstPhotEnInteg = 1.;
		long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
		
		long offset = iz * PerZ + ix * PerX;
		long offsetExIntens = offset / PerX;

		float* pEx_StartForX = pEx0 + offset;
		float* pEz_StartForX = pEz0 + offset;
		if (pI != 0)
		{
			pI += offsetExIntens;
			if (allStokesReq)
			{
				pI1 += offsetExIntens;
				pI2 += offsetExIntens;
				pI3 += offsetExIntens;
			}
		} 

		if (pId != 0)
		{
			pId += offsetExIntens;
			if (allStokesReq)
			{
				pI1d += offsetExIntens;
				pI2d += offsetExIntens;
				pI3d += offsetExIntens;
			}
		} 
		
		//long ixPerX = 0;

		float* pEx_St = pEx_StartForX + Two_ie0;
		float* pEz_St = pEz_StartForX + Two_ie0;
		float* pEx_Fi = pEx_StartForX + Two_ie1;
		float* pEz_Fi = pEz_StartForX + Two_ie1;

		if (intOverEnIsRequired) //OC140813
		{//integrate over photon energy / time
			//float* pEx_StAux = pEx_St;
			//float* pEz_StAux = pEz_St;

			if (!allStokesReq) //OC17042020
			{
				
				//for (ie = 0; ie < ne; ie++) //OC18042020
				//for(int ie=0; ie<RadAccessData.ne; ie++)
				//{
				//	*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, PolCom, Int_or_ReE);
				//	pEx_StAux += 2;
				//	pEz_StAux += 2;
				//}
				//resInt = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep); //OC18042020
				//resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
				resInt = ConstPhotEnInteg * Integ_Intensity<PolCom, NpIsEven>(obj, pEx_St, pEz_St, Int_or_ReE, ne, pRadAccessData->eStep); //HG31072024
			}
			else
			{
				//for (ie = 0; ie < ne; ie++)
				//{
				//	*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -1, Int_or_ReE);
				//	pEx_StAux += 2; pEz_StAux += 2;
				//}
				//resInt = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				//tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				//for (ie = 0; ie < ne; ie++)
				//{
				//	*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -2, Int_or_ReE);
				//	pEx_StAux += 2; pEz_StAux += 2;
				//}
				//resInt1 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				//tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				//for (ie = 0; ie < ne; ie++)
				//{
				//	*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -3, Int_or_ReE);
				//	pEx_StAux += 2; pEz_StAux += 2;
				//}
				//resInt2 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				//tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				//for (ie = 0; ie < ne; ie++)
				//{
				//	*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -4, Int_or_ReE);
				//	pEx_StAux += 2; pEz_StAux += 2;
				//}
				//resInt3 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);
				resInt = ConstPhotEnInteg * Integ_Intensity<-1, NpIsEven>(obj, pEx_St, pEz_St, Int_or_ReE, ne, pRadAccessData->eStep); //HG31072024
				resInt1 = ConstPhotEnInteg * Integ_Intensity<-2, NpIsEven>(obj, pEx_St, pEz_St, Int_or_ReE, ne, pRadAccessData->eStep);
				resInt2 = ConstPhotEnInteg * Integ_Intensity<-3, NpIsEven>(obj, pEx_St, pEz_St, Int_or_ReE, ne, pRadAccessData->eStep);
				resInt3 = ConstPhotEnInteg * Integ_Intensity<-4, NpIsEven>(obj, pEx_St, pEz_St, Int_or_ReE, ne, pRadAccessData->eStep);
			}
		}
		else
		{
			if (!allStokesReq) //OC18042020
			{
				resInt = obj->IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, PolCom, Int_or_ReE);
			}
			else //OC18042020
			{
				resInt = obj->IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -1, Int_or_ReE);
				resInt1 = obj->IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -2, Int_or_ReE);
				resInt2 = obj->IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -3, Int_or_ReE);
				resInt3 = obj->IntensityComponentSimpleInterpol(pEx_St, pEx_Fi, pEz_St, pEz_Fi, InvStepRelArg, -4, Int_or_ReE);
			}
		}
		//OC140813
		if (pI != 0) *pI = (float)resInt;
		if (pId != 0) *pId = resInt; //OC18042020
		//if(pId != 0) *(pId++) = (double)resInt;
		if (allStokesReq) //OC18042020
		{
			if (RadExtract.pExtractedData != 0)
			{
				*pI1 = (float)resInt1; *pI2 = (float)resInt2; *pI3 = (float)resInt3;
			}
			else
			{
				*pI1d = resInt1; *pI2d = resInt2; *pI3d = resInt3;
			}
		}
	}
}

//int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ_GPU(srTRadExtract& RadExtract, double* arAuxInt, long long ie0, long long ie1, double InvStepRelArg, TGPUUsageArg* pGPU)
int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ_GPU(srTRadExtract& RadExtract, long long ie0, long long ie1, double InvStepRelArg, TGPUUsageArg* pGPU) //HG31072024
{
#define GEN_MEMBERS(i) \
	ExtractSingleElecIntensity2DvsXZ_Kernel<false, false, i, false>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<false, false, i, true>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<false, true,  i, false>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<false, true,  i, true>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<true,  false, i, false>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<true,  false, i, true>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<true,  true,  i, false>, \
	ExtractSingleElecIntensity2DvsXZ_Kernel<true,  true,  i, true>,

	decltype(ExtractSingleElecIntensity2DvsXZ_Kernel<false, false, 0, false>) *ExtractSingleElecIntensity2DvsXZ_tbl[] = {
		GEN_MEMBERS(-4)
		GEN_MEMBERS(-3)
		GEN_MEMBERS(-2)
		GEN_MEMBERS(-1)
		GEN_MEMBERS(0)
		GEN_MEMBERS(1)
		GEN_MEMBERS(2)
		GEN_MEMBERS(3)
		GEN_MEMBERS(4)
		GEN_MEMBERS(5)
		GEN_MEMBERS(-5)
	};
#undef GEN_MEMBERS

	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

    dim3 blocks(RadAccessData.nx, RadAccessData.nz);
    dim3 threads(1);

    if (RadAccessData.pBaseRadX != NULL)
	{
		RadAccessData.pBaseRadX = CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadX, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz, CAuxGPU::DISCARD_HOST);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadX);
	}
	if (RadAccessData.pBaseRadZ != NULL)
	{
		RadAccessData.pBaseRadZ = CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadZ, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz, CAuxGPU::DISCARD_HOST);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadZ);
	}

	srTRadGenManip *local_copy = CAuxGPU::ToDevice(pGPU, this, 1, CAuxGPU::DISCARD_HOST);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, local_copy);

	srTSRWRadStructAccessData* pRadAccessData_dev = CAuxGPU::ToDevice(pGPU, &RadAccessData, 1, CAuxGPU::DISCARD_HOST);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData_dev);

	bool allStokesReq = (RadExtract.PolarizCompon == -5);
	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1);

	int Int_or_ReE = RadExtract.Int_or_Phase;
	if (Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here
	
	if (Int_or_ReE != 2) //HG13012024 Fixed bug: Output array was not allocated properly
	{
		if (allStokesReq)
		{
			RadExtract.pExtractedData = CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedData, 4*RadAccessData.nx*RadAccessData.nz, CAuxGPU::DONT_COPY);
			CAuxGPU::Memset(pGPU, RadExtract.pExtractedData, 0.0f, 4*RadAccessData.nx*RadAccessData.nz);
		}
		else
		{
			RadExtract.pExtractedData = CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedData, RadAccessData.nx*RadAccessData.nz, CAuxGPU::DONT_COPY);
			CAuxGPU::Memset(pGPU, RadExtract.pExtractedData, 0.0f, RadAccessData.nx*RadAccessData.nz);
		}
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedData);
	}
	else
	{
		if (allStokesReq)
		{
			RadExtract.pExtractedDataD = CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedDataD, 4*RadAccessData.nx*RadAccessData.nz, CAuxGPU::DONT_COPY);
			CAuxGPU::Memset(pGPU, RadExtract.pExtractedDataD, 0.0, 4*RadAccessData.nx*RadAccessData.nz);
		}
		else
		{
			RadExtract.pExtractedDataD = CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedDataD, RadAccessData.nx*RadAccessData.nz, CAuxGPU::DONT_COPY);
			CAuxGPU::Memset(pGPU, RadExtract.pExtractedDataD, 0.0, RadAccessData.nx*RadAccessData.nz);
		}
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedDataD);
	}

	bool NpIsEven = ((RadAccessData.ne % 2) == 0);
	int idx = RadExtract.PolarizCompon;
	idx = (((idx < -4 || idx > 5) ? 10 : (idx + 4)) << 3) | ((allStokesReq & 1) << 2) | ((intOverEnIsRequired & 1) << 1) | (NpIsEven & 1);
	
	CAuxGPU::CalcLaunchDims(ExtractSingleElecIntensity2DvsXZ_tbl[idx], blocks, blocks, threads);
	ExtractSingleElecIntensity2DvsXZ_tbl[idx]<<<blocks, threads>>>(RadExtract, pRadAccessData_dev, local_copy, ie0, ie1, InvStepRelArg, Int_or_ReE);
	
	if(Int_or_ReE != 2) //HG13012024 Fixed bug: Output array was not allocated properly
	{
		if(RadExtract.pExtractedData != NULL)
		{
			CAuxGPU::MarkUpdated(pGPU, RadExtract.pExtractedData, CAuxGPU::DEVICE);
			RadExtract.pExtractedData = CAuxGPU::GetHostPtr(pGPU, RadExtract.pExtractedData);
		}
	}
	else
	{
		if(RadExtract.pExtractedDataD != NULL)
		{
			CAuxGPU::MarkUpdated(pGPU, RadExtract.pExtractedDataD, CAuxGPU::DEVICE);
			RadExtract.pExtractedDataD = CAuxGPU::GetHostPtr(pGPU, RadExtract.pExtractedDataD);
		}
	}

	CAuxGPU::MarkUpdatedBatch(pGPU, CAuxGPU::DEVICE, RadAccessData.pBaseRadX, RadAccessData.pBaseRadZ, pRadAccessData_dev, local_copy);
	CAuxGPU::ToHostAndFree(pGPU, pRadAccessData_dev); //HG27072024
//HG26022024 (commented out)
//#ifdef _DEBUG
//	if(Int_or_ReE != 2)
//	{
//		if (RadExtract.pExtractedData != NULL)
//			RadExtract.pExtractedData = (float*)CAuxGPU::ToHostAndFree(pGPU, RadExtract.pExtractedData, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
//	}
//	else
//	{
//		if (RadExtract.pExtractedDataD != NULL)
//			RadExtract.pExtractedDataD = (double*)CAuxGPU::ToHostAndFree(pGPU, RadExtract.pExtractedDataD, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(double));
//	}
//#endif

    CAuxGPU::ToHostAndFree(pGPU, local_copy);
	//CAuxGPU::ToHostAndFree(pGPU, arAuxInt, RadAccessData.ne*sizeof(double), true); //HG31072024
    //CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadX, true, false);
	//CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadZ, true, false);

	//if (RadAccessData.pBaseRadX != NULL)
	//	RadAccessData.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadX, 2 * RadAccessData.ne * RadAccessData.nx * RadAccessData.nz * sizeof(float), true); //HG13012024 Original wavefront data does not need to be copied back to CPU
	//if (RadAccessData.pBaseRadZ != NULL)
	//	RadAccessData.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadZ, 2 * RadAccessData.ne * RadAccessData.nx * RadAccessData.nz * sizeof(float), true); //HG13012024 Original wavefront data does not need to be copied back to CPU

	if (RadAccessData.pBaseRadX != NULL)
		RadAccessData.pBaseRadX = CAuxGPU::GetHostPtr(pGPU, RadAccessData.pBaseRadX); //HG13012024 Original wavefront data does not need to be copied back to CPU
	if (RadAccessData.pBaseRadZ != NULL)
		RadAccessData.pBaseRadZ = CAuxGPU::GetHostPtr(pGPU, RadAccessData.pBaseRadZ); //HG13012024 Original wavefront data does not need to be copied back to CPU

//HG26022024 (commented out)
//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
	return 0;
}

template <int PolCom, bool EhOK, bool EvOK, int gt1_iter, int itPerBlk>
__global__ void ExtractSingleElecMutualIntensityVsXZ_Kernel(const float* __restrict__ pEx0, const float* __restrict__ pEz0, float* __restrict__ pMI0, long nxnz, long itStart, long itEnd, long PerX, long iter0)
{
	//Calculate coordinates as the typical triangular matrix
	int i0 = (blockIdx.x * blockDim.x + threadIdx.x); //<=nxnz range
	int it0_0 = (blockIdx.y * blockDim.y + threadIdx.y); //nxnz/(2*itPerBlk) range
	long iter = iter0;

	if (i0 > nxnz) return;
	if (it0_0 > nxnz / 2) return;

	for (int it0 = it0_0 * itPerBlk; it0 < it0_0 * itPerBlk + itPerBlk; it0++)
	{
		long it = it0;
		long i = i0;
		if (i0 > it0) //If the coordinates are past the triangular bounds, switch to the lower half of the triangle
		{
			it = nxnz - it0 - 1;
			i = i0 - (it0 + 1);
		}

		if (it >= itEnd) {
			return;
		}

		//float* pMI = pMI0 + it0 * (nxnz << 1) + (i0 << 1); //Compact representation coordinates
		float* pMI = pMI0 + (it - itStart) * (nxnz << 1) + (i << 1); //Full representation coordinates
		const float* pEx = pEx0 + i * PerX;
		const float* pEz = pEz0 + i * PerX;
		const float* pExT = pEx0 + (it - itStart) * PerX;
		const float* pEzT = pEz0 + (it - itStart) * PerX;

		float ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
		float ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;

		{
			if (EhOK)
			{
				ExRe = *pEx; ExIm = *(pEx + 1);
				if (i != (it - itStart)) {
					ExReT = *pExT; ExImT = *(pExT + 1);
				}
				else {
					ExReT = ExRe;
					ExImT = ExIm;
				}
			}
			if (EvOK) {
				EzRe = *pEz; EzIm = *(pEz + 1);
				if (i != (it - itStart)) {
					EzReT = *pEzT; EzImT = *(pEzT + 1);
				}
				else {
					EzReT = EzRe;
					EzImT = EzIm;
				}
			}
		}
		float ReMI = 0., ImMI = 0.;

		switch (PolCom)
		{
		case 0: // Lin. Hor.
		{
			ReMI = ExRe * ExReT + ExIm * ExImT;
			ImMI = ExIm * ExReT - ExRe * ExImT;
			break;
		}
		case 1: // Lin. Vert.
		{
			ReMI = EzRe * EzReT + EzIm * EzImT;
			ImMI = EzIm * EzReT - EzRe * EzImT;
			break;
		}
		case 2: // Linear 45 deg.
		{
			float ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
			float ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
			ReMI = 0.5f * (ExRe_p_EzRe * ExRe_p_EzReT + ExIm_p_EzIm * ExIm_p_EzImT);
			ImMI = 0.5f * (ExIm_p_EzIm * ExRe_p_EzReT - ExRe_p_EzRe * ExIm_p_EzImT);
			break;
		}
		case 3: // Linear 135 deg.
		{
			float ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
			float ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
			ReMI = 0.5f * (ExRe_mi_EzRe * ExRe_mi_EzReT + ExIm_mi_EzIm * ExIm_mi_EzImT);
			ImMI = 0.5f * (ExIm_mi_EzIm * ExRe_mi_EzReT - ExRe_mi_EzRe * ExIm_mi_EzImT);
			break;
		}
		case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 4: // Circ. Right
		{
			float ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
			float ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
			ReMI = 0.5f * (ExRe_mi_EzIm * ExRe_mi_EzImT + ExIm_p_EzRe * ExIm_p_EzReT);
			ImMI = 0.5f * (ExIm_p_EzRe * ExRe_mi_EzImT - ExRe_mi_EzIm * ExIm_p_EzReT);
			break;
		}
		case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 5: // Circ. Left
		{
			float ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
			float ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
			ReMI = 0.5f * (ExRe_p_EzIm * ExRe_p_EzImT + ExIm_mi_EzRe * ExIm_mi_EzReT);
			ImMI = 0.5f * (ExIm_mi_EzRe * ExRe_p_EzImT - ExRe_p_EzIm * ExIm_mi_EzReT);
			break;
		}
		case -1: // s0
		{
			ReMI = ExRe * ExReT + ExIm * ExImT + EzRe * EzReT + EzIm * EzImT;
			ImMI = ExIm * ExReT - ExRe * ExImT + EzIm * EzReT - EzRe * EzImT;
			break;
		}
		case -2: // s1
		{
			ReMI = ExRe * ExReT + ExIm * ExImT - (EzRe * EzReT + EzIm * EzImT);
			ImMI = ExIm * ExReT - ExRe * ExImT - (EzIm * EzReT - EzRe * EzImT);
			break;
		}
		case -3: // s2
		{
			ReMI = ExImT * EzIm + ExIm * EzImT + ExReT * EzRe + ExRe * EzReT;
			ImMI = ExReT * EzIm - ExRe * EzImT - ExImT * EzRe + ExIm * EzReT;
			break;
		}
		case -4: // s3
		{
			ReMI = ExReT * EzIm + ExRe * EzImT - ExImT * EzRe - ExIm * EzReT;
			ImMI = ExIm * EzImT - ExImT * EzIm - ExReT * EzRe + ExRe * EzReT;
			break;
		}
		default: // total mutual intensity, same as s0
		{
			ReMI = ExRe * ExReT + ExIm * ExImT + EzRe * EzReT + EzIm * EzImT;
			ImMI = ExIm * ExReT - ExRe * ExImT + EzIm * EzReT - EzRe * EzImT;
			break;
			//return CAN_NOT_EXTRACT_MUT_INT;
		}
		}

		if (gt1_iter > 0)
		{
			pMI[0] = (pMI[0] * iter + (float)ReMI) / (float)(iter + 1.);
			pMI[1] = (pMI[1] * iter + (float)ImMI) / (float)(iter + 1.);
		}
		else if (gt1_iter == 0)
		{
			pMI[0] = (float)ReMI;
			pMI[1] = (float)ImMI;
		}
		else
		{
			pMI[0] += (float)ReMI;
			pMI[1] += (float)ImMI;
		}
	}
}

//template <int PolCom, int gt1_iter>
//int ExtractSingleElecMutualIntensityVsXZ_GPUSub(float* pEx, float* pEz, float* pMI0, long nx, long nz, long ne, long itStart, long itEnd, long PerX, long iter, int PolCom, bool EhOK, bool EvOK, TGPUUsageArg* pGPU)
int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ_GPU(float* pEx, float* pEz, float* pMI0, long nx, long nz, long ne, long itStart, long itEnd, long PerX, long iter, int PolCom, bool EhOK, bool EvOK, TGPUUsageArg* pGPU)
{
#define GEN_MEMBERS0(i, a, b) \
	ExtractSingleElecMutualIntensityVsXZ_Kernel<i, a, b, 0, 1>, \
	ExtractSingleElecMutualIntensityVsXZ_Kernel<i, a, b, 1, 1>, \
	ExtractSingleElecMutualIntensityVsXZ_Kernel<i, a, b, -1, 1>, \
	NULL,

#define GEN_MEMBERS(i) \
	GEN_MEMBERS0(i, false, false) \
	GEN_MEMBERS0(i, false, true) \
	GEN_MEMBERS0(i, true, false) \
	GEN_MEMBERS0(i, true, true)

	decltype(ExtractSingleElecMutualIntensityVsXZ_Kernel<0, false, false, 0, 1>) *ExtractSingleElecMutualIntensityVsXZ_tbl[] = {
		GEN_MEMBERS(-5)
		GEN_MEMBERS(-4)
		GEN_MEMBERS(-3)
		GEN_MEMBERS(-2)
		GEN_MEMBERS(-1)
		GEN_MEMBERS(0)
		GEN_MEMBERS(1)
		GEN_MEMBERS(2)
		GEN_MEMBERS(3)
		GEN_MEMBERS(4)
		GEN_MEMBERS(5)
	};
#undef GEN_MEMBERS0
#undef GEN_MEMBERS

	long long nxnz = ((long long)nx) * ((long long)nz); //HG26022024 NOTE: GPU implementation is only called for nxnz < UINT_MAX to avoid integer overflows

	const int itPerBlk = 1;
	dim3 threads = dim3(48, 16, 1);
	dim3 grid = dim3((unsigned int)((nxnz + 1) / threads.x + (threads.x > 1)), (unsigned int)((nxnz / 2) / (threads.y * itPerBlk) + (threads.y > 1)), 1); //OC19022024 (cast to remove warning)
	//dim3 grid = dim3((nxnz + 1) / threads.x + (threads.x > 1), (nxnz / 2) / (threads.y * itPerBlk) + (threads.y > 1), 1);

	pEx = CAuxGPU::ToDevice(pGPU, pEx, nxnz*2, CAuxGPU::DISCARD_HOST);
	pEz = CAuxGPU::ToDevice(pGPU, pEz, nxnz*2, CAuxGPU::DISCARD_HOST);
	pMI0 = CAuxGPU::ToDevice(pGPU, pMI0, (itEnd - itStart)*nxnz*2);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pEx, pEz, pMI0); //HG31072024 Bug-fix

	int idx = ((PolCom + 5) << 4) | ((EhOK & 1) << 3) | ((EvOK & 1) << 2);
	if (iter > 0) idx |= 1;
	else if (iter < 0) idx |= 2;

	ExtractSingleElecMutualIntensityVsXZ_tbl[idx]<<<grid, threads >>> (pEx, pEz, pMI0, (long)nxnz, itStart, itEnd, PerX, iter);

	CAuxGPU::MarkUpdatedBatch(pGPU, CAuxGPU::DEVICE, pEx, pEz, pMI0);
	pEx = CAuxGPU::ToHostAndFree(pGPU, pEx);
	pEz = CAuxGPU::ToHostAndFree(pGPU, pEz);

//HG26022024 (commented out)
//#ifdef _DEBUG
//	if (pMI0 != NULL)
//		pMI0 = (float*)CAuxGPU::ToHostAndFree(pGPU, pMI0, (itEnd - itStart)*ne*nx*nz*2*sizeof(float));
//
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
	return 0;
}

#endif