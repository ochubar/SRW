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
#include <string>
#include <assert.h>
#include <math.h>
#include "srradmnp.h"
#include "gmmeth.h"

template <bool allStokesReq, bool intOverEnIsRequired, int PolCom>
__global__ void ExtractSingleElecIntensity2DvsXZ_Kernel(srTRadExtract RadExtract, srTSRWRadStructAccessData RadAccessData, srTRadGenManip *obj, double* arAuxInt, long long ie0, long long ie1, double InvStepRelArg, int Int_or_ReE)
{
	int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range
    
	if (ix < RadAccessData.nx && iz < RadAccessData.nz) 
    {
		//int PolCom = RadExtract.PolarizCompon;
			
		//bool allStokesReq = (PolCom == -5); //OC18042020

		float* pI = 0, * pI1 = 0, * pI2 = 0, * pI3 = 0; //OC17042020
		double* pId = 0, * pI1d = 0, * pI2d = 0, * pI3d = 0;
		long ne = RadAccessData.ne, nx = RadAccessData.nx, nz = RadAccessData.nz;
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

		float* pEx0 = RadAccessData.pBaseRadX;
		float* pEz0 = RadAccessData.pBaseRadZ;

		//long PerX = RadAccessData.ne << 1;
		//long PerZ = PerX*RadAccessData.nx;
		//long long PerX = RadAccessData.ne << 1;
		//long long PerZ = PerX*RadAccessData.nx;
		long long PerX = ((long long)ne) << 1; //OC18042020
		long long PerZ = PerX * nx;

		//bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (ne > 1); //OC18042020
		double resInt, resInt1, resInt2, resInt3;
		double ConstPhotEnInteg = 1.;
		long long Two_ie0 = ie0 << 1, Two_ie1 = ie1 << 1; //OC26042019
		long ie;

		long offset = iz * PerZ + ix * PerX;
		long offsetDiv2 = offset >> 1;

		float* pEx_StartForX = pEx0 + offset;
		float* pEz_StartForX = pEz0 + offset;
		if (pI != 0)
		{
			pI += offsetDiv2;
			if (allStokesReq)
			{
				pI1 += offsetDiv2;
				pI2 += offsetDiv2;
				pI3 += offsetDiv2;
			}
		} 

		if (pId != 0)
		{
			pId += offsetDiv2;
			if (allStokesReq)
			{
				pI1d += offsetDiv2;
				pI2d += offsetDiv2;
				pI3d += offsetDiv2;
			}
		} 
		
		//long ixPerX = 0;

		float* pEx_St = pEx_StartForX + Two_ie0;
		float* pEz_St = pEz_StartForX + Two_ie0;
		float* pEx_Fi = pEx_StartForX + Two_ie1;
		float* pEz_Fi = pEz_StartForX + Two_ie1;

		if (intOverEnIsRequired) //OC140813
		{//integrate over photon energy / time
			double* tInt = arAuxInt;
			float* pEx_StAux = pEx_St;
			float* pEz_StAux = pEz_St;

			if (!allStokesReq) //OC17042020
			{
				for (ie = 0; ie < ne; ie++) //OC18042020
				//for(int ie=0; ie<RadAccessData.ne; ie++)
				{
					*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, PolCom, Int_or_ReE);
					pEx_StAux += 2;
					pEz_StAux += 2;
				}
				resInt = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep); //OC18042020
				//resInt = ConstPhotEnInteg*CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, RadAccessData.ne, RadAccessData.eStep);
			}
			else
			{
				for (ie = 0; ie < ne; ie++)
				{
					*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -1, Int_or_ReE);
					pEx_StAux += 2; pEz_StAux += 2;
				}
				resInt = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				for (ie = 0; ie < ne; ie++)
				{
					*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -2, Int_or_ReE);
					pEx_StAux += 2; pEz_StAux += 2;
				}
				resInt1 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				for (ie = 0; ie < ne; ie++)
				{
					*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -3, Int_or_ReE);
					pEx_StAux += 2; pEz_StAux += 2;
				}
				resInt2 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);

				tInt = arAuxInt; pEx_StAux = pEx_St; pEz_StAux = pEz_St;
				for (ie = 0; ie < ne; ie++)
				{
					*(tInt++) = obj->IntensityComponent(pEx_StAux, pEz_StAux, -4, Int_or_ReE);
					pEx_StAux += 2; pEz_StAux += 2;
				}
				resInt3 = ConstPhotEnInteg * CGenMathMeth::Integ1D_FuncDefByArray(arAuxInt, ne, RadAccessData.eStep);
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

template <bool allStokesReq, bool intOverEnIsRequired>
static inline void ExtractSingleElecIntensity2DvsXZ_GPUSub(dim3 &blocks, dim3 &threads, srTRadExtract RadExtract, srTSRWRadStructAccessData RadAccessData, srTRadGenManip *local_copy, double* arAuxInt, long long ie0, long long ie1, double InvStepRelArg, int Int_or_ReE)
{
	switch(RadExtract.PolarizCompon)
	{
		case 5: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 5><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case 4: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 4><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case 3: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 3><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case 2: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 2><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case 1: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 1><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case 0: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, 0><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case -1: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, -1><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case -2: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, -2><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case -3: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, -3><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		case -4: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, -4><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
		default: ExtractSingleElecIntensity2DvsXZ_Kernel<allStokesReq, intOverEnIsRequired, -5><<<blocks, threads>>>(RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE); break;
	}
}

int srTRadGenManip::ExtractSingleElecIntensity2DvsXZ_GPU(srTRadExtract& RadExtract, double* arAuxInt, long long ie0, long long ie1, double InvStepRelArg, TGPUUsageArg* pGPU)
{
	srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

    const int bs = 256;
    dim3 blocks(RadAccessData.nx / bs + ((RadAccessData.nx & (bs - 1)) != 0), RadAccessData.nz);
    dim3 threads(bs, 1);

    if (RadAccessData.pBaseRadX != NULL)
	{
		RadAccessData.pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadX, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadX);
	}
	if (RadAccessData.pBaseRadZ != NULL)
	{
		RadAccessData.pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadZ, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadZ);
	}

	srTRadGenManip *local_copy = (srTRadGenManip*)CAuxGPU::ToDevice(pGPU, this, sizeof(srTRadGenManip));
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, local_copy);

    arAuxInt = (double*)CAuxGPU::ToDevice(pGPU, arAuxInt, RadAccessData.ne*sizeof(double));
    CAuxGPU::EnsureDeviceMemoryReady(pGPU, arAuxInt);

	bool allStokesReq = (RadExtract.PolarizCompon == -5);
	bool intOverEnIsRequired = (RadExtract.Int_or_Phase == 7) && (RadAccessData.ne > 1);

	int Int_or_ReE = RadExtract.Int_or_Phase;
	if (Int_or_ReE == 7) Int_or_ReE = 0; //OC150813: time/phot. energy integrated single-e intensity requires "normal" intensity here
	
	if (Int_or_ReE != 2) //HG13012024 Fixed bug: Output array was not allocated properly
	{
		if (allStokesReq)
		{
			RadExtract.pExtractedData = (float*)CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedData, 4*RadAccessData.nx*RadAccessData.nz*sizeof(float), true);
			//CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedData);
		}
		else
		{
			RadExtract.pExtractedData = (float*)CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedData, RadAccessData.nx*RadAccessData.nz*sizeof(float), true);
			//CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedData);
		}
	}
	else
	{
		if (allStokesReq)
		{
			RadExtract.pExtractedDataD = (double*)CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedDataD, 4*RadAccessData.nx*RadAccessData.nz*sizeof(double), true);
			//CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedDataD);
		}
		else
		{
			RadExtract.pExtractedDataD = (double*)CAuxGPU::ToDevice(pGPU, RadExtract.pExtractedDataD, RadAccessData.nx*RadAccessData.nz*sizeof(double), true);
			//CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadExtract.pExtractedDataD);
		}
	}

	if (allStokesReq)
		if (intOverEnIsRequired)
			ExtractSingleElecIntensity2DvsXZ_GPUSub<true, true> (blocks, threads, RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE);
		else
			ExtractSingleElecIntensity2DvsXZ_GPUSub<true, false> (blocks, threads, RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE);
	else
		if (intOverEnIsRequired)
			ExtractSingleElecIntensity2DvsXZ_GPUSub<false, true> (blocks, threads, RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE);
		else
			ExtractSingleElecIntensity2DvsXZ_GPUSub<false, false> (blocks, threads, RadExtract, RadAccessData, local_copy, arAuxInt, ie0, ie1, InvStepRelArg, Int_or_ReE);
	
	if(Int_or_ReE != 2) //HG13012024 Fixed bug: Output array was not allocated properly
	{
		if(RadExtract.pExtractedData != NULL)
			CAuxGPU::MarkUpdated(pGPU, RadExtract.pExtractedData, true, false);
	}
	else
	{
		if(RadExtract.pExtractedDataD != NULL)
			CAuxGPU::MarkUpdated(pGPU, RadExtract.pExtractedDataD, true, false);
	}

//#ifndef _DEBUG //HG26022024 (commented out)
	if(Int_or_ReE != 2)
	{
		if (RadExtract.pExtractedData != NULL)
			RadExtract.pExtractedData = (float*)CAuxGPU::GetHostPtr(pGPU, RadExtract.pExtractedData);
	}
	else
	{
		if (RadExtract.pExtractedDataD != NULL)
			RadExtract.pExtractedDataD = (double*)CAuxGPU::GetHostPtr(pGPU, RadExtract.pExtractedDataD);
	}
//#endif

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

    CAuxGPU::ToHostAndFree(pGPU, local_copy, sizeof(srTRadGenManip), true);
    CAuxGPU::ToHostAndFree(pGPU, arAuxInt, RadAccessData.ne*sizeof(double), true);
	//CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadX, true, false);
	//CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadZ, true, false);
	if (RadAccessData.pBaseRadX != NULL)
		RadAccessData.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadX, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float), true); //HG13012024 Original wavefront data does not need to be copied back to CPU
	if (RadAccessData.pBaseRadZ != NULL)
		RadAccessData.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadZ, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float), true); //HG13012024 Original wavefront data does not need to be copied back to CPU

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

template <int PolCom, int gt1_iter>
int ExtractSingleElecMutualIntensityVsXZ_GPUSub(float* pEx, float* pEz, float* pMI0, long nx, long nz, long ne, long itStart, long itEnd, long PerX, long iter, bool EhOK, bool EvOK, TGPUUsageArg* pGPU)
{
	long long nxnz = ((long long)nx) * ((long long)nz); //HG26022024 NOTE: GPU implementation is only called for nxnz < UINT_MAX to avoid integer overflows

	const int itPerBlk = 1;
	dim3 threads = dim3(48, 16, 1);
	dim3 grid = dim3((unsigned int)((nxnz + 1) / threads.x + (threads.x > 1)), (unsigned int)((nxnz / 2) / (threads.y * itPerBlk) + (threads.y > 1)), 1); //OC19022024 (cast to remove warning)
	//dim3 grid = dim3((nxnz + 1) / threads.x + (threads.x > 1), (nxnz / 2) / (threads.y * itPerBlk) + (threads.y > 1), 1);

	pEx = (float*)CAuxGPU::ToDevice(pGPU, pEx, nxnz*2*sizeof(float));
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pEx);

	pEz = (float*)CAuxGPU::ToDevice(pGPU, pEz, nxnz*2*sizeof(float));
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pEz);

	pMI0 = (float*)CAuxGPU::ToDevice(pGPU, pMI0, (itEnd - itStart)*nxnz*2*sizeof(float));

	if (EhOK)
	{
		//HG26022024
		if (EvOK) ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, true, true, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, (long)nxnz, itStart, itEnd, PerX, iter);
		else ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, true, false, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, (long)nxnz, itStart, itEnd, PerX, iter);
		//if (EvOK) ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, true, true, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, nxnz, itStart, itEnd, PerX, iter);
		//else ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, true, false, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, nxnz, itStart, itEnd, PerX, iter);
	}
	else
	{
		//HG26022024
		if (EvOK) ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, false, true, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, (long)nxnz, itStart, itEnd, PerX, iter);
		else ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, false, false, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, (long)nxnz, itStart, itEnd, PerX, iter);
		//if (EvOK) ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, false, true, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, nxnz, itStart, itEnd, PerX, iter);
		//else ExtractSingleElecMutualIntensityVsXZ_Kernel<PolCom, false, false, gt1_iter, itPerBlk> << <grid, threads >> > (pEx, pEz, pMI0, nxnz, itStart, itEnd, PerX, iter);
	}

	pEx = (float*)CAuxGPU::ToHostAndFree(pGPU, pEx, nxnz * 2 * sizeof(float), true);
	pEz = (float*)CAuxGPU::ToHostAndFree(pGPU, pEz, nxnz * 2 * sizeof(float), true);
	
	CAuxGPU::MarkUpdated(pGPU, pMI0, true, false);

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

int srTRadGenManip::ExtractSingleElecMutualIntensityVsXZ_GPU(float* pEx, float* pEz, float* pMI0, long nx, long nz, long ne, long itStart, long itEnd, long PerX, long iter, int PolCom, bool EhOK, bool EvOK, TGPUUsageArg* pGPU)
{
	if(iter > 0)
	{
		switch (PolCom)
		{
		case  1: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  1, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  2: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  2, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  0: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  0, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  3: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  3, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  4: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  4, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  5: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  5, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -1: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -1, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -2: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -2, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -3: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -3, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -4: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -4, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		default: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -5, 1>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		}
	}
	else if(iter == 0)
	{
		switch (PolCom)
		{
		case  0: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  0, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  1: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  1, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  2: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  2, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  3: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  3, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  4: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  4, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case  5: return ExtractSingleElecMutualIntensityVsXZ_GPUSub<  5, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -1: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -1, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -2: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -2, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -3: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -3, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		case -4: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -4, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		default: return ExtractSingleElecMutualIntensityVsXZ_GPUSub< -5, 0>(pEx, pEz, pMI0, nx, nz, ne, itStart, itEnd, PerX, iter, EhOK, EvOK, pGPU);
		}
	}
	return 0; //OC19022024
}

#endif