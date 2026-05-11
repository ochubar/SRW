/************************************************************************//**
 * File: sroptelm_gpu.h
 * Description: Optical element (general CUDA header)
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
#ifndef __SROPTELMGPU_H
#define __SROPTELMGPU_H

#include "cuda_runtime.h"
#include <sroptelm.h>
#include <srradstr.h>
#include <srstraux.h>

#ifdef __CUDACC__
template<class T, bool combinedE> 
__global__ void 
__launch_bounds__(256)
//RadPointModifierParallel_Kernel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, T* tgt_obj, int xStart, int xFin, int zStart, int zFin) //HG27072024 Redesigned entire function
TraverseRadZXEParallel_Kernel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, T* tgt_obj, int xStart, int xFin, int zStart, int zFin) //HG14042026
{
	int ie = (blockIdx.x * blockDim.x + threadIdx.x); //ne range
	int ix = (blockIdx.y * blockDim.y + threadIdx.y) + xStart; //nx range
	int iz = (blockIdx.z * blockDim.z + threadIdx.z) + zStart; //nz range
	
	int ne = 1;
	if (combinedE)
	{
		ne = pRadAccessData->ne;
		ie = 0;
	} 

	if (ix < xFin && iz < zFin && ie < pRadAccessData->ne) //HG27072024 changed RadAccessData to pRadAccessData
	{
		srTEFieldPtrs EPtrs;
		srTEXZ EXZ;
		EXZ.z = pRadAccessData->zStart + iz * pRadAccessData->zStep;
		EXZ.x = pRadAccessData->xStart + ix * pRadAccessData->xStep;
		EXZ.e = pRadAccessData->eStart + ie * pRadAccessData->eStep;
		EXZ.aux_offset = pRadAccessData->ne * pRadAccessData->nx * 2 * iz + pRadAccessData->ne * 2 * ix + ie * 2;
		if (pRadAccessData->pBaseRadX != 0)
		{
			EPtrs.pExRe = pRadAccessData->pBaseRadX + EXZ.aux_offset;
			EPtrs.pExIm = EPtrs.pExRe + 1;
		}
		else
		{
			EPtrs.pExRe = 0;
			EPtrs.pExIm = 0;
		}
		if (pRadAccessData->pBaseRadZ != 0)
		{
			EPtrs.pEzRe = pRadAccessData->pBaseRadZ + EXZ.aux_offset;
			EPtrs.pEzIm = EPtrs.pEzRe + 1;
		}
		else
		{
			EPtrs.pEzRe = 0;
			EPtrs.pEzIm = 0;
		}

		tgt_obj->RadPointModifierPortable(EXZ, EPtrs, pBufVars);

		for (ie=1; ie < ne; ie++)
		{
			EXZ.e += pRadAccessData->eStep;
			EXZ.aux_offset += 2;
			if (pRadAccessData->pBaseRadX != 0)
			{
				EPtrs.pExRe += 2;
				EPtrs.pExIm += 2;
			}
			if (pRadAccessData->pBaseRadZ != 0)
			{
				EPtrs.pEzRe += 2;
				EPtrs.pEzIm += 2;
			}
			tgt_obj->RadPointModifierPortable(EXZ, EPtrs, pBufVars);
		}
	}
}

template<class T> 
//int RadPointModifierParallelImpl(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, T* tgt_obj, TGPUUsageArg* pGPU, int *pRegion=0, bool combinedE=false) //HG29072024
int TraverseRadZXEParallelImpl(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, T* tgt_obj, TGPUUsageArg* pGPU, int *pRegion=0, bool combinedE=false) //HG14042026
{
	if (pRadAccessData->pBaseRadX != NULL)
	{
		pRadAccessData->pBaseRadX = CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadX, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadX);
	}
	if (pRadAccessData->pBaseRadZ != NULL)
	{
		pRadAccessData->pBaseRadZ = CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadZ, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadZ);
	}
	
	srTSRWRadStructAccessData* pRadAccessData_dev = CAuxGPU::ToDevice(pGPU, pRadAccessData, 1, CAuxGPU::DISCARD_HOST);
    T* local_copy = CAuxGPU::ToDevice(pGPU, tgt_obj, 1, CAuxGPU::DISCARD_HOST);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData_dev, local_copy);
	
	void* pBufVars_dev = NULL;
	if (pBufVarsSz > 0)
	{
		pBufVars_dev = CAuxGPU::ToDevice(pGPU, (char*)pBufVars, pBufVarsSz, CAuxGPU::DISCARD_HOST);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pBufVars_dev);
	}
	
	int xStart = 0;
	int xFin = pRadAccessData->nx;
	int zStart = 0;
	int zFin = pRadAccessData->nz;
	
	//HG30072024 Allow for specifying a region to skip or to only process within the region, reduces extra operations for propagators like apertures and obstacles
	bool HandleInSingleLaunch = (pRegion == 0);
	if (!HandleInSingleLaunch)
	{
		if (pRegion[4] == 0)
		{
			xStart = pRegion[0];
			xFin = pRegion[1];
			zStart = pRegion[2];
			zFin = pRegion[3];
			HandleInSingleLaunch = true;
		}
	}

	void (*kern)(srTSRWRadStructAccessData*, void*, T*, int, int, int, int) = NULL;
	//if (combinedE) kern = RadPointModifierParallel_Kernel<T, true>;
	if (combinedE) kern = TraverseRadZXEParallel_Kernel<T, true>; //HG14042026
	//else kern = RadPointModifierParallel_Kernel<T, false>;
	else kern = TraverseRadZXEParallel_Kernel<T, false>; //HG14042026
	if (HandleInSingleLaunch)
	{
		dim3 blocks(combinedE ? 1 : pRadAccessData->ne, xFin - xStart, zFin - zStart);
		dim3 threads(1);
		CAuxGPU::CalcLaunchDims(kern, blocks, blocks, threads);
		kern<<<blocks, threads >>> (pRadAccessData_dev, pBufVars_dev, local_copy, xStart, xFin, zStart, zFin);
	}
	else
	{
		//Have to split into 4 kernel launches to skip the specified region, run them in parallel
		cudaStream_t stream1 = (cudaStream_t)CAuxGPU::GetComputeStream(pGPU, 0);
		cudaStream_t stream2 = (cudaStream_t)CAuxGPU::GetComputeStream(pGPU, 1);
		cudaStream_t stream3 = (cudaStream_t)CAuxGPU::GetComputeStream(pGPU, 2);

		CAuxGPU::SyncComputeStream(pGPU, 0, (long long)stream1);
		CAuxGPU::SyncComputeStream(pGPU, 0, (long long)stream2);
		CAuxGPU::SyncComputeStream(pGPU, 0, (long long)stream3);

		dim3 blocks0(combinedE ? 1 : pRadAccessData->ne, pRegion[0], pRadAccessData->nz);
		dim3 blocks1(combinedE ? 1 : pRadAccessData->ne, pRadAccessData->nx - pRegion[1], pRadAccessData->nz);
		dim3 blocks2(combinedE ? 1 : pRadAccessData->ne, pRegion[1] - pRegion[0], pRadAccessData->nz - pRegion[3]);
		dim3 blocks3(combinedE ? 1 : pRadAccessData->ne, pRegion[1] - pRegion[0], pRegion[2]);
		dim3 threads0(1);
		dim3 threads1(1);
		dim3 threads2(1);
		dim3 threads3(1);

		if (blocks0.y > 0 && blocks0.z > 0)
		{
			CAuxGPU::CalcLaunchDims(kern, blocks0, blocks0, threads0);
			kern<<<blocks0, threads0 >>> (pRadAccessData_dev, pBufVars_dev, local_copy, 0, pRegion[0], 0, pRadAccessData->nz);
		}
		if (blocks1.y > 0 && blocks1.z > 0)
		{
			CAuxGPU::CalcLaunchDims(kern, blocks1, blocks1, threads1);
			kern<<<blocks1, threads1, 0, stream1 >>> (pRadAccessData_dev, pBufVars_dev, local_copy, pRegion[1], pRadAccessData->nx, 0, pRadAccessData->nz);
		}
		if (blocks2.y > 0 && blocks2.z > 0)
		{
			CAuxGPU::CalcLaunchDims(kern, blocks2, blocks2, threads2);
			kern<<<blocks2, threads2, 0, stream2 >>> (pRadAccessData_dev, pBufVars_dev, local_copy, pRegion[0], pRegion[1], pRegion[3], pRadAccessData->nz);
		}
		if (blocks3.y > 0 && blocks3.z > 0)
		{
			CAuxGPU::CalcLaunchDims(kern, blocks3, blocks3, threads3);
			kern<<<blocks3, threads3, 0, stream3 >>> (pRadAccessData_dev, pBufVars_dev, local_copy, pRegion[0], pRegion[1], 0, pRegion[2]);
		}
		
		CAuxGPU::SyncComputeStream(pGPU, (long long)stream1, 0);
		CAuxGPU::SyncComputeStream(pGPU, (long long)stream2, 0);
		CAuxGPU::SyncComputeStream(pGPU, (long long)stream3, 0);
	}

	CAuxGPU::MarkUpdatedBatch(pGPU, CAuxGPU::DEVICE, pRadAccessData->pBaseRadX, pRadAccessData->pBaseRadZ, pRadAccessData_dev, pBufVars_dev, local_copy);
	if (pBufVarsSz > 0) CAuxGPU::ToHostAndFree(pGPU, (char*)pBufVars_dev);
	CAuxGPU::ToHostAndFree(pGPU, pRadAccessData_dev);
	CAuxGPU::ToHostAndFree(pGPU, local_copy);
	
	if (pRadAccessData->pBaseRadX != NULL)
		pRadAccessData->pBaseRadX = CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadX);
	if (pRadAccessData->pBaseRadZ != NULL)
		pRadAccessData->pBaseRadZ = CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadZ);
	return 0;
}

//template<class T> __global__ void RadPointModifierParallel_Kernel(srTSRWRadStructAccessData RadAccessData, void* pBufVars, T* tgt_obj)
//{
//	int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
//	int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range
//
//	if (ix < RadAccessData.nx && iz < RadAccessData.nz)
//	{
//		srTEFieldPtrs EPtrs;
//		srTEXZ EXZ;
//		EXZ.z = RadAccessData.zStart + iz * RadAccessData.zStep;
//		EXZ.x = RadAccessData.xStart + ix * RadAccessData.xStep;
//
//		for (int ie = 0; ie < RadAccessData.ne; ie++) {
//			EXZ.e = RadAccessData.eStart + ie * RadAccessData.eStep;
//			EXZ.aux_offset = RadAccessData.ne * RadAccessData.nx * 2 * iz + RadAccessData.ne * 2 * ix + ie * 2;
//			if (RadAccessData.pBaseRadX != 0)
//			{
//				EPtrs.pExRe = RadAccessData.pBaseRadX + EXZ.aux_offset;
//				EPtrs.pExIm = EPtrs.pExRe + 1;
//			}
//			else
//			{
//				EPtrs.pExRe = 0;
//				EPtrs.pExIm = 0;
//			}
//			if (RadAccessData.pBaseRadZ != 0)
//			{
//				EPtrs.pEzRe = RadAccessData.pBaseRadZ + EXZ.aux_offset;
//				EPtrs.pEzIm = EPtrs.pEzRe + 1;
//			}
//			else
//			{
//				EPtrs.pEzRe = 0;
//				EPtrs.pEzIm = 0;
//			}
//
//			tgt_obj->RadPointModifierPortable(EXZ, EPtrs, pBufVars);
//		}
//	}
//}
//
//template<class T> int RadPointModifierParallelImpl(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, T* tgt_obj, TGPUUsageArg* pGPU)
//{
//	const int bs = 256;
//	dim3 blocks(pRadAccessData->nx / bs + ((pRadAccessData->nx & (bs - 1)) != 0), pRadAccessData->nz);
//	dim3 threads(bs, 1);
//	
//	if (pRadAccessData->pBaseRadX != NULL)
//	{
//		pRadAccessData->pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadX, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
//		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadX);
//	}
//	if (pRadAccessData->pBaseRadZ != NULL)
//	{
//		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadZ, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
//		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadZ);
//	}
//
//    T* local_copy = (T*)CAuxGPU::ToDevice(pGPU, tgt_obj, sizeof(T));
//	CAuxGPU::EnsureDeviceMemoryReady(pGPU, local_copy);
//    //cudaMalloc(&local_copy, sizeof(T));
//    //cudaMemcpy(local_copy, tgt_obj, sizeof(T), cudaMemcpyHostToDevice);
//	
//	void* pBufVars_dev = NULL;
//	if (pBufVarsSz > 0){
//		pBufVars_dev = CAuxGPU::ToDevice(pGPU, pBufVars, pBufVarsSz);
//		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pBufVars_dev);
//	}
//	RadPointModifierParallel_Kernel<T> << <blocks, threads >> > (*pRadAccessData, pBufVars_dev, local_copy);
//    //cudaDeviceSynchronize();
//    //cudaFreeAsync(local_copy, 0);
//	if (pBufVarsSz > 0) CAuxGPU::ToHostAndFree(pGPU, pBufVars_dev, pBufVarsSz, true);
//	CAuxGPU::ToHostAndFree(pGPU, local_copy, sizeof(T), true);
//
//	CAuxGPU::MarkUpdated(pGPU, pRadAccessData->pBaseRadX, true, false);
//	CAuxGPU::MarkUpdated(pGPU, pRadAccessData->pBaseRadZ, true, false);
//
////#ifndef _DEBUG //HG26022024 (commented-out)
//	if (pRadAccessData->pBaseRadX != NULL)
//		pRadAccessData->pBaseRadX = (float*)CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadX);
//	if (pRadAccessData->pBaseRadZ != NULL)
//		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadZ);
////#endif
//
////HG26022024 (commented-out)
////#ifdef _DEBUG
////	if (pRadAccessData->pBaseRadX != NULL)
////		pRadAccessData->pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, pRadAccessData->pBaseRadX, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
////	if (pRadAccessData->pBaseRadZ != NULL)
////		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, pRadAccessData->pBaseRadZ, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
////	cudaStreamSynchronize(0);
////	auto err = cudaGetLastError();
////	printf("%s\r\n", cudaGetErrorString(err));
////#endif
//
//	return 0;
//}
#endif

#endif //__SROPTELMGPU_H
#endif