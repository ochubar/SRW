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
template<class T> __global__ void RadPointModifierParallel_Kernel(srTSRWRadStructAccessData RadAccessData, void* pBufVars, T* tgt_obj)
{
	int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
	int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range

	if (ix < RadAccessData.nx && iz < RadAccessData.nz)
	{
		srTEFieldPtrs EPtrs;
		srTEXZ EXZ;
		EXZ.z = RadAccessData.zStart + iz * RadAccessData.zStep;
		EXZ.x = RadAccessData.xStart + ix * RadAccessData.xStep;

		for (int ie = 0; ie < RadAccessData.ne; ie++) {
			EXZ.e = RadAccessData.eStart + ie * RadAccessData.eStep;
			EXZ.aux_offset = RadAccessData.ne * RadAccessData.nx * 2 * iz + RadAccessData.ne * 2 * ix + ie * 2;
			if (RadAccessData.pBaseRadX != 0)
			{
				EPtrs.pExRe = RadAccessData.pBaseRadX + EXZ.aux_offset;
				EPtrs.pExIm = EPtrs.pExRe + 1;
			}
			else
			{
				EPtrs.pExRe = 0;
				EPtrs.pExIm = 0;
			}
			if (RadAccessData.pBaseRadZ != 0)
			{
				EPtrs.pEzRe = RadAccessData.pBaseRadZ + EXZ.aux_offset;
				EPtrs.pEzIm = EPtrs.pEzRe + 1;
			}
			else
			{
				EPtrs.pEzRe = 0;
				EPtrs.pEzIm = 0;
			}

			tgt_obj->RadPointModifierPortable(EXZ, EPtrs, pBufVars);
		}
	}
}

template<class T> int RadPointModifierParallelImpl(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, T* tgt_obj, TGPUUsageArg* pGPU)
{
	const int bs = 256;
	dim3 blocks(pRadAccessData->nx / bs + ((pRadAccessData->nx & (bs - 1)) != 0), pRadAccessData->nz);
	dim3 threads(bs, 1);
	
	if (pRadAccessData->pBaseRadX != NULL)
	{
		pRadAccessData->pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadX, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadX);
	}
	if (pRadAccessData->pBaseRadZ != NULL)
	{
		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, pRadAccessData->pBaseRadZ, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pRadAccessData->pBaseRadZ);
	}

    T* local_copy = (T*)CAuxGPU::ToDevice(pGPU, tgt_obj, sizeof(T));
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, local_copy);
    //cudaMalloc(&local_copy, sizeof(T));
    //cudaMemcpy(local_copy, tgt_obj, sizeof(T), cudaMemcpyHostToDevice);
	
	void* pBufVars_dev = NULL;
	if (pBufVarsSz > 0){
		pBufVars_dev = CAuxGPU::ToDevice(pGPU, pBufVars, pBufVarsSz);
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pBufVars_dev);
	}
	RadPointModifierParallel_Kernel<T> << <blocks, threads >> > (*pRadAccessData, pBufVars_dev, local_copy);
    //cudaDeviceSynchronize();
    //cudaFreeAsync(local_copy, 0);
	if (pBufVarsSz > 0) CAuxGPU::ToHostAndFree(pGPU, pBufVars_dev, pBufVarsSz, true);
	CAuxGPU::ToHostAndFree(pGPU, local_copy, sizeof(T), true);

	CAuxGPU::MarkUpdated(pGPU, pRadAccessData->pBaseRadX, true, false);
	CAuxGPU::MarkUpdated(pGPU, pRadAccessData->pBaseRadZ, true, false);

//#ifndef _DEBUG //HG26022024 (commented-out)
	if (pRadAccessData->pBaseRadX != NULL)
		pRadAccessData->pBaseRadX = (float*)CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadX);
	if (pRadAccessData->pBaseRadZ != NULL)
		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::GetHostPtr(pGPU, pRadAccessData->pBaseRadZ);
//#endif

//HG26022024 (commented-out)
//#ifdef _DEBUG
//	if (pRadAccessData->pBaseRadX != NULL)
//		pRadAccessData->pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, pRadAccessData->pBaseRadX, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
//	if (pRadAccessData->pBaseRadZ != NULL)
//		pRadAccessData->pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, pRadAccessData->pBaseRadZ, 2*pRadAccessData->ne*pRadAccessData->nx*pRadAccessData->nz*sizeof(float));
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif

	return 0;
}
#endif

#endif //__SROPTELMGPU_H
#endif