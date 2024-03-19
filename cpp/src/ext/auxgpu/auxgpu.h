/************************************************************************//**
 * File: auxgpu.h
 * Description: Auxiliary utilities to manage GPU usage
 * Project: Synchrotron Radiation Workshop
 * First release: 2023
 *
 * Copyright (C) Brookhaven National Laboratory
 * All Rights Reserved
 *
 * @author H.Goel
 * @version 1.0
 ***************************************************************************/

#ifndef __UTIGPU_H
#define __UTIGPU_H

#include <cstdlib>
#include <stdio.h>

#ifdef _OFFLOAD_GPU
#include <cuda_runtime.h>
#include <map>
//#if CUDART_VERSION < 11020
//#error CUDA version too low, need at least 11.2
//#endif
#endif

//typedef struct
struct TGPUUsageArg //OC18022024
{
	int deviceIndex; // -1 means no device, TODO

	TGPUUsageArg(void* pvGPU=0) //OC18022024
	{
		deviceIndex = -1;
		if(pvGPU == 0) return;
		double *arParGPU = (double*)pvGPU;
		int nPar = (int)arParGPU[0];
		if(nPar > 0) deviceIndex = (int)arParGPU[1];
		//continue here for future params
	}
}; 
//} TGPUUsageArg; //OC18022024 (commented-out)

#ifdef _OFFLOAD_GPU
#define GPU_COND(arg, code) if (arg && CAuxGPU::GPUEnabled((TGPUUsageArg*)arg)) { code }
//#define GPU_COND(arg, code) if (arg && CAuxGPU::GPUEnabled(arg)) { code }
#define GPU_PORTABLE __device__ __host__
#else
#define GPU_COND(arg, code) if(0) { }
#define GPU_PORTABLE 
#endif

 //*************************************************************************
class CAuxGPU
{
private:
public:
	static void Init();
	static void Fini();
	static bool GPUAvailable(); //CheckGPUAvailable etc
	static bool GPUEnabled(TGPUUsageArg *arg);
	static void SetGPUStatus(bool enabled);
	static int GetDevice(TGPUUsageArg* arg);
	static void* ToDevice(TGPUUsageArg* arg, void* hostPtr, size_t size, bool dontCopy = false);
	static void* GetHostPtr(TGPUUsageArg* arg, void* devicePtr);
	static void* ToHostAndFree(TGPUUsageArg* arg, void* devicePtr, size_t size, bool dontCopy = false);
	static void EnsureDeviceMemoryReady(TGPUUsageArg* arg, void* devicePtr);
	static void FreeHost(void* ptr);
	static void MarkUpdated(TGPUUsageArg* arg, void* ptr, bool devToHost, bool hostToDev);
};

//*************************************************************************
#endif