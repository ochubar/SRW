/************************************************************************//**
 * File: utidev.h
 * Description: GPU offloading detection and control
 * Project: Synchrotron Radiation Workshop (and possibly others)
 * First release: 2022
 *
 * @author H. Goel
 * @version 0.1
 ***************************************************************************/

#ifndef __UTIGPU_H
#define __UTIGPU_H

#include <cstdlib>
#include <stdio.h>

#ifdef _OFFLOAD_GPU
#include <cuda_runtime.h>
#endif

typedef int gpuUsageArg_t;

#define ALLOC_ARRAY(type, size) (type *)UtiDev::malloc(sizeof(type)*(size))
#define FREE_ARRAY(x) UtiDev::free(x); x=NULL
#define ALLOC_STRUCT(type) (type *)UtiDev::malloc(sizeof(type))
#define FREE_STRUCT(x) UtiDev::free(x); x=NULL

#ifdef _OFFLOAD_GPU
#define GPU_ENABLED(arg) UtiDev::GPUEnabled(arg)
#define GPU_COND(arg, code) if (GPU_ENABLED(arg)) { code }
#define GPU_PORTABLE __device__ __host__
#else
#define GPU_COND(arg, code) if(0) { }
#define GPU_ENABLED(arg) 0
#define GPU_PORTABLE 
#endif

 //*************************************************************************
class UtiDev
{
public:
	static void Init();
	static void Fini();
	static bool GPUAvailable(); //CheckGPUAvailable etc
	static bool GPUEnabled(gpuUsageArg_t *arg);
	static void SetGPUStatus(bool enabled);
	static int GetDevice(gpuUsageArg_t* arg);
	
	static inline void* malloc(size_t sz) {
#ifdef _OFFLOAD_GPU
			void *ptr;
			auto err = cudaMallocManaged(&ptr, sz);
			if (err != cudaSuccess)
				printf("Allocation Failure\r\n");
			return ptr;
#else
			return std::malloc(sz);
#endif
	}

	static inline void free(void* ptr) {
#ifdef _OFFLOAD_GPU
		cudaFree(ptr);
#else
		std::free(ptr);
#endif
	}
};

//*************************************************************************
#endif