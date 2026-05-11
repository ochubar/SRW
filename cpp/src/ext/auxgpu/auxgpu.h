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

#include <cstdarg>
#include <cstdlib>
#include <stdio.h>
#include <typeinfo>
//#include <type_traits> //OCTEST05122025

#ifdef _OFFLOAD_GPU
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <map>
#include <initializer_list>
//#if CUDART_VERSION < 11020
//#error CUDA version too low, need at least 11.2
//#endif
//#endif //HG23102025

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

//#ifdef _OFFLOAD_GPU //HG23102025
#define GPU_COND(arg, code) if (arg && CAuxGPU::GPUEnabled((TGPUUsageArg*)arg)) { code }
//#define GPU_COND(arg, code) if (arg && CAuxGPU::GPUEnabled(arg)) { code }
#define GPU_PORTABLE __device__ __host__

//#include <execinfo.h>
//#include <unistd.h>
//#include <link.h>
//#include <string.h>
//static int callback(struct dl_phdr_info *info, size_t size, void *data){
//	//If shared library name contains "srwlpy.so" print the base address
//	if(strstr(info->dlpi_name, "srwlpy.so") != NULL) printf("srwlpy.so loaded at address %llx\n", (unsigned long long)info->dlpi_addr);
//    return 0;
//}
//inline void backtrace_printer()
//{
//    void *array[10];
//    size_t size;
//    dl_iterate_phdr(callback, NULL);
//    size = backtrace(array, 10);
//    fprintf(stderr, "Error: CUDA error occurred\n");
//    backtrace_symbols_fd(array, size, STDERR_FILENO);
//	exit(1);
//}
//
//#define CUDA_SAFE(x) { cudaError_t err = x; if (err != cudaSuccess) { printf("CUDA error %d - %s\n", err, cudaGetErrorString(err)); backtrace_printer(); } }
#define CUDA_SAFE(x) { cudaError_t err = x; if (err != cudaSuccess) { printf("CUDA error %d - %s at %s:%d\n", err, cudaGetErrorString(err), __FILE__, __LINE__); } } //HG23102025 Add easy way to find and report CUDA errors
//#else //HG23102025 Not needed
//#define GPU_COND(arg, code) if(0) { }
//#define GPU_PORTABLE 
//#endif

const int PerThread = 16;
#ifdef __CUDACC__
	template<typename T> __global__ void Memset_Kernel(T* p, T val, long long n)
	{
		//Memset kernel, uses coalesced loads across warps
		long long i = blockIdx.x * blockDim.x + threadIdx.x;
		i = i * PerThread;
		if (i >= n) return;
		long long iFin = min(i + PerThread - 1, n - 1);
		for (; i <= iFin; i++) p[i] = val;
	}
#else
    template<typename T> __global__ void Memset_Kernel(T* p, T val, long long n);
#endif

 //*************************************************************************
class CAuxGPU
{
private:

#ifdef _OFFLOAD_GPU
	typedef struct //HG02082024 Move from cpp file to here to make the map a class member
	{
		void *devicePtr;
		void *hostPtr;
		size_t size;
		bool HostToDevUpdated;
		bool DevToHostUpdated;
		cudaEvent_t h2d_event;
		cudaEvent_t d2h_event;
		//bool pinned; //HG26072024
		int flags; //HG23102025
	} memAllocInfo_t;
	static int current_device;
	static std::map<void*, memAllocInfo_t> gpuMap;
	//static bool memcpy_stream_initialized = false; //HG02082024 (commented-out)
	static std::map<int, cudaStream_t*> streams; //HG02082024
	static cudaStream_t memcpy_stream;
#endif

	//static void* ToDevice(TGPUUsageArg* arg, void* hostPtr, size_t size, bool dontCopy = false); //HG26072024
	static void* _ToDevice(TGPUUsageArg* arg, void* hostPtr, size_t size, int flags=0); //HG26072024 Make private

	static void* _ToHostAndFree(TGPUUsageArg* arg, void* devicePtr, int flags=0, size_t size=0); //HG26072024
	//static void* ToHostAndFree(TGPUUsageArg* arg, void* devicePtr, size_t size, bool dontCopy = false);
	
	static void* _GetHostPtr(TGPUUsageArg* arg, void* devicePtr); //HG26072024
	//static void* GetHostPtr(TGPUUsageArg* arg, void* devicePtr);
public:
	/**
	 * Flags used by this class
	 */
	static constexpr int DONT_COPY = (1 << 0);
	static constexpr int PIN_ON_HOST = (1 << 1);
	static constexpr int HOST = (1 << 2);
	static constexpr int DEVICE = (1 << 3);
	static constexpr int DISCARD_HOST = (1 << 4); //HG13102025

	/**
	* Initialize GPU/device functionality
	*/
	//static void Init();
	static void Init(TGPUUsageArg *arg); //HG02082024

	/**
	* Call when returning to the client layer to ensure all memory is accessible on CPU/host again
	*/
	//static void Fini();
	static void Fini(TGPUUsageArg *arg); //HG02082024

	static bool GPUAvailable(); //CheckGPUAvailable etc
	static bool GPUEnabled(TGPUUsageArg *arg);
	static void SetGPUStatus(bool enabled);

	/**
	*  Get the GPU/device number associated with arg
	*  @param [in] arg pointer to a GPU usage argument structure
	*  @return integer number of the GPU/device, -1 if CPU/host
	*/
	static int GetDevice(TGPUUsageArg* arg);

	/**
	*  Associate the specified region of host memory with memory on the device, copies the memory to device by default
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] hostPtr pointer to the region of host memory
	* @param [in] size size in bytes of the memory region
	* @param [in] flags flags to control the memory transfer (DONT_COPY, PIN_ON_HOST, DISCARD_HOST)
	* @return pointer to device memory, NULL on error
	*/
	template <typename T>
	static T* ToDevice(TGPUUsageArg* arg, T* hostPtr, size_t elemCount, int flags=0) //HG30042025
	{
#ifdef _OFFLOAD_GPU
		const int typeSize = sizeof(T);
		return (T*)_ToDevice(arg, (void*)hostPtr, elemCount * typeSize, flags);
#endif
		return hostPtr;
	}
//	//OCTEST05122025
//	template <typename T>
//	static T* ToDevice(TGPUUsageArg* arg, T* hostPtr, size_t elemCount, int flags = 0)
//	{
//#ifdef _OFFLOAD_GPU
//		size_t byteCount;
//
//		if constexpr(std::is_void_v<T>) {
//			// For T = void, elemCount is already in bytes
//			byteCount = elemCount;
//		}
//		else {
//			byteCount = elemCount * sizeof(T);
//		}
//
//		return static_cast<T*>(_ToDevice(arg, static_cast<void*>(hostPtr), byteCount, flags));
//#else
//		(void)arg; (void)elemCount; (void)flags; // avoid unused warnings if you like
//		return hostPtr;
//#endif
//	}

	template <typename T>
	static void Memset(TGPUUsageArg* arg, T* devicePtr, T value, size_t elemCount) //HG30042025
	{
#if defined(_OFFLOAD_GPU)
		if (arg == NULL)
			return;
		if (arg->deviceIndex == 0)
			return;
		if (!GPUEnabled(arg))
			return;
		if (devicePtr == NULL)
			return;
		if (elemCount == 0)
			return;
		//if (gpuMap.find(devicePtr) != gpuMap.end()){
		//	void* devPtr = devicePtr;
		//	if (gpuMap[devPtr].DevToHostUpdated){
		//		cudaStreamWaitEvent(memcpy_stream, gpuMap[devPtr].d2h_event);
		//		gpuMap[devPtr].DevToHostUpdated = false;
		//		if (gpuMap[devPtr].hostPtr != NULL) gpuMap[gpuMap[devPtr].hostPtr].DevToHostUpdated = false;
		//	}
//
		//	int minGridSize = 0;
		//	int bs = 256;
		//	size_t elemCount_orig = elemCount;
		//	cudaOccupancyMaxPotentialBlockSize(&minGridSize, &bs, Memset_Kernel<T>, 0, (elemCount + PerThread - 1) / PerThread);
		//	elemCount = ((elemCount + PerThread - 1) / PerThread + bs - 1) / bs;
		//	Memset_Kernel<T> <<<elemCount, bs, 0, memcpy_stream>>> ((T*)devPtr, value, elemCount_orig);
		//}
		T* ptr = CAuxGPU::ToDevice(arg, devicePtr, elemCount, CAuxGPU::DONT_COPY);
		CAuxGPU::EnsureDeviceMemoryReady(arg, ptr);
		dim3 blocks((unsigned)(elemCount / PerThread + !!(elemCount % PerThread)), 1, 1), threads(1); //HG09122025
		//dim3 blocks(elemCount / PerThread + !!(elemCount % PerThread), 1, 1), threads(1);
		CAuxGPU::CalcLaunchDims(Memset_Kernel<T>, blocks, blocks, threads);
		void* args[] = { &ptr, &value, &elemCount };
		cudaLaunchKernel((void*)Memset_Kernel<T>, blocks, threads, args, 0, 0);
		CAuxGPU::MarkUpdated(arg, ptr, CAuxGPU::DEVICE);
#endif
	}
	
	/**
	* Retrieve the host memory address for a given device or host pointer
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] devicePtr pointer for which the host pointer is desired
	* @return the corresponding host pointer, NULL on errror
	*/
	template<typename T>
	static T* GetHostPtr(TGPUUsageArg* arg, T* devicePtr)
	{
		return (T*)_GetHostPtr(arg, devicePtr);
	}

	/**
	* Transfer memory back to the host if necessary and free the associated device memory. Does not return until the latest copy of the data is on the host.
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] devicePtr device pointer to the memory to be freed, if a host pointer is provided, the corresponding device pointer is freed
	* @param [in] size size of the block to be freed
	* @param [in] flags flags to control the memory transfer (DONT_COPY)
	* @return The corresponding host pointer, NULL on error
	*/
	template <typename T>
	static T* ToHostAndFree(TGPUUsageArg* arg, T* devicePtr, size_t elemCount=0, int flags=0) //HG30042025
	{
#ifdef _OFFLOAD_GPU
		const int typeSize = (typeid(T) == typeid(void)) ? 1 : sizeof(T);
		return (T*)_ToHostAndFree(arg, (void*)devicePtr, flags, elemCount * typeSize);
#endif
		return devicePtr;
	}

	/**
	* Ensure that the device memory has the latest data, used prior to kernel launches
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] devicePtr device pointer to the memory block to be operated on
	*/
	static void EnsureDeviceMemoryReady(TGPUUsageArg* arg, void* devicePtr=0);

	// Makes it possible to pass multiple pointers to EnsureDeviceMemoryReady in a single call
	template <typename First, typename... T> 
	static void EnsureDeviceMemoryReady(TGPUUsageArg* arg, First* devicePtr, T*... ptrs) //HG30042025	
	{
#ifdef _OFFLOAD_GPU
		if (arg == NULL)
			return;
		if (arg->deviceIndex == 0)
			return;
		if (!GPUEnabled(arg))
			return;
		EnsureDeviceMemoryReady(arg, (void*)devicePtr);
		EnsureDeviceMemoryReady(arg, ptrs...);
#endif
	}

	//static void FreeHost(void* ptr); //HG26072024 (Commented out) Unused and potentially breaks this memory management model

	/**
	* If origPtr is a host pointer that has a corresponding device memory block, reassign that block to correspond to newPtr instead, otherwise copy the data from origPtr to newPtr on host
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] origPtr original host pointer
	* @param [in] newPtr host pointer to replace it with
	* @param [in] size size of this memory region
	* @param [in] flags flags to control the memory transfer (DISCARD_HOST)
	* @return 0 on success, -1 on error
	*/
	static int SetHostPtr(TGPUUsageArg* arg, void* origPtr, void* newPtr, size_t size, int flags=0); //HG26072024

	/**
	* Mark the region as having been updated.
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] ptr pointer to the memory region, can be a host or device pointer
	* @param [in] flags flags to control the memory transfer (HOST: host to device, DEVICE: device to host)
	*/
	static void MarkUpdated(TGPUUsageArg* arg, void* ptr, int flags=0); //HG26072024
	//static void MarkUpdated(TGPUUsageArg* arg, void* ptr, bool devToHost, bool hostToDev);

	static void MarkUpdatedBatch(TGPUUsageArg* arg, int flags) //HG30042025 Base case for recursion
	{
		return;
	}

	template<typename First, typename... T> 
	static void MarkUpdatedBatch(TGPUUsageArg* arg, int flags, First* ptr, T*... ptrs) //HG30042025
	{
#ifdef _OFFLOAD_GPU
		if (arg == NULL)
			return;
		if (arg->deviceIndex == 0)
			return;
		if (!GPUEnabled(arg))
			return;
		MarkUpdated(arg, (void*)ptr, flags);
		MarkUpdatedBatch(arg, flags, ptrs...);
#endif
	}

	/**
	* Retrieve a compute stream index to run kernels simultaneously on one GPU.
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] idx the 0-based index of the desired compute stream
	* @return The cudaStream ID associated with the requested compute stream index
	*/
	static long long GetComputeStream(TGPUUsageArg* arg, int idx); //HG26072024

	/**
	* Ensure that the specified compute stream is synchronized with the target compute stream
	* @param [in] arg pointer to a GPU usage argument structure
	* @param [in] targetStreamIdx the 0-based index of the target compute stream
	* @param [in] streamIdx the 0-based index of the compute stream to be synchronized
	*/
	static void SyncComputeStream(TGPUUsageArg* arg, long long targetStreamIdx, long long streamIdx); //HG24042025

	/**
	* Determine a good distribution of threads within a block
	* @param [in] func The kernel to calculate for.
	* @param [in] grid The grid size.
	* @param [out] blocks The resulting block count.
	* @param [out] threads The resulting thread count.
	* @param [in] max_x_bs The maximum block size in the X dimension.
	* @param [in] max_bs_total The maximum block size the kernel is designed to work with.
	* @param [in] warpMultiple Whether the block size should be a multiple of the warp size (32).
	*/
	template<typename T>
	static void CalcLaunchDims(T* kern, dim3 grid, dim3& blocks, dim3& threads, int max_x_bs = 0, int max_bs_total = 1024, bool warpMultiple = false)
	{
		int minGridSize = 0; //HG05082024
    	int bs = max_bs_total;

		//OCTEST (attempted to synchronize or eliminate)
		//cudaDeviceSynchronize(); //HG05082024 Ensure all previous kernels have completed before calculating occupancy for this kernel
		CUDA_SAFE(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &bs, (void*)kern, 0, max_bs_total));
		//bs = 32;
		//OCTEST END

		if (warpMultiple && (bs % 32 != 0)) bs -= (bs % 32);
		if ((max_x_bs == 0) || (grid.x < (unsigned)bs && grid.x < (unsigned)max_x_bs)) max_x_bs = grid.x; //OC09122025
		//if ((max_x_bs == 0) || (grid.x < bs && grid.x < max_x_bs)) max_x_bs = grid.x;
		if (warpMultiple && (max_x_bs % 32 != 0)) max_x_bs -= (max_x_bs % 32);

		threads.x = (max_x_bs < bs) ? max_x_bs : bs;
		threads.y = 1;
		threads.z = 1;

		blocks.x = grid.x / threads.x + !!(grid.x % threads.x); //round up the division result
		blocks.y = grid.y;
		blocks.z = grid.z;

		int y_v = bs / max_x_bs;
		if (y_v > 1 && grid.y > 1)
		{
			//Calculate y grid
			if ((unsigned)y_v >= grid.y) threads.y = grid.y; //OC09122025
			//if (y_v >= grid.y) threads.y = grid.y;
			else
			{
				//Check the remainder if threads.y where set to y_v
				int blky = grid.y / y_v + !!(grid.y % y_v); //round up the division result
				int rem = (blky * y_v) % grid.y;
				if (rem / (float)y_v > 0.1)
				{
					//Adjust threads.y to minimize the remainder
					threads.y = grid.y / blky + !!(grid.y % blky);
					if (threads.y > 32 && threads.y % 32 != 0) threads.y -= (threads.y % 32);
				}
				else threads.y = y_v;
			}
			blocks.y = grid.y / threads.y + !!(grid.y % threads.y); //round up the division result

			int z_v = y_v / threads.y;
			if (z_v > 1 && grid.z > 1)
			{
				threads.z = ((unsigned)z_v > grid.z) ? grid.z : (unsigned)z_v; //OC09122025
				//threads.z = (z_v > grid.z) ? grid.z : z_v;
				blocks.z = grid.z / threads.z + !!(grid.z % threads.z); //round up the division result
			}
		}
	}
};

template <>
inline void* CAuxGPU::ToDevice<void>(TGPUUsageArg* arg, void* hostPtr, size_t elemCount, int flags) //HG05122025
{
	return ToDevice<char>(arg, (char*)hostPtr, elemCount, flags);
}

template<>
inline void* CAuxGPU::ToHostAndFree<void>(TGPUUsageArg* arg, void* devicePtr, size_t elemCount, int flags) //HG05122025
{
	return (void*)ToHostAndFree(arg, (char*)devicePtr, elemCount, flags);
}

//*************************************************************************
#endif //HG23102025
#endif

