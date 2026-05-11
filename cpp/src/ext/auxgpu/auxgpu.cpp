/************************************************************************//**
 * File: auxgpu.cpp
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

#ifdef _OFFLOAD_GPU //HG09122025

#include <cstdio>
#include <cstdlib>
#include <new>
#include <cstring> //HG26072024

//#ifdef _OFFLOAD_GPU
#include <cuda_runtime.h>
//#endif

#include "auxgpu.h"

static bool isGPUAvailable = false;
static bool isGPUEnabled = false;
static bool GPUAvailabilityTested = false;
static bool deviceOffloadInitialized = false;
static int deviceCount = 0;

//#ifdef _OFFLOAD_GPU
//typedef struct
//{
//	void *devicePtr;
//	void *hostPtr;
//	size_t size;
//	bool HostToDevUpdated;
//	bool DevToHostUpdated;
//	cudaEvent_t h2d_event;
//	cudaEvent_t d2h_event;
//} memAllocInfo_t;
//static std::map<void*, memAllocInfo_t> gpuMap;
//static cudaStream_t memcpy_stream;
//static bool memcpy_stream_initialized = false;
//static int current_device = -1;
int CAuxGPU::current_device = -1; //HG02082024 Make into class members
std::map<void*, CAuxGPU::memAllocInfo_t> CAuxGPU::gpuMap; //HG02082024
std::map<int, cudaStream_t*> CAuxGPU::streams; //HG02082024
cudaStream_t CAuxGPU::memcpy_stream; //HG02082024
//#endif

static void CheckGPUAvailability() 
{
//#ifdef _OFFLOAD_GPU
	if (!GPUAvailabilityTested)
	{
		isGPUAvailable = false;
		GPUAvailabilityTested = true;
		int deviceCount = 0;
		if (cudaGetDeviceCount(&deviceCount) != cudaSuccess)
			return;

		if (deviceCount < 1)
			return;

		isGPUAvailable = true;
	}
//#else
//	isGPUAvailable = false;
//	isGPUEnabled = false;
//	GPUAvailabilityTested = true;
//#endif
}

bool CAuxGPU::GPUAvailable()
{
	CheckGPUAvailability();
	return isGPUAvailable;
}

bool CAuxGPU::GPUEnabled(TGPUUsageArg *arg) 
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return false;
	if (arg->deviceIndex > 0) {
		if (arg->deviceIndex <= deviceCount)
		{
			//if (memcpy_stream_initialized && current_device != arg->deviceIndex) //HG02082024 (commented-out)
			//{
			//	cudaStreamDestroy(memcpy_stream);
			//	memcpy_stream_initialized = false;
			//}
			//cudaSetDevice(arg->deviceIndex - 1);
			//if (!memcpy_stream_initialized)
			//	cudaStreamCreateWithFlags(&memcpy_stream, cudaStreamNonBlocking);
			//current_device = arg->deviceIndex;
			//memcpy_stream_initialized = true;

			if (current_device != arg->deviceIndex) //HG02082024
			{
				if (streams.find(arg->deviceIndex) == streams.end()) Init(arg);
				memcpy_stream = streams[arg->deviceIndex][0];
				current_device = arg->deviceIndex;
			}
		}
		//TODO: Add warning that GPU isn't available
		return GPUAvailable();
	}
//#endif
	return false;
}

void CAuxGPU::SetGPUStatus(bool enabled)
{
	isGPUEnabled = enabled && GPUAvailable();
}

int CAuxGPU::GetDevice(TGPUUsageArg* arg)
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return cudaCpuDeviceId;

	int curDevice = 0;
	cudaGetDevice(&curDevice);
	return curDevice;
//#else
//	return 0;
//#endif
}

//void* CAuxGPU::ToDevice(TGPUUsageArg* arg, void* hostPtr, size_t size, bool dontCopy)
void* CAuxGPU::_ToDevice(TGPUUsageArg* arg, void* hostPtr, size_t size, int flags) //HG26072024
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return hostPtr;
	if (arg->deviceIndex == 0)
		return hostPtr;
	//if (hostPtr == NULL)
	//	return hostPtr;
	if (size == 0)
		return hostPtr;
	if (!GPUEnabled(arg))
		return hostPtr;
	bool dontCopy = flags & DONT_COPY; //HG26072024
	bool pinOnHost = flags & PIN_ON_HOST; //HG26072024
	if (hostPtr == NULL) dontCopy = true; //HG26072024
	//if (gpuMap.find(hostPtr) != gpuMap.end()){
	if (hostPtr != NULL) //HG21042025
	{
		auto close_l = gpuMap.lower_bound(hostPtr);
		if (close_l != gpuMap.end())
		{
			memAllocInfo_t info = close_l->second;
			bool devMatchFound = (info.devicePtr <= hostPtr && (char*)info.devicePtr + info.size > hostPtr);
			bool matchFound = (info.hostPtr <= hostPtr && (char*)info.hostPtr + info.size > hostPtr) | devMatchFound;
			if (!matchFound)
			{
				info = std::prev(close_l)->second;
				devMatchFound = (info.devicePtr <= hostPtr && (char*)info.devicePtr + info.size > hostPtr);
				matchFound = (info.hostPtr <= hostPtr && (char*)info.hostPtr + info.size > hostPtr) | devMatchFound;
			}
			if (matchFound)
			{
				if (devMatchFound)
					return (char*)hostPtr; //Already on device
				size_t offset = (char*)hostPtr - (char*)info.hostPtr;
				void* devPtr = info.devicePtr;
				hostPtr = info.hostPtr;
				//if (gpuMap[devPtr].HostToDevUpdated && !dontCopy){
				if (gpuMap[devPtr].HostToDevUpdated && devPtr != NULL && hostPtr != NULL && !dontCopy){
					//cudaMemcpyAsync(devPtr, hostPtr, size, cudaMemcpyHostToDevice, memcpy_stream);
					CUDA_SAFE(cudaMemcpyAsync((char*)devPtr + offset, (char*)hostPtr + offset, size, cudaMemcpyDefault, memcpy_stream)); //HG26072024
					CUDA_SAFE(cudaEventRecord(gpuMap[devPtr].h2d_event, memcpy_stream));
				}
				//#if _DEBUG
				//		printf("ToDevice: %p -> %p, %d, D2H: %d, H2D: %d\n", hostPtr, devPtr, size, gpuMap[devPtr].DevToHostUpdated, gpuMap[devPtr].HostToDevUpdated); //HG28072023
				//#endif
				//gpuMap[devPtr].HostToDevUpdated = false; //HG23102025 Commented-out
				return (char*)devPtr + offset;
			}
		}
	}

	size_t free_mem = 0, total_mem = 0;  //HG26072024 If the memory request is very large, it may be more optimal to pin to host memory
	void *devicePtr = NULL;
	cudaError_t err = cudaSuccess;
	CUDA_SAFE(cudaMemGetInfo(&free_mem, &total_mem));
	if(size >= total_mem * 0.3 && !pinOnHost)
	{
		pinOnHost = true;
	}
	if (!pinOnHost) //Try to use asynchronous allocations first
	{
		err = cudaMallocAsync(&devicePtr, size, memcpy_stream); //Try asynchronous allocation
		if (err != cudaSuccess) // Try again after freeing up some memory HG24072023
		{
			cudaDeviceSynchronize();
			err = cudaMalloc(&devicePtr, size);
			if (err != cudaSuccess) pinOnHost = true; //HG26072024 If allocation still fails, try pinning on host
		}
	}
	if (hostPtr != NULL && pinOnHost) //Fallback to pinning host memory directly
	{
		err = cudaHostRegister(hostPtr, size, cudaHostRegisterDefault | cudaHostRegisterMapped);
		cudaHostGetDevicePointer(&devicePtr, hostPtr, 0);
	}
	if (err != cudaSuccess)
		return NULL;

	//void *devicePtr = NULL; //HG26072024 (commented-out)
	//cudaError_t err = cudaMalloc(&devicePtr, size);
	//if (err != cudaSuccess) // Try again after freeing up some memory HG24072023
	//{
	//	cudaStreamSynchronize(0);
	//	err = cudaMalloc(&devicePtr, size);
	//}
	//if (err != cudaSuccess)
	//	return NULL;
//#if _DEBUG
//	printf("ToDevice: %p -> %p, %d\n", hostPtr, devicePtr, size); //HG28072023
//#endif
	memAllocInfo_t info;
	info.devicePtr = devicePtr;
	info.hostPtr = hostPtr;
	info.DevToHostUpdated = false;
	info.HostToDevUpdated = false;
	//info.pinned = pinOnHost; //HG26072024
	info.flags = pinOnHost ? PIN_ON_HOST : 0; //HG13102025
	if(flags & DISCARD_HOST) info.flags |= DISCARD_HOST; //HG13102025
	cudaEventCreateWithFlags(&info.h2d_event, cudaEventDisableTiming);
	cudaEventCreateWithFlags(&info.d2h_event, cudaEventDisableTiming);
	//if (!dontCopy){ //HG26072024 (commented-out)
	//	cudaMemcpyAsync(devicePtr, hostPtr, size, cudaMemcpyHostToDevice, memcpy_stream);
	//	cudaEventRecord(info.h2d_event, memcpy_stream);
	//}
	info.HostToDevUpdated = true; //HG26072024
	if (hostPtr != NULL && !dontCopy) CUDA_SAFE(cudaMemcpyAsync(devicePtr, hostPtr, size, cudaMemcpyDefault, memcpy_stream)); //HG27072024 Add memset options
	//if (hostPtr != NULL && !dontCopy) printf("Memcpy: %llx %llx %d\r\n", hostPtr, devicePtr, size);
	CUDA_SAFE(cudaEventRecord(info.h2d_event, memcpy_stream));
	info.size = size;
	//gpuMap[hostPtr] = info;
	if(hostPtr != NULL) gpuMap[hostPtr] = info;
	gpuMap[devicePtr] = info;
	//if (hostPtr != NULL) printf("Stored hostPtr entry in gpuMap\r\n");
	//printf("Stored devicePtr entry in gpuMap\r\n");
	return devicePtr;
//#else
//	return hostPtr;
//#endif
}

void CAuxGPU::EnsureDeviceMemoryReady(TGPUUsageArg* arg, void* hostPtr)
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return;
	if (arg->deviceIndex == 0)
		return;
	if (hostPtr == NULL)
		return;
	if (!GPUEnabled(arg))
		return;
	//printf("%s %llx\r\n", __func__, hostPtr);
	if (gpuMap.find(hostPtr) != gpuMap.end()){
		void* devPtr = gpuMap[hostPtr].devicePtr;
		if (gpuMap[devPtr].HostToDevUpdated){
			CUDA_SAFE(cudaStreamWaitEvent(0, gpuMap[devPtr].h2d_event)); //HG23102025 Use CUDA_SAFE to catch errors
			gpuMap[devPtr].HostToDevUpdated = false; //HG26072024 Bug fix: After this event the latest data is known to be on device
			//gpuMap[gpuMap[devPtr].hostPtr].HostToDevUpdated = false; //HG26072024
			if (gpuMap[devPtr].hostPtr != NULL) gpuMap[gpuMap[devPtr].hostPtr].HostToDevUpdated = false; //HG26072024
		}
//#if _DEBUG
//		printf("EnsureDeviceMemoryReady: %p -> %p, %d, D2H: %d, H2D: %d\n", hostPtr, devPtr, gpuMap[devPtr].size, gpuMap[devPtr].DevToHostUpdated, gpuMap[devPtr].HostToDevUpdated); //HG28072023
//#endif
	}
//#endif
}

void* CAuxGPU::_GetHostPtr(TGPUUsageArg* arg, void* devicePtr)
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return devicePtr;
	if (arg->deviceIndex == 0)
		return devicePtr;
	if (devicePtr == NULL)
		return devicePtr;
	if (!GPUEnabled(arg))
		return devicePtr;
	memAllocInfo_t info;
	if (gpuMap.find(devicePtr) == gpuMap.end())
		return devicePtr;
	info = gpuMap[devicePtr];
//#if _DEBUG
//	printf("GetHostPtr: %p -> %p\n", devicePtr, info.hostPtr); //HG28072023
//#endif
	return info.hostPtr;
//#else
//	return devicePtr;
//#endif
}

//void* CAuxGPU::ToHostAndFree(TGPUUsageArg* arg, void* devicePtr, size_t size, bool dontCopy)
void* CAuxGPU::_ToHostAndFree(TGPUUsageArg* arg, void* devicePtr, int flags, size_t size) //HG26072024
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return devicePtr;
	if (arg->deviceIndex == 0)
		return devicePtr;
	if (devicePtr == NULL)
		return devicePtr;
	if (!GPUEnabled(arg))
		return devicePtr;
	memAllocInfo_t info;
	if (gpuMap.find(devicePtr) == gpuMap.end())
		return devicePtr;
	info = gpuMap[devicePtr];
	devicePtr = info.devicePtr;
	void *hostPtr = info.hostPtr;
	bool dontCopy = false; //HG26072024
	if (flags & DONT_COPY) dontCopy = true; //HG13102025
	if (hostPtr == NULL) dontCopy = true;
	if (size == 0) size = info.size;
	if (size == 0) dontCopy = true;
	if (info.flags & DISCARD_HOST) dontCopy = true; //HG13102025
	//if (!dontCopy && info.DevToHostUpdated) //HG26072024 (commented-out)
	//{
	//	cudaStreamWaitEvent(memcpy_stream, info.d2h_event, 0);
	//	cudaMemcpyAsync(hostPtr, devicePtr, size, cudaMemcpyDeviceToHost, memcpy_stream);
	//	cudaEventRecord(info.d2h_event);
	//	cudaEventSynchronize(info.d2h_event); // we can't treat host memory as valid until the copy is complete
	//}
//#if _DEBUG
//	printf("ToHostAndFree: %p -> %p, %d\n", devicePtr, hostPtr, size); //HG28072023
//#endif
	//cudaStreamWaitEvent(0, info.h2d_event); //HG26072024 (commented-out)
	//cudaStreamWaitEvent(0, info.d2h_event);
	//cudaFreeAsync(devicePtr, 0);
	if (hostPtr != NULL && !dontCopy && info.DevToHostUpdated) //HG26072024 Handle free-ing the different tiers of memory we can handle
	{
		CUDA_SAFE(cudaStreamWaitEvent(memcpy_stream, info.d2h_event, 0));
		if(!(info.flags & PIN_ON_HOST)) //HG30072024 Properly handle pinned memory
		{
			CUDA_SAFE(cudaMemcpyAsync(hostPtr, devicePtr, size, cudaMemcpyDefault, memcpy_stream)); //HG26072024 only copy if not using pinned memory
			CUDA_SAFE(cudaFreeAsync(devicePtr, memcpy_stream)); //HG26072024 Doing the async free here is  slightly more efficient and eliminates a potential use-after-free
			CUDA_SAFE(cudaEventSynchronize(info.d2h_event)); // we can't treat host memory as valid until the copy is complete
		}
		else
		{
			CUDA_SAFE(cudaMemcpyAsync(hostPtr, devicePtr, size, cudaMemcpyDefault, memcpy_stream));
			CUDA_SAFE(cudaEventSynchronize(info.d2h_event)); // we can't treat host memory as valid until the copy is complete
			CUDA_SAFE(cudaHostUnregister(devicePtr));
		}
	}
	else //HG26072024
	{
		if(!(info.flags & PIN_ON_HOST)) //HG30072024 Properly handle pinned memory
		{
			CUDA_SAFE(cudaStreamWaitEvent(0, info.d2h_event)); //HG26072024 H2D events are meaningless when the memory is on host
			CUDA_SAFE(cudaFreeAsync(devicePtr, 0));
		} 
		else
		{
			CUDA_SAFE(cudaEventSynchronize(info.d2h_event));
			CUDA_SAFE(cudaHostUnregister(devicePtr));
		} 
	}
    CUDA_SAFE(cudaEventDestroy(info.h2d_event));
	CUDA_SAFE(cudaEventDestroy(info.d2h_event));
	gpuMap.erase(devicePtr);
	//gpuMap.erase(hostPtr);
	if (hostPtr != NULL) gpuMap.erase(hostPtr); //HG21042025
	return hostPtr;
//#else
//	return devicePtr;
//#endif
}

//void CAuxGPU::FreeHost(void* ptr) //HG26072024 (commented-out)
//{
//#ifdef _OFFLOAD_GPU
//	if (ptr == NULL)
//		return;
//	if (gpuMap.find(ptr) == gpuMap.end())
//		return;
//	memAllocInfo_t info = gpuMap[ptr];
//	void *hostPtr = info.hostPtr;
//	void *devicePtr = info.devicePtr;
////#if _DEBUG
////	printf("FreeHost: %p, %p\n", devicePtr, hostPtr);
////#endif
//    //cudaStreamWaitEvent(0, info.h2d_event);
//	//cudaStreamWaitEvent(0, info.d2h_event);
//	cudaFreeAsync(devicePtr, 0);
//	//cudaEventDestroy(info.h2d_event);
//	//cudaEventDestroy(info.d2h_event);
//	std::free(hostPtr); //OC02082023
//	//CAuxGPU::free(hostPtr);
//	gpuMap.erase(devicePtr);
//	gpuMap.erase(hostPtr);
//#endif
//	return;
//}

int CAuxGPU::SetHostPtr(TGPUUsageArg* arg, void* origPtr, void* newPtr, size_t size, int flags) //HG26072024 Add function
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return -1;
	if (arg->deviceIndex == 0)
		return -1;
	if (origPtr == NULL)
		return -1;
	if (newPtr == NULL)
		return -1;
	if (!GPUEnabled(arg))
		return -1;
//#if _DEBUG
//	printf("SetHostPtr: %p -> %p\n", origPtr, newPtr);
//#endif
	if (gpuMap.find(origPtr) == gpuMap.end())
	{
		memcpy(newPtr, origPtr, size);
		ToDevice(arg, newPtr, size, flags); //Register the new pointer
		return 0;
	}
	memAllocInfo_t info = gpuMap[origPtr];
	if (gpuMap.find(newPtr) != gpuMap.end())
		return -1;	//The new pointer should not already be known to the GPU memory map, else we will have a memory leak

	if (info.DevToHostUpdated)
	{
		gpuMap.erase(origPtr);
		info.hostPtr = newPtr;
		if (flags & DISCARD_HOST) //HG13102025
		{
			info.flags |= DISCARD_HOST;
		}
		gpuMap[info.hostPtr] = info;
		gpuMap[info.devicePtr] = info;
	}
	else 
	{
		memcpy(newPtr, origPtr, size);
	}
//#endif
	return 0;
}

//void CAuxGPU::MarkUpdated(TGPUUsageArg* arg, void* ptr, bool devToHost, bool hostToDev)
void CAuxGPU::MarkUpdated(TGPUUsageArg* arg, void* ptr, int flags) //HG26072024
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return;
	if (arg->deviceIndex == 0)
		return;
	if (ptr == NULL)
		return;
	if (!GPUEnabled(arg))
		return;
	if (gpuMap.find(ptr) == gpuMap.end())
		return;
	bool devToHost = flags & DEVICE; //HG26072024
	bool hostToDev = flags & HOST; //HG26072024
	void* devPtr = gpuMap[ptr].devicePtr;
	void* hostPtr = gpuMap[ptr].hostPtr;
	if ((devToHost | gpuMap[devPtr].DevToHostUpdated) && (hostToDev | gpuMap[devPtr].HostToDevUpdated)) //HG26072024 Trying to perform both a copy to host and then a copy back to device doesn't make sense
		return;
	gpuMap[devPtr].DevToHostUpdated = devToHost;
	gpuMap[devPtr].HostToDevUpdated = hostToDev;
	//printf("%s %llx D2H:%d H2D:%d\r\n", __func__, devPtr, devToHost, hostToDev);
	if (hostPtr != NULL)
	{
		gpuMap[hostPtr].DevToHostUpdated = devToHost;
		gpuMap[hostPtr].HostToDevUpdated = hostToDev;
	}
	if (devToHost)
		CUDA_SAFE(cudaEventRecord(gpuMap[devPtr].d2h_event, 0));
	if (hostPtr != NULL && !(gpuMap[hostPtr].flags & PIN_ON_HOST) && hostToDev) //HG26072024 If host data has been updated, copy it over
	{
		CUDA_SAFE(cudaMemcpyAsync(devPtr, hostPtr, gpuMap[devPtr].size, cudaMemcpyDefault, memcpy_stream));
		CUDA_SAFE(cudaEventRecord(gpuMap[devPtr].h2d_event, memcpy_stream));
	}
//#if _DEBUG
//	printf("MarkUpdated: %p -> %p, D2H: %d, H2D: %d\n", ptr, devPtr, devToHost, hostToDev);
//#endif
//#endif
}

long long CAuxGPU::GetComputeStream(TGPUUsageArg* arg, int idx) //HG02082024 Add the ability to run multiple computations simultaneously on a GPU
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return 0;
	if (arg->deviceIndex == 0)
		return 0;
	if (idx < 0)
		return 0;
	if (idx > 8) idx = 8; //Max 9 streams (0 is memcpy stream)
	if (streams.find(arg->deviceIndex) == streams.end())
		return 0;
	return (long long)streams[arg->deviceIndex][idx + 1];
//#endif
//	return -1;
}

void CAuxGPU::SyncComputeStream(TGPUUsageArg* arg, long long targetStreamIdx, long long streamIdx) //HG24042025
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL)
		return;
	if (arg->deviceIndex == 0)
		return;
	if (targetStreamIdx < 0)
		return;
	if (streamIdx < 0)
		return;
	if (targetStreamIdx == streamIdx) //No need to sync the same stream
		return;

	cudaStream_t targetStream = (cudaStream_t)targetStreamIdx;
	cudaStream_t stream = (cudaStream_t)streamIdx;
	cudaEvent_t event;
	CUDA_SAFE(cudaEventCreateWithFlags(&event, cudaEventDisableTiming));
	CUDA_SAFE(cudaEventRecord(event, targetStream));
	CUDA_SAFE(cudaStreamWaitEvent(stream, event, 0));
	CUDA_SAFE(cudaEventDestroy(event));
//#endif
}

//void CAuxGPU::Init() 
void CAuxGPU::Init(TGPUUsageArg* arg) //HG02082024
{
	deviceOffloadInitialized = true;
//#ifdef _OFFLOAD_GPU
	if (arg == NULL) //HG02082024
		return;
	if (arg->deviceIndex <= 0)
		return;
	cudaGetDeviceCount(&deviceCount);
	if (arg->deviceIndex > deviceCount) //HG02082024
		return;
	if (streams.find(arg->deviceIndex) == streams.end())
	{
		//cudaInitDevice(arg->deviceIndex - 1, cudaDeviceMapHost, cudaInitDeviceFlagsAreValid);
		cudaSetDevice(arg->deviceIndex - 1);
		unsigned int flags = 0;
		cudaGetDeviceFlags(&flags);
		if (!(flags & cudaDeviceMapHost))cudaSetDeviceFlags(cudaDeviceMapHost);
		cudaStream_t *cur_streams = new cudaStream_t[10]; //HG23102025 Increase number of streams to 10 (1 memcpy + 9 compute)
		for (int i = 0; i < 10; i++) cudaStreamCreateWithFlags(&cur_streams[i], cudaStreamNonBlocking); //HG23102025 Create multiple streams
		streams[arg->deviceIndex] = cur_streams;
	}
	else
	{
		cudaSetDevice(arg->deviceIndex - 1);
		unsigned int flags = 0;
		cudaGetDeviceFlags(&flags);
		if (!(flags & cudaDeviceMapHost))cudaSetDeviceFlags(cudaDeviceMapHost);
	}
	memcpy_stream = streams[arg->deviceIndex][0]; //HG23102025 Move out of the if-block
	current_device = arg->deviceIndex;
	cudaDeviceSynchronize();
//#endif
}

//void CAuxGPU::Fini() 
void CAuxGPU::Fini(TGPUUsageArg* arg) //HG02082024
{
//#ifdef _OFFLOAD_GPU
	if (arg == NULL) //HG02082024
		return;
	if (arg->deviceIndex == 0)
		return;
	SetGPUStatus(false); //HG30112023 Disable GPU

	// Copy back all updated data
	bool updated = false;
	bool freed = false;
	for (std::map<void*, memAllocInfo_t>::const_iterator it = gpuMap.cbegin(); it != gpuMap.cend(); it++)
	{
		if (it->second.DevToHostUpdated){
			CUDA_SAFE(cudaStreamWaitEvent(memcpy_stream, it->second.d2h_event, 0));
			if (it->second.hostPtr != NULL) CUDA_SAFE(cudaMemcpyAsync(it->second.hostPtr, it->second.devicePtr, it->second.size, cudaMemcpyDeviceToHost, memcpy_stream));
//#if _DEBUG
//			printf("Fini: %p -> %p, %d\n", it->second.devicePtr, it->second.hostPtr, it->second.size);
//#endif
			updated = true;
			if (it->second.hostPtr != NULL) gpuMap[it->second.hostPtr].DevToHostUpdated = false;
			gpuMap[it->second.devicePtr].DevToHostUpdated = false;
		}
	}
	for (std::map<void*, memAllocInfo_t>::const_iterator it = gpuMap.cbegin(); it != gpuMap.cend(); it++)
	{
		if (it->first == it->second.devicePtr)
		{
			CUDA_SAFE(cudaStreamWaitEvent(0, it->second.h2d_event));
			CUDA_SAFE(cudaStreamWaitEvent(0, it->second.d2h_event));
			CUDA_SAFE(cudaFreeAsync(it->second.devicePtr, 0));
			freed = true;
			CUDA_SAFE(cudaEventDestroy(it->second.h2d_event));
			CUDA_SAFE(cudaEventDestroy(it->second.d2h_event));
		}
	}
	if (updated | freed)
		CUDA_SAFE(cudaStreamSynchronize(0));
	gpuMap.clear();
//#if _DEBUG
//	printf("Fini: %d\n", gpuMap.size());
//#endif
//#endif
}

#endif // _OFFLOAD_GPU
