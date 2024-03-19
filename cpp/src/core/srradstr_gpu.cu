/************************************************************************//**
 * File: srradstr_gpu.cu
 * Description: Auxiliary structures for various SR calculation methods (CUDA implementation)
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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_constants.h"
#include <stdio.h>
#include <iostream>
#include <chrono>
#include "srradstr.h"

__global__ void MultiplyElFieldByPhaseLin_Kernel(double xMult, double zMult, float* pBaseRadX, float* pBaseRadZ, int nx, int nz, int ne, float xStart, float zStart, float xStep, float zStep) {
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range
    
    if (ix < nx && iz < nz) 
    {
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);

		double z = zStart + iz * zStep;
		double x = xStart + ix * xStep;
		double dPhZ = zMult * z;
		double dPh = dPhZ + xMult * x;
		double cosPh, sinPh;
		sincos(dPh, &sinPh, &cosPh);

		long long offset = iz * nx * ne * 2 + ix * ne * 2;
		float* tEx = pBaseRadX + offset;
		float* tEz = pBaseRadZ + offset;
		for (int ie = 0; ie < ne; ie++)
		{
			if (RadXisDefined)
			{
				//*(tEx++) *= a; *(tEx++) *= a;
				double newReEx = (*tEx) * cosPh - (*(tEx + 1)) * sinPh;
				double newImEx = (*tEx) * sinPh + (*(tEx + 1)) * cosPh;
				*(tEx++) = (float)newReEx; *(tEx++) = (float)newImEx;
			}
			if (RadZisDefined)
			{
				//*(tEz++) *= a; *(tEz++) *= a;
				double newReEz = (*tEz) * cosPh - (*(tEz + 1)) * sinPh;
				double newImEz = (*tEz) * sinPh + (*(tEz + 1)) * cosPh;
				*(tEz++) = (float)newReEz; *(tEz++) = (float)newImEz;
			}
		}
    }
}

void srTSRWRadStructAccessData::MultiplyElFieldByPhaseLin_GPU(double xMult, double zMult, TGPUUsageArg* pGPU) //OC03082023
//void srTSRWRadStructAccessData::MultiplyElFieldByPhaseLin_GPU(double xMult, double zMult, void* pGpuUsage)
{
	//TGPUUsageArg *pGpuUsage_ = (TGPUUsageArg*)pGpuUsage; //OC03082023 (commented-out)
	if (pBaseRadX != NULL)
	{
		pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, pBaseRadX, nz * nx * ne * 2 * sizeof(float)); //OC03082023
		//pBaseRadX = (float*)CAuxGPU::ToDevice(pGpuUsage_, pBaseRadX, nz * nx * ne * 2 * sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pBaseRadX); //OC03082023
		//CAuxGPU::EnsureDeviceMemoryReady(pGpuUsage_, pBaseRadX);
	}
	if (pBaseRadZ != NULL)
	{
		pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, pBaseRadZ, nz * nx * ne * 2 * sizeof(float)); //OC03082023
		//pBaseRadZ = (float*)CAuxGPU::ToDevice(pGpuUsage_, pBaseRadZ, nz * nx * ne * 2 * sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pBaseRadZ); //OC03082023
		//CAuxGPU::EnsureDeviceMemoryReady(pGpuUsage_, pBaseRadZ);
	}

    const int bs = 256;
    dim3 blocks(nx / bs + ((nx & (bs - 1)) != 0), nz);
    dim3 threads(bs, 1);
    MultiplyElFieldByPhaseLin_Kernel<<<blocks, threads>>> (xMult, zMult, pBaseRadX, pBaseRadZ, nx, nz, ne, (float)xStart, (float)zStart, (float)xStep, (float)zStep);
    //MultiplyElFieldByPhaseLin_Kernel<<<blocks, threads>>> (xMult, zMult, pBaseRadX, pBaseRadZ, nz, nx, ne, zStart, zStep, xStart, xStep);

	if (pBaseRadX != NULL)
		CAuxGPU::MarkUpdated(pGPU, pBaseRadX, true, false); //OC03082023
		//CAuxGPU::MarkUpdated(pGpuUsage_, pBaseRadX, true, false);
	if (pBaseRadZ != NULL)
		CAuxGPU::MarkUpdated(pGPU, pBaseRadZ, true, false); //OC03082023
		//CAuxGPU::MarkUpdated(pGpuUsage_, pBaseRadZ, true, false);

//HG26022024 (commented out)
//#ifdef _DEBUG
//	if (pBaseRadX != NULL)
//		pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, pBaseRadX, nz * nx * ne * 2 * sizeof(float)); //OC03082023
//		//pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGpuUsage_, pBaseRadX, nz * nx * ne * 2 * sizeof(float));
//	if (pBaseRadZ != NULL)
//		pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, pBaseRadZ, nz * nx * ne * 2 * sizeof(float)); //OC03082023
//		//pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGpuUsage_, pBaseRadZ, nz * nx * ne * 2 * sizeof(float));
//	cudaStreamSynchronize(0);
//    //auto err = cudaGetLastError();
//    //printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

template<int mode> __global__ void MirrorFieldData_Kernel(long nx, long nz, long ne, float* pEX0, float* pEZ0) {
	int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
	int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range

	if (ix < nx && iz < nz)
	{
		long long PerX = ne << 1;
		long long PerZ = PerX * nx;
		float buf;

		if (mode == 0)
		{
			if (ix >= (nx >> 1))
				return;

			long long nx_mi_1 = nx - 1; //OC26042019
			for (long long ie = 0; ie < ne; ie++)
			{
				//long Two_ie = ie << 1;
				long long Two_ie = ie << 1; //OC26042019
				
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz * PerZ;
				float* pEX_StartForX = pEX0 + izPerZ;
				float* pEZ_StartForX = pEZ0 + izPerZ;

				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix * PerX + Two_ie;
				float* pEX = pEX_StartForX + ixPerX_p_Two_ie;
				float* pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

				//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
				long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix) * PerX + Two_ie;
				float* rev_pEX = pEX_StartForX + rev_ixPerX_p_Two_ie;
				float* rev_pEZ = pEZ_StartForX + rev_ixPerX_p_Two_ie;

				if (pEX0 != 0)
				{
					buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
					buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
				}
				if (pEZ0 != 0)
				{
					buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
					buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
				}
			}
		}
		else if (mode == 1)
		{
			if (iz >= (nz >> 1))
				return;

			long long nz_mi_1 = nz - 1; //OC26042019
			for (long long ie = 0; ie < ne; ie++)
			{
				//long Two_ie = ie << 1;
				long long Two_ie = ie << 1;
				
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz * PerZ;
				float* pEX_StartForX = pEX0 + izPerZ;
				float* pEZ_StartForX = pEZ0 + izPerZ;

				//long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				long long rev_izPerZ = (nz_mi_1 - iz) * PerZ;
				float* rev_pEX_StartForX = pEX0 + rev_izPerZ;
				float* rev_pEZ_StartForX = pEZ0 + rev_izPerZ;

				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix * PerX + Two_ie;
				float* pEX = pEX_StartForX + ixPerX_p_Two_ie;
				float* pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

				float* rev_pEX = rev_pEX_StartForX + ixPerX_p_Two_ie;
				float* rev_pEZ = rev_pEZ_StartForX + ixPerX_p_Two_ie;

				if (pEX0 != 0)
				{
					buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
					buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
				}
				if (pEZ0 != 0)
				{
					buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
					buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
				}
			}
		}
		else if (mode == 2)
		{
			if (iz >= (nz >> 1))
				return;

			long long nx_mi_1 = nx - 1; //OC26042019
			long long nz_mi_1 = nz - 1;
			for (long long ie = 0; ie < ne; ie++) //OC26042019
				//for(long ie=0; ie<ne; ie++)
			{
				//long Two_ie = ie << 1;
				//for(long iz=0; iz<(nz >> 1); iz++)
				long long Two_ie = ie << 1; //OC26042019
				
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz * PerZ;
				float* pEX_StartForX = pEX0 + izPerZ;
				float* pEZ_StartForX = pEZ0 + izPerZ;

				//long rev_izPerZ = (nz_mi_1 - iz)*PerZ;
				long long rev_izPerZ = (nz_mi_1 - iz) * PerZ;
				float* rev_pEX_StartForX = pEX0 + rev_izPerZ;
				float* rev_pEZ_StartForX = pEZ0 + rev_izPerZ;

				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix * PerX + Two_ie;
				float* pEX = pEX_StartForX + ixPerX_p_Two_ie;
				float* pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

				//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
				long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix) * PerX + Two_ie;
				float* rev_pEX = rev_pEX_StartForX + rev_ixPerX_p_Two_ie;
				float* rev_pEZ = rev_pEZ_StartForX + rev_ixPerX_p_Two_ie;

				if (pEX0 != 0)
				{
					buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
					buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
				}
				if (pEZ0 != 0)
				{
					buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
					buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
				}

				if (((nz >> 1) << 1) != nz)
				{
					//long izPerZ = ((nz >> 1) + 1)*PerZ;
					long long izPerZ = ((nz >> 1) + 1) * PerZ;
					float* pEX_StartForX = pEX0 + izPerZ;
					float* pEZ_StartForX = pEZ0 + izPerZ;

					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix * PerX + Two_ie;
					float* pEX = pEX_StartForX + ixPerX_p_Two_ie;
					float* pEZ = pEZ_StartForX + ixPerX_p_Two_ie;

					//long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix)*PerX + Two_ie;
					long long rev_ixPerX_p_Two_ie = (nx_mi_1 - ix) * PerX + Two_ie;
					float* rev_pEX = pEX_StartForX + rev_ixPerX_p_Two_ie;
					float* rev_pEZ = pEZ_StartForX + rev_ixPerX_p_Two_ie;

					if (pEX0 != 0)
					{
						buf = *rev_pEX; *(rev_pEX++) = *pEX; *(pEX++) = buf;
						buf = *rev_pEX; *rev_pEX = *pEX; *pEX = buf;
					}
					if (pEZ0 != 0)
					{
						buf = *rev_pEZ; *(rev_pEZ++) = *pEZ; *(pEZ++) = buf;
						buf = *rev_pEZ; *rev_pEZ = *pEZ; *pEZ = buf;
					}
				}
			}
		}
	}
}

void srTSRWRadStructAccessData::MirrorFieldData_GPU(int sx, int sz, TGPUUsageArg* pGPU) //OC03082023
//void srTSRWRadStructAccessData::MirrorFieldData_GPU(int sx, int sz, void* pGpuUsage)
{
	//TGPUUsageArg *pGpuUsage_ = (TGPUUsageArg*)pGpuUsage; //OC03082023 (commented-out)
	float *pEX0 = pBaseRadX;
	float *pEZ0 = pBaseRadZ;

	if (pEX0 != NULL)
	{
		pEX0 = (float*)CAuxGPU::ToDevice(pGPU, pEX0, nz * nx * ne * 2 * sizeof(float)); //OC03082023
		//pEX0 = (float*)CAuxGPU::ToDevice(pGpuUsage_, pEX0, nz * nx * ne * 2 * sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pEX0); //OC03082023
		//CAuxGPU::EnsureDeviceMemoryReady(pGpuUsage_, pEX0);
	}
	if (pEZ0 != NULL)
	{
		pEZ0 = (float*)CAuxGPU::ToDevice(pGPU, pEZ0, nz * nx * ne * 2 * sizeof(float)); //OC03082023
		//pEZ0 = (float*)CAuxGPU::ToDevice(pGpuUsage_, pEZ0, nz * nx * ne * 2 * sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, pEZ0); //OC03082023
		//CAuxGPU::EnsureDeviceMemoryReady(pGpuUsage_, pEZ0);
	}

	const int bs = 256;
	dim3 blocks(nx / bs + ((nx & (bs - 1)) != 0), nz);
	dim3 threads(bs, 1);

	if ((sx > 0) && (sz > 0))
		return;
	else if ((sx < 0) && (sz > 0))
		MirrorFieldData_Kernel<0> <<<blocks, threads>>>(nx, nz, ne, pEX0, pEZ0);
	else if ((sx > 0) && (sz < 0))
		MirrorFieldData_Kernel<1> <<<blocks, threads >>> (nx, nz, ne, pEX0, pEZ0);
	else
		MirrorFieldData_Kernel<2> <<<blocks, threads >>> (nx, nz, ne, pEX0, pEZ0);

	if (pEX0 != NULL)
		CAuxGPU::MarkUpdated(pGPU, pEX0, true, false); //OC03082023
		//CAuxGPU::MarkUpdated(pGpuUsage_, pEX0, true, false);
	if (pEZ0 != NULL)
		CAuxGPU::MarkUpdated(pGPU, pEZ0, true, false); //OC03082023
		//CAuxGPU::MarkUpdated(pGpuUsage_, pEZ0, true, false);

//HG26022024 (commented out)
//#ifdef _DEBUG
//	if (pEX0 != NULL)
//		pEX0 = (float*)CAuxGPU::ToHostAndFree(pGPU, pEX0, nz * nx * ne * 2 * sizeof(float)); //OC03082023
//		//pEX0 = (float*)CAuxGPU::ToHostAndFree(pGpuUsage_, pEX0, nz * nx * ne * 2 * sizeof(float));
//	if (pEZ0 != NULL)
//		pEZ0 = (float*)CAuxGPU::ToHostAndFree(pGPU, pEZ0, nz * nx * ne * 2 * sizeof(float));
//		//pEZ0 = (float*)CAuxGPU::ToHostAndFree(pGpuUsage_, pEZ0, nz * nx * ne * 2 * sizeof(float));
//	cudaStreamSynchronize(0);
//	//auto err = cudaGetLastError();
//	//printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

#endif