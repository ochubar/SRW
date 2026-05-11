/************************************************************************//**
 * File: sroptgtr_gpu.cu
 * Description: Optical element: Transmission (CUDA implementation)
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
#include "sroptgtr.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_constants.h"

#include <stdio.h>
#include <iostream>
#include <chrono>

//int srTGenTransmission::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU)
int srTGenTransmission::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU) //HG14042026
{
    //GenTransNumData.pData = (char*)CAuxGPU::ToDevice(pGPU, GenTransNumData.pData, GenTransNumData.DimSizes[0] * (int)GenTransNumData.DimSizes[1] * (int)GenTransNumData.DimSizes[2] * sizeof(double) * 2);
    GenTransNumData.pData = CAuxGPU::ToDevice(pGPU, GenTransNumData.pData, GenTransNumData.DimSizes[0] * (int)GenTransNumData.DimSizes[1] * (int)GenTransNumData.DimSizes[2] * 2 * sizeof(double), CAuxGPU::DISCARD_HOST); //HG01082025
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, GenTransNumData.pData);
    //int retCode = RadPointModifierParallelImpl<srTGenTransmission>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU);
    int retCode = TraverseRadZXEParallelImpl<srTGenTransmission>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU); //HG14042026
    CAuxGPU::MarkUpdated(pGPU, GenTransNumData.pData, CAuxGPU::DEVICE); //HG23102025
    //GenTransNumData.pData = (char*)CAuxGPU::ToHostAndFree(pGPU, GenTransNumData.pData, GenTransNumData.DimSizes[0] * (int)GenTransNumData.DimSizes[1] * (int)GenTransNumData.DimSizes[2] * sizeof(double) * 2, true);
	GenTransNumData.pData = CAuxGPU::ToHostAndFree(pGPU, GenTransNumData.pData); //HG01082025
    return retCode;
} //HG03092022
#endif