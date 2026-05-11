/************************************************************************//**
 * File: sroptzp_gpu.cu
 * Description: Optical element: Zone Plate (CUDA implementation)
 * Project: Synchrotron Radiation Workshop
 * First release: 2024
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
#include "sroptzp.h"

//Implementation of the RadPointModifier's GPU function for the srTRectAperture class
//int srTZonePlate::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) 
int srTZonePlate::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) //HG14042026
{ 
    //return RadPointModifierParallelImpl<srTZonePlate>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); 
    return TraverseRadZXEParallelImpl<srTZonePlate>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); //HG14042026
}

#endif