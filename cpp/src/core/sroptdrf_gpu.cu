/************************************************************************//**
 * File: sroptdrf_gpu.cu
 * Description: Optical element: Drift space (CUDA implementation)
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
#include "sroptdrf.h"

//Implementation of the RadPointModifier's GPU function for the srTDriftSpace class
int srTDriftSpace::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) 
{ 
    return RadPointModifierParallelImpl<srTDriftSpace>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU); 
} //HG03092022
#endif