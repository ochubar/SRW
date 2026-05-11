/************************************************************************//**
 * File: sroptfoc_gpu.cu
 * Description: Optical element: Focusing (CUDA implementation)
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
#include "sroptfoc.h"

//Implementation of the RadPointModifier's GPU function for the srTDriftSpace class
//int srTThinLens::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) 
int srTThinLens::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) //HG14042026
{ 
    //return RadPointModifierParallelImpl<srTThinLens>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU); 
    return TraverseRadZXEParallelImpl<srTThinLens>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU); //HG14042026
} //HG03092022
#endif