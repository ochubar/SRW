/************************************************************************//**
 * File: sroptapt_gpu.cu
 * Description: Optical element: Aperture (CUDA implementation)
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

#include "sroptapt.h"

//Implementation of the RadPointModifier's GPU function for the srTRectAperture class
//int srTRectAperture::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) 
int srTRectAperture::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg *pGPU) //HG14042026
{ 
    if (TransHndl.rep != 0)
        //return RadPointModifierParallelImpl<srTRectAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); 
        return TraverseRadZXEParallelImpl<srTRectAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); //HG14042026
    else
    {
        //Calculate the bounds of the aperture
        const double SmallOffset = 1.E-10;
        const int Margin = 5;

        double AptXStart = TransvCenPoint.x - HalfDx + SmallOffset;
        double AptXEnd = TransvCenPoint.x + HalfDx + SmallOffset;
        double AptZStart = TransvCenPoint.y - HalfDz + SmallOffset;
        double AptZEnd = TransvCenPoint.y + HalfDz + SmallOffset;

        long xStart = (long)floor((AptXStart - pRadAccessData->xStart) / pRadAccessData->xStep + Margin); //HG09122025
        long xEnd = (long)ceil((AptXEnd - pRadAccessData->xStart) / pRadAccessData->xStep - Margin);
        long zStart = (long)floor((AptZStart - pRadAccessData->zStart) / pRadAccessData->zStep + Margin);
        long zEnd = (long)ceil((AptZEnd - pRadAccessData->zStart) / pRadAccessData->zStep - Margin);
        //long xStart = floor((AptXStart - pRadAccessData->xStart) / pRadAccessData->xStep + Margin);
        //long xEnd = ceil((AptXEnd - pRadAccessData->xStart) / pRadAccessData->xStep - Margin);
        //long zStart = floor((AptZStart - pRadAccessData->zStart) / pRadAccessData->zStep + Margin);
        //long zEnd = ceil((AptZEnd - pRadAccessData->zStart) / pRadAccessData->zStep - Margin);

        xStart = min(max(xStart, 0l), pRadAccessData->nx - 1);
        xEnd = min(max(xEnd, 0l), pRadAccessData->nx - 1);
        zStart = min(max(zStart, 0l), pRadAccessData->nz - 1);
        zEnd = min(max(zEnd, 0l), pRadAccessData->nz - 1);

        if (xStart > xEnd) xStart = xEnd;
        if (zStart > zEnd) zStart = zEnd;

        if (xStart == 0 && xEnd == pRadAccessData->nx - 1 && zStart == 0 && zEnd == pRadAccessData->nz - 1)
        {
            return 0;
        }
        
        int region_params[5] = { (int)xStart, (int)xEnd, (int)zStart, (int)zEnd, 1 }; //HG09122025
        //int region_params[5] = { xStart, xEnd, zStart, zEnd, 1 };
        //return RadPointModifierParallelImpl<srTRectAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true);
        return TraverseRadZXEParallelImpl<srTRectAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true); //HG14042026
    }
}

//Implementation of the RadPointModifier's GPU function for the srTRectAperture class
//int srTRectObstacle::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU)
int srTRectObstacle::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU) //HG14042026
{
    if (TransHndl.rep != 0)
        //return RadPointModifierParallelImpl<srTRectObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU,0, true); 
        return TraverseRadZXEParallelImpl<srTRectObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU,0, true); //HG14042026
    else
    {
        //Calculate the bounds of the aperture
        const double SmallOffset = 1.E-10;
        const int Margin = 5;

        double AptXStart = TransvCenPoint.x - HalfDx + SmallOffset;
        double AptXEnd = TransvCenPoint.x + HalfDx + SmallOffset;
        double AptZStart = TransvCenPoint.y - HalfDz + SmallOffset;
        double AptZEnd = TransvCenPoint.y + HalfDz + SmallOffset;

        long xStart = (long)ceil((AptXStart - pRadAccessData->xStart) / pRadAccessData->xStep - Margin); //HG09122025
        long xEnd = (long)floor((AptXEnd - pRadAccessData->xStart) / pRadAccessData->xStep + Margin);
        long zStart = (long)ceil((AptZStart - pRadAccessData->zStart) / pRadAccessData->zStep - Margin);
        long zEnd = (long)floor((AptZEnd - pRadAccessData->zStart) / pRadAccessData->zStep + Margin);
        //long xStart = ceil((AptXStart - pRadAccessData->xStart) / pRadAccessData->xStep - Margin);
        //long xEnd = floor((AptXEnd - pRadAccessData->xStart) / pRadAccessData->xStep + Margin);
        //long zStart = ceil((AptZStart - pRadAccessData->zStart) / pRadAccessData->zStep - Margin);
        //long zEnd = floor((AptZEnd - pRadAccessData->zStart) / pRadAccessData->zStep + Margin);

        xStart = min(max(xStart, 0l), pRadAccessData->nx - 1);
        xEnd = min(max(xEnd, 0l), pRadAccessData->nx - 1);
        zStart = min(max(zStart, 0l), pRadAccessData->nz - 1);
        zEnd = min(max(zEnd, 0l), pRadAccessData->nz - 1);

		int region_params[5] = { (int)xStart, (int)xEnd, (int)zStart, (int)zEnd, 0 }; //HG09122025
        //int region_params[5] = { xStart, xEnd, zStart, zEnd, 0 };
        //return RadPointModifierParallelImpl<srTRectObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true);
        return TraverseRadZXEParallelImpl<srTRectObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true); //HG14042026
    }
}

//Implementation of the RadPointModifier's GPU function for the srTRectAperture class
//int srTCircAperture::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU)
int srTCircAperture::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU) //HG14042026
{
    if (TransHndl.rep != 0)
        //return RadPointModifierParallelImpl<srTCircAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); 
        return TraverseRadZXEParallelImpl<srTCircAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); //HG14042026
    else
    {
        //Calculate the bounds of the aperture
        const double SmallOffset = 1.E-10;
        const int Margin = 5;

        //Calculate the bounds of the aperture as the square fully embedded in the circle
        double Side = (1.414213562373095 * R)/2 + SmallOffset;

        int xStart = (int)(((TransvCenPoint.x-Side) - pRadAccessData->xStart) / pRadAccessData->xStep); //OC09122025
        int xEnd = (int)(((TransvCenPoint.x+Side)  - pRadAccessData->xStart) / pRadAccessData->xStep) + 1;
        int zStart = (int)(((TransvCenPoint.y-Side)  - pRadAccessData->zStart) / pRadAccessData->zStep);
        int zEnd = (int)(((TransvCenPoint.y+Side)  - pRadAccessData->zStart) / pRadAccessData->zStep) + 1;
        //int xStart = ((TransvCenPoint.x-Side) - pRadAccessData->xStart) / pRadAccessData->xStep;
        //int xEnd = ((TransvCenPoint.x+Side)  - pRadAccessData->xStart) / pRadAccessData->xStep + 1;
        //int zStart = ((TransvCenPoint.y-Side)  - pRadAccessData->zStart) / pRadAccessData->zStep;
        //int zEnd = ((TransvCenPoint.y+Side)  - pRadAccessData->zStart) / pRadAccessData->zStep + 1;

        if (xStart > 0) xStart += Margin; 
        if (xEnd < pRadAccessData->nx) xEnd -= Margin;
        if (zStart > 0) zStart += Margin; 
        if (zEnd < pRadAccessData->nz) zEnd -= Margin;

        if (xStart < 0) xStart = 0;
        if (xEnd > pRadAccessData->nx) xEnd = pRadAccessData->nx;
        if (zStart < 0) zStart = 0;
        if (zEnd > pRadAccessData->nz) zEnd = pRadAccessData->nz;

        if (xStart == 0 && xEnd == pRadAccessData->nx - 1 && zStart == 0 && zEnd == pRadAccessData->nz - 1)
        {
            return 0;
        }

        int region_params[5] = { xStart, xEnd, zStart, zEnd, 1 };
        //return RadPointModifierParallelImpl<srTCircAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true);
        return TraverseRadZXEParallelImpl<srTCircAperture>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true); //HG14042026
    }
}

//Implementation of the RadPointModifier's GPU function for the srTRectAperture class
//int srTCircObstacle::RadPointModifierParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU)
int srTCircObstacle::TraverseRadZXEParallel(srTSRWRadStructAccessData* pRadAccessData, void* pBufVars, long pBufVarsSz, TGPUUsageArg* pGPU) //HG14042026
{
    if (TransHndl.rep != 0)
        //return RadPointModifierParallelImpl<srTCircObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); 
        return TraverseRadZXEParallelImpl<srTCircObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, 0, true); //HG14042026
    else
    {
        //Calculate the bounds of the aperture
        const double SmallOffset = 1.E-10;
        const int Margin = 5;

        //Calculate the bounds of the aperture as the square fully embedded in the circle
        double Side = R + SmallOffset;

        int xStart = (int)(((TransvCenPoint.x-Side) - pRadAccessData->xStart) / pRadAccessData->xStep); //OC09122025
        int xEnd = (int)(((TransvCenPoint.x+Side)  - pRadAccessData->xStart) / pRadAccessData->xStep) + 1;
        int zStart = (int)(((TransvCenPoint.y-Side)  - pRadAccessData->zStart) / pRadAccessData->zStep);
        int zEnd = (int)(((TransvCenPoint.y+Side)  - pRadAccessData->zStart) / pRadAccessData->zStep) + 1;
        //int xStart = ((TransvCenPoint.x-Side) - pRadAccessData->xStart) / pRadAccessData->xStep;
        //int xEnd = ((TransvCenPoint.x+Side)  - pRadAccessData->xStart) / pRadAccessData->xStep + 1;
        //int zStart = ((TransvCenPoint.y-Side)  - pRadAccessData->zStart) / pRadAccessData->zStep;
        //int zEnd = ((TransvCenPoint.y+Side)  - pRadAccessData->zStart) / pRadAccessData->zStep + 1;

        xStart -= Margin; xEnd += Margin;
        zStart -= Margin; zEnd += Margin;
        if (xStart < 0) xStart = 0;
        if (xEnd > pRadAccessData->nx) xEnd = pRadAccessData->nx;
        if (zStart < 0) zStart = 0;
        if (zEnd > pRadAccessData->nz) zEnd = pRadAccessData->nz;

        int region_params[5] = { xStart, xEnd, zStart, zEnd, 0 };
        //return RadPointModifierParallelImpl<srTCircObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true);
        return TraverseRadZXEParallelImpl<srTCircObstacle>(pRadAccessData, pBufVars, pBufVarsSz, this, pGPU, region_params, true); //HG14042026
    }
}
#endif
