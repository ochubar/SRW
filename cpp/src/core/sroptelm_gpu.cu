/************************************************************************//**
 * File: sroptelm_gpu.cu
 * Description: Optical element (general CUDA functions)
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
#include "sroptelm.h"
#include "sroptelm_gpu.h"


__global__ void TreatStronglyOscillatingTerm_Kernel(srTSRWRadStructAccessData RadAccessData, bool TreatPolCompX, bool TreatPolCompZ, double ConstRx, double ConstRz, int ieStart) {
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range
    int ie = (blockIdx.z * blockDim.z + threadIdx.z) + ieStart; //ne range
    
    if (ix < RadAccessData.nx && iz < RadAccessData.nz && ie < RadAccessData.ne + ieStart) 
    {
        double ePh = RadAccessData.eStart + RadAccessData.eStep * (ie - ieStart);
        if (RadAccessData.PresT == 1)
        {
            ePh = RadAccessData.avgPhotEn; //?? OC041108
        }

        double ConstRxE = ConstRx * ePh;
        double ConstRzE = ConstRz * ePh;
        if (RadAccessData.Pres == 1)
        {
            //double Lambda_m = 1.239854e-06/ePh;
            double Lambda_m = 1.239842e-06 / ePh;
            if (RadAccessData.PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

            double Lambda_me2 = Lambda_m * Lambda_m;
            ConstRxE *= Lambda_me2;
            ConstRzE *= Lambda_me2;
        }

        double z = (RadAccessData.zStart - RadAccessData.zc) + (iz * RadAccessData.zStep);
        double PhaseAddZ = 0;
        if (RadAccessData.WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE * z * z;

        double x = (RadAccessData.xStart - RadAccessData.xc) + (ix * RadAccessData.xStep);
        double Phase = PhaseAddZ;
        if (RadAccessData.WfrQuadTermCanBeTreatedAtResizeX) Phase += ConstRxE * x * x;

        float SinPh, CosPh;
        sincosf(Phase, &SinPh, &CosPh);

        long long PerX = RadAccessData.ne << 1;
        long long PerZ = PerX * RadAccessData.nx;
        long long offset = ie * 2 + iz * PerZ + ix * PerX;
        
		if (TreatPolCompX)
		{
			float* pExRe = RadAccessData.pBaseRadX + offset;
			float* pExIm = pExRe + 1;
			double ExReNew = (*pExRe) * CosPh - (*pExIm) * SinPh;
			double ExImNew = (*pExRe) * SinPh + (*pExIm) * CosPh;
			*pExRe = (float)ExReNew; *pExIm = (float)ExImNew;
		}
		if (TreatPolCompZ)
		{
			float* pEzRe = RadAccessData.pBaseRadZ + offset;
			float* pEzIm = pEzRe + 1;
			double EzReNew = (*pEzRe) * CosPh - (*pEzIm) * SinPh;
			double EzImNew = (*pEzRe) * SinPh + (*pEzIm) * CosPh;
			*pEzRe = (float)EzReNew; *pEzIm = (float)EzImNew;
		}
    }
}

void srTGenOptElem::TreatStronglyOscillatingTerm_GPU(srTSRWRadStructAccessData& RadAccessData, bool TreatPolCompX, bool TreatPolCompZ, double ConstRx, double ConstRz, int ieStart, int ieBefEnd, TGPUUsageArg* pGPU)
{
	if (RadAccessData.pBaseRadX != NULL)
	{
		RadAccessData.pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadX, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadX);
	}
	if (RadAccessData.pBaseRadZ != NULL)
	{
		RadAccessData.pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, RadAccessData.pBaseRadZ, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
		CAuxGPU::EnsureDeviceMemoryReady(pGPU, RadAccessData.pBaseRadZ);
	}

    const int bs = 256;
    dim3 blocks(RadAccessData.nx / bs + ((RadAccessData.nx & (bs - 1)) != 0), RadAccessData.nz, ieBefEnd - ieStart);
    dim3 threads(bs, 1);
    TreatStronglyOscillatingTerm_Kernel<< <blocks, threads >> > (RadAccessData, TreatPolCompX, TreatPolCompZ, ConstRx, ConstRz, ieStart);
	
	CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadX, true, false);
	CAuxGPU::MarkUpdated(pGPU, RadAccessData.pBaseRadZ, true, false);

#ifndef _DEBUG
	if (RadAccessData.pBaseRadX != NULL)
		RadAccessData.pBaseRadX = (float*)CAuxGPU::GetHostPtr(pGPU, RadAccessData.pBaseRadX);
	if (RadAccessData.pBaseRadZ != NULL)
		RadAccessData.pBaseRadZ = (float*)CAuxGPU::GetHostPtr(pGPU, RadAccessData.pBaseRadZ);
#endif

#ifdef _DEBUG
	if (RadAccessData.pBaseRadX != NULL)
		RadAccessData.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadX, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
	if (RadAccessData.pBaseRadZ != NULL)
		RadAccessData.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, RadAccessData.pBaseRadZ, 2*RadAccessData.ne*RadAccessData.nx*RadAccessData.nz*sizeof(float));
	cudaStreamSynchronize(0);
	auto err = cudaGetLastError();
	printf("%s\r\n", cudaGetErrorString(err));
#endif
}

__global__ void MakeWfrEdgeCorrection_Kernel(srTSRWRadStructAccessData RadAccessData, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr DataPtrs, float dxSt, float dxFi, float dzSt, float dzFi)
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y); //nz range

    if (ix < RadAccessData.nx && iz < RadAccessData.nz)
    {
		//float dxSt = (float)DataPtrs.dxSt;
		//float dxFi = (float)DataPtrs.dxFi;
		//float dzSt = (float)DataPtrs.dzSt;
		//float dzFi = (float)DataPtrs.dzFi;
		float dxSt_dzSt = dxSt * dzSt;
		float dxSt_dzFi = dxSt * dzFi;
		float dxFi_dzSt = dxFi * dzSt;
		float dxFi_dzFi = dxFi * dzFi;

		//long TwoNz = RadAccessData.nz << 1; //OC25012024 (commented-out)
		long PerX = 2;
		long PerZ = PerX * RadAccessData.nx;

        float fSSExRe = DataPtrs.fxStzSt[0];
        float fSSExIm = DataPtrs.fxStzSt[1];
        float fSSEzRe = DataPtrs.fxStzSt[2];
        float fSSEzIm = DataPtrs.fxStzSt[3];
        
        float fFSExRe = DataPtrs.fxFizSt[0];
        float fFSExIm = DataPtrs.fxFizSt[1];
        float fFSEzRe = DataPtrs.fxFizSt[2];
        float fFSEzIm = DataPtrs.fxFizSt[3];
        
        float fSFExRe = DataPtrs.fxStzFi[0];
        float fSFExIm = DataPtrs.fxStzFi[1];
        float fSFEzRe = DataPtrs.fxStzFi[2];
        float fSFEzIm = DataPtrs.fxStzFi[3];
        
        float fFFExRe = DataPtrs.fxFizFi[0];
        float fFFExIm = DataPtrs.fxFizFi[1];
        float fFFEzRe = DataPtrs.fxFizFi[2];
        float fFFEzIm = DataPtrs.fxFizFi[3];

		float bRe, bIm, cRe, cIm;

		long long Two_iz = iz << 1;
		long long Two_iz_p_1 = Two_iz + 1;
		long long Two_ix = ix << 1;
		long long Two_ix_p_1 = Two_ix + 1;

		float* tEx = pDataEx + iz * PerZ + ix * PerX, * tEz = pDataEz + iz * PerZ + ix * PerX;
		float ExRe = *tEx, ExIm = *(tEx + 1);
		float EzRe = *tEz, EzIm = *(tEz + 1);

		if (dxSt != 0.f)
		{
			float ExpXStRe = DataPtrs.ExpArrXSt[Two_ix], ExpXStIm = DataPtrs.ExpArrXSt[Two_ix_p_1];

			bRe = DataPtrs.FFTArrXStEx[Two_iz]; bIm = DataPtrs.FFTArrXStEx[Two_iz_p_1];
			ExRe += (float)(dxSt * (ExpXStRe * bRe - ExpXStIm * bIm));
			ExIm += (float)(dxSt * (ExpXStRe * bIm + ExpXStIm * bRe));

			bRe = DataPtrs.FFTArrXStEz[Two_iz]; bIm = DataPtrs.FFTArrXStEz[Two_iz_p_1];
			EzRe += (float)(dxSt * (ExpXStRe * bRe - ExpXStIm * bIm));
			EzIm += (float)(dxSt * (ExpXStRe * bIm + ExpXStIm * bRe));

			if (dzSt != 0.f)
			{
				bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
				cRe = ExpXStRe * bRe - ExpXStIm * bIm; cIm = ExpXStRe * bIm + ExpXStIm * bRe;

				ExRe += (float)(dxSt_dzSt * (fSSExRe * cRe - fSSExIm * cIm));
				ExIm += (float)(dxSt_dzSt * (fSSExRe * cIm + fSSExIm * cRe));
				EzRe += (float)(dxSt_dzSt * (fSSEzRe * cRe - fSSEzIm * cIm));
				EzIm += (float)(dxSt_dzSt * (fSSEzRe * cIm + fSSEzIm * cRe));
			}
			if (dzFi != 0.f)
			{
				bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
				cRe = ExpXStRe * bRe - ExpXStIm * bIm; cIm = ExpXStRe * bIm + ExpXStIm * bRe;

				ExRe -= (float)(dxSt_dzFi * (fSFExRe * cRe - fSFExIm * cIm));
				ExIm -= (float)(dxSt_dzFi * (fSFExRe * cIm + fSFExIm * cRe));
				EzRe -= (float)(dxSt_dzFi * (fSFEzRe * cRe - fSFEzIm * cIm));
				EzIm -= (float)(dxSt_dzFi * (fSFEzRe * cIm + fSFEzIm * cRe));
			}
		}
		if (dxFi != 0.f)
		{
			float ExpXFiRe = DataPtrs.ExpArrXFi[Two_ix], ExpXFiIm = DataPtrs.ExpArrXFi[Two_ix_p_1];

			bRe = DataPtrs.FFTArrXFiEx[Two_iz]; bIm = DataPtrs.FFTArrXFiEx[Two_iz_p_1];
			ExRe -= (float)(dxFi * (ExpXFiRe * bRe - ExpXFiIm * bIm));
			ExIm -= (float)(dxFi * (ExpXFiRe * bIm + ExpXFiIm * bRe));

			bRe = DataPtrs.FFTArrXFiEz[Two_iz]; bIm = DataPtrs.FFTArrXFiEz[Two_iz_p_1];
			EzRe -= (float)(dxFi * (ExpXFiRe * bRe - ExpXFiIm * bIm));
			EzIm -= (float)(dxFi * (ExpXFiRe * bIm + ExpXFiIm * bRe));

			if (dzSt != 0.f)
			{
				bRe = DataPtrs.ExpArrZSt[Two_iz], bIm = DataPtrs.ExpArrZSt[Two_iz_p_1];
				cRe = ExpXFiRe * bRe - ExpXFiIm * bIm; cIm = ExpXFiRe * bIm + ExpXFiIm * bRe;

				ExRe -= (float)(dxFi_dzSt * (fFSExRe * cRe - fFSExIm * cIm));
				ExIm -= (float)(dxFi_dzSt * (fFSExRe * cIm + fFSExIm * cRe));
				EzRe -= (float)(dxFi_dzSt * (fFSEzRe * cRe - fFSEzIm * cIm));
				EzIm -= (float)(dxFi_dzSt * (fFSEzRe * cIm + fFSEzIm * cRe));
			}
			if (dzFi != 0.f)
			{
				bRe = DataPtrs.ExpArrZFi[Two_iz], bIm = DataPtrs.ExpArrZFi[Two_iz_p_1];
				cRe = ExpXFiRe * bRe - ExpXFiIm * bIm; cIm = ExpXFiRe * bIm + ExpXFiIm * bRe;

				ExRe += (float)(dxFi_dzFi * (fFFExRe * cRe - fFFExIm * cIm));
				ExIm += (float)(dxFi_dzFi * (fFFExRe * cIm + fFFExIm * cRe));
				EzRe += (float)(dxFi_dzFi * (fFFEzRe * cRe - fFFEzIm * cIm));
				EzIm += (float)(dxFi_dzFi * (fFFEzRe * cIm + fFFEzIm * cRe));
			}
		}
		if (dzSt != 0.f)
		{
			float ExpZStRe = DataPtrs.ExpArrZSt[Two_iz], ExpZStIm = DataPtrs.ExpArrZSt[Two_iz_p_1];

			bRe = DataPtrs.FFTArrZStEx[Two_ix]; bIm = DataPtrs.FFTArrZStEx[Two_ix_p_1];
			ExRe += (float)(dzSt * (ExpZStRe * bRe - ExpZStIm * bIm));
			ExIm += (float)(dzSt * (ExpZStRe * bIm + ExpZStIm * bRe));

			bRe = DataPtrs.FFTArrZStEz[Two_ix]; bIm = DataPtrs.FFTArrZStEz[Two_ix_p_1];
			EzRe += (float)(DataPtrs.dzSt * (ExpZStRe * bRe - ExpZStIm * bIm));
			EzIm += (float)(DataPtrs.dzSt * (ExpZStRe * bIm + ExpZStIm * bRe));
		}
		if (dzFi != 0.f)
		{
			float ExpZFiRe = DataPtrs.ExpArrZFi[Two_iz], ExpZFiIm = DataPtrs.ExpArrZFi[Two_iz_p_1];

			bRe = DataPtrs.FFTArrZFiEx[Two_ix]; bIm = DataPtrs.FFTArrZFiEx[Two_ix_p_1];
			ExRe -= (float)(dzFi * (ExpZFiRe * bRe - ExpZFiIm * bIm));
			ExIm -= (float)(dzFi * (ExpZFiRe * bIm + ExpZFiIm * bRe));

			bRe = DataPtrs.FFTArrZFiEz[Two_ix]; bIm = DataPtrs.FFTArrZFiEz[Two_ix_p_1];
			EzRe -= (float)(dzFi * (ExpZFiRe * bRe - ExpZFiIm * bIm));
			EzIm -= (float)(dzFi * (ExpZFiRe * bIm + ExpZFiIm * bRe));
		}

		*tEx = ExRe; *(tEx + 1) = ExIm;
		*tEz = EzRe; *(tEz + 1) = EzIm;
    }
}

void srTGenOptElem::MakeWfrEdgeCorrection_GPU(srTSRWRadStructAccessData* RadAccessData, float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrs, TGPUUsageArg* pGPU)
{
	pDataEx = (float*)CAuxGPU::ToDevice(pGPU, pDataEx, 2*RadAccessData->ne*RadAccessData->nx*RadAccessData->nz*sizeof(float));
	pDataEz = (float*)CAuxGPU::ToDevice(pGPU, pDataEz, 2*RadAccessData->ne*RadAccessData->nx*RadAccessData->nz*sizeof(float));
	DataPtrs.FFTArrXStEx = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrXStEx, 2*RadAccessData->nz*sizeof(float));
	DataPtrs.FFTArrXStEz = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrXStEz, 2*RadAccessData->nz*sizeof(float));
	DataPtrs.FFTArrXFiEx = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrXFiEx, 2*RadAccessData->nz*sizeof(float));
	DataPtrs.FFTArrXFiEz = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrXFiEz, 2*RadAccessData->nz*sizeof(float));
	DataPtrs.FFTArrZStEx = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrZStEx, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.FFTArrZStEz = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrZStEz, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.FFTArrZFiEx = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrZFiEx, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.FFTArrZFiEz = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.FFTArrZFiEz, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.ExpArrXSt = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.ExpArrXSt, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.ExpArrXFi = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.ExpArrXFi, 2*RadAccessData->nx*sizeof(float));
	DataPtrs.ExpArrZSt = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.ExpArrZSt, 2*RadAccessData->nz*sizeof(float));
	DataPtrs.ExpArrZFi = (float*)CAuxGPU::ToDevice(pGPU, DataPtrs.ExpArrZFi, 2*RadAccessData->nz*sizeof(float));

	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pDataEx);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, pDataEz);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrXStEx);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrXStEz);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrXFiEx);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrXFiEz);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrZStEx);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrZStEz);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrZFiEx);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.FFTArrZFiEz);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.ExpArrXSt);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.ExpArrXFi);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.ExpArrZSt);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, DataPtrs.ExpArrZFi);

	const int bs = 256;
	dim3 blocks(RadAccessData->nx / bs + ((RadAccessData->nx & (bs - 1)) != 0), RadAccessData->nz);
	dim3 threads(bs, 1);
	MakeWfrEdgeCorrection_Kernel << <blocks, threads >> > (*RadAccessData, pDataEx, pDataEz, DataPtrs, (float)DataPtrs.dxSt, (float)DataPtrs.dxFi, (float)DataPtrs.dzSt, (float)DataPtrs.dzFi);

	DataPtrs.FFTArrXStEx = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrXStEx, 2*RadAccessData->nz*sizeof(float), true);
	DataPtrs.FFTArrXStEz = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrXStEz, 2*RadAccessData->nz*sizeof(float), true);
	DataPtrs.FFTArrXFiEx = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrXFiEx, 2*RadAccessData->nz*sizeof(float), true);
	DataPtrs.FFTArrXFiEz = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrXFiEz, 2*RadAccessData->nz*sizeof(float), true);
	DataPtrs.FFTArrZStEx = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrZStEx, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.FFTArrZStEz = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrZStEz, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.FFTArrZFiEx = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrZFiEx, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.FFTArrZFiEz = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.FFTArrZFiEz, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.ExpArrXSt = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.ExpArrXSt, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.ExpArrXFi = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.ExpArrXFi, 2*RadAccessData->nx*sizeof(float), true);
	DataPtrs.ExpArrZSt = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.ExpArrZSt, 2*RadAccessData->nz*sizeof(float), true);
	DataPtrs.ExpArrZFi = (float*)CAuxGPU::ToHostAndFree(pGPU, DataPtrs.ExpArrZFi, 2*RadAccessData->nz*sizeof(float), true);

	CAuxGPU::MarkUpdated(pGPU, pDataEx, true, false);
	CAuxGPU::MarkUpdated(pGPU, pDataEz, true, false);

#ifdef _DEBUG
	CAuxGPU::ToHostAndFree(pGPU, pDataEx, 2*RadAccessData->ne*RadAccessData->nx*RadAccessData->nz*sizeof(float));
	CAuxGPU::ToHostAndFree(pGPU, pDataEz, 2*RadAccessData->ne*RadAccessData->nx*RadAccessData->nz*sizeof(float));
	cudaStreamSynchronize(0);
	auto err = cudaGetLastError();
	printf("%s\r\n", cudaGetErrorString(err));
#endif
}

template<bool TreatPolCompX, bool TreatPolCompZ> __global__ void RadResizeCore_Kernel(srTSRWRadStructAccessData OldRadAccessData, srTSRWRadStructAccessData NewRadAccessData)
{
	int ixStart = int(NewRadAccessData.AuxLong1);
	int ixEnd = int(NewRadAccessData.AuxLong2);
	int izStart = int(NewRadAccessData.AuxLong3);
	int izEnd = int(NewRadAccessData.AuxLong4);

    int ix = (blockIdx.x * blockDim.x + threadIdx.x) + ixStart; //nx range
    int iz = (blockIdx.y * blockDim.y + threadIdx.y) + izStart; //nz range
    int ie = (blockIdx.z * blockDim.z + threadIdx.z); //ne range

	if (ix > ixEnd) return;
	if (iz > izEnd) return;

	const double DistAbsTol = 1.E-10;
	double xStepInvOld = 1./OldRadAccessData.xStep;
	double zStepInvOld = 1./OldRadAccessData.zStep;
	int nx_mi_1Old = OldRadAccessData.nx - 1;
	int nz_mi_1Old = OldRadAccessData.nz - 1;
	int nx_mi_2Old = nx_mi_1Old - 1;
	int nz_mi_2Old = nz_mi_1Old - 1;

	//OC31102018: moved by SY at parallelizing SRW via OpenMP
	//srTInterpolAux01 InterpolAux01;
	//srTInterpolAux02 InterpolAux02[4], InterpolAux02I[2];
	//srTInterpolAuxF AuxF[4], AuxFI[2];
	//int ixStOld, izStOld, ixStOldPrev = -1000, izStOldPrev = -1000;

	//long PerX_New = NewRadAccessData.ne << 1;
	//long PerZ_New = PerX_New*NewRadAccessData.nx;
	long long PerX_New = NewRadAccessData.ne << 1;
	long long PerZ_New = PerX_New*NewRadAccessData.nx;

	//long PerX_Old = PerX_New;
	//long PerZ_Old = PerX_Old*OldRadAccessData.nx;
	long long PerX_Old = PerX_New;
	long long PerZ_Old = PerX_Old*OldRadAccessData.nx;

	float *pEX0_New = 0, *pEZ0_New = 0;
	pEX0_New = NewRadAccessData.pBaseRadX;
	pEZ0_New = NewRadAccessData.pBaseRadZ;

	float* pEX0_Old = 0, * pEZ0_Old = 0;
	pEX0_Old = OldRadAccessData.pBaseRadX;
	pEZ0_Old = OldRadAccessData.pBaseRadZ;

	
	int ixStOld, izStOld; //OC25012024 //ixStOldPrev = -1000, izStOldPrev = -1000;
	//int ixStOld, izStOld, ixStOldPrev = -1000, izStOldPrev = -1000;
	//SY: do we need this (always returns 0, updates some clock)
	//if(result = srYield.Check()) return result;

	double zAbs = NewRadAccessData.zStart + iz * NewRadAccessData.zStep;

	char FieldShouldBeZeroedDueToZ = 0;
	if (NewRadAccessData.WfrEdgeCorrShouldBeDone)
	{
		if ((zAbs < NewRadAccessData.zWfrMin - DistAbsTol) || (zAbs > NewRadAccessData.zWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToZ = 1;
	}

	int izcOld = int((zAbs - OldRadAccessData.zStart) * zStepInvOld + 1.E-06);

	double zRel = zAbs - (OldRadAccessData.zStart + izcOld * OldRadAccessData.zStep);

	if (izcOld == nz_mi_1Old) { izStOld = izcOld - 3; zRel += 2. * OldRadAccessData.zStep; }
	else if (izcOld == nz_mi_2Old) { izStOld = izcOld - 2; zRel += OldRadAccessData.zStep; }
	else if (izcOld == 0) { izStOld = izcOld; zRel -= OldRadAccessData.zStep; }
	else izStOld = izcOld - 1;

	zRel *= zStepInvOld;

	int izcOld_mi_izStOld = izcOld - izStOld;
	//long izPerZ_New = iz*PerZ_New;
	long long izPerZ_New = iz * PerZ_New;

	double xAbs = NewRadAccessData.xStart + ix * NewRadAccessData.xStep;

	char FieldShouldBeZeroedDueToX = 0;
	if (NewRadAccessData.WfrEdgeCorrShouldBeDone)
	{
		if ((xAbs < NewRadAccessData.xWfrMin - DistAbsTol) || (xAbs > NewRadAccessData.xWfrMax + DistAbsTol)) FieldShouldBeZeroedDueToX = 1;
	}
	char FieldShouldBeZeroed = (FieldShouldBeZeroedDueToX || FieldShouldBeZeroedDueToZ);

	int ixcOld = int((xAbs - OldRadAccessData.xStart) * xStepInvOld + 1.E-06);
	double xRel = xAbs - (OldRadAccessData.xStart + ixcOld * OldRadAccessData.xStep);

	if (ixcOld == nx_mi_1Old) { ixStOld = ixcOld - 3; xRel += 2. * OldRadAccessData.xStep; }
	else if (ixcOld == nx_mi_2Old) { ixStOld = ixcOld - 2; xRel += OldRadAccessData.xStep; }
	else if (ixcOld == 0) { ixStOld = ixcOld; xRel -= OldRadAccessData.xStep; }
	else ixStOld = ixcOld - 1;

	xRel *= xStepInvOld;

	int ixcOld_mi_ixStOld = ixcOld - ixStOld;

	//or (int ie = 0; ie < NewRadAccessData.ne; ie++)
	{
		//OC31102018: modified by SY at OpenMP parallelization
		//ixStOldPrev = -1000; izStOldPrev = -1000;

		//OC31102018: moved by SY at OpenMP parallelization
		srTInterpolAux01 InterpolAux01;
		srTInterpolAux02 InterpolAux02[4], InterpolAux02I[2];
		srTInterpolAuxF AuxF[4], AuxFI[2];
		//ixStOldPrev = -1000; izStOldPrev = -1000; //OC25012024 (commented-out: never used?)
		float BufF[4], BufFI[2];
		char UseLowOrderInterp_PolCompX = 0, UseLowOrderInterp_PolCompZ = 0;

		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;

		float* pEX_StartForX_New = 0, * pEZ_StartForX_New = 0;
		pEX_StartForX_New = pEX0_New + izPerZ_New;
		pEZ_StartForX_New = pEZ0_New + izPerZ_New;

		//long ixPerX_New_p_Two_ie = ix*PerX_New + Two_ie;
		long long ixPerX_New_p_Two_ie = ix * PerX_New + Two_ie;
		float* pEX_New = 0, * pEZ_New = 0;
		pEX_New = pEX_StartForX_New + ixPerX_New_p_Two_ie;
		pEZ_New = pEZ_StartForX_New + ixPerX_New_p_Two_ie;

		//long TotOffsetOld = izStOld*PerZ_Old + ixStOld*PerX_Old + Two_ie;
		long long TotOffsetOld = izStOld * PerZ_Old + ixStOld * PerX_Old + Two_ie;

		if (TreatPolCompX)
		{
			float* pExSt_Old = pEX0_Old + TotOffsetOld;
			srTGenOptElem::GetCellDataForInterpol(pExSt_Old, PerX_Old, PerZ_Old, AuxF);

			srTGenOptElem::SetupCellDataI(AuxF, AuxFI);
			UseLowOrderInterp_PolCompX = srTGenOptElem::CheckForLowOrderInterp(AuxF, AuxFI, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02, InterpolAux02I);

			if (!UseLowOrderInterp_PolCompX)
			{
				for (int i = 0; i < 2; i++)
				{
					srTGenOptElem::SetupInterpolAux02(AuxF + i, &InterpolAux01, InterpolAux02 + i);
				}
				srTGenOptElem::SetupInterpolAux02(AuxFI, &InterpolAux01, InterpolAux02I);
			}

			if (UseLowOrderInterp_PolCompX)
			{
				srTGenOptElem::InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 0);
				srTGenOptElem::InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 0);
			}
			else
			{
				srTGenOptElem::InterpolF(InterpolAux02, xRel, zRel, BufF, 0);
				srTGenOptElem::InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 0);
			}

			(*BufFI) *= AuxFI->fNorm;
			srTGenOptElem::ImproveReAndIm(BufF, BufFI);

			if (FieldShouldBeZeroed)
			{
				*BufF = 0.; *(BufF + 1) = 0.;
			}

			*pEX_New = *BufF;
			*(pEX_New + 1) = *(BufF + 1);
		}
		if (TreatPolCompZ)
		{
			float* pEzSt_Old = pEZ0_Old + TotOffsetOld;
			srTGenOptElem::GetCellDataForInterpol(pEzSt_Old, PerX_Old, PerZ_Old, AuxF + 2);

			srTGenOptElem::SetupCellDataI(AuxF + 2, AuxFI + 1);
			UseLowOrderInterp_PolCompZ = srTGenOptElem::CheckForLowOrderInterp(AuxF + 2, AuxFI + 1, ixcOld_mi_ixStOld, izcOld_mi_izStOld, &InterpolAux01, InterpolAux02 + 2, InterpolAux02I + 1);

			if (!UseLowOrderInterp_PolCompZ)
			{
				for (int i = 0; i < 2; i++)
				{
					srTGenOptElem::SetupInterpolAux02(AuxF + 2 + i, &InterpolAux01, InterpolAux02 + 2 + i);
				}
				srTGenOptElem::SetupInterpolAux02(AuxFI + 1, &InterpolAux01, InterpolAux02I + 1);
			}
			
			if (UseLowOrderInterp_PolCompZ)
			{
				srTGenOptElem::InterpolF_LowOrder(InterpolAux02, xRel, zRel, BufF, 2);
				srTGenOptElem::InterpolFI_LowOrder(InterpolAux02I, xRel, zRel, BufFI, 1);
			}
			else
			{
				srTGenOptElem::InterpolF(InterpolAux02, xRel, zRel, BufF, 2);
				srTGenOptElem::InterpolFI(InterpolAux02I, xRel, zRel, BufFI, 1);
			}

			(*(BufFI + 1)) *= (AuxFI + 1)->fNorm;
			srTGenOptElem::ImproveReAndIm(BufF + 2, BufFI + 1);

			if (FieldShouldBeZeroed)
			{
				*(BufF + 2) = 0.; *(BufF + 3) = 0.;
			}

			*pEZ_New = *(BufF + 2);
			*(pEZ_New + 1) = *(BufF + 3);
		}
	}
}

int srTGenOptElem::RadResizeCore_GPU(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, char PolComp, TGPUUsageArg* pGPU)
{
	char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x'));
	char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z'));

	int nx = NewRadAccessData.AuxLong2 - NewRadAccessData.AuxLong1 + 1;
	int nz = NewRadAccessData.AuxLong4 - NewRadAccessData.AuxLong3 + 1;
	int ne = NewRadAccessData.ne;
	OldRadAccessData.pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, OldRadAccessData.pBaseRadX, 2*OldRadAccessData.ne*OldRadAccessData.nx*OldRadAccessData.nz*sizeof(float));
	OldRadAccessData.pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, OldRadAccessData.pBaseRadZ, 2*OldRadAccessData.ne*OldRadAccessData.nx*OldRadAccessData.nz*sizeof(float));
	NewRadAccessData.pBaseRadX = (float*)CAuxGPU::ToDevice(pGPU, NewRadAccessData.pBaseRadX, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float), true);
	NewRadAccessData.pBaseRadZ = (float*)CAuxGPU::ToDevice(pGPU, NewRadAccessData.pBaseRadZ, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float), true);
	
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, OldRadAccessData.pBaseRadX);
	CAuxGPU::EnsureDeviceMemoryReady(pGPU, OldRadAccessData.pBaseRadZ);
	//CAuxGPU::EnsureDeviceMemoryReady(pGPU, NewRadAccessData.pBaseRadX);
	//CAuxGPU::EnsureDeviceMemoryReady(pGPU, NewRadAccessData.pBaseRadZ);

	const int bs = 32;
	dim3 blocks(nx / bs + ((nx & (bs - 1)) != 0), nz, ne);
	dim3 threads(bs, 1);
	
	if (TreatPolCompX && TreatPolCompZ) RadResizeCore_Kernel<true, true> << <blocks, threads >> > (OldRadAccessData, NewRadAccessData);
	else if (TreatPolCompX) RadResizeCore_Kernel<true, false> << <blocks, threads >> > (OldRadAccessData, NewRadAccessData);
	else if (TreatPolCompZ) RadResizeCore_Kernel<false, true> << <blocks, threads >> > (OldRadAccessData, NewRadAccessData);

	OldRadAccessData.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, OldRadAccessData.pBaseRadX, 2*OldRadAccessData.ne*OldRadAccessData.nx*OldRadAccessData.nz*sizeof(float), true);
	OldRadAccessData.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, OldRadAccessData.pBaseRadZ, 2*OldRadAccessData.ne*OldRadAccessData.nx*OldRadAccessData.nz*sizeof(float), true);
	//NewRadAccessData.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, NewRadAccessData.pBaseRadX, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float));
	//NewRadAccessData.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, NewRadAccessData.pBaseRadZ, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float));
	CAuxGPU::MarkUpdated(pGPU, NewRadAccessData.pBaseRadX, true, false);
	CAuxGPU::MarkUpdated(pGPU, NewRadAccessData.pBaseRadZ, true, false);
#ifndef _DEBUG
	NewRadAccessData.pBaseRadX = (float*)CAuxGPU::GetHostPtr(pGPU, NewRadAccessData.pBaseRadX);
	NewRadAccessData.pBaseRadZ = (float*)CAuxGPU::GetHostPtr(pGPU, NewRadAccessData.pBaseRadZ);
#endif

#ifdef _DEBUG
	NewRadAccessData.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, NewRadAccessData.pBaseRadX, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float), false);
	NewRadAccessData.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, NewRadAccessData.pBaseRadZ, 2*NewRadAccessData.ne*NewRadAccessData.nx*NewRadAccessData.nz*sizeof(float), false);
	cudaStreamSynchronize(0);
	auto err = cudaGetLastError();
	printf("%s\r\n", cudaGetErrorString(err));

#endif

	return 0;
}

#endif