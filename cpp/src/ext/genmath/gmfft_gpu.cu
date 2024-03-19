/************************************************************************//**
 * File: gmfft_gpu.cu
 * Description: Auxiliary utilities to work with FFTW library (CUDA implementation)
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
#include "gmfft.h"

#define GMFFT_BLOCK_SIZE 256

template <typename T> __global__ void RepairSignAfter1DFFT_Kernel(T* pAfterFFT, long HowMany, long Nx2) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 4 + 2; //Nx range
    
    if (ix < Nx2) 
    {
        for (long k = 0; k < HowMany; k++)
        {
            pAfterFFT[ix + k * Nx2] = -pAfterFFT[ix + k * Nx2];
            pAfterFFT[ix + k * Nx2 + 1] = -pAfterFFT[ix + k * Nx2 + 1];
        }
    }
}

template <typename T> __global__ void RotateDataAfter1DFFT_Kernel(T* pAfterFFT, long HowMany, long Nx2, long Nx) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 2; //HalfNx range
    
    if (ix < Nx) 
    {
        for (long k = 0; k < HowMany; k++)
        {
            T t1_0 = pAfterFFT[ix + Nx2 * k];
            T t1_1 = pAfterFFT[ix + Nx2 * k + 1];

            pAfterFFT[ix + Nx2 * k] = pAfterFFT[ix + Nx + Nx2 * k];
            pAfterFFT[ix + Nx2 * k + 1] = pAfterFFT[ix + Nx + Nx2 * k + 1];
            pAfterFFT[ix + Nx + Nx2 * k] = t1_0;
            pAfterFFT[ix + Nx + Nx2 * k + 1] = t1_1;
        }
    }
}

template <typename T> __global__ void RepairAndRotateAfter1DFFT_Kernel(T* pAfterFFT, long HowMany, long Nx, float Mult) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //HalfNx range
    
    long HalfNx = Nx / 2;
    long Nx2 = Nx * 2;
    if (ix < HalfNx) 
    {
        float sx0 = 1 - 2 * (ix % 2);
        float sx1 = 1 - 2 * ((HalfNx + ix) % 2);
        
        float s1 = sx0 * Mult;
        float s2 = sx1 * Mult;

        int idx = ix * 2;
        for (long i = 0; i < HowMany; i++){
            T* t1 = pAfterFFT + i * Nx2, *t2 = pAfterFFT + (HalfNx) * 2 + i * Nx2;
            
            T buf_r = t1[idx] * s1;
            T buf_im = t1[idx + 1] * s1;

            t1[idx] = t2[idx] * s2;
            t1[idx + 1] = t2[idx + 1] * s2;

            t2[idx] = buf_r;
            t2[idx + 1] = buf_im;
        }
    }
}

template <typename T> __global__ void NormalizeDataAfter1DFFT_Kernel(T* pAfterFFT, long HowMany, long Nx2, T Mult) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 2; //Nx range
    
    if (ix < Nx2) 
    {
        for (long i = 0; i < HowMany; i++) {
            pAfterFFT[ix + i * Nx2] *= Mult;
            pAfterFFT[ix + i * Nx2 + 1] *= Mult;
        }
    }
}

template <typename T> __global__ void FillArrayShift_Kernel(double t0, double tStep, long N, T* arShiftX) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //HalfNx range

    double t0TwoPi = t0 * 2 * CUDART_PI;
    double q = tStep * ix;

    if (ix < N) 
    {
        if (ix == 0) {
            arShiftX[N] = 1.0;
            arShiftX[N + 1] = 0.0;
        }

        ix *= 2;
        if (ix < N - 2) 
        {
            sincos(q * t0TwoPi, &arShiftX[N + 2 + 1 + ix], &arShiftX[N + 2 + ix]);
            arShiftX[N - 2 - ix] = arShiftX[N + 2 + ix];
            arShiftX[N - 1 - ix] = -arShiftX[N + 2 + 1 + ix];
        }

        if (ix == N - 2) 
        {
            sincos(-q * t0TwoPi, &arShiftX[1], &arShiftX[0]);
        }
    }
}

template <typename T> __global__ void TreatShift_Kernel(T* pData, long HowMany, long Nx2, T* tShiftX) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 2; //Nx range
    
    if (ix < Nx2) 
    {
        T MultX_Re = tShiftX[ix];
        T MultX_Im = tShiftX[ix + 1];

        for (long k = 0; k < HowMany; k++)
        {
            T buf_r = pData[ix + k * Nx2];
            T buf_im = pData[ix + k * Nx2 + 1];

            T NewRe = buf_r * MultX_Re - buf_im * MultX_Im;
            T NewIm = buf_r * MultX_Im + buf_im * MultX_Re;
            pData[ix + k * Nx2] = NewRe;
            pData[ix + k * Nx2 + 1] = NewIm;
        }
    }
}

void CGenMathFFT1D::RepairSignAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAfter1DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::RotateDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx) 
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RotateDataAfter1DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2, Nx);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::RepairAndRotateDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx, float Mult) 
{

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif


    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + (((Nx / 2) & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairAndRotateAfter1DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, HowMany, Nx, Mult);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::NormalizeDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx, double Mult) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    NormalizeDataAfter1DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2, (float)Mult); //OC06092023
    //NormalizeDataAfter1DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2, Mult);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::FillArrayShift_GPU(double t0, double tStep, long Nx, float* tShiftX) 
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    FillArrayShift_Kernel<float> << <blocks, threads >> > (t0, tStep, Nx, tShiftX);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::TreatShift_GPU(float* pData, long HowMany, long Nx, float* tShiftX) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    TreatShift_Kernel<float> << <blocks, threads >> > (pData, HowMany, Nx * 2, tShiftX);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::RepairSignAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAfter1DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::RotateDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx) 
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx & (2 * GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RotateDataAfter1DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2, Nx);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::RepairAndRotateDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx, double Mult) 
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + (((Nx / 2) & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairAndRotateAfter1DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, HowMany, Nx, (float)Mult); //OC06092023 (check why it's not ..T Mult..)
    //RepairAndRotateAfter1DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, HowMany, Nx, Mult);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::NormalizeDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx, double Mult) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    NormalizeDataAfter1DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, HowMany, Nx * 2, Mult);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::FillArrayShift_GPU(double t0, double tStep, long Nx, double* tShiftX) 
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx & (2 * GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    FillArrayShift_Kernel<double> << <blocks, threads >> > (t0, tStep, Nx, tShiftX);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}

void CGenMathFFT1D::TreatShift_GPU(double* pData, long HowMany, long Nx, double* tShiftX) 
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0));
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    TreatShift_Kernel<double> << <blocks, threads >> > (pData, HowMany, Nx * 2, tShiftX);

//#ifdef _DEBUG
//	cudaStreamSynchronize(0);
//	auto err = cudaGetLastError();
//	printf("%s\r\n", cudaGetErrorString(err));
//#endif
}


template <typename T> __global__ void RepairSignAfter2DFFT_Kernel(T* pAfterFFT, long Nx, long Ny) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //Nx range
    int iy = (blockIdx.y * blockDim.y + threadIdx.y); //Ny range

    float sx0 = 1 - 2 * (ix % 2);
    float sy0 = 1 - 2 * (iy % 2);
    float s = sx0 * sy0;

    if (ix < Nx && iy < Ny) 
    {
        pAfterFFT[(ix + iy * Nx) * 2] *= s;
        pAfterFFT[(ix + iy * Nx) * 2 + 1] *= s;
    }
}

template <typename T> __global__ void RotateDataAfter2DFFT_Kernel(T* pAfterFFT, long HalfNx, long Nx, long HalfNy, long Ny) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //HalfNx range
    int iy = (blockIdx.y * blockDim.y + threadIdx.y); //HalfNy range

    if (ix < HalfNx && iy < HalfNy) 
    {
        long long idx = ((long long)ix + (long long)iy * Nx) * 2;
        //int idx = (ix + iy * Nx) * 2;
        long long HalfNyNx = ((long long)HalfNy) * ((long long)Nx); //HG26022024
        T* t1 = pAfterFFT, *t2 = pAfterFFT + (HalfNyNx + HalfNx) * 2;
        T* t3 = pAfterFFT + HalfNx * 2, *t4 = pAfterFFT + HalfNyNx * 2;

        T buf_r = t1[idx];
        T buf_im = t1[idx + 1];
        t1[idx] = t2[idx];
        t1[idx + 1] = t2[idx + 1];

        t2[idx] = buf_r;
        t2[idx + 1] = buf_im;

        buf_r = t3[idx];
        buf_im = t3[idx + 1];
        t3[idx] = t4[idx];
        t3[idx + 1] = t4[idx + 1];

        t4[idx] = buf_r;
        t4[idx + 1] = buf_im;
    }
}

template <typename T, typename T2> __global__ void RepairSignAndRotateDataAfter2DFFT_Kernel(T* pAfterFFT, long HalfNx, long Nx, long HalfNy, long Ny, T2 Mult) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x); //HalfNx range
    int iy = (blockIdx.y * blockDim.y + threadIdx.y); //HalfNy range

    if (ix < HalfNx) 
    {
        float sx0 = 1.f - 2.f * (ix % 2);
        float sy0 = 1.f - 2.f * (iy % 2);
        float sx1 = 1.f - 2.f * ((HalfNx + ix) % 2);
        float sy1 = 1.f - 2.f * ((HalfNy + iy) % 2);
        
        float s1 = sx0 * sy0 * Mult;
        float s2 = sx1 * sy1 * Mult;
        float s3 = sx1 * sy0 * Mult;
        float s4 = sx0 * sy1 * Mult;

        long long idx = ((long long)ix + (long long)iy * Nx); //HG26022024
        //int idx = (ix + iy * Nx);
        
        long long HalfNyNx = ((long long)HalfNy) * ((long long)Nx);
        T* t1 = pAfterFFT, *t2 = pAfterFFT + (HalfNyNx + HalfNx);
        T* t3 = pAfterFFT + HalfNx, *t4 = pAfterFFT + HalfNyNx;

        T buf1 = t1[idx];
        buf1.x *= s1;
        buf1.y *= s1;

        T buf2 = t2[idx];
        buf2.x *= s2;
        buf2.y *= s2;

        t1[idx] = buf2;
        t2[idx] = buf1;

        buf1 = t3[idx];
        buf1.x *= s3;
        buf1.y *= s3;

        buf2 = t4[idx];
        buf2.x *= s4;
        buf2.y *= s4;

        t3[idx] = buf2;
        t4[idx] = buf1;
    }
}

template <typename T> __global__ void NormalizeDataAfter2DFFT_Kernel(T* pAfterFFT, long Nx, long Ny, T Mult) //HG26022024
//template <typename T> __global__ void NormalizeDataAfter2DFFT_Kernel(T* pAfterFFT, long Nx2Ny2, long n, T Mult) 
{
    long ix = (blockIdx.x * blockDim.x + threadIdx.x); //Nx range //HG26022024
    //int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 2; //Nx range
    long iy = (blockIdx.y * blockDim.y + threadIdx.y); //Ny range //HG26022024

    if(ix < Nx && iy < Ny) //HG26022024
    {
        long long i = ((long long)iy * Nx + (long long)ix) * 2;
        pAfterFFT[i] *= Mult;
        pAfterFFT[i + 1] *= Mult;
    }

    //HG26022024 (commented out)
    //if (ix < Nx2Ny2) 
    //{
    //    pAfterFFT[ix] *= Mult;
    //    pAfterFFT[ix + 1] *= Mult;
    //}
}

template <typename T, bool NeedsShiftX, bool NeedsShiftY> __global__ void TreatShift2D_Kernel(T* pData, long Nx2, long Ny, T* tShiftX, T* tShiftY) 
{
    int ix = (blockIdx.x * blockDim.x + threadIdx.x) * 2; //Nx range
    int iy = (blockIdx.y * blockDim.y + threadIdx.y); //Ny range

    if (ix < Nx2) 
    {
        T MultRe = 1;
        T MultIm = 0;

        T MultX_Re = 1; 
        T MultX_Im = 0;

        T MultY_Re = 1;
        T MultY_Im = 0;

        if (NeedsShiftY)
        {
            MultY_Re = tShiftY[iy * 2];
            MultY_Im = tShiftY[iy * 2 + 1];
        }
        if (NeedsShiftX)
        {
            MultX_Re = tShiftX[ix];
            MultX_Im = tShiftX[ix + 1];

            if (NeedsShiftY) 
            {
                MultRe = MultX_Re * MultY_Re - MultX_Im * MultY_Im;
                MultIm = MultX_Re * MultY_Im + MultX_Im * MultY_Re;
            }
            else 
            {
                MultRe = MultX_Re;
                MultIm = MultX_Im;
            }
        }
        else 
        {
            MultRe = MultY_Re;
            MultIm = MultY_Im;
        }

        T* pData_tmp = pData + (long long)iy * Nx2 + ix; //HG26022024
        T buf_r = pData_tmp[0];
        T buf_im = pData_tmp[1];
        //long offset = iy * Nx2 + ix;
        //T buf_r = pData[offset];
        //T buf_im = pData[offset + 1];
        T NewRe = buf_r * MultRe - buf_im * MultIm;
        T NewIm = buf_r * MultIm + buf_im * MultRe;
        pData_tmp[0] = NewRe; //HG26022024
        pData_tmp[1] = NewIm;
        //pData[offset] = NewRe;
        //pData[offset + 1] = NewIm;
    }
}

void CGenMathFFT2D::RepairSignAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny)
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAfter2DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, Nx, Ny);
}

void CGenMathFFT2D::RotateDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny)
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RotateDataAfter2DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, Nx / 2, Nx, Ny / 2, Ny);
}

void CGenMathFFT2D::RepairSignAndRotateDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, float Mult)
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny/2);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAndRotateDataAfter2DFFT_Kernel<float2, float> << <blocks, threads >> > ((float2*)pAfterFFT, Nx / 2, Nx, Ny / 2, Ny, Mult);
}

void CGenMathFFT2D::NormalizeDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, double Mult)
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny); //HG26022024
    //dim3 blocks((Nx * Ny) / GMFFT_BLOCK_SIZE + (((Nx * Ny) & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    NormalizeDataAfter2DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, Nx, Ny, (float)Mult); //HG26022024
    //NormalizeDataAfter2DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, Nx * Ny * 2, 1, (float)Mult); //OC06092023
    //NormalizeDataAfter2DFFT_Kernel<float> << <blocks, threads >> > (pAfterFFT, Nx * Ny * 2, howMany,1, Mult);
}

void CGenMathFFT2D::TreatShifts2D_GPU(float* pData, long Nx, long Ny, bool NeedsShiftX, bool NeedsShiftY, float* m_ArrayShiftX, float* m_ArrayShiftY)
{

    dim3 blocks((Nx) / GMFFT_BLOCK_SIZE + (((Nx) & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    
    if (NeedsShiftX && NeedsShiftY) TreatShift2D_Kernel<float, true, true> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
    else if (NeedsShiftX) TreatShift2D_Kernel<float, true, false> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
    else if (NeedsShiftY) TreatShift2D_Kernel<float, false, true> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
}

void CGenMathFFT2D::RepairSignAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny)
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAfter2DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, Nx, Ny);
}

void CGenMathFFT2D::RotateDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny)
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RotateDataAfter2DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, Nx / 2, Nx, Ny / 2, Ny);
}

void CGenMathFFT2D::RepairSignAndRotateDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, double Mult)
{

    dim3 blocks(Nx / (2 * GMFFT_BLOCK_SIZE) + ((Nx / 2 & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny/2);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    RepairSignAndRotateDataAfter2DFFT_Kernel<double2, double> << <blocks, threads >> > ((double2*)pAfterFFT, Nx / 2, Nx, Ny / 2, Ny, Mult);
}

void CGenMathFFT2D::NormalizeDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, double Mult)
{

    dim3 blocks(Nx / GMFFT_BLOCK_SIZE + ((Nx & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny); //HG26022024
    //dim3 blocks((Nx * Ny) / GMFFT_BLOCK_SIZE + (((Nx * Ny) & (GMFFT_BLOCK_SIZE - 1)) != 0), 1);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);
    NormalizeDataAfter2DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, Nx, Ny, Mult); //HG26022024
    //NormalizeDataAfter2DFFT_Kernel<double> << <blocks, threads >> > (pAfterFFT, Nx * Ny * 2, 1, Mult);
}

void CGenMathFFT2D::TreatShifts2D_GPU(double* pData, long Nx, long Ny, bool NeedsShiftX, bool NeedsShiftY, double* m_ArrayShiftX, double* m_ArrayShiftY)
{

    dim3 blocks((Nx) / GMFFT_BLOCK_SIZE + (((Nx) & (GMFFT_BLOCK_SIZE - 1)) != 0), Ny);
    dim3 threads(GMFFT_BLOCK_SIZE, 1);

    if (NeedsShiftX && NeedsShiftY) TreatShift2D_Kernel<double, true, true> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
    else if (NeedsShiftX) TreatShift2D_Kernel<double, true, false> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
    else if (NeedsShiftY) TreatShift2D_Kernel<double, false, true> << <blocks, threads >> > (pData, Nx * 2, Ny, m_ArrayShiftX, m_ArrayShiftY);
}
#endif