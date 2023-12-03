/************************************************************************//**
 * File: gmfft_gpu.h
 * Description: Auxiliary utilities to work with FFTW library (CUDA header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2023
 *
 * Copyright (C) Brookhaven National Laboratory
 * All Rights Reserved
 *
 * @author H.Goel
 * @version 1.0
 ***************************************************************************/

#ifndef __GMFFTGPU0_H
#define __GMFFTGPU0_H

void RepairSignAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx);
void RotateDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx);
void RepairAndRotateDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx, float Mult=1.f);
void NormalizeDataAfter1DFFT_GPU(float* pAfterFFT, long HowMany, long Nx, double Mult);
void FillArrayShift_GPU(double t0, double tStep, long Nx, float* tShiftX);
void TreatShift_GPU(float* pData, long HowMany, long Nx, float* tShiftX);

void RepairSignAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx);
void RotateDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx);
void RepairAndRotateDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx, double Mult=1.);
void NormalizeDataAfter1DFFT_GPU(double* pAfterFFT, long HowMany, long Nx, double Mult);
void FillArrayShift_GPU(double t0, double tStep, long Nx, double* tShiftX);
void TreatShift_GPU(double* pData, long HowMany, long Nx, double* tShiftX);

void RepairSignAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, long howMany);
void RotateDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, long howMany);
void RepairSignAndRotateDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, long howMany, float Mult=1.f); //to check
void NormalizeDataAfter2DFFT_GPU(float* pAfterFFT, long Nx, long Ny, long howMany, double Mult);
void TreatShifts2D_GPU(float* pData, long Nx, long Ny, long howMany, bool NeedsShiftX, bool NeedsShiftY, float* m_ArrayShiftX, float* m_ArrayShiftY);

void RepairSignAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, long howMany);
void RotateDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, long howMany);
void RepairSignAndRotateDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, long howMany, double Mult=1.);
void NormalizeDataAfter2DFFT_GPU(double* pAfterFFT, long Nx, long Ny, long howMany, double Mult);
void TreatShifts2D_GPU(double* pData, long Nx, long Ny, long howMany, bool NeedsShiftX, bool NeedsShiftY, double* m_ArrayShiftX, double* m_ArrayShiftY);

#endif // __GMFFTGPU0_H