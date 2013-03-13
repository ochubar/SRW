/************************************************************************//**
 * File: srremflp.cpp
 * Description: Auxiliary (obsolete or rarely used) class for processing Radiation data
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srremflp.h"

//#ifdef __IGOR_PRO__
//#ifndef __SRIGINTR_H
//#include "srigintr.h"
//#endif
//#else
//#ifndef __SRIGORRE_H
//#include "srigorre.h"
//#endif
//#endif

//*************************************************************************

int srTAuxRemoveFlips::GenRemoveFlips(srTWaveAccessData& WaveData)
{
	int result;
	if(WaveData.AmOfDims == 1)
	{
		if(*(WaveData.WaveType) == 'd')
		{
			DOUBLE* pData = (DOUBLE*)(WaveData.pWaveData);
			RemoveFlips1D(pData, *(WaveData.DimSizes), 0, double(*pData));
		}
		else if(*(WaveData.WaveType) == 'f')
		{
			float* pData = (float*)(WaveData.pWaveData);
			RemoveFlips1D(pData, *(WaveData.DimSizes), 0, double(*pData));
		}
		else return NT_FP32_OR_NT_FP64_WAVE_REQUIRED;
	}
	else if(WaveData.AmOfDims == 2)
	{
		if(*(WaveData.WaveType) == 'd') 
		{
			if(result = RemoveFlips2D_D(WaveData)) return result;
		}
		else if(*(WaveData.WaveType) == 'f')
		{
			if(result = RemoveFlips2D_F(WaveData)) return result;
		}
		else return NT_FP32_OR_NT_FP64_WAVE_REQUIRED;
	}
	else return NEEDS_1D_OR_2D_WAVE;

	return 0;
}

//*************************************************************************

void srTAuxRemoveFlips::RemoveFlips1D(DOUBLE* Slice, long Np, long i0, double Phi0)
{
	const double TwoPi = 6.2831853071796;
	const double cFlip = TwoPi - 2.5; // To steer
	double PhToAdd0 = (i0 != -1)? (Phi0 - Slice[i0]) : 0.;

	long HalfNp = Np >> 1;
	long OtherHalfNp = Np - HalfNp;

	double PhToAdd = PhToAdd0;
	DOUBLE *t = Slice + HalfNp - 1; 
	*t += PhToAdd;
	double PrevPh = *(t--);
	for(long i=0; i<(HalfNp - 1); i++)
	{
		*t += PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t--);
	}

	PhToAdd = PhToAdd0;
	t = Slice + HalfNp - 1;
	PrevPh = *(t++);
	for(long j=0; j<OtherHalfNp; j++)
	{
		*t += PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t++);
	}
}

//*************************************************************************

void srTAuxRemoveFlips::RemoveFlips1D(float* Slice, long Np, long i0, double Phi0)
{
	const double TwoPi = 6.2831853071796;
	const double cFlip = TwoPi - 2.5; // To steer
	double PhToAdd0 = (i0 != -1)? (Phi0 - Slice[i0]) : 0.;

	long HalfNp = Np >> 1;
	long OtherHalfNp = Np - HalfNp;

	double PhToAdd = PhToAdd0;
	float *t = Slice + HalfNp - 1; 
	*t += (float)PhToAdd;
	double PrevPh = *(t--);
	for(long i=0; i<(HalfNp - 1); i++)
	{
		*t += (float)PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += (float)TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= (float)TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t--);
	}

	PhToAdd = PhToAdd0;
	t = Slice + HalfNp - 1;
	PrevPh = *(t++);
	for(long j=0; j<OtherHalfNp; j++)
	{
		*t += (float)PhToAdd;
		if(::fabs(*t - PrevPh) > cFlip)
		{
			if(*t < PrevPh)
			{
				*t += (float)TwoPi;
				PhToAdd += TwoPi;
			}
			else
			{
				*t -= (float)TwoPi;
				PhToAdd -= TwoPi;
			}
		}
		PrevPh = *(t++);
	}
}

//*************************************************************************

int srTAuxRemoveFlips::RemoveFlips2D_D(srTWaveAccessData& WaveData)
{
	long Nx = WaveData.DimSizes[0];
	long Nz = WaveData.DimSizes[1];

	DOUBLE* CenterSlice = new DOUBLE[Nx];
	if(CenterSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	//long ixMid = Nx >> 1, izMid = Nz >> 1;
	long izMid = Nz >> 1;
	long ix, iz;

	DOUBLE *pData = (DOUBLE*)(WaveData.pWaveData);
	DOUBLE *tm = pData + izMid*Nx;
	DOUBLE *t = CenterSlice;
	for(ix=0; ix<Nx; ix++) *(t++) = *(tm++);
	RemoveFlips1D(CenterSlice, Nx, -1, 0.);

	DOUBLE* AuxSlice = new DOUBLE[Nz];
	if(AuxSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	for(ix=0; ix<Nx; ix++)
	{
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *(t++) = *tm; tm += Nx;}
		RemoveFlips1D(AuxSlice, Nz, izMid, CenterSlice[ix]);
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *tm = *(t++); tm += Nx;}
	}

	if(AuxSlice != 0) delete[] AuxSlice;
	if(CenterSlice != 0) delete[] CenterSlice;
	return 0;
}

//*************************************************************************

int srTAuxRemoveFlips::RemoveFlips2D_F(srTWaveAccessData& WaveData)
{
	long Nx = WaveData.DimSizes[0];
	long Nz = WaveData.DimSizes[1];

	float* CenterSlice = new float[Nx];
	if(CenterSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	//long ixMid = Nx >> 1, izMid = Nz >> 1;
	long izMid = Nz >> 1;
	long ix, iz;

	float *pData = (float*)(WaveData.pWaveData);
	float *tm = pData + izMid*Nx;
	float *t = CenterSlice;
	for(ix=0; ix<Nx; ix++) *(t++) = *(tm++);
	RemoveFlips1D(CenterSlice, Nx, -1, 0.);

	float* AuxSlice = new float[Nz];
	if(AuxSlice == 0) return MEMORY_ALLOCATION_FAILURE;

	for(ix=0; ix<Nx; ix++)
	{
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *(t++) = *tm; tm += Nx;}
		RemoveFlips1D(AuxSlice, Nz, izMid, double(CenterSlice[ix]));
		t = AuxSlice; tm = pData + ix;
		for(iz=0; iz<Nz; iz++) { *tm = *(t++); tm += Nx;}
	}

	if(AuxSlice != 0) delete[] AuxSlice;
	if(CenterSlice != 0) delete[] CenterSlice;
	return 0;
}

//*************************************************************************
