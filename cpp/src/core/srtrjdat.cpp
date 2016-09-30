/************************************************************************//**
 * File: srtrjdat.cpp
 * Description: Electron Trajectory (and relevant characteristics) calculation for different types of Magnetic Fields
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srtrjdat.h"
#include "srmagfld.h"

//#include <cmath>
#include <math.h>

//*************************************************************************

//int srTFieldBasedArrays::AllocateArrays(long InNs, srTFieldBasedArrayKeys& Keys)
int srTFieldBasedArrays::AllocateArrays(long long InNs, srTFieldBasedArrayKeys& Keys)
{
	if(InNs <= 0) return 0;
	DisposeArrays();

	if(Keys.Bx_) { BxArr = new double[InNs]; if(BxArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Bz_) { BzArr = new double[InNs]; if(BzArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Btx_) { BtxArr = new double[InNs]; if(BtxArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Btz_) { BtzArr = new double[InNs]; if(BtzArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.X_) { XArr = new double[InNs]; if(XArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Z_) { ZArr = new double[InNs]; if(ZArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntBtxE2_) { IntBtxE2Arr = new double[InNs]; if(IntBtxE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntBtzE2_) { IntBtzE2Arr = new double[InNs]; if(IntBtzE2Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.dBxds_) { dBxdsArr = new double[InNs]; if(dBxdsArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.dBzds_) { dBzdsArr = new double[InNs]; if(dBzdsArr == 0) return MEMORY_ALLOCATION_FAILURE;}

	if(Keys.X1p_) { X1pArr = new double[InNs]; if(X1pArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Z1p_) { Z1pArr = new double[InNs]; if(Z1pArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.X2p_) { X2pArr = new double[InNs]; if(X2pArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Z2p_) { Z2pArr = new double[InNs]; if(Z2pArr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.X1_) { X1Arr = new double[InNs]; if(X1Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Z1_) { Z1Arr = new double[InNs]; if(Z1Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.X2_) { X2Arr = new double[InNs]; if(X2Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.Z2_) { Z2Arr = new double[InNs]; if(Z2Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntX01_) { IntX01Arr = new double[InNs]; if(IntX01Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntX02_) { IntX02Arr = new double[InNs]; if(IntX02Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntX11_) { IntX11Arr = new double[InNs]; if(IntX11Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntX12_) { IntX12Arr = new double[InNs]; if(IntX12Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntX22_) { IntX22Arr = new double[InNs]; if(IntX22Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntZ01_) { IntZ01Arr = new double[InNs]; if(IntZ01Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntZ02_) { IntZ02Arr = new double[InNs]; if(IntZ02Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntZ11_) { IntZ11Arr = new double[InNs]; if(IntZ11Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntZ12_) { IntZ12Arr = new double[InNs]; if(IntZ12Arr == 0) return MEMORY_ALLOCATION_FAILURE;}
	if(Keys.IntZ22_) { IntZ22Arr = new double[InNs]; if(IntZ22Arr == 0) return MEMORY_ALLOCATION_FAILURE;}

	Ns = InNs;
	return 0;
}

//*************************************************************************

srTGenTrjDat* srTGenTrjDat::CreateAndSetupNewTrjDat(srTEbmDat* pEbmDat, srTMagElem* pMagElem)
{
	if(pMagElem == 0) return 0;

	srTGenTrjDat* pOut = pMagElem->CreateAndSetupNewTrjDat(pEbmDat);
	if(pOut == 0) return 0;
	if(pEbmDat != 0) pOut->EbmDat = *pEbmDat;
	return pOut;
}

//*************************************************************************

int srTTrjDat::AllocateMemoryForCfs()
{
	DeallocateMemoryForCfs();

	//int LenFieldData_m_1 = LenFieldData - 1;
	long long LenFieldData_m_1 = LenFieldData - 1;

	BxPlnCf = new double*[LenFieldData_m_1];
	if(BxPlnCf == 0) return MEMORY_ALLOCATION_FAILURE;

	BzPlnCf = new double*[LenFieldData_m_1];
	if(BzPlnCf == 0) { delete[] BxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtxPlnCf = new double*[LenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtzPlnCf = new double*[LenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	xPlnCf = new double*[LenFieldData_m_1];
	if(xPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	zPlnCf = new double*[LenFieldData_m_1];
	if(zPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtx2PlnCf = new double*[LenFieldData_m_1];
	if(IntBtx2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtz2PlnCf = new double*[LenFieldData_m_1];
	if(IntBtz2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	AllCf = new double[LenFieldData_m_1*50];
	if(AllCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; delete[] IntBtz2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	double* tAllCf = AllCf;
	//for(int i=0; i<LenFieldData_m_1; i++)
	for(long long i=0; i<LenFieldData_m_1; i++)
	{
		BxPlnCf[i] = tAllCf; tAllCf += 4;
		BzPlnCf[i] = tAllCf; tAllCf += 4;

		BtxPlnCf[i] = tAllCf; tAllCf += 5;
		BtzPlnCf[i] = tAllCf; tAllCf += 5;

		xPlnCf[i] = tAllCf; tAllCf += 6;
		zPlnCf[i] = tAllCf; tAllCf += 6;

		IntBtx2PlnCf[i] = tAllCf; tAllCf += 10;
		IntBtz2PlnCf[i] = tAllCf; tAllCf += 10;
	}
	return 0;
}

//*************************************************************************

int srTTrjDat::AllocateMemoryForCfs_FromTrj()
{
	if(HorFieldIsNotZero && (xTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;
	if(VerFieldIsNotZero && (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;

	DeallocateMemoryForCfs();

	//int xLenFieldData_m_1 = xTrjInData.np - 1;
	//int zLenFieldData_m_1 = zTrjInData.np - 1;
	long long xLenFieldData_m_1 = xTrjInData.np - 1;
	long long zLenFieldData_m_1 = zTrjInData.np - 1;

	BxPlnCf = new double*[zLenFieldData_m_1];
	if(BxPlnCf == 0) return MEMORY_ALLOCATION_FAILURE;

	BzPlnCf = new double*[xLenFieldData_m_1];
	if(BzPlnCf == 0) { delete[] BxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtxPlnCf = new double*[xLenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtzPlnCf = new double*[zLenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	xPlnCf = new double*[xLenFieldData_m_1];
	if(xPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	zPlnCf = new double*[zLenFieldData_m_1];
	if(zPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtx2PlnCf = new double*[xLenFieldData_m_1];
	if(IntBtx2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtz2PlnCf = new double*[zLenFieldData_m_1];
	if(IntBtz2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	AllCf = new double[(xLenFieldData_m_1 + zLenFieldData_m_1)*21];
	if(AllCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; delete[] IntBtz2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	double* tAllCf = AllCf;
	//for(long ix=0; ix<xLenFieldData_m_1; ix++)
	for(long long ix=0; ix<xLenFieldData_m_1; ix++)
	{
		BzPlnCf[ix] = tAllCf; tAllCf += 4;
		BtxPlnCf[ix] = tAllCf; tAllCf += 5;
		xPlnCf[ix] = tAllCf; tAllCf += 6;
		IntBtx2PlnCf[ix] = tAllCf; tAllCf += 6;
	}
	//for(long iz=0; iz<zLenFieldData_m_1; iz++)
	for(long long iz=0; iz<zLenFieldData_m_1; iz++)
	{
		BxPlnCf[iz] = tAllCf; tAllCf += 4;
		BtzPlnCf[iz] = tAllCf; tAllCf += 5;
		zPlnCf[iz] = tAllCf; tAllCf += 6;
		IntBtz2PlnCf[iz] = tAllCf; tAllCf += 6;
	}

	return 0;
}

//*************************************************************************

//int srTTrjDat::AllocateMemoryForCfsFromTrj(int npLoc)
int srTTrjDat::AllocateMemoryForCfsFromTrj(long long npLoc)
{
	DeallocateMemoryForCfs();

	//int xLenFieldData_m_1 = npLoc - 1;
	//int zLenFieldData_m_1 = npLoc - 1;
	long long xLenFieldData_m_1 = npLoc - 1;
	long long zLenFieldData_m_1 = npLoc - 1;

	BxPlnCf = new double*[zLenFieldData_m_1];
	if(BxPlnCf == 0) return MEMORY_ALLOCATION_FAILURE;

	BzPlnCf = new double*[xLenFieldData_m_1];
	if(BzPlnCf == 0) { delete[] BxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtxPlnCf = new double*[xLenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	BtzPlnCf = new double*[zLenFieldData_m_1];
	if(BtxPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	xPlnCf = new double*[xLenFieldData_m_1];
	if(xPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	zPlnCf = new double*[zLenFieldData_m_1];
	if(zPlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtx2PlnCf = new double*[xLenFieldData_m_1];
	if(IntBtx2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; return MEMORY_ALLOCATION_FAILURE;}

	IntBtz2PlnCf = new double*[zLenFieldData_m_1];
	if(IntBtz2PlnCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	AllCf = new double[(xLenFieldData_m_1 + zLenFieldData_m_1)*21];
	if(AllCf == 0) { delete[] BxPlnCf; delete[] BzPlnCf; delete[] BtxPlnCf; delete[] BtzPlnCf; delete[] xPlnCf; delete[] zPlnCf; delete[] IntBtx2PlnCf; delete[] IntBtz2PlnCf; return MEMORY_ALLOCATION_FAILURE;}

	double* tAllCf = AllCf;
	//for(long ix=0; ix<xLenFieldData_m_1; ix++)
	for(long long ix=0; ix<xLenFieldData_m_1; ix++)
	{
		BzPlnCf[ix] = tAllCf; tAllCf += 4;
		BtxPlnCf[ix] = tAllCf; tAllCf += 5;
		xPlnCf[ix] = tAllCf; tAllCf += 6;
		IntBtx2PlnCf[ix] = tAllCf; tAllCf += 6;
	}
	//for(long iz=0; iz<zLenFieldData_m_1; iz++)
	for(long long iz=0; iz<zLenFieldData_m_1; iz++)
	{
		BxPlnCf[iz] = tAllCf; tAllCf += 4;
		BtzPlnCf[iz] = tAllCf; tAllCf += 5;
		zPlnCf[iz] = tAllCf; tAllCf += 6;
		IntBtz2PlnCf[iz] = tAllCf; tAllCf += 6;
	}

	return 0;
}

//*************************************************************************

int srTTrjDat::DeallocateMemoryForCfs()
{
	//int LenFieldData_m_1 = LenFieldData-1;
	if(BxPlnCf != 0)
	{
		delete[] BxPlnCf; BxPlnCf = 0;
	}
	if(BzPlnCf != 0)
	{
		delete[] BzPlnCf; BzPlnCf = 0;
	}
	if(BtxPlnCf != 0)
	{
		delete[] BtxPlnCf; BtxPlnCf = 0;
	}
	if(BtzPlnCf != 0)
	{
		delete[] BtzPlnCf; BtzPlnCf = 0;
	}
	if(xPlnCf != 0)
	{
		delete[] xPlnCf; xPlnCf = 0;
	}
	if(zPlnCf != 0)
	{
		delete[] zPlnCf; zPlnCf = 0;
	}
	if(IntBtx2PlnCf != 0)
	{
		delete[] IntBtx2PlnCf; IntBtx2PlnCf = 0;
	}
	if(IntBtz2PlnCf != 0)
	{
		delete[] IntBtz2PlnCf; IntBtz2PlnCf = 0;
	}

	if(AllCf != 0)
	{
		delete[] AllCf; AllCf = 0;
	}
	return 0;
}

//*************************************************************************

int srTTrjDat::CompDerivForFieldData(srTFunDer* InitialFieldData)
{
	double SubArray[5];

	for(int k=0; k<5; k++) 
	{
		double& fR = InitialFieldData[k].f;
		SubArray[k] = fR;
	}

	InitialFieldData[0].dfds = Derivative(SubArray, sStep, 0);
	InitialFieldData[1].dfds = Derivative(SubArray, sStep, 1);
	InitialFieldData[2].dfds = Derivative(SubArray, sStep, 2);

	//int LenFieldData_m_2 = LenFieldData - 2;
	long long LenFieldData_m_2 = LenFieldData - 2;
	//for(int i=3; i<LenFieldData_m_2; i++)
	for(long long i=3; i<LenFieldData_m_2; i++)
	{
		//int i_m_2 = i - 2;
		long long i_m_2 = i - 2;
		for(int k=0; k<5; k++) 
		{
			double& fR = InitialFieldData[i_m_2 + k].f;
			SubArray[k] = fR;
		}
		InitialFieldData[i].dfds = Derivative(SubArray, sStep, 2);
	}

	InitialFieldData[LenFieldData_m_2].dfds = Derivative(SubArray, sStep, 3);
	InitialFieldData[LenFieldData - 1].dfds = Derivative(SubArray, sStep, 4);

	return 1;
}

//*************************************************************************

void srTTrjDat::SetupIntegrPlnCfs(char XorZ)
{
	//int i, j;
	int j;
	double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;
	double BtPr = 0., C_Pr = 0., IntBt2Pr = 0.;

	int angZeroCount = 0; //OC220112
	double BtPrPr = 0.;

	//for(i=0; i<LenFieldData-1; i++)
	for(long long i=0; i<LenFieldData-1; i++)
	{
		if(XorZ=='x') // Vertical field is responsible for horizontal dynamics
		{
			B_CfP = BzPlnCf[i]; Bt_CfP = BtxPlnCf[i]; C_CfP = xPlnCf[i]; IntBt2_CfP = IntBtx2PlnCf[i];
		}
		else
		{
			B_CfP = BxPlnCf[i]; Bt_CfP = BtzPlnCf[i]; C_CfP = zPlnCf[i]; IntBt2_CfP = IntBtz2PlnCf[i];
		}

		for(j=0; j<4; j++) Bt_CfP[j+1] = B_CfP[j]/(j+1);
		Bt_CfP[0] = BtPr;

		for(j=0; j<5; j++) C_CfP[j+1] = Bt_CfP[j]/(j+1);
		C_CfP[0] = C_Pr;

		double Bt0 = Bt_CfP[0], Bt1 = Bt_CfP[1], Bt2 = Bt_CfP[2], Bt3 = Bt_CfP[3], Bt4 = Bt_CfP[4];
		IntBt2_CfP[1] = Bt0*Bt0;
		IntBt2_CfP[2] = Bt0*Bt1;
		IntBt2_CfP[3] = (2.*Bt0*Bt2 + Bt1*Bt1)/3.;
		IntBt2_CfP[4] = (Bt0*Bt3 + Bt1*Bt2)/2.;
		IntBt2_CfP[5] = (2.*Bt0*Bt4 + 2.*Bt1*Bt3 + Bt2*Bt2)/5.;
		IntBt2_CfP[6] = (Bt1*Bt4 + Bt2*Bt3)/3.;
		IntBt2_CfP[7] = (2.*Bt2*Bt4 + Bt3*Bt3)/7.;
		IntBt2_CfP[8] = Bt3*Bt4/4.;
		IntBt2_CfP[9] = Bt4*Bt4/9.;
		IntBt2_CfP[0] = IntBt2Pr;

		BtPr = Pol04(sStep, &(Bt_CfP[1])) + Bt0;
		C_Pr = Pol05(sStep, &(C_CfP[1])) + C_CfP[0];
		IntBt2Pr = Pol09(sStep, &(IntBt2_CfP[1])) + IntBt2_CfP[0];

		if(BtPrPr*BtPr < 0) angZeroCount++;
		BtPrPr = BtPr;
	}

	if(angZeroCount > 0) //OC220112
	{
		if(m_estimMinNpForRadInteg < angZeroCount) m_estimMinNpForRadInteg = angZeroCount;
	}
}

//*************************************************************************
/**
void srTTrjDat::CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{// Attention: length units are m at the input/output here !
	if(CompFromTrj) { CompTotalTrjData_FromTrj(Keys, FieldBasedArrays); return;}
	
	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	double dxds0E2 = dxds0*dxds0;
	double dzds0E2 = dzds0*dzds0;
	double s = FieldBasedArrays.sStart;
	double sStp = FieldBasedArrays.sStep;

	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	double *pX = FieldBasedArrays.XArr, *pZ = FieldBasedArrays.ZArr;
	double *pIntBtxE2 = FieldBasedArrays.IntBtxE2Arr, *pIntBtzE2 = FieldBasedArrays.IntBtzE2Arr;
	double *pdBxds = FieldBasedArrays.dBxdsArr, *pdBzds = FieldBasedArrays.dBzdsArr;

	double BufCrd;
	for(int i=0; i<FieldBasedArrays.Ns; i++)
	{
		int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;

		double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;

		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];

			if(Keys.dBzds_) *(pdBzds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
			if(Keys.Bz_) *(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			if(Keys.Btx_) *(pBtx++) = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			if(Keys.X_ || Keys.IntBtxE2_) BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			if(Keys.X_) *(pX++) = xCorr + BtxCorrForX*s + BufCrd;
			if(Keys.IntBtxE2_) *(pIntBtxE2++) = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{
			if(Keys.dBzds_) *(pdBzds++) = 0.;
			if(Keys.Bz_) *(pBz++) = 0.;
			double s_mi_s0 = s - s0;
			if(Keys.Btx_) *(pBtx++) = dxds0;
			if(Keys.X_) *(pX++) = x0 + dxds0*s_mi_s0;
			if(Keys.IntBtxE2_) *(pIntBtxE2++) = dxds0E2*s_mi_s0;
		}
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];

			if(Keys.dBxds_) *(pdBxds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
			if(Keys.Bx_) *(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			if(Keys.Btz_) *(pBtz++) = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			if(Keys.Z_ || Keys.IntBtzE2_) BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			if(Keys.Z_) *(pZ++) = zCorr + BtzCorrForZ*s + BufCrd;
			if(Keys.IntBtzE2_) *(pIntBtzE2++) = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{
			if(Keys.dBzds_) *(pdBxds++) = 0.;
			if(Keys.Bx_) *(pBx++) = 0.;
			double s_mi_s0 = s - s0;
			if(Keys.Btz_) *(pBtz++) = dzds0;
			if(Keys.Z_) *(pZ++) = z0 + dzds0*s_mi_s0;
			if(Keys.IntBtzE2_) *(pIntBtzE2++) = dzds0E2*s_mi_s0;
		}
		s += sStp;
	}
}
**/

void srTTrjDat::CompTotalTrjData(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{// Attention: length units are m at the input/output here !
	if(CompFromTrj) { CompTotalTrjData_FromTrj(Keys, FieldBasedArrays); return;}
	
	char AllEqData_ = Keys.AllEqData_ShouldBeSet();
	char AllNonEqTrajData = Keys.AllNonEqTrajData_ShouldBeSet();
	char AllNonEqIntegData = Keys.AllNonEqIntegData_ShouldBeSet();
	char CalcIsForThickBeam = Keys.CheckIfCalcIsForThickBeam();

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	double dxds0E2 = dxds0*dxds0;
	double dzds0E2 = dzds0*dzds0;
	double s = FieldBasedArrays.sStart;
	double sStp = FieldBasedArrays.sStep;

	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	double *pX = FieldBasedArrays.XArr, *pZ = FieldBasedArrays.ZArr;
	double *pIntBtxE2 = FieldBasedArrays.IntBtxE2Arr, *pIntBtzE2 = FieldBasedArrays.IntBtzE2Arr;
	double *pdBxds = FieldBasedArrays.dBxdsArr, *pdBzds = FieldBasedArrays.dBzdsArr;

    double *pX1p = FieldBasedArrays.X1pArr, *pZ1p = FieldBasedArrays.Z1pArr;
    double *pX2p = FieldBasedArrays.X2pArr, *pZ2p = FieldBasedArrays.Z2pArr;
    double *pX1 = FieldBasedArrays.X1Arr, *pZ1 = FieldBasedArrays.Z1Arr;
	double *pX2 = FieldBasedArrays.X2Arr, *pZ2 = FieldBasedArrays.Z2Arr;
	double *pIntX01 = FieldBasedArrays.IntX01Arr, *pIntX02 = FieldBasedArrays.IntX02Arr;
	double *pIntX11 = FieldBasedArrays.IntX11Arr, *pIntX12 = FieldBasedArrays.IntX12Arr, *pIntX22 = FieldBasedArrays.IntX22Arr;
	double *pIntZ01 = FieldBasedArrays.IntZ01Arr, *pIntZ02 = FieldBasedArrays.IntZ02Arr;
	double *pIntZ11 = FieldBasedArrays.IntZ11Arr, *pIntZ12 = FieldBasedArrays.IntZ12Arr, *pIntZ22 = FieldBasedArrays.IntZ22Arr;

	double BtxCorrForThickBeamCalc = 0., BtzCorrForThickBeamCalc = 0.;
    if(CalcIsForThickBeam) { BtxCorrForThickBeamCalc = -dxds0; BtzCorrForThickBeamCalc = -dzds0;}
	double CurX, CurIntBtxE2, CurZ, CurIntBtzE2, BtxPrev, BtxCur, BtzPrev, BtzCur, BxCur, BzCur, absBxCur, absBzCur;

	double &refAbsBxMax = FieldBasedArrays.absBxMax;
	double &refAbsBzMax = FieldBasedArrays.absBzMax;
	refAbsBxMax = 0; refAbsBzMax = 0;

	double HalfStp = 0.5*sStp;

	double BufCrd;
	//for(int i=0; i<FieldBasedArrays.Ns; i++)
	for(long long i=0; i<FieldBasedArrays.Ns; i++)
	{
		//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;
		double s_mi_s0 = s - s0;

		double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;

        double xCorrForThickBeamCalc = 0., zCorrForThickBeamCalc = 0.;
        if(CalcIsForThickBeam) 
		{ 
			xCorrForThickBeamCalc = -x0 - dxds0*s_mi_s0;
            zCorrForThickBeamCalc = -z0 - dzds0*s_mi_s0;
		}

		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];

			if(AllEqData_)
			{
                *(pdBzds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);

                //*(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
				BzCur = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
				absBzCur = ::fabs(BzCur);
				if(refAbsBzMax < absBzCur) refAbsBzMax = absBzCur;
				*(pBz++) = BzCur;

				BtxCur = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP) + BtxCorrForThickBeamCalc;
                *(pBtx++) = BtxCur;
                BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
				CurX = xCorr + BtxCorrForX*s + BufCrd + xCorrForThickBeamCalc;
                *(pX++) = CurX;
				CurIntBtxE2 = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
				if(CalcIsForThickBeam) CurIntBtxE2 += (-dxds0*(2*CurX + dxds0*s_mi_s0));
                *(pIntBtxE2++) = CurIntBtxE2;
			}
			else
			{
				if(Keys.dBzds_) *(pdBzds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
				if(Keys.Bz_) 
				{
					//*(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
					BzCur = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
					absBzCur = ::fabs(BzCur);
					if(refAbsBzMax < absBzCur) refAbsBzMax = absBzCur;
					*(pBz++) = BzCur;
				}
                if(Keys.Btx_ || Keys.IntX01_ || Keys.IntX02_) BtxCur = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP) + BtxCorrForThickBeamCalc;
				if(Keys.Btx_) *(pBtx++) = BtxCur;
				if(Keys.X_ || Keys.IntBtxE2_ || Keys.X2p_) 
				{
                    BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
                    CurX = xCorr + BtxCorrForX*s + BufCrd + xCorrForThickBeamCalc;
				}
				if(Keys.X_) *(pX++) = CurX;
				if(Keys.IntBtxE2_) 
				{
					CurIntBtxE2 = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
                    if(CalcIsForThickBeam) CurIntBtxE2 += (-dxds0*(2.*CurX + dxds0*s_mi_s0));
					*(pIntBtxE2++) = CurIntBtxE2;
				}
			}
		}
		else 
		{
			if(AllEqData_)
			{
                *(pdBzds++) = 0.;
                *(pBz++) = 0.;
                //double s_mi_s0 = s - s0;
                BtxCur = dxds0 + BtxCorrForThickBeamCalc;
                *(pBtx++) = BtxCur;
				CurX = x0 + dxds0*s_mi_s0 + xCorrForThickBeamCalc;
                *(pX++) = CurX;
				*(pIntBtxE2++) = (CalcIsForThickBeam? 0 : dxds0E2*s_mi_s0);
			}
			else
			{
				if(Keys.dBzds_) *(pdBzds++) = 0.;
				if(Keys.Bz_) *(pBz++) = 0.;
				//double s_mi_s0 = s - s0;
				if(Keys.Btx_ || Keys.IntX01_ || Keys.IntX02_) BtxCur = dxds0 + BtxCorrForThickBeamCalc;
				if(Keys.Btx_) *(pBtx++) = BtxCur;
				if(Keys.X_ || Keys.X2p_) CurX = x0 + dxds0*s_mi_s0 + xCorrForThickBeamCalc;
				if(Keys.X_) *(pX++) = CurX;
				if(Keys.IntBtxE2_) *(pIntBtxE2++) = (CalcIsForThickBeam? 0 : dxds0E2*s_mi_s0);
			}
		}
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];

			if(AllEqData_)
			{
                *(pdBxds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);

                //*(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
				BxCur = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
				absBxCur = ::fabs(BxCur);
				if(refAbsBxMax < absBxCur) refAbsBxMax = absBxCur;
				*(pBx++) = BxCur;

                BtzCur = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP) + BtzCorrForThickBeamCalc;
                *(pBtz++) = BtzCur;
                BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
				CurZ = zCorr + BtzCorrForZ*s + BufCrd + zCorrForThickBeamCalc;
                *(pZ++) = CurZ;
				CurIntBtzE2 = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
                if(CalcIsForThickBeam) CurIntBtzE2 += (-dzds0*(2*CurZ + dzds0*s_mi_s0));
                *(pIntBtzE2++) = CurIntBtzE2;
			}
			else
			{
				if(Keys.dBxds_) *(pdBxds++) = B_CfP[1] + smsb*(2.*B_CfP[2] + 3.*smsb*B_CfP[3]);
				if(Keys.Bx_) 
				{
					//*(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
					BxCur = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
					absBxCur = ::fabs(BxCur);
					if(refAbsBxMax < absBxCur) refAbsBxMax = absBxCur;
					*(pBx++) = BxCur;
				}
				if(Keys.Btz_ || Keys.IntZ01_ || Keys.IntZ02_) BtzCur = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP) + BtzCorrForThickBeamCalc;
				if(Keys.Btz_) *(pBtz++) = BtzCur;
				if(Keys.Z_ || Keys.IntBtzE2_ || Keys.X2p_) 
				{
					BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
					CurZ = zCorr + BtzCorrForZ*s + BufCrd + zCorrForThickBeamCalc;
				}
				if(Keys.Z_) *(pZ++) = CurZ;
				if(Keys.IntBtzE2_) 
				{
					CurIntBtzE2 = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
					if(CalcIsForThickBeam) CurIntBtzE2 += (-dzds0*(2*CurZ + dzds0*s_mi_s0));
					*(pIntBtzE2++) = CurIntBtzE2;
				}
			}
		}
		else 
		{
			if(AllEqData_)
			{
                *(pdBxds++) = 0.;
                *(pBx++) = 0.;
                //double s_mi_s0 = s - s0;
                BtzCur = dzds0 + BtzCorrForThickBeamCalc;
                *(pBtz++) = BtzCur;
				CurZ = z0 + dzds0*s_mi_s0 + zCorrForThickBeamCalc;
                *(pZ++) = CurZ;
                *(pIntBtzE2++) = (CalcIsForThickBeam? 0 : dzds0E2*s_mi_s0);
			}
			else
			{
				if(Keys.dBzds_) *(pdBxds++) = 0.;
				if(Keys.Bx_) *(pBx++) = 0.;
				//double s_mi_s0 = s - s0;
				if(Keys.Btz_ || Keys.IntZ01_ || Keys.IntZ02_) BtzCur = dzds0 + BtzCorrForThickBeamCalc;
				if(Keys.Btz_) *(pBtz++) = BtzCur;
				if(Keys.Z_ || Keys.Z2p_) CurZ = z0 + dzds0*s_mi_s0 + zCorrForThickBeamCalc;
				if(Keys.Z_) *(pZ++) = CurZ;
				if(Keys.IntBtzE2_) *(pIntBtzE2++) = (CalcIsForThickBeam? 0 : dzds0E2*s_mi_s0);
			}
		}

		if(AllNonEqTrajData)
		{
			*(pX1++) = 1;
			*(pX2++) = s_mi_s0;
            *(pX1p++) = 0; 
			*(pX2p++) = 1;
			*(pZ1++) = 1;
			*(pZ2++) = s_mi_s0;
            *(pZ1p++) = 0; 
			*(pZ2p++) = 1;
		}
		else
		{
			if(Keys.X1_) *(pX1++) = 1;
			else pX1++;
			if(Keys.X2_) *(pX2++) = s_mi_s0;
			else pX2++;
            if(Keys.X1p_) *(pX1p++) = 0; 
			else pX1p++;
			if(Keys.X2p_) *(pX2p++) = 1;
			else pX2p++;
			if(Keys.Z1_) *(pZ1++) = 1;
			else pZ1++;
			if(Keys.Z2_) *(pZ2++) = s_mi_s0;
			else pZ2++;
            if(Keys.Z1p_) *(pZ1p++) = 0; 
			else pZ1p++;
			if(Keys.Z2p_) *(pZ2p++) = 1;
			else pZ2p++;
		}

		if(AllNonEqIntegData)
		{
			double X1pPrev = *(pX1p-2), X1pCur = *(pX1p-1);
			double X2pPrev = *(pX2p-2), X2pCur = *(pX2p-1);
			double Z1pPrev = *(pZ1p-2), Z1pCur = *(pZ1p-1);
			double Z2pPrev = *(pZ2p-2), Z2pCur = *(pZ2p-1);

			if(Keys.X1p_ || (i == 0)) 
			{
				*(pIntX01++) = 0;
				*(pIntX11++) = 0;
				*(pIntX12++) = 0;
			}
			else
			{
                *(pIntX01++) = *pIntX01 + HalfStp*(BtxPrev*X1pPrev + BtxCur*X1pCur);
                *(pIntX11++) = *pIntX11 + HalfStp*(X1pPrev*X1pPrev + X1pCur*X1pCur);
                *(pIntX12++) = *pIntX12 + HalfStp*(X1pPrev*X2pPrev + X1pCur*X2pCur);
			}
			if(Keys.X2p_) 
			{
				*(pIntX22++) = s_mi_s0;
				*(pIntX02++) = CurX;
			}
			else
			{
				if(i == 0)
				{
                    *(pIntX22++) = 0;
                    *(pIntX02++) = 0;
				}
				else
				{
                    *(pIntX22++) = *pIntX22 + HalfStp*(X2pPrev*X2pPrev + X2pCur*X2pCur);
                    *(pIntX02++) = *pIntX02 + HalfStp*(BtxPrev*X2pPrev + BtxCur*X2pCur);
				}
			}

			if(Keys.Z1p_ || (i == 0)) 
			{
				*(pIntZ01++) = 0;
				*(pIntZ11++) = 0;
				*(pIntZ12++) = 0;
			}
			else
			{
                *(pIntZ01++) = *pIntZ01 + HalfStp*(BtzPrev*Z1pPrev + BtzCur*Z1pCur);
                *(pIntZ11++) = *pIntZ11 + HalfStp*(Z1pPrev*Z1pPrev + Z1pCur*Z1pCur);
                *(pIntZ12++) = *pIntZ12 + HalfStp*(Z1pPrev*Z2pPrev + Z1pCur*Z2pCur);
			}
			if(Keys.Z2p_) 
			{
				*(pIntZ22++) = s_mi_s0;
				*(pIntZ02++) = CurZ;
			}
			else
			{
				if(i == 0)
				{
                    *(pIntZ22++) = 0;
                    *(pIntZ02++) = 0;
				}
				else
				{
                    *(pIntZ22++) = *pIntZ22 + HalfStp*(Z2pPrev*Z2pPrev + Z2pCur*Z2pCur);
                    *(pIntZ02++) = *pIntZ02 + HalfStp*(BtzPrev*Z2pPrev + BtzCur*Z2pCur);
				}
			}
		}
		else
		{
			if(Keys.IntX01_) 
			{
                if(Keys.X1p_) *(pIntX01++) = 0;
				else 
				{
                    double X1pPrev = *(pX1p-2), X1pCur = *(pX1p-1);
					*(pIntX01++) = ((i == 0)? 0 : (*pIntX01 + HalfStp*(BtxPrev*X1pPrev + BtxCur*X1pCur)));
				}
			}
			if(Keys.IntX11_) 
			{
                if(Keys.X1p_) *(pIntX11++) = 0;
				else 
				{
                    double X1pPrev = *(pX1p-2), X1pCur = *(pX1p-1);
                    *(pIntX11++) = ((i == 0)? 0 : (*pIntX11 + HalfStp*(X1pPrev*X1pPrev + X1pCur*X1pCur)));
				}
			}
			if(Keys.IntX12_) 
			{
                if(Keys.X1p_) *(pIntX12++) = 0;
				else 
				{
                    double X1pPrev = *(pX1p-2), X1pCur = *(pX1p-1);
                    double X2pPrev = *(pX2p-2), X2pCur = *(pX2p-1);
                    *(pIntX12++) = ((i == 0)? 0 : (*pIntX12 + HalfStp*(X1pPrev*X2pPrev + X1pCur*X2pCur)));
				}
			}
			if(Keys.IntX22_) 
			{
				if(Keys.X2p_) *(pIntX22++) = s_mi_s0;
				else 
				{
                    double X2pPrev = *(pX2p-2), X2pCur = *(pX2p-1);
					*(pIntX22++) = ((i == 0)? 0 : (*pIntX22 + HalfStp*(X2pPrev*X2pPrev + X2pCur*X2pCur)));
				}
			}
			if(Keys.IntX02_) 
			{
				if(Keys.X2p_) *(pIntX02++) = CurX;
				else 
				{
                    double X2pPrev = *(pX2p-2), X2pCur = *(pX2p-1);
					*(pIntX02++) = ((i == 0)? 0 : (*pIntX02 + HalfStp*(BtxPrev*X2pPrev + BtxCur*X2pCur)));
				}
			}

			if(Keys.IntZ01_) 
			{
                if(Keys.Z1p_) *(pIntZ01++) = 0;
				else 
				{
                    double Z1pPrev = *(pZ1p-2), Z1pCur = *(pZ1p-1);
					*(pIntZ01++) = ((i == 0)? 0 : (*pIntZ01 + HalfStp*(BtzPrev*Z1pPrev + BtzCur*Z1pCur)));
				}
			}
			if(Keys.IntZ11_) 
			{
                if(Keys.Z1p_) *(pIntZ11++) = 0;
				else 
				{
                    double Z1pPrev = *(pZ1p-2), Z1pCur = *(pZ1p-1);
                    *(pIntZ11++) = ((i == 0)? 0 : (*pIntZ11 + HalfStp*(Z1pPrev*Z1pPrev + Z1pCur*Z1pCur)));
				}
			}
			if(Keys.IntZ12_) 
			{
                if(Keys.Z1p_) *(pIntZ12++) = 0;
				else 
				{
                    double Z1pPrev = *(pZ1p-2), Z1pCur = *(pZ1p-1);
                    double Z2pPrev = *(pZ2p-2), Z2pCur = *(pZ2p-1);
                    *(pIntZ12++) = ((i == 0)? 0 : (*pIntZ12 + HalfStp*(Z1pPrev*Z2pPrev + Z1pCur*Z2pCur)));
				}
			}
			if(Keys.IntZ22_) 
			{
                if(Keys.Z2p_) *(pIntZ22++) = s_mi_s0;
				else 
				{
                    double Z2pPrev = *(pZ2p-2), Z2pCur = *(pZ2p-1);
					*(pIntZ22++) = ((i == 0)? 0 : (*pIntZ22 + HalfStp*(Z2pPrev*Z2pPrev + Z2pCur*Z2pCur)));
				}
			}
			if(Keys.IntZ02_) 
			{
                if(Keys.Z2p_) *(pIntZ02++) = CurZ;
				else 
				{
                    double Z2pPrev = *(pZ2p-2), Z2pCur = *(pZ2p-1);
					*(pIntZ02++) = ((i == 0)? 0 : (*pIntZ02 + HalfStp*(BtzPrev*Z2pPrev + BtzCur*Z2pCur)));
				}
			}
		}
		BtxPrev = BtxCur; BtzPrev = BtzCur;
		s += sStp;
	}

	if(Keys.CheckIfCalcCorrIntIsNecessary())
	{
		//int IndS0;
		long long IndS0;
		double dS0;
		double dIndS0 = (s0 - FieldBasedArrays.sStart)/sStp;
		if(dIndS0 <= 0) { IndS0 = 0; dS0 = 0;}
		else if(dIndS0 >= (FieldBasedArrays.Ns - 1)) { IndS0 = (FieldBasedArrays.Ns - 1); dS0 = 0;}
		else { IndS0 = int(dIndS0); dS0 = s0 - (FieldBasedArrays.sStart + IndS0*sStp);}
		
		if(Keys.IntX01toS0_ && Keys.IntX01_) FieldBasedArrays.IntX01toS0 = EstimateIntermedVal(FieldBasedArrays.IntX01Arr, IndS0, sStp, dS0);
		if(Keys.IntX02toS0_ && Keys.IntX02_) FieldBasedArrays.IntX02toS0 = EstimateIntermedVal(FieldBasedArrays.IntX02Arr, IndS0, sStp, dS0);
		if(Keys.IntX11toS0_ && Keys.IntX11_) FieldBasedArrays.IntX11toS0 = EstimateIntermedVal(FieldBasedArrays.IntX11Arr, IndS0, sStp, dS0);
		if(Keys.IntX12toS0_ && Keys.IntX12_) FieldBasedArrays.IntX12toS0 = EstimateIntermedVal(FieldBasedArrays.IntX12Arr, IndS0, sStp, dS0);
		if(Keys.IntX22toS0_ && Keys.IntX22_) FieldBasedArrays.IntX22toS0 = EstimateIntermedVal(FieldBasedArrays.IntX22Arr, IndS0, sStp, dS0);
		if(Keys.IntZ01toS0_ && Keys.IntZ01_) FieldBasedArrays.IntZ01toS0 = EstimateIntermedVal(FieldBasedArrays.IntZ01Arr, IndS0, sStp, dS0);
		if(Keys.IntZ02toS0_ && Keys.IntZ02_) FieldBasedArrays.IntZ02toS0 = EstimateIntermedVal(FieldBasedArrays.IntZ02Arr, IndS0, sStp, dS0);
		if(Keys.IntZ11toS0_ && Keys.IntZ11_) FieldBasedArrays.IntZ11toS0 = EstimateIntermedVal(FieldBasedArrays.IntZ11Arr, IndS0, sStp, dS0);
		if(Keys.IntZ12toS0_ && Keys.IntZ12_) FieldBasedArrays.IntZ12toS0 = EstimateIntermedVal(FieldBasedArrays.IntZ12Arr, IndS0, sStp, dS0);
		if(Keys.IntZ22toS0_ && Keys.IntZ22_) FieldBasedArrays.IntZ22toS0 = EstimateIntermedVal(FieldBasedArrays.IntZ22Arr, IndS0, sStp, dS0);
	}
	if(Keys.CheckIfNonEqNumIntegIsNecessary())
	{
		char CorrIntX01_IsNecessary = (Keys.IntX01_ && (!Keys.X1p_));
		char CorrIntX02_IsNecessary = (Keys.IntX02_ && (!Keys.X2p_));
		char CorrIntX11_IsNecessary = (Keys.IntX11_ && (!Keys.X1p_));
		char CorrIntX12_IsNecessary = (Keys.IntX12_ && ((!Keys.X1p_) || (!Keys.X2p_)));
		char CorrIntX22_IsNecessary = (Keys.IntX22_ && (!Keys.X2p_));
		char CorrIntZ01_IsNecessary = (Keys.IntZ01_ && (!Keys.Z1p_));
		char CorrIntZ02_IsNecessary = (Keys.IntZ02_ && (!Keys.Z2p_));
		char CorrIntZ11_IsNecessary = (Keys.IntZ11_ && (!Keys.Z1p_));
		char CorrIntZ12_IsNecessary = (Keys.IntZ12_ && ((!Keys.Z1p_) || (!Keys.Z2p_)));
		char CorrIntZ22_IsNecessary = (Keys.IntZ22_ && (!Keys.Z2p_));

        pIntX01 = FieldBasedArrays.IntX01Arr; pIntX02 = FieldBasedArrays.IntX02Arr;
        pIntX11 = FieldBasedArrays.IntX11Arr; pIntX12 = FieldBasedArrays.IntX12Arr; pIntX22 = FieldBasedArrays.IntX22Arr;
        pIntZ01 = FieldBasedArrays.IntZ01Arr; pIntZ02 = FieldBasedArrays.IntZ02Arr;
        pIntZ11 = FieldBasedArrays.IntZ11Arr; pIntZ12 = FieldBasedArrays.IntZ12Arr; pIntZ22 = FieldBasedArrays.IntZ22Arr;

		//for(int k=0; k<FieldBasedArrays.Ns; k++)
		for(long long k=0; k<FieldBasedArrays.Ns; k++)
		{
			if(CorrIntX01_IsNecessary) *(pIntX01++) -= FieldBasedArrays.IntX01toS0;
			if(CorrIntX02_IsNecessary) *(pIntX02++) -= FieldBasedArrays.IntX02toS0;
			if(CorrIntX11_IsNecessary) *(pIntX11++) -= FieldBasedArrays.IntX11toS0;
			if(CorrIntX12_IsNecessary) *(pIntX12++) -= FieldBasedArrays.IntX12toS0;
			if(CorrIntX22_IsNecessary) *(pIntX22++) -= FieldBasedArrays.IntX22toS0;
			if(CorrIntZ01_IsNecessary) *(pIntZ01++) -= FieldBasedArrays.IntZ01toS0;
			if(CorrIntZ02_IsNecessary) *(pIntZ02++) -= FieldBasedArrays.IntZ02toS0;
			if(CorrIntZ11_IsNecessary) *(pIntZ11++) -= FieldBasedArrays.IntZ11toS0;
			if(CorrIntZ12_IsNecessary) *(pIntZ12++) -= FieldBasedArrays.IntZ12toS0;
			if(CorrIntZ22_IsNecessary) *(pIntZ22++) -= FieldBasedArrays.IntZ22toS0;
		}
	}
}

//*************************************************************************

void srTTrjDat::CompTotalTrjData_FromTrj(srTFieldBasedArrayKeys& Keys, srTFieldBasedArrays& FieldBasedArrays)
{
	double s = FieldBasedArrays.sStart;
	double sStp = FieldBasedArrays.sStep;

	double *pBx = FieldBasedArrays.BxArr, *pBz = FieldBasedArrays.BzArr;
	double *pBtx = FieldBasedArrays.BtxArr, *pBtz = FieldBasedArrays.BtzArr;
	double *pX = FieldBasedArrays.XArr, *pZ = FieldBasedArrays.ZArr;
	double *pIntBtxE2 = FieldBasedArrays.IntBtxE2Arr, *pIntBtzE2 = FieldBasedArrays.IntBtzE2Arr;
	double *pdBxds = FieldBasedArrays.dBxdsArr, *pdBzds = FieldBasedArrays.dBzdsArr;

	//int Indx;
	long long Indx;
	double sr;
	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//for(int i=0; i<FieldBasedArrays.Ns; i++)
	for(long long i=0; i<FieldBasedArrays.Ns; i++)
	{
		//Indx = int((s - xTrjInData.Start)/xTrjInData.Step);
		Indx = (long long)((s - xTrjInData.Start)/xTrjInData.Step);
		if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*xTrjInData.Step;
		else if(Indx < xTrjInData.np - 3) ;
		else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
		else sr += 2*xTrjInData.Step;

		pB_Cf = *(BzPlnCf+Indx); pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pIntBtE2_Cf = *(IntBtx2PlnCf+Indx);
		if(Keys.dBzds_) *(pdBzds++) = *(pB_Cf+1) + sr*(*(pB_Cf+2)*2 + sr*(*(pB_Cf+3)*3));
		if(Keys.Bz_) *(pBz++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		if(Keys.Btx_) *(pBtx++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		if(Keys.X_) *(pX++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		if(Keys.IntBtxE2_) *(pIntBtxE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));

		//Indx = int((s - zTrjInData.Start)/zTrjInData.Step);
		Indx = (long long)((s - zTrjInData.Start)/zTrjInData.Step);
		if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*zTrjInData.Step;
		else if(Indx < zTrjInData.np - 3) ;
		else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
		else sr += 2*zTrjInData.Step;
	
		pB_Cf = *(BxPlnCf+Indx); pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pIntBtE2_Cf = *(IntBtz2PlnCf+Indx);
		if(Keys.dBxds_) *(pdBxds++) = *(pB_Cf+1) + sr*(*(pB_Cf+2)*2 + sr*(*(pB_Cf+3)*3));
		if(Keys.Bx_) *(pBx++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		if(Keys.Btz_) *(pBtz++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		if(Keys.Z_) *(pZ++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		if(Keys.IntBtzE2_) *(pIntBtzE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));

		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
{// Attention: length units are m at the input/output here !
	if(CompFromTrj) { CompTotalTrjData_FromTrj(sSt, sEn, Np, pBtx, pBtz, pX, pZ, pBx, pBz); return;}
	
	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	//double dxds0E2 = dxds0*dxds0;
	//double dzds0E2 = dzds0*dzds0;
	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;

	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;

		double *B_CfP, *Bt_CfP, *C_CfP;
		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx];

			*(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtx++) = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pX++) = (xCorr + BtxCorrForX*s + BufCrd);
		}
		else 
		{ 
			*(pBz++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtx++) = dxds0; 
			*(pX++) = (x0 + dxds0*s_mi_s0);
		}
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx];

			*(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtz++) = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pZ++) = (zCorr + BtzCorrForZ*s + BufCrd);
		}
		else 
		{ 
			*(pBx++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtz++) = dzds0; 
			*(pZ++) = (z0 + dzds0*s_mi_s0);
		}
		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pBx, double* pBz)
{
	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//int Indx;
	long long Indx;
	double sr;
	double *pB_Cf, *pBt_Cf, *pCrd_Cf;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//Indx = int((s - xTrjInData.Start)/xTrjInData.Step);
		Indx = (long long)((s - xTrjInData.Start)/xTrjInData.Step);
		if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*xTrjInData.Step;
		else if(Indx < xTrjInData.np - 3) ;
		else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
		else sr += 2*xTrjInData.Step;

		pB_Cf = *(BzPlnCf+Indx); pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx);
		*(pX++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtx++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		*(pBz++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));

		//Indx = int((s - zTrjInData.Start)/zTrjInData.Step);
		Indx = (long long)((s - zTrjInData.Start)/zTrjInData.Step);
		if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*zTrjInData.Step;
		else if(Indx < zTrjInData.np - 3) ;
		else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
		else sr += 2*zTrjInData.Step;

		pB_Cf = *(BxPlnCf+Indx); pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx);
		*(pZ++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtz++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		*(pBx++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));

		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz)
void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz)
{
	if(CompFromTrj) { CompTotalTrjData_FromTrj(sSt, sEn, Np, pBtx, pBtz, pX, pZ, pIntBtxE2, pIntBtzE2, pBx, pBz); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;
	double dxds0E2 = dxds0*dxds0, dzds0E2 = dzds0*dzds0;

	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;

		double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;
		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];

			*(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtx++) = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pX++) = xCorr + BtxCorrForX*s + BufCrd;
			*(pIntBtxE2++) = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{ 
			*(pBz++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtx++) = dxds0; *(pX++) = x0 + dxds0*s_mi_s0; *(pIntBtxE2++) = dxds0E2*s_mi_s0;
		}
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];

			*(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtz++) = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pZ++) = zCorr + BtzCorrForZ*s + BufCrd;
			*(pIntBtzE2++) = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{ 
			*(pBx++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtz++) = dzds0; *(pZ++) = z0 + dzds0*s_mi_s0; *(pIntBtzE2++) = dzds0E2*s_mi_s0;
		}
		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz)
void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz)
{
	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;
	double dxds0E2 = dxds0*dxds0, dzds0E2 = dzds0*dzds0;

	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//int Indx;
	long long Indx;
	double sr;
	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		if(VerFieldIsNotZero)
		{
			//Indx = int((s - xTrjInData.Start)/xTrjInData.Step);
			Indx = (long long)((s - xTrjInData.Start)/xTrjInData.Step);
			if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
			if(Indx < 0) Indx = 0;
			sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
			if(Indx < 2) sr -= (2 - Indx)*xTrjInData.Step;
			else if(Indx < xTrjInData.np - 3) ;
			else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
			else sr += 2*xTrjInData.Step;
			
			pB_Cf = *(BzPlnCf+Indx); pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pIntBtE2_Cf = *(IntBtx2PlnCf+Indx);
			*(pIntBtxE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
			*(pX++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
			*(pBtx++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
			*(pBz++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		}
		else
		{
			*(pBz++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtx++) = dxds0; *(pX++) = x0 + dxds0*s_mi_s0; *(pIntBtxE2++) = dxds0E2*s_mi_s0;
		}
		if(HorFieldIsNotZero)
		{
			//Indx = int((s - zTrjInData.Start)/zTrjInData.Step);
			Indx = (long long)((s - zTrjInData.Start)/zTrjInData.Step);
			if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
			if(Indx < 0) Indx = 0;
			sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
			if(Indx < 2) sr -= (2 - Indx)*zTrjInData.Step;
			else if(Indx < zTrjInData.np - 3) ;
			else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
			else sr += 2*zTrjInData.Step;
			
			pB_Cf = *(BxPlnCf+Indx); pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pIntBtE2_Cf = *(IntBtz2PlnCf+Indx);
			*(pIntBtzE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
			*(pZ++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
			*(pBtz++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
			*(pBx++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		}
		else
		{
			*(pBx++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtz++) = dzds0; *(pZ++) = z0 + dzds0*s_mi_s0; *(pIntBtzE2++) = dzds0E2*s_mi_s0;
		}
		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
void srTTrjDat::CompTotalTrjData(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
{
	if(CompFromTrj) { CompTotalTrjData_FromTrj(sSt, sEn, Np, pBtx, pBtz, pX, pZ, pIntBtxE2, pIntBtzE2, pBx, pBz, pdBxds, pdBzds); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	double dxds0E2 = dxds0*dxds0;
	double dzds0E2 = dzds0*dzds0;

	double s = sSt, sStp = (Np > 1)? (sEn - sSt)/(Np - 1) : 0.;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;

		double *B_CfP, *Bt_CfP, *C_CfP, *IntBt2_CfP;

		if(VerFieldIsNotZero)
		{
			B_CfP = BzPlnCf[Indx]; Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx]; IntBt2_CfP = IntBtx2PlnCf[Indx];

			*(pdBzds++) = *(B_CfP+1) + smsb*(2.*(*(B_CfP+2)) + 3.*smsb*(*(B_CfP+3)));
			*(pBz++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtx++) = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pX++) = xCorr + BtxCorrForX*s + BufCrd;
			*(pIntBtxE2++) = IntBtxE2Corr + BtxCorrForXe2*s + 2.*BtxCorrForX*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{ 
			*(pdBzds++) = 0.; *(pBz++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtx++) = dxds0; *(pX++) = x0 + dxds0*s_mi_s0; *(pIntBtxE2++) = dxds0E2*s_mi_s0;
		}
		if(HorFieldIsNotZero)
		{
			B_CfP = BxPlnCf[Indx]; Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx]; IntBt2_CfP = IntBtz2PlnCf[Indx];

			*(pdBxds++) = *(B_CfP+1) + smsb*(2.*(*(B_CfP+2)) + 3.*smsb*(*(B_CfP+3)));
			*(pBx++) = *B_CfP + smsb*(*(B_CfP+1) + smsb*(*(B_CfP+2) + smsb*(*(B_CfP+3))));
			*(pBtz++) = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pZ++) = zCorr + BtzCorrForZ*s + BufCrd;
			*(pIntBtzE2++) = IntBtzE2Corr + BtzCorrForZe2*s + 2.*BtzCorrForZ*BufCrd + BetaNormConstE2*(Pol09(smsb, IntBt2_CfP+1) + *IntBt2_CfP);
		}
		else 
		{ 
			*(pdBxds++) = 0.; *(pBx++) = 0.;
			double s_mi_s0 = s - s0;
			*(pBtz++) = dzds0; *(pZ++) = z0 + dzds0*s_mi_s0; *(pIntBtzE2++) = dzds0E2*s_mi_s0;
		}
		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
void srTTrjDat::CompTotalTrjData_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, double* pIntBtxE2, double* pIntBtzE2, double* pBx, double* pBz, double* pdBxds, double* pdBzds)
{
	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//int Indx;
	long long Indx;
	double sr;
	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//Indx = int((s - xTrjInData.Start)/xTrjInData.Step);
		Indx = (long long)((s - xTrjInData.Start)/xTrjInData.Step);
		if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*xTrjInData.Step;
		else if(Indx < xTrjInData.np - 3) ;
		else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
		else sr += 2*xTrjInData.Step;

		pB_Cf = *(BzPlnCf+Indx); pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx); pIntBtE2_Cf = *(IntBtx2PlnCf+Indx);
		*(pIntBtxE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
		*(pX++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtx++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		*(pBz++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		*(pdBzds++) = *(pB_Cf+1) + sr*(*(pB_Cf+2)*2 + sr*(*(pB_Cf+3)*3));

		//Indx = int((s - zTrjInData.Start)/zTrjInData.Step);
		Indx = (long long)((s - zTrjInData.Start)/zTrjInData.Step);
		if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*zTrjInData.Step;
		else if(Indx < zTrjInData.np - 3) ;
		else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
		else sr += 2*zTrjInData.Step;

		pB_Cf = *(BxPlnCf+Indx); pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx); pIntBtE2_Cf = *(IntBtz2PlnCf+Indx);
		*(pIntBtzE2++) = *pIntBtE2_Cf + sr*(*(pIntBtE2_Cf+1) + sr*(*(pIntBtE2_Cf+2) + sr*(*(pIntBtE2_Cf+3) + sr*(*(pIntBtE2_Cf+4) + sr*(*(pIntBtE2_Cf+5))))));
		*(pZ++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtz++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));
		*(pBx++) = *pB_Cf + sr*(*(pB_Cf+1) + sr*(*(pB_Cf+2) + sr*(*(pB_Cf+3))));
		*(pdBxds++) = *(pB_Cf+1) + sr*(*(pB_Cf+2)*2 + sr*(*(pB_Cf+3)*3));

		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjDataTrjDisp(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn)
void srTTrjDat::CompTotalTrjDataTrjDisp(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn)
{
	if(CompFromTrj) { CompTotalTrjDataTrjDisp_FromTrj(sSt, sEn, Np, pBtx, pBtz, pX, pZ, DistUn); return;}

	double &dxds0 = EbmDat.dxds0, &x0 = EbmDat.x0, &dzds0 = EbmDat.dzds0, &z0 = EbmDat.z0, &s0 = EbmDat.s0;

	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//int Indx = int((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		long long Indx = (long long)((s - sStart)/sStep); if(Indx >= LenFieldData - 1) Indx = LenFieldData - 2;
		double sb = sStart + Indx*sStep;
		double smsb = s - sb;

		double *Bt_CfP, *C_CfP;
		if(VerFieldIsNotZero)
		{
			Bt_CfP = BtxPlnCf[Indx]; C_CfP = xPlnCf[Indx];

			*(pBtx++) = BtxCorr + BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pX++) = xCorr + BtxCorrForX*s + BufCrd;
		}
		else 
		{ 
			double s_mi_s0 = s - s0;
			*(pBtx++) = dxds0; 
			*(pX++) = x0 + dxds0*s_mi_s0;
		}
		if(HorFieldIsNotZero)
		{
			Bt_CfP = BtzPlnCf[Indx]; C_CfP = zPlnCf[Indx];

			*(pBtz++) = BtzCorr - BetaNormConst*(Pol04(smsb, Bt_CfP+1) + *Bt_CfP);
			double BufCrd = -BetaNormConst*(Pol05(smsb, C_CfP+1) + *C_CfP);
			*(pZ++) = zCorr + BtzCorrForZ*s + BufCrd;
		}
		else 
		{ 
			double s_mi_s0 = s - s0;
			*(pBtz++) = dzds0; 
			*(pZ++) = z0 + dzds0*s_mi_s0;
		}
		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTotalTrjDataTrjDisp_FromTrj(double sSt, double sEn, long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn)
void srTTrjDat::CompTotalTrjDataTrjDisp_FromTrj(double sSt, double sEn, long long Np, double* pBtx, double* pBtz, double* pX, double* pZ, char DistUn)
{
	double s = sSt, sStp = (Np>1)? (sEn - sSt)/(Np - 1) : 0.;
	//int Indx;
	long long Indx;
	double sr;
	double *pBt_Cf, *pCrd_Cf;
	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		//Indx = int((s - xTrjInData.Start)/xTrjInData.Step);
		Indx = (long long)((s - xTrjInData.Start)/xTrjInData.Step);
		if(Indx >= xTrjInData.np - 1) Indx = xTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (xTrjInData.Start + xTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*xTrjInData.Step;
		else if(Indx < xTrjInData.np - 3) ;
		else if(Indx < xTrjInData.np - 2) sr += xTrjInData.Step;
		else sr += 2*xTrjInData.Step;

		pBt_Cf = *(BtxPlnCf+Indx); pCrd_Cf = *(xPlnCf+Indx);
		*(pX++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtx++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));

		//Indx = int((s - zTrjInData.Start)/zTrjInData.Step);
		Indx = (long long)((s - zTrjInData.Start)/zTrjInData.Step);
		if(Indx >= zTrjInData.np - 1) Indx = zTrjInData.np - 2;
		if(Indx < 0) Indx = 0;
		sr = s - (zTrjInData.Start + zTrjInData.Step*Indx);
		if(Indx < 2) sr -= (2 - Indx)*zTrjInData.Step;
		else if(Indx < zTrjInData.np - 3) ;
		else if(Indx < zTrjInData.np - 2) sr += zTrjInData.Step;
		else sr += 2*zTrjInData.Step;

		pBt_Cf = *(BtzPlnCf+Indx); pCrd_Cf = *(zPlnCf+Indx);
		*(pZ++) = *pCrd_Cf + sr*(*(pCrd_Cf+1) + sr*(*(pCrd_Cf+2) + sr*(*(pCrd_Cf+3) + sr*(*(pCrd_Cf+4) + sr*(*(pCrd_Cf+5))))));
		*(pBtz++) = *pBt_Cf + sr*(*(pBt_Cf+1) + sr*(*(pBt_Cf+2) + sr*(*(pBt_Cf+3) + sr*(*(pBt_Cf+4)))));

		s += sStp;
	}
}

//*************************************************************************

//void srTTrjDat::CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep)
void srTTrjDat::CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, long long ns, double sStart, double sStep)
{
	if((pOutBtxData == 0) && (pOutXData == 0) && (pOutBtzData == 0) && (pOutZData == 0)) return;

	int res = 0;
	if(res = ComputeInterpolatingStructure()) throw res;

	double sFi = sStart + sStep*(LenFieldData - 1);
	//in the future: use 
	//int ns, double sStart, double sStep

	char DistUnits = 1; // m for coord.
	CompTotalTrjDataTrjDisp(sStart, sFi, LenFieldData, pOutBtxData, pOutBtzData, pOutXData, pOutZData, DistUnits);
}

//*************************************************************************

//int srTTrjDat::SetupLimitsByAnalizingField(char LongIntType, double& sStartFld, double& sStepFld, long& Ns, int& OutNperTot, int& OutNperLeft)
int srTTrjDat::SetupLimitsByAnalizingField(char LongIntType, double& sStartFld, double& sStepFld, long long& Ns, long long& OutNperTot, long long& OutNperLeft)
{
	const double RelTolForField = 1.E-07; // To steer
	const double MinStepSteerCoef = 0.5; // To steer
	const int NsMin = 11;

	int AmOfExtrem;
	double AbsMax;
	CountFieldExtrem(AmOfExtrem, AbsMax);
	if(AmOfExtrem <= 0) AmOfExtrem = 1;
	double sStartLoc, sFinLoc;
	double AbsTolForField = RelTolForField*AbsMax;
	FindFieldLimitsBasedOnTolValue(AbsTolForField, sStartLoc, sFinLoc);

	double Rmin = 3.3*EbmDat.Energy/AbsMax;
	double dsMin = MinStepSteerCoef*Rmin/EbmDat.Gamma;
	double sRange = sFinLoc - sStartLoc;
	//Ns = long(sRange/dsMin);
	Ns = (long long)(sRange/dsMin);
	if(Ns < NsMin) Ns = NsMin;
	if(((Ns >> 1) << 1) == Ns) Ns++; // Ensure odd

	sStartFld = sStartLoc;
	sStepFld = sRange/double(Ns - 1);

	OutNperTot = NperTot;
	OutNperLeft = NperLeft;
	return 0;
}

//*************************************************************************

void srTTrjDat::CountFieldExtrem(int& AmOfExtrem, double& AbsMax)
{
	AbsMax = 0.;
	int ExtremCountBx = 0, ExtremCountBz = 0;

	srTFunDer *tBxInData = BxInData;
	srTFunDer *tBzInData = BzInData;

	char DerSignBx = 0, DerSignBz = 0;
	if(HorFieldIsNotZero) 
	{
		DerSignBx = (tBxInData->dfds < 0.)? -1 : 1;
		double AbsFld = ::fabs(tBxInData->f);
		if(AbsMax < AbsFld) AbsMax = AbsFld;
		tBxInData++;
	}
	if(VerFieldIsNotZero) 
	{
		DerSignBz = (tBzInData->dfds < 0.)? -1 : 1;
		double AbsFld = ::fabs(tBzInData->f);
		if(AbsMax < AbsFld) AbsMax = AbsFld;
		tBzInData++;
	}

	//for(long i=1; i<LenFieldData; i++)
	for(long long i=1; i<LenFieldData; i++)
	{
		if(HorFieldIsNotZero) 
		{
			char NewDerSignBx = (tBxInData->dfds < 0.)? -1 : 1;
			if(NewDerSignBx != DerSignBx) ExtremCountBx++;
			DerSignBx = NewDerSignBx;
			double AbsFld = ::fabs(tBxInData->f);
			if(AbsMax < AbsFld) AbsMax = AbsFld;
			tBxInData++;
		}
		if(VerFieldIsNotZero) 
		{
			char NewDerSignBz = (tBzInData->dfds < 0.)? -1 : 1;
			if(NewDerSignBz != DerSignBz) ExtremCountBz++;
			DerSignBz = NewDerSignBz;
			double AbsFld = ::fabs(tBzInData->f);
			if(AbsMax < AbsFld) AbsMax = AbsFld;
			tBzInData++;
		}
	}
	AmOfExtrem = (ExtremCountBz > ExtremCountBx)? ExtremCountBz : ExtremCountBx;
}

//*************************************************************************

void srTTrjDat::FindFieldLimitsBasedOnTolValue(double AbsTolForField, double& sStartLoc, double& sFinLoc)
{
	//long LenFieldData_mi_1 = LenFieldData - 1;
	long long LenFieldData_mi_1 = LenFieldData - 1;
	srTFunDer *tBxInData = BxInData, *trBxInData = BxInData + LenFieldData_mi_1;
	srTFunDer *tBzInData = BzInData, *trBzInData = BzInData + LenFieldData_mi_1;

	//long iSt = -1, iFi = -1;
	long long iSt = -1, iFi = -1;
	//for(long i=0; i<LenFieldData; i++)
	for(long long i=0; i<LenFieldData; i++)
	{
		if(HorFieldIsNotZero) 
		{
			if(iSt < 0)
			{
				double AbsFld = ::fabs(tBxInData->f);
				if(AbsFld > AbsTolForField) iSt = i - 1;
			}
			if(iFi < 0)
			{
				double AbsFld = ::fabs(trBxInData->f);
				if(AbsFld > AbsTolForField) iFi = LenFieldData_mi_1 - i + 1;
			}
			tBxInData++; trBxInData--;
		}
		if(VerFieldIsNotZero) 
		{
			if(iSt < 0)
			{
				double AbsFld = ::fabs(tBzInData->f);
				if(AbsFld > AbsTolForField) iSt = i - 1;
			}
			if(iFi < 0)
			{
				double AbsFld = ::fabs(trBzInData->f);
				if(AbsFld > AbsTolForField) iFi = LenFieldData_mi_1 - i + 1;
			}
			tBzInData++; trBzInData--;
		}
		if((iSt >= 0) && (iFi >= 0)) break;
	}
	if(iSt < 0) iSt = 0;
	if((iFi < 0) || (iFi > LenFieldData_mi_1)) iFi = LenFieldData_mi_1;

	sStartLoc = sStart + iSt*sStep;
	sFinLoc = sStart + iFi*sStep;
}

//*************************************************************************

void srTTrjDat::AnalizeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ)
{
	FieldIsSymOverX = FieldIsSymOverZ = 0;
	if(!HorFieldIsNotZero) FieldIsSymOverZ = 1;
	if(!VerFieldIsNotZero) FieldIsSymOverX = 1;
}


//*************************************************************************

int srTTrjDat::CheckAndSetupTrajectoryLimits()
{// LenFieldData, sStart, sStep, Inv_Step;
	if((xTrjInData.pData == 0) || (zTrjInData.pData == 0)) return TRJ_CMPN_WERE_NOT_SETUP;

	double xTrjEnd = xTrjInData.Start + xTrjInData.Step*(xTrjInData.np - 1);
	double zTrjEnd = zTrjInData.Start + zTrjInData.Step*(zTrjInData.np - 1);
	if((xTrjEnd < zTrjInData.Start) || (zTrjEnd < xTrjInData.Start)) return CONTRADICTORY_TRJ_CMPN_DEFINITION;

	double sEnd;
	if(xTrjInData.Start >= zTrjInData.Start)
	{
		sStart = xTrjInData.Start;
		sStep = xTrjInData.Step;

		sEnd = (xTrjEnd <= zTrjEnd)? xTrjEnd : zTrjEnd;
	}
	else
	{
		sStart = zTrjInData.Start;
		sStep = zTrjInData.Step;

		sEnd = (zTrjEnd <= xTrjEnd)? zTrjEnd : xTrjEnd;
	}
	//LenFieldData = long((sEnd - sStart)/sStep + 1.e-04) + 1;
	LenFieldData = (long long)((sEnd - sStart)/sStep + 1.e-04) + 1;
	Inv_Step = 1./sStep;

	return 0;
}

//*************************************************************************

int srTTrjDat::SetupSourcePointFromTrajectory()
{// sets up s0, x0, dxds0, z0, dzds0, ...
	int result;
	if(result = CheckAndSetupTrajectoryLimits()) return result;

	EbmDat.s0 = sStart + sStep*(LenFieldData >> 1); // Make more analysis

	double Bz, Bx;
	TrjCoordAngField(EbmDat.s0, 'x', EbmDat.x0, EbmDat.dxds0, Bz);
	TrjCoordAngField(EbmDat.s0, 'z', EbmDat.z0, EbmDat.dzds0, Bx);

	EbmDat.Mxx = EbmDat.Mxxp = EbmDat.Mxpxp = 0.;
	EbmDat.Mzz = EbmDat.Mzzp = EbmDat.Mzpzp = 0.;
	EbmDat.Mxz = EbmDat.Mxpz = EbmDat.Mxzp = EbmDat.Mxpzp = 0.;
	EbmDat.Mee = EbmDat.SigmaRelE = 0.;
	return 0;
}

//*************************************************************************

int srTTrjDat::ComputeOneQuadPhaseTermFromTrj(char x_or_z)
{
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData;
	DOUBLE *pIntBtE2Arr = (x_or_z == 'x')? IntBtxE2Arr : IntBtzE2Arr;
	DOUBLE *tIntBtE2Arr = pIntBtE2Arr;

	if(TrjInData.pData == 0) return TRJ_CMPN_WERE_NOT_SETUP;

	double sStep = TrjInData.Step;
	double sStepSmall = 0.5*sStep;
	double StepMult = 0.333333333333*sStepSmall;
	double Crd, Ang, B;
	double f0 = 0., f1, f2;

	*(tIntBtE2Arr++) = 0.;
	double CurVal = 0.;
	double sMid = TrjInData.Start + sStepSmall;
	double sLast = TrjInData.Start + TrjInData.Step;
	//for(long i=1; i<TrjInData.np; i++)
	for(long long i=1; i<TrjInData.np; i++)
	{
		TrjCoordAngField(sMid, x_or_z, Crd, Ang, B);
		f1 = Ang*Ang;
		TrjCoordAngField(sLast, x_or_z, Crd, Ang, B);
		f2 = Ang*Ang;

		CurVal += StepMult*(f0 + 4.*f1 + f2);
		*(tIntBtE2Arr++) = CurVal;
		f0 = f2; 
		sMid += sStep; sLast += sStep;
	}

	double dBdsDummy, BDummy, BtDummy, CrdDummy, IntBtE2;
	CompTrjDataAndFieldWithDerAtPoint_FromTrjInitial(x_or_z, EbmDat.s0, dBdsDummy, BDummy, BtDummy, CrdDummy, IntBtE2);
	double AddCorr = -IntBtE2;
	tIntBtE2Arr = pIntBtE2Arr;
	//for(long j=0; j<TrjInData.np; j++) 
	for(long long j=0; j<TrjInData.np; j++) 
	{
		*tIntBtE2Arr += AddCorr;
		tIntBtE2Arr++;
	}
	return 0;
}

//*************************************************************************

int srTTrjDat::ComputeQuadPhaseTermsFromTrj(const SRWLPrtTrj& trj)
{
	bool xIsDefined = (trj.arXp != 0) && (trj.np > 0);
	bool yIsDefined = (trj.arYp != 0) && (trj.np > 0);
	if(!(xIsDefined || yIsDefined)) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;

	//long startInd = 0;
	long long startInd = 0;
	double sRel = 0;
	FindOffestAndRelArgFromTrj(EbmDat.s0, trj, startInd, sRel);

	//int np_mi_1 = trj.np - 1;
	//int np_mi_2 = np_mi_1 - 1;
	long long np_mi_1 = trj.np - 1;
	long long np_mi_2 = np_mi_1 - 1;
	double sStep = (np_mi_1 > 0)? (trj.ctEnd - trj.ctStart)/np_mi_1 : 0;
	double sStepSmall = 0.5*sStep;
	double stepMult = 0.333333333333*sStepSmall;
	double f0 = 0., f1, f2;
	double curVal = 0.;
	double ang, dummyDer1, dummyDer2, dummyDer3, IntBtE2;

	if(xIsDefined)
	{
		*IntBtxE2Arr = 0;
		DOUBLE *tIntBtE2 = IntBtxE2Arr + 1;
		double *p = trj.arXp;
		double *t = p + 1;
		int iMidCase = -1;
		f0 = 0.; curVal = 0.;
		//for(long i=1; i<=np_mi_1; i++)
		for(long long i=1; i<=np_mi_1; i++)
		{
			ang = CGenMathInterp::InterpCubHalfStep(p, iMidCase);
			f1 = ang*ang;
			ang = *(t++);
			f2 = ang*ang;
			curVal += stepMult*(f0 + 4.*f1 + f2);
			*(tIntBtE2++) = curVal;
			f0 = f2; 

			if((i > 1) && (i < np_mi_2)) p++;
			iMidCase = (i < np_mi_1)? 0 : 1;
		}

		InterpFuncAndDerivs(sStep, sRel, IntBtxE2Arr + startInd, IntBtE2, dummyDer1, dummyDer2, dummyDer3);
		tIntBtE2 = IntBtxE2Arr;
		//for(int j=0; j<=np_mi_1; j++) *(tIntBtE2++) -= IntBtE2;
		for(long long j=0; j<=np_mi_1; j++) *(tIntBtE2++) -= IntBtE2;
	}
	if(yIsDefined)
	{
		*IntBtzE2Arr = 0;
		DOUBLE *tIntBtE2 = IntBtzE2Arr + 1;
		double *p = trj.arYp;
		double *t = p + 1;
		int iMidCase = -1;
		f0 = 0.; curVal = 0.;
		//for(long i=1; i<=np_mi_1; i++)
		for(long long i=1; i<=np_mi_1; i++)
		{
			ang = CGenMathInterp::InterpCubHalfStep(p, iMidCase);
			f1 = ang*ang;
			ang = *(t++);
			f2 = ang*ang;
			curVal += stepMult*(f0 + 4.*f1 + f2);
			*(tIntBtE2++) = curVal;
			f0 = f2; 

			if((i > 1) && (i < np_mi_2)) p++;
			iMidCase = (i < np_mi_1)? 0 : 1;
		}

		InterpFuncAndDerivs(sStep, sRel, IntBtzE2Arr + startInd, IntBtE2, dummyDer1, dummyDer2, dummyDer3);
		tIntBtE2 = IntBtzE2Arr;
		//for(int j=0; j<=np_mi_1; j++) *(tIntBtE2++) -= IntBtE2;
		for(long long j=0; j<=np_mi_1; j++) *(tIntBtE2++) -= IntBtE2;
	}
	return 0;
}

//*************************************************************************

int srTTrjDat::ComputeInterpolatingStructure_FromTrj1D(char x_or_z)
{
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData;
	if(TrjInData.pData == 0) return TRJ_CMPN_WERE_NOT_SETUP;
	DOUBLE *pIntBtE2Arr = (x_or_z == 'x')? IntBtxE2Arr : IntBtzE2Arr;

	double BMult = InvBetaNormConst;
	if(x_or_z != 'x') BMult = -BMult;

	int angZeroCount = 0; //OC220112
	double *tAuxCrd = TrjInData.pData, varCrd = (*(TrjInData.pData + 1)) - (*(TrjInData.pData));

	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//for(long i=0; i<(TrjInData.np - 1); i++) 
	for(long long i=0; i<(TrjInData.np - 1); i++) 
	{
		if(x_or_z == 'x')
		{
			pB_Cf = *(BzPlnCf+i); pBt_Cf = *(BtxPlnCf+i); pCrd_Cf = *(xPlnCf+i); pIntBtE2_Cf = *(IntBtx2PlnCf+i);
		}
		else
		{
			pB_Cf = *(BxPlnCf+i); pBt_Cf = *(BtzPlnCf+i); pCrd_Cf = *(zPlnCf+i); pIntBtE2_Cf = *(IntBtz2PlnCf+i);
		}
		
		//long Offset;
		long long Offset;
		if(i < 2) Offset = 0;
		else if(i < TrjInData.np - 3) Offset = i - 2;
		else if(i < TrjInData.np - 2) Offset = i - 3;
		else Offset = i - 4;
		
		DOUBLE *pf0 = pIntBtE2Arr + Offset;
		CoefsPol5thOrder(TrjInData.Step, pf0, pIntBtE2_Cf);
		
		pf0 = TrjInData.pData + Offset;
		CoefsPol5thOrder(TrjInData.Step, pf0, pCrd_Cf);
		
		*pBt_Cf = *(pCrd_Cf+1); *(pBt_Cf+1) = 2.*(*(pCrd_Cf+2)); *(pBt_Cf+2) = 3.*(*(pCrd_Cf+3)); *(pBt_Cf+3) = 4.*(*(pCrd_Cf+4)); *(pBt_Cf+4) = 5.*(*(pCrd_Cf+5));
		*pB_Cf = (*(pBt_Cf+1))*BMult; *(pB_Cf+1) = 2.*(*(pBt_Cf+2))*BMult; *(pB_Cf+2) = 3.*(*(pBt_Cf+3))*BMult; *(pB_Cf+3) = 4.*(*(pBt_Cf+4))*BMult;

		double newVarCrd = (*(tAuxCrd + 1)) - (*tAuxCrd); //OC220112
		if(newVarCrd*varCrd < 0) angZeroCount++;
		tAuxCrd++;
		varCrd = newVarCrd;
	}

	if(angZeroCount > 0) //OC220112
	{
		if(m_estimMinNpForRadInteg < angZeroCount) m_estimMinNpForRadInteg = angZeroCount;
	}

	return 0;
}

//*************************************************************************

int srTTrjDat::ComputeInterpolatingStructureFromTrj1D(char x_or_z, const SRWLPrtTrj& trj)
{//Assumes right-hand frame with Z being longitudinal direction
	//srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData;
	
	double *pCrd = (x_or_z == 'x')? trj.arX : trj.arY;
	if(pCrd == 0) return TRJ_CMPN_WERE_NOT_SETUP;
	DOUBLE *pIntBtE2Arr = (x_or_z == 'x')? IntBtxE2Arr : IntBtzE2Arr;

	double BMult = InvBetaNormConst;
	if(x_or_z == 'x') BMult = -BMult;

	int angZeroCount = 0; //OC220112
	double *tAuxCrd = pCrd, varCrd = (*(pCrd + 1)) - (*pCrd);

	double *pB_Cf, *pBt_Cf, *pCrd_Cf, *pIntBtE2_Cf;
	//for(long i=0; i<(trj.np - 1); i++) 
	for(long long i=0; i<(trj.np - 1); i++) 
	{
		if(x_or_z == 'x')
		{
			pB_Cf = *(BzPlnCf+i); pBt_Cf = *(BtxPlnCf+i); pCrd_Cf = *(xPlnCf+i); pIntBtE2_Cf = *(IntBtx2PlnCf+i);
		}
		else
		{
			pB_Cf = *(BxPlnCf+i); pBt_Cf = *(BtzPlnCf+i); pCrd_Cf = *(zPlnCf+i); pIntBtE2_Cf = *(IntBtz2PlnCf+i);
		}

		//long Offset;
		long long Offset;
		if(i < 2) Offset = 0;
		else if(i < trj.np - 3) Offset = i - 2;
		else if(i < trj.np - 2) Offset = i - 3;
		else Offset = i - 4;
		
		double sStep = (trj.ctEnd - trj.ctStart)/(trj.np - 1);
		DOUBLE *pf0 = pIntBtE2Arr + Offset;
		CoefsPol5thOrder(sStep, pf0, pIntBtE2_Cf);

		pf0 = pCrd + Offset;
		CoefsPol5thOrder(sStep, pf0, pCrd_Cf);

		*pBt_Cf = *(pCrd_Cf+1); *(pBt_Cf+1) = 2.*(*(pCrd_Cf+2)); *(pBt_Cf+2) = 3.*(*(pCrd_Cf+3)); *(pBt_Cf+3) = 4.*(*(pCrd_Cf+4)); *(pBt_Cf+4) = 5.*(*(pCrd_Cf+5));
		*pB_Cf = (*(pBt_Cf+1))*BMult; *(pB_Cf+1) = 2.*(*(pBt_Cf+2))*BMult; *(pB_Cf+2) = 3.*(*(pBt_Cf+3))*BMult; *(pB_Cf+3) = 4.*(*(pBt_Cf+4))*BMult;

		double newVarCrd = (*(tAuxCrd + 1)) - (*tAuxCrd); //OC220112
		if(newVarCrd*varCrd < 0) angZeroCount++;
		tAuxCrd++;
		varCrd = newVarCrd;
	}

	if(angZeroCount > 0) //OC220112
	{
		if(m_estimMinNpForRadInteg < angZeroCount) m_estimMinNpForRadInteg = angZeroCount;
	}

	return 0;
}

//*************************************************************************

int srTTrjDat::FieldComponIsZero_FromTrj(char x_or_z)
{// x_or_z defines Trajectory component !!!
	const double AbsTrjZeroTol = 1.e-10; // m; To steer
	
	srTWaveAccessDataD1D &TrjInData = (x_or_z == 'x')? xTrjInData : zTrjInData;
	if(TrjInData.pData == 0) return TRJ_CMPN_WERE_NOT_SETUP;

	//long iBuf = TrjInData.np - 2;
	long long iBuf = TrjInData.np - 2;
	double f1 = *(TrjInData.pData + 1), f2 = *(TrjInData.pData + iBuf);
	double s1 = TrjInData.Start + TrjInData.Step, s2 = TrjInData.Start + TrjInData.Step*iBuf;
	double Inv_ds = 1./(s1 - s2);
	double a = (f1 - f2)*Inv_ds, b = (f2*s1 - f1*s2)*Inv_ds;

	char FieldIsNotZero = 0;
	DOUBLE *tData = TrjInData.pData;
	double s = TrjInData.Start;
	//for(long i=0; i<TrjInData.np; i++)
	for(long long i=0; i<TrjInData.np; i++)
	{
		double fPredict = s*a + b;
		if(::fabs(fPredict - *tData) > AbsTrjZeroTol) { FieldIsNotZero = 1; break;}

		s += TrjInData.Step; tData++;
	}

	return !FieldIsNotZero;
}

//*************************************************************************

void srTTrjDat::CountFieldExtremums()
{// Fills long AmOfMaxInBx, AmOfMaxInBz; //Auxiliary information on magn. field
	if(HorFieldIsNotZero)
	{
		//long ExtremCount = 0;
		long long ExtremCount = 0;
		srTFunDer* tBxInData = BxInData + 1;
		//for(int i=1; i<LenFieldData; i++)
		for(long long i=1; i<LenFieldData; i++)
		{
			if(((tBxInData-1)->dfds)*(tBxInData->dfds) < 0) ExtremCount++;
			tBxInData++;
		}
		AmOfExtremInBx = ExtremCount;
	}
	else AmOfExtremInBx = 0;
	if(VerFieldIsNotZero)
	{
		//long ExtremCount = 0;
		long long ExtremCount = 0;
		srTFunDer* tBzInData = BzInData + 1;
		//for(int i=1; i<LenFieldData; i++)
		for(long long i=1; i<LenFieldData; i++)
		{
			if(((tBzInData-1)->dfds)*(tBzInData->dfds) < 0) ExtremCount++;
			tBzInData++;
		}
		AmOfExtremInBz = ExtremCount;
	}
	else AmOfExtremInBz = 0;
}

//*************************************************************************

//srTTrjDat::srTTrjDat(srTEbmDat* pEbmDat, srTMagElem* pMagElem)
srTTrjDat::srTTrjDat(srTEbmDat* pEbmDat, srTMagFldTrUnif* pMagElem)
{
	Initialize();
	if(pEbmDat != 0) EbmDat = *pEbmDat;
	pMagElem->SetupTrjDat(this);
}

//*************************************************************************

srTTrjDat::srTTrjDat(SRWLPrtTrj* pTrj)
{
	if(pTrj == 0) throw SRWL_INCORRECT_PARAM_FOR_SR_COMP;

	Initialize();

	const double elecEn0 = 0.51099890221e-03; //[GeV]
	SRWLParticle &part = pTrj->partInitCond;
	double arMom1[] = {(part.gamma)*(part.relE0)*elecEn0, part.x, part.xp, part.y, part.yp, part.z};
	srTEbmDat elecBeam(1., 1., arMom1, 6, 0, 0, part.z, part.nq); //this doesn't set correct current!
	EbmDat = elecBeam;

	LenFieldData = pTrj->np;
	sStart = pTrj->ctStart + pTrj->partInitCond.z;
	sStep = (pTrj->ctEnd - pTrj->ctStart)/(LenFieldData - 1);
	//OC150815 (commented-out)
	//if(pTrj->arZ != 0)
	//{//use tabulated longitudinal position
	//	sStart = pTrj->arZ[0];
	//	sStep = (pTrj->arZ[LenFieldData - 1] - sStart)/(LenFieldData - 1);
	//}

	CheckFromTrjIfFieldCompAreZero(*pTrj, HorFieldIsNotZero, VerFieldIsNotZero);

	int result = ComputeInterpolatingStructureFromTrj(pTrj);
	CompFromTrj = 1; //?

	xTrjInData.np = pTrj->np;
	xTrjInData.Start = sStart;
	xTrjInData.Step = sStep;
	xTrjInData.InvStep = 1./sStep;

	zTrjInData.np = pTrj->np;
	zTrjInData.Start = sStart;
	zTrjInData.Step = sStep;
	zTrjInData.InvStep = xTrjInData.InvStep;

	if(result) throw result;
}

//*************************************************************************

void srTTrjDat::CheckFromTrjIfFieldCompAreZero(SRWLPrtTrj& trj, short& horFieldIsNotZero, short& verFieldIsNotZero)
{
	horFieldIsNotZero = 0; verFieldIsNotZero = 0;

	double *t_arX = trj.arX, *t_arXp = trj.arXp, *t_arY = trj.arY, *t_arYp = trj.arYp;
	//for(int i=0; i<trj.np; i++)
	for(long long i=0; i<trj.np; i++)
	{
		//if((*(t_arX++) != 0.) || (*(t_arXp++) != 0.)) horFieldIsNotZero = 1;
		//if((*(t_arY++) != 0.) || (*(t_arYp++) != 0.)) verFieldIsNotZero = 1;
		if((*(t_arX++) != 0.) || (*(t_arXp++) != 0.)) verFieldIsNotZero = 1; //OC_GIANLUCA
		if((*(t_arY++) != 0.) || (*(t_arYp++) != 0.)) horFieldIsNotZero = 1; //OC_GIANLUCA
		if(horFieldIsNotZero && verFieldIsNotZero) break;
	}
}

//*************************************************************************

int srTTrjDat::ComputeInterpolatingStructureFromTrj(SRWLPrtTrj* pTrj)
{
	if(pTrj == 0) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;
	bool trjIsDefined = (((pTrj->arX != 0) && (pTrj->arXp != 0)) || ((pTrj->arY != 0) && (pTrj->arYp != 0))) && (pTrj->np > 0);
	if(!trjIsDefined) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;

	int result = 0;
	if(result = AllocateQuadPhaseTermsArrFromTrj(*pTrj)) return result;
	if(result = ComputeQuadPhaseTermsFromTrj(*pTrj)) return result;
	if(result = AllocateMemoryForCfsFromTrj(pTrj->np)) return result;

	CompBetaNormConst();
	//xTrjInData.SetupBufVars(); zTrjInData.SetupBufVars();

	m_estimMinNpForRadInteg = -1;

	if(result = ComputeInterpolatingStructureFromTrj1D('x', *pTrj)) return result;
	if(result = ComputeInterpolatingStructureFromTrj1D('y', *pTrj)) return result;

	DeallocateQuadPhaseTermsArr();
	LastCompWasNotOK = 0; 
	return 0;
}

//*************************************************************************
