/************************************************************************//**
 * File: srprdint.cpp
 * Description: Auxiliary (obsolete?) SR calculation method from ~Arbitrary Transversely-Uniform Magnetic Field, in Near-Field observation conditions
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifdef __IGOR_PRO__
#include "srigintr.h"
#endif
#include "srprdint.h"

//*************************************************************************

//srTPartAutoRadInt::srTPartAutoRadInt(srTTrjDat* InTrjDatPtr, srTWfrSmp* InDistrInfoDatPtr, srLambXYZ* InObsCoorPtr, double In_sStart, double In_sFin, double In_sRelPrec, int InNpOnZeroLevel) 
srTPartAutoRadInt::srTPartAutoRadInt(srTTrjDat* InTrjDatPtr, srTWfrSmp* InDistrInfoDatPtr, srLambXYZ* InObsCoorPtr, double In_sStart, double In_sFin, double In_sRelPrec, long long InNpOnZeroLevel) 
{
	TrjDatPtr = InTrjDatPtr; DistrInfoDatPtr = InDistrInfoDatPtr; ObsCoorPtr = InObsCoorPtr;
	sIntegStart = In_sStart; sIntegFin = In_sFin; sIntegRelPrec = In_sRelPrec;
	NpOnZeroLevel = InNpOnZeroLevel;

	AuxDataPtr = 0;
	SomethingIsWrong = 0;

	//int TwoNpOnZeroLevel_mi_1 = 2*NpOnZeroLevel - 1;
	long long TwoNpOnZeroLevel_mi_1 = 2*NpOnZeroLevel - 1;
	double *BufStPtr, *BufPtr;
	if(TwoNpOnZeroLevel_mi_1*16 < 1000) BufStPtr = AuxDataArray;
	else
	{
		AuxDataPtr = new double[TwoNpOnZeroLevel_mi_1*16];
		//if(AuxDataPtr == 0) { Send.ErrorMessage("SR::Error900"); SomethingIsWrong = 1; return;}
		if(AuxDataPtr == 0) { SomethingIsWrong = 1; return;}
		BufStPtr = AuxDataPtr;
	}
	BufPtr = BufStPtr + TwoNpOnZeroLevel_mi_1*8;
	double *BufBtxArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufXArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufIntBtxE2ArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufBxArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufBtzArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufZArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufIntBtzE2ArrP = BufPtr; BufPtr += TwoNpOnZeroLevel_mi_1;
	double *BufBzArrP = BufPtr;

	TrjDatPtr->CompTotalTrjData(sIntegStart, sIntegFin, TwoNpOnZeroLevel_mi_1, BufBtxArrP, BufBtzArrP, BufXArrP, BufZArrP, BufIntBtxE2ArrP, BufIntBtzE2ArrP, BufBxArrP, BufBzArrP);

	//int NpOnZeroLevel_mi_1 = NpOnZeroLevel - 1;
	long long NpOnZeroLevel_mi_1 = NpOnZeroLevel - 1;
	double *TrTotDat = BufStPtr + TwoNpOnZeroLevel_mi_1*8;
	for(int k=0; k<8; k++) 
	{
		double *TrLev0 = BufStPtr + NpOnZeroLevel*k, *TrLev1 = BufStPtr + NpOnZeroLevel*8 + NpOnZeroLevel_mi_1*k;
		//for(int i=0; i<NpOnZeroLevel_mi_1; i++)
		for(long long i=0; i<NpOnZeroLevel_mi_1; i++)
		{
			*(TrLev0++) = *(TrTotDat++);
			*(TrLev1++) = *(TrTotDat++);
		}
		*(TrLev0++) = *(TrTotDat++);
	}
	BufPtr = BufStPtr;
	*BtxArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*XArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*IntBtxE2ArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*BxArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*BtzArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*ZArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*IntBtzE2ArrP = BufPtr; BufPtr += NpOnZeroLevel;
	*BzArrP = BufPtr; BufPtr += NpOnZeroLevel;

	*(BtxArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(XArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(IntBtxE2ArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(BxArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(BtzArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(ZArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(IntBtzE2ArrP+1) = BufPtr; BufPtr += NpOnZeroLevel_mi_1;
	*(BzArrP+1) = BufPtr;

	*AmOfPointsOnLevel = NpOnZeroLevel;
	*(AmOfPointsOnLevel+1) = NpOnZeroLevel_mi_1;
	for(int j=2; j<50; j++) *(AmOfPointsOnLevel+j) = 0;
	NumberOfLevelsFilled = 2;

	MaxFluxDensVal = CurrentAbsPrec = 0.;

	wfe = 7./15.; wf1 = 16./15.; wf2 = 14./15.; wd = 1./15.;
	PI = 3.141592653590;
	TwoPI = 2.*PI;
	ThreePIdTwo = 1.5*PI;

	FivePIdFour = 1.25*PI;

	HalfPI = 0.5* PI;
	One_dTwoPI = 1./TwoPI;
	PIm10e6 = PI*1.E+06;
	PIm10e6dEnCon = PIm10e6*0.80654658;
	TenEm6dTwoPI = 1.E-06/TwoPI;
	a2c = -0.5; a4c = 1./24.; a6c = -1./720.; a8c = 1./40320.; a10c = -1./3628800.; a12c = 1./479001600.;
	a3s = -1./6.; a5s = 1./120.; a7s = -1./5040.; a9s = 1./362880.; a11s = -1./39916800.; a13s = 1./6227020800.;
}

//*************************************************************************

int srTPartAutoRadInt::Integrate(double& OutIntXRe, double& OutIntXIm, double& OutIntZRe, double& OutIntZIm, double* InEdgeFunDer)
{
	double ActNormConst = (DistrInfoDatPtr->TreatLambdaAsEnergyIn_eV)? DistrInfoDatPtr->NormalizingConst*ObsCoorPtr->Lamb*0.80654658E-03 : DistrInfoDatPtr->NormalizingConst/ObsCoorPtr->Lamb;
	double PIm10e9_d_Lamb = (DistrInfoDatPtr->TreatLambdaAsEnergyIn_eV)? PIm10e6dEnCon*ObsCoorPtr->Lamb : PIm10e6*1000./ObsCoorPtr->Lamb;
	char NearField = (DistrInfoDatPtr->CoordOrAngPresentation == CoordPres);

	double EdgeFunDer[16];
	char EdgeFunDerNotComputed = 1;
	if(InEdgeFunDer != 0)
	{
		for(int k=0; k<16; k++) EdgeFunDer[k] = InEdgeFunDer[k];
		EdgeFunDerNotComputed = 0;
	}

	//int NpOnLevel = NpOnZeroLevel - 1;
	long long NpOnLevel = NpOnZeroLevel - 1;

	double sStep = (sIntegFin - sIntegStart)/NpOnLevel;
	double Ax, Az, Ph, CosPh, SinPh, PhPrev, PhInit;

	double xObs = ObsCoorPtr->x, yObs = ObsCoorPtr->y, zObs = ObsCoorPtr->z, GmEm2 = TrjDatPtr->EbmDat.GammaEm2;
	double One_d_ymis, xObs_mi_x, zObs_mi_z, Nx, Nz, Btx_mi_Nx, Btz_mi_Nz, dPhds, dAxds, dAzds;
	double AngPhConst, Two_xObs, Two_zObs;
	
	double Sum1XRe=0., Sum1XIm=0., Sum1ZRe=0., Sum1ZIm=0., Sum2XRe=0., Sum2XIm=0., Sum2ZRe=0., Sum2ZIm=0.;
	double wFxRe, wFxIm, wFzRe, wFzIm;
	double *pBtx = *BtxArrP, *pBtz = *BtzArrP, *pX = *XArrP, *pZ = *ZArrP, *pIntBtxE2 = *IntBtxE2ArrP, *pIntBtzE2 = *IntBtzE2ArrP;

	int result;

	if(NearField) 
	{
		One_d_ymis = 1./(yObs - sIntegStart); 
		xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
		Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
		Ph = PIm10e9_d_Lamb*(sIntegStart*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + xObs_mi_x*Nx + zObs_mi_z*Nz);
		Btx_mi_Nx = *(pBtx++) - Nx; Btz_mi_Nz = *(pBtz++) - Nz;
		Ax = Btx_mi_Nx*One_d_ymis; Az = Btz_mi_Nz*One_d_ymis;
	}
	else
	{
		AngPhConst = GmEm2 + xObs*xObs + zObs*zObs; Two_xObs = 2.*xObs; Two_zObs = 2.*zObs;
		Ph = PIm10e9_d_Lamb*(sIntegStart*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
	}
	PhInit = Ph;

	if(EdgeFunDerNotComputed)
	{
		if(NearField)
		{
			dPhds = PIm10e9_d_Lamb*(GmEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
			dAxds = (2.*Ax + (TrjDatPtr->BetaNormConst)*(**BzArrP))*One_d_ymis;
			dAzds = (2.*Az + (-TrjDatPtr->BetaNormConst)*(**BxArrP))*One_d_ymis;
		}
		else
		{
			dPhds = PIm10e9_d_Lamb*(GmEm2 + Ax*Ax + Az*Az);
			dAxds = (TrjDatPtr->BetaNormConst)*(**BzArrP);
			dAzds = (-TrjDatPtr->BetaNormConst)*(**BxArrP);
		}
		CosAndSin(Ph, CosPh, SinPh);
		wFxRe = Ax*CosPh; wFxIm = Ax*SinPh; wFzRe = Az*CosPh; wFzIm = Az*SinPh;

		double dPhdsSinPh = dPhds*SinPh, dPhdsCosPh = dPhds*CosPh;
		EdgeFunDer[4] = dAxds*CosPh - Ax*dPhdsSinPh;
		EdgeFunDer[5] = dAxds*SinPh + Ax*dPhdsCosPh;
		EdgeFunDer[6] = dAzds*CosPh - Az*dPhdsSinPh;
		EdgeFunDer[7] = dAzds*SinPh + Az*dPhdsCosPh;
	}
	else
	{
		wFxRe = *EdgeFunDer; wFxIm = EdgeFunDer[1]; wFzRe = EdgeFunDer[2]; wFzIm = EdgeFunDer[3];
	}

	//int AmOfLoops = (NpOnZeroLevel - 3) >> 1;
	long long AmOfLoops = (NpOnZeroLevel - 3) >> 1;
	double s = sIntegStart + sStep;

	//for(int i=0; i<AmOfLoops; i++)
	for(long long i=0; i<AmOfLoops; i++)
	{
		if(NearField)
		{
			One_d_ymis = 1./(yObs - s); 
			xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
			Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + xObs_mi_x*Nx + zObs_mi_z*Nz);
			Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
		}
		CosAndSin(Ph, CosPh, SinPh);
		Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

		if(NearField) 
		{
			One_d_ymis = 1./(yObs - s); 
			xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
			Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + xObs_mi_x*Nx + zObs_mi_z*Nz);
			Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
			Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
		}
		CosAndSin(Ph, CosPh, SinPh);
		Sum2XRe += Ax*CosPh; Sum2XIm += Ax*SinPh; Sum2ZRe += Az*CosPh; Sum2ZIm += Az*SinPh; s += sStep;
	}

	if(NearField) 
	{
		One_d_ymis = 1./(yObs - s); 
		xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
		Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
		Ph = PIm10e9_d_Lamb*(s*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + xObs_mi_x*Nx + zObs_mi_z*Nz);
		Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
	}
	else
	{
		Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
		Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
	}
	CosAndSin(Ph, CosPh, SinPh);
	Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

	if(EdgeFunDerNotComputed)
	{
		if(NearField)
		{
			One_d_ymis = 1./(yObs - s); 
			xObs_mi_x = xObs - *pX; zObs_mi_z = zObs - *pZ;
			Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
			Ph = PIm10e9_d_Lamb*(s*GmEm2 + *pIntBtxE2 + *pIntBtzE2 + xObs_mi_x*Nx + zObs_mi_z*Nz);
			Btx_mi_Nx = *pBtx - Nx; Btz_mi_Nz = *pBtz - Nz;
			Ax = Btx_mi_Nx*One_d_ymis; Az = Btz_mi_Nz*One_d_ymis;

			dPhds = PIm10e9_d_Lamb*(GmEm2 + Btx_mi_Nx*Btx_mi_Nx + Btz_mi_Nz*Btz_mi_Nz);
			dAxds = (2.*Ax + (TrjDatPtr->BetaNormConst)*((*BzArrP)[NpOnLevel]))*One_d_ymis;
			dAzds = (2.*Az + (-TrjDatPtr->BetaNormConst)*((*BxArrP)[NpOnLevel]))*One_d_ymis;
		}
		else
		{
			Ph = PIm10e9_d_Lamb*(s*AngPhConst + *pIntBtxE2 + *pIntBtzE2 - (Two_xObs*(*pX) + Two_zObs*(*pZ)));
			Ax = *pBtx - xObs; Az = *pBtz - zObs;

			dPhds = PIm10e9_d_Lamb*(GmEm2 + Ax*Ax + Az*Az);
			dAxds = (TrjDatPtr->BetaNormConst)*((*BzArrP)[NpOnLevel]);
			dAzds = (-TrjDatPtr->BetaNormConst)*((*BxArrP)[NpOnLevel]);
		}
		CosAndSin(Ph, CosPh, SinPh);
		wFxRe += Ax*CosPh; wFxIm += Ax*SinPh; wFzRe += Az*CosPh; wFzIm += Az*SinPh;

		double dPhdsSinPh = dPhds*SinPh, dPhdsCosPh = dPhds*CosPh;
		EdgeFunDer[12] = dAxds*CosPh - Ax*dPhdsSinPh;
		EdgeFunDer[13] = dAxds*SinPh + Ax*dPhdsCosPh;
		EdgeFunDer[14] = dAzds*CosPh - Az*dPhdsSinPh;
		EdgeFunDer[15] = dAzds*SinPh + Az*dPhdsCosPh;
	}
	else
	{
		wFxRe += EdgeFunDer[8]; wFxIm += EdgeFunDer[9]; wFzRe += EdgeFunDer[10]; wFzIm += EdgeFunDer[11];
	}
	wFxRe *= wfe; wFxIm *= wfe; wFzRe *= wfe; wFzIm *= wfe;

	double wDifDerXRe = wd*(EdgeFunDer[4] - EdgeFunDer[12]), wDifDerXIm = wd*(EdgeFunDer[5] - EdgeFunDer[13]);
	double wDifDerZRe = wd*(EdgeFunDer[6] - EdgeFunDer[14]), wDifDerZIm = wd*(EdgeFunDer[7] - EdgeFunDer[15]);

	double ActNormConstStep = ActNormConst*sStep;
	double IntXRe = ActNormConstStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + sStep*wDifDerXRe);
	double IntXIm = ActNormConstStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + sStep*wDifDerXIm);
	double IntZRe = ActNormConstStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + sStep*wDifDerZRe);
	double IntZIm = ActNormConstStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + sStep*wDifDerZIm);
	double SqNorm = IntXRe*IntXRe + IntXIm*IntXIm + IntZRe*IntZRe + IntZIm*IntZIm;

	int LevelNo = 0;
	char NotFinishedYet = 1;
	while(NotFinishedYet)
	{
		Sum2XRe += Sum1XRe; Sum2XIm += Sum1XIm; Sum2ZRe += Sum1ZRe; Sum2ZIm += Sum1ZIm; 
		Sum1XRe = Sum1XIm = Sum1ZRe = Sum1ZIm = 0.;
		char ThisMayBeTheLastLoop = 1;
		PhPrev = PhInit;
		LevelNo++;

		double HalfStep = 0.5*sStep;
		s = sIntegStart + HalfStep;
		if(NumberOfLevelsFilled <= LevelNo) if(result = FillNextLevel(LevelNo, s, sIntegFin - HalfStep, NpOnLevel)) return result;
		pBtx = BtxArrP[LevelNo], pBtz = BtzArrP[LevelNo], pX = XArrP[LevelNo], pZ = ZArrP[LevelNo], pIntBtxE2 = IntBtxE2ArrP[LevelNo], pIntBtzE2 = IntBtzE2ArrP[LevelNo];
		//for(int i=0; i<NpOnLevel; i++)
		for(long long i=0; i<NpOnLevel; i++)
		{
			if(NearField)
			{
				One_d_ymis = 1./(yObs - s);
				xObs_mi_x = xObs - *(pX++); zObs_mi_z = zObs - *(pZ++);
				Nx = xObs_mi_x*One_d_ymis, Nz = zObs_mi_z*One_d_ymis;
				Ph = PIm10e9_d_Lamb*(s*GmEm2 + *(pIntBtxE2++) + *(pIntBtzE2++) + xObs_mi_x*Nx + zObs_mi_z*Nz);
				Ax = (*(pBtx++) - Nx)*One_d_ymis; Az = (*(pBtz++) - Nz)*One_d_ymis;
			}
			else
			{
				Ph = PIm10e9_d_Lamb*(s*AngPhConst + *(pIntBtxE2++) + *(pIntBtzE2++) - (Two_xObs*(*(pX++)) + Two_zObs*(*(pZ++))));
				Ax = *(pBtx++) - xObs; Az = *(pBtz++) - zObs;
			}
			CosAndSin(Ph, CosPh, SinPh);
			Sum1XRe += Ax*CosPh; Sum1XIm += Ax*SinPh; Sum1ZRe += Az*CosPh; Sum1ZIm += Az*SinPh; s += sStep;

			if(Ph - PhPrev > FivePIdFour) ThisMayBeTheLastLoop = 0;
			PhPrev = Ph;
		}
		double ActNormConstHalfStep = ActNormConst*HalfStep;
		double LocIntXRe = ActNormConstHalfStep*(wFxRe + wf1*Sum1XRe + wf2*Sum2XRe + HalfStep*wDifDerXRe);
		double LocIntXIm = ActNormConstHalfStep*(wFxIm + wf1*Sum1XIm + wf2*Sum2XIm + HalfStep*wDifDerXIm);
		double LocIntZRe = ActNormConstHalfStep*(wFzRe + wf1*Sum1ZRe + wf2*Sum2ZRe + HalfStep*wDifDerZRe);
		double LocIntZIm = ActNormConstHalfStep*(wFzIm + wf1*Sum1ZIm + wf2*Sum2ZIm + HalfStep*wDifDerZIm);
		double LocSqNorm = LocIntXRe*LocIntXRe + LocIntXIm*LocIntXIm + LocIntZRe*LocIntZRe + LocIntZIm*LocIntZIm;
		if(ThisMayBeTheLastLoop)
		{
			double TestVal = ::fabs(LocSqNorm - SqNorm);
			if(MaxFluxDensVal > 0) NotFinishedYet = (TestVal > CurrentAbsPrec);
			else NotFinishedYet = (TestVal/LocSqNorm > sIntegRelPrec);
		}
		IntXRe = LocIntXRe; IntXIm = LocIntXIm; IntZRe = LocIntZRe; IntZIm = LocIntZIm; SqNorm = LocSqNorm;
		sStep = HalfStep; NpOnLevel *= 2;
	}
	OutIntXRe += IntXRe; OutIntXIm += IntXIm; OutIntZRe += IntZRe; OutIntZIm += IntZIm;

	if(MaxFluxDensVal < SqNorm) 
	{
		MaxFluxDensVal = SqNorm; CurrentAbsPrec = sIntegRelPrec*MaxFluxDensVal;
	}
	return 1;
}

//*************************************************************************

