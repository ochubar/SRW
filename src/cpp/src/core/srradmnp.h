/************************************************************************//**
 * File: srradmnp.h
 * Description: Various "manipulations" with Radiation data (e.g. "extraction" of Intensity from Electric Field, etc.) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRRADMNP_H
#define __SRRADMNP_H

#include "sroptelm.h"
#include "srstraux.h"
#include "srerror.h"
#include "gminterp.h"

#include <complex>

//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

struct srTAuxInt2DIntegOverAzim {
	double x0, y0, r;
	double xMin, xStep, yMin, yStep;
	//long nx, ny;
	long long nx, ny; //OC26042019
	float *pfI2D;
	double *pdI2D;
	char ordInterp;
	double *arRectToSkipPar; //OC09032019
	int nRectToSkip;
};

//*************************************************************************

class srTRadGenManip {
// Various manipulations with computed Radiation
	bool EhOK, EvOK; //OC111111

public:
	//srTSRWRadStructAccessData RadAccessData;
	CHGenObj hRadAccessData;

	//srTRadGenManip(srTSRWRadStructAccessData& InRadAccessData)
	//{
	//	RadAccessData = InRadAccessData;
	//}
	srTRadGenManip(CHGenObj& In_hRadAccessData)
	{
		//RadAccessData = InRadAccessData;
		hRadAccessData = In_hRadAccessData;

		EhOK = EvOK = false; //OC111111
		srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
		EhOK = (RadAccessData.pBaseRadX != 0);
		EvOK = (RadAccessData.pBaseRadZ != 0);
	}
	srTRadGenManip() 
	{
		EhOK = EvOK = false; //OC111111
	}

	void ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData, double* pMeth=0, srTTrjDat* pTrjDat=0) //OC23022020
	//void ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData, double* pMeth=0) //OC16122019
	//void ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData, int* pMeth=0) //OC13122019
	//void ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData);
	{
		if(pData == 0) throw INCORRECT_PARAMS_WFR_COMPON_EXTRACT;

		srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));
		srTGenOptElem GenOptElem;
		srTRadExtract RadExtract(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData, pMeth); //OC13122019
		//srTRadExtract RadExtract(PolarizCompon, Int_or_Phase, SectID, TransvPres, e, x, z, pData);

		//DEBUG
		//srTIgorSend::WarningMessage("ExtractRadiation #1");
		//END DEBUG

		int res;
		if(TransvPres != RadAccessData.Pres)
			if(res = GenOptElem.SetRadRepres(&RadAccessData, char(TransvPres))) throw res;

		if(RadExtract.Int_or_Phase == 1)
		{//1- Multi-Elec Intensity
			if(res = ComputeConvolutedIntensity(RadExtract)) throw res;
		}
		else if(RadExtract.Int_or_Phase == 4) //OC
		{//4- Sigle-Elec. Flux
			if(res = ExtractFluxFromWfr(RadExtract, 's')) throw res;
		}
		else if(RadExtract.Int_or_Phase == 5) //OC
		{//5- Multi-Elec. Flux
			if(res = ExtractFluxFromWfr(RadExtract, 'm')) throw res;
		}
		else if(RadExtract.Int_or_Phase == 8) //OC06092018
		{
			if(res = ExtractSingleElecMutualIntensity(RadExtract)) throw res;
		}
		//else if(RadExtract.Int_or_Phase == 9) //OC23022020 (under development - OC23032020)
		//{
		//	if(res = ExtractMultiElecMutualIntensity(RadExtract, pTrjDat)) throw res;
		//}
		else
		{
			if(res = ExtractSingleElecIntensity(RadExtract)) throw res;
		}

		//OCTEST17082019
		//if((RadExtract.Int_or_Phase == 2) && (RadExtract.PlotType == 3))
		//{
		//	srTWaveAccessData ExtrWaveData;
		//	if(res = SetupExtractedWaveData(RadExtract, ExtrWaveData)) throw res;
		//	if(res = TryToMakePhaseContinuous(ExtrWaveData)) throw res;
		//}
				//DEBUG
				//srTIgorSend::WarningMessage("ExtractRadiation: ready out");
				//END DEBUG
		//if(res = SetupExtractedWaveData(RadExtract, ExtractedWaveData)) return res;
	}
	
	//OC13122019 (commented-out, to reduce number of function calls)
	//void ExtractRadiationSRWL(char polar, char intType, char depType, char transvPres, double e, double x, double z, char* pData) //OC19082018
	//{
	//	//Re-defining intType
	//	//from SRWL convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Flux; 3- Multi-Elec. Flux; 4- Single-Elec. Rad. Phase; 5- Re(E); 6- Im(E); 7- Time or Photon Energy Integrated Intensity
	//	//to old SRW convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Rad. Phase; 3- Re(E); 4- Single-Elec. Flux; 5- Multi-Elec. Flux; 6- Im(E); 7- Time or Photon Energy Integrated Intensity
	//	if(intType == 2) intType = 4;
	//	else if(intType == 3) intType = 5;
	//	else if(intType == 4) intType = 2;
	//	else if(intType == 5) intType = 3;
	//	ExtractRadiation((int)polar, (int)intType, (int)depType, (int)transvPres, e, x, z, pData);
	//}

	int ExtractRadiation(srTRadExtract& RadExtract, srTWaveAccessData& ExtractedWaveData)
	{
		int result;
        srTSRWRadStructAccessData& RadAccessData = *((srTSRWRadStructAccessData*)(hRadAccessData.ptr()));

		srTGenOptElem GenOptElem;
		if(RadExtract.TransvPres != RadAccessData.Pres)
			if(result = GenOptElem.SetRadRepres(&RadAccessData, char(RadExtract.TransvPres))) return result;
		if(RadExtract.Int_or_Phase == 1)
		{
			if(result = ComputeConvolutedIntensity(RadExtract)) return result;
		}
		else
		{
			if(result = ExtractSingleElecIntensity(RadExtract)) return result;
		}
		if(result = SetupExtractedWaveData(RadExtract, ExtractedWaveData)) return result;

		//For tests only
			//if((RadExtract.Int_or_Phase == 2) && (RadExtract.PlotType == 3)) TryToMakePhaseContinuous(ExtractedWaveData);
		return 0;
	}

	int ExtractSingleElecIntensity(srTRadExtract& RadExtract)
	{
		if(RadExtract.PlotType == 0) return ExtractSingleElecIntensity1DvsE(RadExtract);
		else if(RadExtract.PlotType == 1) return ExtractSingleElecIntensity1DvsX(RadExtract);
		else if(RadExtract.PlotType == 2) return ExtractSingleElecIntensity1DvsZ(RadExtract);
		else if(RadExtract.PlotType == 3) return ExtractSingleElecIntensity2DvsXZ(RadExtract);
		else if(RadExtract.PlotType == 4) return ExtractSingleElecIntensity2DvsEX(RadExtract);
		else if(RadExtract.PlotType == 5) return ExtractSingleElecIntensity2DvsEZ(RadExtract);
		else return ExtractSingleElecIntensity3D(RadExtract);
	}

	int ExtractSingleElecMutualIntensity(srTRadExtract& RadExtract) 
	{//OC06092018
		int PolCom = RadExtract.PolarizCompon;
		int Int_or_ReE = RadExtract.Int_or_Phase;
		//if((PolCom == 6) || (Int_or_ReE != 8)) return CAN_NOT_EXTRACT_MUT_INT;
		if(Int_or_ReE != 8) return CAN_NOT_EXTRACT_MUT_INT;

		//if(RadExtract.PlotType == 0) return ExtractSingleElecMutualIntensityVsE(RadExtract);
		//else if(RadExtract.PlotType == 1) return ExtractSingleElecMutualIntensityVsX(RadExtract);
		if(RadExtract.PlotType == 1) return ExtractSingleElecMutualIntensityVsX(RadExtract);
		else if(RadExtract.PlotType == 2) return ExtractSingleElecMutualIntensityVsZ(RadExtract);
		else if(RadExtract.PlotType == 3) return ExtractSingleElecMutualIntensityVsXZ(RadExtract);
		//else if(RadExtract.PlotType == 4) return ExtractSingleElecMutualIntensityVsEX(RadExtract);
		//else if(RadExtract.PlotType == 5) return ExtractSingleElecMutualIntensityVsEZ(RadExtract);
		//else return ExtractSingleElecMutualIntensityEXZ(RadExtract);
		else return CAN_NOT_EXTRACT_MUT_INT;
	}

/**
	int ExtractMultiElecMutualIntensity(srTRadExtract& RadExtract, srTTrjDat* pTrjDat=0)
	{//OC23022020 (under development - OC23032020)
		int PolCom = RadExtract.PolarizCompon;
		int Int_or_ReE = RadExtract.Int_or_Phase;
		if(Int_or_ReE != 9) return CAN_NOT_EXTRACT_MUT_INT;

		//if(RadExtract.PlotType == 0) return ExtractMultiElecMutualIntensityVsE(RadExtract);
		//else if(RadExtract.PlotType == 1) return ExtractMultiElecMutualIntensityVsX(RadExtract);
		//if(RadExtract.PlotType == 1) return ExtractMultiElecMutualIntensityVsX(RadExtract);
		//else if(RadExtract.PlotType == 2) return ExtractMultiElecMutualIntensityVsZ(RadExtract);
		if(RadExtract.PlotType == 3) return ComputeMultiElecMutualIntensityVsXZ(RadExtract, pTrjDat);
		//else if(RadExtract.PlotType == 4) return ExtractMultiElecMutualIntensityVsEX(RadExtract);
		//else if(RadExtract.PlotType == 5) return ExtractMultiElecMutualIntensityVsEZ(RadExtract);
		//else return ExtractMultiElecMutualIntensityEXZ(RadExtract);
		else return CAN_NOT_EXTRACT_MUT_INT;
	}
**/
	int TryToMakePhaseContinuous(srTWaveAccessData& WaveData);
	void TryToMakePhaseContinuous1D(double* pOutPhase, long long Np, long long i0, float Phi0); //OC26042019
	//void TryToMakePhaseContinuous1D(double* pOutPhase, long Np, long i0, float Phi0);

	int ExtractSingleElecIntensity1DvsE(srTRadExtract&);
	int ExtractSingleElecIntensity1DvsX(srTRadExtract&);
	int ExtractSingleElecIntensity1DvsZ(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsXZ(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsEX(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsEZ(srTRadExtract&);
	int ExtractSingleElecIntensity3D(srTRadExtract&);

	int ExtractSingleElecMutualIntensityVsX(srTRadExtract&); //OC06092018
	int ExtractSingleElecMutualIntensityVsZ(srTRadExtract&);
	int ExtractSingleElecMutualIntensityVsXZ(srTRadExtract&);

	int ComputeMultiElecMutualIntensityVsXZ(srTRadExtract&, srTTrjDat* pTrjDat=0); //23022020

	int SetupExtractedWaveData(srTRadExtract&, srTWaveAccessData&);
	void SetupIntCoord(char, double, long long&, long long&, double&); //OC26042019
	//void SetupIntCoord(char, double, long&, long&, double&);

	int ExtractFluxFromWfr(srTRadExtract& RadExtract, char s_or_m)
	{
		if(RadExtract.PlotType != 0) 
		{
			CErrWarn::AddWarningMessage(&gVectWarnNos, CAN_ONLY_EXTRACT_FLUX_VS_PHOTON_ENERGY);
			if(RadExtract.pExtractedData != 0) *(RadExtract.pExtractedData) = 0.;
		}
		else
		{
			if((s_or_m == 's') || (s_or_m == 'S')) return ExtractSingleElecFlux1DvsE(RadExtract);
			else return ExtractMultiElecFlux1DvsE(RadExtract);
		}
		return 0;
	}

	int ExtractSingleElecFlux1DvsE(srTRadExtract&);
	int ExtractMultiElecFlux1DvsE(srTRadExtract&);

	int ComputeConvolutedIntensity(srTRadExtract&);
	//int ConvoluteWithElecBeamOverTransvCoord(float*, long long, long long);
	int ConvoluteWithElecBeamOverTransvCoord(float*, long, long);
	void PutConstPhotEnergySliceInExtractPlace(long long, long long, long long, srTRadExtract&, srTRadExtract&); //OC26042019
	//void PutConstPhotEnergySliceInExtractPlace(long, long, long, srTRadExtract&, srTRadExtract&);
	void PropagateElecBeamMoments(srTElecBeamMoments&);
	static void PropagateElecBeamMoments(srTElecBeamMoments& ElecBeamMom, double* p4x4PropMatr, double* p4Vect);

	void PadImZerosToRealData(float*, long long, long long); //OC26042019
	//void PadImZerosToRealData(float*, long, long);
	//void ShiftData(float* pStart, long LenData, long ShiftLen)
	void ShiftData(float* pStart, long long LenData, long long ShiftLen)
	{
		//long LenData_mi_1 = LenData - 1;
		long long LenData_mi_1 = LenData - 1;
		float* pOrig = pStart + LenData_mi_1;
		float* pFin = pOrig + ShiftLen;
		for(long i=0; i<LenData; i++) *(pFin--) = *(pOrig--);
	}
	//void SetDataToZero(float* pStart, long LenData)
	void SetDataToZero(float* pStart, long long LenData)
	{
		float* p = pStart;
		//for(long i=0; i<LenData; i++) *(p++) = 0.;
		for(long long i=0; i<LenData; i++) *(p++) = 0.;
	}

	void SetupPolarizVect(complex<float>* PolVect, int PolCom)
	{
		for(int i=0; i<2; i++) PolVect[i] = 0.;
		const float Norm = (float)(1./sqrt(2.));
		complex<float> OneRe(1., 0.), OneIm(0., 1.), OneReN(Norm, 0.), OneImN(0., Norm);

		if(PolCom==0) *PolVect = OneRe;
		else if(PolCom==1) *(PolVect+1) = OneRe;
		else if(PolCom==2) { *PolVect = OneReN; *(PolVect+1) = OneReN;}
		else if(PolCom==3) { *PolVect = OneReN; *(PolVect+1) = -OneReN;}
		else if(PolCom==4) { *PolVect = OneReN; *(PolVect+1) = OneImN;}
		else if(PolCom==5) { *PolVect = OneReN; *(PolVect+1) = -OneImN;}
	}
	float IntensityComponentSimpleInterpol(float* pEx_St, float* pEx_Fi, float* pEz_St, float* pEz_Fi, double InvStepRelArg, int PolCom, int Int_or_ReE)
	{
		float I_St = IntensityComponent(pEx_St, pEz_St, PolCom, Int_or_ReE);
		if(Int_or_ReE == 2) return I_St;

		float I_Fi = IntensityComponent(pEx_Fi, pEz_Fi, PolCom, Int_or_ReE);
		return (float)((I_Fi - I_St)*InvStepRelArg + I_St);
	}
	float IntensityComponentSimpleInterpol2D(float** ExPtrs, float** EzPtrs, double Arg1, double Arg2, int PolCom, int Int_or_ReE)
	{
		float I00 = IntensityComponent(*ExPtrs, *EzPtrs, PolCom, Int_or_ReE);
		if(Int_or_ReE == 2) return I00;

		float I10 = IntensityComponent(*(ExPtrs + 1), *(EzPtrs + 1), PolCom, Int_or_ReE);
		float I01 = IntensityComponent(*(ExPtrs + 2), *(EzPtrs + 2), PolCom, Int_or_ReE);
		float I11 = IntensityComponent(*(ExPtrs + 3), *(EzPtrs + 3), PolCom, Int_or_ReE);
		double Arg1Arg2 = Arg1*Arg2;
		return (float)((I00 - I01 - I10 + I11)*Arg1Arg2 + (I10 - I00)*Arg1 + (I01 - I00)*Arg2 + I00);
	}
	float IntensityComponent(float* pEx, float* pEz, int PolCom, int Int_or_ReE)
	{
		//float ExRe = *pEx, ExIm = *(pEx + 1), EzRe = *pEz, EzIm = *(pEz + 1);
		float ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.; //OC111111
		if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1);}
		if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1);}

		switch(PolCom)
		{
			//case 0: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? ExRe*ExRe + ExIm*ExIm : ((Int_or_ReE == 2)? FormalPhase(ExRe, ExIm) : ExRe))); // Lin. Hor.
			case 0: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? ExRe*ExRe + ExIm*ExIm : ((Int_or_ReE == 2)? FormalPhase(ExRe, ExIm) : ((Int_or_ReE == 3)? ExRe : ExIm)))); // Lin. Hor. //OC031208

			//case 1: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? EzRe*EzRe + EzIm*EzIm : ((Int_or_ReE == 2)? FormalPhase(EzRe, EzIm) : EzRe))); // Lin. Vert. 
			case 1: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? EzRe*EzRe + EzIm*EzIm : ((Int_or_ReE == 2)? FormalPhase(EzRe, EzIm) : ((Int_or_ReE == 3)? EzRe : EzIm)))); // Lin. Vert. //OC031208

			case 2: // Linear 45 deg.
			{
				float ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
				//return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_p_EzRe*ExRe_p_EzRe + ExIm_p_EzIm*ExIm_p_EzIm) : ((Int_or_ReE == 2)? FormalPhase(ExRe_p_EzRe, ExIm_p_EzIm) : 0.70710678*ExRe_p_EzRe)));
				return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_p_EzRe*ExRe_p_EzRe + ExIm_p_EzIm*ExIm_p_EzIm) : ((Int_or_ReE == 2)? FormalPhase(ExRe_p_EzRe, ExIm_p_EzIm) : 0.70710678*((Int_or_ReE == 3)? ExRe_p_EzRe : ExIm_p_EzIm)))); //OC031208
			}
			case 3: // Linear 135 deg.
			{
				float ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
				//return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzRe*ExRe_mi_EzRe + ExIm_mi_EzIm*ExIm_mi_EzIm) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzRe, ExIm_mi_EzIm) : 0.70710678*ExRe_mi_EzRe)));
				return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzRe*ExRe_mi_EzRe + ExIm_mi_EzIm*ExIm_mi_EzIm) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzRe, ExIm_mi_EzIm) : 0.70710678*((Int_or_ReE == 3)? ExRe_mi_EzRe : ExIm_mi_EzIm)))); //OC031208
			}
			case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 4: // Circ. Right
			{
				float ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
				//return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzIm*ExRe_mi_EzIm + ExIm_p_EzRe*ExIm_p_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzIm, ExIm_p_EzRe) : 0.70710678*ExRe_mi_EzIm)));
				return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzIm*ExRe_mi_EzIm + ExIm_p_EzRe*ExIm_p_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzIm, ExIm_p_EzRe) : 0.70710678*((Int_or_ReE == 3)? ExRe_mi_EzIm : ExIm_p_EzRe)))); //OC031208
			}
			case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 5: // Circ. Left
			{
				float ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
				//return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_p_EzIm*ExRe_p_EzIm + ExIm_mi_EzRe*ExIm_mi_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_p_EzIm, ExIm_mi_EzRe) : 0.70710678*ExRe_p_EzIm)));
				return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_p_EzIm*ExRe_p_EzIm + ExIm_mi_EzRe*ExIm_mi_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_p_EzIm, ExIm_mi_EzRe) : 0.70710678*((Int_or_ReE == 3)? ExRe_p_EzIm : ExIm_mi_EzRe)))); //OC031208
			}
			case -1: // s0
			{
				return (float)(ExRe*ExRe + ExIm*ExIm + EzRe*EzRe + EzIm*EzIm);
			}
			case -2: // s1
			{
				return (float)(ExRe*ExRe + ExIm*ExIm - (EzRe*EzRe + EzIm*EzIm));
			}
			case -3: // s2
			{
				return (float)(2.*(ExRe*EzRe + ExIm*EzIm)); //OC07092019 (in line with literature, for righ-hand frame (x,z,s), i.e. Ez meaning Ey in the frame (x,y,z))
				//return (float)(-2.*(ExRe*EzRe + ExIm*EzIm));
			}
			case -4: // s3
			{
				return (float)(2.*(ExRe*EzIm - ExIm*EzRe)); //OC07092019 (in line with literature, for righ-hand frame (x,z,s), i.e. Ez meaning Ey in the frame (x,y,z))
				//return (float)(2.*(-ExRe*EzIm + ExIm*EzRe));
			}
			//default: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? ExRe*ExRe + ExIm*ExIm + EzRe*EzRe + EzIm*EzIm : ((Int_or_ReE == 2)? FormalPhase(ExRe, ExIm) : ExRe)));
			default: return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? ExRe*ExRe + ExIm*ExIm + EzRe*EzRe + EzIm*EzIm : ((Int_or_ReE == 2)? FormalPhase(ExRe, ExIm) : ((Int_or_ReE == 3)? ExRe : ExIm)))); //OC031208
		}
		//return (float)(ExRe*ExRe + ExIm*ExIm + EzRe*EzRe + EzIm*EzIm);
	}
	double FormalPhase(float Re, float Im)
	{
		const double HalhPi = 1.5707963267949;
		const double Pi = 3.1415926535898;
		if(Re != 0.) 
		{
			if(Im <= 0.)
			{
				if(Re < 0.) return atan(double(Im/Re)) - Pi;
				else return atan(double(Im/Re));
			}
			else
			{
				if(Re < 0.) return atan(double(Im/Re)) + Pi;
				else return atan(double(Im/Re));
			}
		}
		else
		{
			if(Im == 0.) return  0.;
			else return (Im > 0.)? HalhPi : -HalhPi;
		}
	}
	
	int MutualIntensityComponentSimpleInterpol(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double InvStepRelArg, int PolCom, double iter, float* pResMI) //OC16122019
	//int MutualIntensityComponentSimpleInterpol(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double InvStepRelArg, int PolCom, int iter, float* pResMI) //OC14122019
	//int MutualIntensityComponentSimpleInterpol(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double InvStepRelArg, int PolCom, float* pResMI)
	{//OC12092018
		float MI0[2], MI1[2];
		int res = 0;
		if(res = MutualIntensityComponent(*ExPtrs, *ExPtrsT, *EzPtrs, *EzPtrsT, PolCom, iter, MI0)) return res; //OC14122019
		if(res = MutualIntensityComponent(*(ExPtrs + 1), *(ExPtrsT + 1), *(EzPtrs + 1), *(EzPtrsT + 1), PolCom, iter, MI1)) return res;
		//if(res = MutualIntensityComponent(*ExPtrs, *ExPtrsT, *EzPtrs, *EzPtrsT, PolCom, MI0)) return res;
		//if(res = MutualIntensityComponent(*(ExPtrs + 1), *(ExPtrsT + 1), *(EzPtrs + 1), *(EzPtrsT + 1), PolCom, MI1)) return res;
		double MI0Re = *MI0, MI0Im = *(MI0 + 1);
		double MI1Re = *MI1, MI1Im = *(MI1 + 1);
		*pResMI = (float)((MI1Re - MI0Re)*InvStepRelArg + MI0Re);
		*(pResMI + 1) = (float)((MI1Im - MI0Im)*InvStepRelArg + MI0Im);
		return 0;
	}

	int MutualIntensityComponentSimpleInterpol2D(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double Arg1, double Arg2, int PolCom, double iter, float* pResMI) //OC16122019
	//int MutualIntensityComponentSimpleInterpol2D(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double Arg1, double Arg2, int PolCom, int iter, float* pResMI) //OC14122019
	//int MutualIntensityComponentSimpleInterpol2D(float** ExPtrs, float** ExPtrsT, float** EzPtrs, float** EzPtrsT, double Arg1, double Arg2, int PolCom, float* pResMI)
	{//OC09092018
		float MI00[2], MI10[2], MI01[2], MI11[2];
		int res = 0;
		if(res = MutualIntensityComponent(*ExPtrs, *ExPtrsT, *EzPtrs, *EzPtrsT, PolCom, iter, MI00)) return res; //OC14122019
		if(res = MutualIntensityComponent(*(ExPtrs + 1), *(ExPtrsT + 1), *(EzPtrs + 1), *(EzPtrsT + 1), PolCom, iter, MI10)) return res;
		if(res = MutualIntensityComponent(*(ExPtrs + 2), *(ExPtrsT + 2), *(EzPtrs + 2), *(EzPtrsT + 2), PolCom, iter, MI01)) return res;
		if(res = MutualIntensityComponent(*(ExPtrs + 3), *(ExPtrsT + 3), *(EzPtrs + 3), *(EzPtrsT + 3), PolCom, iter, MI11)) return res;
		//if(res = MutualIntensityComponent(*ExPtrs, *ExPtrsT, *EzPtrs, *EzPtrsT, PolCom, MI00)) return res;
		//if(res = MutualIntensityComponent(*(ExPtrs + 1), *(ExPtrsT + 1), *(EzPtrs + 1), *(EzPtrsT + 1), PolCom, MI10)) return res;
		//if(res = MutualIntensityComponent(*(ExPtrs + 2), *(ExPtrsT + 2), *(EzPtrs + 2), *(EzPtrsT + 2), PolCom, MI01)) return res;
		//if(res = MutualIntensityComponent(*(ExPtrs + 3), *(ExPtrsT + 3), *(EzPtrs + 3), *(EzPtrsT + 3), PolCom, MI11)) return res;
		double Arg1Arg2 = Arg1*Arg2;
		double MI00Re = *MI00, MI00Im = *(MI00 + 1);
		double MI10Re = *MI10, MI10Im = *(MI10 + 1);
		double MI01Re = *MI01, MI01Im = *(MI01 + 1);
		double MI11Re = *MI11, MI11Im = *(MI11 + 1);
		*pResMI = (float)((MI00Re - MI01Re - MI10Re + MI11Re)*Arg1Arg2 + (MI10Re - MI00Re)*Arg1 + (MI01Re - MI00Re)*Arg2 + MI00Re);
		*(pResMI + 1) = (float)((MI00Im - MI01Im - MI10Im + MI11Im)*Arg1Arg2 + (MI10Im - MI00Im)*Arg1 + (MI01Im - MI00Im)*Arg2 + MI00Im);
		return 0;
	}

	int MutualIntensityComponent(float* pEx, float* pExT, float* pEz, float* pEzT, int PolCom, double iter, float* pResMI) //OC16122019
	//int MutualIntensityComponent(float* pEx, float* pExT, float* pEz, float* pEzT, int PolCom, long long iter, float* pResMI) //OC13122019
	//int MutualIntensityComponent(float* pEx, float* pExT, float* pEz, float* pEzT, int PolCom, float* pResMI)
	{//OC09092018
	 //NOTE: This is based on M.I. definition as: E(x)E*(x'), which differs from existing definition in literature: E*(x)E(x')
	 //The two definitions are related by complex conjugation: E*(x)E(x') = (E(x)E*(x'))*
		double ExRe = 0., ExIm = 0., EzRe = 0., EzIm = 0.;
		double ExReT = 0., ExImT = 0., EzReT = 0., EzImT = 0.;
		if(EhOK) { ExRe = *pEx; ExIm = *(pEx + 1); ExReT = *pExT; ExImT = *(pExT + 1);}
		if(EvOK) { EzRe = *pEz; EzIm = *(pEz + 1); EzReT = *pEzT; EzImT = *(pEzT + 1);}
		double ReMI = 0., ImMI = 0.;

		switch(PolCom)
		{
			case 0: // Lin. Hor.
			{
				ReMI = ExRe*ExReT + ExIm*ExImT;
				ImMI = ExIm*ExReT - ExRe*ExImT;
				break;
			}
			case 1: // Lin. Vert.
			{
				ReMI = EzRe*EzReT + EzIm*EzImT;
				ImMI = EzIm*EzReT - EzRe*EzImT;
				break;
			}
			case 2: // Linear 45 deg.
			{
				double ExRe_p_EzRe = ExRe + EzRe, ExIm_p_EzIm = ExIm + EzIm;
				double ExRe_p_EzReT = ExReT + EzReT, ExIm_p_EzImT = ExImT + EzImT;
				ReMI = 0.5*(ExRe_p_EzRe*ExRe_p_EzReT + ExIm_p_EzIm*ExIm_p_EzImT);
				ImMI = 0.5*(ExIm_p_EzIm*ExRe_p_EzReT - ExRe_p_EzRe*ExIm_p_EzImT);
				break;
			}
			case 3: // Linear 135 deg.
			{
				double ExRe_mi_EzRe = ExRe - EzRe, ExIm_mi_EzIm = ExIm - EzIm;
				double ExRe_mi_EzReT = ExReT - EzReT, ExIm_mi_EzImT = ExImT - EzImT;
				ReMI = 0.5*(ExRe_mi_EzRe*ExRe_mi_EzReT + ExIm_mi_EzIm*ExIm_mi_EzImT);
				ImMI = 0.5*(ExIm_mi_EzIm*ExRe_mi_EzReT - ExRe_mi_EzRe*ExIm_mi_EzImT);
				break;
			}
			case 5: // Circ. Left //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 4: // Circ. Right
			{
				double ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
				double ExRe_mi_EzImT = ExReT - EzImT, ExIm_p_EzReT = ExImT + EzReT;
				ReMI = 0.5*(ExRe_mi_EzIm*ExRe_mi_EzImT + ExIm_p_EzRe*ExIm_p_EzReT);
				ImMI = 0.5*(ExIm_p_EzRe*ExRe_mi_EzImT - ExRe_mi_EzIm*ExIm_p_EzReT);
				break;
			}
			case 4: // Circ. Right //OC08092019: corrected to be in compliance with definitions for right-hand frame (x,z,s) and with corresponding definition and calculation of Stokes params
			//case 5: // Circ. Left
			{
				double ExRe_p_EzIm = ExRe + EzIm, ExIm_mi_EzRe = ExIm - EzRe;
				double ExRe_p_EzImT = ExReT + EzImT, ExIm_mi_EzReT = ExImT - EzReT;
				ReMI = 0.5*(ExRe_p_EzIm*ExRe_p_EzImT + ExIm_mi_EzRe*ExIm_mi_EzReT);
				ImMI = 0.5*(ExIm_mi_EzRe*ExRe_p_EzImT - ExRe_p_EzIm*ExIm_mi_EzReT);
				break;
			}
			case -1: // s0
			{
				ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
				ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
				break;
			}
			case -2: // s1
			{
				ReMI = ExRe*ExReT + ExIm*ExImT - (EzRe*EzReT + EzIm*EzImT);
				ImMI = ExIm*ExReT - ExRe*ExImT - (EzIm*EzReT - EzRe*EzImT);
				break;
			}
			case -3: // s2
			{
				ReMI = ExImT*EzIm + ExIm*EzImT + ExReT*EzRe + ExRe*EzReT;
				ImMI = ExReT*EzIm - ExRe*EzImT - ExImT*EzRe + ExIm*EzReT;
				break;
			}
			case -4: // s3
			{
				ReMI = ExReT*EzIm + ExRe*EzImT - ExImT*EzRe - ExIm*EzReT;
				ImMI = ExIm*EzImT - ExImT*EzIm - ExReT*EzRe + ExRe*EzReT;
				break;
			}
			default: // total mutual intensity, same as s0
			{
				ReMI = ExRe*ExReT + ExIm*ExImT + EzRe*EzReT + EzIm*EzImT;
				ImMI = ExIm*ExReT - ExRe*ExImT + EzIm*EzReT - EzRe*EzImT;
				break;
				//return CAN_NOT_EXTRACT_MUT_INT;
			}
		}
		if(iter == 0)
		{
			*pResMI = (float)ReMI;
			*(pResMI+1) = (float)ImMI;
		}
		else if(iter > 0)
		{
			double iter_p_1 = iter + 1; //OC20012020
			//long long iter_p_1 = iter + 1;
			*pResMI = (float)(((*pResMI)*iter + ReMI)/iter_p_1);
			*(pResMI+1) = (float)(((*(pResMI+1))*iter + ImMI)/iter_p_1);
		}
		else
		{
			*pResMI += (float)ReMI;
			*(pResMI+1) += (float)ImMI;
		}
		return 0;
	}

	static double IntCylCrd(double ph, void* par)
	{//Function to be called by 1D integration routine
		srTAuxInt2DIntegOverAzim *p = (srTAuxInt2DIntegOverAzim*)par;
		double r = p->r;
		double xx = p->x0 + r*cos(ph), yy = p->y0 + r*sin(ph);
		bool pointToBeUsed = true;
		if(p->arRectToSkipPar != 0) pointToBeUsed = CheckIfPointIsOutsideRectangles(xx, yy, p->arRectToSkipPar, p->nRectToSkip);
		if(pointToBeUsed)
		{
			if(p->pfI2D) return CGenMathInterp::InterpOnRegMesh2d(xx, yy, p->xMin, p->xStep, p->nx, p->yMin, p->yStep, p->ny, p->pfI2D, p->ordInterp);
			else return CGenMathInterp::InterpOnRegMesh2d(xx, yy, p->xMin, p->xStep, p->nx, p->yMin, p->yStep, p->ny, p->pdI2D, p->ordInterp);
		}
		else return 0.;
		//if(p->pfI2D) return CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), p->xMin, p->xStep, p->nx, p->yMin, p->yStep, p->ny, p->pfI2D, p->ordInterp);
		//else return CGenMathInterp::InterpOnRegMesh2d(x0 + r*cos(ph), y0 + r*sin(ph), p->xMin, p->xStep, p->nx, p->yMin, p->yStep, p->ny, p->pdI2D, p->ordInterp);
	}

	static bool CheckIfPointIsOutsideRectangles(double x, double y, double* arRectPar, int nRect)
	{
		double x0, halfDx, y0, halfDy;
		double *tRectPar = arRectPar;
		for(int i=0; i<nRect; i++)
		{
			x0 = *(tRectPar++); halfDx = *(tRectPar++);
			y0 = *(tRectPar++); halfDy = *(tRectPar++);
			if(((x0 - halfDx) <= x) && (x <= (x0 + halfDx)) && ((y0 - halfDy) <= y) && (y <= (y0 + halfDy))) return false;
		}
		return true;
	}

	static void IntProc(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar, int nPar); //OC09032019
	//static void IntProc(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar);
	static void Int2DIntegOverAzim(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar, int nPar); //OC09032019
	//static void Int2DIntegOverAzim(srTWaveAccessData* pwI1, srTWaveAccessData* pwI2, double* arPar);

	static void ComponInteg(srTDataMD* pIntensOrigData, srTDataMD* pIntegParData, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEn(srTDataMD* pIntensOrigData, double eMin, double eMax, srTDataMD* pIntegResData);
	static void ComponIntegVsHorOrVertPos(srTDataMD* pIntensOrigData, char x_or_z, double xMin, double xMax, srTDataMD* pIntegResData);
	static void ComponIntegVsHorAndVertPos(srTDataMD* pIntensOrigData, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEnAndPos(srTDataMD* pIntensOrigData, char x_or_z, double zMin, double zMax, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEnAndHorAndVertPos(srTDataMD* pIntensOrigData, double eMin, double eMax, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData);
};

//*************************************************************************

#endif
