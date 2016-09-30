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

#include <complex>

//*************************************************************************

extern srTIntVect gVectWarnNos;

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

	void ExtractRadiation(int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z, char* pData);
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
		//	if((RadExtract.Int_or_Phase == 2) && (RadExtract.PlotType == 3)) TryToMakePhaseContinuous(ExtractedWaveData);

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

	int TryToMakePhaseContinuous(srTWaveAccessData& WaveData);
	void TryToMakePhaseContinuous1D(double* pOutPhase, long Np, long i0, float Phi0);

	int ExtractSingleElecIntensity1DvsE(srTRadExtract&);
	int ExtractSingleElecIntensity1DvsX(srTRadExtract&);
	int ExtractSingleElecIntensity1DvsZ(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsXZ(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsEX(srTRadExtract&);
	int ExtractSingleElecIntensity2DvsEZ(srTRadExtract&);
	int ExtractSingleElecIntensity3D(srTRadExtract&);
	int SetupExtractedWaveData(srTRadExtract&, srTWaveAccessData&);
	void SetupIntCoord(char, double, long&, long&, double&);

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
	int ConvoluteWithElecBeamOverTransvCoord(float*, long, long);
	void PutConstPhotEnergySliceInExtractPlace(long, long, long, srTRadExtract&, srTRadExtract&);
	void PropagateElecBeamMoments(srTElecBeamMoments&);
	static void PropagateElecBeamMoments(srTElecBeamMoments& ElecBeamMom, double* p4x4PropMatr, double* p4Vect);

	void PadImZerosToRealData(float*, long, long);
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
			case 4: // Circ. Right
			{
				float ExRe_mi_EzIm = ExRe - EzIm, ExIm_p_EzRe = ExIm + EzRe;
				//return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzIm*ExRe_mi_EzIm + ExIm_p_EzRe*ExIm_p_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzIm, ExIm_p_EzRe) : 0.70710678*ExRe_mi_EzIm)));
				return (float)((((Int_or_ReE == 0) || (Int_or_ReE == 1))? 0.5*(ExRe_mi_EzIm*ExRe_mi_EzIm + ExIm_p_EzRe*ExIm_p_EzRe) : ((Int_or_ReE == 2)? FormalPhase(ExRe_mi_EzIm, ExIm_p_EzRe) : 0.70710678*((Int_or_ReE == 3)? ExRe_mi_EzIm : ExIm_p_EzRe)))); //OC031208
			}
			case 5: // Circ. Left
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
				return (float)(-2.*(ExRe*EzRe + ExIm*EzIm));
			}
			case -4: // s3
			{
				return (float)(2.*(-ExRe*EzIm + ExIm*EzRe));
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

	static void ComponInteg(srTDataMD* pIntensOrigData, srTDataMD* pIntegParData, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEn(srTDataMD* pIntensOrigData, double eMin, double eMax, srTDataMD* pIntegResData);
	static void ComponIntegVsHorOrVertPos(srTDataMD* pIntensOrigData, char x_or_z, double xMin, double xMax, srTDataMD* pIntegResData);
	static void ComponIntegVsHorAndVertPos(srTDataMD* pIntensOrigData, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEnAndPos(srTDataMD* pIntensOrigData, char x_or_z, double zMin, double zMax, srTDataMD* pIntegResData);
	static void ComponIntegVsPhotEnAndHorAndVertPos(srTDataMD* pIntensOrigData, double eMin, double eMax, double xMin, double xMax, double zMin, double zMax, srTDataMD* pIntegResData);
};

//*************************************************************************

#endif
