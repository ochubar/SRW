/************************************************************************//**
 * File: sroptzp.h
 * Description: Optical element: Zone Plate (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTZP_H
#define __SROPTZP_H

#include "sroptfoc.h"

//*************************************************************************

class srTZonePlate : public srTFocusingElem {

	int Nzones;
	double RnMax;
	double AttenLen1, AttenLen2;
	double RefrDelta1, RefrDelta2;
	double Thickness;

	double m_ZoneHeightRatioExtToCen, m_ZoneHeightRatioIntermedToCen1, m_ZoneHeightRatioIntermedToCen2;
	int m_ZoneIntermedNum1, m_ZoneIntermedNum2;

	double RnMaxe2;
	double m_aModH, m_bModH, m_cModH, m_dModH;
	bool m_ModH_IsDefined;

public:

	srTZonePlate(srTStringVect* pElemInfo)
	{
		Nzones = atoi((*pElemInfo)[1]);
		RnMax = atof((*pElemInfo)[2]); // input in m
		Thickness = atof((*pElemInfo)[3]); // input in m
		AttenLen1 = atof((*pElemInfo)[4]); // input in 1/m
		RefrDelta1 = atof((*pElemInfo)[5]);
		AttenLen2 = atof((*pElemInfo)[6]); // input in 1/m
		RefrDelta2 = atof((*pElemInfo)[7]);

		if(pElemInfo->size() > 8)
		{
            TransvCenPoint.x = atof((*pElemInfo)[8]); // input in m
			TransvCenPoint.y = atof((*pElemInfo)[9]); // input in m
		}

        m_ZoneHeightRatioExtToCen = m_ZoneHeightRatioIntermedToCen1 = m_ZoneHeightRatioIntermedToCen2 = -1;
        m_ZoneIntermedNum1 = m_ZoneIntermedNum2 = -1;

		int AmOfStr = (int)pElemInfo->size();

		if(AmOfStr > 10) m_ZoneHeightRatioExtToCen = atof((*pElemInfo)[10]); // ratio
		if(AmOfStr > 11) m_ZoneIntermedNum1 = atoi((*pElemInfo)[11]); // number
		if(AmOfStr > 12) m_ZoneHeightRatioIntermedToCen1 = atof((*pElemInfo)[12]); // ratio
		if(AmOfStr > 13) m_ZoneIntermedNum2 = atoi((*pElemInfo)[13]); // number
		if(AmOfStr > 14) m_ZoneHeightRatioIntermedToCen2 = atof((*pElemInfo)[14]); // ratio

		RnMaxe2 = RnMax*RnMax;

		DefineAttenModulConstants();
	}
	srTZonePlate(int _nZones, double _rn, double _thick, double _atLen1, double _atLen2, double _delta1, double _delta2, double _x=0, double _y=0)
	{
		Nzones = _nZones;
		RnMax = _rn; // input in m
		Thickness = _thick; // input in m
		AttenLen1 = _atLen1; // input in m
		RefrDelta1 = _delta1;
		AttenLen2 = _atLen2; // input in m
		RefrDelta2 = _delta2;

		TransvCenPoint.x = _x; // input in m
		TransvCenPoint.y = _y; // input in m

        m_ZoneHeightRatioExtToCen = m_ZoneHeightRatioIntermedToCen1 = m_ZoneHeightRatioIntermedToCen2 = -1;
        m_ZoneIntermedNum1 = m_ZoneIntermedNum2 = -1;

		RnMaxe2 = RnMax*RnMax;

		DefineAttenModulConstants();
	}
	srTZonePlate() {}

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			pRadAccessData->CheckAndSubtractPhaseTermsLin(TransvCenPoint.x, TransvCenPoint.y);
		//}

		char &MethNo = ParPrecWfrPropag.MethNo;
		
		int result = 0;

		if(MethNo == 0) result = PropagateRadiationMeth_0(pRadAccessData);
		//else return PropagateRadiationMeth_2(pRadAccessData, ResBeforeAndAfterVect);
		else result = PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);

		//if(ParPrecWfrPropag.AnalTreatment == 1)
		//{// Treating linear terms analytically
			if(!ParPrecWfrPropag.DoNotResetAnalTreatTermsAfterProp) pRadAccessData->CheckAndResetPhaseTermsLin();
		//}

		return result;
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;
		return 0;
	}
  	int PropagateRadiationSimple1D(srTRadSect1D* pSect1D)
	{
		int result;
		if(pSect1D->Pres != 0) if(result = SetRadRepres1D(pSect1D, 0)) return result;
		if(result = TraverseRad1D(pSect1D)) return result;
		return 0;
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{
		SetupFocalDistForPhotonEnergy(pRadAccessData->eStart);
		return srTFocusingElem::PropagateRadMoments(pRadAccessData, MomRatArray);
	}

	void SetupFocalDistForPhotonEnergy(double ePh)
	{
        FocDistX = FocDistZ = 1.E+23;
		if(ePh <= 0.) return;

		//double Wavelength_m = (1.239842E-06)/ePh; // assuming ePh in eV
		//FocDistX = RnMax*RnMax/(Nzones*Wavelength_m); 
		//FocDistZ = FocDistX;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
		double re2 = xRel*xRel + zRel*zRel;
		//int par = ZoneParity(re2);
		//if(par == 0) return;

		double AttenLen = AttenLen1, RefrDelta = RefrDelta1;
		double AttenLenComplem = AttenLen2, RefrDeltaComplem = RefrDelta2;
		double AmpAtten = 1, OptPathDiff = 0;
		float *PExRe = EPtrs.pExRe, *PExIm = EPtrs.pExIm, *PEzRe = EPtrs.pEzRe, *PEzIm = EPtrs.pEzIm;

        if(re2 > RnMaxe2) 
		{
			bool LastZoneNumIsEven = (Nzones == ((Nzones >> 1) << 1));
			if(LastZoneNumIsEven)
			{
                AttenLen = AttenLen1; RefrDelta = RefrDelta1;
			}
			else
			{
                AttenLen = AttenLen2; RefrDelta = RefrDelta2;
			}
            AmpAtten = exp(-0.5*Thickness/AttenLen);
            OptPathDiff = RefrDelta*Thickness;
		}
		else
		{
			int CurZoneNum = int(re2*Nzones/RnMaxe2) + 1;

			bool CurZoneNumIsEven = (CurZoneNum == ((CurZoneNum >> 1) << 1));
			int par = (CurZoneNumIsEven? 2 : 1);

			if(par == 2) 
			{
				AttenLen = AttenLen2; RefrDelta = RefrDelta2;
				AttenLenComplem = AttenLen1; RefrDeltaComplem = RefrDelta1;
			}

			double CurHeight = Thickness;
			double CurHeightComplem = 0;
			if(m_ModH_IsDefined && CurZoneNumIsEven)
			{
				double rnm1 = sqrt(RnMaxe2*(CurZoneNum - 1)/Nzones), rn = sqrt(RnMaxe2*CurZoneNum/Nzones);
				double rZoneCen = 0.5*(rnm1 + rn);
				CurHeight = ((m_aModH*rZoneCen + m_bModH)*rZoneCen + m_cModH)*rZoneCen + m_dModH;
				CurHeightComplem = Thickness - CurHeight;
			}

			//if(AttenLen < 1.e-13) { *PExRe = 0; *PExIm = 0; *PEzRe = 0; *PEzIm = 0; return;}

			//double AmpAtten = exp(-0.5*Thickness/AttenLen);
			//double AmpAtten = exp(-0.5*CurHeight/AttenLen);
			AmpAtten = exp(-0.5*(CurHeight/AttenLen + CurHeightComplem/AttenLenComplem));
			OptPathDiff = RefrDelta*CurHeight + RefrDeltaComplem*CurHeightComplem;
		}

		double k_inv_m = EXZ.e*(5.067681604e+06);
		//double PhaseShift = -k_inv_m*RefrDelta*Thickness;
		//double PhaseShift = -k_inv_m*RefrDelta*CurHeight;
		double PhaseShift = -k_inv_m*OptPathDiff;

		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);
		if((PExRe != 0) && (PExIm != 0))
		{
			float NewExRe = (float)(AmpAtten*((*PExRe)*CosPh - (*PExIm)*SinPh));
			float NewExIm = (float)(AmpAtten*((*PExRe)*SinPh + (*PExIm)*CosPh));
			*PExRe = NewExRe; *PExIm = NewExIm; 
		}
		if((PEzRe != 0) && (PEzIm))
		{
			float NewEzRe = (float)(AmpAtten*((*PEzRe)*CosPh - (*PEzIm)*SinPh));
			float NewEzIm = (float)(AmpAtten*((*PEzRe)*SinPh + (*PEzIm)*CosPh));
			*PEzRe = NewEzRe; *PEzIm = NewEzIm; 
		}
	}

  	void RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coord. side !!!
		double ArgRel = (EXZ.VsXorZ == 'x')? (EXZ.x - TransvCenPoint.x) : (EXZ.z - TransvCenPoint.y);
		double re2 = ArgRel*ArgRel;
		//int par = ZoneParity(re2);
		//if(par == 0) return;

		double AttenLen = AttenLen1, RefrDelta = RefrDelta1;
		double AttenLenComplem = AttenLen2, RefrDeltaComplem = RefrDelta2;
		double AmpAtten = 1, OptPathDiff = 0;
		float *PExRe = EPtrs.pExRe, *PExIm = EPtrs.pExIm, *PEzRe = EPtrs.pEzRe, *PEzIm = EPtrs.pEzIm;

        if(re2 > RnMaxe2)
		{
			bool LastZoneNumIsEven = (Nzones == ((Nzones >> 1) << 1));
			if(LastZoneNumIsEven)
			{
                AttenLen = AttenLen1; RefrDelta = RefrDelta1;
			}
			else
			{
                AttenLen = AttenLen2; RefrDelta = RefrDelta2;
			}
            AmpAtten = exp(-0.5*Thickness/AttenLen);
            OptPathDiff = RefrDelta*Thickness;
		}
		else
		{
			int CurZoneNum = int(re2*Nzones/RnMaxe2) + 1;
			bool CurZoneNumIsEven = (CurZoneNum == ((CurZoneNum >> 1) << 1));
			int par = (CurZoneNumIsEven? 2 : 1);

			//double AttenLen = AttenLen1, RefrDelta = RefrDelta1;
			//double AttenLenComplem = AttenLen2, RefrDeltaComplem = RefrDelta2;
			if(par == 2) 
			{ 
				AttenLen = AttenLen2; RefrDelta = RefrDelta2;
				AttenLenComplem = AttenLen1; RefrDeltaComplem = RefrDelta1;
			}

			double CurHeight = Thickness;
			double CurHeightComplem = 0;
			if(m_ModH_IsDefined && CurZoneNumIsEven)
			{
				double rnm1 = sqrt(RnMaxe2*(CurZoneNum - 1)/Nzones), rn = sqrt(RnMaxe2*CurZoneNum/Nzones);
				double rZoneCen = 0.5*(rnm1 + rn);
				CurHeight = ((m_aModH*rZoneCen + m_bModH)*rZoneCen + m_cModH)*rZoneCen + m_dModH;
				CurHeightComplem = Thickness - CurHeight;
			}

			//float *PExRe = EPtrs.pExRe, *PExIm = EPtrs.pExIm, *PEzRe = EPtrs.pEzRe, *PEzIm = EPtrs.pEzIm;
			//if(AttenLen < 1.e-13) { *PExRe = 0; *PExIm = 0; *PEzRe = 0; *PEzIm = 0; return;}

			//double AmpAtten = exp(-0.5*Thickness/AttenLen);
			//double AmpAtten = exp(-0.5*CurHeight/AttenLen);
			AmpAtten = exp(-0.5*(CurHeight/AttenLen + CurHeightComplem/AttenLenComplem));
			OptPathDiff = CurHeight/AttenLen + CurHeightComplem/AttenLenComplem;
		}

		double k_inv_m = EXZ.e*(5.067681604e+06);
		//double PhaseShift = -k_inv_m*RefrDelta*Thickness;
		//double PhaseShift = -k_inv_m*RefrDelta*CurHeight;
		double PhaseShift = -k_inv_m*OptPathDiff;

		float CosPh, SinPh; CosAndSin(PhaseShift, CosPh, SinPh);

		if((PExRe != 0) && (PExIm != 0))
		{
			float NewExRe = (float)(AmpAtten*((*PExRe)*CosPh - (*PExIm)*SinPh));
			float NewExIm = (float)(AmpAtten*((*PExRe)*SinPh + (*PExIm)*CosPh));
			*PExRe = NewExRe; *PExIm = NewExIm; 
		}
		if((PEzRe != 0) && (PEzIm != 0))
		{
			float NewEzRe = (float)(AmpAtten*((*PEzRe)*CosPh - (*PEzIm)*SinPh));
			float NewEzIm = (float)(AmpAtten*((*PEzRe)*SinPh + (*PEzIm)*CosPh));
			*PEzRe = NewEzRe; *PEzIm = NewEzIm; 
		}
	}

	double RadOptPathDiff(srTEXZ& EXZ) //virtual 
	{// e in eV; Length in m !!!
		double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
		double re2 = xRel*xRel + zRel*zRel;

		double RefrDelta = RefrDelta1, RefrDeltaCen = RefrDelta1;
		double RefrDeltaComplem = RefrDelta2; //, RefrDeltaCenComplem = RefrDelta2;

		double CurHeight = Thickness, HeightCen = Thickness;
		double CurHeightComplem = 0, HeightCenComplem = 0;

        if(re2 > RnMaxe2) 
		{
			bool LastZoneNumIsEven = (Nzones == ((Nzones >> 1) << 1));
			if(LastZoneNumIsEven)
			{
                RefrDelta = RefrDelta1;
			}
			else
			{
                RefrDelta = RefrDelta2;
			}
			return -RefrDelta*Thickness + RefrDeltaCen*HeightCen;
		}

        int CurZoneNum = int(re2*Nzones/RnMaxe2) + 1;
        bool CurZoneNumIsEven = (CurZoneNum == ((CurZoneNum >> 1) << 1));
		int par = (CurZoneNumIsEven? 2 : 1);

		if(par == 2) 
		{
			RefrDelta = RefrDelta2;
			RefrDeltaComplem = RefrDelta1;
		}

		if(m_ModH_IsDefined && CurZoneNumIsEven)
		{
			double rnm1 = sqrt(RnMaxe2*(CurZoneNum - 1)/Nzones), rn = sqrt(RnMaxe2*CurZoneNum/Nzones);
			double rZoneCen = 0.5*(rnm1 + rn);
            CurHeight = ((m_aModH*rZoneCen + m_bModH)*rZoneCen + m_cModH)*rZoneCen + m_dModH;
			CurHeightComplem = Thickness - CurHeight;
		}
		//return -RefrDelta*CurHeight + RefrDeltaCen*HeightCen;
		return -RefrDelta*CurHeight - RefrDeltaComplem*CurHeightComplem + RefrDeltaCen*HeightCen;
	}

	int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
	{
        MinNx = MinNz = 4*Nzones;
		return 0;
	}

	int ZoneParity(double re2)
	{
		if(re2 > RnMaxe2) return 0;
		int n = int(re2*Nzones/RnMaxe2) + 1;
		bool nIsEven = (n == ((n >> 1) << 1));
		return nIsEven? 2 : 1;
	}

	int RangeShouldBeAdjustedAtPropag() { return 0;} // Or switch it On
	int ResolutionShouldBeAdjustedAtPropag() { return 1;}
	int AllowAutoSwitchToUndersamplingMode() { return 0;}

	void DefineAttenModulConstants();
};

//*************************************************************************

#endif

