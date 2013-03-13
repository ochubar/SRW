/************************************************************************//**
 * File: sroptgrat.cpp
 * Description: Optical element: Grating
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptgrat.h"

//*************************************************************************

void srTGrating::SetupPropBufVars_Gen(srTSRWRadStructAccessData* pWfr)
{// Compute any buf. variables necessary for propagation

	double Lambda0 = 1.239842e-06/pWfr->avgPhotEn;
	double Lambda0_d_d0 = m_Order*Lambda0/m_Period;

	m_PropBufVars.ThetaIt = HalfPI - Theta;
	m_PropBufVars.SinThetaI = sin(Theta);
	m_PropBufVars.ThetaM0 = asin(Lambda0_d_d0 - m_PropBufVars.SinThetaI);

	double ThetaM0_p_ThetaIt = m_PropBufVars.ThetaM0 + m_PropBufVars.ThetaIt;

	m_PropBufVars.Sin_ThetaM0_p_ThetaIt = sin(ThetaM0_p_ThetaIt);
	m_PropBufVars.Cos_ThetaM0_p_ThetaIt = cos(ThetaM0_p_ThetaIt);
	m_PropBufVars.Tg_ThetaIt = tan(m_PropBufVars.ThetaIt);

	//double WfrRange = (RotPlane == 'h')? (pWfr->xWfrMax - pWfr->xWfrMin) : (pWfr->zWfrMax - pWfr->zWfrMin);
	//double SpotSizeOnGrating = WfrRange/sin(m_PropBufVars.ThetaIt);

	m_PropBufVars.wfrR = (RotPlane == 'h')? pWfr->RobsX : pWfr->RobsZ;

	//m_PropBufVars.L2 = ::fabs(SpotSizeOnGrating*sin(m_PropBufVars.ThetaM0));
	m_PropBufVars.L2 = 0.;
	//or make m_PropBufVars.L2 = 0; ?

	m_PropBufVars.td_NominTermL2 = m_PropBufVars.L2*(m_PropBufVars.Tg_ThetaIt*m_PropBufVars.Sin_ThetaM0_p_ThetaIt + m_PropBufVars.Cos_ThetaM0_p_ThetaIt);
	m_PropBufVars.td_NominMultX2 = m_PropBufVars.Tg_ThetaIt*m_PropBufVars.Cos_ThetaM0_p_ThetaIt - m_PropBufVars.Sin_ThetaM0_p_ThetaIt;

	double dAux = -Lambda0_d_d0 + cos(m_PropBufVars.ThetaIt);
	
	m_PropBufVars.OptPathCorTiltMultX2 = Lambda0_d_d0/sqrt(1. - dAux*dAux);
	m_PropBufVars.ReflectAmp = sqrt(m_ReflectAvgInt);

	//to check:
	m_PropBufVars.AnamorphMagn = ::fabs(sin(HalfPI + m_PropBufVars.ThetaM0)/sin(m_PropBufVars.ThetaIt));
	m_PropBufVars.PowerConservMultE = m_PropBufVars.ReflectAmp/sqrt(m_PropBufVars.AnamorphMagn);
}

//*************************************************************************

void srTGrating::SetupPropBufVars_SingleE(double PhotEn)
{//call only after SetupPropBufVars_Gen !!!

	m_PropBufVars.CurPhotEn = PhotEn;
	m_PropBufVars.CurWaveNumb = (5.06773065E+06)*PhotEn; //(TwoPI/1.239842e-06)*PhotEn;

	m_PropBufVars.Lambda = 1.239842e-06/PhotEn;
	m_PropBufVars.ThetaM = asin(m_Order*m_PropBufVars.Lambda/m_Period - m_PropBufVars.SinThetaI);
	double ThetaM_p_ThetaIt = m_PropBufVars.ThetaM + m_PropBufVars.ThetaIt;

	m_PropBufVars.Sin_ThetaM_p_ThetaIt = sin(ThetaM_p_ThetaIt);
	m_PropBufVars.Cos_ThetaM_p_ThetaIt = cos(ThetaM_p_ThetaIt);
	//CosAndSin(ThetaM_p_ThetaIt, m_PropBufVars.Cos_ThetaM_p_ThetaIt, m_PropBufVars.Sin_ThetaM_p_ThetaIt);

	m_PropBufVars.td_MultInvDenom = 1./(m_PropBufVars.Tg_ThetaIt*m_PropBufVars.Sin_ThetaM_p_ThetaIt + m_PropBufVars.Cos_ThetaM_p_ThetaIt);

	m_PropBufVars.WaveFrontTermWasTreated = false;
	//if(WaveFrontTermCanBeTreated(*m_pPrevWfr)) //OC120908
	//{
	//	TreatStronglyOscillatingTerm(*m_pPrevWfr, 'r');
	//	m_PropBufVars.WaveFrontTermWasTreated = true;

	//	//double Const = PI*1.E+06/1.239842; // Assumes m and eV
	//	//m_PropBufVars.ConstRxE = Const*PhotEn/(m_pPrevWfr->RobsX);
	//	//m_PropBufVars.ConstRzE = Const*PhotEn/(m_pPrevWfr->RobsZ); //OC171008
	//}
}

//*************************************************************************
/**Old / original version
void srTGrating::RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
{// e in eV; Length in m !!!
 // Operates on Coordinate side !!!

	if(EXZ.e != m_PropBufVars.CurPhotEn) SetupPropBufVars_SingleE(EXZ.e);

	double x2 = (RotPlane == 'h')? EXZ.x : EXZ.z;

	double td = (m_PropBufVars.td_NominTermL2 + x2*m_PropBufVars.td_NominMultX2)*m_PropBufVars.td_MultInvDenom;
	//double rgz = -m_PropBufVars.L2*m_PropBufVars.Sin_ThetaM0_p_ThetaIt - x2*m_PropBufVars.Cos_ThetaM0_p_ThetaIt + td*m_PropBufVars.Sin_ThetaM_p_ThetaIt;
	double rgx = m_PropBufVars.L2*m_PropBufVars.Cos_ThetaM0_p_ThetaIt - x2*m_PropBufVars.Sin_ThetaM0_p_ThetaIt - td*m_PropBufVars.Cos_ThetaM_p_ThetaIt;

	//double OptPathDif = td + rgz - m_PropBufVars.L2 + x2*m_PropBufVars.OptPathCorTiltMultX2;
	//double PhaseShift = m_PropBufVars.CurWaveNumb*OptPathDif;

	//Simplified version: just an angle
	double AngDisp = m_PropBufVars.ThetaM - m_PropBufVars.ThetaM0;
	double PhaseShift = m_PropBufVars.CurWaveNumb*AngDisp*x2;

		//from prop. multi-e:
		//double WaveNumb = EXZ.e*kMult;
		//double Ph = WaveNumb*(AngX*EXZ.x + AngZ*EXZ.z);

	float CosPh, SinPh;
	CosAndSin(PhaseShift, CosPh, SinPh);

	double x1=0, z1=0, rArg=0; //for vertical dispersion
	int nx_mi_1 = m_pPrevWfr->nx - 1, nz_mi_1 = m_pPrevWfr->nz - 1;
	int ie = 0;
	if((m_pPrevWfr->ne > 1) && (m_pPrevWfr->eStep != 0))
	{
		ie = (int)((EXZ.e - m_pPrevWfr->eStart)/(m_pPrevWfr->eStep) + 0.000001);
	}
	long two_ie = ie << 1;

	long PerX = (m_pPrevWfr->ne) << 1;
	long PerZ = (m_pPrevWfr->nx)*PerX;
	long ofstPrev0 = 0, ofstPrev1 = 0;

	if(RotPlane == 'h')
	{
		x1 = rgx; z1 = EXZ.z;
		int ix0 = (int)((x1 - m_pPrevWfr->xStart)/(m_pPrevWfr->xStep) + 0.000001);
		if(ix0 < 0) ix0 = 0;
		else if(ix0 >= nx_mi_1) ix0 = nx_mi_1 - 1;
		int ix1 = ix0 + 1;
		int izc = (int)((z1 - m_pPrevWfr->zStart)/(m_pPrevWfr->zStep) + 0.000001);
		if(izc < 0) izc = 0;
		else if(izc >= m_pPrevWfr->nz) izc = m_pPrevWfr->nz - 1;

		long izc_PerZ_p_two_ie = izc*PerZ + two_ie;
		ofstPrev0 = izc_PerZ_p_two_ie + ix0*PerX;
		ofstPrev1 = izc_PerZ_p_two_ie + ix1*PerX;
		rArg = (x1 - (m_pPrevWfr->xStart + (m_pPrevWfr->xStep)*ix0))/(m_pPrevWfr->xStep);
	}
	else
	{
		x1 = EXZ.x; z1 = rgx;
		int ixc = (int)((x1 - m_pPrevWfr->xStart)/(m_pPrevWfr->xStep) + 0.000001);
		if(ixc < 0) ixc = 0;
		else if(ixc >= m_pPrevWfr->nx) ixc = m_pPrevWfr->nx - 1;
		int iz0 = (int)((z1 - m_pPrevWfr->zStart)/(m_pPrevWfr->zStep) + 0.000001);
		if(iz0 < 0) iz0 = 0;
		else if(iz0 >= nz_mi_1) iz0 = nz_mi_1 - 1;
		int iz1 = iz0 + 1;

		long ixc_PerX_p_two_ie = ixc*PerX + two_ie;
		ofstPrev0 = iz0*PerZ + ixc_PerX_p_two_ie;
		ofstPrev1 = iz1*PerZ + ixc_PerX_p_two_ie;
		rArg = (z1 - (m_pPrevWfr->zStart + (m_pPrevWfr->zStep)*iz0))/(m_pPrevWfr->zStep);
	}
	float *pPrevEx0 = m_pPrevWfr->pBaseRadX + ofstPrev0;
	float *pPrevEx1 = m_pPrevWfr->pBaseRadX + ofstPrev1;
	float *pPrevEz0 = m_pPrevWfr->pBaseRadZ + ofstPrev0;
	float *pPrevEz1 = m_pPrevWfr->pBaseRadZ + ofstPrev1;

	double QuadPhaseTerm = 0;
	float CosQuadPhaseTerm, SinQuadPhaseTerm;
	if(m_PropBufVars.WaveFrontTermWasTreated)
	{//add back quad. phase term
		double relX = x1 - m_pPrevWfr->xc, relZ = z1 - m_pPrevWfr->zc;
		QuadPhaseTerm = relX*relX*m_PropBufVars.ConstRxE + relZ*relZ*m_PropBufVars.ConstRzE;
		CosAndSin(QuadPhaseTerm, CosQuadPhaseTerm, SinQuadPhaseTerm);
	}

	double PrevExRe=0, PrevExIm=0, PrevEzRe=0, PrevEzIm;
	if(m_pPrevWfr->pBaseRadX != 0)
	{
		PrevExRe = InterpLin(rArg, *pPrevEx0, *pPrevEx1);
		PrevExIm = InterpLin(rArg, *(pPrevEx0+1), *(pPrevEx1+1));

		if(m_PropBufVars.WaveFrontTermWasTreated)
		{//add back quad. phase term
			double AuxExRe = PrevExRe*CosQuadPhaseTerm - PrevExIm*SinQuadPhaseTerm;
			PrevExIm = PrevExRe*SinQuadPhaseTerm + PrevExIm*CosQuadPhaseTerm;
			PrevExRe = AuxExRe;
		}

		*(EPtrs.pExRe) = (PrevExRe*CosPh - PrevExIm*SinPh)*m_PropBufVars.PowerConservMultE;
		*(EPtrs.pExIm) = (PrevExRe*SinPh + PrevExIm*CosPh)*m_PropBufVars.PowerConservMultE;
	}
	if(m_pPrevWfr->pBaseRadZ != 0)
	{
		PrevEzRe = InterpLin(rArg, *pPrevEz0, *pPrevEz1);
		PrevEzIm = InterpLin(rArg, *(pPrevEz0+1), *(pPrevEz1+1));

		if(m_PropBufVars.WaveFrontTermWasTreated)
		{//add back quad. phase term
			double AuxEzRe = PrevEzRe*CosQuadPhaseTerm - PrevEzIm*SinQuadPhaseTerm;
			PrevEzIm = PrevEzRe*SinQuadPhaseTerm + PrevEzIm*CosQuadPhaseTerm;
			PrevEzRe = AuxEzRe;
		}

		*(EPtrs.pEzRe) = (PrevEzRe*CosPh - PrevEzIm*SinPh)*m_PropBufVars.PowerConservMultE;
		*(EPtrs.pEzIm) = (PrevEzRe*SinPh + PrevEzIm*CosPh)*m_PropBufVars.PowerConservMultE;
	}
}
**/
//*************************************************************************
