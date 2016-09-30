/************************************************************************//**
 * File: sroptcryst.h
 * Description: Optical element: Crystal (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2013
 *
 * Copyright (C) Brookhaven National Laboratory, Upton, NY, USA
 * Copyright (C) Diamond Light Source, UK
 * All Rights Reserved
 *
 * @author J.Sutter, A.Suvorov, O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTCRYST_H
#define __SROPTCRYST_H

#include "sroptelm.h"
#include "srwlib.h"

#undef max //to avoid name conflict of numeric_limits<double>::max() and #define max(a, b) to allow for compilation with VC++2013
#include <limits>
#include <complex>

//*************************************************************************

struct srTOptCrystMeshTrf {
	double xStart, xStep;
	double zStart, zStep;
	double matrKLabRef[2][3];
	bool crossTermsAreLarge, xMeshTrfIsReq, zMeshTrfIsReq; //, dataTranspIsReq;
};

//*************************************************************************

class srTOptCryst : public srTGenOptElem {

	double m_dA; /* crystal reflecting planes d-spacing (units?) */
	//double m_psi0r, m_psi0i; /* real and imaginary parts of 0-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psi0c;
	//double m_psiHr, m_psiHi; /* real and imaginary parts of H-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psihc;
	//double m_psiHbr, m_psiHbi; /* real and imaginary parts of -H-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psimhc;
	//double m_hMilND, m_kMilND, m_lMilND; /* 1st, 2nd and 3rd  indexes of diffraction vector (Miller indices) */
	//OC180314: Miller indices are removed after discussion with A. Suvorov
	double m_thicum; /* crystal thickness [microns] */
	//double m_alphrd; /* asymmetry angle [rad] */

	char m_itrans; // Input #7: Whether to calculate the transmitted beam as well as the diffracted beam. itrans = 0 for NO, itrans = 1 for YES. 

	//TVector3d m_nv; // horizontal, vertical and longitudinal coordinates of outward normal to crystal surface in the frame of incident beam
	TVector3d m_tv; // horizontal, vertical and longitudinal coordinates of central tangential vector [m] in the frame of incident beam
	TVector3d m_sv; // horizontal, vertical and longitudinal coordinates of central saggital vector [m] in the frame of incident beam
	
	char m_uc; // crystal use case: 1- Bragg Reflection, 2- Bragg Transmission (Laue cases to be added)

	// TRANSFORMATION MATRIX
	double m_RXtLab[3][3]; // RXtLab: 3x3 orthogonal matrix that converts components of a 3x1 vector in crystal coordinates to components in lab coordinates. 
	double m_RLabXt[3][3]; // RLabXt: transpose of RXtLab: converts components of a 3x1 vector in lab coordinates to components in crystal coordinates. 
	double m_RXtRef[3][3]; // RXtRef: 3x3 orthogonal matrix that converts components of a 3x1 vector in crystal coordinates to components in diffracted beam coordinates. 

	//Reflected "horizontal" and "longitudinal" beam reference frame base vectors in the Lab (i.e. the input beam) frame:
	//TVector3d m_vX1, m_vZ1;

	//double m_KLabRef[2][3]; // matrix to transform Kin to Kout

	double m_PolTrn[2][2]; // 2x2 transformation matrix of the polarizations from (e1X,e2X) to (sg0X,pi0X).
	double m_InvPolTrn[2][2]; //OC06092016
	double m_HXAi[3]; // Reciprocal lattice vector coordinates

	double m_sg0X[3];
	double m_cos2t; //cos(2.*thBrd); 

	srTOptCrystMeshTrf *m_pMeshTrf;
	double m_eStartAux, m_eStepAux, m_ne;
	//bool m_xStartIsConstInSlicesE, m_zStartIsConstInSlicesE;

public:

	srTOptCryst(const SRWLOptCryst& srwlCr)
	{
		m_dA = srwlCr.dSp; //crystal reflecting planes d-spacing (units: [A] ?)

		//m_psi0r = srwlCr.psi0r;
		//m_psi0i = srwlCr.psi0i;
		m_psi0c = complex<double>(srwlCr.psi0r, srwlCr.psi0i);
		//m_psiHr = srwlCr.psiHr;
		//m_psiHi = srwlCr.psiHi;
		m_psihc = complex<double>(srwlCr.psiHr, srwlCr.psiHi);
		//m_psiHbr = srwlCr.psiHbr;
		//m_psiHbi = srwlCr.psiHbi;
		m_psimhc = complex<double>(srwlCr.psiHbr, srwlCr.psiHbi);

		//m_hMilND = srwlCr.h1;
		//m_kMilND = srwlCr.h2;
		//m_lMilND = srwlCr.h3;
		//OC180314: the Miller indices are removed after discussion with A. Suvorov (because these are only required for the m_dA, and it is used as input parameter)

		//m_tc = srwlCr.tc;
		m_thicum = srwlCr.tc*1e+06; //assuming srwlCr.tc in [m]
		double alphrd = srwlCr.angAs; //asymmetry angle [rad]

		// From Input #5: reciprocal lattice vector coordinates 
		m_HXAi[0] = 0.;
		m_HXAi[1] = cos(alphrd) / m_dA;
		m_HXAi[2] = -sin(alphrd) / m_dA;

		if(srwlCr.nvz == 0) throw IMPROPER_OPTICAL_COMPONENT_ORIENT;
		
		TVector3d v_nv(srwlCr.nvx, srwlCr.nvy, srwlCr.nvz); /* horizontal, vertical and longitudinal coordinates of outward normal to crystal surface in the frame of incident beam */
		v_nv.Normalize();
		double nv[] = { v_nv.x, v_nv.y, v_nv.z };

		if((srwlCr.tvx == 0) && (srwlCr.tvy == 0)) throw IMPROPER_OPTICAL_COMPONENT_ORIENT;

		//TVector3d v_tv(srwlCr.tvx, srwlCr.tvy, 0); /* horizontal, vertical  and vertical coordinates of central tangential vector [m] in the frame of incident beam */
		//v_tv.z = (-v_nv.x*v_tv.x - v_nv.y*v_tv.y) / v_nv.z;
		//v_tv.Normalize();
		//double tv[] = { v_tv.x, v_tv.y, v_tv.z };
		m_tv.x = srwlCr.tvx; m_tv.y = srwlCr.tvy; /* horizontal and vertical coordinates of central tangential vector [m] in the frame of incident beam */
		m_tv.z = (-v_nv.x*m_tv.x - v_nv.y*m_tv.y) / v_nv.z;
		m_tv.Normalize();
		double tv[] = { m_tv.x, m_tv.y, m_tv.z };

		//sv[0] = nv[1] * tv[2] - nv[2] * tv[1];
		//sv[1] = nv[2] * tv[0] - nv[0] * tv[2];
		//sv[2] = nv[0] * tv[1] - nv[1] * tv[0];
		double sv[] = { nv[1] * tv[2] - nv[2] * tv[1], nv[2] * tv[0] - nv[0] * tv[2], nv[0] * tv[1] - nv[1] * tv[0] };
		m_sv.x = sv[0]; m_sv.y = sv[1]; m_sv.z = sv[2]; 

		m_uc = srwlCr.uc; //OC05092016
		if((m_uc < 1) || (m_uc > 2)) throw IMPROPER_OPTICAL_COMPONENT_MIRROR_USE_CASE;

		// Input #7: Whether to calculate the transmitted beam as well as the diffracted 
		// beam. itrans = 0 for NO, itrans = 1 for YES. 
		m_itrans = 0; //OC: make it input variable
		if(m_uc == 2) m_itrans = 1; //OC05092016

		// RXtLab: 3x3 orthogonal matrix that converts components of a 3x1 vector 
		// in crystal coordinates to components in lab coordinates. 
		for(int i = 0; i < 3; i++)
		{
			m_RXtLab[i][0] = sv[i];
			m_RXtLab[i][1] = nv[i];
			m_RXtLab[i][2] = tv[i];
		}
		// RLabXt: transpose of RXtLab: converts components of a 3x1 vector in lab 
		// coordinates to components in crystal coordinates. 
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++) m_RLabXt[i][j] = m_RXtLab[j][i];

		// Dynamical Bragg angle (angle between incident optical axis and lattice planes)
		double uBrd = m_RLabXt[0][2] * m_HXAi[0] + m_RLabXt[1][2] * m_HXAi[1] + m_RLabXt[2][2] * m_HXAi[2];
		double uHX = sqrt(m_HXAi[0] * m_HXAi[0] + m_HXAi[1] * m_HXAi[1] + m_HXAi[2] * m_HXAi[2]);
		double uRL = sqrt(m_RLabXt[0][2] * m_RLabXt[0][2] + m_RLabXt[1][2] * m_RLabXt[1][2] + m_RLabXt[2][2] * m_RLabXt[2][2]);
		const double pi = 4.*atan(1.);
		double thBrd = acos(uBrd / uHX / uRL) - pi / 2.;
		//printf("Bragg angle = %f\n",thBrd*180./pi); 
		m_cos2t = cos(2.*thBrd);

		// Conversion of polarization vector components to crystal frame: 
		//e1X[0] = RLabXt[0][0];
		//e1X[1] = RLabXt[1][0];
		//e1X[2] = RLabXt[2][0];
		//e2X[0] = RLabXt[0][1];
		//e2X[1] = RLabXt[1][1];
		//e2X[2] = RLabXt[2][1];
		double e1X[] = { m_RLabXt[0][0], m_RLabXt[1][0], m_RLabXt[2][0] };
		double e2X[] = { m_RLabXt[0][1], m_RLabXt[1][1], m_RLabXt[2][1] };

		// Calculation of new polarization vectors in the crystal frame, with the 
		// sigma vector parallel to the diffracting atomic planes of the crystal. 
		//...
		//ucn0X[0] = RLabXt[0][2];
		//ucn0X[1] = RLabXt[1][2];
		//ucn0X[2] = RLabXt[2][2];
		double ucn0X[] = { m_RLabXt[0][2], m_RLabXt[1][2], m_RLabXt[2][2] };

		m_sg0X[0] = m_HXAi[1] * ucn0X[2] - m_HXAi[2] * ucn0X[1];
		m_sg0X[1] = m_HXAi[2] * ucn0X[0] - m_HXAi[0] * ucn0X[2];
		m_sg0X[2] = m_HXAi[0] * ucn0X[1] - m_HXAi[1] * ucn0X[0];

		double sgmag = sqrt(m_sg0X[0] * m_sg0X[0] + m_sg0X[1] * m_sg0X[1] + m_sg0X[2] * m_sg0X[2]);
		//sg0X[0] = sg0X[0] / sgmag;
		//sg0X[1] = sg0X[1] / sgmag;
		//sg0X[2] = sg0X[2] / sgmag;
		for (int i = 0; i < 3; i++) m_sg0X[i] /= sgmag;

		//pi0X[0] = ucn0X[1] * sg0X[2] - ucn0X[2] * sg0X[1];
		//pi0X[1] = ucn0X[2] * sg0X[0] - ucn0X[0] * sg0X[2];
		//pi0X[2] = ucn0X[0] * sg0X[1] - ucn0X[1] * sg0X[0];
		double pi0X[] = {
			ucn0X[1] * m_sg0X[2] - ucn0X[2] * m_sg0X[1],
			ucn0X[2] * m_sg0X[0] - ucn0X[0] * m_sg0X[2],
			ucn0X[0] * m_sg0X[1] - ucn0X[1] * m_sg0X[0]
		};

		// Calculate the 2x2 transformation matrix of the polarizations from 
		// (e1X,e2X) to (sg0X,pi0X). 
		m_PolTrn[0][0] = e1X[0] * m_sg0X[0] + e1X[1] * m_sg0X[1] + e1X[2] * m_sg0X[2];
		m_PolTrn[0][1] = e2X[0] * m_sg0X[0] + e2X[1] * m_sg0X[1] + e2X[2] * m_sg0X[2];
		m_PolTrn[1][0] = e1X[0] * pi0X[0] + e1X[1] * pi0X[1] + e1X[2] * pi0X[2];
		m_PolTrn[1][1] = e2X[0] * pi0X[0] + e2X[1] * pi0X[1] + e2X[2] * pi0X[2];

		double invDetPolTrn = 1./(m_PolTrn[0][0]*m_PolTrn[1][1] - m_PolTrn[1][0]*m_PolTrn[0][1]); //OC06092016
		m_InvPolTrn[0][0] = invDetPolTrn*m_PolTrn[1][1]; m_InvPolTrn[0][1] = -invDetPolTrn*m_PolTrn[0][1];
		m_InvPolTrn[1][0] = -invDetPolTrn*m_PolTrn[1][0]; m_InvPolTrn[1][1] = invDetPolTrn*m_PolTrn[0][0];

		// Calculate polarization vectors in the crystal frame for the diffracted 
		// beam, ignoring the dependence on (kx,ky):
		//... moved to RadPointModifier

		//Finding default reflected beam reference frame base vectors in the Lab frame
		//This should be re-calculated in a separate function before propagation, 
		//if the cooredinates of these vectors will be defines among the propagation parameters

/**
		//Moved to PropagateRadiation to be able to use avgPhotEn

		double lam0Br = 2.*m_dA*sin(thBrd); //check units!?
		//double milNDsq = sqrt(m_hMilND*m_hMilND + m_kMilND*m_kMilND + m_lMilND*m_lMilND);
		//double hn = m_HXAi[1] * milNDsq;
		//double ht = m_HXAi[2] * milNDsq;
		double hn = m_HXAi[1]; //OC180314: fix after conversation with A. Suvorov
		double ht = m_HXAi[2];

		double tz_p_lam0Br_hTau = tv[2] + lam0Br*ht;
		TVector3d vZ1p(sv[2], sqrt(1. - sv[2]*sv[2] - tz_p_lam0Br_hTau*tz_p_lam0Br_hTau), tz_p_lam0Br_hTau);
		vZ1p.Normalize();

		TVector3d vH(0, hn, ht);
		vH.Normalize();

		TVector3d vX1p = vH^vZ1p;
		double absX1p = vX1p.Abs();
		const double tolX1BackRefl = 1e-09; //to tune?
		if (absX1p > tolX1BackRefl) vX1p.Normalize();
		else { vX1p.x = 1.;  vX1p.y = 0.;  vX1p.z = 0.; }

		TVector3d vSt0_Rxl(sv[0], nv[0], tv[0]);
		TVector3d vSt1_Rxl(sv[1], nv[1], tv[1]);
		TVector3d vSt2_Rxl(sv[2], nv[2], tv[2]);
		TMatrix3d mRxl(vSt0_Rxl, vSt1_Rxl, vSt2_Rxl);

		//Default reflected beam reference frame base vectors in the Lab (i.e. input beam) frame:
		TVector3d vX1 = mRxl*vX1p, vZ1 = mRxl*vZ1p;
		vX1.Normalize(); vZ1.Normalize(); //just in case

		//Coordinates of the output Optical Axis vector in the frame of input beam
		//double vLxOut = vZ1.x, vLyOut = vZ1.y, vLzOut = vZ1.z; 
		//double vHxOut = vX1.x, vHyOut = vX1.y; //Coordinates of the Horizontal Base vector of the output frame

		FindRXtRef(vZ1, vX1, m_RLabXt, m_RXtRef);
		//FindRXtRef(vZ1, vX1, m_RLabXt, m_RXtRef, m_KLabRef);
**/
	}

	void FindDefOutFrameVect(double phEn, double hn, double ht, TVector3d& tv, TVector3d& sv, TVector3d& vZ1, TVector3d& vX1)
	{
		const double eAconv = 12398.4193009;
		double lam0Br = eAconv/phEn; //wavelength in [A]

		double tz_p_lam0Br_hTau = tv.z + lam0Br*ht; //tv[2] + lam0Br*ht;
		//TVector3d vZ1p(sv[2], sqrt(1. - sv[2]*sv[2] - tz_p_lam0Br_hTau*tz_p_lam0Br_hTau), tz_p_lam0Br_hTau);
		TVector3d vZ1p(sv.z, sqrt(1. - sv.z*sv.z - tz_p_lam0Br_hTau*tz_p_lam0Br_hTau), tz_p_lam0Br_hTau);
		vZ1p.Normalize();

		TVector3d vH(0, hn, ht);
		vH.Normalize();

		TVector3d vX1p = vH^vZ1p;
		double absX1p = vX1p.Abs();
		const double tolX1BackRefl = 1e-09; //to tune?
		if (absX1p > tolX1BackRefl) vX1p.Normalize();
		else { vX1p.x = 1.;  vX1p.y = 0.;  vX1p.z = 0.; }

		TVector3d nv = tv^sv;
		TVector3d vSt0_Rxl(sv.x, nv.x, tv.x);
		TVector3d vSt1_Rxl(sv.y, nv.y, tv.y);
		TVector3d vSt2_Rxl(sv.z, nv.z, tv.z);
		TMatrix3d mRxl(vSt0_Rxl, vSt1_Rxl, vSt2_Rxl);

		//Default reflected beam reference frame base vectors in the Lab (i.e. input beam) frame:
		vZ1 = mRxl*vZ1p; vZ1.Normalize(); //just in case
		
		//vX1 = mRxl*vX1p; vX1.Normalize(); 
		//OC150514
		//Selecting "Horizontal" and "Vertical" directions of the Output beam frame
		//trying to use "minimum deviation" from the corresponding "Horizontal" and "Vertical" directions of the Input beam frame
		TVector3d vEz0(0, 0, 1);
		TVector3d vE1 = vZ1^vEz0;

		double absE1x = ::fabs(vE1.x), absE1y = ::fabs(vE1.y);
		if(absE1x >= absE1y)
		{
			if(vE1.x > 0) vX1 = vE1;
			else vX1 = (-1)*vE1;
		}
		else
		{
			TVector3d vY1(0, 1, 0);
			if(vE1.y > 0) vY1 = vE1;
			else vY1 = (-1)*vE1;
			vX1 = vY1^vZ1;
		}
		vX1.Normalize(); 
	}

	void FindRXtRef(TVector3d& vZ1, TVector3d& vX1, double RLabXt[3][3], double resRXtRef[3][3])
	//void FindRXtRef(TVector3d& vZ1, TVector3d& vX1, double RLabXt[3][3], double resRXtRef[3][3], double resKLabRef[2][3])
	{//This is called from Ctor and/or eventually before Propagation (if the orientation of the Output Frame is defined among Propagation Parameters)
		//vZ1, vX1 are assumed to be normalized
		TVector3d vY1 = vZ1^vX1;

		// Input #6: reflected beam coordinate vectors in the Lab system
		double rx[] = { vX1.x, vX1.y, vX1.z };
		double ry[] = { vY1.x, vY1.y, vY1.z };
		double rz[] = { vZ1.x, vZ1.y, vZ1.z };

		// transform rx, ry, and rz into crystal coordinate system 
		double rxXt[] = {
			RLabXt[0][0] * rx[0] + RLabXt[0][1] * rx[1] + RLabXt[0][2] * rx[2],
			RLabXt[1][0] * rx[0] + RLabXt[1][1] * rx[1] + RLabXt[1][2] * rx[2],
			RLabXt[2][0] * rx[0] + RLabXt[2][1] * rx[1] + RLabXt[2][2] * rx[2]
		};
		double ryXt[] = {
			RLabXt[0][0] * ry[0] + RLabXt[0][1] * ry[1] + RLabXt[0][2] * ry[2],
			RLabXt[1][0] * ry[0] + RLabXt[1][1] * ry[1] + RLabXt[1][2] * ry[2],
			RLabXt[2][0] * ry[0] + RLabXt[2][1] * ry[1] + RLabXt[2][2] * ry[2]
		};
		double rzXt[] = {
			RLabXt[0][0] * rz[0] + RLabXt[0][1] * rz[1] + RLabXt[0][2] * rz[2],
			RLabXt[1][0] * rz[0] + RLabXt[1][1] * rz[1] + RLabXt[1][2] * rz[2],
			RLabXt[2][0] * rz[0] + RLabXt[2][1] * rz[1] + RLabXt[2][2] * rz[2]
		};

		// build RXtRef
		for(int i = 0; i < 3; i++)
		{
			resRXtRef[0][i] = rxXt[i];
			resRXtRef[1][i] = ryXt[i];
			resRXtRef[2][i] = rzXt[i];
		}
	}

	//void FindGridTransformMatrAngRepres(TVector3d& vtv, TVector3d& vsv, double HXAi[3], double RXtRef[3][3], srTSRWRadStructAccessData* pRad, double resKLabRef[2][3])
	//void FindAngMeshTransformMatr(TVector3d& vtv, TVector3d& vsv, double HXAi[3], double RXtRef[3][3], double phEn, double absMaxAng, double resKLabRef[2][3])
	void FindAngMeshTransformMatr(TVector3d& vtv, TVector3d& vsv, double HXAi[3], double RXtRef[3][3], double phEn, double resKLabRef[2][3])
	{//Transform Kin to Kout (Kout = resKLabRef*Kin)
	 //The output parameter is double resKLabRef[2][3]
	 //This transformation is photon energy dependent; and, strictly speaking, is different for each photon energy "layer";
	 //however, we assume that the photon energy interval is small, and the transformation can be calculated for an average photon energy.

		const double eAconv = 12398.4193009;
		double alamA = eAconv/phEn; //alamA = eAconv / EceV; //OC: wavelength in [A]?
		
		double tZ = vtv.z + alamA*HXAi[2]; //tv[2] + alamA*HXAi[2]; 
		double nY = sqrt(1. - vsv.z*vsv.z - tZ*tZ); //sqrt(1. - sv[2]*sv[2] - tZ*tZ); 

		double aX[] = {vsv.x, -(vsv.z*vsv.x + tZ*vtv.x)/nY, vtv.x};
		//aX[0] = sv[0]; 
		//aX[1] = -(sv[2]*sv[0] + tZ*tv[0])/nY;
		//aX[2] = tv[0]; 
		double aY[] = {vsv.y, -(vsv.z*vsv.y + tZ*vtv.y)/nY, vtv.y};
		//aY[0] = sv[1];
		//aY[1] = -(sv[2]*sv[1] + tZ*tv[1])/nY;
		//aY[2] = tv[1]; 
		double aZ[] = {vsv.z, nY, tZ};
		//aZ[0] = sv[2];
		//aZ[1] = nY;
		//aZ[2] = tZ;
		resKLabRef[0][0] = RXtRef[0][0]*aX[0] + RXtRef[0][1]*aX[1] + RXtRef[0][2]*aX[2]; 
		resKLabRef[1][0] = RXtRef[1][0]*aX[0] + RXtRef[1][1]*aX[1] + RXtRef[1][2]*aX[2]; 
		resKLabRef[0][1] = RXtRef[0][0]*aY[0] + RXtRef[0][1]*aY[1] + RXtRef[0][2]*aY[2]; 
		resKLabRef[1][1] = RXtRef[1][0]*aY[0] + RXtRef[1][1]*aY[1] + RXtRef[1][2]*aY[2]; 
		resKLabRef[0][2] = RXtRef[0][0]*aZ[0] + RXtRef[0][1]*aZ[1] + RXtRef[0][2]*aZ[2]; 
		resKLabRef[1][2] = RXtRef[1][0]*aZ[0] + RXtRef[1][1]*aZ[1] + RXtRef[1][2]*aZ[2]; 

/**
//OC130514: replaced by the block above, following suggestion of A.Suvorov

		double k0Ai = 1./alamA;
		double k0XAi[3], k0HAi[3], k1HAi[3]; //,KLabRef[2][3]; 

		//double k0xAi,k0yAi,k0zAi; 
		k0XAi[0] = vsv.z*k0Ai; //sv[2]*k0Ai; 
		k0XAi[2] = vtv.z*k0Ai + HXAi[2]; // tv[2]*k0Ai + HXAi[2]; 
		k0XAi[1] = sqrt(k0Ai*k0Ai - k0XAi[0]*k0XAi[0] - k0XAi[2]*k0XAi[2]); 

		for(int i = 0; i < 3; i++)
		{ 
			k0HAi[i] = RXtRef[i][0]*k0XAi[0] + RXtRef[i][1]*k0XAi[1] + RXtRef[i][2]*k0XAi[2]; 
		}

		resKLabRef[0][2] = k0HAi[0]/k0Ai; 
		resKLabRef[1][2] = k0HAi[1]/k0Ai; 

		//double k0xAi = 1.e-6*k0Ai; //????????????????????
		double k0xAi = absMaxAng*k0Ai; //????????????????????
		double k0zAi = sqrt(k0Ai*k0Ai - k0xAi*k0xAi); 

		k0XAi[0] = vsv.x*k0xAi + vsv.z*k0zAi; //sv[0]*k0xAi + sv[2]*k0zAi; 
		k0XAi[2] = vtv.x*k0xAi + vtv.z*k0zAi + HXAi[2]; //tv[0]*k0xAi + tv[2]*k0zAi + HXAi[2]; 
		k0XAi[1] = sqrt(k0Ai*k0Ai - k0XAi[0]*k0XAi[0] - k0XAi[2]*k0XAi[2]); 

		for(int i = 0; i < 3; i++)
		{
			k1HAi[i] = RXtRef[i][0]*k0XAi[0] + RXtRef[i][1]*k0XAi[1] + RXtRef[i][2]*k0XAi[2]; 
		}

		resKLabRef[0][0] = (k1HAi[0] - k0HAi[0])/k0xAi; 
		resKLabRef[1][0] = (k1HAi[1] - k0HAi[1])/k0xAi; 

		k0XAi[0] = vsv.y*k0xAi + vsv.z*k0zAi; //sv[1]*k0xAi + sv[2]*k0zAi; 
		k0XAi[2] = vtv.y*k0xAi + vtv.z*k0zAi + HXAi[2]; //tv[1]*k0xAi + tv[2]*k0zAi + HXAi[2]; 
		k0XAi[1] = sqrt(k0Ai*k0Ai - k0XAi[0]*k0XAi[0] - k0XAi[2]*k0XAi[2]); 

		for(int i = 0; i < 3; i++)
		{
			k1HAi[i] = RXtRef[i][0]*k0XAi[0] + RXtRef[i][1]*k0XAi[1] + RXtRef[i][2]*k0XAi[2];
		}

		resKLabRef[0][1] = (k1HAi[0] - k0HAi[0])/k0xAi;
		resKLabRef[1][1] = (k1HAi[1] - k0HAi[1])/k0xAi;
**/
	}

	int FindAngMeshTrf(srTSRWRadStructAccessData* pRad, srTOptCrystMeshTrf* pMeshTrf)
	{
		if((pRad == 0) || (pMeshTrf == 0)) return 0;
		if(pRad->avgPhotEn <= 0.) pRad->SetAvgPhotEnergyFromLimits();

		double xEndOld = (pRad->xStart) + (pRad->xStep)*(pRad->nx - 1);
		double zEndOld = (pRad->zStart) + (pRad->zStep)*(pRad->nz - 1);

		double absMaxAng = ::abs(pRad->xStart);
		double absMaxAngTest = ::abs(xEndOld);
		if(absMaxAng < absMaxAngTest) absMaxAng = absMaxAngTest;
		
		absMaxAngTest = ::abs(pRad->zStart);
		if(absMaxAng < absMaxAngTest) absMaxAng = absMaxAngTest;

		absMaxAngTest = ::abs(zEndOld);
		if(absMaxAng < absMaxAngTest) absMaxAng = absMaxAngTest;

		int nMesh = 1;
		if(pRad->ne > 1) nMesh = pRad->ne + 1;

		const double tolCrossTerm = 1.e-03; //1.e-04; //to steer
		const double tolFractStepMeshChange = 0.1; //to steer

		srTOptCrystMeshTrf *tMeshTrf = pMeshTrf;
		double phEn = pRad->avgPhotEn;
		for(int i=0; i<nMesh; i++)
		{
			if(i == 1) phEn = pRad->eStart;

			double waveLength_m = 1.23984193009e-06/phEn;
			//double absMaxAngRadians = absMaxAng*waveLength_m;

			//FindAngMeshTransformMatr(m_tv, m_sv, m_HXAi, m_RXtRef, phEn, absMaxAngRadians, tMeshTrf->matrKLabRef);
			FindAngMeshTransformMatr(m_tv, m_sv, m_HXAi, m_RXtRef, phEn, tMeshTrf->matrKLabRef); //OC130514

			double &a11 = tMeshTrf->matrKLabRef[0][0], &a12 = tMeshTrf->matrKLabRef[0][1], &a13 =  tMeshTrf->matrKLabRef[0][2];
			double &a21 = tMeshTrf->matrKLabRef[1][0], &a22 = tMeshTrf->matrKLabRef[1][1], &a23 =  tMeshTrf->matrKLabRef[1][2];

			tMeshTrf->crossTermsAreLarge = (::fabs(a12) > ::fabs(a11*tolCrossTerm)) || (::fabs(a21) > ::fabs(a22*tolCrossTerm));
			if(tMeshTrf->crossTermsAreLarge)
			{
				if((::fabs(a12*tolCrossTerm) > ::fabs(a11)) && (::fabs(a21*tolCrossTerm) > ::fabs(a22)))
				{//The case of rotation, i.e. the cross-terms strongly dominate over the diagonal terms; it also doesn't need the interpolations
					//tMeshTrf->dataTranspIsReq = true;
					tMeshTrf->crossTermsAreLarge = false;
				}
			}
			//tMeshTrf->crossTermsAreLarge = ((::fabs(a12) > ::fabs(a11*tolCrossTerm)) && (::fabs(a12*tolCrossTerm) < ::fabs(a11))) //OC130514 (i.e. a12 is comparable to a11)
			//							|| ((::fabs(a21) > ::fabs(a22*tolCrossTerm)) && (::fabs(a21*tolCrossTerm) < ::fabs(a22))); //(i.e. a21 is comparable to a22)

			double k0 = 1./waveLength_m;

			//tMeshTrf->xStart = (pRad->xStart)*a11 + k0*a13;
			//double xEndNew = xEndOld*a11 + k0*a13;
			tMeshTrf->xStart = (pRad->xStart)*a11 + (pRad->zStart)*a12 + k0*a13; //OC130514 (?)
			double xEndNew = xEndOld*a11 + zEndOld*a12 + k0*a13;
			tMeshTrf->xStep = 0; if(pRad->nx > 1) tMeshTrf->xStep = (xEndNew - (tMeshTrf->xStart))/(pRad->nx - 1);

			//OCTEST
			//tMeshTrf->xStart = (pRad->xStart)*a11; //+ k0*a13;
			//double xEndNew = xEndOld*a11; //+ k0*a13;
			//END OCTEST

			//tMeshTrf->zStart = (pRad->zStart)*a22 + k0*a23;
			//double zEndNew = zEndOld*a22 + k0*a23;
			tMeshTrf->zStart = (pRad->xStart)*a21 + (pRad->zStart)*a22 + k0*a23; //OC130514 (?)
			double zEndNew = xEndOld*a21 + zEndOld*a22 + k0*a23;
			tMeshTrf->zStep = 0; if(pRad->nz > 1) tMeshTrf->zStep = (zEndNew - (tMeshTrf->zStart))/(pRad->nz - 1);

			//OCTEST
			//tMeshTrf->zStart = (pRad->zStart)*a22; //+ k0*a23;
			//double zEndNew = zEndOld*a22; //+ k0*a23;
			//END OCTEST

			//tMeshTrf->dataTranspIsReq = (!tMeshTrf->crossTermsAreLarge) && (::fabs(a12) > ::fabs(a11*tolCrossTerm)) && (::fabs(a21) > ::fabs(a22*tolCrossTerm));

			double absTolFractStepMeshChangeX = ::fabs(pRad->xStep)*tolFractStepMeshChange;
			double absTolFractStepMeshChangeZ = ::fabs(pRad->zStep)*tolFractStepMeshChange;

			tMeshTrf->xMeshTrfIsReq = (::fabs((tMeshTrf->xStart) - (pRad->xStart)) > absTolFractStepMeshChangeX) || (::fabs(xEndNew - xEndOld) > absTolFractStepMeshChangeX);
			tMeshTrf->zMeshTrfIsReq = (::fabs((tMeshTrf->zStart) - (pRad->zStart)) > absTolFractStepMeshChangeZ) || (::fabs(zEndNew - zEndOld) > absTolFractStepMeshChangeZ);

			phEn += pRad->eStep;
			tMeshTrf++;
		}
		return 0;
	}

	int CorrectAngMesh(srTSRWRadStructAccessData* pRad, srTOptCrystMeshTrf* pMeshTrf, double*& ar_xStartValuesInSlicesE, double*& ar_zStartValuesInSlicesE)
	{//Check and, if necessary, apply corrections to the angular mesh to take into account eventual anamorphic magnification ("stretching", "rotation") occurring at asymmetric reflections
	 //Should be used before making FFT to Coordinate representation
		ar_xStartValuesInSlicesE = 0;
		ar_zStartValuesInSlicesE = 0;

		if((pRad == 0) || (pMeshTrf == 0)) return 0;

		if(pRad->ne == 1)
		{
			if(pMeshTrf->crossTermsAreLarge)
			{//do interpolation (flip is not required)
				//ddddddddddddddddddddddddddddddd
			}
			else
			{//simply change scale
				bool flipOverX = (pMeshTrf->xMeshTrfIsReq) && (pMeshTrf->xStep < 0);
				bool flipOverZ = (pMeshTrf->zMeshTrfIsReq) && (pMeshTrf->zStep < 0);
				if(flipOverX || flipOverZ) pRad->FlipFieldData(flipOverX, flipOverZ);

				if(pMeshTrf->xMeshTrfIsReq)
				{
					if(flipOverX)
					{
						double xEndNew = (pMeshTrf->xStart) + (pMeshTrf->xStep)*(pRad->nx - 1);
						pRad->xStart = xEndNew; pRad->xStep = -(pMeshTrf->xStep);
					}
					else
					{
						pRad->xStart = pMeshTrf->xStart; pRad->xStep = pMeshTrf->xStep;
					}
				}
				if(pMeshTrf->zMeshTrfIsReq)
				{
					if(flipOverZ)
					{
						double zEndNew = (pMeshTrf->zStart) + (pMeshTrf->zStep)*(pRad->nz - 1);
						pRad->zStart = zEndNew; pRad->zStep = -(pMeshTrf->zStep);
					}
					else
					{
						pRad->zStart = pMeshTrf->zStart; pRad->zStep = pMeshTrf->zStep;
					}
				}

				//if(pMeshTrf->dataTranspIsReq)
				//{
				//	pRad->TransposeFieldData();
				//}
			}
		}
		else
		{//In the case of multiple photon energies, each const-energy "slice" has to be brought to the "average" mesh by interpolation(?)
			//OC151014

			srTOptCrystMeshTrf *tMeshTrf = pMeshTrf + 1;
			bool crossTermsAreLargeTot = pMeshTrf->crossTermsAreLarge;
			bool flipOverXtot = (pMeshTrf->xMeshTrfIsReq) && (pMeshTrf->xStep < 0);
			bool flipOverZtot = (pMeshTrf->zMeshTrfIsReq) && (pMeshTrf->zStep < 0);

			const double relTolArg = 0.05; //0.01; //To steer

			//double xStartAvg = 0, xStepAvg = 0;
			//double zStartAvg = 0, zStepAvg = 0;
			bool xMeshTrfIsReqTot = pMeshTrf->xMeshTrfIsReq;
			bool zMeshTrfIsReqTot = pMeshTrf->zMeshTrfIsReq;
			bool xMeshTrfIsReqHolds = true, zMeshTrfIsReqHolds = true;
			bool xStartIsAvg = true, xStepIsAvg = true;
			bool zStartIsAvg = true, zStepIsAvg = true;

			for(long ie = 0; ie < pRad->ne; ie++)
			{
				crossTermsAreLargeTot |= tMeshTrf->crossTermsAreLarge;

				bool flipOverX = (tMeshTrf->xMeshTrfIsReq) && (tMeshTrf->xStep < 0);
				bool flipOverZ = (tMeshTrf->zMeshTrfIsReq) && (tMeshTrf->zStep < 0);
				flipOverXtot &= flipOverX;
				flipOverZtot &= flipOverZ;

				if(pMeshTrf->xMeshTrfIsReq != tMeshTrf->xMeshTrfIsReq) xMeshTrfIsReqHolds = true;
				if(pMeshTrf->zMeshTrfIsReq != tMeshTrf->zMeshTrfIsReq) zMeshTrfIsReqHolds = true;

				double inv_iep1 = 1./(ie + 1);
				//xStartAvg = (xStartAvg*ie + tMeshTrf->xStart)*inv_iep1;
				//xStepAvg = (xStepAvg*ie + tMeshTrf->xStep)*inv_iep1;
				//zStartAvg = (zStartAvg*ie + tMeshTrf->zStart)*inv_iep1;
				//zStepAvg = (zStepAvg*ie + tMeshTrf->zStep)*inv_iep1;

				xStartIsAvg &= (::fabs(pMeshTrf->xStart - tMeshTrf->xStart) <= ::fabs(pMeshTrf->xStart)*relTolArg);
				xStepIsAvg &= (::fabs(pMeshTrf->xStep - tMeshTrf->xStep) <= ::fabs(pMeshTrf->xStep)*relTolArg);
				zStartIsAvg &= (::fabs(pMeshTrf->zStart - tMeshTrf->zStart) <= ::fabs(pMeshTrf->zStart)*relTolArg);
				zStepIsAvg &= (::fabs(pMeshTrf->zStep - tMeshTrf->zStep) <= ::fabs(pMeshTrf->zStep)*relTolArg);

				tMeshTrf++;
			}

			//bool xStartIsAvg = (::fabs(pMeshTrf->xStart - xStartAvg) <= ::fabs(pMeshTrf->xStart)*relTolArg);
			//bool xStepIsAvg = (::fabs(pMeshTrf->xStep - xStepAvg) <= ::fabs(pMeshTrf->xStep)*relTolArg);
			//bool zStartIsAvg = (::fabs(pMeshTrf->zStart - zStartAvg) <= ::fabs(pMeshTrf->zStart)*relTolArg);
			//bool zStepIsAvg = (::fabs(pMeshTrf->zStep - zStepAvg) <= ::fabs(pMeshTrf->zStep)*relTolArg);

			//if((!crossTermsAreLargeTot) && xStartIsAvg && xStepIsAvg && zStartIsAvg && zStepIsAvg && xMeshTrfIsReqHolds && zMeshTrfIsReqHolds)
			if((!crossTermsAreLargeTot) && xStepIsAvg && zStepIsAvg && xMeshTrfIsReqHolds && zMeshTrfIsReqHolds)
			{
				if(flipOverXtot || flipOverZtot) pRad->FlipFieldData(flipOverXtot, flipOverZtot);

				if(pMeshTrf->xMeshTrfIsReq)
				{
					if(flipOverXtot)
					{
						double xEndNew = (pMeshTrf->xStart) + (pMeshTrf->xStep)*(pRad->nx - 1);
						pRad->xStart = xEndNew; pRad->xStep = -(pMeshTrf->xStep);
					}
					else
					{
						pRad->xStart = pMeshTrf->xStart; pRad->xStep = pMeshTrf->xStep;
					}
				}
				if(pMeshTrf->zMeshTrfIsReq)
				{
					if(flipOverZtot)
					{
						double zEndNew = (pMeshTrf->zStart) + (pMeshTrf->zStep)*(pRad->nz - 1);
						pRad->zStart = zEndNew; pRad->zStep = -(pMeshTrf->zStep);
					}
					else
					{
						pRad->zStart = pMeshTrf->zStart; pRad->zStep = pMeshTrf->zStep;
					}
				}

				//These two switches will allow for taking into account variable start values in slices at doing back FFT
				if((pMeshTrf->xMeshTrfIsReq) && (!xStartIsAvg)) 
				{
					ar_xStartValuesInSlicesE = new double[pRad->ne];
					double *t_ar_xStartValues = ar_xStartValuesInSlicesE;
					srTOptCrystMeshTrf *tMeshTrf = pMeshTrf + 1;
					long nx_mi_1 = pRad->nx - 1;
					for(long ie=0; ie<(pRad->ne); ie++)
					{
						if(flipOverXtot) *(t_ar_xStartValues++) = (tMeshTrf->xStart) + (tMeshTrf->xStep)*nx_mi_1;
						else *(t_ar_xStartValues++) = tMeshTrf->xStart;
						tMeshTrf++;
					}
				}
				if((pMeshTrf->zMeshTrfIsReq) && (!zStartIsAvg))
				{
					ar_zStartValuesInSlicesE = new double[pRad->ne];
					double *t_ar_zStartValues = ar_zStartValuesInSlicesE;
					srTOptCrystMeshTrf *tMeshTrf = pMeshTrf + 1;
					long nz_mi_1 = pRad->nz - 1;
					for(long ie=0; ie<(pRad->ne); ie++)
					{
						if(flipOverZtot) *(t_ar_zStartValues++) = (tMeshTrf->zStart) + (tMeshTrf->zStep)*nz_mi_1;
						else *(t_ar_zStartValues++) = tMeshTrf->zStart;
						tMeshTrf++;
					}
				}
			}
			else
			{//do interpolation of all slices to the Average mesh (flip is not required)
				int res = WfrInterpolOnRegGrid(pRad, pMeshTrf);
				if(res != 0) return res;
			}
		}
		return 0;
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect) //virtual in srTGenOptElem
	{
		m_eStartAux = pRadAccessData->eStart; m_eStepAux = pRadAccessData->eStep; m_ne = pRadAccessData->ne; //required for RadPointModifier

		if((ParPrecWfrPropag.vLxOut != 0) || (ParPrecWfrPropag.vLyOut != 0) || (ParPrecWfrPropag.vLzOut != 0)) 
		{//Process the case when the based vectors of the output beam frame are defined "manually"
			double vHz = FindLongCoordOfHorBaseVect(ParPrecWfrPropag.vLxOut, ParPrecWfrPropag.vLyOut, ParPrecWfrPropag.vLzOut, ParPrecWfrPropag.vHxOut, ParPrecWfrPropag.vHyOut);
			TVector3d vZ1(ParPrecWfrPropag.vLxOut, ParPrecWfrPropag.vLyOut, ParPrecWfrPropag.vLzOut);
			vZ1.Normalize();
			TVector3d vX1(ParPrecWfrPropag.vHxOut, ParPrecWfrPropag.vHyOut, vHz);
			vX1.Normalize();

			FindRXtRef(vZ1, vX1, m_RLabXt, m_RXtRef);
		}
		else
		{
			TVector3d vZ1, vX1;
			//double hn = m_HXAi[1]; //OC180314: fix after conversation with A. Suvorov
			//double ht = m_HXAi[2];
			FindDefOutFrameVect(pRadAccessData->avgPhotEn, m_HXAi[1], m_HXAi[2], m_tv, m_sv, vZ1, vX1);
			FindRXtRef(vZ1, vX1, m_RLabXt, m_RXtRef);
		}

		//return PropagateRadiationMeth_0(pRadAccessData);
		return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0);
	}

	//int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData) //virtual in srTGenOptElem
	//{//because for Crystal (the same way as for the Drift) the following works for many photon energies too!
	//	return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0);
	//}

	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{//It works for many photon energies too (as in the case of Drift)
	 //The "in-place" processing involving FFT for many photon energies greatly improves efficiency of the code for Time-/Frequency-Dependent simulations for FEL and pulsed lasers.
		int result;
		if (result = PropagateRadiationSimple_AngRepres(pRadAccessData)) return result;
		if (result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if (result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if (result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRadAccessData)
	{//Similar to the corresponding function for the Drift (?)
	 //It works for many photon energies too (as in the case of Drift).
	 //The "in-place" processing involving FFT for many photon energies greatly improves efficiency of the code for Time-/Frequency-Dependent simulations for FEL and pulsed lasers.
		int result;

		//double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
		//pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
		//pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
		//double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;
		//pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
		//pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;
		//Special processing of Shift was removed, because in the case of presence of anamorphic magnification it is more complicated than for drift(?)

		pRadAccessData->WfrEdgeCorrShouldBeDone = 0;

		if(pRadAccessData->Pres != 1)
		{//Switch to Angular representation
			if(result = SetRadRepres(pRadAccessData, 1)) return result;
		}

		srTOptCrystMeshTrf meshTrf;
		m_pMeshTrf = &meshTrf;
		int nMesh = 1;
		if(pRadAccessData->ne > 1) 
		{
			nMesh = pRadAccessData->ne + 1;
			m_pMeshTrf = new srTOptCrystMeshTrf[nMesh];
		}
		if(result = FindAngMeshTrf(pRadAccessData, m_pMeshTrf)) return result;

		if(result = TraverseRadZXE(pRadAccessData)) return result;

		//if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		//{
		//	pRadAccessData->xStartTr += xShift;
		//	pRadAccessData->zStartTr += zShift;
		//}

		double *ar_xStartInSlicesE=0,  *ar_zStartInSlicesE=0;
		if(result = CorrectAngMesh(pRadAccessData, m_pMeshTrf, ar_xStartInSlicesE,  ar_zStartInSlicesE)) return result;

		//Switch to Coordinate representation -- make optional?
		if(result = SetRadRepres(pRadAccessData, 0, ar_xStartInSlicesE, ar_zStartInSlicesE)) return result;

		//pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;
		//if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		//{
		//	pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
		//	pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
		//}

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();

		if(ar_xStartInSlicesE != 0) delete[] ar_xStartInSlicesE;
		if(ar_zStartInSlicesE != 0) delete[] ar_zStartInSlicesE;

		if(nMesh > 1) delete[] m_pMeshTrf;
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	//void RadPointModifier_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// E in eV; Length in m !!!
		// Operates on Angles side !!!
		// EXZ.x = angX_rad/wavelength_m, EXZ.z = angY_rad/wavelength_m

		int ie = 0;
		if(m_ne > 1) ie = int((EXZ.e - m_eStartAux)/m_eStepAux + 1e-06) + 1;
		srTOptCrystMeshTrf &meshTrf = m_pMeshTrf[ie];
		double &a11 = (meshTrf.matrKLabRef)[0][0], &a12 = (meshTrf.matrKLabRef)[0][1];
		double &a21 = (meshTrf.matrKLabRef)[1][0], &a22 = (meshTrf.matrKLabRef)[1][1];
		double asymFact = ::fabs(a11*a22 - a12*a21);
		if(asymFact > 0) asymFact = 1./asymFact;
		else asymFact = 1.;

		// Conversion factor: wavelength [A] <-> energy [eV] ?
		const double eAconv = 12398.4193009;
		const double pi = 4.*atan(1.);
		const complex<double> iC(0., 1.);

		const double DBMaxFact = 0.01; //0.1; //to tune
		const double logDBMax = log(DBMaxFact*(numeric_limits<double>::max()));

		double alamA = eAconv/EXZ.e; //alamA = eAconv / EceV; //OC: wavelength in [A]?
		double k0Ai = 1./alamA;

		// Wave vector components: hard-coded for this test only! 
		//kxAi = 0.0; 
		//kyAi = (-50.e-6 + 100.e-6*(double)(ik)/(double)(nk-1))/alamA; 

		double angXrad = EXZ.x*alamA*1.e-10, angYrad = EXZ.z*alamA*1.e-10;
		//double kxAi = angXrad/alamA, kyAi = angYrad/alamA;
		double kxAi = EXZ.x*1.e-10; //OC: ??
		double kyAi = EXZ.z*1.e-10;

		// Electric field amplitudes in k-space: hard-coded for this test only! 
		//Eox = complex<double>(1.,0.); 
		//Eoy = complex<double>(1.,0.); 
		complex<double> Eox = complex<double>(*(EPtrs.pExRe), *(EPtrs.pExIm));
		complex<double> Eoy = complex<double>(*(EPtrs.pEzRe), *(EPtrs.pEzIm));

		//OCTEST111014
		//double testI0 = Eox.real()*Eox.real() + Eox.imag()*Eox.imag() + Eoy.real()*Eoy.real() + Eoy.imag()*Eoy.imag();

		// k0wvAi = incident beam wave vector (kx,ky,kz). 
		double kzAi = sqrt(k0Ai*k0Ai - kxAi*kxAi - kyAi*kyAi); //OC: possible loss of precision?

		//k0wAi[0] = kxAi; 
		//k0wAi[1] = kyAi; 
		//k0wAi[2] = kzAi; 
		double k0wAi[] = { kxAi, kyAi, kzAi };

		//printf("Incident wave vector components in Lab frame: \n"); 
		//printf("%e %e %e\n",k0wAi[0],k0wAi[1],k0wAi[2]); 
		//u0[0] = kxAi/k0Ai; 
		//u0[1] = kyAi/k0Ai; 
		//u0[2] = kzAi/k0Ai;
		double u0[] = { kxAi / k0Ai, kyAi / k0Ai, kzAi / k0Ai };

		// Conversion of wave vector components to crystal frame: 
		//k0wXAi[0] = RLabXt[0][0]*k0wAi[0] + RLabXt[0][1]*k0wAi[1] + RLabXt[0][2]*k0wAi[2]; 
		//k0wXAi[1] = RLabXt[1][0]*k0wAi[0] + RLabXt[1][1]*k0wAi[1] + RLabXt[1][2]*k0wAi[2]; 
		//k0wXAi[2] = RLabXt[2][0]*k0wAi[0] + RLabXt[2][1]*k0wAi[1] + RLabXt[2][2]*k0wAi[2]; 
		double k0wXAi[] = {
			m_RLabXt[0][0] * k0wAi[0] + m_RLabXt[0][1] * k0wAi[1] + m_RLabXt[0][2] * k0wAi[2],
			m_RLabXt[1][0] * k0wAi[0] + m_RLabXt[1][1] * k0wAi[1] + m_RLabXt[1][2] * k0wAi[2],
			m_RLabXt[2][0] * k0wAi[0] + m_RLabXt[2][1] * k0wAi[1] + m_RLabXt[2][2] * k0wAi[2]
		};

		//printf("Incident wave vector components in Crystal frame: \n"); 
		//printf("%e %e %e\n",k0wXAi[0],k0wXAi[1],k0wXAi[2]); 
		//u0X[0] = k0wXAi[0]/k0Ai; 
		//u0X[1] = k0wXAi[1]/k0Ai; 
		//u0X[2] = k0wXAi[2]/k0Ai; 
		double u0X[] = { k0wXAi[0] / k0Ai, k0wXAi[1] / k0Ai, k0wXAi[2] / k0Ai };

		// Convert components of incident beam polarization from (e1X,e2X) to (sg0X,pi0X): 
		//EInSPs = PolTrn[0][0]*Eox + PolTrn[0][1]*Eoy; 
		//EInSPp = PolTrn[1][0]*Eox + PolTrn[1][1]*Eoy; 
		complex<double> EInSPs = m_PolTrn[0][0] * Eox + m_PolTrn[0][1] * Eoy;
		complex<double> EInSPp = m_PolTrn[1][0] * Eox + m_PolTrn[1][1] * Eoy;

		// Calculate the wave vector for the diffracted beam. Note that here the index of 
		// refraction corrections are included. 
		//kHtXAi[0] = k0wXAi[0] + HXAi[0]; 
		//kHtXAi[1] = 0.; 
		//kHtXAi[2] = k0wXAi[2] + HXAi[2]; 
		double kHtXAi[] = { k0wXAi[0] + m_HXAi[0], 0., k0wXAi[2] + m_HXAi[2] };
		double kHt2 = kHtXAi[0] * kHtXAi[0] + kHtXAi[2] * kHtXAi[2];

		//kmHXAi[0] = kHtXAi[0]; 
		//kmHXAi[1] = sqrt(k0Ai*k0Ai-kHt2); 
		//kmHXAi[2] = kHtXAi[2]; 
		double kmHXAi[] = { kHtXAi[0], sqrt(k0Ai*k0Ai - kHt2), kHtXAi[2] }; //OC: can it happen (k0Ai*k0Ai - kHt2) < 0?

		// Calculate direction cosine gamma0, reflection asymmetry parameter bee, 
		// deviation parameter Adev and normalized deviation parameter zeeC (as 
		// defined by Zachariasen gamma0, b, alpha and z.) 
		double gamma0 = -u0X[1];
		double bee = 1. / (1. + m_HXAi[1] / k0wXAi[1]);
		double dotkH = k0wXAi[0] * m_HXAi[0] + k0wXAi[1] * m_HXAi[1] + k0wXAi[2] * m_HXAi[2];
		double Adev = (2.*dotkH + 1. / (m_dA*m_dA)) / (k0Ai*k0Ai);
		complex<double> zeeC = 0.5*((1. - bee)*m_psi0c + bee*Adev);

		// Calculate the complex reflectivity DHsgC for sigma polarization: 
		complex<double> queC = bee*m_psihc*m_psimhc;
		complex<double> sqrqzC = sqrt(queC + zeeC*zeeC);
		complex<double> del1C = 0.5*(m_psi0c - zeeC + sqrqzC);
		complex<double> del2C = 0.5*(m_psi0c - zeeC - sqrqzC);
		complex<double> x1C = (-zeeC + sqrqzC) / m_psimhc;
		complex<double> x2C = (-zeeC - sqrqzC) / m_psimhc;
		//complex<double> ph1C = 2.*pi*iC*k0Ai*del1C*(m_thicum*1.E+04)/gamma0;
		//complex<double> ph2C = 2.*pi*iC*k0Ai*del2C*(m_thicum*1.E+04)/gamma0;
		complex<double> aux_phC = 2.*pi*iC*k0Ai*(m_thicum*1.E+04) / gamma0;
		complex<double> ph1C = aux_phC*del1C;
		complex<double> ph2C = aux_phC*del2C;

		complex<double> DHsgC, Cph1C, Cph2C, D0trsC;
		if(real(ph1C) > logDBMax) DHsgC = x2C;
		else if(real(ph2C) > logDBMax) DHsgC = x1C;
		else
		{
			Cph1C = exp(ph1C);
			Cph2C = exp(ph2C);
			DHsgC = x1C*x2C*(Cph2C - Cph1C)/(Cph2C*x2C - Cph1C*x1C);
		}

		//double re_ph1C = real(ph1C), re_ph2C = real(ph2C);
		//if(re_ph1C > logDBMax) DHsgC = x2C;
		//else if(re_ph2C > logDBMax) DHsgC = x1C;
		//else
		//{
		//	if(re_ph1C < -logDBMax) 
		//	{
		//		Cph1C = complex<double>(0, 0);
		//	}
		//	else Cph1C = exp(ph1C);

		//	if(re_ph2C < -logDBMax) 
		//	{
		//		Cph2C = complex<double>(0, 0);
		//	}
		//	else Cph2C = exp(ph2C);

		//	if((re_ph1C < -logDBMax) && (re_ph2C < -logDBMax)) 
		//	{
		//		DHsgC = complex<double>(0, 0); //?
		//	}
		//	else DHsgC = x1C*x2C*(Cph2C - Cph1C)/(Cph2C*x2C - Cph1C*x1C);
		//}

		if(m_itrans == 1)
		{
			// calculate the complex reflectivity D0trsC of the transmitted beam. 
			if(real(ph1C) > logDBMax)
			{
				Cph2C = exp(ph2C);
				D0trsC = -Cph2C*(x2C - x1C)/x1C;
			}
			else if(real(ph2C) > logDBMax)
			{
				Cph1C = exp(ph1C);
				D0trsC = +Cph1C*(x2C - x1C)/x2C;
			}
			else D0trsC = Cph1C*Cph2C*(x2C - x1C)/(Cph2C*x2C - Cph1C*x1C);
		}

		// Calculate the complex reflectivity DHpiC for pi polarization: 
		//cos2t = cos(2.*thBrd); //OC: moved to members m_cos2t
		queC = bee*m_psihc*m_psimhc*m_cos2t*m_cos2t;
		sqrqzC = sqrt(queC + zeeC*zeeC);
		del1C = 0.5*(m_psi0c - zeeC + sqrqzC);
		del2C = 0.5*(m_psi0c - zeeC - sqrqzC);
		x1C = (-zeeC + sqrqzC) / (m_psimhc*m_cos2t);
		x2C = (-zeeC - sqrqzC) / (m_psimhc*m_cos2t);
		ph1C = 2.*pi*iC*k0Ai*del1C*(m_thicum*1.E+04) / gamma0;
		ph2C = 2.*pi*iC*k0Ai*del2C*(m_thicum*1.E+04) / gamma0;

		complex<double> DHpiC, D0trpC;
		if(real(ph1C) > logDBMax) DHpiC = x2C;
		else if(real(ph2C) > logDBMax) DHpiC = x1C;
		else
		{
			Cph1C = exp(ph1C);
			Cph2C = exp(ph2C);
			DHpiC = x1C*x2C*(Cph2C - Cph1C)/(Cph2C*x2C - Cph1C*x1C);
		}

		if(m_itrans == 1)
		{
			// calculate the complex reflectivity D0trpC of the transmitted beam. 
			if (real(ph1C) > logDBMax)
			{
				Cph2C = exp(ph2C);
				D0trpC = -Cph2C*(x2C - x1C)/x1C;
			}
			else if (real(ph2C) > logDBMax)
			{
				Cph1C = exp(ph1C);
				D0trpC = +Cph1C*(x2C - x1C)/x2C;
			}
			else D0trpC = Cph1C*Cph2C*(x2C - x1C)/(Cph2C*x2C - Cph1C*x1C);
		}

		// Calculate the diffracted amplitudes: 
		complex<double> EHSPCs = DHsgC*EInSPs;
		complex<double> EHSPCp = DHpiC*EInSPp;
		complex<double> E0tSPs, E0tSPp;
		if(m_itrans == 1)
		{
			E0tSPs = D0trsC*EInSPs;
			E0tSPp = D0trpC*EInSPp;
		}

		//kc0XAi[0] = k0Ai*RLabXt[0][2]; 
		//kc0XAi[1] = k0Ai*RLabXt[1][2]; 
		//kc0XAi[2] = k0Ai*RLabXt[2][2];
		double kc0XAi[] = { k0Ai*m_RLabXt[0][2], k0Ai*m_RLabXt[1][2], k0Ai*m_RLabXt[2][2] };

		// Calculate polarization vectors in the crystal frame for the diffracted 
		// beam, ignoring the dependence on (kx,ky): 
		//kcHtXA[0] = kc0XAi[0] + HXAi[0];
		//kcHtXA[1] = 0.;
		//kcHtXA[2] = kc0XAi[2] + HXAi[2];
		double kcHtXA[] = { kc0XAi[0] + m_HXAi[0], 0, kc0XAi[2] + m_HXAi[2] };
		double kcHt2 = kcHtXA[0] * kcHtXA[0] + kcHtXA[2] * kcHtXA[2];
		//kcHXAi[0] = kcHtXA[0];
		//kcHXAi[1] = sqrt(k0Ai*k0Ai - kcHt2);
		//kcHXAi[2] = kcHtXA[2];
		double kcHXAi[] = { kcHtXA[0], sqrt(k0Ai*k0Ai - kcHt2), kcHtXA[2] }; //to check if (k0Ai*k0Ai - kcHt2) can be < 0
		//ucnHX[0] = kcHXAi[0] / k0Ai;
		//ucnHX[1] = kcHXAi[1] / k0Ai;
		//ucnHX[2] = kcHXAi[2] / k0Ai;
		double ucnHX[] = { kcHXAi[0]/k0Ai, kcHXAi[1]/k0Ai, kcHXAi[2]/k0Ai };
		double sgHX[] = { m_sg0X[0], m_sg0X[1], m_sg0X[2] };
		//piHX[0] = ucnHX[1]*sgHX[2] - ucnHX[2]*sgHX[1]; 
		//piHX[1] = ucnHX[2]*sgHX[0] - ucnHX[0]*sgHX[2]; 
		//piHX[2] = ucnHX[0]*sgHX[1] - ucnHX[1]*sgHX[0]; 
		double piHX[] = {
			ucnHX[1] * sgHX[2] - ucnHX[2] * sgHX[1],
			ucnHX[2] * sgHX[0] - ucnHX[0] * sgHX[2],
			ucnHX[0] * sgHX[1] - ucnHX[1] * sgHX[0]
		};

		// Transform the diffracted beam polarization vectors into the diffracted beam frame: 
		double sgH[3], piH[3];
		for(int i = 0; i < 3; i++)
		{
			sgH[i] = m_RXtRef[i][0]*sgHX[0] + m_RXtRef[i][1]*sgHX[1] + m_RXtRef[i][2]*sgHX[2];
			piH[i] = m_RXtRef[i][0]*piHX[0] + m_RXtRef[i][1]*piHX[1] + m_RXtRef[i][2]*piHX[2];
		}

		// Calculate diffracted amplitudes in the diffracted beam frame
		//complex<double> Ehx = sgH[0]*EHSPCs + piH[0]*EHSPCp;
		//complex<double> Ehy = sgH[1]*EHSPCs + piH[1]*EHSPCp;
		//complex<double> Ehz = sgH[2]*EHSPCs + piH[2]*EHSPCp; //Longitudinal component not used (?)

		// Transverse components of the output electric field in the frame of the output beam:
		//*(EPtrs.pExRe) = (float)(asymFact*Ehx.real());
		//*(EPtrs.pExIm) = (float)(asymFact*Ehx.imag());
		//*(EPtrs.pEzRe) = (float)(asymFact*Ehy.real());
		//*(EPtrs.pEzIm) = (float)(asymFact*Ehy.imag());

		//OC05092016
		if(m_uc == 1) //Bragg Reflection
		{
			complex<double> Ehx = sgH[0]*EHSPCs + piH[0]*EHSPCp;
			complex<double> Ehy = sgH[1]*EHSPCs + piH[1]*EHSPCp;
			complex<double> Ehz = sgH[2]*EHSPCs + piH[2]*EHSPCp; //Longitudinal component not used (?)
		
			//Transverse components of the output electric field in the frame of the output beam:
			*(EPtrs.pExRe) = (float)(asymFact*Ehx.real());
			*(EPtrs.pExIm) = (float)(asymFact*Ehx.imag());
			*(EPtrs.pEzRe) = (float)(asymFact*Ehy.real());
			*(EPtrs.pEzIm) = (float)(asymFact*Ehy.imag());
		}
		else if(m_uc == 2) //Bragg Transmission (i.e. Forward Bragg Diffraction (FBD))
		{
			//DEBUG
			//complex<double> DHsgC_p_D0trsC = DHsgC + D0trsC;
			//complex<double> DHsgCae2_p_D0trsCae2 = DHsgC.real()*DHsgC.real() + DHsgC.imag()*DHsgC.imag() + D0trsC.real()*D0trsC.real() + D0trsC.imag()*D0trsC.imag();
			//END DEBUG

			complex<double> Etx = m_InvPolTrn[0][0]*E0tSPs + m_InvPolTrn[0][1]*E0tSPp;
			complex<double> Ety = m_InvPolTrn[1][0]*E0tSPs + m_InvPolTrn[1][1]*E0tSPp;
		
			//Transverse components of the output electric field in the frame of the output beam:
			*(EPtrs.pExRe) = (float)(Etx.real());
			*(EPtrs.pExIm) = (float)(Etx.imag());
			*(EPtrs.pEzRe) = (float)(Ety.real());
			*(EPtrs.pEzIm) = (float)(Ety.imag());
		}

		//OCTEST111014
		//double testI1 = (*(EPtrs.pExRe))*(*(EPtrs.pExRe)) + (*(EPtrs.pExIm))*(*(EPtrs.pExIm)) + 
		//				(*(EPtrs.pEzRe))*(*(EPtrs.pEzRe)) + (*(EPtrs.pEzIm))*(*(EPtrs.pEzIm));

		//complex<double> Etx, Ety, Etz;
		//if(m_itrans == 1)
		//{
		//	Etx = sgH[0]*E0tSPs + piH[0]*E0tSPp;
		//	Ety = sgH[1]*E0tSPs + piH[1]*E0tSPp;
		//	Etz = sgH[2]*E0tSPs + piH[2]*E0tSPp;
		//}

		// Convert the diffracted beam's wave vector from the crystal frame to the diffracted beam frame.
		//kmHAi[0] = RXtRef[0][0]*kmHXAi[0] + RXtRef[0][1]*kmHXAi[1] + RXtRef[0][2]*kmHXAi[2];
		//kmHAi[1] = RXtRef[1][0]*kmHXAi[0] + RXtRef[1][1]*kmHXAi[1] + RXtRef[1][2]*kmHXAi[2];
		//kmHAi[2] = RXtRef[2][0]*kmHXAi[0] + RXtRef[2][1]*kmHXAi[1] + RXtRef[2][2]*kmHXAi[2];
		//double kmHAi[] = {
		//	m_RXtRef[0][0]*kmHXAi[0] + m_RXtRef[0][1]*kmHXAi[1] + m_RXtRef[0][2]*kmHXAi[2],
		//	m_RXtRef[1][0]*kmHXAi[0] + m_RXtRef[1][1]*kmHXAi[1] + m_RXtRef[1][2]*kmHXAi[2],
		//	m_RXtRef[2][0]*kmHXAi[0] + m_RXtRef[2][1]*kmHXAi[1] + m_RXtRef[2][2]*kmHXAi[2]
		//};

		/**
		//printf("Reflected wave vector components in Lab frame: \n");
		//printf("%e %e %e\n",kmHAi[0],kmHAi[1],kmHAi[2]);
		//
		// For test output only:
		dangur = -asin(kyAi/k0Ai)*1.E+06;
		//uref[0] = 0.;
		//uref[1] = sin(2.*thBrd);
		//uref[2] = cos(2.*thBrd);
		//urkmH = uref[0]*kmHAi[0] + uref[1]*kmHAi[1] + uref[2]*kmHAi[2];
		//if (urkmH > k0Ai)
		//kangur = 0.;
		//else
		//kangur = acos(urkmH/k0Ai)*1.E+06;
		//ukcrAi[0] = uref[1]*kmHAi[2] - uref[2]*kmHAi[1];
		//ukcrAi[1] = uref[2]*kmHAi[0] - uref[0]*kmHAi[2];
		//ukcrAi[2] = uref[0]*kmHAi[1] - uref[1]*kmHAi[0];
		//if (ukcrAi[0] > 0.)
		//kangur = -kangur;
		//mgkmk0 = sqrt( kmHAi[0]*kmHAi[0] + kmHAi[1]*kmHAi[1] + kmHAi[2]*kmHAi[2] )/k0Ai;
		//mgDHsg = abs(DHsgC);
		//mgDHpi = abs(DHpiC);
		//mgEHs = abs(EHSPCs);
		//mgEHp = abs(EHSPCp);
		//fprintf(TestK,"%e %e %e\n",dangur,mgkmk0,kangur);
		//fprintf(TestRef,"%e %e %e %e %e\n",dangur,mgDHsg*mgDHsg,arg(DHsgC),mgDHpi*mgDHpi,arg(DHpiC));
		fprintf(TestAmpH,"%e %e %e %e %e %e %e\n",dangur,real(Ehx),imag(Ehx),real(Ehy),imag(Ehy),1.e6*kyAi,1.e6*kmHAi[1]);
		//
		if (itrans == 1){
		//mgE0ts = abs(E0tSPs);
		//mgE0tp = abs(E0tSPp);
		fprintf(TestAmpO,"%e %e %e %e %e\n",dangur,real(E0tSPs),imag(E0tSPs),real(E0tSPp),imag(E0tSPp));
		//fprintf(TestTrns,"%e %e %e %e %e\n",dangur,mgE0ts*mgE0ts,arg(E0tSPs),mgE0tp*mgE0tp,arg(E0tSPp));
		}
		**/
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRadAccessData, srTMomentsRatios* MomRatArray)
	{//Do nothing for Crystal(?), or recalculate, because angular divergence may be changed by Crystal
		//double aStr0[] = { 1., Length };
		//double aStr1[] = { 0., 1. };
		//double* a[] = { aStr0, aStr1 };
		//return AuxPropagateRadMoments(pRadAccessData, a, a, MomRatArray);
		return 0;
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{//Do nothing for Crystal(?)
		// RobsX, RobsZ are treated here as distances to waists 
		return 0;
	}

	int Propagate4x4PropMatr(srTSRWRadStructAccessData* pRadAccessData)
	{//Do nothing for Crystal(?)
		//double Drift4x4Matr[] = { 1., Length, 0., 0.,
		//						  0., 1., 0., 0.,
		//						  0., 0., 1., Length,
		//						  0., 0., 0., 1.};
		//double Drift4Vect[] = { 0., 0., 0., 0.};
		//return GenAuxPropagate4x4PropMatr(pRadAccessData, Drift4x4Matr, Drift4Vect);
		return 0;
	}

	void CollectArgStartValuesInSlices(srTOptCrystMeshTrf* arMeshTrf, long ne, char h_or_v, double*& arStartValuesInSlicesE)
	{
		arStartValuesInSlicesE = 0;
		if((arMeshTrf == 0) || (ne <= 0)) return;

		arStartValuesInSlicesE = new double[ne];
		double *t_arStartValuesInSlicesE = arStartValuesInSlicesE;
		srTOptCrystMeshTrf *t_arMeshTrf = arMeshTrf;

		if((h_or_v == 'h') || (h_or_v == 'H') || (h_or_v == 'x') || (h_or_v == 'X'))
		{
			for(long ie=0; ie<ne; ie++) *(t_arStartValuesInSlicesE++) = (t_arMeshTrf++)->xStart;
		}
		else if((h_or_v == 'v') || (h_or_v == 'V') || (h_or_v == 'z') || (h_or_v == 'Z'))
		{
			for(long ie=0; ie<ne; ie++) *(t_arStartValuesInSlicesE++) = (t_arMeshTrf++)->zStart;
		}
		else
		{
			for(long ie=0; ie<ne; ie++) *(t_arStartValuesInSlicesE++) = 0.;
		}
	}

	int WfrInterpolOnRegGrid(srTSRWRadStructAccessData* pRad, srTOptCrystMeshTrf* pMeshTrf);
};

//*************************************************************************

#endif
