/************************************************************************//**
 * File: sroptcryst.h
 * Description: Optical element: Angle (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2013
 *
 * Copyright (C) Brookhaven National Laboratory, Upton, NY, USA
 * Copyright (C) Diamond Light Source, UK
 * All Rights Reserved
 *
 * @author J.Sutter, A. Suvorov, O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTCRYST_H
#define __SROPTCRYST_H

#include "sroptelm.h"
#include "srwlib.h"

//*************************************************************************

class srTOptCryst : public srTGenOptElem {

	double m_dA; /* crystal reflecting planes d-spacing (units?) */
	//double m_psi0r, m_psi0i; /* real and imaginary parts of 0-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psi0c; 
	//double m_psiHr, m_psiHi; /* real and imaginary parts of H-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psihc;
	//double m_psiHbr, m_psiHbi; /* real and imaginary parts of -H-th Fourier component of crystal polarizability (units?) */
	complex<double> m_psimhc;
	double m_hMilND, m_kMilND, m_lMilND; /* 1st, 2nd and 3rd  indexes of diffraction vector (Miller indices) */
	double m_thicum; /* crystal thickness [microns] */
	double m_alphrd; /* asymmetry angle [rad] */

	TVector3d m_nv; /* horizontal, vertical and longitudinal coordinates of outward normal to crystal surface in the frame of incident beam */
	TVector3d m_tv; /* horizontal, vertical  and vertical coordinates of central tangential vector [m] in the frame of incident beam */

public:

	srTOptCryst(const SRWLOptCryst& srwlCr)
	{
		m_dA = srwlCr.dSp;
		//m_psi0r = srwlCr.psi0r;
		//m_psi0i = srwlCr.psi0i;
		m_psi0c = complex<double>(srwlCr.psi0r, srwlCr.psi0i);
		//m_psiHr = srwlCr.psiHr;
		//m_psiHi = srwlCr.psiHi;
		m_psihc = complex<double>(srwlCr.psiHr, srwlCr.psiHi);
		//m_psiHbr = srwlCr.psiHbr;
		//m_psiHbi = srwlCr.psiHbi;
		m_psimhc = complex<double>(srwlCr.psiHbr, srwlCr.psiHbi);

		m_hMilND = srwlCr.h1;
		m_kMilND = srwlCr.h2;
		m_lMilND = srwlCr.h3;

		//m_tc = srwlCr.tc;
		m_thicum = srwlCr.tc*1e+06; //assuming srwlCr.tc in [m]
		m_alphrd = srwlCr.angAs;

		m_nv.x = srwlCr.nvx;
		m_nv.y = srwlCr.nvy;
		m_nv.z = srwlCr.nvz;
		m_nv.Normalize();

		m_tv.x = srwlCr.tvx;
		m_tv.y = srwlCr.tvy;
		m_tv.z = (-m_nv.x*m_tv.x - m_nv.y*m_tv.y)/m_nv.z;
		m_tv.Normalize();
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect) //virtual
	{
		return PropagateRadiationMeth_0(pRadAccessData);
	}

	int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData) //virtual in srTGenOptElem
	{//because for Crystal (the same way as for the Drift) the following works for many photon energies too!
		return PropagateRadiationSingleE_Meth_0(pRadAccessData, 0);
	}

	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{//it works for many photon energies too!
		int result; 
		if(result = PropagateRadiationSimple_AngRepres(pRadAccessData)) return result;
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRadAccessData)
	{//same as for Drift (?)
		int result;

		double xStartOld = pRadAccessData->xStart, zStartOld = pRadAccessData->zStart;
		pRadAccessData->xStart = -(pRadAccessData->nx >> 1)*pRadAccessData->xStep;
		pRadAccessData->zStart = -(pRadAccessData->nz >> 1)*pRadAccessData->zStep;
		double xShift = pRadAccessData->xStart - xStartOld, zShift = pRadAccessData->zStart - zStartOld;

		pRadAccessData->xWfrMin += xShift; pRadAccessData->xWfrMax += xShift;
		pRadAccessData->zWfrMin += zShift; pRadAccessData->zWfrMax += zShift;

		pRadAccessData->WfrEdgeCorrShouldBeDone = 0;

		if(pRadAccessData->Pres != 1) 
		{
			if(result = SetRadRepres(pRadAccessData, 1)) return result;
		}

		if(result = TraverseRadZXE(pRadAccessData)) return result;

		if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		{
			pRadAccessData->xStartTr += xShift;
			pRadAccessData->zStartTr += zShift;
		}

		if(result = SetRadRepres(pRadAccessData, 0)) return result;

		pRadAccessData->xStart = xStartOld; pRadAccessData->zStart = zStartOld;

		if(pRadAccessData->UseStartTrToShiftAtChangingRepresToCoord)
		{
			pRadAccessData->xStart = pRadAccessData->xStartTr - xShift;
			pRadAccessData->zStart = pRadAccessData->zStartTr - zShift;
		}

		pRadAccessData->SetNonZeroWavefrontLimitsToFullRange();
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	//void RadPointModifier_AngRepres(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Angles side !!!


/**
		// Wave vector components: hard-coded for this test only! 
		kxAi = 0.0; 
		kyAi = (-50.e-6 + 100.e-6*(double)(ik)/(double)(nk-1))/alamA; 
		// Electric field amplitudes in k-space: hard-coded for this test only! 
		Eox = complex<double>(1.,0.); 
		Eoy = complex<double>(1.,0.); 
		// k0wvAi = incident beam wave vector (kx,ky,kz). 
		kzAi = sqrt(k0Ai*k0Ai - kxAi*kxAi - kyAi*kyAi); 
		k0wAi[0] = kxAi; 
		k0wAi[1] = kyAi; 
		k0wAi[2] = kzAi; 
		//
		//printf("Incident wave vector components in Lab frame: \n"); 
		//printf("%e %e %e\n",k0wAi[0],k0wAi[1],k0wAi[2]); 
		//
		u0[0] = kxAi/k0Ai; 
		u0[1] = kyAi/k0Ai; 
		u0[2] = kzAi/k0Ai; 
		// Conversion of wave vector components to crystal frame: 
		k0wXAi[0] = RLabXt[0][0]*k0wAi[0] + RLabXt[0][1]*k0wAi[1] + RLabXt[0][2]*k0wAi[2]; 
		k0wXAi[1] = RLabXt[1][0]*k0wAi[0] + RLabXt[1][1]*k0wAi[1] + RLabXt[1][2]*k0wAi[2]; 
		k0wXAi[2] = RLabXt[2][0]*k0wAi[0] + RLabXt[2][1]*k0wAi[1] + RLabXt[2][2]*k0wAi[2]; 
		//
		//printf("Incident wave vector components in Crystal frame: \n"); 
		//printf("%e %e %e\n",k0wXAi[0],k0wXAi[1],k0wXAi[2]); 
		//
		u0X[0] = k0wXAi[0]/k0Ai; 
		u0X[1] = k0wXAi[1]/k0Ai; 
		u0X[2] = k0wXAi[2]/k0Ai; 
		// Convert components of incident beam polarization from (e1X,e2X) to (sg0X,pi0X): 
		EInSPs = PolTrn[0][0]*Eox + PolTrn[0][1]*Eoy; 
		EInSPp = PolTrn[1][0]*Eox + PolTrn[1][1]*Eoy; 
		// Calculate the wave vector for the diffracted beam. Note that here the index of 
		// refraction corrections are included. 
		kHtXAi[0] = k0wXAi[0] + HXAi[0]; 
		kHtXAi[1] = 0.; 
		kHtXAi[2] = k0wXAi[2] + HXAi[2]; 
		kHt2 = kHtXAi[0]*kHtXAi[0] + kHtXAi[2]*kHtXAi[2]; 
		kmHXAi[0] = kHtXAi[0]; 
		kmHXAi[1] = sqrt(k0Ai*k0Ai-kHt2); 
		kmHXAi[2] = kHtXAi[2]; 
		// Calculate direction cosine gamma0,reflection asymmetry parameter bee, 
		// deviation parameter Adev and normalized deviation parameter zeeC (as 
		// defined by Zachariasen gamma0, b, alpha and z.) 
		gamma0 = -u0X[1]; 
		bee = 1./(1. + HXAi[1]/k0wXAi[1]); 
		dotkH = k0wXAi[0]*HXAi[0] + k0wXAi[1]*HXAi[1] + k0wXAi[2]*HXAi[2]; 
		Adev = ( 2.*dotkH + 1./(dA*dA) )/(k0Ai*k0Ai); 
		zeeC = 0.5*( (1.-bee)*psi0c + bee*Adev ); 
		// Calculate the complex reflectivity DHsgC for sigma polarization: 
		queC = bee*psihc*psimhc; 
		sqrqzC = sqrt(queC + zeeC*zeeC); 
		del1C = 0.5*( psi0c - zeeC + sqrqzC ); 
		del2C = 0.5*( psi0c - zeeC - sqrqzC ); 
		x1C = (-zeeC+sqrqzC)/psimhc; 
		x2C = (-zeeC-sqrqzC)/psimhc; 
		ph1C = 2.*pi*iC*k0Ai*del1C*(thicum*1.E+04)/gamma0; 
		ph2C = 2.*pi*iC*k0Ai*del2C*(thicum*1.E+04)/gamma0; 
		if (real(ph1C) > logDBMax) 
			DHsgC = x2C; 
		else if (real(ph2C) > logDBMax) 
			DHsgC = x1C; 
		else { 
			Cph1C = exp(ph1C); 
			Cph2C = exp(ph2C); 
			DHsgC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C); 
		} 
		if (itrans == 1) { 
			// calculate the complex reflectivity D0trsC of the transmitted beam. 
			if (real(ph1C) > logDBMax){ 
				Cph2C = exp(ph2C); 
				D0trsC = -Cph2C*(x2C-x1C)/x1C; 
			} 
			else if (real(ph2C) > logDBMax){ 
				Cph1C = exp(ph1C); 
				D0trsC = +Cph1C*(x2C-x1C)/x2C; 
			} 
			else 
				D0trsC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C); 
		} 
		// Calculate the complex reflectivity DHpiC for pi polarization: 
		cos2t = cos(2.*thBrd); 
		queC = bee*psihc*psimhc*cos2t*cos2t; 
		sqrqzC = sqrt(queC + zeeC*zeeC); 
		del1C = 0.5*( psi0c - zeeC + sqrqzC ); 
		del2C = 0.5*( psi0c - zeeC - sqrqzC ); 
		x1C = (-zeeC+sqrqzC)/(psimhc*cos2t); 
		x2C = (-zeeC-sqrqzC)/(psimhc*cos2t); 
		ph1C = 2.*pi*iC*k0Ai*del1C*(thicum*1.E+04)/gamma0; 
		ph2C = 2.*pi*iC*k0Ai*del2C*(thicum*1.E+04)/gamma0; 
		if (real(ph1C) > logDBMax) 
			DHpiC = x2C; 
		else if (real(ph2C) > logDBMax) 
			DHpiC = x1C; 
		else { 
			Cph1C = exp(ph1C); 
			Cph2C = exp(ph2C); 
			DHpiC = x1C*x2C*(Cph2C-Cph1C)/(Cph2C*x2C-Cph1C*x1C); 
		} 
		if (itrans == 1) { 
			// calculate the complex reflectivity D0trpC of the transmitted beam. 
			if (real(ph1C) > logDBMax){ 
				Cph2C = exp(ph2C); 
				D0trpC = -Cph2C*(x2C-x1C)/x1C; 
			} 
			else if (real(ph2C) > logDBMax){ 
				Cph1C = exp(ph1C); 
				D0trpC = +Cph1C*(x2C-x1C)/x2C; 
			} 
			else 
				D0trpC = Cph1C*Cph2C*(x2C-x1C)/(Cph2C*x2C-Cph1C*x1C); 
		} 
		// Calculate the diffracted amplitudes: 
		EHSPCs = DHsgC*EInSPs; 
		EHSPCp = DHpiC*EInSPp; 
		// 
		if (itrans == 1){ 
			E0tSPs = D0trsC*EInSPs; 
			E0tSPp = D0trpC*EInSPp; 
		} 
		// Calculate diffracted amplitudes in the diffracted beam frame
		Ehx = sgH[0]*EHSPCs + piH[0]*EHSPCp; 
		Ehy = sgH[1]*EHSPCs + piH[1]*EHSPCp; 
		Ehz = sgH[2]*EHSPCs + piH[2]*EHSPCp; 
		// 
		if (itrans == 1){ 
			Etx = sgH[0]*E0tSPs + piH[0]*E0tSPp; 
			Ety = sgH[1]*E0tSPs + piH[1]*E0tSPp; 
			Etz = sgH[2]*E0tSPs + piH[2]*E0tSPp; 
		} 
		// Convert the diffracted beam's wave vector from the crystal frame to the diffracted beam frame. 
		kmHAi[0] = RXtRef[0][0]*kmHXAi[0] + RXtRef[0][1]*kmHXAi[1] + RXtRef[0][2]*kmHXAi[2]; 
		kmHAi[1] = RXtRef[1][0]*kmHXAi[0] + RXtRef[1][1]*kmHXAi[1] + RXtRef[1][2]*kmHXAi[2]; 
		kmHAi[2] = RXtRef[2][0]*kmHXAi[0] + RXtRef[2][1]*kmHXAi[1] + RXtRef[2][2]*kmHXAi[2]; 
		//
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
	{//Do nothing (?), or recalculate, because angular divergence may be changed by Crystal
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

};

//*************************************************************************

#endif
