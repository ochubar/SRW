/************************************************************************//**
 * File: sropthck.h
 * Description: Optical element: "Thick" Mirror (header)
 * Project: Synchrotron Radiation Workshop
 * First release: October 2012
 *
 * Copyright (C) Brookhaven National Laboratory, Upton, NY, USA
 * Portions Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTHCK_H
#define __SROPTHCK_H

#include "sroptfoc.h"
#include "gmmeth.h"
#include "srwlib.h"

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include <stdio.h>

//*************************************************************************

//class srTOptThickElem : public srTFocusingElem {
//
//protected:
//
//	char m_axRot1, m_axRot2, m_axRot3;
//	double m_angRot1, m_angRot2, m_angRot3;
//	char m_apertShape; //1- rectangular, 2- elliptical
//
//	//TVector3d CenPointVect, CenNormVect;
//	//double RotAng;
//	//double D1, D2; //dimensions
//	//srTrans OptElemTrans;
//
//public:
//	//void NormalizeCenNormVect() {}
//	void SetupNativeTransformation();
//	void CalcOutputFrame(); //virtual
//
//};

//*************************************************************************

class srTMirror : public srTFocusingElem {
//Base class for all Mirrors and Gratings

	srTParPrecWfrPropag m_ParPrecWfrPropag;
	srTSRWRadStructAccessData* m_pRadAux;

protected:

	srTDataMD m_reflData;
	char m_propMeth;
	int m_numPartsProp;
	int m_npt, m_nps;
	TVector3d m_vCenNorm, m_vCenTang, m_vInLoc, m_vOutLoc, m_vHorOutIn, m_vVerOutIn; //Used in all main propagation methods for most mirror types
	TVector3d m_vPtOutLoc; //Used only in Test PropagateRadiationSimple_FourierByParts

	double m_longPosStartPropPart, m_longPosEndPropPart;
	double m_inWfrRh, m_inWfrRv, m_inWfrCh, m_inWfrCv;

	//Grating parameters:
	int m_grM; //Output (diffraction) order to be used
	double m_grDen, m_grDen1, m_grDen2, m_grDen3, m_grDen4; //Grove density (coefficients in the polynomial: m_grDen + m_grDen1*y + m_grDen2*y^2 + m_grDen3*y^3 + m_grDen4*y^4), in [lines/m], [lines/m^2], [lines/m^3], [lines/m^4], and [lines/m^5] respectively 
	double m_grAng; //Angle between the grove direction and the saggital direction of the substrate [rad]
	bool m_isGrating; //Switch between Mirror and Grating Reflection Laws

	double m_grAuxCosAng, m_grAuxSinAng, m_grAuxEphAvg, m_grAuxAnamorphMagnH, m_grAuxAnamorphMagnV, m_grAuxElecFldAnamorphMagnFact; //Auxiliary variables for Grating

public:

	srTMirror(const SRWLOptMir& mir);
	srTMirror(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	static srTMirror* DefineMirror(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	static srTMirror* DefineMirror(char* sType, void* pvData);
	static srTMirror* DefineGrating(char* sType, void* pvData);

	void SetupNativeTransFromLocToBeamFrame(TVector3d& vCenNorm, TVector3d& vCenTang, TVector2d& vCenP);
	int FindBasisVectorTransAndExtents();
	void FindElemExtentsAlongOptAxes(gmTrans& trfMir, TVector3d& vCenNorm, double halfDim1, double halfDim2, double& extIn, double& extOut);
	
	int PropagateRadiationSimple_LocRayTracing(srTSRWRadStructAccessData* pRadAccessData);
	int PropagateRadiationSimple_FourierByParts(srTSRWRadStructAccessData* pRadAccessData); //Test of propagation by Fourier method in steps (failed?)

	void RadPointModifier_ThinElem(srTEXZ& EXZ, srTEFieldPtrs& EPtrs);
	void RadPointModifier_FourierByParts(srTEXZ& EXZ, srTEFieldPtrs& EPtrs); //Test of propagation by Fourier method in steps (failed?)
	void EstimateFocalLengths(double radTan, double radSag); //to make it virtual in srTFocusingElem?

	//int WfrInterpolOnOrigGrid(srTSRWRadStructAccessData* pWfr, float* arRayTrCoord, float* arEX, float* arEY, float xRelOutMin, float xRelOutMax, float yRelOutMin, float yRelOutMax);
	int WfrInterpolOnOrigGrid(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, float* arEX, float* arEY, double xRelOutMin, double xRelOutMax, double yRelOutMin, double yRelOutMax);
	//int WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax);
	//int WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax);
	int WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax, double dxMax, double dzMax); //OC20082018

	virtual void FindSurfNormalInLocFrame(double x, double y, TVector3d& vN) {}
	virtual double SurfHeightInLocFrame(double x, double y) { return 0;}

	virtual bool FindRayIntersectWithSurfInLocFrame(TVector3d& inP, TVector3d& inV, TVector3d& resP, TVector3d* pResN=0) 
	{//find the intersection numerically; to imrove!!
		const int maxIt = 15;
		const double relSurfHeightTol = 1.E-15;
		const double minAbsSurfHeightTol = 1.E-18; //[m]
		
		TVector3d vPpl(0.,0.,0.), vNpl(0.,0.,1.); //Initial Point on Plane and Normal
		double t, x, y, z, zPl, absSurfHeightTol;
		int i;
		for(i=0; i<maxIt; i++)
		{
			//Finding Intersection Point with Tangential Plane
			t = ((vPpl - inP)*vNpl)/(inV*vNpl);
			x = inP.x + inV.x*t;
			y = inP.y + inV.y*t;
			zPl = inP.z + inV.z*t;

			if(i == 0)
			{
				double dx = x - inP.x, dy = y - inP.y, dz = zPl - inP.z;
				double maxDistEstim = sqrt(dx*dx + dy*dy + dz*dz);
				absSurfHeightTol = maxDistEstim*relSurfHeightTol;
				if(absSurfHeightTol < minAbsSurfHeightTol) absSurfHeightTol = minAbsSurfHeightTol;
			}

			//Finding Surface Height for (x, y) coordinates of Intersection Point with Tangential Plane
			z = SurfHeightInLocFrame(x, y);

			if(::fabs(z - zPl) < absSurfHeightTol) break;

			//Finding Normal to the new Tangential Plane
			FindSurfNormalInLocFrame(x, y, vNpl);
			//Setting the Point on the new Tangential Plane
			vPpl.x = x; vPpl.y = y; vPpl.z = z;
		}
		resP.x = x; resP.y = y; resP.z = z;
		if(pResN != 0)
		{
			FindSurfNormalInLocFrame(x, y, *pResN);
		}
		return true;

		//const double badHeightThresh = -1.E+20;
		//double ax = inV.x/inV.z, ay = inV.y/inV.z;
		//double x0 = inP.x, y0 = inP.y, z0 = inP.z;
		//double z = 0., y, x;
		//for(int i=0; i<maxIt; i++)
		//{
		//	x = ax*(z - z0) + x0;
		//	y = ay*(z - z0) + y0;
		//	z = SurfHeightInLocFrame(x, y);
		//	if(z < badHeightThresh) return false;
		//}
		//resP.x = ax*(z - z0) + x0;
		//resP.y = ay*(z - z0) + y0;
		//resP.z = z;

		//if(pResN != 0)
		//{
		//	FindSurfNormalInLocFrame(resP.x, resP.y, *pResN);
		//}
		//return true;
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect) //virtual in srTGenOptElem
	{
		m_ParPrecWfrPropag = ParPrecWfrPropag; //store for use in a composite prapagator (through drif space, etc.)
		
		if(m_isGrating)
		{
			double eFin = pRadAccessData->eStart + (pRadAccessData->ne - 1)*(pRadAccessData->eStep);
			m_grAuxEphAvg = 0.5*(pRadAccessData->eStart + eFin); //required for Grating basis vecotrs setup (in FindBasisVectorTransAndExtents)
		}

		int res = FindBasisVectorTransAndExtents(); //main reason for putting this here (and not in ctor) is that some (intersection with surfaces) functions from derived classes should be called in it
		if(res) return res;

		char &MethNo = ParPrecWfrPropag.MethNo;
		int result = 0;

		if(MethNo == 0) result = PropagateRadiationMeth_0(pRadAccessData); //int srTGenOptElem::PropagateRadiationMeth_0
		//else if(MethNo == 1) result = PropagateRadiationMeth_1(pRadAccessData);
		//else if(MethNo == 2) result = PropagateRadiationMeth_2(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
		return result;
	}

	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadDataSingleE) //virtual in srTGenOptElem
	{
		int result = 0;
		m_wfrRadWasProp = false;
		if(result = PropagateRadiationSimple(pRadAccessData)) return result; //in first place because previous wavefront radius may be required for some derived classes
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		if(!m_wfrRadWasProp) { if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;} //already propagated
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(m_propMeth == 1) return PropagateRadiationSimple_ThinElem(pRadAccessData);
		else if(m_propMeth == 2) return PropagateRadiationSimple_LocRayTracing(pRadAccessData);
		//else if(m_propMeth == 3) return PropagateRadiationSimple_LocRayTracingWithDiffr(pRadAccessData);
		else return 0;
	}

	int PropagateRadiationSimple_ThinElem(srTSRWRadStructAccessData* pRadAccessData)
	{//Thin Optical Element, however with varying reflectivity
		int res = 0;
		if(pRadAccessData->Pres != 0) if(res = SetRadRepres(pRadAccessData, 0)) return res; //to coordinate repres
		if(res = TraverseRadZXE(pRadAccessData)) return res;
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{
		if(m_propMeth == 1) { RadPointModifier_ThinElem(EXZ, EPtrs); return;}
		//else if(m_propMeth == 2) { RadPointModifier_LocRayTracing(EXZ, EPtrs); return;}
		//else if(m_propMeth == 3) { RadPointModifier_LocRayTracingWithDiffr(EXZ, EPtrs); return;}
	}

	void GetComplexReflectCoefFromTable(double phEn, double angInc, double& RsigRe, double& RsigIm, double& RpiRe, double& RpiIm)
	{//to be used only if m_reflData.pData != 0
	 //Getting complex reflecivity coefficients for Sigma and Pi components of the electric field
	 //Log scale case yet to implement

		int ne = m_reflData.DimSizes[0];
		double eStart = m_reflData.DimStartValues[0];
		double eStep = m_reflData.DimSteps[0];
		int nAng = m_reflData.DimSizes[1];
		double angStart = m_reflData.DimStartValues[1];
		double angStep = m_reflData.DimSteps[1];
		int nComp = m_reflData.DimSizes[2];

		const long perPhotEn = 2;
		//long perAng = perPhotEn*ne;
		//const long perSigPi = perAng*nAng;
		long long perAng = perPhotEn*ne;
		long long perSigPi = perAng*nAng;

		int ie = (int)((phEn - eStart)/eStep + 0.00001);
		if((phEn - (eStart + ie*eStep)) > 0.5*eStep) ie++;
		if(ie < 0) ie = 0;
		if(ie >= ne) ie = ne - 1;

		int iAng = (int)((angInc - angStart)/angStep + 0.00001);
		if((angInc - (angStart + iAng*angStep)) > 0.5*angStep) iAng++;
		if(iAng < 0) iAng = 0;
		if(iAng >= nAng) iAng = nAng - 1;

		//long ofstSig = perPhotEn*ie + perAng*iAng;
		long long ofstSig = perPhotEn*ie + perAng*iAng;
		//setting appropriate pointer type 
		if(m_reflData.DataType[1] == 'f')
		{
			float *pRsig = ((float*)(m_reflData.pData)) + ofstSig;
			if(nComp > 1)
			{
				float *pRpi = pRsig + perSigPi;
				RsigRe = *(pRsig++); RsigIm = *pRsig;
				RpiRe = *(pRpi++); RpiIm = *pRpi;
			}
			else
			{
				RsigRe = *(pRsig++); RsigIm = *pRsig;
				RpiRe = RsigRe; RpiIm = RsigIm;
			}
		}
		else
		{
			double *pRsig = ((double*)(m_reflData.pData)) + ofstSig;
			if(nComp > 1)
			{
				double *pRpi = pRsig + perSigPi;
				RsigRe = *(pRsig++); RsigIm = *pRsig;
				RpiRe = *(pRpi++); RpiIm = *pRpi;
			}
			else
			{
				RsigRe = *(pRsig++); RsigIm = *pRsig;
				RpiRe = RsigRe; RpiIm = RsigIm;
			}
		}
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{
		m_wfrRadWasProp = true; //OC190414

		if(m_isGrating)
		{
			double anaMagnHe2 = m_grAuxAnamorphMagnH*m_grAuxAnamorphMagnH;
			double RhOld = pRadAccessData->RobsX;
			double bufErrH = m_grAuxAnamorphMagnH*FocDistX/(FocDistX - anaMagnHe2*RhOld);
			double MagnH = bufErrH*m_grAuxAnamorphMagnH;
			pRadAccessData->RobsX = RhOld*MagnH;
			pRadAccessData->RobsXAbsErr *= (bufErrH*bufErrH);
			pRadAccessData->xc = (pRadAccessData->xc - TransvCenPoint.x)*MagnH;

			double anaMagnVe2 = m_grAuxAnamorphMagnV*m_grAuxAnamorphMagnV;
			double RvOld = pRadAccessData->RobsZ;
			double bufErrV = m_grAuxAnamorphMagnV*FocDistZ/(FocDistZ - anaMagnVe2*RvOld);
			double MagnV = bufErrV*m_grAuxAnamorphMagnV;
			pRadAccessData->RobsZ = RvOld*MagnV;
			pRadAccessData->RobsZAbsErr *= (bufErrV*bufErrV);
			pRadAccessData->zc = (pRadAccessData->zc - TransvCenPoint.y)*MagnV;
			return 0;
		}
		else return srTFocusingElem::PropagateWaveFrontRadius(pRadAccessData);
	}
};

//*************************************************************************

class srTMirrorPlane : public srTMirror {
	
public:

	//srTMirrorPlane(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	srTMirrorPlane(const SRWLOptMirPl& srwlMirPl) : srTMirror(srwlMirPl.baseMir) 
	{
		FocDistX = FocDistZ = 1.e+30;
	}

	bool FindRayIntersectWithSurfInLocFrame(TVector3d& inP, TVector3d& inV, TVector3d& resP, TVector3d* pResN=0) //virtual in srTMirror
	{//Returns coordinates of the intersection point in the Local frame (resP), and, optionally, components of the surface normal at the intersection point (always in the Local frame)
	 //inP and inV are respectively point and vector identifying the input ray
	 //In the Local frame: tangential direction is X, saggital Y, mirror normal is along Z, and plane equation is: z = 0
		
		//test
		//if(inV.z == 0.) inV.z = 1e-13;

		double t = -inP.z/inV.z; //inV.z can't be 0 (!?)
		resP.x = inP.x + t*inV.x;
		resP.y = inP.y + t*inV.y;
		resP.z = 0;

		//resP.x = 0;
		//resP.y = 0;
		//resP.z = 0;

		if(pResN != 0)
		{
			pResN->x = 0; pResN->y = 0; pResN->z = 1;
		}
		return true;
	}
};

//*************************************************************************

class srTMirrorEllipsoid : public srTMirror {
	
	double m_p, m_q, m_angGraz, m_radSag; //input

	double m_ax, m_ay, m_az; //ellipsoid parameters in Local frame (derived)
	double m_axE2, m_ayE2, m_azE2; //ellipsoid parameters in Local frame (derived)
	double m_xcLocNorm, m_zcLocNorm; //coordinates of mirror center in the "Local Normal" frame, where the elipse is described by x^2/m_ax^2 + y^2/m_ay^2 + z^2/m_az^2 = 1
	double m_ellPhiMin, m_ellPhiMax; //angle coordinate of mirror edges in the "Local Normal" frame
	//double m_cosAngGraz, m_sinAngGraz;
	double m_cosAngRotNormLoc, m_sinAngRotNormLoc;

public:

	//srTMirrorEllipsoid(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	srTMirrorEllipsoid(const SRWLOptMirEl& mirEl);

	void DetermineEllipsoidParamsInLocFrame()
	{//In the Local frame: tangential direction is X, saggital Y
	 //Assumes that m_p, m_q, m_angGraz, m_radSag are set !
	 //Determines m_ax, m_ay, m_az
		m_ax = 0.5*(m_p + m_q);
		m_axE2 = m_ax*m_ax;
		double twoTheta = m_angGraz*2.;
		double alp = atan(sin(twoTheta)/(m_p/m_q + cos(twoTheta)));

		if(m_vCenTang.z >= 0.) //OC170116
		{
			if(alp < 0.) alp = -alp;
		}
		else
		{
			if(alp >= 0.) alp = -alp;
		}

		double sinAlp = sin(alp);
		double sinAlpE2 = sinAlp*sinAlp;
		double q_plus_p_sinAlpE2 = m_q + m_p*sinAlpE2;
		m_azE2 = 0.5*m_p*(q_plus_p_sinAlpE2 - sqrt(q_plus_p_sinAlpE2*q_plus_p_sinAlpE2 - 4.*m_axE2*sinAlpE2));
		m_az = sqrt(m_azE2);
		double e2 = (m_axE2 - m_azE2)/m_axE2; //make general
		double x0E2 = (m_axE2 - m_p*m_q)/e2;
		double x0 = sqrt(x0E2);
		if(m_p > m_q) x0 = -x0;
		
		//double z0 = -m_p*sinAlp;
		double z0 = m_p*sinAlp; //OC170116

		double tgBet = -m_az*x0/sqrt(1. - x0E2/m_axE2);
		double tgBetE2 = tgBet*tgBet;
		double aux1 = m_axE2 + m_azE2*tgBetE2;
		double aux2 = x0 + z0*tgBet;
		double dd = sqrt((aux1 - aux2*aux2)/aux1);
		double azt = m_ax*m_az*dd/sqrt(aux1);
		m_ay = sqrt(m_radSag*azt)/dd;
		m_ayE2 = m_ay*m_ay;

		m_xcLocNorm = x0; //coordinates of mirror center in the "Local Normal" frame, where the elipse is described by x^2/m_ax^2 + y^2/m_ay^2 + z^2/m_az^2 = 1
		m_zcLocNorm = z0;
		//m_cosAngGraz = cos(m_angGraz);
		//m_sinAngGraz = sin(m_angGraz);

		double xnLocNorm = -x0/m_axE2, znLocNorm = -z0/m_azE2;
		double invNorm = 1./sqrt(xnLocNorm*xnLocNorm + znLocNorm*znLocNorm);
		m_cosAngRotNormLoc = znLocNorm*invNorm; //cos and sin of rotation angle between Local and "Local Normal" frames
		m_sinAngRotNormLoc = xnLocNorm*invNorm;

		//Coordinates of mirror edges in the Local frame:
		double x1Loc = m_halfDim1, x2Loc = -m_halfDim1;
		double z1Loc = 0., z2Loc = 0.; //to correct?

		const double Pi = 3.141592653589793;
		const double twoPi = 2.*Pi;
		const double tol = 1.E-12;
		//Coordinates of mirror edges in the "Local Normal" frame:
		double x1LocNorm = m_xcLocNorm + x1Loc*m_cosAngRotNormLoc + z1Loc*m_sinAngRotNormLoc;
		double z1LocNorm = m_zcLocNorm - x1Loc*m_sinAngRotNormLoc + z1Loc*m_cosAngRotNormLoc;
		double auxAsin = asin(x1LocNorm/m_ax);
		if(z1LocNorm >= 0.)
		{
			if(x1LocNorm >= 0.) m_ellPhiMin = auxAsin;
			else m_ellPhiMin = twoPi + auxAsin;
		}
		else m_ellPhiMin = Pi - auxAsin; //z1LocNorm < 0.

		double x2LocNorm = m_xcLocNorm + x2Loc*m_cosAngRotNormLoc + z2Loc*m_sinAngRotNormLoc;
		double z2LocNorm = m_zcLocNorm - x2Loc*m_sinAngRotNormLoc + z2Loc*m_cosAngRotNormLoc;
		auxAsin = asin(x2LocNorm/m_ax);
		if(z2LocNorm >= 0.)
		{
			if(x2LocNorm >= 0.) m_ellPhiMax = auxAsin;
			else m_ellPhiMax = twoPi + auxAsin;
		}
		else m_ellPhiMax = Pi - auxAsin; //z2LocNorm < 0.

		double dPhi = fabs(m_ellPhiMax - m_ellPhiMin);
		if(dPhi > Pi) dPhi = twoPi - dPhi;

		double testPhiMax = m_ellPhiMin + dPhi;
		if(::fabs(m_ellPhiMax - testPhiMax) < tol)
		{
			testPhiMax += twoPi;
			if(::fabs(m_ellPhiMax - testPhiMax) < tol) m_ellPhiMin += twoPi;
			else
			{
				testPhiMax = m_ellPhiMin + dPhi - twoPi;
				if(::fabs(m_ellPhiMax - testPhiMax) < tol) m_ellPhiMin -= twoPi;
				else
				{
					testPhiMax = m_ellPhiMax + dPhi;
					if(::fabs(m_ellPhiMin - testPhiMax) < tol)
					{
						m_ellPhiMin = m_ellPhiMax;
						m_ellPhiMax = testPhiMax;
					}
					else
					{
						testPhiMax += twoPi;
						if(::fabs(m_ellPhiMin - testPhiMax) < tol)
						{
							m_ellPhiMin = m_ellPhiMax + twoPi;
							m_ellPhiMax = testPhiMax;
						}
						else
						{
							testPhiMax = m_ellPhiMax + dPhi - twoPi;
							if(::fabs(m_ellPhiMin - testPhiMax) < tol)
							{
								m_ellPhiMin = m_ellPhiMax - twoPi;
								m_ellPhiMax = testPhiMax;
							}
						}
					}
				}
			}
		}
	}

	bool FindRayIntersectWithSurfInLocFrame(TVector3d& inP, TVector3d& inV, TVector3d& resP, TVector3d* pResN=0) //virtual in srTMirror
	{//returns coordinates of the intersection point in the Local frame (resP), and, optionally, components of the surface normal at the intersection point (always in the Local frame)
		
		//Coordinates of all points and vectors in the frame where the elipse is described by x^2/m_ax^2 + y^2/m_ay^2 + z^2/m_az^2 = 1:
		//double x0 = m_xcLocNorm + inP.x*m_cosAngGraz + inP.z*m_sinAngGraz;
		double x0 = m_xcLocNorm + inP.x*m_cosAngRotNormLoc + inP.z*m_sinAngRotNormLoc;
		double y0 = inP.y;
		//double z0 = m_zcLocNorm - inP.x*m_sinAngGraz + inP.z*m_cosAngGraz;
		double z0 = m_zcLocNorm - inP.x*m_sinAngRotNormLoc + inP.z*m_cosAngRotNormLoc;
		//double vx = inV.x*m_cosAngGraz + inV.z*m_sinAngGraz;
		double vx = inV.x*m_cosAngRotNormLoc + inV.z*m_sinAngRotNormLoc;
		double vy = inV.y;
		//double vz = -inV.x*m_sinAngGraz + inV.z*m_cosAngGraz;
		double vz = -inV.x*m_sinAngRotNormLoc + inV.z*m_cosAngRotNormLoc;
		double vxE2 = vx*vx, vyE2 = vy*vy, vzE2 = vz*vz;

		double vy_x0_mi_vx_y0 = vy*x0 - vx*y0;
		double vz_x0_mi_vx_z0 = vz*x0 - vx*z0;
		double vz_y0_mi_vy_z0 = vz*y0 - vy*z0;
		double argRoot = -m_azE2*vy_x0_mi_vx_y0*vy_x0_mi_vx_y0 + m_ayE2*(m_azE2*vxE2 + m_axE2*vzE2 - vz_x0_mi_vx_z0*vz_x0_mi_vx_z0) + m_axE2*(m_azE2*vyE2 - vz_y0_mi_vy_z0*vz_y0_mi_vy_z0);
		if(argRoot < 0) return false;
		
		double ax_ay_az = m_ax*m_ay*m_az;
		double a = m_ayE2*m_azE2*vx*x0 + m_axE2*m_azE2*vy*y0 + m_axE2*m_ayE2*vz*z0;
		double b = m_axE2*m_azE2*vyE2 + m_ayE2*(m_azE2*vxE2 + m_axE2*vzE2);

		//double t0 = (m_p >= m_q)? -(a + sqrt(argRoot))/b : -(a - sqrt(argRoot))/b; //to check
		//double t0 = (m_p < m_q)? -(a + ax_ay_az*sqrt(argRoot))/b : -(a - ax_ay_az*sqrt(argRoot))/b; //to check
		//OC06072017: it looks like the second option of t0 should be used at any m_p,  m_q

			//double testTerm = ax_ay_az*sqrt(argRoot); //OCTEST
			//if(fabs(a - testTerm) < fabs(testTerm)*1e-13)
			//{
			//	int aha = 1;
			//}

		double t0 = -(a - ax_ay_az*sqrt(argRoot))/b; //OC06072017

		//Coordinates of the Intersection Point in the "Local Normal" frame:
		double xi = vx*t0 + x0;
		double yi = vy*t0 + y0;
		double zi = vz*t0 + z0;

		//Verification if angular coordinate of the Intersection Point is withing allowable limits
		const double Pi = 3.141592653589793;
		const double twoPi = 2.*Pi;
		double phi = asin(xi/m_ax);
		if(zi >= 0.)
		{
			if(xi < 0.) phi += twoPi;
		}
		else phi = Pi - phi;

		bool phiIsInside = false;

		if(m_ellPhiMax < m_ellPhiMin) //OC09072017
		{
			const double twoPi = 2.*3.141592653589793;

			if((m_ellPhiMin - twoPi <= phi) && (phi <= m_ellPhiMax)) phiIsInside = true;
			else
			{
				double testPhi = phi + twoPi;
				if((m_ellPhiMin - twoPi <= testPhi) && (testPhi <= m_ellPhiMax)) phiIsInside = true;
				else
				{
					testPhi = phi - twoPi;
					if((m_ellPhiMin - twoPi <= testPhi) && (testPhi <= m_ellPhiMax)) phiIsInside = true;
				}
			}

			if(!phiIsInside)
			{
				if((m_ellPhiMin <= phi) && (phi <= m_ellPhiMax + twoPi)) phiIsInside = true;
				else
				{
					double testPhi = phi + twoPi;
					if((m_ellPhiMin <= testPhi) && (testPhi <= m_ellPhiMax + twoPi)) phiIsInside = true;
					else
					{
						testPhi = phi - twoPi;
						if((m_ellPhiMin - twoPi <= testPhi) && (testPhi <= m_ellPhiMax + twoPi)) phiIsInside = true;
					}
				}
			}
		}
		else
		{
			if((m_ellPhiMin <= phi) && (phi <= m_ellPhiMax)) phiIsInside = true;
			else
			{
				double testPhi = phi + twoPi;
				if((m_ellPhiMin <= testPhi) && (testPhi <= m_ellPhiMax)) phiIsInside = true;
				else
				{
					testPhi = phi - twoPi;
					if((m_ellPhiMin <= testPhi) && (testPhi <= m_ellPhiMax)) phiIsInside = true;
				}
			}
		}

		if(!phiIsInside) 
		{
			return false;
		}

			//test: instant radius of curvature
			//double auxInvSqrt01 = 1./sqrt(1. - xi*xi/m_axE2);
			//double ziP = m_az*xi*auxInvSqrt01/m_axE2;
			//double auxSqrt = sqrt(1. + ziP*ziP);
			//double ziP2 = (m_az*auxInvSqrt01/m_axE2)*(1. + xi*xi*auxInvSqrt01*auxInvSqrt01/m_axE2);
			//double radCurv = auxSqrt*auxSqrt*auxSqrt/ziP2;

		//Transforming coordinates back to the Local frame
		double xi_mi_m_xcLocNorm = xi - m_xcLocNorm;
		double zi_mi_m_zcLocNorm = zi - m_zcLocNorm;
		//resP.x = xi_mi_m_xcLocNorm*m_cosAngGraz - zi_mi_m_zcLocNorm*m_sinAngGraz;
		resP.x = xi_mi_m_xcLocNorm*m_cosAngRotNormLoc - zi_mi_m_zcLocNorm*m_sinAngRotNormLoc;
		resP.y = yi;
		//resP.z = xi_mi_m_xcLocNorm*m_sinAngGraz + zi_mi_m_zcLocNorm*m_cosAngGraz;
		resP.z = xi_mi_m_xcLocNorm*m_sinAngRotNormLoc + zi_mi_m_zcLocNorm*m_cosAngRotNormLoc;

		if(pResN != 0)
		{   //Components of the normal vector in the frame where the elipse is described by x^2/m_ax^2 + y^2/m_ay^2 + z^2/m_az^2 = 1:
			double xnLocNorm = -xi/m_axE2, ynLocNorm = -yi/m_ayE2, znLocNorm = -zi/m_azE2;
			double invNorm = 1./sqrt(xnLocNorm*xnLocNorm + ynLocNorm*ynLocNorm + znLocNorm*znLocNorm);
			xnLocNorm *= invNorm; ynLocNorm *= invNorm; znLocNorm *= invNorm;

			//Same components in the Local frame:
			//pResN->x = xnLocNorm*m_cosAngGraz - znLocNorm*m_sinAngGraz;
			pResN->x = xnLocNorm*m_cosAngRotNormLoc - znLocNorm*m_sinAngRotNormLoc;
			pResN->y = ynLocNorm;
			//pResN->z = xnLocNorm*m_sinAngGraz + znLocNorm*m_cosAngGraz;
			pResN->z = xnLocNorm*m_sinAngRotNormLoc + znLocNorm*m_cosAngRotNormLoc;
			pResN->Normalize();

			//test
			//TVector3d inVaux = inV;
			//inVaux.Normalize();
			//double auxInstGrazAng = 1.5707963267948966 - acos(-((*pResN)*inVaux));
			//double instFocLen = 0.5*radCurv*sin(auxInstGrazAng);
			//int aha = 1;
		}
		return true;
	}

/**
	void FindSurfNormalInLocFrame(double x, double y, TVector3d& vN) //virtual in srTMirror
	{//In Local frame: tangential direction is X, saggital Y; output vector is normalized to 1
		//Coordinates of all points and vectors in the frame where the elipse is described by x^2/m_ax^2 + y^2/m_ay^2 + z^2/m_az^2 = 1:
		double x0 = m_xcLocNorm + inP.x*m_cosAngGraz + inP.z*m_sinAngGraz;
		double y0 = inP.y;

		double sqrt1 = sqrt(m_Rs*m_Rs - y*y);
		double R_mi_r_p_radS = m_Rt - m_Rs + sqrt1;
		double inv_sqrt2 = 1./sqrt(R_mi_r_p_radS*R_mi_r_p_radS - x*x);
		double nx = -x*inv_sqrt2;
		double ny = -y*R_mi_r_p_radS*inv_sqrt2/sqrt1;
		double inv_norm = 1./sqrt(nx*nx + ny*ny + 1.);
		vN.x = nx*inv_norm;
		vN.y = ny*inv_norm;
		vN.z = inv_norm;
	}

	double SurfHeightInLocFrame(double x, double y) //virtual in srTMirror
	{
		const double badRes = -1.E+23; //to return in case inanything goes wrong
		double ry = y/m_Rs;
		double rye2 = ry*ry;
		if(rye2 > 1.) return badRes;
		double a1 =(CGenMathMeth::radicalOnePlusSmallMinusOne(-rye2))*m_Rs/m_Rt;
		double rx = x/m_Rt;
		double a2 = a1*(a1 + 2.) - rx*rx;
		if(a2 < -1.) return badRes;
		return -m_Rt*(CGenMathMeth::radicalOnePlusSmallMinusOne(a2));
	}

**/
};

//*************************************************************************

class srTMirrorSphere : public srTMirror {

	double m_rad; //input

public:

	srTMirrorSphere(const SRWLOptMirSph& mirSph);

	bool FindRayIntersectWithSurfInLocFrame(TVector3d& inP, TVector3d& inV, TVector3d& resP, TVector3d* pResN = 0) //virtual in srTMirror
	{//Returns coordinates of the intersection point in the Local frame (resP), and, optionally, components of the surface normal at the intersection point ((*pResN), always in the Local frame)
	 //inP and inV are respectively point and vector identifying the input ray
	 //In the Local frame: tangential direction is X, saggital Y, mirror normal is along Z, and the sphere equation is: x^2 + y^2 + (z - r)^2 = r^2

		double ax = inV.x/inV.z, ay = inV.y/inV.z;
		double axe2 = ax*ax, aye2 = ay*ay;
		double axe2_p_aye2_p_1 = axe2 + aye2 + 1.;
		double A = m_rad - ax*inP.x - ay*inP.y + (axe2 + aye2)*inP.z;
		double x0_mi_ax_z0 = inP.x - ax*inP.z;
		double y0_mi_ay_z0 = inP.y - ay*inP.z;
		double argR = A*A - axe2_p_aye2_p_1*(x0_mi_ax_z0*x0_mi_ax_z0 + y0_mi_ay_z0*y0_mi_ay_z0);
		if(argR < 0.) return false;

		double R = sqrt(argR);
		double inv_axe2_p_aye2_p_1 = 1./axe2_p_aye2_p_1;

		//To check this solution:
		resP.z = (m_rad > 0.)? (A - R)*inv_axe2_p_aye2_p_1 : (A + R)*inv_axe2_p_aye2_p_1;
		double dzs = resP.z - inP.z;
		resP.x = inP.x + ax*dzs;
		resP.y = inP.y + ay*dzs;

		if(pResN != 0)
		{
			if(m_rad > 0.)
			{
				pResN->x = -resP.x;
				pResN->y = -resP.y;
				pResN->z = m_rad - resP.z;
			}
			else
			{
				pResN->x = resP.x;
				pResN->y = resP.y;
				pResN->z = resP.z - m_rad; //?
			}
			pResN->Normalize();
		}
		return true;
	}

	double SurfHeightInLocFrame(double x, double y) //virtual in srTMirror
	{
		double aSmall = -(x*x + y*y) / (m_rad*m_rad);
		return -m_rad*CGenMathMeth::radicalOnePlusSmallMinusOne(aSmall);
	}

	void FindSurfNormalInLocFrame(double x, double y, TVector3d& vN) //virtual in srTMirror
	{//In Local frame: tangential direction is X, saggital Y; output vector is normalized to 1
	 //In the Local frame: tangential direction is X, saggital Y, mirror normal is along Z, and the sphere equation is: x^2 + y^2 + (z - r)^2 = r^2

		double aSmall = -(x*x + y*y)/(m_rad*m_rad);
		double zi = -m_rad*CGenMathMeth::radicalOnePlusSmallMinusOne(aSmall);

		if(m_rad > 0.)
		{
			vN.x = -x;
			vN.y = -y;
			vN.z = m_rad - zi;
		}
		else
		{
			vN.x = x;
			vN.y = y;
			vN.z = zi - m_rad; //?
		}
		vN.Normalize();
	}
};

//*************************************************************************

class srTMirrorToroid : public srTMirror {
	
	double m_Rt, m_Rs;

public:

	srTMirrorToroid(srTStringVect* pElemInfo, srTDataMD* pExtraData);
	srTMirrorToroid(const SRWLOptMirTor& mirTor);

	//bool FindRayIntersectWithSurfInLocFrame(TVector3d& inP, TVector3d& inV, TVector3d& resP, TVector3d* pResN = 0) 
	//use the iterative version of srTMirror

	double SurfHeightInLocFrame(double x, double y) //virtual in srTMirror
	{
		const double badRes = -1.E+23; //to return in case if anything goes wrong, i.e. if (x, y) don't belong to torus
		double ry = y/m_Rs;
		double rye2 = ry*ry;
		if(rye2 > 1.) return badRes;
		double a1 =(CGenMathMeth::radicalOnePlusSmallMinusOne(-rye2))*m_Rs/m_Rt;
		double rx = x/m_Rt;
		double a2 = a1*(a1 + 2.) - rx*rx;
		if(a2 < -1.) return badRes;
		return -m_Rt*(CGenMathMeth::radicalOnePlusSmallMinusOne(a2));
	}

	void FindSurfNormalInLocFrame(double x, double y, TVector3d& vN) //virtual in srTMirror
	{//In the Local frame: tangential direction is X, saggital Y; output vector is normalized to 1
		//const double badRes = -1.E+23; //to return in case if anything goes wrong, i.e. if (x, y) don't belong to torus
		vN.x = 0.; vN.y = 0.; vN.y = 0.;

		double ry = y/m_Rs;
		double rye2 = ry*ry;
		if(rye2 > 1.) return;

		double radSmi1 = CGenMathMeth::radicalOnePlusSmallMinusOne(-rye2);
		double a1 = radSmi1*m_Rs/m_Rt;
		double rx = x/m_Rt;
		double a2 = a1*(a1 + 2.) - rx*rx;
		if(a2 < -1.) return;
		double invRad = 1./(CGenMathMeth::radicalOnePlusSmallMinusOne(a2) + 1.);

		vN.x = -rx*invRad;
		vN.y = -ry*invRad*(a1 + 1.)/(radSmi1 + 1.);
		vN.z = 1.;
		vN.Normalize();

/**
		double sqrt1 = sqrt(m_Rs*m_Rs - y*y);
		double R_mi_r_p_radS = m_Rt - m_Rs + sqrt1;
		double inv_sqrt2 = 1./sqrt(R_mi_r_p_radS*R_mi_r_p_radS - x*x);
		double nx = -x*inv_sqrt2;
		double ny = -y*R_mi_r_p_radS*inv_sqrt2/sqrt1;
		double inv_norm = 1./sqrt(nx*nx + ny*ny + 1.);
		vN.x = nx*inv_norm;
		vN.y = ny*inv_norm;
		vN.z = inv_norm;
**/
	}
};

//*************************************************************************
//OBSOLETE?
class srTThickMirrorGen : public srTFocusingElem {
//Perhaps this should be a Base class for all mirrors(?)
	srTDataMD m_surfData;
	double m_ampReflectPerp, m_phaseShiftPerp;
	double m_ampReflectPar, m_phaseShiftPar;

public:

	srTThickMirrorGen(srTStringVect* pElemInfo, srTDataMD* pExtraData);

	void SetupPreOrient(gmTrans& tr) //virtual in srTShapedOptElem
	{//find original space transformation that should be applied before any parametrized rotations are applied
		//tr.SetupIdent();
		const double Pi = 3.141592653589793;
		TVector3d Zero(0.,0.,0.), eVert(0.,1.,0.);
		tr.SetupRotation(Zero, eVert, Pi);
	}
	void GetElemDimsInLocFrame(double& horDim, double& vertDim) //virtual in srTShapedOptElem
	{
		//double xLocStart = m_surfData.DimStartValues[0];
		double xLocStep = m_surfData.DimSteps[0];
		long xNp = m_surfData.DimSizes[0];
		horDim = xLocStep*(xNp - 1);

		double yLocStep = m_surfData.DimSteps[1];
		long yNp = m_surfData.DimSizes[1];
		vertDim = yLocStep*(yNp - 1);
	}
	//inline void CalcOutputFrame(); //virtual in srTShapedOptElem
};

//*************************************************************************
/**
inline void srTThickMirrorGen::CalcOutputFrame() //virtual in srTShapedOptElem
{//uses srTransHndl TransHndl; and TVector3d m_eHorOut, m_eVertOut; defined in srTShapedOptElem
 //To call only after TransHndl has been set up !

	TVector3d &OutPlaneCenP = m_arOutFrame[0], &eHorOut = m_arOutFrame[1], &eVertOut = m_arOutFrame[2];

	OutPlaneCenP.Zero();
	eHorOut.x = 1; eHorOut.y = eHorOut.z = 0;
	eVertOut.x = eVertOut.z = 0; eVertOut.y = 1;
	if(TransHndl.rep == 0) return;

	TVector3d vNorm0(0,0,1); //central normal to reflecting surface (in local frame first)
	vNorm0 = TransHndl.rep->TrBiPoint(vNorm0); //in the frame of incident wfr

	ReflectVect(vNorm0, eHorOut);
	ReflectVect(vNorm0, eVertOut);
	//TVector3d vNormOut = eHorOut^eVertOut; //normal to the plane of output wavefront

	TVector3d RayCenIn[2];
	TVector3d &RayCenIn_P = RayCenIn[0], &RayCenIn_V = RayCenIn[1];
	RayCenIn_P.Zero();
	RayCenIn_V.x = RayCenIn_V.y = 0.; RayCenIn_V.z = 1.;
	RayCenIn_P = TransHndl.rep->TrPoint_inv(RayCenIn_P); //to local frame
	RayCenIn_V = TransHndl.rep->TrBiPoint_inv(RayCenIn_V);

	//find central point in output plane:

	TVector3d IntersectP;
	char SurfNoDummy = 1;
	FindRayIntersectWithSurface(RayCenIn, SurfNoDummy, IntersectP);
}
**/
//*************************************************************************
/**
class srTThickMirrorToroid : public srTOptThickElem {
	double Rt, Rs;
	double ReflCoefInt;

public:

	srTThickMirrorToroid(srTStringVect* pElemInfo)
	{
        ReflCoefInt = atof((*pElemInfo)[1]); // intencity reflectivity
		Rt = atof((*pElemInfo)[2]); // tangential radius
		Rs = atof((*pElemInfo)[3]); // sagittal radius
		D1 = atof((*pElemInfo)[4]); // size in tangential plane
		D2 = atof((*pElemInfo)[5]); // size in sagittal plane
		ApertShape = atoi((*pElemInfo)[6]); // size in sagittal plane

		CenPointVect.y = atof((*pElemInfo)[7]); // longitudinal coordinate of the center point in laboratory frame
		CenPointVect.x = atof((*pElemInfo)[8]); // horizontal coordinate of the center point in laboratory frame
		CenPointVect.z = atof((*pElemInfo)[9]); // vertical coordinate of the center point in laboratory frame

		CenNormVect.x = atof((*pElemInfo)[10]); // horizontal coordinate of the central normal vector in laboratory frame
		CenNormVect.y = atof((*pElemInfo)[11]); // longitudinal coordinate of the central normal vector in laboratory frame
		CenNormVect.z = atof((*pElemInfo)[12]); // vertical coordinate of the central normal vector in laboratory frame
		RotAng = atof((*pElemInfo)[13]); // vertical coordinate of the central normal vector in laboratory frame
	}
};
**/
//*************************************************************************
/**
class srTSpherMirror : public srTFocusingElem {
	TVector3d OutPlaneInLocFrame[2];
	TVector3d ExRefInLabFrameBeforeProp, EzRefInLabFrameBeforeProp;

public:
	double CurvRad, Dx, Dy;
	char UseSpherMirrorApprox;

	srTGenOptElemHndl NativeApertureHndl;
	srTGenOptElemHndl SpherMirrorApproxHndl;

	srTSpherMirror(srTStringVect* pElemInfo) 
	{
		CurvRad = atof((*pElemInfo)[1]); // input in m
		Dx = atof((*pElemInfo)[2]); // input in m
		Dy = atof((*pElemInfo)[3]); // input in m

		char* BufString = (*pElemInfo)[4];

		if((!strcmp(BufString, "Horizontal")) || (!strcmp(BufString, "Hor")) || (!strcmp(BufString, "Hor.")) || (!strcmp(BufString, "hor")) || (!strcmp(BufString, "hor.")))
			RotPlane = 'h';
		else if((!strcmp(BufString, "Vertical")) || (!strcmp(BufString, "Ver")) || (!strcmp(BufString, "Ver.")) || (!strcmp(BufString, "ver")) || (!strcmp(BufString, "ver.")))
			RotPlane = 'v';
		else { ErrorCode = ERROR_IN_OPTICAL_ELEMENT_DEFINITION; return;}

		Theta = atof((*pElemInfo)[5]); // input in r

		TransvCenPoint.x = atof((*pElemInfo)[6]); // input in m
		TransvCenPoint.y = atof((*pElemInfo)[7]); // input in m
	
		const double Pi = 3.1415926535898;
		double HalfCurvRad = 0.5*CurvRad;
		double CosTheta = cos(Theta);
		FocDistX = (RotPlane == 'h')? HalfCurvRad*CosTheta : HalfCurvRad/CosTheta; // input in m
		FocDistZ = (RotPlane == 'v')? HalfCurvRad*CosTheta : HalfCurvRad/CosTheta; // input in m

		SetupNativeTransformation();
		SetupNativeAperture();
		SetupSpherMirrorApprox();

		UseSpherMirrorApprox = 0;
	}
	srTSpherMirror() {}

	void SetupNativeAperture();
	void SetupSpherMirrorApprox();

	//int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, int MethNo, srTRadResizeVect& ResBeforeAndAfterVect)
	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect)
	{
		//if(UseSpherMirrorApprox) return ((srTGenOptElem*)(SpherMirrorApproxHndl.rep))->PropagateRadiation(pRadAccessData, MethNo, ResBeforeAndAfterVect);
		if(UseSpherMirrorApprox) return ((srTGenOptElem*)(SpherMirrorApproxHndl.rep))->PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResBeforeAndAfterVect);
		else
		{
            char &MethNo = ParPrecWfrPropag.MethNo;

			if(MethNo == 0) return NativePropagateRadiationMeth_0(pRadAccessData);
			else if(MethNo == 1) return NativePropagateRadiationMeth_1(pRadAccessData);
			return 0;
		}
	}
	int NativePropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		//int MethNo = 0;
		srTRadResizeVect RadResizeVect;

        srTParPrecWfrPropag ManParPrecWfrPropag(0, 0, 0, 1., 0.5);

		//if(result = ((srTGenOptElem*)(NativeApertureHndl.rep))->PropagateRadiation(pRadAccessData, MethNo, RadResizeVect)) return result;
		if(result = ((srTGenOptElem*)(NativeApertureHndl.rep))->PropagateRadiation(pRadAccessData, ManParPrecWfrPropag, RadResizeVect)) return result;
		
		if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = PropagateByRays(pRadAccessData)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}
	int NativePropagateRadiationMeth_1(srTSRWRadStructAccessData* pRadAccessData)
	{
		int result;
		//int MethNo = 1;
		srTRadResizeVect RadResizeVect;

        srTParPrecWfrPropag ParPrecWfrPropag(1, 1, 1, 1., 0.5);

		//if(result = ((srTGenOptElem*)(NativeApertureHndl.rep))->PropagateRadiation(pRadAccessData, MethNo, RadResizeVect)) return result;
		if(result = ((srTGenOptElem*)(NativeApertureHndl.rep))->PropagateRadiation(pRadAccessData, ParPrecWfrPropag, RadResizeVect)) return result;

		srTRadResize PostResize;
		PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;
		if(result = TuneRadForPropMeth_1(pRadAccessData, PostResize)) return result;

		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = PropagateByRays(pRadAccessData)) return result;
		if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;

		const double ResizeTol = 0.15;
		char PostResizeNeeded = (::fabs(PostResize.pxm - 1.) || ::fabs(PostResize.pzm - 1.) || ::fabs(PostResize.pxd - 1.) || ::fabs(PostResize.pzd - 1.));
		if(PostResizeNeeded) if(result = RadResizeGen(*pRadAccessData, PostResize)) return result;

		if(result = Propagate4x4PropMatr(pRadAccessData)) return result;
		return 0;
	}

	double SurfaceFunction(double xLoc, double yLoc, char SurfNo) 
	{
		double SqRt = sqrt(CurvRad*CurvRad - xLoc*xLoc - yLoc*yLoc);
		return (CurvRad >= 0.)? (CurvRad - SqRt) : (CurvRad + SqRt);
	}
	void SurfaceNormalAtPoint(double xLoc, double yLoc, char SurfNo, TVector3d& N)
	{// xLoc, yLoc - Local coord.;
	 // N is in Local frame
		double zLoc = SurfaceFunction(xLoc, yLoc, SurfNo);
		N.x = -xLoc; N.y = -yLoc; N.z = CurvRad - zLoc;
		N = (1./CurvRad)*N;
	}
	void FindRayIntersectWithSurface(TVector3d* Ray, char SurfNo, TVector3d& LocP)
	{// Assumes Ray in Local frame
	 // Assumes V to be unit vector
		TVector3d &R0 = *Ray, &V = Ray[1];
		double VR0 = V*R0;
		double x0x0 = R0.x*R0.x, y0y0 = R0.y*R0.y;
		double Buf1 = -CurvRad*V.z + VR0;
		double SqRt = sqrt(Buf1*Buf1 - (x0x0 + y0y0 + R0.z*(R0.z - 2.*CurvRad)));
		char SignR = (CurvRad >= 0.)? 1 : -1;
		double RVz_mi_VR0_p_SignRSqRt = CurvRad*V.z - VR0 + SignR*SqRt;
		LocP.x = V.x*RVz_mi_VR0_p_SignRSqRt + R0.x;
		LocP.y = V.y*RVz_mi_VR0_p_SignRSqRt + R0.y;
		LocP.z = V.z*RVz_mi_VR0_p_SignRSqRt + R0.z;
	}
	void OneRayTrace(TVector3d* InRay, double& OptPath, double& xOut, double& zOut)
	{// Assumes Ray in Natural Lab frame
	 // Gives xOut, zOut in Natural Lab frame
		TVector3d Ray[2]; *Ray = *InRay; Ray[1] = InRay[1];
		TVector3d &R0 = *Ray, &V = Ray[1];
		FromLabToLocFrame_Point(R0);
		FromLabToLocFrame_Vector(V);
		TVector3d LocP;
		FindRayIntersectWithSurface(Ray, 0, LocP);
		TVector3d LocP_mi_R0 = LocP - R0, LocN;
		OptPath = sqrt(LocP_mi_R0.x*LocP_mi_R0.x + LocP_mi_R0.y*LocP_mi_R0.y + LocP_mi_R0.z*LocP_mi_R0.z);
		SurfaceNormalAtPoint(LocP.x, LocP.y, 0, LocN);
		R0 = LocP;
		ReflectVect(LocN, V);
		FindLineIntersectWithPlane(OutPlaneInLocFrame, Ray, LocP);
		LocP_mi_R0 = LocP - R0;
		OptPath += sqrt(LocP_mi_R0.x*LocP_mi_R0.x + LocP_mi_R0.y*LocP_mi_R0.y + LocP_mi_R0.z*LocP_mi_R0.z);
		FromLocToLabFrame_Point(LocP);
		xOut = LocP*ExRefInLabFrameBeforeProp;
		zOut = LocP*EzRefInLabFrameBeforeProp;
	}

	void SetupInAndOutPlanes(TVector3d*, TVector3d*);

	int RangeShouldBeAdjustedAtPropag() { return 0;} // Or switch it On
	int ResolutionShouldBeAdjustedAtPropag() { return 1;}
};
**/
//*************************************************************************

#endif
