/************************************************************************//**
 * File: sroptshp.h
 * Description: Optical element: "Shaped" (i.e. possessing some shape characteristics) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTSHP_H
#define __SROPTSHP_H

#include "sroptelm.h"

//*************************************************************************

class srTShapedOptElem : public srTGenOptElem {

protected:
	srTSRWRadStructAccessData* m_pPrevWfr;

	char m_axRot1, m_axRot2, m_axRot3;
	double m_angRot1, m_angRot2, m_angRot3;

	char m_apertShape; //1- rectangular, 2- elliptical
	double m_halfDim1, m_halfDim2; //dimensions

	srTransHndl TransHndl; //(auxiliary) space transformation of optical element in the frame of (or with respect to) the incident wavefront
	//TVector3d m_arOutFrame[3]; //center point, "horizontal" and "vertical" vectors of the output wavefront in the frame of the incident wavefront
	
	TVector3d m_eHorOut, m_eVertOut;
	double m_extAlongOptAxIn, m_extAlongOptAxOut; //positive "extents" of real input and output planes with respect to the virtual planes passing through element center

	char m_treatInOut;
	char m_wfrInterpMode; //wavefront interpolation mode at re-sampling from irregular mesh (used at propagation through "thick" mirror): 1- bilinear (based on 4 points); 2- bi-quadratic (based on 5 points).

public:
	TVector2d TransvCenPoint; // in the frame of the incident wavefront
	double Theta;  // Angle bw. Central Normal and Opt. Axis
	char RotPlane; // 'h' or 'v'

	srTShapedOptElem() 
	{
		TransvCenPoint.x = TransvCenPoint.y = 0.;
		Theta = 0.;
		RotPlane = 'h';
		m_pPrevWfr = 0;

		m_extAlongOptAxIn = m_extAlongOptAxOut = 0.;
		m_wfrInterpMode = 2; //1; //to make input parameter(?)
	}

	virtual void SetupNativeTransformation();
	virtual void SetupPreOrient(gmTrans& tr) 
	{//find original space transformation that should be applied before any parametrized rotations are applied
		tr.SetupIdent();
	}

	void FromLabToLocFrame_Point(TVector3d& P)
	{
		if(TransHndl.rep != 0) P = TransHndl.rep->TrPoint_inv(P);
	}
	void FromLabToLocFrame_Vector(TVector3d& V)
	{
		if(TransHndl.rep != 0) V = TransHndl.rep->TrBiPoint_inv(V);
	}
	void FromLocToLabFrame_Point(TVector3d& P)
	{
		if(TransHndl.rep != 0) P = TransHndl.rep->TrPoint(P);
	}
	void FromLocToLabFrame_Vector(TVector3d& V)
	{
		if(TransHndl.rep != 0) V = TransHndl.rep->TrBiPoint(V);
	}

	bool CheckIfPointIsWithinOptElem(double xLoc, double yLoc)
	{//assumes local frame of the optical element
		if((xLoc < -m_halfDim1) || (xLoc > m_halfDim1) || (yLoc < -m_halfDim2) || (yLoc > m_halfDim2)) return false; 
		
		if(m_apertShape == 2) //elliptical
		{//check if it is inside ellipse
			double xr = xLoc/m_halfDim1, yr = yLoc/m_halfDim2;
			if((xr*xr + yr*yr) > 1) return false;
		}
		return true;
	}

	virtual void CalcOutputFrame()
	{//uses srTransHndl TransHndl; and TVector3d m_arOutFrame[3]; defined in srTShapedOptElem
		//TVector3d &OutPlaneCenP = m_arOutFrame[0], &eHorOut = m_arOutFrame[1], &eVertOut = m_arOutFrame[2];
		//OutPlaneCenP.Zero();
		m_eHorOut.x = 1; m_eHorOut.y = m_eHorOut.z = 0;
		m_eVertOut.x = m_eVertOut.z = 0; m_eVertOut.y = 1;
	}

	int PropagateByRays(srTSRWRadStructAccessData*);
	int WfrInterpolOnOrigGrid(srTSRWRadStructAccessData* pRadAccessData, int ie, double* arRayOutTrCoord);

	//void FindElemExtentsAlongOptAxes(gmTrans& ElemRot, double dim1, double dim2, double& extIn, double& extOut);

	inline void FindInWaveFrontNormal(srTEFieldPtrsX&, TVector3d&);

	virtual double SurfaceFunction(double xLoc, double yLoc, char SurfNo) { return 0.;}
	virtual void SurfaceNormalAtPoint(double xLoc, double yLoc, char SurfNo, TVector3d& N) {}
	virtual void FindRayIntersectWithSurface(TVector3d* Ray, char SurfNo, TVector3d& LocPoint) {}
	virtual void SetupInAndOutPlanes(TVector3d*, TVector3d*) {}
	virtual void OneRayTrace(TVector3d* Ray, double& OptPath, double& xOut, double& zOut) {}

	virtual void TraceRayInLocFrame(TVector3d* Ray, TVector3d* arIntersectP, double& extraOptPath) {}
	virtual void ModifElecFieldAtSurf(srTEFieldPtrs& EPtrs, double PhaseShift, TVector3d* Ray, TVector3d* RayOut) {}
	virtual void GetElemDimsInLocFrame(double& horDim, double& vertDim) { horDim = vertDim = 0;}
};

//*************************************************************************

inline void srTShapedOptElem::FindInWaveFrontNormal(srTEFieldPtrsX& EFieldPtrsX, TVector3d& WaveFrN)
{//to improve !!??
	if(EFieldPtrsX.WaveFrontTermCanBeTreated)
	{
		double InvRx = 1./EFieldPtrsX.RobsX;
		double InvRy = 1./EFieldPtrsX.RobsY;

		WaveFrN.x = EFieldPtrsX.x*InvRx;
		//WaveFrN.z = EFieldPtrsX.z*InvRz;
		WaveFrN.y = EFieldPtrsX.y*InvRy;

		double nxe2_p_nye2 = WaveFrN.x*WaveFrN.x + WaveFrN.y*WaveFrN.y;
		if(nxe2_p_nye2 < 0.01)
		{
			WaveFrN.z = 1. - 0.5*nxe2_p_nye2*(1. + 0.25*nxe2_p_nye2*(1. + 0.5*nxe2_p_nye2*(1. + 0.625*nxe2_p_nye2)));
		}
		else
		{
			WaveFrN.z = sqrt(1. - nxe2_p_nye2);
		}
	}
	else
	{
		WaveFrN.x = WaveFrN.y = 0.;
		WaveFrN.z = 1.;
	}
}

//*************************************************************************

#endif
