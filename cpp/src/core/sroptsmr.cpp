/************************************************************************//**
 * File: sroptsmr.cpp
 * Description: Optical element: Spherical Mirror - Obsolete (to be moved to sropthck.h, sropthck.cpp)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptsmr.h"
#include "sroptapt.h"

//*************************************************************************

void srTSpherMirror::SetupNativeAperture()
{
	double CosTheta = cos(Theta);
	double ApDx = (RotPlane == 'h')? Dx*CosTheta : Dx;
	double ApDz = (RotPlane == 'v')? Dy*CosTheta : Dy;
	NativeApertureHndl = srTGenOptElemHndl(new srTRectAperture(ApDx, ApDz, TransvCenPoint.x, TransvCenPoint.y));
}

//*************************************************************************

void srTSpherMirror::SetupSpherMirrorApprox()
{
	SpherMirrorApproxHndl = srTGenOptElemHndl(new srTCompositeOptElem());
	srTGenOptElemHndl ThinLensHndl = srTGenOptElemHndl(new srTThinLens(FocDistX, FocDistZ, TransvCenPoint.x, TransvCenPoint.y));
	((srTCompositeOptElem*)(SpherMirrorApproxHndl.rep))->AddOptElemBack(NativeApertureHndl);
	((srTCompositeOptElem*)(SpherMirrorApproxHndl.rep))->AddOptElemBack(ThinLensHndl);
}

//*************************************************************************

void srTSpherMirror::SetupInAndOutPlanes(TVector3d* InPlane, TVector3d* OutPlane)
{// Assumes InPlaneNorm and transverse part InPlaneCenPo already defined !!!
	TVector3d &InPlaneCenPo = *InPlane, &InPlaneNorm = InPlane[1];
	TVector3d &OutPlaneCenPo = *OutPlane, &OutPlaneNorm = OutPlane[1];

	TVector3d LocInPlaneNorm = InPlaneNorm;
	FromLabToLocFrame_Vector(LocInPlaneNorm);

	double xP = -0.5*Dx, yP = -0.5*Dy;
	TVector3d P0(xP, yP, SurfaceFunction(xP, yP, 0));
	TVector3d P1(-xP, yP, SurfaceFunction(-xP, yP, 0));
	TVector3d P2(xP, -yP, SurfaceFunction(xP, -yP, 0));
	TVector3d P3(-xP, -yP, SurfaceFunction(-xP, -yP, 0));
	TVector3d EdgePoints[] = { P0, P1, P2, P3 };

	int LowestInd, UppestInd;
	FindLowestAndUppestPoints(LocInPlaneNorm, EdgePoints, 4, LowestInd, UppestInd);

	TVector3d PointForInPlane = EdgePoints[LowestInd];
	FromLocToLabFrame_Point(PointForInPlane);
	InPlaneCenPo.y = PointForInPlane.y;

	TVector3d LocInPlaneCenPo = InPlaneCenPo;
	FromLabToLocFrame_Point(LocInPlaneCenPo);
	TVector3d LocInPlane[] = { LocInPlaneCenPo, LocInPlaneNorm };

	TVector3d LocP, LocN;
	FindRayIntersectWithSurface(LocInPlane, 0, LocP);
	SurfaceNormalAtPoint(LocP.x, LocP.y, 0, LocN);

	TVector3d LocOutPlaneNorm = LocInPlaneNorm;
	ReflectVect(LocN, LocOutPlaneNorm);
	FindLowestAndUppestPoints(LocOutPlaneNorm, EdgePoints, 4, LowestInd, UppestInd);
	OutPlaneCenPo = (EdgePoints[UppestInd]*LocOutPlaneNorm)*LocOutPlaneNorm;
	*OutPlaneInLocFrame = OutPlaneCenPo; OutPlaneInLocFrame[1] = LocOutPlaneNorm;
	FromLocToLabFrame_Point(OutPlaneCenPo);
	OutPlaneNorm = LocOutPlaneNorm;
	FromLocToLabFrame_Vector(OutPlaneNorm);

	TVector3d LabN = LocN;
	FromLocToLabFrame_Vector(LabN);
	ExRefInLabFrameBeforeProp = TVector3d(1.,0.,0.);
	ReflectVect(LabN, ExRefInLabFrameBeforeProp);
	EzRefInLabFrameBeforeProp = TVector3d(0.,0.,1.);
	ReflectVect(LabN, EzRefInLabFrameBeforeProp);
}

//*************************************************************************
