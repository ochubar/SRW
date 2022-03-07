/************************************************************************//**
 * File: sroptsmr.h
 * Description: Optical element: Spherical Mirror header - Obsolete (to be moved to sropthck.h, sropthck.cpp)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTSMR_H
#define __SROPTSMR_H

#include "sroptcnt.h"
#include "sroptfoc.h"

//*************************************************************************

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
	
		//const double Pi = 3.1415926535898;
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

		//const double ResizeTol = 0.15;
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

//*************************************************************************

#endif
