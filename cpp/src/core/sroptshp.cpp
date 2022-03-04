/************************************************************************//**
 * File: sroptshp.cpp
 * Description: Optical element: "Shaped" (i.e. possessing some shape characteristics)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptshp.h"

//*************************************************************************
/** Obsolette
void srTShapedOptElem::SetupNativeTransformation()
{// Attention: this was written assuming Mirror !!!
	const double HalfPi = 1.5707963267949;
	TVector3d Ex(1.,0.,0.), Ez(0.,0.,1.), Zero(0.,0.,0.);
	gmTrans Rotation;

	if(RotPlane == 'h')
	{
		gmTrans AuxRotation1;
		AuxRotation1.SetupRotation(Zero, Ex, HalfPi);
		gmTrans AuxRotation2;
		AuxRotation2.SetupRotation(Zero, Ez, Theta);
		TrProduct(&AuxRotation2, &AuxRotation1, Rotation);
	}
	else if(RotPlane == 'v')
	{
		Rotation.SetupRotation(Zero, Ex, HalfPi - Theta);
	}

	TVector3d PositionVect(TransvCenPoint.x, 0., TransvCenPoint.y);
	gmTrans Translation;
	Translation.SetupTranslation(PositionVect);

	gmTrans* pTotalTrans = new gmTrans();
	TrProduct(&Translation, &Rotation, *pTotalTrans);
	TransHndl = srTransHndl(pTotalTrans);
}
**/
//*************************************************************************

void srTShapedOptElem::SetupNativeTransformation()
{//ATTENTION: from fere on, "optical" frame is assumed, i.e:
 //x-hor., y-vert., z-longitudinal
 //the transformation found here orients the optical element with respect to incident beam
	//const double Pi = 3.141592653589793;
	//double RelTol = 1.E-10;
	//TVector3d N0Loc(0.,0.,1.), Zero(0.,0.,0.);
	//TVector3d eVert(0.,1.,0.), eHor(1.,0.,0.), eLong(0.,1.,0.);
	TVector3d Zero(0.,0.,0.);

	gmTrans *pTotalTrans = new gmTrans();
	SetupPreOrient(*pTotalTrans); //virtual
	//pTotalTrans->SetupRotation(Zero, eVert, Pi);

	int arAxRot[] = {m_axRot1, m_axRot2, m_axRot3};
	double arAngRot[] = {m_angRot1, m_angRot2, m_angRot3};
	for(int i=0; i<3; i++)
	{
		int axRot = arAxRot[i];
		double angRot = arAngRot[i];
		if((axRot <= 0) || (angRot == 0)) continue;

		TVector3d vAxRot(0.,0.,0.);
		bool axIsCorrect = false;
		if(axRot == 1) 
		{
			vAxRot.x = 1.; axIsCorrect = true;
		}
		else if(axRot == 2)
		{
			vAxRot.y = 1.; axIsCorrect = true;
		}
		else if(axRot == 3)
		{
			vAxRot.z = 1.; axIsCorrect = true;
		}

		if(!axIsCorrect) continue;

		gmTrans prevTrf(*pTotalTrans), nextRot;
		nextRot.SetupRotation(Zero, vAxRot, angRot);
		TrProduct(&nextRot, &prevTrf, *pTotalTrans);
	}

	double horSizeLoc, vertSizeLoc;
	GetElemDimsInLocFrame(horSizeLoc, vertSizeLoc);

	//FindElemExtentsAlongOptAxes(*pTotalTrans, horSizeLoc, vertSizeLoc, m_extAlongOptAxIn, m_extAlongOptAxOut); //virtual
	//if(ofstLong != 0.)
	//{
	//	TVector3d vOfstTrf(0., 0., ofstLong); //Check this: can be frame-dependent!
	//	gmTrans prevTrf(*pTotalTrans), ofstTrf;
	//	ofstTrf.SetupTranslation(vOfstTrf);
	//	TrProduct(&ofstTrf, &prevTrf, *pTotalTrans);
	//}

	if((TransvCenPoint.x != 0.) || (TransvCenPoint.y != 0.))
	{
		TVector3d vOfstTrf(TransvCenPoint.x, TransvCenPoint.y, 0.); //Check this: can be frame-dependent!
		gmTrans prevTrf(*pTotalTrans), ofstTrf;
		ofstTrf.SetupTranslation(vOfstTrf);
		TrProduct(&ofstTrf, &prevTrf, *pTotalTrans);
	}

	TransHndl = srTransHndl(pTotalTrans);
}

//*************************************************************************

int srTShapedOptElem::PropagateByRays(srTSRWRadStructAccessData* pRadAccessData)
{//ATTENTION: from fere on, "optical" frame is assumed, i.e:
 //x-hor., y-vert., z-longitudinal (along optical path)
 //This function would be eventually called by implementations of "virtual int PropagateRadiationSimple"
 //no wavefront resizing here.

	//TVector3d BufInPlaneCenP, BufInPlaneNorm(0.,0.,1.); //, BufInPlaneNorm(0.,1.,0.);
	//TVector3d InPlane[] = { BufInPlaneCenP, BufInPlaneNorm };
	//TVector3d &InPlaneCenP = *InPlane, &InPlaneNorm = InPlane[1];

	//TVector3d BufOutPlaneCenP, BufOutPlaneNorm;
	//TVector3d OutPlane[] = { BufOutPlaneCenP, BufOutPlaneNorm };
	//TVector3d &OutPlaneCenP = *OutPlane, &OutPlaneNorm = OutPlane[1];
	//SetupInAndOutPlanes(InPlane, OutPlane);

	//TVector3d arOutFrame[3];
	//TVector3d &OutPlaneCenP = m_arOutFrame[0], &eHorOut = m_arOutFrame[1], &eVertOut = m_arOutFrame[2];
	//CalcOutputFrame(arOutFrame); //moved to constructor //virtual: calculates parameters of the output frame in coordinates of the input frame 

	TVector3d WaveFrN;
	TVector3d Ray[2], RayLocFr[2], RayOut[2], arIntersectP[3];
	TVector3d &RayP = Ray[0], &RayV = Ray[1];
	TVector3d &RayOutP = RayOut[0], &RayOutV = RayOut[1];
	TVector3d &RayLocFrP = RayLocFr[0], &RayLocFrV = RayLocFr[1];
	TVector3d &InSurfIntersectP = arIntersectP[0], &OutSurfIntersectP = arIntersectP[1], &OutP = arIntersectP[2];
	TVector3d vPathBef, vPathAft, vDif;
	TVector3d PlaneBefore[2], PlaneCenOut[2], vAuxIntersectP, vAuxDif;
	TVector3d &PlaneBeforeP = PlaneBefore[0], &PlaneBeforeV = PlaneBefore[1];
	TVector3d &PlaneCenOutP = PlaneCenOut[0], &PlaneCenOutV = PlaneCenOut[1];
	
	PlaneCenOutP.Zero();
	PlaneCenOutV = m_eHorOut^m_eVertOut;

	RayP.z = 0; //InPlaneCenP.z;
	double extraOptPath;

	char LocWaveFrontTermCanBeTreated = WaveFrontTermCanBeTreated(*pRadAccessData);
	gmTrans *pTrans = TransHndl.rep;

	double *arAuxRayTrCoord = new double[((pRadAccessData->nx)*(pRadAccessData->nz)) << 1];
	if(arAuxRayTrCoord == 0) return NOT_ENOUGH_MEMORY_FOR_SR_COMP;

	float *pEX0 = pRadAccessData->pBaseRadX;
	float *pEZ0 = pRadAccessData->pBaseRadZ;
	//long PerX = pRadAccessData->ne << 1;
	//long PerY = PerX*pRadAccessData->nx;
	long long PerX = pRadAccessData->ne << 1;
	long long PerY = PerX*pRadAccessData->nx;
	double ePh = pRadAccessData->eStart, x, y;

	//long nx_mi_1 = pRadAccessData->nx - 1;
	//long ny_mi_1 = pRadAccessData->nz - 1;
	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		double TwoPi_d_LambdaM = ePh*5.067681604E+06;
		//long Two_ie = ie << 1;
		long long Two_ie = ie << 1;
		double *t_arAuxRayTrCoord = arAuxRayTrCoord;

		y = pRadAccessData->zStart;

		for(long iy=0; iy<pRadAccessData->nz; iy++)
		{
			//long iyPerY = iy*PerY;
			long long iyPerY = iy*PerY;
			float *pEX_StartForX = pEX0 + iyPerY;
			float *pEZ_StartForX = pEZ0 + iyPerY;

			x = pRadAccessData->xStart;
			for(long ix=0; ix<pRadAccessData->nx; ix++)
			{
				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				float *pExRe = pEX_StartForX + ixPerX_p_Two_ie;
				float *pExIm = pExRe + 1;
				float *pEzRe = pEZ_StartForX + ixPerX_p_Two_ie;
				float *pEzIm = pEzRe + 1;

				srTEFieldPtrsX EFieldPtrsX;

				if((*pExRe != 0) || (*pExIm != 0) || (*pEzRe != 0) || (*pEzIm != 0))
				{
					EFieldPtrsX.WaveFrontTermCanBeTreated = LocWaveFrontTermCanBeTreated;
					EFieldPtrsX.x = x; EFieldPtrsX.y = y; EFieldPtrsX.ePh = ePh;
					EFieldPtrsX.xStep = pRadAccessData->xStep; EFieldPtrsX.yStep = pRadAccessData->zStep;
					EFieldPtrsX.RobsX = pRadAccessData->RobsX; EFieldPtrsX.RobsY = pRadAccessData->RobsZ;

					FindInWaveFrontNormal(EFieldPtrsX, WaveFrN);

					RayP.x = x; RayP.y = y;
					RayV = WaveFrN;

					double OptPathCorThick = 0.;
					if(m_extAlongOptAxIn != 0.) //check sign?
					{//propagate back to a plane before optical element
						PlaneBeforeP.x = PlaneBeforeP.y = 0.; PlaneBeforeP.z = -m_extAlongOptAxIn; //check sign?
						PlaneBeforeV.x = PlaneBeforeV.y = 0.; PlaneBeforeV.z = 1.;
						FindLineIntersectWithPlane(PlaneBefore, Ray, vAuxIntersectP);
						vAuxDif = vAuxIntersectP - RayP;
						OptPathCorThick = -vAuxDif.Abs();
						RayP = vAuxIntersectP;
					}

					if(pTrans != 0) //from input beam frame to local frame
					{
						RayLocFrP = pTrans->TrPoint_inv(RayP);
						RayLocFrV = pTrans->TrBiPoint_inv(RayV);
					}
					else
					{
						RayLocFrP = RayP;
						RayLocFrV = RayV;
					}

					TraceRayInLocFrame(RayLocFr, arIntersectP, extraOptPath); //virtual;
					//at return, arIntersectP should contain: 
					//[0]- intersection point with optical element surface on the Input side
					//[1]- intersection point with optical element surface on the Output side
					//[2]- intersection point with outer plane (perpendicular to optical axis after the element)
					//extraOptPath - opt. path between arIntersectP[0] and arIntersectP[1] (i.e. inside optical element)
					//input Ray should not be modified!
					//propagation is done to the plane just after the optical element!

					vPathBef = arIntersectP[0] - RayLocFrP;
					vPathAft = arIntersectP[2] - arIntersectP[1];
					double OptPath = sqrt(vPathBef.x*vPathBef.x + vPathBef.y*vPathBef.y + vPathBef.z*vPathBef.z);
					OptPath += extraOptPath + sqrt(vPathAft.x*vPathAft.x + vPathAft.y*vPathAft.y + vPathAft.z*vPathAft.z);

					if(pTrans != 0) //from local frame back to the input beam frame
					{
						OutP = pTrans->TrPoint(OutP);
						OutSurfIntersectP = pTrans->TrPoint(OutSurfIntersectP);
						InSurfIntersectP = pTrans->TrPoint(InSurfIntersectP);
					}

					RayV = InSurfIntersectP - RayP;
					RayOutP = OutSurfIntersectP;
					//RayOut[1] = OutP - OutSurfIntersectP;
					vAuxDif = OutP - OutSurfIntersectP;

					if(m_extAlongOptAxOut != 0.) //check sign?
					{//propagate back form OutP to a virtual plane in the center of optical element
						RayOutV = vAuxDif;
						RayOutV.Normalize();
						FindLineIntersectWithPlane(PlaneCenOut, RayOut, vAuxIntersectP);
						vAuxDif = vAuxIntersectP - OutP;
						OptPathCorThick += -vAuxDif.Abs();
						OutP = vAuxIntersectP;
						RayOutV = vAuxDif;
					}
					OptPath += OptPathCorThick;
					double PhaseShift = TwoPi_d_LambdaM*OptPath;

					//Transverse coordinates in the frame of the output wavefront (in "virtual" plane passing through element center):
					vDif = OutP - PlaneCenOutP; // "- PlaneCenOutP" might be unnecessary
					*(t_arAuxRayTrCoord++) = vDif*m_eHorOut; //xOut
					*(t_arAuxRayTrCoord++) = vDif*m_eVertOut; //yOut
					
					srTEFieldPtrs EPtrs(pExRe, pExIm, pEzRe, pEzIm);
					ModifElecFieldAtSurf(EPtrs, PhaseShift, Ray, RayOut); //Update field, taking into account reflection/transmission coeficients !!!
				}
				x += pRadAccessData->xStep;
			}
			y += pRadAccessData->zStep;
		}
		//Re-interpolate the output wavefront (at fixed photon energy) on the initial grid:
		WfrInterpolOnOrigGrid(pRadAccessData, ie, arAuxRayTrCoord);

		ePh += pRadAccessData->eStep;
	}

	if(arAuxRayTrCoord != 0) delete[] arAuxRayTrCoord;
	return 0;
}

//*************************************************************************

int srTShapedOptElem::WfrInterpolOnOrigGrid(srTSRWRadStructAccessData* pRadAccessData, int ie, double* arRayOutTrCoord)
{//To delete ?
	return 0;
}

//*************************************************************************

//void srTShapedOptElem::FindElemExtentsAlongOptAxes(gmTrans& ElemRot, double dim1, double dim2, double& extIn, double& extOut)
//{
	//double horSizeLoc=0, vertSizeLoc=0;
	//GetElemDimsInLocFrame(horSizeLoc, vertSizeLoc);
	//Estimate extents taking extremity points with some security margins
//}

//*************************************************************************
