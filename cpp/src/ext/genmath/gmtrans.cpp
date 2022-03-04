/************************************************************************//**
 * File: gmtrans.cpp
 * Description: Space transformation utilities
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "gmtrans.h"
#include <ctype.h>

//-------------------------------------------------------------------------

void gmTrans::SetupRotation(const TVector3d& PoiOnAxVect, const TVector3d& InAxVect, double Angle)
{
	double NormFact = 1./sqrt(InAxVect.x*InAxVect.x+InAxVect.y*InAxVect.y+InAxVect.z*InAxVect.z);
	TVector3d AxVect = NormFact*InAxVect;
	double VxVx, VyVy, VzVz;
	VxVx=AxVect.x*AxVect.x; VyVy=AxVect.y*AxVect.y; VzVz=AxVect.z*AxVect.z;

	double cosAng, sinAng, One_m_cos;
	cosAng = cos(Angle); sinAng = sin(Angle); One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	M = TMatrix3d(St0, St1, St2);
	M_inv = Matrix3d_inv(M);

	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);
	V = M0*PoiOnAxVect;
	detM = s = 1.;
	ID_No = 2;
}

//-------------------------------------------------------------------------

void gmTrans::SetupRotation(const TVector3d& vPointOnAxis, const TVector3d& vUnit1, const TVector3d& vUnit2)
{//sets up rotation about axis defined by vector product vUnit1^vUnit2 and point vPointOnAxis;
 //vectors vUnit1 and vUnit2 are assumed to have unit lengths
 //rotation angle is always: 0<angle<Pi

	const double minSinAng = 1.E-13; //to tune

	TVector3d AxVect = vUnit1^vUnit2;
	double sinAng = AxVect.Abs(); //??
	if(sinAng < minSinAng)
	{
		SetupIdent(); return;
	}
	AxVect *= (1./sinAng);

	double cosAng = vUnit1*vUnit2;
	double VxVx=AxVect.x*AxVect.x, VyVy=AxVect.y*AxVect.y, VzVz=AxVect.z*AxVect.z;

	double One_m_cos = 1. - cosAng;
	double One_m_cosVxVy, One_m_cosVxVz, One_m_cosVyVz, sinVx, sinVy, sinVz;
	One_m_cosVxVy = One_m_cos*AxVect.x*AxVect.y;
	One_m_cosVxVz = One_m_cos*AxVect.x*AxVect.z;
	One_m_cosVyVz = One_m_cos*AxVect.y*AxVect.z;
	sinVx = sinAng*AxVect.x; sinVy = sinAng*AxVect.y; sinVz = sinAng*AxVect.z;

	TVector3d St0(VxVx+cosAng*(VyVy+VzVz), One_m_cosVxVy-sinVz, One_m_cosVxVz+sinVy);
	TVector3d St1(One_m_cosVxVy+sinVz, VyVy+cosAng*(VxVx+VzVz), One_m_cosVyVz-sinVx);
	TVector3d St2(One_m_cosVxVz-sinVy, One_m_cosVyVz+sinVx, VzVz+cosAng*(VxVx+VyVy));
	M = TMatrix3d(St0, St1, St2);
	M_inv = Matrix3d_inv(M);

	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);
	V = M0*vPointOnAxis;
	detM = s = 1.;
	ID_No = 2;
}

//-------------------------------------------------------------------------

void gmTrans::SetupRotationToPermutAxes(const TVector3d& vCenPoint, char DefOrient, char Orient)
{//copied from int radTApplication::FindSpaceTransToOrientObjAlongMainAxis(double* CPoi, char DefOrient, char Orient)

	DefOrient = (char)toupper(DefOrient);
	Orient = (char)toupper(Orient);

	if(Orient == DefOrient) 
	{
		SetupIdent(); return;
	}

	const double Pi = 3.141592653589793238;
	double RotAng = 0.;
	TVector3d vAxis(1,1,1);

	if(DefOrient == 'X')
	{
		if(Orient == 'Y') 
		{
			RotAng = 2.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
		else if(Orient == 'Z') 
		{
			RotAng = 4.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
	}
	else if(DefOrient == 'Y')
	{
		if(Orient == 'X') 
		{
			RotAng = 4.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
		else if(Orient == 'Z') 
		{
			RotAng = 2.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
	}
	else if(DefOrient == 'Z')
	{
		if(Orient == 'X') 
		{
			RotAng = 2.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
		else if(Orient == 'Y') 
		{
			RotAng = 4.*Pi/3.;
			//TrfInd = SetRotation(CPoi, 3, AxisVect, 3, RotAng);
		}
	}
	else 
	{
		RotAng = 0.;
	}

	if(RotAng == 0.) SetupIdent();
	else SetupRotation(vCenPoint, vAxis, RotAng);
}

//-------------------------------------------------------------------------

void gmTrans::SetupPlaneSym(const TVector3d& PoiOnPlaneVect, const TVector3d& inN)
{
	double InvNormFact = sqrt(inN.x*inN.x+inN.y*inN.y+inN.z*inN.z);

	TVector3d N(0,0,1);
	if(InvNormFact != 0) N = (1./InvNormFact)*inN;

	TVector3d St0(1.-2.*N.x*N.x, -2.*N.x*N.y, -2.*N.x*N.z);
	TVector3d St1(St0.y, 1.-2.*N.y*N.y, -2.*N.y*N.z);
	TVector3d St2(St0.z, St1.z, 1.-2.*N.z*N.z);
	M.Str0 = St0; M.Str1 = St1; M.Str2 = St2;
	M_inv = Matrix3d_inv(M);

	TVector3d St00(1.-St0.x, -St0.y, -St0.z);
	TVector3d St01(-St1.x, 1.-St1.y, -St1.z);
	TVector3d St02(-St2.x, -St2.y, 1.-St2.z);
	TMatrix3d M0(St00, St01, St02);

	V = M0*PoiOnPlaneVect;
	detM = -1.;
	s = 1.;
	ID_No = 3; 
}

//-------------------------------------------------------------------------
