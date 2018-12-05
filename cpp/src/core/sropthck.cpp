/************************************************************************//**
 * File: sropthck.cpp
 * Description: Optical element: "Thick" Mirror
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

#include "sropthck.h"
#include "sroptdrf.h"
#include "gminterp.h"

//*************************************************************************

srTMirror::srTMirror(srTStringVect* pMirInf, srTDataMD* pExtraData) 
{
	if((pMirInf == 0) || (pMirInf->size() < 30)) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}
	if(pExtraData != 0) m_reflData = *pExtraData;

	const char* mirID = (*pMirInf)[1];

	m_halfDim1 = 0.5*atof((*pMirInf)[10]); //dimensions
	m_halfDim2 = 0.5*atof((*pMirInf)[11]);

	m_apertShape = 1; //1- rectangular, 2- elliptical 
	int iShape = atoi((*pMirInf)[12]);
	if((iShape > 0) && (iShape < 3)) m_apertShape = (char)iShape; //keep updated!

	m_vCenNorm.x = atof((*pMirInf)[16]); //central normal in the frame of incident beam
	m_vCenNorm.y = atof((*pMirInf)[17]);
	m_vCenNorm.z = atof((*pMirInf)[18]);
	if(m_vCenNorm.z == 0) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_ORIENT; return;}
	m_vCenNorm.Normalize();

	m_vCenTang.x = atof((*pMirInf)[19]);
	m_vCenTang.y = atof((*pMirInf)[20]);
	if((m_vCenTang.x == 0) && (m_vCenTang.y == 0)) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_ORIENT; return;}

	m_vCenTang.z = (-m_vCenNorm.x*m_vCenTang.x - m_vCenNorm.y*m_vCenTang.y)/m_vCenNorm.z;
	m_vCenTang.Normalize();

	TransvCenPoint.x = atof((*pMirInf)[23]);
	TransvCenPoint.y = atof((*pMirInf)[24]);

	m_propMeth = (char)atoi((*pMirInf)[26]);
	if((m_propMeth < 1) || (m_propMeth > 1)) //to keep updated
	{ ErrorCode = IMPROPER_OPTICAL_COMPONENT_SIM_METH; return;}

	//to program:
	int npT = atoi((*pMirInf)[28]); //number of points for representing the element in Tangential direction (for "thin" approx., etc.)
	int npS = atoi((*pMirInf)[29]); //number of points for representing the element in Sagital direction (for "thin" approx., etc.)

	//m_numPartsProp = (char)atoi((*pMirInf)[23]);
	//if(m_numPartsProp < 1) { ErrorCode = SRWL_INCORRECT_PARAM_FOR_WFR_PROP; return;}

	SetupNativeTransFromLocToBeamFrame(m_vCenNorm, m_vCenTang, TransvCenPoint);
	//FindElemExtentsAlongOptAxes(*(TransHndl.rep), m_vCenNorm, m_halfDim1, m_halfDim2, m_extAlongOptAxIn, m_extAlongOptAxOut); //virtual

	m_pRadAux = 0;
}

//*************************************************************************

srTMirror::srTMirror(const SRWLOptMir& srwlMir) 
{
	m_halfDim1 = 0.5*srwlMir.dt; //dimensions: tangential
	m_halfDim2 = 0.5*srwlMir.ds; //dimensions: sagital

	m_apertShape = 1; //1- rectangular, 2- elliptical 
	if(srwlMir.apShape == 'e') m_apertShape = 2;

	m_propMeth = srwlMir.meth;
	if((m_propMeth < 1) || (m_propMeth > 2)) //to keep updated
	{ ErrorCode = IMPROPER_OPTICAL_COMPONENT_SIM_METH; return;}

	m_npt = srwlMir.npt;
	m_nps = srwlMir.nps;

	m_treatInOut = srwlMir.treatInOut;
	//m_treatOut = srwlMir.treatOut;
	m_extAlongOptAxIn = srwlMir.extIn;
	m_extAlongOptAxOut = srwlMir.extOut;

	m_reflData.pData = 0; //OC12082018
	if(srwlMir.arRefl != 0)
	{
		m_reflData.pData = (char*)srwlMir.arRefl;
		m_reflData.DataType[0] = 'c';
		m_reflData.DataType[1] = 'd'; //?
		m_reflData.AmOfDims = 3;
		m_reflData.DimSizes[0] = srwlMir.reflNumPhEn;
		m_reflData.DimSizes[1] = srwlMir.reflNumAng;
		m_reflData.DimSizes[2] = srwlMir.reflNumComp;
		m_reflData.DimStartValues[0] = srwlMir.reflPhEnStart;
		m_reflData.DimStartValues[1] = srwlMir.reflAngStart;
		m_reflData.DimStartValues[2] = 1;
		
		m_reflData.DimSteps[0] = 0;
		//if(srwlMir.reflNumPhEn > 1) m_reflData.DimSteps[0] = (srwlMir.reflPhEnFin - srwlMir.reflPhEnStart)/(srwlMir.reflNumPhEn - 1);
		m_reflData.DimSteps[1] = 0;
		//if(srwlMir.reflNumAng > 1) m_reflData.DimSteps[1] = (srwlMir.reflAngFin - srwlMir.reflAngStart)/(srwlMir.reflNumAng - 1);
		m_reflData.DimSteps[2] = 0;

		if(strcmp(srwlMir.reflPhEnScaleType, "lin\0") == 0)
		{
			strcpy(m_reflData.DimScales[0], "lin\0");
			if(srwlMir.reflNumPhEn > 1) m_reflData.DimSteps[0] = (srwlMir.reflPhEnFin - srwlMir.reflPhEnStart)/(srwlMir.reflNumPhEn - 1);
		}
		else if(strcmp(srwlMir.reflPhEnScaleType, "log\0") == 0)
		{
			strcpy(m_reflData.DimScales[0], "log\0");
			if(srwlMir.reflNumPhEn > 1) m_reflData.DimSteps[0] = (log10(srwlMir.reflPhEnFin) - log10(srwlMir.reflPhEnStart))/(srwlMir.reflNumPhEn - 1);
		}
		//if(srwlMir.reflNumPhEn > 1) m_reflData.DimSteps[0] /= (srwlMir.reflNumPhEn - 1);
		if(strcmp(srwlMir.reflAngScaleType, "lin\0") == 0)
		{
			strcpy(m_reflData.DimScales[1], "lin\0");
			if(srwlMir.reflNumAng > 1) m_reflData.DimSteps[1] = (srwlMir.reflAngFin - srwlMir.reflAngStart)/(srwlMir.reflNumAng - 1);
		}
		else if(strcmp(srwlMir.reflAngScaleType, "log\0") == 0)
		{
			strcpy(m_reflData.DimScales[1], "log\0");
			if(srwlMir.reflNumAng > 1) m_reflData.DimSteps[1] = (log10(srwlMir.reflAngFin) - log10(srwlMir.reflAngStart))/(srwlMir.reflNumAng - 1);
		}
		//if(srwlMir.reflNumAng > 1) m_reflData.DimSteps[1] /= (srwlMir.reflNumAng - 1);

		strcpy(m_reflData.DimUnits[0], "eV");
		strcpy(m_reflData.DimUnits[1], "rad");
		m_reflData.DimUnits[2][0] = '\0';
		m_reflData.DataUnits[0] = '\0';
		m_reflData.DataName[0] = '\0';
		m_reflData.hState = 1;
	}

	m_vCenNorm.x = srwlMir.nvx; //central normal in the frame of incident beam
	m_vCenNorm.y = srwlMir.nvy;
	m_vCenNorm.z = srwlMir.nvz;
	//if(m_vCenNorm.z == 0) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_ORIENT; return;}
	if(m_vCenNorm.z == 0) { throw IMPROPER_OPTICAL_COMPONENT_ORIENT;}
	m_vCenNorm.Normalize();

	m_vCenTang.x = srwlMir.tvx;
	m_vCenTang.y = srwlMir.tvy;
	//if((m_vCenTang.x == 0) && (m_vCenTang.y == 0)) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_ORIENT; return;}
	if((m_vCenTang.x == 0) && (m_vCenTang.y == 0)) { throw IMPROPER_OPTICAL_COMPONENT_ORIENT;}
	m_vCenTang.z = (-m_vCenNorm.x*m_vCenTang.x - m_vCenNorm.y*m_vCenTang.y)/m_vCenNorm.z;
	m_vCenTang.Normalize();

	TransvCenPoint.x = srwlMir.x;
	TransvCenPoint.y = srwlMir.y;

	//This only calculates the transformation to the local frame
	SetupNativeTransFromLocToBeamFrame(m_vCenNorm, m_vCenTang, TransvCenPoint);
	//Other calculations (transformation of base vectors, finding extents of optical elements along optical axes, etc., will happen just before propagation)
	//FindElemExtentsAlongOptAxes(*(TransHndl.rep), m_vCenNorm, m_halfDim1, m_halfDim2, m_extAlongOptAxIn, m_extAlongOptAxOut); //virtual

	m_pRadAux = 0;
	m_wfrRadWasProp = false;

	m_grAuxAnamorphMagnH = 1.;
	m_grAuxAnamorphMagnV = 1.;
	m_grAuxElecFldAnamorphMagnFact = 1.;
}

//*************************************************************************

srTMirror* srTMirror::DefineMirror(srTStringVect* pMirInf, srTDataMD* pExtraData)
{
	//if((pMirInf == 0) || (pMirInf->size() < 24)) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return 0;}
	if((pMirInf == 0) || (pMirInf->size() < 3)) return 0;
	
	srTMirror *pOutMir = 0;

	//const char* mirID = (*pMirInf)[2];
	const char* mirID = (*pMirInf)[1];
	if(strcmp(mirID, "Toroid") == 0) pOutMir = new srTMirrorToroid(pMirInf, pExtraData);
	//else if(strcmp(mirID, "Paraboloid") == 0) return new srTMirrorToroid(pMirInf, pExtraData);
	//else return 0;

	pOutMir->m_isGrating = false;

	return pOutMir;
}

//*************************************************************************

srTMirror* srTMirror::DefineMirror(char* sType, void* pvData)
{
	if((sType == 0) || (pvData == 0)) throw IMPROPER_OPTICAL_COMPONENT_STRUCTURE;

	srTMirror *pOutMir = 0;

	if(strcmp(sType, "mirror: plane") == 0) pOutMir = new srTMirrorPlane(*((SRWLOptMirPl*)pvData));
	else if(strcmp(sType, "mirror: ellipsoid") == 0) pOutMir = new srTMirrorEllipsoid(*((SRWLOptMirEl*)pvData));
	else if(strcmp(sType, "mirror: toroid") == 0) pOutMir = new srTMirrorToroid(*((SRWLOptMirTor*)pvData));
	else if(strcmp(sType, "mirror: sphere") == 0) pOutMir = new srTMirrorSphere(*((SRWLOptMirSph*)pvData));
	else throw UNKNOWN_OPTICAL_ELEMENT;

	pOutMir->m_isGrating = false;

	return pOutMir;
}

//*************************************************************************

srTMirror* srTMirror::DefineGrating(char* sType, void* pvData)
{
	if((sType == 0) || (pvData == 0)) throw IMPROPER_OPTICAL_COMPONENT_STRUCTURE;

	SRWLOptG *pInGrat = (SRWLOptG*)pvData;
	char *sMirSubType = pInGrat->mirSubType;
	void *pvMirSub = pInGrat->mirSub;
	srTMirror *pOutMir = 0;

	if(strcmp(sMirSubType, "mirror: plane") == 0) pOutMir = new srTMirrorPlane(*((SRWLOptMirPl*)pvMirSub));
	else if(strcmp(sMirSubType, "mirror: ellipsoid") == 0) pOutMir = new srTMirrorEllipsoid(*((SRWLOptMirEl*)pvMirSub));
	else if(strcmp(sMirSubType, "mirror: toroid") == 0) pOutMir = new srTMirrorToroid(*((SRWLOptMirTor*)pvMirSub));
	else throw UNKNOWN_OPTICAL_ELEMENT;

	pOutMir->m_grM = pInGrat->m;
	pOutMir->m_grDen = pInGrat->grDen*1e+03; //[lines/mm] -> [lines/m]?
	pOutMir->m_grDen1 = pInGrat->grDen1*1e+06; //[lines/mm^2] -> [lines/m^2]?
	pOutMir->m_grDen2 = pInGrat->grDen2*1e+09; //[lines/mm^3] -> [lines/m^3]?
	pOutMir->m_grDen3 = pInGrat->grDen3*1e+12; //[lines/mm^4] -> [lines/m^4]?
	pOutMir->m_grDen4 = pInGrat->grDen4*1e+15; //[lines/mm^5] -> [lines/m^5]?
	pOutMir->m_grAng = pInGrat->grAng;

	pOutMir->m_grAuxCosAng = cos(pOutMir->m_grAng);
	pOutMir->m_grAuxSinAng = sin(pOutMir->m_grAng);

	pOutMir->m_isGrating = true;

	return pOutMir;
}

//*************************************************************************

void srTMirror::SetupNativeTransFromLocToBeamFrame(TVector3d& vCenNorm, TVector3d& vCenTang, TVector2d& vCenP2d)
{//In the Local frame, tangential direction is X, sagital Y
	TVector3d mRow1(vCenTang.x, vCenNorm.y*vCenTang.z - vCenNorm.z*vCenTang.y, vCenNorm.x);
	TVector3d mRow2(vCenTang.y, vCenNorm.z*vCenTang.x - vCenNorm.x*vCenTang.z, vCenNorm.y);
	TVector3d mRow3(vCenTang.z, vCenNorm.x*vCenTang.y - vCenNorm.y*vCenTang.x, vCenNorm.z);
	TMatrix3d M(mRow1, mRow2, mRow3);
	TVector3d vCen(vCenP2d.x, vCenP2d.y, 0);

	gmTrans *pTrans = new gmTrans(M, vCen);
	TransHndl = srTransHndl(pTrans);

/**
	TVector3d vUz(0, 0, 1), vUx(1, 0, 0), vUy(0, 1, 0);
	m_vInLoc = pTrans->TrBiPoint_inv(vUz); //direction of input optical axis in the local frame of opt. elem.

	m_vOutLoc = m_vInLoc - (2.*(m_vInLoc*vUz))*vUz; //direction of output optical axis in the local frame of opt. elem.
	//To modify the above: find m_vOutLocintersection with surface


	//Defining the basis vectors of the output beam frame.
	//The Beam frame should stay "right-handed" even after the reflection;
	//therefore the new basis vectors should be obtained by rotation from the previous basis vectors.
	//The rotation should be around the axis perpendicular to the plane of incidence (i.e. plane of reflection)
	//and the rotation angle is the angle between the incident and the reflected central beams (i.e. optical axes before and aftre the reflection).

	const double relTolZeroVect = 1.e-10;
	const double relTolZeroVectE2 = relTolZeroVect*relTolZeroVect;

	//double absCenNormTrE2 = vCenNorm.x*vCenNorm.x + vCenNorm.y*vCenNorm.y;
	//double absCenNormE2 = absCenNormTrE2 + vCenNorm.z*vCenNorm.z;
	//if(absCenNormTrE2 < absCenNormE2*relTolZeroVectE2)

	TVector3d vDifOutIn = m_vOutLoc - m_vInLoc;
	double absDifE2 = vDifOutIn.AmpE2();
	//Special case: central normal parallel to the input (and output) optical axis
	if(absDifE2 < relTolZeroVectE2)
	{
		m_vHorOutIn.x = -1.; m_vHorOutIn.y = m_vHorOutIn.z = 0.; //= -vUx;
		m_vVerOutIn.y = 1.; m_vHorOutIn.x = m_vHorOutIn.z = 0.; //= vUy;
	}

	//General case: rotation about this axis:

	TVector3d vZero(0,0,0), vRotAxis(-vCenNorm.y, vCenNorm.x, 0.); //= vCenNorm^vUz
	double rotAng = acos(m_vInLoc*m_vOutLoc);

	gmTrans auxRot;
	auxRot.SetupRotation(vZero, vRotAxis, rotAng);
	m_vHorOutIn = auxRot.TrBiPoint(vUx); //output horizontal vector in the frame of input beam
	m_vVerOutIn = auxRot.TrBiPoint(vUy); //output vertical vector in the frame of input beam

	//m_vHorOutIn = vUx - (2.*(vUx*vCenNorm))*vCenNorm; //output horizontal vector in the frame of input beam
	//m_vVerOutIn = vUy - (2.*(vUy*vCenNorm))*vCenNorm; 
**/
}

//*************************************************************************

int srTMirror::FindBasisVectorTransAndExtents()
{//Setting up auxiliary vectors before the propagation: m_vInLoc, m_vOutLoc, m_vHorOutIn, m_vVerOutIn;
 //To be called after the native transformation to the local frame has been already set up!
 //Assumes TransHndl, m_ParPrecWfrPropag, m_grAux... have been set in advance

	gmTrans *pTrans = TransHndl.rep;

	TVector3d vUz(0, 0, 1), vUx(1, 0, 0), vUy(0, 1, 0), vZero(0, 0, 0);
	//It is assumed that the optical axis in the frame of input beam is defined by point {0,0,0} and vector {0,0,1}
	m_vInLoc = pTrans->TrBiPoint_inv(vUz); //direction of input optical axis in the Local frame of opt. elem.
	TVector3d vCenPtLoc = pTrans->TrPoint_inv(vZero); //point through which the iput optical axis passes in the local frame
	TVector3d vIntersPtLocFr, vNormAtIntersPtLoc;
	if(!FindRayIntersectWithSurfInLocFrame(vCenPtLoc, m_vInLoc, vIntersPtLocFr, &vNormAtIntersPtLoc)) return FAILED_DETERMINE_OPTICAL_AXIS;
	m_vInLoc.Normalize();

	//OC021213
	//bool OutFrameBaseVectAreDefined = false;
	if((m_ParPrecWfrPropag.vLxOut != 0) || (m_ParPrecWfrPropag.vLyOut != 0) || (m_ParPrecWfrPropag.vLzOut != 0)) 
	{//Defines m_vOutLoc, m_vHorOutIn, m_vVerOutIn in the case if these vectors are defined in the propagation parameters
		TVector3d vOutIn(m_ParPrecWfrPropag.vLxOut, m_ParPrecWfrPropag.vLyOut, m_ParPrecWfrPropag.vLzOut);
		vOutIn.Normalize();

		TVector3d vTestHorOutIn(m_ParPrecWfrPropag.vHxOut, m_ParPrecWfrPropag.vHyOut, 0.);
		if(vOutIn.z == 0)
		{
			if((vTestHorOutIn.x != 0) || (vTestHorOutIn.y != 0))
			{
				vTestHorOutIn.Normalize();
				const double relTol = 1e-12;
				double testScalProd = ::fabs(vOutIn*vTestHorOutIn);
				if(testScalProd > relTol) return FAILED_DETERMINE_OPTICAL_AXIS;
			}
			else vTestHorOutIn.z = 1.;
		}
		else
		{
			vTestHorOutIn.z = (-vOutIn.x*vTestHorOutIn.x - vOutIn.y*vTestHorOutIn.y)/vOutIn.z;
			vTestHorOutIn.Normalize();
		}

		m_vHorOutIn = vTestHorOutIn;
		m_vVerOutIn = vOutIn^vTestHorOutIn;
		m_vOutLoc = pTrans->TrBiPoint_inv(vOutIn);
	}
	else
	{//Defines m_vOutLoc, m_vHorOutIn, m_vVerOutIn by default

		//m_vOutLoc = m_vInLoc - ((2.*(m_vInLoc*vNormAtIntersPtLoc))*vNormAtIntersPtLoc); //direction of output optical axis in the local frame of opt. elem.
		//m_vOutLoc.Normalize();

		if(m_isGrating) 
		{
			//Tangential vector perpendicular to grooves in the center
			TVector3d vTang(vNormAtIntersPtLoc.z*m_grAuxCosAng, vNormAtIntersPtLoc.z*m_grAuxSinAng, -(vNormAtIntersPtLoc.x*m_grAuxCosAng + vNormAtIntersPtLoc.y*m_grAuxSinAng));
			vTang.Normalize();
			double xGr = vIntersPtLocFr.x;
			double locGrDen = m_grDen + xGr*(xGr*(xGr*(xGr*m_grDen4 + m_grDen3) + m_grDen2) + m_grDen1); //Calculate local Groove Density
			//vTang *= (locGrDen*m_grM/(806554.3835*m_grAuxEphAvg));
			vTang *= (-locGrDen*m_grM/(806554.3835*m_grAuxEphAvg)); //OC280214

			TVector3d vInLocTang = m_vInLoc - ((m_vInLoc*vNormAtIntersPtLoc)*vNormAtIntersPtLoc);
			TVector3d vOutLocTang = vInLocTang + vTang;
			double absE2_vOutLocTang = vOutLocTang.AmpE2();
			double abs_vOutLocNorm = sqrt(::fabs(1. - absE2_vOutLocTang));
			m_vOutLoc = vOutLocTang + (abs_vOutLocNorm*vNormAtIntersPtLoc);

			//m_vOutLoc += vTang; //Check the sign!
			//m_vOutLoc.Normalize(); //required here?
		}
		else
		{
			m_vOutLoc = m_vInLoc - ((2.*(m_vInLoc*vNormAtIntersPtLoc))*vNormAtIntersPtLoc); //direction of output optical axis in the local frame of opt. elem.
			m_vOutLoc.Normalize();
		}

		//Defining coordinates of basis vectors of the output beam frame in the input beam frame (m_vHorOutIn, m_vVerOutIn).
		//The Beam frame should stay "right-handed" even after the reflection;
		//therefore the new basis vectors should be obtained by rotation from the previous basis vectors.
		//The rotation should be around the axis perpendicular to the plane of incidence (i.e. plane of reflection)
		//and the rotation angle is the angle between the incident and the reflected central beams (i.e. optical axes before and aftre the reflection).

		//Checking for Special case: central normal parallel to the input (and output) optical axis
		const double relTolZeroVect = 1.e-10;
		double absRotAng = acos(m_vInLoc*m_vOutLoc);
		if(absRotAng < relTolZeroVect)
		{//Special case: central normal parallel to the input (and output) optical axis
			m_vHorOutIn.x = -1.; m_vHorOutIn.y = m_vHorOutIn.z = 0.; //= -vUx;
			m_vVerOutIn.y = 1.; m_vHorOutIn.x = m_vHorOutIn.z = 0.; //= vUy;
		}

		//General case: rotation about vRotAxis:
		TVector3d vRotAxis = m_vInLoc^m_vOutLoc; //rotation axis in the Local opt. elem. frame
		vRotAxis = pTrans->TrBiPoint(vRotAxis); //now in the frame of input beam

		gmTrans auxRot;
		auxRot.SetupRotation(vZero, vRotAxis, absRotAng);
		m_vHorOutIn = auxRot.TrBiPoint(vUx); //output horizontal vector in the frame of input beam
		m_vVerOutIn = auxRot.TrBiPoint(vUy); //output vertical vector in the frame of input beam
	}

	if(m_isGrating)
	{//Estimating Anamorphic Magnification in different planes

		TVector3d vNormAtIntersPtIn = pTrans->TrBiPoint(vNormAtIntersPtLoc); //central normal vector in the frame of input beam
		TVector3d vOutIn = pTrans->TrBiPoint(m_vOutLoc);

		double grFocStrCoef = (m_grDen1*m_grM/(806554.3835*m_grAuxEphAvg));
		TVector3d vIn(0,0,1);
		TVector3d vSagIn = vOutIn^vIn;
		vSagIn.Normalize();

		//In the Horizontal Plane:
		TVector3d vNormAtIntersPtInHorPl(vNormAtIntersPtIn.x, 0., vNormAtIntersPtIn.z);
		vNormAtIntersPtInHorPl.Normalize();
		TVector3d vOutInHorPl(vOutIn.x, 0, vOutIn.z), vInInHorPl(0, 0, 1.);
		vOutInHorPl.Normalize();
		double absCosOutHorPl = ::fabs(vNormAtIntersPtInHorPl*vOutInHorPl);
		double absCosInHorPl = ::fabs(vNormAtIntersPtInHorPl*vInInHorPl);
		m_grAuxAnamorphMagnH = 1.;
		if(absCosInHorPl > 1.e-10) 
		{
			m_grAuxAnamorphMagnH = absCosOutHorPl/absCosInHorPl;
			if(grFocStrCoef != 0) 
			{
				TVector3d vHorOutInPr = m_vHorOutIn - ((m_vHorOutIn*vSagIn)*vSagIn);
				double absHorPr = vHorOutInPr.Abs();
				//double grFocStr = absHorPr*grFocStrCoef/absCosOutHorPl; //focal strength of variable-groove density grating in hor. plane
				double grFocStr = absHorPr*grFocStrCoef/(absCosOutHorPl*absCosOutHorPl); //OC280214, to check
				double totForStr = grFocStr + 1./FocDistX;
				FocDistX = 1./totForStr;
			}
		}

		//In the Vertical Plane:
		TVector3d vNormAtIntersPtInVerPl(0., vNormAtIntersPtIn.y, vNormAtIntersPtIn.z);
		vNormAtIntersPtInVerPl.Normalize();
		TVector3d vOutInVerPl(0, vOutIn.y, vOutIn.z), vInInVerPl(0, 0, 1.);
		vOutInVerPl.Normalize();
		double absCosOutVerPl = ::fabs(vNormAtIntersPtInVerPl*vOutInVerPl);
		double absCosInVerPl = ::fabs(vNormAtIntersPtInVerPl*vInInVerPl);
		m_grAuxAnamorphMagnV = 1.;
		if(absCosInVerPl > 1.e-10) 
		{
			m_grAuxAnamorphMagnV = absCosOutVerPl/absCosInVerPl;
			if(grFocStrCoef != 0) 
			{
				TVector3d vVerOutInPr = m_vVerOutIn - ((m_vVerOutIn*vSagIn)*vSagIn);
				double absVerPr = vVerOutInPr.Abs();
				//double grFocStr = absVerPr*grFocStrCoef/absCosOutVerPl; //focal strength of variable-groove density grating in hor. plane
				double grFocStr = absVerPr*grFocStrCoef/(absCosOutVerPl*absCosOutVerPl); //OC280214, to check
				double totForStr = grFocStr + 1./FocDistZ;
				FocDistZ = 1./totForStr;
			}
		}

		m_grAuxElecFldAnamorphMagnFact = 1./sqrt(m_grAuxAnamorphMagnH*m_grAuxAnamorphMagnV); //to check
	}

	if((m_extAlongOptAxIn == 0.) && (m_extAlongOptAxOut == 0.))
	{
		//Calculate "extents": m_extAlongOptAxIn, m_extAlongOptAxOut, using other member variables (which are assumed to be already defined)
		//Mirror cormers in local frame:
		TVector3d r1(-m_halfDim1, -m_halfDim2, 0), r2(m_halfDim1, -m_halfDim2, 0), r3(-m_halfDim1, m_halfDim2, 0), r4(m_halfDim1, m_halfDim2, 0); 
		//r1 = pTrans->TrBiPoint(r1);
		//r2 = pTrans->TrBiPoint(r2);
		//r3 = pTrans->TrBiPoint(r3);
		//r4 = pTrans->TrBiPoint(r4);
		//Mirror cormers in Input beam frame:
		r1 = pTrans->TrPoint(r1);
		r2 = pTrans->TrPoint(r2);
		r3 = pTrans->TrPoint(r3);
		r4 = pTrans->TrPoint(r4);
		//Longitudinal positions of Mirror cormers in the Input beam frame:
		double arLongCoordIn[] = {r1.z, r2.z, r3.z, r4.z};

		TVector3d m_vLongOutIn = m_vHorOutIn^m_vVerOutIn;
		double arLongCoordOut[] = {r1*m_vLongOutIn, r2*m_vLongOutIn, r3*m_vLongOutIn, r4*m_vLongOutIn};

		double *t_arLongCoordIn = arLongCoordIn, *t_arLongCoordOut = arLongCoordOut;
		double extMin = *(t_arLongCoordIn++), extMax = *(t_arLongCoordOut++);
		for(int i=0; i<3; i++)
		{
			if(extMin > *t_arLongCoordIn) extMin = *t_arLongCoordIn;
			if(extMax < *t_arLongCoordOut) extMax = *t_arLongCoordOut;
			t_arLongCoordIn++; t_arLongCoordOut++;
		}

		m_extAlongOptAxIn = ::fabs(extMin);
		m_extAlongOptAxOut = extMax;
	}

	return 0;
}

//*************************************************************************

void srTMirror::FindElemExtentsAlongOptAxes(gmTrans& trfMir, TVector3d& vCenNorm, double halfDim1, double halfDim2, double& extIn, double& extOut)
{
	TVector3d r1(-halfDim1, -halfDim2, 0), r2(halfDim1, -halfDim2, 0), r3(-halfDim1, halfDim2, 0), r4(halfDim1, halfDim2, 0);
	//TVector3d r1(-halfDim1 + TransvCenPoint.x, -halfDim2 + TransvCenPoint.y, 0);
	//TVector3d r2(halfDim1 + TransvCenPoint.x, -halfDim2 + TransvCenPoint.y, 0);
	//TVector3d r3(-halfDim1 + TransvCenPoint.x, halfDim2 + TransvCenPoint.y, 0);
	//TVEctor3d r4(halfDim1 + TransvCenPoint.x, halfDim2 + TransvCenPoint.y, 0);

	gmTrans *pTrans = TransHndl.rep;
	//r1 = pTrans->TrBiPoint(r1);
	//r2 = pTrans->TrBiPoint(r2);
	//r3 = pTrans->TrBiPoint(r3);
	//r4 = pTrans->TrBiPoint(r4);
	r1 = pTrans->TrPoint(r1);
	r2 = pTrans->TrPoint(r2);
	r3 = pTrans->TrPoint(r3);
	r4 = pTrans->TrPoint(r4);
	double arLongCoordIn[] = {r1.z, r2.z, r3.z, r4.z};

	TVector3d vOptAxis(0, 0, 1);
	vOptAxis -= (2.*(vOptAxis*vCenNorm))*vCenNorm;
	double arLongCoordOut[] = {r1*vOptAxis, r2*vOptAxis, r3*vOptAxis, r4*vOptAxis};

	double *t_arLongCoordIn = arLongCoordIn, *t_arLongCoordOut = arLongCoordOut;
	double extMin = *(t_arLongCoordIn++), extMax = *(t_arLongCoordOut++);
	for(int i=0; i<3; i++)
	{
		if(extMin > *t_arLongCoordIn) extMin = *t_arLongCoordIn;
		if(extMax < *t_arLongCoordOut) extMax = *t_arLongCoordOut;
		t_arLongCoordIn++; t_arLongCoordOut++;
	}

	extIn = ::fabs(extMin);
	extOut = extMax;
}

//*************************************************************************

int srTMirror::WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax, double dxMax, double dzMax) //OC20082018
//int srTMirror::WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax)
//int srTMirror::WfrInterpolOnOrigGrid2(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, long* arIndRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax)
{//OC18032016
	if((pWfr == 0) || (arRayTrCoord == 0) || ((arEX == 0) && (arEZ == 0))) return FAILED_INTERPOL_ELEC_FLD;
	//if((pWfr == 0) || (arRayTrCoord == 0) || (arOptPathDif == 0) || ((arEX == 0) && (arEZ == 0))) return FAILED_INTERPOL_ELEC_FLD;

	bool isCoordRepres = (pWfr->Pres == 0);
	bool isFreqRepres = (pWfr->PresT == 0);
	bool waveFrontTermWasTreated = false;

	if(isFreqRepres && isCoordRepres && WaveFrontTermCanBeTreated(*pWfr, false)) //OC300914
	{
		float *pPrevBaseRadX = pWfr->pBaseRadX; 
		float *pPrevBaseRadZ = pWfr->pBaseRadZ;
		pWfr->pBaseRadX = arEX;
		pWfr->pBaseRadZ = arEZ;
		
		//if(!((fabs(pWfr->RobsX + 0.99) < 0.1) && (fabs(pWfr->RobsZ + 0.99) < 0.1)))
		//{//OCTEST
		
		TreatStronglyOscillatingTermIrregMesh(*pWfr, arRayTrCoord, xMin, xMax, zMin, zMax, 'r');
		
		//}//OCTEST
		
		pWfr->pBaseRadX = pPrevBaseRadX;
		pWfr->pBaseRadZ = pPrevBaseRadZ;
		waveFrontTermWasTreated = true;
	}

	float *t_ExRes = pWfr->pBaseRadX;
	float *t_EzRes = pWfr->pBaseRadZ;

	//long HalfPerX = (pWfr->ne);
	//long PerX = HalfPerX << 1;
	//long HalfPerZ = HalfPerX*(pWfr->nx);
	//long PerZ = PerX*(pWfr->nx);
	//long half_nTot = HalfPerZ*(pWfr->nz);
	//long nTot = PerZ*(pWfr->nz);

	long long HalfPerX = (pWfr->ne);
	long long PerX = HalfPerX << 1;
	long long HalfPerZ = HalfPerX*(pWfr->nx);
	long long PerZ = PerX*(pWfr->nx);
	long long half_nTot = HalfPerZ*(pWfr->nz);
	long long nTot = PerZ*(pWfr->nz);

	//OCTEST
	//if(fabs(pWfr->RobsZ + 18.079) < 0.1)
	//{
		//float *t_arEX = arEX, *t_arEZ = arEZ;
		//for(int k=0; k<nTot; k++)
		//{
		//	*(t_ExRes++) = *(t_arEX++);
		//	*(t_EzRes++) = *(t_arEZ++);
		//}
		//return 0;
		//int aha = 1;
	//}
	//END OCTEST
	
	double arRelCoordXZ[6]; //Coordinates of points to be used for interpolation
	double arReEx[4], arImEx[4], arReEz[4], arImEz[4], arIx[4], arIz[4]; //Aux. arrays to be used for interpolation
	double dx, dz, ReE, ImE;

	//const int maxSearchRad = 3; //To tune
	//OC11082018 (the above didn't allow to find indCloseRayTrCoord for the case of Grating in Example #12)
	int maxSearchRad = 20; //10; //100; //1000; //To tune
	int halfNx_mi_2 = (pWfr->nx - 1) >> 1;
	if(maxSearchRad > halfNx_mi_2) maxSearchRad = halfNx_mi_2;
	int halfNz_mi_2 = (pWfr->nz - 1) >> 1;
	if(maxSearchRad > halfNz_mi_2) maxSearchRad = halfNz_mi_2;

	const double interpSafeFact = 2.5; //OC20082018

	long ixMin = (long)((xMin - (pWfr->xStart))/(pWfr->xStep));
	long ixMax = (long)((xMax - (pWfr->xStart))/(pWfr->xStep));
	long izMin = (long)((zMin - (pWfr->zStart))/(pWfr->zStep));
	long izMax = (long)((zMax - (pWfr->zStart))/(pWfr->zStep));

	double z = pWfr->zStart;
	for(long iz=0; iz<(pWfr->nz); iz++)
	{
		bool pointIsWithinVertLim = (zMin < z) && (z < zMax);

		//long izHalfPerZ = iz*HalfPerZ;
		long long izHalfPerZ = iz*HalfPerZ;
		double x = pWfr->xStart;
		for(long ix=0; ix<(pWfr->nx); ix++)
		{
			bool pointIsWithinTransvLim = (xMin < x) && (x < xMax) && pointIsWithinVertLim;

			//long izHalfPerZ_p_ixHalfPerX = izHalfPerZ + ix*HalfPerX;
			long long izHalfPerZ_p_ixHalfPerX = izHalfPerZ + ix*HalfPerX;
			for(long ie=0; ie<(pWfr->ne); ie++)
			{
				*t_ExRes = 0.; *(t_ExRes+1) = 0.;
				*t_EzRes = 0.; *(t_EzRes+1) = 0.;

						//OCTEST
						//if(m_isGrating && (ix==207) && (iz==1425))
						//if(m_isGrating && (ix==287) && (iz==1427))
						//{
						//	int aha = 1;
						//}
						//END OCTEST

				//long indCloseRayTrCoord = arIndRayTrCoord[izHalfPerZ_p_ixHalfPerX + ie];
				long long indCloseRayTrCoord = arIndRayTrCoord[izHalfPerZ_p_ixHalfPerX + ie];

				if((indCloseRayTrCoord < 0) && pointIsWithinTransvLim) //OC19082018
				//if((indCloseRayTrCoord < 0) && ((xMin < x) && (x < xMax) && (zMin < z) && (z < zMax)))
				{//Try to find nearest non-negative value of the index

					bool thereAreDataXL=false, thereAreDataXU=false;
					bool thereAreDataYL=false, thereAreDataYU=false;
					long long indFirstFound=-1;

					for(int ic=1; ic<=maxSearchRad; ic++)
					{
						int two_ic = ic << 1;

						if((indFirstFound < 0) || (!thereAreDataYL) || (!thereAreDataYU)) //OC12082018
						//if(indFirstFound < 0) //OC12082018
						{
							for(int icz=-ic; icz<=ic; icz+=two_ic)
							{
								long jz = iz + icz;
								if((jz < izMin) || (jz > izMax)) continue;

								//long jzHalfPerZ_p_ie = jz*HalfPerZ + ie;
								long long jzHalfPerZ_p_ie = jz*HalfPerZ + ie;
								for(int icx=-ic; icx<=ic; icx++)
								{
									//long jx = ix + icx;
									long icxCor = icx + ic; //OC21082018
									if(icxCor > 0)
									{
										int icxCor_d_2 = icxCor >> 1;
										if((icxCor_d_2 << 1) != icxCor) icxCor = -(icxCor_d_2 + 1);
										else icxCor = icxCor_d_2;
									}
									long jx = ix + icxCor; //OC21082018
									if((jx < ixMin) || (jx > ixMax)) continue;

									//long testOfst = jzHalfPerZ_p_ie + jx*HalfPerX;
									long long testOfst = jzHalfPerZ_p_ie + jx*HalfPerX;
									//long testOfst = jzHalfPerZ_p_ie + ix*HalfPerX;
									if((testOfst < 0) || (testOfst >= half_nTot)) continue;

									//indCloseRayTrCoord = arIndRayTrCoord[testOfst];
									//if(indCloseRayTrCoord >= 0) 
									//long indCur = arIndRayTrCoord[testOfst];
									long long indCur = arIndRayTrCoord[testOfst];
									if(indCur >= 0)
									{
										if(indFirstFound < 0) 
										{
											indFirstFound = indCur;
										}
										//else
										//{//OC12082018
										//	if(icz >= 0) thereAreDataYU = true;
										//}
										//if(icz < 0) thereAreDataYL = true;

										if(icz < 0) thereAreDataYL = true;
										else thereAreDataYU = true;
										break;
									}
								}
								//if(indCloseRayTrCoord >= 0) break;
								//OC11082018
								if((indFirstFound >= 0) && thereAreDataYL && thereAreDataYU) 
								{
									break;
								}
								//if(indFirstFound >= 0)
								//{
								//	//indCloseRayTrCoord = indFirstFound;
								//	break;
								//}
							}
						}
						//if(indCloseRayTrCoord >= 0) break;
						if((indFirstFound >= 0) && thereAreDataXL && thereAreDataXU && thereAreDataYL && thereAreDataYU)
						//OC12082018: problems with Example #12 showed "negative" impact of thereAreDataXL, thereAreDataXU, thereAreDataYL, thereAreDataYU
						//if(indFirstFound >= 0)
						{
							indCloseRayTrCoord = indFirstFound;
							break;
						}

						if((indFirstFound < 0) || (!thereAreDataXL) || (!thereAreDataXU)) //OC12082018
						//if(indFirstFound < 0) //OC12082018
						{
							for(int icx=-ic; icx<=ic; icx+=two_ic)
							{
								long jx = ix + icx;
								if((jx < ixMin) || (jx > ixMax)) continue;

								//long jxHalfPerX_p_ie = jx*HalfPerX + ie;
								long long jxHalfPerX_p_ie = jx*HalfPerX + ie;
								for(int icz=-ic+1; icz<=(ic-1); icz++)
								{
									//long jz = iz + icz;
									long iczCor = icz + ic - 1; //OC21082018
									if(iczCor > 0)
									{
										int iczCor_d_2 = iczCor >> 1;
										if((iczCor_d_2 << 1) != iczCor) iczCor = -(iczCor_d_2 + 1);
										else iczCor = iczCor_d_2;
									}
									long jz = iz + iczCor; //OC21082018
									if((jz < izMin) || (jz > izMax)) continue;

									//long testOfst = jxHalfPerX_p_ie + jz*HalfPerZ;
									long long testOfst = jxHalfPerX_p_ie + jz*HalfPerZ;
									//long testOfst = jxHalfPerX_p_ie + iz*HalfPerZ;
									if((testOfst < 0) || (testOfst >= half_nTot)) continue;

									//indCloseRayTrCoord = arIndRayTrCoord[testOfst];
									//if(indCloseRayTrCoord >= 0) break;
									//long indCur = arIndRayTrCoord[testOfst];
									long long indCur = arIndRayTrCoord[testOfst];
									if(indCur >= 0)
									{
										if(indFirstFound < 0) 
										{
											indFirstFound = indCur;
										}
										//else
										//{//OC12082018
										//	if(icx >= 0) thereAreDataXU = true;
										//}
										//if(icx < 0) thereAreDataXL = true;

										if(icx < 0) thereAreDataXL = true;
										else thereAreDataXU = true;
										break;
									}
								}
								//if(indCloseRayTrCoord >= 0) break;
								//OC11082018
								if((indFirstFound >= 0) && thereAreDataXL && thereAreDataXU) 
								{
									break;
								}
								//if(indFirstFound >= 0)
								//{
								//	//indCloseRayTrCoord = indFirstFound;
								//	break;
								//}
							}
						}
					
						//if(indCloseRayTrCoord >= 0) break;
						if((indFirstFound >= 0) && thereAreDataXL && thereAreDataXU && thereAreDataYL && thereAreDataYU)
						//OC12082018: problems with Example #12 showed "negative" impact of thereAreDataXL, thereAreDataXU, thereAreDataYL, thereAreDataYU
						//if(indFirstFound >= 0)
						{
							indCloseRayTrCoord = indFirstFound;
							break;
						}

						//if(((fabs(ic*(pWfr->xStep)) > dxMax*interpSafeFact) && ((!thereAreDataXL) || (!thereAreDataXU))) ||
						//   ((fabs(ic*(pWfr->zStep)) > dzMax*interpSafeFact) && ((!thereAreDataYL) || (!thereAreDataYU)))) //OC20082018
						if((indFirstFound < 0) && ((fabs(ic*(pWfr->xStep)) > dxMax*interpSafeFact) || (fabs(ic*(pWfr->zStep)) > dzMax*interpSafeFact))) //OC20082018
						{
							break; //Stop search for a good point over a too large range
						}
					}
				}

				//if(indCloseRayTrCoord >= 0) 
				//if((indCloseRayTrCoord >= 0) && ((xMin <= x) && (x <= xMax) && (zMin <= z) && (z <= zMax))) //OC03072017 (trying to avoid having an aventual "nan" in the resulting field data)
				if((indCloseRayTrCoord >= 0) && pointIsWithinTransvLim) //OC19082018
				{
					double *pRayTrCoord = arRayTrCoord + indCloseRayTrCoord;
					double xCloseRayTr = *pRayTrCoord;
					double zCloseRayTr = *(pRayTrCoord + 1);

					long iz0 = (long)(indCloseRayTrCoord/PerZ);
					long ix0 = (long)((indCloseRayTrCoord - iz0*PerZ)/PerX);
					long ie0 = (long)((indCloseRayTrCoord - iz0*PerZ - ix0*PerX) >> 1);
					long ie0PerE = ie0 << 1;

					//Determine 3 other 2D points "surrounding" the point (x, y)

					long izm1 = iz0 - 1, izp1 = iz0 + 1;
					if(izm1 < 0) 
					{
						izm1 = 0; iz0 = 1; izp1 = 2;
					}
					else if(izp1 >= (pWfr->nz))
					{
						izm1 = (pWfr->nz) - 3; iz0 = izm1 + 1; izp1 = iz0 + 1;
					}
			
					long ixm1 = ix0 - 1, ixp1 = ix0 + 1;
					if(ixm1 < 0) 
					{
						ixm1 = 0; ix0 = 1; ixp1 = 2;
					}
					else if(ixp1 >= (pWfr->nx))
					{
						ixm1 = (pWfr->nx) - 3; ix0 = ixm1 + 1; ixp1 = ix0 + 1;
					}

					//long ixm1PerX_p_ie0PerE = ixm1*PerX + ie0PerE;
					//long izm1PerZ = izm1*PerZ;
					//long ofst_ixm1_izm1 = izm1PerZ + ixm1PerX_p_ie0PerE;
					long long ixm1PerX_p_ie0PerE = ixm1*PerX + ie0PerE;
					long long izm1PerZ = izm1*PerZ;
					long long ofst_ixm1_izm1 = izm1PerZ + ixm1PerX_p_ie0PerE;
					double xCoord_ixm1_izm1 = arRayTrCoord[ofst_ixm1_izm1];
					double zCoord_ixm1_izm1 = arRayTrCoord[ofst_ixm1_izm1 + 1];

					//long iz0PerZ = iz0*PerZ;
					//long ofst_ixm1_iz0 = iz0PerZ + ixm1PerX_p_ie0PerE;
					long long iz0PerZ = iz0*PerZ;
					long long ofst_ixm1_iz0 = iz0PerZ + ixm1PerX_p_ie0PerE;
					double xCoord_ixm1_iz0 = arRayTrCoord[ofst_ixm1_iz0];
					double zCoord_ixm1_iz0 = arRayTrCoord[ofst_ixm1_iz0 + 1];

					//long izp1PerZ = izp1*PerZ;
					//long ofst_ixm1_izp1 = izp1PerZ + ixm1PerX_p_ie0PerE;
					long long izp1PerZ = izp1*PerZ;
					long long ofst_ixm1_izp1 = izp1PerZ + ixm1PerX_p_ie0PerE;
					double xCoord_ixm1_izp1 = arRayTrCoord[ofst_ixm1_izp1];
					double zCoord_ixm1_izp1 = arRayTrCoord[ofst_ixm1_izp1 + 1];

					//long ix0PerX_p_ie0PerE = ix0*PerX + ie0PerE;
					//long ofst_ix0_izm1 = izm1PerZ + ix0PerX_p_ie0PerE;
					long long ix0PerX_p_ie0PerE = ix0*PerX + ie0PerE;
					long long ofst_ix0_izm1 = izm1PerZ + ix0PerX_p_ie0PerE;
					double xCoord_ix0_izm1 = arRayTrCoord[ofst_ix0_izm1];
					double zCoord_ix0_izm1 = arRayTrCoord[ofst_ix0_izm1 + 1];

					//long ofst_ix0_iz0 = iz0PerZ + ix0PerX_p_ie0PerE;
					long long ofst_ix0_iz0 = iz0PerZ + ix0PerX_p_ie0PerE;
					double xCoord_ix0_iz0 = arRayTrCoord[ofst_ix0_iz0];
					double zCoord_ix0_iz0 = arRayTrCoord[ofst_ix0_iz0 + 1];

					//long ofst_ix0_izp1 = izp1PerZ + ix0PerX_p_ie0PerE;
					long long ofst_ix0_izp1 = izp1PerZ + ix0PerX_p_ie0PerE;
					double xCoord_ix0_izp1 = arRayTrCoord[ofst_ix0_izp1];
					double zCoord_ix0_izp1 = arRayTrCoord[ofst_ix0_izp1 + 1];

					//long ixp1PerX_p_ie0PerE = ixp1*PerX + ie0PerE;
					//long ofst_ixp1_izm1 = izm1PerZ + ixp1PerX_p_ie0PerE;
					long long ixp1PerX_p_ie0PerE = ixp1*PerX + ie0PerE;
					long long ofst_ixp1_izm1 = izm1PerZ + ixp1PerX_p_ie0PerE;
					double xCoord_ixp1_izm1 = arRayTrCoord[ofst_ixp1_izm1];
					double zCoord_ixp1_izm1 = arRayTrCoord[ofst_ixp1_izm1 + 1];

					//long ofst_ixp1_iz0 = iz0PerZ + ixp1PerX_p_ie0PerE;
					long long ofst_ixp1_iz0 = iz0PerZ + ixp1PerX_p_ie0PerE;
					double xCoord_ixp1_iz0 = arRayTrCoord[ofst_ixp1_iz0];
					double zCoord_ixp1_iz0 = arRayTrCoord[ofst_ixp1_iz0 + 1];

					//long ofst_ixp1_izp1 = izp1PerZ + ixp1PerX_p_ie0PerE;
					long long ofst_ixp1_izp1 = izp1PerZ + ixp1PerX_p_ie0PerE;
					double xCoord_ixp1_izp1 = arRayTrCoord[ofst_ixp1_izp1];
					double zCoord_ixp1_izp1 = arRayTrCoord[ofst_ixp1_izp1 + 1];

					//Aux. coordinates of centers of quadrants to be used for choosing the most appropriate quadrant for teh interpolation
					double xq00 = 0.25*(xCoord_ixm1_izm1 + xCoord_ix0_izm1 + xCoord_ixm1_iz0 + xCoord_ix0_iz0);
					double zq00 = 0.25*(zCoord_ixm1_izm1 + zCoord_ix0_izm1 + zCoord_ixm1_iz0 + zCoord_ix0_iz0);
					dx = x - xq00; dz = z - zq00;
					double r00 = sqrt(dx*dx + dz*dz);
					double rMin = r00;
					char qID = 0;
						
					double xq10 = 0.25*(xCoord_ix0_izm1 + xCoord_ixp1_izm1 + xCoord_ix0_iz0 + xCoord_ixp1_iz0);
					double zq10 = 0.25*(zCoord_ix0_izm1 + zCoord_ixp1_izm1 + zCoord_ix0_iz0 + zCoord_ixp1_iz0);
					dx = x - xq10; dz = z - zq10;
					double r10 = sqrt(dx*dx + dz*dz);
					if(rMin > r10) { rMin = r10; qID = 1;}
						
					double xq01 = 0.25*(xCoord_ixm1_iz0 + xCoord_ix0_iz0 + xCoord_ixm1_izp1 + xCoord_ix0_izp1);
					double zq01 = 0.25*(zCoord_ixm1_iz0 + zCoord_ix0_iz0 + zCoord_ixm1_izp1 + zCoord_ix0_izp1);
					dx = x - xq01; dz = z - zq01;
					double r01 = sqrt(dx*dx + dz*dz);
					if(rMin > r01) { rMin = r01; qID = 2;}

					double xq11 = 0.25*(xCoord_ix0_iz0 + xCoord_ixp1_iz0 + xCoord_ix0_izp1 + xCoord_ixp1_izp1);
					double zq11 = 0.25*(zCoord_ix0_iz0 + zCoord_ixp1_iz0 + zCoord_ix0_izp1 + zCoord_ixp1_izp1);
					dx = x - xq11; dz = z - zq11;
					double r11 = sqrt(dx*dx + dz*dz);
					if(rMin > r11) { rMin = r11; qID = 3;}

					//Use Bi-Linear interpolation on Irregular mesh
					double relX=0, relZ=0;
					if(qID == 0)
					{
						arRelCoordXZ[0] = xCoord_ix0_izm1 - xCoord_ixm1_izm1;
						arRelCoordXZ[1] = zCoord_ix0_izm1 - zCoord_ixm1_izm1;
						arRelCoordXZ[2] = xCoord_ixm1_iz0 - xCoord_ixm1_izm1;
						arRelCoordXZ[3] = zCoord_ixm1_iz0 - zCoord_ixm1_izm1;
						arRelCoordXZ[4] = xCoord_ix0_iz0 - xCoord_ixm1_izm1;
						arRelCoordXZ[5] = zCoord_ix0_iz0 - zCoord_ixm1_izm1;
						relX = x - xCoord_ixm1_izm1;
						relZ = z - zCoord_ixm1_izm1;
						if(arEX != 0)
						{
							arReEx[0] = arEX[ofst_ixm1_izm1];
							arImEx[0] = arEX[ofst_ixm1_izm1 + 1];
							arReEx[1] = arEX[ofst_ix0_izm1];
							arImEx[1] = arEX[ofst_ix0_izm1 + 1];
							arReEx[2] = arEX[ofst_ixm1_iz0];
							arImEx[2] = arEX[ofst_ixm1_iz0 + 1];
							arReEx[3] = arEX[ofst_ix0_iz0];
							arImEx[3] = arEX[ofst_ix0_iz0 + 1];
						}
						if(arEZ != 0)
						{
							arReEz[0] = arEZ[ofst_ixm1_izm1];
							arImEz[0] = arEZ[ofst_ixm1_izm1 + 1];
							arReEz[1] = arEZ[ofst_ix0_izm1];
							arImEz[1] = arEZ[ofst_ix0_izm1 + 1];
							arReEz[2] = arEZ[ofst_ixm1_iz0];
							arImEz[2] = arEZ[ofst_ixm1_iz0 + 1];
							arReEz[3] = arEZ[ofst_ix0_iz0];
							arImEz[3] = arEZ[ofst_ix0_iz0 + 1];
						}
					}
					else if(qID == 1)
					{
						arRelCoordXZ[0] = xCoord_ixp1_izm1 - xCoord_ix0_izm1;
						arRelCoordXZ[1] = zCoord_ixp1_izm1 - zCoord_ix0_izm1;
						arRelCoordXZ[2] = xCoord_ix0_iz0 - xCoord_ix0_izm1;
						arRelCoordXZ[3] = zCoord_ix0_iz0 - zCoord_ix0_izm1;
						arRelCoordXZ[4] = xCoord_ixp1_iz0 - xCoord_ix0_izm1;
						arRelCoordXZ[5] = zCoord_ixp1_iz0 - zCoord_ix0_izm1;
						relX = x - xCoord_ix0_izm1;
						relZ = z - zCoord_ix0_izm1;
						if(arEX != 0)
						{
							arReEx[0] = arEX[ofst_ix0_izm1];
							arImEx[0] = arEX[ofst_ix0_izm1 + 1];
							arReEx[1] = arEX[ofst_ixp1_izm1];
							arImEx[1] = arEX[ofst_ixp1_izm1 + 1];
							arReEx[2] = arEX[ofst_ix0_iz0];
							arImEx[2] = arEX[ofst_ix0_iz0 + 1];
							arReEx[3] = arEX[ofst_ixp1_iz0];
							arImEx[3] = arEX[ofst_ixp1_iz0 + 1];
						}
						if(arEZ != 0)
						{
							arReEz[0] = arEZ[ofst_ix0_izm1];
							arImEz[0] = arEZ[ofst_ix0_izm1 + 1];
							arReEz[1] = arEZ[ofst_ixp1_izm1];
							arImEz[1] = arEZ[ofst_ixp1_izm1 + 1];
							arReEz[2] = arEZ[ofst_ix0_iz0];
							arImEz[2] = arEZ[ofst_ix0_iz0 + 1];
							arReEz[3] = arEZ[ofst_ixp1_iz0];
							arImEz[3] = arEZ[ofst_ixp1_iz0 + 1];
						}
					}
					else if(qID == 2)
					{
						arRelCoordXZ[0] = xCoord_ix0_iz0 - xCoord_ixm1_iz0;
						arRelCoordXZ[1] = zCoord_ix0_iz0 - zCoord_ixm1_iz0;
						arRelCoordXZ[2] = xCoord_ixm1_izp1 - xCoord_ixm1_iz0;
						arRelCoordXZ[3] = zCoord_ixm1_izp1 - zCoord_ixm1_iz0;
						arRelCoordXZ[4] = xCoord_ix0_izp1 - xCoord_ixm1_iz0;
						arRelCoordXZ[5] = zCoord_ix0_izp1 - zCoord_ixm1_iz0;
						relX = x - xCoord_ixm1_iz0;
						relZ = z - zCoord_ixm1_iz0;
						if(arEX != 0)
						{
							arReEx[0] = arEX[ofst_ixm1_iz0];
							arImEx[0] = arEX[ofst_ixm1_iz0 + 1];
							arReEx[1] = arEX[ofst_ix0_iz0];
							arImEx[1] = arEX[ofst_ix0_iz0 + 1];
							arReEx[2] = arEX[ofst_ixm1_izp1];
							arImEx[2] = arEX[ofst_ixm1_izp1 + 1];
							arReEx[3] = arEX[ofst_ix0_izp1];
							arImEx[3] = arEX[ofst_ix0_izp1 + 1];
						}
						if(arEZ != 0)
						{
							arReEz[0] = arEZ[ofst_ixm1_iz0];
							arImEz[0] = arEZ[ofst_ixm1_iz0 + 1];
							arReEz[1] = arEZ[ofst_ix0_iz0];
							arImEz[1] = arEZ[ofst_ix0_iz0 + 1];
							arReEz[2] = arEZ[ofst_ixm1_izp1];
							arImEz[2] = arEZ[ofst_ixm1_izp1 + 1];
							arReEz[3] = arEZ[ofst_ix0_izp1];
							arImEz[3] = arEZ[ofst_ix0_izp1 + 1];
						}
					}
					else if(qID == 3)
					{
						arRelCoordXZ[0] = xCoord_ixp1_iz0 - xCoord_ix0_iz0;
						arRelCoordXZ[1] = zCoord_ixp1_iz0 - zCoord_ix0_iz0;
						arRelCoordXZ[2] = xCoord_ix0_izp1 - xCoord_ix0_iz0;
						arRelCoordXZ[3] = zCoord_ix0_izp1 - zCoord_ix0_iz0;
						arRelCoordXZ[4] = xCoord_ixp1_izp1 - xCoord_ix0_iz0;
						arRelCoordXZ[5] = zCoord_ixp1_izp1 - zCoord_ix0_iz0;
						relX = x - xCoord_ix0_iz0;
						relZ = z - zCoord_ix0_iz0;
						if(arEX != 0)
						{
							arReEx[0] = arEX[ofst_ix0_iz0];
							arImEx[0] = arEX[ofst_ix0_iz0 + 1];
							arReEx[1] = arEX[ofst_ixp1_iz0];
							arImEx[1] = arEX[ofst_ixp1_iz0 + 1];
							arReEx[2] = arEX[ofst_ix0_izp1];
							arImEx[2] = arEX[ofst_ix0_izp1 + 1];
							arReEx[3] = arEX[ofst_ixp1_izp1];
							arImEx[3] = arEX[ofst_ixp1_izp1 + 1];
						}
						if(arEZ != 0)
						{
							arReEz[0] = arEZ[ofst_ix0_iz0];
							arImEz[0] = arEZ[ofst_ix0_iz0 + 1];
							arReEz[1] = arEZ[ofst_ixp1_iz0];
							arImEz[1] = arEZ[ofst_ixp1_iz0 + 1];
							arReEz[2] = arEZ[ofst_ix0_izp1];
							arImEz[2] = arEZ[ofst_ix0_izp1 + 1];
							arReEz[3] = arEZ[ofst_ixp1_izp1];
							arImEz[3] = arEZ[ofst_ixp1_izp1 + 1];
						}
					}

					double resReEx=0, resImEx=0, resReEz=0, resImEz=0;

					//OC04072017
					bool interpCanBeDone = true;
					double absCoordThres = 1.e+20;
					if((relX < -absCoordThres) || (relX > absCoordThres) || (relZ < -absCoordThres) || (relZ > absCoordThres)) interpCanBeDone = false;
					if(interpCanBeDone)
					{
						for(int jj = 0; jj < 6; jj++)
						{
							double curCoord = arRelCoordXZ[jj];
							if((curCoord < -absCoordThres) || (curCoord > absCoordThres))
							{
								interpCanBeDone = false; break;
							}
						}
					}

					//if(arEX != 0)
					if(interpCanBeDone && (arEX != 0)) //OC04072017
					{
						//static double Interp2dBiLinVar(double x, double y, double* arXY, double* arF)
						//bilinear interpolation on irregular mesh, for relative arguments, first point is x = 0, y = 0
						//arXY is flat array of coordinates of 3 other points {x10, y10, x01, y01, x11, y11}
						//resReEx = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arReEx);
						///resImEx = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arImEx);
						resReEx = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arReEx);
						resImEx = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arImEx);

						//Correction of Electric FIeld values, using Intensity ones
						double appIx = resReEx*resReEx + resImEx*resImEx;
						if(appIx > 0.)
						{
							for(int i=0; i<4; i++) { ReE = arReEx[i]; ImE = arImEx[i]; arIx[i] = ReE*ReE + ImE*ImE;}
							//double resIx = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arIx);
							double resIx = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arIx);
							if(resIx <= 0)
							{
								resReEx = 0.; resImEx = 0.;
							}
							else 
							{
								double fact = sqrt(resIx/appIx);
								resReEx *= fact; resImEx *= fact;
							}
						}
					}
					//if(arEZ != 0)
					if(interpCanBeDone && (arEZ != 0)) //OC04072017
					{
						//resReEz = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arReEz);
						//resImEz = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arImEz);
						resReEz = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arReEz);
						resImEz = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arImEz);

						//Correction of Electric FIeld values, using Intensity ones
						double appIz = resReEz*resReEz + resImEz*resImEz;
						if(appIz > 0.)
						{
							for(int i=0; i<4; i++) { ReE = arReEz[i]; ImE = arImEz[i]; arIz[i] = ReE*ReE + ImE*ImE;}
							//double resIz = CGenMathInterp::Interp2dBiLinVar(x, z, arRelCoordXZ, arIz);
							double resIz = CGenMathInterp::Interp2dBiLinVar(relX, relZ, arRelCoordXZ, arIz);
							if(resIz <= 0)
							{
								resReEz = 0.; resImEz = 0.;
							}
							else
							{
								double fact = sqrt(resIz/appIz);
								resReEz *= fact; resImEz *= fact;
							}
						}
					}

					*t_ExRes = (float)resReEx; *(t_ExRes+1) = (float)resImEx;
					*t_EzRes = (float)resReEz; *(t_EzRes+1) = (float)resImEz;
				}
				t_ExRes += 2;
				t_EzRes += 2;
			}
			x += pWfr->xStep;
		}
		z += pWfr->zStep;
	}
	
	//OCTEST (commented-out)
	//if(fabs(pWfr->RobsZ + 18.079) > 0.1)
	//{//OCTEST

	if(waveFrontTermWasTreated) TreatStronglyOscillatingTerm(*pWfr, 'a');

	//}//OCTEST
	return 0;
}

//*************************************************************************

int srTMirror::WfrInterpolOnOrigGrid(srTSRWRadStructAccessData* pWfr, double* arRayTrCoord, float* arEX, float* arEZ, double xMin, double xMax, double zMin, double zMax)
{
	if((pWfr == 0) || (arRayTrCoord == 0) || ((arEX == 0) && (arEZ == 0))) return FAILED_INTERPOL_ELEC_FLD;
	//if((pWfr == 0) || (arRayTrCoord == 0) || (arOptPathDif == 0) || ((arEX == 0) && (arEZ == 0))) return FAILED_INTERPOL_ELEC_FLD;

	bool isCoordRepres = (pWfr->Pres == 0);
	bool isFreqRepres = (pWfr->PresT == 0);
	bool waveFrontTermWasTreated = false;

	//if(isFreqRepres && isCoordRepres && WaveFrontTermCanBeTreated(*pWfr))
	if(isFreqRepres && isCoordRepres && WaveFrontTermCanBeTreated(*pWfr, false)) //OC300914
	{
		float *pPrevBaseRadX = pWfr->pBaseRadX; 
		float *pPrevBaseRadZ = pWfr->pBaseRadZ;
		pWfr->pBaseRadX = arEX;
		pWfr->pBaseRadZ = arEZ;
		
		//testoc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//TreatStronglyOscillatingTerm(*pWfr, 'r');
		//end testoc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		//if(!((fabs(pWfr->RobsX + 0.99) < 0.1) && (fabs(pWfr->RobsZ + 0.99) < 0.1)))
		//{//OCTEST
		TreatStronglyOscillatingTermIrregMesh(*pWfr, arRayTrCoord, xMin, xMax, zMin, zMax, 'r');
		//}//OCTEST
		
		pWfr->pBaseRadX = pPrevBaseRadX;
		pWfr->pBaseRadZ = pPrevBaseRadZ;
		waveFrontTermWasTreated = true;
	}

	float *t_ExRes = pWfr->pBaseRadX;
	float *t_EzRes = pWfr->pBaseRadZ;

	//long PerX = (pWfr->ne) << 1;
	//long PerZ = PerX*(pWfr->nx);
	//long nTot = PerZ*(pWfr->nz);
	//long nTot_mi_1 = nTot - 1;
	//long nx_mi_1 = pWfr->nx - 1;
	//long nz_mi_1 = pWfr->nz - 1;

	long long PerX = (pWfr->ne) << 1;
	long long PerZ = PerX*(pWfr->nx);
	long long nTot = PerZ*(pWfr->nz);
	long long nTot_mi_1 = nTot - 1;
	long nx_mi_1 = pWfr->nx - 1;
	long nz_mi_1 = pWfr->nz - 1;


	double f0m1, fm10, f00, f10, f01, f11, a10, a01, a11, a20, a02;

	//OCTEST
	//if((fabs(pWfr->RobsX + 20.596) < 0.1) && (fabs(pWfr->RobsZ + 18.079) < 0.1))
	//if(fabs(pWfr->RobsZ + 18.079) < 0.1)
	//{
	//	float *t_arEX = arEX, *t_arEZ = arEZ;
	//	for(int k=0; k<nTot; k++)
	//	{
	//		*(t_ExRes++) = *(t_arEX++);
	//		*(t_EzRes++) = *(t_arEZ++);
	//	}
	//	return 0;
	//	//int aha = 1;
	//}
	//END OCTEST

	long ix0=-1, iz0=-1;
	double phEn, x, z = pWfr->zStart;

	const double dMax = 1.E+20;
	double dx, dz, dTest, dTest0, dTestPrev;

	for(long iz=0; iz<pWfr->nz; iz++)
	{
		x = pWfr->xStart;
		for(long ix=0; ix<pWfr->nx; ix++)
		{
			bool pointIsInsideNonZeroSquare = (x >= xMin) && (x <= xMax) && (z >= zMin) && (z <= zMax);

			phEn = pWfr->eStart;
			for(long ie=0; ie<pWfr->ne; ie++)
			{
				long two_ie = ie << 1;

				if(pointIsInsideNonZeroSquare)
				{//perform interpolation on irregular mesh (bilinear or based on 12 points)
					//find indexes of the relevant point for the interpolation

					//OCTEST
					//if((fabs(pWfr->RobsX + 0.99) < 0.1) && (fabs(pWfr->RobsZ + 0.99) < 0.1))
					//{
					//	//if((fabs(rx_m10) > dMax) || (fabs(rx_10) > dMax) || (fabs(rz_0m1) > dMax) || (fabs(rz_01) > dMax))
					//	//if((iz == 1848) && (ix == 442))// && (fabs(x + 0.0047703) < 0.00001))
					//	if((iz == 748) && (ix == 350))// && (fabs(x + 0.0047703) < 0.00001))
					//	{ 
					//		int aha = 1;
					//	}
					//}

					if(ix0 < 0) ix0 = ix;
					if(iz0 < 0) iz0 = iz;

					bool pointFound = false, candPointFound = false;
					//bool isLeftBordX = false, isRightBordX = false;
					//bool isLeftBordZ = false, isRightBordZ = false;

					long ix0pr = -1, iz0pr = -1;
					while((ix0 != ix0pr) && (iz0 != iz0pr)) 
					{//This while loop is required for a "tilted/rotated" mesh  (to check how ir works!)
						ix0pr = ix0; iz0pr = iz0;

						//long iz0_PerZ_p_2_ie = iz0*PerZ + two_ie;
						long long iz0_PerZ_p_2_ie = iz0*PerZ + two_ie;

						dTestPrev = 1.E+23;
						//long ofst = ix0*PerX + iz0_PerZ_p_2_ie;
						long long ofst = ix0*PerX + iz0_PerZ_p_2_ie;
						if (ofst < nTot_mi_1)
						{
							dx = x - arRayTrCoord[ofst];
							dz = z - arRayTrCoord[ofst + 1];
							dTestPrev = sqrt(dx*dx + dz*dz);
						}
						dTest0 = dTestPrev;

						pointFound = false;
						candPointFound = false;
						//isLeftBordX = false; isRightBordX = false;

						long ix0orig = ix0;
						//for(int iix = ix0 - 1; iix >= 0; iix--)
						for(int iix = ix0orig - 1; iix >= 0; iix--) //OC200414
						{
							ofst = iix*PerX + iz0_PerZ_p_2_ie;
							if(ofst < nTot_mi_1)
							{
								dx = x - arRayTrCoord[ofst];
								dz = z - arRayTrCoord[ofst + 1];
								dTest = sqrt(dx*dx + dz*dz);
								//if((dTest > dMax) && (!candPointFound)) continue;
								if(dTest > dMax)
								{
									if(dTestPrev < dMax) break;
									if(!candPointFound) continue;
								}
								if(dTest < dTestPrev)
								{
									ix0 = iix; dTestPrev = dTest; candPointFound = true;
								}
								else
								{
									if(candPointFound)
									{
										pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isLeftBordX = true;
									}
									break;
								}
							}
						}
						if(!pointFound)
						{
							//dTestPrev = 1.E+23;
							dTestPrev = dTest0;
							candPointFound = false;
							//for(int iix = ix0 + 1; iix < pWfr->nx; iix++)
							for(int iix = ix0orig + 1; iix < pWfr->nx; iix++) //OC200414
							{
								ofst = iix*PerX + iz0_PerZ_p_2_ie;
								if(ofst < nTot_mi_1)
								{
									dx = x - arRayTrCoord[ofst];
									dz = z - arRayTrCoord[ofst + 1];
									dTest = sqrt(dx*dx + dz*dz);
									//if((dTest > dMax) && (!candPointFound)) continue;
									if(dTest > dMax)
									{
										if(dTestPrev < dMax) break;
										if(!candPointFound) continue;
									}
									if(dTest < dTestPrev)
									{
										ix0 = iix; dTestPrev = dTest; candPointFound = true;
									}
									else
									{
										//if(candPointFound) pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isRightBordX = true;
										break;
									}
								}
							}
						}

						//long ix0_PerX_p_2_ie = ix0*PerX + two_ie;
						long long ix0_PerX_p_2_ie = ix0*PerX + two_ie;

						dTestPrev = 1.E+23;
						ofst = iz0*PerZ + ix0_PerX_p_2_ie;
						if (ofst < nTot_mi_1)
						{
							dx = x - arRayTrCoord[ofst];
							dz = z - arRayTrCoord[ofst + 1];
							dTestPrev = sqrt(dx*dx + dz*dz);
						}
						dTest0 = dTestPrev;

						pointFound = false;
						candPointFound = false;
						//isLeftBordZ = false; isRightBordZ = false;

						long iz0orig = iz0;
						//for(int iiz = iz0 - 1; iiz >= 0; iiz--)
						for(int iiz = iz0orig - 1; iiz >= 0; iiz--) //OC200414
						{
							ofst = iiz*PerZ + ix0_PerX_p_2_ie;
							if(ofst < nTot_mi_1)
							{
								dx = x - arRayTrCoord[ofst];
								dz = z - arRayTrCoord[ofst + 1];
								dTest = sqrt(dx*dx + dz*dz);
								//if((dTest > dMax) && (!candPointFound)) continue;
								if(dTest > dMax)
								{
									if(dTestPrev < dMax) break;
									if(!candPointFound) continue;
								}
								if(dTest < dTestPrev)
								{
									iz0 = iiz; dTestPrev = dTest; candPointFound = true;
								}
								else
								{
									if(candPointFound)
									{
										pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) isLeftBordZ = true;
									}
									//OCTEST
									//if((dTest > dMax) || (dTestPrev > dMax))
									//{
									//	int aha = 1;
									//}
									break;
								}
							}
						}
						if(!pointFound)
						{
							//dTestPrev = 1.E+23;
							dTestPrev = dTest0;
							candPointFound = false;
							//for(int iiz = iz0 + 1; iiz < pWfr->nz; iiz++)
							for(int iiz = iz0orig + 1; iiz < pWfr->nz; iiz++) //OC200414
							{
								ofst = iiz*PerZ + ix0_PerX_p_2_ie;
								if(ofst < nTot_mi_1)
								{
									dx = x - arRayTrCoord[ofst];
									dz = z - arRayTrCoord[ofst + 1];
									dTest = sqrt(dx*dx + dz*dz);
									//if((dTest > dMax) && (!candPointFound)) continue;
									if(dTest > dMax)
									{
										if(dTestPrev < dMax) 
										{
											//if(dx < 0) isRightBordZ = true;
											//else if(dx > 0) isLeftBordZ = true;
											break;
										}
										if(!candPointFound) continue;
									}
									if(dTest < dTestPrev)
									{
										iz0 = iiz; dTestPrev = dTest; candPointFound = true;
									}
									else
									{
										//if(candPointFound) pointFound = true;
										//if((dTest > dMax) && (dTestPrev < dMax)) 
										//{
										//	if(dx < 0) isRightBordZ = true;
										//	else if(dx > 0) isLeftBordZ = true;
										//}

										//OCTEST
										//if((dTest > dMax) || (dTestPrev > dMax))
										//{
										//	int aha = 1;
										//}
										break;
									}
								}
							}
						}
					}
					//calculate indexes of other points and interpolate 
					//2 cases are considered: "bi-linear" 2D interpolation based on 4 points and "bi-quadratic" 2D interpolation based on 5 points (mesh can be irregular)
					const double relTolEqualStep = 1.e-04; //to tune
					
					if(ix0 < 0) ix0 = 0;
					else if(ix0 >= nx_mi_1) ix0 = nx_mi_1 - 1;
					if(iz0 < 0) iz0 = 0;
					else if(iz0 >= nz_mi_1) iz0 = nz_mi_1 - 1;

					//long ofst_00 = ix0*PerX + iz0*PerZ + two_ie;
					//long ofst_m10 = ofst_00 - PerX;
					//long ofst_10 = ofst_00 + PerX;
					//long ofst_0m1 = ofst_00 - PerZ;
					//long ofst_01 = ofst_00 + PerZ;
					long long ofst_00 = ix0*PerX + iz0*PerZ + two_ie;
					long long ofst_m10 = ofst_00 - PerX;
					long long ofst_10 = ofst_00 + PerX;
					long long ofst_0m1 = ofst_00 - PerZ;
					long long ofst_01 = ofst_00 + PerZ;

					if(ix0 == 0) ofst_m10 = ofst_00;
					if(ix0 == nx_mi_1) ofst_10 = ofst_00;
					if(iz0 == 0) ofst_0m1 = ofst_00;
					if(iz0 == nz_mi_1) ofst_10 = ofst_00;

					//if((ix0 == 0) || isLeftBordX) ofst_m10 = ofst_00;
					//if((ix0 == nx_mi_1) || isRightBordX) ofst_10 = ofst_00;
					//if((iz0 == 0) || isLeftBordZ) ofst_0m1 = ofst_00;
					//if((iz0 == nz_mi_1) || isRightBordZ) ofst_10 = ofst_00;

					//long ofst_00_p_1 = ofst_00 + 1;
					//long ofst_m10_p_1 = ofst_m10 + 1;
					//long ofst_10_p_1 = ofst_10 + 1;
					//long ofst_0m1_p_1 = ofst_0m1 + 1;
					//long ofst_01_p_1 = ofst_01 + 1;
					long long ofst_00_p_1 = ofst_00 + 1;
					long long ofst_m10_p_1 = ofst_m10 + 1;
					long long ofst_10_p_1 = ofst_10 + 1;
					long long ofst_0m1_p_1 = ofst_0m1 + 1;
					long long ofst_01_p_1 = ofst_01 + 1;

					double x_00 = arRayTrCoord[ofst_00], z_00 = arRayTrCoord[ofst_00_p_1];
					double x_m10 = arRayTrCoord[ofst_m10], z_m10 = arRayTrCoord[ofst_m10_p_1];
					double x_10 = arRayTrCoord[ofst_10], z_10 = arRayTrCoord[ofst_10_p_1];
					double x_0m1 = arRayTrCoord[ofst_0m1], z_0m1 = arRayTrCoord[ofst_0m1_p_1];
					double x_01 = arRayTrCoord[ofst_01], z_01 = arRayTrCoord[ofst_01_p_1];

					double rx_m10 = x_m10 - x_00, rz_m10 = z_m10 - z_00;
					double rx_10 = x_10 - x_00, rz_10 = z_10 - z_00;
					double rx_0m1 = x_0m1 - x_00, rz_0m1 = z_0m1 - z_00;
					double rx_01 = x_01 - x_00, rz_01 = z_01 - z_00;
					double dx_00 = x - x_00, dz_00 = z - z_00;

					//OCTEST
					//if(fabs(pWfr->RobsZ + 0.37) < 0.1)
					//{
					//	//if((fabs(rx_m10) > dMax) || (fabs(rx_10) > dMax) || (fabs(rz_0m1) > dMax) || (fabs(rz_01) > dMax))
					//	if((fabs(z + 0.0054932) < 0.00001)) //&& (fabs(x + 0.0047703) < 0.00001))
					//	{
					//		int aha = 1;
					//	}
					//}

					bool rx_m10_isNotOK = ((fabs(rx_m10) > dMax) || (rx_m10 == 0));
					bool rx_10_isNotOK = ((fabs(rx_10) > dMax) || (rx_10 == 0));
					if(rx_m10_isNotOK && rx_10_isNotOK) goto SetFieldToZero;

					bool rz_0m1_isNotOK = ((fabs(rz_0m1) > dMax) || (rz_0m1 == 0));
					bool rz_01_isNotOK = ((fabs(rz_01) > dMax) || (rz_01 == 0));
					if(rz_0m1_isNotOK && rz_01_isNotOK) goto SetFieldToZero;

					if(rx_m10_isNotOK) 
					{
						rx_m10 = -rx_10;
						ofst_m10 = ofst_00;
					}
					else if(rx_10_isNotOK) 
					{
						rx_10 = -rx_m10;
						ofst_10 = ofst_00;
					}

					bool rx_0m1_isNotOK = (fabs(rx_0m1) > dMax);
					bool rx_01_isNotOK = (fabs(rx_01) > dMax);
					if(rx_0m1_isNotOK) 
					{
						rx_0m1 = 0.; //??
						ofst_0m1 = ofst_00;
					}
					if(rx_01_isNotOK) 
					{
						rx_01 = 0.;
						ofst_01 = ofst_00;
					}

					if(rz_0m1_isNotOK) 
					{
						rz_0m1 = -rz_01;
						ofst_0m1 = ofst_00;
					}
					else if(rz_01_isNotOK) 
					{
						rz_01 = -rz_0m1;
						ofst_01 = ofst_00;
					}

					bool rz_m10_isNotOK = (fabs(rz_m10) > dMax);
					bool rz_10_isNotOK = (fabs(rz_10) > dMax);
					if(rz_m10_isNotOK) 
					{
						rz_m10 = 0.;
						ofst_m10 = ofst_00;
					}
					if(rz_10_isNotOK) 
					{
						rz_10 = 0.;
						ofst_10 = ofst_00;
					}

					const double maxRelArg = 1.5;
					//OC200414
					//if(rx_m10_isNotOK || rx_10_isNotOK || rx_0m1_isNotOK || rx_01_isNotOK)
					//{
					//	double twoRx = ::fabs(rx_m10) + ::fabs(rx_10);
					//	if(::fabs(dx_00/twoRx) > maxRelArg) goto SetFieldToZero;
					//}
					//if(rz_0m1_isNotOK || rz_01_isNotOK || rz_m10_isNotOK || rz_10_isNotOK)
					//{
					//	double twoRz = ::fabs(rz_0m1) + ::fabs(rz_01);
					//	if(::fabs(dz_00/twoRz) > maxRelArg) goto SetFieldToZero;
					//}

					double d_rx_m10 = dx_00 - rx_m10;
					double d_rx_10 = dx_00 - rx_10;
					if(d_rx_m10*d_rx_10 > 0.)
					{
						if((::fabs(rx_m10)*maxRelArg < ::fabs(d_rx_m10)) || (::fabs(rx_10)*maxRelArg < ::fabs(d_rx_10))) goto SetFieldToZero;
					}
					double d_rz_0m1 = dz_00 - rz_0m1;
					double d_rz_01 = dz_00 - rz_01;
					if(d_rz_0m1*d_rz_01 > 0.)
					{
						if((::fabs(rz_0m1)*maxRelArg < ::fabs(d_rz_0m1)) || (::fabs(rz_01)*maxRelArg < ::fabs(d_rz_01))) goto SetFieldToZero;
					}

					if(m_wfrInterpMode == 1)
					{//bi-linear, based on 4 points
						double sp_m10 = rx_m10*dx_00 + rz_m10*dz_00;
						double sp_10 = rx_10*dx_00 + rz_10*dz_00;
						double sp_0m1 = rx_0m1*dx_00 + rz_0m1*dz_00;
						double sp_01 = rx_01*dx_00 + rz_01*dz_00;

						bool initPointMoved = false;
						if((sp_m10 > 0) && (sp_10 < 0))
						{
							if(ix0 > 0) { ix0--; initPointMoved = true;}
						}
						if((sp_0m1 > 0) && (sp_01 < 0))
						{
							if(iz0 > 0) { iz0--; initPointMoved = true;}
						}

						//long ofst0 = initPointMoved? (ix0*PerX + iz0*PerZ + two_ie) : ofst_00;
						long long ofst0 = initPointMoved? (ix0*PerX + iz0*PerZ + two_ie) : ofst_00;

						//long ofst0_p_PerX = ofst0 + PerX;
						//long ofst0_p_PerZ = ofst0 + PerZ;
						//long ofst0_p_PerX_p_PerZ = ofst0_p_PerZ + PerX;
						long long ofst0_p_PerX = ofst0 + PerX;
						long long ofst0_p_PerZ = ofst0 + PerZ;
						long long ofst0_p_PerX_p_PerZ = ofst0_p_PerZ + PerX;
						double x00 = arRayTrCoord[ofst0], x10 = arRayTrCoord[ofst0_p_PerX], x01 = arRayTrCoord[ofst0_p_PerZ], x11 = arRayTrCoord[ofst0_p_PerX_p_PerZ];
						//long ofst0_p_1 = ofst0 + 1;
						//long ofst0_p_PerX_p_1 = ofst0_p_PerX + 1;
						//long ofst0_p_PerZ_p_1 = ofst0_p_PerZ + 1;
						//long ofst0_p_PerX_p_PerZ_p_1 = ofst0_p_PerX_p_PerZ + 1;
						long long ofst0_p_1 = ofst0 + 1;
						long long ofst0_p_PerX_p_1 = ofst0_p_PerX + 1;
						long long ofst0_p_PerZ_p_1 = ofst0_p_PerZ + 1;
						long long ofst0_p_PerX_p_PerZ_p_1 = ofst0_p_PerX_p_PerZ + 1;
						double z00 = arRayTrCoord[ofst0_p_1], z10 = arRayTrCoord[ofst0_p_PerX_p_1], z01 = arRayTrCoord[ofst0_p_PerZ_p_1], z11 = arRayTrCoord[ofst0_p_PerX_p_PerZ_p_1];

						double rX = x - x00, rZ = z - z00;
						double rX10 = x10 - x00, rZ10 = z10 - z00;
						double rX01 = x01 - x00, rZ01 = z01 - z00;
						double rX11 = x11 - x00, rZ11 = z11 - z00;

						double absTolX = relTolEqualStep*fabs(rX10);
						double absTolZ = relTolEqualStep*fabs(rZ01);
						bool isRecX = (fabs(rX01) < absTolX) && (fabs(x11 - x10) < absTolX);
						bool isRecZ = (fabs(rZ10) < absTolZ) && (fabs(z11 - z01) < absTolZ);
						if(isRecX && isRecZ)
						{//regular rectangular mesh
							double xt = rX/rX10, zt = rZ/rZ01;
							if(arEX != 0)
							{
								f00 = arEX[ofst0]; f10 = arEX[ofst0_p_PerX]; f01 = arEX[ofst0_p_PerZ]; f11 = arEX[ofst0_p_PerX_p_PerZ];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*t_ExRes = (float)(xt*(a10 + a11*zt) + a01*zt + f00);

								f00 = arEX[ofst0_p_1]; f10 = arEX[ofst0_p_PerX_p_1]; f01 = arEX[ofst0_p_PerZ_p_1]; f11 = arEX[ofst0_p_PerX_p_PerZ_p_1];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*(t_ExRes + 1) = (float)(xt*(a10 + a11*zt) + a01*zt + f00);
							}
							if(arEZ != 0)
							{
								f00 = arEZ[ofst0]; f10 = arEZ[ofst0_p_PerX]; f01 = arEZ[ofst0_p_PerZ]; f11 = arEZ[ofst0_p_PerX_p_PerZ];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*t_EzRes = (float)(xt*(a10 + a11*zt) + a01*zt + f00);

								f00 = arEZ[ofst0_p_1]; f10 = arEZ[ofst0_p_PerX_p_1]; f01 = arEZ[ofst0_p_PerZ_p_1]; f11 = arEZ[ofst0_p_PerX_p_PerZ_p_1];
								a10 = f10 - f00; a01 = f01 - f00; a11 = f00 - f01 - f10 + f11;
								*(t_EzRes + 1) = (float)(xt*(a10 + a11*zt) + a01*zt + f00);
							}
						}
						else
						{//irregular mesh (general case)
							double arXZ[] = {rX10, rZ10, rX01, rZ01, rX11, rZ11};
							if(arEX != 0)
							{
								double arER[] = {arEX[ofst0], arEX[ofst0_p_PerX], arEX[ofst0_p_PerZ], arEX[ofst0_p_PerX_p_PerZ]};
								*t_ExRes = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arER);
								double arEI[] = {arEX[ofst0_p_1], arEX[ofst0_p_PerX_p_1], arEX[ofst0_p_PerZ_p_1], arEX[ofst0_p_PerX_p_PerZ_p_1]};
								*(t_ExRes + 1) = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arEI);
							}
							if(arEZ != 0)
							{
								double arER[] = {arEZ[ofst0], arEZ[ofst0_p_PerX], arEZ[ofst0_p_PerZ], arEZ[ofst0_p_PerX_p_PerZ]};
								*t_EzRes = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arER);
								double arEI[] = {arEZ[ofst0_p_1], arEZ[ofst0_p_PerX_p_1], arEZ[ofst0_p_PerZ_p_1], arEZ[ofst0_p_PerX_p_PerZ_p_1]};
								*(t_EzRes + 1) = (float)CGenMathInterp::Interp2dBiLinVar(rX, rZ, arXZ, arEI);
							}
						}
					}
					else if(m_wfrInterpMode == 2)
					{//bi-quadratic, based on 5 points
						double absTolX = relTolEqualStep*fabs(rx_10);
						double absTolZ = relTolEqualStep*fabs(rz_01);
						bool isRecX = (fabs(rx_0m1) < absTolX) && (fabs(rx_01) < absTolX);
						bool isRecZ = (fabs(rz_m10) < absTolZ) && (fabs(rz_10) < absTolX);
						if(isRecX && isRecZ)
						{
							bool isEquidistX = (fabs(rx_m10 + rx_10) < absTolX);
							bool isEquidistZ = (fabs(rz_0m1 + rz_01) < absTolZ);
							if(isEquidistX && isEquidistZ)
							{//regular rectangular mesh
								double xt = dx_00/rx_10, zt = dz_00/rz_01;
								if(arEX != 0)
								{
									f0m1 = arEX[ofst_0m1]; fm10 = arEX[ofst_m10]; f00 = arEX[ofst_00]; f10 = arEX[ofst_10]; f01 = arEX[ofst_01];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*t_ExRes = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);

									f0m1 = arEX[ofst_0m1_p_1]; fm10 = arEX[ofst_m10_p_1]; f00 = arEX[ofst_00_p_1]; f10 = arEX[ofst_10_p_1]; f01 = arEX[ofst_01_p_1];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*(t_ExRes + 1) = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);
								}
								if(arEZ != 0)
								{
									f0m1 = arEZ[ofst_0m1]; fm10 = arEZ[ofst_m10]; f00 = arEZ[ofst_00]; f10 = arEZ[ofst_10]; f01 = arEZ[ofst_01];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*t_EzRes = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);

									f0m1 = arEZ[ofst_0m1_p_1]; fm10 = arEZ[ofst_m10_p_1]; f00 = arEZ[ofst_00_p_1]; f10 = arEZ[ofst_10_p_1]; f01 = arEZ[ofst_01_p_1];
									a10 = 0.5*(f10 - fm10); a01 = 0.5*(f01 - f0m1); a20 = 0.5*(f10 + fm10) - f00; a02 = 0.5*(f01 + f0m1) - f00;
									*(t_EzRes + 1) = (float)(xt*(xt*a20 + a10) + zt*(zt*a02 + a01) + f00);
								}
							}
							else
							{//variable-step rectangular mesh
								double DX = 1./((rx_10 - rx_m10)*rx_10*rx_m10), DZ = 1./((rz_01 - rz_0m1)*rz_01*rz_0m1);
								if(arEX != 0)
								{
									f0m1 = arEX[ofst_0m1]; fm10 = arEX[ofst_m10]; f00 = arEX[ofst_00]; f10 = arEX[ofst_10]; f01 = arEX[ofst_01];
									double f00_mi_fm10_x1 = (f00 - fm10)*rx_10, f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									double f00_mi_f10_xm1 = (f00 - f10)*rx_m10, f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*t_ExRes = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);

									f0m1 = arEX[ofst_0m1_p_1]; fm10 = arEX[ofst_m10_p_1]; f00 = arEX[ofst_00_p_1]; f10 = arEX[ofst_10_p_1]; f01 = arEX[ofst_01_p_1];
									f00_mi_fm10_x1 = (f00 - fm10)*rx_10; f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									f00_mi_f10_xm1 = (f00 - f10)*rx_m10; f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*(t_ExRes + 1) = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);
								}
								if(arEZ != 0)
								{
									f0m1 = arEZ[ofst_0m1]; fm10 = arEZ[ofst_m10]; f00 = arEZ[ofst_00]; f10 = arEZ[ofst_10]; f01 = arEZ[ofst_01];
									double f00_mi_fm10_x1 = (f00 - fm10)*rx_10, f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									double f00_mi_f10_xm1 = (f00 - f10)*rx_m10, f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*t_EzRes = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);

									f0m1 = arEZ[ofst_0m1_p_1]; fm10 = arEZ[ofst_m10_p_1]; f00 = arEZ[ofst_00_p_1]; f10 = arEZ[ofst_10_p_1]; f01 = arEZ[ofst_01_p_1];
									f00_mi_fm10_x1 = (f00 - fm10)*rx_10; f00_mi_f0m1_z1 = (f00 - f0m1)*rz_01;
									f00_mi_f10_xm1 = (f00 - f10)*rx_m10; f00_mi_f01_zm1 = (f00 - f01)*rz_0m1;
									a10 = DX*(f00_mi_f10_xm1*rx_m10 - f00_mi_fm10_x1*rx_10); a01 = DZ*(f00_mi_f01_zm1*rz_0m1 - f00_mi_f0m1_z1*rz_01);
									a20 = DX*(f00_mi_fm10_x1 - f00_mi_f10_xm1); a02 = DZ*(f00_mi_f0m1_z1 - f00_mi_f01_zm1);
									*(t_EzRes + 1) = (float)(dx_00*(dx_00*a20 + a10) + dz_00*(dz_00*a02 + a01) + f00);
								}
							}
						}
						else
						{//irregular mesh (general case)
							double arXZ[] = {rx_0m1, rz_0m1, rx_m10, rz_m10, rx_10, rz_10, rx_01, rz_01};
							if(arEX != 0)
							{
								double arER[] = {arEX[ofst_0m1], arEX[ofst_m10], arEX[ofst_00], arEX[ofst_10], arEX[ofst_01]};
								*t_ExRes = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arER);
								double arEI[] = {arEX[ofst_0m1_p_1], arEX[ofst_m10_p_1], arEX[ofst_00_p_1], arEX[ofst_10_p_1], arEX[ofst_01_p_1]};
								*(t_ExRes + 1) = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arEI);
							}
							if(arEZ != 0)
							{
								double arER[] = {arEZ[ofst_0m1], arEZ[ofst_m10], arEZ[ofst_00], arEZ[ofst_10], arEZ[ofst_01]};
								*t_EzRes = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arER);
								double arEI[] = {arEZ[ofst_0m1_p_1], arEZ[ofst_m10_p_1], arEZ[ofst_00_p_1], arEZ[ofst_10_p_1], arEZ[ofst_01_p_1]};
								*(t_EzRes + 1) = (float)CGenMathInterp::Interp2dBiQuad5Var(dx_00, dz_00, arXZ, arEI);
							}
						}
					}
				}
				else
				{
					SetFieldToZero:
					*t_ExRes = 0.; *(t_ExRes+1) = 0.;
					*t_EzRes = 0.; *(t_EzRes+1) = 0.;

					ix0 = -1; iz0 = -1; //OC200414 (forces to restart search)
				}

				//test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				//long ofstAux = iz*PerZ + ix*PerX + ie*2;
				//*t_ExRes = arEX[ofstAux]; *(t_ExRes+1) = arEX[ofstAux+1];
				//*t_EzRes = arEZ[ofstAux]; *(t_EzRes+1) = arEZ[ofstAux+1];
				//end test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				t_ExRes += 2;
				t_EzRes += 2;

				phEn += pWfr->eStep;
			}
			x += pWfr->xStep;
		}
		z += pWfr->zStep;
	}

	if(waveFrontTermWasTreated) TreatStronglyOscillatingTerm(*pWfr, 'a');
	return 0;
}

//*************************************************************************

void srTMirror::RadPointModifier_ThinElem(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
{
	TVector3d inP(EXZ.x, EXZ.z, -m_extAlongOptAxIn), intersP;
	inP = TransHndl.rep->TrPoint_inv(inP); //Point in input transverse plane in local frame

	if(!FindRayIntersectWithSurfInLocFrame(inP, m_vInLoc, intersP)) 
	{
		*(EPtrs.pExIm) = 0.; *(EPtrs.pExRe) = 0.; *(EPtrs.pEzIm) = 0.; *(EPtrs.pEzRe) = 0.;
		return;
	}
	if(!CheckIfPointIsWithinOptElem(intersP.x, intersP.y)) 
	{
		*(EPtrs.pExIm) = 0.; *(EPtrs.pExRe) = 0.; *(EPtrs.pEzIm) = 0.; *(EPtrs.pEzRe) = 0.;
		return;
	}

	//Distance from intersection point with surface to intersection point with output transverse plane of this step
	double distBwIntersPtAndOut = m_vOutLoc*(m_vPtOutLoc - intersP);
	double distBwInAndIntersP = m_vInLoc*(intersP - inP);

	double optPathDif = distBwInAndIntersP + distBwIntersPtAndOut - (m_extAlongOptAxIn + m_extAlongOptAxOut);
	double phShift = 5.067730652e+06*EXZ.e*optPathDif; //to check sign!
	float cosPh, sinPh;
	CosAndSin(phShift, cosPh, sinPh);

	if(m_reflData.pData == 0) //no reflectivity defined
	{
		float NewExRe = (*(EPtrs.pExRe))*cosPh - (*(EPtrs.pExIm))*sinPh;
		float NewExIm = (*(EPtrs.pExRe))*sinPh + (*(EPtrs.pExIm))*cosPh;
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
		float NewEzRe = (*(EPtrs.pEzRe))*cosPh - (*(EPtrs.pEzIm))*sinPh;
		float NewEzIm = (*(EPtrs.pEzRe))*sinPh + (*(EPtrs.pEzIm))*cosPh;
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
		return;
	}

	//Calculate change of the electric field due to reflectivity...
	TVector3d vNormAtP;
	FindSurfNormalInLocFrame(intersP.x, intersP.y, vNormAtP);
	vNormAtP = TransHndl.rep->TrBiPoint(vNormAtP);  //to the frame of incident beam

	//double xp = (EXZ.x - m_inWfrCh)/m_inWfrRh, yp = (EXZ.z - m_inWfrCv)/m_inWfrRv;
	//TVector3d vRay(0, 0, 1.); //in the frame of incident beam
	//TVector3d vSig = vNormAtP^vRay; vSig.Normalize(); //in the frame of incident beam

	double vSigX=1., vSigY=0., vPiX=0., vPiY=1.;
	if((vNormAtP.x != 0.) || (vNormAtP.y != 0.))
	{
		//double multNorm = sqrt(vNormAtP.x*vNormAtP.x + vNormAtP.y*vNormAtP.y);
		double multNorm = 1./sqrt(vNormAtP.x*vNormAtP.x + vNormAtP.y*vNormAtP.y); //?
		vSigX = vNormAtP.y*multNorm; vSigY = -vNormAtP.x*multNorm;
		vPiX = -vSigY; vPiY = vSigX;
	}

	//TVector3d vEr(*(EPtrs.pExRe), *(EPtrs.pEzRe), 0), vEi(*(EPtrs.pExIm), *(EPtrs.pEzIm), 0);
	//Maybe rather this? //TVector3d vEr(-*(EPtrs.pExRe), *(EPtrs.pEzRe), 0), vEi(-*(EPtrs.pExIm), *(EPtrs.pEzIm), 0);
	//double EsigRe = vEr*vSig, EsigIm = vEi*vSig; //in the frame of incident beam
	//double EpiRe = vEr*vPi, EpiIm = vEi*vPi;
	//Sigma and Pi components of input electric field in the frame of incident beam
	double EsigRe = (*(EPtrs.pExRe))*vSigX + (*(EPtrs.pEzRe))*vSigY;
	double EsigIm = (*(EPtrs.pExIm))*vSigX + (*(EPtrs.pEzIm))*vSigY;
	double EpiRe = (*(EPtrs.pExRe))*vPiX + (*(EPtrs.pEzRe))*vPiY;
	double EpiIm = (*(EPtrs.pExIm))*vPiX + (*(EPtrs.pEzIm))*vPiY;

	//double sinAngInc = ::fabs(vRay*vNormAtP);
	double sinAngInc = ::fabs(vNormAtP.z);
	double angInc = asin(sinAngInc);
	double RsigRe=1, RsigIm=0, RpiRe=1, RpiIm=0;
	GetComplexReflectCoefFromTable(EXZ.e, angInc, RsigRe, RsigIm, RpiRe, RpiIm);

	double newEsigRe = -cosPh*(EsigIm*RsigIm - EsigRe*RsigRe) - sinPh*(EsigRe*RsigIm + EsigIm*RsigRe);
	double newEsigIm = cosPh*(EsigRe*RsigIm + EsigIm*RsigRe) - sinPh*(EsigIm*RsigIm - EsigRe*RsigRe);
	double newEpiRe = -cosPh*(EpiIm*RpiIm - EpiRe*RpiRe) - sinPh*(EpiRe*RpiIm + EpiIm*RpiRe);
	double newEpiIm = cosPh*(EpiRe*RpiIm + EpiIm*RpiRe) - sinPh*(EpiIm*RpiIm - EpiRe*RpiRe);

	//vEr = newEsigRe*vSig + newEpiRe*vPi; //in the frame of incident beam
	//vEi = newEsigIm*vSig + newEpiIm*vPi;
	double vErX = newEsigRe*vSigX + newEpiRe*vPiX; //in the frame of incident beam
	double vErY = newEsigRe*vSigY + newEpiRe*vPiY;
	double vEiX = newEsigIm*vSigX + newEpiIm*vPiX;
	double vEiY = newEsigIm*vSigY + newEpiIm*vPiY;

	//electric field components in the frame of output beam 
	//test
	//*(EPtrs.pExRe) = (float)(vEr*m_vHorOutIn);
	//*(EPtrs.pExIm) = (float)(vEi*m_vHorOutIn);
	//*(EPtrs.pEzRe) = (float)(vEr*m_vVerOutIn);
	//*(EPtrs.pEzIm) = (float)(vEi*m_vVerOutIn);
	*(EPtrs.pExRe) = (float)(vErX*m_vHorOutIn.x + vErY*m_vHorOutIn.y);
	*(EPtrs.pExIm) = (float)(vEiX*m_vHorOutIn.x + vEiY*m_vHorOutIn.y);
	*(EPtrs.pEzRe) = (float)(vErX*m_vVerOutIn.x + vErY*m_vVerOutIn.y);
	*(EPtrs.pEzIm) = (float)(vEiX*m_vVerOutIn.x + vEiY*m_vVerOutIn.y);
}

//*************************************************************************

int srTMirror::PropagateRadiationSimple_LocRayTracing(srTSRWRadStructAccessData* pRadAccessData)
{
	int res = 0;
	//char LocWaveFrontTermCanBeTreated = WaveFrontTermCanBeTreated(*pRadAccessData); //checks if quad. term can be treated and set local variables

	const double relTolAstigm = 1.e-06; //to steer (used for calculating directions of inpiut rays)

	if((m_treatInOut == 2) && (m_extAlongOptAxIn != 0.))
	{//Propagate wavefront back (by -m_extAlongOptAxIn) to the beginning of the optical element using Wavefront Propagation through a Drift
		srTRadResizeVect dummyResizeVect; //consider removing this completely
		srTDriftSpace driftIn(-m_extAlongOptAxIn);
		driftIn.PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;
		if(res = driftIn.PropagateRadiation(pRadAccessData, m_ParPrecWfrPropag, dummyResizeVect)) return res;
	}

	double RxInWfr = pRadAccessData->RobsX;
	double xcInWfr = pRadAccessData->xc;
	double RzInWfr = pRadAccessData->RobsZ;
	double zcInWfr = pRadAccessData->zc;

		//OCtest!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//TreatStronglyOscillatingTerm(*pRadAccessData, 'r');
		//return 0;
		//end test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(((m_treatInOut == 0) || (m_treatInOut == 2)) && (m_extAlongOptAxIn != 0))
	{
		srTDriftSpace driftAux(m_extAlongOptAxIn);
		driftAux.PropagateWaveFrontRadius(pRadAccessData);
	}
	if(res = PropagateWaveFrontRadius(pRadAccessData)) return res;
	if(((m_treatInOut == 0) || (m_treatInOut == 2)) && (m_extAlongOptAxOut != 0))
	{
		srTDriftSpace driftAux(m_extAlongOptAxOut);
		driftAux.PropagateWaveFrontRadius(pRadAccessData);
	}

	double RxOutWfr = pRadAccessData->RobsX;
	double RzOutWfr = pRadAccessData->RobsZ;

	//double ampFact = 1.;
	//double RxInCor = RxInWfr + m_extAlongOptAxIn;
	//double RzInCor = RzInWfr + m_extAlongOptAxIn;
	//double RxOutCor = RxOutWfr - m_extAlongOptAxOut;
	//double RzOutCor = RzOutWfr - m_extAlongOptAxOut;
	//if((RxInCor != 0.) && (RzInCor != 0.) && (RxOutWfr != 0.) && (RzOutCor != 0.))
	//{
	//	double ampFactE2 = RxInWfr*RxOutCor/(RxInCor*RxOutWfr);
	//	ampFactE2 *= RzInWfr*RzOutCor/(RzInCor*RzOutWfr);
	//	ampFact = sqrt(fabs(ampFactE2));
	//}
	double ampFact, ampFactE2, RxInCor, RzInCor, RxOutCor, RzOutCor;

	gmTrans *pTrans = TransHndl.rep;

	TVector3d rayLocFr[2]; //ray[2], , RayOut[2], arIntersectP[3];
	TVector3d &rayLocFrP = rayLocFr[0], &rayLocFrV = rayLocFr[1];
	TVector3d vIntersPtLocFr, vSurfNormLocFr;
	TVector3d vAuxIntersectP, vAuxOptPath, vRayIn, vSig, vPi, vTrAux;

	TVector3d planeBeforeLocFr[2]; // vAuxIntersectP, vAuxDif;
	TVector3d &planeBeforeLocFrP = planeBeforeLocFr[0], &planeBeforeLocFrV = planeBeforeLocFr[1];
	
	planeBeforeLocFrP.x = TransvCenPoint.x;
	planeBeforeLocFrP.y = TransvCenPoint.y;
	//planeBeforeLocFrP.z = 0.;
	//if(m_treatInOut == 1) planeBeforeLocFrP.z = -m_extAlongOptAxIn;
	planeBeforeLocFrP.z = -m_extAlongOptAxIn;

	planeBeforeLocFrV.x = planeBeforeLocFrV.y = 0.; 
	planeBeforeLocFrV.z = 1.;
	if(pTrans != 0)
	{
		planeBeforeLocFrP = pTrans->TrPoint_inv(planeBeforeLocFrP);
		planeBeforeLocFrV = pTrans->TrBiPoint_inv(planeBeforeLocFrV);
	}

	TVector3d planeAfterLocFr[2]; //point and normal vector (in the direction of beam propagation) to the plane at the exit of the optical element
	TVector3d &planeAfterLocFrP = planeAfterLocFr[0], &planeAfterLocFrV = planeAfterLocFr[1];

	TVector3d planeCenOutLocFr[2]; //point and normal vector (in the direction of beam propagation) to the plane at the center of the optical element
	TVector3d &planeCenOutLocFrP = planeCenOutLocFr[0], &planeCenOutLocFrV = planeCenOutLocFr[1];

	//Determine Exit (output) plane coordinates in the Local frame
	if(!FindRayIntersectWithSurfInLocFrame(planeBeforeLocFrP, m_vInLoc, planeAfterLocFrP)) return FAILED_DETERMINE_OPTICAL_AXIS;
	planeAfterLocFrV = m_vOutLoc;

			//OCTEST
			//TVector3d vpTestFoc1Loc = -0.2*m_vInLoc;
			//TVector3d vpTestFoc2Loc = 0.2*m_vOutLoc;
			//END OCTEST

	planeCenOutLocFrV = m_vOutLoc;
	planeCenOutLocFrP = planeAfterLocFrP;
	planeAfterLocFrP += m_extAlongOptAxOut*m_vOutLoc;

	float *pEX0 = pRadAccessData->pBaseRadX;
	float *pEZ0 = pRadAccessData->pBaseRadZ;
	double ePh = pRadAccessData->eStart, x, y;

	//long HalfPerX =  pRadAccessData->ne; //OC18032016
	//long PerX = HalfPerX << 1;
	//long HalfPerY = HalfPerX*pRadAccessData->nx; //OC18032016
	//long PerY = PerX*pRadAccessData->nx;
	//long nTot = PerY*pRadAccessData->nz;
	long long HalfPerX =  pRadAccessData->ne; //OC18032016
	long long PerX = HalfPerX << 1;
	long long HalfPerY = HalfPerX*pRadAccessData->nx; //OC18032016
	long long PerY = PerX*pRadAccessData->nx;
	long long nTot = PerY*pRadAccessData->nz;

	//float *arAuxRayTrCoord = new float[nTot];
	double *arAuxRayTrCoord = new double[nTot];
	if(arAuxRayTrCoord == 0) return NOT_ENOUGH_MEMORY_FOR_SR_COMP;

	//long half_nTot = nTot >> 1;
	//long *arAuxIndRayTrCoord = new long[half_nTot]; //OC18032016
	long long half_nTot = nTot >> 1;
	long long *arAuxIndRayTrCoord = new long long[half_nTot]; //OC18032016
	if(arAuxIndRayTrCoord == 0) return NOT_ENOUGH_MEMORY_FOR_SR_COMP;
	//long *t_arAuxIndRayTrCoord = arAuxIndRayTrCoord;
	//for(long j=0; j<half_nTot; j++) *(t_arAuxIndRayTrCoord++) = -1; //OC18032016
	long long *t_arAuxIndRayTrCoord = arAuxIndRayTrCoord;
	for(long long j=0; j<half_nTot; j++) *(t_arAuxIndRayTrCoord++) = -1; //OC18032016

	float *arAuxEX=0, *arAuxEY=0;
	if(pEX0 != 0) 
	{
		arAuxEX = new float[nTot];
		if(arAuxEX == 0) return NOT_ENOUGH_MEMORY_FOR_SR_COMP;
	}
	if(pEZ0 != 0)
	{
		arAuxEY = new float[nTot];
		if(arAuxEY == 0) return NOT_ENOUGH_MEMORY_FOR_SR_COMP;
	}

	double xRelOutMin = 1.E+23, xRelOutMax = -1.E+23;
	double yRelOutMin = 1.E+23, yRelOutMax = -1.E+23;

	long ixMin = pRadAccessData->nx - 1, ixMax = 0;
	long iyMin = pRadAccessData->nz - 1, iyMax = 0;
	//float cosPh, sinPh;
	double cosPh, sinPh;
	double EsigRe, EsigIm, EpiRe, EpiIm;
	double grMult;

	double dxOutMin = 1.e+23*(pRadAccessData->nx)*(pRadAccessData->xStep); //OC20082018
	double dxOutMax = 0.;
	double dyOutMin = 1.e+23*(pRadAccessData->nz)*(pRadAccessData->zStep);
	double dyOutMax = 0.;
	double xRelOutPrev = 1.e+23, yRelOutPrev = 1.e+23;

	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		double TwoPi_d_LambdaM = ePh*5.067730652e+06;
		long Two_ie = ie << 1;

		if(m_isGrating) grMult = m_grM/(806554.3835*ePh);
	
		bool coordOutFoundY = false;
		y = pRadAccessData->zStart;
		for(long iy=0; iy<pRadAccessData->nz; iy++)
		{
			//long iyPerY = iy*PerY;
			//long iyHalfPerY_p_ie = iy*HalfPerY + ie; //OC18032016
			long long iyPerY = iy*PerY;
			long long iyHalfPerY_p_ie = iy*HalfPerY + ie; //OC18032016
			float *pEX_StartForX = pEX0 + iyPerY;
			float *pEZ_StartForX = pEZ0 + iyPerY;

			float *pEX_StartForXres = arAuxEX + iyPerY;
			float *pEY_StartForXres = arAuxEY + iyPerY;

			//float *pAuxRayTrCoord = arAuxRayTrCoord + iyPerY;
			double *pAuxRayTrCoord = arAuxRayTrCoord + iyPerY;

			//bool firstHitForThisY = true; //OC19032016

			bool coordOutFoundX = false;
			x = pRadAccessData->xStart;
			for(long ix=0; ix<pRadAccessData->nx; ix++)
			{
				//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				long long ixPerX_p_Two_ie = ix*PerX + Two_ie;
				float *pExRe = pEX_StartForX + ixPerX_p_Two_ie;
				float *pExIm = pExRe + 1;
				float *pEzRe = pEZ_StartForX + ixPerX_p_Two_ie;
				float *pEzIm = pEzRe + 1;

				float *pExReRes = pEX_StartForXres + ixPerX_p_Two_ie;
				float *pExImRes = pExReRes + 1;
				float *pEyReRes = pEY_StartForXres + ixPerX_p_Two_ie;
				float *pEyImRes = pEyReRes + 1;
				//float *pAuxRayTrCoordX = pAuxRayTrCoord + ixPerX_p_Two_ie;
				//float *pAuxRayTrCoordY = pAuxRayTrCoordX + 1;
				double *pAuxRayTrCoordX = pAuxRayTrCoord + ixPerX_p_Two_ie;
				double *pAuxRayTrCoordY = pAuxRayTrCoordX + 1;

				//*pAuxRayTrCoordX = (float)(-1.E+23); *pAuxRayTrCoordY = (float)(-1.E+23);
				*pAuxRayTrCoordX = -1.E+23; *pAuxRayTrCoordY = -1.E+23;

				//long *pAuxIndRayTrCoord = arAuxIndRayTrCoord + iyHalfPerY_p_ie + ix*HalfPerX;
				//*pAuxIndRayTrCoord = -1;

				bool ExIsNotZero = false, EzIsNotZero = false;
				if(pEX0 != 0)
				{
					*pExReRes = 0.; *pExImRes = 0.;
					if((*pExRe != 0) || (*pExIm != 0)) ExIsNotZero = true;
				}
				if(pEZ0 != 0)
				{
					*pEyReRes = 0.; *pEyImRes = 0.;
					if((*pEzRe != 0) || (*pEzIm != 0)) EzIsNotZero = true;
				}
				if(ExIsNotZero || EzIsNotZero)
				{
					//double tgAngX=0, tgAngY=0;
					//pRadAccessData->GetWaveFrontNormal(x, y, tgAngX, tgAngY);

					//OC20092017 (commented-out)
					//double tgAngX = (x - xcInWfr)/RxInWfr; //check sign
					//double tgAngY = (y - zcInWfr)/RzInWfr; //check sign

					//rayLocFrV.x = tgAngX;
					//rayLocFrV.y = tgAngY;
					//rayLocFrV.z = sqrt(1. - rayLocFrV.x*rayLocFrV.x - rayLocFrV.y*rayLocFrV.y);

					//OC20092017
					double x_mi_xc = x - xcInWfr, y_mi_yc = y - zcInWfr;
					rayLocFrV.x = x_mi_xc;
					rayLocFrV.y = y_mi_yc;
					rayLocFrV.z = RxInWfr;
					rayLocFrV.Normalize();

					if(fabs(RxInWfr - RzInWfr) > fabs(RxInWfr)*relTolAstigm)
					{//astigmatism
						double auxVx = rayLocFrV.x;
						rayLocFrV.x = x_mi_xc;
						rayLocFrV.y = y_mi_yc;
						rayLocFrV.z = RzInWfr;
						rayLocFrV.Normalize();
						double auxVy = rayLocFrV.y;
						rayLocFrV.x = auxVx;
						rayLocFrV.z = sqrt(1. - auxVx*auxVx - auxVy*auxVy);
					}

						//OCTEST
						//TVector3d vInRayInFr(x - xcInWfr, y - zcInWfr, RzInWfr);
						//vInRayInFr.Normalize();
						//double nxTest = sin(atan(tgAngX));
						//double nyTest = sin(atan(tgAngY));
						//double nzTest = sqrt(1. - nxTest*nxTest - nyTest*nyTest);
						//TVector3d vInRayInFrTest(nxTest, nyTest, nzTest);
						//END OCTEST

					rayLocFrP.x = x; rayLocFrP.y = y; rayLocFrP.z = 0.;

					if((m_treatInOut == 0) || (m_treatInOut == 2))
					{
						rayLocFrP.z = -m_extAlongOptAxIn; //?
					}

					if(pTrans != 0)
					{//from input beam frame to local frame
						rayLocFrP = pTrans->TrPoint_inv(rayLocFrP);
						rayLocFrV = pTrans->TrBiPoint_inv(rayLocFrV);
					}
					vRayIn = rayLocFrV;

					//if((m_treatIn == 1) && (m_extAlongOptAxIn != 0.)) //check sign?
					//{//propagate back to a plane before optical element, using geometrical ray-tracing
					//	FindLineIntersectWithPlane(planeBeforeLocFr, rayLocFr, vAuxIntersectP);
					//	rayLocFrP = vAuxIntersectP;
					//}

					//bool intersectHappened = false;
					if(FindRayIntersectWithSurfInLocFrame(rayLocFrP, rayLocFrV, vIntersPtLocFr, &vSurfNormLocFr))
					{
						if(CheckIfPointIsWithinOptElem(vIntersPtLocFr.x, vIntersPtLocFr.y)) 
						{//continue calculating reflected and propagated electric field
									//OCTEST
									//FindRayIntersectWithSurfInLocFrame(rayLocFrP, rayLocFrV, vIntersPtLocFr, &vSurfNormLocFr);
									//intersectHappened = true;
									//END OCTEST
									//OCTEST
									//TVector3d vTestToFoc1 = vpTestFoc1Loc - vIntersPtLocFr;
									//TVector3d vTestToFoc2 = vpTestFoc2Loc - vIntersPtLocFr;
									//double distFoc1 = vTestToFoc1.Abs(), distFoc2 = vTestToFoc2.Abs();
									//double distFoc12 = distFoc1 + distFoc2;
									//vTestToFoc1.Normalize();
									//vTestToFoc2.Normalize();
									//END OCTEST

							vAuxOptPath = vIntersPtLocFr - rayLocFrP;
							//double optPath = vAuxOptPath.Abs();

							double optPath = vAuxOptPath*rayLocFrV;
							double optPathBefore = optPath;

							//Finding the Ray after the reflection (in local frame):
							rayLocFrP = vIntersPtLocFr;
							
							ampFact = 1.; //OC100314
							double phShiftGr = 0.;
							if(m_isGrating)
							{
								//Tangential vector perpendicular to grooves
								TVector3d vTang(vSurfNormLocFr.z*m_grAuxCosAng, vSurfNormLocFr.z*m_grAuxSinAng, -(vSurfNormLocFr.x*m_grAuxCosAng + vSurfNormLocFr.y*m_grAuxSinAng));
								vTang.Normalize();
								double xGr = vIntersPtLocFr.x;
								double locGrDen = m_grDen + xGr*(xGr*(xGr*(xGr*m_grDen4 + m_grDen3) + m_grDen2) + m_grDen1); //Calculate local Groove Density
								//OCTEST
								//double locGrDen = m_grDen; // + xGr*(xGr*(xGr*(xGr*m_grDen4 + m_grDen3) + m_grDen2) + m_grDen1); //Calculate local Groove Density
								//vTang *= (grMult*locGrDen);
								//END OCTEST
								vTang *= (-grMult*locGrDen);

								TVector3d vInLocTang = rayLocFrV - ((rayLocFrV*vSurfNormLocFr)*vSurfNormLocFr);
								TVector3d vOutLocTang = vInLocTang + vTang;
								double absE2_vOutLocTang = vOutLocTang.AmpE2();
								double abs_vOutLocNorm = sqrt(::fabs(1. - absE2_vOutLocTang));
								rayLocFrV = vOutLocTang + (abs_vOutLocNorm*vSurfNormLocFr);

								rayLocFrV.Normalize(); //required here?

								//Number of grooves from center to intersection point
								//double dN = xGr*(xGr*(xGr*(xGr*(xGr*0.2*m_grDen3 + 0.25*m_grDen3) + m_grDen2/3.) + 0.5*m_grDen1) + m_grDen);
								double dN = xGr*(xGr*(xGr*(xGr*(xGr*0.2*m_grDen4 + 0.25*m_grDen3) + m_grDen2/3.) + 0.5*m_grDen1) + m_grDen); //OC08022017
								phShiftGr = -6.283185307179586*dN*m_grM; //Check the sign!

								ampFact = m_grAuxElecFldAnamorphMagnFact;
							}
							else
							{
								rayLocFrV -= (2.*(rayLocFrV*vSurfNormLocFr))*vSurfNormLocFr; //Reflection Law (valid for mirrors only!)
								rayLocFrV.Normalize();
							}

							if((m_treatInOut == 0) || (m_treatInOut == 2))
							{
								FindLineIntersectWithPlane(planeAfterLocFr, rayLocFr, vAuxIntersectP);
							}
							else if(m_treatInOut == 1)
							{
								FindLineIntersectWithPlane(planeCenOutLocFr, rayLocFr, vAuxIntersectP);

								//OCTEST
								//TVector3d vTestAuxIntersectP, rayTestLocFr[2];
								//rayTestLocFr[0] = rayLocFr[0]; rayTestLocFr[1] = m_vOutLoc;
								//FindLineIntersectWithPlane(planeCenOutLocFr, rayTestLocFr, vTestAuxIntersectP);
								//vAuxOptPath = vTestAuxIntersectP - rayLocFrP;
								//auxPathAfter = vAuxOptPath*rayTestLocFr[1];
								//END OCTEST
							}

							vAuxOptPath = vAuxIntersectP - rayLocFrP;
							//optPath += vAuxOptPath.Abs();
							double optPathAfter = vAuxOptPath*rayLocFrV;
							optPath += optPathAfter;

							//double RxInCor = (RxInWfr > 0)? (RxInWfr + optPathBefore) : (RxInWfr - optPathBefore);
							//double RzInCor = (RzInWfr > 0)? (RzInWfr + optPathBefore) : (RzInWfr - optPathBefore);
							//double RxOutCor = (RxOutWfr > 0)? (RxOutWfr + optPathAfter) : (RxOutWfr - optPathAfter);
							//double RzOutCor = (RzOutWfr > 0)? (RzOutWfr + optPathAfter) : (RzOutWfr - optPathAfter);

							//ampFact = 1.; //OC100314
							RxInCor = RxInWfr + optPathBefore; //to check signs
							RzInCor = RzInWfr + optPathBefore;
							RxOutCor = RxOutWfr - optPathAfter;
							RzOutCor = RzOutWfr - optPathAfter;
							if((RxInCor != 0.) && (RzInCor != 0.) && (RxOutWfr != 0.) && (RzOutCor != 0.))
							{//OC: this may require more tests/debugging
								ampFactE2 = RxInWfr*RxOutCor/(RxInCor*RxOutWfr);
								ampFactE2 *= RzInWfr*RzOutCor/(RzInCor*RzOutWfr);
								//ampFact = sqrt(fabs(ampFactE2));
								ampFact *= sqrt(fabs(ampFactE2)); //OC100314
							}

							//Calculating transverse coordinates of intersection point of the ray with the output plane (or central plane) in the frame of the output beam
							vTrAux = vAuxIntersectP - planeCenOutLocFrP;
							if(pTrans != 0)
							{//from local frame to input beam frame
								vTrAux = pTrans->TrBiPoint(vTrAux);
							}
							//float xRelOut = (float)(vTrAux*m_vHorOutIn);
							//float yRelOut = (float)(vTrAux*m_vVerOutIn);
							double xRelOut = vTrAux*m_vHorOutIn; //OC18032016
							double yRelOut = vTrAux*m_vVerOutIn;
							//test!!!!!!!!!!!!!!!!!!!!!
							//float yRelOut = -(float)(vTrAux*m_vVerOutIn);
							//end test!!!!!!!!!!!!!!!!!!!!!

							*pAuxRayTrCoordX = xRelOut;
							*pAuxRayTrCoordY = yRelOut;

							if(coordOutFoundX) //OC20082018
							{
								double dxOut = fabs(xRelOut - xRelOutPrev);
								if(dxOutMin > dxOut) dxOutMin = dxOut;
								else if(dxOutMax < dxOut) dxOutMax = dxOut;
							}
							if(coordOutFoundY) //OC20082018
							{
								double dyOut = fabs(yRelOut - yRelOutPrev);
								if(dyOutMin > dyOut) dyOutMin = dyOut;
								else if(dyOutMax < dyOut) dyOutMax = dyOut;
							}
							xRelOutPrev = xRelOut;
							yRelOutPrev = yRelOut;
							coordOutFoundX = true;
							coordOutFoundY = true;

									//OCTEST
									//if(y > 0)
									//if(y > 0.0043)
									//{
									//	int aha = 1;
									//}
									//if((-0.5*(pRadAccessData->xStep) <= x) && (x <= 0.5*(pRadAccessData->xStep)))
									//{
									//	TVector3d auxTestPlane[2], vTest;
									//	auxTestPlane[0] = 1.8*m_vOutLoc; //point
									//	auxTestPlane[1] = m_vOutLoc; //normal vector
									//	FindLineIntersectWithPlane(auxTestPlane, rayLocFr, vTest);
									//	int aha = 1;
									//}
									//END OCTEST

							if(xRelOutMin > xRelOut) 
							{
								xRelOutMin = xRelOut; ixMin = ix;
							}
							if(xRelOutMax < xRelOut) 
							{
								xRelOutMax = xRelOut; ixMax = ix;
							}
							if(yRelOutMin > yRelOut) 
							{
								yRelOutMin = yRelOut; iyMin = iy;
							}
							if(yRelOutMax < yRelOut) 
							{
								yRelOutMax = yRelOut; iyMax = iy;
							}

							long ixRelOut = (long)((xRelOut - (pRadAccessData->xStart))/(pRadAccessData->xStep)); //OC18032016
							//if(ixRelOut < 0) ixRelOut = 0;
							long iyRelOut = (long)((yRelOut - (pRadAccessData->zStart))/(pRadAccessData->zStep)); //OC18032016
							//if(iyRelOut < 0) iyRelOut = 0;

							if((ixRelOut >= 0) && (ixRelOut < pRadAccessData->nx) && (iyRelOut >= 0) && (iyRelOut < pRadAccessData->nz)) //OC24032016
							{
								//long ofstTrfCoord = iyRelOut*HalfPerY + ixRelOut*HalfPerX + ie; //OC18032016 (to facilitate search of transformed coordinates)
								//long ofstToSet = iyPerY + ix*PerX + Two_ie;
								long long ofstTrfCoord = iyRelOut*HalfPerY + ixRelOut*HalfPerX + ie; //OC18032016 (to facilitate search of transformed coordinates)
								long long ofstToSet = iyPerY + ix*PerX + Two_ie;

								arAuxIndRayTrCoord[ofstTrfCoord] = ofstToSet;

									//OCTEST
									//if(m_isGrating && (ix==207))
									//{
									//	int aha = 1;
									//}
									//END OCTEST

								//last commented:
								//double phShift = TwoPi_d_LambdaM*optPathDif; //to check sign!
								//double phShift = TwoPi_d_LambdaM*optPath; //to check sign!

								double phShift = TwoPi_d_LambdaM*optPath + phShiftGr;

								//CosAndSin(phShift, cosPh, sinPh);
								cosPh = cos(phShift); sinPh = sin(phShift); //OC260114

								if(m_reflData.pData == 0) //no reflectivity defined
								//if(true) //no reflectivity defined
								{
									if(pEX0 != 0)
									{
										//float NewExRe = (float)(ampFact*((*pExRe)*cosPh - (*pExIm)*sinPh));
										//float NewExIm = (float)(ampFact*((*pExRe)*sinPh + (*pExIm)*cosPh));
										double NewExRe = ampFact*((*pExRe)*cosPh - (*pExIm)*sinPh); //OC260114
										double NewExIm = ampFact*((*pExRe)*sinPh + (*pExIm)*cosPh);

										//*pExReRes = NewExRe; *pExImRes = NewExIm;
										*pExReRes = (float)NewExRe; *pExImRes = (float)NewExIm;
									}
									if(pEZ0 != 0)
									{
										//float NewEzRe = (float)(ampFact*((*pEzRe)*cosPh - (*pEzIm)*sinPh));
										//float NewEzIm = (float)(ampFact*((*pEzRe)*sinPh + (*pEzIm)*cosPh));
										double NewEzRe = ampFact*((*pEzRe)*cosPh - (*pEzIm)*sinPh); //OC260114
										double NewEzIm = ampFact*((*pEzRe)*sinPh + (*pEzIm)*cosPh);

										//*pEyReRes = NewEzRe; *pEyImRes = NewEzIm;
										*pEyReRes = (float)NewEzRe; *pEyImRes = (float)NewEzIm;
									}

									//OCTEST
									//double Pi_d_Lambda_m = ePh*2.533840802E+06;
									//double xRel = x - TransvCenPoint.x, zRel = y - TransvCenPoint.y;

									//phShift = -Pi_d_Lambda_m*(xRel*xRel/FocDistX + zRel*zRel/FocDistZ);
									//CosAndSin(phShift, cosPh, sinPh);
									//float NewExRe = (*pExRe)*cosPh - (*pExIm)*sinPh;
									//float NewExIm = (*pExRe)*sinPh + (*pExIm)*cosPh;
									//*pExReRes = NewExRe; *pExImRes = NewExIm; 
									//float NewEzRe = (*pEzRe)*cosPh - (*pEzIm)*sinPh;
									//float NewEzIm = (*pEzRe)*sinPh + (*pEzIm)*cosPh;
									//*pEyReRes = NewEzRe; *pEyImRes = NewEzIm; 
									//END OCTEST

									//*pExRe = phShift; *pExIm = 0; 
									//*pEzRe = phShift; *pEzIm = 0; 
								}
								else
								//if(m_reflData.pData != 0)
								{//Calculate change of the electric field due to reflectivity...
									vRayIn.Normalize();
									vSig = (-1)*(vRayIn^vSurfNormLocFr); //sigma unit vector in Local frame; check sign
									double grazAng = 1.5707963268;
									if(vSig.isZero())
									{//In the frame of incident beam
										vSig.x = 1.; vSig.y = 0.; vSig.z = 0.;
										vPi.x = 0.; vPi.y = 1.; vPi.z = 0.;
									}
									else
									{
										//grazAng = asin(-(vRayIn*vSurfNormLocFr));
										grazAng = acos(vRayIn*vSurfNormLocFr) - 1.5707963267948966;

										vSig.Normalize();
										vPi = vRayIn^vSig;
										if(pTrans != 0)
										{//to the frame of incident beam
											vSig = pTrans->TrBiPoint(vSig);
											vPi = pTrans->TrBiPoint(vPi);
										}
									}

									EsigRe = EsigIm = EpiRe = EpiIm = 0.;
									if(pEX0 != 0)
									{
										EsigRe = (*pExRe)*vSig.x;
										EsigIm = (*pExIm)*vSig.x;
										EpiRe = (*pExRe)*vPi.x;
										EpiIm = (*pExIm)*vPi.x;
									}
									if(pEZ0 != 0)
									{
										EsigRe += (*pEzRe)*vSig.y;
										EsigIm += (*pEzIm)*vSig.y;
										EpiRe += (*pEzRe)*vPi.y;
										EpiIm += (*pEzIm)*vPi.y;
									}
									//double EsigRe = (*pExRe)*vSig.x + (*pEzRe)*vSig.y;
									//double EsigIm = (*pExIm)*vSig.x + (*pEzIm)*vSig.y;
									//double EpiRe = (*pExRe)*vPi.x + (*pEzRe)*vPi.y;
									//double EpiIm = (*pExIm)*vPi.x + (*pEzIm)*vPi.y;

									double RsigRe = 1, RsigIm = 0, RpiRe = 1, RpiIm = 0;
									GetComplexReflectCoefFromTable(ePh, grazAng, RsigRe, RsigIm, RpiRe, RpiIm);

									double newEsigRe = -cosPh*(EsigIm*RsigIm - EsigRe*RsigRe) - sinPh*(EsigRe*RsigIm + EsigIm*RsigRe);
									double newEsigIm = cosPh*(EsigRe*RsigIm + EsigIm*RsigRe) - sinPh*(EsigIm*RsigIm - EsigRe*RsigRe);
									double newEpiRe = -cosPh*(EpiIm*RpiIm - EpiRe*RpiRe) - sinPh*(EpiRe*RpiIm + EpiIm*RpiRe);
									double newEpiIm = cosPh*(EpiRe*RpiIm + EpiIm*RpiRe) - sinPh*(EpiIm*RpiIm - EpiRe*RpiRe);
									//double newEsigRe = -(EsigIm*RsigIm - EsigRe*RsigRe);
									//double newEsigIm = EsigRe*RsigIm + EsigIm*RsigRe;
									//double newEpiRe = -(EpiIm*RpiIm - EpiRe*RpiRe);
									//double newEpiIm = EpiRe*RpiIm + EpiIm*RpiRe;

									//In the frame of incident beam:
									double vErX = newEsigRe*vSig.x + newEpiRe*vPi.x;
									double vErY = newEsigRe*vSig.y + newEpiRe*vPi.y;
									double vEiX = newEsigIm*vSig.x + newEpiIm*vPi.x;
									double vEiY = newEsigIm*vSig.y + newEpiIm*vPi.y;

									//In the frame of output beam:
									if(pEX0 != 0)
									{
										//*pExRe = (float)(vErX*m_vHorOutIn.x + vErY*m_vHorOutIn.y);
										//*pExIm = (float)(vEiX*m_vHorOutIn.x + vEiY*m_vHorOutIn.y);
										*pExReRes = (float)(ampFact*(vErX*m_vHorOutIn.x + vErY*m_vHorOutIn.y));
										*pExImRes = (float)(ampFact*(vEiX*m_vHorOutIn.x + vEiY*m_vHorOutIn.y));
									}
									if(pEZ0 != 0)
									{
										//*pEzRe = (float)(vErX*m_vVerOutIn.x + vErY*m_vVerOutIn.y);
										//*pEzIm = (float)(vEiX*m_vVerOutIn.x + vEiY*m_vVerOutIn.y);
										*pEyReRes = (float)(ampFact*(vErX*m_vVerOutIn.x + vErY*m_vVerOutIn.y));
										*pEyImRes = (float)(ampFact*(vEiX*m_vVerOutIn.x + vEiY*m_vVerOutIn.y));
									}

									//OCTEST!!!!!!!!!!!!!!!!!!!!!
									//if(fabs(pRadAccessData->RobsZ + 18.079) < 0.1)
									//{
										//*pExReRes = *pExRe; *pExImRes = *pExIm;
										//*pEyReRes = *pEzRe; *pEyImRes = *pEzIm;
										//*pExReRes = EsigRe; *pExImRes = EsigIm;
										//*pEyReRes = EpiRe; *pEyImRes = EpiIm;
										//*pExReRes = newEsigRe; *pExImRes = newEsigIm;
										//*pEyReRes = newEpiRe; *pEyImRes = newEpiIm;
										//*pExReRes = phShift; *pExImRes = 0.;
										//*pEyReRes = phShift; *pEyImRes = 0.;
										//*pExReRes = xRelOut; *pExImRes = yRelOut;
										//*pEyReRes = xRelOut; *pEyImRes = yRelOut;
									//}
									//END OCTEST!!!!!!!!!!!!!!!!!!!!!
								}
							}

							//firstHitForThisY = false; //OC19032016
						}
					}
				}
				x += pRadAccessData->xStep;
			}
			y += pRadAccessData->zStep;
		}
		ePh += pRadAccessData->eStep;
	}

	//Re-interpolate the output wavefront (at fixed photon energy) on the initial equidistant grid:
	//if(res = WfrInterpolOnOrigGrid(pRadAccessData, arAuxRayTrCoord, arAuxEX, arAuxEY, xRelOutMin, xRelOutMax, yRelOutMin, yRelOutMax)) return res;
	//OCTEST
	//if(fabs(pRadAccessData->RobsZ + 18.079) > 0.1)
	//if(!m_isGrating)
	//{//OCTEST

	//if(res = WfrInterpolOnOrigGrid2(pRadAccessData, arAuxRayTrCoord, arAuxIndRayTrCoord, arAuxEX, arAuxEY, xRelOutMin, xRelOutMax, yRelOutMin, yRelOutMax)) return res;
	//OC20082018
	if(res = WfrInterpolOnOrigGrid2(pRadAccessData, arAuxRayTrCoord, arAuxIndRayTrCoord, arAuxEX, arAuxEY, xRelOutMin, xRelOutMax, yRelOutMin, yRelOutMax, dxOutMax, dyOutMax)) return res;

	//}
	//else
	//{
	//OCTEST
	//float *t_ExRes = pRadAccessData->pBaseRadX;
	//float *t_EyRes = pRadAccessData->pBaseRadZ;
	//float *t_arEX = arAuxEX, *t_arEY = arAuxEY;
	//for(long k=0; k<nTot; k++)
	//{
	//	*(t_ExRes++) = *(t_arEX++);
	//	*(t_EyRes++) = *(t_arEY++);
	//}
	//END OCTEST
	//}

	if((m_treatInOut == 2) && (m_extAlongOptAxOut != 0.))
	{//Propagate wavefront back (by -m_extAlongOptAxOut) to the center of the optical element using Wavefront Propagation through a Drift
		srTRadResizeVect dummyResizeVect; //consider removing this completely
		srTDriftSpace driftOut(-m_extAlongOptAxOut);
		driftOut.PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;
		if(res = driftOut.PropagateRadiation(pRadAccessData, m_ParPrecWfrPropag, dummyResizeVect)) return res;
	}

	//DEBUG
	//srTMomentsPtrs MomX(pRadAccessData->pMomX, 0);
	//ComputeRadMoments(pRadAccessData);
	//MomX.ComputeCentralMoments();
	//END DEBUG

	if(arAuxEX != 0) delete[] arAuxEX;
	if(arAuxEY != 0) delete[] arAuxEY;
	if(arAuxRayTrCoord != 0) delete[] arAuxRayTrCoord;
	if(arAuxIndRayTrCoord != 0) delete[] arAuxIndRayTrCoord; //OC18032016
	//if(arOptPathDif != 0) delete[] arOptPathDif;

	return 0;
}

//*************************************************************************
//Test of propagation by Fourier method in steps (failed?)
int srTMirror::PropagateRadiationSimple_FourierByParts(srTSRWRadStructAccessData* pRadAccessData)
{
	int res = 0;
	//propagate wavefront back (by -m_extAlongOptAxIn) to the beginning of the optical element
	//to make optional, assuming that the wavefront can be supplied already before the optical element, and not in its middle 
	srTRadResizeVect dummyResizeVect; //consider removing this completely
	srTDriftSpace driftIn(-m_extAlongOptAxIn);
	driftIn.PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;
	if(res = driftIn.PropagateRadiation(pRadAccessData, m_ParPrecWfrPropag, dummyResizeVect)) return res;

	//m_pRadAux = new srTSRWRadStructAccessData(pRadAccessData); //to propagate "old" wavefront part

	//if(res = PropagateWaveFrontRadius(pRadAccessData)) return res;
	//pRadAccessData->AssignElFieldToConst((float)0., (float)0.); //to place and propagate "new" wavefront parts
	//test
	//FocDistZ = 0.153; //test

	//propagate through the optical element by steps
	double stepProp = (m_extAlongOptAxIn + m_extAlongOptAxOut)/m_numPartsProp;
	srTDriftSpace driftStep(stepProp);
	driftStep.PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;

	m_longPosStartPropPart = 0.; m_longPosEndPropPart = stepProp;

	double posOnOutOptAx = m_extAlongOptAxOut - stepProp*(m_numPartsProp - 1);
	m_vPtOutLoc = posOnOutOptAx*m_vOutLoc;
	TVector3d vStepProp = stepProp*m_vOutLoc;

	for(int i=0; i<m_numPartsProp; i++)
	{
		m_inWfrRh = pRadAccessData->RobsX; m_inWfrRv = pRadAccessData->RobsZ;
		m_inWfrCh = pRadAccessData->xc; m_inWfrCv = pRadAccessData->zc;
		//m_inWfrRh = m_pRadAux->RobsX; m_inWfrRv = m_pRadAux->RobsZ;
		//m_inWfrCh = m_pRadAux->xc; m_inWfrCv = m_pRadAux->zc;
		if(res = TraverseRadZXE(pRadAccessData)) return res;

		if(res = driftStep.PropagateRadiation(pRadAccessData, m_ParPrecWfrPropag, dummyResizeVect)) return res;
		//if(res = driftStep.PropagateRadiation(m_pRadAux, m_ParPrecWfrPropag, dummyResizeVect)) return res;
			//if(i == 2) break;
			//break;

		m_longPosStartPropPart = m_longPosEndPropPart; m_longPosEndPropPart += stepProp;
		m_vPtOutLoc += vStepProp;
	}

	//srTDriftSpace driftOutCen(-m_extAlongOptAxOut);
	//driftOutCen.PropBufVars.UseExactRxRzForAnalytTreatQuadPhaseTerm = true;
	//if(res = driftOutCen.PropagateRadiation(pRadAccessData, m_ParPrecWfrPropag, dummyResizeVect)) return res;
	
	//delete m_pRadAux; m_pRadAux = 0; //memory leak is possible!
	return 0;
}

//*************************************************************************
//Test of propagation by Fourier method in steps (failed?)
void srTMirror::RadPointModifier_FourierByParts(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
{//adds optical path difference to simulate propagation through a part of optical element

	TVector3d inP(EXZ.x, EXZ.z, -m_extAlongOptAxIn + m_longPosStartPropPart);
	inP = TransHndl.rep->TrPoint_inv(inP); //point in input transverse plane in local frame
	
	TVector3d intersP;
	FindRayIntersectWithSurfInLocFrame(inP, m_vInLoc, intersP);

	if(!CheckIfPointIsWithinOptElem(intersP.x, intersP.y))
	{
		*(EPtrs.pExIm) = 0.; *(EPtrs.pExRe) = 0.; *(EPtrs.pEzIm) = 0.; *(EPtrs.pEzRe) = 0.;
		//if(m_pRadAux != 0)
		//{
		//	if(m_pRadAux->pBaseRadX != 0)
		//	{
		//		float *pEx = m_pRadAux->pBaseRadX + EXZ.aux_offset;
		//		*(pEx++) = 0; *(pEx++) = 0;
		//	}
		//	if(m_pRadAux->pBaseRadZ != 0)
		//	{
		//		float *pEz = m_pRadAux->pBaseRadZ + EXZ.aux_offset;
		//		*(pEz++) = 0; *(pEz++) = 0;
		//	}
		//}
		return;
	}

	//distance from intersection point with surface to intersection point with output transverse plane of this step
	double distBwIntersPtAndOut = m_vOutLoc*(m_vPtOutLoc - intersP);

	if(distBwIntersPtAndOut < 0)
	{
		//*(EPtrs.pExIm) = 0.; *(EPtrs.pExRe) = 0.; *(EPtrs.pEzIm) = 0.; *(EPtrs.pEzRe) = 0.;
		return;
	}

	//double dx = intersP.x - inP.x, dy = intersP.y - inP.y, dz = intersP.z - inP.z;
	//double distBwInAndIntersP = sqrt(dx*dx + dy*dy + dz*dz);

	double stepProp = m_longPosEndPropPart - m_longPosStartPropPart;

	double distBwInAndIntersP = m_vInLoc*(intersP - inP);
	if(distBwInAndIntersP < 0.) return; //no need to apply opt. path correction at this point, because it has been already applied

	if(distBwInAndIntersP > stepProp) 
	{
		//*(EPtrs.pExIm) = 0.; *(EPtrs.pExRe) = 0.; *(EPtrs.pEzIm) = 0.; *(EPtrs.pEzRe) = 0.;
		return; //no need to apply opt. path correction at this point: it will be applied at next step(s)
	}

	//add optical path difference, taking into account electric field transformation at mirror surface
	double optPathDif = distBwInAndIntersP + distBwIntersPtAndOut - stepProp;
	double phShift = 5.067730652e+06*EXZ.e*optPathDif; //to check sign!
	float cosPh, sinPh;
	CosAndSin(phShift, cosPh, sinPh);
	//test
	//cosPh = 1.; sinPh = 0.;

	//if(m_pRadAux != 0)
	//{
	//	if(m_pRadAux->pBaseRadX != 0)
	//	{
	//		float *pEx = m_pRadAux->pBaseRadX + EXZ.aux_offset;
	//		*(EPtrs.pExRe) = *pEx; *(pEx++) = 0;
	//		*(EPtrs.pExIm) = *pEx; *pEx = 0;
	//	}
	//	if(m_pRadAux->pBaseRadZ != 0)
	//	{
	//		float *pEy = m_pRadAux->pBaseRadZ + EXZ.aux_offset;
	//		*(EPtrs.pEzRe) = *pEy; *(pEy++) = 0;
	//		*(EPtrs.pEzIm) = *pEy; *pEy = 0;
	//	}
	//}

	TVector3d vEr(*(EPtrs.pExRe), *(EPtrs.pEzRe), 0), vEi(*(EPtrs.pExIm), *(EPtrs.pEzIm), 0);
	//Maybe rather this?
	//TVector3d vEr(-*(EPtrs.pExRe), *(EPtrs.pEzRe), 0), vEi(-*(EPtrs.pExIm), *(EPtrs.pEzIm), 0);

	double xp = (EXZ.x - m_inWfrCh)/m_inWfrRh, yp = (EXZ.z - m_inWfrCv)/m_inWfrRv;
	TVector3d vRay(xp, yp, sqrt(1. - xp*xp - yp*yp)); //in the frame of incident beam

	TVector3d vNormAtP;
	FindSurfNormalInLocFrame(intersP.x, intersP.y, vNormAtP);
	vNormAtP = TransHndl.rep->TrBiPoint(vNormAtP);  //to the frame of incident beam

	TVector3d vSig = vNormAtP^vRay; vSig.Normalize(); //in the frame of incident beam
	TVector3d vPi = vRay^vSig;
	double EsigRe = vEr*vSig, EsigIm = vEi*vSig; //in the frame of incident beam
	double EpiRe = vEr*vPi, EpiIm = vEi*vPi;

	//getting complex reflecivity coefficients for Sigma and Pi components of the electric field
	int ne = m_reflData.DimSizes[1];
	double eStart = m_reflData.DimStartValues[1];
	double eStep = m_reflData.DimSteps[1];
	int nAng = m_reflData.DimSizes[2];
	double angStart = m_reflData.DimStartValues[2];
	double angStep = m_reflData.DimSteps[2];

	const long perSigPi = 2;
	const long perPhotEn = perSigPi << 1;
	//long perAng = perPhotEn*ne;
	long long perAng = perPhotEn*ne;

	int ie = (int)((EXZ.e - eStart)/eStep + 0.00001);
	if((EXZ.e - (eStart + ie*eStep)) > 0.5*eStep) ie++;
	if(ie < 0) ie = 0;
	if(ie >= ne) ie = ne - 1;

	double sinAngInc = ::fabs(vRay*vNormAtP);
	double angInc = asin(sinAngInc);

	int iAng = (int)((angInc - angStart)/angStep + 0.00001);
	if((angInc - (angStart + iAng*angStep)) > 0.5*angStep) iAng++;
	if(iAng < 0) iAng = 0;
	if(iAng >= nAng) iAng = nAng - 1;

	//long ofstSig = perPhotEn*ie + perAng*iAng;
	long long ofstSig = perPhotEn*ie + perAng*iAng;
	//long ofstPi = ofstSig + perSigPi;
	double RsigRe=1, RsigIm=0, RpiRe=1, RpiIm=0;

	//setting appropriate pointer type 
	if(m_reflData.pData != 0)
	{
		if(m_reflData.DataType[1] == 'f')
		{
			float *pRsig = ((float*)(m_reflData.pData)) + ofstSig;
			float *pRpi = pRsig + perSigPi;
			RsigRe = *(pRsig++); RsigIm = *pRsig;
			RpiRe = *(pRpi++); RpiIm = *pRpi;
		}
		else
		{
			double *pRsig = ((double*)(m_reflData.pData)) + ofstSig;
			double *pRpi = pRsig + perSigPi;
			RsigRe = *(pRsig++); RsigIm = *pRsig;
			RpiRe = *(pRpi++); RpiIm = *pRpi;
		}
	}

	double newEsigRe = -cosPh*(EsigIm*RsigIm - EsigRe*RsigRe) - sinPh*(EsigRe*RsigIm + EsigIm*RsigRe);
	double newEsigIm = cosPh*(EsigRe*RsigIm + EsigIm*RsigRe) - sinPh*(EsigIm*RsigIm - EsigRe*RsigRe);
	double newEpiRe = -cosPh*(EpiIm*RpiIm - EpiRe*RpiRe) - sinPh*(EpiRe*RpiIm + EpiIm*RpiRe);
	double newEpiIm = cosPh*(EpiRe*RpiIm + EpiIm*RpiRe) - sinPh*(EpiIm*RpiIm - EpiRe*RpiRe);

	vEr = newEsigRe*vSig + newEpiRe*vPi; //in the frame of incident beam
	vEi = newEsigIm*vSig + newEpiIm*vPi;

	//electric field components in the frame of output beam 
	//test
	*(EPtrs.pExRe) = (float)(vEr*m_vHorOutIn);
	*(EPtrs.pExIm) = (float)(vEi*m_vHorOutIn);
	*(EPtrs.pEzRe) = (float)(vEr*m_vVerOutIn);
	*(EPtrs.pEzIm) = (float)(vEi*m_vVerOutIn);
}

//*************************************************************************

void srTMirror::EstimateFocalLengths(double radTan, double radSag) //to make it virtual in srTFocusingElem?
{//Assumes that m_vCenNorm, m_vCenTang are set !
 //Estimates focal lengths (approximately!):
	double cosAng = ::fabs(m_vCenNorm.z);
	if(::fabs(m_vCenTang.x) < ::fabs(m_vCenTang.y))
	{//tangential plane is close to be vertical
		if(::fabs(m_vCenNorm.x) < ::fabs(m_vCenNorm.y))
		{//normal is turned in vertical direction
			//if(FocDistX == 0.) FocDistX = 0.5*radSag/cosAng; //focal length in horizontal plane
			//if(FocDistZ == 0.) FocDistZ = 0.5*radTan*cosAng; //focal length in vertical plane
			FocDistX = 0.5*radSag/cosAng; //focal length in horizontal plane
			FocDistZ = 0.5*radTan*cosAng; //focal length in vertical plane
		}
		else
		{//normal is turned in horizontal direction
			//if(FocDistX == 0.) FocDistX = 0.5*radSag*cosAng; //focal length in horizontal plane
			//if(FocDistZ == 0.) FocDistZ = 0.5*radTan/cosAng; //focal length in vertical plane
			FocDistX = 0.5*radSag*cosAng; //focal length in horizontal plane
			FocDistZ = 0.5*radTan/cosAng; //focal length in vertical plane
		}
	}
	else
	{//tangential plane is close to be horizontal
		if(::fabs(m_vCenNorm.x) < ::fabs(m_vCenNorm.y))
		{//normal is turned in vertical direction
			//if(FocDistX == 0.) FocDistX = 0.5*radTan/cosAng; //focal length in vertical plane
			//if(FocDistZ == 0.) FocDistZ = 0.5*radSag*cosAng; //focal length in vertical plane
			FocDistX = 0.5*radTan/cosAng; //focal length in vertical plane
			FocDistZ = 0.5*radSag*cosAng; //focal length in vertical plane
		}
		else
		{//normal is turned in horizontal direction
			//if(FocDistX == 0.) FocDistX = 0.5*radTan*cosAng; //focal length in vertical plane
			//if(FocDistZ == 0.) FocDistZ = 0.5*radSag/cosAng; //focal length in vertical plane
			FocDistX = 0.5*radTan*cosAng; //focal length in vertical plane
			FocDistZ = 0.5*radSag/cosAng; //focal length in vertical plane
		}
	}
}

//*************************************************************************

srTMirrorEllipsoid::srTMirrorEllipsoid(const SRWLOptMirEl& srwlMirEl) : srTMirror(srwlMirEl.baseMir)
{
	m_p = srwlMirEl.p;
	m_q = srwlMirEl.q;
	m_angGraz = srwlMirEl.angGraz;
	m_radSag = srwlMirEl.radSag;

	//Validate parameters: make sure all are positive
	if((m_p <= 0) || (m_q <= 0) || (m_angGraz <= 0) || (m_radSag <= 0))
	{ ErrorCode = IMPROPER_OPTICAL_COMPONENT_ELLIPSOID; return;} //throw here?

	//Determine ellipsoid parameters in Local frame
	DetermineEllipsoidParamsInLocFrame(); 

	//Estimate focal lengths:
	double pq = m_p*m_q;
	double radTan = sqrt(pq*pq*pq)/(m_ax*m_az);
	EstimateFocalLengths(radTan, m_radSag);
}

//*************************************************************************

srTMirrorSphere::srTMirrorSphere(const SRWLOptMirSph& srwlMirSph) : srTMirror(srwlMirSph.baseMir)
{
	m_rad = srwlMirSph.rad;

	//Validate parameters: make sure all are positive
	if(m_rad == 0) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_MIRROR_SPHERE; return;} //throw here?

/*
	//Determine ellipsoid parameters in Local frame
	DetermineEllipsoidParamsInLocFrame(); 
*/

	//Estimate focal lengths:
	EstimateFocalLengths(m_rad, m_rad);
}

//*************************************************************************

srTMirrorToroid::srTMirrorToroid(srTStringVect* pMirInf, srTDataMD* pExtraData) : srTMirror(pMirInf, pExtraData)
{
	if((pMirInf == 0) || (pMirInf->size() < 5)) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}

	m_Rt = atof((*pMirInf)[2]);
	m_Rs = atof((*pMirInf)[3]);

	FocDistX = atof((*pMirInf)[8]);
	FocDistZ = atof((*pMirInf)[9]);
	if((FocDistX != 0.) && (FocDistZ != 0.)) return;

	//Estimating focal lengths (approximately!):
	EstimateFocalLengths(m_Rt, m_Rs);
}

//*************************************************************************

srTMirrorToroid::srTMirrorToroid(const SRWLOptMirTor& mirTor) : srTMirror(mirTor.baseMir)
{
	m_Rt = mirTor.radTan;
	m_Rs = mirTor.radSag;

	//Estimating focal lengths (approximately!):
	EstimateFocalLengths(m_Rt, m_Rs);
}

//*************************************************************************
//OBSOLETE?
srTThickMirrorGen::srTThickMirrorGen(srTStringVect* pElemInfo, srTDataMD* pExtraData) 
{
	if(pExtraData != 0) m_surfData = *pExtraData;

	char BufStr[256];

	TransvCenPoint.x = 0;
	strcpy(BufStr, (*pElemInfo)[4]); //$name[4]=num2str(xc)
	double aux_xc = atof(BufStr);
	if(::fabs(aux_xc) < 1.e+10) TransvCenPoint.x = aux_xc;

	TransvCenPoint.y = 0;
	strcpy(BufStr, (*pElemInfo)[5]); //$name[5]=num2str(yc)
	double aux_yc = atof(BufStr);
	if(::fabs(aux_yc) < 1.e+10) TransvCenPoint.y = aux_yc;

	m_apertShape = 1; //1- rectangular, 2- elliptical 
	strcpy(BufStr, (*pElemInfo)[6]); //$name[6]=num2str(apertShape)
	int iShape = atoi(BufStr);
	if((iShape > 0) && (iShape < 3)) m_apertShape = (char)iShape; //keep updated!

	m_ampReflectPerp = 1.;
	strcpy(BufStr, (*pElemInfo)[10]); //$name[10]=num2str(ampRefPerp)
	double aux_ampReflectPerp = atof(BufStr);
	if((aux_ampReflectPerp > 0) && (aux_ampReflectPerp < 1.)) m_ampReflectPerp = aux_ampReflectPerp;

	m_phaseShiftPerp = 0.;
	strcpy(BufStr, (*pElemInfo)[11]); //$name[11]=num2str(phShiftPerp)
	double aux_phaseShiftPerp = atof(BufStr);
	if(::fabs(aux_phaseShiftPerp) < 2.*3.141593) m_phaseShiftPerp = aux_phaseShiftPerp;

	m_ampReflectPar = 1.;
	strcpy(BufStr, (*pElemInfo)[12]); //$name[12]=num2str(ampRefPar)
	double aux_ampReflectPar = atof(BufStr);
	if((aux_ampReflectPar > 0) && (aux_ampReflectPar < 1.)) m_ampReflectPar = aux_ampReflectPar;

	m_phaseShiftPar = 0.;
	strcpy(BufStr, (*pElemInfo)[13]); //$name[13]=num2str(phShiftPar)
	double aux_phaseShiftPar = atof(BufStr);
	if(::fabs(aux_phaseShiftPar) < 2.*3.141593) m_phaseShiftPar = aux_phaseShiftPar;

	char m_axRot1 = 0;
	strcpy(BufStr, (*pElemInfo)[14]); //$name[14]=num2str(axRot1 - 1) //"0" means no rotation, "1" means vs "horizontal" axis, ...
	int aux_Rot = atoi(BufStr);
	if((aux_Rot > 0) && (aux_Rot < 4)) m_axRot1 = (char)aux_Rot; 

	strcpy(BufStr, (*pElemInfo)[15]); //$name[15]=num2str(angRot1)
	m_angRot1 = atof(BufStr);

	char m_axRot2 = 0;
	strcpy(BufStr, (*pElemInfo)[16]); //$name[16]=num2str(axRot2 - 1) //"0" means no rotation, "1" means vs "horizontal" axis, ...
	aux_Rot = atoi(BufStr);
	if((aux_Rot > 0) && (aux_Rot < 4)) m_axRot2 = (char)aux_Rot; 

	strcpy(BufStr, (*pElemInfo)[17]); //$name[17]=num2str(angRot2)
	m_angRot2 = atof(BufStr);

	char m_axRot3 = 0;
	strcpy(BufStr, (*pElemInfo)[18]); //$name[18]=num2str(axRot3 - 1) //"0" means no rotation, "1" means vs "horizontal" axis, ...
	aux_Rot = atoi(BufStr);
	if((aux_Rot > 0) && (aux_Rot < 4)) m_axRot3 = (char)aux_Rot; 

	strcpy(BufStr, (*pElemInfo)[19]); //$name[19]=num2str(angRot3)
	m_angRot3 = atof(BufStr);

	SetupNativeTransformation();


	//$name[7]="0" // Setup was finished or not
	//8 - foc. dist. x
	//9 - foc. dist. z
	strcpy(BufStr, (*pElemInfo)[7]); // Setup was completed or not
	int SetupIsCompleted = atoi(BufStr);
	if(SetupIsCompleted) 
	{
		strcpy(BufStr, (*pElemInfo)[8]);
		FocDistX = atof(BufStr);
		if(FocDistX == 0.) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}

		strcpy(BufStr, (*pElemInfo)[9]);
		FocDistZ = atof(BufStr);
		if(FocDistZ == 0.) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}
	}
	else
	{//Complete setup

/**
		if(ErrorCode = EstimateFocalDistancesAndCheckSampling()) return;

		//"erase" part of existing strings
		char *aStr=0;
		int AuxInfStartInd = 7;
		for(int k=AuxInfStartInd; k<(int)(pElemInfo->size()); k++)
		{
			aStr = (*pElemInfo)[k];
			//if(aStr != 0) delete[] aStr;
			if(aStr != 0) *aStr = '\0';
		}
		pElemInfo->erase(pElemInfo->begin() + AuxInfStartInd, pElemInfo->end());

		aStr = (*pElemInfo)[AuxInfStartInd];
		sprintf(aStr, "1");

		aStr = (*pElemInfo)[AuxInfStartInd + 1];
		sprintf(aStr, "%g", FocDistX);

		aStr = (*pElemInfo)[AuxInfStartInd + 2];
		sprintf(aStr, "%g", FocDistZ);

		aStr = (*pElemInfo)[AuxInfStartInd + 3];
		sprintf(aStr, "%d", OptPathOrPhase);
**/
	}

}

//*************************************************************************
