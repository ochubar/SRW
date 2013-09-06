/************************************************************************//**
 * File: sroptang.h
 * Description: Optical element: Angle (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2013
 *
 * Copyright (C) Brookhaven National Laboratory, Upton, NY, USA
 * Copyright (C) European X-ray Free Electron Laser Facility GmbH, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTANG_H
#define __SROPTANG_H

#include "sroptelm.h"

//*************************************************************************

class srTOptAngle : public srTGenOptElem {

public:
	double AngX, AngY;

	srTOptAngle(double InAngX=0, double InAngY=0)
	{
		AngX = InAngX;
		AngY = InAngY;
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect) //virtual
	{
		//return PropagateRadiationMeth_0(pRadAccessData);
		int res = 0;
		if(res = PropagateRadiationMeth_0(pRadAccessData)) return res;

		if(res = PropagateWaveFrontRadius(pRadAccessData)) return res; //because this is the same for all E slices!
		if(res = PropagateRadMoments(pRadAccessData, 0)) return res;
		return 0;
	}

	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		if(result = TraverseRadZXE(pRadAccessData)) return result;
		//consider programming Angle on angular side by simple change of limits
		//however note potential problems for many photon energies!

		//if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		//if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		return 0;
	}

	void RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
	{// e in eV; Length in m !!!
	 // Operates on Coordinates side !!!
		double twoPi_d_Lambda = (5.06773065e+06)*EXZ.e;
		double phaseShift = twoPi_d_Lambda*(AngX*EXZ.x + AngY*EXZ.z); //to check

		double cosPh = cos(phaseShift), sinPh = sin(phaseShift);
		double NewExRe = (*(EPtrs.pExRe))*cosPh - (*(EPtrs.pExIm))*sinPh;
		double NewExIm = (*(EPtrs.pExRe))*sinPh + (*(EPtrs.pExIm))*cosPh;
		*(EPtrs.pExRe) = (float)NewExRe; *(EPtrs.pExIm) = (float)NewExIm; 
		double NewEzRe = (*(EPtrs.pEzRe))*cosPh - (*(EPtrs.pEzIm))*sinPh;
		double NewEzIm = (*(EPtrs.pEzRe))*sinPh + (*(EPtrs.pEzIm))*cosPh;
		*(EPtrs.pEzRe) = (float)NewEzRe; *(EPtrs.pEzIm) = (float)NewEzIm; 
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRad, srTMomentsRatios*) 
	{
		int Offset = 0;
		for(long ie=0; ie<pRad->ne; ie++)
		{
			srTMomentsPtrs MomX(pRad->pMomX + Offset);
			srTMomentsPtrs MomZ(pRad->pMomZ + Offset);

			*(MomX.pXP) += AngX; *(MomX.pZP) += AngY; 
			*(MomZ.pXP) += AngX; *(MomZ.pZP) += AngY; 
			Offset += 11;
		}
		return 0;
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{// Angle modifies xc, zc !?
		//pRadAccessData->xc -= AngX*(pRadAccessData->RobsX);
		//pRadAccessData->zc -= AngY*(pRadAccessData->RobsZ);
		return 0;
	}

	//int PropagateRadiationSimple(srTSRWRadStructAccessData* pRadAccessData)
	//{
	//	int result;
	//	if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
	//	if(result = TraverseRadZXE(pRadAccessData)) return result;
	//	return 0;
	//	//consider programming Angle on angular side by simple change of limits
	//	//however note potential problems for many photon energies!
	//}

	//int PropagateRadiationSimple_CoordRepres(srTSRWRadStructAccessData* pRadAccessData)
	//{
	//	int result;
	//	if(result = TraverseRadZXE(pRadAccessData)) return result;
	//	return 0;
	//}

	//consider programming Angle on angular side by simple change of limits
	//however note potential problems for many photon energies!
	//int PropagateRadiationSimple_AngRepres(srTSRWRadStructAccessData* pRad)
	//{
	//}
};

//*************************************************************************

class srTOptShift : public srTGenOptElem {

public:
	double ShiftX, ShiftY;

	srTOptShift(double InShiftX=0, double InShiftY=0)
	{
		ShiftX = InShiftX;
		ShiftY = InShiftY;
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResBeforeAndAfterVect) //virtual
	{
		//return PropagateRadiationMeth_0(pRadAccessData);
		int res = 0;
		if(res = PropagateRadiationMeth_0(pRadAccessData)) return res;

		if(res = PropagateWaveFrontRadius(pRadAccessData)) return res; //because this is the same for all E slices!
		if(res = PropagateRadMoments(pRadAccessData, 0)) return res;
		return 0;
	}

	int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadAccessData)
	{
		int result;
		if(pRadAccessData->Pres != 0) if(result = SetRadRepres(pRadAccessData, 0)) return result;
		//consider programming Shift on angular side by adding "adding" simusoidal component

		pRadAccessData->ShiftWfrByInterpolVsXZ(ShiftX, ShiftY);

		//if(result = PropagateWaveFrontRadius(pRadAccessData)) return result;
		//if(result = PropagateRadMoments(pRadAccessData, 0)) return result;
		return 0;
	}

	int PropagateRadMoments(srTSRWRadStructAccessData* pRad, srTMomentsRatios*) 
	{
		int Offset = 0;
		for(long ie=0; ie<pRad->ne; ie++)
		{
			srTMomentsPtrs MomX(pRad->pMomX + Offset);
			srTMomentsPtrs MomZ(pRad->pMomZ + Offset);

			*(MomX.pX) += ShiftX; *(MomX.pZ) += ShiftY; 
			*(MomZ.pX) += ShiftX; *(MomZ.pZ) += ShiftY; 
			Offset += 11;
		}
		return 0;
	}

	int PropagateWaveFrontRadius(srTSRWRadStructAccessData* pRadAccessData)
	{// Shift modifies xc, zc !?
		//pRadAccessData->xc += ShiftX;
		//pRadAccessData->zc += ShiftY;
		return 0;
	}
};

//*************************************************************************

#endif
