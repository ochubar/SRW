/************************************************************************//**
 * File: sroptfoc.cpp
 * Description: Optical element: Focusing
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptfoc.h"

//*************************************************************************

int srTFocusingElem::TuneRadForPropMeth_1(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& PostResize)
{
  	srTMomentsRatios* MomRatArray = new srTMomentsRatios[pRadAccessData->ne];
	if(MomRatArray == 0) return MEMORY_ALLOCATION_FAILURE;

	int result;
	if(pRadAccessData->Pres != 0) // Go to spatial...
		if(result = SetRadRepres(pRadAccessData, 0)) return result;

	if(result = PropagateRadMoments(pRadAccessData, MomRatArray)) return result;
	
	srTMomentsRatios* tMomRatArray = MomRatArray;

	//float pxMaxMomX = (float)(1.e-23), pxMinMomX = (float)(1.e+23), pzMaxMomX = (float)(1.e-23), pzMinMomX = (float)(1.e+23);
	//float pxMaxMomZ = (float)(1.e-23), pxMinMomZ = (float)(1.e+23), pzMaxMomZ = (float)(1.e-23), pzMinMomZ = (float)(1.e+23);
	double pxMaxMomX = 1.e-23, pxMinMomX = 1.e+23, pzMaxMomX = 1.e-23, pzMinMomX = 1.e+23; //OC130311
	double pxMaxMomZ = 1.e-23, pxMinMomZ = 1.e+23, pzMaxMomZ = 1.e-23, pzMinMomZ = 1.e+23;

	for(long ie=0; ie<pRadAccessData->ne; ie++)
	{
		if(tMomRatArray->RxpxpMomX > pxMaxMomX) pxMaxMomX = tMomRatArray->RxpxpMomX;
		if(tMomRatArray->RxpxpMomZ > pxMaxMomZ) pxMaxMomZ = tMomRatArray->RxpxpMomZ;
		if(tMomRatArray->RxpxpMomX < pxMinMomX) pxMinMomX = tMomRatArray->RxpxpMomX;
		if(tMomRatArray->RxpxpMomZ < pxMinMomZ) pxMinMomZ = tMomRatArray->RxpxpMomZ;

		if(tMomRatArray->RzpzpMomX > pzMaxMomX) pzMaxMomX = tMomRatArray->RzpzpMomX;
		if(tMomRatArray->RzpzpMomZ > pzMaxMomZ) pzMaxMomZ = tMomRatArray->RzpzpMomZ;
		if(tMomRatArray->RzpzpMomX < pzMinMomX) pzMinMomX = tMomRatArray->RzpzpMomX;
		if(tMomRatArray->RzpzpMomZ < pzMinMomZ) pzMinMomZ = tMomRatArray->RzpzpMomZ;
	
		tMomRatArray++;
	}
	//float pxMax = (pxMaxMomX > pxMaxMomZ)? pxMaxMomX : pxMaxMomZ;
	//float pzMax = (pzMaxMomX > pzMaxMomZ)? pzMaxMomX : pzMaxMomZ;
	double pxMax = (pxMaxMomX > pxMaxMomZ)? pxMaxMomX : pxMaxMomZ; //OC130311
	double pzMax = (pzMaxMomX > pzMaxMomZ)? pzMaxMomX : pzMaxMomZ;

// Only the following is used for tuning
// This is very unsafe !!! Needs Reparation !!!

	//float pxdTot = (float)::fabs((FocDistX - pRadAccessData->RobsX)/FocDistX);
	//float pzdTot = (float)::fabs((FocDistZ - pRadAccessData->RobsZ)/FocDistZ);
	double pxdTot = ::fabs((FocDistX - pRadAccessData->RobsX)/FocDistX); //OC130311
	double pzdTot = ::fabs((FocDistZ - pRadAccessData->RobsZ)/FocDistZ);

	const double pdMinimumAllowed = 0.03;
	if(pxdTot < pdMinimumAllowed) pxdTot = pdMinimumAllowed;
	if(pzdTot < pdMinimumAllowed) pzdTot = pdMinimumAllowed;

	srTRadResize RadResize;
	RadResize.pxm = RadResize.pxd = RadResize.pzm = RadResize.pzd = 1.;

	const double ResizeTol = 0.15;
	char xResizeNeeded = (pxdTot - 1. > ResizeTol);
	char zResizeNeeded = (pzdTot - 1. > ResizeTol);
	if(xResizeNeeded) RadResize.pxd = pxdTot;
	if(zResizeNeeded) RadResize.pzd = pzdTot;
	if(xResizeNeeded || zResizeNeeded) if(result = RadResizeGen(*pRadAccessData, RadResize)) return result;
	
	PostResize.pxm = PostResize.pzm = PostResize.pxd = PostResize.pzd = 1.;
	char xPostResizeNeeded = (1.- ResizeTol > pxdTot);
	char zPostResizeNeeded = (1.- ResizeTol > pzdTot);

	if(xPostResizeNeeded) 
	{
		PostResize.pxd = pxMax;
	}
	if(zPostResizeNeeded) 
	{
		PostResize.pzd = pzMax;
	}

	if(MomRatArray != 0) delete[] MomRatArray;
	return 0;
}

//*************************************************************************
