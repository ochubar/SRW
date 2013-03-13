/************************************************************************//**
 * File: srmagelem.h
 * Description: Magnetic element
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMAGELEM_H
#define __SRMAGELEM_H

#include "srobject.h"
#include "objcont.h"
#include "srercode.h"
#include "gmvect.h"

//*************************************************************************

class srTGenTrjDat;
class srTEbmDat;
class srTWfrSmp;
class srTStokesStructAccessData;
struct srTWigComSASE;
//class srTTrjDat;
//class srTTrjDat3d;
//class srTMagFieldPeriodic;
//class srTMagFldCont;

//*************************************************************************

class srTMagElem : public CGenObject {
protected:
	TVector3d mCenP; //SRWLIB

public:
	int ErrorCode;
	double gsStart;
	double gsEnd;

	srTMagElem(char* Name) : CGenObject(Name) {}
	srTMagElem(const TVector3d& cenP)
	{
		mCenP = cenP;
		ErrorCode = 0;
		gsStart = gsEnd = 0.;
	}
	srTMagElem()
	{
		ErrorCode = 0;
		gsStart = gsEnd = 0.;
	}

	virtual double GetLongExtent() { return 0.;}
	//virtual void SetupTrjDat(srTTrjDat*) { throw SR_COMP_NOT_IMPLEMENTED_FOR_GIVEN_MAG_FLD;}

	virtual srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat*) { return 0;}
	virtual void SetupWigSASE(srTWigComSASE&) {} //sets up SASE wiggler for Genesis
	virtual void ComputeSR_Stokes(srTEbmDat* pElecBeam, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes) { throw SR_COMP_NOT_IMPLEMENTED_FOR_GIVEN_MAG_FLD;}
	virtual void ComputeParticlePropagMatrix(double s, TMatrix2d& Mx, TMatrix2d& Mz) {}
	virtual void compB(TVector3d& inP, TVector3d& outB) {}

	static int FindMagElemWithSmallestLongPos(CObjCont<CGenObject>& AuxCont);
	void GetMagnFieldLongLim(double& sSt, double& sEn) { sSt = gsStart; sEn = gsEnd;} //SRWLIB
};

//*************************************************************************

#endif
