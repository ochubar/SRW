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
//#include "gmvect.h"
#include "gmtrans.h" //OC160615

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
	//TVector3d mCenP; //SRWLIB
	//TVector3d mAxV; //SRWLIB (axis vector) //OC150615
	//double mAng; //SRWLIB (rotation angle about axis)
	gmTrans mTrans; //OC160615

public:
	int ErrorCode;
	double gsStart;
	double gsEnd;

	srTMagElem(char* Name) : CGenObject(Name) {}
	//srTMagElem(const TVector3d& cenP)
	srTMagElem(const TVector3d& cenP, const TVector3d& axV, double ang=0)
	{
		//mCenP = cenP;
		//mAxV = axV; //OC160615
		//mAng = ang;

		SetupOrient(cenP, axV, ang);
		ErrorCode = 0;
		gsStart = gsEnd = 0.;
	}
	srTMagElem()
	{
		mTrans.SetupIdent();

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
	void GetMagnFieldLongLim(double& sSt, double& sEn) //SRWLIB
	{
		//sSt = gsStart; sEn = gsEnd;
		TVector3d Pstart(0., 0., gsStart), Pend(0., 0., gsEnd); //OC170615
		Pstart = mTrans.TrPoint(Pstart); Pend = mTrans.TrPoint(Pend);
		sSt = Pstart.z; sEn = Pend.z;
		if(sSt > sEn)
		{
			double aux = sSt;
			sSt = sEn; sEn = aux;
		}
	}

protected:
	void SetupOrient(const TVector3d& cenP, const TVector3d& axV, double ang=0)
	{
		mTrans.SetupIdent();
		TVector3d vEz(0., 0., 1.), vZero(0., 0., 0.);
		if(ang != 0.)
		{
			mTrans.SetupRotation(vZero, vEz, ang);
		}
		TVector3d inAxV = axV;
		if(!(inAxV.isZero() || ((axV.x == 0.) && (axV.y == 0.) && (axV.z == 1.))))
		{
			inAxV.Normalize();
			gmTrans auxTrans;
			auxTrans.SetupRotation(vZero, vEz, inAxV);
			mTrans.TrMult(auxTrans, 'L');
		}
		if(!((cenP.x == 0.) && (cenP.y == 0.) && (cenP.z == 0.))) 
		{
			gmTrans auxTrans;
			auxTrans.SetupTranslation(cenP);
			mTrans.TrMult(auxTrans, 'L');
		}
	}
};

//*************************************************************************

#endif
