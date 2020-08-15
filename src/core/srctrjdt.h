/************************************************************************//**
 * File: srctrjdt.h
 * Description: Electron trajectory calculation, constant magnetic field case (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRCTRJDT_H
#define __SRCTRJDT_H

#include "srgtrjdt.h"
#include "srmagfld.h"

//*************************************************************************

class srTMagFieldConstant : public srTMagElem {
public:
	double Bx, Bz;

	srTMagFieldConstant(double InBx, double InBz) 
	{
		Bx = InBx; Bz = InBz;
	}
	srTMagFieldConstant() { Bx = Bz = 0.;}

	void AnalyzeFieldSymmetry(char& FieldIsSymOverX, char& FieldIsSymOverZ)
	{
		char zFieldPresent = 0, xFieldPresent = 0;
		if(Bx != 0.) xFieldPresent = 1;
		if(Bz != 0.) zFieldPresent = 1;

		FieldIsSymOverX = FieldIsSymOverZ = 0;
		if(!xFieldPresent) FieldIsSymOverZ = 1;
		if(!zFieldPresent) FieldIsSymOverX = 1;
	}

	//void SetupTrjDat(srTTrjDat*) //virtual
	//{
	//}
    srTGenTrjDat* CreateAndSetupNewTrjDat(srTEbmDat*); //virtual
};

//*************************************************************************

class srTConstTrjDat : public srTGenTrjDat {
public:

	srTMagFieldConstant MagConst;

	srTConstTrjDat() {}

	void CheckIfHorOrVertFieldIsZero()
	{
		HorFieldIsNotZero = (MagConst.Bx != 0.);
		VerFieldIsNotZero = (MagConst.Bz != 0.);
	}
	int MagFieldPeriodicity() { return 1;} // Non-periodic
	char MagFieldIsConstant() { return 1;}

	int ConvertToArbTrjDat(char Periodicity, srTWfrSmp&, srTGenTrjHndl&);
	void DetermineIntegLimitsForArbTrj(srTWfrSmp& DistrInfoDat, double& sStart, double& sEnd);

	int ShowLimitsAndInitInteg(srTWfrSmp& DistrInfoDat, char LongIntType, double& sIntegStart, double& sIntegFin, int& AmOfPer, bool doInit = true) 
	{
		DetermineIntegLimitsForArbTrj(DistrInfoDat, sIntegStart, sIntegFin);
		AmOfPer = 1;
		// Make any initialization, if needed
		return 0;
	}
};

//*************************************************************************

#endif
