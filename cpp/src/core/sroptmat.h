/************************************************************************//**
 * File: sroptmat.h
 * Description: Some material data (to be replaced by a more sophisticated library)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTMAT_H
#define __SROPTMAT_H

//#ifndef __SRSEND_H
//#include "srsend.h"
//#endif
#include "srstraux.h"

//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

class srTOptMater {
	static double AtomicWeightArray[];

	static int LenAttenArgArray1;
	static double AttenArgArray1[];
	static double AttenLengthArrayBe[];
	static double AttenLengthArrayPyroC[];

	int AmOfSpecMaterials;

public:

	srTOptMater()
	{
		AmOfSpecMaterials = 2; // Change when more spec. mater. added
	}

	void CharactOfSpecMater(int SpecMatNo, int& Z, double& Density);

	double AtomicWeight(int Z)
	{
		if((Z < 0) || (Z > 103)) return 0.;
		return AtomicWeightArray[Z-1];
	}
	int CompRefractDelta(int SpecMatNo, int InZ, double InDensity, double PhotEn, double& Delta); // PhotEn in eV, Density in g/cm^3
	double AttenLength(int SpecMatNo, double PhotEn);
	double InterpolFunction(double* ArgArray, double* AttenLenArray, int LenArray, double Arg);
};

//*************************************************************************

#endif
