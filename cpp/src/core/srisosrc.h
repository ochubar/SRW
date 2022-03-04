/************************************************************************//**
 * File: srisosrc.h
 * Description: Isotropic source (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRISOSRC_H
#define __SRISOSRC_H

#include "srebmdat.h"
#include "srstraux.h"
#include "gmmeth.h"

struct SRWLStructPointSource;
typedef struct SRWLStructPointSource SRWLPtSrc;

//*************************************************************************

class srTIsotrSrc {

	double LongDist;
	double NormConstElField;

public:

	srTEbmDat EbmDat;

	char InputWasModified;

	double PhotPerBW;
	double Flux;
	int Polar, UnitFlux;
	double x0, z0, LongPos;
	double SigX, SigZ;

	srTWfrSmp DistrInfoDat;

	srTIsotrSrc() {} //required for compiling for Igor (?)
	srTIsotrSrc(SRWLPtSrc*);

	int CheckInputConsistency();
	void SetupSourceConstants();
	int CreateWavefrontElField(srTSRWRadStructAccessData&);
	void ComputeElectricField(srTSRWRadStructAccessData &wfr); //, double R0, double pow, int polar);

	void SetupProperPolariz(double ReA, double ImA, double x, double z, float* tRadX, float* tRadZ)
	{// Same for Isotropic Source and Gaussian Beam
		const double c = 0.70710678118655;
		switch(Polar) 
		{
		case 1: // Lin. Hor.
			*tRadX = (float)ReA; *(tRadX+1) = (float)ImA; *tRadZ = 0.; *(tRadZ+1) = 0.;
			break;
		case 2: // Lin. Vert.
			*tRadX = 0.; *(tRadX+1) = 0.; *tRadZ = (float)ReA; *(tRadZ+1) = (float)ImA;
			break;
		case 3: // Lin. 45
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(c*ReA); *(tRadZ+1) = (float)(c*ImA);
			break;
		case 4: // Lin. 135
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(-c*ReA); *(tRadZ+1) = (float)(-c*ImA);
			break;
		case 5: // Circ. Right
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(-c*ImA); *(tRadZ+1) = (float)(c*ReA);
			break;
		case 6: // Circ. Left
			*tRadX = (float)(c*ReA); *(tRadX+1) = (float)(c*ImA); *tRadZ = (float)(c*ImA); *(tRadZ+1) = (float)(-c*ReA);
			break;
		case 7: // Radial
			double ang = CGenMathMeth::azimAngFrXY(x, z);
			double cosAng = cos(ang), sinAng = sin(ang);
			*tRadX = (float)(cosAng*ReA); *(tRadX+1) = (float)(cosAng*ImA); *tRadZ = (float)(sinAng*ReA); *(tRadZ+1) = (float)(sinAng*ImA);
			break;
		}
	}
};

//*************************************************************************

#endif
