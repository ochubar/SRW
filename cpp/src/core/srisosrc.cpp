/************************************************************************//**
 * File: srisosrc.cpp
 * Description: Isotropic source
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srisosrc.h"

//*************************************************************************

int srTIsotrSrc::CheckInputConsistency()
{
// Fill-in
	return 0;
}

//*************************************************************************

void srTIsotrSrc::SetupSourceConstants()
{
	LongDist = DistrInfoDat.yStart - EbmDat.s0;

// To give Phot/s/0.1%bw/mm^2
// assuming EbmDat.Current in A, PhotPerBW in Phot/0.1%bw, LongDist in m
	double NormConstSpecFluxPerUnSurf = EbmDat.Current*PhotPerBW/(LongDist*LongDist*((1.602189246E-13)*4.*3.1415926535898));
	NormConstElField = sqrt(NormConstSpecFluxPerUnSurf);

	x0 = EbmDat.x0; z0 = EbmDat.z0; 
	SigX = sqrt(EbmDat.Mxx); SigZ = sqrt(EbmDat.Mzz);
}

//*************************************************************************

int srTIsotrSrc::CreateWavefrontElField(srTSRWRadStructAccessData& RadAccessData)
{// m, eV !
	int result;
	const double TwoPI = 6.28318530717959;
	const double InvTwoPI = 1./TwoPI;

	if(result = CheckInputConsistency()) return result;

	SetupSourceConstants();

	double LongDistE2 = LongDist*LongDist;
	double LongDistE3 = LongDistE2*LongDist;

	double Inv_ye2 = 1./LongDistE2;
	double enMult = (2.53384080189E+06)*LongDist;

	float* tRadX = RadAccessData.pBaseRadX;
	float* tRadZ = RadAccessData.pBaseRadZ;

	double z = RadAccessData.zStart - z0;
	for(int iz=0; iz<RadAccessData.nz; iz++)
	{
		double ze2 = z*z;

		double x = RadAccessData.xStart - x0;
		for(int ix=0; ix<RadAccessData.nx; ix++)
		{
			double xe2 = x*x;
			double a = (xe2 + ze2)*Inv_ye2;
			double CoordPhaseTerm = enMult*a*(1 - 0.25*a + 0.125*a*a);

			double en = RadAccessData.eStart;
			for(int ie=0; ie<RadAccessData.ne; ie++)
			{
				double Phase = en*CoordPhaseTerm;
				//Phase -= TwoPI*long(Phase*InvTwoPI);
				Phase -= TwoPI*((long long)(Phase*InvTwoPI));

				double CosPh = cos(Phase), SinPh = sin(Phase);

				double R2 = LongDistE2 + xe2 + ze2;
				double GeomFact = LongDistE3/(R2*sqrt(R2));
				double NormConst = NormConstElField*GeomFact;

				double ReA = NormConst*CosPh, ImA = NormConst*SinPh;
				SetupProperPolariz(ReA, ImA, tRadX, tRadZ);

				tRadX += 2; tRadZ += 2;
				en += RadAccessData.eStep;
			}
			x += RadAccessData.xStep;
		}
		z += RadAccessData.zStep;
	}
	return 0;
}

//*************************************************************************
