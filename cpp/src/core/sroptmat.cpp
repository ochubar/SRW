/************************************************************************//**
 * File: sroptmat.cpp
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

#include "sroptmat.h"
#include "srerror.h"

//*************************************************************************

double srTOptMater::AtomicWeightArray[] = {
	1.00794, 4.00260, 6.941, 9.01218, 10.81, 12.011, 14.0067, 15.9994, 18.998403, 20.179,
	22.98977, 24.305, 26.98154, 28.0855, 30.97376, 32.06, 35.453, 39.948, 39.0983, 40.08,
	44.9559, 47.88, 50.9415, 51.996, 54.9380, 55.847, 58.9332, 58.69, 63.546, 65.38, 69.72,
	72.59, 74.9216, 78.96, 79.904, 83.80, 85.4678, 87.62, 88.9059, 91.22, 92.9064, 95.94,
	98, 101.07, 102.9055, 106.42, 107.8682, 112.41, 114.82, 118.69, 121.75, 127.60, 126.9045,
	131.29, 132.9054, 137.33, 138.9055, 140.12, 140.9077, 144.24, 145, 150.36, 151.96, 157.25,
	158.9254, 162.50, 164.9304, 167.26, 168.9342, 173.04, 174.967, 178.49, 180.9479, 183.85,
	186.207, 190.2, 192.22, 195.08, 196.9665, 200.59, 204.383, 207.2, 208.9804, 209, 210, 222,
	223, 226.0254, 227.0278, 232.0381, 231.0359, 238.0289, 237.0482, 244, 243, 247, 247, 251,
	252, 257, 258, 259, 260
// Continue starting from Z=43
};

//*************************************************************************

int srTOptMater::LenAttenArgArray1 = 48;
double srTOptMater::AttenArgArray1[] = {
	1000, 1100, 1210, 1330, 1460, 1610, 1770, 1950, 2140, 2360, 2590, 2850, 3140, 3450, 3800, 
	4180, 4590, 5050, 5560, 6120, 6730, 7400, 8140, 8950, 9850, 10830, 11920, 13110, 14420,
	15860, 17450, 19190, 21110, 23230, 25550, 28100, 30910, 34000, 37400, 41140, 45260,
	49790, 54760, 60240, 66260, 72890, 80180, 88200, 97020
};

double srTOptMater::AttenLengthArrayBe[] = {
	9.03E-03, 1.19E-02, 1.58E-02, 2.11E-02, 2.81E-02, 3.77E-02, 5.06E-02, 6.82E-02, 9.22E-02,
	1.25E-01, 1.70E-01, 2.31E-01, 3.16E-01, 4.32E-01, 5.90E-01, 8.07E-01, 1.10E+00, 1.51E+00,
	2.06E+00, 2.80E+00, 3.82E+00, 5.21E+00, 7.11E+00, 9.69E+00, 1.32E+01, 1.80E+01, 2.44E+01,
	3.32E+01, 4.49E+01, 6.06E+01, 0.81E+02, 1.07E+02, 1.40E+02, 1.78E+02, 2.21E+02, 2.66E+02,
	3.09E+02, 3.45E+02, 3.72E+02, 3.88E+02, 3.94E+02, 3.92E+02, 3.84E+02, 3.71E+02, 3.57E+02,
	3.46E+02, 3.32E+02, 3.17E+02, 3.04E+02
};

double srTOptMater::AttenLengthArrayPyroC[] = {
	2.88E-03, 3.77E-03, 4.94E-03, 6.49E-03, 8.54E-03, 1.13E-02, 1.49E-02, 1.97E-02, 2.60E-02,
	3.45E-02, 4.59E-02, 6.10E-02, 8.13E-02, 1.09E-01, 1.45E-01, 1.95E-01, 2.61E-01, 3.52E-01,
	4.74E-01, 6.40E-01, 8.65E-01, 1.17E+00, 1.59E+00, 2.17E+00, 2.96E+00, 4.04E+00, 5.53E+00,
	7.58E+00, 1.04E+01, 1.43E+01, 1.95E+01, 2.67E+01, 3.61E+01, 4.91E+01, 0.66E+02, 0.88E+02,
	1.14E+02, 1.46E+02, 1.80E+02, 2.16E+02, 2.50E+02, 2.80E+02, 3.02E+02, 3.18E+02, 3.27E+02,
	3.22E+02, 3.18E+02, 3.12E+02, 3.03E+02
};

//*************************************************************************

void srTOptMater::CharactOfSpecMater(int SpecMatNo, int& Z, double& Density)
{
	if((SpecMatNo < 1) || (SpecMatNo > AmOfSpecMaterials)) return;
	switch (SpecMatNo) 
	{
	case 1:
		Z = 4; Density = 1.845;
		break; // Be
	case 2:
		Z = 6; Density = 2.20;
		break; // C
	}
}

//*************************************************************************

int srTOptMater::CompRefractDelta(int SpecMatNo, int InZ, double InDensity, double PhotEn, double& Delta) // PhotEn in eV, Density in g/cm^3
{
	const double re = 2.817938070E-13; // class. el. rad. in cm
	const double Na = 6.02204531E+23;
	//const double amu = 1.660565586E-24; // atomic mass unit in g
	const double TwoPI = 6.2831853071796;
	
	const double Con = re*Na/TwoPI;
	double WavelengthIn_cm = 1.239854E-04/PhotEn;
	
	int Z = InZ;
	double Density;
	if((SpecMatNo > 0) && (SpecMatNo <= AmOfSpecMaterials)) 
	{
		CharactOfSpecMater(SpecMatNo, Z, Density);
		if(InDensity > 0.) Density = InDensity;
	}
	else
	{
		if(InDensity <= 0.) return DENSITY_WALUE_NEEDED;
		if((Z < 0) || (Z > 103)) return BAD_MATER_Z_VALUE; // Modify later
		Density = InDensity;
	}
	
	double A = AtomicWeight(Z);

	Delta = Z*Con*WavelengthIn_cm*WavelengthIn_cm*Density/A;
	return 0;
}

//*************************************************************************

double srTOptMater::AttenLength(int SpecMatNo, double PhotEn)
{
	if((SpecMatNo < 1) || (SpecMatNo > AmOfSpecMaterials)) return -1.;
	if((PhotEn < 1000.) || (PhotEn > 97000.))
	{
		//srTSend Send; Send.AddWarningMessage(&gVectWarnNos, POOR_PREC_OF_ATTEN_LEN);
		CErrWarn::AddWarningMessage(&gVectWarnNos, POOR_PREC_OF_ATTEN_LEN);
	}

	double *pArgArray = 0, *pAttenLength;
	int *pLenArray;
	switch (SpecMatNo) 
	{
	case 1:
		pArgArray = AttenArgArray1; pAttenLength = AttenLengthArrayBe;
		pLenArray = &LenAttenArgArray1;
		break; // Be
	case 2:
		pArgArray = AttenArgArray1; pAttenLength = AttenLengthArrayPyroC;
		pLenArray = &LenAttenArgArray1;
		break; // C
	}

	return InterpolFunction(pArgArray, pAttenLength, *pLenArray, PhotEn);
}

//*************************************************************************

double srTOptMater::InterpolFunction(double* ArgArray, double* AttenLenArray, int LenArray, double Arg)
{
	double ArgStart = *ArgArray, ArgFin = *(ArgArray + LenArray - 1);
	if(Arg < ArgStart) return AttenLenArray[0];
	if(Arg > ArgFin) return AttenLenArray[LenArray - 1];

	double *tArgArray = ArgArray + 1;
	int i;
	for(i=1; i<LenArray; i++)
	{
		if(Arg < *(tArgArray++)) break;
	}
	int iSt = i - 1, iFi = i;
	double LinInterpol = AttenLenArray[iSt] + ((Arg - ArgArray[iSt])/(ArgArray[iFi] - ArgArray[iSt]))*(AttenLenArray[iFi] - AttenLenArray[iSt]);
// Interpolate better later
	return LinInterpol;
}

//*************************************************************************
