/************************************************************************//**
 * File: srradind.h
 * Description: Auxiliary SR calculation related structures (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRRADIND_H
#define __SRRADIND_H

//#ifdef __IGOR_PRO__
//#ifndef __SRIGINTR_H
//#include "srigintr.h"
//#endif
//#else
//#ifndef __SRIGORRE_H
//#include "srigorre.h"
//#endif
//#endif

#ifdef __IGOR_PRO__
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#else
#include "srigorre.h"
#endif

//-------------------------------------------------------------------------

struct srTSRWRadInData { 
// Nealy same data as in srTSRWRadStructAccessData
// to be used at the interface only !

	//bool BaseRadWasEmulated;
	float *pBaseRadX, *pBaseRadZ;
	waveHndl wRad, wRadX, wRadZ;
	int hStateRadX, hStateRadZ;
	double eStep, eStart, xStep, xStart, zStep, zStart;
	long ne, nx, nz;

	//double xStartTr, zStartTr;
	//bool UseStartTrToShiftAtChangingRepresToCoord;

	double RobsX, RobsZ;
	double RobsXAbsErr, RobsZAbsErr;
	double xc, zc;
	double xWfrMin, xWfrMax, zWfrMin, zWfrMax; // Exact borders of the Wavefront
	//char WfrEdgeCorrShouldBeDone; // To switch off/on manually

	double UnderSamplingX, UnderSamplingZ;
	char AllowAutoSwitchToPropInUnderSamplingMode;
	double InvUnderSamplingThreshold;

	//bool ResAfterWasEmulated;
	//bool DoNotResizeAfter;
	//srTRadResize* pResAfter;

	char Pres; // 0- Coord, 1- Ang.
	char PresT; // 0- Frequency (Photon Energy), 1- Time
	char LengthUnit; // 0- m; 1- mm; 
	char PhotEnergyUnit; // 0- eV; 1- keV; 
	char ElecFldUnit; // 0- Arb. Units, 1- sqrt(Phot/s/0.1%bw/mm^2)
	double avgPhotEn; //averarage photon energy for time-domain simulations

	//bool WfrQuadTermCanBeTreatedAtResizeX; // is used at the time of one resize only
	//bool WfrQuadTermCanBeTreatedAtResizeZ;

	//char ElectronBeamEmulated; // 0 by def.
	DOUBLE *pElecBeam;
	waveHndl wElecBeam; // Can be this or Trajectory
	int hStateElecBeam;

	waveHndl wTrj; // Can be this or Electron Beam
	int hStateTrj;

	//bool PropMatrWasEmulated;
	DOUBLE *p4x4PropMatr;
	waveHndl w4x4PropMatr;
	int hState4x4PropMatr;

	//bool MomWereEmulated;
	//float *pMomX, *pMomZ;
	DOUBLE *pMomX, *pMomZ; //OC130311
	waveHndl wMomX, wMomZ;
	int hStateMomX, hStateMomZ;

	//bool WfrAuxDataWasEmulated;
	DOUBLE *pWfrAuxData;
	waveHndl wWfrAuxData;
	int hStateWfrAuxData;

	//content of srTSRWRadStructWaveNames
	char NameRad[MAX_OBJ_NAME+1];
	char NameRadX[MAX_OBJ_NAME+1], NameRadZ[MAX_OBJ_NAME+1];
	char NameElecBeam[MAX_OBJ_NAME+1], NameTrj[MAX_OBJ_NAME+1];
	char Name4x4PropMatr[MAX_OBJ_NAME+1];
	char NameMomX[MAX_OBJ_NAME+1], NameMomZ[MAX_OBJ_NAME+1];
	char NameWfrAuxData[MAX_OBJ_NAME+1];

	//content of srTSRWRadStructWaveKeys
	char wRad_;
	char wRadX_, wRadZ_;
	char wElecBeam_, wTrj_;
	char w4x4PropMatr_;
	char wMomX_, wMomZ_;
	char wWfrAuxData_;

	srTSRWRadInData() 
	{
		wRad = NIL;
		pBaseRadX = 0; pBaseRadZ = 0; // This is checked in FinishWorkingWithSRWRadStruct !!!
		wRadX = wRadZ = NIL;

		pElecBeam = 0; wElecBeam = NIL;
		wTrj = NIL;

		p4x4PropMatr = 0; w4x4PropMatr = NIL;

		pMomX = pMomZ = 0;
		wMomX = wMomZ = NIL;

		pWfrAuxData = 0; wWfrAuxData = NIL;

		UnderSamplingX = UnderSamplingZ = 1.;
		AllowAutoSwitchToPropInUnderSamplingMode = 0;

		ElecFldUnit = 1;

		*NameRad = '\0';
		*NameRadX = '\0'; *NameRadZ = '\0';
		*NameElecBeam = '\0'; *NameTrj = '\0';
		*Name4x4PropMatr = '\0';
		*NameMomX = '\0'; *NameMomZ = '\0';
		*NameWfrAuxData = '\0';
	}
	//~srTSRWRadInData() {}
};

//-------------------------------------------------------------------------

struct srTSRWStokesInData {
// Nearly same as srTPowDensStructAccessData
// to be used at the interface only !

	float *pBaseSto;
	waveHndl wSto;
	int hStateSto;
	//bool MemoryWasAllocatedInternally;

	double eStep, eStart, xStep, xStart, zStep, zStart, yStep, yStart;
	long ne, nx, nz, ny;

	srTSRWStokesInData()
	{
		pBaseSto = 0;
        wSto = NIL;
	}
};

//-------------------------------------------------------------------------

struct srTSRWPowDensInData {
// Nearly same as srTPowDensStructAccessData
// to be used at the interface only !
	float *pBasePowDens;
	waveHndl wPowDens;
	int hStatePowDens;

	double xStep, xStart, zStep, zStart;
	long nx, nz;
};

//-------------------------------------------------------------------------

struct srTIgorWaveAccessData {
	char* pWaveData;
	char WaveType[2]; // 'f'|'d'|'cf'|'cd'
	long AmOfDims;
	long DimSizes[10];
	double DimStartValues[10];
	double DimSteps[10];
	char DimUnits[10][255];
	char DataUnits[255];

	waveHndl wHndl;
	int hState;

	char NameOfWave[50];

	srTIgorWaveAccessData()
	{
		pWaveData = 0;
		*WaveType = '\0';
		AmOfDims = -1;
		for(int i=0; i<10; i++) 
		{
			DimSizes[i] = -1;
			DimStartValues[i] = 1.E+23;
			DimSteps[i] = 1.E+23;
			*(DimUnits[i]) = '\0';
		}
		*NameOfWave = '\0';
		*DataUnits = '\0';
	}
};

//-------------------------------------------------------------------------

#endif
