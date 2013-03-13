/************************************************************************//**
 * File: srigintr.h
 * Description: Interface to IGOR Pro (WaveMetrics) header
 * Project: Synchrotron Radiation Workshop (SRW)
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRIGINTR_H
#define __SRIGINTR_H

//#ifdef FIRST_XOP_ERR
//#undef FIRST_XOP_ERR
//#endif

#ifdef __IGOR_PRO__
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#else
#include "srigorre.h"
#endif
//NOTE: this structures can be defined and used without IGOR!

#ifndef __SRERCODE_H
#include "srercode.h"
#endif

//*************************************************************************

/* Prototypes */

//HOST_IMPORT void main(IORecHandle ioRecHandle);

//*************************************************************************

struct srTIgorRadInputStruct {
	waveHndl wRad;
// SRW Radiation structure (text wave)

//	waveHndl wOut2;
// Complex wave for output of Horizontal Fourier-component of the radiation Electric Field (units of sqrt(Ph/(s*mm^2*(.1%BW))))
//	waveHndl wOut1;
// Complex wave for output of Vertical Fourier-component of the radiation Electric Field (units of sqrt(Ph/(s*mm^2*(.1%BW))))
// The output data are the flattened loops, assuming:
// the Outmost (external) loop over Vertical Obs. Coord.
// the next (middle) loop over Horizontal Obs. Coord.
// the Inmost (internal) loop over Photon Energy
// Currently, the output waves are assumed to be created from outside;
// the size of the waves should be (nz*nx*nLambda)

	waveHndl wAuxParam;
// 0: Allow auto choise of Observation points for further Propagation (0- No; 1- Yes)
// 1: Observation oversampling parameter

	waveHndl wIntPar;
// Integration Parameters wave
// 0: Integration Method (0- Manual; 1- AutoUndulator; 2- AutoWiggler)
// 1: Relative Precision or Step [m]
// 2: Initial point for the Integration sStart [m]
// 3: Final point for the Integration sEnd [m]
// 4: Maximum number of points to save in memory
// 5: Flag: Try to apply Near-Field residual term is the far-field one fails (1)

	waveHndl wObservation;
// Observation Wave Format1:
// 0: Flag AtInfinity/InNearField (0/1)
//A) General case
// 1: Initial Longitudinal Obs. Coord. [m]
// 2: Final Longitudinal Obs. Coord. [m]
// 3: Number of points over Long. Obs. Coord.
// 4: Initial Photon Energy [eV]
// 5: Final Photon Energy [eV]
// 6: Number of points over Photon Energy
// 7: Initial Horizontal Obs. Coord. [m]
// 8: Final Horizontal Obs. Coord. [m]
// 9: Number of points over Hor. Obs. Coord.
// 10: Initial Vertical Obs. Coord. [m]
// 11: Final Vertical Obs. Coord. [m]
// 12: Number of points over Vert. Obs. Coord.
//B) Reduced case (one point assumed)
// 1: Photon Energy [eV]
// 2: Horizontal Obs. Coord. [m]
// 3: Vertical Obs. Coord. [m]
// 4: Longitudinal Obs. Coord. [m]

	waveHndl wField;
// Magnetic Field structure (text wave)

	waveHndl wElectronBeam;
// Electron Beam wave:
// 0: Energy [GeV]
// 1: Current [A]
// 2: s0 [m]
// 3: x0 [m]
// 4: dxds0 [r]
// 5: z0 [m]
// 6: dzds0 [r]
};

//*************************************************************************

struct srTIgorWfrFromTrjInputStruct {
	waveHndl wRad;
	waveHndl wAuxParam;
	waveHndl wIntPar;
	waveHndl wObservation;
	waveHndl wTrajectory;
};

//*************************************************************************

struct srTIgorIsotrSrcInputStruct {
	waveHndl wRad;
	waveHndl wAuxParam;
	waveHndl wObservation;
	waveHndl wIsotrSrc;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorGsnBeamInputStruct {
	waveHndl wRad;
	waveHndl wAuxParam;
	waveHndl wObservation;
	waveHndl wGsnBeam;
};

//*************************************************************************

struct srTIgorWfrSASEInputStruct {
	waveHndl wRad;
	waveHndl wControlSASE;
	waveHndl wPrecSASE;
	waveHndl wObservation;
	waveHndl wRadSASE;
	waveHndl wMagGen;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorPerStokesInputStruct {
	waveHndl wStokes;
	waveHndl wPerIntPar;
	waveHndl wObservation;
	waveHndl wPerField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorStokesWigInputStruct {
	waveHndl wStokes;
	waveHndl wWigIntPar;
	waveHndl wObservation;
	waveHndl wGenMagField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorStokesConstInputStruct {
	waveHndl wStokes;
	waveHndl wPrecPar;
	waveHndl wObservation;
	waveHndl wGenMagField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorStokesArbInputStruct {
	waveHndl wStokes;
	waveHndl wPrecPar;
	waveHndl wObservation;
	waveHndl wGenMagField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorPowDensInputStruct {
	waveHndl wPowDens;
	waveHndl wPowDensIntPar;
	waveHndl wObservation;
	waveHndl wGenMagField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorMagArb2Per {
	waveHndl wPerMagField;
	waveHndl wPrecPar;
	waveHndl wTrUnifMagField;
};

//*************************************************************************

struct srTIgorTrjInputStruct {
	waveHndl wOutVerCoor;
// Vertical trajectory coordinate [m] vs longitudinal coordinate [m]
	waveHndl wOutVerAng;
// Vertical trajectory angle [r] vs longitudinal coordinate [m]
	waveHndl wOutHorCoor;
// Horizontal trajectory coordinate [m] vs longitudinal coordinate [m]
	waveHndl wOutHorAng;
// Horizontal trajectory angle [r] vs longitudinal coordinate [m]

	waveHndl wField;
// Magnetic Field structure (text wave)

//	waveHndl wBzField;
// Vertical Magnetic Field wave:
// 0-...: Bz values [T]
// Number of points - in wave Dimensions
// sStart [m], sEnd [m] - in wave Scaling
//	waveHndl wBxField;
// Horizontal Magnetic Field wave:
// 0-...: Bx values [T]
// Number of points - in wave Dimensions
// sStart [m], sEnd [m] - in wave Scaling
// Number of points for Bx and Bz should be the same!

	waveHndl wElectronBeam;
// Electron Beam wave:
// 0: Energy [GeV]
// 1: Current [A]
// 2: s0 [m]
// 3: x0 [m]
// 4: dxds0 [r]
// 5: z0 [m]
// 6: dzds0 [r]
};

//*************************************************************************

struct srTIgorVersionStruct {
	Handle result;
// SRW Version Number
};

//*************************************************************************

struct srTHandleOfSRWRadStruct {
	waveHndl wRad;
// SRW Radiation structure
};

//*************************************************************************

struct srTHandleOfOneWaveStruct {
	DOUBLE p1;
	waveHndl WaveHndl;
// For operations on a single wave
};

//*************************************************************************

struct srTHandleOfReallyOneWaveStruct {
	waveHndl WaveHndl;
// For operations on a single wave
};

//*************************************************************************

struct srTIgorExtractPhaseInputStruct {
	waveHndl Out2WaveHndl;
	waveHndl Out1WaveHndl;
	waveHndl InWaveHndl;
};

//*************************************************************************

struct srTIgorRadResizeInputStruct {
	DOUBLE MethNo;
	DOUBLE pzd;
	DOUBLE pzm;
	DOUBLE pxd;
	DOUBLE pxm;
	waveHndl wRad;
};

//*************************************************************************

struct srTIgorObsSetNxNzInputStruct {
	DOUBLE p1;

	waveHndl wObservation;
	waveHndl wMagnField;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorRadPropagInputStruct {
	DOUBLE AuxPar1;
	DOUBLE MethNo;
	waveHndl wOptElem;
	waveHndl wRad;
};

//*************************************************************************

struct srTIgorWfrPropag {
	waveHndl wPrecPar;
	waveHndl wOptElem;
	waveHndl wRad;
};

//*************************************************************************

struct srTIgorRadPropagTestInputStruct {
	waveHndl wRadOut;
	waveHndl wOptElem;
	waveHndl wRadIn;
};

//*************************************************************************

struct srTIgorRadExtractInputStruct {
	waveHndl wExtractedData;
	waveHndl wExtractParam;
	waveHndl wRad;
};

//*************************************************************************

struct srTHandleOfTextStruct {
	waveHndl wText;
// Text wave
};

//*************************************************************************

struct srTIgorUtiInterTimeInputStruct {
	DOUBLE Delta;
};

//*************************************************************************

//struct srTIgorOptGenTransSetupInputStruct {
//	waveHndl wOptElem;
//};

//*************************************************************************

struct srTIgorOptMatConstInputStruct {
	DOUBLE PhotEn;
	DOUBLE Density;
	DOUBLE AtomNum; // Second arg
	DOUBLE MatNo; // First arg

	DOUBLE Delta; // return Re
	DOUBLE AttenLength; // return Im
};

//*************************************************************************

struct srTIgorOptZonePlateSpecSetupInputStruct {
	waveHndl wRad;
	waveHndl wOptElem;
};

//*************************************************************************

struct srTIgorUtiSpotInfoInputStruct {
	waveHndl wData;
	Handle SpotInfoName; // Return
};


//*************************************************************************

struct srTIgorUtiWrfLimitsInputStruct {
	DOUBLE PowLevel;
	waveHndl wData;
	Handle InfoName; // Return
};

//*************************************************************************

struct srTIgorUtiRemoveFlipsInputStruct {
	waveHndl wData;
};

//*************************************************************************

struct srTIgorUtiMagRadInputStruct {
	DOUBLE OutUnit; // 1- mm; 2- m; 3- km;
	DOUBLE ElecEnergy;
	DOUBLE Bconst;

	DOUBLE Rmag;
};

//*************************************************************************

struct srTIgorUtiMagCritPhotEnInputStruct {
	DOUBLE OutUnit; // 1- ; 2- ; 3- ; ...
	DOUBLE ElecEnergy;
	DOUBLE Bconst;

	DOUBLE CritPhotEn;
};

//*************************************************************************

struct srTIgorUtiUndKInputStruct {
	DOUBLE Period;
	DOUBLE Bpeak;

	DOUBLE K;
};

//*************************************************************************

struct srTIgorUtiUndFundPhotEnInputStruct {
	DOUBLE OutUnit; // 1- ; 2- ; 3- ; ...
	DOUBLE ElecEnergy;
	DOUBLE Period;
	DOUBLE Bpeak;

	DOUBLE FundPhotEn;
};

//*************************************************************************

struct srTIgorKmuInputStruct {
	DOUBLE x;
	DOUBLE mu;
	DOUBLE nInt;

	DOUBLE Kmu;
};

//*************************************************************************

struct srTIgorElecBeamPropag {
	waveHndl wOptElem;
	waveHndl wElectronBeam;
};

//*************************************************************************
//dddddddddddddddddd
struct srTIgorElecEnergyModul {
	//DOUBLE TransvScaleFact;
	//waveHndl wOptElem;
	//waveHndl wObservation;
	waveHndl wMagFld;
	waveHndl wElectronBeam;
};

//*************************************************************************
//TEST
struct srTIgor3DViewInputStruct {
	DOUBLE TransvScaleFact;
	waveHndl wOptElem;
	waveHndl wObservation;
	waveHndl wMagFld;
	waveHndl wElectronBeam;
};

//*************************************************************************
//TEST
struct srTIgorRadEmitPropagStokesMultiElecInputStruct {
	waveHndl wStokes;
	waveHndl wPpecPar;
	waveHndl wObservation;
	waveHndl wOptElem;
	waveHndl wMagField;
	waveHndl wElectronBeam;
};

//*************************************************************************
//TEST
struct srTIgorRadPropagStokesMultiElecInputStruct {
	waveHndl wStokes;
	waveHndl wPpecPar;
	waveHndl wOptElem;
	waveHndl wRad;
	waveHndl wElectronBeam;
};

//*************************************************************************
//TEST
struct srTIgorWfrEmitPropagInputStruct {
	waveHndl wRad; // Radiation structure (text wave) [out]
	waveHndl wPpecPar; // Precision parameters
	waveHndl wObservation; // Observation (rad. sampling)
	waveHndl wOptElem; // Optical element
	waveHndl wMagFld; // Magnetic field
	waveHndl wElectronBeam; // Electron beam
};

//*************************************************************************

struct srTIgorWfrReflectInputStruct {
	waveHndl wOptElem; //reflecting surfaces ???????????
	waveHndl wRad;
	waveHndl wElectronBeam;
};

//*************************************************************************

struct srTIgorRadFocusInputStruct {
	//waveHndl wOut2;
	//waveHndl wOut1;

	waveHndl wIntPar;
// Integration Parameters wave (this only concerns computation of incident SR; the integration over the Optical Element is currently steered by hand within the source...)
// 0: Integration Method (0- Manual; 1- AutoUndulator; 2- AutoWiggler)
// 1: Relative Precision or Step [m]
// 2: Initial point for the Integration sStart [m]
// 3: Final point for the Integration sEnd [m]
// 4: Maximum number of points to save in memory
// 5: Flag: Try to apply Near-Field residual term is the far-field one fails (1)

	waveHndl wObservation;
// Observation Wave Format2:
// 0: Photon Energy [eV]
// 1: X Coord. of the Obs. Centr. Point [m]
// 2: Y Coord. of the Obs. Centr. Point [m]
// 3: Z Coord. of the Obs. Centr. Point [m]
// 4: Initial Transverse Horizontal Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 5: Final Transverse Horizontal Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 6: Number of points over Transverse Horizontal Offset
// 7: Initial Longitudinal Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 8: Final Longitudinal Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 9: Number of points over Longitudinal Offset
// 10: Initial Transverse Vertical Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 11: Final Transverse Vertical Offset (with respect to optical axis) from the Obs. Centr. Point [m]
// 12: Number of points over Transverse Vertical Offset

	waveHndl wOpticalElem;
// Optical Element wave:
// 0: Optical Elem. ID Number (20- Spherical Mirror; 30- Zone Plate)
// 1: X Coord. of the Central Point [m]
// 2: Y Coord. of the Central Point [m]
// 3: Z Coord. of the Central Point [m]
// 4: Nx Compon. of the Central Normal Vector [no]
// 5: Ny Compon. of the Central Normal Vector [no]
// 6: Nz Compon. of the Central Normal Vector [no]
//If ID= 20 (Spherical Mirror):
// 7: Radius of Curvature [m]
// 8: Diameter (used as horizontal and vertical size at computation) [m]
//If ID= 30 (Zone Plate):
// Continue here

	waveHndl wField;
// Magnetic Field structure (text wave)

//	waveHndl wBzField;
// Vertical Magnetic Field wave:
// 0-...: Bz values [T]
// Number of points - in wave Dimensions
// sStart [m], sEnd [m] - in wave Scaling
//	waveHndl wBxField;
// Horizontal Magnetic Field wave:
// 0-...: Bx values [T]
// Number of points - in wave Dimensions
// sStart [m], sEnd [m] - in wave Scaling
// Number of points for Bx and Bz should be the same!

	waveHndl wElectronBeam;
// Electron Beam wave:
// 0: Energy [GeV]
// 1: Current [A]
// 2: s0 [m]
// 3: x0 [m]
// 4: dxds0 [r]
// 5: z0 [m]
// 6: dzds0 [r]
};

//*************************************************************************

#endif
