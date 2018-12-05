/************************************************************************//**
 * File: srerror.cpp
 * Description: Error and Warning Messages
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srerror.h"

//using namespace std;
//#include <vector>

//-------------------------------------------------------------------------

vector<string> CErrWarn::error;
vector<string> CErrWarn::warning;

//-------------------------------------------------------------------------

CErrWarn::CErrWarn()
{
//string CErrWarn::error[] = {
	error.push_back("Wrong error number"); //to check if and how this is used in SRW for Igor

	error.push_back("WaveAccess requires Igor Pro 2.0 or later.\0");		// OLD_IGOR
	error.push_back("Wave does not exist.\0");							// NON_EXISTENT_WAVE
	error.push_back("This function requires 3D wave(s).\0");				// NEEDS_3D_WAVE

	error.push_back("Incorrect Electron Beam Energy value. The Electron Beam Energy should be within 0 and 1 000 000 Gev.\0");
	error.push_back("This function requires wave(s) of type Single Float 32 bit.\0");
	error.push_back("Incorrect Electron Beam wave format.\0");
	error.push_back("Incorrect Magnetic Field wave format.\0");
	error.push_back("Memory allocation failure.\0");
	error.push_back("Inconsistent definition of Horizontal and Vertical Magnetic Field.\0");
	error.push_back("Incorrect Photon Energy or Wavelength value.\0");
	error.push_back("Incorrect Observation wave format.\0");
	error.push_back("Incorrect Precision wave format.\0");
	error.push_back("Incorrect Integration Method number. The Integration Method number can be 0 (manual), 1 or 2 (automatic).\0");
	error.push_back("Incorrect Relalative Precision or Integration Step value.\0");
	error.push_back("Radiation Integration limits are inconsistent with the Magnetic Field definition limits.\0");
	error.push_back("Incorrect Stokes Parameters wave format.\0");
	error.push_back("This function requires wave(s) of type Double Float 64 bit.\0");
	error.push_back("Incorrect Electric Field Fourier Components wave format.\0");
	error.push_back("This function requires wave(s) of type Complex 2 x 64 bit.\0");
	error.push_back("This function requires wave(s) of type Double Float 64 bit or Complex 2 x 32 bit.\0");
	error.push_back("Incorrect dimensions of the Radiation output wave.\0");
	error.push_back("Computation aborted by user.\0");
	error.push_back("The amount of memory allowed for use by SRW is insufficient for SR computation with current input.\0");

	error.push_back("SRW can not compute this case. Try to modify Magnetic Field definition range and/or Longitudinal Integration limits.\0");
	error.push_back("SRW can not compute this case. Electron Trajectory passes too closely to the Observation Point.\0");
	error.push_back("SRW can not compute this case. Longitudinal position of the Observation Plane is within the Integration limits.\0");
	error.push_back("SRW can not compute this case. Electron Beam is not ultra-relativistic.\0");
	error.push_back("Initial Coordinates/Angles of the Electron Trajectory are specified for the Longitudinal position which is out of Magnetic Field limits.\0");
	error.push_back("SRW can not compute this case. The Wavelength is comparable with the Observation Distance.\0");
	error.push_back("SRW can not compute this case. Angle between instant Electron Velocity and direction to the Observation Point is too large.\0");

	error.push_back("This function requires Text wave(s).\0");
	error.push_back("Incorrect wavefront structure.\0");
	error.push_back("The operation needs another presentation of the Radiation data.\0");

	error.push_back("This function requires Numeric wave(s).\0");
	error.push_back("This function requires 2D wave(s).\0");
	error.push_back("The wave dimensions should be even.\0");
	error.push_back("This function requires wave of type Complex 2 x 32 bit.\0");
	error.push_back("This function requires 1D Text wave.\0");
	error.push_back("Unknown Optical Element or incorrect Optical Element definition.\0");
	error.push_back("Internal error in FFT procedure.\0");
	error.push_back("This function requires 1D or 2D wave(s).\0");
	error.push_back("This function requires waves of equal sizes.\0");
	error.push_back("Error in Optical Element definition.\0");
	error.push_back("This function requires 1D wave.\0");
	error.push_back("This operation requires more than one Observation Point in Hor. and Vert. directions.\0");
	error.push_back("This operation accepts no more than one Drift Space, which can be only the last component in a Container.\0");

	error.push_back("Error at reading a number from a text wave.\0");
	error.push_back("Can't find a harmonic wave of the Periodic Magnetic Field structure.\0");
	error.push_back("Improper or corrupt Stokes Parameters structure.\0");
	error.push_back("Improper or corrupt Periodic Magnetic Field Harmonic structure.\0");
	error.push_back("Sorry, this kind of magnet structure is not supported yet.\0");

	error.push_back("Electron beam divergence and slit size should not be zero at the same time for this method of computation.\0");

	error.push_back("Internal error in linear fitting procedure.\0");
	error.push_back("Improper or corrupt Optical Component structure.\0");

	error.push_back("Improper or corrupt Power Density structure.\0");
	error.push_back("SRW can not compute this case in that mode. Try to reduce Magnetic Field limits and/or Observation range and/or Precision level.\0");

	error.push_back("This function requires wave of type Complex 2 x 64 bit.\0");
	error.push_back("Propagation can not be performed correctly due to severe memory limitation.\0");
	error.push_back("For correct propagation, number of points of should be larger than 1 both in horizontal and vertical direction.\0");
	error.push_back("With this method, the propagation can not be done for several photon energies at once.\0");

	error.push_back("For this material, the density value should be specified.\0");
	error.push_back("Material with this atomic number is not supported.\0");

	error.push_back("This function requires wave of type Single (32 bit) or Double Precision (64 bit).\0");
	error.push_back("This function requires 1D or 2D wave.\0");

	error.push_back("Thin electron beam parameters were not properly set up.\0");
	error.push_back("Thick electron beam parameters were not properly set up.\0");

	error.push_back("Number of integrations of Modified Bessel function can not be negative.\0");
	error.push_back("Magnetic field can not be zero for this mode of computation.\0");

	error.push_back("Improper or corrupt Isotropic Source structure.\0");
	error.push_back("Source Power should be positive.\0");
	error.push_back("Number of photons should be positive.\0");
	error.push_back("Improper or corrupt Gaussian Beam structure.\0");
	error.push_back("Gaussian beam waist size can not be negative.\0");
	error.push_back("The order of a Gaussian Beam mode can not be negative.\0");
	error.push_back("Number of photons should be positive.\0");
	error.push_back("Polarization component identifier should be an integer number from 1 to 6.\0");

	error.push_back("Improper or corrupt Trajectory structure.\0");
	error.push_back("Improper or corrupt Trajectory component wave.\0");
	error.push_back("Electron Beam Energy should be positive.\0");
	error.push_back("Trajectory components were not properly set up.\0");
	error.push_back("Contradictory definition of trajectory components.\0");
	error.push_back("Initial longitudinal position is out of trajectory definition limits.\0");

	error.push_back("The number of magnetic field harmonics should not exceed 20.\0");
	error.push_back("Unknown magnet element.\0");

	error.push_back("***** No wiggler field\0");
	error.push_back("***** Invalid wiggler type\0");
	error.push_back("***** Wiggler period must be positive\0");
	error.push_back("***** Electron beam energy is too small\0");
	error.push_back("***** Number of particles is not multiple of 8\0");
	error.push_back("***** Radiation wavelength must be positive\0");
	error.push_back("***** Radiation power must be positive. Set the initial power to small value for SASE.\0");
	error.push_back("***** Integration step (DELZ) must be 0.5 if field errors are included\0");
	error.push_back("***** Parameter ISCAN must be smaller than 7\0");
	error.push_back("***** Parameter scan not defined\0");
	error.push_back("***** Time dependency and parameter scan are selected simultaneously\0");
	error.push_back("***** Step size greater than slice separation\0");
	error.push_back("***** No module defined. NWIG must pe positive.\0");
	error.push_back("***** Focusing length is not an integer of step size\0");
	error.push_back("***** Defocusing length is not an integer of step size\0");
	error.push_back("***** Drift length is not an integer of step size\0");
	error.push_back("***** Start of FODO lattice is not an integer of step size\0");
	error.push_back("***** Slice separation is not an integer of step size\0");
	error.push_back("***** Too many integration steps\0");
	error.push_back("***** Too many slices\0");
	error.push_back("***** Number of grid points of Cartesian mesh is not odd\0");

	error.push_back("Non-compatible Wavefront and Stokes structures.\0");
	error.push_back("Electron beam structure was not set up.\0");
	error.push_back("Distance from the Source to the Observation Plane should be positive.\0");
	error.push_back("Photon Energy should be positive.\0");
	error.push_back("Memory allocation failure. Please note that SASE computation may require ~50 MB of RAM.\0");

	error.push_back("Incorrect arguments for the 3D Viewing function.\0");
	error.push_back("Incorrect function arguments.\0");

	error.push_back("Function arguments: array is not defined correctly.\0");
	error.push_back("Object with this index does not exist.\0");
	error.push_back("Object with this index is not a magnetic field.\0");
	error.push_back("Object with this index is not an electron beam.\0");
	error.push_back("Object with this index is not a wave front sampling structure.\0");
	error.push_back("Incorrect parameters submitted to the SR computation function.\0");
	error.push_back("This type of SR computation is not implemented for given magnetic field.\0");
	error.push_back("Incorrect or no precision parameters were submitted to SR computation function.\0");
	error.push_back("Object with this index is not a radiation wave front.\0");
	error.push_back("Incorrect parameters submitted to the wavefront component extraction function.\0");
	error.push_back("Periodic magnetic field harmonic number should be positive.\0");
	error.push_back("Magnetic field period should be positive.\0");
	error.push_back("Periodic magnetic field length should be positive.\0");
	error.push_back("Number of magnetic field harmonics should be more than 0.\0");
	error.push_back("Incorrect wavefront structure.\0");
	error.push_back("Object with this index is not a gaussian beam.\0");
	error.push_back("Drift space is expected. Object with this index is not a drift space.\0");
	error.push_back("Incorrect relative precision. The relative precision should be a real number between 0 and 1.\0");
	error.push_back("Object with this index is not a transversely uniform magnetic field.\0");
	error.push_back("1D interpolating structure was not set up prperly.\0");
	error.push_back("Object with this index is not a periodic magnetic field.\0");
	error.push_back("Zero pointer is submitted in place of a pointer to an integer number.\0");
	error.push_back("Incorrect grid was specified for the wavefront.\0");
	error.push_back("Can not find magnetic field wave.\0");
	error.push_back("Empty magnetic field container encountered.\0");
	error.push_back("No non-zero magnetic field is defined.\0");
	error.push_back("No Stokes parameters structure found.\0");
	error.push_back("Overlapping magnet elements encountered.\0");
	error.push_back("Incorrect longitudinal step or relative precision value.\0");
	error.push_back("Step of longitudinal integration is too large.\0");
	error.push_back("Incorrect limits of longitudinal integration for power density calculation.\0");
	error.push_back("Incorrect or incompatible data structures were supplied for integration of radiation components.\0");
	error.push_back("Incorrect definition of electron beam.\0");
	error.push_back("Failed to create wave.\0");
	error.push_back("Object with this index is not an optical element.\0");
	error.push_back("Can not extract transmission characteristic of this optical element.\0");
	error.push_back("Incorrect definition of parameters for transmission characteristic.\0");
	error.push_back("Number of points should be positive.\0");
	error.push_back("Maximum number of step subdivisions reached at automatic Runge-Kutta integration.\0");
    error.push_back("Step size is too small in automatic Runge-Kutta integration routine.\0");
    error.push_back("Incorrect definition of Gaussian beam power (energy).\0");
    error.push_back("Time-domain radiation (wavefront) structure is required.\0");
	error.push_back("Character identifying electric field representation should be either 'F' (for frequency domain) or 'T' (for time domain) ot 'C' (for coordinate reperesentation) or 'A' (for angular representation)\0"); //#155

	error.push_back("***** Specified number of slices and/or number of particles in slice is inconsistent with submitted electron distribution.\0");

	error.push_back("Incorrect definition of centers and RMS sizes of multi-dimensional Gaussians.\0");
	error.push_back("The submitted wave has too small number of elements.\0"); //copied from IDB (+13)
	error.push_back("Incorrect (or inconsistent) definition of surface data in radiation sampling structure.\0");

	error.push_back("Incorrect or no data structures were supplied to function.\0");
	error.push_back("Magnetic field supplied is not supported in this type of calculation.\0");
	error.push_back("Incorrect trajectory structure / parameters were supplied to function.\0");
	error.push_back("Incorrect or insufficient parameters for trajectory calculation.\0");
	error.push_back("Incorrect or insufficient parameters for synchrotron radiation calculation.\0");
	error.push_back("Incorrect or insufficient parameters for Gaussian beam definition.\0");
	error.push_back("Incorrect or insufficient parameters for intensity extraction.\0");
	error.push_back("Incorrect or no wavefront structure.\0");
	error.push_back("External (callback) function falied to modify (/ reallocate memory for) wavefront data.\0");
	error.push_back("Incorrect or insufficient parameters for wavefront resizing.\0");
	error.push_back("Incorrect or insufficient parameters for changing electric field (wavefront) representation.\0");
	error.push_back("Incorrect or insufficient parameters for electric field (wavefront) propagation.\0");
	error.push_back("Pointer to external (callback) wavefront modification function is not defined.\0");

	error.push_back("Incorrect optical element orientation.\0");
	error.push_back("Incorrect optical element simulation method number.\0");

	error.push_back("Incorrect or insufficient parameters for synchrotron radiation power density calculation.\0");

	error.push_back("Incorrect ellipsoidal mirror parameters: p, q, grazing angle and sagital radius should be positive.\0"); //#176
	error.push_back("Failed to determine optical axis after reflection from mirror.\0"); //#177
	error.push_back("Failed to interpolate electric field.\0"); //#178

	error.push_back("Incorrect or insufficient parameters for magnetic field calculation.\0"); //#179
	error.push_back("Incorrect input parameters for FFT procedure.\0"); //#180
	error.push_back("Incorrect input parameters for convolution with Gaussian.\0"); //#181
	error.push_back("Incorrect input parameters for conversion of tabulated magnetic field to periodic.\0"); //#182
	error.push_back("Incorrect spherical mirror parameter: radius should not be zero.\0"); //#183
	error.push_back("Incorrect crystal parameter: use case is not supported.\0"); //#184
	error.push_back("Incorrect input parameters for search of indexes of undulator gap and phase values requiired for undulator field interpolation based on gap and phase.\0"); //#185
	error.push_back("Incorrect or insufficient parameters for spherical wave electric field calculation.\0"); //#186
	
	error.push_back("Failed to determine array element index for interpolation.\0"); //#187
	error.push_back("Failed to allocate array in front-end / client environment.\0"); //#188
	error.push_back("Mutual intensity can not be extracted for these input parameters.\0"); //#189
	error.push_back("Incorrect input parameters for calculation of statistical characteristics of intensity.\0"); //#190
	error.push_back("Incorrect input parameters for processing intensity distributions.\0"); //#191

//};

//string CErrWarn::warning[] = {
	warning.push_back("Wrong warning number");
	warning.push_back("Too large taper parameter. Resulting precision may be poor.");
	warning.push_back("Too large Optical Klystron phase shift parameter. Resulting precision may be poor.");
	warning.push_back("Too few periods in the undulator. Resulting precision may be poor.");
	warning.push_back("Electron Beam Emittance was not taken into account.");
	warning.push_back("Phase Shift tabulation is not sufficiently dense. Optical element is not well resolved.");
	warning.push_back("Propagation accuracy may appear reduced. All available RAM was used at the computation.");
	warning.push_back("Propagation precision may be poor.");
	warning.push_back("Attenuation Length precision may be poor.");
	warning.push_back("Emission conditions do not match well the wiggler case. The resulting precision may be poor.");
	warning.push_back("Unable to extract the Multi-electron radiation component. Single-electron component was extracted.");
	warning.push_back("COMPUTATION WAS NOT PERFORMED and the result was set to 0 for some values of input parameters \r due to the following problem: Electron Trajectory passes too closely to the Observation Point.");
	warning.push_back("Computation of terminating terms of radiation integrals WAS NOT PERFORMED for some values of input parameters \r because asymptotic expansion validity criterion was not satisfied. This may influence the accuracy of the computation. \r One can try to fix the problem by modifying Magnetic Field definition range and/or Integration limits."); //don't make it longer!
	warning.push_back("Too many particles requested. The number of particles was set to the default maximum value.");
	warning.push_back("Same base for different loading. Variables are highly correlated.");
	warning.push_back("Too many grid points. The number of grid points was set to the default maximum value.");
	warning.push_back("Too many radial grid points. The number of radial grid points was set to the default maximum value.");
	warning.push_back("Not enough radial grid points for space charge computation. \r The number of radial grid points was set to the default maximum value.");
	warning.push_back("Too many longitudinal modes for space charge computation. \r The number of longitudinal modes was set to 3.");
	warning.push_back("The specified maximum number of magnetic field harmonics is too large. \r The maximum number of harmonics was set to 20.");
	warning.push_back("No magnetic field harmonics suitable for setting up periodic field structure were found.");
	warning.push_back("Spectral flux can only be extracted vs photon energy.");
	warning.push_back("Longitudinal position of the observation point(s) is too close to (or within) the longitudinal integration limits. The resulting precision may be poor.");
	warning.push_back("Electron beam is not ultra-relativistic. The resulting accuracy may be poor.");
	warning.push_back("To calculate emission at harmonics (HGHG), this version of GENESIS requires Electron Distribution data."); //#define GENESIS_RAD_HARM_CALC_NEEDS_ELEC_DISTRIB 24 + SRW_WARNINGS_OFFSET
//};
}

//-------------------------------------------------------------------------

int CErrWarn::ValidateArray(void* Arr, int nElem)
{
	if(Arr == 0) return INCORRECT_ARGUMENTS_ARRAY;
	if(nElem < 0) return INCORRECT_ARGUMENTS_ARRAY;
	return 0;
}

//-------------------------------------------------------------------------
