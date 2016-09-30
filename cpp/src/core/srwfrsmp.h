/************************************************************************//**
 * File: srwfrsmp.h
 * Description: Auxiliary structure for Wavefront Sampling (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRWFRSMP_H
#define __SRWFRSMP_H

#include <string.h>

//#include <cmath>
#include <math.h>
#include <ctype.h>
//#include "gmvect.h"
#include "gmtrans.h"

#include "srobject.h"

//*************************************************************************

enum srTDistanceOrAngle { Distance, Angle };

//*************************************************************************

enum srTDistrValType { FluxDens, FieldFourier, Flux, StokesParam };

//*************************************************************************

enum srTDistrPolariz { HorOnly, VerOnly, HorAndVer, TotOnly };

//*************************************************************************

enum srTIntegrDistrMap { MapVsHor, MapVsVer };

//*************************************************************************

enum srTCoordOrAngPresentation { CoordPres, AngPres };

//*************************************************************************

class srTWfrSmp : public CGenObject {
public:

	double LambStart, LambEnd, xStart, xEnd, yStart, yEnd, zStart, zEnd, tStart, tEnd;
	int nLamb, nx, ny, nz, nt;

	char AllowAutoChoiceOfNxNzForPropagat;
	double NxNzOversamplingParam;
	char DimensionsWereSetAuto;

	int xMult, zMult; // Used at comp. through Coord-Ang FFT
	TVector3d CenP, vHor, vLong; //vNormObsPl; // CenP is used only at Diffraction comp.
	bool obsPlaneIsTransv;

	double NormalizingConst;
	char LoopOrder[5];

	//double *m_pSurfData;
	//use smart pointer!
	CSmartPtr<double> m_spSurfData;

	int PresT; //0- Frequency Domain; 1- Time Domain
	int PhotonEnergyWavelengthUnits; // 0- keV, 1- eV, 2- Ang, 3- nm, 4- micron
	int BandwidthUnits; // 0- 0.1%, 1- keV, 2- eV, 3- Ang, 4- nm, 5- micron

	int CoordUnits; // 0- m, 1- mm
	char FluxComp; // if !=0, Integrated flux is computed

	char TreatLambdaAsEnergyIn_eV;
	char ShowPhaseOnly;
	char PhaseDerOrder;
	char InputWasModified, RadDistrDataContShouldBeRebuild;
	char OnlyOnePoint;
	char AssumeAllPhotonEnergies; //used to signal that power density will be computed

	srTDistanceOrAngle ArgType, DistrNormType;
	srTDistrValType DistrValType;
	srTDistrPolariz DistrPolariz;
	srTIntegrDistrMap IntegrDistrMap;

	srTCoordOrAngPresentation CoordOrAngPresentation;
	char AngPresToSpeedUpCoordPres;

	srTWfrSmp(double s, double hSt, double hFi, int hN, double vSt, double vFi, int vN, double* pSurfData, double eSt, double eFi, int eN, const char* PhotEnUnit, double tSt=0, double tFi=0, int tN=0, int presT =0, double* horOrtObsPlane =0, double* inNormObsPlane =0)
	{
		Initialize();

		LambStart = eSt; LambEnd = eFi; nLamb = eN;
		tStart = tSt; tEnd = tFi; nt = tN; PresT = presT;

		yEnd = yStart = s; ny = 1;
		xStart = hSt; xEnd = hFi; nx = hN;
		zStart = vSt; zEnd = vFi; nz = vN;

		if(PhotEnUnit != 0)
		{
			char cPhotEnBuf[100];
			strcpy(cPhotEnBuf, PhotEnUnit);
			for(int i=0; i<(int)strlen(PhotEnUnit); i++) cPhotEnBuf[i] = toupper(cPhotEnBuf[i]);
			if(strcmp(cPhotEnBuf, "EV") == 0) 
			{
                PhotonEnergyWavelengthUnits = 1; // 0- keV, 1- eV, 2- Ang, 3- nm, 4- micron
                TreatLambdaAsEnergyIn_eV = 1;
			}
			else if(strcmp(cPhotEnBuf, "NM") == 0)
			{
                PhotonEnergyWavelengthUnits = 3; // 0- keV, 1- eV, 2- Ang, 3- nm, 4- micron
                TreatLambdaAsEnergyIn_eV = 0;
			}
		}

        CoordUnits = 0; // 0- m, 1- mm
        BandwidthUnits = 0; // 0- 0.1%, 1- keV, 2- eV, 3- Ang, 4- nm, 5- micron

		strcpy(LoopOrder, "yzxw");
        CoordOrAngPresentation = CoordPres;
        
		//obsPlaneIsTransv = true;
		//vNormObsPl.x = vNormObsPl.z = 0.;
		//vNormObsPl.y = -1.;

		vHor.x = 1; vHor.y = vHor.z = 0;
		vLong.x = 0; vLong.y = 1; vLong.z = 0;

		if(horOrtObsPlane != 0)
		{
			vHor.x = horOrtObsPlane[0]; vHor.y = horOrtObsPlane[1]; vHor.z = horOrtObsPlane[2]; 
			vHor.Normalize();
			//obsPlaneIsTransv = ObsPlaneIsTransverse();
		}
		if(inNormObsPlane != 0)
		{
			vLong.x = inNormObsPlane[0]; vLong.y = inNormObsPlane[1]; vLong.z = inNormObsPlane[2]; 
			vLong.Normalize();
		}

		if(pSurfData != 0)
		{//allocate and copy observ. surface data
			//long totNumPtSurf = hN*vN;
			long long totNumPtSurf = ((long long)hN)*((long long)vN);
			double *pSurfDataLoc = new double[totNumPtSurf];
			double *t_pSurfDataLoc = pSurfDataLoc, *t_pSurfData = pSurfData;
			//for(long i=0; i<totNumPtSurf; i++) 
			for(long long i=0; i<totNumPtSurf; i++) 
			{
				*(t_pSurfDataLoc++) = *(t_pSurfData++);
			}

			CSmartPtr<double> hSurfData(pSurfDataLoc);
			m_spSurfData = hSurfData;
		}

		obsPlaneIsTransv = ObsPlaneIsTransverse(); //to call at the end only

		//if(normObsPlane != 0)
		//{
		//	vNormObsPl.x = normObsPlane[0]; vNormObsPl.y = normObsPlane[1]; vNormObsPl.z = normObsPlane[2]; 
		//	vNormObsPl.Normalize();
		//	obsPlaneIsTransv = ObsPlaneIsTransverse();
		//}
		//RadDistrDataContShouldBeRebuild = 0;
		//PhotonEnergyWavelengthUnits = 1; // We support only eV for the moment
		//TreatLambdaAsEnergyIn_eV = 1;
	}
	srTWfrSmp()
	{
		Initialize();
	}
	~srTWfrSmp()
	{
		//if(m_pSurfData != 0) delete[] m_pSurfData;
	}

	void Initialize()
	{
		CoordOrAngPresentation = CoordPres;
		DistrValType = StokesParam;
		DistrPolariz = HorAndVer;
		TreatLambdaAsEnergyIn_eV = 0;

		LambStart = LambEnd = 0.; nLamb = 1;
		xStart = xEnd = 1.E+23; nx = 1;
		zStart = zEnd = 1.E+23; nz = 1;
		yStart = yEnd = 1.E+23; ny = 1;
		tStart = tEnd = 0; nt = 0;

		PresT = 0;

		InputWasModified = 1;
		ShowPhaseOnly = 0;
		RadDistrDataContShouldBeRebuild = 0;
		OnlyOnePoint = 0;

		AngPresToSpeedUpCoordPres = 0;
		DimensionsWereSetAuto = 0;

		AllowAutoChoiceOfNxNzForPropagat = 0;
		NxNzOversamplingParam = 1;

		CoordUnits = 1;
		FluxComp = 0;

		AssumeAllPhotonEnergies = 0;
		obsPlaneIsTransv = true;
		//m_pSurfData = 0;
	}

	bool IsDefined()
	{//conherent with Initialize
		if((LambStart != 0.) || (LambEnd != 0.) || (nLamb != 1) ||
		   (xStart < 0.99999E+23) || (xStart > 1.00001E+23) || (xEnd < 0.99999E+23) || (xEnd > 1.00001E+23) || (nx != 1) ||
		   (zStart < 0.99999E+23) || (zStart > 1.00001E+23) || (zEnd < 0.99999E+23) || (zEnd > 1.00001E+23) || (nz != 1) ||
		   (yStart < 0.99999E+23) || (yStart > 1.00001E+23) || (yEnd < 0.99999E+23) || (yEnd > 1.00001E+23) || (ny != 1) ||
		   (tStart != 0.) || (tEnd != 0.) || (nt != 0)) return true;
		else return false;
	}

	char SetupLambda(double Start, double End)
	{// Makes eV or nm from any units
		if(PhotonEnergyWavelengthUnits < 2) TreatLambdaAsEnergyIn_eV = 1;

		double OldLambStart = LambStart, OldLambEnd = LambEnd;

		if(PhotonEnergyWavelengthUnits == 0) { LambStart = Start*1000.; LambEnd = End*1000.;}
		else if(PhotonEnergyWavelengthUnits == 1) { LambStart = Start; LambEnd = End;}
		else if(PhotonEnergyWavelengthUnits == 2) { LambStart = Start*0.1; LambEnd = End*0.1;}
		else if(PhotonEnergyWavelengthUnits == 3) { LambStart = Start; LambEnd = End;}
		else if(PhotonEnergyWavelengthUnits == 4) { LambStart = Start*1000.; LambEnd = End*1000.;}

		double AbsPrec = LambStart*1.E-09;
		return ((::fabs(OldLambStart - LambStart) < AbsPrec) && (::fabs(OldLambEnd - LambEnd) < AbsPrec))? 0 : 1;
		// 1 if was modified
	}
	void RetrievePhotonEnergyOrWavelength(double& Start, double& End)
	{
		if(PhotonEnergyWavelengthUnits == 0) Start = LambStart;
		else if(PhotonEnergyWavelengthUnits == 1) Start = LambStart*1000.;
		else if(PhotonEnergyWavelengthUnits == 2) Start = LambStart*10.;
		else if(PhotonEnergyWavelengthUnits == 3) Start = LambStart;
		else if(PhotonEnergyWavelengthUnits == 4) Start = LambStart/1000.;

		if(PhotonEnergyWavelengthUnits == 0) End = LambEnd;
		else if(PhotonEnergyWavelengthUnits == 1) End = LambEnd*1000.;
		else if(PhotonEnergyWavelengthUnits == 2) End = LambEnd*10.;
		else if(PhotonEnergyWavelengthUnits == 3) End = LambEnd;
		else if(PhotonEnergyWavelengthUnits == 4) End = LambEnd/1000.;
	}
	void EnsureZeroTransverseRangesForSinglePoints()
	{
        if(nx == 1)
        {
            //if((!DistrInfoDat.FluxComp) && (!DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat))
            //{
            double xMid = 0.5*(xStart + xEnd);
            xStart = xMid; xEnd = xMid;
            //}
        }
        if(nz == 1)
        {
            //if((!DistrInfoDat.FluxComp) && (!DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat))
            //{
            double zMid = 0.5*(zStart + zEnd);
            zStart = zMid; zEnd = zMid;
            //}
        }
	}
	void OutRangesForFluxComp(double& dx, double& dz)
	{
        if(nx == 1)
        {
            dx = xEnd - xStart;
		}
		else if(nx > 1)
		{
            dx = (xEnd - xStart)/(nx - 1);
		}

        if(nz == 1)
        {
            dz = zEnd - zStart;
		}
		else if(nz > 1)
		{
            dz = (zEnd - zStart)/(nz - 1);
		}
	}
	
	bool ObsPlaneIsTransverse() //and not rotated
	{
		const double relTol = 1.e-10;

		if(m_spSurfData.rep != 0) return false; //make more consistent test of the surface

		//if(vNormObsPl.isZero() || ((::abs(vNormObsPl.x) < 1e-10) && (::abs(vNormObsPl.y + 1) < 1e-10) && (::abs(vNormObsPl.z) < 1e-10)))
		//{
		//	return true;
		//}
		if((vLong.isZero() || ((::fabs(vLong.x) < relTol) && (::fabs(vLong.y - 1) < relTol) && (::fabs(vLong.z) < relTol))) &&
		   (vHor.isZero() || ((::fabs(vHor.x - 1) < relTol) && (::fabs(vHor.y) < relTol) && (::fabs(vHor.z) < relTol))))
		{
			return true;
		}
		else return false;
	}

	bool SetupTrfObsPlaneIfNecessary(gmTrans& trfObsPl)
	{
		obsPlaneIsTransv = ObsPlaneIsTransverse();
		if(obsPlaneIsTransv) return false;

		//TVector3d vEx(1, 0, 0), vEy(0, 1, 0), vEz(0, 0, 1), vEyP(-vNormObsPl.x, -vNormObsPl.y, -vNormObsPl.z), vEzP;
		//TVector3d vExP = vEx - (vEx*vNormObsPl);

		TVector3d vEx(1, 0, 0), vEy(0, 1, 0), vEz(0, 0, 1);
		TVector3d vExP = vHor, vEyP = vLong, vEzP;

		vExP.Normalize();
		vEyP.Normalize();

		if(vExP.AmpE2() < 1e-20) //vExP is zero
		{
			//vEzP = vEz - (vEz*vNormObsPl);
			vEzP = vEz - (vEz*vLong);

			vEzP.Normalize();
			vExP = vEyP^vEzP;
		}
		else
		{
			vEzP = vExP^vEyP;
		}
		
		TVector3d mRow1(vExP*vEx, vEyP*vEx, vEzP*vEx);
		TVector3d mRow2(vExP*vEy, vEyP*vEy, vEzP*vEy);
		TVector3d mRow3(vExP*vEz, vEyP*vEz, vEzP*vEz);
		TMatrix3d M(mRow1, mRow2, mRow3);
		TVector3d vCen(0.5*(xStart + xEnd), yStart, 0.5*(zStart + zEnd));

		//gmTrans trfAux(M, vCen);
		//trfObsPl = trfAux;
		trfObsPl.SetMatrixVector(M, vCen);
		return true;
	}
};

//*************************************************************************

#endif
