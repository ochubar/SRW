/************************************************************************//**
 * File: srinterf.h
 * Description: C/C++ API of Synchrotron Radiation Workshop (SRW) header - OBSOLETE version
 * Project: Synchrotron Radiation Workshop
 * First release: October 2002
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRINTERF_H
#define __SRINTERF_H

#include "alphadll.h"

#ifdef __cplusplus  
extern "C" {
#endif

//-------------------------------------------------------------------------

struct srTParPrecElecFld {
	int IntegMethNo; 
	double RelPrecOrStep;
	double sStartInt, sEndInt;
	double NxNzOversamplingFactor; //active if > 0
	bool ShowProgrIndic;
	char CalcTerminTerms;

	//srTParPrecElecFld(int In_CreateNewWfrObj, int In_IntegMethNo, double In_RelPrecOrStep, double In_sStartInt, double In_sEndInt, double In_NxNzOversamplingFactor)
	//srTParPrecElecFld(int In_IntegMethNo, double In_RelPrecOrStep, double In_sStartInt, double In_sEndInt, double In_NxNzOversamplingFactor, bool In_ShowProgrIndic = true)
	srTParPrecElecFld(int In_IntegMethNo, double In_RelPrecOrStep, double In_sStartInt, double In_sEndInt, double In_NxNzOversamplingFactor, bool In_ShowProgrIndic = true, char In_CalcTerminTerms = 1)
	{
        //CreateNewWfrObj = In_CreateNewWfrObj;
        IntegMethNo = In_IntegMethNo; 
        RelPrecOrStep = In_RelPrecOrStep;
        sStartInt = In_sStartInt; sEndInt = In_sEndInt;
        NxNzOversamplingFactor = In_NxNzOversamplingFactor;
		ShowProgrIndic = In_ShowProgrIndic;
		CalcTerminTerms = In_CalcTerminTerms;
	}
};

//-------------------------------------------------------------------------

struct srTParPrecStokesPer {
	int InitHarm;
	int FinHarm;
	double PrecS;
	double PrecPhi;
	char IntOrFlux;
	double MinPhotEnExtRight; //OC170713

	srTParPrecStokesPer(int In_InitHarm, int In_FinHarm, double In_PrecS, double In_PrecPhi, char In_IntOrFlux, double In_MinPhotEnExtRight=1)
	{
        InitHarm = In_InitHarm;
		FinHarm = In_FinHarm;
		PrecS = In_PrecS;
		PrecPhi = In_PrecPhi;
		IntOrFlux = In_IntOrFlux;
		MinPhotEnExtRight = In_MinPhotEnExtRight;
	}
};

//-------------------------------------------------------------------------

struct srTParPrecStokesArb {
	char IntOrFlux;
	int MethNo; 
	double RelPrecOrStep;
	int NumIter;
};

//-------------------------------------------------------------------------

struct srTParPrecElecFldGaus {
	//int CreateNewWfrObj; //0- don't create, otherwise - create
	double NxNzOversamplingFactor; //active if > 0

	//srTParPrecElecFldGaus(int In_CreateNewWfrObj, double In_NxNzOversamplingFactor)
	srTParPrecElecFldGaus(double In_NxNzOversamplingFactor)
	{
        //CreateNewWfrObj = In_CreateNewWfrObj;
        NxNzOversamplingFactor = In_NxNzOversamplingFactor;
	}
};

//-------------------------------------------------------------------------

struct srTParPrecPowDens {
	//int CreateNewPowDensObj; //0- don't create, otherwise - create
	double PrecFact;
	int MethNo, UseSpecIntLim; 
	double sIntStart, sIntFin;

	//srTParPrecPowDens(int In_CreateNewPowDensObj, int In_MethNo, double In_PrecFact)
	srTParPrecPowDens(int In_MethNo, double In_PrecFact, int In_UseSpecIntLim, double In_sIntStart, double In_sIntFin)
	{
        //CreateNewPowDensObj = In_CreateNewPowDensObj; //0- don't create, otherwise - create
        MethNo = In_MethNo; 
        PrecFact = In_PrecFact;
		UseSpecIntLim = In_UseSpecIntLim;
		sIntStart = In_sIntStart; sIntFin = In_sIntFin;
	}
};

//-------------------------------------------------------------------------

struct srTParIntensExtract {
	int PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
	int Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E)
	int PlotType; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
	int TransvPres; // 0: Spatial; 1: Angular
	double e, x, z;

	srTParIntensExtract(int In_PolarizCompon, int In_Int_or_Phase, int In_PlotType, int In_TransvPres, double In_e, double In_x, double In_z)
	{
		PolarizCompon = In_PolarizCompon; // 0: Linear Hor.; 1: Linear Vert.; 2: Linear 45; 3: Linear 135; 4: Circul. Right; 5: Circul. Left; 6: Total
		Int_or_Phase = In_Int_or_Phase; // 0: 1-e Int; 1: Multi-e Int; 2: Phase; 3: Re(E)
		PlotType = In_PlotType; // vs 0: e; 1: x; 2: z; 3: x&z; 4: e&x; 5: e&z; 6: e&x&z
		TransvPres = In_TransvPres; // 0: Spatial; 1: Angular
		e = In_e; x = In_x; z = In_z;
	}
};

//-------------------------------------------------------------------------

struct srTParPrecWfrPropag {
	char MethNo;
	char UseResBefore;
	char UseResAfter;
    double PrecFact;
    double UnderSampThresh;
	char AnalTreatment; //0- none; 1- linear term; 2- quadratic term
	char DoNotResetAnalTreatTermsAfterProp;

	//OC011213
	double vLxOut, vLyOut, vLzOut; //Coordinates of the output Optical Axis vector
	double vHxOut, vHyOut; //Coordinates of the Horizontal Base vector of the output frame

	srTParPrecWfrPropag(char In_MethNo, char In_UseResBefore, char In_UseResAfter, double In_PrecFact, double In_UnderSampThresh, 
		char In_AnalTreatment =0, char In_DoNotResetAnalTreatTermsAfterProp =0,
		double In_vLxO =0, double In_vLyO =0, double In_vLzO =0, double In_vHxO =0, double In_vHyO =0)
	{
		MethNo = In_MethNo;
        UseResBefore = In_UseResBefore;
        UseResAfter = In_UseResAfter;
        PrecFact = In_PrecFact;
        UnderSampThresh = In_UnderSampThresh;
		AnalTreatment = In_AnalTreatment; //Allow switching to under-sampling mode
		DoNotResetAnalTreatTermsAfterProp = In_DoNotResetAnalTreatTermsAfterProp;

		vLxOut = In_vLxO;
		vLyOut = In_vLyO;
		vLzOut = In_vLzO;
		vHxOut = In_vHxO;
		vHyOut = In_vHyO;
	}
	srTParPrecWfrPropag()
	{
        MethNo = 2;
        UseResBefore = UseResAfter = 1;
        PrecFact = 1.;
        UnderSampThresh = 0.5;
		AnalTreatment = 0;
		DoNotResetAnalTreatTermsAfterProp = 0;
	}
};

//-------------------------------------------------------------------------

struct srTDataMD {
	char* pData;
	char DataType[2]; // 'f'|'d'|'cf'|'cd'
	long AmOfDims;
	long DimSizes[10];
	double DimStartValues[10];
	double DimSteps[10];
	char DimScales[10][4]; //OC040912
	char DimUnits[10][255];
	char DataUnits[255];
	char DataName[255];
	int hState; //auxiliary
};

//-------------------------------------------------------------------------

struct srTMagFldHarm {
	int n;			//Harmonic no.
	char XorZ;		//'z' (vertical) or 'x' (horizontal)
	double K;		//Deflection parameter
	double Phase;	//Inititial phase [r]
};

//-------------------------------------------------------------------------

struct srTSRWRadInData;

//-------------------------------------------------------------------------

//struct srTBeamStatMom1 {
//	double Energy;	//[GeV]
//	double x;		//Horizontal position [m]
//	double dxds;	//Horizontal angle [r]
//	double z;		//Vertical position [m]
//	double dzds;	//Vertical angle [r]
//};

//-------------------------------------------------------------------------

//struct srTBeamStatMom2 {
//	double Mee;		//Squared relative RMS energy spread [r^2]
//	double Mxx;		//Horizontal squared RMS size [m^2]
//	double Mxxp;	//Horizontal mixed second-order central moment [m]
//	double Mxpxp;	//Horizontal squared RMS angular divergence [r^2]
//	double Mzz;		//Vertical squared RMS size [m^2]
//	double Mzzp;	//Vertical mixed second-order central moment [m]
//	double Mzpzp;	//Vertical squared RMS angular divergence [r^2]
//	double Mxz;		//[m^2]
//	double Mxpz;	//[m]
//	double Mxzp;	//[m]
//	double Mxpzp;	//[r^2]
//	double Mss;		//Longitudinal squared RMS size [m^2]
//};

//-------------------------------------------------------------------------

//struct srTElectronBeam {
//	double I;					//Current [A]
//	double s;					//Longitudinal position where stat. moments are defined [m]
//	srTBeamStatMom1 StatMom1;	//First-order statistical moments
//	srTBeamStatMom2 StatMom2;	//Second-order statistical moments
//};

//-------------------------------------------------------------------------

//struct srTMagFieldPer {
//	double Period;				//[m]
//	double Length;				//[m]
//	char Type;					//'c'- conventional, 't'- tapered, 'k'- optical klystron, 'i'- infinite
//	int NumHarm;				//Number of Harmonics
//	srTMagFieldHarm* HarmArr;	//Harmonics array
//};

//-------------------------------------------------------------------------

/* OK = 0 */
//#define OK 0

 	/**.  A teste function which use the database and smart pointer template classes. 
 	@param		No input parameters 
 	@return	An integer  value ( 0 : No Erro,  > 0 : Error Number,  < 0 : Warning Number) 
 	@author	P. Elleaume 
 	@version	1.0 
 	@see		... */ 
//EXP int CALL teste();

//-------------------------------------------------------------------------

EXP void CALL srUtiProgrIndFuncSet(int (*pExtFunc)(double CurVal));
EXP void CALL srUtiWfrExtModifFuncSet(int (*pExtFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp));
EXP void CALL srUtiOptElemGetInfByNameFuncSet(int (*pExtFunc)(const char* sNameOptElem, char** pDescrStr, int* LenDescr, void*));

//-------------------------------------------------------------------------

/** Specifies current API version number.  
@param [out] VerNoStr string specifying current version number of the API
@return	integer error code
@see		... */
EXP int CALL srUtiVerNo(char* VerNoStr);

EXP void CALL srUtiWarnTextGet(char* t, int WarnNo);

/** Creates Electron Beam structure.  
@param [out] i integer reference number of the structure created
@param [in] I current [A]
@param [in] Neb number of electrons in bunch
@param [in] pMom1 first-order moments
@param [in] nMom1 number of first-order central moments
@param [in] pMom2 second-order central moments
@param [in] nMom2 number of second-order central moments
@param [in] s longitudinal position [m]
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srElecBeamSet(int* i, double I, double Neb, double* pMom1, int nMom1, double* pMom2, int nMom2, double s);

EXP int CALL srElecBeamContSet(int* i, int* pElecElem, int nElecElem);

EXP int CALL srElecBeamPropagate(int iElecBeam, int iOptElem);

/** Converts Twiss parameters to second-order central Moments.  
@param [out] pMom2 second-order central moments
@param [in] nMom2 number of second-order central moments
@param [in] pTwissHor horizontal Twiss parameters
@param [in] pTwissVer vertical Twiss parameters
@param [in] RmsEnSpr RMS energy spread
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srElecBeamTwiss2Mom(double* pMom2, int nMom2, double* pTwissHor, double* pTwissVer, double RmsEnSpr);

EXP int CALL srElecBeamMom2EmitAndTwiss(double* pEmit, double* pBeta, double* pAlpha, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double Eta, double EtaPrime);
EXP int CALL srElecBeamMomAndEmit2Twiss(double* pBeta, double* pAlpha, double* pEta, double* pEmit, double Sigma2, double SigmaPrime2, double MixMom, double SigmaE2, double EtaPrime);

EXP int CALL srElecBeamGet(double* pI, double* pMom1, double* pMom2, double* ps0, int iElecBeam);

/** Creates Periodic Magnetic Field Harmonic structure.  
@param [out] i integer reference number of the structure created
@param [in] n harmonic number
@param [in] 'v' (vertical) or 'h' (horizontal)
@param [in] K deflection parameter
@param [in] Ph inititial phase [r]
@return	integer error code
@version	1.0 
@see		... */
//EXP int CALL srMagFldPerHarmSet(int* i, int n, char v_or_h, double K, double Ph);

/** Creates Periodic Magnetic Field structure.  
@param [out] i integer reference number of the structure created
@param [in] Per period [m]
@param [in] L length [m]
@param [in] pHarm array of harmonics ref. numbers
@param [in] nHarm number of harmonics (length of harmonics array)
@param [in] Type 'c'- conventional, 't'- tapered, 'k'- optical klystron, 'i'- infinite
@param [in] SpecPar special parameter related to particular Type of periodic mag. field structure
@return	integer error code
@version	1.0 
@see		... */
//EXP int CALL srMagFldPerSet(int* i, double Per, double L, double s0, int* pHarm, int nHarm, char Type, double SpecPar);
EXP int CALL srMagFldPerSet(int* i, double Per, double L, double s0, int nHarm, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase, char Type, double SpecPar);

/** Attempts to create Periodic Magnetic Field structure from Transversely Uniform Magnetic Field structure.    
@param [out] i integer reference number of the structure created
@param [in] iMagFld integer reference number of the Transversely Uniform Magnetic Field structure
@param [in] RelPrec relative precision (0 < RelPrec < 1)
@param [in] MaxHarm maximal number of field harmonics to create
@param [in] MaxPerLen maximal magnetic field period length to consider
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srMagFldPerSetFromTrUnif(int* i, int iMagFld, double RelPrec, int MaxHarm, double MaxPerLen);

EXP int CALL srMagFldPerGet(int iMagFldPer, double* pPer, double* pL, double* psCen, int* pnHarm, char* pType, double* pSpecPar, double* pkeVperGeV2);

EXP int CALL srMagFldPerHarmGet(int iMagFldPer, int* ArrHarmNo, char* ArrXorZ, double* ArrK, double* ArrPhase);

/** Creates Constant Magnetic Field structure.  
@param [out] i integer reference number of the structure created
@param [in] BH horizontal magnetic field [T]
@param [in] BV vertical magnetic field [T]
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srMagFldConstSet(int* i, double BH, double BV);

/** Creates Transversely Uniform Magnetic Field structure.  
@param [out] i integer reference number of the structure created
@param [in] sStart starting longitudinal position [m]
@param [in] sStep longitudinal step [m]
@param [in] np number of points
@param [in] pBh horizontal magnetic field array [T]
@param [in] pBv vertical magnetic field array [T]
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srMagFldTrUnifSet(int* i, double sStart, double sStep, int np, double* pBh, double* pBv);

/** Creates 3D Magnetic Field structure.  
@param [out] i integer reference number of the structure created
@param [in] xStart starting horizontal position [m]
@param [in] xStep step of horizontal position [m]
@param [in] nx number of points vs horizontal position
@param [in] yStart starting longitudinal position [m]
@param [in] yStep longitudinal step [m]
@param [in] ny number of points vs longitudinal position
@param [in] zStart starting vertical position [m]
@param [in] zStep step of vertical position [m]
@param [in] nz number of points vs vertical position
@param [in] pBx horizontal magnetic field array [T]
@param [in] pBz vertical magnetic field array [T]
@param [in] pBs longitudinal magnetic field array [T]
@return	integer error code
@version	1.0
@see		... */
EXP int CALL srMagFld3dSet(int* i, double xStart, double xStep, int nx, double yStart, double yStep, int ny, double zStart, double zStep, int nz, double* pBx, double* pBy, double* pBz);

EXP int CALL srMagFldContSet(int* i, int* pMagElem, int nMagElem);

/** Creates Wavefront sampling structure.  
@param [out] i integer reference number of the structure created
@param [in] s longitudinal position [m]
@param [in] hSt initial horizontal position [m]
@param [in] hFi final horizontal position [m]
@param [in] hN horizontal number of points
@param [in] vSt initial vertical position [m]
@param [in] vFi final vertical position [m]
@param [in] vN vertical number of points
@param [in] eSt initial photon energy [eV] or [nm]
@param [in] eFi final photon energy [eV] or [nm]
@param [in] eN number of photon energy points
@param [in] eUnit photon energy units ([eV] or [nm])
@param [in] tSt initial time moment [s]
@param [in] tFi final time moment [s]
@param [in] tN number of time moments
@return	integer error code
@version	1.0 
@see		... */
//EXP int CALL srWfrSmpSet(int* i, double s, double hSt, double hFi, int hN, double vSt, double vFi, int vN, double eSt, double eFi, int eN, char* eUnit);
EXP int CALL srWfrSmpSet(int* i, double s, double hSt, double hFi, int hN, double vSt, double vFi, int vN, double* pSurfData, double eSt, double eFi, int eN, char* PhotEnUnit, double tSt, double tFi, int tN, int presT =0, double* horOrtObsPlane =0, double* inNormObsPlane =0);

EXP int CALL srWfrSet(int* i, srTSRWRadInData* pSRWRadInData, int CopyData);

/** Creates Gaussian Beam structure.  
@param [out] i integer reference number of the structure created
@param [in] SpecFlux peak spectral flux [photons/s/0.1%bw]
@param [in] Polar polarization (1- Linear Hor., 2- Linear Vert., 3- Linear 45 deg., 4- Linear 135 deg., 5- Circular Right, 6- Circular Left)
@param [in] SigX horizontal RMS size [m]
@param [in] mx horizontal mode order
@param [in] SigZ vertical RMS size [m]
@param [in] mz vertical mode order
@param [in] SigT RMS pulse duration [s]
@param [in] TypeT pulse form-factor vs time (1- Gaussian, 2- Rectangular, 3- Triangular)
@param [in] pMom1 first-order statistical moments (array of 4 elements: hor. position [m], hor. angle [rad], vert. position [m], vert. angle [rad])
@param [in] s0 longitudinal position [m]
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srGausBeamSet(int* i, double SpecFlux, int Polar, double SigX, int mx, double SigZ, int mz, double SigT, int TypeT, double* pMom1, double s0, double RepRate, double PulseEn, double AvgPhotEn);

/** Sets up an optical element (introduced for compatibility with SRW for Igor)  
@param [out] i integer reference number of the optical element created
@param [in/out] pDescrStr array of strings describing the optical element; ATTENTION may be updated on return
@param [in/out] pLenDescr length of the array of strings describing the optical element; ATTENTION may be updated on return
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srOptElemSet(int* i, char** pDescrStr, int* pLenDescr, srTDataMD* pExtraData, int iWfr);

//EXP int CALL srOptElemGetInf(char** pDescrStr, int* LenDescr, int iOptElem);

EXP int CALL srOptDriftSet(int* piOptElem, double Length);

EXP int CALL srOptContSet(int* piOptElem, int* piMembArr, int nMemb);

/** Computes Stokes Parameters of Synchrotron Radiation.  
@param [out] i integer reference number of the Stokes structure created
@param [in] iElecBeam reference of Electron Beam structure
@param [in] iMagFld reference of Magnetic Field structure
@param [in] iRadSmp reference of Radiation Sampling structure
@param [in] pPrcPar pointer to precision parameters structure
@param [in] nPrcPar length of array of precision parameters
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srStokesComp(int* i, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar);
EXP int CALL srStokesCompExt(void* pStokesInData, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar);

/** Computes Electric Field of Synchrotron Radiation.  
@param [out] i integer reference number or Wavefront structure created
@param [in] iElecBeam reference of Electron Beam structure
@param [in] iMagFld reference of Magnetic Field structure
@param [in] iWfrSmp reference of Radiation Sampling structure
@param [in] pPrcPar precision parameters
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srElecFldComp(int* i, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar);
EXP int CALL srElecFldCompExt(void* pvWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar);

/** Computes Electric Field of Coherent Synchrotron Radiation.  
@param [out] i integer reference number or Wavefront structure created
@param [in] iElecBeam reference of Electron Beam structure
@param [in] iMagFld reference of Magnetic Field structure
@param [in] iWfrSmp reference of Radiation Sampling structure
@param [in] pPrcPar precision parameters
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL csrElecFldCompExt(void* pvWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrcPar);

/** Computes Electric Field of a Gaussian Beam.  
@param [out] i integer reference number of Wavefront structure created
@param [in] iGsnBeam reference of Gaussian Beam structure
@param [in] iMagFld reference of Magnetic Field structure
@param [in] iWfrSmp reference of Radiation Sampling structure
@param [in] pPrecPar precision parameters structure (has impact on the accuracy of wavefront propagation)
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srElecFldGausBeamComp(int* i, int iGsnBeam, int iWfrSmp, void* pPrecPar);
EXP int CALL srElecFldGausBeamCompExt(void* pvWfr, int iGsnBeam, int iWfrSmp, void* pVoidPrecElecFldGaus);

EXP int CALL srPowDensComp(void* pvWfr, int iElecBeam, int iMagFld, int iWfrSmp, void* pPrecPowDens);
EXP int CALL srPowDensCompExt(void* pvPow, int iElecBeam, int iMagFld, int iWfrSmp, void* pvPrecPowDens);

/** Computes Electron Trajectory.  
@param [out] pOutBtxData horizontal angle array data
@param [out] pOutXData horizontal position array data
@param [out] pOutBtyData dy/ds array data
@param [out] pOutYData y array data
@param [out] pOutBtzData vertical angle array data
@param [out] pOutZData vertical position array data
@param [in] ns length of trajectory arrays (number of point vs longitudinal position)
@param [in] sStart initial longitudinal position
@param [in] sStep step of longitudinal position
@param [in] iElecBeam reference of Electron Beam structure
@param [in] iMagFld reference of Magnetic Field structure
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srElecTrjComp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep, int iElecBeam, int iMagFld);

EXP int CALL srWfrComponGet(char* pData, int iWfr, int PolarizID, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z);
EXP int CALL srWfrComponGetExt(char* pData, void* pVoidRadInData, int PolarizCompon, int Int_or_Phase, int SectID, int TransvPres, double e, double x, double z);

EXP int CALL srWfrComponInteg(srTDataMD* pIntensOrigData, srTDataMD* pIntegParData, srTDataMD* pIntegResData);

EXP int CALL srWfrReprSet(int iWfr, char cReprID);
EXP int CALL srWfrReprSetExt(void* pVoidRadInData, char cReprID);

EXP int CALL srWfrSmpGet(double* pHorCenRange, int* nHor, double* pVerCenRange, int* nVer, double* pEnStartEnd, int* nEn, int iRad);

/** Computes definite integrals of (K3/2)^2*x^2 (1) and x*Int(K5/3(x1)) (2) (for SR calculations).
@param [out] ResInt the value computed
@param [in] type calculate (K3/2)^2*x^2 (1) or x*Int(K5/3(x1)) (2)
@param [in] x1 left integration boundary
@param [in] x2 right integration boundary
@param [in] prec precision parameter
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srUtiFuncIntKnXn(double* ResInt, int type, double x1, double x2, double prec);

/** Provides random numbers with n-d Gaussian distribution, n<=6.
@param [out] arGsnND output random numbers
@param [in] arXc centers of the n-d Gaussian distribution
@param [in] arSigma RMS sizes of the n-d Gaussian distribution
@param [in] n number of dimensions (n <= 6)
@param [in] init initialization key (if(init), initialized the random generator)
@param [in] rand_mode type of the random generator (0-standard, 1-LPtau)
@return	integer error code
@version	1.0 
@see		... */
EXP int CALL srUtiRandGsnMD(double* arGsnND, double* arXc, double* arSigma, int n, char init, char rand_mode);

EXP int CALL srUtiTrfOptMirSet(double* pM, double* pV, double* pAngles, double* pP0);

EXP int CALL srUtiConvolComp(srTDataMD* pSrcData, srTDataMD* DstData);

EXP int CALL srDel(int iObj);

/** Extracts transmission (/optical path difference) characteristic of a thin optical element
@param [out] pData pointer to  data buffer where the transmission characteristic should be stored
@param [in] iOptElem reference number of an optical element
@param [in] CharType type of transmission characteristic to extract (1- amplitude transmission, 2- intensity transmissions, 3- opt. path difference)
@param [in] xc horizontal center position [m] for the transmission data
@param [in] xr horizontal range [m]
@param [in] nx number of points vs horizontal position
@param [in] zc vertical center position [m] for the transmission data
@param [in] zr vertical range [m]
@param [in] nz number of points vs vertical position
@return integer error code
@version 1.0
@see ... */
EXP int CALL srOptTransmCharGet(float* pData, int iOptElem, int CharType, double xc, double xr, int nx, double zc, double zr, int nz);

#ifdef __cplusplus  
}
#endif

#endif

