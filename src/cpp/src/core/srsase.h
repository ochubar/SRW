/************************************************************************//**
 * File: srsase.h
 * Description: Wrapper class for calculation of SASE using GENESIS (F2C-ed), with input and output data compatible in SRW formats (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRSASE_H
//#define __SRSTOWIG_H
#define __SRSASE_H

//#ifdef __IGOR_PRO__
#ifndef __SRSEND_H
#include "srsend.h"
#endif
//#endif

#include "srmagfld.h"
#include "srprgind.h"

//*************************************************************************

typedef long int f2c_integer;
typedef float f2c_real;
typedef double f2c_doublereal;
typedef struct { f2c_real r, i; } f2c_complex;
typedef struct { f2c_doublereal r, i; } f2c_doublecomplex;
typedef short f2c_ftnlen;

//*************************************************************************

class srTSASE {

	static double EEV;// = 510999.06E0; //Energy units (mc^2) in eV
	static double VACIMP;// = 376.73E0; //Vacuum impedence in Ohms

	static double PlankConstReduced;// = 1.0545887E-34; //Reduced Plank constant in J*s
	static double SpeedOfLight;// = 2.99792458E08; //in m/s

	static long NPMAX;// = 65536; //# of particles
	static long NDMAX;// = 1250000; //maximum of particle in imported distribution
	static long NHMAX;// = 7; //# of particles
    static long NCMAX;// = 513; //# of gridpoints of cartesian mesh
    static long NRGRID;// = 100; //Number of radial points for s. c.
	static long NZMAX;// = 10000; //# of integration steps
	static long NSMAX;// = 50000; //# of slices

    static double PI;// = 376.73E0; //Vacuum impedence in Ohms
    static double TWOPI;// = 2.E0*PI
/*
     +          NZMAX  = 2501,                  !# of integration steps
     +          NSMAX  = 3300,                  !# of slices
     +          PIHALF = PI/2.D0,               !Pi/2
     +          NTMAX  = 160*160*300)           !dimension of field array 
*/
	bool m_ElecDistribShouldBeUsed;
	double m_dwSeedRad; //difference in cylcic frequencies between seed rad. and FEL resonant freq.

public:

	srTSend* pSend;
	//DLL_IMPLEMENT

	srTCompProgressIndicator CompProgressInd;

	srTEbmDat EbmDat;
	srTMagGroup MagDat;
	srTRadSASE InRad;
	srTSRWRadStructAccessData SeedRad;
	srTWfrSmp DistrInfoDat;
	srTPrecSASE PrecDat;
	srTControlAccessSASE ControlSASE;


	srTSASE() 
	{
		ZeroPtrs_tbunchcom();
		ZeroPtrs_cartcom();
		ZeroPtrs_beamcom();
		ZeroPtrs_wigcom();
		ZeroPtrs_diagcom();
		ZeroPtrs_workspace();
		ZeroPtrs_tslipcom();

		//pSend = 0;
		m_ElecDistribShouldBeUsed = false;
	}

	int InitGenesisStructs_OLD();
	int InitMainGenesisStructs();
	//int InitMainGenesisStructs(int numHarm);

	void ReleaseGenesisStructs();
	int CheckInputConsistency();
	int ConvertInputDataToGenesisFormat(int numHarm);
	int SetupOutputControlStruct(int numHarm);
	int AuxvalGenesisPlus();
	int DiagnoGenesisPlus(long PassCount, long istepz, long islice, int numHarm);

	//int GenesisCompDrive(srTSRWRadStructAccessData& RadAccessData);
	int GenesisCompDrive(srTSRWRadStructAccessData *arRadAccessData, int numHarm);

	int OutResDistrib(long _istepz, long _islice, double _xkw0);
	int OutElecDistrib(long _istepz, long _islice);
    int OutDumpResDistrib(int _islice);
	int EstimateGenesisNSLP();

	//int PrepSRWRadStructTD(srTSRWRadStructAccessData&);
	int PrepSRWRadStructTD(srTSRWRadStructAccessData* arRadAccessData, int numHarm);
	//int CopyRadSliceToSRWRadStructTD(int iSlice, srTSRWRadStructAccessData& RadAccessData);
	int CopyRadSliceToSRWRadStructTD(int iSlice, srTSRWRadStructAccessData* arRadAccessData, int numHarm);
	int FillInSRWRadStruct(srTSRWRadStructAccessData&);
	int InitAuxWfrParams(srTSRWRadStructAccessData *arSRWRadAccessData, int numHarm);
	//int PropagateWavefrontToObservationPlane(srTSRWRadStructAccessData&);
	int PropagateWavefrontToObservationPlane(srTSRWRadStructAccessData *arSRWRadAccessData, int numHarm);
	int ResizeForOtimalSamplingInSpot(srTSRWRadStructAccessData&);
	int ResizeForSpecifiedSampling(srTSRWRadStructAccessData&);
	int ResizeAndPropagateFromWaist(srTSRWRadStructAccessData&);
	int CheckCorrectParamAndResize(srTSRWRadStructAccessData&, srTRadResize&, char);
	int EstimHarmNumFromPhotEn(srTWigComSASE& WigCom);
	//void CalcBasicRadParamsFromTab(srTSRWRadStructAccessData& rad, double& prad0, double& zwaist, double& zrayl);

	double RadMeshRange();
	int UpdateInterface(long);
	double RadFieldMultip();
	double TDElecFieldConvConstSRW2GENESIS();

	int CreateWavefrontElField(srTSRWRadStructAccessData&);
	int CreateWavefrontElField(srTSRWRadStructAccessData*, int);

	int loadbeam_srw(f2c_integer*, f2c_doublereal*);
	int readpart_srw(f2c_integer*);
	int initrun_srw();
	int loadslpfld_srw(f2c_integer* nslp);
	int readfield_srw(f2c_doublecomplex* cin, f2c_integer* irec);
	int loadrad_srw(f2c_integer* islice);
	
	int Alloc_tbunchcom();
	int Alloc_cartcom();
	//int Alloc_cartcom(int numHarm);
	int Alloc_beamcom();
	int Alloc_wigcom();
	int Alloc_diagcom();
	int Alloc_workspace();
	//int Alloc_tslipcom();
	int Alloc_tslipcom(int nWfr, int nSlices);
	
	void Free_tbunchcom();
	void Free_cartcom();
	void Free_beamcom();
	void Free_wigcom();
	void Free_diagcom();
	void Free_workspace();
	void Free_tslipcom();

	void ZeroPtrs_tbunchcom();
	void ZeroPtrs_cartcom();
	void ZeroPtrs_beamcom();
	void ZeroPtrs_wigcom();
	void ZeroPtrs_diagcom();
	void ZeroPtrs_workspace();
	void ZeroPtrs_tslipcom();

	int RoundDoubleToInt(double x)
	{
		int ix = int(x);
		if((x - double(ix)) > 0.5) ix++;
		return ix;
	}
	//void ZeroArr(double *arr, long n)
	void ZeroArr(double *arr, long long n)
	{
		if((arr == 0) || (n <= 0)) return;
		double *t_arr = arr;
		//for(long i=0; i<n; i++) *(t_arr++) = 0;
		for(long long i=0; i<n; i++) *(t_arr++) = 0;
	}
	//void ZeroArr(f2c_integer *arr, long n)
	void ZeroArr(f2c_integer *arr, long long n)
	{
		if((arr == 0) || (n <= 0)) return;
		f2c_integer *t_arr = arr;
		//for(long i=0; i<n; i++) *(t_arr++) = 0;
		for(long long i=0; i<n; i++) *(t_arr++) = 0;
	}
	//void ZeroArr(f2c_doublecomplex *arr, long n)
	void ZeroArr(f2c_doublecomplex *arr, long long n)
	{
		if((arr == 0) || (n <= 0)) return;
		f2c_doublecomplex *t_arr = arr;
		//for(long i=0; i<n; i++) 
		for(long long i=0; i<n; i++) 
		{
			t_arr->r = 0; (t_arr++)->i = 0;
		}
	}
	template<class T> static long FindMaxArrElemInd(T* arToSearch, long len_arToSearch, T& maxElem)
	{
		if((arToSearch == 0) || (len_arToSearch <= 0)) return -1;
		T *t_arToSearch = arToSearch + 1;
		maxElem = arToSearch[0];
		long indMaxElem = 0;
		for(long i=1; i<len_arToSearch; i++)
		{
			if(maxElem < *t_arToSearch) 
			{
				maxElem = *t_arToSearch;
				indMaxElem = i;
			}
			t_arToSearch++;
		}
		return indMaxElem;
	}

	bool CheckIfElecDistrShouldBeUsed()
	{//true if there is any non-zero element in EbmDat.pElecDistr
	 //AND if PrecDat.UseElecDistr != 0
		if((EbmDat.pElecDistr == 0) || (EbmDat.nTotMacroPart <= 0) || (PrecDat.UseElecDistr == 0)) return false;

		float *tElecDistr = EbmDat.pElecDistr;
		for(long i=0; i<EbmDat.nTotMacroPart; i++)
		{
			if(*(tElecDistr++) != 0) return true;
		}
		return false;
	}

	double InterpBilin3D(double rx100, double rx010, double rx001, 
		double f000, double f100, double f010, double f001,
		double f110, double f101, double f011, double f111)
	{
		double a000 = f000;
		double a100 = f100 - f000;
		double a010 = f010 - f000;
		double a001 = f001 - f000;
		double a110 = f110 - f100 - f010 + f000;
		double a101 = f101 - f100 - f001 + f000;
		double a011 = f011 - f010 - f001 + f000;
		double a111 = f111 - f110 - f101 - f011 + f100 + f010 + f001 - f000;
		double res = a000 + rx100*(a100 + rx010*(a110 + rx001*a111) + rx001*a101) + rx010*(a010 + rx001*a011) + rx001*a001;
		return res;
	}
};

//*************************************************************************

#endif
