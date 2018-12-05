/************************************************************************//**
 * File: srradstr.h
 * Description: Auxiliary structures for various SR calculation methods (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRRADSTR_H
#define __SRRADSTR_H

#include <math.h>

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

#include "srobject.h"

//*************************************************************************

struct srTRadResize;
struct srTMinMaxEParam;
class srTWfrSmp;
struct srTMomentsPtrs;
class srTEbmDat;
class srTMagElem;
class srTTrjDat;
class srTGsnBeam;
struct srTSRWRadInData;
struct srTSRWRadStructWaveNames;
struct srTSRWRadStructWaveKeys;
struct srTEXZ;
struct srTParPrecElecFld;
struct srTDataPtrsForWfrEdgeCorr;

struct SRWLStructWaveFront;
typedef struct SRWLStructWaveFront SRWLWfr;
struct SRWLStructParticleBeam;
typedef struct SRWLStructParticleBeam SRWLPartBeam;
struct SRWLStructRadMesh;
typedef struct SRWLStructRadMesh SRWLRadMesh;

//*************************************************************************

class srTSRWRadStructAccessData : public CGenObject {

	void *m_pExtWfr; //OC130311: pointer to external wavefront structure (for eventual call-backs)

public:

	bool BaseRadWasEmulated;
	float *pBaseRadX, *pBaseRadZ;
	float *pBaseRadXaux, *pBaseRadZaux; //OC151115
	waveHndl wRad, wRadX, wRadZ;
	int hStateRadX, hStateRadZ;
	double eStep, eStart, xStep, xStart, zStep, zStart;
	long ne, nx, nz;

	double xStartTr, zStartTr;
	bool UseStartTrToShiftAtChangingRepresToCoord;
	//double tStartTr; //OC091115 //to check!
	//bool UseStartTrToShiftAtChangingRepresToTime; //OC091115 

	double RobsX, RobsZ; //these values should be treated as distances to waists
	double RobsXAbsErr, RobsZAbsErr;
	double xc, zc;
	double xWfrMin, xWfrMax, zWfrMin, zWfrMax; // Exact borders of the Wavefront
	char WfrEdgeCorrShouldBeDone; // To switch off/on manually
	double avgPhotEn; //averarage photon energy for time-domain simulations
	double avgT; //OC101115 //averarage tiem (auxiliary value)

	double UnderSamplingX, UnderSamplingZ;
	char AllowAutoSwitchToPropInUnderSamplingMode;
	double InvUnderSamplingThreshold;

	bool ResAfterWasEmulated;
	bool DoNotResizeAfter;
	srTRadResize* pResAfter;

	char Pres; // 0- Coord, 1- Ang.
	char PresT; // 0- Frequency (Photon Energy), 1- Time Domain (i.e. ne, eStep, eStart contain time parameters)
	char LengthUnit; // 0- m; 1- mm; 
	char PhotEnergyUnit; // 0- eV; 1- keV; 
	char ElecFldUnit; // 0- Arb. Units, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)
	//OC20112017
	char ElecFldAngUnit; //Electric field units in angular representation: 0- sqrt(Wavelength[m]*Phot/s/0.1%bw/mrad^2) vs rad/Wavelength[m], 1- sqrt(Phot/s/0.1%bw/mrad^2) vs rad; [Phot/s/0.1%bw] can be replaced by [J/eV] or [W], depending on ElecFldUnit, PresT and Pres

	bool WfrQuadTermCanBeTreatedAtResizeX; // is used at the time of one resize only
	bool WfrQuadTermCanBeTreatedAtResizeZ;
	double wfrReffX, wfrReffZ; //effective wavefront radii (to be used e.g. at resize) //OC150914

	char ElectronBeamEmulated; // 0 by def.
	DOUBLE *pElecBeam;
	waveHndl wElecBeam; // Can be this or Trajectory
	int hStateElecBeam;

	waveHndl wTrj; // Can be this or Electron Beam
	int hStateTrj;

	bool PropMatrWasEmulated;
	DOUBLE *p4x4PropMatr;
	waveHndl w4x4PropMatr;
	int hState4x4PropMatr;

	bool MomWereEmulated;
	//float *pMomX, *pMomZ;
	double *pMomX, *pMomZ; //OC130311
	waveHndl wMomX, wMomZ;
	int hStateMomX, hStateMomZ;
	bool MomWereCalcNum;

	bool WfrAuxDataWasEmulated;
	DOUBLE *pWfrAuxData;
	waveHndl wWfrAuxData;
	int hStateWfrAuxData;

	long AuxLong1, AuxLong2, AuxLong3, AuxLong4;
	double yStart; // to be used to keep longitudinal position of the initial wavefront

	bool m_xLinOnlyPhaseTermWasSubtracted, m_zLinOnlyPhaseTermWasSubtracted;
	//"xLinOnlyPhaseTermWasSubtracted = true" means that the true electric field is equal to pBaseRadX, pBaseRadZ data
	//multiplied by exp(i*k*dxcSub*x/RobsX)
	bool m_xQuadPhaseTermWasSubtracted, m_zQuadPhaseTermWasSubtracted; 
	//"xQuadPhaseTermWasSubtracted = true" means that the true electric field is equal to pBaseRadX, pBaseRadZ data
	//multiplied by exp(i*k*(x - xc)^2/(2*RobsX))
	double m_dxcSub, m_dzcSub;
	bool m_newExtWfrCreateNotAllowed; //OC130311

	srTSRWRadStructAccessData(srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	srTSRWRadStructAccessData(srTGsnBeam* pGsnBeam, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTEbmDat* pElecBeam, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, srTParPrecElecFld* pPrecElecFld);
	srTSRWRadStructAccessData(srTSRWRadInData* pRadInData, srTGsnBeam* pGsnBeam, srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	srTSRWRadStructAccessData(srTSRWRadInData* pSRWRadInData, int CopyData);
	srTSRWRadStructAccessData(srTSRWRadStructAccessData* pInRadStruct);
	srTSRWRadStructAccessData(const srTSRWRadStructAccessData&, bool createNewEmulStruct=true); //OC140411
	srTSRWRadStructAccessData(srTWfrSmp* pWfrSmp, bool AllocateData); //used for extracting transmission characteristics of optical elements
	srTSRWRadStructAccessData(SRWLWfr*, srTTrjDat* pTrj=0, double* arPrec=0);
	srTSRWRadStructAccessData(SRWLWfr*, srTGsnBeam*, double* arPrec=0);
	srTSRWRadStructAccessData(SRWLWfr*, double longPosSrc, double* arPrec=0); //Used for setting up spherical wave

	srTSRWRadStructAccessData()
	{
		Initialize();
	}

	~srTSRWRadStructAccessData()
	{
		DisposeEmulatedStructs();
	}

	void Initialize();
    void AuxSetupActions(srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, double NxNzOversamplingFactor);
    void AuxSetupActions2(srTSRWRadInData* pRadInData, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, double NxNzOversamplingFactor);
	void AuxSetupActions2SR(srTSRWRadInData* pRadInData, srTTrjDat* pTrjDat, srTWfrSmp* pWfrSmp, double From_s0ToObsPoint, srTParPrecElecFld* pPrecElecFld);
	//void AuxSetupActions2SR(SRWLWfr& srwlWfr, srTTrjDat& trjDat, double* precPar);
	void AuxSetupActionsArbSrc(SRWLWfr& srwlWfr, double Robs, double RobsAbsErr, double xElAtYsrc, double zElAtYsrc, double NxNzOversamplingFactor);
	void DisposeEmulatedStructs();
	void PreserveLogicsOfWfrLimitsAtRangeResizing(srTSRWRadStructAccessData* pOldRadData, char x_or_z);
	void FindMinMaxReE(srTMinMaxEParam& a);
	int EmulateElectronBeamStruct(srTEbmDat& EbmDat);
	int EmulateElectronBeamStruct(const SRWLPartBeam& srwlPartBeam);
	int EmulateElectronBeamStruct(srTGsnBeam& GsnBeam);
	int OutElectronBeamStruct(srTEbmDat& EbmDat);
	int OutElectronBeamStruct(SRWLPartBeam& srwlPartBeam);
	void EstimateAndSetUnderSampling();
	int ReAllocBaseRadAccordingToNeNxNz(char =0);
	int AllocBaseRadAccordingToNeNxNz(char =0);
	void DeAllocBaseRadAccordingToNeNxNz(char =0);
	void ZeroPtrs();
	void UpdateObsParam(srTWfrSmp& DistrInfoDat);
	void SetObsParamFromWfr(srTWfrSmp& smp);
	void SetupRadMomentsPtrs(srTMomentsPtrs& MomPtrsX, srTMomentsPtrs& MomPtrsZ);
	//srTMomentsPtrs OneSetOfMomentsPtrs(float* tMom);
	srTMomentsPtrs OneSetOfMomentsPtrs(double* tMom); //OC130311
	void SetRadSamplingFromObs(srTWfrSmp& DistrInfoDat);

	void ProcessNxNzForPropag(srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	//void ProcessNxNzForPropag(double NxNzOversamplingFactor, long& nx, long& nz);
    void CheckNxNzForSR(srTWfrSmp* pWfrSmp, double NxNzOversamplingFactor);
	void CheckNxNzForSR(double NxNzOversamplingFactor, long& _nx, long& _nz);

	void AllocElectronBeam();
	void Alloc4x4PropMatr();
    void AllocStatMom();
    void AllocWfrAux();

	void CopyBaseRadData(float* pInBaseRadX, float* pInBaseRadZ);
	//void CopyStatMomData(float* pInMomX, float* pInMomZ);
	void CopyStatMomData(double* pInMomX, double* pInMomZ); //OC130311
    void CopyElectronBeamData(DOUBLE* pInElecBeam);
    void Copy4x4PropMatrData(DOUBLE* pIn4x4PropMatr);
    void CopyWfrAuxData(DOUBLE* pInWfrAuxData);

	void InSRWRadPtrs(srTSRWRadInData*, bool DataShouldBeCopied = false);
	void InSRWRadPtrs(SRWLWfr&);
	void OutSRWRadPtrs(srTSRWRadInData*);
	void OutSRWRadPtrs(SRWLWfr&);

	//srTSRWRadInData* CreateCorrespSRWRadInData();
	//int ModifyWfrNeNxNz(char PolarizComp = 0);
	int ModifyWfrNeNxNz(char PolarizComp = 0, bool backupIsReq = false); //OC131115
	int DeleteWfrBackupData(char PolarizComp = 0); //OC151115
	int AllocExtIntArray(char type, char dep, char*& pcAlloc); //OC18082018

	int GetWfrStructNames(srTSRWRadStructWaveNames& RadStructNames);
	int CreateNewWfrStruct(srTSRWRadStructWaveNames& Names);
	int DeleteWfrStructWaves(srTSRWRadStructWaveKeys& RadKeys);
    int RenameWfrStruct(srTSRWRadStructWaveNames& RadStructNames);

	//int FindAverageDistanceToSource(srTTrjDat& TrjDat, srTWfrSmp& WfrSmp, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc);
	int FindAverageDistanceToSource(srTTrjDat& TrjDat, srTWfrSmp& WfrSmp, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc, srTParPrecElecFld* pPrecElecFld=0);
	int FindAverageDistanceToSource(srTTrjDat& TrjDat, double& Robs, double& RobsAbsErr, double& xElAtYsrc, double& zElAtYsrc, double* precPar);
	void AddStokesAtPoint(srTEXZ& EXZ, float* pStokesVal);
	
	void CheckAndSubtractPhaseTermsLin(double newXc, double newZc);
	void CheckAndResetPhaseTermsLin();
	void EstimateOversamplingFactors(double& estimOverSampX, double& estimOverSampZ);
	void MirrorFieldData(int sx, int sz);

	int SetupWfrEdgeCorrData(float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrsForWfrEdgeCorr);
	void MakeWfrEdgeCorrection(float* pDataEx, float* pDataEz, srTDataPtrsForWfrEdgeCorr& DataPtrs);
	int ExtractSliceConstEorT(long ie, float*& pOutEx, float*& pOutEz);
	int SetupSliceConstEorT(long ie, float* pInEx, float* pInEz);

	int ShiftWfrByInterpolVsXZ(double shiftX, double shiftY);
	void FlipFieldData(bool flipOverX, bool flipOverZ);
	void TransposeFieldData();

	int SetRepresCA(char CoordOrAng); //set Coordinate or Angular representation
	int SetRepresFT(char FreqOrTime); //set Frequency or Time representation
	int ComputeRadMoments();

	//void EstimWfrRadCen(double& resR, double& resCen, char cutX_or_Z, char fldX_or_Z=0, double relArgRange=0.2, double relArgCenOther=0.5);
	bool CheckIfQuadTermTreatIsBenefit(char cutX_or_Z, char fldX_or_Z=0);
	void GetIntMesh(char dep, SRWLRadMesh& mesh); //OC23082018

	void SetupSrwWfrAuxData()
	{
		if(pWfrAuxData == 0) return;

		RobsX = *pWfrAuxData;
		RobsZ = *(pWfrAuxData+1);
		RobsXAbsErr = *(pWfrAuxData+2);
		RobsZAbsErr = *(pWfrAuxData+3);
		xWfrMin = *(pWfrAuxData+4);
		xWfrMax = *(pWfrAuxData+5);
		zWfrMin = *(pWfrAuxData+6);
		zWfrMax = *(pWfrAuxData+7);
		UnderSamplingX = *(pWfrAuxData+8);
		UnderSamplingZ = *(pWfrAuxData+9);
		xc = *(pWfrAuxData+10);
		zc = *(pWfrAuxData+11);

		m_xLinOnlyPhaseTermWasSubtracted = (*(pWfrAuxData+12) == 0)? false : true;
		m_zLinOnlyPhaseTermWasSubtracted = (*(pWfrAuxData+13) == 0)? false : true;
		m_xQuadPhaseTermWasSubtracted = (*(pWfrAuxData+14) == 0)? false : true;
		m_zQuadPhaseTermWasSubtracted = (*(pWfrAuxData+15) == 0)? false : true;
		m_dxcSub = *(pWfrAuxData+16);
		m_dzcSub = *(pWfrAuxData+17);
	}
	void UpdateSrwWfrAuxData()
	{
		if(pWfrAuxData == 0) return;

		*pWfrAuxData = RobsX;
		*(pWfrAuxData+1) = RobsZ;
		*(pWfrAuxData+2) = RobsXAbsErr;
		*(pWfrAuxData+3) = RobsZAbsErr;
		*(pWfrAuxData+4) = xWfrMin;
		*(pWfrAuxData+5) = xWfrMax;
		*(pWfrAuxData+6) = zWfrMin;
		*(pWfrAuxData+7) = zWfrMax;
		*(pWfrAuxData+8) = UnderSamplingX;
		*(pWfrAuxData+9) = UnderSamplingZ;
		*(pWfrAuxData+10) =	xc;
		*(pWfrAuxData+11) =	zc;

		*(pWfrAuxData+12) =	m_xLinOnlyPhaseTermWasSubtracted? 1 : 0;
		*(pWfrAuxData+13) = m_zLinOnlyPhaseTermWasSubtracted? 1 : 0;
		*(pWfrAuxData+14) = m_xQuadPhaseTermWasSubtracted? 1 : 0;
		*(pWfrAuxData+15) = m_zQuadPhaseTermWasSubtracted? 1 : 0;
		*(pWfrAuxData+16) = m_dxcSub;
		*(pWfrAuxData+17) = m_dzcSub;
	}
	void SetupNonZeroWavefrontLimitsAtCreation()
	{// Is called ONLY at the Wfr Creation
		xWfrMin = xStart;
		xWfrMax = xStart + xStep*(nx - 1);
		zWfrMin = zStart;
		zWfrMax = zStart + zStep*(nz - 1);
	}
	void SetNonZeroWavefrontLimitsToFullRange()
	{// This switches off the sharp Wfr edges treatment
		xWfrMin = xStart;
		xWfrMax = xStart + nx*xStep;
		zWfrMin = zStart;
		zWfrMax = zStart + nz*zStep;
	}
	bool WfrXLimitsSeemDefined()
	{
		double xEnd = xStart + xStep*nx;
		if((xWfrMin >= (xStart - xStep)) && (xWfrMin < xWfrMax) && (xWfrMax <= (xEnd + xStep))) return true;
		return false;
	}
	bool WfrZLimitsSeemDefined()
	{
		double zEnd = zStart + zStep*nz;
		if((zWfrMin >= (zStart - zStep)) && (zWfrMin < zWfrMax) && (zWfrMax <= (zEnd + zStep))) return true;
		return false;
	}

	void SetupXcZcFromElecData()
	{//this is not used for wavefront, because it is not correct.
		if((pElecBeam == 0) || (p4x4PropMatr == 0)) return;

		double Xc = *(pElecBeam + 2);
		double Xdc = *(pElecBeam + 3);
		double Zc = *(pElecBeam + 4);
		double Zdc = *(pElecBeam + 5);

		DOUBLE *pMStr[4];
		for(int i=0; i<4; i++)
		{
			int i4 = i << 2;
			pMStr[i] = p4x4PropMatr + i4;
		}

		DOUBLE *pStr0 = pMStr[0], *pStr2 = pMStr[2];
		xc = pStr0[0]*Xc + pStr0[1]*Xdc + pStr0[2]*Zc + pStr0[3]*Zdc;
		zc = pStr2[0]*Xc + pStr2[1]*Xdc + pStr2[2]*Zc + pStr2[3]*Zdc;
	}

	DOUBLE* ElemPtrOf4x4PropMatr(int i, int j)
	{
		return p4x4PropMatr + j*4 + i;
	}
	void Setup4x4PropMatrPtrs(DOUBLE** PropMatrPtrs)
	{
		DOUBLE** tPropMatrPtrs = PropMatrPtrs;
		DOUBLE* t4x4PropMatr = p4x4PropMatr;
		for(int i=0; i<16; i++)
		{
			*(tPropMatrPtrs++) = t4x4PropMatr++;
		}
	}
	void InitialSetupOf4x4PropMatr(double Robs)
	{
		if(p4x4PropMatr == 0) return;

		int i;
		DOUBLE* t4x4PropMatr = p4x4PropMatr;
		for(i=0; i<16; i++) *(t4x4PropMatr++) = 0.;
		t4x4PropMatr = p4x4PropMatr;
		for(i=0; i<4; i++)
		{
			*t4x4PropMatr = 1.; t4x4PropMatr += 5;
		}
		*(p4x4PropMatr + 1) = Robs;
		*(p4x4PropMatr + 11) = Robs;
	}

	char ThereIsUnderSamplingX()
	{
		const double UnderSampThresh = 0.2; //0.4; // To steer
		return ((::fabs(UnderSamplingX - 1.) > UnderSampThresh) && (UnderSamplingX != 0))? 1 : 0;
	}
	char ThereIsUnderSamplingZ()
	{
		const double UnderSampThresh = 0.2; //0.4; // To steer
		return ((::fabs(UnderSamplingZ - 1.) > UnderSampThresh) && (UnderSamplingZ != 0))? 1 : 0;
	}

	bool CheckIfSamplingResolvesElecFieldWithQuadPhaseTerm()
	{
		double estimOverSampX, estimOverSampZ;
		EstimateOversamplingFactors(estimOverSampX, estimOverSampZ);
		if((estimOverSampX >= 0.95) && (estimOverSampZ >= 0.95)) return true;
		else return false;
	}

	void NormalizeElFieldToArbUnits()
	{
		if(ElecFldUnit != 0) return;
		if((pBaseRadX == 0) && (pBaseRadZ == 0)) return;
		float MaxAmplitude = DetermineMaxFieldAmplitude();
		if(MaxAmplitude <= 0.) return;

		MultiplyElFieldBy((float)(1./MaxAmplitude));
	}

	float DetermineMaxFieldAmplitude()
	{
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);
		if((!RadXisDefined) && (!RadZisDefined)) return 0.;

		float *tEx = pBaseRadX;
		float *tEz = pBaseRadZ;

		float MaxAmpXe2 = 0., MaxAmpZe2 = 0.;
		for(int iz=0; iz<nz; iz++)
		{
			for(int ix=0; ix<nx; ix++)
			{
				for(int ie=0; ie<ne; ie++)
				{
					if(RadXisDefined) 
					{
						float TestVal = (*tEx)*(*tEx) + (*(tEx+1))*(*(tEx+1));
						if(MaxAmpXe2 < TestVal) MaxAmpXe2 = TestVal;
						tEx += 2;
					}
					if(RadZisDefined) 
					{
						float TestVal = (*tEz)*(*tEz) + (*(tEz+1))*(*(tEz+1));
						if(MaxAmpZe2 < TestVal) MaxAmpZe2 = TestVal;
						tEz += 2;
					}
				}
			}
		}
		float MaxAmpE2 = MaxAmpXe2;
		if(MaxAmpE2 < MaxAmpZe2) MaxAmpE2 = MaxAmpZe2;
		return (float)sqrt(MaxAmpE2);
	}

	void MultiplyElFieldBy(float a)
	{
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);
		if((!RadXisDefined) && (!RadZisDefined)) return;

		float *tEx = pBaseRadX;
		float *tEz = pBaseRadZ;

		for(int iz=0; iz<nz; iz++)
		{
			for(int ix=0; ix<nx; ix++)
			{
				for(int ie=0; ie<ne; ie++)
				{
					if(RadXisDefined) 
					{
						*(tEx++) *= a; *(tEx++) *= a;
					}
					if(RadZisDefined) 
					{
						*(tEz++) *= a; *(tEz++) *= a;
					}
				}
			}
		}
	}

	void MultiplyElFieldByPhaseLin(double xMult, double zMult)
	{
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);
		if((!RadXisDefined) && (!RadZisDefined)) return;

		float *tEx = pBaseRadX;
		float *tEz = pBaseRadZ;
		
		double z = zStart;	
		for(int iz=0; iz<nz; iz++)
		{
			double dPhZ = zMult*z;
			double x = xStart;
			for(int ix=0; ix<nx; ix++)
			{
				double dPh = dPhZ + xMult*x;
				double cosPh = cos(dPh), sinPh = sin(dPh);
				
				for(int ie=0; ie<ne; ie++)
				{
					if(RadXisDefined) 
					{
						//*(tEx++) *= a; *(tEx++) *= a;
						double newReEx = (*tEx)*cosPh - (*(tEx + 1))*sinPh;
						double newImEx = (*tEx)*sinPh + (*(tEx + 1))*cosPh;
						*(tEx++) = (float)newReEx; *(tEx++) = (float)newImEx;
					}
					if(RadZisDefined) 
					{
						//*(tEz++) *= a; *(tEz++) *= a;
						double newReEz = (*tEz)*cosPh - (*(tEz + 1))*sinPh;
						double newImEz = (*tEz)*sinPh + (*(tEz + 1))*cosPh;
						*(tEz++) = (float)newReEz; *(tEz++) = (float)newImEz;
					}
				}
				x += xStep;
			}
			z += zStep;
		}
	}
	
	void AssignElFieldToConst(float E0_Re, float E0_Im)
	{
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);
		if((!RadXisDefined) && (!RadZisDefined)) return;

		if(RadXisDefined)
		{
			float *tEx = pBaseRadX;
			for(int iz=0; iz<nz; iz++)
			{
				for(int ix=0; ix<nx; ix++)
				{
					for(int ie=0; ie<ne; ie++)
					{
						*(tEx++) = E0_Re; *(tEx++) = E0_Im;
					}
				}
			}
		}
		if(RadZisDefined)
		{
            float *tEz = pBaseRadZ;
			for(int iz=0; iz<nz; iz++)
			{
				for(int ix=0; ix<nx; ix++)
				{
					for(int ie=0; ie<ne; ie++)
					{
                        *(tEz++) = E0_Re; *(tEz++) = E0_Im;
					}
				}
			}
		}
	}

	void ExtractElFieldAmplitude(int CharType, float* pData)
	{
		if(pData == 0) return;
		bool RadXisDefined = (pBaseRadX != 0);
		bool RadZisDefined = (pBaseRadZ != 0);
		if((!RadXisDefined) && (!RadZisDefined)) return;
		if((CharType <= 0) || (CharType > 2)) return;

		float *tEx = pBaseRadX,  *tEz = pBaseRadZ;
		float *tData = pData;
		for(int iz=0; iz<nz; iz++)
		{
			for(int ix=0; ix<nx; ix++)
			{
				for(int ie=0; ie<ne; ie++)
				{
					*tData = 0;
					if(RadXisDefined) 
					{ 
						*tData += (*tEx)*(*tEx); tEx++;
						*tData += (*tEx)*(*tEx); tEx++;
					}
					if(RadZisDefined) 
					{ 
						*tData += (*tEz)*(*tEz); tEz++;
						*tData += (*tEz)*(*tEz); tEz++;
					}
					if(CharType == 1) *tData = ::sqrt(*tData); //amplitude transmission
					tData++;
				}
			}
		}
	}

	void CalcBasicRadParamsFromTab(double& prad0, double& zwaist, double& sigRe2)
	{
		prad0 = zwaist = sigRe2 = 0;
		if(!ElecFieldIsDefined()) return;

		bool EXisDefined = pBaseRadX != 0;
		bool EZisDefined = pBaseRadZ != 0;

		double sumEnergy=0, sumSigTe2=0, sumSigRe2=0, sumT=0;
		float *tEX = pBaseRadX, *tEZ = pBaseRadZ;

		double z = zStart;
		for(int iz=0; iz<nz; iz++)
		{
			double zz = z*z;
			double x = xStart;
			for(int ix=0; ix<nx; ix++)
			{
				double xx = x*x;
				double xx_p_zz = xx + zz;
				double t = eStart;
				for(int ie=0; ie<ne; ie++)
				{
					if(EXisDefined) 
					{ 
						sumEnergy += (*tEX)*(*tEX); tEX++;
						sumEnergy += (*tEX)*(*tEX); tEX++;
					}
					if(EZisDefined) 
					{ 
						sumEnergy += (*tEZ)*(*tEZ); tEZ++;
						sumEnergy += (*tEZ)*(*tEZ); tEZ++;
					}
					sumSigTe2 += t*t*sumEnergy;
					sumT += t*sumEnergy;
					sumSigRe2 += xx_p_zz*sumEnergy;

					t += eStep;
				}
				x += xStep;
			}
			z += zStep;
		}
		double tRange = eStep*(ne - 1); 
		double xRange = xStep*(nx - 1);
		double zRange = zStep*(nz - 1);
		double tRange_xRange_zRange = tRange*xRange*zRange;
		double energyTot = (1E+06)*sumEnergy*tRange_xRange_zRange/(ne*nx*nz);

		double sigTe2 = sumSigTe2/sumEnergy;
		double sigT = sqrt(sigTe2);
		double max_sigT = tRange/2.35;
		if(sigT > max_sigT) sigT = max_sigT;

		//double avgT = sumT/sumEnergy; //OC020110
		sigRe2 = sumSigRe2/sumEnergy;
		double max_sigRe2 = sqrt(xRange*xRange + zRange*zRange)/2.35;
		if(sigRe2 > max_sigRe2) sigRe2 = max_sigRe2;

		const double inv_sqrt_two_pi = 1./sqrt(6.28318530717958);
		prad0 = 0.5*energyTot*inv_sqrt_two_pi/sigT; //average power in [W]
		zwaist = 0; //to implement
		//double W0 = InRad.WaistDiam; // or 0.5* ?
		//inputcom_1.zrayl = PI*W0*W0/inputcom_1.xlamds; // "ZRAYL"
	}

	static void CalcStokesFromE(float* tEx, float* tEz, float* Stokes)
	{
		float &EwX_Re = *tEx, &EwX_Im = *(tEx + 1);
		float &EwZ_Re = *tEz, &EwZ_Im = *(tEz + 1);
		double LinHor = EwX_Re*EwX_Re + EwX_Im*EwX_Im;
		double LinVer = EwZ_Re*EwZ_Re + EwZ_Im*EwZ_Im;
		Stokes[0] = (float)(LinHor + LinVer);
		Stokes[1] = (float)(LinHor - LinVer);
		Stokes[2] = (float)(-2.*(EwX_Re*EwZ_Re + EwX_Im*EwZ_Im));
		Stokes[3] = (float)(2.*(-EwX_Re*EwZ_Im + EwX_Im*EwZ_Re));
	}
	//Interp4pRel(tSt00++, tSt10++, tSt01++, tSt11++, xr, zr);

	void CosAndSin(double x, float& Cos, float& Sin)
	{
		const double TwoPI = 6.2831853071796;
		const double ThreePIdTwo = 4.7123889803847;
		const double One_dTwoPI = 0.1591549430919;
		const double HalfPI = 1.5707963267949;
		const double PI = 3.141592653590;
		const double a2c = -0.5, a4c = 0.041666666666667, a6c = -0.0013888888888889, a8c = 0.000024801587301587, a10c = -2.755731922E-07;
		const double a3s = -0.16666666666667, a5s = 0.0083333333333333, a7s = -0.0001984126984127, a9s = 2.755731922E-06, a11s = -2.505210839E-08;

		x -= TwoPI*int(x*One_dTwoPI);
		if(x < 0.) x += TwoPI;

		char ChangeSign=0;
		if(x > ThreePIdTwo) x -= TwoPI;
		else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

		double xe2 = x*x;
		Cos = float(1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c)))));
		Sin = float(x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s))))));
		if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
	}

	//void SetupExpCorrArray(float* pCmpData, long AmOfPt, double x, double qStart, double qStep)
	void SetupExpCorrArray(float* pCmpData, long long AmOfPt, double x, double qStart, double qStep)
	{
		const double TwoPi = 6.28318530717959;
		double TwoPiX = TwoPi*x;
		double q = qStart;
		float *tCmpData = pCmpData;
		//for(long i=0; i<AmOfPt; i++)
		for(long long i=0; i<AmOfPt; i++)
		{
			double Arg = TwoPiX*q;
			float Co, Si;
			CosAndSin(Arg, Co, Si);
			*(tCmpData++) = Co; *(tCmpData++) = -Si;
			q += qStep;
		}
	}

	void SetupRadXorZSectFromSliceConstEorT(float* pInEx, float* pInEz, long _nx, long _nz, char vsX_or_vsZ, long iSect, float* pOutEx, float* pOutEz)
	{
		long Per = (vsX_or_vsZ == 'x')? 2 : (_nx << 1);
		float *tOutEx = pOutEx, *tOutEz = pOutEz;
		long StartOffset = (vsX_or_vsZ == 'x')? iSect*(_nx << 1) : (iSect << 1);
		float *tEx = pInEx + StartOffset, *tEz = pInEz + StartOffset;
		long nPt = (vsX_or_vsZ == 'x')? _nx : _nz;

		for(int i=0; i<nPt; i++)
		{
			*(tOutEx++) = *tEx; *(tOutEx++) = *(tEx + 1);
			*(tOutEz++) = *tEz; *(tOutEz++) = *(tEz + 1);
			tEx += Per; tEz += Per;
		}
	}

	bool QuadPhaseTermCanBeTreated()
	{//Later treat X and Z fully separately here and at removing the corresponding terms from Phase !!!
	 //same as srTGenOptElem::WaveFrontTermCanBeTreated(srTSRWRadStructAccessData& RadAccessData)
		const double CritRatTransvLong = 0.1;
		const double CritRelRobsErr = 0.2; //0.1; //0.2;
		const double Pi = 3.14159265358979;

		char RobsXErrIsSmall = ::fabs(RobsXAbsErr) < CritRelRobsErr*(::fabs(RobsX));
		char RobsZErrIsSmall = ::fabs(RobsZAbsErr) < CritRelRobsErr*(::fabs(RobsZ));

		if(Pres == 0) // Coord
		{
			double xMagn = ::fabs((nx - 1)*xStep);
			double zMagn = ::fabs((nz - 1)*zStep);

			char AnglesXAreSmall = (xMagn < CritRatTransvLong*(::fabs(RobsX)));
			char AnglesZAreSmall = (zMagn < CritRatTransvLong*(::fabs(RobsZ)));

			WfrQuadTermCanBeTreatedAtResizeX = (AnglesXAreSmall && RobsXErrIsSmall);
			WfrQuadTermCanBeTreatedAtResizeZ = (AnglesZAreSmall && RobsZErrIsSmall);

			return (WfrQuadTermCanBeTreatedAtResizeX || WfrQuadTermCanBeTreatedAtResizeZ);
		}
		else // Ang
		{// Not sure about this...
			double ePh = eStart;
			double xRatMax = 1.E-23, zRatMax = 1.E-23, MaxPhaseChange = 1.E-23;
			for(int ie=0; ie<ne; ie++)
			{
				double Lambda_m = 1.239842e-06/ePh;
				double xTetMagn = ::fabs((nx - 1)*xStep*Lambda_m);
				double zTetMagn = ::fabs((nz - 1)*zStep*Lambda_m);
				double PhaseChange = ::fabs((Pi/Lambda_m)*(RobsX*xTetMagn*xTetMagn + RobsZ*zTetMagn*zTetMagn));

				if(xTetMagn > xRatMax) xRatMax = xTetMagn;
				if(zTetMagn > zRatMax) zRatMax = zTetMagn;
				if(PhaseChange > MaxPhaseChange) MaxPhaseChange = PhaseChange;

				ePh += eStep;
			}

			char AnglesXAreSmall = (xRatMax < CritRatTransvLong);
			char AnglesZAreSmall = (zRatMax < CritRatTransvLong);
			char PhaseChangeIsLarge = (MaxPhaseChange > 2.*Pi);

			WfrQuadTermCanBeTreatedAtResizeX = (AnglesXAreSmall && RobsXErrIsSmall);
			WfrQuadTermCanBeTreatedAtResizeZ = (AnglesZAreSmall && RobsZErrIsSmall);

			return ((WfrQuadTermCanBeTreatedAtResizeX || WfrQuadTermCanBeTreatedAtResizeZ) && PhaseChangeIsLarge);
		}
	}

	void TreatQuadPhaseTermTerm(char AddOrRem, char PolComp=0, int ieOnly=-1)
	{//same as srTGenOptElem::TreatStronglyOscillatingTerm(srTSRWRadStructAccessData& RadAccessData, char AddOrRem, char PolComp, int ieOnly)
		//Later treat X and Z coordinates separately here!!!

		char TreatPolCompX = ((PolComp == 0) || (PolComp == 'x')) && (pBaseRadX != 0);
		char TreatPolCompZ = ((PolComp == 0) || (PolComp == 'z')) && (pBaseRadZ != 0);

		const double Pi = 3.14159265358979;
		double Const = Pi*1.E+06/1.239854; // Assumes m and eV

		double ConstRx = (Pres == 0)? Const/RobsX : -Const*RobsX;
		double ConstRz = (Pres == 0)? Const/RobsZ : -Const*RobsZ;

		if(AddOrRem == 'r') { ConstRx = -ConstRx; ConstRz = -ConstRz;}

		double ConstRxE, ConstRzE;
		double ePh = eStart, x, z, zE2;
		double Phase;
		float CosPh, SinPh;

		float *pEX0 = 0, *pEZ0 = 0;
		if(TreatPolCompX) pEX0 = pBaseRadX;
		if(TreatPolCompZ) pEZ0 = pBaseRadZ;

		//long PerX = ne << 1;
		//long PerZ = PerX*nx;
		long long PerX = ne << 1;
		long long PerZ = PerX*nx;

		int ieStart=0, ieBefEnd=ne;
		if((ieOnly >= 0) && (ieOnly < ne))
		{
			ieStart = ieOnly; ieBefEnd = ieOnly + 1;
		}

		for(int ie=ieStart; ie<ieBefEnd; ie++)
		{
			if(PresT == 1)
			{
				ePh = avgPhotEn; //?? OC041108
			}

			//long Two_ie = ie << 1;
			long long Two_ie = ie << 1;

			ConstRxE = ConstRx*ePh;
			ConstRzE = ConstRz*ePh;

			if(Pres == 1)
			{
				double Lambda_m = 1.239842e-06/ePh;
				if(PhotEnergyUnit == 1) Lambda_m *= 0.001; // if keV

				double Lambda_me2 = Lambda_m*Lambda_m;
				ConstRxE *= Lambda_me2;
				ConstRzE *= Lambda_me2;
			}

			z = zStart - zc;

			zE2 = z*z;
			double PhaseAddZ = 0.;
			if(WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE*zE2;

			for(int iz=0; iz<nz; iz++)
			{
				//long izPerZ = iz*PerZ;
				long long izPerZ = iz*PerZ;
				float *pEX_StartForX = pEX0 + izPerZ;
				float *pEZ_StartForX = pEZ0 + izPerZ;

				x = xStart - xc;

				for(int ix=0; ix<nx; ix++)
				{
					//long ixPerX_p_Two_ie = ix*PerX + Two_ie;
					long long ixPerX_p_Two_ie = ix*PerX + Two_ie;

					//Phase = ConstRxE*x*x + ConstRzE*zE2;
					Phase = PhaseAddZ;
					if(WfrQuadTermCanBeTreatedAtResizeX) Phase += ConstRxE*x*x;

					CosAndSin(Phase, CosPh, SinPh);

					if(TreatPolCompX)
					{
						float *pExRe = pEX_StartForX + ixPerX_p_Two_ie;
						float *pExIm = pExRe + 1;
						double ExReNew = (*pExRe)*CosPh - (*pExIm)*SinPh;
						double ExImNew = (*pExRe)*SinPh + (*pExIm)*CosPh;
						*pExRe = (float)ExReNew; *pExIm = (float)ExImNew;
					}
					if(TreatPolCompZ)
					{
						float *pEzRe = pEZ_StartForX + ixPerX_p_Two_ie;
						float *pEzIm = pEzRe + 1;
						double EzReNew = (*pEzRe)*CosPh - (*pEzIm)*SinPh;
						double EzImNew = (*pEzRe)*SinPh + (*pEzIm)*CosPh;
						*pEzRe = (float)EzReNew; *pEzIm = (float)EzImNew;
					}

					x += xStep;
				}
				z += zStep;
				zE2 = z*z;
				PhaseAddZ = 0.;
				if(WfrQuadTermCanBeTreatedAtResizeZ) PhaseAddZ = ConstRzE*zE2;
			}
			ePh += eStep;
		}
	}

	void CopySymEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, bool ChangeSignEx, bool ChangeSignEz)
	{
		float *tOrigEx = pOrigDataEx, *tSymEx = pSymDataEx;
		float *tOrigEz = pOrigDataEz, *tSymEz = pSymDataEz;
		for(int ie = 0; ie < ne; ie++)
		{
			*tSymEx = *(tOrigEx++); *(tSymEx + 1) = *(tOrigEx++);
			if(ChangeSignEx) { *tSymEx = -(*tSymEx); *(tSymEx + 1) *= -(*(tSymEx + 1)); }
			tSymEx += 2;

			*tSymEz = *(tOrigEz++); *(tSymEz + 1) = *(tOrigEz++);
			if(ChangeSignEz) { *tSymEz = -(*tSymEz); *(tSymEz + 1) *= -(*(tSymEz + 1)); }
			tSymEz += 2;
		}
	}

	void SwapDataInEnergySlice(float* pOrigDataEx, float* pOrigDataEz, float* pSymDataEx, float* pSymDataEz, bool treatEx=true, bool treatEz=true)
	{
		float auxE;
		float *tOrigEx = pOrigDataEx, *tSymEx = pSymDataEx;
		float *tOrigEz = pOrigDataEz, *tSymEz = pSymDataEz;
		for(int ie = 0; ie < ne; ie++)
		{
			if(treatEx)
			{
				auxE = *tOrigEx; *(tOrigEx++) = *tSymEx; *(tSymEx++) = auxE;
				auxE = *tOrigEx; *(tOrigEx++) = *tSymEx; *(tSymEx++) = auxE;
			}
			if(treatEz)
			{
				auxE = *tOrigEz; *(tOrigEz++) = *tSymEz; *(tSymEz++) = auxE;
				auxE = *tOrigEz; *(tOrigEz++) = *tSymEz; *(tSymEz++) = auxE;
			}
		}
	}

	//double EstimWaveFrontRadiusOnAxis(char x_or_y)
	//{
	//	double start, step, cen;
	//	long n = 0;

	//	if((x_or_y == 'x') || (x_or_y == 'X'))
	//	{
	//		start = xStart; step = xStep; n = nx; cen = xc;
	//	}
	//	else if((x_or_y == 'y') || (x_or_y == 'Y') || (x_or_y == 'z') || (x_or_y == 'Z'))
	//	{
	//		start = zStart; step = zStep; n = nz; cen = zc;
	//	}
	//	
	//	if(step == 0.) return 0.;
	//	if(n < 3) return 0.;

	//	long ic = long((cen - start)/step + 1.e-13);
	//	if((cen - (ic*step + start)) > 0.5*step) ic++;
	//	long nmi1 = n - 1;
	//	if((ic < 0) || (ic > nmi1)) return 0.;

	//	long np = 5, icm2 = ic - 2, icm1 = ic - 1, icp1 = ic + 1, icp2 = ic + 2;
	//	if(icm2 < 0) np--;
	//	if(icm1 < 0) np--;
	//	if(icp1 > nmi1) np--;
	//	if(icp2 > nmi1) np--;
	//	if(np < 3) return 0.;

	//	//OC: to continue from here
	//
	//	return 0.;
	//}

	void GetWaveFrontNormal(double x, double y, double& nx, double& ny)
	{//to improve!
		nx = (x - xc)/RobsX; //check sign
		ny = (y - zc)/RobsZ; //check sign
	}

	double GetWfrMiddleHor() { return xStart + 0.5*xStep*(nx - 1);}
	double GetWfrMiddleVer() { return zStart + 0.5*zStep*(nz - 1);}

	int EnsureFreqRepres()
	{
		if(PresT) return SetRepresFT('F');
		return 0;
	}

	void SetAvgPhotEnergyFromLimits()
	{
		avgPhotEn = eStart; //averarage photon energy for time-domain simulations
		if(ne > 1)
		{
			avgPhotEn = eStart + 0.5*(ne - 1)*eStep;
		}
	}

	bool ElecFieldIsDefined()
	{
		if((pBaseRadX == 0) && (pBaseRadZ == 0)) return false;
		else return true;
	}

	long long GetIntNumPts(char dep) //OC18082018
	{//To add dependence on type of intensity when/if it is extended to mutual intensity and related charact.
		long long resNp = 0;
		if(dep == 0) resNp = ne;
		else if(dep == 1) resNp = nx;
		else if(dep == 2) resNp = nz;
		else if(dep == 3) resNp = ((long long)nx)*((long long)nz);
		else if(dep == 4) resNp = ((long long)ne)*((long long)nx);
		else if(dep == 5) resNp = ((long long)ne)*((long long)nz);
		else if(dep == 6) resNp = ((long long)ne)*((long long)nx)*((long long)nz);
		return resNp;
	}

	void DeleteElecFieldArrays()
	{
		if(pBaseRadX != 0) delete[] pBaseRadX; pBaseRadX = 0;
		if(pBaseRadZ != 0) delete[] pBaseRadZ; pBaseRadZ = 0;
	}
};

//*************************************************************************

#endif
