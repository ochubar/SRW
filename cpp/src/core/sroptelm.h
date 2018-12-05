/************************************************************************//**
 * File: sroptelm.h
 * Description: Optical element (general functions) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTELM_H
#define __SROPTELM_H

#include <stdlib.h> //required by some (buggy?) version of GCC
#include <cstdlib> //required?

#include "gmtrans.h"
#include "gmvect.h"

//#ifdef __IGOR_PRO__
//#ifndef __SRSEND_H
//#include "srsend.h"
//#endif
//#endif

#include "smartptr.h"
#include "srstraux.h"
#include "srobject.h"
#include "srinterf.h"

#include <list>

#ifdef __MAC__
#include "Memory.h"
#endif

//*************************************************************************

extern srTIntVect gVectWarnNos;

class srTGenOptElem;
//struct srTParPrecWfrPropag;

//typedef CSmartPtr<srTGenOptElem> srTGenOptElemHndl;
typedef CSmartPtr<CGenObject> srTGenOptElemHndl;
typedef list<srTGenOptElemHndl, allocator<srTGenOptElemHndl> > srTGenOptElemHndlList;
typedef list<srTGenOptElem*, allocator<srTGenOptElem*> > srTGenOptElemPtrList;

#ifndef srTransHndl
typedef CSmartPtr<gmTrans> srTransHndl;
#endif

//*************************************************************************

class srTGenOptElem : public CGenObject {

	double a2c, a4c, a6c, a8c, a10c, a12c;
	double a3s, a5s, a7s, a9s, a11s, a13s;

protected:

	srTransvLimits WfrTransmLimits; 
	//Used for estimation of resizing at propagation
	//Should be actually defined for Apertures and Gen. Transmission elements

	bool m_PropWfrInPlace;
	//if "true", no previous electric field is necessary to perform propagation through the optical element
	//"true" for most optical elements;
	//added for Grating, when it is "false" (i.e. previous electric field is necessary) 

	double HalfPI, PI, TwoPI, ThreePIdTwo, One_dTwoPI; // Constants

public:
	int ErrorCode;

	static int SetupOpticalElement(srTStringVect*, srTDataMD*, srTSRWRadStructAccessData*, srTGenOptElemHndl&);

	srTGenOptElem() 
	{
		Initialize();
	}
	virtual ~srTGenOptElem() {}

	void Initialize()
	{
		ErrorCode = 0;

		HalfPI = 1.5707963267949;
		PI = 3.141592653590;
		TwoPI = 6.2831853071796;
		ThreePIdTwo = 4.7123889803847;
		One_dTwoPI = 0.1591549430919;
		a2c = -0.5; a4c = 0.041666666666667; a6c = -0.0013888888888889; a8c = 0.000024801587301587; a10c = -2.755731922E-07;
		a3s = -0.16666666666667; a5s = 0.0083333333333333; a7s = -0.0001984126984127; a9s = 2.755731922E-06; a11s = -2.505210839E-08;
	
		m_PropWfrInPlace = true; //to modify in derived classes, if necessary
	}

	//virtual int PropagateRadiation(srTSRWRadStructAccessData*, int) { return 0;}
	//virtual int PropagateRadiation(srTSRWRadStructAccessData*, int, srTRadResizeVect&) { return 0;}
	virtual int PropagateRadiation(srTSRWRadStructAccessData*, srTParPrecWfrPropag&, srTRadResizeVect&) { return 0;}

	virtual int PropagateRadMoments(srTSRWRadStructAccessData*, srTMomentsRatios*) { return 0;}
	virtual int PropagateWaveFrontRadius(srTSRWRadStructAccessData*) { return 0;}
	virtual int PropagateWaveFrontRadius1D(srTRadSect1D*) { return 0;}
	virtual int Propagate4x4PropMatr(srTSRWRadStructAccessData*) { return 0;}
	virtual int PropagateRadiationSimple(srTSRWRadStructAccessData*) { return 0;}
	virtual int PropagateRadiationSimple1D(srTRadSect1D*) { return 0;}
	virtual int PropagateRadiationSingleE_Meth_0(srTSRWRadStructAccessData* pRadAccessData, srTSRWRadStructAccessData* pPrevRadData) { return 0;}

	virtual int RangeShouldBeAdjustedAtPropag() { return 1;}
	virtual int ResolutionShouldBeAdjustedAtPropag() { return 1;}

	virtual void RadPointModifier(srTEXZ&, srTEFieldPtrs&) {}
	virtual void RadPointModifier1D(srTEXZ&, srTEFieldPtrs&) {}

	virtual int MakePostPropagationProc(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter);
	virtual int EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz) 
	{
		const double MinPo = 40;
		MinNx = MinNz = MinPo;
		return 0;
	}
	virtual int AllowAutoSwitchToUndersamplingMode() { return 1;} //OC17122006
	virtual int PropagateRadiationTest(srTSRWRadStructAccessData* pInRadAccessData, srTSRWRadStructAccessData* pOutRadAccessData);

	virtual void AddPtrOfActualOptElem(srTGenOptElemPtrList& ActOptElemsList)
	{
		ActOptElemsList.push_back(this); //do this for all opt. comp. excepted container
	}

	virtual void ExtractTransmCharact(int CharType, double xc, double xr, int nx, double zc, double zr, int nz, float* pData);
	virtual double RadOptPathDiff(srTEXZ&) { return 0;}
	//virtual void OutOptElemInfo(char** pDescrStr, int* pLenDescr) {} //this functionality is in Optical Element constructors

	int ChangeWfrRepres(srTSRWRadStructAccessData*, int);
	int ChangeWfrRepresMeth_0(srTSRWRadStructAccessData*);
	int ChangeWfrRepresMeth_1(srTSRWRadStructAccessData*);

	int CheckRadStructForPropagation(srTSRWRadStructAccessData*);

	//int GenAuxPropagateRadMoments(srTSRWRadStructAccessData*, float**, float**, srTMomentsRatios*);
	int GenAuxPropagateRadMoments(srTSRWRadStructAccessData*, double**, double**, srTMomentsRatios*); //OC130311
	//void SetupMxxMzzArr(srTSRWRadStructAccessData*, float*, float*);
	void SetupMxxMzzArr(srTSRWRadStructAccessData*, double*, double*); //OC130311

	int GenAuxPropagate4x4PropMatr(srTSRWRadStructAccessData*, double*, double*);

	//int PropagateRadiationMeth_2(srTSRWRadStructAccessData*, srTRadResizeVect&);
	//int PropagateRadiationSingleE_Meth_2(srTSRWRadStructAccessData*, srTRadResizeVect&);
	int PropagateRadiationMeth_2(srTSRWRadStructAccessData*, srTParPrecWfrPropag&, srTRadResizeVect&);
	int PropagateRadiationSingleE_Meth_2(srTSRWRadStructAccessData*, srTParPrecWfrPropag&, srTRadResizeVect&);

	virtual int PropagateRadiationMeth_0(srTSRWRadStructAccessData* pRadAccessData); //moved from derived classes: loops over E, calls derived PropagateRadiationSingleE_Meth_0
	void FindWidestWfrMeshParam(vector<srTSRWRadStructAccessData>& vRadSlices, srTSRWRadStructAccessData* pRad, bool keepConstNumPoints);
	int ReInterpolateWfrDataOnNewTransvMesh(vector<srTSRWRadStructAccessData>& vRadSlices, srTSRWRadStructAccessData* pAuxRadSingleE, srTSRWRadStructAccessData* pRadRes);
	int ReInterpolateWfrSliceSingleE(srTSRWRadStructAccessData& oldRadSingleE, srTSRWRadStructAccessData& newRadMultiE, int ie);
	
	int SetupCharacteristicSections1D(srTSRWRadStructAccessData*, srTRadSect1D*);
	//int DefinePropagScenario(srTSRWRadStructAccessData*, srTPredictedPropagData1D*, srTPropagScenario1D*);
	//int DefinePropagScenario1D(srTRadSect1D&, srTPredictedPropagData1D&, srTPropagScenario1D&);
	int DefinePropagScenario(srTSRWRadStructAccessData*, srTParPrecWfrPropag&, srTPredictedPropagData1D*, srTPropagScenario1D*);
	int DefinePropagScenario1D(srTRadSect1D&, srTParPrecWfrPropag&, srTPredictedPropagData1D&, srTPropagScenario1D&);

	//int EstimateMinimalRangeAllowed1D(srTRadSect1D&, srTRadResize1D&);
	//int TuneRangeNotDegradingPrec1D(srTRadSect1D&, srTRadResize1D&);
	//int TuneAndKeepResolution1D(srTRadSect1D&, srTRadResize1D&, srTFringeInfo&);
	int TuneRangeNotDegradingPrec1D(srTRadSect1D&, srTParPrecWfrPropag&, srTRadResize1D&);
	int TuneAndKeepResolution1D(srTRadSect1D&, srTParPrecWfrPropag&, srTRadResize1D&, srTFringeInfo&);

	int CountFringes(srTRadSect1D&, srTIntVect&, char, srTDoubleVect&);
	int AnalizeFringes(srTRadSect1D&, char, srTFringeInfo&);
	int AnalizeFringes2D(srTSRWRadStructAccessData* pRadAccessData, srTFringeInfo* FringeInfo);
	//int CheckIfScenarioCanBeExecutedOrSuggestReduced(srTSRWRadStructAccessData*, srTRadSect1D*, srTPredictedPropagData1D*, srTPropagScenario1D*);
	//int FindPostResizeForRange1D(srTRadSect1D&, srTRadResize1D&);
	//int FindRelPrecForRangeOverWfr1D(srTRadSect1D&, srTRadResize1D&, double, char, float&);
	int CheckIfScenarioCanBeExecutedOrSuggestReduced(srTSRWRadStructAccessData*, srTRadSect1D*, srTParPrecWfrPropag&, srTPredictedPropagData1D*, srTPropagScenario1D*);
	int FindRelPrecForRangeOverWfr1D(srTRadSect1D&, srTParPrecWfrPropag&, srTRadResize1D&, double, char, float&);
	int FindPostResizeForRange1D(srTRadSect1D&, srTParPrecWfrPropag&, srTRadResize1D&);

	void FindThresholdBorders(srTRadSect1D&, double, char, long&, long&);
	int RecomputeRadMomentsIfPossible(srTSRWRadStructAccessData*);
	void CheckAndCorrectSecondOrderRadAngMoments(srTSRWRadStructAccessData*);
	void CheckRelCenPosAndSetPostResizeParamPmIfNecessary(long np, long iLeftThreshBorder, long iRightThreshBorder, srTRadResize1D& PostResizeParam, bool ModifyPmEvenIfCenPosIsNotSet);

	//int CheckAndSuggestNextValue1D(srTRadSect1D&, char, srTAuxTestValues&, srTRelAndAbsTolerance&);
	//int CheckAndSuggestNextValueWithProp1D(srTRadSect1D&, char, srTAuxTestValues&, srTRelAndAbsTolerance&);
	//char CheckRangeAndSuggestNextValue1D(srTRadSect1D&, double, double&);

	//void FindMaximumAbsReE(srTRadSect1D&, float&, long&, float&, long&);
	void FindMaximumAbsReE(srTRadSect1D&, float&, long long&, float&, long long&);
	//float MaximumAbsReEx(srTRadSect1D&, long&);
	//float ClosestLocMaxAbsReE_FromRight(srTRadSect1D&, long&, char);
	//float ClosestLocMaxAbsReE_FromLeft(srTRadSect1D&, long&, char);
	int EnsureCoordRepres(srTSRWRadStructAccessData* pRadAccessData)
	{
		if(pRadAccessData->Pres != 0) return SetRadRepres(pRadAccessData, 0);
		else return 0;
	}
	//int CheckResolution1D(srTRadSect1D&, char&, srTAuxTestIntegValues&, srTRelAndAbsTolerance&);
	//int CheckPrecAtPropag1D(srTRadSect1D&, char, char&, srTAuxTestIntegValues&, srTRelAndAbsTolerance&);
	//int IntegAfterPropag1D(srTRadSect1D&, char, char, char, float&, float&, srTRelAndAbsTolerance&);
	inline float IntegrateElField1D(srTRadSect1D&, char, float&);

	int FillOutRadFromInRad(srTSRWRadStructAccessData*, srTSRWRadStructAccessData*);

	int TraverseRadZXE(srTSRWRadStructAccessData*);
	int TraverseRad1D(srTRadSect1D*);

	int ExtractRadSliceConstE(srTSRWRadStructAccessData*, long, float*&, float*&, bool forceCopyField=false); //OC120908
	int SetupRadSliceConstE(srTSRWRadStructAccessData*, long, float*, float*);
	inline void SetupRadXorZSectFromSliceConstE(float*, float*, long, long, char, long, float*, float*);

	int ExtractRadSectVsXorZ(srTSRWRadStructAccessData*, long, long, char, float*, float*);
	int SetupSectionArraysVsXandZ(srTSRWRadStructAccessData*, srTRadSect1D&, srTRadSect1D&);

	int SetupNewRadStructFromSliceConstE(srTSRWRadStructAccessData* pRadAccessData, long, srTSRWRadStructAccessData*& pRadDataSingleE);
	//int UpdateGenRadStructFromSlicesConstE(srTSRWRadStructAccessData*, srTSRWRadStructAccessData*);
	//int UpdateGenRadStructSliceConstE_Meth_0(srTSRWRadStructAccessData*, int, srTSRWRadStructAccessData*);
	//OC28102018: modified by S.Yakubov to adopt the code for OpenMP parallelization
	int UpdateGenRadStructSliceConstE_Meth_0(srTSRWRadStructAccessData*, int, srTSRWRadStructAccessData*, int update_mode=0);

	int UpdateGenRadStructSliceConstE_Meth_2(srTSRWRadStructAccessData*, int, srTSRWRadStructAccessData*);
	int RemoveSliceConstE_FromGenRadStruct(srTSRWRadStructAccessData*, long);

	//int SetRadRepres(srTSRWRadStructAccessData*, char);
	int SetRadRepres(srTSRWRadStructAccessData*, char, double* ar_xStartInSlicesE=0, double* ar_zStartInSlicesE=0);
	int SetRadRepres1D(srTRadSect1D*, char);

	int SetupWfrEdgeCorrData(srTSRWRadStructAccessData*, float*, float*, srTDataPtrsForWfrEdgeCorr&);
	//inline void SetupExpCorrArray(float*, long, double, double, double);
	inline void SetupExpCorrArray(float*, long long, double, double, double);
	void MakeWfrEdgeCorrection(srTSRWRadStructAccessData*, float*, float*, srTDataPtrsForWfrEdgeCorr&);

	int SetupWfrEdgeCorrData1D(srTRadSect1D*, float*, float*, srTDataPtrsForWfrEdgeCorr1D&);
	void MakeWfrEdgeCorrection1D(srTRadSect1D*, float*, float*, srTDataPtrsForWfrEdgeCorr1D&);

	int ComputeRadMoments(srTSRWRadStructAccessData*);

	int RadResizeGen(srTSRWRadStructAccessData&, srTRadResize&);
	int RadResizeGenE(srTSRWRadStructAccessData&, srTRadResize&);
	int RadResizeCore(srTSRWRadStructAccessData&, srTSRWRadStructAccessData&, srTRadResize&, char =0);
	int RadResizeCoreE(srTSRWRadStructAccessData&, srTSRWRadStructAccessData&, srTRadResize&, char =0);
	int RadResizeCore_OnlyLargerRange(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp);
	int RadResizeCore_OnlyLargerRangeE(srTSRWRadStructAccessData& OldRadAccessData, srTSRWRadStructAccessData& NewRadAccessData, srTRadResize& RadResizeStruct, char PolComp);

	//inline void GetCellDataForInterpol(float*, long, long, srTInterpolAuxF*);
	inline void GetCellDataForInterpol(float*, long long , long long, srTInterpolAuxF*);
	inline void SetupCellDataI(srTInterpolAuxF*, srTInterpolAuxF*);
	//char WaveFrontTermCanBeTreated(srTSRWRadStructAccessData&);
	//char WaveFrontTermCanBeTreated(srTSRWRadStructAccessData&, bool checkBenefit=true); //OC06012017 (uncommented after some fixes in bool srTSRWRadStructAccessData::CheckIfQuadTermTreatIsBenefit(char, char))
	//char WaveFrontTermCanBeTreated(srTSRWRadStructAccessData&, bool checkBenefit=false); //OC05012017 (changed to checkBenefit=false to resolve problem of resizing in near field at strong under-sampling)
	char WaveFrontTermCanBeTreated(srTSRWRadStructAccessData&, bool checkBenefit=false); //OC29032017 (changed again to checkBenefit=false to resolve problem of resizing of wiggler radiation at strong under-sampling, the ELETTRA SCW case)

	void TreatStronglyOscillatingTerm(srTSRWRadStructAccessData&, char, char =0, int ieOnly =-1);
	//void TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData&, float*, float, float, float, float, char, char =0, int =-1);
	void TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData&, double*, double, double, double, double, char, char =0, int =-1); //OC260114
	//void TreatStronglyOscillatingTermIrregMesh(srTSRWRadStructAccessData&, double*, double, double, double, double, char, char =0, int =-1, double =1, double =1); //OC220214

	inline void SetupInterpolAux02(srTInterpolAuxF*, srTInterpolAux01*, srTInterpolAux02*);
	inline void SetupInterpolAux02_LowOrder(srTInterpolAuxF*, srTInterpolAux01*, srTInterpolAux02*);
	inline void InterpolF(srTInterpolAux02*, double, double, float*, int);
	inline void InterpolFI(srTInterpolAux02*, double, double, float*, int);
	inline void InterpolF_LowOrder(srTInterpolAux02*, double, double, float*, int);
	inline void InterpolFI_LowOrder(srTInterpolAux02*, double, double, float*, int);
	inline double InterpLin(double r, double f1, double f2) { return f1 + r*(f2 - f1);}
	inline void ImproveReAndIm(float*, float*);
	inline int CheckForLowOrderInterp(srTInterpolAuxF*, srTInterpolAuxF*, int, int, srTInterpolAux01*, srTInterpolAux02*, srTInterpolAux02*);

	int RadResizeGen1D(srTRadSect1D&, srTRadResize1D&);
	int RadResizeCore1D(srTRadSect1D&, srTRadSect1D&, srTRadResize1D&);
	char WaveFrontTermCanBeTreated1D(srTRadSect1D&);
	void TreatStronglyOscillatingTerm1D(srTRadSect1D&, char);
	//inline void GetCellDataForInterpol1D(float*, long, srTInterpolAuxF_1D*);
	inline void GetCellDataForInterpol1D(float*, long long, srTInterpolAuxF_1D*);
	inline void SetupCellDataI1D(srTInterpolAuxF_1D*, srTInterpolAuxF_1D*);
	inline int CheckForLowOrderInterp1D(srTInterpolAuxF_1D*, srTInterpolAuxF_1D*, int, srTInterpolAux01_1D*, srTInterpolAux02_1D*, srTInterpolAux02_1D*);
	inline void SetupInterpolAux02_LowOrder1D(srTInterpolAuxF_1D*, srTInterpolAux01_1D*, srTInterpolAux02_1D*);
	inline void SetupInterpolAux02_1D(srTInterpolAuxF_1D*, srTInterpolAux01_1D*, srTInterpolAux02_1D*);
	inline void InterpolF_LowOrder1D(srTInterpolAux02_1D*, double, float*, int);
	inline void InterpolFI_LowOrder1D(srTInterpolAux02_1D*, double, float*, int);
	inline void InterpolF1D(srTInterpolAux02_1D*, double, float*, int);
	inline void InterpolFI1D(srTInterpolAux02_1D*, double, float*, int);

	int RadRearrangeOnRegularMesh(srTSRWRadStructAccessData* pRadAccessData, float* CoordX, float* CoordZ);

	int GenExtractPhase(srTWaveAccessData&, double*, double*, int, int);
	//int ExtractPhase1D(srTWaveAccessData&, double*, double*, long, double);
	int ExtractPhase1D(srTWaveAccessData&, double*, double*, long long, double);
	inline double PredictPhase(double, float, float, float, float);
	inline double FormalPhase(float, float);
	inline double FormalMag(float, float, double);

	//inline void MultSquareMatrByVect(float**, float*, int, float*);
	inline void MultSquareMatrByVect(double**, double*, int, double*); //OC130311
	inline void CosAndSin(double, float&, float&);
	inline void FindLowestAndUppestPoints(TVector3d&, TVector3d*, int, int&, int&);
	inline void ReflectVect(TVector3d& N, TVector3d& V);
	inline void FindLineIntersectWithPlane(TVector3d* Plane, TVector3d* Line, TVector3d& IntersectP);
	inline void TreatPhaseShift(srTEFieldPtrs& EPtrs, double PhShift);

	inline long IntegerOffsetCoord(double xStart, double xStep, double xVal);
	//void FindMinMaxRatio(float*, float*, int, float&, float&);
	void FindMinMaxRatio(double*, double*, int, double&, double&); //OC130311
	int MakeSimpleOversamplingTestAndCorrection(srTSRWRadStructAccessData*);

	inline char ChooseTreatExOrEzBasedOnMax(srTRadSect1D&);
	inline int ResizePropagateAndAnalizeFringes1D(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam, char TreatExOrEz, srTFringeInfo& PropFringeInfo);
	inline int CheckIfOversamplingIsReal(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam, char TreatExOrEz, srTFringeInfo& PropFringeInfo);
	inline double SuggestResolResizeCoef(double PoPerFringeBeforePropag, double PoPerFringeAfterPropag, long AmOfFringesBeforePropag);
	void EstimateMemoryNeededForPropag(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario, double&, double&); // in bytes
	double CheckMemoryAvailable();
	inline int MemoryIsSufficientForProp(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario, double& ExtraSizeBeforeProp, double& ExtraSizeAfterProp);
	int MemoryIsSufficientForTwoResize(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& Resize1, srTRadResize& Resize2);
	int MemoryIsSufficientForResize(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& Resize);
	int SuggestScenarioThatFitsMemory(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario);
	void CorrectResParMinNumPo(long Np, srTRadResize1D& ResizeBefore, srTRadResize1D& ResizeAfter);
	inline void FindCenterWithRespectToThresh(srTRadSect1D& Sect1D, double AbsThresh, char TreatExOrEz, long& ic, double& xc);
	//int FindPostResizeCenterCorrection(srTRadSect1D& Sect1D, srTPropagScenario1D& PropagScenario);
	int FindPostResizeCenterCorrection(srTRadSect1D& Sect1D, srTParPrecWfrPropag&, srTPropagScenario1D& PropagScenario);
	int CheckPostResizeCenterCorrection(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter);
	int PostResizeAndTryToImproveResolInSmallSpot(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter);
	int CheckIfSpotShouldBeResized(srTSRWRadStructAccessData* pRadAccessData, char&, srTRadResize& ImproveRes);
	int AnalizeFringesAroundPoint(srTRadSect1D& Sect1D, char TreatExOrEz, long ic, double& PoPerFr, double& FringeSize);
	int CheckRMS_Sizes1D(srTRadSect1D& Sect1D, char TreatExOrEz, double& Xc, double& SigmaX);
	int CheckWidthMax1D(srTRadSect1D& Sect1D, char TreatExOrEz, double& Xc, double& DeltaX);
	void FindIntensityBorders1D(srTRadSect1D& Sect1D, char TreatExOrEz, double RelZeroTolForIntens, long& iFirst, long& iLast);
	int TuneStepToKeepInterpLimitsTheSameAtResize(srTSRWRadStructAccessData& SRWRadStructAccessData, srTSRWRadStructAccessData& NewSRWRadStructAccessData, srTRadResize&, char, long);
	inline void RejectSmallResize(srTRadResize& Resize);
	void TransferResizeParam(srTPropagScenario1D* PropagScenario, srTRadResize& ResBefore, srTRadResize& ResAfter);
	double ExtraMemSizeForResize(long nxCurRad, long nzCurRad, double pxm, double pxd, double pzm, double pzd, char Mode);
	double ExtraMemSizeForResizeE(long neCurRad, long nxCurRad, long nzCurRad, double pem, double ped, char Mode);
	void SteerPostResizeParam(srTSRWRadStructAccessData* pRadAccessData, srTRadResize& ResAfter);
	inline double MaxIntInHorString(long iz, srTSRWRadStructAccessData* pRadAccessData);
	inline double MaxIntInVertString(long ix, long izSt, long izFi, srTSRWRadStructAccessData* pRadAccessData);

	int TryToSetUnderSamplingMode(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* pSect1D, srTPropagScenario1D* PropagScenario, char& UnderSamplingModeWasSet);
	char SuitableConditionsForUnderSamplingMode(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario);
	char UnderSamplingModeCanBeSuggested(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario);
	int TryToRemoveUndersamplingByResizing(srTSRWRadStructAccessData* pRadAccessData);
	int RemoveUndersamplingByResizingOrStop(srTSRWRadStructAccessData* pRadAccessData);
	void ShowCurrentOverSamplingFactors(srTPropagScenario1D* PropagScenario, double& Fx, double& Fz);
	int EstimateNominalNpForUnderSampling(srTSRWRadStructAccessData* pRadAccessData, srTRadSect1D* Sect1D, double& NomNxForUnderSampl, double& NomNzForUnderSampl);
	long EstimateMinNpForQuadTerm(double en, double R, double xStart, double xEnd);

	int ReduceBiggerResizeParamUntilFitMemory(srTSRWRadStructAccessData& RadAccessData, srTRadResize& RadResize, double& UnderSamplingX, double& UnderSamplingZ);

	double FindLongCoordOfHorBaseVect(double vLx, double vLy, double vLz, double vHx, double vHy);
};

//*************************************************************************

//inline void srTGenOptElem::GetCellDataForInterpol(float* pSt, long PerX_Old, long PerZ_Old, srTInterpolAuxF* tF)
inline void srTGenOptElem::GetCellDataForInterpol(float* pSt, long long PerX_Old, long long PerZ_Old, srTInterpolAuxF* tF)
{// Fills Re and Im parts of Ex or Ez
	float *pf00 = pSt; tF->f00 = *pf00;
	float *pf10 = pf00 + PerX_Old; tF->f10 = *pf10;
	float *pf20 = pf10 + PerX_Old; tF->f20 = *pf20;
	float *pf30 = pf20 + PerX_Old; tF->f30 = *pf30;

	float *pf01 = pf00 + PerZ_Old; tF->f01 = *pf01;
	float *pf11 = pf01 + PerX_Old; tF->f11 = *pf11;
	float *pf21 = pf11 + PerX_Old; tF->f21 = *pf21;
	float *pf31 = pf21 + PerX_Old; tF->f31 = *pf31;

	float *pf02 = pf01 + PerZ_Old; tF->f02 = *pf02;
	float *pf12 = pf02 + PerX_Old; tF->f12 = *pf12;
	float *pf22 = pf12 + PerX_Old; tF->f22 = *pf22;
	float *pf32 = pf22 + PerX_Old; tF->f32 = *pf32;

	float *pf03 = pf02 + PerZ_Old; tF->f03 = *pf03;
	float *pf13 = pf03 + PerX_Old; tF->f13 = *pf13;
	float *pf23 = pf13 + PerX_Old; tF->f23 = *pf23;
	float *pf33 = pf23 + PerX_Old; tF->f33 = *pf33;
	tF++;

	tF->f00 = *(++pf00); tF->f10 = *(++pf10); tF->f20 = *(++pf20); tF->f30 = *(++pf30);
	tF->f01 = *(++pf01); tF->f11 = *(++pf11); tF->f21 = *(++pf21); tF->f31 = *(++pf31);
	tF->f02 = *(++pf02); tF->f12 = *(++pf12); tF->f22 = *(++pf22); tF->f32 = *(++pf32);
	tF->f03 = *(++pf03); tF->f13 = *(++pf13); tF->f23 = *(++pf23); tF->f33 = *(++pf33);
}

//*************************************************************************

//inline void srTGenOptElem::GetCellDataForInterpol1D(float* pSt, long Per_Old, srTInterpolAuxF_1D* tF)
inline void srTGenOptElem::GetCellDataForInterpol1D(float* pSt, long long Per_Old, srTInterpolAuxF_1D* tF)
{// Fills Re and Im parts of Ex or Ez
	float *pf0 = pSt; tF->f0 = *pf0;
	float *pf1 = pf0 + Per_Old; tF->f1 = *pf1;
	float *pf2 = pf1 + Per_Old; tF->f2 = *pf2;
	float *pf3 = pf2 + Per_Old; tF->f3 = *pf3;

	tF++;
	tF->f0 = *(++pf0); tF->f1 = *(++pf1); tF->f2 = *(++pf2); tF->f3 = *(++pf3);
}

//*************************************************************************

inline void srTGenOptElem::SetupCellDataI(srTInterpolAuxF* tF, srTInterpolAuxF* tI)
{
	srTInterpolAuxF* tF1 = tF + 1;

	tI->f00 = (tF->f00)*(tF->f00) + (tF1->f00)*(tF1->f00);
	tI->f10 = (tF->f10)*(tF->f10) + (tF1->f10)*(tF1->f10);
	tI->f20 = (tF->f20)*(tF->f20) + (tF1->f20)*(tF1->f20);
	tI->f30 = (tF->f30)*(tF->f30) + (tF1->f30)*(tF1->f30);

	tI->f01 = (tF->f01)*(tF->f01) + (tF1->f01)*(tF1->f01);
	tI->f11 = (tF->f11)*(tF->f11) + (tF1->f11)*(tF1->f11);
	tI->f21 = (tF->f21)*(tF->f21) + (tF1->f21)*(tF1->f21);
	tI->f31 = (tF->f31)*(tF->f31) + (tF1->f31)*(tF1->f31);

	tI->f02 = (tF->f02)*(tF->f02) + (tF1->f02)*(tF1->f02);
	tI->f12 = (tF->f12)*(tF->f12) + (tF1->f12)*(tF1->f12);
	tI->f22 = (tF->f22)*(tF->f22) + (tF1->f22)*(tF1->f22);
	tI->f32 = (tF->f32)*(tF->f32) + (tF1->f32)*(tF1->f32);

	tI->f03 = (tF->f03)*(tF->f03) + (tF1->f03)*(tF1->f03);
	tI->f13 = (tF->f13)*(tF->f13) + (tF1->f13)*(tF1->f13);
	tI->f23 = (tF->f23)*(tF->f23) + (tF1->f23)*(tF1->f23);
	tI->f33 = (tF->f33)*(tF->f33) + (tF1->f33)*(tF1->f33);

	tI->SetUpAvg(); tI->NormalizeByAvg();
}
//*************************************************************************

inline void srTGenOptElem::SetupCellDataI1D(srTInterpolAuxF_1D* tF, srTInterpolAuxF_1D* tI)
{
	srTInterpolAuxF_1D* tF1 = tF + 1;
	tI->f0 = (tF->f0)*(tF->f0) + (tF1->f0)*(tF1->f0);
	tI->f1 = (tF->f1)*(tF->f1) + (tF1->f1)*(tF1->f1);
	tI->f2 = (tF->f2)*(tF->f2) + (tF1->f2)*(tF1->f2);
	tI->f3 = (tF->f3)*(tF->f3) + (tF1->f3)*(tF1->f3);

	tI->SetUpAvg(); tI->NormalizeByAvg();
}

//*************************************************************************

inline void srTGenOptElem::SetupInterpolAux02(srTInterpolAuxF* pF, srTInterpolAux01* pC, srTInterpolAux02* pA)
{
	pA->Ax0z0 = pF->f11;
	pA->Ax0z1 = (-2*pF->f10 - 3*pF->f11 + 6*pF->f12 - pF->f13)*pC->cAx0z1;
	pA->Ax0z2 = (pF->f10 + pF->f12 - 2*pF->f11)*pC->cAx0z2;
	pA->Ax0z3 = (pF->f13 - pF->f10 + 3*(pF->f11 - pF->f12))*pC->cAx0z3;
	pA->Ax1z0 = (-2*pF->f01 - 3*pF->f11 + 6*pF->f21 - pF->f31)*pC->cAx1z0;
	pA->Ax1z1 = (4*pF->f00 + 6*(pF->f01 + pF->f10 - pF->f23 - pF->f32) - 12*(pF->f02 + pF->f20) + 2*(pF->f03 + pF->f30) + 9*pF->f11 - 18*(pF->f12 + pF->f21) + 3*(pF->f13 + pF->f31) + 36*pF->f22 + pF->f33)*pC->cAx1z1;
	pA->Ax1z2 = (-2*(pF->f00 + pF->f02 - pF->f31) + 4*pF->f01 - 3*(pF->f10 + pF->f12) + 6*(pF->f11 + pF->f20 + pF->f22) - 12*pF->f21 - pF->f30 - pF->f32)*pC->cAx1z2;
	pA->Ax1z3 = (2*(pF->f00 - pF->f03) + 6*(-pF->f01 + pF->f02 - pF->f20 + pF->f23) + 3*(pF->f10 - pF->f13 - pF->f31 + pF->f32) + 9*(pF->f12 - pF->f11) + 18*(pF->f21 - pF->f22) + pF->f30 - pF->f33)*pC->cAx1z3;
	pA->Ax2z0 = (pF->f01 + pF->f21 - 2*pF->f11)*pC->cAx2z0;
	pA->Ax2z1 = (2*(-pF->f00 + pF->f13 - pF->f20) - 3*(pF->f21 + pF->f01) + 6*(pF->f02 + pF->f11 + pF->f22) + 4*pF->f10 - 12*pF->f12 - pF->f23 - pF->f03)*pC->cAx2z1;
	pA->Ax2z2 = (pF->f00 + pF->f02 + pF->f22 + pF->f20 - 2*(pF->f01 + pF->f10 + pF->f12 + pF->f21) + 4*pF->f11)*pC->cAx2z2;
	pA->Ax2z3 = (-pF->f00 + pF->f03 - pF->f20 + pF->f23 + 3*(pF->f01 - pF->f02 + pF->f21 - pF->f22) + 2*(pF->f10 - pF->f13) + 6*(pF->f12 - pF->f11))*pC->cAx2z3;
	pA->Ax3z0 = (pF->f31 - pF->f01 + 3*(pF->f11 - pF->f21))*pC->cAx3z0;
	pA->Ax3z1 = (2*(pF->f00 - pF->f30) + 3*(pF->f01 - pF->f13 + pF->f23 - pF->f31) + 6*(-pF->f02 - pF->f10 + pF->f20 + pF->f32) + 9*(pF->f21 - pF->f11) + 18*(pF->f12 - pF->f22) + pF->f03 - pF->f33)*pC->cAx3z1;
	pA->Ax3z2 = (pF->f30 + pF->f32 - pF->f00 - pF->f02 + 2*(pF->f01 - pF->f31) + 3*(pF->f10 + pF->f12 - pF->f20 - pF->f22) + 6*(pF->f21 - pF->f11))*pC->cAx3z2;
	pA->Ax3z3 = (pF->f00 - pF->f03 - pF->f30 + pF->f33 + 3*(-pF->f01 + pF->f02 - pF->f10 + pF->f13 + pF->f20 - pF->f23 + pF->f31 - pF->f32) + 9*(pF->f11 - pF->f12 - pF->f21 + pF->f22))*pC->cAx3z3;
}

//*************************************************************************

inline void srTGenOptElem::SetupInterpolAux02_1D(srTInterpolAuxF_1D* pF, srTInterpolAux01_1D* pC, srTInterpolAux02_1D* pA)
{
	pA->A0 = pF->f1;
	pA->A1 = (-2.*pF->f0 - 3.*pF->f1 + 6.*pF->f2 - pF->f3)*pC->cA1;
	pA->A2 = (pF->f0 + pF->f2 - 2.*pF->f1)*pC->cA2;
	pA->A3 = (pF->f3 - pF->f0 + 3.*(pF->f1 - pF->f2))*pC->cA3;
}

//*************************************************************************

inline void srTGenOptElem::SetupInterpolAux02_LowOrder(srTInterpolAuxF* pF, srTInterpolAux01* pC, srTInterpolAux02* pA)
{
	pA->Ax0z0 = pF->f00;
	pA->Ax1z0 = pC->cLAx1z0*(pF->f10 - pF->f00);
	pA->Ax0z1 = pC->cLAx0z1*(pF->f01 - pF->f00);
	pA->Ax1z1 = pC->cLAx1z1*(pF->f00 - pF->f01 - pF->f10 + pF->f11);
}

//*************************************************************************

inline void srTGenOptElem::SetupInterpolAux02_LowOrder1D(srTInterpolAuxF_1D* pF, srTInterpolAux01_1D* pC, srTInterpolAux02_1D* pA)
{
	pA->A0 = pF->f0;
	pA->A1 = pC->cLA1*(pF->f1 - pF->f0);
}

//*************************************************************************

inline void srTGenOptElem::InterpolF(srTInterpolAux02* A, double x, double z, float* F, int Offset)
{
	double xE2 = x*x, xz = x*z, zE2 = z*z;
	double xE3 = xE2*x, xE2z = xE2*z, xzE2 = x*zE2, zE3 = zE2*z, xE2zE2 = xE2*zE2;
	double xE3z = xE3*z, xE3zE2 = xE3*zE2, xE3zE3 = xE3*zE3, xE2zE3 = xE2*zE3, xzE3 = x*zE3;
	srTInterpolAux02* tA = A + Offset;
	for(int i=0; i<4-Offset; i++)
	{
		F[i + Offset] = (float)(tA->Ax3z3*xE3zE3 + tA->Ax3z2*xE3zE2 + tA->Ax3z1*xE3z + tA->Ax3z0*xE3 + tA->Ax2z3*xE2zE3 + tA->Ax2z2*xE2zE2 + tA->Ax2z1*xE2z + tA->Ax2z0*xE2 + tA->Ax1z3*xzE3 + tA->Ax1z2*xzE2 + tA->Ax1z1*xz + tA->Ax1z0*x + tA->Ax0z3*zE3 + tA->Ax0z2*zE2 + tA->Ax0z1*z + tA->Ax0z0);
		tA++;
	}
}

//*************************************************************************

inline void srTGenOptElem::InterpolF1D(srTInterpolAux02_1D* A, double x, float* F, int Offset)
{
	double xE2 = x*x;
	double xE3 = xE2*x;
	srTInterpolAux02_1D* tA = A + Offset;
	for(int i=0; i<4-Offset; i++)
	{
		F[i + Offset] = (float)(tA->A3*xE3 + tA->A2*xE2 + tA->A1*x + tA->A0);
		tA++;
	}
}

//*************************************************************************

inline void srTGenOptElem::InterpolFI(srTInterpolAux02* A, double x, double z, float* F, int Offset)
{
	double xE2 = x*x, xz = x*z, zE2 = z*z;
	double xE3 = xE2*x, xE2z = xE2*z, xzE2 = x*zE2, zE3 = zE2*z, xE2zE2 = xE2*zE2;
	double xE3z = xE3*z, xE3zE2 = xE3*zE2, xE3zE3 = xE3*zE3, xE2zE3 = xE2*zE3, xzE3 = x*zE3;
	srTInterpolAux02* tA = A + Offset;
	double Buf = tA->Ax3z3*xE3zE3 + tA->Ax3z2*xE3zE2 + tA->Ax3z1*xE3z + tA->Ax3z0*xE3 + tA->Ax2z3*xE2zE3 + tA->Ax2z2*xE2zE2 + tA->Ax2z1*xE2z + tA->Ax2z0*xE2 + tA->Ax1z3*xzE3 + tA->Ax1z2*xzE2 + tA->Ax1z1*xz + tA->Ax1z0*x + tA->Ax0z3*zE3 + tA->Ax0z2*zE2 + tA->Ax0z1*z + tA->Ax0z0;
	*(F + Offset) = (float)((Buf > 0.)? Buf : 0.);
}

//*************************************************************************

inline void srTGenOptElem::InterpolFI1D(srTInterpolAux02_1D* A, double x, float* F, int Offset)
{
	double xE2 = x*x;
	double xE3 = xE2*x;
	srTInterpolAux02_1D* tA = A + Offset;
	double Buf = tA->A3*xE3 + tA->A2*xE2 + tA->A1*x + tA->A0;
	*(F + Offset) = (float)((Buf > 0.)? Buf : 0.);
}

//*************************************************************************

inline void srTGenOptElem::InterpolF_LowOrder(srTInterpolAux02* A, double x, double z, float* F, int Offset)
{
	double xz = x*z;
	srTInterpolAux02* tA = A + Offset;
	for(int i=0; i<4-Offset; i++)
	{
		F[i + Offset] = (float)(tA->Ax1z1*xz + tA->Ax1z0*x + tA->Ax0z1*z + tA->Ax0z0);
		tA++;
	}
}

//*************************************************************************

inline void srTGenOptElem::InterpolF_LowOrder1D(srTInterpolAux02_1D* A, double x, float* F, int Offset)
{
	srTInterpolAux02_1D* tA = A + Offset;
	for(int i=0; i<4-Offset; i++)
	{
		F[i + Offset] = (float)(tA->A1*x + tA->A0);
		tA++;
	}
}

//*************************************************************************

inline void srTGenOptElem::InterpolFI_LowOrder(srTInterpolAux02* A, double x, double z, float* F, int Offset)
{
	double xz = x*z;
	srTInterpolAux02* tA = A + Offset;
	double Buf = tA->Ax1z1*xz + tA->Ax1z0*x + tA->Ax0z1*z + tA->Ax0z0;
	*(F + Offset) = (float)((Buf > 0.)? Buf : 0.);
}

//*************************************************************************

inline void srTGenOptElem::InterpolFI_LowOrder1D(srTInterpolAux02_1D* A, double x, float* F, int Offset)
{
	srTInterpolAux02_1D* tA = A + Offset;
	double Buf = tA->A1*x + tA->A0;
	*(F + Offset) = (float)((Buf > 0.)? Buf : 0.);
}

//*************************************************************************

inline void srTGenOptElem::ImproveReAndIm(float* pFReIm, float* pFI)
{
	float &FRe = *pFReIm, &FIm = *(pFReIm+1);
	float AppI = FRe*FRe + FIm*FIm;
	if(AppI != 0.)
	{
		float Factor = (float)sqrt(*pFI/AppI);
		FRe *= Factor; FIm *= Factor;
	}
}

//*************************************************************************

inline int srTGenOptElem::CheckForLowOrderInterp(srTInterpolAuxF* CellF, srTInterpolAuxF* CellFI, int ixRel, int izRel, srTInterpolAux01* pC, srTInterpolAux02* pA, srTInterpolAux02* pAI)
{
	if(ixRel < 0) ixRel = 0;
	if(ixRel > 2) ixRel = 2;
	if(izRel < 0) izRel = 0;
	if(izRel > 2) izRel = 2;
	srTInterpolAuxF* t = CellF;
	int iLxLz, iUxLz, iLxUz, iUxUz;
	char LowOrderCaseNoticed = 0;
	for(int i=0; i<2; i++)
	{
		iLxLz = (izRel << 2) + ixRel; iUxLz = iLxLz + 1;
		iLxUz = iLxLz + 4; iUxUz = iLxUz + 1;
		if((t->f00==0.) || (t->f10==0.) || (t->f20==0.) || (t->f30==0.) ||
		   (t->f01==0.) || (t->f11==0.) || (t->f21==0.) || (t->f31==0.) ||
		   (t->f02==0.) || (t->f12==0.) || (t->f22==0.) || (t->f32==0.) ||
		   (t->f03==0.) || (t->f13==0.) || (t->f23==0.) || (t->f33==0.))
		{
			LowOrderCaseNoticed = 1; break;
		}
		t++;
	}
	if(LowOrderCaseNoticed)
	{
		t = CellF;
		srTInterpolAuxF AuxF[2];
		srTInterpolAuxF* tAuxF = AuxF;
		for(int i=0; i<2; i++)
		{
			float BufF[] = { t->f00, t->f10, t->f20, t->f30, t->f01, t->f11, t->f21, t->f31, t->f02, t->f12, t->f22, t->f32, t->f03, t->f13, t->f23, t->f33};
			tAuxF->f00 = BufF[iLxLz]; tAuxF->f10 = BufF[iUxLz]; tAuxF->f01 = BufF[iLxUz]; tAuxF->f11 = BufF[iUxUz];
			SetupInterpolAux02_LowOrder(tAuxF, pC, pA+i);
			t++; tAuxF++;
		}

		t = CellFI;
		srTInterpolAuxF AuxFI;
		float BufFI[] = { t->f00, t->f10, t->f20, t->f30, t->f01, t->f11, t->f21, t->f31, t->f02, t->f12, t->f22, t->f32, t->f03, t->f13, t->f23, t->f33};
		AuxFI.f00 = BufFI[iLxLz]; AuxFI.f10 = BufFI[iUxLz]; AuxFI.f01 = BufFI[iLxUz]; AuxFI.f11 = BufFI[iUxUz];
		SetupInterpolAux02_LowOrder(&AuxFI, pC, pAI);
	}
	return LowOrderCaseNoticed;
}

//*************************************************************************

inline int srTGenOptElem::CheckForLowOrderInterp1D(srTInterpolAuxF_1D* CellF, srTInterpolAuxF_1D* CellFI, int iRel, srTInterpolAux01_1D* pC, srTInterpolAux02_1D* pA, srTInterpolAux02_1D* pAI)
{
	char LowOrderCaseNoticed = 0;
	if(iRel < 0) iRel = 0;
	if(iRel > 2) iRel = 2;
	srTInterpolAuxF_1D* t = CellF;
	int iL = iRel, iU= iRel + 1;
	for(int i=0; i<2; i++)
	{
		if((t->f0==0.) || (t->f1==0.) || (t->f2==0.) || (t->f3==0.))
		{
			LowOrderCaseNoticed = 1; break;
		}
		t++;
	}

	if(LowOrderCaseNoticed)
	{
		t = CellF;
		srTInterpolAuxF_1D AuxF[2];
		srTInterpolAuxF_1D* tAuxF = AuxF;
		for(int i=0; i<2; i++)
		{
			float BufF[] = { t->f0, t->f1, t->f2, t->f3};
			tAuxF->f0 = BufF[iL]; tAuxF->f1 = BufF[iU];
			SetupInterpolAux02_LowOrder1D(tAuxF, pC, pA+i);
			t++; tAuxF++;
		}

		t = CellFI;
		srTInterpolAuxF_1D AuxFI;
		float BufFI[] = { t->f0, t->f1, t->f2, t->f3};
		AuxFI.f0 = BufFI[iL]; AuxFI.f1 = BufFI[iU];
		SetupInterpolAux02_LowOrder1D(&AuxFI, pC, pAI);
	}
	return LowOrderCaseNoticed;
}

//*************************************************************************

inline double srTGenOptElem::PredictPhase(double PhPrev, float RePrev, float ImPrev, float Re, float Im)
{
	double MagE2 = RePrev*RePrev + ImPrev*ImPrev;
	if(MagE2 == 0.) return PhPrev;
	double dPh = (RePrev*(Im - ImPrev) - ImPrev*(Re - RePrev))/MagE2;
	return PhPrev + dPh;
}

//*************************************************************************

inline double srTGenOptElem::FormalPhase(float Re, float Im)
{
	const double HalhPi = 1.5707963267949;
	const double Pi = 3.1415926535898;
	if(Re != 0.) 
	{
		if(Im <= 0.)
		{
			if(Re < 0.) return atan(double(Im/Re)) - Pi;
			else return atan(double(Im/Re));
		}
		else
		{
			if(Re < 0.) return atan(double(Im/Re)) + Pi;
			else return atan(double(Im/Re));
		}
	}
	else
	{
		if(Im == 0.) return  0.;
		else return (Im > 0.)? HalhPi : -HalhPi;
	}
}

//*************************************************************************

inline double srTGenOptElem::FormalMag(float Re, float Im, double Ph)
{
	double CosPh = cos(Ph);
	if(CosPh != 0.) return Re/CosPh;
	else return Im/sin(Ph);
}

//*************************************************************************

//inline void srTGenOptElem::MultSquareMatrByVect(float** b, float* c, int n, float* a)
inline void srTGenOptElem::MultSquareMatrByVect(double** b, double* c, int n, double* a) //OC130311
{// a = b*c
	int i, j;
	//float* ta = a;
	double* ta = a; //OC130311
	for(i=0; i<n; i++) *(ta++) = 0.;
	for(i=0; i<n; i++)
	{
		//float* pai = a + i;
		//float** pbi = b + i;
		double* pai = a + i; //OC130311
		double** pbi = b + i;
		for(j=0; j<n; j++) *pai += (*pbi)[j]*c[j];
	}
}

//*************************************************************************

inline void srTGenOptElem::CosAndSin(double x, float& Cos, float& Sin)
{
	if((x < -1.E+08) || (x > 1.E+08)) { Cos = (float)cos(x); Sin = (float)sin(x); return;} //OC13112011

	//x -= TwoPI*((long)(x*One_dTwoPI));
	x -= TwoPI*((long long)(x*One_dTwoPI));
	if(x < 0.) x += TwoPI;

	char ChangeSign=0;
	if(x > ThreePIdTwo) x -= TwoPI;
	else if(x > HalfPI) { x -= PI; ChangeSign = 1;}

	double xe2 = x*x;
	Cos = float(1. + xe2*(a2c + xe2*(a4c + xe2*(a6c + xe2*(a8c + xe2*a10c)))));
	Sin = float(x*(1. + xe2*(a3s + xe2*(a5s + xe2*(a7s + xe2*(a9s + xe2*a11s))))));
	if(ChangeSign) { Cos = -Cos; Sin = -Sin;}
}

//*************************************************************************

inline void srTGenOptElem::FindLowestAndUppestPoints(TVector3d& Dir, TVector3d* Points, int AmOfPoints, int& LowestInd, int& UppestInd)
{
	TVector3d LowestPo = *Points, UppestPo = *Points;
	LowestInd = UppestInd = 0;
	for(int i=1; i<AmOfPoints; i++)
	{
		TVector3d TestLoV = Points[i] - LowestPo;
		TVector3d TestUpV = Points[i] - UppestPo;
		if(TestLoV*Dir < 0.) { LowestPo = Points[i]; LowestInd = i;}
		if(TestUpV*Dir > 0.) { UppestPo = Points[i]; UppestInd = i;}
	}
}

//*************************************************************************

inline void srTGenOptElem::ReflectVect(TVector3d& N, TVector3d& V)
{
	V = V - ((2.*(N*V))*N);
}

//*************************************************************************

inline void srTGenOptElem::FindLineIntersectWithPlane(TVector3d* Plane, TVector3d* Line, TVector3d& IntersectP)
{// Assumes N, V unit vectors !!!
 // This does not distinguish Zero Intersection and Line belonging to the Plane !!!
	TVector3d &Rp0 = *Plane, &N = Plane[1];
	TVector3d &Rl0 = *Line, &V = Line[1];
	const double RelTol = 1.E-12; //1.E-10;
	double VN = V*N;
	//if(::fabs(VN) > RelTol) IntersectP = Rl0 + ((Rp0 - Rl0)*N/VN)*V;
	if(::fabs(VN) > RelTol) 
	{
		double inv_VN = 1./VN;
		IntersectP = Rl0 + ((inv_VN*((Rp0 - Rl0)*N))*V);
	}
	//if(::fabs(VN) > RelTol) IntersectP = Plus(Rl0, ((Rp0 - Rl0)*N/VN)*V);
	else IntersectP = Rl0;
}

//*************************************************************************

inline void srTGenOptElem::TreatPhaseShift(srTEFieldPtrs& EPtrs, double PhShift)
{
	float CosPh, SinPh;
	CosAndSin(PhShift, CosPh, SinPh);
	float ExReNew = (*EPtrs.pExRe)*CosPh - (*EPtrs.pExIm)*SinPh;
	float ExImNew = (*EPtrs.pExRe)*SinPh + (*EPtrs.pExIm)*CosPh;
	float EzReNew = (*EPtrs.pEzRe)*CosPh - (*EPtrs.pEzIm)*SinPh;
	float EzImNew = (*EPtrs.pEzRe)*SinPh + (*EPtrs.pEzIm)*CosPh;
	*EPtrs.pExRe = ExReNew; *EPtrs.pExIm = ExImNew;
	*EPtrs.pEzRe = EzReNew; *EPtrs.pEzIm = EzImNew;
}

//*************************************************************************

//inline void srTGenOptElem::SetupExpCorrArray(float* pCmpData, long AmOfPt, double x, double qStart, double qStep)
inline void srTGenOptElem::SetupExpCorrArray(float* pCmpData, long long AmOfPt, double x, double qStart, double qStep)
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

//*************************************************************************

inline void srTGenOptElem::SetupRadXorZSectFromSliceConstE(float* pInEx, float* pInEz, long nx, long nz, char vsX_or_vsZ, long iSect, float* pOutEx, float* pOutEz)
{
	//long Per = (vsX_or_vsZ == 'x')? 2 : (nx << 1);
	long long Per = (vsX_or_vsZ == 'x')? 2 : (nx << 1);
	float *tOutEx = pOutEx, *tOutEz = pOutEz;
	//long StartOffset = (vsX_or_vsZ == 'x')? iSect*(nx << 1) : (iSect << 1);
	long long StartOffset = (vsX_or_vsZ == 'x')? iSect*(nx << 1) : (iSect << 1);
	float *tEx = pInEx + StartOffset, *tEz = pInEz + StartOffset;
	//long nPt = (vsX_or_vsZ == 'x')? nx : nz;
	long long nPt = (vsX_or_vsZ == 'x')? nx : nz;

	//for(int i=0; i<nPt; i++)
	for(long long i=0; i<nPt; i++)
	{
		*(tOutEx++) = *tEx; *(tOutEx++) = *(tEx + 1);
		*(tOutEz++) = *tEz; *(tOutEz++) = *(tEz + 1);
		tEx += Per; tEz += Per;
	}
}

//*************************************************************************

inline float srTGenOptElem::IntegrateElField1D(srTRadSect1D& Sect1D, char TreatExOrEz, float& IntegOnHalfInterv)
{
	float* tVal = (TreatExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz;
	//long np_mi_1 = Sect1D.np - 1;
	//long Half_np_mi_1 = np_mi_1 >> 1;
	long long np_mi_1 = Sect1D.np - 1;
	long long Half_np_mi_1 = np_mi_1 >> 1;
	float Sum = (float)(0.5*((*tVal) + (*(tVal + (np_mi_1 << 1)))));
	float SumHalf = (float)(0.5*((*tVal) + np_mi_1));
	//for(long i=1; i<np_mi_1; i++) 
	for(long long i=1; i<np_mi_1; i++) 
	{ 
		Sum += *tVal;
		if(i < Half_np_mi_1) SumHalf += *tVal;
		tVal += 2;
	}
	IntegOnHalfInterv = (float)(SumHalf*Sect1D.ArgStep);
	return (float)(Sum*Sect1D.ArgStep);
}

//*************************************************************************

inline long srTGenOptElem::IntegerOffsetCoord(double xStart, double xStep, double xVal)
{
	long iCoor = long((xVal - xStart)/xStep);
	if(::fabs(xVal - ((iCoor + 1)*xStep + xStart)) < 1.E-05*xStep) iCoor++;
	return iCoor;
}

//*************************************************************************

inline char srTGenOptElem::ChooseTreatExOrEzBasedOnMax(srTRadSect1D& Sect1D)
{
	float MaxAbsEx, MaxAbsEz;
	//long IndMaxAbsEx, IndMaxAbsEz;
	long long IndMaxAbsEx, IndMaxAbsEz;
	FindMaximumAbsReE(Sect1D, MaxAbsEx, IndMaxAbsEx, MaxAbsEz, IndMaxAbsEz);
	return (MaxAbsEx > MaxAbsEz)? 'x' : 'z';
}

//*************************************************************************

inline int srTGenOptElem::ResizePropagateAndAnalizeFringes1D(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam, char TreatExOrEz, srTFringeInfo& PropFringeInfo)
{// Gives PropFringeInfo
	int result;
	srTRadSect1D SectDupl;
	if(result = Sect1D.SetupDupl(SectDupl)) return result;
	if((ResizeParam.pm < 1.) || (ResizeParam.pd < 1.))
	{
		if(result = PropagateRadiationSimple1D(&SectDupl)) return result; // To implement in all elements
		if(result = RadResizeGen1D(SectDupl, ResizeParam)) return result;
	}
	else
	{
		if(result = RadResizeGen1D(SectDupl, ResizeParam)) return result;
		if(result = PropagateRadiationSimple1D(&SectDupl)) return result; // To implement in all elements
	}
	return AnalizeFringes(SectDupl, TreatExOrEz, PropFringeInfo);
}

//*************************************************************************

inline int srTGenOptElem::CheckIfOversamplingIsReal(srTRadSect1D& Sect1D, srTRadResize1D& ResizeParam, char TreatExOrEz, srTFringeInfo& PropFringeInfo)
{// Makes further oversampling and check that the FringeInfo is qualitatively the same after that
	const double RelPointsPerFringeMisfit = 0.15; // To steer
	const double PdIncrCoef = 1.6; // To steer

	const double TolRelNumberOfFringes = 0.01; // To steer
	const int MinTolAbsNumberOfFringes = 3; // To steer

	char LeftIsLarge = (PropFringeInfo.LeftPointsPerFringe > (1. + RelPointsPerFringeMisfit));
	char RightIsLarge = (PropFringeInfo.RightPointsPerFringe > (1. + RelPointsPerFringeMisfit));
	if(!(LeftIsLarge || RightIsLarge)) return 0;

	int TolAbsNumberOfFringes = (int)(TolRelNumberOfFringes*PropFringeInfo.AmOfFringes);
	if(TolAbsNumberOfFringes < MinTolAbsNumberOfFringes) TolAbsNumberOfFringes = MinTolAbsNumberOfFringes;

	int result;
	srTRadResize1D ResizeParamLoc = ResizeParam;
	ResizeParamLoc.pd *= PdIncrCoef;
	srTFringeInfo PropFringeInfoLoc;
	if(result = ResizePropagateAndAnalizeFringes1D(Sect1D, ResizeParamLoc, TreatExOrEz, PropFringeInfoLoc)) return result;

	char NumFringesIsBad = ((::fabs((long double)(PropFringeInfoLoc.AmOfFringes - PropFringeInfo.AmOfFringes))) > TolAbsNumberOfFringes);

	if((LeftIsLarge && (PropFringeInfoLoc.LeftPointsPerFringe*(1. + RelPointsPerFringeMisfit) < PropFringeInfo.LeftPointsPerFringe)) ||
		NumFringesIsBad)
	{
		PropFringeInfo.LeftPointsPerFringe = 1.;
		//double AuxPointsPerFringe = (PropFringeInfoLoc.LeftPointsPerFringe)/PdIncrCoef;
		//if(AuxPointsPerFringe < 1.) AuxPointsPerFringe = 1.;
		//if(AuxPointsPerFringe < PropFringeInfo.LeftPointsPerFringe)
		//{
  //          PropFringeInfo.LeftPointsPerFringe = AuxPointsPerFringe;
		//}
	}
	if((RightIsLarge && (PropFringeInfoLoc.RightPointsPerFringe*(1. + RelPointsPerFringeMisfit) < PropFringeInfo.RightPointsPerFringe)) ||
	   NumFringesIsBad)
	{
		PropFringeInfo.RightPointsPerFringe = 1.;
		//double AuxPointsPerFringe = (PropFringeInfoLoc.RightPointsPerFringe)/PdIncrCoef;
		//if(AuxPointsPerFringe < 1.) AuxPointsPerFringe = 1.;
		//if(AuxPointsPerFringe < PropFringeInfo.RightPointsPerFringe)
		//{
  //          PropFringeInfo.RightPointsPerFringe = AuxPointsPerFringe;
		//}
	}
	return 0;
}

//*************************************************************************

inline double srTGenOptElem::SuggestResolResizeCoef(double PoPerFringeBeforePropag, double PoPerFringeAfterPropag, long AmOfFringesAfterPropag)
{
	const double TransferCoef = 0.4;//0.3; // To choose next value of pd // To steer
	//const int MinAmOfFringesToAllowTuning = 20; //This caused problems for Far IR (JLab)

	double pdLoc = PoPerFringeBeforePropag/PoPerFringeAfterPropag;
	if(pdLoc > 1.) 
	{
		pdLoc = 1. + TransferCoef*(pdLoc - 1.);
	}
	else if(pdLoc < 1.)
	{
		//if(AmOfFringesAfterPropag < MinAmOfFringesToAllowTuning)
		//{
		//	pdLoc = 1.;
		//}
		//else
		//{
        pdLoc = 1. - TransferCoef*(1. - pdLoc);
		//}
	}
	return pdLoc;
}

//*************************************************************************

inline int srTGenOptElem::MemoryIsSufficientForProp(srTSRWRadStructAccessData* pRadAccessData, srTPropagScenario1D* PropagScenario, double& ExtraSizeBeforeProp, double& ExtraSizeAfterProp)
{
	EstimateMemoryNeededForPropag(pRadAccessData, PropagScenario, ExtraSizeBeforeProp, ExtraSizeAfterProp); // in bytes
	double MemNeededForPropag = (ExtraSizeBeforeProp > ExtraSizeAfterProp)? ExtraSizeBeforeProp : ExtraSizeAfterProp;
	double MemAvail = CheckMemoryAvailable();
	return (MemAvail > MemNeededForPropag);
}

//*************************************************************************

inline void srTGenOptElem::FindCenterWithRespectToThresh(srTRadSect1D& Sect1D, double AbsThresh, char TreatExOrEz, long& ic, double& xc)
{
	long iLeftThreshBorder, iRightThreshBorder;
	FindThresholdBorders(Sect1D, AbsThresh, TreatExOrEz, iLeftThreshBorder, iRightThreshBorder);
	ic = iLeftThreshBorder + ((iRightThreshBorder - iLeftThreshBorder) >> 1);
	xc = Sect1D.ArgStart + ic*Sect1D.ArgStep;
}

//*************************************************************************

inline void srTGenOptElem::RejectSmallResize(srTRadResize& Resize)
{
	const double RelDifNotResize = 0.01; //0.05; // To steer
	if(::fabs(Resize.pxm - 1.) < RelDifNotResize) Resize.pxm = 1.;
	if(::fabs(Resize.pxd - 1.) < RelDifNotResize) Resize.pxd = 1.;
	if(::fabs(Resize.pzm - 1.) < RelDifNotResize) Resize.pzm = 1.;
	if(::fabs(Resize.pzd - 1.) < RelDifNotResize) Resize.pzd = 1.;
}

//*************************************************************************

inline double srTGenOptElem::MaxIntInHorString(long iz, srTSRWRadStructAccessData* pRadAccessData)
{
	double IMax = 0.;
	//long Offset = ((pRadAccessData->nx) << 1)*iz;
	long long Offset = ((pRadAccessData->nx) << 1)*iz;
	float *tEx = pRadAccessData->pBaseRadX + Offset, *tEz = pRadAccessData->pBaseRadZ + Offset;
	for(long ix=0; ix<pRadAccessData->nx; ix++)
	{
		double ReEx = *(tEx++), ReEz = *(tEz++);
		double ImEx = *(tEx++), ImEz = *(tEz++);
		double I = ReEx*ReEx + ImEx*ImEx + ReEz*ReEz + ImEz*ImEz;
		if(IMax < I) IMax = I;
	}
	return IMax;
}

//*************************************************************************

inline double srTGenOptElem::MaxIntInVertString(long ix, long izSt, long izFi, srTSRWRadStructAccessData* pRadAccessData)
{
	double IMax = 0.;
	//long Per = (pRadAccessData->nx) << 1;
	//long Offset = Per*izSt + (ix << 1);
	long long Per = (pRadAccessData->nx) << 1;
	long long Offset = Per*izSt + (ix << 1);
	float *tEx = pRadAccessData->pBaseRadX + Offset, *tEz = pRadAccessData->pBaseRadZ + Offset;
	for(long iz=izSt; iz<=izFi; iz++)
	{
		double ReEx = *tEx, ReEz = *tEz;
		double ImEx = *(tEx+1), ImEz = *(tEz+1);
		double I = ReEx*ReEx + ImEx*ImEx + ReEz*ReEz + ImEz*ImEz;
		if(IMax < I) IMax = I;
		tEx += Per; tEz += Per;
	}
	return IMax;
}

//*************************************************************************

inline double srTGenOptElem::FindLongCoordOfHorBaseVect(double vLx, double vLy, double vLz, double vHx, double vHy)
{//Used for processing the case of manual definition of output beam frame for some optical elements
	TVector3d vL(vLx, vLy, vLz);
	vL.Normalize();

	TVector3d vH(vHx, vHy, 0.);
	if(vLz == 0.)
	{
		if((vHx != 0.) || (vHy != 0.))
		{
			vH.Normalize();
			const double relTol = 1.e-12;
			double testScalProd = ::fabs(vL*vH);
			if(testScalProd > relTol) throw FAILED_DETERMINE_OPTICAL_AXIS;
			return 0.;
		}
		else return 1.;
	}
	else
	{
		return (-vL.x*vH.x - vL.y*vH.y)/vL.z;
	}
}

//*************************************************************************

//class srTOptElemSummary {
//public:
//	static int SetupOpticalElement(srTStringVect*, srTGenOptElemHndl&, srTSRWRadStructAccessData*);
//};

//*************************************************************************

#endif
