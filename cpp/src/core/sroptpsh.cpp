/************************************************************************//**
 * File: sroptpsh.cpp
 * Description: Optical element: Phase Shift (obsolete?)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptpsh.h"
#include "gmfft.h"
#include "srradmnp.h"
#include "gmfit.h"

//*************************************************************************

int srTPhaseShift::SetUpPointSourceSect1D(srTRadSect1D* PointSourceSect1D, srTSRWRadStructAccessData* pRad)
{
	const long SmallestN = 10;
	const long OverSampFact = 8;
	CGenMathFFT1D FFT;

	double ePh_keV = pRad->eStart;
	if(pRad->PhotEnergyUnit == 0) ePh_keV *= 0.001;
	double WavelengthIn_m = 1.239854*1.E-09/ePh_keV;

	double HalfLambRx = 0.5*WavelengthIn_m*pRad->RobsX;
	double xStartRel = pRad->xStart - TransvCenPoint.x;
	double xEndRel = pRad->xStart + (pRad->nx - 1)*pRad->xStep - TransvCenPoint.x;
	double dxStart = ::fabs(HalfLambRx/xStartRel);
	double dxEnd = ::fabs(HalfLambRx/xEndRel);
	double dx = (dxStart < dxEnd)? dxStart : dxEnd;
	long Nx = (long(::fabs(xEndRel - xStartRel)/dx) + 1)*OverSampFact;
	if(((Nx>>1)<<1) != Nx) Nx++;
	FFT.NextCorrectNumberForFFT(Nx);
	if(Nx < SmallestN) Nx = SmallestN;

	double HalfLambRz = 0.5*WavelengthIn_m*pRad->RobsZ;
	double zStartRel = pRad->zStart - TransvCenPoint.y;
	double zEndRel = pRad->zStart + (pRad->nz - 1)*pRad->zStep - TransvCenPoint.y;
	double dzStart = ::fabs(HalfLambRz/zStartRel);
	double dzEnd = ::fabs(HalfLambRz/zEndRel);
	double dz = (dzStart < dzEnd)? dzStart : dzEnd;
	long Nz = (long(::fabs(zEndRel - zStartRel)/dz) + 1)*OverSampFact;
	if(((Nz>>1)<<1) != Nz) Nz++;
	FFT.NextCorrectNumberForFFT(Nz);
	if(Nz < SmallestN) Nz = SmallestN;

	PointSourceSect1D[0].pEx = new float[Nx << 1];
	if(PointSourceSect1D[0].pEx == 0) return MEMORY_ALLOCATION_FAILURE;
	PointSourceSect1D[1].pEx = new float[Nz << 1];
	if(PointSourceSect1D[1].pEx == 0) 
	{
		delete[] PointSourceSect1D[0].pEx;
		return MEMORY_ALLOCATION_FAILURE;
	}

	PointSourceSect1D[0].ArgStep = (xEndRel - xStartRel)/(Nx - 1);
	PointSourceSect1D[0].ArgStart = xStartRel;
	PointSourceSect1D[0].np = Nx;
	PointSourceSect1D[0].VsXorZ = 'x';
	PointSourceSect1D[0].Robs = pRad->RobsX;
	PointSourceSect1D[0].RobsAbsErr = pRad->RobsXAbsErr;
	char NameWaveX[] = "AuxOptCompSetupX";
	strcpy(PointSourceSect1D[0].NameOfWave, NameWaveX);

	PointSourceSect1D[1].ArgStep = (zEndRel - zStartRel)/(Nz - 1);
	PointSourceSect1D[1].ArgStart = zStartRel;
	PointSourceSect1D[1].np = Nz;
	PointSourceSect1D[1].VsXorZ = 'z';
	PointSourceSect1D[1].Robs = pRad->RobsZ;
	PointSourceSect1D[1].RobsAbsErr = pRad->RobsZAbsErr;
	char NameWaveZ[] = "AuxOptCompSetupZ";
	strcpy(PointSourceSect1D[1].NameOfWave, NameWaveZ);

	for(int k=0; k<2; k++)
	{
		PointSourceSect1D[k].pEz = 0;
		PointSourceSect1D[k].eVal = pRad->eStart;
		PointSourceSect1D[k].OtherCoordVal = 0.;
		PointSourceSect1D[k].WfrEdgeCorrShouldBeDone = 0;
		PointSourceSect1D[k].Pres = 0;
		PointSourceSect1D[k].PhotEnergyUnit = pRad->PhotEnergyUnit;
		PointSourceSect1D[k].DeleteArraysAtDestruction = 1;
	}

	const double TwoPI = 6.2831853071796;
	double TwoPi_d_Lamb = TwoPI/WavelengthIn_m;

	double Rx = pRad->RobsX;
	double RxRx = Rx*Rx;
	double xStepRel = (xEndRel - xStartRel)/(Nx - 1);
	double xx = xStartRel;
	float *tEx = PointSourceSect1D[0].pEx;
	for(long ix=0; ix<Nx; ix++)
	{
		double Ph = TwoPi_d_Lamb*sqrt(RxRx + xx*xx);
		CosAndSin(Ph, *tEx, *(tEx+1));
		tEx += 2;
		xx += xStepRel;
	}

	double Rz = pRad->RobsZ;
	double RzRz = Rz*Rz;
	double zStepRel = (zEndRel - zStartRel)/(Nz - 1);
	double zz = zStartRel;
	tEx = PointSourceSect1D[1].pEx;
	for(long iz=0; iz<Nz; iz++)
	{
		double Ph = TwoPi_d_Lamb*sqrt(RzRz + zz*zz);
		CosAndSin(Ph, *tEx, *(tEx+1));
		tEx += 2;
		zz += zStepRel;
	}

	return 0;
}

//*************************************************************************

int srTPhaseShift::SetUpPhaseShiftWave1D(srTRadSect1D& Sect1D, srTWaveAccessData& PhShData1D)
{
	char NameOfPhShWaveX[] = "AuxPhaseShiftWaveX";
	char NameOfPhShWaveZ[] = "AuxPhaseShiftWaveZ";

	PhShData1D.pWaveData = 0;
	PhShData1D.WaveType[0] = 'd';
	PhShData1D.AmOfDims = 2;
	if(Sect1D.VsXorZ == 'x')
	{
		PhShData1D.DimSizes[0] = Sect1D.np;
		PhShData1D.DimStartValues[0] = Sect1D.ArgStart;
		PhShData1D.DimSteps[0] = Sect1D.ArgStep;

		PhShData1D.DimSizes[1] = 1;
		PhShData1D.DimStartValues[1] = Sect1D.OtherCoordVal;
		PhShData1D.DimSteps[1] = 1.E-06; // Dummy

		strcpy(PhShData1D.NameOfWave, NameOfPhShWaveX);
	}
	else
	{
		PhShData1D.DimSizes[0] = 1;
		PhShData1D.DimStartValues[0] = Sect1D.OtherCoordVal;
		PhShData1D.DimSteps[0] = 1.E-06; // Dummy

		PhShData1D.DimSizes[1] = Sect1D.np;
		PhShData1D.DimStartValues[1] = Sect1D.ArgStart;
		PhShData1D.DimSteps[1] = Sect1D.ArgStep;

		strcpy(PhShData1D.NameOfWave, NameOfPhShWaveZ);
	}

	char* UnitName = PhShData1D.DimUnits[0];
	UnitName[0] = 'm'; UnitName[1] = '\0';

	//srTSend Send;
	//int result = Send.MakeWaveAccordingToWaveAccessData(PhShData1D);
	//if(result) return result;
	//return Send.SetUpPhaseShiftWave(PhShData1D, FunNo);
	//?????
	return 0;
}

//*************************************************************************

int srTPhaseShift::SetUpPhaseShiftWave(srTSRWRadStructAccessData& RadData, srTWaveAccessData& PhShData)
{
	PhShData.pWaveData = 0;
	PhShData.WaveType[0] = 'd';
	PhShData.AmOfDims = 2;
	
	PhShData.DimSizes[0] = RadData.nx;
	PhShData.DimSizes[1] = RadData.nz;

	PhShData.DimStartValues[0] = RadData.xStart;
	PhShData.DimStartValues[1] = RadData.zStart;

	PhShData.DimSteps[0] = RadData.xStep;
	PhShData.DimSteps[1] = RadData.zStep;

	for(int k=0; k<2; k++)
	{
		char* UnitName = PhShData.DimUnits[k];
		UnitName[0] = 'm'; UnitName[1] = '\0';
	}

	char NameOfPhShWave[] = "AuxPhaseShiftWave";
	strcpy(PhShData.NameOfWave, NameOfPhShWave);

	//srTSend Send;
	//int result = Send.MakeWaveAccordingToWaveAccessData(PhShData);
	//if(result) return result;
	//return Send.SetUpPhaseShiftWave(PhShData, FunNo);
	//?????
	return 0;
}

//*************************************************************************

int srTPhaseShift::SetUpFocalDistByPropag1D(srTRadSect1D& Sect1D)
{
	int result;
	if(result = PropagateRadiationSimple1D(&Sect1D)) return result;

	double* PhaseCont = new double[Sect1D.np];
	if(PhaseCont == 0) return MEMORY_ALLOCATION_FAILURE;

	float *tE = Sect1D.pEx;
	double *tPh = PhaseCont;
	for(long i=0; i<Sect1D.np; i++)
	{
		*(tPh++) = FormalPhase(*tE, *(tE + 1));
		tE += 2;
	}

	srTRadGenManip RadManip;
	RadManip.TryToMakePhaseContinuous1D(PhaseCont, Sect1D.np, -1, 0.);

	float* PhaseContF = new float[Sect1D.np + 1];
	if(PhaseContF == 0) return MEMORY_ALLOCATION_FAILURE;
	float* ArgCont = new float[Sect1D.np + 1];
	if(ArgCont == 0) return MEMORY_ALLOCATION_FAILURE;
	float* Sigm = new float[Sect1D.np + 1];
	if(Sigm == 0) return MEMORY_ALLOCATION_FAILURE;
	const float AbsErrorLevel = (float)0.01; // To steer

	float Arg = (float)Sect1D.ArgStart, Step = (float)Sect1D.ArgStep;
	for(long j=1; j<=Sect1D.np; j++)
	{
		PhaseContF[j] = (float)PhaseCont[j - 1];
		Sigm[j] = AbsErrorLevel;
		ArgCont[j] = Arg;
		Arg += Step;
	}

	float a[4], chisq, qOK;
	int ia[] = {1,1,1,1};
	CGenMathFit Fit;
	if(result = Fit.FitPolynomial(ArgCont, PhaseContF, Sigm, int(Sect1D.np), a, ia, 3, &chisq, &qOK)) return result;

	const double PI = 3.141592653590;
	double ePh_keV = Sect1D.eVal;
	if(Sect1D.PhotEnergyUnit == 0) ePh_keV *= 0.001;
	double WavelengthIn_m = 1.239854*1.E-09/ePh_keV;
	double Rafter = PI/(a[3]*WavelengthIn_m);
	double Rbefore = Sect1D.Robs;

	const double RelTolR = 0.05; // To steer
	if(::fabs(Rbefore - Rafter) > RelTolR*(::fabs(Rbefore)))
	{
		IsFocusing = 1;
		double FocDist = Rbefore*Rafter/(Rafter - Rbefore);
		if(Sect1D.VsXorZ == 'x') FocDistX = FocDist;
		else if(Sect1D.VsXorZ == 'z') FocDistZ = FocDist;
	}

	if(PhaseCont != 0) delete[] PhaseCont;
	if(PhaseContF != 0) delete[] PhaseContF;
	if(ArgCont != 0) delete[] ArgCont;
	if(Sigm != 0) delete[] Sigm;
	return 0;
}

//*************************************************************************
