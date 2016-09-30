/************************************************************************//**
 * File: sroptgtr.cpp
 * Description: Optical element: Transmission
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptgtr.h"
#include "gmfft.h"
#include "gmfit.h"
#include "srradmnp.h"
#include "srerror.h"
#include "srwlib.h"

//*************************************************************************

srTGenTransmission::srTGenTransmission(srTStringVect* pElemInfo, srTDataMD* pExtraData) 
{
	ErrorCode = 0;
	GenTransNumData.pData = 0;

	char NumStructName[256];
	strcpy(NumStructName, (*pElemInfo)[1]);

	//srTSend Send;
	//ErrorCode = Send.FetchNumWave(NumStructName, &GenTransNumData);

	if(pExtraData != 0)
	{
		GenTransNumData = *pExtraData;
	}
	if(GenTransNumData.pData == 0) ErrorCode = NT_FP64_COMPLEX_WAVE_REQUIRED;

	if(ErrorCode) return;
	if((*(GenTransNumData.DataType) != 'c') || (*(GenTransNumData.DataType + 1) != 'd'))
	{
		ErrorCode = NT_FP64_COMPLEX_WAVE_REQUIRED; return;
	}

	char BufStr[256];

	int IndOptPathOrPhase = 10;
	//OptPathOrPhase = 2; // phase shift by default
	OptPathOrPhase = 1; // OptPath by default
	int LenElemInfo = (int)pElemInfo->size();
	if(LenElemInfo > IndOptPathOrPhase)
	{
		strcpy(BufStr, (*pElemInfo)[IndOptPathOrPhase]);
		OptPathOrPhase = (char)atof(BufStr);
		if((OptPathOrPhase <= 0) || (OptPathOrPhase > 2)) OptPathOrPhase = 1;
	}

	strcpy(BufStr, (*pElemInfo)[3]);
	eMid = 0.;
	double Aux_eMid = atof(BufStr);
	if(Aux_eMid > 0) eMid = Aux_eMid;
	//if(eMid == 0.) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}

	TransvCenPoint.x = 0; //???
	strcpy(BufStr, (*pElemInfo)[4]);
	double aux_xc = atof(BufStr);
	if(::fabs(aux_xc) < 1.e+10) TransvCenPoint.x = aux_xc;

	TransvCenPoint.y = 0; //???
	strcpy(BufStr, (*pElemInfo)[5]);
	double aux_yc = atof(BufStr);
	if(::fabs(aux_yc) < 1.e+10) TransvCenPoint.y = aux_yc;

	strcpy(BufStr, (*pElemInfo)[6]);
	OuterTransmIs = atoi(BufStr); // Program this at propagation!

	strcpy(BufStr, (*pElemInfo)[7]); // Setup was completed or not
	int SetupIsCompleted = atoi(BufStr);
	if(SetupIsCompleted) 
	{
		strcpy(BufStr, (*pElemInfo)[8]);
		FocDistX = atof(BufStr);
		if(FocDistX == 0.) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}

		strcpy(BufStr, (*pElemInfo)[9]);
		FocDistZ = atof(BufStr);
		if(FocDistZ == 0.) { ErrorCode = IMPROPER_OPTICAL_COMPONENT_STRUCTURE; return;}
	}
	else
	{
		if(ErrorCode = EstimateFocalDistancesAndCheckSampling()) return;

		//"erase" part of existing strings
		char *aStr=0;
		int AuxInfStartInd = 7;
		for(int k=AuxInfStartInd; k<(int)(pElemInfo->size()); k++)
		{
			aStr = (*pElemInfo)[k];
			//if(aStr != 0) delete[] aStr;
			if(aStr != 0) *aStr = '\0';
		}
		//pElemInfo->erase(pElemInfo->begin() + AuxInfStartInd, pElemInfo->end());

		//fill-in the strings with new information
		//aStr = new char[256]; 
		//if(aStr == 0) { ErrorCode = MEMORY_ALLOCATION_FAILURE; return;}
		aStr = (*pElemInfo)[AuxInfStartInd];
		sprintf(aStr, "1");
		//pElemInfo->push_back(aStr);

		//aStr = new char[256]; 
		//if(aStr == 0) { ErrorCode = MEMORY_ALLOCATION_FAILURE; return;}
		aStr = (*pElemInfo)[AuxInfStartInd + 1];
		sprintf(aStr, "%g", FocDistX);
		//pElemInfo->push_back(aStr);

		//aStr = new char[256]; 
		//if(aStr == 0) { ErrorCode = MEMORY_ALLOCATION_FAILURE; return;}
		aStr = (*pElemInfo)[AuxInfStartInd + 2];
		sprintf(aStr, "%g", FocDistZ);
		//pElemInfo->push_back(aStr);

		//aStr = new char[256]; 
		//if(aStr == 0) { ErrorCode = MEMORY_ALLOCATION_FAILURE; return;}
		aStr = (*pElemInfo)[AuxInfStartInd + 3];
		sprintf(aStr, "%d", OptPathOrPhase);
		//pElemInfo->push_back(aStr);
	}

	if(ErrorCode = EstimateMinimalContinuousIntervals()) return;
}

//*************************************************************************

srTGenTransmission::srTGenTransmission(const SRWLOptT& tr)
{
	OptPathOrPhase = 1; //opt. path dif.
	OuterTransmIs = tr.extTr + 1;
	
	//eMid = 0;
	eMid = 0.5*(tr.mesh.eStart + tr.mesh.eFin);

	//TransvCenPoint.x = tr.x;
	//TransvCenPoint.y = tr.y;
	TransvCenPoint.x = 0.5*(tr.mesh.xStart + tr.mesh.xFin);
	TransvCenPoint.y = 0.5*(tr.mesh.yStart + tr.mesh.yFin);

	FocDistX = tr.Fx;
	FocDistZ = tr.Fy;

	GenTransNumData.pData = (char*)(tr.arTr);
	GenTransNumData.DataType[0] = 'c';
	GenTransNumData.DataType[1] = 'f';
	//GenTransNumData.AmOfDims = 2;
	GenTransNumData.AmOfDims = 3; //OC112312
	
	//GenTransNumData.DimSizes[0] = tr.nx;
	//GenTransNumData.DimSizes[1] = tr.ny;
	GenTransNumData.DimSizes[0] = tr.mesh.ne; //OC112312
	GenTransNumData.DimSizes[1] = tr.mesh.nx;
	GenTransNumData.DimSizes[2] = tr.mesh.ny;

	//GenTransNumData.DimStartValues[0] = tr.x - 0.5*tr.rx;
	//GenTransNumData.DimStartValues[1] = tr.y - 0.5*tr.ry;
	//GenTransNumData.DimSteps[0] = (tr.nx > 1)? tr.rx/(tr.nx - 1) : 0;
	//GenTransNumData.DimSteps[1] = (tr.ny > 1)? tr.ry/(tr.ny - 1) : 0;

	GenTransNumData.DimStartValues[0] = tr.mesh.eStart; //OC112312
	GenTransNumData.DimStartValues[1] = tr.mesh.xStart;
	GenTransNumData.DimStartValues[2] = tr.mesh.yStart;
	GenTransNumData.DimSteps[0] = (tr.mesh.ne > 1)? (tr.mesh.eFin - tr.mesh.eStart)/(tr.mesh.ne - 1) : 0;
	GenTransNumData.DimSteps[1] = (tr.mesh.nx > 1)? (tr.mesh.xFin - tr.mesh.xStart)/(tr.mesh.nx - 1) : 0;
	GenTransNumData.DimSteps[2] = (tr.mesh.ny > 1)? (tr.mesh.yFin - tr.mesh.yStart)/(tr.mesh.ny - 1) : 0;

	//GenTransNumData.DimUnits[10][255];
	//GenTransNumData.DataUnits[255];
	//GenTransNumData.DataName[255];
	//GenTransNumData.hState; //auxiliary
}

//*************************************************************************

int srTGenTransmission::EstimateMinimalContinuousIntervals()
{//	DxContin, DzContin
	const double RelTolForDiscont = 0.6; // To steer
	const double HalfRelTol = 0.5*RelTolForDiscont;

	long Ne=1, Nx=1, Nz=1;
	double xStep=0, zStep=0;

	if(GenTransNumData.AmOfDims == 2)
	{
		Nx = (GenTransNumData.DimSizes)[0]; 
		Nz = (GenTransNumData.DimSizes)[1];
		//double xStart = (GenTransNumData.DimStartValues)[0], zStart = (GenTransNumData.DimStartValues)[1];
		xStep = (GenTransNumData.DimSteps)[0]; 
		zStep = (GenTransNumData.DimSteps)[1];
		//double xRange = (Nx - 1)*xStep, zRange = (Nz - 1)*zStep;
	}
	else if(GenTransNumData.AmOfDims == 3)
	{
		Ne = (GenTransNumData.DimSizes)[0]; 
		Nx = (GenTransNumData.DimSizes)[1]; 
		Nz = (GenTransNumData.DimSizes)[2];

		xStep = (GenTransNumData.DimSteps)[1]; 
		zStep = (GenTransNumData.DimSteps)[2];
	}

	int iDxContinT = Nx - 1, iDxContinP = Nx - 1;
	int iDzContinT = Nz - 1, iDzContinP = Nz - 1;

	//long zPer = Nx << 1;
	//long xPer = Ne << 1;
	//long zPer = xPer*Nx;
	long long xPer = Ne << 1;
	long long zPer = xPer*Nx;

	DOUBLE *pT0 = (DOUBLE*)(GenTransNumData.pData);
	if(pT0 == 0) return IMPROPER_OPTICAL_COMPONENT_STRUCTURE;
	//DOUBLE *pP0 = (DOUBLE*)(GenTransNumData.pData) + 1;

	int *Aux_izStartDiscontT = new int[Nx];
	if(Aux_izStartDiscontT == 0) return MEMORY_ALLOCATION_FAILURE;
	int *Aux_izStartDiscontP = new int[Nx];
	if(Aux_izStartDiscontP == 0) return MEMORY_ALLOCATION_FAILURE;

	for(int ie=0; ie<Ne; ie++)
	{
		int *tAuxT = Aux_izStartDiscontT, *tAuxP = Aux_izStartDiscontP;
		for(int ix=0; ix<Nx; ix++)
		{
			*(tAuxT++) = 0; *(tAuxP++) = 0;
		}

		for(int iz=1; iz<(Nz - 1); iz++)
		{
			int ixStartDiscontT = 0, ixStartDiscontP = 0;
			for(int ix=1; ix<(Nx - 1); ix++)
			{
				//DOUBLE *pT = pT0 + iz*zPer + (ix << 1);
				DOUBLE *pT = pT0 + iz*zPer + ix*xPer + (ie << 1);
				DOUBLE *pP = pT + 1;

				double zPrevT = *(pT - zPer), xPrevT = *(pT - 2);
				double zPrevP = *(pP - zPer), xPrevP = *(pP - 2);

				double xLeftDifT = *pT - xPrevT, xRightDifT = *(pT + 2) - *pT;
				double zLeftDifT = *pT - zPrevT, zRightDifT = *(pT + zPer) - *pT;

				double xLeftDifP = *pP - xPrevP, xRightDifP = *(pP + 2) - *pP;
				double zLeftDifP = *pP - zPrevP, zRightDifP = *(pP + zPer) - *pP;

				if(::fabs(xRightDifT - xLeftDifT) > ::fabs(xLeftDifT + xRightDifT)*HalfRelTol) // Discontinuity
				{
					int Cur_iDxContinT = ix - ixStartDiscontT;
					if(iDxContinT > Cur_iDxContinT) iDxContinT = Cur_iDxContinT;
					ixStartDiscontT = ix;
				}
				if(::fabs(xRightDifP - xLeftDifP) > ::fabs(xLeftDifP + xRightDifP)*HalfRelTol) // Discontinuity
				{
					int Cur_iDxContinP = ix - ixStartDiscontP;
					if(iDxContinP > Cur_iDxContinP) iDxContinP = Cur_iDxContinP;
					ixStartDiscontP = ix;
				}

				if(::fabs(zRightDifT - zLeftDifT) > ::fabs(zLeftDifT + zRightDifT)*HalfRelTol) // Discontinuity
				{
					int *p_izStartDiscontT = Aux_izStartDiscontT + ix;
					int Cur_iDzContinT = iz - *p_izStartDiscontT;
					if(iDzContinT > Cur_iDzContinT) iDzContinT = Cur_iDzContinT;
					*p_izStartDiscontT = iz;
				}
				if(::fabs(zRightDifP - zLeftDifP) > ::fabs(zLeftDifP + zRightDifP)*HalfRelTol) // Discontinuity
				{
					int *p_izStartDiscontP = Aux_izStartDiscontP + ix;
					int Cur_iDzContinP = iz - *p_izStartDiscontP;
					if(iDzContinP > Cur_iDzContinP) iDzContinP = Cur_iDzContinP;
					*p_izStartDiscontP = iz;
				}
			}
		}
	}

	DxContin = ((iDxContinT < iDxContinP)? iDxContinT : iDxContinP)*xStep;
	DzContin = ((iDzContinT < iDzContinP)? iDzContinT : iDzContinP)*zStep;

	if(Aux_izStartDiscontT != 0) delete[] Aux_izStartDiscontT;
	if(Aux_izStartDiscontP != 0) delete[] Aux_izStartDiscontP;
	return 0;
}

//*************************************************************************

void srTGenTransmission::EnsureTransmissionForField()
{
	const char TransUnits[] = "Transmission for Field";
	if(strcmp(TransUnits, GenTransNumData.DataUnits) == 0) return;

	long ne=1, nx=1, nz=1;
	if(GenTransNumData.AmOfDims == 2)
	{
		nx = (GenTransNumData.DimSizes)[0];
		nz = (GenTransNumData.DimSizes)[1];
	}
	else if(GenTransNumData.AmOfDims == 3)
	{
		ne = (GenTransNumData.DimSizes)[0];
		nx = (GenTransNumData.DimSizes)[1];
		nz = (GenTransNumData.DimSizes)[2];
	}

	DOUBLE *t = (DOUBLE*)(GenTransNumData.pData);
	//for(long iz=0; iz<(GenTransNumData.DimSizes)[1]; iz++)
	for(long iz=0; iz<nz; iz++)
	{
		//for(long ix=0; ix<(GenTransNumData.DimSizes)[0]; ix++)
		for(long ix=0; ix<nx; ix++)
		{
			for(long ie=0; ie<ne; ie++)
			{
				double Buf = sqrt(::fabs(*t));
				*t = Buf; t += 2;
			}
		}
	}
	strcpy(GenTransNumData.DataUnits, TransUnits);
}

//*************************************************************************

double srTGenTransmission::DetermineAppropriatePhotEnergyForFocDistTest(double Rx, double Rz)
{
	if((GenTransNumData.AmOfDims == 2) || ((GenTransNumData.AmOfDims == 3) && (GenTransNumData.DimStartValues[0] <= 0.)))
	{
		const double a = 1.239842e-06;
		const int Nm = 256;

		long NpDefX = (GenTransNumData.DimSizes)[0];
		double StartX = (GenTransNumData.DimStartValues)[0] + TransvCenPoint.x;
		double StepX = (GenTransNumData.DimSteps)[0];
		double EndX = StartX + (NpDefX - 1)*StepX;
		double AbsStartX = ::fabs(StartX), AbsEndX = ::fabs(EndX);
		double AbsMaxX = (AbsStartX > AbsEndX)? AbsStartX : AbsEndX;
		double eMidX = a*Rx*Nm/(8*AbsMaxX*AbsMaxX);

		long NpDefZ = (GenTransNumData.DimSizes)[1];
		double StartZ = (GenTransNumData.DimStartValues)[1] + TransvCenPoint.y;
		double StepZ = (GenTransNumData.DimSteps)[1];
		double EndZ = StartZ + (NpDefZ - 1)*StepZ;
		double AbsStartZ = ::fabs(StartZ), AbsEndZ = ::fabs(EndZ);
		double AbsMaxZ = (AbsStartZ > AbsEndZ)? AbsStartZ : AbsEndZ;
		double eMidZ = a*Rz*Nm/(8*AbsMaxZ*AbsMaxZ);

		double res_eMid = (eMidX > eMidZ)? eMidX : eMidZ;
		if(res_eMid < 1.e-04) res_eMid = 1.e-04; //OC071108
		else if(res_eMid > 5.e+04) res_eMid = 5.e+04;

		return res_eMid;
	}
	else //if(GenTransNumData.AmOfDims == 3)
	{
		long ne = (GenTransNumData.DimSizes)[0];
		long iec = ne >> 1;
		return (GenTransNumData.DimStartValues)[0] + iec*((GenTransNumData.DimSteps)[0]); //"central" photon energy
	}
}

//*************************************************************************

int srTGenTransmission::EstimateFocalDistancesAndCheckSampling()
{
	int result;
	const double AcceptedFocError = 0.25; // To steer

	bool eMidIsZero = false;
	if(eMid == 0.)
	{
		eMidIsZero = true;
		eMid = 1;
	}

	double SumFxDir = 0., SumFzDir = 0.;
	double AbsFxErrMax = 0., AbsFzErrMax = 0., FxPrev, FzPrev; 
	for(int k=0; k<3; k++)
	{
		double RelOtherCoord = 0.25*(k + 1);
		double Fx, Fz;
		if(result = DetermineFocalDistDirectAndCheckSampling('x', RelOtherCoord, Fx)) return result;
		if(result = DetermineFocalDistDirectAndCheckSampling('z', RelOtherCoord, Fz)) return result;

		if(((Fx > 1E+22) || (SumFxDir > 1E+22)) && ((Fz > 1E+22) || (SumFzDir > 1E+22))) //OC291009
		{
			SumFxDir = 3E+23; SumFzDir = 3E+23; break;
		}

		SumFxDir += Fx; SumFzDir += Fz;

		if(k > 0)
		{
			double CurrAbsFxErr = ::fabs(Fx - FxPrev);
			if(AbsFxErrMax < CurrAbsFxErr) AbsFxErrMax = CurrAbsFxErr;
			double CurrAbsFzErr = ::fabs(Fz - FzPrev);
			if(AbsFzErrMax < CurrAbsFzErr) AbsFzErrMax = CurrAbsFzErr;
		}
		FxPrev = Fx; FzPrev = Fz;
	}
	double FxDir = SumFxDir/3., FzDir = SumFzDir/3.;
	if((::fabs(FxDir) > 1.E+20) && (::fabs(FzDir) > 1.E+20))
	{
		FocDistX = FocDistZ = 1.E+23; return 0;
	}

	if(::fabs(FxDir)*AcceptedFocError < AbsFxErrMax) FxDir = 1.E+23;
	if(::fabs(FzDir)*AcceptedFocError < AbsFzErrMax) FzDir = 1.E+23;

	double NormalTestRadX[] = {0.5*FxDir, 2.*FxDir};
	double NormalTestRadZ[] = {0.5*FzDir, 2.*FzDir};
	double RescueTestRad[] = { 1., 5.}; // To steer

	double *TestRadX = NormalTestRadX, *TestRadZ = NormalTestRadZ;
	//if(::fabs(FxDir) > 1.E+20) TestRadX = RescueTestRad;
	//if(::fabs(FzDir) > 1.E+20) TestRadZ = RescueTestRad;
	if(::fabs(FxDir) > 1000.) TestRadX = RescueTestRad; //OC071108
	if(::fabs(FzDir) > 1000.) TestRadZ = RescueTestRad;

	if(eMidIsZero) eMid = DetermineAppropriatePhotEnergyForFocDistTest(TestRadX[1], TestRadZ[1]);

	double SumFx = 0., SumFz = 0.;
	AbsFxErrMax = 0.; AbsFzErrMax = 0.;
	int PassCountX = 0, PassCountZ = 0;
	bool PrevFxWasDefined = false, PrevFzWasDefined = false;
	for(int i=0; i<2; i++)
	{
		for(int k=0; k<3; k++)
		{
			double RelOtherCoord = 0.25*(k + 1);
			double Fx = 1.e+23, Fz = 1.e+23;
			srTRadSect1D PointSourceSectX, PointSourceSectZ;
			if(result = SetUpPointSourceSect1D('x', TestRadX[i], RelOtherCoord, PointSourceSectX)) return result;
			if(result = SetUpPointSourceSect1D('z', TestRadZ[i], RelOtherCoord, PointSourceSectZ)) return result;

			if((PointSourceSectX.pEx != 0) && (PointSourceSectX.np > 0)) //OC261009
			{
				if(result = DetermineFocalDistByPropag1D(PointSourceSectX, Fx)) 
				{
					//return result;
					Fx = 1.e+23;
				}
			}
			bool SkipThisX = (Fx == 0);

			if((PointSourceSectZ.pEx != 0) && (PointSourceSectZ.np > 0)) //OC261009
			{
				if(result = DetermineFocalDistByPropag1D(PointSourceSectZ, Fz)) 
				{
					//return result;
					Fz = 1.e+23;
				}
			}
			bool SkipThisZ = (Fz == 0);

			if(((Fx > 1E+22) || (SumFx > 1E+22)) && ((Fz > 1E+22) || (SumFz > 1E+22))) //OC291009
			{
				SumFx = 3E+23; SumFz = 3E+23; 
				PassCountX = PassCountZ = 3; 
				break;
			}

			SumFx += Fx; SumFz += Fz;

			//if((k > 0) || (i > 0))
			if(PrevFxWasDefined)
			{
				if(!SkipThisX)
				{
					double CurrAbsFxErr = ::fabs(Fx - FxPrev);
					if(AbsFxErrMax < CurrAbsFxErr) AbsFxErrMax = CurrAbsFxErr;
				}
			}
			if(PrevFzWasDefined)
			{
				if(!SkipThisZ)
				{
					double CurrAbsFzErr = ::fabs(Fz - FzPrev);
					if(AbsFzErrMax < CurrAbsFzErr) AbsFzErrMax = CurrAbsFzErr;
				}
			}

			if(!SkipThisX) { FxPrev = Fx; PrevFxWasDefined = true; PassCountX++;}
			else PrevFxWasDefined = false;

			if(!SkipThisZ) { FzPrev = Fz; PrevFzWasDefined = true; PassCountZ++;}
			else PrevFzWasDefined = false;
		}
	}
	//double FxProp = SumFx/6., FzProp = SumFz/6.;
	double FxProp = 1.E+23, FzProp = 1.E+23;
	if(PassCountX > 0) FxProp = SumFx/((double)PassCountX);
	if(PassCountZ > 0) FzProp = SumFz/((double)PassCountZ);

	if(::fabs(FxProp)*AcceptedFocError < AbsFxErrMax) FxProp = 1.E+23;
	if(::fabs(FzProp)*AcceptedFocError < AbsFzErrMax) FzProp = 1.E+23;

	if((::fabs(FxProp) > 1.E+20) && (::fabs(FzProp) > 1.E+20))
	{
		FocDistX = FocDistZ = 1.E+23; return 0;
	}

	if(::fabs(FxDir) < 1.E+20) FocDistX = (FxProp + FxDir)*0.5;
	else FocDistX = FxProp;

	if(::fabs(FzDir) < 1.E+20) FocDistZ = (FzProp + FzDir)*0.5;
	else FocDistZ = FzProp;

	if(::fabs(FocDistX) > 1.E+20) FocDistX = 1.E+23;
	if(::fabs(FocDistZ) > 1.E+20) FocDistZ = 1.E+23;

	if(eMidIsZero)
	{
		eMid = 0.;
	}
	return 0;
}

//*************************************************************************

void srTGenTransmission::EstimateEffPointsRange(char x_or_z, long icOtherCoord, long& iFirst, long& iLast, double& ArgFirst, double& ArgLast)
{
	const double AbsTolMagn = 0.1; // To steer

	long Period, InitialOffset, Np, Nx, Ne=1, iec=0, iDimX=0, iDimZ=1;
	//long Nx = (GenTransNumData.DimSizes)[0];
	
	if(GenTransNumData.AmOfDims == 2) 
	{
		Nx = (GenTransNumData.DimSizes)[0];
	}
	else if(GenTransNumData.AmOfDims == 3)
	{
		Ne = (GenTransNumData.DimSizes)[0];
		iec = Ne >> 1;

		iDimX = 1; iDimZ = 2;
		Nx = (GenTransNumData.DimSizes)[iDimX];
	}

	double ArgStep=0, ArgStart=0;

	if(x_or_z == 'x')
	{
		//Period = 2;
		//InitialOffset = icOtherCoord*(Nx << 1);
		
		Period = Ne << 1; //OC241112
		InitialOffset = icOtherCoord*(Nx*Period) + (iec << 1);

		Np = Nx;

		//ArgStep = (GenTransNumData.DimSteps)[0];
		//ArgStart = (GenTransNumData.DimStartValues)[0];
		ArgStep = (GenTransNumData.DimSteps)[iDimX];
		ArgStart = (GenTransNumData.DimStartValues)[iDimX]; // + TransvCenPoint.x;
	}
	else
	{
		//Period = Nx << 1;
		//InitialOffset = icOtherCoord << 1;
		//Np = (GenTransNumData.DimSizes)[1];
		//ArgStep = (GenTransNumData.DimSteps)[1];
		//ArgStart = (GenTransNumData.DimStartValues)[1] + TransvCenPoint.y;

		Period = Nx*(Ne << 1); //OC241112
		InitialOffset = icOtherCoord*(Ne << 1) + (iec << 1);

		Np = (GenTransNumData.DimSizes)[iDimZ];
		
		ArgStep = (GenTransNumData.DimSteps)[iDimZ];
		ArgStart = (GenTransNumData.DimStartValues)[iDimZ]; // + TransvCenPoint.y;
	}
	long Np_mi_1 = Np - 1;

	DOUBLE *tT0 = (DOUBLE*)(GenTransNumData.pData) + InitialOffset;
	DOUBLE *tT = tT0;

	double Tmax = 0., Tmin = 1.E+23;
	long i;
	for(i=0; i<Np; i++)
	{
		if(Tmax < *tT) Tmax = *tT;
		if(Tmin > *tT) Tmin = *tT;
		tT += Period;
	}

	double MagnThresh = Tmin + AbsTolMagn*(Tmax - Tmin);
	tT = tT0;
	DOUBLE *tT_Inv = tT0 + Period*(Np - 1);
	iFirst = -1; iLast = Np;
	for(i=0; i<Np; i++)
	{
		if(iFirst == -1)
		{
			if(*tT >= MagnThresh) iFirst = i;
		}
		if(iLast == Np)
		{
			if(*tT_Inv >= MagnThresh) iLast = Np_mi_1 - i;
		}
		tT += Period;
		tT_Inv -= Period;
	}
	ArgFirst = ArgStart + iFirst*ArgStep;
	ArgLast = ArgStart + iLast*ArgStep;
}

//*************************************************************************

int srTGenTransmission::DetermineFocalDistDirectAndCheckSampling(char x_or_z, double RelOtherCoord, double& FocDist)
{
	int result;
	srTRadSect1D Sect1D;
	double* PhaseCont=0;
	char PhaseIsContinuous;

	if(result = ExtractNumStructSect1DAndCheckSampling(x_or_z, RelOtherCoord, Sect1D, PhaseCont, PhaseIsContinuous)) return result;

	long iFirst, iLast;
	double ArgFirst, ArgLast;
	EstimateEffPointsRange(x_or_z, Sect1D.icOtherCoord, iFirst, iLast, ArgFirst, ArgLast);
	long Np = iLast - iFirst + 1;
	long Np_p_1 = Np + 1;

	if(!PhaseIsContinuous)
	{// This is not at all general: try to improve later...
		srTRadGenManip RadManip;
		RadManip.TryToMakePhaseContinuous1D(PhaseCont + iFirst, Np, -1, 0.);
	}

	float* PhaseContF = new float[Np_p_1];
	if(PhaseContF == 0) return MEMORY_ALLOCATION_FAILURE;
	float* ArgCont = new float[Np_p_1];
	if(ArgCont == 0) return MEMORY_ALLOCATION_FAILURE;
	float* Sigm = new float[Np_p_1];
	if(Sigm == 0) return MEMORY_ALLOCATION_FAILURE;

	const float RelErrorLevel = (float)0.1; // To steer
	//const double qAcceptLevel = 0.0001; // To steer

	double Arg = Sect1D.ArgStart, Step = Sect1D.ArgStep;
	double *tPhaseCont = PhaseCont + iFirst;
	float *tPhaseContF = PhaseContF + 1, *tArgCont = ArgCont + 1;
	double PhAbsMax = 0.;

	for(long j=1; j<=Np; j++)
	{
		double PhAbs = ::fabs(*tPhaseCont);
		if(PhAbsMax < PhAbs) PhAbsMax = PhAbs;

		*(tPhaseContF++) = (float)(*(tPhaseCont++));
		*(tArgCont++) = (float)Arg;
		Arg += Step;
	}

	if(PhAbsMax == 0.) PhAbsMax = 1.;
	float AbsErr = (float)(PhAbsMax*RelErrorLevel);

	float *tSigm = Sigm + 1;
	for(long k=1; k<=Np; k++) *(tSigm++) = AbsErr;

	float a[4], chisq, qOK;
	int ia[] = {1,1,1,1};
	CGenMathFit Fit;
	//if(result = Fit.FitPolynomial(ArgCont, PhaseContF, Sigm, int(Np), a, ia, 3, &chisq, &qOK)) return result;
	result = Fit.FitPolynomial(ArgCont, PhaseContF, Sigm, int(Np), a, ia, 3, &chisq, &qOK);
	if(result != 0) //OC291009
	{
		if(PhaseCont != 0) delete[] PhaseCont;
		if(PhaseContF != 0) delete[] PhaseContF;
		if(ArgCont != 0) delete[] ArgCont;
		if(Sigm != 0) delete[] Sigm;

		FocDist = 1E+23;
		return 0;
	}

	if(a[3] == 0.) a[3] = (float)(1.E-23); // To steer
	FocDist = -2.5338408E+06*eMid/a[3];
	//FocDist = -0.5/a[3];
	//OC

	if(::fabs(FocDist) > 1.E+20) FocDist = 1.E+23;

	if(PhaseCont != 0) delete[] PhaseCont;
	if(PhaseContF != 0) delete[] PhaseContF;
	if(ArgCont != 0) delete[] ArgCont;
	if(Sigm != 0) delete[] Sigm;
	return 0;
}

//*************************************************************************

int srTGenTransmission::ExtractNumStructSect1DAndCheckSampling(char x_or_z, double RelOtherCoord, srTRadSect1D& Sect1D, double*& PhaseCont, char& PhaseIsContinuous)
{
	double OtherStep, OtherStart;
	long OtherNp;

	long iDimX=0, iDimZ=1;
	if(GenTransNumData.AmOfDims == 3)
	{
		iDimX = 1; iDimZ = 2;
	}

	if(x_or_z == 'x')
	{
		//Sect1D.np = (GenTransNumData.DimSizes)[0];
		//Sect1D.ArgStep = (GenTransNumData.DimSteps)[0];
		//Sect1D.ArgStart = (GenTransNumData.DimStartValues)[0];

		//OtherStep = (GenTransNumData.DimSteps)[1];
		//OtherStart = (GenTransNumData.DimStartValues)[1];
		//OtherNp = (GenTransNumData.DimSizes)[1];
		//double OtherRange = ((GenTransNumData.DimSizes)[1] - 1)*OtherStep;
		//Sect1D.OtherCoordVal = (GenTransNumData.DimStartValues)[1] + RelOtherCoord*OtherRange;

		Sect1D.np = (GenTransNumData.DimSizes)[iDimX]; //OC241112
		Sect1D.ArgStep = (GenTransNumData.DimSteps)[iDimX];
		Sect1D.ArgStart = (GenTransNumData.DimStartValues)[iDimX];

		OtherStep = (GenTransNumData.DimSteps)[iDimZ];
		OtherStart = (GenTransNumData.DimStartValues)[iDimZ];
		OtherNp = (GenTransNumData.DimSizes)[iDimZ];
		double OtherRange = ((GenTransNumData.DimSizes)[iDimZ] - 1)*OtherStep;
		Sect1D.OtherCoordVal = (GenTransNumData.DimStartValues)[iDimZ] + RelOtherCoord*OtherRange;
	}
	else
	{
		//Sect1D.np = (GenTransNumData.DimSizes)[1];
		//Sect1D.ArgStep = (GenTransNumData.DimSteps)[1];
		//Sect1D.ArgStart = (GenTransNumData.DimStartValues)[1];

		//OtherStep = (GenTransNumData.DimSteps)[0];
		//OtherStart = (GenTransNumData.DimStartValues)[0];
		//OtherNp = (GenTransNumData.DimSizes)[0];
		//double OtherRange = ((GenTransNumData.DimSizes)[0] - 1)*OtherStep;
		//Sect1D.OtherCoordVal = (GenTransNumData.DimStartValues)[0] + RelOtherCoord*OtherRange;

		Sect1D.np = (GenTransNumData.DimSizes)[iDimZ]; //OC241112
		Sect1D.ArgStep = (GenTransNumData.DimSteps)[iDimZ];
		Sect1D.ArgStart = (GenTransNumData.DimStartValues)[iDimZ];

		OtherStep = (GenTransNumData.DimSteps)[iDimX];
		OtherStart = (GenTransNumData.DimStartValues)[iDimX];
		OtherNp = (GenTransNumData.DimSizes)[iDimX];
		double OtherRange = ((GenTransNumData.DimSizes)[iDimX] - 1)*OtherStep;
		Sect1D.OtherCoordVal = (GenTransNumData.DimStartValues)[iDimX] + RelOtherCoord*OtherRange;
	}
	Sect1D.icOtherCoord = (long)((Sect1D.OtherCoordVal - OtherStart)/OtherStep);
	long OtherNp_mi_2 = OtherNp - 2;
	if(Sect1D.icOtherCoord > OtherNp_mi_2) Sect1D.icOtherCoord = OtherNp_mi_2;

	Sect1D.pEx = 0;
	Sect1D.pEz = 0;
	Sect1D.DeleteArraysAtDestruction = 0;

	PhaseCont = new double[Sect1D.np];
	if(PhaseCont == 0) return MEMORY_ALLOCATION_FAILURE;

	Sect1D.VsXorZ = x_or_z;
	Sect1D.Robs = 1.E+23;
	Sect1D.RobsAbsErr = 1.E+23;
	Sect1D.eVal = eMid;
	Sect1D.WfrEdgeCorrShouldBeDone = 0;
	Sect1D.Pres = 0;
	Sect1D.LengthUnit = 0;
	Sect1D.PhotEnergyUnit = 0;
	char NameWave[] = "AuxOptCompSetup";
	strcpy(Sect1D.NameOfWave, NameWave);

	CopyNumStructValuesToSect1DAndCheckSampling(Sect1D, PhaseCont, PhaseIsContinuous);
	return 0;
}

//*************************************************************************

void srTGenTransmission::CopyNumStructValuesToSect1DAndCheckSampling(srTRadSect1D& Sect1D, double* PhaseCont, char& PhaseIsContinuousOneExtrem)
{
	//long Period, InitialOffset, Np, Nz, Nx, Ne=1, iec=0;
	long long Period, InitialOffset, Np, Nz, Nx, Ne=1, iec=0;

	//Nx = (GenTransNumData.DimSizes)[0];

	if(GenTransNumData.AmOfDims == 2) 
	{
		Nx = (GenTransNumData.DimSizes)[0];
		Nz = (GenTransNumData.DimSizes)[1];
	}
	else if(GenTransNumData.AmOfDims == 3)
	{
		Ne = (GenTransNumData.DimSizes)[0];
		iec = Ne >> 1;
		Nx = (GenTransNumData.DimSizes)[1];
		Nz = (GenTransNumData.DimSizes)[2];
	}

	if(Sect1D.VsXorZ == 'x')
	{
		//Period = 2;
		//InitialOffset = Sect1D.icOtherCoord*(Nx << 1);

		Period = Ne << 1; //OC241112
		InitialOffset = Sect1D.icOtherCoord*Nx*Period + (iec << 1);

		Np = Nx;
	}
	else
	{
		//Period = Nx << 1;
		//InitialOffset = Sect1D.icOtherCoord << 1;

		Period = Nx*(Ne << 1);
		InitialOffset = Sect1D.icOtherCoord*(Ne << 1) + (iec << 1);

		Np = Nz;
	}
	//long Np_mi_1 = Np - 1;
	long long Np_mi_1 = Np - 1;

	DOUBLE *tPh0 = (DOUBLE*)(GenTransNumData.pData) + InitialOffset + 1;
	DOUBLE *tPh = tPh0;

	double *tPhOut = PhaseCont;
	double PhMax = -1.E+23, PhMin = 1.E+23;

	int AmOfExtrem = 0, AmOfPtBwExtrem = 0, MinAmOfPtBwExtrem = 100000;
	char DerSign = (*tPh0 <= *(tPh0 + Period))? 1 : -1;
	const double RelTol = 0.00001; // To steer

	double OptPathMult = eMid*5.0676816042E+06;
	//double InvOptPathMult = 1./OptPathMult;

	double AddPh = -(*tPh0);
	if(OptPathOrPhase == 1) AddPh *= OptPathMult; // TwoPi_d_Lambda_m
	//if(OptPathOrPhase == 2) AddPh *= InvOptPathMult; 
	//OC

	//for(long i=0; i<Np; i++)
	for(long long i=0; i<Np; i++)
	{
		double Ph = *tPh;
		if(OptPathOrPhase == 1) Ph *= OptPathMult; // TwoPi_d_Lambda_m
		//if(OptPathOrPhase == 2) Ph *= InvOptPathMult;
		//OC
		Ph += AddPh;

		if(PhMax < Ph) PhMax = Ph;
		if(PhMin > Ph) PhMin = Ph;

		if(i != Np_mi_1)
		{
			char NewDerSign = DerSign;
			double AbsTol = RelTol*(::fabs(Ph));
			double Dif = *(tPh + Period) - *tPh;
			if(OptPathOrPhase == 1) Dif *= OptPathMult;
			//if(OptPathOrPhase == 2) Dif *= InvOptPathMult;
			//OC

			if(Dif > AbsTol) NewDerSign = 1;
			else if(Dif < -AbsTol) NewDerSign = -1;

			if(NewDerSign != DerSign)
			{
				AmOfExtrem++;
				if(MinAmOfPtBwExtrem > AmOfPtBwExtrem) MinAmOfPtBwExtrem = AmOfPtBwExtrem;
				AmOfPtBwExtrem = 0;
				DerSign = NewDerSign;
			}
			else
			{
				AmOfPtBwExtrem++;
			}
		}

		*(tPhOut++) = Ph;
		tPh += Period;
	}

	PhaseIsContinuousOneExtrem = (AmOfExtrem < 2);
	//if((MinAmOfPtBwExtrem == 0) && (!PhaseIsContinuousOneExtrem))
	//{
		//srTSend Send; Send.AddWarningMessage(&gVectWarnNos, POOR_PHASE_SHIFT_SAMPLING);
		//CErrWarn::AddWarningMessage(&gVectWarnNos, POOR_PHASE_SHIFT_SAMPLING);
		//OC150208 commented out 
	//}

	double Ph1 = *tPh0, Ph2 = *(tPh0 + (Np << 1));
	if((Ph1 <= PhMax) && (Ph2 <= PhMax)) AddPh = -PhMax;
	else AddPh = -PhMin;

	tPhOut = PhaseCont;
	//for(long k=0; k<Np; k++) *(tPhOut++) += AddPh;
	for(long long k=0; k<Np; k++) *(tPhOut++) += AddPh;
}

//*************************************************************************

int srTGenTransmission::SetUpPointSourceSect1D(char x_or_z, double R, double RelOtherCoord, srTRadSect1D& PointSourceSect1D)
{
	const long SmallestN = 10;
	const long LargestN = 1000000;
	const long OverSampFact = 8;
	CGenMathFFT1D FFT;

	double WavelengthIn_m = 1.239854*1.E-06/eMid;
	double HalfLambR = 0.5*WavelengthIn_m*R;

	long iDimX=0, iDimZ=1;
	if(GenTransNumData.AmOfDims == 3)
	{
		iDimX = 1; iDimZ = 2;
	}

	long NpDef, OtherNpDef;
	double Start, End, Step, OtherStart, OtherEnd, OtherStep;
	if(x_or_z == 'x')
	{
		//NpDef = (GenTransNumData.DimSizes)[0];
		//Start = (GenTransNumData.DimStartValues)[0] + TransvCenPoint.x;
		//Step = (GenTransNumData.DimSteps)[0];

		//OtherNpDef = (GenTransNumData.DimSizes)[1];
		//OtherStart = (GenTransNumData.DimStartValues)[1] + TransvCenPoint.y;
		//OtherStep = (GenTransNumData.DimSteps)[1];

		NpDef = (GenTransNumData.DimSizes)[iDimX]; //OC241112
		Start = (GenTransNumData.DimStartValues)[iDimX]; // + TransvCenPoint.x;
		Step = (GenTransNumData.DimSteps)[iDimX];

		OtherNpDef = (GenTransNumData.DimSizes)[iDimZ];
		OtherStart = (GenTransNumData.DimStartValues)[iDimZ]; // + TransvCenPoint.y;
		OtherStep = (GenTransNumData.DimSteps)[iDimZ];
	}
	else
	{
		//NpDef = (GenTransNumData.DimSizes)[1];
		//Start = (GenTransNumData.DimStartValues)[1] + TransvCenPoint.y;
		//Step = (GenTransNumData.DimSteps)[1];

		//OtherNpDef = (GenTransNumData.DimSizes)[0];
		//OtherStart = (GenTransNumData.DimStartValues)[0] + TransvCenPoint.x;
		//OtherStep = (GenTransNumData.DimSteps)[0];

		NpDef = (GenTransNumData.DimSizes)[iDimZ]; //OC241112
		Start = (GenTransNumData.DimStartValues)[iDimZ]; // + TransvCenPoint.y;
		Step = (GenTransNumData.DimSteps)[iDimZ];

		OtherNpDef = (GenTransNumData.DimSizes)[iDimX];
		OtherStart = (GenTransNumData.DimStartValues)[iDimX]; // + TransvCenPoint.x;
		OtherStep = (GenTransNumData.DimSteps)[iDimX];
	}
	End = Start + (NpDef - 1)*Step;
	OtherEnd = OtherStart + (OtherNpDef - 1)*OtherStep;

	double dStart = (Start != 0.)? ::fabs(HalfLambR/Start) : ::fabs(HalfLambR/Step);
	double dEnd = (End != 0.)? ::fabs(HalfLambR/End) : ::fabs(HalfLambR/Step);
	double d = (dStart < dEnd)? dStart : dEnd;
	long Np = (long(::fabs(End - Start)/d) + 1)*OverSampFact;

	if(Np >= LargestN) //OC261009
	{
		PointSourceSect1D.pEx = 0;
		PointSourceSect1D.np = 0;
		return 0;
	}

	if(((Np>>1)<<1) != Np) Np++;
	FFT.NextCorrectNumberForFFT(Np);
	if(Np < SmallestN) Np = SmallestN;

	PointSourceSect1D.pEx = new float[Np << 1];
	if(PointSourceSect1D.pEx == 0) return MEMORY_ALLOCATION_FAILURE;
	PointSourceSect1D.pEz = 0;
	PointSourceSect1D.DeleteArraysAtDestruction = 1;

	PointSourceSect1D.ArgStep = (End - Start)/(Np - 1);
	PointSourceSect1D.ArgStart = Start;
	PointSourceSect1D.np = Np;
	PointSourceSect1D.VsXorZ = x_or_z;
	PointSourceSect1D.Robs = R;
	PointSourceSect1D.RobsAbsErr = 0.01*R; // To steer

	PointSourceSect1D.OtherCoordVal = OtherStart + RelOtherCoord*(OtherEnd - OtherStart);
	PointSourceSect1D.icOtherCoord = (long)((PointSourceSect1D.OtherCoordVal - OtherStart)/OtherStep);

	PointSourceSect1D.eVal = eMid;
	PointSourceSect1D.WfrEdgeCorrShouldBeDone = 0;
	PointSourceSect1D.Pres = 0;
	PointSourceSect1D.LengthUnit = 0;
	PointSourceSect1D.PhotEnergyUnit = 0;

	char NameWave[] = "AuxOptCompSetup";
	strcpy(PointSourceSect1D.NameOfWave, NameWave);

	double Pi_d_Lambda_m = eMid*2.533840802E+06;
	double Pi_d_Lambda_d_R = Pi_d_Lambda_m/R;

	double xStep = PointSourceSect1D.ArgStep;
	double x = Start;
	float *tEx = PointSourceSect1D.pEx;
	for(long i=0; i<Np; i++)
	{
		double Ph = Pi_d_Lambda_d_R*x*x;
		CosAndSin(Ph, *tEx, *(tEx+1));
		tEx += 2;
		x += xStep;
	}

	return 0;
}

//*************************************************************************

int srTGenTransmission::DetermineFocalDistByPropag1D(srTRadSect1D& Sect1D, double& FocDist)
{
	int result;
	if(result = PropagateRadiationSimple1D(&Sect1D)) return result;

	double* PhaseCont = new double[Sect1D.np];
	if(PhaseCont == 0) return MEMORY_ALLOCATION_FAILURE;

	float *tE = Sect1D.pEx;
	double *tPh = PhaseCont;
	bool E_isZero = true;
	for(long i=0; i<Sect1D.np; i++)
	{
		if((*tE != 0) || (*(tE + 1) != 0)) E_isZero = false;
		*(tPh++) = FormalPhase(*tE, *(tE + 1)); tE += 2;
	}
	if(E_isZero)
	{
		if(PhaseCont != 0) delete[] PhaseCont;
		FocDist = 1E+23;
		return 0;
	}

	long iFirst, iLast;
	double ArgFirst, ArgLast;
	EstimateEffPointsRange(Sect1D.VsXorZ, Sect1D.icOtherCoord, iFirst, iLast, ArgFirst, ArgLast);
	iFirst = (long)((ArgFirst - Sect1D.ArgStart)/Sect1D.ArgStep);
	if(iFirst >= Sect1D.np) iFirst = Sect1D.np - 1;
	iLast = (long)((ArgLast - Sect1D.ArgStart)/Sect1D.ArgStep);
	if(iLast >= Sect1D.np) iLast = Sect1D.np - 1;
	long Np = iLast - iFirst + 1;

	const long maxNp = 1000000; 
	if(Np >= maxNp) //OC261009
	{
		if(PhaseCont != 0) delete[] PhaseCont;
		FocDist = 1E+23;
		return 0; //??
	}

	long Np_p_1 = Np + 1;

	srTRadGenManip RadManip;
	RadManip.TryToMakePhaseContinuous1D(PhaseCont + iFirst, Np, -1, 0.);

	float* PhaseContF = new float[Np_p_1];
	if(PhaseContF == 0) return MEMORY_ALLOCATION_FAILURE;
	float* ArgCont = new float[Np_p_1];
	if(ArgCont == 0) return MEMORY_ALLOCATION_FAILURE;
	float* Sigm = new float[Np_p_1];
	if(Sigm == 0) return MEMORY_ALLOCATION_FAILURE;

	const float RelErrorLevel = (float)0.1; // To steer
	//const double qAcceptLevel = 0.01; // To steer
	const double RelTolR = 0.05; // To steer

	double ArgRange = Sect1D.ArgStep*(Np - 1);
	double InvArgRange = 1./ArgRange;
	double Arg = (Sect1D.ArgStart + iFirst*Sect1D.ArgStep)*InvArgRange, Step = Sect1D.ArgStep*InvArgRange;

	double *tPhaseCont = PhaseCont + iFirst;
	float *tPhaseContF = PhaseContF + 1, *tArgCont = ArgCont + 1;
	double PhAbsMax = 0.;

	for(long j=1; j<=Np; j++)
	{
		double PhAbs = ::fabs(*tPhaseCont);
		if(PhAbsMax < PhAbs) PhAbsMax = PhAbs;

		*(tPhaseContF++) = (float)(*(tPhaseCont++));
		*(tArgCont++) = (float)Arg;
		Arg += Step;
	}

	float AbsErr = (float)(PhAbsMax*RelErrorLevel);
	if(AbsErr == 0.) AbsErr = RelErrorLevel;

	float *tSigm = Sigm + 1;
	for(long k=1; k<=Np; k++) *(tSigm++) = AbsErr;

	float a[4], chisq, qOK;
	int ia[] = {1,1,1,1};
	CGenMathFit Fit;
	//if(result = Fit.FitPolynomial(ArgCont, PhaseContF, Sigm, int(Np), a, ia, 3, &chisq, &qOK)) return result;
	result = Fit.FitPolynomial(ArgCont, PhaseContF, Sigm, int(Np), a, ia, 3, &chisq, &qOK);
	if(result != 0) //OC291009
	{
		if(PhaseCont != 0) delete[] PhaseCont;
		if(PhaseContF != 0) delete[] PhaseContF;
		if(ArgCont != 0) delete[] ArgCont;
		if(Sigm != 0) delete[] Sigm;

		FocDist = 1E+23;
		return 0;
	}

	if(a[3] == 0.) a[3] = (float)(1.E-23); // To steer

	double Rafter = (2.5338408E+06)*eMid*ArgRange*ArgRange/a[3];
	double Rbefore = Sect1D.Robs;

	if((::fabs(Rbefore - Rafter) > RelTolR*(::fabs(Rbefore)))) 
		FocDist = Rbefore*Rafter/(Rafter - Rbefore);
	else FocDist = 1.E+23;

	if(PhaseCont != 0) delete[] PhaseCont;
	if(PhaseContF != 0) delete[] PhaseContF;
	if(ArgCont != 0) delete[] ArgCont;
	if(Sigm != 0) delete[] Sigm;
	return 0;
}

//*************************************************************************

void srTGenTransmission::RadPointModifier(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
{// e in eV; Length in m !!!
 // Operates on Coord. side !!!
	//double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
	double xRel = EXZ.x, zRel = EXZ.z; //OC080311

	long Ne = 1, Nemi2 = -1;
	long iDimX = 0, iDimZ = 1;
	if(GenTransNumData.AmOfDims == 3)
	{
		Ne = (GenTransNumData.DimSizes)[0];
		Nemi2 = Ne - 2;
		iDimX = 1; iDimZ = 2;
	}

	//long Nx = (GenTransNumData.DimSizes)[0], Nz = (GenTransNumData.DimSizes)[1];
	long Nx = (GenTransNumData.DimSizes)[iDimX], Nz = (GenTransNumData.DimSizes)[iDimZ]; //OC241112
	long Nxmi2 = Nx - 2, Nzmi2 = Nz - 2;
	
	//double xStart = (GenTransNumData.DimStartValues)[0], zStart = (GenTransNumData.DimStartValues)[1];
	//double xStep = (GenTransNumData.DimSteps)[0], zStep = (GenTransNumData.DimSteps)[1];
	double xStart = (GenTransNumData.DimStartValues)[iDimX], zStart = (GenTransNumData.DimStartValues)[iDimZ];
	double xStep = (GenTransNumData.DimSteps)[iDimX], zStep = (GenTransNumData.DimSteps)[iDimZ];

	double xEnd = xStart + (Nx - 1)*xStep, zEnd = zStart + (Nz - 1)*zStep;

	double AbsTolX = xStep*0.001, AbsTolZ = zStep*0.001; // To steer
	if(OuterTransmIs == 1)
	{
		if((xRel < xStart - AbsTolX) || (xRel > xEnd + AbsTolX) || (zRel < zStart - AbsTolZ) || (zRel > zEnd + AbsTolZ))
		{
			if(EPtrs.pExRe != 0) { *(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;}
			if(EPtrs.pEzRe != 0) { *(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;}
			return;
		}
	}

	double xr = 0., zr = 0.;
	double T = 1., Ph = 0.;
	//char NotExactRightEdgeX = 1, NotExactRightEdgeZ = 1;

	long ix = long((xRel - xStart)/xStep);
	if(::fabs(xRel - ((ix + 1)*xStep + xStart)) < 1.E-05*xStep) ix++;

	//if(ix < 0) { ix = 0; xr = 0.;}
	//else if(ix > Nxmi2) { ix = Nx - 1; xr = 0.; NotExactRightEdgeX = 0;}
	//else xr = (xRel - (ix*xStep + xStart))/xStep;

	if(ix < 0) ix = 0; //OC241112
	//else if(ix > Nxmi2) ix = Nxmi2;
	//xr = (xRel - (ix*xStep + xStart))/xStep;
	else if(ix > Nxmi2) { ix = Nxmi2; xr = 1.;}
	else xr = (xRel - (ix*xStep + xStart))/xStep;

	long iz = long((zRel - zStart)/zStep);
	if(::fabs(zRel - ((iz + 1)*zStep + zStart)) < 1.E-05*zStep) iz++;

	//if(iz < 0) { iz = 0; zr = 0.;}
	//else if(iz > Nzmi2) { iz = Nz - 1; zr = 0.; NotExactRightEdgeZ = 0;}
	//else zr = (zRel - (iz*zStep + zStart))/zStep;

	if(iz < 0) iz = 0;
	//else if(iz > Nzmi2) iz = Nzmi2;
	//zr = (zRel - (iz*zStep + zStart))/zStep;
	else if(iz > Nzmi2) { iz = Nzmi2; zr = 1.;}
	else zr = (zRel - (iz*zStep + zStart))/zStep;
	
	double xrzr = xr*zr;
	if((GenTransNumData.AmOfDims == 2) || ((GenTransNumData.AmOfDims == 3) && (Ne == 1)))
	{
		//long zPer = Nx << 1;
		long long zPer = Nx << 1;

		DOUBLE *p00 = (DOUBLE*)(GenTransNumData.pData) + (iz*zPer + (ix << 1));
		DOUBLE *p10 = p00 + 2, *p01 = p00 + zPer;
		DOUBLE *p11 = p01 + 2;

		DOUBLE *p00p1 = p00+1, *p10p1 = p10+1, *p01p1 = p01+1, *p11p1 = p11+1;

		//double Axz = 0., Ax = 0., Az = 0., Bxz = 0., Bx = 0., Bz = 0.;
		//if(NotExactRightEdgeX && NotExactRightEdgeZ) { Axz = *p00 - *p01 - *p10 + *p11; Bxz = *p00p1 - *p01p1 - *p10p1 + *p11p1;}
		//if(NotExactRightEdgeX) { Ax = (*p10 - *p00); Bx = (*p10p1 - *p00p1);}
		//if(NotExactRightEdgeZ) { Az = (*p01 - *p00); Bz = (*p01p1 - *p00p1);}

		double Axz = *p00 - *p01 - *p10 + *p11, Bxz = *p00p1 - *p01p1 - *p10p1 + *p11p1;
		double Ax = (*p10 - *p00), Bx = (*p10p1 - *p00p1);
		double Az = (*p01 - *p00), Bz = (*p01p1 - *p00p1);

		T = Axz*xrzr + Ax*xr + Az*zr + *p00;
		Ph = Bxz*xrzr + Bx*xr + Bz*zr + *p00p1;
	}
	else if(GenTransNumData.AmOfDims == 3)
	{//bi-linear 3D interpolation
		double eStart = (GenTransNumData.DimStartValues)[0];
		double eStep = (GenTransNumData.DimSteps)[0];

		long ie = long((EXZ.e - eStart)/eStep + 1.e-10);
		if(ie < 0) ie = 0;
		else if(ie > Nemi2) ie = Nemi2;

		double er = (EXZ.e - (ie*eStep + eStart))/eStep;
		//double erxr = er*xr, erzr = er*zr;
		//double erxrzr = erxr*zr;

		//long xPer = Ne << 1;
		//long zPer = Nx*xPer;
		long long xPer = Ne << 1;
		long long zPer = Nx*xPer;
		DOUBLE *p000 = (DOUBLE*)(GenTransNumData.pData) + (iz*zPer + ix*xPer + (ie << 1));
		DOUBLE *p100 = p000 + 2, *p010 = p000 + xPer, *p001 = p000 + zPer;
		DOUBLE *p110 = p100 + xPer, *p101 = p100 + zPer, *p011 = p010 + zPer;
		DOUBLE *p111 = p110 + zPer;

		double one_mi_er = 1.- er, one_mi_xr = 1.- xr, one_mi_zr = 1.- zr;
		double one_mi_er_one_mi_xr = one_mi_er*one_mi_xr, er_one_mi_xr = er*one_mi_xr;
		double one_mi_er_xr = one_mi_er*xr, er_xr = er*xr;
		T = ((*p000)*one_mi_er_one_mi_xr + (*p100)*er_one_mi_xr + (*p010)*one_mi_er_xr + (*p110)*er_xr)*one_mi_zr
		  + ((*p001)*one_mi_er_one_mi_xr + (*p101)*er_one_mi_xr + (*p011)*one_mi_er_xr + (*p111)*er_xr)*zr;
		Ph = ((*(p000+1))*one_mi_er_one_mi_xr + (*(p100+1))*er_one_mi_xr + (*(p010+1))*one_mi_er_xr + (*(p110+1))*er_xr)*one_mi_zr
		   + ((*(p001+1))*one_mi_er_one_mi_xr + (*(p101+1))*er_one_mi_xr + (*(p011+1))*one_mi_er_xr + (*(p111+1))*er_xr)*zr;

	 // inArFunc[] = {f(x0,y0,z0),f(x1,y0,z0),f(x0,y1,z0),f(x0,y0,z1),f(x1,y1,z0),f(x1,y0,z1),f(x0,y1,z1),f(x1,y1,z1)} //function values at the corners of the cube
		//return inArFunc[0]*one_mi_xt*one_mi_yt*one_mi_zt
		//	+ inArFunc[1]*xt*one_mi_yt*one_mi_zt
		//	+ inArFunc[2]*one_mi_xt*yt*one_mi_zt
		//	+ inArFunc[3]*one_mi_xt*one_mi_yt*zt
		//	+ inArFunc[4]*xt*yt*one_mi_zt
		//	+ inArFunc[5]*xt*one_mi_yt*zt
		//	+ inArFunc[6]*one_mi_xt*yt*zt
		//	+ inArFunc[7]*xt*yt*zt;
	}

	if(OptPathOrPhase == 1) Ph *= EXZ.e*5.0676816042E+06; // TwoPi_d_Lambda_m
	float CosPh, SinPh; CosAndSin(Ph, CosPh, SinPh);
	if(EPtrs.pExRe != 0)
	{
		float NewExRe = (float)(T*((*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh));
		float NewExIm = (float)(T*((*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh));
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
	}
	if(EPtrs.pEzRe != 0)
	{
		float NewEzRe = (float)(T*((*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh));
		float NewEzIm = (float)(T*((*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh));
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm; 
	}
}

//*************************************************************************

void srTGenTransmission::RadPointModifier1D(srTEXZ& EXZ, srTEFieldPtrs& EPtrs)
{// e in eV; Length in m !!!
 // Operates on Coord. side !!!
	//double xRel = EXZ.x - TransvCenPoint.x, zRel = EXZ.z - TransvCenPoint.y;
	double xRel = EXZ.x, zRel = EXZ.z; //OC080311

	long Ne = 1, Nemi2 = -1;
	long iDimX = 0, iDimZ = 1;
	if(GenTransNumData.AmOfDims == 3)
	{
		Ne = (GenTransNumData.DimSizes)[0];
		Nemi2 = Ne - 2;
		iDimX = 1; iDimZ = 2;
	}

	//long Nx = (GenTransNumData.DimSizes)[0], Nz = (GenTransNumData.DimSizes)[1];
	long Nx = (GenTransNumData.DimSizes)[iDimX], Nz = (GenTransNumData.DimSizes)[iDimZ]; //OC241112
	long Nxmi2 = Nx - 2, Nzmi2 = Nz - 2;

	//double xStart = (GenTransNumData.DimStartValues)[0], zStart = (GenTransNumData.DimStartValues)[1];
	//double xStep = (GenTransNumData.DimSteps)[0], zStep = (GenTransNumData.DimSteps)[1];
	double xStart = (GenTransNumData.DimStartValues)[iDimX], zStart = (GenTransNumData.DimStartValues)[iDimZ];
	double xStep = (GenTransNumData.DimSteps)[iDimX], zStep = (GenTransNumData.DimSteps)[iDimZ];

	double xEnd = xStart + (Nx - 1)*xStep, zEnd = zStart + (Nz - 1)*zStep;

	double AbsTolX = xStep*0.001, AbsTolZ = zStep*0.001; // To steer
	if(OuterTransmIs == 1)
	{
		if((xRel < xStart - AbsTolX) || (xRel > xEnd + AbsTolX) || (zRel < zStart - AbsTolZ) || (zRel > zEnd + AbsTolZ))
		{
			if(EPtrs.pExRe != 0) { *(EPtrs.pExRe) = 0.; *(EPtrs.pExIm) = 0.;}
			if(EPtrs.pEzRe != 0) { *(EPtrs.pEzRe) = 0.; *(EPtrs.pEzIm) = 0.;}
			return;
		}
	}
	double xr, zr;
	long ix = (long)((xRel - xStart)/xStep + 1.e-08);
	if(ix < 0) { ix = 0; xr = 0.;}
	else if(ix > Nxmi2) { ix = Nxmi2; xr = 1.;}
	else xr = (xRel - (ix*xStep + xStart))/xStep;

	long iz = (long)((zRel - zStart)/zStep + 1.e-08);
	if(iz < 0) { iz = 0; zr = 0.;}
	else if(iz > Nzmi2) { iz = Nzmi2; zr = 1.;}
	else zr = (zRel - (iz*zStep + zStart))/zStep;

	double T=1., Ph=0.;

	if((GenTransNumData.AmOfDims == 2) || ((GenTransNumData.AmOfDims == 3) && (Ne == 1)))
	{
		//long zPer = Nx << 1;
		long long zPer = Nx << 1;
		DOUBLE *p0 = (DOUBLE*)(GenTransNumData.pData) + (iz*zPer + (ix << 1));

		if(EXZ.VsXorZ == 'x')
		{
			T = (*(p0 + 2) - *p0)*xr + *p0; p0++;
			Ph = (*(p0 + 2) - *p0)*xr + *p0;
		}
		else
		{
			T = (*(p0 + zPer) - *p0)*zr + *p0; p0++;
			Ph = (*(p0 + zPer) - *p0)*zr + *p0;
		}
	}
	else if(GenTransNumData.AmOfDims == 3)
	{
		double eStart = (GenTransNumData.DimStartValues)[0];
		double eStep = (GenTransNumData.DimSteps)[0];

		long ie = long((EXZ.e - eStart)/eStep + 1.e-10);
		if(ie < 0) ie = 0;
		else if(ie > Nemi2) ie = Nemi2;

		double er = (EXZ.e - (ie*eStep + eStart))/eStep;

		//long xPer = Ne << 1;
		//long zPer = Nx*xPer;
		long long xPer = Ne << 1;
		long long zPer = Nx*xPer;
		DOUBLE *p00 = (DOUBLE*)(GenTransNumData.pData) + (iz*zPer + ix*xPer + (ie << 1));
		DOUBLE *p10 = p00 + 2;
		DOUBLE *p01=0, *p11=0;
		double one_mi_er = 1.- er, ar = 1., one_mi_ar = 1.;

		if(EXZ.VsXorZ == 'x')
		{
			p01 = p00 + xPer; p11 = p10 + xPer; ar = xr; one_mi_ar = 1.- xr;
		}
		else
		{
			p01 = p00 + zPer; p11 = p10 + zPer; ar = zr; one_mi_ar = 1.- zr;
		}

		double one_mi_er_one_mi_ar = one_mi_er*one_mi_ar, er_one_mi_ar = er*one_mi_ar;
		double one_mi_er_ar = one_mi_er*ar, er_ar = er*ar;

		T = (*p00)*one_mi_er_one_mi_ar + (*p10)*er_one_mi_ar + (*p01)*one_mi_er_ar + (*p11)*er_ar;
		Ph = (*(p00+1))*one_mi_er_one_mi_ar + (*(p10+1))*er_one_mi_ar + (*(p01+1))*one_mi_er_ar + (*(p11+1))*er_ar;
	}

	if(OptPathOrPhase == 1) Ph *= EXZ.e*5.0676816042E+06; // TwoPi_d_Lambda_m
	float CosPh, SinPh; CosAndSin(Ph, CosPh, SinPh);
	if(EPtrs.pExRe != 0)
	{
		float NewExRe = (float)(T*((*(EPtrs.pExRe))*CosPh - (*(EPtrs.pExIm))*SinPh));
		float NewExIm = (float)(T*((*(EPtrs.pExRe))*SinPh + (*(EPtrs.pExIm))*CosPh));
		*(EPtrs.pExRe) = NewExRe; *(EPtrs.pExIm) = NewExIm; 
	}
	if(EPtrs.pEzRe != 0)
	{
		float NewEzRe = (float)(T*((*(EPtrs.pEzRe))*CosPh - (*(EPtrs.pEzIm))*SinPh));
		float NewEzIm = (float)(T*((*(EPtrs.pEzRe))*SinPh + (*(EPtrs.pEzIm))*CosPh));
		*(EPtrs.pEzRe) = NewEzRe; *(EPtrs.pEzIm) = NewEzIm;
	}
}

//*************************************************************************

int srTGenTransmission::EstimateMinNpToResolveOptElem(srTSRWRadStructAccessData* pRadAccessData, double& MinNx, double& MinNz)
{
	const double MinPo = 40; // To steer
	const double MinPoPerContInt = 1.; // To steer
	MinNx = MinNz = MinPo;

	double xRange = pRadAccessData->xStep*(pRadAccessData->nx - 1);
	double MinNxToResolveDiscont = (xRange/DxContin)*MinPoPerContInt;
	if(MinNx <  MinNxToResolveDiscont) MinNx =  MinNxToResolveDiscont;
		
	double zRange = pRadAccessData->zStep*(pRadAccessData->nz - 1);
	double MinNzToResolveDiscont = (zRange/DzContin)*MinPoPerContInt;
	if(MinNz <  MinNzToResolveDiscont) MinNz =  MinNzToResolveDiscont;

	return 0;
}

//*************************************************************************
