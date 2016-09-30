/************************************************************************//**
 * File: srsend.cpp
 * Description: Interface (input-output) functions essentially for IGOR Pro XOP
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srsend.h"

#include <stdio.h>

#ifdef __MATHEMATICA__
extern "C" {
#include <mathlink.h>
}
#endif

#include "srtrjdat.h"
#include "srradint.h"
#include "srebmdat.h"
#include "srpersto.h"
#include "srstowig.h"
#include "srradinc.h"
#include "srpowden.h"
//#include "sroptics.h" //OC221012
//#include "srwarn.h" //OC021112
#include "srctrjdt.h"
#include "srisosrc.h"
#include "srgsnbm.h"
#include "srmagfld.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void srTSend::ErrorMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__
	char err_msg[300];
	sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");
	MLClearError(stdlink);
	MLNewPacket(stdlink);
	MLEvaluate(stdlink, err_msg);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLPutSymbol(stdlink, "$Failed");
#endif
	return;
}

//-------------------------------------------------------------------------

void srTSend::WarningMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__
	char err_msg[300];
	sprintf(err_msg, "%s%s%s", "Message[",MessageString,"]");
	MLNewPacket(stdlink);
	MLEvaluate(stdlink, err_msg);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
#endif
	return;
}

//*************************************************************************

void srTSend::WarningMessages(srTIntVect* pWarnMesNos)
{
#if defined(__IGOR_PRO__)

	if(pWarnMesNos->empty()) return;

	ostringstream DecoratedStream;
	DecoratedStream << "\rS R W   W A R N I N G\r";

	CErrWarn srwlErWar;

	for(srTIntVect::iterator iter = pWarnMesNos->begin(); iter != pWarnMesNos->end(); ++iter)
	{
		//const char* WarningText = srAllWarnings[*iter - SRW_WARNINGS_OFFSET - 1];
		const char* WarningText = srwlErWar.GetWarning(*iter); //OC031112
		if(WarningText != 0)
		{
			DecoratedStream << "- " << WarningText << "\r";
		}
		else DecoratedStream << "Sorry, the warning text is not available.\r";
	}
	DecoratedStream << "\r" << ends;
	basic_string<char> BufDecoratedWarnString = DecoratedStream.str();
	const char* DecoratedWarnString = BufDecoratedWarnString.c_str();

	XOPNotice(DecoratedWarnString);
	pWarnMesNos->erase(pWarnMesNos->begin(), pWarnMesNos->end());

#endif
}

//*************************************************************************

void srTSend::AddWarningMessage(srTIntVect* pWarnMesNos, int WarnNo)
{
	for(srTIntVect::iterator iter = pWarnMesNos->begin(); iter != pWarnMesNos->end(); ++iter)
	{
		if(*iter == WarnNo) return;
	}
	pWarnMesNos->push_back(WarnNo);
}

//*************************************************************************

void srTSend::OrdinaryMessage(const char* MessageString)
{
#ifdef __MATHEMATICA__
	char InfoMessage[500];
	sprintf(InfoMessage, "%s%s%s", "Message[",MessageString,"]");
	MLEvaluate(stdlink, InfoMessage);
	MLNextPacket(stdlink);
	MLNewPacket(stdlink);
	MLPutSymbol(stdlink, "Null");
#endif
}

//-------------------------------------------------------------------------

void srTSend::String(char* MessageString)
{
#ifdef __MATHEMATICA__
	MLPutString(stdlink, MessageString);
#endif
}

//-------------------------------------------------------------------------

void srTSend::Double(double d)
{
#ifdef __MATHEMATICA__
	if(d==0.) d=1.E-17; 
	MLPutDouble(stdlink, d);
#endif
}

//-------------------------------------------------------------------------

void srTSend::MyMLPutDouble(double d)
{
#ifdef __MATHEMATICA__
	if(d==0.) d=1.E-17; 
	MLPutDouble(stdlink, d);
#endif
}

//-------------------------------------------------------------------------

void srTSend::DoubleList(double* ArrayOfDouble, int lenArrayOfDouble)
{
#ifdef __MATHEMATICA__
	InitOutList(lenArrayOfDouble);
	for(int i=0; i<lenArrayOfDouble; i++) MyMLPutDouble(ArrayOfDouble[i]);
#endif
}

//*************************************************************************

void srTSend::DoubleListWithArg(double ArgSt, double ArgEn, double* ArrayOfDouble, int lenArrayOfDouble)
{
#ifdef __MATHEMATICA__
	InitOutList(lenArrayOfDouble);
	double Arg = ArgSt;
	double StepArg = (ArgEn - ArgSt)/(lenArrayOfDouble - 1);
	for(int i=0; i<lenArrayOfDouble; i++)
	{
		InitOutList(2);
		MyMLPutDouble(Arg);
		MyMLPutDouble(ArrayOfDouble[i]);
		Arg += StepArg;
	}
#endif
}

//*************************************************************************

void srTSend::Long(long LongIntValue)
{
#ifdef __MATHEMATICA__
	MLPutLongInteger(stdlink, LongIntValue);
#endif
}

//-------------------------------------------------------------------------

void srTSend::Int(int IntValue)
{
#ifdef __MATHEMATICA__
	MLPutInteger(stdlink, IntValue);
#endif
}

//-------------------------------------------------------------------------

void srTSend::IntList(int* ArrayOfInt, int lenArrayOfInt)
{
#ifdef __MATHEMATICA__
	InitOutList(lenArrayOfInt);
	for(int i=0; i<lenArrayOfInt; i++) Int(ArrayOfInt[i]);
#endif
}

//-------------------------------------------------------------------------

void srTSend::InitOutList(int NumberOfElem)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", NumberOfElem);
#endif
}

//-------------------------------------------------------------------------

void srTSend::SubArbNestedArrays(double* Data, int* Dims, int Depth, int& CntData)
{
#ifdef __MATHEMATICA__
	for(int i=0; i<Dims[Depth-1]; i++)
	{
		if(Depth>1)
		{
			MLPutFunction(stdlink, "List", Dims[Depth-2]);
			SubArbNestedArrays(Data, Dims, Depth-1, CntData);
		}
		else
		{
			double Buf = Data[CntData];
			MyMLPutDouble(Data[CntData++]);
		}
	}
#endif
}

//-------------------------------------------------------------------------

void srTSend::ArbNestedArrays(double* Data, int* Dims, int Depth)
{
#ifdef __MATHEMATICA__
	MLPutFunction(stdlink, "List", Dims[Depth-1]);
	int CntData =0;
	SubArbNestedArrays(Data, Dims, Depth, CntData);
#endif
}

//*************************************************************************

int srTSend::InitTrjOutFormat1(srTTrjDat& TrjDat)
{
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	if(pIgorTrjInputStruct != 0) 
	{
		wavH_OutBtxData = pIgorTrjInputStruct->wOutHorAng;
		wavH_OutXData = pIgorTrjInputStruct->wOutHorCoor;
		wavH_OutBtzData = pIgorTrjInputStruct->wOutVerAng;
		wavH_OutZData = pIgorTrjInputStruct->wOutVerCoor;
	}

	if(wavH_OutBtxData == NIL) return NOWAV;
	if(WaveType(wavH_OutBtxData) != NT_FP64) return NT_FP64_WAVE_REQUIRED;
	if(wavH_OutXData == NIL) return NOWAV;
	if(WaveType(wavH_OutXData) != NT_FP64) return NT_FP64_WAVE_REQUIRED;
	if(wavH_OutBtzData == NIL) return NOWAV;
	if(WaveType(wavH_OutBtzData) != NT_FP64) return NT_FP64_WAVE_REQUIRED;
	if(wavH_OutZData == NIL) return NOWAV;
	if(WaveType(wavH_OutZData) != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH_OutBtxData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutBtxData = 0; //MoveLockHandle(wavH_OutBtxData); //OC180815
	pOutBtxData = (DOUBLE*)((char*)(*wavH_OutBtxData) + dataOffset);
	if(result = MDAccessNumericWaveData(wavH_OutXData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutXData = 0; //MoveLockHandle(wavH_OutXData); //OC180815
	pOutXData = (DOUBLE*)((char*)(*wavH_OutXData) + dataOffset);
	if(result = MDAccessNumericWaveData(wavH_OutBtzData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutBtzData = 0; //MoveLockHandle(wavH_OutBtzData); //OC180815
	pOutBtzData = (DOUBLE*)((char*)(*wavH_OutBtzData) + dataOffset);
	if(result = MDAccessNumericWaveData(wavH_OutZData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutZData = 0; //MoveLockHandle(wavH_OutZData); //OC180815
	pOutZData = (DOUBLE*)((char*)(*wavH_OutZData) + dataOffset);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::InitRadDistrOutFormat1(srTWfrSmp& DistrInfoDat)
{
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	int result = 0;

	if(wavH_OutExData == NIL) return NOWAV;
	int waveType = WaveType(wavH_OutExData);
	if(waveType != (NT_FP32 | NT_CMPLX)) return NT_FP64_OR_NT_CMPLX_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];

	if(result = MDGetWaveDimensions(wavH_OutExData, &numDimensions, dimensionSizes)) return result;
	if (numDimensions != 3) return NEEDS_3D_WAVE;

	long numRowsInWave = dimensionSizes[0];
	long numColumnsInWave = dimensionSizes[1];
	long numLayersInWave = dimensionSizes[2];

	long numRows = long(DistrInfoDat.nLamb);
	long numColumns = long(DistrInfoDat.nx);
	long numLayers = long(DistrInfoDat.nz);

	if((numRowsInWave != numRows) || (numColumnsInWave != numColumns) || (numLayersInWave != numLayers))
		return BAD_RAD_OUTPUT_WAVE_DIMENSIONS;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH_OutExData, kMDWaveAccessMode0, &dataOffset)) return result;

	hStateOutExData = 0; //MoveLockHandle(wavH_OutExData); //OC180815
	pOutExData = (float*)((char*)(*wavH_OutExData) + dataOffset);

	if(wavH_OutEzData == NIL) return NOWAV;
	waveType = WaveType(wavH_OutEzData);
	if(waveType != (NT_FP32 | NT_CMPLX)) return NT_FP64_OR_NT_CMPLX_WAVE_REQUIRED;

	if(result = MDGetWaveDimensions(wavH_OutEzData, &numDimensions, dimensionSizes)) return result;
	if (numDimensions != 3) return NEEDS_3D_WAVE;

	numRowsInWave = dimensionSizes[0];
	numColumnsInWave = dimensionSizes[1];
	numLayersInWave = dimensionSizes[2];

	if((numRowsInWave != numRows) || (numColumnsInWave != numColumns) || (numLayersInWave != numLayers))
		return BAD_RAD_OUTPUT_WAVE_DIMENSIONS;

	if(result = MDAccessNumericWaveData(wavH_OutEzData, kMDWaveAccessMode0, &dataOffset)) return result;

	hStateOutEzData = 0; //MoveLockHandle(wavH_OutEzData); //OC180815
	pOutEzData = (float*)((char*)(*wavH_OutEzData) + dataOffset);
	tOutExData = pOutExData; tOutEzData = pOutEzData;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::InitRadDistrOutFormat2(srTWfrSmp& DistrInfoDat)
{// This is used for diffraction/propagation tests
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	char waveName[MAX_OBJ_NAME+1];
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	strcpy(waveName, "FocusedFluxDensX");
	dimensionSizes[0] = (DistrInfoDat.nx > 1)? DistrInfoDat.nx : 0;
	dimensionSizes[1] = (DistrInfoDat.nz > 1)? DistrInfoDat.nz : 0;
	dimensionSizes[2] = (DistrInfoDat.ny > 1)? DistrInfoDat.ny : 0;
	dimensionSizes[3] = 0;

	//if(result = MDMakeWave(&wavH_OutExData, waveName, NIL, dimensionSizes, NT_FP32, 1)) return result;
	if(result = MDMakeWave(&wavH_OutExData, waveName, NIL, dimensionSizes, NT_FP32, 1)) return FAILED_TO_CREATE_WAVE;

	if(DistrInfoDat.nx > 1)
	{
		DOUBLE xStep = ((DistrInfoDat.xEnd - DistrInfoDat.xStart)/(DistrInfoDat.nx - 1));
		DOUBLE xStart = DistrInfoDat.xStart;
		if(result = MDSetWaveScaling(wavH_OutExData, ROWS, &xStep, &xStart)) return result;

		//if(result = MDSetWaveUnits(wavH_OutExData, ROWS, "m")) return result;
		char strUnits[] = "m\0";
		if(result = MDSetWaveUnits(wavH_OutExData, ROWS, strUnits)) return result;
	}

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH_OutExData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutExData = 0; //MoveLockHandle(wavH_OutExData); //OC180815
	pOutExData = (float*)((char*)(*wavH_OutExData) + dataOffset);

	strcpy(waveName, "FocusedFluxDensZ");
	//if(result = MDMakeWave(&wavH_OutEzData, waveName, NIL, dimensionSizes, NT_FP32, 1)) return result;
	if(result = MDMakeWave(&wavH_OutEzData, waveName, NIL, dimensionSizes, NT_FP32, 1)) return FAILED_TO_CREATE_WAVE;
	if(result = MDAccessNumericWaveData(wavH_OutEzData, kMDWaveAccessMode0, &dataOffset)) return result;
	hStateOutEzData =  0; //MoveLockHandle(wavH_OutEzData); //OC180815
	pOutEzData = (float*)((char*)(*wavH_OutEzData) + dataOffset);

	tOutExData = pOutExData; tOutEzData = pOutEzData;

	return 0;

#else 

	return 0;

#endif
}

//*************************************************************************

int srTSend::InitRadDistrOutFormat3(srTSRWRadStructAccessData& SRWRadStructAccessData, srTWfrSmp& DistrInfoDat)
{
#if defined(__IGOR_PRO__)

	wavH_OutExData = SRWRadStructAccessData.wRadX;
	wavH_OutEzData = SRWRadStructAccessData.wRadZ;

	int result;
	if(result = InitRadDistrOutFormat1(DistrInfoDat)) return result;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::OutRadDistrFormat1(srTRadInt& RadInt)
{
#if defined(__MATHEMATICA__)

	srTWfrSmp& DistrInfoDat = RadInt.DistrInfoDat;
	InitOutList(2);

// Sending header
	InitOutList(5);
	char* Str = 0;

	InitOutList(3);

	if(DistrInfoDat.CoordOrAngPresentation == AngPres) Str = "\"mr\"";
	else if(DistrInfoDat.CoordOrAngPresentation == CoordPres) Str = "\"mm\"";

	if(RadInt.DistrInfoDat.ShowPhaseOnly) Str = (RadInt.DistrInfoDat.PhaseDerOrder = 0)? "\"Ph\"" : ((RadInt.PhaseDerOrder = 1)? "\"dPhds\"" : "\"d2Phds2\"");

	String(Str);

	if(DistrInfoDat.BandwidthUnits == 0) Str = "\".1%\"";
	else if(DistrInfoDat.BandwidthUnits == 1) Str = "\"keV\"";
	else if(DistrInfoDat.BandwidthUnits == 2) Str = "\"eV\"";
	else if(DistrInfoDat.BandwidthUnits == 3) Str = "\"Ang\"";
	else if(DistrInfoDat.BandwidthUnits == 4) Str = "\"nm\"";
	else if(DistrInfoDat.BandwidthUnits == 5) Str = "\"mic\"";
	String(Str);

	if(DistrInfoDat.DistrValType == FieldFourier) Str = "\"Field\"";
	else if(DistrInfoDat.DistrValType == StokesParam) Str = "\"Stokes\"";
	String(Str);

	double LocLambStart, LocLambEnd;

	int *pNArgLoop[4];
	for(int i=0; i<4; i++)
	{
		InitOutList(2);

		double *pArgStart, *pArgEnd;
		int *pNArg;

		if(DistrInfoDat.LoopOrder[i] == 'w')
		{
			DistrInfoDat.RetrievePhotonEnergyOrWavelength(LocLambStart, LocLambEnd);

			pArgStart = &LocLambStart;
			pArgEnd = &LocLambEnd;
			pNArg = &(DistrInfoDat.nLamb);
			if(DistrInfoDat.PhotonEnergyWavelengthUnits == 0) Str = "\"keV\"";
			else if(DistrInfoDat.PhotonEnergyWavelengthUnits == 1) Str = "\"eV\"";
			else if(DistrInfoDat.PhotonEnergyWavelengthUnits == 2) Str = "\"Ang\"";
			else if(DistrInfoDat.PhotonEnergyWavelengthUnits == 3) Str = "\"nm\"";
			else if(DistrInfoDat.PhotonEnergyWavelengthUnits == 4) Str = "\"mic\"";
		}
		else if(DistrInfoDat.LoopOrder[i] == 'x')
		{
			pArgStart = &(DistrInfoDat.xStart);
			pArgEnd = &(DistrInfoDat.xEnd);
			pNArg = &(DistrInfoDat.nx);
			Str = "\"x\"";
		}
		else if(DistrInfoDat.LoopOrder[i] == 'z')
		{
			pArgStart = &(DistrInfoDat.zStart);
			pArgEnd = &(DistrInfoDat.zEnd);
			pNArg = &(DistrInfoDat.nz);
			Str = "\"z\"";
		}
		else if(DistrInfoDat.LoopOrder[i] == 'y')
		{
			pArgStart = &(DistrInfoDat.yStart);
			pArgEnd = &(DistrInfoDat.yEnd);
			pNArg = &(DistrInfoDat.ny);
			Str = "\"long\"";
		}
		String(Str);

		if(*pNArg > 1) 
		{ 
			InitOutList(3); Double(*pArgStart); Double(*pArgEnd); Int(*pNArg);
		}
		else Double(*pArgStart);
		pNArgLoop[i] = pNArg;
	}

// Sending comp. results
	if(RadInt.DistrInfoDat.ShowPhaseOnly) { OutPhaseFormat1(RadInt); return 0;}

	complex<double>* LocHorPolTravers = RadInt.RadDistrFieldFourierHorPol;
	complex<double>* LocVerPolTravers = RadInt.RadDistrFieldFourierVerPol;

	int NArgLoop0 = *(pNArgLoop[0]);
	int NArgLoop1 = *(pNArgLoop[1]);
	int NArgLoop2 = *(pNArgLoop[2]);
	int NArgLoop3 = *(pNArgLoop[3]);

	if(NArgLoop0 > 1) InitOutList(NArgLoop0);
	for(int iLoop0 = 0; iLoop0 < NArgLoop0; iLoop0++)
	{
		if(NArgLoop1 > 1) InitOutList(NArgLoop1);
		for(int iLoop1 = 0; iLoop1 < NArgLoop1; iLoop1++)
		{
			if(NArgLoop2 > 1) InitOutList(NArgLoop2);
			for(int iLoop2 = 0; iLoop2 < NArgLoop2; iLoop2++)
			{
				if(NArgLoop3 > 1) InitOutList(NArgLoop3);
				for(int iLoop3 = 0; iLoop3 < NArgLoop3; iLoop3++)
				{
					if(DistrInfoDat.DistrValType == FieldFourier)
					{
						InitOutList(2);
						InitOutList(2); 
						Double((*LocHorPolTravers).real());
						Double((*(LocHorPolTravers++)).imag());
						InitOutList(2); 
						Double((*LocVerPolTravers).real());
						Double((*(LocVerPolTravers++)).imag());
					}
					else if(DistrInfoDat.DistrValType == StokesParam)
					{
						InitOutList(4);
						Double((*LocHorPolTravers).real());
						Double((*(LocHorPolTravers++)).imag());
						Double((*LocVerPolTravers).real());
						Double((*(LocVerPolTravers++)).imag());
					}
				}
			}
		}
	}
	return 0;

#elif defined(__IGOR_PRO__)

	if(RadInt.DistrInfoDat.ShowPhaseOnly) { OutPhaseFormat1(RadInt); return 0;}

	waveHndl wavH = wavH_OutExData;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if((waveType != NT_FP64) && (waveType != (NT_FP32 | NT_CMPLX))) return NT_FP64_OR_NT_CMPLX_WAVE_REQUIRED;

	if(waveType == NT_FP64)
	{
		complex<double> LocEwX = *(RadInt.RadDistrFieldFourierHorPol);
		complex<double> LocEwZ = *(RadInt.RadDistrFieldFourierVerPol);

		RadInt.ReformatRadDistrCompRes(); // Compute Stokes Par.
		complex<double>* LocHorPolTravers = RadInt.RadDistrFieldFourierHorPol;
		complex<double>* LocVerPolTravers = RadInt.RadDistrFieldFourierVerPol;

		int result;

		long numDimensions;
		long dimensionSizes[MAX_DIMENSIONS+1];
		if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
		long numRows = dimensionSizes[0];
		if(numRows < 4) return BAD_STOKES_WAVE_FORMAT;

		long dataOffset;
		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		int hState = 0; //MoveLockHandle(wavH); //OC180815
		char* dataStartPtr = (char*)(*wavH) + dataOffset;
		DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

		*(dp0++) = (*LocHorPolTravers).real();
		*(dp0++) = (*(LocHorPolTravers++)).imag();
		*(dp0++) = (*LocVerPolTravers).real();
		*dp0 = (*(LocVerPolTravers++)).imag();

		HSetState((Handle)wavH, hState);
		WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.

		wavH = wavH_OutEzData;
		if(wavH == NIL) return NOWAV;
		waveType = WaveType(wavH);
		if(waveType != (NT_FP64 | NT_CMPLX)) return NT_CMPLX_WAVE_REQUIRED;

		if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
		numRows = dimensionSizes[0];
		if(numRows < 2) return BAD_FIELD_FOURIER_WAVE_FORMAT;

		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		hState = 0; //MoveLockHandle(wavH); //OC180815
		dataStartPtr = (char*)(*wavH) + dataOffset;
		dp0 = (DOUBLE*)dataStartPtr;

		*(dp0++) = LocEwX.real();
		*(dp0++) = LocEwX.imag();
		*(dp0++) = LocEwZ.real();
		*dp0 = LocEwZ.imag();

		HSetState((Handle)wavH, hState);
		WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
	}
	else
	{
		int result = 0;
		long numDimensions;
		long dimensionSizes[MAX_DIMENSIONS+1];
		if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
		if(numDimensions != 3) return NEEDS_3D_WAVE;

		long numRowsInWave = dimensionSizes[0];
		long numColumnsInWave = dimensionSizes[1];
		long numLayersInWave = dimensionSizes[2];

		long numRows = long(RadInt.DistrInfoDat.nLamb);
		long numColumns = long(RadInt.DistrInfoDat.nx);
		long numLayers = long(RadInt.DistrInfoDat.nz);

		if((numRowsInWave != numRows) || (numColumnsInWave != numColumns) || (numLayersInWave != numLayers))
			return BAD_RAD_OUTPUT_WAVE_DIMENSIONS;

		//long pointsPerColumn = numRows*2;
		//long pointsPerLayer = pointsPerColumn*numColumns;
		long long pointsPerColumn = numRows*2;
		long long pointsPerLayer = pointsPerColumn*numColumns;

		long dataOffset;
		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		int hState = 0; //MoveLockHandle(wavH); //OC180815
		char* dataStartPtr = (char*)(*wavH) + dataOffset;

		float *fp0 = (float*)dataStartPtr;
		complex<double>* LocRadTravers = RadInt.RadDistrFieldFourierHorPol;

		long layer, column, row;
		for(layer=0; layer<numLayers; layer++) 
		{
			float *flp = fp0 + layer*pointsPerLayer;
			for(column=0; column<numColumns; column++) 
			{
				float *fp = flp + column*pointsPerColumn;
				for(row=0; row<numRows; row++) 
				{
					*(fp++) = float(LocRadTravers->real());
					*(fp++) = float((LocRadTravers++)->imag());
				}
			}
		}
		HSetState((Handle)wavH, hState);
		WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.

		wavH = wavH_OutEzData;
		if(wavH == NIL) return NOWAV;
		waveType = WaveType(wavH);
		if((waveType != NT_FP64) && (waveType != (NT_FP32 | NT_CMPLX))) return NT_FP64_OR_NT_CMPLX_WAVE_REQUIRED;

		if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
		if (numDimensions != 3) return NEEDS_3D_WAVE;

		numRowsInWave = dimensionSizes[0];
		numColumnsInWave = dimensionSizes[1];
		numLayersInWave = dimensionSizes[2];

		if((numRowsInWave != numRows) || (numColumnsInWave != numColumns) || (numLayersInWave != numLayers))
		return BAD_RAD_OUTPUT_WAVE_DIMENSIONS;

		if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
		hState = 0; //MoveLockHandle(wavH); //OC180815
		dataStartPtr = (char*)(*wavH) + dataOffset;

		fp0 = (float*)dataStartPtr;
		LocRadTravers = RadInt.RadDistrFieldFourierVerPol;

		for(layer=0; layer<numLayers; layer++) 
		{
			float *flp = fp0 + layer*pointsPerLayer;
			for(column=0; column<numColumns; column++) 
			{
				float *fp = flp + column*pointsPerColumn;
				for(row=0; row<numRows; row++) 
				{
					*(fp++) = float(LocRadTravers->real());
					*(fp++) = float((LocRadTravers++)->imag());
				}
			}
		}
		HSetState((Handle)wavH, hState);
		WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.
	}

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::OutRadDistrFormat2(srTSRWRadStructAccessData& SRWRadStructAccessData, srTRadInt& RadInt)
{
#if defined(__IGOR_PRO__)

	wavH_OutExData = SRWRadStructAccessData.wRadX;
	wavH_OutEzData = SRWRadStructAccessData.wRadZ;

	int result;
	if(result = OutRadDistrFormat1(RadInt)) return result;

	return 0;

#else 

	return 0;

#endif
}

//*************************************************************************

void srTSend::OutPhaseFormat1(srTRadInt& RadInt)
{
#ifdef __MATHEMATICA__

	int AmOfP = int((RadInt.sIntegFin - RadInt.sIntegStart)/RadInt.sIntegStep) + 1;

	InitOutList(AmOfP);

	double sArg = RadInt.sIntegStart;
	double* tAuxPhaseArray = RadInt.AuxPhaseArray;
	for(int i=0; i<AmOfP; i++)
	{
		InitOutList(2);
		Double(sArg); Double(*(tAuxPhaseArray++)); sArg += RadInt.sIntegStep;
	}
#endif
}

//*************************************************************************

int srTSend::GetCoordData(double& xStart, double& xEnd, int& nx)
{
#ifdef __MATHEMATICA__

	const char *FunName;
	long ArgNum;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { ErrorMessage("SR::Error000"); return -1;}
		MLDisownSymbol(stdlink, FunName);

		if(!MLGetDouble(stdlink, &xStart)) { ErrorMessage("SR::Error000"); return -1;}
		if(!MLGetDouble(stdlink, &xEnd)) { ErrorMessage("SR::Error000"); return -1;}
		if(!MLGetInteger(stdlink, &nx)) { ErrorMessage("SR::Error000"); return -1;}
	}
	else if((Next==MLTKINT) || (Next==MLTKREAL))
	{
		if(!MLGetDouble(stdlink, &xStart)) { ErrorMessage("SR::Error000"); return -1;}
		xEnd = xStart; nx = 1;
	}
	else { ErrorMessage("SR::Error000"); return -1;}

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetPrtclTrjInitData(double& s0, double& x0, double& dxds0, double& z0, double& dzds0)
{
#if defined(__MATHEMATICA__)

	double* DataPtr;
	long LenData;

	int ArrayReadOK = MLGetRealList(stdlink, &DataPtr, &LenData);

	if(!(ArrayReadOK && (LenData==1 || LenData==5)))
	{ 
		ErrorMessage("SR::Error007"); return -1;
	}

	s0 = DataPtr[0];
	if(LenData==1) x0 = dxds0 = z0 = dzds0 = 0.;
	else
	{
		x0 = DataPtr[1]; dxds0 = DataPtr[2];
		z0 = DataPtr[3]; dzds0 = DataPtr[4]; 
	}

	MLDisownRealList(stdlink, DataPtr, LenData);
	return 0;

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetString(char*& String, int MaxSize)
{
#ifdef __MATHEMATICA__

	const char* BufString; 
	if(!MLGetString(stdlink, &BufString)) { ErrorMessage("SR::Error000"); return -1;};
	strncpy(String, BufString, MaxSize);
	MLDisownString(stdlink, BufString);

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetDouble(double& aDouble)
{
#ifdef __MATHEMATICA__

	if(!MLGetDouble(stdlink, &aDouble)) { ErrorMessage("SR::Error000"); return -1;};

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTotalFieldDataFormat1(srTMagFieldAccessData& MagFieldAccessData, srTTrjDat& TrjDat)
{
#if defined(__MATHEMATICA__)

	const char *FunName;
	long ArgNum;
	double DoubleValue;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("SR::Error000"); return -1;}
	if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=4) { ErrorMessage("SR::Error000"); return -1;}
	MLDisownSymbol(stdlink, FunName);

	Next = MLGetNext(stdlink);
	if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) TrjDat.sStart = DoubleValue;}
	else { ErrorMessage("SR::Error000"); return -1;}
	Next = MLGetNext(stdlink);
	double sFin;
	if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) sFin = DoubleValue;}
	else { ErrorMessage("SR::Error000"); return -1;}

	srTFunDer*& BxInData = TrjDat.BxInData;
	srTFunDer*& BzInData = TrjDat.BzInData;
	double& FieldZeroTolerance = TrjDat.FieldZeroTolerance;
	char BxIsNotZero=0, BzIsNotZero=0;

	char BxFilledOK = 0;
	long BxLen = 0;
	double Bx;
	Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &BxLen);
		if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("SR::Error000"); return -1;}
		MLDisownSymbol(stdlink, FunName);

		if(BxLen > 1)
		{
			if(BxInData != 0) { delete[] BxInData; BxInData = 0;}
			BxInData = new srTFunDer[BxLen];
			if(BxInData == 0) { ErrorMessage("SR::Error900"); return -1;}

			for(int k=0; k<BxLen; k++)
			{
				double* pBxVal = &(BxInData[k].f);
				if(!MLGetDouble(stdlink, pBxVal)) { ErrorMessage("SR::Error000"); return -1;}
				if(!BxIsNotZero) if(::fabs(*pBxVal) > FieldZeroTolerance) BxIsNotZero = 1;
			}
			BxFilledOK = 1;
		}
		else
		{
			if(!MLGetDouble(stdlink, &Bx)) { ErrorMessage("SR::Error000"); return -1;}
			if(::fabs(Bx) > FieldZeroTolerance) BxIsNotZero = 1;
		}
	}
	else if((Next==MLTKINT) || (Next==MLTKREAL))
	{
		if(!MLGetDouble(stdlink, &Bx)) { ErrorMessage("SR::Error000"); return -1;}
		if(::fabs(Bx) > FieldZeroTolerance) BxIsNotZero = 1;
	}
	else { ErrorMessage("SR::Error000"); return -1;}

	char BzFilledOK = 0;
	long BzLen = 0;
	double Bz;
	Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC)
	{
		ReadOK = MLGetFunction(stdlink, &FunName, &BzLen);
		if((!ReadOK) || strcmp(FunName, "List")) { ErrorMessage("SR::Error000"); return -1;}
		MLDisownSymbol(stdlink, FunName);

		if(BzLen > 1)
		{
			if(BzInData != 0) { delete[] BzInData; BzInData = 0;}
			BzInData = new srTFunDer[BzLen];
			if(BzInData == 0) { ErrorMessage("SR::Error900"); return -1;}

			for(int k=0; k<BzLen; k++)
			{
				double* pBzVal = &(BzInData[k].f);
				if(!MLGetDouble(stdlink, pBzVal)) { ErrorMessage("SR::Error000"); return -1;}
				if(!BzIsNotZero) if(::fabs(*pBzVal) > FieldZeroTolerance) BzIsNotZero = 1;
			}
			BzFilledOK = 1;
		}
		else
		{
			if(!MLGetDouble(stdlink, &Bz)) { ErrorMessage("SR::Error000"); return -1;}
			if(::fabs(Bz) > FieldZeroTolerance) BzIsNotZero = 1;
		}
	}
	else if((Next==MLTKINT) || (Next==MLTKREAL))
	{
		if(!MLGetDouble(stdlink, &Bz)) { ErrorMessage("SR::Error000"); return -1;}
		if(::fabs(Bz) > FieldZeroTolerance) BzIsNotZero = 1;
	}
	else { ErrorMessage("SR::Error000"); return -1;}

	if(BxFilledOK && BzFilledOK && (BxLen != BzLen))
		TrjDat.RecomputeHorOrVertField((int)BxLen, (int)BzLen);

	if(BxFilledOK || BzFilledOK)
		TrjDat.LenFieldData = (BxLen > BzLen)? BxLen : BzLen;
	else TrjDat.LenFieldData = 10.;

	if(!BxFilledOK)
	{
		if(TrjDat.BxInData != 0) { delete[] TrjDat.BxInData; TrjDat.BxInData = 0;}
		TrjDat.BxInData = new srTFunDer[TrjDat.LenFieldData];
		for(int k=0; k<TrjDat.LenFieldData; k++) (TrjDat.BxInData)[k].f = Bx;
	}
	if(!BzFilledOK)
	{
		if(TrjDat.BzInData != 0) { delete[] TrjDat.BzInData; TrjDat.BzInData = 0;}
		TrjDat.BzInData = new srTFunDer[TrjDat.LenFieldData];
		for(int k=0; k<TrjDat.LenFieldData; k++) (TrjDat.BzInData)[k].f = Bz;
	}

	TrjDat.HorFieldIsNotZero = BxIsNotZero;
	TrjDat.VerFieldIsNotZero = BzIsNotZero;

	TrjDat.sStep = (TrjDat.LenFieldData > 1)? (sFin - TrjDat.sStart)/(TrjDat.LenFieldData - 1) : 0;
	return 0;

#elif defined(__IGOR_PRO__)

	srTFunDer*& BxInData = TrjDat.BxInData;
	srTFunDer*& BzInData = TrjDat.BzInData;
	double& FieldZeroTolerance = TrjDat.FieldZeroTolerance;

	char LocInputWasNotModified = 1;
	waveHndl wavH = MagFieldAccessData.wBx;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long BxLen = dimensionSizes[0];
	if((numDimensions == 0) || (BxLen == 0)) return BAD_MAGN_FIELD_WAVE_FORMAT;

	DOUBLE sBxStep, sBxStart;
	if(result = MDGetWaveScaling(wavH, 0, &sBxStep, &sBxStart)) return result;

	char BxLenWasModified = (TrjDat.AuxBxLen != BxLen)? 1 : 0;
	if(BxLenWasModified)
	{
		if(BxInData != 0) delete[] BxInData;
		BxInData = new srTFunDer[BxLen];
		if(BxInData == 0) return MEMORY_ALLOCATION_FAILURE;

		LocInputWasNotModified = 0;
		TrjDat.AuxBxLen = BxLen;
	}

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) { delete[] BxInData; return result;}
	int hState = 0; //MoveLockHandle(wavH); //OC180815
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	srTFunDer* tBxInData = BxInData;
	double BxVal;
	char BxIsZero = 1;
	for(int k=0; k<BxLen; k++)
	{
		BxVal = *(dp0++);

		if(::fabs(BxVal) > FieldZeroTolerance) BxIsZero = 0;
		else BxVal = 0.;

		if(LocInputWasNotModified) if(tBxInData->f != BxVal) LocInputWasNotModified = 0;
		(tBxInData++)->f = BxVal;
	}

	HSetState((Handle)wavH, hState);

	wavH = MagFieldAccessData.wBz;
	if(wavH == NIL) return NOWAV;
	waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long BzLen = dimensionSizes[0];
	if((numDimensions == 0) || (BzLen == 0)) return BAD_MAGN_FIELD_WAVE_FORMAT;

	DOUBLE sBzStep, sBzStart;
	if(result = MDGetWaveScaling(wavH, 0, &sBzStep, &sBzStart)) return result;

	if((BzLen != BxLen) || (::fabs(sBzStart-sBxStart) > 1.E-05) || (::fabs(sBzStep-sBxStep) > 1.E-05))
		return UNEQUAL_BX_BZ_DEFINITION_LIMITS;

	char BzLenWasModified = (TrjDat.AuxBzLen != BzLen)? 1 : 0;
	if(BzLenWasModified)
	{
		if(BzInData != 0) delete[] BzInData;
		BzInData = new srTFunDer[BzLen];
		if(BzInData == 0) return MEMORY_ALLOCATION_FAILURE;

		LocInputWasNotModified = 0;
		TrjDat.AuxBzLen = BzLen;
	}

	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	hState = 0; //MoveLockHandle(wavH); //OC180815
	dataStartPtr = (char*)(*wavH) + dataOffset;
	dp0 = (DOUBLE*)dataStartPtr;

	srTFunDer* tBzInData = BzInData;
	double BzVal;
	char BzIsZero = 1;
	for(int j=0; j<BzLen; j++)
	{
		BzVal = *(dp0++);

		if(::fabs(BzVal) > FieldZeroTolerance) BzIsZero = 0;
		else BzVal = 0.;

		if(LocInputWasNotModified) if(tBzInData->f != BzVal) LocInputWasNotModified = 0;

		(tBzInData++)->f = BzVal;
	}

	HSetState((Handle)wavH, hState);

	if((TrjDat.LenFieldData != BzLen) || (TrjDat.sStep != sBzStep) || (TrjDat.sStart != sBzStart))
		LocInputWasNotModified = 0;

	TrjDat.InputWasModified = !LocInputWasNotModified;

	TrjDat.LenFieldData = BzLen;
	TrjDat.sStep = sBzStep;
	TrjDat.sStart = sBzStart;

	TrjDat.HorFieldIsNotZero = !BxIsZero;
	TrjDat.VerFieldIsNotZero = !BzIsZero;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalFieldDataFormat2(srTTrjDat& TrjDat)
{// Modify later following General description of the Field wave
#if defined(__IGOR_PRO__)

	waveHndl wField;
	if(pIgorRadInputStruct != 0) wField = pIgorRadInputStruct->wField;
	else if(pIgorTrjInputStruct != 0) wField = pIgorTrjInputStruct->wField;
	else if(pIgorStokesWigInputStruct != 0) wField = pIgorStokesWigInputStruct->wGenMagField;
	else if(pIgorStokesArbInputStruct != 0) wField = pIgorStokesArbInputStruct->wGenMagField;
	else if(pIgorPowDensInputStruct != 0) wField = pIgorPowDensInputStruct->wGenMagField;
	else if(pIgorRadFocusInputStruct != 0) wField = pIgorRadFocusInputStruct->wField;
	else if(pIgorWfrEmitPropagInputStruct != 0) wField = pIgorWfrEmitPropagInputStruct->wMagFld;

	if(wField == NIL) return NOWAV;
	int waveType = WaveType(wField);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	srTMagFieldAccessData MagFieldAccessData;

	int result;
	long RadIndices[MAX_DIMENSIONS];

	Handle textH = NewHandle(0L);
	RadIndices[0] = 0;
	if(result = MDGetTextWavePointValue(wField, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;

	MagFieldAccessData.wBx = FetchWave(*textH);
	DisposeHandle(textH);

	textH = NewHandle(0L);
	RadIndices[0] = 1;
	if(result = MDGetTextWavePointValue(wField, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;

	MagFieldAccessData.wBz = FetchWave(*textH);
	DisposeHandle(textH);

	if(result = GetTotalFieldDataFormat1(MagFieldAccessData, TrjDat)) return result;

// Continue here for the rest field parameters, if needed

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetPeriodicFieldDataFormat1(srTMagFieldPeriodic& MagPer, waveHndl In_wavH)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(In_wavH == 0)
	{
		if(pIgorPerStokesInputStruct != 0) wavH = pIgorPerStokesInputStruct->wPerField;
		else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wGenMagField;
		else if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wGenMagField;
	}
	else wavH = In_wavH;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	int result;
	long numDimensionsM;
	long dimensionSizesM[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensionsM, dimensionSizesM)) return result;
	long NumOfLines = dimensionSizesM[0];

	double BufValue;
	long RadIndices[MAX_DIMENSIONS];
	Handle textH;

	if(result = GetDoubleFromTextWave1D(wavH, 0, MagPer.PerLength)) return result;
	if(result = GetDoubleFromTextWave1D(wavH, 1, MagPer.TotLength)) return result;
	if(result = GetDoubleFromTextWave1D(wavH, 2, BufValue)) return result;

// Type of undulator at the input:
// 1- normal; 2- tapered; 3- optical clystron; 4- infinite; 
// Type of undulator in C code:
// 0- infinite; 1- normal; 2- tapered; 3- optical clystron.
	int TypeOfUndulatorAtInput = int(BufValue + 1.E-12);
	if(TypeOfUndulatorAtInput == 1) MagPer.TypeOfUnd = 1; // Conventional
	else if(TypeOfUndulatorAtInput == 2) MagPer.TypeOfUnd = 2; // Tapered
	else if(TypeOfUndulatorAtInput == 3) MagPer.TypeOfUnd = 3; // Opt. Klystron
	else if(TypeOfUndulatorAtInput == 4) MagPer.TypeOfUnd = 0; // Infinite

	if(MagPer.TypeOfUnd == 2) // Tapered
	{
		if(result = GetDoubleFromTextWave1D(wavH, 3, MagPer.TaperPar_TU)) return result;
	}
	else if(MagPer.TypeOfUnd == 3) // Opt. Klystron
	{
		if(result = GetDoubleFromTextWave1D(wavH, 3, MagPer.PhaseSh_OK)) return result;
	}

	if(result = GetDoubleFromTextWave1D(wavH, 4, MagPer.Fund_keV_per_GeV2)) return result;
	if(result = GetDoubleFromTextWave1D(wavH, 5, BufValue)) return result;
	
	MagPer.AmOfHarm = int(BufValue + 1.E-12);
	if(MagPer.AmOfHarm > 20) return TOO_MANY_MAG_FIELD_HARMONICS;

	MagPer.HarmVect.erase(MagPer.HarmVect.begin(), MagPer.HarmVect.end());

	for(int i=0; i<MagPer.AmOfHarm; i++)
	{
		textH = NewHandle(0L);
		RadIndices[0] = i + 6;
		if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
		waveHndl wHarm = FetchWave(*textH);
		DisposeHandle(textH);
		if(wHarm == NIL) return CAN_NOT_FIND_HARMONIC_WAVE;

		int waveType = WaveType(wHarm);
		if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

		long numDimensions;
		long dimensionSizes[MAX_DIMENSIONS+1];
		if(result = MDGetWaveDimensions(wHarm, &numDimensions, dimensionSizes)) return result;
		if((numDimensions > 1) || (dimensionSizes[0] < 4)) return IMPROPER_FIELD_HARMONIC_STRUCTURE;

		long dataOffset;
		if(result = MDAccessNumericWaveData(wHarm, kMDWaveAccessMode0, &dataOffset)) return result;
		int hState = 0; //MoveLockHandle(wHarm); //OC180815
		char* dataStartPtr = (char*)(*wHarm) + dataOffset;
		DOUBLE* dp = (DOUBLE*)dataStartPtr;

		int HarmNo = int(*(dp++));
		char XorZ = (*(dp++) == 1.)? 'z' : 'x';
		double K = *(dp++);
		double Phase = *dp;

		HSetState((Handle)wHarm, hState);

		srTMagHarm Harm(HarmNo, XorZ, K, Phase);
		MagPer.HarmVect.push_back(Harm);
	}

	if(NumOfLines > 30)
	{
		double BufVal;
		//if(result = GetDoubleFromTextWave1D(wavH, 30, BufVal)) return result; 
		//MagPer.NatFocTypeSASE = int(BufVal + 0.00001 - 1); // 0- natural, 1- other

		if(result = GetDoubleFromTextWave1D(wavH, 30, MagPer.NatFocNxSASE)) return result; 
		if(result = GetDoubleFromTextWave1D(wavH, 31, MagPer.NatFocNySASE)) return result; 

		if(result = GetDoubleFromTextWave1D(wavH, 32, BufVal)) return result; 
		//MagPer.FldErrTypeSASE = int(BufVal + 0.00001 - 1); // 0- No errors; 1- Uniform uncorrelated; 2- Uniform correlated; 3- Gaussian uncorrelated; 4- Gaussian correlated
		MagPer.FldErrTypeSASE = int(BufVal + 0.00001); // 0- No errors; 1- Uniform uncorrelated; 2- Uniform correlated; 3- Gaussian uncorrelated; 4- Gaussian correlated

		if(result = GetDoubleFromTextWave1D(wavH, 33, MagPer.FldErrRMS)) return result;
		if(MagPer.FldErrRMS == 0.) MagPer.FldErrTypeSASE = 0;

		if(result = GetDoubleFromTextWave1D(wavH, 35, BufVal)) return result; 
		MagPer.TaperTypeSASE = int(BufVal + 0.00001); // 0- No taper; 1- Linear; 2- Quadratic

		if(result = GetDoubleFromTextWave1D(wavH, 36, MagPer.TaperStartSASE)) return result; 
		if(result = GetDoubleFromTextWave1D(wavH, 37, MagPer.TaperRelFldChgSASE)) return result;
		if(MagPer.TaperRelFldChgSASE == 0.) MagPer.TaperTypeSASE = 0;
	}

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetConstantFieldDataFormat1(srTMagFieldConstant& MagConst)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wGenMagField;
	else if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wGenMagField;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	int result;
	if(result = GetDoubleFromTextWave1D(wavH, 0, MagConst.Bx)) return result;
	if(result = GetDoubleFromTextWave1D(wavH, 1, MagConst.Bz)) return result;
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::IdentifyMagFieldType(char& MagFieldType)
{ // 1: arbitrary, 2: periodic, 3: constant, 4: mag. optics
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wGenMagField;
	else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wGenMagField;
	else if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wGenMagField;
	else if(pIgorStokesArbInputStruct != 0) wavH = pIgorStokesArbInputStruct->wGenMagField;

	if(wavH == NIL) return NOWAV;

	char NameOfFieldWave[MAX_OBJ_NAME+1];
	WaveName(wavH, NameOfFieldWave);

	MagFieldType = IdentifyMagFieldTypeFromName(NameOfFieldWave);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::IdentifyMagFieldTypeFromName(char* MagFieldName)
{ // 1: arbitrary, 2: periodic, 3: constant, 4: mag. optics
	if(MagFieldName == 0) return 0;
	char* pEnding = strrchr(MagFieldName, '_');
	if(!strcmp(pEnding, "_mag")) return 1;
	else if(!strcmp(pEnding, "_map")) return 2;
	else if(!strcmp(pEnding, "_mac")) return 3;
	else if(!strcmp(pEnding, "_mgo")) return 4;

	else return 0;
}

//*************************************************************************

int srTSend::GetGenMagFieldDataFormat1(srTGenTrjHndl& hGenTrj)
{
#if defined(__IGOR_PRO__)

	char MagFieldType; // 1: arbitrary, 2: periodic, 3: constant
	int result;
	if(result = IdentifyMagFieldType(MagFieldType)) return result;
	if(MagFieldType == 1)
	{
		srTTrjDat* pTrjDat = new srTTrjDat();
		if(pTrjDat == 0) return MEMORY_ALLOCATION_FAILURE;

		if(result = GetTotalFieldDataFormat2(*pTrjDat)) return result;
		hGenTrj = srTGenTrjHndl(pTrjDat);
	}
	else if(MagFieldType == 2)
	{
		srTPerTrjDat* pPerTrjDat = new srTPerTrjDat();
		if(pPerTrjDat == 0) return MEMORY_ALLOCATION_FAILURE;

		if(result = GetPeriodicFieldDataFormat1(pPerTrjDat->MagPer)) return result;
		pPerTrjDat->CheckIfHorOrVertFieldIsZero();
		hGenTrj = srTGenTrjHndl(pPerTrjDat);
	}
	else if(MagFieldType == 3)
	{
		srTConstTrjDat* pConstTrjDat = new srTConstTrjDat();
		if(pConstTrjDat == 0) return MEMORY_ALLOCATION_FAILURE;

		if(result = GetConstantFieldDataFormat1(pConstTrjDat->MagConst)) return result;
		pConstTrjDat->CheckIfHorOrVertFieldIsZero();
		hGenTrj = srTGenTrjHndl(pConstTrjDat);
	}
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetGeneralMagneticFieldDataFormat1(srTMagGroup& MagGen)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wMagGen;
	//else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wGenMagField;
	
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;

	int AmOfElems = int(dimensionSizes[0]*0.5 + 0.0000001);

	for(int i=0; i<AmOfElems; i++)
	{
		srTMagPosAndElem PosAndElem;

		if(result = GetDoubleFromTextWave1D(wavH, i*2, PosAndElem.s)) return result;

		long Indices[MAX_DIMENSIONS];
		Indices[0] = i*2 + 1;
		Handle textH = NewHandle(0L);
		if(result = MDGetTextWavePointValue(wavH, Indices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;

		waveHndl wMagField = FetchWave(*textH);
		int MagTypeNo = IdentifyMagFieldTypeFromName(*textH);
		// 1: arbitrary, 2: periodic, 3: constant, 4: mag. optics
		if(MagTypeNo == 2)
		{
			srTMagElem *pMagElem = new srTMagFieldPeriodic();
			if(pMagElem == 0) return MEMORY_ALLOCATION_FAILURE;
			if(result = GetPeriodicFieldDataFormat1(*((srTMagFieldPeriodic*)pMagElem), wMagField)) return result;
			
			PosAndElem.MagHndl = CHMagFld(pMagElem);
		}
		else if(MagTypeNo == 4)
		{
			srTStringVect OptElemInfo;
			srTMagElemSummary MagElemSummary;
			if(result = GetVectorOfStrings(wMagField, &OptElemInfo)) return result;
			if(result = MagElemSummary.SetupMagElement(&OptElemInfo, PosAndElem.MagHndl)) return result;
		}
		// treat more magnet elements here

		MagGen.PosAndElemVect.push_back(PosAndElem);
		DisposeHandle(textH);
	}
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalElectronBeamDataFormat2(srTTrjDat& TrjDat)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wElectronBeam;
	else if(pIgorTrjInputStruct != 0) wavH = pIgorTrjInputStruct->wElectronBeam;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wElectronBeam;
	else if(pIgorWfrEmitPropagInputStruct != 0) wavH = pIgorWfrEmitPropagInputStruct->wElectronBeam;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 7) return BAD_EL_BEAM_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;

	char LocInputWasNotModified = 1;

	double OldGamma = TrjDat.EbmDat.Gamma;
	double ElectronEnergyInGeV = *(dp++);
	if((ElectronEnergyInGeV < 0.) || (ElectronEnergyInGeV > 1.E+06)) return BAD_ELECTRON_ENERGY;
	TrjDat.EbmDat.SetupGamma(ElectronEnergyInGeV);
	if(TrjDat.EbmDat.Gamma != OldGamma) LocInputWasNotModified = 0;

	double NewElCurrent = *(dp++);
	if(NewElCurrent != TrjDat.EbmDat.Current) LocInputWasNotModified = 0;
	TrjDat.EbmDat.Current = NewElCurrent; 
	
	double New_x0 = *(dp++); // Assuming input in m !

	if(New_x0 != TrjDat.EbmDat.x0) LocInputWasNotModified = 0;
	TrjDat.EbmDat.x0 = New_x0;

	double New_dxds0 = *(dp++); // dxds0 assumed in r at the input !
	if(New_dxds0 != TrjDat.EbmDat.dxds0) LocInputWasNotModified = 0;
	TrjDat.EbmDat.dxds0 = New_dxds0;

	double New_z0 = *(dp++); // Assuming input in m !
	if(New_z0 != TrjDat.EbmDat.z0) LocInputWasNotModified = 0;
	TrjDat.EbmDat.z0 = New_z0;

	double New_dzds0 = *(dp++); // dzds0 assumed in r at the input !
	if(New_dzds0 != TrjDat.EbmDat.dzds0) LocInputWasNotModified = 0;
	TrjDat.EbmDat.dzds0 = New_dzds0;

	double New_s0 = *dp; // Assuming input in m !
	if(New_s0 != TrjDat.EbmDat.s0) LocInputWasNotModified = 0;
	TrjDat.EbmDat.s0 = New_s0;

	TrjDat.InputWasModified = (TrjDat.InputWasModified) || (!LocInputWasNotModified);

	HSetState((Handle)wavH, hState);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalElectronDistribDataFormat1(srTEbmDat& EbmDat)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH = 0;
	//if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wElectronBeam;
	//else if(pIgorTrjInputStruct != 0) wavH = pIgorTrjInputStruct->wElectronBeam;
	//else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wElectronBeam;
	//else if(pIgorPerStokesInputStruct != 0) wavH = pIgorPerStokesInputStruct->wElectronBeam;
	//else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wElectronBeam;
	//else if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wElectronBeam;
	//else if(pIgorStokesArbInputStruct != 0) wavH = pIgorStokesArbInputStruct->wElectronBeam;
	//else if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wElectronBeam;
	//else if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wElectronBeam;
	//else if(pIgorGsnBeamInputStruct != 0) wavH = wavH_AuxElectronBeam;
	//else if(pIgorRadPropagStokesMultiElecInputStruct != 0) wavH = pIgorRadPropagStokesMultiElecInputStruct->wElectronBeam;
	//else if(pIgorWfrReflectInputStruct != 0) wavH = pIgorWfrReflectInputStruct->wElectronBeam;
	if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wElectronBeam;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP32) return NT_FP32_WAVE_REQUIRED;

	EbmDat.wElecDistr = wavH;

	int result = 0;
	long numDimElecDistr;
	long dimSizesElecDistr[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(EbmDat.wElecDistr, &numDimElecDistr, dimSizesElecDistr)) return result;
	EbmDat.nTotMacroPart = dimSizesElecDistr[1];

	long dataOffsetElecDistr;
	if(result = MDAccessNumericWaveData(EbmDat.wElecDistr, kMDWaveAccessMode0, &dataOffsetElecDistr)) return result;
	char* dataStartPtr = (char*)(*(EbmDat.wElecDistr)) + dataOffsetElecDistr;
	EbmDat.pElecDistr = (float*)dataStartPtr;
	EbmDat.hStateElecDistr = 0; //MoveLockHandle(EbmDat.wElecDistr);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalElectronBeamDataFormat3(srTEbmDat& EbmDat)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wElectronBeam;
	else if(pIgorTrjInputStruct != 0) wavH = pIgorTrjInputStruct->wElectronBeam;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wElectronBeam;
	else if(pIgorPerStokesInputStruct != 0) wavH = pIgorPerStokesInputStruct->wElectronBeam;
	else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wElectronBeam;
	else if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wElectronBeam;
	else if(pIgorStokesArbInputStruct != 0) wavH = pIgorStokesArbInputStruct->wElectronBeam;
	else if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wElectronBeam;
	else if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wElectronBeam;
	else if(pIgorGsnBeamInputStruct != 0) wavH = wavH_AuxElectronBeam;
	else if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wElectronBeam;
	else if(pIgorRadPropagStokesMultiElecInputStruct != 0) wavH = pIgorRadPropagStokesMultiElecInputStruct->wElectronBeam;
	else if(pIgorWfrReflectInputStruct != 0) wavH = pIgorWfrReflectInputStruct->wElectronBeam;
	//else if(pIgor3DViewInputStruct != 0) wavH = pIgor3DViewInputStruct->wElectronBeam;

	return GetTotalElectronBeamDataFormat3_FromIgorWave(EbmDat, wavH);

#else

	return 0;

#endif
}

//*************************************************************************

#if defined(__IGOR_PRO__)
int srTSend::GetTotalElectronBeamDataFormat3_FromIgorWave(srTEbmDat& EbmDat, waveHndl wavH)
{
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 7) return BAD_EL_BEAM_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp = dp0;

	char LocInputWasNotModified = 1;

	double OldEnergy = EbmDat.Energy;
	double ElectronEnergyInGeV = *(dp++);
	if((ElectronEnergyInGeV < 0.) || (ElectronEnergyInGeV > 1.E+06)) return BAD_ELECTRON_ENERGY;
	EbmDat.SetupGamma(ElectronEnergyInGeV);
	if(EbmDat.Energy != OldEnergy) LocInputWasNotModified = 0;

	double NewElCurrent = *(dp++);
	if(NewElCurrent != EbmDat.Current) LocInputWasNotModified = 0;
	EbmDat.Current = NewElCurrent; 

// First-order moments
	double New_x0 = *(dp++); // Assuming input in m !
	if(New_x0 != EbmDat.x0) LocInputWasNotModified = 0;
	EbmDat.x0 = New_x0;

	double New_dxds0 = *(dp++); // dxds0 assumed in r at the input !
	if(New_dxds0 != EbmDat.dxds0) LocInputWasNotModified = 0;
	EbmDat.dxds0 = New_dxds0;

	double New_z0 = *(dp++); // Assuming input in m !
	if(New_z0 != EbmDat.z0) LocInputWasNotModified = 0;
	EbmDat.z0 = New_z0;

	double New_dzds0 = *(dp++); // dzds0 assumed in r at the input !
	if(New_dzds0 != EbmDat.dzds0) LocInputWasNotModified = 0;
	EbmDat.dzds0 = New_dzds0;

	double New_s0 = *dp; // Assuming input in m !
	if(New_s0 != EbmDat.s0) LocInputWasNotModified = 0;
	EbmDat.s0 = New_s0;

// Second-order moments
	dp = dp0 + 13;
	double New_Mee = (*dp)*(*dp); 
	if(New_Mee != EbmDat.Mee) LocInputWasNotModified = 0;
	EbmDat.Mee = New_Mee;
	EbmDat.SigmaRelE = *dp;

	dp = dp0 + 20;
	double New_Mxx = *(dp++); // Assuming input in m^2 !
	if(New_Mxx != EbmDat.Mxx) LocInputWasNotModified = 0;
	EbmDat.Mxx = New_Mxx;

	double New_Mxxp = *(dp++); // Assuming input in m !
	if(New_Mxxp != EbmDat.Mxxp) LocInputWasNotModified = 0;
	EbmDat.Mxxp = New_Mxxp;

	double New_Mxpxp = *(dp++); // Assuming input in r^2 !
	if(New_Mxpxp != EbmDat.Mxpxp) LocInputWasNotModified = 0;
	EbmDat.Mxpxp = New_Mxpxp;

	double New_Mzz = *(dp++); // Assuming input in m^2 !
	if(New_Mzz != EbmDat.Mzz) LocInputWasNotModified = 0;
	EbmDat.Mzz = New_Mzz;

	double New_Mzzp = *(dp++); // Assuming input in m !
	if(New_Mzzp != EbmDat.Mzzp) LocInputWasNotModified = 0;
	EbmDat.Mzzp = New_Mzzp;

	double New_Mzpzp = *(dp++); // Assuming input in r^2 !
	if(New_Mzpzp != EbmDat.Mzpzp) LocInputWasNotModified = 0;
	EbmDat.Mzpzp = New_Mzpzp;

	double New_Mxz = *(dp++); // Assuming input in m^2 !
	if(New_Mxz != EbmDat.Mxz) LocInputWasNotModified = 0;
	EbmDat.Mxz = New_Mxz;

	double New_Mxpz = *(dp++); // Assuming input in m !
	if(New_Mxpz != EbmDat.Mxpz) LocInputWasNotModified = 0;
	EbmDat.Mxpz = New_Mxpz;

	double New_Mxzp = *(dp++); // Assuming input in m !
	if(New_Mxzp != EbmDat.Mxzp) LocInputWasNotModified = 0;
	EbmDat.Mxzp = New_Mxzp;

	double New_Mxpzp = *(dp++); // Assuming input in r^2 !
	if(New_Mxpzp != EbmDat.Mxpzp) LocInputWasNotModified = 0;
	EbmDat.Mxpzp = New_Mxpzp;

// More of Transverse & Longitudinal parameters (for SASE, coherent computations)
	if(numRows > 30)
	{
		dp = dp0 + 30;
		double New_TypeDistrTransverse = *(dp++); 
		if(New_TypeDistrTransverse != EbmDat.TypeDistrTransverse) LocInputWasNotModified = 0;
		EbmDat.TypeDistrTransverse = New_TypeDistrTransverse;
		
		double New_TypeDistrLongitudinal = *(dp++); 
		if(New_TypeDistrLongitudinal != EbmDat.TypeDistrLongitudinal) LocInputWasNotModified = 0;
		EbmDat.TypeDistrLongitudinal = New_TypeDistrLongitudinal;
		
		double New_ShotNoiseFactor = *(dp++); 
		if(New_ShotNoiseFactor != EbmDat.ShotNoiseFactor) LocInputWasNotModified = 0;
		EbmDat.ShotNoiseFactor = New_ShotNoiseFactor;
		
		double New_Mss = *(dp++); // Assuming input in m^2 ! //[33]
		if(New_Mss != EbmDat.Mss) LocInputWasNotModified = 0;
		EbmDat.Mss = New_Mss;
	}
	if(numRows > 34)
	{
		dp = dp0 + 34;
		double New_Mse = *(dp++); // Assuming input in m ! //[34]
		if(New_Mse != EbmDat.Mse) LocInputWasNotModified = 0;
		EbmDat.Mse = New_Mse;

		double New_Mxe = *(dp++); // Assuming input in m ! //[35]
		if(New_Mxe != EbmDat.Mxe) LocInputWasNotModified = 0;
		EbmDat.Mxe = New_Mxe;

		double New_Mxpe = *(dp++); // Assuming input in r^2 ! //[36]
		if(New_Mxpe != EbmDat.Mxpe) LocInputWasNotModified = 0;
		EbmDat.Mxpe = New_Mxpe;

		double New_Mze = *(dp++); // Assuming input in m ! //[37]
		if(New_Mze != EbmDat.Mze) LocInputWasNotModified = 0;
		EbmDat.Mze = New_Mze;

		double New_Mzpe = *(dp++); // Assuming input in r^2 ! //[38]
		if(New_Mzpe != EbmDat.Mzpe) LocInputWasNotModified = 0;
		EbmDat.Mzpe = New_Mzpe;

		double New_Mxs = *(dp++); // Assuming input in m^2 ! //[39]
		if(New_Mxs != EbmDat.Mxs) LocInputWasNotModified = 0;
		EbmDat.Mxs = New_Mxs;

		double New_Mxps = *(dp++); // Assuming input in m ! //[40]
		if(New_Mxps != EbmDat.Mxps) LocInputWasNotModified = 0;
		EbmDat.Mxps = New_Mxps;

		double New_Mzs = *(dp++); // Assuming input in m^2 ! //[41]
		if(New_Mzs != EbmDat.Mzs) LocInputWasNotModified = 0;
		EbmDat.Mzs = New_Mzs;

		double New_Mzps = *(dp++); // Assuming input in m ! //[42]
		if(New_Mzps != EbmDat.Mzps) LocInputWasNotModified = 0;
		EbmDat.Mzps = New_Mzps;

		double New_Neb = *(dp++); // Assuming input in m ! //[43]
		if(New_Neb != EbmDat.Neb) LocInputWasNotModified = 0;
		EbmDat.Neb = New_Neb;
	}
	if(numRows > 44)
	{
		dp = dp0 + 44;
		int ElecDistrExists = (int)(*(dp++)); //[44]
		if(ElecDistrExists) 
		{
			const int extWaveElecLen = 4;
			char extWaveElecDistr[] = "_ebd\0"; //to keep in line with srw/igor !!!

			char nameWaveElec[256];
			WaveName(wavH, nameWaveElec);
			int lenNameWaveElec = strlen(nameWaveElec);
			strcpy(nameWaveElec + lenNameWaveElec - extWaveElecLen, extWaveElecDistr);

			EbmDat.wElecDistr = FetchWave(nameWaveElec);
			if(EbmDat.wElecDistr != NULL)
			{
				int waveTypeElecDistr = WaveType(EbmDat.wElecDistr);
				if(waveTypeElecDistr != NT_FP32) return NT_FP32_WAVE_REQUIRED;

				long numDimElecDistr;
				long dimSizesElecDistr[MAX_DIMENSIONS+1];
				if(result = MDGetWaveDimensions(EbmDat.wElecDistr, &numDimElecDistr, dimSizesElecDistr)) return result;
				//EbmDat.nMacroPart = dimSizesElecDistr[0];
				EbmDat.nTotMacroPart = dimSizesElecDistr[1];

				long dataOffsetElecDistr;
				if(result = MDAccessNumericWaveData(EbmDat.wElecDistr, kMDWaveAccessMode0, &dataOffsetElecDistr)) return result;
				char* dataStartPtr = (char*)(*(EbmDat.wElecDistr)) + dataOffsetElecDistr;
				EbmDat.pElecDistr = (float*)dataStartPtr;
				EbmDat.hStateElecDistr = 0; //MoveLockHandle(EbmDat.wElecDistr);
			}
		}
	}
	if(numRows > 45)
	{
		dp = dp0 + 45;
		double NewElCurrentPeak = *(dp++);
		if(NewElCurrentPeak != EbmDat.CurrentPeak) LocInputWasNotModified = 0;
		EbmDat.CurrentPeak = NewElCurrentPeak; 
	}

	EbmDat.InputWasModified = !LocInputWasNotModified;

	HSetState((Handle)wavH, hState);
	return 0;
}
#endif

//*************************************************************************

int srTSend::GetTotalElectronData(srTTrjDat& TrjDat)
{
#ifdef __MATHEMATICA__

	const char *FunName;
	long ArgNum;
	double DoubleValue;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("SR::Error000"); return -1;}
	if((!ReadOK) || strcmp(FunName, "List") || ArgNum!=3) { ErrorMessage("SR::Error000"); return -1;}
	MLDisownSymbol(stdlink, FunName);

	Next = MLGetNext(stdlink);
	double ElectronEnergyInGeV;
	if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) ElectronEnergyInGeV = DoubleValue;}
	else { ErrorMessage("SR::Error000"); return -1;}
	TrjDat.SetupGamma(ElectronEnergyInGeV);

	Next = MLGetNext(stdlink);
	if((Next==MLTKINT) || (Next==MLTKREAL)) { if(MLGetDouble(stdlink, &DoubleValue)) TrjDat.ElCurrent = DoubleValue;}
	else { ErrorMessage("SR::Error000"); return -1;}

	return GetPrtclTrjInitData(TrjDat.s0, TrjDat.x0, TrjDat.dxds0, TrjDat.z0, TrjDat.dzds0);

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetIsotrSrcExtraDataFormat1(srTIsotrSrc& IsotrSrc)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wIsotrSrc;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 5) return BAD_ISOTR_SRC_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp = dp0;

	char LocInputWasNotModified = 1;

	double OldPhotPerBW = IsotrSrc.PhotPerBW;
	double PhotPerBW = *(dp++);
	if(PhotPerBW <= 0.) return BAD_SOURCE_PHOTON_NUM;
	if(PhotPerBW != OldPhotPerBW) LocInputWasNotModified = 0;
	IsotrSrc.PhotPerBW = PhotPerBW;

	int OldPolar = IsotrSrc.Polar;
	int Polar = int(*(dp++));
	if((Polar < 1) || (Polar > 6)) return BAD_POLAR_SPECIFIER;
	if(Polar != OldPolar) LocInputWasNotModified = 0;
	IsotrSrc.Polar = Polar;

	if(!LocInputWasNotModified) IsotrSrc.InputWasModified = 1;
	HSetState((Handle)wavH, hState);

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTotalGsnBeamDataFormat1(srTGsnBeam& GsnBeam)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorGsnBeamInputStruct != 0) wavH = pIgorGsnBeamInputStruct->wGsnBeam;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return IMPROPER_GSN_BEAM_STRUCTURE;

	int result;
	long Indices[MAX_DIMENSIONS];
	Handle textH;

// Electron Beam
	textH = NewHandle(0L);
	Indices[0] = 0;
	if(result = MDGetTextWavePointValue(wavH, Indices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	wavH_AuxElectronBeam = FetchWave(*textH);
	if(wavH_AuxElectronBeam == NIL) return IMPROPER_GSN_BEAM_STRUCTURE;
	if(result = GetTotalElectronBeamDataFormat3(GsnBeam.EbmDat)) return result;
	DisposeHandle(textH);

// Gaussian Beam parameters
	double SigmaX;
	if(result = GetDoubleFromTextWave1D(wavH, 1, SigmaX)) return result;
	if(SigmaX <= 0.) return BAD_WAIST_SIZE;
	GsnBeam.SigmaX = SigmaX;

	double SigmaZ;
	if(result = GetDoubleFromTextWave1D(wavH, 2, SigmaZ)) return result;
	if(SigmaZ <= 0.) return BAD_WAIST_SIZE;
	GsnBeam.SigmaZ = SigmaZ;

	double mx;
	if(result = GetDoubleFromTextWave1D(wavH, 3, mx)) return result;
	if(mx < 0.) return BAD_MODE_ORDER;
	GsnBeam.mx = int(mx);

	double mz;
	if(result = GetDoubleFromTextWave1D(wavH, 4, mz)) return result;
	if(mz < 0.) return BAD_MODE_ORDER;
	GsnBeam.mz = int(mz);

	double Phot;
	if(result = GetDoubleFromTextWave1D(wavH, 5, Phot)) return result;
	if(Phot <= 0.) return BAD_PHOTON_NUMBER;
	GsnBeam.PhotPerBW = Phot;

	double Polar;
	if(result = GetDoubleFromTextWave1D(wavH, 6, Polar)) return result;
	if((Polar < 1.) || (Polar > 6.)) return BAD_POLAR_SPECIFIER;
	GsnBeam.Polar = int(Polar);

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTotalOutTrjParam(short* TypeOfTrjCharact, double* sSt_sEn_Ns, short& Arg)
{
#ifdef __MATHEMATICA__

	const int MaxStrSize = 20;

	char* ComponID = new char[MaxStrSize];
	if(GetString(ComponID, MaxStrSize)) { ErrorMessage("SR::Error006"); return -1;}
	if((!strcmp(ComponID, "Bx")) || (!strcmp(ComponID, "bx")) || (!strcmp(ComponID, "BX"))) TypeOfTrjCharact[0] = 1;
	else if((!strcmp(ComponID, "Bz")) || (!strcmp(ComponID, "bz")) || (!strcmp(ComponID, "BZ"))) TypeOfTrjCharact[1] = 1;
	else if((!strcmp(ComponID, "Btx")) || (!strcmp(ComponID, "btx")) || (!strcmp(ComponID, "BTX"))) TypeOfTrjCharact[2] = 1;
	else if((!strcmp(ComponID, "Btz")) || (!strcmp(ComponID, "btz")) || (!strcmp(ComponID, "BTZ"))) TypeOfTrjCharact[3] = 1;
	else if((!strcmp(ComponID, "x")) || (!strcmp(ComponID, "X"))) TypeOfTrjCharact[4] = 1;
	else if((!strcmp(ComponID, "z")) || (!strcmp(ComponID, "Z"))) TypeOfTrjCharact[5] = 1;
	else if((!strcmp(ComponID, "IntBtxE2")) || (!strcmp(ComponID, "intbtxe2")) || (!strcmp(ComponID, "INTBTXE2"))) TypeOfTrjCharact[6] = 1;
	else if((!strcmp(ComponID, "IntBtzE2")) || (!strcmp(ComponID, "intbtze2")) || (!strcmp(ComponID, "INTBTZE2"))) TypeOfTrjCharact[7] = 1;
	else { ErrorMessage("SR::Error006"); return -1;}

	int Ns;
	if(GetCoordData(sSt_sEn_Ns[0], sSt_sEn_Ns[1], Ns)) return -1;
	sSt_sEn_Ns[2] = Ns;

	if(GetString(ComponID, MaxStrSize)) { ErrorMessage("SR::Error006"); return -1;}
	if((!strcmp(ComponID, "Arg")) || (!strcmp(ComponID, "arg")) || (!strcmp(ComponID, "ARG"))) Arg = 1;
	delete[] ComponID;

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTrajectoryComponentDataFormat1(srTWaveAccessDataD1D& TrjCmpn)
{
#if defined(__IGOR_PRO__)

	int result;
	if(wavH_AuxTrajectory == NIL) return IMPROPER_TRAJECTORY_STRUCTURE;
	int waveType = WaveType(wavH_AuxTrajectory);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	TrjCmpn.wHndl = wavH_AuxTrajectory;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH_AuxTrajectory, &numDimensions, dimensionSizes)) return result;
	if(numDimensions > 1) return IMPROPER_TRAJECTORY_COMPONENT_WAVE;
	TrjCmpn.np = dimensionSizes[0];

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH_AuxTrajectory, kMDWaveAccessMode0, &dataOffset)) return result;
	TrjCmpn.hState = 0; //MoveLockHandle(wavH_AuxTrajectory);
	char* dataStartPtr = (char*)(*wavH_AuxTrajectory) + dataOffset;
	TrjCmpn.pData = (DOUBLE*)dataStartPtr;

	DOUBLE Step, Start;
	if(result = MDGetWaveScaling(wavH_AuxTrajectory, 0, &Step, &Start)) return result;
	TrjCmpn.Step = Step;
	TrjCmpn.Start = Start;

	char Units[MAX_UNIT_CHARS + 1];
	if(result = MDGetWaveUnits(wavH_AuxTrajectory, 0, Units)) return result;
	strcpy(TrjCmpn.ArgUnits, Units);
	if(result = MDGetWaveUnits(wavH_AuxTrajectory, -1, Units)) return result;
	strcpy(TrjCmpn.DataUnits, Units);

	WaveName(wavH_AuxTrajectory, TrjCmpn.NameOfWave);

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTotalTrajectoryDataFormat1(srTTrjDat& TrjDat)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorWfrFromTrjInputStruct != 0) wavH = pIgorWfrFromTrjInputStruct->wTrajectory;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return IMPROPER_TRAJECTORY_STRUCTURE;

	int result;
	long Indices[MAX_DIMENSIONS];
	Handle textH;

// x vs s
	textH = NewHandle(0L);
	Indices[0] = 0;
	if(result = MDGetTextWavePointValue(wavH, Indices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	wavH_AuxTrajectory = FetchWave(*textH);
	if(wavH_AuxTrajectory == NIL) return IMPROPER_TRAJECTORY_STRUCTURE;
	if(result = GetTrajectoryComponentDataFormat1(TrjDat.xTrjInData)) return result;
	DisposeHandle(textH);

// z vs s
	textH = NewHandle(0L);
	Indices[0] = 1;
	if(result = MDGetTextWavePointValue(wavH, Indices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	wavH_AuxTrajectory = FetchWave(*textH);
	if(wavH_AuxTrajectory == NIL) return IMPROPER_TRAJECTORY_STRUCTURE;
	if(result = GetTrajectoryComponentDataFormat1(TrjDat.zTrjInData)) return result;
	DisposeHandle(textH);

// Energy
	double E;
	if(result = GetDoubleFromTextWave1D(wavH, 4, E)) return result;
	if(E <= 0.) return ZERO_OR_NEGATIVE_ELECTRON_ENERGY;
	TrjDat.EbmDat.SetupGamma(E);

// Current
	double ElCurrent;
	if(result = GetDoubleFromTextWave1D(wavH, 5, ElCurrent)) return result;
	TrjDat.EbmDat.Current = ElCurrent; 

// Read more data, if any

// Execute at the end of the input
	if(result = TrjDat.CheckAndSetupTrajectoryLimits()) return result;
	if(result = TrjDat.CheckIfFieldIsZeroFromTrj()) return result;
	if(result = TrjDat.SetupSourcePointFromTrajectory()) return result;

#endif

	return 0;
}

//*************************************************************************

int srTSend::FinishWorkingWithTrajectoryData(srTTrjDat& TrjDat)
{
#if defined(__IGOR_PRO__)

	if(TrjDat.xTrjInData.wHndl != NIL)
	{
		HSetState((Handle)(TrjDat.xTrjInData.wHndl), TrjDat.xTrjInData.hState);
		WaveHandleModified(TrjDat.xTrjInData.wHndl);
	}
	if(TrjDat.zTrjInData.wHndl != NIL)
	{
		HSetState((Handle)(TrjDat.zTrjInData.wHndl), TrjDat.zTrjInData.hState);
		WaveHandleModified(TrjDat.zTrjInData.wHndl);
	}
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalObservationDataFormat1(srTWfrSmp& DistrInfoDat)
{
#if defined(__MATHEMATICA__)

	const char *FunName;
	long ArgNum, InitArgNum;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &InitArgNum);
	else { ErrorMessage("SR::Error000"); return -1;}
	if((!ReadOK) || strcmp(FunName, "List") || (!(InitArgNum==5 || InitArgNum==4))) { ErrorMessage("SR::Error000"); return -1;}
	MLDisownSymbol(stdlink, FunName);

	Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("SR::Error000"); return -1;}
	if((!ReadOK) || strcmp(FunName, "List") || (!(ArgNum==3))) { ErrorMessage("SR::Error000"); return -1;}
	MLDisownSymbol(stdlink, FunName);

	const int MaxStrSize = 50;

	char* DistributionTypeString = new char[MaxStrSize];
	if(GetString(DistributionTypeString, MaxStrSize)) { ErrorMessage("SR::Error000"); return -1;}
	if(!strcmp(DistributionTypeString, "mr"))
	{
		DistrInfoDat.CoordOrAngPresentation = AngPres;
		DistrInfoDat.ArgType = Angle;
		DistrInfoDat.DistrNormType = Angle;
		DistrInfoDat.ShowPhaseOnly = 0;
	}
	else if(!strcmp(DistributionTypeString, "mm"))
	{
		if(InitArgNum==4) { ErrorMessage("SR::Error012"); return -1;}

		DistrInfoDat.CoordOrAngPresentation = CoordPres;
		DistrInfoDat.ArgType = Distance;
		DistrInfoDat.DistrNormType = Distance;
		DistrInfoDat.ShowPhaseOnly = 0;
	}
	else if(!strcmp(DistributionTypeString, "ph") || !strcmp(DistributionTypeString, "Ph"))
	{
		DistrInfoDat.ShowPhaseOnly = 1;
		DistrInfoDat.PhaseDerOrder = 0;
	}
	else if(!strcmp(DistributionTypeString, "dphds") || !strcmp(DistributionTypeString, "dPhds"))
	{
		DistrInfoDat.ShowPhaseOnly = 1;
		DistrInfoDat.PhaseDerOrder = 1;
	}
	else if(!strcmp(DistributionTypeString, "d2hds2") || !strcmp(DistributionTypeString, "d2Phds2"))
	{
		DistrInfoDat.ShowPhaseOnly = 1;
		DistrInfoDat.PhaseDerOrder = 2;
	}
	else { ErrorMessage("SR::Error010"); return -1;}
	delete[] DistributionTypeString;

	char* BandwidthUnitsString = new char[MaxStrSize];
	if(GetString(BandwidthUnitsString, MaxStrSize)) { ErrorMessage("SR::Error000"); return -1;}
	if(!strcmp(BandwidthUnitsString, ".1%")) DistrInfoDat.BandwidthUnits = 0;
	else if(!strcmp(BandwidthUnitsString, "keV")) DistrInfoDat.BandwidthUnits = 1;
	else if(!strcmp(BandwidthUnitsString, "eV")) DistrInfoDat.BandwidthUnits = 2;
	else if(!strcmp(BandwidthUnitsString, "Ang")) DistrInfoDat.BandwidthUnits = 3;
	else if(!strcmp(BandwidthUnitsString, "nm")) DistrInfoDat.BandwidthUnits = 4;
	else if(!strcmp(BandwidthUnitsString, "mic")) DistrInfoDat.BandwidthUnits = 5;
	else { ErrorMessage("SR::Error014"); return -1;}
	delete[] BandwidthUnitsString;

	if(DistrInfoDat.BandwidthUnits > 0) { ErrorMessage("SR::Error101"); return -1;} // To remove later

	DistributionTypeString = new char[MaxStrSize];
	if(GetString(DistributionTypeString, MaxStrSize)) { ErrorMessage("SR::Error000"); return -1;}
	if(!strcmp(DistributionTypeString, "Stokes")) 
		DistrInfoDat.DistrValType = StokesParam;
	else if(!strcmp(DistributionTypeString, "Field")) DistrInfoDat.DistrValType = FieldFourier;
	else { ErrorMessage("SR::Error015"); return -1;}
	delete[] DistributionTypeString;

	if(GetWavelengthAndCoordData(DistrInfoDat, (InitArgNum==5)? 1 : 0)) return -1;
	return 0;

#elif defined(__IGOR_PRO__)

// Trace later DistrInfoDat.RadDistrDataContShouldBeRebuild !

	int LocInputWasModified = 0, LocRadDistrDataContShouldBeRebuild = 0;

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wObservation;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wObservation;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 5) return BAD_OBS_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	int DistrTypeID = int(*(dp0++));
	if(DistrTypeID == 0)
	{
		if(DistrInfoDat.CoordOrAngPresentation != AngPres)
		{
			LocInputWasModified = 1;
		}

		DistrInfoDat.CoordOrAngPresentation = AngPres;
		DistrInfoDat.ArgType = Angle;
		DistrInfoDat.DistrNormType = Angle;
		DistrInfoDat.ShowPhaseOnly = 0;
	}
	else
	{
		if(DistrInfoDat.CoordOrAngPresentation != CoordPres)
		{
			LocInputWasModified = 1;
		}

		DistrInfoDat.CoordOrAngPresentation = CoordPres;
		DistrInfoDat.ArgType = Distance;
		DistrInfoDat.DistrNormType = Distance;
		DistrInfoDat.ShowPhaseOnly = 0;
	}
	if(numRows == 5)
	{
		char* pLoopOrder = DistrInfoDat.LoopOrder;
		*(pLoopOrder++) = 'y';
		*(pLoopOrder++) = 'w';
		*(pLoopOrder++) = 'x';
		*(pLoopOrder++) = 'z';

		DistrInfoDat.DistrValType = StokesParam;

		DistrInfoDat.yStart = DistrInfoDat.yEnd = 1.E+23;
		DistrInfoDat.ny = 1;

		DistrInfoDat.PhotonEnergyWavelengthUnits = 1; // We support only eV for the moment
		double Lambda = *(dp0++); // Assuming input in eV !

		if(Lambda < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;

		if(DistrInfoDat.SetupLambda(Lambda, Lambda)) LocInputWasModified = 1;
		DistrInfoDat.nLamb = 1;

		double X = *(dp0++);
		if((X != DistrInfoDat.xStart) || (X != DistrInfoDat.xEnd)) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.xStart = DistrInfoDat.xEnd = X;
		}
		DistrInfoDat.nx = 1;

		double Z = *(dp0++);
		if((Z != DistrInfoDat.zStart) || (Z != DistrInfoDat.zEnd))
		{
			LocInputWasModified = 1;
			DistrInfoDat.zStart = DistrInfoDat.zEnd = Z;
		}
		DistrInfoDat.nz = 1;

		if(DistrTypeID != 0)
		{
			double Y = *(dp0++); // Assuming input in m !
			if((Y != DistrInfoDat.yStart) || (Y != DistrInfoDat.yEnd))
			{
				LocInputWasModified = 1;
				DistrInfoDat.yStart = DistrInfoDat.yEnd = Y;
			}
			DistrInfoDat.ny = 1;
		}
		DistrInfoDat.OnlyOnePoint = 1; // Modify later if necessary
	}
	else
	{
		char* pLoopOrder = DistrInfoDat.LoopOrder;
		*(pLoopOrder++) = 'y';
		*(pLoopOrder++) = 'z';
		*(pLoopOrder++) = 'x';
		*(pLoopOrder++) = 'w';

		DistrInfoDat.DistrValType = FieldFourier;

		double ySt = *(dp0++); // Assuming input in m !
		if(ySt != DistrInfoDat.yStart) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.yStart = ySt;
		}
		double yFi = *(dp0++); // Assuming input in m !
		if(yFi != DistrInfoDat.yEnd) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.yEnd = yFi;
		}
		int ny = int(*(dp0++));
		if(ny != DistrInfoDat.ny) 
		{
			LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
			DistrInfoDat.ny = ny;
		}

		DistrInfoDat.PhotonEnergyWavelengthUnits = 1; // We support only eV for the moment
		double LambdaSt = *(dp0++); // Assuming input in eV !
		if(LambdaSt < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;
		double LambdaFi = *(dp0++); // Assuming input in eV !
		if(LambdaFi < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;
		if(DistrInfoDat.SetupLambda(LambdaSt, LambdaFi)) LocInputWasModified = 1;
		int nLamb = int(*(dp0++));
		if(nLamb != DistrInfoDat.nLamb) 
		{
			LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
			DistrInfoDat.nLamb = nLamb;
		}

		double xSt = *(dp0++); 
		if(xSt != DistrInfoDat.xStart) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.xStart = xSt;
		}
		double xFi = *(dp0++);
		if(xFi != DistrInfoDat.xEnd) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.xEnd = xFi;
		}
		int nx = int(*(dp0++));
		if(nx != DistrInfoDat.nx) 
		{
			LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
			DistrInfoDat.nx = nx;
		}

		double zSt = *(dp0++);
		if(zSt != DistrInfoDat.zStart) 
		{
			LocInputWasModified = 1;
			DistrInfoDat.zStart = zSt;
		}
		double zFi = *(dp0++);
		if(zFi != DistrInfoDat.zEnd) 
		{

			LocInputWasModified = 1;
			DistrInfoDat.zEnd = zFi;
		}
		int nz = int(*(dp0++));
		if(nz != DistrInfoDat.nz) 
		{
			LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
			DistrInfoDat.nz = nz;
		}

		DistrInfoDat.OnlyOnePoint = (ny*nLamb*nx*nz > 1)? 0 : 1; // Modify later if necessary
		DistrInfoDat.RadDistrDataContShouldBeRebuild = (DistrInfoDat.RadDistrDataContShouldBeRebuild || LocRadDistrDataContShouldBeRebuild);
	}

	HSetState((Handle)wavH, hState);
	DistrInfoDat.InputWasModified = LocInputWasModified;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalObservationDataFormat3(srTWfrSmp& DistrInfoDat)
{
#if defined(__IGOR_PRO__)

	int LocInputWasModified = 0, LocRadDistrDataContShouldBeRebuild = 0;

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wObservation;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wObservation;
	else if(pIgorPerStokesInputStruct != 0) wavH = pIgorPerStokesInputStruct->wObservation;
	else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wObservation;
	else if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wObservation;
	else if(pIgorStokesArbInputStruct != 0) wavH = pIgorStokesArbInputStruct->wObservation;
	else if(pIgorPowDensInputStruct != 0) wavH = pIgorPowDensInputStruct->wObservation;
	else if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wObservation;
	else if(pIgorGsnBeamInputStruct != 0) wavH = pIgorGsnBeamInputStruct->wObservation;
	else if(pIgorWfrFromTrjInputStruct != 0) wavH = pIgorWfrFromTrjInputStruct->wObservation;
	else if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wObservation;
	else if(pIgorWfrEmitPropagInputStruct != 0) wavH = pIgorWfrEmitPropagInputStruct->wObservation;

	//if(wavH == NIL) return NOWAV;
	if(wavH == NIL) return 0; //allow to be absent

	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 5) return BAD_OBS_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp00 = dp0;

	dp0 += 4; // The fields 0 - 3 are not used
	
	char* pLoopOrder = DistrInfoDat.LoopOrder;
	*(pLoopOrder++) = 'y';
	*(pLoopOrder++) = 'z';
	*(pLoopOrder++) = 'x';
	*(pLoopOrder++) = 'w';

	DistrInfoDat.CoordOrAngPresentation = CoordPres;

	double ySt = *(dp0++); 
	if(ySt != DistrInfoDat.yStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.yStart = ySt;
	}
	double yFi = ySt; // Assuming input in m !
	if(yFi != DistrInfoDat.yEnd) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.yEnd = yFi;
	}
	int ny = 1;
	if(ny != DistrInfoDat.ny) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.ny = ny;
	}

	DistrInfoDat.PhotonEnergyWavelengthUnits = 1; // We support only eV for the moment
	DistrInfoDat.TreatLambdaAsEnergyIn_eV = 1;

	double LambdaSt = *(dp0++); // Assuming input in eV !
	if(LambdaSt < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;
	double LambdaFi = *(dp0++); // Assuming input in eV !
	if(LambdaFi < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;
	if(DistrInfoDat.SetupLambda(LambdaSt, LambdaFi)) LocInputWasModified = 1;
	int nLamb = int(*(dp0++));
	if(nLamb != DistrInfoDat.nLamb) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.nLamb = nLamb;
	}

	double xSt = *(dp0++); // Assuming input in m!
	double xFi = *(dp0++); // Assuming input in m!
	int nx = int(*(dp0++));
	if(nx == 1)
	{
		if((!DistrInfoDat.FluxComp) && (!DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat))
		{
			double xMid = 0.5*(xSt + xFi);
			xSt = xMid; xFi = xMid;
		}
	}

	if(xSt != DistrInfoDat.xStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.xStart = xSt;
	}
	if(xFi != DistrInfoDat.xEnd) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.xEnd = xFi;
	}
	if(nx != DistrInfoDat.nx) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.nx = nx;
	}

	double zSt = *(dp0++); // Assuming input in m!
	double zFi = *(dp0++); // Assuming input in m!
	int nz = int(*(dp0++));
	if(nz == 1)
	{
		if((!DistrInfoDat.FluxComp) && (!DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat))
		{
			double zMid = 0.5*(zSt + zFi);
			zSt = zMid; zFi = zMid;
		}
	}

	if(zSt != DistrInfoDat.zStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.zStart = zSt;
	}
	if(zFi != DistrInfoDat.zEnd) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.zEnd = zFi;
	}
	if(nz != DistrInfoDat.nz) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.nz = nz;
	}

	if(numRows > 14)
	{
		DistrInfoDat.PresT = *(dp00 + 14);
		DistrInfoDat.tStart = *(dp00 + 15);
		DistrInfoDat.tEnd = *(dp00 + 16);
		DistrInfoDat.nt = (int)(*(dp00 + 17));
	}

	HSetState((Handle)wavH, hState);
	DistrInfoDat.InputWasModified = LocInputWasModified;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetTotalObservationDataFormat2(srTWfrSmp& DistrInfoDat)
{// For transformed Light; May be changed
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	int LocInputWasModified = 0, LocRadDistrDataContShouldBeRebuild = 0;

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wObservation;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wObservation;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 5) return BAD_OBS_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	DistrInfoDat.PhotonEnergyWavelengthUnits = 1; // We support only eV for the moment
	double Lambda = *(dp0++); // Assuming input in eV !
	if(Lambda < 0.) return WAVELENGTH_SHOULD_BE_POSITIVE;
	if(DistrInfoDat.SetupLambda(Lambda, Lambda)) LocInputWasModified = 1;
	DistrInfoDat.nLamb = 1;

	double xCP = *(dp0++); // Assuming input in m !
	if(xCP != DistrInfoDat.CenP.x)
	{
		LocInputWasModified = 1;
		DistrInfoDat.CenP.x = xCP;
	}
	double yCP = *(dp0++); // Assuming input in m !
	if(yCP != DistrInfoDat.CenP.y)
	{
		LocInputWasModified = 1;
		DistrInfoDat.CenP.y = yCP;
	}
	double zCP = *(dp0++); // Assuming input in m !
	if(zCP != DistrInfoDat.CenP.z)
	{
		LocInputWasModified = 1;
		DistrInfoDat.CenP.z = zCP;
	}

	double xSt = *(dp0++); 
	if(xSt != DistrInfoDat.xStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.xStart = xSt;
	}
	double xFi = *(dp0++);
	if(xFi != DistrInfoDat.xEnd) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.xEnd = xFi;
	}
	int nx = int(*(dp0++));
	if(nx != DistrInfoDat.nx) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.nx = nx;
	}

	double ySt = *(dp0++); 
	if(ySt != DistrInfoDat.yStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.yStart = ySt;
	}
	double yFi = *(dp0++);
	if(yFi != DistrInfoDat.yEnd)
	{
		LocInputWasModified = 1;
		DistrInfoDat.yEnd = yFi;
	}
	int ny = int(*(dp0++));
	if(ny != DistrInfoDat.ny) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.ny = ny;
	}

	double zSt = *(dp0++);
	if(zSt != DistrInfoDat.zStart) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.zStart = zSt;
	}
	double zFi = *(dp0++);
	if(zFi != DistrInfoDat.zEnd) 
	{
		LocInputWasModified = 1;
		DistrInfoDat.zEnd = zFi;
	}
	int nz = int(*(dp0++));
	if(nz != DistrInfoDat.nz) 
	{
		LocInputWasModified = 1; LocRadDistrDataContShouldBeRebuild = 1;
		DistrInfoDat.nz = nz;
	}

	HSetState((Handle)wavH, hState);

	DistrInfoDat.OnlyOnePoint = (nx*nz > 1)? 0 : 1; // Modify later if necessary
	DistrInfoDat.RadDistrDataContShouldBeRebuild = (DistrInfoDat.RadDistrDataContShouldBeRebuild || LocRadDistrDataContShouldBeRebuild);
	DistrInfoDat.InputWasModified = LocInputWasModified;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetWavelengthAndCoordData(srTWfrSmp& DistrInfoDat, char GetLongCoord)
{
#ifdef __MATHEMATICA__

	const char *FunName;
	long ArgNum;
	int ReadOK = 0, Next;

	const int MaxStrSize = 50;
	char* IDString = new char[MaxStrSize];

	char* pLoopOrder = DistrInfoDat.LoopOrder;
	if(!GetLongCoord) 
	{
		*(pLoopOrder++) = 'y';
		DistrInfoDat.yStart = DistrInfoDat.yEnd = 1.E+23;
		DistrInfoDat.ny = 1;
	}

	int AmOfPar = GetLongCoord? 4 : 3;
	for(int i=0; i<AmOfPar; i++)
	{
		Next = MLGetNext(stdlink);
		if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		else { ErrorMessage("SR::Error000"); return -1;}
		if((!ReadOK) || strcmp(FunName, "List") || (!(ArgNum==2))) { ErrorMessage("SR::Error000"); return -1;}
		MLDisownSymbol(stdlink, FunName);

		if(GetString(IDString, MaxStrSize)) { ErrorMessage("SR::Error000"); return -1;}

		if((!strcmp(IDString, "keV")) || (!strcmp(IDString, "eV")) || (!strcmp(IDString, "Ang")) || (!strcmp(IDString, "nm")) || (!strcmp(IDString, "mic"))) 
		{
			*(pLoopOrder++) = 'w';
			double LambStart, LambEnd;
			if(GetCoordData(LambStart, LambEnd, DistrInfoDat.nLamb)) return -1;

			if(!strcmp(IDString, "keV")) DistrInfoDat.PhotonEnergyWavelengthUnits = 0;
			else if(!strcmp(IDString, "eV")) DistrInfoDat.PhotonEnergyWavelengthUnits = 1;
			else if(!strcmp(IDString, "Ang")) DistrInfoDat.PhotonEnergyWavelengthUnits = 2;
			else if(!strcmp(IDString, "nm")) DistrInfoDat.PhotonEnergyWavelengthUnits = 3;
			else if(!strcmp(IDString, "mic")) DistrInfoDat.PhotonEnergyWavelengthUnits = 4;

			DistrInfoDat.SetupLambda(LambStart, LambEnd);
		}
		else if((!strcmp(IDString, "x")) || (!strcmp(IDString, "X")))
		{
			*(pLoopOrder++) = 'x';
			if(GetCoordData(DistrInfoDat.xStart, DistrInfoDat.xEnd, DistrInfoDat.nx)) return -1;
		}
		else if((!strcmp(IDString, "z")) || (!strcmp(IDString, "Z")))
		{
			*(pLoopOrder++) = 'z';
			if(GetCoordData(DistrInfoDat.zStart, DistrInfoDat.zEnd, DistrInfoDat.nz)) return -1;
		}
		else
		{
			if(GetLongCoord)
			{
				if((!strcmp(IDString, "long")) || (!strcmp(IDString, "Long")) || (!strcmp(IDString, "LONG")))
				{
					*(pLoopOrder++) = 'y';
					if(GetCoordData(DistrInfoDat.yStart, DistrInfoDat.yEnd, DistrInfoDat.ny)) return -1;
				}
			}
		}
	}
	*pLoopOrder = '\0';
	delete[] IDString;

#endif

	return 0;
}

//*************************************************************************

int srTSend::GetTotalRadIntegrationParamDataFormat1(srTRadInt& RadInt)
{
#if defined(__MATHEMATICA__)

	const char *FunName;
	long ArgNum;
	int ReadOK = 0;

	int Next = MLGetNext(stdlink);
	if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
	else { ErrorMessage("SR::Error000"); return -1;}
	if((!ReadOK) || strcmp(FunName, "List") || (!(ArgNum==4 || ArgNum==3))) { ErrorMessage("SR::Error000"); return -1;}
	MLDisownSymbol(stdlink, FunName);

	const int MaxStrSize = 50;
	char* RadIntegrMethodString = new char[MaxStrSize];
	if(GetString(RadIntegrMethodString, MaxStrSize)) return -1;
	if(!strcmp(RadIntegrMethodString, "man1") || !strcmp(RadIntegrMethodString, "Man1") || !strcmp(RadIntegrMethodString, "MAN1"))
		RadInt.sIntegMethod = 0;
	else if(!strcmp(RadIntegrMethodString, "man2") || !strcmp(RadIntegrMethodString, "Man2") || !strcmp(RadIntegrMethodString, "MAN2"))
		RadInt.sIntegMethod = 1;
	else if(!strcmp(RadIntegrMethodString, "man3") || !strcmp(RadIntegrMethodString, "Man3") || !strcmp(RadIntegrMethodString, "MAN3"))
		RadInt.sIntegMethod = 2;
	else if(!strcmp(RadIntegrMethodString, "auto1") || !strcmp(RadIntegrMethodString, "Auto1") || !strcmp(RadIntegrMethodString, "AUTO1")) 
		RadInt.sIntegMethod = 10;
	else if(!strcmp(RadIntegrMethodString, "auto2") || !strcmp(RadIntegrMethodString, "Auto2") || !strcmp(RadIntegrMethodString, "AUTO2")) 
		RadInt.sIntegMethod = 11;
	else { ErrorMessage("SR::Error013"); return -1;}
	delete[] RadIntegrMethodString;

	if(!MLGetDouble(stdlink, &(RadInt.sIntegRelPrec))) { ErrorMessage("SR::Error000"); return -1;}
	if(!MLGetDouble(stdlink, &(RadInt.sIntegStep))) { ErrorMessage("SR::Error000"); return -1;}

	if(ArgNum==3)
	{
		RadInt.sIntegStart = RadInt.TrjDatPtr->sStart;
		RadInt.sIntegFin = RadInt.sIntegStart + RadInt.TrjDatPtr->sStep*(RadInt.TrjDatPtr->LenFieldData - 1);
	}
	else
	{
		if(Next==MLTKFUNC) ReadOK = MLGetFunction(stdlink, &FunName, &ArgNum);
		else { ErrorMessage("SR::Error000"); return -1;}
		if((!ReadOK) || strcmp(FunName, "List") || !(ArgNum==2)) { ErrorMessage("SR::Error000"); return -1;}
		MLDisownSymbol(stdlink, FunName);

		if(!MLGetDouble(stdlink, &(RadInt.sIntegStart))) { ErrorMessage("SR::Error000"); return -1;}
		if(!MLGetDouble(stdlink, &(RadInt.sIntegFin))) { ErrorMessage("SR::Error000"); return -1;}
	}

	if(RadInt.sIntegMethod < 10)
	{
		if((::fabs(RadInt.sIntegFin - RadInt.sIntegStart) < RadInt.sIntegStep) || (RadInt.sIntegStep < 0.)) 
		{ ErrorMessage("SR::Error016"); return -1;}
	}
	return 0;

#elif defined(__IGOR_PRO__)

	waveHndl wavH;

	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wIntPar;
	else if(pIgorRadFocusInputStruct != 0) wavH = pIgorRadFocusInputStruct->wIntPar;
	else if(pIgorWfrFromTrjInputStruct != 0) wavH = pIgorWfrFromTrjInputStruct->wIntPar;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 2) return BAD_RAD_INT_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	int IntegMeth = int(*(dp0++));
	if((IntegMeth < 0) || (IntegMeth > 2)) return BAD_RAD_INT_METH_VALUE;

	int Old_sIntegMethod = RadInt.sIntegMethod;
	switch(IntegMeth) 
	{
		case 0:
			RadInt.sIntegMethod = 1; break;
		case 1:
			RadInt.sIntegMethod = 10; break;
		case 2:
			RadInt.sIntegMethod = 11; break;
	}
	if(Old_sIntegMethod != RadInt.sIntegMethod) RadInt.TrjDataContShouldBeRebuild = 1;

	double RelPrecOrStep = *(dp0++);

	if(RelPrecOrStep < 0.) return BAD_PREC_OR_STEP_VALUE;
	if(RadInt.sIntegMethod < 10) 
	{
		if(RadInt.sIntegStep_Input != RelPrecOrStep) RadInt.TrjDataContShouldBeRebuild = 1;
		RadInt.sIntegStep = RelPrecOrStep;
		RadInt.sIntegStep_Input = RelPrecOrStep;
	}
	else
	{
		if(RadInt.sIntegRelPrec != RelPrecOrStep) RadInt.TrjDataContShouldBeRebuild = 1;
		RadInt.sIntegRelPrec = RelPrecOrStep;
	}

	double Old_sIntegStart = RadInt.sIntegStart, Old_sIntegFin = RadInt.sIntegFin;
	double Trj_sFin = RadInt.TrjDatPtr->sStart + RadInt.TrjDatPtr->sStep*(RadInt.TrjDatPtr->LenFieldData - 1);
	if(numRows == 2)
	{
		RadInt.sIntegStart = RadInt.TrjDatPtr->sStart;
		RadInt.sIntegFin = Trj_sFin;
	}
	else
	{
		double sStart = *(dp0++); // Assuming input in m !
		double sEnd = *(dp0++); // Assuming input in m !

		if(::fabs(sEnd - sStart) < 1.E-13)
		{
			RadInt.sIntegStart = RadInt.TrjDatPtr->sStart;
			RadInt.sIntegFin = Trj_sFin;
		}
		else
		{
			if((sStart < RadInt.TrjDatPtr->sStart) || (sEnd > Trj_sFin)) return INCONSISTENT_RAD_INT_LIMITS;
			RadInt.sIntegStart = sStart;
			RadInt.sIntegFin = sEnd;
		}
	}
	if((Old_sIntegStart != RadInt.sIntegStart) || (Old_sIntegFin != RadInt.sIntegFin))
		RadInt.TrjDataContShouldBeRebuild = 1;

	long MaxNumPoToSave = long(*(dp0++));
	if(MaxNumPoToSave != RadInt.MaxNumPoToSave)
	{
		RadInt.MaxNumPoToSave = MaxNumPoToSave;
		RadInt.TrjDataContShouldBeRebuild = 1;
	}

// Testing
	if(numRows == 6) // Switch for Trying to apply Near-Field residual term if the far-field one fails
	{
		RadInt.TryToApplyNearFieldResidual = char(*(dp0++));
	}
// End Testing

	HSetState((Handle)wavH, hState);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetPropagRadStokesMultiElecDataFormat1(double* pPrecPar)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorRadPropagStokesMultiElecInputStruct->wPpecPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp = dp0;

	pPrecPar[0] = *(dp++);
	pPrecPar[1] = *(dp++);

	HSetState((Handle)wavH, hState);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetWfrEmitPropagPrec(double* pPrecPar)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorWfrEmitPropagInputStruct->wPpecPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp = dp0;

	pPrecPar[0] = *(dp++);
	pPrecPar[1] = *(dp++);
	//dddddddddddddddddddddddd
	//add more, or remove some params

	HSetState((Handle)wavH, hState);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetRadIntPeriodicParamDataFormat1(srTRadIntPerStoPrec& IntPerPrec)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorPerStokesInputStruct->wPerIntPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;
	DOUBLE* dp = dp0;

	IntPerPrec.InitHarm = int(*(dp++));
	IntPerPrec.FinHarm = int(*(dp++));
	IntPerPrec.Kns = *(dp++);
	IntPerPrec.Knphi = *(dp++);

	IntPerPrec.IntensityOrFlux = (*dp == 1.)? 'f' : 'i';

	HSetState((Handle)wavH, hState);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetRadIntWigglerParamDataFormat1(srTRadIntWigPrec& IntWigPrec)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorStokesWigInputStruct->wWigIntPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;
	IntWigPrec.PrecFact = *(dp++);

	IntWigPrec.TreatInterf = (char)(*dp) - 1; // "No;Yes" in IGOR

	HSetState((Handle)wavH, hState);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetRadIntConstParamDataFormat1(srTRadIntConstPrec& IntConstPrec)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorStokesConstInputStruct->wPrecPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;
	IntConstPrec.PrecFact = *dp;

	HSetState((Handle)wavH, hState);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetRadIntPowDensParamDataFormat1(srTRadIntPowDenPrec& PowDenPrec)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	wavH = pIgorPowDensInputStruct->wPowDensIntPar;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;
	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;

	PowDenPrec.PrecFact = *(dp++);
	PowDenPrec.Method = *dp;

	HSetState((Handle)wavH, hState);

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************
/**Obsolete
int srTSend::GetTotalOpticalElemDataFormat1(void* OpticsHndlPtr)
{
	srTOpticsHndl& OpticsHndl = *((srTOpticsHndl*)OpticsHndlPtr);

#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	waveHndl wavH = pIgorRadFocusInputStruct->wOpticalElem;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	//long numRows = dimensionSizes[0];

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	int OpticalElemID = int(*(dp0++));

	srTOptics* OldSphMirPtr = OpticsHndl.rep;
	char LocNotFirstInput = (OldSphMirPtr != 0);
	char LocInputWasModified = LocNotFirstInput? 0 : 1;

	double xCP = *(dp0++); // Assuming input in m !
	if(LocNotFirstInput) if(xCP != OldSphMirPtr->CentrPoint.x) LocInputWasModified = 1;
	double yCP = *(dp0++); // Assuming input in m !
	if(LocNotFirstInput) if(yCP != OldSphMirPtr->CentrPoint.y) LocInputWasModified = 1;
	double zCP = *(dp0++); // Assuming input in m !
	if(LocNotFirstInput) if(zCP != OldSphMirPtr->CentrPoint.z) LocInputWasModified = 1;
	TVector3d CP(xCP, yCP, zCP);
	double xSN = *(dp0++);
	if(LocNotFirstInput) if(xSN != OldSphMirPtr->CentrSurfNorm.x) LocInputWasModified = 1;
	double ySN = *(dp0++);
	if(LocNotFirstInput) if(ySN != OldSphMirPtr->CentrSurfNorm.y) LocInputWasModified = 1;
	double zSN = *(dp0++);
	if(LocNotFirstInput) if(zSN != OldSphMirPtr->CentrSurfNorm.z) LocInputWasModified = 1;
	TVector3d SN(xSN, ySN, zSN);

	if(OpticalElemID == 20) // Spherical Mirror
	{
		double CurvRad = *(dp0++); // Assuming input in m !
		if(LocNotFirstInput) if(CurvRad != ((srTSphericalMirror*)OldSphMirPtr)->CurvRad) LocInputWasModified = 1;

		double Diameter = *dp0; // Assuming input in m !
		if(LocNotFirstInput) if(Diameter != ((srTSphericalMirror*)OldSphMirPtr)->Diameter) LocInputWasModified = 1;

		HSetState((Handle)wavH, hState);

		if(LocInputWasModified)
		{
			srTSphericalMirror* SphMirPtr = new srTSphericalMirror(CP, SN, CurvRad, Diameter);
			if(SphMirPtr == 0) return MEMORY_ALLOCATION_FAILURE;

			srTOpticsHndl SphMirHndl(SphMirPtr);
			OpticsHndl = SphMirHndl;
		}
		else OldSphMirPtr->InputWasModified = 0;
	}
	else if(OpticalElemID == 30) // Zone Plate
	{
		// Continue here
	}
	else return UNKNOWN_OPTICAL_ELEMENT;

	return 0;

#else

	return 0;

#endif
}
**/
//*************************************************************************

int srTSend::OutOpticsIncRadDistrFormat1(double* DataPtr, int lenData, double xSt, double xFi)
{
#if defined(__MATHEMATICA__)

	return 0;

#elif defined(__IGOR_PRO__)

	waveHndl wavH;
	char waveName[MAX_OBJ_NAME+1];
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	strcpy(waveName, "IncRadTestWave");
	dimensionSizes[0] = lenData;
	dimensionSizes[1] = dimensionSizes[2] = dimensionSizes[3] = 0;
	//if(result = MDMakeWave(&wavH, waveName, NIL, dimensionSizes, NT_FP64, 1)) return result;
	if(result = MDMakeWave(&wavH, waveName, NIL, dimensionSizes, NT_FP64, 1)) return FAILED_TO_CREATE_WAVE;

	DOUBLE xStart = xSt;
	DOUBLE xStep = (lenData > 1)? (xFi - xSt)/(lenData - 1) : 0.;

	if(result = MDSetWaveScaling(wavH, ROWS, &xStep, &xStart)) return result;
	//if(result = MDSetWaveUnits(wavH, ROWS, "m")) return result;
	char strUnit[] = "m\0";
	if(result = MDSetWaveUnits(wavH, ROWS, strUnit)) return result;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	for(int k=0; k<lenData; k++) *(dp0++) = *(DataPtr++);

	HSetState((Handle)wavH, hState);
	WaveHandleModified(wavH);			// Inform Igor that we have changed the wave.

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetSRWRadStructAccessData(srTSRWRadStructAccessData* pSRWRadStructAccessData)
{
#ifdef __IGOR_PRO__

	waveHndl wavH = HandleOfSRWRadStruct.wRad;
	if(wavH == NIL) 
	{
		if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wRad;
		else if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wRad;
		else if(pIgorGsnBeamInputStruct != 0) wavH = pIgorGsnBeamInputStruct->wRad;
		else if(pIgorWfrFromTrjInputStruct != 0) wavH = pIgorWfrFromTrjInputStruct->wRad;
	}

	if(wavH == NIL) return NOWAV;

	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return IMPROPER_RADIATION_STRUCTURE;

	pSRWRadStructAccessData->wRad = wavH;

	int result;
	long dataOffset, numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long RadIndices[MAX_DIMENSIONS];

// Hor. E data
	Handle textH = NewHandle(0L);
	RadIndices[0] = 0;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wRadX = FetchWave(*textH);
	if(wRadX == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->wRadX = wRadX;
	if(result = MDAccessNumericWaveData(wRadX, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hStateRadX = 0; //MoveLockHandle(wRadX);
	char* dataStartPtr = (char*)(*wRadX) + dataOffset;
	pSRWRadStructAccessData->pBaseRadX = (float*)dataStartPtr;
	if(result = MDGetWaveDimensions(wRadX, &numDimensions, dimensionSizes)) return result;
	if(numDimensions != 3) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->ne = dimensionSizes[0];
	pSRWRadStructAccessData->nx = dimensionSizes[1];
	pSRWRadStructAccessData->nz = dimensionSizes[2];
	DOUBLE eStep, eStart, xStep, xStart, zStep, zStart;
	if(result = MDGetWaveScaling(wRadX, 0, &eStep, &eStart)) return result;
	if(result = MDGetWaveScaling(wRadX, 1, &xStep, &xStart)) return result;
	if(result = MDGetWaveScaling(wRadX, 2, &zStep, &zStart)) return result;
	pSRWRadStructAccessData->eStep = eStep; pSRWRadStructAccessData->eStart = eStart;
	pSRWRadStructAccessData->xStep = xStep; pSRWRadStructAccessData->xStart = xStart;
	pSRWRadStructAccessData->zStep = zStep; pSRWRadStructAccessData->zStart = zStart;
	DisposeHandle(textH);

// Vert. E data
	textH = NewHandle(0L);
	RadIndices[0] = 1;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wRadZ = FetchWave(*textH);
	if(wRadZ == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->wRadZ = wRadZ;
	if(result = MDAccessNumericWaveData(wRadZ, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hStateRadZ = 0; //MoveLockHandle(wRadZ);
	dataStartPtr = (char*)(*wRadZ) + dataOffset;
	pSRWRadStructAccessData->pBaseRadZ = (float*)dataStartPtr;
	if(result = MDGetWaveDimensions(wRadZ, &numDimensions, dimensionSizes)) return result;
	if(numDimensions != 3) return IMPROPER_RADIATION_STRUCTURE;
	if(pSRWRadStructAccessData->ne != dimensionSizes[0]) return IMPROPER_RADIATION_STRUCTURE;
	if(pSRWRadStructAccessData->nx != dimensionSizes[1]) return IMPROPER_RADIATION_STRUCTURE;
	if(pSRWRadStructAccessData->nz != dimensionSizes[2]) return IMPROPER_RADIATION_STRUCTURE;
	if(result = MDGetWaveScaling(wRadZ, 0, &eStep, &eStart)) return result;
	if(result = MDGetWaveScaling(wRadZ, 1, &xStep, &xStart)) return result;
	if(result = MDGetWaveScaling(wRadZ, 2, &zStep, &zStart)) return result;
	if((pSRWRadStructAccessData->eStep != eStep) || (pSRWRadStructAccessData->eStart != eStart)) return IMPROPER_RADIATION_STRUCTURE;
	if((pSRWRadStructAccessData->xStep != xStep) || (pSRWRadStructAccessData->xStart != xStart)) return IMPROPER_RADIATION_STRUCTURE;
	if((pSRWRadStructAccessData->zStep != zStep) || (pSRWRadStructAccessData->zStart != zStart)) return IMPROPER_RADIATION_STRUCTURE;
	DisposeHandle(textH);

// Presentation of rad.
	textH = NewHandle(0L);
	RadIndices[0] = 2;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH == '0') pSRWRadStructAccessData->Pres = 0;
			else if(**textH == '1') pSRWRadStructAccessData->Pres = 1;
			else return IMPROPER_RADIATION_STRUCTURE;
		}
	}
	DisposeHandle(textH);

	pSRWRadStructAccessData->LengthUnit = 0; // Length is in m
	pSRWRadStructAccessData->PhotEnergyUnit = 0; // Photon energy is in eV

// Presentation of rad. (freq. / time)
	textH = NewHandle(0L);
	RadIndices[0] = 10;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH == '0') pSRWRadStructAccessData->PresT = 0;
			else if(**textH == '1') pSRWRadStructAccessData->PresT = 1;
			else return IMPROPER_RADIATION_STRUCTURE;
		}
	}
	DisposeHandle(textH);

// Average photon energy (for time repres.)
	textH = NewHandle(0L);
	RadIndices[0] = 11;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH != '\0') pSRWRadStructAccessData->avgPhotEn = atof(*textH);
		}
	}
	DisposeHandle(textH);

// Electron Beam or Trajectory
	textH = NewHandle(0L);
	RadIndices[0] = 12;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wElecBeamOrTrj = FetchWave(*textH);
	if(wElecBeamOrTrj == NIL) return IMPROPER_RADIATION_STRUCTURE;

	int ElecBeamOrTrjWaveType = WaveType(wElecBeamOrTrj);
	if(ElecBeamOrTrjWaveType == TEXT_WAVE_TYPE) // Assume Trajectory, later modify if necessary
	{
		pSRWRadStructAccessData->wTrj = wElecBeamOrTrj;
		pSRWRadStructAccessData->hStateTrj = 0; //MoveLockHandle(wElecBeamOrTrj);
	}
	else // Assume Trajectory, later modify if necessary
	{
		pSRWRadStructAccessData->wElecBeam = wElecBeamOrTrj;
		if(result = MDAccessNumericWaveData(pSRWRadStructAccessData->wElecBeam, kMDWaveAccessMode0, &dataOffset)) return result;
		pSRWRadStructAccessData->hStateElecBeam = 0; //MoveLockHandle(pSRWRadStructAccessData->wElecBeam);
		dataStartPtr = (char*)(*(pSRWRadStructAccessData->wElecBeam)) + dataOffset;
		pSRWRadStructAccessData->pElecBeam = (DOUBLE*)dataStartPtr;
	}
	DisposeHandle(textH);

// 4x4 Matrix for e-beam Moments propagation
	textH = NewHandle(0L);
	RadIndices[0] = 13;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl w4x4PropMatr = FetchWave(*textH);
	if(w4x4PropMatr == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->w4x4PropMatr = w4x4PropMatr;
	if(result = MDAccessNumericWaveData(w4x4PropMatr, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hState4x4PropMatr = 0; //MoveLockHandle(w4x4PropMatr);
	dataStartPtr = (char*)(*w4x4PropMatr) + dataOffset;
	pSRWRadStructAccessData->p4x4PropMatr = (DOUBLE*)dataStartPtr;
	DisposeHandle(textH);

// Moments of Hor. polarization rad.
	textH = NewHandle(0L);
	RadIndices[0] = 15;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wRadMomX = FetchWave(*textH);
	if(wRadMomX == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->wMomX = wRadMomX;
	if(result = MDAccessNumericWaveData(wRadMomX, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hStateMomX = 0; //MoveLockHandle(wRadMomX);
	dataStartPtr = (char*)(*wRadMomX) + dataOffset;
	//pSRWRadStructAccessData->pMomX = (float*)dataStartPtr;
	pSRWRadStructAccessData->pMomX = (DOUBLE*)dataStartPtr; //130311
	DisposeHandle(textH);

// Moments of Vert. polarization rad.
	textH = NewHandle(0L);
	RadIndices[0] = 16;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wRadMomZ = FetchWave(*textH);
	if(wRadMomZ == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->wMomZ = wRadMomZ;
	if(result = MDAccessNumericWaveData(wRadMomZ, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hStateMomZ = 0; //MoveLockHandle(wRadMomZ);
	dataStartPtr = (char*)(*wRadMomZ) + dataOffset;
	//pSRWRadStructAccessData->pMomZ = (float*)dataStartPtr;
	pSRWRadStructAccessData->pMomZ = (DOUBLE*)dataStartPtr; //130311
	DisposeHandle(textH);

// Auxiliary Wave Front Data
	textH = NewHandle(0L);
	RadIndices[0] = 18;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	waveHndl wWfrAuxData = FetchWave(*textH);
	if(wWfrAuxData == NIL) return IMPROPER_RADIATION_STRUCTURE;
	pSRWRadStructAccessData->wWfrAuxData = wWfrAuxData;
	if(result = MDAccessNumericWaveData(wWfrAuxData, kMDWaveAccessMode0, &dataOffset)) return result;
	pSRWRadStructAccessData->hStateWfrAuxData = 0; //MoveLockHandle(wWfrAuxData);
	dataStartPtr = (char*)(*wWfrAuxData) + dataOffset;
	pSRWRadStructAccessData->pWfrAuxData = (DOUBLE*)dataStartPtr;
	DisposeHandle(textH);

	pSRWRadStructAccessData->SetupSrwWfrAuxData();

// Electric Field Units (to allow arb. units)
	textH = NewHandle(0L);
	RadIndices[0] = 19;
	if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH == '0') pSRWRadStructAccessData->ElecFldUnit = 0;
			else if(**textH == '1') pSRWRadStructAccessData->ElecFldUnit = 1;
			else return IMPROPER_RADIATION_STRUCTURE;
		}
	}
	DisposeHandle(textH);

// Read any new elements of Rad here !!!
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetSRWRadStructArray(srTSRWRadStructAccessData*& arSRWRadStructAccessData, int& numRadStruct)
{
#ifdef __IGOR_PRO__

	waveHndl origWavH = HandleOfSRWRadStruct.wRad;
	if(origWavH == NIL) return NOWAV;
	srTStringVect stringVect;

	int result = 0;
	if(result = GetVectorOfStrings(origWavH, &stringVect)) return result;
	numRadStruct = stringVect.size();

	arSRWRadStructAccessData = new srTSRWRadStructAccessData[numRadStruct];

	for(int i=0; i<numRadStruct; i++)
	{
		//srTHandleOfSRWRadStruct locHandleOfSRWRadStruct;
		//locHandleOfSRWRadStruct.wRad = FetchWave(stringVect[i]);
		//SetHandleOfSRWRadStruct(&locHandleOfSRWRadStruct);
		HandleOfSRWRadStruct.wRad = FetchWave(stringVect[i]);

		if(result = GetSRWRadStructAccessData(arSRWRadStructAccessData + i)) return result;
	}
	HandleOfSRWRadStruct.wRad = origWavH;

	DeleteStringVect(stringVect);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetSRWRadStructAndResizeData(srTIgorRadResizeInputStruct* pIgorRadResizeInputStruct, srTSRWRadStructAccessData* pOutSRWRadStructAccessData, srTRadResize* pOutRadResize)
{
	if(pIgorRadResizeInputStruct->wRad == 0) return IMPROPER_RADIATION_STRUCTURE;

	HandleOfSRWRadStruct.wRad = pIgorRadResizeInputStruct->wRad;
	int result;
	if(result = GetSRWRadStructAccessData(pOutSRWRadStructAccessData)) return result;

	pOutRadResize->pxm = pIgorRadResizeInputStruct->pxm;
	pOutRadResize->pxd = pIgorRadResizeInputStruct->pxd;
	pOutRadResize->pzm = pIgorRadResizeInputStruct->pzm;
	pOutRadResize->pzd = pIgorRadResizeInputStruct->pzd;

	//pOutRadResize->UseOtherSideFFT = (pIgorRadResizeInputStruct->MethNo == 0)? 0 : 1;
	pOutRadResize->useOtherSideFFT((pIgorRadResizeInputStruct->MethNo == 0)? 0 : 1);
	return 0;

//#else

	//todo
	//return 0;

}
#endif

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetSRWRadStructAndExtractData(srTIgorRadExtractInputStruct* pIgorRadExtractInputStruct, srTSRWRadStructAccessData* pRadAccessData, srTRadExtract* pRadExtract)
{
	HandleOfSRWRadStruct.wRad = pIgorRadExtractInputStruct->wRad;
	int result;
	if(result = GetSRWRadStructAccessData(pRadAccessData)) return result;

	long dataOffset, numDimensions, dimensionSizes[MAX_DIMENSIONS+1];
	int waveType = WaveType(pIgorRadExtractInputStruct->wExtractParam);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;
	if(result = MDGetWaveDimensions(pIgorRadExtractInputStruct->wExtractParam, &numDimensions, dimensionSizes)) return result;
	if(numDimensions != 1) return NEEDS_1D_WAVE;

	if(result = MDAccessNumericWaveData(pIgorRadExtractInputStruct->wExtractParam, kMDWaveAccessMode0, &dataOffset)) return result;
	int hStateExtractParam = 0; //MoveLockHandle(pIgorRadExtractInputStruct->wExtractParam);
	char* dataStartPtr = (char*)(*(pIgorRadExtractInputStruct->wExtractParam)) + dataOffset;
	DOUBLE* pD0 = (DOUBLE*)(dataStartPtr);
	DOUBLE* pD = pD0;
	pRadExtract->PolarizCompon = int(*(pD++));
	pRadExtract->Int_or_Phase = int(*(pD++));
	pRadExtract->PlotType = int(*(pD++));
	pRadExtract->TransvPres = int(*(pD++));
	pD = pD0 + 10;
	pRadExtract->ePh = *(pD++);
	pRadExtract->x = *(pD++);
	pRadExtract->z = *pD;

	HSetState((Handle)(pIgorRadExtractInputStruct->wExtractParam), hStateExtractParam);
	WaveHandleModified(pIgorRadExtractInputStruct->wExtractParam);

	waveType = WaveType(pIgorRadExtractInputStruct->wExtractedData);
	if((waveType != NT_FP32) && (waveType != NT_FP64)) return NT_FP32_WAVE_REQUIRED;

	if(result = MDAccessNumericWaveData(pIgorRadExtractInputStruct->wExtractedData, kMDWaveAccessMode0, &dataOffset)) return result;
	pRadExtract->wExtractedData = pIgorRadExtractInputStruct->wExtractedData;
	pRadExtract->hStateExtractedData = 0; //MoveLockHandle(pIgorRadExtractInputStruct->wExtractedData);
	dataStartPtr = (char*)(*(pIgorRadExtractInputStruct->wExtractedData)) + dataOffset;
	if(pRadExtract->Int_or_Phase != 2)
	{
		pRadExtract->pExtractedData = (float*)dataStartPtr;
		pRadExtract->pExtractedDataD = 0;
	}
	else
	{
		pRadExtract->pExtractedData = 0;
		pRadExtract->pExtractedDataD = (DOUBLE*)dataStartPtr;
	}
	return 0;

//#else

	//todo
	//return 0;
}
#endif

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetWaveAccessData(srTHandleOfOneWaveStruct* pHandleOfOneWaveStruct, srTWaveAccessData* pWaveAccessData)
{
	waveHndl wavH = pHandleOfOneWaveStruct->WaveHndl;
	if(wavH == NIL) return NOWAV;

	pWaveAccessData->wHndl = wavH;

	int waveType = WaveType(wavH);
	if(waveType == NT_FP32)
	{
		*(pWaveAccessData->WaveType) = 'f';
		*(pWaveAccessData->WaveType + 1) = '\0';
	}
	else if(waveType == (NT_FP32 | NT_CMPLX))
	{
		*(pWaveAccessData->WaveType) = 'c';
		*(pWaveAccessData->WaveType + 1) = 'f';
	}
	else if(waveType == NT_FP64)
	{
		*(pWaveAccessData->WaveType) = 'd';
		*(pWaveAccessData->WaveType + 1) = '\0';
	}
	else if(waveType == (NT_FP64 | NT_CMPLX))
	{
		*(pWaveAccessData->WaveType) = 'c';
		*(pWaveAccessData->WaveType + 1) = 'd';
	}
	else return NUMERIC_WAVE_REQUIRED;

	int result;
	if(result = MDGetWaveDimensions(wavH, &(pWaveAccessData->AmOfDims), pWaveAccessData->DimSizes)) return result;

	DOUBLE Step, Start;
	char Units[MAX_UNIT_CHARS + 1];
	for(long iDim=0; iDim<pWaveAccessData->AmOfDims; iDim++)
	{
		if(result = MDGetWaveScaling(wavH, iDim, &Step, &Start)) return result;
		pWaveAccessData->DimSteps[iDim] = Step; 
		pWaveAccessData->DimStartValues[iDim] = Start; 

		if(result = MDGetWaveUnits(wavH, iDim, Units)) return result;
		strcpy(pWaveAccessData->DimUnits[iDim], Units);
	}
	if(result = MDGetWaveUnits(wavH, -1, Units)) return result;
	strcpy(pWaveAccessData->DataUnits, Units);

	WaveName(wavH, pWaveAccessData->NameOfWave);

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	pWaveAccessData->hState = 0; //MoveLockHandle(wavH);
	pWaveAccessData->pWaveData = (char*)(*wavH) + dataOffset;

	return 0;

//#else

	//todo
	//return 0;
}
#endif

//*************************************************************************

int srTSend::CreateNewRadStruct(srTSRWRadStructAccessData& RadData, srTSRWRadStructWaveNames& Names)
{//ATTENTION!!! It can only create:
#ifdef __IGOR_PRO__

//- Rad. text wave
//- Ex wave
//- Ez wave
//Program more when mecessary !!!
	int result;

	long RadIndices[MAX_DIMENSIONS];
	long dimensionSizes[MAX_DIMENSIONS+1];
	long dataOffset;
	char CharBuf[256];
	int AmOfBytes;

//Text Wave
	if((RadData.wRad == NIL) && (*(Names.NameRad) != '\0'))
	{
		for(int i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
		dimensionSizes[0] = 22; // Don't forget to increase when needed !!!
		//if(result = MDMakeWave(&(RadData.wRad), Names.NameRad, NIL, dimensionSizes, TEXT_WAVE_TYPE, 1)) return result;
		if(result = MDMakeWave(&(RadData.wRad), Names.NameRad, NIL, dimensionSizes, TEXT_WAVE_TYPE, 1)) return FAILED_TO_CREATE_WAVE;
	}

	RadIndices[0] = 2; // Presentation
	AmOfBytes = sprintf(CharBuf, "%d", char(RadData.Pres));
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;

	DOUBLE eStep = RadData.eStep, eStart = RadData.eStart;
	RadIndices[0] = 3; // eStep
	AmOfBytes = sprintf(CharBuf, "%g", eStep);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 4; // eStart
	AmOfBytes = sprintf(CharBuf, "%g", eStart);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;

	DOUBLE xStep = RadData.xStep, xStart = RadData.xStart;
	DOUBLE zStep = RadData.zStep, zStart = RadData.zStart;
	//DOUBLE RobsX = RadData.RobsX, RobsZ = RadData.RobsZ;
	//DOUBLE RobsXAbsErr = RadData.RobsXAbsErr, RobsZAbsErr = RadData.RobsZAbsErr;
	//DOUBLE xc = RadData.xc, zc = RadData.zc;

	RadIndices[0] = 5; // xStep
	AmOfBytes = sprintf(CharBuf, "%g", xStep);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 6; // xStart
	AmOfBytes = sprintf(CharBuf, "%g", xStart);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 7; // zStep
	AmOfBytes = sprintf(CharBuf, "%g", zStep);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 8; // zStart
	AmOfBytes = sprintf(CharBuf, "%g", zStart);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 9; // Quality of Propag.
	AmOfBytes = sprintf(CharBuf, "%i", 1);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;

	RadIndices[0] = 10; // Representation (Freq. / Time)
	AmOfBytes = sprintf(CharBuf, "%d", char(RadData.PresT));
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;
	
	RadIndices[0] = 11; // Average photon energy (for time repres.)
	AmOfBytes = sprintf(CharBuf, "%g", RadData.avgPhotEn);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;

	RadIndices[0] = 14; // Linearity of Transf.
	AmOfBytes = sprintf(CharBuf, "%i", 1);
	if(result = UpdateNumberPositionInSRWRad(&RadData, RadIndices, CharBuf, AmOfBytes)) return result;

	char TransvUnits[MAX_UNIT_CHARS + 1];
	if(RadData.Pres == 0) *TransvUnits = 'm';
	else *TransvUnits = 'q';
	TransvUnits[1] = '\0';

	char strPhotEnOrTimeUnit[] = "eV\0";
	if(RadData.PresT == 1) strcpy(strPhotEnOrTimeUnit, "s");

//Ex
	if((RadData.wRadX == NIL) && (*(Names.NameRadX) != '\0'))
	{
		for(int i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
		dimensionSizes[0] = RadData.ne;
		dimensionSizes[1] = RadData.nx;
		dimensionSizes[2] = RadData.nz;

		//if(result = MDMakeWave(&(RadData.wRadX), Names.NameRadX, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return result;
		if(result = MDMakeWave(&(RadData.wRadX), Names.NameRadX, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return FAILED_TO_CREATE_WAVE;

		if(result = MDSetWaveUnits(RadData.wRadX, 0, strPhotEnOrTimeUnit)) return result;
		if(result = MDSetWaveUnits(RadData.wRadX, 1, TransvUnits)) return result;
		if(result = MDSetWaveUnits(RadData.wRadX, 2, TransvUnits)) return result;

		if(result = MDSetWaveScaling(RadData.wRadX, 0, &eStep, &eStart)) return result;
		if(result = MDSetWaveScaling(RadData.wRadX, 1, &xStep, &xStart)) return result;
		if(result = MDSetWaveScaling(RadData.wRadX, 2, &zStep, &zStart)) return result;

		if(result = MDAccessNumericWaveData(RadData.wRadX, kMDWaveAccessMode0, &dataOffset)) return result;
		RadData.hStateRadX = 0; //MoveLockHandle(RadData.wRadX);
		RadData.pBaseRadX = (float*)((char*)(*(RadData.wRadX)) + dataOffset);
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 0, Names.NameRadX)) return result;

//Ez
	if((RadData.wRadZ == NIL) && (*(Names.NameRadZ) != '\0'))
	{
		for(int i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
		dimensionSizes[0] = RadData.ne;
		dimensionSizes[1] = RadData.nx;
		dimensionSizes[2] = RadData.nz;

		//if(result = MDMakeWave(&(RadData.wRadZ), Names.NameRadZ, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return result;
		if(result = MDMakeWave(&(RadData.wRadZ), Names.NameRadZ, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return FAILED_TO_CREATE_WAVE;

		if(result = MDSetWaveUnits(RadData.wRadZ, 0, strPhotEnOrTimeUnit)) return result;
		if(result = MDSetWaveUnits(RadData.wRadZ, 1, TransvUnits)) return result;
		if(result = MDSetWaveUnits(RadData.wRadZ, 2, TransvUnits)) return result;

		if(result = MDSetWaveScaling(RadData.wRadZ, 0, &eStep, &eStart)) return result;
		if(result = MDSetWaveScaling(RadData.wRadZ, 1, &xStep, &xStart)) return result;
		if(result = MDSetWaveScaling(RadData.wRadZ, 2, &zStep, &zStart)) return result;

		if(result = MDAccessNumericWaveData(RadData.wRadZ, kMDWaveAccessMode0, &dataOffset)) return result;
		RadData.hStateRadZ = 0; //MoveLockHandle(RadData.wRadZ);
		RadData.pBaseRadZ = (float*)((char*)(*(RadData.wRadZ)) + dataOffset);
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 1, Names.NameRadZ)) return result;

//E-beam
	if((RadData.wElecBeam == NIL) && (*(Names.NameElecBeam) != '\0'))
	{
		// Program when necessary
	}
	if(*(Names.NameElecBeam) != '\0')
	{
		if(result = UpdateTextPositionInSRWRad(&RadData, 12, Names.NameElecBeam)) return result;
	}

// Trajectory
	if((RadData.wTrj == NIL) && (*(Names.NameTrj) != '\0'))
	{
		// Program when necessary
	}
	if(*(Names.NameTrj) != '\0')
	{
		if(result = UpdateTextPositionInSRWRad(&RadData, 12, Names.NameTrj)) return result;
	}

// 4x4 Matrix for e-beam Moments propagation
	if((RadData.w4x4PropMatr == NIL) && (*(Names.Name4x4PropMatr) != '\0'))
	{
		// Program when necessary
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 13, Names.Name4x4PropMatr)) return result;

// Radiation Moments
	if((RadData.wMomX == NIL) && (*(Names.NameMomX) != '\0'))
	{
		// Program when necessary
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 15, Names.NameMomX)) return result;

	if((RadData.wMomZ == NIL) && (*(Names.NameMomZ) != '\0'))
	{
		// Program when necessary
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 16, Names.NameMomZ)) return result;

// Auxiliary Wave Front Data
	if((RadData.wWfrAuxData == NIL) && (*(Names.NameWfrAuxData) != '\0'))
	{
		// Program when necessary
	}
	if(result = UpdateTextPositionInSRWRad(&RadData, 18, Names.NameWfrAuxData)) return result;

// Add here treatment of new Rad data members, if any.

	WaveHandleModified(RadData.wRad);
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::DeleteRadStructWaves(srTSRWRadStructAccessData& RadData, srTSRWRadStructWaveKeys& Keys)
{
#ifdef __IGOR_PRO__

	int result;
	if(Keys.wRadX_ && (RadData.wRadX != NIL))
	{
		HSetState((Handle)(RadData.wRadX), RadData.hStateRadX);
		if(result = KillWave(RadData.wRadX)) return result;
	}
	if(Keys.wRadZ_ && (RadData.wRadZ != NIL))
	{
		HSetState((Handle)(RadData.wRadZ), RadData.hStateRadZ);
		if(result = KillWave(RadData.wRadZ)) return result;
	}
	if(Keys.wElecBeam_ && (RadData.wElecBeam != NIL))
	{
		HSetState((Handle)(RadData.wElecBeam), RadData.hStateElecBeam);
		if(result = KillWave(RadData.wElecBeam)) return result;
	}
	if(Keys.wTrj_ && (RadData.wTrj != NIL))
	{
		HSetState((Handle)(RadData.wTrj), RadData.hStateTrj);
		if(result = KillWave(RadData.wTrj)) return result;
	}
	if(Keys.w4x4PropMatr_ && (RadData.w4x4PropMatr != NIL))
	{
		HSetState((Handle)(RadData.w4x4PropMatr), RadData.hState4x4PropMatr);
		if(result = KillWave(RadData.w4x4PropMatr)) return result;
	}
	if(Keys.wMomX_ && (RadData.wMomX != NIL))
	{
		HSetState((Handle)(RadData.wMomX), RadData.hStateMomX);
		if(result = KillWave(RadData.wMomX)) return result;
	}
	if(Keys.wMomZ_ && (RadData.wMomZ != NIL))
	{
		HSetState((Handle)(RadData.wMomZ), RadData.hStateMomZ);
		if(result = KillWave(RadData.wMomZ)) return result;
	}
	if(Keys.wWfrAuxData_ && (RadData.wWfrAuxData != NIL))
	{
		HSetState((Handle)(RadData.wWfrAuxData), RadData.hStateWfrAuxData);
		if(result = KillWave(RadData.wWfrAuxData)) return result;
	}
	if(Keys.wRad_ && (RadData.wRad != NIL))
	{
		if(result = KillWave(RadData.wRad)) return result;
	}

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::ModifyRadNeNxNz(srTSRWRadStructAccessData& NewRadStructAccessData, char PolarizComp)
{
#ifdef __IGOR_PRO__

	if(NewRadStructAccessData.BaseRadWasEmulated) return 0;

	char TreatPolCompX = ((PolarizComp == 0) || (PolarizComp == 'x'));
	char TreatPolCompZ = ((PolarizComp == 0) || (PolarizComp == 'z'));

	long dimensionSizes[MAX_DIMENSIONS+1];
	int i;
	for(i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
	dimensionSizes[0] = NewRadStructAccessData.ne;
	dimensionSizes[1] = NewRadStructAccessData.nx;
	dimensionSizes[2] = NewRadStructAccessData.nz;

	DOUBLE StepE = NewRadStructAccessData.eStep, StepX = NewRadStructAccessData.xStep, StepZ = NewRadStructAccessData.zStep;
	DOUBLE StartE = NewRadStructAccessData.eStart, StartX = NewRadStructAccessData.xStart, StartZ = NewRadStructAccessData.zStart;

	char TransvUnits[MAX_UNIT_CHARS + 1];
	if(NewRadStructAccessData.Pres == 0) *TransvUnits = 'm';
	else *TransvUnits = 'q';
	TransvUnits[1] = '\0';

	char strPhotEnOrTimeUnit[] = "eV\0";
	if(NewRadStructAccessData.PresT == 1) strcpy(strPhotEnOrTimeUnit, "s");

	char AuxRadWaveName[MAX_OBJ_NAME+1];
	int result;
	long dataOffset;

	if(TreatPolCompX)
	{
		HSetState((Handle)(NewRadStructAccessData.wRadX), NewRadStructAccessData.hStateRadX);

		WaveName(NewRadStructAccessData.wRadX, AuxRadWaveName);
		result = KillWave(NewRadStructAccessData.wRadX);
		//if(result = MDMakeWave(&(NewRadStructAccessData.wRadX), AuxRadWaveName, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return result;
		if(result = MDMakeWave(&(NewRadStructAccessData.wRadX), AuxRadWaveName, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return FAILED_TO_CREATE_WAVE;

		//if(result = MDSetWaveUnits(NewRadStructAccessData.wRadX, 0, "eV")) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadX, 0, strPhotEnOrTimeUnit)) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadX, 1, TransvUnits)) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadX, 2, TransvUnits)) return result;

		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadX, 0, &StepE, &StartE)) return result;
		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadX, 1, &StepX, &StartX)) return result;
		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadX, 2, &StepZ, &StartZ)) return result;

		if(result = MDAccessNumericWaveData(NewRadStructAccessData.wRadX, kMDWaveAccessMode0, &dataOffset)) return result;
		NewRadStructAccessData.hStateRadX = 0; //MoveLockHandle(NewRadStructAccessData.wRadX);
		NewRadStructAccessData.pBaseRadX = (float*)((char*)(*(NewRadStructAccessData.wRadX)) + dataOffset);
	}
	if(TreatPolCompZ)
	{
		HSetState((Handle)(NewRadStructAccessData.wRadZ), NewRadStructAccessData.hStateRadZ);

		WaveName(NewRadStructAccessData.wRadZ, AuxRadWaveName);
		result = KillWave(NewRadStructAccessData.wRadZ);
		//if(result = MDMakeWave(&(NewRadStructAccessData.wRadZ), AuxRadWaveName, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return result;
		if(result = MDMakeWave(&(NewRadStructAccessData.wRadZ), AuxRadWaveName, NIL, dimensionSizes, (NT_FP32 | NT_CMPLX), 1)) return FAILED_TO_CREATE_WAVE;

		//if(result = MDSetWaveUnits(NewRadStructAccessData.wRadZ, 0, "eV")) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadZ, 0, strPhotEnOrTimeUnit)) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadZ, 1, TransvUnits)) return result;
		if(result = MDSetWaveUnits(NewRadStructAccessData.wRadZ, 2, TransvUnits)) return result;

		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadZ, 0, &StepE, &StartE)) return result;
		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadZ, 1, &StepX, &StartX)) return result;
		if(result = MDSetWaveScaling(NewRadStructAccessData.wRadZ, 2, &StepZ, &StartZ)) return result;

		if(result = MDAccessNumericWaveData(NewRadStructAccessData.wRadZ, kMDWaveAccessMode0, &dataOffset)) return result;
		NewRadStructAccessData.hStateRadZ = 0; //MoveLockHandle(NewRadStructAccessData.wRadZ);
		NewRadStructAccessData.pBaseRadZ = (float*)((char*)(*(NewRadStructAccessData.wRadZ)) + dataOffset);
	}

	long RadIndices[MAX_DIMENSIONS];
	int AmOfBytes;
	char CharBuf[15];
	Handle textH;

	RadIndices[0] = 3; // eStep
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.eStep);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);

	RadIndices[0] = 4; // eStart
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.eStart);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);

	RadIndices[0] = 5; // xStep
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.xStep);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);
	
	RadIndices[0] = 6; // xStart
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.xStart);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);

	RadIndices[0] = 7; // zStep
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.zStep);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);
	
	RadIndices[0] = 8; // zStart
	for(i=0; i<15; i++) CharBuf[i] = 0;
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.zStart);
	textH = NewHandle(AmOfBytes);
	strncpy(*textH, CharBuf, AmOfBytes);
	if(result = MDSetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
	DisposeHandle(textH);

	RadIndices[0] = 10; // Representation (Freq. / Time)
	AmOfBytes = sprintf(CharBuf, "%d", char(NewRadStructAccessData.PresT));
	if(result = UpdateNumberPositionInSRWRad(&NewRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;
	
	RadIndices[0] = 11; // Average photon energy (for time repres.)
	AmOfBytes = sprintf(CharBuf, "%g", NewRadStructAccessData.avgPhotEn);
	if(result = UpdateNumberPositionInSRWRad(&NewRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;

	//Updating Statistical Rad. Moments, if necessary (OC250507)
	bool needUpdateMomX = true, needUpdateMomZ = true;
	long numDimMom, dimSizesMom[MAX_DIMENSIONS+1];
	for(i=0; i<=MAX_DIMENSIONS; i++) dimSizesMom[i] = 0;
	dimSizesMom[0] = 11;
	dimSizesMom[1] = NewRadStructAccessData.ne;

	if((NewRadStructAccessData.wMomX == 0) || (NewRadStructAccessData.pMomX == 0)) needUpdateMomX = false;
	//do update only if moments wave already exist
	else
	{
		textH = NewHandle(0L);
		RadIndices[0] = 15;
		if(result = MDGetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
		waveHndl wRadMomX = FetchWave(*textH);
		if(wRadMomX != NIL)
		{
			if(result = MDGetWaveDimensions(wRadMomX, &numDimMom, dimSizesMom)) return result;
			if(numDimMom >= 2)
			{
				if(dimSizesMom[1] == NewRadStructAccessData.ne) needUpdateMomX = false;
			}
		}
		DisposeHandle(textH);
	}

	if((NewRadStructAccessData.wMomZ == 0) || (NewRadStructAccessData.pMomZ == 0)) needUpdateMomZ = false;
	//do update only if moments wave already exist
	else
	{
		textH = NewHandle(0L);
		RadIndices[0] = 16;
		if(result = MDGetTextWavePointValue(NewRadStructAccessData.wRad, RadIndices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
		waveHndl wRadMomZ = FetchWave(*textH);
		if(wRadMomZ != NIL)
		{
			if(result = MDGetWaveDimensions(wRadMomZ, &numDimMom, dimSizesMom)) return result;
			if(numDimMom >= 2)
			{
				if(dimSizesMom[1] == NewRadStructAccessData.ne) needUpdateMomZ = false;
			}
		}
		DisposeHandle(textH);
	}

	if(needUpdateMomX)
	{
		HSetState((Handle)(NewRadStructAccessData.wMomX), NewRadStructAccessData.hStateMomX);
		WaveName(NewRadStructAccessData.wMomX, AuxRadWaveName);
		result = KillWave(NewRadStructAccessData.wMomX);
		dimSizesMom[1] = NewRadStructAccessData.ne;
		if(result = MDMakeWave(&(NewRadStructAccessData.wMomX), AuxRadWaveName, NIL, dimSizesMom, NT_FP32, 1)) return FAILED_TO_CREATE_WAVE;
		if(result = MDAccessNumericWaveData(NewRadStructAccessData.wMomX, kMDWaveAccessMode0, &dataOffset)) return result;
		NewRadStructAccessData.hStateMomX = 0; //MoveLockHandle(NewRadStructAccessData.wMomX);
		//NewRadStructAccessData.pMomX = (float*)((char*)(*(NewRadStructAccessData.wMomX)) + dataOffset);
		NewRadStructAccessData.pMomX = (double*)((char*)(*(NewRadStructAccessData.wMomX)) + dataOffset); //OC130311
	}
	if(needUpdateMomZ)
	{
		HSetState((Handle)(NewRadStructAccessData.wMomZ), NewRadStructAccessData.hStateMomZ);
		WaveName(NewRadStructAccessData.wMomZ, AuxRadWaveName);
		result = KillWave(NewRadStructAccessData.wMomZ);
		dimSizesMom[1] = NewRadStructAccessData.ne;
		if(result = MDMakeWave(&(NewRadStructAccessData.wMomZ), AuxRadWaveName, NIL, dimSizesMom, NT_FP32, 1)) return FAILED_TO_CREATE_WAVE;
		if(result = MDAccessNumericWaveData(NewRadStructAccessData.wMomZ, kMDWaveAccessMode0, &dataOffset)) return result;
		NewRadStructAccessData.hStateMomZ = 0; //MoveLockHandle(NewRadStructAccessData.wMomZ);
		//NewRadStructAccessData.pMomZ = (float*)((char*)(*(NewRadStructAccessData.wMomZ)) + dataOffset);
		NewRadStructAccessData.pMomZ = (double*)((char*)(*(NewRadStructAccessData.wMomZ)) + dataOffset); //OC130311
	}
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::ModifyStokesNeNxNz(srTStokesStructAccessData& Stokes)
{
#ifdef __IGOR_PRO__

	long dimensionSizes[MAX_DIMENSIONS+1];
	int i;
	for(i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
	dimensionSizes[0] = 4;
	dimensionSizes[1] = Stokes.ne;
	dimensionSizes[2] = Stokes.nx;
	dimensionSizes[3] = Stokes.nz;

	DOUBLE StepE = Stokes.eStep, StepX = Stokes.xStep, StepZ = Stokes.zStep;
	DOUBLE StartE = Stokes.eStart, StartX = Stokes.xStart, StartZ = Stokes.zStart;

	//char AuxRadWaveName[MAX_OBJ_NAME+1];
	int result;
	long dataOffset;

	HSetState((Handle)(Stokes.wSto), Stokes.hStateSto);

	if(result = MDChangeWave(Stokes.wSto, -1, dimensionSizes)) return result;
	if(result = MDSetWaveScaling(Stokes.wSto, 1, &StepE, &StartE)) return result;
	if(result = MDSetWaveScaling(Stokes.wSto, 2, &StepX, &StartX)) return result;
	if(result = MDSetWaveScaling(Stokes.wSto, 3, &StepZ, &StartZ)) return result;

	if(result = MDAccessNumericWaveData(Stokes.wSto, kMDWaveAccessMode0, &dataOffset)) return result;
	Stokes.hStateSto = 0; //MoveLockHandle(Stokes.wSto);

	char* dataStartPtr = (char*)(*(Stokes.wSto)) + dataOffset;
	Stokes.pBaseSto = (float*)dataStartPtr;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetRadStructNames(srTSRWRadStructAccessData& RadAccessData, srTSRWRadStructWaveNames& Names)
{
#ifdef __IGOR_PRO__

	WaveName(RadAccessData.wRad, Names.NameRad);
	if(*(Names.NameRad) == '\0') return IMPROPER_RADIATION_STRUCTURE;

	WaveName(RadAccessData.wRadX, Names.NameRadX);
	if(*(Names.NameRadX) == '\0') return IMPROPER_RADIATION_STRUCTURE;
	WaveName(RadAccessData.wRadZ, Names.NameRadZ);
	if(*(Names.NameRadZ) == '\0') return IMPROPER_RADIATION_STRUCTURE;

	if(RadAccessData.wElecBeam != NIL)
	{
		WaveName(RadAccessData.wElecBeam, Names.NameElecBeam);
		if(*(Names.NameElecBeam) == '\0') return IMPROPER_RADIATION_STRUCTURE;
	}
	if(RadAccessData.wTrj != NIL)
	{
		WaveName(RadAccessData.wTrj, Names.NameTrj);
		if(*(Names.NameTrj) == '\0') return IMPROPER_RADIATION_STRUCTURE;
	}

	WaveName(RadAccessData.w4x4PropMatr, Names.Name4x4PropMatr);
	if(*(Names.Name4x4PropMatr) == '\0') return IMPROPER_RADIATION_STRUCTURE;

	WaveName(RadAccessData.wMomX, Names.NameMomX);
	if(*(Names.NameMomX) == '\0') return IMPROPER_RADIATION_STRUCTURE;
	WaveName(RadAccessData.wMomZ, Names.NameMomZ);
	if(*(Names.NameMomZ) == '\0') return IMPROPER_RADIATION_STRUCTURE;

	WaveName(RadAccessData.wWfrAuxData, Names.NameWfrAuxData);
	if(*(Names.NameWfrAuxData) == '\0') return IMPROPER_RADIATION_STRUCTURE;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::RenameRadStruct(srTSRWRadStructAccessData& RadAccessData, srTSRWRadStructWaveNames& Names)
{// Renames text wave and all dependent numerical waves
#ifdef __IGOR_PRO__

	int result;
	char AuxWaveName[MAX_OBJ_NAME+1];

	if(RadAccessData.wRad == NIL) return IMPROPER_RADIATION_STRUCTURE;

	DataFolderHandle dataFolderH;
	if(result = GetWavesDataFolder(RadAccessData.wRad, &dataFolderH)) return result;

//Ex
	if((RadAccessData.wRadX != NIL) && (*(Names.NameRadX) != '\0'))
	{
		WaveName(RadAccessData.wRadX, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameRadX))
		{
			HSetState((Handle)(RadAccessData.wRadX), RadAccessData.hStateRadX);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameRadX)) return result;
			WaveHandleModified(RadAccessData.wRadX);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 0, Names.NameRadX)) return result;
		}
	}
//Ez
	if((RadAccessData.wRadZ != NIL) && (*(Names.NameRadZ) != '\0'))
	{
		WaveName(RadAccessData.wRadZ, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameRadZ))
		{
			HSetState((Handle)(RadAccessData.wRadZ), RadAccessData.hStateRadZ);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameRadZ)) return result;
			WaveHandleModified(RadAccessData.wRadZ);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 1, Names.NameRadZ)) return result;
		}
	}

//Electron Beam or Trajectory
	if((RadAccessData.wElecBeam != NIL) && (*(Names.NameElecBeam) != '\0'))
	{
		WaveName(RadAccessData.wElecBeam, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameElecBeam))
		{
			HSetState((Handle)(RadAccessData.wElecBeam), RadAccessData.hStateElecBeam);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameElecBeam)) return result;
			WaveHandleModified(RadAccessData.wElecBeam);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 12, Names.NameElecBeam)) return result;
		}
	}
	if((RadAccessData.wTrj != NIL) && (*(Names.NameTrj) != '\0'))
	{
		WaveName(RadAccessData.wTrj, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameTrj))
		{
			HSetState((Handle)(RadAccessData.wTrj), RadAccessData.hStateTrj);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameTrj)) return result;
			WaveHandleModified(RadAccessData.wTrj);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 12, Names.NameTrj)) return result;
		}
	}

//4x4PropMatr
	if((RadAccessData.w4x4PropMatr != NIL) && (*(Names.Name4x4PropMatr) != '\0'))
	{
		WaveName(RadAccessData.w4x4PropMatr, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.Name4x4PropMatr))
		{
			HSetState((Handle)(RadAccessData.w4x4PropMatr), RadAccessData.hState4x4PropMatr);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.Name4x4PropMatr)) return result;
			WaveHandleModified(RadAccessData.w4x4PropMatr);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 13, Names.Name4x4PropMatr)) return result;
		}
	}

//MomX
	if((RadAccessData.wMomX != NIL) && (*(Names.NameMomX) != '\0'))
	{
		WaveName(RadAccessData.wMomX, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameMomX))
		{
			HSetState((Handle)(RadAccessData.wMomX), RadAccessData.hStateMomX);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameMomX)) return result;
			WaveHandleModified(RadAccessData.wMomX);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 15, Names.NameMomX)) return result;
		}
	}
//MomZ
	if((RadAccessData.wMomZ != NIL) && (*(Names.NameMomZ) != '\0'))
	{
		WaveName(RadAccessData.wMomZ, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameMomZ))
		{
			HSetState((Handle)(RadAccessData.wMomZ), RadAccessData.hStateMomZ);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameMomZ)) return result;
			WaveHandleModified(RadAccessData.wMomZ);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 16, Names.NameMomZ)) return result;
		}
	}

//WfrAuxData
	if((RadAccessData.wWfrAuxData != NIL) && (*(Names.NameWfrAuxData) != '\0'))
	{
		WaveName(RadAccessData.wWfrAuxData, AuxWaveName);
		if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
		if(strcmp(AuxWaveName, Names.NameWfrAuxData))
		{
			HSetState((Handle)(RadAccessData.wWfrAuxData), RadAccessData.hStateWfrAuxData);
			if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameWfrAuxData)) return result;
			WaveHandleModified(RadAccessData.wWfrAuxData);

			if(result = UpdateTextPositionInSRWRad(&RadAccessData, 18, Names.NameWfrAuxData)) return result;
		}
	}

//Rad
	WaveName(RadAccessData.wRad, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;
	if(strcmp(AuxWaveName, Names.NameRad))
	{
		if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, Names.NameRad)) return result;
		WaveHandleModified(RadAccessData.wRad);
	}

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::RenameRadStruct(srTSRWRadStructAccessData& RadAccessData, char* NewName)
{// Renames text wave and all dependent numerical waves
#ifdef __IGOR_PRO__

	int result;
	char AuxWaveName[MAX_OBJ_NAME+1], AuxNewWaveName[MAX_OBJ_NAME+1];
	char AuxRadExtens[6], AuxRadFieldExtens[6];
	const int RadExtensLen = 4;
	const int RadFieldExtensLen = 5;
	int i;

	DataFolderHandle dataFolderH;
	if(result = GetWavesDataFolder(RadAccessData.wRad, &dataFolderH)) return result;

//Text Wave
	WaveName(RadAccessData.wRad, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	int OldNameLen = strlen(AuxWaveName);
	char *tAuxWaveName = AuxWaveName + (OldNameLen - RadExtensLen);
	for(i=0; i<RadExtensLen; i++) AuxRadExtens[i] = *(tAuxWaveName++);
	AuxRadExtens[RadExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//Num. Wave Ex
	WaveName(RadAccessData.wRadX, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadFieldExtensLen);
	for(i=0; i<RadFieldExtensLen; i++) AuxRadFieldExtens[i] = *(tAuxWaveName++);
	AuxRadFieldExtens[RadFieldExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadFieldExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//Num. Wave Ez
	WaveName(RadAccessData.wRadZ, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadFieldExtensLen);
	for(i=0; i<RadFieldExtensLen; i++) AuxRadFieldExtens[i] = *(tAuxWaveName++);
	AuxRadFieldExtens[RadFieldExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadFieldExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//Electron Beam or Trajectory
	waveHndl wElecOrTraj = (RadAccessData.wElecBeam != NIL)? RadAccessData.wElecBeam : RadAccessData.wTrj;
	WaveName(wElecOrTraj, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadExtensLen);
	for(i=0; i<RadExtensLen; i++) AuxRadExtens[i] = *(tAuxWaveName++);
	AuxRadExtens[RadExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//4x4PropMatr
	WaveName(RadAccessData.w4x4PropMatr, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadExtensLen);
	for(i=0; i<RadExtensLen; i++) AuxRadExtens[i] = *(tAuxWaveName++);
	AuxRadExtens[RadExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//MomX
	WaveName(RadAccessData.wMomX, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadFieldExtensLen);
	for(i=0; i<RadFieldExtensLen; i++) AuxRadFieldExtens[i] = *(tAuxWaveName++);
	AuxRadFieldExtens[RadFieldExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadFieldExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//MomZ
	WaveName(RadAccessData.wMomZ, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadFieldExtensLen);
	for(i=0; i<RadFieldExtensLen; i++) AuxRadFieldExtens[i] = *(tAuxWaveName++);
	AuxRadFieldExtens[RadFieldExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadFieldExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

//WfrAuxData
	WaveName(RadAccessData.wWfrAuxData, AuxWaveName);
	if(*AuxWaveName == '\0') return IMPROPER_RADIATION_STRUCTURE;

	OldNameLen = strlen(AuxWaveName);
	tAuxWaveName = AuxWaveName + (OldNameLen - RadExtensLen);
	for(i=0; i<RadExtensLen; i++) AuxRadExtens[i] = *(tAuxWaveName++);
	AuxRadExtens[RadExtensLen] = '\0';

	strcpy(AuxNewWaveName, NewName);
	strcat(AuxNewWaveName, AuxRadExtens);
	if(result = RenameDataFolderObject(dataFolderH, WAVE_OBJECT, AuxWaveName, AuxNewWaveName)) return result;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetStokesStructAccessData(srTStokesStructAccessData* pStokesAccessData)
{
#ifdef __IGOR_PRO__

	waveHndl wavH;
	if(pIgorPerStokesInputStruct != 0) wavH = pIgorPerStokesInputStruct->wStokes;
	else if(pIgorStokesWigInputStruct != 0) wavH = pIgorStokesWigInputStruct->wStokes;
	else if(pIgorStokesConstInputStruct != 0) wavH = pIgorStokesConstInputStruct->wStokes;
	else if(pIgorRadPropagStokesMultiElecInputStruct != 0) wavH = pIgorRadPropagStokesMultiElecInputStruct->wStokes;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP32) return IMPROPER_STOKES_STRUCTURE;

	int result;
  	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	if(numDimensions != 4) return IMPROPER_STOKES_STRUCTURE;
	if(dimensionSizes[0] != 4) return IMPROPER_STOKES_STRUCTURE;

	pStokesAccessData->wSto = wavH;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	pStokesAccessData->hStateSto = 0; //MoveLockHandle(wavH);

	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	pStokesAccessData->pBaseSto = (float*)dataStartPtr;

	pStokesAccessData->ne = dimensionSizes[1];
	pStokesAccessData->nx = dimensionSizes[2];
	pStokesAccessData->nz = dimensionSizes[3];

	DOUBLE eStep, eStart, xStep, xStart, zStep, zStart;
	if(result = MDGetWaveScaling(wavH, 1, &eStep, &eStart)) return result;
	if(result = MDGetWaveScaling(wavH, 2, &xStep, &xStart)) return result;
	if(result = MDGetWaveScaling(wavH, 3, &zStep, &zStart)) return result;
	pStokesAccessData->eStep = eStep; pStokesAccessData->eStart = eStart;
	pStokesAccessData->xStep = xStep; pStokesAccessData->xStart = xStart;
	pStokesAccessData->zStep = zStep; pStokesAccessData->zStart = zStart;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetPowDensStructAccessData(srTPowDensStructAccessData* pPowDensStructAccessData)
{
#ifdef __IGOR_PRO__

	waveHndl wavH;
	wavH = pIgorPowDensInputStruct->wPowDens;
	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP32) return IMPROPER_POWER_DENSITY_STRUCTURE;

	int result;
  	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	if(numDimensions != 2) return IMPROPER_POWER_DENSITY_STRUCTURE;

	pPowDensStructAccessData->wPowDens = wavH;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	pPowDensStructAccessData->hStatePowDens = 0; //MoveLockHandle(wavH);

	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	pPowDensStructAccessData->pBasePowDens = (float*)dataStartPtr;

	pPowDensStructAccessData->nx = dimensionSizes[0];
	pPowDensStructAccessData->nz = dimensionSizes[1];

	DOUBLE xStep, xStart, zStep, zStart;
	if(result = MDGetWaveScaling(wavH, 0, &xStep, &xStart)) return result;
	if(result = MDGetWaveScaling(wavH, 1, &zStep, &zStart)) return result;

	pPowDensStructAccessData->xStep = xStep; pPowDensStructAccessData->xStart = xStart;
	pPowDensStructAccessData->zStep = zStep; pPowDensStructAccessData->zStart = zStart;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::ChangeObservationData(srTHandleOfOneWaveStruct* pHandleOfOneWaveStruct, srTWfrSmp& DistrInfoDat)
{
	waveHndl wavH = pHandleOfOneWaveStruct->WaveHndl;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	int result;

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];
	if(numRows < 13) return BAD_OBS_WAVE_FORMAT;

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;

	dp += 4; // The fields 0 - 3 are not used
	*(dp++) = DistrInfoDat.yStart; // Output in m

	char LocTreatLambdaAsEnergyIn_eV = DistrInfoDat.TreatLambdaAsEnergyIn_eV || (DistrInfoDat.PhotonEnergyWavelengthUnits == 1);
	if(LocTreatLambdaAsEnergyIn_eV)
	{
		*(dp++) = DistrInfoDat.LambStart; // Output in eV;
		*(dp++) = DistrInfoDat.LambEnd; // Output in eV;
		*(dp++) = DistrInfoDat.nLamb;
	}

	if(DistrInfoDat.CoordOrAngPresentation == CoordPres)
	{
		*(dp++) = DistrInfoDat.xStart; // Output in m;
		*(dp++) = DistrInfoDat.xEnd; // Output in m;
		*(dp++) = DistrInfoDat.nx;

		*(dp++) = DistrInfoDat.zStart; // Output in m;
		*(dp++) = DistrInfoDat.zEnd; // Output in m;
		*dp = DistrInfoDat.nz;
	}

	HSetState((Handle)wavH, hState);
	WaveHandleModified(wavH);

	return 0;

//#else

	//return 0;
}
#endif

//*************************************************************************

int srTSend::GetAuxObsTreatParamFormat1(srTWfrSmp& DistrInfoDat)
{
#if defined(__IGOR_PRO__)

	waveHndl wavH;
	if(pIgorRadInputStruct != 0) wavH = pIgorRadInputStruct->wAuxParam;
	else if(pIgorIsotrSrcInputStruct != 0) wavH = pIgorIsotrSrcInputStruct->wAuxParam;
	else if(pIgorGsnBeamInputStruct != 0) wavH = pIgorGsnBeamInputStruct->wAuxParam;
	else if(pIgorWfrFromTrjInputStruct != 0) wavH = pIgorWfrFromTrjInputStruct->wAuxParam;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	//long numRows = dimensionSizes[0];

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp = (DOUBLE*)dataStartPtr;

	char LocAllowAutoChoiceOfNxNzForPropagat = char(*(dp++));
	if(DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat != LocAllowAutoChoiceOfNxNzForPropagat)
	{
		DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat = LocAllowAutoChoiceOfNxNzForPropagat;
		DistrInfoDat.InputWasModified = 1;
	}
	double LocNxNzOversamplingParam = double(*dp);
	if(DistrInfoDat.NxNzOversamplingParam != LocNxNzOversamplingParam)
	{
		DistrInfoDat.NxNzOversamplingParam = LocNxNzOversamplingParam;
		DistrInfoDat.InputWasModified = 1;
	}

	HSetState((Handle)wavH, hState);

	return 0;

#else

	return 0;

#endif
}


//*************************************************************************

int srTSend::GetSASEInputRadDataFormat1(srTRadSASE& InRadSASE, srTSRWRadStructAccessData& SeedRadSASE)
{
#ifdef __IGOR_PRO__

	waveHndl wRadSASE;
	if(pIgorWfrSASEInputStruct != 0) wRadSASE = pIgorWfrSASEInputStruct->wRadSASE;
	//else if(pIgorTrjInputStruct != 0) wField = pIgorTrjInputStruct->wField;

	if(wRadSASE == NIL) return NOWAV;
	int waveType = WaveType(wRadSASE);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;
	int result = 0;

	long numDimensions, dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wRadSASE, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];

	long RadIndices[MAX_DIMENSIONS];
	RadIndices[0] = 0;
	Handle textH = NewHandle(0L);
	if(result = MDGetTextWavePointValue(wRadSASE, RadIndices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;

	bool isWfr = (numRows >= 20);
	if(textH == 0) isWfr = false;
	else if(*textH == 0) isWfr = false;
	else if(strlen(*textH) <= 0) isWfr = false;
	if(isWfr) 
	{
		InRadSASE.Power = 0;
		
		waveHndl wWfrPrev = HandleOfSRWRadStruct.wRad;
		HandleOfSRWRadStruct.wRad = wRadSASE;
		if(result = GetSRWRadStructAccessData(&SeedRadSASE)) return result;
		HandleOfSRWRadStruct.wRad = wWfrPrev;
		if((SeedRadSASE.PresT != 1) || (SeedRadSASE.ne <= 1)) return TIME_DOMAIN_RAD_STRUCT_REQUIRED;
		else return 0;
	}
	DisposeHandle(textH);

	if(result = GetDoubleFromTextWave1D(wRadSASE, 2, InRadSASE.Power)) return result;
	if(result = GetDoubleFromTextWave1D(wRadSASE, 3, InRadSASE.WaistDiam)) return result;
	if(result = GetDoubleFromTextWave1D(wRadSASE, 4, InRadSASE.WaistLongPos)) return result;
	//if(result = GetDoubleFromTextWave1D(wRadSASE, 5, InRadSASE.PhotonEnergySim)) return result;

	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::GetSASEPrecDataFormat1(srTPrecSASE& PrecSASE)
{
#ifdef __IGOR_PRO__

	waveHndl wavH;
	if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wPrecSASE;
	//else if(pIgorTrjInputStruct != 0) wField = pIgorTrjInputStruct->wField;

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != NT_FP64) return NT_FP64_WAVE_REQUIRED;

	int result;

	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];

	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	//long numRows = dimensionSizes[0];

	long dataOffset;
	if(result = MDAccessNumericWaveData(wavH, kMDWaveAccessMode0, &dataOffset)) return result;
	int hState = 0; //MoveLockHandle(wavH);
	char* dataStartPtr = (char*)(*wavH) + dataOffset;
	DOUBLE* dp0 = (DOUBLE*)dataStartPtr;

	PrecSASE.npart = long(*(dp0));
	PrecSASE.rmax0 = double(*(dp0+1));
	PrecSASE.ncar = long(*(dp0+2));
	PrecSASE.nptr = long(*(dp0+3));
	PrecSASE.nscr = long(*(dp0+4));
	PrecSASE.nscz = long(*(dp0+5));
	PrecSASE.lbc = char(*(dp0+6)) - 1;
	PrecSASE.delz = double(*(dp0+7));
	PrecSASE.zstop = double(*(dp0+8));
	PrecSASE.iorb = char(*(dp0+9)) - 1;

	char TestVal = char(*(dp0+10));
	if(TestVal > 0) PrecSASE.itdp = TestVal - 1;
	else PrecSASE.itdp = 0;

	PrecSASE.nslice = long(*(dp0+11));
	PrecSASE.zsep = double(*(dp0+12));
	PrecSASE.ntail = long(*(dp0+13));
	PrecSASE.UseElecDistr = char(*(dp0+14)); // - 1; //1 subtracted in Igor
	PrecSASE.CreateElecDistr = char(*(dp0+15));
	//field +15 is reserved for "CreateElecDist" param

	PrecSASE.photEn_xlamds = double(*(dp0+16)); //photon energy in [eV] to calc. xlamds
	PrecSASE.alignradf = int(*(dp0+17)) - 1; // "ALIGNRADF" <>0 imported radfile is aligned to electron beam
	PrecSASE.offsetradf = int(*(dp0+18));

	PrecSASE.AllowAutoChoiceOfNxNzForPropagat = char(*(dp0+20)) - 1;
	PrecSASE.NxNzOversamplingParam = double(*(dp0+21));

	//PrecSASE.iPower = int(*(dp0+15)) - 1; // ensure 2:yes, 1:no in Igor
	//PrecSASE.iRadHorSize = int(*(dp0+16)) - 1;
	//PrecSASE.iRadVertSize = int(*(dp0+17)) - 1;

	HSetState((Handle)wavH, hState);
	return 0;

#else

	return 0;

#endif
}

//*************************************************************************

int srTSend::SetupSASEControlStruct(srTControlAccessSASE& ControlAccessSASE, int numHarm)
{
#ifdef __IGOR_PRO__

	waveHndl wavH;
	if(pIgorWfrSASEInputStruct != 0) wavH = pIgorWfrSASEInputStruct->wControlSASE;
	//else if(pIgorTrjInputStruct != 0) wavH = pIgorTrjInputStruct->wField;

	if(wavH == NIL) return 0; // normal exit here

	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	int result = 0;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	long numRows = dimensionSizes[0];

	if(numRows == 0) return 0; // normal exit here

	DOUBLE sStep = (DOUBLE)ControlAccessSASE.sStep;
	DOUBLE sStart = (DOUBLE)ControlAccessSASE.sStart;
	//char* pCharAux = 0;

	ControlAccessSASE.wControl = wavH;

	char* pAuxChar = 0;
	long Indices[MAX_DIMENSIONS];
	Indices[1] = 0;

	long NpAr[] = {ControlAccessSASE.nt, ControlAccessSASE.ns};
	DOUBLE tStep = (DOUBLE)ControlAccessSASE.tStep;
	DOUBLE tStart = (DOUBLE)ControlAccessSASE.tStart;
	DOUBLE StartAr[] = {tStart, sStart};
	DOUBLE StepAr[] = {tStep, sStep};

	for(int iHarm=0; iHarm<numHarm; iHarm++)
	{
		Indices[0] = iHarm; 
		if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.ar_wPower_vs_s[iHarm], pAuxChar, ControlAccessSASE.ar_hStatePower_vs_s[iHarm])) return result;
		ControlAccessSASE.ar_pBasePower_vs_s[iHarm] = (float*)pAuxChar;

		Indices[0] = 7 + iHarm;
		if(result = SetupControlStructMD(wavH, Indices, NpAr, StartAr, StepAr, 2, ControlAccessSASE.ar_wPower_vs_t[iHarm], pAuxChar, ControlAccessSASE.ar_hStatePower_vs_t[iHarm])) return result;
		ControlAccessSASE.ar_pBasePower_vs_t[iHarm] = (float*)pAuxChar;

		Indices[0] = 14 + iHarm;
		if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.ar_wEnergy_vs_s[iHarm], pAuxChar, ControlAccessSASE.ar_hStateEnergy_vs_s[iHarm])) return result;
		ControlAccessSASE.ar_pBaseEnergy_vs_s[iHarm] = (float*)pAuxChar;

		Indices[0] = 21 + iHarm;
		if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.ar_wPeakPower_vs_s[iHarm], pAuxChar, ControlAccessSASE.ar_hStatePeakPower_vs_s[iHarm])) return result;
		ControlAccessSASE.ar_pBasePeakPower_vs_s[iHarm] = (float*)pAuxChar;

		Indices[0] = 32 + iHarm;
		if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.ar_wBunchFact_vs_s[iHarm], pAuxChar, ControlAccessSASE.ar_hStateBunchFact_vs_s[iHarm])) return result;
		ControlAccessSASE.ar_pBaseBunchFact_vs_s[iHarm] = (float*)pAuxChar;
	}

	Indices[0] = 28; 
	if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wRadPhase_vs_s, pAuxChar, ControlAccessSASE.hStateRadPhase_vs_s)) return result;
	ControlAccessSASE.pBaseRadPhase_vs_s = (float*)pAuxChar;

	Indices[0] = 29; 
	if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wRadSize_vs_s, pAuxChar, ControlAccessSASE.hStateRadSize_vs_s)) return result;
	ControlAccessSASE.pBaseRadSize_vs_s = (float*)pAuxChar;

	Indices[0] = 30; 
	if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBmSizeX_vs_s, pAuxChar, ControlAccessSASE.hStateBmSizeX_vs_s)) return result;
	ControlAccessSASE.pBaseBmSizeX_vs_s = (float*)pAuxChar;

	Indices[0] = 31; 
	if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBmSizeZ_vs_s, pAuxChar, ControlAccessSASE.hStateBmSizeZ_vs_s)) return result;
	ControlAccessSASE.pBaseBmSizeZ_vs_s = (float*)pAuxChar;

	//Indices[0] = 5; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBunchFact1_vs_s, pAuxChar, ControlAccessSASE.hStateBunchFact1_vs_s)) return result;
	//ControlAccessSASE.pBaseBunchFact1_vs_s = (float*)pAuxChar;
	//Indices[0] = 6; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBunchFact2_vs_s, pAuxChar, ControlAccessSASE.hStateBunchFact2_vs_s)) return result;
	//ControlAccessSASE.pBaseBunchFact2_vs_s = (float*)pAuxChar;
	//Indices[0] = 7; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBunchFact3_vs_s, pAuxChar, ControlAccessSASE.hStateBunchFact3_vs_s)) return result;
	//ControlAccessSASE.pBaseBunchFact3_vs_s = (float*)pAuxChar;
	//Indices[0] = 8; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBunchFact4_vs_s, pAuxChar, ControlAccessSASE.hStateBunchFact4_vs_s)) return result;
	//ControlAccessSASE.pBaseBunchFact4_vs_s = (float*)pAuxChar;
	//Indices[0] = 9; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wBunchFact5_vs_s, pAuxChar, ControlAccessSASE.hStateBunchFact5_vs_s)) return result;
	//ControlAccessSASE.pBaseBunchFact5_vs_s = (float*)pAuxChar;

	//Indices[0] = 10; 
	//long NpAr[] = {ControlAccessSASE.nt, ControlAccessSASE.ns};
	//DOUBLE StartAr[] = {tStart, sStart};
	//DOUBLE StepAr[] = {tStep, sStep};
	//if(result = SetupControlStructMD(wavH, Indices, NpAr, StartAr, StepAr, 2, ControlAccessSASE.wPower_vs_t, pAuxChar, ControlAccessSASE.hStatePower_vs_t)) return result;
	//ControlAccessSASE.pBasePower_vs_t = (float*)pAuxChar;

	//Indices[0] = 11; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wEnergy_vs_s, pAuxChar, ControlAccessSASE.hStateEnergy_vs_s)) return result;
	//ControlAccessSASE.pBaseEnergy_vs_s = (float*)pAuxChar;

	//Indices[0] = 12; 
	//if(result = SetupControlStruct1D(wavH, Indices, ControlAccessSASE.ns, sStart, sStep, ControlAccessSASE.wPeakPower_vs_s, pAuxChar, ControlAccessSASE.hStatePeakPower_vs_s)) return result;
	//ControlAccessSASE.pBasePeakPower_vs_s = (float*)pAuxChar;

	char NameCntrlWave[MAX_OBJ_NAME+1];
	WaveName(wavH, NameCntrlWave);
	char strToExe[MAX_OBJ_NAME+60];
	strcpy(strToExe, "SrwSASECntrlDisplay(\"");
	strcat(strToExe, NameCntrlWave);
	strcat(strToExe, "\")");

	if(result = XOPSilentCommand(strToExe)) return result;

/**
//old version
	char NameOfWave[MAX_OBJ_NAME+1];
	WaveName(wavH, NameOfWave);
	char* pEnding = strrchr(NameOfWave, '_');
	*pEnding = '\0';
	if(strlen(NameOfWave) > 28) NameOfWave[27] = '\0';

	int result;
	long dataOffset, numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	long Indices[MAX_DIMENSIONS];

	DOUBLE sStep = (DOUBLE)ControlAccessSASE.sStep;
	DOUBLE sStart = (DOUBLE)ControlAccessSASE.sStart;
	char* pCharAux = 0;

	if(PrecSASE.iPower)
	{
		char LocWaveName[32] = "";
		strcat(LocWaveName, NameOfWave);
		strcat(LocWaveName, "0");
		strcat(LocWaveName, ControlWaveNameExtension);

		dimensionSizes[0] = ControlAccessSASE.ns;
		dimensionSizes[1] = dimensionSizes[2] = dimensionSizes[3] = 0;
		if(result = MDMakeWave(&ControlAccessSASE.wPower_vs_s, LocWaveName, NIL, dimensionSizes, NT_FP32, 1)) return result;
		if(result = MDSetWaveScaling(ControlAccessSASE.wPower_vs_s, ROWS, &sStep, &sStart)) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wPower_vs_s, ROWS, "m")) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wPower_vs_s, -1, "W")) return result;

		if(result = GetDataPointerAndStateAndLockNumericWave(ControlAccessSASE.wPower_vs_s, pCharAux, ControlAccessSASE.hStatePower_vs_s)) return result;
		ControlAccessSASE.pBasePower_vs_s = (float*)pCharAux;

		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 0, LocWaveName)) return result;
	}
	else
	{
		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 0, "")) return result;
	}

	if(PrecSASE.iRadHorSize)
	{
		char LocWaveName[32] = "";
		strcat(LocWaveName, NameOfWave);
		strcat(LocWaveName, "1");
		strcat(LocWaveName, ControlWaveNameExtension);

		dimensionSizes[0] = ControlAccessSASE.ns;
		dimensionSizes[1] = dimensionSizes[2] = dimensionSizes[3] = 0;
		if(result = MDMakeWave(&ControlAccessSASE.wRadSizeX_vs_s, LocWaveName, NIL, dimensionSizes, NT_FP32, 1)) return result;
		if(result = MDSetWaveScaling(ControlAccessSASE.wRadSizeX_vs_s, ROWS, &sStep, &sStart)) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wRadSizeX_vs_s, ROWS, "m")) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wRadSizeX_vs_s, -1, "m")) return result;

		if(result = GetDataPointerAndStateAndLockNumericWave(ControlAccessSASE.wRadSizeX_vs_s, pCharAux, ControlAccessSASE.hStateRadSizeX_vs_s)) return result;
		ControlAccessSASE.pBaseRadSizeX_vs_s = (float*)pCharAux;

		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 1, LocWaveName)) return result;
	}
	else
	{
		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 1, "")) return result;
	}

	if(PrecSASE.iRadVertSize)
	{
		char LocWaveName[32] = "";
		strcat(LocWaveName, NameOfWave);
		strcat(LocWaveName, "2");
		strcat(LocWaveName, ControlWaveNameExtension);

		dimensionSizes[0] = ControlAccessSASE.ns;
		dimensionSizes[1] = dimensionSizes[2] = dimensionSizes[3] = 0;

		if(result = MDMakeWave(&ControlAccessSASE.wRadSizeZ_vs_s, LocWaveName, NIL, dimensionSizes, NT_FP32, 1)) return result;
		if(result = MDSetWaveScaling(ControlAccessSASE.wRadSizeZ_vs_s, ROWS, &sStep, &sStart)) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wRadSizeZ_vs_s, ROWS, "m")) return result;
		if(result = MDSetWaveUnits(ControlAccessSASE.wRadSizeZ_vs_s, -1, "m")) return result;

		if(result = GetDataPointerAndStateAndLockNumericWave(ControlAccessSASE.wRadSizeZ_vs_s, pCharAux, ControlAccessSASE.hStateRadSizeZ_vs_s)) return result;
		ControlAccessSASE.pBaseRadSizeZ_vs_s = (float*)pCharAux;

		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 2, LocWaveName)) return result;
	}
	else
	{
		if(result = UpdateTextWave1DPointValue(ControlAccessSASE.wControl, 2, "")) return result;
	}
**/
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__

int srTSend::SetupControlStruct1D(waveHndl wText, long* Indices, long NewNp, DOUBLE NewStart, DOUBLE NewStep, waveHndl& wHndl, char*& pBase, int& hState)
//int srTSend::SetupControlStruct1D(waveHndl wText, long* Indices, long NewNp, double NewStart, double NewStep, waveHndl& wHndl, char*& pBase, int& hState)
{//Resizes, sets new start and step and returns parameters of a numeric wave with the name specified in a text wave
 //there is no validation here
	int result = 0;
	wHndl = NIL; pBase = 0; hState = 0;

	Handle textH = NewHandle(0L);
	if(result = MDGetTextWavePointValue(wText, Indices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH != '\0')
			{
				wHndl = FetchWave(*textH);
				if(wHndl != NIL)
				{
					long dimensionSizes[MAX_DIMENSIONS+1];
					dimensionSizes[0] = NewNp;
					//dimensionSizes[1] = dimensionSizes[2] = dimensionSizes[3] = 0;
					for(int i=1; i<MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;

					if(result = MDChangeWave(wHndl, -1, dimensionSizes)) return result;
					if(result = MDSetWaveScaling(wHndl, ROWS, &NewStep, &NewStart)) return result;

					if(result = GetDataPointerAndStateAndLockNumericWave(wHndl, pBase, hState)) return result;
				}
			}
		}
	}
	DisposeHandle(textH);
	return 0;
}

#endif

//*************************************************************************

#ifdef __IGOR_PRO__

int srTSend::SetupControlStructMD(waveHndl wText, long* Indices, long* NewNpAr, DOUBLE* NewStartAr, DOUBLE* NewStepAr, int numDim, waveHndl& wHndl, char*& pBase, int& hState)
{//Resizes, sets new start and step and returns parameters of a numeric wave with the name specified in a text wave
 // there is no validation here
	int result = 0;
	wHndl = NIL; pBase = 0; hState = 0;

	Handle textH = NewHandle(0L);
	if(result = MDGetTextWavePointValue(wText, Indices, textH)) return result;
	if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
	if(textH != 0)
	{
		if(*textH != 0)
		{
			if(**textH != '\0')
			{
				wHndl = FetchWave(*textH);
				if(wHndl != NIL)
				{
					long dimensionSizes[MAX_DIMENSIONS+1];

					if(numDim >= MAX_DIMENSIONS) numDim = MAX_DIMENSIONS;
					for(int i=0; i<numDim; i++) dimensionSizes[i] = NewNpAr[i];
					for(int i=numDim; i<MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;

					if(result = MDChangeWave(wHndl, -1, dimensionSizes)) return result;

					if(numDim >= 1) if(result = MDSetWaveScaling(wHndl, ROWS, NewStepAr, NewStartAr)) return result;
					if(numDim >= 2) if(result = MDSetWaveScaling(wHndl, COLUMNS, NewStepAr + 1, NewStartAr + 1)) return result;
					if(numDim >= 3) if(result = MDSetWaveScaling(wHndl, LAYERS, NewStepAr + 2, NewStartAr + 2)) return result;
					if(numDim >= 4) if(result = MDSetWaveScaling(wHndl, CHUNKS, NewStepAr + 3, NewStartAr + 3)) return result;

					if(result = GetDataPointerAndStateAndLockNumericWave(wHndl, pBase, hState)) return result;
				}
			}
		}
	}
	DisposeHandle(textH);
	return 0;
}

#endif

//*************************************************************************

int srTSend::FinishWorkingWithControlSASEStruct(srTControlAccessSASE& ControlAccessSASE)
{
#ifdef __IGOR_PRO__

	bool ThereAreControlStructures = false;

	for(int i=0; i<7; i++)
	{
		if((ControlAccessSASE.ar_pBasePower_vs_s[i] != 0) ||
		   (ControlAccessSASE.ar_pBasePower_vs_t[i] != 0) ||
		   (ControlAccessSASE.ar_pBaseEnergy_vs_s[i] != 0) ||
		   (ControlAccessSASE.ar_pBasePeakPower_vs_s[i] != 0) ||
		   (ControlAccessSASE.ar_pBaseBunchFact_vs_s[i] != 0))
		{
			ThereAreControlStructures = true; break;
		}
	}

	if(!ThereAreControlStructures)
	{
		if((ControlAccessSASE.pBaseRadPhase_vs_s != 0) ||
		   (ControlAccessSASE.pBaseRadSize_vs_s != 0) ||
		   (ControlAccessSASE.pBaseBmSizeX_vs_s != 0) ||
		   (ControlAccessSASE.pBaseBmSizeZ_vs_s != 0)) ThereAreControlStructures = true;
	}
	//continue here for more stuctures

	if(!ThereAreControlStructures) return 0;

	int result = 0;
	for(int i=0; i<7; i++)
	{
		if(result = ReleaseWave(ControlAccessSASE.ar_wPower_vs_s[i], ControlAccessSASE.ar_hStatePower_vs_s[i])) return result;
		if(result = ReleaseWave(ControlAccessSASE.ar_wPower_vs_t[i], ControlAccessSASE.ar_hStatePower_vs_t[i])) return result;
		if(result = ReleaseWave(ControlAccessSASE.ar_wEnergy_vs_s[i], ControlAccessSASE.ar_hStateEnergy_vs_s[i])) return result;
		if(result = ReleaseWave(ControlAccessSASE.ar_wPeakPower_vs_s[i], ControlAccessSASE.ar_hStatePeakPower_vs_s[i])) return result;
		if(result = ReleaseWave(ControlAccessSASE.ar_wBunchFact_vs_s[i], ControlAccessSASE.ar_hStateBunchFact_vs_s[i])) return result;
	}
	if(result = ReleaseWave(ControlAccessSASE.wRadPhase_vs_s, ControlAccessSASE.hStateRadPhase_vs_s)) return result;
	if(result = ReleaseWave(ControlAccessSASE.wRadSize_vs_s, ControlAccessSASE.hStateRadSize_vs_s)) return result;
	if(result = ReleaseWave(ControlAccessSASE.wBmSizeX_vs_s, ControlAccessSASE.hStateBmSizeX_vs_s)) return result;
	if(result = ReleaseWave(ControlAccessSASE.wBmSizeZ_vs_s, ControlAccessSASE.hStateBmSizeZ_vs_s)) return result;
	
	//if(result = ReleaseWave(ControlAccessSASE.wBunchFact1_vs_s, ControlAccessSASE.hStateBunchFact1_vs_s)) return result;
	//if(result = ReleaseWave(ControlAccessSASE.wBunchFact2_vs_s, ControlAccessSASE.hStateBunchFact2_vs_s)) return result;
	//if(result = ReleaseWave(ControlAccessSASE.wBunchFact3_vs_s, ControlAccessSASE.hStateBunchFact3_vs_s)) return result;
	//if(result = ReleaseWave(ControlAccessSASE.wBunchFact4_vs_s, ControlAccessSASE.hStateBunchFact4_vs_s)) return result;
	//if(result = ReleaseWave(ControlAccessSASE.wBunchFact5_vs_s, ControlAccessSASE.hStateBunchFact5_vs_s)) return result;
	//add more, if necessary
	
	WaveHandleModified(ControlAccessSASE.wControl);
	return result;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::FinishWorkingWithSRWRadStruct(srTSRWRadStructAccessData* pSRWRadStructAccessData)
{
#ifdef __IGOR_PRO__

	if((pSRWRadStructAccessData->pBaseRadX == 0) && (pSRWRadStructAccessData->pBaseRadZ == 0)) return 0;

	long RadIndices[MAX_DIMENSIONS];
	int AmOfBytes, i, result;
	char CharBuf[15];
	for(i=0; i<15; i++) CharBuf[i] = 0;

	RadIndices[0] = 2; // Presentation
	AmOfBytes = sprintf(CharBuf, "%d", char(pSRWRadStructAccessData->Pres));
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;

	DOUBLE eStep = pSRWRadStructAccessData->eStep, eStart = pSRWRadStructAccessData->eStart;
	RadIndices[0] = 3; // eStep
	AmOfBytes = sprintf(CharBuf, "%g", eStep);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 4; // eStart
	AmOfBytes = sprintf(CharBuf, "%g", eStart);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;

	DOUBLE xStep = pSRWRadStructAccessData->xStep, xStart = pSRWRadStructAccessData->xStart;
	DOUBLE zStep = pSRWRadStructAccessData->zStep, zStart = pSRWRadStructAccessData->zStart;
	//DOUBLE RobsX = pSRWRadStructAccessData->RobsX, RobsZ = pSRWRadStructAccessData->RobsZ;
	//DOUBLE RobsXAbsErr = pSRWRadStructAccessData->RobsXAbsErr, RobsZAbsErr = pSRWRadStructAccessData->RobsZAbsErr;
	//DOUBLE xc = pSRWRadStructAccessData->xc, zc = pSRWRadStructAccessData->zc;

	RadIndices[0] = 5; // xStep
	AmOfBytes = sprintf(CharBuf, "%g", xStep);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 6; // xStart
	AmOfBytes = sprintf(CharBuf, "%g", xStart);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 7; // zStep
	AmOfBytes = sprintf(CharBuf, "%g", zStep);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;
	RadIndices[0] = 8; // zStart
	AmOfBytes = sprintf(CharBuf, "%g", zStart);
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;

	RadIndices[0] = 19; // Elec. fld. units
	AmOfBytes = sprintf(CharBuf, "%d", char(pSRWRadStructAccessData->ElecFldUnit));
	if(result = UpdateNumberPositionInSRWRad(pSRWRadStructAccessData, RadIndices, CharBuf, AmOfBytes)) return result;

	char TransvUnits[MAX_UNIT_CHARS + 1];
	if(pSRWRadStructAccessData->Pres == 0) *TransvUnits = 'm';
	else *TransvUnits = 'q';
	TransvUnits[1] = '\0';

	if(pSRWRadStructAccessData->wRadX != NIL)
	{
		HSetState((Handle)(pSRWRadStructAccessData->wRadX), pSRWRadStructAccessData->hStateRadX);
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadX, 0, &eStep, &eStart)) return result;
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadX, 1, &xStep, &xStart)) return result;
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadX, 2, &zStep, &zStart)) return result;

		if(result = MDSetWaveUnits(pSRWRadStructAccessData->wRadX, COLUMNS, TransvUnits)) return result;
		if(result = MDSetWaveUnits(pSRWRadStructAccessData->wRadX, LAYERS, TransvUnits)) return result;

		WaveHandleModified(pSRWRadStructAccessData->wRadX);
	}

	if(pSRWRadStructAccessData->wRadZ != NIL)
	{
		HSetState((Handle)(pSRWRadStructAccessData->wRadZ), pSRWRadStructAccessData->hStateRadZ);
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadZ, 0, &eStep, &eStart)) return result;
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadZ, 1, &xStep, &xStart)) return result;
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wRadZ, 2, &zStep, &zStart)) return result;

		if(result = MDSetWaveUnits(pSRWRadStructAccessData->wRadZ, COLUMNS, TransvUnits)) return result;
		if(result = MDSetWaveUnits(pSRWRadStructAccessData->wRadZ, LAYERS, TransvUnits)) return result;

		WaveHandleModified(pSRWRadStructAccessData->wRadZ);
	}
// Don't work with El. field waves after this point.

// E-beam
	if(pSRWRadStructAccessData->wElecBeam != NIL)
	{
		HSetState((Handle)(pSRWRadStructAccessData->wElecBeam), pSRWRadStructAccessData->hStateElecBeam);
		WaveHandleModified(pSRWRadStructAccessData->wElecBeam);
	}

// Trajectory
	if(pSRWRadStructAccessData->wTrj != NIL)
	{
		HSetState((Handle)(pSRWRadStructAccessData->wTrj), pSRWRadStructAccessData->hStateTrj);
		WaveHandleModified(pSRWRadStructAccessData->wTrj);
	}

// 4x4 Matrix for e-beam Moments propagation
	if(pSRWRadStructAccessData->w4x4PropMatr != NIL)
	{
		HSetState((Handle)(pSRWRadStructAccessData->w4x4PropMatr), pSRWRadStructAccessData->hState4x4PropMatr);
		WaveHandleModified(pSRWRadStructAccessData->w4x4PropMatr);
	}

// Radiation Moments
	if(pSRWRadStructAccessData->wMomX != NIL)
	{
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wMomX, 1, &eStep, &eStart)) return result;
		HSetState((Handle)(pSRWRadStructAccessData->wMomX), pSRWRadStructAccessData->hStateMomX);
		WaveHandleModified(pSRWRadStructAccessData->wMomX);
	}
	if(pSRWRadStructAccessData->wMomZ != NIL)
	{
		if(result = MDSetWaveScaling(pSRWRadStructAccessData->wMomZ, 1, &eStep, &eStart)) return result;
		HSetState((Handle)(pSRWRadStructAccessData->wMomZ), pSRWRadStructAccessData->hStateMomZ);
		WaveHandleModified(pSRWRadStructAccessData->wMomZ);
	}

// Auxiliary Wave Front Data
	if(pSRWRadStructAccessData->wWfrAuxData != NIL)
	{
		pSRWRadStructAccessData->UpdateSrwWfrAuxData();
		
		HSetState((Handle)(pSRWRadStructAccessData->wWfrAuxData), pSRWRadStructAccessData->hStateWfrAuxData);
		WaveHandleModified(pSRWRadStructAccessData->wWfrAuxData);
	}

// Add here treatment of new Rad data members, if any.

	WaveHandleModified(pSRWRadStructAccessData->wRad);

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetSRWRadStructAndOptElemNames(srTIgorRadPropagInputStruct* pRadPropagInputStruct, srTSRWRadStructAccessData* pRadStructAccessData, srTStringVect* pStringVect)
{
	int result;

	if(pRadPropagInputStruct->wRad == 0) return IMPROPER_RADIATION_STRUCTURE;

	srTHandleOfSRWRadStruct HandleOfSRWRadStruct;
	HandleOfSRWRadStruct.wRad = pRadPropagInputStruct->wRad;
	SetHandleOfSRWRadStruct(&HandleOfSRWRadStruct);

	if(result = GetSRWRadStructAccessData(pRadStructAccessData)) return result;
	if(result = GetVectorOfStrings(pRadPropagInputStruct->wOptElem, pStringVect)) return result;

	return 0;

//#else

	//todo
	//return 0;
}
#endif

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetSRWRadStructAndOptElemNames(srTIgorRadPropagStokesMultiElecInputStruct* pRadPropagInputStruct, srTSRWRadStructAccessData* pRadStructAccessData, srTStringVect* pStringVect)
{
	int result;

	if(pRadPropagInputStruct->wRad == 0) return IMPROPER_RADIATION_STRUCTURE;

	srTHandleOfSRWRadStruct HandleOfSRWRadStruct;
	HandleOfSRWRadStruct.wRad = pRadPropagInputStruct->wRad;
	SetHandleOfSRWRadStruct(&HandleOfSRWRadStruct);

	if(result = GetSRWRadStructAccessData(pRadStructAccessData)) return result;
	if(result = GetVectorOfStrings(pRadPropagInputStruct->wOptElem, pStringVect)) return result;

	return 0;

//#else

	//todo
	//return 0;
}
#endif

//*************************************************************************

int srTSend::GetOptElemNames(srTStringVect& StringVect)
{
#ifdef __IGOR_PRO__

	int result;

	waveHndl wavH;
	if(pIgorWfrEmitPropagInputStruct != 0) wavH = pIgorWfrEmitPropagInputStruct->wOptElem;
	//else if(pIgorWfrEmitPropagInputStruct != 0) wavH = ...;

	if(result = GetVectorOfStrings(wavH, &StringVect)) return result;

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetVectorOfStrings(waveHndl& wavH, srTStringVect* pStringVect)
{
#ifdef __IGOR_PRO__

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	int result;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	if(numDimensions > 1) return TEXT_WAVE_1D_REQUIRED;

	int AmOfRecords = int(*dimensionSizes);
	long RadIndices[MAX_DIMENSIONS];
	for(int k=0; k<AmOfRecords; k++)
	{
		Handle textH = NewHandle(0L);
		*RadIndices = k;
		if(result = MDGetTextWavePointValue(wavH, RadIndices, textH)) return result;
		if(PtrAndHand("\0", textH, 1)) return MEMORY_ALLOCATION_FAILURE;
		if(textH != 0)
		{
			if(*textH != 0)
			{
				char* aString = new char[256];
				if(aString == 0) return MEMORY_ALLOCATION_FAILURE;
				char *tString = aString, *tH = *textH;
				for(int i=0; i<255; i++) 
				{
					*(tString++) = *tH;
					if(*(tH++) == '\0') break;
				}
				pStringVect->push_back(aString);
			}
		}

		DisposeHandle(textH);
	}

	WaveHandleModified(wavH);
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::UpdateTextWave(waveHndl& wavH, srTStringVect* pStringVect)
{
#ifdef __IGOR_PRO__

	if(wavH == NIL) return NOWAV;
	int waveType = WaveType(wavH);
	if(waveType != TEXT_WAVE_TYPE) return TEXT_WAVE_REQUIRED;

	int result;
	long numDimensions;
	long dimensionSizes[MAX_DIMENSIONS+1];
	if(result = MDGetWaveDimensions(wavH, &numDimensions, dimensionSizes)) return result;
	if(numDimensions > 1) return TEXT_WAVE_1D_REQUIRED;

	long RadIndices[MAX_DIMENSIONS+1];
	int AmOfNewRecords = pStringVect->size();
	int AmOfOldRecords = dimensionSizes[0];
	int AmOfRecords = (AmOfOldRecords > AmOfNewRecords)? AmOfOldRecords : AmOfNewRecords;
	char EmptyStr[] = " ";
	for(int k=0; k<AmOfRecords; k++)
	{
		char *pStr;
		if(k < AmOfNewRecords) pStr = (*pStringVect)[k];
		else pStr = EmptyStr;

		int AmOfBytes = strlen(pStr);
		Handle textH = NewHandle(AmOfBytes);
		strncpy(*textH, pStr, AmOfBytes);
		*RadIndices = k;
		if(result = MDSetTextWavePointValue(wavH, RadIndices, textH)) return result;
		DisposeHandle(textH);
	}

	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::GetVectorOfStrings(char* StructName, srTStringVect* pStringVect)
{
#ifdef __IGOR_PRO__

	waveHndl wavH = FetchWave(StructName);
	return GetVectorOfStrings(wavH, pStringVect);

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

int srTSend::FinishWorkingWithWave(srTWaveAccessData* pWaveAccessData)
{
#ifdef __IGOR_PRO__

	int result;
	if(pWaveAccessData->pWaveData == 0) return 0;

	HSetState((Handle)(pWaveAccessData->wHndl), pWaveAccessData->hState);
	if(pWaveAccessData->AmOfDims > -1)
	{
		DOUBLE xStart = pWaveAccessData->DimStartValues[0];
		DOUBLE xStep = pWaveAccessData->DimSteps[0];
		if((xStart < 1.E+23) && (xStep < 1.E+23))
			if(result = MDSetWaveScaling(pWaveAccessData->wHndl, ROWS, &xStep, &xStart)) return result;
		if(*(pWaveAccessData->DimUnits[0]) != '\0')
			if(result = MDSetWaveUnits(pWaveAccessData->wHndl, ROWS, pWaveAccessData->DimUnits[0])) return result;

		if(pWaveAccessData->AmOfDims > 1)
		{
			DOUBLE yStart = pWaveAccessData->DimStartValues[1];
			DOUBLE yStep = pWaveAccessData->DimSteps[1];
			if((yStart < 1.E+23) && (yStep < 1.E+23))
				if(result = MDSetWaveScaling(pWaveAccessData->wHndl, COLUMNS, &yStep, &yStart)) return result;
			if(*(pWaveAccessData->DimUnits[1]) != '\0')
				if(result = MDSetWaveUnits(pWaveAccessData->wHndl, COLUMNS, pWaveAccessData->DimUnits[1])) return result;
		}
		if(pWaveAccessData->AmOfDims > 2)
		{
			DOUBLE zStart = pWaveAccessData->DimStartValues[2];
			DOUBLE zStep = pWaveAccessData->DimSteps[2];
			if((zStart < 1.E+23) && (zStep < 1.E+23))
				if(result = MDSetWaveScaling(pWaveAccessData->wHndl, LAYERS, &zStep, &zStart)) return result;
			if(*(pWaveAccessData->DimUnits[2]) != '\0')
				if(result = MDSetWaveUnits(pWaveAccessData->wHndl, LAYERS, pWaveAccessData->DimUnits[2])) return result;
		}
		if(pWaveAccessData->AmOfDims > 3)
		{
			DOUBLE tStart = pWaveAccessData->DimStartValues[3];
			DOUBLE tStep = pWaveAccessData->DimSteps[3];
			if((tStart < 1.E+23) && (tStep < 1.E+23))
				if(result = MDSetWaveScaling(pWaveAccessData->wHndl, CHUNKS, &tStep, &tStart)) return result;
			if(*(pWaveAccessData->DimUnits[3]) != '\0')
				if(result = MDSetWaveUnits(pWaveAccessData->wHndl, CHUNKS, pWaveAccessData->DimUnits[3])) return result;
		}
	}
	if(*(pWaveAccessData->DataUnits) != '\0')
		if(result = MDSetWaveUnits(pWaveAccessData->wHndl, -1, pWaveAccessData->DataUnits)) return result;

	WaveHandleModified(pWaveAccessData->wHndl);

	pWaveAccessData->pWaveData = 0;
	return 0;

#else

	return 0;

#endif
}


//*************************************************************************

#ifdef __IGOR_PRO__

int srTSend::ShowRadSect1D(srTRadSect1D& Sect1D, char ExOrEz)
{
	srTWaveAccessData ShowingWave;
	ShowingWave.pWaveData = (char*)((ExOrEz == 'x')? Sect1D.pEx : Sect1D.pEz);
	(ShowingWave.WaveType)[0] = 'c';
	(ShowingWave.WaveType)[1] = 'f';
	ShowingWave.AmOfDims = 1;
	(ShowingWave.DimSizes)[0] = Sect1D.np;

	(ShowingWave.DimStartValues)[0] = Sect1D.ArgStart;
	(ShowingWave.DimSteps)[0] = Sect1D.ArgStep;
	
	strcpy((ShowingWave.DimUnits)[0], "m");
	strcpy(ShowingWave.NameOfWave, Sect1D.NameOfWave);

	return MakeWaveToShowData(ShowingWave);
}

#endif

//*************************************************************************

int srTSend::MakeWaveAccordingToWaveAccessData(srTWaveAccessData& ShowingWave)
{
#ifdef __IGOR_PRO__

	int result;
	long dimensionSizes[MAX_DIMENSIONS+1];
	for(int i=0; i<=MAX_DIMENSIONS; i++) dimensionSizes[i] = 0;
	for(int k=0; k<ShowingWave.AmOfDims; k++) dimensionSizes[k] = ShowingWave.DimSizes[k];

	int type;
	if(ShowingWave.WaveType[0] == 'c')
	{
		if(ShowingWave.WaveType[1] == 'f') type = (NT_FP32 | NT_CMPLX);
		else if(ShowingWave.WaveType[1] == 'd') type = (NT_FP64 | NT_CMPLX);
	}
	else if(ShowingWave.WaveType[0] == 'f') type = NT_FP32;
	else if(ShowingWave.WaveType[0] == 'd') type = NT_FP64;

	//if(result = MDMakeWave(&(ShowingWave.wHndl), ShowingWave.NameOfWave, NIL, dimensionSizes, type, 1)) return result;
	if(result = MDMakeWave(&(ShowingWave.wHndl), ShowingWave.NameOfWave, NIL, dimensionSizes, type, 1)) return FAILED_TO_CREATE_WAVE;

	for(int j=0; j<ShowingWave.AmOfDims; j++)
	{
		DOUBLE StepJ = ShowingWave.DimSteps[j];
		DOUBLE StartJ = ShowingWave.DimStartValues[j];
		if(result = MDSetWaveUnits(ShowingWave.wHndl, j, (ShowingWave.DimUnits)[j])) return result;
		if(result = MDSetWaveScaling(ShowingWave.wHndl, j, &StepJ, &StartJ)) return result;
	}

	if(ShowingWave.pWaveData == 0)
	{
		long dataOffset;
		if(result = MDAccessNumericWaveData(ShowingWave.wHndl, kMDWaveAccessMode0, &dataOffset)) return result;
		ShowingWave.hState = 0; //MoveLockHandle(ShowingWave.wHndl);
		ShowingWave.pWaveData = (char*)(*(ShowingWave.wHndl)) + dataOffset;
	}
	return 0;

#else

	//todo
	return 0;

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__

int srTSend::MakeWaveToShowData(srTWaveAccessData& ShowingWave)
{
	int result;
	if(result = MakeWaveAccordingToWaveAccessData(ShowingWave)) return result;

	long dataOffset;
	if(result = MDAccessNumericWaveData(ShowingWave.wHndl, kMDWaveAccessMode0, &dataOffset)) return result;
	ShowingWave.hState = 0; //MoveLockHandle(ShowingWave.wHndl);

	long TotAmOfData = (ShowingWave.DimSizes)[0];
	for(int jj=1; jj<ShowingWave.AmOfDims; jj++) TotAmOfData *= (ShowingWave.DimSizes)[jj];

	char IsComplex = (ShowingWave.WaveType[0] == 'c');
	if((IsComplex && (ShowingWave.WaveType[1] == 'f')) || (ShowingWave.WaveType[0] == 'f'))
	{
		float* tOutF = (float*)((char*)(*(ShowingWave.wHndl)) + dataOffset);
		float* tF = (float*)(ShowingWave.pWaveData);

		for(long ii=0; ii<TotAmOfData; ii++)
		{
			*(tOutF++) = *(tF++);
			if(IsComplex) *(tOutF++) = *(tF++);
		}
	}
	if((IsComplex && (ShowingWave.WaveType[1] == 'd')) || (ShowingWave.WaveType[0] == 'd'))
	{
		DOUBLE* tOutD = (DOUBLE*)((char*)(*(ShowingWave.wHndl)) + dataOffset);
		DOUBLE* tD = (DOUBLE*)(ShowingWave.pWaveData);

		for(long ii=0; ii<TotAmOfData; ii++)
		{
			*(tOutD++) = *(tD++);
			if(IsComplex) *(tOutD++) = *(tD++);
		}
	}
	HSetState((Handle)(ShowingWave.wHndl), ShowingWave.hState);
	WaveHandleModified(ShowingWave.wHndl);

	return 0;
}

#endif

//*************************************************************************

int srTSend::FetchNumWave(char* WaveName, srTWaveAccessData* pWaveAccessData)
{
#ifdef __IGOR_PRO__

	waveHndl wHndl = FetchWave(WaveName);
	if(wHndl == 0) return IMPROPER_OPTICAL_COMPONENT_STRUCTURE;
	srTHandleOfOneWaveStruct HandleOfOneWaveStruct;
	HandleOfOneWaveStruct.WaveHndl = wHndl;
	return GetWaveAccessData(&HandleOfOneWaveStruct, pWaveAccessData);

#else

	return 0;

#endif
}

//*************************************************************************

void srTSend::UpdateInterface()
{
#ifdef __IGOR_PRO__

	DoUpdate();

#else

	//?????

#endif
}

//*************************************************************************

#ifdef __IGOR_PRO__
int srTSend::GetRadSourceData(void* pVoidIn, char* StructNameIn, CHGenObj& hSrcOut)
{
	if((pVoidIn == 0) || (StructNameIn == 0))  return INCORRECT_ARGUMENTS_3DVIEWING;

	waveHndl wavH_ElBeam = 0;
	waveHndl wavH_MagFld = 0;

	if(strcmp(StructNameIn, "srTIgor3DViewInputStruct") == 0)
	{
		srTIgor3DViewInputStruct* pIn = (srTIgor3DViewInputStruct*)pVoidIn;
		if(pIn->wElectronBeam != 0) wavH_ElBeam = pIn->wElectronBeam;
		if(pIn->wMagFld != 0) wavH_MagFld = pIn->wMagFld;
	}
	//else...
	//Add more checking here if necessary

	int res = 0;

	srTEbmDat* pEbmDat = new srTEbmDat();
	if(res = GetTotalElectronBeamDataFormat3_FromIgorWave(*pEbmDat, wavH_ElBeam)) return res;
	CHGenObj hElBeam(pEbmDat);

	//srTMagElem* pMagElem = new srTMagElem();
	//if(res = GetTotalElectronBeamDataFormat3_FromIgorWave(*pMagElem, wavH_MagFld)) return res;

	//dddddddddddddddddddddddddddddddddd

/*
Setup Magnetic Field; distinguish types using data units field of (text) waves
	if(result = MDGetWaveUnits(wavH, -1, Units)) return result;
	strcpy(pWaveAccessData->DataUnits, Units);

int srTSend::GetGenMagFieldDataFormat1(srTGenTrjHndl& hGenTrj)
{
#if defined(__IGOR_PRO__)

	char MagFieldType; // 1: arbitrary, 2: periodic, 3: constant
	int result;
	if(result = IdentifyMagFieldType(MagFieldType)) return result;


*/


	//CSourceSR(CHGenObj In_hElBeam, CHGenObj In_hMagFld, char* Name = 0)

	//CHGenObj hSrcLoc;
/*
	DOUBLE TransvScaleFact;
	waveHndl wOptElem;
	waveHndl wObservation;
	waveHndl wMagFld;
	waveHndl wElectronBeam;
*/
	return res;
}
#endif

//*************************************************************************
