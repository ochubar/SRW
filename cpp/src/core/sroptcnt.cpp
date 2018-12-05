/************************************************************************//**
 * File: sroptcnt.cpp
 * Description: Optical element: Container
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "sroptcnt.h"
#include "sroptdrf.h"
#include "sroptapt.h"
#include "sroptfoc.h"
#include "sroptzp.h"
#include "sroptwgr.h"
#include "sroptgrat.h"
#include "sroptgtr.h"
#include "sropthck.h"
#include "sroptang.h"
#include "sroptcryst.h"
#include "srradmnp.h"
#include "auxparse.h"
#include "srwlib.h"

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include "stdio.h"

//*************************************************************************

extern int (*pgOptElemGetInfByNameFunc)(const char* sNameOptElem, char** pDescrStr, int* pLenDescr, void*);

//*************************************************************************

srTCompositeOptElem::srTCompositeOptElem(srTStringVect* pElemInfo, srTSRWRadStructAccessData* pRad)
{
	int AmOfMembers = (int)pElemInfo->size() - 1;

	//srTSend Send;
	int result = 0;

	const int maxStrLen = 256; //to tune
	const int maxNumParam = 200; //to tune
	char strAuxBuf[maxStrLen*maxNumParam];
	char *pDescrStr[maxNumParam];
	char *tStrAuxBuf = strAuxBuf; //**tDescrStr = pDescrStr;
	for(int i=0; i<maxNumParam; i++)
	{
		*tStrAuxBuf = '\0';
		pDescrStr[i] = tStrAuxBuf;
		tStrAuxBuf += maxStrLen;
	}
	int LenDescr = 0;
	srTDataMD OptElemNumData;

	for(int i=1; i<=AmOfMembers; i++)
	{
		char* MemberID = (*pElemInfo)[i];
		srTStringVect MemberInfo;

		//if(result = Send.GetVectorOfStrings(MemberID, &MemberInfo)) { ErrorCode = result; return;}

		OptElemNumData.pData = 0; //add other marks?
		if(result = (*pgOptElemGetInfByNameFunc)(MemberID, pDescrStr, &LenDescr, &OptElemNumData)) { ErrorCode = result; return;}
		CAuxParse::StringArr2VectCStr(pDescrStr, LenDescr, MemberInfo);

		srTGenOptElemHndl OptElemHndl;
		//if(result = OptElemSummary.SetupOpticalElement(&MemberInfo, OptElemHndl, pRad)) { ErrorCode = result; return;}
		//if(result = srTOptElemSummary::SetupOpticalElement(&MemberInfo, OptElemHndl, pRad)) { ErrorCode = result; return;}
		if(result = SetupOpticalElement(&MemberInfo, &OptElemNumData, pRad, OptElemHndl)) { ErrorCode = result; return;}

		AddOptElemBack(OptElemHndl);

		//for(int k=0; k<(int)(MemberInfo.size()); k++)
		//{
		//	char* aStr = MemberInfo[k];
		//	if(aStr != 0) delete[] aStr;
		//}
	}
}

//*************************************************************************

srTCompositeOptElem::srTCompositeOptElem(const SRWLOptC& opt)
{//to add more optical elements
	void **t_arOpt = opt.arOpt;
	char **t_arOptTypes = opt.arOptTypes;
	if((opt.nElem <= 0) || (t_arOpt == 0) || (t_arOptTypes == 0)) throw UNKNOWN_OPTICAL_ELEMENT;

	double **t_arProp = opt.arProp;
	//srTRadResize propResLast;
	//bool propResWasSet = false;

	for(int i=0; i<=opt.nElem; i++)
	{
		srTGenOptElem *pOptElem=0;
		if(i < opt.nElem)
		{
			if((*t_arOpt) == 0) throw UNKNOWN_OPTICAL_ELEMENT;
			char *sType = *t_arOptTypes;
			if(strcmp(sType, "drift") == 0)
			{
				//pOptElem = new srTDriftSpace(((SRWLOptD*)(*t_arOpt))->L); 
				SRWLOptD *p = (SRWLOptD*)(*t_arOpt);
				pOptElem = new srTDriftSpace(p->L, p->treat); 
			}
			else if((strcmp(sType, "aperture") == 0) || (strcmp(sType, "obstacle") == 0))
			{
				SRWLOptA *p = (SRWLOptA*)(*t_arOpt);
				if((p->ap_or_ob == 'a') || (p->ap_or_ob == 'A'))
				{
					if(p->shape == 'r') pOptElem = new srTRectAperture(p->Dx, p->Dy, p->x, p->y);
					else if(p->shape == 'c') pOptElem = new srTCircAperture(p->Dx, p->x, p->y);
				}
				else if((p->ap_or_ob == 'o') || (p->ap_or_ob == 'O'))
				{
					if(p->shape == 'r') pOptElem = new srTRectObstacle(p->Dx, p->Dy, p->x, p->y);
					else if(p->shape == 'c') pOptElem = new srTCircObstacle(p->Dx, p->x, p->y);
				}
				else throw UNKNOWN_OPTICAL_ELEMENT;
			}
			else if(strcmp(sType, "lens") == 0)
			{
				SRWLOptL *p = (SRWLOptL*)(*t_arOpt);
				pOptElem = new srTThinLens(p->Fx, p->Fy, p->x, p->y);
			}
			else if(strcmp(sType, "angle") == 0)
			{
				SRWLOptAng *p = (SRWLOptAng*)(*t_arOpt);
				pOptElem = new srTOptAngle(p->AngX, p->AngY);
			}
			else if(strcmp(sType, "shift") == 0)
			{
				SRWLOptShift *p = (SRWLOptShift*)(*t_arOpt);
				pOptElem = new srTOptShift(p->ShiftX, p->ShiftY);
			}
			else if((strcmp(sType, "zp") == 0) || (strcmp(sType, "ZP") == 0))
			{
				SRWLOptZP *p = (SRWLOptZP*)(*t_arOpt);
				pOptElem = new srTZonePlate(p->nZones, p->rn, p->thick, p->atLen1, p->atLen2, p->delta1, p->delta2, p->x, p->y);
			}
			else if(strcmp(sType, "waveguide") == 0)
			{
				SRWLOptWG *p = (SRWLOptWG*)(*t_arOpt);
				pOptElem = new srTWaveguideRect(p->L, p->Dx, p->Dy, p->x, p->y); 
			}
			else if(strcmp(sType, "transmission") == 0) pOptElem = new srTGenTransmission(*((SRWLOptT*)(*t_arOpt)));
			else if(strncmp(sType, "mirror", 6) == 0) pOptElem = srTMirror::DefineMirror(sType, *t_arOpt);
			else if(strcmp(sType, "grating") == 0) pOptElem = srTMirror::DefineGrating(sType, *t_arOpt);

			//else if(strcmp(sType, "mirror: plane") == 0)
			//{
			//	pOptElem = new srTMirrorPlane(*((SRWLOptMirPl*)(*t_arOpt)));
			//}
			//else if(strcmp(sType, "mirror: ellipsoid") == 0)
			//{
			//	pOptElem = new srTMirrorEllipsoid(*((SRWLOptMirEl*)(*t_arOpt)));
			//}
			//else if(strcmp(sType, "mirror: toroid") == 0)
			//{
			//	pOptElem = new srTMirrorToroid(*((SRWLOptMirTor*)(*t_arOpt)));
			//}
			//else if(strcmp(sType, "grating") == 0)
			//{
			//	SRWLOptG *p = (SRWLOptG*)(*t_arOpt);
			//	//pOptElem = new srTGrating(p->grDen, p->disPl, p->ang, p->m, p->refl);
			//	//pOptElem = new srTGrating(p->grDen, p->disPl, p->ang, p->m, p->refl, p->grDen1, p->grDen2, p->grDen3, p->grDen4);
			//}

			else if(strcmp(sType, "crystal") == 0) pOptElem = new srTOptCryst(*((SRWLOptCryst*)(*t_arOpt)));
			else if(strcmp(sType, "container") == 0) pOptElem = new srTCompositeOptElem(*((SRWLOptC*)(*t_arOpt)));
			else throw UNKNOWN_OPTICAL_ELEMENT;
		}
		if(pOptElem != 0)
		{
			CSmartPtr<CGenObject> hObj(pOptElem);
			AddOptElemBack(hObj);
		}
		if((opt.arProp != 0) && ((pOptElem != 0) || (i == opt.nElem)))
		{
			srTRadResize propRes; //with all default parameters
			//if((i > 0) && (i < opt.nElem) && (i >= opt.nProp) && propResWasSet) propRes = propResLast;
			if(i < opt.nProp)
			{
				char curNumPropPar = 9;
				if(opt.arPropN != 0) curNumPropPar = opt.arPropN[i];

				double *t_pr = *t_arProp;
				propRes.propAutoResizeBefore((int)(t_pr[0]));
				propRes.propAutoResizeAfter((int)(t_pr[1]));
				propRes.PropAutoPrec = t_pr[2];
				propRes.propAllowUnderSamp((int)(t_pr[3]));
				propRes.useOtherSideFFT((int)(t_pr[4]));
				propRes.pxm = t_pr[5];
				propRes.pxd = t_pr[6];
				propRes.pzm = t_pr[7];
				propRes.pzd = t_pr[8];
				if(curNumPropPar > 9) propRes.ShiftTypeBeforeRes = (char)t_pr[9];
				if(curNumPropPar > 10) propRes.xCenShift = t_pr[10];
				if(curNumPropPar > 11) propRes.zCenShift = t_pr[11];

				if(curNumPropPar > 12) propRes.vLxOut = t_pr[12]; //Default coordinates of the output Optical Axis vector
				if(curNumPropPar > 13) propRes.vLyOut = t_pr[13];
				if(curNumPropPar > 14) propRes.vLzOut = t_pr[14];
				if(curNumPropPar > 15) propRes.vHxOut = t_pr[15]; //Default coordinates of the Horizontal Base vector of the output frame
				if(curNumPropPar > 16) propRes.vHyOut = t_pr[16];
			}

			GenOptElemPropResizeVect.push_back(propRes); //define instructions for propagation/resizing
			//propResLast = propRes;
			//propResWasSet = true;
		}
		t_arOpt++;
		t_arOptTypes++;
		t_arProp++;
	}
}

//*************************************************************************

int srTCompositeOptElem::PropagateRadiationTest(srTSRWRadStructAccessData* pInRadAccessData, srTSRWRadStructAccessData* pOutRadAccessData)
{
	int result;
	int AmOfDrifts = 0;
	for(srTGenOptElemHndlList::iterator iter = GenOptElemList.begin(); iter != GenOptElemList.end(); ++iter)
	{
		srTDriftSpace* pDrift = dynamic_cast<srTDriftSpace*>((*iter).rep);
		if(pDrift != NULL) AmOfDrifts++;
	}
	
	if(AmOfDrifts > 1) return PROP_TEST_CONSTRAINTS;
	else 
	{
		if(AmOfDrifts == 1)
		{
			srTGenOptElemHndlList::iterator pLast = GenOptElemList.end(); pLast--;
			srTDriftSpace* pDrift = dynamic_cast<srTDriftSpace*>((*pLast).rep);
			if(pDrift == NULL) return PROP_TEST_CONSTRAINTS;
		}

		for(srTGenOptElemHndlList::iterator iter = GenOptElemList.begin(); iter != GenOptElemList.end(); ++iter)
			//if(result = ((*iter).rep)->PropagateRadiationTest(pInRadAccessData, pOutRadAccessData)) return result;
			if(result = ((srTGenOptElem*)((*iter).rep))->PropagateRadiationTest(pInRadAccessData, pOutRadAccessData)) return result;
	}
	return 0;
}

//*************************************************************************

int srTCompositeOptElem::PropagateRadiationGuided(srTSRWRadStructAccessData& wfr, int nInt, char** arID, SRWLRadMesh* arIM, char** arI) //OC15082018
//int srTCompositeOptElem::PropagateRadiationGuided(srTSRWRadStructAccessData& wfr)
{
	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start,start1;
	//get_walltime(&start);
	//get_walltime(&start1);

	int numElem = (int)GenOptElemList.size();
	int numResizeInst = (int)GenOptElemPropResizeVect.size();
	const double tolRes = 1.e-04;
	int res = 0, elemCount = 0;

	bool propIntIsNeeded = (nInt != 0) && (arID != 0) && (arI != 0); //OC27082018

	for(srTGenOptElemHndlList::iterator it = GenOptElemList.begin(); it != GenOptElemList.end(); ++it)
	{
		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("PropagateRadiationGuided: start iteration",&start1);
		//srwlPrintTime("PropagateRadiationGuided: start iteration",&start);

		int methNo = 0;
		int useResizeBefore = 0;
		int useResizeAfter = 0;
		double precFact = 1.;
		double underSampThresh = 0.5; //not user
		char analTreatment = 0;

		double vLxO=0, vLyO=0, vLzO=0; //Coordinates of the output Optical Axis vector
		double vHxO=0, vHyO=0; //Default coordinates of the Horizontal Base vector of the output frame

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("Iteration: set params",&start);

		if(elemCount < numResizeInst)
		{
			srTRadResize &curPropResizeInst = GenOptElemPropResizeVect[elemCount];
			useResizeBefore = curPropResizeInst.propAutoResizeBefore();

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime("Iteration: propAutoResizeBefore",&start);

			useResizeAfter = curPropResizeInst.propAutoResizeAfter();

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime("Iteration: propAutoResizeAfter",&start);

			if(useResizeBefore || useResizeAfter) methNo = 2;

			precFact = curPropResizeInst.PropAutoPrec;
			analTreatment = curPropResizeInst.propAllowUnderSamp();

			//TO IMPLEMENT: eventual shift of wavefront before resizing!!!

			if((::fabs(curPropResizeInst.pxd - 1.) > tolRes) || (::fabs(curPropResizeInst.pxm - 1.) > tolRes) ||
			   (::fabs(curPropResizeInst.pzd - 1.) > tolRes) || (::fabs(curPropResizeInst.pzm - 1.) > tolRes))
				if(res = RadResizeGen(wfr, curPropResizeInst)) return res;

			//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
			//srwlPrintTime("Iteration: RadResizeGen",&start);

			vLxO = curPropResizeInst.vLxOut; //OC021213
			vLyO = curPropResizeInst.vLyOut;
			vLzO = curPropResizeInst.vLzOut;
			vHxO = curPropResizeInst.vHxOut;
			vHyO = curPropResizeInst.vHyOut;
		}

		srTParPrecWfrPropag precParWfrPropag(methNo, useResizeBefore, useResizeAfter, precFact, underSampThresh, analTreatment, (char)0, vLxO, vLyO, vLzO, vHxO, vHyO);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("Iteration: precParWfrPropag",&start);

		srTRadResizeVect auxResizeVect;
		if(res = ((srTGenOptElem*)(it->rep))->PropagateRadiation(&wfr, precParWfrPropag, auxResizeVect)) return res;
		//maybe to use "PropagateRadiationGuided" for srTCompositeOptElem?

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("Iteration: PropagateRadiation",&start);

		if(propIntIsNeeded) ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount);

		elemCount++;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//char str[256];
		//sprintf(str,"%s %d","PropagateRadiationGuided: Iteration :",elemCount);
		//srwlPrintTime(str,&start1);
	}
	if(elemCount < numResizeInst)
	{//post-resize
		//TO IMPLEMENT: eventual shift of wavefront before resizing!!!

		srTRadResize &postResize = GenOptElemPropResizeVect[elemCount];

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("PropagateRadiationGuided: GenOptElemPropResizeVect",&start);

		if((::fabs(postResize.pxd - 1.) > tolRes) || (::fabs(postResize.pxm - 1.) > tolRes) ||
		   (::fabs(postResize.pzd - 1.) > tolRes) || (::fabs(postResize.pzm - 1.) > tolRes))
			if(res = RadResizeGen(wfr, postResize)) return res;

		if(propIntIsNeeded) ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount); //OC29082018
		//if(propIntIsNeeded) ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount, nInt - 1);
	}
	return 0;
}

//*************************************************************************

int srTCompositeOptElem::ExtractPropagatedIntensity(srTSRWRadStructAccessData& wfr, int nInt, char** arID, SRWLRadMesh* arIM, char** arI, int elCnt, int indIntSartSearch) //27082018
{
	if((nInt == 0) || (arID == 0) || (arI == 0)) return 0;
	int res = 0;
	int indInt = -1;
	char *pID0 = *arID;
	//char *tID0 = pID0;
	char *tID0 = pID0 + indIntSartSearch;

	for(int ii=indIntSartSearch; ii<nInt; ii++) 
	{
		if(elCnt == (int)(*(tID0++)) - 1)
		{
			char *&arCurI = *(arI + ii);
			char pol = *(arID[1] + ii);
			char type = *(arID[2] + ii);
			char dep = *(arID[3] + ii);
			char pres = *(arID[4] + ii);
			SRWLRadMesh &mesh = *(arIM + ii);

			if((nInt > 1) && (ii > 0))
			{
				if(pol < 0) pol = *(arID[1]);
				if(type < 0) type = *(arID[2]);
				if(dep < 0) dep = *(arID[3]);
				if(pres < 0) pres = *(arID[4]);
				if(mesh.ne < 0) mesh = *(arIM);
			}

			if(arCurI == 0)
			{//Allocate memory for the intensity in the front-end via a function in srTSRWRadStructAccessData (which calculates amount of necessary memory based on type of intensity)
				if(res = wfr.AllocExtIntArray(type, dep, arCurI)) return res; //OC18082018
			}
			//else //?
				
			//Extract the intensity (repeating how this is done in srwlCalcIntFromElecField)
			CHGenObj hWfr(&wfr, true);
			srTRadGenManip radGenManip(hWfr);
			radGenManip.ExtractRadiationSRWL(pol, type, dep, pres, mesh.eStart, mesh.xStart, mesh.yStart, arCurI);

			//Updating mesh for intensity with data from wavefront
			wfr.GetIntMesh(dep, mesh);
			//break; //OC29082018 (to enable intesity indexes in arbitrary order)
		}
	}
	return res;
}

//*************************************************************************
