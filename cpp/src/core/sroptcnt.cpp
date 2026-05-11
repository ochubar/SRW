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
				pOptElem = new srTZonePlate(p->nZones, p->rn, p->thick, p->atLen1, p->atLen2, p->delta1, p->delta2, p->x, p->y, p->e0); //OC22062019
				//pOptElem = new srTZonePlate(p->nZones, p->rn, p->thick, p->atLen1, p->atLen2, p->delta1, p->delta2, p->x, p->y);
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
			else if(strcmp(sType, "interferometer") == 0) pOptElem = new srTOptInterferometer(*((SRWLOptI*)(*t_arOpt))); //OC01102025

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
				if(propRes.pxm > 0.) //OC17102025
				{
					propRes.pxd = t_pr[6];
					propRes.pzm = t_pr[7];
					propRes.pzd = t_pr[8];
					if(curNumPropPar > 9) propRes.ShiftTypeBeforeRes = (char)t_pr[9];
					if(curNumPropPar > 10) propRes.xCenShift = t_pr[10];
					if(curNumPropPar > 11) propRes.zCenShift = t_pr[11];

					propRes.eStart = propRes.eStep = propRes.xStart = propRes.xStep = propRes.zStart = propRes.zStep = 0.; //OC17102025
					propRes.ne = propRes.nx = propRes.nz = 0; //OC17102025
				}
				else //OC17102025
				{
					propRes.nx = (long)(t_pr[8]);
					propRes.xStart = t_pr[6];
					propRes.xStep = (propRes.nx <= 1)? 0. : (t_pr[7] - propRes.xStart)/(propRes.nx - 1);
					propRes.nz = (long)(t_pr[11]);
					propRes.zStart = t_pr[9];
					propRes.zStep = (propRes.nz <= 1)? 0. : (t_pr[10] - propRes.zStart)/(propRes.nz - 1);
					//In Python interface, this option is used for specifying resize to a given mesh, e.g.
					//['op_FinMesh_pp', 'f',  [0, 0, 1.0, 0, 0, -1, -0.0147, 0.0144, 960, -0.0158, 0.01458, 972, 0.0, 0.0, 0.0, 0.0, 0.0, 0., 0., 0., 0.], 'final post-propagation resize parameters to a given mesh'], #OC30092025 (changed/simplified optical component names)

					propRes.eStart = propRes.eStep = 0.;
					propRes.ne = 0;
				}

				if(curNumPropPar > 12) propRes.vLxOut = t_pr[12]; //Default coordinates of the output Optical Axis vector
				if(curNumPropPar > 13) propRes.vLyOut = t_pr[13];
				if(curNumPropPar > 14) propRes.vLzOut = t_pr[14];
				if(curNumPropPar > 15) propRes.vHxOut = t_pr[15]; //Default coordinates of the Horizontal Base vector of the output frame
				if(curNumPropPar > 16) propRes.vHyOut = t_pr[16];

				//OC01102025
				if(curNumPropPar > 17) propRes.Rx = t_pr[17];
				if(curNumPropPar > 18) propRes.xc = t_pr[18];
				if(curNumPropPar > 19) propRes.Rz = t_pr[19];
				if(curNumPropPar > 20) propRes.zc = t_pr[20];
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

int srTCompositeOptElem::PropagateRadiationGuided(srTSRWRadStructAccessData& wfr, int nInt, char** arID, SRWLRadMesh* arIM, char** arI, void* pvGPU) //HG30112023
//int srTCompositeOptElem::PropagateRadiationGuided(srTSRWRadStructAccessData& wfr, int nInt, char** arID, SRWLRadMesh* arIM, char** arI) //OC15082018
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
#ifdef _OFFLOAD_GPU //HG30112023
	bool dataOnDevice = false;
	TGPUUsageArg parGPU(pvGPU); //OC18022024
	TGPUUsageArg *pGPU = &parGPU; //OC18022024
#endif

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

			//TO IMPLEMENT: eventual shift of wavefront before resizing!!!?

			//OC20102025: Distinguish when resize should be made to a precize mesh and call wfr.Resize(..) in that case
			if((curPropResizeInst.nx > 0) && (curPropResizeInst.nz > 0))
			{
				SRWLRadMesh mesh;
				mesh.nx = curPropResizeInst.nx; mesh.ny = curPropResizeInst.nz; //mesh.ne = curPropResizeInst.ne;
				//mesh.eStart = curPropResizeInst.eStart; mesh.eFin = curPropResizeInst.eStart + (curPropResizeInst.ne - 1)*curPropResizeInst.eStep;
				mesh.xStart = curPropResizeInst.xStart; mesh.xFin = curPropResizeInst.xStart + (curPropResizeInst.nx - 1)*curPropResizeInst.xStep;
				mesh.yStart = curPropResizeInst.zStart; mesh.yFin = curPropResizeInst.zStart + (curPropResizeInst.nz - 1)*curPropResizeInst.zStep;

				mesh.ne = wfr.ne; //keep energy mesh unchanged
				mesh.eStart = wfr.eStart; mesh.eFin = wfr.eStart + (wfr.ne - 1)*wfr.eStep;
					
				double arParResizeMesh[] = {0.,1.,1.}; //default
				arParResizeMesh[0] = (double)curPropResizeInst.useOtherSideFFT();

				if(curPropResizeInst.doNotTreatSpherTerm()) arParResizeMesh[1] = 0.;
				else
				{
					if(curPropResizeInst.Rx != 0.)
					{
						wfr.RobsX = curPropResizeInst.Rx;
						wfr.RobsXAbsErr = 0.01*::fabs(wfr.RobsX); //to enable treatment of the quadratic phase terms at resizing
						wfr.xc = curPropResizeInst.xc;
					}
					if(curPropResizeInst.Rz != 0.)
					{
						wfr.RobsZ = curPropResizeInst.Rz;
						wfr.RobsZAbsErr = 0.01*::fabs(wfr.RobsZ);
						wfr.zc = curPropResizeInst.zc;
					}
				}
				//To make arParResizeMesh[2] input variable? : allow or not correction of Re and Im parts of the E-field based on intensity ratio (0- don't allow, 1- allow)

#ifdef _OFFLOAD_GPU //HG191102025
				if(CAuxGPU::GPUEnabled(pGPU) && dataOnDevice)
				{
					if(wfr.pBaseRadX != NULL)
						wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2 * wfr.ne * wfr.nx * wfr.nz);
					if(wfr.pBaseRadZ != NULL)
						wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2 * wfr.ne * wfr.nx * wfr.nz);
					dataOnDevice = false;
				}
#endif
				wfr.Resize(mesh, arParResizeMesh);

			}//OC20102025
			else if((::fabs(curPropResizeInst.pxd - 1.) > tolRes) || (::fabs(curPropResizeInst.pxm - 1.) > tolRes) ||
			//if((::fabs(curPropResizeInst.pxd - 1.) > tolRes) || (::fabs(curPropResizeInst.pxm - 1.) > tolRes) ||
				//(::fabs(curPropResizeInst.pzd - 1.) > tolRes) || (::fabs(curPropResizeInst.pzm - 1.) > tolRes))
				(::fabs(curPropResizeInst.pzd - 1.) > tolRes) || (::fabs(curPropResizeInst.pzm - 1.) > tolRes) || (curPropResizeInst.ShiftTypeBeforeRes > 0)) //OC11072019
			{
				//if(res = RadResizeGen(wfr, curPropResizeInst)) return res;
				if(res = RadResizeGen(wfr, curPropResizeInst, pvGPU)) return res; //HG30112023

#ifdef _OFFLOAD_GPU //HG30112023
				//if(CAuxGPU::GPUEnabled((TGPUUsageArg*)pvGPU)) { //OC18022024 (commented-out)
				if(CAuxGPU::GPUEnabled(pGPU)) { //OC18022024
					dataOnDevice = true; //HG26022024 Add explanation: If GPU is enabled, the resized wavefront is already on the GPU, so mark it appropriately so that the data can be relocated if necessary for the first optical element
				}
#endif
			}

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

#ifdef _OFFLOAD_GPU //HG30112023
		//TGPUUsageArg* pGPU = (TGPUUsageArg*)pvGPU; //OC18022024 (commented-out)
		if(CAuxGPU::GPUEnabled(pGPU)) {
			//if(dataOnDevice && (((srTGenOptElem*)it->rep)->SupportedFeatures() & 1) == 0)
			if(dataOnDevice && (((srTGenOptElem*)it->rep)->GPUImplFeatures() & 1) == 0) //HG07022024 The optical element does not support GPU acceleration, so if necessary, copy the data back to CPU
			{
				//#if DEBUG
				//				printf("Element does not support GPU, transferring to CPU.\r\n");
				//#endif
				if(wfr.pBaseRadX != NULL)
					//wfr.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2 * wfr.ne * wfr.nx * wfr.nz * sizeof(float));
					wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2*wfr.ne*wfr.nx*wfr.nz); //HG29092025
				if(wfr.pBaseRadZ != NULL)
					//wfr.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2 * wfr.ne * wfr.nx * wfr.nz * sizeof(float));
					wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2*wfr.ne*wfr.nx*wfr.nz); //HG29092025
				dataOnDevice = false;
			}
			//else if(!dataOnDevice && (((srTGenOptElem*)it->rep)->SupportedFeatures() & 1) == 1)
			else if(!dataOnDevice && (((srTGenOptElem*)it->rep)->GPUImplFeatures() & 1) == 1) //HG07022024 The optical element supports GPU acceleration, so after it is done running, the data will be on GPU
			{
				dataOnDevice = true;
				//#if DEBUG
				//					printf("Element supports GPU, transferring...\r\n");
				//#endif
			}
		}
#endif

		srTRadResizeVect auxResizeVect;
		//if(res = ((srTGenOptElem*)(it->rep))->PropagateRadiation(&wfr, precParWfrPropag, auxResizeVect)) return res;
		//if(res = ((srTGenOptElem*)(it->rep))->PropagateRadiation(&wfr, precParWfrPropag, auxResizeVect, pvGPU)) return res; //HG30112023
		if(res = ((srTGenOptElem*)(it->rep))->PropagateRadiation(&wfr, precParWfrPropag, auxResizeVect, (((srTGenOptElem*)it->rep)->GPUImplFeatures() & 1) == 0 ? 0 : pvGPU)) return res; //HG26072024 Don't pass pvGPU to propagators that don't support it
		//maybe to use "PropagateRadiationGuided" for srTCompositeOptElem?

		//OC_DEBUG
		//std::cout << "   DEBUG: PropagateRadiationGuided: PropagateRadiation done for element:" << elemCount << "\n";
		//std::cout.flush();
		//END OC_DEBUG

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("Iteration: PropagateRadiation",&start);

		if(propIntIsNeeded)
		{
#ifdef _OFFLOAD_GPU //HG09112022 If the data is on the GPU, transfer it to CPU and synchronize before extracting the intensity
			//TGPUUsageArg* pGPU = (TGPUUsageArg*)pvGPU; //OC18022024 (commented-out)
			if(CAuxGPU::GPUEnabled(pGPU)) {
				if(dataOnDevice)
				{
					if(wfr.pBaseRadX != NULL)
						//wfr.pBaseRadX = (float*)CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2 * wfr.ne * wfr.nx * wfr.nz * sizeof(float));
						wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2*wfr.ne*wfr.nx*wfr.nz); //HG29092025
					if(wfr.pBaseRadZ != NULL)
						//wfr.pBaseRadZ = (float*)CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2 * wfr.ne * wfr.nx * wfr.nz * sizeof(float));
						wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2*wfr.ne*wfr.nx*wfr.nz); //HG29092025

					dataOnDevice = false;
				}
			}
#endif
			ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount);
		}

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

		//OC20102025: Distinguish when resize should be made to a precize mesh and call wfr->Resize(..) in that case
		if((postResize.nx > 0) && (postResize.nz > 0))
		{
			SRWLRadMesh mesh;
			mesh.nx = postResize.nx; mesh.ny = postResize.nz; //mesh.ne = postResize.ne;
			//mesh.eStart = curPropResizeInst.eStart; mesh.eFin = curPropResizeInst.eStart + (curPropResizeInst.ne - 1)*curPropResizeInst.eStep;
			mesh.xStart = postResize.xStart; mesh.xFin = postResize.xStart + (postResize.nx - 1)*postResize.xStep;
			mesh.yStart = postResize.zStart; mesh.yFin = postResize.zStart + (postResize.nz - 1)*postResize.zStep;

			mesh.ne = wfr.ne; //keep energy mesh unchanged
			mesh.eStart = wfr.eStart; mesh.eFin = wfr.eStart + (wfr.ne - 1)*wfr.eStep;

			double arParResizeMesh[] = { 0.,1.,1. }; //default
			arParResizeMesh[0] = (double)postResize.useOtherSideFFT();

			if(postResize.doNotTreatSpherTerm()) arParResizeMesh[1] = 0.;
			else
			{
				if(postResize.Rx != 0.)
				{
					wfr.RobsX = postResize.Rx;
					wfr.RobsXAbsErr = 0.01*::fabs(wfr.RobsX); //to enable treatment of the quadratic phase terms at resizing
					wfr.xc = postResize.xc;
				}
				if(postResize.Rz != 0.)
				{
					wfr.RobsZ = postResize.Rz;
					wfr.RobsZAbsErr = 0.01*::fabs(wfr.RobsZ);
					wfr.zc = postResize.zc;
				}
			}
			//To make arParResizeMesh[2] input variable? : allow or not correction of Re and Im parts of the E-field based on intensity ratio (0- don't allow, 1- allow)

#ifdef _OFFLOAD_GPU //HG191102025
			if(CAuxGPU::GPUEnabled(pGPU) && dataOnDevice)
			{
				if(wfr.pBaseRadX != NULL)
					wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2 * wfr.ne * wfr.nx * wfr.nz);
				if(wfr.pBaseRadZ != NULL)
					wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2 * wfr.ne * wfr.nx * wfr.nz);
				dataOnDevice = false;
			}
#endif
			wfr.Resize(mesh, arParResizeMesh);

		}//OC20102025
		else if((::fabs(postResize.pxd - 1.) > tolRes) || (::fabs(postResize.pxm - 1.) > tolRes) || //OC21102025
		//if((::fabs(postResize.pxd - 1.) > tolRes) || (::fabs(postResize.pxm - 1.) > tolRes) ||
		   (::fabs(postResize.pzd - 1.) > tolRes) || (::fabs(postResize.pzm - 1.) > tolRes))
			//if(res = RadResizeGen(wfr, postResize)) return res;
		{
			if(res = RadResizeGen(wfr, postResize, pvGPU)) return res; //HG26072024 make this resize able to use GPU
#ifdef _OFFLOAD_GPU
			if(CAuxGPU::GPUEnabled(pGPU)) dataOnDevice = true;
#endif
		}

#ifdef _OFFLOAD_GPU //HG26072025 Make sure the data is returned to CPU
		if(CAuxGPU::GPUEnabled(pGPU) && dataOnDevice)
		{
			if(wfr.pBaseRadX != NULL)
				wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2*wfr.ne*wfr.nx*wfr.nz);
			if(wfr.pBaseRadZ != NULL)
				wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2*wfr.ne*wfr.nx*wfr.nz);
			dataOnDevice = false;
		}
#endif
		if(propIntIsNeeded) ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount); //OC29082018
		//if(propIntIsNeeded) ExtractPropagatedIntensity(wfr, nInt, arID, arIM, arI, elemCount, nInt - 1);
	}

#ifdef _OFFLOAD_GPU //HG26072025 Make sure the data is returned to CPU
	if(CAuxGPU::GPUEnabled(pGPU) && dataOnDevice)
	{
		if(wfr.pBaseRadX != NULL)
			wfr.pBaseRadX = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadX, 2*wfr.ne*wfr.nx*wfr.nz);
		if(wfr.pBaseRadZ != NULL)
			wfr.pBaseRadZ = CAuxGPU::ToHostAndFree(pGPU, wfr.pBaseRadZ, 2*wfr.ne*wfr.nx*wfr.nz);
		dataOnDevice = false;
	}
#endif
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

			//char type = *(arID[2] + ii);
			char arIntTypeConv[] = {0,1,4,5,2,3,6,7,8}; //OC14122019
			char type = arIntTypeConv[*(arID[2] + ii)];

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
			radGenManip.ExtractRadiation(pol, type, dep, pres, mesh.eStart, mesh.xStart, mesh.yStart, arCurI); //OC13122019
			//radGenManip.ExtractRadiationSRWL(pol, type, dep, pres, mesh.eStart, mesh.xStart, mesh.yStart, arCurI);

			//Updating mesh for intensity with data from wavefront
			wfr.GetIntMesh(dep, mesh);
			//break; //OC29082018 (to enable intesity indexes in arbitrary order)
		}
	}
	return res;
}

//*************************************************************************

srTOptInterferometer::srTOptInterferometer(const SRWLOptI& opt) //OC01102025
{
	void **t_arOptC = opt.arOptC;
	if((opt.nElem <= 0) || (t_arOptC == 0)) throw UNKNOWN_OPTICAL_ELEMENT;

	for(int i=0; i<opt.nElem; i++)
	{
		if((*t_arOptC) == 0) throw UNKNOWN_OPTICAL_ELEMENT;
		srTGenOptElem *pOptElem = new srTCompositeOptElem(*((SRWLOptC*)(*t_arOptC)));

		if(pOptElem != 0)
		{
			CSmartPtr<CGenObject> hObj(pOptElem);
			GenOptCntList.push_back(hObj);
		}
		else throw UNKNOWN_OPTICAL_ELEMENT; //Maybe throw a better message?
		t_arOptC++;
	}

	//m_irec = opt.irec;
	////OC09102025: this ensures that recombination of wavefronts is always done (to be possibly changed in the future?)
	//int nBranches = (int)GenOptCntList.size();
	//if((m_irec < 0) && (nBranches > 0))
	//{
	//	m_irec = nBranches - 1;
	//}

	//OC12102025
	for(int i=0; i<3; i++) //to increase number of parameters, if necessary
	{
		m_arPar[i] = opt.arPar[i];
	}
	int nBranches = (int)GenOptCntList.size();
	if((m_arPar[0] < 0) && (nBranches > 0))
	{
		m_arPar[0] = 0; //nBranches - 1;
	}

}

//*************************************************************************

int srTOptInterferometer::PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResizeBeforeAndAfterVect, void* pvGPU) 
{//OC08102025
	if(pRadAccessData == 0) return 0;
	int nBranches = (int)GenOptCntList.size();
	if(nBranches < 1) return 0; //nothing to do

	//If number of branches > 1, duplicate wavefront for each branch, except for the "reference" one, indicated by irec parameter
	//For each branch do PropagateRadiationGuided and keep resulting wavefronts
	//If wavefront recombination is necessary, sum up electric fields from all branches with interpolation, if necessary
	//then delete wavefronts of all the branches except the reference one

	int result = 0;
	if(nBranches < 2)
	{//Only one branch - just propagate
		if(result = ((srTCompositeOptElem*)(GenOptCntList.begin()->rep))->PropagateRadiationGuided(*pRadAccessData, 0, 0, 0, 0, pvGPU)) return result;
	}

	//Copy initial wavefront (for local use only, for branches other than the "main" one)
	srTSRWRadStructAccessData origRadAccessData(pRadAccessData); //OC11121025
	//srTSRWRadStructAccessData* pOrigRadAccessData = new srTSRWRadStructAccessData(pRadAccessData);
	//if(pOrigRadAccessData == 0) return MEMORY_ALLOCATION_FAILURE;

	//First, propagate the "main" branch (m_irec)
	int irec = (int)m_arPar[0]; //OC12102025
	srTCompositeOptElem *pMainBranch = (srTCompositeOptElem*)(std::next(GenOptCntList.begin(), irec))->rep;
	//srTCompositeOptElem *pMainBranch = (srTCompositeOptElem*)(std::next(GenOptCntList.begin(), m_irec))->rep;
	if(result = pMainBranch->PropagateRadiationGuided(*pRadAccessData, 0, 0, 0, 0, pvGPU))
	{
		//if(pOrigRadAccessData != 0) delete pOrigRadAccessData; //OC11121025 (commented-out)
		return result;
	}

	int branchInd = -1;
	for(srTGenOptElemHndlList::iterator iter = GenOptCntList.begin(); iter != GenOptCntList.end(); ++iter)
	{
		branchInd++;
		if(branchInd == irec) continue; //OC12102025
		//if(branchInd == m_irec) continue;

		srTSRWRadStructAccessData curRadAccessData(&origRadAccessData); //OC11122025
		//srTSRWRadStructAccessData* pCurRadAccessData = new srTSRWRadStructAccessData(pOrigRadAccessData);
		//if(pCurRadAccessData == 0)
		//{
		//	if(pOrigRadAccessData != 0) delete pOrigRadAccessData;
		//	return MEMORY_ALLOCATION_FAILURE;
		//}

		//Propagate wavefront through the current branch
		if(result = ((srTCompositeOptElem*)((*iter).rep))->PropagateRadiationGuided(curRadAccessData, 0, 0, 0, 0, pvGPU)) //OC11122025
		//if(result = ((srTCompositeOptElem*)((*iter).rep))->PropagateRadiationGuided(*pCurRadAccessData, 0, 0, 0, 0, pvGPU))
		{
			//if(pOrigRadAccessData != 0) delete pOrigRadAccessData;
			//if(pCurRadAccessData != 0) delete pCurRadAccessData;
			return result;
		}

		//Add current wavefront to the "main" one
		pRadAccessData->AddElFieldData(curRadAccessData, m_arPar+1); //OC11122025
		//pRadAccessData->AddElFieldData(*pCurRadAccessData, m_arPar+1); //OC12102025
		//pRadAccessData->AddElFieldDataViaResize(*pCurRadAccessData);

		//if(pCurRadAccessData != 0) delete pCurRadAccessData; //OC11122025 (commented-out)
	}

	//if(pOrigRadAccessData != 0) delete pOrigRadAccessData; //OC11122025 (commented-out)
	return 0;
}

//*************************************************************************
