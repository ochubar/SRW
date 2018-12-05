/************************************************************************//**
 * File: sroptcnt.h
 * Description: Optical element: Container (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SROPTCNT_H
#define __SROPTCNT_H

#include "sroptelm.h"

struct SRWLStructOpticsContainer;
typedef struct SRWLStructOpticsContainer SRWLOptC;
struct SRWLStructRadMesh;
typedef struct SRWLStructRadMesh SRWLRadMesh;

//*************************************************************************

class srTCompositeOptElem : public srTGenOptElem {

public:
	srTGenOptElemHndlList GenOptElemList;
	srTRadResizeVect GenOptElemPropResizeVect; //OC090311

	srTCompositeOptElem(srTStringVect* pElemInfo, srTSRWRadStructAccessData* pRad);
	srTCompositeOptElem(const SRWLOptC& opt);
	srTCompositeOptElem() {}

	int PropagateRadiationTest(srTSRWRadStructAccessData*, srTSRWRadStructAccessData*);
	int PropagateRadiationGuided(srTSRWRadStructAccessData& wfr, int nInt=0, char** arID=0, SRWLRadMesh* arIM=0, char** arI=0); //OC15082018
	//int PropagateRadiationGuided(srTSRWRadStructAccessData& wfr);
	int ExtractPropagatedIntensity(srTSRWRadStructAccessData& wfr, int nInt, char** arID, SRWLRadMesh* arIM, char** arI, int elCnt, int indIntSartSearch=0); //27082018

	void AddOptElemFront(srTGenOptElemHndl& OptElemHndl)
	{
		GenOptElemList.push_front(OptElemHndl);
	}
	void AddOptElemBack(srTGenOptElemHndl& OptElemHndl)
	{
		GenOptElemList.push_back(OptElemHndl);
	}

	int PropagateRadiation(srTSRWRadStructAccessData* pRadAccessData, srTParPrecWfrPropag& ParPrecWfrPropag, srTRadResizeVect& ResizeBeforeAndAfterVect)
	{
		int AmOfElem = (int)GenOptElemList.size(); //OC110104
		int ElemCount = 0; //OC110104
		char GenUseResAfter = ParPrecWfrPropag.UseResAfter; //OC110104
		ParPrecWfrPropag.UseResAfter = 1; //OC110104
		//OC141011: to reconsider this !!!

		int result;
		for(srTGenOptElemHndlList::iterator iter = GenOptElemList.begin(); iter != GenOptElemList.end(); ++iter)
		{
			ElemCount++; //OC110104
			if(ElemCount == AmOfElem) //OC110104
			{//if element is the last
				if(!GenUseResAfter) ParPrecWfrPropag.UseResAfter = GenUseResAfter;
			}

			//if(result = ((srTGenOptElem*)((*iter).rep))->PropagateRadiation(pRadAccessData, MethNo, ResizeBeforeAndAfterVect)) return result;
			if(result = ((srTGenOptElem*)((*iter).rep))->PropagateRadiation(pRadAccessData, ParPrecWfrPropag, ResizeBeforeAndAfterVect)) return result;
		}
		ParPrecWfrPropag.UseResAfter = GenUseResAfter; //OC110104
		return 0;
	}

	void AddPtrOfActualOptElem(srTGenOptElemPtrList& ActOptElemsList)
	{
		for(srTGenOptElemHndlList::iterator iter = GenOptElemList.begin(); iter != GenOptElemList.end(); ++iter)
		{
			//((*iter).rep)->AddPtrOfActualOptElem(ActOptElemsList);
			((srTGenOptElem*)((*iter).rep))->AddPtrOfActualOptElem(ActOptElemsList);
		}
	}

};

//*************************************************************************

#endif
