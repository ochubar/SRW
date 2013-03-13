/************************************************************************//**
 * File: srmagcnt.cpp
 * Description: Magnetic field container
 * Project: Synchrotron Radiation Workshop
 * First release: 2010
 *
 * Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 0.06
 ***************************************************************************/

#include "srmagcnt.h"
#include "srwlib.h"

//*************************************************************************

srTMagFldCont::srTMagFldCont(const SRWLMagFldC& inMagCnt, const TVector3d& inCenP) : srTMagElem(inCenP) //SRWLIB
{
	int nMagElem = inMagCnt.nElem;

	if((inMagCnt.arMagFld == 0) || (inMagCnt.arMagFldTypes == 0)) throw SRWL_NO_FUNC_ARG_DATA;

	void **t_arMagFld = inMagCnt.arMagFld;
	char *t_arMagFldTypes = inMagCnt.arMagFldTypes;
	double *t_arXc = inMagCnt.arXc;
	double *t_arYc = inMagCnt.arYc;
	double *t_arZc = inMagCnt.arZc;
	bool cenPointIsDefined = ((t_arXc != 0) && (t_arYc != 0) && (t_arZc != 0));

	gsStart = 1.e+23; gsEnd = -1.e+23;
	double cur_sStart = 1.e+23, cur_sEnd = -1.e+23;
	bool longFieldLimitsNotSet = true;

	int indMagForCnt = 1;
	for(int iMag=0; iMag<nMagElem; iMag++)
	{
		srTMagElem *pMagElem=0;
		TVector3d vCenP(0, 0, 0);
		if(cenPointIsDefined)
		{
			//vCenP.x = *t_arXc; vCenP.y = *t_arYc; vCenP.z = *t_arZc;
			vCenP.x = *t_arXc + inCenP.x; vCenP.y = *t_arYc + inCenP.y; vCenP.z = *t_arZc + inCenP.z; //OC30012011 shift entire group (?)
		}

		switch(*t_arMagFldTypes) 
		{
			case 'a': //arbitrary 3D
			{
				SRWLMagFld3D *pB = (SRWLMagFld3D*)(*t_arMagFld);
				//pMagElem = new srTMagFld3d(pB->rx, pB->nx, pB->ry, pB->ny, pB->rz, pB->nz, pB->arBx, pB->arBy, pB->arBz, pB->nRep, 0, vCenP);
				pMagElem = new srTMagFld3d(pB->rx, pB->nx, pB->ry, pB->ny, pB->rz, pB->nz, pB->arX, pB->arY, pB->arZ, pB->arBx, pB->arBy, pB->arBz, pB->nRep, pB->interp, 0, vCenP);
				break;
			}
			case 'm': //multipole magnet
			{
				SRWLMagFldM *pB = (SRWLMagFldM*)(*t_arMagFld);
				//pMagElem = new srTMagQuad(pB->G, pB->n_or_s, pB->Leff, pB->Ledge, vCenP);
				pMagElem = new srTMagMult(pB->G, pB->m, pB->n_or_s, pB->Leff, pB->Ledge, vCenP);
				break;
			}
			case 's': //solenoid
			{
				SRWLMagFldS *pB = (SRWLMagFldS*)(*t_arMagFld);
				pMagElem = new srTMagSol(pB->B, pB->Leff, 0, vCenP);
				break;
			}
			case 'u': //undulator
			{
				pMagElem = new srTMagFieldPeriodic(*((SRWLMagFldU*)(*t_arMagFld)), vCenP);
				break;
			}
			case 'c': //container
			{
				pMagElem = new srTMagFldCont(*((SRWLMagFldC*)(*t_arMagFld)), vCenP);
				break;
			}
		}
		if(pMagElem != 0)
		{
			CSmartPtr<CGenObject> hObj(pMagElem);
			gMagElems.insert(indMagForCnt++, hObj);

			if(longFieldLimitsNotSet) 
			{
				pMagElem->GetMagnFieldLongLim(gsStart, gsEnd);
				longFieldLimitsNotSet = false;
			}
			else
			{
				pMagElem->GetMagnFieldLongLim(cur_sStart, cur_sEnd);
				if(gsStart > cur_sStart) gsStart = cur_sStart;
				if(gsEnd < cur_sEnd) gsEnd = cur_sEnd;
			}
		}
		t_arMagFld++;
		t_arMagFldTypes++;
		t_arXc++; t_arYc++; t_arZc++;
	}
}

//*************************************************************************
