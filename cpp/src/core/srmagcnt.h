/************************************************************************//**
 * File: srmagcnt.h
 * Description: Magnetic field container (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2010
 *
 * Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMAGCNT_H
#define __SRMAGCNT_H

#include "srobject.h"
#include "objcont.h"
#include "srmagfld.h"

//-------------------------------------------------------------------------

extern CObjCont<CGenObject> gSRObjects;

struct SRWLStructMagneticFieldContainer;
typedef struct SRWLStructMagneticFieldContainer SRWLMagFldC;

//-------------------------------------------------------------------------

class srTMagFldCont : public srTMagElem {

	CObjCont<CGenObject> gMagElems;
	//vector<TVector3d> mVectCenP;

public:

	srTMagFldCont(int* pMagElem, int nMagElem)
	{
		if((pMagElem == 0) || (nMagElem <= 0)) return;
		AddElements(pMagElem, nMagElem);
	}
	srTMagFldCont(const SRWLMagFldC& inMagCnt, const TVector3d&); //SRWLIB
	srTMagFldCont() {}

	void AddElements(int* pMagElemInd, int nMagElem)
	{
		if((pMagElemInd == 0) || (nMagElem <= 0)) return;

		int* tInd = pMagElemInd;
		for(int i=0; i<nMagElem; i++)
		{
			if(!gSRObjects.exists(*tInd)) throw SR_OBJECT_DOES_NOT_EXIST;
			CSmartPtr<CGenObject> hObj = gSRObjects.get(*tInd);
			if(dynamic_cast<srTMagElem*>(hObj.ptr()) == 0) throw OBJECT_IS_NOT_MAG;
			gMagElems.insert(*tInd, hObj);
		}
	}

	void AddElement(const CHGenObj& hMagElem)
	{
		gMagElems.insert(hMagElem);
	}

	int AmountOfMembers() { return gMagElems.size();}

	void ComputeParticlePropagMatrix(double s, TMatrix2d& Mx, TMatrix2d& Mz) // virtual
	{
		Mx.Str0.x = Mx.Str1.y = 1;
		Mx.Str0.y = Mx.Str1.x = 0;
		Mz.Str0.x = Mz.Str1.y = 1;
		Mz.Str0.y = Mz.Str1.x = 0;

		if(gMagElems.size() <= 0) return;
		for(CMHGenObj::const_iterator iter = gMagElems.data.begin(); iter != gMagElems.data.end(); ++iter)
        {
			srTMagElem* pMag = (srTMagElem*)(((*iter).second).rep);
            TMatrix2d MxNext(1,0,1,0), MzNext(1,0,1,0);
			pMag->ComputeParticlePropagMatrix(s, MxNext, MzNext);
			Mx *= MxNext; Mz *= MzNext; 
		}
	}

	void compB(TVector3d& inP, TVector3d& outB) //virtual in srTMagElem
	{
		if(gMagElems.data.empty()) return;

		TVector3d relP = inP - mCenP;
		for(CMHGenObj::const_iterator iter = gMagElems.data.begin(); iter != gMagElems.data.end(); ++iter)
		{
			((srTMagElem*)((*iter).second.rep))->compB(relP, outB);
		}
	}

    void ComputeSR_Stokes(srTEbmDat* pElecBeam, srTWfrSmp* pWfrSmp, void* pPrcPar, srTStokesStructAccessData* pStokes); //virtual
	void FilterOutTrUnifMagFld(srTMagFldCont& MagTrUnifCont, srTMagFldCont& MagOptCont);
	void PrepareContForParticlePropag();
	void SortContVsStartPos();
	void DetermineLongStartAndEndPos();
};

//-------------------------------------------------------------------------

#endif
