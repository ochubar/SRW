/************************************************************************//**
 * File: srwlib.cpp
 * Description: C/C++ API of Synchrotron Radiation Workshop (SRW) Library
 * Project: Synchrotron Radiation Workshop
 * First release: October 2010
 *
 * SRW is Copyright (c) European Synchrotron Radiation Facility, Grenoble, France
 * SRW C/C++ API (SRWLIB) is Copyright (c) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, G.Geloni, L.Samoylova
 * @version see srwlUtiVerNo
 ***************************************************************************/

#include "srwlib.h"
#include "srerror.h"
#include "srmagfld.h"
#include "srmagcnt.h"
#include "srmlttsk.h"
#include "srtrjdat3d.h"
#include "srradind.h"
#include "srradint.h"
#include "srradmnp.h"
#include "sroptcnt.h"
#include "srgsnbm.h"
#include "srpersto.h"
#include "srpowden.h"
#include "srisosrc.h"
#include "srmatsta.h"

//#include <time.h> //Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP

//-------------------------------------------------------------------------
// Global Variables (used in SRW/SRWLIB, some may be obsolete)
//-------------------------------------------------------------------------

#ifdef __IGOR_PRO__
extern vector<int> gVectWarnNos;
#else
vector<int> gVectWarnNos;
srTYield srYield;
int gCallSpinProcess = 1;
int SpinProcess() { return 0;}
int (*pgWfrExtModifFunc)(int Action, srTSRWRadInData* pWfrIn, char PolComp); //to be replaced by gpWfrModifFunc
int (*pgOptElemGetInfByNameFunc)(const char* sNameOptElem, char** pDescrStr, int* LenDescr, void*); //to be replaced or removed
#endif

//-------------------------------------------------------------------------
// Global Function Pointers
//-------------------------------------------------------------------------

int (*gpWfrModifFunc)(int action, SRWLWfr* pWfrIn, char pol) = 0;
char* (*gpAllocArrayFunc)(char type, long long len) = 0; //OC15082018
int (*gpCompProgressIndicFunc)(double curVal) = 0;

//-------------------------------------------------------------------------
// Auxiliary Functions
//-------------------------------------------------------------------------

void UtiWarnCheck()
{
	if(!gVectWarnNos.empty())
	{
		int CurWarnNo = gVectWarnNos[0];
        gVectWarnNos.erase(gVectWarnNos.begin(), gVectWarnNos.end());
		throw CurWarnNo;
	}
}

//-------------------------------------------------------------------------
// SRW C API (SRWLIB) Functions
//-------------------------------------------------------------------------

EXP int CALL srwlUtiVerNo(char* verNoStr, int code)
{//to modify at each new release!
	if(verNoStr == 0) return SRWL_NO_FUNC_ARG_DATA;

	const char strCurrenVersionSRW[] = "3.964";
	const char strCurrenVersionSRWLIB[] = "0.055";

	const char *pStr=0;
	switch(code) {
		case 1:
			pStr = strCurrenVersionSRW; break;
		case 2:
			pStr = strCurrenVersionSRWLIB; break;
	}

	strcpy(verNoStr, pStr);
	return 0;
}

//-------------------------------------------------------------------------

EXP void CALL srwlUtiSetWfrModifFunc(int (*pExtFunc)(int action, SRWLWfr* pWfrIn, char pol))
{
	//if(pExtFunc != 0) gpWfrModifFunc = pExtFunc;
	gpWfrModifFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP void CALL srwlUtiSetAllocArrayFunc(char* (*pExtFunc)(char type, long long len)) //OC15082018
{
	gpAllocArrayFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP void CALL srwlUtiSetProgrIndFunc(int (*pExtFunc)(double curVal))
{
	//if(pExtFunc != 0) gpCompProgressIndicFunc = pExtFunc;
	gpCompProgressIndicFunc = pExtFunc;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiGetErrText(char* t, int errNo)
{
	CErrWarn srwlErWar;

	if(t == 0) return 0;
	if(errNo > 0) 
	{
		//const char* sErr = CErrWarn::GetError(errNo);
		const char* sErr = srwlErWar.GetError(errNo);
		if(sErr != 0) strcpy(t, sErr);
	}
	//else if(errNo < 0) strcpy(t, CErrWarn::GetWarning(errNo));
	else if(errNo < 0) 
	{
		//strcpy(t, CErrWarn::GetWarning(warNo));
		//const char* sWar = CErrWarn::GetWarning(errNo);
		const char* sWar = srwlErWar.GetWarning(errNo);
		if(errNo != 0) strcpy(t, sWar); //OC130312
	}
	return 0;
}

//-------------------------------------------------------------------------

bool TryToCopyMagFldInsteadOfInterp(SRWLMagFldC* pDispMagFld, SRWLMagFldC* pMagFld)
{
	if((pDispMagFld == 0) || (pMagFld == 0)) throw SRWL_NO_FUNC_ARG_DATA;
	if((pDispMagFld->nElem != 1) || (pDispMagFld->arMagFldTypes[0] != 'a')) throw SRWL_INCORRECT_PARAM_FOR_MAG_FLD_COMP;

	if((pMagFld->nElem != 1) ||
	   (pMagFld->arMagFldTypes[0] != 'a') ||
	   (pDispMagFld->arMagFldTypes[0] != 'a') ||
	   (pMagFld->arXc[0] != pDispMagFld->arXc[0]) ||
	   (pMagFld->arYc[0] != pDispMagFld->arYc[0]) ||
	   (pMagFld->arZc[0] != pDispMagFld->arZc[0])) return false;
	   
	if(((pMagFld->arVx != 0) && (pDispMagFld->arVx == 0)) ||
	   ((pMagFld->arVx == 0) && (pDispMagFld->arVx != 0))) return false;
	else if((pMagFld->arVx != 0) && (pDispMagFld->arVx != 0))
	{
		if(pMagFld->arVx[0] != pDispMagFld->arVx[0]) return false;
	}
	if(((pMagFld->arVy != 0) && (pDispMagFld->arVy == 0)) ||
	   ((pMagFld->arVy == 0) && (pDispMagFld->arVy != 0))) return false;
	else if((pMagFld->arVy != 0) && (pDispMagFld->arVy != 0))
	{
		if(pMagFld->arVy[0] != pDispMagFld->arVy[0]) return false;
	}
	if(((pMagFld->arVz != 0) && (pDispMagFld->arVz == 0)) ||
	   ((pMagFld->arVz == 0) && (pDispMagFld->arVz != 0))) return false;
	else if((pMagFld->arVz != 0) && (pDispMagFld->arVz != 0))
	{
		if(pMagFld->arVz[0] != pDispMagFld->arVz[0]) return false;
	}
		
	if(((pMagFld->arAng != 0) && (pDispMagFld->arAng == 0)) ||
	   ((pMagFld->arAng == 0) && (pDispMagFld->arAng != 0))) return false;
	else if((pMagFld->arAng != 0) && (pDispMagFld->arAng != 0))
	{
		if(pMagFld->arAng[0] != pDispMagFld->arAng[0]) return false;
	}

	SRWLMagFld3D *pOutFld3D = (SRWLMagFld3D*)(*(pDispMagFld->arMagFld));
	SRWLMagFld3D *pInFld3D = (SRWLMagFld3D*)(*(pMagFld->arMagFld));

	if((pOutFld3D->nx == pInFld3D->nx) &&
	   (pOutFld3D->ny == pInFld3D->ny) &&
	   (pOutFld3D->nz == pInFld3D->nz) &&
	   (pOutFld3D->nRep == pInFld3D->nRep))
	{
		if((pInFld3D->arX != 0) && (pOutFld3D->arX != 0))
		{
			double *arXin = pInFld3D->arX, *arXout = pOutFld3D->arX;
			for(int i = 0; i < (pInFld3D->nx); i++)
			{
				if(arXin[i] != arXout[i]) return false;
			}
		}
		else if(((pInFld3D->arX != 0) && (pOutFld3D->arX == 0)) ||
				((pInFld3D->arX == 0) && (pOutFld3D->arX != 0)) ||
				 (pOutFld3D->rx != pInFld3D->rx)) return false;

		if((pInFld3D->arY != 0) && (pOutFld3D->arY != 0))
		{
			double *arYin = pInFld3D->arY, *arYout = pOutFld3D->arY;
			for(int i = 0; i < (pInFld3D->ny); i++)
			{
				if(arYin[i] != arYout[i]) return false;
			}
		}
		else if(((pInFld3D->arY != 0) && (pOutFld3D->arY == 0)) ||
				((pInFld3D->arY == 0) && (pOutFld3D->arY != 0)) ||
				 (pOutFld3D->ry != pInFld3D->ry)) return false;

		if((pInFld3D->arZ != 0) && (pOutFld3D->arZ != 0))
		{
			double *arZin = pInFld3D->arZ, *arZout = pOutFld3D->arZ;
			for(int i = 0; i < (pInFld3D->nz); i++)
			{
				if(arZin[i] != arZout[i]) return false;
			}
		}
		else if(((pInFld3D->arZ != 0) && (pOutFld3D->arZ == 0)) ||
				((pInFld3D->arZ == 0) && (pOutFld3D->arZ != 0)) ||
				 (pOutFld3D->rz != pInFld3D->rz)) return false;

		if(((pInFld3D->arBx != 0) && (pOutFld3D->arBx == 0)) ||
		   ((pInFld3D->arBx == 0) && (pOutFld3D->arBx != 0))) return false;

		if(((pInFld3D->arBy != 0) && (pOutFld3D->arBy == 0)) ||
		   ((pInFld3D->arBy == 0) && (pOutFld3D->arBy != 0))) return false;

		if(((pInFld3D->arBz != 0) && (pOutFld3D->arBz == 0)) ||
		   ((pInFld3D->arBz == 0) && (pOutFld3D->arBz != 0))) return false;

		bool BxIsDef = ((pInFld3D->arBx != 0) && (pOutFld3D->arBx != 0));
		bool ByIsDef = ((pInFld3D->arBy != 0) && (pOutFld3D->arBy != 0));
		bool BzIsDef = ((pInFld3D->arBz != 0) && (pOutFld3D->arBz != 0));
		if(!(BxIsDef || ByIsDef || BzIsDef)) return false;

		long long nTot = (pInFld3D->nx)*(pInFld3D->ny)*(pInFld3D->nz);
		double *tBxIn = pInFld3D->arBx, *tBxOut = pOutFld3D->arBx;
		double *tByIn = pInFld3D->arBy, *tByOut = pOutFld3D->arBy;
		double *tBzIn = pInFld3D->arBz, *tBzOut = pOutFld3D->arBz;

		double cx = 1., cy = 1.; //OC25012018
		if(pMagFld->arPar3 != 0) cx = pMagFld->arPar3[0]; 
		if(pMagFld->arPar4 != 0) cy = pMagFld->arPar4[0];
		for(long long i = 0; i < nTot; i++)
		{
			//if(BxIsDef) *(tBxOut++) = *(tBxIn++);
			if(BxIsDef) *(tBxOut++) = (*(tBxIn++))*cx; //OC25012018
			//if(ByIsDef) *(tByOut++) = *(tByIn++);
			if(ByIsDef) *(tByOut++) = (*(tByIn++))*cy; //OC25012018
			if(BzIsDef) *(tBzOut++) = *(tBzIn++);
		}
		return true;
	}
	return false;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcMagFld(SRWLMagFldC* pDispMagFld, SRWLMagFldC* pMagFld, double* precPar)
{
	if((pDispMagFld == 0) || (pMagFld == 0)) return SRWL_NO_FUNC_ARG_DATA;
	if((pDispMagFld->nElem != 1) || (pDispMagFld->arMagFldTypes[0] != 'a')) return SRWL_INCORRECT_PARAM_FOR_MAG_FLD_COMP;
	try 
	{
		int typeCalc = int(precPar[0]);
		if((typeCalc > 0) && (typeCalc < 3))
		{//testing if parameter(s) coinside with one of mesh values, then the field has to be just copied without any interpolation
			bool skipIndSearch = (bool)precPar[5];
			if(skipIndSearch)
			{
				if(TryToCopyMagFldInsteadOfInterp(pDispMagFld, pMagFld))
				{
					UtiWarnCheck(); return 0;
				}
			}
		}

		//TVector3d vZeroCenP(0,0,0);
		//srTMagFldCont magCont(*pMagFld, vZeroCenP);
		TVector3d vZero(0,0,0);
		srTMagFldCont magCont(*pMagFld, vZero, vZero); //OC170615

		SRWLMagFld3D *pFld3D = (SRWLMagFld3D*)(pDispMagFld->arMagFld[0]);

		TVector3d vDispCenP(0.,0.,0.);
		if((pDispMagFld->arXc != 0) && (pDispMagFld->arYc != 0) && (pDispMagFld->arZc != 0))
		{
			vDispCenP.x = pDispMagFld->arXc[0]; 
			vDispCenP.y = pDispMagFld->arYc[0]; 
			vDispCenP.z = pDispMagFld->arZc[0];
		}
		TVector3d vDispAxV(0.,0.,0.);
		if((pDispMagFld->arVx != 0) && (pDispMagFld->arVy != 0) && (pDispMagFld->arVz != 0))
		{
			vDispAxV.x = pDispMagFld->arVx[0]; 
			vDispAxV.y = pDispMagFld->arVy[0]; 
			vDispAxV.z = pDispMagFld->arVz[0];
		}
		double ang = 0.;
		if(pDispMagFld->arAng != 0) ang = pDispMagFld->arAng[0];

		//srTMagFld3d magFld3d(pFld3D->rx, pFld3D->nx, pFld3D->ry, pFld3D->ny, pFld3D->rz, pFld3D->nz, pFld3D->arX, pFld3D->arY, pFld3D->arZ, pFld3D->arBx, pFld3D->arBy, pFld3D->arBz, pFld3D->nRep, pFld3D->interp, 0, vDispCenP);
		srTMagFld3d magFld3d(pFld3D->rx, pFld3D->nx, pFld3D->ry, pFld3D->ny, pFld3D->rz, pFld3D->nz, pFld3D->arX, pFld3D->arY, pFld3D->arZ, pFld3D->arBx, pFld3D->arBy, pFld3D->arBz, pFld3D->nRep, pFld3D->interp, 0, vDispCenP, vDispAxV, ang); //OC170615
		
		//int typeCalc = int(precPar[0]);
		if(typeCalc == 0) magFld3d.tabulateB(&magCont);
		else if((typeCalc > 0) && (typeCalc < 3)) magFld3d.tabInterpB(magCont, precPar, pMagFld->arPar1, pMagFld->arPar2, pMagFld->arPar3, pMagFld->arPar4);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{ 
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcPartTraj(SRWLPrtTrj* pTrj, SRWLMagFldC* pMagFld, double* precPar)
{//may modify pTrj->ctStart, pTrj->ctEnd !
	if((pTrj == 0) || (pMagFld == 0)) return SRWL_NO_FUNC_ARG_DATA;
	if((pTrj->arX == 0) || (pTrj->arXp == 0) || (pTrj->arY == 0) || (pTrj->arYp == 0) || (pTrj->np <= 0)) return SRWL_INCORRECT_TRJ_STRUCT;

	try 
	{
		const double elecEn0 = 0.51099890221e-03; //[GeV]
		SRWLParticle &part = pTrj->partInitCond;
		double arMom1[] = {(part.gamma)*(part.relE0)*elecEn0, part.x, part.xp, part.y, part.yp, part.z};
		srTEbmDat elecBeam(1., 1., arMom1, 6, 0, 0, part.z, part.nq); //to check what is more appropriate for s0: part.z or pTrj->ctStart

		//TVector3d vZeroCenP(0,0,0);
		//CSmartPtr<srTMagElem> hMagElem(new srTMagFldCont(*pMagFld, vZeroCenP));
		TVector3d vZero(0,0,0);
		CSmartPtr<srTMagElem> hMagElem(new srTMagFldCont(*pMagFld, vZero, vZero)); //OC170615

		double &sSt = pTrj->ctStart, &sEn = pTrj->ctEnd;
		if(sSt >= sEn) 
		{//set sSt and sEn to full range of magnetic field defined
			hMagElem.rep->GetMagnFieldLongLim(sSt, sEn);

			//OC150815 (commented-out)
			//double totLengthFract = (sEn - sSt)*0.01;
			//sSt -= totLengthFract; sEn += totLengthFract;

			if(sSt > part.z) sSt = part.z; //OC130312
			if(sEn < part.z) sEn = part.z;

			sSt -= part.z; sEn -= part.z; //to make sure that ct = 0 at z = part.z
		}

		srTGenTrjDat genTrjDat(&elecBeam, hMagElem);
		genTrjDat.CompTrjCrdVel(pTrj->ctStart, pTrj->ctEnd, pTrj->np, precPar, pTrj->arXp, pTrj->arX, pTrj->arYp, pTrj->arY, pTrj->arZp, pTrj->arZ, pTrj->arBx, pTrj->arBy, pTrj->arBz);

		//hMagElem.rep->DeallocAuxData(); //could be moved to CompTrjCrdVel?
		//not necessary; should be called from destructor

		UtiWarnCheck();
	}
	catch(int erNo) 
	{ 
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcPartTrajFromKickMatr(SRWLPrtTrj* pTrj, SRWLKickM* arKickM, int nKickM, double* precPar)
{
	if((pTrj == 0) || (arKickM == 0) || (nKickM <= 0)) return SRWL_NO_FUNC_ARG_DATA;
	if((pTrj->arX == 0) || (pTrj->arXp == 0) || (pTrj->arY == 0) || (pTrj->arYp == 0) || (pTrj->np <= 0)) return SRWL_INCORRECT_TRJ_STRUCT;

	try 
	{
		const double elecEn0 = 0.51099890221e-03; //[GeV]
		SRWLParticle &part = pTrj->partInitCond;
		double arMom1[] = {(part.gamma)*(part.relE0)*elecEn0, part.x, part.xp, part.y, part.yp, part.z};
		srTEbmDat elecBeam(1., 1., arMom1, 6, 0, 0, part.z, part.nq); //to check what is more appropriate for s0: part.z or pTrj->ctStart

		srTGenTrjDat genTrjDat(&elecBeam);
		genTrjDat.CompTrjKickMatr(arKickM, nKickM, pTrj->ctStart, pTrj->ctEnd, pTrj->np, precPar, pTrj->arXp, pTrj->arX, pTrj->arYp, pTrj->arY, pTrj->arZp, pTrj->arZ);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{ 
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcElecFieldSR(SRWLWfr* pWfr, SRWLPrtTrj* pTrj, SRWLMagFldC* pMagFld, double* precPar, int nPrecPar)
{
	if((pWfr == 0) || (precPar == 0)) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;

	bool trjIsDefined = true;
	if(pTrj == 0) trjIsDefined = false;
	else if((((pTrj->arX == 0) || (pTrj->arXp == 0)) && ((pTrj->arY == 0) || (pTrj->arYp == 0))) || (pTrj->np <= 0)) trjIsDefined = false;

	bool fldIsDefined = true;
	if(pMagFld == 0) fldIsDefined = false;
	else if((pMagFld->arMagFld == 0) || (pMagFld->arMagFldTypes == 0) || (pMagFld->nElem <= 0)) fldIsDefined = false;

	if((!trjIsDefined) && (!fldIsDefined)) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;
	int locErNo = 0;

	try 
	{
		if(!trjIsDefined)
		{
			pTrj = new SRWLPrtTrj();

			int npTraj = 100000;
			//int npTraj = 10000; //OCTEST160815
			if((nPrecPar <= 0) || (nPrecPar > 4)) npTraj = (int)precPar[4];
			pTrj->arX = new double[npTraj];
			pTrj->arXp = new double[npTraj];
			pTrj->arY = new double[npTraj];
			pTrj->arYp = new double[npTraj];
			pTrj->arZ = new double[npTraj]; //required?
			pTrj->arZp = new double[npTraj]; //required?
			pTrj->partInitCond = pWfr->partBeam.partStatMom1;

			pTrj->np = npTraj;

			double sStartInt = 0.; //longitudinal position to start integration
			if((nPrecPar <= 0) || (nPrecPar > 2)) sStartInt = precPar[2];

			double sEndInt = 0.; //longitudinal position to end integration
			if((nPrecPar <= 0) || (nPrecPar > 3)) sEndInt = precPar[3];

			//pTrj->ctStart = precPar[2] - pTrj->partInitCond.z; //precPar[2]; //OC_GIANLUCA
			//pTrj->ctEnd = precPar[3] - pTrj->partInitCond.z;//precPar[3]; //OC_GIANLUCA
			//processing of the case pTrj->ctStart >= pTrj->ctEnd takes place in srwlCalcPartTraj
			pTrj->ctStart = sStartInt - pTrj->partInitCond.z; //precPar[2]; //OC_GIANLUCA
			pTrj->ctEnd = sEndInt - pTrj->partInitCond.z;//precPar[3]; //OC_GIANLUCA

			double* precParForTrj = 0; //assume default for the moment; to update later (?)
			if(locErNo = srwlCalcPartTraj(pTrj, pMagFld, precParForTrj)) throw locErNo;
		}
		else pWfr->partBeam.partStatMom1 = pTrj->partInitCond;

		srTTrjDat trjData(pTrj); //this calculates interpolating structure required for SR calculation
		trjData.EbmDat.SetCurrentAndMom2(pWfr->partBeam.Iavg, pWfr->partBeam.arStatMom2, 21);

		srTSRWRadStructAccessData wfr(pWfr, &trjData, precPar); //ATTENTION: this may request for changing numbers of points in the wavefront mesh

		srTWfrSmp auxSmp;
		wfr.SetObsParamFromWfr(auxSmp);
		//precPar = array('d', [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp])
		
		char calcTerminTerms = 1; //by default do calculate two terminating terms 
		if((nPrecPar <= 0) || (nPrecPar > 5)) calcTerminTerms = (char)precPar[5];

		//srTParPrecElecFld precElecFld((int)precPar[0], precPar[1], precPar[2], precPar[3], precPar[6]);
		//srTParPrecElecFld(int In_IntegMethNo, double In_RelPrecOrStep, double In_sStartInt, double In_sEndInt, double In_NxNzOversamplingFactor, bool In_ShowProgrIndic = true, char In_CalcTerminTerms = 1)
		srTParPrecElecFld precElecFld((int)precPar[0], precPar[1], precPar[2], precPar[3], precPar[6], false, calcTerminTerms);

        srTRadInt RadInt;
		RadInt.ComputeElectricFieldFreqDomain(&trjData, &auxSmp, &precElecFld, &wfr, 0);
		wfr.OutSRWRadPtrs(*pWfr);
		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		locErNo = erNo;
		//return erNo;
	}
	if(!trjIsDefined)
	{
		if(pTrj->arX != 0) { delete[] pTrj->arX; pTrj->arX = 0;}
		if(pTrj->arXp != 0) { delete[] pTrj->arXp; pTrj->arXp = 0;}
		if(pTrj->arY != 0) { delete[] pTrj->arY; pTrj->arY = 0;}
		if(pTrj->arYp != 0) { delete[] pTrj->arYp; pTrj->arYp = 0;}
		if(pTrj->arZ != 0) { delete[] pTrj->arZ; pTrj->arZ = 0;}
		if(pTrj->arZp != 0) { delete[] pTrj->arZp; pTrj->arZp = 0;}
		delete pTrj; 
	}
	return locErNo;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcElecFieldGaussian(SRWLWfr* pWfr, SRWLGsnBm* pGsnBm, double* precPar)
{
	if((pWfr == 0) || (pGsnBm == 0)) return SRWL_INCORRECT_PARAM_FOR_GAUS_BEAM_COMP;

	int locErNo = 0;
	try 
	{
		double arMom1[] = {pGsnBm->x, pGsnBm->xp, pGsnBm->y, pGsnBm->yp};
		srTGsnBeam GsnBm(-1, pGsnBm->polar, pGsnBm->sigX, pGsnBm->mx, pGsnBm->sigY, pGsnBm->my, pGsnBm->sigT, 1, arMom1, pGsnBm->z, pGsnBm->repRate, pGsnBm->pulseEn, pGsnBm->avgPhotEn);
		
		srTSRWRadStructAccessData wfr(pWfr, &GsnBm, precPar); //ATTENTION: this may request for changing numbers of points in the wavefront mesh
		srTWfrSmp auxSmp;
		wfr.SetObsParamFromWfr(auxSmp);

		GsnBm.ComputeElectricField(&auxSmp, &wfr);
		wfr.OutSRWRadPtrs(*pWfr);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		locErNo = erNo;
		//return erNo;
	}
	return locErNo;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcElecFieldPointSrc(SRWLWfr* pWfr, SRWLPtSrc* pPtSrc, double* precPar)
{
	if(pWfr == 0) return SRWL_INCORRECT_PARAM_FOR_SPHER_WAVE_COMP;

	int locErNo = 0;
	try 
	{
		srTIsotrSrc IsotrSrc(pPtSrc);
		srTSRWRadStructAccessData wfr(pWfr, pPtSrc->z, precPar); //ATTENTION: this may request for changing numbers of points in the wavefront mesh

		//srTWfrSmp auxSmp;
		//wfr.SetObsParamFromWfr(auxSmp);

		IsotrSrc.ComputeElectricField(wfr);
		wfr.OutSRWRadPtrs(*pWfr);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		locErNo = erNo;
		//return erNo;
	}
	return locErNo;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcStokesUR(SRWLStokes* pStokes, SRWLPartBeam* pElBeam, SRWLMagFldU* pUnd, double* precPar)
{
	if((pStokes == 0) || (pElBeam == 0) || (pUnd == 0)) return SRWL_INCORRECT_PARAM_FOR_SR_COMP;

	int locErNo = 0;
	try 
	{
		const double elecEn0 = 0.51099890221e-03; //[GeV]
		SRWLParticle &part = pElBeam->partStatMom1;
		double arMom1[] = {(part.gamma)*(part.relE0)*elecEn0, part.x, part.xp, part.y, part.yp, part.z};
		double s0 = part.z; //??? //this corresponds to ct = 0
		srTEbmDat eBeam(pElBeam->Iavg, pElBeam->nPart, arMom1, 6, pElBeam->arStatMom2, 21, s0, part.nq);

		//TVector3d inUndCenP(0,0,0);
		//srTMagFieldPeriodic und(*pUnd, inUndCenP);
		TVector3d vZero(0,0,0);
		srTMagFieldPeriodic und(*pUnd, vZero, vZero); //OC170615

		//srTWfrSmp wfrSmp(pStokes->zStart, pStokes->xStart, pStokes->xFin, pStokes->nx, pStokes->yStart, pStokes->yFin, pStokes->ny, 0, pStokes->eStart, pStokes->eFin, pStokes->ne, "eV");
		SRWLStructRadMesh &mesh = pStokes->mesh;
		srTWfrSmp wfrSmp(mesh.zStart, mesh.xStart, mesh.xFin, mesh.nx, mesh.yStart, mesh.yFin, mesh.ny, 0, mesh.eStart, mesh.eFin, mesh.ne, "eV");

		//Default precision parameters:
		int initHarm = 1;
		int finHarm = 31;
		double precS = 1.;
		double precPhi = 1.;
		char int_or_flux = 'f';
		if(precPar != 0)
		{
			initHarm = (int)precPar[0];
			finHarm = (int)precPar[1];
			precS = precPar[2];
			precPhi = precPar[3];
			int normType = (int)precPar[4];
			if(normType == 1) int_or_flux = 'f';
			else if(normType == 2) int_or_flux = 'i';
			else if(normType == 3) int_or_flux = 'a';
		}
		srTParPrecStokesPer auxPrecPar(initHarm, finHarm, precS, precPhi, int_or_flux);

		srTRadIntPeriodic radInt(&eBeam, &und, &wfrSmp, &auxPrecPar);
		if(locErNo = radInt.ComputeTotalStokesDistr(0, pStokes)) throw locErNo;

		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		locErNo = erNo;
	}

	return locErNo;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcPowDenSR(SRWLStokes* pStokes, SRWLPartBeam* pElBeam, SRWLPrtTrj* pTrj, SRWLMagFldC* pMagFld, double* precPar)
{
	if((pStokes == 0) || (pElBeam == 0)) return SRWL_INCORRECT_PARAM_FOR_SR_POW_COMP;

	bool trjIsDefined = true;
	if(pTrj == 0) trjIsDefined = false;
	else if((((pTrj->arX == 0) || (pTrj->arXp == 0)) && ((pTrj->arY == 0) || (pTrj->arYp == 0))) || (pTrj->np <= 0)) trjIsDefined = false;

	bool fldIsDefined = true;
	if(pMagFld == 0) fldIsDefined = false;
	else if((pMagFld->arMagFld == 0) || (pMagFld->arMagFldTypes == 0) || (pMagFld->nElem <= 0)) fldIsDefined = false;

	if((!trjIsDefined) && (!fldIsDefined)) return SRWL_INCORRECT_PARAM_FOR_SR_POW_COMP;

	int locErNo = 0;
	try 
	{
		if(!trjIsDefined)
		{
			int npTraj = 100000; //default
			if(precPar != 0) npTraj = (int)precPar[4];

			pTrj = new SRWLPrtTrj();
			pTrj->arX = new double[npTraj];
			pTrj->arXp = new double[npTraj];
			pTrj->arY = new double[npTraj];
			pTrj->arYp = new double[npTraj];
			pTrj->arZ = new double[npTraj]; //required?
			pTrj->arZp = new double[npTraj]; //required?
			pTrj->np = npTraj;
			pTrj->partInitCond = pElBeam->partStatMom1;

			pTrj->ctStart = -pTrj->partInitCond.z; //?
			pTrj->ctEnd = -pTrj->partInitCond.z; //?
			if(precPar != 0)
			{
				if(precPar[2] < precPar[3])
				{
					pTrj->ctStart += precPar[2];
					pTrj->ctEnd += precPar[3];
					//processing of the case pTrj->ctStart >= pTrj->ctEnd takes place in srwlCalcPartTraj
				}
			}
			double* precParForTrj = 0; //assume default for the moment; to update later (?)
			if(locErNo = srwlCalcPartTraj(pTrj, pMagFld, precParForTrj)) throw locErNo;
		}

		srTTrjDat trjData(pTrj); //this calculates interpolating structure required for SR calculation
		trjData.EbmDat.SetCurrentAndMom2(pElBeam->Iavg, pElBeam->arStatMom2, 21);

		//Default precision parameters:
		double precFact = 1.;
		int meth = 1; //"near field"
		int useSpecIntLim = 0;
		double sIntStart = 0;
		double sIntFin = 0;
		if(precPar != 0)
		{
			precFact = precPar[0];
			meth = (int)precPar[1];
			useSpecIntLim = (precPar[2] < precPar[3])? 1 : 0;
			sIntStart = precPar[2];
			sIntFin = precPar[3];
		}
		srTParPrecPowDens precPowDens(meth, precFact, useSpecIntLim, sIntStart, sIntFin);
		//srTWfrSmp wfrSmp(pStokes->zStart, pStokes->xStart, pStokes->xFin, pStokes->nx, pStokes->yStart, pStokes->yFin, pStokes->ny, 0, pStokes->eStart, pStokes->eFin, pStokes->ne, "eV");
		SRWLStructRadMesh &mesh = pStokes->mesh;
		//srTWfrSmp wfrSmp(mesh.zStart, mesh.xStart, mesh.xFin, mesh.nx, mesh.yStart, mesh.yFin, mesh.ny, 0, mesh.eStart, mesh.eFin, mesh.ne, "eV");
		double inNormObsPlane[] = {mesh.nvx, mesh.nvz,  mesh.nvy}; //Order of the coordinates is in accord with SRW for Igor Pro; check sign!
		double horOrtObsPlane[] = {mesh.hvx, mesh.hvz, mesh.hvy}; //Order of the coordinates is in accord with SRW for Igor Pro; check sign!
		srTWfrSmp wfrSmp(mesh.zStart, mesh.xStart, mesh.xFin, mesh.nx, mesh.yStart, mesh.yFin, mesh.ny, mesh.arSurf, mesh.eStart, mesh.eFin, mesh.ne, "eV", 0, 0, 0, 0, horOrtObsPlane, inNormObsPlane);

		srTPowDensStructAccessData pow;
		pow.pBasePowDens = (float*)(pStokes->arS0); //consider calculating all Stokes parameters of the power density
		pStokes->numTypeStokes = 'f';

		srTRadIntPowerDensity RadIntPowDens;
		RadIntPowDens.ComputePowerDensity(&trjData, &wfrSmp, &precPowDens, &pow);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		locErNo = erNo;
	}
	if(!trjIsDefined)
	{
		if(pTrj->arX != 0) { delete[] pTrj->arX; pTrj->arX = 0;}
		if(pTrj->arXp != 0) { delete[] pTrj->arXp; pTrj->arXp = 0;}
		if(pTrj->arY != 0) { delete[] pTrj->arY; pTrj->arY = 0;}
		if(pTrj->arYp != 0) { delete[] pTrj->arYp; pTrj->arYp = 0;}
		if(pTrj->arZ != 0) { delete[] pTrj->arZ; pTrj->arZ = 0;}
		if(pTrj->arZp != 0) { delete[] pTrj->arZp; pTrj->arZp = 0;}
		delete pTrj; 
	}
	return locErNo;
}

//-------------------------------------------------------------------------

EXP int CALL srwlCalcIntFromElecField(char* pInt, SRWLWfr* pWfr, char polar, char intType, char depType, double e, double x, double y)
{
	//double start; //Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP
	//get_walltime (&start);

	if((pWfr == 0) || (pInt == 0)) return SRWL_INCORRECT_PARAM_FOR_INT_EXTR;

	try 
	{
		srTSRWRadStructAccessData wfr(pWfr);
		CHGenObj hWfr(&wfr, true);
		srTRadGenManip radGenManip(hWfr);
		
		////Re-defining intType
		////from SRWL convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Flux; 3- Multi-Elec. Flux; 4- Single-Elec. Rad. Phase; 5- Re(E); 6- Im(E); 7- Time or Photon Energy Integrated Intensity
		////to old SRW convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Rad. Phase; 3- Re(E); 4- Single-Elec. Flux; 5- Multi-Elec. Flux; 6- Im(E); 7- Time or Photon Energy Integrated Intensity
		//if(intType == 2) intType = 4;
		//else if(intType == 3) intType = 5;
		//else if(intType == 4) intType = 2;
		//else if(intType == 5) intType = 3;
		//radGenManip.ExtractRadiation((int)polar, (int)intType, (int)depType, wfr.Pres, e, x, y, pInt);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlCalcIntFromElecField : before ExtractRadiation",&start);

		radGenManip.ExtractRadiationSRWL(polar, intType, depType, wfr.Pres, e, x, y, pInt); //OC19082018

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlCalcIntFromElecField : ExtractRadiation",&start);

		//wfr.OutSRWRadPtrs(*pWfr); //not necessary?
		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlResizeElecField(SRWLWfr* pWfr, char type, double* par)
{
	if((pWfr == 0) || (par == 0)) return SRWL_INCORRECT_PARAM_FOR_RESIZE;
	bool isCorA = (type == 'c') || (type == 'C') || (type == 'a') || (type == 'A');
	bool isForT = (type == 'f') || (type == 'F') || (type == 't') || (type == 'T');
	if(!(isCorA || isForT)) return SRWL_INCORRECT_PARAM_FOR_RESIZE;
	//if((type != 'c') && (type != 'C') && (type != 'a') && (type != 'A') && (type != 'f') && (type != 'F') && (type != 't') && (type != 'T')) return SRWL_INCORRECT_PARAM_FOR_RESIZE;
	//if((type != 'c') && (type != 'C') && (type != 'a') && (type != 'A') /*&& (type != 'f') && (type != 'F') && (type != 't') && (type != 'T')*/) return SRWL_INCORRECT_PARAM_FOR_RESIZE;
	//resizing vs photon energy / time yet to be implemented!

	try 
	{
		int locErNo = 0;
		srTGenOptElem GenOptElem;
		srTSRWRadStructAccessData wfr(pWfr);

		srTRadResize resPar;
		//resPar.UseOtherSideFFT = (char)par[0];
		resPar.useOtherSideFFT((int)par[0]);

		if(isCorA)
		{
			resPar.pxm = (double)par[1]; //horizontal range
			resPar.pxd = (double)par[2]; //horizontal resolution
			resPar.pzm = (double)par[3]; //vertical range
			resPar.pzd = (double)par[4]; //vertical resolution
			//OC071014
			resPar.RelCenPosX = (double)par[5]; //relative horizontal wavefront center position (default is 0.5)
			resPar.RelCenPosZ = (double)par[6]; //relative vertical wavefront center position (default is 0.5)

			if(locErNo = GenOptElem.RadResizeGen(wfr, resPar)) return locErNo; //to make it static or put it into srTSRWRadStructAccessData eventually?
		}
		else if(isForT)
		{
			resPar.pem = (double)par[1]; //photon energy / time range
			resPar.ped = (double)par[2]; //photon energy / time resolution
			//OC071014
			resPar.RelCenPosE = (double)par[3]; //relative photon energy / time center position (default is 0.5)

			if(locErNo = GenOptElem.RadResizeGenE(wfr, resPar)) return locErNo; //to make it static or put it into srTSRWRadStructAccessData eventually?
		}

		wfr.OutSRWRadPtrs(*pWfr);
		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlSetRepresElecField(SRWLWfr* pWfr, char repr)
{
	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime (&start);

	if(pWfr == 0) return SRWL_INCORRECT_PARAM_FOR_CHANGE_REP;
	
	char reprCoordOrAng=0, reprFreqOrTime=0;
	if((repr == 'c') || (repr == 'C') || (repr == 'a') || (repr == 'A')) reprCoordOrAng = repr;
	if((repr == 'f') || (repr == 'F') || (repr == 't') || (repr == 'T')) reprFreqOrTime = repr;
	if((!reprCoordOrAng) && (!reprFreqOrTime)) return SRWL_INCORRECT_PARAM_FOR_CHANGE_REP;

	try 
	{
		srTSRWRadStructAccessData wfr(pWfr);

		int locErNo = 0;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlSetRepresElecField : before fft",&start);

		if(reprCoordOrAng) locErNo = wfr.SetRepresCA(reprCoordOrAng); //set Coordinate or Angular representation
		else if(reprFreqOrTime) locErNo = wfr.SetRepresFT(reprFreqOrTime); //set Frequency or Time representation
		
		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlSetRepresElecField : wfr.SetRepresFT",&start);
		
		if(locErNo) return locErNo;

		wfr.OutSRWRadPtrs(*pWfr);

		UtiWarnCheck();
	}
	catch(int erNo) 
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlPropagElecField(SRWLWfr* pWfr, SRWLOptC* pOpt, int nInt, char** arID, SRWLRadMesh* arIM, char** arI) //OC15082018
//EXP int CALL srwlPropagElecField(SRWLWfr* pWfr, SRWLOptC* pOpt)
{
	if((pWfr == 0) || (pOpt == 0)) return SRWL_INCORRECT_PARAM_FOR_WFR_PROP;
	int locErNo = 0;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime (&start);

	try 
	{
		srTCompositeOptElem optCont(*pOpt);
		srTSRWRadStructAccessData wfr(pWfr);
		if(locErNo = optCont.CheckRadStructForPropagation(&wfr)) return locErNo;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("srwlPropagElecField: CheckRadStructForPropagation",&start);

		//if(locErNo = optCont.PropagateRadiationGuided(wfr)) return locErNo;
		if(locErNo = optCont.PropagateRadiationGuided(wfr, nInt, arID, arIM, arI)) return locErNo; //OC15082018

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("srwlPropagElecField: PropagateRadiationGuided",&start);

		wfr.OutSRWRadPtrs(*pWfr);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime("srwlPropagElecField: PropagateRadiationGuided",&start);

		UtiWarnCheck();
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiFFT(char* pcData, char typeData, double* arMesh, int nMesh, int dir)
{
	if((pcData == 0) || (arMesh == 0) || (typeData != 'f') || (nMesh < 3) || (dir == 0)) return SRWL_INCORRECT_PARAM_FOR_FFT;
	//to support typeData == 'd' later
	int locErNo = 0;
	try 
	{
		long nx = (long)arMesh[2];
		if(nx <= 1) return SRWL_INCORRECT_PARAM_FOR_FFT;
		long ny = 1;
		if(nMesh >= 6) ny = (long)arMesh[5];

		int dimFFT = 1;
		if(ny > 1) dimFFT = 2;

		float *pfData = (float*)pcData;

		if(dimFFT == 1)
		{
			CGenMathFFT1DInfo FFT1DInfo;
			FFT1DInfo.pInData = pfData;
			FFT1DInfo.pOutData = pfData; //does this ensure in-place FFT?
			FFT1DInfo.Dir = (char)dir;
			FFT1DInfo.xStart = arMesh[0];
			FFT1DInfo.xStep = arMesh[1];
			FFT1DInfo.Nx = nx;
			FFT1DInfo.HowMany = 1;
			FFT1DInfo.UseGivenStartTrValue = 0;

			CGenMathFFT1D FFT1D;
			if(locErNo = FFT1D.Make1DFFT(FFT1DInfo)) return locErNo;

			arMesh[0] = FFT1DInfo.xStartTr;
			arMesh[1] = FFT1DInfo.xStepTr;
		}
		else
		{
			CGenMathFFT2DInfo FFT2DInfo;
			FFT2DInfo.pData = pfData;
			FFT2DInfo.Dir = (char)dir;
			FFT2DInfo.xStart = arMesh[0];
			FFT2DInfo.xStep = arMesh[1];
			FFT2DInfo.Nx = nx;
			FFT2DInfo.yStart = arMesh[3];
			FFT2DInfo.yStep = arMesh[4];
			FFT2DInfo.Ny = ny;
			FFT2DInfo.UseGivenStartTrValues = 0;

			CGenMathFFT2D FFT2D;
			if(locErNo = FFT2D.Make2DFFT(FFT2DInfo)) return locErNo;

			arMesh[0] = FFT2DInfo.xStartTr;
			arMesh[1] = FFT2DInfo.xStepTr;
			arMesh[3] = FFT2DInfo.yStartTr;
			arMesh[4] = FFT2DInfo.yStepTr;
		}

		UtiWarnCheck();
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiConvWithGaussian(char* pcData, char typeData, double* arMesh, int nMesh, double* arSig)
{
	if((pcData == 0) || (arMesh == 0) || (typeData != 'f') || (nMesh < 3)) return SRWL_INCORRECT_PARAM_FOR_CONV_WITH_GAUS;
	//to support typeData == 'd' later
	int locErNo = 0;
	try 
	{
		long nx = (long)arMesh[2];
		if(nx <= 1) return SRWL_INCORRECT_PARAM_FOR_FFT;
		double xStart = arMesh[0], xStep = arMesh[1];

		long ny = 1;
		double yStart = 0., yStep = 0.;
		if(nMesh >= 6) 
		{
			yStart = arMesh[3];
			yStep = arMesh[4];
			ny = (long)arMesh[5];
		}

		int dimFFT = 1;
		if(ny > 1) dimFFT = 2;

		//Create auxiliary complex array
		//long nTot = nx*ny;
		long long nTot = ((long long)nx)*((long long)ny);

		long nxAdj = nx;
		if(((nx >> 1) << 1) != nx) nxAdj++;
		long nyAdj = ny;
		if((ny != 1) && (((ny >> 1) << 1) != ny)) nyAdj++;

		//long nTotC = (nxAdj*nyAdj) << 1;
		long long nTotC = (((long long)nxAdj)*((long long)nyAdj)) << 1;
		float *arDataC = new float[nTotC];

		float *pfData = (float*)pcData;
		float *t_arDataC = arDataC, *tfData = pfData;

		if((nxAdj == nx) && (nyAdj == ny))
		{
			//for(long i=0; i<nTot; i++) { *(t_arDataC++) = *(tfData++); *(t_arDataC++) = 0.;}
			for(long long i=0; i<nTot; i++) { *(t_arDataC++) = *(tfData++); *(t_arDataC++) = 0.;}
		}
		else
		{
			for(long iy=0; iy<ny; iy++)
			{
				//long iy_nx = iy*nx;
				long long iy_nx = iy*nx;
				for(long ix=0; ix<nx; ix++)
				{
					//long ofst = iy_nx + ix;
					//long ofst2 = ofst << 1;
					long long ofst = iy_nx + ix;
					long long ofst2 = ofst << 1;
					arDataC[ofst2] = pfData[ofst]; arDataC[ofst2+1] = 0.;
				}
				if(nxAdj > nx) 
				{
					//long ofst = iy_nx + nx;
					//long ofst2 = ofst << 1;
					long long ofst = iy_nx + nx;
					long long ofst2 = ofst << 1;
					arDataC[ofst2] = 0.; arDataC[ofst2+1] = 0.;
				}
			}
			if(nyAdj > ny)
			{
				//long ny_nx = ny*nx;
				long long ny_nx = ((long long)ny)*((long long)nx);
				for(long ix=0; ix<nx; ix++)
				{
					//long ofst = ny_nx + ix;
					//long ofst2 = ofst << 1;
					long long ofst = ny_nx + ix;
					long long ofst2 = ofst << 1;
					arDataC[ofst2] = 0.; arDataC[ofst2+1] = 0.;
				}
				if(nxAdj > nx) 
				{
					//long ofst = ny_nx + nx;
					//long ofst2 = ofst << 1;
					long long ofst = ny_nx + nx;
					long long ofst2 = ofst << 1;
					arDataC[ofst2] = 0.; arDataC[ofst2+1] = 0.;
				}
			}
		}

		arMesh[2] = nxAdj; arMesh[5] = nyAdj;
		if(locErNo = srwlUtiFFT((char*)arDataC, typeData, arMesh, nMesh, 1)) return locErNo;

		double xStartTr = arMesh[0];
		double xStepTr = arMesh[1];
		double yStartTr = arMesh[3];
		double yStepTr = arMesh[4];

		double sigX = *arSig, sigY = 0., alp = 0.;
		if(dimFFT > 1) 
		{
			sigY = *(arSig + 1);
			alp = *(arSig + 2);
		}

		const double pi = 3.14159265358979;
		//const double c0 = 2.*pi*pi;
		double c0 = 2.*pi*pi;
		double sigXe2 = sigX*sigX, sigYe2 = sigY*sigY;
		
		double cxy = 0.;
		if(alp != 0.) 
		{
			double alp_sigXe2_sigYe2 = alp*sigXe2*sigYe2;
			c0 /= (1. - alp*alp_sigXe2_sigYe2);
			cxy = c0*2*alp_sigXe2_sigYe2;
		}

		double cx = c0*sigXe2;
		double cy = c0*sigYe2;

		double y = yStartTr; //yStart;
		double argTot;
		float factExp;
		t_arDataC = arDataC;
		for(long iy=0; iy<nyAdj; iy++)
		{
			double argTermY = cy*y*y;
			double x = xStartTr; //xStart;
			for(long ix=0; ix<nxAdj; ix++)
			{
				argTot = -argTermY - cx*x*x;
				if(alp != 0.) argTot += cxy*x*y;

				factExp = (float)exp(argTot);
				*(t_arDataC++) *= factExp;
				*(t_arDataC++) *= factExp;
				x += xStepTr; //xStep;
			}
			y += yStepTr; //yStep;
		}

		if(locErNo = srwlUtiFFT((char*)arDataC, typeData, arMesh, nMesh, -1)) return locErNo;
		arMesh[2] = nx; arMesh[5] = ny;

		if((nxAdj == nx) && (nyAdj == ny))
		{
			t_arDataC = arDataC; tfData = pfData;
			//for(long i=0; i<nTot; i++) { *(tfData++) = *(t_arDataC++); t_arDataC++;}
			for(long long i=0; i<nTot; i++) { *(tfData++) = *(t_arDataC++); t_arDataC++;}
		}
		else
		{
			for(long iy=0; iy<ny; iy++)
			{
				//long iy_nx = iy*nx;
				long long iy_nx = iy*nx;
				for(long ix=0; ix<nx; ix++)
				{
					//long ofst = iy_nx + ix;
					long long ofst = iy_nx + ix;
					pfData[ofst] = arDataC[ofst << 1];
				}
			}
		}
		delete[] arDataC;
		UtiWarnCheck();
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiIntInf(double* arInf, char* pcData, char typeData, SRWLRadMesh* pMesh)
{//OC24092018
	if((arInf == 0) || (pcData == 0) || ((typeData != 'f') && (typeData != 'd')) || (pMesh == 0)) return SRWL_INCORRECT_PARAM_FOR_INT_STAT;
	try 
	{
		srTWaveAccessData InData(pcData, typeData, pMesh); //OC13112018
		if(InData.AmOfDims > 2) throw SRWL_INCORRECT_PARAM_FOR_INT_STAT; //to update in the future

/**
		srTWaveAccessData InData;
		InData.pWaveData = pcData;
		InData.WaveType[0] = typeData; InData.WaveType[1] = '\0';

		int nDims = 0;
		int n1 = 0, n2 = 0, n3 = 0;
		double start1 = 0, start2 = 0, start3 = 0;
		double step1 = 0, step2 = 0, step3 = 0;
		if(pMesh->ne > 1) 
		{
			nDims++; 
			n1 = pMesh->ne;
			start1 = pMesh->eStart;
			step1 = (pMesh->eFin - start1)/(n1 - 1);
		}
		if(pMesh->nx > 1) 
		{
			nDims++;
			if(n1 == 0) 
			{
				n1 = pMesh->nx;
				start1 = pMesh->xStart;
				step1 = (pMesh->xFin - start1)/(n1 - 1);
			}
			else 
			{
				n2 = pMesh->nx;
				start2 = pMesh->xStart;
				step2 = (pMesh->xFin - start2)/(n2 - 1);
			}
		}
		if(pMesh->ny > 1) 
		{
			nDims++;
			if(n1 == 0) 
			{
				n1 = pMesh->ny;
				start1 = pMesh->yStart;
				step1 = (pMesh->yFin - start1)/(n1 - 1);
			}
			else if(n2 == 0) 
			{
				n2 = pMesh->ny;
				start2 = pMesh->yStart;
				step2 = (pMesh->yFin - start2)/(n2 - 1);
			}
			else 
			{
				n3 = pMesh->ny;
				start3 = pMesh->yStart;
				step3 = (pMesh->yFin - start3)/(n3 - 1);
			}
		}
		if(nDims > 2) throw SRWL_INCORRECT_PARAM_FOR_INT_STAT; //to update in the future
		InData.AmOfDims = nDims;

		InData.DimSizes[0] = n1;
		InData.DimSizes[1] = n2;
		InData.DimSizes[2] = n3;
		InData.DimStartValues[0] = start1;
		InData.DimStartValues[1] = start2;
		InData.DimStartValues[2] = start3;
		InData.DimSteps[0] = step1;
		InData.DimSteps[1] = step2;
		InData.DimSteps[2] = step3;
**/

		float arInfAux[5];
		srTWaveAccessData OutData;
		OutData.pWaveData = (char*)arInfAux;
		OutData.WaveType[0] = 'f';
		OutData.AmOfDims = 1;
		OutData.DimSizes[0] = 5; //to update to 3D data case in the future
		OutData.DimSizes[1] = 0;
		OutData.DimStartValues[0] = 0;
		OutData.DimSteps[0] = 1;

		srTAuxMatStat AuxMatStat;
		int res = 0;
		if(res = AuxMatStat.FindSimplestStat(InData, OutData)) throw res;

		//re-arranging res. data for eventual 3D case
		for(int i=0; i<3; i++) arInf[i] = arInfAux[i];
		arInf[3] = 0;
		arInf[4] = arInfAux[3];
		arInf[5] = arInfAux[4];
		arInf[6] = 0;
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiIntProc(char* pcI1, char typeI1, SRWLRadMesh* pMesh1, char* pcI2, char typeI2, SRWLRadMesh* pMesh2, double* arPar)
{//OC13112018
	if((pcI1 == 0) || ((typeI1 != 'f') && (typeI1 != 'd')) || (pMesh1 == 0) || 
	   (pcI2 == 0) || ((typeI2 != 'f') && (typeI2 != 'd')) || (pMesh2 == 0) || (arPar == 0)) return SRWL_INCORRECT_PARAM_FOR_INT_PROC;

	try 
	{
		srTWaveAccessData wI1(pcI1, typeI1, pMesh1), wI2(pcI2, typeI2, pMesh2);
		srTRadGenManip::IntProc(&wI1, &wI2, arPar);
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiUndFromMagFldTab(SRWLMagFldC* pUndCnt, SRWLMagFldC* pMagCnt, double* arPrecPar)
{
	if((pUndCnt == 0) || (pMagCnt == 0) || (arPrecPar == 0)) return SRWL_INCORRECT_PARAM_FOR_CONV_MAG_2_PER;
	if((pUndCnt->nElem != 1) || (pMagCnt->nElem != 1)) return SRWL_INCORRECT_PARAM_FOR_CONV_MAG_2_PER;

	try 
	{
		SRWLMagFld3D *pFld3D = (SRWLMagFld3D*)(pMagCnt->arMagFld[0]);
		double sStart = pMagCnt->arZc[0] - 0.5*pFld3D->rz;
		double sStep = (pFld3D->nz <= 1)? 0. : pFld3D->rz/(pFld3D->nz - 1);
		srTMagFldTrUnif magFldTrUnif(sStart, sStep, pFld3D->nz, pFld3D->arBx, pFld3D->arBy, 0);

		srTMagFieldPeriodic *pMagFldPer = magFldTrUnif.CreateAndSetupMagFieldPeriodic(arPrecPar[0], (int)arPrecPar[1], arPrecPar[2]);
		SRWLMagFldU *pMagFldU = (SRWLMagFldU*)(pUndCnt->arMagFld[0]);
		pMagFldPer->SetupExtMagFldU(*pMagFldU, pUndCnt->arZc[0]);
		pUndCnt->arXc[0] = 0.; pUndCnt->arYc[0] = 0.;

		if(pMagFldPer != 0) delete pMagFldPer;
		UtiWarnCheck();
	}
	catch(int erNo) 
	{ 
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlUtiUndFindMagFldInterpInds(int* arResInds, int* pnResInds, double* arGaps, double* arPhases, int nVals, double arPrecPar[5])
{
	if((arResInds == 0) || (pnResInds == 0) || ((arGaps == 0) && (arPhases == 0)) || (nVals <= 0)) return SRWL_INCORRECT_PARAM_FOR_UND_FLD_INTERP_IND_SEARCH;

	try 
	{
		CGenMathInterp::SelectPointsForInterp1d2d(arGaps, arPhases, nVals, arResInds, *pnResInds, arPrecPar);
	}
	catch(int erNo) 
	{ 
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------

EXP int CALL srwlPropagRadMultiE(SRWLStokes* pStokes, SRWLWfr* pWfr0, SRWLOptC* pOpt, double* precPar, int (*pExtFunc)(int action, SRWLStokes* pStokesInst))
{
	if((pStokes == 0) || (pWfr0 == 0) || (pOpt == 0) || (precPar == 0)) return SRWL_INCORRECT_PARAM_FOR_WFR_PROP;
	int locErNo = 0;

	try 
	{
		srTSRWRadStructAccessData wfr(pWfr0);
		srTCompositeOptElem optCont(*pOpt);

		//to continue


		//if(locErNo = optCont.CheckRadStructForPropagation(&wfr)) return locErNo;
		//if(locErNo = optCont.PropagateRadiationGuided(wfr)) return locErNo;
		//wfr.OutSRWRadPtrs(*pWfr);


		UtiWarnCheck();
	}
	catch(int erNo)
	{
		return erNo;
	}
	return 0;
}

//-------------------------------------------------------------------------
//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
/*
void get_walltime_(double* wcTime) {
    clock_t tp;
    tp = clock();
    *wcTime = (double)(((float)tp)/CLOCKS_PER_SEC);
    // cout << "=== clock = " << tp << "   CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << "\n";
    // cout << "=== *wcTime = " << *wcTime << "\n";
}

EXP void CALL get_walltime(double* wcTime) {
  get_walltime_(wcTime);
}
*/
//-------------------------------------------------------------------------
//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
/*
EXP void CALL srwlPrintTime(const char* str, double* start){
#ifdef MANUAL_PROFILING
	double end;
	get_walltime (&end);
	double dif= end-*start;
	if (dif > 0.1)
	{
		printf ("Elapsed: %80s %5.2f s\n",str,dif);
		fflush(stdout);
	}
	*start=end;
#endif
}
*/
//-------------------------------------------------------------------------
