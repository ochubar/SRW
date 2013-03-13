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

EXP int CALL srwlCalcMagFld(SRWLMagFldC* pDispMagFld, SRWLMagFldC* pMagFld)
{
	if((pDispMagFld == 0) || (pMagFld == 0)) return SRWL_NO_FUNC_ARG_DATA;
	if((pDispMagFld->nElem != 1) || (pDispMagFld->arMagFldTypes[0] != 'a')) return SRWL_INCORRECT_PARAM_FOR_MAG_FLD_COMP;
	try 
	{
		TVector3d vZeroCenP(0,0,0);
		srTMagFldCont magCont(*pMagFld, vZeroCenP);

		SRWLMagFld3D *pFld3D = (SRWLMagFld3D*)(pDispMagFld->arMagFld[0]);
		TVector3d vDispCenP(pDispMagFld->arXc[0], pDispMagFld->arYc[0], pDispMagFld->arZc[0]);
		srTMagFld3d magFld3d(pFld3D->rx, pFld3D->nx, pFld3D->ry, pFld3D->ny, pFld3D->rz, pFld3D->nz, pFld3D->arX, pFld3D->arY, pFld3D->arZ, pFld3D->arBx, pFld3D->arBy, pFld3D->arBz, pFld3D->nRep, pFld3D->interp, 0, vDispCenP);
		magFld3d.tabulateB(&magCont);

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

		TVector3d vZeroCenP(0,0,0);
		CSmartPtr<srTMagElem> hMagElem(new srTMagFldCont(*pMagFld, vZeroCenP));

		double &sSt = pTrj->ctStart, &sEn = pTrj->ctEnd;
		if(sSt >= sEn) 
		{//set sSt and sEn to full range of magnetic field defined
			hMagElem.rep->GetMagnFieldLongLim(sSt, sEn);
			double totLengthFract = (sEn - sSt)*0.01;
			sSt -= totLengthFract; sEn += totLengthFract;

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

		TVector3d inUndCenP(0,0,0);
		srTMagFieldPeriodic und(*pUnd, inUndCenP);
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
		srTWfrSmp wfrSmp(mesh.zStart, mesh.xStart, mesh.xFin, mesh.nx, mesh.yStart, mesh.yFin, mesh.ny, 0, mesh.eStart, mesh.eFin, mesh.ne, "eV");

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
	if((pWfr == 0) || (pInt == 0)) return SRWL_INCORRECT_PARAM_FOR_INT_EXTR;

	try 
	{
		srTSRWRadStructAccessData wfr(pWfr);
		CHGenObj hWfr(&wfr, true);
		srTRadGenManip radGenManip(hWfr);
		
		//Re-defining intType
		//from SRWL convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Flux; 3- Multi-Elec. Flux; 4- Single-Elec. Rad. Phase; 5- Re(E); 6- Im(E)
		//to old SRW convention: 0- Single-Elec. Intensity; 1- Multi-Elec. Intensity; 2- Single-Elec. Rad. Phase; 3- Re(E); 4- Single-Elec. Flux; 5- Multi-Elec. Flux; 6- Im(E)
		if(intType == 2) intType = 4;
		else if(intType == 3) intType = 5;
		else if(intType == 4) intType = 2;
		else if(intType == 5) intType = 3;

		radGenManip.ExtractRadiation((int)polar, (int)intType, (int)depType, wfr.Pres, e, x, y, pInt);
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

			if(locErNo = GenOptElem.RadResizeGen(wfr, resPar)) return locErNo; //to make it static or put it into srTSRWRadStructAccessData eventually?
		}
		else if(isForT)
		{
			resPar.pem = (double)par[1]; //photon energy / time range
			resPar.ped = (double)par[2]; //photon energy / time resolution

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
	if(pWfr == 0) return SRWL_INCORRECT_PARAM_FOR_CHANGE_REP;
	
	char reprCoordOrAng=0, reprFreqOrTime=0;
	if((repr == 'c') || (repr == 'C') || (repr == 'a') || (repr == 'A')) reprCoordOrAng = repr;
	if((repr == 'f') || (repr == 'F') || (repr == 't') || (repr == 'T')) reprFreqOrTime = repr;
	if((!reprCoordOrAng) && (!reprFreqOrTime)) return SRWL_INCORRECT_PARAM_FOR_CHANGE_REP;

	try 
	{
		srTSRWRadStructAccessData wfr(pWfr);

		int locErNo = 0;
		if(reprCoordOrAng) locErNo = wfr.SetRepresCA(reprCoordOrAng); //set Coordinate or Angular representation
		else if(reprFreqOrTime) locErNo = wfr.SetRepresFT(reprFreqOrTime); //set Frequency or Time representation
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

EXP int CALL srwlPropagElecField(SRWLWfr* pWfr, SRWLOptC* pOpt)
{
	if((pWfr == 0) || (pOpt == 0)) return SRWL_INCORRECT_PARAM_FOR_WFR_PROP;
	int locErNo = 0;
	try 
	{
		srTCompositeOptElem optCont(*pOpt);
		srTSRWRadStructAccessData wfr(pWfr);
		if(locErNo = optCont.CheckRadStructForPropagation(&wfr)) return locErNo;
		if(locErNo = optCont.PropagateRadiationGuided(wfr)) return locErNo;

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
