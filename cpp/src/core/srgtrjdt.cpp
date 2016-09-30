/************************************************************************//**
 * File: srgtrjdt.cpp
 * Description: Electron trajectory calculation in 3D magnetic field
 * Project: Synchrotron Radiation Workshop
 * First release: 1998
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * Portions Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srgtrjdt.h"
#include "srwlib.h"
#include "auxparse.h"
#include "gminterp.h"
#include <algorithm>

//*************************************************************************

//void srTGenTrjDat::CompTrjCrdVelRK(double sStart, double sEnd, long ns, double* pPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData)
void srTGenTrjDat::CompTrjCrdVelRK(double sStart, double sEnd, long long ns, double* pPrecPar, double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, double* pOutBxData, double* pOutByData, double* pOutBzData)
{//Corrected for SRWL:
 //Independent variable is s = c*t; initial conditions are assumed to be defined for s = c*t = 0
 //3D trajectory is calculated in Laboratory Frame
 //Z is assumed to be longitudinal variable here!
	if(ns <= 0) throw SRWL_INCORRECT_PARAM_FOR_TRJ_COMP;
	if(m_hMagElem.rep == 0) throw SRWL_INCORRECT_PARAM_FOR_TRJ_COMP;
	//if(sStart*sEnd > 0) throw SRWL_INCORRECT_INIT_COND_FOR_TRJ_COMP;

	//bool x1_is_needed = (pOutBtxData != 0);
	//bool x_is_needed = (pOutXData != 0);
	//bool z1_is_needed = (pOutBtzData != 0);
	//bool z_is_needed = (pOutZData != 0);

	double sStep = (ns <= 1)? 0 : (sEnd - sStart)/(ns - 1);

	const double chElec = 1.602176462e-19; //[C]
	const double mElec = 9.10938188e-31; //[kg]
	const double cLight = 2.99792458e+08; //[m/s]
	//double gamElec = EbmDat.Gamma;
	//m_Mult2ndDer = -chElec/(mElec*cLight*EbmDat.Gamma*sqrt(1. - 1./(EbmDat.Gamma*EbmDat.Gamma)));
	m_Mult2ndDerRK = EbmDat.nQ*chElec/(mElec*cLight*EbmDat.Gamma);

	//double sEnd = sStart + sStep*(ns - 1);
	double sMin, sMax;
	bool dataShouldBeRotatedAfterSolve = false;
	if(sStart <= sEnd)
	{
		sMin = sStart; sMax = sEnd;
	}
	else
	{
		sMin = sEnd; sMax = sStart;
		dataShouldBeRotatedAfterSolve = true;
	}
	const double sResEdgeToler = 1.E-12;
	double sAbsEdgeToler = (sMax - sMin)*sResEdgeToler;

	bool integOnLeftIsNeeded = (sMin < -sAbsEdgeToler);
	bool integOnRightIsNeeded = (sMax > sAbsEdgeToler);

	const int numEq = 5;
	
	bool onPrcRK = false;
	double *pAbsPrecParLoc = 0;
	double epsTol = 1.;
	int maxAutoSteps = 5000;
	if(pPrecPar != 0) 
	{
		int numPrecPar = (int)pPrecPar[0];
		if(numPrecPar > 0)
		{
			int methNo = (int)pPrecPar[1];
			if(methNo > 1) //to keep updated when new methods will be added
			{
				onPrcRK = true;
				pAbsPrecParLoc = pPrecPar + 2;
			}
			if(numPrecPar > 6) epsTol = pPrecPar[7];
			if(numPrecPar > 7) maxAutoSteps = (int)pPrecPar[8];
		}
	}

	CGenMathIntRungeKutta<srTGenTrjDat> gmIntRK(this, &srTGenTrjDat::funcDerivRK, numEq, onPrcRK, pAbsPrecParLoc, epsTol, maxAutoSteps); // manual mode, based on number of points
	//CGenMathIntRungeKutta<srTGenTrjDat> gmIntRK(this, &srTGenTrjDat::funcDerivRK, numEq, onPrcRK, pAbsPrecParLoc); // manual mode, based on number of points
	//CGenMathIntRungeKutta(T* ptrT, void (T::*pFuncDerivF)(double, double*, double*), int amOfEq, bool onPrc, double* precArray, double epsTol =1., int maxAutoStp =5000)
	//check how to define epsTol, maxAutoStp
	
	//double initCond[] = {EbmDat.x0, EbmDat.dxds0, EbmDat.z0, EbmDat.dzds0, EbmDat.sc};
	double initCond[] = {EbmDat.x0, EbmDat.dxds0, EbmDat.z0, EbmDat.dzds0, EbmDat.s0}; //?
	
	//long is0 = 0;
	long long is0 = 0;
	double s0Act = 0.;
	//double gamEm2 = 1./(EbmDat.Gamma*EbmDat.Gamma), btx, bty;
	double gamEm2 = EbmDat.GammaEm2, btx, bty;

	if(integOnLeftIsNeeded)
	{
		double *tOutBtxData = pOutBtxData, *tOutXData = pOutXData;
		double *tOutBtyData = pOutBtyData, *tOutYData = pOutYData;
		double *tOutBtzData = pOutBtzData, *tOutZData = pOutZData;
		
		if(sMax < -sAbsEdgeToler)
		{//need to "arrive" to sMax without storing trajectory data; step can be re-adjusted
			int auxNp = (int)fabs(sMax/sStep) + 1;
			if(auxNp <= 1) auxNp = 2;

			double *auxTrjRes = new double[(auxNp + ns)*(numEq + 1)];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			gmIntRK.solve(initCond, 0., sMax, auxNp, auxTrjRes);

			double *pResEnd = auxTrjRes + auxNp*(numEq + 1) - numEq;
			for(int i=0; i<numEq; i++) initCond[i] = pResEnd[i];

			//arrived to sMax; now solve for the entire main trajectory:
			if(ns <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				*(t_auxTrjRes++) = sMax;
				for(int i=0; i<numEq; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else gmIntRK.solve(initCond, sMax, sMin, ns, auxTrjRes);

			if(dataShouldBeRotatedAfterSolve)
			{
				double *t_auxTrjRes = auxTrjRes;
				//for(int j=0; j<ns; j++)
				for(long long j=0; j<ns; j++)
				{
					t_auxTrjRes++; //may need to be put into another place (to check)
					*(tOutXData++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtxData++) = btx;
					*(tOutYData++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtyData++) = bty;
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty)); //sign may be wrong in non-relativistic case?
				}
			}
			else
			{
				double *t_auxTrjRes = auxTrjRes + ns*(numEq + 1) - 1;
				//for(int j=0; j<ns; j++)
				for(long long j=0; j<ns; j++)
				{
					//if(pOutBtzData) *(tOutBtzData++) = *(t_auxTrjRes--);
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtyData++) = bty;
					*(tOutYData++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtxData++) = btx;
					*(tOutXData++) = *(t_auxTrjRes--);
					t_auxTrjRes--;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
		else
		{
			is0 = (int)fabs((-sMin + sAbsEdgeToler)/sStep);
			if(is0 >= ns) is0 = ns - 1;
			s0Act = sMin + is0*sStep;
			if(s0Act < -sAbsEdgeToler)
			{//one small step to s0Act
				double twoPtTrjRes[2*(numEq + 1)];
				gmIntRK.solve(initCond, 0., s0Act, 2, twoPtTrjRes);

				double *pResEnd = twoPtTrjRes + numEq + 2;
				for(int i=0; i<5; i++) initCond[i] = pResEnd[i];
			}

			//arrived to s0Act; now solve for the left part of the trajectory:
			//int nsLeft = is0 + 1;
			long long nsLeft = is0 + 1;
			double *auxTrjRes = new double[nsLeft*(numEq + 1)];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			if(nsLeft <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				*(t_auxTrjRes++) = s0Act;
				for(int i=0; i<numEq; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else gmIntRK.solve(initCond, s0Act, sMin, nsLeft, auxTrjRes);

			if(dataShouldBeRotatedAfterSolve)
			{
				double *t_auxTrjRes = auxTrjRes;
				//for(int j=0; j<nsLeft; j++)
				for(long long j=0; j<nsLeft; j++)
				{
					t_auxTrjRes++; //may need to be put into another place (to check)
					*(tOutXData++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtxData++) = btx;
					*(tOutYData++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtyData++) = bty;
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			else
			{
				double *t_auxTrjRes = auxTrjRes + nsLeft*(numEq + 1) - 1;
				//for(int j=0; j<nsLeft; j++)
				for(long long j=0; j<nsLeft; j++)
				{
					//if(pOutBtzData) *(tOutBtzData++) = *(t_auxTrjRes--);
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtyData++) = bty;
					*(tOutYData++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtxData++) = btx;
					*(tOutXData++) = *(t_auxTrjRes--);
					t_auxTrjRes--;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
	}

	if(integOnRightIsNeeded)
	{
		double *tOutBtxData = pOutBtxData + is0, *tOutXData = pOutXData + is0;
		double *tOutBtyData = pOutBtyData + is0, *tOutYData = pOutYData + is0;
		double *tOutBtzData = pOutBtzData + is0, *tOutZData = pOutZData + is0;

		if(sMin > sAbsEdgeToler)
		{//need to "arrive" to sMin without storing trajectory data; step can be re-adjusted
			int auxNp = (int)fabs(sMin/sStep) + 1;
			if(auxNp <= 1) auxNp = 2;

			double *auxTrjRes = new double[(auxNp + ns)*(numEq + 1)];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			gmIntRK.solve(initCond, 0., sMin, auxNp, auxTrjRes);

			double *pResEnd = auxTrjRes + auxNp*(numEq + 1) - numEq;
			for(int i=0; i<5; i++) initCond[i] = pResEnd[i];

			//arrived to sMax; now solve for the entire main trajectory:
			if(ns <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				*(t_auxTrjRes++) = sMin;
				for(int i=0; i<numEq; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else gmIntRK.solve(initCond, sMin, sMax, ns, auxTrjRes);

			if(dataShouldBeRotatedAfterSolve)
			{
				double *t_auxTrjRes = auxTrjRes + ns*(numEq + 1) - 1;
				//for(int j=0; j<ns; j++)
				for(long long j=0; j<ns; j++)
				{
					//if(pOutBtzData) *(tOutBtzData++) = *(t_auxTrjRes--);
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtyData++) = bty;
					*(tOutYData++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtxData++) = btx;
					*(tOutXData++) = *(t_auxTrjRes--);
					t_auxTrjRes--;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			else
			{
				double *t_auxTrjRes = auxTrjRes;
				//for(int j=0; j<ns; j++)
				for(long long j=0; j<ns; j++)
				{
					t_auxTrjRes++;
					*(tOutXData++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtxData++) = btx;
					*(tOutYData++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtyData++) = bty;
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
		else
		{
			//normally, initial conditions should be already set here
			//int nsRight = ns - is0;
			long long nsRight = ns - is0;
			if(nsRight < 1) nsRight = 1;
			double *auxTrjRes = new double[nsRight*(numEq + 1)];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			if(nsRight <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				*(t_auxTrjRes++) = s0Act;
				for(int i=0; i<numEq; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else gmIntRK.solve(initCond, s0Act, sMax, nsRight, auxTrjRes);

			if(dataShouldBeRotatedAfterSolve)
			{
				double *t_auxTrjRes = auxTrjRes + nsRight*(numEq + 1) - 1;
				//for(int j=0; j<nsRight; j++)
				for(long long j=0; j<nsRight; j++)
				{
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtyData++) = bty;
					*(tOutYData++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtxData++) = btx;
					*(tOutXData++) = *(t_auxTrjRes--);
					t_auxTrjRes--;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			else
			{
				double *t_auxTrjRes = auxTrjRes;
				//for(int j=0; j<nsRight; j++)
				for(long long j=0; j<nsRight; j++)
				{
					t_auxTrjRes++;
					*(tOutXData++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtxData++) = btx;
					*(tOutYData++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtyData++) = bty;
					if(pOutZData) *(tOutZData++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtzData) *(tOutBtzData++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
	}
	bool BxIsReq = (pOutBxData != 0);
	bool ByIsReq = (pOutByData != 0);
	bool BzIsReq = (pOutBzData != 0);
	if(BxIsReq || ByIsReq || BzIsReq)
	{
		double *tOutXData = pOutXData, *tOutYData = pOutYData, *tOutZData = pOutZData;
		double *tOutBxData = pOutBxData, *tOutByData = pOutByData, *tOutBzData = pOutBzData;
		//for(int i=0; i<ns; i++)
		for(long long i=0; i<ns; i++)
		{
			TVector3d P(*(tOutXData++), *(tOutYData++), *(tOutZData++)), B;
			m_hMagElem.rep->compB(P, B);
			if(BxIsReq) *(tOutBxData++) = B.x;
			if(ByIsReq) *(tOutByData++) = B.y;
			if(BzIsReq) *(tOutBzData++) = B.z;
		}
	}
}

//*************************************************************************

//void srTGenTrjDat::CompTrjKickMatr(SRWLKickM* arKickM, int nKickM, double sStart, double sEnd, long ns, double* pInPrecPar, double* pOutBtX, double* pOutX, double* pOutBtY, double* pOutY, double* pOutBtZ, double* pOutZ)
void srTGenTrjDat::CompTrjKickMatr(SRWLKickM* arKickM, int nKickM, double sStart, double sEnd, long long ns, double* pInPrecPar, double* pOutBtX, double* pOutX, double* pOutBtY, double* pOutY, double* pOutBtZ, double* pOutZ)
{
	if((arKickM == 0) || (nKickM <= 0)) throw SRWL_INCORRECT_PARAM_FOR_TRJ_COMP;
	if(ns <= 0) throw SRWL_INCORRECT_PARAM_FOR_TRJ_COMP;

	//sort arKickM to find groups of kick-matrices with non-overlapping longitudinal intervals
	vector<pair<int, pair<double, double> > > vKickStInd;
	SRWLKickM *t_arKickM = arKickM;
	for(int i=0; i<nKickM; i++)
	{
		double curHalfRange = 0.5*t_arKickM->rz;
		pair<double, double> curRange(t_arKickM->z - curHalfRange, t_arKickM->z + curHalfRange);
		pair<int, pair<double, double> > curPair(i, curRange);
		vKickStInd.push_back(curPair);
	}

	sort(vKickStInd.begin(), vKickStInd.end(), CAuxParse::LessInPairBasedOnFirstInNestedPair<int, double, double>);

	vector<vector<int> > vIndNonOverlapKickGroups;
	vector<pair<double, double> > vIndNonOverlapKickGroupRanges;
	vector<int> vIndKickM;
	vIndKickM.push_back(vKickStInd[0].first);
	double curGroupStart = vKickStInd[0].second.first;
	double curGroupEnd = vKickStInd[0].second.second;
	for(int j=1; j<nKickM; j++)
	{
		double newStart = vKickStInd[j].second.first;
		double newEnd = vKickStInd[j].second.second;

		if((newEnd <= curGroupStart) || (newStart >= curGroupEnd))
		{//end current group
			vIndNonOverlapKickGroups.push_back(vIndKickM);
			pair<double, double> curRange(curGroupStart, curGroupEnd);
			vIndNonOverlapKickGroupRanges.push_back(curRange);

			vIndKickM.erase(vIndKickM.begin(), vIndKickM.end());
			curGroupStart = newStart;
			curGroupEnd = newEnd;
		}
		else
		{//continue stuffing current group
			vIndKickM.push_back(vKickStInd[j].first);
			if(curGroupStart > newStart) curGroupStart = newStart;
			if(curGroupEnd < newEnd) curGroupEnd = newEnd;
		}
	}
	if(!vIndKickM.empty()) 
	{
		vIndNonOverlapKickGroups.push_back(vIndKickM);
		pair<double, double> curRange(curGroupStart, curGroupEnd);
		vIndNonOverlapKickGroupRanges.push_back(curRange);

		vIndKickM.erase(vIndKickM.begin(), vIndKickM.end());
	}
	//sorting completed

	bool trjShouldBeAdded = (pInPrecPar[0] == 1);
	const double sResEdgeToler = 1.E-12;
	double sAbsEdgeToler = (sEnd - sStart)*sResEdgeToler;
	bool integOnLeftIsNeeded = (sStart < -sAbsEdgeToler);
	bool integOnRightIsNeeded = (sEnd > sAbsEdgeToler);

	double sStep = (ns <= 1)? 0 : (sEnd - sStart)/(ns - 1);
	double initCond[] = {EbmDat.x0, EbmDat.dxds0, EbmDat.z0, EbmDat.dzds0, EbmDat.s0}; //?
	double inv_B_pho = EbmDat.Inv_B_rho(); //[1/(T*m)]

	double gamEm2 = EbmDat.GammaEm2, btx, bty;
	//long is0 = 0;
	long long is0 = 0;
	double s0Act = 0.;

	if(integOnLeftIsNeeded)
	{
		double *tOutBtX = pOutBtX, *tOutX = pOutX;
		double *tOutBtY = pOutBtY, *tOutY = pOutY;
		double *tOutBtZ = pOutBtZ, *tOutZ = pOutZ;
		
		if(sEnd < -sAbsEdgeToler) //i.e. < 0
		{//need to "arrive" to sEnd without storing trajectory data; step can be re-adjusted
			long auxNp = (long)fabs(sEnd/sStep) + 1;
			if(auxNp <= 1) auxNp = 2;

			double *auxTrjRes = new double[(auxNp + ns)*5];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, 0., sEnd, auxNp, auxTrjRes);
			//gmIntRK.solve(initCond, 0., sMax, auxNp, auxTrjRes);

			double *pResEnd = auxTrjRes + (auxNp - 1)*5;
			for(int i=0; i<5; i++) initCond[i] = pResEnd[i];

			//arrived to sMax; now solve for the entire main trajectory:
			if(ns <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				for(int i=0; i<5; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, sEnd, sStart, ns, auxTrjRes);
			//else gmIntRK.solve(initCond, sMax, sMin, ns, auxTrjRes);

			double *t_auxTrjRes = auxTrjRes + (ns*5 - 1);
			//for(int j=0; j<ns; j++)
			for(long long j=0; j<ns; j++)
			{
				if(trjShouldBeAdded)
				{
					//if(pOutZ) *(tOutZ++) += *t_auxTrjRes; //longitudinal position is not modified in this case
					t_auxTrjRes--;
					*tOutBtY += *(t_auxTrjRes--); bty = *(tOutBtY++);
					*(tOutY++) += *(t_auxTrjRes--);
					*tOutBtX += *(t_auxTrjRes--); btx = *(tOutBtX++);
					*(tOutX++) += *(t_auxTrjRes--);
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
				else
				{
					if(pOutZ) *(tOutZ++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtY++) = bty;
					*(tOutY++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtX++) = btx;
					*(tOutX++) = *(t_auxTrjRes--);
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
		else
		{
			is0 = (int)fabs((-sStart + sAbsEdgeToler)/sStep);
			if(is0 >= ns) is0 = ns - 1;

			s0Act = sStart + is0*sStep;
			if(s0Act < -sAbsEdgeToler)
			{//one small step to s0Act
				double twoPtTrjRes[2*5];
				IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, 0., s0Act, 2, twoPtTrjRes);
				//gmIntRK.solve(initCond, 0., s0Act, 2, twoPtTrjRes);

				double *pResEnd = twoPtTrjRes + 5;
				for(int i=0; i<5; i++) initCond[i] = pResEnd[i];
			}

			//arrived to s0Act; now solve for the left part of the trajectory:
			//int nsLeft = is0 + 1;
			long long nsLeft = is0 + 1;
			double *auxTrjRes = new double[nsLeft*5];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;
	
			if(nsLeft <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				//*(t_auxTrjRes++) = s0Act;
				for(int i=0; i<5; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, s0Act, sStart, nsLeft, auxTrjRes);
			//else gmIntRK.solve(initCond, s0Act, sMin, nsLeft, auxTrjRes);

			double *t_auxTrjRes = auxTrjRes + nsLeft*5 - 1;
			//for(int j=0; j<nsLeft; j++)
			for(long long j=0; j<nsLeft; j++)
			{
				if(trjShouldBeAdded)
				{
					//if(pOutZ) *(tOutZ++) += *t_auxTrjRes; //longitudinal position is not modified in this case
					t_auxTrjRes--;
					*tOutBtY += *(t_auxTrjRes--); bty = *(tOutBtY++);
					*(tOutY++) += *(t_auxTrjRes--);
					*tOutBtX += *(t_auxTrjRes--); btx = *(tOutBtX++);
					*(tOutX++) += *(t_auxTrjRes--);
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
				else
				{
					if(pOutZ) *(tOutZ++) = *t_auxTrjRes;
					t_auxTrjRes--;
					bty = *(t_auxTrjRes--); *(tOutBtY++) = bty;
					*(tOutY++) = *(t_auxTrjRes--);
					btx = *(t_auxTrjRes--); *(tOutBtX++) = btx;
					*(tOutX++) = *(t_auxTrjRes--);
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
	}

	if(integOnRightIsNeeded)
	{
		double *tOutBtX = pOutBtX + is0, *tOutX = pOutX + is0;
		double *tOutBtY = pOutBtY + is0, *tOutY = pOutY + is0;
		double *tOutBtZ = pOutBtZ + is0, *tOutZ = pOutZ + is0;

		if(sStart > sAbsEdgeToler)
		{//need to "arrive" to sMin without storing trajectory data; step can be re-adjusted
			int auxNp = (int)fabs(sStart/sStep) + 1;
			if(auxNp <= 1) auxNp = 2;

			double *auxTrjRes = new double[(auxNp + ns)*5];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, 0., sStart, auxNp, auxTrjRes);
			//gmIntRK.solve(initCond, 0., sMin, auxNp, auxTrjRes);

			double *pResEnd = auxTrjRes + (auxNp - 1)*5;
			for(int i=0; i<5; i++) initCond[i] = pResEnd[i];

			//arrived to sMax; now solve for the entire main trajectory:
			if(ns <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				//*(t_auxTrjRes++) = sMin;
				for(int i=0; i<5; i++) *(t_auxTrjRes++) = initCond[i];
			}
			else IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, sStart, sEnd, ns, auxTrjRes);
			//else gmIntRK.solve(initCond, sMin, sMax, ns, auxTrjRes);

			double *t_auxTrjRes = auxTrjRes;
			//for(int j=0; j<ns; j++)
			for(long long j=0; j<ns; j++)
			{
				if(trjShouldBeAdded)
				{
					*(tOutX++) += *(t_auxTrjRes++);
					*tOutBtX += *(t_auxTrjRes++); btx = *(tOutBtX++);
					*(tOutY++) += *(t_auxTrjRes++);
					*tOutBtY += *(t_auxTrjRes++); bty = *(tOutBtY++);
					//if(pOutZ) *(tOutZ++) += *t_auxTrjRes; //longitudinal position is not modified in this case
					t_auxTrjRes++;
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
				else
				{
					*(tOutX++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtX++) = btx;
					*(tOutY++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtY++) = bty;
					if(pOutZ) *(tOutZ++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
		else
		{
			//normally, initial conditions should be already set here
			//int nsRight = ns - is0;
			long long nsRight = ns - is0;
			if(nsRight < 1) nsRight = 1;
			double *auxTrjRes = new double[nsRight*5];
			if(auxTrjRes == 0) throw MEMORY_ALLOCATION_FAILURE;

			if(nsRight <= 1)
			{
				double *t_auxTrjRes = auxTrjRes;
				//*(t_auxTrjRes++) = s0Act;
				for(int i=0; i<5; i++) *(t_auxTrjRes++) = initCond[i];
			}
			//s0Act has been defined above!
			else IntegrateKicks(arKickM, vIndNonOverlapKickGroups, vIndNonOverlapKickGroupRanges, inv_B_pho, initCond, s0Act, sEnd, nsRight, auxTrjRes);
			//else gmIntRK.solve(initCond, s0Act, sMax, nsRight, auxTrjRes);

			double *t_auxTrjRes = auxTrjRes;
			//for(int j=0; j<nsRight; j++)
			for(long long j=0; j<nsRight; j++)
			{
				if(trjShouldBeAdded)
				{
					*(tOutX++) += *(t_auxTrjRes++);
					*tOutBtX += *(t_auxTrjRes++); btx = *(tOutBtX++);
					*(tOutY++) += *(t_auxTrjRes++);
					*tOutBtY += *(t_auxTrjRes++); bty = *(tOutBtY++);
					//if(pOutZ) *(tOutZ++) += *t_auxTrjRes; //longitudinal position is not modified in this case
					t_auxTrjRes++;
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
				else
				{
					*(tOutX++) = *(t_auxTrjRes++);
					btx = *(t_auxTrjRes++); *(tOutBtX++) = btx;
					*(tOutY++) = *(t_auxTrjRes++);
					bty = *(t_auxTrjRes++); *(tOutBtY++) = bty;
					if(pOutZ) *(tOutZ++) = *t_auxTrjRes;
					t_auxTrjRes++;
					if(pOutBtZ) *(tOutBtZ++) = CGenMathMeth::radicalOnePlusSmall(-(gamEm2 + btx*btx + bty*bty));
				}
			}
			delete[] auxTrjRes;
		}
	}
}

//*************************************************************************
// Performs integration of trajectory based on kicks
// initCond[5] - array of initial conditions, defined at s = sStart
// sStart - initial argument
// sEnd - final argument
// ns - number of points
// pTrjRes - resulting flat trajectory array, length is equal to 5*ns
//-------------------------------------------------------------------------
//void srTGenTrjDat::IntegrateKicks(SRWLKickM* arKickM, vector<vector<int> >& vIndNonOverlapKickGroups, vector<pair<double, double> >& vIndNonOverlapKickGroupRanges, double inv_B_pho, double* initCond, double sStart, double sEnd, int ns, double* pTrjRes)
void srTGenTrjDat::IntegrateKicks(SRWLKickM* arKickM, vector<vector<int> >& vIndNonOverlapKickGroups, vector<pair<double, double> >& vIndNonOverlapKickGroupRanges, double inv_B_pho, double* initCond, double sStart, double sEnd, long long ns, double* pTrjRes)
{
	int nGroups = (int)vIndNonOverlapKickGroups.size();
	if((arKickM == 0) || (nGroups <= 0) || (initCond == 0) || (pTrjRes == 0)) return;

	double *tX = pTrjRes, *tBtX = pTrjRes + 1, *tY = pTrjRes + 2, *tBtY = pTrjRes + 3, *tZ = pTrjRes + 4;
	*tX = initCond[0]; *tBtX = initCond[1]; *tY = initCond[2]; *tBtY = initCond[3]; *tZ = initCond[4];

	map<int, pair<double, pair<double, double> > >  mAuxPrevKick; //to keep previous kick values and longitudinal positions at which these values were obtained by interpolation
	double arF[12];
	double sStep = (ns > 1)? (sEnd - sStart)/(ns - 1) : 0.;
	//s is longitudinal position from here on!
	double s = initCond[4]; //sStart; // + sStep;
	//for(int is=0; is<(ns-1); is++)
	for(long long is=0; is<(ns-1); is++)
	{
		//double Xprev = *tX, BtXprev = *tBtX, Yprev = *tY, BtYprev = *tBtY, Zprev = *tZ;
		double dX = 0., dBtX = 0., dY = 0., dBtY = 0., dZ = 0.;

		for(int i=0; i<nGroups; i++)
		{//"non-overlapping groups"
			double curGroupStartS = vIndNonOverlapKickGroupRanges[i].first;
			double curGroupEndS = vIndNonOverlapKickGroupRanges[i].second;
			if((curGroupStartS <= s) && (s < curGroupEndS))
			{
				vector<int>& curVectIndKicks = vIndNonOverlapKickGroups[i];
				int nKickInGroup = (int)curVectIndKicks.size();

				double kx_ds = 0., ky_ds = 0., kx = 0., ky = 0.;
				for(int j=0; j<nKickInGroup; j++)
				{//sum-up current kicks from this group
					int indCurKick = curVectIndKicks[j];
					SRWLKickM *pCurKickM = arKickM + indCurKick;

					double sHalfRange = 0.5*pCurKickM->rz;
					double sStartCurKick = pCurKickM->z - sHalfRange;
					double sEndCurKick = pCurKickM->z + sHalfRange;

					double dsTest = s + sStep - sStartCurKick;
					//if((sStartCurKick < (s + sStep)) && (s < sEndCurKick))
					if((0. < dsTest) && (s < sEndCurKick))
					{
						double dsKick = sStep;
						if(dsKick > dsTest) dsKick = dsTest;

						double yHalfRange = 0.5*pCurKickM->ry;
						double yStartCurKick = pCurKickM->y - yHalfRange;
						double yEndCurKick = pCurKickM->y + yHalfRange;
						//if((yStartCurKick <= *tY) && (*tY < yEndCurKick))
						//{
							int iy0 = 0;
							double yStep = 0., yt = 0.;
							int Ny = pCurKickM->ny;
							if(Ny > 1) 
							{
								yStep = pCurKickM->ry/(Ny - 1);
								iy0 = (int)((*tY - yStartCurKick)/yStep + 1.e-09);
								yt = (*tY - (yStartCurKick + iy0*yStep))/yStep;
								if(iy0 < 0) 
								{
									iy0 = 0; yt = 0.;
								}
								else if(iy0 >= Ny) 
								{
									iy0 = Ny - 1; yt = 0.;
								}
							}

							double xHalfRange = 0.5*pCurKickM->rx;
							double xStartCurKick = pCurKickM->x - xHalfRange;
							double xEndCurKick = pCurKickM->x + xHalfRange;
							//if((xStartCurKick <= *tX) && (*tX < xEndCurKick))
							//{
								double mult = inv_B_pho;
								if(pCurKickM->order == 2) mult *= mult;

								mult *= dsKick/pCurKickM->rz;

								bool calcNewKickVals = false;
								map<int, pair<double, pair<double, double> > >::iterator itPrevKick = mAuxPrevKick.find(indCurKick);
								if(itPrevKick != mAuxPrevKick.end())
								{//to save time: don't interpolate (use previous interpolated kick value) if the current kick step is not exceeded
									pair<double, pair<double, double> > &pairKickInf = itPrevKick->second;
									double sPrevKick = pairKickInf.first;
									double sStepKick = (pCurKickM->nz > 1)? pCurKickM->rz/(pCurKickM->nz - 1) : 0.;
									if((s - sPrevKick) <= sStepKick)
									{
										if(pCurKickM->arKickMx != 0) 
										{
											double dkx = mult*pairKickInf.second.first; 
											kx += dkx;
											kx_ds += dkx*dsKick; 
										}
										if(pCurKickM->arKickMy != 0) 
										{
											double dky = mult*pairKickInf.second.second; 
											ky += dky;
											ky_ds += dky*dsKick; 
										}
									}
									else calcNewKickVals = true;
								}
								else calcNewKickVals = true;

								if(calcNewKickVals)
								{
									int ix0 = 0;
									double xStep = 0., xt = 0.;
									int Nx = pCurKickM->nx;
									if(Nx > 1) 
									{
										xStep = pCurKickM->rx/(Nx - 1);
										ix0 = (int)((*tX - xStartCurKick)/xStep + 1.e-09);
										xt = (*tX - (xStartCurKick + ix0*xStep))/xStep;
										if(ix0 < 0) 
										{
											ix0 = 0; xt = 0.;
										}
										else if(ix0 >= Nx) 
										{
											ix0 = Nx - 1; xt = 0.;
										}
									}
									//int ixm1 = ix0 - 1; if(ixm1 < 0) ixm1 = 0;
									//int ix1 = ix0 + 1; if(ix1 >= Nx) ix1 = Nx - 1;
									//int ix2 = ix1 + 1; if(ix2 >= Nx) ix2 = Nx - 1;
									//int iym1 = iy0 - 1; if(iym1 < 0) iym1 = 0;
									//int iy1 = iy0 + 1; if(iy1 >= Ny) iy1 = Ny - 1;
									//int iy2 = iy1 + 1; if(iy2 >= Ny) iy2 = Ny - 1;
									long long ixm1 = ix0 - 1; if(ixm1 < 0) ixm1 = 0;
									long long ix1 = ix0 + 1; if(ix1 >= Nx) ix1 = Nx - 1;
									long long ix2 = ix1 + 1; if(ix2 >= Nx) ix2 = Nx - 1;
									long long iym1 = iy0 - 1; if(iym1 < 0) iym1 = 0;
									long long iy1 = iy0 + 1; if(iy1 >= Ny) iy1 = Ny - 1;
									long long iy2 = iy1 + 1; if(iy2 >= Ny) iy2 = Ny - 1;

									//find current contributions to kx, ky by interpolation, taking into account kick order
									//long Nx_iym1 = Nx*iym1;
									//long i0m1 = ix0 + Nx_iym1; //f0m1
									//long i1m1 = ix1 + Nx_iym1; //f1m1
									//long Nx_iy0 = Nx*iy0;
									//long im10 = ixm1 + Nx_iy0; //fm10
									//long i00 = ix0 + Nx_iy0; //f00
									//long i10 = ix1 + Nx_iy0; //f10
									//long i20 = ix2 + Nx_iy0; //f20
									//long Nx_iy1 = Nx*iy1;
									//long im11 = ixm1 + Nx_iy1; //fm11
									//long i01 = ix0 + Nx_iy1; //f01
									//long i11 = ix1 + Nx_iy1; //f11
									//long i21 = ix2 + Nx_iy1; //f21
									//long Nx_iy2 = Nx*iy2;
									//long i02 = ix0 + Nx_iy2; //f02
									//long i12 = ix1 + Nx_iy2; //f12
									long long Nx_iym1 = Nx*iym1;
									long long i0m1 = ix0 + Nx_iym1; //f0m1
									long long i1m1 = ix1 + Nx_iym1; //f1m1
									long long Nx_iy0 = Nx*iy0;
									long long im10 = ixm1 + Nx_iy0; //fm10
									long long i00 = ix0 + Nx_iy0; //f00
									long long i10 = ix1 + Nx_iy0; //f10
									long long i20 = ix2 + Nx_iy0; //f20
									long long Nx_iy1 = Nx*iy1;
									long long im11 = ixm1 + Nx_iy1; //fm11
									long long i01 = ix0 + Nx_iy1; //f01
									long long i11 = ix1 + Nx_iy1; //f11
									long long i21 = ix2 + Nx_iy1; //f21
									long long Nx_iy2 = Nx*iy2;
									long long i02 = ix0 + Nx_iy2; //f02
									long long i12 = ix1 + Nx_iy2; //f12

									double newKickX = 0., newKickY = 0.;
									double *pM = pCurKickM->arKickMx;
									if(pM != 0)
									{
										double *tF = arF;
										*(tF++) = *(pM + i0m1); //f0m1
										*(tF++) = *(pM + i1m1); //f1m1
										*(tF++) = *(pM + im10); //fm10
										*(tF++) = *(pM + i00); //f00
										*(tF++) = *(pM + i10); //f10
										*(tF++) = *(pM + i20); //f20
										*(tF++) = *(pM + im11); //fm11
										*(tF++) = *(pM + i01); //f01
										*(tF++) = *(pM + i11); //f11
										*(tF++) = *(pM + i21); //f21
										*(tF++) = *(pM + i02); //f02
										*(tF++) = *(pM + i12); //f12
										newKickX = CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arF);
										double dkx = mult*newKickX;
										kx += dkx;
										kx_ds += dkx*dsKick;
									}
									pM = pCurKickM->arKickMy;
									if(pM != 0)
									{
										double *tF = arF;
										*(tF++) = *(pM + i0m1); //f0m1
										*(tF++) = *(pM + i1m1); //f1m1
										*(tF++) = *(pM + im10); //fm10
										*(tF++) = *(pM + i00); //f00
										*(tF++) = *(pM + i10); //f10
										*(tF++) = *(pM + i20); //f20
										*(tF++) = *(pM + im11); //fm11
										*(tF++) = *(pM + i01); //f01
										*(tF++) = *(pM + i11); //f11
										*(tF++) = *(pM + i21); //f21
										*(tF++) = *(pM + i02); //f02
										*(tF++) = *(pM + i12); //f12
										newKickY = CGenMathInterp::Interp2dBiCubic12pRel(xt, yt, arF);
										double dky = mult*newKickY;
										ky += dky;
										ky_ds += dky*dsKick;
									}
									pair<double, pair<double, double> > pairNewKickInf(s, pair<double, double>(newKickX, newKickY));
									mAuxPrevKick[indCurKick] = pairNewKickInf;
								}
							//}
						//}
					}
					else if(s >= sEndCurKick)
					{
						mAuxPrevKick.erase(indCurKick); //to minimize search in this lookup table
					}
				}
				//to add kicks to trajectory here:
				//Xprev += kx_ds; BtXprev += kx; Yprev += ky_ds; BtYprev += ky;
				dX += kx_ds; dBtX += kx; dY += ky_ds; dBtY += ky;
			}
			else if(s < curGroupStartS) break;
		}
		double Xprev = *tX, BtXprev = *tBtX, Yprev = *tY, BtYprev = *tBtY, Zprev = *tZ;
		tX += 5; tBtX += 5; tY += 5; tBtY += 5; tZ += 5;
		//*tX = Xprev; *tBtX = BtXprev; *tY = Yprev; *tBtY = BtYprev; *tZ = Zprev + sStep;
		*tX = Xprev + dX + BtXprev*sStep; 
		*tBtX = BtXprev + dBtX; 
		*tY = Yprev + dY + BtYprev*sStep; 
		*tBtY = BtYprev + dBtY; 
		*tZ = Zprev + sStep;

		s += sStep;
	}
}

//*************************************************************************
