/************************************************************************//**
 * File: srtrjdat3d.cpp
 * Description: Electron Trajectory (and relevant characteristics) calculation in 3D Magnetic Fields
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srtrjdat3d.h"
#include "gmintrk.h"

//*************************************************************************
//Not Used? (RK moved to srTGenTrjDat)

//void srTTrjDat3d::CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, int ns, double sStart, double sStep)
void srTTrjDat3d::CompTrjDataForDisp(double* pOutBtxData, double* pOutXData, double* pOutBtyData, double* pOutYData, double* pOutBtzData, double* pOutZData, long long ns, double sStart, double sStep)
{//Corrected for SRWL:
 //Independent variable is s = c*t
 //3D trajectory is calculated in Laboratory Frame
 //Z is assumed to be longitudinal variable here!
	if(ns <= 0) return;
	//bool x1_is_needed = (pOutBtxData != 0);
	//bool x_is_needed = (pOutXData != 0);
	//bool z1_is_needed = (pOutBtzData != 0);
	//bool z_is_needed = (pOutZData != 0);

	const double chElec = 1.602176462e-19; //[C]
	const double mElec = 9.10938188e-31; //[kg]
	const double cLight = 2.99792458e+08; //[m/s]
	//double gamElec = EbmDat.Gamma;
	//m_Mult2ndDer = -chElec/(mElec*cLight*EbmDat.Gamma*sqrt(1. - 1./(EbmDat.Gamma*EbmDat.Gamma)));
	m_Mult2ndDer = EbmDat.nQ*chElec/(mElec*cLight*EbmDat.Gamma);

	double sEnd = sStart + sStep*(ns - 1);
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
	CGenMathIntRungeKutta<srTTrjDat3d> gmIntRK(this, &srTTrjDat3d::funcDerivF, numEq, false, 0); // manual mode, based on number of points
	double initCond[] = {EbmDat.x0, EbmDat.dxds0, EbmDat.z0, EbmDat.dzds0, EbmDat.sc};
	
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
			//is0 = (int)fabs((-sMin + sAbsEdgeToler)/sStep);
			is0 = (long long)fabs((-sMin + sAbsEdgeToler)/sStep);
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
				double *t_auxTrjRes = auxTrjRes + (nsLeft*(numEq + 1) - 1);
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
				double *t_auxTrjRes = auxTrjRes + (nsRight*(numEq + 1) - 1);
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

/**

	int is0 = (int)((s0_mi_sStart + sAbsEdgeToler)/sStep);
	if(is0 >= ns) is0 = ns - 1;

		//double test = sStart + is0*sStep;

	double dsRight = 0.;
	double dsLeft = EbmDat.s0 - (sStart + is0*sStep);
	if(dsLeft > sAbsEdgeToler) 
	{
		dsRight = sStep - dsLeft;
	}
	else
	{
		dsLeft = 0; dsRight = 0.;
	}

	//if(is0 < (ns - 1))
	//{
	//	dsRight = sStart + (is0 + 1)*sStep - EbmDat.s0;
	//	if(dsRight < sAbsEdgeToler) dsRight = 0;
	//}

	const int numEq = 6;
	CGenMathIntRungeKutta<srTTrjDat3d> gmIntRK(this, &srTTrjDat3d::funcDerivF, numEq, false, 0); // manual mode, based on number of points
	//CGenMathIntRungeKutta(T* ptrT, void (T::*pFuncDerivF)(double, double*, double*), int amOfEq, bool onPrc, double* precArray, double epsTol =1., int maxAutoStp =5000)

	double y0 = EbmDat.s0, dyds0 = 1.; //to check
	double initCond[] = {EbmDat.x0, EbmDat.dxds0, y0, dyds0, EbmDat.z0, EbmDat.dzds0};

	if(integOnRightIsNeeded)
	{
		double initCondRight[numEq];
		int isStartRight = is0;
		double sStartRight = EbmDat.s0;
		for(int i=0; i<numEq; i++) initCondRight[i] = initCond[i];

		if(dsRight > 0.)
		{
            double auxResRight[2*(numEq + 1)];
			gmIntRK.solve(initCond, sStartRight, sStartRight + dsRight, 2, auxResRight);
			sStartRight += dsRight;
			double *tAuxResRight = auxResRight + numEq + 2, *tInitCondRight = initCondRight;
            for(int j=0; j<numEq; j++) *(tInitCondRight++) = *(tAuxResRight++);
			isStartRight++;
		}

		int npRight = ns - isStartRight;
		double *auxTotResRight = new double[npRight*(numEq + 1)];
        if(auxTotResRight == 0) throw MEMORY_ALLOCATION_FAILURE;

		gmIntRK.solve(initCondRight, sStartRight, sEnd, npRight, auxTotResRight);

		long ofst = isStartRight*(numEq + 1);
        double *tOutBtxData = pOutBtxData + ofst, *tOutXData = pOutXData + ofst;
		double *tOutBtyData = pOutBtyData + ofst, *tOutYData = pOutYData + ofst;
		double *tOutBtzData = pOutBtzData + ofst, *tOutZData = pOutZData + ofst;
		double *tAuxTotResRight = auxTotResRight;
		for(int j=0; j<npRight; j++)
		{
			tAuxTotResRight++; //may need to be put into another place (to check)
            *(tOutXData++) = *(tAuxTotResRight++);
            *(tOutBtxData++) = *(tAuxTotResRight++);
            *(tOutYData++) = *(tAuxTotResRight++);
            *(tOutBtyData++) = *(tAuxTotResRight++);
            *(tOutZData++) = *(tAuxTotResRight++);
            *(tOutBtzData++) = *(tAuxTotResRight++);
		}
		if(auxTotResRight != 0) delete[] auxTotResRight;
	}

	if(integOnLeftIsNeeded)
	{
		double initCondLeft[numEq];
		int isStartLeft = is0;
		double sStartLeft = EbmDat.s0;
		for(int i=0; i<numEq; i++) initCondLeft[i] = initCond[i];

		if(dsLeft > 0.)
		{
            double auxResLeft[2*(numEq + 1)];
			gmIntRK.solve(initCond, sStartLeft, sStartLeft - dsLeft, 2, auxResLeft);
			sStartLeft -= dsLeft;
			double *tAuxResLeft = auxResLeft + numEq + 2, *tInitCondLeft = initCondLeft;
			for(int j=0; j<numEq; j++) *(tInitCondLeft++) = *(tAuxResLeft++);
			isStartLeft--;
		}

		int npLeft = isStartLeft + 1;
		double *auxTotResLeft = new double[npLeft*(numEq + 1)];
        if(auxTotResLeft == 0) throw MEMORY_ALLOCATION_FAILURE;

		gmIntRK.solve(initCondLeft, sStartLeft, sStart, npLeft, auxTotResLeft);

		long ofst = isStartLeft*(numEq + 1);
        double *tOutBtxData = pOutBtxData + ofst, *tOutXData = pOutXData + ofst;
		double *tOutBtyData = pOutBtyData + ofst, *tOutYData = pOutYData + ofst;
		double *tOutBtzData = pOutBtzData + ofst, *tOutZData = pOutZData + ofst;
		double *tAuxTotResLeft = auxTotResLeft;
		for(int j=0; j<npLeft; j++)
		{
			tAuxTotResLeft++; //may	need to	be put into	another	place (to check)
			*(tOutXData--) = *(tAuxTotResLeft++);
			*(tOutBtxData--) = *(tAuxTotResLeft++);
			*(tOutYData--) = *(tAuxTotResLeft++);
			*(tOutBtyData--) = *(tAuxTotResLeft++);
			*(tOutZData--) = *(tAuxTotResLeft++);
			*(tOutBtzData--) = *(tAuxTotResLeft++);
		}
		if(auxTotResLeft != 0) delete[] auxTotResLeft;
	}
**/
}

//*************************************************************************
