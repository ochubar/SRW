/************************************************************************//**
 * File: gmintrk.h
 * Description: Template class for Runge-Kutta integration (mostly taken from "Numerical Recipes in C")
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMINTRK_H
#define __GMINTRK_H

#ifndef __GMERCODE_H
#include "gmercode.h"
#endif

//-------------------------------------------------------------------------
// Integration using Runge-Kutta method(s)
//
// The following integer constants can be thrown 
// and therefore must be defined:
// GM_RK_MAX_NUM_STEPS_REACHED
//-------------------------------------------------------------------------

template<class T> class CGenMathIntRungeKutta {

	int m_AmOfEq, m_MaxAutoStp;
	bool m_OnPrc;
	double m_EpsTol;

	double *m_dym_rk4, *m_dyt_rk4, *m_yt_rk4;
	double *m_Y, *m_dYdx;
	double *m_PrecArr;
	double *m_dysav_rks5, *m_ysav_rks5, *m_ytemp_rks5;
	double *m_y_ap, *m_dydx_ap;
	double *m_xp_ap, **m_yp_ap, m_dxsav_ap; 
	int m_count_ap, m_kmax_ap;

	T* m_PtrT;
	void (T::*m_pFuncDerivF)(double in_x, double* in_F, double* out_dFdx);  

public:

	CGenMathIntRungeKutta(T* ptrT, void (T::*pFuncDerivF)(double, double*, double*), int amOfEq, bool onPrc, double* precArray, double epsTol =1., int maxAutoStp =5000)
	{
		m_PtrT = ptrT;
		m_pFuncDerivF = pFuncDerivF;
		m_AmOfEq = amOfEq;
		m_OnPrc = onPrc; 

		m_dym_rk4 = new double[amOfEq];
		m_dyt_rk4 = new double[amOfEq];
		m_yt_rk4 = new double[amOfEq];
        m_Y = new double[amOfEq];
        m_dYdx = new double[amOfEq];

		m_PrecArr = 0;
		m_dysav_rks5 = m_ysav_rks5 = m_ytemp_rks5 = 0;
		m_y_ap = m_dydx_ap = 0;
		if(onPrc)
		{
            m_PrecArr = new double[amOfEq];
			for(int i=0; i<amOfEq; i++) 
			{
                m_PrecArr[i] = precArray[i];
                if(m_PrecArr[i] == 0.) m_PrecArr[i] = 1.E+23;
			}
			m_EpsTol = epsTol;
			m_MaxAutoStp = maxAutoStp;

			m_dysav_rks5 = new double[amOfEq];
			m_ysav_rks5	= new double[amOfEq];
			m_ytemp_rks5 = new double[amOfEq];
			m_y_ap = new double[amOfEq];
			m_dydx_ap =	new	double[amOfEq];

            m_count_ap = 0;
            m_kmax_ap = 0;								// To allow storage of intermediate results in AutoPropagate,
			m_xp_ap = 0; m_yp_ap = 0; m_dxsav_ap = 0.;	// set this to desired values and allocate memory
		}
	}
	~CGenMathIntRungeKutta()
	{
		if(m_dym_rk4 != 0) delete[] m_dym_rk4;
		if(m_dyt_rk4 != 0) delete[] m_dyt_rk4;
		if(m_yt_rk4 != 0) delete[] m_yt_rk4;
		if(m_Y != 0) delete[] m_Y; 
		if(m_dYdx != 0) delete[] m_dYdx;

		if(m_OnPrc)
		{
			if(m_PrecArr != 0) delete[] m_PrecArr;
			if(m_dysav_rks5 != 0) delete[] m_dysav_rks5;
			if(m_ysav_rks5 != 0) delete[] m_ysav_rks5;
			if(m_ytemp_rks5 != 0) delete[] m_ytemp_rks5;
			if(m_y_ap != 0) delete[] m_y_ap;
			if(m_dydx_ap != 0) delete[] m_dydx_ap;
		}
	}

	//void solve(double* initCond, double xmin, double xmax, int np, double* res);
	void solve(double* initCond, double xmin, double xmax, long long np, double* res);
	void stepRungeKutta4(double* y, double* dydx, double x, double h);
	void stepRungeKutta5(double* y, double* dydx, double* x, double htry, double* hdid, double* hnext);
	void autoPropagate(double* ystart, double x1, double x2, double h1, double hmin, int* nok, int* nbad);


	//static bool VectCheckIfCollinear(double xV1, double yV1, double zV1, double xV2, double yV2, double zV2, double RelPrec);

	//static double Integ1D_FuncWithEdgeDer(double (*pF)(double), double x1, double x2, double dFdx1, double dFdx2, double RelPrec);
	//static double Integ1D_FuncDefByArray(double* FuncArr, long Np, double Step);
	//static double Integ1D_FuncDefByArray(float* FuncArr, long Np, double Step);
};

//-------------------------------------------------------------------------
// Performs integration of ordinary differential equations using Runge-Kutta 
// method of the 4th or 5th order
// initCond - array of initial conditions (defined at xmin), length is equal to m_AmOfEq
// xmin - initial argument
// xmax - final argument
// np - number of points
// resArr - resulting flat array, length is equal to (m_AmOfEq + 1)*np
//-------------------------------------------------------------------------
//template <class T> void CGenMathIntRungeKutta<T>::solve(double* initCond, double xmin, double xmax, int np, double* resArr)
template <class T> void CGenMathIntRungeKutta<T>::solve(double* initCond, double xmin, double xmax, long long np, double* resArr)
{
	double step_x = (xmax - xmin)/double(np - 1);
	double x = xmin;

	double minStepAllowed = (1.E-12) * step_x;
	int nok = 0, nbad = 0;

	for(int k=0; k<m_AmOfEq; k++) m_Y[k] = initCond[k];

	int amOfEq_p_1 = m_AmOfEq + 1;
	//int np_mi_1 = np - 1;
	long long np_mi_1 = np - 1;

	//for(int i=0; i<np; i++)
	for(long long i=0; i<np; i++)
	{
		if(!m_OnPrc) (m_PtrT->*m_pFuncDerivF)(x, m_Y, m_dYdx);

		//int i_amOfEq_p_1 = i*amOfEq_p_1; //, i_amOfEq_p_1_p_1 = i_amOfEq_p_1 + 1;
		long long i_amOfEq_p_1 = i*amOfEq_p_1; //, i_amOfEq_p_1_p_1 = i_amOfEq_p_1 + 1;
		
		double *t_res =  resArr + i_amOfEq_p_1;
		*(t_res++) = x;
		double *tY = m_Y;
		for(int k=0; k<m_AmOfEq; k++) *(t_res++) = *(tY++);

		if(i != np_mi_1)
		{
			double x1 = x + step_x;
			if(m_OnPrc) autoPropagate(m_Y, x, x1, step_x, minStepAllowed, &nok, &nbad);
			else stepRungeKutta4(m_Y, m_dYdx, x, step_x);
			//else RungeKutta4(Y, dYdx, x, Step_x, Y);
			x = x1;
		}
	}
}

//-------------------------------------------------------------------------

template <class T> void CGenMathIntRungeKutta<T>::stepRungeKutta4(double* y, double* dydx, double x, double h)
{
	double xh, hh = 0.5*h, h6 = h/6.;
	xh = x + hh;

	double *t_yt_rk4 = m_yt_rk4, *t_y = y, *t_dydx = dydx;
	int i;
	for(i=0; i<m_AmOfEq; i++) *(t_yt_rk4++) = *(t_y++) + hh*(*(t_dydx++));

	(m_PtrT->*m_pFuncDerivF)(xh, m_yt_rk4, m_dyt_rk4);

	t_yt_rk4 = m_yt_rk4; t_y = y;
	double *t_dyt_rk4 = m_dyt_rk4;
	//for(i=0; i<m_AmOfEq; i++) *(t_dyt_rk4++) = *(t_y++) + hh*(*(t_dyt_rk4++));
	for(i=0; i<m_AmOfEq; i++) *(t_yt_rk4++) = *(t_y++) + hh*(*(t_dyt_rk4++)); //OC281110

	(m_PtrT->*m_pFuncDerivF)(xh, m_yt_rk4, m_dym_rk4);

    t_yt_rk4 = m_yt_rk4; t_y = y; t_dyt_rk4 = m_dyt_rk4;
	double *t_dym_rk4 = m_dym_rk4;
	for(i=0; i<m_AmOfEq; i++)
	{
		*(t_yt_rk4++) = *(t_y++) + h*(*t_dym_rk4);
		*(t_dym_rk4++) += *(t_dyt_rk4++);
	}

	(m_PtrT->*m_pFuncDerivF)(x+h, m_yt_rk4, m_dyt_rk4);
	
	t_y = y; t_dydx = dydx; t_dyt_rk4 = m_dyt_rk4; t_dym_rk4 = m_dym_rk4;
	for(i=0; i<m_AmOfEq; i++) *(t_y++) += h6*((*(t_dydx++)) + (*(t_dyt_rk4++)) + 2.*(*(t_dym_rk4++)));
	//for(i=0; i<AmOfEq; i++) yout[i] = y[i] + h6*(dydx[i]+dyt_rk4[i]+2.*dym_rk4[i]);
}

//-------------------------------------------------------------------------

template <class T> void CGenMathIntRungeKutta<T>::stepRungeKutta5(double* y, double* dydx, double* x, double htry, double* hdid, double* hnext)
{
	const double PGROW = -0.2;
	const double PSHRNK = -0.25;
	const double FCOR = 1./15.;
	const double SAFETY = 0.9;
	const double ERRCON = 6.0E-04;

	int i;
	double xsav = *x, hh, h = htry, temp, errmax;

	for(i=0; i<m_AmOfEq; i++)
	{
		m_ysav_rks5[i] = y[i];
		m_dysav_rks5[i] = dydx[i];
	}

	for(;;)
	{
		hh = 0.5*h;

		double *t_ytemp_rks5 = m_ytemp_rks5, *t_ysav_rks5 = m_ysav_rks5;
        for(i=0; i<m_AmOfEq; i++) *(t_ytemp_rks5++) = *(t_ysav_rks5++);
		stepRungeKutta4(m_ytemp_rks5, m_dysav_rks5, xsav, hh);
		*x = xsav + hh;
		(m_PtrT->*m_pFuncDerivF)(*x, m_ytemp_rks5, dydx);

		t_ytemp_rks5 = m_ytemp_rks5;
		double *t_y = y;
        for(i=0; i<m_AmOfEq; i++) *(t_y++) = *(t_ytemp_rks5++);
		stepRungeKutta4(y, dydx, *x, hh);
		*x = xsav + h;
		if(*x == xsav) throw GM_RK_STEP_SIZE_TOO_SMALL;
		//if(*x == xsav) { Send.ErrorMessage("Radia::Error200"); return 0;}

		t_ysav_rks5 = m_ysav_rks5; t_ytemp_rks5 = m_ytemp_rks5;
        for(i=0; i<m_AmOfEq; i++) *(t_ytemp_rks5++) = *(t_ysav_rks5++);
		stepRungeKutta4(m_ytemp_rks5, m_dysav_rks5, xsav, h);

		errmax = 0.;
		for(i=0; i<m_AmOfEq; i++)
		{
			m_ytemp_rks5[i] = y[i] - m_ytemp_rks5[i];
			temp = fabs(m_ytemp_rks5[i]/m_PrecArr[i]);
			if(errmax < temp) errmax = temp;
		}

		errmax /= m_EpsTol;
		if(errmax <= 1.)
		{
			*hdid = h;
			*hnext = (errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4.*h);
			break;
		}
		h *= SAFETY*exp(PSHRNK*log(errmax));
	}
	for(i=0; i<m_AmOfEq; i++) y[i] += m_ytemp_rks5[i]*FCOR;
	//return 1;
}

//-------------------------------------------------------------------------

template <class T> void CGenMathIntRungeKutta<T>::autoPropagate(double* ystart, double x1, double x2, double h1, double hmin, int* nok, int* nbad)
{
	int nstp, i;
	double xsav, x = x1, hnext, hdid;
	double h = (x2 > x1)? fabs(h1) : -fabs(h1);

	*nok = *nbad = m_count_ap = 0;

	for(i=0; i<m_AmOfEq; i++) m_y_ap[i] = ystart[i];

	if(m_kmax_ap > 0) xsav = x - 2.*m_dxsav_ap;

	for(nstp=1; nstp<=m_MaxAutoStp; nstp++)
	{
		(m_PtrT->*m_pFuncDerivF)(x, m_y_ap, m_dydx_ap);

		if(m_kmax_ap > 0)
		{
			if(fabs(x - xsav) > fabs(m_dxsav_ap))
			{
				if(m_count_ap < (m_kmax_ap - 1))
				{
					m_xp_ap[++m_count_ap] = x;
					for(i=0; i<m_AmOfEq; i++) m_yp_ap[i][m_count_ap] = m_y_ap[i];
					xsav = x;
				}
			}
		}

		if((x + h - x2)*(x + h - x1) > 0.) h = x2 - x;

		stepRungeKutta5(m_y_ap, m_dydx_ap, &x, h, &hdid, &hnext);
		
		if(hdid == h) ++(*nok);
		else ++(*nbad);

		if((x-x2)*(x2-x1) >= 0.)
		{
			for(i=0; i<m_AmOfEq; i++) ystart[i] = m_y_ap[i];
			if(m_kmax_ap > 0)
			{
                m_xp_ap[++m_count_ap] = x;
				for(i=0; i<m_AmOfEq; i++) m_yp_ap[i][m_count_ap] = m_y_ap[i];
			}
		}
		if(fabs(hnext) <= hmin) throw GM_RK_STEP_SIZE_TOO_SMALL;
		//if(fabs(hnext) <= hmin) { Send.ErrorMessage("Radia::Error200"); return 0;}
		h = hnext;
	}
	throw GM_RK_MAX_NUM_STEPS_REACHED;
	//Send.ErrorMessage("Radia::Error201"); return 0;
}

//-------------------------------------------------------------------------

#endif


