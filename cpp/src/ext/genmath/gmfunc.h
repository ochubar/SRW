/************************************************************************//**
 * File: gmfunc.h
 * Description: Some mathematical (special) functions header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author P.Elleaume, O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __GMFUNC_H
#define __GMFUNC_H

#ifndef _GM_WITHOUT_BASE
#include "gmobj.h"
#endif

//*************************************************************************

//class srTMathFunctions {
class CGenMathFunc
#ifndef _GM_WITHOUT_BASE
	: public CGenMathObj
#endif
{
public:

	/** Computes the mofified Bessel Function Kn(x)
	@see Method from V.O. Kostrum, NIM 172 (1980) p371-374
    ni=0 : computes k= Kn(x)
    ni=1 : computes k= ¤ Kn(x) dx from x to infinity
    ni=p : computes k= ¤¤..¤ Kn(x) dx from x to infinity p times
    */
	static int Kmu(int ni, double mu, double x, double& f)
	{
		const double h = 0.5; //arbitrary parameter related to the cpu time (0.5 is advised)
		//const double prec = 1.E-05; //precision parameter (1e-5 is advised)
		//const double prec = 1.E-13;
		const double relPrec = 1.E-14; //OC221011
		double absPrec = 0.; //OC221011

		if(ni < 0) return NEGATIVE_NUM_INTEG_OF_MODIF_BESSEL_FUNC;
		//if((x > 100.) && (ni == 0)) return 0; //commented out OC010609

		double k = 0., c3 = 1E+10;
		//long rr = 0;
		long long rr = 0;

		//while(c3 > prec)
		while(c3 > absPrec) //OC221011
		{
			double rrh = (++rr)*h;
			double c1 = cosh(rrh);
			double c2 = cosh(mu*rrh);

			double c1eni = 1.;
			if(ni > 0) for(int i=0; i<ni; i++) c1eni *= c1;
			
			c3 = exp(-x*c1)*c2/c1eni;
			//double arg = x*c1; //OC221011 to change when double accuracy improves
			//if((arg > -680.) && (arg < 680.)) c3 = exp(-arg)*c2/c1eni;
			//else c3 = 0.;

			k += c3;

			if(rr <= 1) absPrec = c3*relPrec; //OC221011
		}
		f = h*(0.5*exp(-x) + k);
		return 0;
	}

	/** Computes x^2*Kn(x)^2
    */
	static double AuxFuncK2d3e2_xe2(double x, void* pv=0) //OC20112018
	//static double AuxFuncK2d3e2_xe2(double x)
	{
		double KmuVal = 0;
		Kmu(0, 2./3., x, KmuVal);
		return x*x*KmuVal*KmuVal;
	}

	/** Computes x^2*Kn(x)^2
    */
	static double AuxFuncIntK5d3_x(double x, void* pv=0) //OC20112018
	//static double AuxFuncIntK5d3_x(double x)
	{
		double KmuVal = 0;
		Kmu(1, 5./3., x, KmuVal);
		return x*KmuVal;
	}

	/** Computes a definite integral ¤ x^2*Kn(x)^2 dx
    */
	static double Int_xe2_K2d3e2(double x1, double x2, double RelPrec)
	{
		if(x1 == x2) return 0;
		double K2d3_1, K2d3_2, K1d3_1, K1d3_2;
		Kmu(0, 2./3., x1, K2d3_1);
		Kmu(0, 2./3., x2, K2d3_2);
		Kmu(0, 1./3., x1, K1d3_1);
		Kmu(0, 1./3., x2, K1d3_2);
		double dK2d3e2_xe2dx1 = (2./3.)*x1*K2d3_1*K2d3_1 - 2*x1*x1*K1d3_1*K2d3_1;
		double dK2d3e2_xe2dx2 = (2./3.)*x2*K2d3_2*K2d3_2 - 2*x2*x2*K1d3_2*K2d3_2;
        return CGenMathMeth::Integ1D_FuncWithEdgeDer(&AuxFuncK2d3e2_xe2, x1, x2, dK2d3e2_xe2dx1, dK2d3e2_xe2dx2, RelPrec);
	}

	/** Computes a definite integral ¤x¤K5/3(x1)dx1dx
    */
	static double Int_x_IntK5d3(double x1, double x2, double RelPrec)
	{
		if(x1 == x2) return 0;
		double IntK5d3_1, K5d3_1, IntK5d3_2, K5d3_2;
		Kmu(0, 5./3., x1, K5d3_1);
		Kmu(1, 5./3., x1, IntK5d3_1);
		Kmu(0, 5./3., x2, K5d3_2);
		Kmu(1, 5./3., x2, IntK5d3_2);
		double Der_x1 = IntK5d3_1 - x1*K5d3_1;
		double Der_x2 = IntK5d3_2 - x2*K5d3_2;
        return CGenMathMeth::Integ1D_FuncWithEdgeDer(&AuxFuncIntK5d3_x, x1, x2, Der_x1, Der_x2, RelPrec);
	}

	/** Computes a definite integrals related to Kn
    */
	static double IntKnXn(int type, double x1, double x2, double RelPrec)
	{
		if(x1 == x2) return 0;
		if(type == 1) return Int_xe2_K2d3e2(x1, x2, RelPrec);
		else if(type == 2) return Int_x_IntK5d3(x1, x2, RelPrec);
		else return 0;
	}

	/** Calculates principal value of argument of a complex number (-Pi < Phi <= Pi)  
	@param [out] x real part
	@param [out] y imaginary part
 	@return	calculated argument value
 	@see		... */
	static double Argument(double x, double y)
	{
		const double Pi = 3.1415926535897932;
		if(x == 0)
		{
			if(y < 0) return -0.5*Pi;
			else if(y == 0) return 0;
			else return 0.5*Pi;
		}
		if(y == 0)
		{
			if(x >= 0) return 0.;
			else return Pi;
		}
		if(y < 0)
		{
			if(x < 0) return -Pi + atan(y/x);
			else return atan(y/x); // x > 0
		}
		else // y > 0
		{
			if(x < 0) return Pi + atan(y/x);
			else return atan(y/x); // x > 0
		}
	}

	/** Searches for minimum and maximum value of an array */
	static void FindMinMax(double* pData, int Np, double& Min, int& IndMin, double& Max, int& IndMax)
	{
		IndMin = IndMax = -1;
		Min = 1.E+30; Max = -1.E+30;
		if((pData == 0) || (Np <= 0)) return;
		double* tData = pData;
		for(int i=0; i<Np; i++)
		{
			double CurVal = *(tData++);
			if(Min > CurVal) { Min = CurVal; IndMin = i;}
			if(Max < CurVal) { Max = CurVal; IndMax = i;}
		}
	}
};

//*************************************************************************

class srTMathInterpol1D {

	double* AllCf;
	double** PlnCf;
	double Orig_sStart, Orig_sStep, InvStep;
	int OrigNp;

public:

	srTMathInterpol1D(double* OrigF, int InOrigNp, double InOrig_sStart, double InOrig_sStep)
	{
		if((OrigF == 0) || (InOrigNp == 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;

        AllCf = 0;
		PlnCf = 0;
        OrigNp = InOrigNp;
        Orig_sStart = InOrig_sStart;
		Orig_sStep = InOrig_sStep;
        InvStep = 1./Orig_sStep;

		SetupPolinomCfs(OrigF);
	}
	srTMathInterpol1D()
	{
		AllCf = 0;
		PlnCf = 0;
		OrigNp = 0;
	}
	~srTMathInterpol1D()
	{
		DeallocateMemoryForCfs();
	}

	void Interpolate(double sSt, double sStp, int Np, double* pInterpData);
	static double Derivative(double* f, double h, int PoIndx, int AmOfPo=5);
	static void CubicPlnCfs(double f1, double f2, double fpr1, double fpr2, double sStep, double* aa)
	{
		double f1mf2_d_s1ms2 = (f2 - f1)/sStep;
		aa[0] = f1;
		aa[1] = fpr1;
		aa[2] = (3.*f1mf2_d_s1ms2 - 2.*fpr1 - fpr2)/sStep;
		aa[3] = (-2.*f1mf2_d_s1ms2 + fpr1 + fpr2)/(sStep*sStep);
	}

private:

	void CompDerivForOrigData(double* OrigF, double* DerF);

	void SetupPolinomCfs(double* OrigF)
	{
		if((OrigF == 0) || (OrigNp <= 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
		AllocateMemoryForCfs();

		double* DerF = new double[OrigNp];
		CompDerivForOrigData(OrigF, DerF);
		CalcPlnCfs(OrigF, DerF);
		if(DerF != 0) delete[] DerF;
	}
	void CalcPlnCfs(double* OrigF, double* DerF)
	{
		if((OrigF == 0) || (DerF == 0) || (OrigNp <= 0)) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
        //int LenFieldData_m_1 = OrigNp - 1; //OC020110
        double f1 = OrigF[0], f2;
        double fpr1 = DerF[0], fpr2;
        for(int is=1; is<OrigNp; is++)
        {
            f2 = OrigF[is];
            fpr2 = DerF[is];
            CubicPlnCfs(f1, f2, fpr1, fpr2, Orig_sStep, PlnCf[is-1]);
            f1 = f2; fpr1 = fpr2;
		}
	}
    void AllocateMemoryForCfs()
    {
		if(OrigNp <= 0) throw MATH_INTERP_STRUCT_WAS_NOT_SETUP;
        DeallocateMemoryForCfs();
        int LenFieldData_m_1 = OrigNp - 1;

        PlnCf = new double*[LenFieldData_m_1];
        if(PlnCf == 0) throw MEMORY_ALLOCATION_FAILURE;
        AllCf = new double[LenFieldData_m_1*4];
        if(AllCf == 0) { delete[] PlnCf; throw MEMORY_ALLOCATION_FAILURE;}
        double* tAllCf = AllCf;
        for(int i=0; i<LenFieldData_m_1; i++) { PlnCf[i] = tAllCf; tAllCf += 4;}
    }
    void DeallocateMemoryForCfs()
    {
        if(AllCf != 0) { delete[] AllCf; AllCf = 0;}
        if(PlnCf != 0) { delete[] PlnCf; PlnCf = 0;}
    }
};

//*************************************************************************

class TMathInterpolMultiD {

public:

	static double Interp2D4pRel(double f00, double f10, double f01, double f11, double xr, double yr)
	{// assumes 0 <= xr <= 1; 0 <= yr <= 1; f00 := f|xr=0,yr=0; f10 := f|xr=1,yr=0;
		return (f11 + f00 - f10 - f01)*xr*yr + (f10 - f00)*xr + (f01 - f00)*xr + f00;
	}

};

//*************************************************************************

#endif
