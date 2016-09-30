/************************************************************************//**
 * File: gmrand.h
 * Description: Mathematical utilities related to pseudo-random numbers generation and Monte-Carlo integration (header)
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __GMRAND_H
#define __GMRAND_H

#include <time.h>
#include <cstdlib>

//*************************************************************************

class CGenMathRandLPTau {

	//long NR[6][20];
	//long IV[6][20];
	//long NG, IQ;
	//long NA, NB, NC, ND, NT;

	//long iQ, mQ;
	//long MaskArr[30], *MaskArrTrav;

	long long NR[6][20];
	long long IV[6][20];
	long long NG, IQ;
	long long NA, NB, NC, ND, NT;

	long long iQ, mQ;
	long long MaskArr[30], *MaskArrTrav;

public:
	CGenMathRandLPTau()
	{
		Initialize();
	}

	void Initialize()
	{
		//long LocNR1[] = {1,3,5,15,17,51,85,255,257,771,1285,3855,4369,13107,21845,65535,65537,196611,327685,983055};
		//long LocNR2[] = {1,1,7,11,13,61,67,79,465,721,823,4091,4125,4141,28723,45311,53505,250113,276231,326411};
		//long LocNR3[] = {1,3,7,5,7,43,49,147,439,1013,727,987,5889,6915,16647,49925,116487,83243,116529,715667};
		//long LocNR4[] = {1,1,5,3,15,51,125,141,177,759,267,1839,6929,16241,16565,17139,82207,50979,252717,851901};
		//long LocNR5[] = {1,3,1,1,9,59,25,89,321,835,833,4033,3913,11643,18777,35225,102401,45059,36865,299009};
		////long LocNR6[] = {1,1,3,7,31,47,109,173,181,949,471,2515,6211,2147,3169,35873,33841,99889,247315,1032727};

		long long LocNR1[] = {1,3,5,15,17,51,85,255,257,771,1285,3855,4369,13107,21845,65535,65537,196611,327685,983055};
		long long LocNR2[] = {1,1,7,11,13,61,67,79,465,721,823,4091,4125,4141,28723,45311,53505,250113,276231,326411};
		long long LocNR3[] = {1,3,7,5,7,43,49,147,439,1013,727,987,5889,6915,16647,49925,116487,83243,116529,715667};
		long long LocNR4[] = {1,1,5,3,15,51,125,141,177,759,267,1839,6929,16241,16565,17139,82207,50979,252717,851901};
		long long LocNR5[] = {1,3,1,1,9,59,25,89,321,835,833,4033,3913,11643,18777,35225,102401,45059,36865,299009};
		//long long LocNR6[] = {1,1,3,7,31,47,109,173,181,949,471,2515,6211,2147,3169,35873,33841,99889,247315,1032727};
		//dddddddddddddddddd

		//long *IV0Trav = IV[0], *IV1Trav = IV[1], *IV2Trav = IV[2], *IV3Trav = IV[3], *IV4Trav = IV[4], *IV5Trav = IV[5],
		long long *IV0Trav = IV[0], *IV1Trav = IV[1], *IV2Trav = IV[2], *IV3Trav = IV[3], *IV4Trav = IV[4], *IV5Trav = IV[5],
			 *LocNR1Trav = LocNR1, *LocNR2Trav = LocNR2, *LocNR3Trav = LocNR3, *LocNR4Trav = LocNR4, *LocNR5Trav = LocNR5,
			 *NR0Trav = NR[0], *NR1Trav = NR[1], *NR2Trav = NR[2], *NR3Trav = NR[3], *NR4Trav = NR[4], *NR5Trav = NR[5];

		for(int i=0; i<20; i++)
		{
			*(NR0Trav++) = 1; 
			*(NR1Trav++) = *LocNR1Trav;
			*(NR2Trav++) = *LocNR2Trav;
			*(NR3Trav++) = *LocNR3Trav;
			*(NR4Trav++) = *LocNR4Trav;
			*(NR5Trav++) = *LocNR5Trav;

			//long OneBuf = 1 << 21;
			long long OneBuf = 1 << 21;
			//long OneBuf = 1 << 39;
			*(IV0Trav++) = OneBuf; 
			*(IV1Trav++) = (*(LocNR1Trav++)) << 20;
			*(IV2Trav++) = (*(LocNR2Trav++)) << 19;
			*(IV3Trav++) = (*(LocNR3Trav++)) << 18;
			*(IV4Trav++) = (*(LocNR4Trav++)) << 17;
			*(IV5Trav++) = (*(LocNR5Trav++)) << 16; //???
		}
		NG = IQ = 0;
		NA = 436207616;
		NB = 872415232;
		NC = 4194303;
		ND = 536870912;
		NT = 64;

		MaskArrTrav = MaskArr; 
		//long BufNumb = 1;
		long long BufNumb = 1;
		for(int k=0; k<29; k++) { BufNumb <<= 1; *(MaskArrTrav++) = BufNumb;}
		MaskArrTrav = MaskArr;
		iQ = 0; mQ = 1;

		srand(1);
	}

/**
	void LPTauSlow(int n, double* Q)
	{
		if((++iQ) >= *MaskArrTrav) { mQ++; MaskArrTrav++;}
		for(int j=0; j<n; j++)
		{
			long* NRj = NR[j];
			double s = 0.;
			for(int k=0; k<mQ; k++)
			{
				long ns = 0;
				long *NRjL_Trav = NRj, *LocMaskArrTrav1 = &(MaskArr[k]), *LocMaskArrTrav2 = MaskArr;
				for(int L=k; L<mQ; L++)
					ns += int(2*D(double(iQ)/(*(LocMaskArrTrav1++))) + 1E-08)*int(2*D(double(*(NRjL_Trav++))/(*(LocMaskArrTrav2++))) + 1E-08);
				s += D(0.5*ns)/(1 << k);
			}
			Q[j] = s;
		}
	}
**/

	void LPTauSlow(int n, double* Q)
	{
		double a = (double)(++iQ);
		int m = 1 + int(log(a)/0.693147);
		
		for(int j = 1; j <= n; j++)
		{
			double s = 0.;
			for(int k = 1; k <= m; k++)
			{
				//int ns = 0;
				long long ns = 0;
				for(int l = k; l <= m; l++)
				{
					//double b = NR[j - 1][l - 1];
					double b = (double)NR[j - 1][l - 1];
					//ns += int(2*D(a/pow(2.,l)))*int(2*D(b/pow(2.,l+1-k)));
					ns += ((long long)(2*D(a/pow(2.,l))))*((long long)(2*D(b/pow(2.,l+1-k))));
				}
				s += D(0.5*ns)/pow(2.,k-1);
			}
			Q[j - 1] = s;
		}
	}
	
	//double D(double x) { return x - long(x);}
	double D(double x) { return x - ((long long)(x));}

	void LPTauQuick(long i, int n, double q)
	{
	}

	void SimpleRand(int n, double* Q)
	{
		double D_RAND_MAX = double(RAND_MAX);
		for(int i=0; i<n; i++) *(Q++) = double(rand())/D_RAND_MAX;
	}
};

//*************************************************************************

class CGenMathRand {

	double InvRAND_MAX;
	double PI, TwoPI, One_d_SqrtTwoPi; //, Two_d_SqrtPi, One_d_SqrtPi;
	CGenMathRandLPTau LPTau;

public:

	CGenMathRand()
	{
		Initialize();
	}
	
	void Initialize(char randInitMode = 0)
	{
		InvRAND_MAX = 1./double(RAND_MAX);
		PI = 3.141592653590;
		TwoPI = 2.*PI;
		One_d_SqrtTwoPi = 1./sqrt(TwoPI);

		if(randInitMode == 0) srand((unsigned)time(NULL));
		else srand(1);

		LPTau.Initialize();
	}

	double NextRandStd(char rand_mode = 1)
	{// Returns standard random number (>0 and <1)

		if(rand_mode == 0) return InvRAND_MAX*rand();
		else
		{
			double g[1];
			LPTau.LPTauSlow(1, g);
			return g[0];
		}
	}

	double NextRandStdN(int n, double* arrQ, char rand_mode = 1)
	{// Calculates n standard random numbers (>0 and <1)

		if(rand_mode == 0) 
		{
			for(int i=0; i<n; i++) arrQ[i] = InvRAND_MAX*rand();
			return arrQ[0];
		}
		else
		{
			LPTau.LPTauSlow(n, arrQ);
			return arrQ[0];
		}
	}

	//double NextStdRandLPTau()
	//{// Returns standard random number (>0 and <1)
	//	double g[1];
	//	LPTau.LPTauSlow(1, g);
	//	return g[0];
	//}

	void NextRandStd2D(double& x, double& y, char rand_mode = 1)
	{// Returns two standard random numbers (>0 and <1)
		double g[2];
		//LPTau.LPTauSlow(2, g);
		NextRandStdN(2, g, rand_mode);
		x = g[0]; y = g[1];
	}

	double NextRandGauss(double xc, double sigma, char rand_mode = 1)
	{
		//double g1 = NextStdRandP(), g2 = NextStdRandP();
		//double a = Sigma*sqrt(-2.*log(g1)), ph = TwoPI*g2;
		double g[2];
		NextRandStdN(2, g, rand_mode);
		double a = sigma*sqrt(-2.*log(g[0])), ph = TwoPI*g[1];
		return a*cos(ph) + xc;
	}

//	double NextRandGauss(double Xc, double Sigma)
//	{
//		double g1 = NextStdRand(), g2 = NextStdRand();
//		double a = Sigma*sqrt(-2.*log(g1)), ph = TwoPI*g2;
//		return a*cos(ph) + Xc;
///**
//		int nSig = 2;
//		double g[2];
//
//		double g0, g1;
//		while(true)
//		{
//			//LPTau.LPTauSlow(2, g);
//		
//			g[0] = NextStdRand(); g[1] = NextStdRand();
//
//			g0 = nSig*(2*g[0] - 1);
//			g1 = One_d_SqrtTwoPi*g[1];
//			if(g1 <= One_d_SqrtTwoPi*exp(-0.5*g0*g0)) break;
//		}
//		return g0;
//**/
//	}

	double NextRandGauss2(double xc, double sigma, double& gr1, double& gr2, char rand_mode = 1)
	{
		//double g1 = NextStdRand(), g2 = NextStdRand();
		//double a = sigma*sqrt(-2.*log(g1)), ph = TwoPI*g2;

		double g[2];
		NextRandStdN(2, g, rand_mode);
		double a = sigma*sqrt(-2.*log(g[0])), ph = TwoPI*g[1];

		gr1 = a*cos(ph) + xc;
		gr2 = a*sin(ph) + xc;
		return gr1;
	}

	void NextRandGauss2D(double Xc, double SigmaX, double Yc, double SigmaY, double& x, double& y, char rand_mode = 1)
	{
		double g[2];
		//LPTau.LPTauSlow(2, g);
		NextRandStdN(2, g, rand_mode);

		double a = sqrt(-2.*log(*g)), ph = TwoPI*(*(g + 1));
		x = SigmaX*a*cos(ph) + Xc; y = SigmaY*a*sin(ph) + Yc; 
	}

	void NextRandGauss4D(double* XcArr, double* SigmaXArr, double& x1, double& x2, double& x3, double& x4, char rand_mode = 1)
	{
		double g[4];
		//LPTau.LPTauSlow(4, g);
		NextRandStdN(4, g, rand_mode);

		double a = sqrt(-2.*log(g[0])), ph = TwoPI*(g[1]);
		x1 = SigmaXArr[0]*a*cos(ph) + XcArr[0]; x2 = SigmaXArr[1]*a*sin(ph) + XcArr[1]; 

		a = sqrt(-2.*log(g[2])); ph = TwoPI*(g[3]);
		x3 = SigmaXArr[2]*a*cos(ph) + XcArr[2]; x4 = SigmaXArr[3]*a*sin(ph) + XcArr[3]; 
	}

	void NextRandGauss6D(double* XcArr, double* SigmaArr, double* p6d, char rand_mode = 1)
	{
		double g[6];
		NextRandStdN(6, g, rand_mode);

		double a, ph;
		for(int i=0; i<3; i++)
		{
			int two_i = i << 1;
			int two_i_p_1 = two_i + 1;

			a = sqrt(-2.*log(g[two_i])); 
			ph = TwoPI*(g[two_i_p_1]);
            p6d[two_i] = SigmaArr[two_i]*a*cos(ph) + XcArr[two_i]; 
			p6d[two_i_p_1] = SigmaArr[two_i_p_1]*a*sin(ph) + XcArr[two_i_p_1]; 
		}
	}

	//void NextRandGaussND(double* XcArr, double* SigmaArr, double n, double* pnd, char rand_mode = 1)
	void NextRandGaussND(double* XcArr, double* SigmaArr, int n, double* pnd, char rand_mode = 1)
	{//for 1 <= n <= 6 !
		if((n <= 0) || (pnd == 0)) return;

		bool xcIsZero = false, sigmaIsOne = false;
		if(XcArr == 0) 
		{
			xcIsZero = true;
			XcArr = new double[n];
			double *tXcArr = XcArr;
			for(int i=0; i<n; i++) *(tXcArr++) = 0;
		}
		if(sigmaIsOne == 0) 
		{
			sigmaIsOne = true;
			SigmaArr = new double[n];
			double *tSigmaArr = SigmaArr;
			for(int i=0; i<n; i++) *(tSigmaArr++) = 1;
		}

		if(n > 6) n = 6;

		if(n == 6) NextRandGauss6D(XcArr, SigmaArr, pnd, rand_mode);
		else if(n == 4) NextRandGauss4D(XcArr, SigmaArr, pnd[0], pnd[1], pnd[2], pnd[3], rand_mode);
		else if(n == 2) NextRandGauss2D(XcArr[0], SigmaArr[0], XcArr[1], SigmaArr[1], pnd[0], pnd[1], rand_mode);
		else if((n == 1) || (n == 3) || (n == 5))
		{
			//int n_p_1 = n + 1;
			double arXcAux[6], arSigmaAux[6], arNdAux[6];
			double *t_arXcAux = arXcAux, *t_arSigmaAux = arSigmaAux;
			double *t_XcArr = XcArr, *t_SigmaArr = SigmaArr;
			for(int i=0; i<n; i++)
			{
				*(t_arXcAux++) = *(t_XcArr++); *(t_arSigmaAux++) = *(t_SigmaArr++);
			}
			*t_arXcAux = 0; *t_arSigmaAux = 1;

			if(n == 1) NextRandGauss2D(arXcAux[0], arSigmaAux[0], arXcAux[1], arSigmaAux[1], pnd[0], arNdAux[1], rand_mode);
			else if(n == 3) NextRandGauss4D(arXcAux, arSigmaAux, pnd[0], pnd[1], pnd[2], arNdAux[3], rand_mode);
			else 
			{
				NextRandGauss6D(arXcAux, arSigmaAux, arNdAux, rand_mode);
				double *t_arNdAux = arNdAux, *t_pnd = pnd;
				for(int i=0; i<n; i++) *(t_pnd++) = *(t_arNdAux++);
			}
		}

		if(xcIsZero && (XcArr != 0)) delete[] XcArr;
		if(sigmaIsOne && (SigmaArr != 0)) delete[] SigmaArr;
	}

	void NextRandGaussNorm2D(double& x, double& y, char rand_mode = 1)
	{
		double g[2];
		//LPTau.LPTauSlow(2, g);
		NextRandStdN(2, g, rand_mode);
		
		//g[0] = NextStdRand(); g[1] = NextStdRand();
		//x = g[0]; y = g[1];

		double a = sqrt(-2.*log(g[0]));
		double ph = TwoPI*g[1];
		x = a*cos(ph); y = a*sin(ph); 
	}
};

//*************************************************************************

#endif
