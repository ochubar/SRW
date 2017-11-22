/************************************************************************//**
 * File: gmfft.cpp
 * Description: Auxiliary utilities to work with FFTW library
 * Project: 
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @author S. Yakubov (E-XFEL) - noticed issue and suggested fix in FFT1D
 * @version 1.1
 ***************************************************************************/

#include "gmfft.h"

//*************************************************************************

long CGenMathFFT::GoodNumbers[] = {
	2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 42, 44, 
	48, 50, 52, 54, 56, 60, 64, 66, 70, 72, 78, 80, 84, 88, 90, 96, 98, 100, 104, 
	108, 110, 112, 120, 126, 128, 130, 132, 140, 144, 150, 154, 156, 160, 162, 
	168, 176, 180, 182, 192, 196, 198, 200, 208, 210, 216, 220, 224, 234, 240, 
	250, 252, 256, 260, 264, 270, 280, 286, 288, 294, 300, 308, 312, 320, 324, 
	330, 336, 350, 352, 360, 364, 378, 384, 390, 392, 396, 400, 416, 420, 432, 
	440, 448, 450, 462, 468, 480, 486, 490, 500, 504, 512, 520, 528, 540, 546, 
	550, 560, 572, 576, 588, 594, 600, 616, 624, 630, 640, 648, 650, 660, 672, 
	686, 700, 702, 704, 720, 728, 750, 756, 768, 770, 780, 784, 792, 800, 810, 
	832, 840, 858, 864, 880, 882, 896, 900, 910, 924, 936, 960, 972, 980, 990, 
	1000, 1008, 1024, 1040, 1050, 1056, 1078, 1080, 1092, 1100, 1120, 1134, 1144, 
	1152, 1170, 1176, 1188, 1200, 1232, 1248, 1250, 1260, 1274, 1280, 1296, 1300, 
	1320, 1344, 1350, 1372, 1386, 1400, 1404, 1408, 1430, 1440, 1456, 1458, 1470, 
	1500, 1512, 1536, 1540, 1560, 1568, 1584, 1600, 1620, 1638, 1650, 1664, 1680, 
	1716, 1728, 1750, 1760, 1764, 1782, 1792, 1800, 1820, 1848, 1872, 1890, 1920, 
	1944, 1950, 1960, 1980, 2000, 2002, 2016, 2048, 2058, 2080, 2100, 2106, 2112, 
	2156, 2160, 2184, 2200, 2240, 2250, 2268, 2288, 2304, 2310, 2340, 2352, 2376, 
	2400, 2430, 2450, 2464, 2496, 2500, 2520, 2548, 2560, 2574, 2592, 2600, 2640, 
	2646, 2688, 2700, 2730, 2744, 2750, 2772, 2800, 2808, 2816, 2860, 2880, 2912, 
	2916, 2940, 2970, 3000, 3024, 3072, 3080, 3120, 3136, 3150, 3168, 3200, 3234, 
	3240, 3250, 3276, 3300, 3328, 3360, 3402, 3430, 3432, 3456, 3500, 3510, 3520, 
	3528, 3564, 3584, 3600, 3640, 3696, 3744, 3750, 3780, 3822, 3840, 3850, 3888, 
	3900, 3920, 3960, 4000, 4004, 4032, 4050, 4096, 4116, 4158, 4160, 4200, 4212, 
	4224, 4290, 4312, 4320, 4368, 4374, 4400, 4410, 4480, 4500, 4536, 4550, 4576, 
	4608, 4620, 4680, 4704, 4752, 4800, 4802, 4860, 4900, 4914, 4928, 4950, 4992,
	5000, 5040, 5096, 5120, 5148, 5184, 5200, 5250, 5280, 5292, 5346, 5376, 
	5390, 5400, 5460, 5488, 5500, 5544, 5600, 5616, 5632, 5670, 5720, 5760, 5824, 
	5832, 5850, 5880, 5940, 6000, 6006, 6048, 6144, 6160, 6174, 6240, 6250, 6272, 
	6300, 6318, 6336, 6370, 6400, 6468, 6480, 6500, 6552, 6600, 6656, 6720, 6750, 
	6804, 6860, 6864, 6912, 6930, 7000, 7020, 7040, 7056, 7128, 7150, 7168, 7200, 
	7280, 7290, 7350, 7392, 7488, 7500, 7546, 7560, 7644, 7680, 7700, 7722, 7776, 
	7800, 7840, 7920, 7938, 8000, 8008, 8064, 8100, 8190, 8192, 8232, 8250, 8316, 
	8320, 8400, 8424, 8448, 8580, 8624, 8640, 8736, 8748, 8750, 8800, 8820, 8910, 
	8918, 8960, 9000, 9072, 9100, 9152, 9216, 9240, 9360, 9408, 9450, 9504, 9600, 
	9604, 9702, 9720, 9750, 9800, 9828, 9856, 9900, 9984, 10000, 10010, 10080, 
	10192, 10206, 10240, 10290, 10296, 10368, 10400, 10500, 10530, 10560, 10584, 
	10692, 10752, 10780, 10800, 10920, 10976, 11000, 11088, 11200, 11232, 11250, 
	11264, 11340, 11440, 11466, 11520, 11550, 11648, 11664, 11700, 11760, 11880, 
	12000, 12012, 12096, 12150, 12250, 12288, 12320, 12348, 12474, 12480, 12500, 
	12544, 12600, 12636, 12672, 12740, 12800, 12870, 12936, 12960, 13000, 13104, 
	13122, 13200, 13230, 13312, 13440, 13500, 13608, 13650, 13720, 13728, 13750, 
	13824, 13860, 14000, 14014, 14040, 14080, 14112, 14256, 14300, 14336, 14400, 
	14406, 14560, 14580, 14700, 14742, 14784, 14850, 14976, 15000, 15092, 15120, 
	15288, 15360, 15400, 15444, 15552, 15600, 15680, 15750, 15840, 15876, 16000, 
	16016, 16038, 16128, 16170, 16200, 16250, 16380, 16384, 16464, 16500, 16632, 
	16640, 16800, 16848, 16896, 17010, 17150, 17160, 17248, 17280, 17472, 17496, 
	17500, 17550, 17600, 17640, 17820, 17836, 17920, 18000, 18018, 18144, 18200, 
	18304, 18432, 18480, 18522, 18720, 18750, 18816, 18900, 18954, 19008, 19110, 
	19200, 19208, 19250, 19404, 19440, 19500, 19600, 19656, 19712, 19800, 19968, 
	20000, 20020, 20160, 20250, 20384, 20412, 20480, 20580, 20592, 20736, 20790, 
	20800, 21000, 21060, 21120, 21168, 21384, 21450, 21504, 21560, 21600, 21840, 
	21870, 21952, 22000, 22050, 22176, 22400, 22464, 22500, 22528, 22638, 22680, 
	22750, 22880, 22932, 23040, 23100, 23166, 23296, 23328, 23400, 23520, 23760, 
	23814, 24000, 24010, 24024, 24192, 24300, 24500, 24570, 24576, 24640, 24696, 
	24750, 24948, 24960, 25000, 25088, 25200, 25272, 25344, 25480, 25600, 25740, 
	25872, 25920, 26000, 26208, 26244, 26250, 26400, 26460, 26624, 26730, 26754, 
	26880, 26950, 27000, 27216, 27300, 27440, 27456, 27500, 27648, 27720, 28000, 
	28028, 28080, 28160, 28224, 28350, 28512, 28600, 28672, 28800, 28812, 29106, 
	29120, 29160, 29250, 29400, 29484, 29568, 29700, 29952, 30000, 30030, 30184, 
	30240, 30576, 30618, 30720, 30800, 30870, 30888, 31104, 31200, 31250, 31360, 
	31500, 31590, 31680, 31752, 31850, 32000, 32032, 32076, 32256, 32340, 32400, 
	32500, 32760, 32768, 32928, 33000, 33264, 33280, 33600, 33614, 33696, 33750, 
	33792, 34020, 34300, 34320, 34398, 34496, 34560, 34650, 34944, 34992, 35000, 
	35100, 35200, 35280, 35640, 35672, 35750, 35840, 36000, 36036, 36288, 36400, 
	36450, 36608, 36750, 36864, 36960, 37044, 37422, 37440, 37500, 37632, 37730, 
	37800, 37908, 38016, 38220, 38400, 38416, 38500, 38610, 38808, 38880, 39000, 
	39200, 39312, 39366, 39424, 39600, 39690, 39936, 40000, 40040, 40320, 40500, 
	40768, 40824, 40950, 40960, 41160, 41184, 41250, 41472, 41580, 41600, 42000, 
	42042, 42120, 42240, 42336, 42768, 42900, 43008, 43120, 43200, 43218, 43680, 
	43740, 43750, 43904, 44000, 44100, 44226, 44352, 44550, 44590, 44800, 44928, 
	45000, 45056, 45276, 45360, 45500, 45760, 45864, 46080, 46200, 46332, 46592, 
	46656, 46800, 47040, 47250, 47520, 47628, 48000, 48020, 48048, 48114, 48384, 
	48510, 48600, 48750, 49000, 49140, 49152, 49280, 49392, 49500, 49896, 49920, 
	50000, 50050, 50176, 50400, 50544, 50688, 50960, 51030, 51200, 51450, 51480, 
	51744, 51840, 52000, 52416, 52488, 52500, 52650, 52800, 52822, 52920, 53248, 
	53460, 53508, 53760, 53900, 54000, 54054, 54432, 54600, 54880, 54912, 55000, 
	55296, 55440, 55566, 56000, 56056, 56160, 56250, 56320, 56448, 56700, 56862, 
	57024, 57200, 57330, 57344, 57600, 57624, 57750, 58212, 58240, 58320, 58500, 
	58800, 58968, 59136, 59400, 59904, 60000, 60060, 60368, 60480, 60750, 61152, 
	61236, 61250, 61440, 61600, 61740, 61776, 62208, 62370, 62400, 62426, 62500, 
	62720, 63000, 63180, 63360, 63504, 63700, 64000, 64064, 64152, 64350, 64512, 
	64680, 64800, 65000, 65520, 65536, 65610, 65856, 66000, 66150, 66528, 66560, 
	67200, 67228, 67392, 67500, 67584, 67914, 68040, 68250, 68600, 68640, 68750, 
	68796, 68992, 69120, 69300, 69498, 69888, 69984, 70000, 70070, 70200, 70400, 
	70560, 71280, 71344, 71442, 71500, 71680, 72000, 72030, 72072, 72576, 72800, 
	72900, 73216, 73500, 73710, 73728, 73920, 74088, 74250, 74844, 74880, 75000, 
	75264, 75460, 75600, 75816, 76032, 76440, 76800, 76832, 77000, 77220, 77616, 
	77760, 78000, 78400, 78624, 78732, 78750, 78848, 79200, 79380, 79872, 80000, 
	80080, 80190, 80262, 80640, 80850, 81000, 81250, 81536, 81648, 81900, 81920, 
	82320, 82368, 82500, 82944, 83160, 83200, 84000, 84084, 84240, 84480, 84672, 
	85050, 85536, 85750, 85800, 86016, 86240, 86400, 86436, 87318, 87360, 87480, 
	87500, 87750, 87808, 88000, 88200, 88452, 88704, 89100, 89180, 89600, 89856, 
	90000, 90090, 90112, 90552, 90720, 91000, 91520, 91728, 91854, 92160, 92400, 
	92610, 92664, 93184, 93312, 93600, 93750, 94080, 94500, 94770, 95040, 95256, 
	95550, 96000, 96040, 96096, 96228, 96250, 96768, 97020, 97200, 97500, 98000, 
	98098, 98280, 98304, 98560, 98784, 99000, 99792, 99840, 100000 
};
long CGenMathFFT::LenGoodNumbers = 1151; //637;

long CGenMathFFT::GoodNum100s[] = { 0,37,61,79,95,107,120,130,142,151,159 };
long CGenMathFFT::LenGoodNum100s = 11;

long CGenMathFFT::GoodNum1000s[] = { 0,159,228,279,318,354,383,410,435,459,479 };
long CGenMathFFT::LenGoodNum1000s = 11;

long CGenMathFFT::GoodNum10000s[] = { 0,479,636,743,830,900,960,1017,1064,1109,1150 };
long CGenMathFFT::LenGoodNum10000s = 11;

//*************************************************************************

void CGenMathFFT::NextCorrectNumberForFFT(long& n)
{
	if(n < 4)
	{
		n = 4; return;
	}
	if(n < 100001)
	{
		long *pGoodPrev, *pGoodNext;

		long n_d_10000 = long(n*0.0001);
		if(n_d_10000 > 0) pGoodPrev = GoodNumbers + GoodNum10000s[n_d_10000] - 1;
		else
		{
			long n_d_1000 = long(n*0.001);
			if(n_d_1000 > 0) pGoodPrev = GoodNumbers + GoodNum1000s[n_d_1000] - 1;
			else
			{
				long n_d_100 = long(n*0.01);
				if(n_d_100 > 0) pGoodPrev = GoodNumbers + GoodNum100s[n_d_100] - 1;
				else pGoodPrev = GoodNumbers;
			}
		}
		pGoodNext = pGoodPrev + 1;
		for(;;)
		{
			if((n > *(pGoodPrev++)) && (n <= *pGoodNext))
			{
				n = *pGoodNext; return;
			}
			pGoodNext++;
		}
	}
	else
	{
		//long k = 16384;
		long k = 65536;
		for(int j=0; j<100; j++)
		{
			k <<= 1; 
			if(n <= k)
			{
				n = k; break;
			}
		}
	}
}

//*************************************************************************

int CGenMathFFT1D::Make1DFFT_InPlace(CGenMathFFT1DInfo& FFT1DInfo)
{
	//long TotAmOfPo = (FFT1DInfo.Nx << 1)*FFT1DInfo.HowMany;
	long long TotAmOfPo = ((long long)(FFT1DInfo.Nx << 1))*((long long)FFT1DInfo.HowMany);
	float* AuxDataCont = new float[TotAmOfPo];
	if(AuxDataCont == 0) return MEMORY_ALLOCATION_FAILURE;
	FFT1DInfo.pOutData = AuxDataCont;

	int result;
	if(result = Make1DFFT(FFT1DInfo)) return result;

	float *tOut = FFT1DInfo.pInData, *t = AuxDataCont;
	for(int ix=0; ix<TotAmOfPo; ix++) *(tOut++) = *(t++);

	if(AuxDataCont != 0) delete[] AuxDataCont;
	return 0;
}

//*************************************************************************

int CGenMathFFT2D::AuxDebug_TestFFT_Plans()
{//debug function to test why fftw2d_create_plan crashed at Nx=Nz=104

	for(long i=3; i<(CGenMathFFT::LenGoodNumbers); i++)
	{
		int CurN = GoodNumbers[i];
		fftwnd_plan Plan2DFFT;
        //fftwnd_destroy_plan(Plan2DFFT);
		Plan2DFFT = fftw2d_create_plan(CurN, CurN, FFTW_FORWARD, FFTW_IN_PLACE);
        fftwnd_destroy_plan(Plan2DFFT);
	}
	return 0;
}

//*************************************************************************
//Forward FFT (FFT2DInfo.Dir = 1?): Int f(x,y)*exp(-i*2*Pi*(qx*x + qy*y)) dx dy
//Backward FFT (FFT2DInfo.Dir = -1?): Int f(qx,qy)*exp(i*2*Pi*(qx*x + qy*y)) dqx dqy
int CGenMathFFT2D::Make2DFFT(CGenMathFFT2DInfo& FFT2DInfo)
{// Assumes Nx, Ny even !
	const double RelShiftTol = 1.E-06;

		//debug
		//AuxDebug_TestFFT_Plans();
		//end debug

	SetupLimitsTr(FFT2DInfo);

	double xStepNx = FFT2DInfo.Nx*FFT2DInfo.xStep;
	double yStepNy = FFT2DInfo.Ny*FFT2DInfo.yStep;

	double x0_After = FFT2DInfo.xStart + 0.5*xStepNx;
	double y0_After = FFT2DInfo.yStart + 0.5*yStepNy;

	NeedsShiftAfterX = (::fabs(x0_After) > RelShiftTol*xStepNx);
	NeedsShiftAfterY = (::fabs(y0_After) > RelShiftTol*yStepNy);

	double xStartTr = -0.5/FFT2DInfo.xStep;
	double yStartTr = -0.5/FFT2DInfo.yStep;

	NeedsShiftBeforeX = NeedsShiftBeforeY = 0;
	double x0_Before = 0., y0_Before = 0.;
	if(FFT2DInfo.UseGivenStartTrValues)
	{
		x0_Before = (FFT2DInfo.xStartTr - xStartTr); // Sign should be probably reversed here: check!!!
		y0_Before = (FFT2DInfo.yStartTr - yStartTr); // Sign should be probably reversed here: check!!!

		NeedsShiftBeforeX = (::fabs(x0_Before) > RelShiftTol*(::fabs(xStartTr)));
		NeedsShiftBeforeY = (::fabs(y0_Before) > RelShiftTol*(::fabs(yStartTr)));
	}

	ArrayShiftX = 0; ArrayShiftY = 0; 
	if(NeedsShiftBeforeX || NeedsShiftAfterX)
	{
		ArrayShiftX = new float[Nx << 1];
		if(ArrayShiftX == 0) return MEMORY_ALLOCATION_FAILURE;
	}
	if(NeedsShiftBeforeY || NeedsShiftAfterY)
	{
		ArrayShiftY = new float[Ny << 1];
		if(ArrayShiftY == 0) return MEMORY_ALLOCATION_FAILURE;
	}

	fftwnd_plan Plan2DFFT;
	FFTW_COMPLEX *DataToFFT = (FFTW_COMPLEX*)(FFT2DInfo.pData);

	char t0SignMult = (FFT2DInfo.Dir > 0)? -1 : 1;

	if(NeedsShiftBeforeX) FillArrayShift('x', t0SignMult*x0_Before, FFT2DInfo.xStep);
	if(NeedsShiftBeforeY) FillArrayShift('y', t0SignMult*y0_Before, FFT2DInfo.yStep);
	if(NeedsShiftBeforeX || NeedsShiftBeforeY) TreatShifts(DataToFFT);

	if(FFT2DInfo.Dir > 0)
	{
		Plan2DFFT = fftw2d_create_plan(Ny, Nx, FFTW_FORWARD, FFTW_IN_PLACE);
		if(Plan2DFFT == 0) return ERROR_IN_FFT;
		fftwnd(Plan2DFFT, 1, DataToFFT, 1, 0, DataToFFT, 1, 0);
		RepairSignAfter2DFFT(DataToFFT);
		RotateDataAfter2DFFT(DataToFFT);
	}
	else
	{
		Plan2DFFT = fftw2d_create_plan(Ny, Nx, FFTW_BACKWARD, FFTW_IN_PLACE);
		if(Plan2DFFT == 0) return ERROR_IN_FFT;
		RotateDataAfter2DFFT(DataToFFT);
		RepairSignAfter2DFFT(DataToFFT);
		fftwnd(Plan2DFFT, 1, DataToFFT, 1, 0, DataToFFT, 1, 0);
	}
	
	//double Mult = FFT2DInfo.xStep*FFT2DInfo.yStep;
	double Mult = FFT2DInfo.xStep*FFT2DInfo.yStep*FFT2DInfo.ExtraMult; //OC20112017

	NormalizeDataAfter2DFFT(DataToFFT, Mult);

	if(NeedsShiftAfterX) FillArrayShift('x', t0SignMult*x0_After, FFT2DInfo.xStepTr);
	if(NeedsShiftAfterY) FillArrayShift('y', t0SignMult*y0_After, FFT2DInfo.yStepTr);
	if(NeedsShiftAfterX || NeedsShiftAfterY) TreatShifts(DataToFFT);

	//OC_NERSC: to comment-out the following line for NERSC (to avoid crash with "python-mpi")
	fftwnd_destroy_plan(Plan2DFFT);

	if(ArrayShiftX != 0) 
	{
		delete[] ArrayShiftX; ArrayShiftX = 0;
	}
	if(ArrayShiftY != 0)
	{
		delete[] ArrayShiftY; ArrayShiftY = 0;
	}
	return 0;
}

//*************************************************************************
//Forward FFT: Int f(x)*exp(-i*2*Pi*qx*x)dx
//Backward FFT: Int f(qx)*exp(i*2*Pi*qx*x)dqx
int CGenMathFFT1D::Make1DFFT(CGenMathFFT1DInfo& FFT1DInfo)
{// Assumes Nx, Ny even !
	const double RelShiftTol = 1.E-06;

	SetupLimitsTr(FFT1DInfo);

	double xStepNx = FFT1DInfo.Nx*FFT1DInfo.xStep;
	double x0_After = FFT1DInfo.xStart + 0.5*xStepNx;
	NeedsShiftAfterX = FFT1DInfo.ApplyAutoShiftAfter && (::fabs(x0_After) > RelShiftTol*xStepNx);

	double xStartTr = -0.5/FFT1DInfo.xStep;

	NeedsShiftBeforeX = 0;
	double x0_Before = 0.;

	if(FFT1DInfo.UseGivenStartTrValue)
	{
		x0_Before = (FFT1DInfo.xStartTr - xStartTr);
		NeedsShiftBeforeX = (::fabs(x0_Before) > RelShiftTol*(::fabs(xStartTr)));
	}

	m_ArrayShiftX = 0;
	if(NeedsShiftBeforeX || NeedsShiftAfterX)
	{
		m_ArrayShiftX = new float[Nx << 1];
		if(m_ArrayShiftX == 0) return MEMORY_ALLOCATION_FAILURE;
	}

	fftw_plan Plan1DFFT;
	FFTW_COMPLEX *DataToFFT = (FFTW_COMPLEX*)(FFT1DInfo.pInData);
	FFTW_COMPLEX *OutDataFFT = (FFTW_COMPLEX*)(FFT1DInfo.pOutData);

	FFTW_COMPLEX *pOutDataFFT = OutDataFFT; //OC03092016 to be used solely in fftw call
/**
	Pointed-out by Sergey Yakubov (E-XFEL).
	From FFTW 2.1.5 docs:
	void fftw(fftw_plan plan, int howmany,
          fftw_complex *in, int istride, int idist,
          fftw_complex *out, int ostride, int odist);
	...
	out, ostride and odist describe the output array(s). The format is the same as for the input array. 
	In-place transforms:  If the plan specifies an in-place transform, ostride and odist are always ignored. 
	If out is NULL, out is ignored, too. Otherwise, out is interpreted as a pointer to an array of n complex numbers, 
	that FFTW will use as temporary space to perform the in-place computation. out is used as scratch space and its contents destroyed. 
	In this case, out must be an ordinary array whose elements are contiguous in memory (no striding). 
**/

	char t0SignMult = (FFT1DInfo.Dir > 0)? -1 : 1;
	if(NeedsShiftBeforeX) 
	{
		FillArrayShift(t0SignMult*x0_Before, FFT1DInfo.xStep);
		TreatShift(DataToFFT, FFT1DInfo.HowMany);
	}

	if(FFT1DInfo.Dir > 0)
	{
		int flags = FFTW_ESTIMATE;
		if(DataToFFT == OutDataFFT)
		{
			flags |= FFTW_IN_PLACE;
			pOutDataFFT = 0; //OC03092016 (see FFTW 2.1.5 doc clause above)
		}
		Plan1DFFT = fftw_create_plan(Nx, FFTW_FORWARD, flags);
		if(Plan1DFFT == 0) return ERROR_IN_FFT;

		//fftw(Plan1DFFT, FFT1DInfo.HowMany, DataToFFT, 1, Nx, OutDataFFT, 1, Nx);
		fftw(Plan1DFFT, FFT1DInfo.HowMany, DataToFFT, 1, Nx, pOutDataFFT, 1, Nx); //OC03092016

		RepairSignAfter1DFFT(OutDataFFT, FFT1DInfo.HowMany);
		RotateDataAfter1DFFT(OutDataFFT, FFT1DInfo.HowMany);
	}
	else
	{
		int flags = FFTW_ESTIMATE;
		if(DataToFFT == OutDataFFT)
		{
			flags |= FFTW_IN_PLACE;
			pOutDataFFT = 0; //OC03092016 (see FFTW 2.1.5 doc clause above)
		}
		Plan1DFFT = fftw_create_plan(Nx, FFTW_BACKWARD, flags);
		if(Plan1DFFT == 0) return ERROR_IN_FFT;

		RotateDataAfter1DFFT(DataToFFT, FFT1DInfo.HowMany);
		RepairSignAfter1DFFT(DataToFFT, FFT1DInfo.HowMany);

		//fftw(Plan1DFFT, FFT1DInfo.HowMany, DataToFFT, 1, Nx, OutDataFFT, 1, Nx);
		fftw(Plan1DFFT, FFT1DInfo.HowMany, DataToFFT, 1, Nx, pOutDataFFT, 1, Nx); //OC03092016
	}
	//double Mult = FFT1DInfo.xStep;
	double Mult = FFT1DInfo.xStep*FFT1DInfo.MultExtra;
	NormalizeDataAfter1DFFT(OutDataFFT, FFT1DInfo.HowMany, Mult);

	if(NeedsShiftAfterX)
	{
		FillArrayShift(t0SignMult*x0_After, FFT1DInfo.xStepTr);
		TreatShift(OutDataFFT, FFT1DInfo.HowMany);
	}

	if(FFT1DInfo.TreatSharpEdges)
	{
		int result = ProcessSharpEdges(FFT1DInfo);
		if(result) return result;
	}

	//OC_NERSC: to comment-out the following line for NERSC (to avoid crash with "python-mpi")
	fftw_destroy_plan(Plan1DFFT);

	if(m_ArrayShiftX != 0) 
	{
		delete[] m_ArrayShiftX; m_ArrayShiftX = 0;
	}
	return 0;
}

//*************************************************************************

int CGenMathFFT1D::SetupAuxDataForSharpEdgeCorr(CGenMathFFT1DInfo& FFT1DInfo, CGenMathAuxDataForSharpEdgeCorr1D& AuxDataForSharpEdgeCorr)
{
	double Step = FFT1DInfo.xStep, Start = FFT1DInfo.xStart;
	double AbsTol = 0.05*Step;

	double EdgeMinOffsetFromStart = FFT1DInfo.LeftSharpEdge - Start;
	long iEdgeMinLower = long(EdgeMinOffsetFromStart/Step + 1.E-04); // Steer: threr was a bug at 1.E-08 and less!
	double EdgeMinLowerMisfit = EdgeMinOffsetFromStart - iEdgeMinLower*Step;

	double EdgeMaxOffsetFromStart = FFT1DInfo.RightSharpEdge - Start;
	long iEdgeMaxLower = long(EdgeMaxOffsetFromStart/Step + 1.E-04); // Steer: threr was a bug at 1.E-08 and less!
	double EdgeMaxLowerMisfit = EdgeMaxOffsetFromStart - iEdgeMaxLower*Step;

	char EdgeMinIsBetweenMeshPoints = (EdgeMinLowerMisfit > AbsTol);
	char EdgeMaxIsBetweenMeshPoints = (EdgeMaxLowerMisfit > AbsTol);
	char EdgeMaxIsSmallerThanDataEnd = (::fabs((Start + FFT1DInfo.Nx*Step) - FFT1DInfo.RightSharpEdge) > AbsTol);
	char EdgeCorrNeeded = (EdgeMinIsBetweenMeshPoints || EdgeMaxIsBetweenMeshPoints || EdgeMaxIsSmallerThanDataEnd);

	float dSt = 0.;
	if(EdgeMinIsBetweenMeshPoints) dSt = (float)(Step - EdgeMinLowerMisfit);
	float dFi = 0.;
	if(EdgeMaxIsBetweenMeshPoints) dFi = (float)(Step - EdgeMaxLowerMisfit);
	else if(EdgeMaxIsSmallerThanDataEnd) dFi = (float)(0.5*Step);

	CGenMathFFT1DInfo FFT1DInfoLoc = FFT1DInfo;
	FFT1DInfoLoc.UseGivenStartTrValue = 0;
	CGenMathFFT1D FFT1D;
	FFT1D.SetupLimitsTr(FFT1DInfoLoc);

	if(EdgeCorrNeeded)
	{
		AuxDataForSharpEdgeCorr.d = Step;
		long TwoN = FFT1DInfo.Nx << 1;

		if(dSt != 0.)
		{
			AuxDataForSharpEdgeCorr.ExpArrSt = new float[TwoN];
			if(AuxDataForSharpEdgeCorr.ExpArrSt == 0) return MEMORY_ALLOCATION_FAILURE;

			AuxDataForSharpEdgeCorr.dSt = dSt;
			long jSt = iEdgeMinLower + 1;
			AuxDataForSharpEdgeCorr.iSt = jSt;

			double ArgjSt = Start + jSt*Step;
			SetupSharpEdgeExpCorrArray(AuxDataForSharpEdgeCorr.ExpArrSt, FFT1DInfoLoc.Nx, ArgjSt, FFT1DInfoLoc.xStartTr, FFT1DInfoLoc.xStepTr);
		}
		if(dFi != 0.)
		{
			AuxDataForSharpEdgeCorr.ExpArrFi = new float[TwoN];
			if(AuxDataForSharpEdgeCorr.ExpArrFi == 0) return MEMORY_ALLOCATION_FAILURE;

			AuxDataForSharpEdgeCorr.dFi = dFi;
			double ArgjFi = Start + iEdgeMaxLower*Step;
			AuxDataForSharpEdgeCorr.iFi = iEdgeMaxLower;

			SetupSharpEdgeExpCorrArray(AuxDataForSharpEdgeCorr.ExpArrFi, FFT1DInfoLoc.Nx, ArgjFi, FFT1DInfoLoc.xStartTr, FFT1DInfoLoc.xStepTr);
		}
		AuxDataForSharpEdgeCorr.WasSetUp = 1;
	}
	return 0;
}

//*************************************************************************

void CGenMathFFT1D::MakeSharpEdgeCorr(CGenMathFFT1DInfo& FFT1DInfo, CGenMathAuxDataForSharpEdgeCorr1D& AuxData)
{
	float *t = FFT1DInfo.pOutData;

	float *tSt = FFT1DInfo.pInData + (AuxData.iSt << 1);
	float fSRe = *tSt, fSIm = *(tSt + 1);

	float *tFi = FFT1DInfo.pInData + (AuxData.iFi << 1);
	float fFRe = *tFi, fFIm = *(tFi + 1);

	for(long i=0; i<FFT1DInfo.Nx; i++)
	{
		long Two_i = i << 1;
		long Two_i_p_1 = Two_i + 1;

		float Re = *t, Im = *(t+1);
		if(AuxData.dSt != 0.)
		{
			float ExpStRe = AuxData.ExpArrSt[Two_i], ExpStIm = AuxData.ExpArrSt[Two_i_p_1];
			Re += (float)(AuxData.dSt*(ExpStRe*fSRe - ExpStIm*fSIm));
			Im += (float)(AuxData.dSt*(ExpStRe*fSIm + ExpStIm*fSRe));
		}
		if(AuxData.dFi != 0.)
		{
			float ExpFiRe = AuxData.ExpArrFi[Two_i], ExpFiIm = AuxData.ExpArrFi[Two_i_p_1];
			Re -= (float)(AuxData.dFi*(ExpFiRe*fFRe - ExpFiIm*fFIm));
			Im -= (float)(AuxData.dFi*(ExpFiRe*fFIm + ExpFiIm*fFRe));
		}
		*t = Re; *(t+1) = Im;
		t += 2;
	}
}

//*************************************************************************
