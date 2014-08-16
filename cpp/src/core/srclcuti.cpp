/************************************************************************//**
 * File: srclcuti.cpp
 * Description: Auxiliary SR calculation utilities
 * Project: Synchrotron Radiation Workshop
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srclcuti.h"

//*************************************************************************

double srTCalcUtils::ChargeEl = 1.60217646263E-19; //1.602189246E-19;
double srTCalcUtils::MassEl_kg = 9.1093818872E-31; //9.10953447E-31;
double srTCalcUtils::MassEl_MeV = 0.51099890221; //0.511003414;
double srTCalcUtils::SpeedLight = 2.9979245812E+08;
double srTCalcUtils::Pi = 3.14159265358979;

//*************************************************************************

double srTCalcUtils::ConvertWavelengthM_ToPhotEnergy(double WavelengthM, int iOutUnit)
{
	//iOutUnit = 1- keV; 2- eV; 3- 1/cm; 4- Å; 5- nm; 6- µm; 7- mm
	//const double Mult = 1.239854;
	const double Mult = 1.239842; //following "X-ray Data Booklet"
	double OutRes = WavelengthM;

	if(iOutUnit == 1) 
	{
		OutRes = Mult*(1.e-9)/WavelengthM;
	}
	else if(iOutUnit == 2) 
	{
		OutRes = Mult*(1.e-6)/WavelengthM;
	}
	else if(iOutUnit == 3) 
	{
		OutRes = 1./(WavelengthM*100.); // 1/cm ???
	}
	else if(iOutUnit == 4) 
	{
		OutRes *= 1.e+10;
	}
	else if(iOutUnit == 5) 
	{
		OutRes *= 1.e+9;
	}
	else if(iOutUnit == 6) 
	{
		OutRes *= 1.e+6;
	}
	else if(iOutUnit == 7) 
	{
		OutRes *= 1000.;
	}
	return OutRes;
}

//*************************************************************************

double srTCalcUtils::MagnetRadius(double Bconst, double ElecEnergy, int iOutUnit)
{
	//iOutUnit = 1- mm; 2- m; 3- km
	double R_m = 3.335640952*ElecEnergy/Bconst;
	double OutR = R_m;

	if(iOutUnit == 1) OutR *= 1000.;
	else if(iOutUnit == 3) OutR *= 0.001;
	return OutR;
}

//*************************************************************************

double srTCalcUtils::MagnetCritPhotEn(double Bconst, double ElecEnergy, int iOutUnit)
{
	const double Pi = 3.141592653589793238;
	double R_m = 3.335640952*ElecEnergy/Bconst;

	double Gamma = ElecEnergy*1000./0.511003414;
	double Lambda_c_m = 4.*Pi*R_m/(3.*Gamma*Gamma*Gamma);
	//double OutRes = Lambda_c_m;

	return ConvertWavelengthM_ToPhotEnergy(Lambda_c_m, iOutUnit);
}

//*************************************************************************

double srTCalcUtils::UndK(double Bpeak, double Period)
{
	const double Mult = ChargeEl/(2.*Pi*MassEl_kg*SpeedLight);
	return Mult*Bpeak*Period;
}

//*************************************************************************

double srTCalcUtils::UndFundPhotEn(double Bpeak, double Period, double ElecEnergy, int iOutUnit)
{
	double K = UndK(Bpeak, Period);
	double Gamma = ElecEnergyFromGevToGamma(ElecEnergy);
	double WavelengthM = Period*(1. + 0.5*K*K)/(2.*Gamma*Gamma);
	return ConvertWavelengthM_ToPhotEnergy(WavelengthM, iOutUnit);
}

//*************************************************************************

double srTCalcUtils::ElecEnergyFromGevToGamma(double ElecEnergyGeV)
{
	return ElecEnergyGeV*1000./MassEl_MeV;
}

//*************************************************************************
