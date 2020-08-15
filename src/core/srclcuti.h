/************************************************************************//**
 * File: srclcuti.h
 * Description: Auxiliary SR calculation utilities (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRCLCUTI_H
#define __SRCLCUTI_H

//*************************************************************************

class srTCalcUtils {

	static double ChargeEl;
	static double MassEl_kg;
	static double MassEl_MeV;
	static double SpeedLight;
	static double Pi;

	static double ConvertWavelengthM_ToPhotEnergy(double WavelengthM, int OutUnit);

public:

	static double ElecEnergyFromGevToGamma(double ElecEnergyGeV);
	static double MagnetRadius(double Bconst, double ElecEnergy, int OutUnit);
	static double MagnetCritPhotEn(double Bconst, double ElecEnergy, int OutUnit);
	static double UndK(double Bpeak, double Period);
	static double UndFundPhotEn(double Bpeak, double Period, double ElecEnergy, int OutUnit);
};

//*************************************************************************

#endif
