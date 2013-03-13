/************************************************************************//**
 * File: srremflp.cpp
 * Description: Auxiliary (obsolete or rarely used) class for processing Radiation data (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRREMFLP_H
#define __SRREMFLP_H

#include "srstraux.h"

//*************************************************************************

class srTAuxRemoveFlips {
public:

	int GenRemoveFlips(srTWaveAccessData& WaveData);
	void RemoveFlips1D(DOUBLE* pData, long Np, long i0, double Phi0);
	void RemoveFlips1D(float* pData, long Np, long i0, double Phi0);
	int RemoveFlips2D_D(srTWaveAccessData& WaveData);
	int RemoveFlips2D_F(srTWaveAccessData& WaveData);
};

//*************************************************************************

#endif
