/************************************************************************//**
 * File: srsysuti.h
 * Description: Auxiliary OS related Utilities (probably work on Windows only) header
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRSYSUTI_H
#define __SRSYSUTI_H

//#ifdef __IGOR_PRO__
//#ifndef __SRIGINTR_H
//#include "srigintr.h"
//#endif
//#else
//#ifndef __SRIGORRE_H
//#include "srigorre.h"
//#endif
//#endif

//*************************************************************************

class srTSystemUtils {
public:

	static double CheckMemoryAvailable();
};

//*************************************************************************

#endif
