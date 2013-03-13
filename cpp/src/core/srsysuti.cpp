/************************************************************************//**
 * File: srsysuti.cpp
 * Description: Auxiliary OS related Utilities (probably work on Windows only)
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
#include "srsysuti.h"
#endif

#ifdef WIN32
#include <windows.h>
#endif
#ifdef __MAC__
#include "Memory.h"
#endif
//#ifdef LINUX
//#include <proc/readproc.h> 
//#endif

//*************************************************************************

double srTSystemUtils::CheckMemoryAvailable()
{
#ifdef WIN32

	MEMORYSTATUS CurMemStatus; // Consider checking Windows version and modifying this by MEMORYSTATUSEX and GlobalMemoryStatusEx
	CurMemStatus.dwLength = sizeof(MEMORYSTATUS);
	GlobalMemoryStatus(&CurMemStatus);

	if(CurMemStatus.dwTotalPhys < 60.E+06)
	{
		const double AbsMemMin = 3.E+06; // To steer
		double PhysMemForSys = 22.E+06; // To steer
 
		double CanBeAvail = double(CurMemStatus.dwTotalPhys) - PhysMemForSys;
		if(CanBeAvail < AbsMemMin) CanBeAvail = AbsMemMin;
		double PhysMemActuallyAvail = double(CurMemStatus.dwAvailPhys);
		if(CanBeAvail < PhysMemActuallyAvail) CanBeAvail = PhysMemActuallyAvail;
		return CanBeAvail;
	}
	else
	{
		const double OverPhysMemoryFact = 0.9; //0.5; //0.9; // To steer //igor does not want to allocate big waves (which require up to 1GB)...
		return OverPhysMemoryFact*double(CurMemStatus.dwAvailPhys);
	}

#endif
#ifdef __MAC__

	Size MaxHeapCanGrow, LagrestBlock;
	LagrestBlock = MaxMem(&MaxHeapCanGrow);
	const double OverPhysMemoryFact = 0.7; // To steer
	return OverPhysMemoryFact*double(LagrestBlock);

#endif
#ifdef LINUX

	return 8.e+09; //Is there any "standard" API for this on LINUX?

#endif
}

//*************************************************************************
