/************************************************************************//**
 * File: srsysuti.cpp
 * Description: Auxiliary OS related Utilities (probably work on Windows only)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume, R.Li
 * @version 1.0
 ***************************************************************************/

#include "srsysuti.h"

#ifdef WIN32
#include <windows.h>
#endif
#ifdef __MAC__
//#include "Memory.h" //RL26112019
#include <sys/sysctl.h>
#endif
#ifdef LINUX
//#include <proc/readproc.h>
#include <sys/sysinfo.h> //RL26112019
#endif

//*************************************************************************

double srTSystemUtils::CheckMemoryAvailable()
{
	const double OverPhysMemoryFact = 0.9; //OC28112019 //0.5; //0.7;

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
		//const double OverPhysMemoryFact = 0.9; //0.5; //0.9;
		return OverPhysMemoryFact*double(CurMemStatus.dwAvailPhys);
	}

#endif
#ifdef __MAC__

//#if 0  //RL06112019
//    // Undefined data type: Size
//	Size MaxHeapCanGrow, LagrestBlock;
//    LagrestBlock = MaxMem(&MaxHeapCanGrow);
//	const double OverPhysMemoryFact = 0.7; // To steer
//	return OverPhysMemoryFact*double(LagrestBlock);
//#endif
	//RL26112019
	int mib[2] = {CTL_HW, HW_MEMSIZE};
	unsigned long long LagrestBlock;
	size_t memlength = sizeof(LagrestBlock);
	sysctl(mib, 2, &LagrestBlock, &memlength, NULL, 0);
	//const double OverPhysMemoryFact = 0.7; // To steer
	return OverPhysMemoryFact*double(LagrestBlock);

#endif
#ifdef LINUX

	//return 8.e+09; //Is there any "standard" API for this on LINUX?
	//RL26112019
	struct sysinfo sinfo;
	sysinfo(&sinfo);
	//const double OverPhysMemoryFact = 0.7; // To steer
	return OverPhysMemoryFact*double(sinfo.totalram);

#endif
}

//*************************************************************************
