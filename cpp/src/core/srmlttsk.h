/************************************************************************//**
 * File: srmlttsk.h
 * Description: Yield (obsolete, was used on Mac before OSX)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author P.Elleaume, O.Chubar
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMLTTSK_H
#define __SRMLTTSK_H

//#ifdef __IGOR_PRO__
//#ifndef __SRIGINTR_H
//#include "srigintr.h"
//#endif
//#else
//#ifndef __SRIGORRE_H
//#include "srigorre.h"
//#endif
//#endif

#include <time.h>

//*************************************************************************

extern int gCallSpinProcess;

//*************************************************************************

class srTYield {
private:
	clock_t oldtime;
	clock_t delta;
public:
	srTYield() { delta=0;}

	inline void Init(double t);
	inline int Check();
};

//*************************************************************************

inline int srTYield::Check() 
{
	if(delta <= 0) return 0;
	if((clock() > oldtime) && gCallSpinProcess) 
	{
		//try
		//{
		if(SpinProcess()) 
		{
			return SR_COMP_PROC_ABORTED;
		}
		//}
		//catch(int)
		//{
		//	return SR_COMP_PROC_ABORTED;
		//}

		oldtime = clock() + delta;
	}
	return 0; // normal
}

//*************************************************************************

inline void srTYield::Init(double t)
{
	if(t<=0)
	{
		delta=0; return;
	}
	delta = (clock_t)(CLOCKS_PER_SEC*t);
	if(delta == 0) delta = 1;

	oldtime = clock() + delta;
	return;
}

//*************************************************************************

#endif
