/************************************************************************//**
 * File: srobject.h
 * Description: Base class for SR related "objects"
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 0.06
 ***************************************************************************/

#ifndef __SROBJECT_H
#define __SROBJECT_H

#include <string>
#include <map>
#include "smartptr.h"

//using namespace std;
#include "stlstart.h"

//-------------------------------------------------------------------------

class CGenObject {
	string Name;

public:
	CGenObject() {}
	CGenObject(char* InName)
	{
		if(InName != 0) Name += InName;
	}
	virtual ~CGenObject() {}

	double sign(double x) { return (x < 0)? -1 : 1;}
};

//-------------------------------------------------------------------------

typedef CSmartPtr<CGenObject> CHGenObj;

#ifdef __GCC__
typedef map <int, CHGenObj, less<int> > CMHGenObj;
#else
typedef map <int, CHGenObj, less<int> > CMHGenObj;
#endif

//-------------------------------------------------------------------------

#endif
