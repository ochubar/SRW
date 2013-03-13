/************************************************************************//**
 * File: srigorre.h
 * Description: "Replacement" of some IGOR Pro / Mac types (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRIGORRE_H
#define __SRIGORRE_H

//-------------------------------------------------------------------------

#ifndef __IGOR_PRO__

#define MAX_OBJ_NAME 255

#ifdef _WINDOWS
#include "windows.h"
#endif

typedef double DOUBLE;
typedef char* Ptr;
typedef Ptr* Handle;

#define waveHndl Handle
#define NIL 0

int SpinProcess(void);
//HOST_EXPORT Handle NewHandle(long size);
Handle NewHandle(long size);

#endif

//-------------------------------------------------------------------------
#endif

