/************************************************************************//**
 * File: stlstart.h
 * Description: Auxiliary header to "load" Standard C++ Libraries (obsolete?)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __STLSTART_H
#define __STLSTART_H

//-------------------------------------------------------------------------

#if defined _MSC_VER 
#pragma warning(disable : 4786) // to suppress annoying warning from STL
#include <map>
//#include <strstream> //obsolette
#include <sstream>
#include <iostream>
using namespace std;

#elif defined  __HP_aCC
#include <map>
#include <strstream.h>
#include <iostream.h>

#elif defined __MWERKS__
#include <map.h>
#include <strstream.h>
#include <iostream.h>
using namespace std;

#else

#include <map>
//#include <strstream>
#include <sstream>
#include <iostream>
using namespace std;

#endif

//-------------------------------------------------------------------------
#endif
