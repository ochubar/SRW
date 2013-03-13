/************************************************************************//**
 * File: srmemory.h
 * Description: Auxiliary memory leak search (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author P.Elleaume, O.Chubar
 * @version 0.06
 ***************************************************************************/

#ifndef __SRMEMORY_H
#define __SRMEMORY_H

//*************************************************************************

#ifdef __MAC__

void* operator new(size_t size)
{
	return NewPtr(size);
}

void operator delete(void* p)
{
    DisposePtr((char*)p);
}

#endif

//*************************************************************************
// Use this to trace and repair memory leaks
/**
#ifdef __VC__
#ifdef _DEBUG

extern long UseOfNew;
extern long UseOfDelete;

//Repair
//	UseOfNew = 0;
//	UseOfDelete = 0;
//End Repair

void* operator new(unsigned int cb)
{
	if(UseOfNew > 10000000) UseOfNew = 0;
	UseOfNew++;

	return malloc(cb);
}

void operator delete(void* p)
{
	if(UseOfDelete > 10000000) UseOfDelete = 0;
	UseOfDelete++;

    free(p);
}

#endif
#endif
**/
//*************************************************************************

#endif
