/************************************************************************//**
 * File: smartptr.h
 * Description: Smart pointer (referece counting) class, taken from Stroustrup, "The C++ Programming Language"
 * Project: 
 * First release: 1997
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author P.Elleaume, O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __SMARTPTR_H
#define __SMARTPTR_H

//-------------------------------------------------------------------------

template<class T> class CSmartPtr {
public:
	T* rep;
	int* pcount;
	bool dontDelPtr; //OC13112010

	CSmartPtr () { rep=0; pcount=0; dontDelPtr=false;}
	CSmartPtr (T* pp, bool dontDel =false) : rep(pp), pcount(new int) { (*pcount)=1; dontDelPtr=dontDel;}
	CSmartPtr (const CSmartPtr& r) : rep(r.rep), pcount(r.pcount), dontDelPtr(r.dontDelPtr)
	{ 
		if(pcount != 0) (*pcount)++;
	}

	void destroy()
	{
		if(pcount!=0)
		{
			if(--(*pcount)==0)
			{
				if((!dontDelPtr) && (rep!=0)) delete rep; //OC13112010
				delete pcount;
				rep=0; pcount=0;
			}
		}
	}

	T* operator->() { return rep;}
	T* ptr() { return rep;}
	void bind(const CSmartPtr& r)
	{
		if(rep!=r.rep)
		{
			if(r.rep!=0)
			{
				destroy();
				rep = r.rep;
				pcount = r.pcount;
				(*pcount)++;
			}
			else
			{
				rep = 0;
				pcount = 0;
			}
			dontDelPtr = r.dontDelPtr; //OC13112010
		}
	}

	CSmartPtr& operator=(const CSmartPtr& r) 
	{ 
		bind(r); return *this;
	}

	int operator<(const CSmartPtr& r)
	{
		if(rep<r.rep) return 1;
		else return 0;
	}
	int operator==(const CSmartPtr& r)
	{
		if(rep==r.rep) return 1;
		else return 0;
	}

	~CSmartPtr()
	{
		destroy();
	}
};

//-------------------------------------------------------------------------

template<class T> inline int operator <(const CSmartPtr<T>& h1, const CSmartPtr<T>& h2)
{
	return (h1.rep < h2.rep); 
}

//-------------------------------------------------------------------------

template<class T> inline int operator ==(const CSmartPtr<T>& h1, const CSmartPtr<T>& h2)
{
	return (h1.rep == h2.rep); 
}

//-------------------------------------------------------------------------

#endif
