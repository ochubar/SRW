/* *********************************************************************
#
#
# Project:       Alpha
#
# Description:   Database class  to gather object 
#
# Author(s):     Oleg Chubar, Pascal Elleaume, Laurent Farvacque
#
# Original:	 May 2000
#
# 
# Copyright (c) 2000 by European Synchrotron Radiation Facility,
#                       Grenoble, France
#
#                       All Rights Reserved
#
#********************************************************************** */

#ifndef __OBJCONT_H
#define __OBJCONT_H

#ifndef __SMARTPTR_H
#include "smartptr.h"
#endif

//-------------------------------------------------------------------------

template<class T> class CObjCont {
	
	int pos;

public:

	map<int, CSmartPtr<T>, less<int> > data; // to enable simple iteration
	
	CObjCont(CObjCont<T>& InObjCont) { copy(InObjCont);}

	CObjCont() { pos=0;}
	~CObjCont() { erase();}
	
	int insert(T* i) 
	{
		CSmartPtr<T> p(i);
		data[++pos]=p;
		return pos;
	}
	int insert(const CSmartPtr<T>& p) 
	{
		data[++pos]=p;
		return pos;
	}
	int insert(int in_pos, const CSmartPtr<T>& p) 
	{
		pos = in_pos;
		data[pos]=p;
		return in_pos;
	}

	void copy(CObjCont<T>& InObjCont)
	{
		if(InObjCont.size() <= 0) return;
		//OC020110: modified because "map<int, CSmartPtr<T>, less<int> >::const_iterator iter;" doesn't compile on GCC 4.2 (on MAC OSX)
		//for(map<int, CSmartPtr<T>, less<int> >::const_iterator iter = InObjCont.data.begin(); iter != InObjCont.data.end(); ++iter)
		//{
		//	insert((*iter).second);
		//}
		data.insert(InObjCont.data.begin(), InObjCont.data.end());
	}

	void erase(int j) 
	{
		if(exists(j)) 
		{
			data.erase(data.find(j));
			if(pos == j) pos--;
		}
	}
	void erase() 
	{
		data.erase(data.begin(), data.end());
		pos = 0;
	}
	
	int exists(int j) 
	{
		if(data.find(j) == data.end()) return 0;
		return 1;
	}

	int size() { return (int)data.size();}
	int getPos() { return pos;}

	CSmartPtr<T> get(int j) { return (data.find(j))->second;}
	T* getPtr(int j) { return ((data.find(j))->second).ptr();}
};

//-------------------------------------------------------------------------

#endif
