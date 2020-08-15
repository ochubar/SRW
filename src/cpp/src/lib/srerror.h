/************************************************************************//**
 * File: srerror.h
 * Description: Error and Warning Messages (header)
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#ifndef __SRERROR_H
#define __SRERROR_H

#include <string>
#include <vector>
#include "srercode.h"

//using namespace std;
#include "stlstart.h"

#ifdef _MSC_VER 
#pragma warning(disable : 4786) // to suppress annoying warning from STL
#endif

//-------------------------------------------------------------------------

class CErrWarn {
protected:
	//static string error[];
	//static string warning[];
	static vector<string> error;
	static vector<string> warning;

public:

	CErrWarn();
	
	static int GetErrorSize(int ErrNo)
	{
		int i = ErrNo - FIRST_XOP_ERR;
		if((i < 0) || (i >= (int)error.size())) return 0; //(int)error[0].size();
		try { return (int)error[i].size();}
		catch(exception e) { return 0;} //(int)error[0].size();}
	}
	static int GetWarningSize(int WarnNo)
	{
		int i = WarnNo - SRW_WARNINGS_OFFSET;
		if((i < 0) || (i >= (int)warning.size())) return 0; //(int)warning[0].size();
		try { return (int)warning[i].size();}
		catch(exception e) { return 0;} //(int)warning[0].size();}
	}

	static const char* GetError(int ErrNo)
	{
		int i = ErrNo - FIRST_XOP_ERR;
		if((i < 0) || (i >= (int)error.size())) return 0; //error[0].c_str();
		try { return error[i].c_str();}
		catch(exception e) { return 0;} //error[0].c_str();}
	}
	static const char* GetWarning(int WarnNo)
	{
		int i = WarnNo - SRW_WARNINGS_OFFSET;
		if((i < 0) || (i >= (int)warning.size())) return 0; //warning[0].c_str();
		try { return warning[i].c_str();}
		catch(exception e) { return 0;} //warning[0].c_str();}
	}

	//static void AddWarningMessage(srTIntVect* pWarnMesNos, int WarnNo)
	static void AddWarningMessage(vector<int>* pWarnMesNos, int WarnNo)
	{
		//for(srTIntVect::iterator iter = pWarnMesNos->begin(); iter != pWarnMesNos->end(); ++iter)
		for(vector<int>::iterator iter = pWarnMesNos->begin(); iter != pWarnMesNos->end(); ++iter)
		{
			if(*iter == WarnNo) return;
		}
		pWarnMesNos->push_back(WarnNo);
	}

	static int ValidateArray(void* Arr, int nElem);
};

//-------------------------------------------------------------------------

#endif
