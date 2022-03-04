// SetupPrep.h : main header file for the SETUPPREP application
//

#if !defined(AFX_SETUPPREP_H__B9F7C588_788D_4EAC_8FD3_C62D8B2B481A__INCLUDED_)
#define AFX_SETUPPREP_H__B9F7C588_788D_4EAC_8FD3_C62D8B2B481A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

#include <sstream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// CSetupPrepApp:
// See SetupPrep.cpp for the implementation of this class
//

class CSetupPrepApp : public CWinApp
{
	static const char gSRWFlatDirName[];
	static const char gSRWStructMapFileName[];
	static char gCharSepar;
	static char gPathPackagerDir[MAX_PATH];
	static char gPackagerFilename[MAX_PATH];
	static const char gPathSetupExe[];
	static const char gDefaultPathPackagerDir[];
	static const char gDefaultPackagerFilename[];

public:
	CSetupPrepApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSetupPrepApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	bool GetPathToFile(char* Explan, char* PathSRW);
	bool GetPathToFolder(char* DialogTitle, char* pWhereToPut);
	bool CreateSRWDirStructDescrFileAndFlatDir(char* PathSRW, char* PathSetupExe, char* PathTempDir);
	int TreatFilesWithinDirectory(char* CurrentDir, ostringstream& ostr, char* PathTempDir);
	int DeleteDirectoryNested(char* DirName);
	bool LaunchPackager();


	//{{AFX_MSG(CSetupPrepApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SETUPPREP_H__B9F7C588_788D_4EAC_8FD3_C62D8B2B481A__INCLUDED_)
