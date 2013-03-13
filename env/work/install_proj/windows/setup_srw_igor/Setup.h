// Setup.h : main header file for the SETUP application
//

#if !defined(AFX_SETUP_H__79681F56_85D3_11D2_B308_0060088DEA10__INCLUDED_)
#define AFX_SETUP_H__79681F56_85D3_11D2_B308_0060088DEA10__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols
#include "InstProgrDlg.h"

#include <sstream>
#include <vector>

/////////////////////////////////////////////////////////////////////////////

using namespace std;
typedef vector<char*, allocator<char*> > TCharPtrVect;

/////////////////////////////////////////////////////////////////////////////
// CSetupApp:
// See Setup.cpp for the implementation of this class
//

class CSetupApp : public CWinApp
{
	static char gSRWStructMapFileName[];
	static char gCharSepar;

	//char StartDirNameG[1024];
	//DWORD LenStartDirNameG;
	int ArbSearchStartedG;

	//CInstProgrDlg* gpProgrDlg;

public:
	CSetupApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSetupApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	int DoAutoInstallSRW();
	//int DoManInstallSRW();

	int FindPathGlobal(char* ToLookFor, DWORD FileType, char* StartName, char* GlobPath, int& LenGlobPath);
	int SearchWithinDirectory(char *CurrentDir, char* ToLookFor, DWORD FileType, char* StartName, char* GlobPath, int& LenGlobPath);
	int SearchWithinDirectoryNotNested(char* CurrentDir, char* ToLookFor, DWORD FileType, char* GlobPath, int& LenGlobPath);
	int DeinstallPreviousVersionOfSRW(char* GlobPath, char ReportOtherVersionFound);
	int DeleteDirectoryNested(char* DirName);

	void ReportUnsuccessfulInstallation();
	void ReportUnsuccessfulInstallationSRWTempNotFound();
	void ReportUnsuccessfulInstallationIgorNotFound();
	void ReportSuccessfulInstallation();
	void ReportDeliveryAndNeedCompletingTheInstallation(char*);
	int ReportConfirmDeletePrevVersionSRW(char* PathPrevVersion);

	int MakeSRWDirectory(char* PathIgorProFolder, char* PathSRWTempFolder);

	//void SRWDirectoryStructure(ostringstream& SRW_DirStruct);
	int DeleteFilePatient(char* FileName);
	int RemoveDirectoryPatient(char* DirName);
	int RemoveDirectoryNamedIfExists(const char* LocName, const char* WhereItIs);
	int CreateDirectoryPatient(const char* DirName, SECURITY_ATTRIBUTES* pSecAttr);
	int SetCurrentDirectoryPatient(const char* DirName);
	int CopyFilePatient(const char* ExistingFilePath, const char* NewFilePath, BOOL Par);
	int PlaceShortCuts(char* PathIgorProFolder);
	int MakeShortCutTo(char* PathOrigFile, char* NameShortCut);
	int ResolveIt(HWND hwnd, LPCSTR lpszLinkFile, LPSTR lpszPath, bool GetShortPath);

	int GetFolderFromUser(char* DialogTitle, char* pWhereToPut);

	int FindSRWTempFolder(char* PathSRWTempFolder);
	int FindIgorProFolder(char* PathIgorProFolder);
	int CheckFoldersInsideIgorFolder(char* PathIgorProFolder);
	int AskUserAboutIgorProFolder(char* PathIgorProFolder);
	int ReportIncorrectIgorFolderAndAskToTryAgain();
	int ReportCanNotUninstallSuggestRetry();
	int ReportCanNotCreateDirectorySuggestRetry();
	int ReportCanNotPlaceShortcutsSuggestRetry();
	void ReportSRW_FolderCreatedButInstallationWasNotCompleted();

	void EraseTerminatingBackslash(char* str);
	void EraseAfterLastBackslashInclusive(char* str);
	int RemoveShortcutsToSRWFilesFromCurrentDir();
	int CheckIfNameIsPointedBySRWShortcut(char* FileNameToCheck);
	BOOL CheckIfInstallationShouldBeCancelled();


	//{{AFX_MSG(CSetupApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SETUP_H__79681F56_85D3_11D2_B308_0060088DEA10__INCLUDED_)
