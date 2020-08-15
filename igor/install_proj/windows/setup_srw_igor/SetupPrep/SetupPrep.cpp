// SetupPrep.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "SetupPrep.h"
#include "SetupPrepDlg.h"
#include <process.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////

const char CSetupPrepApp::gSRWFlatDirName[] = "SRW_Flat_Temp\0";
const char CSetupPrepApp::gSRWStructMapFileName[] = "SRW_Temp_Directory_Locator.txt\0";
char CSetupPrepApp::gCharSepar = '*';
char CSetupPrepApp::gPathPackagerDir[MAX_PATH]; // = "D:\\SRW_Dev\\distribution\\installers\\GkSetup Self Extractor";
char CSetupPrepApp::gPackagerFilename[MAX_PATH]; // = "SFXGEN.EXE";

const char CSetupPrepApp::gPathSetupExe[] = "C:\\SoftwareDevelopments\\SRW_Dev\\distribution\\windows\\setup\\Release\\Setup.exe"; //"E:\\SoftwareDevelopments\\SRW_Dev\\distribution\\windows\\setup\\Release\\Setup.exe";
const char CSetupPrepApp::gDefaultPathPackagerDir[] = "E:\\SoftwareDevelopments\\SRW_Dev\\distribution\\windows\\installers\\GkSetup Self Extractor"; //"D:\\SRW_Dev\\distribution\\installers\\GkSetup Self Extractor";
const char CSetupPrepApp::gDefaultPackagerFilename[] = "SFXGEN.EXE";

/////////////////////////////////////////////////////////////////////////////
// CSetupPrepApp

BEGIN_MESSAGE_MAP(CSetupPrepApp, CWinApp)
	//{{AFX_MSG_MAP(CSetupPrepApp)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSetupPrepApp construction

CSetupPrepApp::CSetupPrepApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CSetupPrepApp object

CSetupPrepApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CSetupPrepApp initialization

BOOL CSetupPrepApp::InitInstance()
{
	AfxEnableControlContainer();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

#ifdef _AFXDLL
	Enable3dControls();			// Call this when using MFC in a shared DLL
#else
	//Enable3dControlsStatic();	// Call this when linking to MFC statically
#endif

	char PathSRW[MAX_PATH], PathSetupExe[MAX_PATH], PathTempDir[MAX_PATH];

	char Explan1[] = "Please specify path of the SRW folder to pack:";
	char Explan2[] = "Please specify path of the SRW Setup.exe file:";
	char Explan3[] = "Specify where a flat SRW directory should be created:";

	if(!GetPathToFolder(Explan1, PathSRW)) return FALSE;
    strcpy(PathSetupExe, gPathSetupExe);
	if(!GetPathToFile(Explan2, PathSetupExe)) return FALSE;
	if(!GetPathToFolder(Explan3, PathTempDir)) return FALSE;

	if(!CreateSRWDirStructDescrFileAndFlatDir(PathSRW, PathSetupExe, PathTempDir)) return FALSE;
	if(!LaunchPackager()) return FALSE;

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return FALSE;
}

/////////////////////////////////////////////////////////////////////////////
bool CSetupPrepApp::LaunchPackager()
{
	char MesTitle[] = "SRW flat directory prepared";
	char Str[] = "SRW flat directory was successfully prepared.\nDistribution pack creation procedure will start now.";
	UINT DlgStyle = MB_OK | MB_ICONASTERISK | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, Str, MesTitle, DlgStyle); 

	strcpy(gPackagerFilename, gDefaultPackagerFilename);
	strcpy(gPathPackagerDir, gDefaultPathPackagerDir);

PackagerStuff:
	if(!SetCurrentDirectory(gPathPackagerDir))
	{
		char PackagerDir[MAX_PATH];
		char Explan4[] = "Specify location of the packager (GkSetup SFXGEN.EXE):";
		if(!GetPathToFile(Explan4, PackagerDir)) return false;
		char* pLastBackSlash = strrchr(PackagerDir, '\\');
		*pLastBackSlash = '\0';
		strcpy(gPackagerFilename, pLastBackSlash + 1);
		strcpy(gPathPackagerDir, PackagerDir);
		goto PackagerStuff;
	}
	char *args[2];
	args[0] = gPackagerFilename;
	args[1] = NULL;
	if(_execv(gPackagerFilename, args) < 0) return false;

	return true;
}

/////////////////////////////////////////////////////////////////////////////
bool CSetupPrepApp::GetPathToFile(char* Explan, char* Path)
{
	CSetupPrepDlg dlg;

	dlg.m_Explanation = Explan;
	if(Path != 0) dlg.m_Path = Path;
	//strcpy(dlg.m_Explanation, Explan);

	//m_pMainWnd = &dlg;
	int nResponse = dlg.DoModal();
	if(nResponse == IDOK)
	{
		strcpy(Path, dlg.m_Path);
	}
	else if(nResponse == IDCANCEL) return false;

	return true;
}

/////////////////////////////////////////////////////////////////////////////
bool CSetupPrepApp::GetPathToFolder(char* DialogTitle, char* pWhereToPut)
{
	BROWSEINFO BrowseForFolderInfo;
	BrowseForFolderInfo.hwndOwner = NULL;
	BrowseForFolderInfo.pidlRoot = NULL;
	BrowseForFolderInfo.pszDisplayName = pWhereToPut;
	BrowseForFolderInfo.lpszTitle = DialogTitle;
	BrowseForFolderInfo.ulFlags = BIF_DONTGOBELOWDOMAIN;
	BrowseForFolderInfo.lpfn = NULL;
	BrowseForFolderInfo.lParam = 0;

	ITEMIDLIST* pItemIDList = SHBrowseForFolder(&BrowseForFolderInfo);
	if(pItemIDList == NULL) return false;

	SHGetPathFromIDList(pItemIDList, pWhereToPut); 
	// In a good style, one should dismiss the *pItemIDList...

	return true;
}

/////////////////////////////////////////////////////////////////////////////
bool CSetupPrepApp::CreateSRWDirStructDescrFileAndFlatDir(char* PathSRW, char* PathSetupExe, char* PathTempDir)
{
	const char* StructFileName = gSRWStructMapFileName;
	ostringstream oStrm;

	char FullPathTempDir[MAX_PATH];
	strcpy(FullPathTempDir, PathTempDir);
	if(strcmp((FullPathTempDir + (strlen(FullPathTempDir) - 1)), "\\") != 0) strcat(FullPathTempDir, "\\");
	strcat(FullPathTempDir, gSRWFlatDirName);

	DeleteDirectoryNested(FullPathTempDir);
	if(!CreateDirectory(FullPathTempDir, NULL)) return false;

	oStrm << "dir:SRW" << gCharSepar;
	if(!TreatFilesWithinDirectory(PathSRW, oStrm, FullPathTempDir)) return false;
	oStrm << "enddir" << gCharSepar << gCharSepar;

	char OutFilePath[MAX_PATH];
	strcpy(OutFilePath, FullPathTempDir);
	strcat(OutFilePath, "\\");
	strcat(OutFilePath, StructFileName);

	CFile outFile;
	outFile.Open(OutFilePath, CFile::modeCreate | CFile::modeWrite);
	basic_string<char>& BufStr = oStrm.str();
	const char* TotStr = BufStr.c_str();
	outFile.Write(TotStr, BufStr.size()); 
	outFile.Close();

	CFile auxFile(PathSetupExe, CFile::modeRead);
	strcpy(OutFilePath, FullPathTempDir);
	strcat(OutFilePath, "\\");
	strcat(OutFilePath, auxFile.GetFileName());
	auxFile.Close();
	if(!CopyFile(PathSetupExe, OutFilePath, FALSE)) return 0;

	return true;
}

/////////////////////////////////////////////////////////////////////////////
//int CSetupApp::TreatFilesWithinDirectory(char* CurrentDir, char* ToLookFor, DWORD FileType, char* StartName, char* GlobPath, int& LenGlobPath)
int CSetupPrepApp::TreatFilesWithinDirectory(char* CurrentDir, ostringstream& ostr, char* PathTempDir)
{
	int IsFound = 0;
	WIN32_FIND_DATA FindData;
	BOOL NewDirIsSet = SetCurrentDirectory(CurrentDir);
	if(!NewDirIsSet) return 0; // not found
	//Debug
		char CurDirName[MAX_PATH];
		DWORD LenCurDirName = GetCurrentDirectory(MAX_PATH, CurDirName);
	//End Debug

	HANDLE SearchHndl;
	BOOL SomethingFound;
		
	//char NeedsFindClose = 0;
		
	SearchHndl = FindFirstFile("*\0", &FindData);
	SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);
	//NeedsFindClose = 1;

	while(SomethingFound)
	{
		IsFound = 1;
		if(*(FindData.cFileName) != '.')
		{
			if(FindData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
			{
				ostr << "dir:" << FindData.cFileName << gCharSepar;

				char LocCurDirName[MAX_PATH];
				DWORD LenLocCurDirName = GetCurrentDirectory(MAX_PATH, LocCurDirName);

				if(!TreatFilesWithinDirectory(FindData.cFileName, ostr, PathTempDir)) return 0;

				BOOL OldDirIsReset = SetCurrentDirectory(LocCurDirName);
				if(!OldDirIsReset) return 0; // not found

				ostr << "enddir" << gCharSepar;
			}
			else
			{
				ostr << FindData.cFileName << gCharSepar;

				char NewFilePath[MAX_PATH];
				strcpy(NewFilePath, PathTempDir);
				strcat(NewFilePath, "\\");
				strcat(NewFilePath, FindData.cFileName);

				if(!CopyFile(FindData.cFileName, NewFilePath, FALSE)) return 0;
			}
		}
		SomethingFound = FindNextFile(SearchHndl, &FindData);
	}

	//if(NeedsFindClose) FindClose(SearchHndl); 
	FindClose(SearchHndl); 
	return IsFound;
}

/////////////////////////////////////////////////////////////////////////////
int CSetupPrepApp::DeleteDirectoryNested(char* DirName)
{
	char LocCurDirName[1024];
	DWORD LenLocCurDirName = GetCurrentDirectory(1024, LocCurDirName);

	WIN32_FIND_DATA FindData;

	BOOL DirIsSet = SetCurrentDirectory(DirName);
	if(!DirIsSet) return 0;

	HANDLE SearchHndl = FindFirstFile("*\0", &FindData);
	BOOL SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);

	while(SomethingFound)
	{
		if(*(FindData.cFileName) != '.')
		{
			//if(FindData.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY)
			if(!(FindData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			{
				if(!DeleteFile(FindData.cFileName)) return 0;
			}
			else
			{
				if(!DeleteDirectoryNested(FindData.cFileName)) return 0;
			}
		}
		SomethingFound = FindNextFile(SearchHndl, &FindData);
	}
	FindClose(SearchHndl);

	BOOL OldDirIsReset = SetCurrentDirectory(LocCurDirName);
	if(!OldDirIsReset) return 0;

	if(!RemoveDirectory(DirName)) return 0; 
	return 1;
}
