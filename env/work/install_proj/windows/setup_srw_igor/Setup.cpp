// Setup.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Setup.h"
#include "SetupDlg.h"
#include "InstProgrDlg.h"
#include "LocateIgorDlg.h"
#include "DirDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <time.h>
#include <shlobj.h>
#include <process.h>

/////////////////////////////////////////////////////////////////////////////

char CSetupApp::gSRWStructMapFileName[] = "SRW_Temp_Directory_Locator.txt\0";
char CSetupApp::gCharSepar = '*';

/////////////////////////////////////////////////////////////////////////////
// CSetupApp

BEGIN_MESSAGE_MAP(CSetupApp, CWinApp)
	//{{AFX_MSG_MAP(CSetupApp)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG
	ON_COMMAND(ID_HELP, CWinApp::OnHelp)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSetupApp construction

CSetupApp::CSetupApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CSetupApp object

CSetupApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CSetupApp initialization

BOOL CSetupApp::InitInstance()
{
	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

#ifdef _AFXDLL
	Enable3dControls();			// Call this when using MFC in a shared DLL
#else
	//Enable3dControlsStatic();	// Call this when linking to MFC statically
#endif

	CSetupDlg dlg;
	//m_pMainWnd = &dlg;
	int nResponse = dlg.DoModal();
	if(nResponse == IDOK)
	{
		DoAutoInstallSRW();
	}
	else if(nResponse == IDCANCEL)
	{
		ReportUnsuccessfulInstallation();
	}

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.

	return FALSE;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::DoAutoInstallSRW()
{
	char PathSRWTempFolder[MAX_PATH], PathIgorProFolder[MAX_PATH];
	
	bool PrevVersionShouldBeUninstalled = true, SRW_FolderShouldBeCreated = true;
	bool ShortcutsShouldBePlaced = true;

	bool SRWTempFolderWasFound = false, IgorFolderWasFound = false, OldSRW_WasDeinstalled = false;
	bool SRW_FolderWasCreated = false, ShortCutsWerePlaced = false;

	CInstProgrDlg ProgrDlg(6);

	//gpProgrDlg = &ProgrDlg;
	//hWndProgrDlg = ProgrDlg.m_hWnd;
	
//LookForPlaceWhereSRWwasExtracted:
	
	ProgrDlg.SetPosition(1, "Looking for temporary directory ...");
	Sleep(500);

	int SRWTempFolderFound = FindSRWTempFolder(PathSRWTempFolder);
	if(!SRWTempFolderFound)
	{
		ReportUnsuccessfulInstallationSRWTempNotFound(); return 0;
	}
	else SRWTempFolderWasFound = true;
	
//LookForIgorProFolder:
	
	ProgrDlg.SetPosition(2, "Looking for Igor Pro ...");
	FindIgorProFolder(PathIgorProFolder);
	int IgorFolderFoundInfo = AskUserAboutIgorProFolder(PathIgorProFolder);
	if(IgorFolderFoundInfo == 0)
	{
		ReportUnsuccessfulInstallation(); return 0;
	}
	else if(IgorFolderFoundInfo < 0) 
	{
		PrevVersionShouldBeUninstalled = false;
		ShortcutsShouldBePlaced = false;
	}
	else IgorFolderWasFound = true;
	
UninstallOlderVersion:
	
	if(PrevVersionShouldBeUninstalled)
	{
		ProgrDlg.SetPosition(3, "Looking for / removing previous SRW version ...");
		char ShowWarningIfPrevVersionFound = 1;
		int OldSRWDeinstalled = DeinstallPreviousVersionOfSRW(PathIgorProFolder, ShowWarningIfPrevVersionFound);
		if(OldSRWDeinstalled == 0)
		{
			ReportUnsuccessfulInstallation(); return 0;
		}
		else if(OldSRWDeinstalled == -1)
		{//suggest to retry
			if(ReportCanNotUninstallSuggestRetry()) goto UninstallOlderVersion;
			else 
			{
				ReportUnsuccessfulInstallation(); return 0;
			}
		}
		else OldSRW_WasDeinstalled = true;
	}
	
CreateFolderSRW:
	
	if(SRW_FolderShouldBeCreated)
	{
		ProgrDlg.SetPosition(4, "Creating SRW directory ...");
		int SRW_DirMadeOK = MakeSRWDirectory(PathIgorProFolder, PathSRWTempFolder);
		if(!SRW_DirMadeOK)
		{
			if(ReportCanNotCreateDirectorySuggestRetry()) goto CreateFolderSRW;
			else 
			{
				ReportUnsuccessfulInstallation(); return 0;
			}
		}
		else SRW_FolderWasCreated = true;
	}
	
PlaceShortcuts:
	
	if(ShortcutsShouldBePlaced)
	{
		ProgrDlg.SetPosition(5, "Creating shortcuts ...");
		int ShortCutsArePlaced = PlaceShortCuts(PathIgorProFolder);
		if(!ShortCutsArePlaced)
		{
			if(ReportCanNotPlaceShortcutsSuggestRetry()) goto PlaceShortcuts;
			else 
			{
				ReportUnsuccessfulInstallation(); return 0;
			}
		}
		else ShortCutsWerePlaced = true;
	}

	if(SRWTempFolderWasFound && IgorFolderWasFound && SRW_FolderWasCreated && OldSRW_WasDeinstalled && ShortCutsWerePlaced)
	{
		ProgrDlg.SetPosition(6, "Installation completed");
		ReportSuccessfulInstallation();
	}
	else if(SRWTempFolderWasFound && IgorFolderWasFound && OldSRW_WasDeinstalled && SRW_FolderWasCreated &&
			(!(OldSRW_WasDeinstalled && ShortCutsWerePlaced)))
	{
		ReportSRW_FolderCreatedButInstallationWasNotCompleted();
	}

	return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Sets the PathSRWTempFolder to the path without ending "\", if succeeds
int CSetupApp::FindSRWTempFolder(char* PathSRWTempFolder)
{
	char PathTempFolder[MAX_PATH];
	//char ToLookFor[] = "SRW_Temp_Directory_Locator.txt\0";
	char* ToLookFor = gSRWStructMapFileName;

	DWORD LocatorFileType = 0;
	int SRWLocatorFound = 0;
	int LenPathSRWTempFolder = MAX_PATH;

	int LenPathTempFolder = GetTempPath(MAX_PATH, PathTempFolder);
	if(LenPathTempFolder > 0)
	{
		//EraseTerminatingBackslash(PathTempFolder);
		SRWLocatorFound = SearchWithinDirectory(PathTempFolder, ToLookFor, LocatorFileType, ToLookFor, PathSRWTempFolder, LenPathSRWTempFolder);
	}

	if(!SRWLocatorFound)
	{
		LenPathTempFolder = GetEnvironmentVariable("TMP", PathTempFolder, MAX_PATH);
		if(LenPathTempFolder > 0)
		{
			//EraseTerminatingBackslash(PathTempFolder);
			SRWLocatorFound = SearchWithinDirectory(PathTempFolder, ToLookFor, LocatorFileType, ToLookFor, PathSRWTempFolder, LenPathSRWTempFolder);
		}
	}
	if(!SRWLocatorFound)
	{
		LenPathTempFolder = GetEnvironmentVariable("TEMP", PathTempFolder, MAX_PATH);
		if(LenPathTempFolder > 0)
		{
			//EraseTerminatingBackslash(PathTempFolder);
			SRWLocatorFound = SearchWithinDirectory(PathTempFolder, ToLookFor, LocatorFileType, ToLookFor, PathSRWTempFolder, LenPathSRWTempFolder);
		}
	}
	if(!SRWLocatorFound)
	{
		SRWLocatorFound = FindPathGlobal(ToLookFor, LocatorFileType, "TMP", PathSRWTempFolder, LenPathSRWTempFolder);
	}
	if(!SRWLocatorFound)
	{
		SRWLocatorFound = FindPathGlobal(ToLookFor, LocatorFileType, "TEMP", PathSRWTempFolder, LenPathSRWTempFolder);
	}

	if(!SRWLocatorFound) return 0;

	EraseAfterLastBackslashInclusive(PathSRWTempFolder);
	return 1;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::FindIgorProFolder(char* PathIgorProFolder)
{
	char ProfileDirPath[MAX_PATH];
	char PathLinkToIgor[MAX_PATH];
	//char IgorLinkFilename[MAX_PATH];

	LPITEMIDLIST pidl;
	DWORD AuxFileType = 0;
	int AuxLenPath = MAX_PATH;

	int LinkToIgorFound = 0;
	int CSIDL_StartMenuArr[] = {CSIDL_COMMON_STARTMENU, CSIDL_STARTMENU};
	char *IgorLinkFilenames[] = {"Igor.lnk\0", "Igor Pro.lnk\0"};
	//char *pIgorLinkFilenameFound = 0;

	for(int k=0; k<2; k++)
	{
		HRESULT hres = SHGetSpecialFolderLocation(NULL, CSIDL_StartMenuArr[k], &pidl);
		if(!SUCCEEDED(hres)) continue;
		if(!SHGetPathFromIDList(pidl, ProfileDirPath)) continue; 
		// In a good style, one should dismiss the *pItemIDList...

		//strcat(ProfileDirPath, "\\Programs");

		for(int j=0; j<2; j++)
		{
			LinkToIgorFound = SearchWithinDirectory(ProfileDirPath, IgorLinkFilenames[j], AuxFileType, "Igor Pro\0", PathLinkToIgor, AuxLenPath);
			if(LinkToIgorFound) 
			{
				//pIgorLinkFilenameFound = IgorLinkFilenames[j];
				break;
			}
		}
		if(LinkToIgorFound) break;
	}

	if(LinkToIgorFound)
	{
		//*(PathLinkToIgor + (strlen(PathLinkToIgor) - 9)) = '\0';
		//if(SetCurrentDirectoryPatient(PathLinkToIgor))
		//{
		int LinkToIgorResolved = ResolveIt(NULL, PathLinkToIgor, PathIgorProFolder, false);
		if(LinkToIgorResolved)
		{
			*(PathIgorProFolder + (strlen(PathIgorProFolder) - 9)) = '\0';
			
			if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return 0;
			//strcpy(PathLinkToIgor, PathIgorProFolder);
			//GetLongPathName(PathLinkToIgor, PathIgorProFolder, MAX_PATH);
			return 1;
		}
		//}
	}
	return 0;

/**
	char ToLookFor[] = "Igor Pro Folder\0";
	char StartNameProgramFiles[] = "Program Files\0";
	//DWORD FileType = FILE_ATTRIBUTE_DIRECTORY;
	DWORD FileType = 0;

	int LenPath = MAX_PATH;
	int IgorProFolderFound = FindPathGlobal("Igor.exe\0", FileType, StartNameProgramFiles, PathIgorProFolder, LenPath);
	if(!IgorProFolderFound)
	{
		IgorProFolderFound = FindPathGlobal("igor.exe\0", FileType, StartNameProgramFiles, PathIgorProFolder, LenPath);
	}
	if(!IgorProFolderFound)
	{
		IgorProFolderFound = FindPathGlobal("IGOR.EXE\0", FileType, StartNameProgramFiles, PathIgorProFolder, LenPath);
	}
	if(!IgorProFolderFound) 
	{
		*PathIgorProFolder = '\0'; return 0;
	}

	EraseAfterLastBackslashInclusive(PathIgorProFolder);
	if(!CheckFoldersInsideIgorFolder(PathIgorProFolder))
	{
		*PathIgorProFolder = '\0'; return 0;
	}
	else return 1;
**/
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::CheckFoldersInsideIgorFolder(char* PathIgorProFolder)
{
	if(PathIgorProFolder == 0) return 0;
	if(*PathIgorProFolder == '\0') return 0;

	if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return 0;

	DWORD FileType = 0;
	char AuxPath[MAX_PATH];
	int LenAuxPath = MAX_PATH;
	char ToLookForIgorExt[] = "Igor Extensions\0";
	if(!SearchWithinDirectoryNotNested(PathIgorProFolder, ToLookForIgorExt, FileType, AuxPath, LenAuxPath)) return 0;

	char ToLookForHelpFiles[] = "Igor Help Files\0";
	LenAuxPath = MAX_PATH;
	if(!SearchWithinDirectoryNotNested(PathIgorProFolder, ToLookForHelpFiles, FileType, AuxPath, LenAuxPath)) return 0;

	char ToLookForIgorProcedures[] = "Igor Procedures\0";
	LenAuxPath = MAX_PATH;
	if(!SearchWithinDirectoryNotNested(PathIgorProFolder, ToLookForIgorProcedures, FileType, AuxPath, LenAuxPath)) return 0;

	return 1;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::AskUserAboutIgorProFolder(char* PathIgorProFolder)
{
	char AuxFileName[MAX_PATH];
	bool NameWasCopied = false;
	if(PathIgorProFolder != 0)
	{
		if(*PathIgorProFolder != '\0')
		{
			strcpy(AuxFileName, PathIgorProFolder);
			strcat(AuxFileName, "\\Igor.exe");
			NameWasCopied = true;
		}
	}
	if(!NameWasCopied)
	{
		*AuxFileName = '\0';
	}

	CLocateIgorDlg aLocateIgorDlg(AuxFileName);

LabelLocateIgor:

	int nResponse = aLocateIgorDlg.DoModal();
	if(nResponse == IDOK)
	{
		strcpy(PathIgorProFolder, aLocateIgorDlg.m_IgorPath);
	}
	else if(nResponse == IDCANCEL)
	{
		return 0;
	}

	char NameOfIgorExeFileLC[] = "igor.exe";
	int LenPathIgorProFolder = strlen(PathIgorProFolder);
	int LenNameOfIgorExeFileLC = strlen(NameOfIgorExeFileLC);
	if(LenPathIgorProFolder > LenNameOfIgorExeFileLC)
	{
		char TestSubStr[MAX_PATH];
		char* pSubStr = PathIgorProFolder + LenPathIgorProFolder - LenNameOfIgorExeFileLC;
		strcpy(TestSubStr, pSubStr);
		_strlwr(TestSubStr); // to lower case a-la MS
		
		if(strcmp(TestSubStr, NameOfIgorExeFileLC) == 0)
		{
			EraseAfterLastBackslashInclusive(PathIgorProFolder);
		}
		else if(*(PathIgorProFolder + LenPathIgorProFolder - 1) == '\\')
		{
			EraseTerminatingBackslash(PathIgorProFolder);
		}
	}

	int IgorIsIdentified = CheckFoldersInsideIgorFolder(PathIgorProFolder);
	if(!IgorIsIdentified)
	{
		EraseAfterLastBackslashInclusive(PathIgorProFolder);
		IgorIsIdentified = CheckFoldersInsideIgorFolder(PathIgorProFolder);
	}

	if(!IgorIsIdentified)
	{
		int Res = ReportIncorrectIgorFolderAndAskToTryAgain();
		if(Res == 0) return 0; // cancel installation
		else if(Res < 0) goto LabelLocateIgor; // try once more
		else return -1; // copy SRW without further installation
	}
	return 1;
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::EraseTerminatingBackslash(char* str)
{
	if(str == NULL) return;

	int len = strlen(str);
	if(len <= 0) return;

	char* pLastBackSlash = strrchr(str, '\\');
	if(pLastBackSlash != NULL)
	{
		if((str + (len - 1)) == pLastBackSlash) *pLastBackSlash = '\0';
	}
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::EraseAfterLastBackslashInclusive(char* str)
{
	if(str == NULL) return;

	int len = strlen(str);
	if(len <= 0) return;

	char* pLastBackSlash = strrchr(str, '\\');
	if(pLastBackSlash != NULL)
	{
		*pLastBackSlash = '\0';
	}
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::FindPathGlobal(char* ToLookFor, DWORD FileType, char* StartName, char* GlobPath, int& LenGlobPath)
{
	DWORD nBufferLength = 512;
	char DriveNamesString[512];
	DWORD LenDriveNames = GetLogicalDriveStrings(nBufferLength, DriveNamesString); 
  
	char StopDrivesLoop = 0;
	char *tDriveNames = DriveNamesString;
	char *CurrentDriveName;
	while(!StopDrivesLoop)
	{
		CurrentDriveName = tDriveNames;
		if(*CurrentDriveName == '\0') break;

		ArbSearchStartedG = (*StartName == '\0')? 1 : 0;

		char CurrentDir[1024];
		UINT DriveType = GetDriveType(CurrentDriveName);
		if((DriveType == DRIVE_FIXED) || (DriveType == DRIVE_REMOTE) || (DriveType == DRIVE_RAMDISK))
		{
			strcpy(CurrentDir, CurrentDriveName);
			
			int IsFound = SearchWithinDirectory(CurrentDir, ToLookFor, FileType, StartName, GlobPath, LenGlobPath);
			if(IsFound) return 1; //found
		}

		while(*tDriveNames != '\0') tDriveNames++;
		tDriveNames++;
	}
	return 0; // not found
}

/////////////////////////////////////////////////////////////////////////////

BOOL CSetupApp::CheckIfInstallationShouldBeCancelled()
{
//	if(gpProgrDlg == 0) return FALSE;
//	MSG aMes;
	//aMes.message = BN_CLICKED;

	//gpProgrDlg->RunModalLoop(MLF_NOKICKIDLE);
	//UpdateWindow(hWndProgrDlg);
    
	//BOOL PeekRes = GetMessage(&aMes, hWndProgrDlg, 0, 0);
	//BOOL PeekRes = PeekMessage(&aMes, hWndProgrDlg, WM_COMMAND, WM_COMMAND, PM_NOREMOVE);
	//BOOL PeekRes = PeekMessage(&aMes, gpProgrDlg->m_hWnd, 0, 0, PM_NOREMOVE);

/*
	//BOOL PeekRes = PeekMessage(&aMes, NULL, 0, 0, PM_NOREMOVE);
//  UINT wMsgFilterMin,  // first message
//  UINT wMsgFilterMax,  // last message
*/
/*
	if(PeekRes)
	{
		PreTranslateMessage(&aMes);
		DispatchMessage(&aMes);

    //LONG lIdle = 0;
    //while ( AfxGetApp()->OnIdle(lIdle++ ) ) 
	//{
	//
	//	int Aha = 1;
	//}
		//HANDLE tstH = GetCurrentProcess();

		//DWORD test = WaitForInputIdle(GetCurrentProcess(), 500);
	}
	//UpdateWindow(hWndProgrDlg);
*/
	return 0;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::SearchWithinDirectory(char* CurrentDir, char* ToLookFor, DWORD FileType, char* StartName, char* GlobPath, int& LenGlobPath)
{
	WIN32_FIND_DATA FindData;

	BOOL NewDirIsSet = SetCurrentDirectoryPatient(CurrentDir);
	if(!NewDirIsSet) return 0; // not found
	int IsFound = 0;
	//Debug
	//	char CurDirName[1024];
	//	DWORD LenCurDirName = GetCurrentDirectory(1024, CurDirName);
	//End Debug
	
	HANDLE SearchHndl;
	BOOL SomethingFound;

	if(CheckIfInstallationShouldBeCancelled()) return 0;

	char NeedsFindClose = 0;
		
	if(!ArbSearchStartedG) 
	{
		SearchHndl = FindFirstFile(StartName, &FindData);

		SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);
		FindClose(SearchHndl); 
		ArbSearchStartedG = 1;

		if(!SomethingFound)
		{
			SearchHndl = FindFirstFile("*\0", &FindData);
			SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);
			NeedsFindClose = 1;
		}
	}
	else
	{
		SearchHndl = FindFirstFile("*\0", &FindData);
		SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);
		NeedsFindClose = 1;
	}

	while(SomethingFound)
	{
		if(*(FindData.cFileName) != '.')
		{
			//if((FileType == 0) || (FindData.dwFileAttributes == FileType))
			if((FileType == 0) || (FindData.dwFileAttributes & FileType))
			{
				if(!strcmp(FindData.cFileName, ToLookFor))
				{
					IsFound = 1; // found
					char **lpFilePartDummy = 0;
					LenGlobPath = GetFullPathName(ToLookFor, LenGlobPath, GlobPath, lpFilePartDummy);
					break; 
				}
			}
			//if((FindData.dwFileAttributes == FILE_ATTRIBUTE_DIRECTORY) || (FindData.dwFileAttributes == 17))
			if(FindData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
			{
				char LocCurDirName[1024];
				DWORD LenLocCurDirName = GetCurrentDirectory(1024, LocCurDirName);

				IsFound = SearchWithinDirectory(FindData.cFileName, ToLookFor, FileType, StartName, GlobPath, LenGlobPath);
				if(IsFound) break; // found

				BOOL OldDirIsReset = SetCurrentDirectoryPatient(LocCurDirName);
				if(!OldDirIsReset) return 0; // not found
			}
		}
		if(ArbSearchStartedG && NeedsFindClose)
		{
			SomethingFound = FindNextFile(SearchHndl, &FindData);
		}
		else
		{
			SearchHndl = FindFirstFile("*\0", &FindData);
			SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);
			ArbSearchStartedG = 1;
			NeedsFindClose = 1;
		}
	}

	if(NeedsFindClose) FindClose(SearchHndl); 
 	return IsFound; // not found
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::SearchWithinDirectoryNotNested(char* CurrentDir, char* ToLookFor, DWORD FileType, char* GlobPath, int& LenGlobPath)
{
	WIN32_FIND_DATA FindData;

	BOOL NewDirIsSet = SetCurrentDirectoryPatient(CurrentDir);
	if(!NewDirIsSet) return 0; // not found

	int IsFound = 0;
	
	HANDLE SearchHndl = FindFirstFile("*\0", &FindData);
	BOOL SomethingFound = (SearchHndl != INVALID_HANDLE_VALUE);

	while(SomethingFound)
	{
		if(*(FindData.cFileName) != '.')
		{
			//if((FileType == 0) || (FindData.dwFileAttributes == FileType))
			if((FileType == 0) || (FindData.dwFileAttributes & FileType))
			{
				if(!strcmp(FindData.cFileName, ToLookFor))
				{
					IsFound = 1; // found
					char **lpFilePartDummy = 0;
					LenGlobPath = GetFullPathName(ToLookFor, LenGlobPath, GlobPath, lpFilePartDummy);
					break; 
				}
			}
		}
		SomethingFound = FindNextFile(SearchHndl, &FindData); 
	}

	FindClose(SearchHndl); 
 	return IsFound; // not found
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::RemoveDirectoryNamedIfExists(const char* LocName, const char* WhereItIs)
{
	char LocCurDirName[1024];
	DWORD LenLocCurDirName = GetCurrentDirectory(1024, LocCurDirName);

 	BOOL WhereItIsOK = SetCurrentDirectoryPatient(WhereItIs);
	if(!WhereItIsOK) return 0;

	WIN32_FIND_DATA FindData;
	HANDLE SearchHndl = FindFirstFile(LocName, &FindData);
	BOOL IsFound = (SearchHndl != INVALID_HANDLE_VALUE);
	FindClose(SearchHndl); 
	if(IsFound) 
	{
		if(!DeleteDirectoryNested(FindData.cFileName)) return 0;
	}

	BOOL OldDirIsReset = SetCurrentDirectoryPatient(LocCurDirName);
	if(!OldDirIsReset) return 0; // not found

	return 1;
}

/////////////////////////////////////////////////////////////////////////////
// returns: 0- failed, cancel further installation; -1- failed, suggest to retry; 1- succeeded
int CSetupApp::DeinstallPreviousVersionOfSRW(char* PathIgorProFolder, char ReportOtherVersionFound)
{
	WIN32_FIND_DATA FindData;

// SRW folder
 	BOOL IgorDirIsSet = SetCurrentDirectoryPatient(PathIgorProFolder);
	if(!IgorDirIsSet) return 0;
	HANDLE SRW_SearchHndl = FindFirstFile("SRW\0", &FindData);
	BOOL SRW_Found = (SRW_SearchHndl != INVALID_HANDLE_VALUE);
	FindClose(SRW_SearchHndl); 
	if(SRW_Found && ReportOtherVersionFound) 
	{
		char AuxPath[MAX_PATH];
		strcpy(AuxPath, PathIgorProFolder);
		strcat(AuxPath, "\\SRW");
		if(!ReportConfirmDeletePrevVersionSRW(AuxPath)) return 0; //cancel installation
	}

// Shortcut to SRW.xop
	if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return -1;
	if(!SetCurrentDirectoryPatient("Igor Extensions\0")) return -1;
	RemoveShortcutsToSRWFilesFromCurrentDir();

// Shortcut to SRW Macro Help.ifn
	if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return -1;
	if(!SetCurrentDirectoryPatient("Igor Help Files\0")) return -1;
	RemoveShortcutsToSRWFilesFromCurrentDir();

// Shortcut to SRW Procedures
	if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return -1;
	if(!SetCurrentDirectoryPatient("Igor Procedures\0")) return -1;
	RemoveShortcutsToSRWFilesFromCurrentDir();

	if(SRW_Found) 
	{// deleting SRW directory after all shortcuts !!!
		if(!SetCurrentDirectoryPatient(PathIgorProFolder)) return 0;
		SRW_SearchHndl = FindFirstFile("SRW\0", &FindData);
		if(SRW_SearchHndl != INVALID_HANDLE_VALUE)
		{
			if(!DeleteDirectoryNested(FindData.cFileName)) return -1;
		}
		FindClose(SRW_SearchHndl); 
	}
	return 1;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::RemoveShortcutsToSRWFilesFromCurrentDir()
{
    CoInitialize(NULL);

	char PathOrigFile[MAX_PATH];
	WIN32_FIND_DATA FindData;
	HANDLE aHndl = FindFirstFile("*", &FindData);
	BOOL SomethingFound = (aHndl != INVALID_HANDLE_VALUE);
	while(SomethingFound)
	{
		int DeletionIsNeeded = 0;
		if(ResolveIt(NULL, FindData.cFileName, PathOrigFile, true) && (*PathOrigFile != '\0'))
		{
			char* pLastBackSlash = strrchr(PathOrigFile, '\\');
			if(*(pLastBackSlash + 1) != '\0')
			{
				char* pFileNameToCheck = pLastBackSlash + 1;
				DeletionIsNeeded = CheckIfNameIsPointedBySRWShortcut(pFileNameToCheck);
			}
		}
		if(DeletionIsNeeded) 
		{
			if(!DeleteFilePatient(FindData.cFileName)) return -1;
		}
		SomethingFound = FindNextFile(aHndl, &FindData);
	}

	CoUninitialize();
	return 1;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::CheckIfNameIsPointedBySRWShortcut(char* FileNameToCheck)
{
	if((FileNameToCheck == 0) || (*FileNameToCheck == '\0')) return 0;
	
	char* LinkedNamesSRW_UC[] = {
		"SRW.XOP", "SRW PROCEDURES", "SRW MACRO HELP.IFN", "SRW MACRO HELP.IHF",
		"SRWPRO~1", "SRWMAC~1.IFN", "SRWMAC~1.IHF"
	};
	int AmOfLinkedNames = 7;

	char AuxName[MAX_PATH];
	strcpy(AuxName, FileNameToCheck);
	_strupr(AuxName); // toUpperCase 

	for(int i=0; i<AmOfLinkedNames; i++)
	{
		if(strcmp(AuxName, LinkedNamesSRW_UC[i]) == 0) return 1;
	}
	return 0;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::DeleteDirectoryNested(char* DirName)
{
	char LocCurDirName[1024];
	DWORD LenLocCurDirName = GetCurrentDirectory(1024, LocCurDirName);

	WIN32_FIND_DATA FindData;

	BOOL DirIsSet = SetCurrentDirectoryPatient(DirName);
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
				if(!DeleteFilePatient(FindData.cFileName)) return 0;
			}
			else
			{
				if(!DeleteDirectoryNested(FindData.cFileName)) return 0;
			}
		}
		SomethingFound = FindNextFile(SearchHndl, &FindData);
	}
	FindClose(SearchHndl);

	BOOL OldDirIsReset = SetCurrentDirectoryPatient(LocCurDirName);
	if(!OldDirIsReset) return 0;

	if(!RemoveDirectoryPatient(DirName)) return 0; 
	return 1;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::ReportIncorrectIgorFolderAndAskToTryAgain()
{
	char ErrorMesTitle[] = "Igor Pro Folder not identified";
	char ErrorStr[] = "Igor Pro Folder was not identified. If SRW will be placed into the folder\nyou have specified, it will not function properly.\n\nAre you sure you want to place SRW into this folder?\n\nPress \"Yes\" if so. Otherwise, press \"No\" to specify another folder,\nor \"Cancel\" to stop the installation.";
	
	UINT DlgStyle = MB_YESNOCANCEL | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int TestRes = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);

	if(TestRes == IDYES) return 1;
	else if(TestRes == IDNO) return -1;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportUnsuccessfulInstallation()
{
	char ErrorMesTitle[] = "SRW was not installed";
	char ErrorStr[] = "SRW was not installed on your system.\nYou can try the installation later.";

	UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportUnsuccessfulInstallationSRWTempNotFound()
{
	char ErrorMesTitle[] = "SRW Installation Error";
	char ErrorStr[] = "SRW installation failed. One possible reason can be \nthat no environmental variable (TEMP, TMP) describing \nthe path to a temporary directory is defined on your system. \nIn case of a serious problem, please contact: chubar@esrf.fr.";
	UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::ReportConfirmDeletePrevVersionSRW(char* PathPrevVersion)
{
	char ErrorMesTitle[] = "SRW Installation Warning";
	
	char ErrorStr[2048];
	char ErrorMainPartStr[] = "Another version of SRW was found.\nThis version will be deleted. Make sure that no your working files \nare located in: ";
	strcpy(ErrorStr, ErrorMainPartStr);

	if((PathPrevVersion == 0) || (*PathPrevVersion == 0))
	{
		strcat(ErrorStr, "\"...\\Igor Pro Folder\\SRW\".");
	}
	else 
	{
		strcat(ErrorStr, "\"");
		strcat(ErrorStr, PathPrevVersion);
		strcat(ErrorStr, "\".");
	}

	strcat(ErrorStr, "\n\nPress \"OK\" to proceed removing the other version and installing \nthe current one, or \"Cancel\" to stop the installation.");
	
	UINT DlgStyle = MB_OKCANCEL | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
	if(MesBoxInf == IDOK) return 1;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::ReportCanNotUninstallSuggestRetry()
{
	char ErrorStr[] = "Older version of SRW can not be uninstalled.\nMake sure you have closed Igor Pro and\nany other applications that may use SRW files.";
	char ErrorMesTitle[] = "SRW Installation Error";
	
	UINT DlgStyle = MB_RETRYCANCEL | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
	
	if(MesBoxInf == IDRETRY) return 1;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::ReportCanNotCreateDirectorySuggestRetry()
{
	char ErrorMesTitle[] = "SRW Installation Error";
	char ErrorStr[] = "Can not create directory structure\nfor the new version of SRW.";
	
	UINT DlgStyle = MB_RETRYCANCEL | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
	if(MesBoxInf == IDRETRY) return 1;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::ReportCanNotPlaceShortcutsSuggestRetry()
{
	char ErrorMesTitle[] = "SRW Installation Error";
	char ErrorStr[] = "Can not place shortcuts to SRW files\nin Igor Pro folder.";

	UINT DlgStyle = MB_RETRYCANCEL | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
	if(MesBoxInf == IDRETRY) return 1;
	else return 0;
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportUnsuccessfulInstallationIgorNotFound()
{
	char ErrorStr[] = "SRW can not be installed correctly.\nNo Igor Pro 3.1x or greater was found.\n";
	char ErrorMesTitle[] = "SRW Installation Error";
	
	UINT DlgStyle = MB_OK | MB_ICONWARNING | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
	
	ReportUnsuccessfulInstallation();
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportSuccessfulInstallation()
{
	char MesTitle[] = "SRW Installation Completed";
	char Str[] = "SRW was successfully installed.\n\nTry to start Igor Pro now. At the first start,\nyou may be asked to compile SRW help file.\nChoose to compile it. Then read the \"Introduction\"\nand try examples from the menu SRWE or SRWP.\nGood luck!";

	UINT DlgStyle = MB_OK | MB_ICONASTERISK | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, Str, MesTitle, DlgStyle); 
}

/////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportSRW_FolderCreatedButInstallationWasNotCompleted()
{
	char MesTitle[] = "SRW Was Copied";
	char Str[] = "SRW directory was copied on your hard drive.\n\nHowever, the installation was not completed.\nTherefore the SRW may not function properly.\n\nYou can try the installation later.";

	UINT DlgStyle = MB_OK | MB_ICONASTERISK | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, Str, MesTitle, DlgStyle); 
}

/////////////////////////////////////////////////////////////////////////////

int CSetupApp::MakeSRWDirectory(char* PathIgorProFolder, char* PathSRWTempFolder)
{
//Debug, to remove!!!
//	strcpy(StartDirNameG, "C:\\TEMP\\SRW_InstallTest");
//End Debug
	//ostringstream o;
	//SRWDirectoryStructure(o);
	//basic_string<char>& BufStr = o.str();
	//const char* TotStr = BufStr.c_str();

	char* StructFileName = gSRWStructMapFileName;
	char StructFilePath[MAX_PATH];
	strcpy(StructFilePath, PathSRWTempFolder);
	if(strcmp((PathSRWTempFolder + (strlen(PathSRWTempFolder) - 1)), "\\") != 0) strcat(StructFilePath, "\\");
	strcat(StructFilePath, gSRWStructMapFileName);

	char TotStr[10000];
	CFile inFile(StructFilePath, CFile::modeRead);
	int AmOfBytes = inFile.Read(TotStr, 10000);
	inFile.Close();
	if(AmOfBytes <= 0) return 0;

	char* tAux = TotStr;
	for(int k=0; k<AmOfBytes; k++) 
	{
		if((*tAux) == gCharSepar) 
		{
			*tAux = '\0';
		}
		tAux++;
	}

	const char* tTotStr = TotStr;

	BOOL IgorDirIsSet = SetCurrentDirectoryPatient(PathIgorProFolder);
	if(!IgorDirIsSet) return 0;

	char DirStart = ((tTotStr[0] == 'd') && (tTotStr[1] == 'i') && (tTotStr[2] == 'r') && (tTotStr[3] == ':'));
	if(DirStart)
	{
		if(!RemoveDirectoryNamedIfExists(tTotStr + 4, PathIgorProFolder)) return 0;
	}

	SECURITY_ATTRIBUTES SecAttr;
	SecAttr.nLength = sizeof(SECURITY_ATTRIBUTES);
	SecAttr.lpSecurityDescriptor = NULL;
	SecAttr.bInheritHandle = FALSE;

	TCharPtrVect NestingDirNames;
	DWORD LenNameOfNestingDir;

	char StopMainLoop = 0;
	while(!StopMainLoop)
	{
		if(*tTotStr == '\0') break;

		char NewDirStart = ((tTotStr[0] == 'd') && (tTotStr[1] == 'i') && (tTotStr[2] == 'r') && (tTotStr[3] == ':'));
		if(NewDirStart)
		{
			tTotStr += 4;

			char* CurNestingDirName = new char[1024];
			LenNameOfNestingDir = GetCurrentDirectory(1024, CurNestingDirName);
			NestingDirNames.push_back(CurNestingDirName);

			BOOL DirCreatedOK = CreateDirectoryPatient(tTotStr, &SecAttr);
			if(!DirCreatedOK) return 0;
			else
			{
				BOOL DirIsSet = SetCurrentDirectoryPatient(tTotStr);
				if(!DirIsSet) return 0;
			}
 		}
		else if(!strcmp(tTotStr, "enddir"))
		{
			char* CurNestingDirName = *(NestingDirNames.end() - 1);
			BOOL DirIsSet = SetCurrentDirectoryPatient(CurNestingDirName);
			if(!DirIsSet) return 0;
			NestingDirNames.pop_back();
			delete[] CurNestingDirName;
		}
		else
		{
			char ExistingFilePath[1024];
			//strcpy(ExistingFilePath, StartDirNameG);
			strcpy(ExistingFilePath, PathSRWTempFolder);
			int CurStrLen = strlen(ExistingFilePath);
			if(*(ExistingFilePath + CurStrLen - 1) != '\\') strcat(ExistingFilePath, "\\");
			strcat(ExistingFilePath, tTotStr);

			BOOL FileCopied = CopyFilePatient(ExistingFilePath, tTotStr, FALSE);
			if(!FileCopied) return 0;
		}

		while(*tTotStr != '\0') tTotStr++;
		tTotStr++;
	}

	return 1;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::PlaceShortCuts(char* PathIgorProFolder)
{
    CoInitialize(NULL);

	BOOL DirIsSet = SetCurrentDirectoryPatient(PathIgorProFolder);
	if(!DirIsSet) return 0;
	DirIsSet = SetCurrentDirectoryPatient("SRW\0");
	if(!DirIsSet) return 0;

	char PathOrigFile[1024], PathShortCut[1024], NewPathShortCut[1024];

// SRW.xop
	strcpy(PathOrigFile, PathIgorProFolder);
	int CurStrLen = strlen(PathOrigFile);
	if(*(PathOrigFile + CurStrLen - 1) != '\\') strcat(PathOrigFile, "\\");
	strcat(PathOrigFile, "SRW\\SRW.xop");

	int ShortCutDone = MakeShortCutTo(PathOrigFile, PathShortCut);
	if(!ShortCutDone) return 0;

	strcpy(NewPathShortCut, PathIgorProFolder);
	CurStrLen = strlen(NewPathShortCut);
	if(*(NewPathShortCut + CurStrLen - 1) != '\\') strcat(NewPathShortCut, "\\");
	strcat(NewPathShortCut, "Igor Extensions");
	char* pShortCutNameBackSlashPrep = strrchr(PathShortCut, '\\');
	strcat(NewPathShortCut, pShortCutNameBackSlashPrep);

	BOOL MovedOK = MoveFile(PathShortCut, NewPathShortCut); 
	if(!MovedOK) return 0;

// SRW Proc
	strcpy(PathOrigFile, PathIgorProFolder);
	CurStrLen = strlen(PathOrigFile);
	if(*(PathOrigFile + CurStrLen - 1) != '\\') strcat(PathOrigFile, "\\");
	strcat(PathOrigFile, "SRW\\SRW Procedures");

	ShortCutDone = MakeShortCutTo(PathOrigFile, PathShortCut);
	if(!ShortCutDone) return 0;

	strcpy(NewPathShortCut, PathIgorProFolder);
	CurStrLen = strlen(NewPathShortCut);
	if(*(NewPathShortCut + CurStrLen - 1) != '\\') strcat(NewPathShortCut, "\\");
	strcat(NewPathShortCut, "Igor Procedures");
	pShortCutNameBackSlashPrep = strrchr(PathShortCut, '\\');
	strcat(NewPathShortCut, pShortCutNameBackSlashPrep);

	MovedOK = MoveFile(PathShortCut, NewPathShortCut); 
	if(!MovedOK) return 0;

// SRW Macro Help
	DirIsSet = SetCurrentDirectoryPatient("SRW Help\0");
	if(!DirIsSet) return 0;

	strcpy(PathOrigFile, PathIgorProFolder);
	CurStrLen = strlen(PathOrigFile);
	if(*(PathOrigFile + CurStrLen - 1) != '\\') strcat(PathOrigFile, "\\");
	//strcat(PathOrigFile, "SRW\\SRW Help\\SRW Macro Help.ifn");
	strcat(PathOrigFile, "SRW\\SRW Help\\SRW Macro Help.ihf");

	ShortCutDone = MakeShortCutTo(PathOrigFile, PathShortCut);
	if(!ShortCutDone) return 0;

	strcpy(NewPathShortCut, PathIgorProFolder);
	CurStrLen = strlen(NewPathShortCut);
	if(*(NewPathShortCut + CurStrLen - 1) != '\\') strcat(NewPathShortCut, "\\");
	strcat(NewPathShortCut, "Igor Help Files");
	pShortCutNameBackSlashPrep = strrchr(PathShortCut, '\\');
	strcat(NewPathShortCut, pShortCutNameBackSlashPrep);

	MovedOK = MoveFile(PathShortCut, NewPathShortCut); 
	if(!MovedOK) return 0;

	CoUninitialize();
	return 1;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::MakeShortCutTo(char* PathOrigFile, char* NameShortCut)
{
	WIN32_FIND_DATA FindData;
	HANDLE SearchHndl = FindFirstFile(PathOrigFile, &FindData);
	BOOL Found = (SearchHndl != INVALID_HANDLE_VALUE);
	if(!Found) return 0;
	FindClose(SearchHndl); 
	//char* pszShortcutFile = FindData.cFileName;

	char pszShortcutFile[1024];
	GetCurrentDirectory(1024, pszShortcutFile);
	strcat(pszShortcutFile, "\\");
	strcat(pszShortcutFile, FindData.cFileName);

    IShellLink *psl;

	// Create an IShellLink object and get a pointer to the IShellLink
	// interface (returned from CoCreateInstance).
    HRESULT hres = CoCreateInstance(CLSID_ShellLink, NULL, CLSCTX_INPROC_SERVER, IID_IShellLink, (void**)&psl);
    if(SUCCEEDED(hres))
    {
		IPersistFile *ppf;
		// Query IShellLink for the IPersistFile interface for
		// saving the shortcut in persistent storage.
		hres = psl->QueryInterface(IID_IPersistFile, (void**)&ppf);
		if(SUCCEEDED(hres))
		{
			//WORD wsz[MAX_PATH];   // buffer for Unicode string
			WCHAR wsz[MAX_PATH];   // buffer for Unicode string
			//OC310806

			// Set the path to the shortcut target.
			hres = psl->SetPath(pszShortcutFile);
			//if(!SUCCEEDED(hres)) AfxMessageBox("SetPath failed!");
			if(!SUCCEEDED(hres)) return 0;

			char pszDesc[512], pszLink[1024];
			strcpy(pszDesc, "Shortcut to ");
			strcat(pszDesc, FindData.cFileName);

			GetCurrentDirectory(1024, pszLink);
			int CurStrLen = strlen(pszLink);
			if(*(pszLink + CurStrLen - 1) != '\\') strcat(pszLink, "\\");
			strcat(pszLink, pszDesc);
			strcat(pszLink, ".LNK");
			// Ensure that the string consists of ANSI characters.
			MultiByteToWideChar(CP_ACP, 0, pszLink, -1, wsz, MAX_PATH);

            // Set the description of the shortcut.
			hres = psl->SetDescription(pszDesc);
			//if(!SUCCEEDED(hres)) AfxMessageBox("SetDescription failed!");
			if(!SUCCEEDED(hres)) return 0;

			// Save the shortcut via the IPersistFile::Save member function.
			hres = ppf->Save(wsz, TRUE);
			//if(!SUCCEEDED(hres)) AfxMessageBox("Save failed!");
			if(!SUCCEEDED(hres)) return 0;

			strcpy(NameShortCut, pszLink);

            // Release the pointer to IPersistFile.
			ppf->Release();
		}
       // Release the pointer to IShellLink.
       psl->Release();
    }
	return 1;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::ResolveIt(HWND hwnd, LPCSTR lpszLinkFile, LPSTR lpszPath, bool GetShortPath) 
{ 
    HRESULT hres; 
    IShellLink* psl; 
    char szGotPath[MAX_PATH]; 
    char szDescription[MAX_PATH]; 
    WIN32_FIND_DATA wfd; 
 
    *lpszPath = 0; // assume failure 
 
	CoInitialize(NULL);

    // Get a pointer to the IShellLink interface. 
    hres = CoCreateInstance(CLSID_ShellLink, NULL, CLSCTX_INPROC_SERVER, IID_IShellLink, (LPVOID *) &psl); 
    if(SUCCEEDED(hres)) 
	{ 
        IPersistFile* ppf; 
 
        // Get a pointer to the IPersistFile interface. 
        hres = psl->QueryInterface(IID_IPersistFile, (void**)&ppf); 
        if(SUCCEEDED(hres)) 
		{ 
            WCHAR wsz[MAX_PATH]; 
 
            // Ensure that the string is Unicode. 
            MultiByteToWideChar(CP_ACP, 0, lpszLinkFile, -1, wsz, MAX_PATH); 
 
            // Load the shortcut. 
            hres = ppf->Load(wsz, STGM_READ); 
            if(SUCCEEDED(hres)) 
			{ 
                // Resolve the link. 
                //hres = psl->Resolve(hwnd, 0); 
                hres = psl->Resolve(hwnd, SLR_NO_UI); 
                if(SUCCEEDED(hres)) 
				{ 
                    // Get the path to the link target. 
					if(GetShortPath)
						hres = psl->GetPath(szGotPath, MAX_PATH, (WIN32_FIND_DATA *)&wfd, SLGP_SHORTPATH ); 
					else hres = psl->GetPath(szGotPath, MAX_PATH, (WIN32_FIND_DATA *)&wfd, SLGP_UNCPRIORITY ); 
                    if(FAILED(hres)) { CoUninitialize(); return 0;}
 
                    // Get the description of the target. 
                    hres = psl->GetDescription(szDescription, MAX_PATH); 
                    if (FAILED(hres)) { CoUninitialize(); return 0;}
                    lstrcpy(lpszPath, szGotPath); 
                }
				else { CoUninitialize(); return 0;}
            }
			else { CoUninitialize(); return 0;}
        // Release the pointer to the IPersistFile interface. 
			ppf->Release(); 
        }
		else { CoUninitialize(); return 0;}
    // Release the pointer to the IShellLink interface. 
		psl->Release(); 
    }
	else { CoUninitialize(); return 0;}
	CoUninitialize();
    return 1; 
} 

////////////////////////////////////////////////////////////////////////////

int CSetupApp::DeleteFilePatient(char* FileName)
{
	const double TrendTime = 5.; //s
	clock_t Delta = (clock_t)(CLOCKS_PER_SEC*TrendTime);
	clock_t EndTrendTime = clock() + Delta;
	int FileDelOK = 0;
	while(clock() < EndTrendTime)
	{
		FileDelOK = DeleteFile(FileName);
		if(FileDelOK) break;
	}
	return FileDelOK;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::RemoveDirectoryPatient(char* DirName)
{
	const double TrendTime = 5.; //s
	clock_t Delta = (clock_t)(CLOCKS_PER_SEC*TrendTime);
	clock_t EndTrendTime = clock() + Delta;
	int DirDelOK = 0;
	while(clock() < EndTrendTime)
	{
		DirDelOK = RemoveDirectory(DirName);
		if(DirDelOK) break;
	}
	return DirDelOK;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::CreateDirectoryPatient(const char* DirName, SECURITY_ATTRIBUTES* pSecAttr)
{
	const double TrendTime = 5.; //s
	clock_t Delta = (clock_t)(CLOCKS_PER_SEC*TrendTime);
	clock_t EndTrendTime = clock() + Delta;
	int DirCrOK = 0;
	while(clock() < EndTrendTime)
	{
		DirCrOK = CreateDirectory(DirName, pSecAttr);
		if(DirCrOK) break;

		//Debug
		//if(!DirCrOK)
		//{
		//	int Aga = 1;
		//}
		//EndDebug
	}
	return DirCrOK;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::SetCurrentDirectoryPatient(const char* DirName)
{
	const double TrendTime = 5.; //s
	clock_t Delta = (clock_t)(CLOCKS_PER_SEC*TrendTime);
	clock_t EndTrendTime = clock() + Delta;
	int DirCrOK = 0;
	while(clock() < EndTrendTime)
	{
		DirCrOK = SetCurrentDirectory(DirName);
		if(DirCrOK) break;
	}
	return DirCrOK;
}

////////////////////////////////////////////////////////////////////////////

int CSetupApp::CopyFilePatient(const char* ExistingFilePath, const char* NewFilePath, BOOL Par)
{
	const double TrendTime = 5.; //s
	clock_t Delta = (clock_t)(CLOCKS_PER_SEC*TrendTime);
	clock_t EndTrendTime = clock() + Delta;
	int CpOK = 0;
	while(clock() < EndTrendTime)
	{
		CpOK = CopyFile(ExistingFilePath, NewFilePath, Par);
		if(CpOK) break;
	}
	return CpOK;
}

////////////////////////////////////////////////////////////////////////////
/**
int CSetupApp::DoManInstallSRW()
{
	//LenStartDirNameG = GetCurrentDirectory(1024, StartDirNameG);
	char ToLookForTemp[] = "SRW Temp Directory Locator.txt\0";
	//DWORD LocatorFileType = FILE_ATTRIBUTE_NORMAL;
	DWORD LocatorFileType = 0;
	char StartName[] = "TEMP\0";
	int LenSRWTempPath = 1024;
	int SRWTempFolderFound = FindPathGlobal(ToLookForTemp, LocatorFileType, StartName, StartDirNameG, LenSRWTempPath);
	if(!SRWTempFolderFound)
	{
		char AnotherStartName[] = "TMP\0";
		LenSRWTempPath = 1024;
		SRWTempFolderFound = FindPathGlobal(ToLookForTemp, LocatorFileType, AnotherStartName, StartDirNameG, LenSRWTempPath);
	}
	if(!SRWTempFolderFound)
	{
		char ErrorMesTitle[] = "SRW Installation Error";
		char ErrorStr[] = "SRW installation failed. \nPlease try another SRW distribution pack\nor send e-mail to: chubar@esrf.fr";
		
		UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
		int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle);
		return 0;
	}
	LenStartDirNameG = LenSRWTempPath;
	char* pLastBackSlash = strrchr(StartDirNameG, '\\');
	*pLastBackSlash = '\0';

	char TitleString[] = "Select a folder where to put SRW (try to find \"Igor Pro Folder\")";
	char WhereToPutSRW[1024];
	if(!GetFolderFromUser(TitleString, WhereToPutSRW))
	{
		ReportUnsuccessfulInstallation(); return 0;
	}

	int SRW_DirMade = MakeSRWDirectory(WhereToPutSRW);
	if(!SRW_DirMade)
	{
		ReportUnsuccessfulInstallation(); return 0;
	}

	ReportDeliveryAndNeedCompletingTheInstallation(WhereToPutSRW);

	return 1;
}
**/
/////////////////////////////////////////////////////////////////////////////

int CSetupApp::GetFolderFromUser(char* DialogTitle, char* pWhereToPut)
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
	if(pItemIDList == NULL) return 0;

	SHGetPathFromIDList(pItemIDList, pWhereToPut); 
	// In a good style, one should dismiss the *pItemIDList...

	return 1;
}

////////////////////////////////////////////////////////////////////////////

void CSetupApp::ReportDeliveryAndNeedCompletingTheInstallation(char* WhereHePutIt)
{
	char MesTitle[] = "SRW Delivered";

	char PathReadMe[1024];
	strcpy(PathReadMe, WhereHePutIt);
	strcat(PathReadMe, "\\SRW\\ReadMe.txt");

	char StrTot[1024];
	char Str1[] = "SRW was delivered to the folder:\n";
	strcpy(StrTot, Str1);
	strcat(StrTot, WhereHePutIt);
	char Str2[] = "\n\nTo complete the installation, please\nread the instructions in the file:\n";
	strcat(StrTot, Str2);
	strcat(StrTot, PathReadMe);

	UINT DlgStyle = MB_OK | MB_ICONASTERISK | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
	int MesBoxInf = MessageBox(NULL, StrTot, MesTitle, DlgStyle); 
}

/////////////////////////////////////////////////////////////////////////////
/**
void CSetupApp::SRWDirectoryStructure(ostringstream& o)
{
	o << "dir:SRW" << ends;
	o << "CopyRight.txt" << ends;
	o << "ReadMe.txt" << ends;
	o << "SRW.xop" << ends;

	o << "dir:SRW Help" << ends;
	o << "SRW Help.ifn" << ends;
	o << "SRW Macro Help.ifn" << ends;
	o << "SRW Macro Help.ifn.igr" << ends;
	o << "enddir" << ends;

	o << "dir:SRW Procedures" << ends;
	o << "dir:Development" << ends;
	o << "SRW Misc.ipf" << ends;
	o << "SRW Genesis.ipf" << ends;
	o << "SRW Radiation Stokes Arb.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Examples" << ends;
	o << "SRW Example BM SR Focusing.ipf" << ends;
	o << "SRW Example BM Std.ipf" << ends;
	o << "SRW Example Brilliance.ipf" << ends;
	o << "SRW Example ER Diffraction.ipf" << ends;
	o << "SRW Example ER Power.ipf" << ends;
	o << "SRW Example ER.ipf" << ends;
	o << "SRW Example GsnBm.ipf" << ends;
	o << "SRW Example UR AngPatThkE.ipf" << ends;
	o << "SRW Example UR Ell.ipf" << ends;
	o << "SRW Example UR Focusing.ipf" << ends;
	o << "SRW Example UR Off-Axis.ipf" << ends;
	o << "SRW Example UR Power.ipf" << ends;
	o << "SRW Example UR SpecThkE.ipf" << ends;
	o << "SRW Example WigEll.ipf" << ends;
	o << "SRW Example WigPlan.ipf" << ends;
	o << "SRW Example X-rayLensCirc.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Globals" << ends;
	o << "SRW Globals.ipf" << ends;
	o << "SRW Menu.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Input" << ends;
	o << "SRW Electron.ipf" << ends;
	o << "SRW MagField Const.ipf" << ends;
	o << "SRW MagField Gen.ipf" << ends;
	o << "SRW MagField Periodic.ipf" << ends;
	o << "SRW Observation.ipf" << ends;
	o << "SRW Other Sources.ipf" << ends;
	o << "SRW Visualize.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Optics" << ends;
	o << "SRW Optics ThinGen.ipf" << ends;
	o << "SRW Optics.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Radiation" << ends;
	o << "SRW Radiation Brilliance.ipf" << ends;
	o << "SRW Radiation Manip.ipf" << ends;
	o << "SRW Radiation NearField.ipf" << ends;
	o << "SRW Radiation Power.ipf" << ends;
	o << "SRW Radiation Und.ipf" << ends;
	o << "SRW Radiation Wiggler.ipf" << ends;
	o << "enddir" << ends;

	o << "dir:Utilities" << ends;
	o << "SRW Utilities.ipf" << ends;
	o << "SRW View Struct.ipf" << ends;
	o << "enddir" << ends;

	o << "enddir" << ends;
	o << "enddir" << ends;
	o << ends;
}
**/
/////////////////////////////////////////////////////////////////////////////
