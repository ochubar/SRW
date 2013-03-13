
#include <afxdlgs.h>

class CDirectoryDialog : public CFileDialog {

public:

	CDirectoryDialog(BOOL bOpenFileDialog, // TRUE for FileOpen, FALSE for FileSaveAs
		LPCTSTR lpszDefExt = NULL,
		LPCTSTR lpszFileName = NULL,
		DWORD dwFlags = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		LPCTSTR lpszFilter = NULL,
		CWnd* pParentWnd = NULL)
		: CFileDialog(bOpenFileDialog, lpszDefExt, lpszFileName, dwFlags, lpszFilter, pParentWnd)
	{
		//CFileDialog


	}

	void OnFolderChange()
	{

		CString FolderPath = GetFolderPath();
		int Aga = 1;

	}
	BOOL OnFileNameOK()
	{
		int Aga = 1;

		return 0;
	}

};