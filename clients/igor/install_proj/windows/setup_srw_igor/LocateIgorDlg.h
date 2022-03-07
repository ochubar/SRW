#if !defined(AFX_LOCATEIGORDLG_H__9B65977B_E4FB_4270_8E40_6D5A3CB5945A__INCLUDED_)
#define AFX_LOCATEIGORDLG_H__9B65977B_E4FB_4270_8E40_6D5A3CB5945A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// LocateIgorDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CLocateIgorDlg dialog

class CLocateIgorDlg : public CDialog
{
// Construction
public:
	//CLocateIgorDlg(CWnd* pParent = NULL);   // standard constructor
	CLocateIgorDlg(char* pPathIgorExe, CWnd* pParent =NULL);

// Dialog Data
	//{{AFX_DATA(CLocateIgorDlg)
	enum { IDD = IDD_LOCATEIGORDLG_DIALOG };
	CString	m_IgorPath;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CLocateIgorDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CLocateIgorDlg)
	afx_msg void OnBrowseButton();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_LOCATEIGORDLG_H__9B65977B_E4FB_4270_8E40_6D5A3CB5945A__INCLUDED_)
