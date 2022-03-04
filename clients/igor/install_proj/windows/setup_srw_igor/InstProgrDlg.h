#if !defined(AFX_INSTPROGRDLG_H__054856BB_8D9C_4A24_88EC_681F865C3B60__INCLUDED_)
#define AFX_INSTPROGRDLG_H__054856BB_8D9C_4A24_88EC_681F865C3B60__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// InstProgrDlg.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CInstProgrDlg dialog

class CInstProgrDlg : public CDialog
{
	BOOL WasInstalled;
// Construction
public:
	CInstProgrDlg(int AmOfSteps, CWnd* pParent = NULL);   // standard constructor

	void SetPosition(int CurStep, char* CurDescr);

// Dialog Data
	//{{AFX_DATA(CInstProgrDlg)
	enum { IDD = IDD_INSTPROGRDLG_DIALOG };
	CProgressCtrl	m_ProgrInd;
	CString	m_Descr;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CInstProgrDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CInstProgrDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_INSTPROGRDLG_H__054856BB_8D9C_4A24_88EC_681F865C3B60__INCLUDED_)
