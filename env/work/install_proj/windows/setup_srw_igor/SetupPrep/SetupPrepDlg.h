// SetupPrepDlg.h : header file
//

#if !defined(AFX_SETUPPREPDLG_H__F74F9E67_832A_4AEC_9DC5_4DFA8E2A4CCE__INCLUDED_)
#define AFX_SETUPPREPDLG_H__F74F9E67_832A_4AEC_9DC5_4DFA8E2A4CCE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CSetupPrepDlg dialog

class CSetupPrepDlg : public CDialog
{
// Construction
public:
	CSetupPrepDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CSetupPrepDlg)
	enum { IDD = IDD_SETUPPREP_DIALOG };
	CString	m_Path;
	CString	m_Explanation;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSetupPrepDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CSetupPrepDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnBrowseButton();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SETUPPREPDLG_H__F74F9E67_832A_4AEC_9DC5_4DFA8E2A4CCE__INCLUDED_)
