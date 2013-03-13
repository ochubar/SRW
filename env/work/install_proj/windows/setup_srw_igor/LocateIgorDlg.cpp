// LocateIgorDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Setup.h"
#include "LocateIgorDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CLocateIgorDlg dialog

//CLocateIgorDlg::CLocateIgorDlg(CWnd* pParent /*=NULL*/)
//	: CDialog(CLocateIgorDlg::IDD, pParent)
//{
//	//{{AFX_DATA_INIT(CLocateIgorDlg)
//	m_IgorPath = _T("");
//	//}}AFX_DATA_INIT
//}

CLocateIgorDlg::CLocateIgorDlg(char* pPathIgorExe, CWnd* pParent /*=NULL*/)
	: CDialog(CLocateIgorDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CLocateIgorDlg)
	m_IgorPath = _T(pPathIgorExe);
	//}}AFX_DATA_INIT
}

void CLocateIgorDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CLocateIgorDlg)
	DDX_Text(pDX, IDC_EDIT1, m_IgorPath);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CLocateIgorDlg, CDialog)
	//{{AFX_MSG_MAP(CLocateIgorDlg)
	ON_BN_CLICKED(IDC_BUTTON1, OnBrowseButton)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CLocateIgorDlg message handlers

void CLocateIgorDlg::OnBrowseButton() 
{
	CFileDialog aFileDlg(TRUE, "", m_IgorPath);
	int nResponse = aFileDlg.DoModal();

	if(nResponse == IDOK)
	{
		m_IgorPath = _T(aFileDlg.m_ofn.lpstrFile);
		UpdateData(FALSE);
	}
}
