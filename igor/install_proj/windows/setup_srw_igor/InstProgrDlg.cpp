// InstProgrDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Setup.h"
#include "InstProgrDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CInstProgrDlg dialog

CInstProgrDlg::CInstProgrDlg(int AmOfSteps, CWnd* pParent /*=NULL*/)
	: CDialog(CInstProgrDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CInstProgrDlg)
	m_Descr = _T("");
	//}}AFX_DATA_INIT

	WasInstalled = Create(IDD_INSTPROGRDLG_DIALOG);
	if(WasInstalled) 
	{
		ShowWindow(SW_SHOW);
		m_ProgrInd.SetRange(0, AmOfSteps);
	}
}

void CInstProgrDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CInstProgrDlg)
	DDX_Control(pDX, IDC_PROGRESS1, m_ProgrInd);
	DDX_Text(pDX, IDC_EDIT2, m_Descr);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CInstProgrDlg, CDialog)
	//{{AFX_MSG_MAP(CInstProgrDlg)
	ON_BN_CLICKED(IDC_CANCEL, OnCancel)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CInstProgrDlg message handlers

/////////////////////////////////////////////////////////////////////////////

void CInstProgrDlg::SetPosition(int CurStep, char* CurDescr)
{
	if(!WasInstalled) return;

	const int MaxStrLen = 1024;
	char TextToSet[MaxStrLen];
	if((CurDescr == 0) || (*CurDescr == '\0')) strcpy(TextToSet, "");
	else if(strlen(CurDescr) > MaxStrLen) strncpy(TextToSet, CurDescr, MaxStrLen);
	else strcpy(TextToSet, CurDescr);

	m_ProgrInd.SetPos(CurStep);
	m_Descr = _T(TextToSet);
	UpdateData(FALSE);
	UpdateWindow();
}

/////////////////////////////////////////////////////////////////////////////
/*
void CInstProgrDlg::OnCancel() 
{
    AfxMessageBox("SRW installation cancelled by User.");
    AfxThrowUserException();
}
*/
