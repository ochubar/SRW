// This file contains equates and prototypes that are needed on Windows only.

#ifdef __cplusplus
extern "C" {						/* This allows C++ to call the XOPSupport routines */
#endif

/* Windows-specific routines (in XOPWinSupport.c) */
HMODULE XOPModule(void);
HMODULE IgorModule(void);
HWND IgorClientHWND(void);
void debugstr(const char *text);
int SendWinMessageToIgor(HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam, int beforeOrAfter);
void PositionWinDialogWindow(HWND theDialog, HWND refWindow);
int IsWinDialogItemHitMessage(XOP_DIALOG_REF theDialog, int itemID, int notificationMessage);
void HideDialogItem(XOP_DIALOG_REF theDialog, int itemID);
void ShowDialogItem(XOP_DIALOG_REF theDialog, int itemID);

#ifdef __cplusplus
}
#endif
