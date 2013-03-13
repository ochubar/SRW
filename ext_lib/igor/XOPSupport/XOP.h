/*
 *	XOP.h
 *		external operations equates, data structures and globals
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

/*	Development System Codes - Use DEV_SYS_CODE in the second field of the
	XOPI resource and stored in the system field of the XOPStuff structure.
	These codes are also defined in XOPResources.h for use in .r and .rc files.
*/
#define DEV_SYS_CODE_ALL_OTHERS 0		// Use this for any type of development not listed below.
#define DEV_SYS_CODE_MAC_MPW 1			// Obsolete.
#define DEV_SYS_CODE_MAC_MACH 2			// Use this for Macintosh Mach-O XOPs.

/*	XOP version history
	
	XOP_VERSION identifies the XOP Toolkit protocol. It is changed only when
	a major change in protocol occurs (that is, almost never).
	
	XOP_VERSION = 1: Released with Igor 1.1
	XOP_VERSION = 3: Released with Igor 1.2
	XOP_VERSION = 4: Released with Igor 2.0 (Igor Pro 2.0D82)	10/1/93
	XOP_VERSION = 5: Released with Igor Pro 5.0, XOP Toolkit 5 beta May, 2003
	XOP_VERSION = 4: Released with Igor Pro 5.01, XOP Toolkit 5 beta December, 2003
	
	HR, 031212: I realized that changing XOP_VERSION to 5 meant that XOPs compiled with
	XOP Toolkit 5 would be considered incompatible by Igor Pro 4. This is not necessary.
	Therefore, with Igor Pro 5.01 and XOP Toolkit 5 Beta from December, 2003, I rolled
	XOP_VERSION back to 4. However, I left XOP_HIGH_VERSION at 5 so that XOPs compiled
	with early beta versions of the toolkit would continue to work with Igor Pro 5.
*/
#define XOP_VERSION 4						/* current XOP version; NOTE: this is in XOPResources.h also. */
#define XOP_LOW_VERSION 2					/* lowest supported XOP version */
#define XOP_HIGH_VERSION 5					/* highest supported XOP version */

#define MAX_XOP_NAME MAX_OBJ_NAME			// This relationship is assumed.
#define MAX_XOP_FILE_NAME (MAX_XOP_NAME+4)	// +4 to allow for extension on Windows.

/*	Structures used for XOPs  */

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

struct XOPStuff {					/* structure of private bookkeeping info */
	struct XOPStuff **nextXOPHandle;/* for linked list of open XOPs */
	char XOPName[MAX_XOP_NAME+1];	/* name of XOP */
	struct IORec **ioRecHandle;		/* handle to main ioRec for this XOP */
	unsigned char flags;			/* private flags used by Igor; 10/28/93 */
	unsigned char system;			/* code for development system */

#ifdef _MACINTOSH_					// Mac-specific fields [
	Handle mapHandle;				/* handle to XOPs resource map */
	unsigned long oldA4;			/* Igor's A4 */
	unsigned long oldA5;			/* Igor's A5 */
	short oldApRefNum;				/* Igor's file reference number */
	unsigned long newA4;			/* XOP's A4 */
	unsigned long newA5;			/* XOP's A5 */
	short newApRefNum;				/* XOP's file reference number */
	short curResFile;				/* XOP's curResFile */
	Ptr worldPtr;					/* pointer to XOP's A5 world */
	
	/* HR: 7/29/94: These fields were added for PowerMac and code fragment manager support. */
	int ISA;						/* instruction set architecture: kM68kISA or kPowerPCISA */
	int isCodeFragmentXOP;			/* true if we are this is a code fragment XOP */
	long connID;					/* connectionID for code fragment or zero */
#endif								// End Mac-specific fields ]

#ifdef _WINDOWS_					// Windows-specific fields [
	HMODULE hModule;				// XOP DLL module handle.
#endif								// End Windows-specific fields ]
};
typedef struct XOPStuff XOPStuff;
typedef XOPStuff *XOPStuffPtr;
typedef XOPStuffPtr *XOPStuffHandle;

/*	On Macintosh, we use a file reference number to identify an XOP (the newApRefNum
	field in XOPStuff. On Windows, we use a handle to a DLL module (the hModule field).
	We define the symbol XOP_MODULE_REF such that it can be used on either platform for
	the type of the thing that identifies a particular XOP.
*/	
#ifdef _MACINTOSH_
	#define XOP_MODULE_REF int
	#define NULL_XOP_MODULE 0
#endif
#ifdef _WINDOWS_
	#define XOP_MODULE_REF HMODULE
	#define NULL_XOP_MODULE NULL
#endif

/*	On Macintosh an XOP_WINDOW_REF is a WindowPtr. On Windows, it is an HWND.
	Igor passes a XOP_WINDOW_REF to the XOP when it sends the XOP window-related
	messages (e.g. ACTIVATE, UPDATE, COPY, PRINT).
*/
#ifdef _MACINTOSH_
	#define XOP_WINDOW_REF WindowPtr
#endif
#ifdef _WINDOWS_
	#define XOP_WINDOW_REF HWND
#endif

// On Macintosh an XOP_DIALOG_REF is a DialogPtr. On Windows, it is an HWND.
#ifdef _MACINTOSH_
	#define XOP_DIALOG_REF DialogPtr
#endif
#ifdef _WINDOWS_
	#define XOP_DIALOG_REF HWND
#endif

/*	We are now using standard C file I/O routines for file I/O.
	Using this #define instead of FILE* makes is easy to switch
	to native routines if desired in the future.
*/
#define XOP_FILE_REF FILE*

// For flags field in XOPStuff.
#define XOP_RUNNING 1				// Used by Igor to detect recursion.
#define XOP_MUST_STAY_RESIDENT 2	// Used by Igor to prevent unloading an XOP that adds a direct XFUNC or a direct operation.

/*	Bit-mask codes for system field in XOPStuff structure.
	These are used internally by Igor and are of no concern to the XOP programmer.
*/
#define MULTI_SEG 1					// XOP is MPW style multi-segment with standard Macintosh CODE 0.
#define NO_DOUBLE 2					// OBSOLETE: XOP uses extended instead of double. This was used in the days of THINK C 3.0.

struct IORec {
	short version;							/* XOP protocol version */
	long XOPType;							/* transient/resident, idle/no idle and other info */
	long result;							/* result code from XOP or from callback -- inited to 0 */
	long status;							/* various status info depending on operation */
#ifdef _MACINTOSH_
	pascal void (*callBackProc)(void *);	/* address of routine for calling Igor back */
#endif
#ifdef _WINDOWS_
	void (*callBackProc)(void *);			/* address of routine for calling Igor back */
#endif
	void (*XOPEntry)(void);					/* address for calling XOP from host */
	XOPStuffHandle stuffHandle;				/* used by Igor for bookkeeping purposes */
	long refCon;							/* for use by XOP --initially zero */
	short message;							/* bidirectional message passed between host and XOP */
	short menuID;							/* ID of menu for XOP or 0 if none */
	short itemID;							/* number of menu item for XOP of 0 if none */
	short subMenuID;						/* ID of submenu for XOP or 0 if none */
	short numItems;							/* total number of items in list */
	long items[1];							/* list of items to operate on -- variable length array */
};
typedef struct IORec IORec;
typedef IORec *IORecPtr;
typedef IORecPtr *IORecHandle;

/*	HR, 3/2/95, for Igor Pro 3.0:

	Previously, both Igor and the XOP would expand or shrink the
	IORec handle depending on the number of items that needed to
	fit in it. 
	
	I have changed Igor Pro 2.1 to allocate a number of IORec items
	which will be sufficient for all time. This is to avoid the
	need to do a SetHandleSize which can fail with some (perhaps all)
	versions of the Mac memory manager if the handle is trapped between
	locked blocks.
	
	I have made a similar change in XOPSupport.c. However, old XOPs
	will still make the handle smaller. Igor Pro will cope with this.
*/
#define NUM_IOREC_ITEMS 16			/* Used in SetRecHandle() */


// XOP type codes identify capabilities of XOP -- used to set XOPType field.
#define TRANSIENT 0					// Default for XOPType.
#define RESIDENT 1					// XOP wants to stick around indefinitely.
#define IDLES 2						// XOP has idle routine to be called periodically.
#define ADDS_TARGET_WINDOWS 4		// XOP has the capability of adding one or more target windows to Igor. For Igor Pro 3.13B03.


/* XOP status codes used by Igor to inform XOP */
#define INFOREGROUND 0				/* Igor is in foreground */
#define INBACKGROUND 1				/* Igor is in background */
#define XOPINITING 2				/* XOP is being inited */


/* XOP error #defines */
#define XOP_ERRS_ID 1100			/* resource ID for STR# resource for custom XOP errors */
#define FIRST_XOP_ERR 10000
#define LAST_XOP_ERR 11999


/* Flags used in XOP menus and XOP menu item resources. */
/* A copy of these is in XOPTypes.r so they can be used in .r files */
#define SHOW_MENU_AT_LAUNCH 1				/* for XMN1 resource menuFlags field */
#define SHOW_MENU_WHEN_ACTIVE 2				/* for XMN1 resource menuFlags field */
#define ITEM_REQUIRES_WAVES 1				/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_GRAPH 2				/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_TABLE 4				/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_LAYOUT 8				/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_GRAPH_OR_TABLE 16		/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_TARGET 32				/* for XMI1 resource itemFlags field */
#define ITEM_REQUIRES_PANEL 64				/* for XMI1 resource itemFlags field; added for Igor version 2.0 */
#define ITEM_REQUIRES_NOTEBOOK 128			/* for XMI1 resource itemFlags field; added for Igor version 2.0 */
#define ITEM_REQUIRES_GRAPH_OR_PANEL 256	/* for XMI1 resource itemFlags field; added for Igor version 2.0 */
#define ITEM_REQUIRES_DRAW_WIN 512			/* for XMI1 resource itemFlags field; added for Igor version 2.0 */
#define ITEM_REQUIRES_PROC_WIN 1024			/* for XMI1 resource itemFlags field; added for Igor version 2.0 */


/* XOP Message Codes -- codes passed from host to XOP */
#define INIT 1
#define IDLE 2
#define CMD 3
#define NEW 4
#define LOAD 5							// Obsolete - use LOADSETTINGS instead.
#define SAVE 6							// Obsolete - use SAVESETTINGS instead.
#define SAVESETTINGS 7
#define LOADSETTINGS 8
#define MODIFIED 9
#define MENUITEM 10
#define MENUENABLE 11
#define CLEANUP 12
#define OBJINUSE 13
#define FUNCTION 14						// Added in Igor Pro 2.0D54.		
#define FUNCADDRS 15					// Added in Igor Pro 2.0D54.
#define SAVE_WINDOW_MACRO 16			// Added in Igor Pro 3.13B03. For target windows only.
#define GET_TARGET_WINDOW_NAME 17		// Added in Igor Pro 3.13B03. For target windows only.
#define SET_TARGET_WINDOW_NAME 18		// Added in Igor Pro 3.13B03. For target windows only.
#define GET_TARGET_WINDOW_REF 19		// Added in Igor Pro 3.13B03. For target windows only.
#define SET_TARGET_WINDOW_TITLE 20		// Added in Igor Pro 3.13B03. For target windows only.
#define CLEAR_MODIFIED 21				// Added in Igor Pro 3.13B03.
#define EXECUTE_OPERATION 22			// Added in Igor Pro 5.00. For testing only. Do not use.


/*	Windowed XOP routine codes -- event codes passed from host to XOP */
#define XOPWINDOWCODE 512
#define EVENT 1 + XOPWINDOWCODE					/* This is obsolete */
#define ACTIVATE 2 + XOPWINDOWCODE
#define UPDATE 3 + XOPWINDOWCODE
#define GROW 4 + XOPWINDOWCODE
#define SETGROW 5 + XOPWINDOWCODE
#define CLICK 6 + XOPWINDOWCODE
#define KEY 7 + XOPWINDOWCODE
#define NULLEVENT 8 + XOPWINDOWCODE
#define UNDO 9 + XOPWINDOWCODE
#define COPY 10 + XOPWINDOWCODE
#define CUT 11 + XOPWINDOWCODE
#define PASTE 12 + XOPWINDOWCODE
#define CLEAR 13 + XOPWINDOWCODE
#define FIND 14 + XOPWINDOWCODE
#define DISPLAYSELECTION 15 + XOPWINDOWCODE
#define REPLACE 16 + XOPWINDOWCODE
#define INDENTLEFT 17 + XOPWINDOWCODE
#define INDENTRIGHT 18 + XOPWINDOWCODE
#define PAGESETUP 19 + XOPWINDOWCODE
#define PRINT 20 + XOPWINDOWCODE
#define CLOSE 21 + XOPWINDOWCODE
#define INSERTFILE 22 + XOPWINDOWCODE
#define SAVEFILE 23 + XOPWINDOWCODE
#define DUPLICATE 24 + XOPWINDOWCODE			/* added in Igor Pro 2.0D82 */
#define SELECT_ALL 25 + XOPWINDOWCODE			/* added in Igor Pro 2.0D82 */
#define EXPORT_GRAPHICS 26 + XOPWINDOWCODE		/* added in Igor Pro 2.0D82 */
#define SAVE_WIN 27 + XOPWINDOWCODE				/* added in Igor Pro 2.0D82 */
#define SAVE_WIN_AS 28 + XOPWINDOWCODE			/* added in Igor Pro 2.0D82 */
#define SAVE_WIN_COPY 29 + XOPWINDOWCODE		/* added in Igor Pro 2.0D82 */
#define REVERT_WIN 30 + XOPWINDOWCODE			/* added in Igor Pro 2.0D82 */
#define MOVE_TO_PREFERRED_POSITION 31 + XOPWINDOWCODE	/* added in Igor Pro 3.10B01 */
#define MOVE_TO_FULL_POSITION 32 + XOPWINDOWCODE		/* added in Igor Pro 3.10B01 */
#define RETRIEVE 33 + XOPWINDOWCODE						/* added in Igor Pro 3.10B01 */
#define WINDOW_MOVED 34 + XOPWINDOWCODE					/* added in Igor Pro 3.10B01 */
#define XOP_HIDE_WINDOW 35 + XOPWINDOWCODE				// Added in Igor Pro 6.00D00
#define XOP_SHOW_WINDOW 36 + XOPWINDOWCODE				// Added in Igor Pro 6.00D00


/* Callback operation codes -- codes passed from XOP to host */
#define WAVEMODIFIED 1
#define WAVEHANDLEMODIFIED 2
#define FETCHNUMVAR 3
#define FETCHSTRVAR 4
#define STORENUMVAR 5
#define STORESTRVAR 6
#define FETCHWAVE 7
#define NOTICE 8
#define COMMAND 9
#define NEXTSYMB 10
#define GETSYMB 11
#define GETFLAG 12
#define GETFLAGNUM 13
#define GETNAME 14
#define GETVAR 15
#define GETNUM 16
#define GETLONG 17
#define GETSTRING 18
#define GETWAVE 19
#define GETWAVELIST 20
#define CHECKTERM 21
#define UNIQUENAME 22
#define GETPATH 23
#define GETFORMAT 24
#define GETWRITESTRING 25		/* OBSOLETE -- use CHIO instead */
#define GETWRITEWAVES 26		/* OBSOLETE -- use CHIO instead */
#define WRITEWAVES 27			/* OBSOLETE -- use CHIO instead */
#define GETREAD 28				/* OBSOLETE -- use CHIO instead */
#define PAUSEUPDATE 29
#define RESUMEUPDATE 30
#define SILENT 31				/* OBSOLETE -- does nothing */
#define UNSILENT 32				/* OBSOLETE -- does nothing */
#define SETCURSOR 33
#define WAVETYPE 34
#define WAVEPOINTS 35
#define WAVENAME 36
#define WAVEDATA 37
#define DEFAULTMENUS 38
#define CHIO 39
#define GETWAVERANGE 40
#define CALCWAVERANGE 41
#define MAKEWAVE 42
#define CHANGEWAVE 43
#define KILLWAVE 44
#define VARIABLE 45
#define IGORERROR 46
#define WAVESMODIFIED 47
#define SPINPROCESS 48
#define PUTCMDLINE 49
#define WAVESCALING 50
#define SETWAVESCALING 51		/* added in Igor 1.20 */
#define WAVEUNITS 52
#define SETWAVEUNITS 53
#define WAVENOTE 54
#define SETWAVENOTE 55
#define WAVELIST 56				/* added in Igor 1.24 */
#define WINLIST 57
#define PATHLIST 58
#define FETCHSTRHANDLE 59
#define GETNUM2 60
#define ACTUALMENUID 61
#define RESOURCEMENUID 62
#define ACTUALITEMNUM 63
#define RESOURCEITEMNUM 64
#define WININFO 65
#define PATHINFO 66				/* added in Igor 1.25 */
#define GETNAMEDFIFO 67			/* added in Igor Pro 2.0D54 */
#define MARKFIFOUPDATED 68
#define SETIGORMENUITEM 69
#define SILENT_COMMAND 70
#define DOUPDATE 71
#define XOPMENUHANDLE 72
#define GETSTRINGINHANDLE 73
#define VARIABLELIST 74
#define STRINGLIST 75
#define UNIQUENAME2 76
#define MD_MAKEWAVE 77			/* added in Igor Pro 3.0 */
#define MD_GETWAVEDIMENSIONS 78
#define MD_CHANGEWAVE 79
#define MD_GETWAVESCALING 80
#define MD_SETWAVESCALING 81
#define MD_GETWAVEUNITS 82
#define MD_SETWAVEUNITS 83
#define MD_GETDIMLABELS 84
#define MD_SETDIMLABELS 85
#define MD_ACCESSNUMERICWAVEDATA 86
#define MD_GETWAVEPOINTVALUE 87
#define MD_SETWAVEPOINTVALUE 88
#define MD_GETDPDATAFROMNUMERICWAVE 89
#define MD_STOREDPDATAINNUMERICWAVE 90
#define MD_GETTEXTWAVEPOINTVALUE 91
#define MD_SETTEXTWAVEPOINTVALUE 92
#define WAVEMODDATE 93						// Added in Igor Pro 3.01
#define WAVEMODSTATE 94						// Added in Igor Pro 4.0D08 but works with all versions of Igor.
#define WAVEMODCOUNT 95						// Added in Igor Pro 4.0D08
#define MD_CHANGEWAVE2 96					// Added in Igor Pro 5.04B06
/* 97 - 99 are reserved. */

#define GET_DATAFOLDER_NAMEORPATH 100		/* added in Igor Pro 3.0 */
#define GET_DATAFOLDER_IDNUMBER 101
#define GET_DATAFOLDER_PROPERTIES 102
#define SET_DATAFOLDER_PROPERTIES 103
#define GET_DATAFOLDER_LISTING 104
#define GETROOT_DATAFOLDER 105
#define GETCURRENT_DATAFOLDER 106
#define SETCURRENT_DATAFOLDER 107
#define GETNAMED_DATAFOLDER 108
#define GET_DATAFOLDER_BYIDNUMBER 109
#define GETPARENT_DATAFOLDER 110
#define GETNUMCHILD_DATAFOLDERS 111
#define GETINDEXEDCHILD_DATAFOLDER 112
#define GETWAVES_DATAFOLDER 113
#define NEW_DATAFOLDER 114
#define KILL_DATAFOLDER 115
#define DUPLICATE_DATAFOLDER 116
#define MOVE_DATAFOLDER 117
#define RENAME_DATAFOLDER 118
#define GETNUM_DATAFOLDER_OBJECTS 119
#define GETINDEXED_DATAFOLDER_OBJECT 120
#define KILL_DATAFOLDER_OBJECT 121
#define MOVE_DATAFOLDER_OBJECT 122
#define RENAME_DATAFOLDER_OBJECT 123
#define GET_DATAFOLDER_AND_NAME 124
#define CLEAR_DATAFOLDER_FLAGS 125
#define GET_DATAFOLDER_CHANGESCOUNT 126
#define GET_DATAFOLDER_CHANGEFLAGS 127
#define DUPLICATE_DATAFOLDER_OBJECT 128

#define CHECKNAME 150						/* Added in Igor Pro 3.0. */
#define POSSIBLY_QUOTE_NAME 151				/* Added in Igor Pro 3.0. */
#define CLEANUP_NAME 152					/* Added in Igor Pro 3.0. */
#define PREPARE_LOAD_IGOR_DATA 153			/* For internal WaveMetrics use only (used by WaveMetrics Browser). */
#define DO_LOAD_IGOR_DATA 154				/* For internal WaveMetrics use only (used by WaveMetrics Browser). */
#define END_LOAD_IGOR_DATA 155				/* For internal WaveMetrics use only (used by WaveMetrics Browser). */
#define IS_STRING_EXPRESSION 156			/* Added in Igor Pro 3.0. */
#define FETCHWAVE_FROM_DATAFOLDER 157		/* Added in Igor Pro 3.0. */
#define GET_DATAFOLDER 158					/* Added in Igor Pro 3.0. */
#define GET_DATAFOLDER_OBJECT 159			/* Added in Igor Pro 3.0. */
#define SET_DATAFOLDER_OBJECT 160			/* Added in Igor Pro 3.0. */

#define SAVE_XOP_PREFS 161					// Added in IGOR Pro 3.10.
#define GET_XOP_PREFS 162					// Added in IGOR Pro 3.10.
#define GET_PREFS_STATE 163					// Added in IGOR Pro 3.10.

#define GET_INDEXED_IGORCOLORTABLE_NAME 164	// Added in IGOR Pro 3.10.
#define GET_NAMED_IGORCOLORTABLE_HANDLE 165	// Added in IGOR Pro 3.10.
#define GET_IGORCOLORTABLE_INFO 166			// Added in IGOR Pro 3.10.
#define GET_IGORCOLORTABLE_VALUES 167		// Added in IGOR Pro 3.10.

#define PATHINFO2 168						// Added in Igor Pro 3.13.
#define GET_DIR_AND_FILE_FROM_FULL_PATH 169	// Added in Igor Pro 3.13.
#define CONCATENATE_PATHS 170				// Added in Igor Pro 3.13.
#define WIN_TO_MAC_PATH 171					// Added in Igor Pro 3.13.
#define MAC_TO_WIN_PATH 172					// Added in Igor Pro 3.13.
#define STRCHR2 173							// Added in Igor Pro 3.13.
#define STRRCHR2 174						// Added in Igor Pro 3.13.

#define DISPLAY_HELP_TOPIC 175				// Added in Igor Pro 3.13.

#define WINDOW_RECREATION_DIALOG 176		// Added in Igor Pro 3.13.
#define GET_IGOR_PROCEDURE_LIST 177			// Added in Igor Pro 3.13.
#define GET_IGOR_PROCEDURE 178				// Added in Igor Pro 3.13.
#define SET_IGOR_PROCEDURE 179				// Added in Igor Pro 3.13.
#define GET_IGOR_ERROR_MESSAGE 180			// Added in Igor Pro 3.14.

#define COMMAND2 181						// Added in Igor Pro 4.0.
#define AT_END_OF_COMMAND 182				// Added in Igor Pro 4.0.

#define SET_CONTEXTUAL_HELP_DIALOG_ID 183	// Added in Igor Pro 5.0. This is a NOP on Windows.
#define SHOW_HIDE_CONTEXTUAL_HELP 184		// Added in Igor Pro 5.0. This is a NOP on Windows.
#define SET_CONTEXTUAL_HELP_MESSAGE 185		// Added in Igor Pro 5.0. Works on Macintosh and Windows.

#define DO_SAVE_IGOR_DATA 186				// For internal WaveMetrics use only (used by WaveMetrics Browser).

#define REGISTER_OPERATION 187					// Added in Igor Pro 5.00.
#define SET_RUNTIME_NUMERIC_VARIABLE 188		// Added in Igor Pro 5.0.
#define SET_RUNTIME_STRING_VARIABLE 189			// Added in Igor Pro 5.0.
#define VAR_NAME_TO_DATA_TYPE 190				// Added in Igor Pro 5.0.
#define STORE_NUMERIC_DATA_USING_VAR_NAME 191	// Added in Igor Pro 5.0.
#define STORE_STRING_DATA_USING_VAR_NAME 192	// Added in Igor Pro 5.0.

#define GET_FUNCTION_INFO 193					// Added in Igor Pro 5.0.
#define CHECK_FUNCTION_FORM 194					// Added in Igor Pro 5.0.
#define CALL_FUNCTION 195						// Added in Igor Pro 5.0.

#define WAVELOCK 196							// Added in Igor Pro 5.0.
#define SETWAVELOCK 197							// Added in Igor Pro 5.0.
#define DATE_TO_IGOR_DATE 198					// Added in Igor Pro 5.0.
#define IGOR_DATE_TO_DATE 199					// Added in Igor Pro 5.0.

#define GET_NVAR 200							// Added in Igor Pro 5.03.
#define SET_NVAR 201							// Added in Igor Pro 5.03.
#define GET_SVAR 202							// Added in Igor Pro 5.03.
#define SET_SVAR 203							// Added in Igor Pro 5.03.
#define GET_FUNCTION_INFO_FROM_FUNCREF 204		// Added in Igor Pro 5.03.

#define GET_TEXT_WAVE_DATA 205					// Added in Igor Pro 5.04.
#define SET_TEXT_WAVE_DATA 206					// Added in Igor Pro 5.04.
#define GET_WAVE_DIMENSION_LABELS 207			// Added in Igor Pro 5.04.
#define SET_WAVE_DIMENSION_LABELS 208			// Added in Igor Pro 5.04.

#define SET_OPERATION_WAVE_REF 209				// Added in Igor Pro 5.04.

#define TELL_IGOR_WINDOW_STATUS 210				// Added in Igor Pro 6.00D00

/*	Text utility callback operation codes -- callback codes passed from XOP to host  */
#define XOPTUCODE 512
#define TUNEW 1 + XOPTUCODE
#define TUDISPOSE 2 + XOPTUCODE
#define TUDISPLAYSELECTION 3 + XOPTUCODE
#define TUGROW 4 + XOPTUCODE
#define TUZOOM 5 + XOPTUCODE
#define TUDRAWWINDOW 6 + XOPTUCODE
#define TUUPDATE 7 + XOPTUCODE
#define TUFIND 8 + XOPTUCODE
#define TUREPLACE 9 + XOPTUCODE
#define TUINDENTLEFT 10 + XOPTUCODE
#define TUINDENTRIGHT 11 + XOPTUCODE
#define TUCLICK 12 + XOPTUCODE
#define TUACTIVATE 13 + XOPTUCODE
#define TUIDLE 14 + XOPTUCODE
#define TUNULL 15 + XOPTUCODE
#define TUCOPY 16 + XOPTUCODE
#define TUCUT 17 + XOPTUCODE
#define TUPASTE 18 + XOPTUCODE
#define TUCLEAR 19 + XOPTUCODE
#define TUKEY 20 + XOPTUCODE
#define TUINSERT 21 + XOPTUCODE
#define TUDELETE 22 + XOPTUCODE
#define TUSETSELECT 23 + XOPTUCODE
#define TUFIXEDITMENU 24 + XOPTUCODE
#define TUFIXFILEMENU 25 + XOPTUCODE
#define TUUNDO 26 + XOPTUCODE
#define TUPRINT 27 + XOPTUCODE
#define TULENGTH 28 + XOPTUCODE
#define TULINES 29 + XOPTUCODE
#define TUSELSTART 30 + XOPTUCODE
#define TUSELEND 31 + XOPTUCODE
#define TUSELLENGTH 32 + XOPTUCODE
#define TUGETTEXT 33 + XOPTUCODE
#define TUFETCHTEXT 34 + XOPTUCODE
#define TUINSERTFILE 35 + XOPTUCODE
#define TUWRITEFILE 36 + XOPTUCODE
#define TUSFINSERTFILE 37 + XOPTUCODE
#define TUSFWRITEFILE 38 + XOPTUCODE
#define TUPAGESETUPDIALOG 39 + XOPTUCODE	/* added in Igor Pro 2.0D82 */
#define TUSELECTALL 40 + XOPTUCODE			/* added in Igor Pro 2.0D82 */
#define TUGETDOCINFO 41 + XOPTUCODE			/* added in Igor Pro 2.0D83 */
#define TUGETSELLOCS 42 + XOPTUCODE			/* added in Igor Pro 2.0D83 */
#define TUSETSELLOCS 43 + XOPTUCODE			/* added in Igor Pro 2.0D83 */
#define TUFETCHPARAGRAPHTEXT 44 + XOPTUCODE	/* added in Igor Pro 2.0D83 */
#define TUFETCHSELECTEDTEXT 45 + XOPTUCODE	/* added in Igor Pro 2.0D83 */
#define TUSETSTATUSAREA 46 + XOPTUCODE		/* added in Igor Pro 2.0D83 */
#define TUMOVE_TO_PREFERRED_POSITION 47 + XOPTUCODE	// Added in Igor Pro 3.10B03.
#define TUMOVE_TO_FULL_POSITION 48 + XOPTUCODE		// Added in Igor Pro 3.10B03.
#define TURETRIEVE 49 + XOPTUCODE					// Added in Igor Pro 3.10B03.
#define TUNEW2 50 + XOPTUCODE						// Added in Igor Pro 3.13B01.


/* Standard file loader flag bit definitions */
#define FILE_LOADER_OVERWRITE 1				/* /O means overwrite */
#define FILE_LOADER_DOUBLE_PRECISION 2		/* /D means double precision */
#define FILE_LOADER_COMPLEX 4				/* /C means complex */
#define FILE_LOADER_INTERACTIVE 8			/* /I means interactive -- use open dialog */
#define FILE_LOADER_AUTONAME 16				/* /A means autoname wave (/N is equivalent to /O/A) */
#define FILE_LOADER_PATH 32					/* /P means use symbolic path */
#define FILE_LOADER_QUIET 64				/* /Q means quiet -- no messages in history */
#define FILE_LOADER_LAST_FLAG 64			/* This allows an XOP to use other bits in flag for its own purposes */


/* Miscellaneous #defines */
#define OPEN_FILE_DLOG 1256					/* DLOG for standard open file dialog */
#define SAVE_FILE_DLOG 1257					/* DLOG for standard save file dialog */
#define MENU_STRINGS 1101					/* STR# defining XOPs menu entry if any */
#define XOP_INFO 1100						/* XOPI resource of various XOP info */
#define XOP_WIND 1100						/* start XOP window resource IDs from here */
#define XOP_SUBMENU 100						/* start XOP menu resource IDs from here */
#define XOP_CMDS 1100						/* ID for XOPC resource describing commands */
#define XOP_MENUS 1100						/* ID for XOPM resource describing menus */
#define XSET_ID 1100						/* ID for XSET resource containing settings */

// These relate to an optional STR# 1101 resource in which IGOR may look for certain strings.
#define XOP_MISC_STR_ID 1101
#define XOP_MISC_STR_MENU_ID_INDEX 1		// Menu ID number if XOP is adding menu item via STR# 1101 method.
#define XOP_MISC_STR_MENU_ITEM_INDEX 2		// Menu item text if XOP is adding menu item via STR# 1101 method.
#define XOP_MISC_STR_HELPFILE_NAME_INDEX 3	// Name of XOP's help file (including ".ihf" index on Windows).

// This relates to the optional STR# 1160 resource which defines a target window type.
#define XOP_TARGET_WINDOW_TYPE_ID 1160
#define XOP_TARGET_WINDOW_SINGULAR_INDEX 1	// Human friendly singular name of the window type (e.g., Surface Plot).
#define XOP_TARGET_WINDOW_PLURAL_INDEX 2	// Human friendly plural name of the window type (e.g., Surface Plots).
#define XOP_TARGET_WINDOW_KEYWORD_INDEX 3	// Keyword used for window recreation macro declaration (e.g., SurfacePlot) or "" if none.
#define XOP_TARGET_STYLE_KEYWORD_INDEX 4	// Keyword used for window style macro declaration (e.g., SurfacePlotStyle) or "" if none.

// Status codes used with TellIgorWindowStatus
#define WINDOW_STATUS_DID_HIDE 1
#define WINDOW_STATUS_DID_SHOW 2
#define WINDOW_STATUS_ACTIVATED 3
#define WINDOW_STATUS_DEACTIVATED 4
#define WINDOW_STATUS_ABOUT_TO_KILL 5

#define TOPMAPHANDLE *(Handle *)0xA50		/* TopMapHandle global in Mac memory */

#include "XOPStructureAlignmentReset.h"		// Reset structure alignment to default.

#ifdef __cplusplus
}
#endif
