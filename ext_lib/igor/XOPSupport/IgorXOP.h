// IgorXOP.h -- Miscellaneous equates for interfacing XOP to Igor.

/*	These equates come from .h files used in compiling Igor and include
	various information that an XOP might need in communicating with Igor.
*/

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

// Miscellaneous
#define TUStuffHandle Handle
#define waveHndl Handle
#define DataFolderHandle Handle

// This is left over from the days when double's were unpredictable sizes. Some XOPs may still use DOUBLE.
#define DOUBLE double

#define CtoToolStr(str) c2pstr(str)						// Converts C to Pascal string.
#define TooltoCStr(str) p2cstr((unsigned char*) str)	// Converts C to Pascal string.


// From DataFolderBits.h. Used by Data Browser to determine events of interest to it.
#define kDFCB_NewChildFolder	1
#define kDFCB_KillChildFolder	2
#define kDFCB_RenameChildFolder	4
#define kDFCB_NewWave			8
#define kDFCB_KillWave			0x10
#define kDFCB_RenameWave		0x20
#define kDFCB_NewVariable		0x40
#define kDFCB_KillVariable		0x80
#define kDFCB_RenameVariable	0x100
#define kDFCB_NewString			0x200
#define kDFCB_KillString		0x400
#define kDFCB_RenameString		0x800
#define kDFCD_LockWave			0x1000					// AG27MAY03, Igor Pro 5


// From WinGenMacs.c. Used for DoWindowRecreationDialog callback.
enum CloseWinAction {
	kCloseWinCancel = 0,
	kCloseWinSave = 1,
	kCloseWinReplace = 2,
	kCloseWinNoSave = 3
};


// From Igor.h
#ifndef NIL
	#define NIL 0L
#endif
#ifndef FALSE				// Conditional compilation is needed because Metrowerks
	#define FALSE 0			// defines TRUE and FALSE in their MacHeaders.c.
#endif
#ifndef TRUE
	#define TRUE -1
#endif

#define MAXCMDLEN 400		// HR, 10/2/93 -- changed from 200 to 400 for Igor 2.0.

#define WAVE_OBJECT 1
#define WAVEARRAY_OBJECT 2
#define VAR_OBJECT 3
#define STR_OBJECT 4
#define STRUCT_OBJECT 5
#define XOPTARGWIN_OBJECT 5		// HR, 980714. Igor Pro 3.13B03.
#define GRAPH_OBJECT 6
#define TABLE_OBJECT 7
#define LAYOUT_OBJECT 8
#define PANEL_OBJECT 9			// HR, 10/2/93. Igor Pro 2.0.
#define NOTEBOOK_OBJECT 10
#define DATAFOLDER_OBJECT 11	// HR, 7/7/95. Igor Pro 3.0.
#define PATH_OBJECT 12			// HR, 7/28/95. Igor Pro 3.0. Symbolic path.
#define PICT_OBJECT 13			// HR, 7/28/95. Igor Pro 3.0. Picture.
#define ANNOTATION_OBJECT 14	// JP, 4/24/98. Igor Pro 3.13.
#define CONTROL_OBJECT 15		// JP, 4/24/98. Igor Pro 3.13.

//	#defines for identifying windows.
#define GRAF_MASK 1
#define SS_MASK 2
#define PL_MASK 4
#define PICT_MASK 8
#define MW_MASK 16
#define TEXT_MASK 32
#define PANEL_MASK 64
#define PROC_MASK 128
#define MOVIE_MASK 256
#define HELP_MASK 512
#define XOP_TARGET_MASK 4096			// HR, 980706: Added for XOP target windows.
#define ALL_MASK -1

#define FIRSTWINCODE CMDWIN
#define CMDWIN 1
#define WMDIALOGWIN 2				// HR, 11/19/92 -- was 10.
#define OLD_PROCWIN 2				// HR, 11/19/92 -- PROCWIN is now 10.
#define GRAFWIN 3
#define SSWIN 4
#define PLWIN 5
#define PICTWIN 6
#define MWWIN 7						// Notebook window (may be plain or formatted text).
#define TEXTWIN 8
#define PANELWIN 9
#define PROCWIN 10					// HR, 11/19/92 -- was 2.
#define MOVIEWIN 11
#define HELPWIN 12
#define HELPDLOGWIN 13
#define XOPWIN 14
#define XOPTARGWIN 15				// To group all XOP target windows together in the windows menu.
#define LASTWINCODE	XOPTARGWIN


/*	Name space codes.
	For use with UniqueName2 callback (Igor Pro 2.0 or later)
	and CheckName callback (Igor Pro 3.0 or later).
*/
#define MAIN_NAME_SPACE 1			// Igor's main name space (waves, variables, windows).
#define PATHS_NAME_SPACE 2			// Igor's symbolic path name space.
#define PICTS_NAME_SPACE 3			// Igor's picture name space.
#define WINDOWS_NAME_SPACE 4		// Igor's windows name space. Added in Igor Pro 3.0.
#define DATAFOLDERS_NAME_SPACE 5	// Igor's data folders name space. Added in Igor Pro 3.0.


// From IgorErrors.h.

#define BAD_EXPR_TERM			2001
#define ILLEGAL_POINT_VAL		2002
#define UNKNOWN_DEST_OR_BAD_CMD	2003
#define BAD_EQN_OR_CMD			2004
#define BAD_STR_EXPR			2005
#define STR_TOO_LONG			2006
#define BLD_MAC_ERR				2007
#define X_MAC_ERR				2008
#define MACRO_NESTING_OVERFLOW	2009
#define MACRO_XS_PARAMS			2010
#define MACRO_NO_LPAREN			2011
#define MAC_END_TOO_SOON		2012
#define MAC_PARM_NOT_DEFINED	2013
#define MAC_BAD_PARM_NAME		2014
#define F_BAD_KEYWORD			2015	// unused
#define FIF_ENDIF_MISSING		2016
#define FITTER_LOOP_MISSING		2017
#define F_NO_LPAREN				2018
#define FIF_IF_MISSING			2019
#define FITTER_ITTER_MISSING	2020
#define F_BREAK_NOITTER			2021
#define FDO_DO_MISSING			2022
#define REDUNDANT_PARAM_DEF		2023	// parameters names must be unique
#define NOT_PARAM_NAME			2024	// not a parameter name
#define EXPECT_WAVEPARAM_NAME	2025	// expected a wave name as a function parameter
#define EXPECT_KW_OR_OBJ_NAME	2026	// expected a keyword or an object name
#define EXPECTED_EQU			2027	// expected assignment operator: =, +=, -=, *= or /=
#define NOT_ALLOWED_IN_USER_FCTN 2028	// this function or operation is not available in user functions
#define UNUSED2029				2029	// "unused"
#define EXPECT_COLON_SEP		2030	// "expected : separator"
#define EXPECT_SUBTYPE			2031	// "expected subtype"
#define EXPECTED_EQU_ONLY		2032	// "expected '='"
#define FUNC_END_TOO_SOON		2033	// "unexpected end of function definition"
#define COMP_FUNC_ERROR			2034	// "Function compilation error"
#define COMP_MENU_ERROR			2035	// "Menu compilation error"
#define NO_MACROS_IN_FUNCTIONS	2036	// "Sorry, you can't invoke a macro . . ."
#define CMD_LINE_TOO_LONG		2037	// "The line is too long. Igor command lines are limited to 400 characters."
#define NO_PROMPT_IN_FUNCTIONS	2038	// "Sorry, you can't use Prompt in a function..."

#define NOMEM 1							//  out of memory
#define NOWAV 2							// expected wave name
#define SYNERR 3						// syntax error 
#define NOCOMMA 4						//  expected comma
#define BADAXIS 5						// expected name of active axis
#define BADUNITS 6						//  expected axis units
#define BADNUM 7						//  expected number
#define NOGRAF 8						// there are no graphs
#define BADFLG 9						// unknown flag 
#define BADNPNTS 10						// "the number points in a wave must be between 0 and 20 million."
#define NOEQUALS 11						// missing equal sign in flag (/N=val)
#define NOTPOW2 12						// wave length must be power of 2 for this operation
#define CANT_OVERWRITE_FILE 13			// "file of wrong type -- can not be overwritten"
#define NO_NLFUNC 14					// expected fitting function
#define PNTS_INCOMPATIBLE 15			//  incompatible waveform size
#define NT_INCOMPATIBLE 16				//  incompatible number types
#define NT_FNOT_AVAIL 17				//  can't do function on this number type
#define BADRANGE 18						//  insufficient range
#define NOREALIFFT 19					// can't do IFFT on real data 
#define WNAME_SYNERR 20					// expected wave name
#define NAME_TOO_LONG 21				// name or string too long
#define NO_TERM_QUOTE 22				// expected terminating quote
#define BAD_NAME 23						// ill-formed name
#define NO_SVAR_OR_FCTN 24				// expected string variable or string function
#define NAME_VAR_CONFLICT 25			// name already exists as a variable
#define NAME_CMD_CONFLICT 26			// name already exists as a command
#define NAME_WAV_CONFLICT 27			// name already exists as a waveform
#define NAME_FUNC_CONFLICT 28			// name already exists as a function
#define NAME_MAC_CONFLICT 29			// name already exists as a macro
#define NO_LOG_AXIS 30					// can't do log axis on this axis type
#define NO_SYMFLDR 31					// no symbolic folder of that name
#define NO_SF_DELETE 32					// can't delete special folder
#define TOO_MANY_SFVARS 33				// No more symbolic folders are available.
#define NO_SF_CHANGE 34					// Can't change a special symbolic path.
#define NO_MULT_BIN_SAVE 35				// Can't save multiple waves to a single binary file.
#define BAD_BINARY_FILE 36				// bad or damaged binary file
#define BAD_SPECIAL_HEADER 37			// <<unused>> -- can't find keyword in special file
#define BAD_BIN_VERSION 38				// Bad or incompatible binary version
#define NOKILL_WAVE_IN_USE 39			// can't kill a wave that is part of a graph or table
#define TOO_MANY_PARAMS 40				// Too many parameters on command line.
#define PF_NOCONV 41					// failure to converge during fit
#define BAD_ATTACH_CODE 42				// Improper letter pair for attach code
#define ORN_OFFSET_OUT_OT_RANGE 43 		// Offset out of range
#define NO_ORNS 44						// there are no tags or textboxes
#define NO_NAMED_ORN 45					// there are no tags or textboxes of that name
#define NO_SUCH_FONT 46					// Font not available
#define SSBADCOLUMNNAME 47				// Invalid column name
#define SSNOCOLUMNS 48					// expected column name
#define NO_TARG_WIN 49					// No target window.
#define EXPECT_MORE_DATA 50				// Expecting more data
#define EXPECT_POS_NUM 51 				// expecting a positive non-zero number
#define HIST_APPEND_BAD_DEST 52			// Destination wave must be single precision real for histogram append
#define BAD_MIRROR_AXIS 53				// Mirror axis only available for left or bottom axies when right or top are unused
#define WAVE_NOT_ON_GRAPH 54			// Wave is not on graph
#define NO_REMOVE_CONTROL_WAVE 55 		// can't remove a wave that controls an axis*/
#define MISMATCH_LISTASSIGNMENT 56 		// "mismatch points in wave[a,b]={n,...,n} expression"
#define USER_ABORT 57					// user abort
#define SINGULAR_MATRIX 58				// singular matrix or other numeric error.
#define NO_DATA_IN_FILE 59				// no data was found in file.
#define LINE_TOO_LONG_IN_FILE 60		// file contains a line that is too long.
#define TOOMANYCOLS_IN_FILE 61			// file contains too many columns.
#define WAVE_NAME_TOO_LONG 62			// Wave names can't exceed 31 characters
#define BAD_CHAR_IN_NAME 63				// names must start with a letter and contain letters, digits or '_'
#define TOO_MANY_MAKE_WAVES 64			// More than 100 waves specified in one make command
#define BAD_IGOR_TEXT_FILE 65			// Bad Igor text file. (missing keyword IGOR)
#define BEGIN_BUT_NO_WAVES 66			// Bad Igor text file. (BEGIN keyword but no WAVES)
// The following LA_XXX's are on-the-spot error reports from loadspecial.c.
#define LA_FLAG_PLACEMENT 67			// Flag in wrong place
#define LA_MAX_WAVES 68					// Limit of 12 waves per group exceeded
#define LA_BAD_SYMB_WAVES 69			// Bad symbol on WAVES line
#define LA_NO_NAMES 70					// no wave names found
#define LA_NO_END 71					// missing END or unknown symbol
#define LA_NO_IMAG 72					// bad or missing imaginary number
#define LA_TOO_FEW_PNTS 73				// too little data (less than 2 points)
#define EXPECT_AS_OR_END 74				// Expected 'as', <cr> or ';'
#define DISP_SYNERR 75					// Expected wave name, 'vs', 'as', <cr> or ';'
#define APPEND_SYNERR 76				// Expected wave name, 'vs', <cr> or ';'
#define EXPECT_END 77					// Expected ';' or <cr>
#define EXPECT_WAVE_OR_END 78			// Expected wave name, ';' or <cr>
#define EXPECT_RPAREN 79				// Expected ')' 
#define EXPECT_KEYWORD_OR_NAME 80		// Expected key word or name 
#define DUP_NO_OVERWRITE 81				// Can't overwrite self 
#define EXPECT_LPAREN_OR_EQUAL 82		// Expected '(' or '='
#define UNKNOWN_KEYWORD 83				// Unknown keyword
#define BAD_IRANGE 84					// Expected number between ^2 and ^3
#define PROC_TOO_BIG 85					// The procedure file is getting too big
#define TEXT_TOO_BIG 86					// The window can't hold any more text
#define BAD_FONT_NAME 87				// Font not available
#define NO_PARAMETRIC_X 88				// Can't use a parametric X wave here
#define BAD_CHECKSUM 89					// bad checksum
#define SSMAXDOCSERR 90					// max table documents already opened
#define NO_MORE_PARAMS 91				// no more parameters are expected
#define EXPECTED_XY 92					// expected 'x' or 'y'
#define ONLY_IN_MACRO 93				// this command is only used in macros
#define EXPECTED_VARNAME 94				// expected variable name
#define EXPECTED_STRINGVARNAME 95		// expected string variable name
#define NUM_EXPECTED_NORM 96			// expected comma or end of command (<cr>, ';', or |)
#define NUM_EXPECTED_CRP 97				// expected comma or ')'
#define EXPECT_LPAREN 98				// expected '('
#define EXPECTED_TARG_NAME 99			// expected name of a target window
#define NO_STACK 100					// out of stack space (too much macro recursion)
#define NO_COMPLEX_WAVE 101				// can't use a complex wave here
#define BAD_INDEX 102					// illegal index value
#define NO_LEVEL_FOUND 103				// level crossing not found
#define GRAFWIN_TOO_BIG 104				// graph dimensions too big

#define EXPECT_GRAFWIN 105				// expected name of graph window
#define TOOMANY_OBJECTS 106				// too many graphs, tables or other objects
#define EXPECTED_XOP 107				// expected XOP
#define EXPECTED_XOP_PARAM 108			// expected XOP parameter
#define UNKNOWN_XOP 109					// unknown XOP
#define XOP_OBSOLETE 110				// XOP is incompatible with this version of Igor
#define XOP_HAS_NO_CMD 111				// this XOP has no command associated with it
#define XOP_CANT_INIT 112				// XOP can't initialize itself
#define XOP_INCOMPATIBLE 113			// XOP incompatible with system software or hardware
#define EXPECT_LPAREN_BRACKET 114		// expected '(' or '['
#define INVALID_FORMAT 115				// format string is invalid
#define XOP_NOT_FOUND 116				// can't find file containing XOP
#define EXPECTED_OBJECT 117				// expected graph, table, picture, or textbox name
#define OBJECT_NOT_IN_LAYOUT 118		// object is not in layout
#define NAME_USER_FUNC_CONFLICT 119		// name already exists as user function
#define BAD_HOLD_STRING 120				// hold string length doesn't match number of parameters
#define NO_USER_FCTN 121				// no user defined function of that name exists
#define BAD_FIT_FCTN_NTYPE 122			// functions for fitting must return a single or double precision scalar result
#define BAD_FIT_FCTN_FORMAT 123			// fitting function does not have required form
#define USER_FUC_TOOMANY_PNTS 124		// coefficient wave implys too many fit coefficients (>50)
#define BAD_GROUT 125					// grout must be between 2 and 72 points
#define BAD_TILE_RECT 126				// the tiling window is too small
#define EXPECTED_LAYOUT_OBJECT 127		// expected name of object in layout
#define NO_LAYOUTS 128					// there are no page layouts
#define EBAR_WAVES_MISSING 129			// both positive & negative error bar waves are missing
#define BAD_LAYOUT_EXPAND 130			// magnification must be .25, .5, 1, or 2
#define EXPECTED_LAYOUT 131				// expected name of page layout window
#define NO_PRINT_RECORD 132				// can't open print record (check print driver)
#define ONLY_GRAF 133					// this operation is for graphs only
#define TIME_OUT_READ 134				// timeout while reading data
#define TIME_OUT_WRITE 135				// timeout while writing data

#define BAD_REFNUM 136					// there is no open file with this reference number
#define EXPECTED_TIME 137				// expected time in hh:mm:ss format
#define NOT_TEXT_FILE 138				// non-text files can not be overwritten
#define TOO_MANY_PARAMETERS 139			// too many parameters for this command
#define BAD_WIN_TITLE 140				// expected window title
#define NUM_EXPECTED_RBRACK 141			// expected ']'
#define EXPECT_VERT_AXIS_NAME 142		// "expected vertical axis keyword"
#define EXPECT_HORIZ_AXIS_NAME 143		// "expected horizontal axis keyword"
#define EXPECT_SIZE_KEYWORD 144			// "expected graph size keyword"
#define INDEX_OUT_OF_RANGE 145			// index out of range
#define EXPECTED_FONT 146				// expected font name
#define EXPECTED_MACRO_NAME 147			// expected macro name
#define EXPECT_0_OR_LCBRACE 148			// "expected 0 or '{'"
#define EXPECT_LBRACK 149				// "expected '['

#define EXPECT_LCBRACE 150				// "expected '{'
#define NUM_EXPECTED_RCBRACE 151		// expected '}'

#define EXPECTED_MENU_NAME 152			// expected menu name
#define EXPECTED_MENU_KEYWORD 153		// expected "menu item", SubMenu or End
#define TOO_MANY_SUBMENUS 154			// there are too many submenus
#define CANT_REMOVE_WAVE_TWICE 155		// can't remove the same wave twice
#define FUNCTION_CHANGED_BACK 156		// "A function name cannot be changed in this dialog. It has been changed back to Ò^0Ó."
#define NEED_ANNOTATION_NAME 157		// expected /N=name
#define FLAG_MUST_BE_FIRST 158			// this must be first flag in command
#define NUM_EXPECTED_CRB 159			// expected comma or '}'
#define NAME_WIN_CONFLICT 160			// name already exists as a window name

#define NUM_EXPECTED_CRBRACK 161		// "expected comma or ']'"
#define EXPECTED_GRAPH_OR_TABLE 162		// "expected name of a graph or table"
#define IGOR_OBSOLETE 163				// XOP requires a later version of Igor
#define EXPECT_LINK 164					// "expected link field in options string"
#define EXPECT_LINK_ARG 165				// "expected link argument field in options string"
#define NAME_STRFUNC_CONFLICT 166		// "name already exists as a string function"
#define EXPECT_WAVELIST_LINK 167		// "expected wavelist option link: 'WIN:'"
#define NO_RENAME_BIVARS 168			// "can't rename built-in variables"
#define EXPECT_WAVE_OR_VAR 169			// "expected wave or variable"

#define BAD_OPTIONS_STRING 170			// "inappropriate options string"
#define	DUPLICATE_NAMES 171				// "two or more names are identical"
#define CANT_LOAD_BINARY_FROM_CLIP 172 	// "can't load binary from clipboard"
#define DEMO_ERR 173					// "This is demo..."

#define FIRST_IGOR2_ERR 174				// Start of Igor Pro 2.0 error codes.

#define USR_BKG_ERROR 174				// "Background user function returned error"
#define BKG_NO_TASK_RUNNING 175			// "No background task to control"

#define DASHSPEC_TOO_LONG 176			// "Dash specification too long"

#define BAD_BOUND_EQN 177				// "Wrong format for dependency formula"
#define EXPECT_EQN_DEF 178				// "Expected ':='"
#define NO_PEAK_FOUND 179				// "Peak center not found"
#define EQN_TOO_LONG 180				// "Formula is too long"
#define NO_LOCAL_EQN_VARS 181			// "Local variables can't have dependency formula"

#define ONLY_GRAPH_OR_LAYOUT 182		// this operation is for graphs or layouts only

#define BAD_FPOS_SET 183				// "can't set file mark past end-of-file"
#define NOGRAF_OR_PANEL 184				// "there are no graphs or panels"

#define USING_NULL_STRVAR 185			// "attempt to use a null string variable"
#define FEATURE_NOT_AVAIL 186			// "This feature is not yet available because I don't know how to do it!"

#define INCOMPATIBLE_FLAGS 187			// "incompatible flags"
#define EXPECT_TWO_EXPRS 188			// "expected two expressions"

#define NO_LOCAL_IN_BOUND_EQN 189		// "Can't use local variables in dependency formula"
#define WAVE_TYPE_INCONSISTENT 190		// "Inconsistent type for a wave variable"

#define BAD_FLAG_NUM 191				// "Flag usage is '/N=number or /N=(expression)"

#define BAD_IRANGE_ODD 192				// "Expected odd number between ^2 and ^3"
#define BAD_EXPECTED_RANGE 193			// "expected ^1 between ^2 and ^3" 
#define BAD_ODD_EXPECTED_RANGE 194		// "expected odd ^1 between ^2 and ^3"

#define BAD_XWave 195					// "X data does not match Y data (length or number type)"
#define BAD_WWave 196					// "Weight data does not match Y data (length or number type)"
#define BAD_DWave 197					// "Destination wave does not match Y wave (length or number type)"

#define BAD_LEVELS_SPEC 198				// "missing level specification value(s)"
#define NO_KILL_VAR_OR_STR 199			// "Can't kill variables or strings from a user function"
#define SIZE_OUT_OF_RANGE 200			// "Size out of range"
#define POSITION_OUT_OF_RANGE 201		// "Position out of range"
#define FONTSIZE_OUT_OF_RANGE 202		// "Font size out of range"
#define SHORTER_THAN_COEFS 203			// "Length of wave ^0 can not be less than the coefficient wave"
#define CANT_MIX_CMPLX_WITH_REAL 204	// "Can't combine complex wave ^0 with real wave ^1"
#define CANT_MIX_REAL_WITH_CMPLX 205	// "Can't combine real wave ^0 with complex wave ^1"
#define NAME_ALREADY_IN_USE		 206	// "Name ^0 is already in use"
#define EXPECTED_CONTROL_NAME	 207	// "Expected control name"
#define NO_CHANGE				 208	// "(no change)"
#define EXPECT_MACRO_PROC_FUNC	 209	// "Expected Macro, Proc, or Function keyword"
#define EXPECTED_FUNCTION_NAME	 210	// "Expected function name"
#define TOO_FEW_PARAMETERS		 211	// "Too few parameters"
#define EXPECTED_SUBTYPE		 212	// "Expected subtype"
#define EXPECTED_END_OR_COLON	 213	// "Expected ':' or end of line"
#define WARN_WRONG_CTL_SUBTYPE	 214	// "Warning: wrong subtype for this control, should be Ò : ^0 Ó "
#define WARN_NO_CTLSUBTYPE		 215	// "(Note: optional subtype Ò : ^0 Ó is missing)"	
#define WARN_SUBTYPE_FIXED		 216	// "(Note: optional subtype Ò : ^0 Ó has been added to this procedure)"
#define WARN_PROC_NOT_FOUND		 217	// "Warning: can't find Ò^0Ó in any procedure window"
#define WARN_PROC_EDITED		 218	// "(no change to control, procedure already changed)"
#define NOTE_PROC_EDITED		 219	// "(Note: procedure ^0 has previously been changed)"
#define EXPECTED_NUMEXPR		 220	// "Expected numeric expression"
#define NOQUICKTIME		 		 221	// "QuickTime not present"
#define ERR_MOVIE_ALREADY_OPEN	 222	// "A movie file is already open"
#define ERR_MOVIE_NOT_OPEN		 223	// "No movie file is open"
#define BADSOUNDWAVE			 224	// "Bad sample rate or amplidude for audio wave"
#define NO_WIND_STYLE_MACRO		 225	// "No style macro for this type of window"
#define WRONG_CONTROL			 226	// "^0 is not a ^1 control"
#define EXPECTED_NAME			 227	// "Expected name"
#define RFLAG_NEEDS_CFLAG		 228	// "/C flag must preceed /R flag"
#define ORN_RENAMED				 229	// "(no changes, except ornament already renamed to Ò^0Ó.)"
#define CROSS_AXIS_MISSING		 230	// "Crossing axis not found"
#define CROSS_AXIS_NOT_PERP		 231	// "Crossing axis is not perpendicular"
#define EXPECTED_FIFO			 232	// "expected name of FIFO"
#define FIFO_IN_USE				 233	// "FIFO in use by XOP"
#define NOT_WHILE_FIFO_RUNNING	 234	// "operation not allowed while FIFO is running"
#define NO_FIFO_CHANS			 235	// "no FIFO channels have been defined"
#define FIFO_STOPPED			 236	// "FIFO is not running"
#define NO_SUCH_FIFO_CHAN		 237	// "no FIFO channel of that name"
#define FIFO_OVERFLOW			 238	// "FIFO overflow (disk did not keep up)"
#define WRONG_NUM_CHANS			 239	// "FIFO has a different number of channels"
#define NO_SUCH_ChartCHAN		 240	// "no chart channel of that number"
#define PATH_BUT_NO_FILE		 241	// "/P flag requires file name argument"
#define FILE_NOT_FOUND			 242	// "File not found"
#define EXPECTED_COMPLEX_NUM	 243	// "Expected complex number"
#define EXPECTED_FUNCTION_KEY	 244	// "Expected Function keyword"
#define EXTRA_MACRO_TEXT		 245	// "Extra text in macro, proc, or function"

#define BAD_FILE_TYPE_SPEC		 246	// "A file type specifier must have exactly 4 characters"
#define BAD_FILE_CREATOR_SPEC	 247	// "A file creator specifier must have exactly 4 characters"
#define PATH_TOO_LONG			 248	// "The path to file is too long"
#define FILE_OPEN_READ_ONLY		 249	// "The file ^1 is already open read-only"
#define FILE_OPEN_WRONG_NAME	 250	// "The file ^1 is already open but with the window name ^2"
#define FILE_OPEN_WRONG_TYPE	 251	// "The file ^1 is already open but as a ^2 file"
#define MENU_ITEM_HAS_NO_SUBMENU 252	// "This menu item has no submenu"
#define NO_MACRO_OR_FUNCTION	 253	// "There is no procedure named ^0"
#define CANT_APPEND_TO_THIS_MENU 254	// "Can't add items to this menu"

#define EXPECTED_PICT			 255	// "Expected picture name"
#define CANT_DRAW_HERE			 256	// "Can't draw here"
#define LAYER_NOT_AVAIL			 257	// "Layer not available"
#define NO_DRAW_OBJ				 258	// "No draw object"
#define NO_DOLLAR_HERE			 259	// "Can't use $ here (compiling)"
#define NO_OPEQU_LIST			 260	// "Can't use op= with {n,..,n} assignments"

#define EXPECTED_NOTEBOOK		261		// "Expected name of a notebook window"
#define NB_LOCS_BACKWARDS 		262		// "Invalid notebook selection: the end location is before the start location"
#define NB_STARTLOC_INVALID		263		// "Invalid notebook selection: the start location is out of bounds"
#define NB_ENDLOC_INVALID		264		// "Invalid notebook selection: the end location is out of bounds"

#define EXPECTED_KEYWORD_OR_END 265		// "Expected ',keyword' or end of command (<cr>, ';', or |)"
#define EXPECTED_GRAPHIC_OBJ_NAME 266	// "Expected name of graph, table, layout or picture"
#define NB_BAD_GRAPH_DIMENSIONS 267		// "The graph width and height must be between 50 and 8200"
#define NB_BAD_LAYOUT_RECT		268		// "The layout rectangle is unreasonable"
#define EXPECTED_STRING_EXPR	269		// "Expected string expression"
#define NB_NO_RULER_SPECIFIED	270		// "No ruler specified (use ruler=rulerName or newRuler=rulerName)"
#define TOO_MANY_TABS			271		// "No more than 20 tabs are allowed"
#define BAD_TAB					272		// "Illegal tab value"
#define EXPECTED_SELECTION_KW	273		// "Expected notebook selection keyword"
#define EXPECTED_PATH			274		// "Expected name of a symbolic path"
#define NB_UNKNOWN_SPECIAL_CHAR_TYPE 275 // "Unknown special character code"
#define FUNCTION_TOO_LONG		276		// "The function is too long"
#define BAD_GRAPH_PREFS			277		// "Bad Graph Preferences"
#define BAD_TABLE_PREFS			278		// "Bad Table Preferences"
#define BAD_LAYOUT_PREFS		279		// "Bad Layout Preferences"
#define BAD_PANEL_PREFS			280		// "Bad Panel Preferences"
#define BAD_NOTEBOOK_PREFS		281		// "Bad Notebook Preferences"
#define BAD_PROCEDURE_PREFS		282		// "Bad Procedure Preferences"
#define BAD_COMMAND_PREFS		283		// "Bad Command Window Preferences"
#define BAD_WINDOWS_PREFS		284		// "Bad Windows Preferences"
#define BAD_PALETTE_PREFS		285		// "Bad Color Palette Preferences"
#define BAD_DASHED_PREFS		286		// "Bad Dashed Lines Preferences"
#define BAD_MISC_PREFS			287		// "Bad Miscellaneous Preferences"
#define BAD_HEADER_FOOTER_PREFS	288		// "Bad Header/Footer Preferences"

#define NO_DATA_IN_CLIP			289		// "no data was found in the clipboard"

#define PICT_ALREADY_EXISTS		290		// "a picture by that name already exists"
#define EXPECTED_PICT_NAME		291		// "expected the name of a picture"
#define NO_PICT_IN_CLIP			292		// "there is no picture in the clipboard"
#define RESOURCE_NOT_FOUND		293		// "resource not found"

#define TOOMANY_WLASSIGN_PNTS 	294		// "too many points in wave={n,..,n} expression"
#define NO_KILL_SELF			295		// "A window hook procedure tried to kill its own window"
#define FIT_COEFF_DIFFERENT		296		// "The coefficient wave is not the same number type as the data. Use Redimension."
#define KN_NO_EQN				297		// "Kn variables can't have dependency formulas."
#define STD_AXISNAME_WRONG_EDGE	298		// "Can't use standard axis name on different edge."
#define EXPECT_VERT_AXIS		299		// "Expected vertical axis."
#define EXPECT_HORIZ_AXIS		300		// "Expected horizontal axis."
#define NOT_FIFO_FILE			301		// "Not a FIFO file."
#define BAD_FIFO_FILE_VERS		302		// "Bad FIFO file version."
#define CORRUPT_FIFO_FILE		303		// "Corrupt FIFO file."

#define CANT_KILL_PICT			304		// "can't kill a picture that is used in a graph, layout or panel"
#define EXPECT_WAVE_VAR_STR		305		// "expected name of wave, variable, or string"
#define NAME_PICT_CONFLICT		306		// "name already exists as a picture"
#define BAD_PICT_VERSION		307		// "picture version is not valid"

#define CMD_CAN_NOT_BE_COMPILED	308		// "Sorry, this operation is not allowed in a function. It is allowed in a macro."
#define EXPECTED_TABLE			309		// "expected name of a table"
#define DELIM_RADIX_CONFLICT	310		// "comma can not be both a delimiter and the decimal character"
#define IT_EXPECTED_BEGIN		311		// "expected BEGIN keyword after WAVES in Igor Text file
#define IT_UNKNOWN_KEYWORD		312		// "unknown keyword in Igor Text file"
#define LF_BAD_COLUMN_NAME_LINE	313		// "column names must be before data"
#define LF_BAD_LINE_OR_COL		314		// "line and column numbers must be ³ 0"

#define MENU_HELP_MISPLACED 	315		// "the help must appear on the line after the menu item string expression"
#define MENU_KEYWORD_UNKNOWN 	316		// "the only keyword allowed after 'Menu' is 'dynamic'"

#define NO_IGOR_OBJECT_INFO_PICCOMMENT 317	// "The picture does not contain Igor object information"
#define BAD_IGOR_OBJECT_INFO_PICCOMMENT 318	// "The picture contains Igor object information that this version of Igor does not understand"

#define NO_TABLES				319		// "There are no tables",
#define NO_XFUNC				320		// "No external function of that name exists"
#define NO_USERFUNC_OR_XFUNC	321		// "No user or external function of that name exists"

#define NO_COEFS_WAVE			322		// "You need a coefficients wave for fitting to a user-defined function"
#define BAD_POLY_NTERMS			323		// "Expected a number of polynomial terms from 3 to 20"
#define TOO_MANY_WAVES			324		// "Too many waves - 100 is the maximum"

#define NAME_PATH_CONFLICT 		325		// "Name already exists as symbolic path"
#define RENAME_CONFLICT			326		// "You have already renamed Ò^0Ó to Ò^1Ó."

#define NO_CHANGE_WAVE_INUSE	327		// "Can't change a wave in the middle of a wave assignment."
#define NO_WAVE_X_DEST			328		// "Can't use wave(x)= in a function. Use x2point and wave[p]= instead"
#define EXPECT_MARGIN			329		// "Expected margin keyword (left,right,top or bottom)"
#define NULL_WAVE_OP			330		// "Attempt to operate on a NULL or missing wave"
#define NAME_IS_RESERVED_KW		331		// "Name is a reserved keyword"
#define NOCOMPILE_APPEND		332		// "Can't compile Append. Use AppendToGraph or AppendToLayout or AppendToTable"
#define NOCOMPILE_REMOVE		333		// "Can't compile Remove. Use RemoveFromGraph or RemoveFromLayout or RemoveFromTable"
#define AXISENAB_RANGE_ERR		334		// "Axis enable settings out of range: must be between 0 and 1 and start < stop"
#define NEED_COMPILE			335		// "The procedure window(s) need to be compiled. Perhaps auto-compile is off."
#define NOKILL_OBJ_IN_FCTN		336		// "Can't kill a wave or variable that is part a user function or dependency expression."
#define TAG_FUNCTION_ERROR		337		// "A tag access function is only valid while a tag is being drawn."
#define TRACE_SPECIFED_TWICE	338		// "A trace was specifed twice."

#define WIN_TITLE_BAD_LENGTH	339		// "Window titles can be 1 to 40 characters long"
#define UNKNOWN_LUMP_REC_VERSION 340	// "This version of Igor can't handle this packed experiment file"
#define CANT_UNTIL_CHOOSE_FROM_LIST 341	// "Select an item in the list"
#define XOP_RESOURCES_MISSING	342		// "The XOP file Ò^1Ó is missing required resources. . ."
#define XOP_NEEDS_FPU			343		// "This XOP can't be loaded because it requires a math coprocessor."

#define NUM_BINS_MUST_BE_TWO 344			// "Histogram requires at least two bins"
#define SRC_AND_DEST_MUST_BE_DIFFERENT 345	// "Source and destination waves must be different"
#define BIN_WIDTH_CANT_BE_ZERO 346			// "Histogram bin width can't be zero"
#define BIN_PARAMS_NAN_OR_INF 347			// "Histogram bin start or bin width is a NaN or an INF"

#define RULERS_MUST_BE_DIFFERENT 348		// "The rulers must be different"
#define EXPECTED_GRAPH_OR_LAYOUT 349		// "expected name of a graph or layout"
#define SAFE_SAVE_DISK_FULL 350				// "The save could not be done because the disk is full"
#define DIRECT_XFUNC_CANT_DO_CALLBACK 351	// "XFUNC programmer error: a direct XFUNC is not allowed to do a callback to Igor."

#define INCLUDE_BAD_FILE_TYPE 352			// "#included files must be of type TEXT"
#define BAD_INCLUDE_SPEC 353				// "the #include file specification is bad"
#define INCLUDE_FILE_NOT_FOUND 354			// "include file not found"
#define READONLY_PROCEDURE 355				// "This is a read-only procedure. Change the name to create a new procedure."

#define TU_BAD_PARAGRAPH 356				// "The TU paragraph is out of range."
#define TU_BAD_LOC_CHAR_POS 357				// "The TU location character position is out of range."
#define TU_BAD_LOC_ORDER 358				// "The first TU location is AFTER the second TU location."

#define XOP_RECURSION_ATTEMPTED 359			// "The XOP has attempted recursion"
#define INCLUDED_FILE_OUTDATED 360			// "The included procedure file Ò^1Ó needs to be updated"

#define BUTTON_NEEDS_COMPILED 361			// "You need to compile the procedure windows before you use controls."
#define NO_BUTTON_PROC 362					// "Cannot find control's action procedure Ò^0Ó

#define TOO_LONG_FOR_CMD_LINE 363			// "Too long to fit on command line"
#define CLICK_AUTO_DUMMY 364				// "Click ÒSet to Auto ValuesÓ button"
#define NEED_N_DIGITS 365					// "You need at least ^2 digits to properly label each tick increment."
#define START_EXCEEDS_END 366				// "Start value must be less than end value"
#define PATH_IN_USE 367						// "The symbolic path is in use and can't be killed."

#define XOP_REQUIRES_PPC 368				// "This XOP will run on a PowerMac but not on this 680x0 Mac."
#define CODE_FRAGMENT_LOAD_ERROR 369		// "A system error (%d) occurred while loading the code fragment \"%s\"."

#define XFUNC_BAD_NT 370					// "This XFUNC was compiled with an illegal number type (must be double)."

#define DO_WINDOW_FROM_FUNCTION 371			// "DoWindow/R requires recompiling functions and thus can't be called from a function."

#define FIRST_IGOR3_ERR 372					// Start of Igor Pro 3.0 error codes.

#define NO_TEXTNUMERIC_WAVE_OVERWRITE 372	// "Can't overwrite a text wave with a numeric wave or vice versa."
#define NEED2PNTS_FOR_THIS_OP 373			// "This operation requires a wave with at least two points."
#define NO_TEXT_OP 374						// "This operation does not work on text waves."
#define NODATA_IN_DIM 375					// "There is no data allocated in this dimension."
#define DIMENSION_MISMATCH 376				// "Mismatch between actual and specified number of dimensions."
#define BAD_INCREMENT 377					// "Increment value is less than 1."
#define ZERO_DATA_IN_WAVE 378				// "The wave has zero data allocated."
#define INCONSISTENT_DIMENSIONS 379			// "Inconsistent number of items in a dimension."
#define BAD_DIMENSION_NUMBER 380			// "Dimension number out of range (0-3)."
#define EXPECT_DIM_LABEL 381				// "Expected a dimension item label (literal, not a string)."
#define NOSUCH_DIM_LABEL 382				// "Couldn't find the given dimension item label."
#define NO_PARENT_DATAFOLDER	383			// "Tried to access the parent data folder of the root or of a non-existent data folder."
#define NO_CHILD_DATAFOLDER		384			// "No child data folder of that name exists."
#define CANT_APPEND_DIF_X_TO_CAT 385		// "Can't append different x-wave to category axis."
#define CANT_APPEND_TEXT_X_TO_NONCAT 386	// "Can't append text x-wave to non-category axis."
#define TOO_MANY_SORT_SRC_WAVES	387			// "Too many sort key waves specified. Maximum is 10."
#define CANT_PUT_DF_IN_SELF		388			// "Can't move a data folder into itself or into a child folder."
#define CANT_RENAME_ROOT		389			// "Can't rename the root data folder."
#define EXPECT_DATAFOLDER_NAME		390		// "Expected name of a data folder"
#define NO_GLOBALS_IN_FUNCTIONS		391		// "Must use WAVE, NVAR & SVAR to access global waves, variables & strings with #pragma rtGlobals=1."
#define EXPECT_SQUARE_MAT		392			// "Expected a square matrix."

#define NO_ROOT_DATAFOLDER		393			// "A non-existent data folder was referenced while accessing a child data folder."
#define CANT_FIND_FOLDER		394			// "Can't find data folder."
#define CANT_KILL_PARENTFOLDER_OF_CURRENT 395	// "Can't kill a data folder that contains the current data folder."
#define FOLDER_NAME_EXISTS		396			// "A data folder of that name already exists at this level."
#define EXPECT_COMPILER_DIRECTIVE	397		// "Expected a compiler directive."
#define UNKNOWN_COMPILER_DIRECTIVE	398		// "Unknown compiler directive."
#define EXPECT_PRAGMA_KW		399			// "Expected a pragma keyword."
#define UNKNOWN_PRAGMA_KW		400			// "Unknown pragma keyword."
#define GVAR_TYPE_INCONSISTENT	401			// "Inconsistent type for a global variable reference."

#define OBSOLETE_SCRAP 402					// "The clipboard contents are not compatible with this version of Igor."
#define NUMERIC_WAVE_CANT_HAVE_TEXT_FORMAT 403	// "A numeric wave can not have a text format."
#define TEXT_WAVE_CANT_HAVE_NUMERIC_FORMAT 404	// "A text wave can not have a numeric or date/time format."

#define BAD_CHAR_IN_WAVE_NAME 405			// "A wave name can't contain any of the following: ' \" ; or : "*/
#define BAD_COLORINDEX_WAVE 406				// "Expected a matrix wave containing 3 columns with red, green, and blue values."
#define EXPECT_IMAGE_NAME 407				// "Expected the name of an image in the top graph."
#define EXPECT_MATRIX 408					// "Expected a 2D (matrix) wave."
#define EXPECT_VS_OR_END 409				// "Expected 'vs' keyword or end of command."
#define EXPECT_COLORTABLE_NAME 410			// "Expected name of a color table."

// These errors would normally be generated by a buggy XOP but could also be generated by a buggy Igor.
#define UNKNOWN_WAVE_ACCESS_MODE 411		// "An attempt was made to access a wave using an incompatible access mode."
#define NUMERIC_ACCESS_ON_TEXT_WAVE 412		// "An attempt was made to treat a text wave as if it were a numeric wave."
#define TEXT_ACCESS_ON_NUMERIC_WAVE 413		// "An attempt was made to treat a numeric wave as if it were a text wave."
#define MD_WAVE_BAD_INDEX 414				// "An invalid index was used to access a wave."

#define CONTOUR_EXPECTED_XY_WAVES 415		// "Expected \"vs {xwave,ywave}\"."
#define EXPECT_1DZ_WAVE 416					// "Expected a 1D (single column) contour data (z) wave."
#define EXPECTED_CONTOUR_XYZMATRIX 417		// "Expected a matrix wave containing 3 columns with X, Y, and Z values, or \"zwave vs {xwave,ywave}\"."
#define EXPECTED_1D_XYWAVE 418				// "Expected a 1D (single column) xwave or ywave in \"vs {xwave,ywave}\"."
#define CONTOUR_SHORT_XWAVE 419				// "xwave in \"vs {xwave,ywave}\" has too few rows for contour data rows."
#define CONTOUR_SHORT_YWAVE 420				// "ywave in \"vs {xwave,ywave}\" has too few rows for contour data columns."
#define CONTOUR_MISMATCHED_WAVES 421		// "XYZ waves must all have the same length."

// These errors are used by GetSelection. Prior to Igor Pro 2.5, they were equated to NOMEM because HR forgot to create real error messages for them.
#define EXPECTED_WINTYPE 422						// "Expected a window type keyword: graph, table, layout, notebook or panel."
#define EXPECTED_GRAPH_TABLE_LAYOUT_NOTEBOOK 423	// "Expected one of these window type keywords: graph, table, layout, notebook."
#define EXPECTED_TABLE_WIN 424						// "Expected name of table window"
#define EXPECTED_PANEL 425							// "Expected name of panel window"

#define NO_ROW_COL_LABELS_POSITIONS 426		// "Reading of row or column labels and positions is supported only for delimited-text matrix loading."
#define TOO_FEW_COLUMNS_FOR_MATRIX_LOAD 427	// "There are not enough columns in the file to load the matrix using the specified parameters."

#define LA_ONE_MD_WAVE_PER_BLOCK 428		// "WAVES declaration error: Each multi-dimensional wave must be in its own data block."
#define LA_BAD_DIMENSIONS 429				// "WAVES declaration error: Improper specification of wave dimension sizes."
#define FORMAT_STR_TOO_LONG 430				// "Format string too long."
#define FORMAT_STR_NO_PCT 431				// "Format string missing '%' character."
#define ILLEGAL_CONTOUR_FORMAT 432			// "Format string needs one '%f' or '%g', may have at most one '%d' and one '%s'."

#define MUST_BE_TWO_FREE_DIMENSIONS 433		// "There must be two free dimensions: one vertical (-2) and one horizontal (-3)."
#define ONE_OF_COLOR_COLORINDEX_OR_CTAB 434 // "Expected only one of color=, ctab= or cindex= keywords."
#define SUBDIVISIONS_REQUIRES_XYZ 435		// "The interpolate= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."
#define ZNULL_REQUIRES_XYZ 436				// "The nullValue= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."
#define ONLY_ONE_LEVELS_ARG_ALLOWED 437		// "Expected only one autoLevels= or manLevels= keyword."
#define EXPECT_1D_LEVEL_WAVE 438			// "Expected a 1D (single column) wave in \"manLevels=waveName\"."
#define EXPECT_CONTOUR_NAME 439				// "Expected the name of a contour matrix or z wave in the top graph."

#define BAD_OBJECT_TYPE 440					// "Illegal object type."
#define BAD_VAR_INDEX 441					// "Numeric variable index out of range."
#define BAD_STR_INDEX 442					// "String variable index out of range."
#define CONTOURXYZ_TOO_FEW_ROWS 443			// "Need at least 4 rows of X, Y, and Z to generate a contour"

#define INCOMPATIBLE_DIMENSIONING	444		// "Incompatible dimensions."
#define EXPECT_MATRIX_OR_VECTOR		445		// "Expected a matrix or vector and got a higher dimensioned object."
#define OSACompileErr				446		// "Got an error when compiling an OSA Script (AppleScript)."
#define OSAExecuteErr				447		// "Got an error when executing an OSA Script (AppleScript)."
#define EXPECT_STR_EXPR_NO_FUNC		448		// "Expected string expression NOT involving local variables or user functions."
#define CANT_FIND_SCRIPTING_SYSTEM	449		// "Can't connect to the Scripting System (AppleScript). Make sure it was installed correctly."
#define BAD_CHAR_IN_DF_NAME		450			// "A data folder name can't contain any of the following: ' \" ; : or any control chars."
#define AUTO_VAR_CONFLICT		451			// "Conflict creating V_ or S_ local variable. Probably user's NVAR, SVAR or wrong type."
#define FFT_ROWS_EVEN			452			// "rows must be even for FFT on real data"
#define BADSOUNDWAVEFREQ		453			// "Invalid sample rate for sound wave (use SetScale)"
#define OLD_SOUNDMGR			454			// "Sound manager version 3.0 or later is not present"
#define NO_TEXT_YTRACE			455			// "Can't use a text wave as y trace. Use ModifyGraph swapXY to get horizontal bars."
#define NVARorSVARfailed		456			// "Failed to resolve a local reference to a global variable or string (NVAR or SVAR)."
#define NOT_ON_WINDOWS			457			// "Feature not available on Windows."

#define CMD_ENDED_WITHOUT_NAME	458			// "Expected name"
#define CANT_CUT_MIX_OF_ROWS_AND_COLUMNS 459 // "Can't cut a mix of rows from one wave and columns from another."
#define ZNULL_MANUAL_OR_AUTO	460			// "Expected nullValueAuto or nullValue=value, but not both"
#define TOO_MANY_LEVELS			461			// "Too many levels (Max ^3)"

#define WAVE_DIMENSION_OR_VIEW_CONFLICT 462	// "All waves must have the same number of dimensions (at least two dimensions) and same viewed dimensions."
#define TRIANGULATION_REQUIRES_XYZ 463		// "The triangulation= keyword applies only to XYZ contour data, yet the named contour data is a Z matrix."

#define NOT_PACKED_EXP_FILE 464				// "This is not a packed Igor experiment file or it has been corrupted."
#define CONFLICT_DIFFERENT_TYPES 465		// "Can't create '^0' because that name is in use for a different type object."
#define INCOMPATIBLE_DATA_BROWSER 466		// "This version of the Data Browser is not compatible with this version of Igor."
#define SUB_DATA_FOLDER_NOT_FOUND 467		// "The specified sub data folder was not found in the file."
#define CONTOURXYZ_ZERO_AREA 468			// "X and Y values define a zero-area boundary: there is nothing to contour."
#define BAD_CONTOUR_LBL_FMT 469				// "Expected labelFormat=0, =1, =3, or =5."

#define XOP_EMPTY_CODE_FRAGMENT 470			// "The XOP's code fragment has nothing in it. It needs to be recompiled."
#define CANT_MODIFY_RGB_IMAGE 471			// "This is an RGB image; there is nothing to modify."
#define CANT_REPLACE_CONTOUR_TRACE 472		// "You can't replace a contour trace; it would just come back when the contours update."

#define OBJECT_DOES_NOT_EXIST 473			// "The named object does not exist."
#define WRONG_OBJECT_TYPE 474				// "The object is not of the specified type."

// End of Igor Pro 3.0 error codes.

#define FIRST_IGOR31_ERR 475				// Start of Igor Pro 3.1 error codes.

#define FUNCFIT_IND_VAR_MISMATCH	475		// "Number of independent variables in fit function does not match number of independent variable waves or columns."
#define FUNCFIT_DIMENSION_MISMATCH 476		// "Number of dimensions must match data wave."
#define FUNCFIT_ROW_MISMATCH 477			// "X waves must have the same number of rows as the data wave."
#define FUNCFIT_TOO_MANY_X_DIMS 478			// "Independent variable array must have no more than 2 dimensions."
#define FUNCFIT_ONE_X_MATRIX 479			// "When you use a matrix for independent variables, only one wave is allowed."
#define FUNCFITXYZ_XWAVE_REQUIRED 480		// "Your fitting function has more than one independent variable; you must specify at least one X wave using the /X flag."
#define FUNCFIT_TOO_MANY_IND_VARS 481		// "Too many independent variables; the limit is 16."
#define BAD_EWave 482						// "Epsilon wave does not match parameter wave (length or number type)"
#define FITFUNC_RETURNED_NAN 483			// "The fitting function returned NaN for at least one X value."
#define FITFUNC_RETURNED_INF 484			// "The fitting function returned INF for at least one X value."
#define MDFUNCFIT_TOO_MANY_IND_VARS 485		// "The fitting function has more than 4 independent variables."
#define MDFUNCFIT_TOO_FEW_IND_VARS 486		// "FitFuncMD requires at least 2 independent variables; use FitFunc instead."
#define MDFUNCFIT_IND_VAR_MISMATCH 487		// "The data wave must have one dimension for each independent variables."
#define BAD_RWave 488						// "Residual wave does not match parameter wave (length or number type)"
#define CONSTRAINT_HOLD_CONFLICT 489		// "Fitting parameter ^3 is being held and also constrained."
#define CONF_WRONG_WAVES	490				// "Wrong number of confidence band waves for the options selected."
#define CRVFIT_CONFLEVEL_NOT_IN_RANGE 491	// "Confidence level must be between 0 and 1, corresponding to 0 to 100% confidence levels."
#define CONSTRAINT_ILLEGAL 492				// "When parsing constraint expression:\n  \"^2\"\nreceived error message \"^3\""
#define CONSTRAINT_MULTIPLE_K 493			// "The constraint expression\n  \"^2\"\n has more than one fit parameter in a single term"  
#define CONSTRAINT_NONLINEAR 494			// "The constraint expression\n  \"^2\"\n lacks a conditional operator (< or >)."
#define CONSTRAINT_K_OUT_OF_RANGE 495		// "Fit parameter ^0 should be in range K0 to ^1 in constraint expression \"^2\"."
#define CONSTRAINT_NO_CONDITIONAL 496		// "Fit parameter ^1 is out of range in constraint expression\n  \"^2\"."
#define CONSTRAINT_ILLEGAL_OP_BEFORE 497	// "Illegal operator \'^0\' before fit parameter ^1 in constraint expression \"^2\"."
#define CONSTRAINT_ILLEGAL_OP_AFTER 498		// "Illegal operator \'^0\' after fit parameter ^1 in constraint expression \"^2\"."
#define MD_MISSING_XWAVE 499				// "The X wave was null or missing"
#define MD_MISSING_YWAVE 500				// "The Y wave was null or missing"
#define MD_MISSING_ZWAVE 501				// "The Z wave was null or missing"
#define MD_MISSING_TWAVE 502				// "The T wave was null or missing"
#define CURVEFIT_MISSING_XWAVE 503			// "X wave was null or missing"
#define CURVEFIT_MISSING_PWAVE 504			// "The parameter wave was null or missing"
#define CURVEFIT_MISSING_EWAVE 505			// "The Epsilon wave was null or missing"
#define CURVEFIT_MISSING_WWAVE 506			// "The weighting wave was null or missing"
#define CURVEFIT_MISSING_DWAVE 507			// "The destination wave was null or missing"
#define CURVEFIT_MISSING_RWAVE 508			// "The residual wave was null or missing"
#define CURVEFIT_MISSING_CMATRIX 509		// "The constraint matrix wave was null or missing"
#define CURVEFIT_MISSING_CVECTOR 510		// "The constraint vector wave was null or missing"
#define CURVEFIT_MISSING_CONFWAVE 511		// "A confidence band wave was null or missing"
#define CONSTRAINT_REQUIRES_TEXT_WAVE 512	// "Expected text wave containing constraint expressions"
#define CONSTRAINT_TEXT_WAVE_MISSING 513	// "Text wave for constraints was null or missing"

#define BAD_MISC_RECORD_VERSION 514			// "A miscellaneous data record is too new or corrupted."
#define CANT_CREATE_HOME_PATH 515			// "An error occurred while creating the home path."

#define BAD_PAGE_ORIENTATION 516			// "Expected 'Portrait' or 'Landscape'."
#define BAD_PRINT_PREFS	517					// "Page-Setup "

#define CANT_REMOVE_NORMAL_RULER 518		// "The Normal ruler can not be removed."

#define CANT_HANDLE_NON_NATIVE_PICT 519		// "This operation can't handle a <other platform> picture."
#define ERR_520 520							// Error 520 is available for duty.

#define SIN_FIFO_EXPECT_1CHAN 521			// "FIFO not setup properly (no channel info)"
#define SIN_FIFO_BAD_NUM_TYPE 522			// "FIFO number type not valid for sound input (8 and 16 bit ints only)"
#define NO_SIN_AVAIL 523					// "Sound input not available."
#define SIN_BAD_WAV_TYPE 524				// "Wave number type not valid for sound input (8 and 16 bit ints only)"
#define SIN_BAD_WAV_DIMS 525				// "Wave used for sound input can't have more than 2 dimensions."
#define SIN_FIFO_ALREADY_STARTED 526		// "Sound input to FIFO already started."
#define SIN_FIFO_ALREADY_STOPPED 527		// "Sound input to FIFO already stopped."
#define SIN_FIFO_ERROR  528					// "An error occurred during sound input to FIFO (probably overflow)."
#define TOOMANYCHANS 529					// "Sound input does not support the number of channels specified."
#define TOOMANYBITS 530						// "Sound input does not support the number of bits specified."
#define FREQNOTAVAIL 531					// "Sound input does not support the sampling rate specified."
#define SI_ERROR_SET_GAIN 532				// "Can't set the sound input gain."
#define SI_ERROR_SET_AGC 533				// "Can't set the sound input AGC."
#define AUDIO_FORMAT_NOT_AVAIL 534			// "The sound device does not support the specified sample rate, sample bits and/or channels."
#define AUDIO_BAD_VIBS 535					// "An audio related error occurred."

#define FAIL_READING_ENHMF 536				// "An error occurred while attempting to read an enhanced metafile."
#define FAIL_READING_AldusMF 537			// "An error occurred while attempting to read a placeable metafile."
#define FAIL_READING_WindowsMF 538			// "An error occurred while attempting to read a windows metafile."
#define FAIL_READING_DIB 539				// "An error occurred while attempting to read a device independent bitmap."

#define CANT_READ_THAT_GRAPHICS_TYPE 540	// "Can't read the graphics format of the specifed file."
#define NO_AUDIO_DEVICE 541					// "Could not find an audio device."
#define AUDIO_SYNCH_ONLY 542				// "Your audio setup does not support asynchronous output."
#define CLIP_NOT_AVAIL 543					// "Clipboard in use by another application."
#define PNG_WRITE_ERROR 544					// "An error occured while writing a PNG file."
#define CLIP_ERROR 545						// "A clipboard error occured."
#define EXPECT_SVAR 546						// "Expected name of string variable reference (SVAR)."
#define EXPECT_NVAR 547						// "Expected name of numeric variable reference (NVAR)."
#define EXPECT_0_OR_LPAREN 548				// "expected 0 or '('" */
#define EXPECT_KEYWORD 549					// "Expected keyword."

#define BAD_COLORTABLE_INDEX 550			// "IGOR color table index out of range"
#define BAD_COLORTABLE_HANDLE 551			// "Invalid IGOR color table handle"
#define BAD_COLORTABLE_PARAM 552			// "IGOR color table parameter out of range"

#define INCLUDE_FILE_ALREADY_OPEN_AS_NOTEBOOK 553	// "Included file is already open as a notebook"
#define INCLUDE_FILE_ALREADY_OPEN_AS_ANOTHER 555	// "Included file is already open as another type of window (e.g., help window)"

#define NOT_EXPERIMENT_FILE 555				// "This does not appear to be a valid IGOR experiment file"

#define FUNCLIST_UNKNOWN_OPTION 556			// "FunctionList does not recognize the option \"^0\""
#define CRVFIT_CONF_NO_MV 557				// "Confidence bands are not supported for multi-variate curve fits."
#define FIFO_ERR_SWAP_INFO 558				// "Probable bug while swapping FIFO info; please contact WaveMetrics tech support"
#define FIFO_ERR_SWAP_CHANINFO 559			// "Probable bug while swapping FIFO channel info; please contact WaveMetrics tech support"

#define ODE_STEP_SIZE_TOO_SMALL 560			// "The integration step size has gotten too small. You might try relaxing the error tolerance."
#define ODE_BAD_METHOD 561					// "The method keyword is set to a value that does not correspond to a supported method."
#define ODE_DERIVE_FUNC_MISSING 562			// "You must specify a function to calculate the derivatives."
#define ODE_PARAM_WAVE_MISSING 563			// "The pWave keyword is missing, or the wave specified by pWave keyword does not exist."
#define ODE_RESWAVE_MISSING 564				// "The resWave keyword is missing, or a wave or waves specified by the resWave keyword does not exist."
#define ODE_XWAVE_WRONG_NPNTS 565			// "The X wave (specified by the xvalues keyword) must have the same number of rows as the result waves (specified by the resWave keyword)."
#define ODE_SCALEWAVE_WRONG_NPNTS 566		// "The error scaling wave (errscale keyword) must have the same number of points as the parameter wave (pwave keyword)."
#define ODE_ZERO_SCALE 567					// "One or more of the values in the error scaling wave (errscale keyword) is zero."
#define ODE_MISSING_SCALE_WAVE 568			// "You have selected an error checking method (errorMethod keyword) that requires a wave to specify error scaling (errscale keyword)."
#define ODE_MISMATCHED_YWAVES 569			// "The lengths of the result waves (reswave keyword) must all be the same."
#define ODE_BAD_FUNCTION 570				// "The function ^0 has the wrong parameter types or the wrong number of parameters to use as an ODE function."
#define CURVEFIT_MISSING_YWAVE 571			// "Y wave was null or missing"
#define BAD_NOTEBOOK_PICTURE_MODE 572		// "Invalid notebook picture mode."
#define kHistEqHistWaveMustBe1D 573			// "The histogram wave must be a 1D wave"
#define kRequiresImageData 574				// "This operation works on image data (i.e., integer) waves only"
#define kImageHistSourceBadWaveType 575		// "Image histogram source wave must be 2D or 3D. If 3D, it must have exactly 3 layers."
#define kNumPointsMustBeEven	576			// "Histogram levels must be in pairs."
#define	kNoSuchThresholdMethod	577			// "No such threshold method."
#define kWavesMustBeSameSize	578			// "Both images must be the same size."
#define kWantsUnsignedCharData	579			// "This image operation supports only unsigned char (/B/U) data"
#define kMissingClipValue		580			// "Adaptive Histogram requires clip value (/C flag)"
#define kBadRegionSpecifier		581			// "Bad region specifier. Check that the image can be evenly divided into the number of specified regions."
#define kBadClipValue			582			// "Clipping value must be positive. /C=1 returns the original image."
#define kBadNumberOfBins		583			// "Bad value for the number of bins."
#define kNumHRegions2Big		584			// "The number of horizontal regions is too big."
#define kNumVRegions2Big		585			// "The number of vertical regions is too big."
#define kRequiresEqualRegions	586			// "Image size must be a multiple of the number of regions."
#define kRequiresMin4Regions	587			// "Adaptive histogram requires at least 4 sub regions."
#define kBadWaveInitialization	588			// "Bad wave initialization"
#define kBadWaveForMask			589			// "Bad wave for mask"
#define kNoSuchFilter			590			// "No such filter"
#define kNoSuchStructElement	591			// "Structure element is undefined."
#define kNeeds3DWave			592			// "This operation requires a 3D wave."	
#define kUnspecifiedThreshold	593			// "Threshold has not been specified"
#define kMethod3NotSupported	594			// "Method 3 is not supported for this operation."
#define kMustHaveROIWave		595			// "An ROI wave is required for this operation."
#define BAD_RECENT_FILES_PREFS	596			// "Recent Files "
#define INCLUDE_EXTENSION_NOT_ALLOWED 597	// "You must omit the \".ipf\" extension in #include statements."
#define PRINTER_DRIVER_SCREWUP 598			// "The printer driver returned an unexpected error. Check the driver and the printer."

#define FIRST_IGOR312_ERR 599				// Start of Igor Pro 3.12 error codes.
#define XOP_LOAD_LIBRARY_FAILED 599			// "Can not load XOP. It may be incorrectly constructed by the development system or may be corrupted."
#define CANT_FIND_XOP_MAIN 600				// "Can't find XOP's main function. (It must be declared HOST_IMPORT.)"
#define kBadROIDimensions		601			// "ROI dimensions do not match the target image."
#define BAD_FILE_TYPE_EXT_SPEC	602			// "A filetype/extension specifier must have exactly 4 characters or start with '.'"

#define FIRST_IGOR313_ERR 603				// Start of Igor Pro 3.13 error codes.

// These were created as platform-independent error codes that can be used when standard C file I/O routines return errors.
#define FILE_READ_ERROR 603					// "A file read error occurred."
#define FILE_WRITE_ERROR 604				// "A file write error occurred."
#define FILE_POS_ERROR 605					// "A file position error occurred."
#define FILE_OPEN_ERROR 606					// "The file can not be opened."
#define FILE_CLOSE_ERROR 607				// "An error occurred while closing the file."
#define FILE_PERM_ERROR 608					// "The file can not be opened because it is already open."
#define FILE_CREATE_ERROR 609				// "Error creating file. The file name or path may be invalid or the file may already exist."
#define FILE_EOF_ERROR 610					// "While reading a file, the end-of-file was reached."

#define XOP_LINK_FAILED 611					// "XOP dynamic link failed. The XOP may require a different version of Igor or may require additional DLLs."
#define SOME_LEVELS_MISSING 612				// "Some level crossings were not found."

// Root finder errors
#define ROOTS_TOO_MANY_EVALS 613			// "The root finder has performed more than the maximum allowed number of iterations."
#define ROOTS_NO_PROGRESS 614				// "The root finder not making progress. It may have gotten trapped at a point that is not a root, or your system may not have any roots."
#define ROOTS_NO_BRACKET 615				// "The root finder was unable to bracket a root before starting; for 1D roots it must find two X values where your function has opposite signs."
#define ROOTS_MISSING_X_WAVE 616			// "The X wave was not found."
#define ROOTS_MISSING_POLY_WAVE 617			// "The wave containing polynomial coefficients was not found."
#define ROOTS_MISSING_PWAVE 618				// "The parameter wave associated with the function ^0 was not found."
#define ROOTS_FUNCTION_TOO_MANY_PARAMS 619	// "Your function parameter wave had just one column, but your function has more than one independent variable. Your parameter wave must have a column for each independent variable."
#define ROOTS_FUNCTION_PARAMS_MISMATCH 620	// "The number of columns in the parameter wave must match the number of independent variables in the function."
#define ROOTS_WRONG_FUNCTION_PARAMS 621		// "You list ^0 functions; the functions must have ^1 parameters: a parameter wave and ^0 independent variables."
#define ROOTS_X_MISMATCH 622				// "Number of X values (length of X wave) does not match the number equations in your system."
#define ROOTS_FUNCTION_BAD_PARAMS 623		// "The function ^0 has the wrong parameter types or the wrong number of parameters."

#define CREATE_PROCESS_ERROR 624			// "CreateProcess failed with error \"^2\"."
#define CONTOUR_NODATA 625					// "Contour data is entirely blank (NaN). There is nothing to contour."

// End of Igor Pro 3.1 error codes.

#define IGORMENUMODE_BADMENUNAME 626		// "The menu \"^1\" is not recognized by SetIgorMenuMode."
#define IGORMENUMODE_BADITEMTEXT 627		// "The menu item text \"^1\" is not recognized by SetIgorMenuMode."
#define IGORMENUMODE_CANTDOTHATMENU 628		// "You can't enable/disable the items in that menu. Try disabling the item it is attached to in the parent menu."
#define IGORMENUMODE_BADACTION 629			// "The action to perform must be one of EnableItem, DisableItem, EnableAllItems or DisableAllItems."

#define CURVEFIT_BAD_POLY2D_ORDER 630		// "Poly 2D order must be at least 1."
#define CURVEFIT_MV_MISSING_XWAVES 631		// "You are fitting a multivariate fit function to a one-column data set.\nYou must select an X wave for each independent variable."
#define CURVEFIT_MISSING_CONF_LEVEL 632		// "You have checked the Error Analysis checkbox. You must enter a confidence level."
#define CURVEFIT_NO_Y_WAVE 633				// "You have not selected data to fit in the Y Data menu."
#define CURVEFIT_AMBIGUOUS_COEFS 634		// "Igor can't tell how many fit coefficients are required by the fit function. You must select a wave from the Coefficient Wave menu",
#define CURVEFIT_NEWWAVE_CONFLICT 635		// "You have selected 'New' from the %s menu, but a wave with the name %s already exists."
#define CURVEFIT_MISSING_NEWWAVE 636		// "You have selected 'New' from the %s menu, but you have not entered a name in the New Wave box."
#define CURVEFIT_MISSING_GUESS_USER 637		// "You have selected a user-defined fit function so you must enter an initial guess for every fit coefficient."
#define CURVEFIT_MISSING_GUESS_BUILTIN 638	// "You have selected manual guesses so you must enter an initial guess for every fit coefficient."
#define CURVEFIT_MISSING_GUESS_HOLD 639		// "You have checked a box to hold the value of a fit coefficient, but you have not entered a value for that coefficient. Go to the Coefficients Tab to do this."
#define CURVEFIT_MISSING_EPSILON_VALUE 640	// "You have selected an Epsilon wave so you must enter values in the Epsilon column of the Coefficients list."
#define CF_PLOTIT_IMSORRY 641				// "I'm sorry- it was impossible to plot your fitting function because"
#define CF_PLOTIT_MultiVariate 642			// "this service is not offered for multivariate fitting functions."

#define EXPECTED_GRAPH 643					// "expected name of a graph"

#define CF_PLOTIT_NoYWave 644				// "you have not selected a data wave from the 'Fit To' menu on the Input Data tab."
#define CF_PLOTIT_NoGraph 645				// "you need to make a graph."
#define CF_PLOTIT_WaveNotOnGraph 646		// "the fit data is not on the top graph."
#define CF_PLOTIT_YWaveIsXWave 647			// "on the top graph the fit data is used as an X wave."
#define CF_PLOTIT_ErrorMakingWave 648		// "there was an error making the destination wave. Igor is probably out of memory."
#define CF_PLOTIT_CouldntAddDisplay 649 	// "there was an error adding the destination wave to the graph. Igor is probably out of memory."
#define CF_PLOTIT_EquationNotAccepted 650	// "there was an error involving creating the function expression. Igor is probably out of memory."
#define CF_PLOTIT_NoCoefs 651				// "Igor couldn't determine the number if fit coefficients.\nTry selecting an appropriate coefficients wave from the Coefficients Wave menu."
#define CF_PLOTIT_ErrorMakingCWave 652 		// "Igor was unable to create a coefficients wave for the operation. Igor is probably out of memory."
#define CF_PLOTIT_NoGuess 653				// "you must enter initial guesses for every coefficient."
#define CURVEFIT_BADNCOEF 654				// "The fitting function requires a coefficient wave with ^0 points."
#define BAD_MWave 655						// "You have provided a mask wave that does not match your data wave."
#define CURVEFIT_MISSING_MWAVE 656			// "Mask wave was null or missing"

#define kMustHaveRowsAndCols	657			// "You must specify /N={extraRows,extraCols} for this operation" 	19MAR99		4.00D01
#define kBadCMAPWaveType		658			// "This wave is not appropriate as a CTAB."						19MAR99		4.00D01
#define kRowAndColsMustBeReasonable	659		// "The requested padding is not appropriate for this wave."		19MAR99		4.00D01
#define kBadWidthOrHeight			660		// "Improper width or height specification."
#define kBadXYWaves					661		// "Mismatch or inappropriate x, y waves."
#define kBadSeedPixel				662		// "Bad seed pixel specification."									09AUG99
#define kBadIPParamSpec				663		// "Bad IP Parameter specification."								09AUG99
#define kExpectedImageTypeName		664		// "Expected name of image file type."
#define kBadSrcWave					665		// "Source wave is bad."								AG	07JUN01
#define kXYZWaveMustBeSP			666		// "XYZ waves must be of type SP (NT_FP32)"				AG	28JUL99
#define kBadMultipleImageCount		667		// "Bad count for multiple images."						AG 	20SEP99
#define kBadValueForFirstImage		668		// "Bad value for first image."							AG	20SEP99
#define kWantsNewerQuickTime		669		// "Operation requires a newer version of QuickTime."	AG 	20SEP99

#define kMustSpecifyDataWave		671		// "Data wave must be specified (see /D flag)."			AG 	11OCT99 
#define kIncompatibleFlagOptions	672		// "Flag options appear to be incompatible."			AG 	29OCT99
#define kAGBloatedFortranCrap		673		// "AG's bloated fortran crap is not compiled."			AG	01NOV99
#define kExpectSquareMatrix			674		// "Expected square matrix."							AG 	05NOV99
#define kInsufficientPointsInWave	675		// "Insufficient number of points in the input wave."	AG	12NOV99
#define kAllPointsColinear			676		// "All points are co-linear."							AG 23NOV99
#define kAllPointsCoPlanar			677		// "All points are co-planar."							AG 23NOV99
#define kFailedConsistency			678		// "Failed consistency test."							AG 23NOV99
#define kNotConvex					679		// "Data does not represent convex set."				AG 23NOV99
#define kFailedEulerTest			680		// "Data failed Euler test."							AG 23NOV99
#define kExpectedTripletWave		681		// "Expected 3 column (Triplet) wave."					AG 03MAY00
#define kExpectedNT_FP32			682		// "Expected an SP (single precision float) wave."		AG 24MAY00
#define kExpectedSPorDP				683		// "Expected SP or DP wave"
#define kInputTooLarge				684		// "Input is too large"
#define kExpectNoNaNorInf			685		// "Source wave should not contain NAN's or INF's"		AG 18MAY01
#define kFailedInternalConsistencyTest	686	// "Failed internal consistency text"

#define ButtonRecursionDetected	687		// "Recursion prevented. The expression \"^0\" causes itself to be executed again, and needs to be changed."
#define NO_SUCH_TOOL_NAME	688			// "Expected \"normal\", \"arrow\", \"text\", \"line\", \"rect\", \"rrect\", \"oval\", or \"poly\"."

#define NO_MOVIE 689					// "no movie"
#define FAILED_TO_PLAY_MOVIE 690		// "failed to play movie"
#define EXPECT_WATERFALL 691			// "expected a waterfall plot"
#define X_VECTOR_MISMATCH 692			// "x vector mismatch"
#define Y_VECTOR_MISMATCH 693			// "y vector mismatch"
#define EXPECTED_INSTANCE 694			// "expected instance"
#define UNMATCHED_CONDITIONAL 695		// "unmatched ?: conditional"
#define NOT_IN_MACROS 696				// "this syntax is not allowed in macros -- use functions"
#define LINK_NO_XFUNC 697				// "During link, couldn't find external function."
#define LINK_TYPE_XFUNCMISMATCH 698		// "During link, external function did not match."
#define LINK_NO_FUNC 699				// "During link, couldn't find user function."
#define LINK_TYPE_MISMATCH 700			// "During link, user function did not match."
#define LINK_NO_CONST 701				// "During link, couldn't find constant."
#define SEMICOLON_EXPECTED 702			// "Expected semicolon."
#define EXPECT_OBJ_NAME 703				// "Expected object name."
#define EXPECT_CONSTANT_OR_LITERAL 704	// "Expected symbolic constant or literal"
#define EXPECT_FUNC_NAME 705			// "Expected function name."
#define FUNCREF_TYPE_INCONSISTENT 706	// "Function reference type inconsistent."
#define CANT_USE_FUNCREF_HERE 707		// "Can't use a function reference here."
#define EXPECT_LOCALVAR_NAME 708		// "Expected a local variable name."
#define REF_VAR_DIFFERENT_TYPE 709		// "Reference variable is of a different type."
#define COULD_NOT_FIND_PROTO_FUNC 710	// "Couln't find prototype function."
#define EXPECT_FUNC_REF 711		// "Expected function reference."
#define NO_STATIC_FUNC_HERE 712	// "Can't use a static function here."
#define NO_PROMPT_THIS_TYPE 713			// "Can't prompt for this type of variable."
#define EXPECT_POPUP 714				// "Expected popup keyword."
#define NO_PROMPT_DEFINED 715			// "No prompt defined for this variable."
#define FASTOP_TOO_MANY_PRODS 716		// "Too many product terms."
#define FASTOP_SYNERR 717				// "Syntax does not conform to FastOp requirements."
#define WAVE_LENGTH_MISMATCH 718		// "Wave lenght mismatch."
#define WAVE_TYPE_MISMATCH 719			// "Wave type mismatch."
#define COMPLEX_INT_NOT_AVAIL 720		// "Complex integers not supported here."
#define COMPLEX_TO_REAL_LOSS 721		// "Complex wave used in real expression."
#define DUP_CONST 722					// "Duplicate constant."
#define WRONG_IGOR_VERS 723				// "A more recent version of Igor is required."
#define NOGRAF_OR_PANEL_OR_TABLE 724	// "Expected graph, panel or table."
#define EXPECT_WIN_NAME 725				// "Expected window name."

#define NOSUBRANGECATWAVE 726				// "Subranged category waves is not currently supported."
#define ONLY_ONE_RANGE 727					// "Only one dimension may have a range."
#define DUPLICATE_RESOURCES 728				// "This experiment can not be loaded in Carbon Igor because of duplicate resources. Open and save it in pre-Carbon Igor to fix the problem."
#define EXPECT_GUDE_NAME 729				// "Expected guide name."
#define NO_REDEFINE_BI_GUIDE 730			// "Can't redefine a built-in guide."
#define NO_SUCH_GUIDE 731					// "Specified guide does not exist."
#define NO_MIX_GUIDE_ORIENTATION 732		// "Guide orientation mismatch."
#define NO_SWITCH_GUIDE_ORIENTATION 733		// "Can't switch guide orientation."
#define ILLEGAL_DUAL_GUIDE_POS 734			// "Illegal dual guide position."
#define NO_CONTROLS_IN_SUBGRAPH 735			// "Can't put controls in a subgraph."
#define NO_PANEL_IN_SUBGRAPH 736			// "Can't put panels in a subgraph."
#define WRONG_GUIDE_ORIENTATION 737			// "Wrong guide orientation."
#define GUIDE_IN_GRAPH_ONLY 738				// "Guide is for graphs only."
#define NOT_VALID_GUIDE_NAME 739			// "Invalid guide name."
#define NO_SUCH_SUBWINDOW 740				// "Specified subwindow does not exist."
#define NOT_AVAILABLE_ON_SUBRANGE 741		// "Action not available on subrange."
#define COMPILE_ONLY 742					// "This feature can only be used in user functions."

#define MACROLIST_UNKNOWN_OPTION 743		// "MacroList does not recognize the option \"^0\""
#define EXPECTED_GRAF_OR_PANEL 744			// "Expected graph or panel name."

#define SEARCH_CANT_FIND_WINDOW	745			// "Can't find the referenced window in the experiment. Try doing the search again."
#define SEARCH_TOO_MANY_TERMS 746			// "Too many search terms - only 8 are allowed."
#define SEARCH_FILE_MODIFIED 747			// "The file was modified after the search. Try doing the search again."

#define WAVE_REF_FAILED	748					// "Failed to resolve a local WAVE \"^0\" reference to a global wave."
#define NVAR_REF_FAILED 749					// "Failed to resolve a local NVAR \"^0\" reference to a global variable."
#define SVAR_REF_FAILED 750					// "Failed to resolve a local SVAR \"^0\" reference to a global string."

#define DIM_LABEL_TOO_LONG 751				// "Dimension labels are limited to 31 characters."

#define CF_PLOTIT_NoGuessBuiltin 752		// "you must enter function coefficient values. Select Manual from the Guess Method menu first."
#define CF_PLOTIT_NoGuessPolyLine 753		// "you must enter function coefficient values. For poly and line fits, check the Hold box in order to enter values."

#define COLUMN_INFO_UNKNOWN_SPECIFIER 754	// "Expected 'C', 'F', 'W', or 'N' in column info specifier."
#define COLUMN_INFO_EXPECTED_NAME 755		// "Expected a wave name in column info specifier."
#define COLUMN_INFO_EXPECTED_NUMBER 756		// "Expected a number in column info specifier."
#define COLUMN_INFO_BAD_NAME 757			// "A name in the column info specifier contained illegal characters."
#define COLUMN_INFO_BAD_NAME_TERM 758		// "Missing comma or semicolon after a name in the column info specifier."
#define COLUMN_INFO_NAME_TOO_LONG 759		// "A name in the column info specifier can not exceed 31 characters in length."
#define COLUMN_INFO_BAD_NUMBER 760			// "Bad number in the column info specifier."
#define COLUMN_INFO_BAD_NUMBER_TERM 761		// "Missing comma or semicolon after a number in the column info specifier."
#define BAD_FIXED_FIELD_NUMBER_OF_COLUMNS 762	// "The number of columns in a fixed field file must be between 1 and 10000."
#define BAD_FIXED_FIELD_FIELD_WIDTH 763			// "The field width in a fixed field file must be >= 1."
#define EXPECT_IMAGE_CONTOUR_TRACE_NAME 764		// "Expected name of image, contour, or trace in top graph."
#define CL_LOOKUP_REQUIRES_CTAB 765				// "ColorScale lookup requires ctab keyword."
#define CURVEFIT_NOSTATICFUNCTIONS 766			// "Static function references are not allowed as curve fitting functions."
#define EXPECTED_GRAPH_TABLE_LAYOUT_PROCEDURE_NOTEBOOK_CMDWIN_BUT_NOT_PANEL 767	 // "This operation is for graphs, tables, layouts, notebooks, procedure windows, or the command/history window."
#define EXPECTED_TARG_PROC_CMDWIN_NAME 768		// "expected name of a target window, procedure window, or \"kwCmdHist\"."
#define CURVEFIT_NOTENOUGHPOINTS 769			// "You must have at least as many data points as fit parameters."
#define BAD_TIMEUNIT 770						// "Expected the name of a time unit like sec, week, year, etc."
#define MANDATE_INCMUSTBEINTEGER 771			// "The manual tick increment for a date/time axis must be an integer."
#define DATEAXIS_NOSWAPENDS 772					// "Date/time axes do not support reversed scaling."
#define FLAG_MUST_BE_FIRST_OR_AFTER_SLASH_W 773	// "this must be first flag in command, or immediately follow /W."
#define EXPECTED_TARG_PROC_NAME 774				// "expected name of a target window or procedure window."

#define EXPECTED_GRAPH_MARGIN_VAL 775			// "Expected a graph margin value."
#define EXPECTED_GRAPH_PLOT_AREA_VAL 776		// "Expected a graph plot area width or height value."

#define LAYOUT_CAN_NOT_BE_COMPILED 777			// "The Layout operation can not be used in a function. Use NewLayout instead."
#define APPENDTOLAYOUT_CAN_NOT_BE_COMPILED 778	// "The AppendToLayout operation can not be used in a function. Use AppendLayoutObject instead."
#define BAD_FRAME_VALUE 779						// "Expected a frame value: 0=none, 1=single, 2=double, 3=triple, 4=shadow."
#define BAD_LAYOUT_OBJECT_TYPE 780				// "Expected a page layout object type keyword: graph, table, picture or textbox."
#define EXPECTED_LAYOUT_NAME 781				// "Expected the name of a page layout window."
#define LAYOUT_USE_TEXTBOX_CMD 782				// "Can't append a textbox via AppendLayoutObject. Use Textbox or Legend instead."
#define FLAG_ALLOWED_JUST_ONCE 783				// "This flag can be used only once per command."

#define OPTIMIZE_NOXFLAG 784					// "When optimizing a multivariate function, you must provide a starting guess using the /X flag."
#define OPTIMIZE_NAN 785						// "The user function you are trying to optimize has returned NaN (Not a Number)."
#define OPTIMIZE_NOBRACKET 786					// "The Optimize operation was unable to find a pair of X values that bracket the minimum (or maximum). Use /L and /H to specify bracketing values."
#define OPTIMIZE_NOPROGRESS 787					// "The Optimize operation could not find a better solution than your starting guess. Your function may be too non-linear, the stopping tolerance may be too large, or your starting guess is a solution."
#define OPTIMIZE_TOOMANYITERATIONS 788			// "The Optimize operation has performed more ^0 iterations, which is more than the maximum allowed."
#define OPTIMIZE_MAXSTEPSIZE 789				// "The Optimize operation  has exceded the maximum step size. It may be that your function is unbounded or approaches a value assymptotically."
#define OPTIMIZE_NTYPSIZEMISMATCH 790			// "The number of values used with the /R flag must match the number of X values."
#define OPTIMIZE_CRITICALPOINT 791				// "Your starting guess is too near a critical point and Optimize can't procede. Try a different starting guess."

#define EXPECTED_LOCAL_NUM_OR_STR_VAR_NAME 792	// "Expected the name of a local numeric variable or NVAR, or a local string variable or SVAR."
#define SSCANF_RAN_OUT_OF_CONVERSIONS 793		// "The sscanf format string does not have enough conversion specifiers or there are too many output variables."
#define SSCANF_TOO_MANY_CONVERSIONS 794			// "The sscanf format string has too many conversion specifiers or there are too few output variables."
#define SSCANF_REQUIRES_L 795					// "\"%e\", \"%f\", and \"%g\" are not allowed in sscanf. Use "\"%le\", \"%lf\", and \"%lg\" instead."
#define SSCANF_L_NOT_ALLOWED 796				// "Do not put an 'l' after the '%' in an sscanf format string."
#define SSCANF_UNKNOWN_FORMAT 797				// "sscanf unsupported format."
#define SSCANF_IN_FUNCTIONS_ONLY 798			// "sscanf can be used in user functions only."
#define SSCANF_EXPECTED_NUM_GOT_STR 799			// "sscanf expected a numeric variable and got a string variable. Check the format string and the variable list."
#define SSCANF_EXPECTED_STR_GOT_NUM 800			// "sscanf expected a string variable and got a numeric variable. Check the format string and the variable list."
#define SSCANF_SCANSET_NOT_TERMINATED 801		// "A scanset (\"%[...]\") was used in sscanf but there was no trailing \"]\" character."
#define SSCANF_TOO_MANY_PARAMS 802				// "sscanf is limited to 100 localVar parameters."
#define CONSTRAINTS_NOT_AVAILABLE 803			// "This version of Igor does not support curve fitting with constraints."
#define EXPECTED_NUMERIC_WAVE 804				// "Expected numeric wave."
#define EXPECTED_TEXT_WAVE 805					// "Expected text wave."

#define CURVEFIT_BADAUTODESTLEN 806				// "The Destination length must be 2 or greater, or \"Auto\" or zero."
#define CF_PLOTIT_NoAllAtOnce 807				// "this service is not offered for all-at-once fitting functions."

// End of Igor Pro 4-compatible error codes.

#define FIRST_IGOR5_ERR 808

#define CANT_OPEN_SPECIFIED_PRINTER_DRIVER 808	// "Unable to open the printer driver for this page setup record."
#define INCOMPATIBLE_PRINT_RECORD 809			// "The page setup record is not compatible with the current printer driver."

#define GETPATH_CALLBACK_OBSOLETE 810			// "The GetPath XOP callback is no longer supported. The XOP needs to be updated for this version of Igor."
#define TU_CALLBACKS_OBSOLETE 811				// "The TUWriteFile and TUInsertFile XOP callbacks are no longer supported. The XOP needs to be updated for this version of Igor."

#define NOT_IMPLEMENTED 812						// "This feature is not yet implemented."
#define FCMD_BAD_INTERACTIVE 813				// "In /I=i, expected i between 0 and 3."
#define FCMD_CANT_COPY_SELF 814					// "Can't copy a file onto itself."
#define FCMD_CANT_MOVE_SELF 815					// "Can't move a file onto itself."

#define XOP_68K_NOT_SUPPORTED 816				// "This is a 68K XOP. It is not compatible with this version of Igor."

#define EXPECTED_MENU_ITEM_FLAG 817				// "Expected a menu item flag (/Q)."

#define XOP_CANT_RUN_ON_OS_X 818				// "The XOP '^3' can not run on Mac OS X."

#define NEGATIVE_FIFO_SIZE 819					// "The fifo size is negative, which is not allowed. Perhaps you used a null variable to set the size."
#define FIFO_SIZE_TOO_BIG 820					// "The fifo size is too big."

#define CANT_TARGET_ROOT_DIRECTORY 821			// "This operation does not permit targeting the root directory of a volume. Target a sub-directory."

#define kBadParameters	822						// "One or more parameters are inappropriate."
#define kParameterOutOfRange	823				// "Parameter is out of range"
#define kBadStructureElement	824				// "Bad Structure Element"
#define kExpectedComplexWave	825				// "Expected Complex Wave"
#define kDivideByZero			826				// "Operation failed beacue of a divide by zero condition"
#define kBadParam				827				// "Bad parameter."
#define kBadWaveletSpecification	828			// "Bad Wavelet Specification."
#define kBadOutputSpecification		829			// "Bad Output Specification."
#define kBadOffsetsSpecification	830			// "Bad Offset Specification."
#define kBadNumCoeff				831			// "Bad Number of Coefficients."
#define kBadRangeSpecification		832			// "Bad range specification."
#define kDestinationMustBeSpecified 833			// "Destination wave must be specified."
#define kRequirePositiveParameter   834			// "Positive number required."
#define kBadTriangulationWave		835			// "Bad Triangulation Wave."
#define kBadDestinationWave			836			// "Bad Destination Wave."
#define kDoesNotSupport4D			837			// "Does not support 4D waves."
#define kBadUserFunctionFormat		838			// "Bad user function format."
#define FFT_COLS_EVEN				839			// "The number of columns must be even."
#define kWave_Scaling_Mismatch		840			// "Wave Scaling Mismatch"
#define kBadMatrixOPToken			841			// "Bad MatrixOPs token."
#define kMatrixDimMismatch			842			// "Matrix dimensions mismatch."
#define kNeedSquareMatrix			843			// "Expected square matrix."
#define kRequireLeftHand			844			// "Left Hand Side is required."
#define kLeftSideMustBeWaveOrVar	845			// "Left Hand Side must be a wave."
#define kBadTokenIndex				846			// "Bad matrix token index."
#define kMatrixStackOutOfSync		847			// "Matrix stack is out of sync."
#define kUnbalancedParenthesis		848			// "Unbalanced parenthesis"
#define kBadProcessRange			849			// "Bad process range."
#define kExpectedOperator			850			// "Expected operator."
#define kCouldNotFindData			851			// "Could not find data."
#define kNaNsNotAllowed				852			// "NaNs are not allowed in this operation."
#define kExpectedPrefixOperator		853			// "Expected prefix operator."
#define kUnknownPrefixOperator		854			// "Unknown prefix operator."
#define kBadDataType				855			// "Bad data type."
#define kExpectRealMatrix			856			// "Expected real matrix."
#define kCantSubtractMatrixFromScalar 	857		// "Can't subtract matrix from a scalar."
#define kBadROIWaveType					858		// "Bad ROI wave type."
#define kBadNumberOfClasses				859		// "Bad number of classes."
#define kLAPACKError					860		// "LAPACK returned an error.  Check your input matrix."
#define kNoComplexSupport				861		// "Does not support complex waves or variables."
#define kAPMMemoryError					862		// "Arbitrary Precision Math memory error."
#define kDoesNotSupportNaNorINF			863		// "Does not support NaN or INF."
#define AGerror864								// "Reserved"
#define AGerror865								// "Reserved"
#define AGerror866								// "Reserved"
#define AGerror867								// "Reserved"
#define AGerror868								// "Reserved"
#define AGerror869								// "Reserved"
#define kExpectedNT_FP64				870		// "Expected an DP (double precision float) wave."		

#define SA_NEED_3_COLUMN_WAVE	871				// "A 3-column wave is required."
#define SA_NEED_2_COLUMN_WAVE	872				// "A 2-column wave is required."
#define SA_MISSING_STEP_SIZE_WAVE 873			// "The simulated annealing step size wave is null or missing."
#define SA_MISSING_MINMAXXWAVE 874				// "The simulated annealing X limit wave is null or missing."
#define PROGRAMMED_STOP 875						// "Termination requested by function."
#define MISMATCHED_NUMBER_OF_POINTS 876			// "Different number of points for X and Y waves." used by New Graph dialog, DisplayWaves.cpp
#define EXPECT_1D_WAVE_FROM_LIST 877			// "2D waves must have a single row or column selected. Click the Add button and select the row or column in the list below."

#define OH_UNKNOWN_PARAM_TYPE 878				// "Valid parameter types are number, string, name, wave."
#define OH_EXPECTED_POSTFIX_CHAR 879			// "Expected trailing ), ] or }."
#define OH_TOO_MANY_PARAMS 880					// "The template contains too many parameters."
#define OH_RUNTIME_STRUCT_TOO_SMALL 881			// "The template requires a runtime structure larger than the specified maximum."

#define OH_TOO_MANY_FLAG_PARAM_GROUPS 882		// "The template contains too many flag parameter groups."
#define OH_TOO_MANY_MAIN_PARAM_GROUPS 883		// "The template contains too many main parameter groups."
#define OH_OPTIONAL_PARAMS_MUST_BE_AT_END 884	// "Optional flag or keyword parameters must appear at the end of a group introduced by (, [ or {."
#define OH_OPTIONAL_SIMPLE_MAIN_PARAMS_MUST_BE_AT_END 885	// "Optional simple main parameters must be the last parameter group in the template."
#define OH_BAD_NUMBER_OF_OPTIONAL_PARAMS 886	// "Number of optional parameters must be between 1 and 100."
#define OH_EXPECTED_EQUALS 887					// "Expected '='."
#define OH_EXPECTED_COMMA 888					// "Expected comma between parameters."
#define OH_EXPECTED_LEFT_PAREN 889				// "Expected '('."
#define OH_EXPECTED_RIGHT_PAREN 890				// "Expected ')'."
#define OH_EXPECTED_LEFT_BRACKET 891			// "Expected ']'."
#define OH_EXPECTED_RIGHT_BRACKET 892			// "Expected '['."
#define OH_EXPECTED_LEFT_BRACE 893				// "Expected '{'."
#define OH_EXPECTED_RIGHT_BRACE 894				// "Expected '}'."
#define OH_CANT_MIX_SIMPLE_AND_KEYWORD_PARAMS 895	// "You can't mix simple and keyword main parameters in the same operation."
#define OH_BAD_XOP_OPERATION_NAME 896			// "The XOP does not add an operation with this name."
#define OH_OPERATION_NOT_REGISTERED 897			// "The operation is not registered with Igor."
#define OH_COMPILABLE_BIT_MUST_BE_SET 898		// "The operation's compilable bit must be set in the XOPC resource."
#define OH_BAD_RUNTIME_PARAM_STRUCT_SIZE 899	// "The runtime parameter structure size is too small."
#define OH_OUT_OF_MEMORY 900					// "Operation Handler ran out of memory while parsing command template."
#define OH_SYNTAX_REQUIRES_PREFIX_CHAR 901		// "This type of optional parameter syntax works only with {}, [] or () syntax."
#define OH_TOO_MANY_LEFT_BRACKETS 902			// "Use only one pair of brackets to indicate optional parameters within {}, [] or () syntax."
#define OH_EXPECTED_NUMBER 903					// "Expected a number or numeric expression"
#define OH_EXPECTED_STRING 904					// "Expected a string or string expression"
#define OH_EXPECTED_NAME 905					// "Expected a name"
#define OH_EXPECTED_VARNAME 906					// "Expected name of a variable, NVAR or SVAR"
#define OH_EXPECTED_WAVE 907					// "Expected a wave"
#define OH_EXPECTED_WAVE_TYPE 908				// "Expected wave type (real, complex, text)"

#define NO_SUCH_OPERATION 909					// "There is no operation with this name."
#define BAD_OPERATION_VAR_NAME 910				// "Names of variables created by operations must start with V_ (numeric) or S_ (string)."

#define EXPECTED_NUM_WAVE 911					// "Expected a numeric wave."
#define EXPECTED_NUM_VAR_OR_NVAR 912			// "Expected the name of a numeric variable or an NVAR."

// These are used for XOP's calling user functions.
#define BAD_COMPILATION_INDEX 913				// "A call to an internal Igor routine (CallUserFunction) was made using stale information."
#define REQUIRES_NUMERIC_PARAMETER 914
#define REQUIRES_PASS_BY_REFERENCE_NUMERIC_PARAMETER 915
#define REQUIRES_COMPLEX_NUMERIC_PARAMETER 916
#define REQUIRES_PASS_BY_REFERENCE_COMPLEX_NUMERIC_PARAMETER 917
#define REQUIRES_STRING_PARAMETER 918
#define REQUIRES_PASS_BY_REFERENCE_STRING_PARAMETER 919
#define REQUIRES_NUMERIC_WAVE_PARAMETER 920
#define REQUIRES_COMPLEX_NUMERIC_WAVE_PARAMETER 921
#define REQUIRES_TEXT_WAVE_PARAMETER 922
#define REQUIRES_FUNCTION_REFERENCE_PARAMETER 923
#define UNKNOWN_PARAMETER_TYPE 924
#define REQUIRES_NUMERIC_RETURN_TYPE 925
#define REQUIRES_COMPLEX_NUMERIC_RETURN_TYPE 926
#define REQUIRES_STRING_RETURN_TYPE 927
#define UNKNOWN_RETURN_TYPE 928
#define FUNCTION_HAS_TOO_FEW_PARAMETERS 929
#define FUNCTION_HAS_TOO_MANY_PARAMETERS 930
#define PASS_BY_REF_BAD_PARAM_TYPE 931			// "Only numeric and string parameters can be pass-by-reference"

#define XOP_MACH_CANT_FIND_EXECUTABLE 932		// "Can't find the executable file in the XOP package."
#define XOP_MACH_CANT_LOAD_EXECUTABLE 933		// "Can't load the executable file in the XOP package."
#define XOP_MACH_CANT_FIND_MAIN 934				// "Can't find the main function in the XOP package."
#define XOP_MACH_CANT_RUN_ON_OS9 935			// "Mach XOPs run on OS X only, they can not execute on Mac OS 9."

#define STOP_RECURSING_THROUGH_FOLDERS 936		// This is not really an error. It is a signal to RecurseThroughFolders.

#define EXPECTED_WSIZE 937

#define WARN_OVERRIDE_NOT_IN_MAIN_PROC 938		// "Override functions should be in the main Procedure window to work reliably."
#define WARN_FUNCTION_IS_OVERRIDDEN 939			// "This function will be overridden by a pre-existing Override function."
#define STATIC_NEEDS_MODULE_NAME 940			// "This procedure window lacks the #pragma Module statement static functions require."

#define EXPECT_STRUCT_NAME 941					// "Expected structure name."
#define EXPECT_STRUCT 942						// "Expected structure."
#define EXPECT_COMPAT_STRUCT 943				// "Expected compatible structure."
#define EXPECT_WAVE_FIELD 944					// "Expected WAVE field."
#define EXPECT_NVAR_FIELD 945					// "Expected NVAR field."
#define EXPECT_SVAR_FIELD 946					// "Expected SVAR field."
#define EXPECT_FUNCREF_FIELD 947				// "Expected FUNCREF field."
#define EXPECT_LOCALCONSTANT_OR_LITERAL_EXPR 948	// "Expected locally defined constant or literal expression."
#define NOT_OPTPARAM_NAME 949					// "Expected name of optional parameter."
#define NO_OPTPARAM_FUNCS 950					// "Functions with optional parameters not allowed here."
#define DUP_PROCPICT 951						// "Duplicate procedure Picture."
#define EXPECTED_END 952						// "Expected End."
#define EXPECT_ASCII85Begin 953					// "Expected ASCII85Begin."
#define DUP_PROCSTRUCT 954						// "Duplicate procedure Structure."
#define STRUCT_FIELD_ARRAY_OUTOFBOUNDS 955		// "Illegal structure field array size."
#define NO_OVERRIDE 956							// "Can't use Override here."
#define NO_SUCH_STRUCT 957						// "No such structure exists."
#define STRUCTS_TOO_DEEP 958					// "Structure nesting too deep."
#define STRUCT_FIELD_INDEX_OB 959				// "Structure field index out of bounds."
#define EXPECT_STRUCT_FIELD 960					// "Expected structure field."
#define ILLEGAL_FIELD_FOR_FBIN 961				// "Illegal field in structure (String, NVAR etc.) for this use."
#define DUP_FIELD_NAME 962						// "Duplicate field name."
#define CANT_CHANGE_LOCKED_WAVE 963				// "Can't change a locked wave."
#define AXIS_NAME_USED 964						// "An axis of that name already exists."
#define BAD_MASTER_AXIS 965						// "Specified master axis not found."
#define BAD_AXIS_HOOK 966						// "Axis hook function not found."
#define WRONG_FUNC_PARAMS 967					// "Function input parameters or output not valid for this use."
#define AXIS_IN_USE 968							// "Axis is in use."
#define STRUCT_REF_ONLY 969						// "Structure input parameters must be pass-by-reference ('&' needed.)"
#define USING_NULL_REFVAR 970					// "attempt to use uninitialized pass-by-reference variable"
#define FONT_ERR 971							// "General font error."
#define FONT_CURVE_EXTRACT_ERR 972				// "Error extracting font outline curves."
#define INCOMPATIBLE_STRUCT_VERSION 973			// "Incompatible structure. The calling function is too old or too new for the called function."
#define STRUCT_ONLY_IN_FUNCTION 974				// "Structure parameters can be used only in user-defined functions"
#define BAD_FUNCREF 975							// "FUNCREF does not reference valid function."
#define BAD_WAVE_LIST 976						// "There was a problem in a list of waves."
#define REQUIRES_STRUCTURE_PARAMETER 977		// "Requires structure parameter". Used for XOP's calling user functions.
#define LH_RES_ERR_978							// "LH Reserved error"
#define LH_RES_ERR_979							// "LH Reserved error"
#define LH_RES_ERR_980							// "LH Reserved error"
#define LH_RES_ERR_981							// "LH Reserved error"
#define LH_RES_ERR_982							// "LH Reserved error"
#define LH_RES_ERR_983							// "LH Reserved error"
#define LH_RES_ERR_984							// "LH Reserved error"

#define OH_BAD_STRUCT_SIZE 985					// "The structure is the wrong size."
#define OH_BAD_STRUCT_TYPE_NAME 986				// "This is the wrong type of structure."
#define OH_BAD_STRUCTURE_TYPE 987				// "Expected structure type, 0 or 1."
#define OH_EXPECTED_LITERAL_INTEGER 988			// "Expected a literal integer number."
#define FILE_STRUCT_MISMATCH 989				// "The file size does not match the structure size."

#define EXPECT_GREATER_THAN 990					// "Expected number greater than ^2"
#define FCMD_AS_NOT_ALLOWED 991					// "This command doesn't accept an \"as\" parameter."

#define NO_COMPLEX_TEXT_WAVES 992				// "Igor does not support complex text waves."

#define STRING_ACCESS_ON_NUMERIC_VARIABLE 993	// "An attempt was made to treat a numeric variable as if it were a string variable."
#define NUMERIC_ACCESS_ON_STRING_VARIABLE 994	// "An attempt was made to treat a string variable as if it were a numeric variable."
#define BAD_VARIABLE_DATA_TYPE 995				// "A variable data type must double, double complex, or string."
#define COMPILE_FAILED 996						// "Procedure compilation failed."
#define BAD_WIN_TYPE 997						// "This operation does not apply to this type of window."
#define BAD_PAGESETUP_SCALE 998					// "Expected a page setup scaling value in percent between 5 and 5000."
#define FLAG_ALLOWED_ONLY_ONCE 999				// "This operation allows each flag to be used only once in a single command."
#define KEYWORD_ALLOWED_ONLY_ONCE 1000			// "This operation allows each keyword to be used only once in a single command."

#define	BAD_NUM 1001							// "expected number"
#define	NAM_SYMB_BAD 1002						// "unknown/inappropriate name or symbol"
#define	LPAREN_MISSING 1003						// "expected left parenthesis"
#define	RPAREN_MISSING 1004						// "expected right parenthesis"
#define	OPAND_MISMATCH 1005						// "expected operand"
#define	OPTOR_MISMATCH 1006						// "expected operator"
#define	OPTOR_OPAND_MISMATCH 1007				// "operator/operand mismatch"
#define WRONG_NO_PARAMS 1008					// "wrong number of parameters"
#define BI_NO_LPAREN 1009						// "expected left parenthesis"
#define RPAREN_EXPECTED 1010					// "expected right parenthesis"
#define NON_NULL_PARAMS 1011					// "this function takes no parameters"
#define BI_NO_FCTN_THIS_NT 1012					// "function not available for this number type"
#define AMBIGUOUS_WAVE_POINT 1013				// "ambiguous wave point number"
#define BI_BAD_WAVEFORM 1014					// "expected wave name"
#define NO_SUCH_CURSOR 1015						// "cursor ^0 is not on the graph \"^1\""
#define NAM_SYMB_BAD_NO_CSR 1016				// "expected cursor name (A or B)"
#define LINE_TOO_LONG 1017						// "line too long"
#define RBRAKET_EXPECTED 1018					// "expected ']'"
#define LBRAKET_EXPECTED 1019					// expected '['
#define NO_USER_VARS_IN_FUNC 1020				// user variables are not allowed in functions

// Automation errors
#define AUTOMATION_RESERVED_1021 1021			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1022 1022			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1023 1023			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1024 1024			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1025 1025			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1026 1026			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1027 1027			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1028 1028			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1029 1029			// Reserved for Automation error.
#define IGORPRO_COM_INIT_FAILED 1030			// "Initialization of COM failed."
#define IGORPRO_COM_UNEXPECTED_ERROR 1031		// "An unexpected error was occurred in an internal Automation routine."
#define IGORPRO_COM_PARAM_ERROR 1032			// "Automation client invalid parameter."
#define WAVE_USED_BY_AUTOMATION_IWave 1033						// "Wave is in use by an Automation client (IWave object)."
#define VARIABLE_USED_BY_AUTOMATION_IVariable 1034				// "Variable is in use by an Automation client (IVariable object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumWaves 1035			// "Data folder is in use by an Automation client (IEnumWaves object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IWaves 1036				// "Data folder is in use by an Automation client (IWaves object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumVariables 1037		// "Data folder is in use by an Automation client (IEnumVariables object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IVariables 1038			// "Data folder is in use by an Automation client (IVariables object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IDataFolder 1039			// "Data folder is in use by an Automation client (IDataFolder object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IEnumDataFolders 1040	// "Data folder is in use by an Automation client (IEnumDataFolders object)."
#define DATA_FOLDER_USED_BY_AUTOMATION_IDataFolders 1041		// "Data folder is in use by an Automation client (IDataFolders object)."
#define AUTOMATION_WAVE_KILLED_IWave 1042						// "A wave referenced by an Automation client was killed (IWave)."
#define AUTOMATION_VARIABLE_KILLED_IVariable 1043				// "A variable referenced by an Automation client was killed (IVariable)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumWaves 1044			// "A data folder referenced by an Automation client was killed (IEnumWaves)."
#define AUTOMATION_DATA_FOLDER_KILLED_IWaves 1045				// "A data folder referenced by an Automation client was killed (IWaves)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumVariables 1046		// "A data folder referenced by an Automation client was killed (IEnumVariables)."
#define AUTOMATION_DATA_FOLDER_KILLED_IVariables 1047			// "A data folder referenced by an Automation client was killed (IVariables)."
#define AUTOMATION_DATA_FOLDER_KILLED_IDataFolder 1048			// "A data folder referenced by an Automation client was killed (IDataFolder)."
#define AUTOMATION_DATA_FOLDER_KILLED_IEnumDataFolders 1049		// "A data folder referenced by an Automation client was killed (IEnumDataFolders)."
#define AUTOMATION_DATA_FOLDER_KILLED_IDataFolders 1050			// "A data folder referenced by an Automation client was killed (IDataFolders)."
#define AUTOMATION_RESERVED_1051 1051			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1052 1052			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1053 1053			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1054 1054			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1055 1055			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1056 1056			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1057 1057			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1058 1058			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1059 1059			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1060 1060			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1061 1061			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1062 1062			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1063 1063			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1064 1064			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1065 1065			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1066 1066			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1067 1067			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1068 1068			// Reserved for Automation error.
#define AUTOMATION_RESERVED_1069 1069			// Reserved for Automation error.

#define SUBRANGE_DLOG_BADLABEL 1070				// "Found what appears to be a dimension label, but it is not a label in the selected wave."
#define SUBRANGE_DLOG_WRONGLABEL 1071			// "Found a dimension label, but it is the overall dimension label. You must use a label for a single element."
#define SUBRANGE_DLOG_EXPECTEDNUM 1072			// "Expected a number."
#define SUBRANGE_DLOG_MISSINGRANGE 1073			// "You must specify a range of elements for one dimension."
#define SUBRANGE_DLOG_BADEND 1074				// "Only the range dimension should specify the End."
#define SUBRANGE_DLOG_BADINCREMENT 1075			// "Only the range dimension should specify the Increment."
#define SUBRANGE_DLOG_WRONGNUMBEROFPOINTS 1076	// "Needed ^2 points in range, have ^3 points."
#define SUBRANGE_DLOG_EXPECTNUMORDIMLABEL 1077	// "Expected number or dimension label.",
#define SUBRANGE_DLOG_EXPECTNUM_GE_1 1078		// "Expected a number greater than or equal to 1.",
#define SUBRANGE_DLOG_REQUIREDPNTS_TOOBIG 1079	// "No dimension in the wave is large enough."
#define SUBRANGE_DLOG_BADRANGEDIM 1080			// "Required points too large for selected range dimension."
#define SUBRANGE_DLOG_LABEL_NOT_ALLOWED 1081	// "A dimension label is not allowed in the range dimension."
#define SUBRANGE_DLOG_INC_TOO_BIG 1082			// "This increment causes the range to excede the wave dimension size. It must be less than ^3."

#define CVODE_ILLEGAL_INPUT 1083				// "IntegrateODE reports illegal input for method 2 or 3"
#define CVODE_SETUP_FAILURE 1084				// "IntegrateODE reports setup failure for method 2 or 3"
#define CVODE_SOLVER_FAILURE 1085				// "IntegrateODE reports solver failure for method 2 or 3"
#define CVODE_BUG_ERROR 1086					// "IntegrateODE returned an unknown error code. Please contact WaveMetrics support.

#define CURVEFIT_CONF_NOALLATONCE 1087			// "It is not possible to calculate confidence bands for all-at-once fit functions."
#define DANGEROUS_FILE_COMMAND_DISABLED 1088	// "Command is disabled. It can be enabled in the Miscellaneous Setting Dialog's Misc Settings category."

#define MISMATCHED_NUMBER_OF_POINTS_ADD_BUTTON 1089		// "Different number of points for X and Y waves. Click Add button and edit range in the list below."
#define MISMATCHED_NUMBER_OF_POINTS_MORE_BUTTON 1090	// "Different number of points for X and Y waves. Click More Choices button."
#define EXPECTED_WAVE_SELECTION 1091					// "You must select a wave."

#define BAD_OBJECT_COORDINATES 1092				// Bad object coordinates.
#define CF_PLOTIT_NOEXPOFFSET 1093				// "The Graph Now button does not yet support the built-in functions exp_Xoffset or dbl_Xoffset."
#define EXPECTED_FILE_NAME 1094					// "Expected file name"	

#define CF_WRONGNUMBEROFCONSTANTS 1095			// "Wrong number of constants specified. The chosen fit function takes ^0 constants."
#define EXPECTNUMBERORAUTO 1096					// "Expected a number or \"Auto\""

#define E_FLAG_REQUIRED	1097					// "The /E flag is required for this operation."

#define INCOMPATIBLE_PACKAGE_PREFS_FILE 1098	// "The package preference file is incompatible with this version of Igor."
#define PACKAGE_PREFS_RECORD_NOT_FOUND 1099		// "Package preference record not found."
#define CORRUPT_PACKAGE_PREFS_FILE 1100			// "The package preference file is corrupt. You should delete it."
#define EXPECTED_NONNEGATIVE_INTEGER 1101		// "Expected a non-negative integer."

// From IgorMenus.h
// These are the menu IDs that you can use in XMI1 resources to attach menu items to Igor menus.
#define APPLEID 1
#define FILEID 2
#define EDITID 3
#define WAVEFORMID 4
#define DATAID 4					// HR, 10/2/93 -- old "Waves" menu is now called "Data"
#define ANALYSISID 5
#define MACROID 6
#define WINDOWSID 7
#define MISCID 8
#define LAYOUTID 10					// HR, 10/2/93 -- this was 9 prior to Igor 2.0
#define GRAPHID 12					// HR, 10/2/93 -- added for Igor Pro 2.0
#define PANELID 13					// HR, 10/2/93 -- added for Igor Pro 2.0
#define TABLEID 14					// HR, 10/2/93 -- added for Igor Pro 2.0
#define PROCEDUREID 15				// HR, 10/2/93 -- added for Igor Pro 2.0
#define NOTEBOOKID 16				// HR, 10/2/93 -- added for Igor Pro 2.0
#define LOAD_SUB_ID 50
#define SAVE_SUB_ID 51
// #define SAVEGRAPHICS_SUB_ID 52	// HR, 981105: The Save Graphics submenu was removed in Igor Pro 3.1.
#define OPEN_FILE_SUB_ID 55
#define CONTROL_WIN_SUB_ID 56
#define NEW_WIN_SUB_ID 58
#define MISC_OPS_SUB_ID 59
#define APPEND_TO_GRAPH_SUBID 89	// HR, 3/2/96 -- added for Igor Pro 3.0

// From Parse.h
#define MAXWAVES 8				// maximum number of waves in wave list for some commands

#define T_COMMA		1			// terminator codes
#define T_RPAREN 	2
#define T_SEMI		4
#define T_RBRACK	8
#define T_RCBRACE	16
#define T_NORM		(T_COMMA | T_SEMI)
#define T_CRP		(T_COMMA | T_RPAREN)
#define T_CRB		(T_COMMA | T_RCBRACE)
#define T_CRBRACK	(T_COMMA | T_RBRACK)	// LH, 3/18/90


// From IgorMath.h
#define NT_CMPLX 1				// complex numbers
#define NT_FP32 2				// 32 bit fp numbers
#define NT_FP64 4				// 64 bit fp numbers
#define NT_I8 8					// 8 bit signed integer (changed 1/21/91)
#define NT_I16 	0x10			// 16 bit integer numbers
#define NT_I32 	0x20			// 32 bit integer numbers
#define NT_UNSIGNED 0x40		// Makes above signed integers unsigned. NOTE: Requires Igor Pro 3.0 or later.


// From wave.h
#define MAX_WAVE_NAME 31		// maximum length of wave name -- not including the null
								//	NOTE: Prior to Igor 3.0, this was 18 and we recommended that you use MAX_OBJ_NAME (31) instead of MAX_WAVE_NAME.
							
#define OLD_MAX_WAVE_NAME 18	// MAX_WAVE_NAME prior to Igor Pro 3.0.

#define TEXT_WAVE_TYPE 0		// The wave type code for text waves. Added in Igor Pro 3.0.

#define MAX_DIMENSIONS 10		// Maximum number of dimensions in a multi-dimension object.
								// In Igor 3.0, the max is actually 4 but this may increase in the future.
#define ROWS 0					// Dimension 0 is the row dimension.
#define COLUMNS 1				// Dimension 1 is the column dimension.
#define LAYERS 2				// Dimension 2 is the layer dimension.
#define CHUNKS 3				// Dimension 3 is the chunk dimension.

#define MAX_UNIT_CHARS 49		// Max number of characters in a units string, not including trailing null, in Igor Pro 3.0 or later.
								// Prior to Igor Pro 3.0, the maximum was 3 characters.

#define MAX_DIM_LABEL_CHARS 31	// Max chars in a dimension label, not including trailing null.
								
#define kMDWaveAccessMode0 0	// Access code for MDAccessNumericWaveData. Used by Igor for future compatibility check.

// From WM.h
#define UNKNOWNCURSOR 0
#define ARROWCURSOR 1
#define WATCHCURSOR 2
#define IBEAMCURSOR 3
#define HANDCURSOR 4
#define SPINNINGCURSOR 5
#define CROSSHAIRCURSOR 6
#define MAX_OBJ_NAME 31			// maximum length of: variables,macros,annotations.
#define MAX_LONG_NAME 255		// Added in 6.00D00. Used for double names.
#define WNAMESIZE 40			// maximum length of window name.

#ifdef _MACINTOSH_
	#define MAX_VOLUMENAME_LEN 27			// Maximum length of volume name
	#define MAX_DIRNAME_LEN 31				// Maximum length of directory name
	#define MAX_FILENAME_LEN 31				// Maximum length of file name
	#define MAX_PATH_LEN 511				// Maximum length of path name. This was 511 in Igor 3.0 so I am leaving it as 511 even though it is not clear whether the Mac OS support more than 255.
#endif
#ifdef _WINDOWS_
	#define MAX_VOLUMENAME_LEN 255			// Maximum length of volume name (e.g., "C:")
	#define MAX_DIRNAME_LEN 255				// Maximum length of directory name
	#define MAX_FILENAME_LEN 255			// maximum length of file name
	#define MAX_PATH_LEN 259				// maximum length of path name
#endif

#define F32_NAN			0x7fffffffL
#define F32_PLUS_INF	0x7f800000L
#define F32_MINUS_INF	0xff800000L
#define F32_SIGN_MASK	0x80000000L
#define F32_NOSIGN_MASK	0x7fffffffL
#define F32_INF 		0x7f800000L		// 32 bit infinity
#define F64_INF 		0x7ff00000L		// 64 bit infinity
#define FSIGN_MASK	(long)0x80000000L	// sign bit

/*	This is used to select one of two ways of doing something or to allow both.
	For an example, see VolumeNameLength() in WMFileUtils.c.
*/
typedef enum PlatformCode {
	kMacPlatform=1,				// This is stored on disk. The value must not be changed.
	kWinPlatform,				// This is stored on disk. The value must not be changed.
	kMacOrWinPlatform,
#ifdef _MACINTOSH_
	kCurrentPlatform=kMacPlatform
#endif
#ifdef _WINDOWS_
	kCurrentPlatform=kWinPlatform
#endif
} PlatformCode;

#define CR_STR "\015"			// Can be used as follows: XOPNotice("Test"CR_STR);
#define CR_CHAR '\015'
#define LF_STR "\012"
#define LF_CHAR '\012'

#define LEFT_ARROW_CHAR_CODE 0x1C				// Low byte of message field of EventRecord.
#define RIGHT_ARROW_CHAR_CODE 0x1D
#define UP_ARROW_CHAR_CODE 0x1E
#define DOWN_ARROW_CHAR_CODE 0x1F
#define PAGE_UP_CHAR_CODE 0x0B
#define PAGE_DOWN_CHAR_CODE 0x0C
#define HOME_CHAR_CODE 0x01
#define END_CHAR_CODE 0x04
#define FORWARD_DELETE_CHAR_CODE 0x7F
#define HELP_CHAR_CODE 0x05

#define LEFT_ARROW_KEY_CODE 0x7B				// Second byte of message field of EventRecord.
#define RIGHT_ARROW_KEY_CODE 0x7C
#define UP_ARROW_KEY_CODE 0x7E
#define DOWN_ARROW_KEY_CODE 0x7D
#define PAGE_UP_KEY_CODE 0x74
#define PAGE_DOWN_KEY_CODE 0x79
#define HOME_KEY_CODE 0x73
#define END_KEY_CODE 0x77
#define FORWARD_DELETE_KEY_CODE 0x75
#define HELP_KEY_CODE 0x72

#define MODIFIER_KEYS_MASK (cmdKey | optionKey | controlKey | shiftKey)


// From Functions.h

// These are used to identify parameters to external functions as waves, strings or names
#define WAVE_TYPE 0x4000		// added to number types above to signify parameter is wave
#define HSTRING_TYPE 0x2000		// signifies parameter is a handle to a string
#define HNAME_TYPE 0x1000		// signifies parameter is a handle to a name

// These are used to test parameter types returned by GetUserFunctionInfo.
#define FV_REF_TYPE 0x1000		// Signifies pass-by-reference
#define FV_FUNC_TYPE 0x0400		// Signifies a function reference
#define WAVE_Z_TYPE	0x8000		// Identifies WAVE/Z argument
#define FV_STRUCT_TYPE 0x0200	// Requires Igor Pro 5.03 or later

struct NVARRec {				// Used for NVAR structure fields.
	long urH;
	long index;
};
typedef struct NVARRec NVARRec;

struct SVARRec {				// Used for SVAR structure fields.
	long urH;
	long index;
};
typedef struct SVARRec SVARRec;

typedef long FUNCREF;			// Used for FUNCREF structure fields.

struct FunctionInfo {			// Used by GetUserFunctionInfo.
	char name[MAX_OBJ_NAME+1];
	int compilationIndex;
	int functionID;
	int subType;
	int isExternalFunction;
	int returnType;
	int reserved[25];					// Do not use. Reserved for future use.
	int numOptionalParameters;
	int numRequiredParameters;
	int totalNumParameters;
	int parameterTypes[100];
};
typedef struct FunctionInfo FunctionInfo;
typedef FunctionInfo* FunctionInfoPtr;


// From TextUtils.h
#define TUMAXBYTES 32000	// max # of bytes allowed in TU document
							// 10/2/93: TUMAXBYTES is obsolete. There is no longer any limit.
							//			We leave TUMAXBYTES here because existing XOP source may reference it.


// structure for getting info about text utility document
struct TUDocInfo {			// 10/23/93: added for Igor Pro 2.0D83
	short version;						// version number of this structure
	short permission;					// 0 = read only, 1 = read/write
	short fileType;						// for future use
	long paragraphs;					// total number of paragraphs in document
	char reserved[256];					// for future use
};
typedef struct TUDocInfo TUDocInfo;
typedef struct TUDocInfo *TUDocInfoPtr;
#define TUDOCINFO_VERSION 1

struct TULoc {							// identifies a location in a text utility document
	long paragraph;						// location's paragraph
	unsigned short pos;					// character offset in paragraph for text paragraph
};
typedef struct TULoc TULoc;
typedef struct TULoc *TULocPtr;

/*	When to erase message in status area.
	This is a bitwise parameter used with the TUSetStatus() callback.
	The status is always changed if a new message comes along.
	This controls if and when it will be erased before a new message comes.
*/
#define TU_ERASE_STATUS_NEVER 0
#define TU_ERASE_STATUS_WHEN_SELECTION_CHANGES 1
#define TU_ERASE_STATUS_WHEN_WINDOW_ACTIVATED 2
#define TU_ERASE_STATUS_WHEN_WINDOW_DEACTIVATED 4
#define TU_ERASE_STATUS_WHEN_DOC_MODIFIED 8
#define TU_ERASE_STATUS_WHEN_ANYTHING_HAPPENS -1


// From CmdWin.h
// modes for PutCmdLine()
#define INSERTCMD 1					// insert text at current insertion point
#define FIRSTCMD 2					// insert text in front of cmd buffer
#define FIRSTCMDCRHIT 3				// insert text in front of cmd buffer and set crHit
#define REPLACEFIRSTCMD 4			// replace first line of cmd buffer with text
#define REPLACEALLCMDSCRHIT 5 		// replace all lines of cmd buffer with text and set crHit
#define REPLACEALLCMDS 6			// replace all lines of cmd buffer with text


// From ColorTables.h
typedef Handle IgorColorTableHandle;
typedef struct IgorColorSpec {
	unsigned long value;			// index or other value
	RGBColor rgb;					// true color
} IgorColorSpec;


// From CommandUtils.h
// structure for getting and passing wave range information
struct WaveRangeRec {
	// the following fields are set by GetWaveRange based on command line
	double x1, x2;	// *** 4/21/90 -- changed from float to double
	long rangeMode;					// bit 0 set if start specified, bit 1 set if end specified
	long isBracket;					// true if range specified by [] instead of ()
	long gotRange;					// true if /R=range was present

	// next, you setup these fields
	Handle waveHandle;
	long minPoints;					// min number of points in acceptable range
	
	// Then, following fields are setup by CalcWaveRange
	long p1, p2;
	long wasBackwards;				// truth p1 > p2 before being swapped
};
typedef struct WaveRangeRec WaveRangeRec;
typedef struct WaveRangeRec *WaveRangeRecPtr;


// From CHIO.h
#define CHIOREAD 1
#define CHIOREADWAVE 2
#define CHIOWRITE 3
#define CHIOWRITEWAVE 4

struct CHIORec {
	short operation;
	long info;							// 0 or XOPRecHandle for XOPs
	short isXOP;						// TRUE if read request is from XOP
	long refCon;						// for use by calling program
	long igorPrivate;					// for use by Igor
	long (*checkProc)(struct CHIORec*);	// address of routine to check # chars available
	long (*readProc)(struct CHIORec*);	// address of routine to read characters
	long (*writeProc)(struct CHIORec*);	// address of routine to write characters
	long (*returnProc)(struct CHIORec*);// address of routine to return characters
	char tf[256];						// terminator string for read operations, format for write
	char *bufPtr;						// pointer to buffer to read into or write from
	long *lengthPtr;					// pointer to long containing # of chars to read or write
	short debug;						// if true, reports debugging info
	short maxChars;						// max chars for input from character input device
	long timeout;						// max ticks allowed for each read or zero for no timeout
	long timeoutTicks;					// ticks at which current read or write should timeout
	short quiet;						// if non-zero, does not return error if timeout occurs
	long flags;							// see flags #defines below
	long items;							// # of items read
};
typedef struct CHIORec CHIORec;
typedef CHIORec *CHIORecPtr;
typedef CHIORecPtr *CHIORecHandle;

// CHIO flags #defines
#define CHIO_WRITE_HEADER 1				// not yet implemented


// From Variables.h
#define VAR_GLOBAL	0x4000				// bit flag for type parameter of Variable XOP callback

struct NumVarValue{						// Used in Get/SetDataFolderObject call.
	long numType;			// NT_FP64 possibly ORed with NT_CMPLX (if variable is complex).
	long spare;				// For future use - set to zero.
	double realValue;
	double imagValue;
};
typedef struct NumVarValue NumVarValue;
typedef NumVarValue* NumVarValuePtr;

union DataObjectValue {
	waveHndl wavH;						// Use this if the object is a wave.
	NumVarValue nv;						// Use this if the object is a numeric variable.
	Handle strH;						// Use this if the object is a string variable.
	DataFolderHandle dfH;				// Use this if the object is a data folder.
	char spare[64];						// For possible future use.
};
typedef union DataObjectValue DataObjectValue;
typedef DataObjectValue* DataObjectValuePtr;


// From Save.h
#define SAVE_TYPE_SAVE 1				// experiment save type codes
#define SAVE_TYPE_SAVEAS 2
#define SAVE_TYPE_SAVEACOPY 3
#define SAVE_TYPE_STATIONERY 4

#define LOAD_TYPE_NEW 1					// experiment load type codes
#define LOAD_TYPE_OPEN 2
#define LOAD_TYPE_REVERT 3
#define LOAD_TYPE_STATIONERY 4
#define LOAD_TYPE_MERGE 5				// Added in Igor Pro 5.

#define EXP_UNPACKED 0					// experiment file type codes
#define EXP_PACKED 1


// From OperationHandler.h
#define kUseCMDMessageForInterpreting 1

struct DataFolderAndName {
	DataFolderHandle dfH;
	char name[MAX_OBJ_NAME+1];
};
typedef struct DataFolderAndName DataFolderAndName;
typedef struct DataFolderAndName *DataFolderAndNamePtr;

struct WaveRange {
	waveHndl waveH;
	double startCoord;					// Start point number or x value
	double endCoord;					// End point number or x value
	int rangeSpecified;					// 0 if user specified no range. 1 if user specified range.
	int isPoint;						// 0 if user specified range using X values. 1 if user specified range using points.
};
typedef struct WaveRange WaveRange;
typedef struct WaveRange *WaveRangePtr;

// This is what is in the runtime parameter structure for a mode 1 structure parameter.
// For a mode 0 structure parameter, the runtime parameter structure contains just a pointer to the structure.
struct IgorStructInfo {
	void* structPtr;					// Pointer to the structure.
	unsigned long structSize;			// Size of structure in bytes.
	char structTypeName[MAX_OBJ_NAME+1];
	int moduleSerialNumber;				// Used by Igor to determine the procedure file in which structure is defined.
	unsigned char reserved[32];			// Reserved for future use.
};
typedef struct IgorStructInfo IgorStructInfo;
typedef struct IgorStructInfo *IgorStructInfoPtr;

#include "XOPStructureAlignmentReset.h"	// Reset structure alignment to default.
