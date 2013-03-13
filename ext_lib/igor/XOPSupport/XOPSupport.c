/*	XOPSupport.c
		Support routines for Igor XOPs
*/

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

// Global Variables.
IORecHandle XOPRecHandle;	// The XOP's ioRecHandle.
int XOPModified = 0;
int igorVersion = 0;		// Set by XOPInit. Example: 128 for version 1.28.


// *** Utility Routines ***

/*	Capitalize(name)

	Changes any lower case letters in name to upper case.
*/
void
Capitalize(char *name)
{
	int len, i;
	char *p;
	
	len = strlen(name);
	p = name;
	for (i = 0;i < len; i++) {
		if (*p >= 'a' && *p <= 'z')
			*p -= 'a' - 'A';
		p++;
	}
}

/*	CmpStr(char *str1, char *str2)

	Compares two C strings case insensitively.
	
	Result is the same as for the strcmp function.
	
	HR, 980331: Removed the constraint that the strings be 255 characters or less.
*/
int
CmpStr(const char *s1, const char *s2)
{
	unsigned const char *str1= (unsigned const char *) s1;
	unsigned const char *str2= (unsigned const char *) s2;
	unsigned int c1,c2;
	int result= 0;						/* assume s1 == s2 */
	
	do {
		c1= *str1++;
		if( c1 >= 'a' && c1 <= 'z' )	/* if (islower(c1))     */
			c1 -= 'a' - 'A';			/*     c1= toupper(c1); */
	
		c2= *str2++;
		if( c2 >= 'a' && c2 <= 'z' )	/* if (islower(c2))     */
			c2 -= 'a' - 'A';			/*     c2= toupper(c2); */

		if( c1 != c2 ) {
			if( c1 < c2 )
				result= -1;				/* s1 < s2 */
			else
				result= 1;				/* s1 > s2 */
			break;
		}

	} while ( c1 );	/* falls through if c1 == 0 (and c2 == 0) because of above if(); s1 == s2 */
	
	return result;
}

/*	strchr2(str, ch)

	strchr2 is like the standard C strchr function except that it is
	Asian-language-aware. 
	
	Returns a pointer to the first occurrence of ch in the null-terminated string
	str or NULL if there is no such occurrence. 
	
	On a system that uses an Asian script system as the default script, strchr2
	knows about two-byte characters. For example, if you are searching for a
	backslash in a full path and if the path contains Asian characters, and if the
	second byte of an Asian character has the same code as the backslash character,
	strchr will mistakenly find this second byte while strchr2 will not. 
	
	On a system that does not use an Asian script system as the default script,
	strchr2 is just like strchr. 
	
	Added in Igor Pro 3.13. If you call this when running with an earlier version,
	it will behave just like strchr. 
*/
char*
strchr2(const char* str, int ch)
{
	if (igorVersion < 313)
		return strchr(str, ch);
	return (char*)CallBack2(STRCHR2, (void*)str, (void*)ch);
}

/*	strrchr2(str, ch)

	strrchr2 is like the standard C strrchr function except that it is
	Asian-language-aware. 
	
	Returns a pointer to the last occurrence of ch in the null-terminated string
	str or NULL if there is no such occurrence. 
	
	On a system that uses an Asian script system as the default script, strrchr2
	knows about two-byte characters. For example, if you are searching for a
	backslash in a full path and if the path contains Asian characters, and if the
	second byte of an Asian character has the same code as the backslash character,
	strrchr will mistakenly find this second byte while strrchr2 will not.
	
	On a system that does not use an Asian script system as the default script,
	strrchr2 is just like strrchr. 
	
	Added in Igor Pro 3.13. If you call this when running with an earlier version,
	it will behave just like strrchr. 
*/
char*
strrchr2(const char* str, int ch)
{
	if (igorVersion < 313)
		return strrchr(str, ch);
	return (char*)CallBack2(STRRCHR2, (void*)str, (void*)ch);
}

/*	MemClear(void *p, long n)

	p points to start of memory to clear.
	n is number of bytes to clear.
*/
void
MemClear(void *p, long n)
{
	long i;
	char *pp;
	
	pp = p;
	for (i = 0; i < n; i++)
		*pp++ = 0;
}

/*	MoveLockHandle(h)

	h is a handle.
	Moves handle to top of heap and locks it.
	Returns state of handle before locking.
*/
int
MoveLockHandle(void *h)
{
	int result;
	
	result = HGetState((Handle)h);
	MoveHHi((Handle)h);
	HLock((Handle)h);
	return(result);
}

/*	GetCStringFromHandle(h, str, maxChars)

	h is a handle containing a string.

	str is a C string (null-terminated character array).

	maxChars is the maximum number of characters that str can hold, not including the
	null terminator byte.
	
	GetCStringFromHandle transfers the characters from h to str.
	
	If h is NULL, it returns USING_NULL_STRVAR. This is typically a programmer error.
	
	If the characters in h will not fit in str, it returns STR_TOO_LONG.
	
	If the characters fit, it returns 0.
*/
int
GetCStringFromHandle(Handle h, char* str, int maxChars)
{
	int numBytesInString;
	int err;
	
	err = 0;
	
	*str = 0;

	if (h == NULL)
		return USING_NULL_STRVAR;

	numBytesInString = GetHandleSize(h);
	if (numBytesInString > maxChars) {
		numBytesInString = maxChars;
		err = STR_TOO_LONG;
	}
	
	memcpy(str, *h, numBytesInString);
	str[numBytesInString] = 0;
	
	return err;
}

/*	PutCStringInHandle(str, h)

	str is a C string (null-terminated character array).

	h is a handle in which the C string data is to be stored.
	
	PutCStringInHandle transfers the characters from str to h. Note that
	the trailing null from the C string is not stored in the handle.
	
	If h is NULL, it returns USING_NULL_STRVAR. This is typically a programmer error.
	
	If an out-of-memory occurs when resizing the handle, it returns NOMEM.
	
	If the operation succeeds, it returns 0.
*/
int
PutCStringInHandle(const char* str, Handle h)
{
	int numBytesInString;

	if (h == NULL)
		return USING_NULL_STRVAR;
		
	numBytesInString = strlen(str);
	SetHandleSize(h, numBytesInString);
	if (MemError())
		return NOMEM;
	
	memcpy(*h, str, numBytesInString);
	
	return 0;
}

/*	CheckAbort(timeoutTicks)

	Returns -1 if user is now pressing cmd-dot (Macintosh) or Ctrl-Break (Windows).
	Returns 1 if TickCount > timeoutTicks.
	Returns 0 otherwise.
	However, if timeoutTicks == 0, does not check for timeout.
	
	Actually does check only every .1 second.
*/
int
CheckAbort(long timeOutTicks)
{
	long ticks;
	static unsigned long lastTicks = 0;
	
	ticks = TickCount();
	if (ticks < lastTicks+6)
		return(0);
	lastTicks = ticks;	
	if (timeOutTicks && ticks>timeOutTicks)
		return 1;					// Timeout.

	#ifdef _MACINTOSH_					// Test for cmd-dot.
	{
		KeyMap keys;
		UInt32 theKeys[4];
	
		/*	The KeyMap data type is defined weird on MacIntel which makes it hard to use.
			We copy the data to an array of 4 UInt32s here. The data will be big-endian
			even when running on MacIntel.
			
			If running on MacIntel we convert to little-endian because we are treating
			the data as an array of UInt32s. If we treated it as an array of bytes then
			we would not want to byte swap it.
		*/
		GetKeys(keys);
		memcpy(theKeys, &keys, sizeof(theKeys));
		#ifdef __LITTLE_ENDIAN__
			FixByteOrder(theKeys, 4, 4);
		#endif
	
		// HR, 12/5/93. Added check so that user procedures can use option-cmd-dot for abort signal.
		// HR, 2/22/94. Changed to allow cmd-dot or shift-cmd-dot because of European keyboards.
		if (theKeys[1]==0x00808000 || theKeys[1]==0x00808001) {		// Cmd-dot or shift-cmd-dot ?
			if (theKeys[0]==0 && theKeys[2]==0 && theKeys[3]==0)	// No other keys pressed ?
				return -1;
		}
	}
	#endif	
	#ifdef _WINDOWS_					// Test for Ctrl-Break.
	{
		/*	Control-break is treated like cmd-dot on Mac.
			Alt-control-break is reserved for user procedures.
			
			HR, 9/25/97: For some reason, GetAsyncKeyState does not always return
			a value with the high bit set if I press and hold Ctrl-Break. Thus,
			I have changed the mask from 0x8000 (key down now) to 0x8001 (key down now
			or key was pressed since last call to GetAsyncKeyState).
		*/
		if ((GetAsyncKeyState(VK_CANCEL) & 0x8001) != 0) {		// This tests for Ctrl-Break.
			if ((GetAsyncKeyState(VK_MENU) & 0x8001) == 0)		// Alt key must not be pressed.
				return -1;
		}
	}
	#endif
	return 0;
}

int
IsNaN32(float *floatPtr)		// Returns truth that floatPtr points to a NaN.
{
	return(((*(long *)floatPtr) & 0x7fffffff) > F32_INF);
}

int
IsNaN64(double *doublePtr)		// Returns truth that doublePtr points to a NaN.
{
	#ifdef XOP_LITTLE_ENDIAN
		return(((*((long *)doublePtr+1)) & 0x7fffffff) > F64_INF);
	#else
		return(((*(long *)doublePtr) & 0x7fffffff) > F64_INF);
	#endif
}

void
SetNaN32(float* fPtr)
{
	*fPtr = SINGLE_NAN;
}

void
SetNaN64(double* dPtr)
{
	*dPtr = DOUBLE_NAN;
}

int
IsINF32(float *floatPtr)		// Returns truth that floatPtr points to a +/- INF.
{
	return ((*(long *)floatPtr) & 0x7FFFFFFF) == F32_INF;
}

int
IsINF64(double *doublePtr)		// Returns truth that doublePtr points to a +/- INF.
{
	#ifdef XOP_LITTLE_ENDIAN
		return ((*((long *)doublePtr+1)) & 0x7FFFFFFF) == F64_INF;
	#else
		return ((*(long *)doublePtr) & 0x7FFFFFFF) == F64_INF;
	#endif
}

int
IgorVersion(void)	// Returns Igor version number times 100 (e.g. 1.13 is returned as 113).
{
	#ifdef _WINDOWS_
		HMODULE igorModule;
		char igorName[MAX_FILENAME_LEN+1];
		char* versionBuffer;
		DWORD versionInfoSize;
		DWORD dwHandle;							// A dummy variable for GetFileVersionInfoSize.
		VS_FIXEDFILEINFO* vsp;
		UINT vsLen;
		int units, tenths, hundredths;
		int version;
		
		igorModule = IgorModule();
		if (igorModule == NULL)					// Should not happen.
			return 0;
		if (GetModuleFileName(igorModule, igorName, MAX_FILENAME_LEN+1) == 0)
			return 0;
		versionInfoSize = GetFileVersionInfoSize(igorName, &dwHandle);
		if (versionInfoSize <= 0)
			return 0;
		versionBuffer = NewPtr(versionInfoSize);
		if (versionBuffer == NULL)
			return 0;
		if (GetFileVersionInfo(igorName, 0L, versionInfoSize, versionBuffer) == 0) {
			DisposePtr(versionBuffer);
			return 0;
		}
		if (VerQueryValue(versionBuffer, "\\", (void**)&vsp, &vsLen) == 0) {
			DisposePtr(versionBuffer);
			return 0;
		}
	
		units = vsp->dwFileVersionMS >> 16;
		tenths = vsp->dwFileVersionMS & 0xFFFF;
		hundredths = vsp->dwFileVersionLS >> 16;
		version = 100*units + 10*tenths + hundredths;
		
		DisposePtr(versionBuffer);
		return version;
	#endif
	
	#ifdef _MACINTOSH_
		int curResFile;
		Handle vHandle;
		char *vPtr;
		int tens, units, tenths, hundredths;
		int version = 0;
		
		curResFile = CurResFile();
		UseResFile((*(*XOPRecHandle)->stuffHandle)->oldApRefNum);	// Use Igor's resource fork.
		vHandle = Get1Resource('vers', 1);
		UseResFile(curResFile);
		if (vHandle) {							// Pick version out of BCD code.
			vPtr = *vHandle;
			#ifdef __LITTLE_ENDIAN__
				tens = (vPtr[3] & 0xF0) >> 4;
				units = vPtr[3] & 0x0F;
				tenths = (vPtr[2] & 0xF0) >> 4;
				hundredths = vPtr[2] & 0x0F;
			#else
				tens = (vPtr[0] & 0xF0) >> 4;
				units = vPtr[0] & 0x0F;
				tenths = (vPtr[1] & 0xF0) >> 4;
				hundredths = vPtr[1] & 0x0F;
			#endif
			version = 1000*tens + 100*units + 10*tenths + hundredths;
		}
		return(version);
	#endif
}

/*	XOPInit(ioRecHandle)

	Does initialization common to all XOPs. ioRecHandle is the parameter passed by host
	application to XOP.
	
	NOTE: XOPInit() must be called before any other call to XOPSupport.
*/
void
XOPInit(IORecHandle ioRecHandle)
{
	XOPRecHandle = ioRecHandle;		// Set global rec handle.

	igorVersion = IgorVersion();	// igorVersion global added 10/9/93.
}

/*	SetXOPType(type)

	Sets XOPType field of XOP's ioRecHandle. This informs host of capabilities and/or mode of
	XOP.
*/
void
SetXOPType(long type)
{
	(*XOPRecHandle)->XOPType = type;
}

/*	SetXOPEntry(entryPoint)

	Sets XOPEntry field of XOP's ioRecHandle. This informs host of routine to call to pass
	messages to XOP after the INIT message.
*/
void
SetXOPEntry(void (*entryPoint)(void))
{
	(*XOPRecHandle)->XOPEntry = entryPoint;
}

/*	SetXOPResult(result)

	Sets the result field of the XOP's ioRecHandle.
*/
void
SetXOPResult(long result)
{
	(*XOPRecHandle)->result = result;
}

/*	GetXOPResult()

	Returns the result field of the XOP's ioRecHandle.
*/
long
GetXOPResult(void)
{
	return((*XOPRecHandle)->result);
}

void
SetXOPMessage(int message)
{
	(*XOPRecHandle)->message = message;
}

/*	GetXOPMessage()

	Returns the message field of the XOP's ioRecHandle.
*/
long
GetXOPMessage(void)
{
	return((*XOPRecHandle)->message);
}

/*	SetXOPRefCon(refCon)

	Sets the refCon field of the XOP's ioRecHandle. The XOP can use the refCon field to store
	any thing it wants.
*/
void
SetXOPRefCon(long refCon)
{
	(*XOPRecHandle)->refCon = refCon;
}

/*	GetXOPRefCon()

	Returns the refCon field of the XOP's ioRecHandle.
*/
long
GetXOPRefCon(void)
{
	return((*XOPRecHandle)->refCon);
}

/*	GetXOPStatus()

	Returns status field of XOP's ioRecHandle.
*/
long
GetXOPStatus(void)
{
	return((*XOPRecHandle)->status);
}

/*	GetXOPItem(itemNumber)

	Returns an item from the XOP's ioRecHandle.
	itemNumber is the number of the item to return starting from zero.
*/
long
GetXOPItem(int itemNumber)
{
	long item = 0;
	
	if ((*XOPRecHandle)->numItems > itemNumber)				// Make sure item exists.
		item = (*XOPRecHandle)->items[itemNumber];
	return(item);
}

/*	SetRecHandle(numItems)

	Given a handle to an IORec, SetRecHandle sets the number of items in the IORec to numItems.
	
	HR, 3/2/95, for Igor Pro 2.1.
		We no longer make the handle smaller. See NUM_IOREC_ITEMS in XOP.h for
		an explanation.
*/
static void
SetRecHandle(int numItems)
{
	SetHandleSize((Handle)XOPRecHandle, (sizeof(IORec)+(NUM_IOREC_ITEMS-1)*sizeof(long)));
	if (MemError()) {
		// We will never get here unless the heap is fragmented.
		// We try again with just the number of items that we need.
		SetHandleSize((Handle)XOPRecHandle, (sizeof(IORec)+(numItems-1)*sizeof(long)));
		if (MemError())
			debugstr("SetRecHandle error");
	}
	(*XOPRecHandle)->numItems = numItems;
}


// *** Callback Routines --- for requesting service from Igor ***

/*	CallBack(message)

	Does a callback to Igor passing it the XOPRecHandle after clearing the result field.
	
	messages specifies the requested operation from Igor.
	
	returns result field from XOPRecHandle after callback.
*/
static long
CallBack(int message)
{
	#ifdef _MACINTOSH_
		pascal void (*CallBack)(void*);			// Address of Igor's callback routine.
	#endif
	#ifdef _WINDOWS_
		void (*CallBack)(void*);				// Address of Igor's callback routine.
	#endif

	SetXOPResult(0L);							// Start with no error.
	SetXOPMessage(message);
	CallBack = (*XOPRecHandle)->callBackProc;
	(*CallBack)(XOPRecHandle);
	return(GetXOPResult());
}

/*	CallBack0(message)

	Does a callback with the specified message which has 0 parameters.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack0(int message)
{
	SetRecHandle(0);
	return(CallBack(message));
}

/*	CallBack1(message, item)

	Does a callback with the specified message passing the item as the only
	parameter in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack1(int message, void *item)
{
	SetRecHandle(1);
	(*XOPRecHandle)->items[0] = (long)item;
	return(CallBack(message));
}

/*	CallBack2(message, item0, item1)

	Does a callback with the specified message passing the item0 and item1 as
	parameters in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack2(int message, void *item0, void *item1)
{
	SetRecHandle(2);
	(*XOPRecHandle)->items[0] = (long)item0;
	(*XOPRecHandle)->items[1] = (long)item1;
	return(CallBack(message));
}

/*	CallBack3(message, item0, item1, item2)

	Does a callback with the specified message passing item0, item1 and item1 as
	parameters in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack3(int message, void *item0, void *item1, void *item2)
{
	SetRecHandle(3);
	(*XOPRecHandle)->items[0] = (long)item0;
	(*XOPRecHandle)->items[1] = (long)item1;
	(*XOPRecHandle)->items[2] = (long)item2;
	return(CallBack(message));
}

/*	CallBack4(message, item0, item1, item2, item3)

	Does a callback with the specified message passing item0, item1, item2 and item3 as
	parameters in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack4(int message, void *item0, void *item1, void *item2, void *item3)
{
	SetRecHandle(4);
	(*XOPRecHandle)->items[0] = (long)item0;
	(*XOPRecHandle)->items[1] = (long)item1;
	(*XOPRecHandle)->items[2] = (long)item2;
	(*XOPRecHandle)->items[3] = (long)item3;
	return(CallBack(message));
}

/*	CallBack5(message, item0, item1, item2, item3, item4)

	Does a callback with the specified message passing item0 through item4 as
	parameters in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack5(int message, void *item0, void *item1, void *item2, void *item3, void *item4)
{
	SetRecHandle(5);
	(*XOPRecHandle)->items[0] = (long)item0;
	(*XOPRecHandle)->items[1] = (long)item1;
	(*XOPRecHandle)->items[2] = (long)item2;
	(*XOPRecHandle)->items[3] = (long)item3;
	(*XOPRecHandle)->items[4] = (long)item4;
	return(CallBack(message));
}

/*	CallBack6(message, item0, item1, item2, item3, item4, item5)

	Does a callback with the specified message passing item0 through item5 as
	parameters in the thingList.
	
	returns result field from XOPRecHandle after callback.
*/
long
CallBack6(int message, void *item0, void *item1, void *item2, void *item3, void *item4, void *item5)
{
	SetRecHandle(6);
	(*XOPRecHandle)->items[0] = (long)item0;
	(*XOPRecHandle)->items[1] = (long)item1;
	(*XOPRecHandle)->items[2] = (long)item2;
	(*XOPRecHandle)->items[3] = (long)item3;
	(*XOPRecHandle)->items[4] = (long)item4;
	(*XOPRecHandle)->items[5] = (long)item5;
	return(CallBack(message));
}

/*	IgorError(title, errCode)

	Displays an error alert appropriate for the specified error code.
	
	Title is a short string that identifies what generated the error.
	
	errCode may be an Igor error code (defined in IgorXOP.h), an XOP-defined error code,
	or, when running on Macintosh, a Mac OS error code. To display a message for a Windows
	OS error, convert the code to an Igor code by calling WindowsErrorToIgorError.
*/
void
IgorError(const char *title, int errCode)
{
	CallBack2(IGORERROR, (void*)title, (void *)errCode);
}

/*	GetIgorErrorMessage(errCode, errorMessage)

	Returns via errorMessage the message corresponding to the specified error code.
	
	errCode may be an Igor error code (defined in IgorXOP.h), an XOP-defined error code,
	or, when running on Macintosh, a Mac OS error code. To display a message for a Windows
	OS error, convert the code to an Igor code by calling WindowsErrorToIgorError.

	Do not pass 0 for errCode. There is no error message corresponding to 0.

	The function result is 0 if OK, IGOR_OBSOLETE if the current version of Igor
	is earlier than 3.14, or another non-zero error code if the errCode parameter
	is invalid. If GetIgorErrorMessage fails to get a message, it sets *errorMessage to 0.
	
	Added in Igor Pro 3.14. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
*/
int
GetIgorErrorMessage(long errCode, char errorMessage[256])
{
	int err;
	
	err = CallBack2(GET_IGOR_ERROR_MESSAGE, (void*)errCode, errorMessage);
	if (err != 0)
		*errorMessage = 0;
	return err;
}

int
WinInfo(int index, int typeMask, char *name, XOP_WINDOW_REF* windowRefPtr)
{
	return(CallBack4(WININFO, (void *)index, (void *)typeMask, name, windowRefPtr));
}


// *** Notice Routines -- for displaying messages in the history area ***

/*	XOPNotice(noticePtr)

	XOPNotice does a callback to Igor to get the notice identified by noticePtr displayed.
	
	noticePtr is a C string (null byte at the end)
	Typically notices should be terminated by carriage return (use \r at end of string).
*/
void
XOPNotice(const char *noticePtr)
{
	CallBack1(NOTICE, (void*)noticePtr);
}

/*	XOPResNotice(strListID, index)

	XOPResNotice does a callback to Igor to get a notice displayed. The notice is fetched
	from a resource which must be in the resource fork of the XOP.
	
	The resource must be of type 'STR#'. The resource ID should be between 1100 and 1199.
	These resource IDs are reserved for XOPs.
	
	strListID is the resource ID of the STR# containing the string.
	index is the number of the string in the STR# resource.
	
	The strings in the STR# resource must be 255 characters or less.
*/
void
XOPResNotice(int strListID, int index)
{
	char theString[256];

	GetXOPIndString(theString, strListID, index);
	XOPNotice(theString);
}


// *** Command Routines -- for executing Igor commands from an XOP ***

/*	XOPCommand(cmdPtr)

	XOPCommand does a callback to Igor to execute the command identified by cmdPtr.
	The command appears in the command line while Igor is executing it.
	
	cmdPtr is a C string (null byte at the end)
	XOPCommand returns the result from the command execution.
*/
int
XOPCommand(const char *cmdPtr)
{
	return(CallBack1(COMMAND, (void*)cmdPtr));
}

/*	XOPSilentCommand(cmdPtr)

	XOPSilentCommand does a callback to Igor to execute the command identified
	by cmdPtr. The command does not appear in the command line while Igor is
	executing it.
	
	cmdPtr is a C string (null byte at the end)
	XOPSilentCommand returns the result from the command execution.
*/
int
XOPSilentCommand(const char *cmdPtr)
{
	return(CallBack1(SILENT_COMMAND, (void*)cmdPtr));
}

/*	XOPCommand2(cmdPtr, silent, sendToHistory)

	XOPCommand2 does a callback to Igor to execute the command identified by cmdPtr.
	
	cmdPtr is a C string (null byte at the end). The string must consist of one
	line of text not longer than MAXCMDLEN and with no carriage return characters.
	
	If silent is non-zero, the command appears in the command line while Igor is
	executing it.
	
	If sendToHistory is non-zero and if the result from command execution is 0 (no error),
	Igor appends the command to the history with a bullet character before it.

	XOPCommand2 returns the result from the command execution.

	Added in Igor Pro 4.0. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
*/
int
XOPCommand2(const char *cmdPtr, int silent, int sendToHistory)
{
	return(CallBack3(COMMAND2, (void*)cmdPtr, (void*)silent, (void*)sendToHistory));
}

/*	PutCmdLine(cmd, mode)

	PutCmdLine puts the specified text into Igor's command line using the specified mode.
	See IgorXOP.h for list and description of modes.
*/
void
PutCmdLine(const char *cmd, int mode)
{
	CallBack2(PUTCMDLINE, (void*)cmd, (void *)mode);
}

/*	FinishDialogCmd(char *cmd, int mode)

	Called at the end of an Igor style dialog.
	cmd is a C string containing the command generated by the Igor style dialog.
	If mode is 1, puts the command in Igor's command line and starts execution.
	If mode is 2, puts the command in Igor's command line but does not start execution.
	If mode is 3, puts the command in the clipboard.
*/
void
FinishDialogCmd(const char *cmd, int mode)
{
	switch (mode) {
		case 1:
			PutCmdLine(cmd, FIRSTCMDCRHIT);
			break;
		case 2:
			PutCmdLine(cmd, INSERTCMD);
			break;
		case 3:
			#ifdef _MACINTOSH_
			{
				ScrapRef scrap;
				int err;
				
				if (err = ClearCurrentScrap())
					return;
				if (err = GetCurrentScrap(&scrap))
					return;
				if (err = PutScrapFlavor(scrap, 'TEXT', kScrapFlavorMaskNone, strlen(cmd), cmd))
					return;
			}
			#endif
			#ifdef _WINDOWS_
				SetScrapData(NULL, strlen(cmd), 'TEXT', (char*)cmd);	// This routine is in IGOR.
			#endif
			break;
	}
}

//  *** Variable Access Routines ***

/*	FetchNumVar(varName, doublePtr1, doublePtr2)

	FetchNumVar returns the value of a named variable via the double pointers doublePtr1
	and doublePtr2. The real part is returned via doublePtr1 and the imaginary part via
	doublePtr2. If the numeric variable is not complex then the *doublePtr2 is meaningless.
	However, doublePtr2 is required anyway.
	
	It returns -1 if the variable does not exist or the numeric type of the variable if it does.
*/
int
FetchNumVar(const char *varName, double *doublePtr1, double *doublePtr2)
{
	return(CallBack3(FETCHNUMVAR, (void*)varName, doublePtr1, doublePtr2));
}

/*	StoreNumVar(varName, doublePtr1, doublePtr2)

	StoreNumVar stores the value in double1 and double2 into a variable whose name is specified
	by varName. double1 contains the real part and double2 contains the imaginary part.
	double2 is required even if the variable is not complex.
	
	It returns -1 if the variable does not exist or the numeric type of the variable if it does.
	
	This routine was changed to take double pointer arguments for release 2 of the XOP Toolkit.
	See Igor XOP Tookit Tech Note #1 for details.
*/
int
StoreNumVar(const char *varName, double *doublePtr1, double *doublePtr2)
{
	return(CallBack3(STORENUMVAR, (void*)varName, doublePtr1, doublePtr2));
}

/*	FetchStrVar(varName, stringPtr)

	FetchStrVar returns the value of a named string variable via the pointer stringPtr.
	It returns 0 if it was able to fetch the string or an error code if it was not able.
	
	stringPtr should be big enough to hold a 255 byte string.
*/
int
FetchStrVar(const char *varName, char *stringPtr)
{
	return(CallBack2(FETCHSTRVAR, (void*)varName, stringPtr));
}

/*	FetchStrHandle(varName)

	FetchStrHandle returns the handle containing the text for the named string
	variable or NIL if no such string variable exists.
	
	The text is not null terminated. Use GetHandleSize to determine the
	number of characters in the string.
	
	You should not dispose of or otherwise modify this handle since it belongs
	to Igor.
*/
Handle
FetchStrHandle(const char *varName)
{
	return((Handle)CallBack1(FETCHSTRHANDLE, (void*)varName));
}

/*	StoreStrVar(varName, stringPtr)

	StoreStrVar stores the value in stringPtr into a string variable whose name is specified by
	varName.
	It returns 0 if it was able to store the string or an error code if it was not able.
*/
int
StoreStrVar(const char *varName, const char *stringPtr)
{
	return(CallBack2(STORESTRVAR, (void*)varName, (void*)stringPtr));
}

/*	Variable(varName, varType)

	Variable creates an Igor variable with the specified name and type.
	varType is
		0 for a string variable
		NT_FP32 for a single precision floating point variable (obsolete)
		NT_FP64 for a double precision floating point variable
		
		To make a numeric variable complex, use one of the following:
		(NT_FP32 | NT_CMPLX) (obsolete)
		(NT_FP64 | NT_CMPLX)
	
	HR, 9/29/93; Igor Pro (2.0) addition:
		Use VAR_GLOBAL in type to force the variable or string to be global.
		Example: (NT_FP64 | VAR_GLOBAL)
		
		If you don't use VAR_GLOBAL then the variable or string will be
		global if executed from the command line or local if executed from
		a macro.
	
	HR, 9/28/95; Igor Pro 3.0 change:
		In Igor Pro 3.0, all numeric variables are double-precision. If you pass
		NT_FP32 for varType, Igor will create a double-precision variable anyway.
	
	It returns 0 if it was able to make the variable or an error code
	if it was not able.
*/
int
Variable(const char *varName, int varType)
{
	if (igorVersion < 200 && (varType&VAR_GLOBAL))
		return IGOR_OBSOLETE;		// Requires Igor 2.00 or later.
	return(CallBack2(VARIABLE, (void*)varName, (void *)varType));
}

int
VariableList(Handle listHandle, const char *match, const char *sep, long varTypeCode)
{
	return(CallBack4(VARIABLELIST, listHandle, (void*)match, (void*)sep, (void*)varTypeCode));
}

int
StringList(Handle listHandle, const char *match, const char *sep)
{
	return(CallBack3(STRINGLIST, listHandle, (void*)match, (void*)sep));
}

/*	SetIgorIntVar(numVarName, value, forceGlobal)

	If the named numeric variable already exists, just set it.
	
	If it does not exist, it creates it and then sets it.
	In this case, it will create a global variable if forceGlobal is non-zero.
	If forceGlobal is zero, the variable will be local if a macro is running
	or global if not.
	
	The forceGlobal feature is available in Igor 2.0 or later only.
	You should not call this with forceGlobal!=0 if igorVersion<200.
	
	Returns 0 or error code.
*/
int
SetIgorIntVar(const char* numVarName, long value, int forceGlobal)
{
	int varType;
	double d1, d2;
	int result = 0;

	if (forceGlobal && igorVersion<200)
		return IGOR_OBSOLETE;		// forceGlobal is available in Igor 2.0 or later.

	d1 = value;
	d2 = 0.0;
	if (StoreNumVar(numVarName, &d1, &d2) == -1) {			// Does it exist ?
		varType = NT_FP64;									// HR, 10/22/95: Was NT_FP32.
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))			// Create new variable.
			return result;
		StoreNumVar(numVarName, &d1, &d2);
	}
	
	return result;
}

/*	SetIgorFloatingVar(numVarName, valuePtr, int forceGlobal)

	If the named numeric variable already exists, just set it.
	
	If it does not exist, it creates it and then sets it.
	In this case, it will create a global variable if forceGlobal is non-zero.
	If forceGlobal is zero, the variable will be local if a macro is running
	or global if not.
	
	The forceGlobal feature is available in Igor 2.0 or later only.
	You should not call this with forceGlobal!=0 if igorVersion<200.
	
	Returns 0 or error code.
*/
int
SetIgorFloatingVar(const char* numVarName, double* valuePtr, int forceGlobal)
{
	int varType;
	double d1, d2;
	int result = 0;

	if (forceGlobal && igorVersion<200)
		return IGOR_OBSOLETE;		// forceGlobal is available in Igor 2.0 or later.

	d1 = *valuePtr;
	d2 = 0.0;
	if (StoreNumVar(numVarName, &d1, &d2) == -1) {			// Does it exist ?
		varType = NT_FP64;
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))			// Create new variable.
			return result;
		StoreNumVar(numVarName, &d1, &d2);
	}
	
	return result;
}

/*	SetIgorComplexVar(numVarName, realValuePtr, imagValuePtr, int forceGlobal)

	If the named numeric variable already exists, just set it.
	
	If it does not exist, it creates it and then sets it.
	In this case, it will create a global variable if forceGlobal is non-zero.
	If forceGlobal is zero, the variable will be local if a macro is running
	or global if not.
	
	The forceGlobal feature is available in Igor 2.0 or later only.
	You should not call this with forceGlobal!=0 if igorVersion<200.
	
	Returns 0 or error code.
*/
int
SetIgorComplexVar(const char* numVarName, double* realValuePtr, double* imagValuePtr, int forceGlobal)
{
	int varType;
	int result = 0;

	if (forceGlobal && igorVersion<200)
		return IGOR_OBSOLETE;		// forceGlobal is available in Igor 2.0 or later.

	if (StoreNumVar(numVarName, realValuePtr, imagValuePtr) == -1) {	// Does it exist ?
		varType = NT_FP64 | NT_CMPLX;
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(numVarName, varType))						// Create new variable.
			return result;
		StoreNumVar(numVarName, realValuePtr, imagValuePtr);
	}
	
	return result;
}

/*	SetIgorStringVar(stringVarName, stringVarValue, forceGlobal)

	If the named string variable already exists, just set it.
	
	If it does not exist, it creates it and then sets it.
	In this case, it will create a global variable if forceGlobal is non-zero.
	If forceGlobal is zero, the variable will be local if a macro is running
	or global if not.
	
	The forceGlobal feature is available in Igor 2.0 or later only.
	You should not call this with forceGlobal!=0 if igorVersion<200.
	
	Returns 0 or error code.
*/
int
SetIgorStringVar(const char* stringVarName, const char* stringVarValue, int forceGlobal)
{
	int varType;
	int result;
	
	if (forceGlobal && igorVersion<200)
		return IGOR_OBSOLETE;		// forceGlobal is available in Igor 2.0 or later.

	if (result = StoreStrVar(stringVarName, stringVarValue)) {	// Error if string does not exist.
		varType = 0;											// 0 means string variable.
		if (forceGlobal)
			varType |= VAR_GLOBAL;
		if (result = Variable(stringVarName, varType))
			return result;
		result = StoreStrVar(stringVarName, stringVarValue);
	}
	return result;
}


// *** Command Parsing Routines ***
// These routines for getting parameters from the command line are note needed
// if you use Operation Handler, which is recommended.

/*	IsStringExpression(assumeString)

	Returns the truth that the next item in the command currently being parsed
	is a string expression (e.g. "hello", myStrVar, myStrVar+"hello", etc.).
	
	Use this when the next thing might be a string expression or some other item,
	such as a number or a wave name. This routine is rarely needed since it is usually
	possible to design your operation's syntax so that the type of the next thing
	is pre-determined.
	
	You should normally pass 1 for the assumeString parameter. Igor uses it only
	in a rare case in which the next item in the command line is ambiguous. You
	don't need to understand the following explanation.
	
	Igor uses the assumeString parameter if the next item is something like:
		$varName
	where varName contains the name of yet another variable which might be a
	numeric variable or might be a string variable. In this case, IsStringExpression
	returns whatever value you pass for assumeString.
	
	Returns one of the following:
		0				The next item is not a string expression
	   -1				The next item is a string expression
		IGOR_OBSOLETE	The version of Igor that is running is too old.
	
	If IsStringExpression returns -1, use GetAString or GetAStringInHandle to
	get the value of the string expression. If it returns 0, use NextSymb to find
	what the next symbol is.

	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
*/
int
IsStringExpression(int assumeString)
{
	return(CallBack1(IS_STRING_EXPRESSION, (void*)assumeString));
}

/*	NextSymb()

	NextSymb() reads the next symbol from the command line but does not skip over it.
	This is useful for testing what the next thing is.
	
	The result can be:
		ASCII code for single character found in command line
		two packed ASCII codes for double symbol like +=
		0 meaning end of command line
		-1 meaning there was a name in command line
		-2 meaning there was a number in the command line
*/
int
NextSymb(void)
{
	return(CallBack0(NEXTSYMB));
}

/*	GetSymb()

	GetSymb() reads the next symbol from the command line.
	
	The result can be:
		ASCII code for single character found in command line
		two packed ASCII codes for double symbol like +=
		0 meaning end of command line
		-1 meaning there was a name in command line
		-2 meaning there was a number in the command line
*/
int
GetSymb(void)
{
	return(CallBack0(GETSYMB));
}

/*	GetFlag(flagsPtr)

	GetFlag() tries to read a flag from the command line.
	
	flagsPtr is a pointer to a C string containing the flags that the calling command handles.
	
	Returns the flag (ASCII code) if a flag was gotten.
	Returns 0 if there are no more flags to get.
	Returns -1 if there was a flag in command line but it is not a flag we're looking for. 
*/
int
GetFlag(char *flagsPtr)
{
	int result;
	
	result = CallBack1(GETFLAG, flagsPtr);
	if (result > 0)
		result = flagsPtr[result-1];		// Convert from offset into string to character.
	return(result);
}

/*	GetName(namePtr)

	GetName tries to read a name (of anything) from the command line.
	
	namePtr is a pointer to a C string to receive the the name from the command line.
	string should be long enough to hold MAX_OBJ_NAME+1 bytes.
	
	Returns
		0 if it was able to get a valid name
		; if reached end of command line without getting a name
	or	an error code.
	
	8/10/93
		Changed internally in Igor Pro 2.0 so that $<str> where <str> evaluates
		to "" returns 0 (no error). This allows syntax like /N=$<str> where,
		if <str> == "", it behaves as if there were no /N at all.
*/
int
GetName(char *namePtr)
{
	return(CallBack1(GETNAME, namePtr));
}

/*	GetDataFolderAndName(dataFolderHandlePtr, namePtr, beLiberal)

	GetDataFolderAndName() tries to read a name of any data object from the command line
	and returns the both the name and a handle to the data folder referenced by
	the command.
	
	Use this for command line operations that need to get the name of numeric variable
	or string variable. This will allow you to operate on objects that are not in the
	current data folder.
	
	For getting a handle to an existing wave, you can use GetWave. However
	if the wave may or may not yet exist, you cannot use GetWave. Use this routine
	instead.
	
	For example, if your operation creates an output wave, you may want to let
	the user specify the name and location (parent data folder) of the output wave.
	The alternative is to use a fixed name and always create the wave in the
	current data folder.

	You can also use this routine to get the data folder and name of an existing wave.
	However, if all you need is a handle to the wave, GetWave will suffice.
	
	For getting a handle to an existing data folder, you can use GetDataFolder. However
	if the data folder may or may not yet exist, you cannot use GetDataFolder. Use this
	routine instead.
	
	namePtr is a pointer to a C string to receive the the name from the command line.
	string should be long enough to hold MAX_OBJ_NAME+1 bytes. The name will be
	unquoted.
	
	If beLiberal is 1, it will accept "liberal" names (see CleanupName for a
	discussion of liberal names). If beLiberal is 0, it will return an error
	if the name is liberal. You should pass 1 for beLiberal if the name that you
	are getting is a data folder name or wave name and if your XOP is capable
	of dealing with liberal names (which it should be). Otherwise, pass 0.
	
	Returns
		0 if it was able to get a valid name
	or	an error code.
	
	If the command ended (CR, semicolon or comment symbol encountered) without a name,
	the error code will be CMD_ENDED_WITHOUT_NAME. You can test for this code if
	you want the name parameter at the end of the command to be optional, as in the
	case of a list of waves names.
		
	Added in Igor Pro 3.0. Prior to Igor Pro 3.0, there were no data folders.
	If you call this when running with an earlier version, it will return
	NIL in *dataFolderHandlePtr and then call the old routine, GetName, which
	will return a name via namePtr.
*/
int
GetDataFolderAndName(DataFolderHandle* dataFolderHandlePtr, char *namePtr, int beLiberal)
{
	if (igorVersion < 300) {
		int result;
		*dataFolderHandlePtr = NIL;
		result = GetName(namePtr);
		if (result == ';')
			result = CMD_ENDED_WITHOUT_NAME;
		return result;
	}
	return(CallBack3(GET_DATAFOLDER_AND_NAME, dataFolderHandlePtr, namePtr, (void*)beLiberal));
}

/*	GetDataFolder(dataFolderHandlePtr)

	Returns via *dataFolderHandlePtr a handle to the data folder referenced
	in the command currently being parsed.
	
	Returns 0 or an error code.
	
	Call this routines when the next item in the command should be the name
	of an existing data folder. If the item is the name of a data folder that
	does not or may not yet exist, call GetDataFolderAndName instead.	
	
	This routine handles any kind of reference to a data folder (a full data
	folder path, a partial data folder path, a simple data folder name, just a
	colon for the current data folder).
		
	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
*/
int
GetDataFolder(DataFolderHandle* dataFolderHandlePtr)
{
	return(CallBack1(GET_DATAFOLDER, dataFolderHandlePtr));
}

/*	Keyword(keywords, keyword)

	Keyword checks to see if the string keyword is in the list of array pointed to by
	keywords. If so, it returns the index into the array of strings. If not it returns -1.
	
	keywords is a pointer to the start of an array of strings pointers.
	Each pointer in the array points to an acceptable keyword except that the
	last pointer points to a null string.
	
	This routine is used by GetKeyword. You are more likely to want to use
	GetKeyword than Keyword.
*/
int
Keyword(const char *keywords[], char *keyword)
{
	const char* keywordPtr;
	int i;
	
	Capitalize(keyword);							// Keywords are case insensitive.
	for (i = 0; ;i++) {
		keywordPtr = keywords[i];
		if (*keywordPtr == '\0')					// Found end of list of keywords ?
			return(-1);
		if (strcmp(keyword, keywordPtr) == 0) {
			return(i);
		}
	}
}

/*	GetKeyword(keywords, keyPtr)

	GetKeyword() tries to read a keyword from the command line.
	
	keywords is a pointer to the start of an array of strings pointers.
	each pointer in the array points to an acceptable keyword.
	
	If a keyword is found, GetKeyword returns the index into the array via keyPtr.
	That is, if the command line input matches the first string in the keywords array,
	GetKeyword returns 0 via keyPtr. If it matches the second string it returns 1 via keyPtr
	and so on.
	
	The function result is as follows:
		Returns 0 if it was able to get a valid keyword.
		Returns ; if reached end of command line.
		Returns UNKNOWN_KEYWORD otherwise.
*/
int
GetKeyword(const char *keywords[], int *keyPtr)
{
	char name[MAX_OBJ_NAME+1];
	int result;
	
	result = CallBack1(GETNAME, name);
	if (result) {
		if (result != ';')
			result = UNKNOWN_KEYWORD;
		return(result);
	}
	
	// Got a name, see if it's a keyword.
	result = Keyword(keywords, name);
	if (result < 0)									// No such keyword ?
		return(UNKNOWN_KEYWORD);
	*keyPtr = result;								// Got valid keyword.
	return(0);
}

/*	GetWaveName(namePtr)

	GetWaveName() tries to read the name of an existing wave from the command line.
	
	namePtr is a pointer to a C string to receive the the name from the command line.
	string should be long enough to hold MAX_OBJ_NAME+1 bytes.
	
	Returns 0 if it was able to get a valid wave name.
	Returns NOWAV if got something but it wasn't a wave name.
*/
int
GetWaveName(char *namePtr)
{
	int result;
	
	result = GetName(namePtr);				// Try to get a name.
	if (result)
		return(NOWAV);
	if (FetchWave(namePtr))					// Check if it's name of a wave.
		return(0);
	else
		return(NOWAV);
}

/*	GetWave()

	GetWave() tries to read a wave name from the command line and returns a handle to the
	wave or NIL if there was no wave name.
*/
waveHndl
GetWave(void)
{
	return((waveHndl)CallBack0(GETWAVE));
}

/*	GetWaveList(waves, numWavesPtr)

	GetWaveList() tries to read a list of waves from the command line.
	
	waves is a pointer to an array of wave handles (holds MAXWAVES wave handles)
	numWavesPtr is a pointer to a long to hold the number of waves in list.
	
	Returns 0 if it was able to get one or more waves.
	Returns error code if error getting wavelist.
*/
int
GetWaveList(waveHndl *waves, long *numWavesPtr)
{
	return(CallBack2(GETWAVELIST, waves, numWavesPtr));
}

/*	GetWaveRange(wrp)

	Call to parse /R=(x1, x2) or /R=[p1, p2] type flag after getting the /R= part.
	wrp is a pointer to a WaveRangeRec.
	See declaration of WaveRangeRec structure for details.
	
	Returns error code if error getting range or 0.
*/
int
GetWaveRange(WaveRangeRecPtr wrp)
{
	return(CallBack1(GETWAVERANGE, wrp));
}

/*	CalcWaveRange(wrp)

	Call after getting /R= flag to convert the specified range to point numbers
	relative to a given wave.
	wrp is a pointer to a WaveRangeRec.
	See declaration of WaveRangeRec structure for details.
	
	Returns error code if error calculating range or 0.
*/
int
CalcWaveRange(WaveRangeRecPtr wrp)
{
	return(CallBack1(CALCWAVERANGE, wrp));
}

/*	GetNumVarName(namePtr)

	GetNumVarName() tries to read a numeric variable name from the command line.
	
	namePtr is a pointer to a C string to receive the the name from the command line.
	string should be long enough to hold MAX_OBJ_NAME+1 bytes.
	
	It returns -1 if it did not get the name of a numeric variable or the numeric type of the
	variable if it did.
*/
int
GetNumVarName(char *namePtr)
{
	int result;
	double d1, d2;
	
	result = GetName(namePtr);				// Try to get a name.
	if (result)
		return(-1);
	return(FetchNumVar(namePtr, &d1, &d2));
}

/*	GetStrVarName(namePtr)

	GetStrVarName() tries to read a string variable name from the command line.
	
	namePtr is a pointer to a C string to receive the the name from the command line.
	string should be long enough to hold MAX_OBJ_NAME+1 bytes.
	
	It returns 0 if it did get the name of a string variable or an error code if not.
*/
int
GetStrVarName(char *namePtr)
{
	int result;
	
	result = GetName(namePtr);				// Try to get a name.
	if (result)
		return(EXPECTED_STRINGVARNAME);
	return(FetchStrHandle(namePtr) ? 0:EXPECTED_STRINGVARNAME);
}

/*	GetNum(doublePtr)

	GetNum() tries to read a number from the command line in to the double pointed to by doublePtr.
	The number can also be a numeric expression.
	
	Returns 0 if it was able to get a valid number.
	Returns error code if not able.
*/
int
GetNum(double *doublePtr)
{
	return(CallBack1(GETNUM, doublePtr));
}

/*	GetNum2(doublePtr, terminator)

	GetNum2() tries to read a number from the command line in to the double pointed to by doublePtr.
	The number can also be a numeric expression.
	
	terminator is a code identifying the characters allowed to end the number of expression.
	See IgorXOP.h for codes.
	
	Returns 0 if it was able to get a valid number.
	Returns error code if not able.
*/
int
GetNum2(double *doublePtr, int terminator)
{
	return(CallBack2(GETNUM2, doublePtr, (void *)terminator));
}

/*	GetFlagNum(doublePtr)

	GetFlagNum() tries to read a number from the command line in to the double pointed to by
	doublePtr. The number can also be a numeric expression.
	
	This is like GetNum but should be used when getting a number associated with a flag, 
	as in XXX /N=32 .... (used to get the 32).
	
	Unlike most command parsing routines, GetFlagNum does not automatically skip a leading
	comma. This is generally what you would want for this type of syntax.

	Returns 0 if it was able to get a valid number.
	Returns error code if not able.
*/
int
GetFlagNum(double *doublePtr)
{
	return(CallBack1(GETFLAGNUM, doublePtr));
}

/*	GetTrueOrFalseFlag(flagMask, flagsPtr)

	For parsing a flag of the form
		/X=x
	where the =x part is optional.
	
	flagMask is a bit mask with one bit set, the bit that represents /X.
	This routine is called after you have parsed the /X part.
	It sets the specified bit of *flagsPtr. Then, it checks to see if
	the =x optional part is there. If it is, it parses it and sets
	or clears the bit in *flagsPtr.
	
	Returns 0 if OK or an error code.
	
	Added 10/14/93.
*/
int
GetTrueOrFalseFlag(long flagMask, long* flagsPtr)
{
	double d1;
	int result;
	
	*flagsPtr |= flagMask;
	if (NextSymb() == '=') {				// Have /X=<boolean> ?
		GetSymb();							// Suck up =.
		if (result = GetFlagNum(&d1))
			return(result);
		if (d1 == 0)
			*flagsPtr &= ~flagMask;
	}
	return 0;
}

/*	GetLong(longPtr)

	GetLong() tries to read a number from the command line in to the long pointed to by longPtr.
	The number can also be a numeric expression.
	
	Returns 0 if it was able to get a valid number.
	Returns error code if not able.
*/
int
GetLong(long *longPtr)
{
	return(CallBack1(GETLONG, longPtr));
}

/*	GetAString(stringPtr)

	GetAString() tries to read a string literal or string expression
	from the command line.

	stringPtr is a pointer to a C string to receive the the string
	from the command line.
	It should be able to handle 256 bytes.
	
	Returns 0 if it was able to get a valid string.
	Returns error code if not able -- string is set to "".
*/
int
GetAString(char *stringPtr)
{
	*stringPtr = '\0';				// In case can't get string.
	return(CallBack1(GETSTRING, stringPtr));
}

/*	GetAStringInHandle(void)

	GetAStringInHandle() tries to read a string literal or string expression
	from the command line.

	Returns a handle if it was able to get a valid string.
	Returns NIL if there was an error in getting the string.
	The error will have already been reported to the user.
	
	You should use this if the string you expect may be greater than 255 characters.
	Otherwise, use GetAString, above.
	
	The handle belongs to the calling routine. You should dispose of
	it when you are finished with it.
	
	The handle contains just the characters in the string, no count
	byte or trailing null. Use GetHandleSize to find count.
	
	Unlike most command parsing routines including GetString, GetAStringInHandle
	does not automatically skip a leading comma. If you want to require or tolerate
	a comma before the string, use GetSymb().
*/
Handle
GetAStringInHandle(void)
{
	return((Handle)CallBack0(GETSTRINGINHANDLE));
}

/*	CheckTerm()

	CheckTerm returns 0 if there is no more stuff in the current command or an error code
	otherwise. It should be called after all XOP parameters are satisfied to make sure that
	there's no extra garbage.
	
	Call CheckTerm when you have finished parsing all items in a command with a fixed
	number of parameters.
*/
int
CheckTerm(void)
{
	return CallBack0(CHECKTERM);
}

/*	AtEndOfCommand()

	AtEndOfCommand returns 1 if there is no more stuff in the current command, meaning
	that you should stop parsing, or 0 if there is more stuff to be parsed. Basically,
	it checks for a semicolon, comment symbol, or carriage return, marking the end
	of a command.
	
	Call AtEndOfCommand to see if you have finished parsing all items in a command with a
	variable number of parameters.
	
	Added in Igor Pro 4.00. If you call this when running with an earlier version,
	it will attempt to simulate the desired behavior. 
*/
int
AtEndOfCommand(void)
{
	if (igorVersion < 400) {
		int symb;
		symb = NextSymb();
		return (symb==0 || symb==CR_CHAR || symb==';' || symb=='|');
	}
	return CallBack0(AT_END_OF_COMMAND);
}

/*	UniqueName(baseName, finalName)

	Given a base name (like "wave") UniqueName returns a name (like "wave3")
	via finalName that does not conflict with any existing names.
	
	Returns the number used to make the name unique.
*/
int
UniqueName(const char *baseName, char *finalName)
{
	return(CallBack2(UNIQUENAME, (void*)baseName, finalName));
}

/*	UniqueName2(nameSpaceCode, baseName, finalName, suffixNumPtr)

	Given a base name (like "wave") UniqueName2 returns a name (like "wave3")
	via finalName that does not conflict with any existing names. The number
	appended to make the name unique will be *suffixNumPtr or greater.
	Igor sets *suffixNumPtr to the number Igor used to make the name unique.
	
	nameSpaceCode is:
		MAIN_NAME_SPACE			for Igor's main name space (waves, variables, windows)
		DATAFOLDER_NAME_SPACE
		See IgorXOP.h for other less frequently-used name space codes.
	
	Returns 0 if OK, -1 for a bad nameSpaceCode or some other non-zero error code.
	
	Support for DATAFOLDER_NAME_SPACE was added in Igor Pro 3.0. Previously, Igor did
	not support data folders. If you pass DATAFOLDER_NAME_SPACE to UniqueName2 when
	running with a pre-3.0 version of Igor, you will get a -1 error code.
	
	You should use suffixNumPtr as follows:
		long suffixNum = 0;
		for(i = 0; i < numObjectsToBeCreated; i++) {
			if (err = UniqueName2(nameSpaceCode, baseName, finalName, &suffixNum))
				break;
			MakeObject(finalName);
		}
*/
int
UniqueName2(int nameSpaceCode, const char *baseName, char *finalName, long* suffixNumPtr)
{
	return(CallBack4(UNIQUENAME2, (void*)nameSpaceCode, (void*)baseName, finalName, suffixNumPtr));
}

/*	SanitizeWaveName(waveName, column)

	NOTE: If your XOP requires IGOR Pro 3 or later, you should use CleanupName
		  instead of the older SanitizeWaveName.

	Given a pointer to a C string containing a proposed wave name,
	SanitizeWaveName() changes it to make it a valid wave name if necessary.
	
	Returns 1 if it had to make a change, 0 if name was OK to begin with.
	
	First, it chops the string off if it is too long.
	Then, it makes sure that the first character is alphabetic.
	Then it replaces any subsequent characters that are not alphanumeric with underscore.
	
	Added for Igor 2.0D83.
	
	HR, 11/22/95: Added ch>0 check to match what Igor does internally.
*/
int
SanitizeWaveName(char *waveName, long column)
{
	int len, i, ch;
	int result = 0;
	
	len = strlen(waveName);
	if (len==0) {
		sprintf(waveName, "NoName%ld", column);
		return 1;
	}

	ch = waveName[0];
	if (!(ch>0 && isalpha(ch))) {					// First char must be alphabetic.
		memmove(waveName+1, waveName, len+1);		// Shift chars over.
		waveName[0] = 'X';
		len += 1;
		result = 1;
	}
	
	if (len > MAX_WAVE_NAME) {
		waveName[MAX_WAVE_NAME] = '\0';				// Truncate name.
		len = MAX_WAVE_NAME;
		result = 1;
	}
	
	for (i = 1; i < len; i++) {						// Subsequent characters must be.
		ch = waveName[i];							// Alphanumeric or underscore.
		if (!(ch>0 && isalnum(ch))) {
			waveName[i] = '_';
			result = 1;
		}
	}
	
	return result;
}

/*	CheckName(dataFolderH, objectType, name)

	Checks the name for legality and uniqueness.
	
	If dataFolderH is NIL, it looks for conflicts with objects in the current data folder.
	If it is not NIL, it looks for conflicts in the folder specified by dataFolderH.
	
	objectType is one of the following which are defined in IgorXOP.h:
		WAVE_OBJECT
		VAR_OBJECT					(numeric variable)
		STR_OBJECT					(string variable)
		GRAPH_OBJECT
		TABLE_OBJECT
		LAYOUT_OBJECT
		PANEL_OBJECT
		NOTEBOOK_OBJECT
		DATAFOLDER_OBJECT
		PATH_OBJECT
		PICT_OBJECT

	Returns 0 if the name is legal and is not in conflict with an existing object.
	Returns an Igor error code otherwise.
		
	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
CheckName(DataFolderHandle dataFolderH, int objectType, const char* name)
{
	return(CallBack3(CHECKNAME, dataFolderH, (void*)objectType, (void*)name));
}

/*	PossiblyQuoteName(name)
	
	name contains an Igor object name.
	PossiblyQuoteName puts single quotes around the name if they would be
	needed to use the name in Igor's command line.

	Igor Pro 3.0 and later allows wave and data folder names to contain characters,
	such as space and dot, that were previously illegal in names. We call this
	"liberal" name rules.
	
	If an object has such a name, you must single-quote the name to use it in Igor's
	command line. This includes using it in the Execute operation or in the XOPCommand,
	XOPSilentCommand, or FinishDialogCmd XOPSupport routines. Thus, if you are going
	to use a wave or data folder name for this purpose, you should call PossiblyQuoteName
	to add the quotes if needed.
	
	NOTE: name must be able to hold two additional characters. Thus, you
		  should declare name: char name[MAX_OBJ_NAME+2+1];
		  
	NOTE: Liberal rules are still not allowed for string and numeric variable names.

	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	it will do nothing and return 0.
*/
int
PossiblyQuoteName(char *name)
{
	if (igorVersion < 300)
		return 0;									// Quoted names added in Igor Pro 3.0.
	return(CallBack1(POSSIBLY_QUOTE_NAME, name));
}

/*	CatPossiblyQuotedName(char* str, char* name)

	Adds the specified Igor object name to the end of the string.
	If necessary, puts single quotes around the name so that it can be
	used in the Igor command line.
	
	Use this to concatenate a wave name to the end of a command string
	when the wave name may be a liberal name that needs to be quoted to
	be used in the command line. 
	
	Example:
		char waveName[MAX_OBJ_NAME+1];				// This contains a wave name.
		char cmd[256];
		strcpy(cmd, "Redimension/N=1000 ");
		CatPossiblyQuotedName(cmd, waveName);
		XOPSilentCommand(cmd);

	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	it will do the concatenation but with no quoting.
*/
void
CatPossiblyQuotedName(char* str, const char* name)
{
	int len;
	
	len = strlen(str);
	strcpy(str + len, name);
	PossiblyQuoteName(str + len);
}

/*	CleanupName(beLiberal, name, maxNameChars)
	
	name contains an Igor object name.
	CleanupName changes it, if necessary, to make it a legal name.
	
	Returns 0 or IGOR_OBSOLETE if running with a version of Igor prior to 3.0.
	
	For most uses, pass MAX_OBJ_NAME for the maxNameChars parameter.

	Igor Pro 3.0 and later allows wave and data folder names to contain characters,
	such as space and dot, that were previously illegal in names. We call this
	"liberal" name rules.
	
	If beLiberal is non-zero, CleanupName uses liberal name rules. Liberal rules are still not
	allowed for string and numeric variable names so pass zero for beLiberal for these objects.
	
	If you are going to use the name in Igor's command line (via an Execute operation
	or via the XOPCommand or XOPSilentCommand callbacks), and if the name uses liberal
	rules, the name may need to be single-quoted. In these cases, you should call
	PossiblyQuoteName after calling CleanupName.

	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	it will do nothing and return IGOR_OBSOLETE.
*/
int
CleanupName(int beLiberal, char *name, long maxNameChars)
{
	if (igorVersion < 300)
		return IGOR_OBSOLETE;						// This call added in Igor Pro 3.0.
	return(CallBack3(CLEANUP_NAME, (void*)beLiberal, name, (void*)maxNameChars));
}

/*	CreateValidDataObjectName(...)
	
	This routine is designed to do all of the nasty work needed to get
	a legal name for a given object. It cleans up illegal names and resolves
	name conflicts.

	It returns in outName a name that can be safely used to create a data object
	of a given type (wave, string, variable, numeric variable).
	
	inName is the proposed name for the object.

	outName is the name after possible cleanup and uniquification.
	
	suffixNumPtr is used to speed up the search for unique names.
	See UniqueName2 for details.
	
	dataFolderH is a handle to a data folder or NIL to use the current data folder.
	
	objectType is one of the following:
		WAVE_OBJECT, VAR_OBJECT, STR_OBJECT, DATAFOLDER_OBJECT.

	beLiberal is 1 if you want to allow liberal names or 0 if not.

	allowOverwrite is 1 if it is OK to for outName to be the name of an
	existing object of the same type.
	
	inNameIsBaseName is 1 if inName is a base name (e.g., "wave") to which a
	suffix (e.g., "0") must always be added to produce the actual name (e.g., "wave0").
	If inNameIsBaseName is 0, then no suffix will be added to inName unless it
	is needed to make the name unique.
	
	printMessage is 1 if you want CreateValidDataObjectName to print a message
	in Igor's history area if an unexpected name conflict occurs. A message is
	printed if you are not using a base name and not allowing overwriting and
	there is a name conflict. A message is also printed if a conflict with
	an object of a different type prevents the normal name from being used.
	
	CreateValidDataObjectName sets *nameChangedPtr to the truth that outName
	is different from inName.
	
	It sets *doOverwritePtr to 1 if outName is the name of an existing object
	of the same type and allowOverwrite is 1.
	
	inName and outName can point to the same array if you don't want to
	preserve the original name. Both must be big enough to hold MAX_OBJ_NAME+1 bytes.
	
	If the object type is VAR_OBJECT or STR_OBJECT, the name will not be
	liberal, even if beLiberal is 1. Igor allows only wave and data folder
	names to be liberal.
	
	Returns 0 if OK or an error code.
		
	Added in Igor Pro 3.0. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
CreateValidDataObjectName(
	DataFolderHandle dataFolderH,
	const char* inName, char* outName, long* suffixNumPtr, int objectType,
	int beLiberal, int allowOverwrite, int inNameIsBaseName, int printMessage,
	int* nameChangedPtr, int* doOverwritePtr)
{
	DataFolderHandle oldCurrentDataFolderH;
	char name[MAX_OBJ_NAME+1];
	char inName2[MAX_OBJ_NAME+1];
	int needToUniquify;
	int objectType2;
	char temp[256];
	int needMessage;
	int err;
	
	err = 0;
	*nameChangedPtr = *doOverwritePtr = 0;
	needMessage = 0;
	
	/*	HR, 030701: Set current data folder based on dataFolderH because some of the calls
		below (e.g., UniqueName2), are not data-folder-aware.
	*/
	oldCurrentDataFolderH = NULL;
	if (dataFolderH != NULL) {
		if (GetCurrentDataFolder(&oldCurrentDataFolderH) != 0)
			oldCurrentDataFolderH = NULL;					// Should never happen.
		if (err = SetCurrentDataFolder(dataFolderH))
			return err;										// dataFolderH is invalid?
	}

	// HR, 050110, 5.04: Avoid crash if inName is too long.
	if (inNameIsBaseName)
		strncpy(inName2, inName, MAX_OBJ_NAME-5);
	else
		strncpy(inName2, inName, MAX_OBJ_NAME);
	inName2[MAX_OBJ_NAME] = 0;								// Needed to guarantee null terminator if inName >= MAX_OBJ_NAME characters. 

	strcpy(name, inName2);	// Assume that inName2 is the desired name, not a base name to which  we must add a suffix.
	
	// If inName2 is a base name to which we must add a suffix, add the suffix here.
	if (inNameIsBaseName) {
		if (allowOverwrite) {
			sprintf(name, "%s%ld", inName2, *suffixNumPtr);				// Use next suffix even if it is already used.
			*suffixNumPtr += 1;
		}
		else {
			UniqueName2(MAIN_NAME_SPACE, inName2, name, suffixNumPtr);	// Choose a suffix that is not already used.
		}
	}

	if (objectType!=WAVE_OBJECT && objectType!=DATAFOLDER_OBJECT)
		beLiberal = 0;										// Liberal names not supported for variables.
	if (err = CleanupName(beLiberal, name, MAX_OBJ_NAME))	// Remove illegal characters, etc.
		goto done;		// This will be IGOR_OBSOLETE if we are running with pre-3.0 Igor.
	
	needToUniquify = 0;
	do {
		if (GetDataFolderObject(dataFolderH, name, &objectType2, NIL) == 0) {
			// If here, an object exists with the specified name.
			
			if (allowOverwrite) {
				if (objectType2 == objectType) {
					*doOverwritePtr = 1;
					break;							// OK to overwrite object.
				}
			}
			
			/*	If here, we must choose another name because the name is in use
				for a different type of object or it is in use for a the same type but
				we are not allowed to overwrite it.
			*/
			needToUniquify = 1;
			needMessage = printMessage;
		}
		else {
			/*	Name is not a name of an existing object. Make sure it is not in conflict
				with a function, operation or some other reserved name.
			*/
			if (CheckName(dataFolderH, objectType, name)) {
				needToUniquify = 1;					// There is some kind of conflict.
				needMessage = printMessage;
			}
		}

		if (!needToUniquify)
			break;

		strcpy(temp, name);
		if (err = UniqueName2(MAIN_NAME_SPACE, temp, name, suffixNumPtr))
			goto done;
	} while(0);
	
	*nameChangedPtr = strcmp(inName, name)!=0;
	if (needMessage) {
		char objectTypeStr[64];
		switch(objectType) {
			case WAVE_OBJECT:
				strcpy(objectTypeStr, "Wave");
				break;
			case VAR_OBJECT:
				strcpy(objectTypeStr, "Variable");
				break;
			case STR_OBJECT:
				strcpy(objectTypeStr, "String");
				break;
			case DATAFOLDER_OBJECT:
				strcpy(objectTypeStr, "Data folder");
				break;
			default:
				sprintf(objectTypeStr, "BUG: CreateValidDataObjectName, objectType=%d", objectType);
				break;
		}
		sprintf(temp, "%s \'%s\' changed to \'%s\' because of a conflict.\015", objectTypeStr, inName2, name);
		XOPNotice(temp);
	}
	
	strcpy(outName, name);

done:
	if (oldCurrentDataFolderH != NULL)
		SetCurrentDataFolder(oldCurrentDataFolderH);
	
	return err;
}

/*	GetFormat(format)

	GetFormat is called to get the =<string expression> part of a /F="..." format specifier.
	The only parameter is a pointer to the format string to contain the format.

	Unlike most command parsing routines, GetFormat does not automatically skip a leading
	comma. This is generally what you would want for this type of syntax.
	
	GetFormat returns 0 or an error code.
*/
int
GetFormat(char *format)
{
	return(CallBack1(GETFORMAT, format));	
}

int
DoCHIO(CHIORecPtr CHIOPtr)
{
	return(CallBack1(CHIO, CHIOPtr));
}


// *** IGOR Color Table Routines ***

/*	GetIndexedIgorColorTableName(index, name)

	Returns via name the name of a color table indicated by the index
	or "" if the index is invalid. Valid indices start from zero. You can find
	the maximum valid index by calling this routine with increasing indices
	until it returns an error.

	The function result is 0 if OK or a non-zero error code.
	
	Added in Igor Pro 3.1. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
GetIndexedIgorColorTableName(int index, char name[MAX_OBJ_NAME+1])
{
	*name = 0;
	
	return CallBack2(GET_INDEXED_IGORCOLORTABLE_NAME, (void*)index, name);
}

/*	GetNamedIgorColorTableHandle(name, ictHPtr)
	
	Returns via *ictHPtr a handle to an IGOR color table or NULL in case of error.
	
	The returned handle belongs to Igor. Do not modify or delete it.
	
	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	The name parameter is case insensitive.
	
	As of this writing, the following IGOR color tables exist:
		Rainbow, Grays, YellowHot, BlueHot
		BlueRedGreen, RedWhiteBlue, PlanetEarth, Terrain
	You can find the names of all color tables using GetIndexedIgorColorTableName.
	
	The function result is 0 if OK or a non-zero error code.
	
	Added in Igor Pro 3.1. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
GetNamedIgorColorTableHandle(const char *name, IgorColorTableHandle* ictHPtr)
{
	*ictHPtr = NULL;

	return CallBack2(GET_NAMED_IGORCOLORTABLE_HANDLE, (void*)name, ictHPtr);
}

/*	GetIgorColorTableInfo(ictH, name, numColorsPtr)

	Provides access to the name and the number of colors in the IGOR
	color table specified by ictH.

	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	If you don't want to know the name, pass NULL for name.
	If you don't want to know the number of colors, pass NULL for numColorsPtr.
	
	The function result is 0 if OK or a non-zero error code.
	
	Added in Igor Pro 3.1. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
GetIgorColorTableInfo(IgorColorTableHandle ictH, char name[MAX_OBJ_NAME+1], int* numColorsPtr)
{
	// In case IGOR_OBSOLETE.
	if (name != NULL)
		*name = 0;
	if (numColorsPtr != NULL)
		*numColorsPtr = 0;
	
	return CallBack3(GET_IGORCOLORTABLE_INFO, ictH, name, numColorsPtr);
}

/*	GetIgorColorTableValues(ictH, startColorIndex, endColorIndex, updatePixelValues, csPtr)
	
	Returns via csPtr a description of the colors associated with the IGOR color
	table specified by ictH.

	The IgorColorTableHandle type is defined in IgorXOP.h.
	
	startColorIndex and endColorIndex specify the indices of the colors for
	which you want to get a description. startColorIndex must be between 0 and
	the number of colors in the table minus one. endColorIndex must be between
	startColorIndex and the number of colors in the table minus one. You can
	find the number of colors in the table using GetIgorColorTableInfo.
	
	The IgorColorSpec structure contains an RGBColor field which identifies the
	RGB color for a color table entry with a particular index.
	
	The value field of the IgorColorSpec structure tells you the pixel value that
	would need to be written to video RAM to get the associated color to appear on
	the screen when the monitor is in 16 or 256 color mode. It is typically used by
	advanced programmers who are writing directly to offscreen bitmap memory. 
	
	However, when a monitor is running in 16 or 256 color mode, this value is
	invalidated whenever the system changes the hardware color lookup table, which
	can happen at any time. If you pass non-zero for the updateVals parameter, then
	Igor will update the value field for each color before returning it to you and
	it will be accurate until the next time the system changes the hardware color
	lookup table. If you pass zero for the updateVals parameter, then Igor will not
	update the value field and it is likely to be stale. 
	
	Updating the value fields takes time so you should pass non-zero for the
	updateVals parameter only if you really need accurate pixel values. For
	example, if you just want to know what RGB colors appear in a particular color
	table then you don't need the pixel values and should pass 0 for the updateVals
	parameter. On the other hand, if you are writing into an offscreen bitmap in
	preparation for blasting it to the screen, then you need accurate pixel values
	and you should pass 1 for updateVals. 
	
	The function result is 0 if OK or a non-zero error code.
	
	Added in Igor Pro 3.1. If you call this when running with an earlier version,
	you will receive the IGOR_OBSOLETE error code as the function result.
*/
int
GetIgorColorTableValues(IgorColorTableHandle ictH, int startColorIndex, int endColorIndex, int updatePixelValues, IgorColorSpec* csPtr)
{
	return CallBack5(GET_IGORCOLORTABLE_VALUES, ictH, (void*)startColorIndex, (void*)endColorIndex,  (void*)updatePixelValues, csPtr);
}


// *** Cross-Platform Utilities ***

void
WinRectToMacRect(const RECT* wr, Rect* mr)
{
	// In principle, the Windows coordinates could exceed 16 bits. In practice, this should not happen.

	mr->left = wr->left;
	mr->right = wr->right;
	mr->top = wr->top;
	mr->bottom = wr->bottom;
}

void
MacRectToWinRect(const Rect *mr, RECT *wr)
{
	wr->left = mr->left;
	wr->top = mr->top;
	wr->right = mr->right;
	wr->bottom = mr->bottom;
}


// *** Miscellaneous Routines ***

/*	XOPBeep()

	Emits a short beep.
*/
void
XOPBeep(void)
{
	#ifdef _MACINTOSH_
		SysBeep(5);
	#endif
	#ifdef _WINDOWS_
		MessageBeep(MB_ICONEXCLAMATION);
	#endif
}

/*	GetXOPIndString(text, strID, item)

	Tries to get string from a STR# resource in the XOP's resource fork.
	Does not search any other resource forks and does not change current
	resource fork on Macintosh.
	
	text is returned as a C string.
*/
void
GetXOPIndString(
	char* text,									// Should hold up to 256 bytes.
	int strID,									// Resource ID of STR# resource.
	int index)									// String number starting from 1.
{
	*text = 0;					// HR, 981022: For XOP Toolkit 3.1.
	
	#ifdef _MACINTOSH_			// [
	{	int curResFile;
		
		if (GetXOPResource('STR#', strID) == NULL) 	// No such STR# resource?
			return;
	
		curResFile = CurResFile();
		UseResFile(XOPRefNum());					// Search XOP's resource fork.
		GetIndString((unsigned char*)text, strID, index);
		CopyPascalStringToC((ConstStr255Param)text, text);
		UseResFile(curResFile);
	}
	#endif						// ]

	#ifdef _WINDOWS_			// [
		GetWinIndString(XOPModule(), text, strID, index);
	#endif						// ]
}

void
ArrowCursor(void)
{
	CallBack1(SETCURSOR, (void *)ARROWCURSOR);
}

void
IBeamCursor(void)
{
	CallBack1(SETCURSOR, (void *)IBEAMCURSOR);
}

void
WatchCursor(void)
{
	CallBack1(SETCURSOR, (void *)WATCHCURSOR);
}

void
HandCursor(void)
{
	CallBack1(SETCURSOR, (void *)HANDCURSOR);
}

void
SpinCursor(void)
{
	static unsigned long lastTime=0;			// Prevents wasting too much time spinning.
	
	if (TickCount() > lastTime+6) {
		CallBack1(SETCURSOR, (void *)SPINNINGCURSOR);
		lastTime = TickCount();
	}
}

/*	SpinProcess()

	SpinProcess spins the beach ball cursor AND gives MultiFinder (if active) a chance
	to do background processing. Also, Igor may be moved from the background to the
	foreground or vice versa when SpinProcess is called.
	
	It returns non-zero if the user has pressed cmd-dot recently or zero otherwise.
*/
int
SpinProcess(void)
{
	return(CallBack0(SPINPROCESS));
}

/*	DoUpdate()

	Causes Igor to do an immediate update.
	An update consists of
		updating any windows that have been uncovered
		updating any graphs, tables or layouts that need it
		recalculating any dependent objects that need it
	
	The DOUPDATE message was added in Igor 2.0.
	If the version of Igor is earlier than 2.0, DoUpdate triggers
	an update by using XOPCommand. The update is a side effect of
	XOPCommand.
*/
int
DoUpdate(void)
{
	if (igorVersion < 200)
		return XOPCommand("");
	
	CallBack0(DOUPDATE);
	return 0;
}

/*	PauseUpdate(savePtr)

	Used to temporarily suppress updating of host windows.
	savePtr is a pointer to a long to save the previous state of PauseUpdate.
	
	This MUST be balanced by a ResumeUpdate(savePtr).
*/
void
PauseUpdate(long *savePtr)
{
	CallBack1(PAUSEUPDATE, savePtr);
}

/*	ResumeUpdate(savePtr)

	Used to undo temporary suppression of updating of host windows.
	
	This MUST called to balance PauseUpdate(savePtr) calls.
*/
void
ResumeUpdate(long *savePtr)
{
	CallBack1(RESUMEUPDATE, savePtr);
}

int
WaveList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	return(CallBack4(WAVELIST, listHandle, (void*)match, (void*)sep, (void*)options));
}

int
WinList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	return(CallBack4(WINLIST, listHandle, (void*)match, (void*)sep, (void*)options));
}

int
PathList(Handle listHandle, const char *match, const char *sep, const char *options)
{
	return(CallBack4(PATHLIST, listHandle, (void*)match, (void*)sep, (void*)options));
}

/*	GetPathInfo(pathName, vRefNumPtr, dirIDPtr, wdRefNumPtr)

	This routines is implemented for Macintosh only. New XOPs should use
	GetPathInfo2 instead.

	Given the name of an Igor symbolic path in pathName, returns via the other
	parameters the volume reference number, directory ID and working directory
	reference number of the path.
	
	Returns 0 or NO_SYMFLDR if the pathName is not the name of an existing Igor symbolic path.
	
	Added for Igor 1.25.
*/
#ifdef _MACINTOSH_		// [
int
GetPathInfo(const char* pathName, long* vRefNumPtr, long*dirIDPtr, long*wdRefNumPtr)
{
	return(CallBack4(PATHINFO, (void*)pathName, vRefNumPtr, dirIDPtr, wdRefNumPtr));
}
#endif 					// _MACINTOSH_ ]

/*	GetPathInfo2(pathName, fullDirPath)

	pathName is the name of an Igor symbolic path.
	
	Returns via fullDirPath the full native path to the directory referenced by pathName.
	The returned path includes a trailing colon on Macintosh and a trailing backslash on Windows.
	
	Returns 0 if OK or an error code if the pathName is not the name of an existing
	Igor symbolic path.
	
	Added in Igor Pro 3.13. If you call this when running with an earlier version,
	it will return IGOR_OBSOLETE.
*/
int
GetPathInfo2(const char* pathName, char fullDirPath[MAX_PATH_LEN+1])
{
	return CallBack2(PATHINFO2, (void*)pathName, (void*)fullDirPath);
}

/*	GetNamedFIFO

	Returns NamedFIFO handle or NULL if none of that name exists
	NamedFIFO data structure and meaning are documented in a separate file.

	Added for Igor Pro 2.0.
*/
struct NamedFIFO **
GetNamedFIFO(const char *name)
{
	return (void *)CallBack1(GETNAMEDFIFO, (void*)name);
}

/*	MarkFIFOUpdated

	Call this after putting data in a named fifo so chart gadgets will refresh.

	Added for Igor Pro 2.0.
*/
void
MarkFIFOUpdated(struct NamedFIFO **fifo)
{
	CallBack1(MARKFIFOUPDATED, fifo);
}

/*	SaveXOPPrefsHandle(prefsHandle)

	Saves the handle in IGOR's preferences file. You can retrieve the handle
	using GetXOPPrefsHandle.
	
	IGOR makes a copy of the data in the handle, so the handle is still yours
	after you call this. Keep or dispose of it as you wish.
	
	If you pass NULL for the prefsHandle parameter, IGOR removes any existing
	XOP preferences from the IGOR preferences file.
	
	IGOR uses the name of your XOP's file to distinguish your preferences from
	the preferences of other XOPs.
	
	Each time you call this routine, the Igor preferences file is opened and closed.
	Therefore, it is best to call each of it only once. One way to do this is to call
	GetXOPPrefsHandle when your XOPs starts and SaveXOPPrefsHandle when you receive
	the CLEANUP message.
	
	Returns 0 if OK or a non-zero error code.

	Added for Igor Pro 3.1. However, you can call this routine with any
	version of IGOR Pro. If you are running with an IGOR Pro prior to 3.1,
	this routine will save the preferences in your XOP's resource fork
	with a resource type of 'XPRF' and a resource ID of 1100. This applies
	to Macintosh only since there are no IGOR versions prior to 3.1 on Windows.
*/
int
SaveXOPPrefsHandle(Handle prefsHandle)
{
	#ifdef _MACINTOSH_
		if (IgorVersion() < 310) {
			Handle h;
			int curResNum;
			int err;
			
			// Remove resource if it exists.
			h = GetXOPResource('XPRF', 1100);
			if (h != NULL) {
				RemoveResource(h);
				DisposeHandle(h);
			}
	
			// Make a copy of input handle.
			h = prefsHandle;
			err = HandToHand(&h);
			if (err != 0)
				return err;
			
			curResNum = CurResFile();
			UseResFile(XOPRefNum());
			AddResource(h, 'XPRF', 1100, (unsigned char*)"");
			err = ResError();
			UseResFile(curResNum);
		
			return err;
		}
	#endif

	return(CallBack1(SAVE_XOP_PREFS, prefsHandle));
}

/*	GetXOPPrefsHandle(prefsHandlePtr)

	Retrieve your XOP's preference handle from the IGOR preferences file, if you
	have previously stored it there using SaveXOPPrefsHandle. In this case,
	on return, *prefsHandlePtr will be your preferences handle. This handle
	is allocated by IGOR but belongs to you to keep or dispose as you wish.
	
	If the IGOR preferences file does not contain your preferences, on return,
	*prefsHandlePtr will be NULL and the function result will be 0.
	
	IGOR uses the name of your XOP's file to distinguish your preferences from
	the preferences of other XOPs.
	
	Each time you call this routine, the Igor preferences file is opened and closed.
	Therefore, it is best to call each of it only once. One way to do this is to call
	GetXOPPrefsHandle when your XOPs starts and SaveXOPPrefsHandle when you receive
	the CLEANUP message.
	
	Returns 0 if OK or a non-zero error code.

	Added for Igor Pro 3.1. However, you can call this routine with any
	version of IGOR Pro. If you are running with an IGOR Pro prior to 3.1,
	this routine will get the preferences from your XOP's resource fork
	with a resource type of 'XPRF' and a resource ID of 1100. This applies
	to Macintosh only since there are no IGOR versions prior to 3.1 on Windows.
*/
int
GetXOPPrefsHandle(Handle* prefsHandlePtr)
{
	*prefsHandlePtr = NULL;

	#ifdef _MACINTOSH_
		if (IgorVersion() < 310) {
			Handle h;
			int err;
			
			// Get resource if it exists.
			h = GetXOPResource('XPRF', 1100);
			if (h == NULL)
				return 0;								// Resource does not exist.
			
			// Make a copy of handle to give to calling routine.
			err = HandToHand(&h);
			if (err != 0)
				return err;
			
			*prefsHandlePtr = h;
			return 0;
		}
	#endif

	return(CallBack1(GET_XOP_PREFS, prefsHandlePtr));
}

/*	GetPrefsState(prefsStatePtr)

	Returns via bit 0 of prefsStatePtr the truth that preferences are on.
	See the IGOR Pro manual for information about the preferences on/off state.
	Other bits are reserved for future use.

	Function result is 0 if OK or an error code.

	Added for Igor Pro 3.10. The function result will be IGOR_OBSOLETE and
	*prefsStatePtr will be set to 1 if you are running with an older version of IGOR.
*/
int
GetPrefsState(long* prefsStatePtr)
{
	*prefsStatePtr = 1;				// Default if we are running with old IGOR.
	return(CallBack1(GET_PREFS_STATE, prefsStatePtr));
}

/*	XOPDisplayHelpTopic(title, topicStr, flags)

	Displays help for the specified topic.
	
	topicStr is a help topic string that matches a help topic or subtopic in a
	native Igor help file. Igor first searches open help files for the topic.
	If it is not found, Igor then searches all Igor native help files in the folder
	containing the XOP file and subfolders. If it is still not found Igor then
	searches all Igor native help files in the Igor Pro folder and subfolders.
	
	The help file must be compiled in order for Igor to find the topic. Each time
	you open a file as a help file, Igor checks to see if it is compiled and if
	not asks if you want to compile it.
	
	topicStr may have one of the following formats:
		Format							Example
		<topic name>					"GBLoadWave XOP"
		<subtopic name>					"The Load General Binary Dialog"
		<topic name>[<subtopic name>]	"GBLoadWave XOP[The Load General Binary Dialog]"
		
	If the topic that you want to display is a subtopic, you should use the last
	form since it minimizes the chance that Igor will find another help file with
	the same subtopic name. Also, you must choose descriptive topic and
	subtopic names to minimize the chance of a conflict between two help files.
	
	Note that once you reference a help topic or subtopic from your executable code,
	you must be careful to avoid changing the name of the topic or subtopic.
	
	The title parameter is used only for modal help and supplies the title for
	a modal dialog containing the help.
	
	The flags parameter is interpreted bitwise as follows:
		Bit 0			If cleared, Igor displays non-modal help.
						If set, Igor displays modal help.						
						
		Bit 1			If cleared, during modal help Igor displays the entire
						help file (if it is not too big) in the modal help dialog,
						with the specified topic initially in view. If set, during
						modal help Igor displays just the specified topic.
						
		Bit 2			If cleared, if the topic can not be found Igor displays an error
						dialog. If set, if the topic can not be found Igor does not display
						an error dialog.
						
		All other bits	Reserved - set to zero.
						
	You MUST set bit 0 if you call XOPDisplayHelpTopic from a modal dialog. This causes
	Igor do display a dialog containing help on top of your dialog. If you fail to set
	bit zero when calling XOPDisplayHelpTopic from a modal dialog, Igor may behave erratically.
	Unfortunately, links in help files don't work during modal help.
	
	If you are calling XOPDisplayHelpTopic in a non-modal situation, it is appropriate to
	clear bit zero, but not required. If you clear bit zero, Igor displays a normal
	Igor help file. If you set bit zero, Igor displays a modal help dialog.
	
	You must set all other bits to zero.

	Function result is 0 if OK or IGOR_OBSOLETE if you are running with an older
	version of IGOR or some other non-zero code if the topic can not be found.

	Added for Igor Pro 3.13B03.
*/
int
XOPDisplayHelpTopic(const char* title, const char* topicStr, long flags)
{
	return CallBack3(DISPLAY_HELP_TOPIC, (void*)title, (void*)topicStr, (void*)flags);
}

/*	XOPSetContextualHelpMessage(theWindow, message, r)

	Displays a message in the Igor Tips help window on Macintosh or in the
	status bar on Windows. Call this when your window is active and the user
	moves the cursor over an icon or other area of the window about which you
	have something to say.
	
	theWindow is your WindowRef on Macintosh or your HWND on Windows.
	This refers to the window containing the control or icon for which
	you are providing help.
	
	message is a C string containing the message to display.
	
	r is a pointer to a Macintosh rectangle, even on Windows, that indicates the
	area of the window that the icon occupies. When the user moves the cursor
	out of this rectangle, Igor will remove the message. On Macintosh, this
	rectangle is in the local coordinates of the window containing the control
	or icon. On Windows, it is in client coordinates of the window containing
	the control or icon. On Windows, use WinRectToMacRect to translate the Windows
	RECT into a Macintosh Rect.

	Function result is 0 if OK or IGOR_OBSOLETE.
	
	The WindowXOP1 sample XOP illustrates the use of this function.

	Added for Igor Pro Carbon. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
XOPSetContextualHelpMessage(XOP_WINDOW_REF theWindow, const char* message, const Rect* r)
{
	return CallBack3(SET_CONTEXTUAL_HELP_MESSAGE, (void*)theWindow, (void*)message, (void*)r);
}

/*	DoWindowRecreationDialog(procedureName)
	
	This routine is of use only to very advanced XOPs that support window recreation
	macros.

	On input, procedureName is the proposed name for the recreation macro.
	
	DoWindowRecreationDialog displays the standard Igor Close Window dialog that
	allows the user to enter a name for a recreation macro and choose to save the
	macro, replace an existing macro, skip saving the macro, or cancel.
	
	If the user clicks Save or Replace, on output, procedureName is the
	name chosen by the user for the recreation macro.
	
	DoWindowRecreationDialog returns the following codes:
		kCloseWinCancel			The user clicked the cancel button.
								The XOP must cancel the close of the window.
								
		kCloseWinSave			The user clicked the Save button.
								The XOP must save the macro in the Igor procedure
								window and close the window.
								
		kCloseWinReplace		The user clicked the Replace button.
								The XOP must save the macro in the Igor procedure
								window and close the window.
								
		kCloseWinNoSave			The user clicked the No Save button.
								The XOP must close the window without saving any macro.

	Added for Igor Pro 3.13B03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
enum CloseWinAction
DoWindowRecreationDialog(char* procedureName)
{
	return (enum CloseWinAction)CallBack1(WINDOW_RECREATION_DIALOG, procedureName);
}

/*	GetIgorProcedureList(hPtr, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later
	
	The main use for this routine is to check if a particular macro or function
	exists in the Igor procedure windows.
	
	GetIgorProcedureList returns via *hPtr a handle to a semicolon-separated
	list of procedure names. Depending on the flags parameter, the list may contain
	names of macros, names of user functions, or names of both macros and user functions.
	
	Note that the handle is not null terminated. Use GetHandleSize to determine how
	many characters are in the handle. This handle belongs to you, so call
	DisposeHandle to dispose it when you no longer need it.
	
	If Igor can not build the list, it returns a non-zero error code and sets
	*hPtr to NULL.
	
	The flags parameter is defined as follows:
		Bit 0:	If set, GetIgorProcedureList will list all macros.
				If cleared, it will ignore all macros.
		
		Bit 1:	If set, GetIgorProcedure will list all user functions.
				If cleared, it will ignore all user functions.
		
		All other bits are reserved for future use and must be set to 0.
	
	Igor will be unable to build the list if a syntactical error in the procedure
	files prevents Igor from successfully scanning them. In this case, GetIgorProcedureList
	will return NEED_COMPILE.
	
	GetIgorProcedureList can also return NOMEM if it runs out of memory. This is
	unlikely.
	
	Future versions of GetIgorProcedureList may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Added for Igor Pro 3.13B03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetIgorProcedureList(Handle* hPtr, long flags)
{
	return CallBack2(GET_IGOR_PROCEDURE_LIST, hPtr, (void*)flags);
}

/*	GetIgorProcedure(procedureName, hPtr, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later.
	
	The main use for this routine is to check if a particular macro or function
	exists in the Igor procedure windows.

	If Igor can find the procedure (macro or function) specified by procedureName,
	GetIgorProcedure returns via *hPtr a handle to the text for the procedure
	and returns a result of 0. The handle will contain the text for the specified
	procedure with a carriage return at the end of each line.
	
	Note that the handle is not null terminated. Use GetHandleSize to determine how
	many characters are in the handle. This handle belongs to you, so call
	DisposeHandle to dispose it when you no longer need it.
	
	If Igor can not find the procedure, it returns a non-zero error code and sets
	*hPtr to NULL.
	
	The flags parameter is defined as follows:
		Bit 0:	If set, GetIgorProcedure will look for macros with the
				specified name. If cleared, it will ignore all macros.
		
		Bit 1:	If set, GetIgorProcedure will look for user functions with the
				specified name. If cleared, it will ignore all user functions.
		
		All other bits are reserved for future use and must be set to 0.
	
	Igor will be unable to find the procedure if there is no such procedure. In
	this case, GetIgorProcedure will return NO_MACRO_OR_FUNCTION.
	
	Igor will be unable to find the procedure if a syntactical error in the procedure
	files prevents Igor from successfully scanning them. In this case, GetIgorProcedure
	will return NEED_COMPILE.
	
	GetIgorProcedure can also return NOMEM if it runs out of memory. This is
	unlikely.
	
	Future versions of GetIgorProcedure may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Added for Igor Pro 3.13B03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetIgorProcedure(const char* procedureName, Handle* hPtr, long flags)
{
	return CallBack3(GET_IGOR_PROCEDURE, (void*)procedureName, hPtr, (void*)flags);
}

/*	SetIgorProcedure(procedureName, h, flags)

	This routine is intended only for very advanced XOPs, such as XOPs that
	generate their own window recreation macros. An example is the WaveMetrics
	Surface Plotter XOP. It requires Igor Pro 3.13B03 or later
	
	SetIgorProcedure will fail if Igor procedures are running when it is called.
	It is intended to be called in response to a user action, such as when a user
	clicks the close box of an XOP window and you want to create a recreation macro
	for that window.
	
	SetIgorProcedure stores the procedure text, which is contained in the handle
	h, in a procedure window.
	
	You must make sure that procedureName agrees with the name of the procedure
	as found in the handle. Otherwise, Igor may be left in an unstable state.
	For example, if the text in h is this:
		Macro HelloWorld()
			Print "Hello World"
		End
	then procedureName must be "HelloWorld".
	
	The handle h must contain the lines of the procedure with a carriage return
	after each line. Do not include linefeeds. There must be no extraneous lines
	before or after the procedure text. The text in the handle should not be
	null terminated. There must be no garbage characters in the handle.
	
	The handle h belongs to Igor. Once you pass it to SetIgorProcedure, you
	must not modify it, access it, or dispose of it.

	The flags parameter is defined as follows:
		Bit 0:	Set if the text in the handle is a macro (defined by the keywords
				Macro, Proc or Window).
				
		Bit 1:	Set if the text in the handle is a function (defined by the keyword
				Function).
				
		Bit 2:	If cleared and if there is already a procedure with the same name,
				SetIgorProcedure will return an error. If set and if there is
				already a procedure with the same name, SetIgorProcedure will
				replace the existing procedure and return 0.
		
		All other bits are reserved for future use and must be set to 0.

	Thus, the flags parameter must be one of the following values:
		1		Text is a macro. Return error if conflict.
		2		Text is a function. Return error if conflict.
		5		Text is a macro. Overwrite existing procedure if conflict.
		6		Text is a function. Overwrite existing procedure if conflict.
		
	The flags parameter must be consistent with the text in the handle. For example,
	you must not pass 1 for flags if the text in the handle is a function.
	
	SetIgorProcedure normally stores the procedure in the built-in Procedure window.
	However, if there already exists a procedure with the specified name and if
	that procedure is in another procedure window, SetIgorProcedure stores the
	procedure in that procedure window.
	
	After storing the procedure in the procedure window, SetIgorProcedure causes
	the Igor procedure windows to be compiled if auto-compile is on. If auto-compile
	is off, then the user must cause procedures to be compiled by choosing Compile
	from the Macros menu. Because of this, an XOP can not depend on the procedure
	being available for execution as soon as SetIgorProcedure returns. There is
	currently no way for the XOP to know the state of the auto-compile setting.
	
	SetIgorProcedure returns 0 if it successfully stores the text in a procedure
	window. Otherwise, it returns a non-zero error code.
	
	SetIgorProcedure returns a non-zero error code in the following cases:

	1. SetIgorProcedure was called while Igor was compiling procedures.
	   This could happen only  under weird recursion circumstances.
	   SetIgorProcedure returns -1 in this case.
	   
	2. SetIgorProcedure was called while Igor was running a macro or function.
	   Igor can not change procedures while they are running because this could
	   cause a crash. SetIgorProcedure returns -1 in this case.
	   
	3. Igor can not store the text because a serious syntax error in the
	   procedure windows prevents it from successfully scanning the procedures.
	   In this case, SetIgorProcedure returns NEED_COMPILE.
	  
	4. procedureName is not a valid name (too long or contains bad characters).
	
	5. procedureName is in conflict with a built-in or external operation or function
	   or a wave or a variable. The error returned depends on the type of object
	   that is in conflict.
	
	6. procedureName is in conflict with a macro or user function and bit 2 of
	   flags is not set. SetIgorProcedure returns NAME_MAC_CONFLICT or
	   NAME_USER_FUNC_CONFLICT.
	
	7. SetIgorProcedure ran out of memory. In this case, SetIgorProcedure returns NOMEM.
	
	8. You are doing a replace and the procedure window containing the original
	   procedure is open for read-only. In this case, you will receive a non-zero
	   error code. If the procedure window is open for read/write but the write-protect
	   icon is on, SetIgorProcedure will replace the procedure and return 0. The write-
	   protect icon is intended only to prevent inadvertent manual changes.
	
	Future versions of SetIgorProcedure may return other error codes so your
	XOP should not crash or otherwise grossly misbehave if it receives some
	other error code.
	
	Added for Igor Pro 3.13B03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
SetIgorProcedure(const char* procedureName, Handle h, long flags)
{
	return CallBack3(SET_IGOR_PROCEDURE, (void*)procedureName, h, (void*)flags);
}

/*	GetFunctionInfo(name, fip)

	Returns information that you need in order to call an Igor user function
	or an external function from an XOP. This is an advanced feature that most
	XOP programmers will not need.
	
	name is the name of an existing user or external function. If there is no such
	function, GetFunctionInfo returns an error. If the function is a user function
	and procedures are not in a compiled state, GetFunctionInfo returns an
	error. If everything is OK, GetFunctionInfo returns zero.
	
	The information returned by GetFunctionInfo should be used and then discarded.
	If the user does anything to cause procedures to be compiled then the values
	returned by GetFunctionInfo are no longer valid.
	
	GetFunctionInfo returns via fip->compilationIndex a value that you will
	pass to CallFunction. This value is used to make sure that procedures
	have not been changed between the time you call GetFunctionInfo and
	the time you call CallFunction.
	
	GetFunctionInfo returns via fip->functionID a value which you will pass to
	CallFunction to specify which user function you want to execute.
	
	GetFunctionInfo returns via fip->subType a code that identifies certain
	special purpose functions. The value returned currently has no use but may
	be used in the future.
	
	GetFunctionInfo returns via fip->isExternalFunction a value that is non-zero
	for external functions and zero for user functions. This field is for your
	information only. Your code will be the same for external and user functions.
	
	GetFunctionInfo returns via fip->returnType one of the following codes:
		NT_FP64:			Return value is a double-precision number
		NT_FP64 | NT_CMPLX:	Return value is a complex double-precision number
		HSTRING_TYPE:		Return value is a string
	
	GetFunctionInfo returns via fip->numOptionalParameters, fip->numRequiredParameters
	and fip->totalNumParameters the number of optional, required and total parameters for
	the function. Currently, an XOP can call a user function that has optional parameters
	but the XOP can not pass optional parameters to the function. In other words, it
	must pass the required parameters only.
	
	GetFunctionInfo returns via fip->parameterTypes an array of parameter types.
	GetFunctionInfo stores a parameter type value in elements 0 through fip->totalNumParameters-1.
	Elements fip->totalNumParameters and higher are undefined so you must not use them.
	
	You must use the CheckFunctionForm XOPSupport routine to make sure that the
	function is of the form you want. You normally don't need to examine the parameter
	type values directly, but in case you are curious, here is how a parameter type value
	is defined:
		if (parameterType == NT_FP64)
			Parameter is a double-precision number
			
		if (parameterType == (NT_FP64 | NT_CMPLX))
			Parameter is a complex double-precision number
		
		if (parameterType == HSTRING_TYPE)
			Parameter is a string

		if (parameterType == FV_FUNC_TYPE)
			Parameter is a function reference
		
		if (WAVE_TYPE bit is set) {
			Parameter is a numeric or text wave
			if (WAVE_Z_TYPE bit is set)
				Wave declared with /Z flag (used only by Igor debugger)
			
			if (parameterType & 0x01) {
				Parameter is a complex numeric wave
			}
			else {
				if (isExternalFunction) {
					Parameter is a numeric or text wave (can't tell which till runtime)
				}
				else {
					if ((parameterType & 0xFE) == 0)
						Parameter is a text wave
					else
						Parameter is a numeric wave
				}
			
			}
		}
	
	In addition, the parameter type values for numeric, complex numeric and string
	parameters may be ORed with FV_REF_TYPE. This indicates that the corresponding
	parameter is "pass-by-reference", meaning that the function can change the
	value of that parameter.
	
	Added for Igor Pro 5.00. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetFunctionInfo(const char* name, FunctionInfoPtr fip)
{
	return CallBack2(GET_FUNCTION_INFO, (void*)name, fip);
}

/*	GetFunctionInfoFromFuncRef(fref, fip)

	Returns information that you need in order to call an Igor user function
	or an external function from an XOP. This is an advanced feature that most
	XOP programmers will not need.
	
	This function does the same thing as GetFunctionInfo except that the input parameter
	is the value of a FUNCREF field from an Igor structure instead of the name of the
	function.	
	
	Added for Igor Pro 5.03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetFunctionInfoFromFuncRef(FUNCREF fref, FunctionInfoPtr fip)
{
	return CallBack2(GET_FUNCTION_INFO_FROM_FUNCREF, (void*)fref, fip);
}

/*	CheckFunctionForm(fip, requiredNumParameters, requiredParameterTypes, badParameterNumberPtr, returnType)

	Checks the form (number of parameters, types of parameters and return type) of the function
	against the required form.  This is an advanced feature that most XOP programmers will not need.
	
	You must call CheckFunctionForm before calling CallFunction to make sure that the function
	you are calling has the form you expect. Otherwise you may cause a crash.
	
	fip is pointer to a structure set by calling GetFunctionInfo.
	
	requiredNumParameters is the number of parameters you expect the function to have.
	
	requiredParameterTypes is an array in which you have set each value to one of the following:
		NT_FP64					The parameter must be scalar numeric
		NT_FP64 | NT_CMPLX		The parameter must be complex numeric
		HSTRING_TYPE			The parameter must be a string
		WAVE_TYPE				The parameter must be a scalar numeric wave
		WAVE_TYPE | NT_CMPLX	The parameter must be a complex numeric wave
		TEXT_WAVE_TYPE			The parameter must be a text wave
		FV_FUNC_TYPE			The parameter must be a function reference
		
	The number of elements in requiredParameterTypes must be at least equal to requiredNumParameters. 

	If the parameter must be pass-by-reference, use FV_REF_TYPE in addition to the
	values above. For example:
		NT_FP64 | FV_REF_TYPE
		NT_FP64 | NT_CMPLX | FV_REF_TYPE
		HSTRING_TYPE | FV_REF_TYPE
	
	Pass-by-reference is applicable to numeric and string parameters only.
	
	If you do not want CheckFunctionForm to check a particular parameter, pass
	-1 in the corresponding element of the requiredParameterTypes array.
	
	If the function is an external function, if the function takes a wave parameter,
	there is no way to know if the external function expects a numeric wave or a text wave.
	Consequently, CheckFunctionForm does not distinguish between numeric and text
	waves for external functions. The external function itself is supposed to check
	the type of the wave passed in at runtime and return an error if it is the
	wrong type.
	
	returnType is the required return type of the function which must be one of
	the following:
		NT_FP64
		NT_FP64 | NT_CMPLX
		HSTRING_TYPE
		
	If you do not want CheckFunctionForm to check the return type, pass -1
	as the returnType parameter.
	
	Sets *badParameterNumberPtr to the zero-based index of the first parameter that does
	not match the required type or to -1 if all parameters match the required type.
	
	Returns 0 if the form matches the requirements or an error code if not.
	
	If a function parameter type does not match the required parameter type, the error
	code returned will indicate the type of parameter required but not which parameter
	type was bad. If you want to inform the user more specifically, use the value returned
	via badParameterNumberPtr to select your own more specific error code. If the
	error was a parameter type mismatch, *badParameterNumberPtr will contain the
	zero-based index of the bad parameter. Otherwise, *badParameterNumberPtr will
	contain -1.
	
	Added for Igor Pro 5.00. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
CheckFunctionForm(FunctionInfoPtr fip, int requiredNumParameters, int requiredParameterTypes[], int* badParameterNumberPtr, int returnType)
{
	return CallBack5(CHECK_FUNCTION_FORM, fip, (void*)requiredNumParameters, requiredParameterTypes, badParameterNumberPtr, (void*)returnType);
}

/*	CallFunction(fip, parameters, resultPtr)

	Calls the function identified by fip. fip is a pointer to a FunctionInfo structure
	whose values were set by calling GetFunctionInfo. This is an advanced feature that most
	XOP programmers will not need.
	
	fip->compilationIndex is used by CallFunction to make sure that procedures were
	not recompiled after you called GetFunctionInfo. If procedures were recompiled, then
	the information in the structure may be invalid so CallFunction returns BAD_COMPILATION_INDEX.
	
	parameters is a pointer to a structure containing the values that you want
	to pass to the function. These values must agree in type with the
	function's parameter list, as indicated by the parameter information that
	you obtain via GetFunctionInfo. To guarantee this, you must call CheckFunctionForm
	before calling CallFunction.
	
	NOTE: The parameters structure must use standard XOP structure packing, namely,
	two-byte packing. If you don't set the structure packing correctly, a crash
	is likely.
	
	Parameter type values are discussed in detail in the comments for the GetFunctionInfo
	function. Here is the basic correspondence between function parameter types and the
	respective structure field:
		if (parameterType == NT_FP64)
			structure field is double
			
		if (parameterType == (NT_FP64 | NT_CMPLX))
			structure field is double[2]
		
		if (parameterType == HSTRING_TYPE)
			structure field is Handle
		
		if (WAVE_TYPE bit is set)
			structure field is waveHndl
		
		if (FV_FUNC_TYPE bit is set)
			structure field is long
	
	NOTE: For pass-by-value strings parameters, ownership of a handle stored in
	the parameter structure is passed to the function when you call CallFunction.
	The function will dispose the handle and CallFunction will set the field
	to NULL. You must not dispose it or otherwise reference it after calling CallFunction.
	
	NOTE: For pass-by-reference string parameters, the handle stored in the parameter
	structure field may be reused or disposed by the called function. When
	CallFunction returns, you own the handle which may be the same handle
	you passed or a different handle. You own this handle and you must dispose
	it when you no longer need it. If the field is NULL, which could occur in the
	event of an error, you must not dispose of it or otherwise access it.
	
	CallFunction stores the function result at the location indicated by
	resultPtr. Here is the correspondence between the function result type and the
	variable pointed to by resultPtr:
		NT_FP64						double
		NT_FP64 | NT_CMPLX			double[2]
		HSTRING_TYPE				Handle
		
	NOTE: A function that returns a string can return NULL instead of a valid
	handle. You must test the returned value to make sure it is not NULL before using
	it. If the returned handle is not NULL, you own it and you must dispose it when
	you no longer need it.
	
	CallFunction returns 0 as the function result if the function executed.
	If the function could not be executed, it returns a non-zero error code.
	
	Added for Igor Pro 5.00. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
CallFunction(FunctionInfoPtr fip, void* parameters, void* resultPtr)
{
	return CallBack3(CALL_FUNCTION, fip, parameters, resultPtr);
}

/*	Example of Calling a User or External Function

	In this example, we have written our own curve fitting routine, analogous to Igor's
	FuncFit operation, as an external operation. We want to call a user or external function
	from our external operation. The function that we want to call has this form:
		Function FitFunc(w, x)
			Wave w
			Variable x
			
	To simplify the example, we assume that we known the name of the function that
	we want to execute.

	// Define the parameter structure. These are parameters we will pass to user or external function.
	#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.
	struct OurParams {							// Used to pass parameters to the function.
		waveHndl waveH;							// For the first function parameter.
		double x;								// For the second function parameter.
	};
	typedef struct OurParams OurParams;
	typedef struct OurParams* OurParamsPtr;
	#include "XOPStructureAlignmentReset.h"		// Reset structure alignment to default.

	int
	DoOurOperation(waveHndl coefsWaveH)
	{
		FunctionInfo fi;
		OurParams parameters;
		int badParameterNumber;
		int requiredParameterTypes[2];
		double result;
		int i;
		double values[5];
		int err;
		
		// Make sure the function exists and get information about it.
		if (err = GetFunctionInfo("TestFitFunc", &fi))
			return err;
		
		// Make sure the function has the right form.
		requiredParameterTypes[0] = NT_FP64;			// First parameter is numeric
		requiredParameterTypes[1] = WAVE_TYPE;			// Second parameter is a numeric wave.
		if (err = CheckFunctionForm(&fi, 2, requiredParameterTypes, &badParameterNumber, NT_FP64))
			return err;
		
		// We have a valid function. Let's call it.
		
		parameters.x = 0;
		parameters.waveH = coefsWaveH;
		for(i=0; i<5; i+=1) {
			parameters.x = i;
			if (err = CallFunction(&fi, (void*)&parameters, &result))
				return err;
			values[i] = result;
		}
			
		return 0;	
	}
*/

/*	RegisterOperation(cmdTemplate, runtimeNumVarList, runtimeStrVarList, runtimeParamStructSize, runtimeAddress, options)

	Registers an XOP operation with Igor's Operation Handler.
	
	cmdTemplate specifies the operation name and syntax.
	
	runtimeNumVarList is a semicolon-separated list of numeric variables that the operation sets
	at runtime or NULL if it sets no numeric variables.
	
	runtimeStrVarList is a semicolon-separated list of string variables that the operation sets
	at runtime or NULL if it sets no string variables.
	
	runtimeParamStructSize is the size of the runtime parameter structure for the operation.
	
	runtimeAddress is the address of the ExecuteOperation function for this operation.
	
	options is reserved for future use. Pass zero for this parameter.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
RegisterOperation(const char* cmdTemplate, const char* runtimeNumVarList, const char* runtimeStrVarList, int runtimeParamStructSize, void* runtimeAddress, int options)
{
	return CallBack6(REGISTER_OPERATION, (void*)cmdTemplate, (void*)runtimeNumVarList, (void*)runtimeStrVarList, (void*)runtimeParamStructSize, runtimeAddress, (void*)options);
}

/*	SetOperationNumVar(varName, dval)

	This is used only to implement an operation using Igor's Operation Handler.
	It is used to store a value in an external operation output numeric variable such as V_flag.
	
	varName must be the name of a numeric variable that you specified via the runtimeNumVarList
	parameter to RegisterOperation.
	
	dval is the value to be stored in the named variable.
	
	If your operation was invoked from the command line, SetOperationNumVar sets
	a global variable in the current data folder, creating it if necessary.
	
	If your operation was invoked from a macro, SetOperationNumVar sets
	a macro local variable.
	
	If your operation was invoked from a user function, SetOperationNumVar sets
	a function local variable.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
SetOperationNumVar(const char* varName, double dval)
{
	return CallBack2(SET_RUNTIME_NUMERIC_VARIABLE, (void*)varName, (void*)&dval);
}

/*	SetOperationStrVar(varName, str)

	This is used only to implement an operation using Igor's Operation Handler.
	It is used to store a value in an external operation output string variable such as S_fileName.
	
	varName must be the name of a string variable that you specified via the runtimeStrVarList
	parameter to RegisterOperation.
	
	str points to the value to be stored in the named variable.
	
	If your operation was invoked from the command line, SetOperationStrVar sets
	a global variable in the current data folder, creating it if necessary.
	
	If your operation was invoked from a macro, SetOperationStrVar sets
	a macro local variable.
	
	If your operation was invoked from a user function, SetOperationStrVar sets
	a function local variable.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
SetOperationStrVar(const char* varName, const char* str)
{
	return CallBack2(SET_RUNTIME_STRING_VARIABLE, (void*)varName, (void*)str);
}

/*	VarNameToDataType(varName, dataTypePtr)

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a variable
	as a parameter and then set the value of that variable. An example is Igor's
	Open operation which takes a "refNum" parameter.
	
	After calling VarNameToDataType you would then call StoreNumericDataUsingVarName
	or StoreStringDataUsingVarName.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.

	On output, *dataTypePtr will contain one of the following:
		0					Means varName refers to a string variable or SVAR.
		NT_FP64				Means varName refers to a scalar local variable or NVAR.
		NT_FP64 | NT_CMPLX	Means varName refers to a complex local variable or NVAR.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
VarNameToDataType(const char* varName, int* dataTypePtr)
{
	return CallBack2(VAR_NAME_TO_DATA_TYPE, (void*)varName, (void*)dataTypePtr);
}

/*	StoreNumericDataUsingVarName(varName, realPart, imagPart)

	Stores data in a numeric variable which may be local or global.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a numeric variable
	as a parameter and then set the value of that variable. An example is Igor's
	Open operation which takes a "refNum" parameter.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a numeric variable or NVAR. This will be the case if VarNameToDataType returns
	a non-zero number as the data type.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
StoreNumericDataUsingVarName(const char* varName, double realPart, double imagPart)
{
	return CallBack3(STORE_NUMERIC_DATA_USING_VAR_NAME, (void*)varName, (void*)&realPart, (void*)&imagPart);
}

/*	StoreStringDataUsingVarName(varName, buf, len)

	Stores data in a string variable which may be local or global.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which take the name of a string variable
	as a parameter and then set the value of that variable. An example is Igor's
	FReadLine operation which takes a "stringVarName" parameter.
	
	varName must be a pointer to a varName field in your operation's runtime parameter
	structure when a pointer to this structure has been passed by Igor to your ExecuteOperation
	function. Igor relies on the varName field being set correctly. If your operation was
	called from a user function, the varName field will actually contain binary data
	used by Igor to access function local variables, NVARs and SVARs.
	
	You should call this routine only after you have determined that varName refers to
	a string variable or SVAR. This will be the case if VarNameToDataType returns
	zero as the data type.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.0. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
StoreStringDataUsingVarName(const char* varName, const char* buf, long len)
{
	return CallBack3(STORE_STRING_DATA_USING_VAR_NAME, (void*)varName, (void*)buf, (void*)len);
}

/*	SetOperationWaveRef(waveH, waveRefIndentifier)

	Sets a wave reference in an Igor user-defined function to refer to a destination wave
	that your operation created.

	This is used only to implement an operation using Igor's Operation Handler.
	It must be called only from your ExecuteOperation function.
	
	This function is used only for operations which use a DataFolderAndName parameter
	and declare it as a destination wave parameter using this kind of syntax in the
	operation template:
	
		SampleOp DataFolderAndName:{dest,<real> or <complex> or <text>}
		
	When you use this syntax and your operation is compiled in a user-defined function,
	the Igor function compiler may automatically create a wave reference in the function
	for the destination wave. However the automatically created wave reference will be
	NULL until you set it by calling SetOperationWaveRef.
	
	The SetOperationWaveRef callback sets the automatically created wave reference to
	refer to a specific wave, namely the wave that you created in response to the
	DataFolderAndName parameter. You must call it after successfully creating the
	destination wave.
	
	Igor will create an automatic wave reference only if the operation is called
	from a user-defined function and only if the destination wave is specified using
	a simple name (e.g., wave0 but not root:wave0 or $destWave). You have no way to know
	whether an automatic wave reference was created so you must call SetOperationWaveRef
	in all cases. SetOperationWaveRef will do nothing if no automatic wave reference
	exists.	
	
	The waveRefIndentifier parameter allows Igor to determine where in memory the wave
	reference is stored. This information is passed to your XOP in the "ParamSet" field
	of your ExecuteOperation structure.
	
	Your code should look something like this:
	
		// In your RuntimeParams structure
		DataFolderAndName dest;
		int destParamsSet[1];
	
		// In your ExecuteOperation function
		destWaveH = NULL;
		err = MDMakeWave(&destWaveH, p->dest.name, p->dest.dfH, dimensionSizes, type, overwrite);
		if (destWaveH != NULL) {
			int waveRefIndentifier = p->destParamsSet[0];
			err = SetOperationWaveRef(destWaveH, waveRefIndentifier);
		}
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.04B05. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
SetOperationWaveRef(waveHndl waveH, int waveRefIndentifier)
{
	return CallBack2(SET_OPERATION_WAVE_REF, waveH, (void*)waveRefIndentifier);
}

/*	DateToIgorDateInSeconds(numValues, year, month, dayOfMonth, secs)
	
	Converts dates into Igor date format (seconds since 1/1/1904).
	
	numValues is the number of dates to convert.
	
	year, month and dayOfMonth and secs are arrays allocated by the calling routine.
	The size of each array is specified by numValues.
	
	On input, year, month and dayOfMonth hold the input values.
	On return, secs holds the output values.
	
	The function result is zero or an error code.
	
	This routine requires Igor Pro 5 or later. Earlier versions will return IGOR_OBSOLETE.
*/
int
DateToIgorDateInSeconds(int numValues, short* year, short* month, short* dayOfMonth, double* secs)
{
	return CallBack5(DATE_TO_IGOR_DATE, (void*)numValues, year, month, dayOfMonth, secs);
}

/*	IgorDateInSecondsToDate(numValues, secs, dates)
	
	Converts dates in Igor date format (seconds since 1/1/1904) into date records.
	
	numValues is the number of Igor dates to convert.
	
	secs is an array of dates in Igor date format. Its size is specified by numValues.

	dates is an array of shorts. It must hold 7*numValues shorts. For each input value,
	7 values are written to the dates array, in the following order:
		year, month, dayOfMonth, hour, minute, second, dayOfWeek		
	
	The function result is zero or an error code.
	
	Example:
		double secs[2];
		short dates[2*7];
		int err;
		
		secs[0] = 0;			// Represents January 1, 1904, midnight.
		secs[1] = 24*60*60;		// Represents January 2, 1904, midnight.
		err = IgorDateInSecondsToDate(2, secs, dates);
	
	This routine requires Igor Pro 5 or later. Earlier versions will return IGOR_OBSOLETE.
*/
int
IgorDateInSecondsToDate(int numValues, double* secs, short* dates)
{
	return CallBack3(IGOR_DATE_TO_DATE, (void*)numValues, secs, dates);
}

/*	GetNVAR(nvp, realPartPtr, imagPartPtr, numTypePtr)

	Retrieves the data and type of a global numeric variable referenced by an NVAR
	field in an Igor Pro structure.
	
	nvp is a pointer to an NVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	If GetNVAR returns 0 then *realPartPtr will be the contents of the real part of
	the global variable and *imagPartPtr will be the contents of the imaginary part of
	the global variable, if it is complex.
	
	realPartPtr and imagPartPtr must each point to storage for a double whether
	the global variable is complex or not.
	
	*numTypePtr is set to the numeric type of the global. This will be either NT_FP64
	or (NT_FP64 | NT_CMPLX).
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetNVAR(NVARRec* nvp, double* realPartPtr, double* imagPartPtr, int* numTypePtr)
{
	return CallBack4(GET_NVAR, nvp, realPartPtr, imagPartPtr, numTypePtr);
}

/*	SetNVAR(nvp, realPartPtr, imagPartPtr)

	Sets the value of a global numeric variable referenced by an NVAR field
	in an Igor Pro structure.

	nvp is a pointer to an NVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	*realPartPtr is the value to store in the real part of the global variable and
	*imagPartPtr is the value to store in the imaginary part of the global variable,
	if it is complex.
	
	realPartPtr and imagPartPtr must each point to storage for a double whether the global
	variable is complex or not.

	Returns 0 or an error code.

	Added for Igor Pro 5.03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int SetNVAR(NVARRec* nvp, double* realPartPtr, double* imagPartPtr)
{
	return CallBack3(SET_NVAR, nvp, realPartPtr, imagPartPtr);
}

/*	GetSVAR(svp, strHPtr)

	Retrieves the handle for a global string variable referenced by an SVAR field
	in an Igor Pro structure.

	svp is a pointer to an SVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	If GetSVAR returns 0 then *strHPtr will be the handle for an Igor global string variable.

	NOTE:	*strHPtr can be NULL if the global string variable contains no characters.

	NOTE:	*strHPtr belongs to Igor. Do not dispose it or alter it in any way.
			Use this function only to read the contents of the string.
		
	Returns 0 or an error code.
	
	Added for Igor Pro 5.03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int
GetSVAR(SVARRec* svp, Handle* strHPtr)
{
	return CallBack2(GET_SVAR, svp, strHPtr);
}

/*	SetSVAR(svp, strH)

	Sets the value of a global string variable referenced by an SVAR field
	in an Igor Pro structure.

	svp is a pointer to an SVAR field in an Igor structure.
	The Igor structure pointer would be passed into the XOP as a parameter.

	strH is a handle that you own. It can be NULL to set the global string variable
	to empty.

	NOTE:	Igor copies the contents of strH. You retain ownership of it and
			must dispose it.
	
	Returns 0 or an error code.
	
	Added for Igor Pro 5.03. If you call this with an earlier version of Igor,
	it will return IGOR_OBSOLETE and do nothing.
*/
int SetSVAR(SVARRec* svp, Handle strH)
{
	return CallBack2(SET_SVAR, svp, strH);
}

#include "XOPStructureAlignmentReset.h"	// Reset structure alignment to default.
