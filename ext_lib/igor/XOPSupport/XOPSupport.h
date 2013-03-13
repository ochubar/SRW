/*	XOPSupport.h

	Declares routines, global variables and other items needed to use XOPSupport routines.
*/

#ifdef __cplusplus
extern "C" {						// This allows C++ to call the XOPSupport routines.
#endif

#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

// Global variables used by XOPSupport.c.
extern IORecHandle XOPRecHandle;
extern int XOPModified;
extern int igorVersion;

// NaN represents a missing value.

#ifdef _MACINTOSH_
	// Visual C++ does not accept this syntax.
	#define SINGLE_NAN ((float)NAN)
	#define DOUBLE_NAN ((double)NAN)
#endif

#ifdef _WINDOWS_
	// Apple's ProjectBuilder (GNU C) does not accept this syntax.
	static unsigned char NAN_BYTES[] = {0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0x7F};
	#define SINGLE_NAN *(float *)NAN_BYTES
	#define DOUBLE_NAN *(double *)NAN_BYTES
#endif

/*	"Private" routines.
	The XOPSupport files use these to call Igor. You should not call them directly.
*/
long CallBack0(int message);
long CallBack1(int message, void* item0);
long CallBack2(int message, void* item0, void* item1);
long CallBack3(int message, void* item0, void* item1, void* item2);
long CallBack4(int message, void* item0, void* item1, void* item2, void* item3);
long CallBack5(int message, void* item0, void* item1, void* item2, void* item3, void* item4);
long CallBack6(int message, void* item0, void* item1, void* item2, void* item3, void* item4, void* item5);

// Notice Routines (in XOPSupport.c).
void XOPNotice(const char* noticePtr);
void XOPResNotice(int strListID, int index);

// Utility routines (in XOPSupport.c).
void Capitalize(char* name);
int CmpStr(const char* str1, const char* str2);
char* strchr2(const char* str, int ch);
char* strrchr2(const char* str, int ch);
void MemClear(void* p, long n);
int MoveLockHandle(void* h);
int GetCStringFromHandle(Handle h, char* str, int maxChars);
int PutCStringInHandle(const char* str, Handle h);
int CheckAbort(long timeOutTicks);
int IsNaN32(float* floatPtr);
int IsNaN64(double* doublePtr);
void SetNaN32(float* fPtr);
void SetNaN64(double* dPtr);
int IsINF32(float* floatPtr);
int IsINF64(double* doublePtr);
int IgorVersion(void);
void XOPInit(IORecHandle ioRecHandle);
void SetXOPType(long type);
void SetXOPEntry(void (*entryPoint)(void));
void SetXOPResult(long result);
long GetXOPResult(void);
void SetXOPMessage(int message);
long GetXOPMessage(void);
void SetXOPRefCon(long refCon);
long GetXOPRefCon(void);
long GetXOPStatus(void);
long GetXOPItem(int itemNumber);
void IgorError(const char* title, int errCode);
int GetIgorErrorMessage(long errCode, char errorMessage[256]);
int WinInfo(int index, int typeMask, char* name, XOP_WINDOW_REF* windowRefPtr);

// Numeric conversion utilities (in XOPNumericConversion.c).
#define SIGNED_INT 1
#define UNSIGNED_INT 2
#define IEEE_FLOAT 3

void DoubleToFloat(double* inPtr, float* outPtr, long numValues);
void DoubleToLong(double* inPtr, long* outPtr, long numValues);
void DoubleToShort(double* inPtr, short* outPtr, long numValues);
void DoubleToByte(double* inPtr, char* outPtr, long numValues);
void DoubleToUnsignedLong(double* inPtr, unsigned long* outPtr, long numValues);
void DoubleToUnsignedShort(double* inPtr, unsigned short* outPtr, long numValues);
void DoubleToUnsignedByte(double* inPtr, unsigned char* outPtr, long numValues);
int ConvertDouble(double* src, void* dest, long numValues, int destFormat, int destBytes);

void FloatToDouble(float* inPtr, double* outPtr, long numValues);
void FloatToLong(float* inPtr, long* outPtr, long numValues);
void FloatToShort(float* inPtr, short* outPtr, long numValues);
void FloatToByte(float* inPtr, char* outPtr, long numValues);
void FloatToUnsignedLong(float* inPtr, unsigned long* outPtr, long numValues);
void FloatToUnsignedShort(float* inPtr, unsigned short* outPtr, long numValues);
void FloatToUnsignedByte(float* inPtr, unsigned char* outPtr, long numValues);
int ConvertFloat(float* src, void* dest, long numValues, int destFormat, int destBytes);

void LongToDouble(long* inPtr, double* outPtr, long numValues);
void LongToFloat(long* inPtr, float* outPtr, long numValues);
void LongToShort(long* inPtr, short* outPtr, long numValues);
void LongToByte(long* inPtr, char* outPtr, long numValues);
void LongToUnsignedLong(long* inPtr, unsigned long* outPtr, long numValues);
void LongToUnsignedShort(long* inPtr, unsigned short* outPtr, long numValues);
void LongToUnsignedByte(long* inPtr, unsigned char* outPtr, long numValues);
int ConvertLong(long* src, void* dest, long numValues, int destFormat, int destBytes);

void ShortToDouble(short* inPtr, double* outPtr, long numValues);
void ShortToFloat(short* inPtr, float* outPtr, long numValues);
void ShortToLong(short* inPtr, long* outPtr, long numValues);
void ShortToByte(short* inPtr, char* outPtr, long numValues);
void ShortToUnsignedLong(short* inPtr, unsigned long* outPtr, long numValues);
void ShortToUnsignedShort(short* inPtr, unsigned short* outPtr, long numValues);
void ShortToUnsignedByte(short* inPtr, unsigned char* outPtr, long numValues);
int ConvertShort(short* src, void* dest, long numValues, int destFormat, int destBytes);

void ByteToDouble(char* inPtr, double* outPtr, long numValues);
void ByteToFloat(char* inPtr, float* outPtr, long numValues);
void ByteToLong(char* inPtr, long* outPtr, long numValues);
void ByteToShort(char* inPtr, short* outPtr, long numValues);
void ByteToUnsignedLong(char* inPtr, unsigned long* outPtr, long numValues);
void ByteToUnsignedShort(char* inPtr, unsigned short* outPtr, long numValues);
void ByteToUnsignedByte(char* inPtr, unsigned char* outPtr, long numValues);
int ConvertByte(char* src, void* dest, long numValues, int destFormat, int destBytes);

void UnsignedLongToDouble(unsigned long* inPtr, double* outPtr, long numValues);
void UnsignedLongToFloat(unsigned long* inPtr, float* outPtr, long numValues);
void UnsignedLongToLong(unsigned long* inPtr, long* outPtr, long numValues);
void UnsignedLongToShort(unsigned long* inPtr, short* outPtr, long numValues);
void UnsignedLongToByte(unsigned long* inPtr, char* outPtr, long numValues);
void UnsignedLongToUnsignedShort(unsigned long* inPtr, unsigned short* outPtr, long numValues);
void UnsignedLongToUnsignedByte(unsigned long* inPtr, unsigned char* outPtr, long numValues);
int ConvertUnsignedLong(unsigned long* src, void* dest, long numValues, int destFormat, int destBytes);

void UnsignedShortToDouble(unsigned short* inPtr, double* outPtr, long numValues);
void UnsignedShortToFloat(unsigned short* inPtr, float* outPtr, long numValues);
void UnsignedShortToLong(unsigned short* inPtr, long* outPtr, long numValues);
void UnsignedShortToShort(unsigned short* inPtr, short* outPtr, long numValues);
void UnsignedShortToByte(unsigned short* inPtr, char* outPtr, long numValues);
void UnsignedShortToUnsignedLong(unsigned short* inPtr, unsigned long* outPtr, long numValues);
void UnsignedShortToUnsignedByte(unsigned short* inPtr, unsigned char* outPtr, long numValues);
int ConvertUnsignedShort(unsigned short* src, void* dest, long numValues, int destFormat, int destBytes);

void UnsignedByteToDouble(unsigned char* inPtr, double* outPtr, long numValues);
void UnsignedByteToFloat(unsigned char* outPtr, float* fPtr, long numValues);
void UnsignedByteToLong(unsigned char* inPtr, long* outPtr, long numValues);
void UnsignedByteToShort(unsigned char* inPtr, short* outPtr, long numValues);
void UnsignedByteToByte(unsigned char* inPtr, char* outPtr, long numValues);
void UnsignedByteToUnsignedLong(unsigned char* inPtr, unsigned long* outPtr, long numValues);
void UnsignedByteToUnsignedShort(unsigned char* inPtr, unsigned short* outPtr, long numValues);
int ConvertUnsignedByte(unsigned char* src, void* dest, long numValues, int destFormat, int destBytes);

int NumTypeToNumBytesAndFormat(int numType, int* numBytesPerPointPtr, int* dataFormatPtr, int* isComplexPtr);
int NumBytesAndFormatToNumType(int numBytesPerValue, int dataFormat, int* numTypePtr);
void FixByteOrder(void* p, int bytesPerPoint, long count);
int ConvertData(void* src, void* dest, long numValues, int srcBytes, int srcFormat, int destBytes, int destFormat);
int ConvertData2(void* src, void* dest, long numValues, int srcDataType, int destDataType);

void ScaleDouble(double* dPtr, double* offset, double* multiplier, long numValues);
void ScaleFloat(float* fPtr, double* offset, double* multiplier, long numValues);
void ScaleLong(long* iPtr, double* offset, double* multiplier, long numValues);						// Added in release 3.0.
void ScaleShort(short* iPtr, double* offset, double* multiplier, long numValues);					// Added in release 3.0.
void ScaleByte(char* iPtr, double* offset, double* multiplier, long numValues);						// Added in release 3.0.
void ScaleUnsignedLong(unsigned long* iPtr, double* offset, double* multiplier, long numValues);	// Added in release 3.0.
void ScaleUnsignedShort(unsigned short* iPtr, double* offset, double* multiplier, long numValues);	// Added in release 3.0.
void ScaleUnsignedByte(unsigned char* iPtr, double* offset, double* multiplier, long numValues);	// Added in release 3.0.
void ScaleData(int dataType, void* dataPtr, double* offsetPtr, double* multiplierPtr, long numValues);
void ScaleClipAndRoundData(int dataType, void* dataPtr, unsigned long numValues, double offset, double multiplier, double dMin, double dMax, int doRound);

// Wave access routines (in XOPWaveAccess.c).
int FetchNumericValue(int type, const char* dataStartPtr, long index, double value[2]);
int StoreNumericValue(int type, char* dataStartPtr, long index, double value[2]);
void WaveHandleModified(waveHndl waveHandle);
void WaveHandlesModified(waveHndl waves[], int numWaves, long start[], long end[]);
void WaveModified(const char* waveName);
waveHndl FetchWave(const char* waveName);
waveHndl FetchWaveFromDataFolder(DataFolderHandle dataFolderH, const char* waveName);	// Added in Igor Pro 3.0.
int WaveType(waveHndl waveHandle);
long WavePoints(waveHndl waveHandle);
void WaveName(waveHndl waveHandle, char* namePtr);
void* WaveData(waveHndl waveHandle);
void WaveScaling(waveHndl waveHandle, double* hsAPtr, double* hsBPtr, double* topPtr, double* botPtr);
void SetWaveScaling(waveHndl waveHandle, double* hsAPtr, double* hsBPtr, double* topPtr, double* botPtr);
void WaveUnits(waveHndl waveHandle, char* xUnits, char* dataUnits);
void SetWaveUnits(waveHndl waveHandle, const char* xUnits, const char* dataUnits);
Handle WaveNote(waveHndl waveHandle);
void SetWaveNote(waveHndl waveHandle, Handle noteHandle);
unsigned long WaveModDate(waveHndl wavH);
int WaveLock(waveHndl wavH);										// Added in Igor Pro 5.0.
int SetWaveLock(waveHndl wavH, int lockState);						// Added in Igor Pro 5.0
int WaveModState(waveHndl wavH);									// Added in Igor Pro 4.0D08 but works with all versions of Igor.
int WaveModCount(waveHndl wavH);									// Added in Igor Pro 4.0D08.
long GetWavesInfo(waveHndl waves[], int numWaves, int waveTypes[], long wavePoints[], int waveStates[], void* wavePtrs[]);
void SetWavesStates(waveHndl waves[], int numWaves, int waveStates[]);
int MakeWave(waveHndl* waveHandlePtr, const char* waveName, long numPoints, int type, int overwrite);
int ChangeWave(waveHndl waveHandle, long numPoints, int type);
int KillWave(waveHndl waveHandle);

// Data folder access routines (in XOPDataFolderAccess.c).
int GetDataFolderNameOrPath(DataFolderHandle dataFolderH, int flags, char dataFolderPathOrName[MAXCMDLEN+1]);
int GetDataFolderIDNumber(DataFolderHandle dataFolderH, long* IDNumberPtr);
int GetDataFolderProperties(DataFolderHandle dataFolderH, long* propertiesPtr);
int SetDataFolderProperties(DataFolderHandle dataFolderH, long properties);
int GetDataFolderListing(DataFolderHandle dataFolderH, int optionsFlag, Handle h);
int GetRootDataFolder(long refNum, DataFolderHandle* rootFolderHPtr);
int GetCurrentDataFolder(DataFolderHandle* currentFolderHPtr);
int SetCurrentDataFolder(DataFolderHandle dataFolderH);
int GetNamedDataFolder(DataFolderHandle startingDataFolderH, char dataFolderPath[MAXCMDLEN+1], DataFolderHandle* dataFolderHPtr);
int GetDataFolderByIDNumber(long IDNumber, DataFolderHandle* dataFolderHPtr);
int GetParentDataFolder(DataFolderHandle dataFolderH, DataFolderHandle* parentFolderHPtr);
int GetNumChildDataFolders(DataFolderHandle parentDataFolderH, long* numChildDataFolderPtr);
int GetIndexedChildDataFolder(DataFolderHandle parentDataFolderH, long index, DataFolderHandle* childDataFolderHPtr);
int GetWavesDataFolder(waveHndl wavH, DataFolderHandle* dataFolderHPtr);
int NewDataFolder(DataFolderHandle parentFolderH, char newDataFolderName[MAX_OBJ_NAME+1], DataFolderHandle* newDataFolderHPtr);
int KillDataFolder(DataFolderHandle dataFolderH);
int DuplicateDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle parentDataFolderH, char newDataFolderName[MAX_OBJ_NAME+1]);
int MoveDataFolder(DataFolderHandle sourceDataFolderH, DataFolderHandle newParentDataFolderH);
int RenameDataFolder(DataFolderHandle dataFolderH, char newName[MAX_OBJ_NAME+1]);
int GetNumDataFolderObjects(DataFolderHandle dataFolderH, int objectType, long* numObjectsPtr);
int GetIndexedDataFolderObject(DataFolderHandle dataFolderH, int objectType, long index, char objectName[MAX_OBJ_NAME+1], DataObjectValuePtr objectValuePtr);
int KillDataFolderObject(DataFolderHandle dataFolderH, int objectType, char objectName[MAX_OBJ_NAME+1]);
int MoveDataFolderObject(DataFolderHandle sourceDataFolderH, int objectType, char objectName[MAX_OBJ_NAME+1], DataFolderHandle destDataFolderH);
int RenameDataFolderObject(DataFolderHandle dataFolderH, int objectType, char objectName[MAX_OBJ_NAME+1], char newObjectName[MAX_OBJ_NAME+1]);
int DuplicateDataFolderObject(DataFolderHandle dataFolderH, int objectType, char objectName[MAX_OBJ_NAME+1], DataFolderHandle destFolderH, char newObjectName[MAX_OBJ_NAME+1], int overwrite);
void ClearDataFolderFlags(void);
long GetDataFolderChangesCount(void);
int GetDataFolderChangeFlags(DataFolderHandle dataFolderH, long* flagsP);
int GetDataFolderObject(DataFolderHandle dataFolderH, char objectName[MAX_OBJ_NAME+1], int* objectTypePtr, DataObjectValuePtr objectValuePtr);
int SetDataFolderObject(DataFolderHandle dataFolderH, char objectName[MAX_OBJ_NAME+1], int objectType, DataObjectValuePtr objectValuePtr);

// Multi-dimension wave access routines (in XOPWaveAccess.c).
int MDMakeWave(waveHndl* wavHPtr, const char* waveName, DataFolderHandle dataFolderH, long dimensionSizes[MAX_DIMENSIONS+1], int type, int overwrite);
int MDGetWaveDimensions(waveHndl wavH, long* numDimensionsPtr, long dimensionSizes[MAX_DIMENSIONS+1]);
int MDChangeWave(waveHndl wavH, int dataType, long dimensionSizes[MAX_DIMENSIONS+1]);
int MDChangeWave2(waveHndl wavH, int dataType, long dimensionSizes[MAX_DIMENSIONS+1], int mode);	// Added in Igor Pro 5.04B06.
int MDGetWaveScaling(waveHndl wavH, int dimension, double* sfA, double* sfB);
int MDSetWaveScaling(waveHndl wavH, int dimension, double* sfA, double* sfB);
int MDGetWaveUnits(waveHndl wavH, int dimension, char units[MAX_UNIT_CHARS+1]);
int MDSetWaveUnits(waveHndl wavH, int dimension, char units[MAX_UNIT_CHARS+1]);
int MDGetDimensionLabel(waveHndl wavH, int dimension, long element, char label[MAX_DIM_LABEL_CHARS+1]);
int MDSetDimensionLabel(waveHndl wavH, int dimension, long element, char label[MAX_DIM_LABEL_CHARS+1]);
int MDAccessNumericWaveData(waveHndl wavH, int accessMode, long* dataOffsetPtr);
int MDGetNumericWavePointValue(waveHndl wavH, long indices[MAX_DIMENSIONS], double value[2]);
int MDSetNumericWavePointValue(waveHndl wavH, long indices[MAX_DIMENSIONS], double value[2]);
int MDGetDPDataFromNumericWave(waveHndl wavH, double* dPtr);
int MDStoreDPDataInNumericWave(waveHndl wavH, double* dPtr);
int MDGetTextWavePointValue(waveHndl wavH, long indices[MAX_DIMENSIONS], Handle textH);
int MDSetTextWavePointValue(waveHndl wavH, long indices[MAX_DIMENSIONS], Handle textH);

// Command line routines (in XOPSupport.c).
int XOPCommand(const char* cmdPtr);
int XOPSilentCommand(const char* cmdPtr);
int XOPCommand2(const char *cmdPtr, int silent, int sendToHistory);		// Added in Igor Pro 4.0D04.
void PutCmdLine(const char* cmd, int mode);
void FinishDialogCmd(const char* cmd, int mode);

// Variable access routines (in XOPSupport.c).
int FetchNumVar(const char* varName, double* doublePtr1, double* doublePtr2);
int StoreNumVar(const char* varName, double* doublePtr1, double* doublePtr2);
int FetchStrVar(const char* varName, char* stringPtr);
Handle FetchStrHandle(const char* varName);
int StoreStrVar(const char* varName, const char* stringPtr);
int Variable(const char* varName, int varType);
int VariableList(Handle listHandle, const char* match, const char* sep, long varTypeCode);
int StringList(Handle listHandle, const char* match, const char* sep);
int SetIgorIntVar(const char* numVarName, long value, int forceGlobal);
int SetIgorFloatingVar(const char* numVarName, double* valuePtr, int forceGlobal);
int SetIgorComplexVar(const char* numVarName, double* realValuePtr, double* imagValuePtr, int forceGlobal);
int SetIgorStringVar(const char* stringVarName, const char* stringVarValue, int forceGlobal);

// Command parsing routines (in XOPSupport.c).
int NextSymb(void);
int GetSymb(void);
int GetFlag(char* flagsPtr);
int GetName(char* namePtr);
int GetDataFolderAndName(DataFolderHandle* dataFolderHandlePtr, char* namePtr, int beLiberal);	// Added in Igor Pro 3.0.
int GetDataFolder(DataFolderHandle* dataFolderHandlePtr);										// Added in Igor Pro 3.0.
int Keyword(const char* keywords[], char* keyword);
int GetKeyword(const char* keywords[], int* keyPtr);
int GetWaveName(char* namePtr);
waveHndl GetWave(void);
int GetWaveList(waveHndl* waves, long* numWavesPtr);
int GetWaveRange(WaveRangeRecPtr wrp);
int CalcWaveRange(WaveRangeRecPtr wrp);
int GetNumVarName(char* namePtr);
int GetStrVarName(char* namePtr);
int GetNum(double* doublePtr);
int GetNum2(double* doublePtr, int terminator);
int GetFlagNum(double* doublePtr);
int GetTrueOrFalseFlag(long flagMask, long* flagsPtr);
int GetLong(long* longPtr);
int GetAString(char* stringPtr);
Handle GetAStringInHandle(void);
int CheckTerm(void);
int AtEndOfCommand(void);															// Added in Igor Pro 4.0 but works with earlier versions.
int UniqueName(const char* baseName, char* finalName);
int UniqueName2(int nameSpaceCode, const char* baseName, char* finalName, long* suffixNumPtr);
int SanitizeWaveName(char* waveName, long column);
int CheckName(DataFolderHandle dataFolderH, int objectType, const char* name);		// Added in Igor Pro 3.0.
int PossiblyQuoteName(char* name);													// Added in Igor Pro 3.0.
void CatPossiblyQuotedName(char* str, const char* name);							// Added in Igor Pro 3.0.
int CleanupName(int beLiberal, char* name, long maxNameChars);						// Added in Igor Pro 3.0.
int CreateValidDataObjectName(DataFolderHandle dataFolderH, const char* inName, char* outName, long* suffixNumPtr, int objectType, int beLiberal, int allowOverwrite, int inNameIsBaseName, int printMessage, int* nameChangedPtr, int* doOverwritePtr);

// As of Igor XOP Toolkit 5.0 (Carbon), this is no longer supported. Use GetName followed by GetPathInfo2 instead.
// int GetSymbolicPath(void);					// HR, 10/8/96: Changed name from GetPath.

int GetFormat(char* format);
int DoCHIO(CHIORecPtr CHIOPtr);
int IsStringExpression(int assumeString);		// Added in Igor Pro 3.0.

// Utilities for XOPs with menu items (in XOPMenus.c).
#ifdef _MACINTOSH_	// The Windows versions of these routines are implemented inside Igor and are declared in XOPWinMacSupport.h.
	void WMDrawMenuBar(void);
	void WMDeleteMenu(short menuID);
	void WMInsertMenu(MenuHandle menuH, short beforeID);
	// void ReleaseMenu(MenuHandle menuH);		// HR, 010427: This name ReleaseMenu was usurped by Apple. It is no longer supplied by the XOP Toolkit.
#endif
MenuHandle WMGetMenu(short resourceID);

// As of Igor XOP Toolkit 5.0 (Carbon), these are no longer supported. Use the XMI1 resource technique for adding a menu item instead.
// int GetXOPMenuID(void);
// int GetXOPItemID(void);
// int GetXOPSubMenuID(void);
// void GetXOPSubMenu(int resourceID);
// void EnableXOPMenuItem(int flag);				// HR, 10/8/96: Changed name from EnableMenuItem.
// void SetXOPItem(const char* itemString);
// void DisposeXOPSubMenu(void);

// As of Igor XOP Toolkit 5.0 (Carbon), this is no longer supported. See "Enabling and Disabling Menu Items" in the XOP Toolkit manual.
// void DefaultMenus(void);

int	ResourceToActualMenuID(int resourceMenuID);
MenuHandle ResourceMenuIDToMenuHandle(int resourceMenuID);
int	ActualToResourceMenuID(int menuID);
int ActualToResourceItem(int igorMenuID, int actualItemNumber);
int ResourceToActualItem(int igorMenuID, int resourceItemNumber);
int SetIgorMenuItem(int message, int enable, const char* text, long param);
void FillMenu(MenuHandle theMenu, const char* itemList, long itemListLen, int afterItem);
void FillMenuNoMeta(MenuHandle theMenu, const char* itemList, long itemListLen, int afterItem);
int FillWaveMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
int FillPathMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
int FillWinMenu(MenuHandle theMenu, const char* match, const char* options, int afterItem);
void WMDeleteMenuItems(MenuHandle theMenu, int afterItem);	// HR, 010427: This was previously called DeleteMenuItems but Apple usurped that name.

// Utilities for XOPs with windows (in XOPWindows.c).
XOP_WINDOW_REF GetActiveWindowRef(void);
int IsXOPWindowActive(XOP_WINDOW_REF w);
void ShowXOPWindow(XOP_WINDOW_REF w);
void HideXOPWindow(XOP_WINDOW_REF w);
void ShowAndActivateXOPWindow(XOP_WINDOW_REF w);
void HideAndDeactivateXOPWindow(XOP_WINDOW_REF w);
void SetXOPWindowTitle(XOP_WINDOW_REF windowRef, const char* title);
void GetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int* winStatePtr);
void SetXOPWindowPositionAndState(XOP_WINDOW_REF windowRef, Rect* r, int winState);
void TransformWindowCoordinates(int mode, double coords[4]);
void GetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int* winStatePtr);
void SetXOPWindowIgorPositionAndState(XOP_WINDOW_REF windowRef, double coords[4], int winState);
int TellIgorWindowStatus(XOP_WINDOW_REF windowRef, long status, long options);			// Added in Igor Pro 6.00D00.

// Utilities for XOPs with text windows (in XOPWindows.c).
#ifdef _MACINTOSH_
	// This routine is not supported on Windows.
	Handle TUNew(WindowPtr winPtr, Rect* borderRectPtr, int font, int size, int crOnly);
#endif
int TUNew2(const char* winTitle, const Rect* winRectPtr, Handle* TUPtr, XOP_WINDOW_REF* windowRefPtr);		// Added in Igor Pro 3.13.
void TUDispose(TUStuffHandle TU);
void TUDisplaySelection(TUStuffHandle TU);
void TUGrow(TUStuffHandle TU, long size);
void TUDrawWindow(TUStuffHandle TU);
void TUUpdate(TUStuffHandle TU);
void TUFind(TUStuffHandle TU, int code);
void TUReplace(TUStuffHandle TU);
void TUIndentLeft(TUStuffHandle TU);
void TUIndentRight(TUStuffHandle TU);
void TUClick(TUStuffHandle TU, EventRecord* eventPtr);
void TUActivate(TUStuffHandle TU, long flag);
void TUIdle(TUStuffHandle TU);
void TUNull(TUStuffHandle TU, EventRecord* eventPtr);
void TUCopy(TUStuffHandle TU);
void TUCut(TUStuffHandle TU);
void TUPaste(TUStuffHandle TU);
void TUClear(TUStuffHandle TU);
void TUKey(TUStuffHandle TU, EventRecord* eventPtr);
void TUInsert(TUStuffHandle TU, const char* dataPtr, long dataLen);
void TUDelete(TUStuffHandle TU);
void TUSetSelect(TUStuffHandle TU, long start, long end);	// TUSetSelect is obsolescent. Use TUSetSelLocs instead.
void TUSelectAll(TUStuffHandle TU);
void TUUndo(TUStuffHandle TU);
void TUPrint(TUStuffHandle TU);
void TUFixEditMenu(TUStuffHandle TU);
void TUFixFileMenu(TUStuffHandle TU);
Handle TUGetText(TUStuffHandle TU);							// TUGetText is obsolescent. Use TUFetchParagraphText instead.
void TUFetchText(TUStuffHandle TU, long offset, long numChars, char* buffer);	// TUFetchText is obsolescent. Use TUFetchParagraphText instead.
long TULength(TUStuffHandle TU);							// TULength is obsolescent. Use TUGetDocInfo instead.
long TULines(TUStuffHandle TU);
long TUSelStart(TUStuffHandle TU);							// TUSelStart is obsolescent. Use TUGetSelLocs instead.
long TUSelEnd(TUStuffHandle TU);							// TUSelEnd is obsolescent. Use TUGetSelLocs instead.
long TUSelectionLength(TUStuffHandle TU);					// TUSelectionLength is obsolescent. Use TUGetSelLocs instead.
int TUInsertFile(TUStuffHandle TU, const char* fileName, int wdRefNum);
int TUWriteFile(TUStuffHandle TU, const char* fileName, int wdRefNum, int allFlag);
int TUSFInsertFile(TUStuffHandle TU, const char* prompt, OSType fileTypes[], int numTypes);
int TUSFWriteFile(TUStuffHandle TU, const char* prompt, OSType fileType, int allFlag);
void TUPageSetupDialog(TUStuffHandle TU);
int TUGetDocInfo(TUStuffHandle TU, TUDocInfoPtr dip);
int TUGetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr);
int TUSetSelLocs(TUStuffHandle TU, TULocPtr startLocPtr, TULocPtr endLocPtr, long flags);
int TUFetchParagraphText(TUStuffHandle TU, long paragraph,  Ptr* textPtrPtr, long* lengthPtr);
int TUFetchSelectedText(TUStuffHandle TU, Handle* textHandlePtr, void* reservedForFuture, long flags);
int TUSetStatusArea(TUStuffHandle TU, const char* message, int eraseFlags, int statusAreaWidth);
void TUMoveToPreferredPosition(TUStuffHandle TU);
void TUMoveToFullSizePosition(TUStuffHandle TU);
void TURetrieveWindow(TUStuffHandle TU);

// Cross-platform dialog routines (in XOPDialogsMac.c and XOPDialogsWin.c)
void XOPOKAlert(const char* title, const char* message);
int XOPOKCancelAlert(const char* title, const char* message);
int XOPYesNoAlert(const char* title, const char* message);
int XOPYesNoCancelAlert(const char* title, const char* message);
void SetDialogBalloonHelpID(int balloonHelpID);
void GetDBox(XOP_DIALOG_REF theDialog, int itemID, Rect *box);
void HiliteDControl(XOP_DIALOG_REF theDialog, int itemID, int enable);
void DisableDControl(XOP_DIALOG_REF theDialog, int itemID);
void EnableDControl(XOP_DIALOG_REF theDialog, int itemID);
void SetRadBut(XOP_DIALOG_REF theDialog, int firstID, int lastID, int theID);
int ToggleCheckBox(XOP_DIALOG_REF theDialog, int itemID);
int GetCheckBox(XOP_DIALOG_REF theDialog, int itemID);
int SetCheckBox(XOP_DIALOG_REF theDialog, int itemID, int val);
void DisplayDialogCmd(XOP_DIALOG_REF theDialog, int itemID, const char* cmd);
#ifdef _MACINTOSH_
	CGrafPtr SetDialogPort(XOP_DIALOG_REF theDialog);
#endif
#ifdef _WINDOWS_
	XOP_DIALOG_REF SetDialogPort(XOP_DIALOG_REF theDialog);		// This is a NOP on Windows.
#endif
void ShowDialogWindow(XOP_DIALOG_REF theDialog);
int GetRadBut(XOP_DIALOG_REF theDialog, int itemID);
int GetDText(XOP_DIALOG_REF theDialog, int theItem, char *theText);
void SetDText(XOP_DIALOG_REF theDialog, int theItem, const char *theText);
int GetDInt(XOP_DIALOG_REF theDialog, int theItem, int *theInt);
void SetDInt(XOP_DIALOG_REF theDialog, int theItem, int theInt);
int GetDLong(XOP_DIALOG_REF theDialog, int theItem, long* theLong);
void SetDLong(XOP_DIALOG_REF theDialog, int theItem, long theLong);
int GetDDouble(XOP_DIALOG_REF theDialog, int theItem, double *theDouble);
void SetDDouble(XOP_DIALOG_REF theDialog, int theItem, double *theDouble);
void SelEditItem(XOP_DIALOG_REF theDialog, int itemID);
void SelMacEditItem(XOP_DIALOG_REF theDialog, int itemID);
int ItemIsPopMenu(XOP_DIALOG_REF theDialog, int itemID);
void InitPopMenus(XOP_DIALOG_REF theDialog);
int CreatePopMenu(XOP_DIALOG_REF theDialog, int popupItemNumber, int titleItemNumber, const char* itemList, int initialItem);	// Added in Igor Pro 3.1 but works with any version.
#ifdef _MACINTOSH_
	// This routine is not supported on Windows.
	MenuHandle GetPopMenuHandle(XOP_DIALOG_REF theDialog, int itemID);
#endif
void GetPopMenu(XOP_DIALOG_REF theDialog, int itemID, int *selItem, char *selStr);
int SetPopMatch(XOP_DIALOG_REF theDialog, int itemID, char *selStr);
void SetPopItem(XOP_DIALOG_REF theDialog, int itemID, int theItem);
void KillPopMenus(XOP_DIALOG_REF theDialog);
void AddPopMenuItems(XOP_DIALOG_REF theDialog, int itemID, const char *itemList);
void FillPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *itemList, long itemListLen, int afterItem);
int FillWavePopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
int FillPathPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
int FillWindowPopMenu(XOP_DIALOG_REF theDialog, int itemID, const char *match, const char *options, int afterItem);
void DeletePopMenuItems(XOP_DIALOG_REF theDialog, int itemID, int afterItem);

// Cross-platform file handling routines (in XOPFiles.c).
#ifdef _MACINTOSH_
	int HFSToPosixPath(const char* hfsPath, char posixPath[MAX_PATH_LEN+1], int isDirectory);
#endif
int XOPCreateFile(const char* fullFilePath, int overwrite, long macCreator, long macFileType);
int XOPDeleteFile(const char* fullFilePath);
int XOPOpenFile(const char* fullFilePath, int readOrWrite, XOP_FILE_REF* fileRefPtr);
int XOPCloseFile(XOP_FILE_REF fileRef);
int XOPReadFile(XOP_FILE_REF fileRef, unsigned long count, void* buffer, unsigned long* numBytesReadPtr);
int XOPReadFile2(XOP_FILE_REF fileRef, unsigned long count, void* buffer, unsigned long* numBytesReadPtr);
int XOPWriteFile(XOP_FILE_REF fileRef, unsigned long count, const void* buffer, unsigned long* numBytesWrittenPtr);
int XOPGetFilePosition(XOP_FILE_REF fileRef, unsigned long* filePosPtr);
int XOPSetFilePosition(XOP_FILE_REF fileRef, long filePos, int mode);
int XOPAtEndOfFile(XOP_FILE_REF fileRef);
int XOPNumberOfBytesInFile(XOP_FILE_REF fileRef, unsigned long* numBytesPtr);
int XOPReadLine(XOP_FILE_REF fileRef, char* buffer, unsigned long bufferLength, unsigned long* numBytesReadPtr);
int FullPathPointsToFile(const char* fullPath);
int FullPathPointsToFolder(const char* fullPath);
int WinToMacPath(char path[MAX_PATH_LEN+1]);
int MacToWinPath(char path[MAX_PATH_LEN+1]);
int GetNativePath(const char* filePathIn, char filePathOut[MAX_PATH_LEN+1]);
int EscapeBackslashesInUNCVolumeName(char macFilePath[MAX_PATH_LEN+1]);		// Added for XOP Toolkit 5.0, 991007.
int GetDirectoryAndFileNameFromFullPath(const char* fullFilePath, char dirPath[MAX_PATH_LEN+1], char fileName[MAX_FILENAME_LEN+1]);	// Added for Igor Pro 3.13.
int GetLeafName(const char* filePath, char name[MAX_FILENAME_LEN+1]);
int GetFullPathFromSymbolicPathAndFilePath(const char* symbolicPathName, const char filePath[MAX_PATH_LEN+1], char fullFilePath[MAX_PATH_LEN+1]);
int ConcatenatePaths(const char* pathIn1, const char* nameOrPathIn2, char pathOut[MAX_PATH_LEN+1]);
int XOPOpenFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, char filePath[MAX_PATH_LEN+1]);										// Added for IGOR Pro 3.13.
int XOPSaveFileDialog(const char* prompt, const char* fileFilterStr, int* fileIndexPtr, const char* initialDir, const char* defaultExtensionStr, char filePath[MAX_PATH_LEN+1]);	// Added for IGOR Pro 3.13.
// As of Igor XOP Toolkit 5.0 (Carbon), these are no longer supported.
// int GetStandardFilePath(char fullPathToDirectory[MAX_PATH_LEN+1]);
// int SetStandardFilePath(const char* fullPathToDirectory);
	
// File loader utilities (in XOPFiles.c).
// As of XOP Toolkit 5.0 (Carbon), FileLoaderGetOperationFlags is not supported. Use FileLoaderGetOperationFlags2 instead.
int FileLoaderGetOperationFlags2(const char* additionalFlagsChars, long* flagsPtr, char* baseName, char symbolicPathName[MAX_OBJ_NAME+1]);
int FileLoaderMakeWave(long column, char* waveName, long numPoints, int fileLoaderFlags, waveHndl* waveHandlePtr);
int SetFileLoaderOutputVariables(const char* fileName, int numWavesLoaded, const char* waveNames);
int SetFileLoaderOperationOutputVariables(int runningInUserFunction, const char* fileName, int numWavesLoaded, const char* waveNames);

// Data loading and saving utilities for internal WaveMetrics use only (used by WaveMetrics Browser). (In XOPFiles.c).
struct LoadDataInfo;		// CWPro 9 and GNU C requires this.
struct LoadFileInfo;		// CWPro 9 and GNU C requires this.
struct SaveDataInfo;		// CWPro 9 and GNU C requires this.
int PrepareLoadIgorData(struct LoadDataInfo* ldiPtr, long* refNumPtr, struct LoadFileInfo*** topFIHPtr);
int LoadIgorData(struct LoadDataInfo* ldiPtr, long refNum, struct LoadFileInfo** topFIH, DataFolderHandle destDataFolderH);
int EndLoadIgorData(struct LoadDataInfo* ldiPtr, long refNum, struct LoadFileInfo** topFIH);
int SaveIgorData(struct SaveDataInfo* sdiPtr, DataFolderHandle topDataFolderH);

// IGOR color table routines (in XOPSupport.c).
int GetIndexedIgorColorTableName(int index, char name[MAX_OBJ_NAME+1]);								// Added for IGOR Pro 3.1
int GetNamedIgorColorTableHandle(const char* name, IgorColorTableHandle* ictHPtr);					// Added for IGOR Pro 3.1
int GetIgorColorTableInfo(IgorColorTableHandle ictH, char name[MAX_OBJ_NAME+1], int* numColorsPtr);	// Added for IGOR Pro 3.1
int GetIgorColorTableValues(IgorColorTableHandle ictH, int startColorIndex, int endColorIndex, int updatePixelValues, IgorColorSpec* csPtr);	// Added for IGOR Pro 3.1

// Cross-Platform Utilities (in XOPSupport.c).
void WinRectToMacRect(const RECT* wr, Rect* mr);
void MacRectToWinRect(const Rect *mr, RECT *wr);

// Miscellaneous routines (in XOPSupport.c).
void XOPBeep(void);
void GetXOPIndString(char* text, int strID, int index);
void ArrowCursor(void);
void IBeamCursor(void);
void WatchCursor(void);
void HandCursor(void); 
void SpinCursor(void);
int SpinProcess(void);
int DoUpdate(void);
void PauseUpdate(long* savePtr);
void ResumeUpdate(long* savePtr);
int WaveList(Handle listHandle, const char* match, const char* sep, const char* options);
int WinList(Handle listHandle, const char* match, const char* sep, const char* options);
int PathList(Handle listHandle, const char* match, const char* sep, const char* options);
#ifdef _MACINTOSH_
	// This routine is not supported on Windows.
	int GetPathInfo(const char* pathName, long* vRefNumPtr, long*dirIDPtr, long*wdRefNumPtr);
#endif
int GetPathInfo2(const char* pathName, char fullDirPath[MAX_PATH_LEN+1]);	// Added for Igor Pro 3.13.
struct NamedFIFO** GetNamedFIFO(const char* name);
void MarkFIFOUpdated(struct NamedFIFO** fifo);
int SaveXOPPrefsHandle(Handle prefsHandle);				// Added for IGOR Pro 3.1 but can be used with any version.
int GetXOPPrefsHandle(Handle* prefsHandlePtr);			// Added for IGOR Pro 3.1 but can be used with any version.
int GetPrefsState(long* prefsStatePtr);					// Added for IGOR Pro 3.1.
int XOPDisplayHelpTopic(const char* title, const char* topicStr, long flags);		// Added for IGOR Pro 3.13B03.
enum CloseWinAction DoWindowRecreationDialog(char* procedureName);
int GetIgorProcedureList(Handle* hPtr, long flags);
int GetIgorProcedure(const char* procedureName, Handle* hPtr, long flags);
int SetIgorProcedure(const char* procedureName, Handle h, long flags);
int XOPSetContextualHelpMessage(XOP_WINDOW_REF theWindow, const char* message, const Rect* r);

int GetFunctionInfo(const char* name, FunctionInfoPtr fip);									// Added for Igor Pro 5.00.
int CheckFunctionForm(struct FunctionInfo* fip, int requiredNumParameters, int requiredParameterTypes[], int* badParameterNumberPtr, int returnType);	// Added for Igor Pro 5.00.
int CallFunction(struct FunctionInfo* fip, void* parameters, void* resultPtr);				// Added for Igor Pro 5.00.
int GetFunctionInfoFromFuncRef(FUNCREF fref, FunctionInfoPtr fip);							// Added for Igor Pro 5.03.

int RegisterOperation(const char* cmdTemplate, const char* runtimeNumVarList, const char* runtimeStrVarList, int runtimeParamStructSize, void* runtimeAddress, int options);	// Added for Igor Pro 5.00.
int SetOperationNumVar(const char* varName, double dval);									// Added for Igor Pro 5.00.
int SetOperationStrVar(const char* varName, const char* str);								// Added for Igor Pro 5.00.
int VarNameToDataType(const char* varName, int* dataTypePtr);								// Added for Igor Pro 5.00.
int StoreNumericDataUsingVarName(const char* varName, double realPart, double imagPart);	// Added for Igor Pro 5.00.
int StoreStringDataUsingVarName(const char* varName, const char* buf, long len);			// Added for Igor Pro 5.00.
int SetOperationWaveRef(waveHndl waveH, int waveRefIndentifier);							// Added for Igor Pro 5.04B05.

int DateToIgorDateInSeconds(int numValues, short* year, short* month, short* dayOfMonth, double* secs);		// Added for Igor Pro 5.00.
int IgorDateInSecondsToDate(int numValues, double* secs, short* dates);										// Added for Igor Pro 5.00.

int GetNVAR(NVARRec* nvp, double* realPartPtr, double* imagPartPtr, int* numTypePtr);		// Added for Igor Pro 5.03.
int SetNVAR(NVARRec* nvp, double* realPartPtr, double* imagPartPtr);						// Added for Igor Pro 5.03.
int GetSVAR(SVARRec* nvp, Handle* strHPtr);													// Added for Igor Pro 5.03.
int SetSVAR(SVARRec* nvp, Handle strH);														// Added for Igor Pro 5.03.

int GetTextWaveData(waveHndl waveH, int mode, Handle* textDataHPtr);						// Added for Igor Pro 5.04.
int SetTextWaveData(waveHndl waveH, int mode, Handle textDataH);							// Added for Igor Pro 5.04.
int GetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS]);			// Added for Igor Pro 5.04.
int SetWaveDimensionLabels(waveHndl waveH, Handle dimLabelsHArray[MAX_DIMENSIONS]);			// Added for Igor Pro 5.04.

#include "XOPStructureAlignmentReset.h"	// Reset structure alignment to default.

#ifdef __cplusplus
}
#endif
