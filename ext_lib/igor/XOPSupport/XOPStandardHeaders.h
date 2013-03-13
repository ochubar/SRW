
#ifndef XOP_STANDARD_HEADERS						// Skip if XOP standard headers have already been included
	#define XOP_STANDARD_HEADERS 1
	#define _IGORXOP_ 1
	
	#include "XOPStructureAlignmentTwoByte.h"	// All structures passed between Igor and XOP are two-byte aligned.

	// The Windows headers create the _WINDOWS_ symbol if we are compiling for Windows.
	// Here, we create an analogous _MACINTOSH_ symbol if we are compiling for Macintosh.
	#if defined(TARGET_OS_MAC)
		#define _MACINTOSH_ 1
	#endif

	#ifdef _MACINTOSH_			// Compiling for Macintosh [
		#include <ctype.h>
		#include <string.h>
		#include <stdlib.h>
		#include <stdio.h>
		#include <stddef.h>
		// #include <math.h>		// This is included via fp.h in the precompiled headers.
		
		#ifdef __LITTLE_ENDIAN__
			#define XOP_LITTLE_ENDIAN	// We are compiling for Intel Macintosh.
		#endif

		typedef struct RECT {			// Windows RECT required for WinRectToMacRect and MacRectToWinRect in XOPSupport.c.
			long left;
			long top;
			long right;
			long bottom;
		} RECT;

		#define HOST_EXPORT				// Null definition for Macintosh. See below for details.
		#define HOST_IMPORT
	#endif						// Compiling for Macintosh ]

	// We use the WIN32 symbol to detect that we are compiling for Windows because _WINDOWS_ is not defined until we include Windows.h.
	#ifdef WIN32				// Compiling for Windows [
		#include <Windows.h>		// This creates the _WINDOWS_ symbol.
		
		#ifdef SetPort				// SetPort is defined in WinSpool.h
			#undef SetPort			// But we use SetPort in the Macintosh sense.
		#endif
		
		#ifdef _MSC_VER				// Microsoft Visual C++?
			#include "VCExtraIncludes.h"
		#endif
		
		#define XOP_LITTLE_ENDIAN	// HR, 030130: Was LITTLE_ENDIAN but Apple stole that from us (see endian.h).

		#undef DUPLICATE		// This is defined in XOP.h but also in NB30 for use with NetBios 3.0.

		/*	HOST_EXPORT declares a routine or variable exported from the host (IGOR) to one or more DLLs (XOPs)
			In IGOR, HOST_EXPORT is defined as __declspec(dllexport). However, in an XOP, HOST_EXPORT is
			defined as __declspec(dllimport). This allows IGOR and an XOP to share a prototype file
			where routines exported from IGOR are marked as HOST_EXPORT.
			
			HOST_IMPORT marks a routine imported by IGOR from a DLL. Currently, this is not
			actually used.
		*/
		#define HOST_EXPORT __declspec(dllimport)		// Declares function or variable as being exported from IGOR to XOPs.
		#define HOST_IMPORT __declspec(dllexport)		// Declares function or variable as being imported by IGOR from an XOP.

		#include <ctype.h>
		#include <string.h>
		#include <stdlib.h>
		#include <stdio.h>
		#include <stddef.h>
		#include <math.h>
		
		// This provides support for a small number of Macintosh routines that are needed by Windows XOPs.
		#include "XOPWinMacSupport.h"
	#endif						// Compiling for Windows ]
	
	#include "IgorXOP.h"
	#include "XOP.h"
	#include "XOPSupport.h"
	#ifdef _MACINTOSH_
		#include "XOPSupportMac.h"
	#endif
	#ifdef _WINDOWS_
		#include "XOPSupportWin.h"
	#endif

	#include "XOPStructureAlignmentReset.h"	// Reset structure alignment to default.

#endif				// XOP_STANDARD_HEADERS
