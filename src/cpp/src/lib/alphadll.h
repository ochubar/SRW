/* *********************************************************************
#
#
# Project:       Alpha
#
# Description:   defining a number of macro for several compilers 
#
# Author(s):     Oleg Chubar, Pascal Elleaume, Laurent Farvacque
#
# Original:	 May 2000
#
# 
# Copyright (c) 2000 by European Synchrotron Radiation Facility,
#                       Grenoble, France
#
#                       All Rights Reserved
#
#********************************************************************** */

#ifndef ALPHA_INC
#define ALPHA_INC

#ifndef ALPHA_NONE
/* _______________ For CodeWarrior PowerMac  ______________*/
#if defined __POWERPC__
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#ifndef EXP
#define EXP __declspec(export)
#endif
#endif

/* ________ For CodeWarrior PC and Visual C++_________*/
#elif defined __INTEL__ || defined WIN32
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#ifndef EXP
#define EXP __declspec(dllexport)
#endif
#else
#ifndef EXP
#define EXP __declspec(dllimport)
#endif
#endif
#ifndef CALL
#define CALL __stdcall
#endif

/* __________________ For HP-UX, gcc  ________________*/
#else
#endif

/* ___________________________________________________*/
#endif /*ALPHA_NONE*/

#ifndef EXP
#define EXP
#endif

#ifndef CALL
#define CALL
#endif

/* ___________________________________________________*/

#endif /*ALPHA_INC*/
