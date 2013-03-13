#ifndef __GMINTERF_H
#define __GMINTERF_H

//-------------------------------------------------------------------------
// Platform- and compiler-dependent macro definitions
//-------------------------------------------------------------------------

#if !(defined(ALPHA_NONE) || defined(ALPHA_CLIENT) || defined(ALPHA__LIB__))
/*---------------For CodeWarrior PowerMac---------------*/
#if defined __POWERPC__
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(export)
#endif
/*---------------For CodeWarrior PC and Visual C++---------------*/
#elif defined __INTEL__ || defined WIN32
#if defined ALPHA__DLL__ || defined MATLAB_MEX_FILE
#define EXP __declspec(dllexport)
#else
#define EXP __declspec(dllimport)
#endif
#define CALL __stdcall
/*---------------For HP-UX, gcc---------------*/
#else
#endif
#endif /*ALPHA_NONE*/

#ifndef EXP
#define EXP
#endif

#ifndef CALL
#define CALL
#endif

//----------------------------------------------------------------------------
// Functions of the Library
//----------------------------------------------------------------------------

#ifdef __cplusplus  
extern "C" {
#endif

EXP int CALL gmMinInit1D(int MethNo, double xLower, double fLower, double xUpper, double fUpper, double xRelPrec, double* xNext, char* ErrWarnText);

EXP int CALL gmMinFindNextArg1D(int i, double fCur, double* xNext, double* xMin, char* ErrWarnText);

EXP int CALL gmMinInitMultiD(int MethNo, int NumDim, double* pArgsLower, double* pArgsUpper, double ArgRelPrec, double* pArgsNext, char* ErrWarnText);

EXP int CALL gmMinFindNextArgMultiD(int i, double fCur, double* pArgsNext, double* pArgsMin, double* fMin, char* ErrWarnText);

EXP int CALL gmMinFindMultiD(int MethNo, int NumDim, double (*pExtFunc)(double* pArgs, double* pArgsMin, double FMin), double* pArgsLower, double* pArgsUpper, double ArgRelPrec, int MaxIter, double* pArgsMin, double* pfMin, char* ErrWarnText);

EXP void CALL gmInterpCubSplinePrep(double *x, double *y, int n, double *y2, char* ErrWarnText);

EXP double CALL gmInterpCubSpline(double *xa, double *ya, double *y2a, int n, double x, char* ErrWarnText);

/** Creates an Interpolating Structure for a tabulated function.
@param meth [in] interpolation method number (1- cubic spline)
@param x [in] array of tabulated function arguments
@param y [in] array of tabulated function values
@param n [in] number of function values (length of x and y arrays)
@param ErrWarnText [out] error or warning text (to be allocated by calling application)
@return integer number referencying the Interpolating Structure
@author O.C.
*/
EXP int CALL gmInterpInit(int meth, double *x, double *y, int n, char* ErrWarnText);

/** Calculates function value using a given Interpolating Structure.
@param i [in] reference number of the Interpolating Structure (created by the function gmInterpInit)
@param x [in] argument for which the function value should be calculated
@param ErrWarnText [out] error or warning text (to be allocated by calling application)
@return the function value f(x)
@author O.C.
*/
EXP double CALL gmInterp(int i, double x, char* ErrWarnText);

/** Performs sorting of entities referenced by array of arrays of integer.
@param arrLen [in] array of lengths of arrays of references
@param numArr [in] number of arrays
@param pExtFitFunc
@param sizePopulation
@param numGenerations
@param typeMutation
@param typeCrossover
@param pAuxInData
@param arrSortResFlat sorted flat array of refernce numbers
@param ErrWarnText [out] error or warning text (to be allocated by calling application)
@return the minimal value of the fitness function found
@author O.C.
*/
EXP double CALL gmGeneticSort(int* arrLen, int numArr, double (*pExtFitFunc)(int* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int sizePopulation, int numGenerations, int typeMutation, double rateMutation, int typeXover, double rateXover, void* pAuxInData, int* arrSortResFlat, char* ErrWarnText);
//gmGeneticSort is obsolete, to be replaced by gmSortDiscr

/** Performs sorting of entities referenced by array of arrays of integer. Used for module sorting (planar or Apple-II undulators).
@param arLen [in] array of lengths of arrays of references
@param arNumCases [in] array of numbers of cases (or sub-references)
@param numAr [in] number of arrays
@param maxNumModif
@param pExtFitFunc
@param pClassFitFunc
@param methOpt
@param pOptParam
@param arSortResFlat sorted flat array of refernce numbers
@param ErrWarnText [out] error or warning text (to be allocated by calling application)
@return the minimal value of the fitness function found
@author O.C.
*/
EXP double CALL gmSortDiscr(int* arLen, int* arNumCases, int numAr, int maxNumModif, double (*pExtFitFunc)(int* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int methOpt, void* pOptParam, int* arSortResFlat, char* ErrWarnText);

EXP double CALL gmSortDiscrNonOrd(int* arFlatNumElemsInHeaps, int* arNumHeapsInTypes, int numTypes, int maxNumModif, double (*pExtFitFunc)(int* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int methOpt, void* pOptParam, int* arSortResFlat, char* ErrWarnText);

EXP double CALL gmGeneticOptimDiscrAndCont(int* arrLen, int numArr, double* wArr, int maxNumElemToModif, int numDim, double* dimRanges, double (*pExtFitFunc)(double* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int sizePopulation, int numGenerations, int typeMutation, double rateMutation, int typeXover, double rateXover, void* pAuxInData, double* arOptResFlat, char* ErrWarnText);
//gmGeneticOptimDiscrAndCont should be removed

EXP double CALL gmMinDiscrCont(int numElem, int numDim, double* dimRanges, double (*pExtFitFunc)(double* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int methOpt, void* pOptParam, double* arOptResFlat, char* ErrWarnText);

EXP double CALL gmMinDiscr(unsigned* arLen, int* arDiscrValMinMax, int* arDiscrValOrig, unsigned* arStructDescr, double* arWeights, unsigned geneDim, unsigned numAr, unsigned maxNumModif, double (*pExtFitFunc)(int* arInstFlat, void* pClassFitFunc), void* pClassFitFunc, int methOpt, void* pOptParam, int* arOptResFlat, char* ErrWarnText);

EXP void CALL gmLinAlgEqSolveSVD(double** a, long numRows, long numCols, double* b, double tolElim, double* x, char* ErrWarnText);

/** Makes 1D or 2D convolution of pSrc and pDst scaled arrays. The result is store in the pDst array.
@param dim [in] dimensions of the source and destination arrays (can be 1 or 2)
@param pSrc [in] 1D or 2D flat array of source data
@param pnSrc [in] length(s) of the source array
@param pxStartSrc [in] initial argument(s) of the source data 
@param pxStepSrc [in] step(s) of the source data 
@param pDst [in, out] 1D or 2D flat array of destination data (another array used in the convolution)
@param pnDst [in] length(s) of the destination array
@param pxStartDst [in] initial argument(s) of the destination data 
@param pxStepDst [in] step(s) of the destination data 
@param ErrWarnText [out] error or warning text
@author O.C.
*/
//EXP void CALL gmConvolve(int dim, double *pSrc, long *pnSrc, double *pxStartSrc, double *pxStepSrc, double *pDst, long *pnDst, double *pxStartDst, double *pxStepDst, char *ErrWarnText);

EXP int CALL gmObjDel(int i, char* ErrWarnText);

EXP int CALL gmObjDelAll(char* ErrWarnText);

#ifdef __cplusplus  
}
#endif

//----------------------------------------------------------------------------

#endif



