/************************************************************************//**
 * File: srwlpy.cpp
 * Description: Python binding
 * Project: Synchrotron Radiation Workshop Library (SRWLib)
 * First release: October 2010
 *
 * SRW is Copyright (c) European Synchrotron Radiation Facility, Grenoble, France
 * SRW C/C++ API (SRWLIB) is Copyright (c) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, G.Geloni, L.Samoylova
 * @version 0.066
 ***************************************************************************/

//#if defined(_DEBUG) 
//#undef _DEBUG
//#include "Python.h"
//#define _DEBUG
//#else
//#include "Python.h"
//#endif

#include "srwlib.h"
#include <vector>
#include <map>
#include <sstream> //OCTEST_161214

//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
//#include <time.h>

using namespace std;

//Without the following Python.h will enforce usage of python**_d.lib in dbug mode, which may not be always existing
//NOTE: to make it compilable with VC2013 (VC12), that blcock had to be moved down and placed after the previous includes
#if defined(_DEBUG) 
#undef _DEBUG
#include "Python.h"
#define _DEBUG
#else
#include "Python.h"
#endif

/************************************************************************//**
 * Error messages related to Python interface
 ***************************************************************************/
static const char strEr_NoObj[] = "No objects were submitted for parsing";
static const char strEr_BadList[] = "Incorrect or no Python List structure";
static const char strEr_BadArray[] = "Incorrect or no Python Array structure";
static const char strEr_BadListArray[] = "Incorrect or no Python List or Array structure";
static const char strEr_BadStr[] = "Error at parsing / converting Python string";
static const char strEr_BadClassName[] = "Error at retrieving Python class name";
static const char strEr_BadNum[] = "Incorrect or no Python number";
static const char strEr_BadTrj[] = "Incorrect Trajectory structure";
static const char strEr_BadPrt[] = "Incorrect Particle structure";
static const char strEr_BadPrtBm[] = "Incorrect Particle Beam structure";
static const char strEr_BadGsnBm[] = "Incorrect Gaussian Beam structure";
static const char strEr_BadPtSrc[] = "Incorrect Point Source structure";
static const char strEr_BadMagC[] = "Incorrect Magnetic Field Container structure";
static const char strEr_BadMag3D[] = "Incorrect 3D Magnetic Field structure";
static const char strEr_BadMagM[] = "Incorrect Multipole Magnet structure";
static const char strEr_BadMagS[] = "Incorrect Solenoid structure";
static const char strEr_BadMagU[] = "Incorrect Undulator structure";
static const char strEr_BadMagH[] = "Incorrect Undulator Harmonic structure";
static const char strEr_BadKickM[] = "Incorrect Kick Matrix structure";
static const char strEr_BadRadMesh[] = "Incorrect Radiation Mesh structure";
static const char strEr_BadWfr[] = "Incorrect Wavefront structure";
static const char strEr_BadStokes[] = "Incorrect Stokes parameters structure";
static const char strEr_BadOptC[] = "Incorrect Optical Element Container structure";
static const char strEr_BadOptD[] = "Incorrect Optical Drift structure";
static const char strEr_BadOptA[] = "Incorrect Optical Aperture / Obstacle structure";
static const char strEr_BadOptL[] = "Incorrect Optical Lens structure";
static const char strEr_BadOptAng[] = "Incorrect Optical Angle structure";
static const char strEr_BadOptShift[] = "Incorrect Optical Shift structure";
static const char strEr_BadOptZP[] = "Incorrect Optical Zone Plate structure";
static const char strEr_BadOptWG[] = "Incorrect Optical Waveguide structure";
static const char strEr_BadOptG[] = "Incorrect Optical Grating structure";
static const char strEr_BadOptT[] = "Incorrect Optical Generic Transmission structure";
static const char strEr_BadOptMir[] = "Incorrect Optical Mirror structure";
static const char strEr_BadOptCryst[] = "Incorrect Optical Crystal structure";
static const char strEr_BadListIntProp[] = "Incorrect list structure defining intensity distributions to be plotted after propagation";
static const char strEr_FloatArrayRequired[] = "This function can be executed for float array(s) only";
static const char strEr_FailedAllocPyArray[] = "Failed to allocate Python array from C";
static const char strEr_FailedUpdateInt[] = "Failed to update intensity data after propagation";
static const char strEr_FailedCreateList[] = "Failed to create resulting data list";

static const char strEr_BadArg_CalcMagnField[] = "Incorrect arguments for magnetic field calculation/tabulation function";
static const char strEr_BadArg_CalcPartTraj[] = "Incorrect arguments for trajectory calculation function";
static const char strEr_BadArg_CalcPartTrajFromKickMatr[] = "Incorrect arguments for trajectory calculation function from kick matrices";
static const char strEr_BadArg_CalcElecFieldSR[] = "Incorrect arguments for SR electric field calculation function";
static const char strEr_BadPrec_CalcElecFieldSR[] = "Incorrect precision parameters for SR electric field calculation";
static const char strEr_BadArg_CalcStokesUR[] = "Incorrect arguments for UR Stokes parameters calculation function";
static const char strEr_BadArg_CalcPowDenSR[] = "Incorrect arguments for SR power density calculation function";
static const char strEr_BadArg_CalcElecFieldGaussian[] = "Incorrect precision parameters for Gaussian beam electric field calculation";
static const char strEr_BadArg_CalcElecFieldSpherWave[] = "Incorrect precision parameters for spherical wave electric field calculation";
static const char strEr_BadArg_CalcIntFromElecField[] = "Incorrect arguments for intensity extraction function";
static const char strEr_BadArg_ResizeElecField[] = "Incorrect arguments for electric field resizing function";
static const char strEr_BadArg_SetRepresElecField[] = "Incorrect arguments for changing electric field representation function";
static const char strEr_BadArg_PropagElecField[] = "Incorrect arguments for electric field wavefront propagation function";
static const char strEr_BadArg_UtiFFT[] = "Incorrect arguments for FFT function";
static const char strEr_BadArg_UtiConvWithGaussian[] = "Incorrect arguments for convolution function";
static const char strEr_BadArg_UtiUndFromMagFldTab[] = "Incorrect arguments for magnetic field conversion to periodic function";
static const char strEr_BadArg_UtiUndFindMagFldInterpInds[] = "Incorrect arguments for magnetic field interpolaton index search function";
static const char strEr_BadArg_UtiIntInf[] = "Incorrect arguments for function analyzing intensity distributions";
static const char strEr_BadArg_UtiIntProc[] = "Incorrect arguments for function performing misc. operations on intensity distributions";

/************************************************************************//**
 * Global objects to be used across different function calls
 ***************************************************************************/
struct AuxStructPyObjectPtrs {
	PyObject *o_wfr;
	//Py_buffer pbEx, pbEy, pbMomX, pbMomY;
	Py_buffer pbEx, pbEy, pbExAux, pbEyAux, pbMomX, pbMomY; //OC151115
	vector<Py_buffer> *pv_buf;
};

static map<SRWLWfr*, AuxStructPyObjectPtrs> gmWfrPyPtr;
static map<char*, PyObject*> gmBufPyObjPtr; //OC16082018 (was added to enable allocation of intensity arrays in Py at propagation)

/************************************************************************//**
 * Auxiliary function dedicated to process errors reported by Library
 ***************************************************************************/
void ProcRes(int er) //throw(...) 
{
	char ErrorBuf[2048];

	if(er == 0) return;
	else
	{
		srwlUtiGetErrText(ErrorBuf, er);
		if(er < 0) 
		{//Print Warning:
			PyErr_SetString(PyExc_Warning, ErrorBuf);
			PyErr_PrintEx(1); //?
		}
		else throw ErrorBuf;
	}
}

/************************************************************************//**
 * Auxiliary function to erase element from map if key found
 ***************************************************************************/
//template<class Tkey, class T> void EraseElementFromMap(Tkey key, map<Tkey, T>& m) //GCC doesn't support this for some reason
void EraseElementFromMap(SRWLWfr* key, map<SRWLWfr*, AuxStructPyObjectPtrs>& m)
{
	map<SRWLWfr*, AuxStructPyObjectPtrs>::iterator iter = m.find(key);
	//map<SRWLWfr*, AuxStructPyObjectPtrs>::const_iterator iter = m.find(key);

	if(iter == m.end()) return;
	m.erase(iter);
}

/************************************************************************//**
 * Gets access to Py array buffer
 ***************************************************************************/
//char* GetPyArrayBuf(PyObject* obj, vector<Py_buffer>& vBuf, int bufFlag, Py_ssize_t* pSizeBuf) //sizeBuf is out
char* GetPyArrayBuf(PyObject* obj, vector<Py_buffer>* pvBuf, Py_ssize_t* pSizeBuf) //sizeBuf is out
{//for simplicity and uniformity of treatment in Py3 and Py2, only writable buffers are supported
	if(obj == 0) return 0;
	if(PyObject_CheckBuffer(obj))
	{
		Py_buffer pb_tmp;
		//if(PyObject_GetBuffer(obj, &pb_tmp, bufFlag)) return 0;
		if(PyObject_GetBuffer(obj, &pb_tmp, PyBUF_WRITABLE)) return 0;
		if(pSizeBuf != 0) *pSizeBuf = pb_tmp.len;
		if(pvBuf != 0) pvBuf->push_back(pb_tmp);
		return (char*)pb_tmp.buf;
	}
#if PY_MAJOR_VERSION < 3
	//else if(PyBuffer_Check(obj))
	else
	{
		//if(bufFlag != PyBUF_WRITABLE) return 0;
		PyObject *pOldBuf = PyBuffer_FromReadWriteObject(obj, 0, Py_END_OF_BUFFER);
		//if(pOldBuf == 0) return 0;
		if(pOldBuf == 0) //RN010814
		{
			PyErr_Clear(); return 0;
		}

		void *pVoidBuffer = 0;
		Py_ssize_t sizeBuf;
		char *pRes = 0;
		if(!PyObject_AsWriteBuffer(pOldBuf, &pVoidBuffer, &sizeBuf))
		{
			if(pVoidBuffer != 0) pRes = (char*)pVoidBuffer;
		}
		Py_DECREF(pOldBuf);
		if(pSizeBuf != 0) *pSizeBuf = sizeBuf;
		return pRes;
	}
#endif
	return 0;
}

/************************************************************************//**
 * Copies elements of Py list or array to a numerical array
 * ATTENTION: it can allocate T *ar !
 * arType can be 'i', 'f' or 'd'
 * Supports both Py lists and arrays
 * If obj is neither List nor Array - returns without thowing error
 ***************************************************************************/
//template<class T> void CopyPyListElemsToNumArray(PyObject* obj, char arType, T*& ar, int& nElem) //throw(...)
template<class T> char CopyPyListElemsToNumArray(PyObject* obj, char arType, T*& ar, int& nElem) //OC03092016
{
	//if(obj == 0) return;
	//if(!((arType == 'i') || (arType == 'l') || (arType == 'f') || (arType == 'd'))) return;
	if(obj == 0) return 0; //OC03092016
	if(!((arType == 'i') || (arType == 'l') || (arType == 'f') || (arType == 'd'))) return 0; //OC03092016
	//if(!PyList_Check(obj)) throw strEr_BadList;
	bool isList = PyList_Check(obj);
	bool isArray = false;
	if(!isList) isArray = PyObject_CheckBuffer(obj);

#if PY_MAJOR_VERSION >= 3
	//if(!(isList || isArray)) throw strEr_BadListArray;
	//if(!(isList || isArray)) return;
	if(!(isList || isArray)) return 0; //OC03092016
#endif

	Py_buffer pb;
	PyObject *pOldBuf=0;
	int *pIntAr=0;
	long *pLongAr=0;
	float *pFloatAr=0;
	double *pDoubleAr=0;

	int nElemInList = 0;
	if(isList) nElemInList = (int)PyList_Size(obj);
	else
	{
		void *pVoidBuffer = 0;
		Py_ssize_t sizeBuf = 0;

		if(isArray)
		{
			if(PyObject_GetBuffer(obj, &pb, PyBUF_SIMPLE)) throw strEr_BadArray;
			pVoidBuffer = pb.buf;
			sizeBuf = pb.len;
		}
#if PY_MAJOR_VERSION < 3
		else
		{
			//if(PyBuffer_Check(obj))
			//{
			pOldBuf = PyBuffer_FromReadWriteObject(obj, 0, Py_END_OF_BUFFER);
			if(pOldBuf != 0)
			{
				if(PyObject_AsWriteBuffer(pOldBuf, &pVoidBuffer, &sizeBuf)) throw strEr_BadArray;
				isArray = true;
				//if((pVoidBuffer == 0) || (sizeBuf <= 0)) throw strEr_BadArray;
				//Py_DECREF(pOldBuf);
			}
			else
			{
				PyErr_Clear();
				//return;
				return 0; //?
			}
			//}
			//else return; //obj is neither List nor Array
		}
#endif

		if(arType == 'i') 
		{
			//nElemInList = (int)pb.len/sizeof(int);
			//pIntAr = (int*)pb.buf;
			nElemInList = (int)(sizeBuf/sizeof(int));
			pIntAr = (int*)pVoidBuffer;
		}
		else if(arType == 'l') 
		{
			//nElemInList = (int)pb.len/sizeof(long);
			//pLongAr = (long*)pb.buf;
			nElemInList = (int)(sizeBuf/sizeof(long));
			pLongAr = (long*)pVoidBuffer;
		}
		else if(arType == 'f') 
		{
			//nElemInList = (int)pb.len/sizeof(float);
			//pFloatAr = (float*)pb.buf;
			nElemInList = (int)(sizeBuf/sizeof(float));
			pFloatAr = (float*)pVoidBuffer;
		}
		else if(arType == 'd') 
		{
			nElemInList = (int)(sizeBuf/sizeof(double));
			pDoubleAr = (double*)pVoidBuffer;
		}
	}
	if(nElemInList <=  0) throw strEr_BadListArray;

	if(ar == 0)
	{
		ar = new T[nElemInList];
		nElem = nElemInList;
	}
	else
	{
		if(nElem > nElemInList) nElem = nElemInList;
	}

	T *t_ar = ar;
	for(int i=0; i<nElem; i++)
	{
		if(isList)
		{
			PyObject *o = PyList_GetItem(obj, (Py_ssize_t)i);
			if(o == 0) throw strEr_BadNum;
			if(!PyNumber_Check(o)) throw strEr_BadNum;

			if((arType == 'i') || (arType == 'l')) *t_ar = (T)PyLong_AsLong(o);
			else if((arType == 'f') || (arType == 'd')) *t_ar = (T)PyFloat_AsDouble(o);
		}
		else if(isArray)
		{
			if(arType == 'i') *t_ar = (T)(*(pIntAr++));
			else if(arType == 'l') *t_ar = (T)(*(pLongAr++));
			else if(arType == 'f') *t_ar = (T)(*(pFloatAr++));
			else if(arType == 'd') *t_ar = (T)(*(pDoubleAr++));
		}
		t_ar++;
	}
	if(pOldBuf != 0) Py_DECREF(pOldBuf);

	return isList? 'l' : 'a'; //OC03092016
}

/************************************************************************//**
 * Copies elements of Py string to a C string (assumed to be allocated outside)
 ***************************************************************************/
void CopyPyStringToC(PyObject* pObj, char* c_str, int maxLenStr)
{
	if((pObj == 0) || (c_str == 0)) throw strEr_BadStr;

	int len = 0;
	char *pStr = 0;
	PyObject *pObjStr = 0;

	if(PyUnicode_Check(pObj))
	{
		pObjStr = PyUnicode_AsUTF8String(pObj);
		if(pObjStr != 0)
		{//usually works with Py3
			if(!PyBytes_Check(pObjStr)) throw strEr_BadStr;
			len = (int)PyBytes_Size(pObjStr);
			pStr = PyBytes_AsString(pObjStr); 
		}
	}
	else
	{//seems to work with Py2.7 and maybe earlier versions
		//pStr = PyString_AsString(pObj);
		Py_ssize_t lenLoc = 0;
#if PY_MAJOR_VERSION < 3
		int res = PyString_AsStringAndSize(pObj, &pStr, &lenLoc);
#else
		int res = PyBytes_AsStringAndSize(pObj, &pStr, &lenLoc);
#endif
		if(res < 0) throw strEr_BadStr;
		len = (int)lenLoc;
	}

	if((len > 0) && (pStr != 0))
	{
		if(len > maxLenStr) len = maxLenStr;
		strncpy(c_str, pStr, len);
		c_str[len] = '\0';
	}

	if(pObjStr != 0) Py_DECREF(pObjStr);
}

/************************************************************************//**
 * Copies char from C to Py
 ***************************************************************************/
PyObject* Py_BuildValueChar(char inC)
{
#if PY_MAJOR_VERSION >= 3
	return Py_BuildValue("C", inC); //doesn't work with Py2.7
#else
	return Py_BuildValue("c", inC);
#endif
}

/************************************************************************//**
 * Sets up output list (eventually of lists) data from an array
 ***************************************************************************/
template<class T> static PyObject* SetPyListOfLists(T* arB, int nB, int nP, char* cType="d") //OC13092018
{
	if((arB == 0) || (nB <= 0) || (nP <= 0)) return 0;
		
	int nElem = 0, nSubElem = 0;
	if(nP == 1)
	{
		nElem = nB;
	}
	else
	{
		nElem = nP;
		nSubElem = (int)round(nB/nP);
	}

	PyObject *oResB = PyList_New(nElem);
	T *t_arB = arB; //OC13092018
	//double *t_arB = arB;
	for(int i=0; i<nElem; i++)
	{
		PyObject *oElem = 0;
		if(nSubElem > 1)
		{
			oElem = PyList_New(nSubElem);
			for(int j=0; j<nSubElem; j++)
			{
				PyObject *oNum = Py_BuildValue(cType, *(t_arB++)); //OC13092018
				//PyObject *oNum = Py_BuildValue("d", *(t_arB++));
				if(PyList_SetItem(oElem, (Py_ssize_t)j, oNum)) throw strEr_FailedCreateList;
			}
		}
		else
		{
			oElem = Py_BuildValue(cType, *(t_arB++)); //OC13092018
			//oElem = Py_BuildValue("d", *(t_arB++));
		}
		if(PyList_SetItem(oResB, (Py_ssize_t)i, oElem)) throw strEr_FailedCreateList;
	}
	return oResB;
}

/************************************************************************//**
 * Copies elements of Py string to a C string (assumed to be allocated outside)
 ***************************************************************************/
void CopyPyClassNameToC(PyObject* pObj, char* c_str, int maxLenStr)
{
	if((pObj == 0) || (c_str == 0)) throw strEr_BadClassName;

	PyTypeObject *pTypeO = pObj->ob_type;
	if(pTypeO != 0)
	{
		const char *sTypeName = pTypeO->tp_name;
		if(sTypeName != 0)
		{
			if(strcmp(sTypeName, "instance") != 0)
			{
				int len = (int)strlen(sTypeName);
				if(len > maxLenStr) len = maxLenStr;
				strncpy(c_str, sTypeName, len);
				c_str[len] = '\0';
				return;
			}
		}
	}

#if PY_MAJOR_VERSION < 3
	PyObject *pStrReprO = PyObject_Repr(pObj);
	if(pStrReprO != 0)
	{
		char *sRepr = PyString_AsString(pStrReprO);
		if(sRepr != 0)
		{
			char *pFirstDot = strchr(sRepr, '.');
			if(pFirstDot != 0)
			{
				char *pEndName = strchr(pFirstDot, ' ');
				if(pEndName != 0) 
				{
					long len = (long)(pEndName - (pFirstDot + 1));
					if(len > 0)
					{
						if(len > maxLenStr) len = maxLenStr;
						strncpy(c_str, pFirstDot + 1, len);
						c_str[len] = '\0';
					}
				}
			}
		}
		Py_DECREF(pStrReprO);
	}
#endif
}

/************************************************************************//**
 * Parses PyObject* to SRWLParticle*
 ***************************************************************************/
void ParseSructSRWLParticle(SRWLParticle* pPrt, PyObject* oPrt) //throw(...) 
{
	if((pPrt == 0) || (oPrt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oPrt, "x");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "y");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "z");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->z = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "xp");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->xp = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "yp");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->yp = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "gamma");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->gamma = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "relE0");
	if(o_tmp == 0) throw strEr_BadPrt;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrt;
	pPrt->relE0 = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrt, "nq");
	if(o_tmp == 0) throw strEr_BadTrj;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadTrj;
	pPrt->nq = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLPartBeam*
 ***************************************************************************/
void ParseSructSRWLPartBeam(SRWLPartBeam* pPrtBm, PyObject* oPrtBm, vector<Py_buffer>& vBuf)
{
	if((pPrtBm == 0) || (oPrtBm == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	o_tmp = PyObject_GetAttrString(oPrtBm, "Iavg");
	if(o_tmp == 0) throw strEr_BadPrtBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrtBm;
	pPrtBm->Iavg = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrtBm, "nPart");
	if(o_tmp == 0) throw strEr_BadPrtBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPrtBm;
	pPrtBm->nPart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPrtBm, "partStatMom1");
	if(o_tmp == 0) throw strEr_BadPrtBm;
	ParseSructSRWLParticle(&(pPrtBm->partStatMom1), o_tmp);
	Py_DECREF(o_tmp);

	//pPrtBm->arStatMom2 = 0;
	o_tmp = PyObject_GetAttrString(oPrtBm, "arStatMom2");
	//if(o_tmp != 0)
	//{
	//	if(PyObject_CheckBuffer(o_tmp))
	//	{
	//		if(!PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE))
	//		{
	//			vBuf.push_back(pb_tmp);
	//			pPrtBm->arStatMom2 = (double*)pb_tmp.buf;
	//		}
	//	}
	//}
	double *pStatMom2 = pPrtBm->arStatMom2;
	int nMom2 = 21;
	CopyPyListElemsToNumArray(o_tmp, 'd', pStatMom2, nMom2);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLPrtTrj*
 * It fills-in pointers of SRWLPrtTrj without allocating memory and copying data
 * vector<Py_buffer>& vBuf is required to store and release all buffers after the execution
 ***************************************************************************/
void ParseSructSRWLPrtTrj(SRWLPrtTrj* pTrj, PyObject* oTrj, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pTrj == 0) || (oTrj == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	o_tmp = PyObject_GetAttrString(oTrj, "arX");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arX = (double*)pb_tmp.buf;
	//if(!(pTrj->arX = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arX = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "arXp");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arXp = (double*)pb_tmp.buf;
	//if(!(pTrj->arXp = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arXp = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "arY");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arY = (double*)pb_tmp.buf;
	//if(!(pTrj->arY = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arY = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "arYp");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arYp = (double*)pb_tmp.buf;
	//if(!(pTrj->arYp = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arYp = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "arZ");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arZ = (double*)pb_tmp.buf;
	//if(!(pTrj->arZ = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arZ = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "arZp");
	if(o_tmp == 0) throw strEr_BadTrj;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadTrj;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadTrj;
	//vBuf.push_back(pb_tmp);
	//pTrj->arZp = (double*)pb_tmp.buf;
	//if(!(pTrj->arZp = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0))) throw strEr_BadTrj;
	if(!(pTrj->arZp = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
	Py_DECREF(o_tmp);

	pTrj->arBx = 0;
	if(PyObject_HasAttrString(oTrj, "arBx"))
	{
		o_tmp = PyObject_GetAttrString(oTrj, "arBx");
		if(o_tmp != 0)
		{
			if(!(pTrj->arBx = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
			Py_DECREF(o_tmp);
		}
	}

	pTrj->arBy = 0;
	if(PyObject_HasAttrString(oTrj, "arBy"))
	{
		o_tmp = PyObject_GetAttrString(oTrj, "arBy");
		if(o_tmp != 0)
		{
			if(!(pTrj->arBy = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
			Py_DECREF(o_tmp);
		}
	}

	pTrj->arBz = 0;
	if(PyObject_HasAttrString(oTrj, "arBz"))
	{
		o_tmp = PyObject_GetAttrString(oTrj, "arBz");
		if(o_tmp != 0)
		{
			if(!(pTrj->arBz = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadTrj;
			Py_DECREF(o_tmp);
		}
	}

	o_tmp = PyObject_GetAttrString(oTrj, "np");
	if(o_tmp == 0) throw strEr_BadTrj;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadTrj;
	pTrj->np = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "ctStart");
	if(o_tmp == 0) throw strEr_BadTrj;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadTrj;
	pTrj->ctStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "ctEnd");
	if(o_tmp == 0) throw strEr_BadTrj;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadTrj;
	pTrj->ctEnd = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oTrj, "partInitCond");
	if(o_tmp == 0) throw strEr_BadTrj;
	ParseSructSRWLParticle(&(pTrj->partInitCond), o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLKickM*
 ***************************************************************************/
void ParseSructSRWLKickM(SRWLKickM* pKickM, PyObject* oKickM, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pKickM == 0) || (oKickM == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	o_tmp = PyObject_GetAttrString(oKickM, "nx");
	if(o_tmp == 0) throw strEr_BadKickM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->nx = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "ny");
	if(o_tmp == 0) throw strEr_BadKickM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->ny = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "nz");
	if(o_tmp == 0) throw strEr_BadKickM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->nz = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	//long npTot = long(pKickM->nx)*long(pKickM->ny);
	long long npTot = ((long long)(pKickM->nx))*((long long)(pKickM->ny));

	o_tmp = PyObject_GetAttrString(oKickM, "rx");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->rx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "ry");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->ry = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "rz");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->rz = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "x");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "y");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "z");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->z = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "order");
	if(o_tmp == 0) throw strEr_BadKickM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadKickM;
	pKickM->order = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	Py_ssize_t sizeBuf = 0;
	o_tmp = PyObject_GetAttrString(oKickM, "arKickMx");
	if(o_tmp == 0) throw strEr_BadKickM;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadKickM;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadKickM;
	//if(pb_tmp.len <= 0) pKickM->arKickMx = 0; 
	//else
	//{
	//	vBuf.push_back(pb_tmp);
	//	pKickM->arKickMx = (double*)pb_tmp.buf;
	//}
	//char *cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	char *cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf == 0) || (sizeBuf <= 0)) pKickM->arKickMx = 0;
	else
	{
		if((long long)sizeBuf != (long long)(npTot*sizeof(double))) throw strEr_BadKickM;
		pKickM->arKickMx = (double*)cpBuf;
	}
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oKickM, "arKickMy");
	if(o_tmp == 0) throw strEr_BadKickM;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadKickM;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadKickM;
	//if(pb_tmp.len <= 0) pKickM->arKickMy = 0; 
	//else
	//{
	//	vBuf.push_back(pb_tmp);
	//	pKickM->arKickMy = (double*)pb_tmp.buf;
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf == 0) || (sizeBuf <= 0)) pKickM->arKickMy = 0;
	else
	{
		if((long long)sizeBuf != (long long)(npTot*sizeof(double))) throw strEr_BadKickM;
		pKickM->arKickMy = (double*)cpBuf;
	}
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFld3D*
 * May allocate memory!
 ***************************************************************************/
void ParseSructSRWLMagFld3D(SRWLMagFld3D* pMag, PyObject* oMag, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	o_tmp = PyObject_GetAttrString(oMag, "nx");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->nx = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "ny");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->ny = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "nz");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->nz = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	//long npTot = long(pMag->nx)*long(pMag->ny)*long(pMag->nz);
	long long npTot = ((long long)(pMag->nx))*((long long)(pMag->ny))*((long long)(pMag->nz));

	o_tmp = PyObject_GetAttrString(oMag, "rx");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->rx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "ry");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->ry = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "rz");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->rz = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "nRep");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->nRep = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "interp");
	if(o_tmp == 0) throw strEr_BadMag3D;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMag3D;
	pMag->interp = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	Py_ssize_t sizeBuf = 0;
	o_tmp = PyObject_GetAttrString(oMag, "arBx");
	if(o_tmp == 0) throw strEr_BadMag3D;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMag3D;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//if(pb_tmp.len <= 0) pMag->arBx = 0; 
	//else
	//{
	//	vBuf.push_back(pb_tmp);
	//	pMag->arBx = (double*)pb_tmp.buf;
	//}
	//char *cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	char *cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf == 0) || (sizeBuf <= 0)) pMag->arBx = 0;
	else
	{
		if((long long)sizeBuf != (long long)(npTot*sizeof(double))) throw strEr_BadMag3D;
		pMag->arBx = (double*)cpBuf;
	}
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arBy");
	if(o_tmp == 0) throw strEr_BadMag3D;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMag3D;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//if(pb_tmp.len <= 0) pMag->arBy = 0; 
	//else
	//{
	//	vBuf.push_back(pb_tmp);
	//	pMag->arBy = (double*)pb_tmp.buf;
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf == 0) || (sizeBuf <= 0)) pMag->arBy = 0;
	else
	{
		if((long long)sizeBuf != (long long)(npTot*sizeof(double))) throw strEr_BadMag3D;
		pMag->arBy = (double*)cpBuf;
	}
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arBz");
	if(o_tmp == 0) throw strEr_BadMag3D;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMag3D;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//if(pb_tmp.len <= 0) pMag->arBz = 0; 
	//else
	//{
	//	vBuf.push_back(pb_tmp);
	//	pMag->arBz = (double*)pb_tmp.buf;
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf == 0) || (sizeBuf <= 0)) pMag->arBz = 0;
	else
	{
		if((long long)sizeBuf != (long long)(npTot*sizeof(double))) throw strEr_BadMag3D;
		pMag->arBz = (double*)cpBuf;
	}
	Py_DECREF(o_tmp);

	if((pMag->arBx == 0) && (pMag->arBy == 0) && (pMag->arBz == 0)) throw strEr_BadMag3D; //OC170515

	o_tmp = PyObject_GetAttrString(oMag, "arX");
	if(o_tmp == 0) throw strEr_BadMag3D;
	pMag->arX = 0;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//	if(pb_tmp.len > 0)
	//	{
	//		vBuf.push_back(pb_tmp);
	//		pMag->arX = (double*)pb_tmp.buf;
	//	}
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf != 0) && (sizeBuf > 0)) pMag->arX = (double*)cpBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arY");
	if(o_tmp == 0) throw strEr_BadMag3D;
	pMag->arY = 0;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//	if(pb_tmp.len > 0)
	//	{
	//		vBuf.push_back(pb_tmp);
	//		pMag->arY = (double*)pb_tmp.buf;
	//	}
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf != 0) && (sizeBuf > 0)) pMag->arY = (double*)cpBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arZ");
	if(o_tmp == 0) throw strEr_BadMag3D;
	pMag->arZ = 0;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMag3D;
	//	if(pb_tmp.len > 0)
	//	{
	//		vBuf.push_back(pb_tmp);
	//		pMag->arZ = (double*)pb_tmp.buf;
	//	}
	//}
	//cpBuf = GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf);
	cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
	if((cpBuf != 0) && (sizeBuf > 0)) pMag->arZ = (double*)cpBuf;
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFldQ*
 ***************************************************************************/
void ParseSructSRWLMagFldM(SRWLMagFldM* pMag, PyObject* oMag) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oMag, "G");
	if(o_tmp == 0) throw strEr_BadMagM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagM;
	pMag->G = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "m");
	if(o_tmp == 0) throw strEr_BadMagM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagM;
	pMag->m = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "n_or_s");
	if(o_tmp == 0) throw strEr_BadMagM;
	//PyObject *o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(!PyBytes_Check(o_str)) throw strEr_BadMagM;
	//pMag->n_or_s = *PyBytes_AsString(o_str);
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pMag->n_or_s = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "Leff");
	if(o_tmp == 0) throw strEr_BadMagM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagM;
	pMag->Leff = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "Ledge");
	if(o_tmp == 0) throw strEr_BadMagM;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagM;
	pMag->Ledge = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	pMag->R = 0;
	o_tmp = PyObject_GetAttrString(oMag, "R");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadMagM;
		pMag->R = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFldS*
 ***************************************************************************/
void ParseSructSRWLMagFldS(SRWLMagFldS* pMag, PyObject* oMag) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oMag, "B");
	if(o_tmp == 0) throw strEr_BadMagS;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagS;
	pMag->B = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "Leff");
	if(o_tmp == 0) throw strEr_BadMagS;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagS;
	pMag->Leff = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFldS*
 ***************************************************************************/
void ParseSructSRWLMagFldH(SRWLMagFldH* pMag, PyObject* oMag) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oMag, "n");
	if(o_tmp == 0) throw strEr_BadMagH;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagH;
	pMag->n = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "h_or_v");
	if(o_tmp == 0) throw strEr_BadMagH;
	//PyObject *o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(o_str != 0)
	//{//usually works with Py3
	//	if(!PyBytes_Check(o_str)) throw strEr_BadMagH;
	//	pMag->h_or_v = *PyBytes_AsString(o_str); 
	//}
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pMag->h_or_v = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "B");
	if(o_tmp == 0) throw strEr_BadMagH;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagH;
	pMag->B = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "ph");
	if(o_tmp == 0) throw strEr_BadMagH;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagH;
	pMag->ph = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "s");
	if(o_tmp == 0) throw strEr_BadMagH;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagH;
	pMag->s = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "a");
	if(o_tmp == 0) throw strEr_BadMagH;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagH;
	pMag->a = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFldU*
 * ATTENTION: allocates SRWLMagFldH *arHarm and its elements and char *arMagFldTypes
 ***************************************************************************/
void ParseSructSRWLMagFldU(SRWLMagFldU* pMag, PyObject* oMag) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	pMag->nHarm = 0;
	pMag->arHarm = 0;

	o_tmp = PyObject_GetAttrString(oMag, "per");
	if(o_tmp == 0) throw strEr_BadMagU;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagU;
	pMag->per = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "nPer");
	if(o_tmp == 0) throw strEr_BadMagU;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadMagU;
	pMag->nPer = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arHarm");
	if(o_tmp == 0) throw strEr_BadMagU;
	if(!PyList_Check(o_tmp)) throw strEr_BadMagU;
	//PyObject *o_List = o_tmp;
	int nHarm = (int)PyList_Size(o_tmp);
	if(nHarm <=  0) throw strEr_NoObj;
	pMag->nHarm = nHarm;
	pMag->arHarm = new SRWLMagFldH[nHarm];
	for(int i=0; i<nHarm; i++)
	{
		//PyObject *o = PyList_GetItem(o_List, (Py_ssize_t)i);
		PyObject *o = PyList_GetItem(o_tmp, (Py_ssize_t)i);
		ParseSructSRWLMagFldH(pMag->arHarm + i, o);
	}
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLWfr*
 * vector<Py_buffer>& vBuf is required to release all buffers after the end of execution
 ***************************************************************************/
//void ParseSructSRWLRadMesh(SRWLRadMesh* pRadMesh, PyObject* oRadMesh)
void ParseSructSRWLRadMesh(SRWLRadMesh* pRadMesh, PyObject* oRadMesh, vector<Py_buffer>* pvBuf =0)
{
	if((pRadMesh == 0) || (oRadMesh == 0)) throw strEr_NoObj;

	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oRadMesh, "eStart");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->eStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "eFin");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->eFin = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "xStart");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->xStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "xFin");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->xFin = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "yStart");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->yStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "yFin");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->yFin = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "zStart");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->zStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "ne");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->ne = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "nx");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->nx = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oRadMesh, "ny");
	if(o_tmp == 0) throw strEr_BadRadMesh;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
	pRadMesh->ny = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	//Optional parameters
	pRadMesh->nvx = 0.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "nvx");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->nvx = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
	pRadMesh->nvy = 0.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "nvy");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->nvy = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
	pRadMesh->nvz = 1.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "nvz");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->nvz = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
	pRadMesh->hvx = 1.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "hvx");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->hvx = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
	pRadMesh->hvy = 0.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "hvy");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->hvy = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}
	pRadMesh->hvz = 0.;
	o_tmp = PyObject_GetAttrString(oRadMesh, "hvz");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadRadMesh;
		pRadMesh->hvz = PyFloat_AsDouble(o_tmp);
		Py_DECREF(o_tmp);
	}

	pRadMesh->arSurf = 0;
	o_tmp = PyObject_GetAttrString(oRadMesh, "arSurf");
	if((o_tmp != 0) && (pvBuf != 0))
	{
		Py_ssize_t sizeBuf = 0;
		char *cpBuf = GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf);
		if((cpBuf != 0) && (sizeBuf > 0)) pRadMesh->arSurf = (double*)cpBuf;
		Py_DECREF(o_tmp);
	}
}

/************************************************************************//**
 * Parses PyObject* to SRWLMagFldCnt*
 * ATTENTION: allocates void **arMagFld and its elements and char *arMagFldTypes
 * vector<Py_buffer>& vBuf is required to store and release all buffers after the execution
 ***************************************************************************/
void ParseSructSRWLMagFldC(SRWLMagFldC* pMag, PyObject* oMag, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pMag == 0) || (oMag == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;
	Py_ssize_t sizeBuf = 0;

	o_tmp = PyObject_GetAttrString(oMag, "arMagFld");
	if(o_tmp == 0) throw strEr_BadMagC;
	if(!PyList_Check(o_tmp)) throw strEr_BadMagC;

	PyObject *o_List = o_tmp;
	//int nElem = (int)PyList_Size(o_tmp);
	int nElem = (int)PyList_Size(o_List);
	if(nElem <=  0) throw strEr_NoObj;
	//pMag->nElem = nElem;

	o_tmp = PyObject_GetAttrString(oMag, "arXc");
	if(o_tmp == 0) throw strEr_BadMagC;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMagC;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMagC;
	//if(pb_tmp.len != nElem*sizeof(double)) throw strEr_BadMagC;
	//vBuf.push_back(pb_tmp);
	//pMag->arXc = (double*)pb_tmp.buf;
	//if(!(pMag->arXc = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf))) throw strEr_BadMagC;
	if(!(pMag->arXc = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
	if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arYc");
	if(o_tmp == 0) throw strEr_BadMagC;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMagC;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMagC;
	//if(pb_tmp.len != nElem*sizeof(double)) throw strEr_BadMagC;
	//vBuf.push_back(pb_tmp);
	//pMag->arYc = (double*)pb_tmp.buf;
	//if(!(pMag->arYc = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf))) throw strEr_BadMagC;
	if(!(pMag->arYc = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
	if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oMag, "arZc");
	if(o_tmp == 0) throw strEr_BadMagC;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadMagC;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadMagC;
	//if(pb_tmp.len != nElem*sizeof(double)) throw strEr_BadMagC;
	//vBuf.push_back(pb_tmp);
	//pMag->arZc = (double*)pb_tmp.buf;
	//if(!(pMag->arZc = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, &sizeBuf))) throw strEr_BadMagC;
	if(!(pMag->arZc = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
	if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
	Py_DECREF(o_tmp);

	pMag->arVx = 0;
	if(PyObject_HasAttrString(oMag, "arVx"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arVx");
		if(o_tmp != 0) 
		{
			if(!(pMag->arVx = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if(sizeBuf == 0) pMag->arVx = 0;
			else if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
			Py_DECREF(o_tmp);
		}
	}
	pMag->arVy = 0;
	if(PyObject_HasAttrString(oMag, "arVy"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arVy");
		if(o_tmp != 0) 
		{
			if(!(pMag->arVy = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if(sizeBuf == 0) pMag->arVy = 0;
			else if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
			Py_DECREF(o_tmp);
		}
	}
	pMag->arVz = 0;
	if(PyObject_HasAttrString(oMag, "arVz"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arVz");
		if(o_tmp != 0) 
		{
			if(!(pMag->arVz = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if(sizeBuf == 0) pMag->arVz = 0;
			else if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
			Py_DECREF(o_tmp);
		}
	}
	pMag->arAng = 0;
	if(PyObject_HasAttrString(oMag, "arAng"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arAng");
		if(o_tmp != 0) 
		{
			if(!(pMag->arAng = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC;
			Py_DECREF(o_tmp);
		}
	}

	pMag->arPar1 = 0;
	if(PyObject_HasAttrString(oMag, "arPar1"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arPar1");
		if(o_tmp != 0)
		{
			if(!(pMag->arPar1 = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC; //?
			Py_DECREF(o_tmp);
		}
	}
	pMag->arPar2 = 0;
	if(PyObject_HasAttrString(oMag, "arPar2"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arPar2");
		if(o_tmp != 0)
		{
			if(!(pMag->arPar2 = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC; //?
			Py_DECREF(o_tmp);
		}
	}
	pMag->arPar3 = 0;
	if(PyObject_HasAttrString(oMag, "arPar3"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arPar3");
		if(o_tmp != 0)
		{
			if(!(pMag->arPar3 = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC; //?
			Py_DECREF(o_tmp);
		}
	}
	pMag->arPar4 = 0;
	if(PyObject_HasAttrString(oMag, "arPar4"))
	{
		o_tmp = PyObject_GetAttrString(oMag, "arPar4");
		if(o_tmp != 0)
		{
			if(!(pMag->arPar4 = (double*)GetPyArrayBuf(o_tmp, pvBuf, &sizeBuf))) throw strEr_BadMagC;
			if((long long)sizeBuf != (long long)(nElem*sizeof(double))) throw strEr_BadMagC; //?
			Py_DECREF(o_tmp);
		}
	}

	pMag->arMagFld = new void*[nElem];
	pMag->arMagFldTypes = new char[nElem + 1];
	pMag->arMagFldTypes[nElem] = '\0';

	pMag->nElem = 0; //in case if there will be reading error
	for(int i=0; i<nElem; i++)
	{
		PyObject *o = PyList_GetItem(o_List, (Py_ssize_t)i);
		
		//PyTypeObject *pTypeO = o->ob_type;
		//const char* sTypeName = pTypeO->tp_name;
		char sTypeName[1025];
		CopyPyClassNameToC(o, sTypeName, 1024);

		//if((strcmp(sTypeName, "SRWLMagFldC") == 0) || (cFldTypeID == 'c'))
		if(strcmp(sTypeName, "SRWLMagFldC") == 0)
		{
			SRWLMagFldC *pMagElem = new SRWLMagFldC();
			pMag->arMagFldTypes[i] = 'c';
			pMag->arMagFld[i] = (void*)pMagElem;
			ParseSructSRWLMagFldC(pMagElem, o, pvBuf);
		}
		//else if((strcmp(sTypeName, "SRWLMagFld3D") == 0) || (cFldTypeID == 'a'))
		else if(strcmp(sTypeName, "SRWLMagFld3D") == 0)
		{
			SRWLMagFld3D *pMagElem = new SRWLMagFld3D();
			pMag->arMagFldTypes[i] = 'a';
			pMag->arMagFld[i] = (void*)pMagElem;
			ParseSructSRWLMagFld3D(pMagElem, o, pvBuf);
		}
		//else if((strcmp(sTypeName, "SRWLMagFldM") == 0) || (cFldTypeID == 'm'))
		else if(strcmp(sTypeName, "SRWLMagFldM") == 0)
		{
			SRWLMagFldM *pMagElem = new SRWLMagFldM();
			pMag->arMagFldTypes[i] = 'm';
			pMag->arMagFld[i] = (void*)pMagElem;
			ParseSructSRWLMagFldM(pMagElem, o);
		}
		//else if((strcmp(sTypeName, "SRWLMagFldS") == 0) || (cFldTypeID == 's'))
		else if(strcmp(sTypeName, "SRWLMagFldS") == 0)
		{
			SRWLMagFldS *pMagElem = new SRWLMagFldS();
			pMag->arMagFldTypes[i] = 's';
			pMag->arMagFld[i] = (void*)pMagElem;
			ParseSructSRWLMagFldS(pMagElem, o);
		}
		//else if((strcmp(sTypeName, "SRWLMagFldU") == 0) || (cFldTypeID == 'u'))
		else if(strcmp(sTypeName, "SRWLMagFldU") == 0)
		{
			SRWLMagFldU *pMagElem = new SRWLMagFldU();
			pMag->arMagFldTypes[i] = 'u';
			pMag->arMagFld[i] = (void*)pMagElem;
			ParseSructSRWLMagFldU(pMagElem, o);
		}
		//to add more magnetic elements

		(pMag->nElem)++; //in case if there will be reading error
	}
	Py_DECREF(o_List);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptD*
 ***************************************************************************/
void ParseSructSRWLOptD(SRWLOptD* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "L");
	if(o_tmp == 0) throw strEr_BadOptD;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptD;
	pOpt->L = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "treat");
	if(o_tmp != 0)
	{
		if(!PyNumber_Check(o_tmp)) throw strEr_BadOptD;
		pOpt->treat = (char)PyLong_AsLong(o_tmp);
		Py_DECREF(o_tmp);
	}

}

/************************************************************************//**
 * Parses PyObject* to SRWLOptA*
 ***************************************************************************/
void ParseSructSRWLOptA(SRWLOptA* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "shape");
	if(o_tmp == 0) throw strEr_BadOptA;
	//PyObject *o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(!PyBytes_Check(o_str)) throw strEr_BadOptA;
	//pOpt->shape = *PyBytes_AsString(o_str); 
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pOpt->shape = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "ap_or_ob");
	if(o_tmp == 0) throw strEr_BadOptA;
	//o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(!PyBytes_Check(o_str)) throw strEr_BadOptA;
	//pOpt->ap_or_ob = *PyBytes_AsString(o_str); 
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pOpt->ap_or_ob = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Dx");
	if(o_tmp == 0) throw strEr_BadOptA;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptA;
	pOpt->Dx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Dy");
	if(o_tmp == 0) throw strEr_BadOptA;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptA;
	pOpt->Dy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "x");
	if(o_tmp == 0) throw strEr_BadOptA;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptA;
	pOpt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "y");
	if(o_tmp == 0) throw strEr_BadOptA;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptA;
	pOpt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptL*
 ***************************************************************************/
void ParseSructSRWLOptL(SRWLOptL* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "Fx");
	if(o_tmp == 0) throw strEr_BadOptL;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptL;
	pOpt->Fx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Fy");
	if(o_tmp == 0) throw strEr_BadOptL;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptL;
	pOpt->Fy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "x");
	if(o_tmp == 0) throw strEr_BadOptL;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptL;
	pOpt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "y");
	if(o_tmp == 0) throw strEr_BadOptL;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptL;
	pOpt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptAng*
 ***************************************************************************/
void ParseSructSRWLOptAng(SRWLOptAng* pOpt, PyObject* oOpt)
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "AngX");
	if(o_tmp == 0) throw strEr_BadOptAng;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptAng;
	pOpt->AngX = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "AngY");
	if(o_tmp == 0) throw strEr_BadOptAng;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptAng;
	pOpt->AngY = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptAng*
 ***************************************************************************/
void ParseSructSRWLOptShift(SRWLOptShift* pOpt, PyObject* oOpt)
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "ShiftX");
	if(o_tmp == 0) throw strEr_BadOptAng;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptShift;
	pOpt->ShiftX = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "ShiftY");
	if(o_tmp == 0) throw strEr_BadOptAng;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptShift;
	pOpt->ShiftY = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptZP*
 ***************************************************************************/
void ParseSructSRWLOptZP(SRWLOptZP* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "nZones");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->nZones = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "rn");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->rn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "thick");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->thick = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "delta1");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->delta1 = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "delta2");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->delta2 = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "atLen1");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->atLen1 = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "atLen2");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->atLen2 = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "x");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "y");
	if(o_tmp == 0) throw strEr_BadOptZP;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptZP;
	pOpt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptWG*
 ***************************************************************************/
void ParseSructSRWLOptWG(SRWLOptWG* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "L");
	if(o_tmp == 0) throw strEr_BadOptWG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptWG;
	pOpt->L = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Dx");
	if(o_tmp == 0) throw strEr_BadOptWG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptWG;
	pOpt->Dx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Dy");
	if(o_tmp == 0) throw strEr_BadOptWG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptWG;
	pOpt->Dy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "x");
	if(o_tmp == 0) throw strEr_BadOptWG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptWG;
	pOpt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "y");
	if(o_tmp == 0) throw strEr_BadOptWG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptWG;
	pOpt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptT*
 ***************************************************************************/
void ParseSructSRWLOptT(SRWLOptT* pOpt, PyObject* oOpt, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "arTr");
	//pOpt->arTr = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_SIMPLE);
	//pOpt->arTr = (double*)GetPyArrayBuf(o_tmp, vBuf, PyBUF_WRITABLE, 0);
	pOpt->arTr = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0);
	if(pOpt->arTr == 0) throw strEr_BadOptT;
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "nx");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->nx = PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "ny");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->ny = PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "rx");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->rx = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "ry");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->ry = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "mesh");
	if(o_tmp == 0) throw strEr_BadOptT;
	ParseSructSRWLRadMesh(&(pOpt->mesh), o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "extTr");
	if(o_tmp == 0) throw strEr_BadOptT;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	pOpt->extTr = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Fx");
	if(o_tmp == 0) throw strEr_BadOptT;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	pOpt->Fx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Fy");
	if(o_tmp == 0) throw strEr_BadOptT;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	pOpt->Fy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "x");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->x = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "y");
	//if(o_tmp == 0) throw strEr_BadOptT;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptT;
	//pOpt->y = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptMirEl*
 ***************************************************************************/
void ParseSructSRWLOptMir(SRWLOptMir* pOpt, PyObject* oOpt, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	//PyObject *o_tmp = PyObject_GetAttrString(oOpt, "arRefl");
	//pOpt->arRefl = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0);
	//if(pOpt->arRefl == 0) throw strEr_BadOptMir;
	//Py_DECREF(o_tmp);

	//OC12082018
	pOpt->arRefl = 0; //To allow not to process reflectivity if it is not defined
	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "arRefl");
	if(o_tmp != 0)
	{
		pOpt->arRefl = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0);
		Py_DECREF(o_tmp);
	}

	o_tmp = PyObject_GetAttrString(oOpt, "reflNumPhEn");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflNumPhEn = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflNumAng");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflNumAng = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflNumComp");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflNumComp = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflPhEnScaleType");
	if(o_tmp == 0) throw strEr_BadOptMir;
	CopyPyStringToC(o_tmp, pOpt->reflPhEnScaleType, 3);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflAngScaleType");
	if(o_tmp == 0) throw strEr_BadOptMir;
	CopyPyStringToC(o_tmp, pOpt->reflAngScaleType, 3);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflPhEnStart");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflPhEnStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflPhEnFin");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflPhEnFin = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflAngStart");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflAngStart = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "reflAngFin");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->reflAngFin = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "dt");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->dt = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "ds");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->ds = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "apShape");
	if(o_tmp == 0) throw strEr_BadOptMir;
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pOpt->apShape = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "meth");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->meth = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "npt");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->npt = (int)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nps");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->nps = (int)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "treatInOut");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->treatInOut = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "treatOut");
	//if(o_tmp == 0) throw strEr_BadOptMir;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	//pOpt->treatOut = (char)PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "extIn");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->extIn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "extOut");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->extOut = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvx");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->nvx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvy");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->nvy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvz");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->nvz = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "tvx");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->tvx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "tvy");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->tvy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "x");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "y");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Fx");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->Fx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "Fy");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->Fy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptMirPl*
 ***************************************************************************/
//void ParseSructSRWLOptMirPl(SRWLOptMirPl* pOpt, PyObject* oOpt, vector<Py_buffer>* pvBuf) //throw(...) 
//{
//	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
//
//	ParseSructSRWLOptMir(&(pOpt->baseMir), oOpt, pvBuf);
//}

/************************************************************************//**
 * Parses PyObject* to SRWLOptMirEl*
 ***************************************************************************/
void ParseSructSRWLOptMirExtEl(SRWLOptMirEl* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	//ParseSructSRWLOptMir(&(pOpt->baseMir), oOpt, pvBuf);

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "p");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->p = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "q");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->q = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "angGraz");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->angGraz = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "radSag");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->radSag = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptMirTor*
 ***************************************************************************/
void ParseSructSRWLOptMirExtTor(SRWLOptMirTor* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	//ParseSructSRWLOptMir(&(pOpt->baseMir), oOpt, pvBuf);

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "radTan");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->radTan = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "radSag");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->radSag = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptMirTor*
 ***************************************************************************/
void ParseSructSRWLOptMirExtSph(SRWLOptMirSph* pOpt, PyObject* oOpt) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	//ParseSructSRWLOptMir(&(pOpt->baseMir), oOpt, pvBuf);

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "rad");
	if(o_tmp == 0) throw strEr_BadOptMir;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptMir;
	pOpt->rad = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to a SRWLOptMir*
 ***************************************************************************/
void* ParseSructSRWLOptMirAll(PyObject* oOpt, char* sPyTypeName, vector<Py_buffer>* pvBuf, char* srwOptTypeName)
{
	if((oOpt == 0) || (srwOptTypeName == 0)) throw strEr_NoObj;

	bool needTypeName = false;
	if(sPyTypeName == 0) needTypeName = true; 
	else if(*sPyTypeName == '\0') needTypeName = true;

	char sPyTypeNameLoc[1025];
	if(needTypeName)
	{
		sPyTypeName = sPyTypeNameLoc;
		CopyPyClassNameToC(oOpt, sPyTypeName, 1024);
	}

	void *pMir = 0;
	strcpy(srwOptTypeName, "mirror: ");
	if(strcmp(sPyTypeName, "SRWLOptMirPl") == 0)
	{
		pMir = new SRWLOptMirPl();
		strcat(srwOptTypeName, "plane\0");
		ParseSructSRWLOptMir(&(((SRWLOptMirPl*)pMir)->baseMir), oOpt, pvBuf);
	}
	else if(strcmp(sPyTypeName, "SRWLOptMirEl") == 0)
	{
		pMir = new SRWLOptMirEl();
		strcat(srwOptTypeName, "ellipsoid\0");
		ParseSructSRWLOptMir(&(((SRWLOptMirEl*)pMir)->baseMir), oOpt, pvBuf);
		ParseSructSRWLOptMirExtEl((SRWLOptMirEl*)pMir, oOpt);
	}
	else if(strcmp(sPyTypeName, "SRWLOptMirTor") == 0)
	{
		pMir = new SRWLOptMirTor();
		strcat(srwOptTypeName, "toroid\0");
		ParseSructSRWLOptMir(&(((SRWLOptMirTor*)pMir)->baseMir), oOpt, pvBuf);
		ParseSructSRWLOptMirExtTor((SRWLOptMirTor*)pMir, oOpt);
	}
	else if (strcmp(sPyTypeName, "SRWLOptMirSph") == 0)
	{
		pMir = new SRWLOptMirSph();
		strcat(srwOptTypeName, "sphere\0");
		ParseSructSRWLOptMir(&(((SRWLOptMirSph*)pMir)->baseMir), oOpt, pvBuf);
		ParseSructSRWLOptMirExtSph((SRWLOptMirSph*)pMir, oOpt);
	}

	return pMir;
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptG*
 ***************************************************************************/
void ParseSructSRWLOptG(SRWLOptG* pOpt, PyObject* oOpt, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "mirSub");
	if(o_tmp == 0) throw strEr_BadOptG;
	pOpt->mirSub = ParseSructSRWLOptMirAll(o_tmp, 0, pvBuf, pOpt->mirSubType);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "m");
	if(o_tmp == 0) throw strEr_BadOptG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptG;
	pOpt->m = PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "grDen");
	if(o_tmp == 0) throw strEr_BadOptG;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptG;
	pOpt->grDen = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	pOpt->grDen1 = 0;
	o_tmp = PyObject_GetAttrString(oOpt, "grDen1");
	if(o_tmp != 0) 
	{
		if(PyNumber_Check(o_tmp)) 
		{
			pOpt->grDen1 = PyFloat_AsDouble(o_tmp);
			Py_DECREF(o_tmp);
		}
	}

	pOpt->grDen2 = 0;
	o_tmp = PyObject_GetAttrString(oOpt, "grDen2");
	if(o_tmp != 0) 
	{
		if(PyNumber_Check(o_tmp)) 
		{
			pOpt->grDen2 = PyFloat_AsDouble(o_tmp);
			Py_DECREF(o_tmp);
		}
	}

	pOpt->grDen3 = 0;
	o_tmp = PyObject_GetAttrString(oOpt, "grDen3");
	if(o_tmp != 0) 
	{
		if(PyNumber_Check(o_tmp)) 
		{
			pOpt->grDen3 = PyFloat_AsDouble(o_tmp);
			Py_DECREF(o_tmp);
		}
	}

	pOpt->grDen4 = 0;
	o_tmp = PyObject_GetAttrString(oOpt, "grDen4");
	if(o_tmp != 0) 
	{
		if(PyNumber_Check(o_tmp)) 
		{
			pOpt->grDen4 = PyFloat_AsDouble(o_tmp);
			Py_DECREF(o_tmp);
		}
	}

	pOpt->grAng = 0;
	o_tmp = PyObject_GetAttrString(oOpt, "grAng");
	if(o_tmp != 0) 
	{
		if(PyNumber_Check(o_tmp)) 
		{
			pOpt->grAng = PyFloat_AsDouble(o_tmp);
			Py_DECREF(o_tmp);
		}
	}
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptCryst*
 ***************************************************************************/
void ParseSructSRWLOptCryst(SRWLOptCryst* pOpt, PyObject* oOpt)
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;

	PyObject *o_tmp = PyObject_GetAttrString(oOpt, "dSp");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->dSp = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psi0r");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psi0r = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psi0i");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psi0i = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psiHr");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psiHr = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psiHi");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psiHi = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psiHbr");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psiHbr = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "psiHbi");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->psiHbi = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "h1"); //OC180314 (h1, h2, h3 removed)
	//if(o_tmp == 0) throw strEr_BadOptCryst;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	//pOpt->h1 = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "h2");
	//if(o_tmp == 0) throw strEr_BadOptCryst;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	//pOpt->h2 = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oOpt, "h3");
	//if(o_tmp == 0) throw strEr_BadOptCryst;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	//pOpt->h3 = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "tc");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->tc = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "angAs");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->angAs = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvx");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->nvx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvy");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->nvy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "nvz");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->nvz = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "tvx");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->tvx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "tvy");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->tvy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oOpt, "uc");
	if(o_tmp == 0) throw strEr_BadOptCryst;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadOptCryst;
	pOpt->uc = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLOptC*
 * ATTENTION: allocates void **arOpt and its elements and char **arOptTypes and double **arProp
 * vector<Py_buffer>& vBuf is required to store and release all buffers after the execution
 ***************************************************************************/
void ParseSructSRWLOptC(SRWLOptC* pOpt, PyObject* oOpt, vector<Py_buffer>* pvBuf) //throw(...) 
{
	if((pOpt == 0) || (oOpt == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oOpt, "arOpt");
	if(o_tmp == 0) throw strEr_BadOptC;
	if(!PyList_Check(o_tmp)) throw strEr_BadOptC;
	PyObject *o_List = o_tmp;
	//int nElem = (int)PyList_Size(o_tmp);
	int nElem = (int)PyList_Size(o_List);
	if(nElem <=  0) throw strEr_NoObj;

	pOpt->arPropN = 0; //OC031213

	o_tmp = PyObject_GetAttrString(oOpt, "arProp");
	if(o_tmp == 0) throw strEr_BadOptC;
	if(!PyList_Check(o_tmp)) throw strEr_BadOptC;
	int nElemProp = (int)PyList_Size(o_tmp);
	if(nElemProp > 0)
	{
		pOpt->arPropN = new char[nElemProp]; //OC031213
		pOpt->arProp = new double*[nElemProp];
		double **t_arProp = pOpt->arProp;
		for(int i=0; i<nElemProp; i++)
		{
			pOpt->arPropN[i] = 0;
			*t_arProp = 0;
			PyObject *o = PyList_GetItem(o_tmp, (Py_ssize_t)i);
			if(o != 0)
			{
				if(PyList_Check(o))
				{
					int nElemSub = 0;
					CopyPyListElemsToNumArray(o, 'd', *t_arProp, nElemSub);
					pOpt->arPropN[i] = (char)nElemSub;
				}
			}
			t_arProp++;
		}
	}
	Py_DECREF(o_tmp);
	pOpt->nProp = nElemProp;

	pOpt->arOpt = new void*[nElem];
	pOpt->arOptTypes = new char*[nElem];
	pOpt->nElem = 0; //in case if there will be reading error
	for(int i=0; i<nElem; i++)
	{
		PyObject *o = PyList_GetItem(o_List, (Py_ssize_t)i);

		//PyTypeObject *pTypeO = o->ob_type;
		//const char* sTypeName = pTypeO->tp_name;
		char sTypeName[1025];
		CopyPyClassNameToC(o, sTypeName, 1024);

		void *pOptElem=0;
		char *sOptType = new char[256];
		pOpt->arOptTypes[i] = 0;

		if(strcmp(sTypeName, "SRWLOptC") == 0)
		{
			pOptElem = new SRWLOptC();
			strcpy(sOptType, "container\0");
			ParseSructSRWLOptC((SRWLOptC*)pOptElem, o, pvBuf);
		}
		else if(strcmp(sTypeName, "SRWLOptD") == 0)
		{
			pOptElem = new SRWLOptD();
			strcpy(sOptType, "drift\0");
			ParseSructSRWLOptD((SRWLOptD*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptA") == 0)
		{
			pOptElem = new SRWLOptA();
			strcpy(sOptType, "aperture\0");
			ParseSructSRWLOptA((SRWLOptA*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptL") == 0)
		{
			pOptElem = new SRWLOptL();
			strcpy(sOptType, "lens\0");
			ParseSructSRWLOptL((SRWLOptL*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptAng") == 0)
		{
			pOptElem = new SRWLOptAng();
			strcpy(sOptType, "angle\0");
			ParseSructSRWLOptAng((SRWLOptAng*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptShift") == 0)
		{
			pOptElem = new SRWLOptShift();
			strcpy(sOptType, "shift\0");
			ParseSructSRWLOptShift((SRWLOptShift*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptZP") == 0)
		{
			pOptElem = new SRWLOptZP();
			strcpy(sOptType, "zp\0");
			ParseSructSRWLOptZP((SRWLOptZP*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptWG") == 0)
		{
			pOptElem = new SRWLOptWG();
			strcpy(sOptType, "waveguide\0");
			ParseSructSRWLOptWG((SRWLOptWG*)pOptElem, o);
		}
		else if(strcmp(sTypeName, "SRWLOptG") == 0)
		{
			pOptElem = new SRWLOptG();
			strcpy(sOptType, "grating\0");
			ParseSructSRWLOptG((SRWLOptG*)pOptElem, o, pvBuf);
		}
		else if(strcmp(sTypeName, "SRWLOptT") == 0)
		{
			pOptElem = new SRWLOptT();
			strcpy(sOptType, "transmission\0");
			ParseSructSRWLOptT((SRWLOptT*)pOptElem, o, pvBuf);
		}
		else if(strncmp(sTypeName, "SRWLOptMir", 10) == 0) pOptElem = ParseSructSRWLOptMirAll(o, sTypeName, pvBuf, sOptType);
		//else if(strcmp(sTypeName, "SRWLOptMirPl") == 0)
		//{
		//	pOptElem = new SRWLOptMirPl();
		//	strcpy(sOptType, "mirror: plane\0");
		//	ParseSructSRWLOptMirPl((SRWLOptMirPl*)pOptElem, o, pvBuf);
		//}
		//else if(strcmp(sTypeName, "SRWLOptMirEl") == 0)
		//{
		//	pOptElem = new SRWLOptMirEl();
		//	strcpy(sOptType, "mirror: ellipsoid\0");
		//	ParseSructSRWLOptMirEl((SRWLOptMirEl*)pOptElem, o, pvBuf);
		//}
		//else if(strcmp(sTypeName, "SRWLOptMirTor") == 0)
		//{
		//	pOptElem = new SRWLOptMirTor();
		//	strcpy(sOptType, "mirror: toroid\0");
		//	ParseSructSRWLOptMirTor((SRWLOptMirTor*)pOptElem, o, pvBuf);
		//}
		else if(strcmp(sTypeName, "SRWLOptCryst") == 0)
		{
			pOptElem = new SRWLOptCryst();
			strcpy(sOptType, "crystal\0");
			ParseSructSRWLOptCryst((SRWLOptCryst*)pOptElem, o);
		}

		//to add more optical elements here
		if(pOptElem != 0) 
		{
			pOpt->arOpt[i] = pOptElem;
			pOpt->arOptTypes[i] = sOptType;
			(pOpt->nElem)++; //in case if there will be reading error
		}
		else
		{
			if(sOptType != 0) delete[] sOptType;
		}
	}
	Py_DECREF(o_List);
}

/************************************************************************//**
 * Parses PyObject* to auxiliary arrays defining after which opt. elements
 * and what type of intensity distributions have to be calculated.
 * ATTENTION: it allocates arrays.
 ***************************************************************************/
int ParseSructSRWLPropIntDef(char** arIntDescr, SRWLRadMesh*& arIntMesh, PyObject* oInt)
{
	if(oInt == 0) throw strEr_NoObj;

	const int numIntDescr = 5;

	if(!PyList_Check(oInt)) throw strEr_BadList;
	int nElem = (int)PyList_Size(oInt);
	if(nElem <= 0) throw strEr_BadListIntProp;

	PyObject *o = PyList_GetItem(oInt, (Py_ssize_t)0);
	if(o == 0) throw strEr_BadListIntProp;

	PyObject *oSub = 0; //, *oSubListPrev = 0;
	bool assumeFlatListFormat = false;

	int nInt = 1;
	bool elemsMeanLists = false;
	if(PyList_Check(o)) 
	{
		nInt = (int)PyList_Size(o);
		if(nInt < 0) throw strEr_BadListIntProp;

		if(nInt == (numIntDescr + 1)) 
		{
			oSub = PyList_GetItem(o, (Py_ssize_t)numIntDescr);
			if(strcmp(oSub->ob_type->tp_name, "SRWLRadMesh") == 0) assumeFlatListFormat = true;
		}

		elemsMeanLists = true;
	}

	if(assumeFlatListFormat)
	{//Parses this format:
	 //[[4, 6, 0, 3, 0, SRWLRadMesh(_xStart=0, _yStart=0)],
     //[11, 6, 0, 3, 0, SRWLRadMesh(_xStart=0, _yStart=0)],
     //[11, 6, 5, 3, 0, SRWLRadMesh(_xStart=0, _yStart=0)],
     //[11, 6, 6, 3, 0, SRWLRadMesh(_xStart=0, _yStart=0)]]

		for(int j=0; j<numIntDescr; j++) arIntDescr[j] = new char[nElem];
		arIntMesh = new SRWLRadMesh[nElem];

		for(int i=0; i<nElem; i++)
		{
			o = PyList_GetItem(oInt, (Py_ssize_t)i);
			if(o == 0) throw strEr_BadListIntProp;

			for(int j=0; j<=numIntDescr; j++)
			{
				oSub = PyList_GetItem(o, (Py_ssize_t)j);
				if(oSub == 0) throw strEr_BadListIntProp;

				if(j < numIntDescr)
				{
					if(!PyNumber_Check(oSub)) throw strEr_BadListIntProp;
					char *pIntDescr = arIntDescr[j];
					pIntDescr[i] = (char)PyLong_AsLong(oSub);
				}
				else
				{
					ParseSructSRWLRadMesh(arIntMesh + i, oSub);
  				}
			}
		}
		return nElem;
	}
		
	for(int j=0; j<numIntDescr; j++) arIntDescr[j] = new char[nInt];
	arIntMesh = new SRWLRadMesh[nInt];

	//SRWLRadMesh *tMesh = arIntMesh;

	for(int j=0; j<=numIntDescr; j++)
	{
		char *pIntDescr = arIntDescr[j];
		o = PyList_GetItem(oInt, (Py_ssize_t)j);
		if(o == 0) throw strEr_BadListIntProp;

		bool curElemIsList = elemsMeanLists;
		if(elemsMeanLists)
		{
			if(!PyList_Check(o)) 
			{
				if(j == 0) throw strEr_BadListIntProp;
				else curElemIsList = false;
			}
		}
		if(curElemIsList)
		{
			for(int i=0; i<nInt; i++)
			{
				oSub = PyList_GetItem(o, (Py_ssize_t)i);
				if(j < numIntDescr)
				{
					pIntDescr[i] = (char)PyLong_AsLong(oSub);
				}
				else
				{
					ParseSructSRWLRadMesh(arIntMesh + i, oSub);
  				}
			}
		}
		else
		{
			if(j < numIntDescr)
			{
				if(!PyNumber_Check(o)) throw strEr_BadListIntProp;
				*pIntDescr = (char)PyLong_AsLong(o);
				if(nInt > 1)
				{
					*(pIntDescr + 1) = -1; //to mark that only the first element of the array is defined
				}
			}
			else
			{
				ParseSructSRWLRadMesh(arIntMesh, o);
				if(nInt > 1)
				{
					(arIntMesh + 1)->ne = -1; //to mark that only the first element of the array is defined
				}
  			}
		}
	}
	return nInt;
}

/************************************************************************//**
 * Parses PyObject* to SRWLGsnBm*
 ***************************************************************************/
void ParseSructSRWLGsnBm(SRWLGsnBm* pGsnBm, PyObject* oGsnBm) //throw(...) 
{
	if((pGsnBm == 0) || (oGsnBm == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oGsnBm, "x");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "y");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "z");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->z = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "xp");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->xp = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "yp");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->yp = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "avgPhotEn");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->avgPhotEn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "pulseEn");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->pulseEn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "repRate");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->repRate = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oGsnBm, "repRate");
	//if(o_tmp == 0) throw strEr_BadGsnBm;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	//pGsnBm->repRate = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "polar");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->polar = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "sigX");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->sigX = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "sigY");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->sigY = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "sigT");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->sigT = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "mx");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->mx = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oGsnBm, "my");
	if(o_tmp == 0) throw strEr_BadGsnBm;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadGsnBm;
	pGsnBm->my = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLPtSrc*
 ***************************************************************************/
void ParseSructSRWLPtSrc(SRWLPtSrc* pPtSrc, PyObject* oPtSrc) //throw(...) 
{
	if((pPtSrc == 0) || (oPtSrc == 0)) throw strEr_NoObj;
	PyObject *o_tmp = 0;

	o_tmp = PyObject_GetAttrString(oPtSrc, "x");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->x = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPtSrc, "y");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->y = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPtSrc, "z");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->z = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPtSrc, "flux");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->flux = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPtSrc, "unitFlux");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->unitFlux = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oPtSrc, "polar");
	if(o_tmp == 0) throw strEr_BadPtSrc;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadPtSrc;
	pPtSrc->polar = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Parses PyObject* to SRWLWfr*
 * vector<Py_buffer>& vBuf is required to release all buffers after the end of execution
 ***************************************************************************/
void ParseSructSRWLWfr(SRWLWfr* pWfr, PyObject* oWfr, vector<Py_buffer>* pvBuf, map<SRWLWfr*, AuxStructPyObjectPtrs>& mWfrPyPtr)
{
	if((pWfr == 0) || (oWfr == 0)) throw strEr_NoObj;

	AuxStructPyObjectPtrs sPyObjectPtrs;
	sPyObjectPtrs.o_wfr = oWfr;
	sPyObjectPtrs.pv_buf = pvBuf;

	//Py_INCREF(oWfr); //??

	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	pWfr->arEx = pWfr->arEy = 0;
	o_tmp = PyObject_GetAttrString(oWfr, "arEx");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadWfr;
	//	pvBuf->push_back(pb_tmp);
	//	pWfr->arEx = (char*)pb_tmp.buf;
	//	sPyObjectPtrs.pbEx = pb_tmp;
	//}
	//else pWfr->arEx = 0;
	int sizeVectBuf = (int)pvBuf->size();
	if(!(pWfr->arEx = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbEx = (*pvBuf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "arEy");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) throw strEr_BadWfr;
	//	pvBuf->push_back(pb_tmp);
	//	pWfr->arEy = (char*)pb_tmp.buf;
	//	sPyObjectPtrs.pbEy = pb_tmp;
	//}
	//else pWfr->arEy = 0;
	sizeVectBuf = (int)pvBuf->size();
	if(!(pWfr->arEy = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbEy = (*pvBuf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	if((pWfr->arEx == 0) && (pWfr->arEy == 0)) throw strEr_BadWfr;

	pWfr->arExAux = 0; //OC161115
	if(PyObject_HasAttrString(oWfr, "arExAux"))
	{
		o_tmp = PyObject_GetAttrString(oWfr, "arExAux");
		if(o_tmp != 0)
		{
			sizeVectBuf = (int)pvBuf->size();
			//if(pWfr->arExAux = GetPyArrayBuf(o_tmp, pvBuf, 0))
			pWfr->arExAux = GetPyArrayBuf(o_tmp, pvBuf, 0);
			if(pWfr->arExAux)
			{
				if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbExAux = (*pvBuf)[sizeVectBuf];
				Py_DECREF(o_tmp);
			}
		}
	}

	pWfr->arEyAux = 0; //OC161115
	if(PyObject_HasAttrString(oWfr, "arEyAux"))
	{
		o_tmp = PyObject_GetAttrString(oWfr, "arEyAux");
		if(o_tmp != 0)
		{
			sizeVectBuf = (int)pvBuf->size();
			//if(pWfr->arEyAux = GetPyArrayBuf(o_tmp, pvBuf, 0))
			pWfr->arEyAux = GetPyArrayBuf(o_tmp, pvBuf, 0);
			if(pWfr->arEyAux)
			{
				if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbEyAux = (*pvBuf)[sizeVectBuf];
				Py_DECREF(o_tmp);
			}
		}
	}

	//o_tmp = PyObject_GetAttrString(oWfr, "eStart");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->eStart = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "eFin");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->eFin = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "xStart");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->xStart = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "xFin");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->xFin = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "yStart");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->yStart = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "yFin");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->yFin = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "zStart");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->zStart = PyFloat_AsDouble(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "ne");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->ne = PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "nx");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->nx = PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oWfr, "ny");
	//if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	//pWfr->ny = PyLong_AsLong(o_tmp);
	//Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "mesh");
	if(o_tmp == 0) throw strEr_BadWfr;
	ParseSructSRWLRadMesh(&(pWfr->mesh), o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "Rx");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->Rx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "Ry");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->Ry = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "dRx");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->dRx = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "dRy");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->dRy = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "xc");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->xc = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "yc");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->yc = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "avgPhotEn");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->avgPhotEn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "presCA");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->presCA = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "presFT");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->presFT = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "numTypeElFld");
	if(o_tmp == 0) throw strEr_BadWfr;
	//PyObject *o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(!PyBytes_Check(o_str)) throw strEr_BadWfr;
	//pWfr->numTypeElFld = *PyBytes_AsString(o_str);
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pWfr->numTypeElFld = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "unitElFld");
	if(o_tmp == 0) throw strEr_BadWfr;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
	pWfr->unitElFld = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	pWfr->unitElFldAng = 0; //OC19112017
	if(PyObject_HasAttrString(oWfr, "unitElFldAng"))
	{
		o_tmp = PyObject_GetAttrString(oWfr, "unitElFldAng");
		if(o_tmp == 0) throw strEr_BadWfr;
		if(!PyNumber_Check(o_tmp)) throw strEr_BadWfr;
		pWfr->unitElFldAng = (char)PyLong_AsLong(o_tmp);
		Py_DECREF(o_tmp);
	}

	o_tmp = PyObject_GetAttrString(oWfr, "partBeam");
	if(o_tmp == 0) throw strEr_BadWfr;
	ParseSructSRWLPartBeam(&(pWfr->partBeam), o_tmp, *pvBuf);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "arElecPropMatr");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadWfr;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadWfr;
	//pvBuf->push_back(pb_tmp);
	//pWfr->arElecPropMatr = (double*)pb_tmp.buf;
	if(!(pWfr->arElecPropMatr = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "arMomX");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadWfr;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadWfr;
	//pvBuf->push_back(pb_tmp);
	//pWfr->arMomX = (double*)pb_tmp.buf;
	//sPyObjectPtrs.pbMomX = pb_tmp;
	sizeVectBuf = (int)pvBuf->size();
	if(!(pWfr->arMomX = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbMomX = (*pvBuf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "arMomY");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadWfr;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadWfr;
	//pvBuf->push_back(pb_tmp);
	//pWfr->arMomY = (double*)pb_tmp.buf;
	//sPyObjectPtrs.pbMomY = pb_tmp;
	sizeVectBuf = (int)pvBuf->size();
	if(!(pWfr->arMomY = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	if((int)pvBuf->size() > sizeVectBuf) sPyObjectPtrs.pbMomY = (*pvBuf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oWfr, "arWfrAuxData");
	if(o_tmp == 0) throw strEr_BadWfr;
	//if(!PyObject_CheckBuffer(o_tmp)) throw strEr_BadWfr;
	//if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_SIMPLE)) throw strEr_BadWfr;
	//pvBuf->push_back(pb_tmp);
	//pWfr->arWfrAuxData = (double*)pb_tmp.buf;
	if(!(pWfr->arWfrAuxData = (double*)GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadWfr;
	Py_DECREF(o_tmp);

	mWfrPyPtr[pWfr] = sPyObjectPtrs;
}

/************************************************************************//**
 * Parses PyObject* to SRWLStokes*
 * vector<Py_buffer>& vBuf is required to release all buffers after the end of execution
 ***************************************************************************/
void ParseSructSRWLStokes(SRWLStokes* pStokes, PyObject* oStokes, vector<Py_buffer>* pvBuf)
{
	if((pStokes == 0) || (oStokes == 0)) throw strEr_NoObj;

	//Py_INCREF(oStokes); //??
	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;

	pStokes->arS0 = pStokes->arS1 = pStokes->arS2 = pStokes->arS3 = 0;

	o_tmp = PyObject_GetAttrString(oStokes, "arS");
	//o_tmp = PyObject_GetAttrString(oStokes, "arS0");
	if(o_tmp == 0) throw strEr_BadStokes;
	if(!(pStokes->arS0 = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadStokes;
	Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oStokes, "arS1");
	//if(o_tmp == 0) throw strEr_BadStokes;
	//if(!(pStokes->arS1 = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadStokes;
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oStokes, "arS2");
	//if(o_tmp == 0) throw strEr_BadStokes;
	//if(!(pStokes->arS2 = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadStokes;
	//Py_DECREF(o_tmp);

	//o_tmp = PyObject_GetAttrString(oStokes, "arS3");
	//if(o_tmp == 0) throw strEr_BadStokes;
	//if(!(pStokes->arS3 = GetPyArrayBuf(o_tmp, pvBuf, 0))) throw strEr_BadStokes;
	//Py_DECREF(o_tmp);

	//if((pStokes->arS0 == 0) && (pStokes->arS1 == 0) && (pStokes->arS2 == 0) && (pStokes->arS3 == 0)) throw strEr_BadStokes;

	o_tmp = PyObject_GetAttrString(oStokes, "mesh");
	if(o_tmp == 0) throw strEr_BadStokes;
	//ParseSructSRWLRadMesh(&(pStokes->mesh), o_tmp);
	ParseSructSRWLRadMesh(&(pStokes->mesh), o_tmp, pvBuf);
	Py_DECREF(o_tmp);

	//long npTotS0 = (pStokes->mesh.ne)*(pStokes->mesh.nx)*(pStokes->mesh.ny);
	long long npTotS0 = ((long long)(pStokes->mesh.ne))*((long long)(pStokes->mesh.nx))*((long long)(pStokes->mesh.ny));
	pStokes->arS1 = (char*)(((float*)(pStokes->arS0)) + npTotS0);
	pStokes->arS2 = (char*)(((float*)(pStokes->arS1)) + npTotS0);
	pStokes->arS3 = (char*)(((float*)(pStokes->arS2)) + npTotS0);

	o_tmp = PyObject_GetAttrString(oStokes, "avgPhotEn");
	if(o_tmp == 0) throw strEr_BadStokes;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
	pStokes->avgPhotEn = PyFloat_AsDouble(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oStokes, "presCA");
	if(o_tmp == 0) throw strEr_BadStokes;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
	pStokes->presCA = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oStokes, "presFT");
	if(o_tmp == 0) throw strEr_BadStokes;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
	pStokes->presFT = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oStokes, "numTypeStokes");
	if(o_tmp == 0) throw strEr_BadStokes;
	//PyObject *o_str = PyUnicode_AsUTF8String(o_tmp);
	//if(!PyBytes_Check(o_str)) throw strEr_BadStokes;
	//pStokes->numTypeStokes = *PyBytes_AsString(o_str);
	char cStrBuf[2];
	CopyPyStringToC(o_tmp, cStrBuf, 1);
	pStokes->numTypeStokes = *cStrBuf;
	Py_DECREF(o_tmp);

	o_tmp = PyObject_GetAttrString(oStokes, "unitStokes");
	if(o_tmp == 0) throw strEr_BadStokes;
	if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
	pStokes->unitStokes = (char)PyLong_AsLong(o_tmp);
	Py_DECREF(o_tmp);
}

/************************************************************************//**
 * Updates Py List by numbers
 ***************************************************************************/
template<class T> void UpdatePyListNum(PyObject* oList, const T* ar, int n) //OC03092016
{
	if((ar == 0) || (n <= 0)) return;

	if(!PyList_Check(oList)) throw strEr_BadList;
	int nElem = (int)PyList_Size(oList);

	int nExist = n;
	if(nExist > nElem) nExist = nElem;

	for(int i=0; i<nExist; i++)
	{
		PyObject *oElemOld = PyList_GetItem(oList, (Py_ssize_t)i);
		if(oElemOld == 0) throw strEr_BadNum;
		if(PyNumber_Check(oElemOld) != 1) throw strEr_BadNum;

		char arNumType[2];
		arNumType[1] = '\0';
		//if(PyLong_Check(oElemOld)) arNumType[0] = 'i';
		//else if(PyFloat_Check(oElemOld))  arNumType[0] = 'd';
		//if(PyList_SetItem(oList, (Py_ssize_t)i, Py_BuildValue(arNumType, ar[i])) != 0) throw strEr_BadNum;
		//OC02112017

#if PY_MAJOR_VERSION >= 3 //OC21112017
		if(PyLong_Check(oElemOld)) //This compiles, but desn't make a correct test under Py 2.x
#else
		if(PyInt_Check(oElemOld)) //This doesn't compile with Py 3.x
#endif
		{
			arNumType[0] = 'i';
			if(PyList_SetItem(oList, (Py_ssize_t)i, Py_BuildValue(arNumType, (int)(ar[i]))) != 0) throw strEr_BadNum;
		}
		else if(PyFloat_Check(oElemOld))  
		{
			arNumType[0] = 'd';
			if(PyList_SetItem(oList, (Py_ssize_t)i, Py_BuildValue(arNumType, (double)(ar[i]))) != 0) throw strEr_BadNum;
		}
	}

	for(int j=nExist; j<n; j++) //??
	{
		if(PyList_Append(oList, Py_BuildValue("d", ar[j])) != 0) throw strEr_BadNum; //consider analyzing T ?
	}
}

/************************************************************************//**
 * Updates Electric Field Wavefront structure in Py (excluding buffers)
 ***************************************************************************/
void UpdatePyRadMesh(PyObject* oRadMesh, SRWLRadMesh* pMesh)
{//OC19082018
	if((pMesh == 0) || (oRadMesh == 0)) throw strEr_NoObj;

	if(PyObject_SetAttrString(oRadMesh, "eStart", Py_BuildValue("d", pMesh->eStart))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "eFin", Py_BuildValue("d", pMesh->eFin))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "xStart", Py_BuildValue("d", pMesh->xStart))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "xFin", Py_BuildValue("d", pMesh->xFin))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "yStart", Py_BuildValue("d", pMesh->yStart))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "yFin", Py_BuildValue("d", pMesh->yFin))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "zStart", Py_BuildValue("d", pMesh->zStart))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "ne", Py_BuildValue("i", pMesh->ne))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "nx", Py_BuildValue("i", pMesh->nx))) throw strEr_BadRadMesh;
	if(PyObject_SetAttrString(oRadMesh, "ny", Py_BuildValue("i", pMesh->ny))) throw strEr_BadRadMesh;
	//Add more member updates if necessary
}

/************************************************************************//**
 * Updates Electric Field Wavefront structure in Py (excluding buffers)
 ***************************************************************************/
void UpdatePyWfr(PyObject* oWfr, SRWLWfr* pWfr)
{
	if((pWfr == 0) || (oWfr == 0)) throw strEr_NoObj;

	//if(PyObject_SetAttrString(oWfr, "eStart", Py_BuildValue("d", pWfr->eStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "eFin", Py_BuildValue("d", pWfr->eFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "xStart", Py_BuildValue("d", pWfr->xStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "xFin", Py_BuildValue("d", pWfr->xFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "yStart", Py_BuildValue("d", pWfr->yStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "yFin", Py_BuildValue("d", pWfr->yFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "zStart", Py_BuildValue("d", pWfr->zStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "ne", Py_BuildValue("i", pWfr->ne))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "nx", Py_BuildValue("i", pWfr->nx))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "ny", Py_BuildValue("i", pWfr->ny))) throw strEr_BadWfr;

	PyObject *oRadMesh = PyObject_GetAttrString(oWfr, "mesh");
	if(oRadMesh == 0) throw strEr_BadWfr;
	//SRWLRadMesh &mesh = pWfr->mesh;
	//if(PyObject_SetAttrString(oRadMesh, "eStart", Py_BuildValue("d", mesh.eStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "eFin", Py_BuildValue("d", mesh.eFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "xStart", Py_BuildValue("d", mesh.xStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "xFin", Py_BuildValue("d", mesh.xFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "yStart", Py_BuildValue("d", mesh.yStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "yFin", Py_BuildValue("d", mesh.yFin))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "zStart", Py_BuildValue("d", mesh.zStart))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "ne", Py_BuildValue("i", mesh.ne))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "nx", Py_BuildValue("i", mesh.nx))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oRadMesh, "ny", Py_BuildValue("i", mesh.ny))) throw strEr_BadWfr;
	UpdatePyRadMesh(oRadMesh, &(pWfr->mesh)); //OC19082018
	Py_DECREF(oRadMesh);

	if(PyObject_SetAttrString(oWfr, "Rx", Py_BuildValue("d", pWfr->Rx))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "Ry", Py_BuildValue("d", pWfr->Ry))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "dRx", Py_BuildValue("d", pWfr->dRx))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "dRy", Py_BuildValue("d", pWfr->dRy))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "xc", Py_BuildValue("d", pWfr->xc))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "yc", Py_BuildValue("d", pWfr->yc))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "avgPhotEn", Py_BuildValue("d", pWfr->avgPhotEn))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "presCA", Py_BuildValue("i", pWfr->presCA))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "presFT", Py_BuildValue("i", pWfr->presFT))) throw strEr_BadWfr;
	//if(PyObject_SetAttrString(oWfr, "numTypeElFld", Py_BuildValue("C", pWfr->numTypeElFld))) throw strEr_BadWfr; //doesn't work with Py2.7
	if(PyObject_SetAttrString(oWfr, "numTypeElFld", Py_BuildValue("c", pWfr->numTypeElFld))) throw strEr_BadWfr;
	if(PyObject_SetAttrString(oWfr, "unitElFld", Py_BuildValue("i", pWfr->unitElFld))) throw strEr_BadWfr;
}

/************************************************************************//**
 * Updates Stokes parameters structure in Py (excluding buffers)
 ***************************************************************************/
void UpdatePyStokes(PyObject* oStk, SRWLStokes* pStk)
{
	if((pStk == 0) || (oStk == 0)) throw strEr_NoObj;

	//if(PyObject_SetAttrString(oStk, "eStart", Py_BuildValue("d", pStk->eStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "eFin", Py_BuildValue("d", pStk->eFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "xStart", Py_BuildValue("d", pStk->xStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "xFin", Py_BuildValue("d", pStk->xFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "yStart", Py_BuildValue("d", pStk->yStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "yFin", Py_BuildValue("d", pStk->yFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "zStart", Py_BuildValue("d", pStk->zStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "ne", Py_BuildValue("i", pStk->ne))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "nx", Py_BuildValue("i", pStk->nx))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "ny", Py_BuildValue("i", pStk->ny))) throw strEr_BadStokes;

	PyObject *oRadMesh = PyObject_GetAttrString(oStk, "mesh");
	if(oRadMesh == 0) throw strEr_BadStokes;
	//SRWLRadMesh &mesh = pStk->mesh;
	//if(PyObject_SetAttrString(oRadMesh, "eStart", Py_BuildValue("d", mesh.eStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "eFin", Py_BuildValue("d", mesh.eFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "xStart", Py_BuildValue("d", mesh.xStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "xFin", Py_BuildValue("d", mesh.xFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "yStart", Py_BuildValue("d", mesh.yStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "yFin", Py_BuildValue("d", mesh.yFin))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "zStart", Py_BuildValue("d", mesh.zStart))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "ne", Py_BuildValue("i", mesh.ne))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "nx", Py_BuildValue("i", mesh.nx))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oRadMesh, "ny", Py_BuildValue("i", mesh.ny))) throw strEr_BadStokes;
	UpdatePyRadMesh(oRadMesh, &(pStk->mesh)); //OC19082018
	Py_DECREF(oRadMesh);

	if(PyObject_SetAttrString(oStk, "avgPhotEn", Py_BuildValue("d", pStk->avgPhotEn))) throw strEr_BadStokes;
	if(PyObject_SetAttrString(oStk, "presCA", Py_BuildValue("i", pStk->presCA))) throw strEr_BadStokes;
	if(PyObject_SetAttrString(oStk, "presFT", Py_BuildValue("i", pStk->presFT))) throw strEr_BadStokes;
	//if(PyObject_SetAttrString(oStk, "numTypeStokes", Py_BuildValue("C", pStk->numTypeStokes))) throw strEr_BadStokes; //doesn't work with Py2.7
	//if(PyObject_SetAttrString(oStk, "numTypeStokes", Py_BuildValue("c", pStk->numTypeStokes))) throw strEr_BadStokes;
	char sNumTypeStokes[] = {pStk->numTypeStokes, '\0'}; //OC060714 (because the above generates a byte, not a string, in Py)
	if(PyObject_SetAttrString(oStk, "numTypeStokes", Py_BuildValue("s", sNumTypeStokes))) throw strEr_BadStokes;
	if(PyObject_SetAttrString(oStk, "unitStokes", Py_BuildValue("i", pStk->unitStokes))) throw strEr_BadStokes;
}

/************************************************************************//**
 * Updates Py magnetic field Undulator Harmonic structure (primary use: at converting tabulated field to periodic)
 ***************************************************************************/
void UpdatePyMagFldH(PyObject* oMagFldH, SRWLMagFldH* pMagFldH)
{
	if(oMagFldH == 0) throw strEr_NoObj;

	SRWLMagFldH HarmZero = {0,0,0,0,0,0};
	if(pMagFldH == 0) pMagFldH = &HarmZero;

	if(PyObject_SetAttrString(oMagFldH, "n", Py_BuildValue("i", pMagFldH->n))) throw strEr_BadMagH;
	if(PyObject_SetAttrString(oMagFldH, "h_or_v", Py_BuildValueChar(pMagFldH->h_or_v))) throw strEr_BadMagH; //to check with Py 2.7 & 3.x
	if(PyObject_SetAttrString(oMagFldH, "B", Py_BuildValue("d", pMagFldH->B))) throw strEr_BadMagH;
	if(PyObject_SetAttrString(oMagFldH, "ph", Py_BuildValue("d", pMagFldH->ph))) throw strEr_BadMagH;
	if(PyObject_SetAttrString(oMagFldH, "s", Py_BuildValue("i", pMagFldH->s))) throw strEr_BadMagH;
	if(PyObject_SetAttrString(oMagFldH, "a", Py_BuildValue("d", pMagFldH->a))) throw strEr_BadMagH;
}

/************************************************************************//**
 * Updates Py magnetic field Undulator structure (primary use: at converting tabulated field to periodic)
 ***************************************************************************/
void UpdatePyMagFldU(PyObject* oMagFldU, SRWLMagFldU* pMagFldU)
{
	if((oMagFldU == 0) || (pMagFldU == 0)) throw strEr_NoObj;

	if(PyObject_SetAttrString(oMagFldU, "per", Py_BuildValue("d", pMagFldU->per))) throw strEr_BadMagU;
	if(PyObject_SetAttrString(oMagFldU, "nPer", Py_BuildValue("i", pMagFldU->nPer))) throw strEr_BadMagU;

	//NOTE: this doesn't modify the number of harmonics in oMagFldU! (to be updated)

	PyObject* o_ListHarm = PyObject_GetAttrString(oMagFldU, "arHarm");
	if(o_ListHarm == 0) throw strEr_BadMagU;
	if(!PyList_Check(o_ListHarm)) throw strEr_BadMagU;
	int nHarmExt = (int)PyList_Size(o_ListHarm);
	if(nHarmExt <=  0) throw strEr_NoObj;

	SRWLMagFldH *pHarm;
	for(int i=0; i<nHarmExt; i++)
	{
		PyObject *oHarm = PyList_GetItem(o_ListHarm, (Py_ssize_t)i);
		pHarm = 0;
		if(i < pMagFldU->nHarm) pHarm = (pMagFldU->arHarm) + i;

		if(pHarm == 0) break;
		UpdatePyMagFldH(oHarm, pHarm);

		//if(pHarm != 0)
		//{
		//	if(PyObject_SetAttrString(oHarm, "B", Py_BuildValue("d", pHarm->B))) throw strEr_BadMagH;
		//}

		//oHarm = PyList_GetItem(o_ListHarm, (Py_ssize_t)i);
		//PyObject* o_tmp = PyObject_GetAttrString(oHarm, "B");
		//if(o_tmp == 0) throw strEr_BadStokes;
		//if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
		//double testB = PyFloat_AsDouble(o_tmp);
		//int aha = 1;

		////if(PyList_SetItem(o_ListHarm, (Py_ssize_t)i, oHarm)) throw strEr_BadMagU;
		//Py_XINCREF(oHarm);
	}

		//o_ListHarm = PyObject_GetAttrString(oMagFldU, "arHarm");
		//PyObject *oHarm = PyList_GetItem(o_ListHarm, (Py_ssize_t)0);
		//PyObject *o_tmp = PyObject_GetAttrString(oHarm, "B");
		//if(o_tmp == 0) throw strEr_BadStokes;
		//if(!PyNumber_Check(o_tmp)) throw strEr_BadStokes;
		//double testB = PyFloat_AsDouble(o_tmp);
		//Py_DECREF(o_tmp);

	Py_DECREF(o_ListHarm);
}

/************************************************************************//**
 * Updates Py magnetic field Container structure (primary use: at converting tabulated field to periodic)
 ***************************************************************************/
void UpdatePyMagFldC(PyObject* oMagFldC, SRWLMagFldC* pMagFldC)
{
	if((oMagFldC == 0) || (pMagFldC == 0)) throw strEr_NoObj;

	PyObject *o_List = PyObject_GetAttrString(oMagFldC, "arMagFld");
	if(o_List == 0) throw strEr_BadMagC;
	if(!PyList_Check(o_List)) throw strEr_BadMagC;

	int nElem = (int)PyList_Size(o_List);
	if(nElem <=  0) throw strEr_NoObj;

	for(int i=0; i<nElem; i++)
	{
		PyObject *o = PyList_GetItem(o_List, (Py_ssize_t)i);
		char cFldType = pMagFldC->arMagFldTypes[i];
		void *pvMagFld = pMagFldC->arMagFld[i];

		if(cFldType == 'c') //SRWLMagFldC
		{
			UpdatePyMagFldC(o, (SRWLMagFldC*)pvMagFld);
		}
		else if(cFldType == 'a') //SRWLMagFld3D
		{//to implement
			//UpdatePyMagFld3D(o, (SRWLMagFld3D*)pvMagFld);
		}
		else if(cFldType == 'm') //SRWLMagFldM
		{//to implement
			//UpdatePyMagFldM(o, (SRWLMagFldM*)pvMagFld);
		}
		else if(cFldType == 's') //SRWLMagFldS
		{//to implement
			//UpdatePyMagFldS(o, (SRWLMagFldS*)pvMagFld);
		}
		else if(cFldType == 'u') //SRWLMagFldU
		{//to implement
			UpdatePyMagFldU(o, (SRWLMagFldU*)pvMagFld);
		}
		//to add more magnetic elements
	}
	Py_DECREF(o_List);
}

/************************************************************************//**
 * Updates propagated intensity and corresponding meshes in Py 
 ***************************************************************************/
void UpdatePyPropInt(PyObject* oInt, SRWLRadMesh* arIntMesh, char** arInts, int nInt) //OC19082018
{
	if((oInt == 0) || (arIntMesh == 0) || (arInts == 0)) return;

	const int indMesh = 5; //keep updated
	const int indInt = indMesh + 1; //keep updated

	if(!PyList_Check(oInt)) throw strEr_BadListIntProp;
	int nElem = (int)PyList_Size(oInt);
	//if(nElem <= indMesh) throw strEr_BadListIntProp;
	if(nElem <= 0) throw strEr_BadListIntProp; //29082018

	PyObject *o = PyList_GetItem(oInt, (Py_ssize_t)0);
	if(o == 0) throw strEr_BadListIntProp;

	PyObject *oMesh=0, *oAr=0;
	SRWLRadMesh *t_arIntMesh = arIntMesh;

	bool assumeFlatListFormat = false;
	int nSub = 0;
	if(PyList_Check(o)) 
	{
		nSub = (int)PyList_Size(o);
		if(nSub < 0) throw strEr_BadListIntProp;
		if(nSub >= indInt) 
		{
			PyObject *oSub = PyList_GetItem(o, (Py_ssize_t)indMesh);
			if(strcmp(oSub->ob_type->tp_name, "SRWLRadMesh") == 0) assumeFlatListFormat = true;
		}
	}

	if(assumeFlatListFormat)
	{//Returnes this format:
	 //[[4, 6, 0, 3, 0, SRWLRadMesh(..), arI],
     //[11, 6, 0, 3, 0, SRWLRadMesh(..), arI],
     //[11, 6, 5, 3, 0, SRWLRadMesh(..), arI],
     //[11, 6, 6, 3, 0, SRWLRadMesh(..), arI]]

		for(int i=0; i<nElem; i++)
		{
			o = PyList_GetItem(oInt, (Py_ssize_t)i);
			if(o == 0) throw strEr_BadListIntProp;

			oMesh = PyList_GetItem(o, (Py_ssize_t)indMesh);
			if(oMesh == 0) throw strEr_BadListIntProp;

			UpdatePyRadMesh(oMesh, t_arIntMesh++);

			//Updating Intensity data
			char *pCurInt = arInts[i];
			if(pCurInt != 0)
			{
				map<char*, PyObject*>::iterator it = gmBufPyObjPtr.find(pCurInt);
				if(it == gmBufPyObjPtr.end()) throw strEr_FailedUpdateInt;
				oAr = it->second;
				if(oAr == 0) throw strEr_FailedUpdateInt;
			}
			if(oAr != 0)
			{
				if(nSub > indInt)
				{
					if(PyList_SetItem(o, (Py_ssize_t)indInt, oAr)) throw strEr_FailedUpdateInt;
				}
				else
				{
					if(PyList_Append(o, oAr)) throw strEr_FailedUpdateInt;
				}
			}
		}
		return;
	}

	oMesh = PyList_GetItem(oInt, (Py_ssize_t)indMesh);
	if(oMesh == 0) throw strEr_BadListIntProp;

	PyObject *oMeshTrue=0; 
	PyObject *oNewListMesh = 0;
	bool elemsAreLists = false;
	bool meshNeedsToBeCopied = false;
	bool meshListNeedsToBeSet = false;
	if(nInt > 1)
	{
		if(!PyList_Check(oMesh)) 
		{
			oNewListMesh = PyList_New((Py_ssize_t)nInt);
			if(oNewListMesh == 0) throw strEr_FailedUpdateInt;

			if(PyList_SetItem(oNewListMesh, (Py_ssize_t)0, oMesh)) throw strEr_FailedUpdateInt;
			oMeshTrue = oMesh;
			oMesh = oNewListMesh;

			meshNeedsToBeCopied = true;
			meshListNeedsToBeSet = true;
		}

		//if((int)PyList_Size(oMesh) < nInt) throw strEr_BadListIntProp;
		elemsAreLists = true;
	}
	else
	{
		if(PyList_Check(oMesh)) elemsAreLists = true;
		//else if(!PyNumber_Check(oMesh)) throw strEr_BadListIntProp;
	}

	PyObject *oM=0, *oIntAr=0;
	PyObject *oFunc=0, *argList=0;

	if(elemsAreLists)
	{
		oIntAr = PyList_New((Py_ssize_t)nInt);
	}
	for(int i=0; i<nInt; i++)
	{
		if(elemsAreLists)
		{
			if(!(meshNeedsToBeCopied && (i > 0)))
			{
				oM = PyList_GetItem(oMesh, (Py_ssize_t)i);
				if(oM == 0) throw strEr_FailedUpdateInt;
				
				if(oMeshTrue == 0) oMeshTrue = oM;
				else
				{
					if(oM == oMeshTrue)
					{//Request to copy the mesh object in Py (because it appeared to be the same as the the previous one - this happens e.g. 
					 //when the list of meshes is created as: [SRWLRadMesh(_eStart=wfr.mesh.eStart, _xStart=0, _yStart=0)]*3)
						meshNeedsToBeCopied = true;
					}
				}
			}
			if(meshNeedsToBeCopied && (i > 0))
			{
				if(oFunc == 0) 
				{
					oFunc = PyObject_GetAttrString(oMeshTrue, "copy");
					if(oFunc == 0) throw strEr_FailedUpdateInt;
				}
				if(argList == 0) 
				{
					argList = Py_BuildValue("()");
					if(argList == 0) throw strEr_FailedUpdateInt;
				}

				oM = PyObject_CallObject(oFunc, argList); 
				if(oM == 0) throw strEr_FailedUpdateInt;
			}

			if(oM == 0) throw strEr_FailedUpdateInt;

			//This copies Mesh element in Py:
			//PyObject *oFunc = PyObject_GetAttrString(oMesh, "copy");
			//PyObject *argList = Py_BuildValue("()");
			//PyObject *newMesh = PyObject_CallObject(oFunc, argList); 
		}
		else
		{
			oM = oMesh;
		}

		UpdatePyRadMesh(oM, t_arIntMesh++);

		if(elemsAreLists && meshNeedsToBeCopied && (i > 0))
		{
			if(PyList_SetItem(oMesh, (Py_ssize_t)i, oM)) throw strEr_FailedUpdateInt;
		}

		char *pCurInt = arInts[i];
		if(pCurInt != 0)
		{
			map<char*, PyObject*>::iterator it = gmBufPyObjPtr.find(pCurInt);
			if(it == gmBufPyObjPtr.end()) throw strEr_FailedUpdateInt;
			oAr = it->second;
			if(oAr == 0) throw strEr_FailedUpdateInt;
		}
		if(oAr != 0)
		{
			if(elemsAreLists)
			{
				if(PyList_SetItem(oIntAr, (Py_ssize_t)i, oAr)) throw strEr_FailedUpdateInt;
			}
			else
			{
				oIntAr = oAr;
			}
		}
	}

	if(meshListNeedsToBeSet)
	//if(meshNeedsToBeCopied)
	{
		if(PyList_SetItem(oInt, (Py_ssize_t)indMesh, oMesh)) throw strEr_FailedUpdateInt;
	}

	if(nElem > indInt)
	{
		if(PyList_SetItem(oInt, (Py_ssize_t)indInt, oIntAr)) throw strEr_FailedUpdateInt;
	}
	else
	{
		if(PyList_Append(oInt, oIntAr)) throw strEr_FailedUpdateInt;
	}
}

/************************************************************************//**
 * Auxiliary function: releases Python buffers
 ***************************************************************************/
void ReleasePyBuffers(vector<Py_buffer>& vBuf)
{
	if(vBuf.empty()) return;
	int nBuf = (int)vBuf.size();
	for(int i=0; i<nBuf; i++) 
	{
		PyBuffer_Release(&(vBuf[i]));
	}
	vBuf.erase(vBuf.begin(), vBuf.end());
}

/************************************************************************//**
 * Auxiliary function: deallocates arrays of magnetic field elements (which were allocated at parsing)
 ***************************************************************************/
void DeallocMagCntArrays(SRWLMagFldC* pMagCnt)
{
	if(pMagCnt == 0) return;
	if((pMagCnt->arMagFld != 0) && (pMagCnt->arMagFldTypes != 0) && (pMagCnt->nElem > 0))
	{
		void **t_arMagFld = pMagCnt->arMagFld;
		for(int i=0; i<pMagCnt->nElem; i++)
		{
			char cType = pMagCnt->arMagFldTypes[i];
			if(cType == 'a')
			{
				SRWLMagFld3D *pFld3D = (SRWLMagFld3D*)(*t_arMagFld);
				//arBx, arBy, arBz, arX, arY, arZ were not allocated (read via buffer)
				delete pFld3D;
			}
			else if(cType == 'm') delete (SRWLMagFldM*)(*t_arMagFld);
			else if(cType == 's') delete (SRWLMagFldS*)(*t_arMagFld);
			else if(cType == 'h') delete (SRWLMagFldH*)(*t_arMagFld);
			else if(cType == 'u')
			{
				SRWLMagFldU *pUnd = (SRWLMagFldU*)(*t_arMagFld);
				if(pUnd->arHarm != 0) delete[] pUnd->arHarm;
				delete pUnd;
			}
			else if(cType == 'c')
			{
				SRWLMagFldC *pFldC = (SRWLMagFldC*)(*t_arMagFld);
				DeallocMagCntArrays(pFldC);
				delete pFldC;
			}
			//to add new magnetic field types here
			//delete *t_arMagFld;

			*t_arMagFld = 0;
			t_arMagFld++;
		}
		delete[] pMagCnt->arMagFld; pMagCnt->arMagFld = 0;
		delete[] pMagCnt->arMagFldTypes; pMagCnt->arMagFldTypes = 0;
		//arXc, arYc, arZc, arPar1, arPar2 were not allocated (read via buffer)
	}
	else if(pMagCnt->arMagFld != 0) { delete[] pMagCnt->arMagFld; pMagCnt->arMagFld = 0;}
	else if(pMagCnt->arMagFldTypes != 0) { delete[] pMagCnt->arMagFldTypes; pMagCnt->arMagFldTypes = 0;}
}

/************************************************************************//**
 * Auxiliary function: deallocates arrays of optical elements (which were allocated at parsing)
 ***************************************************************************/
void DeallocOptCntArrays(SRWLOptC* pOptCnt)
{
	if(pOptCnt == 0) return;

	for(int i=0; i<pOptCnt->nElem; i++)
	{
		if(pOptCnt->arOpt != 0)
		{
			if(pOptCnt->arOpt[i] != 0)
			{
				if(pOptCnt->arOptTypes != 0)
				{
					char *sType = pOptCnt->arOptTypes[i];
					if(sType != 0)
					{
						if(strcmp(sType, "drift") == 0) delete (SRWLOptD*)(pOptCnt->arOpt[i]);
						else if((strcmp(sType, "aperture") == 0) || (strcmp(sType, "obstacle") == 0)) delete (SRWLOptA*)(pOptCnt->arOpt[i]);
						else if(strcmp(sType, "lens") == 0) delete (SRWLOptL*)(pOptCnt->arOpt[i]);
						else if(strcmp(sType, "zp") == 0) delete (SRWLOptZP*)(pOptCnt->arOpt[i]);
						else if(strcmp(sType, "waveguide") == 0) delete (SRWLOptWG*)(pOptCnt->arOpt[i]);
						else if(strcmp(sType, "grating") == 0) delete (SRWLOptG*)(pOptCnt->arOpt[i]);
						else if(strcmp(sType, "transmission") == 0) delete (SRWLOptT*)(pOptCnt->arOpt[i]);
						//{
						//	SRWLOptT *pT = (SRWLOptT*)(pOptCnt->arOpt[i]);
						//	if(pT->arTr != 0) delete[] (pT->arTr);
						//}
						else if(strcmp(sType, "container") == 0)
						{
							DeallocOptCntArrays((SRWLOptC*)(pOptCnt->arOpt[i]));
						}
					}
				}
				//delete pOptCnt->arOpt[i];
				pOptCnt->arOpt[i] = 0;
			}
		}
		if(pOptCnt->arOptTypes != 0)
		{
			if(pOptCnt->arOptTypes[i] != 0)
			{
				delete[] pOptCnt->arOptTypes[i];
			}
		}
	}
	if(pOptCnt->arOpt != 0) { delete[] pOptCnt->arOpt; pOptCnt->arOpt = 0;}
	if(pOptCnt->arOptTypes != 0) { delete[] pOptCnt->arOptTypes; pOptCnt->arOptTypes = 0;}
	
	if(pOptCnt->arProp != 0)
	{
		for(int i=0; i<pOptCnt->nProp; i++)
		{
			if(pOptCnt->arProp[i] != 0) delete[] pOptCnt->arProp[i];
		}
		delete[] pOptCnt->arProp; pOptCnt->arProp = 0;
	}

	if (pOptCnt->arPropN != 0)
	{
		delete[] pOptCnt->arPropN; pOptCnt->arPropN = 0;
	}
}

/************************************************************************//**
 * Wavefront modification (re-allocation) function; to be called by pointer from SRWLIB
 ***************************************************************************/
int ModifySRWLWfr(int action, SRWLWfr* pWfr, char pol)
{
	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime (&start);

	if(pWfr == 0) return -1; //returning non-zero means Wfr modification did not succeed; no throwing allowed here
	//if((action < 0) || (action > 2)) return -1;
	if(action < 0) return -1; //OC151115

	map<SRWLWfr*, AuxStructPyObjectPtrs>::iterator it = gmWfrPyPtr.find(pWfr);
	//map<SRWLWfr*, AuxStructPyObjectPtrs>::const_iterator it = gmWfrPyPtr.find(pWfr);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":ModifySRWLWfr : find",&start);

	if(it == gmWfrPyPtr.end()) return -1;
	PyObject *oWfr = it->second.o_wfr;
	if(oWfr == 0) return -1;

	//PyObject *oFunc = PyObject_GetAttrString(oWfr, "allocate");
	//if(oFunc == 0) return -1;
	//if(!PyCallable_Check(oFunc)) return -1;

	PyObject *oFunc=0; //OC151115
	PyObject *argList=0;

	//We probably need to release Py buffers of Ex, Ey (?):
	//PyBuffer_Release(&(it->second.pbEx));
	//PyBuffer_Release(&(it->second.pbEy));
	//PyBuffer_Release(&(it->second.pbMomX));
	//PyBuffer_Release(&(it->second.pbMomY));

	int ExNeeded = ((pol == 0) || (pol == 'x') || (pol == 'X')) ? 1 : 0;
	int EyNeeded = ((pol == 0) || (pol == 'y') || (pol == 'Y') || (pol == 'z') || (pol == 'Z')) ? 1 : 0;

	if(action == 0)
	{//delete existing wavefront data
		oFunc = PyObject_GetAttrString(oWfr, "allocate");

	 //trying to call an equivalent of allocate(ne, nx, ny) in Py
#if PY_MAJOR_VERSION >= 3
		argList = Py_BuildValue("(i,i,i,i,i,C)", 0, 0, 0, 1, 1, pWfr->numTypeElFld); //doesn't work with Py2.7
#else
		argList = Py_BuildValue("(i,i,i,i,i,c)", 0, 0, 0, 1, 1, pWfr->numTypeElFld);
#endif
	}
	//else if(action == 1)
	//{//allocate new wavefront data (without checking/deleting any existing data)
	//}
	//else if(action == 2)
	//{//modify wavefront size (numbers of points vs photon energy, horizontal or vertical position)
	//}
	//else
	else if((action == 2) || (action == 12))
	{
		oFunc = PyObject_GetAttrString(oWfr, "allocate");

		//int ExNeeded = ((pol == 0) || (pol == 'x') || (pol == 'X'))? 1 : 0;
		//int EyNeeded = ((pol == 0) || (pol == 'y') || (pol == 'Y') || (pol == 'z') || (pol == 'Z'))? 1 : 0;

		int backupNeeded = 0; //OC141115
		if(action == 12) backupNeeded = 1;

#if PY_MAJOR_VERSION >= 3
		//argList = Py_BuildValue("(i,i,i,i,i,C)", pWfr->ne, pWfr->nx, pWfr->ny, ExNeeded, EyNeeded, pWfr->numTypeElFld); //doesn't work with Py2.7
		//argList = Py_BuildValue("(i,i,i,i,i,C)", pWfr->mesh.ne, pWfr->mesh.nx, pWfr->mesh.ny, ExNeeded, EyNeeded, pWfr->numTypeElFld); //doesn't work with Py2.7
		argList = Py_BuildValue("(i,i,i,i,i,C,i)", pWfr->mesh.ne, pWfr->mesh.nx, pWfr->mesh.ny, ExNeeded, EyNeeded, pWfr->numTypeElFld, backupNeeded); //OC141115
#else
		//argList = Py_BuildValue("(i,i,i,i,i,c)", pWfr->ne, pWfr->nx, pWfr->ny, ExNeeded, EyNeeded, pWfr->numTypeElFld);
		//argList = Py_BuildValue("(i,i,i,i,i,c)", pWfr->mesh.ne, pWfr->mesh.nx, pWfr->mesh.ny, ExNeeded, EyNeeded, pWfr->numTypeElFld);
		argList = Py_BuildValue("(i,i,i,i,i,c,i)", pWfr->mesh.ne, pWfr->mesh.nx, pWfr->mesh.ny, ExNeeded, EyNeeded, pWfr->numTypeElFld, backupNeeded); //OC141115
#endif
	}
	else if(action == 20) //OC151115
	{//delete wavefront backup data
		oFunc = PyObject_GetAttrString(oWfr, "delE");

		int typeData = 2; //backup only
		argList = Py_BuildValue("(i,i,i)", typeData, ExNeeded, EyNeeded); //OC151115
	}

	if(argList == 0) return -1;
	if(oFunc == 0) return -1; //OC151115
	if(!PyCallable_Check(oFunc)) return -1;

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : before PyObject_CallObject",&start);

	PyObject *res = PyObject_CallObject(oFunc, argList); //re-allocate in Py
	Py_DECREF(argList);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : PyObject_CallObject",&start);

	Py_DECREF(oFunc);
	if(res == 0) return -1;
	Py_DECREF(res);

	PyObject *o_tmp = 0;
	//Py_buffer pb_tmp;
	pWfr->arEx = pWfr->arEy = 0;

	o_tmp = PyObject_GetAttrString(oWfr, "arEx");
	if(o_tmp == 0) return -1;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) return -1;
	//	it->second.pv_buf->push_back(pb_tmp);
	//	pWfr->arEx = (char*)pb_tmp.buf;
	//	it->second.pbEx = pb_tmp;
	//}
	int sizeVectBuf = 0;
	if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
	if(!(pWfr->arEx = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))) return -1;
	if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbEx = (*it->second.pv_buf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arEx",&start);

	o_tmp = PyObject_GetAttrString(oWfr, "arEy");
	if(o_tmp == 0) return -1;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) return -1;
	//	it->second.pv_buf->push_back(pb_tmp);
	//	pWfr->arEy = (char*)pb_tmp.buf;
	//	it->second.pbEy = pb_tmp;
	//}
	sizeVectBuf = 0;
	if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
	if(!(pWfr->arEy = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))) return -1;
	if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbEy = (*it->second.pv_buf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arEy",&start);

	pWfr->arExAux = 0; //OC151115
	if(PyObject_HasAttrString(oWfr, "arExAux"))
	{
		o_tmp = PyObject_GetAttrString(oWfr, "arExAux");
		if(o_tmp != 0)
		{
			sizeVectBuf = 0;
			if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
			//if(pWfr->arExAux = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))
			pWfr->arExAux = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0);
			if(pWfr->arExAux)
			{
				if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbExAux = (*it->second.pv_buf)[sizeVectBuf];
				Py_DECREF(o_tmp);
			}
		}
	}

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arExAux",&start);

	//if(pWfr->arExAux == 0)
	//{//This creates problems
	//	if((action == 20) && ExNeeded) PyBuffer_Release(&(it->second.pbExAux));
	//}

	pWfr->arEyAux = 0; //OC151115
	if(PyObject_HasAttrString(oWfr, "arEyAux"))
	{
		o_tmp = PyObject_GetAttrString(oWfr, "arEyAux");
		if(o_tmp != 0)
		{
			sizeVectBuf = 0;
			if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
			//if(pWfr->arEyAux = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))
			pWfr->arEyAux = GetPyArrayBuf(o_tmp, it->second.pv_buf, 0);
			if(pWfr->arEyAux)
			{
				if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbEyAux = (*it->second.pv_buf)[sizeVectBuf];
				Py_DECREF(o_tmp);
			}
		}
	}

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arEyAux",&start);

	//if(pWfr->arEyAux == 0)
	//{//This creates problems
	//	if((action == 20) && EyNeeded) PyBuffer_Release(&(it->second.pbEyAux));
	//}

	o_tmp = PyObject_GetAttrString(oWfr, "arMomX");
	if(o_tmp == 0) return -1;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) return -1;
	//	it->second.pv_buf->push_back(pb_tmp);
	//	pWfr->arMomX = (double*)pb_tmp.buf;
	//	it->second.pbMomX = pb_tmp;
	//}
	sizeVectBuf = 0;
	if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
	if(!(pWfr->arMomX = (double*)GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))) return -1;
	if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbMomX = (*it->second.pv_buf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arMomX",&start);

	o_tmp = PyObject_GetAttrString(oWfr, "arMomY");
	if(o_tmp == 0) return -1;
	//if(PyObject_CheckBuffer(o_tmp))
	//{
	//	if(PyObject_GetBuffer(o_tmp, &pb_tmp, PyBUF_WRITABLE)) return -1;
	//	it->second.pv_buf->push_back(pb_tmp);
	//	pWfr->arMomY = (double*)pb_tmp.buf;
	//	it->second.pbMomY = pb_tmp;
	//}
	sizeVectBuf = 0;
	if(it->second.pv_buf != 0) sizeVectBuf = (int)it->second.pv_buf->size();
	if(!(pWfr->arMomY = (double*)GetPyArrayBuf(o_tmp, it->second.pv_buf, 0))) return -1;
	if((int)it->second.pv_buf->size() > sizeVectBuf) it->second.pbMomY = (*it->second.pv_buf)[sizeVectBuf];
	Py_DECREF(o_tmp);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime("::ModifySRWLWfr : arMomY",&start);

	return 0;
}

//void* (*pExtFunc)(char type, long long len)
/************************************************************************//**
 * Array allocation function; to be called by pointer from SRWLIB
 ***************************************************************************/
char* AllocPyArrayGetBuf(char type, long long len)
{
	if(!((type == 'd') || (type == 'f') || (type == 'i'))) return 0; //returning 0 means allocation did not succeed; no throwing allowed here
	if(len <= 0) return 0;

	PyObject *oSRWLIB = PyImport_AddModule("srwlib");
	PyObject *oFunc = PyObject_GetAttrString(oSRWLIB, "srwl_uti_array_alloc");
	if((oFunc == 0) || (!PyCallable_Check(oFunc))) throw strEr_FailedAllocPyArray;

#if PY_MAJOR_VERSION >= 3
	PyObject *oArgList = Py_BuildValue("(C,l)", type, len);
#else
	PyObject *oArgList = Py_BuildValue("(c,l)", type, len);
#endif

	PyObject *oAr = PyObject_CallObject(oFunc, oArgList); //allocate array in Py
	Py_DECREF(oArgList);

	if(oAr == 0) throw strEr_FailedAllocPyArray;

	Py_ssize_t sizeBuf=0;
	char *resBuf = GetPyArrayBuf(oAr, 0, &sizeBuf);
	if((resBuf == 0) || (sizeBuf <= 0)) throw strEr_FailedAllocPyArray;

	gmBufPyObjPtr[resBuf] = oAr; //to be able to manipulate with Py objects corresponding to arrays
	return resBuf;
}

/************************************************************************//**
 * Tabulates 3D magnetic field (in Cartesian laboratory frame);
 * see help to srwlCalcMagnField
 ***************************************************************************/
static PyObject* srwlpy_CalcMagnField(PyObject *self, PyObject *args)
{
	PyObject *oDispMagCnt=0, *oMagFldCnt=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	//SRWLMagFldC magCnt = {0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	//SRWLMagFldC dispMagCnt = {0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	SRWLMagFldC magCnt = {0,0,0,0,0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	SRWLMagFldC dispMagCnt = {0,0,0,0,0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)

	try
	{
		try
		{//OC190115 //For backwards compatibility
			if(!PyArg_ParseTuple(args, "OOO:CalcMagnField", &oDispMagCnt, &oMagFldCnt, &oPrecPar)) throw strEr_BadArg_CalcMagnField;
		}
		catch(...)
		{//OC190115 //For backwards compatibility
			if(!PyArg_ParseTuple(args, "OO:CalcMagnField", &oDispMagCnt, &oMagFldCnt)) throw strEr_BadArg_CalcMagnField;
		}

		if((oDispMagCnt == 0) || (oMagFldCnt == 0)) throw strEr_BadArg_CalcMagnField;

		ParseSructSRWLMagFldC(&dispMagCnt, oDispMagCnt, &vBuf);
		if((dispMagCnt.nElem != 1) || (dispMagCnt.arMagFldTypes[0] != 'a')) throw strEr_BadArg_CalcMagnField;
		ParseSructSRWLMagFldC(&magCnt, oMagFldCnt, &vBuf);

		double arPrecPar[] = {0,0,0,0,0,0}; //to increase if necessary
		double *pPrecPar = arPrecPar;
		int nPrecPar = 6;
		if(oPrecPar != 0)
		{
			//pPrecPar = arPrecPar;
			CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);
		}

		ProcRes(srwlCalcMagFld(&dispMagCnt, &magCnt, pPrecPar));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		oDispMagCnt = 0;
	}

	DeallocMagCntArrays(&dispMagCnt);
	DeallocMagCntArrays(&magCnt);
	ReleasePyBuffers(vBuf);

	if(oDispMagCnt) Py_XINCREF(oDispMagCnt);
	return oDispMagCnt;
}

/************************************************************************//**
 * Calculates charged particle trajectory in external 3D magnetic field (in Cartesian laboratory frame);
 * see help to srwlCalcPartTraj
 ***************************************************************************/
static PyObject* srwlpy_CalcPartTraj(PyObject *self, PyObject *args)
{
	PyObject *oPartTraj=0, *oMagFldCnt=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;

	//SRWLMagFldC magCnt = {0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	SRWLMagFldC magCnt = {0,0,0,0,0,0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	//SRWLPrtTrj trj = {0,0,0,0,0,0}; //zero pointers
	SRWLPrtTrj trj = {0,0,0,0,0,0,0,0,0}; //zero pointers
	try
	{
		if(!PyArg_ParseTuple(args, "OOO:CalcPartTraj", &oPartTraj, &oMagFldCnt, &oPrecPar)) throw strEr_BadArg_CalcPartTraj;
		if((oPartTraj == 0) || (oMagFldCnt == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcPartTraj;

		ParseSructSRWLPrtTrj(&trj, oPartTraj, &vBuf);
		ParseSructSRWLMagFldC(&magCnt, oMagFldCnt, &vBuf);

		double arPrecPar[9]; //to increase if necessary
		int nPrecPar = 1;
		arPrecPar[1] = 1; //default integration method
		double *pPrecPar = arPrecPar + 1;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);
		arPrecPar[0] = nPrecPar; //!

		ProcRes(srwlCalcPartTraj(&trj, &magCnt, arPrecPar));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oPartTraj = 0;
	}

	DeallocMagCntArrays(&magCnt);
	ReleasePyBuffers(vBuf);

	if(oPartTraj) Py_XINCREF(oPartTraj);
	return oPartTraj;
}

/************************************************************************//**
 * Calculates charged particle trajectory from an array of kick matrices;
 * see help to srwlCalcPartTrajFromKickMatr
 ***************************************************************************/
static PyObject* srwlpy_CalcPartTrajFromKickMatr(PyObject *self, PyObject *args)
{
	PyObject *oPartTraj=0, *oListKickMatr=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;

	//SRWLPrtTrj trj = {0,0,0,0,0,0}; //zero pointers, since SRWL structures are definied in C (no constructors)
	SRWLPrtTrj trj = {0,0,0,0,0,0,0,0,0}; //zero pointers, since SRWL structures are definied in C (no constructors)
	SRWLKickM *arKickM = 0;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:CalcPartTrajFromKickMatr", &oPartTraj, &oListKickMatr, &oPrecPar)) throw strEr_BadArg_CalcPartTrajFromKickMatr;
		if((oPartTraj == 0) || (oListKickMatr == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcPartTrajFromKickMatr;
		
		ParseSructSRWLPrtTrj(&trj, oPartTraj, &vBuf);

		int nKickM = 0;
		if(PyList_Check(oListKickMatr))
		{
			nKickM = (int)PyList_Size(oListKickMatr);
			if(nKickM <= 0) throw strEr_BadArg_CalcPartTrajFromKickMatr;

			arKickM = new SRWLKickM[nKickM];

			for(int i=0; i<nKickM; i++)
			{
				PyObject *oKickMatr = PyList_GetItem(oListKickMatr, (Py_ssize_t)i);
				if(oKickMatr == 0) throw strEr_BadArg_CalcPartTrajFromKickMatr;

				ParseSructSRWLKickM(arKickM + i, oKickMatr, &vBuf);
			}
		}
		else 
		{
			nKickM = 1;
			arKickM = new SRWLKickM[nKickM];
			ParseSructSRWLKickM(arKickM, oListKickMatr, &vBuf);
		}

		double arPrecPar[9]; //to increase if necessary
		int nPrecPar = 1;
		double *pPrecPar = arPrecPar;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcPartTrajFromKickMatr(&trj, arKickM, nKickM, arPrecPar));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		oPartTraj = 0;
	}

	if(arKickM != 0) delete[] arKickM;

	//DeallocMagCntArrays(&magCnt);
	ReleasePyBuffers(vBuf);

	if(oPartTraj) Py_XINCREF(oPartTraj);
	return oPartTraj;
}

/************************************************************************//**
 * Calculates Wavefront (electric field) of Synchrotron Radiation 
 * by a relativistic charged particle traveling in external (3D) magnetic field;
 * see help to srwlCalcElecFieldSR
 ***************************************************************************/
static PyObject* srwlpy_CalcElecFieldSR(PyObject *self, PyObject *args)
{
	PyObject *oWfr=0, *oPartTraj=0, *oMagFldCnt=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	//SRWLMagFldC magCnt = {0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC magCnt = {0,0,0,0,0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC *pMagCnt = &magCnt;
	//SRWLPrtTrj trj = {0,0,0,0,0,0};
	SRWLPrtTrj trj = {0,0,0,0,0,0,0,0,0};
	SRWLPrtTrj *pTrj = &trj;
	SRWLWfr wfr;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOO:CalcElecFieldSR", &oWfr, &oPartTraj, &oMagFldCnt, &oPrecPar)) throw strEr_BadArg_CalcElecFieldSR;
		if((oWfr == 0) || (oPartTraj == 0) || (oMagFldCnt == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcElecFieldSR;

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);

		//if(strcmp(oPartTraj->ob_type->tp_name, "SRWLPrtTrj") == 0) ParseSructSRWLPrtTrj(pTrj, oPartTraj, &vBuf);
		//else pTrj = 0;
		char sTypeName[1025];
		CopyPyClassNameToC(oPartTraj, sTypeName, 1024);
		if(strcmp(sTypeName, "SRWLPrtTrj") == 0) ParseSructSRWLPrtTrj(pTrj, oPartTraj, &vBuf);
		else pTrj = 0;

		//if(strcmp(oMagFldCnt->ob_type->tp_name, "SRWLMagFldC") == 0) ParseSructSRWLMagFldC(pMagCnt, oMagFldCnt, &vBuf);
		//else pMagCnt = 0;
		CopyPyClassNameToC(oMagFldCnt, sTypeName, 1024);
		if(strcmp(sTypeName, "SRWLMagFldC") == 0) ParseSructSRWLMagFldC(pMagCnt, oMagFldCnt, &vBuf);
		else pMagCnt = 0;

		if((pTrj == 0) && (pMagCnt == 0)) throw strEr_BadArg_CalcElecFieldSR;

		//double *arPrecPar = (double*)GetPyArrayBuf(vBuf, oPrecPar, PyBUF_SIMPLE);
		//if(arPrecPar == 0) throw strEr_BadPrec_CalcElecFieldSR;
		double arPrecPar[7];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 7;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcElecFieldSR(&wfr, pTrj, pMagCnt, arPrecPar, nPrecPar));
		UpdatePyWfr(oWfr, &wfr);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oWfr = 0;
	}

	if(pMagCnt != 0) DeallocMagCntArrays(pMagCnt);
	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	if(oWfr) Py_XINCREF(oWfr);
	return oWfr;
}

/************************************************************************//**
 * Calculates Wavefront (electric field) of a Gaussian Beam;
 * see help to srwlCalcElecFieldSR
 ***************************************************************************/
static PyObject* srwlpy_CalcElecFieldGaussian(PyObject *self, PyObject *args)
{
	PyObject *oWfr=0, *oGsnBm=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	SRWLWfr wfr;
	SRWLGsnBm gsnBm;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:CalcElecFieldGaussian", &oWfr, &oGsnBm, &oPrecPar)) throw strEr_BadArg_CalcElecFieldGaussian;
		if((oWfr == 0) || (oGsnBm == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcElecFieldGaussian;

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);
		ParseSructSRWLGsnBm(&gsnBm, oGsnBm);

		double arPrecPar[1];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 1;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcElecFieldGaussian(&wfr, &gsnBm, arPrecPar));
		UpdatePyWfr(oWfr, &wfr);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		oWfr = 0;
	}

	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	if(oWfr) Py_XINCREF(oWfr);
	return oWfr;
}

/************************************************************************//**
 * Calculates Wavefront (electric field) of a Gaussian Beam;
 * see help to srwlCalcElecFieldSR
 ***************************************************************************/
static PyObject* srwlpy_CalcElecFieldPointSrc(PyObject *self, PyObject *args)
{
	PyObject *oWfr=0, *oPtSrc=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	SRWLWfr wfr;
	SRWLPtSrc ptSrc;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:CalcElecFieldSpherWave", &oWfr, &oPtSrc, &oPrecPar)) throw strEr_BadArg_CalcElecFieldSpherWave;
		if((oWfr == 0) || (oPtSrc == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcElecFieldSpherWave;

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);
		ParseSructSRWLPtSrc(&ptSrc, oPtSrc);

		double arPrecPar[1];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 1;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcElecFieldPointSrc(&wfr, &ptSrc, arPrecPar));
		UpdatePyWfr(oWfr, &wfr);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		oWfr = 0;
	}

	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	if(oWfr) Py_XINCREF(oWfr);
	return oWfr;
}

/************************************************************************//**
 * Calculates Stokes parameters of Undulator Radiation by a relativistic 
 * finite-emittance electron beam traveling in periodic magnetic field of an undulator;
 * see help to srwlCalcStokesUR
 ***************************************************************************/
static PyObject* srwlpy_CalcStokesUR(PyObject *self, PyObject *args)
{
	PyObject *oStokes=0, *oElBeam=0, *oUnd=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	SRWLStokes stokes;
	SRWLPartBeam eBeam;
	SRWLMagFldU und;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOO:CalcStokesUR", &oStokes, &oElBeam, &oUnd, &oPrecPar)) throw strEr_BadArg_CalcStokesUR;
		if((oStokes == 0) || (oElBeam == 0) || (oUnd == 0) || (oPrecPar == 0)) throw strEr_BadArg_CalcStokesUR;

		ParseSructSRWLStokes(&stokes, oStokes, &vBuf);
		ParseSructSRWLPartBeam(&eBeam, oElBeam, vBuf);
		ParseSructSRWLMagFldU(&und, oUnd);

		double arPrecPar[5];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 5;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcStokesUR(&stokes, &eBeam, &und, arPrecPar));
		UpdatePyStokes(oStokes, &stokes);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oStokes = 0;
	}

	if(und.arHarm != 0) delete[] und.arHarm;
	ReleasePyBuffers(vBuf);

	if(oStokes) Py_XINCREF(oStokes);
	return oStokes;
}

/************************************************************************//**
 * Calculates Power Density distribution of Synchrotron Radiation 
 * by a relativistic finite-emittance electron beam traveling in arbitrary magnetic field;
 * see help to srwlCalcStokesUR
 ***************************************************************************/
static PyObject* srwlpy_CalcPowDenSR(PyObject *self, PyObject *args)
{
	PyObject *oStokes=0, *oElBeam=0, *oPartTraj=0, *oMagFldCnt=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;
	SRWLStokes stokes;
	SRWLPartBeam eBeam;
	//SRWLMagFldC magCnt = {0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC magCnt = {0,0,0,0,0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC *pMagCnt = &magCnt;
	//SRWLPrtTrj trj = {0,0,0,0,0,0};
	SRWLPrtTrj trj = {0,0,0,0,0,0,0,0,0};
	SRWLPrtTrj *pTrj = &trj;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOOO:CalcPowDenSR", &oStokes, &oElBeam, &oPartTraj, &oMagFldCnt, &oPrecPar)) throw strEr_BadArg_CalcPowDenSR;
		if((oStokes == 0) || (oElBeam == 0) || ((oPartTraj == 0) && (oMagFldCnt == 0)) || (oPrecPar == 0)) throw strEr_BadArg_CalcPowDenSR;

		ParseSructSRWLStokes(&stokes, oStokes, &vBuf);
		ParseSructSRWLPartBeam(&eBeam, oElBeam, vBuf);

		//if(strcmp(oPartTraj->ob_type->tp_name, "SRWLPrtTrj") == 0) ParseSructSRWLPrtTrj(pTrj, oPartTraj, &vBuf);
		//else pTrj = 0;
		char sTypeName[1025];
		CopyPyClassNameToC(oPartTraj, sTypeName, 1024);
		if(strcmp(sTypeName, "SRWLPrtTrj") == 0) ParseSructSRWLPrtTrj(pTrj, oPartTraj, &vBuf);
		else pTrj = 0;

		//if(strcmp(oMagFldCnt->ob_type->tp_name, "SRWLMagFldC") == 0) ParseSructSRWLMagFldC(pMagCnt, oMagFldCnt, &vBuf);
		//else pMagCnt = 0;
		CopyPyClassNameToC(oMagFldCnt, sTypeName, 1024);
		if(strcmp(sTypeName, "SRWLMagFldC") == 0) ParseSructSRWLMagFldC(pMagCnt, oMagFldCnt, &vBuf);
		else pMagCnt = 0;

		if((pTrj == 0) && (pMagCnt == 0)) throw strEr_BadArg_CalcPowDenSR;

		double arPrecPar[5];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 5;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlCalcPowDenSR(&stokes, &eBeam, pTrj, pMagCnt, arPrecPar));
		UpdatePyStokes(oStokes, &stokes);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oStokes = 0;
	}

	if(pMagCnt != 0) DeallocMagCntArrays(pMagCnt);
	ReleasePyBuffers(vBuf);

	if(oStokes) Py_XINCREF(oStokes);
	return oStokes;
}

/************************************************************************//**
 * Calculates/extracts Intensity and/or other characteristics from pre-calculated Electric Field
 * see help to srwlCalcIntFromElecField
 ***************************************************************************/
static PyObject* srwlpy_CalcIntFromElecField(PyObject *self, PyObject *args)
{
	PyObject *oInt=0, *oWfr=0, *oPol=0, *oIntType=0, *oDepType=0, *oE=0, *oX=0, *oY=0;
	vector<Py_buffer> vBuf;
	SRWLWfr wfr;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOOOOOO:CalcIntFromElecField", &oInt, &oWfr, &oPol, &oIntType, &oDepType, &oE, &oX, &oY)) throw strEr_BadArg_CalcIntFromElecField;
		if((oInt == 0) || (oWfr == 0) || (oPol == 0) || (oIntType == 0) || (oDepType == 0) || (oE == 0) || (oX == 0) || (oY == 0)) throw strEr_BadArg_CalcIntFromElecField;

		//char *arInt = (char*)GetPyArrayBuf(oInt, vBuf, PyBUF_WRITABLE, 0);
		char *arInt = (char*)GetPyArrayBuf(oInt, &vBuf, 0);

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);

		if(!PyNumber_Check(oPol)) throw strEr_BadArg_CalcIntFromElecField;
		char pol = (char)PyLong_AsLong(oPol);

		if(!PyNumber_Check(oIntType)) throw strEr_BadArg_CalcIntFromElecField;
		char intType = (char)PyLong_AsLong(oIntType);

		if(!PyNumber_Check(oDepType)) throw strEr_BadArg_CalcIntFromElecField;
		char depType = (char)PyLong_AsLong(oDepType);

		if(!PyNumber_Check(oE)) throw strEr_BadArg_CalcIntFromElecField;
		double e = PyFloat_AsDouble(oE);

		if(!PyNumber_Check(oX)) throw strEr_BadArg_CalcIntFromElecField;
		double x = PyFloat_AsDouble(oX);

		if(!PyNumber_Check(oY)) throw strEr_BadArg_CalcIntFromElecField;
		double y = PyFloat_AsDouble(oY);

		ProcRes(srwlCalcIntFromElecField(arInt, &wfr, pol, intType, depType, e, x, y));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oInt = 0;
	}

	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	if(oInt) Py_XINCREF(oInt);
	return oInt;
}

/************************************************************************//**
 * "Resizes" Electric Field Wavefront vs transverse positions / angles or photon energy / time
 * see help to srwlResizeElecField
 ***************************************************************************/
static PyObject* srwlpy_ResizeElecField(PyObject *self, PyObject *args)
{
	PyObject *oWfr=0, *oType, *oPar;
	vector<Py_buffer> vBuf;
	SRWLWfr wfr;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:ResizeElecField", &oWfr, &oType, &oPar)) throw strEr_BadArg_ResizeElecField;
		if((oWfr == 0) || (oType == 0) || (oPar == 0)) throw strEr_BadArg_ResizeElecField;

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);

		//PyObject *o_str = PyUnicode_AsUTF8String(oType);
		//if(!PyBytes_Check(o_str)) throw strEr_BadArg_ResizeElecField;
		//char typeRes = *PyBytes_AsString(o_str); 
		//Py_INCREF(o_str);
		char cTypeRes[2];
		CopyPyStringToC(oType, cTypeRes, 1);

		//double arPar[] = {0.,1.,1.,1.,1.}; int nPar = 5; double *pPar = arPar;
		double arPar[] = {0.,1.,1.,1.,1.,0.5,0.5}; int nPar = 7; double *pPar = arPar; //OC071014
		CopyPyListElemsToNumArray(oPar, 'd', pPar, nPar);

		if((nPar < 4) && ((cTypeRes[0] == 'f') || (cTypeRes[0] == 't') || (cTypeRes[0] == 'F') || (cTypeRes[0] == 'T'))) 
		{//OC081014
			arPar[3] = 0.5; arPar[4] = 0.5; 
		}

		ProcRes(srwlResizeElecField(&wfr, *cTypeRes, arPar));
		UpdatePyWfr(oWfr, &wfr);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oWfr = 0;
	}

	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	if(oWfr) Py_XINCREF(oWfr);
	return oWfr;
}

/************************************************************************//**
 * Changes Representation of Electric Field: coordinates<->angles, frequency<->time
 * see help to srwlSetRepresElecField
 ***************************************************************************/
static PyObject* srwlpy_SetRepresElecField(PyObject *self, PyObject *args)
{
	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//double start;
	//get_walltime (&start);

	PyObject *oWfr=0, *oRepr;
	vector<Py_buffer> vBuf;
	SRWLWfr wfr;

	try
	{
		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_SetRepresElecField : begin",&start);

		if(!PyArg_ParseTuple(args, "OO:SetRepresElecField", &oWfr, &oRepr)) throw strEr_BadArg_SetRepresElecField;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_SetRepresElecField : PyArg_ParseTuple",&start);

		if((oWfr == 0) || (oRepr == 0)) throw strEr_BadArg_SetRepresElecField;

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_SetRepresElecField : ParseSructSRWLWfr",&start);

		//PyObject *o_str = PyUnicode_AsUTF8String(oRepr);
		//if(!PyBytes_Check(o_str)) throw strEr_BadArg_SetRepresElecField;
		//char repr = *PyBytes_AsString(o_str); 
		//Py_INCREF(o_str);
		char cRepr[2];
		CopyPyStringToC(oRepr, cRepr, 1);

		ProcRes(srwlSetRepresElecField(&wfr, *cRepr));

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_SetRepresElecField : srwlSetRepresElecField",&start);

		UpdatePyWfr(oWfr, &wfr);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_SetRepresElecField : UpdatePyWfr",&start);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oWfr = 0;
	}

	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":srwlpy_SetRepresElecField : EraseElementFromMap",&start);

	if(oWfr) Py_XINCREF(oWfr);

	//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
	//srwlPrintTime(":srwlpy_SetRepresElecField : Py_XINCREF",&start);

	return oWfr;
}

/************************************************************************//**
 * "Propagates" Electric Field Wavefront through Optical Elements and free space
 * see help to srwlPropagElecField
 ***************************************************************************/
static PyObject* srwlpy_PropagElecField(PyObject *self, PyObject *args)
{
	//PyObject *oWfr=0, *oOptCnt=0;
	PyObject *oWfr=0, *oOptCnt=0, *oInt=0; //OC14082018

	vector<Py_buffer> vBuf;
	SRWLWfr wfr;
	SRWLOptC optCnt = {0,0,0,0,0}; //since SRWL structures are definied in C (no constructors)
	
	//char *arIndsInt=0, *arIntType=0, *arPol=0, *arDepType=0; //OC14082018
	//double *arE=0, *arX=0, *arY=0;
	char *arIntDescr[] = {0,0,0,0,0}; //OC14082018
	SRWLRadMesh *arIntMesh=0;
	//float **arInts=0;
	char **arInts=0;

	try
	{
		//if(!PyArg_ParseTuple(args, "OO:PropagElecField", &oWfr, &oOptCnt)) throw strEr_BadArg_PropagElecField;
		if(!PyArg_ParseTuple(args, "OO|O:PropagElecField", &oWfr, &oOptCnt, &oInt)) throw strEr_BadArg_PropagElecField; //OC14082018
		if((oWfr == 0) || (oOptCnt == 0)) throw strEr_BadArg_PropagElecField;

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//double start;
		//get_walltime(&start);

		ParseSructSRWLWfr(&wfr, oWfr, &vBuf, gmWfrPyPtr);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_PropagElecField : ParseSructSRWLWfr", &start);

		ParseSructSRWLOptC(&optCnt, oOptCnt, &vBuf);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_PropagElecField :ParseSructSRWLOptC", &start);

		int nInt = 0;
		if(oInt != 0) //OC14082018
		{
			nInt = ParseSructSRWLPropIntDef(arIntDescr, arIntMesh, oInt);
			if(nInt > 0)
			{
				//arInts = new float*[nInt];
				arInts = new char*[nInt];
				for(int i=0; i<nInt; i++) arInts[i] = 0;
			}
		}
		//OCTEST
		//char *pRes = AllocPyArrayGetBuf('d', 1000);
		//END OCTEST

		//ProcRes(srwlPropagElecField(&wfr, &optCnt));
		ProcRes(srwlPropagElecField(&wfr, &optCnt, nInt, arIntDescr, arIntMesh, arInts)); //OC15082018

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_PropagElecField :srwlPropagElecField", &start);

		if((oInt != 0) && (nInt > 0)) //OC14082018
		{//Find and add objects corresponding to different intensity distributions to the oInt list
			UpdatePyPropInt(oInt, arIntMesh, arInts, nInt);
		}

		UpdatePyWfr(oWfr, &wfr);

		//Added by S.Yakubov (for profiling?) at parallelizing SRW via OpenMP:
		//srwlPrintTime(":srwlpy_PropagElecField :UpdatePyWfr", &start);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oWfr = 0;
	}

	DeallocOptCntArrays(&optCnt);
	ReleasePyBuffers(vBuf);
	EraseElementFromMap(&wfr, gmWfrPyPtr);

	for(int i=0; i<4; i++) if(arIntDescr[i] != 0) delete[] arIntDescr[i];
	if(arIntMesh != 0) delete[] arIntMesh;
	if(arInts != 0) delete[] arInts;
	//arInts[i] should not be deleted, because these should be available in Py

	if(oWfr) Py_XINCREF(oWfr);
	return oWfr;
}

/************************************************************************//**
 * Performs FFT (1D or 2D, depending on dimensionality of input arrays)
 ***************************************************************************/
static PyObject* srwlpy_UtiFFT(PyObject *self, PyObject *args)
{
	PyObject *oData=0, *oMesh=0, *oDir=0;
	vector<Py_buffer> vBuf;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:UtiFFT", &oData, &oMesh, &oDir)) throw strEr_BadArg_UtiFFT;
		if((oData == 0) || (oMesh == 0) || (oDir == 0)) throw strEr_BadArg_UtiFFT;

		//int sizeVectBuf = (int)vBuf.size();
		char *pcData=0;
		Py_ssize_t sizeBuf;
		if(!(pcData = GetPyArrayBuf(oData, &vBuf, &sizeBuf))) throw strEr_BadArg_UtiFFT;
		//if((int)vBuf.size() > sizeVectBuf) sPyObjectPtrs.pbEx = (*pvBuf)[sizeVectBuf];
		//Py_DECREF(o_tmp);

		double arMesh[6];
		double *pMesh = arMesh;
		//int nMesh=0;
		int nMesh=6; //OC03092016 (should match arMesh[6] !)
		//CopyPyListElemsToNumArray(oMesh, 'd', pMesh, nMesh);
		char meshArType = CopyPyListElemsToNumArray(oMesh, 'd', pMesh, nMesh); //OC03092016
		if(nMesh < 3) throw strEr_BadArg_UtiFFT;

		char typeData = 'f';
		long nPt = (long)arMesh[2];
		if(nMesh >= 6) nPt *= (long)arMesh[5];
		long nElemTest = (long)(sizeBuf/sizeof(float));
		if(nElemTest != 2*nPt)
		{
			nElemTest = (long)(sizeBuf/sizeof(double));
			if(nElemTest != 2*nPt) throw strEr_BadArg_UtiFFT;
			else 
			{
				typeData = 'd'; //Yet to implement
				throw strEr_FloatArrayRequired;
			}
		}

		if(!PyNumber_Check(oDir)) throw strEr_BadArg_UtiFFT;
		int dir = (int)PyLong_AsLong(oDir);

		ProcRes(srwlUtiFFT(pcData, typeData, arMesh, nMesh, dir));

		if(meshArType == 'l') UpdatePyListNum(oMesh, arMesh, nMesh); //04092016
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//if(vBuf.size() > 0) ReleasePyBuffers(vBuf);
		oData = 0; oMesh = 0; oDir = 0;
	}

	ReleasePyBuffers(vBuf);

	if(oData) Py_XINCREF(oData);
	return oData;
}

/************************************************************************//**
 * Performs FFT (1D or 2D, depending on dimensionality of input arrays)
 ***************************************************************************/
static PyObject* srwlpy_UtiConvWithGaussian(PyObject *self, PyObject *args)
{
	PyObject *oData=0, *oMesh=0, *oSig=0;
	vector<Py_buffer> vBuf;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:UtiConvWithGaussian", &oData, &oMesh, &oSig)) throw strEr_BadArg_UtiConvWithGaussian;
		if((oData == 0) || (oMesh == 0) || (oSig == 0)) throw strEr_BadArg_UtiConvWithGaussian;

		//int sizeVectBuf = (int)vBuf.size();
		char *pcData=0;
		Py_ssize_t sizeBuf;
		if(!(pcData = GetPyArrayBuf(oData, &vBuf, &sizeBuf))) throw strEr_BadArg_UtiConvWithGaussian;
		//if((int)vBuf.size() > sizeVectBuf) sPyObjectPtrs.pbEx = (*pvBuf)[sizeVectBuf];
		//Py_DECREF(o_tmp);

		double arMesh[8];
		double *pMesh = arMesh;
		int nMesh=8;
		CopyPyListElemsToNumArray(oMesh, 'd', pMesh, nMesh);
		if(nMesh < 3) throw strEr_BadArg_UtiConvWithGaussian;

		char typeData = 'f';
		long nPt = (long)arMesh[2];
		int nDim = 1;
		if(nMesh >= 6) 
		{
			long ny = (long)arMesh[5];
			if(ny > 1)
			{
				nPt *= ny;
				nDim = 2;
			}
		}
		long nElemTest = (long)(sizeBuf/sizeof(float));
		if(nElemTest != nPt)
		{
			nElemTest = (long)(sizeBuf/sizeof(double));
			if(nElemTest != nPt) throw strEr_BadArg_UtiConvWithGaussian;
			else 
			{
				typeData = 'd';
				throw strEr_FloatArrayRequired;
			}
		}

		double arSig[3]; //[2];
		arSig[2] = 0.; //cross-term
		double *pSig = arSig;
		int nSig=3; //2;
		CopyPyListElemsToNumArray(oSig, 'd', pSig, nSig);
		if(nSig < nDim) throw strEr_BadArg_UtiConvWithGaussian;

		ProcRes(srwlUtiConvWithGaussian(pcData, typeData, arMesh, nMesh, arSig));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//if(vBuf.size() > 0) ReleasePyBuffers(vBuf);
		//oData = 0; oMesh = 0; oSig = 0;
	}

	ReleasePyBuffers(vBuf);

	if(oData) Py_XINCREF(oData);
	return oData;
}

/************************************************************************//**
 * Calculates basic statistical characteristics of intensity distribution
 ***************************************************************************/
static PyObject* srwlpy_UtiIntInf(PyObject *self, PyObject *args)
{
	PyObject *oData=0, *oMesh=0, *oRes=0;
	vector<Py_buffer> vBuf;

	try
	{
		if(!PyArg_ParseTuple(args, "OO:UtiIntInf", &oData, &oMesh)) throw strEr_BadArg_UtiIntInf;
		if((oData == 0) || (oMesh == 0)) throw strEr_BadArg_UtiIntInf;

		char *pcData=0;
		Py_ssize_t sizeBuf;
		if(!(pcData = GetPyArrayBuf(oData, &vBuf, &sizeBuf))) throw strEr_BadArg_UtiIntInf;

		SRWLRadMesh mesh;
		ParseSructSRWLRadMesh(&mesh, oMesh);

		//Py_buffer curBuf = vBuf[vBuf.size() - 1]; //Not compatible with old buffer enterf.
		//Py_ssize_t dataItemSize = curBuf.itemsize;

		Py_ssize_t dataItemSize = (Py_ssize_t)round((sizeBuf/(mesh.ne*mesh.nx*mesh.ny)));
		char typeData = 0;
		if(dataItemSize == (Py_ssize_t)sizeof(float)) typeData = 'f';
		else if(dataItemSize == (Py_ssize_t)sizeof(double)) typeData = 'd';
		else throw strEr_BadArg_UtiIntInf;

		const int nInf = 7;
		double resInf[nInf];
		ProcRes(srwlUtiIntInf(resInf, pcData, typeData, &mesh));

		oRes = SetPyListOfLists(resInf, nInf, 1, (char*)"d");
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}

	ReleasePyBuffers(vBuf);

	if(oRes) Py_XINCREF(oRes);
	return oRes;
}

/************************************************************************//**
 * Performs misc. operations on input 
 ***************************************************************************/
static PyObject* srwlpy_UtiIntProc(PyObject *self, PyObject *args)
{
	PyObject *oInt1=0, *oMesh1=0, *oInt2=0, *oMesh2=0, *oPar=0;
	vector<Py_buffer> vBuf;
	double *arPar=0;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOOO:UtiIntProc", &oInt1, &oMesh1, &oInt2, &oMesh2, &oPar)) throw strEr_BadArg_UtiIntProc;
		if((oInt1 == 0) || (oMesh1 == 0) || (oInt2 == 0) || (oMesh2 == 0) || (oPar == 0)) throw strEr_BadArg_UtiIntProc;

		SRWLRadMesh mesh1, mesh2;
		ParseSructSRWLRadMesh(&mesh1, oMesh1);
		ParseSructSRWLRadMesh(&mesh2, oMesh2);

		Py_ssize_t sizeBuf;
		char *pcInt1=0;
		if(!(pcInt1 = GetPyArrayBuf(oInt1, &vBuf, &sizeBuf))) throw strEr_BadArg_UtiIntProc;
		
		char typeInt1=0;
		Py_ssize_t intItemSize = (Py_ssize_t)round((sizeBuf/(mesh1.ne*mesh1.nx*mesh1.ny)));
		if(intItemSize == (Py_ssize_t)sizeof(float)) typeInt1 = 'f';
		else if(intItemSize == (Py_ssize_t)sizeof(double)) typeInt1 = 'd';
		else throw strEr_BadArg_UtiIntProc;
		
		char *pcInt2=0;
		if(!(pcInt2 = GetPyArrayBuf(oInt2, &vBuf, &sizeBuf))) throw strEr_BadArg_UtiIntProc;

		char typeInt2=0;
		intItemSize = (Py_ssize_t)round((sizeBuf/(mesh2.ne*mesh2.nx*mesh2.ny)));
		if(intItemSize == (Py_ssize_t)sizeof(float)) typeInt2 = 'f';
		else if(intItemSize == (Py_ssize_t)sizeof(double)) typeInt2 = 'd';
		else throw strEr_BadArg_UtiIntProc;

		int nPar=0;
		CopyPyListElemsToNumArray(oPar, 'd', arPar, nPar);
		if(nPar < 1) throw strEr_BadArg_UtiIntProc;

		ProcRes(srwlUtiIntProc(pcInt1, typeInt1, &mesh1, pcInt2, typeInt2, &mesh2, arPar));
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
	}

	ReleasePyBuffers(vBuf);
	if(arPar) delete[] arPar;

	if(oInt2) Py_XINCREF(oInt2);
	return oInt2;
}

/************************************************************************//**
 * Attempts to deduce parameters of periodic undulator magnetic field from tabulated field
 ***************************************************************************/
static PyObject* srwlpy_UtiUndFromMagFldTab(PyObject *self, PyObject *args)
{
	PyObject *oUndC=0, *oFld3DC=0, *oPrecPar=0;
	vector<Py_buffer> vBuf;

	SRWLMagFldC undCnt = {0,0,0,0,0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC *pUndCnt = &undCnt;
	SRWLMagFldC magCnt = {0,0,0,0,0,0,0,0,0,0}; //just zero pointers
	SRWLMagFldC *pMagCnt = &magCnt;

	try
	{
		if(!PyArg_ParseTuple(args, "OOO:UtiUndFromMagFldTab", &oUndC, &oFld3DC, &oPrecPar)) throw strEr_BadArg_UtiUndFromMagFldTab;
		if((oUndC == 0) || (oFld3DC == 0) || (oPrecPar == 0)) throw strEr_BadArg_UtiUndFromMagFldTab;

		ParseSructSRWLMagFldC(pUndCnt, oUndC, &vBuf);
		ParseSructSRWLMagFldC(pMagCnt, oFld3DC, &vBuf);

		double arPrecPar[3];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 3;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlUtiUndFromMagFldTab(pUndCnt, pMagCnt, arPrecPar));

			//OCTEST_161214
			//SRWLMagFldU *pMagFldU = (SRWLMagFldU*)(pUndCnt->arMagFld[0]);
			//SRWLMagFldH *pMagFldH = pMagFldU->arHarm;
			//std::stringstream ss;
			//ss << "n=" << pMagFldH->n << ", h_or_v=" << pMagFldH->h_or_v << ", B=" << pMagFldH->B;
			//string s2throw = ss.str();
			//char *sOut = new char[1000];
			//strcpy(sOut, s2throw.c_str());
			//throw sOut;

		UpdatePyMagFldC(oUndC, pUndCnt);

			//SRWLMagFldC undCntTest = {0,0,0,0,0,0,0,0,0,0}; //just zero pointers
			//SRWLMagFldC *pUndCntTest = &undCntTest;
			//ParseSructSRWLMagFldC(pUndCntTest, oUndC, &vBuf);
			//SRWLMagFldU *pTestU = (SRWLMagFldU*)(pUndCntTest->arMagFld[0]);
			//int aha = 1;
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
		oUndC = 0;
	}

	if(pUndCnt != 0) DeallocMagCntArrays(pUndCnt);
	if(pMagCnt != 0) DeallocMagCntArrays(pMagCnt);
	ReleasePyBuffers(vBuf);

	if(oUndC) Py_XINCREF(oUndC);
	return oUndC;
}

/************************************************************************//**
 * Finds indexes of undulator gap and phase values and associated magnetic fields requiired to be used in field interpolation based on gap and phase
 * see help to srwlUtiUndFindMagFldInterpInds
 ***************************************************************************/
static PyObject* srwlpy_UtiUndFindMagFldInterpInds(PyObject *self, PyObject *args)
{
	PyObject *oResInds=0, *oGaps=0, *oPhases=0, *oPrecPar=0;
	int *arResInds=0;
	double *arGaps=0, *arPhases=0;
	int nResInds=0;

	try
	{
		if(!PyArg_ParseTuple(args, "OOOO:UtiUndFindMagFldInterpInds", &oResInds, &oGaps, &oPhases, &oPrecPar)) throw strEr_BadArg_UtiUndFindMagFldInterpInds;
		if((oResInds == 0) || (oGaps == 0) || (oPhases == 0) || (oPrecPar == 0)) throw strEr_BadArg_UtiUndFindMagFldInterpInds;

		CopyPyListElemsToNumArray(oResInds, 'i', arResInds, nResInds);

		int nGaps=0, nPhases=0;
		CopyPyListElemsToNumArray(oGaps, 'd', arGaps, nGaps);
		CopyPyListElemsToNumArray(oPhases, 'd', arPhases, nPhases);

		if((arGaps != 0) && (arPhases != 0))
		{
			if(nGaps != nPhases) throw strEr_BadArg_UtiUndFindMagFldInterpInds;
		}

		double arPrecPar[5];
		double *pPrecPar = arPrecPar;
		int nPrecPar = 5;
		CopyPyListElemsToNumArray(oPrecPar, 'd', pPrecPar, nPrecPar);

		ProcRes(srwlUtiUndFindMagFldInterpInds(arResInds, &nResInds, arGaps, arPhases, nGaps, arPrecPar));

		UpdatePyListNum(oResInds, arResInds, nResInds);
		UpdatePyListNum(oPrecPar, arPrecPar, nPrecPar);
	}
	catch(const char* erText) 
	{
		PyErr_SetString(PyExc_RuntimeError, erText);
		//PyErr_PrintEx(1);
	}
	if(arResInds != 0) delete[] arResInds;
	if(arGaps != 0) delete[] arGaps;
	if(arPhases != 0) delete[] arPhases;

	PyObject *oResNumInds = Py_BuildValue("i", nResInds);
	Py_XINCREF(oResNumInds); //?
	return oResNumInds;
}

/************************************************************************//**
 * Python C API stuff: module & method definition2, etc.
 ***************************************************************************/

//This is from Python C API docs (seems to be a strict minimum):
//struct srwlpy_state {
//    PyObject *error;
//};
//
//#if PY_MAJOR_VERSION >= 3
//#define GETSTATE(m) ((struct srwlpy_state*)PyModule_GetState(m))
//#else
//#define GETSTATE(m) (&_state)
//static struct srwlpy_state _state;
//#endif
//
//static PyObject* error_out(PyObject *m) 
//{
//    struct srwlpy_state *st = GETSTATE(m);
//    PyErr_SetString(st->error, "something bad happened");
//    return NULL;
//}

static PyMethodDef srwlpy_methods[] = {
	{"CalcMagnField", srwlpy_CalcMagnField, METH_VARARGS, "CalcMagnField() Calculates (tabulates) 3D magnetic field created by multiple elements"},
	{"CalcPartTraj", srwlpy_CalcPartTraj, METH_VARARGS, "CalcPartTraj() Calculates charged particle trajectory in external 3D magnetic field (in Cartesian laboratory frame)"},
	{"CalcPartTrajFromKickMatr", srwlpy_CalcPartTrajFromKickMatr, METH_VARARGS, "CalcPartTrajFromKickMatr() Calculates charged particle trajectory from an array of kick matrices"},
	{"CalcElecFieldSR", srwlpy_CalcElecFieldSR, METH_VARARGS, "CalcElecFieldSR() Calculates Electric Field (Wavefront) of Synchrotron Radiation by a relativistic charged particle traveling in external 3D magnetic field"},
	{"CalcElecFieldGaussian", srwlpy_CalcElecFieldGaussian, METH_VARARGS, "CalcElecFieldGaussian() Calculates Electric Field (Wavefront) of a coherent Gaussian Beam"},
	{"CalcElecFieldPointSrc", srwlpy_CalcElecFieldPointSrc, METH_VARARGS, "CalcElecFieldPointSrc() Calculates Electric Field (Wavefront) of a spherical wave"},
	{"CalcStokesUR", srwlpy_CalcStokesUR, METH_VARARGS, "CalcStokesUR() Calculates Stokes parameters of Synchrotron Radiation by a relativistic finite-emittance electron beam traveling in periodic magnetic field of an undulator"},
	{"CalcPowDenSR", srwlpy_CalcPowDenSR, METH_VARARGS, "CalcPowDenSR() Calculates Power Density distribution of Synchrotron Radiation by a relativistic finite-emittance electron beam traveling in arbitrary magnetic field"},
	{"CalcIntFromElecField", srwlpy_CalcIntFromElecField, METH_VARARGS, "CalcIntFromElecField() Calculates/extracts Intensity from pre-calculated Electric Field"},
	{"ResizeElecField", srwlpy_ResizeElecField, METH_VARARGS, "ResizeElecField() \"Resizes\" Electric Field Wavefront vs transverse positions / angles or photon energy / time"},
	{"SetRepresElecField", srwlpy_SetRepresElecField, METH_VARARGS, "SetRepresElecField() Changes Representation of Electric Field: coordinates<->angles, frequency<->time"},
	{"PropagElecField", srwlpy_PropagElecField, METH_VARARGS, "PropagElecField() \"Propagates\" Electric Field Wavefront through Optical Elements and free space"},
	{"UtiFFT", srwlpy_UtiFFT, METH_VARARGS, "UtiFFT() Performs 1D or 2D FFT (as defined by arguments)"},
	{"UtiConvWithGaussian", srwlpy_UtiConvWithGaussian, METH_VARARGS, "UtiConvWithGaussian() Performs convolution of 1D or 2D data wave with 1D or 2D Gaussian (as defined by arguments)"},
	{"UtiIntInf", srwlpy_UtiIntInf, METH_VARARGS, "UtiIntInf() Calculates basic statistical characteristics of intensity distribution"},
	{"UtiIntProc", srwlpy_UtiIntProc, METH_VARARGS, "UtiIntProc() Performs misc. operations on one or two intensity distributions"},
	{"UtiUndFromMagFldTab", srwlpy_UtiUndFromMagFldTab, METH_VARARGS, "UtiUndFromMagFldTab() Attempts to create periodic undulator structure from tabulated magnetic field"},
	{"UtiUndFindMagFldInterpInds", srwlpy_UtiUndFindMagFldInterpInds, METH_VARARGS, "UtiUndFindMagFldInterpInds() Finds indexes of undulator gap and phase values and associated magnetic fields requiired to be used in field interpolation based on gap and phase"},
	{NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef srwlpymodule = {
    PyModuleDef_HEAD_INIT,
    "srwlpy",
    "srwlpy module is Python binding of Synchrotron Radiation Workshop (SRW) Library",
    -1,
    srwlpy_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_srwlpy(void)
{ 
	//setting pointer to function to be eventually called from SRWLIB
	srwlUtiSetWfrModifFunc(&ModifySRWLWfr);
	srwlUtiSetAllocArrayFunc(&AllocPyArrayGetBuf); //OC15082018

    return PyModule_Create(&srwlpymodule);
}

#else

PyMODINIT_FUNC initsrwlpy(void)
{
	//setting pointer to function to be eventually called from SRWLIB
	srwlUtiSetWfrModifFunc(&ModifySRWLWfr);
	srwlUtiSetAllocArrayFunc(&AllocPyArrayGetBuf); //OC15082018

	Py_InitModule("srwlpy", srwlpy_methods);
}

#endif
