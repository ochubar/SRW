/************************************************************************//**
 * File: pyparse.h
 * Description: Python / C parsing interface functions (supporting both Py2.* and 3.*) 
 * Project: 
 * First release: 2018
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __PYPARSE_H
#define __PYPARSE_H

//#include <ctype.h>
#include <vector>
//#include <string>
//#include <cstring>

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

//-------------------------------------------------------------------------
	
static const char strEr_BadArray[] = "Incorrect or no Python Array structure";
static const char strEr_BadList[] = "Incorrect or no Python List structure";
static const char strEr_BadListArray[] = "Incorrect or no Python List or Array structure";
static const char strEr_BadNum[] = "Incorrect or no Python number";
static const char strEr_BadStr[] = "Error at parsing / converting Python string or byte array";
static const char strEr_BadClassName[] = "Error at retrieving Python class name";
static const char strEr_FailedCreateList[] = "Failed to create resulting data list";

//-------------------------------------------------------------------------

class CPyParse {

	//Pointer to function returning text of error message associated with error number after a library function call.
	int (*m_pGetErrText)(char*, int); 

	//Auxiliary vector for storing Py array buffers (maybe not needed?)
	vector<Py_buffer> m_vBuf;

public:

	CPyParse(int (*pInGetErrText)(char*, int))
	{
		m_pGetErrText = pInGetErrText;
	}

	/************************************************************************//**
	 * Auxiliary function dedicated to process errors reported by Library
	 ***************************************************************************/
	void ProcRes(int er) //throw(...) 
	{
		char sErrBuf[2048];

		if(er == 0) return;
		else
		{
			//srwlUtiGetErrText(ErrorBuf, er);
			(*m_pGetErrText)(sErrBuf, er);
			if(er < 0)
			{//Print Warning:
				PyErr_SetString(PyExc_Warning, sErrBuf);
				PyErr_PrintEx(1); //?
			}
			else throw sErrBuf;
		}
	}

	/************************************************************************//**
	 * Gets access to Py array buffer, and returns it pointer to it
	 ***************************************************************************/
	//char* GetPyArrayBuf(PyObject* obj, vector<Py_buffer>* pvBuf, Py_ssize_t* pSizeBuf) //sizeBuf is out
	char* GetPyArrayBuf(PyObject* obj, Py_ssize_t* pSizeBuf =0) //sizeBuf is out
	{//for simplicity and uniformity of treatment in Py3 and Py2, only writable buffers are supported
		if(obj == 0) return 0;
		if(PyObject_CheckBuffer(obj))
		{
			Py_buffer pb_tmp;
			//if(PyObject_GetBuffer(obj, &pb_tmp, bufFlag)) return 0;
			if(PyObject_GetBuffer(obj, &pb_tmp, PyBUF_WRITABLE)) return 0;
			if(pSizeBuf != 0) *pSizeBuf = pb_tmp.len;
			//if(pvBuf != 0) pvBuf->push_back(pb_tmp);
			m_vBuf.push_back(pb_tmp); //Maybe not necessary?
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
				if (pVoidBuffer != 0) pRes = (char*)pVoidBuffer;
			}
			Py_DECREF(pOldBuf);
			if(pSizeBuf != 0) *pSizeBuf = sizeBuf;
			return pRes;
		}
#endif
		return 0;
	}

	/************************************************************************//**
	 * Checks if Py object is List of Array
	 ***************************************************************************/
	static char CheckIfObjIsListOrArray(PyObject* obj, Py_buffer& pb, void*& pvb, Py_ssize_t& len)
	{
		pvb = 0; len = 0;
		if(PyList_Check(obj)) 
		{
			len = PyList_Size(obj);
			return 'l';
		}

		//if(PyObject_CheckBuffer(obj)) return 'a';
		bool isArray = PyObject_CheckBuffer(obj);

#if PY_MAJOR_VERSION >= 3

		if(!isArray) return 0;

#endif

		if(isArray)
		{
			if(PyObject_GetBuffer(obj, &pb, PyBUF_SIMPLE)) return 0;
			pvb = pb.buf;
			len = pb.len;
			return 'a';
		}

#if PY_MAJOR_VERSION < 3

		else
		{
			PyObject *pOldBuf = PyBuffer_FromReadWriteObject(obj, 0, Py_END_OF_BUFFER);
			if(pOldBuf != 0)
			{
				if(PyObject_AsWriteBuffer(pOldBuf, &pvb, &len)) return 0;
				else return 'a';
			}
			else
			{
				PyErr_Clear();
				return 0;
			}
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
	 * If the length of Py list or array is larger than nElem at input, then nElemTooSmall is set to true
	 ***************************************************************************/
	template<class T> static char CopyPyListElemsToNumArray(PyObject* obj, char arType, T*& ar, int& nElem, bool& nElemTooSmall)
	//template<class T> static char CopyPyListElemsToNumArray(PyObject* obj, char arType, T*& ar, int& nElem)
	{
		nElemTooSmall = false; //OC29072018

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
		PyObject *pOldBuf = 0;
		int *pIntAr = 0;
		long *pLongAr = 0;
		float *pFloatAr = 0;
		double *pDoubleAr = 0;

		int nElemInList = 0;
		if(isList) 
		{
			nElemInList = (int)PyList_Size(obj);
		}
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
			//Case of Array:
			if(arType == 'i')
			{
				//nElemInList = (int)pb.len/sizeof(int);
				//pIntAr = (int*)pb.buf;
				nElemInList = (int)(sizeBuf / sizeof(int));
				pIntAr = (int*)pVoidBuffer;
			}
			else if(arType == 'l')
			{
				//nElemInList = (int)pb.len/sizeof(long);
				//pLongAr = (long*)pb.buf;
				nElemInList = (int)(sizeBuf / sizeof(long));
				pLongAr = (long*)pVoidBuffer;
			}
			else if(arType == 'f')
			{
				//nElemInList = (int)pb.len/sizeof(float);
				//pFloatAr = (float*)pb.buf;
				nElemInList = (int)(sizeBuf / sizeof(float));
				pFloatAr = (float*)pVoidBuffer;
			}
			else if(arType == 'd')
			{
				nElemInList = (int)(sizeBuf / sizeof(double));
				pDoubleAr = (double*)pVoidBuffer;
			}
		}
		if(nElemInList < 0) throw strEr_BadListArray;
		else if(nElemInList == 0) return 0; //OC29062018 (empty lists are allowed)

		if(ar == 0)
		{
			ar = new T[nElemInList];
			nElem = nElemInList;
		}
		else
		{
			if(nElem > nElemInList) nElem = nElemInList;
			else if(nElem < nElemInList) nElemTooSmall = true; //OC29072018
		}

		T *t_ar = ar;
		for(int i = 0; i < nElem; i++)
		{
			if(isList)
			{
				PyObject *o = PyList_GetItem(obj, (Py_ssize_t)i);
				if(o == 0) throw strEr_BadNum;
				if(!PyNumber_Check(o)) throw strEr_BadNum;

				if((arType == 'i') || (arType == 'l')) *t_ar = (T)PyLong_AsLong(o);
				else if((arType == 'f') || (arType == 'd')) *t_ar = (T)PyFloat_AsDouble(o);

				//Py_DECREF(o); //Uncommenting possibly causes crash
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

		return isList ? 'l' : 'a'; //OC03092016
	}

	/************************************************************************//**
	 * Copies elements of Py list or array of known length to a numerical array
	 * The output T *ar is expected to be allocated in calling function
	 * arType can be 'i', 'f' or 'd'
	 * Supports both Py lists only
	 ***************************************************************************/
	template<class T> static void CopyPyListElemsToNumArrayKnownLen(PyObject* o, char arType, T* ar, int nElem, const char* sEr=0)
	{
		if(o == 0) return;

		bool lenIsSmall = false;
		int len = nElem;
		CPyParse::CopyPyListElemsToNumArray(o, arType, ar, len, lenIsSmall);
		if(sEr != 0)
		{
			if((len != nElem) || lenIsSmall) throw sEr;
		}
	}

	/************************************************************************//**
	 * Copies elements of Py list of lists or list of arrays or array to a vector
	 * arType can be 'i', 'f' or 'd'
	 * Supports both Py lists and arrays
	 ***************************************************************************/
	template<class T> static char CopyPyNestedListElemsToNumVect(PyObject* obj, char typeElem, vector<T>* pv)
	{
		if(obj == 0) return 0;
		if(!((typeElem == 'i') || (typeElem == 'l') || (typeElem == 'f') || (typeElem == 'd'))) return 0; 

		T num;
		Py_buffer pb;
		void *pVoidBuf = 0;
		Py_ssize_t lenListOrAr = 0;
		char l_or_a = CheckIfObjIsListOrArray(obj, pb, pVoidBuf, lenListOrAr);

		if(l_or_a == 0) 
		{//Try to parse a number
			if(!PyNumber_Check(obj)) throw strEr_BadNum;

			if((typeElem == 'i') || (typeElem == 'l')) num = (T)PyLong_AsLong(obj);
			else if((typeElem == 'f') || (typeElem == 'd')) num = (T)PyFloat_AsDouble(obj);
			pv->push_back(num);
			return 'n';
		}
		else if(l_or_a == 'a')
		{//Collect numbers from array
			int *pIntAr = 0;
			long *pLongAr = 0;
			float *pFloatAr = 0;
			double *pDoubleAr = 0;

			switch(typeElem) 
			{
				case 'd':
					pDoubleAr = (double*)pVoidBuf; break;
				case 'f':
					pFloatAr = (float*)pVoidBuf; break;
				case 'i':
					pIntAr = (int*)pVoidBuf; break;
				case 'l':
					pLongAr = (long*)pVoidBuf; break;
			}

			for(Py_ssize_t i=0; i<lenListOrAr; i++)
			{
				switch(typeElem) 
				{
					case 'd':
						num = (T)(*(pDoubleAr++)); break;
					case 'f':
						num = (T)(*(pFloatAr++)); break;
					case 'i':
						num = (T)(*(pIntAr++)); break;
					case 'l':
						num = (T)(*(pLongAr++)); break;
				}
				pv->push_back(num);
			}
			return 'a';
		}
		else if(l_or_a == 'l')
		{
			for(Py_ssize_t i=0; i<lenListOrAr; i++)
			{
				PyObject *oElem = PyList_GetItem(obj, (Py_ssize_t)i);
				if(oElem == 0) throw strEr_BadListArray;
				if(PyNumber_Check(oElem))
				{
					if((typeElem == 'i') || (typeElem == 'l')) num = (T)PyLong_AsLong(oElem);
					else if((typeElem == 'f') || (typeElem == 'd')) num = (T)PyFloat_AsDouble(oElem);
					pv->push_back(num);
				}
				else
				{
					char res = CopyPyNestedListElemsToNumVect(oElem, typeElem, pv);
					if(res == 0) return 0;
				}
				//Py_DECREF(oElem); //Uncommenting possibly causes crash
			}
			return 'l';
		}
		return 0;
	}

	/************************************************************************//**
	 * Copies elements of Py list of lists or list of arrays or array to a flat C-aligned numerical array
	 * ATTENTION: it can allocate T *ar !
	 * arType can be 'i', 'f' or 'd'
	 * Supports both Py lists and arrays
	 ***************************************************************************/
	template<class T> static char CopyPyNestedListElemsToNumAr(PyObject* obj, char typeElem, T*& ar, int& nElem)
	{
		vector<T> vData;
		char res = CopyPyNestedListElemsToNumVect(obj, typeElem, &vData);
		if(res == 0) return 0;

		int nElemAct = (int)vData.size();
		if(nElemAct <= 0) return 0;

		if(ar == 0)
		{
			ar = new T[nElemAct];
			nElem = nElemAct;
		}
		else
		{
			if(nElem > nElemAct) nElem = nElemAct;
		}

		T *t_ar = ar;
		for(int i=0; i<nElem; i++)
		{
			*(t_ar++) = vData[i];
		}
		return res;
	}

	/************************************************************************//**
	 * Find lengths of Py lists or arrays that are elements of a list
	 * ATTENTION: it can allocate int *arLen !
	 ***************************************************************************/
	static void FindLengthsOfElemListsOrArrays(PyObject* oList, int*& arLens, int& nElem)
	{
		if(oList == 0) throw strEr_BadList;
		if(!PyList_Check(oList)) throw strEr_BadList;

		int nElemAct = (int)(int)PyList_Size(oList);
		if(nElemAct <= 0) return;

		if(arLens == 0)
		{
			arLens = new int[nElemAct];
			nElem = nElemAct;
		}
		else
		{
			if(nElem > nElemAct) nElem = nElemAct;
		}

		int *t_arLens = arLens;
		for(int i=0; i<nElem; i++)
		{
			PyObject *o = PyList_GetItem(oList, (Py_ssize_t)i);
			if(o == 0) throw strEr_BadListArray;

			bool isList = PyList_Check(o);
			bool isArray = false;
			if(!isList) isArray = PyObject_CheckBuffer(o);

#if PY_MAJOR_VERSION >= 3

			if(!(isList || isArray)) throw strEr_BadListArray;

#endif

			int nElemCur = 0;
			if(isList) 
			{
				nElemCur = (int)PyList_Size(o);
			}
			else
			{
				void *pVoidBuffer = 0;
				Py_ssize_t sizeBuf = 0;

				Py_buffer pb;
				PyObject *pOldBuf = 0;

				if(isArray)
				{
					if(PyObject_GetBuffer(o, &pb, PyBUF_SIMPLE)) throw strEr_BadListArray;
					pVoidBuffer = pb.buf;
					sizeBuf = pb.len;
				}

#if PY_MAJOR_VERSION < 3

				else
				{
					pOldBuf = PyBuffer_FromReadWriteObject(o, 0, Py_END_OF_BUFFER);
					if(pOldBuf != 0)
					{
						if(PyObject_AsWriteBuffer(pOldBuf, &pVoidBuffer, &sizeBuf)) throw strEr_BadListArray;
						isArray = true;
					}
					else
					{
						PyErr_Clear();
						throw strEr_BadListArray;
					}
				}

#endif

				nElemCur = (int)sizeBuf;
			}
			*(t_arLens++) = nElemCur;
		}
	}

	/************************************************************************//**
	 * Copies elements of Py string to a C string (assumed to be allocated outside)
	 ***************************************************************************/
	static void CopyPyStringToC(PyObject* pObj, char* c_str, int maxLenStr)
	{
		//const char strEr_BadStr[] = "Error at parsing / converting Python string";
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
	 * Copies elements of Py byte array a C string (assumed to be allocated outside)
	 ***************************************************************************/
//	static void CopyPyByteArrToC(PyObject* pObj, char*& arBytes, int& nBytes)
//	{
//		if(pObj == 0) throw strEr_BadStr;
//
//		Py_ssize_t len = 0;
//
//#if PY_MAJOR_VERSION < 3
//		//if(PyString_AsStringAndSize(pObj, &arBytes, &len) == -1) throw strEr_BadStr;
//		if(PyBytes_AsStringAndSize(pObj, &arBytes, &len) == -1) throw strEr_BadStr;
//#else
//		//nBytes = (int)PyBytes_Size(pObj);
//		//if(nBytes <= 0) throw strEr_BadStr;
//		////arBytes = new char[nBytes];
//		//arBytes = PyBytes_AsString(pObj);
//
//		if(PyBytes_AsStringAndSize(pObj, &arBytes, &len) == -1) throw strEr_BadStr;
//#endif
//		nBytes = (int)len;
//	}

	/************************************************************************//**
	 * Copies char from C to Py
	 ***************************************************************************/
	static PyObject* Py_BuildValueChar(char inC)
	{
#if PY_MAJOR_VERSION >= 3
		return Py_BuildValue("C", inC); //doesn't work with Py2.7
#else
		return Py_BuildValue("c", inC);
#endif
	}

	/************************************************************************//**
	 * Copies char from C to Py
	 ***************************************************************************/
	static PyObject* Py_BuildValueByteStr(char* inStr, int len)
	{
#if PY_MAJOR_VERSION >= 3
		return Py_BuildValue("y#", inStr, len);
#else
		//return Py_BuildValue("c#", inStr, len);
		return Py_BuildValue("s#", inStr, len); //OC04102018
#endif
	}

	/************************************************************************//**
	 * Copies elements of Py string to a C string (assumed to be allocated outside)
	 ***************************************************************************/
	static void CopyPyClassNameToC(PyObject* pObj, char* c_str, int maxLenStr)
	{
		//const char strEr_BadClassName[] = "Error at retrieving Python class name";
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
	* Sets up output list (eventually of lists) data from an array
	***************************************************************************/
	template<class T> static PyObject* SetDataListOfLists(T* arB, int nB, int nP, const char* cType="d") //OC13092018
	//static PyObject* SetDataListOfLists(double* arB, int nB, int nP)
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
	* Sets up output list (eventually of lists) data from an array
	***************************************************************************/
	template<class T> static PyObject* SetDataListsNested(T*& ar, int* arDims, char* cType="d") //OC16092018
	//static PyObject* SetDataListOfLists(double* arB, int nB, int nP)
	{//More testing may be required!
		//arDims[0] has to define number of dimensions
		//arDims[1] is number of points in the in-most dimenstion
		//...
		//arDims[arDims[0] - 1] is number of points in the out-most dimenstion

		if((ar == 0) || (arDims == 0)) return 0;
		int nDims = arDims[0];
		
		int indNumElem = 1;
		for(int ii=nDims; ii>0; ii--)
		{
			if(arDims[ii] > 1)
			{
				indNumElem = ii; break;
			}
		}

		//int nElem = arDims[nDims_mi_1];
		int nElem = arDims[indNumElem];
		PyObject *oRes = PyList_New(nElem);

		//PyObject *oRes = 0; 
		//if(nElem > 1) oRes = PyList_New(nElem);
		//else return 0;
		//arDims[0] = nDims_mi_1;

		arDims[0] = indNumElem - 1;
		PyObject *oElem = 0;

		for(int i=0; i<nElem; i++)
		{
			if(indNumElem > 1) oElem = SetDataListsNested(ar, arDims, cType);
			else oElem = Py_BuildValue(cType, *(ar++));

			//if((oElem == 0) && (nDims == 1) && (arDims[nDims] == 1)) oElem = Py_BuildValue(cType, *(ar++));
			//if((oElem == 0) && (nDims_mi_2 >= 0) && (arDims[nDims_mi_2] == 1)) oElem = Py_BuildValue(cType, *(ar++));

			if(oElem != 0)
			{
				if(PyList_SetItem(oRes, (Py_ssize_t)i, oElem)) throw strEr_FailedCreateList;
			}
		}
		arDims[0] = nDims;
		return oRes;
	}
};

//-------------------------------------------------------------------------

#endif
