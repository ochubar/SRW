/************************************************************************//**
 * File: auxparse.h
 * Description: String manipulation utilities (header)
 * Project: 
 * First release: 2002
 *
 * Copyright (C) Synchrotron SOLEIL, France
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#ifndef __AUXPARSE_H
#define __AUXPARSE_H

#include <ctype.h>
#include <vector>
#include <string>
#include <cstring>

//#ifdef __GCC__
//#define std
//#else
using namespace std;
//#endif

//-------------------------------------------------------------------------

class CAuxParse {

public:

	static void toUpperCase(char* s)
	{
		if(s == 0) return;
		size_t len = strlen(s);
		if(len == 0) return;

		char *ts = s;
		for(unsigned i=0; i<len; i++)
		{
			*ts = (char)toupper(*ts);
			ts++;
		}
	}

	static char** StringArrayAllocate(int lenStr, long numStr)
	{
		if((lenStr <= 0) || (numStr <= 0)) return 0; 

		char** arStr = new char*[numStr];
		char **t_arStr = arStr;
		for(long k=0; k<numStr; k++) 
		{
			*t_arStr = new char[lenStr];
			**t_arStr = '\0';
			t_arStr++;
		}
		return arStr;
	}

	static void StringArrayDeallocate(char** arrStr, long numStr)
	{
		if((arrStr == 0) || (numStr <= 0)) return; 

		char **tArrStr = arrStr;
		for(long i=0; i<numStr; i++)
		{
			if(*tArrStr != 0) { delete[] *tArrStr; *tArrStr = 0;}
			tArrStr++;
		}
		delete[] arrStr; 
	}
	
	static void DoubleVect2Arr(vector<double>& vd, double*& ad)
	{//ATTENTION: this allocates array of double !

		int num = (int)vd.size();
		if(num <= 0) return;

		ad = new double[num];
		double *p_ad = ad;

		for(vector<double>::iterator it = vd.begin(); it != vd.end(); ++it)
		{
			*(p_ad++) = *it;
		}
	}
	
	template<class T> static int FindElemInd(T& cToFind, vector<T>& vcToSearch)
	{
		for(int i=0; i<(int)vcToSearch.size(); i++)
		{
			if(vcToSearch[i] == cToFind) return i;
		}
		return -1;
	}

	static int FindElemInd(char* cToFind, char** arToSearch, int size_arToSearch)
	{
		char **t_arToSearch = arToSearch;
		for(int i=0; i<size_arToSearch; i++)
		{
			if(strcmp(*(t_arToSearch++), cToFind) == 0) return i;
		}
		return -1;
	}

	template<class T> static long FindMaxElemInd(T* arToSearch, long len_arToSearch, T& maxElem)
	{
		if((arToSearch == 0) || (len_arToSearch <= 0)) return -1;
		T *t_arToSearch = arToSearch + 1;
		maxElem = arToSearch[0];
		long indMaxElem = 0;
		for(long i=1; i<len_arToSearch; i++)
		{
			if(maxElem < *t_arToSearch) 
			{
				maxElem = *t_arToSearch;
				indMaxElem = i;
			}
			t_arToSearch++;
		}
		return indMaxElem;
	}
	
	template<class T> static void FindMinMax(T* ar, long len, T& v_min, long& i_min, T& v_max, long& i_max)
	{
		i_min = i_max = -1;
		if((len <= 0) || (ar == 0)) return;

		T *t_ar = ar;
		//v_max = (T)(-1E+23), v_min = (T)(1E+23);
		v_max = v_min = ar[0];
		for(long i=0; i<len; i++) 
		{
			if(v_max < *t_ar)
			{
				v_max = *t_ar; i_max = i;
			}
			if(v_min > *t_ar)
			{
				v_min = *t_ar; i_min = i;
			}
			t_ar++;
		}
	}

	template<class T> static int FindLargestPartIndInVect2D(vector<vector<T> >& v, int& maxSize)
	{
		maxSize = 0;
		int indPart = -1;
		for(int i=0; i<(int)v.size(); i++)
		{
			int curSize = (int)(v[i].size());
			if(maxSize < curSize)
			{
				maxSize = curSize;
				indPart = i;
			}
		}
		return indPart;
	}
	
	template<class T> static void FindPartSizesVect2D(vector<vector<T> >& v2D, vector<int>& vPartSizes)
	{
		if(vPartSizes.size() > 0) vPartSizes.erase(vPartSizes.begin(), vPartSizes.end());
		for(int i=0; i<(int)v2D.size(); i++)
		{
			vPartSizes.push_back((int)v2D[i].size());
		}
	}
	template<class T> static int* FindPartSizesVect2D(vector<vector<T> >& v2D)
	{//ATTENTION: this allocates array !
		vector<int> vPartSizes;
		FindPartSizesVect2D(v2D, vPartSizes);
		return Vect2Ar(vPartSizes);
	}
	
	template<class T> static void ArraySubtractMin(T* ar, long len)
	{
		if((ar == 0) || (len <= 0)) return;

		T v_min, v_max;
		long i_min, i_max;
		FindMinMax(ar, len, v_min, i_min, v_max, i_max);

		if(v_min == 0) return;
		for(int i=0; i<len; i++) ar[i] -= v_min;
	}

	template<class T> static T* Vect2Ar(vector<T>& v)
	{//ATTENTION: this allocates array !
		long size_v = (long)v.size();
		if(size_v <= 0) return 0;
		T* ar = new T[size_v];
		T *t_ar = ar;
		for(long i=0; i<size_v; i++) *(t_ar++) = v[i];
		return ar;
	}

	template<class T> static T** Vect2Ar2D(vector<vector<T> >& vv, int*& arLen)
	{//ATTENTION: this allocates arrays !
		int size_vv = (int)vv.size();
		if(size_vv <= 0) return 0;
		T** ar = new T*[size_vv];
		arLen = new int[size_vv];
		T **t_ar = ar;
		int *t_arLen = arLen;
		for(int i=0; i<size_vv; i++) 
		{
			vector<T>& v = vv[i];
			int size_v = (int)v.size();
			*t_ar = (size_v > 0)? new T[size_v] : 0;
			*(t_arLen++) = size_v;
			T *t_ar2 = *(t_ar++);
			for(int j=0; j<size_v; j++) *(t_ar2++) = v[j];
		}
		return ar;
	}

	template<class T> static void DeallocAr2D(T**& ar2D, int len)
	{
		if((ar2D == 0) || (len <= 0)) return;
		T **t_ar2D = ar2D;
		for(int i=0; i<len; i++)
		{
			if(*t_ar2D != 0) { delete[] *t_ar2D; *t_ar2D = 0;}
			t_ar2D++;
		}
		delete[] ar2D;
		ar2D = 0;
	}

	static long FindSumVectElem(vector<int>& v)
	{
		long sum = 0;
		for(vector<int>::iterator it = v.begin(); it != v.end(); ++it)
		{
			sum += *it;
		}
		return sum;
	}
	
	static long FindSumArElem(int* ar, long lenAr)
	{
		long sum = 0;
		int *t_ar = ar;
		for(long i=0; i<lenAr; i++)
		{
			sum += *(t_ar++);
		}
		return sum;
	}

	template<class T> static long FindSumVectElemSizes(vector<vector<T> >& v2D)
	{
		long sum = 0;
		for(int i=0; i<(int)v2D.size(); i++)
		{
			sum += (int)v2D[i].size();
		}
		return sum;
	}
	
	static int FindVectLargestElem(vector<int>& v)
	{
		int max = v[0];
		for(vector<int>::iterator it = v.begin(); it != v.end(); ++it)
		{
			if(max < *it) max = *it;
		}
		return max;
	}
	
	static int FindIndVectSmallestElem(vector<double>& v, double& min)
	{
		int iMin = -1;
		int v_size = (int)v.size();
		if(v_size <= 0) return iMin;
		iMin = 0;
		min = v[0];
		for(int i=0; i<v_size; i++)
		{
			double curElem = v[i];
			if(min > curElem)
			{
				min = curElem; iMin = i;
			}
		}
		return iMin;
	}

	static int FindLargestElemSmallerThan(double s, double* ar, int lenAr)
	{
		if((ar == 0) || (lenAr <= 0)) return -1;
		
		double *t_ar = ar;
		double prev = *(t_ar++);
		for(int i=1; i<lenAr; i++)
		{
			if((prev < s) && (*t_ar > s)) return (i - 1);
			prev = *(t_ar++);
		}
		return -1;
	}

	static int FindLargestElemSmallerThanOrEqualTo(double s, double* ar, int lenAr)
	{
		if((ar == 0) || (lenAr <= 0)) return -1;
		
		double *t_ar = ar;
		double prev = *(t_ar++);
		for(int i=1; i<lenAr; i++)
		{
			if((prev <= s) && (*t_ar > s)) return (i - 1);
			prev = *(t_ar++);
		}
		if(prev == s) return (lenAr - 1); //?
		return -1;
	}

	static int FindSmallestElemLargerThan(double s, double* ar, int lenAr)
	{
		if((ar == 0) || (lenAr <= 0)) return -1;
	
		double *t_ar = ar + lenAr - 1;
		double prev = *(t_ar--);
		
		for(int i=(lenAr - 2); i>=0; i--)
		{
			if((prev > s) && (*t_ar < s)) return (i - 1);
			prev = *(t_ar--);
		}
		return -1;
	}

	template<class T> static int FindIndVectElemEqualToWithinTol(T valueToFind, vector<T>& v, T& absTol)
	{
		int v_size = (int)v.size();
		if(v_size <= 0) return -1;
		for(int i=0; i<v_size; i++)
		{
			T curElemValDif = v[i] - valueToFind;
			if((curElemValDif >= -absTol) && (curElemValDif <= absTol)) return i;
		}
		return -1;
	}

	static bool CheckIfStringExistsInArray(char* str, char** arStr, int lenAr)
	{
		if((str == 0) || (arStr == 0) || (lenAr <= 0)) return false;
		for(int i=0; i<lenAr; i++)
		{
			if(strcmp(str, arStr[i]) == 0) return true;
		}
		return false;
	}

	template<class T> static bool CheckIfElemIsPresent(const T& elem, const vector<T>& vectOld)
	{
		long numVectOld = (long)vectOld.size();
		for(long i=0; i<numVectOld; i++)
		{
			if(elem == vectOld[i]) return true;
		}
		return false;
	}

	template<class T> static bool CheckIfVectElemArePresent(const vector<T>& vectNew, const vector<vector<T> >& vectVectOld)
	{
		long numVectOld = (long)vectVectOld.size();
		long numElemVectNew = (long)vectNew.size();

		for(long i=0; i<numVectOld; i++)
		{
			const vector<T> &vectOld = vectVectOld[i];
			bool allElemArePresentInCurOldVect = true;
			for(long j=0; j<numElemVectNew; j++)
			{
				if(!CheckIfElemIsPresent(vectNew[j], vectOld))
				{
					allElemArePresentInCurOldVect = false;
					break;
				}
			}
			if(allElemArePresentInCurOldVect) return true;
		}
		return false;
	}

	//static int CalcAbsFlatIndFromPartIndAndRelInd(int* arNumMagDifTypes, int numMagTypes, int indCurMagType, int relMagIndInType)
	static int CalcAbsFlatIndFromPartIndAndRelInd(int indPart, int relInd, int* arPartLengths, int numParts)
	{
		if((arPartLengths == 0) || (numParts <= 0) || (indPart < 0) || (relInd < 0)) return -1;
		if(indPart >= numParts) return -1;
		if(relInd >= arPartLengths[indPart]) return -1;
		int absCount = 0;
		for(int i=0; i<indPart; i++) absCount += arPartLengths[i];
		absCount += relInd;
		return absCount;
	}

	static int CalcPartIndAndRelIndFromAbsFlatInd(int indAbs, int* arPartLengths, int numParts, int& indPart, int& relInd)
	{//to test !!!
		indPart = relInd = -1;
		if((arPartLengths == 0) || (numParts <= 0)) return relInd;
		if(indAbs == 0)
		{
			indPart = relInd = 0; return 0;
		}
		int *t_arPartLengths = arPartLengths;
		int sumPartLengths = *(t_arPartLengths++);
		for(int iPart=1; iPart<numParts; iPart++)
		{
			if(sumPartLengths <= indAbs) 
			{
				sumPartLengths += *t_arPartLengths;
			}
			else 
			{
				indPart = iPart - 1;
				sumPartLengths -= *(t_arPartLengths - 1);
				relInd = indAbs - sumPartLengths;
				return relInd;
			}
			t_arPartLengths++;
		}
		if(sumPartLengths > indAbs)
		{
			indPart = numParts - 1;
			sumPartLengths -= *(t_arPartLengths - 1);
			relInd = indAbs - sumPartLengths;
		}

		return relInd;
	}

	static char* PrependZeros(char* strNum, int reqLength)
	{
		if((strNum == 0) || (reqLength <= 0)) return 0;

		int origStrLen = (int)strlen(strNum);
		if(origStrLen >= reqLength) return strNum;

		char *strBuf = new char[reqLength + 1];
		*strBuf = '\0';
		for(int i=origStrLen; i<reqLength; i++)
		{
			strcat(strBuf, "0");
		}
		strcat(strBuf, strNum);
		strcpy(strNum, strBuf);
		delete[] strBuf;
		return strNum;
	}

	static void StringSplitSimple(const char* c_strTot, char cSep, vector<string>& vsSepRes)
	{
		if(c_strTot == NULL) return;
		StringSplitSimple(c_strTot, (int)strlen(c_strTot), cSep, vsSepRes);
	}

	static void StringArr2Vect(char** as, int numStr, vector<string>& vs)
	{
		if((as == 0) || (numStr <= 0)) return;
		for(int i=0; i<numStr; i++)
		{
			string curStr(as[i]);
			vs.push_back(curStr);
		}
	}

	static void StringArr2VectCStr(char** as, int numStr, vector<char*, allocator<char*> >& vCStr)
	{
		if((as == 0) || (numStr == 0)) return;
		for(int i=0; i<numStr; i++) 
		{
			char* CurStr = as[i];
			vCStr.push_back(CurStr);
		}
	}

	template<class T> static void ZeroArr(T* ar, unsigned len_ar)
	{
		T *t_ar = ar;
		for(unsigned i=0; i<len_ar; i++)
		{
			*(t_ar++) = (T)0;
		}
	}

	template<class T> static int Round(T x)
	{
		T ix = floor(x);
		T dx = x - ix;
		return (dx < 0.5)? (int)ix : (int)(ix + 1);
	}

	template<class T1, class T2> static bool LessInPairBasedOnSecond(pair<T1, T2> p1, pair<T1, T2> p2)
	{
		return p1.second < p2.second;
	}
	template<class T1, class T2, class T3> static bool LessInPairBasedOnFirstInNestedPair(pair<T1, pair<T2, T3> > p1, pair<T1, pair<T2, T3> > p2)
	{
		return p1.second.first < p2.second.first;
	}

	static void StringSplitSimple(const char* c_strTot, int maxNumSymb, char cSep, vector<string>& vsSepRes);
	static void StringSplit(char* TotStr, char** SepStrArr, int AmOfSepStr, char* SymbToCut, vector<string>& ResStrings);
	static void StringSymbolsRemove(const char* OrigStr, char* SymbToCut, char* FinStr);
	static void StringSymbolsRemoveAtBegin(const char* OrigStr, char* SymbToCut, char* FinStr);
	static void StringArrayRemoveEmptyElemFromTable(char** arrPlacesTableFlat, int numRows, int numCols, char**& arrPlacesFlat, int*& arrNonEmpty, long& totNumNonEmptyElem);
	static long StringArrayCountNonEmptyElem(char** arrStr, long len);
	static void StringTokenize(const char* in_str, const char* sepChars, vector<string>& vTokens);
	static void StringFindTokens(const char* in_str, const char* sepChars, vector<pair<int, int> >& vIndTokens);
	static void StringRemoveCharsFromStartAndEnd(char*& strToken, const char* chars);
	static void StringVect2Arr(vector<string>& vs, char**& as); //ATTENTION: this allocates array of c strings !
	static char** StringVect2Arr(vector<string>& vs);
	static void StringVect2ArrNoAlloc(vector<string>& vs, char** as, int& numStr);

	static void StringFindDifFirstLetters(vector<string>& vs, vector<char>& vcFirstLetters);

	static void StringRemoveOutDecor(char* str);
	static void StringRemoveQuotesSepBySpaces(char* strToken, const char* arSep);
	static bool StringCheckIfIncludesOnlyGivenChars(char* ar, int len_ar, char ch);
	static void StringFindIndCharBetweenQuotes(const char* str, char ch, vector<int>& vIndChars);
	static void StringMergeTokensAroundGivenCharInd(vector<pair<int, int> >& vIndTokens, vector<int>& vIndChars);
};

//-------------------------------------------------------------------------
//Auxiliary class for binary data collection / storage / manipulations (added for programming persistence in Radia)
class CAuxBinStr : public string {

private:
	string::iterator itOut;

public:

	void setItOut(long ofst=0)
	{//To call before using ">>" first time!
		itOut = begin() + ofst;
	}

	template<class T> void setFromPos(long ofst, T a)
	{//Assumes that memory has been allocated previously
		int sz = (int)sizeof(T);
		//char *p = reinterpret_cast<char*>(&a);
		unsigned char *p = reinterpret_cast<unsigned char*>(&a);
		string::iterator it = begin() + ofst;
		for(int i=0; i<sz; i++) *(it++) = *(p++); 
	}

	template<class T> inline friend CAuxBinStr& operator<<(CAuxBinStr& str, T a);
	template<class T> inline friend CAuxBinStr& operator>>(CAuxBinStr& str, T& a);
};

template<class T> inline CAuxBinStr& operator<<(CAuxBinStr& str, T a)
{
	int sz = (int)sizeof(T);
	//char *p = reinterpret_cast<char*>(&a);
	unsigned char *p = reinterpret_cast<unsigned char*>(&a);
	for(int i=0; i<sz; i++) str.push_back(*(p++)); //str += *(p++);
	return str;
};

template<class T> inline CAuxBinStr& operator>>(CAuxBinStr& str, T& a)
{//Extracts starting from position specified by itOut (assuming it exists)
	int sz = (int)sizeof(T);
	//char *p = reinterpret_cast<char*>(&a);
	unsigned char *p = reinterpret_cast<unsigned char*>(&a);
	string::iterator &it = str.itOut;
	for(int i=0; i<sz; i++) *(p++) = *(it++); 
	return str;
};

//-------------------------------------------------------------------------

#endif
