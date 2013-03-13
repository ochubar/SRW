/************************************************************************//**
 * File: auxparse.cpp
 * Description: String manipulation utilities
 * Project: 
 * First release: 2002
 *
 * Copyright (C) Synchrotron SOLEIL, France
 * All Rights Reserved
 *
 * @author O.Chubar
 * @version 1.0
 ***************************************************************************/

#include "auxparse.h"

//-------------------------------------------------------------------------

void CAuxParse::StringSplitSimple(const char* c_strTot, int maxNumSymb, char cSep, vector<string>& vsSepRes)
{
	if(c_strTot == 0) return;
	int lenStrTot = (int)strlen(c_strTot);
	if(lenStrTot == 0) return;
	if(lenStrTot > maxNumSymb) lenStrTot = maxNumSymb;
	
	char *sBufStrTot = new char[lenStrTot + 1];
	char *sBuf = new char[lenStrTot + 1];
	strncpy(sBufStrTot, c_strTot, lenStrTot);
	*(sBufStrTot + lenStrTot) = '\0';
	
	char *pSep, *pCurSubStr = sBufStrTot;
	for(;;)
	{
		pSep = strchr(pCurSubStr, cSep);
		if(pSep != NULL)
		{
			int lenNewPart = (int)(pSep - pCurSubStr);
			int ofstNewPart = (int)(pSep - sBufStrTot);

			if((lenNewPart > 0) && (ofstNewPart < maxNumSymb))
			{
				strncpy(sBuf, pCurSubStr, lenNewPart);
				sBuf[lenNewPart] = '\0';
				string sNewPart(sBuf);
				vsSepRes.push_back(sNewPart);
				pCurSubStr = pSep + 1;
			}
			else break;
		}
		else break;
	}
	int len_pCurSubStr = (int)strlen(pCurSubStr);
	if(len_pCurSubStr > 0)
	{
		int maxNumSymbToCopy = maxNumSymb - (int)(pCurSubStr - sBufStrTot);
		if(maxNumSymbToCopy > 0)
		{
			if(maxNumSymbToCopy > len_pCurSubStr) maxNumSymbToCopy = len_pCurSubStr;

			strncpy(sBuf, pCurSubStr, maxNumSymbToCopy);
			*(sBuf + maxNumSymbToCopy) = '\0';

			string sNewPart(sBuf);
			vsSepRes.push_back(sNewPart);
		}
	}

	if(sBufStrTot != NULL) delete[] sBufStrTot;
	if(sBuf != NULL) delete[] sBuf;
}

//-------------------------------------------------------------------------

void CAuxParse::StringSplit(char* TotStr, char** SepStrArr, int AmOfSepStr, char* SymbToCut, vector<string>& ResStrings)
{
	if(TotStr == 0) return;
	long LenTotStr = (long)strlen(TotStr);
	if(LenTotStr <= 0) return;

	if((SepStrArr == 0) || (AmOfSepStr == 0)) return;

	char *StrBuf = new char[LenTotStr + 1];
	char *StrBuf1 = new char[LenTotStr + 1];

	char *pStartStr = TotStr;
	long LenCurStr = (long)strlen(TotStr);
	for(long i=0; i<LenTotStr; i++)
	{
		char StopLoopNow = 0;
		for(int k=0; k<AmOfSepStr; k++)
		{
			char *SepStr = SepStrArr[k];
			if(SepStr == 0) continue;
            long LenSepStr = (long)strlen(SepStr);
            if(LenSepStr <= 0) continue;

			char *pCurSubStr = strstr(pStartStr, SepStr);
			if(pCurSubStr != NULL)
			{
				//char ActChar = *pCurSubStr;
				*pCurSubStr = '\0';
				strcpy(StrBuf, pStartStr);
				StringSymbolsRemove(StrBuf, SymbToCut, StrBuf1);

				string CurSubStr(StrBuf1);
				ResStrings.push_back(CurSubStr);

				long NewOffset = (long)(strlen(pStartStr) + LenSepStr);
				if(NewOffset >= LenCurStr) 
				{
					StopLoopNow = 1; break;
				}

				pStartStr += NewOffset;
				LenCurStr = (long)strlen(pStartStr);
			}
			else 
			{
				StringSymbolsRemove(pStartStr, SymbToCut, StrBuf1);
				string CurSubStr(StrBuf1);
				ResStrings.push_back(CurSubStr);
				StopLoopNow = 1; 
				break;
			}
		}
		if(StopLoopNow) break;
	}
	if(StrBuf != NULL) delete[] StrBuf;
	if(StrBuf1 != NULL) delete[] StrBuf1;
}

//-------------------------------------------------------------------------
//FinStr should be allocated in the calling function.
//To be safe, ensure length of FinStr to be strlen(OrigStr) + 1
//-------------------------------------------------------------------------
void CAuxParse::StringSymbolsRemove(const char* OrigStr, char* SymbToCut, char* FinStr)
{
	if(OrigStr == 0) return;
	long LenOrigStr = (long)strlen(OrigStr);
	if(LenOrigStr == 0) return;

	if(FinStr == 0) return;

	if(SymbToCut == 0) { strcpy(FinStr, OrigStr); return;}
	long LenSymbToCut = (long)strlen(SymbToCut);
	if(LenSymbToCut == 0) { strcpy(FinStr, OrigStr); return;}

	long CopiedSymbCount = 0;
	const char *tOrigStr = OrigStr;
	char *tFinStr = FinStr;
	for(long i=0; i<LenOrigStr; i++)
	{
		char *tSymbToCut = SymbToCut;
		char CopyingIsNotNeeded = 0;
		for(long j=0; j<LenSymbToCut; j++)
		{
			if(*tOrigStr == *tSymbToCut)
			{
				CopyingIsNotNeeded = 1; break;
			}
			tSymbToCut++;
		}
		if(!CopyingIsNotNeeded) 
		{
            *(tFinStr++) = *tOrigStr;
			CopiedSymbCount++;
		}
		tOrigStr++;
	}
	*tFinStr = '\0';
}

//-------------------------------------------------------------------------

void CAuxParse::StringSymbolsRemoveAtBegin(const char* OrigStr, char* SymbToCut, char* FinStr)
{
	if(OrigStr == 0) return;
	long LenOrigStr = (long)strlen(OrigStr);
	if(LenOrigStr == 0) return;
	if(FinStr == 0) return;

	long LenSymbToCut = (SymbToCut == 0)? 0 : (long)strlen(SymbToCut);
	if(LenSymbToCut <= 0) { strcpy(FinStr, OrigStr); return;}

	bool MainCopyingStarted = false;
	long CopiedSymbCount = 0;
	const char *tOrigStr = OrigStr;
	char *tFinStr = FinStr;
	for(long i=0; i<LenOrigStr; i++)
	{
		char CopyingIsNotNeeded = 0;
		if(!MainCopyingStarted)
		{
			char *tSymbToCut = SymbToCut;
			for(long j=0; j<LenSymbToCut; j++)
			{
				if(*tOrigStr == *tSymbToCut)
				{
					CopyingIsNotNeeded = 1; break;
				}
				tSymbToCut++;
			}
		}
		if(!CopyingIsNotNeeded) 
		{
            *(tFinStr++) = *tOrigStr;
			CopiedSymbCount++;
			MainCopyingStarted = true;
		}
		tOrigStr++;
	}
	*tFinStr = '\0';
}

//-------------------------------------------------------------------------

long CAuxParse::StringArrayCountNonEmptyElem(char** arrStr, long len)
{
	if((arrStr == 0) || (len <= 0)) return 0;

	long counter = 0;
	char **tArrStr = arrStr;
	for(long i=0; i<len; i++)
	{
		if((*tArrStr != 0) && (**tArrStr != '\0') && (**tArrStr != ' ')) counter++;
		tArrStr++;
	}
	return counter;
}

//-------------------------------------------------------------------------

void CAuxParse::StringArrayRemoveEmptyElemFromTable(char** arrPlacesTableFlat, int numRows, int numCols, char**& arrNonEmptyFlat, int*& arrLen, long& totNonEmptyElem)
{
	if((arrPlacesTableFlat == 0) || (numRows <= 0) || (numCols <= 0)) return;
	if(arrLen == 0) arrLen = new int[numCols];

	totNonEmptyElem = 0;
	char** tArrPlacesTableFlat = arrPlacesTableFlat;
	for(int i=0; i<numCols; i++)
	{
		long curNonEmptyElem = StringArrayCountNonEmptyElem(tArrPlacesTableFlat, numRows);
		arrLen[i] = curNonEmptyElem;
        totNonEmptyElem += curNonEmptyElem;
		tArrPlacesTableFlat += numRows;
	}
	if(totNonEmptyElem <= 0) return;

    arrNonEmptyFlat = new char*[totNonEmptyElem];
	char** tArrNonEmptyFlat = arrNonEmptyFlat;
	tArrPlacesTableFlat = arrPlacesTableFlat;
	for(int j=0; j<numCols; j++)
	{
        for(int k=0; k<numRows; k++)
		{
			if((*tArrPlacesTableFlat != 0) && (**tArrPlacesTableFlat != '\0') && (**tArrPlacesTableFlat != ' '))
			{
				*tArrNonEmptyFlat = new char[strlen(*tArrPlacesTableFlat) + 10];
				strcpy(*tArrNonEmptyFlat, *tArrPlacesTableFlat);
                tArrNonEmptyFlat++;
			}
			//else *tArrNonEmptyFlat = 0;
			tArrPlacesTableFlat++;
		}
	}
}
// ----------------------------------------------------------------------------

void CAuxParse::StringFindTokens(const char* in_str, const char* sepChars, vector<pair<int, int> >& vIndTokens)
{
	if((in_str == 0) || (sepChars == 0)) return;
	int len_in_str = (int)strlen(in_str);
	int len_sepChars = (int)strlen(sepChars);
	if((len_in_str <= 0) || (len_sepChars <= 0)) return;

	const char *t_in_str = in_str;
	int curStartInd = -1, curEndInd = -1;
	for(int i=0; i<len_in_str; i++)
	{
		bool isSeparChar = false;
		const char *t_sepChars = sepChars;
		for(int j=0; j<len_sepChars; j++)
		{
			if(*t_in_str == *(t_sepChars++))
			{
				isSeparChar = true; break;
			}
		}
		if(isSeparChar)
		{
			if((curStartInd >= 0) && (curEndInd < 0))
			{
				curEndInd = i;
				pair<int, int> newStartEnd(curStartInd, curEndInd - curStartInd);
				vIndTokens.push_back(newStartEnd);
				curStartInd = -1; curEndInd = -1;
			}
		}
		else
		{
			if(curStartInd < 0) curStartInd = i;
		}
		t_in_str++;
	}
	if((curStartInd >= 0) && (curEndInd < 0))
	{
		pair<int, int> newStartEnd(curStartInd, len_in_str - curStartInd);
		vIndTokens.push_back(newStartEnd);
	}
}

// ----------------------------------------------------------------------------

void CAuxParse::StringTokenize(const char* in_str, const char* sepChars, vector<string>& vTokens)
{
	vector<pair<int, int> > vIndTokens;
	if((sepChars == 0) || (strlen(sepChars) <= 0)) 
	{
		if((in_str != 0) && (strlen(in_str) > 0))
		{
			vTokens.push_back(string(in_str));
		}
		return;
	}

	StringFindTokens(in_str, sepChars, vIndTokens);
	int numTokens = (int)vIndTokens.size();
	if(numTokens <= 0) return;

	for(int i=0; i<numTokens; i++)
	{
		pair<int, int>& curTokenInd = vIndTokens[i];
		string newToken(in_str + curTokenInd.first, curTokenInd.second);

		//removing eventual spaces from beginning of string:
		const char *OrigStr = newToken.c_str();
		if(*OrigStr == ' ')
		{
			int lenOrigStr = (int)strlen(OrigStr);
			char *NewStr = new char[lenOrigStr + 1];
			//CAuxParse::StringSymbolsRemoveAtBegin(OrigStr, " \0", NewStr);
			char strAux[] = " \0";
			CAuxParse::StringSymbolsRemoveAtBegin(OrigStr, strAux, NewStr);
			string cleanedNewToken(NewStr);
			newToken = cleanedNewToken;
			delete[] NewStr;
		}

		vTokens.push_back(newToken);
	}
}

// ----------------------------------------------------------------------------

void CAuxParse::StringRemoveCharsFromStartAndEnd(char*& strToken, const char* chars)
{//this doesn't allocate
	if((strToken == 0) || (chars == 0)) return;
	int len_strToken = (int)strlen(strToken);
	int len_chars = (int)strlen(chars);
	if((len_strToken <= 0) || (len_chars <= 0)) return;

	char *t_strToken = strToken;
	int offset = 0;
	for(int i=0; i<len_strToken; i++)
	{
		bool charShouldBeRemoved = false;
		const char *t_chars = chars;
		for(int j=0; j<len_chars; j++)
		{
			if(*t_strToken == *(t_chars++))
			{
				charShouldBeRemoved = true; break;
			}
		}
		if(charShouldBeRemoved) t_strToken++;
		else
		{
			offset = i; break;
		}
	}
	if(offset > 0) 
	{
		strToken += offset;
		len_strToken = (int)strlen(strToken);
	}

	int len_strToken_mi_1 = len_strToken - 1;
	offset = len_strToken_mi_1;
	t_strToken = strToken + len_strToken_mi_1;
	for(int i=0; i<len_strToken; i++)
	{
		bool charShouldBeRemoved = false;
		const char *t_chars = chars;
		for(int j=0; j<len_chars; j++)
		{
			if(*t_strToken == *(t_chars++))
			{
				charShouldBeRemoved = true; break;
			}
		}
		if(charShouldBeRemoved) t_strToken--;
		else
		{
			offset -= i; break;
		}
	}
	if(offset != len_strToken_mi_1)
	{
		*(strToken + offset) = '\0';
	}
}

//-------------------------------------------------------------------------

void CAuxParse::StringVect2Arr(vector<string>& vs, char**& as)
{//ATTENTION: this allocates array of c strings !

	int numStr = (int)vs.size();
	if(numStr <= 0) return;
	
	as = new char*[numStr];
	char **p_as = as;
	
	for(vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	{
		*p_as = NULL;
		int curSize = (int)(it->size());
		if(curSize > 0)
		{
			*p_as = new char[curSize + 1];
			strcpy(*p_as, it->c_str());
		}
		p_as++;
	}
}

//-------------------------------------------------------------------------

void CAuxParse::StringVect2ArrNoAlloc(vector<string>& vs, char** as, int& numStr)
{//ATTENTION: this DOES NOT allocate array of c strings !

	numStr = (int)vs.size();
	if(numStr <= 0) return;
	
	//as = new char*[numStr];
	//char **p_as = as;
	
	//for(vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	for(int i=0; i<numStr; i++)
	{
		//*p_as = NULL;
		//**p_as = '\0';
		*(as[i]) = '\0';
		//int curSize = (int)(it->size());
		string& curString = vs[i];
		int curSize = (int)curString.size();
		if(curSize > 0)
		{
			//*p_as = new char[curSize + 1];
			//const char* curStr = it->c_str();
			const char* curStr = curString.c_str();
			int lenCurStr = (int)strlen(curStr);
			strncpy(as[i], curStr, lenCurStr);
			*(as[i] + lenCurStr) = '\0';
		}
		//p_as++;
	}
}

//-------------------------------------------------------------------------

char** CAuxParse::StringVect2Arr(vector<string>& vs)
{//ATTENTION: this allocates array of c strings !

	int numStr = (int)(vs.size());
	if(numStr <= 0) return 0;
	
	char **as = new char*[numStr];
	char **p_as = as;
	
	for(vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	{
		*p_as = NULL;
		int curSize = (int)(it->size());
		if(curSize > 0)
		{
			*p_as = new char[curSize + 1];
			strcpy(*p_as, it->c_str());
		}
		p_as++;
	}
	return as;
}

//-------------------------------------------------------------------------

void CAuxParse::StringFindDifFirstLetters(vector<string>& vs, vector<char>& vcFirstLetters)
{
	if(vs.size() <= 0) return;
	for(vector<string>::iterator it = vs.begin(); it != vs.end(); ++it)
	{
		char curFirstLetter = (*it)[0];
		bool isAlreadyThere = false;
		for(vector<char>::iterator it0 = vcFirstLetters.begin(); it0 != vcFirstLetters.end(); ++it0)
		{
			if(curFirstLetter == *it0)
			{
				isAlreadyThere = true;
				break;
			}
		}
		if(!isAlreadyThere) vcFirstLetters.push_back(curFirstLetter);
	}
}
// ----------------------------------------------------------------------------

void CAuxParse::StringRemoveOutDecor(char* str)
{//remove quotes and \n at the end from a string
	if(str == 0) return;
	int origStrLen = (int)strlen(str);
	if(origStrLen < 2) return;

	int numCharsToCopy = origStrLen - 1;
	bool thereIsCR = false;
	if(str[origStrLen - 2] == 10) 
	{
		numCharsToCopy--;
		thereIsCR = true;
	}

	char *t_str_new = str, *t_str_orig = str;
	if((*t_str_orig == '\'') || (*t_str_orig == '\"')) 
	{
		t_str_orig++; numCharsToCopy--;
	}

	for(int i=0; i<numCharsToCopy; i++)
	{
		*(t_str_new++) = *(t_str_orig++);
	}
	if((thereIsCR && (*t_str_orig == 10)) || ((*t_str_orig == '\'') || (*t_str_orig == '\"')))
	{
		*t_str_new = '\0';
	}
}

// ----------------------------------------------------------------------------
//removes '    ' or "   " from string
void CAuxParse::StringRemoveQuotesSepBySpaces(char* strToken, const char* arSep)
{
	if(strToken == 0) return;
	int len_strToken = (int)strlen(strToken);
	if(len_strToken <= 0) return;
	int numSepChar = (int)strlen(arSep);
	if(numSepChar <= 0) return;

	int maxNumGaps = 100;

	for(int j=0; j<numSepChar; j++)
	{
		char sepChar = arSep[j];

		char *p_strToken = strToken;
		for(int i=0; i<maxNumGaps; i++)
		{
			char *pFirstQuote = strchr(p_strToken, sepChar);
			if(pFirstQuote == NULL) break;

			int loc_len_strToken = (int)strlen(p_strToken);
			int loc_len_strToken_mi_1 = loc_len_strToken - 1;
			if(pFirstQuote >= (p_strToken + loc_len_strToken_mi_1)) throw "Syntax error: incorrect string format";
			
			char *pSecondQuote = strchr(pFirstQuote + 1, sepChar);
			if((pSecondQuote == NULL) || (pSecondQuote >= (p_strToken + loc_len_strToken_mi_1))) throw "Syntax error: incorrect string format";

			int lenInterv = (int)(pSecondQuote - pFirstQuote - 1);
			if(!StringCheckIfIncludesOnlyGivenChars(pFirstQuote + 1, lenInterv, ' ')) throw "Syntax error: incorrect string format";

			int numCharsAfterGap = (int)strlen(pSecondQuote + 1);
			char *t_main = pFirstQuote, *t_to_move = pSecondQuote + 1;
			for(int k=0; k<numCharsAfterGap; k++) *(t_main++) = *(t_to_move++);
			*t_main = '\0';

			p_strToken = pFirstQuote;
		}
	}
}

// ----------------------------------------------------------------------------

bool CAuxParse::StringCheckIfIncludesOnlyGivenChars(char* ar, int len_ar, char ch)
{
	if((ar == 0) || (len_ar <= 0)) return false;

	char *t_ar = ar;
	for(int i=0; i<len_ar; i++)
	{
		if(*(t_ar++) != ch) return false;
	}
	return true;
}
// ----------------------------------------------------------------------------

void CAuxParse::StringFindIndCharBetweenQuotes(const char* str, char ch, vector<int>& vIndChars)
{
	if((str == 0) || (*str == '\0')) return;

	char quoteChar = '\'', anotherQuoteChar = '\"';
	int len_str = (int)strlen(str);
	
	vector<int> vAuxIndCharsFound;
	bool firstQuoteCharWasFound = false, firstAnotherQuoteCharWasFound = false;
	bool secondQuoteCharWasFound = false, secondAnotherQuoteCharWasFound = false;
	int indFirstQuoteChar = -1, indSecondQuoteChar = -1;
	int indFirstAnotherQuoteChar = -1, indSecondAnotherQuoteChar = -1;
	char firstExternQuoteChar = 0;
	const char *t_str = str;
	for(int i=0; i<len_str; i++)
	{
		if((*t_str == quoteChar) && ((i == 0) || ((i > 0) && (*(t_str - 1) != '\\'))))
		{
			if(!firstQuoteCharWasFound) 
			{
				firstQuoteCharWasFound = true;
				indFirstQuoteChar = i;
				if(!firstAnotherQuoteCharWasFound) firstExternQuoteChar = quoteChar;
			}
			else 
			{
				secondQuoteCharWasFound = true;
				indSecondQuoteChar = i;
			}
		}
		else if((*t_str == anotherQuoteChar) && ((i == 0) || ((i > 0) && (*(t_str - 1) != '\\'))))
		{
			if(!firstAnotherQuoteCharWasFound) 
			{
				firstAnotherQuoteCharWasFound = true;
				indFirstAnotherQuoteChar = i;
				if(!firstQuoteCharWasFound) firstExternQuoteChar = anotherQuoteChar;
			}
			else 
			{
				secondAnotherQuoteCharWasFound = true;
				indSecondAnotherQuoteChar = i;
			}
		}
		else if(*t_str == ch)
		{
			vAuxIndCharsFound.push_back(i);
		}

		bool qouteCharPairWasFound = firstQuoteCharWasFound && secondQuoteCharWasFound && (firstExternQuoteChar == quoteChar);
		bool anotherQouteCharPairWasFound = firstAnotherQuoteCharWasFound && secondAnotherQuoteCharWasFound && (firstExternQuoteChar == anotherQuoteChar);

		if(qouteCharPairWasFound || anotherQouteCharPairWasFound)
		{
			int numCharsFound = (int)vAuxIndCharsFound.size();
			if(numCharsFound > 0)
			{
				for(int j=0; j<numCharsFound; j++) 
				{
					int curIndChar = vAuxIndCharsFound[j];
					if((qouteCharPairWasFound && ((indFirstQuoteChar < curIndChar) && (curIndChar < indSecondQuoteChar))) ||
					   (anotherQouteCharPairWasFound && ((indFirstAnotherQuoteChar < curIndChar) && (curIndChar < indSecondAnotherQuoteChar))))
					{
						vIndChars.push_back(curIndChar);
					}
				}
				vAuxIndCharsFound.erase(vAuxIndCharsFound.begin(), vAuxIndCharsFound.end());
			}
			firstQuoteCharWasFound = false; secondQuoteCharWasFound = false;
			firstAnotherQuoteCharWasFound = false; secondAnotherQuoteCharWasFound = false;
			indFirstQuoteChar = -1; indSecondQuoteChar = -1;
			indFirstAnotherQuoteChar = -1; indSecondAnotherQuoteChar = -1;
			firstExternQuoteChar = 0;
		}
		t_str++;
	}
}
// ----------------------------------------------------------------------------

void CAuxParse::StringMergeTokensAroundGivenCharInd(vector<pair<int, int> >& vIndTokens, vector<int>& vIndChars)
{
	int numTokens = (int)vIndTokens.size(), numChars = (int)vIndChars.size();
	if((numTokens <= 0) || (numChars <= 0)) return;

	for(int i=0; i<numChars; i++)
	{
		int curCharPos = vIndChars[i];
		int i1TokenToBeMerged = -1;
		for(int j=1; j<numTokens; j++)
		{
			pair<int, int>& prevTokenInds = vIndTokens[j - 1];
			pair<int, int>& curTokenInds = vIndTokens[j];

			if(((prevTokenInds.first + prevTokenInds.second - 1) < curCharPos) && (curCharPos < curTokenInds.first))
			{
				i1TokenToBeMerged = j - 1;
				break;
			}
		}

		if(i1TokenToBeMerged >= 0)
		{
			pair<int, int>& tokenToBeExtended = vIndTokens[i1TokenToBeMerged];
			int i1TokenToBeMerged_p_1 = i1TokenToBeMerged + 1;
			pair<int, int>& tokenToBeErased = vIndTokens[i1TokenToBeMerged_p_1];
			tokenToBeExtended.second = (tokenToBeErased.first - tokenToBeExtended.first) + tokenToBeErased.second;
			vIndTokens.erase(vIndTokens.begin() + i1TokenToBeMerged_p_1);
		}
	}
}

//-------------------------------------------------------------------------
