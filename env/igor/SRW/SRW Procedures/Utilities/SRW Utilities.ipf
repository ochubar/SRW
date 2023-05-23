#pragma rtGlobals=1		// Use modern global access method

//+++++++++++++++++++++++++++++++++++++++
//
//General SRW Utilities
//
//+++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++
//
//Returns SRW Path string 
//
//+++++++++++++++++++++++++++++++++++++++
function/T srwUtiPathStr()
string s_srwPathStr = FunctionPath("srwUtiGetValN")
s_srwPathStr = ParseFilePath(1, s_srwPathStr, ":", 1, 2)
return s_srwPathStr
end

//+++++++++++++++++++++++++++++++++++++++
//
//Function used to output various warnings  
//
//+++++++++++++++++++++++++++++++++++++++
function sRWarning(s)
string s
print s
end

//+++++++++++++++++++++++++++++++++++++++
//
//Rounding Function
//
//+++++++++++++++++++++++++++++++++++++++
function sRRound(var,num)
variable var,num

variable lo=10^round(log(var)-num)
return round(var/lo)*lo
end

//+++++++++++++++++++++++++++++++++++++++
//
//Utility proc to convert from various units of photon energies
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiPhotEnConv(val,inputunit,outputunit)
variable val=SrwVal
variable inputunit=SrwInputUnit
variable outputunit=SrwOutputUnit
prompt val,"Input Photon Energy value"
prompt inputunit,"Input Unit of the Photon Energy",popup "keV;eV;1/cm;A;nm;µm;mm;THz"
prompt outputunit,"Output Unit of the Photon Energy",popup "keV;eV;1/cm;A;nm;µm;mm;THz"
Silent 1						|	Energy Conversion  ...
PauseUpdate

SrwVal=val
SrwOutputUnit=outputunit
SrwInputUnit=inputunit

variable outputval=1
string sinp,sout

// convert  to "keV"
if (inputunit==1)
	sinp="keV";outputval=val
endif
if (inputunit==2)
	sinp="eV";outputval=val*1.e-3
endif
if (inputunit==3)
	sinp="1/cm";outputval=val*(12.39842*1e-8)
endif
if (inputunit==4)
	sinp="A";outputval=12.39842/val
endif
if (inputunit==5)
	sinp="nm";outputval=(12.39842*0.1)/val
endif
if (inputunit==6)
	sinp="µm";outputval=(12.39842*1.e-4)/val
endif
if (inputunit==7)
	sinp="mm";outputval=(12.39842*1.e-7)/val
endif
if (inputunit==8)
	//sinp="THz";outputval=val*(4.13285e-06)
	sinp="THz";outputval=val*(4.1356672e-06)
endif

SrwSmpEdep=outputval

// convert  to outputunit
if (outputunit==1)
	sout="keV";outputval=outputval
endif
if (outputunit==2)
	sout="eV";outputval=outputval*1000
endif
if (outputunit==3)
	sout="1/cm";outputval=outputval/(12.39842*1e-8)
endif
if (outputunit==4)
	sout="A";outputval=12.39842/outputval
endif
if (outputunit==5)
	sout="nm";outputval=(12.39842*0.1)/outputval
endif
if (outputunit==6)
	sout="µm";outputval=(12.39842*1.e-4)/outputval
endif
if (outputunit==7)
	sout="mm";outputval=(12.39842*1.e-7)/outputval
endif
if (outputunit==8)
	//sout="THz";outputval=outputval/(4.13285e-06)
	sout="THz";outputval=outputval/(4.1356672e-06)
endif

print  val,sinp, "  correspond(s) to : ",outputval,sout
end

//+++++++++++++++++++++++++++++++++++++++
//
//Utility proc to convert Spectral Flux units
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiSpecFluxConv(val,inputunit,outputunit)
variable val=SrwSpecFluxVal
variable inputunit=SrwSpecFluxInUnit
variable outputunit=SrwSpecFluxOutUnit
prompt val,"Input Spectral Flux value"
prompt inputunit,"Input Unit of the Spectral Flux",popup "Phot/s/.1%bw;W/eV;W/keV;W/cm^(-1)"
prompt outputunit,"Output Unit of the Spectral Flux",popup "Phot/s/.1%bw;W/eV;W/keV;W/cm^(-1)"
Silent 1						|	Spectral Flux Conversion  ...
PauseUpdate

SrwSpecFluxVal=val
SrwSpecFluxOutUnit=outputunit
SrwSpecFluxInUnit=inputunit

variable outputval=1
string sinp,sout

// convert  to "Phot/s/.1%bw"
if(inputunit==1)
	sinp="Phot/s/.1%bw"; outputval=val
endif
if(inputunit==2)
	//sinp="W/eV"; outputval = val/(1.6021892e-16)
	sinp="W/eV"; outputval = val/(1.6021765e-16)
endif
if(inputunit==3)
	//sinp="W/keV"; outputval = val/(1.6021892e-13)
	sinp="W/keV"; outputval = val/(1.6021765e-13)
endif
if(inputunit==4)
	//sinp="W/cm^(-1)"; outputval = val/((1.6021892e-16)*0.0001239854)
	sinp="W/cm^(-1)"; outputval = val/((1.6021765e-16)*0.0001239842)
endif

// convert  to outputunit
if(outputunit==1)
	sout="Phot/s/.1%bw"; //outputval=outputval
endif
if(outputunit==2)
	//sout="W/eV"; outputval *= (1.6021892e-16)
	sout="W/eV"; outputval *= (1.6021765e-16)
endif
if(outputunit==3)
	//sout="W/keV"; outputval *= (1.6021892e-13)
	sout="W/keV"; outputval *= (1.6021765e-13)
endif
if(outputunit==4)
	//sout="W/cm^(-1)"; outputval *= (1.6021892e-16)*0.0001239854
	sout="W/cm^(-1)"; outputval *= (1.6021765e-16)*0.0001239842
endif

print  val,sinp, "  correspond(s) to : ",outputval,sout
end

//+++++++++++++++++++++++++++++++++++++++
//
//Utility proc to estimate focal distance of a mirror at grazing incidence
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiOptMirFocLength(radCurv, incAngle, plane)
variable radCurv = srwUtiGetValN("radCurv", 1, "SrwUtiOptMirFocDist")
variable incAngle = srwUtiGetValN("incAngle", 45, "SrwUtiOptMirFocDist")
variable plane = srwUtiGetValN("insAngle", 3, "SrwUtiOptMirFocDist")
prompt radCurv, "Radius of Curvature [m]"
prompt incAngle, "Incidence Angle [degree]"
prompt plane, "Plane", popup "Perpendicular to Line about which Mirror is Rotated;Passing through Line about which Mirror is Rotated;Both"
Silent 1						|	Spectral Flux Conversion  ...
PauseUpdate

srwUtiSetValN("radCurv", radCurv, "SrwUtiOptMirFocDist")
srwUtiSetValN("incAngle", incAngle, "SrwUtiOptMirFocDist")
srwUtiSetValN("plane", plane, "SrwUtiOptMirFocDist")

variable focLen1, focLen2 = 0
if(plane == 1)
	focLen1 = srwUtiOptMirFocLen(radCurv, incAngle, plane)
	print "Focal Length in Plane Perpendicular to Line about which Mirror is Rotated:", focLen1
else
	if(plane == 2)
		focLen2 = srwUtiOptMirFocLen(radCurv, incAngle, 2)
		print "Focal Length in Plane Passing through Line about which Mirror is Rotated:", focLen2
	else
		focLen1 = srwUtiOptMirFocLen(radCurv, incAngle, 1)
		focLen2 = srwUtiOptMirFocLen(radCurv, incAngle, 2)
		print "Focal Length in Plane Perpendicular to Line about which Mirror is Rotated:", focLen1
		print "Focal Length in Plane Passing through Line about which Mirror is Rotated:", focLen2	
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Function to estimate focal distance of a mirror at grazing incidence
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiOptMirFocLen(radCurv, incAngle, plane)
variable radCurv, incAngle, plane
if((plane > 2) %| (plane <= 0))
	abort "Plane is not defined correctly"
endif

variable halfR = 0.5*radCurv
//variable cosFact = cos((0.5 - incAngle/180)*Pi) 
variable cosFact = cos(Pi*incAngle/180) //OC, ML141108

if(plane == 1)
	return halfR*cosFact
else
	return halfR/cosFact
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Returns Item Number from a Wave List with a given separator
// 1- based
//+++++++++++++++++++++++++++++++++++++++
Function sRWaveListItemNo(AllNames, SepChar, Name)
String AllNames, Name, SepChar;

Variable AllNamesLen = strlen(AllNames);

//if(AllNamesLen==0) 
//	return -1;
//endif

Variable i = 0, k = 0, NameFits = 1, ItemNo = 1;
Variable SemiCol = char2num(SepChar);
	//Print AllNames
do
	Variable chAllNames = char2num(AllNames[i]), chName = char2num(Name[k]);

	if(chAllNames == SemiCol)
		if(NameFits == 1)
			return ItemNo;
		else
			ItemNo += 1;
			k = -1; NameFits = 1;
			if(i == (AllNamesLen - 1))
				NameFits = 0;
			endif
		endif
	else
		if(chAllNames != chName)
			NameFits = 0;
		endif
	endif
	i += 1; k += 1;
while(i < AllNamesLen)

ItemNo -= 1;

if(NameFits == 1)
	return ItemNo;
else 
	return 0; // If no fit
endif

End

//+++++++++++++++++++++++++++++++++++++++
//
//Fills Item String from a Wave List with a given separator
// 1- based
//+++++++++++++++++++++++++++++++++++++++
Function/S sRWaveListItemName(AllNames, SepChar, InItemNo)
String AllNames, SepChar;
Variable InItemNo;

String ItemName;

Variable AllNamesLen = strlen(AllNames);
Variable i = 0, ItemNo = 1;

Variable SemiCol = char2num(SepChar);
ItemName = "";
	//Print AllNames
do
	Variable chAllNames = char2num(AllNames[i]);

	if(chAllNames == SemiCol)
		if(ItemNo == InItemNo)
				//Print ItemName
			return ItemName;
		else
			ItemNo += 1;
			ItemName = "";
		endif
	else
		ItemName += AllNames[i];
	endif
	i += 1;
while(i < AllNamesLen)

return ""; // No fit
End

//+++++++++++++++++++++++++++++++++++++++
//
//Counts how many times a character encounters in a given string
//
//+++++++++++++++++++++++++++++++++++++++
Function sRCountCharCasesInStr(Str, Char)
String Str, Char;

Variable StrLength = strlen(Str);
Variable i = 0, Count = 0;
Variable vChar = char2num(Char);
do
	Variable vCurrChar = char2num(Str[i]);
	if(vCurrChar == vChar)
		Count += 1;
	endif
	i += 1;
while(i < StrLength)
return Count;
End

//+++++++++++++++++++++++++++++++++++++++
//
//Kills all waves of a given type
//
//+++++++++++++++++++++++++++++++++++++++
Function sRKillAllWavesOfType(Type)
String Type;

String WavesToKill = Wavelist("*"+Type,";","");
Variable AmOfWaves = sRCountCharCasesInStr(WavesToKill, ";");

String CurrWaveName;
Variable i = 0;
do
	CurrWaveName = sRWaveListItemName(WavesToKill, ";", i+1);
	 KillWaves/Z $CurrWaveName;
	 	//Print CurrWaveName;

	i += 1;
while(i < AmOfWaves)

End

//+++++++++++++++++++++++++++++++++++++++
//
//Data Windowing function
//
//+++++++++++++++++++++++++++++++++++++++
Function srwUtiDataWindowZero(s, sdep, sfin)
Variable s, sdep, sfin;
if((s > sdep) %& (s < sfin))
	return 1;
else
	return 0;
endif

End

//+++++++++++++++++++++++++++++++++++++++
//
//Show SRW Help Topic from non-formatted Help File notebook.
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiShowHelpTopic(TopicStr)
string TopicStr

PathInfo Igor
newpath/O/Q SRW_Proc, S_path+"SRW:SRW Help"
opennotebook/R/P=SRW_Proc/N=srwHelp "SRW Help.ifn"

notebook srwHelp selection={startOfFile, startOfFile}
notebook srwHelp findText={TopicStr, 1};
end

//+++++++++++++++++++++++++++++++++++++++
//
// Utility proc to trigger printing extra info to history
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiTriggerPrint(OnOrOff)
variable OnOrOff=SrwAllowPrintingExtraInfo
prompt OnOrOff,SrwPAllowPrintingExtraInfo,popup SrwPOPUPAllowPrintingExtraInfo
Silent 1						|	...
PauseUpdate

SrwAllowPrintingExtraInfo=OnOrOff

end

//+++++++++++++++++++++++++++++++++++++++
//
// Utility to steer Interruption Time to allow multitasking on Mac
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiInterTime(Delta)
Variable Delta
srUtiInterTime(Delta)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Print SRW Warning message
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiPrintWarn(WarnStr)
String WarnStr
print SrwPWarnMessageHeader
print SrwPWarnFirstSymb, WarnStr
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiKillAllGraphs()
DoWindow/K Graph0
DoWindow/K Graph1
DoWindow/K Graph2
DoWindow/K Graph3
DoWindow/K Graph4
DoWindow/K Graph5
DoWindow/K Graph6
DoWindow/K Graph7
DoWindow/K Graph8
DoWindow/K Graph9
DoWindow/K Graph10
DoWindow/K Graph11
DoWindow/K Graph12
DoWindow/K Graph13
DoWindow/K Graph14
DoWindow/K Graph15
DoWindow/K Graph16
DoWindow/K Graph17
DoWindow/K Graph18
DoWindow/K Graph19
DoWindow/K Graph20
DoWindow/K Graph21
DoWindow/K Graph22
DoWindow/K Graph23
DoWindow/K Graph24
DoWindow/K Graph25
DoWindow/K Graph26
DoWindow/K Graph27
DoWindow/K Graph28
DoWindow/K Graph29
DoWindow/K Graph30
DoWindow/K Graph31
DoWindow/K Graph32
DoWindow/K Graph33
DoWindow/K Graph34
DoWindow/K Graph35
DoWindow/K Graph36
DoWindow/K Graph37
DoWindow/K Graph38
DoWindow/K Graph39
DoWindow/K Graph40
DoWindow/K Graph41
DoWindow/K Graph42
DoWindow/K Graph43
DoWindow/K Graph44
DoWindow/K Graph45
DoWindow/K Graph46
DoWindow/K Graph47
DoWindow/K Graph48
DoWindow/K Graph49
DoWindow/K Graph50
DoWindow/K Graph51
DoWindow/K Graph52
DoWindow/K Graph53
DoWindow/K Graph54
DoWindow/K Graph55
DoWindow/K Graph56
DoWindow/K Graph57
DoWindow/K Graph58
DoWindow/K Graph59
DoWindow/K Graph60
DoWindow/K Graph61
DoWindow/K Graph62
DoWindow/K Graph63
DoWindow/K Graph64
DoWindow/K Graph65
DoWindow/K Graph66
DoWindow/K Graph67
DoWindow/K Graph68
DoWindow/K Graph69
DoWindow/K Graph70
DoWindow/K Graph71
DoWindow/K Graph72
DoWindow/K Graph73
DoWindow/K Graph74
DoWindow/K Graph75
DoWindow/K Graph76
DoWindow/K Graph77
DoWindow/K Graph78
DoWindow/K Graph79
DoWindow/K Graph80
DoWindow/K Graph81
DoWindow/K Graph82
DoWindow/K Graph83
DoWindow/K Graph84
DoWindow/K Graph85
DoWindow/K Graph86
DoWindow/K Graph87
DoWindow/K Graph88
DoWindow/K Graph89
DoWindow/K Graph90
DoWindow/K Graph91
DoWindow/K Graph92
DoWindow/K Graph93
DoWindow/K Graph94
DoWindow/K Graph95
DoWindow/K Graph96
DoWindow/K Graph97
DoWindow/K Graph98
DoWindow/K Graph99
DoWindow/K Graph100
// Make more smart
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiElecLib()

SrwUtiTriggerPrint(2)

//SrwElecFilament("SPRING8_LowBeta",8.,0.1,0,0,0,0,0);SrwElecThick("SPRING8_LowBeta"+SrwElecType,0.0012,7.3,0.07,0.95,5.5,0,0,0,0)
//SrwElecFilament("SPRING8_HighBeta",8.,0.1,0,0,0,0,0);SrwElecThick("SPRING8_HighBeta"+SrwElecType,0.0012,7.3,0.07,24.,11.9,0,0,0,0)
SrwElecFilament("SPRING8_HighBeta",8,0.1,0,0,0,0,0);SrwElecThick("SPRING8_HighBeta"+SrwElecType,0.0011,3.4,0.0068,22.6,5.6,0,0,0.107,0)
//SrwElecFilament("APS",7.,0.1,0,0,0,0,0);SrwElecThick("APS"+SrwElecType,0.00096,7.7,0.12,16.,3.,0,0,0.13,0)
SrwElecFilament("APS",7,0.1,0,0,0,0,0);SrwElecThick("APS"+SrwElecType,0.00096,2.79,0.0084,22.7,3.1,0,0,0.206,0)
//SrwElecFilament("ESRF_LowBeta",6.04,0.2,0,0,0,0,0);SrwElecThick("ESRF_LowBeta"+SrwElecType,0.001,4,0.03,0.5,2.73,0,0,0,0)
SrwElecFilament("ESRF_LowBeta",6.03,0.2,0,0,0,0,0);SrwElecThick("ESRF_LowBeta"+SrwElecType,0.001,4,0.005,0.35,2.97,0,0,0.031,0) //Dec.2010
//SrwElecFilament("ESRF_HighBeta",6.04,0.2,0,0,0,0,0);SrwElecThick("ESRF_HighBeta"+SrwElecType,0.001,4,0.03,35.6,2.5,0,0,0,0)
SrwElecFilament("ESRF_HighBeta",6.03,0.2,0,0,0,0,0);SrwElecThick("ESRF_HighBeta"+SrwElecType,0.001,4,0.005,37.6,2.95,0,0,0.134,0)
SrwElecFilament("ELETTRA",2.,0.3,0,0,0,0,0);SrwElecThick("ELETTRA"+SrwElecType,0.0028,7.,0.07,8.2,2.6,0,0,0,0)
SrwElecFilament("MAXII",1.5,0.2,0,0,0,0,0);SrwElecThick("MAXII_ebm",0.0007,8.8,0.09,13,2.3,0,0,0.13,0)
//SrwElecFilament("SuperACO",0.8,0.4,0,0,0,0,0);SrwElecThick("SuperACO"+SrwElecType,0.0006,37,3.7,6,11,0,0,0,0)
SrwElecFilament("SOLEIL_LongSect",2.75,0.5,0,0,0,0,0);SrwElecThick("SOLEIL_LongSect"+SrwElecType,0.001016,3.73,0.037,10.09,8.01,0,0,0.2,0)
SrwElecFilament("SOLEIL_MedSect",2.75,0.5,0,0,0,0,0);SrwElecThick("SOLEIL_MedSect"+SrwElecType,0.001016,3.73,0.037,4,1.77,0,0,0.13,0)
SrwElecFilament("SOLEIL_ShortSect",2.75,0.5,0,0,0,0,0);SrwElecThick("SOLEIL_ShortSect"+SrwElecType,0.001016,3.73,0.037,17.78,1.75,0,0,0.28,0)
SrwElecFilament("SOLEIL_BM1DEG",2.75,0.5,0,0,0,0,0);SrwElecThick("SOLEIL_BM1DEG"+SrwElecType,0.001016,3.73,0.037,0.603,16.53,0.776,0.931,0.039,-0.088)
SrwElecFilament("SOLEIL_BM4DEG",2.75,0.5,0,0,0,0,0);SrwElecThick("SOLEIL_BM4DEG"+SrwElecType,0.001016,3.73,0.037,0.375,16.01,0.024,0.899,0.021,-0.037)
SrwElecFilament("DIAMOND_LowBeta",3,0.3,0,0,0,0,0);SrwElecThick("DIAMOND_LowBeta"+SrwElecType,0.001,2.7965,0.0243202,2.28286,2.51447,0,0,0,0)
//SrwElecFilament("NSLSII_LowBeta_Day1",3,0.5,0,0,0,0,0);SrwElecThick("NSLSII_LowBeta_Day1"+SrwElecType,0.00089,0.9,0.008,2.02,1.06,0,0,0,0)
SrwElecFilament("NSLSII_LowBeta_Day1",3,0.5,0,0,0,0,0);SrwElecThick("NSLSII_LowBeta_Day1"+SrwElecType,0.00089,0.9,0.008,1.84,1.17,0,0,0,0)

SrwElecFilament("NSLSII_HighBeta_Day1",3,0.5,0,0,0,0,0);SrwElecThick("NSLSII_HighBeta_Day1"+SrwElecType,0.00089,0.9,0.008,20.85,3.4,0,0,0,0)
SrwElecFilament("NSLSII_BM_Day1",3,0.5,0,0,0,0,0);SrwElecThick("NSLSII_BM_Day1"+SrwElecType,0.00089,0.9,0.008,1.5,22.5,0,-0.9,0.137,-0.1)
SrwElecFilament("NSLSII_TPW_Day1",3,0.5,0,0,0,0,0);SrwElecThick("NSLSII_TPW_Day1"+SrwElecType,0.00089,0.9,0.008,2.956,19.653,1.932,-0.806,0.137,-0.105)

SrwUtiTriggerPrint(1)

//SrwUtiElecLibPrintMes("SPRING8_LowBeta")
SrwUtiElecLibPrintMes("SPRING8_HighBeta")
SrwUtiElecLibPrintMes("APS")
SrwUtiElecLibPrintMes("ESRF_LowBeta")
SrwUtiElecLibPrintMes("ESRF_HighBeta")
SrwUtiElecLibPrintMes("ELETTRA")
SrwUtiElecLibPrintMes("MAXII")
//SrwUtiElecLibPrintMes("SuperACO")
SrwUtiElecLibPrintMes("SOLEIL_LongSect")
SrwUtiElecLibPrintMes("SOLEIL_MedSect")
SrwUtiElecLibPrintMes("SOLEIL_ShortSect")
SrwUtiElecLibPrintMes("SOLEIL_BM1DEG")
SrwUtiElecLibPrintMes("SOLEIL_BM4DEG")
SrwUtiElecLibPrintMes("DIAMOND_LowBeta")
SrwUtiElecLibPrintMes("NSLSII_LowBeta_Day1")
SrwUtiElecLibPrintMes("NSLSII_HighBeta_Day1")
SrwUtiElecLibPrintMes("NSLSII_BM_Day1")
SrwUtiElecLibPrintMes("NSLSII_TPW_Day1")

end

//+++++++++++++++++++++++++++++++++++++++
//Runs all SRW Examples one-by-one
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiAllExam()
SrwUtiTriggerPrint(2)

variable DelT = 10 //[s] Time to keep graphs after execution of each example

SrwExamUR(); SrwUtiKillExamWindows(DelT)
SrwExamBMSRPolar(); SrwUtiKillExamWindows(DelT)
SrwExamER(1); SrwUtiKillExamWindows(DelT)
SrwExamERPowerDens(); SrwUtiKillExamWindows(DelT)
//SrwExamCSR_RotBunch(); SrwUtiKillExamWindows(DelT)

SrwExamURSpecThickEbeam(); SrwUtiKillExamWindows(DelT)
SrwExamURIntDistrThickEbeam(); SrwUtiKillExamWindows(DelT)
SrwExamURPowerDens(); SrwUtiKillExamWindows(DelT)
SrwExamURElliptic(); SrwUtiKillExamWindows(DelT)
SrwExamBrilUR(); SrwUtiKillExamWindows(DelT)

SrwExamWigPlanar(); SrwUtiKillExamWindows(DelT)
SrwExamWigElliptic(); SrwUtiKillExamWindows(DelT)

SrwExamStdBM(); SrwUtiKillExamWindows(DelT)

SrwExamSASE_GENESIS(); SrwUtiKillExamWindows(DelT)

SrwExamImagBm(); SrwUtiKillExamWindows(DelT)
SrwExamImagUnd(); SrwUtiKillExamWindows(DelT)
SrwExamUndRad2Slits(); SrwUtiKillExamWindows(DelT)
SrwExamDiffrER(); SrwUtiKillExamWindows(DelT)
SrwExamXrayLensCirc(); SrwUtiKillExamWindows(DelT)
SrwExamGsnBm(); SrwUtiKillExamWindows(DelT)
//add next here

SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiElecLibPrintMes(BeamName)
String BeamName
if(SrwAllowPrintingExtraInfo==1)
	String StrToPrint="\""+BeamName+"\""+" electron beam is loaded"
	print StrToPrint
endif
end

//+++++++++++++++++++++++++++++++++++++++
function srwUtiCurrentPath()
//???
end

//+++++++++++++++++++++++++++++++++++++++
//Returns index of field in a one-dimensional text wave 
//where the txt string occours for the last time
//+++++++++++++++++++++++++++++++++++++++
function srwUtiLastIndOfTxtFld(TxtW, txt)
wave/T TxtW 
String txt

Variable SizeW = DimSize(TxtW, 0)
String TestStr
Variable i = 0, OutInd = -1
do
	TestStr = TxtW[i]
	
	if(strsearch(TestStr,txt,0) >= 0)
		OutInd = i
	endif
	
	i += 1
while(i < SizeW)
return OutInd
end

//+++++++++++++++++++++++++++++++++++++++
//Composes new wave name
//if(EnsureUnique != 0) and NameMain is not unique,
//digit suffix is appended to the NameMain
//+++++++++++++++++++++++++++++++++++++++
function/S SrwUtiGiveNewName(NameMain, NameEnding, EnsureUnique)
String NameMain, NameEnding
Variable EnsureUnique

String TestName = NameMain + NameEnding
if(EnsureUnique == 0)
	return TestName
endif

if(exists(TestName) == 0)
	return TestName
endif

String/G SrwSeparator
Variable Ind = 1
do
	TestName = NameMain + SrwSeparator + num2str(Ind) + NameEnding
	Ind += 1
while(exists(TestName) != 0)
return TestName
end

//+++++++++++++++++++++++++++++++++++++++
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiKillExamWindows(DelayTime)
Variable DelayTime = 4 // Time delay before killing example graphs
DoUpdate; Sleep/S 0.5*DelayTime; DoWindow/K srwHelp; Sleep/S 0.5*DelayTime; SrwUtiKillAllGraphs()
end

//+++++++++++++++++++++++++++++++++++++++
//Tries to find current value of a global numerical variable
//(is used to facilitate storage/extraction of previous values of arguments in dialogs)
//+++++++++++++++++++++++++++++++++++++++
function srwUtiGetValN(VarName, DefaultVal, FuncName)
string VarName, FuncName
variable DefaultVal
string TotName = srwUtiEnsureShortName(srwUtiComposeName(VarName, "_", FuncName))
//NVAR var1=$TotName
//if(NVAR_Exists(var1)==0)
if(exists(TotName)==0)
	return DefaultVal
endif
NVAR var1=$TotName
return var1
end

//+++++++++++++++++++++++++++++++++++++++
///Tries to find current value of a global string variable
//(is used to facilitate storage/extraction of previous values of arguments in dialogs)
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiGetValS(VarName, DefaultVal, FuncName)
string VarName, DefaultVal, FuncName
string TotName = srwUtiEnsureShortName(srwUtiComposeName(VarName, "_", FuncName))
//SVAR var1=$TotName
//if(SVAR_Exists(var1)==0)
if(exists(TotName)==0)
	return DefaultVal
endif
SVAR var1=$TotName
return var1
end

//+++++++++++++++++++++++++++++++++++++++
//Sets a global numerical variable to a given value
//if the variable does not exist, it is created
//(is used to facilitate storage/extraction of previous values of arguments in dialogs)
//+++++++++++++++++++++++++++++++++++++++
function srwUtiSetValN(VarName, Val, FuncName)
variable Val
string VarName, FuncName
string TotName = srwUtiEnsureShortName(srwUtiComposeName(VarName, "_", FuncName))
//NVAR var1=$TotName
//if(NVAR_Exists(var1)==0)
if(exists(TotName)==0)
	variable/G $TotName = Val
else
	NVAR var1=$TotName
	var1 = Val
endif
return Val
end

//+++++++++++++++++++++++++++++++++++++++
//Sets a global string variable to a given value
//if the variable does not exist, it is created
//(is used to facilitate storage/extraction of previous values of arguments in dialogs)
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiSetValS(VarName, Val, FuncName)
string Val, VarName, FuncName
string TotName = srwUtiEnsureShortName(srwUtiComposeName(VarName, "_", FuncName))
//SVAR var1=$TotName
//if(SVAR_Exists(var1)==0)
if(exists(TotName)==0)
	string/G $TotName = Val
else
	SVAR var1=$TotName
	var1 = Val
endif
return Val
end

//+++++++++++++++++++++++++++++++++++++++
//Composes a string name from two parts separated by 
//a chain of symbols (if any)
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiComposeName(Part1, SepSymb, Part2)
string Part1, SepSymb, Part2
if(strlen(Part2) == 0)
	return Part1
endif
return (Part1 + SepSymb + Part2)
end

//+++++++++++++++++++++++++++++++++++++++
//Trancates a string to 31 symbol
//(Igor does not support names longer than 31 symbol)
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiEnsureShortName(VarName)
string VarName
variable MaxStrLen = 31

if(strlen(VarName) <= MaxStrLen)
	return VarName
endif
variable i=0
string OutName = ""
do
	OutName[i] = VarName[i]
	i += 1
while(i < MaxStrLen)
return OutName
end

//+++++++++++++++++++++++++++++++++++++++
//Trancates a string to MaxStrLen symbols
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiTruncString(VarName, MaxStrLen)
string VarName
variable MaxStrLen

if(strlen(VarName) <= MaxStrLen)
	return VarName
endif
variable i=0
string OutName = ""
do
	OutName[i] = VarName[i]
	i += 1
while(i < MaxStrLen)
return OutName
end

//+++++++++++++++++++++++++++++++++++++++
//Returns the end part of SRW wave name, e.g. "rad", of "x"
//assumes SeparSymb has only one symbol
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiGetNameEnd(sName, SeparSymb)
string sName, SeparSymb

string OutName = ""
variable sNameLen = strlen(sName)
if((sNameLen <= 0) %| (strlen(SeparSymb) <= 0))
	return OutName
endif

string sTestSymb = ""
variable i = sNameLen - 1
do
	sTestSymb = sName[i]
	if(cmpstr(sTestSymb, SeparSymb) == 0)
		break
	endif
	i -= 1
while(i >= 0)
if(i <= 0)
	return OutName
endif
variable i0 = i
if(i0 >= sNameLen)
	return OutName
endif
do
	OutName[i - i0] = sName[i]
	i += 1
while(i < sNameLen)

return OutName
end

//+++++++++++++++++++++++++++++++++++++++
//Step function
//+++++++++++++++++++++++++++++++++++++++
function srwUtiStep(p)
variable p
if(p < 0.)
	return 0.
else
	return 1.
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Gate function
//+++++++++++++++++++++++++++++++++++++++
function srwUtiNonZeroInterval(p, pmin, pmax)
variable p, pmin, pmax
if((p < pmin) %| (p >= pmax))
	return 0.
else
	return 1.
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Gate function including two boundaries
//+++++++++++++++++++++++++++++++++++++++
function srwUtiNonZeroIntervB(p, pmin, pmax)
variable p, pmin, pmax
if((p < pmin) %| (p > pmax))
	return 0.
else
	return 1.
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Trivial auxiliary function for wave operations
//+++++++++++++++++++++++++++++++++++++++
function srwUtiPow2(x)
variable x
return x*x
end

//+++++++++++++++++++++++++++++++++++++++
//Function used at computation of power density distribution 
//after an aperture and a drift space.
//+++++++++++++++++++++++++++++++++++++++
function srwUtiPowDensApertFunc(x, L, R0, R1, dx, x0)
variable x, L, R0, R1, dx, x0

variable R1mR0 = R1 - R0
variable HalfL = 0.5*L
variable FactB = 1. + R1mR0/(R0 - HalfL)
variable FactS = 1. + R1mR0/(R0 + HalfL)

variable x1 = x0 - 0.5*dx
variable x2 = x1 + dx

variable x1m = x1*FactB
variable x1p = x1*FactS
if(x1 >= 0.)
	x1m = x1*FactS
	x1p = x1*FactB
endif
variable x2m = x2*FactS
variable x2p = x2*FactB
if(x2 < 0.)
	x2m = x2*FactB
	x2p = x2*FactS
endif

if((x <= x1m) %| (x >= x2p))
	return 0.
endif
if((x >= x1p) %& (x <= x2m))
	return 1.
endif

variable dL = 0.
if((x > x1m) %& (x < x1p))
	if(x1 == 0.)
		return 1.
	endif
	if(x1 < 0)
		return (HalfL - (R0 - R1mR0/(x/x1 - 1.)))/L
	else
		return ((R0 - R1mR0/(x/x1 - 1.)) + HalfL)/L
	endif
endif
if((x > x2m) %& (x < x2p))
	if(x2 == 0.)
		return 1.
	endif
	if(x2 > 0)
		return (HalfL - (R0 - R1mR0/(x/x2 - 1.)))/L
	else
		return ((R0 - R1mR0/(x/x2 - 1.)) + HalfL)/L
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
function srwUtiNumOrMax(q, qmax)
variable q, qmax
if(q < qmax)
	return q
endif
return qmax
end

//+++++++++++++++++++++++++++++++++++++++
function srwUtiMaxOfTwo(v1, v2)
variable v1, v2
if(v1 > v2)
	return v1
endif
return v2
end

//+++++++++++++++++++++++++++++++++++++++
function srwUtiNumOrDefault(q, qmin, qmax, qdef)
variable q, qmin, qmax, qdef
if((q < qmin) %| (q > qmax))
	return qdef
endif
return q
end

//+++++++++++++++++++++++++++++++++++++++
function SrwUtiSinXdX(x) //Auxiliary sin(x)/x
variable x
if(x == 0)
	return 1
else
	return sin(x)/x
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Transpose 2D wave
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiTransposeWave2D(nmWave, dupl, nmWaveD)
string nmWave
variable dupl
string nmWaveD

//dddddddddddddddddddddddddddddddddddddddddd

end

//+++++++++++++++++++++++++++++++++++++++
//Makes in-place convolution of planes of a 3D wave with 
//2D Interval or 2D Gaussian distribution
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiConvWave3DWithFunc2D(NameW, DelX, DelY, convType, indPlanes)
string NameW = srwUtiGetValS("NameW", "", "SrwUtiConvWave3DWithFunc2D")
variable DelX = srwUtiGetValN("DelX", 1., "SrwUtiConvWave3DWithFunc2D")
variable DelY = srwUtiGetValN("DelY", 1., "SrwUtiConvWave3DWithFunc2D")
variable convType = srwUtiGetValN("convType", 1, "SrwUtiConvWave3DWithFunc2D")
variable indPlanes = srwUtiGetValN("indPlanes", 3, "SrwUtiConvWave3DWithFunc2D")
prompt NameW, "2D wave to convolve with Interval (in place)", popup WaveList("*",";","TEXT:0,DIMS:3")
prompt DelX, "Size of 2D Function in 1st Dimension"
prompt DelY, "Size of 2D Function in 2nd Dimension"
prompt convType, "Function Type to Convolve with", popup "Gaussian;Interval"
prompt indPlanes, "Index of 3D Wave Planes", popup "0;1;2"
Silent 1						|	Making 2D convolution ...
PauseUpdate

srwUtiSetValS("NameW", NameW, "SrwUtiConvWave3DWithFunc2D")
srwUtiSetValN("DelX", DelX, "SrwUtiConvWave3DWithFunc2D")
srwUtiSetValN("DelY", DelY, "SrwUtiConvWave3DWithFunc2D")
srwUtiSetValN("convType", convType, "SrwUtiConvWave3DWithFunc2D")
srwUtiSetValN("indPlanes", indPlanes, "SrwUtiConvWave3DWithFunc2D")

variable indPlanes0based = indPlanes - 1
variable indX = 0, indY = 1
if(indPlanes0based == 0)
	indX = 1; indY = 2
endif
if(indPlanes0based == 1)
	indX = 0; indY = 2
endif

variable dimSizeX = dimsize($NameW, indX)
variable startX = dimoffset($NameW, indX)
variable stepX = dimdelta($NameW, indX)
variable dimSizeY = dimsize($NameW, indY)
variable startY = dimoffset($NameW, indY)
variable stepY = dimdelta($NameW, indY)

string nmPlane = "wPlaneConvWave3DWithFunc2D"
make/O/N=(dimSizeX, dimSizeY) $nmPlane
SetScale/P x startX, stepX, WaveUnits($NameW, indX), $nmPlane
SetScale/P y startY, stepY, WaveUnits($NameW, indY), $nmPlane

variable numPlanes = dimsize($NameW, indPlanes0based)
variable iPlane = 0
do
	if(indPlanes0based == 0)
		$nmPlane = $NameW[iPlane][p][q]
		if(convType == 1)
			SrwUtiConvWaveWithGaus2D(nmPlane, DelX, DelY)
		endif
		if(convType == 2)
			SrwUtiConvWaveWithInterv2D(nmPlane, DelX, DelY)
		endif
		$NameW[iPlane][][] = $nmPlane[q][r]
	endif
	if(indPlanes0based == 1)
		$nmPlane = $NameW[p][iPlane][q]
		if(convType == 1)
			SrwUtiConvWaveWithGaus2D(nmPlane, DelX, DelY)
		endif
		if(convType == 2)
			SrwUtiConvWaveWithInterv2D(nmPlane, DelX, DelY)
		endif
		$NameW[][iPlane][] = $nmPlane[p][r]
	endif
	if(indPlanes0based == 2)
		$nmPlane = $NameW[p][q][iPlane]
		if(convType == 1)
			SrwUtiConvWaveWithGaus2D(nmPlane, DelX, DelY)
		endif
		if(convType == 2)
			SrwUtiConvWaveWithInterv2D(nmPlane, DelX, DelY)
		endif
		$NameW[][][iPlane] = $nmPlane[p][q]
	endif
	iPlane += 1
while(iPlane < numPlanes)
killwaves/Z $nmPlane
end

//+++++++++++++++++++++++++++++++++++++++
//Makes in-place convolution of a 2D wave with 2D Interval 
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiConvWaveWithInterv2D(NameW, DelX, DelY)
string NameW = srwUtiGetValS("NameW", "", "SrwUtiConvWaveWithInterv2D")
variable DelX = srwUtiGetValN("DelX", 1., "SrwUtiConvWaveWithInterv2D")
variable DelY = srwUtiGetValN("DelY", 1., "SrwUtiConvWaveWithInterv2D")
prompt NameW, "2D wave to convolve with Interval (in place)", popup WaveList("*",";","TEXT:0,DIMS:2")
prompt DelX, "Horizontal Size of 2D Interval"
prompt DelY, "Vertical Size of 2D Interval"
Silent 1						|	Making 2D convolution ...
PauseUpdate
srwUtiSetValS("NameW", NameW, "SrwUtiConvWaveWithInterv2D")
srwUtiSetValN("DelX", DelX, "SrwUtiConvWaveWithInterv2D")
srwUtiSetValN("DelY", DelY, "SrwUtiConvWaveWithInterv2D")

string AuxWaveC = "AuxUtiConvWaveWithInterv2D"
variable NxOrig = DimSize($NameW, 0)
variable NyOrig = DimSize($NameW, 1)
variable NxOrigM1 = NxOrig - 1
variable NyOrigM1 = NyOrig - 1
variable NxAct = 2*trunc(0.5*NxOrig + 0.000001)
variable NyAct = 2*trunc(0.5*NyOrig + 0.000001)
variable StepX = DimDelta($NameW, 0)
variable StepY = DimDelta($NameW, 1)

make/O/C/N=(NxAct, NyAct) $AuxWaveC
$AuxWaveC = $NameW[srwUtiNumOrMax(p, NxOrigM1)][srwUtiNumOrMax(q, NyOrigM1)]
SetScale/P x -0.5*StepX*NxAct, StepX, WaveUnits($NameW, 0), $AuxWaveC
SetScale/P y -0.5*StepY*NyAct, StepY, WaveUnits($NameW, 1), $AuxWaveC

srFFT2D($AuxWaveC, 1)

variable piAx = Pi*DelX, piAy = Pi*DelY
$AuxWaveC *= SrwUtiSinXdX(piAx*x)*SrwUtiSinXdX(piAy*y)

srFFT2D($AuxWaveC, -1)
$NameW = real($AuxWaveC[p][q])

KillWaves/Z $AuxWaveC
end

//+++++++++++++++++++++++++++++++++++++++
//Makes in-place convolution of a 2D wave with 2D Gaussian 
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiConvWaveWithGaus2D(NameW, SigX, SigY)
string NameW = srwUtiGetValS("NameW", "", "SrwUtiConvWaveWithGaus2D")
variable SigX = srwUtiGetValN("SigX", 1., "SrwUtiConvWaveWithGaus2D")
variable SigY = srwUtiGetValN("SigY", 1., "SrwUtiConvWaveWithGaus2D")
prompt NameW, "2D wave to convolve with Gaussian (in place)", popup WaveList("*",";","TEXT:0,DIMS:2")
prompt SigX, "Horizontal RMS size of 2D Gaussian"
prompt SigY, "Vertical RMS size of 2D Gaussian"
Silent 1						|	Making 2D convolution ...
PauseUpdate
srwUtiSetValS("NameW", NameW, "SrwUtiConvWaveWithGaus2D")
srwUtiSetValN("SigX", SigX, "SrwUtiConvWaveWithGaus2D")
srwUtiSetValN("SigY", SigY, "SrwUtiConvWaveWithGaus2D")

string AuxWaveC = "AuxUtiConvWaveWithGaus2D"
variable NxOrig = DimSize($NameW, 0)
variable NyOrig = DimSize($NameW, 1)
variable NxOrigM1 = NxOrig - 1
variable NyOrigM1 = NyOrig - 1
variable NxAct = 2*trunc(0.5*NxOrig + 0.000001)
variable NyAct = 2*trunc(0.5*NyOrig + 0.000001)
variable StepX = DimDelta($NameW, 0)
variable StepY = DimDelta($NameW, 1)

make/O/C/N=(NxAct, NyAct) $AuxWaveC
$AuxWaveC = $NameW[srwUtiNumOrMax(p, NxOrigM1)][srwUtiNumOrMax(q, NyOrigM1)]
SetScale/P x -0.5*StepX*NxAct, StepX, WaveUnits($NameW, 0), $AuxWaveC
SetScale/P y -0.5*StepY*NyAct, StepY, WaveUnits($NameW, 1), $AuxWaveC

srFFT2D($AuxWaveC, 1)
variable c0 = 2*Pi*Pi //OC fixed 16082004 //4*Pi*Pi
variable cx = c0*SigX*SigX
variable cy = c0*SigY*SigY
$AuxWaveC *= exp(-cx*x*x - cy*y*y)

srFFT2D($AuxWaveC, -1)
$NameW = real($AuxWaveC[p][q])

KillWaves/Z $AuxWaveC
end

//+++++++++++++++++++++++++++++++++++++++
//Makes in-place convolution of a 1D wave with 1D Gaussian 
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiConvWaveWithGaus1D(NameW, SigX)
string NameW = srwUtiGetValS("NameW", "", "SrwUtiConvWaveWithGaus1D")
variable SigX = srwUtiGetValN("SigX", 1., "SrwUtiConvWaveWithGaus1D")
prompt NameW, "Name of the 1D wave to convolve with Gaussian"
prompt SigX, "RMS size of 1D Gaussian"
Silent 1						|	Making 1D convolution ...
PauseUpdate
srwUtiSetValS("NameW", NameW, "SrwUtiConvWaveWithGaus1D")
srwUtiSetValN("SigX", SigX, "SrwUtiConvWaveWithGaus1D")

string AuxWaveC = "AuxUtiConvWaveWithGaus1D"
variable NxOrig = DimSize($NameW, 0)
variable NxOrigM1 = NxOrig - 1
variable NxAct = 2*trunc(0.5*NxOrig + 0.000001)
variable StepX = DimDelta($NameW, 0)

make/O/C/N=(NxAct) $AuxWaveC
$AuxWaveC = $NameW[srwUtiNumOrMax(p, NxOrigM1)]
SetScale/P x -0.5*StepX*NxAct, StepX, WaveUnits($NameW, 0), $AuxWaveC

srFFT1D($AuxWaveC, 1)

variable c0 = 2*Pi*Pi //OC fixed 16082004 //4*Pi*Pi
variable cx = c0*SigX*SigX
$AuxWaveC *= exp(-cx*x*x)

srFFT1D($AuxWaveC, -1)
$NameW = real($AuxWaveC[p])

KillWaves/Z $AuxWaveC
end

//+++++++++++++++++++++++++++++++++++++++
//Makes in-place convolution of a 1D wave with 1D Gaussian,
//with linearly changing RMS.
//Implemented for treating Energy Spread
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiConvWaveWithGausLinVar1D(wIn, SigX_p, numSigInteg, wOut)
wave wIn //"1D wave to convolve with Gaussian"
wave wOut //"Resulting wave"
variable SigX_p
variable numSigInteg //= 6 by default

if(numSigInteg <= 0)
	numSigInteg = 6 //to tune
endif
variable minNumStepsPerSigma = 7 //to tune
variable numStepsPerSigma

variable np = dimsize(wIn, 0)
variable step = dimdelta(wIn, 0), start = dimoffset(wIn, 0), step_t
variable i = 0, x = start, xt
variable curSigX, halfNumStepsGaus, numStepsGaus, invHalfCurSigXe2
variable jStartInteg, jEndInteg, j, sumExp
variable jArgSpecStart, jArgSpecEnd, jArgSpec
variable invSqrtTwoPi = 1./sqrt(2*Pi)

do
	curSigX = SigX_p*x
	invHalfCurSigXe2 = 0.5/(curSigX*curSigX)
	
	step_t = step
	numStepsPerSigma = curSigX/step_t
	if(numStepsPerSigma < minNumStepsPerSigma)
		step_t = step/round(minNumStepsPerSigma/numStepsPerSigma)
	endif
	
	halfNumStepsGaus = round(numSigInteg*curSigX/step_t)
	if(halfNumStepsGaus < 2)
		halfNumStepsGaus = 2
	endif
	
	jStartInteg = -halfNumStepsGaus
	jArgSpecStart = i - jStartInteg
	if(jArgSpecStart >= np)
		jArgSpecStart = np - 1
		jStartInteg = i - jArgSpecStart
	endif
	
	jEndInteg = halfNumStepsGaus
	jArgSpecEnd = i - jEndInteg
	if(jArgSpecEnd < 0)
		jArgSpecEnd = 0
		jEndInteg = i - jArgSpecEnd
	endif
	
	jArgSpec = i - jStartInteg
	xt = jStartInteg*step_t
	//sumExp = 0.5*(wIn[jArgSpec]*exp(-xt*xt*invHalfCurSigXe2))
	sumExp = 0.5*(wIn(x - xt)*exp(-xt*xt*invHalfCurSigXe2))
	
	xt += step_t
	j = jStartInteg + 1
	do
		jArgSpec = i - j
		//sumExp += wIn[jArgSpec]*exp(-xt*xt*invHalfCurSigXe2)
		sumExp += wIn(x - xt)*exp(-xt*xt*invHalfCurSigXe2)
		
		xt += step_t
		j += 1
	while(j < jEndInteg)
	jArgSpec = i - jEndInteg
	//sumExp += 0.5*(wIn[jArgSpec]*exp(-xt*xt*invHalfCurSigXe2))
	sumExp += 0.5*(wIn(x - xt)*exp(-xt*xt*invHalfCurSigXe2))
	
	wOut[i] = sumExp*step_t*invSqrtTwoPi/curSigX
		//DoUpdate
		
	x += step
	i += 1
while(i < np)
end

//+++++++++++++++++++++++++++++++++++++++
//Makes convolution of two 2D waves 
//places the result into w1
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiConvWaves2D(w1, w2)
wave w1, w2

variable NxOrigW1 = DimSize(w1, 0)
variable NyOrigW1 = DimSize(w1, 1)
//variable NxOrigW1M1 = NxOrigW1 - 1
//variable NyOrigW1M1 = NyOrigW1 - 1
//variable NxActW1 = 2*trunc(0.5*NxOrigW1 + 0.000001)
//variable NyActW1 = 2*trunc(0.5*NyOrigW1 + 0.000001)
variable StepXW1 = DimDelta(w1, 0)
variable StepYW1 = DimDelta(w1, 1)
//make/O/C/N=(NxActW1, NyActW1) AuxWaveUtiConvWavesW1C
//AuxWaveUtiConvWavesW1C = w1[srwUtiNumOrMax(p, NxOrigW1M1)][srwUtiNumOrMax(q, NyOrigW1M1)]
//SetScale/P x -0.5*StepXW1*NxActW1, StepXW1, WaveUnits(w1, 0), AuxWaveUtiConvWavesW1C
//SetScale/P y -0.5*StepYW1*NyActW1, StepYW1, WaveUnits(w1, 1), AuxWaveUtiConvWavesW1C
make/O/C/N=(NxOrigW1, NyOrigW1) AuxWaveUtiConvWavesW1C
AuxWaveUtiConvWavesW1C = w1[p][q]
SetScale/P x DimOffset(w1, 0), StepXW1, WaveUnits(w1, 0), AuxWaveUtiConvWavesW1C
SetScale/P y DimOffset(w1, 1), StepYW1, WaveUnits(w1, 1), AuxWaveUtiConvWavesW1C

variable NxOrigW2 = DimSize(w2, 0)
variable NyOrigW2 = DimSize(w2, 1)
//variable NxOrigW2M1 = NxOrigW2 - 1
//variable NyOrigW2M1 = NyOrigW2 - 1
//variable NxActW2 = 2*trunc(0.5*NxOrigW2 + 0.000001)
//variable NyActW2= 2*trunc(0.5*NyOrigW2 + 0.000001)
variable StepXW2 = DimDelta(w2, 0)
variable StepYW2 = DimDelta(w2, 1)
//make/O/C/N=(NxActW2, NyActW2) AuxWaveUtiConvWavesW2C
//AuxWaveUtiConvWavesW2C = w2[srwUtiNumOrMax(p, NxOrigW2M1)][srwUtiNumOrMax(q, NyOrigW2M1)]
//SetScale/P x -0.5*StepXW2*NxActW2, StepXW2, WaveUnits(w2, 0), AuxWaveUtiConvWavesW2C
//SetScale/P y -0.5*StepYW2*NyActW2, StepYW2, WaveUnits(w2, 1), AuxWaveUtiConvWavesW2C
make/O/C/N=(NxOrigW2, NyOrigW2) AuxWaveUtiConvWavesW2C
AuxWaveUtiConvWavesW2C = w2[p][q]
SetScale/P x DimOffset(w2, 0), StepXW2, WaveUnits(w2, 0), AuxWaveUtiConvWavesW2C
SetScale/P y DimOffset(w2, 1), StepYW2, WaveUnits(w2, 1), AuxWaveUtiConvWavesW2C

srFFT2D(AuxWaveUtiConvWavesW1C, 1)
srFFT2D(AuxWaveUtiConvWavesW2C, 1)
AuxWaveUtiConvWavesW1C *= AuxWaveUtiConvWavesW2C(x)(y)
srFFT2D(AuxWaveUtiConvWavesW1C, -1)
w1 = real(AuxWaveUtiConvWavesW1C[p][q])
KillWaves/Z AuxWaveUtiConvWavesW1C, AuxWaveUtiConvWavesW2C
end

//+++++++++++++++++++++++++++++++++++++++
//Makes convolution of two complex 1D waves 
//places the result into w1
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiConvWavesC1D(w1, w2)
wave/C w1, w2

duplicate/O w2 AuxWaveUtiConvWavesW2C

srFFT1D(w1, 1)
srFFT1D(AuxWaveUtiConvWavesW2C, 1)
w1 *= AuxWaveUtiConvWavesW2C(x)
srFFT1D(w1, -1)
KillWaves/Z AuxWaveUtiConvWavesW2C
end

//+++++++++++++++++++++++++++++++++++++++
//Finds average value of a scaled 2D wave around a given point
//Added to make ~same type of calculation as SrwUtiConvWaveWithInterv2D(),
//but eventually more robust. Can be used e.g. to take into account 
//the effect of finite pixel size of a detector.
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiAvgValNearPtWave2D(w, x, y, dx, dy, [xwStart, xwStep, ywStart, ywStep])
wave w
variable x, y, dx, dy
variable xwStart, xwStep, ywStart, ywStep //mesh parameters of the w wave, to avoid repetitive calls to dimoffset, dimdelta

if(ParamIsDefault(xwStart)) 
	xwStart = dimoffset(w, 0)
endif
if(ParamIsDefault(xwStep)) 
	xwStep = dimdelta(w, 0)
endif
if(ParamIsDefault(ywStart)) 
	ywStart = dimoffset(w, 1)
endif
if(ParamIsDefault(ywStep)) 
	ywStep = dimdelta(w, 1)
endif

variable dxHalf = 0.5*dx
variable xStart = x - dxHalf, xEnd = x + dxHalf

variable fpStart = (xStart - xwStart)/xwStep
variable pStart = trunc(fpStart)
if(fpStart - pStart > 0)
	pStart += 1
endif
variable fpEnd = (xEnd - xwStart)/xwStep
variable pEnd = trunc(fpEnd)
if(fpEnd - pEnd == 0)
	pEnd -= 1
endif
if(pEnd < pStart)
	return 0
endif

variable dyHalf = 0.5*dy
variable yStart = y - dyHalf, yEnd = y + dyHalf

variable fqStart = (yStart - ywStart)/ywStep
variable qStart = trunc(fqStart)
if(fqStart - qStart > 0)
	qStart += 1
endif
variable fqEnd = (yEnd - ywStart)/ywStep
variable qEnd = trunc(fqEnd)
if(fqEnd - qEnd == 0)
	qEnd -= 1
endif
if(qEnd < qStart)
	return 0
endif

//print pStart, pEnd, qStart, qEnd
variable ip, iq, sumVal = 0
for(iq = qStart; iq <= qEnd; iq += 1)
	for(ip = pStart; ip <= pEnd; ip += 1)
		sumVal += w[ip][iq]
	endfor
endfor
variable np = (pEnd - pStart + 1)*(qEnd - qStart + 1)
//print np
return sumVal/np
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates statistical relative standard deviation of one scaled wave from another
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiStdDevWaves1D(w1, w2, [iStart1, iEnd1, iSkipStart1, iSkipEnd1])
wave w1, w2
variable iStart1, iEnd1
variable iSkipStart1, iSkipEnd1

variable npw1 = dimsize(w1, 0)
variable xStart1 = dimoffset(w1, 0)
variable xStep1= dimdelta(w1, 0)

if(ParamIsDefault(iStart1))
	iStart1 = 0
endif
if(ParamIsDefault(iEnd1))
	iEnd1 = npw1 - 1
endif

if(ParamIsDefault(iSkipStart1))
	iSkipStart1 = -1
endif
if(ParamIsDefault(iSkipEnd1))
	iSkipEnd1 = -1
endif

variable treatSkip = 0
if(iSkipStart1 <= iSkipEnd1)
	if((iStart1 <= iSkipStart1) && (iSkipStart1 <= iEnd1))
		if((iStart1 <= iSkipEnd1) && (iSkipEnd1 <= iEnd1))
			treatSkip = 1
		endif
	endif
endif

variable i, s2 = 0, v2, dv, useThisPt
variable xx = xStart1 + iStart1*xStep1

for(i = iStart1; i <= iEnd1; i += 1)
	useThisPt = 1
	if(treatSkip > 0)
		if((iSkipStart1 <= i) && (i <= iSkipEnd1))
			useThisPt = 0
		endif
	endif

	if(useThisPt > 0)
		v2 = w2(xx)
		dv = w1(xx) - v2
		s2 += dv*dv/(v2*v2)
		//s2 += abs(dv/v1)
	endif
	
	xx += xStep1
endfor
variable np = iEnd1 - iStart1 + 1
return sqrt(s2/np)
//return s2/np
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave partially
//by O.Marcouille
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiIntWave2D(wave2Dname,x_min,x_max,y_min,y_max)
String wave2Dname
Variable/D x_min,x_max,y_min,y_max
Variable/D delta_x,delta_y,delta,sum
Variable/D index_x_min,index_x,index_x_max,index_y_min,index_y,index_y_max

if((x_min == x_max) %| (y_min == y_max))
	return 0
endif

Duplicate/O/R=(x_min,x_max)(y_min,y_max) $wave2Dname wave2DtempIn,wave2DtempOut

//print wave2Dname
delta_x=DimDelta(wave2DtempIn,0)
delta_y=DimDelta(wave2DtempIn,1)
delta=delta_x*delta_y

index_x_min=trunc((x_min-DimOffset(wave2DtempIn,0))/delta_x + 0.01)
index_x_max=trunc((x_max-DimOffset(wave2DtempIn,0))/delta_x + 0.01)
index_y_min=trunc((y_min-DimOffset(wave2DtempIn,1))/delta_y + 0.01)
index_y_max=trunc((y_max-DimOffset(wave2DtempIn,1))/delta_y + 0.01)
wave2DtempOut[][]=0
do
	index_x=index_x_min
	do
		if ((index_x>0)%&(index_y>0))
			wave2DtempOut[index_x][index_y]=wave2DtempOut[index_x-1][index_y]
			wave2DtempOut[index_x][index_y]+=0.5*(wave2DtempIn[index_x-1][index_y-1]+wave2DtempIn[index_x][index_y])*delta
		endif
		index_x+=1
	while (index_x<=index_x_max)
		
	wave2DtempOut[0][index_y+1]=wave2DtempOut[index_x_max][index_y]
	index_y+=1
while (index_y<=index_y_max)
sum=wave2DtempOut[index_x_max][index_y_max]	
return sum
End

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave partially
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiIntegWave2D(w2d, x_min, x_max, y_min, y_max)
wave w2d
variable/D x_min,x_max,y_min,y_max

if((x_min == x_max) %| (y_min == y_max))
	return 0
endif
variable xStart = dimoffset(w2d, 0), xStep = dimdelta(w2d, 0), nx = dimsize(w2d, 0)
variable yStart = dimoffset(w2d, 1), yStep = dimdelta(w2d, 1), ny = dimsize(w2d, 1)
if((xStep == 0) %| (yStep == 0))
	return 0
endif

variable ixStart = trunc((x_min - xStart)/xStep)
if(ixStart < 0)
	ixStart = 0
else
	if(ixStart >= nx)
		ixStart = nx - 1
	endif
endif
variable iyStart = trunc((y_min - yStart)/yStep)
if(iyStart < 0)
	iyStart = 0
else
	if(iyStart >= ny)
		iyStart = ny - 1
	endif
endif
variable ixEnd = trunc((x_max - xStart)/xStep - 1e-06) + 1
if(ixEnd < 0)
	ixEnd = 0
else
	if(ixEnd >= nx)
		ixEnd = nx - 1
	endif
endif
variable iyEnd = trunc((y_max - yStart)/yStep - 1e-06) + 1
if(iyEnd < 0)
	iyEnd = 0
else
	if(iyEnd >= ny)
		iyEnd = ny - 1
	endif
endif

if((ixStart >= ixEnd) %| (iyStart >= iyEnd))
	return 0
endif

variable nxi = ixEnd - ixStart + 1
make/O/N=(nxi) wAuxIntegWave2Dx
SetScale/P x, xStart + ixStart*xStep, xStep, "", wAuxIntegWave2Dx
variable nyi = iyEnd - iyStart + 1
make/O/N=(nyi) wAuxIntegWave2Dy
SetScale/P x, yStart + iyStart*yStep, yStep, "", wAuxIntegWave2Dy

variable ix, iy
for(iy = 0; iy < nyi; iy += 1)
	wAuxIntegWave2Dx = w2d[ixStart + p][iyStart + iy]
	integrate/T wAuxIntegWave2Dx
	wAuxIntegWave2Dy[iy] = wAuxIntegWave2Dx(x_max) - wAuxIntegWave2Dx(x_min)
endfor
integrate/T wAuxIntegWave2Dy
variable res = wAuxIntegWave2Dy(y_max) - wAuxIntegWave2Dy(y_min)
killwaves/Z wAuxIntegWave2Dx, wAuxIntegWave2Dy
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave fully
//+++++++++++++++++++++++++++++++++++++++
function srwUtiIntTotWave2D(wave2d)
wave wave2d
wavestats/Q wave2d
return V_avg*dimdelta(wave2d, 0)*dimsize(wave2d, 0)*dimdelta(wave2d, 1)*dimsize(wave2d, 1)
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave fully, using Trapethoidal integration
//+++++++++++++++++++++++++++++++++++++++
function srwUtiIntTotWave2DT(wave2d)
wave wave2d
//wavestats/Q wave2d
//return V_avg*dimdelta(wave2d, 0)*dimsize(wave2d, 0)*dimdelta(wave2d, 1)*dimsize(wave2d, 1)

variable nx = dimsize(wave2d, 0)
variable ny = dimsize(wave2d, 1)
variable ny_mi_1 = ny - 1
variable xStart = dimoffset(wave2d, 0), xStep = dimdelta(wave2d, 0)
variable yStart = dimoffset(wave2d, 1), yStep = dimdelta(wave2d, 1)

make/O/N=(nx) auxWave1D
SetScale/P x xStart, xStep,"", auxWave1D
variable iy = 0, sumVal = 0, mult = 1
do
	auxWave1D = wave2d[p][iy]
	integrate/T auxWave1D
	if((iy != 0) %& (iy != ny_mi_1))
		mult = 1
	else 
		mult = 0.5
	endif
	sumVal += mult*auxWave1D[nx - 1]
	iy += 1
while(iy < ny)
variable res = yStep*sumVal
killwaves/Z auxWave1D
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave fully within given limits
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiIntWave2D_Aux(wave2Dname,x_min,x_max,y_min,y_max)
string wave2Dname
variable x_min,x_max,y_min,y_max
//Variable/D delta_x,delta_y,delta,sum
//Variable/D index_x_min,index_x,index_x_max,index_y_min,index_y,index_y_max

if((x_min == x_max) %| (y_min == y_max))
	return 0
endif

duplicate/O/R=(x_min,x_max)(y_min,y_max) $wave2Dname wave2DtempIn //,wave2DtempOut
variable res = srwUtiIntTotWave2D(wave2DtempIn)
killwaves/Z wave2DtempIn
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 3D wave over 2D
//+++++++++++++++++++++++++++++++++++++++
function srwUtiIntWave3Dvs2D(wave3d, p0)
wave wave3d
variable p0
make/O/N=(dimsize(wave3d, 1), dimsize(wave3d, 2)) auxwaveaa
auxwaveaa = wave3d[p0][p][q]
wavestats/Q auxwaveaa
variable res = V_avg*dimdelta(wave3d, 1)*dimsize(wave3d, 1)*dimdelta(wave3d, 2)*dimsize(wave3d, 2)
killwaves/Z auxwaveaa
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Sets up a scaled 2D wave from 1D wave
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMakeWave2DFrom1D(NameW2D, NameW1D, Np, Nq, pStart, qStart, pStep, qStep)
string NameW2D = srwUtiGetValS("NameW2D", "", "SrwUtiMakeWave2DFrom1D")
string NameW1D = srwUtiGetValS("NameW1D", "", "SrwUtiMakeWave2DFrom1D")
variable Np = srwUtiGetValN("Np", 50, "SrwUtiMakeWave2DFrom1D")
variable pStart = srwUtiGetValN("pStart", 0, "SrwUtiMakeWave2DFrom1D")
variable pStep = srwUtiGetValN("pStep", 1, "SrwUtiMakeWave2DFrom1D")
variable Nq = srwUtiGetValN("Nq", 50, "SrwUtiMakeWave2DFrom1D")
variable qStart = srwUtiGetValN("qStart", 0, "SrwUtiMakeWave2DFrom1D")
variable qStep = srwUtiGetValN("qStep", 1, "SrwUtiMakeWave2DFrom1D")
prompt NameW2D,"Name for 2D wave to produce"
prompt NameW1D,"1D wave",popup Wavelist("*",";", "")
prompt Np,"Number of columns in 2D wave"
prompt pStart,"Initial argument value in each row"
prompt pStep,"Step of argument in each row"
prompt Nq"Number of rows in 2D wave"
prompt qStart,"Initial argument value in each column"
prompt qStep,"Step of argument in each column"
PauseUpdate
Silent 1						|	Importing the Field Component ...

srwUtiSetValS("NameW2D", NameW2D, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValS("NameW1D", NameW1D, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("Np", Np, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("pStart", pStart, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("pStep", pStep, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("Nq", Nq, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("qStart", qStart, "SrwUtiMakeWave2DFrom1D")
srwUtiSetValN("qStep", qStep, "SrwUtiMakeWave2DFrom1D")

variable Ntot = DimSize($NameW1D, 0)
if(Np*Nq != Ntot)
	Abort "The total number of points in the 2D wave should be equal to the number of points in 1D wave"
endif

string WaveUnitsStr = WaveUnits($NameW1D,0)
make/O/N=(Np,Nq) $NameW2D
SetScale/P x pStart, pStep, WaveUnitsStr, $NameW2D
SetScale/P y qStart, qStep, WaveUnitsStr, $NameW2D
$NameW2D = $NameW1D[q*Np + p]
end

//+++++++++++++++++++++++++++++++++++++++
//Applies spatial transformation to a point
//+++++++++++++++++++++++++++++++++++++++
function srwUtiTrfP(wM, wV, wP) 
wave wM, wV, wP
variable xt = wM[0][0]*wP[0] + wM[1][0]*wP[1] + wM[2][0]*wP[2]
variable yt = wM[0][1]*wP[0] + wM[1][1]*wP[1] + wM[2][1]*wP[2]
variable zt = wM[0][2]*wP[0] + wM[1][2]*wP[1] + wM[2][2]*wP[2]

//variable a31 = wM[0][2]
//variable a32 = wM[1][2]
//variable a33 = wM[2][2]
//variable p0 = wP[0]
//variable p1 = wP[1]
//variable p2 = wP[2]

wP[0] = xt + wV[0]
wP[1] = yt + wV[1]
wP[2] = zt + wV[2]
end

//+++++++++++++++++++++++++++++++++++++++
//Applies spatial transformation to a vector
//+++++++++++++++++++++++++++++++++++++++
function srwUtiTrfV(wM, wP) 
wave wM, wP
variable xt = wM[0][0]*wP[0] + wM[1][0]*wP[1] + wM[2][0]*wP[2]
variable yt = wM[0][1]*wP[0] + wM[1][1]*wP[1] + wM[2][1]*wP[2]
variable zt = wM[0][2]*wP[0] + wM[1][2]*wP[1] + wM[2][2]*wP[2]

wP[0] = xt
wP[1] = yt
wP[2] = zt
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculate Norm and eventually Normalize vector
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiVectNorm(v, newNorm)
wave v
variable newNorm //<=0- don't normalize v; >0- make newNorm a new normal of v

variable nDim = dimsize(v, 0)
if(nDim <= 0)
	return 0
endif
variable sumE2 = 0, i
for(i = 0; i < nDim; i += 1)
	sumE2 += v[i]*v[i]
endfor
variable oldNorm = sqrt(sumE2)
if((oldNorm > 0) %& (newNorm > 0))
	v *= newNorm/oldNorm
endif
return oldNorm
end 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculate Norm and eventually Normalize vector
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiVectScalProd(v1, v2)
wave v1, v2
variable nDim1 = dimsize(v1, 0), nDim2 = dimsize(v2, 0), i
variable nMin = nDim1
if(nMin > nDim2)
	nMin = nDim2
endif
variable scalProd = 0
for(i = 0; i < nMIn; i += 1)
	scalProd += v1[i]*v2[i]
endfor
return scalProd
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates Intersection of a straight line with Torus
//At Input: wR0 is the point defining the line
//At Output: wR0 contains coordinates of the intersection point
//+++++++++++++++++++++++++++++++++++++++
function srwUtiOptIntersTorAndLine(rt, rs, wV, wR0) 
variable rt, rs
wave wR0, wV

variable ax = wV[0]/wV[1], az = wV[2]/wV[1]
variable b = rt/rs
variable x0 = wR0[0], y0 = wR0[1], z0 = wR0[2]

variable buf1 = 1 + ax*ax + az*az*b
variable dx = x0 - ax*y0, dz = z0 - az*y0
variable buf2 = rt + ax*dx + az*b*dz
variable buf3 = dx*dx + b*dz*dz

wR0[1] = (-buf2 + sqrt(buf2*buf2 - buf1*buf3))/buf1 // yi 
variable dy = wR0[1] - y0
wR0[0] = ax*dy + x0 // xi
wR0[2] = az*dy + z0 // zi
end

//+++++++++++++++++++++++++++++++++++++++
//Round number
//+++++++++++++++++++++++++++++++++++++++
function srwAuxRoundDecNum(num, digs)
variable num, digs
if(num == 0)
	return 0
endif
return (10^(-digs+trunc(log(num))))*round((10^(digs-trunc(log(num))))*num)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns non-zero number
//+++++++++++++++++++++++++++++++++++++++
function srwAuxRetNonZero(Number, ZeroTol)
variable Number, ZeroTol
if(Number != 0)
	return Number
else
	return ZeroTol
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Returns zero if negative
//+++++++++++++++++++++++++++++++++++++++
function srwAuxRetPosOrZero(Number)
variable Number
if(Number >= 0)
	return Number
else
	return 0
endif
end

//+++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++
function SrwAuxFindWaveLevel(v, w)
variable v
wave w
FindLevel/Q w, v
return V_LevelX
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave vs one dimension
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiIntegWave2Dvs1D(wnamesect, wname, DimNoToInteg, StartIntPos)
string wnamesect=srwUtiGetValS("wnamesect", "", "AuxIntegWave2DAlong1D")
string wname=srwUtiGetValS("wname", "", "AuxIntegWave2DAlong1D")
variable DimNoToInteg=srwUtiGetValN("DimNoToInteg", 1, "AuxIntegWave2DAlong1D")
variable StartIntPos=srwUtiGetValN("StartIntPos", 0, "AuxIntegWave2DAlong1D")
prompt wnamesect, "Name of Profile Wave to Create"
prompt wname, "Name of Wave to Integrate"
prompt DimNoToInteg, "Dimension Along Which to Integrate",popup "1st;2nd"
prompt StartIntPos, "Initial Argument Value in Another Dimension from which to Integrate"//"Initial Position for Integration"

if(DimNoToInteg > 2)
	abort "Incorrect dimension number"
endif

srwUtiSetValS("wnamesect", wnamesect, "AuxIntegWave2DAlong1D")
srwUtiSetValS("wname", wname, "AuxIntegWave2DAlong1D")
srwUtiSetValN("DimNoToInteg", DimNoToInteg, "AuxIntegWave2DAlong1D")
srwUtiSetValN("StartIntPos", StartIntPos, "AuxIntegWave2DAlong1D")

DimNoToInteg -= 1
variable DimNoToLeave = 1
if(DimNoToInteg == 1)
	DimNoToLeave = 0
endif

variable npts = dimsize($wname, DimNoToLeave)
make/O/N=(npts) $wnamesect
SetScale/P x dimoffset($wname,DimNoToLeave),dimdelta($wname,DimNoToLeave),waveunits($wname,DimNoToLeave), $wnamesect
$wnamesect = 1000*srwUtiIntegWave2Dvs1DatPos($wname, DimNoToInteg, p, StartIntPos)
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave vs one dimension interval
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiIntegWave2Dvs1DInterval(wnamesect, wname, DimNoToInteg, StartIntPos, EndIntPos)
string wnamesect=srwUtiGetValS("wnamesect", "", "AuxIntegWave2DAlong1D")
string wname=srwUtiGetValS("wname", "", "AuxIntegWave2DAlong1D")
variable DimNoToInteg=srwUtiGetValN("DimNoToInteg", 1, "AuxIntegWave2DAlong1D")
variable StartIntPos=srwUtiGetValN("StartIntPos", 0, "AuxIntegWave2DAlong1D")
variable EndIntPos=srwUtiGetValN("EndIntPos", 0, "AuxIntegWave2DAlong1D")
prompt wnamesect, "Name of Profile Wave to Create"
prompt wname, "2D Wave to Integrate", popup Wavelist("*",";","TEXT:0,DIMS:2")
prompt DimNoToInteg, "Dimension Along Which to Integrate",popup "1st;2nd"
prompt StartIntPos, "Initial Argument Value in Another Dimension from which to Integrate"//"Initial Position for Integration"
prompt EndIntPos, "Final Argument Value in Another Dimension from which to Integrate"//"Final Position for Integration"
Silent 1						|	 ...
PauseUpdate

if(DimNoToInteg > 2)
	abort "Incorrect dimension number"
endif

srwUtiSetValS("wnamesect", wnamesect, "AuxIntegWave2DAlong1D")
srwUtiSetValS("wname", wname, "AuxIntegWave2DAlong1D")
srwUtiSetValN("DimNoToInteg", DimNoToInteg, "AuxIntegWave2DAlong1D")
srwUtiSetValN("StartIntPos", StartIntPos, "AuxIntegWave2DAlong1D")
srwUtiSetValN("EndIntPos", EndIntPos, "AuxIntegWave2DAlong1D")

DimNoToInteg -= 1
variable DimNoToLeave = 1
if(DimNoToInteg == 1)
	DimNoToLeave = 0
endif

variable npts = dimsize($wname, DimNoToLeave)
make/O/N=(npts) $wnamesect
SetScale/P x dimoffset($wname,DimNoToLeave),dimdelta($wname,DimNoToLeave),waveunits($wname,DimNoToLeave), $wnamesect
$wnamesect = srwUtiIntegWave2Dvs1Dinterv($wname, DimNoToInteg, p, StartIntPos, EndIntPos)
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave vs one dimension at given position
//+++++++++++++++++++++++++++++++++++++++
function srwUtiIntegWave2Dvs1DatPos(inwave, DimNoToInteg, ind, StartPos)
wave inwave
variable DimNoToInteg, ind, StartPos

variable npts = dimsize(inwave, DimNoToInteg)
make/O/N=(npts) auxwave
SetScale/P x dimoffset(inwave,DimNoToInteg),dimdelta(inwave,DimNoToInteg),"", auxwave

if(DimNoToInteg == 0)
	auxwave = inwave[p][ind]*srwUtiStep(x - StartPos)
else 
	auxwave = inwave[ind][p]*srwUtiStep(x - StartPos)
endif
integrate/T auxwave
variable res = auxwave[npts - 1]
killwaves/Z auxwave
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Integrate 2D wave vs one dimension at given position
//+++++++++++++++++++++++++++++++++++++++
function srwUtiIntegWave2Dvs1Dinterv(inwave, DimNoToInteg, ind, StartPos, EndPos)
wave inwave
variable DimNoToInteg, ind, StartPos, EndPos

variable npts = dimsize(inwave, DimNoToInteg)
make/O/N=(npts) auxwave
SetScale/P x dimoffset(inwave,DimNoToInteg),dimdelta(inwave,DimNoToInteg),"", auxwave

if(DimNoToInteg == 0)
	auxwave = inwave[p][ind]*srwUtiNonZeroIntervB(x, StartPos, EndPos)
	//srwUtiStep(x - StartPos)
else 
	auxwave = inwave[ind][p]*srwUtiNonZeroIntervB(x, StartPos, EndPos)
	//auxwave = inwave[ind][p]*srwUtiStep(x - StartPos)
endif
integrate/T auxwave
variable res = auxwave[npts - 1]
killwaves/Z auxwave
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Creates 1D Numerical wave
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiCreateNumWave1D(nmWave, Type, Npts, Start, Step, ArgUnits, ArgLabel, ValUnits)
string nmWave
string Type
variable Npts
variable Start
variable Step
string ArgUnits
string ArgLabel
string ValUnits
string ValLabel

if(cmpstr(Type,"D") == 0)
	Make/D/O/N=(Npts) $nmWave
else
	if(cmpstr(Type,"C") == 0)
		Make/C/O/N=(Npts) $nmWave
	else
		Make/O/N=(Npts) $nmWave
	endif
endif

$nmWave = 0
SetScale/P x, Start, Step, ArgUnits, $nmWave
SetScale d, 0, 0, ValUnits, $nmWave

//string strExe = "SetDimLabel 0, -1, " + ArgLabel + ", '" + nmWave + "'"
//execute strExe
SetDimLabel 0, -1, $ArgLabel, $nmWave
//SetDimLabel -1, -1, $ValLabel, $nmWave
end

//+++++++++++++++++++++++++++++++++++++++
//(Re-)Creates Data Wave Information structure 
//(text wave where auxiliary information about data waves is stored)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiDataWaveInfCreate()
if(exists("SrwDataWaveInf") != 1)
	make/O/T/N=(50, 5) SrwDataWaveInf
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Enters information into Data Wave Information structure 
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiDataWaveInfStore(NameWave, DataType, DataDescr)
string NameWave, DataType, DataDescr
Silent 1
PauseUpdate

if((strlen(NameWave) <= 0) %| (strlen(DataType) <= 0))
	return
endif
if(strlen(DataDescr) <= 0)
	return
endif
SrwUtiDataWaveInfCreate()

if(cmpstr(DataDescr, SrwPUnitSpAngFluxPerUnSurf) == 0)
	DataDescr = SrwPUnitSpAngFluxPerUnSurf1
endif
if(cmpstr(DataDescr, SrwPUnitSpAngFluxPerUnAngle) == 0)
	DataDescr = SrwPUnitSpAngFluxPerUnAngle1
endif
if(cmpstr(DataDescr, SrwPUnitBrilliance) == 0)
	DataDescr = SrwPUnitBrilliance1
endif
if(cmpstr(DataDescr, SrwPUnitPowDen) == 0)
	DataDescr = SrwPUnitPowDen1
endif
if(cmpstr(DataDescr, SrwPUnitElectricField) == 0)
	DataDescr = SrwPUnitElectricField1
endif

variable AmOfLines = dimsize(SrwDataWaveInf, 0)
variable i = 0, NameWaveExists = 0
string ExistingNameWave = ""
do
	ExistingNameWave = SrwDataWaveInf[i][0]
	if(strlen(ExistingNameWave) <= 0)
		break
	endif
	if(cmpstr(ExistingNameWave, NameWave) == 0)
		NameWaveExists = 1
		break
	endif
	i += 1
while(i < AmOfLines)

string SeparSymb = ":"
variable SrwDataWaveInfNumCol = dimsize(SrwDataWaveInf, 1)
string DataTypeWithSep = DataType + SeparSymb

string NewRecordStr = DataType + SeparSymb + DataDescr
if(NameWaveExists == 0) // wave does not exist
	if(i >= AmOfLines)  // no more free lines
		variable NewAmOfLines = AmOfLines + 10
		redimension/N=(NewAmOfLines, SrwDataWaveInfNumCol) SrwDataWaveInf
	endif
	SrwDataWaveInf[i][0] = NameWave
	SrwDataWaveInf[i][1] = NewRecordStr
else
	variable k = 1, NewRecordEntered = 0
	string CurStr = ""
	do
		CurStr = SrwDataWaveInf[i][k]
		if(strlen(CurStr) <= 0)
			break
		endif
		if(strsearch(CurStr, DataTypeWithSep, 0) == 0) // found
			SrwDataWaveInf[i][k] = NewRecordStr
			NewRecordEntered = 1
			break
		endif
		k += 1
	while(k < SrwDataWaveInfNumCol)
	if(NewRecordEntered == 0)
		if(k >= SrwDataWaveInfNumCol)
			variable NewAmOfCols = SrwDataWaveInfNumCol + 1
			redimension/N=(AmOfLines, NewAmOfCols) SrwDataWaveInf
		endif
		SrwDataWaveInf[i][k] = NewRecordStr
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Extracts information from Data Wave Information structure 
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiDataWaveInfGet(NameWave, DataType)
string NameWave, DataType

string NameDataWaveInf = "SrwDataWaveInf"
if(exists(NameDataWaveInf) != 1)
	return ""
endif

wave/T wDataWaveInf = $NameDataWaveInf
variable AmOfLines = dimsize(wDataWaveInf, 0)
variable i = 0, NameWaveExists = 0
string ExistingNameWave
do
	ExistingNameWave = wDataWaveInf[i][0]
	if(strlen(ExistingNameWave) <= 0)
		break
	endif
	if(cmpstr(ExistingNameWave, NameWave) == 0)
		NameWaveExists = 1
		break
	endif
	i += 1
while(i < AmOfLines)

if(NameWaveExists == 0)
	return ""
endif

string SeparSymb = ":"
variable SrwDataWaveInfNumCol = dimsize(wDataWaveInf, 1)
string DataTypeWithSep = DataType + SeparSymb

variable k = 1
string ExistingRecord = ""
do
	ExistingRecord = wDataWaveInf[i][k]
	if(strsearch(ExistingRecord, DataTypeWithSep, 0) == 0) // found
		break
	endif
	k += 1
while(k < SrwDataWaveInfNumCol)

variable LenExistingRecord = strlen(ExistingRecord)
if(LenExistingRecord == 0)
	return ""
endif

variable SepPos = strsearch(ExistingRecord, SeparSymb, 0)
if((SepPos < 0) %| (SepPos >= LenExistingRecord))
	return ""
endif

string OutStr = ""
variable j = SepPos + 1
do
	OutStr += ExistingRecord[j]
	j += 1
while(j < LenExistingRecord)
return OutStr
end

//+++++++++++++++++++++++++++++++++++++++
//Simplest bilinear interpolation in 2D
//+++++++++++++++++++++++++++++++++++++++
function srwUtiInterp2DBilin(x, y, w)
wave w
variable x, y

variable xmin = dimoffset(w, 0)
variable nx = dimsize(w, 0)
variable xstep = dimdelta(w, 0)
variable xmax = xmin + (nx - 1)*xstep

variable ymin = dimoffset(w, 1)
variable ny = dimsize(w, 1)
variable ystep = dimdelta(w, 1)
variable ymax = ymin + (ny - 1)*ystep

//OC13012020 (uncommented)
if(x < xmin)
	x = xmin
endif
if(x > xmax)
	x = xmax
endif
if(y < ymin)
	y = ymin
endif
if(y > ymax)
	y = ymax
endif
//OC23112017
//if((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax))
//	return 0
//endif

variable x0 = xmin + trunc((x - xmin)/xstep)*xstep
if(x0 >= xmax)
	x0 = xmax - xstep
endif
variable x1 = x0 + xstep

variable y0 = ymin + trunc((y - ymin)/ystep)*ystep
if(y0 >= ymax)
	y0 = ymax - ystep
endif
variable y1 = y0 + ystep

variable t = (x - x0)/xstep, u = (y - y0)/ystep
return (1 - t)*(1 - u)*(w(x0)(y0)) + t*(1 - u)*(w(x1)(y0)) + t*u*(w(x1)(y1)) + (1 - t)*u*(w(x0)(y1))
end

//+++++++++++++++++++++++++++++++++++++++
//Line-by-line spline interpolation in 2D
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiInterp2DNonLin(nmWold, nmWnew)
string nmWold, nmWnew
Silent 1						|	 ...
PauseUpdate

variable nxOld = dimsize($nmWold, 0)
variable xStartOld = dimoffset($nmWold, 0)
variable xStepOld = dimdelta($nmWold, 0)
variable nyOld = dimsize($nmWold, 1)
variable yStartOld = dimoffset($nmWold, 1)
variable yStepOld = dimdelta($nmWold, 1)

variable nxNew = dimsize($nmWnew, 0)
variable xStartNew = dimoffset($nmWnew, 0)
variable xStepNew = dimdelta($nmWnew, 0)
variable nyNew = dimsize($nmWnew, 1)
variable yStartNew = dimoffset($nmWnew, 1)
variable yStepNew = dimdelta($nmWnew, 1)

make/O/N=(nxOld) wAuxOldUtiInterp2DNonLin_x
SetScale/P x xStartOld,xStepOld,"", wAuxOldUtiInterp2DNonLin_x
make/O/N=(nyOld) wAuxOldUtiInterp2DNonLin_y
SetScale/P y yStartOld,yStepOld,"", wAuxOldUtiInterp2DNonLin_y

make/O/N=(nxNew) wAuxNewUtiInterp2DNonLin_x
SetScale/P x xStartNew,xStepNew,"", wAuxNewUtiInterp2DNonLin_x
make/O/N=(nyNew) wAuxNewUtiInterp2DNonLin_y
SetScale/P y yStartNew,yStepNew,"", wAuxNewUtiInterp2DNonLin_y

make/O/N=(nxOld, nyNew) wAuxOldNewUtiInterp2DNonLin_xy
SetScale/P x xStartOld,xStepOld,"", wAuxOldNewUtiInterp2DNonLin_xy
SetScale/P y yStartNew,yStepNew,"", wAuxOldNewUtiInterp2DNonLin_xy

make/O/N=(nxNew, nyOld) wAuxNewOldUtiInterp2DNonLin_xy
SetScale/P x xStartNew,xStepNew,"", wAuxNewOldUtiInterp2DNonLin_xy
SetScale/P y yStartOld,yStepOld,"", wAuxNewOldUtiInterp2DNonLin_xy

make/O/N=(nxNew, nyNew) wAuxNewNewUtiInterp2DNonLin_xy
SetScale/P x xStartNew,xStepNew,"", wAuxNewNewUtiInterp2DNonLin_xy
SetScale/P y yStartNew,yStepNew,"", wAuxNewNewUtiInterp2DNonLin_xy

variable ix = 0
do
	wAuxOldUtiInterp2DNonLin_y = $nmWold[ix][p]
	Interpolate2/T=2/N=(nyNew)/E=2/Y=wAuxNewUtiInterp2DNonLin_y wAuxOldUtiInterp2DNonLin_y
	wAuxOldNewUtiInterp2DNonLin_xy[ix][] = wAuxNewUtiInterp2DNonLin_y[q]
	ix += 1
while(ix < nxOld)
variable iy = 0
do
	wAuxOldUtiInterp2DNonLin_x = wAuxOldNewUtiInterp2DNonLin_xy[p][iy]
	Interpolate2/T=2/N=(nxNew)/E=2/Y=wAuxNewUtiInterp2DNonLin_x wAuxOldUtiInterp2DNonLin_x
	wAuxNewNewUtiInterp2DNonLin_xy[][iy] = wAuxNewUtiInterp2DNonLin_x[p]
	iy += 1
while(iy < nyNew)

iy = 0
do
	wAuxOldUtiInterp2DNonLin_x = $nmWold[p][iy]
	Interpolate2/T=2/N=(nxNew)/E=2/Y=wAuxNewUtiInterp2DNonLin_x wAuxOldUtiInterp2DNonLin_x
	wAuxNewOldUtiInterp2DNonLin_xy[][iy] = wAuxNewUtiInterp2DNonLin_x[p]
	iy += 1
while(iy < nyOld)
ix = 0
do
	wAuxOldUtiInterp2DNonLin_y = wAuxNewOldUtiInterp2DNonLin_xy[ix][p]
	Interpolate2/T=2/N=(nyNew)/E=2/Y=wAuxNewUtiInterp2DNonLin_y wAuxOldUtiInterp2DNonLin_y
	$nmWnew[ix][] = wAuxNewUtiInterp2DNonLin_y[q]
	ix += 1
while(ix < nxNew)

$nmWnew = 0.5*(wAuxNewNewUtiInterp2DNonLin_xy[p][q] + $nmWnew[p][q])
killwaves/Z wAuxNewNewUtiInterp2DNonLin_xy, wAuxNewOldUtiInterp2DNonLin_xy, wAuxOldNewUtiInterp2DNonLin_xy
killwaves/Z wAuxOldUtiInterp2DNonLin_x, wAuxOldUtiInterp2DNonLin_y
killwaves/Z wAuxNewUtiInterp2DNonLin_x, wAuxNewUtiInterp2DNonLin_y
end

//+++++++++++++++++++++++++++++++++++++++
//Compare 2 scaled waves and find max values
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiFindMaxValuesIn2Waves(NameWaveRes, NameWave1, NameWave2)
string NameWaveRes, NameWave1, NameWave2
$NameWaveRes = srwUtiMaxOfTwo($NameWave1(x), $NameWave2(x))
end

//+++++++++++++++++++++++++++++++++++++++
//Compare 2 scaled waves and find max values
//+++++++++++++++++++++++++++++++++++++++
function srwUtiWfrMultExpLin(nameWfr, xc_m, rx_m, yc_m, ry_m, phEn_eV)
string nameWfr
variable xc_m, rx_m, yc_m, ry_m, phEn_eV

wave/T wfr = $nameWfr
string nameEx = wfr[0], nameEz = wfr[1]
wave/C wEx = $nameEx
wave/C wEz = $nameEz

variable wavelength_m =1.239854e-06/phEn_eV
variable k_m = 2*Pi/wavelength_m
variable constX = k_m*xc_m/rx_m, constY = k_m*yc_m/ry_m

wEx *= cmplx(cos(constX*y + constY*z), sin(constX*y + constY*z))
wEz *= cmplx(cos(constX*y + constY*z), sin(constX*y + constY*z))
end

//+++++++++++++++++++++++++++++++++++++++
//Find local maximum or minimum of 1D wave.
//The returned complex number contains: 
//	Re - min / max value
//	Im - min / max argument
//min_or_max == 1 -> search for minimum
//min_or_max == 2 -> search for maximum
//+++++++++++++++++++++++++++++++++++++++
function/C srwUtiFindLocalPeak1D(wData, min_or_max, xCen, xRange)
wave wData
variable min_or_max, xCen, xRange

WaveStats/Q/R=(xCen - 0.5*xRange, xCen + 0.5*xRange) wData
variable/C res
if(min_or_max == 1)
	res = cmplx(V_min, V_minloc)
else
	res = cmplx(V_max, V_maxloc)
endif
return res
end

//+++++++++++++++++++++++++++++++++++++++
//Find harmonic peaks in undulator spectrum
//wSpec is radiation spectral flux density vs photon energy in eV
//The results are stored in 2D wave wRes
//+++++++++++++++++++++++++++++++++++++++
function srwUtiFindSpecPeaks(wSpec, e1_eV, de_eV, nHarmBeg, nHarmEnd, wRes)
wave wSpec, wRes
variable e1_eV, de_eV, nHarmBeg, nHarmEnd

variable amOfHarm = nHarmEnd - nHarmBeg + 1
if(amOfHarm <= 0)
	return 0
endif

redimension/N=(amOfHarm, 3) wRes
variable i, ic = 0, ec
variable/C auxC
for(i = nHarmBeg; i <= nHarmEnd; i += 1)
	ec = e1_eV*i
	auxC = srwUtiFindLocalPeak1D(wSpec, 2, ec, de_eV)
	wRes[ic][0] = i
	wRes[ic][1] = imag(auxC)
	wRes[ic][2] = real(auxC)
	ic += 1
endfor
end

//+++++++++++++++++++++++++++++++++++++++
//Generates 2D Gaussian random number (with a cross-term)
//+++++++++++++++++++++++++++++++++++++++
function/C srwUtiRandGauss2D(SigXe2, MXXp, SigXpe2)
variable SigXe2, MXXp, SigXpe2
variable mult = 0.5/(SigXe2*SigXpe2 - MXXp*MXXp)
variable B = SigXe2*mult
variable G = SigXpe2*mult
variable A = MXXp*mult

variable SigP = 1/sqrt(2*G)
variable SigQ = sqrt(G/(2*(B*G - A*A)))
variable p = gnoise(SigP)
variable Xp = gnoise(SigQ)
variable X = p + A*Xp/G
return cmplx(X, Xp)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns sub-string of a given string
//+++++++++++++++++++++++++++++++++++++++
function/S srwUtiSubStr(str, iStart, lenSubStr)
string str
variable iStart, lenSubStr

string strRes = ""
if(lenSubStr <= 0)
	return strRes
endif

variable lenStr = strlen(str)
if(lenStr < lenSubStr)
	lenSubStr = lenStr
endif
variable i = 0
do
	strRes += str[iStart + i]
	i += 1
while(i < lenSubStr)
return strRes
end

//+++++++++++++++++++++++++++++++++++++++
//Returns 1 if (x,y) point is inside a tilted centered ellipse
//(on boundary means still inside)
//+++++++++++++++++++++++++++++++++++++++
function srwUtiCheckIfInsideEllipse(xx,yy,a,b,ang)
variable xx, yy, a, b, ang

variable xt = xx*cos(ang) + yy*sin(ang)
variable yt = -xx*sin(ang) + yy*cos(ang)

if((xt < -a) %| (xt > a))
	return 0
endif
if((yt < -b) %| (yt > b))
	return 0
endif

variable aE2 = a*a, bE2 = b*b, ytE2 = yt*yt
variable ytBordE2 = bE2*(1 - xt*xt/aE2)
if(ytE2 < ytBordE2)
	return 1
else
	return 0
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Finds maximum of cubic polynomial passing through 4 points
//Returns cmplx(iArgMax, fMax)
//TO DEBUG!!!
//+++++++++++++++++++++++++++++++++++++++
function/C srwUtiFindMaxFunc4P(f)
wave f
variable a0 = f[0], a1 = -(11*f[0]/6) + 3*f[1] - 3*f[2]/2 + f[3]/3, a2 = f[0] - 5*f[1]/2 + 2*f[2] - f[3]/2, a3 = -f[0]/6 + f[1]/2 - f[2]/2 + f[3]/6
//print a0, a1, a2, a3
//abort

variable d = a2*a2 - 3*a1*a3
if(d < 0)
	return Nan
endif

variable sqrt_d = sqrt(d), a3_3 =  a3*3
variable x1 = (-a2 - sqrt_d)/a3_3, x2 = (-a2 + sqrt_d)/a3_3

variable f1 = f[0], f2 = f[0]
if((0 <= x1) && (x1 <= 3))
	f1 = a0 + (a1 + (a2 + a3*x1)*x1)*x1
else
	x1 = 0
endif

if((0 <= x2) && (x2 <= 3))
	f2 = a0 + (a1 + (a2 + a3*x2)*x2)*x2
else
	x2 = 0
endif

make/O wAuxUtiFindMaxPtPol3 = {f[0], f1, f2, f[3]}, wAuxArgUtiFindMaxPtPol3 = {0, x1, x2, 3}
variable i, fmax=wAuxUtiFindMaxPtPol3[0], imax=0
for(i=0; i<4; i+=1)
	if(fmax < wAuxUtiFindMaxPtPol3[i])
		fmax = wAuxUtiFindMaxPtPol3[i]
		imax = i
	endif
endfor
//print sqrt_d, x1, x2, imax, fmax
//abort
return cmplx(wAuxArgUtiFindMaxPtPol3[imax], fmax)
end

//+++++++++++++++++++++++++++++++++++++++
//Deletes any SRW structure after obtaining confirmation
//+++++++++++++++++++++++++++++++++++++++
//proc SrwUtiDelStructConfirm(StructName)
//string StructName=SrwRadName+SrwRadType
//prompt RadiniName, "Wavefront structure to delete", popup Wavelist("*"+SrwRadType, ";", "")
//Silent 1						|	Deleting the Radiation  structure  ...
//PauseUpdate
//
//if(cmpstr(RadiniName,"_none_")==0)
//	DoAlert 0, SrwPAlertWavefrontNeeded
//	return
//endif
//
//string PromptStr = "Are you sure you want to delete the wavefront structure " + RadiniName + "?"
//DoAlert 1, PromptStr
//if(V_Flag != 1)
//	return
//endif
//end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates "Argument Extension Coefficient" (a, see below) of one curve (wRel) with respect to another (wBase):
//wRel(x) = wBase(x*(1 + a))
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiFindArgExtCoef(wBase, wRel, maxAbsVal, absTol)
wave wBase, wRel //scaled 1D waves
variable maxAbsVal, absTol

string nmWaveBase = NameOfWave(wBase)
string nmWaveRel = NameOfWave(wRel)
string strExe = "duplicate/O " + nmWaveBase + " gwBaseFindArgExtCoef; duplicate/O " + nmWaveRel + " gwRelFindArgExtCoef"
execute strExe
make/O wDummyFindArgExtCoef = {0}
Optimize/Q/L=(-maxAbsVal)/H=(maxAbsVal)/T=(absTol) srwUtiFuncForFindArgExtCoef, wDummyFindArgExtCoef

wave gwBaseFindArgExtCoef, gwRelFindArgExtCoef 
killwaves/Z wDummyFindArgExtCoef, gwBaseFindArgExtCoef, gwRelFindArgExtCoef
return V_minloc
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates "Argument Extension Coefficient" (a, see below) of one curve (wRel) with respect to another (wBase):
//wRel(x) = wBase(x*(1 + a))
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiFuncForFindArgExtCoef(wDummy, xExtCoef)
wave wDummy
variable xExtCoef

//string nmWaveBase = w[0]
//string nmWaveRel = w[1]
//wave wBase = $nmWaveBase, wRel = $nmWaveRel
wave gwBaseFindArgExtCoef, gwRelFindArgExtCoef //global waves, expected to be created in srwUtiFindArgExtCoef

variable np = dimsize(gwBaseFindArgExtCoef, 0)
variable arg = dimoffset(gwRelFindArgExtCoef, 0), argStep = dimdelta(gwRelFindArgExtCoef, 0)
variable i, s = 0, auxBuf

for(i = 0; i < np; i += 1)
	auxBuf = gwBaseFindArgExtCoef(arg*(1 + xExtCoef)) - gwRelFindArgExtCoef(arg)
	s += auxBuf*auxBuf
	arg += argStep
endfor
return sqrt(s/np)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates index of an argument wave (used at re-arranging wave: from non-uniform to uniform scaling)
//Usage: MultiLayerRefUnifMesh = Reflectivity(srwUtiCalcUniformArgInd(x, PhotonEnergy))
//where Reflectivity and PhotonEnergy are waves (function and argument)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiCalcUniformArgInd(x, wArg)
variable x
wave wArg //argument wave, assumed to be monotone

variable np = dimsize(wArg, 0)
variable i, i0
for(i=0; i<(np-1); i+=1)
	if((wArg[i] <= x) %& (x < wArg[i+1]))
		i0 = i
		break
	endif
endfor
variable fractPart = (x - wArg[i0])/(wArg[i0 + 1] - wArg[i0])
return round(i0 + 1 + fractPart)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//proc SrwUtiCalcFitPolyCoef2D(nmCoefs, nmTableOrig, nmArg1, nmArg2, nmFunc1, nmFunc2, numCoefFunc1, numCoefFunc2, nmWcoef01, nmWcoef02)
proc SrwUtiCalcFitPolyCoef2D(nmCoefs, nmTableOrig, nmArg1, nmArg2, nmFunc1, nmFunc2, nmWeights01, nmWeights02, nmWcoef01, nmWcoef02)
string nmCoefs
string nmTableOrig
string nmArg1
string nmArg2
string nmFunc1
string nmFunc2
//variable numCoefFunc1
//variable numCoefFunc2
string nmWeights01
string nmWeights02
string nmWcoef01
string nmWcoef02

variable numCoefFunc1 = dimsize($nmWcoef01, 0)
variable numCoefFunc2 = dimsize($nmWcoef02, 0)

variable npArg1 = dimsize($nmTableOrig, 0)
variable npArg2 = dimsize($nmTableOrig, 1)

make/O/N=(numCoefFunc1, numCoefFunc2) $nmCoefs
make/O/N=(numCoefFunc1) wAuxCoefs1
make/O/N=(numCoefFunc2) wAuxCoefs2
make/O/N=(npArg1) wAuxColTable
make/O/N=(npArg2) wAuxRowTable

make/O/N=(numCoefFunc1, npArg2) wAuxDemiCoefTable
string strToExe
duplicate/O $nmWcoef01 W_coef
variable i = 0
do
	wAuxColTable  = $nmTableOrig[p][i]
	//duplicate/O $nmWcoef01 W_coef
	
	//strToExe = "FuncFit/Q /H=\"00000\" " + nmFunc1 + " W_coef  wAuxColTable /X=" + nmArg1 + " /D" 
	strToExe = "FuncFit/Q " + nmFunc1 + " W_coef  wAuxColTable /X=" + nmArg1 + " /W=" + nmWeights01 + " /I=1 /D" 
		print strToExe
	execute strToExe
	
	wAuxDemiCoefTable[][i] = W_coef[p]
	
	i += 1
while(i < npArg2)

duplicate/O $nmWcoef02 W_coef
i = 0
do
	wAuxRowTable  = wAuxDemiCoefTable[i][p]
	duplicate/O $nmWcoef02 W_coef
	
	//strToExe = "FuncFit/Q " + nmFunc2 + " W_coef  wAuxRowTable /X=" + nmArg2 + " /D" 
	strToExe = "FuncFit/Q " + nmFunc2 + " W_coef  wAuxRowTable /X=" + nmArg2 + " /W=" + nmWeights02 + " /I=1 /D" 
		print strToExe
	execute strToExe
	
	$nmCoefs[i][] = W_coef[q]
		//DoUpdate
	
	i += 1
while(i < numCoefFunc1)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//CSR Form-Factor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiCSRFormFact(e_eV, sig_s_mm, Ne)
variable e_eV, sig_s_mm, Ne

variable wavelength_mm = (1.239854e-3)/e_eV
variable arg = 2*Pi*sig_s_mm/wavelength_mm
return Ne*exp(-arg*arg)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates Stokes parameters from Complex Electric Field
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiE2Stokes(EX, EZ, StokesNo)
variable/C EX, EZ
variable StokesNo

variable ReEX = real(EX), ImEX = imag(EX)
variable ReEZ = real(EZ), ImEZ = imag(EZ)

if(StokesNo == 0)
	return ReEX*ReEX + ImEX*ImEX + ReEZ*ReEZ + ImEZ*ImEZ
endif
if(StokesNo == 1)
	return ReEX*ReEX + ImEX*ImEX - ReEZ*ReEZ - ImEZ*ImEZ
endif
if(StokesNo == 2)
	return -2*(ReEX*ReEZ + ImEX*ImEZ)
endif
if(StokesNo == 3)
	return 2*(-ReEX*ImEZ + ImEX*ReEZ)
endif
return 0
//from srradinc.h:
//double LinHor = E.EwX_Re*E.EwX_Re + E.EwX_Im*E.EwX_Im;
//double LinVer = E.EwZ_Re*E.EwZ_Re + E.EwZ_Im*E.EwZ_Im;
//Stokes.s0 = LinHor + LinVer;
//Stokes.s1 = LinHor - LinVer;
//Stokes.s2 = -2.*(E.EwX_Re*E.EwZ_Re + E.EwX_Im*E.EwZ_Im);
//Stokes.s3 = 2.*(-E.EwX_Re*E.EwZ_Im + E.EwX_Im*E.EwZ_Re);
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Convert wave from Real to Complex and reverse
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiWaveDuplTypeChange(nmOutWave, nmInWave, outType)
string nmInWave = srwUtiGetValS("nmInWave", "", "SrwUtiWaveDuplTypeChange")
variable outType = srwUtiGetValN("outType", 3, "SrwUtiWaveDuplTypeChange")
string nmOutWave = srwUtiGetValS("nmOutWave", "newWave", "SrwUtiWaveDuplTypeChange")
prompt nmInWave, "Existing Wave", popup Wavelist("*",";","")
prompt outType, "Change Wave Type to:", popup "Real 32 bit;Real 64 bit;Complex 32 bit;Complex 64 bit;Text"
prompt nmOutWave, "Name for New Wave"
silent 1         |       Duplicating wave ...

if(!exists(nmInWave))
	abort "Input wave was not found"
endif
variable lenNmOutWave = strlen(nmOutWave)
if(lenNmOutWave <= 0)
	abort "No new wave name has been provided"
else
	if(lenNmOutWave > 31)
		abort "New wave name is too long"
	endif
endif

srwUtiSetValS("nmOutWave", nmOutWave, "SrwUtiWaveDuplTypeChange")
srwUtiSetValS("nmInWave", nmInWave, "SrwUtiWaveDuplTypeChange")
srwUtiSetValN("outType", outType, "SrwUtiWaveDuplTypeChange")

variable n0 = dimsize($nmInWave, 0), start0 = dimoffset($nmInWave, 0), step0 = dimdelta($nmInWave, 0)
variable n1 = dimsize($nmInWave, 1), start1 = dimoffset($nmInWave, 1), step1 = dimdelta($nmInWave, 1)
variable n2 = dimsize($nmInWave, 2), start2 = dimoffset($nmInWave, 2), step2 = dimdelta($nmInWave, 2)
variable n3 = dimsize($nmInWave, 3), start3 = dimoffset($nmInWave, 3), step3 = dimdelta($nmInWave, 3)
string units0 = WaveUnits($nmInWave,0)
string units1 = WaveUnits($nmInWave,1)
string units2 = WaveUnits($nmInWave,2)
string units3 = WaveUnits($nmInWave,3)
string unitsData = WaveUnits($nmInWave,-1)

variable oldType = WaveType($nmInWave)
variable oldIsText = 0
if(oldType == 0)
	oldIsText = 1
endif
variable oldIsComplex = oldType & 0x01

if(outType == 1)
	make/O/N=(n0,n1,n2,n3) $nmOutWave
endif
if(outType == 2)
	make/O/D/N=(n0,n1,n2,n3) $nmOutWave
endif
if(outType == 3)
	make/O/C/N=(n0,n1,n2,n3) $nmOutWave
endif
if(outType == 4)
	make/O/C/D/N=(n0,n1,n2,n3) $nmOutWave
endif
if(outType == 5)
	make/O/T/N=(n0,n1,n2,n3) $nmOutWave
endif

SetScale/P x start0,step0,units0, $nmOutWave
SetScale/P y start1,step1,units1, $nmOutWave
SetScale/P z start2,step2,units2, $nmOutWave
SetScale/P t start3,step3,units3, $nmOutWave
SetScale d 0,0,unitsData, $nmOutWave

if((outType == 1) %| (outType == 2)) //to real
	if(oldIsText)
		$nmOutWave = real(str2num($nmInWave[p][q][r][s]))
	else
		$nmOutWave = real($nmInWave[p][q][r][s])	
	endif
endif
if((outType == 3) %| (outType == 4)) //to complex
	if(oldIsComplex)
		$nmOutWave = $nmInWave[p][q][r][s]
	else
		if(oldIsText)
			$nmOutWave = cmplx(str2num($nmInWave[p][q][r][s]), 0)
		else
			$nmOutWave = cmplx($nmInWave[p][q][r][s], 0)
		endif
	endif
endif
if(outType == 5) //to text
	if(oldIsText)
		$nmOutWave = $nmInWave[p][q][r][s]
	else
		if(oldIsComplex)
			print "Copying of Complex data to Text wave was not performed"
		else
			$nmOutWave = str2num($nmInWave[p][q][r][s])
		endif
	endif
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Get 1D section (cut) of a MD wave
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiWaveGetSect(nmOutWave, nmInWave, dimNo, x0, y0, z0, t0, disp)
string nmInWave = srwUtiGetValS("nmInWave", "", "SrwUtiWaveGetSect")
variable dimNo = srwUtiGetValN("dimNo", 1, "SrwUtiWaveGetSect")
variable x0 = srwUtiGetValN("x0", 0, "SrwUtiWaveGetSect")
variable y0 = srwUtiGetValN("y0", 0, "SrwUtiWaveGetSect")
variable z0 = srwUtiGetValN("z0", 0, "SrwUtiWaveGetSect")
variable t0 = srwUtiGetValN("t0", 0, "SrwUtiWaveGetSect")
string nmOutWave = srwUtiGetValS("nmOutWave", "newWave", "SrwUtiWaveGetSect")
variable disp = srwUtiGetValN("disp", 2, "SrwUtiWaveGetSect")
prompt nmInWave, "Existing Wave", popup Wavelist("*",";","")
prompt dimNo, "Dimension over which to make Cut:", popup "x;y;z;t"
prompt x0, "First Coordinate Value"
prompt y0, "Second Coordinate Value"
prompt z0, "Third Coordinate Value"
prompt t0, "Fourth Coordinate Value"
prompt nmOutWave, "Name for New Wave"
prompt disp, "Display Extracted Wave?", popup "No;Yes"
silent 1         |       Extracting data ...

if(!exists(nmInWave))
	abort "Input wave was not found"
endif
variable lenNmOutWave = strlen(nmOutWave)
if(lenNmOutWave <= 0)
	abort "No new wave name has been provided"
else
	if(lenNmOutWave > 31)
		abort "New wave name is too long"
	endif
endif

srwUtiSetValS("nmInWave", nmInWave, "SrwUtiWaveGetSect")
srwUtiSetValN("dimNo", dimNo, "SrwUtiWaveGetSect")
srwUtiSetValN("x0", x0, "SrwUtiWaveGetSect")
srwUtiSetValN("y0", y0, "SrwUtiWaveGetSect")
srwUtiSetValN("z0", z0, "SrwUtiWaveGetSect")
srwUtiSetValN("t0", t0, "SrwUtiWaveGetSect")
srwUtiSetValS("nmOutWave", nmOutWave, "SrwUtiWaveGetSect")
srwUtiSetValN("disp", disp, "SrwUtiWaveGetSect")

variable n0 = dimsize($nmInWave, 0), start0 = dimoffset($nmInWave, 0), step0 = dimdelta($nmInWave, 0)
string units0 = WaveUnits($nmInWave,0)
variable n1 = dimsize($nmInWave, 1), start1 = dimoffset($nmInWave, 1), step1 = dimdelta($nmInWave, 1)
string units1 = WaveUnits($nmInWave,1)
variable n2 = dimsize($nmInWave, 2), start2 = dimoffset($nmInWave, 2), step2 = dimdelta($nmInWave, 2)
string units2 = WaveUnits($nmInWave,2)
variable n3 = dimsize($nmInWave, 3), start3 = dimoffset($nmInWave, 3), step3 = dimdelta($nmInWave, 3)
string units3 = WaveUnits($nmInWave,3)
string unitsData = WaveUnits($nmInWave,-1)

if(dimNo == 1)
	make/O/N=(n0) $nmOutWave
	SetScale/P x start0,step0,units0, $nmOutWave
	$nmOutWave = $nmInWave(x)(y0)(z0)(t0)
endif
if(dimNo == 2)
	make/O/N=(n1) $nmOutWave
	SetScale/P x start1,step1,units1, $nmOutWave
	$nmOutWave = $nmInWave(x0)(x)(z0)(t0)
endif
if(dimNo == 3)
	make/O/N=(n2) $nmOutWave
	SetScale/P x start2,step2,units2, $nmOutWave
	$nmOutWave = $nmInWave(x0)(y0)(x)(t0)
endif
if(dimNo == 4)
	make/O/N=(n3) $nmOutWave
	SetScale/P x start3,step3,units3, $nmOutWave
	$nmOutWave = $nmInWave(x0)(y0)(z0)(x)
endif
SetScale d 0,0,unitsData, $nmOutWave

if(disp == 2)
	display $nmOutWave
	SrwUtiGraphAddFrameAndGrid()
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculate RMS of a Radial distribution
//Returns cmplx(RMS, TotalFlux)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function/C srwUtiRadialDistrRMS(wRadialDistr, ord, fluxPortion)
wave wRadialDistr
variable fluxPortion
variable ord //if ord==0: returns 0-order moment; if ord != 0: returns RMS

duplicate/O wRadialDistr wAuxUtiRadialDistrRMS

wAuxUtiRadialDistrRMS *= x
integrate/T wAuxUtiRadialDistrRMS
variable nR = dimsize(wAuxUtiRadialDistrRMS, 0)
variable totRelFlux = wAuxUtiRadialDistrRMS[nR-1], cutRelFlux, rMax
variable momRes = 2*Pi*totRelFlux
variable momRes0 = momRes

if(ord != 0) //2nd order moment
	cutRelFlux = fluxPortion*totRelFlux
	FindLevel/Q wAuxUtiRadialDistrRMS, cutRelFlux
	rMax = V_LevelX

	wAuxUtiRadialDistrRMS = x*x*x*wRadialDistr(x)

	integrate/T wAuxUtiRadialDistrRMS
	momRes = sqrt(wAuxUtiRadialDistrRMS(rMax)/cutRelFlux)
endif

killwaves/Z wAuxUtiRadialDistrRMS
return cmplx(momRes, momRes0)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculate RMS of a Radial distribution, using Gaussian appriximation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiRadialDistrGausFitRMS(wRadialDistr)
wave wRadialDistr

variable nR = dimsize(wRadialDistr, 0)
variable stepR = dimdelta(wRadialDistr, 0)
variable rangeR = stepR*(nR - 1)
make/O/N=(2*nR - 1) wAuxUtiRadialDistrGausRMS
SetScale/P x -rangeR,stepR,"", wAuxUtiRadialDistrGausRMS

wAuxUtiRadialDistrGausRMS = 0
wAuxUtiRadialDistrGausRMS = wRadialDistr[nR - 1 - p]*srwUtiNonZeroIntervB(p, 0, nR - 1)
wAuxUtiRadialDistrGausRMS += wRadialDistr[p - nR + 1]*srwUtiNonZeroIntervB(p, nR, 2*nR - 2)
wave W_coef
K0 = 0; K2 = 0;
//CurveFit/Q/W=0/H="1010"/NTHR=0 gauss wAuxUtiRadialDistrGausRMS
//The option /NTHR=0 doesn't seem to be supported in IGOR versions < 6
CurveFit/Q/W=0/H="1010" gauss wAuxUtiRadialDistrGausRMS

killwaves/Z wAuxUtiRadialDistrGausRMS
return abs(W_coef[3]) //SigmaR = SigmaX*sqrt(2)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Appends names of waves representing traces in a graph to a text wave
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiListGraphTracesWaves(nmTextWave, nmGraph)
string nmTextWave = srwUtiGetValS("nmTextWave", "", "SrwUtiListGraphTracesWaves")
string nmGraph = srwUtiGetValS("nmGraph", "", "SrwUtiListGraphTracesWaves")
prompt nmTextWave, "Name of text wave for trace wave names"
prompt nmGraph, "Name of graph", popup " ;" + WinList("*", ";", "WIN:1")
Silent 1						|	...
PauseUpdate
srwUtiSetValS("nmTextWave", nmTextWave, "SrwUtiListGraphTracesWaves")
srwUtiSetValS("nmGraph", nmGraph, "SrwUtiListGraphTracesWaves")

if(strlen(nmTextWave) <= 0)
	abort "Please specify name for the text wave to keep names of numeric waves"
endif
if(cmpstr(nmGraph, " ") == 0)
	nmGraph = ""
endif
if(exists(nmTextWave) != 1)
	make/O/T/N=(0, 2) $nmTextWave
endif

string strTraceList = TraceNameList(nmGraph, ";", 1)
string curTraceName, curTraceWaveNames, curWaveNameY, curWaveNameX
variable i = 0, iRow = dimsize($nmTextWave, 0) - 1, alreadyExists
do
	curTraceName = StringFromList(i, strTraceList)
	if(strlen(curTraceName) == 0)
		break
	endif
	curTraceWaveNames = SrwUtiGraphTraceToWaveName(nmGraph, curTraceName)
	alreadyExists = 0
	if(strlen(curTraceWaveNames) > 0)
		curWaveNameY = StringFromList(0, curTraceWaveNames)
		if(strlen(curWaveNameY) > 0)
			alreadyExists = SrwUtiCheckIfStringExistsInCol(nmTextWave, curWaveNameY, 0)
			if(alreadyExists == 0)
				redimension/N=(dimsize($nmTextWave, 0) + 1, 2) $nmTextWave
				iRow += 1
				$nmTextWave[iRow][0] = curWaveNameY
			endif
		endif
		if(alreadyExists == 0)
			curWaveNameX = StringFromList(1, curTraceWaveNames)
			if(strlen(curWaveNameX) > 0)
				$nmTextWave[iRow][1] = curWaveNameX
			endif
		endif
	endif
	i += 1
while(1)	// exit is via break statement
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary function; checks if a string exists in a column of a text wave
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SrwUtiCheckIfStringExistsInCol(nmTextWave, strToFind, iCol)
string nmTextWave, strToFind
variable iCol
if(exists(nmTextWave) != 1)
	return 0
endif
wave/T wText = $nmTextWave
variable numRows = dimsize(wText, 0)
variable numCols = dimsize(wText, 1)
if((iCol < 0) %| (iCol >= numCols))
	return 0
endif
variable i
for(i = 0; i < numRows; i += 1)
	if(cmpstr(strToFind, wText[i][iCol]) == 0)
		return 1
	endif
endfor
return 0
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary function; used by SrwUtiListGraphTracesWaves
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function/T SrwUtiGraphTraceToWaveName(nmGraph, strTraceName)
string nmGraph, strTraceName
string resStr = ""
wave wTraceY = TraceNameToWaveRef(nmGraph, strTraceName)
if(WaveExists(wTraceY))
	resStr += NameOfWave(wTraceY)
endif
wave wTraceX = XWaveRefFromTrace(nmGraph, strTraceName)
if(WaveExists(wTraceX))
	resStr += ";" + NameOfWave(wTraceX)
endif
return resStr
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Searches for minimum or for maximum value among list of waves at a point
//Used to quickly analize multiple undulator flux/brightness tuning curves
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiFindWavesExtremAtPoint(nmResWave, nmWaveNames, x0, min_or_max)
string nmResWave = srwUtiGetValS("nmResWave", "", "SrwUtiFindWavesExtremAtPoint")
string nmWaveNames = srwUtiGetValS("nmWaveNames", "", "SrwUtiFindWavesExtremAtPoint")
variable x0 = srwUtiGetValN("x0", 1000, "SrwUtiFindWavesExtremAtPoint")
variable min_or_max = srwUtiGetValN("min_or_max", 1, "SrwUtiFindWavesExtremAtPoint")
prompt nmResWave, "Name for resulting text wave to store extrema information"
prompt nmWaveNames, "Name of text wave storing names of numeric waves to explore", popup Wavelist("*",";","TEXT:1,DIMS:2")
prompt x0, "Argument at which waves extremum should be found"
prompt min_or_max, "Look for:", popup "Minimum;Maximum"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmResWave", nmResWave, "SrwUtiFindWavesExtremAtPoint")
srwUtiSetValS("nmWaveNames", nmWaveNames, "SrwUtiFindWavesExtremAtPoint")
srwUtiSetValN("x0", x0, "SrwUtiFindWavesExtremAtPoint")
srwUtiSetValN("min_or_max", min_or_max, "SrwUtiFindWavesExtremAtPoint")

if(strlen(nmResWave) <= 0)
	abort "Please specify name for the resulting text wave to keep information on extrema"
endif
if(exists(nmWaveNames) != 1)
	abort "Can't find text wave with names of numeric waves to explore"
endif

variable numWaves = dimsize($nmWaveNames, 0)
string nmAuxVals = "wAuxValsUtiFindWavesExtrem"
string nmAuxInf = "wAuxInfUtiFindWavesExtrem"

make/O/N=(numWaves) $nmAuxVals
make/O/T/N=(numWaves) $nmAuxInf
$nmAuxVals = 1e+50
if(min_or_max == 2)
	$nmAuxVals = -1e+50
endif

string curValWaveName, curArgWaveName
variable i = 0, xStartCur, xEndCur
do
	curValWaveName = $nmWaveNames[i][0]
	curArgWaveName = $nmWaveNames[i][1]
	
	if(exists(curValWaveName) == 1)
		if(exists(curArgWaveName) == 1)
			FindLevel/Q/P $curArgWaveName, x0
			if((V_flag == 0) %& (V_LevelX >= 0) %& (V_LevelX < dimsize($curArgWaveName, 0)))
				$nmAuxVals[i] = $curValWaveName[V_LevelX]
			endif
		else
			xStartCur = dimoffset($curValWaveName, 0)
			xEndCur = xStartCur + (dimsize($curValWaveName, 0) - 1)*dimdelta($curValWaveName, 0)
			if((x0 >= xStartCur) %& (x0 <= xEndCur))
				$nmAuxVals[i] = $curValWaveName(x0)
			endif
		endif
		$nmAuxInf[i] = curValWaveName
	endif
	i += 1
while(i < numWaves)

if(min_or_max == 1)
	sort $nmAuxVals, $nmAuxVals, $nmAuxInf
endif
if(min_or_max == 2)
	sort/R $nmAuxVals, $nmAuxVals, $nmAuxInf
endif
make/O/T/N=(numWaves, 2) $nmResWave
$nmResWave[][0] = num2str($nmAuxVals[p])
$nmResWave[][1] = $nmAuxInf[p]

killwaves/Z $nmAuxVals, $nmAuxInf
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Post-processes undulator flux/brightness curves
//Used to find "envelope" of undulator flux/brightness tuning curves
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndTuningCurvesPostProc(nmResWaveCore, nmWaveNames, meth, truncMaxMinRat, disp)
string nmResWaveCore = srwUtiGetValS("nmResWaveCore", "", "SrwUtiUndTuningCurvesPostProc")
string nmWaveNames = srwUtiGetValS("nmWaveNames", "", "SrwUtiUndTuningCurvesPostProc")
variable meth = srwUtiGetValN("meth", 1, "SrwUtiUndTuningCurvesPostProc")
variable truncMaxMinRat = srwUtiGetValN("truncMaxMinRat", 100, "SrwUtiUndTuningCurvesPostProc")
variable disp = srwUtiGetValN("disp", 2, "SrwUtiUndTuningCurvesPostProc")
prompt nmResWaveCore, "Name core for resulting curve wave(s)"
prompt nmWaveNames, "Name of text wave storing names of numeric waves to explore", popup Wavelist("*",";","TEXT:1,DIMS:2")
prompt meth, "Post-processing method", popup "Just truncate individual curves;Make separate envelopes for non-overlapping curves;Make one envelope for all curves"
prompt truncMaxMinRat, "Maximal to minimal value ratio for curve truncation"
prompt disp, "Display new curves?", popup "No;Yes"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmResWaveCore", nmResWaveCore, "SrwUtiUndTuningCurvesPostProc")
srwUtiSetValS("nmWaveNames", nmWaveNames, "SrwUtiUndTuningCurvesPostProc")
srwUtiSetValN("meth", meth, "SrwUtiUndTuningCurvesPostProc")
srwUtiSetValN("truncMaxMinRat", truncMaxMinRat, "SrwUtiUndTuningCurvesPostProc")
srwUtiSetValN("disp", disp, "SrwUtiUndTuningCurvesPostProc")

if(strlen(nmResWaveCore) <= 0)
	abort "Please specify name core for resulting curve wave(s)"
endif
if(exists(nmWaveNames) != 1)
	abort "Can't find text wave with names of numeric waves to explore"
endif

variable numOrigCurves = dimsize($nmWaveNames, 0)
variable i, iCol, j, k, xCur, yCur, numEnvelopes, xStart, xStep, xEnd, xStartEnv, xEndEnv, xTestEnv, yTestEnv, iTestEnv, isOutside, makeNewEnv
variable xPrev, xNext, xAux, yAux, yAuxTestEnv, kEnv, pShift, pShift2, minAccept
variable nCurve,  nEnv, countLeft, countRight, countLeftRight, prevNpEnv
string nmWaveCur, nmWaveArgCur, nmWaveNew, nmWaveArgNew, nmWaveEnv, nmWaveEnvArg, charPref

if(meth == 1) //Just truncate individual curves
	make/O/T/N=(numOrigCurves, 2) $nmResWaveCore
	i = 0
	do
		iCol = 0
		do
			nmWaveCur = $nmWaveNames[i][iCol]
			charPref = "T0"
			if(iCol == 1)
				charPref = "A0"
			endif
			if(exists(nmWaveCur) == 1)
				nmWaveNew = srwUtiEnsureShortName(charPref + nmWaveCur)
				duplicate/O $nmWaveCur $nmWaveNew
				$nmResWaveCore[i][iCol] = nmWaveNew
			else
				if(iCol == 1)
					nmWaveArgCur = srwUtiEnsureShortName(charPref + nmWaveCur)
					nmWaveCur = nmWaveNames[i][0]
					if(exists(nmWaveCur) == 1)
						duplicate/O $nmWaveCur $nmWaveArgCur
						xStart = dimoffset($nmWaveCur, 0)
						xStep = dimdelta($nmWaveCur, 0)
						$nmWaveArgCur = xStart + p*xStep
						nmWaveNames[i][1] = nmWaveArgCur
					endif
				endif
			endif
			iCol += 1
		while(iCol < 2)
		i += 1
	while(i < numOrigCurves)
	numEnvelopes = numOrigCurves
endif
if((meth == 2) %| (meth == 3)) //Make separate envelopes for non-overlapping curves
	//find intervals for overlapping curves
	
	make/O/T/N=(0, 2) $nmResWaveCore
	numEnvelopes = 0
	i = 0
	do
		makeNewEnv = 1
		nmWaveCur = $nmWaveNames[i][0]
		nmWaveArgCur = $nmWaveNames[i][1]
		if(exists(nmWaveCur) == 1)
			if(dimsize($nmResWaveCore, 0) > 0)
				if(exists(nmWaveArgCur) != 1)
					nmWaveArgCur = srwUtiEnsureShortName("A" + nmWaveCur)
					duplicate/O $nmWaveCur $nmWaveArgCur
					xStart = dimoffset($nmWaveCur, 0)
					xStep = dimdelta($nmWaveCur, 0)
					$nmWaveArgCur = xStart + p*xStep
				endif
				xStart = $nmWaveArgCur[0]
				nCurve = dimsize($nmWaveArgCur, 0)
				xEnd = $nmWaveArgCur[nCurve - 1]
				j = 0
				do
					nmWaveEnv = $nmResWaveCore[j][0]
					nmWaveEnvArg = $nmResWaveCore[j][1]
					xStartEnv = $nmWaveEnvArg[0]
					nEnv = dimsize($nmWaveEnvArg, 0)
					xEndEnv = $nmWaveEnvArg[nEnv - 1]
					
					isOutside = 0
					if(xStartEnv <= xEndEnv)
						if((xStart < xStartEnv) %& (xEnd < xStartEnv))
							isOutside = 1
						endif
						if((xStart > xEndEnv) %& (xEnd > xEndEnv))
							isOutside = 1
						endif
					else
						if((xStart < xEndEnv) %& (xEnd < xEndEnv))
							isOutside = 1
						endif
						if((xStart > xStartEnv) %& (xEnd > xStartEnv))
							isOutside = 1
						endif
					endif
					
					if(!isOutside)
						makeNewEnv = 0
						make/O/N=(nCurve) wAuxLeftArgTuningCurvePP, wAuxLeftValTuningCurvePP, wAuxRightArgTuningCurvePP, wAuxRightValTuningCurvePP
						countLeft = 0; countRight = 0
						k = 0
						do
							xCur = $nmWaveArgCur[k]
							yCur = $nmWaveCur[k]
															
							if(xCur < xStartEnv)
								wAuxLeftArgTuningCurvePP[countLeft] = xCur
								wAuxLeftValTuningCurvePP[countLeft] = yCur
								countLeft += 1
							endif
							if(xCur > xEndEnv)
								wAuxRightArgTuningCurvePP[countRight] = xCur
								wAuxRightValTuningCurvePP[countRight] = yCur
								countRight += 1
							endif
							
							if((xCur >= xStartEnv) %& (xCur <= xEndEnv))
								xPrev = xCur
								if(k > 0)
									xPrev = $nmWaveArgCur[k - 1]
								endif
								xNext = xCur
								if(k < (nCurve - 1))
									xNext = $nmWaveArgCur[k + 1]
								endif
								if(xNext < xPrev)
									xAux = xPrev
									xPrev = xNext
									xNext = xAux
								endif
							
								FindLevel/P/Q $nmWaveEnvArg, xCur
								iTestEnv = round(V_LevelX)
								xTestEnv = $nmWaveEnvArg[iTestEnv]
								yTestEnv = $nmWaveEnv[iTestEnv]
							
								if((xStart <= xTestEnv) && (xTestEnv <= xEnd)) //OC17012019
									FindLevel/P/Q $nmWaveArgCur, xTestEnv
									yCur = $nmWaveCur[V_LevelX]
									if(yCur > yTestEnv)
										$nmWaveEnv[iTestEnv] = yCur
									endif
								endif
								
								//make small loops over neighbouring points
								if(iTestEnv > 0)
									kEnv = iTestEnv - 1									
									do
										xAux = $nmWaveEnvArg[kEnv]
										if((xAux >= xStartEnv) %& (xAux <= xEndEnv))
											if((xAux > xPrev) %& (xAux < xNext))
												yAuxTestEnv = $nmWaveEnv[kEnv]
										
												FindLevel/P/Q $nmWaveArgCur, xAux
												yAux = $nmWaveCur[V_LevelX]
												if(yAux > yAuxTestEnv)
													$nmWaveEnv[kEnv] = yAux
												endif
											else
												break
											endif
										else
											break
										endif
										kEnv -= 1
									while(kEnv >= 0)
								endif
								if(iTestEnv < (nEnv - 1))
									kEnv = iTestEnv + 1
									do
										xAux = $nmWaveEnvArg[kEnv]
										if((xAux >= xStartEnv) %& (xAux <= xEndEnv))
											if((xAux > xPrev) %& (xAux < xNext))
												yAuxTestEnv = $nmWaveEnv[kEnv]
										
												FindLevel/P/Q $nmWaveArgCur, xAux
												yAux = $nmWaveCur[V_LevelX]
												if(yAux > yAuxTestEnv)
													$nmWaveEnv[kEnv] = yAux
												endif
											else
												break
											endif
										else
											break
										endif
										kEnv += 1
									while(kEnv < nEnv)
								endif
							endif
							k += 1
						while(k < nCurve)
						countLeftRight = countLeft + countRight
						if(countLeftRight > 0)
							prevNpEnv = dimsize($nmWaveEnv, 0)
							make/O/N=(prevNpEnv + countLeftRight) wAuxNewEnvTuningCurvePP, wAuxNewEnvArgTuningCurvePP
							wAuxNewEnvTuningCurvePP = 0; wAuxNewEnvArgTuningCurvePP = 0
							if(countLeft > 0)
								pShift = countLeft - 1 + 0.001
								wAuxNewEnvArgTuningCurvePP = wAuxLeftArgTuningCurvePP[p]*srwUtiStep(pShift - p)
								wAuxNewEnvTuningCurvePP = wAuxLeftValTuningCurvePP[p]*srwUtiStep(pShift - p)
							endif
							pShift = countLeft - 0.001; pShift2 = countLeft + prevNpEnv - 1 + 0.001
							wAuxNewEnvArgTuningCurvePP += $nmWaveEnvArg[p - countLeft]*srwUtiNonZeroIntervB(p, pShift, pShift2)
							wAuxNewEnvTuningCurvePP += $nmWaveEnv[p - countLeft]*srwUtiNonZeroIntervB(p, pShift, pShift2)
							if(countRight > 0)
								pShift = countLeft + prevNpEnv; pShift2 = countLeft + prevNpEnv - 0.001
								wAuxNewEnvArgTuningCurvePP += wAuxRightArgTuningCurvePP[p - pShift]*srwUtiStep(p - pShift2)
								wAuxNewEnvTuningCurvePP += wAuxRightValTuningCurvePP[p - pShift]*srwUtiStep(p - pShift2)
							endif
							
							duplicate/O wAuxNewEnvArgTuningCurvePP $nmWaveEnvArg
							duplicate/O wAuxNewEnvTuningCurvePP $nmWaveEnv
							killwaves/Z wAuxNewEnvTuningCurvePP, wAuxNewEnvArgTuningCurvePP
						endif
						break
					endif
					j += 1
				while(j < numEnvelopes)
			endif

			if(makeNewEnv)
				redimension/N=(numEnvelopes + 1, 2) $nmResWaveCore
				nmWaveNew = srwUtiEnsureShortName("T" + nmWaveCur)
				duplicate/O $nmWaveCur $nmWaveNew
				$nmResWaveCore[numEnvelopes][0] = nmWaveNew
				if(exists(nmWaveArgCur) == 1)
					nmWaveArgNew = srwUtiEnsureShortName("T" + nmWaveArgCur)
					duplicate/O $nmWaveArgCur $nmWaveArgNew
				else
					nmWaveArgNew = srwUtiEnsureShortName("A" + nmWaveCur)
					duplicate/O $nmWaveCur $nmWaveArgNew
					xStart = dimoffset($nmWaveCur, 0)
					xStep = dimdelta($nmWaveCur, 0)
					$nmWaveArgNew = xStart + p*xStep
				endif
				$nmResWaveCore[numEnvelopes][1] = nmWaveArgNew
				numEnvelopes += 1
			endif
		endif
		i += 1
	while(i < numOrigCurves)
endif

if(meth == 3) //Make one envelope for all curves; i.e. merge all envelopes obtained at previous step
	numEnvelopes = dimsize($nmResWaveCore, 0)
	if(numEnvelopes > 1)
		nmWaveEnv = $nmResWaveCore[0][0]
		nmWaveEnvArg = $nmResWaveCore[0][1]
		xStartEnv = $nmWaveEnvArg[0]
		nEnv = dimsize($nmWaveEnvArg, 0)
		xEndEnv = $nmWaveEnvArg[nEnv - 1]
		i = 1
		do
			nmWaveCur = $nmResWaveCore[i][0]
			nmWaveArgCur = $nmResWaveCore[i][1]
			xStart = $nmWaveArgCur[0]
			nCurve = dimsize($nmWaveArgCur, 0)
			xEnd = $nmWaveArgCur[nCurve - 1]
			
			if(nCurve > 0)
				make/O/N=(dimsize($nmWaveEnvArg, 0) + nCurve) wAuxNewEnvTuningCurvePP, wAuxNewEnvArgTuningCurvePP
				wAuxNewEnvTuningCurvePP = 0; wAuxNewEnvArgTuningCurvePP = 0
				
				if(0.5*(xStart + xEnd) < 0.5*(xStartEnv + xEndEnv)) //assuming that the input envelope curves are not overlapping
					pShift = nCurve - 1 + 0.001
					wAuxNewEnvArgTuningCurvePP = $nmWaveArgCur[p]*srwUtiStep(pShift - p)
					wAuxNewEnvTuningCurvePP = $nmWaveCur[p]*srwUtiStep(pShift - p)
					pShift = nCurve - 0.001; pShift2 = nCurve + nEnv - 1 + 0.001
					wAuxNewEnvArgTuningCurvePP += $nmWaveEnvArg[p - nCurve]*srwUtiNonZeroIntervB(p, pShift, pShift2)
					wAuxNewEnvTuningCurvePP += $nmWaveEnv[p - nCurve]*srwUtiNonZeroIntervB(p, pShift, pShift2)
				else
					pShift = nEnv - 1 + 0.001
					wAuxNewEnvArgTuningCurvePP = $nmWaveEnvArg[p]*srwUtiStep(pShift - p)
					wAuxNewEnvTuningCurvePP = $nmWaveEnv[p]*srwUtiStep(pShift - p)
					pShift = nEnv - 0.001; pShift2 = nCurve + nEnv - 1 + 0.001
					wAuxNewEnvArgTuningCurvePP += $nmWaveArgCur[p - nEnv]*srwUtiNonZeroIntervB(p, pShift, pShift2)
					wAuxNewEnvTuningCurvePP += $nmWaveCur[p - nEnv]*srwUtiNonZeroIntervB(p, pShift, pShift2)
				endif
				duplicate/O wAuxNewEnvArgTuningCurvePP $nmWaveEnvArg
				duplicate/O wAuxNewEnvTuningCurvePP $nmWaveEnv
				killwaves/Z wAuxNewEnvArgTuningCurvePP, wAuxNewEnvTuningCurvePP
			endif
			i += 1
		while(i < numEnvelopes)
		redimension/N=(1, 2) $nmResWaveCore
		$nmResWaveCore[0][0] = nmWaveEnv
		$nmResWaveCore[0][1] = nmWaveEnvArg
	endif
endif

if(truncMaxMinRat >= 0) //curve truncation is necessary
	numEnvelopes = dimsize($nmResWaveCore, 0)
	if(numEnvelopes <= 0)
		return
	endif
	i = 0
	do
		nmWaveEnv = $nmResWaveCore[i][0]
		nmWaveEnvArg = $nmResWaveCore[i][1]
		wavestats/Q $nmWaveEnv
		minAccept = V_min
		
		if(V_min*V_max >= 0)
			if(V_min > 0)
				if(V_max/V_min > truncMaxMinRat)
					minAccept = V_max/truncMaxMinRat
				endif
			else
				if(V_min < 0)
					if(V_max/V_min > 1/truncMaxMinRat)
						minAccept = V_max*truncMaxMinRat
					endif
				else
					if(V_max > 0)
						minAccept = V_max/truncMaxMinRat
					endif
					if(V_max < 0)
						minAccept = V_max*truncMaxMinRat
					endif
				endif
			endif
		endif
	
		if(minAccept != V_min)
			duplicate/O $nmWaveEnv wAuxNewEnvTuningCurvePP
			duplicate/O $nmWaveEnvArg wAuxNewEnvArgTuningCurvePP
			
			nEnv = dimsize($nmWaveEnvArg, 0)
			countLeft = 0
			k = 0
			do
				yCur = $nmWaveEnv[k]
				if(yCur >= minAccept)
					wAuxNewEnvTuningCurvePP[countLeft] = yCur
					wAuxNewEnvArgTuningCurvePP[countLeft] = $nmWaveEnvArg[k]
					countLeft += 1
				endif
				k += 1
			while(k < nEnv)
			redimension/N=(countLeft) $nmWaveEnv, $nmWaveEnvArg
			$nmWaveEnv = wAuxNewEnvTuningCurvePP[p]
			$nmWaveEnvArg = wAuxNewEnvArgTuningCurvePP[p]
			killwaves/Z wAuxNewEnvTuningCurvePP, wAuxNewEnvArgTuningCurvePP
		endif
		i += 1
	while(i < numEnvelopes)
endif

variable someCurvesAppended = 0
variable maxNumResCurves = dimsize($nmResWaveCore, 0)
if(disp == 2) //Display all waves with names from $nmResWaveCore in one graph
	if(maxNumResCurves > 0)
		display
		i = 0
		do
			nmWaveEnv = $nmResWaveCore[i][0]
			nmWaveEnvArg = $nmResWaveCore[i][1]
			if(exists(nmWaveEnv) == 1)
				if(exists(nmWaveEnvArg) == 1)
					appendtograph $nmWaveEnv vs $nmWaveEnvArg
					someCurvesAppended = 1
				else
					appendtograph $nmWaveEnv
				endif
			endif
			i += 1
		while(i < maxNumResCurves)
		if(someCurvesAppended)
			SrwUtiGraphAddFrameAndGrid()
		endif
	endif
endif

killwaves/Z wAuxLeftArgTuningCurvePP, wAuxLeftValTuningCurvePP, wAuxRightArgTuningCurvePP, wAuxRightValTuningCurvePP
killwaves/Z wAuxNewEnvTuningCurvePP, wAuxNewEnvArgTuningCurvePP
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Generate coordinates of ponts randomly filling 3D volume limited by two arbitrary curves (defining base) and two surfaces
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiRandFillVolume(nmResCoord, np, xMin, xMax, nmCurveMinY, nmCurveMaxY, nmSurfMinZ, nmSurfMaxZ)
string nmResCoord = srwUtiGetValS("nmResCoord", "wRandCoord", "SrwUtiRandFillVolume")
variable np = srwUtiGetValN("np", 10, "SrwUtiRandFillVolume")
variable xMin = srwUtiGetValN("xMin", -10, "SrwUtiRandFillVolume")
variable xMax = srwUtiGetValN("xMax", 10, "SrwUtiRandFillVolume")
string nmCurveMinY = srwUtiGetValS("nmCurveMinY", "wCurveMinY", "SrwUtiRandFillVolume")
string nmCurveMaxY = srwUtiGetValS("nmCurveMaxY", "wCurveMaxY", "SrwUtiRandFillVolume")
string nmSurfMinZ = srwUtiGetValS("nmSurfMinZ", "wSurfMinZ", "SrwUtiRandFillVolume")
string nmSurfMaxZ = srwUtiGetValS("nmSurfMaxZ", "wSurfMinZ", "SrwUtiRandFillVolume")
prompt nmResCoord, "Name for res. random coord. wave"
prompt np, "Number of random points"
prompt xMin, "Min. X coordinate"
prompt xMax, "Max. X coordinate"
prompt nmCurveMinY, "Name of lower Y curve wave", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt nmCurveMaxY, "Name of upper Y curve wave", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt nmSurfMinZ, "Name of lower Z surface wave", popup Wavelist("*",";","TEXT:0,DIMS:2")
prompt nmSurfMaxZ, "Name of upper Z surface wave", popup Wavelist("*",";","TEXT:0,DIMS:2")
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmResCoord", nmResCoord, "SrwUtiRandFillVolume")
srwUtiSetValN("np", np, "SrwUtiRandFillVolume")
srwUtiSetValN("xMin", xMin, "SrwUtiRandFillVolume")
srwUtiSetValN("xMax", xMax, "SrwUtiRandFillVolume")
srwUtiSetValS("nmCurveMinY", nmCurveMinY, "SrwUtiRandFillVolume")
srwUtiSetValS("nmCurveMaxY", nmCurveMaxY, "SrwUtiRandFillVolume")
srwUtiSetValS("nmSurfMinZ", nmSurfMinZ, "SrwUtiRandFillVolume")
srwUtiSetValS("nmSurfMaxZ", nmSurfMaxZ, "SrwUtiRandFillVolume")

variable len_nmResCoord = strlen(nmResCoord)

if(len_nmResCoord <= 0)
	abort "Resulting coordinates wave name was not provided"
endif
if(len_nmResCoord > 28)
	abort "Resulting coordinates wave name is too long"
endif
if(np <= 0)
	abort "Number of points should be positive"
endif
if(xMin >= xMax)
	abort "X coordinate limits are not defined correctly"
endif

wavestats/Q $nmCurveMinY
variable yMin = V_min
wavestats/Q $nmCurveMaxY
variable yMax = V_max
if(yMin >= yMax)
	abort "Y coordinate limits are not defined correctly"
endif

wavestats/Q $nmSurfMinZ
variable zMin = V_min
wavestats/Q $nmSurfMaxZ
variable zMax = V_max
if(zMin >= zMax)
	abort "Z coordinate limits are not defined correctly"
endif

make/O/N=(np, 3) $nmResCoord

variable xCen = 0.5*(xMIn + xMax), yCen = 0.5*(yMIn + yMax), zCen = 0.5*(zMIn + zMax)
variable xHalfRange = 0.5*(xMax - xMIn), yHalfRange = 0.5*(yMax - yMIn), zHalfRange = 0.5*(zMax - zMIn)

variable i=0, iTot=0, xx, yy, zz, yTestMin, yTestMax, zTestMin, zTestMax
do
	xx = xCen + enoise(xHalfRange)
	yy = yCen + enoise(yHalfRange)
	
	yTestMin = $nmCurveMinY(xx)
	yTestMax = $nmCurveMaxY(xx)
	
	if((yy >= yTestMin) %& (yy <= yTestMax))
		zz = zCen + enoise(zHalfRange)
		
		zTestMin = $nmSurfMinZ(xx)(yy)
		zTestMax = $nmSurfMaxZ(xx)(yy)
		
		if((zz >= zTestMin) %& (zz <= zTestMax))
			$nmResCoord[i][0] = xx
			$nmResCoord[i][1] = yy
			$nmResCoord[i][2] = zz
			i += 1
		endif
	endif
	iTot += 1
while(i < np)
	print "Number of points in rect. parallelepiped:", iTot
	print "Eslimated volume:", (np/iTot)*(xMax - xMIn)*(yMax - yMIn)*(zMax - zMIn)
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Calculates geom. path in "bubbles", without keeping track of eventual "bubble" overlapping 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SrwUtiRayPathInSpheres(xx, yy, wSphCenCoord, wSphRad)
variable xx, yy
wave wSphCenCoord, wSphRad

variable n = dimsize(wSphCenCoord, 0), nRad = dimsize(wSphRad, 0)
if(n > nRad)
	n = nRad
endif
variable i, sumPath = 0, dx, dy, ue2, rr, re2

for(i = 0; i < n; i += 1)
	dx = xx - wSphCenCoord[i][0]
	dy = yy - wSphCenCoord[i][1]
	rr = wSphRad[i]
	re2 = rr*rr
	ue2 = dx*dx + dy*dy
	if(ue2 < re2)
		sumPath += 2*sqrt(re2 - ue2)
	endif
endfor
return sumPath
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// E.g.: "#5000.0 #Initial Photon Energy [eV]" -> "5000.0 ".
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function/S SrwUtiGetStrBwSep(str, sep, [sep2])
string str, sep, sep2

//print str
//print sep //, sep2
//return str

variable i1stSep = strsearch(str, sep, 0)
if(i1stSep < 0)
	return ""
endif
variable iStartNum = i1stSep + strlen(sep)

if(ParamIsDefault(sep2) == 0)
	sep = sep2
endif

variable i2ndSep = strsearch(str, sep, iStartNum)
if(i2ndSep < 0)
	i2ndSep = strlen(str)
endif
return str[iStartNum, i2ndSep - 1]
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// E.g.: "#5000.0 #Initial Photon Energy [eV]" -> 5000.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SrwUtiParseNumBwSep(str, sep, [separ2])
string str, sep, separ2

//variable i1stSep = strsearch(str, sep, 0)
//if(i1stSep < 0)
//	return 0
//endif
//variable iStartNum = i1stSep + strlen(sep)
//
//if(ParamIsDefault(sep2) == 0)
//	sep = sep2
//endif
//
//variable i2ndSep = strsearch(str, sep, iStartNum)
//if(i2ndSep < 0)
//	i2ndSep = strlen(str)
//endif
//string sNum = str[iStartNum, i2ndSep - 1]

if(ParamIsDefault(separ2))
	return str2num(SrwUtiGetStrBwSep(str, sep))
else
	return str2num(SrwUtiGetStrBwSep(str, sep, sep2=separ2))
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Import a single-column intensity data file and plot the data in a graph or as an image
//Expected header format:
//#C-aligned Intensity (inner loop is vs photon energy, outer loop vs vertical position)
//#5000.0 #Initial Photon Energy [eV]
//#5000.0 #Final Photon Energy [eV]
//#1 #Number of points vs Photon Energy
//#-4.1579932835224167e-07 #Initial Horizontal Position [m]
//#4.157993283522281e-07 #Final Horizontal Position [m]
//#260 #Number of points vs Horizontal Position
//#-5.193669768581371e-07 #Initial Vertical Position [m]
//#5.193669768581371e-07 #Final Vertical Position [m]
//#416 #Number of points vs Vertical Position
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiLoadAndPlotDataCol(nmWaveMD, xc, yc, disp)
string nmWaveMD = srwUtiGetValS("nmWaveMD", "wInt", "SrwUtiLoadAndPlotDataCol")
variable xc = srwUtiGetValN("xc", 0, "SrwUtiLoadAndPlotDataCol")
variable yc = srwUtiGetValN("yc", 0, "SrwUtiLoadAndPlotDataCol")
variable disp = srwUtiGetValN("disp", 2, "SrwUtiLoadAndPlotDataCol")
prompt nmWaveMD, "Name of wave to create"
prompt xc, "Horizontal Position for Vertical Cut"
prompt yc, "Vertical Position for Horizontal Cut"
prompt disp, "New Display?", popup "No;Yes"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmWaveMD", nmWaveMD, "SrwUtiLoadAndPlotDataCol")
srwUtiSetValN("xc", xc, "SrwUtiLoadAndPlotDataCol")
srwUtiSetValN("yc", yc, "SrwUtiLoadAndPlotDataCol")
srwUtiSetValN("disp", disp, "SrwUtiLoadAndPlotDataCol")

//To keep updated!
string sep = "#"
variable iPhEnStart = 1 //0
variable iPhEnEnd = iPhEnStart + 1
variable iPhEnN = iPhEnStart + 2
variable iHorPosStart = 4 //3
variable iHorPosEnd = iHorPosStart + 1
variable iHorPosN = iHorPosStart + 2
variable iVertPosStart = 7 //6
variable iVertPosEnd = iVertPosStart + 1
variable iVertPosN = iVertPosStart + 2
variable iStokesN = iVertPosStart + 3

//Loading Header:
string nmHeader = "wAuxHeadUtiLoadAndPlotData"
//LoadWave/Q/A/J/O/B="N=wAuxHeadUtiLoadAndPlotData;"/L={0, 1, 9, 0, 1}/K=2
//LoadWave/Q/A/J/O/B="N=wAuxHeadUtiLoadAndPlotData;"/L={0, 1, 10, 0, 1}/K=2
//LoadWave/Q/A/J/O/B="N=wAuxHeadUtiLoadAndPlotData;"/L={0, 0, 10, 0, 1}/K=2
LoadWave/Q/A/J/O/B="N=wAuxHeadUtiLoadAndPlotData;"/L={0, 0, 11, 0, 1}/K=2

string sFilePath = S_path + S_fileName

variable phEnStart = SrwUtiParseNumBwSep($nmHeader[iPhEnStart], sep)
variable phEnEnd = SrwUtiParseNumBwSep($nmHeader[iPhEnEnd], sep)
variable phEnN = SrwUtiParseNumBwSep($nmHeader[iPhEnN], sep)
variable horPosStart = SrwUtiParseNumBwSep($nmHeader[iHorPosStart], sep)
variable horPosEnd = SrwUtiParseNumBwSep($nmHeader[iHorPosEnd], sep)
variable horPosN = SrwUtiParseNumBwSep($nmHeader[iHorPosN], sep)
variable vertPosStart = SrwUtiParseNumBwSep($nmHeader[iVertPosStart], sep)
variable vertPosEnd = SrwUtiParseNumBwSep($nmHeader[iVertPosEnd], sep)
variable vertPosN = SrwUtiParseNumBwSep($nmHeader[iVertPosN], sep)

//Using units defined in the input file
string strInUnitPhEn = SrwUtiGetStrBwSep($nmHeader[iPhEnStart], "[", sep2="]")
string strInUnitHorPos = SrwUtiGetStrBwSep($nmHeader[iHorPosStart], "[", sep2="]")
string strInUnitVertPos = SrwUtiGetStrBwSep($nmHeader[iVertPosStart], "[", sep2="]")

variable nStokes = 1
string sTestStokesN = $nmHeader[iStokesN]

//print "Aha, iStokesN=", iStokesN, " sTestStokesN=", sTestStokesN
//variable it=0
//do
//	sTestStokesN = $nmHeader[it]
//	print sTestStokesN
//	it+=1
//while(it<13)

string sTestSymbStokesN = sTestStokesN[0, 0]
variable iDataStart = 10
//print sTestStokesN, sTestSymbStokesN
if(cmpstr(sTestSymbStokesN, "#") == 0)
	nStokes = SrwUtiParseNumBwSep($nmHeader[iStokesN], sep)
	iDataStart += 1
	//print "aha"
endif
//currently, this proc reads only first Stokes component (s0) -- to upgrade

variable isMutual = 0, isCmplx = 0, isDegCoh = 0
string sEnt = $nmHeader[0]
//if(cmpstr(sEnt[0, 6], "#Mutual") == 0)

string strLabelEntity = "Intensity" //OC13112019

if(cmpstr(sEnt[0, 19], "#Degree of Coherence") == 0)
	isDegCoh = 1
	strLabelEntity = "Degree of Coherence" //OC13112019
endif
//if((cmpstr(sEnt[0, 6], "#Mutual") == 0) || (cmpstr(sEnt[0, 19], "#Degree of Coherence") == 0))
if(cmpstr(sEnt[0, 6], "#Mutual") == 0)
	isMutual = 1
	strLabelEntity = "Mutual Intensity" //OC13112019
endif

if(cmpstr(sEnt[0, 14], "#Complex Mutual") == 0)
	isMutual = 1
	isCmplx = 1
	strLabelEntity = "Complex Mutual Intensity" //OC13112019
endif

//print phEnStart, phEnEnd, phEnN, horPosStart, horPosEnd, horPosN

variable nVal = phEnN*horPosN*vertPosN
//if(isMutual != 0)
if((isMutual != 0) || (isDegCoh != 0)) //OC12112019
	if(phEnN == 1)
		if(horPosN == 1)
			nVal = vertPosN*vertPosN
		else
			if(vertPosN == 1)
				nVal = horPosN*horPosN
			endif
		endif
	endif
endif
if(isCmplx != 0)
	nVal *= 2
endif

//print iDataStart

string nmFlatWave = "wAuxFlatUtiLoadAndPlotData"
//LoadWave/Q/A/J/O/B="N=wAuxFlatUtiLoadAndPlotData;"/L={0, 10, phEnN*horPosN*vertPosN, 0, 1}/K=1 sFilePath

//print iDataStart, nVal
LoadWave/Q/A/J/O/B="N=wAuxFlatUtiLoadAndPlotData;"/L={0, iDataStart, nVal, 0, 1}/K=1 sFilePath

variable numDim = 0, arg1Start = 0, arg1End = 0, arg1N = 0, arg2Start = 0, arg2End = 0, arg2N = 0, arg3Start = 0, arg3End = 0, arg3N = 0
string strLabel1, strLabel2, strUnit1, strUnit2, strUnit3
if(phEnN > 1)
	numDim += 1
	arg1Start = phEnStart
	arg1End = phEnEnd
	arg1N = phEnN
	strUnit1 = "eV"
	if(strlen(strInUnitPhEn) > 0)
		strUnit1 = strInUnitPhEn
	endif
	strLabel1 = "Photon Energy"
	strLabel2 = "Intensity"
endif
if(horPosN > 1)
	numDim += 1
	if(phEnN > 1)
		arg2Start = horPosStart
		arg2End = horPosEnd
		arg2N = horPosN
		strUnit2 = "m"
		strLabel2 = "Horizontal Position"
		if(strlen(strInUnitHorPos) > 0)
			strUnit2 = strInUnitHorPos
			if(cmpstr(strInUnitHorPos, "rad") == 0)
				strLabel2 = "Horizontal Angle"
			endif
		endif
	else
		arg1Start = horPosStart
		arg1End = horPosEnd
		arg1N = horPosN
		strUnit1 = "m"
		strLabel1 = "Horizontal Position"
		strLabel2 = "Intensity"		
		if(strlen(strInUnitHorPos) > 0)
			strUnit1 = strInUnitHorPos
			if(cmpstr(strInUnitHorPos, "rad") == 0)
				strLabel1 = "Horizontal Angle"
			endif
		endif
		if(isMutual != 0)
			numDim += 1
			arg2Start = horPosStart
			arg2End = horPosEnd
			arg2N = horPosN
			strUnit2 = "m"
			strLabel2 = "Horizontal Position (conj.)"
			if(strlen(strInUnitHorPos) > 0)
				strUnit2 = strInUnitHorPos
				if(cmpstr(strInUnitHorPos, "rad") == 0)
					strLabel1 = "Horizontal Angle (conj.)"
				endif
			endif
		endif
	endif
endif
if(vertPosN > 1)
	numDim += 1
	if(phEnN > 1)
		if(horPosN > 1)
			arg3Start = vertPosStart
			arg3End = vertPosEnd
			arg3N = vertPosN
			strUnit3 = "m"
			if(strlen(strInUnitVertPos) > 0)
				strUnit3 = strInUnitVertPos
			endif
		else
			arg2Start = vertPosStart
			arg2End = vertPosEnd
			arg2N = vertPosN
			strUnit2 = "m"
			strLabel2 = "Vertical Position"
			if(strlen(strInUnitVertPos) > 0)
				strUnit2 = strInUnitVertPos
				if(cmpstr(strInUnitVertPos, "rad") == 0)
					strLabel2 = "Vertical Angle"
				endif
			endif
		endif
	else
		if(horPosN > 1)
			arg2Start = vertPosStart
			arg2End = vertPosEnd
			arg2N = vertPosN
			strUnit2 = "m"
			strLabel2 = "Vertical Position"
			if(strlen(strInUnitVertPos) > 0)
				strUnit2 = strInUnitVertPos
				if(cmpstr(strInUnitVertPos, "rad") == 0)
					strLabel2 = "Vertical Angle"
				endif
			endif
		else
			arg1Start = vertPosStart
			arg1End = vertPosEnd
			arg1N = vertPosN
			strUnit1 = "m"
			strLabel1 = "Vertical Position"
			strLabel2 = "Intensity"			
			if(strlen(strInUnitVertPos) > 0)
				strUnit1 = strInUnitVertPos
				if(cmpstr(strInUnitVertPos, "rad") == 0)
					strLabel1 = "Vertical Angle"
				endif
			endif
			if(isMutual != 0)
				numDim += 1
				arg2Start = vertPosStart
				arg2End = vertPosEnd
				arg2N = vertPosN
				strUnit2 = "m"
				strLabel2 = "Vertical Position (conj.)"
				if(strlen(strInUnitVertPos) > 0)
					strUnit2 = strInUnitVertPos
					if(cmpstr(strInUnitVertPos, "rad") == 0)
						strLabel2 = "Vertical Angle (conj.)"
					endif
				endif
			endif
		endif	
	endif			
endif

if(isDegCoh != 0) //OC12112019
	numDim = 2
	strLabelEntity = "Degree of Coherence"
	if(horPosN > 1)
		arg1Start = horPosStart
		arg1End = horPosEnd
		arg1N = horPosN
		strUnit1 = "m"
		strLabel1 = "(x1+x2)/2"
		if(strlen(strInUnitHorPos) > 0)
			strUnit1 = strInUnitHorPos
			if(cmpstr(strInUnitHorPos, "rad") == 0)
				strLabel1 = "(x'1+x'2)/2"
			endif
		endif
		arg2Start = (horPosStart - horPosEnd)/2
		arg2End = (horPosEnd - horPosStart)/2
		arg2N = horPosN
		strUnit2 = "m"
		strLabel2 = "(x1-x2)/2"
		if(strlen(strInUnitHorPos) > 0)
			strUnit2 = strInUnitHorPos
			if(cmpstr(strInUnitHorPos, "rad") == 0)
				strLabel2 = "(x'1-x'2)/2"
			endif
		endif
	endif
	if(vertPosN > 1)	
		arg1Start = vertPosStart
		arg1End = vertPosEnd
		arg1N = vertPosN
		strUnit1 = "m"
		strLabel1 = "(y1+y2)/2"
		if(strlen(strInUnitVertPos) > 0)
			strUnit1 = strInUnitVertPos
			if(cmpstr(strInUnitVertPos, "rad") == 0)
				strLabel1 = "(y'1+y'2)/2"
			endif
		endif
		arg2Start = (vertPosStart - vertPosEnd)/2
		arg2End = (vertPosEnd - vertPosStart)/2
		arg2N = vertPosN
		strUnit2 = "m"
		strLabel2 = "(y1-y2)/2"
		if(strlen(strInUnitVertPos) > 0)
			strUnit2 = strInUnitVertPos
			if(cmpstr(strInUnitVertPos, "rad") == 0)
				strLabel2 = "(y'1-y'2)/2"
			endif
		endif
	endif
endif

string nmWaveMDx = nmWaveMD + "x"
string nmWaveMDy = nmWaveMD + "y"
//print numDim, strLabel1, strLabel2

string nmReWaveMD, nmImWaveMD, nmReWaveMDx, nmImWaveMDx, nmReWaveMDy, nmImWaveMDy
if(isCmplx != 0)
	nmReWaveMD = "re" + nmWaveMD
	nmImWaveMD = "im" + nmWaveMD
	nmReWaveMDx = "re" + nmWaveMDx
	nmImWaveMDx = "im" + nmWaveMDx
	nmReWaveMDy = "re" + nmWaveMDy
	nmImWaveMDy = "im" + nmWaveMDy
endif

if(numDim == 1)
	make/O/N=(arg1N) $nmWaveMD
	$nmWaveMD = $nmFlatWave[p]
	SetScale/I x arg1Start, arg1End, strUnit1, $nmWaveMD
	if(disp == 2)
		display $nmWaveMD
		Label bottom strLabel1
		Label left strLabel2
		SrwUtiGraphAddFrameAndGrid()
	endif
endif
if(numDim == 2)
	if(isCmplx == 0)
		make/O/N=(arg1N, arg2N) $nmWaveMD
		$nmWaveMD = $nmFlatWave[arg1N*q + p]
	else
		make/O/N=(arg1N, arg2N) $nmReWaveMD, $nmImWaveMD
		$nmReWaveMD = $nmFlatWave[(2*arg1N)*q + 2*p]
		$nmImWaveMD = $nmFlatWave[(2*arg1N)*q + 2*p + 1]
	endif
	
	if(isCmplx == 0)
		SetScale/I x arg1Start, arg1End, strUnit1, $nmWaveMD
		SetScale/I y arg2Start, arg2End, strUnit2, $nmWaveMD
		if(disp == 2)
			display; appendimage $nmWaveMD
			Label bottom strLabel1
			Label left strLabel2
		endif
	else
		SetScale/I x arg1Start, arg1End, strUnit1, $nmReWaveMD, $nmImWaveMD
		SetScale/I y arg2Start, arg2End, strUnit2, $nmReWaveMD, $nmImWaveMD
		if(disp == 2)
			display; appendimage $nmReWaveMD
			Label bottom strLabel1
			Label left strLabel2
			display; appendimage $nmImWaveMD
			Label bottom strLabel1
			Label left strLabel2
		endif
	endif

	if(isCmplx == 0)	
		make/O/N=(arg1N) $nmWaveMDx
		SetScale/I x arg1Start, arg1End, strUnit1, $nmWaveMDx
		$nmWaveMDx = interp2d($nmWaveMD,x,yc) //OC23092022
		//$nmWaveMDx = $nmWaveMD(x)(yc)
		if(disp == 2)
			display $nmWaveMDx
			Label bottom strLabel1
			if(strlen(strLabelEntity) > 0) //OC13112019
				Label left strLabelEntity
			endif
			SrwUtiGraphAddFrameAndGrid()
		endif
		make/O/N=(arg2N) $nmWaveMDy
		SetScale/I x arg2Start, arg2End, strUnit2, $nmWaveMDy
		$nmWaveMDy = interp2d($nmWaveMD,xc,x) //OC23092022
		//$nmWaveMDy = $nmWaveMD(xc)(x)
		if(disp == 2)	
			display $nmWaveMDy
			Label bottom strLabel2
			if(strlen(strLabelEntity) > 0) //OC13112019
				Label left strLabelEntity
			endif
			SrwUtiGraphAddFrameAndGrid()
		endif
	else
		make/O/N=(arg1N) $nmReWaveMDx, $nmImWaveMDx
		SetScale/I x arg1Start, arg1End, strUnit1, $nmReWaveMDx, $nmImWaveMDx
		$nmReWaveMDx = interp2d($nmReWaveMD,x,yc) //OC23092022
		$nmImWaveMDx = interp2d($nmImWaveMD,x,yc) //OC23092022
		//$nmReWaveMDx = $nmReWaveMD(x)(yc)
		//$nmImWaveMDx = $nmImWaveMD(x)(yc)
		if(disp == 2)
			display $nmReWaveMDx, $nmImWaveMDx
			Label bottom strLabel1
			if(strlen(strLabelEntity) > 0) //OC13112019
				Label left strLabelEntity
			endif
			SrwUtiGraphAddFrameAndGrid()
		endif
		make/O/N=(arg2N) $nmReWaveMDy, $nmImWaveMDy
		SetScale/I x arg2Start, arg2End, strUnit2, $nmReWaveMDy, $nmImWaveMDy
		$nmReWaveMDy = interp2d($nmReWaveMD,xc,x) //OC23092022
		$nmImWaveMDy = interp2d($nmImWaveMD,xc,x) //OC23092022
		//$nmReWaveMDy = $nmReWaveMD(xc)(x)
		//$nmImWaveMDy = $nmImWaveMD(xc)(x)
		if(disp == 2)	
			display $nmReWaveMDy, $nmImWaveMDy
			Label bottom strLabel2
			if(strlen(strLabelEntity) > 0) //OC13112019
				Label left strLabelEntity
			endif
			SrwUtiGraphAddFrameAndGrid()
		endif
	endif
endif
if(numDim == 3)
	make/O/N=(arg1N, arg2N, arg3N) $nmWaveMD
	SetScale/I x arg1Start, arg1End, strUnit1, $nmWaveMD
	SetScale/I y arg2Start, arg2End, strUnit2, $nmWaveMD
	SetScale/I z arg3Start, arg3End, strUnit3, $nmWaveMD
	$nmWaveMD = $nmFlatWave[arg2N*arg1N*r + arg1N*q + p]
	if(disp == 2)	
		display $nmWaveMD[][round((xc - arg2Start)/dimdelta($nmWaveMD,1))][round((yc - arg3Start)/dimdelta($nmWaveMD,2))]
		Label bottom "Photon Energy"; Label left "Intensity"; SrwUtiGraphAddFrameAndGrid()
		display $nmWaveMD[round(0.5*(arg1End - arg1Start)/dimdelta($nmWaveMD,0))][][round((yc - arg3Start)/dimdelta($nmWaveMD,2))]
		Label bottom "Horizontal Position"; Label left "Intensity"; SrwUtiGraphAddFrameAndGrid()
		display $nmWaveMD[round(0.5*(arg1End - arg1Start)/dimdelta($nmWaveMD,0))][round((xc - arg2Start)/dimdelta($nmWaveMD,1))][]
		Label bottom "Vertical Position"; Label left "Intensity"; SrwUtiGraphAddFrameAndGrid()
	endif
endif

killwaves/Z $nmFlatWave
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Auxiliary function to identify and get rid of NaN
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiNan2Num(x, num)
variable x, num
if(numtype(x) != 0)
	return num
else
	return x
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//FInd Transverse Coherence Length from Degree of Coherence (1D cut)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwUtiCohLenFromDegCoh(wDC)
wave wDC

FindLevels wDC, 0.5
wave W_FindLevels
variable nLev = dimsize(W_FindLevels, 0) 
variable i=0, curX, prevX, iFirstPos = -1
do 
	curX = W_FindLevels[i]
	
	if((i > 0) && (prevX < 0) && (curX > 0))
		iFirstPos = i
		break
	endif
	prevX = curX
	i += 1
while(i < nLev)
if(iFirstPos < 0)
	return -1
else
	return curX - prevX
endif

end



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Conversion from Mutial Intensity to (rotated) Degree of Coherence function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiCmplxMutInt2RotDegCoh(nmResDegCoh, nmInReMutInt, nmInImMutInt, rel_thresh, new_disp)
string nmResDegCoh = srwUtiGetValS("nmResDegCoh", "wDegCoh", "SrwUtiMutInt2RotDegCoh")
string nmInReMutInt = srwUtiGetValS("nmInReMutInt", "wInReMutInt", "SrwUtiMutInt2RotDegCoh")
string nmInImMutInt = srwUtiGetValS("nmInImMutInt", "wInImMutInt", "SrwUtiMutInt2RotDegCoh")
variable rel_thresh = srwUtiGetValN("rel_thresh", 1e-04, "SrwUtiMutInt2RotDegCoh")
variable new_disp = srwUtiGetValN("new_disp", 1, "SrwUtiMutInt2RotDegCoh")
prompt nmResDegCoh, "Name of Degree of Coherence wave to create"
prompt nmInReMutInt, "Name of Mutual Intensity wave, Re part", popup Wavelist("*",";","TEXT:0,DIMS:2")
prompt nmInImMutInt, "Name of Mutual Intensity wave, Im part", popup Wavelist("*",";","TEXT:0,DIMS:2")
prompt rel_thresh, "Relative Regularization Threshold"
prompt new_disp, "New Display?", popup "No;Yes"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmResDegCoh", nmResDegCoh, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValS("nmInReMutInt", nmInReMutInt, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValS("nmInImMutInt", nmInImMutInt, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValN("rel_thresh", rel_thresh, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValN("new_disp", new_disp, "SrwUtiMutInt2RotDegCoh")

string nmInMutInt = nmInReMutInt

variable xStart = dimoffset($nmInMutInt, 0)
variable xNp = dimsize($nmInMutInt, 0)
variable xStep = dimdelta($nmInMutInt, 0)
variable xEnd = xStart + (xNp - 1)*xStep
variable xc = 0.5*(xStart + xEnd)
variable xRange = xEnd - xStart

variable yStart = dimoffset($nmInMutInt, 1)
variable yNp = dimsize($nmInMutInt, 1)
variable yStep = dimdelta($nmInMutInt, 1)
variable yEnd = yStart + (yNp - 1)*yStep
variable yc = 0.5*(yStart + yEnd)
variable yRange = yEnd - yStart

variable xNpNew = 2*xNp - 1, yNpNew = 2*yNp - 1

//make/O/D/N=(xNpNew, yNpNew) wInReMutIntAux, wInImMutIntAux
//variable xHalfNp = round(xNp*0.5), yHalfNp = round(yNp*0.5)
//SetScale/P x xStart - xHalfNp*xStep, xStep, "m", wInReMutIntAux, wInImMutIntAux
//SetScale/P y yStart - yHalfNp*yStep, yStep, "m", wInReMutIntAux, wInImMutIntAux
//wInReMutIntAux = $nmInReMutInt(x)(y)*srwUtiNonZeroIntervB(x, xStart, xEnd)*srwUtiNonZeroIntervB(y, yStart, yEnd)
//wInImMutIntAux = $nmInImMutInt(x)(y)*srwUtiNonZeroIntervB(x, xStart, xEnd)*srwUtiNonZeroIntervB(y, yStart, yEnd)
//duplicate/O wInReMutIntAux wMutCohNonRot

make/O/D/N=(xNpNew, yNpNew) wMutCohNonRot
variable xHalfNp = round(xNp*0.5), yHalfNp = round(yNp*0.5)
SetScale/P x xStart - xHalfNp*xStep, xStep, "m", wMutCohNonRot
SetScale/P y yStart - yHalfNp*yStep, yStep, "m", wMutCohNonRot

variable abs_thresh = rel_thresh*abs($nmInMutInt(0)(0))
//print abs_thresh

//wMutCohNonRot = abs(wInMutCohRes(x)(y))/(sqrt(abs(wInMutCohRes(x)(x)*wInMutCohRes(y)(y))) + abs_thresh)
wMutCohNonRot = (sqrt(($nmInReMutInt(x)(y))^2 + ($nmInImMutInt(x)(y))^2)/(sqrt(abs($nmInReMutInt(x)(x)*$nmInReMutInt(y)(y))) + abs_thresh))*srwUtiNonZeroIntervB(x, xStart, xEnd)*srwUtiNonZeroIntervB(y, yStart, yEnd)

duplicate/O $nmInMutInt $nmResDegCoh

$nmResDegCoh = srwUtiInterp2DBilin(x+y, x-y, wMutCohNonRot)
//$nmResDegCoh = srwUtiNan2Num($nmResDegCoh(x)(y), 0)

string nmInMutIntCut = nmInMutInt + "x"
string nmResDegCohCut = nmResDegCoh + "cut"
string nmResDegCohCut0 = nmResDegCoh + "cut_aux"

duplicate/O $nmInMutIntCut $nmResDegCohCut
$nmResDegCohCut = $nmResDegCoh(0)(x)

duplicate/O $nmResDegCohCut $nmResDegCohCut0
$nmResDegCohCut0 = $nmResDegCoh(x)(0)

if(new_disp > 1)
	display; appendimage $nmResDegCoh
	Label bottom "(x + x*)/2"; Label left "(x - x*)/2"
	display $nmResDegCohCut0; SrwUtiGraphAddFrameAndGrid()
	Label bottom "(x + x*)/2"; Label left "Degree of Coherence"
	display $nmResDegCohCut; SrwUtiGraphAddFrameAndGrid()
	Label bottom "(x - x*)/2"; Label left "Degree of Coherence"
endif

//Determining Coherence Length
FindLevels $nmResDegCohCut, 0.5

variable nLev = dimsize(W_FindLevels, 0) 
variable i=0, curX, prevX, iFirstPos = -1
do 
	curX = W_FindLevels[i]
	
	if((i > 0) && (prevX < 0) && (curX > 0))
		iFirstPos = i
		break
	endif
	prevX = curX
	i += 1
while(i < nLev)
if(iFirstPos < 0)
	print "Coherence Length can't be determined"
else
	print "Estimated Coherence Length:", curX - prevX, "m"
endif

KillWaves/Z wMutCohNonRot
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Conversion from Mutial Intensity to (rotated) Degree of Coherence function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwUtiMutInt2RotDegCoh(nmResDegCoh, nmInMutInt, rel_thresh, new_disp)
string nmResDegCoh = srwUtiGetValS("nmResDegCoh", "wDegCoh", "SrwUtiMutInt2RotDegCoh")
string nmInMutInt = srwUtiGetValS("nmInMutInt", "wInMutInt", "SrwUtiMutInt2RotDegCoh")
variable rel_thresh = srwUtiGetValN("rel_thresh", 1e-04, "SrwUtiMutInt2RotDegCoh")
variable new_disp = srwUtiGetValN("new_disp", 1, "SrwUtiMutInt2RotDegCoh")
prompt nmResDegCoh, "Name of Degree of Coherence wave to create"
prompt nmInMutInt, "Name of Mutual Intensity wave", popup Wavelist("*",";","TEXT:0,DIMS:2")
prompt rel_thresh, "Relative Regularization Threshold"
prompt new_disp, "New Display?", popup "No;Yes"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("nmResDegCoh", nmResDegCoh, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValS("nmInMutInt", nmInMutInt, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValN("rel_thresh", rel_thresh, "SrwUtiMutInt2RotDegCoh")
srwUtiSetValN("new_disp", new_disp, "SrwUtiMutInt2RotDegCoh")

variable xStart = dimoffset($nmInMutInt, 0)
variable xNp = dimsize($nmInMutInt, 0)
variable xStep = dimdelta($nmInMutInt, 0)
variable xEnd = xStart + (xNp - 1)*xStep
variable xc = 0.5*(xStart + xEnd)
variable xRange = xEnd - xStart
//variable xHalfRange = 0.5*xRange

variable yStart = dimoffset($nmInMutInt, 1)
variable yNp = dimsize($nmInMutInt, 1)
variable yStep = dimdelta($nmInMutInt, 1)
variable yEnd = yStart + (yNp - 1)*yStep
variable yc = 0.5*(yStart + yEnd)
variable yRange = yEnd - yStart
//variable yHalfRange = 0.5*yRange

variable xNpNew = 2*xNp - 1, yNpNew = 2*yNp - 1
//print xNpNew, yNpNew

make/O/D/N=(xNpNew, yNpNew) wInMutCohRes
//SetScale/I x xc - xRange, xc + xRange, "m", wInMutCohRes
//SetScale/I y yc - yRange, yc + yRange, "m", wInMutCohRes

variable xHalfNp = round(xNp*0.5), yHalfNp = round(yNp*0.5)
SetScale/P x xStart - xHalfNp*xStep, xStep, "m", wInMutCohRes
SetScale/P y yStart - yHalfNp*yStep, yStep, "m", wInMutCohRes
wInMutCohRes = $nmInMutInt(x)(y)*srwUtiNonZeroIntervB(x, xStart, xEnd)*srwUtiNonZeroIntervB(y, yStart, yEnd)

duplicate/O $nmInMutInt wMutCohNonRot
duplicate/O $nmInMutInt $nmResDegCoh

variable abs_thresh = rel_thresh*abs($nmInMutInt(0)(0))
//print abs_thresh

//wMutCohNonRot = abs($nmInMutInt(x)(y))/(sqrt(abs($nmInMutInt(x)(x)*$nmInMutInt(y)(y))) + abs_thresh)
wMutCohNonRot = abs(wInMutCohRes(x)(y))/(sqrt(abs(wInMutCohRes(x)(x)*wInMutCohRes(y)(y))) + abs_thresh)

//TEST
//wMutCohNonRot = abs(wInMutCohRes(x)(y))

//$nmResDegCoh = wMutCohNonRot(x + y)(x - y)
//$nmResDegCoh = Interp2D(wMutCohNonRot, x+y, x-y)
$nmResDegCoh = srwUtiInterp2DBilin(x+y, x-y, wMutCohNonRot)

//$nmResDegCoh = wMutCohNonRot(x)(y)

$nmResDegCoh = srwUtiNan2Num($nmResDegCoh(x)(y), 0)

string nmInMutIntCut = nmInMutInt + "x"
string nmResDegCohCut = nmResDegCoh + "cut"
string nmResDegCohCut0 = nmResDegCoh + "cut_aux"

duplicate/O $nmInMutIntCut $nmResDegCohCut
$nmResDegCohCut = $nmResDegCoh(0)(x)

duplicate/O $nmResDegCohCut $nmResDegCohCut0
$nmResDegCohCut0 = $nmResDegCoh(x)(0)

if(new_disp > 1)
	display; appendimage $nmResDegCoh
	Label bottom "(x + x*)/2"; Label left "(x - x*)/2"
	display $nmResDegCohCut0; SrwUtiGraphAddFrameAndGrid()
	Label bottom "(x + x*)/2"; Label left "Degree of Coherence"
	display $nmResDegCohCut; SrwUtiGraphAddFrameAndGrid()
	Label bottom "(x - x*)/2"; Label left "Degree of Coherence"
endif

//Determining Coherence Length
FindLevels $nmResDegCohCut, 0.5

variable nLev = dimsize(W_FindLevels, 0) 
variable i=0, curX, prevX, iFirstPos = -1
do 
	curX = W_FindLevels[i]
	
	if((i > 0) && (prevX < 0) && (curX > 0))
		iFirstPos = i
		break
	endif
	prevX = curX
	i += 1
while(i < nLev)
if(iFirstPos < 0)
	print "Coherence Length can't be determined"
else
	print "Estimated Coherence Length:", curX - prevX, "m"
endif
end