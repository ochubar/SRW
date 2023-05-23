
//+++++++++++++++++++++++++++++++++++++++
//
//Compute Power Density (General)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowDenCreate(RadName, ElecName, MagName, ObsName, PrecPar, Meth, UseSpecLim, sStart, sFin, Disp)
String RadName=srwUtiTruncString(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType
String MagName=SrwMagGenTotName
String ObsName=SrwSmpGenTotName
Variable PrecPar=SrwPowCompPrec
Variable Meth=SrwPowCompMeth
variable UseSpecLim=srwUtiGetValN("UseSpecLim", 1, "SrwPowDenCreate")
variable sStart=srwUtiGetValN("sStart", -10, "SrwPowDenCreate")
variable sFin=srwUtiGetValN("sFin", 10, "SrwPowDenCreate")
Variable Disp=SrwPowDispImmed
prompt RadName,SrwPPowName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwMagConstType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt PrecPar,SrwPPowCompPrec
prompt Meth,SrwPPowCompMeth,popup "Near Field;Far Field"
prompt UseSpecLim,"Use Spec. Integration Limits",popup "No;Yes"
prompt sStart,"Initial Longit. Position [m]"
prompt sFin,"Final Longit. Position [m]"
prompt Disp,SrwPPowDispImmed,popup "No;Yes"
Silent 1						|	Computing the Power Density  ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0;
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0;
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0;
endif
if(ElecWavePresent==0)
	SrwElecFilament();
	SrwElecThick();
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwPowCreate(); 
		Return;
	endif
endif
if(MagWavePresent==0)
	DoAlert 0, SrwPAlertMagFieldNeeded;
	Abort;
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 4;
	SrwStartMacrosAfterRadSmp2 *= -1;
	SrwRadSamplDialog(SrwSmpPowType);
	Return;
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1];
SrwSmpGenTotName=ObsName;

SrwPowCompPrec=PrecPar;
SrwPowCompMeth=Meth;
SrwPowDispImmed=Disp;

srwUtiSetValN("UseSpecLim", UseSpecLim, "SrwPowDenCreate")
srwUtiSetValN("sStart", sStart, "SrwPowDenCreate")
srwUtiSetValN("sFin", sFin, "SrwPowDenCreate")

string BufMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1];
string MagType = MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1];
if(cmpstr(MagType,SrwFieldType)==0)
	SrwMagName = BufMagName;
endif
if(cmpstr(MagType,SrwUndType)==0)
	SrwUndName = BufMagName;
endif
if(cmpstr(MagType,SrwMagConstType)==0)
	SrwMagConstName = BufMagName;
endif
SrwMagGenTotName = MagName;

if(strlen(RadName)==0)
	RadName=SrwElecName+BufMagName+SrwSmpName;
endif
SrwPowName=RadName;

// Preparing data for C function
Make/D/O/N=6 waveprec;
waveprec[0]=PrecPar;
waveprec[1]=Meth;
waveprec[2]=UseSpecLim
waveprec[3]=sStart
waveprec[4]=sFin

SrwPowPrep(ElecName,MagName,Obsname,RadName);
String TotRadName=RadName+SrwPowType;
SrwRadGenTotName=TotRadName

srPowDens($ElecName, $MagName, $ObsName, waveprec, $TotRadName);
KillWaves/Z  waveprec;

String IntSuffix="I";
if(Disp==2)
	SrwPow2Int(TotRadName,IntSuffix,4,0,0,2);
else
	if(($Obsname[10]>1) %& ($Obsname[13]>1)) // if 2D
		String IntWaveNameXZ=RadName+IntSuffix+SrwRadXZType;
		if(exists(IntWaveNameXZ)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
	if($Obsname[13]<=1) // if hor. section
		String IntWaveNameX = RadName+IntSuffix+SrwSeparator+SrwRadXType;
		if(exists(IntWaveNameX)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
	if($Obsname[10]<=1) // if vert. section
		String IntWaveNameZ = RadName+IntSuffix+SrwSeparator+SrwRadZType;
		if(exists(IntWaveNameZ)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Power Density (General) Old Version
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowCreate(RadName, ElecName, MagName, ObsName, PrecPar, Meth, Disp)
String RadName=srwUtiTruncString(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName, 25)
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwMagGenTotName;
String ObsName=SrwSmpGenTotName;
Variable PrecPar=SrwPowCompPrec;
Variable Meth=SrwPowCompMeth;
Variable Disp=SrwPowDispImmed;
prompt RadName,SrwPPowName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","");
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwMagConstType,";","");
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","");
prompt PrecPar,SrwPPowCompPrec;
prompt Meth,SrwPPowCompMeth,popup "Near Field;Far Field";
prompt Disp,SrwPPowDispImmed,popup "No;Yes";

Silent 1						|	Computing the Power Density  ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0;
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0;
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0;
endif
if(ElecWavePresent==0)
	SrwElecFilament();
	SrwElecThick();
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwPowCreate(); 
		Return;
	endif
endif
if(MagWavePresent==0)
	DoAlert 0, SrwPAlertMagFieldNeeded;
	Abort;
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 4;
	SrwStartMacrosAfterRadSmp2 *= -1;
	SrwRadSamplDialog(SrwSmpPowType);
	Return;
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1];
SrwSmpGenTotName=ObsName;

SrwPowCompPrec=PrecPar;
SrwPowCompMeth=Meth;
SrwPowDispImmed=Disp;

string BufMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1];
string MagType = MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1];
if(cmpstr(MagType,SrwFieldType)==0)
	SrwMagName = BufMagName;
endif
if(cmpstr(MagType,SrwUndType)==0)
	SrwUndName = BufMagName;
endif
if(cmpstr(MagType,SrwMagConstType)==0)
	SrwMagConstName = BufMagName;
endif
SrwMagGenTotName = MagName;

variable StrLenRadName = strlen(RadName);
if(StrLenRadName > 25)
	abort "The name is too long."
endif
if(StrLenRadName==0)
	RadName=SrwElecName+BufMagName+SrwSmpName;
	RadName=srwUtiTruncString(RadName, 25);
endif

SrwPowName=RadName;

// Preparing data for C function
Make/D/O/N=6 waveprec;
waveprec[0]=PrecPar;
waveprec[1]=Meth;

SrwPowPrep(ElecName,MagName,Obsname,RadName);
String TotRadName=RadName+SrwPowType;
SrwRadGenTotName=TotRadName

srPowDens($ElecName, $MagName, $ObsName, waveprec, $TotRadName);
KillWaves/Z  waveprec;

String IntSuffix="I";
if(Disp==2)
	SrwPow2Int(TotRadName,IntSuffix,4,0,0,2);
else
	if(($Obsname[10]>1) %& ($Obsname[13]>1)) // if 2D
		String IntWaveNameXZ=RadName+IntSuffix+SrwRadXZType;
		if(exists(IntWaveNameXZ)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
	if($Obsname[13]<=1) // if hor. section
		String IntWaveNameX = RadName+IntSuffix+SrwSeparator+SrwRadXType;
		if(exists(IntWaveNameX)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
	if($Obsname[10]<=1) // if vert. section
		String IntWaveNameZ = RadName+IntSuffix+SrwSeparator+SrwRadZType;
		if(exists(IntWaveNameZ)==1)
			SrwPow2Int(TotRadName,IntSuffix,4,0,0,1);
		endif
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Prepare Power Density structure (numerical wave)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowPrep(ElecName,MagName,ObsName,PowName)
String ElecName,MagName,ObsName,PowName;

PowName += SrwPowType;

variable xSt=($ObsName[8]), xFi=($ObsName[9]);
if(($ObsName[10]==1)%&(xSt==xFi)) 
	xFi = xSt*1.0001; // Can be anything
endIf
variable zSt=($ObsName[11]), zFi=($ObsName[12]);
if(($ObsName[13]==1)%&(zSt==zFi)) 
	zFi = zSt*1.0001; // Can be anything
endIf

variable nx = $ObsName[10], nz = $ObsName[13];
if(nx < 1)
	nx = 1;
endif
if(nz < 1)
	nz = 1;
endif

Make/O/N=(nx, nz) $PowName;
SetScale/I x xSt, xFi, "m", $PowName;
SetScale/I y zSt, zFi, "m", $PowName;
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract  results
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPow2Int(PowName, SuffixExtract, PlotType, xVal, zVal, dis)
String PowName=SrwPowName+SrwPowType
String SuffixExtract=SrwSuffixExtract
Variable PlotType=SrwViewPlotTypeXZ
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt PowName,SrwPPowName, popup Wavelist("*"+SrwPowType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt PlotType, SrwPViewPlotType, popup SrwPOPUPViewPlotTypeXZ
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,"New Display",popup "No;Yes"
Silent 1						|	  ...
PauseUpdate

SrwPowName=PowName[0,strlen(PowName)-strlen(SrwPowType)-1]
SrwSuffixExtract=SuffixExtract
SrwViewPlotTypeXZ=PlotType
SrwViewX=xVal
SrwViewZ=zVal

if(cmpstr(PowName,"_none_")==0)
	SrwStartMacrosAfterRadSmp2 = 4 // To proceed default chain through RadSampling panel
	SrwPowCreate()
	if(SrwStartMacrosAfterRadSmp2 > 0)
		SrwPow2Int()
		SrwStartMacrosAfterRadSmp2 = 0
	endif
	Return
endif

xVal *= 0.001
zVal *= 0.001

// Treating "Auto"
String st=PowName
if(PlotType==4)
	PlotType=srwPow2IntTreatAuto($st)
endif

String ViewPowName=SrwPowName+SuffixExtract

if(PlotType==1)
	ViewPowName+=SrwSeparator+SrwRadXType
	Make/O/N=(DimSize($st, 0)) $ViewPowName
	SetScale/P x DimOffset($st, 0), DimDelta($st, 0), WaveUnits($st, 0), $ViewPowName
	//SetScale d, 0, 0, SrwPUnitPowDen, $ViewPowName
	$ViewPowName=$st(x)(zVal)
	if(dis==2)
		if(DimSize($st, 0)>1)
			Display;Append $ViewPowName
			Label bottom SrwPLabelHorPos
			Label left SrwPUnitPowDen
		else
			print $ViewPowName[0], SrwPUnitPowDen
		endif
	endif
endif
if(PlotType==2)
	ViewPowName+=SrwSeparator+SrwRadZType
	Make/O/N=(DimSize($st, 1)) $ViewPowName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewPowName
	//SetScale d, 0, 0, SrwPUnitPowDen, $ViewPowName
	$ViewPowName=$st(xVal)(x)
	if(dis==2)
		if(DimSize($st, 1)>1)
			Display;Append $ViewPowName
			Label bottom SrwPLabelVerPos
			Label left SrwPUnitPowDen
		else
			print $ViewPowName[0], SrwPUnitPowDen
		endif
	endif
endif
if(PlotType==3)
	ViewPowName+=SrwSeparator+SrwRadXType+SrwRadZType
	Make/O/N=((DimSize($st, 0)), (DimSize($st, 1))) $ViewPowName
	SetScale/P x DimOffset($st, 0), DimDelta($st, 0), WaveUnits($st, 0), $ViewPowName
	SetScale/P y DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewPowName
	$ViewPowName=$st[p][q]
	if(dis==2)
		if((DimSize($st, 0)>1) %| (DimSize($st, 1)>1))
			Display;AppendImage $ViewPowName
			SrwImageFormat(ViewPowName)
			Label bottom SrwPLabelHorPos
			Label left SrwPLabelVerPos
		else
			print $ViewPowName[0][0], SrwPUnitPowDen
		endif
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
function srwPow2IntTreatAuto(pow)
wave pow
if(DimSize(pow, 0)==1)
	return 2
else
	if(DimSize(pow, 1)==1)
		return 1
	else
		return 3
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Duplicate Power Density structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowDupl(PowName,DuplPowName)
string PowName=SrwPowName
string DuplPowName=SrwPowName+"d"
prompt PowName,SrwPPowName1, popup Wavelist("*"+SrwPowType,";","")
prompt DuplPowName,SrwPPowNameDpl
Silent 1						|	  ...
PauseUpdate

if(cmpstr(PowName,"_none_")==0)
	DoAlert 0, SrwPAlertNoCompResultsFound
	Return
endif

SrwPowName=DuplPowName

String NewPowName=SrwPowName+SrwPowType
Duplicate/O $PowName $NewPowName
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculate Undulator Power and Power Density vs Fund. Photon Energy
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowUndVsFundPhotEn(nmPowCore, nmElec, nmMag, nmObs, type, Kmin, numK, precPar, disp)
string nmPowCore=srwUtiTruncString(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName, 25)
string nmElec=SrwElecName+SrwElecType
string nmMag=SrwMagGenTotName
string nmObs=SrwSmpGenTotName
variable type=srwUtiGetValN("type", 1, "SrwPowUndVsFundPhotEn")
variable Kmin=srwUtiGetValN("Kmin", 0.3, "SrwPowUndVsFundPhotEn")
variable numK=srwUtiGetValN("numK", 100, "SrwPowUndVsFundPhotEn")
variable precPar=SrwPowCompPrec
variable meth=SrwPowCompMeth
variable disp=SrwPowDispImmed
prompt nmPowCore,"Name Core for Power data"
prompt nmElec,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt nmMag,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","")//+Wavelist("*"+SrwUndType,";","")
prompt nmObs,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt type,"Type of Power Density to produce",popup "Per Unit Surface;Per Unit Solid Angle"
prompt Kmin,"Minimal Deflecting Parameter"
prompt numK,"Number of Deflecting Param. values"
prompt precPar,SrwPPowCompPrec
prompt disp,SrwPPowDispImmed,popup "No;Yes"
Silent 1						|	Computing Power Density  ...
//PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(nmElec,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(nmMag,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(nmObs,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwPowUndVsFundPhotEn() 
		Return
	endif
endif
if(MagWavePresent==0)
	DoAlert 0, SrwPAlertMagFieldNeeded
	Abort
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 4
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpPowType)
	Return
endif

if(Kmin <= 0)
	Abort "Minimal Deflecting Parameter shound be positive"
endif
if(numK <= 1)
	Abort "Number of Deflecting Parameter values should be greater than 1"
endif

//Extra parameters (to edit)
variable intMeth = 2 //2- undulator, 3- wiggler
variable intPrec = 0.006
variable relMinPhEn = 0.6
variable relMaxPhEn = 1.4
//End of extra parameters

SrwElecName=nmElec[0,strlen(nmElec)-strlen(SrwElecType)-1]
SrwSmpName=nmObs[0,strlen(nmObs)-strlen(SrwSmpPowType)-1]
SrwSmpGenTotName=nmObs
srwUtiSetValN("type", type, "SrwPowUndVsFundPhotEn")
srwUtiSetValN("Kmin", Kmin, "SrwPowUndVsFundPhotEn")
srwUtiSetValN("numK", numK, "SrwPowUndVsFundPhotEn")
SrwPowCompPrec=precPar
SrwPowCompMeth=meth

string BufMagName = nmMag[0,strlen(nmMag)-strlen(SrwFieldType)-1]
string MagType = nmMag[strlen(nmMag)-strlen(SrwFieldType),strlen(nmMag)-1]

if(cmpstr(MagType,SrwFieldType)==0)
	SrwMagName = BufMagName
endif
//if(cmpstr(MagType,SrwUndType)==0)
//	SrwUndName = BufMagName
//endif
SrwMagGenTotName = nmMag

variable StrLenRadName = strlen(nmPowCore)
if(StrLenRadName > 25)
	abort "The name is too long."
endif
if(StrLenRadName==0)
	nmPowCore=SrwElecName+BufMagName+SrwSmpName
	nmPowCore=srwUtiTruncString(nmPowCore, 25)
endif

SrwPowName=nmPowCore

string nmResK = nmPowCore + "_K"
string nmResPhEn = nmPowCore + "_E"
string nmResPowTot = nmPowCore + "_P"
string nmResPowDensOnAx = nmPowCore + "_PA"
string nmAuxPowDenSuf = ""
string nmResPowDens = nmPowCore + nmAuxPowDenSuf + "_xz"

make/O/N=(numK) $nmResK, $nmResPhEn, $nmResPowTot, $nmResPowDensOnAx
SetScale d 0,0,"eV", $nmResPhEn
SetScale d 0,0,"W", $nmResPowTot

string nmAuxMagCore = "AuxMagPowUndVsFundPhEn"
string nmAuxMag = nmAuxMagCore + "_mag"
string nmAuxMagBX = nmAuxMagCore + "BX_fld"
string nmAuxMagBZ = nmAuxMagCore + "BZ_fld"
string nmAuxMagPerCore = "AuxMagPerPowUndVsFundPhEn"
string nmAuxMagPer = nmAuxMagPerCore + "_map"

string nmAuxSpecCore = "AuxWfrSpecPowUndVsPhEn"
string nmAuxSpecSingleE = nmAuxSpecCore + "_rad"
string nmAuxSpecIntSuf = "Is"
string nmAuxSpecInt = nmAuxSpecCore + nmAuxSpecIntSuf + "_e"
string nmAuxElecDriftCore = "AuxDriftElecPowUndVsPhEn"
string nmAuxElecDrift = nmAuxElecDriftCore + "_bli"
string nmAuxElec = "AuxElecPowUndVsFundPhEn_ebm"

string nmAuxMagPerHarm

string nmObsAux = "AuxObsPowUndVsFundPhotEn_obs"
duplicate/O $nmObs $nmObsAux
srwSetSmpHorPosNp(nmObsAux, 1)
srwSetSmpVertPosNp(nmObsAux, 1)

variable xStart = srwGetSmpHorPosStart(nmObs), xEnd = srwGetSmpHorPosEnd(nmObs)
variable xcPowDens = 0.5*(xStart + xEnd)
variable yStart = srwGetSmpVertPosStart(nmObs), yEnd = srwGetSmpVertPosEnd(nmObs)
variable ycPowDens = 0.5*(yStart + yEnd)
variable Robs = srwGetSmpLongPos(nmObs)

variable xAngStart = xStart/Robs, xAngEnd = xEnd/Robs
variable yAngStart = yStart/Robs, yAngEnd = yEnd/Robs

SrwMagDupl(nmMag, nmAuxMagCore)
variable sStart = dimoffset($nmAuxMagBZ, 0)
variable sStep = dimdelta($nmAuxMagBZ, 0)
variable ns = dimsize($nmAuxMagBZ, 0)
variable sRange = sStep*(ns - 1)
variable maxUndPer = 0.25*sRange*1000

variable elecEn = srwGetElecBeamEnergy(nmElec)
variable elecS0 = srwGetElecBeamLongPos(nmElec)
variable elecDriftLen = (sStart + sStep) - elecS0

if(abs(elecDriftLen) > 10*sStep)
	duplicate/O $nmElec $nmAuxElec
	SrwOptDrift(nmAuxElecDriftCore, elecDriftLen)
	srElecBeamPropag($nmAuxElec, $nmAuxElecDrift)
	nmElec = nmAuxElec
endif

SrwMagArb2Per(nmAuxMagPerCore, nmMag, 0.03, 5, maxUndPer)
variable effUndPer = str2num($nmAuxMagPer[0]) //undulator period in [m]
variable nMagPerHarm = str2num($nmAuxMagPer[5]) 
variable KxE2 = 0, KyE2 = 0
variable iMagPerHarm = 0
do
	nmAuxMagPerHarm = $nmAuxMagPer[6 + iMagPerHarm]
	if($nmAuxMagPerHarm[1] == 1)
		KyE2 += ($nmAuxMagPerHarm[2]/$nmAuxMagPerHarm[0])^2
	else
		KxE2 += ($nmAuxMagPerHarm[2]/$nmAuxMagPerHarm[0])^2	
	endif
	iMagPerHarm += 1
while(iMagPerHarm < nMagPerHarm)
variable KeffMax = sqrt(KxE2 + KyE2)
variable factB = (Kmin/KeffMax)^(1/(numK - 1))
variable Beff = KeffMax/(0.0933729e+03*effUndPer)
variable approxResonPhotEn = srUtiUndFundPhotEn(Beff, effUndPer, elecEn, 2) // [eV]
variable infFundPhotEn, supFundPhotEn

variable iK = 0, Kprev, Kcur = KeffMax, dispPowDens = 1
if(disp == 2)
	dispPowDens = 2
endif
do
	if(iK > 0)
		$nmAuxMagBX *= factB
		$nmAuxMagBZ *= factB
		
		Kcur = factB*Kprev
		approxResonPhotEn *= (2 + Kprev*Kprev)/(2 + Kcur*Kcur)
		
		dispPowDens = 1
		nmAuxPowDenSuf = "aux"
		nmResPowDens = nmPowCore + nmAuxPowDenSuf + "_xz"
	endif
	
	infFundPhotEn = relMinPhEn*approxResonPhotEn
	supFundPhotEn = relMaxPhEn*approxResonPhotEn
	
	srwSetSmpPhotEnStart(nmObsAux, infFundPhotEn)
	srwSetSmpPhotEnEnd(nmObsAux, supFundPhotEn)

	SrwMagPrec(nmAuxMag, intMeth, intPrec, intPrec, 10000, 1, 0, 0)
	SrwWfrCreate(nmAuxSpecCore, nmElec, nmAuxMag, nmObsAux, 1, 1)
	SrwWfr2Int(nmAuxSpecSingleE, nmAuxSpecIntSuf, 7, 1, 1, 1, 1, 0, 0, 1) //on-axis single-e intensity vs photon energy

	WaveStats/Q/R=(infFundPhotEn, supFundPhotEn) $nmAuxSpecInt
	if(V_maxloc <= infFundPhotEn)
		print "WARNING: Lower limit reached when searching harmonic No:", iHarm, "(finite-emittance e-beam)"
	endif
	if(V_maxloc >= supFundPhotEn)
		print "WARNING: Upper limit reached when searching harmonic No:", iHarm, "(finite-emittance e-beam)"
	endif
	approxResonPhotEn = V_maxloc
	$nmResPhEn[iK] = V_maxloc
	$nmResK[iK] = Kcur
	
	SrwPowDenCreate(nmPowCore, nmElec, nmAuxMag, nmObs, precPar, 1,1,-10,10, 1)
	SrwPow2Int(nmPowCore+"_pow", nmAuxPowDenSuf, 3,0,0, dispPowDens)
	
	$nmResPowTot[iK] = srwUtiIntTotWave2DT($nmResPowDens)*1e+06 //total power in [W]
	$nmResPowDensOnAx[iK] = $nmResPowDens(xcPowDens)(ycPowDens)
	
	if(type == 2) //convert to power per unit angle
		$nmResPowDensOnAx[iK] *= Robs*Robs
		
		if(iK == 0)
			$nmResPowDens *= Robs*Robs //to have [W/mrad^2]
			SetScale/I x xAngStart,xAngEnd,"rad", $nmResPowDens
			SetScale/I y yAngStart,yAngEnd,"rad", $nmResPowDens
		
			Label left "Vertical Angle"; 	Label bottom "Horizontal Angle"
		endif
	endif

	if(iK == 0)
		if(disp == 2)
			display $nmResPowTot vs $nmResK; Label bottom "Deflecting Parameter"; Label left "Power [W]"; SrwUtiGraphAddFrameAndGrid()
			display $nmResPowTot vs $nmResPhEn; Label bottom "Fundamental Photon Energy"; Label left "Power [W]"; SrwUtiGraphAddFrameAndGrid()
	
			display $nmResPowDensOnAx vs $nmResK; Label bottom "Deflecting Parameter"; Label left "Power per Unit Surface [W/mm\\S2\\M]"; SrwUtiGraphAddFrameAndGrid()
			if(type == 2)
				Label left "Power per Unit Solid Angle [W/mrad\\S2\\M]"
			endif
			display $nmResPowDensOnAx vs $nmResPhEn; Label bottom "Fundamental Photon Energy"; Label left "Power per Unit Surface [W/mm\\S2\\M]"; SrwUtiGraphAddFrameAndGrid()
			if(type == 2)
				Label left "Power per Unit Solid Angle [W/mrad\\S2\\M]"
			endif
		endif
	endif
	
	Kprev = Kcur
	iK += 1
while(iK < numK)

SrwPowDispImmed=disp
killwaves/Z $nmObsAux
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculate Power Density distributions for a set of Electron Beam Initial Conditions
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPowDenVarElecInitCond(nmPowCore, nmElecBase, nmElecInitCond, nmMag, nmObs, precPar, meth, useSpecLim, sStart, sFin) //, disp)
string nmPowCore = srwUtiGetValS("nmPowCore", "PowDens", "SrwPowDenVarInitCond")
string nmElecBase = srwUtiGetValS("nmElecBase", "Elec", "SrwPowDenVarInitCond")
string nmElecInitCond = srwUtiGetValS("nmElecInitCond", "wElecInitCond", "SrwPowDenVarInitCond")
string nmMag = srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string nmObs = srwUtiGetValS("SrwSmpGenTotName", "Obs", "")
variable precPar = srwUtiGetValN("SrwPowCompPrec", 1, "")
variable meth = srwUtiGetValN("SrwPowCompMeth", 1, "")
variable useSpecLim = srwUtiGetValN("UseSpecLim", 1, "SrwPowDenCreate")
variable sStart = srwUtiGetValN("sStart", -10, "SrwPowDenCreate")
variable sFin = srwUtiGetValN("sFin", 10, "SrwPowDenCreate")
//variable disp = srwUtiGetValN("SrwPowDispImmed", 2, "")
prompt nmPowCore, "Core name for Power Density structures"
prompt nmElecBase, "Base Electron Beam structure", popup Wavelist("*"+SrwElecType,";","")
prompt nmElecInitCond, "2D Num. Wave with Elec. Initial Cond.", popup Wavelist("*",";","TEXT:0,DIMS:2,MINCOLS:4,MAXCOLS:4")
prompt nmMag, SrwPMagName2, popup Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwMagConstType,";","")
prompt nmObs, SrwPSmpName2, popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt precPar,SrwPPowCompPrec
prompt meth,SrwPPowCompMeth,popup "Near Field;Far Field"
prompt useSpecLim,"Use Spec. Integration Limits",popup "No;Yes"
prompt sStart,"Initial Longit. Position [m]"
prompt sFin,"Final Longit. Position [m]"
//prompt Disp,SrwPPowDispImmed,popup "No;Yes"
Silent 1						|	Computing Power Density  ...
PauseUpdate

if(strlen(nmPowCore) > 25)
	abort "Power density core name is too long"
endif
if((strlen(nmElecInitCond) <= 0) %| (exists(nmElecInitCond) != 1))
	abort "2D numerical wave with electron initial conditions is not defined"
endif

srwUtiSetValS("nmPowCore", nmPowCore, "SrwPowDenVarInitCond")
srwUtiSetValS("nmElecBase", nmElecBase, "SrwPowDenVarInitCond")
srwUtiSetValS("nmElecInitCond", nmElecInitCond, "SrwPowDenVarInitCond")
srwUtiSetValS("SrwMagGenTotName", nmMag, "")
srwUtiSetValS("SrwSmpGenTotName", nmObs, "")
srwUtiSetValN("SrwPowCompPrec", precPar, "")
srwUtiSetValN("SrwPowCompMeth", meth, "")
srwUtiSetValN("UseSpecLim", useSpecLim, "SrwPowDenCreate")
srwUtiSetValN("sStart", sStart, "SrwPowDenCreate")
srwUtiSetValN("sFin", sFin, "SrwPowDenCreate")

//indexes of columns in 2D initial conditions wave 
variable indElecPosX = 0 //horizontal position [m]
variable indElecAngX = 1 //horizontal angle [rad]
variable indElecPosY = 2 //vertical position [m]
variable indElecAngY = 3 //vertical angle [rad]

variable disp = 2 //do display all calculations

string nmElecAux = "AuxElecVarInitCond" + SrwElecType
duplicate/O $nmElecBase $nmElecAux

string nmPowCoreAux
variable xe, xpe, ye, ype

variable numElec = dimsize($nmElecInitCond, 0)
make/O/T/N=(numElec, 5) $nmPowCore

variable iElec = 0
do
	xe = $nmElecInitCond[iElec][indElecPosX]
	xpe = $nmElecInitCond[iElec][indElecAngX]
	ye = $nmElecInitCond[iElec][indElecPosY]
	ype = $nmElecInitCond[iElec][indElecAngY]
	
	srwSetElecBeamHorPos(nmElecAux, xe*0.001)
	srwSetElecBeamHorAng(nmElecAux, xpe*0.001)
	srwSetElecBeamVertPos(nmElecAux, ye*0.001)
	srwSetElecBeamVertAng(nmElecAux, ype*0.001)
	
	nmPowCoreAux = nmPowCore + "_" + num2str(iElec + 1)
	SrwPowDenCreate(nmPowCoreAux, nmElecAux, nmMag, nmObs, precPar, meth, useSpecLim, sStart, sFin, disp)
	
	$nmPowCore[iElec][0] = nmPowCoreAux + "_pow"
	$nmPowCore[iElec][indElecPosX + 1] = num2str(xe)
	$nmPowCore[iElec][indElecAngX + 1] = num2str(xpe)
	$nmPowCore[iElec][indElecPosY + 1] = num2str(ye)
	$nmPowCore[iElec][indElecAngY + 1] = num2str(ype)

	iElec += 1
while(iElec < numElec)

edit $nmPowCore
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculates integrated Power Density in [W/mm^2] deposited at transverse position (xx, yy) 
//in material thickness dz, from zz - dz/2 to zz + dz/2 (if z > dz/2), or from 0 to dz (if zz < dz/2),
//using Spectral Stokes parameters calculated on 3D mesh (e, x, y), 
//taking into account Spectral Attenuation Length of the material.
//Assumes Input Stokes Parameters in [Ph/s/0.1%bw/mm^2], Spectral Attenuation Length in [m], positions in [m].
//To obtain volume power density in [W/mm^3], this output must be multiplied by 0.001/dz.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwPowDensFromSpecStokes(xx, yy, zz, dz, nmStokes, nmSpecAttenLen)
variable xx, yy, zz, dz //z is longitudinal coord. (assuming z=0 where absorbtion starts)
string nmStokes, nmSpecAttenLen
PauseUpdate

if((!(xx <= 0)) %& (!(xx >= 0)))
	return NaN
endif
if((!(yy <= 0)) %& (!(yy >= 0)))
	return NaN
endif

if((zz < 0) %| (dz <= 0))
	return 0
endif

variable halfDz = 0.5*dz
variable zSt = zz - halfDz, zFi = zz + halfDz
if(zSt < 0)
	zSt = 0
	zFi = zSt + dz
endif

string nmSpecInt = nmStokes[0,strlen(nmStokes)-strlen("ras")-1] + "e"
variable eN = dimsize($nmStokes, 1)
variable eStart = dimoffset($nmStokes, 1)
variable eStep = dimdelta($nmStokes, 1)
make/O/N=(eN) $nmSpecInt
SetScale/P x eStart, eStep, "eV", $nmSpecInt
wave wSpecInt = $nmSpecInt
wave wStokes = $nmStokes
wave wSpecAttenLen = $nmSpecAttenLen

wSpecInt = (wStokes[0](x)(xx)(yy)*(exp(-zSt/wSpecAttenLen(x)) - exp(-zFi/wSpecAttenLen(x)))*(1.60218e-16))
integrate/T wSpecInt //[W/mm^2]
variable resPowDen = wSpecInt[dimsize(wSpecInt, 0) - 1] //[W/mm^2]
killwaves/Z wSpecInt
return resPowDen
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculates Volume Power Density in [W/mm^3] deposited at position (xx, yy, zz) 
//using Spectral Stokes parameters calculated on 3D mesh (e, x, y), 
//taking into account Spectral Attenuation Length of the material.
//Assumes Input Stokes Parameters in [Ph/s/0.1%bw/mm^2], Spectral Attenuation Length in [m], positions in [m].
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwPowDensVolFromSpecStokes(xx, yy, zz, nmStokes, nmSpecAttenLen)
variable xx, yy, zz //zz is longitudinal coord. (assuming zz=0 where absorbtion starts)
string nmStokes, nmSpecAttenLen
PauseUpdate

if((!(xx <= 0)) %& (!(xx >= 0)))
	return NaN
endif
if((!(yy <= 0)) %& (!(yy >= 0)))
	return NaN
endif

if(zz < 0)
	return 0
endif

string nmSpecInt = nmStokes[0,strlen(nmStokes)-strlen("ras")-1] + "e"
variable eN = dimsize($nmStokes, 1)
variable eStart = dimoffset($nmStokes, 1)
variable eStep = dimdelta($nmStokes, 1)
make/O/N=(eN) $nmSpecInt
SetScale/P x eStart, eStep, "eV", $nmSpecInt
wave wSpecInt = $nmSpecInt
wave wStokes = $nmStokes
wave wSpecAttenLen = $nmSpecAttenLen

//wSpecInt = (wStokes[0](x)(xx)(yy)*(exp(-zSt/wSpecAttenLen(x)) - exp(-zFi/wSpecAttenLen(x)))*(1.60218e-16))
wSpecInt = wStokes[0](x)(xx)(yy)*(1.60218e-19)*exp(-zz/wSpecAttenLen(x))/wSpecAttenLen(x)

integrate/T wSpecInt //[W/mm^3]
variable resPowDen = wSpecInt[dimsize(wSpecInt, 0) - 1] //[W/mm^3]
killwaves/Z wSpecInt
return resPowDen
end
