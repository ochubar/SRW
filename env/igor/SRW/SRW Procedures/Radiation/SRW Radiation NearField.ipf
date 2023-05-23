
//+++++++++++++++++++++++++++++++++++++++
//
//Create Radiation Field Structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreate_(RadName, ElecName, MagName, ObsName)
//String RadName=SrwRadName
String RadName=SrwElecName+SrwMagName+SrwSmpName
String ElecName=SrwElecName+SrwElecType
String MagName=SrwMagName+SrwFieldType
String ObsName=SrwSmpName+SrwSmpType
prompt RadName,SrwPRadName
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
Silent 1						|	Computing the Radiation  ...
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
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwWfrCreate_(); 
		Return;
	endif
endif
if(MagWavePresent==0)
	SrwMagFieldCreate();
	DoAlert 0, SrwPAlertGenMagFieldNeeded;
	Abort;
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 1;
	SrwStartMacrosAfterRadSmp2 *= -1;
	SrwRadSamplDialog(SrwSmpType); 
	Return;
endif

SrwWfrCreate(RadName, ElecName, MagName, ObsName, 1, 1)
end 

//+++++++++++++++++++++++++++++++++++++++
//
//Create Electric Field Wavefront structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreate(RadName, ElecName, MagName, ObsName, ObsNxNzForProp, ObsNxNzSamplFact)
//String RadName=SrwRadName
String RadName=srwUtiTruncString(SrwElecName+SrwMagName+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType
String MagName=SrwMagName+SrwFieldType
String ObsName=SrwSmpName+SrwSmpType
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt RadName,SrwPRadName
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	Computing the Radiation  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwWfrCreate()
		Return;
	endif
endif
if(MagWavePresent==0)
	SrwMagFieldCreate()
	DoAlert 0, SrwPAlertGenMagFieldNeeded
	Abort
endif
if(ObsWavePresent==0)
	if(SrwStartMacrosAfterRadSmp == 0)
		SrwStartMacrosAfterRadSmp = 2
	endif
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

if(ObsNxNzForProp<1)
	ObsNxNzForProp=1
endif
if(ObsNxNzForProp>2)
	ObsNxNzForProp=2
endif
if(ObsNxNzSamplFact<=0.)
	Abort SrwPAlertBadNxNzSamplFact
endif

SrwWfrCrDestroy() // Kills Create Wavefront Panel, if it exists

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]

SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwMagName+SrwSmpName							
endif
SrwRadName=RadName

// Preparing data for srLoop
Make/D/O/N=7 waveprec

//waveprec=str2num($MagName[p+2]) //produces warning with IGOR Pro 6.32
variable j = 0
do
	waveprec[j] =str2num($MagName[j+2])
	j += 1
while(j < 4)

waveprec[4]=str2num($MagName[7])   // Maximum memory [MB] SRW can use (obsolete)
waveprec[5]=1   // TryToApplyNearFieldResidual

waveprec[6]=1 //Do show Progress Indicator by default
if(exists("SrwUtiShowProgrIndic") == 2) //To change default behavior
	waveprec[6]=SrwUtiShowProgrIndic
endif

Make/D/O/N=5 AuxObsTreat; AuxObsTreat[0]=ObsNxNzForProp-1; AuxObsTreat[1]=ObsNxNzSamplFact

SrwWfrPrep(ElecName,ObsName,RadName,0)
RadName += SrwRadType
SrwRadGenTotName=RadName

srLoop($ElecName, $MagName, $ObsName, waveprec, AuxObsTreat, $RadName)

KillWaves/Z  waveprec, AuxObsTreat
end 

//+++++++++++++++++++++++++++++++++++++++
//
//Create Electric Field Wavefront structure from trajectory
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreateFromTrj(RadName,TrjName,ObsName,Mode,Prec,UseDiffLimits,sdep,sfin,ObsNxNzForProp,ObsNxNzSamplFact)
//String RadName=SrwRadName
String RadName=SrwTrjName+SrwSmpName
String TrjName=SrwTrjName+SrwTrjType
String ObsName=SrwSmpName+SrwSmpType
Variable Mode=SrwMode
Variable Prec=SrwPrec
Variable UseDiffLimits=SrwUseDiffRadIntLimits
Variable sdep=SrwSdep
Variable sfin=SrwSfin
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt RadName,SrwPRadName
prompt TrjName,SrwPTrjName2,popup Wavelist("*"+SrwTrjType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt Mode,SrwPMode,popup SrwPOPUPMode
prompt Prec,"Rel. Prec. (for Auto) / Step [m] (for Man.)"
prompt UseDiffLimits,SrwPUseDiffRadIntLimits,popup SrwPOPUPUseDiffRadIntLimits
prompt sdep,"Longitud. Position to Start Integ. [m]"
prompt sfin,"Longitud. Position to Finish Integ. [m]"
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	Computing the Radiation  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

if(cmpstr(TrjName,"_none_")==0)
	Abort SrwPAlertNoTrj
endif
if(cmpstr(ObsName,"_none_")==0)
	Abort SrwPAlertRadSmplNeeded
endif

if(ObsNxNzForProp<1)
	ObsNxNzForProp=1
endif
if(ObsNxNzForProp>2)
	ObsNxNzForProp=2
endif
if(ObsNxNzSamplFact<=0.)
	Abort SrwPAlertBadNxNzSamplFact
endif

SrwWfrCrDestroy() // Kills Create Wavefront Panel, if it exists

SrwRadName=RadName
SrwTrjName=TrjName[0,strlen(TrjName)-strlen(SrwTrjType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwMode=Mode
if(Mode==1)
	SrwRadIntStep=Prec
else
	SrwPrec=Prec
endif
SrwUseDiffRadIntLimits=UseDiffLimits
SrwSdep=sdep
SrwSfin=sfin
SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

if(UseDiffLimits==1) // "No"
	sdep=0; sfin=0
endif

// Preparing data for SR computation routine
Make/D/O/N=6 waveprec
waveprec[0]=Mode-1
waveprec[1]=Prec
waveprec[2]=sdep
waveprec[3]=sfin
waveprec[4]=15000   // Number of point to save
waveprec[5]=1   // TryToApplyNearFieldResidual

Make/D/O/N=5 AuxObsTreat
AuxObsTreat[0]=ObsNxNzForProp-1
AuxObsTreat[1]=ObsNxNzSamplFact

SrwWfrPrep(TrjName,ObsName,RadName,0)
RadName += SrwRadType
SrwRadGenTotName=RadName

srWfrFromTrj($TrjName,$ObsName,waveprec,AuxObsTreat,$RadName)

KillWaves/Z  waveprec, AuxObsTreat
end 

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Electric Field from Trajectory
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreateFromTrj_(RadName,TrjName,ObsName,Mode,Prec,UseDiffLimits,sdep,sfin)
//String RadName=SrwRadName
String RadName=SrwTrjName+SrwSmpName
String TrjName=SrwTrjName+SrwTrjType
String ObsName=SrwSmpName+SrwSmpType
Variable Mode=SrwMode
Variable Prec=SrwPrec
Variable UseDiffLimits=SrwUseDiffRadIntLimits
Variable sdep=SrwSdep
Variable sfin=SrwSfin
prompt RadName,SrwPRadName
prompt TrjName,SrwPTrjName2,popup Wavelist("*"+SrwTrjType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt Mode,SrwPMode,popup SrwPOPUPMode
prompt Prec,"Rel. Prec. (for Auto) / Step [m] (for Man.)"
prompt UseDiffLimits,SrwPUseDiffRadIntLimits,popup SrwPOPUPUseDiffRadIntLimits
prompt sdep,"Longitud. Position to Start Integ. [m]"
prompt sfin,"Longitud. Position to Finish Integ. [m]"
Silent 1						|	Computing the Radiation  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

if(cmpstr(TrjName,"_none_")==0)
	Abort SrwPAlertNoTrj
endif
if(cmpstr(ObsName,"_none_")==0)
	Abort SrwPAlertRadSmplNeeded
endif

SrwWfrCreateFromTrj(RadName,TrjName,ObsName,Mode,Prec,UseDiffLimits,sdep,sfin,1,1)
end 

//+++++++++++++++++++++++++++++++++++++++
//
//Compute CSR Electric Field in frequency domain
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCSRCreate(RadName, ElecName, MagName, ObsName, MethNo, PrecParam, sdep, sfin, ObsNxNzForProp, ObsNxNzSamplFact)
String RadName=srwUtiTruncString(SrwElecName+SrwMagName+SrwSmpName, 27)
String ElecName=srwUtiGetValS("SrwElecLastTot", "", "") //SrwElecLastTot is assumed to have extension
String MagName=SrwMagName+SrwFieldType
String ObsName=SrwSmpName+SrwSmpType
Variable MethNo=srwUtiGetValN("MethNo", 1, "SrwWfrCSRCreate")
Variable PrecParam=srwUtiGetValN("PrecParam", 1, "SrwWfrCSRCreate")
Variable sdep=srwUtiGetValN("sdep", -100, "SrwWfrCSRCreate")
Variable sfin=srwUtiGetValN("sfin", 100, "SrwWfrCSRCreate")
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt RadName,SrwPRadName
prompt ElecName,"Electron Beam structure",popup Wavelist("*"+SrwElecType ,";", "")+Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc","") ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt MethNo,"Method of Longitudinal Integration",popup "Manual" //SrwPOPUPMode
prompt PrecParam,"Prec. Param. (for Auto) / Step [m] (for Man.)"
prompt sdep,"Longitud. Position to Start Integ. [m]"
prompt sfin,"Longitud. Position to Finish Integ. [m]"
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	Computing the Radiation  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwWfrCSRCreate() 
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagFieldCreate()
	DoAlert 0, SrwPAlertGenMagFieldNeeded
	Abort
endif
if(ObsWavePresent==0)
	if(SrwStartMacrosAfterRadSmp == 0)
		SrwStartMacrosAfterRadSmp = 2
	endif
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType) 
	Return
endif

if(ObsNxNzForProp<1)
	ObsNxNzForProp=1
endif
if(ObsNxNzForProp>2)
	ObsNxNzForProp=2
endif
if(ObsNxNzSamplFact<=0.)
	Abort SrwPAlertBadNxNzSamplFact
endif

variable LongitParamWereNotSetUp = 0
if(dimsize($ElecName, 0) < 31)
	LongitParamWereNotSetUp = 1
endif
if(LongitParamWereNotSetUp == 1)
	Abort "CSR can not be computed, because longitudinal bunch parameters were not defined."
endif

srwUtiSetValS("SrwElecLastTot", ElecName, "")  //SrwElecLastTot is assumed to have extension
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
srwUtiSetValN("MethNo", MethNo, "SrwWfrCSRCreate")
srwUtiSetValN("PrecParam", PrecParam, "SrwWfrCSRCreate")
srwUtiSetValN("sdep", sdep, "SrwWfrCSRCreate")
srwUtiSetValN("sfin", sfin, "SrwWfrCSRCreate")

SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwMagName+SrwSmpName							
endif
SrwRadName=RadName

// Preparing data for srWfrCSR
Make/D/O/N=8 waveprec
waveprec = 0
waveprec[0] = MethNo - 1 
waveprec[1] = PrecParam 
waveprec[2] = sdep
waveprec[3] = sfin
waveprec[4] = ObsNxNzForProp-1
waveprec[5] = ObsNxNzSamplFact

//variable/G SrwWfrCSR_MC, SrwWfrCSR_MC_Npart
//if(SrwWfrCSR_MC != 0)
//	waveprec[6] = 1
//	if(SrwWfrCSR_MC_Npart > 0)
//		waveprec[7] = SrwWfrCSR_MC_Npart
//	else
//		waveprec[7] = 100 // default number of macro-particles	
//	endif
//endif

SrwWfrPrep(ElecName,ObsName,RadName,0)
RadName += SrwRadType
SrwRadGenTotName=RadName

srWfrCSR($ElecName, $MagName, $ObsName, waveprec, $RadName)

KillWaves/Z  waveprec
end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute SR Electric Field Panel
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreateDialog(ViewMode)
Variable ViewMode // 0- Reduced version; 1- Full version (for propagation)

SrwWfrCrBufVars(ViewMode)
DoWindow/K SrwWfrCrPanel
SrwWfrCrPanel()

end

//+++++++++++++++++++++++++++++++++++++++
Window SrwWfrCrPanel() : Panel

PauseUpdate; Silent 1		// building window...
NewPanel /W=(410,65,758,497) as "Compute SR Electric Field"
SetDrawLayer UserBack

SrwWfrCrDrawAllContr()

EndMacro

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrDrawAllContr()

Variable VertOffset=-12

// Name
DrawRect 10,22+VertOffset,338,70+VertOffset

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,39+VertOffset,"Name of the Wavefront structure"
SetVariable setvar1WfrCr,pos={30,44+VertOffset},size={290,17},title=" ",fSize=14;
SetVariable setvar1WfrCr,limits={-Inf,Inf,1},value=SrwRadNameBuf

// Input structures
DrawRect 10,75+VertOffset,338,186+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,92+VertOffset,"Input Structures"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,111+VertOffset,"Compute from"
PopupMenu popup1WfrCr,pos={150,94+VertOffset},size={255,21}
PopupMenu popup1WfrCr,value="Magnetic Field;Electron Trajectory",mode=SrwWfrCrFromBuf,proc=SrwWfrCrFromPopupProc

if(SrwWfrCrFromBuf==1) // From Elec. + Mag. Field
	DrawText 30,133+VertOffset,"Electron Beam"
	PopupMenu popup2WfrCr,pos={150,116+VertOffset},size={255,21}
	Variable ElecItemNo=srwWfrCrItemNo(SrwElecNameBuf, SrwElecType)
	SrwElecNameBuf=srwWfrCrBufStrVar(SrwElecType, ElecItemNo)
	PopupMenu popup2WfrCr,value=#"Wavelist(\"*\"+SrwElecType,\";\",\"\")",mode=ElecItemNo,proc=SrwWfrCrElecPopupProc

	DrawText 30,155+VertOffset,"Magnetic Field"
	PopupMenu popup3WfrCr,pos={150,138+VertOffset},size={255,21}
	Variable MagItemNo=srwWfrCrItemNo(SrwMagNameBuf, SrwFieldType)
	SrwMagNameBuf=srwWfrCrBufStrVar(SrwFieldType, MagItemNo)
	PopupMenu popup3WfrCr,value=#"Wavelist(\"*\"+SrwFieldType,\";\",\"\")",mode=MagItemNo,proc=SrwWfrCrMagPopupProc
endif
if(SrwWfrCrFromBuf==2) // Trajectory
	DrawText 30,133+VertOffset,"Electron Trajectory"
	PopupMenu popup2WfrCr,pos={150,116+VertOffset},size={255,21}
	Variable TrjItemNo=srwWfrCrItemNo(SrwTrjNameBuf, SrwTrjType)
	SrwTrjNameBuf=srwWfrCrBufStrVar(SrwTrjType, TrjItemNo)
	PopupMenu popup2WfrCr,value=#"Wavelist(\"*\"+SrwTrjType,\";\",\"\")",mode=TrjItemNo,proc=SrwWfrCrTrjPopupProc
endif

DrawText 30,177+VertOffset,"Radiation Sampling"
PopupMenu popup4WfrCr,pos={150,160+VertOffset},size={255,21}
Variable SmpItemNo=srwWfrCrItemNo(SrwSmpNameBuf, SrwSmpType)
SrwSmpNameBuf=srwWfrCrBufStrVar(SrwSmpType, SmpItemNo)
PopupMenu popup4WfrCr,value=#"Wavelist(\"*\"+SrwSmpType,\";\",\"\")",mode=SmpItemNo,proc=SrwWfrCrSmpPopupProc

// Precision parameters
DrawRect 10,191+VertOffset,338,328+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,208+VertOffset,"Precision Parameters"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 30,228+VertOffset,"Method for Integration"
PopupMenu popup5WfrCr,pos={200,211+VertOffset},size={255,21}
//PopupMenu popup5WfrCr,value="Manual;Auto Undulator;Auto Wiggler",mode=SrwModeBuf,proc=SrwWfrCrModePopupProc
//PopupMenu popup5WfrCr,value="Manual;Auto Undulator;Auto Wiggler;FFT",mode=SrwModeBuf,proc=SrwWfrCrModePopupProc
PopupMenu popup5WfrCr,value="Manual;Auto Undulator;Auto Wiggler",mode=SrwModeBuf,proc=SrwWfrCrModePopupProc

if(SrwModeBuf==1)
	DrawText 30,251+VertOffset,"Longitudinal Step [m]"
	SetVariable setvar2WfrCr,pos={200,234+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar2WfrCr,limits={-Inf,Inf,1},value=SrwRadIntStepBuf
endif
if(SrwModeBuf>1)
	DrawText 30,251+VertOffset,"Relative Precision"
	SetVariable setvar2WfrCr,pos={200,234+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar2WfrCr,limits={-Inf,Inf,1},value=SrwPrecBuf
endif

DrawRect 29,258+VertOffset,320,278+VertOffset
CheckBox check1WfrCr,pos={30,259+VertOffset},size={289,18},value=SrwUseDiffRadIntLimitsBuf-1,proc=SrwWfrCrDiffLimCheckProc
CheckBox check1WfrCr,title="Use Special Limits for Integration"

if(SrwUseDiffRadIntLimitsBuf==2)
	DrawText 30,298+VertOffset,"Initial Position [m]"
	SetVariable setvar3WfrCr,pos={200,281+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar3WfrCr,limits={-Inf,Inf,1},value=SrwSdepBuf

	DrawText 30,318+VertOffset,"Final Position [m]"
	SetVariable setvar4WfrCr,pos={200,301+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar4WfrCr,limits={-Inf,Inf,1},value=SrwSfinBuf
endif

// For Propagation
if(SrwWfrCrViewModeBuf != 0)
	DrawRect 10,333+VertOffset,338,406+VertOffset
	SetDrawEnv fname= "default",fsize=12,fstyle=1
	DrawText 20,350+VertOffset,"For Propagation"

	DrawRect 29,356+VertOffset,320,376+VertOffset
	CheckBox check2WfrCr,pos={30,357+VertOffset},size={289,18},value=SrwSmpNxNzForPropBuf-1,proc=SrwWfrCrNxNzPropCheckProc
	CheckBox check2WfrCr,title=SrwPSmpNxNzForProp

	if(SrwSmpNxNzForPropBuf==2)
		DrawText 30,396+VertOffset,SrwPSmpNxNzSamplFact
		SetVariable setvar5WfrCr,pos={200,379+VertOffset},size={120,14},title=" ",fSize=14
		SetVariable setvar5WfrCr,limits={-Inf,Inf,1},value=SrwSmpNxNzSamplFactBuf
	endif
endif

// OK-Cancel-Help
Button button1WfrCr,pos={38,415+VertOffset},size={70,20},proc=SrwWfrCrCancelButtonProc,title=SrwPPanelCancelButton
Button button2WfrCr,pos={138,415+VertOffset},size={70,20},proc=SrwWfrCrOKButtonProc,title=SrwPPanelOKButton
Button button3WfrCr,pos={238,415+VertOffset},size={70,20},proc=SrwWfrCrHelpButtonProc,title=SrwPPanelHelpButton

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrKillAllContr()

if(CheckName("SrwWfrCrPanel",9) != 0)
	KillControl setvar1WfrCr
	KillControl setvar2WfrCr
	KillControl setvar3WfrCr
	KillControl setvar4WfrCr
	KillControl setvar5WfrCr
	KillControl popup1WfrCr
	KillControl popup2WfrCr
	KillControl popup3WfrCr
	KillControl popup4WfrCr
	KillControl popup5WfrCr
	KillControl check1WfrCr
	KillControl check2WfrCr
endif

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrBufVars(ViewMode)
Variable ViewMode // 0- Reduced version; 1- Full version (for propagation)

Variable/G SrwWfrCrFromBuf=1
Variable/G SrwWfrCrViewModeBuf=ViewMode
Variable/G SrwWfrCrPrecOrStepBuf

Variable/G SrwModeBuf=SrwMode
Variable/G SrwRadIntStepBuf=SrwRadIntStep
Variable/G SrwPrecBuf=SrwPrec
Variable/G SrwUseDiffRadIntLimitsBuf=SrwUseDiffRadIntLimits
Variable/G SrwSdepBuf=SrwSdep
Variable/G SrwSfinBuf=SrwSfin
Variable/G SrwMaxPtsToSaveBuf=SrwMaxPtsToSave
Variable/G SrwSmpNxNzForPropBuf=SrwSmpNxNzForProp
Variable/G SrwSmpNxNzSamplFactBuf=SrwSmpNxNzSamplFact

String/G SrwRadNameBuf=SrwRadName
String/G SrwElecNameBuf=SrwElecName
String/G SrwMagNameBuf=SrwMagName
String/G SrwTrjNameBuf=SrwTrjName
String/G SrwSmpNameBuf=SrwSmpName

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrKillBufVars()

KillVariables/Z SrwWfrCrFromBuf, SrwWfrCrViewModeBuf, SrwWfrCrPrecOrStepBuf
KillVariables/Z SrwModeBuf, SrwRadIntStepBuf, SrwPrecBuf, SrwUseDiffRadIntLimitsBuf, SrwSdepBuf, SrwSfinBuf, SrwMaxPtsToSaveBuf
KillVariables/Z SrwSmpNxNzForPropBuf, SrwSmpNxNzSamplFactBuf
KillStrings/Z SrwRadNameBuf, SrwElecNameBuf, SrwMagNameBuf, SrwTrjNameBuf, SrwSmpNameBuf

end

//+++++++++++++++++++++++++++++++++++++++
function srwWfrCrItemNo(NameBuf, Type)
String NameBuf, Type

String AllWaves=Wavelist("*"+Type,";","")
String NameTot=NameBuf+Type
Variable ItemNo=sRWaveListItemNo(AllWaves, ";", NameTot)
if(ItemNo<1)
	ItemNo=1
endif
return ItemNo
end

//+++++++++++++++++++++++++++++++++++++++
function/S srwWfrCrBufStrVar(Type, ItemNo)
String Type
Variable ItemNo

String AllWaves=Wavelist("*"+Type,";","")
String NameTot=sRWaveListItemName(AllWaves, ";", ItemNo)
String NameBuf
if(cmpstr(NameTot,"") != 0)
	NameBuf=NameTot[0,strlen(NameTot)-strlen(Type)-1]
else
	NameBuf=""
endif
return NameBuf
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrFromPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwWfrCrFromBuf = popNum

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrElecPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwElecNameBuf = popStr[0,strlen(popStr)-strlen(SrwElecType)-1]

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrMagPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwMagNameBuf = popStr[0,strlen(popStr)-strlen(SrwFieldType)-1]

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrTrjPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwTrjNameBuf = popStr[0,strlen(popStr)-strlen(SrwTrjType)-1]

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrSmpPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwSmpNameBuf = popStr[0,strlen(popStr)-strlen(SrwSmpType)-1]

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrModePopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwModeBuf = popNum

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrDiffLimCheckProc(ctrlName,checked) : CheckBoxControl
String ctrlName
Variable checked			// 1 if checked, 0 if not

SrwUseDiffRadIntLimitsBuf = checked+1

if(SrwUseDiffRadIntLimitsBuf==2)
	Variable ItemNo
	String StructName = "", AuxCoreName = ""

	if(SrwWfrCrFromBuf==1) // From Elec. + Mag. Field
		ItemNo=srwWfrCrItemNo(SrwMagNameBuf, SrwFieldType)
		AuxCoreName=srwWfrCrBufStrVar(SrwFieldType, ItemNo)
		if(strlen(AuxCoreName) > 0)
			StructName=AuxCoreName+SrwFieldType
			StructName=$StructName[1]
		endif
	endif
	if(SrwWfrCrFromBuf==2) // Trajectory
		ItemNo=srwWfrCrItemNo(SrwTrjNameBuf, SrwTrjType)
		AuxCoreName=srwWfrCrBufStrVar(SrwTrjType, ItemNo)
		if(strlen(AuxCoreName) > 0)
			StructName=AuxCoreName+SrwTrjType
			StructName=$StructName[0]
		endif
	endif

	if(cmpstr(StructName,"") != 0)
		Variable sStart = DimOffset($StructName, 0)
		Variable sStep = DimDelta($StructName, 0)
		Variable Np = DimSize($StructName, 0)
		Variable sFin = sStart + sStep*(Np - 1)
		
		if(SrwSdepBuf==SrwSfinBuf)
			SrwSdepBuf = sStart
			SrwSfinBuf = sFin
		endif
		if(SrwSdepBuf < sStart)
			SrwSdepBuf = sStart
		endif
		if(SrwSfinBuf > sFin)
			SrwSfinBuf = sFin
		endif
	endif
endif

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrNxNzPropCheckProc(ctrlName,checked) : CheckBoxControl
String ctrlName
Variable checked			// 1 if checked, 0 if not

SrwSmpNxNzForPropBuf = checked+1

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwWfrCrKillAllContr()
SrwWfrCrDrawAllContr()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrCancelButtonProc(ctrlName) : ButtonControl
String ctrlName

SrwWfrCrDestroy()
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrDestroy()

SrwWfrCrKillAllContr()
DoWindow/K SrwWfrCrPanel
SrwWfrCrKillBufVars()
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrValidateGlobVars()

if(strlen(SrwRadNameBuf)==0)
	Abort SrwPAlertBadName
endif
if(strlen(SrwRadNameBuf)>28)
	Abort SrwPAlertTooLongName
endif

if(SrwWfrCrFromBuf==1)
	if(strlen(SrwElecNameBuf)==0)
		Abort SrwPAlertNoElec
	endif
	if(strlen(SrwMagNameBuf)==0)
		Abort SrwPAlertNoMag
	endif
endif
if(SrwWfrCrFromBuf==2)
	if(strlen(SrwTrjNameBuf)==0)
		Abort SrwPAlertNoTrj
	endif
endif
if(strlen(SrwSmpNameBuf)==0)
	Abort SrwPAlertNoSmp
endif

if(SrwModeBuf==1)
	if(SrwRadIntStepBuf <= 0)
		Abort SrwPAlertPrecStep
	endif
	
	SrwWfrCrPrecOrStepBuf=SrwRadIntStepBuf
else
	if(SrwPrecBuf <= 0)
		Abort SrwPAlertPrecRel
	endif
	
	SrwWfrCrPrecOrStepBuf=SrwPrecBuf
endif

if(SrwUseDiffRadIntLimitsBuf==2) // use
	if(SrwSdepBuf >= SrwSfinBuf)
		Abort SrwPAlertPrecIntLim
	endif
else 
	SrwSdepBuf = 0
	SrwSfinBuf = 0
endif

if(SrwSmpNxNzForPropBuf==2) // yes
	if(SrwSmpNxNzSamplFactBuf <= 0)
		Abort SrwPAlertBadNxNzSamplFact
	endif
endif

SrwMaxPtsToSaveBuf = 10000

SrwElecNameBuf+=SrwElecType
SrwMagNameBuf+=SrwFieldType
SrwTrjNameBuf+=SrwTrjType
SrwSmpNameBuf+=SrwSmpType

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrOKButtonProc(ctrlName) : ButtonControl
String ctrlName

SrwWfrCrValidateGlobVars()
String ComLineStr, ComLineStr1, ComLineStr2

if(SrwWfrCrFromBuf==1) // From Elec. + Mag. Field

	sprintf ComLineStr1, "SrwMagPrec(\"%s\",%g,%g,%g,%g,%g,%g,%g)",SrwMagNameBuf,SrwModeBuf,SrwRadIntStepBuf,SrwPrecBuf,SrwMaxPtsToSaveBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf
	if(SrwWfrCrViewModeBuf==0) // Not for propagation
	
		sprintf ComLineStr2, "SrwWfrCreate_(\"%s\",\"%s\",\"%s\",\"%s\")",SrwRadNameBuf,SrwElecNameBuf,SrwMagNameBuf,SrwSmpNameBuf
		ComLineStr=ComLineStr1+";"+ComLineStr2
		print ComLineStr
		SrwMagPrec(SrwMagNameBuf,SrwModeBuf,SrwRadIntStepBuf,SrwPrecBuf,SrwMaxPtsToSaveBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf)
		SrwWfrCreate_(SrwRadNameBuf,SrwElecNameBuf,SrwMagNameBuf,SrwSmpNameBuf)
		Return
		
	else // For Propagation
	
		sprintf ComLineStr2, "SrwWfrCreate(\"%s\",\"%s\",\"%s\",\"%s\",%g,%g)",SrwRadNameBuf,SrwElecNameBuf,SrwMagNameBuf,SrwSmpNameBuf,SrwSmpNxNzForPropBuf,SrwSmpNxNzSamplFactBuf
		ComLineStr=ComLineStr1+";"+ComLineStr2
		print ComLineStr
		SrwMagPrec(SrwMagNameBuf,SrwModeBuf,SrwRadIntStepBuf,SrwPrecBuf,SrwMaxPtsToSaveBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf)
		SrwWfrCreate(SrwRadNameBuf,SrwElecNameBuf,SrwMagNameBuf,SrwSmpNameBuf,SrwSmpNxNzForPropBuf,SrwSmpNxNzSamplFactBuf)
		Return
	
	endif
endif
if(SrwWfrCrFromBuf==2) // From Trajectory
	if(SrwWfrCrViewModeBuf==0) // Not for propagation
	
		sprintf ComLineStr, "SrwWfrCreateFromTrj_(\"%s\",\"%s\",\"%s\",%g,%g,%g,%g,%g)",SrwRadNameBuf,SrwTrjNameBuf,SrwSmpNameBuf,SrwModeBuf,SrwWfrCrPrecOrStepBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf
		print ComLineStr
		SrwWfrCreateFromTrj_(SrwRadNameBuf,SrwTrjNameBuf,SrwSmpNameBuf,SrwModeBuf,SrwWfrCrPrecOrStepBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf)
		Return

	else // For Propagation
	
		sprintf ComLineStr, "SrwWfrCreateFromTrj(\"%s\",\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g)",SrwRadNameBuf,SrwTrjNameBuf,SrwSmpNameBuf,SrwModeBuf,SrwWfrCrPrecOrStepBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf,SrwSmpNxNzForPropBuf,SrwSmpNxNzSamplFactBuf
		print ComLineStr
		SrwWfrCreateFromTrj(SrwRadNameBuf,SrwTrjNameBuf,SrwSmpNameBuf,SrwModeBuf,SrwWfrCrPrecOrStepBuf,SrwUseDiffRadIntLimitsBuf,SrwSdepBuf,SrwSfinBuf,SrwSmpNxNzForPropBuf,SrwSmpNxNzSamplFactBuf)
		Return

	endif
endif

SrwWfrCrDestroy()

end

//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCrHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwWfrCreateDialog     ")
end
