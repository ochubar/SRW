
//+++++++++++++++++++++++++++++++++++++++
//
//Create a Periodic Magnetic field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerCreate2D(UndName,Period,Kz,Kx,Length,Ph0x,UndKind,UndTaper,UndPhaseShift)
string UndName=SrwUndName
Variable Period=SrwPeriod
Variable Kz=SrwKz
Variable Kx=SrwKx
Variable Length=SrwLength
Variable Ph0x=SrwPh0x
Variable UndKind=SrwUndKind
Variable UndTaper=SrwUndTaper
Variable UndPhaseShift=SrwUndOptKlystPhaseShift
prompt UndName,SrwPUndName
prompt Period,SrwPPeriod
prompt Kz,SrwPKz
prompt Kx,SrwPKx
prompt Length,SrwPLength
prompt Ph0x,SrwPPh0x
prompt UndKind,SrwPUndKind,popup "Conventional;Tapered;Optical Klystron;Infinite"
prompt UndTaper,SrwPUndTaper
prompt UndPhaseShift,SrwPUndOptKlystPhaseShift
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(Length <= 0.)
	Abort SrwPAlertUndLength
endif
if((Kx < 0.) %| (Kz < 0.))
	Abort SrwPAlertUndKK
endif

SrwKz=Kz
SrwKx=Kx
SrwPh0x=Ph0x

if(Kz > 0.)
	SrwMagPerCreate(UndName,Period,1,Kz,Length,UndKind,UndTaper,UndPhaseShift)
	SrwUndFieldType = 1
	if(Kx > 0.)
		UndName+=SrwUndType
		SrwMagPerAddHarm(UndName,1,2,Kx,Ph0x)
		SrwUndFieldType = 3
	endif
else
	if(Kx > 0.)
		SrwMagPerCreate(UndName,Period,2,Kx,Length,UndKind,UndTaper,UndPhaseShift)
		SrwUndFieldType = 2
	else
		Abort SrwPAlertUndKK
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create a Periodic Magnetic field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerCreate2D_(UndName,Period,Kz,Kx,Length,Ph0x)
string UndName=SrwUndName
Variable Period=SrwPeriod
Variable Kz=SrwKz
Variable Kx=SrwKx
Variable Length=SrwLength
Variable Ph0x=SrwPh0x
prompt UndName,SrwPUndName
prompt Period,SrwPPeriod
prompt Kz,SrwPKz
prompt Kx,SrwPKx
prompt Length,SrwPLength
prompt Ph0x,SrwPPh0x
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(Length <= 0.)
	Abort SrwPAlertUndLength
endif
if((Kx < 0.) %| (Kz < 0.))
	Abort SrwPAlertUndKK
endif

SrwKz=Kz
SrwKx=Kx
SrwPh0x=Ph0x

if(Kz > 0.)
	SrwMagPerCreate_(UndName,Period,1,Kz,Length)
	if(Kx > 0.)
		UndName+=SrwUndType
		SrwMagPerAddHarm(UndName,1,2,Kx,Ph0x)
	endif
else
	if(Kx > 0.)
		SrwMagPerCreate_(UndName,Period,2,Kx,Length)
	else
		Abort SrwPAlertUndKK
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create a Periodic Magnetic field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerCreate(UndName,Period,Plane,KK,Length,UndKind,UndTaper,UndPhaseShift)
string UndName=SrwUndName
Variable Period=SrwPeriod
Variable Plane=SrwPlane
Variable KK=SrwKK
Variable Length=SrwLength
Variable UndKind=SrwUndKind
Variable UndTaper=SrwUndTaper
Variable UndPhaseShift=SrwUndOptKlystPhaseShift
prompt UndName,SrwPUndName
prompt Period,SrwPPeriod
prompt Plane,SrwPPlane,popup "Vertical;Horizontal"
prompt KK,SrwPKK
prompt Length,SrwPLength
prompt UndKind,SrwPUndKind,popup "Conventional;Tapered;Optical Klystron;Infinite"
prompt UndTaper,SrwPUndTaper
prompt UndPhaseShift,SrwPUndOptKlystPhaseShift
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(KK <= 0.)
	Abort SrwPAlertUndKK
endif
if(Length <= 0.)
	Abort SrwPAlertUndLength
endif

SrwUndName=UndName
SrwMagGenTotName=UndName+SrwUndType

SrwPeriod=Period
SrwLength=Length
SrwUndKind=UndKind
SrwUndTaper=UndTaper
SrwUndOptKlystPhaseShift=UndPhaseShift
SrwPlane=Plane
SrwKK=KK

Period*=0.001
UndName+=SrwUndType
Make/T/O/N=6  $UndName
//SetScale d -1E+23, 1E+23, SrwMagType_Periodic, $UndName

$UndName[0]=num2str(Period)
$UndName[1]=num2str(Length)
$UndName[2]=num2str(UndKind)

$UndName[3]=""
if(UndKind == 2) // Tapered
	$UndName[3]=num2str(UndTaper)
endif
if(UndKind == 3) // Opt. Klystron
	$UndName[3]=num2str(UndPhaseShift)
endif

$UndName[4]=num2str(9.5/Period/1000)
$UndName[5]=num2str(0)
SrwMagPerAddHarm(UndName,1,Plane,KK,0)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create a Periodic Magnetic field structure (Trancated)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerCreate_(UndName,Period,Plane,KK,Length)
string UndName=SrwUndName
Variable Period=SrwPeriod
Variable Plane=SrwPlane
Variable KK=SrwKK
Variable Length=SrwLength
prompt UndName,SrwPUndName
prompt Period,SrwPPeriod
prompt Plane,SrwPPlane,popup "Vertical;Horizontal"
prompt KK,SrwPKK
prompt Length,SrwPLength
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(KK <= 0.)
	Abort SrwPAlertUndKK
endif
if(Length <= 0.)
	Abort SrwPAlertUndLength
endif

SrwUndName=UndName
SrwMagGenTotName=UndName+SrwUndType

SrwPeriod=Period
SrwLength=Length
SrwUndKind=1
SrwPlane=Plane
SrwKK=KK

Period*=0.001
UndName+=SrwUndType
Make/T/O/N=6  $UndName
//SetScale d -1E+23, 1E+23, SrwMagType_Periodic, $UndName

$UndName[0]=num2str(Period)
$UndName[1]=num2str(Length)
$UndName[2]=num2str(SrwUndKind)

$UndName[3]=""

$UndName[4]=num2str(9.5/Period/1000)
$UndName[5]=num2str(0)
SrwMagPerAddHarm(UndName,1,Plane,KK,0)

end

//+++++++++++++++++++++++++++++++++++++++
//
// Append harmonic to a Periodic Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerAddHarm(UndName,Harm,Plane,KK,Phase)
string UndName=SrwUndName+SrwUndType
Variable Harm=SrwHarm
Variable Plane=SrwPlane
Variable KK=SrwKK
Variable Phase=SrwPhase
prompt UndName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt Harm,SrwPHarm
prompt Plane,SrwPPlane,popup "Vertical;Horizontal"
prompt KK,"Deflection Parameter (0.0934 x B x <Period of Undulator> / <Harmonic Number>)"
prompt Phase,SrwPPhase
Silent 1						|	...
PauseUpdate

if(cmpstr(UndName,"_none_")==0)
	SrwMagPerCreate()
	SrwMagPerAddHarm()
	Abort
endif

// Validation of parameters
if((abs(round(Harm) - Harm) > 1.E-08) %| (Harm <= 0))
	Abort SrwPAlertUndHarm
endif
if(KK <= 0.)
	Abort SrwPAlertUndKK
endif

SrwUndName=UndName[0,strlen(UndName)-strlen(SrwUndType)-1]
SrwHarm=Harm
SrwPlane=Plane
SrwPhase=Phase
SrwKK=KK

Variable n=dimsize($undname,0)+1

redimension/N=(n)  $UndName
$UndName[5]=num2str(n-6)
String HarmName=SrwUndName+num2str(n-6)+SrwUndHarmType
Make/O/D/N=(4) $HarmName
$UndName[n-1]=HarmName
$HarmName[0]=Harm
$HarmName[1]=Plane
$HarmName[2]=KK
$HarmName[3]=Phase

Variable Period=str2num($UndName[0])*1000
Variable k2=1/(str2num($UndName[4])/9.5*Period)
Variable En=9.5/Period/(k2+kk*kk/harm/2)

$UndName[4]=num2str(En)

if(SrwAllowPrintingExtraInfo==1)
	print "Bmax = ",kk/(0.0934*Period)," [T]"
	print "1+ K^2/2 = ",k2+kk*kk/harm/2
	print "Fundamental = ", En, " [keV/GeV^2]"
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Compute field and deflection parameter of a pure permanent magnet undulator
//
//+++++++++++++++++++++++++++++++++++++++
function srwMagPerPPM_B(Period,Gap,Height,Air,MagPer,Br)
variable Period //=SrwPeriod
variable Gap //=SrwGap
variable Air //=SrwAir
variable Height //=SrwHeight
variable MagPer //=SrwMagPer
variable Br //=SrwBr

SVAR SrwPAlertUndPeriod, SrwPAlertUndGap, SrwPAlertUndHeight, SrwPAlertUndAir
NVAR SrwPeriod, SrwGap, SrwAir, SrwMagPer, SrwBr, SrwHeight

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(Gap <= 0.)
	Abort SrwPAlertUndGap
endif
if(Height <= 0.)
	Abort SrwPAlertUndHeight
endif
if(Air <= 0.)
	Abort SrwPAlertUndAir
endif

SrwPeriod=Period
SrwGap=Gap
SrwAir=Air
SrwMagPer=MagPer
SrwBr=Br
SrwHeight=Height

variable Bc=2*Br*exp(-Pi*Gap/Period)*sin(Pi/MagPer)/(Pi/MagPer)
Bc*=(1-exp(-2*Pi*Height/Period))
Bc*= cos(Pi/2*MagPer/Period*Air)
return Bc
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerPPM(Period,Gap,Height,Air,MagPer,Br)
Variable Period=SrwPeriod
Variable Gap=SrwGap
Variable Air=SrwAir
Variable Height=SrwHeight
Variable MagPer=SrwMagPer
Variable Br=SrwBr
prompt Period,SrwPPeriod
prompt Gap,SrwPGap
prompt Height,SrwPHeight
prompt Air,SrwPAir
prompt MagPer,SrwPMagPer
prompt Br,SrwPBr
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(Gap <= 0.)
	Abort SrwPAlertUndGap
endif
if(Height <= 0.)
	Abort SrwPAlertUndHeight
endif
if(Air <= 0.)
	Abort SrwPAlertUndAir
endif

SrwPeriod=Period
SrwGap=Gap
SrwAir=Air
SrwMagPer=MagPer
SrwBr=Br
SrwHeight=Height

variable Bc = srwMagPerPPM_B(Period,Gap,Height,Air,MagPer,Br)
variable KK=0.0934*Period*Bc

if(SrwAllowPrintingExtraInfo==1)
	Print "Peak Field =",Bc," [T]"
	Print "K =             ",KK
endif

SrwKK=srRound(KK,3)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Duplicates Periodic Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerDupl(MagName,Name)
String MagName=SrwUndName+SrwUndType
String  Name=SrwUndName+"d"
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType,";","")
prompt Name,"Name of the Duplicated Periodic Field structure"
Silent 1						|	Duplicating the field structure  ....
PauseUpdate

SrwUndName=Name

Name += SrwUndType
duplicate/O $MagName $Name

String FullOldHarmName, MainOldHarmName, MainNewHarmName, FullNewHarmName
Variable AmOfHarm = str2num($MagName[5])
Variable i=0
do
	FullOldHarmName = $MagName[6+i]
	MainOldHarmName = FullOldHarmName[0,strlen(FullOldHarmName)-strlen(SrwUndHarmType)-2]
	//MainNewHarmName = MainOldHarmName + "d"
	MainNewHarmName = SrwUndName
	FullNewHarmName = MainNewHarmName + FullOldHarmName[strlen(FullOldHarmName)-strlen(SrwUndHarmType)-1] + SrwUndHarmType

	duplicate/O $FullOldHarmName $FullNewHarmName
	$Name[6+i] = FullNewHarmName
	
	i += 1
while(i<AmOfHarm)

end

//+++++++++++++++++++++++++++++++++++++++
//
// Periodic Magnetic Field Panel
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPerCreateDialog()

SrwMagPerCreateBufVars()
DoWindow/K SrwMagPerPanel
SrwMagPerPanel()

end

//+++++++++++++++++++++++++++++++++++++++
Window SrwMagPerPanel() : Panel

PauseUpdate; Silent 1		// building window...
NewPanel /W=(410,65,748,414) as "Periodic Magnetic Field"
SetDrawLayer UserBack

SrwMagPerDrawAllContr()

EndMacro

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerDrawAllContr()

Variable VertOffset=0

String AllMagPerWaves=Wavelist("*"+SrwUndType,";","")
String UndNameTot = SrwUndName+SrwUndType
Variable ItemNo=sRWaveListItemNo(AllMagPerWaves, ";", UndNameTot)
//if(ItemNo == 0)
//	ItemNo=1
//endif

// Existing structures
DrawRect 10,16+VertOffset,328,71+VertOffset

SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,40+VertOffset,"Existing structures"
PopupMenu popup0MagPer,pos={140,22+VertOffset},size={255,21}
PopupMenu popup0MagPer,value=#"Wavelist(\"*\"+SrwUndType,\";\",\"\")",mode=ItemNo,proc=SrwMagPerSelectPopupProc

// Name
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,62+VertOffset,"Name"
SetVariable setvar0MagPer,pos={140,44+VertOffset},size={180,17},title=" ",fSize=14;
SetVariable setvar0MagPer,limits={-Inf,Inf,1},value=SrwUndName

// Geometry
DrawRect 10,76+VertOffset,328,144+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,96+VertOffset,"Geometry"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 50,115+VertOffset,SrwPPeriod
SetVariable setvar1MagPer,pos={200,97+VertOffset},size={120,14},title=" ",fSize=14
SetVariable setvar1MagPer,limits={-Inf,Inf,1},value=SrwPeriod

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText 50,135+VertOffset,SrwPLength
SetVariable setvar2MagPer,pos={200,117+VertOffset},size={120,14},title=" ",fSize=14
SetVariable setvar2MagPer,limits={-Inf,Inf,1},value=SrwLength,proc=SrwLengthSetVarProc

// Field
DrawRect 10,149+VertOffset,328,244+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,173+VertOffset,"Field"

PopupMenu popup1MagPer,pos={200,155+VertOffset},size={255,21}
PopupMenu popup1MagPer,value="Vertical;Horizontal;Ellipsoidal",mode=SrwUndFieldType,proc=SrwMagPerFieldPopupProc

if((SrwUndFieldType==1) %| (SrwUndFieldType==3))
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 50,195+VertOffset,SrwPKz
	SetVariable setvar3MagPer,pos={200,177+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar3MagPer,limits={-Inf,Inf,1},value=SrwKz
endif
if((SrwUndFieldType==2) %| (SrwUndFieldType==3))
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 50,215+VertOffset,SrwPKx
	SetVariable setvar4MagPer,pos={200,197+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar4MagPer,limits={-Inf,Inf,1},value=SrwKx
endif
if(SrwUndFieldType==3)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 50,235+VertOffset,SrwPPh0x
	SetVariable setvar5MagPer,pos={200,217+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar5MagPer,limits={-Inf,Inf,1},value=SrwPh0x
endif

// Type
DrawRect 10,249+VertOffset,328,304+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText 20,273+VertOffset,"Type"

PopupMenu popup2MagPer,pos={200,255+VertOffset},size={255,21}
PopupMenu popup2MagPer,value="Conventional;Tapered;Optical Klystron;Infinite",mode=SrwUndKind,proc=SrwMagPerKindPopupProc

if(SrwUndKind==2)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 50,295+VertOffset,"Taper param. = N dE/E"
	SetVariable setvar6MagPer,pos={200,277+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar6MagPer,limits={-Inf,Inf,1},value=SrwUndTaper
endif
if(SrwUndKind==3)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText 50,295+VertOffset,"Phase shift param. = Nd/N"
	SetVariable setvar7MagPer,pos={200,277+VertOffset},size={120,14},title=" ",fSize=14
	SetVariable setvar7MagPer,limits={-Inf,Inf,1},value=SrwUndOptKlystPhaseShift
endif
if(SrwUndKind==4)
	String Str1="Spectral Flux as for Undulator of  ", Str2=" m  Length"
	String StrTot=Str1+num2str(SrwLength)+Str2
	DrawText 50,298+VertOffset,StrTot
endif

// OK-Cancel-Help
Button button0MagPer,pos={35,317+VertOffset},size={70,20},proc=SrwMagPerCancelButtonProc,title=SrwPPanelCancelButton
Button button1MagPer,pos={135,317+VertOffset},size={70,20},proc=SrwMagPerOKButtonProc,title=SrwPPanelOKButton
Button button2MagPer,pos={235,317+VertOffset},size={70,20},proc=SrwMagPerHelpButtonProc,title=SrwPPanelHelpButton

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerKillAllContr()

KillControl popup0MagPer
KillControl popup1MagPer
KillControl popup2MagPer
KillControl setvar0MagPer
KillControl setvar1MagPer
KillControl setvar2MagPer
KillControl setvar3MagPer
KillControl setvar4MagPer
KillControl setvar5MagPer
KillControl setvar6MagPer
KillControl setvar7MagPer

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerSelectPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwUndName=popStr[0,strlen(popStr)-strlen(SrwUndType)-1]
SrwMagPerSetupGlobVars(popStr)

SrwMagPerUpdatePanel()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerSetupGlobVars(UndName)
String UndName

SrwUndName=UndName[0,strlen(UndName)-strlen(SrwUndType)-1]

SrwPeriod=1000.*str2num($UndName[0])
SrwLength=str2num($UndName[1])
SrwUndKind=str2num($UndName[2])
if(SrwUndKind==2) // Tapered
	SrwUndTaper=str2num($UndName[3])
endif
if(SrwUndKind==3) // Opt. Klystron
	SrwUndOptKlystPhaseShift=str2num($UndName[3])
endif

Variable AmOfHarm=str2num($UndName[5]), i=0
Variable Kz=0, Kx=0, Phz=0, Phx=0
String HarmName
do
	HarmName=$UndName[6+i]
	if($HarmName[1]==1)
		if(Kz==0)
			Kz=$HarmName[2]
			Phz=$HarmName[3]
		endif
	endif
	if($HarmName[1]==2)
		if(Kx==0)
			Kx=$HarmName[2]
			Phx=$HarmName[3]
		endif
	endif
	i+=1
while(i<AmOfHarm)

if(Kz != 0)
	if(Kx != 0)
		SrwUndFieldType=3
	else
		SrwUndFieldType=1
	endif
else
	if(Kx != 0)
		SrwUndFieldType=2
	endif
endif
SrwKz=Kz
SrwKx=Kx
SrwPh0x=Phx-Phz

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerValidateGlobVars()

// Validation of parameters

if(cmpstr(SrwUndName,"")==0)
	Abort SrwPAlertMagPerBad
endif
if(strlen(SrwUndName)>28)
	Abort SrwPAlertTooLongName
endif

if(SrwPeriod <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(SrwLength <= 0.)
	Abort SrwPAlertUndLength
endif
if(SrwLength < 0.001*SrwPeriod)
	Abort SrwPAlertUndLengthSm
endif

if((SrwKz < 0.) %| (SrwKx < 0.))
	Abort SrwPAlertUndKK
endif

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerFieldPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwUndFieldType=popNum
SrwMagPerUpdatePanel()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerKindPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwUndKind=popNum
SrwMagPerUpdatePanel()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwLengthSetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
String ctrlName
Variable varNum	// value of variable as number
String varStr		// value of variable as string
String varName	// name of variable

SrwLength=varNum
SrwMagPerUpdatePanel()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerCreateBufVars()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerKillBufVars()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerUpdatePanel()
SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwMagPerKillAllContr()
SrwMagPerDrawAllContr()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerDestroy()
SrwMagPerKillAllContr()
DoWindow/K SrwMagPerPanel
SrwMagPerKillBufVars()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerCancelButtonProc(ctrlName) : ButtonControl
String ctrlName

SrwMagPerDestroy()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerOKButtonProc(ctrlName) : ButtonControl
String ctrlName
String ComLineStr

// Validate parameters
SrwMagPerValidateGlobVars()

if(SrwUndFieldType==1) // Vert.
sprintf ComLineStr, "SrwMagPerCreate2D(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwUndName, SrwPeriod, SrwKz, 0, SrwLength, 0, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift
print ComLineStr
SrwMagPerCreate2D(SrwUndName, SrwPeriod, SrwKz, 0, SrwLength, 0, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift)
endif
if(SrwUndFieldType==2) // Hor.
sprintf ComLineStr, "SrwMagPerCreate2D(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwUndName, SrwPeriod, 0, SrwKx, SrwLength, 0, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift
print ComLineStr
SrwMagPerCreate2D(SrwUndName, SrwPeriod, 0, SrwKx, SrwLength, 0, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift)
endif
if(SrwUndFieldType==3) // Vert.+Hor.
sprintf ComLineStr, "SrwMagPerCreate2D(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwUndName, SrwPeriod, SrwKz, SrwKx, SrwLength, SrwPh0x, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift
print ComLineStr
SrwMagPerCreate2D(SrwUndName, SrwPeriod, SrwKz, SrwKx, SrwLength, SrwPh0x, SrwUndKind, SrwUndTaper, SrwUndOptKlystPhaseShift)
endif

SrwMagPerDestroy()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagPerHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwMagPerCreateDialog     ")
End
