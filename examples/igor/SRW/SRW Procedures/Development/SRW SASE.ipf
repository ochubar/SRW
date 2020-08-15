
//+++++++++++++++++++++++++++++++++++++++
//
// SASE calculations mainly based on the code GENESIS 3D
// converted to C
//
//+++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++
//
//Initialization of relevant globals for SASE and relevant
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEInit()

String/G SrwMagContType=SrwSeparator+"mgg"		// general magnetic field structure

String/G SrwRadSASEType=SrwSeparator+"rss"		// sase rad. structure (temporary?)
String/G SrwPrecSASEType=SrwSeparator+"pss"	// sase precision structure
String/G SrwCntrlSASEType=SrwSeparator+"css"	// sase control output structure
String/G SrwElecDistrType=SrwSeparator+"ebd"

String/G SrwMagContName="M"
String/G SrwMagOptName="L"
String/G SrwRadSASEName="InR"
String/G SrwPrecSASEName="P"
String/G SrwCntrlSASEName="C"

//SrwElecSASEExtra(ElecName,TypeTransv,TypeLong,SigLong,ShotNoise)
Variable/G SrwElecTypeTransv=2
Variable/G SrwElecTypeLong=1
Variable/G SrwElecSigLong=10
Variable/G SrwElecShotNoise=1

//SrwSASEUndCreate(UndName,Period,K,PlanOrEl,FocType,xkx,NperInSec,Nsec,FldErrType,FldErrRMS)
SrwKz=1.26714 //0.896*SQRT(2)
Variable/G SrwMagSASEPlanOrEl=1
Variable/G SrwMagSASEFocType=1
Variable/G SrwMagSASExkx=0.
Variable/G SrwMagSASExky=1.
Variable/G SrwMagSASESecLen=4.477 //m
Variable/G SrwMagSASENsec=3
Variable/G SrwMagSASEIntBwSec=0.15 //m
Variable/G SrwMagSASEFldErrType=1
Variable/G SrwMagSASEFldErrRMS=0

//SrwSASEFODOAdd(MagName,Drift,QuadF,QuadD,LenF,LenD,Dx,Dz,F1St)
Variable/G SrwMagSASEDrift=0.341
Variable/G SrwMagSASEQuadF=18.3
Variable/G SrwMagSASEQuadD=18.3
Variable/G SrwMagSASELenF=0.137
Variable/G SrwMagSASELenD=0.137
Variable/G SrwMagSASEDx=0
Variable/G SrwMagSASEDz=0
Variable/G SrwMagSASEF1St=0.06825

//SrwSASEUndTaperAdd(MagName,TaperType,TaperStart,TaperVal)
Variable/G SrwMagSASETaperType=1
Variable/G SrwMagSASETaperStart=0
Variable/G SrwMagSASETaperVal=0

//SrwSASEInRad(InRadName,InPower,Waist,WaistPos)
String/G SrwRadSASEName="InR"
Variable/G SrwRadSASEPower=100.
Variable/G SrwRadSASEWaist=111.
Variable/G SrwRadSASEWaistPos=0.

//SrwPrecSASEMain(PrecName,npart,rmax0,ncar,nptr,nscr,nscz,delz,zstop,iorb)
Variable/G SrwPrecSASEnpart=4096.
Variable/G SrwPrecSASErmax0=8.
Variable/G SrwPrecSASEncar=129
Variable/G SrwPrecSASEnptr=40
Variable/G SrwPrecSASEnscr=0.
Variable/G SrwPrecSASEnscz=0.
Variable/G SrwPrecSASEdelz=0.5
Variable/G SrwPrecSASEzstop=-1
Variable/G SrwPrecSASEiorb=1

//SrwPrecSASETime(PrecName,itdp,nslice,zsep,ntail)
Variable/G SrwPrecSASEitdp=1
Variable/G SrwPrecSASEnslice=1
Variable/G SrwPrecSASEzsep=1
Variable/G SrwPrecSASEntail=0

//SrwCntrlSASE(CntrlName,iPower,iRadPhase,iRadSize,iBmHorSize,iBmVertSize,iBunchFact,iDisplay)
Variable/G SrwCntrlSASEiPower=1
Variable/G SrwCntrlSASEiRadPhase=1
Variable/G SrwCntrlSASEiRadSize=1
Variable/G SrwCntrlSASEiBmHorSize=1
Variable/G SrwCntrlSASEiBmVertSize=1
Variable/G SrwCntrlSASEiBunchFact=1
Variable/G SrwCntrlSASEiDisplay=1

end

//+++++++++++++++++++++++++++++++++++++++
//
//Create General Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagContCreate(MagName)
String MagName=SrwMagContName
prompt MagName,"Name of the Magnetic Field structure"
Silent 1						|	 ...

SrwMagContName = MagName
String TotMagName = MagName + SrwMagContType

KillWaves/Z $TotMagName
Make/T/O/N=0 $TotMagName
//SetScale d -1E+23, 1E+23, SrwMagType_Group, $TotMagName

end

//+++++++++++++++++++++++++++++++++++++++
//
//Create General Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagContFull(MagName,Pos1,Elem1,Pos2,Elem2,Pos3,Elem3,Pos4,Elem4)
String MagName=SrwMagContName
Variable Pos1=0
String Elem1=""
Variable Pos2=0
String Elem2=""
Variable Pos3=0
String Elem3=""
Variable Pos4=0
String Elem4=""
prompt MagName,"Name of the General Magnetic Field structure"
prompt Pos1,"Long. start pos. of the 1-st magnet element [m]"
prompt Elem1,"1-st magnet element",popup " ;"+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwMagOptType,";","")
prompt Pos2,"Long. start pos. of the 2-nd magnet element [m]"
prompt Elem2,"2-nd magnet element",popup " ;"+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwMagOptType,";","")
prompt Pos3,"Long. start pos. of the 3-rd magnet element [m]"
prompt Elem3,"3-rd magnet element",popup " ;"+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwMagOptType,";","")
prompt Pos4,"Long. start pos. of the 4-th magnet element [m]"
prompt Elem4,"4-th magnet element",popup " ;"+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwMagOptType,";","")
Silent 1						|	 ...

Variable Elem1Defined = 1, Elem2Defined = 1, Elem3Defined = 1, Elem4Defined = 1
if((cmpstr(Elem1," ") == 0) %| (cmpstr(Elem1,"") == 0))
	Elem1Defined = 0
endif
if((cmpstr(Elem2," ") == 0) %| (cmpstr(Elem2,"") == 0))
	Elem2Defined = 0
endif
if((cmpstr(Elem3," ") == 0) %| (cmpstr(Elem3,"") == 0))
	Elem3Defined = 0
endif
if((cmpstr(Elem4," ") == 0) %| (cmpstr(Elem4,"") == 0))
	Elem4Defined = 0
endif

if((Elem1Defined == 0) %& (Elem2Defined == 0) %& (Elem3Defined == 0) %& (Elem3Defined == 0))
	Abort "No magnet elements selected"
endif

Variable AmOfElems = (Elem1Defined + Elem2Defined + Elem3Defined + Elem4Defined)

SrwMagContName = MagName
String TotMagName = MagName + SrwMagContType

KillWaves/Z $TotMagName
Make/T/O/N=(AmOfElems*2) $TotMagName
//SetScale d -1E+23, 1E+23, SrwMagType_Group, $TotMagName

Variable ElCount = 0
if(Elem1Defined != 0)
	$TotMagName[ElCount] = num2str(Pos1) // start long. position
	$TotMagName[ElCount + 1] = Elem1 // elem. name
	ElCount += 2
endif
if(Elem2Defined != 0)
	$TotMagName[ElCount] = num2str(Pos2) // start long. position
	$TotMagName[ElCount + 1] = Elem2 // elem. name
	ElCount += 2
endif
if(Elem3Defined != 0)
	$TotMagName[ElCount] = num2str(Pos3) // start long. position
	$TotMagName[ElCount + 1] = Elem3 // elem. name
	ElCount += 2
endif
if(Elem4Defined != 0)
	$TotMagName[ElCount] = num2str(Pos4) // start long. position
	$TotMagName[ElCount + 1] = Elem4 // elem. name
	ElCount += 2
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Add element to General Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagContAddElem(MagName,Pos1,Elem1)
String MagName=SrwMagContName+SrwMagContName
Variable Pos1=0
String Elem1
prompt MagName,"General Magnetic Field structure",popup Wavelist("*"+SrwMagContName,";","")
prompt Pos1,"Longitudinal (start) position of the next magnet [m]"
prompt Elem1,"Next magnet",popup " ;"+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwMagOptType,";","")
Silent 1						|	 ...

if(cmpstr(Elem1," ") == 0)
	Elem1 = ""
endif
if(cmpstr(Elem1,"") == 0)
	Abort "No magnet selected"
endif

SrwMagContName = MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]

Variable i0 = DimSize($MagName,0)
//if(Indx > 0)
//	Variable PosIndx = (Indx - 1)*2
//	if(PosIndx < i0)
//		$MagName[PosIndx] = num2str(Pos1) // start long. position
//		$MagName[PosIndx + 1] = Elem1 // elem. name
//		Return
//	endif
//endif

Variable AmOfElems = 1
String AuxWave1 = "aux1", AuxWave2 = "aux2"
Make/T/O/N=(AmOfElems) $AuxWave1
Make/T/O/N=(AmOfElems) $AuxWave2

$AuxWave1[0] = Elem1; $AuxWave2[0] = num2str(Pos1)
//$AuxWave1[1] = Elem2; $AuxWave2[1] = num2str(Pos2)

Variable i = 0
do
	if(cmpstr($AuxWave1[i], "") != 0)
		redimension/N=(i0 + i*2 + 2)  $MagName
		$MagName[i*2 + i0] = $AuxWave2[i] // start long. position
		$MagName[i*2 + 1 + i0] = $AuxWave1[i] // elem. name
	endif
	i += 1
while(i < AmOfElems)

KillWaves/Z $AuxWave1, $AuxWave2
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Drift
//not used
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagDriftCreate(LensName,Length)
String LensName="Drift"
Variable Length=0.341
prompt LensName,"Name of the Drift structure"
prompt Length,"Drift length [m]"
Silent 1						|	 ...

SrwMagOptName = LensName
String TotLensName = LensName + SrwMagOptType

KillWaves/Z $TotLensName
Make/T/O/N = 2 $TotLensName
$TotLensName[0] = "Drift" // identificator
$TotLensName[1] = num2str(Length)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Quadrupole lens
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagQuadCreate(LensName,LensType,LensStrength,LensLength,HorPos,VertPos)
String LensName="Quad"
Variable LensType=1
Variable LensStrength=18.3
Variable LensLength=0.137
Variable HorPos=0
Variable VertPos=0
prompt LensName,"Name of the Quadrupole Lens structure"
prompt LensType,"Quadrupole Type",popup "Focusing;Defocusing"
prompt LensStrength,"Quadrupole Strength (absolute value) [T/m]"
prompt LensLength,"Effective Length [m]"
prompt HorPos,"Horizontal Position [mm]"
prompt VertPos,"Vertical Position [mm]"
Silent 1						|	 ...
PauseUpdate

// Validation of parameters
if(LensLength < 0)
	Abort "Effective length should be positive"
endif

SrwMagOptName = LensName
String TotLensName = LensName + SrwMagOptType

if(LensStrength < 0)
	Lenstrength = -LensStrength
endif
if(LensType == 2) // defocusing
	LensStrength = -LensStrength
endif

KillWaves/Z $TotLensName
Make/T/O/N = 10 $TotLensName
//SetScale d -1E+23, 1E+23, SrwMagType_Quad, $MagName

$TotLensName[0] = "Quadrupole" // identificator
$TotLensName[1] = num2str(LensStrength)
$TotLensName[2] = num2str(LensLength)
$TotLensName[3] = num2str(HorPos*0.001)
$TotLensName[4] = num2str(VertPos*0.001)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create simple 4-Magnet Chicane
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagChicaneCreate(ChicaneName, DipoleB, DipoleLen, DipoleNum, TotDriftLen)
string ChicaneName=srwUtiGetValS("SrwMagSASEChicName", "Chicane", "")
variable DipoleB=srwUtiGetValN("SrwMagSASEChicDipoleB", 0.015, "") //IBFIELD in GENESIS 
variable DipoleLen=srwUtiGetValN("SrwMagSASEChicDipoleLen", 0.25, "") //IMAGL in GENESIS
variable DipoleNum=srwUtiGetValN("SrwMagSASEChicDipoleNum", 4, "")
variable TotDriftLen=srwUtiGetValN("SrwMagSASEChicTotDriftLen", 0.13, "") //IDRIL in GENESIS
prompt ChicaneName,"Name for the Chicane structure"
prompt DipoleB,"Magnetic Field in Bending Magnets [T]"
prompt DipoleLen,"Length of each of Bending Magnets [m]"
prompt DipoleNum,"Number of Bending Magnets in the Chicane"
prompt TotDriftLen,"Total Drift Length of the Chicane (sum of 5 drifts) [m]"
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(DipoleLen < 0)
	Abort "Bending magnet length should be positive"
endif
if(DipoleLen == 0)
	DipoleNum = 0
endif

if(DipoleNum < 0)
	Abort "Number of Bending Magnets should be positive"
endif
if(TotDriftLen < 0)
	Abort "Total drift length between bending magnets should be positive"
endif

if((DipoleB == 0) %& (TotDriftLen > 0))
	DipoleB = 0.01
	TotDriftLen += 4*DipoleLen
	DipoleLen = 0
endif

srwUtiSetValS("SrwMagSASEChicName", ChicaneName, "")
srwUtiSetValN("SrwMagSASEChicDipoleB", DipoleB, "")
srwUtiSetValN("SrwMagSASEChicDipoleLen", DipoleLen, "")
srwUtiSetValN("SrwMagSASEChicDipoleNum", DipoleNum, "")
srwUtiSetValN("SrwMagSASEChicTotDriftLen", TotDriftLen, "")

String TotChicaneName = ChicaneName + SrwMagOptType

KillWaves/Z $TotChicaneName
Make/T/O/N = 10 $TotChicaneName

$TotChicaneName[0] = "Chicane" // identificator
$TotChicaneName[1] = num2str(DipoleB)
$TotChicaneName[2] = num2str(DipoleLen)
$TotChicaneName[3] = num2str(DipoleNum)
$TotChicaneName[4] = num2str(TotDriftLen)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Undulator Magnetic field for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEUndCreate(UndName,PlanOrEl,K,Period,SecLen,Nsec,IntBwSec,xkx,xky)
String UndName=SrwMagContName
Variable Period=SrwPeriod
Variable SecLen=SrwMagSASESecLen
Variable PlanOrEl=SrwMagSASEPlanOrEl
Variable K=SrwKz
Variable xkx=SrwMagSASExkx
Variable xky=SrwMagSASExky
Variable Nsec=SrwMagSASENsec
Variable IntBwSec=SrwMagSASEIntBwSec
prompt UndName,"Name of the Magnetic Field structure"
prompt Period,SrwPPeriod
prompt SecLen,"Section Length [m]"
prompt PlanOrEl,"Planar or Helical",popup "Planar;Helical"
prompt K,"K (= 0.0934 x Per. x |B|max)"
prompt xkx,"Horizontal Natural Focusing (XKX)"
prompt xky,"Vertical Natural Focusing (XKY)"
prompt Nsec,"Number of Sections"
prompt IntBwSec,"Interval between Sections [m]"
Silent 1						|	...
PauseUpdate

// Validation of parameters
if(Period <= 0.)
	Abort SrwPAlertUndPeriod
endif
if(K < 0.)
	Abort SrwPAlertUndKK
endif

SrwMagContName=UndName
SrwUndName=UndName

SrwPeriod=Period
SrwKz=K
SrwMagSASExkx=xkx
SrwMagSASExky=xky
SrwMagSASESecLen=SecLen
SrwMagSASEPlanOrEl=PlanOrEl
SrwMagSASENsec=Nsec
SrwMagSASEIntBwSec=IntBwSec

String TotGenMagName=UndName+SrwMagContType
String TotUndName=UndName+SrwUndType

Variable Kx = 0, PhSh = 1.5708
if(PlanOrEl == 2) // Elliptical
	K /= sqrt(2) //OC 140306
	Kx = K
	//if(FocType == 1) // default
	//	xkx = 0.5
	//endif
endif
SrwMagPerCreate2D(UndName,Period,K,Kx,SecLen,PhSh,1,0,0)

redimension/N=(40) $TotUndName
$TotUndName[30]=num2str(xkx)
$TotUndName[31]=num2str(xky)

$TotUndName[32]=num2str(1)
$TotUndName[33]=num2str(0) //No errors

$TotUndName[35] = "1" //No Taper num2str(TaperType)
$TotUndName[36] = "0" //num2str(TaperStart)
$TotUndName[37] = "0" //num2str(TaperVal)

SrwMagContCreate(UndName)

Variable i = 1
Variable CurPos = 0
do
	SrwMagContAddElem(TotGenMagName,CurPos,TotUndName)
	CurPos += (SecLen + IntBwSec)
	i += 1
while(i <= Nsec)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Writes start positions of adjascent undulator segments
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagSetUndStartPos(MagName)
String MagName

Variable AmOfSect = (srwUtiLastIndOfTxtFld($MagName, SrwUndType) + 1)*0.5
Variable StartPos = str2num($MagName[0])
Variable SectLens
Variable i = 1
do
	SectLens = str2num($($MagName[i*2 - 1])[1])
	StartPos += SectLens
	$MagName[i*2] = num2str(StartPos)
	i += 1
while(i < AmOfSect)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Tapering parameters for Magnetic field for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEUndTaperAdd(MagName,TaperType,TaperStart,TaperVal)
String MagName=SrwMagContName+SrwMagContType
Variable TaperType=SrwMagSASETaperType
Variable TaperStart=SrwMagSASETaperStart
Variable TaperVal=SrwMagSASETaperVal
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt TaperType,"Taper Type",popup "Linear;Quadratic"
prompt TaperStart,"Taper Start Position [m]"
prompt TaperVal,"Relative Field Change due to the Taper"
Silent 1						|	...
PauseUpdate

SrwMagContName = MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]
String TotUndName = $MagName[1]
SrwUndName = TotUndName[0,strlen(TotUndName)-strlen(SrwUndType)-1]

SrwMagSASETaperType=TaperType
SrwMagSASETaperStart=TaperStart
SrwMagSASETaperVal=TaperVal

$TotUndName[35] = num2str(TaperType)
$TotUndName[36] = num2str(TaperStart)
$TotUndName[37] = num2str(TaperVal)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Add Undulator Field Errors for SASE computation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEUndErrAdd(MagName,FldErrType,FldErrRMS)
String MagName=SrwMagContName+SrwMagContType
Variable FldErrType=SrwMagSASEFldErrType
Variable FldErrRMS=SrwMagSASEFldErrRMS
prompt FldErrType,"Type of Magnetic Field Errors",popup "Uniform Uncorrelated;Uniform Correlated;Gaussian Uncorrelated;Gaussian Correlated"
prompt FldErrRMS,"RMS Magnetic Field Error [T]"

SrwMagContName = MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]
String TotUndName = $MagName[1]
SrwUndName = TotUndName[0,strlen(TotUndName)-strlen(SrwUndType)-1]

SrwMagSASEFldErrType=FldErrType
SrwMagSASEFldErrRMS=FldErrRMS

$TotUndName[32] = num2str(FldErrType)
$TotUndName[33] = num2str(FldErrRMS)
end

//+++++++++++++++++++++++++++++++++++++++
//
//FODO parameters for Magnetic field for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASE_FODO_Add(MagName,Drift,QuadF,QuadD,LenF,LenD,F1St,Dx,Dz)
string MagName=SrwMagContName+SrwMagContType
Variable Drift=SrwMagSASEDrift
Variable QuadF=SrwMagSASEQuadF
Variable QuadD=SrwMagSASEQuadD
Variable LenF=SrwMagSASELenF
Variable LenD=SrwMagSASELenD
Variable Dx=SrwMagSASEDx
Variable Dz=SrwMagSASEDz
Variable F1St=SrwMagSASEF1St
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt Drift,"Interval between F- and D-sections [m]"
prompt QuadF,"Focusing Strength [T/m]"
prompt QuadD,"Defocusing Strength [T/m]"
prompt LenF,"Length of F-section [m]"
prompt LenD,"Length of D-section [m]"
prompt Dx,"Max. Horizontal Position Error [mm]"
prompt Dz,"Max. Vertical Position Error [mm]"
prompt F1St,"Offset w. resp. to undulator [m]"
Silent 1						|	...
PauseUpdate

//Variable ThereIsFODO = 0
//if((QuadF != 0) %| (QuadD != 0))
//	ThereIsFODO = 1
//endif
//if(ThereIsFODO == 0)
//	Return
//endif

SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]

SrwMagSASEDrift=Drift
SrwMagSASEQuadF=QuadF
SrwMagSASEQuadD=QuadD
SrwMagSASELenF=LenF
SrwMagSASELenD=LenD
SrwMagSASEDx=Dx
SrwMagSASEDz=Dz
SrwMagSASEF1St=F1St
//SrwUtiMagSetUndStartPos(MagName) //sets start positions of segments assuming no gaps

//Re-initializing entire magnetic field structure with global variables (necessary for correct container content)
SrwSASEUndCreate(SrwMagContName,SrwMagSASEPlanOrEl,SrwKz,SrwPeriod,SrwMagSASESecLen,SrwMagSASENsec,SrwMagSASEIntBwSec,SrwMagSASExkx,SrwMagSASExky)

String NameQuadF="", NameQuadD=""
if(QuadF != 0)
	NameQuadF = SrwMagContName + "F"
	SrwMagQuadCreate(NameQuadF,1,QuadF,LenF,Dx,Dz)
	NameQuadF += SrwMagOptType
endif
if(QuadD != 0)
	NameQuadD = SrwMagContName + "D"
	SrwMagQuadCreate(NameQuadD,2,QuadD,LenD,Dx,Dz)
	NameQuadD += SrwMagOptType
endif

//Variable GapLen = LenF + LenD + Drift - 2*F1St
Variable TotFocLen = LenF + LenD + Drift
Variable iFodoStart = srwUtiLastIndOfTxtFld($MagName, SrwUndType) + 1
Variable AmOfSectFODO = iFodoStart*0.5 - 1

Variable IndLengthInUndStruct = 1
Variable CurSegmLen = 0

Variable UndSegmStart, UndSegmEnd, NextUndSegmStart,CurGapLen,CurOffset,FodoSegmStart
Variable i = 1
String UndStructName
do
	UndSegmStart = str2num($MagName[(i - 1)*2])
	UndStructName = $MagName[i*2 - 1]
	CurSegmLen = str2num($UndStructName[IndLengthInUndStruct])
	UndSegmEnd = UndSegmStart + CurSegmLen
	
	//NextUndSegmStart = str2num($MagName[i*2])
	NextUndSegmStart = UndSegmEnd
	if(dimsize($MagName, 0) > i*2) //to program better: wave may be not long enough here
		NextUndSegmStart = str2num($MagName[i*2])
	endif

	CurGapLen = NextUndSegmStart - UndSegmEnd
	CurOffset = (CurGapLen - TotFocLen)*0.5
	//CurOffset = (CurGapLen - TotFocLen)*0.5 + F1St

	//FodoSegmStart = UndSegmEnd + CurOffset
	FodoSegmStart = UndSegmStart - F1St //?
	
	if(QuadF != 0)
		SrwMagContAddElem(MagName,FodoSegmStart,NameQuadF)
	endif
	FodoSegmStart += (LenF + Drift)
	if(QuadD != 0)
		SrwMagContAddElem(MagName,FodoSegmStart,NameQuadD)
	endif
	i += 1
while(i <= AmOfSectFODO)
end

//+++++++++++++++++++++++++++++++++++++++
//
//FODO parameters for Magnetic field for SASE
//Obsolette
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEFODOAdd(MagName,Drift,QuadF,QuadD,LenF,LenD,Dx,Dz)
string MagName=SrwMagContName+SrwMagContType
Variable Drift=SrwMagSASEDrift
Variable QuadF=SrwMagSASEQuadF
Variable QuadD=SrwMagSASEQuadD
Variable LenF=SrwMagSASELenF
Variable LenD=SrwMagSASELenD
Variable Dx=SrwMagSASEDx
Variable Dz=SrwMagSASEDz
//Variable F1St=SrwMagSASEF1St
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt Drift,"Interval between F- and D-sections [m]"
prompt QuadF,"Focusing Strength [T/m]"
prompt QuadD,"Defocusing Strength [T/m]"
prompt LenF,"Length of F-section [m]"
prompt LenD,"Length of D-section [m]"
prompt Dx,"Max. Horizontal Position Error [mm]"
prompt Dz,"Max. Vertical Position Error [mm]"
//prompt F1St,"Offset w. resp. to undulator [m]"
Silent 1						|	...
PauseUpdate

//Variable ThereIsFODO = 0
//if((QuadF != 0) %| (QuadD != 0))
//	ThereIsFODO = 1
//endif
//if(ThereIsFODO == 0)
//	Return
//endif

SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]

SrwMagSASEDrift=Drift
SrwMagSASEQuadF=QuadF
SrwMagSASEQuadD=QuadD
SrwMagSASELenF=LenF
SrwMagSASELenD=LenD
SrwMagSASEDx=Dx
SrwMagSASEDz=Dz
//SrwMagSASEF1St=F1St
//SrwUtiMagSetUndStartPos(MagName) //sets start positions of segments assuming no gaps

//Re-initializing entire magnetic field structure with global variables (necessary for correct container content)
SrwSASEUndCreate(SrwMagContName,SrwMagSASEPlanOrEl,SrwKz,SrwPeriod,SrwMagSASESecLen,SrwMagSASENsec,SrwMagSASEIntBwSec,SrwMagSASExkx,SrwMagSASExky)

String NameQuadF="", NameQuadD=""
if(QuadF != 0)
	NameQuadF = SrwMagContName + "F"
	SrwMagQuadCreate(NameQuadF,1,QuadF,LenF,Dx,Dz)
	NameQuadF += SrwMagOptType
endif
if(QuadD != 0)
	NameQuadD = SrwMagContName + "D"
	SrwMagQuadCreate(NameQuadD,2,QuadD,LenD,Dx,Dz)
	NameQuadD += SrwMagOptType
endif

//Variable GapLen = LenF + LenD + Drift - 2*F1St
Variable TotFocLen = LenF + LenD + Drift
Variable iFodoStart = srwUtiLastIndOfTxtFld($MagName, SrwUndType) + 1
Variable AmOfSectFODO = iFodoStart*0.5 - 1

Variable IndLengthInUndStruct = 1
Variable CurSegmLen = 0

Variable UndSegmStart, UndSegmEnd, NextUndSegmStart,CurGapLen,CurOffset,FodoSegmStart
Variable i = 1
String UndStructName
do
	UndSegmStart = str2num($MagName[(i - 1)*2])
	UndStructName = $MagName[i*2 - 1]
	CurSegmLen = str2num($UndStructName[IndLengthInUndStruct])
	UndSegmEnd = UndSegmStart + CurSegmLen
	
	NextUndSegmStart = UndSegmEnd
	if(dimsize($MagName, 0) > i*2) //to program better: wave may be not long enough here
		NextUndSegmStart = str2num($MagName[i*2])
	endif
	
	CurGapLen = NextUndSegmStart - UndSegmEnd
	CurOffset = (CurGapLen - TotFocLen)*0.5
	//CurOffset = (CurGapLen - TotFocLen)*0.5 + F1St

	FodoSegmStart = UndSegmEnd + CurOffset
	
	if(QuadF != 0)
		SrwMagContAddElem(MagName,FodoSegmStart,NameQuadF)
	endif
	FodoSegmStart += (LenF + Drift)
	if(QuadD != 0)
		SrwMagContAddElem(MagName,FodoSegmStart,NameQuadD)
	endif
	i += 1
while(i <= AmOfSectFODO)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Add 4-Magnet Chicane to the Magnet structure
//for HGHG (SASE) computation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEChicaneAdd(MagName, DipoleB, DipoleLen, TotDriftLen)
string MagName=SrwMagContName+SrwMagContType
variable DipoleB=srwUtiGetValN("SrwMagSASEChicDipoleB", 0.015, "") //IBFIELD in GENESIS 
variable DipoleLen=srwUtiGetValN("SrwMagSASEChicDipoleLen", 0.25, "") //IMAGL in GENESIS
variable TotDriftLen=srwUtiGetValN("SrwMagSASEChicTotDriftLen", 0.13, "") //IDRIL in GENESIS
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt DipoleB,"Magnetic Field in Bending Magnets [T]"
prompt DipoleLen,"Length of each of Bending Magnets [m]"
prompt TotDriftLen,"Total Drift Length of the Chicane (sum of 5 drifts) [m]"
Silent 1						|	...
PauseUpdate

SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]

srwUtiSetValN("SrwMagSASEChicDipoleB", DipoleB, "")
srwUtiSetValN("SrwMagSASEChicDipoleLen", DipoleLen, "")
srwUtiSetValN("SrwMagSASEChicTotDriftLen", TotDriftLen, "")

string NameChicane = SrwMagContName + "C"
SrwMagChicaneCreate(NameChicane, DipoleB, DipoleLen, 4, TotDriftLen)
NameChicane += SrwMagOptType

variable TotChicaneLength = TotDriftLen + 4*DipoleLen
SrwMagContAddElem(MagName, -TotChicaneLength, NameChicane) 
//chicane is located before the first undulator segment, which typically starts from 0 longitudinal position
end

//+++++++++++++++++++++++++++++++++++++++
//
//More e-beam parameters for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEElecExtra(ElecName,TypeTransv,TypeLong,SigLong,ShotNoise)
string ElecName=SrwElecName+SrwElecType
Variable TypeTransv=SrwElecTypeTransv
Variable TypeLong=SrwElecTypeLong
Variable SigLong=SrwElecSigLong
Variable ShotNoise=SrwElecShotNoise
prompt ElecName,"Electron Beam structure",popup Wavelist("*"+SrwElecType,";","")
prompt TypeTransv,"Type of Transverse Particle Distribution",popup "Uniform;Gaussian;Parabolic"
prompt TypeLong,"Type of Longitudinal Particle Distribution",popup "Infinite Uniform;Gaussian"
prompt SigLong,"RMS Bunch Length (for Gaussian Longitudinal Distribution) [mm]"
prompt ShotNoise,"Shot Noise (i.e. phase fluctuations scaling) Factor"
Silent 1						|	...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]

if(DimSize($ElecName, 0) < 40)
	redimension/N=(40) $ElecName
endif

SrwElecTypeTransv=TypeTransv
SrwElecTypeLong=TypeLong
SrwElecSigLong=SigLong
SrwElecShotNoise=ShotNoise

Variable SigLong_m = SigLong*0.001

$ElecName[30] = TypeTransv
$ElecName[31] = TypeLong
$ElecName[32] = ShotNoise
$ElecName[33] = SigLong_m*SigLong_m
end

//+++++++++++++++++++++++++++++++++++++++
//
//Input Radiation parameters for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEInRadGsn(InRadName,Waist,WaistPos,InPower)
string InRadName=SrwRadSASEName
Variable InPower=SrwRadSASEPower
Variable Waist=SrwRadSASEWaist
Variable WaistPos=SrwRadSASEWaistPos
prompt InRadName,"Name of Input Radiation structure"
prompt InPower,"Input Radiation Power [W]"
prompt Waist,"RMS Waist of the Radiation Field [µm]"
prompt WaistPos,"Longitudinal Position of the Waist [m] (0 corresponds to exit of undulator)"
Silent 1						|	...
PauseUpdate

SrwRadSASEName = InRadName
SrwRadSASEPower=InPower
SrwRadSASEWaist=Waist
SrwRadSASEWaistPos=WaistPos

String TotInRadName = InRadName + SrwRadSASEType

Make/T/O/N=10 $TotInRadName

// first fields are reserved for rad. structure name
$TotInRadName[2] = num2str(InPower)
$TotInRadName[3] = num2str(Waist*0.000001)
$TotInRadName[4] = num2str(WaistPos)
//$TotInRadName[5] = num2str(PhEnSim*1000.)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Main Precision parameters for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEPrecGen(PrecName,delz,npart,UseElecDistr,ProduceElecDistr,photEn_xlamds,rmax0,ncar) //,nptr,nscr,nscz) //,iorb)
string PrecName=SrwPrecSASEName
variable delz=SrwPrecSASEdelz
variable npart=SrwPrecSASEnpart
variable UseElecDistr=srwUtiGetValN("UseElecDistr", 1, "SrwSASEPrecGen")
variable ProduceElecDistr=srwUtiGetValN("ProduceElecDistr", 1, "SrwSASEPrecGen")
variable photEn_xlamds=srwUtiGetValN("photEn_xlamds", 1, "SrwSASEPrecGen")
variable rmax0=SrwPrecSASErmax0
variable ncar=SrwPrecSASEncar
//variable nptr=SrwPrecSASEnptr
//variable nscr=SrwPrecSASEnscr
//variable nscz=SrwPrecSASEnscz
//variable zstop=SrwPrecSASEzstop
//variable iorb=SrwPrecSASEiorb
prompt PrecName,"Name of SASE Precision structure"
prompt delz,"Long. Step per Und. Period  (DELZ)"
prompt npart,"Num. of Macro-Particles (NPART)"
prompt rmax0,"Transverse Mesh Extent (RMAX0)"
prompt ncar,"Num. of Transverse Mesh Pts. (NCAR)"
//prompt nptr,"Num. of Rad. Pts. for Sp. Charge (NPTR)"
//prompt nscr,"Num. of Az. Modes for Sp. Charge (NSCR)"
//prompt nscz,"Num. of Long. Modes for Sp. Charge (NSCZ)"
//prompt lbc,"Type of boundary conditions",popup "Dirichlet;Neumann"
//prompt zstop,"Longit. Integration End Point [m]"
//prompt iorb,"Use Orbit Correction (IORB)",popup "No;Yes"
prompt UseElecDistr,"Use Pre-Calculated Particle Distribution?",popup "No;Yes"
prompt ProduceElecDistr,"Produce Particle Distribution in the next run?",popup "No;Yes"
prompt photEn_xlamds,"Reson. Photon Energy [keV] (~1/XLAMDS)"
Silent 1						|	...
PauseUpdate

String TotPrecName = PrecName + SrwPrecSASEType
SrwPrecSASEName = PrecName

SrwPrecSASEnpart=npart
SrwPrecSASErmax0=rmax0
SrwPrecSASEncar=ncar
//SrwPrecSASEnptr=nptr
//SrwPrecSASEnscr=nscr
//SrwPrecSASEnscz=nscz
SrwPrecSASEdelz=delz
SrwPrecSASEzstop=-1
//SrwPrecSASEiorb=iorb
srwUtiSetValN("UseElecDistr", UseElecDistr, "SrwSASEPrecGen")
srwUtiSetValN("ProduceElecDistr", ProduceElecDistr, "SrwSASEPrecGen")
srwUtiSetValN("photEn_xlamds", photEn_xlamds, "SrwSASEPrecGen")

Make/D/O/N=30 $TotPrecName
$TotPrecName[0] = npart
$TotPrecName[1] = rmax0
$TotPrecName[2] = ncar
//$TotPrecName[3] = nptr
//$TotPrecName[4] = nscr
//$TotPrecName[5] = nscz
$TotPrecName[6] = 1 //lbc
$TotPrecName[7] = delz
$TotPrecName[8] = -1 //zstop
//$TotPrecName[9] = iorb
$TotPrecName[14] = UseElecDistr - 1
$TotPrecName[15] = ProduceElecDistr - 1

$TotPrecName[16] = photEn_xlamds*1000

if($TotPrecName[10] <= 0) 
	$TotPrecName[10] = 1 
endif
if($TotPrecName[11] <= 0) 
	$TotPrecName[11] = 1 
endif
if($TotPrecName[12] <= 0) 
	$TotPrecName[12] = delz 
endif
//$TotPrecName[10] = itdp
//$TotPrecName[11] = nslice
//$TotPrecName[12] = zsep
//$TotPrecName[13] = ntail
//SrwSASEPrecTime(TotPrecName,1,1,delz,0) // default values
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extra Precision parameters for SASE (space charge, etc.)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEPrecExtra(PrecName,nptr,nscr,nscz,iorb)
string PrecName=SrwPrecSASEName+SrwPrecSASEType
variable nptr=SrwPrecSASEnptr
variable nscr=SrwPrecSASEnscr
variable nscz=SrwPrecSASEnscz
variable iorb=SrwPrecSASEiorb
prompt PrecName,"SASE Precision Parameters structure",popup Wavelist("*"+SrwPrecSASEType,";","")
prompt nptr,"Number of Transverse Points for Space Charge (NPTR)"
prompt nscr,"Number of Azimuthal Modes for Space Charge (NSCR)"
prompt nscz,"Number of Longitudinal Modes for Space Charge (NSCZ)"
prompt iorb,"Use Orbit Correction (IORB)",popup "No;Yes"
Silent 1						|	...
PauseUpdate

SrwPrecSASEName=PrecName[0,strlen(PrecName)-strlen(SrwPrecSASEType)-1]
string TotPrecName = PrecName

SrwPrecSASEnptr=nptr
SrwPrecSASEnscr=nscr
SrwPrecSASEnscz=nscz
SrwPrecSASEiorb=iorb

if(exists(TotPrecName) == 0)
	Make/D/O/N=30 $TotPrecName
endif

$TotPrecName[3] = nptr
$TotPrecName[4] = nscr
$TotPrecName[5] = nscz
$TotPrecName[9] = iorb
end

//+++++++++++++++++++++++++++++++++++++++
//
//Main Precision parameters for SASE
//Obsolete version, use "SrwSASEPrecGen" instead
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEPrecMain(PrecName,delz,npart,rmax0,ncar,nptr,nscr,nscz,iorb)
string PrecName=SrwPrecSASEName
Variable npart=SrwPrecSASEnpart
Variable rmax0=SrwPrecSASErmax0
Variable ncar=SrwPrecSASEncar
Variable nptr=SrwPrecSASEnptr
Variable nscr=SrwPrecSASEnscr
Variable nscz=SrwPrecSASEnscz
Variable delz=SrwPrecSASEdelz
//Variable zstop=SrwPrecSASEzstop
Variable iorb=SrwPrecSASEiorb
prompt PrecName,"Name of SASE Precision structure"
prompt npart,"Number of Macro-Particles"
prompt rmax0,"Transverse Mesh Extent Factor"
prompt ncar,"Number of Transverse Mesh Pts."
prompt nptr,"Num. of Radial Pts. for Space Charge"
prompt nscr,"Num. of Azim. Modes for Space Charge"
prompt nscz,"Num. of Longit. Modes for Space Charge"
//prompt lbc,"Type of boundary conditions",popup "Dirichlet;Neumann"
prompt delz,"Longit. Step divided by Undulator Period"
//prompt zstop,"Longit. Integration End Point [m]"
prompt iorb,"Use Orbit Correction",popup "No;Yes"
Silent 1						|	...
PauseUpdate

String TotPrecName = PrecName + SrwPrecSASEType
SrwPrecSASEName = PrecName

SrwPrecSASEnpart=npart
SrwPrecSASErmax0=rmax0
SrwPrecSASEncar=ncar
SrwPrecSASEnptr=nptr
SrwPrecSASEnscr=nscr
SrwPrecSASEnscz=nscz
SrwPrecSASEdelz=delz
SrwPrecSASEzstop=-1
SrwPrecSASEiorb=iorb

Make/D/O/N=30 $TotPrecName
$TotPrecName[0] = npart
$TotPrecName[1] = rmax0
$TotPrecName[2] = ncar
$TotPrecName[3] = nptr
$TotPrecName[4] = nscr
$TotPrecName[5] = nscz
$TotPrecName[6] = 1 //lbc
$TotPrecName[7] = delz
$TotPrecName[8] = -1 //zstop
$TotPrecName[9] = iorb

SrwSASEPrecTime(TotPrecName,1,1,delz,0,0,0) // default values
end

//+++++++++++++++++++++++++++++++++++++++
//
//Time-dependency Precision parameters for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASEPrecTime(PrecName,itdp,nslice,zsep_d_delz,ntail,alignradf,offsetradf)
string PrecName=SrwPrecSASEName+SrwPrecSASEType
variable itdp=SrwPrecSASEitdp
variable nslice=SrwPrecSASEnslice
variable zsep_d_delz=trunc(SrwPrecSASEzsep/SrwPrecSASEdelz + 0.00001)
variable ntail=SrwPrecSASEntail
variable alignradf=srwUtiGetValN("SrwPrecSASEalignradf", 1, "")
variable offsetradf=srwUtiGetValN("SrwPrecSASEoffsetradf", 0, "")
prompt PrecName,"SASE Precision Parameters structure",popup Wavelist("*"+SrwPrecSASEType,";","")
prompt itdp,"Run in Time Dependent Mode (ITDP)",popup "No;Yes"
prompt nslice,"Number of Bunch Slices (NSLICE)"
prompt zsep_d_delz,"Slice Separ. per Long. Step (ZSEP/DELZ)", popup "1;2;3;4;5;6;7;8;9;10;15;20;25;30"
prompt ntail,"First Slice Position rel. to Bunch Cen. (NTAIL)" // (*ZSEP*Wavelength, GENESIS: NTAIL)"
prompt alignradf,"Align Seed Field manually? (ALIGNRADF)", popup "No;Yes"
prompt offsetradf,"Seed Field Offset in Slices (OFFSETRADF)"
Silent 1						|	...
PauseUpdate

SrwPrecSASEName=PrecName[0,strlen(PrecName)-strlen(SrwPrecSASEType)-1]
String TotPrecName = PrecName

if(zsep_d_delz == 11)
	zsep_d_delz = 15
else
	if(zsep_d_delz == 12)
		zsep_d_delz = 20
	else
		if(zsep_d_delz == 13)
			zsep_d_delz = 25
		else
			if(zsep_d_delz == 14)
				zsep_d_delz = 30
			endif
		endif
	endif
endif

variable zsep = zsep_d_delz*SrwPrecSASEdelz
SrwPrecSASEitdp=itdp
SrwPrecSASEnslice=nslice
SrwPrecSASEzsep=zsep
SrwPrecSASEntail=ntail
srwUtiSetValN("SrwPrecSASEalignradf", alignradf, "")
srwUtiSetValN("SrwPrecSASEoffsetradf", offsetradf, "")

if(exists(TotPrecName) == 0)
	Make/D/O/N=30 $TotPrecName
endif

$TotPrecName[10] = itdp
$TotPrecName[11] = nslice
$TotPrecName[12] = zsep
$TotPrecName[13] = ntail
$TotPrecName[17] = alignradf
$TotPrecName[18] = offsetradf
end

//+++++++++++++++++++++++++++++++++++++++
//
//Output control parameters for SASE
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASECntrl(CntrlName,iPower,iRadPhase,iRadSize,iBmHorSize,iBmVertSize,iBunchFact,iDisplay)
string CntrlName=SrwCntrlSASEName
variable iPower=SrwCntrlSASEiPower
variable iRadPhase=SrwCntrlSASEiRadPhase
variable iRadSize=SrwCntrlSASEiRadSize
variable iBmHorSize=SrwCntrlSASEiBmHorSize
variable iBmVertSize=SrwCntrlSASEiBmVertSize
variable iBunchFact=SrwCntrlSASEiBunchFact
variable iDisplay=SrwCntrlSASEiDisplay
prompt CntrlName,"Name of SASE Control structure"
prompt iPower,"Radiation Power",popup "No;Yes"
prompt iRadPhase,"On-axis Radiation Phase",popup "No;Yes"
prompt iRadSize,"RMS Size of Power Density Distribution",popup "No;Yes"
prompt iBmHorSize,"Hor. RMS E-Beam Size",popup "No;Yes"
prompt iBmVertSize,"Vertical RMS E-Beam Size",popup "No;Yes"
prompt iBunchFact,"Bunching Factor",popup "No;Yes"
prompt iDisplay,"New Display",popup "No;Yes"
Silent 1						|	...
PauseUpdate

String TotCntrlName = CntrlName + SrwCntrlSASEType
SrwCntrlSASEName = CntrlName

SrwCntrlSASEiPower=iPower
SrwCntrlSASEiRadPhase=iRadPhase
SrwCntrlSASEiRadSize=iRadSize
SrwCntrlSASEiBmHorSize=iBmHorSize
SrwCntrlSASEiBmVertSize=iBmVertSize
SrwCntrlSASEiBunchFact=iBunchFact
SrwCntrlSASEiDisplay=iDisplay

Make/T/O/N=(40,2) $TotCntrlName
variable i=0
do
	$TotCntrlName[i][0] = ""
	$TotCntrlName[i][1] = ""
	i += 1
while(i <= 40)

if(iPower == 2)
	$TotCntrlName[0][0] = CntrlName + "001_aux"
	$TotCntrlName[0][1] = "Radiation Power in Slice at Harm. #1"
	SrwUtiCreateNumWave1D($TotCntrlName[0][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[1][0] = CntrlName + "002_aux"
	$TotCntrlName[1][1] = "Radiation Power in Slice at Harm. #2"
	SrwUtiCreateNumWave1D($TotCntrlName[1][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[2][0] = CntrlName + "003_aux"
	$TotCntrlName[2][1] = "Radiation Power in Slice at Harm. #3"
	SrwUtiCreateNumWave1D($TotCntrlName[2][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[3][0] = CntrlName + "004_aux"
	$TotCntrlName[3][1] = "Radiation Power in Slice at Harm. #4"
	SrwUtiCreateNumWave1D($TotCntrlName[3][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[4][0] = CntrlName + "005_aux"
	$TotCntrlName[4][1] = "Radiation Power in Slice at Harm. #5"
	SrwUtiCreateNumWave1D($TotCntrlName[4][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[5][0] = CntrlName + "006_aux"
	$TotCntrlName[5][1] = "Radiation Power in Slice at Harm. #6"
	SrwUtiCreateNumWave1D($TotCntrlName[5][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[6][0] = CntrlName + "007_aux"
	$TotCntrlName[6][1] = "Radiation Power in Slice at Harm. #7"
	SrwUtiCreateNumWave1D($TotCntrlName[6][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	
	$TotCntrlName[7][0] = CntrlName + "101_aux"
	$TotCntrlName[7][1] = "Power vs Time at Harm. #1"
	make/O/N=(10,10) $($TotCntrlName[7][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[7][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[7][0])
	string strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[7][0] + "'"; execute strExe
	$TotCntrlName[8][0] = CntrlName + "102_aux"
	$TotCntrlName[8][1] = "Power vs Time at Harm. #2"
	make/O/N=(10,10) $($TotCntrlName[8][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[8][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[8][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[8][0] + "'"; execute strExe
	$TotCntrlName[9][0] = CntrlName + "103_aux"
	$TotCntrlName[9][1] = "Power vs Time at Harm. #3"
	make/O/N=(10,10) $($TotCntrlName[9][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[9][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[9][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[9][0] + "'"; execute strExe
	$TotCntrlName[10][0] = CntrlName + "104_aux"
	$TotCntrlName[10][1] = "Power vs Time at Harm. #4"
	make/O/N=(10,10) $($TotCntrlName[10][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[10][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[10][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[10][0] + "'"; execute strExe
	$TotCntrlName[11][0] = CntrlName + "105_aux"
	$TotCntrlName[11][1] = "Power vs Time at Harm. #5"
	make/O/N=(10,10) $($TotCntrlName[11][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[11][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[11][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[11][0] + "'"; execute strExe
	$TotCntrlName[12][0] = CntrlName + "106_aux"
	$TotCntrlName[12][1] = "Power vs Time at Harm. #6"
	make/O/N=(10,10) $($TotCntrlName[12][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[12][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[12][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[12][0] + "'"; execute strExe
	$TotCntrlName[13][0] = CntrlName + "107_aux"
	$TotCntrlName[13][1] = "Power vs Time at Harm. #7"
	make/O/N=(10,10) $($TotCntrlName[13][0]) //to be redimensioned in C
	SetScale/P x, 0, 1, "s", $($TotCntrlName[13][0])
	SetScale d, 0, 0, "W", $($TotCntrlName[13][0])
	strExe = "SetDimLabel 0, -1, Time, '" + $TotCntrlName[13][0] + "'"; execute strExe
	
	$TotCntrlName[14][0] = CntrlName + "201_aux"
	$TotCntrlName[14][1] = "Radiation Energy at Harm. #1"
	SrwUtiCreateNumWave1D($TotCntrlName[14][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[15][0] = CntrlName + "202_aux"
	$TotCntrlName[15][1] = "Radiation Energy at Harm. #2"
	SrwUtiCreateNumWave1D($TotCntrlName[15][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[16][0] = CntrlName + "203_aux"
	$TotCntrlName[16][1] = "Radiation Energy at Harm. #3"
	SrwUtiCreateNumWave1D($TotCntrlName[16][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[17][0] = CntrlName + "204_aux"
	$TotCntrlName[17][1] = "Radiation Energy at Harm. #4"
	SrwUtiCreateNumWave1D($TotCntrlName[17][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[18][0] = CntrlName + "205_aux"
	$TotCntrlName[18][1] = "Radiation Energy at Harm. #5"
	SrwUtiCreateNumWave1D($TotCntrlName[18][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[19][0] = CntrlName + "206_aux"
	$TotCntrlName[19][1] = "Radiation Energy at Harm. #6"
	SrwUtiCreateNumWave1D($TotCntrlName[19][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	$TotCntrlName[20][0] = CntrlName + "207_aux"
	$TotCntrlName[20][1] = "Radiation Energy at Harm. #7"
	SrwUtiCreateNumWave1D($TotCntrlName[20][0], "", 10, 0, 1, "m", "Longitudinal Position", "J")
	
	$TotCntrlName[21][0] = CntrlName + "301_aux"
	$TotCntrlName[21][1] = "Peak Radiation Power at Harm. #1"
	SrwUtiCreateNumWave1D($TotCntrlName[21][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[22][0] = CntrlName + "302_aux"
	$TotCntrlName[22][1] = "Peak Radiation Power at Harm. #2"
	SrwUtiCreateNumWave1D($TotCntrlName[22][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[23][0] = CntrlName + "303_aux"
	$TotCntrlName[23][1] = "Peak Radiation Power at Harm. #3"
	SrwUtiCreateNumWave1D($TotCntrlName[23][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[24][0] = CntrlName + "304_aux"
	$TotCntrlName[24][1] = "Peak Radiation Power at Harm. #4"
	SrwUtiCreateNumWave1D($TotCntrlName[24][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[25][0] = CntrlName + "305_aux"
	$TotCntrlName[25][1] = "Peak Radiation Power at Harm. #5"
	SrwUtiCreateNumWave1D($TotCntrlName[25][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[26][0] = CntrlName + "306_aux"
	$TotCntrlName[26][1] = "Peak Radiation Power at Harm. #6"
	SrwUtiCreateNumWave1D($TotCntrlName[26][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
	$TotCntrlName[27][0] = CntrlName + "307_aux"
	$TotCntrlName[27][1] = "Peak Radiation Power at Harm. #7"
	SrwUtiCreateNumWave1D($TotCntrlName[27][0], "", 10, 0, 1, "m", "Longitudinal Position", "W")
endif

if(iRadPhase == 2)
	$TotCntrlName[28][0]= CntrlName + "008_aux"
	$TotCntrlName[28][1] = "On-Axis Radiation Phase"
	SrwUtiCreateNumWave1D($TotCntrlName[28][0], "", 10, 0, 1, "m", "Longitudinal Position", "r")
endif
if(iRadSize == 2)
	$TotCntrlName[29][0] = CntrlName + "009_aux"
	$TotCntrlName[29][1] = "RMS Size of Radiation Power Density Distribution"
	SrwUtiCreateNumWave1D($TotCntrlName[29][0], "", 10, 0, 1, "m", "Longitudinal Position", "m")
endif
if(iBmHorSize == 2)
	$TotCntrlName[30][0] = CntrlName + "010_aux"
	$TotCntrlName[30][1] = "Horizontal RMS E-Beam Size"
	SrwUtiCreateNumWave1D($TotCntrlName[30][0], "", 10, 0, 1, "m", "Longitudinal Position", "m")
endif
if(iBmVertSize == 2)
	$TotCntrlName[31][0] = CntrlName + "011_aux"
	$TotCntrlName[31][1] = "Vertical RMS E-Beam Size"
	SrwUtiCreateNumWave1D($TotCntrlName[31][0], "", 10, 0, 1, "m", "Longitudinal Position", "m")
endif

if(iBunchFact == 2)
	$TotCntrlName[32][0] = CntrlName + "401_aux"
	$TotCntrlName[32][1] = "Bunching Factor at Harm. #1"
	SrwUtiCreateNumWave1D($TotCntrlName[32][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[33][0] = CntrlName + "402_aux"
	$TotCntrlName[33][1] = "Bunching Factor at Harm. #2"
	SrwUtiCreateNumWave1D($TotCntrlName[33][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[34][0] = CntrlName + "403_aux"
	$TotCntrlName[34][1] = "Bunching Factor at Harm. #3"
	SrwUtiCreateNumWave1D($TotCntrlName[34][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[35][0] = CntrlName + "404_aux"
	$TotCntrlName[35][1] = "Bunching Factor at Harm. #4"
	SrwUtiCreateNumWave1D($TotCntrlName[35][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[36][0] = CntrlName + "405_aux"
	$TotCntrlName[36][1] = "Bunching Factor at Harm. #5"
	SrwUtiCreateNumWave1D($TotCntrlName[36][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[37][0] = CntrlName + "406_aux"
	$TotCntrlName[37][1] = "Bunching Factor at Harm. #6"
	SrwUtiCreateNumWave1D($TotCntrlName[37][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
	$TotCntrlName[38][0] = CntrlName + "407_aux"
	$TotCntrlName[38][1] = "Bunching Factor at Harm. #7"
	SrwUtiCreateNumWave1D($TotCntrlName[38][0], "", 10, 0, 1, "m", "Longitudinal Position", "")
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Zeros SASE control output in graphs
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASECntrlZero(CntrlName)
String CntrlName
Silent 1						|	...
PauseUpdate

String DataLabel, ArgLabel, WaveNameStr
Variable MaxAmOfNumWaves = DimSize($CntrlName, 0)
Variable i = 0
do
	if((cmpstr($CntrlName[i][0],"") != 0))
		WaveNameStr = $CntrlName[i][0]
		$WaveNameStr = 0
	endif
	i += 1
while(i < MaxAmOfNumWaves)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Displays SASE control output in graphs
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSASECntrlDisplay(CntrlName)
string CntrlName
Silent 1						|	...
PauseUpdate

//make validation here

string strAuxBunchFactLabel = "Bunching Factor", strAuxTestDataLabel = ""
string strAuxPowerVsTimeLabel = "Power vs Time"
string strAuxRadPowerInSliceLabel = "Radiation Power in Slice"
string strAuxRadEnLabel = "Radiation Energy"
string strAuxPeakRadPowerLabel = "Peak Radiation Power"

string DataLabel, ArgLabel, WaveNameStr, AllWaveNamesStr, actDataLabel = ""
string DataLabelTestCore
string StrToExe = "", strHarmNo = "", strHarmID = ""

variable MaxAmOfNumWaves = DimSize($CntrlName, 0)
variable i = 0, k = 0, lenStrLabel = 0

//variable lenStrAuxBunchFactLabel = strlen(strAuxBunchFactLabel)
//variable lenStrAuxRadPowerInSliceLabel = strlen(strAuxRadPowerInSliceLabel)
//variable lenStrAuxRadEnLabel = strlen(strAuxRadEnLabel)
//variable lenStrAuxPeakRadPowerLabel = strlen(strAuxPeakRadPowerLabel)

variable isBunchingFact = 0, bunchingFactGraphOpened = 0, newGraphIsNecessary = 0
variable isRadPowerInSlice = 0, radPowerInSliceGraphOpened = 0
variable isRadEn = 0, radEnGraphOpened = 0
variable isPeakRadPower = 0, peakRadPowerGraphOpened = 0
variable isPowerVsTime = 0, powerVsTimeGraphOpened = 0
variable argOffset, argDelta, argNp, argValForTag

do
	if((cmpstr($CntrlName[i][0],"") != 0) %& (SrwCntrlSASEiDisplay == 2))
		WaveNameStr = $CntrlName[i][0]
		DataLabel = $CntrlName[i][1]
		ArgLabel = GetDimLabel($WaveNameStr, 0, -1)
		
		DataLabelTestCore = DataLabel[0, strlen(DataLabel) - strlen(" at Harm. #1") - 1]
		
		isBunchingFact = 0
		isRadPowerInSlice = 0
		isRadEn = 0
		isPeakRadPower = 0
		isPowerVsTime = 0

		if(cmpstr(DataLabelTestCore, strAuxBunchFactLabel) == 0)
			isBunchingFact = 1
		endif
		if(cmpstr(DataLabelTestCore, strAuxRadPowerInSliceLabel) == 0)
			isRadPowerInSlice = 1
		endif
		if(cmpstr(DataLabelTestCore, strAuxRadEnLabel) == 0)
			isRadEn = 1
		endif
		if(cmpstr(DataLabelTestCore, strAuxPeakRadPowerLabel) == 0)
			isPeakRadPower = 1
		endif
		if(cmpstr(DataLabelTestCore, strAuxPowerVsTimeLabel) == 0)
			isPowerVsTime = 1
		endif
		
		newGraphIsNecessary = 1		
		if((isBunchingFact != 0) %& (bunchingFactGraphOpened != 0))
			newGraphIsNecessary = 0
		endif
		if((isRadPowerInSlice != 0) %& (radPowerInSliceGraphOpened != 0))
			newGraphIsNecessary = 0
		endif
		if((isRadEn != 0) %& (radEnGraphOpened != 0))
			newGraphIsNecessary = 0
		endif
		if((isPeakRadPower != 0) %& (peakRadPowerGraphOpened != 0))
			newGraphIsNecessary = 0
		endif
		if((isPowerVsTime != 0) %& (powerVsTimeGraphOpened != 0))
			newGraphIsNecessary = 0
		endif

		//lenStrLabel = strlen(DataLabel)
		//if(lenStrLabel > lenStrAuxBunchFactLabel) 
		//	lenStrLabel = lenStrAuxBunchFactLabel
		//endif
		//strAuxTestDataLabel = srwUtiSubStr(DataLabel, 0, lenStrLabel)
		
		if(newGraphIsNecessary != 0)
		
			if(isPowerVsTime != 0)
				Display $WaveNameStr[][dimsize($WaveNameStr, 1) - 1]
			else
				Display $WaveNameStr
			endif
			Label bottom ArgLabel
			
			actDataLabel = DataLabel

			if(isBunchingFact != 0)
				//actDataLabel = strAuxBunchFactLabel + "s"
				actDataLabel = DataLabelTestCore
				bunchingFactGraphOpened = 1
			endif
			if(isRadPowerInSlice != 0)
				actDataLabel = DataLabelTestCore
				radPowerInSliceGraphOpened = 1
				ModifyGraph log(left)=1
			endif
			if(isRadEn != 0)
				actDataLabel = DataLabelTestCore
				radEnGraphOpened = 1
				ModifyGraph log(left)=1
			endif
			if(isPeakRadPower != 0)
				actDataLabel = DataLabelTestCore
				peakRadPowerGraphOpened = 1
				ModifyGraph log(left)=1
			endif
			if(isPowerVsTime != 0)
				actDataLabel = DataLabelTestCore
				powerVsTimeGraphOpened = 1
			endif
		
			Textbox/A=MT/N=text0/F=0 actDataLabel
		
			ModifyGraph rgb=(0,0,52224)
			ModifyGraph mode=4,marker=19,lsize=1
			
			if(strlen(AllWaveNamesStr) > 0)
				AllWaveNamesStr += ","
			endif
			AllWaveNamesStr += WinName(0,1)
		else
		
			if(isPowerVsTime != 0)
				AppendToGraph $WaveNameStr[][dimsize($WaveNameStr, 1) - 1]
			else
				AppendToGraph $WaveNameStr
			endif
			
			ModifyGraph mode=4,marker=19,lsize=1
			ModifyGraph rgb=(0,0,52224)
		endif
		
		if(isBunchingFact %| isRadPowerInSlice %| isRadEn %| isPeakRadPower %| isPowerVsTime)
			strHarmNo = DataLabel[strlen(DataLabel) - 1]
			//strHarmID = srwUtiSubStr(DataLabel, strlen(DataLabel) - 8, 8)
			strHarmID = "Harm. #" + strHarmNo
			
			argValForTag = 100000
			StrToExe = "Tag/C/N=text" + strHarmNo + "/X=-10/Y=0/Z=0/L=2/D=0.3 " + WaveNameStr + ", " + num2str(argValForTag) + ", \"" + strHarmID + "\""
			Execute StrToExe
		endif
	endif
	i += 1
while(i < MaxAmOfNumWaves)

StrToExe = ""
if(strlen(AllWaveNamesStr) > 0)
	StrToExe = "TileWindows " + AllWaveNamesStr
	Execute StrToExe
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Start SASE computation (main version, to be used in prference)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrSASECreate(WfrName,ElecName,MagName,InRadName,RadSmpName,PrecName,CntrlName,ObsNxNzForProp,ObsNxNzSamplFact)
string WfrName="Wfr"
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagContName+SrwMagContType
string InRadName=SrwRadSASEName+SrwRadSASEType
string RadSmpName=SrwSmpName+SrwSmpType
string PrecName=SrwPrecSASEName+SrwPrecSASEType
string CntrlName=SrwCntrlSASEName+SrwCntrlSASEType
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt WfrName,"Name of the Wavefront structure"
prompt ElecName,"Electron Beam structure",popup Wavelist("*"+SrwElecType,";","")
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt InRadName,"Input Radiation structure",popup Wavelist("*"+SrwRadSASEType,";","") + Wavelist("*"+SrwRadType,";","")
prompt RadSmpName,"Radiation Sampling structure (optional)",popup Wavelist("*"+SrwSmpType,";","")
prompt PrecName,"SASE Precision structure",popup Wavelist("*"+SrwPrecSASEType,";","")
prompt CntrlName,"SASE Control structure",popup Wavelist("*"+SrwCntrlSASEType,";","")
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	...
//PauseUpdate

if((strlen(ElecName) == 0) %| (cmpstr(ElecName, "_none_") == 0))
	abort "Electron Beam structure is not defined."
endif
if((strlen(MagName) == 0) %| (cmpstr(MagName, "_none_") == 0))
	abort "Magnetic Field structure is not defined."
endif
if((strlen(InRadName) == 0) %| (cmpstr(InRadName, "_none_") == 0))
	abort "Input Radiation structure is not defined."
endif
if((strlen(PrecName) == 0) %| (cmpstr(PrecName, "_none_") == 0))
	abort "Precision Parameters structure is not defined."
endif
if((strlen(CntrlName) == 0) %| (cmpstr(CntrlName, "_none_") == 0))
	abort "Control structure is not defined."
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]
SrwRadSASEName=InRadName[0,strlen(InRadName)-strlen(SrwRadSASEType)-1]

//variable RadSamplindIsDefined = 0
if((strlen(RadSmpName) > 0) %& (cmpstr(RadSmpName, "_none_") != 0))
//SASE computation can be performed even without Observation structure (just to the exit of undulator, with the sampling specified in PrecName)
	SrwSmpName=RadSmpName[0,strlen(RadSmpName)-strlen(SrwSmpType)-1]
	//RadSamplindIsDefined = 1
endif

SrwPrecSASEName=PrecName[0,strlen(PrecName)-strlen(SrwPrecSASEType)-1]

Variable CntrIsDefined = 0
if((strlen(CntrlName) > 0) %& (cmpstr(CntrlName, "_none_") != 0))
	CntrIsDefined = 1
endif
if(CntrIsDefined != 1)
	SrwSASECntrl("CntrAux",1,1,1,1,1,1,1)
else
	SrwCntrlSASEName=CntrlName[0,strlen(CntrlName)-strlen(SrwCntrlSASEType)-1]
endif

SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

$PrecName[20] = ObsNxNzForProp
$PrecName[21] = ObsNxNzSamplFact

if(strlen(WfrName)==0)
	WfrName=SrwElecName+SrwMagContName+SrwSmpName+SrwPrecSASEName
endif
SrwRadName=WfrName

//SrwSASEElecExtra(ElecName,2,1,10,1) // default extra electron beam parameters for steady state computation
//SrwSASEPrecTime(PrecName,1,1,0.5,0) // default time-dependent precision parameters for steady state computation

SrwWfrPrep(ElecName,RadSmpName,WfrName,0)
WfrName += SrwRadType
SrwRadGenTotName=WfrName

variable NumMacroPart = $PrecName[0]
variable NumSlices = $PrecName[11]
variable CreateElecDistr = $PrecName[15]
if(CreateElecDistr)
	//SrwElecDistrPrep(ElecName, CreateElecDistr + 1, NumMacroPart, NumSlices)
	SrwElecDistrPrep(ElecName, NumMacroPart, NumSlices)
endif
SrwWfrFldUnitSet(WfrName, 0) // Electric field units (0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...)

variable ii
if(CntrIsDefined == 1)
	//display only fundamental harmonic in this function
	ii = 0
	do
		$CntrlName[ii + 1][0] = ""; $CntrlName[ii + 1][1] = ""
		$CntrlName[ii + 8][0] = ""; $CntrlName[ii + 8][1] = ""
		$CntrlName[ii + 15][0] = ""; $CntrlName[ii + 15][1] = ""
		$CntrlName[ii + 22][0] = ""; $CntrlName[ii + 22][1] = ""
		$CntrlName[ii + 33][0] = ""; $CntrlName[ii + 33][1] = ""
		ii += 1
	while(ii < 6)

	SrwSASECntrlZero(CntrlName)
	//If no calculation in time-dependent mode is planned, then 
	//the resulting power vs time, radiation energy and peak radiation power should not be displayed:
	if($PrecName[10] != 2) //itdp
		if(dimsize($CntrlName, 1) < 1)
			redimension/N=(dimsize($CntrlName, 0), 2) $CntrlName
		endif
		
		ii = 0
		do
			$CntrlName[ii + 7][0] = ""; $CntrlName[ii + 7][1] = ""
			$CntrlName[ii + 14][0] = ""; $CntrlName[ii + 14][1] = ""
			$CntrlName[ii + 21][0] = ""; $CntrlName[ii + 21][1] = ""
			ii += 1
		while(ii < 7)
	endif
	//SrwSASECntrlDisplay(CntrlName) //this is called from C
endif

srWfrSASE($ElecName, $MagName, $InRadName, $RadSmpName, $PrecName, $CntrlName, $WfrName)

$ElecName[44] = CreateElecDistr //"Create Electron Distribution: 0 or 1", setting After the action performed
string nmElecDistrib = SrwElecName + "_ebd"
if(CreateElecDistr == 0)
	killwaves/Z $nmElecDistrib
endif

if(CntrIsDefined != 1)
	String DummyCntrName = "CntrAux" + SrwCntrlSASEType
	KillWaves $DummyCntrName
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Start SASE computation (test version, calculates harmonics)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrSASEHarmCreate(WfrName,ElecName,MagName,InRadName,PrecName,CntrlName,NumHarm)
string WfrName="Wfr"
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagContName+SrwMagContType
string InRadName=SrwRadSASEName+SrwRadSASEType
//string RadSmpName=SrwSmpName+SrwSmpType
string PrecName=SrwPrecSASEName+SrwPrecSASEType
string CntrlName=SrwCntrlSASEName+SrwCntrlSASEType
variable NumHarm=srwUtiGetValN("SrwRadSASENumHarm", 1, "")
//variable ObsNxNzForProp=SrwSmpNxNzForProp
//variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt WfrName,"Name of the Wavefront structure"
prompt ElecName,"Electron Beam structure",popup Wavelist("*"+SrwElecType,";","")
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt InRadName,"Input Radiation structure",popup Wavelist("*"+SrwRadSASEType,";","") + Wavelist("*"+SrwRadType,";","")
//prompt RadSmpName,"Radiation Sampling structure (optional)",popup Wavelist("*"+SrwSmpType,";","")
prompt PrecName,"SASE Precision structure",popup Wavelist("*"+SrwPrecSASEType,";","")
prompt CntrlName,"SASE Control structure",popup Wavelist("*"+SrwCntrlSASEType,";","")
prompt NumHarm,"Number of Harmonics to calculate",popup "1;2;3;4;5;6;7"
//prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
//prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	...
//PauseUpdate

if((strlen(ElecName) == 0) %| (cmpstr(ElecName, "_none_") == 0))
	abort "Electron Beam structure is not defined."
endif
if((strlen(MagName) == 0) %| (cmpstr(MagName, "_none_") == 0))
	abort "Magnetic Field structure is not defined."
endif
if((strlen(InRadName) == 0) %| (cmpstr(InRadName, "_none_") == 0))
	abort "Input Radiation structure is not defined."
endif
if((strlen(PrecName) == 0) %| (cmpstr(PrecName, "_none_") == 0))
	abort "Precision Parameters structure is not defined."
endif
if((strlen(CntrlName) == 0) %| (cmpstr(CntrlName, "_none_") == 0))
	abort "Control structure is not defined."
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]
SrwRadSASEName=InRadName[0,strlen(InRadName)-strlen(SrwRadSASEType)-1]

srwUtiSetValN("SrwRadSASENumHarm", NumHarm, "")

//if((strlen(RadSmpName) > 0) %& (cmpstr(RadSmpName, "_none_") != 0))
////SASE computation can be performed even without Observation structure (just for the exit of undulator, with the sampling specified in PrecName)
//	SrwSmpName=RadSmpName[0,strlen(RadSmpName)-strlen(SrwSmpType)-1]
//endif

SrwPrecSASEName=PrecName[0,strlen(PrecName)-strlen(SrwPrecSASEType)-1]

Variable CntrIsDefined = 0
if((strlen(CntrlName) > 0) %& (cmpstr(CntrlName, "_none_") != 0))
	CntrIsDefined = 1
endif
if(CntrIsDefined != 1)
	SrwSASECntrl("CntrAux",1,1,1,1,1,1,1)
else
	SrwCntrlSASEName=CntrlName[0,strlen(CntrlName)-strlen(SrwCntrlSASEType)-1]
endif

//SrwSmpNxNzForProp=ObsNxNzForProp
//SrwSmpNxNzSamplFact=ObsNxNzSamplFact

//$PrecName[20] = ObsNxNzForProp
//$PrecName[21] = ObsNxNzSamplFact

if(strlen(WfrName)==0)
	WfrName=SrwElecName+SrwMagContName+SrwSmpName+SrwPrecSASEName
endif
SrwRadName=WfrName

string nmRadHarmList = WfrName + "_rhl"
string nmRadHarm, nmRadFin
variable iHarm

string dummyRadSmpName = ""
if(NumHarm <= 1)
	//SrwWfrPrep(ElecName,RadSmpName,WfrName,0)
	SrwWfrPrep(ElecName,dummyRadSmpName,WfrName,0)

	WfrName += SrwRadType
	SrwWfrFldUnitSet(WfrName, 0) // Electric field units (0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...)
	nmRadFin = WfrName
else
	make/O/T/N=(NumHarm) $nmRadHarmList
	iHarm = 0
	do
		nmRadHarm = WfrName + "H" + num2str(iHarm + 1)
		//SrwWfrPrep(ElecName,RadSmpName,nmRadHarm,0)
		SrwWfrPrep(ElecName,dummyRadSmpName,nmRadHarm,0)

		nmRadHarm += SrwRadType
		$nmRadHarmList[iHarm] = nmRadHarm
		SrwWfrFldUnitSet(WfrName, 0) // Electric field units (0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...)
		iHarm += 1
	while(iHarm < NumHarm)
	nmRadFin = nmRadHarmList
endif

SrwRadGenTotName=WfrName

variable NumMacroPart = $PrecName[0]
variable NumSlices = $PrecName[11]
variable CreateElecDistr = $PrecName[15]
if(CreateElecDistr)
	SrwElecDistrPrep(ElecName, NumMacroPart, NumSlices)
endif

variable ii
if(CntrIsDefined == 1)

	if(NumHarm < 7)
		ii = NumHarm
		do
			$CntrlName[ii][0] = ""; $CntrlName[ii][1] = ""
			$CntrlName[ii + 7][0] = ""; $CntrlName[ii + 7][1] = ""
			$CntrlName[ii + 14][0] = ""; $CntrlName[ii + 14][1] = ""
			$CntrlName[ii + 21][0] = ""; $CntrlName[ii + 21][1] = ""
			$CntrlName[ii + 32][0] = ""; $CntrlName[ii + 32][1] = ""
			ii += 1
		while(ii < 7)
	endif

	SrwSASECntrlZero(CntrlName)
	//If no calculation in time-dependent mode is planned, then 
	//the resulting power vs time should not be displayed:
	if($PrecName[10] != 2) //itdp
		if(dimsize($CntrlName, 1) < 1)
			redimension/N=(dimsize($CntrlName, 0), 2) $CntrlName
		endif
		
		ii = 0
		do
			$CntrlName[ii + 7][0] = ""; $CntrlName[ii + 7][1] = ""
			$CntrlName[ii + 14][0] = ""; $CntrlName[ii + 14][1] = ""
			$CntrlName[ii + 21][0] = ""; $CntrlName[ii + 21][1] = ""
			ii += 1
		while(ii < 7)
	endif
	//SrwSASECntrlDisplay(CntrlName) //this is called from C
endif

//srWfrSASE($ElecName, $MagName, $InRadName, $RadSmpName, $PrecName, $CntrlName, $WfrName)
srWfrSASE($ElecName, $MagName, $InRadName, $dummyRadSmpName, $PrecName, $CntrlName, $nmRadFin)

$ElecName[44] = CreateElecDistr //"Create Electron Distribution: 0 or 1", setting After the action performed
string nmElecDistrib = SrwElecName + "_ebd"
if(CreateElecDistr == 0)
	killwaves/Z $nmElecDistrib
endif

if(CntrIsDefined != 1)
	String DummyCntrName = "CntrAux" + SrwCntrlSASEType
	KillWaves $DummyCntrName
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Start SASE computation (with electron distribution calculation)
//Accepts Electron Beam input either via SRW structure, or as 
//6D Phase-Space Distribution of Macro-Particles
//Unnecessary; to be removed
//+++++++++++++++++++++++++++++++++++++++
//proc SrwWfrSASECreate(WfrName,ElecName,MagName,InRadName,RadSmpName,PrecName,CntrlName,ObsNxNzForProp,ObsNxNzSamplFact)
proc SrwWfrSASEBeamCreate(WfrName,ElecName,MagName,InRadName,RadSmpName,PrecName,CntrlName,CreateElecDistr,ObsNxNzForProp,ObsNxNzSamplFact)
string WfrName="Wfr"
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagContName+SrwMagContType
string InRadName=SrwRadSASEName+SrwRadSASEType
string RadSmpName=SrwSmpName+SrwSmpType
string PrecName=SrwPrecSASEName+SrwPrecSASEType
string CntrlName=SrwCntrlSASEName+SrwCntrlSASEType
variable CreateElecDistr=srwUtiGetValN("CreateElecDistr", 1, "SrwWfrSASEBeamCreate")
variable ObsNxNzForProp=SrwSmpNxNzForProp
variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt WfrName,"Name of the Wavefront structure"
prompt ElecName,"Electron Beam structure",popup Wavelist("*"+SrwElecType,";","") //+Wavelist("*"+SrwElecDistrType,";","")
prompt MagName,"Magnetic Field structure",popup Wavelist("*"+SrwMagContType,";","")
prompt InRadName,"Input Radiation structure",popup Wavelist("*"+SrwRadSASEType,";","") + Wavelist("*"+SrwRadType,";","")
prompt RadSmpName,"Radiation Sampling structure",popup Wavelist("*"+SrwSmpType,";","")
prompt PrecName,"SASE Precision structure",popup Wavelist("*"+SrwPrecSASEType,";","")
prompt CntrlName,"SASE Control structure",popup Wavelist("*"+SrwCntrlSASEType,";","")
prompt CreateElecDistr, "Create Electron Distribution",popup "No;Yes"
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagContName=MagName[0,strlen(MagName)-strlen(SrwMagContType)-1]
SrwRadSASEName=InRadName[0,strlen(InRadName)-strlen(SrwRadSASEType)-1]
SrwSmpName=RadSmpName[0,strlen(RadSmpName)-strlen(SrwSmpType)-1]
SrwPrecSASEName=PrecName[0,strlen(PrecName)-strlen(SrwPrecSASEType)-1]

Variable CntrIsDefined = 0
if((strlen(CntrlName) > 0) %& (cmpstr(CntrlName, "_none_") != 0))
	CntrIsDefined = 1
endif
if(CntrIsDefined != 1)
	SrwSASECntrl("CntrAux",1,1,1,1,1,1,1)
else
	SrwCntrlSASEName=CntrlName[0,strlen(CntrlName)-strlen(SrwCntrlSASEType)-1]
endif

srwUtiSetValN("CreateElecDistr", CreateElecDistr, "SrwWfrSASEBeamCreate")
SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

$PrecName[20] = ObsNxNzForProp
$PrecName[21] = ObsNxNzSamplFact

if(strlen(WfrName)==0)
	WfrName=SrwElecName+SrwMagContName+SrwSmpName+SrwPrecSASEName
endif
SrwRadName=WfrName

//SrwSASEElecExtra(ElecName,2,1,10,1) // default extra electron beam parameters for steady state computation
////SrwSASEPrecTime(PrecName,1,1,0.5,0) // default time-dependent precision parameters for steady state computation

SrwWfrPrep(ElecName,RadSmpName,WfrName,0)
WfrName += SrwRadType
SrwRadGenTotName=WfrName

variable NumMacroPart = $PrecName[0]
variable NumSlices = $PrecName[11]
CreateElecDistr -= 1
if(CreateElecDistr)
	SrwElecDistrPrep(ElecName, CreateElecDistr, NumMacroPart, NumSlices)
//else
//	$ElecName[44] = 0
endif

SrwWfrFldUnitSet(WfrName, 0) // Electric field units (0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...)

if(CntrIsDefined == 1)
	SrwSASECntrlZero(CntrlName)
	//If no calculation in time-dependent mode is planned, then 
	//the resulting power vs time should not be displayed:
	if($PrecName[10] != 2) //itdp
		$CntrlName[10][0] = ""
		$CntrlName[10][1] = ""
		$CntrlName[11][0] = ""
		$CntrlName[11][1] = ""
		$CntrlName[12][0] = ""
		$CntrlName[12][1] = ""
	endif
	//SrwSASECntrlDisplay(CntrlName) //this is called from C
endif

srWfrSASE($ElecName, $MagName, $InRadName, $RadSmpName, $PrecName, $CntrlName, $WfrName)

$ElecName[44] = CreateElecDistr //"Create Electron Distribution: 0 or 1", setting After the action performed
string nmElecDistrib = SrwElecName + "_ebd"
if(CreateElecDistr == 0)
	killwaves/Z $nmElecDistrib
endif

if(CntrIsDefined != 1)
	String DummyCntrName = "CntrAux" + SrwCntrlSASEType
	KillWaves $DummyCntrName
endif
end
