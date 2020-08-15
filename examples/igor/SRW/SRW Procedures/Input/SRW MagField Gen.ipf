
//+++++++++++++++++++++++++++++++++++++++
//
//Duplicate the Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagDupl(MagName,Name)
String MagName=SrwMagName+SrwFieldType
String Name=SrwMagName+"d"
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt Name,"Name of the Duplicated Field structure"
Silent 1						|	Duplicating the field structure  ...
PauseUpdate
//SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwMagName=Name

String BxName=Name+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType
String BzName=Name+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType
//SrwMagBname=MagName+SrwSuffixMagField+SrwSuffixZ
SrwMagBname=BzName[0,strlen(BzName)-strlen(SrwFieldWaveType)-1]

Name=Name+SrwFieldType
SrwMagGenTotName=Name

duplicate/O $($MagName[0]) $BxName
duplicate/O $($MagName[1]) $BzName
duplicate/O $MagName $Name

$Name[0]=BxName
$Name[1]=BzName
end

//+++++++++++++++++++++++++++++++++++++++
//
//Set up precision parameters for SR computation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagPrec(MagName,Mode,Step,Prec,MaxPtsToSave,UseDiffLimits,sdep,sfin)
String MagName=SrwMagName+SrwFieldType
Variable Mode=SrwMode
Variable Step=SrwRadIntStep
Variable Prec=SrwPrec
Variable MaxPtsToSave=SrwMaxPtsToSave
Variable UseDiffLimits=SrwUseDiffRadIntLimits
Variable sdep=SrwSdep
Variable sfin=SrwSfin
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt Mode,SrwPMode,popup SrwPOPUPMode
prompt Step,SrwPRadIntStep
prompt Prec,SrwPPrec
prompt MaxPtsToSave,SrwPMaxPtsToSave
prompt UseDiffLimits,SrwPUseDiffRadIntLimits,popup SrwPOPUPUseDiffRadIntLimits
prompt sdep,SrwPSdep
prompt sfin,SrwPsfin

Silent 1						|	Modifying the Field structure  ....
PauseUpdate

SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwMode=Mode
SrwRadIntStep=Step
SrwPrec=Prec
SrwMaxPtsToSave=MaxPtsToSave
SrwUseDiffRadIntLimits=UseDiffLimits
SrwSdep=sdep
SrwSfin=sfin

if(Mode==1) // man.
	Prec = Step
endif
if(UseDiffLimits==1) // "No"
	sdep = 0.
	sfin = 0.
endif

$MagName[2]=num2str(Mode-1)   // Mode of Rad. Integration (0- man; 1- auto und.; 2- auto wigg.)
$MagName[3]=num2str(Prec)   // Rel. Prec. (for "auto" modes) or Step Size in m (for "man" mode)
$MagName[4]=num2str(sdep)   // sdep
$MagName[5]=num2str(sfin)   // sfin
$MagName[7]=num2str(MaxPtsToSave)   // Max. number of points to keep in memory at SR comp.

end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagFieldCreate(MagName,FieldCenter,FieldLength,FieldNpts)
String MagName=SrwMagName
Variable FieldCenter=SrwFieldCenter
Variable FieldLength=SrwFieldLength
Variable FieldNpts=SrwFieldNpts
prompt MagName,SrwPMagName
prompt FieldCenter,SrwPFieldCenter
prompt FieldLength,SrwPFieldLength
prompt FieldNpts,SrwPFieldNpts

Silent 1						|	Creating Magnetic Field Waves  ...
PauseUpdate

SrwFieldCenter=FieldCenter
//SrwMagScent=FieldCenter

SrwFieldLength=FieldLength
SrwFieldNpts=FieldNpts
SrwMagName=MagName
SrwMagGenTotName=MagName+SrwFieldType

String BxName=MagName+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType
String BzName=MagName+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType
SrwMagBname=MagName+SrwSuffixMagField+SrwSuffixZ

Make/N=(FieldNpts)/D/O $BxName, $BzName
SetScale/I x FieldCenter-FieldLength/2, FieldCenter+FieldLength/2, "m", $BxName, $BzName
SetScale d 0, 0, SrwPUnitMagField, $BxName, $BzName

MagName+=SrwFieldType
Make/N=8/T/O  $MagName
//SetScale d 0, 1E+23, SrwMagType_TabTrUnif, $MagName

$MagName[0]=BxName
$MagName[1]=BzName
$MagName[2]=num2str(SrwMode-1)   // Mode of Rad. Integration (0- man; 1- auto und.; 2- auto wigg.)
$MagName[3]=num2str(SrwPrec)   // Rel. Prec. (for "auto" modes) or Step Size in m (for "man" mode)
$MagName[4]=num2str(0)   // sdep
$MagName[5]=num2str(0)   // sfin

//$MagName[6]=num2str(FieldElecS0)   // s0 where particle traj. transverse coordinates and angles are zero (term of the {E,I,{s0,x0,x10,z0,z10}} elec. input wave for C routine)

$MagName[7]=num2str(SrwMaxPtsToSave)   // Max. number of points to keep in memory

end  // proc SrwMagFieldCreate

//+++++++++++++++++++++++++++++++++++++++
//
//Create 3D Magnetic Field structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagField3dCreate(MagName,FieldCenS,FieldLenS,FieldNpS,FieldCenX,FieldLenX,FieldNpX,FieldCenZ,FieldLenZ,FieldNpZ)
String MagName=SrwMagName
Variable FieldCenS=srwUtiGetValN("SrwFieldCenter", 0, "")
Variable FieldLenS=srwUtiGetValN("SrwFieldLength", 2., "")
Variable FieldNpS=srwUtiGetValN("SrwFieldNpts", 1000, "")
Variable FieldCenX=srwUtiGetValN("SrwFieldCenX", 0, "")
Variable FieldLenX=srwUtiGetValN("SrwFieldLenX", 0.1, "")
Variable FieldNpX=srwUtiGetValN("SrwFieldNpX", 10, "")
Variable FieldCenZ=srwUtiGetValN("SrwFieldCenZ", 0, "")
Variable FieldLenZ=srwUtiGetValN("SrwFieldLenZ", 0.1, "")
Variable FieldNpZ=srwUtiGetValN("SrwFieldNpZ", 10, "")
prompt MagName,SrwPMagName
prompt FieldCenS,"Longitudinal Center Position"
prompt FieldLenS,"Range of Longitudinal Position"
prompt FieldNpS,"Number of Points vs Longitudinal Position"
prompt FieldCenX,"Horizontal Center Position"
prompt FieldLenX,"Range of Horizontal Position"
prompt FieldNpX,"Number of Points vs Horizontal Position"
prompt FieldCenZ,"Vertical Center Position"
prompt FieldLenZ,"Range of Vertical Position"
prompt FieldNpZ,"Number of Points vs Vertical Position"
Silent 1						|	Creating Magnetic Field Waves  ...
PauseUpdate

SrwFieldCenter=FieldCenS
SrwFieldLength=FieldLenS
SrwFieldNpts=FieldNpS
srwUtiSetValN("SrwFieldCenX", FieldCenX, "")
srwUtiSetValN("SrwFieldLenX", FieldLenX, "")
srwUtiSetValN("SrwFieldNpX", FieldNpX, "")
srwUtiSetValN("SrwFieldCenZ", FieldCenZ, "")
srwUtiSetValN("SrwFieldLenZ", FieldLenZ, "")
srwUtiSetValN("SrwFieldNpZ", FieldNpZ, "")

String/G SrwMagField3dType = "_m3d", SrwMagField3dWaveType = "_f3d", SrwSuffixS = "S"

SrwMagName=MagName
SrwMagGenTotName=MagName+SrwMagField3dType

String BsName=MagName+SrwSuffixMagField+SrwSuffixS+SrwMagField3dWaveType
String BxName=MagName+SrwSuffixMagField+SrwSuffixX+SrwMagField3dWaveType
String BzName=MagName+SrwSuffixMagField+SrwSuffixZ+SrwMagField3dWaveType
//SrwMagBname=MagName+SrwSuffixMagField+SrwSuffixZ
SrwMagBname=BzName[0,strlen(BzName)-strlen(SrwMagField3dWaveType)-1]

Make/N=(FieldNpS,FieldNpX,FieldNpZ)/D/O $BsName, $BxName, $BzName
SetScale/I x FieldCenS-FieldLenS/2, FieldCenS+FieldLenS/2, "m", $BsName, $BxName, $BzName
SetScale/I y FieldCenX-FieldLenX/2, FieldCenX+FieldLenX/2, "m", $BsName, $BxName, $BzName
SetScale/I z FieldCenZ-FieldLenZ/2, FieldCenZ+FieldLenZ/2, "m", $BsName, $BxName, $BzName
SetScale d 0, 0, SrwPUnitMagField, $BsName, $BxName, $BzName

MagName=SrwMagGenTotName
Make/N=3/T/O  $MagName
//$MagName[0]=BsName
//$MagName[1]=BxName
//$MagName[2]=BzName

$MagName[0]=BxName //for compatibility with transversely-uniform field (*_mag)
$MagName[1]=BzName
$MagName[2]=BsName

end  // proc SrwMagField3dCreate

//+++++++++++++++++++++++++++++++++++++++
//
//Returns longitudinal center position of the magnetic field
//
//+++++++++++++++++++++++++++++++++++++++
function srwGetMagFldCenter(MagName)
string MagName //=srwUtiGetValS("SrwMagGenTotName", "", "")
SVAR SrwFieldType, SrwUndType, SrwMagConstType
string BufMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
string MagType = MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1]
if(cmpstr(MagType,SrwFieldType)==0)
	srwUtiSetValS("SrwMagName", BufMagName, "")
	
	wave/T wMag = $MagName
	string BxName = wMag[0]
	wave wBx = $BxName
	variable sStart = DimOffset(wBx, 0)
	variable sDelta = DimDelta(wBx, 0)
	variable ns = DimSize(wBx, 0)
	return sStart + 0.5*sDelta*(ns - 1)
endif
if(cmpstr(MagType,SrwUndType)==0)
	srwUtiSetValS("SrwUndName", BufMagName, "")
	return 0.
endif
if(cmpstr(MagType,SrwMagConstType)==0)
	srwUtiSetValS("SrwMagConstName", BufMagName, "")
	return 0.
endif
DoAlert 0, "No Magnetic Field structure found"
end

//+++++++++++++++++++++++++++++++++++++++
//
//Add a sinus to a Magnetic Field Wave
//Copied from "B2E"
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagSin(FldName, zero, scent, period, taper, sag, nperiod, bpeak, bnoise, bnoisei)
String FldName=SrwMagBname+SrwFieldWaveType
Variable zero=SrwMagZeroSgn
Variable scent=SrwMagScent
Variable period=SrwPeriod
Variable nperiod=SrwMagNper
Variable bpeak=SrwMagBpeak
Variable taper=SrwMagTaper
Variable sag=SrwMagSag
Variable bnoise=SrwMagBnoise
Variable bnoisei=SrwMagBnoisei
Prompt FldName,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Prompt zero,SrwPMagZero,popup "Yes;No"
Prompt scent,SrwPMagScent
Prompt period,SrwPPeriod
Prompt nperiod,SrwPMagNper
Prompt bpeak,SrwPMagBpeak
Prompt bnoise,SrwPMagBnoise
Prompt bnoisei,SrwPMagBnoisei
Prompt taper,SrwPMagTaper
Prompt sag,SrwPMagSag
PauseUpdate
Silent 1						|	Add a SineWave...

if(cmpstr(FldName,"_none_")==0)
	SrwMagFieldCreate()
	SrwMagSin()
	return
endif

SrwMagBname=FldName[0,strlen(FldName)-strlen(SrwFieldWaveType)-1]
SrwMagZeroSgn=zero
SrwMagScent=scent
SrwPeriod=period
SrwMagNper=nperiod
SrwMagBpeak=bpeak
SrwMagTaper=taper
SrwMagSag=sag
SrwMagBnoise=bnoise
SrwMagBnoisei=bnoisei

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

if(zero==1)
$FldName=0
endif
Variable pinit, pend, n
Variable a,b,c,sfin,sdep,dep,fin
period /= 1000				//[mm] -> [m]
taper /=1000
sag /=1000
if (nperiod<0)
scent=(fin+dep)/2
dep=leftx($FldName);fin=rightx($FldName)
nperiod=floor((fin-dep+2*nperiod/1000)/period)
endif
//lwv=FldName
sdep=scent-nperiod/2*period
sfin=sdep+nperiod*period
a=taper/(sfin-sdep)
b=-taper/2
c=sag*4/(sfin-sdep)/(sfin-sdep)
pinit=x2pnt($FldName,sdep)
pend=x2pnt($FldName,sfin)-1
Duplicate/O $FldName, tamp
tamp=bpeak*Sin(2*Pi*(x-sdep)/period)
tamp *=exp(-Pi*(a*(x-sdep)+b+c*(x-sdep)*(x-sfin))/period)

String AbortMes = "Out of Magnetic Field definition range."

Make/O/N=2 w1crenau
w1crenau = 0

if(pend-pinit-1 <= 0) 
	abort(AbortMes)
endif
//InsertPoints 2,pend-pinit-1,w1crenau
InsertPoints 2,pend-pinit-2,w1crenau //oc

w1crenau=0.5
if(pinit <= 0) 
	abort(AbortMes)
endif
//InsertPoints 0,pinit,w1crenau
InsertPoints 0,pinit+1,w1crenau //oc

if(numpnts($FldName)-pend-1 <= 0) 
	abort(AbortMes)
endif
//InsertPoints pend+1,numpnts($FldName)-pend-1,w1crenau
InsertPoints pend+2,numpnts($FldName)-pend-1,w1crenau //oc

pinit=x2pnt($FldName,sdep+period/2)
pend=x2pnt($FldName,sfin-period/2)-1
Make/O/N=2 w2crenau
//InsertPoints 2,pend-pinit-1,w2crenau
InsertPoints 2,pend-pinit-2,w2crenau //oc

w2crenau=0.5
//InsertPoints 0,pinit,w2crenau
InsertPoints 0,pinit+1,w2crenau //oc

//InsertPoints pend+1,numpnts($FldName)-pend-1,w2crenau
InsertPoints pend+2,numpnts($FldName)-pend-1,w2crenau //oc

w1crenau+=w2crenau
CopyScales $FldName, w1crenau
tamp*=w1crenau[p]
$FldName+= tamp
if (bnoise>0)
$FldName += gnoise(bnoise)*srwUtiDataWindowZero(x, sdep, sfin)
endif

if (bnoisei>0)
Duplicate $FldName a1
a1=gnoise(bnoisei)
Integrate a1
$FldName+=a1
Killwaves a1
endif
Killwaves/Z w1crenau, w2crenau,tamp

End  // proc SrwMagSin(FldName, zero, scent, period, taper, sag, nperiod, bpeak, bnoise, bnoisei)

//+++++++++++++++++++++++++++++++++++++++
//
//Zero Magnetic Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagZero(Bname)
string Bname=SrwMagBname+SrwFieldWaveType
prompt Bname,SrwPMagBname,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Silent 1						|	Iniatialize the Field Component  ....
PauseUpdate

if(cmpstr(Bname,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagZero();
	Return;
endif

SrwMagBname=Bname[0,strlen(Bname)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

$Bname=0
end

//+++++++++++++++++++++++++++++++++++++++
//
//Constant Magnetic Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagConst(Bname, Bconst)
string Bname=SrwMagBname+SrwFieldWaveType;
variable Bconst=SrwMagBconst;
prompt Bname,SrwPMagBname,popup Wavelist("*"+SrwFieldWaveType ,";", "");
prompt Bconst,SrwPMagBconst;
Silent 1						|	Iniatializing Field Component  ...
PauseUpdate;

if(cmpstr(Bname,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagConst();
	Return;
endif

SrwMagBname=Bname[0,strlen(Bname)-strlen(SrwFieldWaveType)-1]
SrwMagBconst=Bconst

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

$Bname=Bconst;
end

//+++++++++++++++++++++++++++++++++++++++
//
//Display Magnetic Field Component
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagDisplayField(Bname)
string  Bname=SrwMagBname+SrwFieldWaveType
prompt Bname,SrwPMagBname,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Silent 1						|	Initialize the Field Component  ....
PauseUpdate

SrwMagBname=Bname[0,strlen(Bname)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

display  $Bname
Label bottom SrwPLabelLongPos
Label left SrwPLabelMagField

end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Fringe Field of a Bending Magnet Edge
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagEdge(FldName, zero, scent, StrSectLen, FringeSize, bconst)
String FldName=SrwMagBname+SrwFieldWaveType
Variable zero=SrwMagZeroSgn
Variable scent=SrwMagScentEdge
Variable StrSectLen=SrwMagStrSectLen
Variable FringeSize=SrwMagFringeSize
Variable  bconst=SrwMagBconst
Prompt FldName,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Prompt zero,SrwPMagZero,popup "Yes;No"
Prompt scent,SrwPMagScent
Prompt StrSectLen,SrwPMagStrSectLen
Prompt FringeSize,SrwPMagFringeSize
Prompt bconst,SrwPMagBconst
PauseUpdate
Silent 1						|	Adding a Fringe Field ...

if(cmpstr(FldName,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagEdge();
	Return;
endif

SrwMagZeroSgn=zero
SrwMagScentEdge=scent
SrwMagStrSectLen=StrSectLen
SrwMagFringeSize=FringeSize
SrwMagBconst=bconst

SrwMagBname=FldName[0,strlen(FldName)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

if(zero==1)
	$FldName=0
endif
Variable Con1=0.22756
Variable ds=Con1*FringeSize*0.001
$FldName += sREdgeFieldExp(x, scent, StrSectLen, ds, bconst)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create a Dipole Magnet Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagDipole(FldName, Zero, Scent, Len, FringeSize, Bconst)
String FldName=SrwMagBname+SrwFieldWaveType
Variable Zero=srwUtiGetValN("SrwMagZeroSgn", 1, "")
Variable Scent=srwUtiGetValN("SrwMagScentEdge", 0, "")
Variable Len=srwUtiGetValN("SrwMagDipoleLen", 2, "")
Variable FringeSize=srwUtiGetValN("SrwMagFringeSize", 50, "")
Variable  Bconst=srwUtiGetValN("SrwMagBconst", 1.56, "")
Prompt FldName,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Prompt Zero,SrwPMagZero,popup "Yes;No"
Prompt Scent,SrwPMagScent
Prompt Len,"Iron Length [m]"
Prompt FringeSize,SrwPMagFringeSize
Prompt Bconst,SrwPMagBconst
PauseUpdate
Silent 1						|	Adding a Dipole Field ...

if(cmpstr(FldName,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagDipole();
	Return;
endif

srwUtiSetValN("SrwMagZeroSgn", Zero, "") //SrwMagZeroSgn=Zero
srwUtiSetValN("SrwMagScentEdge", Scent, "") //SrwMagScentEdge=Scent
srwUtiSetValN("SrwMagDipoleLen", Len, "") //SrwMagDipoleLen=Len
srwUtiSetValN("SrwMagFringeSize", FringeSize, "") //SrwMagFringeSize=FringeSize
srwUtiSetValN("SrwMagBconst", Bconst, "") //SrwMagBconst=Bconst

SrwMagBname=FldName[0,strlen(FldName)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

if(Zero==1)
	$FldName=0
endif
Variable Con1=0.22756
Variable ds=Con1*FringeSize*0.001
$FldName += sRDipoleFieldExp(x, Scent, Len, ds, Bconst)
end

//+++++++++++++++++++++++++++++++++++++++
function sREdgeFieldExp(s, s0, L, ds, b)
variable s, s0, L, ds, b
return b*(1./(1.+exp((s-(s0-0.5*L))/ds)) + 1./(1.+exp(-(s-(s0+0.5*L))/ds)))
end

//+++++++++++++++++++++++++++++++++++++++
function sRDipoleFieldExp(s, s0, L, ds, b)
variable s, s0, L, ds, b
return b-sREdgeFieldExp(s, s0, L, ds, b)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Steerer Field of Gaussian Longitudinal Profile
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagGsnAng(FldName, zero, scent, sigma, FieldInt)
String FldName=SrwMagBname+SrwFieldWaveType
Variable zero=SrwMagZeroSgn
Variable scent=SrwMagScentGsnAng
Variable sigma=SrwMagSigmaGsnAng
Variable FieldInt=SrwFieldIntGsnAng
Prompt FldName,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "")
Prompt zero,SrwPMagZero,popup "Yes;No"
Prompt scent,SrwPMagScent
Prompt sigma,SrwPMagSigmaGsnAng
Prompt FieldInt,SrwPFieldIntGsnAng
PauseUpdate
Silent 1						|	Adding a Fringe Field ...

if(cmpstr(FldName,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagGsnAng();
	Return;
endif

SrwMagZeroSgn=zero
SrwMagScentGsnAng=scent
SrwMagSigmaGsnAng=sigma
SrwFieldIntGsnAng=FieldInt

SrwMagBname=FldName[0,strlen(FldName)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

if(zero==1)
	$FldName=0
endif
$FldName += sRGssn(x, scent, sigma*0.001, FieldInt*0.001)
End

//+++++++++++++++++++++++++++++++++++++++
Function sRGssn(s, s0, sigma, intgr)
Variable s, s0, sigma, intgr
Variable t=(s-s0)/sigma
return (intgr*0.3989422804/sigma)*exp(-0.5*t*t)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Importing Magnetic Field component data
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagImportCmpn(name,units,fname)
String name=SrwMagBname+SrwFieldWaveType;
Variable units=SrwMagImportUnits;
String fname=SrwMagImportFileName;
Prompt name,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "");
Prompt units,SrwPMagImportUnits,popup SrwPOPUPMagImportUnits;
Prompt fname,SrwPMagImportFileName
Silent 1						|	Importing the Field Component ...

if(cmpstr(name,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagImportCmpn();
	Return;
endif

SrwMagBname=name[0,strlen(name)-strlen(SrwFieldWaveType)-1]

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

SrwMagImportUnits=units
SrwMagImportFileName=fname

LoadWave/G/D/O/A=wMagImportAux/N/Q fname;

String FailedStr="Import Failed";
String AdviceStr="Try to modify the input file or import it manually as explained in the SRW Help"

Print "";

if(V_flag != 1)
	Print FailedStr;
	Print "Invalid File Format";
	Print AdviceStr;
	return;
endif
Variable NumPtsImp = numpnts(wMagImportAux0), NumPtsDef = DimSize($name,0);
if(NumPtsImp  != NumPtsDef)
	Print FailedStr;
	Print "Incorrect Number of Points";
	Print "	Expected:", NumPtsDef;
	Print "	Input:", NumPtsImp;
	Print AdviceStr;
	return;
endif

Print "Import Successful";

SetScale/P x DimOffset($name, 0), DimDelta($name, 1), wMagImportAux0;
$name = wMagImportAux0;

if(units == 2) // Gauss
	$name *= 0.0001;
endif

KillWaves/Z wMagImportAux0;
End

//+++++++++++++++++++++++++++++++++++++++
//
//Importing Magnetic Field component data and add it 
//to existing field
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwMagAddImpCmpn(name,zero,s0,step,units,fname)
String name=SrwMagBname+SrwFieldWaveType;
Variable zero=SrwMagZeroSgn
Variable s0=srwUtiGetValN("s0", 0., "SrwMagAddImpCmpn")
Variable step=srwUtiGetValN("step", 1., "SrwMagAddImpCmpn")
Variable units=SrwMagImportUnits;
String fname=SrwMagImportFileName;
Prompt zero,SrwPMagZero,popup "Yes;No"
Prompt name,SrwPFldName,popup Wavelist("*"+SrwFieldWaveType ,";", "");
Prompt units,SrwPMagImportUnits,popup SrwPOPUPMagImportUnits;
Prompt s0,"Longitudinal Center Position [m]"
Prompt step,"Step [mm]"
Prompt fname,SrwPMagImportFileName
PauseUpdate
Silent 1						|	Importing the Field Component ...

if(cmpstr(name,"_none_")==0)
	SrwMagFieldCreate();
	SrwMagImportCmpn();
	Return;
endif

SrwMagBname=name[0,strlen(name)-strlen(SrwFieldWaveType)-1]
SrwMagZeroSgn=zero
SrwMagImportUnits=units
SrwMagImportFileName=fname
srwUtiSetValN("s0", s0, "SrwMagAddImpCmpn")
srwUtiSetValN("step", step, "SrwMagAddImpCmpn")

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwFieldType

LoadWave/G/D/O/A=wMagImportAux/N/Q fname

String FailedStr="Import Failed"
String AdviceStr="Try to modify the input file or import it manually as explained in the SRW Help"

if(V_flag != 1)
	Print FailedStr
	Print "Invalid File Format"
	Print AdviceStr
	return;
endif
Variable NumPtsImp = numpnts(wMagImportAux0)
Variable RangeImp = (NumPtsImp - 1)*step*0.001
Variable OffsetImp = s0 - 0.5*RangeImp

SetScale/P x OffsetImp, step*0.001, wMagImportAux0

if(zero == 1)
	$name = 0.
endif

$name += wMagImportAux0(x)*srwUtiStep(x - OffsetImp)*srwUtiStep(OffsetImp + RangeImp - x)

if(units == 2) // Gauss
	$name *= 0.0001
endif

//Print "Import Successful"
KillWaves/Z wMagImportAux0
End

//+++++++++++++++++++++++++++++++++++++++
//
//Importing Magnetic Field component data and add it 
//to existing field
//For the moment, it imports only 2D field components
//i.e. those depending on s and x
//To improve !!!
//+++++++++++++++++++++++++++++++++++++++
Proc SrwMag3dAddImpCmpn(name,zero,sb,step,xb,xstep,units,fname)
//Proc SrwMag3dAddImpCmpn(name,zero,sb,step,xb,xstep,zb,zstep,units,fname)
String name=SrwMagBname+srwUtiGetValS("SrwMagField3dWaveType", "_f3d", "")
Variable zero=SrwMagZeroSgn
Variable sb=srwUtiGetValN("sb", 0., "SrwMag3dAddImpCmpn")
Variable step=srwUtiGetValN("step", 1., "SrwMag3dAddImpCmpn")
Variable xb=srwUtiGetValN("xb", 0., "SrwMag3dAddImpCmpn")
Variable xstep=srwUtiGetValN("xstep", 1., "SrwMag3dAddImpCmpn")
//Variable zb=srwUtiGetValN("zb", 0., "SrwMag3dAddImpCmpn")
//Variable zstep=srwUtiGetValN("zstep", 1., "SrwMag3dAddImpCmpn")
Variable units=SrwMagImportUnits
String fname=SrwMagImportFileName
Prompt zero,SrwPMagZero,popup "Yes;No"
Prompt name,SrwPFldName,popup Wavelist("*"+srwUtiGetValS("SrwMagField3dWaveType", "_f3d", "") ,";", "")
Prompt sb,"Longitudinal Start Position [m]"
Prompt step,"Longitudinal Step [m]"
Prompt xb,"Horizontal Start Position [m]"
Prompt xstep,"Horizontal Step [m]"
//Prompt zb,"Vertical Start Position [m]"
//Prompt zstep,"Vertical Step [m]"
Prompt units,SrwPMagImportUnits,popup SrwPOPUPMagImportUnits
Prompt fname,SrwPMagImportFileName
PauseUpdate
Silent 1						|	Importing the Field Component ...

if(cmpstr(name,"_none_")==0)
	//SrwMagFieldCreate()
	//SrwMagImportCmpn()
	Return
endif

SrwMagBname=name[0,strlen(name)-strlen(SrwMagField3dWaveType)-1]
SrwMagZeroSgn=zero
SrwMagImportUnits=units
SrwMagImportFileName=fname
srwUtiSetValN("sb", sb, "SrwMag3dAddImpCmpn")
srwUtiSetValN("step", step, "SrwMag3dAddImpCmpn")
srwUtiSetValN("xb", xb, "SrwMag3dAddImpCmpn")
srwUtiSetValN("xstep", xstep, "SrwMag3dAddImpCmpn")
//srwUtiSetValN("zb", zb, "SrwMag3dAddImpCmpn")
//srwUtiSetValN("zstep", zstep, "SrwMag3dAddImpCmpn")

SrwMagName=SrwMagBname[0,strlen(SrwMagBname)-strlen(SrwSuffixMagField)-strlen(SrwSuffixX)-1]
SrwMagGenTotName=SrwMagName+SrwMagField3dType //SrwFieldType

//LoadWave/G/D/O/N=wMagImportAux/Q fname
LoadWave/J/K=1/D/O/N=wMagImportAux/Q/M/U={0, 0, 0, 0}  fname

String FailedStr="Import Failed"
String AdviceStr="Try to modify the input file or import it manually as explained in the SRW Help"

//Print "";

if(V_flag != 1)
	Print FailedStr
	Print "Invalid File Format"
	Print AdviceStr
	return
endif

//Variable NumPtsImp = numpnts(wMagImportAux0)
//Variable RangeImp = (NumPtsImp - 1)*step*0.001
//Variable OffsetImp = s0 - 0.5*RangeImp

Variable nsImp = dimsize(wMagImportAux0, 0)
Variable rsImp = (nsImp - 1)*step

Variable nxImp = dimsize(wMagImportAux0, 1)
Variable rxImp = (nxImp - 1)*xstep

SetScale/P x sb, step, wMagImportAux0
if(nxImp > 0)
	SetScale/P y xb, xstep, wMagImportAux0
endif

if(zero == 1)
	$name = 0.
endif

//$name += wMagImportAux0(x)*srwUtiStep(x - OffsetImp)*srwUtiStep(OffsetImp + RangeImp - x)
if(nxImp > 0)
	$name += wMagImportAux0(x)(y)*srwUtiStep(x - sb)*srwUtiStep(sb + rsImp - x)*srwUtiStep(y - xb)*srwUtiStep(xb + rxImp - y)
else
	$name += wMagImportAux0(x)*srwUtiStep(x - sb)*srwUtiStep(sb + rsImp - x)
endif

if(units == 2) // Gauss
	$name *= 0.0001
endif

//Print "Import Successful"
KillWaves/Z wMagImportAux0
End

//+++++++++++++++++++++++++++++++++++++++
//
// Import Magnetic Field Wave
// The field must be in Tesla
// the data must be equidistant
// the file must be in text with two coulmns, the first is the horizontal
// field, the second is the vertical field
// the delimiter can be a space or a tab or comma
//
// propose two units Gauss or Tesla
//
//+++++++++++++++++++++++++++++++++++++++
Proc SRMagImport(MagName,FieldCenter,FieldLength)
String MagName=SrwMagName
Variable FieldCenter=SrwFieldCenter
Variable FieldLength=SrwFieldLength
prompt MagName,SrwPMagName
prompt FieldCenter,SrwPFieldCenter
prompt FieldLength,SrwPFieldLength

Silent 1      | Creating Magnetic Field Waves  ...
PauseUpdate

SrwMagFieldCreate(MagName,FieldCenter,FieldLength,2)

string aa
aa ="Your Magnetic Field data should be in the proper format (see the SRW documentation)."
aa+="Select a text file containing both the horizontal and vertical components of the mangetic field to be imported"
DoAlert  0  aa

SrwMagImportAlert();

LoadWave/G/D/O/A=Field/N/Q

variable info=0;

if ( v_flag  != 2)
print "Invalid file format"
info=1
endif

if ( numpnts(Field0)  != numpnts(Field1))
print "Horizontal and Vertical Fields must have the same number of points "
info=1
endif


if (info==0)
  print "Successfull"
  
SetScale/I x FieldCenter-FieldLength/2, FieldCenter+FieldLength/2, "m",
Field0, Field1

String BxName=MagName+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType
String BzName=MagName+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType
Redimension/N=(numpnts(Field0)) $BxName $BzName
SetScale/I x FieldCenter-FieldLength/2, FieldCenter+FieldLength/2, "m",
$BxName, $BzName
$BxName=Field0
$BzName=Field1

 else
  print "Unsuccessfull."
  print "Modidy your input file or try importing the Field manually as
explained in the documentation"
 endif

end  

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Equilibrium Trajectory
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagElecTraj(ElecName, MagName, disxang, disxtra, diszang, disztra)
String ElecName=SrwElecName+SrwElecType
String MagName=SrwMagName+SrwFieldType
Variable disxang=SrwDisplayAngX, disxtra=SrwDisplayTraX, diszang=SrwDisplayAngZ, disztra=SrwDisplayTraZ
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt disxang,"Display x' vs s?",popup "No;Yes"
prompt disxtra,"Display x vs s?",popup "No;Yes"
prompt diszang,"Display z' vs s?",popup "No;Yes"
prompt disztra,"Display z vs s?",popup "No;Yes"
Silent 1						|	Computing the Trajectory  ...
PauseUpdate

if(cmpstr(ElecName,"_none_")==0)
	Abort SrwPAlertNoElec
endif
if(cmpstr(MagName,"_none_")==0)
	Abort SrwPAlertNoMag
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwDisplayAngX=disxang
SrwDisplayTraX=disxtra
SrwDisplayAngZ=diszang
SrwDisplayTraZ=disztra

String xAng=SrwElecName+SrwMagName+SrwSuffixX+SrwAngleWaveType
String xTra=SrwElecName+SrwMagName+SrwSuffixX+SrwTrajWaveType
String zAng=SrwElecName+SrwMagName+SrwSuffixZ+SrwAngleWaveType
String zTra=SrwElecName+SrwMagName+SrwSuffixZ+SrwTrajWaveType

Duplicate/O $($MagName[0]) $xAng
SetScale d 0,0,"r", $xAng
Duplicate/O $xAng $xTra
SetScale d 0,0,"m", $xTra
Duplicate/O $($MagName[1]) $zAng
SetScale d 0,0,"r", $zAng
Duplicate/O $zAng $zTra
SetScale d 0,0,"m", $zTra

srTraj($ElecName, $MagName, $xAng, $xTra, $zAng, $zTra)

if(disxang==2)
	Display $xAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorAngle
endif
if(disxtra==2)
	Display $xTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorPos
endif
if(diszang==2)
	Display $zAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerAngle
endif
if(disztra==2)
	Display $zTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerPos
endif
end  // proc SrwMagElecTraj

//+++++++++++++++++++++++++++++++++++++++
//
// Trajectory structure: Init globals (remove when stabilized)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwTrjInit()

String/G SrwTrjType=SrwSeparator+"trj"
String/G SrwTrjName=""; String/G SrwPTrjName="Name of the Trajectory structure", SrwPTrjName2="Trajectory structure"

end

//+++++++++++++++++++++++++++++++++++++++
//
// Create Trajectory structure and compute trajectory in Transversely Uniform field 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwTrjCreateTransvUnif(TrjName, ElecName, MagName, disxang, disxtra, diszang, disztra)
String TrjName=SrwElecName+SrwMagName
String ElecName=SrwElecName+SrwElecType
String MagName=SrwMagName+SrwFieldType
Variable disxang=SrwDisplayAngX, disxtra=SrwDisplayTraX, diszang=SrwDisplayAngZ, disztra=SrwDisplayTraZ
prompt TrjName,SrwPTrjName
prompt ElecName,SrwPElecName2,popup wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup wavelist("*"+SrwFieldType,";","")
prompt disxang,"Display x' vs s?",popup "No;Yes"
prompt disxtra,"Display x vs s?",popup "No;Yes"
prompt diszang,"Display z' vs s?",popup "No;Yes"
prompt disztra,"Display z vs s?",popup "No;Yes"
Silent 1						|	Computing the Trajectory  ...
PauseUpdate

if(strlen(TrjName)==0)
	abort SrwPAlertBadName
endif
if(strlen(TrjName)>28)
	abort SrwPAlertTooLongName
endif

if(cmpstr(ElecName,"_none_")==0)
	abort SrwPAlertNoElec
endif
if(cmpstr(MagName,"_none_")==0)
	abort SrwPAlertNoMag
endif

SrwTrjName=TrjName
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwDisplayAngX=disxang
SrwDisplayTraX=disxtra
SrwDisplayAngZ=diszang
SrwDisplayTraZ=disztra

String xAng=TrjName+SrwSuffixX+SrwAngleWaveType
String xTra=TrjName+SrwSuffixX+SrwTrajWaveType
String zAng=TrjName+SrwSuffixZ+SrwAngleWaveType
String zTra=TrjName+SrwSuffixZ+SrwTrajWaveType

Duplicate/O $($MagName[0]) $xAng
SetScale d 0,0,"r", $xAng
Duplicate/O $xAng $xTra
SetScale d 0,0,"m", $xTra
Duplicate/O $($MagName[1]) $zAng
SetScale d 0,0,"r", $zAng
Duplicate/O $zAng $zTra
SetScale d 0,0,"m", $zTra

srTraj($ElecName, $MagName, $xAng, $xTra, $zAng, $zTra)
SrwTrjPrep(TrjName,ElecName,xTra,zTra)

if(disxang==2)
	Display $xAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorAngle
//else //OC180805
//	KillWaves/Z $xAng
endif
if(disxtra==2)
	Display $xTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorPos
endif
if(diszang==2)
	Display $zAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerAngle
//else
//	KillWaves/Z $zAng
endif
if(disztra==2)
	Display $zTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerPos
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Create Trajectory structure and compute trajectory in 3D magnetic field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwTrjCreate(TrjName, ElecName, MagName, ns, sbeg, send, disxang, disxtra, diszang, disztra)
string TrjName=SrwElecName+SrwMagName
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagGenTotName //SrwMagName+SrwFieldType
variable ns=srwUtiGetValN("ns", 101, "SrwTrjCreate")
variable sbeg=srwUtiGetValN("sbeg", 0, "SrwTrjCreate")
variable send=srwUtiGetValN("send", 0, "SrwTrjCreate")
variable disxang=SrwDisplayAngX, disxtra=SrwDisplayTraX, diszang=SrwDisplayAngZ, disztra=SrwDisplayTraZ
prompt TrjName,SrwPTrjName
prompt ElecName,SrwPElecName2,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+srwUtiGetValS("SrwMagField3dType", "_m3d", ""),";","")
prompt ns,"Number of Points vs Longitudinal Position"
prompt sbeg,"Initial Longitudinal Position [m]"
prompt send,"Final Longitudinal Position [m]"
prompt disxang,"Display x' vs s?",popup "No;Yes"
prompt disxtra,"Display x vs s?",popup "No;Yes"
prompt diszang,"Display z' vs s?",popup "No;Yes"
prompt disztra,"Display z vs s?",popup "No;Yes"
Silent 1						|	Computing the Trajectory  ...
PauseUpdate

if(strlen(TrjName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(TrjName)>28)
	Abort SrwPAlertTooLongName
endif

if(cmpstr(ElecName,"_none_")==0)
	Abort SrwPAlertNoElec
endif
if(cmpstr(MagName,"_none_")==0)
	Abort SrwPAlertNoMag
endif

SrwTrjName=TrjName
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
srwUtiSetValN("ns", ns, "SrwTrjCreate")
srwUtiSetValN("sbeg", sbeg, "SrwTrjCreate")
srwUtiSetValN("send", send, "SrwTrjCreate")
SrwDisplayAngX=disxang
SrwDisplayTraX=disxtra
SrwDisplayAngZ=diszang
SrwDisplayTraZ=disztra

String MagType=MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1]
//String/G SrwMagField3dType = "_m3d", SrwMagField3dWaveType = "_f3d"

String xAng=TrjName+SrwSuffixX+SrwAngleWaveType
String xTra=TrjName+SrwSuffixX+SrwTrajWaveType
String zAng=TrjName+SrwSuffixZ+SrwAngleWaveType
String zTra=TrjName+SrwSuffixZ+SrwTrajWaveType

String/G SrwSuffixY="Y"
String yAng=TrjName+SrwSuffixY+SrwAngleWaveType
String yTra=TrjName+SrwSuffixY+SrwTrajWaveType

if(cmpstr(MagType,SrwFieldType) == 0)
	Duplicate/O $($MagName[0]) $xAng
	Duplicate/O $($MagName[1]) $zAng
endif

Variable NpS, StartS, StepS
Variable aux_sRange
if(cmpstr(MagType,SrwMagField3dType) == 0)
	NpS = dimsize($($MagName[1]), 0)
	StartS = dimoffset($($MagName[1]), 0)
	StepS = dimdelta($($MagName[1]), 0)
	
	aux_sRange = StepS*(NpS - 1)
	
	if(ns > 1)
		NpS = ns
		StepS = aux_sRange/(NpS - 1)
	endif
	if(sbeg != send)
		StartS = sbeg
		StepS = 0
		if(NpS > 1)
			StepS = (send - sbeg)/(NpS - 1)
		endif
	endif
	
	Make/O/D/N=(NpS) $xAng, $zAng, $yAng
	SetScale/P x StartS, StepS, "m", $xAng, $zAng, $yAng
endif

SetScale d 0,0,"r", $xAng, $zAng
Duplicate/O $xAng $xTra
Duplicate/O $zAng $zTra
Duplicate/O $yAng $yTra
SetScale d 0,0,"m", $xTra, $zTra, $yTra

srTraj3d($ElecName, $MagName, $xAng, $xTra, $yAng, $yTra, $zAng, $zTra)
SrwTrj3dPrep(TrjName,ElecName,xTra,yTra,zTra)

if(disxang==2)
	Display $xAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorAngle
//else
//	KillWaves/Z $xAng
endif
if(disxtra==2)
	Display $xTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelHorPos
endif
if(diszang==2)
	Display $zAng
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerAngle
//else
//	KillWaves/Z $zAng
endif
if(disztra==2)
	Display $zTra
	Label bottom SrwPLabelLongPos
	Label left SrwPLabelVerPos
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Prepare Trajectory structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwTrjPrep(TrjName,ElecName,xTraName,zTraName)
String TrjName,ElecName,xTraName,zTraName

TrjName+=SrwTrjType

Make/T/O/N=6 $TrjName // Don't forget to increase when needed !!!
$TrjName[0]=xTraName
$TrjName[1]=zTraName
$TrjName[2]=""
$TrjName[3]=""
$TrjName[4]=num2str($ElecName[0]) 		// Electron Energy in GeV
$TrjName[5]=num2str($ElecName[1]) 		// Electron Current in A
end

//+++++++++++++++++++++++++++++++++++++++
//
// Prepare Trajectory structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwTrj3dPrep(TrjName,ElecName,xTraName,yTraName,zTraName)
String TrjName,ElecName,xTraName,yTraName,zTraName

TrjName+=SrwTrjType

Make/T/O/N=6 $TrjName // Don't forget to increase when needed !!!
$TrjName[0]=xTraName
$TrjName[1]=zTraName
$TrjName[2]=yTraName
$TrjName[3]=""
$TrjName[4]=num2str($ElecName[0]) 		// Electron Energy in GeV
$TrjName[5]=num2str($ElecName[1]) 		// Electron Current in A
end

//+++++++++++++++++++++++++++++++++++++++
//
// Auxiliary function, used to fing longitudinal position 
// of source point in a bending magnet
//
//+++++++++++++++++++++++++++++++++++++++
function srwTrjLongPosForAngle(nameElec, nameMagFld, angRad, strXorZ)
string nameMagFld, nameElec, strXorZ
variable angRad

//variable prevS0 = srwGetElecBeamLongPos(nameMagFld)
//variable prevX0 = srwGetElecBeamHorPos(nameMagFld)
//variable prevX0p = srwGetElecBeamHorAng(nameMagFld)
//variable prevZ0 = srwGetElecBeamVertPos(nameMagFld)
//variable prevZ0p = srwGetElecBeamVertAng(nameMagFld)

wave/T wMagFld = $nameMagFld
wave wMagFldBz = $(wMagFld[1])
variable sStartMagFld = dimoffset(wMagFldBz, 0)

//srwSetElecBeamLongPos(nameMagFld, 0)
//srwSetElecBeamHorPos(nameMagFld, 0)
//srwSetElecBeamHorAng(nameMagFld, 0)
//srwSetElecBeamVertPos(nameMagFld, 0)
//srwSetElecBeamVertAng(nameMagFld, 0)

string nameTrjCore = "AuxTrj"
string strToExe = "SrwTrjCreateTransvUnif(\"" + nameTrjCore + "\",\"" + nameElec + "\",\"" + nameMagFld + "\",1,1,1,1)"
execute/Z/Q strToExe

string nameTrj = nameTrjCore + "_trj"
string nameTrjXang = nameTrjCore + "X_ang", nameTrjXpos = nameTrjCore + "X_tra"
string nameTrjZang = nameTrjCore + "Z_ang", nameTrjZpos = nameTrjCore + "Z_tra"

string nameTrjAng = ""
if(cmpstr(strXorZ, "x") == 0)
	nameTrjAng = nameTrjXang
else
	nameTrjAng = nameTrjZang
endif
wave wTrjAng = $nameTrjAng

FindLevel/Q wTrjAng, angRad

variable resLongitOffset = V_LevelX
variable resWasFound = V_Flag
if(resWasFound == 1) //was not found
	print "WARNING : Longitudinal position was not found"
	resLongitOffset = 0
endif

//srwSetElecBeamLongPos(nameMagFld, prevS0)
//srwSetElecBeamHorPos(nameMagFld, prevX0)
//srwSetElecBeamHorAng(nameMagFld, prevX0p)
//srwSetElecBeamVertPos(nameMagFld, prevZ0)
//srwSetElecBeamVertAng(nameMagFld, prevZ0p)

wave wTrj = $nameTrj, wTrjXang = $nameTrjXang, wTrjXpos = $nameTrjXpos, wTrjZang = $nameTrjZang, wTrjZpos = $nameTrjZpos
killwaves/Z wTrj, wTrjXang, wTrjXpos, wTrjZang, wTrjZpos
return resLongitOffset
end

//+++++++++++++++++++++++++++++++++++++++
//
// Convert Arb. Transversely-Uniform Magnetic field to Periodic field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagArb2Per(MagPerName,MagName,RelPrec,MaxHarm,MaxPerLen)
string MagPerName = srwUtiTruncString(srwUtiGetValS("SrwMagName", "Mag", ""), 24) + "Per"
string MagName = srwUtiGetValS("SrwMagName", "Mag", "") + SrwFieldType
variable RelPrec = srwUtiGetValN("RelPrec", 0.05, "SrwMagArb2Per")
variable MaxHarm = srwUtiGetValN("MaxHarm",10, "SrwMagArb2Per")
variable MaxPerLen = srwUtiGetValN("MaxPerLen",1000, "SrwMagArb2Per")
prompt MagPerName, "Name of Periodic Magnetic Field structure to create"
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt RelPrec, "Relative accuracy threshold (bw 0 and 1)"
prompt MaxHarm, "Maximal number of magnetic field harmonics to create (from 1 to 20)"
prompt MaxPerLen, "Maximal magnetic field period length [mm]"
Silent 1						|	...
PauseUpdate

if(strlen(MagPerName)==0)
	Abort "Please specify name of the periodic magnetic feld structure to be created."
endif
if((cmpstr(MagName,"_none_")==0) %| (cmpstr(MagName,"")==0))
	Abort "Please set up a transversely uniform magnetic field structure to be converted to the periodic magnetic field structure."
endif
if((RelPrec <= 0) %| (RelPrec >= 1))
	Abort "Relative accuracy threshold should be between 0 and 1."
endif
if((MaxHarm <= 0) %| (MaxHarm > 20))
	Abort "Maximum number of field harmonics should be more than 0 and less than 21."
endif

srwUtiSetValS("SrwUndName", MagPerName, "")
srwUtiSetValS("SrwMagName", MagName[0,strlen(MagName)-strlen(SrwFieldType)-1], "")
srwUtiSetValN("RelPrec", RelPrec, "SrwMagArb2Per")
srwUtiSetValN("MaxHarm", MaxHarm, "SrwMagArb2Per")
srwUtiSetValN("MaxPerLen", MaxPerLen, "SrwMagArb2Per")

make/O/D/N=3 AuxPrecWave
AuxPrecWave[0] = RelPrec
AuxPrecWave[1] = MaxHarm
AuxPrecWave[2] = MaxPerLen*0.001

SrwUtiTriggerPrint(2)
SrwMagPerCreate2D(MagPerName,35,1,0,1.6,0,1,0,0)
SrwUtiTriggerPrint(1)
MagPerName += SrwUndType

srMagArb2Per($MagName, AuxPrecWave, $MagPerName)

KillWaves/Z AuxPrecWave
end

//+++++++++++++++++++++++++++++++++++++++
//
// Copies and scales Transversely - Uniform magnetic field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwMagDuplScale(nmMagRes, gapRes, nmMagGapMin, gapMin, nmBxPeakVsGap, nmBzPeakVsGap)
string nmMagRes=srwUtiGetValS("nmMagRes", "MagRes", "SrwMagDuplScale")
string nmMagGapMin=srwUtiGetValS("SrwMagGenTotName", "Mag", "")
string nmBxPeakVsGap=srwUtiGetValS("BxPeakVsGap", "BxPeakVsGap", "SrwMagDuplScale")
string nmBzPeakVsGap=srwUtiGetValS("BzPeakVsGap", "BzPeakVsGap", "SrwMagDuplScale")
variable gapRes=srwUtiGetValN("gapRes", 15.5, "SrwMagDuplScale")
variable gapMin=srwUtiGetValN("gapMin", 15.5, "SrwMagDuplScale")
prompt nmMagRes,"Name for Resulting Magnet Field"
prompt nmMagGapMin,"Magnetic Field at Min. Gap",popup Wavelist("*"+SrwFieldType,";","") + Wavelist("*"+SrwMagContainerType,";","")
prompt gapRes, "Resulting Gap [mm]"
prompt gapMin, "Minimal Gap [mm]"
prompt nmBxPeakVsGap, "Horizontal Field vs Gap"
prompt nmBzPeakVsGap, "Vertical Field vs Gap"
Silent 1						|	Computing Radiation  ...
PauseUpdate

srwUtiSetValS("nmMagRes", nmMagRes, "SrwMagDuplScale")
//srwUtiSetValS("nmMagGapMin", nmMagGapMin[0,strlen(nmMagGapMin)-strlen(SrwFieldType)-1], "")
srwUtiSetValS("SrwMagGenTotName", nmMagGapMin[0,strlen(nmMagGapMin)-strlen(SrwFieldType)-1], "")
srwUtiSetValS("nmBxPeakVsGap", nmBxPeakVsGap, "SrwMagDuplScale")
srwUtiSetValS("nmBzPeakVsGap", nmBzPeakVsGap, "SrwMagDuplScale")
srwUtiSetValN("gapRes", gapRes, "SrwMagDuplScale")
srwUtiSetValN("gapMin", gapMin, "SrwMagDuplScale")

variable BxPeakGapMin = $nmBxPeakVsGap(gapMin)
variable BxPeakGapRes = $nmBxPeakVsGap(gapRes)
variable BzPeakGapMin = $nmBzPeakVsGap(gapMin)
variable BzPeakGapRes = $nmBzPeakVsGap(gapRes)

SrwMagDupl(nmMagGapMin, nmMagRes)
string nmMagResBX = nmMagRes + "BX_fld"
string nmMagResBZ = nmMagRes + "BZ_fld"

$nmMagResBX *= BxPeakGapRes/BxPeakGapMin
$nmMagResBZ *= BzPeakGapRes/BzPeakGapMin
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculates section lengths and period lengths for an "Adaptive-Gap Undulator"
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwMagSectPerLenAGU(nmSectPerK, fundPhotEn, elecEn, nmUndGapVsLen, nmFundEnVsGapPer)
string nmSectPerK=srwUtiGetValS("nmSectPerK", "aguSectPerK", "SrwMagSectPerLenAGU")
variable fundPhotEn=srwUtiGetValN("fundPhotEn", 1000, "SrwMagSectPerLenAGU")
variable elecEn=srwUtiGetValN("elecEn", 3, "SrwMagSectPerLenAGU")
variable minPer=srwUtiGetValN("minPer", 20, "SrwMagSectPerLenAGU")
string nmUndGapVsLen=srwUtiGetValS("nmUndGapVsLen", "aguGapVsLen", "SrwMagSectPerLenAGU")
string nmFundEnVsGapPer=srwUtiGetValS("nmFundEnVsGapPer", "aguFundEnVsGapPer", "SrwMagSectPerLenAGU")
prompt nmSectPerK,"Name for Resulting Period and K wave"
prompt fundPhotEn,"Fundamental Photon Energy [eV]"
prompt elecEn,"Electron Energy [GeV]"
prompt nmUndGapVsLen,"Name of Undulator Gap vs Length wave", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt nmFundEnVsGapPer,"Name of Fund. Energy vs Gap and Period wave", popup Wavelist("*",";","TEXT:0,DIMS:2")
Silent 1						|	Computing Radiation  ...
PauseUpdate

srwUtiSetValS("nmSectPerK", nmSectPerK, "SrwMagSectPerLenAGU")
srwUtiSetValN("fundPhotEn", fundPhotEn, "SrwMagSectPerLenAGU")
srwUtiSetValN("elecEn", elecEn, "SrwMagSectPerLenAGU")
srwUtiSetValS("nmUndGapVsLen", nmUndGapVsLen, "SrwMagSectPerLenAGU")
srwUtiSetValS("nmFundEnVsGapPer", nmFundEnVsGapPer, "SrwMagSectPerLenAGU")

variable numSect = dimsize($nmUndGapVsLen, 0)
variable sStart = dimoffset($nmUndGapVsLen, 0)
variable sStep = dimdelta($nmUndGapVsLen, 0)
make/O/N=(numSect, 2) $nmSectPerK
SetScale/P x sStart, sStep,"m", $nmSectPerK
$nmSectPerK = 0

variable numPer = dimsize($nmFundEnVsGapPer, 1)
variable perStart = dimoffset($nmFundEnVsGapPer, 1)
variable perStep = dimdelta($nmUndGapVsLen, 1)
variable perEnd = perStart + (numPer - 1)*perStep
make/O/N=(numPer) wAuxFundEnVsPer
SetScale/P x perStart, perStep,"mm", wAuxFundEnVsPer

variable numGaps = dimsize($nmFundEnVsGapPer, 0)
variable gapStart = dimoffset($nmFundEnVsGapPer, 0)
variable gapStep = dimdelta($nmFundEnVsGapPer, 0)
make/O/N=(numGaps) wAuxPhotEnVsGap
SetScale/P x gapStart, gapStep,"mm", wAuxPhotEnVsGap

variable nCoefInterpVsGap = 5, nCoefInterpVsPer = numPer - 2 //to tune
if(nCoefInterpVsPer <= 0)
	nCoefInterpVsPer = 1
endif
variable numPtVsPerAfterInterp = 200
make/O/N=(nCoefInterpVsGap, numPer) wCoefInterpVsGap
make/O/N=(nCoefInterpVsGap) wAuxInterpCoef

make/O/N=(numPtVsPerAfterInterp) wAuxFundEnVsPerAfterInterp
SetScale/I x perStart-0.5, perEnd+0.5,"mm", wAuxFundEnVsPerAfterInterp

variable iPer = 0, per
do
	wAuxPhotEnVsGap = $nmFundEnVsGapPer[p][iPer]
	CurveFit/Q/W=0 poly nCoefInterpVsGap, wAuxPhotEnVsGap
	wCoefInterpVsGap[][iPer] = W_coef[p]
	
	iPer += 1
while(iPer < numPer)

variable iSect = 0, gap
do
	gap = $nmUndGapVsLen[iSect]
	iPer = 0
	do
		wAuxInterpCoef = wCoefInterpVsGap[p][iPer]
		wAuxFundEnVsPer[iPer] = poly(wAuxInterpCoef, gap)
		iPer += 1
	while(iPer < numPer)
	
	CurveFit/Q/W=0 poly nCoefInterpVsPer, wAuxFundEnVsPer
	wAuxFundEnVsPerAfterInterp = poly(W_coef, x)
	
	FindLevel/Q wAuxFundEnVsPerAfterInterp, fundPhotEn
	if(V_flag == 0)
		per = V_LevelX
		$nmSectPerK[iSect][0] = per //period
		$nmSectPerK[iSect][1] = sqrt((9496.3421866853*elecEn*elecEn/per/fundPhotEn - 1)*2)
	endif
	iSect += 1
while(iSect < numSect)

edit/K=0 $nmSectPerK

killwaves/Z wAuxFundEnVsPerAfterInterp, wAuxInterpCoef, wAuxPhotEnVsGap, wAuxFundEnVsPer
end
