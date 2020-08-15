
//+++++++++++++++++++++++++++++++++++++++
//
//Create Observation General
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpCreate(ObsName, Dist)
String ObsName=SrwSmpName
Variable Dist=SrwSmpDist
prompt ObsName,SrwPSmpName
prompt Dist,SrwPSmpDist
Silent 1						|	Creating the Observation Structure  ...
PauseUpdate

SrwSmpName=ObsName
SrwSmpDist=Dist

ObsName+=SrwSmpType
SrwSmpGenTotName=ObsName

//Make/N=14/D/O $ObsName
//Make/N=18/D/O $ObsName
Make/N=27/D/O $ObsName

$ObsName[0]=0
$ObsName[4]=Dist

// Energy scan
$ObsName[5]=SrwSmpEdep*1000
$ObsName[6]=SrwSmpEfin*1000
$ObsName[7]=SrwSmpEnpts
// Hor. scan
$ObsName[8]=SrwSmpXmid*0.001
$ObsName[9]=SrwSmpXmagn*0.001
$ObsName[10]=SrwSmpXnpts
// Vert. scan
$ObsName[11]=SrwSmpZmid*0.001
$ObsName[12]=SrwSmpZmagn*0.001
$ObsName[13]=SrwSmpZnpts
// Time scan takes positions [15], [16], [17]; [14] defines Time/Frequency representation

End  

//+++++++++++++++++++++++++++++++++++++++
//
//Sets parameters of Electron Beam structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpParam(Name, IndParam, ValToSet)
string Name
variable IndParam, ValToSet
SVAR SrwSmpType
srwUtiSetValS("SrwSmpName", Name[0,strlen(Name)-strlen(SrwSmpType)-1], "")
wave wSmp = $Name
variable sizeSmp = DimSize(wSmp, 0)
if(sizeSmp < 14)
	abort "This is not a radiation sampling structure."
endif
if(IndParam > 16)
	abort "Parameter index value is too large."
endif
if(IndParam >= sizeSmp)
	redimension/N=(IndParam + 1) wSmp
endif
wSmp[IndParam] = ValToSet
end

//+++++++++++++++++++++++++++++++++++++++
//
//Return parameters of Radiation Sampling structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpParam(Name, IndParam)
string Name
variable IndParam
SVAR SrwSmpType
srwUtiSetValS("SrwSmpName", Name[0,strlen(Name)-strlen(SrwSmpType)-1], "")
wave wSmp = $Name
if(DimSize(wSmp, 0) < 14)
	abort "This is not a radiation sampling structure."
endif
if(IndParam > 16)
	abort "Parameter index value is too large."
endif
return wSmp[IndParam]
end

//+++++++++++++++++++++++++++++++++++++++
//Returns longitudinal position of the observation plane
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpLongPos(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 4)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets longitudinal position of the observation plane
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpLongPos(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 4, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns number of points vs horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpHorPosNp(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 10)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets number of points vs horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpHorPosNp(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 10, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns initial horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpHorPosStart(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 8)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets initial horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpHorPosStart(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 8, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns final horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpHorPosEnd(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 9)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets final horizontal position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpHorPosEnd(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 9, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns number of points vs vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpVertPosNp(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 13)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets number of points vs vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpVertPosNp(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 13, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns initial vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpVertPosStart(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 11)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets initial vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpVertPosStart(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 11, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns final vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpVertPosEnd(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 12)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets final vertical position
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpVertPosEnd(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 12, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns number of points vs photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpPhotEnNp(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 7)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets number of points vs photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpPhotEnNp(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 7, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns initial photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpPhotEnStart(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 5)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets initial photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpPhotEnStart(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 5, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns final photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpPhotEnEnd(ObsName)
string ObsName
return srwGetSmpParam(ObsName, 6)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets final photon energy
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpPhotEnEnd(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 6, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//Returns number of points vs time
//+++++++++++++++++++++++++++++++++++++++
function srwGetSmpTimeNp(ObsName)
string ObsName
variable totNumData = dimsize($ObsName, 0)
variable indTimeNp = 14
if(totNumData <= indTimeNp)
	return 0
endif
return srwGetSmpParam(ObsName, indTimeNp)
end
//+++++++++++++++++++++++++++++++++++++++
//Sets number of points vs time
//+++++++++++++++++++++++++++++++++++++++
function srwSetSmpTimeNp(ObsName, ValToSet)
string ObsName
variable ValToSet
srwSetSmpParam(ObsName, 14, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Modify Horizontal,Vertical,Energy scan values
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpScanXZE(ObsName, Xmid, Xmagn, Xnpts, Zmid, Zmagn, Znpts,dep,fin,npts)
String ObsName=SrwSmpName+SrwSmpType
Variable Xmid=SrwSmpXmid
Variable Xmagn=SrwSmpXmagn
Variable Xnpts=SrwSmpXnpts
Variable Zmid=SrwSmpZmid
Variable Zmagn=SrwSmpZmagn
Variable Znpts=SrwSmpZnpts
Variable dep=SrwSmpEdep
Variable fin=SrwSmpEfin
Variable npts=SrwSmpEnpts
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType ,";", "")
prompt Xmid,SrwPSmpXmid
prompt Xmagn,SrwPSmpXmagn
prompt Xnpts,SrwPSmpXnpts
prompt Zmid,SrwPSmpZmid
prompt Zmagn,SrwPSmpZmagn
prompt Znpts,SrwPSmpZnpts
prompt dep,SrwPSmpEdep
prompt fin,SrwPSmpEfin
prompt npts,SrwPSmpEnpts
Silent 1						|	Updating the Observation Structure  ...
PauseUpdate

// Validation of parameters
if(Xmagn < 0.)
	Abort SrwPAlertSmpRange
endif
if(Zmagn < 0.)
	Abort SrwPAlertSmpRange
endif
if(dep < 0.)
	Abort SrwPAlertSmpEnergy
endif
if(fin < 0.)
	Abort SrwPAlertSmpEnergy
endif
if((abs(round(npts) - npts) > 1.E-08) %| (npts <= 0))
	Abort SrwPAlertSmpEnergyNpts
endif
if((abs(round(Xnpts) - Xnpts) > 1.E-08))
	Abort SrwPAlertSmpNpts
endif
if((abs(round(Znpts) - Znpts) > 1.E-08))
	Abort SrwPAlertSmpNpts
endif

if(Xnpts<1)
	Xnpts=1
endif
if(Znpts<1)
	Znpts=1
endif
if(npts<1)
	npts=1
endif

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwSmpXmid=Xmid
SrwSmpXmagn=Xmagn
SrwSmpXnpts=Xnpts
SrwSmpZmid=Zmid
SrwSmpZmagn=Zmagn
SrwSmpZnpts=Znpts
SrwSmpEdep=dep
SrwSmpEfin=fin
SrwSmpEnpts=npts
SrwSmpGenTotName=ObsName

if(cmpstr(ObsName,"_none_")==0)
abort  SrwNoObsName
endif

if(dimsize($ObsName, 0) < 15)
	redimension/N=15 $ObsName
endif

//Setting required Electric Field representation to FREQUENCY:
$ObsName[14] = 0

$ObsName[5]=dep*1000
$ObsName[6]=fin*1000
$ObsName[7]=npts

//if(Xnpts != 1) //AuxDebugTest
	$ObsName[8]=(Xmid - 0.5*Xmagn)*0.001
	$ObsName[9]=(Xmid + 0.5*Xmagn)*0.001
//else
//	$ObsName[8]=Xmid*0.001
//	$ObsName[9]=Xmid*0.001
//endif
$ObsName[10]=Xnpts

//if(Znpts != 1) //AuxDebugTest
	$ObsName[11]=(Zmid - 0.5*Zmagn)*0.001
	$ObsName[12]=(Zmid + 0.5*Zmagn)*0.001
//else
//	$ObsName[11]=Zmid*0.001
//	$ObsName[12]=Zmid*0.001
//endif
$ObsName[13]=Znpts
end  

//+++++++++++++++++++++++++++++++++++++++
//
//Modify/Set Time scan values
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpScanTime(ObsName, tStart, tEnd, tNp)
string ObsName=SrwSmpName+SrwSmpType
variable tStart=srwUtiGetValN("SrwSmpTstart", 0, "")
variable tEnd=srwUtiGetValN("SrwSmpTend", 0, "")
variable tNp=srwUtiGetValN("SrwSmpTnpts", 1, "")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType ,";", "")
prompt tStart, "Initial Time [fs]"
prompt tEnd, "Final Time [fs]"
prompt tNp, "Number of Time Points"
Silent 1						|	Updating the Observation Structure  ...
PauseUpdate

// Validation of parameters
if(tNp <= 0)
	abort "Number of Time Moments should be positive"
endif
if(cmpstr(ObsName,"_none_") == 0)
	abort  SrwNoObsName
endif
if(tNp < 1)
	tNp = 1
endif

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
srwUtiSetValN("SrwSmpTstart", tStart, "")
srwUtiSetValN("SrwSmpTend", tEnd, "")
srwUtiSetValN("SrwSmpTnpts", tNp, "")
SrwSmpGenTotName=ObsName

if(dimsize($ObsName, 0) < 18)
	//make/O/N=18/D $ObsName
	redimension/N=18 $ObsName
endif

//Setting required Electric Field representation to TIME:
$ObsName[14] = 1

// Time scan takes positions [15], [16], [17]
// Time untits are [s] for internal storage
$ObsName[15] = tStart*1e-15
$ObsName[16] = tEnd*1e-15
$ObsName[17] = tNp
end 

//+++++++++++++++++++++++++++++++++++++++
//
//Create Observation for Power Density (without energy)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpPowCreate(ObsName, Dist)
String ObsName=SrwSmpName
Variable Dist=SrwSmpDist
prompt ObsName,SrwPSmpName
prompt Dist,SrwPSmpDist
Silent 1						|	Creating the Observation Structure ...
PauseUpdate

SrwSmpName=ObsName
SrwSmpDist=Dist

ObsName+=SrwSmpPowType
SrwSmpGenTotName=ObsName

Make/N=14/D/O $ObsName
$ObsName[0]=0
$ObsName[4]=Dist

// Energy scan
$ObsName[5]=0.12345 // Arbitrary number
$ObsName[6]=0.12345 // Arbitrary number
$ObsName[7]=12345 // Arbitrary number
// Hor scan
$ObsName[8]=SrwSmpXmid*0.001
$ObsName[9]=SrwSmpXmagn*0.001
$ObsName[10]=SrwSmpXnpts
// Ver scan
$ObsName[11]=SrwSmpZmid*0.001
$ObsName[12]=SrwSmpZmagn*0.001
$ObsName[13]=SrwSmpZnpts

End  

//+++++++++++++++++++++++++++++++++++++++
//
//Modify Horizontal, Vertical scan values for Power Density
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpPowScanXZ(ObsName, Xmid, Xmagn, Xnpts, Zmid, Zmagn, Znpts)
String ObsName=SrwSmpName+SrwSmpPowType
Variable Xmid=SrwSmpXmid
Variable Xmagn=SrwSmpXmagn
Variable Xnpts=SrwSmpXnpts
Variable Zmid=SrwSmpZmid
Variable Zmagn=SrwSmpZmagn
Variable Znpts=SrwSmpZnpts
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType ,";", "")
prompt Xmid,SrwPSmpXmid
prompt Xmagn,SrwPSmpXmagn
prompt Xnpts,SrwPSmpXnpts
prompt Zmid,SrwPSmpZmid
prompt Zmagn,SrwPSmpZmagn
prompt Znpts,SrwPSmpZnpts
Silent 1						|	Updating the Observation Structure  ...
PauseUpdate

// Validation of parameters
if(Xmagn < 0.)
	Abort SrwPAlertSmpRange
endif
if(Zmagn < 0.)
	Abort SrwPAlertSmpRange
endif
if((abs(round(Xnpts) - Xnpts) > 1.E-08))
	Abort SrwPAlertSmpNpts
endif
if((abs(round(Znpts) - Znpts) > 1.E-08))
	Abort SrwPAlertSmpNpts
endif

if(Xnpts<1)
	Xnpts=1
endif
if(Znpts<1)
	Znpts=1
endif

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1]
SrwSmpXmid=Xmid
SrwSmpXmagn=Xmagn
SrwSmpXnpts=Xnpts
SrwSmpZmid=Zmid
SrwSmpZmagn=Zmagn
SrwSmpZnpts=Znpts
SrwSmpGenTotName=ObsName

if(cmpstr(ObsName,"_none_")==0)
abort  SrwNoObsName
endif

$ObsName[8]=(Xmid - 0.5*Xmagn)*0.001
$ObsName[9]=(Xmid + 0.5*Xmagn)*0.001
$ObsName[10]=Xnpts

$ObsName[11]=(Zmid - 0.5*Zmagn)*0.001
$ObsName[12]=(Zmid + 0.5*Zmagn)*0.001
$ObsName[13]=Znpts

end  

//+++++++++++++++++++++++++++++++++++++++
//
//Extends Orservation ranges to allow for accurate treatment 
//of e-beam emittance by convolution
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpExtForElecBeam(nmObsNew, nmObsOld, nmElec, numExtraSigProj)
string nmObsNew=SrwSmpName
string nmObsOld=srwUtiGetValS("SrwSmpName", "Obs", "")+SrwSmpType
string nmElec=srwUtiGetValS("SrwElecName", "Elec", "")+SrwElecType
variable numExtraSigProj=srwUtiGetValN("numExtraSigProj", 2, "SrwSmpExtForElecBeam")
prompt nmObsNew,"Name for new Radiation Sampling structure"
prompt nmObsOld,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt nmElec,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt numExtraSigProj,"Number of Projected RMS E-Beam Sizes to add"
Silent 1						|	Creating Observation Structure  ...
PauseUpdate

SrwElecName=nmElec[0,strlen(nmElec)-strlen(SrwElecType)-1]
SrwSmpName=nmObsNew

nmObsNew += SrwSmpType
duplicate/O $nmObsOld $nmObsNew

variable xStartOld = srwGetSmpHorPosStart(nmObsOld)
variable xEndOld = srwGetSmpHorPosEnd(nmObsOld)
variable xNpOld = srwGetSmpHorPosNp(nmObsOld)

variable zStartOld = srwGetSmpVertPosStart(nmObsOld)
variable zEndOld = srwGetSmpVertPosEnd(nmObsOld)
variable zNpOld = srwGetSmpVertPosNp(nmObsOld)

variable rObs = srwGetSmpLongPos(nmObsOld)

variable xSigProj = srwGetElecBeamHorSizeProjRMS(nmElec, rObs)
variable zSigProj = srwGetElecBeamVertSizeProjRMS(nmElec, rObs)

variable xCenOld = 0.5*(xStartOld + xEndOld)
variable xStep = 0
if(xNpOld > 1)
	xStep = (xEndOld - xStartOld)/(xNpOld - 1)
endif
variable xNumStepsExtraHalf = trunc(numExtraSigProj*xSigProj/xStep) + 1
variable xRangeExtraHalf = xNumStepsExtraHalf*xStep
variable xRangeNew = (xEndOld - xStartOld) + 2*xRangeExtraHalf
variable xNpNew = xNpOld + 2*xNumStepsExtraHalf

variable zCenOld = 0.5*(zStartOld + zEndOld)
variable zStep = 0
if(zNpOld > 1)
	zStep = (zEndOld - zStartOld)/(zNpOld - 1)
endif
variable zNumStepsExtraHalf = trunc(numExtraSigProj*zSigProj/zStep) + 1
variable zRangeExtraHalf = zNumStepsExtraHalf*zStep
variable zRangeNew = (zEndOld - zStartOld) + 2*zRangeExtraHalf
variable zNpNew = zNpOld + 2*zNumStepsExtraHalf

variable eStartOld = srwGetSmpPhotEnStart(nmObsOld)
variable eEndOld = srwGetSmpPhotEnEnd(nmObsOld)
variable eNpOld = srwGetSmpPhotEnNp(nmObsOld)
variable dEperE = srwGetElecBeamRelEnSprRMS(nmElec)
variable photEnExtCoef = 2*dEperE*numExtraSigProj
variable eNumStepsExtraLeft = 0, eNumStepsExtraRight = 0
variable eStep = 0
if(eNpOld > 1)
	eStep = (eEndOld - eStartOld)/(eNpOld - 1)
	eNumStepsExtraLeft = trunc(eStartOld*photEnExtCoef/eStep) + 1
	eNumStepsExtraRight = trunc(eEndOld*photEnExtCoef/eStep) + 1
endif
variable eStartNew = eStartOld - eNumStepsExtraLeft*eStep
variable eEndNew = eEndOld + eNumStepsExtraRight*eStep
variable eNpNew = eNpOld + eNumStepsExtraLeft + eNumStepsExtraRight

SrwSmpScanXZE(nmObsNew, xCenOld*1000, xRangeNew*1000, xNpNew, zCenOld*1000, zRangeNew*1000, zNpNew, eStartNew*0.001, eEndNew*0.001, eNpNew)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Observation Panel
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplDialog(SamplType)
string SamplType

string/G SrwSmpType
if(strlen(SamplType) <= 1)
	SamplType = SrwSmpType
endif

SrwRadSamplCreateBufVars(SamplType)
DoWindow/K SrwRadSamplPanel

variable ShowEnergy = 1
if(cmpstr(SamplType,SrwSmpPowType)==0)
	ShowEnergy = 0
endif

SrwRadSamplPanel(ShowEnergy)

 End
 
//+++++++++++++++++++++++++++++++++++++++
Proc SrwExecuteAnyMacroAfterRadSmp()

Variable LocStartMacrosAfterRadSmp = SrwStartMacrosAfterRadSmp
Variable LocStartMacrosAfterRadSmp2 = SrwStartMacrosAfterRadSmp2

SrwStartMacrosAfterRadSmp = 0
SrwStartMacrosAfterRadSmp2 = 0

if(LocStartMacrosAfterRadSmp == 1)
	SrwWfrCreate_()
endif
if(LocStartMacrosAfterRadSmp == 2)
	SrwWfrCreate()
endif
if(LocStartMacrosAfterRadSmp == 3)
	SrwPerStoCreate()
endif
if(LocStartMacrosAfterRadSmp == 4)
	SrwPowCreate()
endif
if(LocStartMacrosAfterRadSmp == 5)
	SrwStoWigCreate()
endif

if(abs(LocStartMacrosAfterRadSmp2) == 1)
	SrwWfr2Int_()
endif
if(abs(LocStartMacrosAfterRadSmp2) == 2)
	SrwWfr2Int();
endif
if(abs(LocStartMacrosAfterRadSmp2) == 3)
	SrwSto2Int()
endif
if(abs(LocStartMacrosAfterRadSmp2) == 4)
	SrwPow2Int()
endif
if(abs(LocStartMacrosAfterRadSmp2) == 5)
	SrwSto2Int()
endif

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwSuppressMacrosAfterRadSmp()

SrwStartMacrosAfterRadSmp = 0;
SrwStartMacrosAfterRadSmp2 = 0;

End

//+++++++++++++++++++++++++++++++++++++++
Window SrwRadSamplPanel(ShowEnergy) : Panel
Variable ShowEnergy;

PauseUpdate; Silent 1		// building window...
NewPanel /W=(416,65,725,478) as SrwPRadSamplTitle;
SetDrawLayer UserBack;

SrwRadSamplDrawAllContr(ShowEnergy);

EndMacro

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplDrawAllContr(ShowEnergy)
Variable ShowEnergy;

String/G SrwSmpNameBuf, SrwSmpTypeBuf

String ObsName=SrwSmpNameBuf+SrwSmpTypeBuf;
String AllObsWaves=Wavelist("*"+SrwSmpTypeBuf,";","");

Variable ItemNo=sRWaveListItemNo(AllObsWaves, ";", ObsName);
//if(ItemNo == 0)
//	ItemNo = 1;
//endif

Variable VertOffset = -105;

SetDrawEnv fname= "default",fsize= 12;
SetDrawEnv fstyle= 1;
DrawText 20,140+VertOffset,SrwPRadSamplExist;
PopupMenu popup0RadSampl,pos={20,145+VertOffset},size={255,21};
PopupMenu popup0RadSampl,value= # "Wavelist(\"*\"+SrwSmpTypeBuf,\";\",\"\")",mode=ItemNo,proc=SrwRadSamplSelectPopupProc;

// Name and Long. position
DrawRect 13,180+VertOffset,296,245+VertOffset;

SetDrawEnv textyjust= 2;
SetDrawEnv fstyle= 1;
DrawText 25,187+VertOffset,SrwPRadSamplName;

SetVariable setvar0RadSampl,pos={65,185+VertOffset},size={226,17},title=" ",fSize=14;
SetVariable setvar0RadSampl,limits={-Inf,Inf,1},value= SrwSmpNameBuf;

SetDrawEnv fstyle= 1;
DrawText 25,238+VertOffset,"Longitudinal";
SetDrawEnv textyjust= 2;
DrawText 120,223+VertOffset,SrwPRadSamplLongPos;

SetVariable setvar1RadSampl,pos={205,220+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar1RadSampl,limits={-Inf,Inf,1},value= SrwSmpDistBuf

// Horizontal Sampling
DrawRect 13,250+VertOffset,296,320+VertOffset;

SetDrawEnv textyjust= 2;
DrawText 120,258+VertOffset,SrwPRadSamplHorPosMid;
SetDrawEnv textyjust= 2;
DrawText 120,278+VertOffset,SrwPRadSamplHorPosRng;
SetDrawEnv textyjust= 2;
DrawText 120,298+VertOffset,SrwPRadSamplHorPtNum;
SetDrawEnv fstyle= 1;
DrawText 25,272+VertOffset,"Horizontal";

SetVariable setvar2RadSampl,pos={205,255+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar2RadSampl,limits={-Inf,Inf,1},value= SrwSmpXmidBuf
SetVariable setvar3RadSampl,pos={205,275+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar3RadSampl,limits={-Inf,Inf,1},value= SrwSmpXmagnBuf
SetVariable setvar4RadSampl,pos={205,295+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar4RadSampl,limits={-Inf,Inf,1},value= SrwSmpXnptsBuf

// Vertical Sampling
DrawRect 13,325+VertOffset,296,395+VertOffset;

SetDrawEnv textyjust= 2;
DrawText 120,333+VertOffset,SrwPRadSamplVerPosMid;
SetDrawEnv textyjust= 2;
DrawText 120,353+VertOffset,SrwPRadSamplVerPosRng;
SetDrawEnv textyjust= 2;
DrawText 120,373+VertOffset,SrwPRadSamplVerPtNum;
SetDrawEnv fstyle= 1;
DrawText 25,347+VertOffset,"Vertical";

SetVariable setvar5RadSampl,pos={205,330+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar5RadSampl,limits={-Inf,Inf,1},value= SrwSmpZmidBuf
SetVariable setvar6RadSampl,pos={205,350+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar6RadSampl,limits={-Inf,Inf,1},value= SrwSmpZmagnBuf
SetVariable setvar7RadSampl,pos={205,370+VertOffset},size={85,17},title=" ",fSize=14;
SetVariable setvar7RadSampl,limits={-Inf,Inf,1},value= SrwSmpZnptsBuf

if(ShowEnergy != 0) // Photon Energy Sampling
	DrawRect 13,400+VertOffset,296,470+VertOffset;

	SetDrawEnv textyjust= 2;
	DrawText 120,408+VertOffset,SrwPRadSamplEnInit;
	SetDrawEnv textyjust= 2;
	DrawText 120,428+VertOffset,SrwPRadSamplEnFin;
	SetDrawEnv textyjust= 2;
	DrawText 120,448+VertOffset,SrwPRadSamplEnPtNum;
	SetDrawEnv fstyle= 1;
	DrawText 25,422+VertOffset,"Energy";

	SetVariable setvar8RadSampl,pos={205,405+VertOffset},size={85,17},title=" ",fSize=14;
	SetVariable setvar8RadSampl,limits={-Inf,Inf,1},value= SrwSmpEdepBuf
	SetVariable setvar9RadSampl,pos={205,425+VertOffset},size={85,17},title=" ",fSize=14;
	SetVariable setvar9RadSampl,limits={-Inf,Inf,1},value= SrwSmpEfinBuf
	SetVariable setvar10RadSampl,pos={205,445+VertOffset},size={85,17},title=" ",fSize=14;
	SetVariable setvar10RadSampl,limits={-Inf,Inf,1},value= SrwSmpEnptsBuf
endif

// OK-Cancel-Help
Button button5RadSampl,pos={120,485+VertOffset},size={70,20},proc=SrwRadSamplPanelOKButtonProc,title=SrwPPanelOKButton;
Button button2RadSampl,pos={30,485+VertOffset},size={70,20},proc=SrwRadSamplCancelButtonProc,title=SrwPPanelCancelButton;
Button button6RadSampl,pos={210,485+VertOffset},size={70,20},proc=SrwRadSamplPanelHelpButtonProc,title=SrwPPanelHelpButton;

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplCancelButtonProc(ctrlName) : ButtonControl
String ctrlName;

DoWindow/K SrwRadSamplPanel;
SrwRadSamplKillBufVars();

SrwSuppressMacrosAfterRadSmp();

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplPanelOKButtonProc(ctrlName) : ButtonControl
String ctrlName;
	
DoWindow/K SrwRadSamplPanel;

String ComLineStr;

String TotObsName = SrwSmpNameBuf + SrwSmpTypeBuf;
if(cmpstr(SrwSmpTypeBuf,SrwSmpType)==0)
	SrwSmpCreate(SrwSmpNameBuf, SrwSmpDistBuf); 
	SrwSmpScanXZE(TotObsName, SrwSmpXmidBuf, SrwSmpXmagnBuf, SrwSmpXnptsBuf, SrwSmpZmidBuf, SrwSmpZmagnBuf, SrwSmpZnptsBuf, SrwSmpEdepBuf, SrwSmpEfinBuf, SrwSmpEnptsBuf);
	sprintf ComLineStr, "SrwSmpCreate(\"%s\",%g);SrwSmpScanXZE(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g,%g)", SrwSmpNameBuf, SrwSmpDistBuf, TotObsName, SrwSmpXmidBuf, SrwSmpXmagnBuf, SrwSmpXnptsBuf, SrwSmpZmidBuf, SrwSmpZmagnBuf, SrwSmpZnptsBuf, SrwSmpEdepBuf, SrwSmpEfinBuf, SrwSmpEnptsBuf;
endif
if(cmpstr(SrwSmpTypeBuf,SrwSmpPowType)==0)
	SrwSmpPowCreate(SrwSmpNameBuf, SrwSmpDistBuf); 
	SrwSmpPowScanXZ(TotObsName, SrwSmpXmidBuf, SrwSmpXmagnBuf, SrwSmpXnptsBuf, SrwSmpZmidBuf, SrwSmpZmagnBuf, SrwSmpZnptsBuf);
	sprintf ComLineStr, "SrwSmpPowCreate(\"%s\",%g);SrwSmpPowScanXZ(\"%s\",%g,%g,%g,%g,%g,%g)", SrwSmpNameBuf, SrwSmpDistBuf, TotObsName, SrwSmpXmidBuf, SrwSmpXmagnBuf, SrwSmpXnptsBuf, SrwSmpZmidBuf, SrwSmpZmagnBuf, SrwSmpZnptsBuf;
endif

Print ComLineStr;

SrwRadSamplKillBufVars();
SrwExecuteAnyMacroAfterRadSmp();

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplSelectPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName;
Variable popNum;		// which item is currently selected (1-based)
String popStr;		// contents of current popup item as string
	
SrwSmpNameBuf = popStr[0,strlen(popStr) - strlen(SrwSmpTypeBuf) - 1];
String ObsName = SrwSmpNameBuf+SrwSmpTypeBuf;
	
SrwSmpDistBuf = $ObsName[4];
// Energy scan
SrwSmpEdepBuf = ($ObsName[5])*0.001;
SrwSmpEfinBuf = ($ObsName[6])*0.001;
SrwSmpEnptsBuf = $ObsName[7];
// Hor scan
SrwSmpXmidBuf = 0.5*($ObsName[8] + $ObsName[9])*1000.;
SrwSmpXmagnBuf = ($ObsName[9] - $ObsName[8])*1000.;
SrwSmpXnptsBuf = $ObsName[10];
// Ver scan
SrwSmpZmidBuf = 0.5*($ObsName[11] + $ObsName[12])*1000.;
SrwSmpZmagnBuf = ($ObsName[12] - $ObsName[11])*1000.;
SrwSmpZnptsBuf = $ObsName[13];

ControlUpdate/A/W=SrwRadSamplPanel;

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplPanelHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwRadSamplDialog     ")
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplPanelKillContr()

KillControl popup0RadSampl;
KillControl setvar0RadSampl;
KillControl setvar1RadSampl;
KillControl setvar2RadSampl;
KillControl setvar3RadSampl;
KillControl setvar4RadSampl;
KillControl setvar5RadSampl;
KillControl setvar6RadSampl;
KillControl setvar7RadSampl;
KillControl setvar8RadSampl;
KillControl setvar9RadSampl;
KillControl setvar10RadSampl;

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplCreateBufVars(SamplType)
String SamplType;

String/G SrwSmpNameBuf = SrwSmpName;
//String/G SrwSmpTypeBuf = SrwSmpType;
String/G SrwSmpTypeBuf = SamplType;

Variable/G SrwSmpDistBuf = SrwSmpDist;
Variable/G SrwSmpXmidBuf = SrwSmpXmid;
Variable/G SrwSmpXmagnBuf=SrwSmpXmagn;
Variable/G SrwSmpXnptsBuf=SrwSmpXnpts;
Variable/G SrwSmpZmidBuf=SrwSmpZmid;
Variable/G SrwSmpZmagnBuf=SrwSmpZmagn;
Variable/G SrwSmpZnptsBuf=SrwSmpZnpts;
Variable/G SrwSmpEdepBuf=SrwSmpEdep;
Variable/G SrwSmpEfinBuf=SrwSmpEfin;
Variable/G SrwSmpEnptsBuf=SrwSmpEnpts;
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwRadSamplKillBufVars()

KillStrings/Z SrwSmpNameBuf, SrwSmpTypeBuf;
KillVariables/Z SrwSmpDistBuf, SrwSmpXmidBuf, SrwSmpXmagnBuf, SrwSmpXnptsBuf, SrwSmpZmidBuf, SrwSmpZmagnBuf, SrwSmpZnptsBuf, SrwSmpEdepBuf, SrwSmpEfinBuf, SrwSmpEnptsBuf;

End

//+++++++++++++++++++++++++++++++++++++++
//
//Define Orientation of the Observation Plane
//Test version
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpPlaneOrient(ObsName, nx, ny, nz)
string ObsName=SrwSmpGenTotName
variable nx=srwUtiGetValN("nx", 0, "SrwSmpPlaneOrient")
variable ny=srwUtiGetValN("ny", -1, "SrwSmpPlaneOrient")
variable nz=srwUtiGetValN("nz", 0, "SrwSmpPlaneOrient")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt nx, "Horizontal Coordinate of the Vector Normal to Observation Plane"
prompt ny, "Longitudinal Coordinate of the Vector Normal to Observation Plane"
prompt nz, "Vertical Coordinate of the Vector Normal to Observation Plane"
Silent 1						|	Creating the Observation Structure  ...
PauseUpdate

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1]
SrwSmpGenTotName=ObsName

srwUtiSetValN("nx", nx, "SrwSmpPlaneOrient")
srwUtiSetValN("ny", ny, "SrwSmpPlaneOrient")
srwUtiSetValN("nz", nz, "SrwSmpPlaneOrient")

if(exists(ObsName) != 1)
	abort "Radiation Sampling structure was not found"
endif

if(dimsize($ObsName, 0) < 24)
	redimension/N=24 $ObsName
endif

$ObsName[21] = nx
$ObsName[22] = ny
$ObsName[23] = nz

//Make/N=18/D/O $ObsName
//
//$ObsName[0]=0
//$ObsName[4]=Dist
//
// Energy scan
//$ObsName[5]=SrwSmpEdep*1000
//$ObsName[6]=SrwSmpEfin*1000
//$ObsName[7]=SrwSmpEnpts
// Hor. scan
//$ObsName[8]=SrwSmpXmid*0.001
//$ObsName[9]=SrwSmpXmagn*0.001
//$ObsName[10]=SrwSmpXnpts
// Vert. scan
//$ObsName[11]=SrwSmpZmid*0.001
//$ObsName[12]=SrwSmpZmagn*0.001
//$ObsName[13]=SrwSmpZnpts
// Time scan takes positions [15], [16], [17]; [14] defines Time/Frequency representation
end  

//+++++++++++++++++++++++++++++++++++++++
//
//Define Orientation of the Observation Plane
//Complete Version
//+++++++++++++++++++++++++++++++++++++++
proc SrwSmpFrameOrient(ObsName, axx, axy, axz, ayx, ayy, ayz)
string ObsName=SrwSmpGenTotName
variable axx=srwUtiGetValN("axx", 1, "SrwSmpFrameOrient")
variable axy=srwUtiGetValN("axy", 0, "SrwSmpFrameOrient")
variable axz=srwUtiGetValN("axz", 0, "SrwSmpFrameOrient")
variable ayx=srwUtiGetValN("ayx", 0, "SrwSmpFrameOrient")
variable ayy=srwUtiGetValN("ayy", 1, "SrwSmpFrameOrient")
variable ayz=srwUtiGetValN("ayz", 0, "SrwSmpFrameOrient")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt axx, "Hor. Lab-Frame Coord. of \"Hor.\" Ort"
prompt axy, "Long. Lab-Frame Coord. of \"Hor.\" Ort"
prompt axz, "Vert. Lab-Frame Coord. of \"Hor.\" Ort"
prompt ayx, "Hor. Lab-Frame Coord. of Inner Normal"
prompt ayy, "Long. Lab-Frame Coord. of Inner Normal"
prompt ayz, "Vert. Lab-Frame Coord. of Inner Normal"
Silent 1						|	Creating the Observation Structure  ...
PauseUpdate

SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1]
SrwSmpGenTotName=ObsName

srwUtiSetValN("axx", axx, "SrwSmpFrameOrient")
srwUtiSetValN("axy", axy, "SrwSmpFrameOrient")
srwUtiSetValN("axz", axz, "SrwSmpFrameOrient")
srwUtiSetValN("ayx", ayx, "SrwSmpFrameOrient")
srwUtiSetValN("ayy", ayy, "SrwSmpFrameOrient")
srwUtiSetValN("ayz", ayz, "SrwSmpFrameOrient")

if(exists(ObsName) != 1)
	abort "Radiation Sampling structure was not found"
endif

make/O vExP = {axx, axy, axz}, vEyP = {ayx, ayy, ayz}

variable oldNormExP = srwUtiVectNorm(vExP, 1)
variable oldNormEyP = srwUtiVectNorm(vEyP, 1)

if((oldNormExP <= 0) %| (oldNormEyP <= 0))
	abort "Base vector length can't be zero"
endif

variable relTol = 1e-8
variable scalProd = srwUtiVectScalProd(vExP, vEyP)
if(abs(scalProd) > relTol)
	abort "Base vectors are not perpendicular"
endif

if(dimsize($ObsName, 0) < 27)
	redimension/N=27 $ObsName
endif

$ObsName[21] = vExP[0]
$ObsName[22] = vExP[1]
$ObsName[23] = vExP[2]

$ObsName[24] = vEyP[0]
$ObsName[25] = vEyP[1]
$ObsName[26] = vEyP[2]
end  
