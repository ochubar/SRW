
//+++++++++++++++++++++++++++++++++++++++
//
// ID Simulations
//
//+++++++++++++++++++++++++++++++++++++++
//Imports 2D magnetic field from 2 files 
//(can be used to input field along 1 period)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagFieldFromFile(mname,sst,sfi,np,fnamex,fnamez)
string mname=SrwMagName
variable sst=0
variable sfi=50
variable np=101
string fnamex
string fnamez
prompt mname,"Name of Magn. Field structure to create"
prompt sst,"Initial Longitudinal Position [m]"
prompt sfi,"Final Longitudinal Position [m]"
prompt np,"Number of Points"
prompt fnamex,"File Name with Horizontal Field Data"
prompt fnamez,"File Name with Vertical Field Data"

//to improve: allow only one component
SrwMagFieldCreate(mname,0.5*(sst+sfi),sfi-sst,np)
SrwMagImportCmpn(mname+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType,1,fnamex)
SrwMagImportCmpn(mname+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType,1,fnamez)
end

//+++++++++++++++++++++++++++++++++++++++
//Insert periods in center of a given magnetic field
//(most frequently used version)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagInsertPeriods(resMagCore, inMag, per, nper)
string resMagCore=srwUtiGetValS("resMagCore", "", "SrwUtiMagInsertPeriods")
string inMag=srwUtiGetValS("inMag", "", "SrwUtiMagInsertPeriods") + SrwFieldType
variable per=srwUtiGetValN("per", 0.08, "SrwUtiMagInsertPeriods")
variable nper=srwUtiGetValN("nper", 20, "SrwUtiMagInsertPeriods")
prompt resMagCore,"Name of Resulting Magn. Field structure"
prompt inMag,"Name of Short Magn. Field structure",popup Wavelist("*"+SrwFieldType ,";", "")
prompt per,"Period [m]"
prompt nper,"Number of Periods to Insert"
Silent 1						|	Constructing Undulator Field ...
PauseUpdate

srwUtiSetValS("resMagCore", resMagCore, "SrwUtiMagInsertPeriods")
srwUtiSetValS("inMag", inMag[0,strlen(inMag)-strlen(SrwFieldWaveType)-1], "SrwUtiMagInsertPeriods")
srwUtiSetValN("per", per, "SrwUtiMagInsertPeriods")
srwUtiSetValN("nper", nper, "SrwUtiMagInsertPeriods")

variable npWavePer = 1501
make/O/N=(npWavePer) wAuxOnePerBX, wAuxOnePerBZ
SetScale/I x -0.5*per,0.5*per,"m", wAuxOnePerBX, wAuxOnePerBZ

string nmInMagBX = $inMag[0], nmInMagBZ = $inMag[1]
wAuxOnePerBX = $nmInMagBX(x)
wAuxOnePerBZ = $nmInMagBZ(x)

duplicate/O wAuxOnePerBX wAuxIntOnePerBX
integrate/T wAuxIntOnePerBX
duplicate/O wAuxOnePerBZ wAuxIntOnePerBZ
integrate/T wAuxIntOnePerBZ
variable corBX = -wAuxIntOnePerBX[npWavePer - 1]/per
variable corBZ = -wAuxIntOnePerBZ[npWavePer - 1]/per
wAuxOnePerBX += corBX
wAuxOnePerBZ += corBZ

variable inStart = dimoffset($nmInMagBZ, 0)
variable inStep = dimdelta($nmInMagBZ, 0), inSize = dimsize($nmInMagBZ, 0)
variable inRange = inStep*(inSize - 1)
variable inEnd = inStart + inRange
variable extraRange = per*nper
variable outStart = inStart - 0.5*extraRange
variable outEnd = inEnd + 0.5*extraRange
variable outRange = outEnd - outStart
variable outCen = 0.5*(outStart + outEnd)
SrwMagFieldCreate(resMagCore, outCen, outRange, 15000)
string nmResMagBX = resMagCore + "BX_fld", nmResMagBZ = resMagCore + "BZ_fld"

variable halfPer = 0.5*per
variable extraStart = outCen - 0.5*extraRange
variable extraStartPextraRange = extraStart + extraRange

$nmResMagBX = 0
//$nmResMagBX += $nmInMagBX(x - extraStart + halfPer)*srwUtiNonZeroInterval(x, outStart, extraStart)
//$nmResMagBX += $nmInMagBX(x - extraStartPextraRange - halfPer)*srwUtiNonZeroInterval(x, extraStartPextraRange, outEnd)
$nmResMagBX += $nmInMagBX(x - extraStart)*srwUtiNonZeroInterval(x, outStart, extraStart)
$nmResMagBX += $nmInMagBX(x - extraStartPextraRange)*srwUtiNonZeroInterval(x, extraStartPextraRange, outEnd)

$nmResMagBZ = 0
//$nmResMagBZ += $nmInMagBZ(x - extraStart + halfPer)*srwUtiNonZeroInterval(x, outStart, extraStart)
//$nmResMagBZ += $nmInMagBZ(x - extraStartPextraRange - halfPer)*srwUtiNonZeroInterval(x, extraStartPextraRange, outEnd)
$nmResMagBZ += $nmInMagBZ(x - extraStart)*srwUtiNonZeroInterval(x, outStart, extraStart)
$nmResMagBZ += $nmInMagBZ(x - extraStartPextraRange)*srwUtiNonZeroInterval(x, extraStartPextraRange, outEnd)

variable curCen = extraStart + halfPer
//variable curCen = extraStart //+ per

variable iPer=0
do
	//$nmResMagBX += wAuxOnePerBX(x - curCen)*srwUtiNonZeroInterval(x, curCen - halfPer, curCen + halfPer)
	//$nmResMagBZ += wAuxOnePerBZ(x - curCen)*srwUtiNonZeroInterval(x, curCen - halfPer, curCen + halfPer)
	
	$nmResMagBX -= wAuxOnePerBX(x - curCen)*srwUtiNonZeroInterval(x, curCen - halfPer, curCen + halfPer)
	$nmResMagBZ -= wAuxOnePerBZ(x - curCen)*srwUtiNonZeroInterval(x, curCen - halfPer, curCen + halfPer)

	curCen += per
	iPer += 1
while(iPer < nper)

killwaves/Z wAuxOnePerBX, wAuxIntOnePerBX, wAuxOnePerBZ, wAuxIntOnePerBZ
end

//+++++++++++++++++++++++++++++++++++++++
//Add periods to a given magnetic field by serializing it
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagAddPeriods(mname,nperx,nperz,mnamefin)
string mname=SrwMagName+SrwFieldType
variable nperx=20
variable nperz=20
string mnamefin=""
prompt mname,"Name of Magn. Field structure",popup Wavelist("*"+SrwFieldType ,";", "")
prompt nperx,"Number of Hor. Periods in the final Magn. Field"
prompt nperz,"Number of Vert. Periods in the final Magn. Field"
prompt mnamefin,"Name of Magn. Field structure to create"

Silent 1						|	Creating Periodic Magnetic Field  ...
PauseUpdate

//string VertFldCompName=mname[0,strlen(mname)-strlen(SrwFieldType)-1]+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType
//string HorFldCompName=mname[0,strlen(mname)-strlen(SrwFieldType)-1]+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType
string VertFldCompName=$mname[1]
string HorFldCompName=$mname[0]

variable npermax=nperz, npermin=nperz
if(nperx>nperz)
	npermax=nperx
else
	npermin=nperx
endif

variable InitNp=numpnts($VertFldCompName)
variable InitStartS=DimOffset($VertFldCompName, 0)
variable InitStepS=DimDelta($VertFldCompName, 0)

variable FinNp=round(InitNp*npermax)
variable FinFieldLength=InitStepS*(FinNp - 1)
variable FinFieldCenter=InitStartS+0.5*FinFieldLength

SrwMagFieldCreate(mnamefin,FinFieldCenter,FinFieldLength,FinNp)
string FinMagFieldName=mnamefin+SrwFieldType
string FinVertFldCompName=$FinMagFieldName[1] //mnamefin+SrwSuffixMagField+SrwSuffixZ+SrwFieldWaveType
string FinHorFldCompName=$FinMagFieldName[0] //mnamefin+SrwSuffixMagField+SrwSuffixX+SrwFieldWaveType
SrwMagZero(FinVertFldCompName)
SrwMagZero(FinHorFldCompName)

variable NperInt=floor(npermax)
variable PerCount=0, PointCount
do
	PointCount=0
	do
		if(PerCount<floor(nperz))
			$FinVertFldCompName[PerCount*InitNp+PointCount]=$VertFldCompName[PointCount]
		endif
		if(PerCount<floor(nperx))
			$FinHorFldCompName[PerCount*InitNp+PointCount]=$HorFldCompName[PointCount]
		endif
		PointCount += 1
	while(PointCount<InitNp)
	PerCount += 1
while(PerCount<NperInt)

variable NperIntMin=floor(npermin)
variable nperzInitNp=nperz*InitNp, nperxInitNp=nperx*InitNp
if((npermax-NperIntMin)>0.01)
	variable MoreNp=round((npermax-NperIntMin)*InitNp)
	variable iStart=NperIntMin*InitNp
	variable i=0
	variable iStartpi
	do
		iStartpi=iStart+i
		if(iStartpi<=nperzInitNp)
			$FinVertFldCompName[iStartpi]=$VertFldCompName[i]
		endif
		if(iStartpi<=nperxInitNp)
			$FinHorFldCompName[iStartpi]=$HorFldCompName[i]
		endif
		i += 1
	while(i<MoreNp)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Add terminations to a given magnetic field.
//The terminations try to ensure axis of emission at 0 angle,
//and to compensate 1st and 2nd field integrals.
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagAddTermin(mname,ds,dszero)
string mname=SrwMagName+SrwFieldType
variable ds=20
variable dszero=20
prompt mname,"Name of Magn. Field structure",popup Wavelist("*"+SrwFieldType,";","")
prompt ds,"Termination Length [mm]"
prompt dszero,"Zero Field Length [mm]"

string VertFldCompName=$mname[1]
string HorFldCompName=$mname[0]

variable InitNp=numpnts($VertFldCompName)
variable InitStartS=DimOffset($VertFldCompName, 0)
variable InitStepS=DimDelta($VertFldCompName, 0)

variable HalfMoreNp=floor((ds*0.001)/InitStepS)+1
ds=1000.*HalfMoreNp*InitStepS

variable HalfMoreNpZero=floor((dszero*0.001)/InitStepS)+1
variable dszero_m=HalfMoreNpZero*InitStepS
dszero=1000.*dszero_m

variable HalfMoreNpTot=HalfMoreNp+HalfMoreNpZero

variable MoreNp=2*HalfMoreNpTot
variable NewNp=InitNp+MoreNp

string AuxNameX="aux_x", IntAuxNameX="int_aux_x"
string AuxNameZ="aux_z", IntAuxNameZ="int_aux_z"

Duplicate/O $VertFldCompName, $AuxNameZ
Duplicate/O $HorFldCompName, $AuxNameX
$VertFldCompName=0
$HorFldCompName=0

Make/O/D/N=(NewNp) $VertFldCompName, $HorFldCompName
SetScale/P x, InitStartS-dszero_m, InitStepS, "m", $VertFldCompName, $HorFldCompName
SetScale d 0, 0, SrwPUnitMagField, $VertFldCompName, $HorFldCompName

variable FeStX=$AuxNameX[0], FeFiX=$AuxNameX[InitNp-1]
variable FeDStX=($AuxNameX[1]-FeStX)/InitStepS, FeDFiX=(FeFiX-$AuxNameX[InitNp-2])/InitStepS
WaveStats/Q $AuxNameX
//variable FhiX=-0.5*V_avg*InitStepS*InitNp
variable I1X=V_avg*InitStepS*InitNp

Duplicate/O $AuxNameX, $IntAuxNameX
Integrate/T $IntAuxNameX
WaveStats/Q $IntAuxNameX
variable I2X=V_avg*InitStepS*InitNp

variable FeStZ=$AuxNameZ[0], FeFiZ=$AuxNameZ[InitNp-1]
variable FeDStZ=($AuxNameZ[1]-FeStZ)/InitStepS, FeDFiZ=(FeFiZ-$AuxNameZ[InitNp-2])/InitStepS
WaveStats/Q $AuxNameZ
//variable FhiZ=-0.5*V_avg*InitStepS*InitNp
variable I1Z=V_avg*InitStepS*InitNp
//variable CorI1Z = -InitStepS*FeStZ

Duplicate/O $AuxNameZ, $IntAuxNameZ
Integrate/T $IntAuxNameZ
WaveStats/Q $IntAuxNameZ
variable I2Z=V_avg*InitStepS*InitNp

variable L=0.001*ds
variable SeSt=InitStartS+L, SeFi=SeSt+(InitNp-1)*InitStepS
//L=-L
//variable a0StZ=-((L + SeSt)^2*(-2*FeStZ*L*(L^2 - 2*L*SeSt - 15*SeSt^2) + SeSt*(60*FhiZ*SeSt + FeDStZ*L^2*(2*L + 5*SeSt))))/(2.*L^5)
//variable a1StZ=((L + SeSt)*(FeDStZ*L^2*(L^2 + 8*L*SeSt + 10*SeSt^2) + 12*SeSt*(5*FhiZ*(L + 2*SeSt) + FeStZ*L*(3*L + 5*SeSt))))/L^5
//variable a2StZ=(-3*(20*FhiZ*(L^2 + 6*L*SeSt + 6*SeSt^2) + L*(FeDStZ*L*(3*L^2 + 12*L*SeSt + 10*SeSt^2) + 4*FeStZ*(3*L^2 + 16*L*SeSt + 15*SeSt^2))))/(2.*L^5)
//variable a3StZ=(2*(30*FhiZ*(L + 2*SeSt) + L*(FeDStZ*L*(3*L + 5*SeSt) + 2*FeStZ*(8*L + 15*SeSt))))/L^5
//variable a4StZ=(-5*(12*FhiZ + L*(6*FeStZ + FeDStZ*L)))/(2.*L^5)
variable/D a0StZ=SrwUtiAuxA00(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)
variable/D a1StZ=SrwUtiAuxA10(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)
variable/D a2StZ=SrwUtiAuxA20(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)
variable/D a3StZ=SrwUtiAuxA30(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)
variable/D a4StZ=SrwUtiAuxA40(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)
variable/D a5StZ=SrwUtiAuxA50(SeSt, SeFi, L, FeStZ, FeDStZ, I1Z, I2Z)

//variable a0StX=-((L + SeSt)^2*(-2*FeStX*L*(L^2 - 2*L*SeSt - 15*SeSt^2) + SeSt*(60*FhiX*SeSt + FeDStX*L^2*(2*L + 5*SeSt))))/(2.*L^5)
//variable a1StX=((L + SeSt)*(FeDStX*L^2*(L^2 + 8*L*SeSt + 10*SeSt^2) + 12*SeSt*(5*FhiX*(L + 2*SeSt) + FeStX*L*(3*L + 5*SeSt))))/L^5
//variable a2StX=(-3*(20*FhiX*(L^2 + 6*L*SeSt + 6*SeSt^2) + L*(FeDStX*L*(3*L^2 + 12*L*SeSt + 10*SeSt^2) + 4*FeStX*(3*L^2 + 16*L*SeSt + 15*SeSt^2))))/(2.*L^5)
//variable a3StX=(2*(30*FhiX*(L + 2*SeSt) + L*(FeDStX*L*(3*L + 5*SeSt) + 2*FeStX*(8*L + 15*SeSt))))/L^5
//variable a4StX=(-5*(12*FhiX + L*(6*FeStX + FeDStX*L)))/(2.*L^5)
variable/D a0StX=SrwUtiAuxA00(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)
variable/D a1StX=SrwUtiAuxA10(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)
variable/D a2StX=SrwUtiAuxA20(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)
variable/D a3StX=SrwUtiAuxA30(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)
variable/D a4StX=SrwUtiAuxA40(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)
variable/D a5StX=SrwUtiAuxA50(SeSt, SeFi, L, FeStX, FeDStX, I1X, I2X)

//L=-L
//variable a0FiZ=-((L + SeFi)^2*(-2*FeFiZ*L*(L^2 - 2*L*SeFi - 15*SeFi^2) + SeFi*(-60*FhiZ*SeFi + FeDFiZ*L^2*(2*L + 5*SeFi))))/(2.*L^5)
//variable a1FiZ=((L + SeFi)*(FeDFiZ*L^2*(L^2 + 8*L*SeFi + 10*SeFi^2) + 12*SeFi*(-5*FhiZ*(L + 2*SeFi) + FeFiZ*L*(3*L + 5*SeFi))))/L^5
//variable a2FiZ=(-3*(-20*FhiZ*(L^2 + 6*L*SeFi + 6*SeFi^2) + L*(FeDFiZ*L*(3*L^2 + 12*L*SeFi + 10*SeFi^2) + 4*FeFiZ*(3*L^2 + 16*L*SeFi + 15*SeFi^2))))/(2.*L^5)
//variable a3FiZ=(2*(-30*FhiZ*(L + 2*SeFi) + L*(FeDFiZ*L*(3*L + 5*SeFi) + 2*FeFiZ*(8*L + 15*SeFi))))/L^5
//variable a4FiZ=(-5*(-12*FhiZ + L*(6*FeFiZ + FeDFiZ*L)))/(2.*L^5)
variable/D a0FiZ=SrwUtiAuxA0e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)
variable/D a1FiZ=SrwUtiAuxA1e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)
variable/D a2FiZ=SrwUtiAuxA2e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)
variable/D a3FiZ=SrwUtiAuxA3e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)
variable/D a4FiZ=SrwUtiAuxA4e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)
variable/D a5FiZ=SrwUtiAuxA5e(SeSt, SeFi, L, FeFiZ, FeDFiZ, I1Z, I2Z)

//variable a0FiX=-((L + SeFi)^2*(-2*FeFiX*L*(L^2 - 2*L*SeFi - 15*SeFi^2) + SeFi*(-60*FhiX*SeFi + FeDFiX*L^2*(2*L + 5*SeFi))))/(2.*L^5)
//variable a1FiX=((L + SeFi)*(FeDFiX*L^2*(L^2 + 8*L*SeFi + 10*SeFi^2) + 12*SeFi*(-5*FhiX*(L + 2*SeFi) + FeFiX*L*(3*L + 5*SeFi))))/L^5
//variable a2FiX=(-3*(-20*FhiX*(L^2 + 6*L*SeFi + 6*SeFi^2) + L*(FeDFiX*L*(3*L^2 + 12*L*SeFi + 10*SeFi^2) + 4*FeFiX*(3*L^2 + 16*L*SeFi + 15*SeFi^2))))/(2.*L^5)
//variable a3FiX=(2*(-30*FhiX*(L + 2*SeFi) + L*(FeDFiX*L*(3*L + 5*SeFi) + 2*FeFiX*(8*L + 15*SeFi))))/L^5
//variable a4FiX=(-5*(-12*FhiX + L*(6*FeFiX + FeDFiX*L)))/(2.*L^5)
variable/D a0FiX=SrwUtiAuxA0e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)
variable/D a1FiX=SrwUtiAuxA1e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)
variable/D a2FiX=SrwUtiAuxA2e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)
variable/D a3FiX=SrwUtiAuxA3e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)
variable/D a4FiX=SrwUtiAuxA4e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)
variable/D a5FiX=SrwUtiAuxA5e(SeSt, SeFi, L, FeFiX, FeDFiX, I1X, I2X)

//variable vIntZ0 = TestInt0(SeSt, L, a0StZ, a1StZ, a2StZ, a3StZ, a4StZ, a5StZ)
//variable testHalfInt = -0.5*I1Z
//variable vIntZE = TestIntE(SeFi, L, a0FiZ, a1FiZ, a2FiZ, a3FiZ, a4FiZ, a5FiZ)

$VertFldCompName=0
$HorFldCompName=0

variable sZeroOffset = InitStartS-dszero_m
variable iStart = floor((InitStartS - sZeroOffset)/InitStepS + 0.000001)
variable iEndP1 = floor((SeSt - sZeroOffset)/InitStepS + 0.000001) 
$VertFldCompName += srwUtiMagTermFunc(sZeroOffset + p*InitStepS,a0StZ,a1StZ,a2StZ,a3StZ,a4StZ,a5StZ,sZeroOffset + iStart*InitStepS,sZeroOffset + iEndP1*InitStepS)
$HorFldCompName += srwUtiMagTermFunc(sZeroOffset + p*InitStepS,a0StX,a1StX,a2StX,a3StX,a4StX,a5StX,sZeroOffset + iStart*InitStepS,sZeroOffset + iEndP1*InitStepS)

iStart = floor((SeFi - sZeroOffset)/InitStepS + 0.000001)
iEndP1 = floor((SeFi+abs(L) - sZeroOffset)/InitStepS + 0.000001)
$VertFldCompName += srwUtiMagTermFunc(sZeroOffset + p*InitStepS,a0FiZ,a1FiZ,a2FiZ,a3FiZ,a4FiZ,a5FiZ,sZeroOffset + iStart*InitStepS,sZeroOffset + iEndP1*InitStepS)
$HorFldCompName += srwUtiMagTermFunc(sZeroOffset + p*InitStepS,a0FiX,a1FiX,a2FiX,a3FiX,a4FiX,a5FiX,sZeroOffset + iStart*InitStepS,sZeroOffset + iEndP1*InitStepS)

iStart = floor((SeSt - sZeroOffset)/InitStepS + 0.000001)
iEndP1 = floor((SeFi - sZeroOffset)/InitStepS + 0.000001)
$VertFldCompName += srwUtiWaveShifted(p,iStart,iEndP1,$AuxNameZ)
$HorFldCompName += srwUtiWaveShifted(p,iStart,iEndP1,$AuxNameX)

KillWaves/Z $AuxNameX,  $AuxNameZ, $IntAuxNameX, $IntAuxNameZ
end

//+++++++++++++++++++++++++++++++++++++++
//Auxiliary functions
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA0e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return -(((L + xe)^2*(-(Fe*L^5) + L^4*(2*Fe + FeD*L)*xe + 3*L*(70*I2 + L^2*(19*Fe + 2*FeD*L) - 5*I1*(6*L + 7*R))*xe^2 + 7*(60*I2 + L^2*(12*Fe + FeD*L) - 30*I1*(L + R))*xe^3))/L^7)
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA00(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return ((L - x0)^2*(F0*L^2*(L^3 + 2*L^2*x0 - 57*L*x0^2 + 84*x0^3) - x0*(F0D*L^3*(L^2 - 6*L*x0 + 7*x0^2) + 15*x0*(-14*I2*(L - 2*x0) + I1*(8*L^2 + 7*L*(r - 2*x0) - 14*r*x0)))))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA1e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return ((L + xe)*(FeD*L^6 + 15*L^2*(28*I2 + L^2*(8*Fe + FeD*L) - 2*I1*(6*L + 7*R))*xe + 15*L*(140*I2 - 66*I1*L + 32*Fe*L^2 + 3*FeD*L^3 - 70*I1*R)*xe^2 + 35*(60*I2 + L^2*(12*Fe + FeD*L) - 30*I1*(L + R))*xe^3))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA10(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return ((L - x0)*(30*x0*(L^2*(-14*I2 + 8*I1*L + 4*F0*L^2 + 7*I1*r) - L*(-70*I2 + 37*I1*L + 16*F0*L^2 + 35*I1*r)*x0 - 7*(10*I2 - 2*F0*L^2 - 5*I1*(L + r))*x0^2) + F0D*L^3*(L^3 - 15*L^2*x0 + 45*L*x0^2 - 35*x0^3)))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA2e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return -((L^3*(210*I2 + 4*L^2*(15*Fe + 2*FeD*L) - 15*I1*(6*L + 7*R)) + 30*L^2*(84*I2 - 39*I1*L + 20*Fe*L^2 + 2*FeD*L^3 - 42*I1*R)*xe + 30*L*(210*I2 + L^2*(45*Fe + 4*FeD*L) - 3*I1*(34*L + 35*R))*xe^2 + 70*(60*I2 + L^2*(12*Fe + FeD*L) - 30*I1*(L + R))*xe^3)/L^7)
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA20(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return (L^3*(210*I2 - 60*F0*L^2 + 8*F0D*L^3 - 15*I1*(8*L + 7*r)) - 30*L^2*(84*I2 - 45*I1*L - 20*F0*L^2 + 2*F0D*L^3 - 42*I1*r)*x0 + 30*L*(210*I2 + L^2*(-45*F0 + 4*F0D*L) - 3*I1*(36*L + 35*r))*x0^2 - 70*(60*I2 + L^2*(-12*F0 + F0D*L) - 30*I1*(L + r))*x0^3)/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA3e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return (10*(84*I2*(L^2 + 5*L*xe + 5*xe^2) - 3*I1*(13*L^3 + 70*R*xe^2 + 70*L*xe*(R + xe) + 2*L^2*(7*R + 34*xe)) + L^2*(FeD*L*(2*L^2 + 8*L*xe + 7*xe^2) + Fe*(20*L^2 + 90*L*xe + 84*xe^2))))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA30(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return (-10*(-84*I2*(L^2 - 5*L*x0 + 5*x0^2) + 3*I1*(15*L^3 + 2*L^2*(7*r - 36*x0) - 70*L*(r - x0)*x0 + 70*r*x0^2) + L^2*(F0D*L*(-2*L^2 + 8*L*x0 - 7*x0^2) + F0*(20*L^2 - 90*L*x0 + 84*x0^2))))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA4e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return (-5*(210*I2*(L + 2*xe) - 3*I1*(34*L^2 + 70*R*xe + 35*L*(R + 2*xe)) + L^2*(FeD*L*(4*L + 7*xe) + Fe*(45*L + 84*xe))))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA40(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return (5*(210*I2*(L - 2*x0) - 3*I1*(36*L^2 + 35*L*(r - 2*x0) - 70*r*x0) + L^2*(F0D*L*(4*L - 7*x0) + F0*(-45*L + 84*x0))))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA5e(x0, xe, L, Fe, FeD, I1, I2)
variable/D x0, xe, L, Fe, FeD, I1, I2
variable/D R = xe - x0 - 2*L
variable/D d = xe - x0
variable/D Ap = -I1 + I2/d
I1 = -2*Ap
I2 = -Ap*d
return (7*(60*I2 + L^2*(12*Fe + FeD*L) - 30*I1*(L + R)))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D SrwUtiAuxA50(x0, xe, L, F0, F0D, inI1, I2)
variable/D x0, xe, L, F0, F0D, inI1, I2
variable/D r = xe - x0
variable/D I1 = 2*I2/r
return (7*(60*I2 + L^2*(-12*F0 + F0D*L) - 30*I1*(L + r)))/L^7
end
//+++++++++++++++++++++++++++++++++++++++
function/D srwUtiMagTermFunc(x,a0,a1,a2,a3,a4,a5,xSt,xFi)
variable/D x,a0,a1,a2,a3,a4,a5,xSt,xFi
if((x<=(xSt-0.000001)) %| (x>(xFi-0.000001))) 
	return 0
endif
variable/D xe2=x*x
variable/D xe3=xe2*x
variable/D xe4=xe3*x
return a0+a1*x+a2*xe2+a3*xe3+a4*xe4+a5*xe4*x
end
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiWaveShifted(i,ist,ifip1,w)
variable i,ist,ifip1
wave w
if((i<ist) %| (i>=ifip1))
	return 0
endif
return w[i-ist]
end

//+++++++++++++++++++++++++++++++++++++++
//Creates periodic magnetic field from arb. field within one period,
//adds terminations to compensate field integral, and a kicker magnet
//to deflect particle trajectory after the undulator (to avoid the trajectroy passing
//through obsrvation plane)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagAddPerAndTermin(mname,nperx,nperz,lenterm, lenzero,kickon,mnamefin)
string mname=SrwMagName+SrwFieldType
variable nperx=20
variable nperz=20
variable lenterm=120
variable lenzero=200
variable kickon=1
string mnamefin=""
prompt mname,"Name of 1-per. Magn. Field structure",popup Wavelist("*"+SrwFieldType ,";", "")
prompt nperx,"Number of Hor. Periods"
prompt nperz,"Number of Vert. Periods"
prompt lenterm,"Termination Length [mm]"
prompt lenzero,"Zero Field Length [mm]"
prompt kickon,"Deflect Trajectory after Undulator?",popup "Yes;No"
prompt mnamefin,"Name of the final structure"

SrwUtiMagAddPeriods(mname,nperx,nperz,mnamefin)
SrwUtiMagAddTermin(mnamefin+SrwFieldType,lenterm,lenzero)

if(kickon==1)
	string VertFldCompName=$(mnamefin+SrwFieldType)[1]
	variable dsKick=40. //[mm]
	variable BintKick=-3. //[T*mm]
	variable sKickCen=DimOffset($VertFldCompName,0)+DimSize($VertFldCompName,0)*DimDelta($VertFldCompName,0)-0.001*lenzero+2*0.001*dsKick
	SrwMagGsnAng(VertFldCompName,2,sKickCen,dsKick,BintKick)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates undulator radiation flux through a given aperture vs deflection parameters
//The deflection parameter values should Kx, Kz should be specified in a 2D wave: /N=(n,2)
//where n is the number of Kx(Kz) values.
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndFluxVsK(ebmname, photen, angapx, angapz, undlen, per, wknm, phsh, hrange, wname)
string ebmname=srwUtiGetValS("ebmname", "", "SrwUtiUndFluxVsK")
string wknm=srwUtiGetValS("wknm", "", "SrwUtiUndFluxVsK")
variable undlen=srwUtiGetValN("undlen", 2., "SrwUtiUndFluxVsK")
variable per=srwUtiGetValN("per", 80., "SrwUtiUndFluxVsK")
variable phsh=srwUtiGetValN("phsh", 0., "SrwUtiUndFluxVsK")
variable angapx=srwUtiGetValN("angapx", 0.5, "SrwUtiUndFluxVsK")
variable angapz=srwUtiGetValN("angapz", 0.5, "SrwUtiUndFluxVsK")
variable hrange=srwUtiGetValN("hrange", 6, "SrwUtiUndFluxVsK")
variable photen=srwUtiGetValN("photen", 1000., "SrwUtiUndFluxVsK")
string wname=srwUtiGetValS("wname", "spec", "SrwUtiUndFluxVsK")
prompt ebmname,SrwPElecName,popup Wavelist("*"+SrwElecType,";","");
prompt wknm,"Name of Kx and Kz values wave"
prompt undlen,"Undulator length [m]"
prompt per,"Undulator period [mm]"
prompt phsh,"Phase shift bw vert. and horiz. fields [rad]"
prompt angapx, "Horizontal angular aperture [mrad]"
prompt angapz, "Vertical angular aperture [mrad]"
prompt hrange, "Number or UR harmonics"
prompt photen, "Photon energy [eV]"
prompt wname,"Name of the wave to produce"
srwUtiSetValS("ebmname", ebmname, "SrwUtiUndFluxVsK")
srwUtiSetValS("wknm", wknm, "SrwUtiUndFluxVsK")
srwUtiSetValN("undlen", undlen, "SrwUtiUndFluxVsK")
srwUtiSetValN("per", per, "SrwUtiUndFluxVsK")
srwUtiSetValN("phsh", phsh, "SrwUtiUndFluxVsK")
srwUtiSetValN("angapx", angapx, "SrwUtiUndFluxVsK")
srwUtiSetValN("angapz", angapz, "SrwUtiUndFluxVsK")
srwUtiSetValN("hrange", hrange, "SrwUtiUndFluxVsK")
srwUtiSetValN("photen", photen, "SrwUtiUndFluxVsK")
srwUtiSetValS("wname", wname, "SrwUtiUndFluxVsK")

string AuxObsName = "AuxObs"
string AuxMagFldName = "AuxMagFld"
string AuxStokesName = "AuxStokes"

variable Polarization = 7 // 1- linear hor., 2- linear ver., 3- lin. 45, 4- lin. 135, 5- circ. right, 6- circ. left, 7- total 
variable AmOfExtraEnPtsOnOneSide = 30
variable ObsDist = 100
variable LongIntPar = 1., AzimIntPar = 1.
variable PhotEn_keV = photen*0.001
variable ElecEnergy_GeV = $ebmname[0]

variable PhotEnMin_keV = 0.8*PhotEn_keV
variable PhotEnMax_keV = 1.2*PhotEn_keV
variable TotAmOfEnPts = 2*AmOfExtraEnPtsOnOneSide + 1

SrwSmpCreate(AuxObsName,ObsDist)
SrwSmpScanXZE(AuxObsName + SrwSmpType,0,ObsDist*angapx,1,0,ObsDist*angapz,1,PhotEnMin_keV,PhotEnMax_keV,TotAmOfEnPts)

variable AmOfPt = DimSize($wknm, 0)
variable AmOfKDims = 2
if(DimSize($wknm, 1)==0)
	AmOfKDims = 1
endif

make/O/N=(AmOfPt) $wname

SrwUtiTriggerPrint(2)
variable i = 0
variable Kx = 0, Kz = 0, Ktot = 0, NhArg = 1, NhMin = 1, NhMax = 10
do
	if(AmOfKDims == 1)
		Kz = $wknm[i]
	else 
		Kx = $wknm[i][0]
		Kz = $wknm[i][1]
	endif
	
	Ktot = sqrt(Kx*Kx + Kz*Kz)
	NhArg = srwUtiClosestHarmNum(PhotEn_keV,Ktot,per,ElecEnergy_GeV)
	NhMin = trunc(NhArg - 0.5*hrange)
	if(NhMin <= 0)
		NhMin = 1
	endif
	NhMax = round(NhArg + 0.5*hrange)
	if(NhMax <= 0)
		NhMax = 1
	endif

	SrwMagPerCreate2D(AuxMagFldName,per,Kz,Kx,undlen,phsh,1,0,0)
	SrwPerStoCreate(AuxStokesName,ebmname,AuxMagFldName+SrwUndType,AuxObsName+SrwSmpType,NhMin,NhMax,LongIntPar,AzimIntPar,1)
	SrwSto2Int(AuxStokesName+SrwStoType,"Iaux",Polarization,1,PhotEn_keV,1e-09,1e-09,1)
	$wname[i] = $(AuxStokesName+"Iaux"+SrwSeparator+SrwRadEType)[AmOfExtraEnPtsOnOneSide]
	
	i += 1
while(i < AmOfPt)
SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates power density distribution generated by finite-emittance electron beam 
//at its motion in transversely-uniform magnetic field and propagated through
//a rectangular aperture (e.g. at the exit of vacuum chamber) and a drift space.
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiPowAfterSlitAndDrift(RadName, ElecName, MagName, EffMagLen, ApertName, PosApert, ObsName, PrecPar, Meth, IntNameW)
string RadName=srwUtiGetValS("RadName", srwUtiEnsureShortName(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName), "SrwUtiPowAfterSlitAndDrift")
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string MagName=srwUtiGetValS("SrwMagGenTotName", "", "")
variable EffMagLen=srwUtiGetValN("EffMagLen", 2., "SrwUtiPowAfterSlitAndDrift")
variable PosApert=srwUtiGetValN("PosApert", 6., "SrwUtiPowAfterSlitAndDrift")
string ApertName=srwUtiGetValS("SrwBliRectApert", "RectApert", "") + SrwBeamlineType
string ObsName=srwUtiGetValS("SrwSmpGenTotName", "", "")
variable PrecPar=srwUtiGetValN("SrwPowCompPrec", 1., "")
variable Meth=srwUtiGetValN("SrwPowCompMeth", 1, "")
string IntNameW=srwUtiGetValS("IntNameW", "", "SrwUtiPowAfterSlitAndDrift")
prompt RadName,SrwPPowName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType,";","")+Wavelist("*"+SrwUndType,";","")+Wavelist("*"+SrwMagConstType,";","")
prompt EffMagLen, "Eff. Longitudinal Extent of the Source"
prompt PosApert, "Longitudinal Position of the Aperture"
prompt ApertName, "Rectangular Aperture structure", popup Wavelist("*"+SrwBeamlineType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpPowType,";","")+Wavelist("*"+SrwSmpType,";","")
prompt PrecPar,SrwPPowCompPrec
prompt Meth,SrwPPowCompMeth,popup "Near Field;Far Field"
prompt IntNameW, "Name of wave to store total power", popup Wavelist("*",";","")
Silent 1						|	Computing Power Density  ...
PauseUpdate

srwUtiSetValS("RadName", RadName, "SrwUtiPowAfterSlitAndDrift")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "")
string BufMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
string MagType = MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1]
if(cmpstr(MagType,SrwFieldType)==0)
	srwUtiSetValS("SrwMagName", BufMagName, "")
endif
if(cmpstr(MagType,SrwUndType)==0)
	srwUtiSetValS("SrwUndName", BufMagName, "")
endif
if(cmpstr(MagType,SrwMagConstType)==0)
	srwUtiSetValS("SrwMagConstName", BufMagName, "")
endif
srwUtiSetValS("SrwMagGenTotName", MagName, "")
srwUtiSetValN("EffMagLen", EffMagLen, "SrwUtiPowAfterSlitAndDrift")
srwUtiSetValN("PosApert", PosApert, "SrwUtiPowAfterSlitAndDrift")
srwUtiSetValS("SrwBliRectApert", ApertName[0,strlen(ApertName)-strlen(SrwBeamlineType)-1], "")
srwUtiSetValS("SrwBliLast", ApertName[0,strlen(ApertName)-strlen(SrwBeamlineType)-1], "")
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpPowType)-1], "")
srwUtiSetValS("SrwSmpGenTotName", ObsName, "")
srwUtiSetValN("SrwPowCompPrec", PrecPar, "")
srwUtiSetValN("SrwPowCompMeth", Meth, "")
srwUtiSetValS("IntNameW", IntNameW, "SrwUtiPowAfterSlitAndDrift")

variable s0 = srwGetMagFldCenter(MagName)
variable sMagEnd = s0 + 0.5*EffMagLen
variable sObs = srwGetSmpLongPos(ObsName)

if(sObs < sMagEnd)
	Abort "Longitudinal position of the observation plane should be out of the effective length of the magnetic field."
endif
if(PosApert < sMagEnd)
	Abort "Longitudinal position of the aperture should be out of the effective length of the magnetic field."
endif

variable ApertIsBeforeObsPlane = 0
if(PosApert < sObs)
	ApertIsBeforeObsPlane = 1
endif

SrwPowCreate(RadName, ElecName, MagName, ObsName, PrecPar, Meth, 1)

variable x0Apert = srwGetOptApertRectHorPos(ApertName)
variable dxApert = srwGetOptApertRectHorSize(ApertName)
variable z0Apert = srwGetOptApertRectVertPos(ApertName)
variable dzApert = srwGetOptApertRectVertSize(ApertName)
variable rApert = PosApert - s0
variable rObs = sObs - s0
variable rIm = rObs - rApert
variable rElec = PosApert - srwGetElecBeamLongPos(ElecName)
variable MagnFact = rIm/rElec
variable x0Elec = srwGetElecBeamHorPos(ElecName)
variable z0Elec = srwGetElecBeamVertPos(ElecName)
variable sigxElec = srwGetElecBeamHorSizeRMS(ElecName)
variable sigzElec = srwGetElecBeamVertSizeRMS(ElecName)
variable sigxElecIm = sigxElec*MagnFact
variable sigzElecIm = sigzElec*MagnFact

string ExtrPowDensName = RadName + SrwPowType
string AuxWaveName = "AuxUtiPowAfterSlitAndDrift"
duplicate/O $ExtrPowDensName $AuxWaveName
$AuxWaveName = srwUtiPowDensApertFunc(x, EffMagLen, rApert, rObs, dxApert, x0Apert - x0Elec)*srwUtiPowDensApertFunc(y, EffMagLen, rApert, rObs, dzApert, z0Apert - z0Elec)
//SrwUtiConvWaveWithGaus2D(AuxWaveName, sigxElecIm, sigzElecIm)

variable TotPowWaveExists = 0
if(strlen(IntNameW) > 0)
	if(exists(IntNameW) == 1)
		if(DimSize($IntNameW, 0) > 1)
			TotPowWaveExists = 1
		endif
	endif
endif
variable IntFact = DimDelta($ExtrPowDensName, 0)*DimDelta($ExtrPowDensName, 1)*DimSize($ExtrPowDensName, 0)*DimSize($ExtrPowDensName, 1)*(10^6)
if(TotPowWaveExists == 1)
	WaveStats/Q $ExtrPowDensName
	$IntNameW[0] = V_avg*IntFact
endif

$ExtrPowDensName *= $AuxWaveName[p][q]
if(TotPowWaveExists == 1)
	WaveStats/Q $ExtrPowDensName
	$IntNameW[1] = V_avg*IntFact
endif

KillWaves/Z $AuxWaveName
end

//+++++++++++++++++++++++++++++++++++++++
//Estimates useful aperture for undulator radiation
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndUsefulAngApert(OutWaveName, ElecName, MagName, FluxPortion, PhotEnMin, PhotEnMax, MaxAngApX, MaxAngApZ)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "AngApert", "SrwUtiUndUsefulAngApert"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string MagName = srwUtiGetValS("SrwUndName", "Und", "") + SrwUndType
variable FluxPortion=srwUtiGetValN("FluxPortion", 0.9, "SrwUtiUndUsefulAngApert")
variable PhotEnMin=srwUtiGetValN("PhotEnMin", 0.002, "SrwUtiUndUsefulAngApert")
variable PhotEnMax=srwUtiGetValN("PhotEnMax", 1., "SrwUtiUndUsefulAngApert")
variable MaxAngApX=srwUtiGetValN("MaxAngApX", 0.5, "SrwUtiUndUsefulAngApert")
variable MaxAngApZ=srwUtiGetValN("MaxAngApZ", 0.5, "SrwUtiUndUsefulAngApert")
prompt OutWaveName, "Name of wave to store aperture dimentions"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")
prompt FluxPortion,"Flux portion in the useful aperture [bw 0 and 1]"
prompt PhotEnMin,"Minimal photon energy [keV]"
prompt PhotEnMax,"Maximal photon energy [keV]"
prompt MaxAngApX,"Maximal hor. angular aperture [mrad]"
prompt MaxAngApZ,"Maximal vert. angular aperture [mrad]"
Silent 1						|	...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndUsefulAngApert")
srwUtiGetValN("FluxPortion", FluxPortion, "SrwUtiUndUsefulAngApert")
srwUtiGetValN("PhotEnMin", PhotEnMin, "SrwUtiUndUsefulAngApert")
srwUtiGetValN("PhotEnMax", PhotEnMax, "SrwUtiUndUsefulAngApert")
srwUtiGetValN("MaxAngApX", MaxAngApX, "SrwUtiUndUsefulAngApert")
srwUtiGetValN("MaxAngApZ", MaxAngApZ, "SrwUtiUndUsefulAngApert")

//to implement !!!
//SrwPerStoCreate("ElecMagObs","Elec_ebm","Mag1Per_map","Obs_obs",1,7,1,1,2)
end


//+++++++++++++++++++++++++++++++++++++++
//Function defining dependence of amplitudes Bx, Bz and Phase Shift
//on undulator parameters (gap/shift, currents,...)
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiFuncBxBzPhi_PureCoils(t, p1, p2, p3)
variable t //1- Bz, 2- Bx, 3- Phi
variable p1, p2, p3 //relative parameters varying bw 0 and 1

variable Bz1max = 0.105 //0.10741 // [T]
variable Iz1max = 1 //0.5

variable Bz2max = 0.105 //0.10741 // [T]
variable Iz2max = 1 //0.5

variable Bx1max = 0.1012 //0.10741 // [T]
variable Ix1max = 1 //0.5
// p1, p2  - relative currents creating Bz; vary bw 0 and 1
// p3 - relative current creating Bx; vary bw 0 and 1

variable bz1, bz2
if(t == 1) //Bz
	//bz1 = Bz1max*((p1 - 0.5)/Iz1max)
	//bz2 = Bz2max*((p2 - 0.5)/Iz2max)
	bz1 = Bz1max*(p1/Iz1max)
	bz2 = Bz2max*(p2/Iz2max)
	return sqrt(bz1*bz1 + bz2*bz2)
endif
if(t == 2) //Bx
	//return abs(Bx1max*((p3 - 0.5)/Ix1max))
	return abs(Bx1max*(p3/Ix1max))
endif
if(t == 3) //Phi
	//bz1 = Bz1max*((p1 - 0.5)/Iz1max)
	//bz2 = Bz2max*((p2 - 0.5)/Iz2max)
	bz1 = Bz1max*(p1/Iz1max)
	bz2 = Bz2max*(p2/Iz2max)
	
	variable ph0 = 0
	if(p3 < 0)
		ph0 = Pi
	endif
	
	if(bz2 == 0)
		if(bz1 < 0)
			return -0.5*Pi - ph0
		endif
		if(bz1 == 0)
			return -ph0
		else
			return 0.5*Pi - ph0
		endif
	else
		return atan(bz1/bz2) - ph0
	endif
endif
return 0
end

//+++++++++++++++++++++++++++++++++++++++
//Function defining dependence of amplitudes Bx, Bz and Phase Shift
//on undulator parameters (gap/shift, currents,...)
//+++++++++++++++++++++++++++++++++++++++
function SrwUtiFuncBxBzPhi_CoilsS(t, p1, p2, p3)
variable t //1- Bz, 2- Bx, 3- Phi
variable p1, p2, p3 //relative parameters varying bw 0 and 1

variable Bz1max = 0.105 //0.10741 // [T]
variable Iz1max = 1 //0.5

//variable Bz2max = 0.105 //0.10741 // [T]
//variable Iz2max = 1 //0.5

variable PhiMax = 0.5*Pi

variable Bx1max = 0.1012 //0.10741 // [T]
variable Ix1max = 1 //0.5
//p1 - relative current creating Bz; vary bw 0 and 1
//p2 - phase bw currents creating Bz; vary bw 0 and 1
//p3 - relative current creating Bx; vary bw 0 and 1

//variable bz1, bz2
if(t == 1) //Bz
	//bz1 = Bz1max*((p1 - 0.5)/Iz1max)
	//bz2 = Bz2max*((p2 - 0.5)/Iz2max)
	//bz1 = Bz1max*(p1/Iz1max)
	//bz2 = Bz2max*(p2/Iz2max)
	//return sqrt(bz1*bz1 + bz2*bz2)
	return abs(Bz1max*p1/Iz1max)
endif
if(t == 2) //Bx
	//return abs(Bx1max*((p3 - 0.5)/Ix1max))
	return abs(Bx1max*(p3/Ix1max))
endif
if(t == 3) //Phi
	//bz1 = Bz1max*((p1 - 0.5)/Iz1max)
	//bz2 = Bz2max*((p2 - 0.5)/Iz2max)
	//bz1 = Bz1max*(p1/Iz1max)
	//bz2 = Bz2max*(p2/Iz2max)
	
	//variable ph0 = 0
	//if(p3 < 0)
	//	ph0 = Pi
	//endif
	
	//if(bz2 == 0)
	//	if(bz1 < 0)
	//		return -0.5*Pi - ph0
	//	endif
	//	if(bz1 == 0)
	//		return -ph0
	//	else
	//		return 0.5*Pi - ph0
	//	endif
	//else
	//	return atan(bz1/bz2) - ph0
	//endif
	
	return PhiMax //*p2
endif
return 0
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates max. spectral flux of undulator vs photon energy
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndOptimSpecFuncPar(OutWaveName, ElecName, ObsName, UndPer, UndLen, FuncName, Np1, Np2, Np3, MinPolRate)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "OptSpec", "SrwUtiUndOptimSpecFuncPar"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string ObsName = srwUtiGetValS("SrwSmpName", "Und", "") + SrwSmpType
variable UndPer=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen=srwUtiGetValN("SrwLength", 1.6, "")
string FuncName=srwUtiGetValS("FuncName", "SrwUtiFuncBxBzPhi", "SrwUtiUndOptimSpecFuncPar")
variable Np1=srwUtiGetValN("Np1", 50, "SrwUtiUndOptimSpecFuncPar")
variable Np2=srwUtiGetValN("Np2", 50, "SrwUtiUndOptimSpecFuncPar")
variable Np3=srwUtiGetValN("Np3", 1, "SrwUtiUndOptimSpecFuncPar")
variable MinPolRate=srwUtiGetValN("MinPolRate", 0.9, "SrwUtiUndOptimSpecFuncPar")
prompt OutWaveName, "Name of wave to store spectrum"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType,";","")
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt FuncName,"B vs Parameters Function Name"
prompt Np1,"Number of 1st Parameter Values"
prompt Np2,"Number of 2nd Parameter Values"
prompt Np3,"Number of 3rd Parameter Values"
prompt MinPolRate,"Minimum acceptable polarization rate"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndOptimSpecFuncPar")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "") 
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1], "") 
srwUtiSetValN("SrwPeriod", UndPer, "")
srwUtiSetValN("SrwLength", UndLen, "")
srwUtiSetValS("FuncName", FuncName, "SrwUtiUndOptimSpecFuncPar")
srwUtiSetValN("Np1", Np1, "SrwUtiUndOptimSpecFuncPar")
srwUtiSetValN("Np2", Np2, "SrwUtiUndOptimSpecFuncPar")
srwUtiSetValN("Np3", Np3, "SrwUtiUndOptimSpecFuncPar")
srwUtiSetValN("MinPolRate", MinPolRate, "SrwUtiUndOptimSpecFuncPar")

SrwUtiTriggerPrint(2)

if(Np1 < 1)
	Np1 = 1
endif
if(Np2 < 1)
	Np2 = 1
endif
if(Np3 < 1)
	Np3 = 1
endif

variable NumDims = 3
if(Np1 == 1)
	NumDims -= 1
endif
if(Np2 == 1)
	NumDims -= 1
endif
if(Np3 == 1)
	NumDims -= 1
endif

//variable NumDimB = WaveDims($BzWaveName)
//if(NumDimB != WaveDims($BxWaveName))
//	abort "Dimensions of waves containing horizontal and vertical field components are different"
//endif
//
//variable KsAreDefined = 0
//if(NumDimB == 4)
//	KsAreDefined = 1
//endif
//if(NumDimB == 3)
//	KsAreDefined = 1
//endif
//
//variable BxIsDefined = 0
//if(strlen(BxWaveName) > 0)
//	BxIsDefined = 1
//endif
//variable BzIsDefined = 0
//if(strlen(BzWaveName) > 0)
//	BzIsDefined = 1
//endif

//variable NumMagHarm
//
//variable Np, pStart, pStep
//variable Nq = 1, qStart = 0, qStep = 0
//variable NumDims

//if(KsAreDefined == 1)
//
//	NumMagHarm = DimSize($BzWaveName, 1)
//	
//	if(BxIsDefined > 0)
//		Np = DimSize($BxWaveName, NumDimB-1)
//		pStart = DimOffset($BxWaveName, NumDimB-1)
//		pStep = DimDelta($BxWaveName, NumDimB-1)
//	else
//		Np = DimSize($BzWaveName, NumDimB-1)
//		pStart = DimOffset($BzWaveName, NumDimB-1)
//		pStep = DimDelta($BzWaveName, NumDimB-1)
//	endif
//
//	if(BzIsDefined > 0)
//		Nq = DimSize($BzWaveName, NumDimB-1)
//		qStart = DimOffset($BzWaveName, NumDimB-1)
//		qStep = DimDelta($BzWaveName, NumDimB-1)
//	else
//		Nq = DimSize($BxWaveName, NumDimB-1)
//		qStart = DimOffset($BxWaveName, NumDimB-1)
//		qStep = DimDelta($BxWaveName, NumDimB-1)
//	endif
//	
//	if(NumDimB == 3)
//		NumDims = 1
//	endif
//	if(NumDimB == 4)
//		NumDims = 2
//	endif
//	
//else
//	Np = DimSize($BzWaveName, 0)
//	pStart = DimOffset($BzWaveName, 0)
//	pStep = DimDelta($BzWaveName, 0)
//
//	NumDims = WaveDims($BzWaveName)
//	if(NumDims > 1)
//		Nq = DimSize($BzWaveName, 1)
//		qStart = DimOffset($BzWaveName, 1)
//		qStep = DimDelta($BzWaveName, 1)
//	endif
//endif

//variable Nu = Np*Nq
//if(NumDims == 1)
//	if(Np < Nq)
//		Nu = Np
//	else
//		Nu = Nq
//	endif
//endif

variable Nu = Np1*Np2*Np3

variable Ne = srwGetSmpPhotEnNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
string eUnitsStr = "eV"
variable ElecEnergy = srwGetElecBeamEnergy(ElecName)

variable RegectK = 0.08 //100

variable AmOfExtraHarm = 3 //6 // to make input variable ?
variable PrecPar = 1.5
variable ShowIntermGraphs = 2 // 1- No,  2- Yes
variable AuxShowGraphs = 2

string UndName = "AuxUnd", StokesName = "AuxSto", AuxSpecName, AuxRateName
make/O/N=(Nu, Ne) AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL
make/O/N=(Nu, Ne) AuxRatesLH, AuxRatesLV, AuxRatesCR, AuxRatesCL

AuxSpectraLH = 0
AuxSpectraLV = 0
AuxSpectraCR = 0
AuxSpectraCL = 0

AuxRatesLH = 0
AuxRatesLV = 0
AuxRatesCR = 0
AuxRatesCL = 0

variable PolTypeLH = 1, PolTypeLV = 2, PolTypeCR = 5, PolTypeCL = 6
string AuxSufStr = ""

//variable Bx, Bz, Kx, Kz, MaxHarm
variable/G gAuxBx, gAuxBz, gAuxPhi
variable Kx, Kz, MaxHarm

variable MagHarmCount=0, PhiX, PhiZ, KeffE2, FundPhotEn, Kz1, Kx1
variable ip1 = 0, ip2 = 0, ip3 = 0, iUnd = 0
variable p1 = 0, p2 = 0, p3 = 0

variable sp1 = 0, sp2 = 0, sp3 = 0
if(Np1 > 1)
	sp1 = 1/(Np1 - 1)
endif
if(Np2 > 1)
	sp2 = 1/(Np2 - 1)
endif
if(Np3 > 1)
	sp3 = 1/(Np3 - 1)
endif

variable FieldWasSet = 0, TestMaxHarm, IsRegected = 0

string ComLineStr = ""

do
	p3 = ip3*sp3
	ip2 = 0
	do
		p2 = ip2*sp2
		ip1 = 0
		do
			p1 = ip1*sp1
			
			FieldWasSet = 0
		 
			//if(KsAreDefined == 1)
			//
			//	Kz = $BzWaveName[0][0][ip][iq]
			//	Kx = $BxWaveName[0][0][ip][iq]
			//	PhiX = $BxWaveName[1][0][ip][iq]
			//
			//	Kz1 = Kz
			//	Kx1 = Kx
			//
			//	IsRegected = 1
			//	if((Kx > 0) %| (Kz > 0))
			//		If(Kz <=  RegectK*Kx)
			//			IsRegected = 0
			//		else
			//			if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
			//				IsRegected = 0
			//			else
			//				if(Kz >= (1/RegectK)*Kx)
			//					IsRegected = 0
			//				endif
			//			endif
			//		endif
			//	endif
			//
			//	if(IsRegected != 1)
			//
			//		SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase+PhiX,1,0,0)
			//
			//		KeffE2 = Kx*Kx + Kz*Kz
			//
			//		MagHarmCount=1
			//		do
			//			Kz = $BzWaveName[0][MagHarmCount][ip][iq]
			//			if(Kz > 0)
			//				PhiZ = $BzWaveName[1][MagHarmCount][ip][iq]
			//				SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,1,Kz,PhiZ)
			//				KeffE2 += Kz*Kz
			//			endif
			//			Kx = $BxWaveName[0][MagHarmCount][ip][iq]
			//			if(Kx > 0)
			//				PhiX = $BxWaveName[1][MagHarmCount][ip][iq]
			//				SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,2,Kx,PhiX+Phase)
			//				KeffE2 += Kx*Kx
			//			endif
			//			MagHarmCount += 1
			//		while(MagHarmCount < NumMagHarm)
			//
			//		FundPhotEn = 950*ElecEnergy*ElecEnergy/(1 + 0.5*KeffE2)/(UndPer*0.1)
			//		MaxHarm = round(eEnd/FundPhotEn) //+ AmOfExtraHarm
			//	
			//		TestMaxHarm = round(1.5*MaxHarm)
			//		if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
			//			MaxHarm += AmOfExtraHarm
			//		else 
			//			MaxHarm = TestMaxHarm
			//		endif
			//
			//		FieldWasSet = 1
			//	endif
			//
			//else
			
			//Bx = $BxWaveName[ip][iq]
			//Bz = $BzWaveName[ip][iq]
			
			sprintf ComLineStr, "gAuxBz=%s(1,%g,%g,%g)", FuncName, p1, p2, p3
			execute ComLineStr
			sprintf ComLineStr, "gAuxBx=%s(2,%g,%g,%g)", FuncName, p1, p2, p3
			execute ComLineStr
			sprintf ComLineStr, "gAuxPhi=%s(3,%g,%g,%g)", FuncName, p1, p2, p3
			execute ComLineStr
		
			if(gAuxBz < 0)
				gAuxBz = 0
			endif
			if(gAuxBx < 0)
				gAuxBx = 0
			endif
		
			Kx = srUtiUndK(gAuxBx, UndPer*0.001)
			Kz = srUtiUndK(gAuxBz, UndPer*0.001)
			
			IsRegected = 1
			if((Kx > 0) %| (Kz > 0))
				If(Kz <=  RegectK*Kx)
					IsRegected = 0
				else
					if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
						IsRegected = 0
					else
						if(Kz >= (1/RegectK)*Kx)
							IsRegected = 0
						endif
					endif
				endif
			endif
			
			if(IsRegected != 1)
			
				SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,gAuxPhi,1,0,0)
				MaxHarm = round(eEnd/srUtiUndFundPhotEn(sqrt(gAuxBx*gAuxBx + gAuxBz*gAuxBz), UndPer*0.001, ElecEnergy, 2)) //+ AmOfExtraHarm
				if(MaxHarm == 0)
					MaxHarm = 1
				endif
				
				TestMaxHarm = round(1.5*MaxHarm)
				if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
					MaxHarm += AmOfExtraHarm
				else 
					MaxHarm = TestMaxHarm
				endif
				
				 FieldWasSet = 1
			endif
		//endif
		
			if(FieldWasSet  > 0)
		
				SrwPerStoCreate(StokesName,ElecName,UndName + SrwUndType,ObsName,1,MaxHarm,PrecPar,PrecPar,1)
		
				if(ShowIntermGraphs == 2)
					//if(KsAreDefined == 1)
					//	print ip, iq, "  Kx1 =", Kx1, "   ", "Kz1 =", Kz1
					//else
					print ip1, ip2, ip3, "  Bx =", gAuxBx, "T  ", "Bz =", gAuxBz, "T", "Phi =", gAuxPhi
					//endif
				endif

				AuxSufStr = "I" + num2str(PolTypeLH)
				SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,40,300,150,0,0); DoUpdate
				endif
				AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxSpectraLH[iUnd] = $AuxSpecName[q]
		
				AuxSufStr = "I" + num2str(PolTypeLV)
				SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,180,300,150,0,0); DoUpdate
				endif
				AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxSpectraLV[iUnd] = $AuxSpecName[q]
	
				AuxSufStr = "I" + num2str(PolTypeCR)
				SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,320,300,150,0,0); DoUpdate
				endif
				AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxSpectraCR[iUnd] = $AuxSpecName[q]
	
				AuxSufStr = "I" + num2str(PolTypeCL)
				SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,460,300,150,0,0); DoUpdate
				endif
				AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxSpectraCL[iUnd] = $AuxSpecName[q]
		
				AuxSufStr = "R" + num2str(PolTypeLH)
				SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,40,300,150,0,0); DoUpdate
				endif
				AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxRatesLH[iUnd] = $AuxRateName[q]
	
				AuxSufStr = "R" + num2str(PolTypeLV)
				SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,180,300,150,0,0); DoUpdate
				endif
				AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxRatesLV[iUnd] = $AuxRateName[q]

				AuxSufStr = "R" + num2str(PolTypeCR)
				SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,320,300,150,0,0); DoUpdate
				endif
				AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxRatesCR[iUnd] = $AuxRateName[q]

				AuxSufStr = "R" + num2str(PolTypeCL)
				SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
				if(ShowIntermGraphs == 2)
					SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,460,300,150,0,0); DoUpdate
				endif
				AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
				AuxRatesCL[iUnd] = $AuxRateName[q]
		
				if(ShowIntermGraphs == 2)
					if(AuxShowGraphs == 2)
						AuxShowGraphs = 1
					endif
				endif
		
			endif // FieldWasSet
	
			iUnd += 1
			ip1 += 1
		while(ip1 < Np1)
		ip2 += 1
	while(ip2 < Np2)
	ip3 += 1
while(ip3 < Np3)

//Looking for maximum for each photon energy value

string OutWaveNameLH = OutWaveName + "LH"
string OutWaveNameLV = OutWaveName + "LV"
string OutWaveNameCR = OutWaveName + "CR"
string OutWaveNameCL = OutWaveName + "CL"

string OutWaveNameLHPar1 = OutWaveName + "LH_par1"
string OutWaveNameLVPar1 = OutWaveName + "LV_par1"
string OutWaveNameCRPar1 = OutWaveName + "CR_par1"
string OutWaveNameCLPar1 = OutWaveName + "CL_par1"

string OutWaveNameLHPar2 = OutWaveName + "LH_par2"
string OutWaveNameLVPar2 = OutWaveName + "LV_par2"
string OutWaveNameCRPar2 = OutWaveName + "CR_par2"
string OutWaveNameCLPar2 = OutWaveName + "CL_par2"

string OutWaveNameLHPar3 = OutWaveName + "LH_par3"
string OutWaveNameLVPar3 = OutWaveName + "LV_par3"
string OutWaveNameCRPar3 = OutWaveName + "CR_par3"
string OutWaveNameCLPar3 = OutWaveName + "CL_par3"

make/O/N=(Ne) $OutWaveNameLH, $OutWaveNameLV, $OutWaveNameCR, $OutWaveNameCL
make/O/N=(Ne) $OutWaveNameLHPar1, $OutWaveNameLVPar1, $OutWaveNameCRPar1, $OutWaveNameCLPar1

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLH
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLV
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCR
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCL

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar1

if(NumDims > 1)
	make/O/N=(Ne) $OutWaveNameLHPar2, $OutWaveNameLVPar2, $OutWaveNameCRPar2, $OutWaveNameCLPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar2
endif
if(NumDims > 3)
	make/O/N=(Ne) $OutWaveNameLHPar3, $OutWaveNameLVPar3, $OutWaveNameCRPar3, $OutWaveNameCLPar3
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar3
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar3
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar3
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar3
endif

variable MaxFluxLH, CurFluxLH, CurPolRateLH
variable MaxFluxLV, CurFluxLV, CurPolRateLV
variable MaxFluxCR, CurFluxCR, CurPolRateCR
variable MaxFluxCL, CurFluxCL, CurPolRateCL

variable ie=0, iUndMaxLH=-1, iUndMaxLV=-1, iUndMaxCR=-1, iUndMaxCL=-1
do 
	iUnd = 0
	MaxFluxLH = 0; MaxFluxLV = 0; MaxFluxCR = 0; MaxFluxCL = 0
	do
		CurPolRateLH = AuxRatesLH[iUnd][ie]
		if(CurPolRateLH >= MinPolRate)
			CurFluxLH = AuxSpectraLH[iUnd][ie]
			if(MaxFluxLH < CurFluxLH) 
				MaxFluxLH = CurFluxLH
				iUndMaxLH = iUnd
			endif
		endif
		CurPolRateLV = AuxRatesLV[iUnd][ie]
		if(CurPolRateLV >= MinPolRate)
			CurFluxLV = AuxSpectraLV[iUnd][ie]
			if(MaxFluxLV < CurFluxLV) 
				MaxFluxLV = CurFluxLV
				iUndMaxLV = iUnd
			endif
		endif
		CurPolRateCR = AuxRatesCR[iUnd][ie]
		if(CurPolRateCR >= MinPolRate)
			CurFluxCR = AuxSpectraCR[iUnd][ie]
			if(MaxFluxCR < CurFluxCR) 
				MaxFluxCR = CurFluxCR
				iUndMaxCR = iUnd
			endif
		endif
		CurPolRateCL = AuxRatesCL[iUnd][ie]
		if(CurPolRateCL >= MinPolRate)
			CurFluxCL = AuxSpectraCL[iUnd][ie]
			if(MaxFluxCL < CurFluxCL) 
				MaxFluxCL = CurFluxCL
				iUndMaxCL = iUnd
			endif
		endif
		iUnd += 1
	while(iUnd < Nu)
	
	$OutWaveNameLH[ie] = MaxFluxLH
	$OutWaveNameLV[ie] = MaxFluxLV
	$OutWaveNameCR[ie] = MaxFluxCR
	$OutWaveNameCL[ie] = MaxFluxCL
	
//	ip3 = trunc(iUndMaxLH/(Np1*Np2) + 1e-07)
//	ip = iUndMaxLH - iq*Np
//	$OutWaveNameLHPar1[ie] = pStart + pStep*ip
//	if(NumDims > 1)
//		$OutWaveNameLHPar2[ie] = qStart + qStep*iq
//	endif
//	
//	iq = trunc(iUndMaxLV/Np + 1e-07)
//	ip = iUndMaxLV - iq*Np
//	$OutWaveNameLVPar1[ie] = pStart + pStep*ip
//	if(NumDims > 1)
//		$OutWaveNameLVPar2[ie] = qStart + qStep*iq
//	endif
//	
//	iq = trunc(iUndMaxCR/Np + 1e-07)
//	ip = iUndMaxCR - iq*Np
//	$OutWaveNameCRPar1[ie] = pStart + pStep*ip
//	if(NumDims > 1)
//		$OutWaveNameCRPar2[ie] = qStart + qStep*iq
//	endif
//	
//	iq = trunc(iUndMaxCL/Np + 1e-07)
//	ip = iUndMaxCL - iq*Np
//	$OutWaveNameCLPar1[ie] = pStart + pStep*ip
//	if(NumDims > 1)
//		$OutWaveNameCLPar2[ie] = qStart + qStep*iq
//	endif
	
	ie += 1
while(ie < Ne)

if(ShowIntermGraphs == 2)
	Display $OutWaveNameLH
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Hor. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLV
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Vert. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCR
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Right Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCL
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Left Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
//	Display $OutWaveNameLHPar1
//	Label bottom SrwPLabelPhotEn
//	Label left "Parameter 1 LH"
//	SrwUtiGraphAddFrameAndGrid()
//	
//	Display $OutWaveNameLVPar1
//	Label bottom SrwPLabelPhotEn
//	Label left "Parameter 1 LV"
//	SrwUtiGraphAddFrameAndGrid()
//	
//	Display $OutWaveNameCRPar1
//	Label bottom SrwPLabelPhotEn
//	Label left "Parameter 1 CR"
//	SrwUtiGraphAddFrameAndGrid()
//
//	Display $OutWaveNameCLPar1
//	Label bottom SrwPLabelPhotEn
//	Label left "Parameter 1 CL"
//	SrwUtiGraphAddFrameAndGrid()
//	
//	if(NumDims > 1)
//		Display $OutWaveNameLHPar2
//		Label bottom SrwPLabelPhotEn
//		Label left "Parameter 2 LH"
//		SrwUtiGraphAddFrameAndGrid()
//		
//		Display $OutWaveNameLVPar2
//		Label bottom SrwPLabelPhotEn
//		Label left "Parameter 2 LV"
//		SrwUtiGraphAddFrameAndGrid()
//		
//		Display $OutWaveNameCRPar2
//		Label bottom SrwPLabelPhotEn
//		Label left "Parameter 2 CR"
//		SrwUtiGraphAddFrameAndGrid()
//		
//		Display $OutWaveNameCLPar2
//		Label bottom SrwPLabelPhotEn
//		Label left "Parameter 2 CL"
//		SrwUtiGraphAddFrameAndGrid()
//	endif
endif

KillWaves/Z AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL, AuxPolRateLH, AuxPolRateLV, AuxPolRateCR, AuxPolRateCL, $UndName
SrwUtiTriggerPrint(1)
end


//+++++++++++++++++++++++++++++++++++++++
//Calculates max. spectral flux of undulator vs photon energy
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndOptimSpec(OutWaveName, ElecName, ObsName, UndPer, UndLen, Phase, BxWaveName, BzWaveName, MinPolRate)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "OptSpec", "SrwUtiUndOptimSpec"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string ObsName = srwUtiGetValS("SrwSmpName", "Und", "") + SrwSmpType
variable UndPer=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen=srwUtiGetValN("SrwLength", 1.6, "")
variable Phase=srwUtiGetValN("SrwPh0x", 1.6, "")
string BxWaveName = srwUtiGetValS("BxWaveName", " ", "SrwUtiUndOptimSpec")
string BzWaveName = srwUtiGetValS("BzWaveName", " ", "SrwUtiUndOptimSpec")
variable MinPolRate=srwUtiGetValN("MinPolRate", 0.95, "SrwUtiUndOptimSpec")
prompt OutWaveName, "Name of wave to store spectrum"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType,";","")
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt Phase,"Phase Shift [rad]"
prompt BxWaveName,"Wave of Hor. Mag. Field values", popup Wavelist("*",";","")
prompt BzWaveName,"Wave of Vert. Mag. Field values", popup Wavelist("*",";","")
prompt MinPolRate,"Minimum acceptable polarization rate"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "") 
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1], "") 
srwUtiSetValN("SrwPeriod", UndPer, "")
srwUtiSetValN("SrwLength", UndLen, "")
srwUtiSetValN("SrwPh0x", Phase, "")
srwUtiSetValS("BxWaveName", BxWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("BzWaveName", BzWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValN("MinPolRate", MinPolRate, "SrwUtiUndOptimSpec")

SrwUtiTriggerPrint(2)

variable NumDimB = WaveDims($BzWaveName)
if(NumDimB != WaveDims($BxWaveName))
	abort "Dimensions of waves containing horizontal and vertical field components are different"
endif
//if((NumDimB != 2) &% (NumDimB != 4))
//	abort "The waves containing horizontal and vertical field components are incorrect"
//endif

variable KsAreDefined = 0
if(NumDimB == 4)
	KsAreDefined = 1
endif
if(NumDimB == 3)
	KsAreDefined = 1
endif

variable BxIsDefined = 0
if(strlen(BxWaveName) > 0)
	BxIsDefined = 1
endif
variable BzIsDefined = 0
if(strlen(BzWaveName) > 0)
	BzIsDefined = 1
endif

variable NumMagHarm

variable Np, pStart, pStep
variable Nq = 1, qStart = 0, qStep = 0
variable NumDims

if(KsAreDefined == 1)

	NumMagHarm = DimSize($BzWaveName, 1)
	
	if(BxIsDefined > 0)
		Np = DimSize($BxWaveName, NumDimB-1)
		pStart = DimOffset($BxWaveName, NumDimB-1)
		pStep = DimDelta($BxWaveName, NumDimB-1)
	else
		Np = DimSize($BzWaveName, NumDimB-1)
		pStart = DimOffset($BzWaveName, NumDimB-1)
		pStep = DimDelta($BzWaveName, NumDimB-1)
	endif

	if(BzIsDefined > 0)
		Nq = DimSize($BzWaveName, NumDimB-1)
		qStart = DimOffset($BzWaveName, NumDimB-1)
		qStep = DimDelta($BzWaveName, NumDimB-1)
	else
		Nq = DimSize($BxWaveName, NumDimB-1)
		qStart = DimOffset($BxWaveName, NumDimB-1)
		qStep = DimDelta($BxWaveName, NumDimB-1)
	endif
	
	if(NumDimB == 3)
		NumDims = 1
	endif
	if(NumDimB == 4)
		NumDims = 2
	endif
	
else
	Np = DimSize($BzWaveName, 0)
	pStart = DimOffset($BzWaveName, 0)
	pStep = DimDelta($BzWaveName, 0)

	NumDims = WaveDims($BzWaveName)
	if(NumDims > 1)
		Nq = DimSize($BzWaveName, 1)
		qStart = DimOffset($BzWaveName, 1)
		qStep = DimDelta($BzWaveName, 1)
	endif
endif

variable Nu = Np*Nq
if(NumDims == 1)
	if(Np < Nq)
		Nu = Np
	else
		Nu = Nq
	endif
endif

variable Ne = srwGetSmpPhotEnNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
string eUnitsStr = "eV"

variable ElecEnergy = srwGetElecBeamEnergy(ElecName)

variable RegectK = 100 //0.08

variable AmOfExtraHarm = 6 //6 // to make input variable ?
variable PrecPar = 1.5
variable ShowIntermGraphs = 2 // 1- No,  2- Yes
variable AuxShowGraphs = 2

string UndName = "AuxUnd", StokesName = "AuxSto", AuxSpecName, AuxRateName
make/O/N=(Nu, Ne) AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL
make/O/N=(Nu, Ne) AuxRatesLH, AuxRatesLV, AuxRatesCR, AuxRatesCL

AuxSpectraLH = 0
AuxSpectraLV = 0
AuxSpectraCR = 0
AuxSpectraCL = 0

AuxRatesLH = 0
AuxRatesLV = 0
AuxRatesCR = 0
AuxRatesCL = 0

variable PolTypeLH = 1, PolTypeLV = 2, PolTypeCR = 5, PolTypeCL = 6
string AuxSufStr = ""

variable Bx, Bz, Kx, Kz, MaxHarm
variable MagHarmCount=0, PhiX, PhiZ, KeffE2, FundPhotEn, Kz1, Kx1
variable iq = 0, ip = 0, iUnd = 0
variable FieldWasSet = 0, TestMaxHarm, IsRegected = 0
do
	ip = 0
	do
		 FieldWasSet = 0
		 
		if(KsAreDefined == 1)
		
			Kz = $BzWaveName[0][0][ip][iq]
			Kx = $BxWaveName[0][0][ip][iq]
			PhiX = $BxWaveName[1][0][ip][iq]
			
			Kz1 = Kz
			Kx1 = Kx
			
			IsRegected = 1
			if((Kx > 0) %| (Kz > 0))
				If(Kz <=  RegectK*Kx)
					IsRegected = 0
				else
					if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
						IsRegected = 0
					else
						if(Kz >= (1/RegectK)*Kx)
							IsRegected = 0
						endif
					endif
				endif
			endif
			
			if(IsRegected != 1)
			
				SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase+PhiX,1,0,0)
			
				KeffE2 = Kx*Kx + Kz*Kz
			
				MagHarmCount=1
				do
					Kz = $BzWaveName[0][MagHarmCount][ip][iq]
					if(Kz > 0)
						PhiZ = $BzWaveName[1][MagHarmCount][ip][iq]
						SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,1,Kz,PhiZ)
						KeffE2 += Kz*Kz
					endif
					Kx = $BxWaveName[0][MagHarmCount][ip][iq]
					if(Kx > 0)
						PhiX = $BxWaveName[1][MagHarmCount][ip][iq]
						SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,2,Kx,PhiX+Phase)
						KeffE2 += Kx*Kx
					endif
					MagHarmCount += 1
				while(MagHarmCount < NumMagHarm)
			
				FundPhotEn = 950*ElecEnergy*ElecEnergy/(1 + 0.5*KeffE2)/(UndPer*0.1)
				MaxHarm = round(eEnd/FundPhotEn) //+ AmOfExtraHarm
				
				TestMaxHarm = round(1.5*MaxHarm)
				if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
					MaxHarm += AmOfExtraHarm
				else 
					MaxHarm = TestMaxHarm
				endif

				FieldWasSet = 1
			endif
			
		else
			Bx = $BxWaveName[ip][iq]
			Bz = $BzWaveName[ip][iq]
		
			if(Bx < 0)
				Bx = 0
			endif
			if(Bz < 0)
				Bz = 0
			endif
		
			Kx = srUtiUndK(Bx, UndPer*0.001)
			Kz = srUtiUndK(Bz, UndPer*0.001)
			
			IsRegected = 1
			if((Kx > 0) %| (Kz > 0))
				If(Kz <=  RegectK*Kx)
					IsRegected = 0
				else
					if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
						IsRegected = 0
					else
						if(Kz >= (1/RegectK)*Kx)
							IsRegected = 0
						endif
					endif
				endif
			endif
			
			if(IsRegected != 1)
			
				SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase,1,0,0)
				MaxHarm = round(eEnd/srUtiUndFundPhotEn(sqrt(Bx*Bx + Bz*Bz), UndPer*0.001, ElecEnergy, 2)) //+ AmOfExtraHarm
				
				TestMaxHarm = round(1.5*MaxHarm)
				if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
					MaxHarm += AmOfExtraHarm
				else 
					MaxHarm = TestMaxHarm
				endif
				
				 FieldWasSet = 1
			endif
		endif
		
		if(FieldWasSet  > 0)
		
		SrwPerStoCreate(StokesName,ElecName,UndName + SrwUndType,ObsName,1,MaxHarm,PrecPar,PrecPar,1)
		
		if(ShowIntermGraphs == 2)

			if(KsAreDefined == 1)
				print ip, iq, "  Kx1 =", Kx1, "   ", "Kz1 =", Kz1
			else
				print ip, iq, "  Bx =", Bx, "T  ", "Bz =", Bz, "T"
			endif
			
			//if((ip == 0) %& (iq == 0))
			//	AuxShowGraphs = 2
			//else
			//	AuxShowGraphs = 1
			//endif
		endif

		//SrwSto2Int(StokesName + SrwStoType,"I",PolType,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		
		AuxSufStr = "I" + num2str(PolTypeLH)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,40,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraLH[iUnd] = $AuxSpecName[q]
		
		AuxSufStr = "I" + num2str(PolTypeLV)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,180,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraLV[iUnd] = $AuxSpecName[q]
	
		AuxSufStr = "I" + num2str(PolTypeCR)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,320,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraCR[iUnd] = $AuxSpecName[q]
	
		AuxSufStr = "I" + num2str(PolTypeCL)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,460,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraCL[iUnd] = $AuxSpecName[q]
		
		//SrwSto2PolRate(StokesName + SrwStoType,"R",PolType,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		
		AuxSufStr = "R" + num2str(PolTypeLH)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,40,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesLH[iUnd] = $AuxRateName[q]
	
		AuxSufStr = "R" + num2str(PolTypeLV)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,180,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesLV[iUnd] = $AuxRateName[q]

		AuxSufStr = "R" + num2str(PolTypeCR)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,320,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesCR[iUnd] = $AuxRateName[q]

		AuxSufStr = "R" + num2str(PolTypeCL)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,460,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesCL[iUnd] = $AuxRateName[q]
		
		if(ShowIntermGraphs == 2)
			if(AuxShowGraphs == 2)
				AuxShowGraphs = 1
			endif
		endif
		
		endif // FieldWasSet
	
		iUnd += 1
		ip += 1
	while(ip < Np)
	iq += 1
while(iq < Nq)

//Looking for maximum for each photon energy value

string OutWaveNameLH = OutWaveName + "LH"
string OutWaveNameLV = OutWaveName + "LV"
string OutWaveNameCR = OutWaveName + "CR"
string OutWaveNameCL = OutWaveName + "CL"

string OutWaveNameLHPar1 = OutWaveName + "LH_par1"
string OutWaveNameLVPar1 = OutWaveName + "LV_par1"
string OutWaveNameCRPar1 = OutWaveName + "CR_par1"
string OutWaveNameCLPar1 = OutWaveName + "CL_par1"

string OutWaveNameLHPar2 = OutWaveName + "LH_par2"
string OutWaveNameLVPar2 = OutWaveName + "LV_par2"
string OutWaveNameCRPar2 = OutWaveName + "CR_par2"
string OutWaveNameCLPar2 = OutWaveName + "CL_par2"

make/O/N=(Ne) $OutWaveNameLH, $OutWaveNameLV, $OutWaveNameCR, $OutWaveNameCL
make/O/N=(Ne) $OutWaveNameLHPar1, $OutWaveNameLVPar1, $OutWaveNameCRPar1, $OutWaveNameCLPar1

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLH
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLV
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCR
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCL

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar1

if(NumDims > 1)
	make/O/N=(Ne) $OutWaveNameLHPar2, $OutWaveNameLVPar2, $OutWaveNameCRPar2, $OutWaveNameCLPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar2
endif

variable MaxFluxLH, CurFluxLH, CurPolRateLH
variable MaxFluxLV, CurFluxLV, CurPolRateLV
variable MaxFluxCR, CurFluxCR, CurPolRateCR
variable MaxFluxCL, CurFluxCL, CurPolRateCL

variable ie=0, iUndMaxLH=-1, iUndMaxLV=-1, iUndMaxCR=-1, iUndMaxCL=-1
do 
	iUnd = 0
	MaxFluxLH = 0; MaxFluxLV = 0; MaxFluxCR = 0; MaxFluxCL = 0
	do
		CurPolRateLH = AuxRatesLH[iUnd][ie]
		if(CurPolRateLH >= MinPolRate)
			CurFluxLH = AuxSpectraLH[iUnd][ie]
			if(MaxFluxLH < CurFluxLH) 
				MaxFluxLH = CurFluxLH
				iUndMaxLH = iUnd
			endif
		endif
		CurPolRateLV = AuxRatesLV[iUnd][ie]
		if(CurPolRateLV >= MinPolRate)
			CurFluxLV = AuxSpectraLV[iUnd][ie]
			if(MaxFluxLV < CurFluxLV) 
				MaxFluxLV = CurFluxLV
				iUndMaxLV = iUnd
			endif
		endif
		CurPolRateCR = AuxRatesCR[iUnd][ie]
		if(CurPolRateCR >= MinPolRate)
			CurFluxCR = AuxSpectraCR[iUnd][ie]
			if(MaxFluxCR < CurFluxCR) 
				MaxFluxCR = CurFluxCR
				iUndMaxCR = iUnd
			endif
		endif
		CurPolRateCL = AuxRatesCL[iUnd][ie]
		if(CurPolRateCL >= MinPolRate)
			CurFluxCL = AuxSpectraCL[iUnd][ie]
			if(MaxFluxCL < CurFluxCL) 
				MaxFluxCL = CurFluxCL
				iUndMaxCL = iUnd
			endif
		endif
		iUnd += 1
	while(iUnd < Nu)
	
	$OutWaveNameLH[ie] = MaxFluxLH
	$OutWaveNameLV[ie] = MaxFluxLV
	$OutWaveNameCR[ie] = MaxFluxCR
	$OutWaveNameCL[ie] = MaxFluxCL
	
	iq = trunc(iUndMaxLH/Np + 1e-07)
	ip = iUndMaxLH - iq*Np
	$OutWaveNameLHPar1[ie] = pStart + pStep*ip
	if(NumDims > 1)
		$OutWaveNameLHPar2[ie] = qStart + qStep*iq
	endif
	
	iq = trunc(iUndMaxLV/Np + 1e-07)
	ip = iUndMaxLV - iq*Np
	$OutWaveNameLVPar1[ie] = pStart + pStep*ip
	if(NumDims > 1)
		$OutWaveNameLVPar2[ie] = qStart + qStep*iq
	endif
	
	iq = trunc(iUndMaxCR/Np + 1e-07)
	ip = iUndMaxCR - iq*Np
	$OutWaveNameCRPar1[ie] = pStart + pStep*ip
	if(NumDims > 1)
		$OutWaveNameCRPar2[ie] = qStart + qStep*iq
	endif
	
	iq = trunc(iUndMaxCL/Np + 1e-07)
	ip = iUndMaxCL - iq*Np
	$OutWaveNameCLPar1[ie] = pStart + pStep*ip
	if(NumDims > 1)
		$OutWaveNameCLPar2[ie] = qStart + qStep*iq
	endif
	
	ie += 1
while(ie < Ne)

if(ShowIntermGraphs == 2)
	Display $OutWaveNameLH
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Hor. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLV
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Vert. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCR
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Right Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCL
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Left Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLHPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LH"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLVPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LV"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCRPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CR"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCLPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CL"
	SrwUtiGraphAddFrameAndGrid()
	
	if(NumDims > 1)
		Display $OutWaveNameLHPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 LH"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameLVPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 LV"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameCRPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 CR"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameCLPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 CL"
		SrwUtiGraphAddFrameAndGrid()
	endif
endif

KillWaves/Z AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL, AuxPolRateLH, AuxPolRateLV, AuxPolRateCR, AuxPolRateCL, $UndName
SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//Compares two waves and chooses max. values for each point
//+++++++++++++++++++++++++++++++++++++++
function AuxUtiChooseMaxValuesInWaves(wMax, wNew, wConstr, MinConstrVal, iPar1, iPar2, wPar1, wPar2)
wave wMax, wNew, wPar1, wPar2, wConstr
variable MinConstrVal, iPar1, iPar2

variable i=0, n=dimsize(wMax, 0)
variable CurVal=0
do
	if(wConstr[i] > MinConstrVal)
		CurVal = wNew[i]
		if(wMax[i] < CurVal)
			wMax[i] = CurVal
			wPar1[i] = iPar1
			wPar2[i] = iPar2
		endif
	endif
	i += 1
while(i < n)
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates max. spectral flux of undulator vs photon energy
//with Interpolation for B components
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndOptimSpecInterpB(OutWaveName, ElecName, ObsName, UndPer, UndLen, MinPolRate, BxWaveName, BzWaveName, pStep, qStep)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "OptSpec", "SrwUtiUndOptimSpec"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string ObsName = srwUtiGetValS("SrwSmpName", "Und", "") + SrwSmpType
variable UndPer=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen=srwUtiGetValN("SrwLength", 1.6, "")
variable MinPolRate=srwUtiGetValN("MinPolRate", 0.95, "SrwUtiUndOptimSpec")
string BxWaveName = srwUtiGetValS("BxWaveName", " ", "SrwUtiUndOptimSpec")
string BzWaveName = srwUtiGetValS("BzWaveName", " ", "SrwUtiUndOptimSpec")
variable pStep = srwUtiGetValN("pStep", 0.001, "SrwUtiUndOptimSpec")
variable qStep = srwUtiGetValN("qStep", 0.001, "SrwUtiUndOptimSpec")
prompt OutWaveName, "Name of wave to store spectrum"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType,";","")
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt MinPolRate,"Minimum acceptable polarization rate"
prompt BxWaveName,"Wave of Hor. Mag. Field values", popup Wavelist("*",";","")
prompt BzWaveName,"Wave of Vert. Mag. Field values", popup Wavelist("*",";","")
prompt pStep,"1st Und. Param. Step (units of waves)"
prompt qStep,"2nd Und. Param. Step (units of waves)"
Silent 1						|	...
PauseUpdate

//============EDIT THIS============
variable Phase = 0.5*Pi //phase shift between horizontal and vertical magnetic field components
variable RegectK = 0.08 //0.1 //0.08 //K regection parameter (RegectK -> 0 saves computation time)
variable AmOfExtraHarm = 6 //number of extra higher harmonics to take into account
variable PrecPar = 1.5
variable ShowIntermGraphs = 2 // 1- No,  2- Yes
variable AuxShowGraphs = 2
//================================

srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "") 
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1], "") 
srwUtiSetValN("SrwPeriod", UndPer, "")
srwUtiSetValN("SrwLength", UndLen, "")
srwUtiSetValS("BxWaveName", BxWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("BzWaveName", BzWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValN("MinPolRate", MinPolRate, "SrwUtiUndOptimSpec")
srwUtiSetValN("pStep", pStep, "SrwUtiUndOptimSpec")
srwUtiSetValN("qStep", qStep, "SrwUtiUndOptimSpec")

SrwUtiTriggerPrint(2)

variable NumDimB = WaveDims($BzWaveName)
if(NumDimB != WaveDims($BxWaveName))
	abort "Dimensions of waves containing horizontal and vertical field components are different"
endif

variable BxIsDefined = 0
if(strlen(BxWaveName) > 0)
	BxIsDefined = 1
endif
variable BzIsDefined = 0
if(strlen(BzWaveName) > 0)
	BzIsDefined = 1
endif

variable pStart = DimOffset($BzWaveName, 0)
variable NpOrig = DimSize($BzWaveName, 0)
variable pStepOrig = DimDelta($BzWaveName, 0)
variable pEnd = pStart + (NpOrig - 1)*pStepOrig

variable NumDims = WaveDims($BzWaveName)
variable qStart = DimOffset($BzWaveName, 1)
variable NqOrig = DimSize($BzWaveName, 1)
if(NqOrig == 0)
	NqOrig = 1
endif
variable qStepOrig = DimDelta($BzWaveName, 1)
variable qEnd = qStart + (NqOrig - 1)*qStepOrig

variable Ne = srwGetSmpPhotEnNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
string eUnitsStr = "eV"

string OutWaveNameLH = OutWaveName + "LH"
string OutWaveNameLV = OutWaveName + "LV"
string OutWaveNameCR = OutWaveName + "CR"
string OutWaveNameCL = OutWaveName + "CL"

string OutWaveNameLHPar1 = OutWaveName + "LH_par1"
string OutWaveNameLVPar1 = OutWaveName + "LV_par1"
string OutWaveNameCRPar1 = OutWaveName + "CR_par1"
string OutWaveNameCLPar1 = OutWaveName + "CL_par1"

string OutWaveNameLHPar2 = OutWaveName + "LH_par2"
string OutWaveNameLVPar2 = OutWaveName + "LV_par2"
string OutWaveNameCRPar2 = OutWaveName + "CR_par2"
string OutWaveNameCLPar2 = OutWaveName + "CL_par2"

make/O/N=(Ne) $OutWaveNameLH, $OutWaveNameLV, $OutWaveNameCR, $OutWaveNameCL
make/O/N=(Ne) $OutWaveNameLHPar1, $OutWaveNameLVPar1, $OutWaveNameCRPar1, $OutWaveNameCLPar1

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLH
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLV
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCR
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCL

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar1

$OutWaveNameLH = 0; $OutWaveNameLV = 0; $OutWaveNameCR = 0; $OutWaveNameCL = 0
$OutWaveNameLHPar1 = 0; $OutWaveNameLVPar1 = 0; $OutWaveNameCRPar1 = 0; $OutWaveNameCLPar1 = 0

if(NumDims > 1)
	make/O/N=(Ne) $OutWaveNameLHPar2, $OutWaveNameLVPar2, $OutWaveNameCRPar2, $OutWaveNameCLPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar2
	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar2
	
	$OutWaveNameLHPar2 = 0; $OutWaveNameLVPar2 = 0; $OutWaveNameCRPar2 = 0; $OutWaveNameCLPar2 = 0
endif

variable ElecEnergy = srwGetElecBeamEnergy(ElecName)

string UndName = "AuxUnd", StokesName = "AuxSto", AuxSpecName, AuxRateName

variable PolTypeLH = 1, PolTypeLV = 2, PolTypeCR = 5, PolTypeCL = 6
string AuxSufStr = ""

variable Bx, Bz, Kx, Kz, MaxHarm
variable MagHarmCount=0, PhiX, PhiZ, KeffE2, FundPhotEn, Kz1, Kx1
variable iq = 0, ip = 0
variable TestMaxHarm, IsRegected = 0

variable par1 = pStart, par2 = qStart

do
	par1 = pStart
	
	do
		Bx = srwUtiInterp2DBilin(par1, par2, $BxWaveName)
		Bz = srwUtiInterp2DBilin(par1, par2, $BzWaveName)
		
		if(Bx < 0)
			Bx = 0
		endif
		if(Bz < 0)
			Bz = 0
		endif
		
		Kx = srUtiUndK(Bx, UndPer*0.001)
		Kz = srUtiUndK(Bz, UndPer*0.001)
			
		IsRegected = 1
		if((Kx > 0) %| (Kz > 0))
			If(Kz <=  RegectK*Kx)
				IsRegected = 0
			else
				if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
					IsRegected = 0
				else
					if(Kz >= (1/RegectK)*Kx)
						IsRegected = 0
					endif
				endif
			endif
		endif
			
		if(IsRegected != 1)
		
			SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase,1,0,0)
			MaxHarm = round(eEnd/srUtiUndFundPhotEn(sqrt(Bx*Bx + Bz*Bz), UndPer*0.001, ElecEnergy, 2)) //+ AmOfExtraHarm
				
			TestMaxHarm = round(1.5*MaxHarm)
			if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
				MaxHarm += AmOfExtraHarm
			else 
				MaxHarm = TestMaxHarm
			endif
				
			SrwPerStoCreate(StokesName,ElecName,UndName + SrwUndType,ObsName,1,MaxHarm,PrecPar,PrecPar,1)
		
			if(ShowIntermGraphs == 2)
				print par1, par2, "  Bx =", Bx, "T  ", "Bz =", Bz, "T"
			endif

			AuxSufStr = "I" + num2str(PolTypeLH)
			SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,40,300,130,0,0); DoUpdate
			endif
			AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
			
			AuxSufStr = "R" + num2str(PolTypeLH)
			SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(440,40,300,130,0,0); DoUpdate
			endif
			AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType

			AuxUtiChooseMaxValuesInWaves($OutWaveNameLH, $AuxSpecName, $AuxRateName, MinPolRate, par1, par2, $OutWaveNameLHPar1, $OutWaveNameLHPar2)
		
			AuxSufStr = "I" + num2str(PolTypeLV)
			SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,180,300,130,0,0); DoUpdate
			endif
			AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
			
			AuxSufStr = "R" + num2str(PolTypeLV)
			SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(440,180,300,130,0,0); DoUpdate
			endif
			AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
	
			AuxUtiChooseMaxValuesInWaves($OutWaveNameLV, $AuxSpecName, $AuxRateName, MinPolRate, par1, par2, $OutWaveNameLVPar1, $OutWaveNameLVPar2)
			
			AuxSufStr = "I" + num2str(PolTypeCR)
			SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,320,300,130,0,0); DoUpdate
			endif
			AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
			
			AuxSufStr = "R" + num2str(PolTypeCR)
			SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(440,320,300,130,0,0); DoUpdate
			endif
			AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
			
			AuxUtiChooseMaxValuesInWaves($OutWaveNameCR, $AuxSpecName, $AuxRateName, MinPolRate, par1, par2, $OutWaveNameCRPar1, $OutWaveNameCRPar2)
	
			AuxSufStr = "I" + num2str(PolTypeCL)
			SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,460,300,130,0,0); DoUpdate
			endif
			AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
			
			AuxSufStr = "R" + num2str(PolTypeCL)
			SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
			if(ShowIntermGraphs == 2)
				SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(440,460,300,130,0,0); DoUpdate
			endif
			AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType

			AuxUtiChooseMaxValuesInWaves($OutWaveNameCL, $AuxSpecName, $AuxRateName, MinPolRate, par1, par2, $OutWaveNameCLPar1, $OutWaveNameCLPar2)
				
			if(ShowIntermGraphs == 2)
				if(AuxShowGraphs == 2)
					AuxShowGraphs = 1
				endif
			endif
		
		endif // FieldWasSet
	
		par1 += pStep
	while(par1 <= pEnd)
	
	par2 += qStep
while(par2 <= qEnd)

if(ShowIntermGraphs == 2)
	Display $OutWaveNameLH
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Hor. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLV
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Vert. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCR
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Right Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCL
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Left Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLHPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LH"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLVPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LV"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCRPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CR"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCLPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CL"
	SrwUtiGraphAddFrameAndGrid()
	
	if(NumDims > 1)
		Display $OutWaveNameLHPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 LH"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameLVPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 LV"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameCRPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 CR"
		SrwUtiGraphAddFrameAndGrid()
		
		Display $OutWaveNameCLPar2
		Label bottom SrwPLabelPhotEn
		Label left "Parameter 2 CL"
		SrwUtiGraphAddFrameAndGrid()
	endif
endif

//KillWaves/Z AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL, AuxPolRateLH, AuxPolRateLV, AuxPolRateCR, AuxPolRateCL, $UndName
KillWaves/Z $UndName, $StokesName

SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates max. spectral flux of undulator vs photon energy
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndOptimSpec1D(OutWaveName, ElecName, ObsName, UndPer, UndLen, Phase, BxWaveName, BzWaveName, MinPolRate)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "OptSpec", "SrwUtiUndOptimSpec"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string ObsName = srwUtiGetValS("SrwSmpName", "Und", "") + SrwSmpType
variable UndPer=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen=srwUtiGetValN("SrwLength", 1.6, "")
variable Phase=srwUtiGetValN("SrwPh0x", 1.6, "")
string BxWaveName = srwUtiGetValS("BxWaveName", " ", "SrwUtiUndOptimSpec")
string BzWaveName = srwUtiGetValS("BzWaveName", " ", "SrwUtiUndOptimSpec")
variable MinPolRate=srwUtiGetValN("MinPolRate", 0.95, "SrwUtiUndOptimSpec")
prompt OutWaveName, "Name of wave to store spectrum"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType,";","")
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt Phase,"Phase Shift [rad]"
prompt BxWaveName,"Wave of Hor. Mag. Field values", popup Wavelist("*",";","")
prompt BzWaveName,"Wave of Vert. Mag. Field values", popup Wavelist("*",";","")
prompt MinPolRate,"Minimum acceptable polarization rate"
Silent 1						|	...
PauseUpdate

srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "") 
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1], "") 
srwUtiSetValN("SrwPeriod", UndPer, "")
srwUtiSetValN("SrwLength", UndLen, "")
srwUtiSetValN("SrwPh0x", Phase, "")
srwUtiSetValS("BxWaveName", BxWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("BzWaveName", BzWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValN("MinPolRate", MinPolRate, "SrwUtiUndOptimSpec")

SrwUtiTriggerPrint(2)

variable NumDimB = WaveDims($BzWaveName)
if(NumDimB != WaveDims($BxWaveName))
	abort "Dimensions of waves containing horizontal and vertical field components are different"
endif
//if((NumDimB != 2) &% (NumDimB != 4))
//	abort "The waves containing horizontal and vertical field components are incorrect"
//endif

variable KsAreDefined = 0
if(NumDimB == 3)
	KsAreDefined = 1
endif

variable BxIsDefined = 0
if(strlen(BxWaveName) > 0)
	BxIsDefined = 1
endif
variable BzIsDefined = 0
if(strlen(BzWaveName) > 0)
	BzIsDefined = 1
endif

variable NumMagHarm

variable Np, pStart, pStep
variable Nq = 1, qStart = 0, qStep = 0
variable NumDims

if(KsAreDefined == 1)

	NumMagHarm = DimSize($BzWaveName, 1)
	
	if(BxIsDefined > 0)
		Np = DimSize($BxWaveName, NumDimB-1)
		pStart = DimOffset($BxWaveName, NumDimB-1)
		pStep = DimDelta($BxWaveName, NumDimB-1)
	else
		Np = DimSize($BzWaveName, NumDimB-1)
		pStart = DimOffset($BzWaveName, NumDimB-1)
		pStep = DimDelta($BzWaveName, NumDimB-1)
	endif

	if(BzIsDefined > 0)
		Nq = DimSize($BzWaveName, NumDimB-1)
		qStart = DimOffset($BzWaveName, NumDimB-1)
		qStep = DimDelta($BzWaveName, NumDimB-1)
	else
		Nq = DimSize($BxWaveName, NumDimB-1)
		qStart = DimOffset($BxWaveName, NumDimB-1)
		qStep = DimDelta($BxWaveName, NumDimB-1)
	endif
	
	NumDims = 1
	
else
	Np = DimSize($BzWaveName, 0)
	pStart = DimOffset($BzWaveName, 0)
	pStep = DimDelta($BzWaveName, 0)

	NumDims = WaveDims($BzWaveName)
	if(NumDims > 1)
		Nq = DimSize($BzWaveName, 1)
		qStart = DimOffset($BzWaveName, 1)
		qStep = DimDelta($BzWaveName, 1)
	endif
endif

variable Nu=Nq
//if(NumDims == 1)
//	if(Np < Nq)
//		Nu = Np
//	else
//		Nu = Nq
//	endif
//endif

variable Ne = srwGetSmpPhotEnNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
string eUnitsStr = "eV"

variable ElecEnergy = srwGetElecBeamEnergy(ElecName)

variable RegectK = 0.08

variable AmOfExtraHarm = 10 //6 // to make input variable ?
variable PrecPar = 1.5
variable ShowIntermGraphs = 2 // 1- No,  2- Yes
variable AuxShowGraphs = 2

string UndName = "AuxUnd", StokesName = "AuxSto", AuxSpecName, AuxRateName
make/O/N=(Nu, Ne) AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL
make/O/N=(Nu, Ne) AuxRatesLH, AuxRatesLV, AuxRatesCR, AuxRatesCL

AuxSpectraLH = 0
AuxSpectraLV = 0
AuxSpectraCR = 0
AuxSpectraCL = 0

AuxRatesLH = 0
AuxRatesLV = 0
AuxRatesCR = 0
AuxRatesCL = 0

variable PolTypeLH = 1, PolTypeLV = 2, PolTypeCR = 5, PolTypeCL = 6
string AuxSufStr = ""

variable Bx, Bz, Kx, Kz, MaxHarm
variable MagHarmCount=0, PhiX, PhiZ, KeffE2, FundPhotEn, Kz1, Kx1
variable iq = 0, iUnd = 0
variable FieldWasSet = 0, TestMaxHarm, IsRegected = 0
do
	//ip = 0
	//do
		 FieldWasSet = 0
		 
		if(KsAreDefined == 1)
		
			Kz = $BzWaveName[0][0][iq]
			Kx = $BxWaveName[0][0][iq]
			PhiX = $BxWaveName[1][0][iq]
			
			Kz1 = Kz
			Kx1 = Kx
			
			IsRegected = 1
			if((Kx > 0) %| (Kz > 0))
				IsRegected = 0
			//	If(Kz <=  RegectK*Kx)
			//		IsRegected = 0
			//	else
			//		if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
			//			IsRegected = 0
			//		else
			//			if(Kz >= (1/RegectK)*Kx)
			//				IsRegected = 0
			//			endif
			//		endif
			//	endif
			endif
			
			if(IsRegected != 1)
			
				SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase+PhiX,1,0,0)
			
				KeffE2 = Kx*Kx + Kz*Kz
			
				MagHarmCount=1
				do
					Kz = $BzWaveName[0][MagHarmCount][iq]
					if(Kz > 0)
						PhiZ = $BzWaveName[1][MagHarmCount][iq]
						SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,1,Kz,PhiZ)
						KeffE2 += Kz*Kz
					endif
					Kx = $BxWaveName[0][MagHarmCount][iq]
					if(Kx > 0)
						PhiX = $BxWaveName[1][MagHarmCount][iq]
						SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,2,Kx,PhiX+Phase)
						KeffE2 += Kx*Kx
					endif
					MagHarmCount += 1
				while(MagHarmCount < NumMagHarm)
			
				FundPhotEn = 950*ElecEnergy*ElecEnergy/(1 + 0.5*KeffE2)/(UndPer*0.1)
				MaxHarm = round(eEnd/FundPhotEn) //+ AmOfExtraHarm
				
				TestMaxHarm = round(1.5*MaxHarm)
				if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
					MaxHarm += AmOfExtraHarm
				else 
					MaxHarm = TestMaxHarm
				endif

				FieldWasSet = 1
			endif
			
		else
			Bx = $BxWaveName[iq]
			Bz = $BzWaveName[iq]
		
			if(Bx < 0)
				Bx = 0
			endif
			if(Bz < 0)
				Bz = 0
			endif
		
			Kx = srUtiUndK(Bx, UndPer*0.001)
			Kz = srUtiUndK(Bz, UndPer*0.001)
			
			//IsRegected = 1
			//if((Kx > 0) %| (Kz > 0))
			//	If(Kz <=  RegectK*Kx)
			//		IsRegected = 0
			//	else
			//		if((Kz <= (1+RegectK)*Kx) %& (Kz >= (1-RegectK)*Kx))
			//			IsRegected = 0
			//		else
			//			if(Kz >= (1/RegectK)*Kx)
			//				IsRegected = 0
			//			endif
			//		endif
			//	endif
			//endif
			
			//if(IsRegected != 1)
			
				SrwMagPerCreate2D(UndName,UndPer,Kz,Kx,UndLen,Phase,1,0,0)
				MaxHarm = round(eEnd/srUtiUndFundPhotEn(sqrt(Bx*Bx + Bz*Bz), UndPer*0.001, ElecEnergy, 2)) //+ AmOfExtraHarm
				
				TestMaxHarm = round(1.5*MaxHarm)
				if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
					MaxHarm += AmOfExtraHarm
				else 
					MaxHarm = TestMaxHarm
				endif
				
				 FieldWasSet = 1
			//endif
		endif
		
		if(FieldWasSet  > 0)
		
		SrwPerStoCreate(StokesName,ElecName,UndName + SrwUndType,ObsName,1,MaxHarm,PrecPar,PrecPar,1)
		
		if(ShowIntermGraphs == 2)

			if(KsAreDefined == 1)
				print iq, "  Kx1 =", Kx1, "   ", "Kz1 =", Kz1
			else
				print iq, "  Bx =", Bx, "T  ", "Bz =", Bz, "T"
			endif
			
			//if((ip == 0) %& (iq == 0))
			//	AuxShowGraphs = 2
			//else
			//	AuxShowGraphs = 1
			//endif
		endif

		//SrwSto2Int(StokesName + SrwStoType,"I",PolType,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		
		AuxSufStr = "I" + num2str(PolTypeLH)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,40,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraLH[iUnd] = $AuxSpecName[q]
		
		AuxSufStr = "I" + num2str(PolTypeLV)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,180,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraLV[iUnd] = $AuxSpecName[q]
	
		AuxSufStr = "I" + num2str(PolTypeCR)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,320,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraCR[iUnd] = $AuxSpecName[q]
	
		AuxSufStr = "I" + num2str(PolTypeCL)
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,460,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraCL[iUnd] = $AuxSpecName[q]
		
		//SrwSto2PolRate(StokesName + SrwStoType,"R",PolType,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		
		AuxSufStr = "R" + num2str(PolTypeLH)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLH,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,40,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesLH[iUnd] = $AuxRateName[q]
	
		AuxSufStr = "R" + num2str(PolTypeLV)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeLV,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,180,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesLV[iUnd] = $AuxRateName[q]

		AuxSufStr = "R" + num2str(PolTypeCR)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCR,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,320,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesCR[iUnd] = $AuxRateName[q]

		AuxSufStr = "R" + num2str(PolTypeCL)
		SrwSto2PolRate(StokesName + SrwStoType,AuxSufStr,PolTypeCL,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(450,460,300,150,0,0); DoUpdate
		endif
		AuxRateName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxRatesCL[iUnd] = $AuxRateName[q]
		
		if(ShowIntermGraphs == 2)
			if(AuxShowGraphs == 2)
				AuxShowGraphs = 1
			endif
		endif
		
		endif // FieldWasSet
	
		iUnd += 1
		//ip += 1
	//while(ip < Np)
	iq += 1
while(iq < Nu)

//Looking for maximum for each photon energy value

string OutWaveNameLH = OutWaveName + "LH"
string OutWaveNameLV = OutWaveName + "LV"
string OutWaveNameCR = OutWaveName + "CR"
string OutWaveNameCL = OutWaveName + "CL"

string OutWaveNameLHPar1 = OutWaveName + "LH_par1"
string OutWaveNameLVPar1 = OutWaveName + "LV_par1"
string OutWaveNameCRPar1 = OutWaveName + "CR_par1"
string OutWaveNameCLPar1 = OutWaveName + "CL_par1"

string OutWaveNameLHPar2 = OutWaveName + "LH_par2"
string OutWaveNameLVPar2 = OutWaveName + "LV_par2"
string OutWaveNameCRPar2 = OutWaveName + "CR_par2"
string OutWaveNameCLPar2 = OutWaveName + "CL_par2"

make/O/N=(Ne) $OutWaveNameLH, $OutWaveNameLV, $OutWaveNameCR, $OutWaveNameCL
make/O/N=(Ne) $OutWaveNameLHPar1, $OutWaveNameLVPar1, $OutWaveNameCRPar1, $OutWaveNameCLPar1

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLH
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLV
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCR
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCL

SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar1

//if(NumDims > 1)
//	make/O/N=(Ne) $OutWaveNameLHPar2, $OutWaveNameLVPar2, $OutWaveNameCRPar2, $OutWaveNameCLPar2
//	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar2
//	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLVPar2
//	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCRPar2
//	SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameCLPar2
//endif

variable MaxFluxLH, CurFluxLH, CurPolRateLH
variable MaxFluxLV, CurFluxLV, CurPolRateLV
variable MaxFluxCR, CurFluxCR, CurPolRateCR
variable MaxFluxCL, CurFluxCL, CurPolRateCL

variable ie=0, iUndMaxLH=-1, iUndMaxLV=-1, iUndMaxCR=-1, iUndMaxCL=-1
do 
	iUnd = 0
	MaxFluxLH = 0; MaxFluxLV = 0; MaxFluxCR = 0; MaxFluxCL = 0
	do
		CurPolRateLH = AuxRatesLH[iUnd][ie]
		if(CurPolRateLH >= MinPolRate)
			CurFluxLH = AuxSpectraLH[iUnd][ie]
			if(MaxFluxLH < CurFluxLH) 
				MaxFluxLH = CurFluxLH
				iUndMaxLH = iUnd
			endif
		endif
		CurPolRateLV = AuxRatesLV[iUnd][ie]
		if(CurPolRateLV >= MinPolRate)
			CurFluxLV = AuxSpectraLV[iUnd][ie]
			if(MaxFluxLV < CurFluxLV) 
				MaxFluxLV = CurFluxLV
				iUndMaxLV = iUnd
			endif
		endif
		CurPolRateCR = AuxRatesCR[iUnd][ie]
		if(CurPolRateCR >= MinPolRate)
			CurFluxCR = AuxSpectraCR[iUnd][ie]
			if(MaxFluxCR < CurFluxCR) 
				MaxFluxCR = CurFluxCR
				iUndMaxCR = iUnd
			endif
		endif
		CurPolRateCL = AuxRatesCL[iUnd][ie]
		if(CurPolRateCL >= MinPolRate)
			CurFluxCL = AuxSpectraCL[iUnd][ie]
			if(MaxFluxCL < CurFluxCL) 
				MaxFluxCL = CurFluxCL
				iUndMaxCL = iUnd
			endif
		endif
		iUnd += 1
	while(iUnd < Nu)
	
	$OutWaveNameLH[ie] = MaxFluxLH
	$OutWaveNameLV[ie] = MaxFluxLV
	$OutWaveNameCR[ie] = MaxFluxCR
	$OutWaveNameCL[ie] = MaxFluxCL
	
	//iq = trunc(iUndMaxLH/Np + 1e-07)
	//ip = iUndMaxLH - iq*Np
	//$OutWaveNameLHPar1[ie] = pStart + pStep*ip
	$OutWaveNameLHPar1[ie] = pStart + pStep*iUndMaxLH
	
	//if(NumDims > 1)
	//	$OutWaveNameLHPar2[ie] = qStart + qStep*iq
	//endif
	
	//iq = trunc(iUndMaxLV/Np + 1e-07)
	//ip = iUndMaxLV - iq*Np
	//$OutWaveNameLVPar1[ie] = pStart + pStep*ip
	$OutWaveNameLVPar1[ie] = pStart + pStep*iUndMaxLV
	
	//if(NumDims > 1)
	//	$OutWaveNameLVPar2[ie] = qStart + qStep*iq
	//endif
	
	//iq = trunc(iUndMaxCR/Np + 1e-07)
	//ip = iUndMaxCR - iq*Np
	//$OutWaveNameCRPar1[ie] = pStart + pStep*ip
	$OutWaveNameCRPar1[ie] = pStart + pStep*iUndMaxCR
	
	//if(NumDims > 1)
	//	$OutWaveNameCRPar2[ie] = qStart + qStep*iq
	//endif
	
	//iq = trunc(iUndMaxCL/Np + 1e-07)
	//ip = iUndMaxCL - iq*Np
	//$OutWaveNameCLPar1[ie] = pStart + pStep*ip
	$OutWaveNameCLPar1[ie] = pStart + pStep*iUndMaxCL
	
	//if(NumDims > 1)
	//	$OutWaveNameCLPar2[ie] = qStart + qStep*iq
	//endif
	
	ie += 1
while(ie < Ne)

if(ShowIntermGraphs == 2)
	Display $OutWaveNameLH
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Hor. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLV
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Vert. Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCR
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Right Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCL
	Label bottom SrwPLabelPhotEn
	Label left "Flux at Circ. Left Pol. [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLHPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LH"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameLVPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 LV"
	SrwUtiGraphAddFrameAndGrid()
	
	Display $OutWaveNameCRPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CR"
	SrwUtiGraphAddFrameAndGrid()

	Display $OutWaveNameCLPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1 CL"
	SrwUtiGraphAddFrameAndGrid()
	
	//if(NumDims > 1)
	//	Display $OutWaveNameLHPar2
	//	Label bottom SrwPLabelPhotEn
	//	Label left "Parameter 2 LH"
	//	SrwUtiGraphAddFrameAndGrid()
	//	
	//	Display $OutWaveNameLVPar2
	//	Label bottom SrwPLabelPhotEn
	//	Label left "Parameter 2 LV"
	//	SrwUtiGraphAddFrameAndGrid()
	//	
	//	Display $OutWaveNameCRPar2
	//	Label bottom SrwPLabelPhotEn
	//	Label left "Parameter 2 CR"
	//	SrwUtiGraphAddFrameAndGrid()
	//	
	//	Display $OutWaveNameCLPar2
	//	Label bottom SrwPLabelPhotEn
	//	Label left "Parameter 2 CL"
	//	SrwUtiGraphAddFrameAndGrid()
	//endif
endif

KillWaves/Z AuxSpectraLH, AuxSpectraLV, AuxSpectraCR, AuxSpectraCL, AuxPolRateLH, AuxPolRateLV, AuxPolRateCR, AuxPolRateCL, $UndName
SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//Calculates max. spectral flux of Planar Undulator vs photon energy
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndOptimSpec1DPlan(OutWaveName, ElecName, ObsName, UndPer, UndLen, BzWaveName)
string OutWaveName=srwUtiTruncString(srwUtiGetValS("OutWaveName", "OptSpec", "SrwUtiUndOptimSpec"), 30)
string ElecName=srwUtiGetValS("SrwElecName", "Elec", "") + SrwElecType
string ObsName = srwUtiGetValS("SrwSmpName", "Und", "") + SrwSmpType
variable UndPer=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen=srwUtiGetValN("SrwLength", 1.6, "")
string BzWaveName = srwUtiGetValS("BzWaveName", " ", "SrwUtiUndOptimSpec")
prompt OutWaveName, "Name of wave to store spectrum"
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt ObsName,SrwPSmpName1,popup Wavelist("*"+SrwSmpType,";","")
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt BzWaveName,"Wave of Vert. Mag. Field values", popup Wavelist("*",";","")
Silent 1						|	...
PauseUpdate

srwUtiSetValS("OutWaveName", OutWaveName, "SrwUtiUndOptimSpec")
srwUtiSetValS("SrwElecName", ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1], "") 
srwUtiSetValS("SrwSmpName", ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1], "") 
srwUtiSetValN("SrwPeriod", UndPer, "")
srwUtiSetValN("SrwLength", UndLen, "")
srwUtiSetValS("BzWaveName", BzWaveName, "SrwUtiUndOptimSpec")

SrwUtiTriggerPrint(2)

variable NumDimB = WaveDims($BzWaveName)

variable KsAreDefined = 0
if(NumDimB == 2)
	KsAreDefined = 1
endif

variable NumMagHarm = 1

variable Np, pStart, pStep
//variable Nq = 1, qStart = 0, qStep = 0
variable NumDims = 1

if(KsAreDefined == 1)

	NumMagHarm = DimSize($BzWaveName, 0)
	
	Np = DimSize($BzWaveName, 1)
	pStart = DimOffset($BzWaveName, 1)
	pStep = DimDelta($BzWaveName, 1)
else
	Np = DimSize($BzWaveName, 0)
	pStart = DimOffset($BzWaveName, 0)
	pStep = DimDelta($BzWaveName, 0)
endif

//variable Nu=Nq

variable Ne = srwGetSmpPhotEnNp(ObsName)
variable eStart = srwGetSmpPhotEnStart(ObsName)
variable eEnd = srwGetSmpPhotEnEnd(ObsName)
string eUnitsStr = "eV"

variable ElecEnergy = srwGetElecBeamEnergy(ElecName)

variable AmOfExtraHarm = 5 // to make input variable ?
variable PrecPar = 2 //1.5
variable ShowIntermGraphs = 2 // 1- No,  2- Yes
variable AuxShowGraphs = 2

string UndName = "AuxUnd", StokesName = "AuxSto", AuxSpecName
make/O/N=(Np, Ne) AuxSpectraLH
AuxSpectraLH = 0

//variable PolTypeLH = 1, PolTypeLV = 2, PolTypeCR = 5, PolTypeCL = 6
string AuxSufStr = ""

variable Bz, Kz, MaxHarm
variable MagHarmCount=0, KeffE2, FundPhotEn, Kz1
variable ip = 0, iUnd = 0
variable FieldWasSet = 0, TestMaxHarm, IsRegected = 0
do
	 FieldWasSet = 0
		 
	if(KsAreDefined == 1)
		
		Kz = $BzWaveName[0][ip]
		Kz1 = Kz
			
		//Hack: to remove!!!
		SrwMagPerCreate2D(UndName,UndPer,Kz,0,UndLen,0,1,0,0)
		//Hack: to remove!!!
		//SrwMagPerCreate2D(UndName,UndPer,Kz/sqrt(2),Kz/sqrt(2),UndLen,Pi/2,1,0,0)
		
		KeffE2 = Kz*Kz
			
		MagHarmCount=1
		do
			Kz = $BzWaveName[MagHarmCount][ip]
			if(Kz > 0)
				SrwMagPerAddHarm(UndName + "_map",MagHarmCount+1,1,Kz,0)
				KeffE2 += Kz*Kz
			endif
			MagHarmCount += 1
		while(MagHarmCount < NumMagHarm)
			
		FundPhotEn = 950*ElecEnergy*ElecEnergy/(1 + 0.5*KeffE2)/(UndPer*0.1)
		MaxHarm = round(eEnd/FundPhotEn) //+ AmOfExtraHarm
		TestMaxHarm = round(1.5*MaxHarm)
		if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
			MaxHarm += AmOfExtraHarm
		else 
			MaxHarm = TestMaxHarm
		endif
		FieldWasSet = 1
	else
		Bz = $BzWaveName[ip]
		if(Bz < 0)
			Bz = 0
		endif
		
		Kz = srUtiUndK(Bz, UndPer*0.001)
			
		SrwMagPerCreate2D(UndName,UndPer,Kz,0,UndLen,0,1,0,0)
		
		MaxHarm = round(eEnd/srUtiUndFundPhotEn(Bz, UndPer*0.001, ElecEnergy, 2)) //+ AmOfExtraHarm
		TestMaxHarm = round(1.5*MaxHarm)
		if(TestMaxHarm > (AmOfExtraHarm + MaxHarm))
			MaxHarm += AmOfExtraHarm
		else 
			MaxHarm = TestMaxHarm
		endif
		 FieldWasSet = 1
	endif
		
	if(FieldWasSet  > 0)
		
		SrwPerStoCreate(StokesName,ElecName,UndName + SrwUndType,ObsName,1,MaxHarm,PrecPar,PrecPar,1)
		
		if(ShowIntermGraphs == 2)
			if(KsAreDefined == 1)
				print ip, "Kz1 =", Kz1
			else
				print ip, "Bz =", Bz, "T"
			endif
		endif

		AuxSufStr = "I"
		SrwSto2Int(StokesName + SrwStoType,AuxSufStr,7,1,eEnd,2.5e-09,2.5e-09,AuxShowGraphs)
		if(ShowIntermGraphs == 2)
			SrwUtiGraphAddFrameAndGrid(); SrwUtiGraphWindResize(120,40,300,150,0,0); DoUpdate
		endif
		AuxSpecName = StokesName + AuxSufStr + SrwSeparator + SrwRadEType
		AuxSpectraLH[iUnd] = $AuxSpecName[q]
		
		if(ShowIntermGraphs == 2)
			if(AuxShowGraphs == 2)
				AuxShowGraphs = 1
			endif
		endif
	endif // FieldWasSet
	
	iUnd += 1
	ip += 1
while(ip < Np)

//Looking for maximum for each photon energy value

string OutWaveNameLH = OutWaveName + "LH"
string OutWaveNameLHPar1 = OutWaveName + "LH_par1"

make/O/N=(Ne) $OutWaveNameLH
make/O/N=(Ne) $OutWaveNameLHPar1
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLH
SetScale/I x eStart, eEnd, eUnitsStr, $OutWaveNameLHPar1

variable MaxFluxLH, CurFluxLH, CurPolRateLH

variable ie=0, iUndMaxLH=-1
do 
	iUnd = 0
	MaxFluxLH = 0
	do
		//CurPolRateLH = AuxRatesLH[iUnd][ie]
		//if(CurPolRateLH >= MinPolRate)
		CurFluxLH = AuxSpectraLH[iUnd][ie]
		if(MaxFluxLH < CurFluxLH) 
			MaxFluxLH = CurFluxLH
			iUndMaxLH = iUnd
		endif
		//endif

		iUnd += 1
	while(iUnd < Np)
	
	$OutWaveNameLH[ie] = MaxFluxLH
	$OutWaveNameLHPar1[ie] = pStart + pStep*iUndMaxLH
	
	ie += 1
while(ie < Ne)

if(ShowIntermGraphs == 2)
	Display $OutWaveNameLH
	Label bottom SrwPLabelPhotEn
	Label left "Total Flux [Ph/s/0.1%bw]"
	SrwUtiGraphAddFrameAndGrid()
		
	Display $OutWaveNameLHPar1
	Label bottom SrwPLabelPhotEn
	Label left "Parameter 1"
	SrwUtiGraphAddFrameAndGrid()
endif

KillWaves/Z AuxSpectraLH, $UndName
SrwUtiTriggerPrint(1)
end

//+++++++++++++++++++++++++++++++++++++++
//Estimates Power (heat load) through finite aperture
//+++++++++++++++++++++++++++++++++++++++
function srwUtiUndPowerEstim(ElecName, ObsName, UndPer, UndLen, Kz, Kx, Ph0x, PrecPar, Meth, Disp)
string ElecName //=srwUtiGetValS("SrwElecName", "Elec", "") + "_ebm"
string ObsName //=srwUtiGetValS("SrwSmpName", "Und", "") + "_obs"
variable UndPer //=srwUtiGetValN("SrwPeriod", 80, "")
variable UndLen //=srwUtiGetValN("SrwLength", 1.6, "")
variable Kz //=srwUtiGetValN("SrwKz", 1, "")
variable Kx //=srwUtiGetValN("SrwKx", 1, "")
variable Ph0x //=srwUtiGetValN("SrwPh0x", 0, "")
variable PrecPar //=srwUtiGetValN("SrwPowCompPrec", 1, "")
variable Meth //=srwUtiGetValN("SrwPowCompMeth", 1, "")
variable Disp //=srwUtiGetValN("SrwPowDispImmed", 1, "")
prompt ElecName,"Electron Beam"
prompt ObsName,"Radiation Sampling"
prompt UndPer,"Undulator Period [mm]"
prompt UndLen,"Undulator Length [m]"
prompt Kz,"Vertical Deflecting Parameter"
prompt Kx,"Horizontal Deflecting Parameter"
prompt Ph0x,"Phase Shift [rad]"
prompt PrecPar,"Precision parameter"
prompt Meth,"Computation Method",popup "Near Field;Far Field"
prompt Disp,"Display Results",popup "No;Yes"

execute/Z "SrwUtiTriggerPrint(2)"

string AuxUndName = "AuxUndPowerEstim", AuxPowDensName = "AuxPowDens"
string cmdStr
sprintf cmdStr, "SrwMagPerCreate2D(\"%s\",%g,%g,%g,%g,%g,1,0,0)", AuxUndName, UndPer, Kz, Kx, UndLen, Ph0x
execute/Z cmdStr
AuxUndName += "_map"

sprintf cmdStr, "SrwPowCreate(\"%s\",\"%s\",\"%s\",\"%s\",%g,%g,%g)", AuxPowDensName, ElecName, AuxUndName, ObsName, PrecPar, Meth, Disp
execute/Z cmdStr
AuxPowDensName += "_pow"

string StatWaveName = srUtiWfrLimits($AuxPowDensName, 0.9)
wave StatWave = $StatWaveName
variable OutPow = (StatWave[0])*(1e+06)

execute/Z "SrwUtiTriggerPrint(1)"
KillWaves/Z $AuxUndName, $AuxPowDensName

return OutPow
end

//+++++++++++++++++++++++++++++++++++++++
//Aux. function
//+++++++++++++++++++++++++++++++++++++++
function srwUtiAuxUndPowerVsPhotEn(ElecName, ObsName, UndPer, UndLen, Ph0x, WaveBx, WaveBz, WavePar1, WavePar2, PhotEn)
string ElecName, ObsName
variable UndPer, UndLen, PhotEn, Ph0x
wave  WaveBx, WaveBz, WavePar1, WavePar2

variable PrecPar = 1, Meth = 1, Disp = 1
return srwUtiUndPowerEstim(ElecName, ObsName, UndPer, UndLen, srUtiUndK(abs(WaveBz(WavePar1(PhotEn))(WavePar2(PhotEn))), UndPer*0.001), srUtiUndK(abs(WaveBx(WavePar1(PhotEn))(WavePar2(PhotEn))), UndPer*0.001), Ph0x, PrecPar, Meth, Disp)
end

//+++++++++++++++++++++++++++++++++++++++
//Sets up a Q-periodic undulator (with modulated magnetic field)
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagFldModul(BName, ModWaveName, sStart, sPer, NumPer)
string BName=SrwMagBname+SrwFieldWaveType
string ModWaveName=srwUtiGetValS("ModWaveName", "", "SrwUtiMagFldModul")
variable sStart=srwUtiGetValN("sStart", 0, "SrwUtiMagFldModul")
variable sPer=srwUtiGetValN("sPer", 0.01, "SrwUtiMagFldModul")
variable NumPer=srwUtiGetValN("NumPer", 0, "SrwUtiMagFldModul")
prompt BName,SrwPMagBname,popup Wavelist("*"+SrwFieldWaveType ,";", "")
prompt ModWaveName,"Name of the wave containing modulation values",popup Wavelist("*",";","")
prompt sStart, "Initial longitudinal position to start modulation [m]"
prompt sPer, "Period of modulation [m]"
prompt NumPer, "Number of periods (/modulation values)"

SrwMagBname=BName[0,strlen(BName)-strlen(SrwFieldWaveType)-1]
srwUtiSetValS("ModWaveName", ModWaveName, "SrwUtiMagFldModul")
srwUtiSetValN("sStart", sStart, "SrwUtiMagFldModul")
srwUtiSetValN("sPer", sPer, "SrwUtiMagFldModul")
srwUtiSetValN("NumPer", NumPer, "SrwUtiMagFldModul")

if(cmpstr(BName,"_none_")==0)
	abort "Please supply a magnetic field component wave"
endif
if(sPer<=0)
	abort "Period length should be positive"
endif
if(NumPer<=0)
	abort "Number of periods (/modulation values) should be positive"
endif
if(cmpstr(ModWaveName,"_none_")==0)
	abort "Please supply a modulation values wave"
endif
if(DimSIze($ModWaveName, 0) < NumPer)
	abort "The length of wave containing modulation values is too small"
endif

variable i = 0
variable sCurStartMod, sCurEndMod, CurModVal
do
	sCurStartMod = sStart + i*sPer
	sCurEndMod = sCurStartMod + sPer
	CurModVal = $ModWaveName[i] - 1

	$BName *= (1 + CurModVal*srwUtiStep(x - sCurStartMod)*srwUtiStep(sCurEndMod - x))
	i += 1
while(i < NumPer)
end

//+++++++++++++++++++++++++++++++++++++++
//Sets up a Q-periodic undulator with various positions of periods
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiMagSetupQPer(MagName, MagFldLeft, MagFldPer, MagFldRight, PerLen, NumPer, StartPer, PerModCoefsName, PerShiftsName)
string MagName=srwUtiGetValS("SrwMagBname", "Mag", "")+SrwFieldWaveType
string MagFldLeft=srwUtiGetValS("MagFldLeft", "", "SrwUtiMagSetupQPer")
string MagFldPer=srwUtiGetValS("MagFldPer", "", "SrwUtiMagSetupQPer")
string MagFldRight=srwUtiGetValS("MagFldRight", "", "SrwUtiMagSetupQPer")
variable PerLen=srwUtiGetValN("PerLen", 50, "SrwUtiMagSetupQPer")
variable StartPer=srwUtiGetValN("StartPer", 0, "SrwUtiMagSetupQPer")
variable NumPer=srwUtiGetValN("NumPer", 20, "SrwUtiMagSetupQPer")
string PerModCoefsName=srwUtiGetValS("PerModCoefsName", "", "SrwUtiMagSetupQPer")
string PerShiftsName=srwUtiGetValS("PerShiftsName", "", "SrwUtiMagSetupQPer")
prompt MagName,"Resulting Magnetic Field",popup Wavelist("*"+SrwFieldWaveType ,";", "")
prompt MagFldLeft,"Left Termination",popup Wavelist("*"+SrwFieldWaveType ,";", "")
prompt MagFldPer,"Field over One Q-Period",popup Wavelist("*"+SrwFieldWaveType ,";", "")
prompt MagFldRight,"Right Termination",popup Wavelist("*"+SrwFieldWaveType ,";", "")
prompt PerLen,"Period Length [mm]"
prompt NumPer,"Number of Q-Periods"
prompt StartPer,"Start Position of the Central Part [m]"
prompt PerModCoefsName,"Field Modulation Coef. Wave",popup Wavelist("*",";", "")
prompt PerShiftsName,"Q-Periods Shifts Wave [mm]",popup Wavelist("*",";", "")
Silent 1						|	...
PauseUpdate

SrwMagBname=MagName[0,strlen(MagName)-strlen(SrwFieldWaveType)-1]

srwUtiSetValS("MagFldLeft", MagFldLeft, "SrwUtiMagSetupQPer")
srwUtiSetValS("MagFldPer", MagFldPer, "SrwUtiMagSetupQPer")
srwUtiSetValS("MagFldRight", MagFldRight, "SrwUtiMagSetupQPer")
srwUtiSetValN("PerLen", PerLen, "SrwUtiMagSetupQPer")
srwUtiSetValN("NumPer", NumPer, "SrwUtiMagSetupQPer")
srwUtiSetValN("StartPer", StartPer, "SrwUtiMagSetupQPer")
srwUtiSetValS("PerModCoefsName", PerModCoefsName, "SrwUtiMagSetupQPer")
srwUtiSetValS("PerShiftsName", PerShiftsName, "SrwUtiMagSetupQPer")

if((cmpstr(MagName,"_none_")==0) %| (cmpstr(MagFldLeft,"_none_")==0) %| (cmpstr(MagFldPer,"_none_")==0) %| (cmpstr(MagFldRight,"_none_")==0))
	abort "Please supply a magnetic field component wave"
endif
if(PerLen<=0)
	abort "Period length should be positive"
endif
if(NumPer<=0)
	abort "Number of periods should be positive"
endif
if(cmpstr(PerModCoefsName,"_none_")==0)
	abort "Please supply modulation coefficients wave"
endif
if(cmpstr(PerShiftsName,"_none_")==0)
	abort "Please supply periods shifts wave"
endif

//determine where mag. field of left termination finishes
variable Len
variable iEndFld = dimsize($MagFldLeft, 0) - 1
do
	iEndFld -= 1
while((iEndFld > 0) %& ($MagFldLeft[iEndFld] == 0))
variable xEndFld = dimoffset($MagFldLeft, 0) + dimdelta($MagFldLeft, 0)*iEndFld
variable dx = StartPer - xEndFld

$MagName = 0
$MagName += $MagFldLeft(x - dx)

//determine where mag. field of 1 period starts
variable LenMagFldPer = dimsize($MagFldPer, 0), iStartField = 0
do
	iStartField += 1
while((iStartField < LenMagFldPer) %& ($MagFldPer[iStartField] == 0))
variable xStartField = dimoffset($MagFldPer, 0) + dimdelta($MagFldPer, 0)*iStartField
dx = StartPer - xStartField

variable PerLenM = PerLen*0.001
variable i = 0
do
	$MagName += (($PerModCoefsName[i])*($MagFldPer(x - dx)))
	dx += (PerLenM + ($PerShiftsName[i])*0.001)
	i += 1
while(i < NumPer)

//determine where mag. field of right termination starts
variable LenMagFldRight = dimsize($MagFldRight, 0)
iStartField = 0
do
	iStartField += 1
while((iStartField < LenMagFldRight) %& ($MagFldRight[iStartField] == 0))
xStartField = dimoffset($MagFldRight, 0) + dimdelta($MagFldRight, 0)*iStartField
dx = dx - xStartField
$MagName += $MagFldRight(x - dx)

end