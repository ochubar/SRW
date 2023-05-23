
//+++++++++++++++++++++++++++++++++++++++
//
// Manipulations with Radiation: Extraction of data, Propagation, etc.
//
//+++++++++++++++++++++++++++++++++++++++
//
// Prepare Radiation structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPrep(ElecName,ObsName,RadName,mode)
string ElecName,Obsname,RadName
variable mode

if(strlen(ObsName) <= 0)
	ObsName = "AuxObsWfrPrep" + SrwSmpType
endif
string ObsNameCore = ""
if(exists(ObsName) == 0)
	ObsNameCore = ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
	SrwSmpCreate(ObsNameCore,0) //creating dummy observation structure
	SrwSmpScanXZE(ObsName,0,1,1,0,1,1,1,1,1)
endif

string Radx=RadName+SrwSuffixX+SrwRadElType
string Radz=RadName+SrwSuffixZ+SrwRadElType
RadName += SrwRadType

variable eSt=($ObsName[5]), eFi=($ObsName[6])
if(($ObsName[7]==1)%&(eSt==eFi)) 
	eFi = eSt*1.000001 // Can be anything
endIf
variable xSt=($ObsName[8]), xFi=($ObsName[9])
if(($ObsName[10]==1)%&(xSt==xFi)) 
	if(xSt==0)
		xFi = 1e-12
	else
		xFi = xSt*1.000001 // Can be anything
	endif
endIf
variable zSt=($ObsName[11]), zFi=($ObsName[12])
if(($ObsName[13]==1)%&(zSt==zFi)) 
	if(zSt==0)
		zFi = 1e-12
	else
		zFi = zSt*1.000001 // Can be anything
	endif
endIf

variable ne = $ObsName[7], nx = $ObsName[10], nz = $ObsName[13]
if(ne < 1)
	ne = 1
endif
if(nx < 1)
	nx = 1
endif
if(nz < 1)
	nz = 1
endif

Make/T/O/N=22 $RadName 	// Don't forget to increase when needed !!!

string PhotEnOrTimeUnit = "eV"
$RadName[10] = "0" //frequency domain by default

variable nt = 0, tSt = 0, tFi = 0
if(dimsize($ObsName, 0) > 17)
	tSt = $ObsName[15]
	tFi = $ObsName[16]
	nt = $ObsName[17]
	if((nt > 0) %& ($ObsName[14] > 0))
		PhotEnOrTimeUnit = "s"
		ne = nt
		eSt = tSt
		eFi = tFi
		$RadName[10] = "1" //time domain
	endif
endif

KillWaves/Z $Radx,$Radz

Make/C/O/N=(ne, nx, nz) $Radx,$Radz
SetScale/I x eSt, eFi, PhotEnOrTimeUnit, $Radx,$Radz
SetScale/I y xSt, xFi, "m", $Radx,$Radz
SetScale/I z zSt, zFi, "m", $Radx,$Radz

$RadName[0]=Radx 
$RadName[1]=Radz
$RadName[2]=num2str(mode)
$RadName[3]=num2str(DimDelta($Radx,0))			// Energy step in eV
$RadName[4]=num2str(DimOffset($Radx,0)) 			// Energy start in eV
$RadName[5]=num2str(DimDelta($Radx,1))			// X step in m
$RadName[6]=num2str(DimOffset($Radx,1))			// X start in m
$RadName[7]=num2str(DimDelta($Radx,2))			// Z step in m
$RadName[8]=num2str(DimOffset($Radx,2))			// Z start in m

$RadName[9]=num2str(1)  // Quality of the diffracted wave 

//setting average photon energy:
if(ne > 1)
	$RadName[11]=num2str(0.5*(eSt + eFi))
else
	$RadName[11]=num2str(eSt)
endif

string Mat=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwMat44Type
make/O/D/N=(4,5) $Mat

if(strlen(ElecName) > 0)
	sRMatDrift($Mat,srwWfrSrcLongPos(ElecName))
endif

$RadName[12]=ElecName
$RadName[13]=Mat
$RadName[14]=num2str(1)  // Quality of the electron beam convolution

string MomX=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwSuffixX+SrwMomType
//make/O/N=(11,DimSize($Radx,0)) $MomX
make/O/D/N=(11,DimSize($Radx,0)) $MomX //OC130311
$RadName[15]=MomX
string MomZ=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwSuffixZ+SrwMomType
//make/O/N=(11,DimSize($Radz,0)) $MomZ
make/O/D/N=(11,DimSize($Radz,0)) $MomZ //OC130311
$RadName[16]=MomZ

string WfrAuxData=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwWfrAuxDataType
make/O/D/N=20 $WfrAuxData
$WfrAuxData=0
$RadName[18]=WfrAuxData

SrwWfrFldUnitSet(RadName, 1) // Electric field units (0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...)
if(strlen(ObsNameCore) > 0)
	killwaves/Z $ObsName
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Sets units of Electric Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrFldUnitSet(RadName, UnitNo)
String RadName
Variable UnitNo // 0-arbitrary, 1-sqrt(Phot/s/0.1%bw/mm^2), ...
Variable UnitInd = 19

if(exists(RadName) != 1)
	Make/T/O/N=22  $RadName
endif
if(DimSize($RadName,0) <= UnitInd)
	Redimension/N = (UnitInd + 1) $RadName
endif
$RadName[UnitInd] = num2str(UnitNo)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Sets units of Electric Field
//
//+++++++++++++++++++++++++++++++++++++++
function SrwWfrFldUnitGet(RadName)
String RadName
Variable UnitInd = 19
if(exists(RadName) != 1)
	return 0
endif
wave/T w = $RadName
if(DimSize(w,0) <= UnitInd)
	return 1 // def. units
endif
String AuxStr = w[UnitInd]
if(cmpstr(AuxStr,"")==0)
	return 1 // def. units
endif
return str2num(w[UnitInd])
end

//+++++++++++++++++++++++++++++++++++++++
function srwWfrSrcLongPos(SrcName)
String SrcName
Wave SrcWave=$SrcName

String/G SrwElecType,SrwTrjType,SrwIsotrSrcType
String SrcType=SrcName[strlen(SrcName)-strlen(SrwElecType),strlen(SrcName)-1]
if(cmpstr(SrcType,SrwElecType)==0) // Electron Beam
	return SrcWave[6]
endif
if(cmpstr(SrcType,SrwTrjType)==0) // Trajectory
	return 0. // Actual Source position is set up in C++
endif
if(cmpstr(SrcType,SrwIsotrSrcType)==0) // Isotropic Source
	return SrcWave[2]
endif
// Continue for other sources...
end

//+++++++++++++++++++++++++++++++++++++++
//
// Kill Radiation structure
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwRadKill(rad)
string rad;

KillWaves/Z $($rad[0]);
KillWaves/Z $($rad[1]);
KillWaves/Z $($rad[13]);
KillWaves/Z $($rad[15]);
KillWaves/Z $($rad[16]);
KillWaves/Z $rad;
end

//+++++++++++++++++++++++++++++++++++++++
//
// Propagate Radiation through a Beamline (uses C++ routine)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropagate(rad, blname, PropMeth, DplBeforeProp, RadDplName)
string rad=SrwRadName+SrwRadType;
string blname=SrwBliLast+SrwBeamlineType;
variable PropMeth=SrwPropMeth;
variable DplBeforeProp=SrwDplBeforeProp;
string RadDplName=SrwRadName+"d";
prompt rad,SrwPRadName1,popup Wavelist("*"+SrwRadType ,";", "");
prompt blname,SrwPBliName,popup Wavelist("*"+SrwBeamlineType ,";", "");
prompt PropMeth,SrwPPropMeth,popup "Yes;No";
prompt DplBeforeProp,SrwPDplBeforeProp,popup "No;Yes";
prompt RadDplName,SrwPRadNameDpl;
Silent 1						|	Propagating the Wavefront ...
PauseUpdate;

if(cmpstr(rad,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded;
	Return;
endif
if(cmpstr(blname,"_none_")==0)
	DoAlert 0, SrwPAlertNoOptCompFound;
	Return;
endif

SrwBliLast=blname[0,strlen(blname)-strlen(SrwBeamlineType)-1];
SrwPropMeth=PropMeth;
SrwDplBeforeProp=DplBeforeProp;
SrwRadName=rad[0,strlen(rad)-strlen(SrwRadType)-1];

if(DplBeforeProp==2)
	SrwWfrDupl(rad, RadDplName);
	string RadDpl = RadDplName + SrwRadType;
	rad=RadDpl;
	SrwRadName=RadDplName;
endif

SrwRadGenTotName=rad

variable MethNo = PropMeth;
if(PropMeth == 2)
	MethNo = 0;
else
	if(exists("SrwWfrPropDefaultMethNo") != 0)
		MethNo = SrwWfrPropDefaultMethNo;
	else
		MethNo = 1;
	endif
endif

Variable AuxPropPar1 = 0.5 // Testing

srRadPropag($rad, $blname, MethNo, AuxPropPar1);
end

//+++++++++++++++++++++++++++++++++++++++
//
// Temporary utility to set default propagation method
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrSetDefaultPropMeth(MethNo)
Variable MethNo;

if(exists("SrwWfrPropDefaultMethNo") == 0)
	Variable/G SrwWfrPropDefaultMethNo = 1;
endif

SrwWfrPropDefaultMethNo = MethNo;
if(MethNo < 1)
	SrwWfrPropDefaultMethNo = 1;
endif
if(MethNo > 3)
	SrwWfrPropDefaultMethNo = 3;
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Propagate Radiation through a Beamline (uses C++ routine)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropag(rad, blname, UseResBefore, UseResAfter, PrecParam, DplBeforeProp, RadDplName)
string rad=SrwRadName+SrwRadType
string blname=SrwBliLast+SrwBeamlineType
variable UseResBefore=srwUtiGetValN("UseResBefore", 1, "SrwWfrPropag")
variable UseResAfter=srwUtiGetValN("UseResAfter", 1, "SrwWfrPropag")
variable PrecParam=srwUtiGetValN("PrecParam", 1, "SrwWfrPropag")
variable DplBeforeProp=SrwDplBeforeProp
string RadDplName=SrwRadName+"d"
prompt rad,SrwPRadName1,popup Wavelist("*"+SrwRadType ,";", "")
prompt blname,SrwPBliName,popup Wavelist("*"+SrwBeamlineType ,";", "")
prompt UseResBefore,"Auto-Resize Before Propagation?",popup "Yes;No"
prompt UseResAfter,"Auto-Resize After Propagation?",popup "Yes;No"
prompt PrecParam,"Precision Parameter"
prompt DplBeforeProp,SrwPDplBeforeProp,popup "No;Yes"
prompt RadDplName,SrwPRadNameDpl
Silent 1						|	Propagating the Wavefront ...
PauseUpdate

if(cmpstr(rad,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded
	return
endif
if(cmpstr(blname,"_none_")==0)
	DoAlert 0, SrwPAlertNoOptCompFound
	return
endif

SrwBliLast=blname[0,strlen(blname)-strlen(SrwBeamlineType)-1]
srwUtiSetValN("UseResBefore", UseResBefore, "SrwWfrPropag")
srwUtiSetValN("UseResAfter", UseResAfter, "SrwWfrPropag")
srwUtiSetValN("PrecParam", PrecParam, "SrwWfrPropag")
SrwDplBeforeProp=DplBeforeProp
SrwRadName=rad[0,strlen(rad)-strlen(SrwRadType)-1]

if(DplBeforeProp==2)
	SrwWfrDupl(rad, RadDplName)
	string RadDpl = RadDplName + SrwRadType
	rad=RadDpl
	SrwRadName=RadDplName
endif

SrwRadGenTotName=rad

if(UseResBefore > 1)
	UseResBefore = 0
endif
if(UseResAfter > 1)
	UseResAfter = 0
endif

// Precision parameters
Make/O/N=4 waveprec 
waveprec[0]=UseResBefore
waveprec[1]=UseResAfter
waveprec[2]=PrecParam

waveprec[3]=0.5 //Variable AuxPropPar1 = 0.5 // Testing

srWfrPropag($rad, $blname, waveprec)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Propagate Radiation through a Beamline (uses C++ routine)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrProp(rad, blname, UseResBefore, UseResAfter, PrecParam, AnalTreat, DplBeforeProp, RadDplName)
string rad=SrwRadName+SrwRadType
string blname=SrwBliLast+SrwBeamlineType
variable UseResBefore=srwUtiGetValN("UseResBefore", 1, "SrwWfrPropag")
variable UseResAfter=srwUtiGetValN("UseResAfter", 1, "SrwWfrPropag")
variable PrecParam=srwUtiGetValN("PrecParam", 1, "SrwWfrPropag")
variable AnalTreat=srwUtiGetValN("AnalTreat", 2, "SrwWfrProp")
variable DplBeforeProp=SrwDplBeforeProp
string RadDplName=SrwRadName+"d"
prompt rad,SrwPRadName1,popup Wavelist("*"+SrwRadType ,";", "")
prompt blname,SrwPBliName,popup Wavelist("*"+SrwBeamlineType ,";", "")
prompt UseResBefore,"Auto-Resize Before Propagation?",popup "Yes;No"
prompt UseResAfter,"Auto-Resize After Propagation?",popup "Yes;No"
prompt PrecParam,"Precision Parameter"
prompt AnalTreat,"Allow Under-Sampling Mode?",popup "Yes;No"
prompt DplBeforeProp,SrwPDplBeforeProp,popup "No;Yes"
prompt RadDplName,SrwPRadNameDpl
Silent 1						|	Propagating the Wavefront ...
PauseUpdate

if(cmpstr(rad,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded
	return
endif
if(cmpstr(blname,"_none_")==0)
	DoAlert 0, SrwPAlertNoOptCompFound
	return
endif

if((UseResBefore == 1) %| (UseResAfter == 1))
	if(AnalTreat == 1) //to remove when implemented
		abort "Auto-Resize in the Under-Sampling mode is not implemented yet."
	endif
endif

SrwBliLast=blname[0,strlen(blname)-strlen(SrwBeamlineType)-1]
srwUtiSetValN("UseResBefore", UseResBefore, "SrwWfrPropag")
srwUtiSetValN("UseResAfter", UseResAfter, "SrwWfrPropag")
srwUtiSetValN("PrecParam", PrecParam, "SrwWfrPropag")
srwUtiSetValN("AnalTreat", AnalTreat, "SrwWfrProp")
SrwDplBeforeProp=DplBeforeProp
SrwRadName=rad[0,strlen(rad)-strlen(SrwRadType)-1]

if(DplBeforeProp==2)
	SrwWfrDupl(rad, RadDplName)
	string RadDpl = RadDplName + SrwRadType
	rad=RadDpl
	SrwRadName=RadDplName
endif

SrwRadGenTotName=rad

if(UseResBefore > 1)
	UseResBefore = 0
endif
if(UseResAfter > 1)
	UseResAfter = 0
endif

// Precision parameters
Make/O/N=5 waveprec 
waveprec[0]=UseResBefore
waveprec[1]=UseResAfter
waveprec[2]=PrecParam
waveprec[3]=0.5 //Variable AuxPropPar1 = 0.5 // Testing

if(AnalTreat == 2)
	AnalTreat = 0
endif
waveprec[4]=AnalTreat //AnalTreat - 1

srWfrPropag($rad, $blname, waveprec)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Propagate Radiation through a List of Optical Elements,
// with resizing the wavefront and applying individual 
// propagation options for each Optical Element
// 
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropList(rad, nmListOptElem, dplBeforeProp, radDplName)
string rad=SrwRadName+SrwRadType
string nmListOptElem=srwUtiGetValS("nmListOptElem", "", "SrwWfrPropList")
variable DplBeforeProp=SrwDplBeforeProp
string RadDplName=SrwRadName+"d"
prompt rad,SrwPRadName1,popup Wavelist("*"+SrwRadType ,";", "")
prompt nmListOptElem, "Optical Element and Propagation Param. List",popup WaveList("*",";","TEXT:1,DIMS:2,MINROWS:1,MINCOLS:10")
prompt DplBeforeProp,SrwPDplBeforeProp,popup "No;Yes"
prompt RadDplName,SrwPRadNameDpl
Silent 1						|	Propagating the Wavefront ...
PauseUpdate

if(cmpstr(rad,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded
	return
endif
if(exists(nmListOptElem)!=1)
	DoAlert 0, "List of Optical Elements was not provided"
	return
endif

srwUtiSetValS("nmListOptElem", nmListOptElem, "SrwWfrPropList")

string RadDpl = ""
if(DplBeforeProp==2)
	SrwWfrDupl(rad, RadDplName)
	RadDpl = RadDplName + SrwRadType
	rad=RadDpl
	//SrwRadName=RadDplName
endif
SrwRadName=rad[0,strlen(rad)-strlen(SrwRadType)-1]
SrwRadGenTotName=rad

string nmOptElem
variable useResBefore, useResAfter, precParam, allowUnderSamp
variable resizeMeth, resizeScanX, resizeStepX, resizeScanZ, resizeStepZ, resizeIsReq = 0
variable typeShift = 0, newCenX = 0, newCenZ = 0, newCenE = 1
variable numOptElem = dimsize($nmListOptElem, 0)
variable numParams = dimsize($nmListOptElem, 1)
variable iOptElem = 0
do
	nmOptElem = $nmListOptElem[iOptElem][0]
	useResBefore = str2num($nmListOptElem[iOptElem][1])
	useResAfter = str2num($nmListOptElem[iOptElem][2])
	precParam = str2num($nmListOptElem[iOptElem][3])
	allowUnderSamp = str2num($nmListOptElem[iOptElem][4])

	resizeMeth = str2num($nmListOptElem[iOptElem][5])
	resizeScanX = str2num($nmListOptElem[iOptElem][6])
	resizeStepX = str2num($nmListOptElem[iOptElem][7])
	resizeScanZ = str2num($nmListOptElem[iOptElem][8])
	resizeStepZ = str2num($nmListOptElem[iOptElem][9])
	
	if(numParams > 10)
		typeShift = str2num($nmListOptElem[iOptElem][10])
		newCenX = str2num($nmListOptElem[iOptElem][11])
		newCenZ = str2num($nmListOptElem[iOptElem][12])	
	endif	
	if(numParams > 13)
		newCenE = str2num($nmListOptElem[iOptElem][13])
	endif	
	if(typeShift > 0)
		SrwWfrShiftMesh(rad, typeShift, newCenX, newCenZ, newCenE, 1, "")
	endif

	resizeIsReq = 0
	if((resizeScanX > 0) %& (resizeScanX != 1))
		resizeIsReq = 1
	endif
	if((resizeStepX > 0) %& (resizeStepX != 1))
		resizeIsReq = 1
	endif
	if((resizeScanZ > 0) %& (resizeScanZ != 1))
		resizeIsReq = 1
	endif
	if((resizeStepZ > 0) %& (resizeStepZ != 1))
		resizeIsReq = 1
	endif
	if(resizeIsReq)
		SrwWfrResize(rad, resizeMeth, resizeScanX, resizeStepX, resizeScanZ, resizeStepZ, 1, "")
	endif
	if((strlen(nmOptElem) > 0) %& exists(nmOptElem))
		SrwWfrProp(rad, nmOptElem, useResBefore, useResAfter, precParam, allowUnderSamp, 1, "")
	endif

	iOptElem += 1
while(iOptElem < numOptElem)

SrwDplBeforeProp=DplBeforeProp
end

//+++++++++++++++++++++++++++++++++++++++
//
//Resize Wavefront vs Hor. and Vert. position (interface to C++ routine)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrResize(RadName, RadResizeMeth, RadResizeScanX, RadResizeStepX, RadResizeScanZ, RadResizeStepZ, DplBeforeResize, RadDplName)
string RadName=SrwRadName+SrwRadType;
variable RadResizeMeth=SrwRadResizeMeth;
variable RadResizeScanX=SrwRadResizeScanX;
variable RadResizeStepX=SrwRadResizeStepX;
variable RadResizeScanZ=SrwRadResizeScanZ; 
variable RadResizeStepZ=SrwRadResizeStepZ;
variable DplBeforeResize=SrwDplBeforeResize;
string RadDplName=SrwRadName+"d"
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "");
prompt RadResizeMeth,SrwPRadResizeMeth, popup "Normal;Special";
prompt RadResizeScanX, SrwPRadResizeScanX;
prompt RadResizeStepX, SrwPRadResizeStepX;
prompt RadResizeScanZ, SrwPRadResizeScanZ;
prompt RadResizeStepZ, SrwPRadResizeStepZ;
prompt DplBeforeResize,SrwPDplBeforeResize, popup "No;Yes";
prompt RadDplName,SrwPRadNameDpl;
Silent 1						|	Resizing Wavefront  ...
PauseUpdate

if(cmpstr(RadName,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded;
	Return;
endif

SrwRadResizeScanX=RadResizeScanX;
SrwRadResizeStepX=RadResizeStepX;
SrwRadResizeScanZ=RadResizeScanZ;
SrwRadResizeStepZ=RadResizeStepZ;
SrwRadResizeMeth=RadResizeMeth;

SrwDplBeforeResize=DplBeforeResize;
SrwRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1];

if(DplBeforeResize==2)
	SrwWfrDupl(RadName, RadDplName);
	RadName=RadDplName+SrwRadType;
	SrwRadName=RadDplName
endif

SrwRadGenTotName=RadName

srRadResizeXZ($RadName, RadResizeScanX, RadResizeStepX, RadResizeScanZ, RadResizeStepZ, RadResizeMeth-1);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Resize Wavefront vs Photon Energy or Time (test version; written entirely in Igor)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrResizePhotEn(RadName, RadResizeRangeE, RadResizeStepE, DplBeforeResize, RadDplName)
string RadName=SrwRadName+SrwRadType
variable RadResizeRangeE=srwUtiGetValN("RadResizeRangeE", 1, "SrwWfrResizePhotEn")
variable RadResizeStepE=srwUtiGetValN("RadResizeStepE", 1, "SrwWfrResizePhotEn")
variable DplBeforeResize=srwUtiGetValN("DplBeforeResize", 1, "SrwWfrResizePhotEn")
string RadDplName=srwUtiGetValS("RadDplName", SrwRadName+"d", "")
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "")
prompt RadResizeRangeE, "Photon Energy / Time Range Resizing"
prompt RadResizeStepE, "Photon Energy / Time Resolution Resizing"
prompt DplBeforeResize,SrwPDplBeforeResize, popup "No;Yes"
prompt RadDplName,SrwPRadNameDpl
Silent 1						|	Resizing Wavefront  ...
PauseUpdate

if((strlen(RadName) <= 0) %| (cmpstr(RadName,"_none_") == 0))
	abort SrwPAlertWavefrontNeeded
endif
if((RadResizeRangeE <= 0) %| (RadResizeStepE <= 0))
	abort "Resizing factors should positive."
endif

srwUtiSetValN("RadResizeRangeE", RadResizeRangeE, "SrwWfrResizePhotEn")
srwUtiSetValN("RadResizeStepE", RadResizeStepE, "SrwWfrResizePhotEn")
srwUtiSetValN("DplBeforeResize", DplBeforeResize, "SrwWfrResizePhotEn")

SrwRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]
SrwRadGenTotName=RadName

if(DplBeforeResize==2)
	SrwWfrDupl(RadName, RadDplName)
	RadName=RadDplName+SrwRadType
	SrwRadName=RadDplName
endif

string RadNameCore = SrwRadName
string RadX=RadNameCore+SrwSuffixX+SrwRadElType
string RadZ=RadNameCore+SrwSuffixZ+SrwRadElType

variable neOld = 0, nx = 0, nz = 0
variable eStartOld = 0, eStepOld = 0
variable xStart = 0, xStep = 0, zStart = 0, zStep = 0
string sUnitsE, sUnitsX, sUnitsZ

if(exists(RadX))
	neOld = dimsize($RadX, 0); eStartOld = dimoffset($RadX, 0); eStepOld = dimdelta($RadX, 0)
	nx = dimsize($RadX, 1); xStart = dimoffset($RadX, 1); xStep = dimdelta($RadX, 1)
	nz = dimsize($RadX, 2); zStart = dimoffset($RadX, 2); zStep = dimdelta($RadX, 2)
	sUnitsE = WaveUnits($RadX, 0)
	sUnitsX = WaveUnits($RadX, 1)
	sUnitsZ = WaveUnits($RadX, 2)
else
	neOld = dimsize($RadZ, 0); eStartOld = dimoffset($RadZ, 0); eStepOld = dimdelta($RadZ, 0)
	nx = dimsize($RadZ, 1); xStart = dimoffset($RadZ, 1); xStep = dimdelta($RadZ, 1)
	nz = dimsize($RadZ, 2); zStart = dimoffset($RadZ, 2); zStep = dimdelta($RadZ, 2)
	sUnitsE = WaveUnits($RadZ, 0)
	sUnitsX = WaveUnits($RadZ, 1)
	sUnitsZ = WaveUnits($RadZ, 2)
endif
variable eRangeOld = eStepOld*neOld
variable eHalfRangeOld = 0.5*eRangeOld
variable eEndOld = eStartOld + eRangeOld
variable eCen = eStartOld + eHalfRangeOld

variable ne =  neOld, eStart = eStartOld, eStep = eStepOld
variable eHalfRange = eHalfRangeOld

if(RadResizeStepE != 1)
	eStep = eStepOld/RadResizeStepE
endif
if(RadResizeRangeE != 1)
	eHalfRange = RadResizeRangeE*eHalfRangeOld
endif

string nmAuxRadComp = "wAuxRadCompResize"
string nmAuxRadCompMom = nmAuxRadComp + "_mom"
string wfrMomX = $RadName[15]
string wfrMomZ = $RadName[16]

ne = 2*round(eHalfRange/eStep)
//eStart = eCen - eHalfRange
variable neExtraLeft = round((eHalfRange - eHalfRangeOld)/eStep)
eStart = eStartOld - neExtraLeft*eStep

if((abs(ne - neOld) >= 2) %| (abs(eStart - eStartOld) > eStepOld))

	make/C/O/N=(ne, nx, nz) $nmAuxRadComp
	SetScale/P x eStart, eStep, sUnitsE, $nmAuxRadComp
	SetScale/P y xStart, xStep, sUnitsX, $nmAuxRadComp
	SetScale/P z zStart, zStep, sUnitsZ, $nmAuxRadComp
	
	if(exists(RadX))
		$nmAuxRadComp = $RadX(x)[q][r]*srwUtiNonZeroInterval(x, eStartOld, eEndOld)
		killwaves/Z $RadX
		rename $nmAuxRadComp $RadX
	endif
	if(exists(RadZ))
		if(exists(RadX))
			duplicate/O $RadX $nmAuxRadComp		
		endif
		$nmAuxRadComp = $RadZ(x)[q][r]*srwUtiNonZeroInterval(x, eStartOld, eEndOld)
		killwaves/Z $RadZ
		rename $nmAuxRadComp $RadZ
	endif
	
	$RadName[3] = num2str(eStep)
	$RadName[4] = num2str(eStart)
	$RadName[5] = num2str(xStep)
	$RadName[6] = num2str(xStart)
	$RadName[7] = num2str(zStep)
	$RadName[8] = num2str(zStart)
	
	//Modify Stat. Moments (_mom) structures:
	if(exists(wfrMomX))
		make/O/N=(11, ne) $nmAuxRadCompMom
		SetScale/P y eStart, eStep, sUnitsE, $nmAuxRadCompMom
		$nmAuxRadCompMom = 0
		$nmAuxRadCompMom = $wfrMomX[p](y)*srwUtiNonZeroIntervB(y, eStartOld, eEndOld)
		killwaves/Z $wfrMomX
		rename $nmAuxRadCompMom $wfrMomX
	endif
	if(exists(wfrMomZ))
		make/O/N=(11, ne) $nmAuxRadCompMom
		SetScale/P y eStart, eStep, sUnitsE, $nmAuxRadCompMom
		$nmAuxRadCompMom = 0
		$nmAuxRadCompMom = $wfrMomZ[p](y)*srwUtiNonZeroIntervB(y, eStartOld, eEndOld)
		killwaves/Z $wfrMomZ
		rename $nmAuxRadCompMom $wfrMomZ
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Shift Wavefront Mesh (keeping numbers of points and step size unchanged)
//Written entirely in Igor
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrShiftMesh(RadName, typeShift, newCenX, newCenZ, newCenE, DplBefore, RadDplName)
string RadName=SrwRadName+SrwRadType
variable typeShift=srwUtiGetValN("typeShift", 4, "SrwWfrShiftMesh")
variable newCenX=srwUtiGetValN("newCenX", 0, "SrwWfrShiftMesh")
variable newCenZ=srwUtiGetValN("newCenZ", 0, "SrwWfrShiftMesh")
variable newCenE=srwUtiGetValN("newCenE", 1, "SrwWfrShiftMesh")
variable DplBefore=srwUtiGetValN("DplBefore", 1, "SrwWfrShiftMesh")
string RadDplName=SrwRadName+"D" //srwUtiGetValS("RadDplName", SrwRadName+"d", "")
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "")
prompt typeShift, "Shift Mesh vs:", popup "Photon Energy (or Time);Horizontal Pos.;Vertical Pos.;Hor. + Vert. Pos.;Phot. En. (or Time) + Hor. Pos.;Phot. En. (or Time) + Vert. Pos.;Phot. En. (or Time) + Hor. + Vert. Pos."
prompt newCenX, "New Horizontal Center Position [mm]"
prompt newCenZ, "New Vertical Center Position [mm]"
prompt newCenE, "New Cen. Photon Energy [keV] or Time [fs]"
prompt DplBefore, "Duplicate Wavefront before Shifting?", popup "No;Yes"
prompt RadDplName, SrwPRadNameDpl
Silent 1						|	Shifting Wavefront  ...
PauseUpdate

if((strlen(RadName) <= 0) %| (cmpstr(RadName,"_none_") == 0))
	abort SrwPAlertWavefrontNeeded
endif
if(typeShift <= 0)
	return
endif

srwUtiSetValN("typeShift", typeShift, "SrwWfrShiftMesh")
srwUtiSetValN("newCenX", newCenX, "SrwWfrShiftMesh")
srwUtiSetValN("newCenZ", newCenZ, "SrwWfrShiftMesh")
srwUtiSetValN("newCenE", newCenE, "SrwWfrShiftMesh")
srwUtiSetValN("DplBefore", DplBefore, "SrwWfrShiftMesh")

SrwRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]
SrwRadGenTotName=RadName

if(DplBefore==2)
	SrwWfrDupl(RadName, RadDplName)
	RadName=RadDplName+SrwRadType
	SrwRadName=RadDplName
endif

SrwRadGenTotName=RadName

string RadX = $RadName[0]
string RadZ = $RadName[1]
variable EXexists = exists(RadX), EZexists = exists(RadZ)
variable FreqTimeRep = str2num($RadName[10])

string nmWfrComp = RadX
if(!EXexists)
	nmWfrComp = RadZ
	if(!EZexists)
		abort "Incorrect or Incomplete Wavefront structure"
	endif
endif

string nmAuxWfrComp = "A" + nmWfrComp[1, strlen(nmWfrComp) - 1] //to check
duplicate/O $nmWfrComp $nmAuxWfrComp

variable ne = 0, nx = 0, nz = 0
variable eStartOld = 0, eStep = 0
variable xStartOld = 0, xStep = 0, zStartOld = 0, zStep = 0
string sUnitsE, sUnitsX, sUnitsZ

ne = dimsize($nmWfrComp, 0); eStartOld = dimoffset($nmWfrComp, 0); eStep = dimdelta($nmWfrComp, 0)
nx = dimsize($nmWfrComp, 1); xStartOld = dimoffset($nmWfrComp, 1); xStep = dimdelta($nmWfrComp, 1)
nz = dimsize($nmWfrComp, 2); zStartOld = dimoffset($nmWfrComp, 2); zStep = dimdelta($nmWfrComp, 2)
sUnitsE = WaveUnits($nmWfrComp, 0)
sUnitsX = WaveUnits($nmWfrComp, 1)
sUnitsZ = WaveUnits($nmWfrComp, 2)

variable eRange = eStep*(ne - 1) //to check
variable eEndOld = eStartOld + eRange
variable eHalfRange = 0.5*eRange
variable eCenOld = eStartOld + eHalfRange
variable eCen = newCenE*1000 //to have Photon Energy in [eV]
if(FreqTimeRep == 1) //time domain
	eCen = newCenE*1.e-15 //to have Time in [s]
endif

variable xRange = xStep*(nx - 1)
variable xEndOld = xStartOld + xRange
variable xHalfRange = 0.5*xRange
variable xCenOld = xStartOld + xHalfRange
variable xCen = newCenX*0.001 //to have it in [m]

variable zRange = zStep*(nz - 1)
variable zEndOld = zStartOld + zRange
variable zHalfRange = 0.5*zRange
variable zCenOld = zStartOld + zHalfRange
variable zCen = newCenZ*0.001 //to have it in [m]

variable eStart = eStartOld, xStart = xStartOld, zStart = zStartOld
variable neDif = 0, nxDif = 0, nzDif = 0

if(typeShift == 1) //vs Photon Energy (or Time)
	if(eStep != 0) 
		//neDif = round((eCen - eCenOld)/eStep)
		neDif = (eCen - eCenOld)/eStep
	endif
	eStart = eStartOld + neDif*eStep
	if(eStart != eStartOld)
		SetScale/P x eStart, eStep, sUnitsE, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX(x)[q][r]*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ(x)[q][r]*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[4] = num2str(eStart)
	endif
endif
if(typeShift == 2) //vs Horizontal Pos.
	if(xStep != 0) 
		//nxDif = round((xCen - xCenOld)/xStep)
		nxDif = (xCen - xCenOld)/xStep
	endif
	xStart = xStartOld + nxDif*xStep
	if(xStart != xStartOld)
		SetScale/P y xStart, xStep, sUnitsX, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX[p](y)[r]*srwUtiNonZeroIntervB(y, xStartOld, xEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ[p](y)[r]*srwUtiNonZeroIntervB(y, xStartOld, xEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[6] = num2str(xStart)
	endif
endif
if(typeShift == 3) //vs Vertical Pos.
	if(zStep != 0) 
		//nzDif = round((zCen - zCenOld)/zStep)
		nzDif = (zCen - zCenOld)/zStep
	endif
	zStart = zStartOld + nzDif*zStep
	if(zStart != zStartOld)
		SetScale/P z zStart, zStep, sUnitsZ, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX[p][q](z)*srwUtiNonZeroIntervB(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ[p][q](z)*srwUtiNonZeroIntervB(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[8] = num2str(zStart)
	endif
endif
if(typeShift == 4) //vs Hor. + Vert. Pos.
	if(xStep != 0) 
		//nxDif = round((xCen - xCenOld)/xStep)
		nxDif = (xCen - xCenOld)/xStep
	endif
	if(zStep != 0) 
		//nzDif = round((zCen - zCenOld)/zStep)
		nzDif = (zCen - zCenOld)/zStep
	endif
	xStart = xStartOld + nxDif*xStep
	zStart = zStartOld + nzDif*zStep
	if((xStart != xStartOld) %| (zStart != zStartOld))
		SetScale/P y xStart, xStep, sUnitsX, $nmAuxWfrComp
		SetScale/P z zStart, zStep, sUnitsZ, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX[p](y)(z)*srwUtiNonZeroIntervB(y, xStartOld, xEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ[p](y)(z)*srwUtiNonZeroIntervB(y, xStartOld, xEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[6] = num2str(xStart)
		$RadName[8] = num2str(zStart)
	endif
endif
if(typeShift == 5) //vs Phot. En. (or Time) + Hor. Pos.
	if(eStep != 0) 
		//neDif = round((eCen - eCenOld)/eStep)
		neDif = (eCen - eCenOld)/eStep
	endif
	if(xStep != 0) 
		//nxDif = round((xCen - xCenOld)/xStep)
		nxDif = (xCen - xCenOld)/xStep
	endif
	eStart = eStartOld + neDif*eStep
	xStart = xStartOld + nxDif*xStep
	if((eStart != eStartOld) %| (xStart != xStartOld))
		SetScale/P x eStart, eStep, sUnitsE, $nmAuxWfrComp
		SetScale/P y xStart, xStep, sUnitsX, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX(x)(y)[r]*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(y, xStartOld, xEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ(x)(y)[r]*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(y, xStartOld, xEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[4] = num2str(eStart)
		$RadName[6] = num2str(xStart)
	endif
endif
if(typeShift == 6) //vs Phot. En. (or Time) + Vert. Pos.
	if(eStep != 0) 
		//neDif = round((eCen - eCenOld)/eStep)
		neDif = (eCen - eCenOld)/eStep
	endif
	if(zStep != 0) 
		//nzDif = round((zCen - zCenOld)/zStep)
		nzDif = (zCen - zCenOld)/zStep
	endif
	eStart = eStartOld + neDif*eStep
	zStart = zStartOld + nzDif*zStep
	if((eStart != eStartOld) %| (zStart != zStartOld))
		SetScale/P x eStart, eStep, sUnitsE, $nmAuxWfrComp
		SetScale/P z zStart, zStep, sUnitsZ, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX(x)[q](z)*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ(x)[q](z)*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[4] = num2str(eStart)
		$RadName[8] = num2str(zStart)
	endif
endif
if(typeShift == 7) //vs Phot. En. (or Time) + Hor. + Vert. Pos."
	if(eStep != 0) 
		//neDif = round((eCen - eCenOld)/eStep)
		neDif = (eCen - eCenOld)/eStep
	endif
	if(xStep != 0) 
		//nxDif = round((xCen - xCenOld)/xStep)
		nxDif = (xCen - xCenOld)/xStep
	endif
	if(zStep != 0) 
		//nzDif = round((zCen - zCenOld)/zStep)
		nzDif = (zCen - zCenOld)/zStep
	endif
	eStart = eStartOld + neDif*eStep
	xStart = xStartOld + nxDif*xStep
	zStart = zStartOld + nzDif*zStep
	if((eStart != eStartOld) %| (xStart != xStartOld) %| (zStart != zStartOld))
		SetScale/P x eStart, eStep, sUnitsE, $nmAuxWfrComp
		SetScale/P y xStart, xStep, sUnitsX, $nmAuxWfrComp
		SetScale/P z zStart, zStep, sUnitsZ, $nmAuxWfrComp
		if(EXexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadX(x)(y)(z)*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(y, xStartOld, xEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadX
		endif
		if(EZexists)
			$nmAuxWfrComp = 0
			$nmAuxWfrComp = $RadZ(x)(y)(z)*srwUtiNonZeroIntervB(x, eStartOld, eEndOld)*srwUtiNonZeroInterval(y, xStartOld, xEndOld)*srwUtiNonZeroInterval(z, zStartOld, zEndOld)
			duplicate/O $nmAuxWfrComp $RadZ
		endif
		$RadName[4] = num2str(eStart)
		$RadName[6] = num2str(xStart)
		$RadName[8] = num2str(zStart)
	endif
endif
killwaves/Z $nmAuxWfrComp
end

//+++++++++++++++++++++++++++++++++++++++
//
//Duplicate Radiation structure (general)
//Supports Electric Field (Wavefront), Stokes, Power Density
//+++++++++++++++++++++++++++++++++++++++
proc SrwRadDupl(RadIniName,RadName)
String RadIniName=SrwRadGenTotName
String RadName=SrwRadGenTotName[0,strlen(SrwRadGenTotName)-strlen(SrwRadType)-1]+"d"
prompt RadIniName, SrwPRadGenName, popup Wavelist("*"+SrwRadType,";","")+Wavelist("*"+SrwStoType,";","")+Wavelist("*"+SrwPowType,";","")
prompt RadName, SrwPRadGenNameDpl
Silent 1						|	Duplicating the Radiation structure  ...
PauseUpdate

if(cmpstr(RadIniName,"_none_")==0)
	DoAlert 0, SrwPAlertRadiationNeeded
	Return
endif

String ChosenRadType=SrwRadGenTotName[strlen(SrwRadGenTotName)-strlen(SrwRadType),strlen(SrwRadGenTotName)-1]
if(cmpstr(ChosenRadType,SrwRadType)==0)
	SrwWfrDupl(RadIniName,RadName) // Electric Field
endif
if(cmpstr(ChosenRadType,SrwStoType)==0)
	SrwStoDupl(RadIniName,RadName) // Stokes
endif
if(cmpstr(ChosenRadType,SrwPowType)==0)
	SrwPowDupl(RadIniName,RadName) // Power Density
endif

SrwRadGenTotName=RadName+ChosenRadType

end

//+++++++++++++++++++++++++++++++++++++++
//
//Duplicate Wavefront (electric field)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrDupl(RadiniName,RadName)
String RadiniName=SrwRadName+SrwRadType
String RadName=SrwRadName+"d"
prompt RadiniName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "")
prompt RadName, SrwPRadNameDpl
Silent 1						|	Duplicating the Radiation  structure  ...
PauseUpdate

if(cmpstr(RadiniName,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded;
	Return;
endif

SrwRadName=RadName

string Radx=RadName+SrwSuffixX+SrwRadElType
string Radz=RadName+SrwSuffixZ+SrwRadElType
RadName += SrwRadType

duplicate/O  $RadiniName  $RadName
$RadName[0]=Radx
$RadName[1]=Radz
duplicate/O $($RadiniName[0]) $Radx
duplicate/O $($RadiniName[1]) $Radz

string MatIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwMat44Type
string Mat=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwMat44Type
duplicate/O $MatIni $Mat
$RadName[13]=Mat

string MomXIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwSuffixX+SrwMomType
string MomX=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwSuffixX+SrwMomType
duplicate/O $MomXIni $MomX
$RadName[15]=MomX
string MomZIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwSuffixZ+SrwMomType
string MomZ=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwSuffixZ+SrwMomType
duplicate/O $MomZIni $MomZ
$RadName[16]=MomZ

string WfrAuxDataIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwWfrAuxDataType
string WfrAuxData=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SrwWfrAuxDataType
duplicate/O $WfrAuxDataIni $WfrAuxData
$RadName[18]=WfrAuxData

// Continue Duplication here, if any new elements pointed from the Rad
end

//+++++++++++++++++++++++++++++++++++++++
//
//Delete Wavefront (electric field)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrDel(RadiniName)
string RadiniName=SrwRadName+SrwRadType
prompt RadiniName, "Wavefront structure to delete", popup Wavelist("*"+SrwRadType, ";", "")
Silent 1						|	Deleting Wavefront  structure  ...
PauseUpdate

if(cmpstr(RadiniName,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded
	return
endif

string RadXIni = $RadiniName[0], RadZIni = $RadiniName[1]
string MatIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwMat44Type
string MomXIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwSuffixX+SrwMomType
string MomZIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwSuffixZ+SrwMomType
string WfrAuxDataIni=RadiniName[0,strlen(RadiniName)-strlen(SrwRadType)-1]+SrwWfrAuxDataType

KillWaves/Z $RadXIni, $RadZIni, $MatIni, $MomXIni, $MomZIni, $WfrAuxDataIni, $RadiniName

// Continue Deleting here, if any new elements pointed from the Rad
end

//+++++++++++++++++++++++++++++++++++++++
//
//Delete Wavefront (electric field)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrDelConfirm(RadiniName)
string RadiniName=SrwRadName+SrwRadType
prompt RadiniName, "Wavefront structure to delete", popup Wavelist("*"+SrwRadType, ";", "")
Silent 1						|	Duplicating the Radiation  structure  ...
PauseUpdate

if(cmpstr(RadiniName,"_none_")==0)
	DoAlert 0, SrwPAlertWavefrontNeeded
	return
endif

string PromptStr = "Are you sure you want to delete the wavefront structure  \"" + RadiniName + "\"?"
DoAlert 1, PromptStr
if(V_Flag != 1)
	return
endif
SrwWfrDel(RadiniName)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (uses C++ routine srRadExtract)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfr2Int(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, dis)
String RadName=SrwRadName+SrwRadType;
String SuffixExtract=SrwSuffixExtract;
Variable RadCmpnType=SrwViewRadCmpnType;
Variable CmpnNo=SrwCmpnNo;
Variable PlotType=SrwViewPlotType;
Variable eVal=SrwViewE;
Variable xVal=SrwViewX;
Variable zVal=SrwViewZ;
Variable dis=1;
Variable Repr=SrwRepr;
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "");
prompt SuffixExtract, SrwPViewSuffixExtract;
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total";
prompt CmpnNo,SrwPCmpnNo,popup SrwPOPUPCmpnNo;
prompt PlotType, SrwPViewPlotType, popup SrwPOPUPViewPlotType;
prompt Repr, SrwPRepr, popup "Position;Angle";
prompt eVal, SrwPViewE;
prompt xVal, SrwPViewX;
prompt zVal, SrwPViewZ;
prompt dis,SrwPViewDisplay,popup SrwPOPUPViewDisplay;
Silent 1						|	  ...
PauseUpdate;

if(cmpstr(radname,"_none_")==0)
	//DoAlert 0, SrwPAlertNoCompResultsFound; Return;

	if(SrwStartMacrosAfterRadSmp2 == 0)
		SrwStartMacrosAfterRadSmp2 = 2; // To proceed default chain through RadSampling panel
	endif
	SrwWfrCreate();
	if(SrwStartMacrosAfterRadSmp2 > 0)
		SrwWfr2Int();
		SrwStartMacrosAfterRadSmp2 = 0;
	endif
	Return;
endif

//if((RadCmpnType==7) %& (CmpnNo>2)) // if "Total" and ("Single-e Phase" or "Single-e Re(E)")
//	CmpnNo=1; // "Single-e Intens."
//endif
variable IsPhaseOrE = 0
//if((CmpnNo==3) %| (CmpnNo==4))
if((CmpnNo==3) %| (CmpnNo==4) %| (CmpnNo==7))
	IsPhaseOrE = 1
endif
if((RadCmpnType==7) %& (IsPhaseOrE > 0)) // if "Total" and ("Single-e Phase" or "Single-e Re(E)" or "Single-e Im(E)")
	abort "To extract Phase or Electric Field, a particular Polarization component should be chosen"
endif

SrwRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1];
SrwSuffixExtract=SuffixExtract;

SrwViewRadCmpnType=RadCmpnType;
SrwCmpnNo=CmpnNo;
SrwViewPlotType=PlotType;
SrwRepr=Repr;

SrwViewE=eVal
SrwViewX=xVal
SrwViewZ=zVal
SrwRepr=Repr

// Treating "Auto"
string RadNameX=$RadName[0]
if(PlotType==8) //Auto
	PlotType=srwWfr2IntTreatAuto($RadNameX)
endif

//string ViewRadName=SrwRadName+SuffixExtract;
string ViewRootRadName = SrwRadName
string ViewRadName = ViewRootRadName + SuffixExtract
variable LenViewRadName = strlen(ViewRadName)
variable LenSuff = strlen(SuffixExtract)
if(LenViewRadName > 27)
	if(LenSuff > 10)
		abort "The suffix string is too long."
	endif
	ViewRootRadName = srwUtiTruncString(ViewRootRadName, 27 - LenSuff)
	ViewRadName = ViewRootRadName + SuffixExtract
endif

SrwViewRadName=ViewRadName

string sUnitsEorT = WaveUnits($RadNameX, 0)
if(cmpstr(sUnitsEorT, "eV") == 0)
	eVal *= 1000. //make it in [eV]
endif
if(cmpstr(sUnitsEorT, "s") == 0)
	eVal *= 1e-15 //make it in [s]
endif

xVal *= 0.001 //make it in [m]
zVal *= 0.001 //make it in [m]

// Create necessary waves
Make/O/D/N=20 wAuxRadExtract
wAuxRadExtract[0] = RadCmpnType - 1
wAuxRadExtract[1] = CmpnNo - 1
wAuxRadExtract[2] = PlotType - 1
wAuxRadExtract[3] = Repr - 1
wAuxRadExtract[10] = eVal
wAuxRadExtract[11] = xVal
wAuxRadExtract[12] = zVal

Variable FldUnit = SrwWfrFldUnitGet(RadName)

String UnitsString=" ", UnitsString1=" "
if(Repr==1)
	if((CmpnNo==1) %| (CmpnNo==2)) // Flux per unit surface
		UnitsString = SrwPUnitSpAngFluxPerUnSurf
		UnitsString1 = SrwPUnitSpAngFluxPerUnSurf1
		
		if(FldUnit==0)
			UnitsString = SrwPUnitArbitrary
			UnitsString1 = SrwPUnitArbitrary
		endif
	endif
	if(CmpnNo==3) // Phase
		UnitsString = SrwPUnitPhase
		UnitsString1 = SrwPUnitPhase
	endif
	//if(CmpnNo==4) // Electric field
	if((CmpnNo==4) %| (CmpnNo==7)) // Electric field
		UnitsString = SrwPUnitElectricField
		UnitsString1 = SrwPUnitElectricField1
		
		if(FldUnit==0)
			UnitsString = SrwPUnitArbitrary
			UnitsString1 = SrwPUnitArbitrary
		endif
	endif
	if((CmpnNo==5) %| (CmpnNo==6)) // Flux
		UnitsString = SrwPUnitSpAngFlux
		UnitsString1 = SrwPUnitSpAngFlux
		
		if(FldUnit==0)
			UnitsString = SrwPUnitArbitrary
			UnitsString1 = SrwPUnitArbitrary
		endif
	endif
endif

if(PlotType==1) //vs Energy or Time
	ViewRadName+=SrwSeparator+SrwRadEType;
	if(CmpnNo == 3)
		Make/O/D/N=(DimSize($RadNameX, 0)) $ViewRadName;
	else
		Make/O/N=(DimSize($RadNameX, 0)) $ViewRadName;
	endif
	SetScale/P x DimOffset($RadNameX, 0), DimDelta($RadNameX, 0), WaveUnits($RadNameX, 0), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
	if(dis==2)
		if(DimSize($RadNameX, 0)>1)
			Display; Append $ViewRadName;
			Label left UnitsString;
			SrwUtiGraphAddFrameAndGrid()
			
			if(cmpstr(WaveUnits($RadNameX, 0), "eV") == 0)
				Label bottom SrwPLabelPhotEn;
			endif
			if(cmpstr(WaveUnits($RadNameX, 0), "s") == 0)
				Label bottom "Time";
			endif
		else
			print $ViewRadName[0], UnitsString1;
		endif
	endif
endif
if(PlotType==2) //vs Horizontal Pos.
	ViewRadName+=SrwSeparator+SrwRadXType;
	if(CmpnNo == 3)
		Make/O/D/N=(DimSize($RadNameX, 1)) $ViewRadName;
	else
		Make/O/N=(DimSize($RadNameX, 1)) $ViewRadName;
	endif
	SetScale/P x DimOffset($RadNameX, 1), DimDelta($RadNameX, 1), WaveUnits($RadNameX, 1), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
	if(dis==2)
		if(DimSize($RadNameX, 1)>1)
			Display; 	Append $ViewRadName;
			Label bottom SrwPLabelHorPos;
			Label left UnitsString;
			SrwUtiGraphAddFrameAndGrid()
		else
			print $ViewRadName[0], UnitsString1;
		endif
	endif
endif
if(PlotType==3) //vs Vertical Pos.
	ViewRadName+=SrwSeparator+SrwRadZType;
	if(CmpnNo==3)
		Make/O/D/N=(DimSize($RadNameX, 2))  $ViewRadName;
	else
		Make/O/N=(DimSize($RadNameX, 2))  $ViewRadName;
	endif
	SetScale/P x DimOffset($RadNameX, 2), DimDelta($RadNameX, 2), WaveUnits($RadNameX, 2), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
	if(dis==2)
		if(DimSize($RadNameX, 2)>1)
			Display; 	Append $ViewRadName;
			Label bottom SrwPLabelVerPos;
			Label left UnitsString;
			SrwUtiGraphAddFrameAndGrid()
		else
			print $ViewRadName[0], UnitsString1;
		endif
	endif
endif
if(PlotType==4) //vs Hor. + Vert. Pos.
	ViewRadName+=SrwSeparator+SrwRadXType+SrwRadZType;
	if(CmpnNo==3)
		Make/O/D/N=((DimSize($RadNameX, 1)), (DimSize($RadNameX, 2))) $ViewRadName;
	else
		//variable testNx = (DimSize($RadNameX, 1)
		//variable testNy = (DimSize($RadNameX, 2)
		//Make/O/N=((DimSize($RadNameX, 1)), (DimSize($RadNameX, 2))) $ViewRadName;
		Make/O/N=((DimSize($RadNameX, 1)), (DimSize($RadNameX, 2))) AuxSpecWave
		duplicate/O AuxSpecWave $ViewRadName
		killwaves/Z AuxSpecWave
	endif
	SetScale/P x DimOffset($RadNameX, 1), DimDelta($RadNameX, 1), WaveUnits($RadNameX, 1), $ViewRadName;
	SetScale/P y DimOffset($RadNameX, 2), DimDelta($RadNameX, 2), WaveUnits($RadNameX, 2), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
		
	if(dis==2)
		if((DimSize($RadNameX, 1)>1) %| (DimSize($RadNameX, 2)>1))
			Display; 	AppendImage $ViewRadName;
			SrwImageFormat(ViewRadName);
			Label bottom SrwPLabelHorPos;
			Label left SrwPLabelVerPos;
		else
			print $ViewRadName[0][0], UnitsString1;
		endif
	endif
endif
if(PlotType==5) //vs En. + Hor. Pos.
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType;
	Make/O/N=((DimSize($RadNameX, 0)), (DimSize($RadNameX, 1))) $ViewRadName;
	SetScale/P x DimOffset($RadNameX, 0), DimDelta($RadNameX, 0), WaveUnits($RadNameX, 0), $ViewRadName;
	SetScale/P y DimOffset($RadNameX, 1), DimDelta($RadNameX, 1), WaveUnits($RadNameX, 1), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
	if(dis==2)
		if((DimSize($RadNameX, 0)>1) %| (DimSize($RadNameX, 1)>1))
			Display; AppendImage $ViewRadName;
			SrwImageFormat(ViewRadName);
			Label bottom SrwPLabelPhotEn;
			Label left SrwPLabelHorPos;
		else
			print $ViewRadName[0][0], UnitsString1;
		endif
	endif
endif
if(PlotType==6) //vs En. + Vert. Pos.
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadZType;
	Make/O/N=((DimSize($RadNameX, 0)), (DimSize($RadNameX, 2))) $ViewRadName;
	SetScale/P x DimOffset($RadNameX, 0), DimDelta($RadNameX, 0), WaveUnits($RadNameX, 0), $ViewRadName;
	SetScale/P y DimOffset($RadNameX, 2), DimDelta($RadNameX, 2), WaveUnits($RadNameX, 2), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
	if(dis==2)
		if((DimSize($RadNameX, 0)>1) %| (DimSize($RadNameX, 2)>1))
			Display; 	AppendImage $ViewRadName;
			SrwImageFormat(ViewRadName);
			Label bottom SrwPLabelPhotEn;
			Label left SrwPLabelVerPos;
		else
			print $ViewRadName[0][0], UnitsString1;
		endif
	endif
endif
if(PlotType==7) //vs En. + Hor. + Vert.
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType+SrwRadZType;
	
		//print DimSize($RadNameX, 0),DimSize($RadNameX, 1),DimSize($RadNameX, 2)
	
	Make/O/N=((DimSize($RadNameX, 0)), (DimSize($RadNameX, 1)), (DimSize($RadNameX, 2))) $ViewRadName;
	SetScale/P x DimOffset($RadNameX, 0), DimDelta($RadNameX, 0), WaveUnits($RadNameX, 0), $ViewRadName;
	SetScale/P y DimOffset($RadNameX, 1), DimDelta($RadNameX, 1), WaveUnits($RadNameX, 1), $ViewRadName;
	SetScale/P z DimOffset($RadNameX, 2), DimDelta($RadNameX, 2), WaveUnits($RadNameX, 2), $ViewRadName;
	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
endif
//if(PlotType==9) //integrate over Hor. and Vert. position, display vs Phot. En.
//	ViewRadName+=SrwSeparator+SrwRadEType;
//	if(CmpnNo == 3)
//		Make/O/D/N=(DimSize($RadNameX, 0)) $ViewRadName;
//	else
//		Make/O/N=(DimSize($RadNameX, 0)) $ViewRadName;
//	endif
//	SetScale/P x DimOffset($RadNameX, 0), DimDelta($RadNameX, 0), WaveUnits($RadNameX, 0), $ViewRadName;
//	
//	srRadExtract($RadName, wAuxRadExtract, $ViewRadName);
//	
//	if(dis==2)
//		if(DimSize($RadNameX, 0)>1)
//			Display; 	Append $ViewRadName;
//			Label bottom SrwPLabelPhotEn;
//			Label left UnitsString;
//			//SrwPUnitSpAngFlux
//		else
//			print $ViewRadName[0], UnitsString1;
//		endif
//	endif
//endif

SrwUtiDataWaveInfStore(ViewRadName, "Unit", UnitsString1)

killwaves/Z wAuxRadExtract;
end

//+++++++++++++++++++++++++++++++++++++++
function srwWfr2IntTreatAuto(WaveEx)
wave WaveEx
if(DimSize(WaveEx, 0)==1)
	if(DimSize(WaveEx, 1)==1)
		return 3
	else
		if(DimSize(WaveEx, 2)==1)
			return 2
		else
			return 4
		endif
	endif
else
	if(DimSize(WaveEx, 1)==1)
		if(DimSize(WaveEx, 2)==1)
			return 1
		else
			return 6
		endif
	else
		if(DimSize(WaveEx, 2)==1)
			return 5
		else
			return 7
		endif
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (reduced version, without Representation)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfr2Int_(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, eVal, xVal, zVal, dis)
String RadName=SrwRadName+SrwRadType;
String SuffixExtract=SrwSuffixExtract;
Variable RadCmpnType=SrwViewRadCmpnType;
Variable CmpnNo=SrwCmpnNo;
Variable PlotType=SrwViewPlotType;
Variable eVal=SrwViewE;
Variable xVal=SrwViewX;
Variable zVal=SrwViewZ;
Variable dis=1;
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "");
prompt SuffixExtract, SrwPViewSuffixExtract;
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total";
prompt CmpnNo,SrwPCmpnNo,popup "Single-e Intensity;Multi-e Intensity;Single-e Phase;Single-e Re(E);Single-e Im(E)";
prompt PlotType, SrwPViewPlotType, popup "Energy;Hor.;Vert.;Hor. & Vert.;Energy & Hor.;Energy & Vert.;Energy & Hor. & Vert.;Auto";
prompt eVal, SrwPViewE;
prompt xVal, SrwPViewX;
prompt zVal, SrwPViewZ;
prompt dis,"New Display",popup "No;Yes";
Silent 1						|	  ...
PauseUpdate;

if(cmpstr(radname,"_none_")==0)
	//DoAlert 0, SrwPAlertNoCompResultsFound; Return;

	SrwStartMacrosAfterRadSmp2 = 1; // To proceed default chain through RadSampling panel
	SrwWfrCreate_();
	if(SrwStartMacrosAfterRadSmp2 > 0)
		SrwWfr2Int_();
		SrwStartMacrosAfterRadSmp2 = 0;
	endif
	Return;
endif

Variable Repr = 1 // Coordinate representation
if(CmpnNo == 5) //Single-e Im(E)
	CmpnNo == 7
endif
SrwWfr2Int(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, dis)

end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (polarization rate, calls SrwWfr2Int)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfr2PolRate(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, dis)
String RadName=SrwRadName+SrwRadType;
String SuffixExtract=SrwSuffixExtract;
Variable RadCmpnType=SrwViewRadCmpnType;
Variable CmpnNo=SrwCmpnNo;
Variable PlotType=SrwViewPlotType;
Variable eVal=SrwViewE;
Variable xVal=SrwViewX;
Variable zVal=SrwViewZ;
Variable dis=1;
Variable Repr=SrwRepr;
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "");
prompt SuffixExtract, SrwPViewSuffixExtract;
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total";
prompt CmpnNo,SrwPCmpnNo,popup SrwPOPUPCmpnNo;
prompt PlotType, SrwPViewPlotType, popup SrwPOPUPViewPlotType;
prompt Repr, SrwPRepr, popup "Position;Angle";
prompt eVal, SrwPViewE;
prompt xVal, SrwPViewX;
prompt zVal, SrwPViewZ;
prompt dis,SrwPViewDisplay,popup SrwPOPUPViewDisplay;
Silent 1						|	  ...
PauseUpdate;

if(RadCmpnType==7)
	Abort SrwPAlertBadPolRate
endif
if((CmpnNo != 1) %& (CmpnNo != 2))
	Abort SrwPAlertBadPolRate
endif

Variable CmplmntType
if(RadCmpnType==1)
	CmplmntType=2
endif
if(RadCmpnType==2)
	CmplmntType=1
endif
if(RadCmpnType==3)
	CmplmntType=4
endif
if(RadCmpnType==4)
	CmplmntType=3
endif
if(RadCmpnType==5)
	CmplmntType=6
endif
if(RadCmpnType==6)
	CmplmntType=5
endif

String ViewRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SuffixExtract
if(PlotType==8)
	PlotType=srwWfr2IntTreatAuto($($RadName[0]))
endif

if(PlotType==1)
	ViewRadName+=SrwSeparator+SrwRadEType
endif
if(PlotType==2)
	ViewRadName+=SrwSeparator+SrwRadXType
endif
if(PlotType==3)
	ViewRadName+=SrwSeparator+SrwRadZType
endif
if(PlotType==4)
	ViewRadName+=SrwSeparator+SrwRadXType+SrwRadZType
endif
if(PlotType==5)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType
endif
if(PlotType==6)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadZType
endif
if(PlotType==7)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType+SrwRadZType
endif

SrwWfr2Int(RadName, SuffixExtract, CmplmntType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, 1)
String OrtPolWave="AuxWave"
duplicate/O $ViewRadName $OrtPolWave

Variable ResIsOneVal=0
if((DimSize($ViewRadName,0)<=1) %& (DimSize($ViewRadName,1)<=1) %& (DimSize($ViewRadName,2)<=1))
	ResIsOneVal=1
endif
Variable DisInt=dis
if((ResIsOneVal != 0) %& (dis==2))
	DisInt=1
endif

SrwWfr2Int(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, DisInt)

WaveStats/Q $ViewRadName
Variable ZeroAbsTol = V_max*(1e-04) // To prevent division by zero

$ViewRadName = Abs($ViewRadName[p][q][r])/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)
KillWaves/Z $OrtPolWave
SrwUtiDataWaveInfStore(ViewRadName, "Unit", "")

if((dis==2) %& (PlotType<=3) %& (ResIsOneVal==0))
	Label left SrwPLabelPolRate
endif
if((ResIsOneVal != 0) %& (dis==2))
	print "Polarization rate: ", $ViewRadName[0][0][0]
endif

end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (polarization rate, calls SrwWfr2Int)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfr2PolRateExt(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, RateType, eVal, xVal, zVal, dis)
String RadName=SrwRadName+SrwRadType
String SuffixExtract=SrwSuffixExtract
Variable RadCmpnType=SrwViewRadCmpnType
Variable CmpnNo=SrwCmpnNo
Variable PlotType=SrwViewPlotType
Variable eVal=SrwViewE
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
//Variable Repr=SrwRepr
Variable RateType=srwUtiGetValN("SrwPolRateType", 1, "")
prompt RadName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt CmpnNo,SrwPCmpnNo,popup SrwPOPUPCmpnNo
prompt PlotType, SrwPViewPlotType, popup SrwPOPUPViewPlotType
//prompt Repr, SrwPRepr, popup "Position;Angle"
prompt RateType, "Normalisation", popup "I1/(I1+I2);(I1-I2)/(I1+I2)"
prompt eVal, SrwPViewE
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,SrwPViewDisplay,popup SrwPOPUPViewDisplay
Silent 1						|	  ...
PauseUpdate;

if(RadCmpnType==7)
	Abort SrwPAlertBadPolRate
endif
//if((CmpnNo != 1) %& (CmpnNo != 2))
//	Abort SrwPAlertBadPolRate
//endif
if((CmpnNo == 3) %| (CmpnNo == 4))
	Abort SrwPAlertBadPolRate
endif

srwUtiSetValN("SrwPolRateType", RateType, "")
Variable Repr=1
SrwRepr=Repr

Variable CmplmntType
if(RadCmpnType==1)
	CmplmntType=2
endif
if(RadCmpnType==2)
	CmplmntType=1
endif
if(RadCmpnType==3)
	CmplmntType=4
endif
if(RadCmpnType==4)
	CmplmntType=3
endif
if(RadCmpnType==5)
	CmplmntType=6
endif
if(RadCmpnType==6)
	CmplmntType=5
endif

//String ViewRadName=RadName[0,strlen(RadName)-strlen(SrwRadType)-1]+SuffixExtract
string ViewRootRadName = RadName[0,strlen(RadName)-strlen(SrwRadType)-1]
string ViewRadName = ViewRootRadName + SuffixExtract
variable LenViewRadName = strlen(ViewRadName)
variable LenSuff = strlen(SuffixExtract)
if(LenViewRadName > 27)
	if(LenSuff > 10)
		abort "The suffix string is too long."
	endif
	ViewRootRadName = srwUtiTruncString(ViewRootRadName, 27 - LenSuff)
	ViewRadName = ViewRootRadName + SuffixExtract
endif

if(PlotType==8)
	PlotType=srwWfr2IntTreatAuto($($RadName[0]))
endif

if(PlotType==1)
	ViewRadName+=SrwSeparator+SrwRadEType
endif
if(PlotType==2)
	ViewRadName+=SrwSeparator+SrwRadXType
endif
if(PlotType==3)
	ViewRadName+=SrwSeparator+SrwRadZType
endif
if(PlotType==4)
	ViewRadName+=SrwSeparator+SrwRadXType+SrwRadZType
endif
if(PlotType==5)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType
endif
if(PlotType==6)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadZType
endif
if(PlotType==7)
	ViewRadName+=SrwSeparator+SrwRadEType+SrwRadXType+SrwRadZType
endif

SrwWfr2Int(RadName, SuffixExtract, CmplmntType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, 1)
String OrtPolWave="AuxWave"
duplicate/O $ViewRadName $OrtPolWave

Variable ResIsOneVal=0
if((DimSize($ViewRadName,0)<=1) %& (DimSize($ViewRadName,1)<=1) %& (DimSize($ViewRadName,2)<=1))
	ResIsOneVal=1
endif
Variable DisInt=dis
if((ResIsOneVal != 0) %& (dis==2))
	DisInt=1
endif

SrwWfr2Int(RadName, SuffixExtract, RadCmpnType, CmpnNo, PlotType, Repr, eVal, xVal, zVal, DisInt)

WaveStats/Q $ViewRadName
Variable ZeroAbsTol = V_max*(1e-13) //V_max*(1e-04) // To prevent division by zero

if(RateType == 1)
	//$ViewRadName = Abs($ViewRadName[p][q][r])/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)
	$ViewRadName = Abs($ViewRadName[p][q][r])/srwAuxRetNonZero(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]), ZeroAbsTol)
else
	//$ViewRadName = (Abs($ViewRadName[p][q][r]) - Abs($OrtPolWave[p][q][r]))/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)	
	$ViewRadName = (Abs($ViewRadName[p][q][r]) - Abs($OrtPolWave[p][q][r]))/srwAuxRetNonZero(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]), ZeroAbsTol)
endif
SrwUtiDataWaveInfStore(ViewRadName, "Unit", "")

KillWaves/Z $OrtPolWave

if((dis==2) %& (PlotType<=3) %& (ResIsOneVal==0))
	Label left SrwPLabelPolRate
endif
if((ResIsOneVal != 0) %& (dis==2))
	print "Polarization rate: ", $ViewRadName[0][0][0]
endif

end

//+++++++++++++++++++++++++++++++++++++++
// Integrate radiation intensity
// assumes input data either in Ph/s/0.1%bw or in Ph/s/0.1%bw/mm2 or in W/mm2
// In case of integration over photon energy, the input data should be in Ph/s/0.1%bw or in Ph/s/0.1%bw/mm2,
// then the result is in W or W/mm2
//+++++++++++++++++++++++++++++++++++++++
proc SrwRadIntensInteg(IntegName, RadName, IntegType, dis, eMin, eMax, xMin, xMax, zMin, zMax)
string IntegName=srwUtiGetValS("IntegName", SrwRadName + "d", "SrwRadIntensInteg")
string RadName=srwUtiGetValS("RadName", SrwRadName, "SrwRadIntensInteg")
variable IntegType=srwUtiGetValN("IntegType", 1, "SrwRadIntensInteg")
variable eMin=srwUtiGetValN("eMin", 1, "SrwRadIntensInteg")
variable eMax=srwUtiGetValN("eMax", 1, "SrwRadIntensInteg")
variable xMin=srwUtiGetValN("xMin", 0, "SrwRadIntensInteg")
variable xMax=srwUtiGetValN("xMax", 0, "SrwRadIntensInteg")
variable zMin=srwUtiGetValN("zMin", 0, "SrwRadIntensInteg")
variable zMax=srwUtiGetValN("zMax", 0, "SrwRadIntensInteg")
variable dis=srwUtiGetValN("dis", 1, "SrwRadIntensInteg")
//prompt RadName, "Radiation structure", popup wavelist("*"+SrwRadType, ";", "")+wavelist("*"+SrwStoType, ";", "")+wavelist("*"+SrwPowType, ";", "")+wavelist("*"+"_e", ";", "")+wavelist("*"+"_x", ";", "")+wavelist("*"+"_z", ";", "")+wavelist("*"+"_xz", ";", "")+wavelist("*"+"_ex", ";", "")+wavelist("*"+"_ez", ";", "")+wavelist("*"+"_exz", ";", "")
prompt RadName, "Radiation structure to integrate", popup wavelist("*"+"_e", ";", "")+wavelist("*"+"_x", ";", "")+wavelist("*"+"_z", ";", "")+wavelist("*"+"_xz", ";", "")+wavelist("*"+"_ex", ";", "")+wavelist("*"+"_ez", ";", "")+wavelist("*"+"_exz", ";", "")+wavelist("*",";","TEXT:0,DIMS:2")
prompt IntegName, "Name of integrated structure"
prompt IntegType, "Integrate Over", popup "Photon Energy (or Time);Horizontal Position;Vertical Position;Hor. + Vert. Pos.;Ph. En. (or Time) + Hor. Pos.;Ph. En. (or TIme) + Vert. Pos.;Ph. En.(or Time) + Hor. + Vert. Pos."
prompt eMin, "Initial Photon Energy [keV]"
prompt eMax, "Final Photon Energy [keV]"
prompt xMin, "Initial Horizontal Position [mm]"
prompt xMax, "Final Horizontal Position [mm]"
prompt zMin, "Initial Vertical Position [mm]"
prompt zMax, "Final Vertical Position [mm]"
prompt dis,SrwPViewDisplay,popup SrwPOPUPViewDisplay
Silent 1						|	  ...
PauseUpdate

srwUtiSetValS("IntegName", IntegName, "SrwRadIntensInteg")
srwUtiSetValS("RadName", RadName, "SrwRadIntensInteg")
srwUtiSetValN("IntegType", IntegType, "SrwRadIntensInteg")
srwUtiSetValN("eMin", eMin, "SrwRadIntensInteg")
srwUtiSetValN("eMax", eMax, "SrwRadIntensInteg")
srwUtiSetValN("xMin", xMin, "SrwRadIntensInteg")
srwUtiSetValN("xMax", xMax, "SrwRadIntensInteg")
srwUtiSetValN("zMin", zMin, "SrwRadIntensInteg")
srwUtiSetValN("zMax", zMax, "SrwRadIntensInteg")
srwUtiSetValN("dis", dis, "SrwRadIntensInteg")

eMin *= 1000.; eMax *= 1000.
xMin *= 0.001; xMax *= 0.001
zMin *= 0.001; zMax *= 0.001

variable LenIntegName = strlen(IntegName)
if(LenIntegName > 27)
	IntegName = srwUtiTruncString(IntegName, 27)
endif

//Validation
variable DimSizeS=0, DimSizeE=0, DimSizeX=0, DimSizeZ=0 
variable DimOffsetE=0, DimOffsetX=0, DimOffsetZ=0
variable DimDeltaE=0, DimDeltaX=0, DimDeltaZ=0
string DimUnitE="", DimUnitX="", DimUnitZ=""

string WaveNameEnd = srwUtiGetNameEnd(RadName, SrwSeparator)
if(cmpstr(WaveNameEnd, "_e") == 0)
	DimSizeE = dimsize($RadName, 0); DimOffsetE = dimoffset($RadName, 0); DimDeltaE = dimdelta($RadName, 0); DimUnitE = waveunits($RadName, 0)
	DimSizeX = 0
	DimSizeZ = 0
endif
if(cmpstr(WaveNameEnd, "_x") == 0)
	DimSizeE = 0
	DimSizeX = dimsize($RadName, 0); DimOffsetX = dimoffset($RadName, 0); DimDeltaX = dimdelta($RadName, 0); DimUnitX = waveunits($RadName, 0)
	DimSizeZ = 0
endif
if(cmpstr(WaveNameEnd, "_z") == 0)
	DimSizeE = 0
	DimSizeX = 0
	DimSizeZ = dimsize($RadName, 0); DimOffsetZ = dimoffset($RadName, 0); DimDeltaZ = dimdelta($RadName, 0); DimUnitZ = waveunits($RadName, 0)
endif
if(cmpstr(WaveNameEnd, "_xz") == 0)
	DimSizeE = 0
	DimSizeX = dimsize($RadName, 0); DimOffsetX = dimoffset($RadName, 0); DimDeltaX = dimdelta($RadName, 0); DimUnitX = waveunits($RadName, 0)
	DimSizeZ = dimsize($RadName, 1); DimOffsetZ = dimoffset($RadName, 1); DimDeltaZ = dimdelta($RadName, 1); DimUnitZ = waveunits($RadName, 1)
endif
if(cmpstr(WaveNameEnd, "_ex") == 0)
	DimSizeE = dimsize($RadName, 0); DimOffsetE = dimoffset($RadName, 0); DimDeltaE = dimdelta($RadName, 0); DimUnitE = waveunits($RadName, 0)
	DimSizeX = dimsize($RadName, 1); DimOffsetX = dimoffset($RadName, 1); DimDeltaX = dimdelta($RadName, 1); DimUnitX = waveunits($RadName, 1)
	DimSizeZ = 0
endif
if(cmpstr(WaveNameEnd, "_ez") == 0)
	DimSizeE = dimsize($RadName, 0); DimOffsetE = dimoffset($RadName, 0); DimDeltaE = dimdelta($RadName, 0); DimUnitE = waveunits($RadName, 0)
	DimSizeX = 0
	DimSizeZ = dimsize($RadName, 1); DimOffsetZ = dimoffset($RadName, 1); DimDeltaZ = dimdelta($RadName, 1); DimUnitZ = waveunits($RadName, 1)
endif
if(cmpstr(WaveNameEnd, "_exz") == 0)
	DimSizeE = dimsize($RadName, 0); DimOffsetE = dimoffset($RadName, 0); DimDeltaE = dimdelta($RadName, 0); DimUnitE = waveunits($RadName, 0)
	DimSizeX = dimsize($RadName, 1); DimOffsetX = dimoffset($RadName, 1); DimDeltaX = dimdelta($RadName, 1); DimUnitX = waveunits($RadName, 1)
	DimSizeZ = dimsize($RadName, 2); DimOffsetZ = dimoffset($RadName, 2); DimDeltaZ = dimdelta($RadName, 2); DimUnitZ = waveunits($RadName, 2)
endif

string OrigDataUnitStr = srwUtiDataWaveInfGet(RadName, "Unit")
string ResDataUnitStr = ""

string NotImplementedStr = "Sorry, this type of integration is not implemented yet."

variable IsCompatible = 0
if(IntegType == 1) //Photon Energy
	if(DimSizeE > 1)
		IsCompatible = 1
		DimSizeE = 0
		
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			if((DimSizeX > 1) %| (DimSizeZ > 1))
				ResDataUnitStr =  SrwPUnitPowDen
			else
				ResDataUnitStr =  SrwPUnitPowDen1
			endif
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitPow
		endif
		
	endif
endif
if(IntegType == 2) //Horizontal Position
	if(DimSizeX > 1)
	
		abort NotImplementedStr
		
		IsCompatible = 1
		DimSizeX = 0
		
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			ResDataUnitStr = SrwPUnitSpAngFluxPerUnLen
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitSpAngFlux
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitPowDen) == 0)
			ResDataUnitStr = SrwPUnitPowPerUnLen
		endif
		
	endif
endif
if(IntegType == 3) //Vertical Position
	if(DimSizeZ > 1)
		IsCompatible = 1
		DimSizeZ = 0
		
		abort NotImplementedStr
		
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			ResDataUnitStr = SrwPUnitSpAngFluxPerUnLen
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitSpAngFlux
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitPowDen) == 0)
			ResDataUnitStr = SrwPUnitPowPerUnLen
		endif
		
	endif
endif
if(IntegType == 4) //Hor. + Vert. Pos.
	if((DimSizeX > 1) %& (DimSizeZ > 1))
		IsCompatible = 1
		DimSizeX = 0
		DimSizeZ = 0
		
		abort NotImplementedStr
		
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			ResDataUnitStr = SrwPUnitSpAngFlux
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitSpAngFlux
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitPowDen) == 0)
			ResDataUnitStr = SrwPUnitPow
		endif
		
	endif
endif
if(IntegType == 5) //En. + Hor. Pos.
	if((DimSizeE > 1) %& (DimSizeX > 1))
		IsCompatible = 1
		DimSizeE = 0
		DimSizeX = 0
		
		abort NotImplementedStr

		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			ResDataUnitStr = SrwPUnitPowPerUnLen
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitPow
		endif
		
	endif
endif
if(IntegType == 6) //En. + Vert. Pos.
	if((DimSizeE > 1) %& (DimSizeZ > 1))
		IsCompatible = 1
		DimSizeE = 0
		DimSizeZ = 0
		
		abort NotImplementedStr
		
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
			ResDataUnitStr = SrwPUnitPowPerUnLen
		endif
		if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
			ResDataUnitStr = SrwPUnitPow
		endif
		
	endif
endif
if(IntegType == 7) //En. + Hor. + Vert. Pos.
	if((DimSizeE > 1) %& (DimSizeX > 1))
		if(DimSizeZ > 1)
			IsCompatible = 1
			DimSizeE = 0
			DimSizeX = 0
			DimSizeZ = 0
			
			abort NotImplementedStr
			
			if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFluxPerUnSurf1) == 0)
				ResDataUnitStr = SrwPUnitPow
			endif
			if(cmpstr(OrigDataUnitStr, SrwPUnitSpAngFlux) == 0)
				ResDataUnitStr = SrwPUnitPow
			endif
			
		endif
	endif
endif
if(IsCompatible == 0)
	abort "The selected data structure can not be integrated over selected argument(s)."
endif
variable LenResDataUnitStr = strlen(ResDataUnitStr)

string ResWaveNameEnd = ""
if(DimSizeE > 0)
	if(DimSizeX > 0)
		if(DimSizeZ > 0)
			ResWaveNameEnd = "exz"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeE, DimSizeX, DimSizeZ) $IntegName
			setscale/P x, DimOffsetE, DimDeltaE, DimUnitE, $IntegName
			setscale/P y, DimOffsetX, DimDeltaX, DimUnitX, $IntegName			
			setscale/P z, DimOffsetZ, DimDeltaZ, DimUnitZ, $IntegName	
		else
			ResWaveNameEnd = "ex"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeE, DimSizeX) $IntegName
			setscale/P x, DimOffsetE, DimDeltaE, DimUnitE, $IntegName
			setscale/P y, DimOffsetX, DimDeltaX, DimUnitX, $IntegName	
		endif
	else
		if(DimSizeZ > 0)
			ResWaveNameEnd = "ez"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeE, DimSizeZ) $IntegName
			setscale/P x, DimOffsetE, DimDeltaE, DimUnitE, $IntegName
			setscale/P y, DimOffsetZ, DimDeltaZ, DimUnitZ, $IntegName
		else
			ResWaveNameEnd = "e"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeE) $IntegName
			setscale/P x, DimOffsetE, DimDeltaE, DimUnitE, $IntegName
		endif
	endif
else
	if(DimSizeX > 0)
		if(DimSizeZ > 0)
			ResWaveNameEnd = "xz"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeX, DimSizeZ) $IntegName
			setscale/P x, DimOffsetX, DimDeltaX, DimUnitX, $IntegName
			setscale/P y, DimOffsetZ, DimDeltaZ, DimUnitZ, $IntegName
		else
			ResWaveNameEnd = "x"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeX) $IntegName
			setscale/P x, DimOffsetX, DimDeltaX, DimUnitX, $IntegName
		endif
	else
		if(DimSizeZ > 0)
			ResWaveNameEnd = "z"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=(DimSizeZ) $IntegName
			setscale/P x, DimOffsetZ, DimDeltaZ, DimUnitZ, $IntegName
		else
			ResWaveNameEnd = "e"
			IntegName += SrwSeparator + ResWaveNameEnd
			make/O/N=1 $IntegName
		endif
	endif
endif

make/O wIntegPar = {IntegType, eMin, eMax, xMin, xMax, zMin, zMax}
srRadIntensInteg($RadName, wIntegPar, $IntegName)
killwaves/Z wIntegPar

variable PlotE=0, PlotX=0, PlotZ=0, PlotEX=0, PlotEZ=0, PlotXZ=0, PrintVal=0
if(dis == 2)
	if(cmpstr(ResWaveNameEnd,  "e") == 0)
		if(DimSizeE > 1)
			PlotE = 1
		else
			PrintVal = 1
		endif
	endif
	if(cmpstr(ResWaveNameEnd,  "x") == 0)
		if(DimSizeX > 1)
			PlotX = 1
		else
			PrintVal = 1
		endif
	endif
	if(cmpstr(ResWaveNameEnd,  "z") == 0)
		if(DimSizeZ > 1)
			PlotZ = 1
		else
			PrintVal = 1
		endif
	endif
	if(cmpstr(ResWaveNameEnd,  "ex") == 0)
		if(DimSizeE > 1)
			if(DimSizeX > 1)
				PlotEX = 1
			else
				PlotE = 1
			endif
		else
			if(DimSizeX > 1)
				PlotX = 1
			else
				PrintVal = 1
			endif
		endif
	endif
	if(cmpstr(ResWaveNameEnd,  "ez") == 0)
		if(DimSizeE > 1)
			if(DimSizeZ > 1)
				PlotEZ = 1
			else
				PlotE = 1
			endif
		else
			if(DimSizeZ > 1)
				PlotZ = 1
			else
				PrintVal = 1
			endif
		endif
	endif
	if(cmpstr(ResWaveNameEnd,  "xz") == 0)
		if(DimSizeX > 1)
			if(DimSizeZ > 1)
				PlotXZ = 1
			else
				PlotX = 1
			endif
		else
			if(DimSizeZ > 1)
				PlotZ = 1
			else
				PrintVal = 1
			endif
		endif
	endif
	
	if(PrintVal != 0)
		print $IntegName[0], ResDataUnitStr
	endif
	if(PlotE != 0)
		Display; Append $IntegName
		Label bottom SrwPLabelPhotEn
		if(LenResDataUnitStr >= 0)
			Label left ResDataUnitStr
		endif
		SrwPlotFormat(IntegName)
	endif
	if(PlotX != 0)
		Display; Append $IntegName
		Label bottom SrwPLabelHorPos
		if(LenResDataUnitStr >= 0)
			Label left ResDataUnitStr
		endif
		SrwPlotFormat(IntegName)
	endif
	if(PlotZ != 0)
		Display; Append $IntegName
		Label bottom SrwPLabelVerPos
		if(LenResDataUnitStr >= 0)
			Label left ResDataUnitStr
		endif
		SrwPlotFormat(IntegName)
	endif
	if(PlotEX != 0)
		Display; AppendImage $IntegName
		SrwImageFormat(IntegName)
		Label bottom SrwPLabelPhotEn
		Label left SrwPLabelHorPos
	endif
	if(PlotEZ != 0)
		Display; AppendImage $IntegName
		SrwImageFormat(IntegName)
		Label bottom SrwPLabelPhotEn
		Label left SrwPLabelVerPos
	endif
	if(PlotXZ != 0)
		Display; AppendImage $IntegName
		SrwImageFormat(IntegName)
		Label bottom SrwPLabelHorPos
		Label left SrwPLabelVerPos	
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
// Integrate radiation intensity
// assumes input data in Ph/s/0.1%bw/mm2 or in W/mm2
// Aux. function, to be replaced by SrwRadIntensInteg
//+++++++++++++++++++++++++++++++++++++++
proc SrwRadIntensIntegXZ(IntegName, RadName, xMin, xMax, zMin, zMax, dis)
string IntegName=srwUtiGetValS("IntegName", SrwRadName + "d", "SrwRadIntensInteg")
string RadName=srwUtiGetValS("RadName", SrwRadName, "SrwRadIntensInteg")
variable xMin=srwUtiGetValN("xMin", 0, "SrwRadIntensInteg")
variable xMax=srwUtiGetValN("xMax", 0, "SrwRadIntensInteg")
variable zMin=srwUtiGetValN("zMin", 0, "SrwRadIntensInteg")
variable zMax=srwUtiGetValN("zMax", 0, "SrwRadIntensInteg")
variable dis=srwUtiGetValN("dis", 1, "SrwRadIntensInteg")
prompt RadName, "Radiation structure to integrate", popup wavelist("*",";","TEXT:0,DIMS:2")+wavelist("*",";","TEXT:0,DIMS:3")
prompt IntegName, "Name of integrated structure"
prompt xMin, "Initial Horizontal Position [mm]"
prompt xMax, "Final Horizontal Position [mm]"
prompt zMin, "Initial Vertical Position [mm]"
prompt zMax, "Final Vertical Position [mm]"
prompt dis,SrwPViewDisplay,popup SrwPOPUPViewDisplay
Silent 1						|	  ...
PauseUpdate

srwUtiSetValS("IntegName", IntegName, "SrwRadIntensInteg")
srwUtiSetValS("RadName", RadName, "SrwRadIntensInteg")
srwUtiSetValN("xMin", xMin, "SrwRadIntensInteg")
srwUtiSetValN("xMax", xMax, "SrwRadIntensInteg")
srwUtiSetValN("zMin", zMin, "SrwRadIntensInteg")
srwUtiSetValN("zMax", zMax, "SrwRadIntensInteg")
srwUtiSetValN("dis", dis, "SrwRadIntensInteg")

xMin *= 0.001; xMax *= 0.001
zMin *= 0.001; zMax *= 0.001

variable npFirstDim = dimsize($RadName, 0)
variable npThirdDim = dimsize($RadName, 2)
variable i
string nmAux = "wAuxRadIntensIntegXZ"
if(npThirdDim == 0)
	make/O/N=1 $IntegName
	$IntegName[0] = SrwUtiIntWave2D(RadName, xMin, xMax, zMin, zMax)*1e+06 //assuming data in Ph/mm^2 and arg in m
	if(dis == 2)
		print $IntegName[0], "Ph/s/0.1%bw (or W, depending on units of the input wave)"
	endif
else
	make/O/N=(npFirstDim) $IntegName
	SetScale/P x dimoffset($RadName, 0), dimdelta($RadName, 0), WaveUnits($RadName, 0), $IntegName
	make/O/N=(dimsize($RadName, 1),  dimsize($RadName, 2)) $nmAux
	SetScale/P x dimoffset($RadName, 1), dimdelta($RadName, 1), WaveUnits($RadName, 1), $nmAux
	SetScale/P y dimoffset($RadName, 2), dimdelta($RadName, 2), WaveUnits($RadName, 2), $nmAux

	i = 0
	do
		$nmAux = $RadName[i][p][q]
		$IntegName[i] = SrwUtiIntWave2D(nmAux, xMin, xMax, zMin, zMax)*1e+06 //assuming data in Ph/mm^2 and arg in m
		i += 1
	while(i < npFirstDim)
	if(dis == 2)
		display $IntegName
	endif
	killwaves/Z $nmAux
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Integrates Intensity calculated on 3D mesh (e, x, z) vs photon energy (e), with a given Spectral Absorption.
//Assumes Input Intensity in Ph/s/0.1%bw/mm^2; Spectral Absorption dimensionless, vs [keV];
//Returns intensity in W/mm^2.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwRadIntensIntegVsPhotEnSpec(nmResInt, nmIntensMD, nmSpecTr)
string nmResInt=srwUtiGetValS("nmResInt", "Int", "SrwIntensIntegVsPhotEnWithSpec")
string nmIntensMD=srwUtiGetValS("nmIntensMD", "IntMD", "SrwIntensIntegVsPhotEnWithSpec")
string nmSpecTr=srwUtiGetValS("nmSpecTr", "Spec", "SrwIntensIntegVsPhotEnWithSpec")
prompt nmResInt,"Name for Resulting Intensity wave"
prompt nmIntensMD,"Spectral Intensity (vs e, x, z) to treat",popup Wavelist("*",";","TEXT:0") //,DIMS:3,DIMS:2,DIMS:1")
prompt nmSpecTr,"Spectral Absorption wave",popup Wavelist("*",";","TEXT:0,DIMS:1")
Silent 1						|	Computing Power Density  ...

srwUtiSetValS("nmResInt",nmResInt,"SrwIntensIntegVsPhotEnWithSpec")
srwUtiSetValS("nmIntensMD",nmIntensMD,"SrwIntensIntegVsPhotEnWithSpec")
srwUtiSetValS("nmSpecTr",nmSpecTr,"SrwIntensIntegVsPhotEnWithSpec")

variable obsEnp = dimsize($nmIntensMD, 0)
variable obsEstart_eV = dimoffset($nmIntensMD, 0)
variable obsEstep_eV = dimdelta($nmIntensMD, 0)
variable obsEfin_eV = obsEstart_eV + obsEstep_eV*(obsEnp - 1)
variable obsXnp = dimsize($nmIntensMD, 1)
variable obsZnp = dimsize($nmIntensMD, 2)

if(obsEnp <= 1)
	abort "The number of points of the in-most dimension should be > 1" 
endif

string nmAuxSpec = "wAuxIntegVsPhotEnWithSpec"
make/O/N=(obsEnp) $nmAuxSpec
SetScale/P x obsEstart_eV, obsEstep_eV, WaveUnits($nmIntensMD, 0), $nmAuxSpec

variable unitConvCoef = 1.60219e-16 //1  Phot/s/.1%bw    correspond(s) to :   1.60219e-16  W/eV
variable ix, iz
variable auxVar = 0

if(obsZnp > 0)
	make/O/N=(obsXnp, obsZnp) $nmResInt
	SetScale/P x dimoffset($nmIntensMD, 1), dimdelta($nmIntensMD, 1), WaveUnits($nmIntensMD, 1), $nmResInt
	SetScale/P y dimoffset($nmIntensMD, 2), dimdelta($nmIntensMD, 2), WaveUnits($nmIntensMD, 2), $nmResInt
	//iz = 0
	//do
	//	ix = 0
	//	do
	//		$nmAuxSpec = $nmIntensMD[p][ix][iz]
	//		$nmAuxSpec *= $nmSpecTr(x)
	//		integrate/T $nmAuxSpec
	//		$nmResInt[ix][iz] = unitConvCoef*$nmAuxSpec[obsEnp - 1]
	//		ix += 1
	//	while(ix < obsXnp)
	//	iz += 1
	//while(iz < obsZnp)
	variable phEn = obsEstart_eV
	$nmResInt = 0.5*$nmSpecTr(phEn)*$nmIntensMD(phEn)[p][q]
	do
		auxVar = $nmSpecTr(phEn)
		$nmResInt += auxVar*$nmIntensMD(phEn)[p][q]
			
		phEn += obsEstep_eV
	while(phEn < obsEfin_eV)
	$nmResInt += 0.5*$nmSpecTr(phEn)*$nmIntensMD(phEn)[p][q]
	$nmResInt *= unitConvCoef*obsEstep_eV
else
	if(obsXnp > 0)
		make/O/N=(obsXnp) $nmResInt
		SetScale/P x dimoffset($nmIntensMD, 1), dimdelta($nmIntensMD, 1), WaveUnits($nmIntensMD, 1), $nmResInt
		ix = 0
		do
			$nmAuxSpec = $nmIntensMD[p][ix]
			$nmAuxSpec *= $nmSpecTr(x)
			integrate/T $nmAuxSpec
			$nmResInt[ix] = unitConvCoef*$nmAuxSpec[obsEnp - 1]
			ix += 1
		while(ix < obsXnp)
	else
		make/O/N=(1) $nmResInt
		$nmAuxSpec = $nmIntensMD[p]
		$nmAuxSpec *= $nmSpecTr(x)
		integrate/T $nmAuxSpec
		$nmResInt[0] = unitConvCoef*$nmAuxSpec[obsEnp - 1]
	endif
endif

KillWaves/Z $nmAuxSpec
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Calculates integrated Power Density in [W/mm^2] deposited at transverse position (x,y) 
//in material from longitudinal position zSt to zFi,
//from Spectral Intensity calculated on 3D mesh (e, x, y), 
//taking into account Spectral Attenuation Length of the material.
//Assumes Input Intensity in Ph/s/0.1%bw/mm^2, Spectral Attenuation Length in [m],
//positions in [m].
//To obtain volume power density in [W/mm^3], this output must be multiplied by 0.001/dz.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function srwRadPowDensFromStokes(xx, yy, zz, dz, nmStokes, nmSpecAttenLen)
variable xx, yy, zz, dz //z is longitudinal coord. (assuming z=0 where absorbtion starts)
string nmStokes, nmSpecAttenLen
PauseUpdate

//variable xStart = dimoffset(wSpecInt, 1), yStart = dimoffset(wSpecInt, 2)
//variable xNp = dimsize(wSpecInt, 1), yNp = dimsize(wSpecInt, 2)
//variable xEnd = xStart + (xNp - 1)*dimdelta(wSpecInt, 1), yEnd = yStart + (yNp - 1)*dimdelta(wSpecInt, 1)
//if((xx < xStart) %| (xx > xEnd))
//	return 0
//endif
//if((yy < yStart) %| (yy > yEnd))
//	return 0
//endif

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

string str2exe = "SrwSto2IntF(\"" + nmStokes + "\",\"\",7,1,1,1," + num2str(1000*xx) + "," + num2str(1000*yy) + ",1)"
execute/Q str2exe
string nmSpecInt = nmStokes[0,strlen(nmStokes)-strlen("ras")-1] + "e"
wave wSpecInt = $nmSpecInt

wave wSpecAttenLen = $nmSpecAttenLen
//wAuxSpec001 = Interp3d(wSpecInt, x, xx, yy)

wSpecInt *= (exp(-zSt/wSpecAttenLen(x)) - exp(-zFi/wSpecAttenLen(x)))*(1.60218e-16)
integrate/T wSpecInt //[W/mm^2]
//variable resPowDen = 0.001*wAuxSpec001[eNp - 1]/(zFi - zSt) //[W/mm^3]
variable resPowDen = wSpecInt[dimsize(wSpecInt, 0) - 1] //[W/mm^2]
killwaves/Z wSpecInt
return resPowDen
end

//+++++++++++++++++++++++++++++++++++++++
//
// Representation
//
// Change the representation from spatial to angular
//  mode =0  : go to spatial representation
//  mode =1  : go to angular representation
//
//+++++++++++++++++++++++++++++++++++++++
function SrwRadRep(rad,mode)
wave/T rad
Variable mode

wave/C bx=$rad[0],bz=$rad[1]
variable ne=DimSize(bx,0),nx=DimSize(bx,1),nz=DimSize(bx,2),ie, r=str2num(rad[2])
Variable wavelength,angleoffsetX,angleoffsetZ,angledeltaX,angledeltaZ

if(mode != r)

Make/C/O/N=(nx,nz) tmpx, tmpz //AuxDebugTest

if (ne>1) 
//Make/C/O/N=(nx,nz) tmpx, tmpz //AuxDebugTest
//print "OK"
endif
variable i,j

//variable t0=startmstimer

if (r==0)

if (ne>1) 
ie=ne-1
	do
	
	tmpx=bx[ie][p][q]
	tmpz=bz[ie][p][q]

	//AuxDebugTest
		SetScale/P x DimOffset(bx,1),DimDelta(bx,1), tmpx,tmpz;
		SetScale/P y DimOffset(bx,2),DimDelta(bx,2), tmpx,tmpz;
		
		//srFFT2D(tmpx, 1); srFFT2D(tmpz, 1);
		//print "srFFT2D Used"
	//End AuxDebugTest
	
	j=0
	do
	i=0
	do
	bx[ie][i][j]=tmpx[i][j]
	bz[ie][i][j]=tmpz[i][j]
	i+=1
	while(i<nx)
	j+=1
	while (j<nz)
	ie=ie-1
	while (ie>-1)
	
	//AuxDebugTest
		SetScale/P y DimOffset(tmpx,0),DimDelta(tmpx,0), bx,bz;
		SetScale/P z DimOffset(tmpx,1),DimDelta(tmpx,1), bx,bz;
	//End AuxDebugTest

else

	//AuxDebugTest
		tmpx=bx[0][p][q]
		tmpz=bz[0][p][q]
		SetScale/P x DimOffset(bx,1),DimDelta(bx,1), tmpx,tmpz;
		SetScale/P y DimOffset(bx,2),DimDelta(bx,2), tmpx,tmpz;
		
		//srFFT2D(tmpx, 1); srFFT2D(tmpz, 1);
		//print "srFFT2D Used"
		
		j=0
		do
			i=0
			do
				bx[0][i][j]=tmpx[i][j]
				bz[0][i][j]=tmpz[i][j]
				i+=1
			while(i<nx)
			j+=1
		while (j<nz)
		
		SetScale/P y DimOffset(tmpx,0),DimDelta(tmpx,0), bx,bz;
		SetScale/P z DimOffset(tmpx,1),DimDelta(tmpx,1), bx,bz;
	//End AuxDebugTest

endif  //if (ne>1)  

	rad[2]=num2str(1)

//	wavelength= 12.4/str2num(rad[4])*1e-7    // wavelength in m 
	wavelength=1; 
	
	string unit="q" //  q= teta / Lambda  [m-1]
	angledeltaX=wavelength/str2num(rad[5])/DimSize(bx,1 )			//X Angle  range in mrad
	angledeltaZ=wavelength/str2num(rad[7])/DimSize(bx,2 )			//Z Angle  range in mrad
	angleoffsetX=-angledeltaX*(DimSize(bx,1 ))/2			
	angleoffsetZ=-angledeltaZ*(DimSize(bx,2 ))/2
	
	SetScale/P x str2num(rad[4]),str2num(rad[3]),"eV", bx,bz
	//SetScale/P y angleoffsetX,angledeltaX,unit, bx,bz //AuxDebugTest
	//SetScale/P z  angleoffsetZ,angledeltaZ,unit, bx,bz //AuxDebugTest
	
else  // if (r==0)

if (ne>1) 
ie=ne-1
	do
	tmpx=bx[ie][p][q]
	tmpz=bz[ie][p][q]

	//AuxDebugTest
		SetScale/P x DimOffset(bx,1),DimDelta(bx,1), tmpx,tmpz;
		SetScale/P y DimOffset(bx,2),DimDelta(bx,2), tmpx,tmpz;

		//srFFT2D(tmpx, -1); srFFT2D(tmpz, -1);
		//print "srFFT2D Inverse Used"
	//End AuxDebugTest
	
	j=0
	do
	i=0
	do
	bx[ie][i][j]=tmpx[i][j]
	bz[ie][i][j]=tmpz[i][j]
	i+=1
	while(i<nx)
	j+=1
	while (j<nz)
	ie=ie-1							
	while (ie>-1)
	
	//AuxDebugTest
		SetScale/P y DimOffset(tmpx,0),DimDelta(tmpx,0),"m", bx,bz;
		SetScale/P z DimOffset(tmpx,1),DimDelta(tmpx,1),"m", bx,bz;
	//End AuxDebugTest

else

	//AuxDebugTest
		tmpx=bx[0][p][q]
		tmpz=bz[0][p][q]
		
		SetScale/P x DimOffset(bx,1),DimDelta(bx,1), tmpx,tmpz;
		SetScale/P y DimOffset(bx,2),DimDelta(bx,2), tmpx,tmpz;

		//srFFT2D(tmpx, -1); srFFT2D(tmpz, -1);
		//print "srFFT2D Inverse Used"

		j=0
		do
			i=0
			do
				bx[0][i][j]=tmpx[i][j]
				bz[0][i][j]=tmpz[i][j]
				i+=1
			while(i<nx)
			j+=1
		while (j<nz)
		
		SetScale/P y DimOffset(tmpx,0),DimDelta(tmpx,0),"m", bx,bz;
		SetScale/P z DimOffset(tmpx,1),DimDelta(tmpx,1),"m", bx,bz;
	//End AuxDebugTest

endif  //if (ne>1) 

	rad[2]=num2str(0)
	SetScale/P x str2num(rad[4]),str2num(rad[3]),"eV", bx,bz
	
	//SetScale/P y str2num(rad[6]),str2num(rad[5]),"m", bx,bz  //AuxDebugTest
	//SetScale/P z str2num(rad[8]),str2num(rad[7]),"m", bx,bz  //AuxDebugTest
	
endif  // if (r==0)
//print srround(stopmstimer(t0)/1e6,1)," seconds"

//if (ne>1)  //AuxDebugTest
killwaves/Z tmpx tmpz
//endif  //AuxDebugTest

endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Setting up Stokes Parameters structure from Electric Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfr2Sto(StoName, WfrName)
string WfrName=SrwRadName+SrwRadType
string StoName=WfrName[0,strlen(WfrName)-strlen(SrwRadType)-1]
prompt StoName, "Name Core for the Stokes Parameters structure"
prompt WfrName, SrwPViewRadName, popup Wavelist("*"+SrwRadType, ";", "")
Silent 1						|	Calculating Stokes Parameters ...
PauseUpdate

if(exists(WfrName) != 1)
	abort "Wavefront structure was not provided"
endif

SrwStoName = StoName
SrwRadName = WfrName[0,strlen(WfrName)-strlen(SrwRadType)-1]

string nmEX = $WfrName[0]
string nmEZ = $WfrName[1]

string nmEref = ""
variable EX_exists = 0, EZ_exists = 0
if(exists(nmEX) == 1)
	nmEref = nmEX
	EX_exists = 1
endif
if(exists(nmEZ) == 1)
	if(exists(nmEX) != 1)
		nmEref = nmEZ
	endif
	EZ_exists = 1
endif
if((strlen(nmEref) <= 0) %| (exists(nmEref) != 1))
	abort "Incorrect Wavefront structure"
endif

variable ne = dimsize($nmEref, 0)
variable eSt = dimoffset($nmEref, 0)
variable eStep = dimdelta($nmEref, 0)
variable eFi = eSt + eStep*(ne -1)
string eUnit = waveunits($nmEref, 0)
variable nx = dimsize($nmEref, 1)
variable xSt = dimoffset($nmEref, 1)
variable xStep = dimdelta($nmEref, 1)
variable xFi = xSt + xStep*(nx -1)
string xUnit = waveunits($nmEref, 1)
variable nz = dimsize($nmEref, 2)
variable zSt = dimoffset($nmEref, 2)
variable zStep = dimdelta($nmEref, 2)
variable zFi = zSt + zStep*(nz -1)
string zUnit = waveunits($nmEref, 2)

//SrwStoPrepSimpleStr(StoName, eSt, eFi, ne, xSt, xFi, nx, zSt, zFi, nz)
StoName += SrwStoType

Make/O/N=(4, ne, nx, nz) $StoName
SetScale/I y eSt, eFi,eUnit, $StoName
SetScale/I z xSt, xFi, xUnit, $StoName
SetScale/I t zSt, zFi, zUnit, $StoName
SetScale d 0, 0, SrwPUnitSpAngFluxPerUnSurf, $StoName


variable/C zeroC = cmplx(0, 0)
if(EX_exists)
	if(EZ_exists)
		$StoName = srwUtiE2Stokes($nmEX[q][r][s], $nmEZ[q][r][s], p)
	else
		$StoName = srwUtiE2Stokes($nmEX[q][r][s], zeroC, p)
	endif
else
	if(EZ_exists)
		$StoName = srwUtiE2Stokes(zeroC, $nmEZ[q][r][s], p)
	else
		abort "Incorrect Wavefront structure"
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
// Test of Propagation of the Radiation through a Beamline (uses C++ routine)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPropagTest(rad, RadNewName, blname, SliceType, eVal, xVal, zVal)
string rad=SrwRadName+SrwRadType;
string RadNewName=SrwRadName+"t";
string blname=SrwBliLast+SrwBeamlineType;
Variable SliceType=2;
Variable eVal=SrwViewE;
Variable xVal=SrwViewX;
Variable zVal=SrwViewZ;
prompt rad,SrwPRadName1,popup Wavelist("*"+SrwRadType ,";","");
prompt RadNewName,"Name of Test Wavefront struct.";
prompt blname,SrwPBliName,popup Wavelist("*"+SrwBeamlineType ,";","");
prompt SliceType,SrwPViewPlotType,popup "Energy;Hor.;Vert.;Hor. & Vert.";
prompt eVal, SrwPViewE;
prompt xVal, SrwPViewX;
prompt zVal, SrwPViewZ;
Silent 1						|	 ...
PauseUpdate;

SrwBliLast=blname[0,strlen(blname)-strlen(SrwBeamlineType)-1];
SrwRadName=rad[0,strlen(rad)-strlen(SrwRadType)-1];

string RadxNew=RadNewName+SrwSuffixX+SrwRadElType;
string RadzNew=RadNewName+SrwSuffixZ+SrwRadElType;
RadNewName += SrwRadType;
string ElecName = $rad[12];

string ObsNameNew = "Obs";
variable Dist = 0.;
SrwSmpCreate(ObsNameNew, Dist);

ObsNameNew = SrwSmpName+SrwSmpType;
variable Xmid, Xmagn, Xnpts, Zmid, Zmagn, Znpts, Edep, Efin, Enpts;
if(SliceType==1) // vs Energy
	Xmid = xVal; Xmagn = 0.; Xnpts = 1;
	Zmid = zVal; Zmagn = 0.; Znpts = 1;
	Edep = str2num($rad[4])*0.001; Enpts = DimSize($rad[0], 0); Efin = Edep + (Enpts-1)*str2num($rad[3])*0.001;
endif
if(SliceType==2) // vs Hor.
	variable Xdep = str2num($rad[6])*1000.;
	Xnpts = DimSize($($rad[0]), 1);
	if(Xnpts == 1)
		Xmid = Xdep; Xmagn = 0.;
	else
		variable Xfin = Xdep + (Xnpts-1)*str2num($rad[5])*1000.;
		Xmid = 0.5*(Xdep + Xfin); Xmagn = Xfin - Xdep;
	endif
	Zmid = zVal; Zmagn = 0.; Znpts = 1;
	Edep = eVal; Efin = eVal; Enpts = 1;
endif
if(SliceType==3) // vs Vert.
	Xmid = xVal; Xmagn = 0.; Xnpts = 1;
	variable Zdep = str2num($rad[8])*1000.;
	Znpts = DimSize($rad[0], 2);
	if(Znpts == 1)
		Zmid = Zdep; Zmagn = 0.;
	else
		variable Zfin = Zdep + (Znpts-1)*str2num($rad[7])*1000.;
		Zmid = 0.5*(Zdep + Zfin); Zmagn = Zfin - Zdep;
	endif
	Edep = eVal; Efin = eVal; Enpts = 1;
endif
if(SliceType==4) // vs Hor. & Vert.
	variable Xdep = str2num($rad[6])*1000.;
	Xnpts = DimSize($rad[0], 1);
	if(Xnpts == 1)
		Xmid = Xdep; Xmagn = 0.;
	else
		variable Xfin = Xdep + (Xnpts-1)*str2num($rad[5])*1000.;
		Xmid = 0.5*(Xdep + Xfin); Xmagn = Xfin - Xdep;
	endif
	variable Zdep = str2num($rad[8])*1000.;
	Znpts = DimSize($rad[0], 2);
	if(Znpts == 1)
		Zmid = Zdep; Zmagn = 0.;
	else
		variable Zfin = Zdep + (Znpts-1)*str2num($rad[7])*1000.;
		Zmid = 0.5*(Zdep + Zfin); Zmagn = Zfin - Zdep;
	endif
	Edep = eVal; Efin = eVal; Enpts = 1;
endif
SrwSmpScanXZE(ObsNameNew, Xmid, Xmagn, Xnpts, Zmid, Zmagn, Znpts, Edep, Efin, Enpts);

SrwRadPrep(ElecName, ObsNameNew, RadNewName, RadxNew, RadzNew, 0);
$RadNewName[2] = $rad[2];
$RadNewName[9] = $rad[9];
$($RadNewName[13]) = $($rad[13]);
$RadNewName[14] = $rad[14];
$($RadNewName[15]) = $($rad[15]);
$($RadNewName[16]) = $($rad[16]);
$RadNewName[17] = $rad[17];
$RadNewName[18] = $rad[18];
$RadNewName[19] = $rad[19];
$RadNewName[20] = $rad[20];
$RadNewName[21] = $rad[21];

string RadAuxName=SrwRadName+"Aux";
SrwWfrDupl(rad, RadAuxName);
RadAuxName += SrwRadType;

srRadPropagTest($RadAuxName, $blname, $RadNewName);

SrwRadKill(RadAuxName);
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwImageFormat(str)
string str

//ModifyImage $str ctab= {*,*,Rainbow,0}
// add other settings
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwPlotFormat(str)
string str

SrwUtiGraphAddFrameAndGrid()

end

//+++++++++++++++++++++++++++++++++++++++
function sRMatMult(a,b,c)
wave a,b,c
variable i,j,k
a=0

i=0
do
j=0
do
k=0
do
a[i][j]+=b[i][k]*c[k][j]
k=k+1
while(k<4)
j=j+1
while(j<4)
i=i+1
while(i<4)
end

//+++++++++++++++++++++++++++++++++++++++
function sRMatIni(a)
wave a
a=0
a[0][0]=1
a[1][1]=1
a[2][2]=1
a[3][3]=1
end

//+++++++++++++++++++++++++++++++++++++++
function sRMatDrift(a,length)
wave a;Variable length
sRMatIni(a)
a[0][1]=length
a[2][3]=length
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculates projected position or RMS size of electron beam for given wavefront
//if(m1or2 == 1) returns (x, y)
//if(m1or2 == 2) returns (SigmaX, SigmaY)
//
//+++++++++++++++++++++++++++++++++++++++
function/C srwWfrElecBeamProj(wfr, m1or2)
wave/T wfr
variable m1or2

wave wElec = $(wfr[12])
wave wMatr = $(wfr[13])

variable MX = wElec[2]
variable MXp = wElec[3]
variable MY = wElec[4]
variable MYp = wElec[5]

variable a00 = wMatr[0][0], a01 = wMatr[1][0], a02 = wMatr[2][0], a03 = wMatr[3][0]
variable a10 = wMatr[0][1], a11 = wMatr[1][1], a12 = wMatr[2][1], a13 = wMatr[3][1]
variable a20 = wMatr[0][2], a21 = wMatr[1][2], a22 = wMatr[2][2], a23 = wMatr[3][2]
variable a30 = wMatr[0][3], a31 = wMatr[1][3], a32 = wMatr[2][3], a33 = wMatr[3][3]

if(m1or2 == 1)
	return cmplx(a00*MX + a01*MXp + a02*MY + a03*MYp, a20*MX + a21*MXp + a22*MY + a23*MYp)
endif

variable MXX = wElec[20]
variable MXXp = wElec[21]
variable MXpXp = wElec[22]
variable MYY = wElec[23]
variable MYYp = wElec[24]
variable MYpYp = wElec[25]
variable MXY = wElec[26]
variable MXpY = wElec[27]
variable MXYp = wElec[28]
variable MXpYp = wElec[29]

variable b00=a00*a00, b01=2*a00*a01, b02=a01*a01, b03=a02*a02, b04=2*a02*a03, b05=a03*a03, b06=2*a00*a02, b07=2*a01*a02, b08=2*a00*a03, b09=2*a01*a03
variable b10=a00*a10, b11=a01*a10 + a00*a11, b12=a01*a11, b13=a02*a12, b14=a03*a12 + a02*a13, b15=a03*a13, b16=a02*a10 + a00*a12, b17=a02*a11 + a01*a12, b18=a03*a10 + a00*a13, b19=a03*a11 + a01*a13
variable b20=a10*a10, b21=2*a10*a11, b22=a11*a11, b23=a12*a12, b24=2*a12*a13, b25=a13*a13, b26=2*a10*a12, b27=2*a11*a12, b28=2*a10*a13, b29=2*a11*a13
variable b30=a20*a20, b31=2*a20*a21, b32=a21*a21, b33=a22*a22, b34=2*a22*a23, b35=a23*a23, b36=2*a20*a22, b37=2*a21*a22, b38=2*a20*a23, b39=2*a21*a23
variable b40=a20*a30, b41=a21*a30 + a20*a31, b42=21*a31, b43=a22*a32, b44=a23*a32 + a22*a33, b45=a23*a33, b46=a22*a30 + a20*a32, b47=a22*a31 + a21*a32, b48=a23*a30 + a20*a33, b49=a23*a31 + a21*a33
variable b50=a30*a30, b51=2*a30*a31, b52=a31*a31, b53=a32*a32, b54=2*a32*a33, b55=a33*a33, b56=2*a30*a32, b57=2*a31*a32, b58=2*a30*a33, b59=2*a31*a33
variable b60=a00*a20, b61=a01*a20 + a00*a21, b62=a01*a21, b63=a02*a22, b64=a03*a22 + a02*a23, b65=a03*a23, b66=a02*a20 + a00*a22, b67=a02*a21 + a01*a22, b68=a03*a20 + a00*a23, b69=a03*a21 + a01*a23
variable b70=a10*a20, b71=a11*a20 + a10*a21, b72=a11*a21, b73=a12*a22, b74=a13*a22 + a12*a23, b75=a13*a23, b76=a12*a20 + a10*a22, b77=a12*a21 + a11*a22, b78=a13*a20 + a10*a23, b79=a13*a21 + a11*a23
variable b80=a00*a30, b81=a01*a30 + a00*a31, b82=a01*a31, b83=a02*a32, b84=a03*a32 + a02*a33, b85=a03*a33, b86=a02*a30 + a00*a32, b87=a02*a31 + a01*a32, b88=a03*a30 + a00*a33, b89=a03*a31 + a01*a33
variable b90=a10*a30, b91=a11*a30 + a10*a31, b92=a11*a31, b93=a12*a32, b94=a13*a32 + a12*a33, b95=a13*a33, b96=a12*a30 + a10*a32, b97=a12*a31 + a11*a32, b98=a13*a30 + a10*a33, b99=a13*a31 + a11*a33

variable SigXe2 = 0, SigYe2 = 0
if(m1or2 == 2)
	SigXe2 = b00*MXX + b01*MXXp + b02*MXpXp + b03*MYY + b04*MYYp + b05*MYpYp + b06*MXY + b07*MXpY + b08*MXYp + b09*MXpYp
	SigYe2 = b30*MXX + b31*MXXp + b32*MXpXp + b33*MYY + b34*MYYp + b35*MYpYp + b36*MXY + b37*MXpY + b38*MXYp + b39*MXpYp
endif
return cmplx(sqrt(abs(SigXe2)), sqrt(abs(SigYe2)))
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculates initial wavefronts at different photon energies 
//and propagates them in automatic mode
//SpecSensName is Relative Spectral Sensitivity vs Phot. Energy in [eV]
//Creates Intesity wave named IntName + "_xz", in [W/mm^2]
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwWfrCreatePropagChrom(IntName, ElecName, MagName, ObsName, BLName, ObsNxNzSamplFact, Prec, PolCmpn, MultiElecCmpn, SpecSensName)
proc SrwWfrCreatePropagChrom(IntName, ElecName, MagName, ObsName, BLName, Prec, PolCmpn, MultiElecCmpn, SpecSensName, Disp)
string IntName=srwUtiTruncString(SrwElecName+SrwMagName+SrwSmpName+SrwBliLast, 25)
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagName+SrwFieldType
string ObsName=SrwSmpName+SrwSmpType
string BLName=SrwBliLast+SrwBeamlineType
//variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
variable Prec=srwUtiGetValN("Prec", 1, "SrwWfrCreatePropagChrom")
variable PolCmpn=SrwViewRadCmpnType
variable MultiElecCmpn=srwUtiGetValN("MultiElecCmpn", 1, "SrwWfrCreatePropagChrom")
string SpecSensName=srwUtiGetValS("SpecSensName", "", "SrwWfrCreatePropagChrom")
variable Disp=srwUtiGetValN("Disp", 1, "SrwWfrCreatePropagChrom")
prompt IntName, "Name of the Final Intensity structure"
prompt ElecName,"Electron Beam structure", popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2, popup Wavelist("*"+SrwFieldType,";","")
prompt ObsName,SrwPSmpName2, popup Wavelist("*"+SrwSmpType,";","")
prompt BLName,SrwPBliName, popup Wavelist("*"+SrwBeamlineType,";","")
//prompt ObsNxNzSamplFact, SrwPSmpNxNzSamplFact
prompt Prec, "Overal Precision Factor"
prompt PolCmpn, "Polarization Component to extract", popup SrwPOPUPPolar+";Total"
prompt MultiElecCmpn, "Multi-e or Single-e Intensity to extract", popup "Single-e Intensity;Multi-e Intensity;Single-e Flux;Multi-e Flux"
prompt SpecSensName, "Rel. Spectral Sensitivity vs ph. en. in [eV]", popup Wavelist("*",";","")
prompt Disp, "Display Intermed. and Final Results?", popup "No;Yes"
Silent 1						|	  ...
//PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwBliLast=BLName[0,strlen(BLName)-strlen(SrwBeamlineType)-1]
//SrwSmpNxNzSamplFact=ObsNxNzSamplFact
srwUtiSetValN("Prec", Prec, "SrwWfrCreatePropagChrom")
SrwViewRadCmpnType=PolCmpn
srwUtiSetValN("MultiElecCmpn", MultiElecCmpn, "SrwWfrCreatePropagChrom")
srwUtiSetValS("SpecSensName", SpecSensName, "SrwWfrCreatePropagChrom")
srwUtiSetValN("Disp", Disp, "SrwWfrCreatePropagChrom")

variable minInitNpX = 60 //to steer
variable minInitNpZ = 60
variable maxPropNpX = 2000 //to steer
variable maxPropNpZ = 2000

variable InitWfrCompPrec = 0.01/Prec
variable ObsNxNzSamplFact = 1.2*Prec
SrwMagPrec(MagName, 3, 0.01, InitWfrCompPrec, 10000, 1, 0, 0)

variable IntensCmpnNo = MultiElecCmpn
if(MultiElecCmpn == 3)
	IntensCmpnNo = 5
endif
if(MultiElecCmpn == 4)
	IntensCmpnNo = 6
endif

string AuxObsName = "AuxObs_WfrCreatePropagChrom", AuxWfrName = "Aux_WfrCrPrCh"
SrwSmpCreate(AuxObsName, srwGetSmpLongPos(ObsName))

string FinIntName = IntName + "_xz"
string IntermedIntName = AuxWfrName + "_xz"
string IntermedIntName1 = AuxWfrName + "1_xz"
string IntermedIntNameX = AuxWfrName + "_x"
string IntermedIntNameZ = AuxWfrName + "_z"
string strLabel = "", strWinNameCutX, strWinNameCutZ
string strExe
string nmRadCmpnX = AuxWfrName + "X_rae"
string nmRadCmpnZ = AuxWfrName + "Z_rae"

variable xStart_mm = 1000.*srwGetSmpHorPosStart(ObsName), xEnd_mm = 1000.*srwGetSmpHorPosEnd(ObsName), xNp = srwGetSmpHorPosNp(ObsName)
variable xc_mm = 0.5*(xStart_mm + xEnd_mm), xRange_mm = xEnd_mm - xStart_mm
variable zStart_mm = 1000.*srwGetSmpVertPosStart(ObsName), zEnd_mm = 1000.*srwGetSmpVertPosEnd(ObsName), zNp = srwGetSmpVertPosNp(ObsName)
variable zc_mm = 0.5*(zStart_mm + zEnd_mm), zRange_mm = zEnd_mm - zStart_mm
variable eStart_keV = 0.001*srwGetSmpPhotEnStart(ObsName), eEnd_keV = 0.001*srwGetSmpPhotEnEnd(ObsName), eNp = srwGetSmpPhotEnNp(ObsName)
variable eStep_keV = 0
if(eNp > 0)
	eStep_keV = (eEnd_keV - eStart_keV)/(eNp - 1)
endif

variable resizeInitXD, resizeInitZD
variable e_keV = eStart_keV, RelSpecSensitivity = 1

SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, ObsNxNzSamplFact)
variable curNpX = dimsize($nmRadCmpnX, 1), curNpZ = dimsize($nmRadCmpnX, 2)
resizeInitXD = minInitNpX/curNpX
if(resizeInitXD < 1.2)
	resizeInitXD = 1
endif
resizeInitZD = minInitNpZ/curNpZ
if(resizeInitZD < 1.2)
	resizeInitZD = 1
endif
if((resizeInitXD != 1) %| (resizeInitZD != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
endif
//SrwWfrPropag(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 1, "")
SrwWfrProp(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 2, 1, "")
curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
resizeInitXD = maxPropNpX/curNpX
if(resizeInitXD > 0.8)
	resizeInitXD = 1
endif
resizeInitZD = maxPropNpZ/curNpZ
if(resizeInitZD > 0.8)
	resizeInitZD = 1
endif
if((resizeInitXD != 1) %| (resizeInitZD != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
endif
SrwWfrDupl(AuxWfrName+SrwRadType, AuxWfrName + "D1")

SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
RelSpecSensitivity = $SpecSensName(e_keV*1000)

if(Disp == 2)
	strLabel =  num2str(e_keV) + " keV"
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 2)
	SrwUtiGraphAddFrameAndGrid()
	SrwUtiGraphWindResize(2,242+170,190,150,-1,0)
	$IntermedIntNameX *= RelSpecSensitivity
	strWinNameCutX = WinName(0, 1)
	strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 2)
	SrwUtiGraphAddFrameAndGrid()
	SrwUtiGraphWindResize(200,242+170,190,150,-1,0)
	$IntermedIntNameZ *= RelSpecSensitivity
	strWinNameCutZ = WinName(0, 1)
	strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
endif

duplicate/O $IntermedIntName $IntermedIntName1
$IntermedIntName1 *= (0.5*RelSpecSensitivity)

variable xStartEmin = dimoffset($IntermedIntName1, 0), xStepEmin = dimdelta($IntermedIntName1, 0), xNpEmin = dimsize($IntermedIntName1, 0)
variable xEndEmin = xStartEmin + xStepEmin*(xNpEmin - 1)
variable zStartEmin = dimoffset($IntermedIntName1, 1), zStepEmin = dimdelta($IntermedIntName1, 1), zNpEmin = dimsize($IntermedIntName1, 1)
variable zEndEmin = zStartEmin + zStepEmin*(zNpEmin - 1)

srUtiSpotInfo($IntermedIntName1)
string nmIntermedIntName1inf = IntermedIntName1 + "_inf"
variable xcAuxInt = $nmIntermedIntName1inf[1], zcAuxInt = $nmIntermedIntName1inf[2]
variable xFWHM = $nmIntermedIntName1inf[3], zFWHM = $nmIntermedIntName1inf[3]
variable xStartEminCor = xcAuxInt - 8*xFWHM, zStartEminCor = zcAuxInt - 8*zFWHM
variable xEndEminCor = xcAuxInt + 8*xFWHM, zEndEminCor = zcAuxInt + 8*zFWHM
if(xStartEmin < xStartEminCor)
	xStartEmin = xStartEminCor
endif
if(xEndEmin > xEndEminCor)
	xEndEmin = xEndEminCor
endif
if(xStartEmin > -xEndEmin)
	xStartEmin = -xEndEmin
endif
if(xEndEmin < -xStartEmin)
	xEndEmin = -xStartEmin
endif
xNpEmin = trunc((xEndEmin - xStartEmin)/xStepEmin) + 1
xStepEmin = (xEndEmin - xStartEmin)/(xNpEmin - 1)

if(zStartEmin < zStartEminCor)
	zStartEmin = zStartEminCor
endif
if(zEndEmin > zEndEminCor)
	zEndEmin = zEndEminCor
endif
if(zStartEmin > -zEndEmin)
	zStartEmin = -zEndEmin
endif
if(zEndEmin < -zStartEmin)
	zEndEmin = -zStartEmin
endif
zNpEmin = trunc((zEndEmin - zStartEmin)/zStepEmin) + 1
zStepEmin = (zEndEmin - zStartEmin)/(zNpEmin - 1)

e_keV = eEnd_keV
SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, ObsNxNzSamplFact)
curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
resizeInitXD = minInitNpX/curNpX
if(resizeInitXD < 1.2)
	resizeInitXD = 1
endif
resizeInitZD = minInitNpZ/curNpZ
if(resizeInitZD < 1.2)
	resizeInitZD = 1
endif
if((resizeInitXD != 1) %| (resizeInitZD != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
endif
//SrwWfrPropag(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 1, "")
SrwWfrProp(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 2, 1, "")
curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
resizeInitXD = maxPropNpX/curNpX
if(resizeInitXD > 0.8)
	resizeInitXD = 1
endif
resizeInitZD = maxPropNpZ/curNpZ
if(resizeInitZD > 0.8)
	resizeInitZD = 1
endif
if((resizeInitXD != 1) %| (resizeInitZD != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
endif
SrwWfrDupl(AuxWfrName+SrwRadType, AuxWfrName + "D")

SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
RelSpecSensitivity = $SpecSensName(e_keV*1000)

if(Disp == 2)
	strLabel =  num2str(e_keV) + " keV"
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 1)
	strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 1)
	strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	$IntermedIntNameX *= RelSpecSensitivity
	$IntermedIntNameZ *= RelSpecSensitivity
endif

variable xStartEmax = dimoffset($IntermedIntName, 0), xStepEmax = dimdelta($IntermedIntName, 0), xNpEmax = dimsize($IntermedIntName, 0)
variable xEndEmax = xStartEmax + xStepEmax*(xNpEmax - 1)
variable zStartEmax = dimoffset($IntermedIntName, 1), zStepEmax = dimdelta($IntermedIntName, 1), zNpEmax = dimsize($IntermedIntName, 1)
variable zEndEmax = zStartEmax + zStepEmax*(zNpEmax - 1)

srUtiSpotInfo($IntermedIntName)
string nmIntermedIntNameInf = IntermedIntName + "_inf"
xcAuxInt = $nmIntermedIntNameInf[1]; zcAuxInt = $nmIntermedIntNameInf[2]
xFWHM = $nmIntermedIntNameInf[3]; zFWHM = $nmIntermedIntNameInf[3]
variable xStartEmaxCor = xcAuxInt - 8*xFWHM, zStartEmaxCor = zcAuxInt - 8*zFWHM
variable xEndEmaxCor = xcAuxInt + 8*xFWHM, zEndEmaxCor = zcAuxInt + 8*zFWHM
if(xStartEmax < xStartEmaxCor)
	xStartEmax = xStartEmaxCor
endif
if(xEndEmax > xEndEmaxCor)
	xEndEmax = xEndEmaxCor
endif
if(xStartEmax > -xEndEmax)
	xStartEmax = -xEndEmax
endif
if(xEndEmax < -xStartEmax)
	xEndEmax = -xStartEmax
endif
xNpEmax = trunc((xEndEmax - xStartEmax)/xStepEmax) + 1
xStepEmax = (xEndEmax - xStartEmax)/(xNpEmax - 1)

if(zStartEmax < zStartEmaxCor)
	zStartEmax = zStartEmaxCor
endif
if(zEndEmax > zEndEmaxCor)
	zEndEmax = zEndEmaxCor
endif

if(zStartEmax > -zEndEmax)
	zStartEmax = -zEndEmax
endif
if(zEndEmax < -zStartEmax)
	zEndEmax = -zStartEmax
endif
zNpEmax = trunc((zEndEmax - zStartEmax)/zStepEmax) + 1
zStepEmax = (zEndEmax - zStartEmax)/(zNpEmax - 1)

variable xStepFin = xStepEmax, xStartFin = xStartEmin,  xEndFin = xEndEmin
if(xStepFin > xStepEmin)
	xStepFin = xStepEmin
endif
//xStepFin = 0.2*xStepEmin + 0.8*xStepEmax
if(xStartFin > xStartEmax)
	xStartFin = xStartEmax
endif
if(xEndFin < xEndEmin)
	xEndFin = xEndEmin
endif
//variable xNpFin = round((xEndFin - xStartFin)/(xStepFin))
//xNpFin = 2*trunc(0.5*xNpFin + 0.00001)
variable xNpFin = trunc((xEndFin - xStartFin)/(xStepFin)) + 1
xStepFin = (xEndFin - xStartFin)/(xNpFin - 1)
if(xNpFin > maxPropNpX)
	xNpFin = maxPropNpX
	//xNpFin = 2*trunc(0.5*xNpFin + 0.00001)
	xStepFin = (xEndFin - xStartFin)/(xNpFin - 1)
endif
//if(xNpFin < xNpEmin)
//	xNpFin = xNpEmin
//endif
//if(xNpFin < xNpEmax)
//	xNpFin = xNpEmax
//endif

variable zStepFin = zStepEmax, zStartFin = zStartEmin,  zEndFin = zEndEmin
if(zStepFin > zStepEmin)
	zStepFin = zStepEmin
endif
//zStepFin = 0.2*zStepEmin + 0.8*zStepEmax
if(zStartFin > zStartEmax)
	zStartFin = zStartEmax
endif
if(zEndFin < zEndEmin)
	zEndFin = zEndEmin
endif
//variable zNpFin = round((zEndFin - zStartFin)/(zStepFin))
//zNpFin = 2*trunc(0.5*zNpFin + 0.00001)
variable zNpFin = trunc((zEndFin - zStartFin)/(zStepFin)) + 1
zStepFin = (zEndFin - zStartFin)/(zNpFin - 1)
if(zNpFin > maxPropNpZ)
	zNpFin = maxPropNpZ
	//zNpFin = 2*trunc(0.5*zNpFin + 0.00001)
	zStepFin = (zEndFin - zStartFin)/(zNpFin - 1)
endif
//if(zNpFin < zNpEmin)
//	zNpFin = zNpEmin
//endif
//if(zNpFin < zNpEmax)
//	zNpFin = zNpEmax
//endif

variable resizeXD = abs(xStepEmin/xStepFin)
variable resizeZD = abs(zStepEmin/zStepFin)
if(resizeXD < 1.2)
	resizeXD = 1
endif
if(resizeZD < 1.2)
	resizeZD = 1
endif
if((resizeXD > 1) %| (resizeXD > 1))
	SrwWfrResize(AuxWfrName+"D1"+SrwRadType,1,1,resizeXD,1,resizeZD,2,AuxWfrName)
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, eStart_keV, 0, 0, 1)
	RelSpecSensitivity = $SpecSensName(eStart_keV*1000)
	duplicate/O $IntermedIntName $IntermedIntName1
	$IntermedIntName1 *= (0.5*RelSpecSensitivity)
endif
SrwWfrDel(AuxWfrName+"D1"+SrwRadType)

resizeXD = abs(xStepEmax/xStepFin)
resizeZD = abs(zStepEmax/zStepFin)
if(resizeXD < 1.2)
	resizeXD = 1
endif
if(resizeZD < 1.2)
	resizeZD = 1
endif
if((resizeXD > 1) %| (resizeXD > 1))
	SrwWfrResize(AuxWfrName+"D"+SrwRadType,1,1,resizeXD,1,resizeZD,2,AuxWfrName)
endif
SrwWfrDupl(AuxWfrName+"D"+SrwRadType,AuxWfrName)
SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, eEnd_keV, 0, 0, 1)
RelSpecSensitivity = $SpecSensName(eEnd_keV*1000)
SrwWfrDel(AuxWfrName+"D"+SrwRadType)

make/O/N=(xNpFin, zNpFin) $FinIntName
setscale/P x xStartFin, xStepFin, "m", $FinIntName
setscale/P y zStartFin, zStepFin, "m", $FinIntName

string FinIntNameX = FinIntName[0,strlen(FinIntName)-strlen("_xz")-1] + "_x"
string FinIntNameZ = FinIntName[0,strlen(FinIntName)-strlen("_xz")-1] + "_z"

variable xStartCur = dimoffset($IntermedIntName1, 0), xStepCur = dimdelta($IntermedIntName1, 0), xNpCur = dimsize($IntermedIntName1, 0)
variable xEndCur = xStartCur + (xNpCur - 1)*xStepCur
variable zStartCur = dimoffset($IntermedIntName1, 1), zStepCur = dimdelta($IntermedIntName1, 1), zNpCur = dimsize($IntermedIntName1, 1)
variable zEndCur = zStartCur + (zNpCur - 1)*zStepCur
$FinIntName = $IntermedIntName1(x)(y)*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)

xStartCur = dimoffset($IntermedIntName, 0); xStepCur = dimdelta($IntermedIntName, 0); xNpCur = dimsize($IntermedIntName, 0)
xEndCur = xStartCur + (xNpCur - 1)*xStepCur
zStartCur = dimoffset($IntermedIntName, 1); zStepCur = dimdelta($IntermedIntName, 1); zNpCur = dimsize($IntermedIntName, 1)
zEndCur = zStartCur + (zNpCur - 1)*zStepCur
$FinIntName += (0.5*RelSpecSensitivity)*($IntermedIntName(x)(y))*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)
variable xc = 0.5*(xStartFin + xEndFin), zc = 0.5*(zStartFin + zEndFin)

variable ConstPhotToW = 1.60219e-16 // 1 Phot/s/.1%bw correspond(s) to: 1.60219e-16  W/eV
variable IntMultConst = ConstPhotToW*eStep_keV*1000
$FinIntName *= IntMultConst

if(Disp == 2)
	display; appendimage $FinIntName
	label bottom "Horizontal Position"
	label left "Vertical Position"
	SrwUtiGraphWindResize(2,10,190,170,-2,0)
	display; appendimage $IntermedIntName
	label bottom "Horizontal Position"
	label left "Vertical Position"
	SrwUtiGraphWindResize(200,10,190,170,-2,0)
	
	make/O/N=(xNpFin) $FinIntNameX
	setscale/P x xStartFin, xStepFin, "m", $FinIntNameX
	$FinIntNameX = $FinIntName(x)(zc)
	display $FinIntNameX
	label bottom "Horizontal Position"
	label left "W/mm\\S2\\M"
	SrwUtiGraphAddFrameAndGrid()
	SrwUtiGraphWindResize(2,242,190,150,-1,0)

	make/O/N=(zNpFin) $FinIntNameZ
	setscale/P x zStartFin, zStepFin, "m", $FinIntNameZ
	$FinIntNameZ = $FinIntName(xc)(x)
	display $FinIntNameZ
	label bottom "Vertical Position"
	label left "W/mm\\S2\\M"
	SrwUtiGraphAddFrameAndGrid()
	SrwUtiGraphWindResize(200,242,190,150,-1,0)
endif

variable ie = 1, AuxConst
e_keV = eStart_keV + eStep_keV

do 
	SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
	SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, ObsNxNzSamplFact)
	curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
	resizeInitXD = minInitNpX/curNpX
	if(resizeInitXD < 1.2)
		resizeInitXD = 1
	endif
	resizeInitZD = minInitNpZ/curNpZ
	if(resizeInitZD < 1.2)
		resizeInitZD = 1
	endif
	if((resizeInitXD != 1) %| (resizeInitZD != 1))
		SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
	endif
	//SrwWfrPropag(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 1, "")
	SrwWfrProp(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 2, 1, "")
	curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
	resizeInitXD = maxPropNpX/curNpX
	if(resizeInitXD > 0.8)
		resizeInitXD = 1
	endif
	resizeInitZD = maxPropNpZ/curNpZ
	if(resizeInitZD > 0.8)
		resizeInitZD = 1
	endif
	if((resizeInitXD != 1) %| (resizeInitZD != 1))
		SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
	endif
	
	xStepCur = dimdelta($nmRadCmpnX, 1)
	zStepCur = dimdelta($nmRadCmpnX, 2)
	
	resizeXD = abs(xStepCur/xStepFin)
	resizeZD = abs(zStepCur/zStepFin)
	if(resizeXD < 1.2)
		resizeXD = 1
	endif
	if(resizeZD < 1.2)
		resizeZD = 1
	endif
	if((resizeXD > 1) %| (resizeXD > 1))
		SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeXD,1,resizeZD,1,"")
	endif
	
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
	RelSpecSensitivity = $SpecSensName(e_keV*1000)
	
	if(Disp == 2)
		strLabel =  num2str(e_keV) + " keV"
		SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 1)
		strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
		SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 1)
		strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
		$IntermedIntNameX *= RelSpecSensitivity
		$IntermedIntNameZ *= RelSpecSensitivity
	endif
	
	AuxConst = IntMultConst*RelSpecSensitivity
	
	xStartCur = dimoffset($nmRadCmpnX, 1); xStepCur = dimdelta($nmRadCmpnX, 1); xNpCur = dimsize($nmRadCmpnX, 1)
	xEndCur = xStartCur + (xNpCur - 1)*xStepCur
	zStartCur = dimoffset($nmRadCmpnX, 2); zStepCur = dimdelta($nmRadCmpnX, 2); zNpCur = dimsize($nmRadCmpnX, 2)
	zEndCur = zStartCur + (zNpCur - 1)*zStepCur
	$FinIntName += AuxConst*($IntermedIntName(x)(y))*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)
	
	if(Disp == 2)
		$FinIntNameX = $FinIntName(x)(zc)
		$FinIntNameZ = $FinIntName(xc)(x)
	endif
	
	e_keV += eStep_keV
	ie += 1
while(ie < (eNp - 1))

SrwWfrDel(AuxWfrName+SrwRadType)
killwaves/Z IntermedIntName, IntermedIntName1
end

//+++++++++++++++++++++++++++++++++++++++
//
//Defines precision parameters for SrwWfrCreatePropagChromS
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrPrecCreatePropagChromS(nmPrec, precInit, overSampInit, resizeRangeX, resizeResolX, resizeRangeZ, resizeResolZ)
string nmPrec = srwUtiGetValS("nmPrec", "", "SrwWfrCreatePropagChromS")
variable precInit = srwUtiGetValN("precInit", 0.002, "SrwWfrCreatePropagChromS")
variable overSampInit = srwUtiGetValN("overSampInit", 1, "SrwWfrCreatePropagChromS")
variable resizeRangeX = srwUtiGetValN("resizeRangelX", 1, "SrwWfrCreatePropagChromS")
variable resizeResolX = srwUtiGetValN("resizeResolX", 1, "SrwWfrCreatePropagChromS")
variable resizeRangeZ = srwUtiGetValN("resizeRangelZ", 1, "SrwWfrCreatePropagChromS")
variable resizeResolZ = srwUtiGetValN("resizeResolZ", 1, "SrwWfrCreatePropagChromS")
prompt nmPrec, "Name for Precision Param. structure"
prompt precInit, "Precision Factor for Retarded Potential Integ."
prompt overSampInit, "Oversampling Factor for Initial Wavefront Calc."
prompt resizeRangeX, "Horizontal Range Resizing Factor"
prompt resizeResolX, "Horizontal Resolution Resizing Factor"
prompt resizeRangeZ, "Vertical Range Resizing Factor"
prompt resizeResolZ, "Vertical Resolution Resizing Factor"
Silent 1						|	  ...
PauseUpdate

srwUtiSetValS("nmPrec", nmPrec, "SrwWfrCreatePropagChromS")
srwUtiSetValN("precInit", precInit, "SrwWfrCreatePropagChromS")
srwUtiSetValN("overSampInit", overSampInit, "SrwWfrCreatePropagChromS")
srwUtiSetValN("resizeRangeX", resizeRangeX, "SrwWfrCreatePropagChromS")
srwUtiSetValN("resizeResolX", resizeResolX, "SrwWfrCreatePropagChromS")
srwUtiSetValN("resizeRangeZ", resizeRangeZ, "SrwWfrCreatePropagChromS")
srwUtiSetValN("resizeResolZ", resizeResolZ, "SrwWfrCreatePropagChromS")

make/O/N=6 $nmPrec
$nmPrec[0] = precInit
$nmPrec[1] = overSampInit
$nmPrec[2] = resizeRangeX
$nmPrec[3] = resizeResolX
$nmPrec[4] = resizeRangeZ
$nmPrec[5] = resizeResolZ
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculates initial wavefronts at different photon energies 
//and propagates them in MANUAL mode
//SpecSensName is Relative Spectral Sensitivity vs Phot. Energy in [eV]
//Creates Intesity wave named IntName + "_xz", in [W/mm^2]
//Written to simulate propagation of partially-coherent wavefront through Zone Plate setup
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrCreatePropagChromS(IntName, ElecName, MagName, ObsName, BLName, PolCmpn, SpecSensName, nmPrecParam, ebmPropSigX, ebmPropSigZ)
string IntName=srwUtiGetValS("IntName", "Int", "SrwWfrCreatePropagChromS") //srwUtiTruncString(SrwElecName+SrwMagName+SrwSmpName+SrwBliLast, 25)
string ElecName=SrwElecName+SrwElecType
string MagName=SrwMagName+SrwFieldType
string ObsName=srwUtiGetValS("ObsName", "", "SrwWfrCreatePropagChromS")
string BLName=SrwBliLast+SrwBeamlineType
//variable Prec=srwUtiGetValN("Prec", 1, "SrwWfrCreatePropagChromS")
variable PolCmpn=SrwViewRadCmpnType
//variable MultiElecCmpn=srwUtiGetValN("MultiElecCmpn", 1, "SrwWfrCreatePropagChromS")
//variable Disp=srwUtiGetValN("Disp", 1, "SrwWfrCreatePropagChromS")
string SpecSensName=srwUtiGetValS("SpecSensName", "", "SrwWfrCreatePropagChromS")
string nmPrecParam=srwUtiGetValS("nmPrec", "", "SrwWfrCreatePropagChromS")
variable ebmPropSigX = srwUtiGetValN("ebmPropSigX", 0, "SrwWfrCreatePropagChromS")
variable ebmPropSigZ = srwUtiGetValN("ebmPropSigZ", 0, "SrwWfrCreatePropagChromS")
prompt IntName, "Name of the Final Intensity structure"
prompt ElecName,"Filament Electron Beam structure", popup Wavelist("*"+SrwElecType,";","")
prompt MagName,SrwPMagName2, popup Wavelist("*"+SrwFieldType,";","")
prompt ObsName,SrwPSmpName2, popup Wavelist("*"+SrwSmpType,";","")
prompt BLName,SrwPBliName, popup Wavelist("*"+SrwBeamlineType,";","")
//prompt Prec, "Overal Precision Factor"
prompt PolCmpn, "Polarization Component to extract", popup SrwPOPUPPolar+";Total"
//prompt MultiElecCmpn, "Multi-e or Single-e Intensity to extract", popup "Single-e Intensity;Multi-e Intensity;Single-e Flux;Multi-e Flux"
prompt SpecSensName, "Rel. Spectral Sensitivity vs ph. en. in [eV]", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt nmPrecParam, "Precision Parameters structure", popup Wavelist("*",";","TEXT:0,DIMS:1")
prompt ebmPropSigX, "Projected Horizontal E-Beam Size [µm]"
prompt ebmPropSigZ, "Projected Vertical E-Beam Size [µm]"
//prompt Disp, "Display Intermed. and Final Results?", popup "No;Yes"
Silent 1						|	  ...
//PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagName=MagName[0,strlen(MagName)-strlen(SrwFieldType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwBliLast=BLName[0,strlen(BLName)-strlen(SrwBeamlineType)-1]
//srwUtiSetValN("Prec", Prec, "SrwWfrCreatePropagChromS")
SrwViewRadCmpnType=PolCmpn
//srwUtiSetValN("MultiElecCmpn", MultiElecCmpn, "SrwWfrCreatePropagChromS")
srwUtiSetValS("ObsName", ObsName, "SrwWfrCreatePropagChromS")
srwUtiSetValS("SpecSensName", SpecSensName, "SrwWfrCreatePropagChromS")
srwUtiSetValS("nmPrec", nmPrecParam, "SrwWfrCreatePropagChromS")
srwUtiSetValN("ebmPropSigX", ebmPropSigX, "SrwWfrCreatePropagChromS")
srwUtiSetValN("ebmPropSigZ", ebmPropSigZ, "SrwWfrCreatePropagChromS")
//srwUtiSetValN("Disp", Disp, "SrwWfrCreatePropagChromS")
srwUtiSetValS("IntName", IntName, "SrwWfrCreatePropagChromS")

ebmPropSigX *= 1.e-06
ebmPropSigZ *= 1.e-06

SrwMagPrec(MagName, 3, 0.01, $nmPrecParam[0], 10000, 1, 0, 0)

//variable IntensCmpnNo = MultiElecCmpn
//if(MultiElecCmpn == 3)
//	IntensCmpnNo = 5
//endif
//if(MultiElecCmpn == 4)
//	IntensCmpnNo = 6
//endif
variable IntensCmpnNo = 1

string AuxObsName = "AuxObs_WfrCreatePropagChrom", AuxWfrName = "Aux_WfrCrPrCh"
SrwSmpCreate(AuxObsName, srwGetSmpLongPos(ObsName))

string FinIntName = IntName + "_xz"
string FinIntNameX = FinIntName[0,strlen(FinIntName)-strlen("_xz")-1] + "_x"
string FinIntNameZ = FinIntName[0,strlen(FinIntName)-strlen("_xz")-1] + "_z"
string IntermedIntName = AuxWfrName + "_xz"
//string IntermedIntName1 = AuxWfrName + "1_xz"
string IntermedIntNameX = AuxWfrName + "_x"
string IntermedIntNameZ = AuxWfrName + "_z"
string strLabel = "", strWinNameCutX, strWinNameCutZ
string strExe
string nmRadCmpnX = AuxWfrName + "X_rae"
string nmRadCmpnZ = AuxWfrName + "Z_rae"

variable xStart_mm = 1000.*srwGetSmpHorPosStart(ObsName), xEnd_mm = 1000.*srwGetSmpHorPosEnd(ObsName), xNp = srwGetSmpHorPosNp(ObsName)
variable xc_mm = 0.5*(xStart_mm + xEnd_mm), xRange_mm = xEnd_mm - xStart_mm
variable zStart_mm = 1000.*srwGetSmpVertPosStart(ObsName), zEnd_mm = 1000.*srwGetSmpVertPosEnd(ObsName), zNp = srwGetSmpVertPosNp(ObsName)
variable zc_mm = 0.5*(zStart_mm + zEnd_mm), zRange_mm = zEnd_mm - zStart_mm
variable eStart_keV = 0.001*srwGetSmpPhotEnStart(ObsName), eEnd_keV = 0.001*srwGetSmpPhotEnEnd(ObsName), eNp = srwGetSmpPhotEnNp(ObsName)
variable eStep_keV = 0
if(eNp > 0)
	eStep_keV = (eEnd_keV - eStart_keV)/(eNp - 1)
endif

variable resizeInitRangeX = $nmPrecParam[2]
variable resizeInitResolX = $nmPrecParam[3]
variable resizeInitRangeZ = $nmPrecParam[4]
variable resizeInitResolZ = $nmPrecParam[5]
variable e_keV = eStart_keV, RelSpecSensitivity = 1

variable ConstPhotToW = 1.60219e-16 // 1 Phot/s/.1%bw correspond(s) to: 1.60219e-16  W/eV
variable IntMultConst = ConstPhotToW*eStep_keV*1000

SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, $nmPrecParam[1])
if((resizeInitRangeX != 1) %| (resizeInitResolX != 1) %| (resizeInitRangeZ != 1) %| (resizeInitResolZ != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,resizeInitRangeX,resizeInitResolX,resizeInitRangeZ,resizeInitResolZ,1,"")
endif

SrwWfrProp(AuxWfrName+SrwRadType, BLName, 2, 2, 1, 2, 1, "")

SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
RelSpecSensitivity = $SpecSensName(e_keV*1000)
variable AuxConst = 0.5*IntMultConst*RelSpecSensitivity

//if(Disp == 2)
strLabel =  num2str(e_keV) + " keV"
SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 2)
SrwUtiGraphAddFrameAndGrid()
SrwUtiGraphWindResize(2,242+170,190,150,-1,0)
$IntermedIntNameX *= RelSpecSensitivity
strWinNameCutX = WinName(0, 1)
strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	
SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 2)
SrwUtiGraphAddFrameAndGrid()
SrwUtiGraphWindResize(200,242+170,190,150,-1,0)
$IntermedIntNameZ *= RelSpecSensitivity
strWinNameCutZ = WinName(0, 1)
strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
//endif

//duplicate/O $IntermedIntName $IntermedIntName1
//$IntermedIntName1 *= (0.5*RelSpecSensitivity)

if((ebmPropSigX != 0) %| (ebmPropSigZ != 0))
	SrwUtiConvWaveWithGaus2D(IntermedIntName, ebmPropSigX, ebmPropSigZ)
endif

duplicate/O $IntermedIntName $FinIntName
//$FinIntName *= (0.5*RelSpecSensitivity)
$FinIntName *= AuxConst
killwaves/Z $IntermedIntName

e_keV = eEnd_keV
SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, $nmPrecParam[1])
if((resizeInitRangeX != 1) %| (resizeInitResolX != 1) %| (resizeInitRangeZ != 1) %| (resizeInitResolZ != 1))
	SrwWfrResize(AuxWfrName+SrwRadType,1,resizeInitRangeX,resizeInitResolX,resizeInitRangeZ,resizeInitResolZ,1,"")
endif

SrwWfrProp(AuxWfrName+SrwRadType, BLName, 2, 2, 1, 2, 1, "")

SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
RelSpecSensitivity = $SpecSensName(e_keV*1000)
AuxConst = 0.5*IntMultConst*RelSpecSensitivity

//if(Disp == 2)
strLabel =  num2str(e_keV) + " keV"
SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 1)
strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 1)
strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
$IntermedIntNameX *= RelSpecSensitivity
$IntermedIntNameZ *= RelSpecSensitivity
//endif

if((ebmPropSigX != 0) %| (ebmPropSigZ != 0))
	SrwUtiConvWaveWithGaus2D(IntermedIntName, ebmPropSigX, ebmPropSigZ)
endif

$FinIntName += AuxConst*$IntermedIntName(x)(y)

variable xStartFin = dimoffset($IntermedIntName, 0), xStepFin = dimdelta($IntermedIntName, 0), xNpFin = dimsize($IntermedIntName, 0)
variable xEndFin = xStartFin + xStepFin*(xNpFin - 1)
variable zStartFin = dimoffset($IntermedIntName, 1), zStepFin = dimdelta($IntermedIntName, 1), zNpFin = dimsize($IntermedIntName, 1)
variable zEndFin = zStartFin + zStepFin*(zNpFin - 1)
killwaves/Z $IntermedIntName

//variable xStartCur = dimoffset($IntermedIntName1, 0), xStepCur = dimdelta($IntermedIntName1, 0), xNpCur = dimsize($IntermedIntName1, 0)
//variable xEndCur = xStartCur + (xNpCur - 1)*xStepCur
//variable zStartCur = dimoffset($IntermedIntName1, 1), zStepCur = dimdelta($IntermedIntName1, 1), zNpCur = dimsize($IntermedIntName1, 1)
//variable zEndCur = zStartCur + (zNpCur - 1)*zStepCur
//$FinIntName = $IntermedIntName1(x)(y)*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)

//xStartCur = dimoffset($IntermedIntName, 0); xStepCur = dimdelta($IntermedIntName, 0); xNpCur = dimsize($IntermedIntName, 0)
//xEndCur = xStartCur + (xNpCur - 1)*xStepCur
//zStartCur = dimoffset($IntermedIntName, 1); zStepCur = dimdelta($IntermedIntName, 1); zNpCur = dimsize($IntermedIntName, 1)
//zEndCur = zStartCur + (zNpCur - 1)*zStepCur
//$FinIntName += (0.5*RelSpecSensitivity)*($IntermedIntName(x)(y))*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)
variable xc = 0.5*(xStartFin + xEndFin), zc = 0.5*(zStartFin + zEndFin)

//$FinIntName *= IntMultConst

//if(Disp == 2)
display; appendimage $FinIntName
label bottom "Horizontal Position"
label left "Vertical Position"
SrwUtiGraphWindResize(2,10,190,170,-2,0)
//display; appendimage $IntermedIntName
//label bottom "Horizontal Position"
//label left "Vertical Position"
//SrwUtiGraphWindResize(200,10,190,170,-2,0)
	
make/O/N=(xNpFin) $FinIntNameX
setscale/P x xStartFin, xStepFin, "m", $FinIntNameX
$FinIntNameX = $FinIntName(x)(zc)
display $FinIntNameX
label bottom "Horizontal Position"
label left "W/mm\\S2\\M"
SrwUtiGraphAddFrameAndGrid()
SrwUtiGraphWindResize(2,242,190,150,-1,0)

make/O/N=(zNpFin) $FinIntNameZ
setscale/P x zStartFin, zStepFin, "m", $FinIntNameZ
$FinIntNameZ = $FinIntName(xc)(x)
display $FinIntNameZ
label bottom "Vertical Position"
label left "W/mm\\S2\\M"
SrwUtiGraphAddFrameAndGrid()
SrwUtiGraphWindResize(200,242,190,150,-1,0)
//endif

variable ie = 1
e_keV = eStart_keV + eStep_keV

do 
	SrwSmpScanXZE(AuxObsName+SrwSmpType, xc_mm, xRange_mm, xNp, zc_mm, zRange_mm, zNp, e_keV, e_keV, 1)
	SrwWfrCreate(AuxWfrName, ElecName, MagName, AuxObsName+SrwSmpType, 2, $nmPrecParam[1])
	
	//curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
	//resizeInitXD = minInitNpX/curNpX
	//if(resizeInitXD < 1.2)
	//	resizeInitXD = 1
	//endif
	//resizeInitZD = minInitNpZ/curNpZ
	//if(resizeInitZD < 1.2)
	//	resizeInitZD = 1
	//endif
	//if((resizeInitXD != 1) %| (resizeInitZD != 1))
	//	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
	//endif
	if((resizeInitRangeX != 1) %| (resizeInitResolX != 1) %| (resizeInitRangeZ != 1) %| (resizeInitResolZ != 1))
		SrwWfrResize(AuxWfrName+SrwRadType,1,resizeInitRangeX,resizeInitResolX,resizeInitRangeZ,resizeInitResolZ,1,"")
	endif
	
	//SrwWfrPropag(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 1, "")
	//SrwWfrProp(AuxWfrName+SrwRadType, BLName, 1, 1, Prec, 2, 1, "")
	SrwWfrProp(AuxWfrName+SrwRadType, BLName, 2, 2, 1, 2, 1, "")
	
	//curNpX = dimsize($nmRadCmpnX, 1); curNpZ = dimsize($nmRadCmpnX, 2)
	//resizeInitXD = maxPropNpX/curNpX
	//if(resizeInitXD > 0.8)
	//	resizeInitXD = 1
	//endif
	//resizeInitZD = maxPropNpZ/curNpZ
	//if(resizeInitZD > 0.8)
	//	resizeInitZD = 1
	//endif
	//if((resizeInitXD != 1) %| (resizeInitZD != 1))
	//	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeInitXD,1,resizeInitZD,1,"")
	//endif
	
	//xStepCur = dimdelta($nmRadCmpnX, 1)
	//zStepCur = dimdelta($nmRadCmpnX, 2)
	
	//resizeXD = abs(xStepCur/xStepFin)
	//resizeZD = abs(zStepCur/zStepFin)
	//if(resizeXD < 1.2)
	//	resizeXD = 1
	//endif
	//if(resizeZD < 1.2)
	//	resizeZD = 1
	//endif
	//if((resizeXD > 1) %| (resizeXD > 1))
	//	SrwWfrResize(AuxWfrName+SrwRadType,1,1,resizeXD,1,resizeZD,1,"")
	//endif
	
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 4, 1, e_keV, 0, 0, 1)
	RelSpecSensitivity = $SpecSensName(e_keV*1000)
	
	//if(Disp == 2)
	strLabel =  num2str(e_keV) + " keV"
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 2, 1, e_keV, 0, zc_mm, 1)
	strExe = "TextBox/W=" + strWinNameCutX +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	SrwWfr2Int(AuxWfrName+SrwRadType, "", PolCmpn, IntensCmpnNo, 3, 1, e_keV, xc_mm, 0, 1)
	strExe = "TextBox/W=" + strWinNameCutZ +"/C/N=text0/D=0.3/X=5/Y=9 \"" + strLabel + "\""; execute strExe
	$IntermedIntNameX *= RelSpecSensitivity
	$IntermedIntNameZ *= RelSpecSensitivity
	//endif
	
	AuxConst = IntMultConst*RelSpecSensitivity
	
	//xStartCur = dimoffset($nmRadCmpnX, 1); xStepCur = dimdelta($nmRadCmpnX, 1); xNpCur = dimsize($nmRadCmpnX, 1)
	//xEndCur = xStartCur + (xNpCur - 1)*xStepCur
	//zStartCur = dimoffset($nmRadCmpnX, 2); zStepCur = dimdelta($nmRadCmpnX, 2); zNpCur = dimsize($nmRadCmpnX, 2)
	//zEndCur = zStartCur + (zNpCur - 1)*zStepCur
	//$FinIntName += AuxConst*($IntermedIntName(x)(y))*srwUtiNonZeroInterval(x, xStartCur, xEndCur)*srwUtiNonZeroInterval(y, zStartCur, zEndCur)
	
	if((ebmPropSigX != 0) %| (ebmPropSigZ != 0))
		SrwUtiConvWaveWithGaus2D(IntermedIntName, ebmPropSigX, ebmPropSigZ)
	endif
	$FinIntName += AuxConst*$IntermedIntName(x)(y)
	killwaves/Z $IntermedIntName
	
	//if(Disp == 2)
	$FinIntNameX = $FinIntName(x)(zc)
	$FinIntNameZ = $FinIntName(xc)(x)
	//endif
	
		SaveExperiment
	
	e_keV += eStep_keV
	ie += 1
while(ie < (eNp - 1))

SrwWfrDel(AuxWfrName+SrwRadType)
killwaves/Z IntermedIntName //, IntermedIntName1
end

//+++++++++++++++++++++++++++++++++++++++
//
//Auxiliary proc: application of SrwWfrCreatePropagChrom to pinhole simulations
//Call Example:
//AuxUtiCalcMultiPinholeChrom("Ch_Cu_0_1_", 0.025, "wDimsPinZ", "SOLEIL_BM3-8DEG_comm_ebm", "Mag_mag", "ObsEXZ_obs", 2.2, 7, 1, "wSpecTr_Cu_0_1_dense")
//AuxUtiCalcMultiPinholeChrom("Ch_Cu_0_3_", 0.025, "wDimsPinZ", "SOLEIL_BM3-8DEG_comm_ebm", "Mag_mag", "ObsEXZ_obs", 2.6, 7, 1, "wSpecTr_Cu_0_3_dense")
//AuxUtiCalcMultiPinholeChrom("Ch_Cu_1a_", 0.025, "wDimsPinZa", "SOLEIL_BM3-8DEG_comm_ebm", "Mag_mag", "ObsEXZ_obs", 2.6, 7, 1, "wSpecTr_Cu_1_dense")
//+++++++++++++++++++++++++++++++++++++++
proc AuxUtiCalcMultiPinholeChrom(nmCoreResInt, dimPinX, nmDimsPinZ, nmElec, nmMag, nmObs, prec, componPol, componMultiE, nmSpecTr)
variable dimPinX, prec, componPol, componMultiE
string nmCoreResInt, nmDimsPinZ, nmElec, nmMag, nmObs, nmSpecTr

variable distAfterPin = 5.703 //[m] to edit !!!
variable numberPinDims = dimsize($nmDimsPinZ, 0)
SrwOptDrift("Drift", distAfterPin)
SrwOptApertRect("RectAperture", dimPinX, $nmDimsPinZ[0], 0, 0)
SrwOptContFull("PinholeBL", "RectAperture_bli", "Drift_bli", "", "")
string nmResInt, strSuffix

variable iDim = 0
do
	SrwOptApertRect("RectAperture", dimPinX, $nmDimsPinZ[iDim], 0, 0)
	strSuffix = num2str($nmDimsPinZ[iDim]) + "_mm"
		print "Current Pinhole size:", strSuffix
	nmResInt = nmCoreResInt + strSuffix
	SrwWfrCreatePropagChrom(nmResInt, nmElec, nmMag, nmObs, "PinholeBL_bli", prec, componPol, componMultiE, nmSpecTr, 2)

	iDim += 1
while(iDim < numberPinDims)
end