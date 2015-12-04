
//+++++++++++++++++++++++++++++++++++++++
//
//Compute Stokes from Per. Magn. Field
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPerStoCreate(RadName, ElecName, MagName, ObsName, InitHarm, FinHarm, ns, nphi, FluxType)
String RadName=srwUtiTruncString(SrwElecName+SrwUndName+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType
String MagName=SrwUndName+SrwUndType
String ObsName=SrwSmpName+SrwSmpType
Variable InitHarm=SrwPerInitHarm
Variable FinHarm=SrwPerFinHarm
Variable ns=SrwPerNs
Variable nphi=SrwPerNphi
Variable FluxType=SrwFluxType
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt InitHarm,SrwPPerInitHarm
prompt FinHarm,SrwPPerFinHarm
prompt ns,SrwPPerNs
prompt nphi,SrwPPerNphi
prompt FluxType,SrwPFluxType,popup "Photons/s/.1%bw/pixel;Photons/s/.1%bw/mm^2";
Silent 1						|	Computing the Radiation  ...
PauseUpdate

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
	SrwElecThick()
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwPerStoCreate()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	if(ObsWavePresent == 1)
		SrwPerStoCreate()
		Return
	endif
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 3
	SrwStartMacrosAfterRadSmp2 *= -1
	SrwRadSamplDialog(SrwSmpType)
	Return
endif

// Validation of parameters
if((abs(round(InitHarm) - InitHarm) > 1.E-08) %| (InitHarm <= 0))
	Abort SrwPAlertUndRadHarm
endif
if((abs(round(FinHarm) - FinHarm) > 1.E-08) %| (FinHarm <= 0))
	Abort SrwPAlertUndRadHarm
endif
if(ns <= 0.)
	Abort SrwPAlertUndPrecNs
endif
if(nphi <= 0.)
	Abort SrwPAlertUndPrecNphi
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]

SrwPerInitHarm=InitHarm
SrwPerFinHarm=FinHarm
SrwPerNs=ns
SrwPerNphi=nphi
SrwFluxType=FluxType

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwUndName+SrwSmpName
endif
SrwStoName=RadName

// Preparing data for C function
Make/D/O/N=6 waveprec
waveprec[0]=InitHarm   // Initial harmonic to take into account
waveprec[1]=FinHarm   // Final harmonic to take into account
waveprec[2]=ns   // Parameter of Longitudinal integration
waveprec[3]=nphi   // Parameter of Circular integration
waveprec[4]=FluxType

waveprec[5]=1
if(exists("SrwPerStoPrecEnExtRight") == 2)
	waveprec[5]=SrwPerStoPrecEnExtRight
endif

SrwStoPrep(ElecName,MagName,Obsname,RadName,FluxType)
RadName += SrwStoType
SrwRadGenTotName=RadName

srStokesUnd($ElecName, $MagName, $ObsName, waveprec, $RadName)

KillWaves/Z  waveprec
end 

//+++++++++++++++++++++++++++++++++++++++
//TEST version
//Compute Stokes from Per. Magn. Field (testing extra precision parameter, responsible for accuracy of convolution to take into account finite number of periods and energy spread)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwPerStokesCreate(RadName, ElecName, MagName, ObsName, InitHarm, FinHarm, ns, nphi, nConv, FluxType)
string RadName=srwUtiTruncString(SrwElecName+SrwUndName+SrwSmpName, 27)
string ElecName=SrwElecName+SrwElecType
string MagName=SrwUndName+SrwUndType
string ObsName=SrwSmpName+SrwSmpType
variable InitHarm=SrwPerInitHarm
variable FinHarm=SrwPerFinHarm
variable ns=SrwPerNs
variable nphi=SrwPerNphi
variable FluxType=SrwFluxType
variable nConv=srwUtiGetValN("SrwPerStoPrecEnExtRight", 1, "")
prompt RadName,SrwPStoName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "")
prompt InitHarm,SrwPPerInitHarm
prompt FinHarm,SrwPPerFinHarm
prompt ns,SrwPPerNs
prompt nphi,SrwPPerNphi
prompt FluxType,SrwPFluxType,popup "Photons/s/.1%bw/pixel;Photons/s/.1%bw/mm^2"
prompt nConv,"Length, En. Spread Convol. Param."

Silent 1						|	Computing the Radiation  ...
PauseUpdate

srwUtiSetValN("SrwPerStoPrecEnExtRight", nConv, "")
SrwPerStoCreate(RadName, ElecName, MagName, ObsName, InitHarm, FinHarm, ns, nphi, FluxType)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Prepare Stokes structure
//For the moment, Stokes are only a numerical wave
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoPrep(ElecName,MagName,ObsName,StoName,FluxType)
String ElecName,MagName,ObsName,StoName;
Variable FluxType;

StoName += SrwStoType;

variable ne = $ObsName[7], nx = $ObsName[10], nz = $ObsName[13];

variable eSt=($ObsName[5]), eFi=($ObsName[6]);
if((ne==1)%&(eSt==eFi)) 
	eFi = eSt*1.000001; // Can be anything
endIf

variable xSt=($ObsName[8]), xFi=($ObsName[9]);

if(nx==1)
	if(xSt != xFi)
		xSt = 0.5*(xSt + xFi);
	endif
	
	if(xSt==0)
		xFi = Abs(xFi - xSt)*1.e-08;
	else
		xFi = xSt*1.000001; // Can be anything
	endif
endIf

variable zSt=($ObsName[11]), zFi=($ObsName[12]);
if(nz==1)
	if(zSt != zFi)
		zSt = 0.5*(zSt + zFi);
	endif
	
	if(zSt==0)
		zFi = Abs(zFi - zSt)*1.e-08;
	else
		zFi = zSt*1.000001; // Can be anything
	endif
endIf

if(ne < 1)
	ne = 1;
endif
if(nx < 1)
	nx = 1;
endif
if(nz < 1)
	nz = 1;
endif

Make/O/N=(4, ne, nx, nz) $StoName;
SetScale/I y eSt, eFi,"eV", $StoName;
SetScale/I z xSt, xFi, "m", $StoName;
SetScale/I t zSt, zFi, "m", $StoName;
if(FluxType==1)
	SetScale d 0, 0, SrwPUnitSpAngFlux, $StoName;
endif
if(FluxType==2)
	SetScale d 0, 0, SrwPUnitSpAngFluxPerUnSurf, $StoName;
endif
	//print xSt, xFi

end

//+++++++++++++++++++++++++++++++++++++++
//
//Prepare Stokes structure
//For the moment, Stokes are only a numerical wave
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoPrepSimple(ObsName,StoName,FluxType)
String ObsName,StoName
Variable FluxType

StoName += SrwStoType;

variable ne = $ObsName[7], nx = $ObsName[10], nz = $ObsName[13];

variable eSt=($ObsName[5]), eFi=($ObsName[6]);
if((ne==1)%&(eSt==eFi)) 
	eFi = eSt*1.000001; // Can be anything
endIf

variable xSt=($ObsName[8]), xFi=($ObsName[9]);

if(nx==1)
	if(xSt != xFi)
		xSt = 0.5*(xSt + xFi);
	endif
	
	if(xSt==0)
		xFi = Abs(xFi - xSt)*1.e-08;
	else
		xFi = xSt*1.000001; // Can be anything
	endif
endIf

variable zSt=($ObsName[11]), zFi=($ObsName[12]);
if(nz==1)
	if(zSt != zFi)
		zSt = 0.5*(zSt + zFi);
	endif
	
	if(zSt==0)
		zFi = Abs(zFi - zSt)*1.e-08;
	else
		zFi = zSt*1.000001; // Can be anything
	endif
endIf

if(ne < 1)
	ne = 1;
endif
if(nx < 1)
	nx = 1;
endif
if(nz < 1)
	nz = 1;
endif

Make/O/N=(4, ne, nx, nz) $StoName;
SetScale/I y eSt, eFi,"eV", $StoName;
SetScale/I z xSt, xFi, "m", $StoName;
SetScale/I t zSt, zFi, "m", $StoName;
if(FluxType==1)
	SetScale d 0, 0, SrwPUnitSpAngFlux, $StoName;
endif
if(FluxType==2)
	SetScale d 0, 0, SrwPUnitSpAngFluxPerUnSurf, $StoName;
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Prepare Stokes structure directly.
//For the moment, Stokes are only a numerical wave
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoPrepSimpleStr(StoName, eSt, eFi, ne, xSt, xFi, nx, zSt, zFi, nz)
string StoName
variable eSt, eFi, ne, xSt, xFi, nx, zSt, zFi, nz

StoName += SrwStoType

if(ne < 1)
	ne = 1
endif
if(nx < 1)
	nx = 1
endif
if(nz < 1)
	nz = 1
endif

Make/O/N=(4, ne, nx, nz) $StoName
SetScale/I y eSt, eFi,"eV", $StoName
SetScale/I z xSt, xFi, "m", $StoName
SetScale/I t zSt, zFi, "m", $StoName
SetScale d 0, 0, SrwPUnitSpAngFluxPerUnSurf, $StoName
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract Intensity from Stokes Components
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSto2Int(StoName, SuffixExtract, RadCmpnType, PlotType, eVal, xVal, zVal, dis)
String StoName=SrwStoName+SrwStoType
String SuffixExtract=SrwSuffixExtract
Variable RadCmpnType=SrwViewRadCmpnType
Variable PlotType=SrwViewPlotType
Variable eVal=SrwViewE
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt StoName,SrwPStoName1, popup Wavelist("*"+SrwStoType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt PlotType, SrwPViewPlotType, popup "Energy;Hor.;Vert.;Hor. & Vert.;Energy & Hor.;Energy & Vert.;Energy & Hor. & Vert.;Auto"
prompt eVal, SrwPViewE
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,"New Display",popup "No;Yes"
Silent 1						|	  ...
PauseUpdate

SrwStoName=StoName[0,strlen(StoName)-strlen(SrwStoType)-1]
SrwSuffixExtract=SuffixExtract
SrwViewPlotType=PlotType
SrwViewRadCmpnType=RadCmpnType
SrwViewE=eVal
SrwViewX=xVal
SrwViewZ=zVal

if(cmpstr(StoName,"_none_")==0)
	//DoAlert 0, SrwPAlertNoCompResultsFound; Return;
	
	SrwStartMacrosAfterRadSmp2 = 3; // To proceed default chain through RadSampling panel
	SrwPerStoCreate();
	if(SrwStartMacrosAfterRadSmp2 > 0)
		SrwSto2Int();
		SrwStartMacrosAfterRadSmp2 = 0;
	endif
	Return;
endif

eVal *= 1000.
xVal *= 0.001
zVal *= 0.001

String st=StoName

// Treating "Auto"
if(PlotType==8)
	PlotType=srwSto2IntTreatAuto($st)
endif

String ViewStoName=SrwStoName+SuffixExtract
String DataUn=WaveUnits($st,-1)
	
if(PlotType==1) // Energy
	ViewStoName+=SrwSeparator+SrwRadEType
	Make/O/N=(DimSize($st, 1)) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, xVal, zVal)
	if(dis==2)
		if(DimSize($st, 1)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelPhotEn
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==2) // Hor
	ViewStoName+=SrwSeparator+SrwRadXType
	Make/O/N=(DimSize($st, 2)) $ViewStoName
	SetScale/P x DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, x, zVal)
	if(dis==2)
		if(DimSize($st, 2)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelHorPos
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==3) // Ver
	ViewStoName+=SrwSeparator+SrwRadZType
	Make/O/N=(DimSize($st, 3)) $ViewStoName
	SetScale/P x DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, xVal, x)
	if(dis==2)
		if(DimSize($st, 3)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelVerPos
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==4) // Hor & Ver
	ViewStoName+=SrwSeparator+SrwRadXType+SrwRadZType
	Make/O/N=((DimSize($st, 2)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	SetScale/P y DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, x, y)
	if(dis==2)
		if((DimSize($st, 2)>1) %| (DimSize($st, 3)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelHorPos
			Label left SrwPLabelVerPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==5) // Energy & Hor.
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadXType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 2))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	$ViewStoName=sRComSto($s1, RadCmpnType, x, y, zVal)
	if(dis==2)
		if((DimSize($st, 1)>1) %| (DimSize($st, 2)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelPhotEn
			Label left SrwPLabelHorPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==6) // Energy & Vert.
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadZType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, xVal, y)
	if(dis==2)
		if((DimSize($st, 1)>1) %| (DimSize($st, 3)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelPhotEn
			Label left SrwPLabelVerPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==7)
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadXType+SrwRadZType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 2)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	SetScale/P z DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, y, z)
endif

SrwUtiDataWaveInfStore(ViewStoName, "Unit", DataUn)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract Intensity from Stokes Components
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSto2IntF(StoName, SuffixExtract, RadCmpnType, IntOrFlux, PlotType, eVal, xVal, zVal, dis)
String StoName=SrwStoName+SrwStoType
String SuffixExtract=SrwSuffixExtract
Variable RadCmpnType=SrwViewRadCmpnType
variable IntOrFlux=srwUtiGetValN("IntOrFlux", 1, "SrwSto2IntF")
Variable PlotType=SrwViewPlotType
Variable eVal=SrwViewE
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt StoName,SrwPStoName1, popup Wavelist("*"+SrwStoType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt IntOrFlux, "Intensity or Flux", popup "Flux / Surface;Flux"
prompt PlotType, SrwPViewPlotType, popup "Energy;Hor.;Vert.;Hor. & Vert.;Energy & Hor.;Energy & Vert.;Energy & Hor. & Vert.;Auto"
prompt eVal, SrwPViewE
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,"New Display",popup "No;Yes"
Silent 1						|	  ...
PauseUpdate

SrwStoName=StoName[0,strlen(StoName)-strlen(SrwStoType)-1]
SrwSuffixExtract=SuffixExtract
SrwViewPlotType=PlotType
SrwViewRadCmpnType=RadCmpnType
srwUtiSetValN("IntOrFlux", IntOrFlux, "SrwSto2IntF")
SrwViewE=eVal
SrwViewX=xVal
SrwViewZ=zVal

if(cmpstr(StoName,"_none_")==0)
	//DoAlert 0, SrwPAlertNoCompResultsFound; Return;
	
	SrwStartMacrosAfterRadSmp2 = 3; // To proceed default chain through RadSampling panel
	SrwPerStoCreate();
	if(SrwStartMacrosAfterRadSmp2 > 0)
		SrwSto2Int()
		SrwStartMacrosAfterRadSmp2 = 0
	endif
	return
endif

eVal *= 1000.
xVal *= 0.001
zVal *= 0.001

string st=StoName

variable npe = dimsize($st, 1), npx = dimsize($st, 2), npz = dimsize($st, 3)
variable IntegrationIsNeeded = 0
if(IntOrFlux == 2)
	//if((npx > 0) %& (npz > 0))
	if((npx > 1) %& (npz > 1)) //OC050609
		IntegrationIsNeeded = 1
	endif
endif

variable ArgDisplay = dis
if(IntegrationIsNeeded == 1) //flux
	dis = 1 // switching off display temporarily
	if(npe <= 1)
		PlotType = 4
	else
		PlotType = 7
	endif
endif

// Treating "Auto"
if(PlotType==8)
	PlotType=srwSto2IntTreatAuto($st)
endif

string ViewStoName=SrwStoName+SuffixExtract
variable LenExtrCore = strlen(ViewStoName)
variable SuffLen = strlen(SuffixExtract)
if(LenExtrCore > 27)
	if(SuffLen > 10)
		abort "The suffix string is too long."
	endif
	SrwStoName = srwUtiTruncString(SrwStoName, 27 - SuffLen)
	ViewStoName=SrwStoName+SuffixExtract
endif

string CoreViewStoName = ViewStoName
string DataUn=WaveUnits($st,-1)
	
if(PlotType==1) // Energy
	ViewStoName+=SrwSeparator+SrwRadEType
	Make/O/N=(DimSize($st, 1)) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, xVal, zVal)
	if(dis==2)
		if(DimSize($st, 1)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelPhotEn
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==2) // Hor
	ViewStoName+=SrwSeparator+SrwRadXType
	Make/O/N=(DimSize($st, 2)) $ViewStoName
	SetScale/P x DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, x, zVal)
	if(dis==2)
		if(DimSize($st, 2)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelHorPos
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==3) // Ver
	ViewStoName+=SrwSeparator+SrwRadZType
	Make/O/N=(DimSize($st, 3)) $ViewStoName
	SetScale/P x DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, xVal, x)
	if(dis==2)
		if(DimSize($st, 3)>1)
			Display;Append $ViewStoName
			Label bottom SrwPLabelVerPos
			Label left DataUn
		else
			print $ViewStoName[0], DataUn
		endif
	endif
endif
if(PlotType==4) // Hor & Ver
	ViewStoName+=SrwSeparator+SrwRadXType+SrwRadZType
	Make/O/N=((DimSize($st, 2)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	SetScale/P y DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, eVal, x, y)
	if(dis==2)
		if((DimSize($st, 2)>1) %| (DimSize($st, 3)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelHorPos
			Label left SrwPLabelVerPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==5) // Energy & Hor.
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadXType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 2))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	$ViewStoName=sRComSto($s1, RadCmpnType, x, y, zVal)
	if(dis==2)
		if((DimSize($st, 1)>1) %| (DimSize($st, 2)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelPhotEn
			Label left SrwPLabelHorPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==6) // Energy & Vert.
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadZType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, xVal, y)
	if(dis==2)
		if((DimSize($st, 1)>1) %| (DimSize($st, 3)>1))
			Display;AppendImage $ViewStoName
			SrwImageFormat(ViewStoName)
			Label bottom SrwPLabelPhotEn
			Label left SrwPLabelVerPos
		else
			print $ViewStoName[0][0], DataUn
		endif
	endif
endif
if(PlotType==7)
	ViewStoName+=SrwSeparator+SrwRadEType+SrwRadXType+SrwRadZType
	Make/O/N=((DimSize($st, 1)), (DimSize($st, 2)), (DimSize($st, 3))) $ViewStoName
	SetScale/P x DimOffset($st, 1), DimDelta($st, 1), WaveUnits($st, 1), $ViewStoName
	SetScale/P y DimOffset($st, 2), DimDelta($st, 2), WaveUnits($st, 2), $ViewStoName
	SetScale/P z DimOffset($st, 3), DimDelta($st, 3), WaveUnits($st, 3), $ViewStoName
	$ViewStoName=sRComSto($st, RadCmpnType, x, y, z)
endif

SrwUtiDataWaveInfStore(ViewStoName, "Unit", DataUn)

//processing flux
if(IntegrationIsNeeded == 1) //flux

	string NewVeiwName = CoreViewStoName + SrwSeparator + SrwRadEType
	make/O/N=(npe) $NewVeiwName
	if(npe <= 1)
		$NewVeiwName[0] = srwUtiIntTotWave2D($ViewStoName)*(1e+6)
		if(ArgDisplay==2)
			print $NewVeiwName[0], SrwPUnitSpAngFlux
		endif
	else
		//make/O/N=(npe) auxwavesrwsto2intf
		setscale/P x dimoffset($ViewStoName, 0), dimdelta($ViewStoName, 0), SrwPUnitPhotEn, $NewVeiwName
		$NewVeiwName = srwUtiIntWave3Dvs2D($ViewStoName, p)*(1e+6)
		
		if(ArgDisplay==2)
			Display;Append $NewVeiwName
			Label bottom SrwPLabelPhotEn
			Label left SrwPUnitSpAngFlux
		endif
	endif
	SrwUtiDataWaveInfStore(NewVeiwName, "Unit", SrwPUnitSpAngFlux)

	killwaves/Z $ViewStoName

endif
end

//+++++++++++++++++++++++++++++++++++++++
function srwSto2IntTreatAuto(st)
wave st
if(DimSize(st, 1)==1)
	if(DimSize(st, 2)==1)
		return 3
	else
		if(DimSize(st, 3)==1)
			return 2
		else
			return 4
		endif
	endif
else
	if(DimSize(st, 2)==1)
		if(DimSize(st, 3)==1)
			return 1
		else
			return 6
		endif
	else
		if(DimSize(st, 3)==1)
			return 5
		else
			return 7
		endif
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
function sRComSto(st, cmp, eVal, xVal, zVal)
wave st
Variable eVal, xVal, zVal, cmp

if(cmp==7)
	return st[0](eVal)(xVal)(zVal)
endif
if(cmp==1)
	return 0.5*(st[0](eVal)(xVal)(zVal)+st[1](eVal)(xVal)(zVal))
endif
if(cmp==2)
	return 0.5*(st[0](eVal)(xVal)(zVal)-st[1](eVal)(xVal)(zVal))
endif
if(cmp==3)
	return 0.5*(st[0](eVal)(xVal)(zVal)+st[2](eVal)(xVal)(zVal))
endif
if(cmp==4)
	return 0.5*(st[0](eVal)(xVal)(zVal)-st[2](eVal)(xVal)(zVal))
endif
if(cmp==5)
	return 0.5*(st[0](eVal)(xVal)(zVal)+st[3](eVal)(xVal)(zVal))
endif
if(cmp==6)
	return  0.5*(st[0](eVal)(xVal)(zVal)-st[3](eVal)(xVal)(zVal))
endif

// etc...
return 0
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (polarization rate, calls SrwSto2Int)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSto2PolRate(StoName, SuffixExtract, RadCmpnType, PlotType, eVal, xVal, zVal, dis)
String StoName=SrwStoName+SrwStoType
String SuffixExtract=SrwSuffixExtract
Variable RadCmpnType=SrwViewRadCmpnType
Variable PlotType=SrwViewPlotType
Variable eVal=SrwViewE
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt StoName,SrwPStoName1, popup Wavelist("*"+SrwStoType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt PlotType, SrwPViewPlotType, popup "Energy;Hor.;Vert.;Hor. & Vert.;Energy & Hor.;Energy & Vert.;Energy & Hor. & Vert.;Auto"
prompt eVal, SrwPViewE
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,"New Display",popup "No;Yes"
Silent 1						|	  ...
PauseUpdate

if(RadCmpnType==7)
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

String ViewRadName=StoName[0,strlen(StoName)-strlen(SrwStoType)-1]+SuffixExtract
if(PlotType==8)
	PlotType=srwSto2IntTreatAuto($StoName)
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

SrwSto2Int(StoName, SuffixExtract, CmplmntType, PlotType, eVal, xVal, zVal, 1)
String OrtPolWave="AuxWaveOrtPol"
duplicate/O $ViewRadName $OrtPolWave

Variable ResIsOneVal=0
if((DimSize($ViewRadName,0)<=1) %& (DimSize($ViewRadName,1)<=1) %& (DimSize($ViewRadName,2)<=1))
	ResIsOneVal=1
endif
Variable DisInt=dis
if((ResIsOneVal != 0) %& (dis==2))
	DisInt=1
endif

SrwSto2Int(StoName, SuffixExtract, RadCmpnType, PlotType, eVal, xVal, zVal, DisInt)

WaveStats/Q $ViewRadName
Variable ZeroAbsTol = V_max*(1e-13) //V_max*(1e-04) // To prevent division by zero

$ViewRadName = Abs($ViewRadName[p][q][r])/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)
KillWaves/Z $OrtPolWave

if((dis==2) %& (PlotType<=3) %& (ResIsOneVal==0))
	Label left SrwPLabelPolRate
endif
if((ResIsOneVal != 0) %& (dis==2))
	print "Polarization rate: ", $ViewRadName[0][0][0]
endif

end

//+++++++++++++++++++++++++++++++++++++++
//
//Extract results (polarization rate, calls SrwSto2Int)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwSto2PolRateExt(StoName, SuffixExtract, RadCmpnType, IntOrFlux, PlotType, RateType, eVal, xVal, zVal, dis)
String StoName=SrwStoName+SrwStoType
String SuffixExtract=SrwSuffixExtract
Variable RadCmpnType=SrwViewRadCmpnType
Variable IntOrFlux=srwUtiGetValN("SrwCmpnIntOrFlux", 1, "")
Variable PlotType=SrwViewPlotType
Variable RateType=srwUtiGetValN("SrwPolRateType", 1, "")
Variable eVal=SrwViewE
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt StoName,SrwPStoName1, popup Wavelist("*"+SrwStoType, ";", "")
prompt SuffixExtract, SrwPViewSuffixExtract
prompt RadCmpnType, SrwPViewRadCmpnType, popup SrwPOPUPPolar+";Total"
prompt IntOrFlux,"Via Intensity or Flux",popup "Intensity;Flux"
prompt PlotType, SrwPViewPlotType, popup "Energy;Hor.;Vert.;Hor. & Vert.;Energy & Hor.;Energy & Vert.;Energy & Hor. & Vert.;Auto"
prompt RateType, "Normalisation", popup "I1/(I1+I2);(I1-I2)/(I1+I2)"
prompt eVal, SrwPViewE
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,"New Display",popup "No;Yes"
Silent 1						|	  ...
PauseUpdate

if(RadCmpnType==7)
	Abort SrwPAlertBadPolRate, CmpnNo
endif

srwUtiSetValN("SrwCmpnIntOrFlux", IntOrFlux, "")
srwUtiSetValN("SrwPolRateType", RateType, "")

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

string ViewRootRadName = StoName[0,strlen(StoName)-strlen(SrwStoType)-1]
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
	PlotType=srwSto2IntTreatAuto($StoName)
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

if(IntOrFlux < 1)
	IntOrFlux = 1
endif
if(IntOrFlux > 2)
	IntOrFlux = 2
endif

//SrwSto2Int(StoName, SuffixExtract, CmplmntType, PlotType, eVal, xVal, zVal, 1)
SrwSto2IntF(StoName, SuffixExtract, CmplmntType, IntOrFlux, PlotType, eVal, xVal, zVal, 1)

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

//if(IntOrFlux == 1)
//	SrwSto2Int(StoName, SuffixExtract, RadCmpnType, PlotType, eVal, xVal, zVal, DisInt)
//else
//	if(IntOrFlux == 2)
//		SrwSto2IntF(StoName, SuffixExtract, RadCmpnType, IntOrFlux, PlotType, eVal, xVal, zVal, DisInt)
//	endif
//endif
SrwSto2IntF(StoName, SuffixExtract, RadCmpnType, IntOrFlux, PlotType, eVal, xVal, zVal, DisInt)

WaveStats/Q $ViewRadName
Variable ZeroAbsTol = V_max*(1e-13) //V_max*(1e-04) // To prevent division by zero

if(RateType == 1)
	//$ViewRadName = Abs($ViewRadName[p][q][r])/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)
	$ViewRadName = Abs($ViewRadName[p][q][r])/srwAuxRetNonZero(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]), ZeroAbsTol)
else
	//$ViewRadName = (Abs($ViewRadName[p][q][r]) - Abs($OrtPolWave[p][q][r]))/(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]) + ZeroAbsTol)
	$ViewRadName = (Abs($ViewRadName[p][q][r]) - Abs($OrtPolWave[p][q][r]))/srwAuxRetNonZero(Abs($OrtPolWave[p][q][r]) + Abs($ViewRadName[p][q][r]), ZeroAbsTol)
endif
KillWaves/Z $OrtPolWave

if((dis==2) %& (PlotType<=3) %& (ResIsOneVal==0))
	Label left SrwPLabelPolRate
endif
if((ResIsOneVal != 0) %& (dis==2))
	print "Polarization rate: ", $ViewRadName[0][0][0]
endif

end

//+++++++++++++++++++++++++++++++++++++++
//
// Duplicate Stokes
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoDupl(StoName,DuplStoName)
string StoName=SrwStoName
string DuplStoName=SrwStoName+"d"
prompt StoName,SrwPStoName1, popup Wavelist("*"+SrwStoType,";","")
prompt DuplStoName,SrwPStoNameDpl

Silent 1						|	  ...
PauseUpdate

if(cmpstr(StoName,"_none_")==0)
	DoAlert 0, SrwPAlertNoCompResultsFound;
	Return;
endif

SrwStoName=DuplStoName;

String NewStoName=SrwStoName+SrwStoType
duplicate/O $StoName $NewStoName

end

//+++++++++++++++++++++++++++++++++++++++
//
//Deflecting parameter
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndK(MagField,Period)
Variable MagField=SrwMagBpeak
Variable Period=SrwPeriod
prompt MagField,"On-axis Peak Magnetic Field [T]"
prompt Period,"Period [mm]"
Silent 1						|	...
PauseUpdate

SrwPeriod=Period
SrwMagBpeak=MagField

SrwKK=srUtiUndK(MagField, Period*0.001)
SrwKK=srRound(SrwKK,5)

print SrwKK
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extracts total K value from Undulator structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiUndGetK(wUnd)
wave/T wUnd

variable numHarm = str2num(wUnd[5]), i, sumKxE2 = 0, sumKzE2 = 0, curK, curN
string nmHarm
for(i = 6; i < dimsize(wUnd, 0); i += 1)
	nmHarm = wUnd[i]
	if(cmpstr(nmHarm, "") == 0)
		continue
	endif
	wave wHarm = $nmHarm
	if(dimsize(wHarm, 0) < 3)
		continue
	endif
	
	curN = wHarm[0]
	curK = wHarm[2]
	if(wHarm[1] == 1)
		sumKzE2 += curK*curK/(curN*curN)
	endif
	if(wHarm[1] == 2)
		sumKxE2 += curK*curK/(curN*curN)
	endif
endfor
return sqrt(sumKxE2 + sumKzE2)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Multiplies K values for all harmonics by a given factor
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiUndScaleK(wUnd, facK)
wave/T wUnd
variable facK

variable numHarm = str2num(wUnd[5]), i
string nmHarm
for(i = 6; i < dimsize(wUnd, 0); i += 1)
	nmHarm = wUnd[i]
	if(cmpstr(nmHarm, "") == 0)
		continue
	endif
	wave wHarm = $nmHarm
	if(dimsize(wHarm, 0) < 3)
		continue
	endif
	wHarm[2] *= facK
endfor
end

//+++++++++++++++++++++++++++++++++++++++
//
//Multiplies K values for all harmonics by a given factor
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiUndScaleK1D(wUnd, facK, cmpn)
wave/T wUnd
variable facK
variable cmpn //1- vert., 2- hor.

variable numHarm = str2num(wUnd[5]), i
string nmHarm
for(i = 6; i < dimsize(wUnd, 0); i += 1)
	nmHarm = wUnd[i]
	if(cmpstr(nmHarm, "") == 0)
		continue
	endif
	wave wHarm = $nmHarm
	if(dimsize(wHarm, 0) < 3)
		continue
	endif
	
	if(wHarm[1] == cmpn)
		wHarm[2] *= facK
	endif	
endfor
end

//+++++++++++++++++++++++++++++++++++++++
//
//Extracts period value from Undulator structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiUndGetPer(wUnd)
wave/T wUnd
return str2num(wUnd[0])
end

//+++++++++++++++++++++++++++++++++++++++
//
//Fundamental Photon Energy
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndFundPhotEn(MagField,Period,ElecEnergy,outunit)
Variable MagField=SrwMagBpeak
Variable Period=SrwPeriod
Variable ElecEnergy=SrwElecEn
Variable outunit=1
prompt MagField,"On-axis Peak Magnetic Field [T]"
prompt Period,"Period [mm]"
prompt ElecEnergy,"Electron Energy [GeV]"
prompt outunit,"Output Units",popup "keV;eV;1/cm;Å;nm;µm;mm"
Silent 1						|	...
PauseUpdate

SrwPeriod=Period
SrwMagBpeak=MagField
SrwElecEn=ElecEnergy

String sout = ""
if(outunit == 1)
	sout = "keV"
endif
if(outunit == 2)
	sout = "eV"
endif
if(outunit == 3)
	sout = "1/cm"
endif
if(outunit == 4)
	sout = "Å"
endif
if(outunit == 5)
	sout = "nm"
endif
if(outunit == 6)
	sout = "µm"
endif
if(outunit == 7)
	sout = "mm"
endif

print srUtiUndFundPhotEn(MagField, Period*0.001, ElecEnergy, outunit), sout
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwUtiUndFundPhotEnFromK(K0,Period,ElecEnergy,outunit)
Variable K0=SrwKK
Variable Period=SrwPeriod
Variable ElecEnergy=SrwElecEn
Variable outunit=1
prompt K0,"Deflection Parameter"
prompt Period,"Period [mm]"
prompt ElecEnergy,"Electron Energy [GeV]"
prompt outunit,"Output Units",popup "keV;eV;1/cm;Å;nm;µm;mm"
Silent 1						|	...
PauseUpdate

SrwKK = K0

variable B0 = K0/(0.0933729*Period)
SrwUtiUndFundPhotEn(B0,Period,ElecEnergy,outunit)
end

//+++++++++++++++++++++++++++++++++++++++
function SrwUtiUndKfromFundPhotEn(PhotEn_eV,Period_mm,ElecEn_GeV)
variable PhotEn_eV,Period_mm,ElecEn_GeV
NVAR SrwKK
SrwKK = sqrt((9496.3421866853*ElecEn_GeV*ElecEn_GeV/Period_mm/PhotEn_eV - 1)*2)
return SrwKK
end

//+++++++++++++++++++++++++++++++++++++++
function SrwUtiUndB1fromFundPhotEn(PhotEn_eV,Period_mm,ElecEn_GeV)
variable PhotEn_eV,Period_mm,ElecEn_GeV
NVAR SrwKK
SrwKK = sqrt((9496.3421866853*ElecEn_GeV*ElecEn_GeV/Period_mm/PhotEn_eV - 1)*2)
return SrwKK/(0.0933729*Period_mm)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Harmonic number with on-axis resonant photon energy value closest to the given photon energy
//
//+++++++++++++++++++++++++++++++++++++++
function srwUtiClosestHarmNum(PhotEn_keV,K,Period_mm,ElecEnergy_GeV)
variable PhotEn_keV,K,Period_mm,ElecEnergy_GeV
variable FundPhotEn_keV = 0.94963421866853*ElecEnergy_GeV*ElecEnergy_GeV/(1 + 0.5*K*K)/(0.1*Period_mm)
variable Nh = round(PhotEn_keV/FundPhotEn_keV)
if(Nh <= 0) 
	Nh = 1 
endif
return Nh
end

//+++++++++++++++++++++++++++++++++++++++
//
//Calculates maximal spectral flux at harmonics 
//through a fixed aperture  vs K 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwUndPeakSpecFlux_vsK(nmSpec, nmElec, nmUndMaxK, nmObs, hBeg, hEnd, minK, numK, pol, prec)
string nmSpec=srwUtiTruncString(SrwElecName+SrwUndName+SrwSmpName, 22)
string nmElec=SrwElecName+SrwElecType
string nmUndMaxK=SrwUndName+SrwUndType
string nmObs=SrwSmpName+SrwSmpType
variable hBeg=srwUtiGetValN("hBeg", 1, "SrwUndPeakSpecFlux_vsK")
variable hEnd=srwUtiGetValN("hEnd", 3, "SrwUndPeakSpecFlux_vsK")
variable minK=srwUtiGetValN("minK", 0.03, "SrwUndPeakSpecFlux_vsK")
variable numK=srwUtiGetValN("numK", 3, "SrwUndPeakSpecFlux_vsK")
variable pol=srwUtiGetValN("pol", 1, "SrwUndPeakSpecFlux_vsK")
variable prec=srwUtiGetValN("prec", 1, "SrwUndPeakSpecFlux_vsK")
prompt nmSpec, "Core Name of the Spectrum to be calculated"
prompt nmElec, "Electron Beam structure", popup Wavelist("*"+SrwElecType ,";", "")
prompt nmUndMaxK, "Undulator structure (max. K)", popup Wavelist("*"+SrwUndType ,";", "")
prompt nmObs, "Radiation Sampling structure", popup Wavelist("*"+SrwSmpType ,";", "")
prompt hBeg, "Initial Harmonic Number"
prompt hEnd, "Final Harmonic Number"
prompt minK, "Minimal K value"
prompt numK, "Number of K values"
prompt pol, "Polarization", popup SrwPOPUPPolar+";Total"
prompt prec, "Precision parameter"
Silent 1						|	Calculating spectra ...
PauseUpdate

if((hBeg <= 0) %| (hBeg > hEnd))
	abort "Incorrect harmonic number(s)"
endif
if((minK <= 0) %| (numK <= 0))
	abort "Incorrect deflection parameter value(s)"
endif
if((pol < 1) %| (pol > 7))
	abort "Incorrect polarization type"
endif

srwUtiSetValN("hBeg", hBeg, "SrwUndPeakSpecFlux_vsK")
srwUtiSetValN("hEnd", hEnd, "SrwUndPeakSpecFlux_vsK")
srwUtiSetValN("minK", minK, "SrwUndPeakSpecFlux_vsK")
srwUtiSetValN("numK", numK, "SrwUndPeakSpecFlux_vsK")
srwUtiSetValN("pol", pol, "SrwUndPeakSpecFlux_vsK")
srwUtiSetValN("prec", prec, "SrwUndPeakSpecFlux_vsK")

variable typePolRate = 2 //popup "I1/(I1+I2);(I1-I2)/(I1+I2)"
variable numExtraHarm = 3
variable factRange = 0.2

variable maxK = srwUtiUndGetK($nmUndMaxK)
string auxCoreNameUnd = srwUtiTruncString("Aux" + nmUndMaxK[0,strlen(nmUndMaxK)-strlen(SrwUndType)-1], 25)
string auxNameUnd = auxCoreNameUnd + SrwUndType

variable stepK = 0
if(numK > 1)
	stepK = (maxK - minK)/(numK - 1)
endif

string nmResPhotEn, nmResSpec, nmResPolRate
string nmHarmCore
variable ih = hBeg, ih0
do
	nmHarmCore = nmSpec + num2str(ih)
	nmResSpec = nmHarmCore + ".F"
	nmResPhotEn = nmHarmCore + ".E"
	make/O/N=(numK) $nmResPhotEn, $nmResSpec
	
	if(pol != 7)
		nmResPolRate = nmHarmCore + ".R"
		make/O/N=(numK) $nmResPolRate
	endif
	ih += 1
while(ih <= hEnd)

string nmStoAuxCore = "Aux" + nmSpec
string nmStoAux = nmStoAuxCore + SrwStoType
string sufFlux = "f", sufPolRate = "r"
string nmFluxAux = nmStoAuxCore + sufFlux + "_e"
string nmPolRateAux = nmStoAuxCore + sufPolRate + "_e"

variable per_m =  srwUtiUndGetPer($nmUndMaxK)
variable elEnGeV = srwGetElecBeamEnergy(nmElec)
variable e1_eV, eRange
variable amOfHarm =  hEnd - hBeg + 1

make/O/N=(amOfHarm, 3) wAuxPeakData

variable curK = maxK, facK, ik = 0
do
	SrwMagPerDupl(nmUndMaxK, auxCoreNameUnd)
	facK = curK/maxK
	srwUtiUndScaleK($auxNameUnd, facK)
	
	SrwPerStoCreate(nmStoAuxCore, nmElec, auxNameUnd, nmObs, hBeg, hEnd + numExtraHarm, prec, prec, 1)
	SrwSto2IntF(nmStoAux, sufFlux, 7, 1, 1, 1., 0., 0., 0)
	
	e1_eV = 9.50*elEnGeV*elEnGeV/(per_m*(1 + 0.5*curK*curK))
	
	eRange = factRange*e1_eV	
	
	wAuxPeakData = 0
	srwUtiFindSpecPeaks($nmFluxAux, e1_eV, eRange, hBeg, hEnd, wAuxPeakData)
	if(pol != 7)
		SrwSto2PolRateExt(nmStoAux, sufPolRate, pol, 1, 1, typePolRate, 1., 0., 0., 0)
	endif
	
	ih = hBeg
	ih0 = 0
	do
		nmHarmCore = nmSpec + num2str(ih)
		nmResSpec = nmHarmCore + ".F"
		nmResPhotEn = nmHarmCore + ".E"
		
		$nmResPhotEn[ik] = wAuxPeakData[ih0][1]
		$nmResSpec[ik] = wAuxPeakData[ih0][2]
	
		if(pol != 7)
			nmResPolRate = nmHarmCore + ".R"		
			$nmResPolRate[ik] = $nmPolRateAux($nmResPhotEn[ik])	
		endif

		ih0 += 1
		ih += 1
	while(ih <= hEnd)

	ik += 1
	curK -= stepK
while(ik < numK)

display
ih = hBeg
do
	nmHarmCore = nmSpec + num2str(ih)
	nmResSpec = nmHarmCore + ".F"
	nmResPhotEn = nmHarmCore + ".E"
	AppendToGraph $nmResSpec vs $nmResPhotEn
	ih += 1
while (ih <= hEnd)

if(pol != 7)
	display
	ih = hBeg
	do
		nmHarmCore = nmSpec + num2str(ih)
		nmResPolRate = nmHarmCore + ".R"
		nmResPhotEn = nmHarmCore + ".E"
		AppendToGraph $nmResPolRate vs $nmResPhotEn
		ih += 1
	while (ih <= hEnd)
endif

killwaves/Z $nmStoAux, wAuxPeakData
end