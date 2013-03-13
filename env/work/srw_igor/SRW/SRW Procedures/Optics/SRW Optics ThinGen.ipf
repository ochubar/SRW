
//+++++++++++++++++++++++++++++++++++++++
//
//Optical Element: Thin Generic 
//
//+++++++++++++++++++++++++++++++++++++++
//
//Create Template with empty complex transmission wave
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinTempl(ElemName,TransmOuter,xc,zc,xr,zr,nx,nz)
String ElemName=SrwBliThinGen
Variable TransmOuter=SrwBliThinGenOuter
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
Variable xr=SrwBliThinGenWaveRangeX
Variable zr=SrwBliThinGenWaveRangeZ
Variable nx=SrwBliThinGenWaveNx
Variable nz=SrwBliThinGenWaveNz
prompt ElemName,SrwPBli
prompt TransmOuter,SrwPBliThinGenOuter,popup SrwPOPUPBliThinGenOuter
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
prompt xr,SrwPBliThinGenWaveRangeX
prompt zr,SrwPBliThinGenWaveRangeZ
prompt nx,SrwPBliThinGenWaveNx
prompt nz,SrwPBliThinGenWaveNz
Silent 1						|	 ...
PauseUpdate

// Validation
if((nx <= 0) %| (nz <= 0))
	Abort SrwPAlertBliThinGenWaveNpts
endif
if((xr <= 0) %| (zr <= 0))
	Abort SrwPAlertBliThinGenWaveRange
endif

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliThinGenWavePosX=xc
SrwBliThinGenWavePosZ=zc
SrwBliThinGenWaveRangeX=xr
SrwBliThinGenWaveRangeZ=zr
SrwBliThinGenWaveNx=nx
SrwBliThinGenWaveNz=nz
SrwBliThinGenOuter=TransmOuter

xc*=0.001
zc*=0.001
xr*=0.001
zr*=0.001

String CWaveName=ElemName+SrwOptThinGenWaveType
killwaves/Z $CWaveName
Make/C/D/O/N=(nx,nz) $CWaveName
SetScale/I x (xc-0.5*xr),(xc+0.5*xr),"m",$CWaveName
SetScale/I y (zc-0.5*zr),(zc+0.5*zr),"m",$CWaveName

ElemName+=SrwBeamlineType
Make/T/O/N=11 $ElemName
$ElemName[0]=SrwBliThinGenType
$ElemName[1]=CWaveName

$ElemName[4]=num2str(xc)
$ElemName[5]=num2str(zc)
$ElemName[6]=num2str(TransmOuter)
$ElemName[7]="0" // Setup was finished or not
$ElemName[8]="" //8 - foc. dist. x
$ElemName[9]="" //9 - foc. dist. z
$ElemName[10]="1" //"1"- opt. path; "2"- ph. shift
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create Template with empty complex transmission wave
//Obsolete
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenTempl(ElemName,en,xc,zc,xr,zr,nx,nz,TransmOuter)
String ElemName=SrwBliThinGen
Variable en=SrwBliThinGenE
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
Variable xr=SrwBliThinGenWaveRangeX
Variable zr=SrwBliThinGenWaveRangeZ
Variable nx=SrwBliThinGenWaveNx
Variable nz=SrwBliThinGenWaveNz
Variable TransmOuter=SrwBliThinGenOuter
prompt ElemName,SrwPBli
prompt en,SrwPBliThinGenE
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
prompt xr,SrwPBliThinGenWaveRangeX
prompt zr,SrwPBliThinGenWaveRangeZ
prompt nx,SrwPBliThinGenWaveNx
prompt nz,SrwPBliThinGenWaveNz
prompt TransmOuter,SrwPBliThinGenOuter,popup SrwPOPUPBliThinGenOuter
Silent 1						|	 ...
PauseUpdate

// Validation
if((nx <= 0) %| (nz <= 0))
	Abort SrwPAlertBliThinGenWaveNpts
endif
if((xr <= 0) %| (zr <= 0))
	Abort SrwPAlertBliThinGenWaveRange
endif

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliThinGenE=en
SrwBliThinGenWavePosX=xc
SrwBliThinGenWavePosZ=zc
SrwBliThinGenWaveRangeX=xr
SrwBliThinGenWaveRangeZ=zr
SrwBliThinGenWaveNx=nx
SrwBliThinGenWaveNz=nz
SrwBliThinGenOuter=TransmOuter

en*=1000.
xc*=0.001
zc*=0.001
xr*=0.001
zr*=0.001

String CWaveName=ElemName+SrwOptThinGenWaveType
Make/C/D/O/N=(nx,nz) $CWaveName
SetScale/I x (xc-0.5*xr),(xc+0.5*xr),"m",$CWaveName
SetScale/I y (zc-0.5*zr),(zc+0.5*zr),"m",$CWaveName

ElemName+=SrwBeamlineType
Make/T/O/N=11 $ElemName
$ElemName[0]=SrwBliThinGenType
$ElemName[1]=CWaveName

$ElemName[3]=num2str(en)
$ElemName[4]=num2str(xc)
$ElemName[5]=num2str(zc)
$ElemName[6]=num2str(TransmOuter)
$ElemName[7]="0" // Setup was finished or not
$ElemName[8]="" //8 - foc. dist. x
$ElemName[9]="" //9 - foc. dist. z
//$ElemName[10]="1" //"1"- opt. path; "2"- ph. shift
end

//+++++++++++++++++++++++++++++++++++++++
// Amplitude Transmission and Optical Path Difference Proto-function
//+++++++++++++++++++++++++++++++++++++++
function AmpTrOptPathDifProtoFunc(xx, zz)
variable xx, zz
return 1
end

//+++++++++++++++++++++++++++++++++++++++
//Function to setup 2D complex transmission wave and optical component, if exists. 
//Uses Optical Path
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThinSetupFunc(CWaveName, FunTransName, FunOptPathName)
string CWaveName, FunTransName, FunOptPathName
FUNCREF AmpTrOptPathDifProtoFunc fAmpTransm = $FunTransName
FUNCREF AmpTrOptPathDifProtoFunc fOptPathDif = $FunOptPathName
wave/C wCTr = $CWaveName
wCTr = cmplx(fAmpTransm(x, y), fOptPathDif(x, y))
end

//+++++++++++++++++++++++++++++++++++++++
//
//Procedure Setup 2D complex transmission wave and optical component, if exists. 
//Uses Optical Path
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinSetup(CWaveName,FunTransName,FunOptPathName)
string CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
string FunTransName=SrwBliThinGenWaveFunTrans
string FunOptPathName=SrwBliThinGenWaveFunOptPath
prompt CWaveName,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt FunTransName,SrwPBliThinGenWaveFunTrans,popup FunctionList("*",";","KIND:2,NPARAMS:2,VALTYPE:1")
prompt FunOptPathName,SrwPBliThinGenWaveFunOptPath,popup FunctionList("*",";","KIND:2,NPARAMS:2,VALTYPE:1")
Silent 1						|	Setting up complex transmission ...
PauseUpdate

variable FunTransExistInfo=exists(FunTransName)
variable FunOptPathExistInfo=exists(FunOptPathName)
if(((FunTransExistInfo != 3) %& (FunTransExistInfo != 6)) %| ((FunOptPathExistInfo != 3) %& (FunOptPathExistInfo != 6)))
	Abort SrwPAlertBliThinGenFunAbsent
endif

SrwOptThinGenValidateCmplWave(CWaveName)

SrwBliThinGen=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]
SrwBliThinGenWaveFunTrans=FunTransName
SrwBliThinGenWaveFunOptPath=FunOptPathName

//String ComLineStr
//sprintf ComLineStr, "'%s'=cmplx(%s(x,y),%s(x,y))", CWaveName, FunTransName, FunOptPathName
//Execute ComLineStr
SrwOptThinSetupFunc(CWaveName, FunTransName, FunOptPathName)

string ElemName=SrwBliThinGen+SrwBeamlineType
variable ElemExistInfo=exists(ElemName)
if(ElemExistInfo == 1)
	$ElemName[10]="1" //"1"- opt. path; "2"- ph. shift
	//srOptThinGenSetup($ElemName) // Finish Setup in C (analize focal distances, etc.)
	srOptElemSetup($ElemName) // Finish Setup in C (analize focal distances, etc.)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete: uses Phase Shift function
//Setup 2D complex transmission wave and optical component, if exists
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenSetup(CWaveName,FunTransName,FunPhaseName)
String CWaveName=SrwBliThinGen
String FunTransName=SrwBliThinGenWaveFunTrans
String FunPhaseName=SrwBliThinGenWaveFunPhase
prompt CWaveName,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt FunTransName,SrwPBliThinGenWaveFunTrans
prompt FunPhaseName,SrwPBliThinGenWaveFunPhase
Silent 1						|	Setting up complex transmission ...
PauseUpdate

Variable FunTransExistInfo=exists(FunTransName)
Variable FunPhaseExistInfo=exists(FunPhaseName)
if(((FunTransExistInfo != 3) %& (FunTransExistInfo != 6)) %| ((FunPhaseExistInfo != 3) %& (FunPhaseExistInfo != 6)))
	Abort SrwPAlertBliThinGenFunAbsent
endif

SrwOptThinGenValidateCmplWave(CWaveName)

SrwBliThinGen=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]
SrwBliThinGenWaveFunTrans=FunTransName
SrwBliThinGenWaveFunPhase=FunPhaseName

String ComLineStr
sprintf ComLineStr, "'%s'=cmplx(%s(x,y),%s(x,y))", CWaveName, FunTransName, FunPhaseName
Execute ComLineStr

String ElemName=SrwBliThinGen+SrwBeamlineType
Variable ElemExistInfo=exists(ElemName)
if(ElemExistInfo == 1)
	$ElemName[10]="2" //"1"- opt. path; "2"- ph. shift
	//srOptThinGenSetup($ElemName) // Finish Setup in C (analize focal lengths, etc.)
	srOptElemSetup($ElemName) // Finish Setup in C (analize focal lengths, etc.)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Formal validation of 2D Complex Transmission wave
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenValidateCmplWave(CWaveName)
String CWaveName

Variable ExistanceInfo = exists(CWaveName)
if(ExistanceInfo != 1)
	Abort SrwPAlertBliThinGenWaveAbsent
endif

Variable CWaveInf = WaveType($CWaveName)
Variable WaveIsComplex = CWaveInf %& 0x01
if(WaveIsComplex == 0)
	Abort SrwPAlertBliThinGenImproperWave
endif

Variable WaveIs32BitFloat = CWaveInf %& 0x02
Variable WaveIs64BitFloat = CWaveInf %& 0x04
if((WaveIs32BitFloat == 0) %& (WaveIs64BitFloat == 0))
	Abort SrwPAlertBliThinGenImproperWave
endif

Variable AmOfDims = WaveDims($CWaveName)
if(AmOfDims != 2)
	Abort SrwPAlertBliThinGenImproperWave
endif

if(WaveIs32BitFloat != 0)
	Redimension/C/D $CWaveName
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Display Thin Generic Opt. Component
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinTransmDisplay(CWave,DispWave,Cmpn,PlotType,xVal,zVal,dis)
String CWave=SrwBliThinGen+SrwOptThinGenWaveType
String DispWave=SrwBliThinGenDefDispName
Variable Cmpn=SrwBliThinGenViewCmpn
Variable PlotType=SrwBliThinGenViewPlotType
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=SrwBliThinGenViewNewDisp
prompt CWave,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt DispWave, SrwPThinGenViewWave
prompt Cmpn,SrwPBliThinGenViewCmpn,popup SrwPOPUPBliThinViewCmpn
prompt PlotType,SrwPBliThinGenViewPlotType,popup SrwPOPUPBliThinGenViewPlotType
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,SrwPBliThinGenViewNewDisp,popup SrwPOPUPBliThinGenViewNewDisp
Silent 1						|	Extracting data ...
PauseUpdate

SrwOptThinGenValidateCmplWave(CWave)

SrwBliThinGen=CWave[0,strlen(CWave)-strlen(SrwOptThinGenWaveType)-1]
SrwBliThinGenViewCmpn=Cmpn
SrwBliThinGenViewPlotType=PlotType
SrwBliThinGenViewNewDisp=dis
SrwViewX=xVal
SrwViewZ=zVal
SrwBliThinGenDefDispName=DispWave

xVal *= 0.001
zVal *= 0.001

String LabelString = " "

Variable ModifyNameToEnsureUnique = 0
if(dis == 2)  // if 
	ModifyNameToEnsureUnique = 1 
endif

if(PlotType==1)
	DispWave+=SrwSeparator+SrwRadXType+SrwRadZType
	//DispWave = SrwUtiGiveNewName(DispWave, SrwSeparator+SrwRadXType+SrwRadZType, ModifyNameToEnsureUnique)
	Make/O/D/N=((DimSize($CWave, 0)), (DimSize($CWave, 1))) $DispWave
	SetScale/P x DimOffset($CWave, 0), DimDelta($CWave, 0), WaveUnits($CWave, 0), $DispWave
	SetScale/P y DimOffset($CWave, 1), DimDelta($CWave, 1), WaveUnits($CWave, 1), $DispWave
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave[p][q])
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave[p][q])^2
	endif
	if(Cmpn==3) // Optical Path
		$DispWave=imag($CWave[p][q])
	endif
endif
if(PlotType==2)
	DispWave+=SrwSeparator+SrwRadXType
	//DispWave = SrwUtiGiveNewName(DispWave, SrwSeparator+SrwRadXType, ModifyNameToEnsureUnique)
	Make/O/D/N=(DimSize($CWave, 0)) $DispWave
	SetScale/P x DimOffset($CWave, 0), DimDelta($CWave, 0), WaveUnits($CWave, 0), $DispWave
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave(x)(zVal))
		LabelString = SrwPLabelTransmAmp
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave(x)(zVal))^2
		LabelString = SrwPLabelTransmInt
	endif
	if(Cmpn==3) // Optical Path
		$DispWave=imag($CWave(x)(zVal))
		LabelString = SrwPLabelOptPath
	endif
endif
if(PlotType==3)
	DispWave+=SrwSeparator+SrwRadZType
	//DispWave = SrwUtiGiveNewName(DispWave, SrwSeparator+SrwRadZType, ModifyNameToEnsureUnique)
	Make/O/D/N=(DimSize($CWave, 1)) $DispWave
	SetScale/P x DimOffset($CWave, 1), DimDelta($CWave, 1), WaveUnits($CWave, 1), $DispWave
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave(xVal)(x))
		LabelString = SrwPLabelTransmAmp
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave(xVal)(x))^2
		LabelString = SrwPLabelTransmInt
	endif
	if(Cmpn==3) // Optical Path
		$DispWave=imag($CWave(xVal)(x))
		LabelString = SrwPLabelOptPath
	endif
endif

if(Cmpn==3) // Optical Path
	SetScale d 0, 1.e+23, "m", $DispWave
endif

if(dis==2)
	Display;
	if(PlotType==1)
		AppendImage $DispWave
		SrwImageFormat(DispWave)
		Label bottom SrwPLabelHorPos
		Label left SrwPLabelVerPos
	endif
	if(PlotType==2)
		Append $DispWave
		Label bottom SrwPLabelHorPos
		Label left LabelString
	endif
	if(PlotType==3)
		Append $DispWave
		Label bottom SrwPLabelVerPos
		Label left LabelString
	endif
endif
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete
//Display Thin Generic Opt. Component
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenDisplay(CWave,DispWave,Cmpn,PlotType,xVal,zVal,dis)
String CWave=SrwBliThinGen+SrwOptThinGenWaveType
String DispWave=SrwBliThinGen
Variable Cmpn=SrwBliThinGenViewCmpn
Variable PlotType=SrwBliThinGenViewPlotType
Variable xVal=SrwViewX
Variable zVal=SrwViewZ
Variable dis=1
prompt CWave,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt DispWave, SrwPThinGenViewWave
prompt Cmpn,SrwPBliThinGenViewCmpn,popup SrwPOPUPBliThinGenViewCmpn
prompt PlotType,SrwPBliThinGenViewPlotType,popup SrwPOPUPBliThinGenViewPlotType
prompt xVal, SrwPViewX
prompt zVal, SrwPViewZ
prompt dis,SrwPBliThinGenViewNewDisp,popup SrwPOPUPBliThinGenViewNewDisp
Silent 1						|	Extracting data ...
PauseUpdate

SrwOptThinGenValidateCmplWave(CWave)

SrwBliThinGen=CWave[0,strlen(CWave)-strlen(SrwOptThinGenWaveType)-1]
SrwBliThinGenViewCmpn=Cmpn
SrwBliThinGenViewPlotType=PlotType
SrwViewX=xVal
SrwViewZ=zVal

xVal *= 0.001;
zVal *= 0.001;

String LabelString = " "

if(PlotType==1)
	DispWave+=SrwSeparator+SrwRadXType+SrwRadZType
	Make/O/D/N=((DimSize($CWave, 0)), (DimSize($CWave, 1))) $DispWave
	SetScale/P x DimOffset($CWave, 0), DimDelta($CWave, 0), WaveUnits($CWave, 0), $DispWave
	SetScale/P y DimOffset($CWave, 1), DimDelta($CWave, 1), WaveUnits($CWave, 1), $DispWave;
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave[p][q])
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave[p][q])^2
	endif
	if(Cmpn==3) // Phase Shift
		$DispWave=imag($CWave[p][q])
	endif
	if(Cmpn==4) // Re(T)
		$DispWave=real(p2rect($CWave[p][q]))
	endif
	if(Cmpn==5) // Im(T)
		$DispWave=imag(p2rect($CWave[p][q]))
	endif
	
endif
if(PlotType==2)
	DispWave+=SrwSeparator+SrwRadXType
	Make/O/D/N=(DimSize($CWave, 0)) $DispWave
	SetScale/P x DimOffset($CWave, 0), DimDelta($CWave, 0), WaveUnits($CWave, 0), $DispWave
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave(x)(zVal))
		LabelString = SrwPLabelTransmAmp
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave(x)(zVal))^2
		LabelString = SrwPLabelTransmInt
	endif
	if(Cmpn==3) // Phase Shift
		$DispWave=imag($CWave(x)(zVal))
		LabelString = SrwPLabelPhaseShift
	endif
	if(Cmpn==4) // Re(T)
		$DispWave=real(p2rect($CWave(x)(zVal)))
		LabelString = SrwPLabelTransmCmplxRe
	endif
	if(Cmpn==5) // Im(T)
		$DispWave=imag(p2rect($CWave(x)(zVal)))
		LabelString = SrwPLabelTransmCmplxIm
	endif
	
endif
if(PlotType==3)
	DispWave+=SrwSeparator+SrwRadZType
	Make/O/D/N=(DimSize($CWave, 1)) $DispWave
	SetScale/P x DimOffset($CWave, 1), DimDelta($CWave, 1), WaveUnits($CWave, 1), $DispWave
	
	if(Cmpn==1) // Amp. Transm.
		$DispWave=real($CWave(xVal)(x))
		LabelString = SrwPLabelTransmAmp
	endif
	if(Cmpn==2) // Int. Transm.
		$DispWave=real($CWave(xVal)(x))^2
		LabelString = SrwPLabelTransmInt
	endif
	if(Cmpn==3) // Phase Shift
		$DispWave=imag($CWave(xVal)(x))
		LabelString = SrwPLabelPhaseShift
	endif
	if(Cmpn==4) // Re(T)
		$DispWave=real(p2rect($CWave(xVal)(x)))
		LabelString = SrwPLabelTransmCmplxRe
	endif
	if(Cmpn==5) // Im(T)
		$DispWave=imag(p2rect($CWave(xVal)(x)))
		LabelString = SrwPLabelTransmCmplxIm
	endif
	
endif

if(dis==2)
	Display;
	if(PlotType==1)
		AppendImage $DispWave
		SrwImageFormat(DispWave)
		Label bottom SrwPLabelHorPos
		Label left SrwPLabelVerPos
	endif
	if(PlotType==2)
		Append $DispWave
		Label bottom SrwPLabelHorPos
		Label left LabelString
	endif
	if(PlotType==3)
		Append $DispWave
		Label bottom SrwPLabelVerPos
		Label left LabelString
	endif

endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Estimate index of refraction (1 - n)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptRefrDelta(Z,Dens,PhotEn)
Variable Z=SrwBliMatZ
Variable Dens=SrwBliMatDens
Variable PhotEn=SrwBliThinGenE
prompt Z, SrwPBliMatZ
prompt Dens, SrwPBliMatDens
prompt PhotEn,SrwPBliThinGenE

SrwBliMatZ=Z
SrwBliMatDens=Dens
SrwBliThinGenE=PhotEn

PhotEn*=1000. // in eV

Variable RefrDelta = real(srOptMatConst(-1,Z,Dens,PhotEn))
if(SrwAllowPrintingExtraInfo==1)
	print "1 - n =", RefrDelta
endif

SrwBliZonePlateDeltaRefr=str2num(num2str(RefrDelta)) // to remove extra digits
SrwBliXrayLensDelta=SrwBliZonePlateDeltaRefr
end

//+++++++++++++++++++++++++++++++++++++++
//
//Estimate index of refraction (1 - n), Function version
//
//+++++++++++++++++++++++++++++++++++++++
function SrwOptRefrDeltaFun(Z,A,Dens,PhotEn)
Variable Z,A,Dens,PhotEn

PhotEn*=1000. // in eV

Variable Re = 2.817938070E-13 // class. el. rad. in cm
Variable Na = 6.02204531E+23
Variable TwoPI = 6.2831853071796
Variable Con = Re*Na/TwoPI
Variable Wavelength_cm = 1.239854E-04/PhotEn

return Z*Con*Wavelength_cm*Wavelength_cm*Dens/A
end

//+++++++++++++++++++++++++++++++++++++++
//
//Refractive X-ray Lens. Circular Holes. 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinXrayLensCirc(ElemName,FocPlane,Delta,AttenLen,Diam,Nhol,WallThick,xc,zc)
String ElemName=SrwBliXrayLens
Variable FocPlane=SrwBliXrayLensFocPlane
Variable Delta=SrwBliXrayLensDelta
Variable AttenLen=SrwBliXrayLensAttenLen
Variable Diam=SrwBliXrayLensDiam
Variable Nhol=SrwBliXrayLensNhol
Variable WallThick=SrwBliXrayLensWallThick
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
prompt ElemName,SrwPBli
prompt FocPlane,SrwPBliXrayLensFocPlane,popup SrwPOPUPBliXrayLensFocPlane
prompt Delta,"Refractive Index Decrement (= 1 - n)"
prompt AttenLen,"Intensity Attenuation Length [mm]"
prompt Diam,SrwPBliXrayLensDiam
prompt Nhol,SrwPBliXrayLensNhol
prompt WallThick,SrwPBliXrayLensWallThick
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
Silent 1						|	Setting up complex transmission ...
PauseUpdate

Variable FocDis_m=(0.25*Diam/(Nhol*Delta))*0.001
if(SrwAllowPrintingExtraInfo==1)
	print "Focal Distance:", FocDis_m, "m"
endif

Variable xr = Diam*1.1, zr = Diam*1.1 // default ranges of approximating mesh
Variable nx = 300, nz = 300 // default number of points in approximating mesh
Variable TransmOuter = 2 // outer transmission is the same as at the border
SrwOptThinTempl(ElemName,TransmOuter,xc,zc,xr,zr,nx,nz)

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliXrayLens=ElemName

ElemName=ElemName+SrwBeamlineType
$ElemName[7]="1" // Setup was finished
if(FocPlane==1) // Hor.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(1e+023)
endif
if(FocPlane==2) // Ver.
	$ElemName[8]=num2str(1e+023); $ElemName[9]=num2str(FocDis_m)
endif
if(FocPlane==3) // Hor. + Ver.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(FocDis_m)
endif

SrwBliXrayLensAttenLen=AttenLen
SrwBliXrayLensFocPlane=FocPlane
SrwBliXrayLensDiam=Diam
SrwBliXrayLensNhol=Nhol
SrwBliXrayLensWallThick=WallThick
SrwBliXrayLensDelta=Delta
//SrwBliXrayLensXc=xc
//SrwBliXrayLensZc=zc
SrwBliXrayLensXc=xc*0.001 //OC080311
SrwBliXrayLensZc=zc*0.001 //OC080311

String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
SrwOptThinSetup(CWaveName,"AmpTransm_XrayLensCirc","OptPath_XrayLensCirc")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Refractive X-ray Lens. Parabolic Holes.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinXrayLensParab(ElemName,FocPlane,Delta,AttenLen,Diam,Rmin,Nhol,WallThick,xc,zc)
String ElemName=SrwBliXrayLens
Variable FocPlane=SrwBliXrayLensFocPlane
Variable Delta=SrwBliXrayLensDelta
Variable AttenLen=SrwBliXrayLensAttenLen
Variable Diam=SrwBliXrayLensDiam
Variable Rmin=SrwBliXrayLensRmin
Variable Nhol=SrwBliXrayLensNhol
Variable WallThick=SrwBliXrayLensWallThick
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
prompt ElemName,SrwPBli
prompt FocPlane,SrwPBliXrayLensFocPlane,popup SrwPOPUPBliXrayLensFocPlane
prompt Delta,"Refractive Index Decrement (= 1 - n)"
prompt AttenLen,"Intensity Attenuation Length [mm]"
prompt Diam,"Geometrical Aperture [mm]"
prompt Rmin,"Minimal Radius of Parabola [mm]"
prompt Nhol,SrwPBliXrayLensNhol
prompt WallThick,SrwPBliXrayLensWallThick
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
Silent 1						|	Setting up complex transmission ...
PauseUpdate

Variable FocDis_m=(0.5*Rmin/(Nhol*Delta))*0.001
if(SrwAllowPrintingExtraInfo==1)
	print "Focal Length:", FocDis_m, "m"
endif

Variable xr = Diam*1.1, zr = Diam*1.1 // default ranges of approximating mesh
Variable nx = 800, nz = 800 // default number of points in approximating mesh
Variable TransmOuter = 2 // outer transmission is the same as at the border
SrwOptThinTempl(ElemName,TransmOuter,xc,zc,xr,zr,nx,nz)

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliXrayLens=ElemName

ElemName=ElemName+SrwBeamlineType
$ElemName[7]="1" // Setup was finished
if(FocPlane==1) // Hor.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(1e+023)
endif
if(FocPlane==2) // Ver.
	$ElemName[8]=num2str(1e+023); $ElemName[9]=num2str(FocDis_m)
endif
if(FocPlane==3) // Hor. + Ver.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(FocDis_m)
endif

SrwBliXrayLensAttenLen=AttenLen
SrwBliXrayLensFocPlane=FocPlane
SrwBliXrayLensDiam=Diam
SrwBliXrayLensRmin=Rmin
SrwBliXrayLensNhol=Nhol
SrwBliXrayLensWallThick=WallThick
SrwBliXrayLensDelta=Delta
//SrwBliXrayLensXc=xc
//SrwBliXrayLensZc=zc
SrwBliXrayLensXc=xc*0.001 //OC080311
SrwBliXrayLensZc=zc*0.001 //OC080311

String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
SrwOptThinSetup(CWaveName,"AmpTransm_XrayLensParab","OptPath_XrayLensParab")
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete version
//Refractive X-ray Lens. Circular Holes
//Obsolete version
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenSetupXrayLensCirc(CWaveName,Mat,AtomNum,AttenLen,Density,FocPlane,Diam,Nhol,WallThick)
String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
Variable Mat=SrwBliXrayLensMat
Variable AtomNum=SrwBliXrayLensAtomNum
Variable AttenLen=SrwBliXrayLensAttenLen
Variable Density=SrwBliXrayLensDensity
Variable FocPlane=SrwBliXrayLensFocPlane
Variable Diam=SrwBliXrayLensDiam
Variable Nhol=SrwBliXrayLensNhol
Variable WallThick=SrwBliXrayLensWallThick
prompt CWaveName,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt Mat,SrwPBliXrayLensMat,popup SrwPOPUPBliXrayLensMat
prompt AtomNum,SrwPBliXrayLensAtomNum
prompt AttenLen,SrwPBliXrayLensAttenLen
prompt Density,SrwPBliXrayLensDensity
prompt FocPlane,SrwPBliXrayLensFocPlane,popup SrwPOPUPBliXrayLensFocPlane
prompt Diam,SrwPBliXrayLensDiam
prompt Nhol,SrwPBliXrayLensNhol
prompt WallThick,SrwPBliXrayLensWallThick
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwOptThinGenValidateCmplWave(CWaveName)

if(Mat==1) // Be
	AtomNum=4; Density=1.845
endif
if(Mat==2) // C
	AtomNum=6; Density=2.20
endif // Continue for any more elements supported

String ElemName=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]+SrwBeamlineType
Variable PhotEn=str2num($ElemName[3])
Variable/C MatConst=srOptMatConst(Mat,AtomNum,Density,PhotEn)
Variable AuxDelta=str2num(num2str(real(MatConst))), AuxAttenLen=str2num(num2str(imag(MatConst)))
if(AuxAttenLen > 0.)
	AttenLen=AuxAttenLen
endif

Variable FocDis_m=(0.25*Diam/(Nhol*AuxDelta))*0.001

if(SrwAllowPrintingExtraInfo==1)
	print "Focal Distance:", FocDis_m, "m"
endif

$ElemName[7]="1" // Setup was finished
if(FocPlane==1) // Hor.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(1e+023)
endif
if(FocPlane==2) // Ver.
	$ElemName[8]=num2str(1e+023); $ElemName[9]=num2str(FocDis_m)
endif

SrwBliThinGen=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]
SrwBliXrayLensMat=Mat
SrwBliXrayLensAtomNum=AtomNum
SrwBliXrayLensAttenLen=AttenLen
SrwBliXrayLensDensity=Density
SrwBliXrayLensFocPlane=FocPlane
SrwBliXrayLensDiam=Diam
SrwBliXrayLensNhol=Nhol
SrwBliXrayLensWallThick=WallThick
SrwBliXrayLensDelta=AuxDelta
SrwBliXrayLensPhotEn=PhotEn
SrwBliXrayLensXc=str2num($ElemName[4])
SrwBliXrayLensZc=str2num($ElemName[5])

SrwOptThinGenSetup(CWaveName,"AmpTransm_XrayLensCirc","PhaseShift_XrayLensCirc")
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete version
//Refractive X-ray Lens. Parabolic Holes.
//Obsolete version
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenSetupXrayLensParab(CWaveName,Mat,AtomNum,AttenLen,Density,FocPlane,Diam,Rmin,Nhol,WallThick)
String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
Variable Mat=SrwBliXrayLensMat
Variable AtomNum=SrwBliXrayLensAtomNum
Variable AttenLen=SrwBliXrayLensAttenLen
Variable Density=SrwBliXrayLensDensity
Variable FocPlane=SrwBliXrayLensFocPlane
Variable Diam=SrwBliXrayLensDiam
Variable Rmin=SrwBliXrayLensRmin
Variable Nhol=SrwBliXrayLensNhol
Variable WallThick=SrwBliXrayLensWallThick
prompt CWaveName,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt Mat,SrwPBliXrayLensMat,popup SrwPOPUPBliXrayLensMat
prompt AtomNum,SrwPBliXrayLensAtomNum
prompt AttenLen,SrwPBliXrayLensAttenLen
prompt Density,SrwPBliXrayLensDensity
prompt FocPlane,SrwPBliXrayLensFocPlane,popup SrwPOPUPBliXrayLensFocPlane
prompt Diam,SrwPBliXrayLensDiam1
prompt Rmin,SrwPBliXrayLensRmin
prompt Nhol,SrwPBliXrayLensNhol
prompt WallThick,SrwPBliXrayLensWallThick
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwOptThinGenValidateCmplWave(CWaveName)

if(Mat==1) // Be
	AtomNum=4; Density=1.845
endif
if(Mat==2) // C
	AtomNum=6; Density=2.20
endif // Continue for any more elements supported

String ElemName=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]+SrwBeamlineType
Variable PhotEn=str2num($ElemName[3])
Variable/C MatConst=srOptMatConst(Mat,AtomNum,Density,PhotEn);
Variable AuxDelta=str2num(num2str(real(MatConst))), AuxAttenLen=str2num(num2str(imag(MatConst)))
if(AuxAttenLen > 0.)
	AttenLen=AuxAttenLen
endif

Variable FocDis_m=(0.5*Rmin/(Nhol*AuxDelta))*0.001

if(SrwAllowPrintingExtraInfo==1)
	print "Focal Distance:", FocDis_m, "m"
endif

$ElemName[7]="1" // Setup was finished
if(FocPlane==1) // Hor.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(1e+023)
endif
if(FocPlane==2) // Ver.
	$ElemName[8]=num2str(1e+023); $ElemName[9]=num2str(FocDis_m)
endif
if(FocPlane==3) // Hor. + Ver.
	$ElemName[8]=num2str(FocDis_m); $ElemName[9]=num2str(FocDis_m)
endif

SrwBliThinGen=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]
SrwBliXrayLensMat=Mat
SrwBliXrayLensAtomNum=AtomNum
SrwBliXrayLensAttenLen=AttenLen
SrwBliXrayLensDensity=Density
SrwBliXrayLensFocPlane=FocPlane
SrwBliXrayLensDiam=Diam
SrwBliXrayLensRmin=Rmin
SrwBliXrayLensNhol=Nhol
SrwBliXrayLensWallThick=WallThick
SrwBliXrayLensDelta=AuxDelta
SrwBliXrayLensPhotEn=PhotEn
SrwBliXrayLensXc=str2num($ElemName[4])
SrwBliXrayLensZc=str2num($ElemName[5])

SrwOptThinGenSetup(CWaveName,"AmpTransm_XrayLensParab","PhaseShift_XrayLensParab")
end

//+++++++++++++++++++++++++++++++++++++++
//Validation of X-ray lens parameters
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenValidateXrayLens(CWaveName,Mat,AtomNum,AttenLen,Density,FocPlane,Diam,Nhol,WallThick)
String CWaveName
Variable Mat, AtomNum, AttenLen, Density, FocPlane, Diam, Nhol, WallThick

Variable xStep = DimDelta($CWaveName, 0), zStep = DimDelta($CWaveName, 1)
Variable xStart = DimOffset($CWaveName, 0), zStart = DimOffset($CWaveName, 1)
Variable nx = DimSize($CWaveName, 0), nz = DimSize($CWaveName, 1)
Variable xEnd = xStart + xStep*(nx - 1), zEnd = zStart + zStep*(nz - 1)

// To fill, if any ...

end

//+++++++++++++++++++++++++++++++++++++++
// Amplitude Transmission function for X-ray lens with circular holes.
// Returns the Amplitude Transmission for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function AmpTransm_XrayLensCirc(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensAttenLen
Variable/G SrwBliXrayLensXc, SrwBliXrayLensZc

Variable r = 0.5*SrwBliXrayLensDiam*0.001
Variable SectLen = (SrwBliXrayLensWallThick*0.001 + 2*r)
Variable AttenLen_m = SrwBliXrayLensAttenLen*0.001
Variable PathInBodyPerSect = 0.
Variable u=0, w=0, ue2=0
Variable re2=r*r

if(SrwBliXrayLensFocPlane==1)
	u=xx-SrwBliXrayLensXc
endif
if(SrwBliXrayLensFocPlane==2)
	w=zz-SrwBliXrayLensZc
endif
if(SrwBliXrayLensFocPlane==3)
	u=xx-SrwBliXrayLensXc
	w=zz-SrwBliXrayLensZc
endif
ue2 = u*u + w*w

if(ue2 >= re2)
	PathInBodyPerSect = SectLen
else
	PathInBodyPerSect = SectLen - 2.*sqrt(re2 - ue2)
endif

return exp(-0.5*PathInBodyPerSect*SrwBliXrayLensNhol/AttenLen_m)
end

//+++++++++++++++++++++++++++++++++++++++
// Amplitude Transmission function for X-ray lens with parabolic holes.
// Returns the Amplitude Transmission for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function AmpTransm_XrayLensParab(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensRmin, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensAttenLen
Variable/G SrwBliXrayLensXc, SrwBliXrayLensZc

Variable HalfAp = SrwBliXrayLensDiam*0.0005 // m
Variable a = 1000./SrwBliXrayLensRmin
Variable d = SrwBliXrayLensWallThick*0.001
Variable HalfApe2 = HalfAp*HalfAp
Variable SectLen = d + a*HalfApe2
Variable AttenLen_m = SrwBliXrayLensAttenLen*0.001
Variable PathInBodyPerSect = 0.
Variable u = xx - SrwBliXrayLensXc
Variable xe2, ze2, re2

if(SrwBliXrayLensFocPlane==3) //Foc. in two planes
	xe2 = u*u
	u = zz - SrwBliXrayLensZc
	ze2 = u*u
	re2 = xe2 + ze2
	
	if(re2 > HalfApe2)
		PathInBodyPerSect = SectLen;
	else
		PathInBodyPerSect = a*re2 + d;
	endif
else
	if(SrwBliXrayLensFocPlane==2)
		u=zz-SrwBliXrayLensZc
	endif

	if((u < -HalfAp) %| (u > HalfAp))
		PathInBodyPerSect = SectLen;
	else
		PathInBodyPerSect = a*u*u + d;
	endif
endif

return exp(-0.5*PathInBodyPerSect*SrwBliXrayLensNhol/AttenLen_m)
end

//+++++++++++++++++++++++++++++++++++++++
// Optical Path function for X-ray lens with circular holes.
// Returns the Optical Path in m for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinSetup.
//+++++++++++++++++++++++++++++++++++++++
function OptPath_XrayLensCirc(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensDelta
Variable/G SrwBliXrayLensXc, SrwBliXrayLensZc

Variable r = 0.5*SrwBliXrayLensDiam*0.001
Variable SectLen = (SrwBliXrayLensWallThick*0.001 + 2*r)
Variable re2=r*r

Variable PathInBodyPerSect = 0
Variable u=0, w=0, ue2=0
if(SrwBliXrayLensFocPlane==1)
	u=xx-SrwBliXrayLensXc
	
	if((u < -r) %| (u > r))
		PathInBodyPerSect = SectLen
	else
		ue2 = u*u
		PathInBodyPerSect = SectLen - 2*sqrt(re2 - ue2)
	endif
endif
if(SrwBliXrayLensFocPlane==2)
	//w=zz-SrwBliXrayLensZc
	u=zz-SrwBliXrayLensZc
	
	if((u < -r) %| (u > r))
		PathInBodyPerSect = SectLen
	else
		ue2 = u*u
		PathInBodyPerSect = SectLen - 2*sqrt(re2 - ue2)
	endif
endif
if(SrwBliXrayLensFocPlane==3) //Foc. in two planes
	u=xx-SrwBliXrayLensXc
	w=zz-SrwBliXrayLensZc
	ue2 = u*u + w*w

	if(ue2 > r)
		PathInBodyPerSect = SectLen
	else
		PathInBodyPerSect = SectLen - 2*sqrt(re2 - ue2)
	endif
endif
//ue2 = u*u + w*w

//Variable OptPathPerSect = 0.
//if(ue2 < re2)
//	OptPathPerSect = 2.*SrwBliXrayLensDelta*sqrt(re2 - ue2)
//endif
//return OptPathPerSect*SrwBliXrayLensNhol
return -SrwBliXrayLensDelta*SrwBliXrayLensNhol*PathInBodyPerSect //OC200310
end

//+++++++++++++++++++++++++++++++++++++++
// Optical Path function for X-ray lens with parabolic holes.
// Returns the Optical Path im m for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinSetup.
//+++++++++++++++++++++++++++++++++++++++
function OptPath_XrayLensParab(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensRmin, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensDelta
Variable/G SrwBliXrayLensXc, SrwBliXrayLensZc

Variable HalfAp = SrwBliXrayLensDiam*0.0005 // m
Variable a = 1000./SrwBliXrayLensRmin
Variable d = SrwBliXrayLensWallThick*0.001
Variable HalfApe2 = HalfAp*HalfAp
Variable SectLen = d + a*HalfApe2

Variable PathInBodyPerSect = 0.
Variable u = xx - SrwBliXrayLensXc
Variable xe2, ze2, re2

if(SrwBliXrayLensFocPlane==3) //Foc. in two planes
	xe2 = u*u
	u = zz - SrwBliXrayLensZc
	ze2 = u*u
	re2 = xe2 + ze2
	
	if(re2 > HalfApe2)
		PathInBodyPerSect = SectLen
	else
		PathInBodyPerSect = a*re2 + d
	endif
else
	if(SrwBliXrayLensFocPlane==2)
		u=zz-SrwBliXrayLensZc
	endif

	if((u < -HalfAp) %| (u > HalfAp))
		PathInBodyPerSect = SectLen
	else
		PathInBodyPerSect = a*u*u + d
	endif
endif

//return -SrwBliXrayLensDelta*SrwBliXrayLensNhol*(PathInBodyPerSect-SectLen)
return -SrwBliXrayLensDelta*SrwBliXrayLensNhol*PathInBodyPerSect //OC200310
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete
// Phase Shift function for X-ray lens with circular holes.
// Returns the Phase Shift for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function PhaseShift_XrayLensCirc(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensDelta
Variable/G SrwBliXrayLensPhotEn, SrwBliXrayLensXc, SrwBliXrayLensZc

Variable r = 0.5*SrwBliXrayLensDiam*0.001
Variable SectLen = (SrwBliXrayLensWallThick*0.001 + 2*r)
Variable re2=r*r

Variable u=0, w=0, ue2=0
if(SrwBliXrayLensFocPlane==1)
	u=xx-SrwBliXrayLensXc
endif
if(SrwBliXrayLensFocPlane==2)
	w=zz-SrwBliXrayLensZc
endif
if(SrwBliXrayLensFocPlane==3)
	u=xx-SrwBliXrayLensXc
	w=zz-SrwBliXrayLensZc
endif
ue2 = u*u + w*w

Variable OptPathPerSect = 0.
if(ue2 < re2)
	OptPathPerSect = 2.*SrwBliXrayLensDelta*sqrt(re2 - ue2)
endif

Variable WaveNumb = (5.067681604e+06)*SrwBliXrayLensPhotEn // in 1/m, assuming PhotEn in eV
return WaveNumb*OptPathPerSect*SrwBliXrayLensNhol
end

//+++++++++++++++++++++++++++++++++++++++
//Obselete
// Phase Shift function for X-ray lens with parabolic holes.
// Returns the Phase Shift for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function PhaseShift_XrayLensParab(xx, zz)
Variable xx, zz

Variable/G SrwBliXrayLensFocPlane, SrwBliXrayLensDiam, SrwBliXrayLensRmin, SrwBliXrayLensWallThick, SrwBliXrayLensNhol, SrwBliXrayLensDelta
Variable/G SrwBliXrayLensPhotEn, SrwBliXrayLensXc, SrwBliXrayLensZc

Variable HalfAp = SrwBliXrayLensDiam*0.0005 // m
Variable a = 1000./SrwBliXrayLensRmin
Variable d = SrwBliXrayLensWallThick*0.001
Variable HalfApe2 = HalfAp*HalfAp
Variable SectLen = d + a*HalfApe2

//Variable u=xx-SrwBliXrayLensXc
//Variable OptPathPerSect = 0.
//
//if(SrwBliXrayLensFocPlane==2)
//	u=zz-SrwBliXrayLensZc
//endif
//
//if((u >= -HalfAp) %& (u <= HalfAp))
//	OptPathPerSect = SrwBliXrayLensDelta*(SectLen - (a*u*u + d))
//endif

Variable PathInBodyPerSect = 0.
Variable u = xx - SrwBliXrayLensXc
Variable xe2, ze2, re2

if(SrwBliXrayLensFocPlane==3) //Foc. in two planes
	xe2 = u*u
	u = zz - SrwBliXrayLensZc
	ze2 = u*u
	re2 = xe2 + ze2
	
	if(re2 > HalfApe2)
		PathInBodyPerSect = SectLen;
	else
		PathInBodyPerSect = a*re2 + d;
	endif
else
	if(SrwBliXrayLensFocPlane==2)
		u=zz-SrwBliXrayLensZc
	endif

	if((u < -HalfAp) %| (u > HalfAp))
		PathInBodyPerSect = SectLen;
	else
		PathInBodyPerSect = a*u*u + d;
	endif
endif

Variable WaveNumb = (5.067681604e+06)*SrwBliXrayLensPhotEn // in 1/m, assuming PhotEn in eV
//return WaveNumb*OptPathPerSect*SrwBliXrayLensNhol
return -WaveNumb*SrwBliXrayLensDelta*SrwBliXrayLensNhol*PathInBodyPerSect

end

//+++++++++++++++++++++++++++++++++++++++
//
//Simple Circular Zone Plate
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinSetupZonePlateSmpl(ElemName,Nzones,Rn,AttenLen,DeltaRefr,Thick,xc,zc,NperZone)
String ElemName=SrwBliZonePlate
Variable Rn=SrwBliZonePlateRn
Variable Nzones=SrwBliZonePlateNzones
Variable AttenLen=SrwBliZonePlateAttenLen1
Variable DeltaRefr=SrwBliZonePlateDeltaRefr1
Variable Thick=SrwBliZonePlateThick
Variable NperZone=5
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
prompt ElemName,SrwPBli
prompt Rn,"Auter Zone Radius [mm]"
prompt Nzones,"Number of Zones"
prompt NperZone,"Number of Mesh Points per Outer Zone"
prompt AttenLen,"Attenuation Length [mm]"
prompt DeltaRefr,"Delta (= 1 - n)"
prompt Thick,"Thickness [mm]"
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
Silent 1						|	Setting up complex transmission ...
PauseUpdate

//SrwBliZonePlatePhaseMult = 1
//if(DelPhi == 2)
//	SrwBliZonePlatePhaseMult = 3
//endif
//if(DelPhi == 3)
//	SrwBliZonePlatePhaseMult = 5
//endif

//String ElemName=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]+SrwBeamlineType
//Variable PhotEn=str2num($ElemName[3])
//SrwBliZonePlatePhotEn=PhotEn

//SrwBliZonePlateThick=Thick
//if(IgnThick==1)
//	Variable Lamb_mm = 1.239854e-03/PhotEn // if PhotEn in eV
//	SrwBliZonePlateThick = 0.5*Lamb_mm*SrwBliZonePlatePhaseMult/abs(DeltaRefr1-DeltaRefr2)
//	
//	if(SrwAllowPrintingExtraInfo==1)
//		print "Optimal Thickness:", SrwBliZonePlateThick, "mm"
//	endif
//endif
//SrwOptThinGenValidateZonePlate(CWaveName,FocDist,AttenLen1,DeltaRefr1,AttenLen2,DeltaRefr2,IgnThick,Thick,SrwBliZonePlatePhaseMult,Nzones)

Variable AuxPtsPerOuterZone = NperZone //to steer
Variable AuxRange = 2*Rn*1.1
Variable dRn = 0.5*Rn/Nzones
Variable Np = AuxRange*AuxPtsPerOuterZone/dRn
Variable TransmOuter = 2 // outer transmission is the same as at the border
SrwOptThinTempl(ElemName,TransmOuter,xc,zc,AuxRange,AuxRange,Np,Np)

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliZonePlate=ElemName

ElemName=ElemName+SrwBeamlineType
$ElemName[7]="1" // Setup was finished
$ElemName[8]=num2str(1.)
$ElemName[9]=num2str(1.)

SrwBliZonePlateThick=Thick
SrwBliZonePlateAttenLen1=AttenLen
SrwBliZonePlateDeltaRefr1=DeltaRefr
SrwBliZonePlateNzones=Nzones
SrwBliZonePlateRn=Rn
SrwBliZonePlateXc=str2num($ElemName[4])
SrwBliZonePlateZc=str2num($ElemName[5])

String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
SrwOptThinSetup(CWaveName,"AmpTransm_ZonePlateSmpl","OptPath_ZonePlateSmpl")
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete; use SrwOptThinSetupZonePlateSmpl
//Simple Circular Zone Plate
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenSetupZonePlateSmpl(CWaveName,FocDist,AttenLen1,DeltaRefr1,AttenLen2,DeltaRefr2,IgnThick,Thick,DelPhi,Nzones)
String CWaveName=SrwBliThinGen+SrwOptThinGenWaveType
Variable FocDist=SrwBliZonePlateFocDist
Variable AttenLen1=SrwBliZonePlateAttenLen1
Variable DeltaRefr1=SrwBliZonePlateDeltaRefr1
Variable AttenLen2=SrwBliZonePlateAttenLen2
Variable DeltaRefr2=SrwBliZonePlateDeltaRefr2
Variable IgnThick=SrwBliZonePlateIgnThick
Variable Thick=SrwBliZonePlateThick
Variable DelPhi=SrwBliZonePlateDelPhi
Variable Nzones=SrwBliZonePlateNzones
prompt CWaveName,SrwPBliThinGenWave1,popup Wavelist("*"+SrwOptThinGenWaveType,";","")
prompt FocDist,SrwPBliZonePlateFocDist
prompt AttenLen1,SrwPBliZonePlateAttenLen1
prompt DeltaRefr1,SrwPBliZonePlateDeltaRefr1
prompt AttenLen2,SrwPBliZonePlateAttenLen2
prompt DeltaRefr2,SrwPBliZonePlateDeltaRefr2
prompt IgnThick,SrwPBliZonePlateIgnThick,popup SrwPOPUPBliZonePlateIgnThick
prompt Thick,SrwPBliZonePlateThick
prompt DelPhi,SrwPBliZonePlateDelPhi,popup SrwPOPUPBliZonePlateDelPhi
prompt Nzones,SrwPBliZonePlateNzones
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwOptThinGenValidateCmplWave(CWaveName)

SrwBliZonePlatePhaseMult = 1
if(DelPhi == 2)
	SrwBliZonePlatePhaseMult = 3
endif
if(DelPhi == 3)
	SrwBliZonePlatePhaseMult = 5
endif

String ElemName=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]+SrwBeamlineType
Variable PhotEn=str2num($ElemName[3])
SrwBliZonePlatePhotEn=PhotEn

SrwBliZonePlateThick=Thick
if(IgnThick==1)
	Variable Lamb_mm = 1.239854e-03/PhotEn // if PhotEn in eV
	SrwBliZonePlateThick = 0.5*Lamb_mm*SrwBliZonePlatePhaseMult/abs(DeltaRefr1-DeltaRefr2)
	
	if(SrwAllowPrintingExtraInfo==1)
		print "Optimal Thickness:", SrwBliZonePlateThick, "mm"
	endif
endif

SrwOptThinGenValidateZonePlate(CWaveName,FocDist,AttenLen1,DeltaRefr1,AttenLen2,DeltaRefr2,IgnThick,Thick,SrwBliZonePlatePhaseMult,Nzones)

$ElemName[7]="1" // Setup was finished
$ElemName[8]=num2str(FocDist)
$ElemName[9]=num2str(FocDist)

SrwBliThinGen=CWaveName[0,strlen(CWaveName)-strlen(SrwOptThinGenWaveType)-1]
SrwBliZonePlateFocDist=FocDist
SrwBliZonePlateAttenLen1=AttenLen1
SrwBliZonePlateDeltaRefr1=DeltaRefr1
SrwBliZonePlateAttenLen2=AttenLen2
SrwBliZonePlateDeltaRefr2=DeltaRefr2

SrwBliZonePlateIgnThick=IgnThick
SrwBliZonePlateDelPhi=DelPhi
SrwBliZonePlateNzones=Nzones
SrwBliZonePlateXc=str2num($ElemName[4])
SrwBliZonePlateZc=str2num($ElemName[5])

SrwOptThinGenSetup(CWaveName,"AmpTransm_ZonePlateSmpl","PhaseShift_ZonePlateSmpl")
end

//+++++++++++++++++++++++++++++++++++++++
//Obsolete
//Validation of Zone Plate parameters
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinGenValidateZonePlate(CWaveName,FocDist,AttenLen1,DeltaRefr1,AttenLen2,DeltaRefr2,IgnThick,Thick,PhaseMult,Nzones)
String CWaveName
Variable FocDist,AttenLen1,DeltaRefr1,AttenLen2,DeltaRefr2,IgnThick,Thick,PhaseMult,Nzones

Variable xStep = DimDelta($CWaveName, 0), zStep = DimDelta($CWaveName, 1)
Variable xStart = DimOffset($CWaveName, 0), zStart = DimOffset($CWaveName, 1)
Variable nx = DimSize($CWaveName, 0), nz = DimSize($CWaveName, 1)
Variable xEnd = xStart + xStep*(nx - 1), zEnd = zStart + zStep*(nz - 1)

Variable f = FocDist // in m
Variable Lamb_m = 1.239854e-06/SrwBliZonePlatePhotEn // PhotEn in eV
Variable LambK = Lamb_m*Nzones*PhaseMult
Variable rMax = sqrt((f + 0.25*LambK)*LambK)
Variable drMin = sqrt(0.25*f*Lamb_m*PhaseMult/Nzones)

if(SrwAllowPrintingExtraInfo==1)
	print "Outmost Zone Radius:", rMax*1.e+03, "mm"
	print "Outmost Zone Width:", drMin*1.e+06, "microns"
endif

if((-rMax < xStart) %| (rMax > xEnd) %| (-rMax < zStart) %| (rMax > zEnd))
	SrwUtiPrintWarn(SrwPWarnBliThinGenWaveLim)
endif

if((drMin < 1.5*xStep) %| (drMin < 1.5*zStep))
	SrwUtiPrintWarn(SrwPWarnBliThinGenWaveSmp)
endif

end

//+++++++++++++++++++++++++++++++++++++++
// Amplitude Transmission function for simple Zone Plate.
// Returns the Amplitude Transmission for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function AmpTransm_ZonePlateSmpl(xx, zz)
Variable xx, zz

Variable/G SrwBliZonePlateThick, SrwBliZonePlateAttenLen1 //, SrwBliZonePlateAttenLen2
Variable/G SrwBliZonePlateXc, SrwBliZonePlateZc

Variable xr = xx - SrwBliZonePlateXc, zr = zz - SrwBliZonePlateZc
Variable re2 = xr*xr + zr*zr
Variable ZoneParity = SrwOptZonePlateParity(re2)
if(ZoneParity == 0)
	return 0 // Assume zero transmission in outer parts
endif

if(ZoneParity == 2)
	//ExpArg = SrwBliZonePlateThick/SrwBliZonePlateAttenLen2
	return 1.
endif
Variable ExpArg = SrwBliZonePlateThick/SrwBliZonePlateAttenLen1
return exp(-0.5*ExpArg)
end

//+++++++++++++++++++++++++++++++++++++++
// Optical Path function for simple Zone Plate.
// Returns the Optical Path for any point with transverse coordinates (x,z) 
// in a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinSetup.
//+++++++++++++++++++++++++++++++++++++++
function OptPath_ZonePlateSmpl(xx, zz)
Variable xx, zz

Variable/G SrwBliZonePlateThick, SrwBliZonePlateDeltaRefr1 //, SrwBliZonePlateDeltaRefr2
Variable/G SrwBliZonePlateXc, SrwBliZonePlateZc //, SrwBliZonePlatePhotEn, 

Variable xr = xx - SrwBliZonePlateXc, zr = zz - SrwBliZonePlateZc
Variable re2 = xr*xr + zr*zr
//Variable WaveNumb = (5.067681604e+06)*SrwBliZonePlatePhotEn // in 1/m, assuming PhotEn in eV

Variable ZoneParity = SrwOptZonePlateParity(re2)
if((ZoneParity == 0) %| (ZoneParity == 2))
	return 0
endif
//Variable Buf = -WaveNumb*SrwBliZonePlateThick*0.001
Variable Buf = -SrwBliZonePlateThick*0.001

//if(ZoneParity == 1)
return Buf*SrwBliZonePlateDeltaRefr1
//else
//	//return Buf*SrwBliZonePlateDeltaRefr2
//endif
end

//+++++++++++++++++++++++++++++++++++++++
// Obsolete (Optical Path should be used instead)
// Phase Shift function for simple Zone Plate.
// Returns the Phase Shift for any point with transverse coordinates (x,z) 
// on a plane perpendicular to optical axis.
// To simulate an arbitrary "Thin" optical component, 
// user needs to prepare such a function for it
// and submit the name of this function to the macro SrwOptThinGenSetup.
//+++++++++++++++++++++++++++++++++++++++
function PhaseShift_ZonePlateSmpl(xx, zz)
Variable xx, zz

Variable/G SrwBliZonePlateDeltaRefr1, SrwBliZonePlateDeltaRefr2, SrwBliZonePlateThick
Variable/G SrwBliZonePlatePhotEn, SrwBliZonePlateXc, SrwBliZonePlateZc

Variable xr = xx - SrwBliZonePlateXc, zr = zz - SrwBliZonePlateZc
Variable re2 = xr*xr + zr*zr
Variable WaveNumb = (5.067681604e+06)*SrwBliZonePlatePhotEn // in 1/m, assuming PhotEn in eV

Variable ZoneParity = SrwOptZonePlateParity(re2)
if(ZoneParity == 0)
	return 0
endif
Variable Buf = -WaveNumb*SrwBliZonePlateThick*0.001
if(ZoneParity == 1)
	return Buf*SrwBliZonePlateDeltaRefr1
else
	return Buf*SrwBliZonePlateDeltaRefr2
endif
end

//+++++++++++++++++++++++++++++++++++++++
function SrwOptZonePlateParity(re2) // r in m
Variable re2

//Variable/G SrwBliZonePlateFocDist, SrwBliZonePlatePhotEn
Variable/G SrwBliZonePlateNzones //, SrwBliZonePlatePhaseMult
Variable/G SrwBliZonePlateRn

//dddddddddddddddddddddddddddd
//Variable f = SrwBliZonePlateFocDist // in m
//Variable Lamb_m = 1.239854e-06/SrwBliZonePlatePhotEn // PhotEn in eV

//Variable LambK = Lamb_m*SrwBliZonePlateNzones*SrwBliZonePlatePhaseMult
//Variable re2Max = (f + 0.25*LambK)*LambK

Variable re2Max = SrwBliZonePlateRn*SrwBliZonePlateRn*(1e-6)

if(re2 > re2Max)
	return 0
endif
Variable two_k_real = trunc(2*re2*SrwBliZonePlateNzones/re2Max) + 1
if(abs(two_k_real - trunc(two_k_real*0.5)*2.) < 0.001) // k is even
	return 2
else
	return 1
endif

//Variable Oned2f = 0.5/f
//Variable t1 = re2*Oned2f
//Variable t2 = t1*t1*Oned2f
//Variable t3 = t2*t1*Oned2f*0.5

//Variable a = 2./(Lamb_m*SrwBliZonePlatePhaseMult)
//Variable b = t1 - t2 + t3
//Variable k = trunc(a*b) + 1

//if(abs(k - trunc(k*0.5)*2.) < 0.001) // k is even
//	return 2
//else
//	return 1
//endif
end

//+++++++++++++++++++++++++++++++++++++++
// Optical Path Difference function for 1D Bimorph Lens / Mirror
// Requires:
// make/O/N=(2,numParts+1) wThinMirBimorphPartCirCen 
// [0]: transverse start positions of parts; [1]: focal lengths of parts; [2][3]: centers of circles
//+++++++++++++++++++++++++++++++++++++++
function OptPathDif_MirBimorph(xx, zz)
variable xx, zz

wave wThinMirBimorphPartCirCen
NVAR SrwBliThinMirBimFocPlane

variable numParts = dimsize(wThinMirBimorphPartCirCen, 1) - 1

variable xFoc = xx, xPlane = zz
if(SrwBliThinMirBimFocPlane == 2) //foc. in vert. plane
	 xFoc = zz; xPlane = xx
endif

variable zn = wThinMirBimorphPartCirCen[0][0]
variable R = 2*wThinMirBimorphPartCirCen[1][0]
variable zc = wThinMirBimorphPartCirCen[2][0]
variable yc = wThinMirBimorphPartCirCen[3][0], iPart

if((xFoc < zn) %| (xFoc > wThinMirBimorphPartCirCen[0][numParts]))
	return 0
endif

for(iPart = 1; iPart < numParts; iPart += 1)
	zn = wThinMirBimorphPartCirCen[0][iPart]
	if(xFoc <= zn)
		break
	endif
	R = 2*wThinMirBimorphPartCirCen[1][iPart]
	zc = wThinMirBimorphPartCirCen[2][iPart]
	yc = wThinMirBimorphPartCirCen[3][iPart]
endfor
variable dxFoc = xFoc - zc
return -2*(yc - sqrt(R*R - dxFoc*dxFoc))
end

//+++++++++++++++++++++++++++++++++++++++
function AmpTransm_MirBimorph(xx, zz)
variable xx, zz
return 1
end

//+++++++++++++++++++++++++++++++++++++++
//
// 1D Bimorph Lens / Mirror
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirBimorphInit(ElemName,FocPlane,xr,zr,xc,zc,IndCenPart,f0,x1,f1)
string ElemName=srwUtiGetValS("SrwBliThinMirBim", "BiMir", "")
variable FocPlane=srwUtiGetValN("SrwBliThinMirBimFocPlane", 1, "")
variable xr=srwUtiGetValN("SrwBliThinGenWaveRangeX", 1, "")
variable zr=srwUtiGetValN("SrwBliThinGenWaveRangeZ", 1, "")
variable xc=srwUtiGetValN("SrwBliThinGenWavePosX", 0, "")
variable zc=srwUtiGetValN("SrwBliThinGenWavePosZ", 0, "")
variable IndCenPart=srwUtiGetValN("SrwBliThinMirBimIndCenPart", 1, "")
variable f0=srwUtiGetValN("f0", 1, "SrwOptThinMirBimorphInit")
variable x1=srwUtiGetValN("x1", 0, "SrwOptThinMirBimorphInit")
variable f1=srwUtiGetValN("f1", 1, "SrwOptThinMirBimorphInit")
prompt ElemName,SrwPBli
prompt FocPlane,SrwPBliXrayLensFocPlane,popup "Horizontal;Vertical"
prompt xr,"Horizontal Aperture [mm]"
prompt zr,"Vertical Aperture [mm]"
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
prompt IndCenPart,"Indice of Central Mirror Part"
prompt f0,"Focal Length of Part #0 [m]"
prompt x1,"Initial Position of Part #1 [mm]"
prompt f1,"Focal Length of Part #1 [m]"
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwBliThinGen=ElemName
SrwBliLast=ElemName
SrwBliXrayLens=ElemName

srwUtiSetValN("SrwBliThinMirBimFocPlane", FocPlane, "")
srwUtiSetValN("SrwBliThinGenWaveRangeX", xr, "")
srwUtiSetValN("SrwBliThinGenWaveRangeZ", zr, "")
srwUtiSetValN("SrwBliThinGenWavePosX", xc, "")
srwUtiSetValN("SrwBliThinGenWavePosZ", zc, "")
srwUtiSetValN("SrwBliThinMirBimIndCenPart", IndCenPart, "")
srwUtiSetValN("f0", f0, "SrwOptThinMirBimorphInit")
srwUtiSetValN("x1", x1, "SrwOptThinMirBimorphInit")
srwUtiSetValN("f1", f1, "SrwOptThinMirBimorphInit")

variable nConst = 50, nFoc = 4000 // default number of points in approximating mesh

variable x0, xRange, nx, nz
if(FocPlane == 1)
	x0 = xc - 0.5*xr
	xRange = xr
	nx = nFoc
	nz = nConst
	if((x1 < x0) %| (x1 > (xc + 0.5*xr)))
		abort "Initial position of part is out of aperture"
	endif
else
	x0 = zc - 0.5*zr
	xRange = zr
	nx = nConst
	nz = nFoc
	if((x1 < x0) %| (x1 > (zc + 0.5*zr)))
		abort "Initial position of part is out of aperture"
	endif
endif

variable TransmOuter = 1 // outer transmission is zero
SrwOptThinTempl(ElemName,TransmOuter,xc,zc,xr,zr,nx,nz)

ElemName=ElemName+SrwBeamlineType

$ElemName[7]="0" //"1" // Setup was finished
if(FocPlane==1) // Hor.
	$ElemName[8]=num2str(f0); $ElemName[9]=num2str(1e+023)
endif
if(FocPlane==2) // Ver.
	$ElemName[8]=num2str(1e+023); $ElemName[9]=num2str(f0)
endif

make/O/D wThinMirBimorphPartStartPos = {x0,x1,x0+xRange}
wThinMirBimorphPartStartPos *= 0.001
make/O/D wThinMirBimorphPartFocLen = {f0,f1}
end

//+++++++++++++++++++++++++++++++++++++++
//
// 1D Bimorph Lens / Mirror: add parts and setup
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirBimorphSetup(ElemName,FinishSetup,x1,f1,x2,f2,x3,f3,x4,f4)
string ElemName=srwUtiGetValS("SrwBliThinMirBim", "BiMir", "") + SrwBeamlineType
variable NumPartsAdd=srwUtiGetValN("NumPartsAdd", NumPartsAdd, "SrwOptThinMirBimorphSetup")
variable FinishSetup=srwUtiGetValN("FinishSetup", NumPartsAdd, "SrwOptThinMirBimorphSetup")
variable x1=srwUtiGetValN("x1", x1, "SrwOptThinMirBimorphSetup")
variable f1=srwUtiGetValN("f1", f1, "SrwOptThinMirBimorphSetup")
variable x2=srwUtiGetValN("x2", x2, "SrwOptThinMirBimorphSetup")
variable f2=srwUtiGetValN("f2", f2, "SrwOptThinMirBimorphSetup")
variable x3=srwUtiGetValN("x3", x3, "SrwOptThinMirBimorphSetup")
variable f3=srwUtiGetValN("f3", f3, "SrwOptThinMirBimorphSetup")
variable x4=srwUtiGetValN("x4", x4, "SrwOptThinMirBimorphSetup")
variable f4=srwUtiGetValN("f4", f4, "SrwOptThinMirBimorphSetup")
prompt ElemName,SrwPBli, popup Wavelist("*"+SrwBeamlineType,";","")
prompt FinishSetup,"Finish Setup?",popup "No;Yes"
prompt x1,"Initial Position of Next Part [mm]"
prompt f1,"Focal Length of Next Part [m]"
prompt x2,"Initial Position of Next Part [mm]"
prompt f2,"Focal Length of Next Part [m]"
prompt x3,"Initial Position of Next Part [mm]"
prompt f3,"Focal Length of Next Part [m]"
prompt x4,"Initial Position of Next Part [mm]"
prompt f4,"Focal Length of Next Part [m]"
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwBliThinGen=ElemName[0,strlen(ElemName)-strlen(SrwBeamlineType)-1]
SrwBliLast=ElemName
SrwBliXrayLens=ElemName

srwUtiSetValN("NumPartsAdd", NumPartsAdd, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("FinishSetup", FinishSetup, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("x1", x1, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("f1", f1, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("x2", x2, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("f2", f2, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("x3", x3, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("f3", f3, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("x4", x4, "SrwOptThinMirBimorphSetup")
srwUtiSetValN("f4", f4, "SrwOptThinMirBimorphSetup")

if((!exists("wThinMirBimorphPartStartPos")) %| (!exists("wThinMirBimorphPartFocLen")) %| (SrwBliThinMirBimIndCenPart < 0))
	abort "Bimorph Mirror setup was not initiated. \"SrwOptThinMirBimorphInit\" macro must be called first."
endif

variable numStartPos = dimsize(wThinMirBimorphPartStartPos, 0) - 1
variable numFocLen = dimsize(wThinMirBimorphPartFocLen, 0)

if(numStartPos != numFocLen)
	abort "Inconsistent parts' Focal Lengths and Start Positions waves."
endif

make/O/D wAuxBiMirStartPosFocLen = {{x1*0.001,f1},{x2*0.001,f2},{x3*0.001,f3},{x4*0.001,f4}}
variable iNewPart=0, iPart, xi, fi

do
	fi = wAuxBiMirStartPosFocLen[1][iNewPart]
	if(fi != 0)
		xi = wAuxBiMirStartPosFocLen[0][iNewPart]
		if((xi >= wThinMirBimorphPartStartPos[0]) %& (xi < wThinMirBimorphPartStartPos[numFocLen]))
			iPart = 1
			do
				if((xi > wThinMirBimorphPartStartPos[iPart - 1]) %& (xi < wThinMirBimorphPartStartPos[iPart]))
					InsertPoints iPart, 1, wThinMirBimorphPartStartPos, wThinMirBimorphPartFocLen
					wThinMirBimorphPartStartPos[iPart] = xi
					wThinMirBimorphPartFocLen[iPart] = fi
					numFocLen += 1
					numStartPos += 1
					break
				endif
				iPart += 1
			while(iPart <= numFocLen)
		endif
	endif
	iNewPart += 1
while(iNewPart < 4)
killwaves/Z wAuxBiMirStartPosFocLen

if((FinishSetup != 2) %| (numFocLen <= 0))
	return
endif
//Finishing setup

if(SrwBliThinMirBimIndCenPart < 0)
	SrwBliThinMirBimIndCenPart = 0
endif
if(SrwBliThinMirBimIndCenPart >= numFocLen)
	SrwBliThinMirBimIndCenPart = numFocLen - 1
endif

//make/O/N=(2,numFocLen) wThinMirBimorphPartCirCen
make/D/O/N=(4,numFocLen+1) wThinMirBimorphPartCirCen

variable fn = wThinMirBimorphPartFocLen[SrwBliThinMirBimIndCenPart]
variable zc0 = 0, yc0 = 2*fn

wThinMirBimorphPartCirCen[0][SrwBliThinMirBimIndCenPart] =  wThinMirBimorphPartStartPos[SrwBliThinMirBimIndCenPart]
wThinMirBimorphPartCirCen[1][SrwBliThinMirBimIndCenPart] = fn
wThinMirBimorphPartCirCen[2][SrwBliThinMirBimIndCenPart] = zc0
wThinMirBimorphPartCirCen[3][SrwBliThinMirBimIndCenPart] = yc0

variable zcnPrev = zc0, ycnPrev = yc0, zn, yn, zcn, ycn, Rn, RnPrev = yc0, zn_mi_zcn, zn_mi_zcnPrev
iPart = SrwBliThinMirBimIndCenPart - 1
if(iPart >= 0)
	do
		fn = wThinMirBimorphPartFocLen[iPart]
		Rn = 2*fn
		zn = wThinMirBimorphPartStartPos[iPart + 1]
		
		zcn = (Rn/RnPrev)*zcnPrev + (1 - (Rn/RnPrev))*zn
		zn_mi_zcn = zn - zcn
		zn_mi_zcnPrev = zn - zcnPrev
		ycn = ycnPrev + sqrt(Rn*Rn - zn_mi_zcn*zn_mi_zcn) - sqrt(RnPrev*RnPrev - zn_mi_zcnPrev*zn_mi_zcnPrev)

		wThinMirBimorphPartCirCen[2][iPart] = zcn
		wThinMirBimorphPartCirCen[3][iPart] = ycn
		
		wThinMirBimorphPartCirCen[1][iPart] = fn
		wThinMirBimorphPartCirCen[0][iPart] =  wThinMirBimorphPartStartPos[iPart]

		RnPrev = Rn
		zcnPrev = zcn
		ycnPrev = ycn
		iPart -= 1
	while(iPart >= 0)
endif
iPart = SrwBliThinMirBimIndCenPart + 1
if(iPart < numFocLen)
	zcnPrev = zc0; ycnPrev = yc0; RnPrev = yc0
	do
		fn = wThinMirBimorphPartFocLen[iPart]
		Rn = 2*fn
		zn = wThinMirBimorphPartStartPos[iPart]
 
		zcn = (Rn/RnPrev)*zcnPrev + (1 - (Rn/RnPrev))*zn
		zn_mi_zcn = zn - zcn
		zn_mi_zcnPrev = zn - zcnPrev
		
		ycn = ycnPrev + sqrt(Rn*Rn - zn_mi_zcn*zn_mi_zcn) - sqrt(RnPrev*RnPrev - zn_mi_zcnPrev*zn_mi_zcnPrev)
		//print zn, fn, zcn, ycn, 2*(ycn - sqrt(Rn*Rn - zn_mi_zcn*zn_mi_zcn)), 2*(ycnPrev - sqrt(RnPrev*RnPrev - zn_mi_zcnPrev*zn_mi_zcnPrev))

		wThinMirBimorphPartCirCen[2][iPart] = zcn
		wThinMirBimorphPartCirCen[3][iPart] = ycn
		wThinMirBimorphPartCirCen[1][iPart] = fn
		wThinMirBimorphPartCirCen[0][iPart] = zn

		RnPrev = Rn
		zcnPrev = zcn
		ycnPrev = ycn
		iPart += 1
	while(iPart < numFocLen)
endif

wThinMirBimorphPartCirCen[0][numFocLen] =  wThinMirBimorphPartStartPos[numFocLen]
killwaves/Z wThinMirBimorphPartStartPos, wThinMirBimorphPartFocLen

string CWaveName = SrwBliThinGen + SrwOptThinGenWaveType
SrwOptThinSetup(CWaveName,"AmpTransm_MirBimorph","OptPathDif_MirBimorph")
killwaves/Z wThinMirBimorphPartCirCen
end