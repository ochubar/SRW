
//+++++++++++++++++++++++++++++++++++++++
//
//Optical Element: Thick Generic Mirror
//
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThickMirGenDialogs()

NVAR SrwOptThickMirGenInitPassed,SrwOptThickMirGenSetupPassed

SrwOptThickMirGenInitPassed=0
execute "SrwOptThickMirGenInit()"
if(SrwOptThickMirGenInitPassed == 0) 
	return 0
endif
SVAR SrwBliThickMirGen
NVAR SrwBliThickMirGenApShape,SrwBliThickMirGenSizeH,SrwBliThickMirGenSizeV,SrwBliThickMirGenNpH,SrwBliThickMirGenNpV,SrwBliThickMirAmpRefPerp,SrwBliThickMirPhShiftPerp,SrwBliThickMirAmpRefPar,SrwBliThickMirPhShiftPar
string ComLineStr = ""
sprintf ComLineStr, "SrwOptThickMirGenInit(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g,%g)",SrwBliThickMirGen,SrwBliThickMirGenApShape,SrwBliThickMirGenSizeH*1000,SrwBliThickMirGenSizeV*1000,SrwBliThickMirGenNpH,SrwBliThickMirGenNpV,SrwBliThickMirAmpRefPerp,SrwBliThickMirPhShiftPerp,SrwBliThickMirAmpRefPar,SrwBliThickMirPhShiftPar

SrwOptThickMirGenSetupPassed = 0
execute "SrwOptThickMirGenSetup()"
if(SrwOptThickMirGenSetupPassed == 0) 
	print ComLineStr
	return 0
endif
SVAR SrwBeamlineType, SrwBliThickMirGenSurfFunc
NVAR SrwBliThickMirGenPosH,SrwBliThickMirGenPosV,SrwBliThickMirGenAxRot1,SrwBliThickMirGenAngRot1,SrwBliThickMirGenAxRot2,SrwBliThickMirGenAngRot2,SrwBliThickMirGenAxRot3,SrwBliThickMirGenAngRot3
string ComLineStr2 = ""
sprintf ComLineStr2, "SrwOptThickMirGenSetup(\"%s\",\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)",SrwBliThickMirGen+SrwBeamlineType,SrwBliThickMirGenSurfFunc,SrwBliThickMirGenPosH*1000,SrwBliThickMirGenPosV*1000,SrwBliThickMirGenAxRot1,SrwBliThickMirGenAngRot1,SrwBliThickMirGenAxRot2,SrwBliThickMirGenAngRot2,SrwBliThickMirGenAxRot3,SrwBliThickMirGenAngRot3
ComLineStr += (";" + ComLineStr2)
print ComLineStr
end

//+++++++++++++++++++++++++++++++++++++++
//
//Optical Element: Thick Generic Mirror: Init
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirGenInit(name,apertShape,xr,yr,npx,npy,ampRefPerp,phShiftPerp,ampRefPar,phShiftPar)
string name=srwUtiGetValS("SrwBliThickMirGen", "ThickMir", "")
variable apertShape=srwUtiGetValN("SrwBliThickMirGenApShape", 1, "")
variable xr=srwUtiGetValN("SrwBliThickMirGenSizeH", 10, "")*1000
variable yr=srwUtiGetValN("SrwBliThickMirGenSizeV", 10, "")*1000
variable npx=srwUtiGetValN("SrwBliThickMirGenNpH", 400, "")
variable npy=srwUtiGetValN("SrwBliThickMirGenNpV", 400, "")
variable ampRefPerp=srwUtiGetValN("SrwBliThickMirAmpRefPerp", 1, "")
variable phShiftPerp=srwUtiGetValN("SrwBliThickMirPhShiftPerp", 0, "")
variable ampRefPar=srwUtiGetValN("SrwBliThickMirAmpRefPar", 1, "")
variable phShiftPar=srwUtiGetValN("SrwBliThickMirPhShiftPar", Pi, "")
prompt name,"Name for the Mirror structure"
prompt apertShape,"Mirror Aperture Shape",popup "rectangular;elliptical"
prompt xr,"\"Horizontal\" Mirror Size in Local Frame [mm]"
prompt yr,"\"Vertical\" Mirror Size in Local Frame [mm]"
prompt npx,"Num. of Points to repr. Surface in \"Horiz.\" dir."
prompt npy,"Num. of Points to repr. Surface in \"Vert.\" dir."
prompt ampRefPerp,"Ampl. Reflection for Polar. Perp. to Inc. Plane"
prompt phShiftPerp,"Phase Shift for Polar. Perp. to Inc. Plane"
prompt ampRefPar,"Ampl. Reflection for Polar. Par. to Inc. Plane"
prompt phShiftPar,"Phase Shift for Polar. Parallel to Inc. Plane"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThickMirGen", name, "")
srwUtiSetValN("SrwBliThickMirGenApShape", apertShape, "")
srwUtiSetValN("SrwBliThickMirGenSizeH", xr*0.001, "")
srwUtiSetValN("SrwBliThickMirGenSizeV", yr*0.001, "")
srwUtiSetValN("SrwBliThickMirGenNpH", npx, "")
srwUtiSetValN("SrwBliThickMirGenNpV", npy, "")
srwUtiSetValN("SrwBliThickMirAmpRefPerp", ampRefPerp, "")
srwUtiSetValN("SrwBliThickMirPhShiftPerp", phShiftPerp, "")
srwUtiSetValN("SrwBliThickMirAmpRefPar", ampRefPar, "")
srwUtiSetValN("SrwBliThickMirPhShiftPar", phShiftPar, "")

name += SrwBeamlineType
make/T/O/N=20 $name

variable/G SrwOptThickMirGenInitPassed = 1
end

//+++++++++++++++++++++++++++++++++++++++
//
//Optical Element: Thick Generic Mirror: Orient
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirGenSetup(name,nmSurfFunc,xc,yc,axRot1,angRot1,axRot2,angRot2,axRot3,angRot3)
string name=srwUtiGetValS("SrwBliThickMirGen", "ThickMir", "")+SrwBeamlineType
string nmSurfFunc=srwUtiGetValS("SrwBliThickMirGenSurfFunc", "MirSurfFunc", "")
variable xc=srwUtiGetValN("SrwBliThickMirGenPosH", 0, "")*1000
variable yc=srwUtiGetValN("SrwBliThickMirGenPosV", 0, "")*1000
variable axRot1=srwUtiGetValN("SrwBliThickMirGenAxRot1", 1, "")
variable angRot1=srwUtiGetValN("SrwBliThickMirGenAngRot1", 0, "")
variable axRot2=srwUtiGetValN("SrwBliThickMirGenAxRot2", 1, "")
variable angRot2=srwUtiGetValN("SrwBliThickMirGenAngRot2", 0, "")
variable axRot3=srwUtiGetValN("SrwBliThickMirGenAxRot3", 1, "")
variable angRot3=srwUtiGetValN("SrwBliThickMirGenAngRot3", 0, "")
prompt name,"Name of the Mirror structure",popup WaveList("*"+SrwBeamlineType ,";", "")
prompt nmSurfFunc,"Name of the Mirror Surface function",popup FunctionList("*",";","KIND:2,NPARAMS:2,VALTYPE:1")
prompt xc,"\"Horizontal\" Center Position [mm]"
prompt yc,"\"Vertical\" Center Position [mm]"
prompt axRot1,"Mirror Rotation #1: Axis (resp. to inc. beam)",popup "_none_;\"horizontal\";\"vertical\";\"longitudinal\""
prompt angRot1,"Mirrot Rotation #1: Angle [rad]"
prompt axRot2,"Mirror Rotation #2: Axis (resp. to inc. beam)",popup "_none_;\"horizontal\" (resp. to incident beam);\"vertical\";\"longitudinal\""
prompt angRot2,"Mirror Rotation #2: Angle [rad]"
prompt axRot3,"Mirror Rotation #3: Axis (resp. to inc. beam)",popup "_none_;\"horizontal\" (resp. to incident beam);\"vertical\";\"longitudinal\""
prompt angRot3,"Mirror Rotation #3: Angle [rad]"
Silent 1						|	 ...
PauseUpdate

string nameCore = name[0,strlen(name)-strlen(SrwBeamlineType)-1]

srwUtiSetValS("SrwBliThickMirGen", nameCore, "")
srwUtiSetValS("SrwBliThickMirGenSurfFunc", nmSurfFunc, "")
xc = xc*0.001; yc = yc*0.001
srwUtiSetValN("SrwBliThickMirGenPosH", xc, "")
srwUtiSetValN("SrwBliThickMirGenPosV", yc, "")
srwUtiSetValN("SrwBliThickMirGenAxRot1", axRot1, "")
srwUtiSetValN("SrwBliThickMirGenAngRot1", angRot1, "")
srwUtiSetValN("SrwBliThickMirGenAxRot2", axRot2, "")
srwUtiSetValN("SrwBliThickMirGenAngRot2", angRot2, "")
srwUtiSetValN("SrwBliThickMirGenAxRot3", axRot3, "")
srwUtiSetValN("SrwBliThickMirGenAngRot3", angRot3, "")

variable apertShape=srwUtiGetValN("SrwBliThickMirGenApShape", 1, "")
variable xr=srwUtiGetValN("SrwBliThickMirGenSizeH", 10, "")
variable yr=srwUtiGetValN("SrwBliThickMirGenSizeV", 10, "")
variable npx = srwUtiGetValN("SrwBliThickMirGenNpH", 400, "")
variable npy = srwUtiGetValN("SrwBliThickMirGenNpV", 400, "")
variable ampRefPerp=srwUtiGetValN("SrwBliThickMirAmpRefPerp", 1, "")
variable phShiftPerp=srwUtiGetValN("SrwBliThickMirPhShiftPerp", 0, "")
variable ampRefPar=srwUtiGetValN("SrwBliThickMirAmpRefPar", 1, "")
variable phShiftPar=srwUtiGetValN("SrwBliThickMirPhShiftPar", Pi, "")

string nameExt = srwUtiGetValS("SrwBliThickMirGenSurfType", SrwSeparator+"bms", "")
string SWaveName=nameCore+nameExt
killwaves/Z $SWaveName
make/D/O/N=(npx,npy) $SWaveName
SetScale/I x (xc-0.5*xr),(xc+0.5*xr),"m",$SWaveName
SetScale/I y (yc-0.5*yr),(yc+0.5*yr),"m",$SWaveName

String ComLineStr
sprintf ComLineStr, "'%s'=%s(x,y)", SWaveName, nmSurfFunc
execute ComLineStr

$name[0] = srwUtiGetValS("SrwBliThickMirGenID", "ThickMirrorGen", "")
$name[1] = SWaveName

$name[4]=num2str(xc)
$name[5]=num2str(yc)
$name[6]=num2str(apertShape)
$name[7]="0" // Setup was finished or not
//8 - foc. dist. x
//9 - foc. dist. z
$name[10]=num2str(ampRefPerp)
$name[11]=num2str(phShiftPerp)
$name[12]=num2str(ampRefPar)
$name[13]=num2str(phShiftPar)
$name[14]=num2str(axRot1 - 1) //"0" means no rotation, "1" means vs "horizontal" axis, ...
$name[15]=num2str(angRot1)
$name[16]=num2str(axRot2 - 1)
$name[17]=num2str(angRot2)
$name[18]=num2str(axRot3 - 1)
$name[19]=num2str(angRot3)

//to determine focal distances:
if(exists(name) == 1)
	srOptThickMirGenSetup($name) // Finish Setup in C (analize focal distances, etc.)
endif

variable/G SrwOptThickMirGenSetupPassed = 1
end
