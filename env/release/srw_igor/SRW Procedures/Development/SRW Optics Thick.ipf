
//+++++++++++++++++++++++++++++++++++++++
//
//Initialise constants
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickInit() 
string/G SrwBliThickMirTorType = "ThickMirrorToroid"

end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror, "thick" approximation
//
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThickMirTorDialogs()

NVAR SrwOptThickMirInitPassed, SrwOptThickMirSetupPassed
SrwOptThickMirInitPassed=0
execute "SrwOptThickMirTorInit()"
if(SrwOptThickMirInitPassed == 0) 
	return 0
endif
string ComLineStr = ""
SVAR SrwBliThick
NVAR SrwOptThickMirRefl,SrwOptThickMirRt,SrwOptThickMirRs,SrwOptThickMirSizeT,SrwOptThickMirSizeS,SrwOptThickMirCenPosH,SrwOptThickMirCenPosV,SrwOptThickMirCenPosL
sprintf ComLineStr, "SrwOptThickMirTorInit(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwBliThick,SrwOptThickMirRefl,SrwOptThickMirRt,SrwOptThickMirRs,SrwOptThickMirSizeT,SrwOptThickMirSizeS,SrwOptThickMirCenPosH,SrwOptThickMirCenPosV,SrwOptThickMirCenPosL

SrwOptThickMirSetupPassed = 0
execute "SrwOptThickMirTorSetup()"
if(SrwOptThickMirSetupPassed == 0) 
	print ComLineStr
	return 0
endif

string ComLineStr2 = ""
NVAR SrwOptThickMirNx,SrwOptThickMirNy,SrwOptThickMirNz,SrwOptThickMirRotAng
sprintf ComLineStr2, "SrwOptThickMirTorSetup(\"%s\",%g,%g,%g,%g)", SrwBliThick,SrwOptThickMirNx,SrwOptThickMirNy,SrwOptThickMirNz,SrwOptThickMirRotAng
ComLineStr += (";" + ComLineStr2)
print ComLineStr
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror, "thick" approximation: Init Setup
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirTorInit(name, refl, rt, rs, dt, ds, shape, yc, xc, zc)
string name=srwUtiGetValS("SrwBliThick", "ThickTorMir", "")
variable refl=srwUtiGetValN("SrwOptThickMirRefl", 1, "")
variable rt=srwUtiGetValN("SrwOptThickMirRt", 1, "")
variable rs=srwUtiGetValN("SrwOptThickMirRs", 1, "")
variable dt=srwUtiGetValN("SrwOptThickMirSizeT", 0.05, "")*1000
variable ds=srwUtiGetValN("SrwOptThickMirSizeS", 0.05, "")*1000
variable shape=srwUtiGetValN("SrwOptThickMirShape", 1, "")
variable yc=srwUtiGetValN("SrwOptThickMirCenPosL", 0, "")
variable xc=srwUtiGetValN("SrwOptThickMirCenPosH", 0, "")*1000
variable zc=srwUtiGetValN("SrwOptThickMirCenPosV", 0, "")*1000
prompt name,"Name of the Mirror structure"
prompt refl,"Intensity Reflectivity Coefficient"
prompt rt,"Tangential (or Major) Radius [m]"
prompt rs,"Sagittal (or Minor) Radius [m]"
prompt dt,"Mirror Size in Tangential Plane [mm]"
prompt ds,"Mirror Size in Sagittal Plane [mm]"
prompt shape,"Mirror Aperture Shape in Local Frame",popup "rectangular;elliptical"
prompt yc,"Longitudinal Position of the Mirror Center [mm]"
prompt xc,"Horizontal Position of the Mirror Center [mm]"
prompt zc,"Vertical Position of the Mirror Center [mm]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThick", name, "")
srwUtiSetValN("SrwOptThickMirRefl", refl, "")
srwUtiSetValN("SrwOptThickMirRt", rt, "")
srwUtiSetValN("SrwOptThickMirRs", rs, "")
srwUtiSetValN("SrwOptThickMirSizeT", dt*0.001, "")
srwUtiSetValN("SrwOptThickMirSizeS", ds*0.001, "")
srwUtiSetValN("SrwOptThickMirSizeS", ds*0.001, "")
srwUtiSetValN("SrwOptThickMirShape", shape, "")
srwUtiSetValN("SrwOptThickMirCenPosL", yc, "")
srwUtiSetValN("SrwOptThickMirCenPosH", xc*0.001, "")
srwUtiSetValN("SrwOptThickMirCenPosV", zc*0.001, "")

variable/G SrwOptThickMirInitPassed = 1
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror, "thick" approximation: Finish Setup
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirTorSetup(name, nx, ny, nz, rotang)
string name=srwUtiGetValS("SrwBliThick", "ThickTorMir", "")
variable nx=srwUtiGetValN("SrwOptThickMirNx", 0, "")
variable ny=srwUtiGetValN("SrwOptThickMirNy", -1, "")
variable nz=srwUtiGetValN("SrwOptThickMirNz", 1, "")
variable rotang=srwUtiGetValN("SrwOptThickMirRotAng", 0, "")
prompt name,"Name of the Mirror structure"
prompt nx,"Horizontal Coordinate of the Central Normal Vector [mm]"
prompt ny,"Longitudinal Coordinate of the Central Normal Vector [mm]"
prompt nz,"Vertical Coordinate of the Central Normal Vector [mm]"
prompt rotang,"Rotation Angle about Cetral Normal Vector [rad]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThick", name, "")
srwUtiSetValN("SrwOptThickMirNx", nx, "")
srwUtiSetValN("SrwOptThickMirNy", ny, "")
srwUtiSetValN("SrwOptThickMirNz", nz, "")
srwUtiSetValN("SrwOptThickMirRotAng", rotang, "")

name+=SrwBeamlineType
srwUtiGetValS("SrwBliThickMirTorType", "ThickMirrorToroid", "")

make/T/O/N=15 $name
$name[0]=SrwBliThickMirTorType
$name[1]=num2str(SrwOptThickMirRefl)
$name[2]=num2str(SrwOptThickMirRt)
$name[3]=num2str(SrwOptThickMirRs)
$name[4]=num2str(SrwOptThickMirSizeT)
$name[5]=num2str(SrwOptThickMirSizeS)
$name[6]=num2str(SrwOptThickMirShape)
$name[7]=num2str(SrwOptThickMirCenPosL)
$name[8]=num2str(SrwOptThickMirCenPosH)
$name[9]=num2str(SrwOptThickMirCenPosV)
$name[10]=num2str(SrwOptThickMirNx)
$name[11]=num2str(SrwOptThickMirNy)
$name[12]=num2str(SrwOptThickMirNz)
$name[13]=num2str(SrwOptThickMirRotAng)

variable/G SrwOptThickMirSetupPassed = 1
end
