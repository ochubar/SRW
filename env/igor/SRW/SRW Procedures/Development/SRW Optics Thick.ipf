
//+++++++++++++++++++++++++++++++++++++++
//
//Initialise constants
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwOptThickInit() 
//string/G SrwBliThickMirTorType = "ThickMirrorToroid"
//end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror, "thick" approximation: Init 
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwOptThickMirTorInit(name, rt, rs, shape, dt, ds, xc, zc, nmrs, nmrp)
//proc SrwOptThickMirTorInit(name, rt, rs, shape, dt, ds, nmr)
proc SrwOptThickMirInitTor(name, rt, rs)
string name=srwUtiGetValS("SrwBliThick", "ThickMir", "")
//variable refl=srwUtiGetValN("SrwOptThickMirRefl", 1, "")
variable rt=srwUtiGetValN("SrwOptThickMirRt", 1, "")
variable rs=srwUtiGetValN("SrwOptThickMirRs", 1, "")
//variable dt=srwUtiGetValN("SrwOptThickMirSizeT", 0.05, "")*1000
//variable ds=srwUtiGetValN("SrwOptThickMirSizeS", 0.05, "")*1000
//variable shape=srwUtiGetValN("SrwOptThickMirShape", 1, "")
//variable yc=srwUtiGetValN("SrwOptThickMirCenPosL", 0, "")
//variable xc=srwUtiGetValN("SrwOptThickMirCenPosH", 0, "")*1000
//variable zc=srwUtiGetValN("SrwOptThickMirCenPosV", 0, "")*1000
//variable meth=srwUtiGetValN("SrwOptThickMirPropMeth", 1, "")
//variable npart=srwUtiGetValN("SrwOptThickMirPropPart", 1, "")
//string nmr=srwUtiGetValS("SrwOptRefl", "wRefl", "")
prompt name,"Name for the Mirror structure"
//prompt refl,"Intensity Reflectivity Coefficient"
prompt rt,"Tangential (Major) Radius [m]"
prompt rs,"Sagittal (Minor) Radius [m]"
//prompt dt,"Mirror Size in Tangential Plane [mm]"
//prompt ds,"Mirror Size in Sagittal Plane [mm]"
//prompt shape,"Mirror Aperture Shape in Local Frame",popup "rectangular;elliptical"
//prompt nmr,"Compl. Reflect. of Sigma and Pi Comp.", popup " ;" + Wavelist("*",";","TEXT:0,CMPLX:1,DIMS:3")
//prompt yc,"Longitudinal Position of the Mirror Center [mm]"
//prompt xc,"Horizontal Position of the Mirror Center [mm]"
//prompt zc,"Vertical Position of the Mirror Center [mm]"
//prompt meth,"Simulation Method",popup "Sequence of Thin Elements"
//prompt npart,"Number of Sub-Elements"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThick", name, "")
//srwUtiSetValN("SrwOptThickMirRefl", refl, "")
srwUtiSetValN("SrwOptThickMirRt", rt, "")
srwUtiSetValN("SrwOptThickMirRs", rs, "")

//srwUtiSetValN("SrwOptThickMirShape", shape, "")
//srwUtiSetValN("SrwOptThickMirSizeT", dt*0.001, "")
//srwUtiSetValN("SrwOptThickMirSizeS", ds*0.001, "")
//srwUtiSetValS("SrwOptRefl", nmr, "")

//srwUtiSetValN("SrwOptThickMirCenPosL", yc, "")
//srwUtiSetValN("SrwOptThickMirCenPosH", xc*0.001, "")
//srwUtiSetValN("SrwOptThickMirCenPosV", zc*0.001, "")
//srwUtiSetValN("SrwOptThickMirPropMeth", meth, "")
//srwUtiSetValN("SrwOptThickMirPropPart", npart, "")

//srwUtiSetValS("SrwOptReflPi", nmrp, "")

SrwBliLast=name
name+=SrwBeamlineType
//srwUtiGetValS("SrwBliThickMirTorType", "ThickMirrorToroid", "")

make/T/O/N=30 $name
$name[0]="Mirror" //SrwBliThickMirTorType
$name[1]="Toroid"

$name[2]=num2str(SrwOptThickMirRt)
$name[3]=num2str(SrwOptThickMirRs)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Ellipsoidal mirror, "thick" approximation: Init
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirInitEl(name, p, q, angGraz, rSag)
string name=srwUtiGetValS("SrwBliThick", "ThickMir", "")
variable p=srwUtiGetValN("SrwOptThickMirElP", 1, "")
variable q=srwUtiGetValN("SrwOptThickMirElQ", 1, "")
variable angGraz=srwUtiGetValN("SrwOptThickMirElAngGraz", 1e-02, "")
variable rSag=srwUtiGetValN("SrwOptThickMirElRSag", 1000, "")
prompt name, "Name for the Mirror structure"
prompt p, "Distance from First Focus (\"source\") to Mirror Center [m]"
prompt q, "Distance from Mirror Center to Second Focus (\"image\") [m]"
prompt angGraz, "Grazing Angle at Mirror Center at perfect orientation [rad]"
prompt rSag, "Sagital Radius of Curvature at Mirror Center [rad]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThick", name, "")
srwUtiSetValN("SrwOptThickMirElP", p, "")
srwUtiSetValN("SrwOptThickMirElQ", q, "")
srwUtiSetValN("SrwOptThickMirElAngGraz", angGraz, "")
srwUtiSetValN("SrwOptThickMirElRSag", rSag, "")

SrwBliLast=name
name+=SrwBeamlineType

make/T/O/N=30 $name
$name[0]="Mirror"
$name[1]="Ellipsoid"

$name[2]=num2str(SrwOptThickMirElP)
$name[3]=num2str(SrwOptThickMirElQ)
$name[4]=num2str(angGraz)
$name[5]=num2str(rSag)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Mirror: dimensions, reflectivity, simulation method
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThickMirDimReflMeth(name, shape, dt, ds, nmr, meth, npt, nps)
string name=srwUtiGetValS("SrwBliThick", "ThickMir", "")+SrwBeamlineType
variable shape=srwUtiGetValN("SrwOptThickMirShape", 1, "")
variable dt=srwUtiGetValN("SrwOptThickMirSizeT", 0.05, "")*1000
variable ds=srwUtiGetValN("SrwOptThickMirSizeS", 0.05, "")*1000
string nmr=srwUtiGetValS("SrwOptRefl", "wRefl", "")
variable meth=srwUtiGetValN("SrwOptThickMirPropMeth", 1, "")
variable npt=srwUtiGetValN("SrwOptThinMirNumPointsTan", 0, "")
variable nps=srwUtiGetValN("SrwOptThinMirNumPointsSag", 0, "")
prompt name,"Name of the Mirror structure", popup Wavelist("*"+SrwBeamlineType ,";", "");
prompt dt,"Mirror Size in Tangential Plane [mm]"
prompt ds,"Mirror Size in Sagittal Plane [mm]"
prompt shape,"Mirror Aperture Shape in Local Frame",popup "rectangular;elliptical"
prompt nmr,"Compl. Reflect. of Sigma and Pi Comp.", popup " ;" + Wavelist("*",";","TEXT:0,CMPLX:1,DIMS:3")
prompt meth,"Simulation Method",popup "Thin Element;Local Ray-Tracing"//;Local Ray-Tracing with Diffraction"
prompt npt,"Num. of Points to repr. Surface in Tang. Dir."
prompt nps,"Num. of Points to repr. Surface in Sag. Dir."
Silent 1						|	 ...
PauseUpdate

srwUtiSetValN("SrwOptThickMirShape", shape, "")
srwUtiSetValN("SrwOptThickMirSizeT", dt*0.001, "")
srwUtiSetValN("SrwOptThickMirSizeS", ds*0.001, "")
srwUtiSetValS("SrwOptRefl", nmr, "")
srwUtiSetValN("SrwOptThinMirNumPointsTan", npt, "")
srwUtiSetValN("SrwOptThinMirNumPointsSag", nps, "")

SrwBliLast=name[0,strlen(name)-strlen(SrwBeamlineType)-1]
srwUtiSetValS("SrwBliThick", SrwBliLast, "")

if(dimsize($name,0) < 30)
	redimension/N=30 $name
endif

$name[6]="" //reserved for opt. elem. params
$name[7]="" //reserved for opt. elem. params
$name[8]="" //reserved for Fx
$name[9]="" //reserved for Fy

$name[10]=num2str(SrwOptThickMirSizeT)
$name[11]=num2str(SrwOptThickMirSizeS)
$name[12]=num2str(SrwOptThickMirShape)
$name[13]="" //reserved

$name[14]=SrwOptRefl

$name[26]=num2str(SrwOptThickMirPropMeth)
$name[27]="" //reserved
$name[28]=num2str(npt)
$name[29]=num2str(nps)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Thick mirror: orientation (one function for all mirror shapes: toroid, ellipsoid, paraboloid, spherical,...)
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwOptThickOrient(name, nh, nv, nl, th, tv, meth, npart, endpr)
proc SrwOptThickMirOrient(name, nh, nv, nl, th, tv, xc, zc)//, meth) //, npart)
string name=srwUtiGetValS("SrwBliThick", "ThickMir", "")+SrwBeamlineType
variable nh=srwUtiGetValN("SrwOptThickMirNh", 0, "")
variable nv=srwUtiGetValN("SrwOptThickMirNv", 0.5, "")
variable nl=srwUtiGetValN("SrwOptThickMirNl", -0.5, "")
variable th=srwUtiGetValN("SrwOptThickMirTh", 1, "")
variable tv=srwUtiGetValN("SrwOptThickMirTv", 0, "")
variable xc=srwUtiGetValN("SrwOptThickMirCenPosH", 0, "")*1000
variable zc=srwUtiGetValN("SrwOptThickMirCenPosV", 0, "")*1000
//variable meth=srwUtiGetValN("SrwOptThickMirPropMeth", 1, "")
//variable npart=srwUtiGetValN("SrwOptThickMirPropPart", 1, "")
//variable endpr=srwUtiGetValN("SrwOptThickMirPropEnd", 1, "")
//variable meth=srwUtiGetValN("SrwOptThickMirPropMeth", 1, "")
//variable npart=srwUtiGetValN("SrwOptThickMirPropPart", 1, "")
prompt name,"Name of the Mirror structure", popup Wavelist("*"+SrwBeamlineType ,";", "");
prompt nh,"Horiz. Coordinate of Central Normal"
prompt nv,"Vert. Coordinate of Central Normal"
prompt nl,"Long. Coord. of Central Normal"
prompt th,"Horiz. Coord. of Central Tangent. Vect."
prompt tv,"Vert. Coord. of Central Tangent. Vect."
prompt xc,"Horiz. Position of Mirror Center [mm]"
prompt zc,"Vert. Position of  Mirror Center [mm]"
//prompt meth,"Simulation Method",popup "Thin Element;Local Ray-Tracing;Local Ray-Tracing with Diffraction"
//prompt npart,"Number of Sub-Elements"
//prompt endpr,"After the Propagation:",popup "Stop at This Element End;Propagate to This Element Center"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValN("SrwOptThickMirNh", nh, "")
srwUtiSetValN("SrwOptThickMirNv", nv, "")
srwUtiSetValN("SrwOptThickMirNl", nl, "")
srwUtiSetValN("SrwOptThickMirTh", th, "")
srwUtiSetValN("SrwOptThickMirTv", tv, "")
srwUtiSetValN("SrwOptThickMirCenPosH", xc*0.001, "")
srwUtiSetValN("SrwOptThickMirCenPosV", zc*0.001, "")
//srwUtiSetValN("SrwOptThickMirPropMeth", meth, "")
//srwUtiSetValN("SrwOptThickMirPropPart", npart, "")
//srwUtiSetValN("SrwOptThickMirPropEnd", endpr, "")

SrwBliLast=name[0,strlen(name)-strlen(SrwBeamlineType)-1]
srwUtiSetValS("SrwBliThick", SrwBliLast, "")

if(dimsize($name,0) < 30)
	redimension/N=30 $name
endif

$name[15]="" //reserved

$name[16]=num2str(SrwOptThickMirNh)
$name[17]=num2str(SrwOptThickMirNv)
$name[18]=num2str(SrwOptThickMirNl)
$name[19]=num2str(SrwOptThickMirTh)
$name[20]=num2str(SrwOptThickMirTv)

$name[21]="" //reserved
$name[22]="" //reserved

$name[23]=num2str(SrwOptThickMirCenPosH)
$name[24]=num2str(SrwOptThickMirCenPosV)
$name[25]="" //reserved
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror, "thick" approximation
//
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThickMirTorDialogs_OLD()

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
proc SrwOptThickMirTorInit_OLD(name, refl, rt, rs, dt, ds, shape, yc, xc, zc)
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
proc SrwOptThickMirTorSetup_OLD(name, nx, ny, nz, rotang)
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
