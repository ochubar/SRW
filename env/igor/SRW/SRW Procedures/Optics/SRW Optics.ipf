
//+++++++++++++++++++++++++++++++++++++++
//
//Optical Elements
//
//+++++++++++++++++++++++++++++++++++++++
//
//Create container and add optical elements to it
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptContFull(cntname,cmp1name,cmp2name,cmp3name,cmp4name)
String cntname=SrwBliCont
String cmp1name
String cmp2name
String cmp3name
String cmp4name
prompt cntname,"Name of the Container"
prompt cmp1name,"1st Component to add",popup " ;" + Wavelist("*"+SrwBeamlineType,";","")
prompt cmp2name,"2nd Component to add",popup " ;" + Wavelist("*"+SrwBeamlineType,";","")
prompt cmp3name,"3rd Component to add",popup " ;" + Wavelist("*"+SrwBeamlineType,";","")
prompt cmp4name,"4th Component to add",popup " ;" + Wavelist("*"+SrwBeamlineType,";","")
Silent 1				|	 ...
PauseUpdate
SrwBliCont=cntname
SrwBliLast=cntname

cntname+=SrwBeamlineType

if((cmpstr(cntname,cmp1name)==0) %| (cmpstr(cntname,cmp2name)==0) %| (cmpstr(cntname,cmp3name)==0) %| (cmpstr(cntname,cmp4name)==0))
abort  "Can not add container to itself"
endif

Make/T/O/N=1 $cntname
$cntname[0]=SrwBliContType

if((cmpstr(cmp1name," ") != 0) %& (cmpstr(cmp1name,"") != 0))
	SrwOptContAdd(cmp1name,cntname,1)
endif
if((cmpstr(cmp2name," ") != 0) %& (cmpstr(cmp2name,"") != 0))
	SrwOptContAdd(cmp2name,cntname,1)
endif
if((cmpstr(cmp3name," ") != 0) %& (cmpstr(cmp3name,"") != 0))
	SrwOptContAdd(cmp3name,cntname,1)
endif
if((cmpstr(cmp4name," ") != 0) %& (cmpstr(cmp4name,"") != 0))
	SrwOptContAdd(cmp4name,cntname,1)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Append a Beamline Element to a Container
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptContAdd(radname,groupname,where)
String radname=SrwBliLast+SrwBeamlineType
string groupname=SrwBliCont+SrwBeamlineType
variable where=SrwBliwhere
prompt radname,SrwPBliLast,popup Wavelist("*"+SrwBeamlineType ,";", "")
prompt groupname,SrwPBliCont,popup Wavelist("*"+SrwBeamlineType ,";", "")
prompt where,SrwPBliwhere,popup "Front;End"
Silent 1						|	 ...
PauseUpdate
SrwBliCont=groupname[0,strlen(groupname)-strlen(SrwBeamlineType)-1]
SrwBliLast=SrwBliCont
SrwBliwhere=where

if (cmpstr($groupname[0],SrwBliContType)==1)
abort  "Cannot append to a non Container Structure"
endif

if (cmpstr(groupname,radname)==0)
abort  "Cannot append a container to itself"
endif

variable n= DimSize($groupname,0)
redimension/N=(n+1)  $groupname
$groupname[n]=radname
end

//+++++++++++++++++++++++++++++++++++++++
//
//Container
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptCont(radname)
String radname=SrwBliCont
prompt radname,SrwPBli
Silent 1						|	 ...
PauseUpdate
SrwBliCont=radname;
SrwBliLast=radname

radname+=SrwBeamlineType
Make/T/O/N=1 $radname
$radname[0]=SrwBliContType

end

//+++++++++++++++++++++++++++++++++++++++
//
//Drift space
// length:		distance of the drift in [m]
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptDrift(radname,length)
String radname=SrwBliDrift
variable length=SrwBliDriftLength
prompt radname,SrwPBli
prompt length,SrwPBliDriftLength
Silent 1						|	 ...
PauseUpdate
SrwBliDriftLength=length;
SrwBliDrift=radname;
SrwBliLast=radname

radname+=SrwBeamlineType
Make/T/O/N=2 $radname
$radname[0]=SrwBliDriftType
$radname[1]=num2str(length)

end

//+++++++++++++++++++++++++++++++++++++++
//
//Lens
//
//focalx:		horizontal focal length [m]  ( > 0 if focuising)
//focalz:		vertical focal length [m]  ( > 0 if focuising)
//posx:			horizontal position of the lenth [m]
//posz:			vertical position of the lenth [m]
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinLens(radname,focalx,focalz,posx,posz)
String radname=SrwBliThinLens
variable focalx=SrwBliTlFocalX
variable focalz=SrwBliTlFocalZ
variable posx=SrwBliTlPosX
variable posz=SrwBliTlPosZ
prompt radname,SrwPBli
prompt focalx,SrwPBliTlFocalX
prompt focalz,SrwPBliTlFocalZ
prompt posx,SrwPBliTlPosX
prompt posz,SrwPBliTlPosZ
Silent 1						|	 ...
PauseUpdate
SrwBliThinLens=radname
SrwBliLast=radname
SrwBliTlFocalX=focalx
SrwBliTlFocalZ=focalz
SrwBliTlPosX=posx
SrwBliTlPosZ=posz

posx/=1000
posz/=1000

radname+=SrwBeamlineType
Make/T/O/N=5 $radname
$radname[0]=SrwBliThinLensType
$radname[1]=num2str(focalx)
$radname[2]=num2str(focalz)
$radname[3]=num2str(posx)
$radname[4]=num2str(posz)

end

//+++++++++++++++++++++++++++++++++++++++
//
//Rectangular Aperture
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptApertRect(radname,dx,dz,posx,posz)
String radname=SrwBliRectApert;
variable dx=SrwBliRectApertDx;
variable dz=SrwBliRectApertDz;
variable posx=SrwBliRectApertPosX;
variable posz=SrwBliRectApertPosZ;
prompt radname,SrwPBli;
prompt dx,SrwPBliRectApertDx;
prompt dz,SrwPBliRectApertDz;
prompt posx,SrwPBliRectApertPosX;
prompt posz,SrwPBliRectApertPosZ;
Silent 1						|	 ...
PauseUpdate
SrwBliRectApert=radname;
SrwBliLast=radname;
SrwBliRectApertDx=dx;
SrwBliRectApertDz=dz;
SrwBliRectApertPosX=posx;
SrwBliRectApertPosZ=posz;

dx*=0.001;
dz*=0.001;
posx*=0.001;
posz*=0.001;

radname+=SrwBeamlineType;
Make/T/O/N=5 $radname;
$radname[0]=SrwBliRectApertType;
$radname[1]=num2str(dx);
$radname[2]=num2str(dz);
$radname[3]=num2str(posx);
$radname[4]=num2str(posz);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Return parameters of the Rectangular Aperture structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwGetOptApertRectParam(ApertName, IndParam)
string ApertName
variable IndParam
SVAR SrwBeamlineType, SrwBliRectApertType
srwUtiSetValS("SrwBliRectApert", ApertName[0,strlen(ApertName)-strlen(SrwBeamlineType)-1], "")
wave/T wApert = $ApertName
if(cmpstr(wApert[0], SrwBliRectApertType) != 0)
	abort "This structure is not an optical aperture."
endif
if(IndParam > 4)
	abort "Parameter index value is too big."
endif
return str2num(wApert[IndParam])
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetOptApertRectHorPos(ApertName)
string ApertName
return srwGetOptApertRectParam(ApertName, 3)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetOptApertRectHorSize(ApertName)
string ApertName
return srwGetOptApertRectParam(ApertName, 1)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetOptApertRectVertPos(ApertName)
string ApertName
return srwGetOptApertRectParam(ApertName, 4)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetOptApertRectVertSize(ApertName)
string ApertName
return srwGetOptApertRectParam(ApertName, 2)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Circular Aperture
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptApertCirc(radname,d,posx,posz)
String radname=SrwBliCircApert;
variable d=SrwBliCircApertD;
variable posx=SrwBliCircApertPosX;
variable posz=SrwBliCircApertPosZ;
prompt radname,SrwPBli;
prompt d,SrwPBliCircApertD;
prompt posx,SrwPBliCircApertPosX;
prompt posz,SrwPBliCircApertPosZ;
Silent 1						|	 ...
PauseUpdate
SrwBliCircApert=radname;
SrwBliLast=radname;
SrwBliCircApertD=d;
SrwBliCircApertPosX=posx;
SrwBliCircApertPosZ=posz;

d*=0.001;
posx*=0.001;
posz*=0.001;

radname+=SrwBeamlineType;
Make/T/O/N=5 $radname;
$radname[0]=SrwBliCircApertType;
$radname[1]=num2str(d);
$radname[2]=num2str(posx);
$radname[3]=num2str(posz);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Obstacle: Rectangular
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptObstRect(radname,dx,dz,posx,posz)
string radname=srwUtiGetValS("SrwBliObstRect", "RectObst", "");
variable dx=srwUtiGetValN("SrwBliObstRectDx", 1, "");
variable dz=srwUtiGetValN("SrwBliObstRectDz", 1, "");
variable posx=srwUtiGetValN("SrwBliObstRectPosX", 0, "");
variable posz=srwUtiGetValN("SrwBliObstRectPosZ", 0, "");
prompt radname,"Name of the Rectangular Obstacle Optical Element"  //SrwPBli;
prompt dx,SrwPBliRectApertDx;
prompt dz,SrwPBliRectApertDz;
prompt posx,SrwPBliRectApertPosX;
prompt posz,SrwPBliRectApertPosZ;
Silent 1						|	 ...
PauseUpdate
srwUtiSetValS("SrwBliObstRect", radname, "");
SrwBliLast=radname;
srwUtiSetValN("SrwBliObstRectDx", dx, "");
srwUtiSetValN("SrwBliObstRectDz", dz, "");
srwUtiSetValN("SrwBliObstRectPosX", posx, "");
srwUtiSetValN("SrwBliObstRectPosZ", posz, "");

dx*=0.001;
dz*=0.001;
posx*=0.001;
posz*=0.001;

radname+=SrwBeamlineType;
Make/T/O/N=5 $radname;
$radname[0]=srwUtiGetValS("SrwBliObstRectType", "RectObstacle", "");
$radname[1]=num2str(dx);
$radname[2]=num2str(dz);
$radname[3]=num2str(posx);
$radname[4]=num2str(posz);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Obstacle: Circular
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptObstCirc(radname,d,posx,posz)
string radname=srwUtiGetValS("SrwBliObstCirc", "CircObst", "");
variable d=srwUtiGetValN("SrwBliObstCircD", 1, "");
variable posx=srwUtiGetValN("SrwBliObstCircPosX", 1, "");
variable posz=srwUtiGetValN("SrwBliObstCircPosZ", 1, "");
prompt radname,"Name of the Circular Obstacle Optical Element"; //SrwPBli;
prompt d,SrwPBliCircApertD;
prompt posx,SrwPBliCircApertPosX;
prompt posz,SrwPBliCircApertPosZ;
Silent 1						|	 ...
PauseUpdate
srwUtiSetValS("SrwBliObstCirc", radname, "");
SrwBliLast=radname;
SrwBliCircApertD=d;
srwUtiSetValN("SrwBliObstCircD", d, "");
srwUtiSetValN("SrwBliObstCircPosX", posx, "");
srwUtiSetValN("SrwBliObstCircPosZ", posz, "");

d*=0.001;
posx*=0.001;
posz*=0.001;

radname+=SrwBeamlineType;
Make/T/O/N=5 $radname;
$radname[0]=srwUtiGetValS("SrwBliObstCircType", "CircObstacle", "");
$radname[1]=num2str(d);
$radname[2]=num2str(posx);
$radname[3]=num2str(posz);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Spherical Mirror
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptMirSpher(radname,radius,dx,dz,rotplane,theta,posx,posz)
String radname=SrwBliSpherMirror;
variable radius=SrwBliSpherMirrorR;
variable dx=SrwBliSpherMirrorDx;
variable dz=SrwBliSpherMirrorDz;
variable rotplane=SrwBliSpherMirrorRotPlane;
variable theta=SrwBliSpherMirrorTheta;
variable posx=SrwBliSpherMirrorPosX;
variable posz=SrwBliSpherMirrorPosZ;
prompt radname,SrwPBli;
prompt radius,SrwPBliSpherMirrorR;
prompt dx,SrwPBliSpherMirrorDx;
prompt dz,SrwPBliSpherMirrorDz;
prompt rotplane,SrwPBliSpherMirrorRotPlane,popup "Horizontal;Vertical";
prompt theta,SrwPBliSpherMirrorTheta;
prompt posx,SrwPBliSpherMirrorPosX;
prompt posz,SrwPBliSpherMirrorPosZ;
Silent 1						|	 ...
PauseUpdate
SrwBliSpherMirror=radname;
SrwBliLast=radname;
SrwBliSpherMirrorR=radius;
SrwBliSpherMirrorDx=dx;
SrwBliSpherMirrorDz=dz;
SrwBliSpherMirrorRotPlane=rotplane;
SrwBliSpherMirrorTheta=theta;
SrwBliSpherMirrorPosX=posx;
SrwBliSpherMirrorPosZ=posz;

dx*=0.001;
dz*=0.001;
posx*=0.001;
posz*=0.001;

radname+=SrwBeamlineType;
Make/T/O/N=8 $radname;
$radname[0]=SrwBliSpherMirrorType;
$radname[1]=num2str(radius);
$radname[2]=num2str(dx);
$radname[3]=num2str(dz);
string rotplaneName;
if(rotplane==1)
	rotplaneName="Horizontal"
else
	rotplaneName="Vertical"
endif
$radname[4]=rotplaneName;
$radname[5]=num2str(theta);
$radname[6]=num2str(posx);
$radname[7]=num2str(posz);
end

//+++++++++++++++++++++++++++++++++++++++
//
//Waveguide: Rectangular
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptWgRect(name,len,dx,dz,posx,posz)
String name=SrwBliWgRect
variable len=SrwBliWgRectLen
variable dx=SrwBliWgRectDx
variable dz=SrwBliWgRectDz
variable posx=SrwBliWgRectPosX
variable posz=SrwBliWgRectPosZ
prompt name,SrwPBli
prompt len,SrwPBliWgRectLen
prompt dx,SrwPBliWgRectDx
prompt dz,SrwPBliWgRectDz
prompt posx,SrwPBliWgRectPosX
prompt posz,SrwPBliWgRectPosZ
Silent 1						|	 ...
PauseUpdate

SrwBliWgRect=name
SrwBliLast=name
SrwBliWgRectLen=len
SrwBliWgRectDx=dx
SrwBliWgRectDz=dz
SrwBliWgRectPosX=posx
SrwBliWgRectPosZ=posz

dx*=0.001
dz*=0.001
posx*=0.001
posz*=0.001

name+=SrwBeamlineType
Make/T/O/N=6 $name
$name[0]=SrwBliWgRectType
$name[1]=num2str(len)
$name[2]=num2str(dx)
$name[3]=num2str(dz)
$name[4]=num2str(posx)
$name[5]=num2str(posz)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Zone Plate
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinZonePlate(ElemName,Nzones,Rn,Thick,AttenLen1,RefrDelta1,AttenLen2,RefrDelta2,xc,zc)
String ElemName=SrwBliZonePlate
Variable Nzones=SrwBliZonePlateNzones
Variable Rn=SrwBliZonePlateRn
Variable Thick=SrwBliZonePlateThick
Variable AttenLen1=SrwBliZonePlateAttenLen1
Variable RefrDelta1=SrwBliZonePlateDeltaRefr1
Variable AttenLen2=SrwBliZonePlateAttenLen2
Variable RefrDelta2=SrwBliZonePlateDeltaRefr2
Variable xc=SrwBliThinGenWavePosX
Variable zc=SrwBliThinGenWavePosZ
prompt ElemName,SrwPBli
prompt Thick,"Thickness [mm]"
prompt Nzones,"Total Number of Zones"
prompt Rn,"Outer Zone Radius [mm]"
prompt AttenLen1,"Atten. Length of 1st Material [mm]"
prompt RefrDelta1,"Refr. Index Decr. (1-n) of 1st Material"
prompt AttenLen2,"Atten. Length of 2nd Material [mm]"
prompt RefrDelta2,"Refr. Index Decr. (1-n) of 2nd Material"
prompt xc,SrwPBliThinGenWavePosX
prompt zc,SrwPBliThinGenWavePosZ
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwBliZonePlate=ElemName
SrwBliLast=SrwBliZonePlate
SrwBliZonePlateThick=Thick
SrwBliZonePlateRn=Rn
SrwBliZonePlateNzones=Nzones
SrwBliZonePlateAttenLen1=AttenLen1
SrwBliZonePlateDeltaRefr1=RefrDelta1
SrwBliZonePlateAttenLen2=AttenLen2
SrwBliZonePlateDeltaRefr2=RefrDelta2
SrwBliZonePlateXc=xc
SrwBliZonePlateZc=zc

ElemName+=SrwBeamlineType
Make/T/O/N=10 $ElemName
$ElemName[0]=SrwBliZonePlateType
$ElemName[1]=num2str(Nzones)
$ElemName[2]=num2str(Rn*0.001)
$ElemName[3]=num2str(Thick*0.001)
$ElemName[4]=num2str(AttenLen1*0.001)
$ElemName[5]=num2str(RefrDelta1)
$ElemName[6]=num2str(AttenLen2*0.001)
$ElemName[7]=num2str(RefrDelta2)
$ElemName[8]=num2str(xc*0.001)
$ElemName[9]=num2str(zc*0.001)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Zone Plate Radial Modulation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinZonePlateRadialMod(ZPName,ZoneHeightRatExtToCen,ZoneNum1,ZoneHeightRatZoneNum1ToCen,ZoneNum2,ZoneHeightRatZoneNum2ToCen)
string ZPName=SrwBliLast+SrwBeamlineType
variable ZoneHeightRatExtToCen=srwUtiGetValN("ZoneHeightRatExtToCen", 1, "SrwOptThinZonePlateRadialMod")
variable ZoneNum1=srwUtiGetValN("ZoneNum1", -1, "SrwOptThinZonePlateRadialMod")
variable ZoneHeightRatZoneNum1ToCen=srwUtiGetValN("ZoneHeightRatZoneNum1ToCen", -1, "SrwOptThinZonePlateRadialMod")
variable ZoneNum2=srwUtiGetValN("ZoneNum2", -1, "SrwOptThinZonePlateRadialMod")
variable ZoneHeightRatZoneNum2ToCen=srwUtiGetValN("ZoneHeightRatZoneNum2ToCen", -1, "SrwOptThinZonePlateRadialMod")
prompt ZPName,"Zone Plate Name",popup srwOptElemExistList("ZonePlate") //Wavelist("*"+SrwBeamlineType ,";", "")
prompt ZoneHeightRatExtToCen,"Height Ratio of External Zone to Cen. Zone"
prompt ZoneNum1,"First Intermediate Zone Number"
prompt ZoneHeightRatZoneNum1ToCen,"Height Ratio of 1st Intermed. Zone to Cen."
prompt ZoneNum2,"Second Intermediate Zone Number"
prompt ZoneHeightRatZoneNum2ToCen,"Height Ratio of 2nd Intermed. Zone to Cen."
Silent 1						|	Setting up complex transmission ...
PauseUpdate

SrwBliLast=ZPName[0,strlen(ZPName)-strlen(SrwBeamlineType)-1]
srwUtiSetValN("ZoneHeightRatExtToCen", ZoneHeightRatExtToCen, "SrwOptThinZonePlateRadialMod")
srwUtiSetValN("ZoneNum1", ZoneNum1, "SrwOptThinZonePlateRadialMod")
srwUtiSetValN("ZoneHeightRatZoneNum1ToCen", ZoneHeightRatZoneNum1ToCen, "SrwOptThinZonePlateRadialMod")
srwUtiSetValN("ZoneNum2", ZoneNum2, "SrwOptThinZonePlateRadialMod")
srwUtiSetValN("ZoneHeightRatZoneNum2ToCen", ZoneHeightRatZoneNum2ToCen, "SrwOptThinZonePlateRadialMod")

ReDimension/N=15 $ZPName

$ZPName[10]=num2str(ZoneHeightRatExtToCen)
$ZPName[11]=num2str(ZoneNum1)
$ZPName[12]=num2str(ZoneHeightRatZoneNum1ToCen)
$ZPName[13]=num2str(ZoneNum2)
$ZPName[14]=num2str(ZoneHeightRatZoneNum2ToCen)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Simple Plane Grating
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptGratingPlaneSimple(ElemName,GrooveDens,DispersPlane,AngIns,mOrder,ReflectAvg)
string ElemName=srwUtiGetValS("SrwBliGratPlan", "PlaneGrating", "")
variable GrooveDens=srwUtiGetValN("GrooveDens", 100, "SrwOptGratingPlaneSimple")
variable DispersPlane=srwUtiGetValN("DispersPlane", 1, "SrwOptGratingPlaneSimple")
variable AngIns=srwUtiGetValN("AngIns", 2, "SrwOptGratingPlaneSimple")
variable mOrder=srwUtiGetValN("mOrder", 1, "SrwOptGratingPlaneSimple")
variable ReflectAvg=srwUtiGetValN("ReflectAvg", 1, "SrwOptGratingPlaneSimple")
prompt ElemName,SrwPBli
prompt GrooveDens,"Groove Density [lines/mm]"
prompt DispersPlane,"Dispersion (Deflection) Plane",popup "Horizontal;Vertical"
prompt AngIns,"Angle bw Opt. Axis and Grating Plane [deg.]"
prompt mOrder,"Output Order to be used"
prompt ReflectAvg,"Average Intensity Reflectivity"
Silent 1						|	Setting up complex transmission ...
PauseUpdate

srwUtiSetValS("SrwBliGratPlan", ElemName, "")
SrwBliLast=ElemName

if(GrooveDens <= 0)
	abort "Groove Density should be positive"
endif
srwUtiSetValN("GrooveDens", GrooveDens, "SrwOptGratingPlaneSimple")

if((DispersPlane <= 0) %| (DispersPlane > 2))
	abort "Dispersion Plane can be Horizontal (1) or Vertical (2)"
endif
srwUtiSetValN("DispersPlane", DispersPlane, "SrwOptGratingPlaneSimple")

if((AngIns <= 0) %| (AngIns <= 0))
	abort "Incidence Angle should be between 0 and 90 degrees"
endif
srwUtiSetValN("AngIns", AngIns, "SrwOptGratingPlaneSimple")
srwUtiSetValN("mOrder", mOrder, "SrwOptGratingPlaneSimple")
if(ReflectAvg < 0) %| (ReflectAvg > 1))
	abort "Reflectivity should be between 0 and 1"
endif
srwUtiSetValN("ReflectAvg", ReflectAvg, "SrwOptGratingPlaneSimple")

ElemName += SrwBeamlineType

make/T/O/N=6 $ElemName
$ElemName[0]=srwUtiGetValS("SrwBliGratPlanType", "PlaneGrating", "")
$ElemName[1]=num2str(GrooveDens)
$ElemName[2]=num2str(DispersPlane)
$ElemName[3]=num2str(AngIns)
$ElemName[4]=num2str(mOrder)
$ElemName[5]=num2str(ReflectAvg)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Generate List of Optical Elements of Particular Type
//
//+++++++++++++++++++++++++++++++++++++++
function/S srwOptElemExistList(ElemType)
string ElemType

if(strlen(ElemType) <= 0)
	return ""
endif

SVAR SrwBeamlineType
string AllElemList = Wavelist("*"+SrwBeamlineType,";","")

variable Nc = strlen(AllElemList)
if(Nc <= 1)
	return ""
endif

string NewElemList = "", TestElem = "", CurChar = ""
variable i = 0
do
	CurChar = AllElemList[i]
	if(cmpstr(CurChar, ";") == 0)
		wave/T wTestElem = $TestElem
		if(cmpstr(wTestElem[0], ElemType) == 0)
			NewElemList += TestElem + ";"
		endif
		TestElem = ""
	else
		TestElem += CurChar
	endif
	i += 1
while(i < Nc)
return NewElemList
end

//+++++++++++++++++++++++++++++++++++++++
//
//Returns type name of an optical element
//
//+++++++++++++++++++++++++++++++++++++++
function/S srwOptElemType(wOptElem)
wave/T wOptElem
if(dimsize(wOptElem, 0) > 0)
	return wOptElem[0]
else
	return ""
endif
end

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//Dedicated to the exrtaction of a transmission characteristic of a thin optical element
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proc SrwOptTransmCharDisplay(NameCharact,OptElem,CharType,disp,xc,zc,xr,zr,nx,nz)
string NameCharact=srwUtiGetValS("NameCharact","TransmCharact","SrwOptTransmCharDisplay")
string OptElem=srwUtiGetValS("SrwBliLast","OptElem","")+SrwBeamlineType
variable CharType=srwUtiGetValN("CharType",1,"SrwOptTransmCharDisplay")
variable xc=srwUtiGetValN("xc",0,"SrwOptTransmCharDisplay")
variable xr=srwUtiGetValN("xr",1,"SrwOptTransmCharDisplay")
variable nx=srwUtiGetValN("nx",1,"SrwOptTransmCharDisplay")
variable zc=srwUtiGetValN("zc",0,"SrwOptTransmCharDisplay")
variable zr=srwUtiGetValN("zr",1,"SrwOptTransmCharDisplay")
variable nz=srwUtiGetValN("nz",1,"SrwOptTransmCharDisplay")
variable disp=srwUtiGetValN("disp",1,"SrwOptTransmCharDisplay")
prompt NameCharact,"Name for Characteristic structure"
prompt OptElem,"Optical Element structure",popup wavelist("*"+SrwBeamlineType,";","")
prompt CharType,"Type of Characteristic to extract",popup "Amplitude Transmission;Intensity Transmission;Optical Path Difference"
prompt xc,"Horizontal Center Position [mm]"
prompt xr,"Horizontal Range [mm]"
prompt nx,"Number of Points vs Horizontal Position"
prompt zc,"Vertical Center Position [mm]"
prompt zr,"Vertical Range [mm]"
prompt nz,"Number of Points vs Vertical Position"
prompt disp,"New Display ?",popup "No;Yes"
Silent 1						|	Extracting transmission characteristic of a thin optical element ...
PauseUpdate

if((strlen(NameCharact)==0) %| (cmpstr(NameCharact,"_none_")==0))
	abort "The value of the parameter \"Name for Characteristic structure\" is not defined"
endif
if(strlen(NameCharact) > 28)
	abort "The string parameter \"Name for Characteristic structure\" is too long"
endif
if((strlen(OptElem)==0) %| (cmpstr(OptElem,"_none_")==0))
	abort "The value of the parameter \"Optical Element structure\" is not defined"
endif
if((nx < 1) %| (nz < 1))
	abort "Number of points in horizontal and vertical directions should be positive"
endif
if((cmpstr(srwOptElemType($OptElem), "Drift") == 0) %| (cmpstr(srwOptElemType($OptElem), "Container") == 0))
	abort "Can not extract transmission characteristics of this optical element"
endif
 
srwUtiSetValS("NameCharact",NameCharact,"SrwOptTransmCharDisplay")
//srwUtiSetValS("SrwBliLast",OptElem,"")
SrwBliLast=OptElem[0,strlen(OptElem)-strlen(SrwBeamlineType)-1]
srwUtiSetValN("CharType",CharType,"SrwOptTransmCharDisplay")
srwUtiSetValN("xc",xc,"SrwOptTransmCharDisplay")
srwUtiSetValN("xr",xr,"SrwOptTransmCharDisplay")
srwUtiSetValN("nx",nx,"SrwOptTransmCharDisplay")
srwUtiSetValN("zc",zc,"SrwOptTransmCharDisplay")
srwUtiSetValN("zr",zr,"SrwOptTransmCharDisplay")
srwUtiSetValN("nz",nz,"SrwOptTransmCharDisplay")
srwUtiSetValN("disp",disp,"SrwOptTransmCharDisplay")

string DataUnit = "", CharTypeStr = ""
if(CharType == 1)
	CharTypeStr = "Amplitude Transmission"
endif
if(CharType == 2)
	CharTypeStr = "Intensity Transmission"
endif
if(CharType == 3)
	DataUnit = "m"
	CharTypeStr = "Optical Path Difference"
endif

variable PrintIsNecessary = 0
string LabelBottom = "", LabelLeft = ""
if(nx == 1)
	if(nz == 1)
		NameCharact += "_x"
		make/O/N=(1) $NameCharact
		PrintIsNecessary = 1
	else
		NameCharact += "_z"
		make/O/N=(nz) $NameCharact
		SetScale/I x (zc - 0.5*zr)*0.001, (zc + 0.5*zr)*0.001, "m", $NameCharact
		LabelBottom = "Vertical Position"
		LabelLeft = CharTypeStr
	endif
else
	if(nz == 1)
		NameCharact += "_x"
		make/O/N=(nx) $NameCharact
		SetScale/I x (xc - 0.5*xr)*0.001, (xc + 0.5*xr)*0.001, "m", $NameCharact
		LabelBottom = "Horizontal Position"
		LabelLeft = CharTypeStr
	else
		NameCharact += "_xz"
		make/O/N=(nx, nz) $NameCharact
		SetScale/I x (xc - 0.5*xr)*0.001, (xc + 0.5*xr)*0.001, "m", $NameCharact
		SetScale/I y (zc - 0.5*zr)*0.001, (zc + 0.5*zr)*0.001, "m", $NameCharact
		LabelBottom = "Horizontal Position"
		LabelLeft = "Vertical Position"
	endif
endif
if(strlen(DataUnit) > 0)
	//waveunits($NameCharact, -1)
	SetScale d 0, 1, DataUnit, $NameCharact
endif

if(cmpstr(srwOptElemType($OptElem), "ThinGen") == 0) 

	string NameTransmWave = "Aux_SrwOptTransmCharDisplay"
	string NameOpTranWave = $OptElem[1]
	variable AsFuncOf = 1 //vs x and z
	if((nx == 1) %& (nz == 1))
		AsFuncOf = 2
	endif
	if((nx > 1) %& (nz == 1))
		AsFuncOf = 2
	endif
	if((nx == 1) %& (nz > 1))
		AsFuncOf = 3
	endif
	SrwOptThinTransmDisplay(NameOpTranWave, NameTransmWave, CharType, AsFuncOf, xc, zc, 1)
	
	string NameExtrTot = NameTransmWave
	if((nx == 1) %& (nz == 1))
		NameExtrTot += SrwSeparator+SrwRadXType
		$NameCharact = $NameExtrTot
	endif
	if((nx > 1) %& (nz == 1))
		NameExtrTot += SrwSeparator+SrwRadXType
		$NameCharact = $NameExtrTot(x)
	endif
	if((nx == 1) %& (nz > 1))
		NameExtrTot += SrwSeparator+SrwRadZType
		$NameCharact = $NameExtrTot(x)
	endif
	if((nx > 1) %& (nz > 1))
		NameExtrTot += SrwSeparator+SrwRadXType+SrwRadZType
		$NameCharact = $NameExtrTot(x)(y)
	endif
	killwaves/Z $NameExtrTot
else
	make/O wAuxPrec = {CharType, xc*0.001, xr*0.001, nx, zc*0.001, zr*0.001, nz}
	srOptThinTransmCharExtract($OptElem, wAuxPrec, $NameCharact)
	killwaves/Z wAuxPrec
endif

if(disp == 2)
	if(PrintIsNecessary == 1)
		print $NameCharact[0], DataUnit
	else 
		if((nx > 1) %& (nz > 1))
			display; AppendImage $NameCharact
		else
			display $NameCharact
		endif
		Label bottom LabelBottom
		Label left LabelLeft
	endif
endif
end
