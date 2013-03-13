
//+++++++++++++++++++++++++++++++++++++++
//
//Optical Elements: Mirrors based on Thin Generic 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirInit() 
Variable/G SrwOptThinMirRefl=1
Variable/G SrwOptThinMirAmpTrans
Variable/G SrwOptThinMirRx=1
Variable/G SrwOptThinMirRy=1
Variable/G SrwOptThinMirAng=0
Variable/G SrwOptThinMirAngPl=1
Variable/G SrwOptThinMirHorPos=0 //[m]
Variable/G SrwOptThinMirVerPos=0 //[m]
Variable/G SrwOptThinMirRad=0.01 //[m]
Variable/G SrwOptThinMirApTyp=1
Variable/G SrwOptThinMirHorHalfAp=0.01
Variable/G SrwOptThinMirVerHalfAp=0.01
Variable/G SrwOptThinMirTorOrient=1

Variable/G SrwOptThinMirRh = 1
Variable/G SrwOptThinMirRv = 1
Variable/G SrwOptThinMirPhi0 = 0
Variable/G SrwOptThinMirSizeH = 0.05 //[m]
Variable/G SrwOptThinMirSizeV = 0.05 //[m]
Variable/G SrwOptThinMirCenPosH = 0
Variable/G SrwOptThinMirCenPosV = 0

Variable/G SrwOptThinMirInitPassed = 0
Variable/G SrwOptThinMirSetupPassed = 0

end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror
//
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThinMirTorDialogs()

NVAR SrwOptThinMirInitPassed,SrwOptThinMirSetupPassed

SrwOptThinMirInitPassed=0
execute "SrwOptThinMirTorInit()"
if(SrwOptThinMirInitPassed == 0) 
	return 0
endif
SVAR SrwBliThinGen
NVAR SrwBliThinGenWaveNx, SrwBliThinGenWaveNz, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV
string ComLineStr = ""
sprintf ComLineStr, "SrwOptThinMirTorInit(\"%s\",%g,%g,%g,%g)", SrwBliThinGen, SrwBliThinGenWaveNx, SrwBliThinGenWaveNz, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV

SrwOptThinMirSetupPassed = 0
execute "SrwOptThinMirTorSetup()"
if(SrwOptThinMirSetupPassed == 0) 
	print ComLineStr
	return 0
endif
NVAR SrwOptThinMirRefl, SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirPhi0, SrwOptThinMirTheta0, SrwOptThinMirKsi0, SrwOptThinMirSizeH, SrwOptThinMirSizeV
string ComLineStr2 = ""
sprintf ComLineStr2, "SrwOptThinMirTorSetup(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g)", SrwBliThinGen, SrwOptThinMirRefl, SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirPhi0, SrwOptThinMirTheta0, SrwOptThinMirKsi0, SrwOptThinMirSizeH*1000, SrwOptThinMirSizeV*1000
ComLineStr += (";" + ComLineStr2)
print ComLineStr
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror: Init Setup
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirTorInit(name, npx, npz, xc, yc)
string name=srwUtiGetValS("SrwBliThinGen", "TorMir", "")
variable npx=srwUtiGetValN("SrwBliThinGenWaveNx", 300, "")
variable npz=srwUtiGetValN("SrwBliThinGenWaveNz", 300, "")
variable xc=srwUtiGetValN("SrwOptThinMirCenPosH", 0, "")*1000
variable yc=srwUtiGetValN("SrwOptThinMirCenPosV", 0, "")*1000
prompt name,"Name of the Mirror structure"
prompt npx,"Number of Horizontal Points to represent Transmission"
prompt npz,"Number of Vertical Points to represent Transmission"
prompt xc,"Horizontal Center Position [mm]"
prompt yc,"Vertical Center Position [mm]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThinGen", name, "")
srwUtiSetValN("SrwBliThinGenWaveNx", npx, "")
srwUtiSetValN("SrwBliThinGenWaveNz", npz, "")
srwUtiSetValN("SrwOptThinMirCenPosH", xc*0.001, "")
srwUtiSetValN("SrwOptThinMirCenPosV", yc*0.001, "")
variable/G SrwOptThinMirInitPassed = 1
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirTorSetup(name, refl, rh, rv, phi, theta, ksi, dh, dv)
string name=srwUtiGetValS("SrwBliThinGen", "TorMir", "")
variable refl=srwUtiGetValN("SrwOptThinMirRefl", 1, "")
variable rh=srwUtiGetValN("SrwOptThinMirRh", 1, "")
variable rv=srwUtiGetValN("SrwOptThinMirRv", 1, "")
variable phi=srwUtiGetValN("SrwOptThinMirPhi0", 0, "")
variable theta=srwUtiGetValN("SrwOptThinMirTheta0", 0, "")
variable ksi=srwUtiGetValN("SrwOptThinMirKsi0", 0, "")
variable dh=srwUtiGetValN("SrwOptThinMirSizeH", 0.05, "")*1000
variable dv=srwUtiGetValN("SrwOptThinMirSizeV", 0.05, "")*1000
prompt name,"Name of the Mirror structure"
prompt refl,"Reflectivity"
prompt rh,"Tangential Radius [m]"
prompt rv,"Sagittal (or Minor) Radius [m]"
prompt phi,"Horizontal Angle of Cen. Normal [r]"
prompt theta,"Vertical Angle of Cen. Normal [r]"
prompt ksi,"Rotation Angle about Cen. Normal [r]"
prompt dh,"Mirror Size in Tangential Plane [mm]"
prompt dv,"Mirror Size in Sagittal Plane [mm]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThinGen", name, "")
srwUtiSetValN("SrwOptThinMirRefl", refl, "")
srwUtiSetValN("SrwOptThinMirRh", rh, "")
srwUtiSetValN("SrwOptThinMirRv", rv, "")
srwUtiSetValN("SrwOptThinMirPhi0", phi, "")
srwUtiSetValN("SrwOptThinMirTheta0", theta, "")
srwUtiSetValN("SrwOptThinMirKsi0", ksi, "")
srwUtiSetValN("SrwOptThinMirSizeH", dh*0.001, "")
srwUtiSetValN("SrwOptThinMirSizeV", dv*0.001, "")

variable/G SrwBliThinGenWaveNx, SrwBliThinGenWaveNz, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV
variable Nx = SrwBliThinGenWaveNx, Nz = SrwBliThinGenWaveNz, xc = SrwOptThinMirCenPosH*1000, zc = SrwOptThinMirCenPosV*1000

if(Nx == 0)
	Nx = 400
endif
if(Nz == 0)
	Nz = 400
endif

variable SafetyFact = 1.
variable dd = SafetyFact*sqrt(dh*dh + dv*dv)

string FuncNameAmpTrans = "SrwOptThinAmpTransMirTor", FuncNameOptPath = "SrwOptThinPathMirTor"

killwaves/Z AngAndP, SrwOptThinTrfMatr, SrwOptThinTrfVect, SrwOptThinAuxP
make/O AngAndP = {theta, phi, ksi, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV}, SrwOptAxisVectMirFr = {0, 1, 0}, SrwOptMirNormVect = {0, -1, 0}, SrwOptAxisReflVectMirFr = {0, 0, 0}
make/O/N=(3,3) SrwOptThinTrfMatr
make/O/N=3 SrwOptThinTrfVect, SrwOptThinAuxP

srUtiTrfOptMir(AngAndP, SrwOptThinTrfMatr, SrwOptThinTrfVect)
srwUtiTrfV(SrwOptThinTrfMatr, SrwOptAxisVectMirFr)

variable TwoScalProd = 2*(SrwOptAxisVectMirFr[0]*SrwOptMirNormVect[0] + SrwOptAxisVectMirFr[1]*SrwOptMirNormVect[1] + SrwOptAxisVectMirFr[2]*SrwOptMirNormVect[2])
SrwOptAxisReflVectMirFr[0] = SrwOptAxisVectMirFr[0] - TwoScalProd*SrwOptMirNormVect[0]
SrwOptAxisReflVectMirFr[1] = SrwOptAxisVectMirFr[1] - TwoScalProd*SrwOptMirNormVect[1]
SrwOptAxisReflVectMirFr[2] = SrwOptAxisVectMirFr[2] - TwoScalProd*SrwOptMirNormVect[2]

variable/G SrwOptThinMirTorHalfSizeH = 0.5*SrwOptThinMirSizeH, SrwOptThinMirTorHalfSizeV = 0.5*SrwOptThinMirSizeV, SrwOptThinMirReflAmp = sqrt(abs(SrwOptThinMirRefl))

SrwOptThinTempl(name,1,xc,zc,dd,dd,Nx,Nz)
SrwOptThinSetup(name + SrwOptThinGenWaveType,FuncNameAmpTrans,FuncNameOptPath)

killwaves/Z AngAndP, SrwOptThinTrfMatr, SrwOptThinTrfVect, SrwOptThinAuxP
variable/G SrwOptThinMirSetupPassed = 1
end

//+++++++++++++++++++++++++++++++++++++++
//Toroidal mirror: Amplitude Transmission
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThinAmpTransMirTor(x, z) 
variable x, z
string AuxName = "SrwOptThinAuxP", MatrName = "SrwOptThinTrfMatr", VectName = "SrwOptThinTrfVect", OptAxVectName = "SrwOptAxisVectMirFr"
wave wP = $AuxName, wM = $MatrName, wV = $VectName, wAx = $OptAxVectName

variable/G SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirTorHalfSizeH, SrwOptThinMirTorHalfSizeV, SrwOptThinMirReflAmp

wP[0] = x; wP[1] = 0; wP[2] = z
srwUtiTrfP(wM, wV, wP) 
srwUtiOptIntersTorAndLine(SrwOptThinMirRh, SrwOptThinMirRv, wAx, wP) 
variable xt = wP[0], zt = wP[2]

if((xt < -SrwOptThinMirTorHalfSizeH) %| (xt > SrwOptThinMirTorHalfSizeH))
	return 0
endif
if((zt < -SrwOptThinMirTorHalfSizeV) %| (zt > SrwOptThinMirTorHalfSizeV))
	return 0
endif
return SrwOptThinMirReflAmp
end

//+++++++++++++++++++++++++++++++++++++++
//Toroidal mirror: Optical Path Difference
//+++++++++++++++++++++++++++++++++++++++
function SrwOptThinPathMirTor(x, z) 
variable x, z

string AuxName = "SrwOptThinAuxP", MatrName = "SrwOptThinTrfMatr", VectName = "SrwOptThinTrfVect", OptAxVectName = "SrwOptAxisVectMirFr", OptAxReflVectName = "SrwOptAxisReflVectMirFr"
wave wP = $AuxName, wM = $MatrName, wV = $VectName, wAx = $OptAxVectName, wAxR = $OptAxReflVectName

variable/G SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirTorHalfSizeH, SrwOptThinMirTorHalfSizeV, SrwOptThinMirReflAmp

wP[0] = x; wP[1] = 0; wP[2] = z
srwUtiTrfP(wM, wV, wP) 
variable x0 = wP[0], y0 = wP[1], z0 = wP[2]
srwUtiOptIntersTorAndLine(SrwOptThinMirRh, SrwOptThinMirRv, wAx, wP) 
variable ddx = wP[0] - x0, ddy = wP[1] - y0, ddz = wP[2] - z0
variable Path1 = ddx*wAx[0] + ddy*wAx[1] + ddz*wAx[2]
variable Path2 = -wP[0]*wAxR[0] - wP[1]*wAxR[1] - wP[2]*wAxR[2]

return Path1 + Path2
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror OLD VERSION
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwOptThinMirTor_OLD(name, refl, rh, rv, phi0, orient, dh, dv, xc, yc)
string name=srwUtiGetValS("SrwBliThinGen", "TorMir", "")
variable refl=srwUtiGetValN("SrwOptThinMirRefl", 1, "")
variable rh=srwUtiGetValN("SrwOptThinMirRh", 1, "")
variable rv=srwUtiGetValN("SrwOptThinMirRv", 1, "")
variable phi0=srwUtiGetValN("SrwOptThinMirPhi0", 0, "")
variable orient=srwUtiGetValN("SrwOptThinMirTorOrient", 1, "")
variable dh=srwUtiGetValN("SrwOptThinMirSizeH", 0.05, "")*1000
variable dv=srwUtiGetValN("SrwOptThinMirSizeV", 0.05, "")*1000
variable xc=srwUtiGetValN("SrwOptThinMirCenPosH", 0, "")*1000
variable yc=srwUtiGetValN("SrwOptThinMirCenPosV", 0, "")*1000
prompt name,"Name of the Mirror structure"
prompt refl,"Reflectivity"
prompt rh,"Tangential Radius [m]"
prompt rv,"Sagittal (or Minor) Radius [m]"
prompt phi0,"Angle bw Opt. Axis and Cen. Normal [r]" //and Plane passing through Torus Axis and Central Normal to the reflecting surface [r]"
prompt orient, "Orientation", popup "Tangential Plane Horizontal;Tangential Plane Vertical;Sagittal Plane Horizontal;Sagittal Plane Vertical" //Tilted in Vertical Plane;Tilted in Horizontal Plane"
prompt dh,"Mirror Size in Tangential Plane [mm]"
prompt dv,"Mirror Size in Sagittal Plane [mm]"
prompt xc,"Horizontal Center Position [mm]"
prompt yc,"Vertical Center Position [mm]"
Silent 1						|	 ...
PauseUpdate

srwUtiSetValS("SrwBliThinGen", name, "")
srwUtiSetValN("SrwOptThinMirRefl", refl, "")
srwUtiSetValN("SrwOptThinMirRh", rh, "")
srwUtiSetValN("SrwOptThinMirRv", rv, "")
srwUtiSetValN("SrwOptThinMirPhi0", phi0, "")
srwUtiSetValN("SrwOptThinMirTorOrient", orient, "")
srwUtiSetValN("SrwOptThinMirSizeH", dh*0.001, "")
srwUtiSetValN("SrwOptThinMirSizeV", dv*0.001, "")
srwUtiSetValN("SrwOptThinMirCenPosH", xc*0.001, "")
srwUtiSetValN("SrwOptThinMirCenPosV", yc*0.001, "")

variable Nx = 400, Ny = 400 
variable dx, dz
variable SafetyFact = 1.2
string FuncNameAmpTrans = "AmpTransMirTor", FuncNameOptPath = ""

variable/G SrwOptThinMirTorX0=0, SrwOptThinMirTorY0=0, SrwOptThinMirTorZ0=0, SrwOptThinMirTorLrX=0, SrwOptThinMirTorLrY=0, SrwOptThinMirTorLrZ=0
variable/G SrwOptThinMirTorHalfApX=0, SrwOptThinMirTorHalfApZ=0, SrwOptThinMirReflAmp = sqrt(abs(SrwOptThinMirRefl))

variable CosPhi0 = cos(SrwOptThinMirPhi0), SinPhi0 = sin(SrwOptThinMirPhi0)
if(orient==1) //Torus Axis is oriented VERTICALLY, Optical Axis is deflected in HORIZONTAL plane
	dx = SafetyFact*dh
	dz = SafetyFact*dv
	FuncNameOptPath = "OptPathMirTorOrient1"
	SrwOptThinMirTorX0 = SrwOptThinMirRh*SinPhi0
	SrwOptThinMirTorY0 = -SrwOptThinMirRh*CosPhi0
	SrwOptThinMirTorZ0 = 0
	SrwOptThinMirTorLrX = 2*CosPhi0*SinPhi0
	SrwOptThinMirTorLrY = 1 - 2*CosPhi0*CosPhi0
	SrwOptThinMirTorLrZ = 0
	SrwOptThinMirTorHalfApX = 0.5*SrwOptThinMirSizeH*CosPhi0
	SrwOptThinMirTorHalfApZ = 0.5*SrwOptThinMirSizeV
endif
if(orient==2) //Torus Axis is oriented HORIZONTALLY, Optical Axis is deflected in VERTICAL plane
	dx = SafetyFact*dv
	dz = SafetyFact*dh
	FuncNameOptPath = "OptPathMirTorOrient2"
	SrwOptThinMirTorX0 = 0
	SrwOptThinMirTorY0 = -SrwOptThinMirRh*CosPhi0
	SrwOptThinMirTorZ0 = SrwOptThinMirRh*SinPhi0
	SrwOptThinMirTorLrX = 0
	SrwOptThinMirTorLrY = 1 - 2*CosPhi0*CosPhi0
	SrwOptThinMirTorLrZ = 2*CosPhi0*SinPhi0
	SrwOptThinMirTorHalfApX = 0.5*SrwOptThinMirSizeV
	SrwOptThinMirTorHalfApZ = 0.5*SrwOptThinMirSizeH*CosPhi0
endif
if(orient==3) //Sagittal Plane Horizontal, i.e. Torus Axis is Tilted in HORIZONTAL plane, Optical Axis is deflected in HORIZONTAL plane
	dx = SafetyFact*dv
	dz = SafetyFact*dh
	FuncNameOptPath = "OptPathMirTorOrient4"
	SrwOptThinMirTorX0 = SrwOptThinMirRh*SinPhi0
	SrwOptThinMirTorY0 = -SrwOptThinMirRh*CosPhi0
	SrwOptThinMirTorZ0 = 0
	SrwOptThinMirTorLrX = 2*CosPhi0*SinPhi0
	SrwOptThinMirTorLrY = 1 - 2*CosPhi0*CosPhi0
	SrwOptThinMirTorLrZ = 0
	SrwOptThinMirTorHalfApX = 0.5*SrwOptThinMirSizeV*CosPhi0
	SrwOptThinMirTorHalfApZ = 0.5*SrwOptThinMirSizeH
endif
if(orient==4) //Sagittal Plane Vertical, i.e. Torus Axis is Tilted in VERTICAL plane, Optical Axis is deflected in VERTICAL plane
	dx = SafetyFact*dv
	dz = SafetyFact*dh
	FuncNameOptPath = "OptPathMirTorOrient3"
	SrwOptThinMirTorX0 = 0
	SrwOptThinMirTorY0 = -SrwOptThinMirRh*CosPhi0
	SrwOptThinMirTorZ0 = SrwOptThinMirRh*SinPhi0
	SrwOptThinMirTorLrX = 0
	SrwOptThinMirTorLrY = 1 - 2*CosPhi0*CosPhi0
	SrwOptThinMirTorLrZ = 2*CosPhi0*SinPhi0
	SrwOptThinMirTorHalfApX = 0.5*SrwOptThinMirSizeH
	SrwOptThinMirTorHalfApZ = 0.5*SrwOptThinMirSizeV*CosPhi0
endif

SrwOptThinTempl(name,1,xc,yc,dx,dz,Nx,Ny)
SrwOptThinSetup(name + SrwOptThinGenWaveType,FuncNameAmpTrans,FuncNameOptPath)
end

//+++++++++++++++++++++++++++++++++++++++
function OptPathMirTorOrient1(x, z)  //Torus Axis is oriented VERTICALLY, Optical Axis is deflected in HORIZONTAL plane
variable x, z

Variable/G SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV
Variable/G SrwOptThinMirTorLrX, SrwOptThinMirTorLrY, SrwOptThinMirTorLrZ
Variable/G SrwOptThinMirTorX0, SrwOptThinMirTorY0, SrwOptThinMirTorZ0

variable RelX = x - SrwOptThinMirCenPosH - SrwOptThinMirTorX0
variable RelZ = z - SrwOptThinMirCenPosV
variable rve2 = SrwOptThinMirRv*SrwOptThinMirRv
variable AuxR = (SrwOptThinMirRh - SrwOptThinMirRv) + sqrt(rve2 - RelZ*RelZ)
variable y = SrwOptThinMirTorY0 + sqrt(AuxR*AuxR - RelX*RelX)
variable dsr = -x*SrwOptThinMirTorLrX - y*SrwOptThinMirTorLrY - z*SrwOptThinMirTorLrZ
variable OptPath = y + dsr
return OptPath
end

//+++++++++++++++++++++++++++++++++++++++
function OptPathMirTorOrient2(x, z)  //Torus Axis is oriented HORIZONTALLY, Optical Axis is deflected in VERTICAL plane
variable x, z

NVAR SrwOptThinMirRh, SrwOptThinMirRv, SrwOptThinMirCenPosH, SrwOptThinMirCenPosV
NVAR SrwOptThinMirTorLrX, SrwOptThinMirTorLrY, SrwOptThinMirTorLrZ
NVAR SrwOptThinMirTorX0, SrwOptThinMirTorY0, SrwOptThinMirTorZ0

//variable RelX = x - SrwOptThinMirCenPosH - SrwOptThinMirTorX0
variable RelZ = z - SrwOptThinMirCenPosV - SrwOptThinMirTorZ0
//variable RelZ = z - SrwOptThinMirCenPosV
variable RelX = x - SrwOptThinMirCenPosH
variable rve2 = SrwOptThinMirRv*SrwOptThinMirRv
variable AuxR = (SrwOptThinMirRh - SrwOptThinMirRv) + sqrt(rve2 - RelX*RelX)
variable y = SrwOptThinMirTorY0 + sqrt(AuxR*AuxR - RelZ*RelZ)
variable dsr = -x*SrwOptThinMirTorLrX - y*SrwOptThinMirTorLrY - z*SrwOptThinMirTorLrZ
variable OptPath = y + dsr
return OptPath
end

//+++++++++++++++++++++++++++++++++++++++
//
//Toroidal mirror @ Normal incidence
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwOptThinMirTorNorm(name, refl, rx, ry, aptyp, d, dx, dy, xc, yc)
proc SrwOptThinMirTorNorm(name, refl, rx, ry, dx, dy, orient, xc, yc)
String name=SrwBliThinGen
Variable refl=SrwOptThinMirRefl
Variable rx=SrwOptThinMirRx
Variable ry=SrwOptThinMirRy
Variable dx=SrwOptThinMirHorHalfAp*2*1000
Variable dy=SrwOptThinMirVerHalfAp*2*1000
Variable orient=SrwOptThinMirTorOrient
Variable xc=SrwOptThinMirHorPos*1000
Variable yc=SrwOptThinMirVerPos*1000
prompt name,"Name of the Mirror structure"
prompt refl,"Reflection Coefficient"
prompt rx,"Horizontal Radius of Curvature [m]"
prompt ry,"Vertical Radius of Curvature [m]"
prompt dx,"Hor. Apert. Size [mm]"
prompt dy,"Ver. Apert. Size [mm]"
prompt orient,"Median Plane Orientation",popup "Horizontal;Vertical"
prompt xc,"Horizontal Center Position [m]"
prompt yc,"Vertical Center Position [m]"
Silent 1						|	 ...
PauseUpdate

SrwBliThinGen=name
SrwBliLast=name
SrwOptThinMirRefl=refl
SrwOptThinMirAmpTrans=sqrt(refl)
SrwOptThinMirRx=rx
SrwOptThinMirRy=ry
SrwOptThinMirHorPos=xc*0.001
SrwOptThinMirVerPos=yc*0.001
SrwOptThinMirTorOrient = orient
SrwOptThinMirHorHalfAp = 0.5*dx*0.001
SrwOptThinMirVerHalfAp = 0.5*dy*0.001

Variable re2 = SrwOptThinMirHorHalfAp*SrwOptThinMirHorHalfAp + SrwOptThinMirVerHalfAp*SrwOptThinMirVerHalfAp

//Variable dz1 = -ry + sqrt(ry*ry - re2)
//Variable dz2 = -rx + sqrt(rx*rx - re2)

//if(rx < ry)
//	SrwOptThinMirPlaneOffset = -rx + sqrt(rx*rx - re2)
//else
//	SrwOptThinMirPlaneOffset = -ry + sqrt(ry*ry - re2)
//endif

Variable BufVar
if(orient == 1)
	BufVar = rx - ry + sqrt(ry*ry - SrwOptThinMirVerHalfAp*SrwOptThinMirVerHalfAp)
	SrwOptThinMirPlaneOffset = -rx + sqrt(BufVar*BufVar - SrwOptThinMirHorHalfAp*SrwOptThinMirHorHalfAp)
else
	BufVar = ry - rx + sqrt(rx*rx - SrwOptThinMirHorHalfAp*SrwOptThinMirHorHalfAp)
	SrwOptThinMirPlaneOffset = -ry + sqrt(BufVar*BufVar - SrwOptThinMirVerHalfAp*SrwOptThinMirVerHalfAp)
endif

Variable Nx = 300, Ny = 300
Variable AddRngX = 2*abs(SrwOptThinMirHorPos), AddRngY = 2*abs(SrwOptThinMirVerPos)

SrwOptThinTempl(name,2,0,0,dx*1.1 + AddRngX,dy*1.1 + AddRngY,Nx,Ny)
SrwOptThinSetup(name + SrwOptThinGenWaveType,"AmpTransMirTorNorm","OptPathMirTorNorm")

end

//+++++++++++++++++++++++++++++++++++++++
function OptPathMirTorNorm(x, y)
Variable x, y

Variable/G SrwOptThinMirRx, SrwOptThinMirRy, SrwOptThinMirPlaneOffset
Variable/G SrwOptThinMirHorPos, SrwOptThinMirVerPos, SrwOptThinMirTorOrient

Variable RelX = x - SrwOptThinMirHorPos, RelY = y - SrwOptThinMirVerPos

Variable Aux1, CurPath
if(SrwOptThinMirTorOrient == 1)
	Aux1 = SrwOptThinMirRx - SrwOptThinMirRy + sqrt(SrwOptThinMirRy*SrwOptThinMirRy - RelY*RelY)
	CurPath = 2*((-SrwOptThinMirRx + sqrt(Aux1*Aux1 - RelX*RelX)) - SrwOptThinMirPlaneOffset)
else // median plane oriented vertically (check this !)
	Aux1 = SrwOptThinMirRy - SrwOptThinMirRx + sqrt(SrwOptThinMirRx*SrwOptThinMirRx - RelX*RelX)
	CurPath = 2*((-SrwOptThinMirRy + sqrt(Aux1*Aux1 - RelY*RelY)) - SrwOptThinMirPlaneOffset)
endif
// Phase flip Pi is not taken into account
return CurPath
end

//+++++++++++++++++++++++++++++++++++++++
function AmpTransMirTorNorm(x, y)
Variable x, y
//Variable/G SrwOptThinMirRad, SrwOptThinMirAmpTrans
Variable/G SrwOptThinMirAmpTrans
Variable/G SrwOptThinMirHorPos, SrwOptThinMirVerPos
//Variable/G SrwOptThinMirApTyp, SrwOptThinMirHorHalfAp, SrwOptThinMirVerHalfAp
Variable/G SrwOptThinMirHorHalfAp, SrwOptThinMirVerHalfAp

Variable RelX = x - SrwOptThinMirHorPos, RelY = y - SrwOptThinMirVerPos
Variable xe2 = RelX*RelX, ye2 = RelY*RelY

//if(SrwOptThinMirApTyp == 1) //Circular
//	if(xe2 + ye2 > SrwOptThinMirRad*SrwOptThinMirRad)
//		return 0
//	endif
//else //Rectangular

if((RelX < -SrwOptThinMirHorHalfAp) %| (RelX > SrwOptThinMirHorHalfAp))
	return 0
endif
if((RelY < -SrwOptThinMirVerHalfAp) %| (RelY > SrwOptThinMirVerHalfAp))
	return 0
endif
//endif
return SrwOptThinMirAmpTrans
end