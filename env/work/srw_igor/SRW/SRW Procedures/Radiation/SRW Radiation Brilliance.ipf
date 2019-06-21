
//+++++++++++++++++++++++++++++++++++++++
//
// Routines to estimate Brilliance
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwBrilInit()

String/G SrwBrilName=""; String/G SrwPBrilName="Name of the Brillance curve"
Variable/G SrwBrilHarm=1; String/G SrwPBrilHarm="Harmonic"
Variable/G SrwBrilHarmMin=1; String/G SrwPBrilHarmMin="Initial Harmonic"
Variable/G SrwBrilHarmMax=5; String/G SrwPBrilHarmMax="Final Harmonic"

Variable/G SrwKmin=0.3; String/G SrwPKmin="Minimum Deflection Parameter" 
Variable/G SrwNbEnpts=100; String/G SrwPNbEnpts="Number of Energy Points" 

Variable/G SrwType=3; String/G SrwPSrwType="Type" 
Variable/G SrwPlot=2; String/G SrwPPlot="Display" 
String/G SrwBrStr="B",SrwFlStr="F",SrwAfStr="A",SrwEnStr="E",SrwSepStr=".",SrwGlStr="H"

String/G SrwPFIeld="Vertical Field [T]" 

//Variable/G SrwFIeld=0.85
//Variable/G SrwEnMin=1; String/G SrwPEnMin="Initial Energy [keV]" 
//Variable/G SrwEnMax=100; String/G SrwPEnMax="Final Energy [keV]" 
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwBrilUndHarm(BrilName,ElecName,MagName,Kmin,Harm,NbEnpts,Type,Plot)
String BrilName=SrwElecName+SrwUndName;
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwUndName+SrwUndType;
Variable Kmin=SrwKmin;
Variable Harm=SrwBrilHarm;
Variable NbEnpts=SrwNbEnpts;
Variable Type=SrwType;
Variable Plot=SrwPlot;
prompt BrilName,SrwPBrilName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "");
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "");
prompt Kmin,SrwPKmin;
prompt Harm,SrwPBrilHarm;
prompt NbEnpts,SrwPNbEnpts;
prompt Type,SrwPSrwType,popup "Phot/s/.1%;Phot/s/.1%/mr2;Phot/s/.1%/mr2/mm2"
prompt Plot,SrwPPlot,popup "No;Yes"
Silent 1						|	Computing the Brilliance or Flux or ...
PauseUpdate

Harm=round(Harm*0.50001)*2-1

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwKmin=Kmin
SrwBrilHarm=Harm
SrwNbEnpts=NbEnpts
SrwType=Type
SrwPlot=Plot
SrwBrilName=BrilName

//  Get Kx and Kz 
Variable Kz=0,Kx=0,Per,N,Phase=0
if ($($MagName[6])[1]==1)
Kz=$($MagName[6])[2]
else
Kx=$($MagName[6])[2]
endif

if (str2num($MagName[5])==2)
if ($($MagName[7])[1]==1)
Kz=$($MagName[7])[2]
Phase=Phase+$($MagName[7])[3]
else
Kx=$($MagName[7])[2]
Phase=Phase-$($MagName[7])[3]
endif
endif
Per=str2num($MagName[0])
N=str2num($MagName[1])/Per

//Print Kx,Kz,Per,N,Phase
// Get electron Beam Parameters

Variable en,cur,sigx,sigpx,sigz,sigpz,mx,mz
en=$ElecName[0];cur=$ElecName[1]
sigx=$ElecName[20];sigpx=$ElecName[22];mx=$ElecName[21]
sigz=$ElecName[23];sigpz=$ElecName[25];mz=$ElecName[24]

String/G Yy,Xx

if(Type==1)
	Yy=SrwFlStr
endif
if(Type==2)
	Yy=SrwAfStr
endif
if(Type==3)
	Yy=SrwBrStr
endif
	
Yy=BrilName+num2str(Harm)+SrwSepStr+Yy
Xx=BrilName+num2str(Harm)+SrwSepStr+SrwEnStr
Make/N=(NbEnpts)/D/O $Yy $Xx temp
Variable K2=Kx*Kx+Kz*Kz
//Variable d=sqrt(Kx^4+Kz^4+2*Kx^2*Kz^2*cos(Phase))*Harm
Variable d=sqrt(Kx^4+Kz^4+2*Kx^2*Kz^2*sin(Phase))*Harm
Variable h1=(Harm-1)/2
Variable h2=(Harm+1)/2
Variable L=N*Per,cst
SetScale/I x 1,kmin^2/K2,"", temp
temp=x
//print d,harm,h1,h2,k2

// Energy and Flux
$Xx=9.5*en*en/Per/(1+temp*k2/2)*Harm
$Yy=1.431E+14*N*cur*temp*k2/(1+temp*k2/2)*Harm*(bessJ(h1,temp*d/(4+temp*2*K2))-bessJ(h2,temp*d/(4+2*temp*K2)))^2

// Angular Flux
if(Type==2)
	$Yy*=1.744E*N*en^2/1.431*N*harm/(1+temp*k2/2)
	$Yy/=sqrt(1+sigpx/(12.4/L/$Xx*1e-7))*sqrt(1+sigpz/(12.4/L/$Xx*1e-7))
endif
// Brilliance
if(Type==3)
	cst=(2*pi)^2*1e12
	$Yy/=cst*sqrt((sigx+12.4/16/pi/pi*L/$Xx*1e-7)*(sigpx+12.4/L/$Xx*1e-7)-mx*mx)
	$Yy/=sqrt((sigz+12.4/16/pi/pi*L/$Xx*1e-7)*(sigpz+12.4/L/$Xx*1e-7)-mz*mz)
endif
SetScale/P y 0,1,"eV", $Xx

//print sqrt(sigpx)*1e6
//print sqrt(sigpz)*1e6

if(Plot==2)
	Display $Yy vs $Xx
endif
killwaves/Z temp
end 

//+++++++++++++++++++++++++++++++++++++++
proc SrwBrilUnd(BrilName,ElecName,MagName,Kmin,HarmMin,HarmMax,NbEnpts,Type,Plot)
String BrilName=SrwElecName+SrwUndName;
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwUndName+SrwUndType;
Variable Kmin=SrwKmin;
Variable HarmMin=SrwBrilHarmMin;
Variable HarmMax=SrwBrilHarmMax;
Variable NbEnpts=SrwNbEnpts;
Variable Type=SrwType;
Variable Plot=SrwPlot;
prompt BrilName,SrwPBrilName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "");
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "");
prompt Kmin,SrwPKmin;
prompt HarmMin,SrwPBrilHarmMin;
prompt HarmMax,SrwPBrilHarmMax;
prompt NbEnpts,SrwPNbEnpts;
prompt Type,SrwPSrwType,popup "Phot/s/.1%;Phot/s/.1%/mr2;Phot/s/.1%/mr2/mm2"
prompt Plot,SrwPPlot,popup "No;Yes"

variable ElecWavePresent = 1, MagWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if(MagWavePresent == 1)
		SrwBrilUnd()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	SrwBrilUnd()
	Return
endif
if(SrwUndIsEllips(MagName)==0)
	Abort "Sorry, this type of computation supports only sinusoidal magnetic field."
endif
if(HarmMin>HarmMax)
	Abort "Wrong Harmonic Numbers"
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwKmin=Kmin
SrwBrilHarmMin=HarmMin
SrwBrilHarmMax=HarmMax
SrwNbEnpts=NbEnpts
SrwType=Type
SrwBrilName=BrilName

HarmMin=Max(HarmMin,1)
variable h=HarmMin

if(Plot==2)
	Display 
endif
do
	SrwBrilUndHarm(BrilName, ElecName, MagName, Kmin,h,NbEnpts,Type,1)
	if(Plot==2)
		AppendToGraph $Yy vs $Xx
	endif
	h+=2
while (h<harmMax+1)					// as long as condition is true
if(Plot==2)
	ModifyGraph log(left)=1
	ModifyGraph log=1
	Label bottom SrwPLabelPhotEn

	if(Type==1)
		Label left SrwPUnitSpAngFlux
	endif
	if(Type==2)
		Label left SrwPUnitSpAngFluxPerUnAngle
	endif
	if(Type==3)
		Label left SrwPUnitBrilliance
	endif
endif

SrwPlot=Plot
end

//+++++++++++++++++++++++++++++++++++++++
function SrwUndIsEllips(MagName)
String MagName
Wave/T w = $MagName	

if(exists(MagName)==0)
	return 0
endif
Variable Nharm = str2num(w[5])
if(Nharm > 2)
	return 0
endif
if(Nharm==2) 
	Wave wh1=$(w[6])
	Wave wh2=$(w[7])
	
	Variable PlaneInd1 = wh1[1]
	Variable PlaneInd2 = wh2[1]
	if(PlaneInd1 == PlaneInd2) // two harmonics are only accepted if they describe field in different planes
		return 0
	endif
	
	//Variable DelPhi=Abs(wh1[3] - wh2[3])
	//if(Abs(DelPhi - 1.571) > 0.05)
	//	return 0
	//endif
endif
return 1
end

//+++++++++++++++++++++++++++++++++++++++
function SrwUndIsPlanar(MagName)
String MagName
Wave/T w = $MagName

if(exists(MagName)==0)
	return 0
endif
Variable Nharm = str2num(w[5])
if(Nharm > 1)
	return 0
endif
return 1
end

//+++++++++++++++++++++++++++++++++++++++
// Returns estimate of RMS Angular Divergence of Bending Magnet Radiation
//+++++++++++++++++++++++++++++++++++++++
function srwUtiAngDivBM(en_d_enc, inv_gam)
variable en_d_enc, inv_gam

variable b0 = -0.57988 //coefficients found by fitting
variable b1 = -0.45258
variable b2 = -0.012217
variable b3 = 0.00027521
variable b4 = 6.1404e-005

variable xt = ln(en_d_enc)
return inv_gam*exp(b0 + (b1 + (b2 + (b3 + b4*xt)*xt)*xt)*xt)
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwBrilBM(BrilName,ElecName,Field,NbEnpts,EnMin,EnMax,Type,Plot)
String BrilName=SrwElecName+num2str(SrwMagConstBz)
String ElecName=SrwElecName+SrwElecType
Variable Field=SrwMagConstBz
Variable EnMin=SrwSmpEdep
Variable EnMax=SrwSmpEfin
Variable NbEnpts=SrwNbEnpts
Variable Type=SrwType
Variable Plot=SrwPlot
prompt BrilName,SrwPBrilName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt Field,SrwPField
prompt NbEnpts,SrwPNbEnpts
prompt EnMin,SrwPSmpEdep
prompt EnMax,SrwPSmpEfin
prompt Type,SrwPSrwType,popup "Phot/s/.1%/mr;Phot/s/.1%/mr2;Phot/s/.1%/mr2/mm2;W/mr (integ. over photon energy)" //;W/mr2 (--//--);W/mr2/mm2 (--//--)"
prompt Plot,SrwPPlot,popup "No;Yes"

Silent 1						|	Computing Brilliance ...
PauseUpdate

if(cmpstr(ElecName,"_none_")==0)
	SrwElecFilament()
	SrwElecThick()
	SrwBrilBM()
	Return
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwMagConstBz=Field
SrwSmpEdep=EnMin
SrwSmpEfin=EnMax
SrwNbEnpts=NbEnpts
SrwType=Type
SrwBrilName=BrilName

Variable EnMin_=1000*EnMin
Variable EnMax_=1000*EnMax
String/G Yy

if(Type==1)
	Yy=SrwGlStr
endif
if(Type==2)
	Yy=SrwAfStr
endif
if(Type==3)
	Yy=SrwBrStr
endif
if(Type==4)
	NbEnpts = 1
	Yy = SrwGlStr
endif
if(Type==5)
	NbEnpts = 1
	Yy = SrwAfStr
endif
if(Type==6)
	NbEnpts = 1
	Yy = SrwBrStr
endif

Yy = BrilName+SrwSepStr+Yy
Make/N=(NbEnpts)/D/O $Yy

if(NbEnpts > 1)
	SetScale/I x EnMin_,EnMax_,"eV", $Yy
else
	SetScale/P x EnMin_,0.01*EnMin_,"eV", $Yy
endif

Variable en,cur,sigx,sigpx,sigz,sigpz
en=$ElecName[0] // E [GeV]
cur=$ElecName[1] // I [A]
sigx=$ElecName[20] // <(x-<x>)^2> [m^2]
sigpx=$ElecName[22] // <(x'-<x'>)^2> [r^2]
sigz=$ElecName[23] // <(z-<z>)^2> [m^2]
sigpz=$ElecName[25] // <(z'-<z'>)^2> [r^2]

if((sigx <= 0.) %| (sigz <= 0))
	Abort SrwPAlertElecThickPar
endif

variable cst,ec,cst1,cst2,cst3
ec=0.665*en*en*field*1000
cst=1.327E13*en*en*cur
cst1=0.107/en/en*1e-6
cst2=1/(2*pi*sqrt(sigx*sigz))*1e-6
cst3=2.457E13*en*cur

variable cwpkev=8*1.60219e-13
variable eckev=0.001*ec
variable invGamma = 0.510998902e-03/en

variable RelPrecNumInt = 0.001

if(Type==1) // Phot/s/.1%/mr
	$Yy = cst3*(x/ec)*srKn(1,5/3,x/ec)
endif
if(Type==2) // Phot/s/.1%/mr2
	//$Yy=cst*(x/ec)*srKn(0,2/3 ,(x/2/ec))/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))
	//$Yy=cst*((x/ec)*srKn(0,2/3,(x/2/ec)))^2/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))
	$Yy = cst*((x/ec)*srKn(0,2/3,(x/2/ec)))^2/sqrt(1 + sigpz/(srwUtiAngDivBM(x/ec,invGamma))^2) //OC010609	
endif
if(Type==3) // Phot/s/.1%/mr2/mm2
	//$Yy=cst*(x/ec)*srKn(0,2/3 ,(x/2/ec))/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))*cst2
	//$Yy=cst*((x/ec)*srKn(0,2/3,(x/2/ec)))^2/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))*cst2
	$Yy = cst*((x/ec)*srKn(0,2/3,(x/2/ec)))^2/sqrt(1 + sigpz/(srwUtiAngDivBM(x/ec,invGamma))^2)*(1e-6/(2*Pi))/sqrt((sigx + (1.239842e-06/x/(4*Pi)/srwUtiAngDivBM(x/ec,invGamma))^2)*(sigz + (1.239842e-06/x/(4*Pi)/srwUtiAngDivBM(x/ec,invGamma))^2))  //OC010609
endif

if(Type==4) // W/mr (integ. over photon energy)
	$Yy = cwpkev*cst3*srUtiIntKnXn(2, EnMin/eckev, EnMax/eckev, RelPrecNumInt)
endif
//Needs checking/correction:
//if(Type==5) // W/mr2
//	$Yy = cwpkev*cst*eckev*srUtiIntKnXn(1, 0.5*EnMin/eckev, 0.5*EnMax/eckev, RelPrecNumInt)/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))
//endif
//if(Type==6) // W/mr2/mm2
//	$Yy = cwpkev*cst*cst2*eckev*srUtiIntKnXn(1, 0.5*EnMin/eckev, 0.5*EnMax/eckev, RelPrecNumInt)/sqrt(1+sigpz/(cst1*(x/ec)^(-1.1)))
//endif

if(Plot==2)
	SrwBrilWigPlot(Yy,Type)
endif

SrwPlot=Plot
end

//+++++++++++++++++++++++++++++++++++++++
proc SrwBrilWigPlot(wav,Type)
string wav
variable Type

variable Np = dimsize($Yy, 0)
string BrilUnits = ""
if(Type==1)
	BrilUnits = SrwPUnitSpAngFluxPerUnPlAng
endif
if(Type==2)
	BrilUnits = "Phot/s/.1%/mr2" //SrwPUnitSpAngFluxPerUnAngle
endif
if(Type==3)
	BrilUnits = "Phot/s/.1%/mr2/mm2" //SrwPUnitBrilliance
endif
if(Type==4)
	BrilUnits = "W/mr"
endif
if(Type==5)
	BrilUnits = "W/mr2"
endif
if(Type==6)
	BrilUnits = "W/mr2/mm2"
endif

if(Np > 1)
	display $Yy
	Label bottom SrwPLabelPhotEn
	Label left BrilUnits
else
	print $Yy[0], BrilUnits
endif
end

//++++++++++++++++++++++++++++++++++++++
proc SrwBrilWig(BrilName,ElecName,MagName,NbEnpts,EnMin,EnMax,Type,Plot)
String BrilName=SrwElecName+SrwUndName
String ElecName=SrwElecName+SrwElecType
String MagName=SrwUndName+SrwUndType
Variable EnMin=SrwSmpEdep
Variable EnMax=SrwSmpEfin
Variable NbEnpts=SrwNbEnpts
Variable Type=SrwType
Variable Plot=SrwPlot
prompt BrilName,SrwPBrilName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt MagName,SrwPUndName2,popup Wavelist("*"+SrwUndType ,";", "")
prompt EnMin,SrwPSmpEdep
prompt EnMax,SrwPSmpEfin
prompt NbEnpts,SrwPNbEnpts
prompt Type,SrwPSrwType,popup "Phot/s/.1%/mr;Phot/s/.1%/mr2;Phot/s/.1%/mr2/mm2;W/mr (integ. vs photon energy);W/mr2 (--//--);W/mr2/mm2 (--//--)"
prompt Plot,SrwPPlot,popup "No;Yes"

Silent 1						|	Computing the Brilliance or Flux or  ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0
endif
if(ElecWavePresent==0)
	SrwElecFilament()
	SrwElecThick()
	if(MagWavePresent==1)
		SrwBrilWig()
		Return
	endif
endif
if(MagWavePresent==0)
	SrwMagPerCreate2D()
	SrwBrilWig()
	Return
endif

if(SrwUndIsPlanar(MagName)==0)
	Abort "Sorry, this type of computation supports only planar sinusoidal magnetic field."
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwUndName=MagName[0,strlen(MagName)-strlen(SrwUndType)-1]
SrwSmpEdep=EnMin
SrwSmpEfin=EnMax
SrwNbEnpts=NbEnpts
SrwType=Type
SrwBrilName=BrilName

//get field and number of periods
Variable Kz=0,per,N,field,sigx,en
if($($MagName[6])[1]==1)
	Kz=$($MagName[6])[2]
endif
sigx=$ElecName[20]
en=$ElecName[0]

per = str2num($MagName[0])
N = str2num($MagName[1])/per

field = Kz/0.0934/per/1000
SrwBrilBM(BrilName,ElecName,field,NbEnpts,EnMin,EnMax,Type,1)

$Yy *= 2*N

variable Buf1
if(Type==3)
	Buf1 = sqrt(1 + (per/2/pi/1957/en*Kz)^2/sigx)
	$Yy /= Buf1
endif

if(Plot==2)
	SrwBrilWigPlot(Yy,Type)
endif

SrwPlot=Plot
end
