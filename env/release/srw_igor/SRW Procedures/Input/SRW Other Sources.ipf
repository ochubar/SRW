
//+++++++++++++++++++++++++++++++++++++++
//
// Isotropic Light Source: Create Wavefront
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrIsotrSrcCreate(RadName,ElecName,Phot,Polar,ObsName,ObsNxNzForProp,ObsNxNzSamplFact)
String RadName=SrwElecName+SrwSmpName
String ElecName=SrwElecName+SrwElecType
Variable Phot=SrwIsotrSrcPhot
Variable Polar=SrwGsnBeamPolar
String ObsName=SrwSmpName+SrwSmpType
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt RadName,SrwPRadName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "")
prompt Phot,SrwPIsotrSrcPhot
prompt Polar,SrwPGsnBeamPolar,popup SrwPOPUPPolar
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	Creating the Wavefront  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

if(cmpstr(ElecName,"_none_")==0)
	Abort SrwPAlertElecBeamNeeded
endif
if(cmpstr(ObsName,"_none_")==0)
	Abort SrwPAlertRadSmplNeeded
endif
if(Phot<=0)
	Abort SrwPAlertIsotrSrcPhot
endif
if((Polar<1) %| (Polar>6))
	Abort SrwPAlertGsnBeamPolar
endif

if(ObsNxNzForProp<1)
	ObsNxNzForProp=1
endif
if(ObsNxNzForProp>2)
	ObsNxNzForProp=2
endif
if(ObsNxNzSamplFact<=0.)
	Abort SrwPAlertBadNxNzSamplFact
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwIsotrSrcPhot=Phot
SrwGsnBeamPolar=Polar
SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact
if(strlen(RadName)==0)
	RadName=SrwElecName+SrwMagName+SrwSmpName
endif
SrwRadName=RadName

SrwWfrPrep(ElecName,ObsName,RadName,0)
RadName += SrwRadType
SrwRadGenTotName=RadName

Make/D/O/N=5 AuxObsTreat; AuxObsTreat[0]=ObsNxNzForProp-1; AuxObsTreat[1]=ObsNxNzSamplFact
Make/D/O/N=5 AuxSrcParam; AuxSrcParam[0]=Phot; AuxSrcParam[1]=Polar

srWfrIsotrSrc($ElecName,AuxSrcParam,$ObsName,AuxObsTreat,$RadName)
KillWaves/Z AuxSrcParam,AuxObsTreat

end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: Init globals
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeamInit()

String/G SrwGsnBeamType=SrwSeparator+"gsb"
String/G SrwGsnBeamName="Gsn"; String/G SrwPGsnBeamName="Name of the Source structure", SrwPGsnBeamName2="Gaussian Beam structure"
Variable/G SrwGsnBeamWaistX=1; String/G SrwPGsnBeamWaistX="Horizontal RMS Waist Size [µm]"
Variable/G SrwGsnBeamWaistZ=1; String/G SrwPGsnBeamWaistZ="Vertical RMS Waist Size [µm]"
Variable/G SrwGsnBeamMx=0; String/G SrwPGsnBeamMx="Horizontal Mode Order"
Variable/G SrwGsnBeamMz=0; String/G SrwPGsnBeamMz="Vertical Mode Order"
Variable/G SrwGsnBeamPolar=1; String/G SrwPGsnBeamPolar="Polarization"
Variable/G SrwGsnBeamPhot=1

end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: Setup Source
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeam(SrcName,ElecName,WaistX,WaistZ,mx,mz,Phot,Polar)
String SrcName=SrwGsnBeamName
String ElecName=SrwElecName+SrwElecType
Variable WaistX=SrwGsnBeamWaistX
Variable WaistZ=SrwGsnBeamWaistZ
Variable mx=SrwGsnBeamMx
Variable mz=SrwGsnBeamMz
Variable Phot=SrwGsnBeamPhot
Variable Polar=SrwGsnBeamPolar
prompt SrcName,SrwPGsnBeamName
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","")
prompt WaistX,SrwPGsnBeamWaistX
prompt WaistZ,SrwPGsnBeamWaistZ
prompt mx,SrwPGsnBeamMx
prompt mz,SrwPGsnBeamMz
prompt Phot,SrwPIsotrSrcPhot
prompt Polar,SrwPGsnBeamPolar,popup SrwPOPUPPolar
Silent 1						|	Creating Gaussian Beam structure  ...
PauseUpdate

if(strlen(SrcName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(SrcName)>28)
	Abort SrwPAlertTooLongName
endif
if(cmpstr(ElecName,"_none_")==0)
	Abort SrwPAlertElecBeamNeeded
endif
if((WaistX<=0.) %| (WaistZ<=0.))
	Abort SrwPAlertGsnBeamWaist
endif
if((mx<0) %| (mz<0))
	Abort SrwPAlertGsnBeamOrder
endif
if(Phot<=0)
	Abort SrwPAlertIsotrSrcPhot
endif
if((Polar<1) %| (Polar>6))
	Abort SrwPAlertGsnBeamPolar
endif

SrwGsnBeamName=SrcName
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
SrwGsnBeamWaistX=WaistX
SrwGsnBeamWaistZ=WaistZ
SrwGsnBeamMx=mx
SrwGsnBeamMz=mz
SrwGsnBeamPhot=Phot
SrwGsnBeamPolar=Polar

SrcName+=SrwGsnBeamType

Make/T/O/N=10 $SrcName
$SrcName[0]=ElecName
$SrcName[1]=num2str(WaistX*1.e-06)
$SrcName[2]=num2str(WaistZ*1.e-06)
$SrcName[3]=num2str(mx)
$SrcName[4]=num2str(mz)
$SrcName[5]=num2str(Phot)
$SrcName[6]=num2str(Polar)

end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: create source
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeamCreate(SrcName,SpecFlux,WaistX,WaistZ,mx,mz,Polar)
String SrcName=srwUtiGetValS("SrwGsnBeamName","GausBeam","")
Variable SpecFlux=srwUtiGetValN("SpecFlux", 1., "SrwGsnBeamCreate")
Variable WaistX=srwUtiGetValN("SrwGsnBeamWaistX", 1., "")
Variable WaistZ=srwUtiGetValN("SrwGsnBeamWaistZ", 1., "")
Variable mx=srwUtiGetValN("SrwGsnBeamMx", 0, "")
Variable mz=srwUtiGetValN("SrwGsnBeamMz", 0, "")
Variable Polar=srwUtiGetValN("SrwGsnBeamPolar", 1, "")
prompt SrcName,SrwPGsnBeamName
prompt SpecFlux,"Spectral Flux [x 1.E+15 Phot/s/0.1%bw]"
prompt WaistX,"Horizontal RMS Waist Size [µm]"
prompt WaistZ,"Vertical RMS Waist Size [µm]"
prompt mx,SrwPGsnBeamMx
prompt mz,SrwPGsnBeamMz
prompt Polar,SrwPGsnBeamPolar,popup SrwPOPUPPolar
Silent 1						|	Creating Gaussian Beam structure  ...
PauseUpdate

if(strlen(SrcName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(SrcName)>28)
	Abort SrwPAlertTooLongName
endif
if(SpecFlux<=0)
	Abort SrwPAlertIsotrSrcPhot
endif
if((WaistX<=0.) %| (WaistZ<=0.))
	Abort SrwPAlertGsnBeamWaist
endif
if((mx<0) %| (mz<0))
	Abort SrwPAlertGsnBeamOrder
endif
if((Polar<1) %| (Polar>6))
	Abort SrwPAlertGsnBeamPolar
endif

srwUtiSetValS("SrwGsnBeamName", SrcName, "")
srwUtiSetValN("SpecFlux", SpecFlux, "SrwGsnBeamCreate")
srwUtiSetValN("SrwGsnBeamWaistX", WaistX, "")
srwUtiSetValN("SrwGsnBeamWaistZ", WaistZ, "")
srwUtiSetValN("SrwGsnBeamMx", mx, "")
srwUtiSetValN("SrwGsnBeamMz", mz, "")
srwUtiSetValN("SrwGsnBeamPolar", Polar, "")

string ElecName = SrcName+"Elec"

string TotElecName = ElecName + SrwElecType
if(exists(TotElecName) == 0)
	SrwElecFilament(ElecName,1,-1,0,0,0,0,0)
endif

SrcName+=SrwGsnBeamType

Make/T/O/N=20 $SrcName
$SrcName[0]=TotElecName
$SrcName[1]=num2str(WaistX*1.e-06)
$SrcName[2]=num2str(WaistZ*1.e-06)
$SrcName[3]=num2str(mx)
$SrcName[4]=num2str(mz)
$SrcName[5]=num2str(SpecFlux)
$SrcName[6]=num2str(Polar)
end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: create source
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeamPulsedLaser(SrcName,SpecFlux,WaistX,WaistZ,mx,mz,Wavelength,Polar,SigT,RepRate)
string SrcName=srwUtiGetValS("SrwGsnBeamName","GausBeam","")
variable SpecFlux=srwUtiGetValN("SpecFlux", 1., "SrwGsnBeamCreate")
variable WaistX=srwUtiGetValN("SrwGsnBeamWaistX", 100., "")
variable WaistZ=srwUtiGetValN("SrwGsnBeamWaistZ", 100., "")
variable mx=srwUtiGetValN("SrwGsnBeamMx", 0, "")
variable mz=srwUtiGetValN("SrwGsnBeamMz", 0, "")
variable Polar=srwUtiGetValN("SrwGsnBeamPolar", 1, "")
variable Wavelength=srwUtiGetValN("Wavelength", 800, "SrwGsnBeamPulsedLaser")
variable SigT=srwUtiGetValN("SigT", 50, "SrwGsnBeamPulsedLaser")
variable RepRate=srwUtiGetValN("RepRate", 1, "SrwGsnBeamPulsedLaser")
prompt SrcName,SrwPGsnBeamName
prompt SpecFlux,"Energy in Pulse [mJ]"
prompt WaistX,"Horiz. RMS Size (Intens.) at Waist [µm]"
prompt WaistZ,"Vert. RMS Size (Intens.) at Waist [µm]"
prompt mx,SrwPGsnBeamMx
prompt mz,SrwPGsnBeamMz
prompt Wavelength,"Average Wavelength [nm]"
prompt Polar,SrwPGsnBeamPolar,popup SrwPOPUPPolar
prompt SigT,"RMS Pulse Duration [fs]"
prompt RepRate,"Repetition Rate [kHz]"
Silent 1						|	Creating Gaussian Beam structure  ...
PauseUpdate
srwUtiSetValS("SrwGsnBeamName", SrcName, "")
srwUtiSetValN("SpecFlux", SpecFlux, "SrwGsnBeamCreate")
srwUtiSetValN("SrwGsnBeamWaistX", WaistX, "")
srwUtiSetValN("SrwGsnBeamWaistZ", WaistZ, "")
srwUtiSetValN("SrwGsnBeamMx", mx, "")
srwUtiSetValN("SrwGsnBeamMz", mz, "")
srwUtiSetValN("SrwGsnBeamPolar", Polar, "")
srwUtiSetValN("Wavelength", Wavelength, "SrwGsnBeamPulsedLaser")
srwUtiSetValN("SigT", SigT, "SrwGsnBeamPulsedLaser")
srwUtiSetValN("RepRate", RepRate, "SrwGsnBeamPulsedLaser")

if(strlen(SrcName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(SrcName)>28)
	Abort SrwPAlertTooLongName
endif
if(SpecFlux<=0)
	Abort SrwPAlertIsotrSrcPhot
endif
if((WaistX<=0.) %| (WaistZ<=0.))
	Abort SrwPAlertGsnBeamWaist
endif
if((mx<0) %| (mz<0))
	Abort SrwPAlertGsnBeamOrder
endif
if((Polar<1) %| (Polar>6))
	Abort SrwPAlertGsnBeamPolar
endif

string ElecName = SrcName+"Elec"

string TotElecName = ElecName + SrwElecType
if(exists(TotElecName) == 0)
	SrwElecFilament(ElecName,1,-1,0,0,0,0,0)
endif

SrcName+=SrwGsnBeamType

Make/T/O/N=20 $SrcName
$SrcName[0]=TotElecName
$SrcName[1]=num2str(WaistX*1.e-06)
$SrcName[2]=num2str(WaistZ*1.e-06)
$SrcName[3]=num2str(mx)
$SrcName[4]=num2str(mz)
$SrcName[5]=num2str(0) //set 0 for spectral flux
$SrcName[6]=num2str(Polar)

$SrcName[14]=num2str(SigT*1.e-15) //[s]
$SrcName[16]=num2str(RepRate*1000) //[Hz]
$SrcName[17]=num2str(SpecFlux*0.001) //[J]

variable avgPhotEn_eV = 1.239842*1000/Wavelength
string strAvgPhotEn
sprintf strAvgPhotEn, "%16.9g", avgPhotEn_eV
$SrcName[18]=strAvgPhotEn //[eV]
//$SrcName[18]=num2str(1.239842*1000/Wavelength) //[eV]

end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: set positions and angles
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeamMom1(SrcName,s0,x0,z0,xp0,zp0)
String SrcName=srwUtiGetValS("SrwGsnBeamName","GausBeam","")+SrwGsnBeamType
Variable s0=srwUtiGetValN("s0", 0, "SrwGsnBeamMom1")
Variable x0=srwUtiGetValN("x0", 0, "SrwGsnBeamMom1")
Variable xp0=srwUtiGetValN("xp0", 0, "SrwGsnBeamMom1")
Variable z0=srwUtiGetValN("z0", 0, "SrwGsnBeamMom1")
Variable zp0=srwUtiGetValN("zp0", 0, "SrwGsnBeamMom1")
prompt SrcName,SrwPGsnBeamName2,popup Wavelist("*"+SrwGsnBeamType,";","")
prompt s0,"Initial Longitudinal Position [m]"
prompt x0,"Average Horizontal Position [mm]"
prompt xp0,"Average Horizontal Angle [mrad]"
prompt z0,"Average Vertical Position [mm]"
prompt zp0,"Average Vertical Angle [mrad]"
Silent 1						|	Creating Gaussian Beam structure  ...
PauseUpdate
srwUtiSetValS("SrwGsnBeamName", SrcName[0,strlen(SrcName)-strlen(SrwGsnBeamType)-1], "")
srwUtiSetValN("s0", s0, "SrwGsnBeamMom1")
srwUtiSetValN("x0", x0, "SrwGsnBeamMom1")
srwUtiSetValN("xp0", xp0, "SrwGsnBeamMom1")
srwUtiSetValN("z0", z0, "SrwGsnBeamMom1")
srwUtiSetValN("zp0", zp0, "SrwGsnBeamMom1")

if(strlen(SrcName)==0)
	Abort SrwPAlertBadName
endif

$SrcName[8]=num2str(s0)
$SrcName[9]=num2str(x0*1.e-03)
$SrcName[10]=num2str(xp0*1.e-03)
$SrcName[11]=num2str(z0*1.e-03)
$SrcName[12]=num2str(zp0*1.e-03)

string ElecName = SrwGsnBeamName+"Elec"
SrwElecFilament(ElecName,1,-1,s0,x0,z0,xp0,zp0)
ElecName += SrwElecType
$SrcName[0]=ElecName
end

//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: set time parameters
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwGsnBeamTime(SrcName,sigt,type)
String SrcName=srwUtiGetValS("SrwGsnBeamName","GausBeam","")+SrwGsnBeamType
Variable sigt=srwUtiGetValN("sigt", 100, "SrwGsnBeamTime")
Variable type=srwUtiGetValN("type", 1, "SrwGsnBeamTime")
prompt SrcName,SrwPGsnBeamName2,popup Wavelist("*"+SrwGsnBeamType,";","")
prompt sigt,"RMS Pulse Duration [fs]"
prompt type,"Pulse Form-Factor",popup "Gaussian;Rectangular;Triangular"
Silent 1						|	Creating Gaussian Beam structure  ...
PauseUpdate
srwUtiSetValS("SrwGsnBeamName", SrcName[0,strlen(SrcName)-strlen(SrwGsnBeamType)-1], "")
srwUtiSetValN("sigt", sigt, "SrwGsnBeamTime")
srwUtiSetValN("type", type, "SrwGsnBeamTime")

if(strlen(SrcName)==0)
	Abort SrwPAlertBadName
endif

$SrcName[14]=num2str(sigt*1.e-15)
$SrcName[15]=num2str(type)
end


//+++++++++++++++++++++++++++++++++++++++
//
// Gaussian Beam: Create Wavefront
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwWfrGsnBeamCreate(RadName,SrcName,ObsName,ObsNxNzForProp,ObsNxNzSamplFact)
String RadName=srwUtiTruncString(SrwGsnBeamName+SrwSmpName, 27)
String SrcName=SrwGsnBeamName+SrwGsnBeamType
String ObsName=SrwSmpName+SrwSmpType
Variable ObsNxNzForProp=SrwSmpNxNzForProp
Variable ObsNxNzSamplFact=SrwSmpNxNzSamplFact
prompt RadName,SrwPRadName
prompt SrcName,SrwPGsnBeamName2,popup Wavelist("*"+SrwGsnBeamType,";","")
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","")
prompt ObsNxNzForProp,SrwPSmpNxNzForProp,popup "No;Yes"
prompt ObsNxNzSamplFact,SrwPSmpNxNzSamplFact
Silent 1						|	Creating the Wavefront  ...
PauseUpdate

if(strlen(RadName)==0)
	Abort SrwPAlertBadName
endif
if(strlen(RadName)>28)
	Abort SrwPAlertTooLongName
endif

if(cmpstr(SrcName,"_none_")==0)
	Abort SrwPAlertGsnBeamNeeded
endif
if(cmpstr(ObsName,"_none_")==0)
	Abort SrwPAlertRadSmplNeeded
endif

if(ObsNxNzForProp<1)
	ObsNxNzForProp=1
endif
if(ObsNxNzForProp>2)
	ObsNxNzForProp=2
endif
if(ObsNxNzSamplFact<=0.)
	Abort SrwPAlertBadNxNzSamplFact
endif

SrwGsnBeamName=SrcName[0,strlen(SrcName)-strlen(SrwGsnBeamType)-1]
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1]
SrwSmpNxNzForProp=ObsNxNzForProp
SrwSmpNxNzSamplFact=ObsNxNzSamplFact

if(strlen(RadName)==0)
	RadName=SrwGsnBeamName+SrwSmpName							
endif
SrwRadName=RadName

variable numObsPhotEnVals = srwGetSmpPhotEnNp(ObsName)
variable numObsTimeVals = srwGetSmpTimeNp(ObsName)
if((numObsPhotEnVals > 1) %& (numObsTimeVals > 1))
	abort "Can't compute electric field for many photon energy values and many time moments."
endif

Make/D/O/N=5 AuxObsTreat; AuxObsTreat[0]=ObsNxNzForProp-1; AuxObsTreat[1]=ObsNxNzSamplFact

SrwWfrPrep($SrcName[0],ObsName,RadName,0)
RadName += SrwRadType
SrwRadGenTotName=RadName

srWfrGsnBeam($SrcName,$ObsName,AuxObsTreat,$RadName)
KillWaves/Z AuxObsTreat

end
