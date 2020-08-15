
//+++++++++++++++++++++++++++++++++++++++
//
//Compute Stokes parameters of Wiggler radiation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoWigglerCreate(RadName, ElecName, MagName, ObsName, PrecPar, TreatInterf)
String RadName=srwUtiTruncString(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwMagGenTotName;
String ObsName=SrwSmpName+SrwSmpType;
Variable PrecPar=SrwStoWigCompPrec;
Variable TreatInterf=srwUtiGetValN("TreatInterf", 1, "SrwStoWigCreate")
prompt RadName,SrwPStoName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "");
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")+Wavelist("*"+SrwUndType ,";", "");
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "");
prompt PrecPar,SrwPStoWigCompPrec;
prompt TreatInterf,"Include Interference Effects?",popup "No;Yes";
Silent 1						|	Computing the Radiation ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0;
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0;
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0;
endif
if(ElecWavePresent==0)
	SrwElecFilament();
	SrwElecThick();
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwStoWigCreate(); 
		Return;
	endif
endif
if(MagWavePresent==0)
	DoAlert 0, SrwPAlertMagFieldNeeded;
	Abort;
endif
if(ObsWavePresent==0)
	SrwStartMacrosAfterRadSmp = 5;
	SrwStartMacrosAfterRadSmp2 *= -1;
	SrwRadSamplDialog(SrwSmpType);
	Return;
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1];
SrwStoWigCompPrec=PrecPar;

srwUtiSetValN("TreatInterf", TreatInterf, "SrwStoWigCreate");

string BufMagName = MagName[0,strlen(MagName)-strlen(SrwFieldType)-1];
string MagType = MagName[strlen(MagName)-strlen(SrwFieldType),strlen(MagName)-1];
if(cmpstr(MagType,SrwFieldType)==0)
	SrwMagName = BufMagName;
endif
if(cmpstr(MagType,SrwUndType)==0)
	SrwUndName = BufMagName;
endif
SrwMagGenTotName = MagName;

if(strlen(RadName)==0)
	RadName=SrwElecName+BufMagName+SrwSmpName;
endif
SrwStoName=RadName
SrwRadGenTotName=RadName

// Preparing data for C function
Make/D/O/N=6 waveprec;
waveprec[0]=PrecPar;
waveprec[1]=TreatInterf;

Variable FluxUn=2; // To set Flux per untt surface units
SrwStoPrep(ElecName,MagName,Obsname,RadName,FluxUn)
RadName += SrwStoType;
SrwRadGenTotName=RadName

srStokesWig($ElecName, $MagName, $ObsName, waveprec, $RadName);
KillWaves/Z  waveprec;
end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Stokes parameters of Wiggler radiation
//Obsolete version
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoWigCreate(RadName, ElecName, MagName, ObsName, PrecPar)
String RadName=srwUtiTruncString(SrwElecName+SrwMagGenTotName[0,strlen(SrwMagGenTotName)-strlen(SrwFieldType)-1]+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwMagGenTotName;
String ObsName=SrwSmpName+SrwSmpType;
Variable PrecPar=SrwStoWigCompPrec;
prompt RadName,SrwPStoName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType ,";", "");
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwFieldType ,";", "")+Wavelist("*"+SrwUndType ,";", "");
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType ,";", "");
prompt PrecPar,SrwPStoWigCompPrec;
Silent 1						|	Computing the Radiation ...
PauseUpdate

SrwStoWigglerCreate(RadName, ElecName, MagName, ObsName, PrecPar, 2)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute Stokes parameters of dipole magnet radiation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwStoConstCreate(RadName, ElecName, MagName, ObsName)
String RadName=srwUtiTruncString(SrwElecName+SrwMagConstName+SrwSmpName, 27)
String ElecName=SrwElecName+SrwElecType;
String MagName=SrwMagConstName+SrwMagConstType;
String ObsName=SrwSmpName+SrwSmpType;
prompt RadName,SrwPStoName;
prompt ElecName,SrwPElecName1,popup Wavelist("*"+SrwElecType,";","");
prompt MagName,SrwPMagName2,popup Wavelist("*"+SrwMagConstType,";","");
prompt ObsName,SrwPSmpName2,popup Wavelist("*"+SrwSmpType,";","");
Silent 1						|	Computing the Radiation ...
PauseUpdate

variable ElecWavePresent = 1, MagWavePresent = 1, ObsWavePresent = 1
if(cmpstr(ElecName,"_none_")==0)
	ElecWavePresent = 0;
endif
if(cmpstr(MagName,"_none_")==0)
	MagWavePresent = 0;
endif
if(cmpstr(ObsName,"_none_")==0)
	ObsWavePresent = 0;
endif
if(ElecWavePresent==0)
	SrwElecFilament();
	SrwElecThick();
	if((MagWavePresent == 1) %& (ObsWavePresent == 1))
		SrwStoConstCreate(); 
		Return;
	endif
endif
if(MagWavePresent==0)
	SrwMagConstCreate();
	if(ObsWavePresent == 1)
		SrwStoConstCreate(); 
		Return;
	endif
endif
if(ObsWavePresent==0)
	DoAlert 0, SrwPAlertRadSmplNeeded;
	Abort;
endif

Variable PrecPar=1;

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];
SrwMagConstName=MagName[0,strlen(MagName)-strlen(SrwMagConstType)-1];
SrwSmpName=ObsName[0,strlen(ObsName)-strlen(SrwSmpType)-1];
SrwStoConstCompPrec=PrecPar;

if(strlen(RadName)==0)
	RadName=SrwElecName+SrwMagConstName+SrwSmpName
endif
SrwStoName=RadName

// Preparing data for C function
Make/D/O/N=6 waveprec;
waveprec[0]=PrecPar;

Variable FluxUn=2; // To set Flux per unit surface units
SrwStoPrep(ElecName,MagName,Obsname,RadName,FluxUn);
RadName += SrwStoType;
SrwRadGenTotName=RadName

srStokesConst($ElecName, $MagName, $ObsName, waveprec, $RadName);

KillWaves/Z  waveprec;

end

//+++++++++++++++++++++++++++++++++++++++
//
//Compute maximal spectral flux of wiggler radiation in the on-axis direction
//taking into account constraints on Total Power and Wiggler Length
//
//Call examples:
//make/O/N=(101, 101) wFluxPerUnitAngleWhiteSol;SetScale/I x 10,300,"", wFluxPerUnitAngleWhiteSol;SetScale/I y 0.1,10,"", wFluxPerUnitAngleWhiteSol; wFluxPerUnitAngleWhiteSol=0
//make/O/N=(101, 101) wFluxPerUnitAngleSol;SetScale/I x 10,300,"", wFluxPerUnitAngleSol;SetScale/I y 0.1,10,"", wFluxPerUnitAngleSol; wFluxPerUnitAngleSol=0
//wFluxPerUnitAngleWhiteSol = SrwWigSpecFluxAtPowerConstr("SOL_SS_ebm", 25, 1.8, 5, 20, 100, y, x)
//wFluxPerUnitAngleDia = SrwWigSpecFluxAtPowerConstr(3., 0.3, 20, 2, 50, 50, y, x)
//+++++++++++++++++++++++++++++++++++++++
function SrwWigSpecFluxAtPowerConstr(nmElecBeam, pow_kW, Lmax, f_type, ph_en1_keV, ph_en2_keV, b_T, per_mm)
string nmElecBeam
variable pow_kW, Lmax, f_type, ph_en1_keV, ph_en2_keV, b_T, per_mm
//f_type: 
//	1- "Phot/s/.1%/mr"
//	2- "Phot/s/.1%/mr2"
//	3- "Phot/s/.1%/mr2/mm2"
//	4- "W/mr (integ. vs photon energy)"
//	5- "W/mr2 (integ. vs photon energy)"
//	6- "W/mr2/mm2 (integ. vs photon energy)"

//variable EmitX = 2.7 //[nm] //3.73 // [nm]
//variable Lmax = 2 //[m]

variable en_GeV = srwGetElecBeamEnergy(nmElecBeam)
variable cur_A = srwGetElecBeamCurrent(nmElecBeam)

variable Np = pow_kW/((0.633e-03)*en_GeV*en_GeV*cur_A*b_T*b_T*per_mm)
if(Np < 1) 
	return 0
endif
if(0.001*Np*per_mm > Lmax) 
	Np = Lmax/(0.001*per_mm)
endif

variable wigL_m = 0.001*Np*per_mm
variable wigK = 0.0934*per_mm*b_T

variable/G SrwAllowPrintingExtraInfo
string/G SrwGlStr, SrwAfStr, SrwBrStr

string ExtStr = ""
if((f_type==1) %| (f_type==4))
	ExtStr=SrwGlStr
endif
if((f_type==2) %| (f_type==5))
	ExtStr=SrwAfStr
endif
if((f_type==3) %| (f_type==6))
	ExtStr=SrwBrStr
endif
string BrilWaveName = "ElecW." + ExtStr

variable PrevSrwAllowPrintingExtraInfo = SrwAllowPrintingExtraInfo
SrwAllowPrintingExtraInfo = 0

string StrToExe = ""
//sprintf StrToExe, "SrwElecFilament(\"Elec\",%g,%g,0,0,0,0,0);SrwElecThick(\"Elec_ebm\",0.001016,%g,0.037,17.78,1.75,0,0,0.28,0)", en_GeV, cur_A, EmitX //for diamond

execute StrToExe
sprintf StrToExe, "SrwMagPerCreate2D(\"W\",%g,%g,0,%g,0,1,0,0)", per_mm, wigK, wigL_m
execute StrToExe
sprintf StrToExe, "SrwBrilWig(\"ElecW\",\"Elec_ebm\",\"W_map\",1,%g,%g,%g,1)", ph_en1_keV, ph_en2_keV, f_type
	//print StrToExe
execute StrToExe

SrwAllowPrintingExtraInfo = PrevSrwAllowPrintingExtraInfo

wave wBrilWaveName = $BrilWaveName
return wBrilWaveName[0]
end
