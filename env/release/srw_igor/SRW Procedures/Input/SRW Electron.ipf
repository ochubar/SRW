
//+++++++++++++++++++++++++++++++++++++++
//
//Modify e-beam structure. Thick e-beam parameters. Based on Twiss.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecThick(ElecName,sige,emx,emz,betax,betaz,alphax,alphaz,etax,etaxpr)
string ElecName=SrwElecName+SrwElecType;
variable emx=SrwElecEmx;
variable betax=SrwElecBetax;
variable alphax=SrwElecAlphax;
variable emz=SrwElecEmz;
variable betaz=SrwElecBetaz;
variable alphaz=SrwElecAlphaz;
variable sige=SrwElecSige;
variable etax=SrwElecEtax;
variable etaxpr=SrwElecEtaxPr;
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","");
prompt sige,SrwPElecSige;
prompt emx,SrwPElecEmx;
prompt emz,SrwPElecEmz;
prompt betax,SrwPElecBetax;
prompt betaz,SrwPElecBetaz;
prompt alphax,SrwPElecAlphax;
prompt alphaz,SrwPElecAlphaz;
prompt etax,SrwPElecEtax;
prompt etaxpr,SrwPElecEtaxPr;
Silent 1						|	 ...
PauseUpdate;

if(cmpstr(ElecName,"_none_")==0)
	SrwElecFilament();
	SrwElecThick();
	Return;
endif

// Validation of parameters
if(sige < 0.)
	Abort SrwPAlertElecEnSpr
endif
if(emx < 0.)
	Abort SrwPAlertElecEmittance
endif
if(emx < 0.)
	Abort SrwPAlertElecEmittance
endif

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];
SrwElecEmx=emx;
SrwElecBetax=betax;
SrwElecAlphax=alphax;
SrwElecEmz=emz;
SrwElecBetaz=betaz;
SrwElecAlphaz=alphaz;
SrwElecSige=sige;
SrwElecEtax=etax;
SrwElecEtaxPr=etaxpr;
SrwElecEtaz=0;
SrwElecEtazPr=0;

$ElecName[7]=emx;
$ElecName[8]=betax;
$ElecName[9]=alphax;
$ElecName[10]=emz; 
$ElecName[11]=betaz;
$ElecName[12]=alphaz;
$ElecName[13]=sige;
$ElecName[14]=etax;
$ElecName[15]=etaxpr;
$ElecName[16]=0;
$ElecName[17]=0;
$ElecName[18]=0; // Set up from Twiss

variable emx_m = emx*1.e-9, emz_m = emz*1.e-9, sigeE2 = sige*sige;
 
$ElecName[20]=emx_m*betax + sigeE2*etax*etax;   // <(x-<x>)^2>
$ElecName[21]=-emx_m*alphax + sigeE2*etax*etaxpr;   // <(x-<x>)(x'-<x'>)>
$ElecName[22]=emx_m*(1 + alphax*alphax)/betax + sigeE2*etaxpr*etaxpr;   // <(x'-<x'>)^2>

$ElecName[23]=emz_m*betaz;   // <(z-<z>)^2>
$ElecName[24]=-emz_m*alphaz;   // <(z-<z>)(z'-<z'>)>
$ElecName[25]=emz_m*(1 + alphaz*alphaz)/betaz;   // <(z'-<z'>)^2>

$ElecName[26]=0.;   // <(x-<x>)(z-<z>)>
$ElecName[27]=0.;   // <(x'-<x'>)(z-<z>)>
$ElecName[28]=0.;   // <(x-<x>)(z'-<z'>)>
$ElecName[29]=0.;   // <(x'-<x'>)(z'-<z'>)>

variable s0 = ($ElecName[6]);
variable sigx = sqrt($ElecName[20])*1.e+3;
variable sigxpr = sqrt($ElecName[22])*1.e+3;
variable sxWaist = s0 - ($ElecName[21])/($ElecName[22]);
variable sigz = sqrt($ElecName[23])*1.e+3;
variable sigzpr = sqrt($ElecName[25])*1.e+3;
variable szWaist = s0 - ($ElecName[24])/($ElecName[25]);

if(SrwAllowPrintingExtraInfo==1)
	sRElecPrintSizes(s0, sigx, sigxpr, sxWaist, sigz, sigzpr, szWaist);
endif

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
function sRElecPrintSizes(s0, sigx, sigxpr, swx, sigz, sigzpr, swz)
variable s0, sigx, sigxpr, swx, sigz, sigzpr, swz;
string s0Str, sigxStr, sigxprStr, swxStr, sigzStr, sigzprStr, swzStr;
sprintf s0Str "%7.3g" s0;
sprintf sigxStr "%9.4g" sigx;
sprintf sigzStr "%9.4g" sigz;
sprintf sigxprStr "%9.4g" sigxpr;
sprintf sigzprStr "%9.4g" sigzpr;
sprintf swxStr "%9.4g" swx;
sprintf swzStr "%9.4g" swz;
print ""
print "Electron Beam Dimensions at Position :", s0Str,"m";
print "\t\t\t\t\t   Horizontal\t\t   Vertical";
print "RMS Size [mm]\t\t", sigxStr, "\t\t", sigzStr;
print "RMS Diverg. [mr]\t", sigxprStr, "\t\t", sigzprStr;
print "Pos. of Waist [m]\t", swxStr, "\t\t", swzStr;
end

//+++++++++++++++++++++++++++++++++++++++
//
//Modify e-beam structure. Thick e-beam parameters. Based on Moments.
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecThickMom(ElecName,sige,sigx,sigz,sigxp,sigzp,mxxp,mzzp)
string ElecName=SrwElecName+SrwElecType;
variable sige=SrwElecSige;
variable sigx=SrwElecSigx;
variable sigz=SrwElecSigz;
variable sigxp=SrwElecSigxp;
variable sigzp=SrwElecSigzp;
variable mxxp=SrwElecMxxp;
variable mzzp=SrwElecMzzp;
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","");
prompt sige,SrwPElecSige;
prompt sigx,SrwPElecSigx2;
prompt sigz,SrwPElecSigz2;
prompt sigxp,SrwPElecSigxp;
prompt sigzp,SrwPElecSigzp;
prompt mxxp,SrwPElecMxxp;
prompt mzzp,SrwPElecMzzp;
Silent 1						|	 ...
PauseUpdate;

//if(cmpstr(ElecName,"_none_")==0)
//	SrwElecFilament();
//	SrwElecThick();
//	Return;
//endif

// Validation of parameters
//...

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1];
SrwElecSige=sige;
SrwElecSigx=sigx;
SrwElecSigz=sigz;
SrwElecSigxp=sigxp;
SrwElecSigzp=sigzp;
SrwElecMxxp=mxxp;
SrwElecMzzp=mzzp;

variable sigx_m = sigx*1.e-6, sigz_m = sigz*1.e-6;
variable sigxp_r = sigxp*1.e-6, sigzp_r = sigzp*1.e-6;
variable mxxp_m = mxxp*1.e-9, mzzp_m = mzzp*1.e-9;
variable sigx_mE2 = sigx_m*sigx_m, sigz_mE2 = sigz_m*sigz_m;
variable sigxp_rE2 = sigxp_r*sigxp_r, sigzp_rE2 = sigzp_r*sigzp_r;
variable emx_m = sqrt(sigx_mE2*sigxp_rE2 - mxxp_m*mxxp_m);
variable emz_m = sqrt(sigz_mE2*sigzp_rE2 - mzzp_m*mzzp_m);

	//print sigx_mE2*sigxp_rE2, mxxp_m*mxxp_m
	//print sigz_mE2

$ElecName[7]=(1.e+09)*emx_m; // emx in nm
$ElecName[8]=sigx_mE2/emx_m; // betax in m
$ElecName[9]=-mxxp_m/emx_m; // alphax in r
$ElecName[10]=(1.e+09)*emz_m; // emz in nm
$ElecName[11]=sigz_mE2/emz_m; // betaz in m
$ElecName[12]=-mzzp_m/emz_m; // alphaz in r
$ElecName[13]=sige;
$ElecName[14]=0;
$ElecName[15]=0;
$ElecName[16]=0;
$ElecName[17]=0;
$ElecName[18]=1; // Set up from Moments

$ElecName[20]=sigx_mE2;   // <(x-<x>)^2>
$ElecName[21]=mxxp_m;   // <(x-<x>)(x'-<x'>)>
$ElecName[22]=sigxp_rE2;   // <(x'-<x'>)^2>
$ElecName[23]=sigz_mE2;   // <(z-<z>)^2>
$ElecName[24]=mzzp_m;   // <(z-<z>)(z'-<z'>)>
$ElecName[25]=sigzp_rE2;   // <(z'-<z'>)^2>
$ElecName[26]=0.;   // <(x-<x>)(z-<z>)>
$ElecName[27]=0.;   // <(x'-<x'>)(z-<z>)>
$ElecName[28]=0.;   // <(x-<x>)(z'-<z'>)>
$ElecName[29]=0.;   // <(x'-<x'>)(z'-<z'>)>

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Creates the Electron beam structure and fills it with electron energy and current
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecFilament(ElecName,ElecEn,ElecCur,s0,xx,xp,zz,zp)
String  ElecName=SrwElecName
Variable ElecEn=SrwElecEn
Variable ElecCur=SrwElecCur
Variable s0=SrwElecs0
Variable xx=SrwElecxx
Variable xp=SrwElecxp
Variable zz=SrwEleczz
Variable zp=SrwEleczp
prompt ElecName,SrwPElecName
prompt ElecEn,SrwPElecEn
prompt ElecCur,SrwPElecCur
prompt s0,SrwPElecs0
prompt xx,SrwPElecxx
prompt xp,SrwPElecxp
prompt zz,SrwPEleczz
prompt zp,SrwPEleczp
Silent 1						|	Creating Electron Beam structure ...
PauseUpdate

// Validation of parameters
if(ElecEn <= 0.)
	Abort SrwPAlertElecEnergy
endif

SrwElecEn=ElecEn
SrwElecCur=ElecCur
SrwElecName=ElecName
SrwElecxx=xx
SrwElecxp=xp
SrwEleczz=zz
SrwEleczp=zp
SrwElecs0=s0

ElecName+=SrwElecType
//Make/N=30/D/O $ElecName
//Make/N=44/D/O $ElecName
Make/N=46/D/O $ElecName
$ElecName[0]=ElecEn
$ElecName[1]=ElecCur
$ElecName[2]=xx*0.001
$ElecName[3]=xp*0.001
$ElecName[4]=zz*0.001
$ElecName[5]=zp*0.001
$ElecName[6]=s0

//$ElecName[7]=SrwElecEmx;
//$ElecName[8]=SrwElecBetax;
//$ElecName[9]=SrwElecAlphax;
//$ElecName[10]=SrwElecEmz; 
//$ElecName[11]=SrwElecBetaz;
//$ElecName[12]=SrwElecAlphaz;
//$ElecName[13]=SrwElecSige;
//$ElecName[14]=SrwElecEtax;
//$ElecName[15]=SrwElecEtaxPr;

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
End

//+++++++++++++++++++++++++++++++++++++++
//
//Sets number of electrons in bunch and longitudinal bunch parameters
//Obsolete version
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecLong(ElecName,Neb,Sc,SigS)
string ElecName=SrwElecName+SrwElecType
variable Neb=srwUtiGetValN("SrwElecNeb", 1, "")
variable Sc=srwUtiGetValN("SrwElecSc", 0, "")
variable SigS=srwUtiGetValN("SrwElecSigS", 1, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
prompt Neb, "Number of Electrons in Bunch [x 1.E+09]"
prompt Sc, "Longitudinal Position of the Bunch Center [mm]"
prompt SigS, "RMS Bunch Length [mm]"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
srwUtiSetValN("SrwElecNeb", Neb, "")
srwUtiSetValN("SrwElecSc", Sc, "")
srwUtiSetValN("SrwElecSigS", SigS, "")

if(dimsize($ElecName, 0) < 46)
	Redimension/N=46 $ElecName
endif
//$ElecName[18]=(1.e+09)*Neb // Neb
$ElecName[19]=(1.e-03)*Sc // Sc in m
$ElecName[33]=(1.e-06)*SigS*SigS // SigS in m
$ElecName[43]=(1.e+09)*Neb // Neb

variable CurPeak=$ElecName[45]
variable ElCharge = 1.602189246e-19 // [Coulomb]
variable LightSpeed = 2.99792458e+08 // [m/s]
if(CurPeak > 0)
	//update this field only if the peak current value was set previously
	$ElecName[45] = Neb*(1.e+09)*ElCharge*LightSpeed/(sqrt(2*Pi)*SigS*(1.e-03))
endif

$ElecName[30] = 2 // Gaussian Transverse distribution of particles (TypeTransv)
$ElecName[31] = 2 // Gaussian Longitudinal distribution of particles (TypeLong)
$ElecName[32] = 1 // ShotNoise

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets number of electrons in bunch and longitudinal bunch parameters
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecLongCur(ElecName,CurPeak,Sc,SigS,ShotNoise)
string ElecName=SrwElecName+SrwElecType
variable CurPeak=srwUtiGetValN("SrwElecCurPeak", 1, "")
variable Sc=srwUtiGetValN("SrwElecSc", 0, "")
variable SigS=srwUtiGetValN("SrwElecSigS", 1, "")
variable ShotNoise=srwUtiGetValN("SrwElecShotNoise", 1, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
prompt CurPeak, "Peak Electron Current [A]"
prompt Sc, "Longitudinal Position of Bunch Center [mm]"
prompt SigS, "RMS Bunch Length [mm]"
prompt ShotNoise,"Shot Noise Factor (i.e. phase fluctuations scaling for SASE calc.)"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
srwUtiSetValN("SrwElecCurPeak", CurPeak, "")
srwUtiSetValN("SrwElecSc", Sc, "")
srwUtiSetValN("SrwElecSigS", SigS, "")
srwUtiSetValN("SrwElecShotNoise", ShotNoise, "")

if(dimsize($ElecName, 0) < 46)
	Redimension/N=46 $ElecName
endif
//$ElecName[18]=(1.e+09)*Neb // Neb
$ElecName[19]=(1.e-03)*Sc // Sc in m

$ElecName[30] = 2 // Gaussian Transverse distribution of particles (TypeTransv)
$ElecName[31] = 2 // Gaussian Longitudinal distribution of particles (TypeLong)

$ElecName[32] = ShotNoise

$ElecName[33]=(1.e-06)*SigS*SigS // SigS in m

variable ElCharge = 1.602189246e-19 // [Coulomb]
variable LightSpeed = 2.99792458e+08 // [m/s]
variable Neb = sqrt(2*Pi)*SigS*(1.e-03)*CurPeak/(ElCharge*LightSpeed)
$ElecName[43]=Neb // number of electrons in a bunch, assuming Gaussian distribution in longitudinal direction
//$ElecName[44] is taken by a sign whether the partcle distribution exists or not
$ElecName[45]=CurPeak // this value to be used in SASE calculations

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets number of electrons in bunch and longitudinal bunch parameters
//
//+++++++++++++++++++++++++++++++++++++++
//proc SrwElecDistrPrep(ElecName, CreateElecDistr, NumMacroPart, NumSlices)
proc SrwElecDistrPrep(ElecName, NumMacroPart, NumSlices)
string ElecName=SrwElecName+SrwElecType
//variable CreateElecDistr=srwUtiGetValN("CreateElecDistr", 1, "SrwElecDistrPrep")
variable NumMacroPart=srwUtiGetValN("SrwPrecSASEnpart", 16384, "")
variable NumSlices=srwUtiGetValN("SrwPrecSASEnslice", 1, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
//prompt CreateElecDistr, "Create Electron Distribution",popup "No;Yes"
prompt NumMacroPart,"Number of Macro-Particles"
prompt NumSlices,"Number of Bunch Slices"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwPrecSASEnpart=NumMacroPart
SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
//srwUtiSetValN("CreateElecDistr", CreateElecDistr, "SrwElecDistrPrep")

string/G SrwElecDistrType = "_ebd" //in line with definition in Init !!!
string ElecNameExtStr = ElecName[strlen(SrwElecName), strlen(ElecName)]
if(cmpstr(ElecNameExtStr, SrwElecDistrType) == 0)
	return //if Electron Distribution file is submitted instead of SRW Electron Beam, nothing should be done here
endif

if(dimsize($ElecName, 0) < 46)
	Redimension/N=46 $ElecName
endif

variable PrevDistrExisted = $ElecName[44]
//$ElecName[44] = CreateElecDistr - 1 //meaning that there is a particle distribution wave with the same core name and different extension
$ElecName[44] = 1
//the above should be set AFTER execution of srSASE

string ElecDistrName = SrwElecName + SrwElecDistrType
variable TotNumPartRecords = NumMacroPart*NumSlices

//if(CreateElecDistr == 2)
	//create numerical wave of particle distribution
if(exists(ElecDistrName) == 0)
	make/O/N=(6, TotNumPartRecords) $ElecDistrName
	$ElecDistrName = 0
else
	redimension/N=(6, TotNumPartRecords) $ElecDistrName
endif
if(PrevDistrExisted == 0)
	$ElecDistrName = 0
endif
	//if non-zero electron distribution will be found, it will be used for SASE computation
//endif

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets second-order statistical moments of particle distribution in e-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecThickMomCrossTransv(ElecName,mxz,mxzp,mxpz,mxpzp)
string ElecName=SrwElecName+SrwElecType
variable mxz=srwUtiGetValN("SrwElecMxz", 0, "")
variable mxzp=srwUtiGetValN("SrwElecMxzp", 0, "")
variable mxpz=srwUtiGetValN("SrwElecMxpz", 0, "")
variable mxpzp=srwUtiGetValN("SrwElecMxpzp", 0, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
prompt mxz,"<(x-<x>)(z-<z>)> [µm**2]"
prompt mxzp,"<(x-<x>)(z'-<z'>)> [nm]"
prompt mxpz,"<(x'-<x'>)(z-<z>)> [nm]"
prompt mxpzp,"<(x'-<x'>)(z'-<z'>)> [µrad**2]"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
srwUtiSetValN("SrwElecMxz", mxz, "")
srwUtiSetValN("SrwElecMxzp", mxzp, "")
srwUtiSetValN("SrwElecMxpz", mxpz, "")
srwUtiSetValN("SrwElecMxpzp", mxpzp, "")

if(dimsize($ElecName, 0) < 30)
	Redimension/N=30 $ElecName
endif
$ElecName[26]=(1.e-12)*mxz // mxz in m^2
$ElecName[27]=(1.e-09)*mxpz // mxpz in m
$ElecName[28]=(1.e-09)*mxzp // mxzp in m
$ElecName[29]=(1.e-12)*mxpzp // mxpzp in rad^2

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets second-order statistical moments of particle distribution in e-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecThickMomCrossEnLong(ElecName,mse,mxe,mxpe,mze,mzpe,mxs,mxps,mzs,mzps)
string ElecName=SrwElecName+SrwElecType
variable mse=srwUtiGetValN("SrwElecMse", 0, "")
variable mxe=srwUtiGetValN("SrwElecMxe", 0, "")
variable mxpe=srwUtiGetValN("SrwElecMxpe", 0, "")
variable mze=srwUtiGetValN("SrwElecMze", 0, "")
variable mzpe=srwUtiGetValN("SrwElecMzpe", 0, "")
variable mxs=srwUtiGetValN("SrwElecMxs", 0, "")
variable mxps=srwUtiGetValN("SrwElecMxps", 0, "")
variable mzs=srwUtiGetValN("SrwElecMzs", 0, "")
variable mzps=srwUtiGetValN("SrwElecMzps", 0, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
prompt mse,"<(s-<s>)(E-<E>)>/<E> [nm]"
prompt mxe,"<(x-<x>)(E-<E>)>/<E> [nm]"
prompt mxpe,"<(x'-<x'>)(E-<E>)>/<E> [µrad**2]"
prompt mze,"<(z-<z>)(E-<E>)>/<E> [nm]"
prompt mzpe,"<(z'-<z'>)(E-<E>)>/<E> [µrad**2]"
prompt mxs,"<(x-<x>)(s-<s>)> [µm**2]"
prompt mxps,"<(x'-<x'>)(s-<s>)> [nm]"
prompt mzs,"<(z-<z>)(s-<s>)> [µm**2]"
prompt mzps,"<(z'-<z'>)(s-<s>)> [nm]"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
srwUtiSetValN("SrwElecMse", mse, "")
srwUtiSetValN("SrwElecMxe", mxe, "")
srwUtiSetValN("SrwElecMxpe", mxpe, "")
srwUtiSetValN("SrwElecMze", mze, "")
srwUtiSetValN("SrwElecMzpe", mzpe, "")
srwUtiSetValN("SrwElecMxs", mxs, "")
srwUtiSetValN("SrwElecMxps", mxps, "")
srwUtiSetValN("SrwElecMzs", mzs, "")
srwUtiSetValN("SrwElecMzps", mzps, "")

if(dimsize($ElecName, 0) < 43)
	Redimension/N=43 $ElecName
endif
$ElecName[34]=(1.e-09)*mse // mse in m
$ElecName[35]=(1.e-09)*mxe // mxe in m
$ElecName[36]=(1.e-12)*mxpe // mxpe in rad^2
$ElecName[37]=(1.e-09)*mze // mze in m
$ElecName[38]=(1.e-12)*mzpe // mzpe in rad^2
$ElecName[39]=(1.e-12)*mxs // mxs in m^2
$ElecName[40]=(1.e-09)*mxps // mxps in m
$ElecName[41]=(1.e-12)*mzs // mzs in m^2
$ElecName[42]=(1.e-09)*mzps // mzpe in m

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets "hard" limits on e-beam phase-space coordinates
//Used by some types of calculations (e.g. "macro-particles" method)
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecThickPhaseCoordRange(ElecName,rE,rX,rXp,rZ,rZp,rS)
string ElecName=SrwElecName+SrwElecType
variable rE=srwUtiGetValN("SrwElec_rE", 1000, "")
variable rX=srwUtiGetValN("SrwElec_rX", 1000, "")
variable rXp=srwUtiGetValN("SrwElec_rXp", 1000, "")
variable rZ=srwUtiGetValN("SrwElec_rZ", 1000, "")
variable rZp=srwUtiGetValN("SrwElec_rZp", 1000, "")
variable rS=srwUtiGetValN("SrwElec_rS", 1000, "")
prompt ElecName,SrwPElecName,popup Wavelist("*"+SrwElecType,";","")
prompt rE,"Energy Range [MeV]"
prompt rX,"Horizontal Coord. Range [mm]"
prompt rXp,"Horizontal Angle Range [mrad]"
prompt rZ,"Vertical Coord. Range [mm]"
prompt rZp,"Vertical Angle Range [mrad]"
prompt rS,"Longitudinal Coord. Range [mm]"
Silent 1						|	Setting-up Electron Beam structure ...
PauseUpdate

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
srwUtiSetValN("SrwElec_rE", rE, "")
srwUtiSetValN("SrwElec_rX", rX, "")
srwUtiSetValN("SrwElec_rXp", rXp, "")
srwUtiSetValN("SrwElec_rZ", rZ, "")
srwUtiSetValN("SrwElec_rZp", rZp, "")
srwUtiSetValN("SrwElec_rS", rS, "")

if(dimsize($ElecName, 0) < 52)
	Redimension/N=52 $ElecName
endif
$ElecName[46]=(1.e-03)*rE // Energy Range [GeV]
$ElecName[47]=(1.e-03)*rX // Horizontal Coord. Range [m]
$ElecName[48]=(1.e-03)*rXp // Horizontal Angle Range [rad]
$ElecName[49]=(1.e-03)*rZ // Vertical Coord. Range [m]
$ElecName[50]=(1.e-03)*rZp // Vertical Angle Range [rad]
$ElecName[51]=(1.e-03)*rS // Longitudinal Coord. Range [mm]

srwUtiSetValS("SrwElecLastTot", SrwElecName+SrwElecType, "")
end

//+++++++++++++++++++++++++++++++++++++++
//
//Create e-beam container and add elements to it
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecCont(cntname,cmp1name,cmp2name,cmp3name,cmp4name)
String cntname=srwUtiGetValS("SrwElecBeamCont", "ElCont", "")
String cmp1name
String cmp2name
String cmp3name
String cmp4name
prompt cntname,"Name of the Electron Beam Container structure"
prompt cmp1name,"1st Beam (/Bunch) to add",popup " ;" + Wavelist("*"+SrwElecType,";","") + Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc",""),";","")
prompt cmp2name,"2nd Beam (/Bunch) to add",popup " ;" + Wavelist("*"+SrwElecType,";","") + Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc",""),";","")
prompt cmp3name,"3rd Beam (/Bunch) to add",popup " ;" + Wavelist("*"+SrwElecType,";","") + Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc",""),";","")
prompt cmp4name,"4th Beam (/Bunch) to add",popup " ;" + Wavelist("*"+SrwElecType,";","") + Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc",""),";","")
Silent 1				|	 ...
PauseUpdate
srwUtiSetValS("SrwElecBeamCont", cntname, "")
srwUtiSetValS("SrwElecLastTot", cntname+srwUtiGetValS("SrwElecContType","_ebc",""), "")

srwUtiSetValS("SrwElecContType","_ebc","")
cntname+=SrwElecContType

if((cmpstr(cntname,cmp1name)==0) %| (cmpstr(cntname,cmp2name)==0) %| (cmpstr(cntname,cmp3name)==0) %| (cmpstr(cntname,cmp4name)==0))
abort  "Can not add container to itself"
endif

Make/T/O/N=1 $cntname
//$cntname[0]=SrwBliContType

if((cmpstr(cmp1name," ") != 0) %& (cmpstr(cmp1name,"") != 0))
	$cntname[0]=cmp1name
	//SrwElecContAdd(cmp1name,cntname)
endif
if((cmpstr(cmp2name," ") != 0) %& (cmpstr(cmp2name,"") != 0))
	SrwElecContAdd(cmp2name,cntname)
endif
if((cmpstr(cmp3name," ") != 0) %& (cmpstr(cmp3name,"") != 0))
	SrwElecContAdd(cmp3name,cntname)
endif
if((cmpstr(cmp4name," ") != 0) %& (cmpstr(cmp4name,"") != 0))
	SrwElecContAdd(cmp4name,cntname)
endif
end

//+++++++++++++++++++++++++++++++++++++++
//
//Append a Beamline Element to a Container
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecContAdd(bmname,groupname)
string bmname=srwUtiGetValS("SrwElecLastTot", "", "")
string groupname=srwUtiGetValS("SrwElecBeamCont", "", "")+srwUtiGetValS("SrwElecContType","_ebc","")
prompt bmname,"Electron Beam structure to be added to the Container",popup Wavelist("*"+SrwElecType ,";", "")
prompt groupname,"Container structure",popup Wavelist("*"+srwUtiGetValS("SrwElecContType","_ebc","") ,";", "")
Silent 1						|	 ...
PauseUpdate
srwUtiSetValS("SrwElecBeamCont", groupname[0,strlen(groupname)-strlen(srwUtiGetValS("SrwElecContType","_ebc",""))-1], "")
srwUtiSetValS("SrwElecLastTot", groupname, "")

if(cmpstr($groupname[0],SrwElecContType)==1)
	abort  "Can't append to a non-container structure"
endif
if(cmpstr(groupname,bmname)==0)
	abort  "Can't append a container to itself"
endif

variable n = DimSize($groupname,0)
redimension/N=(n+1)  $groupname
$groupname[n]=bmname
end

//+++++++++++++++++++++++++++++++++++++++
//
//Creates the Electron beam structure and fills it with electron energy and current
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwElecEnsureZeroEmit(ElecName)
string  ElecName=SrwElecName+SrwElecType
prompt ElecName, "Electron Beam structure"
Silent 1						|	Modifying Electron Beam structure  ...
PauseUpdate
SrwElecName = ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]
$ElecName *= srwUtiStep(6 - p)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Sets parameters of Electron Beam structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamParam(Name, IndParam, ValToSet)
string Name
variable IndParam, ValToSet
SVAR SrwElecType

srwUtiSetValS("SrwElecName", Name[0,strlen(Name)-strlen(SrwElecType)-1], "")
wave wElec = $Name
if(DimSize(wElec, 0) < 7)
	abort "This is not an electron beam structure."
endif
if(IndParam > 30)
	abort "Parameter index value is too big."
endif
wElec[IndParam] = ValToSet
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamEnergy(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 0, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamLongPos(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 6, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamHorPos(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 2, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamHorAng(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 3, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamVertPos(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 4, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamVertAng(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 5, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamRelEnSprRMS(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 13, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamEmitX(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 7, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamBetaX(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 8, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamAlphaX(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 9, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamEmitZ(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 10, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamBetaZ(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 11, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamAlphaZ(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 12, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamDispX(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 14, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamDispDerivX(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 15, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamDispZ(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 16, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
function srwSetElecBeamDispDerivZ(ElecName, ValToSet)
string ElecName
variable ValToSet
srwSetElecBeamParam(ElecName, 17, ValToSet)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Return parameters of Electron Beam structure
//
//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamParam(Name, IndParam)
string Name
variable IndParam
SVAR SrwElecType

srwUtiSetValS("SrwElecName", Name[0,strlen(Name)-strlen(SrwElecType)-1], "")
wave wElec = $Name
if(DimSize(wElec, 0) < 7)
	abort "This is not an electron beam structure."
endif
if(IndParam > 30)
	abort "Parameter index value is too big."
endif
return wElec[IndParam]
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamEnergy(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 0)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamCurrent(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 1)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorPos(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 2)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorAng(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 3)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertPos(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 4)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertAng(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 5)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamLongPos(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 6)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorSizeRMS(ElecName)
string ElecName
return sqrt(srwGetElecBeamParam(ElecName, 20))
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorMixedMom(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 21)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorDivergRMS(ElecName)
string ElecName
return sqrt(srwGetElecBeamParam(ElecName, 22))
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamHorSizeProjRMS(ElecName, rObs)
string ElecName
variable rObs 
rObs += srwGetElecBeamLongPos(ElecName)
variable elecSigXe2 = srwGetElecBeamParam(ElecName, 20)
variable elecSigXpe2 = srwGetElecBeamParam(ElecName, 22)
variable elecMXXp = srwGetElecBeamHorMixedMom(ElecName)
variable elecSigXeffE2 = elecSigXe2 + rObs*rObs*elecSigXpe2 + 2*rObs*elecMXXp
return sqrt(elecSigXeffE2)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertSizeRMS(ElecName)
string ElecName
return sqrt(srwGetElecBeamParam(ElecName, 23))
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertMixedMom(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 24)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertDivergRMS(ElecName)
string ElecName
return sqrt(srwGetElecBeamParam(ElecName, 25))
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamVertSizeProjRMS(ElecName, rObs)
string ElecName
variable rObs
rObs += srwGetElecBeamLongPos(ElecName)
variable elecSigZe2 = srwGetElecBeamParam(ElecName, 23)
variable elecSigZpe2 = srwGetElecBeamParam(ElecName, 25)
variable elecMZZp = srwGetElecBeamHorMixedMom(ElecName)
variable elecSigZeffE2 = elecSigZe2 + rObs*rObs*elecSigZpe2 + 2*rObs*elecMZZp
return sqrt(elecSigZeffE2)
end

//+++++++++++++++++++++++++++++++++++++++
function srwGetElecBeamRelEnSprRMS(ElecName)
string ElecName
return srwGetElecBeamParam(ElecName, 13)
end

//+++++++++++++++++++++++++++++++++++++++
//
//Update second-order Moments from Twiss
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecUpdateMomentsFromTwiss()

Variable emx_m = SrwElecEmx*1.e-9, emz_m = SrwElecEmz*1.e-9, sigeE2 = SrwElecSige*SrwElecSige

SrwElecSigx=(1.e+06)*sqrt(emx_m*SrwElecBetax + sigeE2*SrwElecEtax*SrwElecEtax) //microns
SrwElecSigz=(1.e+06)*sqrt(emz_m*SrwElecBetaz + sigeE2*SrwElecEtaz*SrwElecEtaz) //microns
SrwElecMxxp=(1.e+09)*(-emx_m*SrwElecAlphax + sigeE2*SrwElecEtax*SrwElecEtaxPr) //nm
SrwElecMzzp=(1.e+09)*(-emz_m*SrwElecAlphaz + sigeE2*SrwElecEtaz*SrwElecEtazPr) //nm
SrwElecSigxp=(1.e+06)*sqrt(emx_m*(1 + SrwElecAlphax*SrwElecAlphax)/SrwElecBetax + sigeE2*SrwElecEtaxPr*SrwElecEtaxPr) //micro-r
SrwElecSigzp=(1.e+06)*sqrt(emz_m*(1 + SrwElecAlphaz*SrwElecAlphaz)/SrwElecBetaz + sigeE2*SrwElecEtazPr*SrwElecEtazPr) //micro-r
End

//+++++++++++++++++++++++++++++++++++++++
//
//Update Twiss param from second-order Moments
//Attention: dispersion/disp. deriv. is taken from old values !
//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecUpdateTwissFromMoments()

Variable sigx_m = SrwElecSigx*1.e-6, sigz_m = SrwElecSigz*1.e-6
Variable sigxp_r = SrwElecSigxp*1.e-6, sigzp_r = SrwElecSigzp*1.e-6
Variable mxxp_m = SrwElecMxxp*1.e-9, mzzp_m = SrwElecMzzp*1.e-9

Variable sigeE2 = SrwElecSige*SrwElecSige
Variable sigx_mE2 = sigx_m*sigx_m, sigz_mE2 = sigz_m*sigz_m
Variable sigxp_rE2 = sigxp_r*sigxp_r, sigzp_rE2 = sigzp_r*sigzp_r

//Variable emx_m = sqrt(sigx_mE2*sigxp_rE2 - mxxp_m*mxxp_m)
Variable emx_mE2 = sigx_mE2*sigxp_rE2 - mxxp_m*mxxp_m - sigeE2*(SrwElecEtaxPr*SrwElecEtaxPr*sigx_mE2 + SrwElecEtax*SrwElecEtax*sigxp_rE2 - 2.*mxxp_m*SrwElecEtax*SrwElecEtaxPr)
//double Emit2 = Sigma2*SigmaPrime2 - MixMom2 - SigmaE2*(EtaPrime2*Sigma2 + Eta2*SigmaPrime2 - 2*MixMom*Eta*EtaPrime);
Variable emx_m = 0
if(emx_mE2 > 0.)
	emx_m = sqrt(emx_mE2)
endif

//Variable emz_m = sqrt(sigz_mE2*sigzp_rE2 - mzzp_m*mzzp_m)
Variable emz_mE2 = sigz_mE2*sigzp_rE2 - mzzp_m*mzzp_m - sigeE2*(SrwElecEtazPr*SrwElecEtazPr*sigz_mE2 + SrwElecEtaz*SrwElecEtaz*sigzp_rE2 - 2.*mzzp_m*SrwElecEtaz*SrwElecEtazPr)
//double Emit2 = Sigma2*SigmaPrime2 - MixMom2 - SigmaE2*(EtaPrime2*Sigma2 + Eta2*SigmaPrime2 - 2*MixMom*Eta*EtaPrime);
Variable emz_m = 0
if(emz_mE2 > 0.)
	emz_m = sqrt(emz_mE2)
endif

Variable emx_mOld = SrwElecEmx*(1.e-9), emz_mOld = SrwElecEmz*(1.e-9)

if(abs(emx_m - emx_mOld) > 0.0001*emx_m)
	Variable TermSigX = (1.e-12)*SrwElecSigx*SrwElecSigx
	//SrwElecEtax = 0
	//SrwElecEtaxPr = 0
	SrwElecBetax = TermSigX/emx_m
	SrwElecAlphax = -(1.e-09)*SrwElecMxxp/emx_m
	SrwElecEmx = (1.e+09)*emx_m
endif
if(abs(emz_m - emz_mOld) > 0.0001*emz_m)
	Variable TermSigZ = (1.e-12)*SrwElecSigz*SrwElecSigz
	//SrwElecEtaz = 0
	//SrwElecEtazPr = 0
	SrwElecBetaz = TermSigZ/emz_m
	SrwElecAlphaz = -(1.e-09)*SrwElecMzzp/emz_m
	SrwElecEmz = (1.e+09)*emz_m
endif
End

//+++++++++++++++++++++++++++++++++++++++
//
//Electron Beam Panel
//
//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecDialog(SrcType)
Variable SrcType // 0- for SR, 1- for Gaussian Beam, 2- for Isotropic Source

SrwElecCreateBufVars(SrcType)
DoWindow/K SrwElecPanel
SrwElecPanel()

End

//+++++++++++++++++++++++++++++++++++++++
Window SrwElecPanel() : Panel

String/G SrwPElecTitle
PauseUpdate; Silent 1		// building window...
NewPanel /W=(210,65,655,502) as SrwPElecTitle
SetDrawLayer UserBack

SrwElecDrawAllContr()

EndMacro

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecDrawAllContr()

Variable VertOffset=12
Variable xStartNames1=20, xStartVars1=135
Variable xStartNames2=230, xStartVars2=345

String AllElecWaves=Wavelist("*"+SrwElecType,";","")
String ElecNameTot=SrwElecName+SrwElecType
Variable ItemNo=sRWaveListItemNo(AllElecWaves, ";", ElecNameTot)

//if(ItemNo == 0)
//	ItemNo=1
//	ElecNameTot=sRWaveListItemName(AllElecWaves, ";", 1)
//	SrwElecName=ElecNameTot[0,strlen(ElecNameTot)-strlen(SrwElecType)-1]
//else
//	SrwElecSetupGlobVars(ElecNameTot)
//endif

DrawRect 10,VertOffset-2,435,55+VertOffset

// Existing structures
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText xStartNames1,20+VertOffset,SrwPElecExist
PopupMenu popup0Elec,pos={xStartVars1,4+VertOffset},size={255,21}
PopupMenu popup0Elec,value=#"Wavelist(\"*\"+SrwElecType,\";\",\"\")",mode=ItemNo,proc=SrwElecSelectPopupProc

// Name
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText xStartNames1,48+VertOffset,"Name"
SetVariable setvar18Elec,pos={xStartVars1,30+VertOffset},size={180,17},title=" ",fSize=14;
SetVariable setvar18Elec,limits={-Inf,Inf,1},value=SrwElecName

// Thin beam
DrawRect 10,60+VertOffset,435,195+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText xStartNames1-5,77+VertOffset,"\"Filament\" beam"

if(SrwElecSrcTypeBuf==0)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1,97+VertOffset,SrwPElecEn
	SetVariable setvar0Elec,pos={xStartVars1,80+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar0Elec,limits={-Inf,Inf,1},value=SrwElecEn
endif

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames2,97+VertOffset,SrwPElecCur
SetVariable setvar1Elec,pos={xStartVars2,80+VertOffset},size={80,14},title=" ",fSize=14
SetVariable setvar1Elec,limits={-Inf,Inf,1},value=SrwElecCur

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames1,121+VertOffset,"Longitudinal Position where Initial Conditions are defined [m]"
SetVariable setvar2Elec,pos={xStartVars2,105+VertOffset},size={80,14},title=" ",fSize=14
SetVariable setvar2Elec,limits={-Inf,Inf,1},value=SrwElecs0,proc=SrwElecS0SetVarProc

DrawRect xStartNames1-5,130+VertOffset,220,190+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames1,145+VertOffset,"Horizontal"

DrawRect xStartNames2-5,130+VertOffset,430,190+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames2,145+VertOffset,"Vertical"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames1+10,162+VertOffset,"Position [mm]"
SetVariable setvar3Elec,pos={xStartVars1,145+VertOffset},size={80,14},title=" ",fSize=14
SetVariable setvar3Elec,limits={-Inf,Inf,1},value=SrwElecxx

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames2+10,162+VertOffset,"Position [mm]"
SetVariable setvar4Elec,pos={xStartVars2,145+VertOffset},size={80,14},title=" ",fSize=14
SetVariable setvar4Elec,limits={-Inf,Inf,1},value=SrwEleczz

if(SrwElecSrcTypeBuf != 2)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,182+VertOffset,"Angle [mr]"
	SetVariable setvar5Elec,pos={xStartVars1,165+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar5Elec,limits={-Inf,Inf,1},value=SrwElecxp
endif

if(SrwElecSrcTypeBuf != 2)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,182+VertOffset,"Angle [mr]"
	SetVariable setvar6Elec,pos={xStartVars2,165+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar6Elec,limits={-Inf,Inf,1},value=SrwEleczp
endif

// Thick beam
DrawRect 10,200+VertOffset,435,390+VertOffset
SetDrawEnv fname= "default",fsize=12,fstyle=1
DrawText xStartNames1-5,217+VertOffset,"\"Thick\" beam"

SetDrawEnv fname= "default",fsize=12,fstyle=0
DrawText xStartNames1,237+VertOffset,"Definition by"
PopupMenu popup1Elec,pos={xStartVars1,220+VertOffset},size={255,21}
if(SrwElecSrcTypeBuf != 2)
	PopupMenu popup1Elec,value="Twiss;Moments",mode=SrwElecThickDefBy,proc=SrwElecDefinedByPopupProc
endif
if(SrwElecSrcTypeBuf == 2)
	PopupMenu popup1Elec,value="Moments",mode=1
	SrwElecThickDefBy=2
endif

if(SrwElecSrcTypeBuf == 0)
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2,237+VertOffset,SrwPElecSige
	SetVariable setvar7Elec,pos={xStartVars2,220+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar7Elec,limits={-Inf,Inf,1},value=SrwElecSige
endif

if(SrwElecThickDefBy==1) // Twiss

	DrawRect xStartNames1-5,245+VertOffset,220,365+VertOffset
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1,260+VertOffset,"Horizontal"

	DrawRect xStartNames2-5,245+VertOffset,430,365+VertOffset
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2,260+VertOffset,"Vertical"

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,277+VertOffset,"Emittance [nm]"
	SetVariable setvar8Elec,pos={xStartVars1,260+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar8Elec,limits={-Inf,Inf,1},value=SrwElecEmx
	
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,277+VertOffset,"Emittance [nm]"
	SetVariable setvar9Elec,pos={xStartVars2,260+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar9Elec,limits={-Inf,Inf,1},value=SrwElecEmz

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,297+VertOffset,"Beta [m]"
	SetVariable setvar10Elec,pos={xStartVars1,280+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar10Elec,limits={-Inf,Inf,1},value=SrwElecBetax

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,297+VertOffset,"Beta [m]"
	SetVariable setvar11Elec,pos={xStartVars2,280+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar11Elec,limits={-Inf,Inf,1},value=SrwElecBetaz

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,317+VertOffset,"Alpha [r]"
	SetVariable setvar12Elec,pos={xStartVars1,300+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar12Elec,limits={-Inf,Inf,1},value=SrwElecAlphax

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,317+VertOffset,"Alpha [r]"
	SetVariable setvar13Elec,pos={xStartVars2,300+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar13Elec,limits={-Inf,Inf,1},value=SrwElecAlphaz

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,337+VertOffset,"Dispersion [m]"
	SetVariable setvar14Elec,pos={xStartVars1,320+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar14Elec,limits={-Inf,Inf,1},value=SrwElecEtax

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,337+VertOffset,"Dispersion [m]"
	SetVariable setvar15Elec,pos={xStartVars2,320+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar15Elec,limits={-Inf,Inf,1},value=SrwElecEtaz

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,357+VertOffset,"Dispers. Deriv. [r]"
	SetVariable setvar16Elec,pos={xStartVars1,340+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar16Elec,limits={-Inf,Inf,1},value=SrwElecEtaxPr

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,357+VertOffset,"Dispers. Deriv. [r]"
	SetVariable setvar17Elec,pos={xStartVars2,340+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar17Elec,limits={-Inf,Inf,1},value=SrwElecEtazPr

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	//String StrS0="NOTE:"
	//DrawText xStartNames1,385+VertOffset,StrS0
	String StrS1="Twiss parameters refer to the Longitudinal Position:  "+num2str(SrwElecs0)+"  m"
	DrawText xStartNames1,385+VertOffset,StrS1

endif
if(SrwElecThickDefBy==2) // Moments

	DrawRect xStartNames1-5,245+VertOffset,220,325+VertOffset
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1,260+VertOffset,"Horizontal"

	DrawRect xStartNames2-5,245+VertOffset,430,325+VertOffset
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2,260+VertOffset,"Vertical"

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames1+10,277+VertOffset,"RMS Size [µm]"
	SetVariable setvar8Elec,pos={xStartVars1,260+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar8Elec,limits={-Inf,Inf,1},value=SrwElecSigx
	
	SetDrawEnv fname= "default",fsize=12,fstyle=0
	DrawText xStartNames2+10,277+VertOffset,"RMS Size [µm]"
	SetVariable setvar9Elec,pos={xStartVars2,260+VertOffset},size={80,14},title=" ",fSize=14
	SetVariable setvar9Elec,limits={-Inf,Inf,1},value=SrwElecSigz

	if(SrwElecSrcTypeBuf != 2)

		SetDrawEnv fname= "default",fsize=12,fstyle=0
		DrawText xStartNames1+10,297+VertOffset,"RMS Diverg. [µr]"
		SetVariable setvar10Elec,pos={xStartVars1,280+VertOffset},size={80,14},title=" ",fSize=14
		SetVariable setvar10Elec,limits={-Inf,Inf,1},value=SrwElecSigxp

		SetDrawEnv fname= "default",fsize=12,fstyle=0
		DrawText xStartNames2+10,297+VertOffset,"RMS Diverg. [µr]"
		SetVariable setvar11Elec,pos={xStartVars2,280+VertOffset},size={80,14},title=" ",fSize=14
		SetVariable setvar11Elec,limits={-Inf,Inf,1},value=SrwElecSigzp

		SetDrawEnv fname= "default",fsize=11,fstyle=0
		DrawText xStartNames1-3,317+VertOffset,"<(x-<x>)(x'-<x'>)>[nm]"
		SetVariable setvar12Elec,pos={xStartVars1,300+VertOffset},size={80,14},title=" ",fSize=14
		SetVariable setvar12Elec,limits={-Inf,Inf,1},value=SrwElecMxxp

		SetDrawEnv fname= "default",fsize=11,fstyle=0
		DrawText xStartNames2-3,317+VertOffset,"<(x-<x>)(x'-<x'>)>[nm]"
		SetVariable setvar13Elec,pos={xStartVars2,300+VertOffset},size={80,14},title=" ",fSize=14
		SetVariable setvar13Elec,limits={-Inf,Inf,1},value=SrwElecMzzp
		
		SetDrawEnv fname= "default",fsize=12,fstyle=0
		String StrS1="NOTE:"
		DrawText xStartNames1,370+VertOffset,StrS1
		String StrS2="RMS Size = <(x-<x>)^2>^(1/2);  RMS Divergence = <(x'-<x'>)^2>^(1/2)"
		DrawText xStartNames1+40,370+VertOffset,StrS2
		String StrS3="<(x-<x>)(x'-<x'>)>^2 = <(x-<x>)^2><(x'-<x'>)^2> - Emittance^2"
		DrawText xStartNames1+40,385+VertOffset,StrS3
		
	endif

	SetDrawEnv fname= "default",fsize=12,fstyle=0
	String StrS5="Second-order Moments refer to the Longitudinal Position:  "+num2str(SrwElecs0)+"  m"
	DrawText xStartNames1,347+VertOffset,StrS5

endif

// OK-Cancel-Help
Button button0Elec,pos={85,398+VertOffset},size={70,20},proc=SrwElecCancelButtonProc,title=SrwPPanelCancelButton
Button button1Elec,pos={185,398+VertOffset},size={70,20},proc=SrwElecOKButtonProc,title=SrwPPanelOKButton
Button button2Elec,pos={285,398+VertOffset},size={70,20},proc=SrwElecHelpButtonProc,title=SrwPPanelHelpButton

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecKillAllContr()

KillControl popup0Elec
KillControl popup1Elec
KillControl setvar0Elec
KillControl setvar1Elec
KillControl setvar2Elec
KillControl setvar3Elec
KillControl setvar4Elec
KillControl setvar5Elec
KillControl setvar6Elec
KillControl setvar7Elec
KillControl setvar8Elec
KillControl setvar9Elec
KillControl setvar10Elec
KillControl setvar11Elec
KillControl setvar12Elec
KillControl setvar13Elec
KillControl setvar14Elec
KillControl setvar15Elec
KillControl setvar16Elec
KillControl setvar17Elec
KillControl setvar18Elec

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecCreateBufVars(SrcType)
Variable SrcType

Variable/G SrwElecSrcTypeBuf=SrcType

Variable/G SrwElecSigeBuf=SrwElecSige

Variable/G SrwElecSigxBuf=SrwElecSigx
Variable/G SrwElecSigzBuf=SrwElecSigz
Variable/G SrwElecSigxpBuf=SrwElecSigxp
Variable/G SrwElecSigzpBuf=SrwElecSigzp
Variable/G SrwElecMxxpBuf=SrwElecMxxp
Variable/G SrwElecMzzpBuf=SrwElecMzzp

Variable/G SrwElecEmxBuf=SrwElecEmx
Variable/G SrwElecEmzBuf=SrwElecEmz
Variable/G SrwElecBetaxBuf=SrwElecBetax
Variable/G SrwElecBetazBuf=SrwElecBetaz
Variable/G SrwElecAlphaxBuf=SrwElecAlphax
Variable/G SrwElecAlphazBuf=SrwElecAlphaz
Variable/G SrwElecEtaxBuf=SrwElecEtax
Variable/G SrwElecEtazBuf=SrwElecEtaz
Variable/G SrwElecEtaxPrBuf=SrwElecEtaxPr
Variable/G SrwElecEtazPrBuf=SrwElecEtazPr

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecKillBufVars()

KillVariables/Z SrwElecSrcTypeBuf

KillVariables/Z SrwElecSigeBuf
KillVariables/Z SrwElecSigxBuf, SrwElecSigzBuf, SrwElecSigxpBuf, SrwElecSigzpBuf, SrwElecMxxpBuf, SrwElecMzzpBuf
KillVariables/Z SrwElecEmxBuf, SrwElecEmzBuf, SrwElecBetaxBuf, SrwElecBetazBuf, SrwElecAlphaxBuf, SrwElecAlphazBuf, SrwElecEtaxBuf, SrwElecEtazBuf, SrwElecEtaxPrBuf, SrwElecEtazPrBuf

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecSelectPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

SrwElecSetupGlobVars(popStr)

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwElecKillAllContr()
SrwElecDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecSetupGlobVars(ElecName)
String ElecName

SrwElecName=ElecName[0,strlen(ElecName)-strlen(SrwElecType)-1]

SrwElecEn=$ElecName[0]
SrwElecCur=$ElecName[1]
SrwElecxx=($ElecName[2])*1000.
SrwElecxp=($ElecName[3])*1000.
SrwEleczz=($ElecName[4])*1000.
SrwEleczp=($ElecName[5])*1000.
SrwElecs0=$ElecName[6]

SrwElecEmx=$ElecName[7]
SrwElecBetax=$ElecName[8]
SrwElecAlphax=$ElecName[9]
SrwElecEmz=$ElecName[10]
SrwElecBetaz=$ElecName[11]
SrwElecAlphaz=$ElecName[12]
SrwElecSige=$ElecName[13]
SrwElecEtax=$ElecName[14]
SrwElecEtaxPr=$ElecName[15]
SrwElecEtaz=$ElecName[16]
SrwElecEtazPr=$ElecName[17]

Variable MomSetup=$ElecName[18]
if(MomSetup==0)
	SrwElecThickDefBy=1
else
	SrwElecThickDefBy=2
endif

SrwElecSigx=(1.e+06)*sqrt(abs($ElecName[20]))
SrwElecSigz=(1.e+06)*sqrt(abs($ElecName[23]))
SrwElecSigxp=(1.e+06)*sqrt(abs($ElecName[22]))
SrwElecSigzp=(1.e+06)*sqrt(abs($ElecName[25]))
SrwElecMxxp=(1.e+09)*($ElecName[21])
SrwElecMzzp=(1.e+09)*($ElecName[24])

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecRestoreGlobVars()

SrwElecSige=SrwElecSigeBuf

SrwElecSigx=SrwElecSigxBuf
SrwElecSigz=SrwElecSigzBuf
SrwElecSigxp=SrwElecSigxpBuf
SrwElecSigzp=SrwElecSigzpBuf
SrwElecMxxp=SrwElecMxxpBuf
SrwElecMzzp=SrwElecMzzpBuf

SrwElecEmx=SrwElecEmxBuf
SrwElecEmz=SrwElecEmzBuf
SrwElecBetax=SrwElecBetaxBuf
SrwElecBetaz=SrwElecBetazBuf
SrwElecAlphax=SrwElecAlphaxBuf
SrwElecAlphaz=SrwElecAlphazBuf
SrwElecEtax=SrwElecEtaxBuf
SrwElecEtaz=SrwElecEtazBuf
SrwElecEtaxPr=SrwElecEtaxPrBuf
SrwElecEtazPr=SrwElecEtazPrBuf

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecDefinedByPopupProc(ctrlName,popNum,popStr) : PopupMenuControl
String ctrlName
Variable popNum		// which item is currently selected (1-based)
String popStr			// contents of current popup item as string

Variable OldDefBy=SrwElecThickDefBy
SrwElecThickDefBy=popNum

if((OldDefBy==1) %& (SrwElecThickDefBy==2)) // from Twiss to Moments
	SrwElecUpdateMomentsFromTwiss()
endif
if((OldDefBy==2) %& (SrwElecThickDefBy==1)) // from Moments to Twiss
	SrwElecUpdateTwissFromMoments()
endif

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwElecKillAllContr()
SrwElecDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecS0SetVarProc(ctrlName,varNum,varStr,varName) : SetVariableControl
String ctrlName
Variable varNum	// value of variable as number
String varStr		// value of variable as string
String varName	// name of variable

SetDrawLayer/K UserBack
SetDrawLayer UserBack
SrwElecKillAllContr()
SrwElecDrawAllContr()

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecCancelButtonProc(ctrlName) : ButtonControl
String ctrlName

SrwElecKillAllContr()
DoWindow/K SrwElecPanel
SrwElecRestoreGlobVars()
SrwElecKillBufVars()
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecOKButtonProc(ctrlName) : ButtonControl
String ctrlName

// Validate parameters
SrwElecValidateGlobVars()

String ComLineStr1, ComLineStr2
sprintf ComLineStr1, "SrwElecFilament(\"%s\",%g,%g,%g,%g,%g,%g,%g)", SrwElecName, SrwElecEn, SrwElecCur, SrwElecs0, SrwElecxx, SrwElecxp, SrwEleczz, SrwEleczp
String ElecName=SrwElecName+SrwElecType
if(SrwElecThickDefBy==1)
	sprintf ComLineStr2, "SrwElecThick(\"%s\",%g,%g,%g,%g,%g,%g,%g,%g,%g)", ElecName, SrwElecSige, SrwElecEmx, SrwElecEmz, SrwElecBetax, SrwElecBetaz, SrwElecAlphax, SrwElecAlphaz, SrwElecEtax, SrwElecEtaxPr
	String ComLineStr=ComLineStr1+";"+ComLineStr2
	print ComLineStr
	SrwElecFilament(SrwElecName, SrwElecEn, SrwElecCur, SrwElecs0, SrwElecxx, SrwElecxp, SrwEleczz, SrwEleczp)
	SrwElecThick(ElecName, SrwElecSige, SrwElecEmx, SrwElecEmz, SrwElecBetax, SrwElecBetaz, SrwElecAlphax, SrwElecAlphaz, SrwElecEtax, SrwElecEtaxPr)
endif
if(SrwElecThickDefBy==2)
	sprintf ComLineStr2, "SrwElecThickMom(\"%s\",%g,%g,%g,%g,%g,%g,%g)", ElecName, SrwElecSige, SrwElecSigx, SrwElecSigz, SrwElecSigxp, SrwElecSigzp, SrwElecMxxp, SrwElecMzzp
	String ComLineStr=ComLineStr1+";"+ComLineStr2
	print ComLineStr
	SrwElecFilament(SrwElecName, SrwElecEn, SrwElecCur, SrwElecs0, SrwElecxx, SrwElecxp, SrwEleczz, SrwEleczp)
	SrwElecThickMom(ElecName, SrwElecSige, SrwElecSigx, SrwElecSigz, SrwElecSigxp, SrwElecSigzp, SrwElecMxxp, SrwElecMzzp)
endif

SrwElecKillBufVars()
SrwElecKillAllContr()
DoWindow/K SrwElecPanel

End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecHelpButtonProc(ctrlName) : ButtonControl
String ctrlName
srwUtiShowHelpTopic("SrwElecDialog     ")
End

//+++++++++++++++++++++++++++++++++++++++
Proc SrwElecValidateGlobVars()

// Filament
if(SrwElecEn <= 0.)
	Abort SrwPAlertElecEnergy
endif

// Thick
if(SrwElecSige < 0.)
	Abort SrwPAlertElecEnSpr
endif
if(SrwElecThickDefBy==1)
	if((SrwElecEmx < 0.) %| (SrwElecEmz < 0.))
		Abort SrwPAlertElecEmittance
	endif
	if((SrwElecBetax < 0.) %| (SrwElecBetaz < 0.))
		Abort SrwPAlertElecBeta
	endif
endif
if(SrwElecThickDefBy==2)
	if((SrwElecSigx < 0.) %| (SrwElecSigz < 0.))
		Abort SrwPAlertElecSize
	endif
	if((SrwElecSigxp < 0.) %| (SrwElecSigzp < 0.))
		Abort SrwPAlertElecDiverg
	endif
	Variable MaxAbsMxxp = (1.e-03)*SrwElecSigx*SrwElecSigxp;
	Variable MaxAbsMzzp = (1.e-03)*SrwElecSigz*SrwElecSigzp
	if(abs(SrwElecMxxp) > MaxAbsMxxp)
		Abort SrwPAlertElecMxxp
	endif
	if(abs(SrwElecMzzp) > MaxAbsMzzp)
		Abort SrwPAlertElecMzzp
	endif
endif

End
