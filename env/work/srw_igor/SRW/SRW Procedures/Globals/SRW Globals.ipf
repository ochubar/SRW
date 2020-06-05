
//+++++++++++++++++++++++++++++++++++++++
//
// SRW Globals
//
//++++++++++++++ Remarks++++++++++++++++++
//
// All globals and main procs written in Igor start with Srw
// All global functions start with srw
// All functions written in C start with sr
//
//+++++++++++++++++++++++++++++++++++++++

#pragma rtGlobals=1		// Use modern global access method.
//#include <Keyword-Value>		// For NumByKey

//+++++++++++++++++++++++++++++++++++++++
//
//Create and Intialize all global variables
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwInit(v)
Variable v=1
Silent 1						|	Initializing SRW  ...

String/G SrwVerProc="3.97"
String VerExt=srVerNo()
String SrwURL="http://ftp.esrf.fr/pub/InsertionDevices/"

Variable p0=strsearch(SrwVerProc,".",0)+1
Variable p1=strsearch(SrwVerProc,".",p0)+1
String pv0=SrwVerProc[0,p0-2]
String pv1=SrwVerProc[p0,p1-2]

Variable e0=strsearch(VerExt,".",0)+1
Variable e1=strsearch(VerExt,".",e0)+1
String ev0=VerExt[0,e0-2]
String ev1=VerExt[e0,e1-2]

DefaultFont "Times"

if(v==1)
//String info=IgorInfo(0)
//Variable Vers=NumberByKey("IGORVERS", info)
//Print "Igor Pro Version : ", vers

Print "SRW Proc. Version : ", SrwVerProc
Print "SRW Ext. Version : ", VerExt

//Print "Memory Available : ", NumByKey("FREEMEM", info)/1000000, " M"
endif

//print p0,p1, pv0,pv1,SrwVerProc[p1,strlen(SrwVerProc)]
if(cmpstr(ev0,pv0)+cmpstr(ev1,pv1)!=0)
DoAlert  0  "Incompatibility between the SRW Procedure and the SRW Extension files"
Print "Download a new SRW from : "+SrwURL
abort
endif

if(v==1)
print " "
Print "SRWE Examples are accessible from \"Help\" sub-menus of individual types of computation, e.g.:\r      SRWE->Undulator->Help->Example: Spectrum through a Slit"
print " "
Print "SRWP Examples are accessible from the general SRWP \"Help\" sub-menu, e.g.:\r      SRWP->Help->Example: Focusing the Undulator Radiation"
endif

//   Path to SRW folder
//String SrwPathStr = FunctionPath("srwUtiGetValN")
//SrwPathStr = ParseFilePath(1, SrwPathStr, ":", 1, 2)
string SrwPathStr = srwUtiPathStr()
newpath/O/Q SrwUtiPath, SrwPathStr

//   All types of waves
String/G SrwSeparator="_"
String/G SrwElecType=SrwSeparator+"ebm"				// filament electron beams structure
String/G SrwElecContType=SrwSeparator+"ebc"			// container of electron beam structures
String/G SrwFieldWaveType=SrwSeparator+"fld"			// field wave
String/G SrwAngleWaveType=SrwSeparator+"ang"		// angle wave
String/G SrwTrajWaveType=SrwSeparator+"tra"			// trjectory wave
String/G SrwFieldType=SrwSeparator+"mag"				// arbitrary magnetic field structure
String/G SrwUndType=SrwSeparator+"map"				// periodic magnetic field structure
String/G SrwUndHarmType=SrwSeparator+"fha"			// harmonic of periodic magnetic field structure
String/G SrwMagConstType=SrwSeparator+"mac"		// constant magnetic field structure
String/G SrwMagOptType=SrwSeparator+"mgo"			// magnetic lens structure
String/G SrwMagContainerType=SrwSeparator+"mgc"		// container to keep various magnetic fields

//String/G SrwMagType_TabTrUnif="MagTabTransvUnif"
//String/G SrwMagType_Periodic="MagPeriodic"
//String/G SrwMagType_Quad="MagQuad"
//String/G SrwMagType_Const="MagConst"
//String/G SrwMagType_Group="MagGroup"

String/G SrwPowType=SrwSeparator+"pow"				// power density structure

String/G SrwSmpType=SrwSeparator+"obs"				// radiation sampling structure
String/G SrwSmpPowType=SrwSeparator+"obp"			// radiation sampling structure for power density

String/G SrwRadType=SrwSeparator+"rad"				// radiation field structure
String/G SrwRadElType=SrwSeparator+"rae"				// radiation field 3D complex wave
String/G SrwRadXType="x"								// Wave vs Horizontal Position
String/G SrwRadZType="z"								// Wave vs Vertical Position
String/G SrwRadEType="e"								// Wave vs Energy
String/G SrwRadXZType=SrwSeparator+SrwRadXType+SrwRadZType

String/G SrwSuffixExtract="I"							// For Stokes Waves
String/G SrwSuffixMagField="B"						// For Stokes Waves
String/G SrwSuffixX="X"									
String/G SrwSuffixZ="Z"
String/G SrwSuffixY="Y"
String/G SrwSuffixE="E"

String/G SrwAngType="A"								// Angle Type
String/G SrwPosType="P"								// Position Type
String/G SrwMat44Type=SrwSeparator+"m44"				// 4x4 matrices
String/G SrwStoType=SrwSeparator+"ras"				// Stokes parameters

String/G SrwMomType=SrwSeparator+"mom"				// first-second order Moments
String/G SrwWfrAuxDataType=SrwSeparator+"wfx"		// Auxiliary Data far Wave Front Manupulations
String/G SrwBeamlineType=SrwSeparator+"bli"			// BeamLine type
String/G SrwOptThinGenWaveType=SrwSeparator+"bgt"			// Generic transmission numeric wave type

Variable/G SrwVal=1;
Variable/G SrwInputUnit=1,SrwOutputUnit=4

Variable/G SrwSpecFluxVal=1;
Variable/G SrwSpecFluxInUnit=1,SrwSpecFluxOutUnit=2

Variable/G SrwAllowPrintingExtraInfo=1;String/G SrwPAllowPrintingExtraInfo="Printing Extra Information to History";
String/G SrwPOPUPAllowPrintingExtraInfo="Allow;Suppress";

//   SrwElec
String/G SrwElecName="Elec",SrwPElecName="Name of the Electron Beam structure"
String/G SrwPElecName1="Electron Beam structure"
Variable/G SrwElecEn=6;String/G SrwPElecEn="Electron Energy [GeV]";
Variable/G SrwElecCur=0.2;String/G SrwPElecCur="Electron Current [A]";

String/G SrwElecName="Elec"
Variable/G SrwElecxx=0;String/G SrwPElecxx="Horizontal Pos. [mm]"
Variable/G SrwElecxp=0;String/G SrwPElecxp="Horizontal Angle [mr]"
Variable/G SrwEleczz=0;String/G SrwPEleczz="Vertical Pos. [mm]"
Variable/G SrwEleczp=0;String/G SrwPEleczp="Vertical Angle [mr]"
Variable/G SrwElecs0=0;String/G SrwPElecs0="Longitudinal Position [m]"

Variable/G SrwElecSige=1e-3;String/G SrwPElecSige="Rel. RMS En. Spread"
Variable/G SrwElecEmx=3.9;String/G SrwPElecEmx="Hor. Emittance [nm]"
Variable/G SrwElecBetax=35.6;String/G SrwPElecBetax="Horizontal Beta [m]"
Variable/G SrwElecAlphax=0.;String/G SrwPElecAlphax="Horizontal Alpha [r]"
Variable/G SrwElecEmz=0.039;String/G SrwPElecEmz="Vert. Emittance [nm]"
Variable/G SrwElecBetaz=2.5;String/G SrwPElecBetaz="Vertical Beta [m]"
Variable/G SrwElecAlphaz=0.;String/G SrwPElecAlphaz="Vertical Alpha [r]"
Variable/G SrwElecEtax=0.;String/G SrwPElecEtax="Hor. Dispersion [m]"
Variable/G SrwElecEtaxPr=0.;String/G SrwPElecEtaxPr="Hor. Disp. Deriv. [r]"
Variable/G SrwElecEtaz=0.;String/G SrwPElecEtaz="Vert. Dispersion [m]"
Variable/G SrwElecEtazPr=0.;String/G SrwPElecEtazPr="Vert. Disp. Deriv. [r]"

Variable/G SrwElecSigx=373;String/G SrwPElecSigx="RMS Size [µm]";String/G SrwPElecSigx2="Horizontal RMS Size [µm]"
Variable/G SrwElecSigz=9.9;String/G SrwPElecSigz="RMS Size [µm]";String/G SrwPElecSigz2="Vertical RMS Size [µm]"
Variable/G SrwElecSigxp=10.5;String/G SrwPElecSigxp="RMS Diverg. [µr]";String/G SrwPElecSigxp2="Horizontal RMS Diverg. [µr]"
Variable/G SrwElecSigzp=3.95;String/G SrwPElecSigzp="RMS Diverg. [µr]";String/G SrwPElecSigzp2="Vertical RMS Diverg. [µr]"
Variable/G SrwElecMxxp=0;String/G SrwPElecMxxp="<(x-<x>)(x'-<x'>)> [nm]"
Variable/G SrwElecMzzp=0;String/G SrwPElecMzzp="<(z-<z>)(z'-<z'>)> [nm]"

Variable/G SrwElecThickDefBy=1

//   SrwMagFieldCreate and SrwMagPrec
String/G SrwMagName="Mag";String/G SrwPMagName="Name of the Magnetic Field structure"
String/G SrwBxName="Bx";String/G SrwPBxName="Horizontal Component of the Field"
String/G SrwBzName="Bz";String/G SrwPBzName="Vertical Component of the Field"
Variable/G SrwFieldCenter=0;String/G SrwPFieldCenter="Center of the Field [m] "
Variable/G SrwFieldLength=6.6;String/G SrwPFieldLength="Length of the Field [m] "
Variable/G SrwFieldNpts=1500;String/G SrwPFieldNpts="Number of Points "

//   SrwMagPrec
String/G SrwPMagName1="Name of the Magnetic Field structure to be modified"
Variable/G SrwMode=3;String/G SrwPMode="Method for Longitudinal Integration";String/G SrwPOPUPMode="Manual;Auto Undulator;Auto Wiggler"
Variable/G SrwRadIntStep=0.01; String/G SrwPRadIntStep="Long. Step [m] (for \"Man\" only)"
Variable/G SrwPrec=0.01; String/G SrwPPrec="Rel. Prec. (for \"Auto\" only)"
Variable/G SrwMaxPtsToSave=10000; String/G SrwPMaxPtsToSave="Num. of points to keep in Memory"
Variable/G SrwUseDiffRadIntLimits=1; String/G SrwPUseDiffRadIntLimits="Use Special Integration Limits";String/G SrwPOPUPUseDiffRadIntLimits="No;Yes"
Variable/G SrwSdep=0; String/G SrwPSdep="Initial Longitudinal Coordinate [m]"
Variable/G SrwSfin=0; String/G  SrwPsfin="Final Longitudinal Coordinate [m]"
Variable/G SrwFieldElecS0=(SrwFieldCenter - 0.5*SrwFieldLength)/1.2; String/G SrwPFieldElecS0="Long. Pos. [m]  with  Zero Angle and Pos."

//   SrwMagSin
String/G SrwPFldName="Name of the Field Component"
Variable/G SrwMagZeroSgn=2;String/G SrwPMagZero="Zero the Field before ?"
Variable/G SrwMagScent=0;String/G SrwPMagScent="Center Point [m]"
Variable/G SrwPeriod=35;String/G SrwPPeriod="Period [mm]"
Variable/G SrwMagNper=46.5;String/G SrwPMagNper="N. of Per. [Integer/2]"
Variable/G SrwMagBpeak=0.7;String/G SrwPMagBpeak="Peak Field [T]"
Variable/G SrwMagBnoise=0;String/G SrwPMagBnoise="Rms Noise [T]"
Variable/G SrwMagBnoisei=0;String/G SrwPMagBnoisei="Rms Integrated Noise [T]"
Variable/G SrwMagTaper=0;String/G SrwPMagTaper="Taper = Gap diff. bw. Entry and Exit [mm]"
Variable/G SrwMagSag=0;String/G SrwPMagSag="Sag = Gap diff. bw. Center and Extrem. [mm]"

//   SrwMagEdge
Variable/G SrwMagScentEdge=0
Variable/G SrwMagStrSectLen=SrwFieldLength/1.1; String/G SrwPMagStrSectLen="Length of Straight Section [m]"
Variable/G SrwMagFringeSize=50; String/G SrwPMagFringeSize="Fringe Field Length (10-90%) [mm]"
Variable/G SrwMagBconst=1.56;String/G SrwPMagBconst="Uniform Field [T]"
Variable/G SrwMagDipoleLen=2 //[m]

//   SrwMagGsnAng
Variable/G SrwMagScentGsnAng=0
Variable/G SrwMagSigmaGsnAng=50; String/G SrwPMagSigmaGsnAng="Rms Length [mm]"
Variable/G SrwFieldIntGsnAng=1; String/G SrwPFieldIntGsnAng="Field Integral [T*mm]"

//   SrwMagImportCmpn
Variable/G SrwMagImportUnits=1; String/G SrwPMagImportUnits="Units of Data to Import", SrwPOPUPMagImportUnits="Tesla;Gauss"
String/G SrwMagImportFileName=""; String/G SrwPMagImportFileName="Full File Name with Path (if known)"; 

String/G SrwMagBname=SrwBzName,SrwPMagBname="Name of the Field Component"

//   Constant Magnetic Field
String/G SrwMagConstName="BM"
String/G SrwPMagConstName="Name of the Constant Magnetic Field structure"
Variable/G SrwMagConstBz=0.85; String/G SrwPMagConstBz="Vertical Magnetic Field [T]"
Variable/G SrwMagConstBx=0.; String/G SrwPMagConstBx="Horizontal Magnetic Field [T]"

//   SrwMagElecTraj
Variable/G SrwDisplayAngX=2, SrwDisplayTraX=2, SrwDisplayAngZ=2, SrwDisplayTraZ=2

String/G SrwSmpName="Obs";String/G SrwPSmpName="Name of the Radiation Sampling structure"
Variable/G SrwSmpDist=30;String/G SrwPSmpDist="Distance between Source and Observation Plane [m]"
Variable/G SrwSmpEdep=0.1; String/G SrwPSmpEdep="Initial Photon Energy [keV]"
Variable/G SrwSmpEfin=20; String/G SrwPSmpEfin="Final Photon Energy [keV]"
Variable/G SrwSmpEnpts=1000; String/G SrwPSmpEnpts="Number of Energy Points"

String/G SrwPSmpName1="Radiation Sampling structure"
Variable/G SrwSmpXmid=0; String/G SrwPSmpXmid="Middle Horizontal Position [mm]"
Variable/G SrwSmpXmagn=1; String/G SrwPSmpXmagn="Range of Horizontal Position [mm]"
Variable/G SrwSmpXnpts=1; String/G SrwPSmpXnpts="Number of Horizontal Points"
Variable/G SrwSmpZmid=0; String/G SrwPSmpZmid="Middle Vertical Position [mm]"
Variable/G SrwSmpZmagn=1; String/G SrwPSmpZmagn="Range of Vertical Position [mm]"
Variable/G SrwSmpZnpts=1; String/G SrwPSmpZnpts="Number of Vertical Points"

//   SrwWfr
String/G SrwRadName="Wfr",SrwPRadName="Name of the Wavefront structure"
String/G SrwPElecName2="Filament Electron Beam structure"
String/G SrwPSmpName2="Radiation Sampling structure"
String/G SrwPMagName2="Magnetic Field structure"
Variable/G SrwModeField=1; 
String/G SrwPModeField="Type of Computation "
String/G SrwPSmpNxNzForProp="Use Automatic Radiation Sampling"
String/G SrwPSmpNxNzSamplFact="Oversampling Factor"

Variable/G SrwSmpNxNzForProp=1;
Variable/G SrwSmpNxNzSamplFact=1;

//   SrwViewRes
Variable/G SrwViewRadCmpnType=7
String/G SrwPOPUPPolar="Linear Hor.;Linear Vert.;Linear 45;Linear 135;Circular Right;Circular Left"

Variable/G SrwViewPlotType=8
Variable/G SrwViewPlotTypeXZ=4
Variable/G SrwViewE=(SrwSmpEdep+SrwSmpEfin)/2
Variable/G SrwViewX=SrwSmpXmid
Variable/G SrwViewZ=SrwSmpZmid
Variable/G SrwViewDisplay=2

String/G SrwPViewRadName="Wavefront structure"
String/G SrwPViewSuffixExtract="Suffix"
String/G SrwPViewPlotType="As a function of"
String/G SrwPOPUPViewPlotType="Phot. Energy (or Time);Horizontal Pos.;Vertical Pos.;Hor. + Vert. Pos.;Phot. En. (or Time) + Hor. Pos.;Phot. En. (or Time) + Vert. Pos.;Phot. En. (or Time) + Hor. + Vert.;Auto"
String/G SrwPOPUPViewPlotTypeXZ="Horizontal Pos.;Vertical Pos.;Hor. + Vert. Pos.;Auto"
String/G SrwPViewE="Const. Photon Energy [keV] (or Time [fs])"
String/G SrwPViewX= "Const. Horizontal Pos. [mm]"
String/G SrwPViewZ= "Const. Vertical Pos. [mm]"
String/G SrwViewRadName=""

Variable/G SrwRepr=1;String/G SrwPRepr="Tranverse Representation"

Variable/G SrwCmpnNo=1;String/G SrwPCmpnNo="Single-e or Multi-e";String/G SrwPOPUPCmpnNo="Single-e Intens.;Multi-e Intens.;Single-e Phase;Single-e Re(E);Single-e Flux;Multi-e Flux;Single-e Im(E)" //"Single-e Intens.;Multi-e Intens.;Single-e Phase;Single-e Re(E);Single-e Flux;Multi-e Flux"
String/G SrwPViewRadCmpnType="Polarization Component"
String/G SrwPViewDisplay="New Display";String/G SrwPOPUPViewDisplay="No;Yes"
Variable/G SrwViewRadConv2PolRate=1

//   Other Light Sources
Variable/G SrwIsotrSrcPhot=1000;String/G SrwPIsotrSrcPhot="Phot./0.1%bw emitted by one electron"

//   Optics
String/G SrwPBli="Name of the Optical Element"
String/G SrwBliDrift="Drift";String/G SrwBliDriftType="Drift";
Variable/G SrwBliDriftLength=1;String/G SrwPBliDriftLength="Length of the Drift Space [m]"

String/G SrwBliThinLens="ThinLens";String/G SrwBliThinLensType="ThinLens"
Variable/G SrwBliTlFocalX=1;String/G SrwPBliTlFocalX="Horizontal Focal Length [m]"
Variable/G SrwBliTlFocalZ=1;String/G SrwPBliTlFocalZ="Vertical Focal Length [m]"
Variable/G SrwBliTlPosX=0;String/G SrwPBliTlPosX="Horizontal Position [mm]"
Variable/G SrwBliTlPosZ=0;String/G SrwPBliTlPosZ="Vertical Position [mm]"

String/G SrwBliRectApert="RectAperture";String/G SrwBliRectApertType="RectAperture"
Variable/G SrwBliRectApertDx=10;String/G SrwPBliRectApertDx="Horizontal Size [mm]"
Variable/G SrwBliRectApertDz=10;String/G SrwPBliRectApertDz="Vertical Size [mm]"
Variable/G SrwBliRectApertPosX=0;String/G SrwPBliRectApertPosX="Horizontal Position [mm]"
Variable/G SrwBliRectApertPosZ=0;String/G SrwPBliRectApertPosZ="Vertical Position [mm]"

String/G SrwBliCircApert="CircAperture";String/G SrwBliCircApertType="CircAperture"
Variable/G SrwBliCircApertD=10;String/G SrwPBliCircApertD="Diameter [mm]"
Variable/G SrwBliCircApertPosX=0;String/G SrwPBliCircApertPosX="Horizontal Position [mm]"
Variable/G SrwBliCircApertPosZ=0;String/G SrwPBliCircApertPosZ="Vertical Position [mm]"

String/G SrwBliObstRectType="RectObstacle"
String/G SrwBliObstCircType="CircObstacle"

String/G SrwBliSpherMirror="SpherMirror";String/G SrwBliSpherMirrorType="SpherMirror";
Variable/G SrwBliSpherMirrorR=5;String/G SrwPBliSpherMirrorR="Radius of Curvature [m]";
Variable/G SrwBliSpherMirrorDx=20;String/G SrwPBliSpherMirrorDx="Horizontal Size [mm]";
Variable/G SrwBliSpherMirrorDz=20;String/G SrwPBliSpherMirrorDz="Vertical Size [mm]";
Variable/G SrwBliSpherMirrorRotPlane=1;String/G SrwPBliSpherMirrorRotPlane="Rotation Plane";
Variable/G SrwBliSpherMirrorTheta=2.3562;String/G SrwPBliSpherMirrorTheta="Angle bw. Norm. and Opt. Axis [r]";
Variable/G SrwBliSpherMirrorPosX=0;String/G SrwPBliSpherMirrorPosX="Horizontal Position [mm]"
Variable/G SrwBliSpherMirrorPosZ=0;String/G SrwPBliSpherMirrorPosZ="Vertical Position [mm]"

String/G SrwBliAngle="Angle";String/G SrwBliAngType="Angle";
Variable/G SrwBliAngleX=0;String/G SrwPBliAngleX="Horizontal Angle [mrad]"
Variable/G SrwBliAngleZ=0;String/G SrwPBliAngleZ="Vertical Angle [mrad]"

String/G SrwBliCont="Container",SrwBliContType="Container"
String/G SrwBliLast="",SrwPBliLast="Element to be Added"
String/G SrwPBliCont="Container"
Variable/G SrwBliwhere=2;String/G SrwPBliwhere="Where?"

String/G SrwBliWgRect="RectWaveguide";String/G SrwBliWgRectType="WaveguideRect"
Variable/G SrwBliWgRectLen=1;String/G SrwPBliWgRectLen="Length [m]"
Variable/G SrwBliWgRectDx=10;String/G SrwPBliWgRectDx="Horizontal Size [mm]"
Variable/G SrwBliWgRectDz=10;String/G SrwPBliWgRectDz="Vertical Size [mm]"
Variable/G SrwBliWgRectPosX=0;String/G SrwPBliWgRectPosX="Horizontal Position [mm]"
Variable/G SrwBliWgRectPosZ=0;String/G SrwPBliWgRectPosZ="Vertical Position [mm]"

//   Generic Transmission Optical Element
String/G SrwBliThinGenType="ThinGen"
String/G SrwBliThinGen="ThinGen"
String/G SrwPBliThinGenWave="Name of 2D Complex Transmission wave",SrwPBliThinGenWave1="2D Complex Transmission wave"
//Variable/G SrwBliOptPathOrPhase=1;String/G SrwPBliOptPathOrPhase="Im Part of Transmission is";String/G SrwPOPUPBliOptPathOrPhase="Optical Path Difference;Phase Shift"
Variable/G SrwBliThinGenOuter=1;String/G SrwPBliThinGenOuter="Transmission in Outer Region is";String/G SrwPOPUPBliThinGenOuter="Zero;Same as at the Border"
Variable/G SrwBliThinGenE=20.;String/G SrwPBliThinGenE="Photon Energy [keV]"
Variable/G SrwBliThinGenPosX=0.;String/G SrwPBliThinGenPosX="Horizontal Position [mm]"
Variable/G SrwBliThinGenPosZ=0.;String/G SrwPBliThinGenPosZ="Vertical Position [mm]"
Variable/G SrwBliThinGenWaveNx=100;String/G SrwPBliThinGenWaveNx="Horizontal Number of Points"
Variable/G SrwBliThinGenWavePosX=0.;String/G SrwPBliThinGenWavePosX="Horizontal Position [mm]"
Variable/G SrwBliThinGenWaveRangeX=1.;String/G SrwPBliThinGenWaveRangeX="Horizontal Range [mm]"
Variable/G SrwBliThinGenWaveNz=100;String/G SrwPBliThinGenWaveNz="Vertical Number of Points"
Variable/G SrwBliThinGenWavePosZ=0.;String/G SrwPBliThinGenWavePosZ="Vertical Position [mm]"
Variable/G SrwBliThinGenWaveRangeZ=1.;String/G SrwPBliThinGenWaveRangeZ="Vertical Range [mm]"
String/G SrwBliThinGenWaveFunTrans="TransmFunc",SrwPBliThinGenWaveFunTrans="Name of Amplitude Transmission function"
String/G SrwBliThinGenWaveFunPhase="PhaseFunc",SrwPBliThinGenWaveFunPhase="Name of Phase Shift Function"
String/G SrwBliThinGenWaveFunOptPath="OptPathFunc",SrwPBliThinGenWaveFunOptPath="Name of Optical Path Difference function"

String/G SrwPThinGenViewWave="Name of the wave to create"
Variable/G SrwBliThinGenViewCmpn=1;String/G SrwPBliThinGenViewCmpn="Component to extract";String/G SrwPOPUPBliThinViewCmpn="Amplitude Transmission;Intensity Transmission;Optical Path Difference";String/G SrwPOPUPBliThinGenViewCmpn="Amplitude Transmission;Intensity Transmission;Phase Shift;Re(T);Im(T)"
Variable/G SrwBliThinGenViewPlotType=1;String/G SrwPBliThinGenViewPlotType="As a function of";String/G SrwPOPUPBliThinGenViewPlotType="Hor. and Vert. Position;Horizontal Position;Vertical Position"
Variable/G SrwBliThinGenViewNewDisp=2;String/G SrwPBliThinGenViewNewDisp="New Display";String/G SrwPOPUPBliThinGenViewNewDisp="No;Yes"
String/G SrwBliThinGenDefDispName="Transm"

//   X-ray Lens
String/G SrwBliXrayLensType="XrayLens"
String/G SrwBliXrayLens="XrayLens"
Variable/G SrwBliXrayLensMat=1;String/G SrwPBliXrayLensMat="Material";String/G SrwPOPUPBliXrayLensMat="Beryllium;PyroCarbon;Other"
Variable/G SrwBliXrayLensFocPlane=1;String/G SrwPBliXrayLensFocPlane="Focusing Plane";String/G SrwPOPUPBliXrayLensFocPlane="Horizontal;Vertical;Horizontal + Vertical"
Variable/G SrwBliXrayLensShape=1;String/G SrwPBliXrayLensShape="Hole Shape";String/G SrwPOPUPBliXrayLensShape="Circular;Parabolic"
Variable/G SrwBliXrayLensDiam=1.;String/G SrwPBliXrayLensDiam="Hole Diameter [mm]";String/G SrwPBliXrayLensDiam1="Aperture [mm]"
Variable/G SrwBliXrayLensRmin=1.;String/G SrwPBliXrayLensRmin="Min. Radius [mm]"
Variable/G SrwBliXrayLensWallThick=0.1;String/G SrwPBliXrayLensWallThick="Wall Thickness [mm]"
Variable/G SrwBliXrayLensNhol=10;String/G SrwPBliXrayLensNhol="Number of Holes"
Variable/G SrwBliXrayLensPos=0;String/G SrwPBliXrayLensPos="Transverse Position [mm]"
Variable/G SrwBliXrayLensAttenLen=1.083;String/G SrwPBliXrayLensAttenLen="Inten. Atten. Len. [mm] (for \"Other\" mat.)"
Variable/G SrwBliXrayLensDensity=1.845;String/G SrwPBliXrayLensDensity="Density [g/cm3]"
Variable/G SrwBliXrayLensAtomNum=4;String/G SrwPBliXrayLensAtomNum="Atomic Number (for \"Other\" mat.)"
Variable/G SrwBliXrayLensDelta=5.44E-07
Variable/G SrwBliXrayLensPhotEn
Variable/G SrwBliXrayLensXc, SrwBliXrayLensZc

//   Simple Circular Zone Plate
String/G SrwBliZonePlate="ZonePlate";String/G SrwBliZonePlateType="ZonePlate"
Variable/G SrwBliZonePlateFocDist=1.;String/G SrwPBliZonePlateFocDist="Focal Distance [m]"
Variable/G SrwBliZonePlateAttenLenSb=1.;String/G SrwPBliZonePlateAttenLenSb="Substrate Intens. Atten. Length [mm]"
Variable/G SrwBliZonePlateThickSb=0.1;String/G SrwPBliZonePlateThickSb="Substrate Thickness [mm]"
Variable/G SrwBliZonePlateAttenLenMn=1.;String/G SrwPBliZonePlateAttenLenMn="Zone Layer Intens. Atten. Length [mm]"
Variable/G SrwBliZonePlateDeltaRefrMn=1.e-06;String/G SrwPBliZonePlateDeltaRefrMn="Zone Layer Refraction (1 - n)"
Variable/G SrwBliZonePlateIgnThickMn=1.;String/G SrwPBliZonePlateIgnThickMn="Choice of Zone Layer Thickness";String/G SrwPOPUPBliZonePlateIgnThickMn="Set for Best Focusing;Define Explicitly"
Variable/G SrwBliZonePlateThickMn=1.;String/G SrwPBliZonePlateThickMn="Zone Layer Thickness [mm]"
Variable/G SrwBliZonePlateDelPhi=1;String/G SrwPBliZonePlateDelPhi="Phase Shift in successive zones";String/G SrwPOPUPBliZonePlateDelPhi="Pi;3 Pi;5 Pi"
Variable/G SrwBliZonePlateNzones=50;String/G SrwPBliZonePlateNzones="Total Number of Zones"
Variable/G SrwBliZonePlatePhotEn
Variable/G SrwBliZonePlatePhaseMult

Variable/G SrwBliZonePlateAttenLen1=100.;String/G SrwPBliZonePlateAttenLen1="Material 1 Atten. Length [mm]"
Variable/G SrwBliZonePlateDeltaRefr1=1.e-06;String/G SrwPBliZonePlateDeltaRefr1="Material 1 Delta (= 1 - n)"
Variable/G SrwBliZonePlateAttenLen2=1.e-03;String/G SrwPBliZonePlateAttenLen2="Material 2 Atten. Length [mm]"
Variable/G SrwBliZonePlateDeltaRefr2=1.e-06;String/G SrwPBliZonePlateDeltaRefr2="Material 2 Delta (= 1 - n)"
Variable/G SrwBliZonePlateIgnThick=1;String/G SrwPBliZonePlateIgnThick="Choice of Thickness";String/G SrwPOPUPBliZonePlateIgnThick="Best Efficiency;Entered Value"
Variable/G SrwBliZonePlateThick=1.;String/G SrwPBliZonePlateThick="Thickness [mm]"
Variable/G SrwBliZonePlateXc, SrwBliZonePlateZc
Variable/G SrwBliZonePlateRn=0.1
Variable/G SrwBliZonePlatedRn=0.001

//   Estimate Refraction
Variable/G SrwBliMatZ=4; String/G SrwPBliMatZ="Atomic Number"
Variable/G SrwBliMatDens=1.845; String/G SrwPBliMatDens="Density [g/cm3]"
Variable/G SrwBliZonePlateDeltaRefr;

//   SrwWfrResize
Variable/G SrwRadResizeMeth=1; String/G SrwPRadResizeMeth="Resizing Method"
Variable/G SrwRadResizeScanX=1; String/G SrwPRadResizeScanX="Horizontal Range Resizing"
Variable/G SrwRadResizeStepX=1; String/G SrwPRadResizeStepX="Horizontal Resolution Resizing"
Variable/G SrwRadResizeScanZ=1; String/G SrwPRadResizeScanZ="Vertical Range Resizing"
Variable/G SrwRadResizeStepZ=1; String/G SrwPRadResizeStepZ="Vertical Resolution Resizing"
Variable/G SrwDplBeforeResize=1; String/G SrwPDplBeforeResize="Duplicate Wavefront before Resizing?"

String/G SrwBliResize="Resize",SrwBliResType="Resize"
Variable/G SrwBliResmode=1;String/G SrwPBliResmode="Mode ?"
Variable/G SrwBliResdepx=0;String/G SrwPBliResdepx="Horizontal Initial Coordinate [m or q]"
Variable/G SrwBliResfinx=0;String/G SrwPBliResfinx="Horizontal Final Coordinate [m or q]"
Variable/G SrwBliResnpx=-1;String/G SrwPBliResnpx="Number of Horizontal Points"
Variable/G SrwBliResdepz=0;String/G SrwPBliResdepz="Vertical Initial Coordinate [m or q]"
Variable/G SrwBliResfinz=0;String/G SrwPBliResfinz="Vertical Final Coordinate [m or q]"
Variable/G SrwBliResnpz=-1;String/G SrwPBliResnpz="Number of Vertical Points"

//   SrwRadProp
String/G SrwPRadName1="Wavefront Structure"
String/G SrwPBliName="Beamline / Optical Component"
Variable/G SrwPropMeth=1;String/G SrwPPropMeth="Auto-Resize Wavefront?"
Variable/G SrwDplBeforeProp=1;String/G SrwPDplBeforeProp="Duplicate Wavefront Before Propagation?"
String/G SrwPRadNameDpl="Name of Duplicated Wavefront structure"

String/G SrwSuffixStokes="S",SrwPSuffixStokes="Suffix"
String/G SrwStoName="",SrwPStoName="Name of the Stokes structure",SrwPStoName1="Stokes structure"
Variable/G SrwFluxType=1;String/G SrwPFluxType="Output Units"
String/G SrwPStoNameDpl="Name of Duplicated Stokes structure"

String/G SrwUndName="U",SrwPUndName="Name of the Periodic Field structure",SrwPUndName2="Periodic Field structure"
Variable/G SrwLength=3.2;String/G SrwPLength="Length [m]"
Variable/G SrwUndKind=1;String/G SrwPUndKind="Type of Undulator"
Variable/G SrwUndTaper=0;String/G SrwPUndTaper="N dE/E (for Tapered Undulator only)"
Variable/G SrwUndOptKlystPhaseShift=0;String/G SrwPUndOptKlystPhaseShift="Nd/N (for Optical Klystron only)"

Variable/G SrwPlane=1;String/G SrwPPlane="Plane of Field"
Variable/G SrwKK=2.2;String/G SrwPKK="Deflection Parameter (Per. x B x 0.0934)"
Variable/G SrwHarm=1;String/G SrwPHarm="Harmonic Number"
Variable/G SrwPhase=Pi/2;String/G SrwPPhase="Phase [radian]"

Variable/G SrwKz=2.2;String/G SrwPKz="Kz = 0.0934 x Per. x Bz"
Variable/G SrwKx=0.;String/G SrwPKx="Kx = 0.0934 x Per. x Bx"
Variable/G SrwPh0x=1.571;String/G SrwPPh0x="Phase bw. Vert. and Hor. [r]"

Variable/G SrwPerInitHarm=1;String/G SrwPPerInitHarm="Initial Harmonic of Spectrum"
Variable/G SrwPerFinHarm=7;String/G SrwPPerFinHarm="Final Harmonic of Spectrum"
Variable/G SrwPerNs=1;String/G SrwPPerNs="Longitudinal Integ. Param."
Variable/G SrwPerNphi=1;String/G SrwPPerNphi="Azimuthal Integ. Param."

Variable/G SrwGap=11;String/G SrwPGap="Magnetic Gap [mm]"
Variable/G SrwAir=0.05;String/G SrwPAir="Air Space between Blocks [mm]"
Variable/G SrwMagPer=4;String/G SrwPMagPer="Number of Magnets per Period"
Variable/G SrwBr=1.17;String/G SrwPBr="Magnetization of the Blocks [T]"
Variable/G SrwHeight=12;String/G SrwPHeight="Height of the Blocks [mm]"

Variable/G SrwUndFieldType=1

String/G SrwPowName="",SrwPPowName="Name of the Power Density structure",SrwPPowName1="Power Density structure"
String/G SrwPPowNameDpl="Name of duplicated Power Density structure"
Variable/G SrwPowCompPrec=1.;String/G SrwPPowCompPrec="Precision Parameter"
Variable/G SrwPowCompMeth=1;String/G SrwPPowCompMeth="Computation Method"
Variable/G SrwPowDispImmed=2;String/G SrwPPowDispImmed="Display Results"

String/G SrwMagGenTotName=""
String/G SrwRadGenTotName=""
String/G SrwPRadGenName="Radiation structure"
String/G SrwPRadGenNameDpl="Name of duplicated Radiation structure"

String/G SrwSmpGenTotName=""

String/G SrwStoWigName="",SrwPStoWigName="Name of the Wiggler structure"
Variable/G SrwStoWigCompPrec=1.;String/G SrwPStoWigCompPrec="Precision Parameter"

String/G SrwStoConstName="",SrwPStoConstName="Name of the Wiggler structure"

Variable/G SrwStoConstCompPrec=1.;String/G SrwPStoConstCompPrec="Precision Parameter"

//   General for Panels
String/G SrwPPanelCancelButton="Quit";
String/G SrwPPanelOKButton="Continue";
String/G SrwPPanelHelpButton="Help";

//   SrwRadSamplPanel
String/G SrwPRadSamplTitle="Radiation Sampling";
String/G SrwPRadSamplExist="Existing structures";
String/G SrwPRadSamplName="Name";
String/G SrwPRadSamplLongPos="Position [m]";
String/G SrwPRadSamplHorPosMid="Center [mm]";
String/G SrwPRadSamplHorPosRng="Range [mm]";
String/G SrwPRadSamplHorPtNum="Number of Pts.";
String/G SrwPRadSamplVerPosMid="Center [mm]";
String/G SrwPRadSamplVerPosRng="Range [mm]";
String/G SrwPRadSamplVerPtNum="Number of Pts.";
String/G SrwPRadSamplEnInit="Initial [keV]";
String/G SrwPRadSamplEnFin="Final [keV]";
String/G SrwPRadSamplEnPtNum="Number of Pts.";

Variable/G SrwStartMacrosAfterRadSmp=0; // 0- no, 1- SrwWfrCreate_(), 2- SrwWfrCreate(), 3- SrwPerStoCreate()
Variable/G SrwStartMacrosAfterRadSmp2=0; // 0- no, 1- SrwWfr2Int_(), 2- SrwWfr2Int(), 3- SrwSto2Int()

//   SrwUtiViewStructPanel
Variable/G SrwUtiViewStructListTypeNo=1;
String/G SrwPViewStructTitle="Structure to View / Modify";
String/G SrwPViewStructType="Type";
String/G SrwPViewStructExist="Existing structures";

//   SrwVisualizePanel
String/G SrwPVisualizeTitle="Visualize"
String/G SrwPRadAllKind="Computed Radiation structure"

//   Titles of Radiation structures
String/G SrwRadTitleStokes="Stokes"
String/G SrwRadTitleElecField="Electric Field"
String/G SrwRadTitlePowDens="Power Density"

//   SrwElecPanel
String/G SrwPElecTitle="Electron Beam"
String/G SrwPElecExist="Existing structures"

//   Temporary 
Variable/G SrwWfrPropDefaultMethNo=2

//   Alert Messages
String/G SrwNoRadName="Please first create a Radiation Field structure"
String/G SrwNoElecName="Please first create an Electron Beam structure"
String/G SrwNoMagName="Please first create a Magnetic Field structure"
String/G SrwNoObsName="Please first create a Radiation Sampling structure"
String/G SrwNoFieldName="Please first create a Field wave"

String/G SrwPAlertNoElec="No Electron Beam structure found."
String/G SrwPAlertNoMag="No Magnetic Field structure found."
String/G SrwPAlertNoTrj="No Trajectory structure found."
String/G SrwPAlertNoSmp="No Radiation Sampling structure found."

String/G SrwPAlertNoCompResultsFound="No relevant SR computation results found. Make sure you have done any SR computation."
String/G SrwPAlertMagFieldNeeded="You need to set up magnetic field now. After this, you may proceed to the SR computation."
String/G SrwPAlertGenMagFieldNeeded="You need to set up magnetic field components now. After this, you may proceed to the SR computation."
String/G SrwPAlertRadSmplNeeded="You need to set up Radiation Sampling structure now. After this, you may proceed to the SR computation."
String/G SrwPAlertWavefrontNeeded="No Wavefront structure found. Make sure you have done any SR computation."
String/G SrwPAlertRadiationNeeded="No Radiation structure found. Make sure you have done any SR computation."
String/G SrwPAlertNoOptCompFound="No Optical Component structure found. Make sure you have defined any Optical Component."
String/G SrwPAlertBadPolRate="Can not compute polarization rate for this radiation component."

String/G SrwPAlertElecBeamNeeded="You need to set up Electron Beam structure now. After this you may proceed to computation of the radiation."
String/G SrwPAlertElecEnergy="Electron beam energy should be positive."
String/G SrwPAlertElecEnSpr="Electron beam energy spread can not be negative."
String/G SrwPAlertElecEmittance="Electron beam emittance can not be negative."
String/G SrwPAlertElecBeta="Beta function can not be negative."
String/G SrwPAlertElecSize="Electron beam size can not be negative."
String/G SrwPAlertTransvSize="Transverse size can not be negative."
String/G SrwPAlertElecDiverg="Electron beam divergence can not be negative."
String/G SrwPAlertElecMxxp="Violation of the requirement:\r<(x-<x>)(x'-<x'>)>^2 <= <(x-<x>)^2><(x'-<x'>)^2> for horizontal moments."
String/G SrwPAlertElecMzzp="Violation of the requirement:\r<(x-<x>)(x'-<x'>)>^2 <= <(x-<x>)^2><(x'-<x'>)^2> for vertical moments."
String/G SrwPAlertElecThickPar="\"Thick\" electron beam parameters were not properly set up."

String/G SrwPAlertGsnBeamNeeded="You need to set up Gaussian Beam structure now. After this, you may proceed to computation of the radiation."
String/G SrwPAlertIsotrSrcPhot="Number of photons should be positive."
String/G SrwPAlertGsnBeamWaist="Gaussian beam waist size can not be negative."
String/G SrwPAlertGsnBeamOrder="The order of a gaussian beam mode can not be negative."
String/G SrwPAlertGsnBeamPolar="Polarization component identifier should be an integer number from 1 to 6."

String/G SrwPAlertPrecStep="Step of integration should be positive."
String/G SrwPAlertPrecRel="Relative precision should be positive."
String/G SrwPAlertPrecIntLim="Improper definition of the integration limits."

String/G SrwPAlertUndPeriod="Undulator period should be positive."
String/G SrwPAlertUndKK="Undulator deflection parameter should be positive."
String/G SrwPAlertUndLength="Undulator length should be positive."
String/G SrwPAlertUndLengthSm="Undulator length should be larger than period."
String/G SrwPAlertUndHarm="Periodic field harmonic number should be positive integer."
String/G SrwPAlertUndGap="Undulator gap should be positive."
String/G SrwPAlertUndHeight="Height of magnet blocks should be positive."
String/G SrwPAlertUndAir="Air space between blocks should be positive."
String/G SrwPAlertUndMagPer="Number of magnet blocks per period should be positive integer."
String/G SrwPAlertUndRadHarm="Undulator radiation harmonic number should be positive integer."
String/G SrwPAlertUndPrecNs="Longitudinal integration parameter should be positive."
String/G SrwPAlertUndPrecNphi="Azimuthal integration parameter should be positive."
String/G SrwPAlertMagPerBad="Improper definition of the Periodic Magnetic Field structure."

String/G SrwPAlertBadName="Improper name of the structure."
String/G SrwPAlertTooLongName="Too long name of the structure. The name should not exceed 28 characters."

String/G SrwPAlertSmpRange="Sampling range can not be negative."
String/G SrwPAlertSmpEnergy="Photon energy can not be negative."
String/G SrwPAlertSmpEnergyNpts="Number of energy points should be positive integer."
String/G SrwPAlertSmpNpts="Number of points should be integer."
String/G SrwPAlertBadNxNzSamplFact="Oversampling factor should be positive real."

String/G SrwPAlertBliThinGenWaveAbsent="2D Complex Transmission wave was not found.\rTo set up this optical component, one needs first to prepare and submit a 2D Complex Transmission wave."
String/G SrwPAlertBliThinGenImproperWave="Improper format of Complex Transmission wave. The Complex Transmission should be a 2D complex 32 bit wave scaled in meters in each of its two dimensions."
String/G SrwPAlertBliThinGenWaveNpts="Number of points should be positive."
String/G SrwPAlertBliThinGenWaveRange="Transverse range should be positive."
String/G SrwPAlertBliThinGenFunAbsent="Can not find user-defined function(s) with the specified name(s)."

String/G SrwPAlertVisMultiE="Extraction of multi-electron intensity requires electric field computed on a grid with many points in horizontal and vertical directions."

//   Units & Graph Labels
String/G SrwPUnitSpAngFluxPerUnSurf="Ph/s/0.1%bw/mm\\S2\\M", SrwPUnitSpAngFluxPerUnSurf1="Ph/s/0.1%bw/mm^2"
String/G SrwPUnitSpAngFluxPerUnAngle="Ph/s/0.1%bw/mr\\S2\\M", SrwPUnitSpAngFluxPerUnAngle1="Ph/s/0.1%bw/mr^2"
String/G SrwPUnitBrilliance="Ph/s/0.1%bw/mm\\S2\\M/mr\\S2\\M", SrwPUnitBrilliance1="Ph/s/0.1%bw/mm^2/mr^2"
String/G SrwPUnitSpAngFlux="Ph/s/0.1%bw"
String/G SrwPUnitSpAngFluxPerUnPlAng="Ph/s/0.1%bw/mr"
String/G SrwPUnitSpAngFluxPerUnLen="Ph/s/0.1%bw/mm"

String/G SrwPUnitArbitrary="Arb. Units"

String/G SrwPUnitPow="W"
String/G SrwPUnitPowDen="W/mm\\S2\\M", SrwPUnitPowDen1="W/mm^2"
String/G SrwPUnitPowPerUnLen="W/mm"

String/G SrwPUnitPhase="radian"
String/G SrwPUnitElectricField="(Ph/s/0.1%bw/mm\\S2\\M)\\S1/2\\M", SrwPUnitElectricField1="(Ph/s/0.1%bw/mm^2)^1/2"
String/G SrwPUnitMagField="T"
String/G SrwPUnitLength="m"
String/G SrwPUnitPhotEn="eV"

String/G SrwPLabelPhotEn="Photon Energy"
String/G SrwPLabelLongPos="Longitudinal Position"
String/G SrwPLabelHorPos="Horizontal Position"
String/G SrwPLabelVerPos="Vertical Position"
String/G SrwPLabelMagField="Magnetic Field"
String/G SrwPLabelHorAngle="Horizontal Angle"
String/G SrwPLabelVerAngle="Vertical Angle"
String/G SrwPLabelPowDen="Power Density"
String/G SrwPLabelPolRate="Polarization Rate"
String/G SrwPLabelTransmAmp="Amplitude Transmission"
String/G SrwPLabelTransmCmplxRe="Re of Complex Transmission"
String/G SrwPLabelTransmCmplxIm="Im of Complex Transmission"

String/G SrwPLabelTransmInt="Intensity Transmission"
String/G SrwPLabelPhaseShift="Phase Shift"
String/G SrwPLabelOptPath="Optical Path Difference"

//   Warning Messages
String/G SrwPWarnMessageHeader="S R W   W A R N I N G"
String/G SrwPWarnFirstSymb="-"

String/G SrwPWarnBliThinGenWaveLim="Dimensions of the optical component are larger than the range of Complex Transmission wave."
String/G SrwPWarnBliThinGenWaveSmp="Number of points in the Complex Transmission wave may be too small to resolve the optical component."

//   Local parts initialization
SrwBrilInit() // Brilliance
SrwGsnBeamInit() // Gaussian Beam
SrwTrjInit() // Trajectory (remove when stabilized)
SrwSASEInit() // SASE
SrwOptThinMirInit() // "Thin" Optics
//SrwOptThickInit() // "Thick" Optics

SrwUtiDataWaveInfCreate() // (Re-)Creates Data Wave Information structure 

End  // proc SrwInit()
