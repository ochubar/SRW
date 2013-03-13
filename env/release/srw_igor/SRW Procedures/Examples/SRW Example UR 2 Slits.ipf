
//+++++++++++++++++++++++++++++++++++++++
//
//Focusing the central cone of the undulator radiation 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamUndRad2Slits()

srwUtiShowHelpTopic("Propagating Undulator Radiation through a Two-Slit Interferometer     ")

Print "Callculating Undulator Radiation. Be patient, It may take from a few 10s of seconds to a few minutes."

// Initialization
SrwInit(2)

variable t0=stopmstimer(0)
t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
variable eBeamEnergy = 7. //[GeV]
variable eBeamVertPos = 0.02 //[mm]
variable eBeamVertAng = 0. //[mrad]
SrwElecFilament("eBeam",eBeamEnergy,0.1,0,0,0,eBeamVertPos,eBeamVertAng)
variable eBeamVertEmit = 0.05 //[nm]
variable eBeamVertBeta = 3 //[m]
SrwElecThick("eBeam_ebm",0.00096,7.7,eBeamVertEmit,16,eBeamVertBeta,0,0,0.13,0)

// Create Undulator Magnetic Field 
SrwMagFieldCreate("U33",0,2.5,5000)
SrwMagZero("U33BX_fld")
variable Nper = 70.5
variable Period = 33 //[mm]
SrwMagSin("U33BZ_fld",1,0,Period,0,0,Nper,0.7,0,0)

// Set Mode of Longitudinal Integration
SrwMagPrec("U33_mag",3,0.001,0.005,10000,1,0,0)

// Display Magnetic Field
SrwMagDisplayField("U33BZ_fld")

// Define Observation parameters for On-Axis Spectrum computation
variable Robs = 30
SrwSmpCreate("O",Robs)
SrwSmpScanXZE("O_obs",0,0,1,0,0,1,0.1,8,1500)

// Compute On-Axis Spectrum and plot it
SrwWfrCreate("S","eBeam_ebm","U33_mag","O_obs",1,1)
SrwWfr2Int("S_rad","I",1,1,8,1,5,0,0,2)
Textbox/N=text1/F=0 " Single-Electron Spectral Flux \r per Unit Surface \r On the Axis "
TextBox/C/N=text1/F=2/D=0.3

// Find harmonic peak of the On-Axis Spectrum
wavestats/Q SI_e
variable en = V_maxloc*0.001 //[keV]

// Define Observation parameters for Electric Field computation vs transverse positions at one photon energy
SrwSmpCreate("OF",Robs)
variable RangeObs = Robs*3000*sqrt(1.23984e-06/(en*Period*Nper)) //[mm]
SrwSmpScanXZE("OF_obs",0,RangeObs,100,0,RangeObs,100,en,en,1)

// Compute Electric Field (in mode appropriate for further propagation) and plot Intensity vs transverse position
SrwWfrCreate("W","eBeam_ebm","U33_mag","OF_obs",2,1.6)
SrwWfr2Int("W_rad","beforeSlits",7,1,3,1,en,0,0,1)
SrwWfr2Int("W_rad","beforeSlits",1,1,4,1,en,0,0,2)
string strAnnotBeforeSlits = " UR Intensity Before Slits \r Single-Electron Emission; Harmonic No: 1 \r Photon Energy: " + num2str(en) + " keV \r Observation Distance: " + num2str(Robs) + " m "
Textbox/N=text1/F=0 strAnnotBeforeSlits

// Create Optical Elements, Propagate the Wavefront and plot intermediate Intensity distributions
//	Ideal Lens:
SrwOptThinLens("ThinLens",Robs/2,Robs/2,0,0)
SrwWfrPropagate("W_rad","ThinLens_bli",1,2,"Wd") //in "automatic" mode, with re-sampling / resizing
//SrwWfrProp("W_rad","ThinLens_bli",1,1,1,2,2,"Wd")

//	Two Slits:
variable SlitSize = 0.07 //[mm]
variable SlitSepar = 0.2 //[mm] 
variable SlitVertCen = 0. //[mm]
SrwOptApertRect("RectAperture",10,SlitSepar + 2*SlitSize,0,SlitVertCen)
SrwOptObstRect("RectObst",10,SlitSepar,0,SlitVertCen)
SrwOptContFull("TwoSlits","RectAperture_bli","RectObst_bli","","") 
SrwWfrPropagate("Wd_rad","TwoSlits_bli",2,1,"") //in "manual" mode, without re-sampling
//SrwWfrProp("Wd_rad","TwoSlits_bli",2,2,1,2,1,"")
SrwWfr2Int("Wd_rad","afterSlits",7,1,4,1,en,0,0,2)
string strAnnotAfterSlits = " UR Intensity After Slits \r Single-Electron Emission "
Textbox/N=text1/F=0 strAnnotAfterSlits
SrwWfr2Int("Wd_rad","afterSlits",7,1,3,1,en,0,0,2)
string strAnnotAfterSlitsCut = " UR Intensity Before \r and After Slits \r Vertical Cuts "
Textbox/N=text1/F=0 strAnnotAfterSlitsCut
TextBox/C/N=text1/F=2/D=0.3
AppendToGraph WbeforeSlits_z
ModifyGraph rgb(WdafterSlits_z)=(0,0,0)

//	Drift space between Slits (Lens) and Image plane:
SrwOptDrift("Drift",Robs)
SrwWfrPropagate("Wd_rad","Drift_bli",1,1,"") //in "automatic" mode, with re-sampling / resizing
//SrwWfrProp("Wd_rad","Drift_bli",1,1,1,2,1,"")

//	Single-e Intensity in the Image plane
SrwWfr2Int("Wd_rad","atImage",7,1,4,1,en,0,0,2)
SetAxis left -0.0003,0.0003
string strAnnotImage = " Single-Electron Intensity in the Image Plane "
Textbox/N=text1/F=0 strAnnotImage

SrwWfr2Int("Wd_rad","atImage",7,1,3,1,en,0,0,2)
SetAxis bottom -0.0003,0.0003
string strAnnotImageCut = " Single-Electron Intensity \r in the Image Plane \r Vertical Cut"
Textbox/N=text1/F=0 strAnnotImageCut
TextBox/C/N=text1/F=2/D=0.3

//	Multi-e Intensity in the Image plane
SrwWfrResize("Wd_rad",1,15,1,1,1,2,"Wdd") //resizing horizontal range
SrwWfr2Int("Wdd_rad","atImageM",7,2,4,1,en,0,0,2)
SetAxis left -0.0003,0.0003
string strAnnotImageM = " Multi-Electron Intensity in the Image Plane \r (estimation by convolution, perhaps not very accurate...) "
Textbox/N=text1/F=0 strAnnotImageM

SrwWfr2Int("Wdd_rad","atImageM",7,2,3,1,en,0,0,2)
SetAxis bottom -0.0003,0.0003
string strAnnotImageCutM = " Multi-Electron Intensity \r in the Image Plane \r Vertical Cut \r (estimation by convolution) "
Textbox/N=text1/F=0 strAnnotImageCutM
TextBox/C/N=text1/F=2/D=0.3

TileWindows/O=1/C

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

srwUtiShowHelpTopic("Propagating Undulator Radiation through a Two-Slit Interferometer     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example UR 2 Slits.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
