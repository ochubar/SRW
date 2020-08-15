
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Focusing X-rays from Bending Magnet by 
// a Refractive Lens with circular holes
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamXrayLensCirc()

// Initialization
SrwInit(2)

srwUtiShowHelpTopic("Focusing X-rays by a Refractive Lens     ")

Print "Focusing the X-rays from Bending Magnet by a Refractive Lens."
Print "Be patient, it may take from a few seconds to a few minutes."

variable t0=stopmstimer(0)
t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
SrwElecFilament("Elec",6,0.2,0,0,0,0,0)

// Create Magnetic Field
SrwMagFieldCreate("Mag",0,0.2,1500)
SrwMagConst("MagBZ_fld",1.6)

// Setup Mode of Longitudinal Integration
SrwMagPrec("Mag_mag",3,0.03,0.05,3000,1,0,0)

// Radiation Sampling
Variable PhotEnergy = 25 // [keV]
Variable S = 4.5
SrwSmpCreate("Obs",S)
SrwSmpScanXZE("Obs_obs",0,0.05,10,0,0.4,10,PhotEnergy,PhotEnergy,1)

// Estimate distances for desired magnification M
Variable M, F, S1
M = 2.6
F = S/(1+1/M)
S1 = S*F/(S-F) // Distance from lens to geometrical image point

// Refraction and Attenuation taken from
// LBNL Center for X-ray Optics
// http://www-cxro.lbl.gov/optical_constants/
Variable DeltaBe = 5.44841E-07 // refraction index decrement
Variable AttenLenBe = 31.7516 // [mm] at 25 keV

// Hole diameter in mm
Variable D = 0.8 

// Estimate necessary number of holes in the lens for a given material, hole diameter and focal distance
Variable Nhol = round(0.25*D*0.001/(F*DeltaBe)) // Number of holes
Variable Wall = 0.1 // [mm]

// Create Refractive Lens with necessary parameters
SrwOptThinXrayLensCirc("XrayLensCirc",2,DeltaBe,AttenLenBe,D,Nhol,Wall,0,0)

// Display Transmission of the Lens
//SrwOptThinGenDisplay("XrayLensCirc_bgt","XrayLensCirc",2,1,0,0,2)
SrwOptThinTransmDisplay("XrayLensCirc_bgt","XrayLensCirc",2,1,0,0,2)
Textbox/A=RT/N=text0/F=0 " Intensity Transmission \r of the Refractive Lens \r with circular holes. "

//SrwOptThinGenDisplay("XrayLensCirc_bgt","XrayLensCirc",3,3,0,0,2)
SrwOptThinTransmDisplay("XrayLensCirc_bgt","XrayLensCirc",3,3,0,0,2)
Textbox/A=RT/N=text0/F=0 " Optical Path \r in the Refractive Lens \r with circular holes. "

// Create drift space 
// (due to spherical aberration, better focusing is expected a bit closer to the lens than the geometrical image plane)
SrwOptDrift("Drift",S1*0.975)

// Create Rectangular Aperture (to simulate a slit/pinhole in front of the lens)
SrwOptApertRect("RectAperture",0.05,0.4,0,0)

// Create a beamline (Container) and put there all the optical components in proper order
SrwOptContFull("Container","RectAperture_bli","XrayLensCirc_bli","Drift_bli","")

// Compute SR before the lens (at the input aperture)
SrwWfrCreate("ElecMagObs","Elec_ebm","Mag_mag","Obs_obs",2,1)

// Visualize SR at the input
SrwWfr2Int("ElecMagObs_rad","I",7,1,4,1,5.05,0,0,2)
Textbox/A=RT/N=text0/F=0 " Spectral Flux / Surface before the Lens."

// Resize in order to increase precision
SrwWfrResize("ElecMagObs_rad",1,1,2,1,1.2,2,"ElecMagObsd")

// Perform Propagation through the beamline
SrwWfrPropagate("ElecMagObsd_rad","Container_bli",1,2,"ElecMagObsdd")

// Visualize results
SrwWfr2Int("ElecMagObsdd_rad","I",7,1,3,1,5.05,0,0,2)
Textbox/A=RT/N=text0/F=0 " Spectral Flux / Surface in the Image plane. \r Vertical Profile. "

SrwWfr2Int("ElecMagObsdd_rad","I",7,1,2,1,5.05,0,0,2)
Textbox/A=RT/N=text0/F=0 " Spectral Flux / Surface in the Image plane. \r Horizontal Profile. "

SrwWfr2Int("ElecMagObsdd_rad","I",7,1,4,1,5.05,0,0,2)
Textbox/A=RT/N=text0/F=0 " Spectral Flux / Surface in the Image plane. \r Photon Energy 25 keV. "

TileWindows/O=1/C

srwUtiShowHelpTopic("Focusing X-rays by a Refractive Lens     ")
print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example X-rayLensCirc.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end