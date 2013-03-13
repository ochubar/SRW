#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//
// Example: Ellipsoidal Wiggler
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamWigElliptic()

srwUtiShowHelpTopic("Ellipsoidal Wiggler     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
Print "Computing Wiggler Radiation. Be patient, It may take from a few 10s of seconds to a few minutes."
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=6, current=0.2
SrwElecFilament("Elec",energy,current,0,0,0,0,0)

Variable enspread=0.001
Variable emx=3.9, betax=35.6
Variable emz=0.039, betaz=2.5
SrwElecThick("Elec_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Elliptic Wiggler
Variable Per=150, Len=0.45, Kz=20, Kx=1.2
//SrwMagPerCreate_("Wig",Per,1,Kz,Len)
//SrwMagPerAddHarm("Wig_map",1,2,Kx,1.571)
//SrwMagPerCreate2D_("Wig",Per,Kz,Kx,Len,1.571)
SrwMagPerCreate2D("Wig",Per,Kz,Kx,Len,1.571,1,0,0)

// Set up Observation
SrwSmpCreate("Obs",30)
SrwSmpScanXZE("Obs_obs",0,110,60,0,20,40,30,30,1)

// Compute SR
SrwStoWigCreate("ElecWigObs","Elec_ebm","Wig_map","Obs_obs",1)

// Extract data
SrwSto2Int("ElecWigObs_ras","Icl",6,4,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r 30 m from wiggler Kz=20, Kx=1.2 \r Photon energy: 30 keV \r Circular left polarization "

SrwSto2Int("ElecWigObs_ras","Icr",5,4,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r 30 m from wiggler Kz=20, Kx=1.2 \r Photon energy: 30 keV \r Circular right polarization "

SrwSto2Int("ElecWigObs_ras","Iv",2,4,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r 30 m from wiggler Kz=20, Kx=1.2 \r Photon energy: 30 keV \r Vertical polarization "

SrwSto2Int("ElecWigObs_ras","Ih",1,4,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r 30 m from wiggler Kz=20, Kx=1.2 \r Photon energy: 30 keV \r Horizontal polarization "

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("Ellipsoidal Wiggler     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example WigEll.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end