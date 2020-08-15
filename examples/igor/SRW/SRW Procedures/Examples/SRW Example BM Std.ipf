#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//
// Example: Bending Magnet
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamStdBM()

srwUtiShowHelpTopic("Standard Bending Magnet Radiation     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=6, current=0.2, s0=-0.49
SrwElecFilament("Elec",energy,current,s0,0,0,0,0)

Variable enspread=0.001
Variable emx=3.9, betax=35.6
Variable emz=0.039, betaz=2.5
SrwElecThick("Elec_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Constant Magnetic Field
Variable Bconst=0.85
SrwMagConstCreate("BM",Bconst,0)

// Set up Observation
SrwSmpCreate("ObsE",30);SrwSmpScanXZE("ObsE_obs",0,1,1,0,1,1,0.1,100,200)
SrwSmpCreate("ObsXZ",30);SrwSmpScanXZE("ObsXZ_obs",0,20,100,0,20,100,10,10,1)

// Compute Spectrum
SrwStoConstCreate("ElecBMObsE","Elec_ebm","BM_mac","ObsE_obs")
// Extract data
SrwSto2Int("ElecBMObsE_ras","I",1,8,5.05,0,0,2)
Textbox/N=text0/F=0 " On-axis spectrum\r 30 m from BM "

// Compute Spectral Angular Distribution
SrwStoConstCreate("ElecBMObsXZ","Elec_ebm","BM_mac","ObsXZ_obs")
// Extract data
SrwSto2Int("ElecBMObsXZ_ras","Iv",2,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux / Surface \r 30 m from BM; Phot. energy: 10 keV \r Vertical polarization "
SrwSto2Int("ElecBMObsXZ_ras","Ih",1,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux / Surface \r 30 m from BM; Phot. energy: 10 keV \r Horizontal polarization "
SrwSto2Int("ElecBMObsXZ_ras","J",7,3,10,0,0,2)
Textbox/N=text0/F=0 " Spectral flux / Surface \r 30 m from BM; Phot. energy: 10 keV "

// Compute Power Density
SrwPowCreate("ElecBMObsXZ","Elec_ebm","BM_mac","ObsXZ_obs",1,2,2)
Textbox/N=text0/F=0 " Power per unit surface 30 m from BM "
SrwPow2Int("ElecBMObsXZ_pow","Jp",2,0,0,2)
Textbox/N=text0/F=0 " Power per unit surface 30 m from BM "

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("Standard Bending Magnet Radiation     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example BM Std.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end