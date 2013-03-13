#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//
// Example: Elliptic Wiggler
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamWigPlanar()

srwUtiShowHelpTopic("Planar Wiggler     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=6.04, current=0.2, s0=-0.8
SrwElecFilament("Elec",energy,current,s0,0,0,0,0)

Variable enspread=0.001
Variable emx=4.0, betax=0.5
Variable emz=0.03, betaz=2.73
SrwElecThick("Elec_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Planar Wiggler
Variable Bmax=1.8, Per=150
SrwMagFieldCreate("Mag",0,2,2500)
SrwMagSin("MagBZ_fld",1,0,150,0,0,8.5,1.8,0,0)

// Display Magnetic Field
SrwMagDisplayField("MagBZ_fld")

// Set up Observation
SrwSmpCreate("ObsE1",30);SrwSmpScanXZE("ObsE1_obs",0,1,1,0,1,1,1,130,800)
SrwSmpCreate("ObsE2",30);SrwSmpScanXZE("ObsE2_obs",25,1,1,0,1,1,1,130,100)
SrwSmpCreate("ObsX",30);SrwSmpScanXZE("ObsX_obs",0,150,300,0,1,1,39,39,1)
SrwSmpCreate("ObsZ",30);SrwSmpScanXZE("ObsZ_obs",0,1,1,0,20,250,39,39,1)

// Compute SR
SrwStoWigCreate("ElecMagObsE1","Elec_ebm","Mag_mag","ObsE1_obs",1)
SrwStoWigCreate("ElecMagObsE2","Elec_ebm","Mag_mag","ObsE2_obs",1)

// Extract data
SrwSto2Int("ElecMagObsE1_ras","I",7,8,5.05,0,0,2)
SrwSto2Int("ElecMagObsE2_ras","I",7,8,5.05,0,0,1)
AppendToGraph ElecMagObsE2I_e
ModifyGraph rgb(ElecMagObsE2I_e)=(0,0,0)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r \"Thick\" e-beam (ESRF Low Beta) \r Observation 30 m from wiggler \r   red: on-axis \r   black: off-axis (hor. pos. 25 mm) "

// Compute SR
SrwStoWigCreate("ElecMagObsZ","Elec_ebm","Mag_mag","ObsZ_obs",1)
// Extract data
SrwSto2Int("ElecMagObsZ_ras","I",7,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r in vertical plane 30 m from wiggler \r Photon energy: 39 keV "

// Compute SR
SrwStoWigCreate("ElecMagObsX","Elec_ebm","Mag_mag","ObsX_obs",1)
// Extract data
SrwSto2Int("ElecMagObsX_ras","I",7,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r in horizontal plane 30 m from wiggler\r Photon energy: 39 keV "

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("Planar Wiggler     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example WigPlan.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end