#pragma rtGlobals=1		// Use modern global access method.

//+++++++++++++++++++++++++++++++++++++++
//
// Example: Asymmetric Wiggler
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamWigAsym()

//srwUtiShowHelpTopic("Asymmetric Wiggler     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=2.75, current=0.5, s0=-1
SrwElecFilament("Elec",energy,current,s0,0,0,0,0)

Variable enspread=0.001
Variable emx=3.73, betax=17.78
Variable emz=0.037, betaz=1.75
SrwElecThick("Elec_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Asymmetric Wiggler (from Gaussian kicks)
Variable Bint=140 //[T*mm]
Variable Per=0.2, sm0=-0.8 //[m] 
Variable PoleSig=20 //[mm]
Variable Nper=7
SrwMagFieldCreate("AsymWig",0,2,3500)
SrwMagGsnAng("AsymWigBZ_fld",1,sm0,PoleSig,0.5*Bint); sm0+=0.25*Per
Variable i=0
do
	SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,-0.5*Bint); sm0+=0.5*Per
	SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,-0.5*Bint); sm0+=0.25*Per
	SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,Bint); sm0+=0.25*Per
	i+=1
while(i<Nper)
SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,-0.5*Bint); sm0+=0.5*Per
SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,-0.5*Bint); sm0+=0.25*Per
SrwMagGsnAng("AsymWigBZ_fld",2,sm0,PoleSig,0.5*Bint)

// Display Magnetic Field
SrwMagDisplayField("AsymWigBZ_fld")
SrwUtiGraphAddFrameAndGrid()

// Setup Observation
SrwSmpCreate("ObsE1",10);SrwSmpScanXZE("ObsE1_obs",0,1,1,0,1,1,0.1,60,4000)
SrwSmpCreate("ObsX",10);SrwSmpScanXZE("ObsX_obs",0,150,301,0,1,1,30,30,1)
SrwSmpCreate("ObsZ",10);SrwSmpScanXZE("ObsZ_obs",0,1,1,0,30,101,30,30,1)
//SrwSmpCreate("ObsE2",30);SrwSmpScanXZE("ObsE2_obs",25,1,1,0,1,1,1,130,100)
//SrwSmpCreate("ObsX",30);SrwSmpScanXZE("ObsX_obs",0,150,300,0,1,1,39,39,1)
//SrwSmpCreate("ObsZ",30);SrwSmpScanXZE("ObsZ_obs",0,1,1,0,20,250,39,39,1)

// Compute SR
SrwStoWigCreate("SOLEIL_MedSectAsymWigObsE1","SOLEIL_MedSect_ebm","AsymWig_mag","ObsE1_obs",2)
//SrwStoWigCreate("ElecMagObsE1","Elec_ebm","Mag_mag","ObsE1_obs",1)
//SrwStoWigCreate("ElecMagObsE2","Elec_ebm","Mag_mag","ObsE2_obs",1)

// Extract data
SrwSto2IntF("SOLEIL_MedSectAsymWigObsE1_ras","I1",7,1,8,10.05,0,0,2)
//SrwSto2Int("ElecMagObsE1_ras","I",7,8,5.05,0,0,2)
//SrwSto2Int("ElecMagObsE2_ras","I",7,8,5.05,0,0,1)
//AppendToGraph ElecMagObsE2I_e
//ModifyGraph rgb(ElecMagObsE2I_e)=(0,0,0)
//Textbox/N=text0/F=0 " Spectral flux per unit surface \r \"Thick\" e-beam (ESRF Low Beta) \r Observation 30 m from wiggler \r   red: on-axis \r   black: off-axis (hor. pos. 25 mm) "

// Compute SR
SrwStoWigCreate("SOLEIL_MedSectAsymWigObsX","SOLEIL_MedSect_ebm","AsymWig_mag","ObsX_obs",2)
SrwSto2IntF("SOLEIL_MedSectAsymWigObsX_ras","I1",7,1,2,30,0,2.5e-09,2)
//SrwStoWigCreate("ElecMagObsZ","Elec_ebm","Mag_mag","ObsZ_obs",1)
// Extract data
//SrwSto2Int("ElecMagObsZ_ras","I",7,8,5.05,0,0,2)
//Textbox/N=text0/F=0 " Spectral flux per unit surface \r in vertical plane 30 m from wiggler \r Photon energy: 39 keV "

// Compute SR
SrwStoWigCreate("SOLEIL_MedSectAsymWigObsZ","SOLEIL_MedSect_ebm","AsymWig_mag","ObsZ_obs",2)
SrwSto2IntF("SOLEIL_MedSectAsymWigObsZ_ras","I1",7,1,3,30,2.5e-09,2.5e-09,2)
SrwSto2PolRateExt("SOLEIL_MedSectAsymWigObsZ_ras","Rcr",5,1,3,2,30,2.5e-09,2.5e-09,2)
// Extract data
//SrwSto2Int("ElecMagObsX_ras","I",7,8,5.05,0,0,2)
//Textbox/N=text0/F=0 " Spectral flux per unit surface \r in horizontal plane 30 m from wiggler\r Photon energy: 39 keV "

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
//srwUtiShowHelpTopic("Asymmetric Wiggler     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example WigAsym.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end