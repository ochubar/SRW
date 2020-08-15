
//+++++++++++++++++++++++++++++++++++++++
//
//Focusing the central cone of the undulator radiation 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamImagUnd()

srwUtiShowHelpTopic("Focusing the central cone of Undulator Radiation     ")

Print "Focusing the Undulator Radiation. Be patient, It may take from a few 10s of seconds to a few minutes."

// Initialization
SrwInit(2)

variable t0=stopmstimer(0)
t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
SrwElecFilament("ESRF",6,0.2,-0.89,0,0,0,0)
SrwElecThick("ESRF_ebm",0.001,4,0.03,0.5,2.73,0,0,0,0)

// Create The Magnetic Field 
SrwMagFieldCreate("U35",0,1.8,4000)
SrwMagZero("U35BX_fld")
SrwMagSin("U35BZ_fld",1,0,35,0,0,46,0.7,0,0)

// Setup Mode of Longitudinal Integration
SrwMagPrec("U35_mag",1,0.0018,0.1,10000,1,0,0)

// Display Magnetic Field
SrwMagDisplayField("U35BZ_fld")

// Compute the Spectrum and plot it
SrwSmpCreate("O",10)
SrwSmpScanXZE("O_obs",0,0,1,0,0,1,0.1,9,600)
SrwWfrCreate("W","ESRF_ebm","U35_mag","O_obs",1,1)
SrwWfr2Int("W_rad","I",1,1,8,1,5,0,0,2)
Textbox/N=text1/F=0 " Spectral Flux / Surface on the axis "

// Create another observation range centered on the fundamental on axis
variable en=8.1
SrwSmpCreate("OF",10)
SrwSmpScanXZE("OF_obs",0,0.7/1.7,30,0,0.7/1.7,30,en,en,1)

SrwWfrCreate("Z","ESRF_ebm","U35_mag","OF_obs",2,1)
SrwWfr2Int("Z_rad","I",1,1,4,1,5,0,0,2)
Textbox/N=text1/F=0 " Spectral Flux / Surface \r Central Cone \r Harmonic : 3 \r 10 m from the Source "

// Create a Beamline 
SrwOptThinLens("ThinLens",5,5,0,0)
SrwOptDrift("Drift",8)
SrwOptCont("Container")
SrwOptContAdd("ThinLens_bli","Container_bli",2)
SrwOptContAdd("Drift_bli","Container_bli",2)
SrwOptDrift("D",2)

// Propagate the Wavefront through the Beamline
SrwWfrPropagate("Z_rad","Container_bli",1,2,"Zd")
SrwWfrResize("Zd_rad",2,1,4,1,4,2,"Zdd")
SrwWfr2Int("Zdd_rad","I_1",1,1,3,1,5,0,0,1)

SrwWfrPropagate("Zd_rad","D_bli",2,1,"")
SrwWfrResize("Zd_rad",2,1,4,1,4,2,"Zdd")
SrwWfr2Int("Zdd_rad","I0",1,1,3,1,5,0,0,1)

SrwWfrPropagate("Zd_rad","D_bli",2,1,"")
SrwWfrResize("Zd_rad",2,1,4,1,4,2,"Zdd")
SrwWfr2Int("Zdd_rad","I1",1,1,3,1,5,0,0,1)

Display ZddI_1_z,ZddI0_z,ZddI1_z
SetAxis bottom -3.0e-05,3.0e-05 
Tag/N=text0/F=0/A=RT ZddI0_z, 0," 10 m From Lens "
Tag/N=text1/F=0/A=RT ZddI1_z, -1e-5," 12 m From Lens "
Tag/N=text2/F=0/A=RT ZddI_1_z, 1e-5," 8 m From Lens "
Label bottom "Vertical Position "
Label left "Phot/sec/0.1%bw/mm^2"

SrwOptDrift("D",-2)
SrwWfrPropagate("Zd_rad","D_bli",2,1,"")
SrwWfrResize("Zd_rad",2,1,2,1,2,1,"")
SrwWfrResize("Zd_rad",1,3,1,3,1,1,"")
SrwWfr2Int("Zd_rad","I2",1,1,8,1,5,0,0,2)
Textbox/N=text1/F=0 " Spectral Flux / Surface \r After Lens \r 20 m from the Source \r 1:1 Imaging \r Filament Beam "

SrwWfr2Int("Zd_rad","I2m",1,2,4,1,5,0,0,2)
Textbox/N=text1/F=0 " Spectral Flux / Surface \r After Lens \r 20 m from the Source \r 1:1 Imaging \r Thick Beam "

TileWindows/O=1/C

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

srwUtiShowHelpTopic("Focusing the central cone of Undulator Radiation     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example UR Focusing.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
