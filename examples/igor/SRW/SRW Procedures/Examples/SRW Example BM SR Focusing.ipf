
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Focusing of visible light from a Bending Magnet 
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamImagBm()

srwUtiShowHelpTopic("Focusing the Bending Magnet Radiation     ")
Print "Focusing the Bending Magnet Radiation. Be patient, It may take from a few 10s of seconds to a few minutes."

// Initialization
SrwInit(2)

variable t0=stopmstimer(0)
t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
SrwElecFilament("E",6,0.2,0,0,0,0,0)
SrwElecThick("E_ebm",0.001,4,0.1,3,35,1.77,-0.022,0.0349,-0.0262)

// Create Magnetic Field
SrwMagFieldCreate("BM",0,3,1000)
SrwMagConst("BMBZ_fld",0.85)
SrwMagConst("BMBX_fld",0.0)

// Setup Mode of Longitudinal Integration
SrwMagPrec("BM_mag",3,0.01,0.01,3000,1,0,0)

// Display the Trajectory
SrwMagElecTraj("E_ebm","BM_mag",1,2,1,1)
Textbox/N=text0/F=0 " Electron Trajectory \r Energy : 6 GeV \r Field : 0.85 Tesla "

// Create the range of Observation 
SrwSmpCreate("O",20)
SrwSmpScanXZE("O_obs",0,10,2,0,30,2,0.003,0.003,1)

// Compute the SR
SrwWfrCreate("W","E_ebm","BM_mag","O_obs",2,1)

// Visualize Intensity
SrwWfr2Int("W_rad","I",1,1,8,1,5,0,0,2)
Textbox/N=text0/F=0 " Spectral Flux / Surface Horiz. Polar. \r Photon Energy : 3 eV  \r Dist. from Bending Magnet : 20 m"

// Create a Beamline with optical magnification ratio R
Variable R, F, D
R=0.1
F=20/(1+1/R)
D=20*F/(20-F)

SrwOptThinLens("ThinLens",F,F,0,-0.05)
SrwOptDrift("Drift",D)
SrwOptCont("Container")
SrwOptContAdd("ThinLens_bli","Container_bli",2)
SrwOptContAdd("Drift_bli","Container_bli",2)

// Resize the SR Wavefront in order to chek/increase the accuracy
SrwWfrResize("W_rad",1,1,1.3,1,1.3,2,"Wd")
//SrwWfrDupl("W_rad","Wd")

// Propagate the Wavefront through the Beamline
SrwWfrPropagate("Wd_rad","Container_bli",1,1,"")

// Reduce the range (for better viewing)
SrwWfrResize("Wd_rad",1,0.3,1,0.6,1,1,"Wd")

// Extract and Plot the results 
SrwWfr2Int("Wd_rad","I",1,1,8,1,5,0,0,2)
Textbox/N=text0/F=0 " Spectral Flux / Surface Horiz. Polar. \r In conjugate plane of 10:1 imaging \r Filament Beam "

SrwWfr2Int("Wd_rad","J",2,1,8,1,5,0,0,2)
Textbox/N=text0/F=0 " Spectral Flux / Surface Vert. Polar. \r In conjugate plane of 10:1 imaging \r Filament Beam "

SrwWfr2Int("Wd_rad","J",2,1,3,1,5,0,0,2)
SrwWfr2Int("Wd_rad","Jm",2,2,3,1,5,0,0,1)
AppendToGraph WdJm_z
ModifyGraph rgb(WdJm_z)=(0,12800,52224)
Textbox/A=RT/N=text0/F=0 " Spectral Flux / Surface \r Vertically Polarized \r In the Image Plane \r  red- Filament Beam \r  blue- Thick Beam"

TileWindows/O=1/C

srwUtiShowHelpTopic("Focusing the Bending Magnet Radiation     ")

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example BM SR Focusing.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher) "

end
