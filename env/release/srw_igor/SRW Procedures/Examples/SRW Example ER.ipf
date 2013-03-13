
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Edge Radiation at SOLEIL
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamER(v)
variable v=1
Silent 1						|	Computing  ...
PauseUpdate

srwUtiShowHelpTopic("Near-Field Edge Radiation     ")

Print "Bending Magnet Edge Radiation at SOLEIL"
variable t0=startmstimer
SrwInit(2)

// Create Electron Beam
SrwElecFilament("SOL",2.5,0.5,0,0,0,0,0)

// Create Magnetic Field
SrwMagFieldCreate("Edge",0,6,1000)
SrwMagEdge("EdgeBZ_fld",1,0,5,50,1.56)

SrwMagDisplayField("EdgeBZ_fld")

// Display Trajectory
SrwMagElecTraj("SOL_ebm","Edge_mag",2,2,1,1)

// Setup Longitudinal Integration parameters
SrwMagPrec("Edge_mag",3,0.01,0.1,10000,1,-3,3)

// Create Observation
SrwSmpCreate("O3",3.2)
SrwSmpScanXZE("O3_obs",0,20,200,0,0,1,0.0001,9,1)

// Compute
SrwWfrCreate_("SOLEdgeO3","SOL_ebm","Edge_mag","O3_obs")

// Visualize Intensity
SrwWfr2Int_("SOLEdgeO3_rad","I",1,1,8,5,0,0,2)

TileWindows/O=1/C

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

srwUtiShowHelpTopic("Near-Field Edge Radiation     ")

// Error Testing
if (abs(SOLEdgeO3I_x(-3e-3)/1.38582e12 - 1) > 5e-2)
print "Error"
abort
else
print "Test Successful"
print""
endif

print "The macro generating this computation can be found in the file \"SRW Example ER.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
