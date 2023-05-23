
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Off-Axis UR
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamUR()

srwUtiShowHelpTopic("Off-Axis Undulator Radiation     ")

variable t0=startmstimer

// Initialization
SrwInit(2)

// Create Electron Beam
SrwElecFilament("Elec",6,0.2,0,0,0,0,0)

// Create Magnetic Field
SrwMagFieldCreate("U35",0,1.7,2000)
SrwMagSin("U35BZ_fld",1,0,35,0,0,46.5,0.7,0,0)

// Set Mode of Longitudinal Integration
SrwMagPrec("U35_mag",2,0.01,0.01,3000,1,0,0)

// Setup Radiation Sampling
SrwSmpCreate("Obs",30); SrwSmpScanXZE("Obs_obs",0,0,1,0,0,1,11,14,400)
// Compute Electric Field
SrwWfrCreate_("ElecU35Obs","Elec_ebm","U35_mag","Obs_obs")
// Extract Intensity
SrwWfr2Int_("ElecU35Obs_rad","I",1,1,8,5.05,0,0,1)

// Setup Radiation Sampling
SrwSmpCreate("Obs1",30); SrwSmpScanXZE("Obs1_obs",0,0,1,1,0,1,11,14,400)
// Compute Electric Field
SrwWfrCreate_("ElecU35Obs1","Elec_ebm","U35_mag","Obs1_obs")
// Extract Intensity
SrwWfr2Int_("ElecU35Obs1_rad","I",1,1,8,5.05,0,0,1)

// Setup Radiation Sampling
SrwSmpCreate("Obs2",30); SrwSmpScanXZE("Obs2_obs",0,0,1,1.5,0,1,11,14,400)
// Compute Electric Field
SrwWfrCreate_("ElecU35Obs2","Elec_ebm","U35_mag","Obs2_obs")
// Extract Intensity
SrwWfr2Int_("ElecU35Obs2_rad","I",1,1,8,5.05,0,0,1)

// Setup Radiation Sampling
SrwSmpCreate("Obs3",30); SrwSmpScanXZE("Obs3_obs",0,0,1,2,0,1,11,14,400)
// Compute Electric Field
SrwWfrCreate_("ElecU35Obs3","Elec_ebm","U35_mag","Obs3_obs")
// Extract Intensity
SrwWfr2Int_("ElecU35Obs3_rad","I",1,1,8,5.05,0,0,1)

// Display Magnetic Field
SrwMagDisplayField("U35BZ_fld")

// Display Trajectory
SrwMagElecTraj("Elec_ebm","U35_mag",2,1,1,1)

// Display Extracted Intensity Spectra
Display /W=(78,93,534,382) ElecU35ObsI_e,ElecU35Obs1I_e,ElecU35Obs2I_e,ElecU35Obs3I_e
ModifyGraph lStyle(ElecU35Obs1I_e)=1,lStyle(ElecU35Obs2I_e)=2,lStyle(ElecU35Obs3I_e)=4
ModifyGraph rgb=(0,0,0)
ModifyGraph tick=2
ModifyGraph mirror=1
ModifyGraph lblMargin(left)=13
ModifyGraph lblLatPos(left)=2
Label left "Photons/sec/0.1%bw/mm^2"
Label bottom "Photon Energy"
SetAxis left -31300900000000,1.26788e+15
Textbox/N=text0/A=MC/X=-20.71/Y=28.63 "Electron Beam : ESRF\rUndulator : 46.5 periods of 35 mm B = 0.7 T\rHarmonic : 5"
AppendText "\\s(ElecU35ObsI_e) : On axis\r\\s(ElecU35Obs1I_e) : Displaced 1 mm Vertically"
AppendText "\\s(ElecU35Obs2I_e) : Displaced 1.5 mm Vertically\r\\s(ElecU35Obs3I_e) : Displaced 2 mm Vertically"

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

// Error Testing
if(abs(ElecU35ObsI_e(13503.8)/(1.23575e15) - 1.) > 5e-2)
	print "Error"
	abort
else
	print "Test Successful"
	print""
	
	TileWindows/O=1/C
	srwUtiShowHelpTopic("Off-Axis Undulator Radiation     ")
endif

print "The macro generating this computation can be found in the file \"SRW Example UR Off-Axis.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end