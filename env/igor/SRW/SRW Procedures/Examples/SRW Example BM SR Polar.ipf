
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Polarization of Bending Magnet Radiation
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamBMSRPolar()
Silent 1						|	Computing  ...

srwUtiShowHelpTopic("Polarization of Bending Magnet Synchrotron Radiation     ")

Print "Polarization of Bending Magnet Radiation"
variable t0=startmstimer
SrwInit(2)

// Create Electron Beam
SrwElecFilament("SOL_BM",2.75,0.5,0,0,0,0,0)

// Create Observation structures
SrwSmpCreate("ObsE",10);SrwSmpScanXZE("ObsE_obs",0,1,1,0,1,1,1,20,500)
SrwSmpCreate("ObsXZ",10);SrwSmpScanXZE("ObsXZ_obs",0,8,31,0,8,51,7,7,1)

// Create Magnetic Field
SrwMagFieldCreate("MagBM",0,1,1500)
SrwMagConst("MagBMBZ_fld",1.72)

// Compute Electric Field of emitted radiation vs photon energy
SrwMagPrec("MagBM_mag",3,0.01,0.01,10000,1,0,0)
SrwWfrCreate("Wfr","SOL_BM_ebm","MagBM_mag","ObsE_obs",1,1)

// Extract Intensity spectrum
SrwWfr2Int("Wfr_rad","Itot",7,1,1,1,7,0,0,2)
ModifyGraph lsize=1.5
SrwUtiGraphAddFrameAndGrid()
TextBox/C/N=text0/F=0/A=MC "Spectrum of Bending Magnet SR\rSpectral Flux per Unit Surface\rin Median Plane at 10 m from Source"
AppendText "E = 2.75 GeV, I = 0.5 A, B = 1.72 T"
TextBox/C/N=text0/X=2.45/Y=-23.2

// Compute Electric Field of emitted radiation on transverse grid
SrwMagPrec("MagBM_mag",3,0.01,0.01,10000,1,0,0)
SrwWfrCreate("Wfr","SOL_BM_ebm","MagBM_mag","ObsXZ_obs",1,1)

// Extract Intensities at Various Polarizations
SrwWfr2Int("Wfr_rad","Itot",7,1,4,1,7,0,0,2)
TextBox/C/N=text0/F=0/A=MC " Total Spectral Flux / Surface \r at 7 keV, at 10 m from Source "
AppendText " E = 2.75 GeV, I = 0.5 A, B = 1.72 T "
TextBox/C/N=text0/X=24.50/Y=40.00

SrwWfr2Int("Wfr_rad","Ilh",1,1,4,1,7,0,0,2)
TextBox/C/N=text0/F=0/A=MC " Spectral Flux / Surface at \r Linear Horizontal Polarization "
TextBox/C/N=text0/X=24.50/Y=40.00

SrwWfr2Int("Wfr_rad","Ilv",2,1,4,1,7,0,0,2)
TextBox/C/N=text0/F=0/A=MC " Spectral Flux / Surface at \r Linear Vertical Polarization "
TextBox/C/N=text0/X=24.50/Y=40.00

SrwWfr2Int("Wfr_rad","Icr",5,1,4,1,7,0,0,2)
TextBox/C/N=text0/F=0/A=MC " Spectral Flux / Surface at \r Circular Right Polarization "
TextBox/C/N=text0/X=24.50/Y=40.00

SrwWfr2Int("Wfr_rad","Icl",6,1,4,1,7,0,0,2)
TextBox/C/N=text0/F=0/A=MC " Spectral Flux / Surface at \r Circular Left Polarization "
TextBox/C/N=text0/X=24.50/Y=40.00

SrwWfr2Int("Wfr_rad","Itot",7,1,3,1,7,0,0,2)
SrwUtiGraphAddFrameAndGrid()
ModifyGraph rgb(WfrItot_z)=(0,0,0)
SrwWfr2Int("Wfr_rad","Ilh",1,1,3,1,7,0,0,1)
AppendToGraph WfrIlh_z
SrwWfr2Int("Wfr_rad","Ilv",2,1,3,1,7,0,0,1)
AppendToGraph WfrIlv_z
ModifyGraph rgb(WfrIlv_z)=(0,0,52224)
SrwWfr2Int("Wfr_rad","Icr",5,1,3,1,7,0,0,1)
AppendToGraph WfrIcr_z
ModifyGraph lstyle(WfrIcr_z)=3
SrwWfr2Int("Wfr_rad","Icl",6,1,3,1,7,0,0,1)
AppendToGraph WfrIcl_z
ModifyGraph rgb(WfrIcl_z)=(0,0,52224), lstyle(WfrIcl_z)=3, lsize=1.5
Legend/C/N=text0/J/F=0/A=MC "Spectral Flux / Surface\rat Different Polarizations\r\\s(WfrIlh_z) Linear Horizontal\r\\s(WfrIlv_z) Linear Vertical"
AppendText "\\s(WfrIcr_z) Circular Right\r\\s(WfrIcl_z) Circular Left\r\\s(WfrItot_z) Total"
TextBox/C/N=text0/X=27.45/Y=-12.00

// Extract Polarization Rates
SrwWfr2PolRateExt("Wfr_rad","Rlh",1,1,3,2,7,0,0,2)
SrwUtiGraphAddFrameAndGrid()
SrwWfr2PolRateExt("Wfr_rad","Rcr",5,1,3,2,7,0,0,1)
AppendToGraph WfrRcr_z
ModifyGraph lstyle(WfrRcr_z)=3, lsize=1.5
Legend/C/N=text0/J/F=0/A=MC "Polarization Rates (I1-I2)/(I1+I2)\r\\s(WfrRlh_z) Linear Horizontal\r\\s(WfrRcr_z) Circular Right"
TextBox/C/N=text0/X=27.45/Y=-30.00

TileWindows/O=1/C

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

srwUtiShowHelpTopic("Polarization of Bending Magnet Synchrotron Radiation     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example BM SR Polar.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
