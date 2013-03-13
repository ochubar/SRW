
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Ellipsoidal Undulator
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamURElliptic()

srwUtiShowHelpTopic("Ellipsoidal Undulator     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=6, current=0.2
SrwElecFilament("E",energy,current,0,0,0,0,0)

Variable enspread=0.001
Variable emx=3.9, betax=35.6
Variable emz=0.039, betaz=2.5
SrwElecThick("E_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Undulator
Variable period=42, length=3.2, Kz=2, Kx=1
//SrwMagPerCreate("U",period,1,Kz,length,1,0,0)
//SrwMagPerAddHarm("U_map",1,2,Kx,1.571)
SrwMagPerCreate2D("U",period,Kz,Kx,length,1.571,1,0,0)

// Create Observation structures
Variable Wx=0.5,Posx=0
Variable Wz=0.5,Posz=0
Variable EnDepart=1,EnFinal=17,NpEn=2000
SrwSmpCreate("Oe",30); SrwSmpScanXZE("Oe_obs",Posx,Wx,1,Posz,Wz,1,EnDepart,EnFinal,NpEn)
SrwSmpCreate("Oxz",30);SrwSmpScanXZE("Oxz_obs",0,6,60,0,6,60,1,1,1)

// Compute Power Density
SrwPowCreate("EUOxz","E_ebm","U_map","Oxz_obs",1,2,2)
Textbox/A=RT/N=text0/F=0 "\\Z12 Power per unit surface \r 30 m from the undulator \r Kz = 2, Kx = 1, Phase = Pi/2 "

// Compute Spectrum
Variable HarmDepart=1,HarmFinal=7
Variable nlong=1,nangle=1
SrwPerStoCreate("EUOe","E_ebm","U_map","Oe_obs",HarmDepart,HarmFinal,nlong,nangle,1)

// Extract Intensity
SrwSto2Int("EUOe_ras","tot",7,8,5.05,0,0,2)
SrwSto2Int("EUOe_ras","cl",6,8,5.05,0,0,1)
SrwUtiGraphAddFrameAndGrid()
AppendToGraph EUOecl_e
ModifyGraph rgb('EUOecl_e')=(0,0,65280)
Textbox/A=RT/N=text0/F=0 "\\Z12 Spectral Flux observed through a 0.5 mm x 0.5 mm slit \r centered on the electron beam and placed \r at a distance of 30 m from the undulator \r Kz = 2, Kx = 1, Phase = Pi/2 \r   red curve: total flux \r   blue curve: circular-left polarized flux. "

// Extract Polarization Rate
SrwSto2PolRateExt("EUOe_ras","rcl",6,1,8,2,5.05,0,0,2)
SrwUtiGraphAddFrameAndGrid()
ModifyGraph rgb=(0,0,65280)
Textbox/A=RT/N=text0/F=0 "\\Z12 Circular (Left) Polarization Rate (Fcl - Fcr)/(Fcl + Fcr) \r calculated from flux at Kz = 2, Kx = 1, Phase = Pi/2. "
TextBox/C/N=text0/X=5.0/Y=20.0

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"
SrwUtiTriggerPrint(1)

srwUtiShowHelpTopic("Ellipsoidal Undulator     ")

TileWindows/O=1/C

print "The macro generating this computation can be found in the file \"SRW Example UR Ell.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
