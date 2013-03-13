
//+++++++++++++++++++++++++++++++++++++++
//
// Example: UR Spectrum from Thick E-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamURSpecThickEbeam()

srwUtiShowHelpTopic("UR Spectrum through a Slit     ")

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
Variable period=35, length=3.2, k=2.2
//SrwMagPerCreate("U",period,1,k,length,1,0,0)
SrwMagPerCreate2D("U",period,k,0,length,0,1,0,0)

// Create Observation Structure
Variable Wx=0.5,Posx=0
Variable Wz=0.5,Posz=0
Variable EnDepart=1,EnFinal=21,NpEn=2000
SrwSmpCreate("O",30); SrwSmpScanXZE("O_obs",Posx,Wx,1,Posz,Wz,1,EnDepart,EnFinal,NpEn)
 
// Compute Spectrum
Variable HarmDepart=1,HarmFinal=7
Variable nlong=1,nangle=1
SrwPerStoCreate("EUO","E_ebm","U_map","O_obs",HarmDepart,HarmFinal,nlong,nangle,1)

// Extract Intensity
SrwSto2Int("EUO_ras","I",7,8,5.05,0,0,2)
Textbox/A=RT/N=text0/F=0 "Spectrum through a 0.5 mm x 0.5 mm aperture\ron the axis of the electron beam\rat 30 m distance from the undulator"

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"
SrwUtiTriggerPrint(1)

// Error Testing
if(abs(EUOI_e(17028.)/(6.3395e13) - 1.) > 5e-2)
	Print "Error executing Test"
	abort
else
	Print "Test Successful"
endif
Print ""

TileWindows/O=17/C

print "The macro generating this computation can be found in the file \"SRW Example UR SpecThkE.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
