
//+++++++++++++++++++++++++++++++++++++++
//
// Example: UR Intensity Distribution from Thick E-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamURPowerDens()

srwUtiShowHelpTopic("UR Power Density     ")

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
//SrwMagPerCreate("U",period,1,k,length,4,0,0)
SrwMagPerCreate2D("U",period,k,0,length,0,4,0,0)

// Create Observation Structure
Variable Wx=20,Posx=0,npx=80;
Variable Wz=20,Posz=0,npz=80;
SrwSmpPowCreate("O",30); SrwSmpPowScanXZ("O_obp",Posx,Wx,npx,Posz,Wz,npz)

// Compute Power Density Distribution
Variable prec=1
SrwPowCreate("EUO","E_ebm","U_map","O_obp",prec,2,2)
Textbox/N=text0/F=0 " Power per unit surface \r at 30 m distance from the undulator. "

// Extract Horizontal and Vertical Intensity profiles
SrwPow2Int("EUO_pow","I",2,0,0,2)
Textbox/N=text0/F=0 " Power per unit surface \r at 30 m distance from the undulator. \r Vertical profile at zero hor. pos. "
SrwPow2Int("EUO_pow","I",1,0,0,2)
Textbox/N=text0/F=0 " Power per unit surface \r at 30 m distance from the undulator. \r Horizontal profile at zero vert. pos. "

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"
SrwUtiTriggerPrint(1)

TileWindows/O=1/C
srwUtiShowHelpTopic("UR Power Density     ")

print " "
print "The macro generating this computation can be found in the file \"SRW Example UR Power.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end