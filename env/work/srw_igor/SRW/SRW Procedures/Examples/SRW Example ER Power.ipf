
//+++++++++++++++++++++++++++++++++++++++
//
// Example: ER Intensity Distribution from Thick E-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamERPowerDens()

srwUtiShowHelpTopic("ER Power Density     ")

// Initialization
SrwInit(2)

Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Create Electron Beam
Variable energy=2.5, current=0.2
SrwElecFilament("E",energy,current,0,0,0,0,0)

Variable enspread=0.001
Variable emx=3.9, betax=35.6
Variable emz=0.039, betaz=2.5
SrwElecThick("E_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Create Magnetic Field
SrwMagFieldCreate("M",0,6,1000)
SrwMagEdge("MBZ_fld",1,0,5,50,1.56)

SrwMagDisplayField("MBZ_fld")

// Create Observation Structure
Variable Wx=80,Posx=0,npx=60;
Variable Wz=10,Posz=0,npz=30;
SrwSmpPowCreate("O",20); SrwSmpPowScanXZ("O_obp",Posx,Wx,npx,Posz,Wz,npz)

// Compute Power Density Distribution
Variable prec=1
SrwPowCreate("EMO","E_ebm","M_mag","O_obp",prec,1,2)
Textbox/N=text0/F=0 " Power per unit surface at 20 m from \r the middle of straight section."

// Extract Horizontal and Vertical Intensity profiles
SrwPow2Int("EMO_pow","I",1,0,0,2)

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("ER Power Density     ")
SrwUtiTriggerPrint(1)

print "The macro generating this computation can be found in the file \"SRW Example ER Power.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
 
 
 
  