
//+++++++++++++++++++++++++++++++++++++++
//
// Example: UR Intensity Distribution from Thick E-beam
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamURIntDistrThickEbeam()

srwUtiShowHelpTopic("UR Angular Pattern     ")

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
Variable Wx=6,Posx=0,npx=80;
Variable Wz=6,Posz=0,npz=80;
Variable en=13.5,NpEn=1;
SrwSmpCreate("O",30); SrwSmpScanXZE("O_obs",Posx,Wx,npx,Posz,Wz,npz,en,en,NpEn)
  
// Compute Intensity Distribution
Variable HarmDepart=3,HarmFinal=7;
Variable nlong=1,nangle=1
SrwPerStoCreate("EUO","E_ebm","U_map","O_obs",HarmDepart,HarmFinal,nlong,nangle,1)

// Extract Intensity
SrwSto2Int("EUO_ras","I1",7,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r at a distance of 30 m from an ESRF undulator \r for a photon energy of 13.5 keV \r taking into account 4 (0.04) nm emittance \r in the horizontal (vertical) plane."

// Reduce emittances
emx=0.039; emz=0.00039
enspread=0.00001
SrwElecThick("E_ebm",enspread,emx,emz,betax,betaz,0,0,0,0)

// Compute Intensity Distribution
SrwPerStoCreate("EUO","E_ebm","U_map","O_obs",HarmDepart,HarmFinal,nlong,nangle,1)

// Extract Intensity
SrwSto2Int("EUO_ras","I2",7,8,5.05,0,0,2)
Textbox/N=text0/F=0 " Spectral flux per unit surface \r at a distance of 30 m from an ESRF undulator \r for a photon energy of 13.5 keV \r assuming a filament electron beam."

Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("UR Intensity Distribution from Thick Electron Beam     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example UR AngPatThkE.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
 
 
 
  