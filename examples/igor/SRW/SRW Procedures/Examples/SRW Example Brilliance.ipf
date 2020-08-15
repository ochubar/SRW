
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Brilliance of UR source
//
//+++++++++++++++++++++++++++++++++++++++

proc SrwExamBrilUR()

srwUtiShowHelpTopic("Brilliance at ESRF     ")

// Initialization
SrwInit(2)
Variable t0=startmstimer
SrwUtiTriggerPrint(2)

// Load Electron Beam Library
SrwUtiElecLib()

// Magnetic Field
SrwMagPerCreate2D("U42",42,2.2,0,5,1.571,1,0,0)

// Compute Brilliance
SrwBrilUnd("BrU42ESRF","ESRF_HighBeta_ebm","U42_map",0.2,1,5,100,3,2)
    
ModifyGraph tick(bottom)=2,mirror(bottom)=1
ModifyGraph tick=2,mirror=1
SetAxis left 4.90066e+17,5e+20 
Textbox/N=text0/A=MC "Brilliance\rESRF\rBeamline : ID26\rU42\r Length : 5m\r Harmonic 1 to 5"

TileWindows/O=17/C

SrwUtiTriggerPrint(1)
Print "CPU Time : ",srround(stopmstimer(t0)/1e6,2)," seconds"

print " "
print "The macro generating this computation can be found in the file \"SRW Example Brilliance.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
