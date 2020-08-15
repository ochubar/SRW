
//+++++++++++++++++++++++++++++++++++++++
//
//SASE Computation Example
//based on GENESIS 1.3
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamSASE_GENESIS()

srwUtiShowHelpTopic("SASE Wavefront     ")

//------------------------Initialization
SrwInit(2)
SrwUtiTriggerPrint(2) //Suppress printing of aux. SRW messages to the History window

//------------------------Electron Beam
Variable ElecEnergy = 25.0187 //[GeV]
Variable ElecPeakCurrent = 5000 //[A]
SrwElecFilament("E",ElecEnergy,ElecPeakCurrent,0,0,0,0,0)
SrwElecThick("E_ebm",8.e-05,0.02,0.02,20,20,0,0,0,0)

SrwElecLongCur("E_ebm",5000,0,1,0)

//------------------------Magnetic Field
Variable Period = 48.5 //Undulator Period [mm]
Variable UndK = 4.2128 //Deflecting Parameter
Variable NumOfSect = 1 //Number of Sections
Variable IntervalBwSect = 0 //Distance between Sections [m]
Variable NumOfPeriodsInSect = 3000 //Number of Periods per Section
Variable SectLen = NumOfPeriodsInSect*Period*0.001 //Section Length [m]
Variable NatFoc = 39.5 //Natural Focusing (exagerated value, to speed-up computation by avoiding quads)
SrwSASEUndCreate("M",1,UndK,Period,SectLen,NumOfSect,IntervalBwSect,NatFoc,NatFoc)
SrwSASEUndTaperAdd("M_mgg",1,0,0)
SrwSASEUndErrAdd("M_mgg",1,0)
SrwSASEFODOAdd("M_mgg",0.34,0,0,0.137,0.137,0,0)

//------------------------Input Radiation
Variable InRadPower = 7060. // Power per Unit Bandwidth dP/(dw/w),  [W]
Variable InRadWaist = 28. //Waist Diameter [microns]
Variable InRadWaistPos = -SectLen*NumOfSect - (NumOfSect - 1)*IntervalBwSect 
						//Position of the Waist at the beginning of the undulator [m]
SrwSASEInRadGsn("InR",InRadWaist,InRadWaistPos,InRadPower)

//------------------------Precision, Control
Variable IntegStepDivByPeriod = 25.//Longitudinal Integration Step div. by Und. Period
Variable NumMacroPart = 2^14 //Number of Macro-Particles
Variable PhotEnComp = 12.398 //Resonant Photon Energy value for which the computation is performed [keV]
Variable NumTransvMeshPts = 161 //Number of Transverse Mesh Points; this defines the number of radiation points vs horizontal and vertical positions in the wavefront at undulator exit
Variable MeshExtent = 10 //Transverse Mesh Extent Factor
//SrwSASEPrecMain("P",IntegStepDivByPeriod,NumOfMacroParticles,MeshExtent,NumOfTransvMeshPoints,40,0,0,1)
SrwSASEPrecGen("P",IntegStepDivByPeriod,NumMacroPart,1,1,PhotEnComp,MeshExtent,NumTransvMeshPts)
SrwSASEPrecExtra("P_pss",40,0,0,1)
SrwSASECntrl("C",2,1,2,1,1,2,2)

//------------------------SASE Computation
//SrwWfrSASECreate("Wfr","E_ebm","M_mgg","InR_rss","O_obs","P_pss","C_css",1,1)
SrwWfrSASEHarmCreate("Wfr","E_ebm","M_mgg","InR_rss","P_pss","C_css",5)

//Extracting intensity at 1st and 3rd harmonics:
SrwWfr2Int("WfrH1_rad","I",7,1,2,1,PhotEnComp,0,0,2) 
SrwWfr2Int("WfrH3_rad","I",7,1,2,1,PhotEnComp,0,0,1) 
TextBox/C/N=text1/F=0 "Intensity at the Exit of the Undulator"
AppendToGraph WfrH3I_x
ModifyGraph rgb(WfrH3I_x)=(0,0,0)
//Tag/C/N=text2/D=0.3/TL=0 WfrH1I_x, 0,"Harm. #1"
//Tag/C/N=text3/D=0.3/TL=0 WfrH3I_x, 0,"Harm. #3"
Tag/C/N=text2/D=0.3 WfrH1I_x, 0,"Harm. #1" //option /TL=0 doesn't work with Igor 6.0.0 (tested on Mac)
Tag/C/N=text3/D=0.3 WfrH3I_x, 0,"Harm. #3" //option /TL=0 doesn't work with Igor 6.0.0 (tested on Mac)

TileWindows/O=1/C

//------------------------Wavefront Propagation through free space
Variable LenDrift = 500 //[m]
SrwOptDrift("Drift",LenDrift)

SrwWfrResize("WfrH1_rad",1,1,0.3,1,0.3,1,"") //Reducing resolution to save memory for the propagation
SrwWfrPropagate("WfrH1_rad","Drift_bli",1,2,"WfrH1d")
SrwWfr2Int("WfrH1d_rad","I",7,1,2,1,PhotEnComp,0,0,2)
TextBox/C/N=text1/F=0 " Intensity at " +  num2str(LenDrift) + " m from \r the Exit of the Undulator "

SrwWfrResize("WfrH3_rad",1,1,0.7,1,0.7,1,"") //Reducing resolution to save memory for the propagation
SrwWfrPropagate("WfrH3_rad","Drift_bli",1,2,"WfrH3d")
SrwWfr2Int("WfrH3d_rad","I",7,1,2,1,PhotEnComp,0,0,1)
AppendToGraph WfrH3dI_x
ModifyGraph rgb(WfrH3dI_x)=(0,0,0)
//Tag/C/N=text2/D=0.3/TL=0 WfrH1dI_x, 0,"Harm. #1"
//Tag/C/N=text3/D=0.3/TL=0 WfrH3dI_x, 0,"Harm. #3"
Tag/C/N=text2/D=0.3/TL=0 WfrH1dI_x, 0,"Harm. #1" //option /TL=0 doesn't work with Igor 6.0.0 (tested on Mac)
Tag/C/N=text3/D=0.3/TL=0 WfrH3dI_x, 0,"Harm. #3" //option /TL=0 doesn't work with Igor 6.0.0 (tested on Mac)

SrwWfr2Int("WfrH1d_rad","I",7,1,4,1,PhotEnComp,0,0,2)
TextBox/C/N=text1/F=0 " Intensity at " +  num2str(LenDrift) + " m from \r the Exit of the Undulator (Harm. #1) "

TileWindows/O=1/C
srwUtiShowHelpTopic("SASE Wavefront     ")
SrwUtiTriggerPrint(1)

print " "
print "The macro generating this computation can be found in the file \"SRW Example SASE.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end
