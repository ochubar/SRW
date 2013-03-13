
//+++++++++++++++++++++++++++++++++++++++
//
// Example: Coherent Synchrotron Radiation from Bending Magnet
//
//+++++++++++++++++++++++++++++++++++++++
proc SrwExamCSR_RotBunch()

srwUtiShowHelpTopic("Coherent Synchrotron Radiation from Rotated Electron Bunches     ")

variable t0=startmstimer

//Initialization
SrwInit(2)

//Create Electron Beam:
//     - zero-order and first-order statistical moments:
variable E = 2.75 //electron energy [GeV]
variable Iavg = 0.02 //average electron current [A]
SrwElecFilament("ElecR",E,Iavg,0,0,0,0,0)

//     - second-order transverse statistical moments and energy spread:
variable SigmaX = 43 //horizontal RMS bunch size [microns]
variable SigmaZ = 2000 //vertical RMS bunch size [microns], assuming "rotated" bunch
variable SigmaXp = 106 //horizontal RMS bunch divergence [micro-rad]
variable SigmaZp = 204 //vertical RMS bunch divergence [micro-rad]
variable MXXp = -0.89 //horizontal mixed moment [nm*rad]
variable MZZp = -332 //vertical mixed moment [nm*rad]
variable SigmaEdE = 0.001016 //RMS relative energy spread
SrwElecThickMom("ElecR_ebm",SigmaEdE,SigmaX,SigmaZ,SigmaXp,SigmaZp,MXXp,MZZp)

//     - longitudinal bunch parameters:
variable NumElec = 7.3 // x 1e+09 number of electrons in one bunch
variable SigmaS = 0.03 //RMS bunch length [mm]
SrwElecLong("ElecR_ebm",NumElec,0,SigmaS)

//Create Magnetic Field
variable Bmax = 1.72 //constant magnetic field in bending magnet
variable RangeB = 1.5 //magnetic field definition range [m]
SrwMagFieldCreate("BM",0,RangeB,1000)
SrwMagConst("BMBZ_fld",Bmax)

//Define Radiation Sampling
variable ObsDist = 3 //observation distance [m]

//     - scan vs photon energy
variable WavelengthMax_mic = 1000 //max. wavelength [microns]
variable WavelengthMinCSR_mic = 30 //min. wavelength for CSR calculation [microns]
variable WavelengthMinISR_mic = 10 //min. wavelength for incoherent SR calculation [microns]
variable NpE = 100 //number of point vs wavelength / photon energy
SrwSmpCreate("ObsECSR",ObsDist)
SrwSmpScanXZE("ObsECSR_obs",0,0,1,0,0,1,(12.3985e-4)/WavelengthMax_mic,(12.3985e-4)/WavelengthMinCSR_mic,NpE)
SrwSmpCreate("ObsEISR",ObsDist)
SrwSmpScanXZE("ObsEISR_obs",0,0,1,0,0,1,(12.3985e-4)/WavelengthMax_mic,(12.3985e-4)/WavelengthMinISR_mic,NpE)

//     - scan vs vertical position
variable WavelengthCen_mic = 100 //constant wavelength for intensity profile calculation [microns]
variable RangeZ = 180 //range of vertical position of the observation point [mm]
variable NpZ = 61 //number of point vs vertical position
SrwSmpCreate("ObsZ",ObsDist); SrwSmpScanXZE("ObsZ_obs",0,1,1,0,RangeZ,NpZ,(12.3985e-4)/WavelengthCen_mic,(12.3985e-4)/WavelengthCen_mic,1)

//Compute CSR and incoherent SR
variable StartIntegS = -0.5 //initial point of longitudinal integration [m]
variable FinIntegS = 0.5 //initial point of longitudinal integration [m]
variable StepIntegS = 0.001 //step of longitudinal integration [m]

//     - calculate CSR electric field vs photon energy
SrwWfrCSRCreate("ElecRBMObsECSR","ElecR_ebm","BM_mag","ObsECSR_obs",1,StepIntegS,StartIntegS,FinIntegS,1,1)
//     - extract CSR intensity spectrum vs photon energy
SrwWfr2Int("ElecRBMObsECSR_rad","I",7,1,1,1,1,0,0,2)
SrwUtiGraphAddFrameAndGrid()
ModifyGraph log(left)=1

//     - calculate single-electron SR electric field vs photon energy
SrwMagPrec("BM_mag",3,0.01,0.01,10000,1,0,0)
SrwWfrCreate("ElecRBMObsEISR","ElecR_ebm","BM_mag","ObsEISR_obs",1,1)
//     - extract incoherent SR intensity spectrum vs photon energy
SrwWfr2Int("ElecRBMObsEISR_rad","I",7,1,1,1,1,0,0,1)
AppendToGraph 'ElecRBMObsEISRI_e'
ModifyGraph rgb('ElecRBMObsEISRI_e')=(0,0,0)
ModifyGraph lsize=1.2
Legend/C/N=text0/J/X=29.4/Y=33.7/A=MC " Intensity Spectrum in the Median Plane \r\\s(ElecRBMObsECSRI_e) CSR \r\\s(ElecRBMObsEISRI_e) Incoherent SR "
AppendText " Number of Electrons per Bunch: " + num2str(NumElec*1e+9) + " \r Average Current: " + num2str(Iavg) + " A \r Distance from Source: " + num2str(ObsDist) + " m "

//     - calculate incoherent SR electric field vs vertical position
SrwMagPrec("BM_mag",3,0.01,0.01,10000,1,0,0);SrwWfrCreate("ElecRBMObsZ","ElecR_ebm","BM_mag","ObsZ_obs",1,1)
//     - extract incoherent SR intensity vs vertical position
SrwWfr2Int("ElecRBMObsZ_rad","ISRtot",7,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,2) //total
ModifyGraph rgb('ElecRBMObsZISRtot_z')=(0,0,0)
SrwWfr2Int("ElecRBMObsZ_rad","ISRh",1,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,1) //linear horizontal polarization
AppendToGraph 'ElecRBMObsZISRh_z'
SrwWfr2Int("ElecRBMObsZ_rad","ISRv",2,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,1) //linear vertical polarization
AppendToGraph 'ElecRBMObsZISRv_z'
ModifyGraph rgb('ElecRBMObsZISRv_z')=(0,0,65280)
ModifyGraph lsize=1.2
Legend/C/N=text0/J/A=MC " Incoherent SR Intensity vs Vertical Position \r at " + num2str(WavelengthCen_mic) + " \\F'Symbol'm\\F'Times New Roman'm Wavelength; Distance from Source: " + num2str(ObsDist) + " m \r\\s(ElecRBMObsZISRh_z) Linear Horizontal Polarization "
AppendText "\\s(ElecRBMObsZISRv_z) Linear Vertical Polarization \r\\s(ElecRBMObsZISRtot_z) Total "
SrwUtiGraphAddFrameAndGrid()

//     - calculate CSR electric field vs vertical position
SrwWfrCSRCreate("ElecRBMObsZ","ElecR_ebm","BM_mag","ObsZ_obs",1,StepIntegS,StartIntegS,FinIntegS,1,1)
//     - extract CSR intensity vs vertical position
SrwWfr2Int("ElecRBMObsZ_rad","CSRtot",7,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,2) //total
ModifyGraph rgb('ElecRBMObsZCSRtot_z')=(0,0,0)
SrwWfr2Int("ElecRBMObsZ_rad","CSRh",1,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,1) //linear horizontal polarization
AppendToGraph 'ElecRBMObsZCSRh_z'
SrwWfr2Int("ElecRBMObsZ_rad","CSRv",2,1,3,1,(12.3985e-4)/WavelengthCen_mic,0,0,1) //linear vertical polarization
AppendToGraph 'ElecRBMObsZCSRv_z'
ModifyGraph rgb('ElecRBMObsZCSRv_z')=(0,0,65280)
ModifyGraph lsize=1.2
Legend/C/N=text0/J/A=MC " CSR Intensity vs Vertical Position at " + num2str(WavelengthCen_mic) + " \\F'Symbol'm\\F'Times New Roman'm Wavelength \r Distance from Source: " + num2str(ObsDist) + " m \r\\s(ElecRBMObsZCSRh_z) Linear Horizontal Polarization "
AppendText "\\s(ElecRBMObsZCSRv_z) Linear Vertical Polarization \r\\s(ElecRBMObsZCSRtot_z) Total "
AppendText " RMS Sizes of (rotated) Electron Bunch: \r    Longitud.: " + num2str(SigmaS) + " mm; Horiz.: " + num2str(SigmaX*0.001) + " mm; Vert.: " + num2str(SigmaZ*0.001) + " mm "
SrwUtiGraphAddFrameAndGrid()

print "CPU Time : ",srround(stopmstimer(t0)/1e6,1)," seconds"

TileWindows/O=1/C
srwUtiShowHelpTopic("Coherent Synchrotron Radiation from Rotated Electron Bunches     ")

print "The macro generating this computation can be found in the file \"SRW Example CSR.ipf\""
print "that can be accessed through the menu \"Windows->Procedure Windows\" (assuming IGOR Pro version 6 or higher)"

end