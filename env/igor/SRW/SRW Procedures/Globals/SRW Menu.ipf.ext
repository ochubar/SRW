
//+++++++++++++++++++++++++++++++++++++++
//
// Menu SRW:Emission
//
//+++++++++++++++++++++++++++++++++++++++
Menu "SRWE"

SubMenu "Help"
"Introduction",srwUtiShowHelpTopic("Introduction     ")
"Arbitrary SR Source",srwUtiShowHelpTopic("Near Field Computation     ")
"Undulator",srwUtiShowHelpTopic("Fast Computation of Undulator Radiation Spectra through a Slit     ")
"Wiggler",srwUtiShowHelpTopic("Computation of Wiggler Radiation     ")
"Bending Magnet",srwUtiShowHelpTopic("Computation of Bending Magnet Radiation     ")
"Gaussian Beam",srwUtiShowHelpTopic("Gaussian Beam     ")
"SASE FEL",srwUtiShowHelpTopic("SASE FEL: Wavefront Amplification     ")
"Power Density",srwUtiShowHelpTopic("Computation of SR Power Density     ")
"-"
"Run All SRW Examples",SrwUtiAllExam()
"FAQ",srwUtiShowHelpTopic("Frequently Asked Questions     ")
"About SRW",srwUtiShowHelpTopic("About SRW     ")
end

SubMenu "Utilities"
"Initialize",SrwInit(1)
"View | Modify SRW Structures...",SrwUtiViewStructDialog()
"Printing to History...",SrwUtiTriggerPrint()
"-"
"Graph: Add Frame and Grid", SrwUtiGraphAddFrameAndGrid()
"Graph: Add Axes Labels...", SrwUtiGraphAddLabels()
"Graph: Resize Window...", SrwUtiGraphWindResize()
"Graph: Save to File...", SrwUtiGraphSaveFile()
"Graphs: Close All",SrwUtiKillAllGraphs()
"-"
"Conversion: Photon Energy...",SrwUtiPhotEnConv()
"Conversion: Spectral Flux...",SrwUtiSpecFluxConv()
"-"
//"Electron Beam: Load Library",SrwUtiElecLib()
//"-"
"Bending Magnet: Angle from Field Integral...",SrwUtiMagFieldInt2Ang()
"Bending Magnet: Radius...",SrwUtiMagRad()
"Bending Magnet: Crit. Photon Energy...",SrwUtiMagCritPhotEn()
"-"
"Undulator: Deflection Parameter [K]...",SrwUtiUndK()
"Undulator: Fund. Photon Energy from K...",SrwUtiUndFundPhotEnFromK()
"Undulator: Fund. Photon Energy from B...",SrwUtiUndFundPhotEn()
"Undulator: PPM...",SrwMagPerPPM()
"-"
"Load Electron Beam Library",SrwUtiElecLib()
//"Load SR Source Data...",SrwUtiSrcSelect()
end

"-"
//"Electron Beam...",SrwElecDialog(0)
SubMenu "Electron Beam"
"Main Parameters...",SrwElecDialog(0)
"-"
"Longitudinal Dependence...",SrwElecLongCur()
"Cross-Moments: Transverse...",SrwElecThickMomCrossTransv()
"Cross-Moments: Long. Position && Energy...",SrwElecThickMomCrossEnLong()
"-"
"Container...",SrwElecCont()
"Add to Container...",SrwElecContAdd()
end
"Radiation Sampling...",SrwRadSamplDialog("") //,SrwRadSamplDialog(SrwSmpType)

"-"
SubMenu "Arbitrary SR Source"
SubMenu "Help"
"Near Field",srwUtiShowHelpTopic("Near Field Computation     ")
"Power Density",srwUtiShowHelpTopic("Computation of SR Power Density     ")
"-"
"Example: Off-Axis UR",SrwExamUR()
"Example: SR Polarization",SrwExamBMSRPolar()
"Example: ER Spectral Flux",SrwExamER(1)
"Example: ER Power Density",SrwExamERPowerDens()
"Example: CSR from Rotated Bunch",SrwExamCSR_RotBunch()
end
"-"
SubMenu "Magnetic Field"
"Create and Modify...",SrwMagFieldCreate()
"-"
"Zero...",SrwMagZero()
"Constant...",SrwMagConst()
"Sinusoidal...",SrwMagSin()
"BM Edges...",SrwMagEdge()
"Gaussian Angle...",SrwMagGsnAng()
"Dipole Magnet...",SrwMagDipole()
"-"
"Import Component...",SrwMagAddImpCmpn() 
//SrwMagImportCmpn()
"-"
"Display...",SrwMagDisplayField()
"Duplicate...",SrwMagDupl()
"Convert to Periodic...",SrwMagArb2Per()
"-"
//"Create Trajectory...",SrwTrjCreate()
"Create Trajectory...",SrwTrjCreateTransvUnif()
//"Trajectory: Display ...",SrwMagElecTraj()
end
"-"
//"Compute Electric Field...",SrwWfrCreateDialog(1)
"Compute Electric Field: Single-Electron SR...",SrwWfrCreateDialog(1)
"Compute Electric Field: CSR [test]...",SrwWfrCSRCreate()
//"Compute Stokes...",SrwStoArbCreate()
//"Compute Power Density...",SrwPowCreate()
"Compute Power Density...",SrwPowDenCreate()
end  // Arbitrary Source

SubMenu "Undulator"
SubMenu "Help"
"Spectrum through a Slit",srwUtiShowHelpTopic("Fast Computation of Undulator Radiation Spectra through a Slit     ")
"Power Density",srwUtiShowHelpTopic("Computation of SR Power Density     ")
"-"
"Example: Spectrum through a Slit",SrwExamURSpecThickEbeam()
"Example: Angular Pattern",SrwExamURIntDistrThickEbeam()
"Example: Power Density",SrwExamURPowerDens()
"Example: Ellipsoidal",SrwExamURElliptic()
"Example: Brilliance",SrwExamBrilUR()
end
"-"
SubMenu "Periodic Magnetic Field"
"Create and Modify...",SrwMagPerCreateDialog()
"Add Harmonic...",SrwMagPerAddHarm()
"-"
"Duplicate...",SrwMagPerDupl()
end
"-"
//"Compute Stokes...",SrwPerStoCreate()
"Compute Stokes...",SrwPerStokesCreate()
"Compute Power Density...",SrwPowCreate()
SubMenu "Estimate Brilliance"
"Following \"X-Ray Data Booklet\"...",SrwBrilUnd()
"Treating Energy Spread and \"Detuning\"...",SrwBrilUndEnDet()
end
end  // Undulator Radiation

SubMenu "Wiggler"
SubMenu "Help"
"Spectral Angular Distribution",srwUtiShowHelpTopic("Computation of Wiggler Radiation     ")
"Power Density",srwUtiShowHelpTopic("Computation of SR Power Density     ")
"-"
"Example: Planar Wiggler",SrwExamWigPlanar()
"Example: Ellipsoidal Wiggler",SrwExamWigElliptic()
end
"-"
SubMenu "Magnetic Field"
"Create and Modify...",SrwMagFieldCreate()
"-"
"Zero...",SrwMagZero()
"Constant...",SrwMagConst()
"Sinusoidal...",SrwMagSin()
"BM Edges...",SrwMagEdge()
"Gaussian Angle...",SrwMagGsnAng()
"-"
"Import Component...",SrwMagAddImpCmpn() 
//SrwMagImportCmpn()
"-"
"Display...",SrwMagDisplayField()
"Duplicate...",SrwMagDupl()
"Convert to Periodic...",SrwMagArb2Per()
"-"
"Create Trajectory...",SrwTrjCreateTransvUnif()
//"Trajectory: Display ...",SrwMagElecTraj()
end
SubMenu "Periodic Magnetic Field"
"Create and Modify...",SrwMagPerCreateDialog()
"Add Harmonic...",SrwMagPerAddHarm()
"PPM Undulator...",SrwMagPerPPM()
"-"
"Duplicate...",SrwMagPerDupl()
end
"-"
"Compute Stokes...",SrwStoWigglerCreate()
"Compute Power Density...",SrwPowCreate()
"Estimate Brilliance...",SrwBrilWig()
end  // Wiggler

SubMenu "Bending Magnet"
SubMenu "Help"
"Spectral Angular Distribution",srwUtiShowHelpTopic("Computation of Bending Magnet Radiation     ")
"Power Density",srwUtiShowHelpTopic("Computation of SR Power Density     ")
"-"
"Example: ESRF Bending Magnet",SrwExamStdBM()
"Example: SR Polarization",SrwExamBMSRPolar()
end
"-"
"Constant Magnetic Field...",SrwMagConstCreate()
"-"
"Compute Stokes...",SrwStoConstCreate()
"Compute Power Density...",SrwPowCreate()
"Estimate Brilliance...",SrwBrilBM()
end  // Bending Magnet

SubMenu "Gaussian Beam"
SubMenu "Help"
"Basics...",srwUtiShowHelpTopic("Gaussian Beam     ")
"-"
"Example: Definition and Propagation",SrwExamGsnBm()
end
"-"
//"TEM Mode...",SrwGsnBeam()
//"Create...",SrwGsnBeamCreate()
"Main Parameters...",SrwGsnBeamPulsedLaser()
"Position and Angle...",SrwGsnBeamMom1()
//"Pulse Duration...",SrwGsnBeamTime()
"-"
"Compute Electric Field...",SrwWfrGsnBeamCreate()
end

SubMenu "SASE FEL"
SubMenu "Help"
"Wavefront Amplification...",srwUtiShowHelpTopic("SASE FEL: Wavefront Amplification     ")
"-"
"Example: SASE Wavefront",SrwExamSASE_GENESIS()
end
"-"
"Undulator: General...",SrwSASEUndCreate()
"Undulator: Taper...",SrwSASEUndTaperAdd()
"Undulator: Field Errors...",SrwSASEUndErrAdd()
"FODO...",SrwSASEFODOAdd()
//"Electron Beam Extras...",SrwElecSASEExtra()
"-"
"Input Radiation...",SrwSASEInRadGsn()
"-"
"Precision: General...",SrwSASEPrecGen()
"Precision: Time...",SrwSASEPrecTime()
"Precision: Extra...",SrwSASEPrecExtra()
"Control...",SrwSASECntrl()
"-"
"Compute SASE...",SrwWfrSASEHarmCreate()
end
"-"
"Visualize...",SrwVisualizeDialog()
"Duplicate...",SrwRadDupl()

end  // Menu SRWE

//+++++++++++++++++++++++++++++++++++++++
//
// Menu SRW:Propagation
//
//+++++++++++++++++++++++++++++++++++++++
Menu "SRWP"

SubMenu "Help"
"Introduction",srwUtiShowHelpTopic("Introduction     ")
"Wavefront Propagation",srwUtiShowHelpTopic("Wavefront Propagation     ")
"-"
"FAQ",srwUtiShowHelpTopic("Frequently Asked Questions     ")
"About SRW",srwUtiShowHelpTopic("About SRW     ")
"-"
"Example: Focusing Bending Magnet Radiation",SrwExamImagBm()
"Example: Focusing Undulator Radiation",SrwExamImagUnd()
"Example: Undulator Radiation in 2 Slit Interferometer",SrwExamUndRad2Slits()
"Example: Diffraction of Infra-Red Edge Radiation",SrwExamDiffrER()
"Example: Focusing X-rays by a Refractive Lens",SrwExamXrayLensCirc()
"Example: Gaussian Beam Definition and Propagation",SrwExamGsnBm()
end
SubMenu "Utilities"
"Initialize",SrwInit(1)
"View | Modify SRW Structures...",SrwUtiViewStructDialog()
"Printing to History...",SrwUtiTriggerPrint()
"-"
"Graph: Add Frame and Grid", SrwUtiGraphAddFrameAndGrid()
"Graph: Add Axes Labels...", SrwUtiGraphAddLabels()
"Graph: Resize Window...", SrwUtiGraphWindResize()
"Graph: Save to File...", SrwUtiGraphSaveFile()
"Graphs: Close All",SrwUtiKillAllGraphs()
"-"
"Conversion: Photon Energy...",SrwUtiPhotEnConv()
"Conversion: Spectral Flux...",SrwUtiSpecFluxConv()
"-"
"Optics: Mirror Focal Distance...",SrwUtiOptMirFocLength()
"-"
"Refraction Index...",SrwOptRefrDelta()
end
"-"
SubMenu "Optical Element"
"Aperture: Rectangular...",SrwOptApertRect()
"Aperture: Circular...",SrwOptApertCirc()
"Obstacle: Rectangular...",SrwOptObstRect()
"Obstacle: Circular...",SrwOptObstCirc()
"-"
"Lens: Ideal...",SrwOptThinLens()
"Lens: X-Ray CRL: Parabolic...",SrwOptThinXrayLensParab()
"Lens: X-Ray CRL: Circular...",SrwOptThinXrayLensCirc()
"-"
"Zone Plate...",SrwOptThinZonePlate()
"Zone Plate: Radial Dependence...",SrwOptThinZonePlateRadialMod()
//"Zone Plate...",SrwOptThinGenSetupZonePlateSmpl()
"-"
"Mirror: Toroidal...",SrwOptThinMirTorDialogs()
//"Mirror: Spherical...",SrwOptThinMirSph()
"-"
"Transmission: Template...",SrwOptThinTempl()
"Transmission: Setup...",SrwOptThinSetup()
"Transmission: Display...",SrwOptThinTransmDisplay()
"-"
"Waveguide [test]...",SrwOptWgRect()
end
//SubMenu "Optics Element: Thick"
//"Mirror: Generic...",SrwOptThickMirGenDialogs()
//end

"Drift Space...",SrwOptDrift()
//"Waveguide [test]...",SrwOptWgRect()
"-"
"Container...",SrwOptContFull()
"Add to Container...",SrwOptContAdd()
"-"
"Propagate...",SrwWfrProp()
"Resize...",SrwWfrResize()
"Visualize...",SrwVisualizeDialog()
"Duplicate...",SrwWfrDupl()
"Delete...",SrwWfrDelConfirm()
end // Propagation

//+++++++++++++++++++++++++++++++++++++++
//
// Add-Ons Menu
//
//+++++++++++++++++++++++++++++++++++++++
//Menu "SRW_Add-Ons"
////SubMenu "Report Generation"
//"Add Frame and Grid to Graph", SrwUtiGraphAddFrameAndGrid()
//"Add Axes Labels to Graph...", SrwUtiGraphAddLabels()
//"Resize Graph Window...", SrwUtiGraphWindResize()
//"Save Graph As PNG File...", SrwUtiGraphSave()
//"Close all Open Graphs",SrwUtiKillAllGraphs()
////end
//end

