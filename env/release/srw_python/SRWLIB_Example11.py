# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # 11: Calculating spectral flux of undulator radiation by finite-emittance electron beam
# and performing partially-coherent wavefront propagation through a simple optical system containing dispersive CRL
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os
#import sys
import time

print('SRWLIB Example # 11:')
print('Calculating spectral flux of undulator radiation by finite-emittance electron beam and performing partially-coherent wavefront propagation')
print('through a simple optical system containing dispersive CRL')

#**********************Input Parameters:
strDataFolderName = 'data_example_11' #example data sub-folder name
strOptPathDifCRL_Be = 'ex11_res_opt_path_dif_crl_be.dat' #file name for optical path difference of Be CRL
strOptPathDifCRL_Al = 'ex11_res_opt_path_dif_crl_al.dat' #file name for optical path difference of Al CRL
strFluxOutFileName = 'ex11_res_multi_e_in_spec_flux.dat' #file name for output UR flux data
strTrjOutFileName = 'ex11_res_trj.dat' #file name for output trajectory data
strSingleElSpecOutFileName = 'ex11_res_single_e_spec.dat' #file name for output single-e spectrum data
strSingleElIntOutFileName = 'ex11_res_single_e_int.dat' #file name for output single-e intensity data
strSingleElPropIntOutFileName = 'ex11_res_single_e_prop_int.dat' #file name for output propagated single-e intensity data
strMultiElPropIntOutFileName = 'ex11_res_multi_e_prop_int.dat' #file name for output propagated multi-e intensity data
strSingleGausBeamIntOutFileName = 'ex11_res_single_g_int.dat' #file name for output single-e intensity data
strSingleGausBeamPropIntOutFileName = 'ex11_res_single_g_prop_int.dat' #file name for output propagated single-e intensity data
strMultiGausBeamPropIntOutFileName = 'ex11_res_multi_g_prop_int.dat' #file name for output propagated multi-e intensity data

#***********Undulator
harmB = SRWLMagFldH() #magnetic field harmonic
harmB.n = 1 #harmonic number
harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB.B = 0.391486 #0.920902 #1. #magnetic field amplitude [T]
harmB.s = 0 #assuming anti-symmetric magnetic field vs long. pos.
und = SRWLMagFldU([harmB])
und.per = 0.022 #0.02 #period length [m]
und.nPer = 92 #150 #number of periods (will be rounded to integer)

xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 1.3695 #Longitudinal Coordinate of Undulator Center [m]

magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all magnetic field elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.1 #0.5 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = 6.04/0.5109989e-03 #relative energy
sigEperE = 0.001 #relative RMS energy spread
sigX = 100e-06 #33.33e-06 #horizontal RMS size of e-beam [m]
sigXp = 48.23e-06 #16.5e-06 #horizontal RMS angular divergence [rad]
sigY = 9.525e-06 #2.912e-06 #vertical RMS size of e-beam [m]
sigYp = 3.1623e-06 #2.7472e-06 #vertical RMS angular divergence [rad]
#2nd order stat. moments:
eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2> 
eBeam.arStatMom2[1] = 2.695e-09 #<(x-<x>)(x'-<x'>)> [m]
eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2> 
eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
eBeam.arStatMom2[4] = 0.002695e-09 #<(y-<y>)(y'-<y'>)> [m]
eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

#***********Gaussian Beam (source option)
gBeam = SRWLGsnBm() #Gaussian Beam structure (just parameters)
gBeam.x = 0 #Transverse Coordinates of Gaussian Beam Center at Waist [m]
gBeam.y = 0
gBeam.z = 0 #Longitudinal Coordinate of Waist [m]
gBeam.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
gBeam.yp = 0
gBeam.avgPhotEn = 35.66e+03 #Photon Energy [eV]
gBeam.pulseEn = 0.001 #Energy per Pulse [J] - dummy
gBeam.repRate = 1 #Rep. Rate [Hz] - dummy
gBeam.polar = 1 #1- linear hoirizontal
gBeamWavelength = 1.239842e-06/gBeam.avgPhotEn
gBeamSigDiv = sqrt(gBeamWavelength/(2*und.per*und.nPer)) #Horiz. RMS angular divergence (close to the single-electron UR central cone size value)
gBeamSigSize = gBeamWavelength/(4*3.1416*gBeamSigDiv) #RMS size at Waist [m] (close to the single-electron UR central cone size value)
gBeam.sigX = gBeamSigSize #Horiz. RMS size at Waist [m]
gBeam.sigY = gBeamSigSize #Vert. RMS size at Waist [m]
#print('RMS size of Gaussian Beam at Waist: ', gBeam.sigX)
gBeam.sigT = 1e-27 #Pulse duration [fs] (to have large broadband)
gBeam.mx = 0 #Transverse Gauss-Hermite Mode Orders
gBeam.my = 0

#To edit:
srcToUse = 'e' #'e'- e-beam in undulator; 'g'- Gaussian Beam

#***********Auxiliary Electron Trajectory structure (for test)
partTraj = SRWLPrtTrj() #defining auxiliary trajectory structure
partTraj.partInitCond = eBeam.partStatMom1
partTraj.allocate(20001) 
partTraj.ctStart = 0.2695 #Start "time" for the calculation
partTraj.ctEnd = 2.5

#***********Precision Parameters
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = 9 #final UR harmonic to take into account
arPrecF[2] = 1. #longitudinal integration precision parameter
arPrecF[3] = 1. #azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux within aperture (1) or flux per unit surface (2)

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

arPrecS = [0]*7 #for electric field and single-electron intensity spectrum
arPrecS[0] = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
arPrecS[1] = 0.01 #relative precision
arPrecS[2] = 0 #longitudinal position to start integration (effective if < zEndInteg)
arPrecS[3] = 0 #longitudinal position to finish integration (effective if > zStartInteg)
arPrecS[4] = 20000 #Number of points for intermediate trajectory calculation 
arPrecS[5] = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
arPrecS[6] = -1 #sampling factor for adjusting nx, ny (effective if > 0)

sampFactNxNyForProp = 0.15 #0.1 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecI = deepcopy(arPrecS) #for electric field and single-electron intensity (initial wavefront)
arPrecI[6] = sampFactNxNyForProp

arPrecG = [sampFactNxNyForProp] #for Gaussian Beam
#sys.exit(0)

#***********Defining mesh and allocating memory for UR Stokes parameters and Power Density distributions
stkF = SRWLStokes() #for spectral flux vs photon energy
stkF.allocate(2000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.zStart = 32.87 #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = 10. #initial photon energy [eV]
stkF.mesh.eFin = 40.e+03 #20000. #final photon energy [eV]
stkF.mesh.xStart = -0.0005 #initial horizontal position of the collection aperture [m]
stkF.mesh.xFin = 0.0005 #final horizontal position of the collection aperture [m]
stkF.mesh.yStart = -0.0005 #initial vertical position of the collection aperture [m]
stkF.mesh.yFin = 0.0005 #final vertical position of the collection aperture [m]

wfrS = SRWLWfr() #For on-axis single-electron intensity vs photon energy
wfrS.allocate(10000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfrS.mesh.zStart = stkF.mesh.zStart #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated
wfrS.mesh.eStart = stkF.mesh.eStart #Initial Photon Energy [eV]
wfrS.mesh.eFin = stkF.mesh.eFin #Final Photon Energy [eV]
wfrS.mesh.xStart = 0 #Initial Horizontal Position [m]
wfrS.mesh.xFin = 0 #Final Horizontal Position [m]
wfrS.mesh.yStart = 0 #Initial Vertical Position [m]
wfrS.mesh.yFin = 0 #Final Vertical Position [m]
wfrS.partBeam = eBeam

wfrI = SRWLWfr() #For single-electron intensity at fixed photon energy
wfrI.allocate(1, 201, 201) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfrI.mesh.zStart = stkF.mesh.zStart #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated
wfrI.mesh.eStart = 35.66e+03 #35.3e+03 #Initial Photon Energy [eV]
wfrI.mesh.eFin = wfrI.mesh.eStart #Final Photon Energy [eV]
wfrI.mesh.xStart = stkF.mesh.xStart #Initial Horizontal Position [m]
wfrI.mesh.xFin = stkF.mesh.xFin #Final Horizontal Position [m]
wfrI.mesh.yStart = stkF.mesh.yStart #Initial Vertical Position [m]
wfrI.mesh.yFin = stkF.mesh.yFin #Final Vertical Position [m]
wfrI.partBeam = eBeam
auxMeshP = deepcopy(wfrI.mesh)

#***********Defining optical elements 

#Defining a perfect 2D parabolic CRL (with dispersion):
#Be CRL array, 5 m focusing distance 905 m
print('Setting-up CRL (with dispersion)...')
eStartCRL = 34.01e+03 # Initial Photon Energy [eV]
eFinCRL = 36.11e+03 #Final Photon Energy [eV]
#Refractive Index Decrement tabulated (with constant step) for the photon energy from eStartCRL to eFinCRL
arDeltaBe = [2.93914e-07, 2.92193e-07, 2.90487e-07, 2.88796e-07, 2.8712e-07, 2.85458e-07, 2.8381e-07, 2.82177e-07, 2.80558e-07, 2.78953e-07, 2.77362e-07,
             2.75784e-07, 2.74219e-07, 2.72668e-07, 2.7113e-07, 2.69605e-07, 2.68093e-07, 2.66594e-07, 2.65107e-07, 2.63632e-07, 2.6217e-07, 2.60719e-07]
arDeltaAl = [4.66441e-07, 4.63705e-07, 4.60993e-07, 4.58305e-07, 4.55641e-07, 4.53e-07, 4.50382e-07, 4.47786e-07, 4.45213e-07, 4.42661e-07, 4.40132e-07,
             4.37624e-07, 4.35138e-07, 4.32673e-07, 4.30228e-07, 4.27804e-07, 4.254e-07, 4.23018e-07, 4.20654e-07, 4.18311e-07, 4.15987e-07, 4.13683e-07]
#Attenuatio Length tabulated (with constant step) for the photon energy from eStartCRL to eFinCRL
arAttenLenBe = [31.6869, 31.7151, 31.7433, 31.7714, 31.7994, 31.8271, 31.8537, 31.8802, 31.9066, 31.9329, 31.9592,
                31.9853, 32.0114, 32.0375, 32.0634, 32.0893, 32.1151, 32.1408, 32.1665, 32.1921, 32.2176, 32.2421]
arAttenLenAl = [4.48624, 4.51822, 4.54987, 4.58177, 4.61393, 4.64627, 4.67821, 4.71039, 4.74283, 4.77546, 4.80771,
                4.84021, 4.87296, 4.9059, 4.93842, 4.97119, 5.00425, 5.03741, 5.07024, 5.10331, 5.13662, 5.1696]
diamCRL = 1.e-03 #CRL diameter
rMinCRL = 0.2e-03 #CRL radius at the tip of parabola [m] #0.3482e-03 for 15000 eV ?
wallThickCRL = 50e-06 #CRL wall thickness [m]
nCRL_Be = 16 #number of Be lenses
nCRL_Al = 21 #number of Al lenses

opCRL_Be = srwl_opt_setup_CRL(3, arDeltaBe, arAttenLenBe, 1, diamCRL, diamCRL, rMinCRL, nCRL_Be, wallThickCRL, 0, 0, None, eStartCRL, eFinCRL, 501, 501)
opCRL_Al = srwl_opt_setup_CRL(3, arDeltaAl, arAttenLenAl, 1, diamCRL, diamCRL, rMinCRL, nCRL_Al, wallThickCRL, 0, 0, None, eStartCRL, eFinCRL, 501, 501)

opPathDifCRL_Be = opCRL_Be.get_data(3, 3, _e = wfrI.mesh.eStart)
auxMeshViewCRL_Be = deepcopy(opCRL_Be.mesh)
auxMeshViewCRL_Be.ne = 1
auxMeshViewCRL_Be.eStart = wfrI.mesh.eStart
auxMeshViewCRL_Be.eFin = wfrI.mesh.eStart
srwl_uti_save_intens_ascii(opPathDifCRL_Be, auxMeshViewCRL_Be, os.path.join(os.getcwd(), strDataFolderName, strOptPathDifCRL_Be), 0,
                           ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
opPathDifCRL_Al = opCRL_Al.get_data(3, 3, _e = wfrI.mesh.eStart)
auxMeshViewCRL_Al = deepcopy(opCRL_Al.mesh)
auxMeshViewCRL_Al.ne = 1
auxMeshViewCRL_Al.eStart = wfrI.mesh.eStart
auxMeshViewCRL_Al.eFin = wfrI.mesh.eStart
srwl_uti_save_intens_ascii(opPathDifCRL_Al, auxMeshViewCRL_Al, os.path.join(os.getcwd(), strDataFolderName, strOptPathDifCRL_Al), 0,
                           ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
print('done')

opApCRL = SRWLOptA('c', 'a', diamCRL, diamCRL) #CRL Circular Aperture before CRL

opDriftCRL_Obs = SRWLOptD(9.9519) #Drift space from CRL to ~Waist

#***********Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#                [0] [1] [2] [3] [4] [5] [6]  [7]   [8] [9][10][11]
ppApCRL =        [0., 0., 1., 1,  0, 7.0, 1.0, 3.0, 1.0, 0, 0, 0]
ppCRL_Be =       [0,  0,  1., 1,  0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppCRL_Al =       [0,  0,  1., 1,  0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppDriftCRL_Obs = [0,  0,  1., 2,  0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
ppFinal =        [0,  0,  1., 1,  0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
opBL = SRWLOptC([opApCRL, opCRL_Be, opCRL_Al, opDriftCRL_Obs],
                [ppApCRL, ppCRL_Be, ppCRL_Al, ppDriftCRL_Obs])#, ppFinal])
#sys.exit(0)

#**********************Calculation (SRWLIB function calls)
#print(srwl_uti_proc_is_master())

if(srwl_uti_proc_is_master()):
    print('   Performing Spectral Flux (Stokes parameters) vs Photon Energy calculation ... ', end='')
    srwl.CalcStokesUR(stkF, eBeam, und, arPrecF)
    print('done')

    print('   Performing Electron Trajectory calculation (for test) ... ', end='')
    srwl.CalcPartTraj(partTraj, magFldCnt, [1])
    print('done')

    print('   Performing Single-E Electric Field calculation (vs photon energy) ... ', end='')
    srwl.CalcElecFieldSR(wfrS, 0, magFldCnt, arPrecS)
    print('done')
    print('   Extracting Intensity from the Calculated Electric Field ... ', end='')
    arIS = array('f', [0]*wfrS.mesh.ne) #array to take intensity data
    srwl.CalcIntFromElecField(arIS, wfrS, 6, 0, 0, wfrS.mesh.eStart, 0, 0)
    print('done')

    print('   Performing Single-E Electric Field calculation (vs hor. and vert. positions) ... ', end='')
    if srcToUse == 'e':
        srwl.CalcElecFieldSR(wfrI, 0, magFldCnt, arPrecI)
    elif srcToUse == 'g':
        srwl.CalcElecFieldGaussian(wfrI, gBeam, arPrecG)
    print('done')
    print('   Extracting Intensity from the Calculated Electric Field ... ', end='')
    auxMeshI0 = deepcopy(wfrI.mesh)
    arI = array('f', [0]*auxMeshI0.nx*auxMeshI0.ny) #array to take intensity data
    srwl.CalcIntFromElecField(arI, wfrI, 6, 0, 3, auxMeshI0.eStart, 0, 0)
    print('done')

    print('   Performing Simulation of Single-E Wavefront Propagation ... ', end='')
    startTime = time.time()
    srwl.PropagElecField(wfrI, opBL)
    print('done (lasted', round((time.time() - startTime)/6.)/10., 'min)')

    print('   Extracting Intensity from the Propagated Electric Field ... ', end='')
    auxMeshIP = deepcopy(wfrI.mesh)
    arIP = array('f', [0]*auxMeshIP.nx*auxMeshIP.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIP, wfrI, 6, 0, 3, auxMeshIP.eStart, 0, 0) #extracts intensity
    print('done')

    #***********Saving some intermediate results
    print('   Saving spectral flux [ph/s/.1%bw], flux per unit surface [ph/s/.1%bw/mm^2], and electron trajectory data to files ... ', end='')
    srwl_uti_save_intens_ascii(stkF.arS, stkF.mesh, os.path.join(os.getcwd(), strDataFolderName, strFluxOutFileName))
    partTraj.save_ascii(os.path.join(os.getcwd(), strDataFolderName, strTrjOutFileName))
    srwl_uti_save_intens_ascii(arIS, wfrS.mesh, os.path.join(os.getcwd(), strDataFolderName, strSingleElSpecOutFileName))

    strI0 = strSingleElIntOutFileName
    strIP = strSingleElPropIntOutFileName
    if srcToUse == 'g':
        strI0 = strSingleGausBeamIntOutFileName
        strIP = strSingleGausBeamPropIntOutFileName
    srwl_uti_save_intens_ascii(arI, auxMeshI0, os.path.join(os.getcwd(), strDataFolderName, strI0))
    srwl_uti_save_intens_ascii(arIP, auxMeshIP, os.path.join(os.getcwd(), strDataFolderName, strIP))
    print('done')
#sys.exit(0)
    
print('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... ')
nMacroElec = 50000 #Total number of Macro-Electrons (Wavefronts)
nMacroElecAvgOneProc = 5 #Number of Macro-Electrons (Wavefronts) to average on each node (for MPI calculations)
nMacroElecSavePer = 5 #Saving periodicity (in terms of Macro-Electrons) for the Resulting Intensity
srCalcMeth = 1 #SR calculation method (1- undulator)
srCalcPrec = 0.01 #SR calculation rel. accuracy

if srcToUse == 'e':
    radStokesProp = srwl_wfr_emit_prop_multi_e(eBeam, magFldCnt, auxMeshP, srCalcMeth, srCalcPrec, nMacroElec, nMacroElecAvgOneProc, nMacroElecSavePer, 
                                               os.path.join(os.getcwd(), strDataFolderName, strMultiElPropIntOutFileName), sampFactNxNyForProp, opBL) 
elif srcToUse == 'g':
    radStokesProp = srwl_wfr_emit_prop_multi_e(eBeam, gBeam, auxMeshP, None, None, nMacroElec, nMacroElecAvgOneProc, nMacroElecSavePer, 
                                               os.path.join(os.getcwd(), strDataFolderName, strMultiGausBeamPropIntOutFileName), sampFactNxNyForProp, opBL) 
print('done')
