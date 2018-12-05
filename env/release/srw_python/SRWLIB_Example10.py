# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # 10: Simulating emission and propagation of partially-coherent undulator radiation 
# through a microscopy beamline with a secondary source aperture and ellopsoidal K-B mirrors used for final focusing
# v 0.07
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os
import time

print('SRWLIB Python Example # 10:')
print('!!!!!Under testing!!!!!')
print('Simulating emission and propagation of undulator radiation (UR) wavefront through a microscopy beamline with a secondary source aperture and ellopsoidal K-B mirrors used for final focusing')
print('')
print('First, single-electron UR wavefront at a fixed photon energy is calculated and propagated through the optical scheme. ', end='')
print('After this, calculation of partially-coherent UR from entire electron beam is started as a loop over "macro-electrons", using "srwl_wfr_emit_prop_multi_e" function. ', end='')
print('This function can run either in "normal" sequential mode, or in parallel mode under "mpi4py".', end='')
print('For this, an MPI2 package and the "mpi4py" Python package have to be installed and configured, and this example has to be started e.g. as:')
print('    mpiexec -n 5 python SRWLIB_Example10.py')
print('For more information on parallel calculations under "mpi4py" please see documentation to the "mpi4py" and MPI.')
print('Note that the long-lasting partially-coherent UR calculation saves from time to time instant average intensity to an ASCII file, ', end='')
print('so the execution of the long loop over "macro-electrons" can be aborted after some time without the danger that all results will be lost.')
print('')

#**********************Input Parameters:
strExDataFolderName = 'data_example_10' #example data sub-folder name
strIntSE_OutFileName = 'ex10_res_int_se.dat' #file name for output initial single-electron SR intensity data
strIntPropSE_OutFileName = 'ex10_res_int_prop_se.dat' #file name for output propagated single-electron SR intensity data
strIntPropME_OutFileName = 'ex10_res_int_prop_me.dat' #file name for output propagated multi-electron SR intensity data

#***********Undulator
numPer = 71.5 #Number of ID Periods (without counting for terminations
undPer = 0.021 #Period Length [m]
By = 0.80371 #Peak Vertical field [T]
phBy = 0 #Initial Phase of the Vertical field component
sBy = -1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 1.25 #0 #Longitudinal Coordinate of Undulator Center [m]

und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1)], undPer, numPer) #Planar Undulator
magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements

#***********Electron Beam
elecBeam = SRWLPartBeam()
elecBeam.Iavg = 0.5 #Average Current [A]
elecBeam.partStatMom1.x = 0. #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
elecBeam.partStatMom1.y = 0. #-0.00025
elecBeam.partStatMom1.z = 0. #-0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
elecBeam.partStatMom1.xp = 0. #Initial Relative Transverse Velocities
elecBeam.partStatMom1.yp = 0.
elecBeam.partStatMom1.gamma = 3./0.51099890221e-03 #Relative Energy
#2nd order statistical moments
elecBeam.arStatMom2[0] = (33.3317e-06)**2 #<(x-x0)^2>
elecBeam.arStatMom2[1] = 0
elecBeam.arStatMom2[2] = (16.5008e-06)**2 #<(x'-x'0)^2>
elecBeam.arStatMom2[3] = (2.91204e-06)**2 #<(y-y0)^2>
elecBeam.arStatMom2[4] = 0
elecBeam.arStatMom2[5] = (2.74721e-06)**2 #<(y'-y'0)^2>
elecBeam.arStatMom2[10] = (0.89e-03)**2 #<(E-E0)^2>/E0^2

#***********Precision Parameters for SR calculation
meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.01 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 0.2 #0.6 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecParSpec = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, 0]

#****************** Initial Wavefront
wfr = SRWLWfr() #For intensity distribution at fixed photon energy
wfr.allocate(1, 101, 101) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr.mesh.zStart = 20 #Longitudinal Position [m] from Center of Straight Section at which SR has to be calculated
wfr.mesh.eStart = 12700 #Initial Photon Energy [eV]
wfr.mesh.eFin = wfr.mesh.eStart #Final Photon Energy [eV]
wfr.mesh.xStart = -0.0004 #-0.00015 #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.0004 #0.00015 #Final Horizontal Position [m]
wfr.mesh.yStart = -0.00025 #-0.00015 #Initial Vertical Position [m]
wfr.mesh.yFin = 0.00025 #0.00015 #Final Vertical Position [m]
meshInitPartCoh = deepcopy(wfr.mesh)
wfr.partBeam = elecBeam

#***************** Optical Elements and Propagation Parameters

Drift_Slits_HFM = SRWLOptD(22.) #Drift from first Slits to Horizontally-Focusing Mirror (HFM)

HFM = SRWLOptL(9.92727) #HFM as Thin Lens

Drift_HFM_SSA = SRWLOptD(13.) #Drift from HFM to Secondary Source Aperture (SSA)

SSA = SRWLOptA('r', 'a', 30e-06, 10e-03)

Drift_SSA_VKB = SRWLOptD(11.) #Drift from SSA to Center of Vertically Focusing K-B Mirror (VKB)

twoSlitPlane = None #'h' #'v'  #plane for the two slit interference scheme (should be 'h' or 'v' for calculating interference)
#twoSlitInterfAngRepres = False #calculate interference in angular representation
twoSlitSep = 0.2e-03 #separation distance between two slits
twoSlitSize = 0.002e-03 #size of each slit
twoSlitObstaclePos = 0 #0.0021e-03 #center position of obstacle (e.g. for calculating diffraction on 1 slit)

TwoSlit_Drift = SRWLOptD(20.) #Drift space after 2 slits

twoSlitSizeAp = twoSlitSep + twoSlitSize
twoSlitSizeOb = twoSlitSep - twoSlitSize
TwoSlit_Ap = None; TwoSlit_Ob = None;
if twoSlitPlane == 'h':
    TwoSlit_Ap = SRWLOptA('r', 'a', twoSlitSizeAp, 10e-03)
    TwoSlit_Ob = SRWLOptA('r', 'o', twoSlitSizeOb, 10e-03, twoSlitObstaclePos, 0)
elif twoSlitPlane == 'v':
    TwoSlit_Ap = SRWLOptA('r', 'a', 10e-03, twoSlitSizeAp)
    TwoSlit_Ob = SRWLOptA('r', 'o', 10e-03, twoSlitSizeOb, 0, twoSlitObstaclePos)

ApKB = SRWLOptA('r', 'a', 1.25e-03, 1.25e-03) #Aperture Before K-B

angVKB = 2.5e-03 #grazing angle at VKB center [rad]
VKB = SRWLOptMirEl(_p=64.75, _q=1., _ang_graz=angVKB, _r_sag=1.e+40, _size_tang=0.5, _nvx=0, _nvy=cos(angVKB), _nvz=-sin(angVKB), _tvx=0, _tvy=-sin(angVKB), _x=0, _y=0, _treat_in_out=1) #VKB Ellipsoidal Mirror

Drift_VKB_HKB = SRWLOptD(0.5) #Distance between centers of Vertically and Horizontally focusing K-B mirrors 

angHKB = 2.5e-03 #[rad]
HKB = SRWLOptMirEl(_p=11.5, _q=0.5, _ang_graz=angHKB, _r_sag=1.e+40, _size_tang=0.5, _nvx=cos(angHKB), _nvy=0, _nvz=-sin(angHKB), _tvx=-sin(angHKB), _tvy=0, _x=0, _y=0, _treat_in_out=1) #HKB Ellipsoidal Mirror
#HKB = SRWLOptL(_Fx=11.5*0.5/(11.5 + 0.5))

Drift_HKB_Sample = SRWLOptD(0.5) #Drift from HKB Center to Sample

#Wavefront Propagation Parameters:
#                    [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
ppDrift_Slits_HFM =  [ 0,  0, 1.0,  1,  0, 4.0, 4.0, 3.0, 3.0,  0,  0,   0]
ppHFM =              [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_HFM_SSA =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
#ppSSA =              [ 0,  0, 1.0,  0,  0, 1.0, 3.0, 1.0, 1.0,  0,  0,   0]
ppSSA =              [ 0,  0, 1.0,  0,  0, 1.0, 6.0, 1.0, 1.0,  0,  0,   0]
ppDrift_SSA_VKB =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppTwoSlit_Ap = None
if twoSlitPlane == 'h':
    ppTwoSlit_Ap =   [ 0,  0, 1.0,  1,  0, 0.5, 7.0, 0.5, 0.5,  0,  0,   0]
    #ppTwoSlit_Ap =   [ 0,  0, 1.0,  1,  0, 1.0, 30.0, 1.0, 1.0,  0,  0,   0]
elif twoSlitPlane == 'v':
    ppTwoSlit_Ap =   [ 0,  0, 1.0,  1,  0, 0.5, 0.04, 0.8, 40.0,  0,  0,   0]
    #ppTwoSlit_Ap =   [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 40.0,  0,  0,   0]
ppTwoSlit_Ob =       [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppTwoSlit_Drift =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

ppApKB =             [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppVKB =              [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_VKB_HKB =    [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppHKB =              [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppDrift_HKB_Sample = [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]
ppFinal =            [ 0,  0, 1.0,  0,  1, 0.2, 2.0, 0.2, 4.0,  0,  0,   0]

#[ 0]: Auto-Resize (1) or not (0) Before propagation
#[ 1]: Auto-Resize (1) or not (0) After propagation
#[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[ 6]: Horizontal Resolution modification factor at Resizing
#[ 7]: Vertical Range modification factor at Resizing
#[ 8]: Vertical Resolution modification factor at Resizing
#[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
if (twoSlitPlane == 'h') or (twoSlitPlane == 'v'):
    optBL = SRWLOptC([Drift_Slits_HFM,   HFM,   Drift_HFM_SSA,   SSA,   Drift_SSA_VKB,   TwoSlit_Ap,   TwoSlit_Ob,   TwoSlit_Drift],#   ApKB,   VKB,   Drift_VKB_HKB,   HKB,   Drift_HKB_Sample], 
                     [ppDrift_Slits_HFM, ppHFM, ppDrift_HFM_SSA, ppSSA, ppDrift_SSA_VKB, ppTwoSlit_Ap, ppTwoSlit_Ob, ppTwoSlit_Drift])#, ppApKB, ppVKB, ppDrift_VKB_HKB, ppHKB, ppDrift_HKB_Sample, ppFinal]) 
else:
    optBL = SRWLOptC([Drift_Slits_HFM,   HFM,   Drift_HFM_SSA,   SSA,   Drift_SSA_VKB,   ApKB,   VKB,   Drift_VKB_HKB,   HKB,   Drift_HKB_Sample], 
                     [ppDrift_Slits_HFM, ppHFM, ppDrift_HFM_SSA, ppSSA, ppDrift_SSA_VKB, ppApKB, ppVKB, ppDrift_VKB_HKB, ppHKB, ppDrift_HKB_Sample, ppFinal])

#**********************Calculation (SRWLIB function calls)
if(srwl_uti_proc_is_master()):
#if(False):
    print('   Performing Initial Single-E Electric Field calculation ... ', end='')
    arPrecParSpec[6] = sampFactNxNyForProp #sampling factor for adjusting nx, ny (effective if > 0)
    srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecParSpec)
    print('done')
    print('   Extracting Intensity from the Calculated Initial Electric Field ... ', end='')
    arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    print('done')
    print('   Saving the Initial Wavefront Intensity into a file ... ', end='')
    srwl_uti_save_intens_ascii(arI, wfr.mesh, os.path.join(os.getcwd(), strExDataFolderName, strIntSE_OutFileName))
    print('done')

    print('   Simulating Electric Field Wavefront Propagation ... ', end='')
    t0 = time.time();
    srwl.PropagElecField(wfr, optBL)
    print('done; lasted', round(time.time() - t0), 's')

    #if twoSlitInterfAngRepres: srwl.SetRepresElecField(wfr, 'a')
    
    print('   Extracting Intensity from the Propagated Electric Field  ... ', end='')
    arI1 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
    srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    #srwl.CalcIntFromElecField(arI1, wfr, 0, 6, 3, wfr.mesh.eStart, 0, 0) #test (E-field extraction)
    print('done')
    print('   Saving the Propagated Wavefront Intensity data to a file ... ', end='')
    srwl_uti_save_intens_ascii(arI1, wfr.mesh, os.path.join(os.getcwd(), strExDataFolderName, strIntPropSE_OutFileName))
    print('done')

#sys.exit()

print('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... ')
nMacroElec = 50000 #Total number of Macro-Electrons (Wavefronts)
nMacroElecAvgOneProc = 5 #Number of Macro-Electrons (Wavefronts) to average on each node (for MPI calculations)
nMacroElecSavePer = 5 #Saving periodicity (in terms of Macro-Electrons) for the Resulting Intensity
srCalcMeth = 1 #SR calculation method (1- undulator)
srCalcPrec = 0.01 #SR calculation rel. accuracy
radStokesProp = srwl_wfr_emit_prop_multi_e(elecBeam, magFldCnt, meshInitPartCoh, srCalcMeth, srCalcPrec, nMacroElec, nMacroElecAvgOneProc, nMacroElecSavePer, 
                                           os.path.join(os.getcwd(), strExDataFolderName, strIntPropME_OutFileName), sampFactNxNyForProp, optBL)
print('done')

