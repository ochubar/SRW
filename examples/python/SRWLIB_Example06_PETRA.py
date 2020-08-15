#############################################################################
# SRWLIB Example#6: Calculating spectral flux of undulator radiation 
# by finite-emittance electron beam collected through a finite aperture
# and power density distribution of this radiation (integrated over all photon energies)
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
#from uti_plot import * #required for plotting
import os
import sys

print('SRWLIB Extended Example # 6:')
print('Calculating spectral flux of undulator radiation by finite-emittance electron beam collected through a finite aperture and power density distribution of this radiation (integrated over all photon energies)')

#**********************Input Parameters:
strExDataFolderName = 'data_example_06' #example data sub-folder name
strFluxOutFileName = 'ex06_res_flux.dat' #file name for output UR flux data
strTrjOutFileName = 'ex06_res_trj.dat' #file name for output trajectory data
strSingleElSpecOutFileName = 'ex06_res_single_e_spec.dat' #file name for output single-e spectrum data
strSingleElIntOutFileName = 'ex06_res_single_e_int.dat' #file name for output single-e intensity data
strPowOutFileName = 'ex06_res_pow.dat' #file name for output power density data
strPowTrjOutFileName = 'ex06_res_pow_trj.dat' #file name for output power density data
strSpecIntOutFileName = 'ex06_res_spec_int.dat' #file name for output UR spectral intensity data

#***********Undulator
harmB = SRWLMagFldH() #magnetic field harmonic
harmB.n = 1 #harmonic number
harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
harmB.B = 0.920902 #1. #magnetic field amplitude [T]
harmB.s = 1 #assuming symmetric magnetic field vs long. pos.
und = SRWLMagFldU([harmB])
und.per = 0.0314 #0.02 #period length [m]
und.nPer = 63 #150 #number of periods (will be rounded to integer)
magFldCnt = SRWLMagFldC([und], array('d', [0]), array('d', [0]), array('d', [0])) #Container of all magnetic field elements

#***********Electron Beam
eBeam = SRWLPartBeam()
eBeam.Iavg = 0.1 #0.5 #average current [A]
eBeam.partStatMom1.x = 0. #initial transverse positions [m]
eBeam.partStatMom1.y = 0.
eBeam.partStatMom1.z = 0. #initial longitudinal positions (set in the middle of undulator)
eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
eBeam.partStatMom1.yp = 0
eBeam.partStatMom1.gamma = 6.08/0.5109989e-03 #3./0.51099890221e-03 #relative energy
sigEperE = 0.0011 #0.00089 #relative RMS energy spread
sigX = 141.4e-06 #33.33e-06 #horizontal RMS size of e-beam [m]
sigXp = 7.071e-06 #16.5e-06 #horizontal RMS angular divergence [rad]
sigY = 4.879e-06 #2.912e-06 #vertical RMS size of e-beam [m]
sigYp = 2.05e-06 #2.7472e-06 #vertical RMS angular divergence [rad]
#2nd order stat. moments:
eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2> 
eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2> 
eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

#***********Auxiliary Electron Trajectory structure (for test)
partTraj = SRWLPrtTrj() #defining auxiliary trajectory structure
partTraj.partInitCond = eBeam.partStatMom1
partTraj.allocate(20001) 
partTraj.ctStart = -1.6 #Start "time" for the calculation
partTraj.ctEnd = 1.6

#***********Precision Parameters
arPrecF = [0]*5 #for spectral flux vs photon energy
arPrecF[0] = 1 #initial UR harmonic to take into account
arPrecF[1] = 21 #final UR harmonic to take into account
arPrecF[2] = 1. #longitudinal integration precision parameter
arPrecF[3] = 1. #azimuthal integration precision parameter
arPrecF[4] = 1 #calculate flux within aperture (1) or flux per unit surface (2)

arPrecP = [0]*5 #for power density
arPrecP[0] = 1.5 #precision factor
arPrecP[1] = 1 #power density computation method (1- "near field", 2- "far field")
arPrecP[2] = 0 #initial longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[3] = 0 #final longitudinal position (effective if arPrecP[2] < arPrecP[3])
arPrecP[4] = 20000 #number of points for (intermediate) trajectory calculation

arPrecS = [0]*7 #for electric field and single-electron intensity
arPrecS[0] = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
arPrecS[1] = 0.01 #relative precision
arPrecS[2] = 0 #longitudinal position to start integration (effective if < zEndInteg)
arPrecS[3] = 0 #longitudinal position to finish integration (effective if > zStartInteg)
arPrecS[4] = 20000 #Number of points for intermediate trajectory calculation 
arPrecS[5] = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
arPrecS[6] = -1 #0.1 #sampling factor for adjusting nx, ny (effective if > 0)

#***********Defining mesh and allocating memory for UR Stokes parameters and Power Density distributions
stkF = SRWLStokes() #for spectral flux vs photon energy
stkF.allocate(10000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
stkF.mesh.zStart = 50. #30. #longitudinal position [m] at which UR has to be calculated
stkF.mesh.eStart = 10. #initial photon energy [eV]
stkF.mesh.eFin = 40.e+03 #20000. #final photon energy [eV]
stkF.mesh.xStart = -0.0015 #initial horizontal position of the collection aperture [m]
stkF.mesh.xFin = 0.0015 #final horizontal position of the collection aperture [m]
stkF.mesh.yStart = -0.001 #initial vertical position of the collection aperture [m]
stkF.mesh.yFin = 0.001 #final vertical position of the collection aperture [m]

wfrS = SRWLWfr() #For on-axis single-electron intensity vs photon energy
wfrS.allocate(20000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
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
wfrI.mesh.eStart = 7.208e+03 #Initial Photon Energy [eV]
wfrI.mesh.eFin = wfrI.mesh.eStart #Final Photon Energy [eV]
wfrI.mesh.xStart = stkF.mesh.xStart #Initial Horizontal Position [m]
wfrI.mesh.xFin = stkF.mesh.xFin #Final Horizontal Position [m]
wfrI.mesh.yStart = stkF.mesh.yStart #Initial Vertical Position [m]
wfrI.mesh.yFin = stkF.mesh.yFin #Final Vertical Position [m]
wfrI.partBeam = eBeam

stkP = SRWLStokes() #for power density
stkP.allocate(1, 101, 101) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
stkP.mesh.zStart = 50. #30. #longitudinal position [m] at which power density has to be calculated
stkP.mesh.xStart = -0.02 #initial horizontal position [m]
stkP.mesh.xFin = 0.02 #final horizontal position [m]
stkP.mesh.yStart = -0.015 #initial vertical position [m]
stkP.mesh.yFin = 0.015 #final vertical position [m]

stkFS = SRWLStokes() #for flux per unit surface vs photon energy, horizontal and vertical position
stkFS.allocate(9, 21, 21) #numbers of points vs photon energy, horizontal and vertical positions
stkFS.mesh.zStart = 50. #longitudinal position [m] at which UR has to be calculated
stkFS.mesh.eStart =  7.14e+03 #10. #initial photon energy [eV]
stkFS.mesh.eFin = 7.26e+03 #final photon energy [eV]
stkFS.mesh.xStart = -0.002 #initial horizontal position of the collection aperture [m]
stkFS.mesh.xFin = 0.002 #final horizontal position of the collection aperture [m]
stkFS.mesh.yStart = -0.0015 #initial vertical position of the collection aperture [m]
stkFS.mesh.yFin = 0.0015 #final vertical position of the collection aperture [m]

#sys.exit(0)

#**********************Calculation (SRWLIB function calls)

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
srwl.CalcElecFieldSR(wfrI, 0, magFldCnt, arPrecS)
print('done')
print('   Extracting Intensity from the Calculated Electric Field ... ', end='')
arI = array('f', [0]*wfrI.mesh.nx*wfrI.mesh.ny) #array to take intensity data
srwl.CalcIntFromElecField(arI, wfrI, 6, 0, 3, wfrI.mesh.eStart, 0, 0)
print('done')

print('   Performing Power Density calculation (from magnetic field) ... ', end='')
srwl.CalcPowDenSR(stkP, eBeam, 0, magFldCnt, arPrecP)
print('done')

print('   Performing Power Density calculation (from trajectory) ... ', end='')
stkP1 = deepcopy(stkP)
srwl.CalcPowDenSR(stkP1, eBeam, partTraj, 0, arPrecP)
print('done')

print('   Performing calculation of Spectral Intensity (Stokes parameters) vs photon energy, hor. and vert. positions (takes a bit of time) ... ', end='')
arPrecF[0] = 3 #initial UR harmonic to take into account
arPrecF[1] = 4 #final UR harmonic to take into account
arPrecF[2] = 1. #longitudinal integration precision parameter
arPrecF[3] = 1. #azimuthal integration precision parameter
arPrecF[4] = 2 #calculate flux within aperture (1) or flux per unit surface (2)
#srwl.CalcStokesUR(stkFS, eBeam, und, arPrecF)
print('done')

#**********************Saving results
#Auxiliary function to write tabulated Trajectory data to ASCII file:
def AuxSaveTrajData(traj, filePath):
    f = open(filePath, 'w')
    resStr = '#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [rad]'
    if(hasattr(traj, 'arBx')):
        resStr += ', Bx [T]'
    if(hasattr(traj, 'arBy')):
        resStr += ', By [T]'
    if(hasattr(traj, 'arBz')):
        resStr += ', Bz [T]'
    f.write(resStr + '\n')
    ctStep = 0
    if traj.np > 0:
        ctStep = (traj.ctEnd - traj.ctStart)/(traj.np - 1)
    ct = traj.ctStart
    for i in range(traj.np):
        resStr = str(ct) + '\t' + repr(traj.arX[i]) + '\t' + repr(traj.arXp[i]) + '\t' + repr(traj.arY[i]) + '\t' + repr(traj.arYp[i]) + '\t' + repr(traj.arZ[i]) + '\t' + repr(traj.arZp[i])
        if(hasattr(traj, 'arBx')):
            resStr += '\t' + repr(traj.arBx[i])
        if(hasattr(traj, 'arBy')):
            resStr += '\t' + repr(traj.arBy[i])
        if(hasattr(traj, 'arBz')):
            resStr += '\t' + repr(traj.arBz[i])
        f.write(resStr + '\n')        
        ct += ctStep
    f.close()

print('   Saving spectral flux [ph/s/.1%bw], flux per unit surface [ph/s/.1%bw/mm^2], poower density [W/mm^2], and electron trajectory data to files ... ', end='')

srwl_uti_save_intens_ascii(stkF.arS, stkF.mesh, os.path.join(os.getcwd(), strExDataFolderName, strFluxOutFileName), 0, ['Photon Energy', '', '', 'Flux'], _arUnits=['eV', '', '', 'ph/s/.1%bw'])

AuxSaveTrajData(partTraj, os.path.join(os.getcwd(), strExDataFolderName, strTrjOutFileName))
srwl_uti_save_intens_ascii(arIS, wfrS.mesh, os.path.join(os.getcwd(), strExDataFolderName, strSingleElSpecOutFileName))

srwl_uti_save_intens_ascii(arI, wfrI.mesh, os.path.join(os.getcwd(), strExDataFolderName, strSingleElIntOutFileName))

srwl_uti_save_intens_ascii(stkP.arS, stkP.mesh, os.path.join(os.getcwd(), strExDataFolderName, strPowOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Power Density'], _arUnits=['', 'm', 'm', 'W/mm^2'])

srwl_uti_save_intens_ascii(stkP1.arS, stkP1.mesh, os.path.join(os.getcwd(), strExDataFolderName, strPowTrjOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Power Density'], _arUnits=['', 'm', 'm', 'W/mm^2'])

#AuxSaveS0Data(stkFS, os.path.join(os.getcwd(), strExDataFolderName, strSpecIntOutFileName))

print('done')
