#############################################################################
# SRWLIB Example#3: Calculating synchrotron (undulator) radiation emitted by an electron travelling in ellipsoidal undulator
# v 0.04
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os

print('SRWLIB Python Example # 3:')
print('Calculating synchrotron (undulator) radiation emitted by an electron travelling in a helical undulator')

#**********************Input Parameters:
strExDataFolderName = 'data_example_03' #example data sub-folder name
strTrajOutFileName = 'ex03_res_traj.dat' #file name for output trajectory data
strIntOutFileName1 = 'ex03_res_int1.dat' #file name for output SR intensity data
strIntOutFileName2 = 'ex03_res_int2.dat' #file name for output SR intensity data

#***********Undulator
numPer = 40.5 #Number of ID Periods (without counting for terminations
undPer = 0.049 #Period Length [m]
Bx = 0.57/3. #Peak Horizontal field [T]
By = 0.57 #Peak Vertical field [T]
phBx = 0 #Initial Phase of the Horizontal field component
phBy = 0 #Initial Phase of the Vertical field component
sBx = -1 #Symmetry of the Horizontal field component vs Longitudinal position
sBy = 1 #Symmetry of the Vertical field component vs Longitudinal position
xcID = 0 #Transverse Coordinates of Undulator Center [m]
ycID = 0
zcID = 0 #Longitudinal Coordinate of Undulator Center [m]

und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1), SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer, numPer) #Ellipsoidal Undulator
magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field Elements

#***********Electron Beam
elecBeam = SRWLPartBeam()
elecBeam.Iavg = 0.5 #Average Current [A]
elecBeam.partStatMom1.x = 0. #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
elecBeam.partStatMom1.y = 0.
elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
elecBeam.partStatMom1.xp = 0 #Initial Relative Transverse Velocities
elecBeam.partStatMom1.yp = 0
elecBeam.partStatMom1.gamma = 3./0.51099890221e-03 #Relative Energy

#***********Precision
meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.01 #relative precision
zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000
sampFactNxNyForProp = 0 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp]

#***********Wavefront
wfr1 = SRWLWfr() #For spectrum vs photon energy

wfr1.allocate(10000, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr1.mesh.zStart = 20. #Longitudinal Position [m] at which SR has to be calculated
wfr1.mesh.eStart = 10. #Initial Photon Energy [eV]
wfr1.mesh.eFin = 3000. #Final Photon Energy [eV]
wfr1.mesh.xStart = 0. #Initial Horizontal Position [m]
wfr1.mesh.xFin = 0 #Final Horizontal Position [m]
wfr1.mesh.yStart = 0 #Initial Vertical Position [m]
wfr1.mesh.yFin = 0 #Final Vertical Position [m]
wfr1.partBeam = elecBeam

wfr2 = SRWLWfr() #For intensity distribution at fixed photon energy
wfr2.allocate(1, 101, 101) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
wfr2.mesh.zStart = 20. #Longitudinal Position [m] at which SR has to be calculated
wfr2.mesh.eStart = 1090. #Initial Photon Energy [eV]
wfr2.mesh.eFin = 1090. #Final Photon Energy [eV]
wfr2.mesh.xStart = -0.001 #Initial Horizontal Position [m]
wfr2.mesh.xFin = 0.001 #Final Horizontal Position [m]
wfr2.mesh.yStart = -0.001 #Initial Vertical Position [m]
wfr2.mesh.yFin = 0.001 #Final Vertical Position [m]
wfr2.partBeam = elecBeam

#**********************Auxiliary function to write tabulated resulting Intensity data to ASCII file:
def AuxSaveIntData(arI, mesh, filePath):
    f = open(filePath, 'w')
    f.write('#C-aligned Intensity (inner loop is vs photon energy, outer loop vs vertical position)\n')
    f.write('#' + repr(mesh.eStart) + ' #Initial Photon Energy [eV]\n')
    f.write('#' + repr(mesh.eFin) + ' #Final Photon Energy [eV]\n')
    f.write('#' + repr(mesh.ne) + ' #Number of points vs Photon Energy\n')
    f.write('#' + repr(mesh.xStart) + ' #Initial Horizontal Position [m]\n')
    f.write('#' + repr(mesh.xFin) + ' #Final Horizontal Position [m]\n')
    f.write('#' + repr(mesh.nx) + ' #Number of points vs Horizontal Position\n')
    f.write('#' + repr(mesh.yStart) + ' #Initial Vertical Position [m]\n')
    f.write('#' + repr(mesh.yFin) + ' #Final Vertical Position [m]\n')
    f.write('#' + repr(mesh.ny) + ' #Number of points vs Vertical Position\n')
    for i in range(mesh.ne*mesh.nx*mesh.ny): #write all data into one column using "C-alignment" as a "flat" 1D array
        f.write(' ' + repr(arI[i]) + '\n')
    f.close()

#**********************Calculation (SRWLIB function calls)
print('   Performing Electric Field (spectrum vs photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
print('done')
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI1 = array('f', [0]*wfr1.mesh.ne)
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 0, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
print('done')
print('   Performing Electric Field (wavefront at fixed photon energy) calculation ... ', end='')
srwl.CalcElecFieldSR(wfr2, 0, magFldCnt, arPrecPar)
print('done')
print('   Extracting Intensity from calculated Electric Field ... ', end='')
arI2 = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0)
print('done')

#**********************Saving results
print('   Saving intensity data to files ... ', end='')
AuxSaveIntData(arI1, wfr1.mesh, os.path.join(os.getcwd(), strExDataFolderName, strIntOutFileName1))
AuxSaveIntData(arI2, wfr2.mesh, os.path.join(os.getcwd(), strExDataFolderName, strIntOutFileName2))
print('done')
