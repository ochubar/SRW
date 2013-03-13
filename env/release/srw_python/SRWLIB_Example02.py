#############################################################################
# SRWLIB Example#2: Calculating electron trajectory in magnetic field of a segmented planar undulator
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os

print('SRWLIB Python Example # 2:')
print('Calculating electron trajectory in magnetic field of a segmented planar undulator with FODO lattice')

#**********************Input Parameters:
strExDataFolderName = 'data_example_02' #example data sub-folder name
strTrajOutFileName = 'ex02_res_traj.dat' #file name for output trajectory data

numSegm = 5 #Number of ID Segments
numPer = 100 #Number of Periods in one Segment (without counting for terminations)
undPer = 0.02 #Period Length [m]
xcID = 0 #Transverse Coordinates of ID Center [m]
ycID = 0
zcID = 0 #Longitudinal Coordinate of ID Center [m]

part = SRWLParticle()
part.x = 0.0001 #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
part.y = 0.0001
part.xp = 0 #Initial Transverse Velocities
part.yp = 0
part.gamma = 3/0.51099890221e-03 #Relative Energy
part.relE0 = 1 #Electron Rest Mass
part.nq = -1 #Electron Charge

npTraj = 20001 #Number of Points for Trajectory calculation

arPrecPar = [1] #General Precision parameters for Trajectory calculation:
#[0]: integration method No:
    #1- fourth-order Runge-Kutta (precision is driven by number of points)
    #2- fifth-order Runge-Kutta
#[1],[2],[3],[4],[5]: absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] (yet to be tested!!) - to be taken into account only for R-K fifth order or higher
#[6]: tolerance (default = 1) for R-K fifth order or higher
#[7]: max. number of auto-steps for R-K fifth order or higher (default = 5000)

und = SRWLMagFldU([SRWLMagFldH(1, 'v', 1.05, 0, 1)], undPer, numPer) #Undulator Segment
qf = SRWLMagFldM(0.5, 2, 'n', 0.2) #Focusing Quad
qd = SRWLMagFldM(-0.5, 2, 'n', 0.2) #Defocusing Quad

arZero = array('d', [0]*9)
undLen = (numPer + 2)*undPer
distBwSegm = 0.4 #Distance between Undulator Segments
undLenExt = undLen + distBwSegm
arZc = array('d', [-2*undLenExt, -1.5*undLenExt, -undLenExt, -0.5*undLenExt, 0, 0.5*undLenExt, undLenExt, 1.5*undLenExt, 2*undLenExt])
magFldCnt = SRWLMagFldC([und, qf, und, qd, und, qf, und, qd, und], arZero, arZero, arZc) #Container of all Field Elements

part.z = arZc[0] - 0.5*undLenExt #Initial Longitudinal Coordinate (set before the ID)

#**********************Auxiliary function to write tabulated resulting Trajectory data to ASCII file:
def AuxSaveTrajData(traj, filePath):
    f = open(filePath, 'w')
    f.write('#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [rad]\n')
    #print('file opened: ', filePath)
    ctStep = 0
    if traj.np > 0:
        ctStep = (traj.ctEnd - traj.ctStart)/(traj.np - 1)
    ct = traj.ctStart
    for i in range(traj.np):
        f.write(str(ct) + '\t' + repr(traj.arX[i]) + '\t' + repr(traj.arXp[i]) + '\t' + repr(traj.arY[i]) + '\t' + repr(traj.arYp[i]) + '\t' + repr(traj.arZ[i]) + '\t' + repr(traj.arZp[i]) + '\n')        
        ct += ctStep
    f.close()

#**********************Trajectory structure, where the results will be stored
partTraj = SRWLPrtTrj()
partTraj.partInitCond = part
partTraj.allocate(npTraj)
partTraj.ctStart = 0 #"Start Time" (c*t) for the calculation (0 corresponds to the time moment for which the initial conditions are defined)
partTraj.ctEnd = partTraj.ctStart + 5*undLenExt #End Time

#**********************Calculation (SRWLIB function call)
print('   Performing calculation ... ', end='')
partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, arPrecPar)
print('done')

#**********************Saving results
print('   Saving trajectory data to a file ... ', end='')
#raw_input()
AuxSaveTrajData(partTraj, os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName))
#raw_input()
print('done')
