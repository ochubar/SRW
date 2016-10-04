# -*- coding: utf-8 -*-

"""
The example was created by Timur Shaftan (BNL) for RadTrack project (https://github.com/radiasoft/radtrack).
Adapted by Maksim Rakitin (BNL).
The purpose of the example is to demonstrate good agreement of the SRW simulation
of intensity distribution after diffraction on a circular aperture with an analytical distribution.

The example requires SciPy library to perform comparison.
"""
from __future__ import print_function  # Python 2.7 compatibility
import uti_plot
from srwlib import *

print('SRWLIB Python Example # 16:')
print('Comparison of intensity distribution after diffraction on a circular aperture with an analytical distribution')

#************************************* Create examples directory if it does not exist
example_folder = 'data_example_16'  # example data sub-folder name
if not os.path.isdir(example_folder):
    os.mkdir(example_folder)
strIntOutFileName1 = 'ex16_res_int1.dat' # file name for output SR intensity data
strIntOutFileName2 = 'ex16_res_int2.dat' # file name for output SR intensity data

#************************************* Perform SRW calculations
# Gaussian beam definition:
GsnBm = SRWLGsnBm()
GsnBm.x = 0  # Transverse Coordinates of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0.0  # Longitudinal Coordinate of Waist [m]
GsnBm.xp = 0  # Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = 0.5  # Photon Energy [eV]
GsnBm.pulseEn = 1.0E7  # Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1  # Rep. Rate [Hz] - to be corrected
GsnBm.polar = 6  # 0- Linear Horizontal /  1- Linear Vertical 2- Linear 45 degrees / 3- Linear 135 degrees / 4- Circular Right /  5- Circular /  6- Total
GsnBm.sigX = 2.0E-3  # Horiz. RMS size at Waist [m]
GsnBm.sigY = 2.0E-3  # Vert. RMS size at Waist [m]
GsnBm.sigT = 1E-12  # Pulse duration [fs] (not used?)
GsnBm.mx = 0  # Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

# Wavefront definition:
wfr = SRWLWfr()
NEnergy = 1  # Number of points along Energy
Nx = 300  # Number of points along X
Ny = 300  # Number of points along Y
wfr.allocate(NEnergy, Nx, Ny)  # Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 0.35  # Longitudinal Position [m] at which Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = 0.5  # Initial Photon Energy [eV]
wfr.mesh.eFin = 0.5  # Final Photon Energy [eV]
firstHorAp = 2.0E-3  # First Aperture [m]
firstVertAp = 2.0E-3  # [m]
wfr.mesh.xStart = -0.01  # Initial Horizontal Position [m]
wfr.mesh.xFin = 0.01  # Final Horizontal Position [m]
wfr.mesh.yStart = -0.01  # Initial Vertical Position [m]
wfr.mesh.yFin = 0.01  # Final Vertical Position [m]

# Precision parameters for SR calculation:
meth = 2  # SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
npTraj = 1  # number of points for trajectory calculation (not needed)
relPrec = 0.01  # relative precision
zStartInteg = 0.0  # longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 0.0  # longitudinal position to finish integration (effective if > zStartInteg)
useTermin = 1  # Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
sampFactNxNyForProp = 1  # sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]

# Calculating initial wavefront:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
meshIn = deepcopy(wfr.mesh)
wfrIn = deepcopy(wfr)

arIin = array('f', [0] * wfrIn.mesh.nx * wfrIn.mesh.ny)
srwl.CalcIntFromElecField(arIin, wfrIn, 0, 0, 3, wfr.mesh.eStart, 0, 0)

arIinY = array('f', [0] * wfrIn.mesh.ny)
srwl.CalcIntFromElecField(arIinY, wfrIn, 0, 0, 2, wfrIn.mesh.eStart, 0, 0)  # extracts intensity

# Plotting initial wavefront:
plotMeshInX = [1000 * wfrIn.mesh.xStart, 1000 * wfrIn.mesh.xFin, wfrIn.mesh.nx]
plotMeshInY = [1000 * wfrIn.mesh.yStart, 1000 * wfrIn.mesh.yFin, wfrIn.mesh.ny]
srwl_uti_save_intens_ascii(arIin, wfrIn.mesh,
                           '{}/{}'.format(example_folder, strIntOutFileName1),
                           0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''],
                           _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
uti_plot.uti_plot2d(arIin, plotMeshInX, plotMeshInY,
                    ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation [a.u.]'])
uti_plot.uti_plot1d(arIinY, plotMeshInY, ['Vertical Position [mm]', 'Intensity [a.u.]',
                                          'Intensity Before Propagation\n(cut vs vertical position at x=0)'])

# Element definition:
apertureSize = 0.00075  # Aperture radius, m
driftLength = 1.0  # Drift length, m
OpElement = []
ppOpElement = []

OpElement.append(SRWLOptA('c', 'a', apertureSize, apertureSize))
ppOpElement.append([0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0])
OpElement.append(SRWLOptD(driftLength))  # Drift space
ppOpElement.append([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0])

opBL = SRWLOptC(OpElement, ppOpElement)
srwl.PropagElecField(wfr, opBL)  # Propagate Electric Field

polarization = 6  # 0- Linear Horizontal /  1- Linear Vertical 2- Linear 45 degrees / 3- Linear 135 degrees / 4- Circular Right /  5- Circular /  6- Total
intensity = 0  # 0=Single-e I/1=Multi-e I/2=Single-e F/3=Multi-e F/4=Single-e RadPhase/5=Re single-e Efield/6=Im single-e Efield
dependArg = 3  # 0 - vs e, 1 - vs x, 2 - vs y, 3- vs x&y, 4-vs x&e, 5-vs y&e, 6-vs x&y&e

# Calculating output wavefront:
arIout = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D intensity data
arII = arIout
arIE = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)
srwl.CalcIntFromElecField(arII, wfr, polarization, intensity, dependArg, wfr.mesh.eStart, 0, 0)

arI1y = array('f', [0] * wfr.mesh.ny)
arRe = array('f', [0] * wfr.mesh.ny)
arIm = array('f', [0] * wfr.mesh.ny)
srwl.CalcIntFromElecField(arI1y, wfr, polarization, intensity, 2, wfr.mesh.eStart, 0, 0)  # extracts intensity

# Normalize intensities:
arI1ymax = max(arI1y)
arIinymax = max(arIinY)
for i in range(len(arI1y)):
    arI1y[i] /= arI1ymax
for i in range(len(arIinY)):
    arIinY[i] /= arIinymax

# Plotting output wavefront:
plotNum = 1000
plotMeshx = [plotNum * wfr.mesh.xStart, plotNum * wfr.mesh.xFin, wfr.mesh.nx]
plotMeshy = [plotNum * wfr.mesh.yStart, plotNum * wfr.mesh.yFin, wfr.mesh.ny]
srwl_uti_save_intens_ascii(arII, wfr.mesh,
                           '{}/{}'.format(example_folder, strIntOutFileName2),
                           0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''],
                           _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
uti_plot.uti_plot2d(arII, plotMeshx, plotMeshy,
                    ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intenisty After Propagation [a.u.]'])
uti_plot.uti_plot1d(arI1y, plotMeshy, ['Vertical Position [mm]', 'Intensity [a.u.]',
                                       'Intensity After Propagation\n(cut vs vertical position at x=0)'])

srwl.CalcIntFromElecField(arRe, wfr, polarization, 5, 2, wfr.mesh.eStart, 0, 0)
srwl.CalcIntFromElecField(arIm, wfr, polarization, 6, 2, wfr.mesh.eStart, 0, 0)

parameters = [
    ['Max Inten B and A propagation', arIinymax, arI1ymax],
    ['Num of y mesh pts B and A propagation', len(arIinY), len(arI1y)],
    ['Num of x mesh pts B and A propagation', wfrIn.mesh.nx, wfr.mesh.nx],
    ['Num of y mesh pts B and A propagation', wfrIn.mesh.ny, wfr.mesh.ny],
    ['Wfr xS size [mm] B and A propagation', wfrIn.mesh.xStart, wfr.mesh.xStart],
    ['Wfr xE size [mm] B and A propagation', wfrIn.mesh.xFin, wfr.mesh.xFin],
    ['Wfr yS size [mm] B and A propagation', wfrIn.mesh.yStart, wfr.mesh.yStart],
    ['Wfr yE size [mm] B and A propagation', wfrIn.mesh.yFin, wfr.mesh.yFin],
]
for i in range(len(parameters)):
    print('{}{}: [{}, {}]'.format('    ', parameters[i][0], parameters[i][1], parameters[i][2]))

#************************************* 2. Defining parameters for analytic calculation
lam = 2.4796e-6  # 1.2398/0.5 eV
numPointsIn = len(arIinY)
numPointsOut = len(arI1y)
meshSize = float(wfr.mesh.xFin)

#************************************* 3. Computing intensity distribution as per Born & Wolf, Principles of Optics
th = []
sIn = []
sOut = []
analyticalIntens = []
for i in range(numPointsOut):
    thx = 2.0 * (i - numPointsOut / 2.0 + 0.5) * meshSize / numPointsOut / driftLength
    th.append(thx)
    sOut.append(thx * driftLength * 1000)
try:
    from scipy.special import jv
    for i in range(numPointsOut):
        x = 3.1415 * apertureSize * sin(th[i]) / lam
        analyticalIntens.append((2 * jv(1, x) / x) ** 2)
except:
    pass
for i in range(numPointsIn):
    sIn.append(2000.0 * (i - numPointsIn / 2.0) * float(wfrIn.mesh.xFin) / numPointsIn)

#************************************* 4. Plotting
try:
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sIn, arIinY, '--g.', label='SRW (before aperture)')
    ax.plot(sOut, arI1y, '-r.', label='SRW (after aperture)')
    if analyticalIntens:
        ax.plot(sOut, analyticalIntens, '-b.', label='Analytical estimation')
    ax.legend()
    ax.set_xlabel('Vertical Position [mm]')
    ax.set_ylabel('Normalized Intensity, [a.u.]')
    ax.set_title('Intensity After Propagation\n(cut vs vertical position at x=0)')
    ax.grid()
    plt.savefig('{}/compare.png'.format(example_folder))
    plt.show()
except:
    pass

print('done')
