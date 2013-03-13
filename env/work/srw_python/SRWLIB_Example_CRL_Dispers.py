###########################################################################
# SRWLIB Example#*: Simulating propagation of a Time-Dependent Gaussian X-ray beam through a simple optical scheme containing CRL (with dispersion)
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os
import sys
import math
import random
import copy
import time

print('SRWLIB Python Example:')
print('Simulating propagation of a Time-Dependent Gaussian X-ray beam through a simple optical scheme containing CRL (with dispersion)')
print('Requires ~10 GB of RAM or a bit more !')

#**********************Auxiliary Functions

#Read data comumns from ASCII file:
def AuxReadInDataColumns(filePath, nCol, strSep):
    f = open(filePath, 'r')
    resCols = []
    for iCol in range(nCol):
        resCols.append([])

    curLine = f.readline()
    while len(curLine) > 0:
        curLineParts = curLine.split(strSep)
        for iCol in range(nCol):
            if(iCol < len(curLineParts)):
                resCols[iCol].append(float(curLineParts[iCol]))
        curLine = f.readline()
    f.close()
    return resCols #attn: returns lists, not arrays!

#Write Optical Transmission characteristic data to ASCII file:
def AuxSaveOpTransmData(optTr, t, filePath):
    f = open(filePath, 'w')
    f.write('#C-aligned optical Transmission characteristic (inner loop is vs horizontal position, outer loop vs vertical position)\n')
    f.write('#' + repr(optTr.mesh.eStart) + ' #Initial Photon Energy [eV]\n')
    f.write('#' + repr(optTr.mesh.eFin) + ' #Final Photon Energy [eV]\n')
    f.write('#' + repr(optTr.mesh.ne) + ' #Number of points vs Photon Energy\n')
    f.write('#' + repr(optTr.mesh.xStart) + ' #Initial Horizontal Position [m]\n')
    f.write('#' + repr(optTr.mesh.xFin) + ' #Final Horizontal Position [m]\n')
    f.write('#' + repr(optTr.mesh.nx) + ' #Number of points vs Horizontal Position\n')
    f.write('#' + repr(optTr.mesh.yStart) + ' #Initial Vertical Position [m]\n')
    f.write('#' + repr(optTr.mesh.yFin) + ' #Final Vertical Position [m]\n')
    f.write('#' + repr(optTr.mesh.ny) + ' #Number of points vs Vertical Position\n')
    neLoc = 1
    if(optTr.mesh.ne > 1):
        neLoc = optTr.mesh.ne
    for i in range(neLoc*optTr.mesh.nx*optTr.mesh.ny): #write all data into one column using "C-alignment" as a "flat" 1D array
        tr = 0
        if((t == 1) or (t == 2)): #amplitude or intensity transmission
            tr = optTr.arTr[i*2]
            if(t == 2): #intensity transmission
                tr *= tr
        else: #optical path difference
            tr = optTr.arTr[i*2 + 1]
        f.write(' ' + repr(tr) + '\n')
    f.close()

#Setup Transmission optical element with 1D heght profile data
def AuxTransmAddSurfHeightProfile(optSlopeErr, heightProfData, dim, ang):
    argHeightProfData = heightProfData[0]
    valHeightProfData = heightProfData[1]
    sinAng = math.sin(ang)
    npData = len(heightProfData[0])
    
    #xStep = optSlopeErr.rx/(optSlopeErr.nx - 1)
    #yStep = optSlopeErr.ry/(optSlopeErr.ny - 1)
    #y = optSlopeErr.y - 0.5*optSlopeErr.ry

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    y = auxMesh.yStart
    hApprox = 0
    ipStart = 0
    #for iy in range(optSlopeErr.ny):
    for iy in range(auxMesh.ny):
        if('y' in dim):
            hApprox = 0
            y1 = argHeightProfData[ipStart]*sinAng
            for i in range(ipStart + 1, npData):
                y2 = argHeightProfData[i]*sinAng
                if((y1 <= y) and (y < y2)):
                    hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAng))*(y - y1) + valHeightProfData[i-1]
                    #print(ipStart, i, iy, y1, y, y2, argHeightProfData[i-1], argHeightProfData[i], valHeightProfData[i-1], valHeightProfData[i], hApprox)
                    ipStart = i - 1
                    break
                y1 = y2

        #x = optSlopeErr.x - 0.5*optSlopeErr.rx
        x = auxMesh.xStart
        
        #for ix in range(optSlopeErr.nx):
        for ix in range(auxMesh.nx):
            if('x' in dim):
                if(ix == 0): ipStart = 0
                hApprox = 0
                x1 = argHeightProfData[ipStart]*sinAng
                for i in range(ipStart + 1, npData):
                    x2 = argHeightProfData[i]*sinAng
                    if((x1 <= x) and (x < x2)):
                        hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAng))*(x - x1) + valHeightProfData[i-1]
                        ipStart = i - 1
                        break
                    x1 = x2
            #ofst = 2*ix + (2*optSlopeErr.nx)*iy
            ofst = 2*ix + (2*auxMesh.nx)*iy

            optSlopeErr.arTr[ofst] = 1. #Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0. #Optical Path Difference
            if(hApprox != 0):
                optSlopeErr.arTr[ofst + 1] = -2*sinAng*hApprox #Optical Path Difference (to check sign!)
                #print(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep


#**********************Input Parameters and Structures:
    
strDataFolderName = 'data_example_crl_dispers' #data sub-folder name

GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Coordinates of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Coordinate of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0

#Coherence time (~ Gaussian pulse duration)  0.12 fs @ 15 keV and 0.17 fs @ 7 keV
#Far field angular divergence: 14.1e-6 ./ (ekev) .^0.75

GsnBm.avgPhotEn = 7000. #15000. #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear hoirizontal
#Far field angular divergence: 14.1e-6 ./ (ekev) .^0.75
GsnBm.sigX = 4.30192e-06 #0.17712e-09/(4*Pi)/(14.1e-06/((7)^0.75)) for 7 keV, 3.55561e-06 = 0.0826561e-09/(4*Pi)/(14.1e-06/((15)^0.75)) for 15 keV #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]
#Coherence time (~ Gaussian pulse duration)           0.12 fs @ 15 keV and 0.17 fs @ 7 keV
GsnBm.sigT = 0.12e-15 #0.17e-15 for 15 keV #Pulse duration [s] #To check: Is it 0.12 fs or 12 fs ? 
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

wfr = SRWLWfr() #Initial Electric Field Wavefront

wfr.allocate(200, 200, 200) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.presFT = 1 #Defining Initial Wavefront in Time Domain
#TEST
#wfr.allocate(41, 200, 200) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
#wfr.allocate(1, 200, 200) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
#wfr.presFT = 0 #Defining Initial Wavefront in Time Domain

wfr.avgPhotEn = GsnBm.avgPhotEn

wfr.mesh.eStart = -20*GsnBm.sigT #Initial Time [s]
wfr.mesh.eFin = 20*GsnBm.sigT #Final Time [s]
#TEST
#wfr.mesh.eStart = GsnBm.avgPhotEn - 20.
#wfr.mesh.eFin = GsnBm.avgPhotEn + 20.

wfr.mesh.zStart = 213. #Longitudinal Position [m] at which Electric Field has to be calculated, i.e. the position of the first optical element

firstHorAp = 0.5e-03 #First Aperture [m]
firstVertAp = 0.5e-03 #[m] 
wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

sampFactNxNyForProp = -1 #5 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

#**********************Calculating Initial Wavefront, Changing Domain, Resizing and Extracting Intensity:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)

auxMeshXY = deepcopy(wfr.mesh)
auxMeshXY.ne = 1
arIxy = array('f', [0]*auxMeshXY.nx*auxMeshXY.ny) #"flat" array to take 2D intensity data
#TEST
srwl.CalcIntFromElecField(arIxy, wfr, 6, 0, 3, 0, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIxy, auxMeshXY, os.path.join(os.getcwd(), strDataFolderName, "res_int_t_xy_in.dat"))

#arP = array('d', [0]*auxMeshXY.nx*auxMeshXY.ny) #"flat" array to take 2D phase data (note it should be 'd')
#srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, 0, 0, 0) #extracts radiation phase
#srwl_uti_save_intens_ascii(arP, auxMeshXY, os.path.join(os.getcwd(), strDataFolderName, "res_phase_in.dat"))

#TEST
auxMeshT = deepcopy(wfr.mesh)
auxMeshT.nx = auxMeshT.ny = 1
arIt = array('f', [0]*auxMeshT.ne) #"flat" array to take 1D intensity data (vs Time / Phot. En.)
srwl.CalcIntFromElecField(arIt, wfr, 6, 0, 0, 0, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIt, auxMeshT, os.path.join(os.getcwd(), strDataFolderName, "res_int_t_in.dat"))

#Switching domain from Time to Frequency:
srwl.SetRepresElecField(wfr, 'f');

#Resizing: decreasing Range in Frequency domain to save memory:
srwl.ResizeElecField(wfr, 'f', [0, 0.2, 1]);

srwl.CalcIntFromElecField(arIt, wfr, 6, 0, 0, 0, 0, 0) #extracts intensity
auxMeshF = deepcopy(wfr.mesh)
auxMeshF.nx = auxMeshF.ny = 1
srwl_uti_save_intens_ascii(arIt, auxMeshF, os.path.join(os.getcwd(), strDataFolderName, "res_int_f_in.dat"))

srwl.CalcIntFromElecField(arIxy, wfr, 6, 0, 3, GsnBm.avgPhotEn, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIxy, auxMeshXY, os.path.join(os.getcwd(), strDataFolderName, "res_int_f_xy_in.dat"))

#sys.exit(0)

#**********************Optical Elements

#Defining a perfect 2D parabolic CRL (collimator):
#Be CRL (0.5 mm aperture, collimating)      213  m   # focal distance 213 m @ 7 keV and 593.6 @ 15 keV
print('Setting-up weakly-focusing collimating CRL (without dispersion)...')
delta = 6.9593134e-06 #Be @ 7000 eV #1.51381016e-06 #Be @ 15000 eV
attenLen = 3573.44e-06 #[m] #Be @ 7000 eV #21106.8e-06 #Be @ 15000 eV 
diamCRL = 0.5e-03 #CRL diameter
rMinCRL = 2.9646e-03 #CRL radius at the tip of parabola [m] #1.797e-03 for 15000 eV ?
nCRL = 1 #number of lenses
wallThickCRL = 50e-06 #CRL wall thickness [m]
opColCRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0)
print('done')
#Saving transmission data to file:
AuxSaveOpTransmData(opColCRL, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_col_crl.dat"))

#CRL Aperture (circular)
opApColCRL = SRWLOptA('c', 'a', diamCRL, diamCRL)

#Drift space from collimating CRL to M1
opDrCRL_M1 = SRWLOptD(260 - 213) 

#M1 (mirror1.dat), first Horiz Offset Mirror (HOM) 260 m
#read mirror slope arror data from file and setup the corresponding optical element
print('Defining Trnansmission element (to simulate M1 surface slope error)...', end='')
horApM1 = 1.44e-03 #[m]
vertApM1 = 5e-03 #[m] (dummy?)
angM1 = 1.8e-03 #? [rad]
opTrErM1 = SRWLOptT(1500, 100, horApM1, vertApM1)
heightProfData = AuxReadInDataColumns(os.path.join(os.getcwd(), strDataFolderName, "mirror1.dat"), 2, '\t')
AuxTransmAddSurfHeightProfile(opTrErM1, heightProfData, 'x', angM1)
print('done')
#Saving transmission data to file:
AuxSaveOpTransmData(opTrErM1, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_m1.dat"))

#Drift space from M1 to M2
opDrM1_M2 = SRWLOptD(270 - 260) 

#M2 (mirror2.dat), second HOM 270 m
#read mirror slope arror data from file and setup the corresponding optical element
print('Defining Trnansmission element (to simulate M2 surface slope error)...', end='')
opTrErM2 = SRWLOptT(1500, 100, horApM1, vertApM1)
#can't find good "mirror2.dat", so using "mirror1.dat"
heightProfData = AuxReadInDataColumns(os.path.join(os.getcwd(), strDataFolderName, "mirror1.dat"), 2, '\t')
AuxTransmAddSurfHeightProfile(opTrErM2, heightProfData, 'x', angM1)
print('done')
#Saving transmission data to file:
AuxSaveOpTransmData(opTrErM2, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_m2.dat"))

#Drift space from M2 to Be CRL array
opDrM2_CRL = SRWLOptD(905 - 270) 

#Defining a perfect 2D parabolic CRL (strongly focusing):
#Be CRL array, 5 m focusing distance 905 m
print('Setting-up strongly-focusing CRL (with dispersion)...')

eStartCRL = 6980. # Initial Photon Energy [eV]
eFinCRL = 7020. #Final Photon Energy [eV]
#Refractive Index Decrement tabulated (with equal step) for the photon energy from eStartCRL to eFinCRL
arDelta = [
    6.99931e-06,6.9953e-06,6.9913e-06,6.9873e-06,6.98331e-06,6.97931e-06,6.97532e-06,6.97133e-06,6.96734e-06,6.96335e-06,6.95937e-06,
    6.95539e-06,6.95141e-06,6.94743e-06,6.94346e-06,6.93949e-06,6.93552e-06,6.93155e-06,6.92758e-06,6.92362e-06,6.91966e-06
]

#TEST: Exagerating the dependence of Delta on E:
exFact = 10
yc = arDelta[int(0.5*len(arDelta))]
#print('yc=',yc)
print('Delta data (exaggerated dispersion):')
for ie in range(len(arDelta)):
    arDelta[ie] *= exFact
    arDelta[ie] += (1 - exFact)*yc
    print(arDelta[ie])

#sys.exit(0)

#Attenuation Length tabulated (with equal step) for the photon energy from eStartCRL to eFinCRL
arAttenLen = [
    0.00354301,0.00354604,0.00354908,0.00355212,0.00355516,0.0035582,0.00356125,0.0035643,0.00356734,0.00357039,0.00357344,
    0.00357649,0.00357954,0.00358259,0.00358565,0.00358871,0.00359177,0.00359483,0.00359789,0.00360095,0.00360402
]

diamCRL = 0.5e-03 #CRL diameter
rMinCRL = 0.3479e-03 #CRL radius at the tip of parabola [m] #0.3482e-03 for 15000 eV ?
nCRL = 5 #number of lenses #23 for 15000 eV ?
wallThickCRL = 50e-06 #CRL wall thickness [m]
opFocCRL = srwl_opt_setup_CRL(3, arDelta, arAttenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0, None, eStartCRL, eFinCRL)
#TEST
#opFocCRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0)

print('done')
#Saving transmission data to file (produces very large file):
#AuxSaveOpTransmData(opFocCRL, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_foc_crl.dat"))

#CRL Aperture
opApFocCRL = SRWLOptA('c', 'a', diamCRL, diamCRL)

#Drift space from CRL to Waist
opDrCRL_Waist = SRWLOptD(5.02)  #~5 m for 7 keV

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
#             [0]  [1] [2] [3] [4][5]  [6]  [7]  [8] [9][10][11]
prParRes1 =   [0., 0., 1., 1,  0, 10.0, 1.0, 10.0, 1.0, 0, 0, 0]
prParRes0 =   [0,  0,  1., 1,  0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
prParRes2 =   [0,  0,  1., 1,  0, 0.25, 3.0, 0.25, 3.0, 0, 0, 0]
prParResFin = [0,  0,  1., 1,  0, 0.1, 2.0, 0.1, 2.0, 0, 0, 0]

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
#optBL = SRWLOptC([opApColCRL, opColCRL, opDrCRL_M1, opTrErM1, opDrM1_M2, opTrErM2, opDrM2_CRL, opApFocCRL, opFocCRL, opDrCRL_Waist],
#                     [prParRes1, prParRes0, prParRes0, prParRes0, prParRes0, prParRes0, prParRes0, prParRes0, prParRes0, prParRes0, prParResFin])
optBL = SRWLOptC([opApColCRL, opColCRL, opDrCRL_M1
                  #, opTrErM1
                  , opDrM1_M2
                  #, opTrErM2
                  , opDrM2_CRL
                  , opApFocCRL
                  , opFocCRL
                  , opDrCRL_Waist
                 ],
                 [prParRes1, prParRes0, prParRes0
                  #, prParRes0
                  , prParRes0
                  #, prParRes0
                  , prParRes0
                  , prParRes2
                  , prParRes0
                  , prParRes0
                  , prParResFin
                 ])

#***********Wavefront Propagation

print('Propagating Wavefront...', end='')
startTime = time.time()
srwl.PropagElecField(wfr, optBL)
print('done (lasted', round((time.time() - startTime)/6.)/10., 'min)')

print('Saving resulting data to files...', end='')
auxMeshXY = deepcopy(wfr.mesh)
auxMeshXY.ne = 1
arIxy = array('f', [0]*auxMeshXY.nx*auxMeshXY.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arIxy, wfr, 6, 0, 3, GsnBm.avgPhotEn, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIxy, auxMeshXY, os.path.join(os.getcwd(), strDataFolderName, "res_int_f_xy_out.dat"))

#sys.exit(0)

auxMeshT = deepcopy(wfr.mesh)
auxMeshT.nx = auxMeshT.ny = 1
arIt = array('f', [0]*auxMeshT.ne) #"flat" array to take 1D intensity data (vs Time / Phot. En.)
srwl.CalcIntFromElecField(arIt, wfr, 6, 0, 0, 0, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIt, auxMeshT, os.path.join(os.getcwd(), strDataFolderName, "res_int_f_out.dat"))

#sys.exit(0)

#Resizing: increasing Range in Frequency domain to reduce step size in Time domain:
srwl.ResizeElecField(wfr, 'f', [0, 4., 1]);

#Switching domain from Frequency to Time:
srwl.SetRepresElecField(wfr, 't');

print('Repres. changed to T')

#Resizing: decreasing Range in Time domain to save memory:
#srwl.ResizeElecField(wfr, 'f', [0, 0.5, 1]);

auxMeshXY = deepcopy(wfr.mesh)
auxMeshXY.ne = 1
srwl.CalcIntFromElecField(arIxy, wfr, 6, 0, 3, 0, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIxy, auxMeshXY, os.path.join(os.getcwd(), strDataFolderName, "res_int_t_xy_out.dat"))

auxMeshT = deepcopy(wfr.mesh)
auxMeshT.nx = auxMeshT.ny = 1
arIt = array('f', [0]*auxMeshT.ne) #"flat" array to take 1D intensity data (vs Time / Phot. En.)
srwl.CalcIntFromElecField(arIt, wfr, 6, 0, 0, 0, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arIt, auxMeshT, os.path.join(os.getcwd(), strDataFolderName, "res_int_t_out.dat"))
print('done')
