###########################################################################
# SRWLIB Example#7: Simulating propagation of a Gaussian X-ray beam through a simple optical scheme containing CRL
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
import os
#import sys
#import math
import random
import copy

print('SRWLIB Python Example # 7:')
print('Simulating propagation of a Gaussian X-ray beam through a simple optical scheme containing CRL')

#**********************Auxiliary Functions

#Write tabulated resulting Intensity data to ASCII file:
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

#def AuxSaveOpTransmData(optTr, t, filePath):
#    f = open(filePath, 'w')
#    f.write('#C-aligned optical Transmission characteristic (inner loop is vs horizontal position, outer loop vs vertical position)\n')
#    f.write('#' + repr(1) + ' #Reserved for Initial Photon Energy [eV]\n')
#    f.write('#' + repr(1) + ' #Reserved for Final Photon Energy [eV]\n')
#    f.write('#' + repr(1) + ' #Reserved for Number of points vs Photon Energy\n')
#    f.write('#' + repr(optTr.x - 0.5*optTr.rx) + ' #Initial Horizontal Position [m]\n')
#    f.write('#' + repr(optTr.x + 0.5*optTr.rx) + ' #Final Horizontal Position [m]\n')
#    f.write('#' + repr(optTr.nx) + ' #Number of points vs Horizontal Position\n')
#    f.write('#' + repr(optTr.y - 0.5*optTr.ry) + ' #Initial Vertical Position [m]\n')
#    f.write('#' + repr(optTr.y + 0.5*optTr.ry) + ' #Final Vertical Position [m]\n')
#    f.write('#' + repr(optTr.ny) + ' #Number of points vs Vertical Position\n')
#    for i in range(optTr.nx*optTr.ny): #write all data into one column using "C-alignment" as a "flat" 1D array
#        tr = 0
#        if((t == 1) or (t == 2)): #amplitude or intensity transmission
#            tr = optTr.arTr[i*2]
#            if(t == 2): #intensity transmission
#                tr *= tr
#        else: #optical path difference
#            tr = optTr.arTr[i*2 + 1]
#        f.write(' ' + repr(tr) + '\n')
#    f.close()

#**********************Input Parameters and structures:
    
strDataFolderName = 'data_example_07' #data sub-folder name

GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Coordinates of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Coordinate of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = 12400 #5000 #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear hoirizontal
GsnBm.sigX = 23e-06/2.35 #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]
GsnBm.sigT = 10e-15 #Pulse duration [fs] (not used?)
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, 100, 100) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 300 #Longitudinal Position [m] at which Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]
firstHorAp = 1.e-03 #First Aperture [m]
firstVertAp = 1.e-03 #[m] 
wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

sampFactNxNyForProp = 5 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

#**********************Calculating Initial Wavefront and extracting Intensity:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
AuxSaveIntData(arI, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_int_in.dat"))

arP = array('d', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0) #extracts radiation phase
AuxSaveIntData(arP, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_phase_in.dat"))

#***********Optical Elements

delta = 2.21555115e-06 #Be @ 12400 eV #1.3657359E-05 #Be @ 5000 eV
attenLen = 15453.7e-06 #[m] #1268.65e-06
diamCRL = 1.e-03 #CRL diameter
rMinCRL = 250e-06 #CRL radius at the tip of parabola [m]
nCRL = 6 #number of lenses
wallThickCRL = 50e-06 #CRL wall thickness [m]

opCRLperf = 1 #set to 1 if simulation of perfect CRL is required, otherwise set to None
opCRLdist = None #set to 1 if simulation of distorted CRL is required, otherwise set to None

#Generating a perfect 2D parabolic CRL:
if(opCRLperf != None):
    print('Setting-up Perfect CRL...')
    opCRLperf = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0)
    print('done')
    #Saving transmission data to file
    AuxSaveOpTransmData(opCRLperf, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_perf_crl.dat"))

if(opCRLdist != None):
    print('Setting-up CRL Distorted by \"voids\" in volume (takes time)...')
    #Generating array of Void Centers and Radii, and a CRL with these (spherical) Voids:
    nVoidInRectPar = 1000 #parameter controlling density of voids
    baseNx = 201 #(auxiliary) numbers of points in tabulated functions to be used (in rejectrion method)
    baseNy = 201
    rMinVoid = 2e-06 #min. void radius
    rMaxVoid = 25e-06 #max. void radius
    baseYvsXmin = array('d', [0]*baseNx)
    baseYvsXmax = array('d', [0]*baseNx)
    rCRL = 0.5*diamCRL
    rCRLe2 = rCRL*rCRL
    dx = diamCRL/(baseNx - 1)
    x = -rCRL
    for ix in range(baseNx):
        baseYvsXmin[ix] = -sqrt(abs(rCRLe2 - x*x))
        baseYvsXmax[ix] = -baseYvsXmin[ix]
        x += dx
    baseNpSurf = baseNx*baseNy
    baseZvsXYmin = array('d', [0]*baseNpSurf)
    baseZvsXYmax = array('d', [0]*baseNpSurf)
    a = 0.5/rMinCRL
    halfWallThickCRL = 0.5*wallThickCRL
    y = -rCRL
    it = 0
    for iy in range(baseNy):
        ye2 = y*y
        x = -rCRL
        for ix in range(baseNx):
            baseZvsXYmin[it] = -halfWallThickCRL - a*(x*x + ye2)
            baseZvsXYmax[it] = -baseZvsXYmin[it]
            it += 1
            x += dx
        y += dx
    arVoidCenCoordInCRL = srwl_uti_rand_fill_vol(nVoidInRectPar, -rCRL, rCRL, baseNx, baseYvsXmin, baseYvsXmax, -rCRL, rCRL, baseNy, baseZvsXYmin, baseZvsXYmax)
    #Replacing longitudinal coordinate by random Radii:
    random.seed()
    for i in range(int(len(arVoidCenCoordInCRL)/3)):
        arVoidCenCoordInCRL[3*i + 2] = rMinVoid + (rMaxVoid - rMinVoid)*random.random()
        #print('Void Coord. and Rad. (x, y, r):', arVoidCenCoordInCRL[3*i], arVoidCenCoordInCRL[3*i + 1], arVoidCenCoordInCRL[3*i + 2])
    #Generating distorted CRL, with voids
    opCRLdist = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0, arVoidCenCoordInCRL)
    print('done')
    AuxSaveOpTransmData(opCRLdist, 3, os.path.join(os.getcwd(), strDataFolderName, "res_opt_path_dif_dist_crl.dat"))

opApCRL = SRWLOptA('r', 'a', diamCRL, diamCRL) #Aperture at CRL

opDrCRL_CDI = SRWLOptD(0.5) #Drift space from CRL to a plane where a CDI image of imperfections could be observed
opDrCRL_Waist = SRWLOptD(9.7075) #Drift space from CRL to Waist (for 12.4 keV)
#opDrCRL_Waist = SRWLOptD(1.53323) #Drift space from CRL to Waist (for 5 keV)

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
prParRes1 = [0, 0, 1., 1, 0, 1.1, 8., 1.1, 8., 0, 0, 0]
prParRes0 = [0, 0, 1., 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

#"Beamline" - Container of Optical Elements (together with the corresponding wavefront propagation instructions)
optBLperf = SRWLOptC([opApCRL, opCRLperf, opDrCRL_Waist],
                     [prParRes1, prParRes0, prParRes0])
optBLdist = SRWLOptC([opApCRL, opCRLdist, opDrCRL_Waist],
                     [prParRes1, prParRes0, prParRes0])
optBLdistCDI = SRWLOptC([opApCRL, opCRLdist, opDrCRL_CDI],
                     [prParRes1, prParRes0, prParRes0])
#sys.exit(0)

#***********Wavefront Propagation
if(opCRLdist != None):
    #Duplicating wavefront (by re-creating all objects/arrays):
    wfr1 = copy.deepcopy(wfr)
    wfr2 = copy.deepcopy(wfr)

    print('Propagating Wavefront (through Distorted CRL and short Drift)...', end='')
    srwl.PropagElecField(wfr1, optBLdistCDI)
    print('done')
    print('Saving resulting Intensity and Phase data to files...', end='')
    arI = array('f', [0]*wfr1.mesh.nx*wfr1.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr1, 6, 0, 3, wfr1.mesh.eStart, 0, 0) #extracts intensity
    AuxSaveIntData(arI, wfr1.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_int_dist_crl_short_drift.dat"))
    arP = array('d', [0]*wfr1.mesh.nx*wfr1.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arP, wfr1, 0, 4, 3, wfr1.mesh.eStart, 0, 0) #extracts radiation phase
    AuxSaveIntData(arP, wfr1.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_phase_dist_crl_short_drift.dat"))
    del wfr1
    print('done')

    print('Propagating Wavefront (through Distorted CRL and Drift to waist)...', end='')
    srwl.PropagElecField(wfr2, optBLdist)
    print('done')
    print('Saving resulting Intensity and Phase data to files...', end='')
    arI = array('f', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr2, 6, 0, 3, wfr2.mesh.eStart, 0, 0) #extracts intensity
    AuxSaveIntData(arI, wfr2.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_int_dist_crl_waist.dat"))
    arP = array('d', [0]*wfr2.mesh.nx*wfr2.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arP, wfr2, 0, 4, 3, wfr2.mesh.eStart, 0, 0) #extracts radiation phase
    AuxSaveIntData(arP, wfr2.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_phase_dist_crl_waist.dat"))
    del wfr2
    print('done')

if(opCRLperf != None):
    print('Propagating Wavefront (through Perfect CRL and Drift to waist)...', end='')
    srwl.PropagElecField(wfr, optBLperf)
    print('done')
    print('Saving resulting data to files...', end='')
    arI = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
    AuxSaveIntData(arI, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_int_perf_crl_waist.dat"))
    arP = array('d', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arP, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0) #extracts radiation phase
    AuxSaveIntData(arP, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, "res_phase_perf_crl_waist.dat"))
    print('done')
