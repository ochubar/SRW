#############################################################################
# SRWLIB Example#18: Simulating Coherent X-ray (Gaussian beam) Scattering from from ensemble of 3D Nano-Particles randomly distributed over a volume
# Authors: H. Goel (SBU/ECE), O.C. (BNL/NSLS-II)
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

try: #OC15112022
    import sys
    sys.path.append('../')
    from srwlib import *
    from srwl_uti_smp import *
    import srwl_uti_smp_rnd_obj3d
    from uti_plot import *
except:
    from srwpy.srwlib import *
    from srwpy.srwl_uti_smp import *
    from srwpy import srwl_uti_smp_rnd_obj3d
    from srwpy.uti_plot import *
#from srwlib import *
#from srwl_uti_smp import *
#import srwl_uti_smp_rnd_obj3d
#from uti_plot import * #required for plotting
import copy
import os
import time

print('SRWLIB Python Example # 18:')
print('Simulating Coherent X-ray (Gaussian beam) Scattering on randomly-distibuted 3D spheres modeling an experimental sample (e.g. a colloidal solution with nano-objects).')
print('The propagation of X-ray beam through the sample is simulated in several steps (allowing to take into account multiple scattering).')

#***********Folder and Data File Names
strDataFolderName = 'data_example_18' #Data sub-folder name
strSampDefFolderName = 'samples' #Sub-folder name for storing sample data
strSampDefInFileName = 'smp_sph3d_list_%d.dat' #Sample definition file name core
strSampOptPathDifOutFileName = 'ex18_samp_opt_path_dif_%d.dat' #File name of Optical Path Difference corresponding to selected sample
strIntInitOutFileName = 'ex18_res_int_in.dat' #initial wavefront intensity distribution output file name
strIntPropOutFileName = 'ex18_res_int_prop.dat' #propagated wavefront intensity distribution output file name

#***********Gaussian Beam Source
GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Position of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = 8000 #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear horizontal
GsnBm.sigX = 7.e-06/2.35 #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

constConvRad = 1.23984186e-06/(4*3.1415926536)
rmsAngDiv = constConvRad/(GsnBm.avgPhotEn*GsnBm.sigX) #RMS angular divergence [rad]
print('RMS Source Size:', round(GsnBm.sigX*1.e+06, 3), 'um; RMS Divergence:', round(rmsAngDiv*1.e+06, 3), 'urad')

GsnBm.sigT = 10e-15 #Pulse duration [fs] (not used?)
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

#***********Initial Wavefront Parameters
wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, 100, 100) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 0.5 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]
wfr.unitElFld = 1 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

#Horizontal and Vertical Position Range for the Initial Wavefront calculation
firstHorAp = 40.e-06 #10.*sqrt((rmsAngDiv*wfr.mesh.zStart)**2 + (GsnBm.sigX)**2) #[m]
firstVertAp = firstHorAp #[m]
wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]

sampFactNxNyForProp = 3 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp

#***********Defining Sample, Beamline, Propagation Parameters

#Sample (/ Nano-Objects) Material Parameters (Au at 8 keV)
matDelta = 4.773e-05 #Refractive Index Decrement
matAttenLen = 2.48644e-06 #Attenuation Length [m]

nObj = 200 #Total number of Nano-Objects in the Sample
dAvgObj = 100e-09 #Average size of objects [m]
dSigObj = 25e-09 #Standard deviation of object sizes [m]

#Volume dimensions within which Sample is defined
rx = 10e-06
ry = 10e-06
rz = 100e-06

nSlices = 5 #10 #Number of "slices" of the Sample separated by small drifts on the beam propagation direction
nObjPart = int(nObj/nSlices) #Number of 3D Nano-Objects in one Slice
dRz = rz/nSlices #Small drift length (if nSlices > 1)

wfr.mesh.zStart -= 0.5*dRz*(nSlices - 1) #Slightly reducing the distance to the source for "thick" samples (nSlices > 1)

#Small Drift between "slices" (thin Transmission elements) for "thick" Sample
opSmpDr = SRWLOptD(dRz)

#Drift from Sample to Detector
lenDriftSmp_Det = 10. - 0.5*dRz*(nSlices - 1)
opSmp_Det = SRWLOptD(lenDriftSmp_Det)

#Wavefront Propagation Parameters:
#[0]: Auto-Resize (1) or not (0) Before propagation
#[1]: Auto-Resize (1) or not (0) After propagation
#[2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
#[3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
#[4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
#[5]: Horizontal Range modification factor at Resizing (1. means no modification)
#[6]: Horizontal Resolution modification factor at Resizing
#[7]: Vertical Range modification factor at Resizing
#[8]: Vertical Resolution modification factor at Resizing
#[9]: Type of wavefront Shift before Resizing
#[10]: New Horizontal wavefront Center position after Shift
#[11]: New Vertical wavefront Center position after Shift
#           [0][1][2] [3][4] [5] [6] [7]  [8]  [9][10][11] 
ppSmp =     [0, 0, 1., 0, 0, 1., 55., 1., 55.,  0, 0, 0]
ppSmpDr =   [0, 0, 1., 0, 0, 1., 1.,  1.,  1.,  0, 0, 0] #?
ppSmp_Det = [0, 0, 1., 3, 0, 1., 1.,  1.,  1.,  0, 0, 0]
ppFin =     [0, 0, 1., 0, 0, 1., 1.,  1.,  1.,  0, 0, 0]

arOpEl = [] #List of all Optical Elements (to be updated further)
arPP = [ppSmp] #List of the corresponding Propagation Parameters (to be updated further)

seed = 10 #Random generator seed number
zc = -0.5*dRz*(nSlices - 1) #Local longitudinal position of the first slice

print('   Generating Sample Slices and all Beamline ... ', end='')
t = time.time()

for i in range(nSlices):

    curSeed = seed if(i == 0) else None #For deterministic distribution of Nano-Particles

    #Defining list of Nano-Particles composing the Slice of the Sample 
    curSmpData = srwl_uti_smp_rnd_obj3d.setup_list_obj3d( #Initial list of 3D object (sphere) parameters
        _n = nObjPart, #Number of 3D nano-particles
        _ranges = [0.95*rx, 0.95*ry, dRz], #Ranges of horizontal, vertical and longitudinal position within which the 3D objects are defined
        _cen = [0., 0., zc], #Horizontal, Vertical and Longitudinal coordinates of center position around which the 3D objects are defined
        _dist = 'uniform', #Type (and eventual parameters) of distributions of 3D objects
        _obj_shape = ['S', 'gauss', dAvgObj, dSigObj], #Type of 3D objects, their distribution type and parameters (min. and max. diameter for the 'uniform' distribution)
        _allow_overlap = False, #Allow or not the 3D objects to overlap
        _seed = curSeed,
        _fp = os.path.join(os.getcwd(), strDataFolderName, strSampDefFolderName, strSampDefInFileName%(i)))

    #Defining thin Transmission object corresponding to the Slice
    curOpTr = srwl_opt_setup_transm_from_obj3d(shape_defs=curSmpData, delta=matDelta, atten_len=matAttenLen, rx=rx, ry=ry, nx=2000, ny=2000, xc=0, yc=0, extTr=1)

    #Extracting Optical Path Difference for the Slice and saving it to a file
    opPathDif = curOpTr.get_data(_typ = 3, _dep = 3)
    meshS = curOpTr.mesh
    srwl_uti_save_intens_ascii(opPathDif, meshS, os.path.join(os.getcwd(), strDataFolderName, strSampOptPathDifOutFileName%(i)), 0, 
        ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Optical Path Difference'], _arUnits=['eV', 'm', 'm', 'm'])

    #Creating Plots of Optical Path Difference in Sample Slices (to be shown at the end of the example)
    plotMeshSx = [meshS.xStart, meshS.xFin, meshS.nx]
    plotMeshSy = [meshS.yStart, meshS.yFin, meshS.ny]
    uti_plot2d(opPathDif, plotMeshSx, plotMeshSy, labels=['Horizontal Position', 'Vertical Position', 'Optical Path Diff. in Sample Slice #%d'%(i+1)], units=['m', 'm', 'm'])

    arOpEl.append(curOpTr) #Adding thin Transmission optical element corresponding to this Slice to the Container (/ BEamline)
    
    if(i > 0): arPP.append(ppSmpDr) #Same P.P. for Sample Slice and Drift bw Slices
   
    if(i < nSlices - 1): 
        arOpEl.append(opSmpDr) #Adding small Drift Space between Slices
        arPP.append(ppSmpDr) #Same P.P. for sample Slice and Drift bw Slices

    zc += dRz

print('done in', round(time.time() - t), 's')

#Adding Drift from Sample to Detector and corresponding Propagation Params
arOpEl.append(opSmp_Det)
arPP.append(ppSmp_Det)
arPP.append(ppFin) #Adding Final Resizing Params

opBL = SRWLOptC(arOpEl, arPP) #Defining Beamline / Container of Optical Elements

#***********Calculating Initial Wavefront of Gaussian Beam
#Calculating Initial Wavefront and extracting Intensity:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
mesh0 = deepcopy(wfr.mesh)
arI0 = array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, mesh0.eStart, 0, 0) #Extracting Intensity
srwl_uti_save_intens_ascii(
    arI0, mesh0, os.path.join(os.getcwd(), strDataFolderName, strIntInitOutFileName), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])

#***********Wavefront Propagation
print('   Propagating wavefront ... ', end='')
t = time.time()
srwl.PropagElecField(wfr, opBL)
print('done in', round(time.time() - t), 's')

print('   Extracting and projecting propagated wavefront intensity on detector and saving it to file ... ', end='')
t = time.time()
mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0) #Extracting Intensity

#***********Detector
nxDet = 2048 #Detector Number of Pixels in Horizontal direction
nyDet = 2048 #Detector Number of Pixels in Vertical direction
pSize = 75.e-06 #Detector Pizel Size
xrDet = (nxDet - 1)*pSize
yrDet = (nyDet - 1)*pSize
det = SRWLDet(_xStart = -0.5*xrDet, _xFin = 0.5*xrDet, _nx = nxDet, _yStart = -0.5*yrDet, _yFin = 0.5*yrDet, _ny = nyDet)

stkDet = det.treat_int(arI1, _mesh = mesh1) #"Projecting" Intensity on Detector (by interpolation)
mesh1 = stkDet.mesh; arI1 = stkDet.arS

srwl_uti_save_intens_ascii( #Saving "Projected" Intensity to a file
    arI1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Spectral Fluence'], _arUnits=['eV', 'm', 'm', 'J/eV/mm^2'])
print('done in', round(time.time() - t), 's')

#**********************Plotting Results (requires 3rd party graphics package)
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
plotMesh0x = [mesh0.xStart, mesh0.xFin, mesh0.nx]
plotMesh0y = [mesh0.yStart, mesh0.yFin, mesh0.ny]
uti_plot2d1d(arI0, plotMesh0x, plotMesh0y, x=0, y=0, labels=['Horizontal Position', 'Vertical Position', 'Intensity at Sample'], units=['m', 'm', 'ph/s/.1%bw/mm^2'])

plotMesh1x = [mesh1.xStart, mesh1.xFin, mesh1.nx]
plotMesh1y = [mesh1.yStart, mesh1.yFin, mesh1.ny]
arLogI1 = copy.copy(arI1)
nTot = mesh1.ne*mesh1.nx*mesh1.ny
for i in range(nTot):
    curI = arI1[i]
    if(curI <= 0.): arLogI1[i] = 0
    else: arLogI1[i] = log(curI, 10)
uti_plot2d1d(arLogI1, plotMesh1x, plotMesh1y, x=0, y=0, labels=['Horizontal Position', 'Vertical Position', 'Log of Intensity at Detector'], units=['m', 'm', ''])

uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')
