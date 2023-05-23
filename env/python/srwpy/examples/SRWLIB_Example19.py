#############################################################################
# SRWLIB Example#19: Simulating Coherent X-ray (Gaussian beam) Scattering from ensemble
# of 3D Nano-Particles modifying their positions due to Brownian motion
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
import os
import time

npIsAvail = False
try:
    import numpy as np
    npIsAvail = True
except:
    print('NumPy can not be loaded. You may need to install numpy, otherwise some functionality of this example will not be available. If you are using pip, you can use the following command to install it: \npip install numpy')

print('SRWLIB Python Example # 19:')
print('Simulating Coherent X-ray (Gaussian beam) Scattering from ensemble of 3D particles modifying their positions due to Brownian motion')

#***********Folder and Data File Names
strDataFolderName = 'data_example_19' #Data sub-folder name
strSampleSubFolderName = 'samples' #Sub-folder name for storing Sample data
strListSampObjFileName = 'ex19_smp_obj_list_%d.dat' #List of 3D Nano-Objects / Sample file name
strSampOptPathDifOutFileName = 'ex19_smp_opt_path_dif_%d.dat' #optical path difference corresponding to selected sample
strIntInitOutFileName = 'ex19_res_int_in.dat' #initial wavefront intensity distribution output file name
strIntPropOutFileName = 'ex19_res_int_prop_%d.dat' #propagated wavefront intensity distribution output file name
strIntPropOutFileNameDet = 'ex19_res_int_det.h5' #intensity distribution regisgtered by detector output file name

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

#***********Initial Wavefront
wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, 100, 100) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 0.5 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]
wfr.unitElFld = 1 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

#Horizontal and Vertical Position Range for the Initial Wavefront calculation
firstHorAp = 40.e-06 #[m]
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

#Calculate Initial Wavefront, extract Intensity and save it to a file:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
mesh0 = deepcopy(wfr.mesh)
arI0 = array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, mesh0.eStart, 0, 0) #Extract intensity
srwl_uti_save_intens_ascii( #Save Intensity of to a file
    arI0, mesh0, os.path.join(os.getcwd(), strDataFolderName, strIntInitOutFileName), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])

#Plot the Initial Wavefront (without showing it yet)
plotMesh0x = [mesh0.xStart, mesh0.xFin, mesh0.nx]
plotMesh0y = [mesh0.yStart, mesh0.yFin, mesh0.ny]
uti_plot2d1d(arI0, plotMesh0x, plotMesh0y, 0, 0, ['Horizontal Position', 'Vertical Position', 'Intensity at Sample'], ['m', 'm', 'ph/s/.1%bw/mm^2'])

#************Defining Samples (lists of 3D objects (spheres))
#Initial set of 3D objects
rx = 20.e-06 #Range of Horizontal position [m] within which 3D Objects constituing Sample are defined
ry = 20.e-06 #Range of Vertical position [m]
rz = 20.e-06 #Range of Longitudinal position [m]
xc = 0 #Horizontal Center position of the Sample
yc = 0 #Vertical Center position of the Sample
zc = 0 #Longitudinal Center position of the Sample

listObjInit = srwl_uti_smp_rnd_obj3d.setup_list_obj3d( #Initial list of 3D object (sphere) parameters
    _n = 100, #Number of 3D nano-objects
    _ranges = [0.95*rx, 0.95*ry, rz], #Ranges of horizontal, vertical and longitudinal position within which the 3D objects are defined
    #_ranges = [rx, ry, rz], #Ranges of horizontal, vertical and longitudinal position within which the 3D objects are defined
    _cen = [xc, yc, zc], #Horizontal, Vertical and Longitudinal coordinates of center position around which the 3D objects are defined
    _dist = 'uniform', #Type (and eventual parameters) of distributions of 3D objects
    _obj_shape = ['S', 'uniform', 25.e-09, 250.e-09], #Type of 3D objects, their distribution type and parameters (min. and max. diameter for the 'uniform' distribution)
    _allow_overlap = False, #Allow or not the 3D objects to overlap
    _fp = os.path.join(os.getcwd(), strDataFolderName, strSampleSubFolderName, strListSampObjFileName%(0)))
    
#Generate timesteps of Brownian motion of the 3D nano-objects (spheres) simulating particles suspended in water at room temperature
timeStep = 0.1 #Time step between different Sample "snapshots" / scattering patterns
timeInterv = 0.5 #Total time interval covered by the "snapshots"
listObjBrownian = srwl_uti_smp_rnd_obj3d.brownian_motion3d(
    _obj_crd = listObjInit, #Initial list of 3D objects
    _viscosity = 1.e-3, #[Pa*s]
    _temperature = 293, #[K]
    _timestep = timeStep, #[s]
    _duration = timeInterv, #[s]
    _seed = 0,
    _fp = os.path.join(os.getcwd(), strDataFolderName, strSampleSubFolderName, strListSampObjFileName))

#Sample Material Characteristics (Au at 8 keV)
matDelta = 4.773e-05 #Refractive Index Decrement
matAttenLen = 2.48644e-06 #Attenuation Length [m]

#***********Detector
nxDet = 2048 #Detector Number of Pixels in Horizontal direction
nyDet = 2048 #Detector Number of Pixels in Vertical direction
pSize = 75e-06 #Detector Pixel Size
xrDet = nxDet*pSize
yrDet = nyDet*pSize
det = SRWLDet(_xStart = -0.5*xrDet, _xFin = 0.5*xrDet, _nx = nxDet, _yStart = -0.5*yrDet, _yFin = 0.5*yrDet, _ny = nyDet, _eStart=GsnBm.avgPhotEn, _eFin=GsnBm.avgPhotEn) #OC20092021
#det = SRWLDet(_xStart = -0.5*xrDet, _xFin = 0.5*xrDet, _nx = nxDet, _yStart = -0.5*yrDet, _yFin = 0.5*yrDet, _ny = nyDet)

arDetFrames = None #Array to store all detector frames data
if(npIsAvail): arDetFrames = np.zeros((len(listObjBrownian), nxDet, nyDet))

#***********Defining Drift from Sample to Detector and Propagation Parameters
distSmp_Det = 10.
opSmp_Det = SRWLOptD(distSmp_Det)

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
ppSmp_Det = [0, 0, 1., 3, 0, 1., 1.,  1.,  1.,  0, 0, 0]
ppFin =     [0, 0, 1., 0, 0, 1., 1.,  1.,  1.,  0, 0, 0]

#***********Wavefront Propagation / Scattering calculation for different instances of Sample created by Brownnian motion
for it in range(len(listObjBrownian)):

    if(it == 0): print('   Performing Simulaitons for the Initial Nano-Particle Distribution ***********')
    else: print('   Performing Simulaitons for Brownian Motion Step:', it, '***********')

    print('   Setting up Transmission optical element from input Sample data ... ', end='')
    t = time.time()
    opSmp = srwl_opt_setup_transm_from_obj3d( #Defining Sample (Transmission object)
        shape_defs = listObjBrownian[it], #List of 3D Nano-Object params for the current step
        delta = matDelta, atten_len = matAttenLen, #3D Nano-Object Material params
        rx = rx, ry = ry, #Range of Horizontal and Vertical position [m] within which Nano-Objects constituing the Sample are defined
        nx = 2000, ny = 2000, #Numbers of points vs Horizontal and Vertical position for the Transmission
        xc = xc, yc = yc, #Horizontal and Vertical Center positions of the Sample
        extTr = 1) #Transmission outside the grid/mesh is zero (0), or the same as on boundary (1)
    print('done in', round(time.time() - t, 3), 's')

    print('   Extracting Optical Path Difference data from Sample Transmission optical element ... ', end='')
    t = time.time()
    opPathDif = opSmp.get_data(_typ = 3, _dep = 3)

    print('done in', round(time.time() - t, 3), 's')

    print('   Saving Optical Path Difference data from Sample Transmission optical element ... ', end='')
    t = time.time()
    srwl_uti_save_intens_ascii(
        opPathDif, opSmp.mesh, os.path.join(os.getcwd(), strDataFolderName, strSampOptPathDifOutFileName%(it)), 0,
        ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Optical Path Difference'], _arUnits=['eV', 'm', 'm', 'm'])
    print('done in', round(time.time() - t, 3), 's')

    #Defining "Beamline" to Propagate the Wavefront through
    opBL = SRWLOptC([opSmp, opSmp_Det], 
                    [ppSmp, ppSmp_Det, ppFin])

    #Duplicating Initial Wavefront to perform its Propagaton
    wfrP = deepcopy(wfr)

    print('   Propagating Wavefront ... ', end='')
    t = time.time()
    srwl.PropagElecField(wfrP, opBL)
    print('done in', round(time.time() - t), 's')

    print('   Extracting, Projecting the Propagated Wavefront Intensity on Detector and Saving it to file ... ', end='')
    t = time.time()
    mesh1 = deepcopy(wfrP.mesh)
    arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI1, wfrP, 6, 0, 3, mesh1.eStart, 0, 0) #extracts intensity

    stkDet = det.treat_int(arI1, _mesh = mesh1) #"Projecting" intensity on detector (by interpolation)
    mesh1 = stkDet.mesh
    arI1 = stkDet.arS
    srwl_uti_save_intens_ascii(
        arI1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName%(it)), 0,
        ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Spectral Fluence'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])

    if(arDetFrames is not None): arDetFrames[it] = np.reshape(arI1, (mesh1.ny, mesh1.nx)).transpose()
    print('done in', round(time.time() - t), 's')

    #Plotting the Results (requires 3rd party graphics package)
    print('   Plotting the results (i.e. creating plots without showing them yet) ... ', end='')

    #Sample Optical Path Diff.
    meshS = opSmp.mesh
    plotMeshSx = [meshS.xStart, meshS.xFin, meshS.nx]
    plotMeshSy = [meshS.yStart, meshS.yFin, meshS.ny]
    uti_plot2d(opPathDif, plotMeshSx, plotMeshSy, ['Horizontal Position', 'Vertical Position', 'Optical Path Diff. in Sample (Time = %.3fs)' % (it*timeStep)], ['m', 'm', 'm'])
        
    #Scattered Radiation Intensity Distribution in Log Scale
    plotMesh1x = [mesh1.xStart, mesh1.xFin, mesh1.nx]
    plotMesh1y = [mesh1.yStart, mesh1.yFin, mesh1.ny]
    arLogI1 = copy(arI1)
    nTot = mesh1.ne*mesh1.nx*mesh1.ny
    for i in range(nTot):
        curI = arI1[i]
        if(curI <= 0.): arLogI1[i] = 0 #?
        else: arLogI1[i] = log(curI, 10)

    uti_plot2d1d(arLogI1, plotMesh1x, plotMesh1y, 0, 0, ['Horizontal Position', 'Vertical Position', 'Log of Intensity at Detector (Time = %.3f s)' % (it*timeStep)], ['m', 'm', ''])

    print('done')

if(arDetFrames is not None): #Saving simulated Detector data file
    print('   Saving all Detector data to another file (that can be used in subsequent processing) ... ', end='')
    srwl_uti_save_intens_hdf5_exp(arDetFrames, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileNameDet), 
        _exp_type = 'XPCS', _dt = timeStep, _dist_smp = distSmp_Det, _bm_size_x = GsnBm.sigX*2.35, _bm_size_y = GsnBm.sigY*2.35)
    print('done')

uti_plot_show() #Show all plots created
