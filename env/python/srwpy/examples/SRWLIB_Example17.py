#############################################################################
# SRWLIB Example#17: Simulating Coherent X-ray (Gaussian beam) Scattering
# on Experimental Samples defined using different methods
# Authors: O.C., M. Rakitin, R. Coles (BNL, NSLS-II)
# SEM images of nano-fabricated samples from J. Lhermitte, K. Yager (BNL, CFN)
# v 0.02
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

try: #OC15112022
    import sys
    sys.path.append('../')
    from srwlib import *
    from srwl_uti_smp import *
    from uti_plot import *
except:
    from srwpy.srwlib import *
    from srwpy.srwl_uti_smp import *
    from srwpy.uti_plot import *
#from srwlib import *
#from srwl_uti_smp import *
#from uti_plot import * #required for plotting
import copy
import os
import time

print('SRWLIB Python Example # 17:')
print('Simulating Coherent X-ray (Gaussian beam) Scattering on Experimental Sample modelled by 2D (nano-)objects')

#**********************Input Parameters and Structures
#***********Folder and Data File Names
strDataFolderName = 'data_example_17' #data sub-folder name
strSampImgFolderName = 'samples'
strSampImgInFileName01 = 'R5.tif' #SEM image of nano-fabricated sample: Rings
strSampImgInFileName02 = 'H5.tif' #SEM image of nano-fabricated sample: Dots
strSampImgInFileName03 = 'H5R5.tif' #SEM image of nano-fabricated sample: Rings and Dots
strSampImgOutFileName01 = 'ex17_samp_img_proc.tif' #Processes (output) sample file name
strSampOptPathDifOutFileName01 = 'ex17_samp_opt_path_dif.dat' #Optical path difference corresponding to selected sample
strIntInitOutFileName01 = 'ex17_res_int_in.dat' #initial wavefront intensity distribution output file name
strIntPropOutFileName01 = 'ex17_res_int_prop.dat' #propagated wavefront intensity distribution output file name

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

#***********Detector
npx = 2070 #Detector Number of Pixels in Horizontal direction
npy = 2167 #Detector Number of Pixels in Vertical direction
pSize = 75e-06 #Detector Pizel Size
xrDet = npx*pSize
yrDet = npy*pSize
det = SRWLDet(_xStart = -0.5*xrDet, _xFin = 0.5*xrDet, _nx = npx, _yStart = -0.5*yrDet, _yFin = 0.5*yrDet, _ny = npy)

#***********Defining Sample and Propagation Parameters
print('   Setting up Transmission optical element from input Sample data ... ', end='')
t = time.time()

smpType = 'img' #to choose between 2 options of sample definition
#smpType = 'rnd' #to choose

#Object thickness along X-ray propagation direction
objThickn = 50.e-09
#Object material characteristics (Au at 8 keV)
matDelta = 4.773e-05 #Refractive Index Decrement
matAttenLen = 2.48644e-06 #Attenuation Length [m]

opSmp = None; imgProc = None
if(smpType == 'img'): #from SEM image
    opSmp, imgProc = srwl_opt_setup_transm_from_file(
        file_path = os.path.join(os.getcwd(), strDataFolderName, strSampImgFolderName,
        #Select (uncomment) one of 3 following lines, coresponding to 3 different images:
        strSampImgInFileName01), resolution = 2.481e-09, area = [100, 1179, 130, 830], #useful area of original image in pixels: [x_start, x_end, y_start, y_end]
        #strSampImgInFileName02), resolution = 1.416e-09, area = [0, 1279, 160, 835],
        #strSampImgInFileName03), resolution = 4.95e-09, area = [150, 1279, 160, 800],
        thickness = objThickn, delta = matDelta, atten_len = matAttenLen, 
        extTr = 1, xc = 0, yc = 0, shift_x = 0, shift_y = 0, 
        rotate_angle = None, invert = None,
        cutoff_background_noise = 0.4, background_color = 0,
        _cutoff_max_fact = 0.5, _max_color = 255,
        _ret = 'all')
    
elif(smpType == 'rnd'):
    opSmp, imgProc = srwl_opt_setup_smp_rnd_obj2d(
        _thickness = objThickn, _delta = matDelta, _atten_len = matAttenLen,
        _rx = 10.e-06, _ry = 10.e-06, _xc = 0, _yc = 0, _nx = 4001, _ny = 4001,
        _dens = 20.e+06, _edge_frac = 0.02,
        _obj_type = 2, #1- recatngles, 2- ellipses, 3- triangles, 4- polygons, 5- mixed shapes
        _r_min_bw_obj = 1.e-09, _obj_size_min = 100e-09, _obj_size_max = 120e-09, _size_dist = 1, 
        _ang_min = 30., _ang_max = 60., _ang_dist = 1, _rand_alg = 1, #2, _sim_step_size = 0.5e-06,
        _obj_par1 = 0.5, _obj_par2 = None,
        _ext_tr = 1, _ret='all')

print('done in', round(time.time() - t, 3), 's')

print('   Saving auxiliary processed image file ... ', end='')
t = time.time()
#Saving Processed Image for eventual tests
imgProc.save(os.path.join(os.getcwd(), strDataFolderName, strSampImgOutFileName01))
print('done in', round(time.time() - t, 3), 's')
print('   Extracting Optical Path Difference data from Sample Transmission optical element ... ', end='')
t = time.time()
#Extracting and saving Optical Path Difference for eventual tests
opPathDif = opSmp.get_data(_typ = 3, _dep = 3)
print('done in', round(time.time() - t, 3), 's')
print('   Saving Optical Path Difference data from Sample Transmission optical element ... ', end='')
t = time.time()
srwl_uti_save_intens_ascii(
    opPathDif, opSmp.mesh, os.path.join(os.getcwd(), strDataFolderName, strSampOptPathDifOutFileName01), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Optical Path Difference'], _arUnits=['eV', 'm', 'm', 'm'])
print('done in', round(time.time() - t, 3), 's')

#Drift from Sample to Detector
opSmp_Det = SRWLOptD(10.)

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

opBL = SRWLOptC([opSmp, opSmp_Det],
                [ppSmp, ppSmp_Det, ppFin])

#**********************Main Calculations
#***********Initial Wavefront of Gaussian Beam
#Initial Wavefront and extracting Intensity:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
mesh0 = deepcopy(wfr.mesh)
arI0 = array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, mesh0.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(
    arI0, mesh0, os.path.join(os.getcwd(), strDataFolderName, strIntInitOutFileName01), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])

#***********Wavefront Propagation
print('   Propagating wavefront ... ', end='')
tryUsingGPU = 1 #0 #Set to 1 if GPU should be used, 0 otherwise #OC21032024
if(tryUsingGPU): print('trying to use GPU ... ', end='')
t = time.time()
srwl.PropagElecField(wfr, opBL, None, tryUsingGPU)
#srwl.PropagElecField(wfr, opBL)
print('done in', round(time.time() - t), 's')

print('   Extracting, projecting propagated wavefront intensity on detector and saving it to file ... ', end='')
t = time.time()
mesh1 = deepcopy(wfr.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0, None, None, tryUsingGPU) #extracts intensity (eventually using GPU)
#srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0) #extracts intensity

stkDet = det.treat_int(arI1, _mesh = mesh1) #"Projecting" intensity on detector (by interpolation)
mesh1 = stkDet.mesh; arI1 = stkDet.arS
srwl_uti_save_intens_ascii(
    arI1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName01), 0,
    ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Spectral Fluence'], _arUnits=['eV', 'm', 'm', 'J/eV/mm^2'])
print('done in', round(time.time() - t), 's')

#**********************Plotting Results (requires 3rd party graphics package)
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
plotMesh0x = [mesh0.xStart, mesh0.xFin, mesh0.nx]
plotMesh0y = [mesh0.yStart, mesh0.yFin, mesh0.ny]
uti_plot2d1d(arI0, plotMesh0x, plotMesh0y, x=0, y=0, labels=['Horizontal Position', 'Vertical Position', 'Intensity at Sample'], units=['m', 'm', 'ph/s/.1%bw/mm^2'])

uti_plot2d(imgProc, labels=['Horizontal Pixel #', 'Vertical Pixel #', 'Processed Image of Sample'])

meshS = opSmp.mesh
plotMeshSx = [meshS.xStart, meshS.xFin, meshS.nx]
plotMeshSy = [meshS.yStart, meshS.yFin, meshS.ny]
uti_plot2d1d(opPathDif, plotMeshSx, plotMeshSy, x=0, y=0, labels=['Horizontal Position', 'Vertical Position', 'Optical Path Diff. in Sample'], units=['m', 'm', 'm'])

plotMesh1x = [mesh1.xStart, mesh1.xFin, mesh1.nx]
plotMesh1y = [mesh1.yStart, mesh1.yFin, mesh1.ny]
arLogI1 = copy.copy(arI1)
nTot = mesh1.ne*mesh1.nx*mesh1.ny
for i in range(nTot):
    curI = arI1[i]
    if(curI <= 0.): arLogI1[i] = 0
    else: arLogI1[i] = log(curI)
uti_plot2d1d(arLogI1, plotMesh1x, plotMesh1y, x=0, y=0, labels=['Horizontal Position', 'Vertical Position', 'Log of Intensity at Detector'], units=['m', 'm', ''])

uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')

