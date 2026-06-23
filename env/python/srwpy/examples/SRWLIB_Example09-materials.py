#############################################################################
# SRWLIB Example#9_Materials: Simulating propagation of a Gaussian X-ray beam
# through a Beamline Containing Mirrors made of different materials
# v 0.04
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

try: #OC15112022
    import sys
    sys.path.append('../')
    from srwlib import *
    from uti_plot import *
    from uti_mtrl import *
except:
    from srwpy.srwlib import *
    from srwpy.uti_plot import *
    from srwpy.uti_mtrl import *
import os
import copy

print('SRWLIB Python Example # 9 with material properties:')
print('Simulating propagation of a coherent Gaussian X-ray beam through Pt and B4C mirrors')



#******************* Reflectivity Parameters

_n_ph_en=100
_n_ang=1000
_n_comp=2
_ph_en_start=12000
_ph_en_fin=13000
_ph_en_scale_type="lin"
_ang_start=0.001
_ang_fin=0.01
_ang_scale_type="lin"
avgEn = 12400


Pt_refl=calc_refl_arr('Pt', _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                  _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)
B4C_refl=calc_refl_arr('B4C', _dens=2.5, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                  _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)



#**********************Input Parameters and Structures
#***********Folder and Data File Names
strDataFolderName = 'data_example_09-materials' #data sub-folder name
strIntOutFileName = 'ex09_res_int_in.dat' #initial wavefront intensity distribution output file name
strPhOutFileName = 'ex09_res_phase_in.dat' #initial wavefront phase output file name
strIntPropOutFileName = 'ex09_res_int_prop.dat' #propagated wavefront intensity distribution output file name
strPhasePropOutFileName  = 'ex09_res_phase_prop.dat' #propagated wavefront phase output file name
strIntPropOutFileName2 = 'ex09_res_int_prop_Pt.dat' #propagated wavefront intensity distribution output file name
strPhPropOutFileName2  = 'ex09_res_phase_prop_Pt.dat' #propagated wavefront phase output file name
strIntPropOutFileName3 = 'ex09_res_int_prop_B4C.dat' #propagated wavefront intensity distribution output file name
strPhPropOutFileName3  = 'ex09_res_phase_prop_B4C.dat' #propagated wavefront phase output file name


#***********Gaussian Beam Source
GsnBm = SRWLGsnBm() #Gaussian Beam structure (just parameters)
GsnBm.x = 0 #Transverse Positions of Gaussian Beam Center at Waist [m]
GsnBm.y = 0
GsnBm.z = 0 #Longitudinal Position of Waist [m]
GsnBm.xp = 0 #Average Angles of Gaussian Beam at Waist [rad]
GsnBm.yp = 0
GsnBm.avgPhotEn = avgEn #Photon Energy [eV]
GsnBm.pulseEn = 0.001 #Energy per Pulse [J] - to be corrected
GsnBm.repRate = 1 #Rep. Rate [Hz] - to be corrected
GsnBm.polar = 1 #1- linear horizontal
GsnBm.sigX = 23e-06/2.35 #Horiz. RMS size at Waist [m]
GsnBm.sigY = GsnBm.sigX #Vert. RMS size at Waist [m]

constConvRad = 1.23984186e-06/(4*3.1415926536)
rmsAngDiv = constConvRad/(GsnBm.avgPhotEn*GsnBm.sigX) #RMS angular divergence [rad]
print('RMS Source Size:', round(GsnBm.sigX*1.e+06, 3), 'micro-m; RMS Divergence:', round(rmsAngDiv*1.e+06, 3), 'micro-rad')

GsnBm.sigT = 10e-15 #Pulse duration [fs] (not used?)
GsnBm.mx = 0 #Transverse Gauss-Hermite Mode Orders
GsnBm.my = 0

#***********Initial Wavefront
wfr = SRWLWfr() #Initial Electric Field Wavefront
wfr.allocate(1, 100, 100) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
wfr.mesh.zStart = 25 #Longitudinal Position [m] at which initial Electric Field has to be calculated, i.e. the position of the first optical element
wfr.mesh.eStart = GsnBm.avgPhotEn #Initial Photon Energy [eV]
wfr.mesh.eFin = GsnBm.avgPhotEn #Final Photon Energy [eV]

wfr.unitElFld = 2 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

distSrc_VFM = wfr.mesh.zStart - GsnBm.z
#Horizontal and Vertical Position Range for the Initial Wavefront calculation
#can be used to simulate the First Aperture
firstHorAp = 8.*rmsAngDiv*distSrc_VFM #[m]
firstVertAp = firstHorAp #[m]

wfr.mesh.xStart = -0.5*firstHorAp #Initial Horizontal Position [m]
wfr.mesh.xFin = 0.5*firstHorAp #Final Horizontal Position [m]
wfr.mesh.yStart = -0.5*firstVertAp #Initial Vertical Position [m]
wfr.mesh.yFin = 0.5*firstVertAp #Final Vertical Position [m]s

sampFactNxNyForProp = 4 #sampling factor for adjusting nx, ny (effective if > 0)
arPrecPar = [sampFactNxNyForProp]

wfr.partBeam.partStatMom1.x = GsnBm.x #Some information about the source in the Wavefront structure
wfr.partBeam.partStatMom1.y = GsnBm.y
wfr.partBeam.partStatMom1.z = GsnBm.z
wfr.partBeam.partStatMom1.xp = GsnBm.xp
wfr.partBeam.partStatMom1.yp = GsnBm.yp


#**********************Calculating Initial Wavefront and extracting Intensity:
srwl.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
arI0 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI0, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, strIntOutFileName), 0)

arI0x = array('f', [0]*wfr.mesh.nx) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI0x, wfr, 6, 0, 1, wfr.mesh.eStart, 0, 0) #extracts intensity

arP0 = array('d', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP0, wfr, 0, 4, 3, wfr.mesh.eStart, 0, 0) #extracts radiation phase
srwl_uti_save_intens_ascii(arP0, wfr.mesh, os.path.join(os.getcwd(), strDataFolderName, strPhOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Phase'], _arUnits=['', 'm', 'm', 'rad'])

mesh0 = deepcopy(wfr.mesh)



#***********Optical Elements and Propagation Parameters
#Sequence of Optical Elements:
#   <Drift to VFM>
#   <Aperture of KB>
#   <VFM>
#   <Drift to HFM>
#   <HFM>
#   <Drift to Sample>



distVFM_HFM = 1 #Distance from VFM to HFM [m]
distHFM_Samp = 1 #Distance from HFM to Sample [m]


lenKB = 0.1 #Length of VFM and HFM (same for each) [m]
angKB = 0.004 #b4c_ca #2.521e-3 #Incident Angle of VFM and HFM [rad]



#Aperture of KB system
opApKB = SRWLOptA('r', 'a', lenKB*angKB, lenKB*angKB)


#VFM simulated by Ideal Lens:
#opVFM = SRWLOptL(_Fy=(distSrc_M1 + distM1_VFM)*(distVFM_HFM + distHFM_Samp)/(distSrc_M1 + distM1_VFM + distVFM_HFM + distHFM_Samp))
#VFM simulated by Extended Elliptical Mirror:
opVFM = SRWLOptMirEl(_p=(distSrc_VFM), _q=(distVFM_HFM + distHFM_Samp), _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=0, _nvy=cos(angKB), _nvz=-sin(angKB), _tvx=0, _tvy=-sin(angKB))

#Drift from VFM to HFM
opDrVFM_HFM = SRWLOptD(distVFM_HFM)

#VFM simulated by Extended Elliptical Mirror:
opHFM = SRWLOptMirEl(_p=(distSrc_VFM + distVFM_HFM), _q=distHFM_Samp, _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=cos(angKB), _nvy=0, _nvz=-sin(angKB), _tvx=-sin(angKB), _tvy=0)

#Drift from HFM to Sample
opDrHFM_Samp = SRWLOptD(distHFM_Samp)


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
#[9]: Type of wavefront Shift before Resizing (not yet implemented)
#[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
#[11]: New Vertical wavefront Center position after Shift (not yet implemented)
#           [0][1][2] [3][4] [5]   [6] [7]  [8] [9][10][11]
#prParInit = [0, 0, 1., 1, 0, 2.,   5., 2.,  3.,  0, 0, 0]
prParInit = [0, 0, 1., 1, 0, 2.,   2., 2.,   2.,  0, 0, 0]
prPar0 =    [0, 0, 1., 1, 0, 1.,   2,  1.,   2,  0, 0, 0]
prPar1 =    [0, 0, 1., 1, 0, 1,   2,  0.1,   2,  0, 0, 0]
#prParPost = [0, 0, 1., 1, 0, 0.06, 3., 0.1, 2.,  0, 0, 0]

#NOTE: in this case of simulation, it can be enough to define the precision parameters only Before and After
#the propagation through the entire Beamline. However, if necessary, different propagation parameters can be specified
#for each optical element.
#The optimal values of propagation parameters may depend on photon energy and optical layout.

#"Beamline" - a sequenced Container of Optical Elements (together with the corresponding wavefront propagation parameters,
#and the "post-propagation" wavefront resizing parameters for better viewing).
# optBL = SRWLOptC([opVFM, opTrErVFM, opDrVFM_HFM, opHFM, opTrErHFM, opDrHFM_Samp],
#                  [prPar0, prPar0,   prPar0,     prPar0,   prPar0,  prParPost])
optBL = SRWLOptC([opApKB, opVFM, opDrVFM_HFM, opHFM, opDrHFM_Samp],
                 [prParInit, prPar0, prPar0, prPar0, prPar1])


# # Beamline with material properties
# #VFM simulated by Ideal Lens:
# #opVFM = SRWLOptL(_Fy=(distSrc_M1 + distM1_VFM)*(distVFM_HFM + distHFM_Samp)/(distSrc_M1 + distM1_VFM + distVFM_HFM + distHFM_Samp))
# #VFM simulated by Extended Elliptical Mirror:
opVFM_Pt = SRWLOptMirEl(_p=(distSrc_VFM), _q=(distVFM_HFM + distHFM_Samp), _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=0, _nvy=cos(angKB), _nvz=-sin(angKB), _tvx=0, _tvy=-sin(angKB),
                     _refl=Pt_refl, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                     _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)
# # opTrErVFM_Pt = srwl_opt_setup_surf_height_1d(heightProfData, _dim='y', _ang=angKB, _amp_coef=1)

# #HFM simulated by Ideal Lens:
# #opHFM = SRWLOptL(_Fx=(distSrc_M1 + distM1_VFM + distVFM_HFM)*distHFM_Samp/(distSrc_M1 + distM1_VFM + distVFM_HFM + distHFM_Samp))
# #VFM simulated by Extended Elliptical Mirror:
opHFM_Pt = SRWLOptMirEl(_p=(distSrc_VFM + distVFM_HFM), _q=distHFM_Samp, _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=cos(angKB), _nvy=0, _nvz=-sin(angKB), _tvx=-sin(angKB), _tvy=0,
                     _refl=Pt_refl, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                     _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)
# opTrErHFM_Pt = srwl_opt_setup_surf_height_1d(heightProfData, _dim='x', _ang=angKB, _amp_coef=1)

# #"Beamline" - a sequenced Container of Optical Elements (together with the corresponding wavefront propagation parameters,
# #and the "post-propagation" wavefront resizing parameters for better viewing).
# # optBL_Pt = SRWLOptC([opVFM_Pt, opTrErVFM_Pt, opDrVFM_HFM, opHFM_Pt, opTrErHFM_Pt, opDrHFM_Samp],
# #                     [prPar0,    prPar0,     prPar0,         prPar0,   prPar0,  prParPost])
optBL_Pt = SRWLOptC([opApKB, opVFM_Pt, opDrVFM_HFM, opHFM_Pt, opDrHFM_Samp],
                 [prParInit, prPar0, prPar0, prPar0, prPar1])



# Boron Carbide Mirrors
opVFM_B4C = SRWLOptMirEl(_p=(distSrc_VFM), _q=(distVFM_HFM + distHFM_Samp), _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=0, _nvy=cos(angKB), _nvz=-sin(angKB), _tvx=0, _tvy=-sin(angKB),
                     _refl=B4C_refl, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                     _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)

opHFM_B4C = SRWLOptMirEl(_p=(distSrc_VFM + distVFM_HFM), _q=distHFM_Samp, _ang_graz=angKB, _size_tang=lenKB, _size_sag=10.e-03,
                     _nvx=cos(angKB), _nvy=0, _nvz=-sin(angKB), _tvx=-sin(angKB), _tvy=0,
                     _refl=B4C_refl, _n_ang=_n_ang, _n_ph_en=_n_ph_en, _n_comp=_n_comp, _ph_en_start=_ph_en_start, _ph_en_fin=_ph_en_fin, _ph_en_scale_type=_ph_en_scale_type,
                     _ang_start=_ang_start, _ang_fin=_ang_fin, _ang_scale_type=_ang_scale_type)

optBL_B4C = SRWLOptC([opApKB, opVFM_B4C, opDrVFM_HFM, opHFM_B4C, opDrHFM_Samp],
                 [prParInit, prPar0, prPar0, prPar0, prPar1])




#***********Wavefront Propagation
arI1 = None; arI1x = None; arI1y = None
arP1 = None
mesh1 = None


#Duplicating wavefront (by re-creating all objects/arrays):
wfr1 = copy.deepcopy(wfr)
wfr2 = copy.deepcopy(wfr)
wfr3 = copy.deepcopy(wfr)



print('   Propagating Wavefront (through two elliptical mirrors)... ', end='')
srwl.PropagElecField(wfr1, optBL)
print('done')
print('   Saving resulting Intensity and Phase data to files ... ', end='')
mesh1 = deepcopy(wfr1.mesh)
arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 3, mesh1.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName), 0)

arI1x = array('f', [0]*mesh1.nx) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI1x, wfr1, 6, 0, 1, mesh1.eStart, 0, 0) #extracts intensity
arI1y = array('f', [0]*mesh1.ny) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI1y, wfr1, 6, 0, 2, mesh1.eStart, 0, 0) #extracts intensity
arP1 = array('d', [0]*mesh1.nx*mesh1.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP1, wfr1, 0, 4, 3, mesh1.eStart, 0, 0) #extracts radiation phase
srwl_uti_save_intens_ascii(arP1, mesh1, os.path.join(os.getcwd(), strDataFolderName, strPhasePropOutFileName), 0, ['', 'Horizontal Position', 'Vertical Position', 'Phase'], _arUnits=['', 'm', 'm', 'rad'])
del wfr1
print('done')


print('   Propagating Wavefront (through 2 Platinum Mirrors) ... ', end='')
srwl.PropagElecField(wfr2, optBL_Pt)
print('done')
print('   Saving resulting Intensity and Phase data to files ... ', end='')
mesh2 = deepcopy(wfr2.mesh)
arI2 = array('f', [0]*mesh2.nx*mesh2.ny) #"flat" array to take 2D intensity data
srwl.CalcIntFromElecField(arI2, wfr2, 6, 0, 3, mesh2.eStart, 0, 0) #extracts intensity
srwl_uti_save_intens_ascii(arI2, mesh2, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName2), 0)

arI2x = array('f', [0]*mesh2.nx) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI2x, wfr2, 6, 0, 1, mesh2.eStart, 0, 0) #extracts intensity
arI2y = array('f', [0]*mesh2.ny) #array to take 1D intensity data
srwl.CalcIntFromElecField(arI2y, wfr2, 6, 0, 2, mesh2.eStart, 0, 0) #extracts intensity
arP2 = array('d', [0]*mesh2.nx*mesh2.ny) #"flat" array to take 2D phase data (note it should be 'd')
srwl.CalcIntFromElecField(arP2, wfr2, 0, 4, 3, mesh2.eStart, 0, 0) #extracts radiation phase
srwl_uti_save_intens_ascii(arP2, mesh2, os.path.join(os.getcwd(), strDataFolderName, strPhPropOutFileName2), 0, ['', 'Horizontal Position', 'Vertical Position', 'Phase'], _arUnits=['', 'm', 'm', 'rad'])
del wfr2
print('done')

print('   Propagating Wavefront (through 2 B4C Mirrors) ... ', end='')
srwl.PropagElecField(wfr3, optBL_B4C)
print('done')
print('   Saving resulting Intensity and Phase data to files ... ', end='')
mesh3 = deepcopy(wfr3.mesh)
arI3 = array('f', [0]*mesh3.nx*mesh3.ny)
srwl.CalcIntFromElecField(arI3, wfr3, 6, 0, 3, mesh3.eStart, 0, 0)
srwl_uti_save_intens_ascii(arI3, mesh3, os.path.join(os.getcwd(), strDataFolderName, strIntPropOutFileName3), 0)

arP3 = array('d', [0]*mesh3.nx*mesh3.ny)
srwl.CalcIntFromElecField(arP3, wfr3, 0, 4, 3, mesh3.eStart, 0, 0)
srwl_uti_save_intens_ascii(arP3, mesh3, os.path.join(os.getcwd(), strDataFolderName, strPhPropOutFileName3), 0, ['', 'Horizontal Position', 'Vertical Position', 'Phase'], _arUnits=['', 'm', 'm', 'rad'])
del wfr3
print('done')

#**********************Plotting results (requires 3rd party graphics package)
print('   Plotting the results (blocks script execution; close any graph windows to proceed) ... ', end='')
plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
uti_plot2d1d(arI0, plotMesh0x, plotMesh0y, labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])


plotMesh1x = [1e+06*mesh1.xStart, 1e+06*mesh1.xFin, mesh1.nx]
plotMesh1y = [1e+06*mesh1.yStart, 1e+06*mesh1.yFin, mesh1.ny]
uti_plot2d1d(arI1, plotMesh1x, plotMesh1y, labels=['Horizontal Position [microns]', 'Vertical Position [microns]', 'Intensity After Propagation through 2 Elliptical Mirrors'])

plotMesh2x = [1e+06*mesh2.xStart, 1e+06*mesh2.xFin, mesh2.nx]
plotMesh2y = [1e+06*mesh2.yStart, 1e+06*mesh2.yFin, mesh2.ny]
uti_plot2d1d(arI2, plotMesh2x, plotMesh2y, labels=['Horizontal Position [microns]', 'Vertical Position [microns]', 'Intensity After Propagation through 2 Platinum Elliptical Mirrors'])

plotMesh3x = [1e+06*mesh3.xStart, 1e+06*mesh3.xFin, mesh3.nx]
plotMesh3y = [1e+06*mesh3.yStart, 1e+06*mesh3.yFin, mesh3.ny]
uti_plot2d1d(arI3, plotMesh3x, plotMesh3y, labels=['Horizontal Position [microns]', 'Vertical Position [microns]', 'Intensity After Propagation through 2 Boron Carbide Elliptical Mirrors'])



# Calculate relative intensity after propagation through the material mirrors
print()
print("maximum intensity for perfect reflectivity: ", max(arI1))
print("maximum intensity for Pt reflectivity: ", max(arI2))
print("maximum intensity for B4C reflectivity: ", max(arI3))
print("simulated relative intensity for Pt mirrors: ", max(arI2)/max(arI1))
print("simulated relative intensity for B4C mirrors: ", max(arI3)/max(arI1))

uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
print('done')
