#############################################################################
# SRWLIB Example # 22: Material-dependent reflective and refractive optics
# Simulating a Gaussian X-ray beam passing through an Aluminum filter,
# a B4C mirror, a 20 nm Pt-coated Si mirror, and an array of diamond CRLs
# v 0.01
#############################################################################

from __future__ import print_function

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
from copy import deepcopy


print('SRWLIB Python Example # 22:')
print('Simulating material-dependent reflective and refractive X-ray optics')


#******************* Photon Energy and Material Parameters
photon_energy = 8000. # [eV]

n_ph_en = 101
n_ang = 401
n_comp = 2
ph_en_start = 7500.
ph_en_fin = 8500.
ang_start = 0.001
ang_fin = 0.005
ang_graz = 0.003

# Aluminum filter
filter_thick = 50.e-06 # [m]
al_delta, al_atten_len = calc_delta_atten_len('Al', photon_energy)
opFilter = srwl_opt_setup_transm_from_material(
    'Al', filter_thick, photon_energy,
    _rx=2.e-03, _ry=2.e-03, _nx=2, _ny=2, _ext_tr=1
)

# B4C mirror; B4C is a formula not registered with a default xraydb density
b4c_density = 2.5 # [g/cm^3]
b4c_refl = calc_refl(
    'B4C', _dens=b4c_density,
    _n_ph_en=n_ph_en, _n_ang=n_ang, _n_comp=n_comp,
    _ph_en_start=ph_en_start, _ph_en_fin=ph_en_fin,
    _ang_start=ang_start, _ang_fin=ang_fin,
)

# 20 nm Pt coating on a Si substrate; xraydb expects thickness in Angstroms
pt_coating_thick = 200. # [Angstrom]
pt_si_refl = calc_refl_coated(
    'Pt', pt_coating_thick, 'Si',
    _n_ph_en=n_ph_en, _n_ang=n_ang, _n_comp=n_comp,
    _ph_en_start=ph_en_start, _ph_en_fin=ph_en_fin,
    _ang_start=ang_start, _ang_fin=ang_fin,
)

# Diamond CRL material: elemental C at diamond density
diamond_density = 3.52 # [g/cm^3]
diamond_delta, diamond_atten_len = calc_delta_atten_len(
    'C', photon_energy, _dens=diamond_density
)

print('Al filter: delta = {:.6g}, attenuation length = {:.6g} m'.format(
    al_delta, al_atten_len
))
print('Diamond: delta = {:.6g}, attenuation length = {:.6g} m'.format(
    diamond_delta, diamond_atten_len
))


#********************** Output Files
strDataFolderName = 'data_example_22'
strIntInFileName = 'ex22_res_int_in.dat'
strIntOutFileName = 'ex22_res_int_out.dat'

if not os.path.exists(strDataFolderName):
    os.makedirs(strDataFolderName)


#********************** Gaussian Beam Source
gsnBm = SRWLGsnBm()
gsnBm.x = 0
gsnBm.y = 0
gsnBm.z = 0
gsnBm.xp = 0
gsnBm.yp = 0
gsnBm.avgPhotEn = photon_energy
gsnBm.pulseEn = 0.001
gsnBm.repRate = 1
gsnBm.polar = 1 # linear horizontal polarization
gsnBm.sigX = 20.e-06
gsnBm.sigY = 20.e-06
gsnBm.sigT = 10.e-15
gsnBm.mx = 0
gsnBm.my = 0

wfr = SRWLWfr()
wfr.allocate(1, 151, 151)
wfr.mesh.zStart = 10.
wfr.mesh.eStart = photon_energy
wfr.mesh.eFin = photon_energy
wfr.mesh.xStart = -0.3e-03
wfr.mesh.xFin = 0.3e-03
wfr.mesh.yStart = -0.3e-03
wfr.mesh.yFin = 0.3e-03
wfr.unitElFld = 2
wfr.partBeam.partStatMom1.x = gsnBm.x
wfr.partBeam.partStatMom1.y = gsnBm.y
wfr.partBeam.partStatMom1.z = gsnBm.z
wfr.partBeam.partStatMom1.xp = gsnBm.xp
wfr.partBeam.partStatMom1.yp = gsnBm.yp

srwl.CalcElecFieldGaussian(wfr, gsnBm, [2])

meshIn = deepcopy(wfr.mesh)
arIIn = array('f', [0]*meshIn.nx*meshIn.ny)
srwl.CalcIntFromElecField(arIIn, wfr, 6, 0, 3, meshIn.eStart, 0, 0)
srwl_uti_save_intens_ascii(
    arIIn, meshIn, os.path.join(strDataFolderName, strIntInFileName), 0
)


#********************** Optical Elements
# Uniform Aluminum filter.
filter_amp = opFilter.arTr[0]

mirror_len = 0.2
mirror_width = 5.e-03

# B4C plane mirror, deflecting vertically
opMirB4C = SRWLOptMirPl(
    _size_tang=mirror_len, _size_sag=mirror_width,
    _nvx=0, _nvy=cos(ang_graz), _nvz=-sin(ang_graz),
    _tvx=0, _tvy=-sin(ang_graz),
    _refl=b4c_refl,
    _n_ph_en=n_ph_en, _n_ang=n_ang, _n_comp=n_comp,
    _ph_en_start=ph_en_start, _ph_en_fin=ph_en_fin,
    _ang_start=ang_start, _ang_fin=ang_fin,
)

# Pt-coated Si plane mirror, deflecting horizontally
opMirPtSi = SRWLOptMirPl(
    _size_tang=mirror_len, _size_sag=mirror_width,
    _nvx=cos(ang_graz), _nvy=0, _nvz=-sin(ang_graz),
    _tvx=-sin(ang_graz), _tvy=0,
    _refl=pt_si_refl,
    _n_ph_en=n_ph_en, _n_ang=n_ang, _n_comp=n_comp,
    _ph_en_start=ph_en_start, _ph_en_fin=ph_en_fin,
    _ang_start=ang_start, _ang_fin=ang_fin,
)

# Array of ten 2D parabolic diamond CRLs
crl_apert = 1.e-03
crl_r_min = 100.e-06
crl_number = 10
crl_wall_thick = 20.e-06
if diamond_delta > 0:
    opCRL = srwl_opt_setup_CRL(
        3, diamond_delta, diamond_atten_len, 1,
        crl_apert, crl_apert, crl_r_min, crl_number, crl_wall_thick,
        0, 0, None, 0, 0, 201, 201
    )
    crl_focal_len = opCRL.Fx
else:
    print('Warning: diamond refractive properties are unavailable; the CRL is disabled.')
    opCRL = SRWLOptT(
        _nx=2, _ny=2, _rx=crl_apert, _ry=crl_apert,
        _extTr=1, _alloc_base=[1, 0]
    )
    crl_focal_len = 1.

opDrFilter_M1 = SRWLOptD(1.)
opDrM1_M2 = SRWLOptD(1.)
opDrM2_CRL = SRWLOptD(1.)
opDrCRL_Obs = SRWLOptD(crl_focal_len)


#********************** Propagation
ppOpt = [0, 0, 1., 1, 0, 1., 1., 1., 1., 0, 0, 0]
ppDrift = [0, 0, 1., 1, 0, 1., 1., 1., 1., 0, 0, 0]
ppCRL = [0, 0, 1., 1, 0, 1.1, 2., 1.1, 2., 0, 0, 0]
ppFinal = [0, 0, 1., 1, 0, 0.2, 2., 0.2, 2., 0, 0, 0]

optBL = SRWLOptC(
    [opFilter, opDrFilter_M1, opMirB4C, opDrM1_M2,
     opMirPtSi, opDrM2_CRL, opCRL, opDrCRL_Obs],
    [ppOpt, ppDrift, ppOpt, ppDrift, ppOpt, ppDrift, ppCRL, ppFinal]
)

print('   Propagating wavefront through filter, mirrors, and diamond CRLs ... ', end='')
srwl.PropagElecField(wfr, optBL)
print('done')

meshOut = deepcopy(wfr.mesh)
arIOut = array('f', [0]*meshOut.nx*meshOut.ny)
srwl.CalcIntFromElecField(arIOut, wfr, 6, 0, 3, meshOut.eStart, 0, 0)
srwl_uti_save_intens_ascii(
    arIOut, meshOut, os.path.join(strDataFolderName, strIntOutFileName), 0
)


def integrated_intensity(arI, mesh):
    dx = (mesh.xFin - mesh.xStart)/(mesh.nx - 1) if mesh.nx > 1 else 1.
    dy = (mesh.yFin - mesh.yStart)/(mesh.ny - 1) if mesh.ny > 1 else 1.
    return sum(arI)*dx*dy


intIn = integrated_intensity(arIIn, meshIn)
intOut = integrated_intensity(arIOut, meshOut)
print('CRL focal length: {:.6g} m'.format(crl_focal_len))
print('Al filter intensity transmission: {:.6g}'.format(filter_amp*filter_amp))
print('Integrated beamline transmission: {:.6g}'.format(intOut/intIn))


#********************** Plotting
uti_plot2d1d(
    arIIn,
    [1.e+03*meshIn.xStart, 1.e+03*meshIn.xFin, meshIn.nx],
    [1.e+03*meshIn.yStart, 1.e+03*meshIn.yFin, meshIn.ny],
    labels=['Horizontal Position [mm]', 'Vertical Position [mm]', 'Input Intensity']
)
uti_plot2d1d(
    arIOut,
    [1.e+06*meshOut.xStart, 1.e+06*meshOut.xFin, meshOut.nx],
    [1.e+06*meshOut.yStart, 1.e+06*meshOut.yFin, meshOut.ny],
    labels=['Horizontal Position [microns]', 'Vertical Position [microns]', 'Focused Intensity']
)
uti_plot_show()
print('done')
