# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example: Virtual Beamline: a set of utilities and functions allowing to simulate
# operation of an SR Beamline.
# The standard use of this script is from command line, with some optional arguments,
# e.g. for calculation (with default parameter values) of:
# UR Spectrum Through a Slit (Flux within a default aperture):
#    python SRWLIB_ExampleVirtBL_01.py --sm
# Single-Electron UR Spectrum (Flux per Unit Surface):
#    python SRWLIB_ExampleVirtBL_01.py --ss
# UR Power Density (at the first optical element):
#    python SRWLIB_ExampleVirtBL_01.py --pw
# Input Single-Electron UR Intensity Distribution (at the first optical element):
#    python SRWLIB_ExampleVirtBL_01.py --si
# Single-Electron Wavefront Propagation:
#    python SRWLIB_ExampleVirtBL_01.py --ws
# Multi-Electron Wavefront Propagation:
#  Sequential Mode:
#    python SRWLIB_ExampleVirtBL_01.py --wm
#  Parallel Mode (using MPI / mpi4py), e.g.:
#    mpiexec -n 6 python SRWLIB_ExampleVirtBL_01.py --wm
# For changing parameters of all these calculaitons from the default valuse, see the definition
# of all options in the list at the end of the script.
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwl_bl import *
#import time

#*********************************Setting Up Optical Elements and Propagation Parameters
def set_optics(_v):
    """This function describes optical layout of the Coherent Hoard X-ray (CHX) beamline of NSLS-II.
    Such function has to be written for every beamline to be simulated; it is specific to a particular beamline.
    :param _v: structure containing all parameters allowed to be varied for that particular beamline
    """
    
#---Nominal Positions of Optical Elements [m] (with respect to straight section center)
    zS0 = 20.5 #S0 (primary slit)
    zHDM = 27.4 #Horizontally-Deflecting Mirror (HDM)
    zS1 = 29.9 #S1 slit
    zS2 = 34.3 #S2 slit
    zBPM = 34.6 #BPM for beam visualization
    zCRL = 35.4 #+tzCRL*1e-3 #CRL transfocator (corrected by translation)
    zKL = 44.5 #+tzKL*1e-3 #Kinoform Lens for horizontal focusing (corrected by translation)
    zS3 = 48. #S3 slit ('pinhole', waist position)
    zSample = 48.7 #Sample position, COR of diffractometer

#---Instantiation of the Optical Elements
    arElNamesAll = ['S0', 'S0_HDM', 'HDM_S1', 'S1', 'S1_S2', 'S2', 'S2_BPM', 'BPM_CRL', 'CRL1', 'CRL2', 'CRL_KL', 'KLA', 'KL', 'KL_S3', 'S3', 'S3_SMP']

    #Treat beamline sub-cases / alternative configurations
    if(len(v.op_HDM_ifn) > 0): arElNamesAll.insert(2, 'HDM')

    if(len(v.op_fin) > 0):
        if(v.op_fin not in arElNamesAll): raise Exception('Optical element with the name specified in the "op_fin" option is not present in this beamline')
        #Could be made more general

    arElNames = []; 
    for i in range(len(arElNamesAll)):
        arElNames.append(arElNamesAll[i])
        if(len(v.op_fin) > 0):
            if(arElNamesAll[i] == v.op_fin): break

    el = []; pp = [] #lists of SRW optical element objects and their corresponding propagation parameters

    #S0 (primary slit)
    if('S0' in arElNames): 
        el.append(SRWLOptA('r', 'a', v.op_S0_dx, v.op_S0_dy)); pp.append(v.op_S0_pp)

    #Drift S0 -> HDM
    if('S0_HDM' in arElNames): 
        el.append(SRWLOptD(zHDM - zS0)); pp.append(v.op_S0_HDM_pp)
        
    #HDM (Height Profile Error)
    if('HDM' in arElNames): 
        horApHDM = 0.94e-03 #Projected dimensions
        verApHDM = 1.e-03
        angHDM = 3.1415926e-03 #? grazing angle
        ifnHDM = os.path.join(v.fdir, v.op_HDM_ifn) if len(v.op_HDM_ifn) > 0 else ''
        if(len(ifnHDM) > 0):
            hProfDataHDM = srwl_uti_read_data_cols(ifnHDM, '\t', 0, 1)
            opHDM = srwl_opt_setup_surf_height_1d(hProfDataHDM, 'x', _ang=angHDM, _amp_coef=v.op_HDM_amp, _nx=1000, _ny=200, _size_x=horApHDM, _size_y=verApHDM)
            ofnHDM = os.path.join(v.fdir, v.op_HDM_ofn) if len(v.op_HDM_ofn) > 0 else ''
            if(len(ofnHDM) > 0):
                pathDifHDM = opHDM.get_data(3, 3)
                srwl_uti_save_intens_ascii(pathDifHDM, opHDM.mesh, ofnHDM, 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Dif.'], _arUnits=['', 'm', 'm', 'm'])
            el.append(opHDM); pp.append(v.op_HDM_pp)
            
    #Drift HDM -> S1
    if('HDM_S1' in arElNames): 
        el.append(SRWLOptD(zS1 - zHDM + v.op_S1_dz)); pp.append(v.op_HDM_S1_pp)

    #S1 slit
    if('S1' in arElNames): 
        el.append(SRWLOptA('r', 'a', v.op_S1_dx, v.op_S1_dy)); pp.append(v.op_S1_pp)

    #Drift S1 -> S2
    if('S1_S2' in arElNames): 
        el.append(SRWLOptD(zS2 - zS1 + v.op_S2_dz)); pp.append(v.op_S1_S2_pp)

    #S2 slit
    if('S2' in arElNames): 
        el.append(SRWLOptA('r', 'a', v.op_S2_dx, v.op_S2_dy)); pp.append(v.op_S2_pp)

    #Drift S2 -> BPM
    if('S2_BPM' in arElNames): 
        el.append(SRWLOptD(zBPM - zS2 + v.op_BPM_dz)); pp.append(v.op_S2_BPM_pp)

    #Drift BPM -> CRL
    if('BPM_CRL' in arElNames): 
        el.append(SRWLOptD(zCRL - zBPM + v.op_CRL_dz)); pp.append(v.op_BPM_CRL_pp)

    #CRL1 (1D, vertically-focusing)
    if('CRL1' in arElNames):
        el.append(srwl_opt_setup_CRL(2, v.op_CRL1_delta, v.op_CRL1_atnl, 1, v.op_CRL1_apnf, v.op_CRL1_apf, v.op_CRL1_rmin, v.op_CRL1_n, v.op_CRL1_thck, 0, 0))
        pp.append(v.op_CRL1_pp)

    #CRL2 (1D, vertically-focusing)
    if('CRL2' in arElNames):
        el.append(srwl_opt_setup_CRL(2, v.op_CRL2_delta, v.op_CRL2_atnl, 1, v.op_CRL2_apnf, v.op_CRL2_apf, v.op_CRL2_rmin, v.op_CRL2_n, v.op_CRL2_thck, 0, 0))
        pp.append(v.op_CRL2_pp)

    #Drift CRL -> KL
    if('CRL_KL' in arElNames):
        el.append(SRWLOptD(zKL - zCRL + v.op_KL_dz)); pp.append(v.op_CRL_KL_pp)

    #KL Aperture
    if('KLA' in arElNames): 
        el.append(SRWLOptA('r', 'a', v.op_KLA_dx, v.op_KLA_dy)); pp.append(v.op_KLA_pp)

    #KL (1D, horizontally-focusing)
    if('KL' in arElNames):
        el.append(SRWLOptL(v.op_KL_fx, v.op_KL_fy)) #KL as Ideal Lens; to make it a transmission element with a profile read from a file
        pp.append(v.op_KL_pp)

    #Drift KL -> S3
    if('KL_S3' in arElNames):
        el.append(SRWLOptD(zS3 - zKL + v.op_S3_dz)); pp.append(v.op_KL_S3_pp)

    #S3 slit
    if('S3' in arElNames):
        el.append(SRWLOptA('r', 'a', v.op_S3_dx, v.op_S3_dy)); pp.append(v.op_S3_pp)

    #Drift S3 -> Sample
    if('S3_SMP' in arElNames):
        el.append(SRWLOptD(zSample - zS3 + v.op_SMP_dz)); pp.append(v.op_S3_SMP_pp)

    pp.append(v.op_fin_pp)

    return SRWLOptC(el, pp)

#*********************************List of Parameters allowed to be varied
#---List of supported options / commands / parameters allowed to be varied for this Beamline (comment-out unnecessary):   
varParam = [
#---Data Folder
    ['fdir', 's', os.path.join(os.getcwd(), 'data_example_virt_bl_01'), 'folder (directory) name for reading-in input and saving output data files'],

#---Electron Beam
    ['ebm_nm', 's', 'NSLS-II Low Beta ', 'standard electron beam name'],
    ['ebm_nms', 's', 'Day1', 'standard electron beam name suffix: e.g. can be Day1, Final'],
    ['ebm_i', 'f', 0.5, 'electron beam current [A]'],
    #['ebeam_e', 'f', 3., 'electron beam avarage energy [GeV]'],
    ['ebm_de', 'f', 0., 'electron beam average energy deviation [GeV]'],
    ['ebm_x', 'f', 0., 'electron beam initial average horizontal position [m]'],
    ['ebm_y', 'f', 0., 'electron beam initial average vertical position [m]'],
    ['ebm_z', 'f', 0., 'electron beam initial average longitudinal position [m]'],
    ['ebm_dr', 'f', -1.65, 'electron beam longitudinal drift [m] to be performed before a required calculation'],
        
#---Undulator
    #['und_per', 'f', 0.02, 'undulator period [m]'],
    #['und_len', 'f', 3., 'undulator length [m]'],
    ['und_b', 'f', 0.88770981, 'undulator vertical peak magnetic field [T]'],
    #['und_bx', 'f', 0., 'undulator horizontal peak magnetic field [T]'],
    #['und_by', 'f', 1., 'undulator vertical peak magnetic field [T]'],
    #['und_phx', 'f', 1.5708, 'undulator horizontal magnetic field phase [rad]'],
    #['und_phy', 'f', 0., 'undulator vertical magnetic field phase [rad]'],
    #['und_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
    #['und_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
    #['und_zc', 'f', 0., 'undulator center longitudinal position [m]'],

    ['und_b2e', '', '', 'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by', 'store_true'],
    ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e', 'store_true'],

#---Calculation Types
    #Single-Electron Spectrum vs Photon Energy
    ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
    ['ss_ei', 'f', 100., 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ef', 'f', 20000., 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
    ['ss_x', 'f', 0., 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_y', 'f', 0., 'vertical position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_meth', 'i', 1, 'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['ss_prec', 'f', 0.01, 'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
    ['ss_pol', 'i', 6, 'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['ss_mag', 'i', 1, 'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
    ['ss_fn', 's', 'res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
    ['ss_pl', 's', 'e', 'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    #Multi-Electron Spectrum vs Photon Energy (taking into account e-beam emittance, energy spread and collection aperture size)
    ['sm', '', '', 'calculate multi-e spectrum vs photon energy', 'store_true'],
    ['sm_ei', 'f', 100., 'initial photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ef', 'f', 20000., 'final photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ne', 'i', 10000, 'number of points vs photon energy for multi-e spectrum vs photon energy calculation'],
    ['sm_x', 'f', 0., 'horizontal center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_rx', 'f', 0.001, 'range of horizontal position / horizontal aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_nx', 'i', 1, 'number of points vs horizontal position for multi-e spectrum vs photon energy calculation'],
    ['sm_y', 'f', 0., 'vertical center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ry', 'f', 0.001, 'range of vertical position / vertical aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ny', 'i', 1, 'number of points vs vertical position for multi-e spectrum vs photon energy calculation'],
    ['sm_mag', 'i', 1, 'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
    ['sm_hi', 'i', 1, 'initial UR spectral harmonic to be taken into accountfor multi-e spectrum vs photon energy calculation'],
    ['sm_hf', 'i', 15, 'final UR spectral harmonic to be taken into accountfor multi-e spectrum vs photon energy calculation'],
    ['sm_prl', 'f', 1., 'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_pra', 'f', 1., 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
    ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
    ['sm_pl', 's', 'e', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
    #to add options for the multi-e calculation from "accurate" magnetic field

    #Power Density Distribution vs horizontal and vertical position
    ['pw', '', '', 'calculate SR power density distribution', 'store_true'],
    ['pw_x', 'f', 0., 'central horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_rx', 'f', 0.015, 'range of horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_nx', 'i', 100, 'number of points vs horizontal position for calculation of power density distribution'],
    ['pw_y', 'f', 0., 'central vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ry', 'f', 0.015, 'range of vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ny', 'i', 100, 'number of points vs vertical position for calculation of power density distribution'],
    ['pw_pr', 'f', 1., 'precision factor for calculation of power density distribution'],
    ['pw_meth', 'i', 1, 'power density computation method (1- "near field", 2- "far field")'],
    ['pw_zst', 'f', 0., 'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_zfi', 'f', 0., 'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
    ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
    ['pw_pl', 's', 'xy', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],
    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],
    
    ['w_e', 'f', 9000., 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1., 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0., 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.4e-03, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0., 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.6e-03, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', 1., 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],

    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', 'xy', 'plot the propagated radiaiton intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['si_pl', 's', 'xy', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 100000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y'],
    ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_fni', 's', 'res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],

    #['ws_fn', 's', '', 'file name for saving single-e (/ fully coherent) wavefront data'],
    #['wm_fn', 's', '', 'file name for saving multi-e (/ purtially coherent) wavefront data'],
    #to add options
    
    ['op_r', 'f', 20.5, 'longitudinal position of the first optical element [m]'],
    ['op_fin', 's', 'S3_SMP', 'name of the final optical element wavefront has to be propagated through'],
    
    #NOTE: the above option/variable names (fdir, ebm*, und*, ss*, sm*, pw*, is*, ws*, wm*) should be the same in all beamline scripts
    #on the other hand, the beamline optics related options below (op*) are specific to a particular beamline (and can be differ from beamline to beamline).
    #However, the default values of all the options/variables (above and below) can differ from beamline to beamline.
    
#---Beamline Optics
    ['op_S0_dx', 'f', 0.2e-03, 'slit S0: horizontal size [m]'],
    ['op_S0_dy', 'f', 1.0e-03, 'slit S0: vertical size [m]'],
    ['op_HDM_ifn', 's', 'CHX_HDM_height_prof_1d.dat', 'mirror HDM: input file name of height profile data'],
    ['op_HDM_amp', 'f', 1., 'mirror HDM: amplification coefficient for height profile data'],
    ['op_HDM_ofn', 's', 'res_CHX_HDM_opt_path_dif.dat', 'mirror HDM: output file name of optical path difference data'],
    ['op_S1_dz', 'f', 0., 'S1: offset of longitudinal position [m]'],
    ['op_S1_dx', 'f', 0.2e-03, 'slit S1: horizontal size [m]'],
    ['op_S1_dy', 'f', 1.0e-03, 'slit S1: vertical size [m]'],
    ['op_S2_dz', 'f', 0., 'S2: offset of longitudinal position [m]'],
    ['op_S2_dx', 'f', 0.05e-03, 'slit S2: horizontal size [m]'],
    ['op_S2_dy', 'f', 1.0e-03, 'slit S2: vertical size [m]'],
    ['op_BPM_dz', 'f', 0., 'BPM: offset of longitudinal position [m]'],
    ['op_CRL_dz', 'f', 0., 'CRL: offset of longitudinal position [m]'],
    ['op_CRL1_delta', 'f', 4.20756805e-06, 'CRL1: refractive index decrements of material'],
    ['op_CRL1_atnl', 'f', 7312.94e-06, 'CRL1: attenuation length of material [m]'],
    ['op_CRL1_apnf', 'f', 1.e-03, 'CRL1: geometrical aparture of 1D CRL in the plane where there is no focusing'],
    ['op_CRL1_apf', 'f', 2.4e-03, 'CRL1: geometrical aparture of 1D CRL in the focusing plane'],
    ['op_CRL1_rmin', 'f', 1.5e-03, 'CRL1: radius of curface curvature at the tip of parabola [m]'],
    ['op_CRL1_n', 'i', 1, 'CRL1: number of individual lenses'],
    ['op_CRL1_thck', 'f', 80.e-06, 'CRL1: wall thickness (at the tip of parabola) [m]'],
    ['op_CRL2_delta', 'f', 4.20756805e-06, 'CRL2: refractive index decrements of material'],
    ['op_CRL2_atnl', 'f', 7312.94e-06, 'CRL2: attenuation length of material [m]'],
    ['op_CRL2_apnf', 'f', 1.e-03, 'CRL2: geometrical aparture of 1D CRL in the plane where there is no focusing'],
    ['op_CRL2_apf', 'f', 1.4e-03, 'CRL2: geometrical aparture of 1D CRL in the focusing plane'],
    ['op_CRL2_rmin', 'f', 0.5e-03, 'CRL2: radius of curface curvature at the tip of parabola [m]'],
    ['op_CRL2_n', 'i', 6, 'CRL2: number of individual lenses'],
    ['op_CRL2_thck', 'f', 80.e-06, 'CRL2: wall thickness (at the tip of parabola) [m]'],
    ['op_KLA_dx', 'f', 1.4e-03, 'KL Aperture: horizontal size [m]'],
    ['op_KLA_dy', 'f', 0.2e-03, 'KL Aperture: vertical size [m]'],
    ['op_KL_dz', 'f', 0., 'KL: offset of longitudinal position [m]'],
    ['op_KL_fx', 'f', 3.24479, 'KL: horizontal focal length [m]'],
    ['op_KL_fy', 'f', 1.e+23, 'KL: vertical focal length [m]'],
    ['op_S3_dz', 'f', 0., 'S3: offset of longitudinal position [m]'],
    ['op_S3_dx', 'f', 10.e-06, 'slit S3: horizontal size [m]'],
    ['op_S3_dy', 'f', 10.e-06, 'slit S3: vertical size [m]'],
    ['op_SMP_dz', 'f', 0., 'Sample: offset of longitudinal position [m]'],

    #to add options for different beamline cases, etc.

    #Propagation Param.:   [0][1][2][3][4] [5]  [6]  [7]  [8] [9][10][11]
    #['op_S0_pp', 'f',      [0, 0, 1, 0, 0, 4.5, 5.0, 1.5, 2.5, 0, 0, 0], 'slit S0: propagation parameters'],
    ['op_S0_pp', 'f',      [0, 0, 1, 0, 0, 2.5, 5.0, 1.5, 2.5, 0, 0, 0], 'slit S0: propagation parameters'],
    ['op_S0_HDM_pp', 'f',  [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift S0 -> HDM: propagation parameters'],
    ['op_HDM_pp', 'f',     [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'mirror HDM: propagation parameters'],
    ['op_HDM_S1_pp', 'f',  [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift HDM -> S1: propagation parameters'],
    ['op_S1_pp', 'f',      [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'slit S1: propagation parameters'],
    ['op_S1_S2_pp', 'f',   [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift S1 -> S2: propagation parameters'],
    ['op_S2_pp', 'f',      [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'slit S2: propagation parameters'],
    ['op_S2_BPM_pp', 'f',  [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift S2 -> BPM: propagation parameters'],
    ['op_BPM_CRL_pp', 'f', [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift BPM -> CRL: propagation parameters'],
    ['op_CRL1_pp', 'f',    [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'CRL1: propagation parameters'],
    ['op_CRL2_pp', 'f',    [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'CRL2: propagation parameters'],
    ['op_CRL_KL_pp', 'f',  [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift CRL -> KL: propagation parameters'],
    ['op_KLA_pp', 'f',     [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'KL Aperture: propagation parameters'],
    #['op_KL_pp', 'f',      [0, 0, 1, 0, 0, 1.0, 5.0, 1.0, 7.0, 0, 0, 0], 'KL: propagation parameters'],
    ['op_KL_pp', 'f',      [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'KL: propagation parameters'],
    ['op_KL_S3_pp', 'f',   [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift KL -> S3: propagation parameters'],
    #['op_S3_pp', 'f',      [0, 0, 1, 0, 0, 0.3, 3.0, 0.3, 3.0, 0, 0, 0], 'slit S3: propagation parameters'],
    ['op_S3_pp', 'f',      [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'slit S3: propagation parameters'],
    #['op_S3_SMP_pp', 'f',  [0, 0, 1, 1, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift S3 -> Sample: propagation parameters'],
    ['op_S3_SMP_pp', 'f',  [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'drift S3 -> Sample: propagation parameters'],
    #['op_fin_pp', 'f',     [0, 0, 1, 0, 1, 0.1, 5.0, 1.0, 1.5, 0, 0, 0], 'final post-propagation (resize) parameters'],
    ['op_fin_pp', 'f',     [0, 0, 1, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0], 'final post-propagation (resize) parameters'],

    #[ 0]: Auto-Resize (1) or not (0) Before propagation
    #[ 1]: Auto-Resize (1) or not (0) After propagation
    #[ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    #[ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
    #[ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    #[ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
    #[ 6]: Horizontal Resolution modification factor at Resizing
    #[ 7]: Vertical Range modification factor at Resizing
    #[ 8]: Vertical Resolution modification factor at Resizing
    #[ 9]: Type of wavefront Shift before Resizing (not yet implemented)
    #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
    #[12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
    #[13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
    #[14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
    #[15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
    #[16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
]
    
#*********************************Entry
if __name__ == "__main__":

#---Parse options, defining Beamline elements and running calculations
    v = srwl_uti_parse_options(varParam)
    
#---Add some constant "parameters" (not allowed to be varied) for the beamline
    v.und_per = 0.02 #['und_per', 'f', 0.02, 'undulator period [m]'],
    v.und_len = 3. #['und_len', 'f', 3., 'undulator length [m]'],
    v.und_zc = 0. #['und_zc', 'f', 0., 'undulator center longitudinal position [m]'],
    v.und_sy = -1 #['und_sy', 'i', -1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
    v.und_sx = 1 #['und_sx', 'i', 1, 'undulator vertical magnetic field symmetry vs longitudinal position'],

#---Setup optics only if Wavefront Propagation is required:
    op = set_optics(v) if(v.ws or v.wm) else None

#---Run all requested calculations
    SRWLBeamline('Coherent Hard X-ray beamline').calc_all(v, op)
