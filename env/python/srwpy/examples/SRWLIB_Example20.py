# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # 20: Calculating 4D cross-spectral density (/mutual intensity) of partially-coherent
# undulator radiation at a distance from the undulator, performing its coherent mode decomposition, 
# and illustrating propagation of the coherent modes through a high-resolution microscopy beamline.
# The example makes use of the SRW Virtual Beamline (srwl_bl.py) module, and is developed based on
# script generated authomatically by the Sirepo web interface to SRW (https://www.sirepo.com/srw#/)
# Based on developments by R. Li
# v 0.04
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import sys
try:
    __IPYTHON__
    del sys.argv[1:]
except:
    pass

try: #OC15112022
    import sys
    sys.path.append('../')
    import srwl_bl
    import srwlib
except:
    from srwpy import srwl_bl
    from srwpy import srwlib
#import srwl_bl
#import srwlib

#*********************************Story
help_str='''SRWLIB Python Example # 20:
Calculating 4D cross-spectral density (CSD) of partially-coherent undulator radiation, performing its coherent mode decomposition
(CMD) and simulating the propagation of the coherent modes through a high-resolution microscopy beamline. The electron beam,
undulator, and optical layout parameters are defined in the list below (see comments next to the corresponding parameters).
!!!!!Under testing!!!!!

To calculate (on a multi-core server computer) the 4D cross-spectral density at a distance from undulator, using MPI-parallelization, and perform
the coherent mode decomposition right after this, execute something like:

    mpiexec -n 40 python SRWLIB_Example20.py --wm --wm_nop --wm_ch=61 --wm_nm=300000 --wm_nmm=10 --wm_ncm=800 --wm_na=100 --wm_ns=100 --wm_ff=h5 --wm_fni=ex20_res

The exact meaning of the options used explicitly (as well as those used implicitly) for this calculation, can be found in the varParam list below. 
Depending on available resources (CPU and memory), the number of processes ("-n *") can be larger or smaller, as well as the number of MPI 
"masters"/groups to be used ("--wm_nmm=*"). The total number of processes and the number of MPI groups are related. We found it close to optimum
to have ~4 (or 3) MPI processes in each "MPI group", i.e. if one uses "-n 40" (or "-n 30") then one would preferably use "--wm_nmm=10".

Instead of calculating the CSD and making its CMD in one run, it can be more practical to calculate CSD first, save it to a file,
and then perform CMD after reading-in the CSD from that file, e.g. as follows:

    mpiexec -n 40 python SRWLIB_Example20.py --wm --wm_nop --wm_ch=6 --wm_nm=300000 --wm_nmm=10 --wm_na=100 --wm_ns=100 --wm_ff=h5 --wm_fni=ex20_res
    python SRWLIB_Example20.py --wm_ch=7 --wm_fnmi=ex20_res_mi.h5 --wm_ncm=800 --wm_fni=ex20_res
    
In the above example, the file core name is defined by "--wm_fni=ex20_res", which results in creation of the 4D CSD / mutual ntensity file named
"ex20_res_mi.h5" and the coherent modes file named "ex20_res_cm.h5" (i.e. suffixes "_mi" and "_cm" are added to the file name core, as well as
the file extension ".h5", corresponding to the specified HDF5 file type, "--wm_ff=h5").
Depending on optical conditions, the initial CSD / CMD calculations can be done for smaller transverse ranges and at smaller numbers
of points. For this one needs to change (/reduce values of) parameters w_rx, w_ry, w_nx, w_ny, either by using the mechanism of optins
or by changing the values in the list below. This may result in a faster calculation.

To simulate propagation of the coherent modes (after these were calculated by running this script as described above) through
the beamline defined in the function "set_optics", and to derive intensity and degree of coherence after the propagation, execute e.g.:

    mpiexec -n 41 python SRWLIB_Example20.py --wm --wm_ch=41 --wm_fncm=ex20_res_cm.h5 --op_S1_pp=[0,0,1.,0,0,5.0,5.0,1.5,6.0,0.,0.,0.,0.,0.,0.,0.,0.] --wm_fni=ex20_res_pr.dat

In the above example, propagation parameters are specified by the option "--op_S1_pp=*". Note that these parameters can be different for
coherent modes, as compared to "direct calculations" without the CMD, see below.

The radiation intensity and degree of coherence can be also calculated "directly" in this example, without using the CMD:

    mpiexec -n 41 python SRWLIB_Example20.py --wm --wm_ch=41 --w_smpf=0.08 --wm_nm=300000 --op_S1_pp=[0,0,1.,0,0,1.5,15.0,1.5,6.0,0.,0.,0.,0.,0.,0.,0.,0.] --wm_na=20 --wm_ns=20 --wm_fni=ex20_res_pr_dir.dat

Further on, the 4D CSD, followed by CMD, can also be performed at any location of a beamline; to make these calculations for the sample position,
try running this:

    mpiexec -n 40 python SRWLIB_Example20.py --wm --wm_ch=61 --w_smpf=0.08 --wm_nm=300000 --wm_nmm=10 --wm_ns=100 --op_S1_pp=[0,0,1.,0,0,1.5,15.0,1.5,6.0,0.,0.,0.,0.,0.,0.,0.,0.] --op_fin_pp=[0,0,1.,0,0,0.06,0.5,0.17,1.,0.,0.,0.,0.,0.,0.,0.,0.] --d_rx=200.e-09 --d_nx=100 --d_ry=200.e-09 --d_ny=100 --wm_ncm=500 --wm_ff=h5 --wm_fni=ex20_res_pr_dir

In the above call example, the transverse mesh of the final wavefront (and the CSD and Coherent Modes) is defined via "detector" params
(options "--d_rx=200.e-09 --d_nx=100 --d_ry=200.e-09 --d_ny=100").

For more information on parallel calculations under "mpi4py" please see documentation to the "mpi4py" and MPI.
Note that some long-lasting partially-coherent UR calculations save from time to time instant average intensity / degree of coherence 
or other characteristics to file(s), so the execution of the long loop over "macro-electrons" can be aborted after some time without 
the danger that all results will be lost.
'''

#*********************************Beamline
def set_optics(v=None):
    el = []
    pp = []
    names = ['S1', 'S1_HCM', 'HCM', 'HCM_HFM', 'HFM', 'HFM_CRL1', 'CRL1', 'CRL2', 'CRL2_SSA', 'SSA', 'SSA_AFFO', 'AFFO', 'FFO', 'FFO_Sample']
    for el_name in names:
        if el_name == 'S1':
            # S1: aperture 20.5m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_S1_shape,
                _ap_or_ob='a',
                _Dx=v.op_S1_Dx,
                _Dy=v.op_S1_Dy,
                _x=v.op_S1_x,
                _y=v.op_S1_y,
            ))
            pp.append(v.op_S1_pp)
        elif el_name == 'S1_HCM':
            # S1_HCM: drift 20.5m
            el.append(srwlib.SRWLOptD(
                _L=v.op_S1_HCM_L,
            ))
            pp.append(v.op_S1_HCM_pp)
        elif el_name == 'HCM':
            # HCM: sphericalMirror 28.35m
            el.append(srwlib.SRWLOptMirSph(
                _r=v.op_HCM_r,
                _size_tang=v.op_HCM_size_tang,
                _size_sag=v.op_HCM_size_sag,
                _nvx=v.op_HCM_nvx,
                _nvy=v.op_HCM_nvy,
                _nvz=v.op_HCM_nvz,
                _tvx=v.op_HCM_tvx,
                _tvy=v.op_HCM_tvy,
                _x=v.op_HCM_x,
                _y=v.op_HCM_y,
            ))
            pp.append(v.op_HCM_pp)
        elif el_name == 'HCM_HFM':
            # HCM_HFM: drift 28.35m
            el.append(srwlib.SRWLOptD(
                _L=v.op_HCM_HFM_L,
            ))
            pp.append(v.op_HCM_HFM_pp)
        elif el_name == 'HFM':
            # HFM: sphericalMirror 32.64m
            el.append(srwlib.SRWLOptMirSph(
                _r=v.op_HFM_r,
                _size_tang=v.op_HFM_size_tang,
                _size_sag=v.op_HFM_size_sag,
                _nvx=v.op_HFM_nvx,
                _nvy=v.op_HFM_nvy,
                _nvz=v.op_HFM_nvz,
                _tvx=v.op_HFM_tvx,
                _tvy=v.op_HFM_tvy,
                _x=v.op_HFM_x,
                _y=v.op_HFM_y,
            ))
            pp.append(v.op_HFM_pp)
        elif el_name == 'HFM_CRL1':
            # HFM_CRL1: drift 32.64m
            el.append(srwlib.SRWLOptD(
                _L=v.op_HFM_CRL1_L,
            ))
            pp.append(v.op_HFM_CRL1_pp)
        elif el_name == 'CRL1':
            # CRL1: crl 34.15m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL1_foc_plane,
                _delta=v.op_CRL1_delta,
                _atten_len=v.op_CRL1_atten_len,
                _shape=v.op_CRL1_shape,
                _apert_h=v.op_CRL1_apert_h,
                _apert_v=v.op_CRL1_apert_v,
                _r_min=v.op_CRL1_r_min,
                _n=v.op_CRL1_n,
                _wall_thick=v.op_CRL1_wall_thick,
                _xc=v.op_CRL1_x,
                _yc=v.op_CRL1_y,
            ))
            pp.append(v.op_CRL1_pp)
        elif el_name == 'CRL2':
            # CRL2: crl 34.15m
            el.append(srwlib.srwl_opt_setup_CRL(
                _foc_plane=v.op_CRL2_foc_plane,
                _delta=v.op_CRL2_delta,
                _atten_len=v.op_CRL2_atten_len,
                _shape=v.op_CRL2_shape,
                _apert_h=v.op_CRL2_apert_h,
                _apert_v=v.op_CRL2_apert_v,
                _r_min=v.op_CRL2_r_min,
                _n=v.op_CRL2_n,
                _wall_thick=v.op_CRL2_wall_thick,
                _xc=v.op_CRL2_x,
                _yc=v.op_CRL2_y,
            ))
            pp.append(v.op_CRL2_pp)
        elif el_name == 'CRL2_SSA':
            # CRL2_SSA: drift 34.15m
            el.append(srwlib.SRWLOptD(
                _L=v.op_CRL2_SSA_L,
            ))
            pp.append(v.op_CRL2_SSA_pp)
        elif el_name == 'SSA':
            # SSA: aperture 94.5m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_SSA_shape,
                _ap_or_ob='a',
                _Dx=v.op_SSA_Dx,
                _Dy=v.op_SSA_Dy,
                _x=v.op_SSA_x,
                _y=v.op_SSA_y,
            ))
            pp.append(v.op_SSA_pp)
        elif el_name == 'SSA_AFFO':
            # SSA_AFFO: drift 94.5m
            el.append(srwlib.SRWLOptD(
                _L=v.op_SSA_AFFO_L,
            ))
            pp.append(v.op_SSA_AFFO_pp)
        elif el_name == 'AFFO':
            # AFFO: aperture 109.0m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_AFFO_shape,
                _ap_or_ob='a',
                _Dx=v.op_AFFO_Dx,
                _Dy=v.op_AFFO_Dy,
                _x=v.op_AFFO_x,
                _y=v.op_AFFO_y,
            ))
            pp.append(v.op_AFFO_pp)
        elif el_name == 'FFO':
            # FFO: lens 109.0m
            el.append(srwlib.SRWLOptL(
                _Fx=v.op_FFO_Fx,
                _Fy=v.op_FFO_Fy,
                _x=v.op_FFO_x,
                _y=v.op_FFO_y,
            ))
            pp.append(v.op_FFO_pp)
        elif el_name == 'FFO_Sample':
            # FFO_Sample: drift 109.0m
            el.append(srwlib.SRWLOptD(
                _L=v.op_FFO_Sample_L,
            ))
            pp.append(v.op_FFO_Sample_pp)
    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)

#*********************************Parameters
varParam = srwl_bl.srwl_uti_ext_options([
    ['name', 's', 'NSLS-II HXN beamline for CMD', 'simulation name'],

#---Data Folder
    ['fdir', 's', 'data_example_20', 'folder (directory) name for reading-in input and saving output data files'],

#---Electron Beam
    ['ebm_nm', 's', '', 'standard electron beam name'],
    ['ebm_nms', 's', '', 'standard electron beam name suffix: e.g. can be Day1, Final'],
    ['ebm_i', 'f', 0.5, 'electron beam current [A]'],
    ['ebm_e', 'f', 3.0, 'electron beam avarage energy [GeV]'],
    ['ebm_de', 'f', 0.0, 'electron beam average energy deviation [GeV]'],
    ['ebm_x', 'f', 0.0, 'electron beam initial average horizontal position [m]'],
    ['ebm_y', 'f', 0.0, 'electron beam initial average vertical position [m]'],
    ['ebm_xp', 'f', 0.0, 'electron beam initial average horizontal angle [rad]'],
    ['ebm_yp', 'f', 0.0, 'electron beam initial average vertical angle [rad]'],
    ['ebm_z', 'f', 0.0, 'electron beam initial average longitudinal position [m]'],
    ['ebm_dr', 'f', -1.8, 'electron beam longitudinal drift [m] to be performed before a required calculation'],
    ['ebm_ens', 'f', 0.00089, 'electron beam relative energy spread'],
    ['ebm_emx', 'f', 7.6e-10, 'electron beam horizontal emittance [m]'],
    ['ebm_emy', 'f', 8e-12, 'electron beam vertical emittance [m]'],
    # Definition of the beam through Twiss params:
    ['ebm_betax', 'f', 1.84, 'horizontal beta-function [m]'],
    ['ebm_betay', 'f', 1.17, 'vertical beta-function [m]'],
    ['ebm_alphax', 'f', 0.0, 'horizontal alpha-function [rad]'],
    ['ebm_alphay', 'f', 0.0, 'vertical alpha-function [rad]'],
    ['ebm_etax', 'f', 0.0, 'horizontal dispersion function [m]'],
    ['ebm_etay', 'f', 0.0, 'vertical dispersion function [m]'],
    ['ebm_etaxp', 'f', 0.0, 'horizontal dispersion function derivative [rad]'],
    ['ebm_etayp', 'f', 0.0, 'vertical dispersion function derivative [rad]'],

#---Undulator
    ['und_bx', 'f', 0.0, 'undulator horizontal peak magnetic field [T]'],
    ['und_by', 'f', 0.88770981, 'undulator vertical peak magnetic field [T]'],
    
    ['und_phx', 'f', 0.0, 'initial phase of the horizontal magnetic field [rad]'],
    ['und_phy', 'f', 0.0, 'initial phase of the vertical magnetic field [rad]'],
    ['und_b2e', '', '', 'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by', 'store_true'],
    ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e', 'store_true'],
    ['und_per', 'f', 0.02, 'undulator period [m]'],
    ['und_len', 'f', 3.0, 'undulator length [m]'],
    ['und_zc', 'f', 0.0, 'undulator center longitudinal position [m]'],
    ['und_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
    ['und_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
    ['und_g', 'f', 0.0, 'undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
    ['und_ph', 'f', 0.0, 'shift of magnet arrays [mm] for which the field should be set up'],
    ['und_mdir', 's', '', 'name of magnetic measurements sub-folder'],
    ['und_mfs', 's', '', 'name of magnetic measurements for different gaps summary file'],

#---Calculation Types
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '1', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', 9000.0, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1.0, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],

    ['w_rx', 'f', 0.0025, 'range of horizontal position [m] for calculation of intensity distribution'],
    #['w_rx', 'f', 0.003, 'range of horizontal position [m] for calculation of intensity distribution'],

    ['w_nx', 'i', 624, 'number of points vs horizontal position for calculation of intensity distribution'],
    #['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],

    ['w_ry', 'f', 0.001, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    #['w_ry', 'f', 0.0007, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],

    ['w_ny', 'i', 60, 'number of points vs vertical position for calculation of intensity distribution'],
    #['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],

    ['w_smpf', 'f', -1, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    #['w_smpf', 'f', 0.05, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    #['w_smpf', 'f', 0.15, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],

    ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 1000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    #['wm_nm', 'i', 280000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 10, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 10, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    #['wm_na', 'i', 100, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    #['wm_ns', 'i', 100, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],

    #['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    ['wm_ch', 'i', 61, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    #['wm_ch', 'i', 6, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    #['wm_ch', 'i', 41, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],

    ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0.0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0.0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],

    ['wm_nmm', 'i', 14, 'number of MPI masters to use'],
    #['wm_nmm', 'i', 2, 'number of MPI masters to use'],
    ['wm_ncm', 'i', 1000, 'number of Coherent Modes to calculate'],

    ['wm_nop', '', '', 'switch forcing to do calculations ignoring any optics defined (by set_optics function)', 'store_true'],

    ['wm_fni', 's', 'ex20_res.h5', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],
    #['wm_fni', 's', 'ex20_res_pr.h5', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],

    ['wm_fnmi', 's', '', 'file name of input cross-spectral density / mutual intensity; if this file name is supplied, the initial cross-spectral density (for such operations as coherent mode decomposition) will not be calculated, but rathre it will be taken from that file.'],
    ['wm_fncm', 's', '', 'file name of input cross-spectral density / mutual intensity; if this file name is supplied, the initial cross-spectral density (for such operations as coherent mode decomposition) will not be calculated, but rather it will be taken from that file.'],
    ['wm_ff', 's', 'ascii', 'format of data file for saving propagated multi-e intensity distribution vs horizontal and vertical position (ascii and hdf5 supported)'],

    #['wm_fbk', '', '1', 'create backup file(s) with propagated multi-e intensity distribution vs horizontal and vertical position and other radiation characteristics', 'store_true'],

    #to add options
    ['op_r', 'f', 20.5, 'longitudinal position of the first optical element [m]'],
    # Former appParam:
    ['rs_type', 's', 'u', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # S1: aperture
    ['op_S1_shape', 's', 'r', 'shape'],
    ['op_S1_Dx', 'f', 0.0025, 'horizontalSize'],
    ['op_S1_Dy', 'f', 0.0007, 'verticalSize'],
    ['op_S1_x', 'f', 0.0, 'horizontalOffset'],
    ['op_S1_y', 'f', 0.0, 'verticalOffset'],

    # S1_HCM: drift
    ['op_S1_HCM_L', 'f', 7.85, 'length'],

    # HCM: sphericalMirror
    ['op_HCM_hfn', 's', 'None', 'heightProfileFile'],
    ['op_HCM_dim', 's', 'x', 'orientation'],
    ['op_HCM_r', 'f', 17718.8, 'radius'],
    ['op_HCM_size_tang', 'f', 1.0, 'tangentialSize'],
    ['op_HCM_size_sag', 'f', 0.006, 'sagittalSize'],
    ['op_HCM_ang', 'f', 0.0032, 'grazingAngle'],
    ['op_HCM_nvx', 'f', 0.9999948800043691, 'normalVectorX'],
    ['op_HCM_nvy', 'f', 0.0, 'normalVectorY'],
    ['op_HCM_nvz', 'f', -0.003199994538669463, 'normalVectorZ'],
    ['op_HCM_tvx', 'f', 0.003199994538669463, 'tangentialVectorX'],
    ['op_HCM_tvy', 'f', 0.0, 'tangentialVectorY'],
    ['op_HCM_amp_coef', 'f', 0.001, 'heightAmplification'],
    ['op_HCM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_HCM_y', 'f', 0.0, 'verticalOffset'],

    # HCM_HFM: drift
    ['op_HCM_HFM_L', 'f', 4.29, 'length'],

    # HFM: sphericalMirror
    ['op_HFM_hfn', 's', 'None', 'heightProfileFile'],
    ['op_HFM_dim', 's', 'x', 'orientation'],
    ['op_HFM_r', 'f', 38660.0, 'radius'],
    ['op_HFM_size_tang', 'f', 1.0, 'tangentialSize'],
    ['op_HFM_size_sag', 'f', 0.06, 'sagittalSize'],
    ['op_HFM_ang', 'f', 0.0032, 'grazingAngle'],
    ['op_HFM_nvx', 'f', -0.9999948800043691, 'normalVectorX'],
    ['op_HFM_nvy', 'f', 0.0, 'normalVectorY'],
    ['op_HFM_nvz', 'f', -0.003199994538669463, 'normalVectorZ'],
    ['op_HFM_tvx', 'f', -0.003199994538669463, 'tangentialVectorX'],
    ['op_HFM_tvy', 'f', 0.0, 'tangentialVectorY'],
    ['op_HFM_amp_coef', 'f', 0.001, 'heightAmplification'],
    ['op_HFM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_HFM_y', 'f', 0.0, 'verticalOffset'],

    # HFM_CRL1: drift
    ['op_HFM_CRL1_L', 'f', 1.51, 'length'],

    # CRL1: crl
    ['op_CRL1_foc_plane', 'f', 2, 'focalPlane'],
    ['op_CRL1_delta', 'f', 5.326453e-06, 'refractiveIndex'],
    ['op_CRL1_atten_len', 'f', 0.005276, 'attenuationLength'],
    ['op_CRL1_shape', 'f', 1, 'shape'],
    ['op_CRL1_apert_h', 'f', 0.003, 'horizontalApertureSize'],
    ['op_CRL1_apert_v', 'f', 0.0015, 'verticalApertureSize'],
    ['op_CRL1_r_min', 'f', 0.001, 'tipRadius'],
    ['op_CRL1_wall_thick', 'f', 5e-05, 'tipWallThickness'],
    ['op_CRL1_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL1_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL1_n', 'i', 3, 'numberOfLenses'],

    # CRL2: crl
    ['op_CRL2_foc_plane', 'f', 2, 'focalPlane'],
    ['op_CRL2_delta', 'f', 5.326453e-06, 'refractiveIndex'],
    ['op_CRL2_atten_len', 'f', 0.005276, 'attenuationLength'],
    ['op_CRL2_shape', 'f', 1, 'shape'],
    ['op_CRL2_apert_h', 'f', 0.003, 'horizontalApertureSize'],
    ['op_CRL2_apert_v', 'f', 0.0015, 'verticalApertureSize'],
    ['op_CRL2_r_min', 'f', 0.0015, 'tipRadius'],
    ['op_CRL2_wall_thick', 'f', 5e-05, 'tipWallThickness'],
    ['op_CRL2_x', 'f', 0.0, 'horizontalOffset'],
    ['op_CRL2_y', 'f', 0.0, 'verticalOffset'],
    ['op_CRL2_n', 'i', 2, 'numberOfLenses'],

    # CRL2_SSA: drift
    ['op_CRL2_SSA_L', 'f', 60.35, 'length'],

    # SSA: aperture
    ['op_SSA_shape', 's', 'r', 'shape'],
    ['op_SSA_Dx', 'f', 2e-05, 'horizontalSize'],
    ['op_SSA_Dy', 'f', 0.0001, 'verticalSize'],
    ['op_SSA_x', 'f', 0.0, 'horizontalOffset'],
    ['op_SSA_y', 'f', 0.0, 'verticalOffset'],

    # SSA_AFFO: drift
    ['op_SSA_AFFO_L', 'f', 14.5, 'length'],

    # AFFO: aperture
    ['op_AFFO_shape', 's', 'r', 'shape'],
    ['op_AFFO_Dx', 'f', 0.00015, 'horizontalSize'],
    ['op_AFFO_Dy', 'f', 0.00015, 'verticalSize'],
    ['op_AFFO_x', 'f', 0.0, 'horizontalOffset'],
    ['op_AFFO_y', 'f', 0.0, 'verticalOffset'],

    # FFO: lens
    ['op_FFO_Fx', 'f', 0.01814, 'horizontalFocalLength'],
    ['op_FFO_Fy', 'f', 0.01814, 'verticalFocalLength'],
    ['op_FFO_x', 'f', 0.0, 'horizontalOffset'],
    ['op_FFO_y', 'f', 0.0, 'verticalOffset'],

    # FFO_Sample: drift
    ['op_FFO_Sample_L', 'f', 0.018163, 'length'],

#---Propagation parameters
    ['op_S1_pp', 'f',         [0, 0, 1.0, 0, 0, 5.0, 5.0, 1.5, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'S1'],
    #['op_S1_pp', 'f',         [0, 0, 1.0, 0, 0, 1.2, 1.0, 1.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'S1'],

    ['op_S1_HCM_pp', 'f',     [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'S1_HCM'],
    ['op_HCM_pp', 'f',        [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HCM'],
    ['op_HCM_HFM_pp', 'f',    [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HCM_HFM'],
    ['op_HFM_pp', 'f',        [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HFM'],
    ['op_HFM_CRL1_pp', 'f',   [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HFM_CRL1'],
    ['op_CRL1_pp', 'f',       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL1'],
    ['op_CRL2_pp', 'f',       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL2'],
    ['op_CRL2_SSA_pp', 'f',   [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'CRL2_SSA'],
    ['op_SSA_pp', 'f',        [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'SSA'],
    ['op_SSA_AFFO_pp', 'f',   [0, 0, 1.0, 3, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'SSA_AFFO'],
    ['op_AFFO_pp', 'f',       [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'AFFO'],
    ['op_FFO_pp', 'f',        [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'FFO'],
    ['op_FFO_Sample_pp', 'f', [0, 0, 1.0, 4, 0, 1.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'FFO_Sample'],

    ['op_fin_pp', 'f',        [0, 0, 1.0, 0, 0, 0.07, 1.0, 0.25, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],
    #['op_fin_pp', 'f',        [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

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
])

#*********************************Run All
def main():

    if srwlib.srwl_uti_proc_is_master(): print(help_str)
    
    v = srwl_bl.srwl_uti_parse_options(varParam, use_sys_argv=True)

    op = set_optics(v)
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v, op)

main()
