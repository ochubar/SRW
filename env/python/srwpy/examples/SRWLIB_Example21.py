# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # 21: Simulating wavefront propagation through an optical layout
# containing a Hartmann sensor dedicated for analyzing aberrations of optical elements. 
# The example makes use of the SRW Virtual Beamline (srwl_bl.py) module, and
# is developed based on script generated authomatically by the Sirepo web
# interface to SRW (https://www.sirepo.com/srw#/).
# Based on developments by L. Huang (BNL / NSLS-II)
# v 0.02
#############################################################################

try:
    __IPYTHON__
    import sys
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
import math
import time

#*********************************Story
help_str='''SRWLIB Python Example # 21:
Simulating wavefront propagation through an optical layout containing a Hartmann sensor dedicated for analyzing aberrations of optical elements. 
The simulation is set up for the basic parameters of the ISR beamline of NSLS-II.
This example demonstrates the use of the function 'srwlib.srwl_opt_setup_Hartmann_sensor()' to simulate the generation of a Hartmanngram.

By default, this script runs in the single-electron mode, generating a fully coherent undulator radiation beam passing through the sensor,
with the Hartmanngram observed at some distance from the sensor.
To make it runing in the multi-electron mode for partially coherent undulator radiation beam, one can change the following line:
['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],
to
['wm', '', '1', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],
or add the line
v.wm = True
in the main(), before 'srwl_bl.SRWLBeamline(..).calc_all(..)'.
Note that the calculation in the multi-electron mode may be very time-consuming. 
The script can run in parallel mode (via MPI). For this, it should be started via mpiexec, e.g.:
mpiexec -n 11 python SRWLIB_Example21.py
or via slurm on an HPC cluster.

The resulting intensity distribution data are saved to the folder 'data_example_21'.
The single-electron data is saved to 'ex21_res_int_pr_se.dat' file; the multi-electron data to 'ex21_res_int_pr_me.dat' file.

The upsampling factor 'upsampling_factor' (before the varParam list) is used for controlling the number of data points in the final intensity distribution.
For high-accuracy wavefront reconstruction from this intensity distribution of the Hartmanngram, increasing of the upsampling factor may be required.
'''

#*********************************Beamline
def set_optics(v):
    el = []
    pp = []
    names = ['Ap_HFM', 'HFM', 'HFM_Before_VFM', 'Before_VFM', 'Ap_VFM', 'VFM', 'VFM_Before_SSA', 'Before_SSA', 'SSA', 'After_SSA',
             'After_SSA_Before_Mask', 'Before_Mask', 'HM', 'HM_HD']
    for el_name in names:
        if el_name == 'Ap_HFM':
            # Ap_HFM: aperture 31.3m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Ap_HFM_shape,
                _ap_or_ob='a',
                _Dx=v.op_Ap_HFM_Dx,
                _Dy=v.op_Ap_HFM_Dy,
                _x=v.op_Ap_HFM_x,
                _y=v.op_Ap_HFM_y,
            ))
            pp.append(v.op_Ap_HFM_pp)
        elif el_name == 'HFM':
            # HFM: sphericalMirror 31.3m
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
        elif el_name == 'HFM_Before_VFM':
            # HFM_Before_VFM: drift 31.3m
            el.append(srwlib.SRWLOptD(
                _L=v.op_HFM_Before_VFM_L,
            ))
            pp.append(v.op_HFM_Before_VFM_pp)
        elif el_name == 'Before_VFM':
            # Before_VFM: watch 37.9m
            pass
        elif el_name == 'Ap_VFM':
            # Ap_VFM: aperture 37.9m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_Ap_VFM_shape,
                _ap_or_ob='a',
                _Dx=v.op_Ap_VFM_Dx,
                _Dy=v.op_Ap_VFM_Dy,
                _x=v.op_Ap_VFM_x,
                _y=v.op_Ap_VFM_y,
            ))
            pp.append(v.op_Ap_VFM_pp)
        elif el_name == 'VFM':
            # VFM: sphericalMirror 37.9m
            el.append(srwlib.SRWLOptMirSph(
                _r=v.op_VFM_r,
                _size_tang=v.op_VFM_size_tang,
                _size_sag=v.op_VFM_size_sag,
                _nvx=v.op_VFM_nvx,
                _nvy=v.op_VFM_nvy,
                _nvz=v.op_VFM_nvz,
                _tvx=v.op_VFM_tvx,
                _tvy=v.op_VFM_tvy,
                _x=v.op_VFM_x,
                _y=v.op_VFM_y,
            ))
            pp.append(v.op_VFM_pp)
        elif el_name == 'VFM_Before_SSA':
            # VFM_Before_SSA: drift 37.9m
            el.append(srwlib.SRWLOptD(
                _L=v.op_VFM_Before_SSA_L,
            ))
            pp.append(v.op_VFM_Before_SSA_pp)
        elif el_name == 'Before_SSA':
            # Before_SSA: watch 49.2m
            pass
        elif el_name == 'SSA':
            # SSA: aperture 49.2m
            el.append(srwlib.SRWLOptA(
                _shape=v.op_SSA_shape,
                _ap_or_ob='a',
                _Dx=v.op_SSA_Dx,
                _Dy=v.op_SSA_Dy,
                _x=v.op_SSA_x,
                _y=v.op_SSA_y,
            ))
            pp.append(v.op_SSA_pp)
        elif el_name == 'After_SSA':
            # After_SSA: watch 49.2m
            pass
        elif el_name == 'After_SSA_Before_Mask':
            # After_SSA_Before_Mask: drift 49.2m
            el.append(srwlib.SRWLOptD(
                _L=v.op_After_SSA_Before_Mask_L,
            ))
            pp.append(v.op_After_SSA_Before_Mask_pp)
        elif el_name == 'Before_Mask':
            # Before_Mask: watch 64.2m
            pass

        # Hartmann sensor
        elif el_name == 'HM':
            # HM:
            print('\nSetting up Hartmann mask ... ', end='')
            t0 = time.time()
            opTr = srwlib.srwl_opt_setup_Hartmann_sensor(
                _delta=1.0, _atten_len=1.0,
                _thick=v.op_HM_thickness,
                _grid_sh=v.op_HM_grid_shape,
                _grid_dx=v.op_HM_grid_dx,
                _grid_dy=v.op_HM_grid_dy,
                _pitch_x=v.op_HM_grid_pitch_x,
                _pitch_y=v.op_HM_grid_pitch_y,
                _grid_nx=v.op_HM_grid_nx,
                _grid_ny=v.op_HM_grid_nx,
                _grid_angle=v.op_HM_grid_angle,
                _mask_nx=v.op_HM_mask_nx,
                _mask_ny=v.op_HM_mask_ny,
                _hx=v.op_HM_mask_hx,
                _hy=v.op_HM_mask_hy,
                _mask_x0=v.op_HM_mask_x,
                _mask_y0=v.op_HM_mask_y,
            )
            el.append(opTr)

            print('completed (lasted', round(time.time() - t0, 6), 's)')
            pp.append(v.op_HM_pp)
        elif el_name == 'HM_HD':
            # HM_HD:
            el.append(srwlib.SRWLOptD(
                _L=v.op_HM_HD_L,
            ))
            pp.append(v.op_HM_HD_pp)

    pp.append(v.op_fin_pp)
    return srwlib.SRWLOptC(el, pp)

#*********************************Parameters
upsampling_factor = 2 #Upsampling factor used for controlling the number of data points in the final intensity distribution.

varParam = [
    ['name', 's', 'Hartmanngram simulation under NSLS-II 4-ID ISR beamline layout', 'simulation name'],

#---Data Folder
    ['fdir', 's', 'data_example_21', 'folder (directory) name for reading-in input and saving output data files'],

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
    ['ebm_z', 'f', 0., 'electron beam initial average longitudinal position [m]'],
    ['ebm_dr', 'f', -1.5, 'electron beam longitudinal drift [m] to be performed before a required calculation'],
    ['ebm_ens', 'f', 0.0008, 'electron beam relative energy spread'],
    ['ebm_emx', 'f', 7.6e-10, 'electron beam horizontal emittance [m]'],
    ['ebm_emy', 'f', 7e-12, 'electron beam vertical emittance [m]'],
    # Definition of the beam through Twiss:
    ['ebm_betax', 'f', 20.1, 'horizontal beta-function [m]'],
    ['ebm_betay', 'f', 3.4, 'vertical beta-function [m]'],
    ['ebm_alphax', 'f', 0.0, 'horizontal alpha-function [rad]'],
    ['ebm_alphay', 'f', 0.0, 'vertical alpha-function [rad]'],
    ['ebm_etax', 'f', 0.0, 'horizontal dispersion function [m]'],
    ['ebm_etay', 'f', 0.0, 'vertical dispersion function [m]'],
    ['ebm_etaxp', 'f', 0.0, 'horizontal dispersion function derivative [rad]'],
    ['ebm_etayp', 'f', 0.0, 'vertical dispersion function derivative [rad]'],

#---Undulator
    ['und_bx', 'f', 0.0, 'undulator horizontal peak magnetic field [T]'],
    ['und_by', 'f', 0.7512, 'undulator vertical peak magnetic field [T]'],
    ['und_phx', 'f', 0.0, 'initial phase of the horizontal magnetic field [rad]'],
    ['und_phy', 'f', 0.0, 'initial phase of the vertical magnetic field [rad]'],
    ['und_b2e', '', '', 'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by', 'store_true'],
    ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e', 'store_true'],
    ['und_per', 'f', 0.023, 'undulator period [m]'],
    ['und_len', 'f', 2.8, 'undulator length [m]'],
    ['und_zc', 'f', 0.0, 'undulator center longitudinal position [m]'],
    ['und_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
    ['und_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
    ['und_g', 'f', 0.0, 'undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
    ['und_ph', 'f', 0.0, 'shift of magnet arrays [mm] for which the field should be set up'],
    ['und_mdir', 's', '', 'name of magnetic measurements sub-folder'],
    ['und_mfs', 's', '', 'name of magnetic measurements for different gaps summary file'],

#---Calculation Types
    # Electron Trajectory
    ['tr', '', '', 'calculate electron trajectory', 'store_true'],
    ['tr_cti', 'f', 0.0, 'initial time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_ctf', 'f', 0.0, 'final time moment (c*t) for electron trajectory calculation [m]'],
    ['tr_np', 'f', 10000, 'number of points for trajectory calculation'],
    ['tr_mag', 'i', 1, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
    ['tr_fn', 's', 'ex21_res_trj.dat', 'file name for saving calculated trajectory data'],
    ['tr_pl', 's', '', 'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

    #Single-Electron Spectrum vs Photon Energy
    ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
    ['ss_ei', 'f', 100.0, 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ef', 'f', 20000.0, 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
    ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
    ['ss_x', 'f', 0.0, 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_y', 'f', 0.0, 'vertical position [m] for single-e spectrum vs photon energy calculation'],
    ['ss_meth', 'i', 1, 'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['ss_prec', 'f', 0.01, 'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
    ['ss_pol', 'i', 6, 'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['ss_mag', 'i', 1, 'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
    ['ss_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
    ['ss_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['ss_fn', 's', 'ex21_res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
    ['ss_pl', 's', '', 'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    #Multi-Electron Spectrum vs Photon Energy (taking into account e-beam emittance, energy spread and collection aperture size)
    ['sm', '', '', 'calculate multi-e spectrum vs photon energy', 'store_true'],
    ['sm_ei', 'f', 100.0, 'initial photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ef', 'f', 20000.0, 'final photon energy [eV] for multi-e spectrum vs photon energy calculation'],
    ['sm_ne', 'i', 10000, 'number of points vs photon energy for multi-e spectrum vs photon energy calculation'],
    ['sm_x', 'f', 0.0, 'horizontal center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_rx', 'f', 0.001, 'range of horizontal position / horizontal aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_nx', 'i', 1, 'number of points vs horizontal position for multi-e spectrum vs photon energy calculation'],
    ['sm_y', 'f', 0.0, 'vertical center position [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ry', 'f', 0.001, 'range of vertical position / vertical aperture size [m] for multi-e spectrum vs photon energy calculation'],
    ['sm_ny', 'i', 1, 'number of points vs vertical position for multi-e spectrum vs photon energy calculation'],
    ['sm_mag', 'i', 1, 'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
    ['sm_hi', 'i', 1, 'initial UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_hf', 'i', 15, 'final UR spectral harmonic to be taken into account for multi-e spectrum vs photon energy calculation'],
    ['sm_prl', 'f', 1.0, 'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_pra', 'f', 1.0, 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
    ['sm_meth', 'i', -1, 'method to use for spectrum vs photon energy calculation in case of arbitrary input magnetic field: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler", -1- dont use this accurate integration method (rather use approximate if possible)'],
    ['sm_prec', 'f', 0.01, 'relative precision for spectrum vs photon energy calculation in case of arbitrary input magnetic field (nominal value is 0.01)'],
    ['sm_nm', 'i', 1, 'number of macro-electrons for calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_na', 'i', 5, 'number of macro-electrons to average on each node at parallel (MPI-based) calculation of spectrum in case of arbitrary input magnetic field'],
    ['sm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons) for intermediate intensity at calculation of multi-electron spectrum in case of arbitrary input magnetic field'],
    ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
    ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['sm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['sm_fn', 's', 'ex21_res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
    ['sm_pl', 's', '', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
    #to add options for the multi-e calculation from "accurate" magnetic field

    #Power Density Distribution vs horizontal and vertical position
    ['pw', '', '', 'calculate SR power density distribution', 'store_true'],
    ['pw_x', 'f', 0.0, 'central horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_rx', 'f', 0.015, 'range of horizontal position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_nx', 'i', 100, 'number of points vs horizontal position for calculation of power density distribution'],
    ['pw_y', 'f', 0.0, 'central vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ry', 'f', 0.015, 'range of vertical position [m] for calculation of power density distribution vs horizontal and vertical position'],
    ['pw_ny', 'i', 100, 'number of points vs vertical position for calculation of power density distribution'],
    ['pw_pr', 'f', 1.0, 'precision factor for calculation of power density distribution'],
    ['pw_meth', 'i', 1, 'power density computation method (1- "near field", 2- "far field")'],
    ['pw_zst', 'f', 0., 'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_zfi', 'f', 0., 'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
    ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
    ['pw_fn', 's', 'ex21_res_pow.dat', 'file name for saving calculated power density distribution'],
    ['pw_pl', 's', '', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #Single-Electron Intensity distribution vs horizontal and vertical position
    ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],
    #Single-Electron Wavefront Propagation
    ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],
    #Multi-Electron (partially-coherent) Wavefront Propagation
    ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    ['w_e', 'f', 11300.0, 'photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ef', 'f', -1.0, 'final photon energy [eV] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ne', 'i', 1, 'number of points vs photon energy for calculation of intensity distribution'],
    ['w_x', 'f', 0.0, 'central horizontal position [m] for calculation of intensity distribution'],
    ['w_rx', 'f', 0.0021, 'range of horizontal position [m] for calculation of intensity distribution'],
    ['w_nx', 'i', 100, 'number of points vs horizontal position for calculation of intensity distribution'],
    ['w_y', 'f', 0.0, 'central vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ry', 'f', 0.0015, 'range of vertical position [m] for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_ny', 'i', 100, 'number of points vs vertical position for calculation of intensity distribution'],
    ['w_smpf', 'f', 0.2, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
    ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
    ['w_u', 'i', 1, 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
    ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
    ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
    ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],

    ['si_fn', 's', 'ex21_res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
    ['si_pl', 's', '', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
    ['ws_fni', 's', 'ex21_res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
    ['ws_pl', 's', '', 'plot the resulting intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    ['wm_nm', 'i', 100000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
    ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node for parallel (MPI-based) calculation of multi-electron wavefront propagation'],
    ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
    ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y; 40- intensity(s0), mutual intensity cuts and degree of coherence vs X & Y'],
    ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
    ['wm_x0', 'f', 0.0, 'horizontal center position for mutual intensity cut calculation'],
    ['wm_y0', 'f', 0.0, 'vertical center position for mutual intensity cut calculation'],
    ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
    ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
    ['wm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],
    ['wm_fni', 's', 'ex21_res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],
    ['wm_ff', 's', 'ascii', 'format of file name for saving propagated multi-e intensity distribution vs horizontal and vertical position (ascii and hdf5 supported)'],

    # Optics parameters
    ['op_r', 'f', 31.3, 'longitudinal position of the first optical element [m]'],
    # Former appParam:
    ['rs_type', 's', 'u', 'source type, (u) idealized undulator, (t), tabulated undulator, (m) multipole, (g) gaussian beam'],

#---Beamline optics:
    # Ap_HFM: aperture
    ['op_Ap_HFM_shape', 's', 'r', 'shape'],
    ['op_Ap_HFM_Dx', 'f', 0.00208, 'horizontalSize'],
    ['op_Ap_HFM_Dy', 'f', 0.01, 'verticalSize'],
    ['op_Ap_HFM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Ap_HFM_y', 'f', 0.0, 'verticalOffset'],

    # HFM: sphericalMirror
    ['op_HFM_hfn', 's', 'None', 'heightProfileFile'],
    ['op_HFM_dim', 's', 'x', 'orientation'],
    ['op_HFM_r', 'f', 8759.69355847405, 'radius'],
    ['op_HFM_size_tang', 'f', 1.0, 'tangentialSize'],
    ['op_HFM_size_sag', 'f', 0.1, 'sagittalSize'],
    ['op_HFM_ang', 'f', 0.0026, 'grazingAngle'],
    ['op_HFM_nvx', 'f', 0.999996620001904, 'normalVectorX'],
    ['op_HFM_nvy', 'f', 0.0, 'normalVectorY'],
    ['op_HFM_nvz', 'f', -0.0025999970706676568, 'normalVectorZ'],
    ['op_HFM_tvx', 'f', 0.0025999970706676568, 'tangentialVectorX'],
    ['op_HFM_tvy', 'f', 0.0, 'tangentialVectorY'],
    ['op_HFM_amp_coef', 'f', 0., 'heightAmplification'],
    ['op_HFM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_HFM_y', 'f', 0.0, 'verticalOffset'],

    # HFM_Before_VFM: drift
    ['op_HFM_Before_VFM_L', 'f', 6.6, 'length'],

    # Ap_VFM: aperture
    ['op_Ap_VFM_shape', 's', 'r', 'shape'],
    ['op_Ap_VFM_Dx', 'f', 0.01, 'horizontalSize'],
    ['op_Ap_VFM_Dy', 'f', 0.00156, 'verticalSize'],
    ['op_Ap_VFM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_Ap_VFM_y', 'f', 0.0, 'verticalOffset'],

    # VFM: sphericalMirror
    ['op_VFM_hfn', 's', 'None', 'heightProfileFile'],
    ['op_VFM_dim', 's', 'y', 'orientation'],
    ['op_VFM_r', 'f', 6695.9036898, 'radius'],
    ['op_VFM_size_tang', 'f', 1.0, 'tangentialSize'],
    ['op_VFM_size_sag', 'f', 0.1, 'sagittalSize'],
    ['op_VFM_ang', 'f', 0.0026, 'grazingAngle'],
    ['op_VFM_nvx', 'f', 0.0, 'normalVectorX'],
    ['op_VFM_nvy', 'f', 0.999996620001904, 'normalVectorY'],
    ['op_VFM_nvz', 'f', -0.0025999970706676568, 'normalVectorZ'],
    ['op_VFM_tvx', 'f', 0.0, 'tangentialVectorX'],
    ['op_VFM_tvy', 'f', 0.0025999970706676568, 'tangentialVectorY'],
    ['op_VFM_amp_coef', 'f', 0., 'heightAmplification'],
    ['op_VFM_x', 'f', 0.0, 'horizontalOffset'],
    ['op_VFM_y', 'f', 0.0, 'verticalOffset'],

    # VFM_Before_SSA: drift
    ['op_VFM_Before_SSA_L', 'f', 11.3, 'length'],

    # SSA: aperture
    ['op_SSA_shape', 's', 'r', 'shape'],
    ['op_SSA_Dx', 'f', 1e-03, 'horizontalSize'],
    ['op_SSA_Dy', 'f', 1e-03, 'verticalSize'],
    ['op_SSA_x', 'f', 0.0, 'horizontalOffset'],
    ['op_SSA_y', 'f', 0.0, 'verticalOffset'],

    # After_SSA_Before_Mask: drift
    ['op_After_SSA_Before_Mask_L', 'f', 15.0, 'length'],

    # HM (Hartmann sensor)
    ['op_HM_thickness', 'f', 1.e-03, 'Thickness of the mask plate'],
    ['op_HM_grid_shape', 'i', 1, 'Shape of grid: 0: circle hartmann; 1: rectangle hartmann; 2: 2D grating'],
    ['op_HM_grid_dx', 'f', 5.e-06, 'Size of the grid hole in x # The design is 6 um, but the measurement is 5 um'],
    ['op_HM_grid_dy', 'f', 5.e-06, 'Size of the grid hole in y # The design is 6 um, but the measurement is 5 um'],
    ['op_HM_grid_pitch_x', 'f', 20.e-06, 'Grid pitch in x'],
    ['op_HM_grid_pitch_y', 'f', 20.e-06, 'Grid pitch in y'],
    ['op_HM_grid_nx', 'i', 151, 'Number of the grid holes in x # ideally 151 for the HXR sensor'],
    ['op_HM_grid_ny', 'i', 151, 'Number of the grid holes in y # ideally 151 for the HXR sensor'],
    ['op_HM_grid_angle', 'f', (90.0-25.0)*math.pi/180., 'Grid holes rotation angle [rad]'],
    ['op_HM_mask_nx', 'i', 2001, 'Number of of Transmission data points in x'],
    ['op_HM_mask_ny', 'i', 2001, 'Number of of Transmission data points in y'],
    ['op_HM_mask_hx', 'f', 0.5e-06, 'Transmission data step size in x'],
    ['op_HM_mask_hy', 'f', 0.5e-06, 'Transmission data step size in y'],
    ['op_HM_mask_x', 'f', 0., 'Center of the mask in x'],
    ['op_HM_mask_y', 'f', 0., 'Center of the mask in y'],

    # HM_HD: drift
    ['op_HM_HD_L', 'f', 0.2, 'length'],

    # Detector (determining mesh of the final intensity distribution)
    ['d_x', 'f', 0., 'central horizontal position [m] of detector active area'],
    ['d_rx', 'f', 1.48e-06*512, 'range of horizontal position [m] of detector active area (should be >0 to be taken into account)'],
    ['d_nx', 'i', 512*upsampling_factor, 'number of pixels vs horizontal position (should be >0 to be taken into account)'],
    ['d_dx', 'f', 0., 'horizontal size of pixel [m]'],
    ['d_y', 'f', 0., 'central vertical position [m] of detector active area'],
    ['d_ry', 'f', 1.48e-06*512, 'range of vertical position [m] of detector active area (should be >0 to be taken into account)'],
    ['d_ny', 'i', 512*upsampling_factor, 'number of pixels vs vertical position (should be >0 to be taken into account)'],
    ['d_dy', 'f', 0., 'vertical size of pixel [m]'],
    ['d_or', 'f', 1, 'interpolation order (i.e. order of polynomials to be used at 2D interpolation)'],
    ['d_ifn', 's', '', 'file name with detector spectral efficiency data (on 1D mesh vs photon energy: _eStart, _eFin, _ne)'],

#---Propagation parameters
    ['op_Ap_HFM_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.5, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Ap_HFM'],
    ['op_HFM_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HFM'],
    ['op_HFM_Before_VFM_pp', 'f',        [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HFM_Before_VFM'],
    ['op_Ap_VFM_pp', 'f',                [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'Ap_VFM'],
    ['op_VFM_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VFM'],
    ['op_VFM_Before_SSA_pp', 'f',        [0, 0, 1.0, 4, 0, 2.0, 1.3, 1.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'VFM_Before_SSA'],
    ['op_SSA_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'SSA'],
    ['op_After_SSA_Before_Mask_pp', 'f', [0, 0, 1.0, 3, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'After_SSA_Before_Mask'],
    ['op_HM_pp', 'f',                    [0, 0, 1.0, 0, 0, 1.0, 5.0, 1.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HM'],
    ['op_HM_HD_pp', 'f',                 [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'HM_HD'],
    ['op_fin_pp', 'f',                   [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 'final post-propagation (resize) parameters'],

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

#*********************************Run All
def main():

    if srwlib.srwl_uti_proc_is_master(): print(help_str)

    v = srwl_bl.srwl_uti_parse_options(srwl_bl.srwl_uti_ext_options(varParam), use_sys_argv=True)
    op = set_optics(v)

    if srwlib.srwl_uti_proc_is_master(): 
        v.ws = True
        v.ws_pl = 'xy'
    #v.wm = True #Uncomment this to run multi-electron partially coherent calculation
   
    srwl_bl.SRWLBeamline(_name=v.name).calc_all(v, op)

main()
