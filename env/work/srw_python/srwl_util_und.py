#############################################################################
# SRWLib for Python: Undulator Utilities v 0.0
#############################################################################

from srwlib import *

#****************************************************************************
def srwl_und_max_flux_spec(_e_beam, _mesh, _per, _len, _bx_vs_par, _by_vs_par, [_par1_start, _par1_step], [_par1_start, _par1_step], _par1_to_use, _par2_step_to_use, _min_pol_rate):
    """
    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR
    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _per: undulator period
    :param _len: undulator length
    :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
    :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
    :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
    :param _file_path: path to file for saving intermediate average Stokes data by master process
    :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
    :param _opt_bl: optical beamline (container) to propagate the radiation through (SRWLOptC type)
    :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1)
    :param _char: radiation characteristic to calculate: 0- Intensity (s0); 1- Four Stokes components; 2- Mutual Intensity Cut vs X; 3- Mutual Intensity Cut vs Y; 4- Mutual Intensity Cut vs X & Y
    :param _x0: horizontal center position for mutual intensity calculation
    :param _y0: vertical center position for mutual intensity calculation
    """

