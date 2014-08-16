# -*- coding: utf-8 -*-
#############################################################################
# SRWLib SR Beamline Base Class
# Contains a set of member objects and functions for simulating basic operation and characteristics
# of a complete user beamline in a synchrotron radiation source.
# Under development!!!
# v 0.02
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import optparse

#****************************************************************************
class SRWLBeamline(object):
    """Basic Synchrotron Radiation Beamline (in a storage ring based source)
    Several different types of simulations are planned to be supported, including:
    - single-electron UR intensity spectrum;
    - UR spectrum through a slit (finite aperture), taking into account e-beam emittance and energy spread;
    - single-electron (fully-coherent) and multi-electron (partially-coherent) wavefront propagation calculations;
    - simulation of (coherent scattering type) experiments for different types of samples and beamline settings;
    - misc. auxiliary calculations/estimations facilitating beamline operation (e.g. estimation of undulator gap and crystal monochromator angle for a required photon energy, etc.)
    """

    #------------------------------------------------------------------------
    def __init__(self, _name='', _e_beam=None, _mag_approx=None, _mag=None, _gsn_beam=None, _op=None):
        """
        :param _name: beamline name string
        :param _e_beam: electron beam (SRWLPartBeam instance)
        :param _mag_approx: approximate magnetic field container (SRWLMagFldC instance)
        :param _mag: accurate magnetic field container (SRWLMagFldC instance)
        :param _gsn_beam: coherent Gaussian beam (SRWLGsnBm instance)
        :param _op: optics sequence (SRWLOptC instance)
        """

        if _e_beam != None:
            if(isinstance(_e_beam, SRWLPartBeam) == False):
                raise Exception("Incorrect Electron Beam structure")
        if _mag_approx != None:
            if(isinstance(_mag_approx, SRWLMagFldC) == False):
                raise Exception("Incorrect Magnetic Field Container structure")
        if _mag != None:
            if(isinstance(_mag, SRWLMagFldC) == False):
                raise Exception("Incorrect Magnetic Field Container structure")
        if _gsn_beam != None:
            if(isinstance(_gsn_beam, SRWLGsnBm) == False):
                raise Exception("Incorrect Gaussian Beam structure") 
        if _op != None:
            if(isinstance(_op, SRWLOptC) == False):
                raise Exception("Incorrect Ptical Container structure") 

        self.name = _name
        self.eBeam = SRWLPartBeam() if _e_beam is None else _e_beam
        self.mag_approx = _mag_approx
        self.mag = _mag
        self.gsnBeam = _gsn_beam
        self.optics = _op

    #------------------------------------------------------------------------
    def set_e_beam(self, _e_beam_name='', _e_beam=None, _i=-1, _drift=0, _x=0, _y=0, _dE=0):
        """Setup Electron Beam
        :param _e_beam_name: e-beam unique name, e.g. 'NSLS-II Low Beta Day 1' (see srwl_uti_src.py)
        :param _e_beam: e-beam structure (SRWLPartBeam instance)
        :param _i: e-beam current [A]
        :param _drift: e-beam drift length in [m] from center of straight section
        :param _x0: initial average horizontal position [m]
        :param _y0: initial average vertical position [m]
        :param _dE0: average energy deviation [GeV]
        """
        #add more parameters when/if necessary

        sIncInpElecBeam = 'Incorrect input for setting up Electron Beam structure'
        if(len(_e_beam_name) > 0):
            self.eBeam = srwl_uti_src_e_beam(_e_beam_name)
            if(self.eBeam == None):
                if((_e_beam == None) or (isinstance(_e_beam, SRWLPartBeam) == False)):
                    raise Exception(sIncInpElecBeam)
                else: self.eBeam = _e_beam
        else:
            if((_e_beam == None) or (isinstance(_e_beam, SRWLPartBeam) == False)):
                raise Exception(sIncInpElecBeam)
            else: self.eBeam = _e_beam

        if(_i > 0): self.eBeam.Iavg = _i
        if(_drift != 0): self.eBeam.drift(_drift)
        self.eBeam.partStatMom1.x = _x
        self.eBeam.partStatMom1.y = _y
        if(_dE != 0):
            elRestMassGeV = 0.51099890221e-03
            curE0 = self.eBeam.partStatMom1.gamma*self.eBeam.partStatMom1.relE0*elRestMassGeV
            self.eBeam.partStatMom1.gamma = (curE0 + _dE)/(self.eBeam.partStatMom1.relE0*elRestMassGeV)

    #------------------------------------------------------------------------
    def set_und_sin(self, _per=0.02, _len=1, _bx=0, _by=0, _phx=0, _phy=0, _sx=1, _sy=1, _zc=0, _add=0):
        """Setup magnetic field container with basic undulator having sinusoidal magnetic field
        :param _per: period length [m]
        :param _len: undulator length [m]
        :param _bx: horizontal magnetic field amplitude [m]
        :param _by: vertical magnetic field amplitude [m]
        :param _phx: initial phase of the horizontal magnetic field [rad]
        :param _phx: initial phase of the vertical magnetic field [rad]
        :param _sx: symmetry of the horizontal magnetic field vs longitudinal position 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
        :param _sy: symmetry of the vertical magnetic field vs longitudinal position
        :param _zc: longitudinal position of the undulator center
        :param _reset: add (=1) or reset (=0) the new undulator structure to the existing approximate magnetic field container
        """

        if(_add == 0):
            if(self.mag_approx != None):
                del self.mag_approx

        if(self.mag_approx == None): self.mag_approx = SRWLMagFldC()

        und = SRWLMagFldU()
        und.set_sin(_per, _len, _bx, _by, _phx, _phy, _sx, _sy)
        self.mag_approx.arMagFld.append(und)
        self.mag_approx.arXc.append(0)
        self.mag_approx.arYc.append(0)
        self.mag_approx.arZc.append(_zc)

    #------------------------------------------------------------------------
    #def set_mag_fld(self):

    #------------------------------------------------------------------------
    #def set_gsn_beam(self):

    #------------------------------------------------------------------------
    def set_optics(self, _op):
        """Setup optical element container
        :param _op: optical element container (SRWLOptC instance)
        """

        if((_op == None) or (isinstance(_op, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')
        if(self.optics != None): del self.optics
        self.optics = _op

    #------------------------------------------------------------------------
    def calc_sr_se(self, _mesh, _samp_fact=-1, _meth=2, _rel_prec=0.01, _pol=6, _int_type=0, _mag_type=1, _fname=''):
        """Calculates single-electron intensity
        :param _mesh: mesh on which the intensity has to be calculated (SRWLRadMesh instance)
        :param _samp_fact: sampling factor for adjusting nx, ny (effective if > 0)
        :param _meth: SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
        :param _rel_prec: relative precision
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
        :param _int_type: "type" of a characteristic to be extracted:
           -1- No Intensity / Electric Field components extraction is necessary (only Wavefront will be calculated)
            0- "Single-Electron" Intensity; 
            1- "Multi-Electron" Intensity; 
            2- "Single-Electron" Flux; 
            3- "Multi-Electron" Flux; 
            4- "Single-Electron" Radiation Phase; 
            5- Re(E): Real part of Single-Electron Electric Field;
            6- Im(E): Imaginary part of Single-Electron Electric Field;
            7- "Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence);
        :param _mag_type: "type" of magnetic field to use: 
            1- "Approximate", referenced by self.mag_approx; 
            2- "Accurate" (tabulated), referenced by self.mag; 
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: 1D array with (C-aligned) resulting intensity data
        """

        if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
            raise Exception('Incorrect SRWLRadMesh structure')

        depType = -1
        if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6
        if(depType < 0): Exception('Incorrect numbers of points in the mesh structure')

        if(self.eBeam == None): Exception('Electron Beam structure is not defined')

        if(_mag_type == 1):
            if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag == None): Exception('Magnetic Field is not defined')
        else: Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2): magToUse = self.mag

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)
        wfr.partBeam = self.eBeam

        zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
        zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
        npTraj = 20000 #Number of points for trajectory calculation 
        useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
        arPrecPar = [_meth, _rel_prec, zStartInteg, zEndInteg, npTraj, useTermin, _samp_fact]

        srwl.CalcElecFieldSR(wfr, 0, magToUse, arPrecPar) #calculate SR

        arI = None
        if(_int_type >= 0): 
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            arI = array(sNumTypeInt, [0]*wfr.mesh.ne*wfr.mesh.nx*wfr.mesh.ny)
            srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
            if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, wfr.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
        return wfr, arI

    #------------------------------------------------------------------------
    def calc_ur_spec_me(self, _mesh, _harm_init=1, _harm_fin=15, _prec_long=1., _prec_azim=1., _type=1, _pol=6, _fname=''):
        """Calculates multi-electron flux of undulator radiation (within fixed aperture of per unit surface), using approximate periodic magnetic field
        :param _mesh: mesh on which the intensity has to be calculated (SRWLRadMesh instance)
        :param _harm_init: initial UR spectral harmonic to be taken into account
        :param _harm_fin: final UR spectral harmonic to be taken into account
        :param _prec_long: longitudinal integration precision parameter
        :param _prec_azim: azimuthal integration precision parameter
        :param _type: calculate flux (=1) or flux per unit surface (=2)
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: 1D array with (C-aligned) resulting intensity data
        """

        if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
            raise Exception('Incorrect SRWLRadMesh structure')

        depType = -1
        if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6
        if(depType < 0): Exception('Incorrect numbers of points in the mesh structure')

        if(self.eBeam == None): Exception('Electron Beam structure is not defined')
        if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')

        stk = SRWLStokes()
        stk.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #numbers of points vs photon energy, horizontal and vertical positions
        stk.mesh = deepcopy(_mesh)

        #Precision Parampeters for Spectral Flux through a Slit:
        arPrecPar = [_harm_init, _harm_fin, _prec_long, _prec_azim, _type]

        und = self.mag_approx.arMagFld[0]
        if(isinstance(und, SRWLMagFldU) == False): raise Exception('Incorrect SRWLMagFldU (i.e. undulator) structure')

        eBeamAux = self.eBeam
        zc = self.mag_approx.arZc[0]
        if(zc != 0.):
            #since the srwl.CalcStokesUR method assumes undulator located at zc = 0,
            #we have to "move" observation plane and the longitudinal position at which e-beam parameters are defined
            stk.mesh.zStart -= zc
            eBeamAux = deepcopy(self.eBeam)
            eBeamAux.drift(-zc)

        srwl.CalcStokesUR(stk, eBeamAux, und, arPrecPar)

        arI = stk.to_int(_pol)
        if(len(_fname) > 0):
            sValName = 'Flux'
            sValUnitName = 'ph/s/.1%bw'
            if(_type == 2):
                sValName = 'Intensity'
                sValUnitName = 'ph/s/.1%bw/mm^2'
            srwl_uti_save_intens_ascii(arI, stk.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
        return arI

    #------------------------------------------------------------------------
    def calc_pow_den(self, _mesh, _prec=1, _meth=1, _z_start=0., _z_fin=0., _mag_type=1, _fname=''):
        """Calculates multi-electron flux of undulator radiation (within fixed aperture of per unit surface), using approximate periodic magnetic field
        :param _mesh: mesh (grid) on which the power density has to be calculated (SRWLRadMesh instance)
        :param _prec: precision factor for calculation of power density distribution
        :param _meth: power density computation method (1- "near field", 2- "far field")
        :param _z_start: initial longitudinal position along electron trajectory of power density distribution (effective if _s_start < _s_fin)
        :param _z_fin: final longitudinal position along electron trajectory of power density distribution (effective if _s_start < _s_fin)
        :param _mag_type: "type" of magnetic field to use: 
            1- "Approximate", referenced by self.mag_approx; 
            2- "Accurate" (tabulated), referenced by self.mag; 
        :param _fname: name of file to save the resulting data to
        :return: 1D array with (C-aligned) resulting power density data
        """

        if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
            raise Exception('Incorrect SRWLRadMesh structure')

        if(self.eBeam == None): Exception('Electron Beam structure is not defined')

        if(_mag_type == 1):
            if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag == None): Exception('Magnetic Field is not defined')
        else: Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2): magToUse = self.mag

        arPrecPar = [_prec, _meth, _z_start, _z_fin, 20000]

        stkP = SRWLStokes() #for power density vs x and y
        stkP.allocate(1, _mesh.nx, _mesh.ny) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
        stkP.mesh = deepcopy(_mesh)
        srwl.CalcPowDenSR(stkP, self.eBeam, 0, magToUse, arPrecPar) #calculate SR

        if(len(_fname) > 0):
            srwl_uti_save_intens_ascii(stkP.arS, _mesh, _fname, 0,
                                       ['', 'Horizontal Position', 'Vertical Position', 'Power Density'],
                                       _arUnits=['', 'm', 'm', 'W/mm^2'])
        return stkP.arS#, arSx, arSy

    #------------------------------------------------------------------------
    def calc_wfr_prop(self, _wfr, _pol=6, _int_type=0, _dep_type=3, _fname=''):
        """Calculates single-electron (/ fully coherent) wavefront propagation
        :param _wfr: wavefront (instance of SRWLWfr) to be propagated (and modified in place!)
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
        :param _int_type: "type" of a characteristic to be extracted:
           -1- No Intensity / Electric Field components extraction is necessary (only Wavefront will be calculated)
            0- "Single-Electron" Intensity; 
            1- "Multi-Electron" Intensity; 
            2- "Single-Electron" Flux; 
            3- "Multi-Electron" Flux; 
            4- "Single-Electron" Radiation Phase; 
            5- Re(E): Real part of Single-Electron Electric Field;
            6- Im(E): Imaginary part of Single-Electron Electric Field;
            7- "Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence);
        :param _dep_type: "type" of dependence to be extracted:
            0- vs e (photon energy or time);
            1- vs x (horizontal position or angle);
            2- vs y (vertical position or angle);
            3- vs x&y (horizontal and vertical positions or angles);
            4- vs e&x (photon energy or time and horizontal position or angle);
            5- vs e&y (photon energy or time and vertical position or angle);
            6- vs e&x&y (photon energy or time, horizontal and vertical positions or angles);
            
        :param _fname: name of file to save the resulting data to
        :return: 1D array with (C-aligned) resulting intensity data; it also modified _wfr in place
        """

        if((hasattr(self, 'optics') == False) or (isinstance(self.optics, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')
        if(isinstance(_wfr, SRWLWfr) == False):
            raise Exception('Incorrect wavefront (SRWLWfr) structure')

        srwl.PropagElecField(_wfr, self.optics)
 
        arI = None
        if(_int_type >= 0): 
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            arI = array(sNumTypeInt, [0]*_wfr.mesh.ne*_wfr.mesh.nx*_wfr.mesh.ny)
            srwl.CalcIntFromElecField(arI, _wfr, _pol, _int_type, _dep_type, _wfr.mesh.eStart, _wfr.mesh.xStart, _wfr.mesh.yStart)
            if(len(_fname) > 0):
                sValUnitName = 'ph/s/.1%bw/mm^2' #consider allowing for other units (for FEL applications)
                srwl_uti_save_intens_ascii(arI, _wfr.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', sValUnitName])
        return arI

    #------------------------------------------------------------------------
    def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None):
        """Calculates multi-electron (/ partially coherent) SR emission and wavefront propagation
        :param _mesh: mesh (grid) on which the initial wavefront has to be calculated (SRWLRadMesh instance)
        :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
        :param _sr_meth: SR Electric Field calculation method to be used (0- "manual", 1- "auto-undulator", 2- "auto-wiggler")
        :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, smaller the more accurate)
        :param _mag_type: "type" of magnetic field to use: 
            1- "Approximate", referenced by self.mag_approx; 
            2- "Accurate" (tabulated), referenced by self.mag; 
        :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
        :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
        :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
        :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1)
        :param _char: radiation characteristic to calculate:
            0- Intensity (s0);
            1- Four Stokes components;
            2- Mutual Intensity Cut vs X;
            3- Mutual Intensity Cut vs Y;
            4- Mutual Intensity Cut vs X & Y
        :param _x0: horizontal center position for mutual intensity cut calculation
        :param _y0: vertical center position for mutual intensity cut calculation
        :param _e_ph_integ: integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from _mesh
        :param _rand_meth: method for generation of pseudo-random numbers for e-beam phase-space integration:
            1- standard pseudo-random number generator
            2- Halton sequences
            3- LPtau sequences (to be implemented)
        :param _fname: name of file to save the resulting data to
        :return: 1D array with (C-aligned) resulting intensity data
        """

        if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
            raise Exception('Incorrect SRWLRadMesh structure')

        if((hasattr(self, 'eBeam') == False) or (isinstance(self.eBeam, SRWLPartBeam) == False)):
            raise Exception('Incorrect electron beam (SRWLPartBeam) structure')

        if((hasattr(self, 'optics') == False) or (isinstance(self.optics, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')

        if(_mag_type == 1):
            if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag == None): Exception('Magnetic Field is not defined')
        else: Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2): magToUse = self.mag

        return srwl_wfr_emit_prop_multi_e(
            _e_beam = self.eBeam, _mag = magToUse, _mesh = _mesh, _sr_samp_fact = _sr_samp_fact,
            _sr_meth = _sr_meth, _sr_rel_prec = _sr_rel_prec,
            _n_part_tot = _n_part_tot, _n_part_avg_proc = _n_part_avg_proc, _n_save_per = _n_save_per,
            _file_path = _fname,
            _opt_bl = self.optics,
            _pres_ang = _pres_ang, _char = _char, _x0 = _x0, _y0 = _y0,
            _e_ph_integ = _e_ph_integ, _rand_meth = _rand_meth)

    #------------------------------------------------------------------------
    def calc_all(self, _v, _op):
        """Performs setup of electron beam, magnetic field, and performs calculations according to options specified in _v
        :param _v: an object containing set of variables / options defining SR source and required calculations
        :param _op: optical element container (SRWLOptC instance) that is assumed to be set up before calling this function and eventually used for wavefront propagation
        """

        #---setup electron beam
        self.set_e_beam(
            _e_beam_name = (_v.ebm_nm + _v.ebm_nms),
            _i = _v.ebm_i,
            _drift = _v.ebm_dr,
            _x = _v.ebm_x,
            _y = _v.ebm_y,
            _dE = _v.ebm_de)

        #---setup magnetic field: undulator
        if hasattr(_v, 'und_b'):
            if hasattr(_v, 'und_bx') == False: _v.und_bx = 0
            if hasattr(_v, 'und_by') == False: _v.und_by = _v.und_b
            if hasattr(_v, 'und_phx') == False: _v.und_phx = 0
            if hasattr(_v, 'und_phy') == False: _v.und_phy = 0
            if hasattr(_v, 'und_zc') == False: _v.und_zc = 0
            self.set_und_sin(#setup approximate undulator field parameters
                _per = _v.und_per,
                _len = _v.und_len,
                _bx = _v.und_bx,
                _by = _v.und_by,
                _phx = _v.und_phx,
                _phy = _v.und_phy,
                _sx = _v.und_sx,
                _sy = _v.und_sy,
                _zc = _v.und_zc)
            
            if _v.und_b2e:
                e1 = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].get_E1(_en_elec=self.eBeam.partStatMom1.get_E(), _unit='eV')
                print('Fundamental Photon Energy:', srwl_uti_num_round(e1), 'eV') #check how it will work under IPython

            if _v.und_e2b:
                b = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].E1_2_B(_e1=_v.w_e, _en_elec=self.eBeam.partStatMom1.get_E())
                print('Magnetic Field Amplitude:', srwl_uti_num_round(b), 'T') #check how it will work under IPython

        #to program setting-up tabulated magnetic field (calculated or measured)

        #---calculate single-e spectrum vs photon energy
        if(_v.ss): 
            mesh_ss = SRWLRadMesh(
                _v.ss_ei, _v.ss_ef, _v.ss_ne,
                _v.ss_x, _v.ss_x, 1,
                _v.ss_y, _v.ss_y, 1,
                _v.op_r)
            wfr_ss, int_ss = self.calc_sr_se(
                _mesh = mesh_ss,
                _meth = _v.ss_meth,
                _rel_prec = _v.ss_prec,
                _pol = _v.ss_pol,
                _int_type = 0,
                _mag_type = _v.ss_mag,
                _fname = os.path.join(_v.fdir, _v.ss_fn) if(len(_v.ss_fn) > 0) else '')
            
        #---calculate multi-e spectrum vs photon energy
        if(_v.sm):
            mesh_sm = SRWLRadMesh(
                _v.sm_ei, _v.sm_ef, _v.sm_ne,
                _v.sm_x - 0.5*_v.sm_rx, _v.sm_x + 0.5*_v.sm_rx, _v.sm_nx,
                _v.sm_y - 0.5*_v.sm_ry, _v.sm_y + 0.5*_v.sm_ry, _v.sm_ny,
                _v.op_r)
            if(_v.sm_mag == 1):
                int_sm = self.calc_ur_spec_me(
                    _mesh = mesh_sm,
                    _harm_init = _v.sm_hi,
                    _harm_fin = _v.sm_hf,
                    _prec_long = _v.sm_prl,
                    _prec_azim = _v.sm_pra,
                    _type = _v.sm_type,
                    _pol = _v.sm_pol,
                    _fname = os.path.join(_v.fdir, _v.sm_fn) if(len(_v.sm_fn) > 0) else '')
            #else: #to program for "accurate" magnetic field

        #---calculate SR power density distribution
        if(_v.pw):
            mesh_pw = SRWLRadMesh(
                1000, 1000, 1, #dummy (not used for power density)
                _v.pw_x - 0.5*_v.pw_rx, _v.pw_x + 0.5*_v.pw_rx, _v.pw_nx,
                _v.pw_y - 0.5*_v.pw_ry, _v.pw_y + 0.5*_v.pw_ry, _v.pw_ny,
                _v.op_r)
            #intxy_pw, intx_pw, inty_pw = self.calc_pow_den(
            int_pw = self.calc_pow_den(
                _mesh = mesh_pw,
                _prec = _v.pw_pr,
                _meth = _v.pw_meth,
                _z_start = _v.pw_zst,
                _z_fin = _v.pw_zfi,
                _mag_type = _v.pw_mag,
                _fname= os.path.join(_v.fdir, _v.pw_fn) if(len(_v.pw_fn) > 0) else '')
            
        #---calculate single-e and multi-e intensity distributions (before and after wavefront propagation through a beamline)
        if(_v.si or _v.ws or _v.wm):
            if(_v.ws or _v.wm): self.set_optics(_op)
            
            ef = _v.w_e
            if(_v.w_ef > 0): ef = _v.w_ef
            mesh_w = SRWLRadMesh(
                _v.w_e, ef, _v.w_ne,
                _v.w_x - 0.5*_v.w_rx, _v.w_x + 0.5*_v.w_rx, _v.w_nx,
                _v.w_y - 0.5*_v.w_ry, _v.w_y + 0.5*_v.w_ry, _v.w_ny,
                _v.op_r)
        #---calculate single-e electric field and intensity (before wavefront propagation through a beamline)
            if(_v.si or _v.ws):
                wfr, int_w0 = self.calc_sr_se(
                    _mesh = deepcopy(mesh_w),
                    _samp_fact = _v.w_smpf,
                    _meth = _v.w_meth,
                    _rel_prec = _v.w_prec,
                    _pol = _v.si_pol,
                    _int_type = _v.si_type,
                    _mag_type = _v.w_mag,
                    _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '')
                mesh_si = deepcopy(wfr.mesh)
        #---calculate single-e electric field and intensity (after wavefront propagation through a beamline)
                if(_v.ws):
                    int_ws = self.calc_wfr_prop(
                        _wfr = wfr,
                        _pol = _v.si_pol,
                        _int_type = _v.si_type,
                        _dep_type=3, #consider adding other cases (e.g. for TD FEL calculations)
                        _fname = os.path.join(_v.fdir, _v.ws_fni) if(len(_v.ws_fni) > 0) else '')
                    mesh_ws = wfr.mesh
                    #if(len(_v.ws_fn) > 0): to implement saving single-e (/ fully coherent) wavefront data (wfr) to a file

        #---calculate multi-electron (/ partially coherent) wavefront propagation
            if(_v.wm):
                res_ipm = self.calc_wfr_emit_prop_me(
                    _mesh = mesh_w,
                    _sr_samp_fact = _v.w_smpf,
                    _sr_meth = _v.w_meth,
                    _sr_rel_prec = _v.w_prec,
                    _mag_type = _v.w_mag,
                    _n_part_tot = _v.wm_nm,
                    _n_part_avg_proc = _v.wm_na,
                    _n_save_per = _v.wm_ns,
                    _pres_ang = _v.wm_ap,
                    _char = _v.wm_ch,
                    _x0 = _v.wm_x0,
                    _y0 = _v.wm_y0,
                    _e_ph_integ = _v.wm_ei,
                    _rand_meth = _v.wm_rm,
                    _fname = os.path.join(_v.fdir, _v.wm_fni) if(len(_v.wm_fni) > 0) else None)
                #if(len(_v.wp_fn) > 0): to implement saving single-e (/ fully coherent) wavefront data (wfr) to a file

        #---plot results of all calculatiopns here (because the plotting "from the middle of the script" may hang up script execution)
        #uti_plot_init('TkAgg') #make the backend name an input option or move this to uti_plot ?
        plotOK = False
        if (_v.ss == True) and (len(_v.ss_pl) > 0):
            uti_plot1d(int_ss, [mesh_ss.eStart, mesh_ss.eFin, mesh_ss.ne], ['Photon Energy', 'Intensity', 'Intensity'], ['eV', 'ph/s/.1%bw/mm^2'])
            plotOK = True

        if (_v.sm == True) and (len(_v.sm_pl) > 0):
            sValType = 'Flux'; sValUnit = 'ph/s/.1%bw'
            if(_v.sm_type == 2):
                sValType = 'Flux per Unit Surface'; sValUnit = 'ph/s/.1%bw/mm^2'
            uti_plot1d(int_sm, [mesh_sm.eStart, mesh_sm.eFin, mesh_sm.ne], ['Photon Energy', sValType, sValType], ['eV', sValUnit])
            plotOK = True

        if (_v.pw == True) and (len(_v.pw_pl) > 0):
            if ((_v.pw_pl == 'xy') or (_v.pw_pl == 'yx') or (_v.pw_pl == 'XY') or (_v.pw_pl == 'YX')) and (_v.pw_nx > 1) and (_v.pw_ny > 1):
                uti_plot2d1d(
                    int_pw,
                    [mesh_pw.xStart, mesh_pw.xFin, mesh_pw.nx],
                    [mesh_pw.yStart, mesh_pw.yFin, mesh_pw.ny],
                    0.5*(mesh_pw.xStart + mesh_pw.xFin),
                    0.5*(mesh_pw.yStart + mesh_pw.yFin),
                    ['Horizontal Position', 'Vertical Position', 'Power Density'],
                    ['m', 'm', 'W/mm^2'], True)
            elif ((_v.pw_pl == 'x') or (_v.pw_pl == 'X')) and (_v.pw_nx > 1):
                uti_plot1d(int_pw, [mesh_pw.xStart, mesh_pw.xFin, mesh_pw.nx],
                           ['Horizontal Position', 'Power Density', 'Power Density'],
                           ['m', 'W/mm^2'])
            elif ((_v.pw_pl == 'y') or (_v.pw_pl == 'Y')) and (_v.pw_ny > 1):
                uti_plot1d(int_pw, [mesh_pw.yStart, mesh_pw.yFin, mesh_pw.ny],
                           ['Vertical Position', 'Power Density', 'Power Density'],
                           ['m', 'W/mm^2'])
            plotOK = True

        if _v.si and (len(_v.si_pl) > 0):
            if _v.si_pl in ['xy', 'yx', 'XY', 'YX']:
                #print('testing', _v.si_pl)
                uti_plot2d1d(
                    int_w0,
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', 'Intensity Before Propagation'],
                    ['m', 'm', 'ph/s/.1%bw/mm^2'], #add other units for FEL
                    True)
                plotOK = True

        if _v.ws and (len(_v.ws_pl) > 0):
            if (_v.ws_pl == 'xy') or (_v.ws_pl == 'yx') or (_v.ws_pl == 'XY') or (_v.ws_pl == 'YX'):
                uti_plot2d1d(
                    int_w0,
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', 'Intensity Before Propagation'],
                    ['m', 'm', 'ph/s/.1%bw/mm^2'], #add other units for FEL
                    True)
                uti_plot2d1d(
                    int_ws,
                    [mesh_ws.xStart, mesh_ws.xFin, mesh_ws.nx],
                    [mesh_ws.yStart, mesh_ws.yFin, mesh_ws.ny],
                    0, #0.5*(mesh_ws.xStart + mesh_ws.xFin),
                    0, #0.5*(mesh_ws.yStart + mesh_ws.yFin),
                    ['Horizontal Position', 'Vertical Position', 'Intensity After Propagation'],
                    ['m', 'm', 'ph/s/.1%bw/mm^2'], #add other units for FEL
                    True)

            #to continue here
                
            plotOK = True

        if(plotOK): uti_plot_show()

#****************************************************************************
def srwl_uti_parse_str2list(_str):
    """Parse _str, such as '1 2 3' or '1,2,3', to list
    """
    sLoc = copy(_str)
    sLoc = sLoc.replace('[', ''); sLoc = sLoc.replace(']', '')
    sLoc = sLoc.replace('(', ''); sLoc = sLoc.replace(')', '')
    sLoc = sLoc.replace('{', ''); sLoc = sLoc.replace('}', '')
    
    resList = []
    if(',' in sLoc): resList = sLoc.split(',')
    elif(';' in sLoc): resList = sLoc.split(';')
    else: resList = sLoc.split(' ')

    for i in range(len(resList)):
        resList[i] = float(resList[i])

    return resList

#****************************************************************************
def srwl_uti_parse_options(_descr):
    """Set and parse command-prompt options from a compact description provided in _descr
    :param _descr: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    """

    p = optparse.OptionParser()
    nOpt = len(_descr)

    listOptNamesPostParse = []
    for i in range(nOpt):
        curOpt = _descr[i]
        
        sTypeShort = curOpt[1]
        sType = 'string'
        if(sTypeShort == 'f'): sType = 'float'
        elif(sTypeShort == 'i'): sType = 'int'        
        #elif(sTypeShort == 's'): sType = 'string'

        sAct = 'store'
        if(len(curOpt) > 4): sAct = curOpt[4]

        defVal = curOpt[2]
        
        optIsList = False
        if(isinstance(defVal, list) or isinstance(defVal, array)): optIsList = True

        if(optIsList):
            sType = 'string'
            listOptNamesPostParse.append(curOpt[0])

        if(len(sTypeShort) <= 0):
            p.add_option('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
        else:
            p.add_option('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)

    v, args = p.parse_args()

    #"post-parsing" list-type options
    for i in range(len(listOptNamesPostParse)):
        curOptName = listOptNamesPostParse[i]
        valCurOpt = getattr(v, curOptName)

        if((isinstance(valCurOpt, list) == False) and (isinstance(valCurOpt, array) == False)):
            parsedVal = srwl_uti_parse_str2list(valCurOpt)
            setattr(v, curOptName, parsedVal)
    
    return v
