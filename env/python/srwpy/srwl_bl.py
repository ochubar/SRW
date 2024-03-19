# -*- coding: utf-8 -*-
#############################################################################
# SRWLib SR Beamline Base Class
# Contains a set of member objects and functions for simulating basic operation and characteristics
# of a complete user beamline in a synchrotron radiation source.
# Under development!!!
# Authors/Contributors: O.C., Maksim Rakitin, An He
# v 0.13
#############################################################################

from __future__ import print_function #Python 2.7 compatibility

try: #OC15112022
    from .srwlib import *
    from .srwl_uti_mag import *
    from .srwl_uti_und import *
    from .uti_plot import *
    from . import uti_math
    from . import uti_io
except:
    from srwlib import *
    from srwl_uti_mag import *
    from srwl_uti_und import *
    from uti_plot import *
    import uti_math
    import uti_io
#from srwlib import *
#from srwl_uti_mag import *
#from srwl_uti_und import *
#from uti_plot import *
#import uti_math
#import uti_io

#import optparse #MR081032016 #Consider placing import argparse here
import time
import re #OC23052020
import zipfile #OC23052020

try: #OC05042018
    import cPickle as pickle
except:
    import pickle

#****************************************************************************
# Global Constants

_Pi = 3.14159265358979
_ElCh = 1.60217646263E-19 #1.602189246E-19 #Electron Charge [Q]
_ElMass_kg = 9.1093818872E-31 #9.10953447E-31 #Electron Mass in [kg]
_ElMass_MeV = 0.51099890221 #Electron Mass in [MeV]
_LightSp = 2.9979245812E+08 #Speed of Light [m/c]
_Light_eV_mu = 1.23984186 #Wavelength <-> Photon Energy conversion constant ([um] <-> [eV])
_PlanckConst_eVs = 4.13566766225E-15 #Planck constant in [eV*s]

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
    def __init__(self, _name='', _e_beam=None, _mag_approx=None, _mag=None, _gsn_beam=None, _op=None): #, _det=None):
        """
        :param _name: beamline name string
        :param _e_beam: electron beam (SRWLPartBeam instance)
        :param _mag_approx: approximate magnetic field container (SRWLMagFldC instance)
        :param _mag: accurate magnetic field container (SRWLMagFldC instance)
        :param _gsn_beam: coherent Gaussian beam (SRWLGsnBm instance)
        :param _op: optics sequence (SRWLOptC instance)
        :param _det: detector (SRWLDet instance)
        """

        if _e_beam is not None:
            if(isinstance(_e_beam, SRWLPartBeam) == False):
                raise Exception("Incorrect Electron Beam structure")
        if _mag_approx is not None:
            if(isinstance(_mag_approx, SRWLMagFldC) == False):
                raise Exception("Incorrect Magnetic Field Container structure")
        if _mag is not None:
            if(isinstance(_mag, SRWLMagFldC) == False):
                raise Exception("Incorrect Magnetic Field Container structure")
        if _gsn_beam is not None:
            if(isinstance(_gsn_beam, SRWLGsnBm) == False):
                raise Exception("Incorrect Gaussian Beam structure") 
        if _op is not None:
            if(isinstance(_op, SRWLOptC) == False):
                raise Exception("Incorrect Optical Container structure") 
        #if _det != None:
        #    if(isinstance(_det, SRWLDet) == False):
        #        raise Exception("Incorrect Detector structure") 

        self.name = _name
        self.eBeam = SRWLPartBeam() if _e_beam is None else _e_beam
        self.mag_approx = _mag_approx
        self.mag = _mag
        self.gsnBeam = _gsn_beam
        self.optics = _op
        #self.detector = _det

    #------------------------------------------------------------------------
##    def set_e_beam(self, _e_beam_name='', _e_beam=None, _i=-1, _sig_e=-1, _emit_x=-1, _emit_y=-1, _drift=0, _x=0, _y=0, _xp=0, _yp=0, _dE=0):
##        """Setup Electron Beam.
##        NOTE: The beam is assumed to be first set up at z = 0 longitudinal position, and then it is propagated according to _drift (if _drift != 0.)
##        :param _e_beam_name: e-beam unique name, e.g. 'NSLS-II Low Beta Day 1' (see srwl_uti_src.py)
##        :param _e_beam: e-beam structure (SRWLPartBeam instance)
##        :param _i: e-beam current [A]
##        :param _sig_e: e-beam relative energy spread
##        :param _emit_x: e-beam horizontal emittance
##        :param _emit_y: e-beam vertical emittance
##        :param _drift: e-beam drift length in [m] from center of straight section
##        :param _x: initial average horizontal position [m]
##        :param _y: initial average vertical position [m]
##        :param _xp: initial average horizontal angle [m]
##        :param _yp: initial average vertical angle [m]
##        :param _dE0: average energy deviation [GeV]
##        """
##        #add more parameters when/if necessary
##
##        if(_sig_e < 0.): _sig_e = None
##        if(_emit_x < 0.): _emit_x = None
##        if(_emit_y < 0.): _emit_y = None
##
##        sIncInpElecBeam = 'Incorrect input for setting up Electron Beam structure'
##        if(len(_e_beam_name) > 0):
##            self.eBeam = srwl_uti_src_e_beam(_e_beam_name, _sig_e=_sig_e, _emit_x=_emit_x, _emit_y=_emit_y)
##            if(self.eBeam == None):
##                if((_e_beam == None) or (isinstance(_e_beam, SRWLPartBeam) == False)):
##                    raise Exception(sIncInpElecBeam)
##                else: self.eBeam = _e_beam
##        else:
##            if((_e_beam == None) or (isinstance(_e_beam, SRWLPartBeam) == False)):
##                raise Exception(sIncInpElecBeam)
##            else: self.eBeam = _e_beam
##
##        #OC: Add Twiss parameters and 2nd order moments to function arguments and program logic of switching bw these definitions
##        #OC: consider treating _sig_e=-1, _emit_x=-1, _emit_y=-1 is defined, in all cases!
##
##        if(_i > 0): self.eBeam.Iavg = _i
##        if(_drift != 0): self.eBeam.drift(_drift)
##        self.eBeam.partStatMom1.x = _x
##        self.eBeam.partStatMom1.y = _y
##        self.eBeam.partStatMom1.xp = _xp
##        self.eBeam.partStatMom1.yp = _yp
##        
##        if(_dE != 0):
##            elRestMassGeV = 0.51099890221e-03
##            curE0 = self.eBeam.partStatMom1.gamma*self.eBeam.partStatMom1.relE0*elRestMassGeV
##            self.eBeam.partStatMom1.gamma = (curE0 + _dE)/(self.eBeam.partStatMom1.relE0*elRestMassGeV)

    #------------------------------------------------------------------------
    def set_e_beam(self, _e_beam_name='', _e_beam=None, _i=-1, _ens=-1, _emx=-1, _emy=-1, _dr=0, _x=0, _y=0, _xp=0, _yp=0, _e=None, _de=0,
                   _betax=None, _alphax=None, _etax=None, _etaxp=None, _betay=None, _alphay=None, _etay=0, _etayp=0,
                   _sigx=None, _sigxp=None, _mxxp=None, _sigy=None, _sigyp=None, _myyp=None):
        """Setup Electron Beam

        :param _e_beam_name: e-beam unique name, e.g. 'NSLS-II Low Beta Day 1' (see srwl_uti_src.py)
        :param _e_beam: e-beam structure (SRWLPartBeam instance)
        :param _i: e-beam current [A]
        :param _ens: e-beam relative energy spread
        :param _emx: e-beam horizontal emittance
        :param _emy: e-beam vertical emittance
        :param _dr: e-beam drift length in [m] from center of straight section
        :param _x: initial average horizontal position [m]
        :param _y: initial average vertical position [m]
        :param _xp: initial average horizontal angle [m]
        :param _yp: initial average vertical angle [m]
        :param _e: energy [GeV]
        :param _de: average energy deviation [GeV]

        #MR28092016 - added parameters to define the beam explicitly:
        # Parameters for SRWLPartBeam.from_Twiss():
        # def from_Twiss(self, _Iavg=0, _e=0, _sig_e=0, _emit_x=0, _beta_x=0, _alpha_x=0, _eta_x=0, _eta_x_pr=0, _emit_y=0, _beta_y=0, _alpha_y=0, _eta_y=0, _eta_y_pr=0):
        :param _betax: horizontal beta-function [m]
        :param _alphax: horizontal alpha-function [rad]
        :param _etax: horizontal dispersion function [m]
        :param _etaxp: horizontal dispersion function derivative [rad]
        :param _betay: vertical beta-function [m]
        :param _alphay: vertical alpha-function [rad]
        :param _etay: vertical dispersion function [m]
        :param _etayp: vertical dispersion function derivative [rad]

        # Parameters for SRWLPartBeam.from_RMS():
        # def from_RMS(self, _Iavg=0, _e=0, _sig_e=0, _sig_x=0, _sig_x_pr=0, _m_xx_pr=0, _sig_y=0, _sig_y_pr=0, _m_yy_pr=0):
        :param _sigx: horizontal RMS size [m]
        :param _sigxp: horizontal RMS divergence [rad]
        :param _mxxp: <(x-<x>)(x'-<x'>)> [m]
        :param _sigy: vertical RMS size [m]
        :param _sigyp: vertical RMS divergence [rad]
        :param _myyp: <(y-<y>)(y'-<y'>)> [m]
        """
        #add more parameters when/if necessary

        #print('In set_e_beam: x=', _x, ' xp=', _xp)

        varParamStd = srwl_uti_std_options()
        help_dict = {}
        for v in varParamStd:
            help_dict[v[0]] = '{}{}'.format(v[3][0].upper(), v[3][1:])

        def check_positive(d):
            for k, v in d.items():
                if v is None or v < 0:
                    return False, k
            return True, None

        if(_ens < 0.): _ens = None
        if(_emx < 0.): _emx = None
        if(_emy < 0.): _emy = None

        sIncInpElecBeam = 'Incorrect input for setting up Electron Beam structure'

        eBeamWasSetFromDB = False #OC28092016
        if len(_e_beam_name) > 0:
            self.eBeam = srwl_uti_src_e_beam(_e_beam_name, _sig_e=_ens, _emit_x=_emx, _emit_y=_emy)
            eBeamWasSetFromDB = True
            if self.eBeam is None:
                if (_e_beam is None) or (isinstance(_e_beam, SRWLPartBeam) is False):
                    # raise ValueError(sIncInpElecBeam)
                    raise ValueError('The beam name "{}" was not found in the database and _e_beam is empty'.format(_e_beam_name))
                else:
                    self.eBeam = _e_beam

        eBeamWasSetFromInObj = False
        if((eBeamWasSetFromDB is False) and (isinstance(_e_beam, SRWLPartBeam) is True)):
            self.eBeam = _e_beam
            eBeamWasSetFromInObj = True

        if((eBeamWasSetFromDB is False) and (eBeamWasSetFromInObj is False)): #Try to set up e-beam from input Twiss params or Moments
            
            #beam_inputs = {'_i': _i, '_e': _e, '_ens': _ens}
            #beam_check, beam_bad_var = check_positive(beam_inputs)
            #if beam_check: 
            ##OC15092017 removed the above check because for partially-coherent calculations with Gaussian beam _i, _e, _ens can be 0 or undefined
            twiss_inputs = {'_emx': _emx, '_betax': _betax, '_etax': _etax,
                            '_emy': _emy, '_betay': _betay, '_etay': _etay}
            moments_inputs = {'_sigx': _sigx, '_sigxp': _sigxp,
                              '_sigy': _sigy, '_sigyp': _sigyp}
            twiss_check, twiss_bad_var = check_positive(twiss_inputs)
            if((_alphax is None) or (_alphay is None)): twiss_check = False #OC28092016
            if((_etaxp is None) or (_etayp is None)): twiss_check = False

            moments_check, moments_bad_var = check_positive(moments_inputs)
            if((_mxxp is None) or (_myyp is None)): moments_check = False #OC28092016
                    
            if twiss_check:
                # Defined by Twiss parameters:
                self.eBeam.from_Twiss(
                    _Iavg=_i, _e=_e, _sig_e=_ens,
                    _emit_x=_emx, _beta_x=_betax, _alpha_x=_alphax, _eta_x=_etax, _eta_x_pr=_etaxp,
                    _emit_y=_emy, _beta_y=_betay, _alpha_y=_alphay, _eta_y=_etay, _eta_y_pr=_etayp,
                )
            elif moments_check:
                # Defined by Moments:
                self.eBeam.from_RMS(
                    _Iavg=_i, _e=_e, _sig_e=_ens,
                    _sig_x=_sigx, _sig_x_pr=_sigxp, _m_xx_pr=_mxxp,
                    _sig_y=_sigy, _sig_y_pr=_sigyp, _m_yy_pr=_myyp,
                )
            else:
                # raise ValueError(sIncInpElecBeam)
                err_msg = 'Twiss and/or Moments parameters are not set correctly:\n  - {} ({}): {}\n  - {} ({}): {}\n'
                raise ValueError(err_msg.format(
                    help_dict['ebm{}'.format(twiss_bad_var)], twiss_bad_var, twiss_inputs[twiss_bad_var],
                    help_dict['ebm{}'.format(moments_bad_var)], moments_bad_var, moments_inputs[moments_bad_var],
                ))
            #else:
            #    err_msg = 'Beam parameters are not set correctly:\n  - {} ({}): {}\n'
            #    raise ValueError(err_msg.format(
            #        help_dict['ebm{}'.format(beam_bad_var)], beam_bad_var, beam_inputs[beam_bad_var],
            #    ))

        else: #Allow overriding some 2nd order moments if eBeamWasSetFromDB or eBeamWasSetFromInObj
           
            if((_ens is not None) and (_ens > 0)): self.eBeam.arStatMom2[10] = _ens*_ens
            if((_sigx is not None) and (_sigx > 0)): self.eBeam.arStatMom2[0] = _sigx*_sigx
            if((_sigxp is not None) and (_sigxp > 0)): self.eBeam.arStatMom2[2] = _sigxp*_sigxp
            if((_mxxp is not None) and (_mxxp != 1.e+23)): self.eBeam.arStatMom2[1] = _mxxp
            if((_sigy is not None) and (_sigy > 0)): self.eBeam.arStatMom2[3] = _sigy*_sigy
            if((_sigyp is not None) and (_sigyp > 0)): self.eBeam.arStatMom2[5] = _sigyp*_sigyp
            if((_myyp is not None) and (_myyp != 1.e+23)): self.eBeam.arStatMom2[4] = _myyp

        #Allow applying drift and overriding 1st order moments in any case
        #if((_dr is not None) and (_dr != 0)): self.eBeam.drift(_dr) #OC16032017: moved to after the following lines

        if((_i is not None) and (_i > 0)): self.eBeam.Iavg = _i
        if(_x is not None): self.eBeam.partStatMom1.x = _x
        if(_y is not None): self.eBeam.partStatMom1.y = _y
        if(_xp is not None): self.eBeam.partStatMom1.xp = _xp
        if(_yp is not None): self.eBeam.partStatMom1.yp = _yp
        
        if((_dr is not None) and (_dr != 0)): self.eBeam.drift(_dr) #OC16032017 (see above)

        if((_de is not None) and (_de != 0)):
            elRestMassGeV = 0.51099890221e-03
            curE0 = self.eBeam.partStatMom1.gamma*self.eBeam.partStatMom1.relE0*elRestMassGeV
            self.eBeam.partStatMom1.gamma = (curE0 + _de)/(self.eBeam.partStatMom1.relE0*elRestMassGeV)

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
            if(hasattr(self, 'mag_approx')): #OC16122018
                if(self.mag_approx is not None):
                    del self.mag_approx

        #if(self.mag_approx is None): self.mag_approx = SRWLMagFldC()
        if((not hasattr(self, 'mag_approx')) or (self.mag_approx is None)): #OC16122018
            #print('Setting up mag_approx')
            self.mag_approx = SRWLMagFldC()

        und = SRWLMagFldU()
        und.set_sin(_per, _len, _bx, _by, _phx, _phy, _sx, _sy)
        self.mag_approx.arMagFld.append(und)
        self.mag_approx.arXc.append(0)
        self.mag_approx.arYc.append(0)
        self.mag_approx.arZc.append(_zc)

        #print(self.mag_approx.arVx)
        
        return self.mag_approx

    #------------------------------------------------------------------------
    def set_mag_multipole(self, _bx=0, _by=1., _gn=0, _gs=0, _len=1.5, _led=0, _r=0, _zc=0, _add=0):
        """Setup magnetic field container with basic dipole / quadrupole magnet
        :param _bx: horizontal magnetic field [m]
        :param _by: vertical magnetic field [m]
        :param _gn: magnetic field gradient of normal quad [m]
        :param _gs: magnetic field gradient of skew quad [m]
        :param _len: magnet length [m]
        :param _led: "soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed [m]
        :param _zc: longitudinal position of the magnet center
        :param _add: add (=1) or reset (=0) the new magnet structure to the existing approximate magnetic field container
        """
        
        if(_add == 0):
            if(self.mag_approx is not None):
                del self.mag_approx

        if(self.mag_approx is None): self.mag_approx = SRWLMagFldC()

        if(_bx != 0):
            dipBx = SRWLMagFldM(_G=_bx, _m=1, _n_or_s='s', _Leff=_len, _Ledge=_led, _R=_r) #?
            self.mag_approx.arMagFld.append(dipBx)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        if(_by != 0):
            dipBy = SRWLMagFldM(_G=_by, _m=1, _n_or_s='n', _Leff=_len, _Ledge=_led, _R=_r)
            self.mag_approx.arMagFld.append(dipBy)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        if(_gn != 0):
            quadN = SRWLMagFldM(_G=_gn, _m=2, _n_or_s='n', _Leff=_len, _Ledge=_led, _R=_r)
            self.mag_approx.arMagFld.append(quadN)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        if(_gs != 0):
            quadS = SRWLMagFldM(_G=_gs, _m=2, _n_or_s='s', _Leff=_len, _Ledge=_led, _R=_r)
            self.mag_approx.arMagFld.append(quadS)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        #print('At the end of set_mag_multipole, mag_approx:', self.mag_approx) #DEBUG
        return self.mag_approx

    #------------------------------------------------------------------------
    def set_mag_kick(self, _angx=1., _angy=0, _len=1., _led=0, _zc=0, _add=0):
        """Setup magnetic field container with basic dipole "kick" magnet
        :param _angx: horizontal kick angle [rad]
        :param _angy: vertical kick angle [rad]
        :param _len: magnet length [m]
        :param _led: "soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed [m]
        :param _zc: longitudinal position of the magnet center
        :param _add: add (=1) or reset (=0) the new magnet structure to the existing approximate magnetic field container
        """
        
        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined (it is required for defining kick magnet parameters from angle)')
        elEn_GeV = self.eBeam.partStatMom1.get_E(_unit='GeV')

        if(_add == 0):
            if(self.mag_approx is not None): del self.mag_approx

        if(_angx != 0.): 
            kickMag = srwl_mag_kick(_el_en=elEn_GeV, _ang=_angx, _x_or_y='x', _len=_len, _led=_led)
            self.mag_approx.arMagFld.append(kickMag)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        if(_angy != 0.): 
            kickMag = srwl_mag_kick(_el_en=elEn_GeV, _ang=_angy, _x_or_y='y', _len=_len, _led=_led)
            self.mag_approx.arMagFld.append(kickMag)
            self.mag_approx.arXc.append(0)
            self.mag_approx.arYc.append(0)
            self.mag_approx.arZc.append(_zc)

        return self.mag_approx

    #------------------------------------------------------------------------
    def set_und_tab_total(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m', _per=0.02, _c1=0, _c2=0, _a=0, _dg_by_len=0, _y0=0, _yp=0):
    #def set_und_tab(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m', _per=0.02, _c1=0, _c2=0, _a=0, _dg_by_len=0, _y0=0, _yp=0):
        """Setup magnetic field container with magnetic measurements or calculation data interpolated for given gap and phase
        :param _gap: magnetic gap [mm] for which the field should be set up
        :param _ph_mode: type of phase (shift of magnet arrays) motion
        :param _phase: shift of magnet arrays [mm] for which the field should be set up
        :param _zc: center position [m]
        :param _interp_ord: order of interpolation: 1- (bi-)linear, 2- (bi-)quadratic, 3- (bi-)cubic
        :param _meas_or_calc: use magnetic measurements ('m') or calculation ('c') data
        :param _per: undulator period [m]
        :param _c1: constant defining (approximate) undulator field dependence on gap (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
        :param _c2: constant defining (approximate) undulator field dependence on gap (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
        :param _a: constant defining (approximate) undulator field dependence on vertical position (i.e. a in cosh(2*Pi*a*y/per)
        :param _dg_by_len: gap taper (exit minus entrance) divided by undulator length
        :param _y0: vertical electron position in the center of undulator relative to undulator median plane [m]
        :param _dydz: vertical electron angle in the center of undulator relative to undulator median plane [rad]
        """

        fPathSum = ''
        if(_meas_or_calc == 'm'):
            if(hasattr(self, 'dir_magn_meas') and hasattr(self, 'fn_magn_meas_sum')):
                fPathSum = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
            else: raise Exception('No magnetic measurements data are supplied')
        elif(_meas_or_calc == 'c'):
            raise Exception('No magnetic calculation data are supplied')
        
        f = open(fPathSum, 'r')
        lines = f.readlines() #read-in all lines
        nRows = len(lines)

        strSep = '\t'
        arGaps = []; arPhases = []; arMagFld3D = []
        arXc = []; arYc = []; arZc = []
        arCoefBx = []; arCoefBy = []

        phaseIsVar = False
        phasePrev = None
        #print('Setting up tabulated magnetic field')
        #print('Starting to read-in', nRows, 'files')
        
        for i in range(nRows):
            curLine = lines[i]
            curLineParts = curLine.split(strSep)
            curLenLineParts = len(curLineParts)
            if(curLenLineParts >= 4):
                curPhaseMode = curLineParts[1]
                if(curPhaseMode != _ph_mode): continue

                arGaps.append(float(curLineParts[0]))

                curPhase = float(curLineParts[2])
                if((phasePrev is not None) and (curPhase != phasePrev)): phaseIsVar = True
                arPhases.append(curPhase)
                phasePrev = curPhase
                
                #curFileName = curLineParts[3]
                curFileName = curLineParts[3].strip() #MR13012017
                #print(curFileName)
                
                if(len(curFileName) > 0):
                    curFilePath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, curFileName)
                    #print(curFilePath)
                    curFldCnt = srwl_uti_read_mag_fld_3d(curFilePath, '#')
                    arMagFld3D.append(curFldCnt.arMagFld[0])
                    arXc.append(curFldCnt.arXc[0])
                    arYc.append(curFldCnt.arYc[0])
                    arZc.append(curFldCnt.arZc[0] + _zc)

                if(curLenLineParts >= 6):
                    arCoefBx.append(float(curLineParts[4]))
                    arCoefBy.append(float(curLineParts[5]))
        f.close()
        #print('Done reading-in magnetic field files')

        fldCnt = SRWLMagFldC(arMagFld3D, array('d', arXc), array('d', arYc), array('d', arZc))
        nElem = len(arMagFld3D)
        if((nElem != len(arGaps)) or (nElem != len(arPhases))):
            raise Exception('Inconsistent magnetic field data summary file')
        
        fldCnt.arPar1 = array('d', arGaps)
        fldCnt.arPar2 = array('d', arPhases)
        if(len(arCoefBx) == nElem): fldCnt.arPar3 = array('d', arCoefBx)
        if(len(arCoefBy) == nElem): fldCnt.arPar4 = array('d', arCoefBy)

        #print('  Gaps:', fldCnt.arPar1)
        #print('  Phase Mode:', _ph_mode)
        #print('  Phases:', fldCnt.arPar2)
        #if(phaseIsVar is True): print('  Phase is variable')
        #else: print('  Phase is constant')

        fldCntRes = SRWLMagFldC(arMagFld3D[0], arXc[0], arYc[0], arZc[0])

        numDims = 1
        if(phaseIsVar): numDims += 1
        precPar = [numDims, _gap, _phase, _interp_ord]
        #precPar = [1, _gap, _phase, _interp_ord]
        self.mag = srwl.CalcMagnField(fldCntRes, fldCnt, precPar)

        if((_dg_by_len != 0.) or (_y0 != 0.) or (_yp != 0.)):
            self.mag.arMagFld[0] = srwl_und_fld_1d_mis(self.mag.arMagFld[0], _per, _dg_by_len, _c1, _c2, 0.001*_gap, _a, _y0, _yp)
        
        return self.mag

    #------------------------------------------------------------------------
    def set_und_tab(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m', _per=0.02, _c1=0, _c2=0, _a=0, _dg_by_len=0, _y0=0, _yp=0, _fn_magn_sum=None):
    #def set_und_tab(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m', _per=0.02, _c1=0, _c2=0, _a=0, _dg_by_len=0, _y0=0, _yp=0):
    #def set_und_tab_faster(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m', _per=0.02, _c1=0, _c2=0, _a=0, _dg_by_len=0, _y0=0, _yp=0):
        """Setup magnetic field container with magnetic measurements or calculation data interpolated for given gap and phase. Faster version, reading only requiired field files.
        :param _gap: magnetic gap [mm] for which the field should be set up
        :param _ph_mode: type of phase (shift of magnet arrays) motion
        :param _phase: shift of magnet arrays [mm] for which the field should be set up
        :param _zc: center position [m]
        :param _interp_ord: order of interpolation: 1- (bi-)linear, 2- (bi-)quadratic, 3- (bi-)cubic
        :param _meas_or_calc: use magnetic measurements ('m') or calculation ('c') data
        :param _per: undulator period [m]
        :param _c1: constant defining (approximate) undulator field dependence on gap (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
        :param _c2: constant defining (approximate) undulator field dependence on gap (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))
        :param _a: constant defining (approximate) undulator field dependence on vertical position (i.e. a in cosh(2*Pi*a*y/per)
        :param _dg_by_len: gap taper (exit minus entrance) divided by undulator length
        :param _y0: vertical electron position in the center of undulator relative to undulator median plane [m]
        :param _yp: vertical electron angle in the center of undulator relative to undulator median plane [rad]
        :param _fn_magn_sum: tabulated magnetic field data summary file name (optional)
        """

        fPathSum = ''
        if(_meas_or_calc == 'm'):
            #errStrNoMagnMeasData = 'No magnetic measurements data are supplied'
            #if(hasattr(self, 'dir_magn_meas') and hasattr(self, 'fn_magn_meas_sum')):
            #    fPathSum = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
            #else: raise Exception('No magnetic measurements data are supplied')
            #OC21012023
            if hasattr(self, 'dir_magn_meas'):
                if _fn_magn_sum is not None:
                    fPathSum = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, _fn_magn_sum)
                elif hasattr(self, 'fn_magn_meas_sum'):
                    fPathSum = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
            if not os.path.exists(fPathSum): raise Exception('No magnetic measurements data are supplied')

        elif(_meas_or_calc == 'c'):
            raise Exception('No magnetic calculation data are supplied') #??
        
        f = open(fPathSum, 'r')
        lines = f.readlines() #read-in all lines
        nRows = len(lines)

        strSep = '\t'
        arGaps = []; arPhases = []; arFieldFileNames = []; arIndReqFiles = []
        #arXc = []; arYc = []; arZc = []
        arCoefBx = []; arCoefBy = []

        phaseIsVar = False
        phasePrev = None
        #print('Setting up tabulated magnetic field')

        #Read-in summary file and setup aux. arrays
        for i in range(nRows):
            curLine = lines[i]
            curLineParts = curLine.split(strSep)
            curLenLineParts = len(curLineParts)
            if(curLenLineParts >= 4):
                curPhaseMode = curLineParts[1]
                if(curPhaseMode != _ph_mode): continue

                curFileName = curLineParts[3].strip() #MR13012017
                #print(curFileName)

                if(len(curFileName) > 0):
                    arFieldFileNames.append(curFileName)

                    arGaps.append(float(curLineParts[0]))

                    curPhase = float(curLineParts[2])
                    if((phasePrev is not None) and (curPhase != phasePrev)): phaseIsVar = True
                    arPhases.append(curPhase)
                    phasePrev = curPhase

                    arIndReqFiles.append(-1) #Array to be filled-out from C

                if(curLenLineParts >= 6):
                    arCoefBx.append(float(curLineParts[4]))
                    arCoefBy.append(float(curLineParts[5]))
        f.close()

        numDims = 1
        if(phaseIsVar): numDims += 1
        meshIsRect = 0
        precPar = [numDims, _gap, _phase, _interp_ord, meshIsRect] #NOTE: these values can be modified by srwl.UtiUndFindMagFldInterpInds

        #OCTEST
        #print('precPar for UtiUndFindMagFldInterpInds:', precPar)
        
        nIndReq = srwl.UtiUndFindMagFldInterpInds(arIndReqFiles, arGaps, arPhases, precPar) #to implement

        #OCTEST
        #print('precPar after UtiUndFindMagFldInterpInds:', precPar)
        #print(nIndReq)
        #print(arIndReqFiles)
        #print(arFieldFileNames)

        if((nIndReq <= 0) or (nIndReq > nRows)):
            raise Exception('Inconsistent magnetic field data summary file')

        arResMagFld3D = []
        arResGaps = []; arResPhases = []
        arResXc = []; arResYc = []; arResZc = []
        arResCoefBx = []; arResCoefBy = []
        for j in range(nIndReq):
            curInd = arIndReqFiles[j]

            curFileName = arFieldFileNames[curInd]
            curFilePath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, curFileName)

            #OCTEST
            #print(curFilePath)
            
            curFldCnt = srwl_uti_read_mag_fld_3d(curFilePath, '#')
            if(curFldCnt is not None):
                arResMagFld3D.append(curFldCnt.arMagFld[0])
                arResXc.append(curFldCnt.arXc[0])
                arResYc.append(curFldCnt.arYc[0])
                arResZc.append(curFldCnt.arZc[0] + _zc)

                arResGaps.append(arGaps[curInd])
                arResPhases.append(arPhases[curInd])
                
                if(len(arCoefBx) > curInd): arResCoefBx.append(arCoefBx[curInd])
                if(len(arCoefBy) > curInd): arResCoefBy.append(arCoefBy[curInd])

        #OCTEST
        #print('Gaps:', arResGaps)
        #print('Phases:', arResPhases)
        
        fldResCnt = SRWLMagFldC(arResMagFld3D, array('d', arResXc), array('d', arResYc), array('d', arResZc))
        #nElem = len(arMagFld3D)
        #if((nElem != len(arGaps)) or (nElem != len(arPhases))):
        #    raise Exception('Inconsistent magnetic field data summary file')
        
        fldResCnt.arPar1 = array('d', arResGaps)
        fldResCnt.arPar2 = array('d', arResPhases)
        if(len(arResCoefBx) == nIndReq): fldResCnt.arPar3 = array('d', arResCoefBx)
        if(len(arResCoefBy) == nIndReq): fldResCnt.arPar4 = array('d', arResCoefBy)

        #print('  Gaps:', fldCnt.arPar1)
        #print('  Phase Mode:', _ph_mode)
        #print('  Phases:', fldCnt.arPar2)
        #if(phaseIsVar is True): print('  Phase is variable')
        #else: print('  Phase is constant')

        #fldCntFinRes = SRWLMagFldC(arResMagFld3D[0], arResXc[0], arResYc[0], arResZc[0])
        fldCntFinRes = SRWLMagFldC(copy(arResMagFld3D[0]), arResXc[0], arResYc[0], arResZc[0])

        precPar.append(1) #this means that search for necessary subset of indexes should not be done, i.e. the field container is assumed to include only the fields necessary for the interpolaton

        #OCTEST
        #print('precPar for CalcMagnField:', precPar)

        #OC21012023
        curMagCnt = srwl.CalcMagnField(fldCntFinRes, fldResCnt, precPar)
        curCntWasNotAdded = True
        if hasattr(self, 'mag'):
            if self.mag is not None: 
                if isinstance(self.mag, SRWLMagFldC): 
                    self.mag.add(curMagCnt)
                    curCntWasNotAdded = False
        if curCntWasNotAdded: self.mag = curMagCnt
        #self.mag = srwl.CalcMagnField(fldCntFinRes, fldResCnt, precPar)

        if((_dg_by_len != 0.) or (_y0 != 0.) or (_yp != 0.)):
            self.mag.arMagFld[0] = srwl_und_fld_1d_mis(self.mag.arMagFld[0], _per, _dg_by_len, _c1, _c2, 0.001*_gap, _a, _y0, _yp)
            #Add mis-steering to other und. segments
        
        return self.mag

    #------------------------------------------------------------------------
    def set_und_per_from_tab(self, _rel_ac_thr=0.05, _max_nh=5, _max_per=0.1):
        """Setup periodic Magnetic Field from Tabulated one
        :param _rel_ac_thr: relative accuracy threshold
        :param _max_nh: max. number of harmonics to create
        :param _max_per: max. period length to consider
        """

        sErMes = 'Magnetic Field is not defined'
        if(self.mag is None): raise Exception(sErMes)
        if(isinstance(self.mag, SRWLMagFldC) == False): raise Exception(sErMes)

        arHarm = []
        for i in range(_max_nh): arHarm.append(SRWLMagFldH())

        self.mag_approx = SRWLMagFldC(SRWLMagFldU(arHarm))
        srwl.UtiUndFromMagFldTab(self.mag_approx, self.mag, [_rel_ac_thr, _max_nh, _max_per])
        return self.mag_approx

    #------------------------------------------------------------------------
    def set_mag_tab(self, _fpath, _zc=0, _interp_ord=1):
        """Setup magnetic field container with tabulated magnetic field data
        :param _fpath: path to magnetic field data file
        :param _zc: center position[m]
        :param _interp_ord: order of interpolation: 1- (bi-)linear, 2- (bi-)quadratic, 3- (bi-)cubic
        """

        fpath = _fpath.strip() #MR13012017 #OC13012017

        #if(os.path.exists(_fpath) == False):
        if(os.path.exists(fpath) == False): #OC13012017
            raise Exception('No magnetic field data are supplied')

        #self.mag = srwl_uti_read_mag_fld_3d(_fpath)
        self.mag = srwl_uti_read_mag_fld_3d(fpath) #OC13012017
        self.mag.arZc[0] += _zc #?

        #print('mag was set up in srwl_bl')
        return self.mag

    #------------------------------------------------------------------------
    def set_gsn_beam(self, _x=0, _y=0, _z=0, _xp=0, _yp=0, _avgPhotEn=1, _pulseEn=1, _repRate=1, _polar=1,
                     _sigX=10e-06, _sigY=10e-06, _sigT=1e-15, _mx=0, _my=0, _presCA='c', _presFT='t'):
        """Setup Gaussian beam source
        :param _x: average horizontal coordinates of waist [m]
        :param _y: average vertical coordinates of waist [m]
        :param _z: average longitudinal coordinate of waist [m]
        :param _xp: average horizontal angle at waist [rad]
        :param _yp: average verical angle at waist [rad]
        :param _avgPhotEn: average photon energy [eV]
        :param _pulseEn: energy per pulse [J]
        :param _repRate: rep. rate [Hz]
        :param _polar: polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left
        :param _sigX: RMS beam size vs horizontal position [m] at waist (for intensity)
        :param _sigY: RMS beam size vs vertical position [m] at waist (for intensity)
        :param _sigT: RMS pulse duration [s] (for intensity)
        :param _mx: transverse Gauss-Hermite mode order in horizontal direction
        :param _my: transverse Gauss-Hermite mode order in vertical direction
        :param _presCA: treat _sigX, _sigY as sizes in [m] in coordinate representation (_presCA="c") or as angular divergences in [rad] in angular representation (_presCA="a")
        :param _presFT: treat _sigT as pulse duration in [s] in time domain/representation (_presFT="t") or as bandwidth in [eV] in frequency domain/representation (_presFT="f")
        """

        if(self.gsnBeam is not None): del self.gsnBeam

        sigX = _sigX
        sigY = _sigY
        if(_presCA == 'a'):
            convConstCA = _LightSp*_PlanckConst_eVs/(4*_Pi*_avgPhotEn)
            sigX = convConstCA/sigX
            sigY = convConstCA/sigY

        sigT = _sigT
        if(_presFT == 'f'):
            convConstFT = _PlanckConst_eVs/(4*_Pi)
            sigT = convConstFT/sigT
        
        self.gsnBeam = SRWLGsnBm(_x, _y, _z, _xp, _yp, _avgPhotEn, _pulseEn, _repRate, _polar, sigX, sigY, sigT, _mx, _my)
        return self.gsnBeam

    #------------------------------------------------------------------------
    def set_pt_src(self, _x=0, _y=0, _z=0, _flux=1, _unitFlux=1, _polar=1):
        """Setup Gaussian beam source
        :param _x: average horizontal coordinates of waist [m]
        :param _y: average vertical coordinates of waist [m]
        :param _z: average longitudinal coordinate of waist [m]
        :param _flux: spectral flux
        :param _unitFlux: spectral flux units: 1- ph/s/.1%bw, 2- W/eV
        :param _polar: polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left, 7- radial
        """

        if(hasattr(self, 'ptSrc')):
            if(self.ptSrc is not None): del self.ptSrc
        
        self.ptSrc = SRWLPtSrc(_x, _y, _z, _flux, _unitFlux, _polar)
        return self.ptSrc

    #------------------------------------------------------------------------
    def set_optics(self, _op, _v=None): #OC28042018
    #def set_optics(self, _op):
        """Setup optical element container
        :param _op: optical element container (SRWLOptC instance)
        :param _v: additional propagation-related options
        """

        #if((_op is None) or (isinstance(_op, SRWLOptC) == False)):
        #if((_op is None) or (not isinstance(_op, SRWLOptC))):
        #    raise Exception('Incorrect optics container (SRWLOptC) structure')

        if(_op is None): #OC16122018
            if(_v is not None):
                if(hasattr(_v, 'op_func')): #OC28032020
                    if(_v.op_func is not None): #OC28092022
                    #if(_v.op_func is None):
                        _op = _v.op_func(_v) #assuming _v.op_func is reference to set_optics function of virtual beamline script
                #if(_v.op_func is None):
                #    raise Exception('Optics container structure (SRWLOptC) or function for setting it up is not defined')
                #else: _op = _v.op_func(_v) #assuming _v.op_func is reference to set_optics function of virtual beamline script

        #print(_v)
        if((_v is not None) and (_op is not None)): #OC28032020
        #if(_v is not None): #OC28042018
            if(_v.op_dp != 0):
                _op.append_drift(_v.op_dp)
        
        if(self.optics is not None): del self.optics
        self.optics = _op

    #------------------------------------------------------------------------
    def set_detector(self, _x=0, _rx=0, _nx=0, _dx=0, _y=0, _ry=0, _ny=0, _dy=0, _ord=1, _fname=''):
        """Setup detector
        :param _x: horizontal center position of active area [m]
        :param _rx: horizontal size of active area [m]
        :param _nx: number of pixels in horizontal direction
        :param _dx: horizontal pixel size [m]
        :param _y: vertical center position of active area [m]
        :param _ry: vertical size of active area [m]
        :param _ny: number of pixels in vertical direction
        :param _dy: vertical pixel size [m]
        :param _ord: interpolation order (i.e. order of polynomials to be used at 2D interpolation)
        :param _fname: file name with detector spectral efficiency data
        """

        if((_rx <= 0) and (_nx <= 0) and (_ry <= 0) and (_ny <= 0)):
            raise Exception('Incorrect detector parameters')

        detSpecEff = 1
        eStartDetSpecEff = 0
        eFinDetSpecEff = 0
        #if(len(_fname) > 0):
            #Read-in detector spectral efficiency
            #detSpecEff = ...

        xHalfRangeDet = 0.5*_rx
        yHalfRangeDet = 0.5*_ry
        #self.detector = SRWLDet(
        return SRWLDet(
            _xStart = _x - xHalfRangeDet, _xFin = _x + xHalfRangeDet, _nx = _nx, 
            _yStart = _y - yHalfRangeDet, _yFin = _y + yHalfRangeDet, _ny = _ny, 
            _dx = _dx, _dy = _dy, 
            _spec_eff = detSpecEff, _eStart = eStartDetSpecEff, _eFin = eFinDetSpecEff)

    #------------------------------------------------------------------------
    def calc_mag_fld(self, _z, _rz, _nz, _x=0, _rx=0, _nx=1, _y=0, _ry=0, _ny=1, _mag_type=1, _ior=3, _fname=''):
        """Calculates / tabulates magnetic field (from other magnetic field inputs)
        :param _z: longitudinal center position [m] of the magnetic field calculation mesh
        :param _rz: range of longitudinal position of the magnetic field calculation mesh
        :param _nz: number of points vs longitudinal position of the magnetic field calculation mesh
        :param _x: horizontal center position [m] of the magnetic field calculation mesh
        :param _rx: range of horizontal position of the magnetic field calculation mesh
        :param _nx: number of points vs horizontal position of the magnetic field calculation mesh
        :param _y: vertical center position [m] of the magnetic field calculation mesh
        :param _ry: range of vertical position of the magnetic field calculation mesh
        :param _ny: number of points vs vertical position of the magnetic field calculation mesh
        :param _mag_type: "type" of magnetic field to use: 
            1- "Approximate", referenced by self.mag_approx; 
            2- "Accurate" (tabulated), referenced by self.mag; 
        :param _ior: interpolaton order of tabulated magnetic field
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: container with one 3D magnetic field structure
        """

        if(_mag_type == 1):
            if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag is None): raise Exception('Magnetic Field is not defined')
        else: raise Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2):
            magToUse = self.mag
            #print('Using tabulated magnetic field...')

        nTot = _nx*_ny*_nz
        aux = [0]*nTot
        arBx = array('d', aux)
        arBy = array('d', aux)
        arBz = array('d', aux)
        resFld3D = SRWLMagFld3D(_arBx=arBx, _arBy=arBy, _arBz=arBz, _nx=_nx, _ny=_ny, _nz=_nz, _rx=_rx, _ry=_ry, _rz=_rz, _nRep=1, _interp=_ior)
        resFld3DC = SRWLMagFldC(_arMagFld=[resFld3D], _arXc=[_x], _arYc=[_y], _arZc=[_z])
        prec = [0]

        print('Magnetic field calculation ... ', end='')
        #srwl.CalcMagnField(resFld3DC, 0, prec)
        srwl.CalcMagnField(resFld3DC, magToUse, prec)
        print('completed')

        if(len(_fname) > 0):
            print('Saving magnetic field data to a file ... ', end='')
            resFld3D.save_ascii(_fname, _x, _y, _z)
            print('completed')

        return resFld3DC

    #------------------------------------------------------------------------
    def calc_el_trj(self, _ctst, _ctfi, _np=50000, _mag_type=1, _fname=''):
        """Calculates electron trajectory
        :param _ctst: initial time (ct) for trajectory calculation 
        :param _ctfi: final time (ct) for trajectory calculation
        :param _np: number of points for trajectory calculation
        :param _mag_type: "type" of magnetic field to use: 
            1- "Approximate", referenced by self.mag_approx; 
            2- "Accurate" (tabulated), referenced by self.mag; 
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: trajectory structure
        """

        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined')

        if(_mag_type == 1):
            if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag is None): raise Exception('Magnetic Field is not defined')
        else: raise Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2):
            magToUse = self.mag
            #print('Using tabulated magnetic field...')

        partTraj = SRWLPrtTrj()
        partTraj.partInitCond = self.eBeam.partStatMom1
        partTraj.allocate(_np, True)
        partTraj.ctStart = _ctst #Start Time (ct) for the calculation
        partTraj.ctEnd = _ctfi

        arPrecPar = [1] #General Precision parameters for Trajectory calculation:
        #[0]: integration method No:
            #1- fourth-order Runge-Kutta (precision is driven by number of points)
            #2- fifth-order Runge-Kutta
        #[1],[2],[3],[4],[5]: absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] (yet to be tested!!) - to be taken into account only for R-K fifth order or higher
        #[6]: tolerance (default = 1) for R-K fifth order or higher
        #[7]: max. number of auto-steps for R-K fifth order or higher (default = 5000)

        print('Electron trajectory calculation ... ', end='')
        #print('Magnetic Field Object:', magToUse)
        srwl.CalcPartTraj(partTraj, magToUse, arPrecPar)
        print('completed')

        if(len(_fname) > 0):
            print('Saving trajectory data to a file ... ', end='')
            partTraj.save_ascii(_fname)
            print('completed')

        return partTraj

    #------------------------------------------------------------------------
    def calc_int_from_wfr(self, _wfr, _pol=6, _int_type=0, _det=None, _fname='', _pr=True): #OC06042018
        """Calculates intensity from electric field and saving it to a file
        :param _wfr: electric field wavefront (instance of SRWLWfr)
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
        :param _det: detector (instance of SRWLDet)
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :param _pr: switch specifying if printing tracing the execution should be done or not
        :return: 1D array with (C-aligned) resulting intensity data
        """

        if _pr:
            print('Extracting intensity and saving it to a file ... ', end='')
            t0 = time.time();
            
        sNumTypeInt = 'f'
        if(_int_type == 4): sNumTypeInt = 'd' #Phase?

        resMeshI = deepcopy(_wfr.mesh)

        depType = resMeshI.get_dep_type()
        if(depType < 0): raise Exception('Incorrect numbers of points in the mesh structure')
        
        arI = array(sNumTypeInt, [0]*resMeshI.ne*resMeshI.nx*resMeshI.ny)
        srwl.CalcIntFromElecField(arI, _wfr, _pol, _int_type, depType, resMeshI.eStart, resMeshI.xStart, resMeshI.yStart)

        if(_det is not None):
            resStkDet = _det.treat_int(arI, resMeshI)
            arI = resStkDet.arS
            resMeshI = resStkDet.mesh

        if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, resMeshI, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
        #if _pr: print('completed (lasted', round(time.time() - t0, 2), 's)')
        if _pr: print('completed (lasted', round(time.time() - t0, 3), 's)')

        return arI, resMeshI

    #------------------------------------------------------------------------
    #def calc_sr_se(self, _mesh, _samp_fact=-1, _meth=2, _rel_prec=0.01, _pol=6, _int_type=0, _mag_type=1, _fname=''):
    #def calc_sr_se(self, _mesh, _samp_fact=-1, _meth=2, _rel_prec=0.01, _pol=6, _int_type=0, _mag_type=1, _fname='', _det=None): #OC06122016
    def calc_sr_se(self, _mesh, _samp_fact=-1, _meth=2, _rel_prec=0.01, _pol=6, _int_type=0, _mag_type=1, _fname='', _det=None, _zi=0, _zf=0, _pr=True): #OC27122016
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
        :param _det: detector (instance of SRWLDet)
        :param _zi: initial lonngitudinal position [m] of electron trajectory for SR calculation
        :param _zf: final lonngitudinal position [m] of electron trajectory for SR calculation
        :return: resulting wavefront (instance of SRWLWfr) and 1D array with (C-aligned) intensity data
        """

        #if((_mesh is None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #depType = -1
        #if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        #elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        #elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6

        #depType = _mesh.get_dep_type() #OC06042018
        #if(depType < 0): raise Exception('Incorrect numbers of points in the mesh structure')

        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined')

        #print('In the beginning of calc_sr_se, mag_approx:', self.mag_approx)

        if(_mag_type == 1):
            if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag is None): raise Exception('Magnetic Field is not defined')
        else: raise Exception('Incorrect Magnetic Field type identificator')

        magToUse = self.mag_approx
        if(_mag_type == 2):
            magToUse = self.mag
            #print('Using tabulated magnetic field...')

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)
        wfr.partBeam = self.eBeam

        #zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
        #zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
        zStartInteg = _zi #longitudinal position to start integration (effective if < zEndInteg)
        zEndInteg = _zf #longitudinal position to finish integration (effective if > zStartInteg)
        npTraj = 50000 #Number of points for trajectory calculation 
        useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
        arPrecPar = [_meth, _rel_prec, zStartInteg, zEndInteg, npTraj, useTermin, _samp_fact]

        #print('calc_sr_se: magToUse=', magToUse)
        #print('calc_sr_se: magToUse.arMagFld[0]=', magToUse.arMagFld[0])

        if _pr:
            print('Single-electron SR calculation ... ', end='')
            t0 = time.time();

        #DEBUG
        #print('       e-beam z=', wfr.partBeam.partStatMom1.z)
        #print('Eel=', wfr.partBeam.partStatMom1.get_E(), ' GeV')
        #print('xe=', wfr.partBeam.partStatMom1.x, ' xpe=', wfr.partBeam.partStatMom1.xp, ' ye=', wfr.partBeam.partStatMom1.y, ' ype=', wfr.partBeam.partStatMom1.yp)
        #print('nx=', wfr.mesh.nx, ' xStart=', wfr.mesh.xStart, ' xFin=', wfr.mesh.xFin)
        #print('ny=', wfr.mesh.ny, ' yStart=', wfr.mesh.yStart, ' yFin=', wfr.mesh.yFin)
        #END DEBUG
        
        srwl.CalcElecFieldSR(wfr, 0, magToUse, arPrecPar) #calculate SR
        #if _pr: print('completed (lasted', round(time.time() - t0, 2), 's)')
        if _pr: print('completed (lasted', round(time.time() - t0, 3), 's)')

        arI = None
        resMeshI = None
        if(_int_type >= 0):
            arI, resMeshI = self.calc_int_from_wfr(wfr, _pol, _int_type, _det, _fname)
            
            #print('Extracting intensity and saving it to a file ... ', end='')
            #t0 = time.time();
            #sNumTypeInt = 'f'
            #if(_int_type == 4): sNumTypeInt = 'd' #Phase?

            #resMeshI = deepcopy(wfr.mesh)
            #arI = array(sNumTypeInt, [0]*resMeshI.ne*resMeshI.nx*resMeshI.ny)
            #srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, resMeshI.eStart, resMeshI.xStart, resMeshI.yStart)

            #if(_det is not None): #OC06122016
            #    resStkDet = _det.treat_int(arI, resMeshI) #OC11012017
            #    arI = resStkDet.arS
            #    resMeshI = resStkDet.mesh

            #if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, resMeshI, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
            #print('completed (lasted', round(time.time() - t0, 6), 's)')

        return wfr, arI, resMeshI #OC06122016

    #------------------------------------------------------------------------
    #def calc_rad_gsn(self, _mesh, _samp_fact=-1, _pol=6, _int_type=0, _presFT='f', _unitE=2, _fname=''):
    def calc_rad_gsn(self, _mesh, _samp_fact=-1, _pol=6, _int_type=0, _presFT='f', _unitE=2, _fname='', _det=None, _pr=True): #OC06122016
        """Calculates Gaussian beam wavefront (electric field) and intensity
        :param _mesh: mesh on which the intensity has to be calculated (SRWLRadMesh instance)
        :param _samp_fact: sampling factor for adjusting nx, ny (effective if > 0)
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
            0- "Single-Electron" / Coherent Beam Intensity; 
            1- "Multi-Electron" / Partially-Coherent Beam Intensity; 
            2- "Single-Electron" / Coherent Beam Flux; 
            3- "Multi-Electron" / Partially-Coherent Beam Flux; 
            4- "Single-Electron" / Coherent Beam Radiation Phase; 
            5- Re(E): Real part of Single-Electron / Coherent Beam Electric Field;
            6- Im(E): Imaginary part of Single-Electron / Coherent Beam Electric Field;
            7- "Single-Electron" / Coherent Beam Intensity, integrated over Time or Photon Energy (i.e. Fluence);
        :param _presFT: calculate electric field (and intensity) in time domain/representation (="t") or in frequency domain/representation (="f")
        :param _unitE: #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: 1D array with (C-aligned) resulting intensity data
        """

        #if((_mesh is None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #depType = -1
        #if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        #elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        #elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6

        #depType = _mesh.get_dep_type() #OC11102017
        #if(depType < 0): raise Exception('Incorrect numbers of points in the mesh structure')

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)

        wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        if(_presFT == "t"): wfr.presFT = 1

        wfr.unitElFld = _unitE;

        locGsnBeam = copy(self.gsnBeam) #OC16102017
        if(self.eBeam is not None):
            locGsnBeam.x += self.eBeam.partStatMom1.x
            locGsnBeam.y += self.eBeam.partStatMom1.y
            locGsnBeam.z = self.eBeam.partStatMom1.z #?
            locGsnBeam.xp += self.eBeam.partStatMom1.xp
            locGsnBeam.yp += self.eBeam.partStatMom1.yp

        #wfr.partBeam.partStatMom1.x = self.gsnBeam.x #Some information about the source in the Wavefront structure
        #wfr.partBeam.partStatMom1.y = self.gsnBeam.y
        #wfr.partBeam.partStatMom1.z = self.gsnBeam.z
        #wfr.partBeam.partStatMom1.xp = self.gsnBeam.xp
        #wfr.partBeam.partStatMom1.yp = self.gsnBeam.yp
        #OC16102017
        wfr.partBeam.partStatMom1.x = locGsnBeam.x #Some information about the source in the Wavefront structure
        wfr.partBeam.partStatMom1.y = locGsnBeam.y
        wfr.partBeam.partStatMom1.z = locGsnBeam.z
        wfr.partBeam.partStatMom1.xp = locGsnBeam.xp
        wfr.partBeam.partStatMom1.yp = locGsnBeam.yp

        if(self.eBeam is not None): #OC16102017
            #print('Debug: 2nd order stat moments of incoherent beam:')
            for i in range(len(wfr.partBeam.arStatMom2)):
                wfr.partBeam.arStatMom2[i] = self.eBeam.arStatMom2[i]
                #print(i, wfr.partBeam.arStatMom2[i])

        if _pr:
            print('Gaussian beam electric field calculation ... ', end='')
            t0 = time.time();
        
        #TEST
        #print(wfr.mesh.zStart)
        
        #srwl.CalcElecFieldGaussian(wfr, self.gsnBeam, [_samp_fact])
        srwl.CalcElecFieldGaussian(wfr, locGsnBeam, [_samp_fact]) #OC16102017
        #if _pr: print('completed (lasted', round(time.time() - t0, 2), 's)')
        if _pr: print('completed (lasted', round(time.time() - t0, 3), 's)')

        arI = None
        resMeshI = None
        if(_int_type >= 0):
            arI, resMeshI = self.calc_int_from_wfr(wfr, _pol, _int_type, _det, _fname)
            
            #print('Extracting intensity and saving it to a file ... ', end='')
            #t0 = time.time();
            #sNumTypeInt = 'f'
            #if(_int_type == 4): sNumTypeInt = 'd'

            #resMeshI = deepcopy(wfr.mesh)
            #arI = array(sNumTypeInt, [0]*resMeshI.ne*resMeshI.nx*resMeshI.ny)
            #srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, resMeshI.eStart, resMeshI.xStart, resMeshI.yStart)

            #if(_det is not None): #OC06122016
            #    resStkDet = _det.treat_int(arI, resMeshI) #OC11012017
            #    arI = resStkDet.arS
            #    resMeshI = resStkDet.mesh

            #if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, resMeshI, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
            #print('completed (lasted', round(time.time() - t0, 6), 's)')

        return wfr, arI, resMeshI #OC06122016

    #------------------------------------------------------------------------
    def calc_rad_pt_src(self, _mesh, _samp_fact=-1, _pol=6, _int_type=0, _presFT='f', _unitE=2, _fname='', _det=None, _pr=True): #OC11102017
        """Calculates Spherical Wave electric field and intensity
        :param _mesh: mesh on which the intensity has to be calculated (SRWLRadMesh instance)
        :param _samp_fact: sampling factor for adjusting nx, ny (effective if > 0)
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
            0- "Single-Electron" / Coherent Beam Intensity; 
            1- "Multi-Electron" / Partially-Coherent Beam Intensity; 
            2- "Single-Electron" / Coherent Beam Flux; 
            3- "Multi-Electron" / Partially-Coherent Beam Flux; 
            4- "Single-Electron" / Coherent Beam Radiation Phase; 
            5- Re(E): Real part of Single-Electron / Coherent Beam Electric Field;
            6- Im(E): Imaginary part of Single-Electron / Coherent Beam Electric Field;
            7- "Single-Electron" / Coherent Beam Intensity, integrated over Time or Photon Energy (i.e. Fluence);
        :param _presFT: calculate electric field (and intensity) in time domain/representation (="t") or in frequency domain/representation (="f")
        :param _unitE: #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: 1D array with (C-aligned) resulting intensity data
        """

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #depType = -1
        #if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        #elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        #elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6

        #depType = _mesh.get_dep_type()
        #if(depType < 0): raise Exception('Incorrect numbers of points in the mesh structure')

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)

        wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        if(_presFT == "t"): wfr.presFT = 1

        wfr.unitElFld = _unitE;

        locPtSrc = copy(self.ptSrc) #OC16102017
        if(self.eBeam is not None): 
            locPtSrc.x += self.eBeam.partStatMom1.x
            locPtSrc.y += self.eBeam.partStatMom1.y
            locPtSrc.z = self.eBeam.partStatMom1.z #?

        #wfr.partBeam.partStatMom1.x = self.ptSrc.x #Some information about the source in the Wavefront structure
        #wfr.partBeam.partStatMom1.y = self.ptSrc.y
        #wfr.partBeam.partStatMom1.z = self.ptSrc.z
        #wfr.partBeam.partStatMom1.xp = 0
        #wfr.partBeam.partStatMom1.yp = 0
        #OC16102017
        wfr.partBeam.partStatMom1.x = locPtSrc.x #Some information about the source in the Wavefront structure
        wfr.partBeam.partStatMom1.y = locPtSrc.y
        wfr.partBeam.partStatMom1.z = locPtSrc.z
        wfr.partBeam.partStatMom1.xp = 0.
        wfr.partBeam.partStatMom1.yp = 0.

        if(self.eBeam is not None): #OC16102017
            for i in range(len(wfr.partBeam.arStatMom2)):
                wfr.partBeam.arStatMom2[i] = self.eBeam.arStatMom2[i]

        if _pr:
            print('Spherical wave electric field calculation ... ', end='')
            t0 = time.time()
            
        #srwl.CalcElecFieldPointSrc(wfr, self.ptSrc, [_samp_fact])
        srwl.CalcElecFieldPointSrc(wfr, locPtSrc, [_samp_fact])
        #if _pr: print('completed (lasted', round(time.time() - t0, 2), 's)')
        if _pr: print('completed (lasted', round(time.time() - t0, 3), 's)')

        arI = None
        resMeshI = None
        if(_int_type >= 0):
            arI, resMeshI = self.calc_int_from_wfr(wfr, _pol, _int_type, _det, _fname)
            
            #print('Extracting intensity and saving it to a file ... ', end='')
            #t0 = time.time();
            #sNumTypeInt = 'f'
            #if(_int_type == 4): sNumTypeInt = 'd'

            #resMeshI = deepcopy(wfr.mesh)
            #arI = array(sNumTypeInt, [0]*resMeshI.ne*resMeshI.nx*resMeshI.ny)
            #srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, resMeshI.eStart, resMeshI.xStart, resMeshI.yStart)

            #if(_det is not None): #OC06122016
            #    resStkDet = _det.treat_int(arI, resMeshI) #OC11012017
            #    arI = resStkDet.arS
            #    resMeshI = resStkDet.mesh

            #if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, resMeshI, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
            #print('completed (lasted', round(time.time() - t0, 6), 's)')

        #print('wfr.mesh.eStart=', wfr.mesh.eStart, 'wfr.mesh.eFin=', wfr.mesh.eFin, 'wfr.mesh.ne=', wfr.mesh.ne)
        #print('resMeshI.eStart=', resMeshI.eStart, 'resMeshI.eFin=', resMeshI.eFin, 'resMeshI.ne=', resMeshI.ne)
        #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.nx=', wfr.mesh.nx)
        #print('resMeshI.xStart=', resMeshI.xStart, 'resMeshI.xFin=', resMeshI.xFin, 'resMeshI.nx=', resMeshI.nx)
        #print('wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin, 'wfr.mesh.ny=', wfr.mesh.ny)
        #print('resMeshI.yStart=', resMeshI.yStart, 'resMeshI.yFin=', resMeshI.yFin, 'resMeshI.ny=', resMeshI.ny)

        return wfr, arI, resMeshI

    #------------------------------------------------------------------------
    def prep_stokes_fname(self, _type, _pol=6, _fname=None): #OC10102023
        """Adjusts / modifies (aux.) Stokes file name"""

        charMultiE = 0 #Calculate intensity (flux per unit surface by default)
        if(_type == 1): charMultiE = 10 #Calculate flux

        fnStk = _fname #copy(_fname)
        
        if(((_pol >= 0) and (_pol < 6)) or (_pol == 7)): #OC10102023 (requiring 4 Stokes components)
        #if((_pol >= 0) and (_pol < 6)): #OC05102023 (requiring 4 Stokes components)
            charMultiE += 1
            
            #OC10102023 (adding "_stk" before ".dat" here)
            len_fname = 0 if _fname is None else len(_fname)
            indLastDot = fnStk.rfind('.')
            if((indLastDot >= 0) and (indLastDot < len_fname)):
                nameCore = fnStk[:indLastDot]
                sExt = fnStk[indLastDot:len_fname] #extension with '.'
                #sExt = fnStk[indLastDot+1:len_fname]
                fnStk = nameCore + '_stk' + sExt

            #if(charMultiE == 0): charMultiE = 1
            #elif(charMultiE == 10): charMultiE = 11
        return fnStk, charMultiE

    #------------------------------------------------------------------------
    #OC27122023: moved to srwlib.SRWLStokes
    # def extr_save_spec_comp(self, stk, _type=2, _pol=6, _fname=None, _fname_stk=None): #OC10102023
    #     """Extracts Flux / Intensity at different polarizations from Stokes data"""
        
    #     arI = None
    #     if(_pol < 7):
    #         arI = stk.to_int(_pol)
    #     else:
    #         arI = []
    #         for i in range(6):
    #             arI.append(stk.to_int(i))
            
    #     if(srwl_uti_proc_is_master()):
    #         sValName = 'Flux'
    #         sValUnitName = 'ph/s/.1%bw'
    #         if(_type == 2):
    #             sValName = 'Intensity'
    #             sValUnitName = 'ph/s/.1%bw/mm^2'

    #         if(_pol != 6): 
    #             if(_fname_stk is not None):
    #                 if(len(_fname_stk) > 0):
    #                     srwl_uti_save_intens_ascii(stk.arS, stk.mesh, _fname_stk, _n_stokes=4, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
                
    #         if(_fname is not None):
    #             if(len(_fname) > 0):
    #                 #sValName = 'Flux'
    #                 #sValUnitName = 'ph/s/.1%bw'
    #                 #if(_type == 2):
    #                 #    sValName = 'Intensity'
    #                 #    sValUnitName = 'ph/s/.1%bw/mm^2'
                        
    #                 if(_pol < 7): #OC10102023
    #                     srwl_uti_save_intens_ascii(arI, stk.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
    #                 else: #OC10102023 (adding "_0" .. "_5" before ".dat" here)
    #                     len_fname = len(_fname)
    #                     indLastDot = _fname.rfind('.')
    #                     nameCore = ''; sExt = ''
    #                     if((indLastDot >= 0) and (indLastDot < len_fname)):
    #                         nameCore = _fname[:indLastDot]
    #                         sExt = _fname[indLastDot:len_fname] #extension with '.'
    #                     for i in range(6):
    #                         fnPol = nameCore + '_' + repr(i) + sExt
    #                         srwl_uti_save_intens_ascii(arI[i], stk.mesh, fnPol, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
    #     return arI

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
            7- All Polarization Components (0-5)
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :return: 1D array with (C-aligned) resulting intensity data
        """

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #depType = -1
        #if((_mesh.ne >= 1) and (_mesh.nx == 1) and (_mesh.ny == 1)): depType = 0
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 1
        #elif((_mesh.ne == 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 2
        #elif((_mesh.ne == 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 3
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny == 1)): depType = 4
        #elif((_mesh.ne > 1) and (_mesh.nx == 1) and (_mesh.ny > 1)): depType = 5
        #elif((_mesh.ne > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)): depType = 6

        depType = _mesh.get_dep_type() #OC06042018
        if(depType < 0): raise Exception('Incorrect numbers of points in the mesh structure')

        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined')
        if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')

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

        fnStk, charMultiE = self.prep_stokes_fname(_type, _pol, _fname) #OC10102023

        srwl.CalcStokesUR(stk, eBeamAux, und, arPrecPar)
        #Consider treating detector here?

        return stk.extr_save_pol_comp(_type, _pol, _fname, fnStk) #OC15022024
        #return stk.extr_save_spec_comp(_type, _pol, _fname, fnStk) #OC27122023
        #return self.extr_save_spec_comp(stk, _type, _pol, _fname, fnStk) #OC10102023

        #arI = stk.to_int(_pol)
        #if(len(_fname) > 0):
        #    sValName = 'Flux'
        #    sValUnitName = 'ph/s/.1%bw'
        #    if(_type == 2):
        #        sValName = 'Intensity'
        #        sValUnitName = 'ph/s/.1%bw/mm^2'
        #    srwl_uti_save_intens_ascii(arI, stk.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
        #return arI
    
    #------------------------------------------------------------------------
    #def calc_arb_spec_me(self, _mesh, _meth=2, _rel_prec=0.01, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=10, _type=2, _mag=2, _pol=0, _rand_meth=1, _fname=None):
    #def calc_arb_spec_me(self, _mesh, _meth=2, _rel_prec=0.01, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=10, _type=2, _mag=2, _pol=0, _rand_meth=1, _fname=None, _sr_samp_fact=-1, _det=None, _me_approx=0): #OC13042018
    def calc_arb_spec_me(self, _mesh, _meth=2, _rel_prec=0.01, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=10, _type=2, _mag=2, _pol=0, _rand_meth=1, _fname=None, _sr_samp_fact=-1, _det=None, _me_approx=0, _fbk=False): #OC14082018
        """Calculates multi-electron flux of undulator radiation (within fixed aperture of per unit surface), using approximate periodic magnetic field
        :param _mesh: mesh on which the intensity has to be calculated (SRWLRadMesh instance)
        :param _meth: SR Electric Field calculation method to be used (0- "manual", 1- "auto-undulator", 2- "auto-wiggler")
        :param _rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, smaller the more accurate)
        :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
        :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
        :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
        :param _type: calculate flux (=1) or flux per unit surface (=2)
        :param _mag: magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
        :param _rand_meth: method for generation of pseudo-random numbers for e-beam phase-space integration:
            1- standard pseudo-random number generator
            2- Halton sequences
            3- LPtau sequences (to be implemented)
        :param _fname: name of file to save the resulting data to (for the moment, in ASCII format)
        :param _det: detector structure ensuring a given final mesh on which the calculated intensity (or other characteristic) will be interpolated
        :param _me_approx: multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance
        :return: 1D array with (C-aligned) resulting intensity data
        """

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined')

        mag2use = self.mag
        if(_mag == 1):
            if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
            mag2use = self.mag_approx
        else:
            if(self.mag is None): raise Exception('Accurate Magnetic Field is not defined')

        fnStk, charMultiE = self.prep_stokes_fname(_type, _pol, _fname) #OC10102023
        
        #charMultiE = 0 #Calculate intensity (flux per unit surface by default)
        #if(_type == 1): charMultiE = 10 #Calculate flux

        # fnStk = copy(_fname) #OC10102023
        
        # if(((_pol >= 0) and (_pol < 6)) or (_pol == 7)): #OC10102023 (requiring 4 Stokes components)
        # #if((_pol >= 0) and (_pol < 6)): #OC05102023 (requiring 4 Stokes components)
        #     charMultiE += 1
            
        #     #OC10102023 (adding "_stk" before ".dat" here)
        #     len_fname = len(_fname)
        #     indLastDot = fnStk.rfind('.')
        #     if((indLastDot >= 0) and (indLastDot < len_fname)):
        #         nameCore = fnStk[:indLastDot]
        #         sExt = fnStk[indLastDot:len_fname] #extension with '.'
        #         #sExt = fnStk[indLastDot+1:len_fname]
        #         fnStk = nameCore + '_stk' + sExt

        #     #if(charMultiE == 0): charMultiE = 1
        #     #elif(charMultiE == 10): charMultiE = 11
        
        #print(_fname)
        stk = srwl_wfr_emit_prop_multi_e(
            _e_beam = self.eBeam, _mag = mag2use, _mesh = _mesh,
            _sr_meth = _meth, _sr_rel_prec = _rel_prec,
            _n_part_tot = _n_part_tot, _n_part_avg_proc = _n_part_avg_proc, _n_save_per = _n_save_per, _rand_meth = _rand_meth,
            #_file_path = _fname, _char = charMultiE)
            #_file_path = _fname, _sr_samp_fact = _sr_samp_fact, _char = charMultiE, _det = _det, _me_approx = _me_approx, #) #OC14042018
            _file_path = fnStk, _sr_samp_fact = _sr_samp_fact, _char = charMultiE, _det = _det, _me_approx = _me_approx, #) #OC14042018
            _file_bkp = True if(_fbk == True) else False) #OC14082018
            
        #Consider treating detector here?
        
        if(srwl_uti_proc_is_master()): #OC22122023
            return stk.extr_save_pol_comp(_type, _pol, _fname) #OC15022024
            #return stk.extr_save_spec_comp(_type, _pol, _fname) #OC27122023
            #return self.extr_save_spec_comp(stk, _type, _pol, _fname)
        else:
            return None 

        ##arI = None
        ##if(stk is not None):
        #return self.extr_save_spec_comp(stk, _type, _pol, _fname)
        ##return arI

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

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #if(self.eBeam is None): Exception('Electron Beam structure is not defined')
        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined') #AH03022019 (identified raise missing in many places)

        if(_mag_type == 1):
            if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag is None): raise Exception('Magnetic Field is not defined')
        else: raise Exception('Incorrect Magnetic Field type identificator')

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
        return stkP.arS #, arSx, arSy

    #------------------------------------------------------------------------
    def calc_und_oper_tab(self, _mesh, _pol=0, _hi=1, _hf=1, _meas_or_calc='m', _zc=0, _fname=''):
        """Calculate undulator "operation table", i.e. dependence of gap (and phase) on photon energy (for a given polarization)
        :param _mesh: mesh (grid) for which the operation table has to be calculated (SRWLRadMesh instance)

        """
        #['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],

        #print('Calculating undulator operation table')

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #print('_mesh.xStart=', _mesh.xStart, '_mesh.xFin=', _mesh.xFin)

        if(self.eBeam is None): raise Exception('Electron Beam structure is not defined')

        fPathSum = ''
        if(_meas_or_calc == 'm'):
            if(hasattr(self, 'dir_magn_meas') and hasattr(self, 'fn_magn_meas_sum')):
                fPathSum = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
            else: raise Exception('No magnetic measurements data are supplied')
        elif(_meas_or_calc == 'c'):
            raise Exception('No magnetic calculation data are supplied')

        f = open(fPathSum, 'r')
        lines = f.readlines() #read-in all lines
        nRows = len(lines)

        strSep = '\t'
        arGaps = []; arPhases = []; arMagFld3D = []
        arXc = []; arYc = []; arZc = []
        arCoefBx = []; arCoefBy = []

        arHarmUR = []; arEnMaxIntSE = []; arEnMaxFluxEst = []; arMaxFlux = [];  arPowTot = []; arPowInAp = []

        phaseIsVar = False
        phasePrev = None
        #print('Setting up tabulated magnetic field')

        phModeReq = 'p1'
        if((_pol == 2) or (_pol == 3)): #2- Linear 45 degrees, 3- Linear 135 degrees
            phModeReq = 'p2'

        powNumX = 101; powNumY = 101
        #powNumX = 11; powNumY = 11
        nxPartIntegPowDens = 101; nyPartIntegPowDens = 101

        relAcThrConv2Per = 0.05
        maxNumHarmConv2Per = 7
        maxPerConv2Per = 0.2 #to steer?

        curMeshF = deepcopy(_mesh)
        curMeshF.nx = 1; curMeshF.ny = 1

        curMeshI = deepcopy(_mesh)

        for i in range(nRows):
            curLine = lines[i]
            curLineParts = curLine.split(strSep)
            curLenLineParts = len(curLineParts)
            
            curGap = None
            curPhase = 0
            if(curLenLineParts >= 4):
                curPhaseMode = curLineParts[1]
                if(curPhaseMode != phModeReq): continue

                curGap = float(curLineParts[0])
                arGaps.append(curGap)
                print('Gap:', curGap, 'mm')

                curPhase = float(curLineParts[2])
                if((phasePrev is not None) and (curPhase != phasePrev)): phaseIsVar = True
                arPhases.append(curPhase)
                phasePrev = curPhase

                #curFileName = curLineParts[3]
                curFileName = curLineParts[3].strip() #MR13012017
                print('Magnetic Field Data File:', curFileName)

                curFldCnt = None
                curZc = 0.
                curCoefBx = 1.
                curCoefBy = 1.
                
                if(len(curFileName) > 0):
                    curFilePath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, curFileName)
                    curFldCnt = srwl_uti_read_mag_fld_3d(curFilePath, '#')
                    #arMagFld3D.append(curFldCnt.arMagFld[0])
                    #arXc.append(curFldCnt.arXc[0])
                    #arYc.append(curFldCnt.arYc[0])
                    #arZc.append(curFldCnt.arZc[0] + _zc)
                    curFldCnt.arZc[0] += _zc
                    
                if(curLenLineParts >= 6):
                    #arCoefBx.append(float(curLineParts[4]))
                    #arCoefBy.append(float(curLineParts[5]))
                    curCoefBx = float(curLineParts[4])
                    curCoefBy = float(curLineParts[5])

                    if(((curCoefBx != 1.) or (curCoefBy != 1.)) and (curFldCnt is not None)):
                        curFld3D = curFldCnt.arMagFld[0]
                        iif = 0
                        for iz in range(curFld3D.nz):
                            for iy in range(curFld3D.ny):
                                for ix in range(curFld3D.nx):
                                    curFld3D.arBx[iif] *= curCoefBx
                                    curFld3D.arBy[iif] *= curCoefBy
                                    iif += 1

                if(curFldCnt is not None):
                    
                    #Convert field to periodic
                    arHarm = []
                    for iih in range(7): arHarm.append(SRWLMagFldH())
                    undMagApprox = SRWLMagFldC(SRWLMagFldU(arHarm))

                    srwl.UtiUndFromMagFldTab(undMagApprox, curFldCnt, [relAcThrConv2Per, maxNumHarmConv2Per, maxPerConv2Per])

                    undApprox = undMagApprox.arMagFld[0]
                    elEnGeV = self.eBeam.partStatMom1.get_E('GeV')
                    e1Approx = undApprox.get_E1(elEnGeV, 'eV')
                    
                    #print('Undulator period:', undApprox.per, 'm')
                    #print('Electron energy:', elEnGeV, 'GeV')
                    #print('Undulator e1 =', e1Approx, 'eV')

                    #Define spectral range and mesh for single-e & multi-e calculation
                    #curMesh = deepcopy(_mesh)
                    curMeshF.eStart = e1Approx*(_hi - 0.5)
                    curMeshF.eFin = e1Approx*(_hf + 0.5)
                    curMeshF.ne = 1000*(_hf - _hi + 1)
                    eStep = (curMeshF.eFin - curMeshF.eStart)/(curMeshF.ne - 1)

                    #Calculate multi-e spectrum
                    stkSpF = SRWLStokes() #for spectral flux vs photon energy
                    stkSpF.allocate(curMeshF.ne, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
                    stkSpF.mesh = curMeshF
                    longPrecStkF = 1.5 #longitudinal integration precision parameter
                    azPrecStkF = 1.5 #azimuthal integration precision parameter
                    srwl.CalcStokesUR(stkSpF, self.eBeam, undMagApprox.arMagFld[0], [1, _hf+5, longPrecStkF, azPrecStkF, 1])

                    #Calculate single-e spectrum
                    wfrSp = SRWLWfr()
                    curMeshI.eStart = curMeshF.eStart
                    curMeshI.eFin = curMeshF.eFin
                    curMeshI.ne = curMeshF.ne
                    wfrSp.allocate(curMeshI.ne, 1, 1) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions (may be modified by the library!)
                    
                    wfrSp.mesh = curMeshI
                    wfrSp.partBeam = self.eBeam

                    methSR = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
                    relPrecSR = 0.01 #relative precision
                    srwl.CalcElecFieldSR(wfrSp, 0, curFldCnt, [methSR, relPrecSR, 0, 0, 20000, 1, -1])
                    
                    arSpecI = array('f', [0]*wfrSp.mesh.ne) #"flat" array to take 2D intensity data
                    srwl.CalcIntFromElecField(arSpecI, wfrSp, _pol, 0, 0, wfrSp.mesh.eStart, 0, 0)

                    #Calculate multi-e power-density distribution
                    stkP = SRWLStokes() #for power density
                    stkP.allocate(1, powNumX, powNumY) #numbers of points vs horizontal and vertical positions (photon energy is not taken into account)
                    stkP.mesh.zStart = curMeshF.zStart #longitudinal position [m] at which power density has to be calculated
                    powHalfRangeX = 2.*curMeshF.zStart*undApprox.get_K()/self.eBeam.partStatMom1.gamma
                    powHalfRangeY = 2.*curMeshF.zStart/self.eBeam.partStatMom1.gamma
                    stkP.mesh.xStart = -powHalfRangeX #initial horizontal position [m]
                    stkP.mesh.xFin = powHalfRangeX #final horizontal position [m]
                    stkP.mesh.yStart = -powHalfRangeY #initial vertical position [m]
                    stkP.mesh.yFin = powHalfRangeY #final vertical position [m]

                    #print('powHalfRangeX=', powHalfRangeX, 'powHalfRangeY=', powHalfRangeY)
                    #print('curMeshF.xStart=', curMeshF.xStart, 'curMeshF.xFin=', curMeshF.xFin)
                    #print('curMeshF.yStart=', curMeshF.yStart, 'curMeshF.yFin=', curMeshF.yFin)
                    #print('curMeshI.xStart=', curMeshI.xStart, 'curMeshI.xFin=', curMeshI.xFin)
                    #print('curMeshI.yStart=', curMeshI.yStart, 'curMeshI.yFin=', curMeshI.yFin)

                    precFactPow = 1.5 #precision factor
                    methPow = 1 #power density computation method (1- "near field", 2- "far field")
                    srwl.CalcPowDenSR(stkP, self.eBeam, 0, curFldCnt, [precFactPow, methPow, 0, 0, 20000])
                    powTot = uti_math.integ_ar_2d(stkP.arS, 1, [stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx], [stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny])*1.e+06

                    powInAp = uti_math.integ_ar_2d(stkP.arS, 1, [stkP.mesh.xStart, stkP.mesh.xFin, stkP.mesh.nx], [stkP.mesh.yStart, stkP.mesh.yFin, stkP.mesh.ny],
                                                   [curMeshF.xStart, curMeshF.xFin, nxPartIntegPowDens], [curMeshF.yStart, curMeshF.yFin, nyPartIntegPowDens])*1.e+06

                    print('Power ~total:', powTot, 'W')
                    print('Power within work aperture:', powInAp, 'W')
                    
                    #Determine energy shifts of different harmonics and required harmonic positions
                    #arEnShift = []
                    for iHarmUR in range(_hi, _hf + 1, 1):
                        curEnStart = (iHarmUR - 0.5)*e1Approx
                        ieStartSearch = 0
                        if((curEnStart > curMeshF.eStart) and (curEnStart <= curMeshF.eFin)):
                            ieStartSearch = int((curEnStart - curMeshF.eStart)/eStep + 1.e-10)
                        if(ieStartSearch < 0): ieStartSearch = 0

                        curEnEnd = (iHarmUR + 0.5)*e1Approx
                        ieEndSearch = 0
                        if((curEnEnd > curMeshF.eStart) and (curEnEnd <= curMeshF.eFin)):
                            ieEndSearch = int((curEnEnd - curMeshF.eStart)/eStep + 1.e-10)

                        if(ieEndSearch >= curMeshF.ne): ieEndSearch = curMeshF.ne - 1
                        if(ieEndSearch < ieStartSearch): ieEndSearch = ieStartSearch

                        curMaxFlux, iCurMaxFlux = uti_math.find_ar_max(stkSpF.arS, ieStartSearch, ieEndSearch)
                        enMaxFlux = curMeshF.eStart + eStep*iCurMaxFlux
                        enShift = enMaxFlux - e1Approx*iHarmUR
                        #arEnShift.append(enShift)

                        curMaxIntSE, iCurMaxIntSE = uti_math.find_ar_max(arSpecI, ieStartSearch, ieEndSearch)
                        enMaxIntSE = curMeshI.eStart + eStep*iCurMaxIntSE
                        enMaxFluxEst = enMaxIntSE + enShift
                        
                        print('Harmonic Number:', iHarmUR)
                        print('Approx. Photon Energy for Max. Flux:', enMaxFlux, 'eV')
                        print('Estimated Photon Energy for Max. Flux:', enMaxFluxEst, 'eV')
                        print('Photon Energy (Red) Shift:', enShift, 'eV')

                        #To improve:
                        #Treat different polarizations
                        #Determine polarization rate of required polarization

                        #arResForHarm.append([iHarmUR, curGap, curMaxIntSE, enMaxFluxEst, curMaxFlux, powInAp])
                        arHarmUR.append(iHarmUR)
                        arEnMaxIntSE.append(round(enMaxIntSE, 5))
                        arEnMaxFluxEst.append(round(enMaxFluxEst, 5))
                        arMaxFlux.append(round(curMaxFlux, 4))
                        arPowTot.append(round(powTot, 4))
                        arPowInAp.append(round(powInAp, 4))

        f.close()
        arResForHarm = [arHarmUR, arGaps, arPhases, arEnMaxIntSE, arEnMaxFluxEst, arMaxFlux, arPowTot, arPowInAp]

        if(len(_fname) > 0):
            strHeader = '#harm, gap, phase, en_res, en_max, flux, pow_tot, pow_in_ap'
            srwl_uti_write_data_cols(_fname, arResForHarm, '\t', strHeader)

        return arResForHarm

    #------------------------------------------------------------------------
    #def calc_wfr_prop(self, _wfr, _pres_ang=0, _pol=6, _int_type=0, _dep_type=3, _fname=''):
    #def calc_wfr_prop(self, _wfr, _pres_ang=0, _pol=6, _int_type=0, _dep_type=3, _fname='', _det=None): #OC06122016
    def calc_wfr_prop(self, _wfr, _pres_ang=0, _pol=6, _int_type=0, _dep_type=3, _fname='', _det=None, _rad_view=None): #OC08022021
        """Calculates single-electron (/ fully coherent) wavefront propagation
        :param _wfr: wavefront (instance of SRWLWfr) to be propagated (and modified in place!)
        :param _pres_ang: switch specifying whether the result of the propagation should be shown in angular presentation (1) or not (0)
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
        :param _det: detector object defining mesh and possibly other characteristics of radiation to be extracted after the propagation
        :param _rad_view: list of indexes (1-based) of optical elements in self.optics container after which the intensity distributions or other radiation characteristics shoudl be extracted
        :return: 1D array with (C-aligned) resulting intensity data; it also modified _wfr in place
        """

        if((hasattr(self, 'optics') == False) or (isinstance(self.optics, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')
        if(isinstance(_wfr, SRWLWfr) == False):
            raise Exception('Incorrect wavefront (SRWLWfr) structure')

        from datetime import datetime #OCTEST

        print('Propagation ... ', end='')
        t0 = time.time();
        #dt = datetime.now() #OCTEST
        #t0 = dt.now(); #OCTEST
        #print(t0) #OCTEST

        radView = None #OC08022021
        if(_rad_view is not None):
            nInd = len(_rad_view)
            if(nInd > 0):
                radView = []
                for i in range(nInd):
                    intType = 0 if(_int_type < 0) else _int_type
                    radView.append([_rad_view[i], _pol, intType, _dep_type, 0, SRWLRadMesh()])
            
        srwl.PropagElecField(_wfr, self.optics, radView) #OC08022021
        #srwl.PropagElecField(_wfr, self.optics)

        #dt = datetime.now() #OCTEST
        #t = dt.now(); #OCTEST
        #print(t) #OCTEST
        
        #print('completed (lasted', t - t0, 's)')#OCTEST
        #print('completed (lasted', round(time.time() - t0, 2), 's)')
        print('completed (lasted', round(time.time() - t0, 3), 's)')

        #print('_wfr.Rx=',  _wfr.Rx, '   _wfr.Ry=',  _wfr.Ry)

        if(_pres_ang != 0): srwl.SetRepresElecField(_wfr, 'a')

        sValUnitName = 'ph/s/.1%bw/mm^2' #OC08022021 (moved up) #consider allowing for other units (for FEL applications)

        if(radView is not None): #OC08022021
            print('Saving Intermediate Propagation Results ... ', end='')
            t0 = time.time();

            nInd = len(radView)
            if(nInd > 0):
                intNameCore = 'intermed_rad_data'
                sExt = 'dat'
                len_fname = len(_fname)
                if((_fname is not None) and (len_fname > 0)):
                    indLastDot = _fname.rfind('.')
                    if((indLastDot >= 0) and (indLastDot < len_fname)):
                        intNameCore = _fname[:indLastDot]
                        sExt = _fname[indLastDot+1:len_fname]

                for i in range(nInd):
                    curRadData = radView[i]
                    if((len(curRadData) > 6) and (curRadData[6] is not None)):
                        curIntFileName = intNameCore + '_intermed_' + repr(curRadData[0]) + '.' + sExt
                        #DEBUG
                        #print(curIntFileName)
                        #END DEBUG
                        srwl_uti_save_intens_ascii(curRadData[6], curRadData[5], curIntFileName, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', sValUnitName])

            print('completed (lasted', round(time.time() - t0, 3), 's)')
 
        arI = None
        if(_int_type >= 0): 
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            #resMeshI = _wfr.mesh
            resMeshI = deepcopy(_wfr.mesh)
            arI = array(sNumTypeInt, [0]*resMeshI.ne*resMeshI.nx*resMeshI.ny)
            srwl.CalcIntFromElecField(arI, _wfr, _pol, _int_type, _dep_type, resMeshI.eStart, resMeshI.xStart, resMeshI.yStart)
            
            if(_det is not None): #OC06122016
                #resStkDet = _det.treat_int(arI, resMeshI, _ord_interp=1)
                resStkDet = _det.treat_int(arI, resMeshI) #OC11012017
                arI = resStkDet.arS
                resMeshI = resStkDet.mesh

            if(len(_fname) > 0):
                #sValUnitName = 'ph/s/.1%bw/mm^2' #consider allowing for other units (for FEL applications)

                print('Saving Propagation Results ... ', end='')
                t0 = time.time();
                srwl_uti_save_intens_ascii(arI, resMeshI, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', sValUnitName])
                #print('completed (lasted', round(time.time() - t0, 2), 's)')
                print('completed (lasted', round(time.time() - t0, 3), 's)')

        #return arI
        #return arI, resMeshI #OC06122016
        return arI, resMeshI, radView #OC08022021

    #------------------------------------------------------------------------
    #def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None):
    #def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _in_wr=0., _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None):
    #def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _in_wr=0., _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None): #OC06122016
    #def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _in_wr=0., _in_wre=0., _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0): #OC05042017
    #def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _in_wr=0., _in_wre=0., _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False): #OC14082018
    def calc_wfr_emit_prop_me(self, _mesh, _sr_samp_fact=1, _sr_meth=2, _sr_rel_prec=0.01, _in_wr=0., _in_wre=0., _mag_type=1, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=50,
                              _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii', _no_opt=False, _nmm=1, _ncm=100, _cm_wfr=None, _ms=0, _pol=None): #OC27122023
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii', _no_opt=False, _nmm=1, _ncm=100, _cm_wfr=None, _ms=0): #OC22112022
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii', _no_opt=False, _nmm=1, _ncm=100, _cm_wfr=None): #OC02072021
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii', _nmm=1, _ncm=1000): #OC27062021
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii', _nmm=1): #OC16042021
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False, _fform='ascii'): #OC25022021
                              #_pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1, _fname=None, _det=None, _me_approx=0, _fbk=False, _op_rnd=False): #OC24042020
        """Calculates multi-electron (/ partially coherent) SR emission and wavefront propagation
        :param _mesh: mesh (grid) on which the initial wavefront has to be calculated (SRWLRadMesh instance)
        :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
        :param _sr_meth: SR Electric Field calculation method to be used (0- "manual", 1- "auto-undulator", 2- "auto-wiggler")
        :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, smaller the more accurate)
        :param _in_wr: initial wavefront radius [m] to assume at wavefront propagation (is taken into account if != 0)
        :param _in_wre: initial wavefront radius error [m] to assume at wavefront propagation (is taken into account if != 0)
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
        :param _det: detector structure ensuring a given final mesh on which the calculated intensity (or other characteristic) will be interpolated
        :param _me_approx: multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance
        :param _fbk: switch specifying if backup intermadiate intensity files should be produced in the course of calculation or not
        :param _op_rnd: switch specifying whether some optical element parameters, as specified in elements self.optics, should be randomized in the course of the calculation or not
        :param _fform: switch specifying file format of some input/output data files (not fully supported yet)
        :param _nmm: number of MPI "masters" to be used at 4D cross-spectral density (CSD) calculation
        :param _ncm: number of coherent modes to be produced by coherent mode decomposition (CMD) of 4D CSD
        :param _cm_wfr: list of wavefronts (objects of SRWLWfr type) describing coherent modes to be used in partially-coherent radiation propagation calculations; if this list is supplied, it is used instead of any other input radiaiton data
        :param _ms: index of the first wavefront / coherent mode to start computation
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
            7- All Polarization Components (0-5)
        :return: data structure(s) of resulting intensity
        """

        #if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
        if((_mesh is None) or (not isinstance(_mesh, SRWLRadMesh))):
            raise Exception('Incorrect SRWLRadMesh structure')

        #if((hasattr(self, 'eBeam') == False) or (isinstance(self.eBeam, SRWLPartBeam) == False)):
        if((not hasattr(self, 'eBeam')) or (not isinstance(self.eBeam, SRWLPartBeam))):
            raise Exception('Incorrect electron beam (SRWLPartBeam) structure')

        #OC25022021 (commented-out to enable calculaiton of P-C SR without propagation)
        #if((not hasattr(self, 'optics')) or (not isinstance(self.optics, SRWLOptC))):
        #    raise Exception('Incorrect optics container (SRWLOptC) structure')

        #OC26042020 (commented-out to make Gaussian Schell or Zernike source work)
        #if(_mag_type == 1):
        #    if(self.mag_approx is None): raise Exception('Approximate Magnetic Field is not defined')
        #elif(_mag_type == 2):
        #    if(self.mag is None): raise Exception('Magnetic Field is not defined')
        #else: raise Exception('Incorrect Magnetic Field type identificator')

        #print(self.mag_approx) #DEBUG
        #print(self.mag) #DEBUG

        magToUse = self.mag_approx
        if(_mag_type == 2): magToUse = self.mag

        #if((magToUse is None) and (self.gsnBeam is not None)): magToUse = self.gsnBeam #OC15092017 (because _mag is used for SRWLGsnBm when doing partially-coherent simulations in the scope of Gaussian-Schell model)
        if(magToUse is None):
            if(self.gsnBeam is not None): magToUse = self.gsnBeam #OC15092017 (because _mag is used for SRWLGsnBm when doing partially-coherent simulations in the scope of Gaussian-Schell model)
            elif(self.ptSrc is not None): magToUse = self.ptSrc #OC16102017 (because _mag is used for SRWLPtSrc when doing partially-coherent simulations in the scope of Van-Cittert / Zernike model)

        optBL = self.optics if(not _no_opt) else None #OC03072021

        #DEBUG
        #print('_no_opt=', _no_opt)
        #sys.stdout.flush()
        #END DEBUG

        if(_cm_wfr is not None): magToUse = _cm_wfr #OC02072021 (note it overrides other cases - should it be like this?)

        res = srwl_wfr_emit_prop_multi_e( #OC27122023
        #return srwl_wfr_emit_prop_multi_e(
            _e_beam = self.eBeam, _mag = magToUse, _mesh = _mesh, _sr_samp_fact = _sr_samp_fact,
            #_sr_meth = _sr_meth, _sr_rel_prec = _sr_rel_prec,
            #_sr_meth = _sr_meth, _sr_rel_prec = _sr_rel_prec, _w_wr = _in_wr, #OC26032016
            #_sr_meth = _sr_meth, _sr_rel_prec = _sr_rel_prec, _wr = _in_wr, #OC07092016
            _sr_meth = _sr_meth, _sr_rel_prec = _sr_rel_prec, _wr = _in_wr, _wre = _in_wre, #OC05012017
            _n_part_tot = _n_part_tot, _n_part_avg_proc = _n_part_avg_proc, _n_save_per = _n_save_per,
            _file_path = _fname,
            _opt_bl = optBL, #OC03072021
            #_opt_bl = self.optics,
            _pres_ang = _pres_ang, _char = _char, _x0 = _x0, _y0 = _y0,
            #_e_ph_integ = _e_ph_integ, _rand_meth = _rand_meth)
            #_e_ph_integ = _e_ph_integ, _rand_meth = _rand_meth, _det = _det) #OC06122016
            #_e_ph_integ = _e_ph_integ, _rand_meth = _rand_meth, _det = _det, _me_approx = _multi_e_approx) #OC05042017
            _e_ph_integ = _e_ph_integ, _rand_meth = _rand_meth, _det = _det, _me_approx = _me_approx, 
            #_file_bkp = True if(_fbk == True) else False) #OC14082018
            #_file_bkp = _fbk, _rand_opt = _op_rnd) #OC24042020
            #_file_bkp = _fbk, _rand_opt = _op_rnd, _file_form = _fform) #OC25022021
            #_file_bkp = _fbk, _rand_opt = _op_rnd, _file_form = _fform, _n_mpi=_nmm) #OC16042021
            #_file_bkp = _fbk, _rand_opt = _op_rnd, _file_form = _fform, _n_mpi=_nmm, _n_cm=_ncm) #OC27062021
            _file_bkp = _fbk, _rand_opt = _op_rnd, _file_form = _fform, _n_mpi = _nmm, _n_cm = _ncm, _ms = _ms) #OC22112022
        
        if(_pol is not None):
            if(res is not None):
                if(isinstance(res, SRWLStokes)):
                    res.extr_save_pol_comp(_type=2, _pol=_pol, _fname=_fname, _fname_stk=None) #Save intensity at a given polarization
        return res

    #------------------------------------------------------------------------
    def cost_func_aux_int_distr(self, _x, *_aux):

        x_lim = _aux[0]
        if(x_lim is not None):
            if(not self.uti_math_opt.check_limits(_x, x_lim)): return 1.e+23 #to tune

        v = _aux[1] #all params for optical simulation
        rStart = v.om_aux_rStart
        rFin = v.om_aux_rFin
        irStart = v.om_aux_irStart
        irFinP1 = v.om_aux_irFinP1
        arI1d = v.om_aux_arI1d

        rStep = (self.meshIntFit.xFin - self.meshIntFit.xStart)/(self.meshIntFit.nx - 1)
        nr = irFinP1 - irStart

        scaleMult = _x[0]

        useLogOnlyForF = False
        useLogForFandArg = False
        rMult = 1; 

        if(v.om_fs[1] == 1): useLogOnlyForF = True #do fitting in normal scale or apply log to function only
        elif(v.om_fs[1] == 2): #do fitting in log scale vs argument in log scale
            useLogForFandArg = True
            rMult = 10.**((log10(rFin) - log10(rStart))/(nr - 1))
        
        useNormalScale = not (useLogOnlyForF or useLogForFandArg)
        sumDifE2 = 0
        sumMeas = 0; sumCalc = 0
        r = rStart
        for ir in range(irStart, irFinP1):
            if useNormalScale or useLogOnlyForF:
                fMeas = self.arIntFit[ir]
                fCalc = scaleMult*arI1d[ir - irStart]
                if(useLogOnlyForF):
                    fMeas = log10(fMeas)
                    fCalc = log10(fCalc)
            elif useLogForFandArg:
                fMeas = log10(uti_math.interp_1d(r, self.meshIntFit.xStart, rStep, self.meshIntFit.nx, self.arIntFit, _ord=v.om_fp[4])) #interpolate using same order as at azimuthal avaraging
                fCalc = log10(scaleMult*uti_math.interp_1d(r, rStart, rStep, nr, arI1d, _ord=v.om_fp[4]))
                r *= rMult
            dF = fMeas - fCalc
            sumDifE2 += dF*dF
            sumMeas += fMeas
            sumCalc += fCalc

        v.om_aux_avg_fMeas = sumMeas/nr
        v.om_aux_avg_fCalc = sumCalc/nr
        #_aux[1] = v #not necessary?

        res = sqrt(sumDifE2/nr)
        #OCTEST
        #print('cost_func_aux_int_distr: res=', res)
        
        return res #maybe normalize it by something / apply weight?

    #------------------------------------------------------------------------
    def cost_func(self, _x, *_aux):
        """Standard Cost Function used for optimization of optical elements and possibly sources"""

        #OC28092022
        procIsMaster = True
        sizeMPI = 0
        rank = 0
        if(self.comMPI is not None):
            rank = self.comMPI.Get_rank()
            if(rank > 0): procIsMaster = False
            sizeMPI = self.comMPI.Get_size()

        if(sizeMPI > 1): #Master has to send / broadcast _x to all workers to make sure that they will be doing calculations for the same source and optics params
            self.comMPI.Barrier() #synchronizing all processes in case if MPI is used (e.g. for partially-coherent calculations)
            arX = array('d', _x)
            if procIsMaster:
                for iRank in range(sizeMPI - 1):
                    #dst = iRank + 1 + rankMaster
                    dst = iRank + 1
                    self.comMPI.Send([arX, self.MPI.DOUBLE], dest=dst)

                #self.comMPI.Bcast([arX, self.MPI.DOUBLE], root=self.MPI.ROOT)
            else:
                self.comMPI.Recv([arX, self.MPI.DOUBLE], source=self.MPI.ANY_SOURCE)
                for i in range(len(arX)): _x[i] = arX[i]

        x_lim = _aux[0]
        #print(x_lim)
        
        v = _aux[1] #all params for optical simulation
        cw = v.om_cw #optimization weights for different criterions
        ct = v.om_ct #target values of different criterions
        cn = v.om_cn #nominal values of different criterions
        #[0]- horizontal spot size,
        #[1]- vertucal spot size,
        #[2]- horizontal angular divergence,
        #[3]- vertical angular divergence,
        #[4]- horizontal coherence length,
        #[5]- vertical coherence length,
        #[6]- peak intensity,
        #[7]- given intensity distribution
        
        nmVars = v.om_pn
        nVars = len(_x)
        if(nVars != len(nmVars)): #OC08122022
            raise Exception("Inconsistent numbers of optimization parameter values and their listed names")
    
        if(v.om_mt > 0): #OC08122022 (don't check limits for the "param. scan" method) 
            if(nVars != len(x_lim)):
                raise Exception("Inconsistent numbers of optimization parameter values and their listed limits")
            if(not self.uti_math_opt.check_limits(_x, x_lim)): return 1.e+23 #to tune
            
        #if((nVars != len(nmVars)) or (nVars != len(x_lim))):
        #    raise Exception("Inconsistent numbers of optimization parameters / names / limits")
        #if(not self.uti_math_opt.check_limits(_x, x_lim)): return 1.e+23 #to tune
        
        for i in range(nVars): 
            #curNameVar = nmVars[i]
            #curVar = getattr(v, curNameVar)
            #print(curNameVar, '= ', curVar)
            setattr(v, nmVars[i], _x[i]) #setting instant values of optimization variables

        #DEBUG
        #print('     rank=', rank, 'starting calc_all()')
        #sys.stdout.flush()
        #END DEBUG

        self.calc_all(v) #run forward-simulation

        infIsReq = v.om_pr or v.om_fl #OC28092022
        s2pr = '  ' #OC28092022
        frm = '{:04.6g}' #OC28092022
        def critInf(_parDescr, _parName, _subCost): #OC28092022
            return (' '+_parDescr+'='+frm).format(_parName) + (' (sub-cost:'+frm+')').format(_subCost)

        cost = 0
        
        #Calculate cost based of weights, target and nominal values of radiation characteristics
        if(v.si or v.ws or v.wg): #fully-coherent / single-electron calculations

            #To uncomment when / if MPI will be used for fully-coherent (e.g. time-dependent) calculations:
            #if(self.comMPI is not None): self.comMPI.Barrier() #synchronizing all processes in case if MPI is used (e.g. for partially-coherent calculations)

            #s2pr = '  ' #OC28092022 (moved up)
            #frm = '{:04.6g}'
            #def critInf(_parDescr, _parName, _subCost):
            #    return (' '+_parDescr+'='+frm).format(_parName) + (' (sub-cost:'+frm+')').format(_subCost)
                
            wfr = v.w_res            
            #infIsReq = v.om_pr or v.om_fl #OC28092022 (moved up)

            if((cw['xFWHM'] > 0) or (cw['yFWHM'] > 0) or (cw['iM'] > 0) or (cw['xFWFM'] > 0) or (cw['yFWFM'] > 0)): #horizontal, vertical sizes, and peak intensity
                
                meshI = wfr.mesh
                arI = array('f', [0]*meshI.nx*meshI.ny) #"flat" 2D array to take intensity data
                srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, meshI.eStart, 0, 0) #total intensity vs x and y; consider adding different polarizations

                arParInf = [0]
                if((cw['xFWFM'] > 0) or (cw['yFWFM'] > 0)):
                    arParInf = [0, 0, 0]
                    if(cw['xFWFM'] > 0): arParInf[1] = v.om_ce['xFWFM']
                    if(cw['yFWFM'] > 0): arParInf[2] = v.om_ce['yFWFM']
                
                infI = srwl.UtiIntInf(arI, meshI, arParInf)
                
                if((cw['xFWHM'] > 0) and (cn['xFWHM'] > 0)): #horizontal size
                #if((cw[0] > 0) and (cn[0] > 0)): #horizontal size
                    xFWHM = infI[4]
                    rx = (xFWHM - ct['xFWHM'])/cn['xFWHM']
                    dCost = cw['xFWHM']*rx*rx
                    #if(infIsReq): s2pr += (' xFWHM='+frm).format(xFWHM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('xFWHM', xFWHM, dCost)
                    cost += dCost

                if((cw['yFWHM'] > 0) and (cn['yFWHM'] > 0)): #vertical size
                    yFWHM = infI[5]
                    ry = (yFWHM - ct['yFWHM'])/cn['yFWHM']
                    dCost = cw['yFWHM']*ry*ry
                    #if(infIsReq): s2pr += (' yFWHM='+frm).format(yFWHM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('yFWHM', yFWHM, dCost)
                    cost += dCost

                if((cw['xFWFM'] > 0) and (cn['xFWFM'] > 0)): #additional horizontal size
                    xFWFM = infI[7]
                    rx = (xFWFM - ct['xFWFM'])/cn['xFWFM']
                    dCost = cw['xFWFM']*rx*rx
                    #if(infIsReq): s2pr += (' xFWFM='+frm).format(xFWFM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('xFWHM', xFWHM, dCost)
                    cost += dCost

                if((cw['yFWFM'] > 0) and (cn['yFWFM'] > 0)): #additional vertical size
                    yFWFM = infI[8]
                    ry = (yFWFM - ct['yFWFM'])/cn['yFWFM']
                    dCost = cw['yFWFM']*ry*ry
                    #if(infIsReq): s2pr += (' yFWFM='+frm).format(yFWFM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('yFWFM', yFWFM, dCost)
                    cost += dCost

                if((cw['iM'] > 0) and (cn['iM'] > 0)): #peak intensity
                    maxI = infI[0]
                    rI = (maxI - ct['iM'])/cn['iM']
                    dCost = cw['iM']*rI*rI
                    #if(infIsReq): s2pr += (' iM='+frm).format(maxI) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('iM', maxI, dCost)
                    cost += dCost

            if((cw['xpFWHM'] > 0) or (cw['ypFWHM'] > 0)): #horizontal and vertical angular divergences
                srwl.SetRepresElecField(wfr, 'a')
                meshI = wfr.mesh
                arI = array('f', [0]*meshI.nx*meshI.ny) #"flat" 2D array to take intensity data
                srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, meshI.eStart, 0, 0) #total intensity vs x and y; consider adding different polarizations

                arParInf = [0]
                if((cw['xpFWFM'] > 0) or (cw['ypFWFM'] > 0)):
                    arParInf = [0, 0, 0]
                    if(cw['xpFWFM'] > 0): arParInf[1] = v.om_ce['xpFWFM']
                    if(cw['ypFWFM'] > 0): arParInf[2] = v.om_ce['ypFWFM']

                infI = srwl.UtiIntInf(arI, meshI, arParInf)
                
                if((cw['xpFWHM'] > 0) and (cn['xpFWHM'] > 0)): #horizontal angular divergence
                    xpFWHM = infI[4]
                    rxp = (xpFWHM - ct['xpFWHM'])/cn['xpFWHM']
                    dCost = cw['xpFWHM']*rxp*rxp
                    #if(infIsReq): s2pr += (' xpFWHM='+frm).format(xpFWHM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('xpFWHM', xpFWHM, dCost)
                    cost += dCost

                if((cw['ypFWHM'] > 0) and (cn['ypFWHM'] > 0)): #vertical angular divergence
                    ypFWHM = infI[5]
                    ryp = (ypFWHM - ct['ypFWHM'])/cn['ypFWHM']
                    dCost = cw['ypFWHM']*ryp*ryp
                    #if(infIsReq): s2pr += (' ypFWHM='+frm).format(ypFWHM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('ypFWHM', ypFWHM, dCost)
                    cost += dCost

                if((cw['xpFWFM'] > 0) and (cn['xpFWFM'] > 0)): #additional horizontal angular divergence
                    xpFWFM = infI[7]
                    rxp = (xpFWFM - ct['xpFWFM'])/cn['xpFWFM']
                    dCost = cw['xpFWFM']*rxp*rxp
                    #if(infIsReq): s2pr += (' xpFWFM='+frm).format(xpFWFM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('xpFWFM', xpFWFM, dCost)
                    cost += dCost

                if((cw['ypFWFM'] > 0) and (cn['ypFWFM'] > 0)): #additional vertical angular divergence
                    ypFWFM = infI[8]
                    ryp = (ypFWFM - ct['ypFWFM'])/cn['ypFWFM']
                    dCost = cw['ypFWFM']*ryp*ryp
                    #if(infIsReq): s2pr += (' ypFWFM='+frm).format(ypFWFM) + (' (sub-cost:'+frm+')').format(dCost)
                    if(infIsReq): s2pr += critInf('ypFWFM', ypFWFM, dCost)
                    cost += dCost

            if(cw['iD'] > 0): #arbitrary intensity distribution to fit; this optimization should run separately with others (const values don't accumulate)
                if(hasattr(self, 'arIntFit') and hasattr(self, 'meshIntFit')):

                    #To move this section to a separate function(?)
                    meshI2d = wfr.mesh
                    if(v.om_fp[0] == 3): #azimuthally-averaged flux per unit surface area
                        arI2d = array('f', [0]*meshI2d.nx*meshI2d.ny) #"flat" 2D array to take intensity data
                        srwl.CalcIntFromElecField(arI2d, wfr, 6, 0, 3, meshI2d.eStart, 0, 0) #total intensity vs x and y

                        xcAzAvg = v.om_fp[7]; ycAzAvg = v.om_fp[8]
                        rMaxFromI2d = xcAzAvg - meshI2d.xStart
                        rMaxTest = meshI2d.xFin - xcAzAvg
                        if(rMaxFromI2d > rMaxTest): rMaxFromI2d = rMaxTest
                        rMaxTest = ycAzAvg - meshI2d.yStart
                        if(rMaxFromI2d > rMaxTest): rMaxFromI2d = rMaxTest
                        rMaxTest = meshI2d.yFin - ycAzAvg
                        if(rMaxFromI2d > rMaxTest): rMaxFromI2d = rMaxTest

                        rStep = (self.meshIntFit.xFin - self.meshIntFit.xStart)/(self.meshIntFit.nx - 1)
                        rStart = self.meshIntFit.xStart
                        rFin = self.meshIntFit.xFin
                        irStart = 0
                        irFinP1 = self.meshIntFit.nx 

                        if(v.om_fr is not None):
                            if(v.om_fr[0] > self.meshIntFit.xStart): 
                                rStart = v.om_fr[0]
                                irStart = int((rStart - self.meshIntFit.xStart)/rStep)
                                rStart = self.meshIntFit.xStart + irStart*rStep

                            if(v.om_fr[1] < self.meshIntFit.xFin): 
                                rFin = v.om_fr[1]
                                irFinP1 = int((rFin - self.meshIntFit.xStart)/rStep) + 1
                                rFin = self.meshIntFit.xStart + (irFinP1 - 1)*rStep

                        if(rStart > rMaxFromI2d): raise Exception("Inconsistent argument limits of calculated and external / measured intensity distributions")

                        if(rFin > rMaxFromI2d): 
                            irFinP1 = int((rMaxFromI2d - self.meshIntFit.xStart)/rStep) + 1
                            rFin = self.meshIntFit.xStart + (irFinP1 - 1)*rStep

                        nr = irFinP1 - irStart
                        meshI1d = SRWLRadMesh(_eStart=meshI2d.eStart, _eFin=meshI2d.eFin, _ne=meshI2d.ne, 
                                              _xStart=rStart, _xFin=rFin, _nx=nr, 
                                              _yStart=0, _yFin=0, _ny=1, _zStart=meshI2d.zStart)
                        arI1d = array('f', [0]*nr)
                        srwl.UtiIntProc(arI1d, meshI1d, arI2d, meshI2d, v.om_fp)
                        #print('Azimuthally-averaged intensity calculated')

                        v.om_aux_rStart = rStart
                        v.om_aux_rFin = rFin
                        v.om_aux_irStart = irStart
                        v.om_aux_irFinP1 = irFinP1
                        v.om_aux_arI1d = arI1d

                        #OCTEST
                        print('cost_func: rStart=', rStart, ' rFin=', rFin)

                        if(v.om_fs[0] == 0): #do not scale data at fitting (i.e. do not multiply it by a best-fit constant)

                            dCost = self.cost_func_aux_int_distr([1.], None, v)

                        elif(v.om_fs[0] == 1): #do scale data at fitting (i.e. multiply it by a best-fit constant)

                            self.cost_func_aux_int_distr([1.], None, v)
                            ivScaleMult = v.om_aux_avg_fMeas/v.om_aux_avg_fCalc
                            #OCTEST
                            #print('cost_func: ivScaleMult=', ivScaleMult)

                            optRes = self.uti_math_opt.minimize(self.cost_func_aux_int_distr, 
                                                                _x=[ivScaleMult], _x_lim=[[0.01*ivScaleMult,100*ivScaleMult]], #to tune?
                                                                _meth=v.om_fs[0], _opt=v.om_mp, _aux=v)
                            dCost = optRes.fun
                            #OCTEST
                            #print('cost_func: optRes.fun=', optRes.fun)
                            #return optRes.fun #?

                    #if(infIsReq): s2pr += critInf('yFWHM', yFWHM, dCost)
                    cost += dCost

        if(v.wm): #partially-coherent / multi-electron calculations (add treatment of the Gaussian-Schell case later)

            if(self.comMPI is not None): self.comMPI.Barrier() #synchronizing all processes in case if MPI is used (e.g. for partially-coherent calculations)

            if(len(v.wm_fni) > 0): #load Intensity / Degree of Coherence from files and analyze the data

                if(((cw['xFWHM'] > 0) or (cw['yFWHM'] > 0) or (cw['iM'] > 0)) and ((v.wm_ap == 0) or (v.wm_ap == 2))): #horizontal, vertical sizes and peak intensity
                    if((v.wm_ch == 0) or (v.wm_ch == 40)):
                        arI, meshI = srwl_uti_read_intens_ascii(os.path.join(os.getcwd(), v.fdir, v.wm_fni), 'f')

                        arParInf = [0]
                        if((cw['xFWFM'] > 0) or (cw['yFWFM'] > 0)):
                            arParInf = [0, 0, 0]
                            if(cw['xFWFM'] > 0): arParInf[1] = v.om_ce['xFWFM']
                            if(cw['yFWFM'] > 0): arParInf[2] = v.om_ce['yFWFM']

                        infI = srwl.UtiIntInf(arI, meshI, arParInf)
            
                        if((cw['xFWHM'] > 0) and (cn['xFWHM'] > 0)): #horizontal size
                            xFWHM = infI[4]
                            rx = (xFWHM - ct['xFWHM'])/cn['xFWHM']
                            dCost = cw['xFWHM']*rx*rx
                            if(infIsReq): s2pr += critInf('xFWHM', xFWHM, dCost)
                            cost += dCost

                        if((cw['yFWHM'] > 0) and (cn['yFWHM'] > 0)): #vertical size
                            yFWHM = infI[5]
                            ry = (yFWHM - ct['yFWHM'])/cn['yFWHM']
                            dCost = cw['yFWHM']*ry*ry
                            if(infIsReq): s2pr += critInf('yFWHM', yFWHM, dCost)
                            cost += dCost

                        if((cw['xFWFM'] > 0) and (cn['xFWFM'] > 0)): #additional horizontal size
                            xFWFM = infI[7]
                            rx = (xFWFM - ct['xFWFM'])/cn['xFWFM']
                            dCost = cw['xFWFM']*rx*rx
                            if(infIsReq): s2pr += critInf('xFWFM', xFWFM, dCost)
                            cost += dCost

                        if((cw['yFWFM'] > 0) and (cn['yFWFM'] > 0)): #additional vertical size
                            yFWFM = infI[8]
                            ry = (yFWFM - ct['yFWFM'])/cn['yFWFM']
                            dCost = cw['yFWFM']*ry*ry
                            if(infIsReq): s2pr += critInf('yFWFM', yFWFM, dCost)
                            cost += dCost

                        if((cw['iM'] > 0) and (cn['iM'] > 0)): #peak intensity
                            maxI = infI[0]
                            rI = (maxI - ct['iM'])/cn['iM']
                            dCost = cw['iM']*rI*rI
                            if(infIsReq): s2pr += critInf('iM', maxI, dCost)
                            cost += dCost
            
                if(((cw['xpFWHM'] > 0) or (cw['ypFWHM'] > 0)) and ((v.wm_ap == 1) or (v.wm_ap == 2))): #horizontal and vertical angular divergences
                    #if((v.wm_ch == 0) or (v.wm_ch == 40)):
                    fnAngI = v.wm_fni
                    if(v.wm_ap == 2): fnAngI = srwl_wfr_fn(v.wm_fni, 2)
                    arI, meshI = srwl_uti_read_intens_ascii(os.path.join(os.getcwd(), v.fdir, fnAngI), 'f')

                    arParInf = [0]
                    if((cw['xpFWFM'] > 0) or (cw['ypFWFM'] > 0)):
                        arParInf = [0, 0, 0]
                        if(cw['xpFWFM'] > 0): arParInf[1] = v.om_ce['xpFWFM']
                        if(cw['ypFWFM'] > 0): arParInf[2] = v.om_ce['ypFWFM']

                    infI = srwl.UtiIntInf(arI, meshI, arParInf) #estimate FWHMs, starting search from extremities of the distribution

                    if((cw['xpFWHM'] > 0) and (cn['xpFWHM'] > 0)): #horizontal divergence
                        xpFWHM = infI[4]
                        rx = (xpFWHM - ct['xpFWHM'])/cn['xpFWHM']
                        dCost = cw['xpFWHM']*rx*rx
                        if(infIsReq): s2pr += critInf('xpFWHM', xpFWHM, dCost)
                        cost += dCost

                    if((cw['ypFWHM'] > 0) and (cn['ypFWHM'] > 0)): #vertical divergence
                        ypFWHM = infI[5]
                        ry = (ypFWHM - ct['ypFWHM'])/cn['ypFWHM']
                        dCost = cw['ypFWHM']*ry*ry
                        if(infIsReq): s2pr += critInf('ypFWHM', ypFWHM, dCost)
                        cost += dCost

                    if((cw['xpFWFM'] > 0) and (cn['xpFWFM'] > 0)): #additional horizontal angular divergence
                        xpFWFM = infI[7]
                        rxp = (xpFWFM - ct['xpFWFM'])/cn['xpFWFM']
                        dCost = cw['xpFWFM']*rxp*rxp
                        if(infIsReq): s2pr += critInf('xpFWFM', xpFWFM, dCost)
                        cost += dCost

                    if((cw['ypFWFM'] > 0) and (cn['ypFWFM'] > 0)): #additional vertical angular divergence
                        ypFWFM = infI[8]
                        ryp = (ypFWFM - ct['ypFWFM'])/cn['ypFWFM']
                        dCost = cw['ypFWFM']*ryp*ryp
                        if(infIsReq): s2pr += critInf('ypFWFM', ypFWFM, dCost)
                        cost += dCost

                if(((cw['xCL'] > 0) or (cw['yCL'] > 0)) and ((v.wm_ch == 4) or (v.wm_ch == 40))): #horizontal and vertical coherence lengths

                    if((cw['xCL'] > 0) and (cn['xCL'] > 0)): #horizontal coherence length
                        fnDC = srwl_wfr_fn(v.wm_fni, 41)
                        arDC, meshDC = srwl_uti_read_intens_ascii(os.path.join(os.getcwd(), v.fdir, fnDC), 'f')
                        infDC = srwl.UtiIntInf(arDC, meshDC, [1]) #estimate FWHMs, starting search from center of the distribution
                        xCoh = infDC[5] #cut vs (x1+x2)/2
                        rx = (xCoh - ct['xCL'])/cn['xCL']
                        dCost = cw['xCL']*rx*rx
                        if(infIsReq): s2pr += critInf('xCL', xCoh, dCost)
                        cost += dCost
                        
                    if((cw['yCL'] > 0) and (cn['yCL'] > 0)): #vertical coherence length
                        fnDC = srwl_wfr_fn(v.wm_fni, 42)
                        arDC, meshDC = srwl_uti_read_intens_ascii(os.path.join(os.getcwd(), v.fdir, fnDC), 'f')
                        infDC = srwl.UtiIntInf(arDC, meshDC, [1]) #estimate FWHMs, starting search from center of the distribution
                        yCoh = infDC[5] #cut vs (y1+y2)/2
                        ry = (yCoh - ct['yCL'])/cn['yCL']
                        dCost = cw['yCL']*ry*ry
                        if(infIsReq): s2pr += critInf('yCL', yCoh, dCost)
                        cost += dCost

            #if(cw['iD'] > 0): #arbitrary intensity distribution, to implement!

        if(v.om_res is not None): #OC08122022
            nPar_p_2 = nVars + 2
            lenStat = len(v.om_res)
            if(lenStat < nPar_p_2):
                extStat = [0]*(nPar_p_2 - lenStat)
                v.om_res.extend(extStat)

        if procIsMaster:
            if(v.om_pr or v.om_fl):

                statStr = self.uti_math_opt.status_str(v.om_pn, _x, cost, frm, v.om_res)
                if(v.om_pr): 
                    print(statStr)
                    sys.stdout.flush()
                if(v.om_fl): self.uti_math_opt.log_update(v.om_fnl, statStr)
            
        v.om_res = self.uti_math_opt.status_update(v.om_res, _x, cost)
        return cost

    #------------------------------------------------------------------------
    def cancel_calc_req(self, _v):
        """Calcel requests of any calculations in the parameter list _v
        :param _v: an object containing set of variables / options defining SR source and required calculations
        """
        _v.om = False
        _v.tr = False
        _v.ss = False
        _v.sm = False
        _v.pw = False
        _v.si = False
        _v.ws = False
        _v.wm = False
        return _v

    #------------------------------------------------------------------------
    def calc_all(self, _v, _op=None): #16122018
    #def calc_all(self, _v, _op):
        """Performs setup of electron beam, magnetic field, and performs calculations according to options specified in _v
        :param _v: an object containing set of variables / options defining SR source and required calculations
        :param _op: optical element container (SRWLOptC instance) that is assumed to be set up before calling this function and eventually used for wavefront propagation
        """
        
        #try: #OC23052020
        #    from mpi4py import MPI
        #    self.comMPI = MPI.COMM_WORLD
        #except:
        #    self.comMPI = None

        #---perform optimization of beamline and possibly source parameters (this should be processed before any other things)
        if(_v.om): #since the optimization modifies params in _v and calls calc_all, it goes first
            try: #16122018
                import uti_math_opt
                self.uti_math_opt = uti_math_opt
                
            except:
                traceback.print_exc()
                print('Optimization can not be performed (probably because required libraries are missing).')
                return

            try: #OC23052020
                from mpi4py import MPI
                self.MPI = MPI #OC28092022
                self.comMPI = MPI.COMM_WORLD
            except:
                self.comMPI = None

            _v.om = False #to avoid infinite nested loop

            #_v.om_cw = uti_math_opt.norm_weights(_v.om_cw)
            #Updating weight values in the hashtable (is there a simpler war to do the above manipulation?):
            om_cw_vals = uti_math_opt.norm_weights([_v.om_cw[k] for k in _v.om_cw.keys()])
            ik = 0
            for k in _v.om_cw.keys():
                _v.om_cw[k] = om_cw_vals[ik]
                ik += 1
            #print(_v.om_cw)

            if(_v.om_pr): 
                print('Optimization started ...')
                sys.stdout.flush()

            if(_v.om_fl): _v.om_fnl = self.uti_math_opt.log_init('__srwl_logs__')

            _v.om_res = [0]*(len(_v.om_iv) + 2) #array for storing instant optimization results in the process of the optimization
            for ip in range(len(_v.om_iv)): _v.om_res[ip] = _v.om_iv[ip]

            if(_v.om_cw['iD'] > 0):
                if((_v.om_fn is not None) and (len(_v.om_fn) > 0)):
                    fPathOptIntDistr = os.path.join(os.getcwd(), _v.fdir, _v.om_fn)
                    #print(self.meshIntFit)
                    self.arIntFit, self.meshIntFit = srwl_uti_read_intens_ascii(fPathOptIntDistr)
            
            uti_math_opt.minimize(self.cost_func, _v.om_iv, _x_lim=_v.om_lm, _meth=_v.om_mt, _opt=_v.om_mp, _aux=_v)
            _v = self.cancel_calc_req(_v)

        #---main folder
        if hasattr(_v, 'fdir'): self.dir_main = _v.fdir

        #---defining radiation source type (moved from def setup_source(v))
        if hasattr(_v, 'rs_type'): #OC29012019

            #_v.rs_tp == 'c' means combined SR source that may include several magnetic fields
            
            if((_v.rs_type != 'u') and (_v.rs_type != 'c')): #not idealized undulator / multi-pole wiggler (and e-beam)
                if hasattr(_v, 'und_b'): _v.und_b = 0
                if hasattr(_v, 'und_by'): _v.und_by = 0
                if hasattr(_v, 'und_bx'): _v.und_bx = 0

            if((_v.rs_type != 't') and (_v.rs_type != 'c')): #not ID with tabulated magnetic field, including dependence on gap / phase (and e-beam)
                if hasattr(_v, 'und_g'): del _v.und_g

            if((_v.rs_type != 'd') and (_v.rs_type != 'm') and (_v.rs_type != 'c')): #not idealized dipole magnet with eventual gradient (and e-beam)
            #if((_v.rs_type != 'd') and (_v.rs_type != 'c')): #not idealized dipole magnet with eventual gradient (and e-beam)
                if hasattr(_v, 'mag_bx'): _v.mag_bx = 0
                if hasattr(_v, 'mag_by'): _v.mag_by = 0
                if hasattr(_v, 'mag_gn'): _v.mag_gn = 0
                if hasattr(_v, 'mag_gs'): _v.mag_gs = 0

            if((_v.rs_type != 'a') and (_v.rs_type != 'c')): #not tabulated magnetic field not depending on gap / phase, i.e. arbitrary magnetic field source (and e-beam)
                #if hasattr(_v, 'mag_ifn'): del _v.mag_ifn
                if hasattr(_v, 'mag_fn'): del _v.mag_fn #OC29052019

            if(_v.rs_type != 'g'): #fully-coherent Gaussian radiation beam
                if hasattr(_v, 'gbm_pen'): del _v.gbm_pen

            if(_v.rs_type != 'p'): #point radiation source
                if hasattr(_v, 'psc_fl'): del _v.psc_fl

            #if(_v.rs_tp != 'gm'): #partially-coherent Gaussian Schell model radiation beam, defined through e-beam and coherent Gaussian beam
            #    pass
            #if(_v.rs_tp != 'f'): #FEL radiation beam in the format of GENESIS code output file(s)
            #    pass
            #if(_v.rs_tp != 'ws'): #fully-coherent SRW wavefront in pickle (or other) file format
            #    pass
            #if(_v.rs_tp != 'wm'): #partially-coherent SRW radiation represented as a set of fully-coherent wavefronts (format TBD)
            #    pass

        #---setup electron beam
        #if(hasattr(_v, 'ebm_nm')): #To improve
        #    self.set_e_beam(
        #        _e_beam_name = (_v.ebm_nm + _v.ebm_nms),
        #        _i = _v.ebm_i,
        #        _sig_e = _v.ebm_ens,
        #        _emit_x = _v.ebm_emx,
        #        _emit_y = _v.ebm_emy,
        #        _drift = _v.ebm_dr,
        #        _x = _v.ebm_x,
        #        _y = _v.ebm_y,
        #        _xp = _v.ebm_xp,
        #        _yp = _v.ebm_yp,
        #        _dE = _v.ebm_de)
        #    #Re-define some 2-nd order moments, if necessary:
        #    if(_v.ebm_sigx > 0): self.eBeam.arStatMom2[0] = (_v.ebm_sigx)*(_v.ebm_sigx)
        #    if(_v.ebm_mxxp != 1.e+23): self.eBeam.arStatMom2[1] = _v.ebm_mxxp
        #    if(_v.ebm_sigxp > 0): self.eBeam.arStatMom2[2] = (_v.ebm_sigxp)*(_v.ebm_sigxp)
        #    if(_v.ebm_sigy > 0): self.eBeam.arStatMom2[3] = (_v.ebm_sigy)*(_v.ebm_sigy)
        #    if(_v.ebm_myyp != 1.e+23): self.eBeam.arStatMom2[4] = _v.ebm_myyp
        #    if(_v.ebm_sigyp > 0): self.eBeam.arStatMom2[5] = (_v.ebm_sigyp)*(_v.ebm_sigyp)
        
        if hasattr(_v, 'ebm_nm'): #MR28092016
            #OC: to check if the above is the right condition
            self.set_e_beam(
                _e_beam_name=(_v.ebm_nm + _v.ebm_nms),
                _e_beam=None,
                _i=_v.ebm_i,
                _ens=_v.ebm_ens,
                _emx=_v.ebm_emx,
                _emy=_v.ebm_emy,
                _dr=_v.ebm_dr,
                _x=_v.ebm_x,
                _y=_v.ebm_y,
                _xp=_v.ebm_xp,
                _yp=_v.ebm_yp,
                _e=_v.ebm_e,
                _de=_v.ebm_de,
                # Twiss parameters:
                _betax=_v.ebm_betax,
                _alphax=_v.ebm_alphax,
                _etax=_v.ebm_etax,
                _etaxp=_v.ebm_etaxp,
                _betay=_v.ebm_betay,
                _alphay=_v.ebm_alphay,
                _etay=_v.ebm_etay,
                _etayp=_v.ebm_etayp,
                # Moments:
                _sigx=_v.ebm_sigx,
                _sigxp=_v.ebm_sigxp,
                _mxxp=_v.ebm_mxxp,
                _sigy=_v.ebm_sigy,
                _sigyp=_v.ebm_sigyp,
                _myyp=_v.ebm_myyp)
            
        #print('e-beam was set up')

        #---setup magnetic field: undulator, sinusoidal approximation
        #if hasattr(_v, 'und_b'):
        if(hasattr(_v, 'und_b') or hasattr(_v, 'und_by') or hasattr(_v, 'und_bx')): #OC25052016
            if hasattr(_v, 'und_bx') == False: _v.und_bx = 0
            if hasattr(_v, 'und_by') == False: _v.und_by = _v.und_b
            if hasattr(_v, 'und_phx') == False: _v.und_phx = 0
            if hasattr(_v, 'und_phy') == False: _v.und_phy = 0
            if hasattr(_v, 'und_zc') == False: _v.und_zc = 0
            if((_v.und_bx != 0) or (_v.und_by != 0)): #OC01062016
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
            self.mag = None
            
            if _v.und_b2e:
                e1 = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].get_E1(_en_elec=self.eBeam.partStatMom1.get_E(), _unit='eV')
                print('Fundamental Photon Energy:', srwl_uti_num_round(e1), 'eV') #check how it will work under IPython

            if _v.und_e2b:
                b = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].E1_2_B(_e1=_v.w_e, _en_elec=self.eBeam.partStatMom1.get_E())
                print('Magnetic Field Amplitude:', srwl_uti_num_round(b), 'T') #check how it will work under IPython

        if(hasattr(_v, 'und2_b') or hasattr(_v, 'und2_by') or hasattr(_v, 'und2_bx')): #OC03122016
            if not hasattr(_v, 'und2_bx'): _v.und2_bx = 0
            if not hasattr(_v, 'und2_by'): _v.und2_by = _v.und2_b
            if not hasattr(_v, 'und2_phx'): _v.und2_phx = 0
            if not hasattr(_v, 'und2_phy'): _v.und2_phy = 0
            if not hasattr(_v, 'und2_zc'): _v.und2_zc = 0
            if((_v.und2_bx != 0) or (_v.und2_by != 0)):
                self.set_und_sin(#setup second approximate undulator field parameters
                    _per = _v.und2_per,
                    _len = _v.und2_len,
                    _bx = _v.und2_bx,
                    _by = _v.und2_by,
                    _phx = _v.und2_phx,
                    _phy = _v.und2_phy,
                    _sx = _v.und2_sx,
                    _sy = _v.und2_sy,
                    _zc = _v.und2_zc,
                    _add = 1)

                if((_v.und2_cma != 0) and (_v.und2_cml > 0)): #setting-up canting magnet (to move to a separate function and make more general; treat soft edges ensuring corrct kick angle)
                    self.set_mag_kick(
                        _angx = _v.und2_cma,
                        _angy = 0,
                        _len = _v.und2_cml,
                        _led = _v.und2_cmd,
                        _zc = _v.und2_cmz,
                        _add = 1)

        #---setup magnetic field: undulator, tabulated (e.g. measured) magnetic field (a set of 3D magnetic field data files for different gaps / phases with a summary file)
        if hasattr(_v, 'und_mfz'): #OC23052020
            if len(_v.und_mfz) > 0:
                testPath = os.path.join(os.getcwd(), self.dir_main, _v.und_mfz)
                if os.path.exists(testPath):
                    z = zipfile.ZipFile(testPath) #OC23052020: from Sirepo script-generation code by P. Moeller, R. Nagner
                    rank = 0
                    comMPI = None
                    try:
                        from mpi4py import MPI
                        comMPI = MPI.COMM_WORLD
                        rank = comMPI.Get_rank()
                    except: pass
                    #print(z) #DEBUG
                    
                    #rank = 0 if self.comMPI is None else self.comMPI.Get_rank()
                    if rank == 0:
                        z.extractall(path=os.path.join(os.getcwd(), self.dir_main))
                        #z.extractall()
                        #print('Zip file unzipped') #DEBUG
                    
                    if comMPI is not None: comMPI.Barrier()
                    magMeasSummaryFileFound = False
                    for f in z.namelist():
                        if re.search(r'\.txt', f): #summary file
                            _v.und_mfs = os.path.basename(f)
                            _v.und_mdir = os.path.dirname(f) or './'
                            magMeasSummaryFileFound = True
                            break
                    if not magMeasSummaryFileFound:
                        raise RuntimeError('Magnetic measurement / simulation summary *.txt file can not be found')

        magnMeasDirExists = False
        if hasattr(_v, 'und_mdir'):
            self.dir_magn_meas = _v.und_mdir
            if(len(self.dir_magn_meas) > 0): magnMeasDirExists = True
            #print('self.dir_magn_meas = ', self.dir_magn_meas) #DEBUG

        #OC21012023 (commented-out the following, moved some vars to args of set_und_tab)
        #magnMeasSumFileExists = False
        if hasattr(_v, 'und_mfs'):
           self.fn_magn_meas_sum = _v.und_mfs
           if(magnMeasDirExists):
               testPath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
               magnMeasSumFileExists = os.path.exists(testPath)
               if not magnMeasSumFileExists: raise Exception('Tabulated undulator magnetic field summary file can not be found') #OC21012023 ??
               #print(testPath) #DEBUG

        if magnMeasDirExists: #OC21012023
        #if magnMeasSumFileExists and hasattr(_v, 'und_g'):

            if(_v.und_g > 0.):
                phase = 0.
                if hasattr(_v, 'und_ph'): phase = _v.und_ph
                phase_mode = 'p1'
                if hasattr(_v, 'und_phm'): phase_mode = _v.und_phm

                #print('Before set_und_tab()') #DEBUG
                self.set_und_tab(#setup undulator source from measured magnetic field data
                    _gap = _v.und_g,
                    _ph_mode = phase_mode,
                    _phase = phase,
                    _zc = _v.und_zc,
                    _interp_ord = _v.und_ior, #3, #1,
                    _meas_or_calc = 'm',
                    _per = _v.und_per,
                    _c1 = _v.und_c1,
                    _c2 = _v.und_c2,
                    _a = _v.und_a,
                    _dg_by_len = _v.und_dg/_v.und_len,
                    _y0 = _v.ebm_y + _v.ebm_yp*_v.und_zc - _v.und_dy, #this assumes that e-beam parameters are defined at z=0
                    _yp = _v.ebm_yp - _v.und_yp, #)
                    _fn_magn_sum = _v.und_mfs) #OC21012023

                if(_v.und2_g > 0.):
                    phase2 = 0.
                    if hasattr(_v, 'und2_ph'): phase2 = _v.und2_ph
                    phase2_mode = 'p1'
                    if hasattr(_v, 'und2_phm'): phase2_mode = _v.und2_phm

                    #print('Before 2nd set_und_tab()') #DEBUG
                    self.set_und_tab(#setup undulator source from measured magnetic field data
                        _gap = _v.und2_g,
                        _ph_mode = phase2_mode,
                        _phase = phase2,
                        _zc = _v.und2_zc,
                        _interp_ord = _v.und2_ior, #3, #1,
                        _meas_or_calc = 'm',
                        _per = _v.und2_per,
                        _c1 = _v.und2_c1,
                        _c2 = _v.und2_c2,
                        _a = _v.und2_a,
                        _dg_by_len = _v.und2_dg/_v.und2_len,
                        _y0 = _v.ebm_y + _v.ebm_yp*_v.und2_zc - _v.und2_dy, #this assumes that e-beam parameters are defined at z=0
                        _yp = _v.ebm_yp - _v.und2_yp, #)
                        _fn_magn_sum = _v.und2_mfs) #OC21012023

                    if(_v.und3_g > 0.):
                        phase3 = 0.
                        if hasattr(_v, 'und3_ph'): phase3 = _v.und3_ph
                        phase3_mode = 'p1'
                        if hasattr(_v, 'und3_phm'): phase3_mode = _v.und3_phm

                        #print('Before 3rd set_und_tab()') #DEBUG
                        self.set_und_tab(#setup undulator source from measured magnetic field data
                            _gap = _v.und3_g,
                            _ph_mode = phase3_mode,
                            _phase = phase3,
                            _zc = _v.und3_zc,
                            _interp_ord = _v.und3_ior, #3, #1,
                            _meas_or_calc = 'm',
                            _per = _v.und3_per,
                            _c1 = _v.und3_c1,
                            _c2 = _v.und3_c2,
                            _a = _v.und3_a,
                            _dg_by_len = _v.und3_dg/_v.und3_len,
                            _y0 = _v.ebm_y + _v.ebm_yp*_v.und3_zc - _v.und3_dy, #this assumes that e-beam parameters are defined at z=0
                            _yp = _v.ebm_yp - _v.und3_yp, #)
                            _fn_magn_sum = _v.und3_mfs) #OC21012023

                #if((_v.ss_mag == 1) or (_v.ss_mag == 1) or (_v.w_mag == 1) or (_v.tr_mag == 1)):
                if((_v.ss and (_v.ss_mag == 1)) or (_v.sm and (_v.sm_mag == 1)) or ((_v.ws or _v.wm) and (_v.w_mag == 1)) or (_v.tr and (_v.tr_mag == 1))):
                    #print('test field conversion')
                    maxPer = _v.und_per + 0.01
                    self.set_und_per_from_tab(
                        _rel_ac_thr=0.05,
                        _max_nh=7,
                        _max_per=maxPer)
                
                #self.mag_approx = None
                ##forcing using tabulated field for whatever calculaitons (?):
                #_v.ss_mag = 2
                #_v.w_mag = 2
                #_v.tr_mag = 2

        #---setup magnetic field of a dipole magnet
        if hasattr(_v, 'mag_bx') == False: _v.mag_bx = 0
        if hasattr(_v, 'mag_by') == False: _v.mag_by = 0
        if hasattr(_v, 'mag_gn') == False: _v.mag_gn = 0
        if hasattr(_v, 'mag_gs') == False: _v.mag_gs = 0
        if hasattr(_v, 'mag_len') == False: _v.mag_len = 1.5 #?
        if hasattr(_v, 'mag_led') == False: _v.mag_led = 0
        if hasattr(_v, 'mag_r') == False: _v.mag_r = 0
        if hasattr(_v, 'mag_zc') == False: _v.mag_zc = 0
        if((_v.mag_bx != 0) or (_v.mag_by != 0) or (_v.mag_gn != 0) or (_v.mag_gs != 0)):
            print('Before setting up multipole') #DEBUG
            self.set_mag_multipole(#setup dipole / quad magnet parameters
                _bx = _v.mag_bx,
                _by = _v.mag_by,
                _gn = _v.mag_gn,
                _gs = _v.mag_gs,
                _len = _v.mag_len,
                _led = _v.mag_led,
                _r = _v.mag_r,
                _zc = _v.mag_zc)
            #self.mag = None #OC16122016 (commented-out)
            
        #---setup arbitrary magnetic field source from one file of tabulated 3D magnetic field data
        if hasattr(_v, 'mag_fn'):
            if(len(_v.mag_fn) > 0):
                magPath = os.getcwd()
                if hasattr(self, 'dir_main'):
                    if(len(self.dir_main) > 0): magPath = os.path.join(magPath, self.dir_main)
                if hasattr(_v, 'mag_mdir'):
                    if(len(_v.mag_mdir) > 0): magPath = os.path.join(magPath, _v.mag_mdir)
                magPath = os.path.join(magPath, _v.mag_fn)

                #if(os.path.exists(magPath)):
                #    #print('')
                #    #print(magPath)
                #    #print('')
                #    ordInterp = 1
                #    if hasattr(_v, 'mag_ior'): ordInterp = _v.mag_ior
                #    self.set_mag_tab(#setup magnet from tabulated magnetic field data file
                #        _fpath = magPath,
                #        _zc = _v.mag_zc,
                #        _interp_ord = ordInterp)
                    
                if(not os.path.exists(magPath)): #AH03022019
                    raise Exception('Magnetic field file {} can not be found'.format(_v.mag_fn))
                ordInterp = 1
                if hasattr(_v, 'mag_ior'): ordInterp = _v.mag_ior
                self.set_mag_tab(#setup magnet from tabulated magnetic field data file
                    _fpath = magPath,
                    _zc = _v.mag_zc,
                    _interp_ord = ordInterp)

        #magnMeasDirExists = False
        #if hasattr(_v, 'mag_mdir'):
        #    #self.dir_magn_meas = _v.mag_mdir
        #    if(len(_v.mag_mdir) > 0):
        #        if hasattr(_v, 'mag_ifn'):
        #            magPath = os.path.join(os.getcwd(), self.dir_main, _v.mag_mdir, _v.mag_ifn)
        #           
        #            if(os.path.exists(magPath)):
        #                #print("")
        #                #print(magPath)
        #                #print("")
        #                self.set_mag_tab(#setup magnet from tabulated magnetic field data file
        #                    _fpath = magPath,
        #                    _zc = _v.mag_zc,
        #                    _interp_ord = 3)
        #                
        #                #self.mag_approx = None #OC16122016 (commented-out)
        #                #forcing using tabulated field for whatever calculaitons (?):
        #                #_v.ss_mag = 2
        #                #_v.w_mag = 2
        #                #_v.tr_mag = 2

        #---setup Gaussian beam
        if hasattr(_v, 'gbm_pen'):
            if(_v.gbm_pen > 0): #OC11102017
                self.set_gsn_beam(#setup Gaussiam beam (i.e. only define parameters, without calculating wavefront)
                    _x = _v.gbm_x,
                    _y = _v.gbm_y,
                    _z = _v.gbm_z,
                    _xp = _v.gbm_xp,
                    _yp = _v.gbm_yp,
                    _avgPhotEn = _v.gbm_ave,
                    _pulseEn = _v.gbm_pen,
                    _repRate = _v.gbm_rep,
                    _polar = _v.gbm_pol,
                    _sigX = _v.gbm_sx,
                    _sigY = _v.gbm_sy,
                    _sigT = _v.gbm_st,
                    _mx = _v.gbm_mx,
                    _my = _v.gbm_my,
                    _presCA = _v.gbm_ca,
                    _presFT = _v.gbm_ft)

        #---setup Point Source (i.e. only define parameters, without calculating spherical wavefront)
        if hasattr(_v, 'psc_fl'):
            if(_v.psc_fl > 0):
                self.set_pt_src(#setup Point Source (i.e. only define parameters, without calculating wavefront)
                    _x = _v.psc_x,
                    _y = _v.psc_y,
                    _z = _v.psc_z,
                    _flux = _v.psc_fl,
                    _unitFlux = _v.psc_ufl,
                    _polar = _v.psc_pol)

        #---setup detector (that may be used at different calculations)
        detector = None
        if((_v.d_rx > 0.) and (_v.d_nx > 0) and (_v.d_ry > 0.) and (_v.d_ny > 0)): #OC06122016
            detector = self.set_detector(
                _x = _v.d_x,
                _rx = _v.d_rx,
                _nx = _v.d_nx,
                _dx = _v.d_dx,
                _y = _v.d_y,
                _ry = _v.d_ry,
                _ny = _v.d_ny,
                _dy = _v.d_dy,
                _ord = _v.d_or,
                _fname = os.path.join(_v.fdir, _v.d_ifn) if(len(_v.d_ifn) > 0) else '')

        #---calculate / tabulate magnetic field
        if(_v.mf): #OC26012023
            _v.mf_res = self.calc_mag_fld(
                _z = _v.mf_z, _rz = _v.mf_rz, _nz = _v.mf_nz,
                _x = _v.mf_x, _rx = _v.mf_rx, _nx = _v.mf_nx,
                _y = _v.mf_y, _ry = _v.mf_ry, _ny = _v.mf_ny,
                _mag_type = _v.mf_mag,
                _ior = _v.mf_ior,
                _fname = os.path.join(_v.fdir, _v.mf_fn) if(len(_v.mf_fn) > 0) else '')

        #---calculate electron trajectory
        if(_v.tr):
            #print(self.eBeam.partStatMom1.z)
            #trj = self.calc_el_trj(
            _v.tr_res = self.calc_el_trj( #OC15112017
                _ctst = _v.tr_cti, _ctfi = _v.tr_ctf, _np = _v.tr_np,
                _mag_type = _v.tr_mag,
                _fname = os.path.join(_v.fdir, _v.tr_fn) if(len(_v.tr_fn) > 0) else '')

        ##---setup detector (that may be used at different calculations)
        #detector = None
        #if((_v.d_rx > 0.) and (_v.d_nx > 0) and (_v.d_ry > 0.) and (_v.d_ny > 0)): #OC06122016
        #    detector = self.set_detector(
        #        _x = _v.d_x,
        #        _rx = _v.d_rx,
        #        _nx = _v.d_nx,
        #        _dx = _v.d_dx,
        #        _y = _v.d_y,
        #        _ry = _v.d_ry,
        #        _ny = _v.d_ny,
        #        _dy = _v.d_dy,
        #        _ord = _v.d_or,
        #        _fname = os.path.join(_v.fdir, _v.d_ifn) if(len(_v.d_ifn) > 0) else '')

        #---calculate single-e spectrum vs photon energy
        if(_v.ss or _v.gs): 
            mesh_ss = SRWLRadMesh(
                _v.ss_ei, _v.ss_ef, _v.ss_ne,
                _v.ss_x, _v.ss_x, 1,
                _v.ss_y, _v.ss_y, 1,
                _v.op_r)
            
            srCanBeCalc = (self.eBeam is not None) and ((self.mag_approx is not None) or (self.mag is not None))
            #print("                    _v.ss=", _v.ss)
            #print("                    _v.gs=", _v.gs)
            
            if((_v.gs != True) and (srCanBeCalc == True)):
                #print("                                  Before calc_sr_se")
                #print("                                  self.mag=", self.mag)
                #wfr_ss, int_ss = self.calc_sr_se(
                #wfr_ss, int_ss, mesh_dummy = self.calc_sr_se( #OC06122016
                _v.w_res, _v.ss_res, mesh_dummy = self.calc_sr_se( #OC16102017
                    _mesh = mesh_ss,
                    _meth = _v.ss_meth,
                    _rel_prec = _v.ss_prec,
                    _zi = _v.ss_zi,
                    _zf = _v.ss_zf,
                    _pol = _v.ss_pol,
                    _int_type = 0,
                    _mag_type = _v.ss_mag,
                    _fname = os.path.join(_v.fdir, _v.ss_fn) if(len(_v.ss_fn) > 0) else '')

            if((_v.gs == True) or ((self.gsnBeam is not None) and (srCanBeCalc == False))):
                #wfr_ss, int_ss = self.calc_rad_gsn(
                #wfr_ss, int_ss, mesh_dummy = self.calc_rad_gsn( #OC06122016
                _v.w_res, _v.ss_res, mesh_dummy = self.calc_rad_gsn( #OC16102017
                    _mesh = mesh_ss,
                    _pol = _v.ss_pol,
                    _int_type = 0,
                    _presFT = _v.ss_ft,
                    _unitE = _v.ss_u,
                    _fname = os.path.join(_v.fdir, _v.ss_fn) if(len(_v.ss_fn) > 0) else '')
            
        #---calculate multi-e spectrum vs photon energy
        if(_v.sm):
            mesh_sm = SRWLRadMesh(
                _v.sm_ei, _v.sm_ef, _v.sm_ne,
                _v.sm_x - 0.5*_v.sm_rx, _v.sm_x + 0.5*_v.sm_rx, _v.sm_nx,
                _v.sm_y - 0.5*_v.sm_ry, _v.sm_y + 0.5*_v.sm_ry, _v.sm_ny,
                _v.op_r)
            if((_v.sm_mag == 1) and (_v.sm_meth < 0)):
                #int_sm = self.calc_ur_spec_me(
                _v.sm_res = self.calc_ur_spec_me( #OC16102017
                    _mesh = mesh_sm,
                    _harm_init = _v.sm_hi,
                    _harm_fin = _v.sm_hf,
                    _prec_long = _v.sm_prl,
                    _prec_azim = _v.sm_pra,
                    _type = _v.sm_type,
                    _pol = _v.sm_pol,
                    _fname = os.path.join(_v.fdir, _v.sm_fn) if(len(_v.sm_fn) > 0) else '')
            else: #for "accurate" magnetic field
                #int_sm = self.calc_arb_spec_me(
                _v.sm_res = self.calc_arb_spec_me( #OC16102017
                    _mesh = mesh_sm,
                    _meth = _v.sm_meth,
                    _rel_prec = _v.sm_prec,
                    _n_part_tot = _v.sm_nm,
                    _n_part_avg_proc = _v.sm_na,
                    _n_save_per = _v.sm_ns,
                    _type = _v.sm_type,
                    _mag = _v.sm_mag,
                    _pol = _v.sm_pol,
                    _rand_meth = _v.sm_rm,
                    _fname = os.path.join(_v.fdir, _v.sm_fn) if(len(_v.sm_fn) > 0) else '',
                    _sr_samp_fact = _v.sm_smpf,
                    _det = detector,
                    _me_approx = _v.sm_am, #) #OC13042018
                    _fbk = True if(_v.sm_fbk) else False) #OC14082018

        #---calculate undulator "operation table", i.e. dependence of gap (and phase) on photon energy (for a given polarization)
        if(_v.ut):
            #phase_mode = 'p1'
            #if hasattr(_v, 'und_phm'): phase_mode = _v.und_phm
            mesh_sm = SRWLRadMesh(
                _v.sm_ei, _v.sm_ef, _v.sm_ne,
                _v.sm_x - 0.5*_v.sm_rx, _v.sm_x + 0.5*_v.sm_rx, _v.sm_nx,
                _v.sm_y - 0.5*_v.sm_ry, _v.sm_y + 0.5*_v.sm_ry, _v.sm_ny,
                _v.op_r)
            #print('_v.sm_rx=', _v.sm_rx, '_v.sm_ry=', _v.sm_ry)
            
            self.calc_und_oper_tab(
                _mesh = mesh_sm,
                _pol = _v.sm_pol,
                _hi = _v.sm_hi,
                _hf = _v.sm_hf,
                _meas_or_calc = 'm',
                _zc = _v.und_zc,
                _fname = os.path.join(_v.fdir, _v.ut_fn) if(len(_v.ut_fn) > 0) else '')

        #---calculate SR power density distribution
        if(_v.pw):
            mesh_pw = SRWLRadMesh(
                1000, 1000, 1, #dummy (not used for power density)
                _v.pw_x - 0.5*_v.pw_rx, _v.pw_x + 0.5*_v.pw_rx, _v.pw_nx,
                _v.pw_y - 0.5*_v.pw_ry, _v.pw_y + 0.5*_v.pw_ry, _v.pw_ny,
                _v.op_r)
            #int_pw = self.calc_pow_den(
            _v.pw_res = self.calc_pow_den( #OC16102017
                _mesh = mesh_pw,
                _prec = _v.pw_pr,
                _meth = _v.pw_meth,
                #_z_start = _v.pw_zst,
                #_z_fin = _v.pw_zfi,
                _z_start = _v.pw_zi,
                _z_fin = _v.pw_zf,
                _mag_type = _v.pw_mag,
                _fname= os.path.join(_v.fdir, _v.pw_fn) if(len(_v.pw_fn) > 0) else '')

        #---calculate beamline optics orientations and save it to file
        if(len(_v.op_fno) > 0): #OC23032020
            avgPhEn = _v.w_e
            if(_v.w_ef > 0): avgPhEn = 0.5*(_v.w_e + _v.w_ef)
            self.set_optics(_op, _v)
            if(self.optics is not None): #OC28032020
                opOrientData = self.optics.get_orient_lab_fr(_e=avgPhEn, _r=_v.op_r)
                if(opOrientData is not None): #OC29012020
                    if(len(opOrientData) > 0):
                        uti_io.write_ascii_data_rows(_file_path=os.path.join(_v.fdir, _v.op_fno), _rows=opOrientData, _str_sep='\t',
                                                     _str_head='#Types of optical elements and Cartesian coordinates of their center positions and base vectors (t, s, n) in the Lab frame')
            
        #---calculate single-e and multi-e intensity distributions (before and after wavefront propagation through a beamline)
        if(_v.si or _v.ws or _v.gi or _v.wg or _v.wm):
            #if(_v.ws or _v.wg or _v.wm): self.set_optics(_op)
            if(_v.ws or _v.wg or _v.wm):

                #self.set_optics(_op, _v) #OC28042018
                needSetOptics = False #OC23032020
                if not hasattr(self, 'optics'): needSetOptics = True
                elif self.optics is None: needSetOptics = True
                elif hasattr(_v, 'op_func'): #OC28092022 (to make sure that Optics is updated when optimization is running)
                    if(_v.op_func is not None): needSetOptics = True 

                if needSetOptics: self.set_optics(_op, _v)
                
            ef = _v.w_e
            if(_v.w_ef > 0): ef = _v.w_ef
            mesh_w = SRWLRadMesh(
                _v.w_e, ef, _v.w_ne,
                _v.w_x - 0.5*_v.w_rx, _v.w_x + 0.5*_v.w_rx, _v.w_nx,
                _v.w_y - 0.5*_v.w_ry, _v.w_y + 0.5*_v.w_ry, _v.w_ny,
                _v.op_r)
        #---calculate single-e electric field and intensity (before wavefront propagation through a beamline)
            if(_v.si or _v.ws or _v.gi or _v.wg):

                srCanBeCalc = (self.eBeam is not None) and ((self.mag_approx is not None) or (self.mag is not None))
                gsnBeamCanBeCalc = (self.gsnBeam is not None)
                
                ptSrcSphWaveCanBeCalc = False
                if(hasattr(self, 'ptSrc')): ptSrcSphWaveCanBeCalc = (self.ptSrc is not None)

                #print(self.gsnBeam)
                detForSI = None if((_v.wg is True) or (_v.ws is True)) else detector
                #print('detForSI=', detForSI)

                #OC05042018
                wfrWasNotLoaded = True
                if(len(_v.ws_fnei) > 0):
                    print('Loading initial wavefront data from file ... ', end='')
                    t0 = time.time()
                    fnWfr = os.path.join(_v.fdir, _v.ws_fnei)
                    in_s = open(fnWfr, 'rb')
                    _v.w_res = pickle.load(in_s)
                    in_s.close()
                    if(_v.w_res is not None):
                        _v.si_res, mesh_si = self.calc_int_from_wfr(_v.w_res, _v.si_pol, _v.si_type, detForSI, _pr=False)
                        wfrWasNotLoaded = False
                        #print('completed (lasted', round(time.time() - t0, 2), 's)')
                        print('completed (lasted', round(time.time() - t0, 3), 's)')
                    else: print('failed to load wavefront file')

                if wfrWasNotLoaded:
                    #if((_v.gi == False) and (_v.wg == False) and (srCanBeCalc == True)):
                    if((_v.gi != True) and (_v.wg != True) and (srCanBeCalc == True)):

                        #print('Before wfr, int_w0 = self.calc_sr_se')
                
                        #wfr, int_w0 = self.calc_sr_se(
                        #wfr, int_w0, mesh_si = self.calc_sr_se( #OC06122016
                        _v.w_res, _v.si_res, mesh_si = self.calc_sr_se( #OC16102017
                            _mesh = deepcopy(mesh_w),
                            _samp_fact = _v.w_smpf,
                            _meth = _v.w_meth,
                            _rel_prec = _v.w_prec,
                            _zi = _v.w_zi,
                            _zf = _v.w_zf,
                            _pol = _v.si_pol,
                            _int_type = _v.si_type,
                            _mag_type = _v.w_mag,
                            _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '',
                            _det = detForSI)

                    #if((_v.gs == True) or ((gsnBeamCanBeCalc == True) and (srCanBeCalc == False))):
                    if((_v.gs == True) or (_v.wg == True) or ((gsnBeamCanBeCalc == True) and (srCanBeCalc == False))): #OC01062016
                        #wfr, int_w0 = self.calc_rad_gsn(
                        #wfr, int_w0, mesh_si = self.calc_rad_gsn( #OC06122016
                        _v.w_res, _v.si_res, mesh_si = self.calc_rad_gsn( #OC16102017
                            _mesh = deepcopy(mesh_w),
                            _samp_fact = _v.w_smpf,
                            _pol = _v.si_pol,
                            _int_type = _v.si_type,
                            _presFT = _v.w_ft,
                            _unitE = _v.w_u,
                            _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '',
                            _det = detForSI)

                    if((ptSrcSphWaveCanBeCalc is True) and (gsnBeamCanBeCalc is not True) and (srCanBeCalc is not True)): #OC11102017
                        #wfr, int_w0, mesh_si = self.calc_rad_pt_src(
                        _v.w_res, _v.si_res, mesh_si = self.calc_rad_pt_src( #OC16102017
                            _mesh = deepcopy(mesh_w),
                            _samp_fact = _v.w_smpf,
                            _pol = _v.si_pol,
                            _int_type = _v.si_type,
                            _presFT = _v.w_ft,
                            _unitE = _v.w_u,
                            _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '',
                            _det = detForSI)

                    if((len(_v.ws_fne) > 0) and (_v.w_res is not None)): #OC05042018
                        #Dumping initial wavefront
                        print('Saving initial wavefront data to a file ... ', end='')
                        t0 = time.time();
                        fnWfr = os.path.join(_v.fdir, _v.ws_fne)
                        out_s = open(fnWfr, 'wb')
                        pickle.dump(_v.w_res, out_s)
                        out_s.flush()
                        out_s.close()
                        #print('completed (lasted', round(time.time() - t0, 2), 's)')
                        print('completed (lasted', round(time.time() - t0, 3), 's)')

                #mesh_si = deepcopy(wfr.mesh) #OC06122016 (commented-out)
                #if(detForSI is None): mesh_si = deepcopy(wfr.mesh)

                #print('mesh_si.eStart=', mesh_si.eStart, 'mesh_si.eFin=', mesh_si.eFin, 'mesh_si.ne=', mesh_si.ne)
                #print('mesh_si.xStart=', mesh_si.xStart, 'mesh_si.xFin=', mesh_si.xFin, 'mesh_si.nx=', mesh_si.nx)
                #print('mesh_si.yStart=', mesh_si.yStart, 'mesh_si.yFin=', mesh_si.yFin, 'mesh_si.ny=', mesh_si.ny)
                
        #---calculate single-e electric field and intensity (after wavefront propagation through a beamline)
                if(_v.ws or _v.wg):
                #if(_v.ws or _v.wg or _v.wsm): #OC10052016 (commented-out)

                    if(_v.w_wr != 0.): #OC26032016
                        #wfr.Rx = _v.w_wr
                        #wfr.Ry = _v.w_wr
                        _v.w_res.Rx = _v.w_wr #OC16102017
                        _v.w_res.Ry = _v.w_wr

                    if(_v.w_wre > 0.): #OC05012017
                        #wfr.dRx = _v.w_wre
                        #wfr.dRy = _v.w_wre
                        _v.w_res.dRx = _v.w_wre #OC16102017
                        _v.w_res.dRy = _v.w_wre

                    #int_ws = self.calc_wfr_prop(
                    #int_ws, mesh_ws = self.calc_wfr_prop( #OC06122016
                    #_v.ws_res, mesh_ws = self.calc_wfr_prop( #OC16102017
                    _v.ws_res, mesh_ws, rad_intermed_ws = self.calc_wfr_prop( #OC08022021
                        #_wfr = wfr,
                        _wfr = _v.w_res, #OC16102017
                        _pres_ang = _v.ws_ap,
                        _pol = _v.si_pol,
                        _int_type = _v.si_type,
                        _dep_type=3, #consider adding other cases (e.g. for TD FEL calculations)
                        _fname = os.path.join(_v.fdir, _v.ws_fni) if(len(_v.ws_fni) > 0) else '',
                        _det = detector,
                        _rad_view = None if(not hasattr(_v, 'op_rv')) else _v.op_rv) #OC08022021
                        #_det = detector)
                    #mesh_ws = wfr.mesh #OC06122016 (commented-out)
                    #if(len(_v.ws_fn) > 0): to implement saving single-e (/ fully coherent) wavefront data (wfr) to a file

                    if((len(_v.ws_fnep) > 0) and (_v.w_res is not None)): #OC05042018
                        #Dumping propagated wavefront
                        print('Saving propagated wavefront data to a file ... ', end='')
                        t0 = time.time()
                        fnWfr = os.path.join(_v.fdir, _v.ws_fnep)
                        out_s = open(fnWfr, 'wb')
                        pickle.dump(_v.w_res, out_s)
                        out_s.flush()
                        out_s.close()
                        #print('completed (lasted', round(time.time() - t0, 2), 's)')
                        print('completed (lasted', round(time.time() - t0, 3), 's)')

        #---calculate multi-electron (/ partially coherent) wavefront propagation
            if(_v.wm):
                #wmResFileName = os.path.join(_v.fdir, _v.wm_fni) if(len(_v.wm_fni) > 0) else None
                #print(wmResFileName)
                #print('_v.wm_fbk=', _v.wm_fbk)

                wmCalcReq = True
                if(((_v.wm_ch == 7) or (_v.wm_ch == 71)) and (_v.wm_fnmi is not None)): #OC11102021 (calculation of CMD, or averaging CSD, based on (partial) CSD stored in an input file(s))
                #if((_v.wm_ch == 7) and (_v.wm_fnmi is not None)): #OC02072021 (calculation of CMD based on CSD stored in an input file)
                    if(len(_v.wm_fnmi) > 0):

                        fpCSD = os.path.join(_v.fdir, _v.wm_fnmi)

                        if(_v.wm_nmm > 1): #OC10102021 (CSD averaging is required)

                            fi_st = 0
                            fi_en = _v.wm_nmm - 1
                            if(_v.wm_fmis >= 0): 
                                fi_st = _v.wm_fmis
                                fi_en += fi_st

                            fpCore = fpCSD
                            fpExt = ''
                            iDot = fpCSD.rfind('.')
                            lenFP = len(fpCSD)
                            if((iDot >= 0) and (iDot >= lenFP - 5)):
                                fpCore = fpCSD[0:iDot]
                                fpExt = fpCSD[iDot:lenFP]

                            sTestMI = ''
                            lenCore = len(fpCore)
                            if(lenCore > 3):
                                sTestMI = fpCSD[lenCore-3:lenCore]
                                if(sTestMI == '_mi'): fpCore = fpCSD[0:lenCore-3]

                            fin_fpCore = None #OC11102021: processing the case _v.wm_ch == 71, when CSD just has to be averaged and saved to a new file
                            if(_v.wm_ch == 71): fin_fpCore = os.path.join(_v.fdir, _v.wm_fni) 

                            CSD, wfr0 = srwl_wfr_csd_avg(fpCore, fpExt, fi_st, fi_en, _fin_fp_core=fin_fpCore)
                            mesh = CSD.mesh

                        elif(_v.wm_ch == 7): #OC11102021
                        #else: #OC10102021 (made it after else)
                            data, mesh, dattr, wfr0 = srwl_uti_read_intens_hdf5(fpCSD)
                            CSD = SRWLStokes(_arS=data, _typeStokes='f', _mutual=1, _n_comp=1) #OC02072021: to cover different polarizations / Stokes components in the future
                            CSD.mesh = mesh

                        cohModes = None; eigVals = None #OC11102021
                        if(_v.wm_ch == 7): #OC11102021
                            cohModes, eigVals = srwl_wfr_cmd(CSD, _n_modes=_v.wm_ncm, _awfr=wfr0, _alg=_v.wm_acm)
                            fp_cm = os.path.join(_v.fdir, srwl_wfr_fn(_v.wm_fni, _type=7, _form='hdf5'))
                            srwl_uti_save_wfr_cm_hdf5(cohModes, None, _awfr=wfr0, _file_path=fp_cm)

                        _v.wm_res = cohModes, eigVals, mesh #?
                        wmCalcReq = False
                
                if(wmCalcReq): #OC02072021

                    lstWfrCM = None  #OC02072021
                    if(_v.wm_fncm is not None):
                        if(len(_v.wm_fncm) > 0):
                            fpCM = os.path.join(_v.fdir, _v.wm_fncm)
                            lstWfrCM = srwl_uti_read_wfr_cm_hdf5(fpCM)

                    #DEBUG
                    #print('_v.wm_nop=', _v.wm_nop)
                    #sys.stdout.flush()
                    #END DEBUG

                    _v.wm_res = self.calc_wfr_emit_prop_me( #OC19012021
                    #res_ipm = self.calc_wfr_emit_prop_me(
                        _mesh = mesh_w,
                        _sr_samp_fact = _v.w_smpf,
                        _sr_meth = _v.w_meth,
                        _sr_rel_prec = _v.w_prec,
                        _in_wr = _v.w_wr,
                        _in_wre = _v.w_wre, #OC05012017
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
                        _fname = os.path.join(_v.fdir, _v.wm_fni) if(len(_v.wm_fni) > 0) else None,
                        _det = detector,
                        #_multi_e_approx = _v.wm_am)
                        _me_approx = _v.wm_am, #) #OC13042018
                        _fbk = True if(_v.wm_fbk) else False, #) #OC14082018
                        _op_rnd = True if(_v.op_rnd) else False, #) #OC24042020
                        _fform = _v.wm_ff, #) #OC25022021
                        _no_opt = _v.wm_nop, #OC03072021
                        _nmm = _v.wm_nmm, #) #OC16042021
                        _ncm = _v.wm_ncm, #) #OC27062021
                        _cm_wfr = lstWfrCM, #) #OC02072021
                        _ms = _v.wm_ms, #) #OC22112022
                        _pol = _v.wm_pol) #OC28122023

        if(not srwl_uti_proc_is_master()): return #OC26122023 
        #ATTENTION: only plot by Master process after this point!

        #---plot results of all calculatiopns here (because the plotting "from the middle of the script" may hang up script execution)
        #uti_plot_init('TkAgg') #make the backend name an input option or move this to uti_plot ?
        plotOK = False

        #if (_v.tr == True) and (len(_v.tr_pl) > 0):
        if (_v.tr) and (len(_v.tr_pl) > 0): #OC06042018 (to make sure that it starts when _v.tr = 1)

            #args = []
            #kwargs = {
            #    'labels': ['Longitudinal Position'],
            #    'units': ['m', 'm'],
            #}
            ##OC: add more options: 'xpz', 'ypz', 'xy',...
            #traj_rep_allowed_values = ['xz', 'yz']
            #traj_rep_allowed_values += [x.upper() for x in traj_rep_allowed_values]

            #if _v.tr_pl.lower() == 'xz':
            #    args.append(trj.arX)
            #    kwargs['labels'].append('Horizontal Position')
            #elif _v.tr_pl.lower() == 'yz':
            #    args.append(trj.arY)
            #    kwargs['labels'].append('Vertical Position')
            #else:
            #    raise ValueError('No such option allowed: {}. Allowed values: {}'.format(_v.tr_pl, traj_rep_allowed_values))
            #kwargs['labels'].append('Electron Trajectory')
            #args.append([min(trj.arZ), max(trj.arZ), len(trj.arZ)])
            #uti_plot1d(*args, **kwargs)

            #OC15112017
            strSep = ',' #The possible separators are actually ',' and ' '
            strAuxOpt = copy(_v.tr_pl)
            strAuxOpt = strAuxOpt.replace(strSep, ' ')
            arOpt = strAuxOpt.split(' ')
            trj = _v.tr_res

            for i in range(len(arOpt)):
                curOpt = copy(arOpt[i])
                curOpt = curOpt.replace(' ', '')
                if(len(curOpt) > 0):
                    arPlotAbsc = None #reference to Abscissa array
                    arPlotOrd = None #reference to Ordinate array
                    gridIsReg = False
                    labels = None
                    units = None
                    curOpt = curOpt.lower()
                    if(curOpt == 'xz'):
                        arPlotAbsc = trj.arZ; arPlotOrd = trj.arX; gridIsReg = True
                        labels = ['Longitudinal Position', 'Horizontal Position', 'Horizontal Trajectory']
                        units = ['m', 'm']
                    elif(curOpt == 'xpz'):
                        arPlotAbsc = trj.arZ; arPlotOrd = trj.arXp; gridIsReg = True
                        labels = ['Longitudinal Position', 'Horizontal Angle', 'Horizontal Trajectory']
                        units = ['m', 'rad']
                    elif(curOpt == 'yz'):
                        arPlotAbsc = trj.arZ; arPlotOrd = trj.arY; gridIsReg = True
                        labels = ['Longitudinal Position', 'Vertical Position', 'Vertical Trajectory']
                        units = ['m', 'm']
                    elif(curOpt == 'ypz'):
                        arPlotAbsc = trj.arZ; arPlotOrd = trj.arYp; gridIsReg = True
                        labels = ['Longitudinal Position', 'Vertical Angle', 'Vertical Trajectory']
                        units = ['m', 'rad']
                    elif(curOpt == 'yx'):
                        arPlotAbsc = trj.arX; arPlotOrd = trj.arY
                        labels = ['Horizontal Position', 'Vertical Position', 'Transverse Trajectory']
                        units = ['m', 'm']
                    elif(curOpt == 'xy'):
                        arPlotAbsc = trj.arY; arPlotOrd = trj.arX
                        labels = ['Vertical Position', 'Horizontal Position', 'Transverse Trajectory']
                        units = ['m', 'm']
                    elif(curOpt == 'ypxp'):
                        arPlotAbsc = trj.arXp; arPlotOrd = trj.arYp
                        labels = ['Horizontal Angle', 'Vertical Angle', 'Transverse Trajectory']
                        units = ['rad', 'rad']
                    elif(curOpt == 'xpyp'):
                        arPlotAbsc = trj.arYp; arPlotOrd = trj.arXp
                        labels = ['Vertical Angle', 'Horizontal Angle', 'Transverse Trajectory']
                        units = ['rad', 'rad']
                    elif(curOpt == 'bxz'):
                        bxExists = False
                        if(hasattr(trj, 'arBx')): 
                            if(trj.arBx is not None):
                                arPlotAbsc = trj.arZ; arPlotOrd = trj.arBx; gridIsReg = True
                                labels = ['Longitudinal Position', 'Horizontal Magnetic Field', 'Magnetic Field']
                                units = ['m', 'T']
                                bxExists = True
                        if(not bxExists): raise ValueError('Horizontal magnetic field seen by electron was not calculated')
                    elif(curOpt == 'byz'):
                        byExists = False
                        if(hasattr(trj, 'arBy')): 
                            if(trj.arBy is not None):
                                arPlotAbsc = trj.arZ; arPlotOrd = trj.arBy; gridIsReg = True
                                labels = ['Longitudinal Position', 'Vertical Magnetic Field', 'Magnetic Field']
                                units = ['m', 'T']
                                byExists = True
                        if(not byExists): raise ValueError('Vertical magnetic field seen by electron was not calculated')

                    if((arPlotAbsc is not None) and (arPlotOrd is not None)):
                        if(gridIsReg): uti_plot1d(arPlotOrd, [min(arPlotAbsc), max(arPlotAbsc), len(arPlotAbsc)], labels, units)
                        else: uti_plot1d_ir(arPlotOrd, arPlotAbsc, labels, units)
                    else: raise ValueError('This trajectory plot option value is not supported: {}. The supported values are: xz,xpz,yz,ypz,yx,ypxp,bxz,byz'.format(curOpt))

            plotOK = True
        
        if((_v.ss or _v.gs) and (len(_v.ss_pl) > 0)):

            sArgLabel = 'Photon Energy'
            sArgUnit = 'eV'
            if(_v.ss_ft == 't'):
                sArgLabel = 'Time'
                sArgUnit = 's'
            
            sValLabel = 'Flux per Unit Surface'
            sValUnit = 'ph/s/.1%bw/mm^2'
            if(_v.w_u == 0):
                sValLabel = 'Intensity'
                sValUnit = 'a.u.'
            elif(_v.w_u == 2):
                if(_v.ss_ft == 't'):
                    sValLabel = 'Power Density'
                    sValUnit = 'W/mm^2'
                elif(_v.ss_ft == 'f'):
                    sValLabel = 'Spectral Fluence'
                    sValUnit = 'J/eV/mm^2'
            
            #uti_plot1d(int_ss, [mesh_ss.eStart, mesh_ss.eFin, mesh_ss.ne], [sArgLabel, sValLabel, sValLabel], [sArgUnit, sValUnit])
            uti_plot1d(_v.ss_res, [mesh_ss.eStart, mesh_ss.eFin, mesh_ss.ne], [sArgLabel, sValLabel, sValLabel], [sArgUnit, sValUnit]) #OC16102017

            plotOK = True

        #print(_v.sm)
        #print(_v.sm_pl)
        
        if _v.sm and (len(_v.sm_pl) > 0): #OC05082019
        #if (_v.sm == True) and (len(_v.sm_pl) > 0):

            #print('printing multi-e flux')
            sValType = 'Flux'; sValUnit = 'ph/s/.1%bw'
            if(_v.sm_type == 2):
                sValType = 'Flux per Unit Surface'; sValUnit = 'ph/s/.1%bw/mm^2'
            #uti_plot1d(int_sm, [mesh_sm.eStart, mesh_sm.eFin, mesh_sm.ne], ['Photon Energy', sValType, sValType], ['eV', sValUnit])
            
            if(isinstance(_v.sm_res, list)): #OC10102023
                arPolName = ['LH', 'LV', 'LT45', 'LT135', 'CR', 'CL']
                for i in range(6):
                    sValTypeP = sValType + ' at Polarization: ' + arPolName[i]
                    uti_plot1d(_v.sm_res[i], [mesh_sm.eStart, mesh_sm.eFin, mesh_sm.ne], ['Photon Energy', sValType, sValTypeP], ['eV', sValUnit])
            else:
                uti_plot1d(_v.sm_res, [mesh_sm.eStart, mesh_sm.eFin, mesh_sm.ne], ['Photon Energy', sValType, sValType], ['eV', sValUnit])

            plotOK = True

        if _v.pw and (len(_v.pw_pl) > 0): #OC01102022
        #if (_v.pw == True) and (len(_v.pw_pl) > 0):
            if ((_v.pw_pl == 'xy') or (_v.pw_pl == 'yx') or (_v.pw_pl == 'XY') or (_v.pw_pl == 'YX')) and (_v.pw_nx > 1) and (_v.pw_ny > 1):
                uti_plot2d1d(
                    _v.pw_res, #int_pw, #OC16102017
                    [mesh_pw.xStart, mesh_pw.xFin, mesh_pw.nx],
                    [mesh_pw.yStart, mesh_pw.yFin, mesh_pw.ny],
                    0.5*(mesh_pw.xStart + mesh_pw.xFin),
                    0.5*(mesh_pw.yStart + mesh_pw.yFin),
                    ['Horizontal Position', 'Vertical Position', 'Power Density'],
                    ['m', 'm', 'W/mm^2'], True)
            elif ((_v.pw_pl == 'x') or (_v.pw_pl == 'X')) and (_v.pw_nx > 1):
                #uti_plot1d(int_pw, [mesh_pw.xStart, mesh_pw.xFin, mesh_pw.nx],
                uti_plot1d(_v.pw_res, [mesh_pw.xStart, mesh_pw.xFin, mesh_pw.nx], #OC16102017
                           ['Horizontal Position', 'Power Density', 'Power Density'],
                           ['m', 'W/mm^2'])
            elif ((_v.pw_pl == 'y') or (_v.pw_pl == 'Y')) and (_v.pw_ny > 1):
                #uti_plot1d(int_pw, [mesh_pw.yStart, mesh_pw.yFin, mesh_pw.ny],
                uti_plot1d(_v.pw_res, [mesh_pw.yStart, mesh_pw.yFin, mesh_pw.ny], #OC16102017
                           ['Vertical Position', 'Power Density', 'Power Density'],
                           ['m', 'W/mm^2'])
            plotOK = True

        if (_v.si or _v.gs) and (len(_v.si_pl) > 0):
            if _v.si_pl in ['xy', 'yx', 'XY', 'YX']:

                sValLabel = 'Flux per Unit Surface'
                sValUnit = 'ph/s/.1%bw/mm^2'
                if(_v.w_u == 0):
                    sValLabel = 'Intensity'
                    sValUnit = 'a.u.'
                elif(_v.w_u == 2):
                    if(_v.w_ft == 't'):
                        sValLabel = 'Power Density'
                        sValUnit = 'W/mm^2'
                    elif(_v.w_ft == 'f'):
                        sValLabel = 'Spectral Fluence'
                        sValUnit = 'J/eV/mm^2'
                
                #print('testing', _v.si_pl)
                uti_plot2d1d(
                    #int_w0,
                    _v.si_res, #OC16102017
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel],
                    ['m', 'm', sValUnit], #add other units for FEL
                    True)
                plotOK = True

        #if _v.ws and (len(_v.ws_pl) > 0):
        if (_v.ws or _v.wg) and (len(_v.ws_pl) > 0): #OC01062016
            if (_v.ws_pl == 'xy') or (_v.ws_pl == 'yx') or (_v.ws_pl == 'XY') or (_v.ws_pl == 'YX'):
                #print('2D plot panel is to be prepared')
                
                sValLabel = 'Flux per Unit Surface'
                sValUnit = 'ph/s/.1%bw/mm^2'
                if(_v.w_u == 0):
                    sValLabel = 'Intensity'
                    sValUnit = 'a.u.'
                elif(_v.w_u == 2):
                    if(_v.w_ft == 't'):
                        sValLabel = 'Power Density'
                        sValUnit = 'W/mm^2'
                    elif(_v.w_ft == 'f'):
                        sValLabel = 'Spectral Fluence'
                        sValUnit = 'J/eV/mm^2'

                uti_plot2d1d(
                    #int_w0,
                    _v.si_res, #OC16102017
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel + ' Before Propagation'],
                    ['m', 'm', sValUnit],
                    True)
                uti_plot2d1d(
                    #int_ws,
                    _v.ws_res, #OC16102017
                    [mesh_ws.xStart, mesh_ws.xFin, mesh_ws.nx],
                    [mesh_ws.yStart, mesh_ws.yFin, mesh_ws.ny],
                    0, #0.5*(mesh_ws.xStart + mesh_ws.xFin),
                    0, #0.5*(mesh_ws.yStart + mesh_ws.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel + ' After Propagation'],
                    ['m', 'm', sValUnit],
                    True)

                if(rad_intermed_ws is not None): #OC08022021
                    for i in range(len(rad_intermed_ws)):
                        curIntData = rad_intermed_ws[i]
                        if(len(curIntData) > 6):
                            curMesh = curIntData[5]
                            uti_plot2d1d(
                                curIntData[6],
                                [curMesh.xStart, curMesh.xFin, curMesh.nx],
                                [curMesh.yStart, curMesh.yFin, curMesh.ny],
                                0, #0.5*(curMesh.xStart + curMesh.xFin),
                                0, #0.5*(curMesh.yStart + curMesh.yFin),
                                ['Horizontal Position', 'Vertical Position', sValLabel + ' After Elem. #' + repr(curIntData[0])],
                                ['m', 'm', sValUnit],
                                True)

        #if _v.wm and (len(_v.wm_pl) > 0): #OC19012021
        #    if (_v.wm_pl == 'xy') or (_v.wm_pl == 'yx') or (_v.wm_pl == 'XY') or (_v.wm_pl == 'YX'):
        #        #print('2D plot panel is to be prepared')
                
        #        sValLabel = 'Flux per Unit Surface'
        #        sValUnit = 'ph/s/.1%bw/mm^2'
        #        if(_v.w_u == 0):
        #            sValLabel = 'Intensity'
        #            sValUnit = 'a.u.'
        #        elif(_v.w_u == 2):
        #            if(_v.w_ft == 't'):
        #                sValLabel = 'Power Density'
        #                sValUnit = 'W/mm^2'
        #            elif(_v.w_ft == 'f'):
        #                sValLabel = 'Spectral Fluence'
        #                sValUnit = 'J/eV/mm^2'

        #        uti_plot2d1d(
        #            #int_w0,
        #            _v.wm_res.arS,
        #            [_v.wm_res.mesh.xStart, _v.wm_res.mesh.xFin, _v.wm_res.mesh.nx],
        #            [_v.wm_res.mesh.yStart, _v.wm_res.mesh.yFin, _v.wm_res.mesh.ny],
        #            0, #0.5*(mesh_si.xStart + mesh_si.xFin),
        #            0, #0.5*(mesh_si.yStart + mesh_si.yFin),
        #            ['Horizontal Position', 'Vertical Position', sValLabel + ' After Propagation'],
        #            ['m', 'm', sValUnit],
        #            True)

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
def srwl_uti_std_options():
    """Defines sets of standard default options (applicable to any beamline) for general types of calculation
    :returns: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    """
    varParamStd = [
#---Electron Beam
        ['ebm_nm', 's', 'NSLS-II Low Beta ', 'standard electron beam name'],
        ['ebm_nms', 's', 'Day1', 'standard electron beam name suffix: e.g. can be Day1, Final'],
        ['ebm_i', 'f', 0.5, 'electron beam current [A]'],
        ['ebm_e', 'f', 3., 'electron beam avarage energy [GeV]'],
        ['ebm_de', 'f', 0., 'electron beam average energy deviation [GeV]'],
        ['ebm_x', 'f', 0., 'electron beam initial average horizontal position [m]'],
        ['ebm_y', 'f', 0., 'electron beam initial average vertical position [m]'],
        ['ebm_xp', 'f', 0., 'electron beam initial average horizontal angle [rad]'],
        ['ebm_yp', 'f', 0., 'electron beam initial average vertical angle [rad]'],
        #['ebm_z', 'f', 0., 'electron beam initial average longitudinal position [m]'], #it is always assumed to be 0.
        ['ebm_dr', 'f', 0., 'electron beam longitudinal drift [m] to be performed before a required calculation'],
        ['ebm_ens', 'f', -1, 'electron beam relative energy spread'],
        ['ebm_emx', 'f', -1, 'electron beam horizontal emittance [m]'],
        ['ebm_emy', 'f', -1, 'electron beam vertical emittance [m]'],
        # Definition of the beam through Twiss parameters: #MR28092016
        ['ebm_betax', 'f', None, 'horizontal beta-function [m]'], #OC: re-check the default values
        ['ebm_alphax', 'f', None, 'horizontal alpha-function [rad]'],
        ['ebm_etax', 'f', None, 'horizontal dispersion function [m]'],
        ['ebm_etaxp', 'f', None, 'horizontal dispersion function derivative [rad]'],
        ['ebm_betay', 'f', None, 'vertical beta-function [m]'],
        ['ebm_alphay', 'f', None, 'vertical alpha-function [rad]'],
        ['ebm_etay', 'f', None, 'vertical dispersion function [m]'],
        ['ebm_etayp', 'f', None, 'vertical dispersion function derivative [rad]'],
        # Definition of the beam through Moments:
        ['ebm_sigx', 'f', None, 'horizontal RMS size of electron beam [m]'],
        ['ebm_sigy', 'f', None, 'vertical RMS size of electron beam [m]'],
        ['ebm_sigxp', 'f', None, 'horizontal RMS angular divergence of electron beam [rad]'],
        ['ebm_sigyp', 'f', None, 'vertical RMS angular divergence of electron beam [rad]'],
        ['ebm_mxxp', 'f', None, 'horizontal position-angle mixed 2nd order moment of electron beam [m]'],
        ['ebm_myyp', 'f', None, 'vertical position-angle mixed 2nd order moment of electron beam [m]'],
        #['ebm_sigx', 'f', -1, 'horizontal RMS size of electron beam [m] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],
        #['ebm_sigy', 'f', -1, 'vertical RMS size of electron beam [m] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],
        #['ebm_sigxp', 'f', -1, 'horizontal RMS angular divergence of electron beam [rad] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],
        #['ebm_sigyp', 'f', -1, 'vertical RMS angular divergence of electron beam [rad] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],
        #['ebm_mxxp', 'f', 1.e+23, 'horizontal position-angle mixed 2nd order moment of electron beam [m] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],
        #['ebm_myyp', 'f', 1.e+23, 'vertical position-angle mixed 2nd order moment of electron beam [m] (is taken into account if > 0, in that case it overrides Emittance and Twiss parameters)'],

#---Undulator
        ['und_per', 'f', 0.02, 'undulator period [m]'],
        ['und2_per', 'f', 0.02, 'undulator period [m]'],
        ['und3_per', 'f', 0.02, 'undulator period [m]'],
        ['und_len', 'f', 3., 'undulator length [m]'],
        ['und2_len', 'f', 3., 'undulator length [m]'],
        ['und3_len', 'f', 3., 'undulator length [m]'],
        ['und_b', 'f', 0., 'undulator vertical peak magnetic field [T]'], #Keeping it 0 here is important for parsing calculation options!
        ['und2_b', 'f', 0., 'undulator vertical peak magnetic field [T]'], #Keeping it 0 here is important for parsing calculation options!
        ['und3_b', 'f', 0., 'undulator vertical peak magnetic field [T]'], #Keeping it 0 here is important for parsing calculation options!
        ['und_bx', 'f', 0., 'undulator horizontal peak magnetic field [T]'],
        ['und2_bx', 'f', 0., 'undulator horizontal peak magnetic field [T]'],
        ['und3_bx', 'f', 0., 'undulator horizontal peak magnetic field [T]'],
        #['und_by', 'f', 0., 'undulator vertical peak magnetic field [T]'],
        #['und2_by', 'f', 0., 'undulator vertical peak magnetic field [T]'],
        #['und3_by', 'f', 0., 'undulator vertical peak magnetic field [T]'],
        ['und_g', 'f', 0., 'undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und2_g', 'f', 0., 'second undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und3_g', 'f', 0., 'third undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],

        ['und_ph', 'f', 0., 'undulator phase, i.e. longitudinal shift of magnet arrays [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und2_ph', 'f', 0., 'second undulator phase, i.e. longitudinal shift of magnet arrays [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und3_ph', 'f', 0., 'third undulator phase, i.e. longitudinal shift of magnet arrays [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und_phm', 's', 'p1', 'undulator phase move mode'],
        ['und2_phm', 's', 'p1', 'second undulator phase move mode'],
        ['und3_phm', 's', 'p1', 'third undulator phase move mode'],

        #['und_ior', 'f', 1, 'interpolaton order of tabulated undulator magnetic field'],
        ['und_ior', 'f', 3, 'interpolaton order of tabulated undulator magnetic field'], #OC03082019
        #['und2_ior', 'f', 1, 'interpolaton order of second tabulated undulator magnetic field'],
        ['und2_ior', 'f', 3, 'interpolaton order of second tabulated undulator magnetic field'], #OC03082019
        ['und3_ior', 'f', 3, 'interpolaton order of third tabulated undulator magnetic field'], #OC20012023

        ['und_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
        ['und2_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
        ['und3_sx', 'i', 1, 'undulator horizontal magnetic field symmetry vs longitudinal position'],
        ['und_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
        ['und2_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
        ['und3_sy', 'i', -1, 'undulator vertical magnetic field symmetry vs longitudinal position'],
        ['und_zc', 'f', 0., 'undulator center longitudinal position [m]'],
        ['und2_zc', 'f', 0., 'second undulator center longitudinal position [m]'],
        ['und3_zc', 'f', 0., 'third undulator center longitudinal position [m]'],

        ['und2_cma', 'f', 0., 'canting magnet angle [rad]'],
        ['und3_cma', 'f', 0., 'canting magnet angle [rad]'],
        ['und2_cmz', 'f', 0., 'canting magnet longitudinal position [m]'],
        ['und3_cmz', 'f', 0., 'canting magnet longitudinal position [m]'],
        ['und2_cml', 'f', 0., 'canting magnet effective length [m]'],
        ['und3_cml', 'f', 0., 'canting magnet effective length [m]'],
        ['und2_cmd', 'f', 0., 'canting magnet edge length [m]'],
        ['und3_cmd', 'f', 0., 'canting magnet edge length [m]'],

        ['und_b0', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. b0 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und2_b0', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. b0 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und3_b0', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. b0 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und_c1', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und2_c1', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und3_c1', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c1 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und_c2', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und2_c2', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und3_c2', 'f', 0., 'constant defining (approximate) undulator field dependence on gap (i.e. c2 in b0*exp(-c1*gap/per + c2*(gap/per)^2))'],
        ['und_a', 'f', 0., 'constant defining (approximate) undulator field dependence on vertical position (i.e. a in cosh(2*Pi*a*y/per)'],
        ['und2_a', 'f', 0., 'constant defining (approximate) undulator field dependence on vertical position (i.e. a in cosh(2*Pi*a*y/per)'],
        ['und3_a', 'f', 0., 'constant defining (approximate) undulator field dependence on vertical position (i.e. a in cosh(2*Pi*a*y/per)'],

        ['und_dg', 'f', 0., 'undulator gap taper, i.e. gap difference between exit and entrance [m]'],
        ['und2_dg', 'f', 0., 'undulator gap taper, i.e. gap difference between exit and entrance [m]'],
        ['und3_dg', 'f', 0., 'undulator gap taper, i.e. gap difference between exit and entrance [m]'],
        ['und_dy', 'f', 0., 'undulator elevation in vertical direction over the median plane [m]'],
        ['und2_dy', 'f', 0., 'undulator elevation in vertical direction over the median plane [m]'],
        ['und3_dy', 'f', 0., 'undulator elevation in vertical direction over the median plane [m]'],
        ['und_yp', 'f', 0., 'undulator vertical angular misalignment over the median plane [rad]'],
        ['und2_yp', 'f', 0., 'undulator vertical angular misalignment over the median plane [rad]'],
        ['und3_yp', 'f', 0., 'undulator vertical angular misalignment over the median plane [rad]'],
        
        ['und_mdir', 's', '', 'name of magnetic measurements sub-folder'], #OC23072021
        #['und_mdir', 's', 'magn_meas', 'name of magnetic measurements sub-folder'],
        ['und_mfs', 's', '', 'name of undulator magnetic measurements for different gaps summary file'],
        ['und2_mfs', 's', '', 'name of second undulator magnetic measurements for different gaps summary file'],
        ['und3_mfs', 's', '', 'name of third undulator magnetic measurements for different gaps summary file'],
        ['und_mfz', 's', '', 'name of zip-file of directory with magnetic measurement files for different gaps + summary file (if it is defined, it overrides the values of und_mdir and und_mfs)'],

        ['und_b2e', '', '', 'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by', 'store_true'],
        ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e', 'store_true'],

#---Dipole Magnet (with eventual field gradient)
        ['mag_by', 'f', 0., 'vertical magnetic field [T]'],
        ['mag_bx', 'f', 0., 'horizontal magnetic field [T]'],
        ['mag_gn', 'f', 0., 'gradient of normal magnetic field quadrupole [T/m]'],
        ['mag_gs', 'f', 0., 'gradient of skew magnetic field quadrupole [T/m]'],
        ['mag_len', 'f', 3.0, 'magnetic length [m]'],
        ['mag_led', 'f', 0., 'magnetic edge length bw 10% and 90% of peak field [m]; G/(1 + ((z-zc)/d)^2)^2 field dependence is assumed'],
        ['mag_r', 'f', 0., 'radius of curvature of central trajectory [m] (for simulating e.g. quadrupole component integrated to a bending magnet; effective if > 0)'],
        ['mag_zc', 'f', 0., 'center longitudinal position [m]'],

#---Arbitrary Magnetic Field
        ['mag_fn', 's', '', 'input 3D magnetic field file name'],
        ['mag_ior', 'f', 3, 'interpolation order to use for tabulated magnetic field (1<= mag_ior <= 3)'],

#---Coherent Gaussian Beam
        ['gbm_pen', 'f', 0., 'energy per pulse [J] (needs to be >0 for the Gaussian beam to be set up)'],
        ['gbm_x', 'f', 0.0, 'average horizontal coordinates of waist [m]'],
        ['gbm_y', 'f', 0.0, 'average vertical coordinates of waist [m]'],
        ['gbm_z', 'f', 0.0, 'average longitudinal coordinate of waist [m]'],
        ['gbm_xp', 'f', 0.0, 'average horizontal angle at waist [rad]'],
        ['gbm_yp', 'f', 0.0, 'average verical angle at waist [rad]'],
        ['gbm_ave', 'f', 9000.0, 'average photon energy [eV]'],
        ['gbm_rep', 'f', 1, 'rep. rate [Hz]'],
        ['gbm_pol', 'f', 1, 'polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left'],
        ['gbm_sx', 'f', 9.78723e-06, 'rms beam size vs horizontal position [m] at waist (for intensity)'],
        ['gbm_sy', 'f', 9.78723e-06, 'rms beam size vs vertical position [m] at waist (for intensity)'],
        ['gbm_st', 'f', 1e-13, 'rms pulse duration [s] (for intensity)'],
        ['gbm_mx', 'f', 0, 'transverse Gauss-Hermite mode order in horizontal direction'],
        ['gbm_my', 'f', 0, 'transverse Gauss-Hermite mode order in vertical direction'],
        ['gbm_ca', 's', 'c', 'treat _sigX, _sigY as sizes in [m] in coordinate representation (_presCA="c") or as angular divergences in [rad] in angular representation (_presCA="a")'],
        ['gbm_ft', 's', 't', 'treat _sigT as pulse duration in [s] in time domain/representation (_presFT="t") or as bandwidth in [eV] in frequency domain/representation (_presFT="f")'],

#---Calculation Types
    #Magnetic Field
        ['mf', '', '', 'calculate/tabulate magnetic field (from other magnetic field input)', 'store_true'],
        ['mf_x', 'f', 0., 'horizontal center position [m] of the magnetic field calculation mesh'],
        ['mf_rx', 'f', 0., 'range of horizontal position of the magnetic field calculation mesh'],
        ['mf_nx', 'i', 1, 'number of points vs horizontal position of the magnetic field calculation mesh'],
        ['mf_y', 'f', 0., 'vertical center position [m] of the magnetic field calculation mesh'],
        ['mf_ry', 'f', 0., 'range of vertical position of the magnetic field calculation mesh'],
        ['mf_ny', 'i', 1, 'number of points vs vertical position of the magnetic field calculation mesh'],
        ['mf_z', 'f', 0., 'longitudinal center position [m] of the magnetic field calculation mesh'],
        ['mf_rz', 'f', 0., 'range of longitudinal position of the magnetic field calculation mesh'],
        ['mf_nz', 'i', 1, 'number of points vs longitudinal position of the magnetic field calculation mesh'],
        ['mf_ior', 'f', 3, 'interpolaton order to use for the resulting magnetic field'],
        ['mf_mag', 'i', 2, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
        ['mf_fn', 's', 'res_mag_fld.dat', 'file name for saving calculated magnetic field data'],
        ['mf_pl', 's', 'bxz,byz,bxby', 'plot the resulting magnetic field in graph(s): ""- dont plot, otherwise the string should list the magnetic field components to plot'],

    #Electron Trajectory
        ['tr', '', '', 'calculate electron trajectory', 'store_true'],
        ['tr_cti', 'f', 0., 'initial time moment (c*t) for electron trajectory calculation [m]'],
        ['tr_ctf', 'f', 0., 'final time moment (c*t) for electron trajectory calculation [m]'],
        ['tr_np', 'f', 50000, 'number of points for trajectory calculation'],
        ['tr_mag', 'i', 2, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
        ['tr_fn', 's', 'res_trj.dat', 'file name for saving calculated trajectory data'],
        ['tr_pl', 's', 'xz,xpz,yz,ypz,xy,xpyp', 'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

    #Single-Electron Spectrum vs Photon Energy
        ['ss', '', '', 'calculate single-e spectrum vs photon energy', 'store_true'],
        ['ss_ei', 'f', 100., 'initial photon energy [eV] for single-e spectrum vs photon energy calculation'],
        ['ss_ef', 'f', 20000., 'final photon energy [eV] for single-e spectrum vs photon energy calculation'],
        ['ss_ne', 'i', 10000, 'number of points vs photon energy for single-e spectrum vs photon energy calculation'],
        ['ss_x', 'f', 0., 'horizontal position [m] for single-e spectrum vs photon energy calculation'],
        ['ss_y', 'f', 0., 'vertical position [m] for single-e spectrum vs photon energy calculation'],
        ['ss_meth', 'i', 1, 'method to use for single-e spectrum vs photon energy calculation: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
        ['ss_prec', 'f', 0.01, 'relative precision for single-e spectrum vs photon energy calculation (nominal value is 0.01)'],
        ['ss_zi', 'f', 0., 'initial longitudinal position along electron trajectory for SR spectrum calculation (effective if ss_zi < ss_zf)'],
        ['ss_zf', 'f', 0., 'final longitudinal position along electron trajectory for SR spectrum calculation (effective if ss_zi < ss_zf)'],
        ['ss_pol', 'i', 6, 'polarization component to extract after spectrum vs photon energy calculation: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        ['ss_mag', 'i', 1, 'magnetic field to be used for single-e spectrum vs photon energy calculation: 1- approximate, 2- accurate'],
        ['ss_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
        ['ss_u', 'i', '1', 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
        ['ss_fn', 's', 'res_spec_se.dat', 'file name for saving calculated single-e spectrum vs photon energy'],
        ['ss_pl', 's', 'e', 'plot the resulting single-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

    #Coherent Gaussian Beam Spectrum vs Photon Energy or Time
        ['gs', '', '', 'calculate Gaussian beam spectrum vs photon energy or time (has priority over "ss" if both Gaussian beam and e-beam + magnetic field are defined)', 'store_true'],

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
        ['sm_smpf', 'f', -1, 'sampling factor for calculation of intensity distribution vs horizontal and vertical position (active if >0)'],
        ['sm_mag', 'i', 1, 'magnetic field to be used for calculation of multi-e spectrum spectrum or intensity distribution: 1- approximate, 2- accurate'],
        ['sm_hi', 'i', 1, 'initial UR spectral harmonic to be taken into accountfor multi-e spectrum vs photon energy calculation'],
        ['sm_hf', 'i', 15, 'final UR spectral harmonic to be taken into accountfor multi-e spectrum vs photon energy calculation'],
        ['sm_prl', 'f', 1., 'longitudinal integration precision parameter for multi-e spectrum vs photon energy calculation'],
        ['sm_pra', 'f', 1., 'azimuthal integration precision parameter for multi-e spectrum vs photon energy calculation'],
        ['sm_meth', 'i', 1, 'method to use for spectrum vs photon energy calculation in case of arbitrary input magnetic field: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler", -1- dont use this accurate integration method (rather use approximate if possible)'],
        ['sm_prec', 'f', 0.01, 'relative precision for spectrum vs photon energy calculation in case of arbitrary input magnetic field (nominal value is 0.01)'],
        ['sm_nm', 'i', 1000000, 'number of macro-electrons for calculation of spectrum in case of arbitrary input magnetic field'],
        ['sm_na', 'i', 10, 'number of macro-electrons to average on each node at parallel (MPI-based) calculation of spectrum in case of arbitrary input magnetic field'],
        ['sm_ns', 'i', 10, 'saving periodicity (in terms of macro-electrons) for intermediate intensity at calculation of multi-electron spectrum in case of arbitrary input magnetic field'],
        ['sm_type', 'i', 1, 'calculate flux (=1) or flux per unit surface (=2)'],
        ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total, 7- All Polarizations (0..5) in separate files'],
        #['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        ['sm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
        ['sm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],
        ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
        ['sm_pl', 's', 'e', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],
        ['sm_fbk', '', '', 'create backup file(s) with multi-e spectrum (only is it is calculated using the macro-particle method)', 'store_true'],

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
        ['pw_zi', 'f', 0., 'initial longitudinal position along electron trajectory for power density calculation (effective if pw_zi < pw_zf)'],
        ['pw_zf', 'f', 0., 'final longitudinal position along electron trajectory for power density calculation (effective if pw_zi < pw_zf)'],
        ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
        ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
        ['pw_pl', 's', 'xy', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

    #Undulator "operation table", i.e. dependence of gap (and phase) on photon energy (for a given polarization)
        ['ut', '', '', 'calculate undulator "operation table", i.e. dependence of gap (and phase) on photon energy (for a given polarization)', 'store_true'],
        ['ut_fn', 's', 'und_oper_table.dat', 'file name for saving calculated undulator operation table data'],

    #Single-Electron Intensity distribution vs horizontal and vertical position
        ['si', '', '', 'calculate single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position', 'store_true'],

    #Coherent Gaussian Beam Intensity distribution vs horizontal and vertical position
        ['gi', '', '', 'calculate coherent Gaussian beam intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position (has priority over "si" if both Gaussian beam and e-beam + magnetic field are defined)', 'store_true'],

    #Single-Electron Wavefront Propagation
        ['ws', '', '', 'calculate single-electron (/ fully coherent) wavefront propagation', 'store_true'],

    #Coherent Gaussian Beam Wavefront Propagation
        ['wg', '', '', 'calculate coherent Gaussian beam wavefront propagation (has priority over "si" if both Gaussian beam and e-beam + magnetic field are defined)', 'store_true'],
        
    #Multi-Electron (partially-coherent) Wavefront Propagation
        ['wm', '', '', 'calculate multi-electron (/ partially coherent) wavefront propagation', 'store_true'],

    #General Wavefront parameters (used at several different calculations)
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
        ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"'],
        ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
        ['w_zi', 'f', 0., 'initial longitudinal position [m] along electron trajectory for SR calculation (effective if w_zi < w_zf)'],
        ['w_zf', 'f', 0., 'final longitudinal position [m] along electron trajectory for SR calculation (effective if w_zi < w_zf)'],
        ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],
        ['w_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
        ['w_u', 'i', '1', 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
        ['w_wr', 'f', 0., 'wavefront radius to set (is taken into account if != 0) [m]; this parameter may be important for subsequent wavefront propagation simulations; by default, it is set by a function calculating the initial wavefront; however, it can also be set manually using this variable'],
        ['w_wre', 'f', 0., 'wavefront radius error (is taken into account if != 0) [m]; this parameter may be important for subsequent wavefront propagation simulations; by default, it is set by a function calculating the initial wavefront; however, it can also be set manually using this variable'],
        
        ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
        ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],

        #['mi_pol', 'i', 6, 'polarization component: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        #['mi_fn', 's', 'res_int_me.dat', 'file name for saving calculated multi-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position, for one or several photon energies'],

        ['ws_fne', 's', '', 'file name for saving initial electric field wavefront (before propagation)'],
        ['ws_fnep', 's', '', 'file name for saving electric field wavefront after propagation'],
        ['ws_fnei', 's', '', 'file name for input electric field for further eventual propagation; if this file name is supplied, the initial electric field will not be calculated, but loaded from this file'],
        ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],

        ['ws_pl', 's', 'xy', 'plot the propagated radiaiton intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
        ['ws_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters (/ Intensity distribution): coordinate (0) or angular (1)'],
        ['si_pl', 's', 'xy', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

        ['wm_nm', 'i', 100000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
        ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node at parallel (MPI-based) calculation of multi-electron wavefront propagation'],
        ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
        ['wm_ms', 'i', 0, 'index of initial coherent wavefront (/ mode) to start partially-coherent calculation from'],

        ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y'],
        ['wm_pol', 'i', 6, 'polarization component of to extract after calculation of multi-electron intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total, 7- All Main Polarizations (0-5)'], #OC28122023

        ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
        ['wm_x0', 'f', 0, 'horizontal center position for mutual intensity cut calculation'],
        ['wm_y0', 'f', 0, 'vertical center position for mutual intensity cut calculation'],
        ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
        ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
        ['wm_am', 'i', 0, 'multi-electron integration approximation method: 0- no approximation (use the standard 5D integration method), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance'],
        ['wm_fni', 's', 'res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],
        ['wm_ff', 's', 'ascii', 'format of file name for saving propagated multi-e intensity distribution vs horizontal and vertical position (ascii and hdf5 supported)'], #OC25022021

        ['wm_nmm', 'i', 1, 'number of MPI Masters to use'], #OC16042021
        ['wm_ncm', 'i', 100, 'number of coherent modes to calculate'], #OC02072021
        ['wm_acm', 's', 'SP', 'coherent mode decomposition algorithm to be used (supported algorithms are: "SP" for SciPy, "SPS" for SciPy Sparse, "PM" for Primme, based on names of software packages)'], #OC02072021
        ['wm_nop', '', '', 'switch forcing to do calculations ignoring any optics defined (by set_optics function)', 'store_true'], #OC03072021

        ['wm_fnmi', 's', '', 'file name of input cross-spectral density / mutual intensity; if this file name is supplied, the initial cross-spectral density (for such operations as coherent mode decomposition) will not be calculated, instead, it will be taken from that file (or a set of files).'], #OC02072021
        ['wm_fmis', 'i', -1, 'initial intex of partial cross-spectral density / mutual intensity file (is taken into account if >=0)'], #OC10102021
        #['wm_fmin', 'i', 1, 'number of partial cross-spectral density / mutual intensity files (is taken into account if >1)'], #OC10102021

        ['wm_fncm', 's', '', 'file name of input coherent modes; if this file name is supplied, the eventual partially-coherent radiation propagation simulation will be done based on propagation of the coherent modes from that file.'], #OC02072021

        ['wm_fbk', '', '', 'create backup file(s) with propagated multi-e intensity distribution vs horizontal and vertical position and other radiation characteristics', 'store_true'],
        ['wm_pl', 's', 'xy', 'plot the propagated radiaiton intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

        #['ws_fn', 's', '', 'file name for saving single-e (/ fully coherent) wavefront data'],
        #['wm_fn', 's', '', 'file name for saving multi-e (/ partially coherent) wavefront data'],

    #Optics parameters
        ['op_r', 'f', 30., 'longitudinal position of the first optical element [m]'],
        ['op_dp', 'f', 0., 'length of drift space to be applied after propagation through a beamline [m]'],
        ['op_fno', 's', '', 'file name for saving orientations of optical elements in the lab frame'],
        ['op_rnd', '', '', 'randomize parameters of optical elements or not (the randomizaiton params have to be defined separately)', 'store_true'], #OC24042020

    #Detector parameters
        ['d_x', 'f', 0., 'central horizontal position [m] of detector active area'],
        ['d_rx', 'f', 0., 'range of horizontal position [m] of detector active area (should be >0 to be taken into account)'],
        ['d_nx', 'i', 0, 'number of pixels vs horizontal position (should be >0 to be taken into account)'],
        ['d_dx', 'f', 0., 'horizontal size of pixel [m]'],
        ['d_y', 'f', 0., 'central vertical position [m] of detector active area'],
        ['d_ry', 'f', 0., 'range of vertical position [m] of detector active area (should be >0 to be taken into account)'],
        ['d_ny', 'i', 0, 'number of pixels vs vertical position (should be >0 to be taken into account)'],
        ['d_dy', 'f', 0., 'vertical size of pixel [m]'],
        ['d_or', 'f', 1, 'interpolation order (i.e. order of polynomials to be used at 2D interpolation)'],
        ['d_ifn', 's', '', 'file name with detector spectral efficiency data (on 1D mesh vs photon energy: _eStart, _eFin, _ne)'],
        #['d_ifn', 's', 'in_det_spec_eff.dat', 'file name with detector spectral efficiency data (on 1D mesh vs photon energy: _eStart, _eFin, _ne)'],

    #Optimization
        ['om', '', '', 'perform optimization of beamline parameters', 'store_true'],
        ['om_pn', '', ['op_VFM_r', 'op_HFM_r', 'op_S0_Dy', 'op_S0H_Dx'], 'names of parameters to be optimized'],
        ['om_iv', 'f', [20000., 20000., 1e-05, 1e-05], 'initial values of parameters to be optimized'],
        ['om_lm', 'f', [[10000.,30000.], [10000.,30000.], [1.e-06,0.001], [1.e-06,0.001]], 'limiting (i.e. min. amd max. possible) values of parameters to be optimized'],
        #['om_ls', 'f', [[10000.,10000.,1.e-06,0.001], [20000.,10000.,1.e-06,0.001]], 'explicit list of parameter values to try / use (for basic parameter scan method)'],

        ['om_cw', 'f', {'xFWHM':1, 'yFWHM':1, 'xpFWHM':0, 'ypFWHM':0, 'xCL':0, 'yCL':0, 'iM':0, 'iD':0, 'xFWFM':1.5, 'yFWFM':3., 'xpFWFM':0, 'ypFWFM':0}, 'weights of different criterions the optimization should be performed for: "xFWHM"- horizontal FWHM spot size, "yFWHM"- vertical FWHM spot size, "xpFWHM"- horizontal FWHM angular divergence, "ypFWHM"- vertical FWHM angular divergence, "xCL"- horizontal coherence length, "yCL"- vertical coherence length, "iM"- peak intensity, "iD"- given intensity distribution, "xFWFM"- full horizontal width at fraction of maximum, "yFWFM"- full vertical width at fraction of maximum, "xpFWFM"- full horizontal angular divergence at fraction of maximum, , "vpFWFM"- full vertical angular divergence at fraction of maximum'],
        ['om_ct', 'f', {'xFWHM':0.7e-06, 'yFWHM':0.7e-06, 'xpFWHM':0, 'ypFWHM':0, 'xCL':0, 'yCL':0, 'iM':1.e+17, 'iD':0, 'xFWFM':(0.7e-06*1.52379), 'yFWFM':(0.7e-06*1.52379), 'xpFWFM':0, 'ypFWFM':0}, 'target values of different criterions / figures of merit for the optimization, in the same order as defined by om_cw'],
        ['om_cn', 'f', {'xFWHM':0.7e-06, 'yFWHM':0.7e-06, 'xpFWHM':0, 'ypFWHM':0, 'xCL':0, 'yCL':0, 'iM':1.e+17, 'iD':0, 'xFWFM':(0.7e-06*1.52379), 'yFWFM':(0.7e-06*1.52379), 'xpFWFM':0, 'ypFWFM':0}, 'nominal values of different criterions / figures of merit for the optimization, in the same order as defined by om_cw'],
        ['om_ce', 'f', {'xFWFM':0.2, 'yFWFM':0.2, 'xpFWFM':0.2, 'ypFWFM':0.2}, 'extra data related to some optimization criterions, e.g. for *FWFM this is the fraction of maximum at which full width of intensity distribution should be considered'],

        ['om_fn', 's', '', 'input intensity distribution file name'], #?
        ['om_ft', 'i', 2, 'type of input intensity distribution: 1- flux per unit surface area, 2- azimuthally-averaged flux per unit surface area'],

        ['om_mt', 'i', 1, 'optimization method to use: 0- basic scan on parameter values, 1- ..., 2- ...'],

        ['om_mp', 'f', {'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 1000}, 'method-dependent optimization parameters (dictionary)'],
        #['om_mp', 'i', {'xtol': 1e-3, 'ftol': 1e-3, 'maxiter': 1000}, 'method-dependent optimization parameters (dictionary)'],

        ['om_pr', '', '1', 'print-out auxiliary information in the course of the optimization procedure', 'store_true'],
        ['om_fl', '', '1', 'save auxiliary information during the optimization procedure to a listing file', 'store_true'],
        
        #add more options here
    ]
    return varParamStd

#****************************************************************************
def srwl_uti_merge_options(_arOpt1, _arOpt2):
    """Merge arrays of options specified in _arOpt1 and _arOpt2, eliminating duplicates
    :param _arOpt1: list providing compact description of options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    :param _arOpt2: list providing compact description of extra options
    """
    nOpt1 = len(_arOpt1)
    nOpt2 = len(_arOpt2)
    arOptRes = copy(_arOpt1)
    
    for i in range(nOpt2):
        curOpt2 = _arOpt2[i]
        curOpt2nm = curOpt2[0]
        iOpt1 = nOpt1
        for j in range(nOpt1):
            curOpt1 = arOptRes[j]
            curOpt1nm = curOpt1[0]
            if(curOpt2nm == curOpt1nm):
                iOpt1 = j
                break
        curOpt2c = copy(curOpt2)
        if((iOpt1 >= 0) and (iOpt1 < nOpt1)):
            arOptRes[iOpt1] = curOpt2c
        else:
            arOptRes.append(curOpt2c)
    return arOptRes

#****************************************************************************
def srwl_uti_ext_options(_arOpt):
    """Extend array of standard options with the options specified in _arOpt, eliminating duplicates
    :param _arOpt: list providing compact description of extra options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    """
    return srwl_uti_merge_options(srwl_uti_std_options(), _arOpt)

#****************************************************************************
##def srwl_uti_parse_options(_descr): #OC08032016 (commented-out)
##    """Set and parse command-prompt options from a compact description provided in _descr
##    :param _descr: list providing compact description of all options; every element of this list is supposed to contain:
##        [0]: string containing option (/ variable) name
##        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
##        [2]: default value
##        [3]: string containing help / explanation of the option / variable
##        [4]: optional string describing formal action to be taken if option is fired
##    """
##
##    p = optparse.OptionParser()
##    nOpt = len(_descr)
##
##    listOptNamesPostParse = []
##    for i in range(nOpt):
##        curOpt = _descr[i]
##        
##        sTypeShort = curOpt[1]
##        sType = 'string'
##        if(sTypeShort == 'f'): sType = 'float'
##        elif(sTypeShort == 'i'): sType = 'int'        
##        #elif(sTypeShort == 's'): sType = 'string'
##
##        sAct = 'store'
##        if(len(curOpt) > 4): sAct = curOpt[4]
##
##        defVal = curOpt[2]
##        
##        optIsList = False
##        if(isinstance(defVal, list) or isinstance(defVal, array)): optIsList = True
##
##        if(optIsList):
##            sType = 'string'
##            listOptNamesPostParse.append(curOpt[0])
##
##        if(len(sTypeShort) <= 0):
##            p.add_option('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
##        else:
##            p.add_option('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)
##
##    v, args = p.parse_args()
##
##    #"post-parsing" list-type options
##    for i in range(len(listOptNamesPostParse)):
##        curOptName = listOptNamesPostParse[i]
##        valCurOpt = getattr(v, curOptName)
##
##        if((isinstance(valCurOpt, list) == False) and (isinstance(valCurOpt, array) == False)):
##            parsedVal = srwl_uti_parse_str2list(valCurOpt)
##            setattr(v, curOptName, parsedVal)
##    
##    return v

#****************************************************************************
#def _optparse(_descr, use_sys_argv=True, args=None):  #MR26022016, MR04032016
def srwl_uti_parse_options_obs(_descr, use_sys_argv=True, args=None): #OC08032016 #MR26022016, MR04032016
    """Set and parse command-prompt options from a compact description provided in _descr using optparse (OBSOLETE: deprecated since Python 2.7).
    :param _descr: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    :param use_sys_argv: a flag which manages use of sys.argv values in optparse.
    :param args: arbitrary arguments to be parsed, used when use_sys_argv is set to False.
    """
    import optparse

    p = optparse.OptionParser(None if use_sys_argv else __name__)
    nOpt = len(_descr)

    listOptNamesPostParse = []
    for i in range(nOpt):
        curOpt = _descr[i]

        sTypeShort = curOpt[1]
        sType = 'string'
        if (sTypeShort == 'f'):
            sType = 'float'
        elif (sTypeShort == 'i'):
            sType = 'int'
        # elif(sTypeShort == 's'): sType = 'string'

        sAct = 'store'
        if (len(curOpt) > 4): sAct = curOpt[4]

        defVal = curOpt[2]

        optIsList = False
        if (isinstance(defVal, list) or isinstance(defVal, array)): optIsList = True

        if (optIsList):
            sType = 'string'
            listOptNamesPostParse.append(curOpt[0])

        if (len(sTypeShort) <= 0):
            p.add_option('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
        else:
            p.add_option('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)

    #MR07032016:
    if use_sys_argv:
        v, args = p.parse_args() #MR07032016
    else:
        try:
            v, args = p.parse_args(args if args else []) #MR07032016
        except SystemExit as e:
            raise ValueError('Exit code: {}'.format(e))

    # "post-parsing" list-type options
    for i in range(len(listOptNamesPostParse)):
        curOptName = listOptNamesPostParse[i]
        valCurOpt = getattr(v, curOptName)

        if ((isinstance(valCurOpt, list) == False) and (isinstance(valCurOpt, array) == False)):
            parsedVal = srwl_uti_parse_str2list(valCurOpt)
            setattr(v, curOptName, parsedVal)

    return v

#****************************************************************************
#def _argparse(_descr, use_sys_argv=True, args=None): #MR26022016, MR04032016
def srwl_uti_parse_options(_descr, use_sys_argv=True, args=None): #OC08032016 #MR26022016, MR04032016
    """Set and parse command-prompt options from a compact description provided in _descr using argparse (recommended since Python 2.7).
    :param _descr: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    :param use_sys_argv: a flag which manages use of sys.argv values in argparse.
    :param args: arbitrary arguments to be parsed, used when use_sys_argv is set to False.
    """
    import argparse

    #_descr = srwl_uti_ext_options(deepcopy(_descr)) #OCTEST: consider adding this?

    p = argparse.ArgumentParser(None if use_sys_argv else __name__) #MR07032016
    nOpt = len(_descr)

    listOptNamesPostParse = []
    for i in range(nOpt):
        curOpt = _descr[i]

        sTypeShort = curOpt[1]
        sType = str
        if sTypeShort == 'f':
            sType = float
        elif sTypeShort == 'i':
            sType = int

        sAct = 'store'
        if (len(curOpt) > 4): sAct = curOpt[4]

        defVal = curOpt[2]

        optIsList = False
        if (isinstance(defVal, list) or isinstance(defVal, array)): optIsList = True

        if (optIsList):
            sType = str
            listOptNamesPostParse.append(curOpt[0])

        #curOpt[3] = curOpt[3].replace('%', '%%') # screen special '%' symbol
        if '%%' not in curOpt[3]: #MR14042018 (via GitHub pull req.)
            curOpt[3] = curOpt[3].replace('%', '%%') # screen special '%' symbol 
        
        if (len(sTypeShort) <= 0):
            p.add_argument('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
        else:
            p.add_argument('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)

    #MR07032016:
    if use_sys_argv:
        v = p.parse_args() #MR07032016
    else:
        try:
            v = p.parse_args(args if args else []) #MR07032016
        except SystemExit as e:
            raise ValueError('Exit code: {}'.format(e))

    # "post-parsing" list-type options
    for i in range(len(listOptNamesPostParse)):
        curOptName = listOptNamesPostParse[i]
        valCurOpt = getattr(v, curOptName)

        if not isinstance(valCurOpt, list) and not isinstance(valCurOpt, array): #MR07032016
            parsedVal = srwl_uti_parse_str2list(valCurOpt)
            setattr(v, curOptName, parsedVal)

    return v

'''
MR26022016: Here we can specify which parser to use for parsing the options. Argparse is used by default, but if
a user wants to execute optparse, the following environment variable has to be set up:
    SRWL_OPTPARSE=1

On Linux systems it can be set either in ~/.bashrc using:
    export SRWL_OPTPARSE=1
or just before the executed command:
    SRWL_OPTPARSE=1 python script.py -h

On Windows systems it can be set using "Environment Variables..." button in System Properties or via command line:
    set SRWL_OPTPARSE=1

If you wish to use a particular parser, edit the statement below as follows:
    srwl_uti_parse_options = _optparse  # for optparse
    srwl_uti_parse_options = _argparse  # for argparse

Various options can be specified including lists, e.g.:
    python script.py --op_S0_pp="[0, 0, 1, 0, 0, 5.0, 8.0, 2.5, 3.5, 0, 0, 0]"
'''

#srwl_uti_parse_options = srwl_uti_parse_options_obs if os.getenv('SRWL_OPTPARSE') else _argparse  # MR07032016
if(os.getenv('SRWL_OPTPARSE')): #OC08032016
    srwl_uti_parse_options = srwl_uti_parse_options_obs

#****************************************************************************
#OC25012019: This function is obsolete; this processing is taking place in SRWLBeamline.calc_all().
#The source type is defined by v.rs_tp string variable (instead of v.source_type).
#OC (previously): This function function may need to be moved to SRWLBeamline class or renamed (to have name with prefixes / "decorations")
def setup_source(v):  #MR20160617 - moved from Sirepo .jinja template
    mag = None
    if v.source_type in ['u', 't']:
        if v.source_type == 'u' or (not v.und_g or v.und_g == 0):
            v.und_b = 1
            if hasattr(v, 'und_g'):
                del v.und_g
            if hasattr(v, 'gbm_pen'):
                del v.gbm_pen
        elif v.source_type == 't' or (v.und_g and v.und_g > 0):
            if hasattr(v, 'gbm_pen'):
                del v.gbm_pen
            v.pw_mag = 2
            v.w_mag = 2
    elif v.source_type == 'g':
        pass
    elif v.source_type == 'm':
        mag = SRWLMagFldC()
        mag.arXc.append(0)
        mag.arYc.append(0)
        mag.arMagFld.append(SRWLMagFldM(
            v.mp_field,
            v.mp_order,
            v.mp_distribution,
            v.mp_len
        ))
        mag.arZc.append(v.mp_zc)
        if hasattr(v, 'gbm_pen'):
            del v.gbm_pen
    else:
        raise AssertionError('{}: unknown source_type'.format(v.source_type))

    return v.source_type, mag




